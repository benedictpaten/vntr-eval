#!/usr/bin/env python3
"""mafft_panel.py -- the full-panel arm of the mafft realigners.

For every region and mafft method this aligns the deduplicated FULL panel (work/panel/union/<id>.fa,
from `tools/panel.py union`: every distinct sequence of hprc.fa.gz and hap32.fa), keeps the full MSA
and its full-panel graph, and projects the MSA onto the 34 hap32 rows:

    tools/panel.py union            work/panel/union/<id>.fa + <id>.map.tsv (made here when missing)
    realign.align_fasta(METHOD)     the realigner's MSA-only mode: the same method definition, dedup,
                                    N masking and output checks as the hap32 arm (tools/realign.py)
    -> work/panel/<METHOD>/<id>.msa.fa.gz   the full MSA, rows named by distinct id (u0001..; the map
                                    TSV lists the haplotypes that carry each)
    -> work/panel/<METHOD>/<id>.gfa         its graph, one P line per hprc.fa.gz record (panel.py
                                    panel-graph expands the distinct rows back to member names;
                                    every path's spelling is re-read and checked)
    -> candidates/<METHOD>__all/<id>.msa.fa the hap32 rows (panel.py project; all-gap columns dropped)
    -> candidates/<METHOD>__all/<id>.gfa    msa_graph.py, one P line per hap32.fa name; every path is
                                    re-read from the written file and asserted to spell hap32.fa
    -> candidates/<METHOD>__all/<id>.realign.json  how it was made (also written for timeouts,
                                    memouts, skips and failures, so a re-run skips them)
    -> results/realign_runtime.all.mafft.tsv   one row per (method, region); the latest run wins

Methods (tools/realign.py METHODS; every call adds --nuc --thread T --threadit 0):
    mafft_fftns2   --retree 2 --maxiterate 0                 progressive; the method run on every region
    mafft_fftnsi   --retree 2 --maxiterate 1000              FFT-NS-i: progressive + iterative refinement,
                                                            no pairwise stage (all-arm method; the mafft
                                                            script caps iterations at 16, and its own
                                                            `fftnsi` wrapper uses 2)
    mafft_linsi    --localpair --maxiterate 1000             L-INS-i  \\  every pair of distinct sequences
    mafft_einsi    --genafpair --maxiterate 1000 --ep 0      E-INS-i   > is aligned by full DP first:
    mafft_ginsi    --globalpair --maxiterate 1000            G-INS-i  /  O(N^2 L^2) time
Threads: 3 for the pairwise methods (dropped to 1 when the longest pair's DP would not fit the cap),
2 for the progressive ones (their progressive and refinement steps are single-threaded anyway).

Usage (from the repository root):

    python3 tools/mafft_panel.py run [--methods M1,M2,..] [--regions ID,..|all] [--stratum S]
        [--jobs 2] [--timeout 1800] [--mem-mb 12000] [--budget-mb 18000]
        [--threads-pairwise 3] [--threads-progressive 2] [--skip-after 2] [--skip-factor 3]
        [--no-skip-methods M,..] [--priority M=W,..] [--prefer-regions ID,.. --prefer-methods M,..]
        [--deadline-h H] [--busy-pids PID,..] [--exclude M:ID,..]
        [--force | --retry STATUS[,STATUS]] [--workdir DIR]
    python3 tools/mafft_panel.py one METHOD REGION [--threads T --timeout S --mem-mb MB]
    python3 tools/mafft_panel.py predict [--methods ..]          # the cost/memory model per job
    python3 tools/mafft_panel.py status                          # counts per method, status, stratum
    python3 tools/mafft_panel.py mark-skipped --methods .. --reason TEXT   # record every job that
                                                  # has no result yet as 'skipped'

Scheduling. Jobs are separate processes, at most --jobs at a time, started in increasing order of
rank = predicted aligner seconds x the method's --priority weight (default fftns2 0.25, linsi and
fftnsi 1, einsi and ginsi 1.5), x --prefer-factor for --prefer-regions (a first look's regions).
So within one method regions always start smallest first, and a slow region never starves the rest.
A job starts only when the predicted memory of the running jobs plus its own stays within
--budget-mb and the OS reports that much available; when nothing is running the next job starts
regardless (its RSS cap still applies). A running job holds <id>.realign.json.running (its pid), so
a restarted scheduler neither repeats nor double-starts it; --busy-pids makes the new scheduler
count jobs an older one started as occupying slots. After --deadline-h no job starts, and every
pending job is recorded as skipped ("batch wall budget").

Not attempted, and recorded as such (status 'skipped', reason in `note`):
  * a job whose method already timed out on --skip-after (2) regions with a smaller cost measure
    (pair cells for the pairwise methods, N x L1 x L2 for the progressive ones);
  * a job predicted to need more than --skip-factor (3) x --timeout seconds;
  * a pairwise job predicted to need more than --mem-mb even with one thread (status 'memout',
    "predicted" in the note);
  * whatever was still pending at --deadline-h.
A job the aligner runs past --timeout is 'timeout'; past --mem-mb (sampled RSS of its process
group) 'memout'.

Cost model (predicted aligner seconds; distinct sequences after N masking, lengths L_i, L1 >= L2
the two longest):
    pairwise    : 8 ns x pair_cells / threads, pair_cells = sum_{i<j} L_i L_j. Measured 6.3-8.2 ns x
                  thread on the full panel's small regions; it underestimates regions with hundreds
                  of distinct sequences, where the refinement adds (E-INS-i L004989: 813 s vs 131).
    fftns2      : 1e-8 s x N x L1 x L2 (fitted on the first 51 runs: median 3.3e-9, hotspots to
                  2.1e-8; an upper-side value so the --skip-factor rule refuses only the largest).
    fftnsi      : 2e-8 s x N x L1 x L2.
Memory model (MB), used to reserve memory (and, for pairwise jobs, to refuse or drop to 1 thread):
    pairwise    : 400 + 20 B x L1 x L2 x threads (the DP matrices of the longest pair, one per
                  thread; hap32 arm: 10.5-14.5 GB at 1 thread and L1 20-34 kb);
    progressive : min(400 + 55 B x L1 x L2, 8500) (see predict_mb).
"""
import argparse
import collections
import csv
import datetime
import fcntl
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402
import msa_graph  # noqa: E402
import panel  # noqa: E402
import realign  # noqa: E402

METHODS = ('mafft_fftns2', 'mafft_fftnsi', 'mafft_linsi', 'mafft_einsi', 'mafft_ginsi')
PAIRWISE = ('mafft_linsi', 'mafft_einsi', 'mafft_ginsi')
RUNTIME_TSV = os.path.join(config.RESULTS_DIR, 'realign_runtime.all.mafft.tsv')
RUNTIME_COLS = ['region_id', 'method', 'panel', 'n_distinct', 'seconds', 'status',
                'stratum', 'span_bp', 'n_panel_records', 'n_hap32_only', 'max_len', 'len2', 'total_bp',
                'pair_cells', 'predicted_s', 'predicted_mb', 'threads', 'peak_rss_mb', 'mem_cap_mb',
                'timeout_s', 'align_s', 'panel_graph_s', 'project_s', 'graph_s', 'columns_full',
                'nodes_full', 'edges_full', 'columns_proj', 'nodes_proj', 'nodes_per_kb_proj',
                'frac_nodes_1bp_proj', 'masked_runs', 'masked_bp', 'host_load', 'finished', 'note']
DEFAULT_TIMEOUT = 1800
DEFAULT_MEM_MB = 12000
DEFAULT_BUDGET_MB = 20000
C_PAIRWISE = 8e-9          # s per pair cell per thread (all of L-/E-/G-INS-i; hap32 arm)
PROGRESSIVE_MAX_MB = 8500.0
C_PROGRESSIVE = {'mafft_fftns2': 1e-8, 'mafft_fftnsi': 2e-8}    # s per N x L1 x L2 (first 51 full-panel runs)
T0 = time.time()


def log(msg):
    sys.stderr.write('[mafft_panel %7.1fs] %s\n' % (time.time() - T0, msg))
    sys.stderr.flush()


# ---------------------------------------------------------------- regions and paths

def union_files(rid):
    return (os.path.join(panel.UNION_DIR, rid + '.fa'), os.path.join(panel.UNION_DIR, rid + '.map.tsv'))


def ensure_union(rd, rid):
    fa, mp = union_files(rid)
    if not (os.path.exists(fa) and os.path.exists(mp)):
        panel.union_to_files(rd, fa, mp)
    return fa, mp


def out_paths(method, rid, cand_root=None, panel_root=None):
    cd = os.path.join(cand_root or config.CANDIDATES_DIR, method + '__all')
    pd = os.path.join(panel_root or panel.PANEL_DIR, method)
    return {'cand_dir': cd, 'panel_dir': pd,
            'msa': os.path.join(cd, rid + '.msa.fa'), 'gfa': os.path.join(cd, rid + '.gfa'),
            'json': os.path.join(cd, rid + '.realign.json'),
            'full_msa': os.path.join(pd, rid + '.msa.fa.gz'), 'full_gfa': os.path.join(pd, rid + '.gfa')}


def pid_alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


def previous_status(method, rid, cand_root=None, panel_root=None):
    """Status of the last attempt ('running' while a live job holds the region's marker)."""
    o = out_paths(method, rid, cand_root, panel_root)
    try:
        with open(o['json'] + '.running') as f:
            if pid_alive(f.read().strip()):
                return 'running'
    except OSError:
        pass
    try:
        with open(o['json']) as f:
            j = json.load(f)
    except (OSError, ValueError):
        return None
    st = j.get('status')
    if st == 'ok' and not all(os.path.exists(o[k]) for k in ('msa', 'gfa', 'full_msa', 'full_gfa')):
        return None
    return st


def region_info(rid):
    """Sizes of a region's union: raw lengths, and the aligner's (N runs masked) lengths."""
    rd = realign.region_dir(rid)
    fa, mp = ensure_union(rd, rid)
    lens, mlens, n_only, n_rec = [], [], 0, 0
    for r in panel.read_map(mp):
        lens.append(int(r['length']))
        mlens.append(int(r['length']) - int(r['n_bp_N']))
        if int(r['weight']) == 0:
            n_only += 1
        n_rec += int(r['weight'])
    stratum, span = realign.region_meta(rd)
    s = sorted(mlens, reverse=True)
    tot = sum(mlens)
    return {'region_id': rid, 'rd': rd, 'stratum': stratum, 'span_bp': span, 'n_distinct': len(lens),
            'n_hap32_only': n_only, 'n_panel_records': n_rec, 'max_len': max(lens),
            'len2': s[1] if len(s) > 1 else 0, 'total_bp': sum(lens), 'mlens': mlens,
            'pair_cells': (tot * tot - sum(x * x for x in mlens)) // 2}


# ---------------------------------------------------------------- models

def predict_s(method, info, threads):
    s = sorted(info['mlens'], reverse=True)
    l1, l2 = s[0], (s[1] if len(s) > 1 else s[0])
    if method in PAIRWISE:
        return C_PAIRWISE * info['pair_cells'] / max(1, threads)
    return C_PROGRESSIVE[method] * len(s) * l1 * l2


def cost_of(method, info):
    """The size measure a method's run time grows with: pair cells for the all-pairs methods,
    N x L1 x L2 (distinct sequences x the two longest) for the progressive ones."""
    if method in PAIRWISE:
        return info['pair_cells']
    s = sorted(info['mlens'], reverse=True)
    return len(s) * s[0] * (s[1] if len(s) > 1 else s[0])


COST_UNIT = {True: 'pair cells', False: 'N x L1 x L2'}


def predict_mb(method, info, threads):
    """Predicted peak RSS (MB). Pairwise methods: the DP matrices of the longest pair, one per
    thread (well founded on the hap32 arm; a job predicted over the cap is not run). Progressive
    methods: 55 B x L1 x L2 (L005990 full panel: 7.4 GB at 11.7 kb, the profiles being ~2.7x
    longer than the longest sequence), at most 8.5 GB (the largest of 228 FFT-NS-2 / FFT-NS-i
    full-panel runs peaked at 7.95 GB, up to 39.6 kb sequences: mafft's own memory saving takes
    over); used only to reserve memory in the scheduler, never to refuse a job -- the RSS cap
    decides."""
    s = sorted(info['mlens'], reverse=True)
    l1, l2 = s[0], (s[1] if len(s) > 1 else s[0])
    if method in PAIRWISE:
        return 400 + 20.0 * l1 * l2 * threads / 2 ** 20
    return min(400 + 55.0 * l1 * l2 / 2 ** 20, PROGRESSIVE_MAX_MB)


def plan_threads(method, info, a):
    """Threads for a job, and its predicted MB: pairwise jobs drop to 1 thread when the DP of the
    longest pair would not fit the cap at --threads-pairwise (as the hap32 arm's re-runs did)."""
    if method in PAIRWISE:
        t = a.threads_pairwise
        if predict_mb(method, info, t) > a.mem_mb:
            t = 1
    else:
        t = a.threads_progressive
    return t, predict_mb(method, info, t)


# ---------------------------------------------------------------- runtime table

def update_runtime(row, path=RUNTIME_TSV):
    """Replace the (method, region_id) row of the runtime table, under a file lock."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + '.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        rows = []
        if os.path.exists(path):
            with open(path) as f:
                rows = list(csv.DictReader(f, delimiter='\t'))
        rows = [r for r in rows if (r.get('method'), r.get('region_id')) != (row['method'], row['region_id'])]
        rows.append(row)
        order = dict((m, i) for i, m in enumerate(METHODS))
        rows.sort(key=lambda r: (order.get(r['method'], 99), r['method'], r['region_id']))
        tmp = path + '.tmp%d' % os.getpid()
        with open(tmp, 'w', newline='') as f:
            w = csv.DictWriter(f, RUNTIME_COLS, delimiter='\t', extrasaction='ignore', lineterminator='\n')
            w.writeheader()
            for r in rows:
                w.writerow(dict((c, r.get(c, '')) for c in RUNTIME_COLS))
        os.replace(tmp, path)
        fcntl.flock(lk, fcntl.LOCK_UN)


def read_runtime(path=RUNTIME_TSV):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


# ---------------------------------------------------------------- one job

def check_paths(gfa, fasta_recs):
    """Re-read gfa and assert that each named path spells its sequence (and every name has a path)."""
    seqs, _, paths = panel.read_gfa_paths(gfa)
    want = collections.OrderedDict((n, s.upper()) for n, s in fasta_recs)
    bad = []
    for n, s in want.items():
        w = paths.get(n)
        if w is None or ''.join(seqs[x] if x > 0 else panel.revcomp(seqs[-x]) for x in w) != s:
            bad.append(n)
    extra = [n for n in paths if n not in want]
    if bad or extra:
        raise RuntimeError('%s: %d paths do not spell their sequence (%s), %d unexpected paths'
                           % (gfa, len(bad), ', '.join(bad[:3]), len(extra)))
    return len(want)


def base_row(method, info_r, status, **kw):
    row = {'region_id': info_r['region_id'], 'method': method, 'panel': 'all',
           'n_distinct': info_r['n_distinct'], 'status': status, 'stratum': info_r['stratum'],
           'span_bp': info_r['span_bp'], 'n_panel_records': info_r['n_panel_records'],
           'n_hap32_only': info_r['n_hap32_only'], 'max_len': info_r['max_len'], 'len2': info_r['len2'],
           'total_bp': info_r['total_bp'], 'pair_cells': info_r['pair_cells'],
           'finished': datetime.datetime.now().isoformat(timespec='seconds'),
           'host_load': round(os.getloadavg()[0], 1)}
    row.update(kw)
    return row


def write_json(path, info):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + '.tmp%d' % os.getpid()
    with open(tmp, 'w') as f:
        json.dump(info, f, indent=1)
    os.replace(tmp, path)


def record_not_run(method, rid, status, note, predicted_s='', predicted_mb='', threads='', a=None,
                   record_runtime=True):
    """A job that is not run (skipped or predicted memout): .realign.json + runtime row."""
    info_r = region_info(rid)
    o = out_paths(method, rid, getattr(a, 'cand_root', None), getattr(a, 'panel_root', None))
    for k in ('msa', 'gfa', 'full_msa', 'full_gfa'):
        if os.path.exists(o[k]):
            os.remove(o[k])
    info = collections.OrderedDict([
        ('method', method), ('panel', 'all'), ('region_id', rid), ('stratum', info_r['stratum']),
        ('span_bp', info_r['span_bp']), ('status', status), ('message', note),
        ('n_distinct', info_r['n_distinct']), ('max_len', info_r['max_len']), ('len2', info_r['len2']),
        ('total_bp', info_r['total_bp']), ('pair_cells', info_r['pair_cells']),
        ('predicted_s', predicted_s), ('predicted_mb', predicted_mb), ('threads', threads),
        ('seconds', 0.0), ('finished', datetime.datetime.now().isoformat(timespec='seconds')),
        ('mafft_panel_version', 1)])
    write_json(o['json'], info)
    row = base_row(method, info_r, status, seconds=0.0, predicted_s=predicted_s, predicted_mb=predicted_mb,
                   threads=threads, mem_cap_mb=getattr(a, 'mem_mb', ''), timeout_s=getattr(a, 'timeout', ''),
                   note=note[:300])
    if record_runtime:
        update_runtime(row)
    return row


def run_job(method, rid, threads=2, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB, cand_root=None,
            panel_root=None, workroot=None, record_runtime=True, panel_engine='native'):
    """Align region rid's full-panel union with method; write every output; return the runtime row."""
    m = realign.get_method(method)
    info_r = region_info(rid)
    rd = info_r['rd']
    fa, mp = union_files(rid)
    hap_fa = os.path.join(rd, 'hap32.fa')
    o = out_paths(method, rid, cand_root, panel_root)
    os.makedirs(o['cand_dir'], exist_ok=True)
    os.makedirs(o['panel_dir'], exist_ok=True)
    workroot = workroot or os.path.join(config.WORK_DIR, 'realign')
    os.makedirs(workroot, exist_ok=True)
    for k in ('msa', 'gfa', 'full_msa', 'full_gfa'):     # never leave an older result beside a new attempt
        if os.path.exists(o[k]):
            os.remove(o[k])
    t0 = time.time()
    marker = o['json'] + '.running'
    with open(marker, 'w') as f:
        f.write(str(os.getpid()))
    ps, pm = predict_s(method, info_r, threads), predict_mb(method, info_r, threads)
    info = collections.OrderedDict([
        ('method', method), ('panel', 'all'), ('region_id', rid),
        ('stratum', info_r['stratum']), ('span_bp', info_r['span_bp']), ('status', None),
        ('description', m.description), ('union_fa', realign._rel(fa)), ('union_map', realign._rel(mp)),
        ('n_distinct', info_r['n_distinct']), ('n_hap32_only', info_r['n_hap32_only']),
        ('n_panel_records', info_r['n_panel_records']), ('max_len', info_r['max_len']),
        ('len2', info_r['len2']), ('total_bp', info_r['total_bp']), ('pair_cells', info_r['pair_cells']),
        ('predicted_s', round(ps, 1)), ('predicted_mb', round(pm)), ('threads', threads),
        ('mem_cap_mb', mem_mb), ('timeout_s', timeout)])
    times = {}
    wd = tempfile.mkdtemp(prefix='mafft_panel.%s.%s.' % (method, rid), dir=workroot)
    try:
        full_msa = os.path.join(wd, 'full.msa.fa')
        al = realign.align_fasta(m, fa, full_msa, threads=threads, timeout=timeout, mem_mb=mem_mb, workdir=wd,
                                 dedup=True)
        times['align_s'] = al.get('align_s')
        info['align'] = al
        info['status'] = al['status']
        if al.get('message'):
            info['message'] = al['message']
        if al['status'] == 'ok':
            try:
                tg = time.time()
                tmp = o['full_msa'] + '.tmp%d' % os.getpid()
                with open(full_msa, 'rb') as fi, gzip.open(tmp, 'wb', compresslevel=6) as fo:
                    shutil.copyfileobj(fi, fo)
                os.replace(tmp, o['full_msa'])
                pg = panel.panel_graph(full_msa, mp, rd, o['full_gfa'], engine=panel_engine)
                hprc = [(n, s) for n, _, s in panel.read_fasta(os.path.join(rd, 'hprc.fa.gz'))]
                pg['paths_verified'] = check_paths(o['full_gfa'], hprc)
                pg['gfa'] = realign._rel(o['full_gfa'])
                pg['msa'] = realign._rel(o['full_msa'])
                info['panel_graph'] = pg
                times['panel_graph_s'] = round(time.time() - tg, 2)
                tp = time.time()
                tmp_msa = o['msa'] + '.tmp%d' % os.getpid()
                pr = panel.project(full_msa, mp, hap_fa, tmp_msa)
                info['projection'] = pr
                times['project_s'] = round(time.time() - tp, 2)
                th = time.time()
                tmp_gfa = o['gfa'] + '.tmp%d' % os.getpid()
                st = msa_graph.msa_to_gfa(tmp_msa, hap_fa, tmp_gfa, workdir=wd)
                hap = msa_graph.read_fasta(hap_fa)
                st['paths_verified'] = check_paths(tmp_gfa, hap)
                os.replace(tmp_msa, o['msa'])
                os.replace(tmp_gfa, o['gfa'])
                st['msa'] = realign._rel(o['msa'])
                st['gfa'] = realign._rel(o['gfa'])
                st['hap32'] = realign._rel(hap_fa)
                info['graph'] = st
                times['graph_s'] = round(time.time() - th, 2)
            except Exception as e:  # noqa: BLE001 -- recorded, the batch goes on
                info['status'] = 'graph_error'
                info['message'] = '%s: %s' % (type(e).__name__, e)
                for k in ('msa', 'gfa', 'full_msa', 'full_gfa'):
                    if os.path.exists(o[k]):
                        os.remove(o[k])
    finally:
        shutil.rmtree(wd, ignore_errors=True)
        for p in os.listdir(o['cand_dir']):
            if p.startswith(rid + '.') and '.tmp' in p:
                os.remove(os.path.join(o['cand_dir'], p))
    info.update(times)
    info['seconds'] = round(time.time() - t0, 2)
    info['finished'] = datetime.datetime.now().isoformat(timespec='seconds')
    info['host_load'] = round(os.getloadavg()[0], 1)
    info['mafft_panel_version'] = 1
    write_json(o['json'], info)
    if os.path.exists(marker):
        os.remove(marker)
    al = info.get('align') or {}
    aln = al.get('aligner') or {}
    pg = info.get('panel_graph') or {}
    g = info.get('graph') or {}
    pr = info.get('projection') or {}
    row = base_row(method, info_r, info['status'], seconds=info['seconds'], predicted_s=round(ps, 1),
                   predicted_mb=round(pm), threads=threads, peak_rss_mb=aln.get('peak_rss_mb', ''),
                   mem_cap_mb=mem_mb, timeout_s=timeout, align_s=times.get('align_s', ''),
                   panel_graph_s=times.get('panel_graph_s', ''), project_s=times.get('project_s', ''),
                   graph_s=times.get('graph_s', ''), columns_full=al.get('columns', ''),
                   nodes_full=pg.get('nodes', ''), edges_full=pg.get('edges', ''),
                   columns_proj=pr.get('columns_out', ''), nodes_proj=g.get('nodes', ''),
                   nodes_per_kb_proj=g.get('nodes_per_kb', ''), frac_nodes_1bp_proj=g.get('frac_nodes_1bp', ''),
                   masked_runs=al.get('masked_runs', ''), masked_bp=al.get('masked_bp', ''),
                   note=(info.get('message') or '')[:300].replace('\t', ' ').replace('\n', ' '))
    if record_runtime:
        update_runtime(row)
    return row


# ---------------------------------------------------------------- scheduler

def pick_regions(arg, stratum=None):
    rows = realign.region_list()
    if arg and arg != 'all':
        want = [x.strip() for x in arg.split(',') if x.strip()]
        known = set(r['region_id'] for r in rows)
        missing = [x for x in want if x not in known]
        if missing:
            raise SystemExit('mafft_panel: no region package for %s' % ', '.join(missing))
        rows = [r for r in rows if r['region_id'] in set(want)]
    if stratum:
        keep = set(stratum.split(','))
        rows = [r for r in rows if r.get('stratum') in keep]
    return [r['region_id'] for r in rows]


def job_argv(j, a):
    cmd = [sys.executable, os.path.abspath(__file__), 'one', j['method'], j['rid'], '--threads', str(j['threads']),
           '--timeout', str(a.timeout), '--mem-mb', str(a.mem_mb), '--panel-engine', a.panel_engine]
    for opt in ('cand_root', 'panel_root', 'workdir'):
        v = getattr(a, opt)
        if v:
            cmd += ['--' + opt.replace('_', '-'), v]
    if a.no_runtime:
        cmd.append('--no-runtime')
    return cmd


def skip_reason(j, timeouts, a):
    """Why a pending job is not attempted (None = run it)."""
    if a.skip_after > 0 and j['method'] not in a.no_skip:
        smaller = [(c, r) for c, r in timeouts.get(j['method'], []) if c <= j['cost']]
        if len(smaller) >= a.skip_after:
            smaller.sort()
            unit = COST_UNIT[j['method'] in PAIRWISE]
            return ('not attempted: %s timed out (%.0f s) on %d regions with a smaller %s, e.g. %s; '
                    'this region has %.3g' % (j['method'], a.timeout, len(smaller), unit,
                                             ', '.join('%s (%.3g)' % (r, c) for c, r in smaller[:3]), j['cost']))
    if a.skip_factor > 0 and j['pred_s'] > a.skip_factor * a.timeout and j['method'] not in a.no_skip:
        return 'not attempted: predicted %.0f s > %.0f x timeout %.0f s' % (j['pred_s'], a.skip_factor, a.timeout)
    return None


def parse_priority(text):
    out = {}
    for item in (text or '').split(','):
        if '=' in item:
            k, v = item.split('=', 1)
            out[k.strip()] = float(v)
    return out


def run_batch(a):
    methods = [x for x in a.methods.split(',') if x]
    for mname in methods:
        realign.get_method(mname)
    rids = pick_regions(a.regions, a.stratum)
    infos = dict((rid, region_info(rid)) for rid in rids)
    retry = set(x for x in (a.retry or '').split(',') if x)
    jobs = []
    timeouts = collections.defaultdict(list)       # method -> [(cost, rid)] of timeouts (past and new)
    for rid in rids:
        for mname in methods:
            if '%s:%s' % (mname, rid) in set((a.exclude or '').split(',')):
                continue
            prev = previous_status(mname, rid, a.cand_root, a.panel_root)
            if prev == 'timeout' and not (a.force or 'timeout' in retry):
                timeouts[mname].append((cost_of(mname, infos[rid]), rid))
            if prev and not a.force and prev not in retry:
                continue
            t, pm = plan_threads(mname, infos[rid], a)
            jobs.append({'method': mname, 'rid': rid, 'threads': t, 'pred_mb': pm,
                         'pred_s': predict_s(mname, infos[rid], t), 'cost': cost_of(mname, infos[rid]),
                         'reserve': min(pm, a.mem_mb)})
    prio = parse_priority(a.priority)
    a.no_skip = set(x for x in (a.no_skip_methods or '').split(',') if x)
    pref_r = set(x for x in (a.prefer_regions or '').split(',') if x)
    pref_m = set(x for x in (a.prefer_methods or '').split(',') if x) or set(methods)
    for j in jobs:
        j['rank'] = j['pred_s'] * prio.get(j['method'], 1.0) * (
            a.prefer_factor if (j['rid'] in pref_r and j['method'] in pref_m) else 1.0)
    jobs.sort(key=lambda j: (j['rank'], j['cost'], j['rid'], METHODS.index(j['method'])
                             if j['method'] in METHODS else 99))
    log('%d jobs (%d regions x %s); budget %d MB, cap %d MB, %d at a time, timeout %d s'
        % (len(jobs), len(rids), ','.join(methods), a.budget_mb, a.mem_mb, a.jobs, a.timeout))
    logdir = os.path.join(a.workdir or os.path.join(config.WORK_DIR, 'realign'), 'mafft_panel_logs')
    os.makedirs(logdir, exist_ok=True)
    running = {}
    counts = collections.Counter()
    pending = list(jobs)
    last_wait = 0
    stop_file = os.path.join(logdir, 'STOP')
    a.busy = [x for x in (a.busy_pids or '').split(',') if x]
    while pending or running:
        for p in list(running):
            if p.poll() is not None:
                j = running.pop(p)
                try:
                    with open(j['result']) as f:
                        row = json.load(f)
                except (OSError, ValueError):
                    row = {'status': 'crash', 'note': 'no result; exit %s; log %s' % (p.returncode, j['log'])}
                    if not a.no_runtime:
                        update_runtime(base_row(j['method'], infos[j['rid']], 'crash',
                                                seconds=round(time.time() - j['t0'], 1), threads=j['threads'],
                                                note=row['note']))
                if row.get('status') == 'timeout':
                    timeouts[j['method']].append((j['cost'], j['rid']))
                counts[row.get('status')] += 1
                log('%-13s %-9s %-11s %8ss align %8ss (pred %6.0f) rss %7s MB (pred %6.0f) t%d nodes full %-6s '
                    'proj %-6s %s' % (j['method'], j['rid'], row.get('status'), row.get('seconds', ''),
                                      row.get('align_s', ''), j['pred_s'], row.get('peak_rss_mb', ''),
                                      j['pred_mb'], j['threads'], row.get('nodes_full', ''),
                                      row.get('nodes_proj', ''), (row.get('note') or '')[:100]))
            elif time.time() - running[p]['t0'] > a.timeout * 2 + 1800:
                log('killing %s %s: job wall clock over %d s' % (running[p]['method'], running[p]['rid'],
                                                                 a.timeout * 2 + 1800))
                realign._kill_group(p)
        deadline_hit = a.deadline_h and time.time() - T0 > a.deadline_h * 3600
        if os.path.exists(stop_file) or deadline_hit:
            if pending:
                if deadline_hit:
                    log('deadline %.2f h reached: %d pending jobs recorded as skipped' % (a.deadline_h, len(pending)))
                    for j in pending:
                        why = ('not attempted: batch wall budget of %.1f h spent (jobs start in order of predicted '
                               'seconds x method priority; this one was predicted %.0f s)' % (a.deadline_h, j['pred_s']))
                        record_not_run(j['method'], j['rid'], 'skipped', why, round(j['pred_s'], 1),
                                       round(j['pred_mb']), j['threads'], a, not a.no_runtime)
                        counts['skipped'] += 1
                else:
                    log('STOP file present: %d pending jobs left unstarted (no record)' % len(pending))
                pending = []
            if not running:
                break
        # skip what the timeouts so far rule out, and refuse predicted memouts
        keep = []
        for j in pending:
            why = skip_reason(j, timeouts, a)
            if why:
                record_not_run(j['method'], j['rid'], 'skipped', why, round(j['pred_s'], 1), round(j['pred_mb']),
                               j['threads'], a, not a.no_runtime)
                counts['skipped'] += 1
                log('%-13s %-9s skipped     %s' % (j['method'], j['rid'], why[:150]))
            elif j['method'] in PAIRWISE and j['pred_mb'] > a.mem_mb:
                why = 'not run: predicted %.0f MB > cap %.0f MB at %d thread(s)' % (j['pred_mb'], a.mem_mb, j['threads'])
                record_not_run(j['method'], j['rid'], 'memout', why, round(j['pred_s'], 1), round(j['pred_mb']),
                               j['threads'], a, not a.no_runtime)
                counts['memout'] += 1
                log('%-13s %-9s memout      %s' % (j['method'], j['rid'], why))
            else:
                keep.append(j)
        pending = keep
        started = False
        busy = [pid for pid in a.busy if pid_alive(pid)]
        if pending and len(running) + len(busy) < a.jobs:
            reserved = sum(j['reserve'] for j in running.values()) + a.busy_reserve_mb * len(busy)
            avail = available_mb()
            for k, j in enumerate(pending[:20]):
                need = j['reserve']
                fits = reserved + need <= a.budget_mb and (avail is None or need <= 2000 or avail >= need + 1500)
                if fits or not (running or busy):
                    pending.pop(k)
                    prev = previous_status(j['method'], j['rid'], a.cand_root, a.panel_root)
                    if prev and not a.force and prev not in retry:      # done or started elsewhere meanwhile
                        log('%-13s %-9s already %s elsewhere; not started' % (j['method'], j['rid'], prev))
                        break
                    j['result'] = os.path.join(logdir, '%s.%s.result.json' % (j['method'], j['rid']))
                    j['log'] = os.path.join(logdir, '%s.%s.log' % (j['method'], j['rid']))
                    if os.path.exists(j['result']):
                        os.remove(j['result'])
                    cmd = job_argv(j, a) + ['--result', j['result']]
                    j['t0'] = time.time()
                    with open(j['log'], 'w') as lf:
                        p = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT, start_new_session=True,
                                             cwd=config.REPO)
                    running[p] = j
                    started = True
                    break
            if not started and time.time() - last_wait > 300:
                log('waiting for memory: %d pending, next needs %.0f MB, reserved %.0f, available %s'
                    % (len(pending), pending[0]['reserve'], reserved, avail))
                last_wait = time.time()
        if not started:
            time.sleep(0.5 if running else 3)
    log('done: %s' % ', '.join('%s %d' % kv for kv in sorted(counts.items(), key=lambda kv: str(kv[0]))))
    return 0


def available_mb():
    """Memory the OS could hand out now (free + inactive + speculative + purgeable), MB; None if unknown."""
    try:
        out = subprocess.run(['vm_stat'], capture_output=True, text=True, timeout=20).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    page, tot = 4096, 0
    for line in out.splitlines():
        if 'page size of' in line:
            page = int(line.split('page size of')[1].split()[0])
        for key in ('Pages free', 'Pages inactive', 'Pages speculative', 'Pages purgeable'):
            if line.startswith(key + ':'):
                tot += int(line.split(':')[1].strip().rstrip('.'))
    return tot * page / 2 ** 20 if tot else None


# ---------------------------------------------------------------- reports

def status_table(methods):
    rows = read_runtime()
    strata = ['hotspot_vntr', 'control_vntr_matched', 'control_vntr_correct', 'hotspot_other', 'control_nontr_sv']
    n_reg = collections.Counter(r['stratum'] for r in (json.load(open(os.path.join(config.REGIONS_DIR, rid,
                                                                                    'region.json')))
                                                       for rid in pick_regions('all')))
    out = ['method\tstatus\t' + '\t'.join(strata) + '\ttotal']
    for m in methods:
        mr = [r for r in rows if r['method'] == m]
        for st in sorted(set(r['status'] for r in mr)):
            c = collections.Counter(r['stratum'] for r in mr if r['status'] == st)
            out.append('%s\t%s\t%s\t%d' % (m, st, '\t'.join(str(c.get(s, 0)) for s in strata), sum(c.values())))
        c = collections.Counter(r['stratum'] for r in mr)
        miss = dict((s, n_reg[s] - c.get(s, 0)) for s in strata)
        if sum(miss.values()):
            out.append('%s\t(no row)\t%s\t%d' % (m, '\t'.join(str(miss[s]) for s in strata), sum(miss.values())))
    out.append('regions\t\t%s\t%d' % ('\t'.join(str(n_reg[s]) for s in strata), sum(n_reg.values())))
    return '\n'.join(out)


# ---------------------------------------------------------------- CLI

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    def common(p):
        p.add_argument('--timeout', type=float, default=DEFAULT_TIMEOUT, help='aligner wall clock (s)')
        p.add_argument('--mem-mb', type=float, default=DEFAULT_MEM_MB, help='aligner RSS cap (MB)')
        p.add_argument('--cand-root', help='candidates root (default %s)' % config.CANDIDATES_DIR)
        p.add_argument('--panel-root', help='full-panel outputs root (default %s)' % panel.PANEL_DIR)
        p.add_argument('--workdir', help='aligner scratch (default $VNTR_WORK/realign)')
        p.add_argument('--no-runtime', action='store_true', help='do not write the runtime table')
        p.add_argument('--panel-engine', default='native', choices=['auto', 'vg', 'native'],
                       help='msa_graph engine for the full-panel graph (vg construct -M is ~50x slower on '
                            'hundreds of rows; default native, which msa_graph tests against vg)')
    p = sub.add_parser('run', help='a batch of regions x methods')
    p.add_argument('--methods', default=','.join(METHODS))
    p.add_argument('--regions', default='all')
    p.add_argument('--stratum')
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--budget-mb', type=float, default=DEFAULT_BUDGET_MB)
    p.add_argument('--threads-pairwise', type=int, default=3)
    p.add_argument('--threads-progressive', type=int, default=2)
    p.add_argument('--skip-after', type=int, default=2,
                   help='skip a pairwise job once its method timed out on this many regions with fewer pair '
                        'cells (0 = never)')
    p.add_argument('--skip-factor', type=float, default=3.0,
                   help='skip a job predicted to need more than this x --timeout (0 = never)')
    p.add_argument('--no-skip-methods', default='mafft_fftns2',
                   help='methods never skipped by --skip-after or --skip-factor (default mafft_fftns2: every '
                        'region is attempted)')
    p.add_argument('--priority', default='mafft_fftns2=0.25,mafft_linsi=1,mafft_fftnsi=1,mafft_einsi=1.5,'
                                         'mafft_ginsi=1.5',
                   help='jobs start in order of predicted seconds x this per-method weight')
    p.add_argument('--prefer-regions', default='',
                   help='comma list of regions whose jobs (of --prefer-methods) jump the queue: their rank is '
                        'multiplied by --prefer-factor (e.g. the regions of a first look)')
    p.add_argument('--prefer-methods', default='', help='methods --prefer-regions applies to (default all)')
    p.add_argument('--prefer-factor', type=float, default=0.02)
    p.add_argument('--busy-pids', default='',
                   help='comma list of pids of jobs started by an earlier scheduler: each occupies a job slot '
                        '(and --busy-reserve-mb) while alive')
    p.add_argument('--busy-reserve-mb', type=float, default=7000)
    p.add_argument('--exclude', default='', help='comma list of METHOD:REGION jobs not to schedule')
    p.add_argument('--deadline-h', type=float, default=0,
                   help='start no job after this many hours; record the rest as skipped (0 = none)')
    p.add_argument('--force', action='store_true')
    p.add_argument('--retry', help='re-run jobs whose last status is one of these (e.g. timeout,memout,skipped)')
    common(p)
    p = sub.add_parser('one', help='one job in this process')
    p.add_argument('method')
    p.add_argument('region')
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--result', help='also write the runtime row here (JSON)')
    common(p)
    p = sub.add_parser('predict', help='predicted seconds and memory per region and method')
    p.add_argument('--methods', default=','.join(METHODS))
    p.add_argument('--regions', default='all')
    p.add_argument('--mem-mb', type=float, default=DEFAULT_MEM_MB)
    p.add_argument('--threads-pairwise', type=int, default=3)
    p.add_argument('--threads-progressive', type=int, default=2)
    p = sub.add_parser('status', help='counts per method, status and stratum')
    p.add_argument('--methods', default=','.join(METHODS))
    p = sub.add_parser('mark-skipped', help="record every job without a result as 'skipped'")
    p.add_argument('--methods', default=','.join(METHODS))
    p.add_argument('--regions', default='all')
    p.add_argument('--reason', required=True)
    common(p)
    a = ap.parse_args(argv)

    if a.cmd == 'predict':
        w = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
        w.writerow(['region_id', 'stratum', 'method', 'n_distinct', 'max_len', 'len2', 'total_bp', 'pair_cells',
                    'threads', 'predicted_s', 'predicted_mb'])
        rows = []
        for rid in pick_regions(a.regions):
            r = region_info(rid)
            for m in a.methods.split(','):
                t, pm = plan_threads(m, r, a)
                rows.append([rid, r['stratum'], m, r['n_distinct'], r['max_len'], r['len2'], r['total_bp'],
                             r['pair_cells'], t, round(predict_s(m, r, t), 1), round(pm)])
        rows.sort(key=lambda x: (x[2], x[9]))
        w.writerows(rows)
        return 0
    if a.cmd == 'status':
        print(status_table(a.methods.split(',')))
        return 0
    if a.cmd == 'mark-skipped':
        n = 0
        for rid in pick_regions(a.regions):
            for m in a.methods.split(','):
                if previous_status(m, rid, a.cand_root, a.panel_root) is None:
                    record_not_run(m, rid, 'skipped', 'not attempted: ' + a.reason, a=a,
                                   record_runtime=not a.no_runtime)
                    n += 1
        log('%d jobs recorded as skipped' % n)
        return 0
    if a.cmd == 'one':
        row = run_job(a.method, a.region, threads=a.threads, timeout=a.timeout, mem_mb=a.mem_mb,
                      cand_root=a.cand_root, panel_root=a.panel_root, workroot=a.workdir,
                      record_runtime=not a.no_runtime, panel_engine=a.panel_engine)
        if a.result:
            with open(a.result + '.tmp', 'w') as f:
                json.dump(row, f)
            os.replace(a.result + '.tmp', a.result)
        log('%s %s: %s in %.1f s %s' % (a.method, a.region, row['status'], row['seconds'], row.get('note', '')))
        return 0 if row['status'] == 'ok' else (4 if row['status'] in ('timeout', 'memout') else 2)
    return run_batch(a)


if __name__ == '__main__':
    sys.exit(main())
