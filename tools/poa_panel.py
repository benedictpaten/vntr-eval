#!/usr/bin/env python3
"""poa_panel.py -- the full-panel arm of the POA realigners (poa_abpoa, poa_spoa).

For every region and method this aligns the deduplicated FULL panel (work/panel/union/<id>.fa,
written by `tools/panel.py union`: every distinct sequence of hprc.fa.gz and hap32.fa), keeps the
full MSA and its full-panel graph, and projects the MSA onto the 34 hap32 rows:

    tools/panel.py union            work/panel/union/<id>.fa + <id>.map.tsv (made here when missing)
    realign.align_fasta(METHOD)     the realigner's MSA-only mode: same dedup, N masking, order rule
                                    (longest first), timeout and RSS cap as the hap32 arm
    -> work/panel/<METHOD>/<id>.msa.fa.gz   the full MSA, rows named by distinct id (u0001..; the
                                    map TSV lists who carries each)
    -> work/panel/<METHOD>/<id>.gfa         its graph, one P line per hprc.fa.gz record (panel.py
                                    panel-graph: distinct rows expanded back to member names;
                                    every path's spelling checked)
    -> candidates/<METHOD>__all/<id>.msa.fa the hap32 rows (panel.py project; all-gap columns dropped)
    -> candidates/<METHOD>__all/<id>.gfa    msa_graph.py, one P line per hap32.fa name; every path is
                                    re-read from the written file and asserted to spell hap32.fa
    -> candidates/<METHOD>__all/<id>.realign.json  how it was made (also written for skips,
                                    timeouts and failures, so a re-run skips them)
    -> results/realign_runtime.all.poa.tsv  one row per (method, region); the latest run wins

Usage (from the repository root):

    python3 tools/poa_panel.py run [--methods poa_abpoa,poa_spoa] [--regions ID,...|all]
        [--stratum S] [--jobs 2] [--timeout 1800] [--mem-mb 12000] [--budget-mb 14000]
        [--force | --retry-failed] [--order longest|random:SEED|...] [--name NAME]
        [--cand-root DIR --panel-root DIR --no-runtime]      (pilot variants go to scratch)
    python3 tools/poa_panel.py one METHOD REGION [same options]   # one job in this process
    python3 tools/poa_panel.py predict [--regions ...]            # predicted memory table

Scheduling: jobs run as separate processes (the graph building is Python), at most --jobs at a
time, regions smallest first (total distinct bp) with the two methods interleaved. A job starts
only when the predicted memory of the running jobs plus its own stays within --budget-mb and the
OS reports that much available; when nothing is running the next job starts regardless (its RSS
cap still applies). A later, smaller job may overtake one that waits for memory.

Memory prediction (MB), from the two longest sequences len1 >= len2 of the union:
    spoa   50 + 80 B x len1^2        (its preallocation for 4 x len1 graph nodes; tools/POA.md. The
                                     full-panel graphs stay far below 4 x len1 nodes: ~1.2 x len1)
    abPOA  50 + c x len1 x len2,     c = 11 B (int16 DP) when len1 <= 13000, else 22 B (int32:
                                     abPOA widens its cells once the graph passes ~16k nodes)
           Longest first, the second sequence is aligned to a graph of >= len1 nodes. Full panel,
           /usr/bin/time -l: L014297 692 MB (predicted 745), L005990 1,403 MB (1,490).
This replaces realign_poa.predict_mb (the hap32 arm's 22-24 B x len1^2, which refuses regions
whose longest sequence is a lone outlier) inside a job process only. A job predicted over
--mem-mb is not run (status memout, "predicted" in the note); a running aligner is killed when
its sampled RSS passes --mem-mb (status memout) or its wall clock passes --timeout (status
timeout). The scheduler reserves 1.25 x the prediction (at most --mem-mb) per running job.
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
import realign_poa  # noqa: E402

BASE_METHODS = ('poa_abpoa', 'poa_spoa')
RUNTIME_TSV = os.path.join(config.RESULTS_DIR, 'realign_runtime.all.poa.tsv')
RUNTIME_COLS = ['region_id', 'method', 'panel', 'n_distinct', 'seconds', 'status',
                'stratum', 'span_bp', 'n_panel_records', 'n_hap32_only', 'max_len', 'len2', 'total_bp',
                'order', 'predicted_mb', 'peak_rss_mb', 'mem_cap_mb', 'timeout_s', 'align_s', 'panel_graph_s',
                'project_s', 'graph_s', 'columns_full', 'nodes_full', 'edges_full', 'columns_proj',
                'nodes_proj', 'nodes_per_kb_proj', 'frac_nodes_1bp_proj', 'masked_runs', 'masked_bp',
                'finished', 'note']
DEFAULT_TIMEOUT = 1800
DEFAULT_MEM_MB = 12000
DEFAULT_BUDGET_MB = 14000
T0 = time.time()


def log(msg):
    sys.stderr.write('[poa_panel %7.1fs] %s\n' % (time.time() - T0, msg))
    sys.stderr.flush()


def tool_of(method):
    return 'spoa' if method.startswith('poa_spoa') else 'abpoa'


def predict_mb(tool, lens):
    """Predicted peak RSS (MB) of the aligner on sequences of these lengths (module docstring)."""
    if not lens:
        return 50.0
    s = sorted(lens, reverse=True)
    l1, l2 = s[0], (s[1] if len(s) > 1 else 0)
    if tool == 'spoa':
        return round(50 + 80.0 * l1 * l1 / 2 ** 20, 1)
    c = 11.0 if l1 <= 13000 else 22.0
    return round(50 + c * l1 * l2 / 2 ** 20, 1)


def _aligner_predict(tool, lens, args=()):
    """Stand-in for realign_poa.predict_mb inside a job (same signature)."""
    return predict_mb(tool, lens)


# ---------------------------------------------------------------- paths

def union_files(rid):
    fa = os.path.join(panel.UNION_DIR, rid + '.fa')
    mp = os.path.join(panel.UNION_DIR, rid + '.map.tsv')
    return fa, mp


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


def previous_status(method, rid, cand_root=None, panel_root=None):
    o = out_paths(method, rid, cand_root, panel_root)
    try:
        j = json.load(open(o['json']))
    except (OSError, ValueError):
        return None
    st = j.get('status')
    if st == 'ok' and not all(os.path.exists(o[k]) for k in ('msa', 'gfa', 'full_msa', 'full_gfa')):
        return None
    return st


def region_info(rid):
    """Sizes of a region's union (for ordering and prediction)."""
    rd = realign.region_dir(rid)
    fa, mp = ensure_union(rd, rid)
    lens, n_only, n_rec = [], 0, 0
    for r in panel.read_map(mp):
        lens.append(int(r['length']))
        if int(r['weight']) == 0:
            n_only += 1
        n_rec += int(r['weight'])
    stratum, span = realign.region_meta(rd)
    s = sorted(lens, reverse=True)
    return {'region_id': rid, 'rd': rd, 'stratum': stratum, 'span_bp': span, 'n_distinct': len(lens),
            'n_hap32_only': n_only, 'n_panel_records': n_rec, 'max_len': s[0], 'len2': s[1] if len(s) > 1 else 0,
            'total_bp': sum(lens), 'lens': lens}


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
        rows.sort(key=lambda r: (r['method'], r['region_id']))
        tmp = path + '.tmp%d' % os.getpid()
        with open(tmp, 'w', newline='') as f:
            w = csv.DictWriter(f, RUNTIME_COLS, delimiter='\t', extrasaction='ignore', lineterminator='\n')
            w.writeheader()
            for r in rows:
                w.writerow(dict((c, r.get(c, '')) for c in RUNTIME_COLS))
        os.replace(tmp, path)
        fcntl.flock(lk, fcntl.LOCK_UN)


# ---------------------------------------------------------------- one job

def register(method, base, order=None):
    """Make `method` known to realign's registry (a variant of base with another input order)."""
    ms = realign.all_methods()
    if method in ms and order is None:
        return ms[method]
    name, spec = realign_poa.variant(base, order=order, name=method)
    ms[name] = realign._as_method(name, spec)
    return ms[name]


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


def run_job(method, rid, base=None, order=None, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB,
            cand_root=None, panel_root=None, workroot=None, record_runtime=True, panel_engine='native',
            keep_uncompressed=False):
    """Align region rid's full-panel union with method; write every output; return the runtime row."""
    base = base or method
    m = register(method, base, order)
    realign_poa.predict_mb = _aligner_predict          # the full-panel memory model (docstring)
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
    tool = tool_of(base)
    pred = predict_mb(tool, info_r['lens'])
    info = collections.OrderedDict([
        ('method', method), ('base_method', base), ('panel', 'all'), ('region_id', rid),
        ('stratum', info_r['stratum']), ('span_bp', info_r['span_bp']), ('status', None),
        ('description', m.description), ('order', m.params.get('order')),
        ('union_fa', realign._rel(fa)), ('union_map', realign._rel(mp)),
        ('n_distinct', info_r['n_distinct']), ('n_hap32_only', info_r['n_hap32_only']),
        ('n_panel_records', info_r['n_panel_records']), ('max_len', info_r['max_len']),
        ('len2', info_r['len2']), ('total_bp', info_r['total_bp']), ('predicted_mb', pred),
        ('mem_cap_mb', mem_mb), ('timeout_s', timeout)])
    row = {}
    times = {}
    wd = tempfile.mkdtemp(prefix='poa_panel.%s.%s.' % (method, rid), dir=workroot)
    try:
        if mem_mb and pred > mem_mb:
            info['status'] = 'memout'
            info['message'] = 'not run: predicted %.0f MB > cap %.0f MB (len1 %d, len2 %d)' % (
                pred, mem_mb, info_r['max_len'], info_r['len2'])
        else:
            full_msa = os.path.join(wd, 'full.msa.fa')
            al = realign.align_fasta(m, fa, full_msa, threads=1, timeout=timeout, mem_mb=mem_mb, workdir=wd,
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
                    if keep_uncompressed:
                        shutil.copy(full_msa, o['full_msa'][:-3])
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
    info['poa_panel_version'] = 1
    tmpj = o['json'] + '.tmp%d' % os.getpid()
    with open(tmpj, 'w') as f:
        json.dump(info, f, indent=1)
    os.replace(tmpj, o['json'])
    al = info.get('align') or {}
    aln = al.get('aligner') or {}
    pg = info.get('panel_graph') or {}
    g = info.get('graph') or {}
    pr = info.get('projection') or {}
    row = {'region_id': rid, 'method': method, 'panel': 'all', 'n_distinct': info_r['n_distinct'],
           'seconds': info['seconds'], 'status': info['status'], 'stratum': info_r['stratum'],
           'span_bp': info_r['span_bp'], 'n_panel_records': info_r['n_panel_records'],
           'n_hap32_only': info_r['n_hap32_only'], 'max_len': info_r['max_len'], 'len2': info_r['len2'],
           'total_bp': info_r['total_bp'], 'order': info.get('order'), 'predicted_mb': pred,
           'peak_rss_mb': aln.get('peak_rss_mb', ''), 'mem_cap_mb': mem_mb, 'timeout_s': timeout,
           'align_s': times.get('align_s', ''), 'panel_graph_s': times.get('panel_graph_s', ''),
           'project_s': times.get('project_s', ''), 'graph_s': times.get('graph_s', ''),
           'columns_full': al.get('columns', ''), 'nodes_full': pg.get('nodes', ''), 'edges_full': pg.get('edges', ''),
           'columns_proj': pr.get('columns_out', ''), 'nodes_proj': g.get('nodes', ''),
           'nodes_per_kb_proj': g.get('nodes_per_kb', ''), 'frac_nodes_1bp_proj': g.get('frac_nodes_1bp', ''),
           'masked_runs': al.get('masked_runs', ''), 'masked_bp': al.get('masked_bp', ''),
           'finished': info['finished'],
           'note': (info.get('message') or '')[:200].replace('\t', ' ').replace('\n', ' ')}
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
            raise SystemExit('poa_panel: no region package for %s' % ', '.join(missing))
        rows = [r for r in rows if r['region_id'] in set(want)]
    if stratum:
        keep = set(stratum.split(','))
        rows = [r for r in rows if r.get('stratum') in keep]
    return [r['region_id'] for r in rows]


def job_argv(method, rid, a):
    cmd = [sys.executable, os.path.abspath(__file__), 'one', method, rid, '--timeout', str(a.timeout),
           '--mem-mb', str(a.mem_mb), '--panel-engine', a.panel_engine]
    if a.base:
        cmd += ['--base', a.base]
    if a.order:
        cmd += ['--order', a.order]
    for opt in ('cand_root', 'panel_root', 'workdir'):
        v = getattr(a, opt)
        if v:
            cmd += ['--' + opt.replace('_', '-'), v]
    if a.no_runtime:
        cmd.append('--no-runtime')
    return cmd


def run_batch(a):
    methods = [x for x in a.methods.split(',') if x]
    if a.name:
        if len(methods) != 1:
            raise SystemExit('--name takes one method')
        a.base, names = methods[0], [a.name]
    elif a.order:
        a.base = methods[0]
        if len(methods) != 1:
            raise SystemExit('--order takes one method (the variant gets its own name)')
        names = [realign_poa.variant(methods[0], order=a.order)[0]]
    else:
        a.base, names = None, methods
    rids = pick_regions(a.regions, a.stratum)
    infos = {}
    for rid in rids:
        infos[rid] = region_info(rid)
    rids.sort(key=lambda r: (infos[r]['total_bp'], infos[r]['max_len'], r))
    jobs = []
    for rid in rids:
        for mname in names:
            base = a.base or mname
            prev = previous_status(mname, rid, a.cand_root, a.panel_root)
            if prev and not a.force and not (a.retry_failed and prev != 'ok'):
                continue
            pred = predict_mb(tool_of(base), infos[rid]['lens'])
            jobs.append({'method': mname, 'rid': rid, 'pred': pred, 'reserve': min(1.25 * pred, a.mem_mb)})
    log('%d jobs (%d regions x %s); budget %d MB, cap %d MB, %d at a time, timeout %d s'
        % (len(jobs), len(rids), ','.join(names), a.budget_mb, a.mem_mb, a.jobs, a.timeout))
    logdir = os.path.join(a.workdir or os.path.join(config.WORK_DIR, 'realign'), 'poa_panel_logs')
    os.makedirs(logdir, exist_ok=True)
    running = {}
    counts = collections.Counter()
    pending = list(jobs)
    last_wait = 0
    while pending or running:
        # reap
        for p in list(running):
            if p.poll() is not None:
                j = running.pop(p)
                row = {}
                try:
                    row = json.load(open(j['result']))
                except (OSError, ValueError):
                    row = {'status': 'crash', 'note': 'no result; exit %s; log %s' % (p.returncode, j['log'])}
                    if not a.no_runtime:
                        update_runtime({'region_id': j['rid'], 'method': j['method'], 'panel': 'all',
                                        'status': 'crash', 'n_distinct': infos[j['rid']]['n_distinct'],
                                        'seconds': round(time.time() - j['t0'], 1), 'note': row['note']})
                counts[row.get('status')] += 1
                log('%-22s %-9s %-11s %7ss align %7ss rss %7s MB  nodes full %-6s proj %-6s %s' % (
                    j['method'], j['rid'], row.get('status'), row.get('seconds', ''), row.get('align_s', ''),
                    row.get('peak_rss_mb', ''), row.get('nodes_full', ''), row.get('nodes_proj', ''),
                    (row.get('note') or '')[:120]))
            elif time.time() - running[p]['t0'] > a.timeout * 2 + 1800:
                log('killing %s %s: job wall clock over %d s' % (running[p]['method'], running[p]['rid'],
                                                                 a.timeout * 2 + 1800))
                realign._kill_group(p)
        # start
        started = False
        if pending and len(running) < a.jobs:
            reserved = sum(j['reserve'] for j in running.values())
            avail = realign_poa.available_mb()
            for k, j in enumerate(pending):
                need = j['reserve']
                fits = reserved + need <= a.budget_mb and (avail is None or need <= 2000 or
                                                            avail >= need + 1500)
                if j['pred'] > a.mem_mb:       # refused by prediction: costs nothing
                    fits = True
                if fits or not running:
                    pending.pop(k)
                    j['result'] = os.path.join(logdir, '%s.%s.result.json' % (j['method'], j['rid']))
                    j['log'] = os.path.join(logdir, '%s.%s.log' % (j['method'], j['rid']))
                    if os.path.exists(j['result']):
                        os.remove(j['result'])
                    cmd = job_argv(j['method'], j['rid'], a) + ['--result', j['result']]
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


# ---------------------------------------------------------------- CLI

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    def common(p):
        p.add_argument('--timeout', type=float, default=DEFAULT_TIMEOUT, help='aligner wall clock (s)')
        p.add_argument('--mem-mb', type=float, default=DEFAULT_MEM_MB, help='aligner RSS cap (MB)')
        p.add_argument('--order', help='input order for a variant (given, longest, shortest, guide, random:SEED)')
        p.add_argument('--cand-root', help='candidates root (default %s)' % config.CANDIDATES_DIR)
        p.add_argument('--panel-root', help='full-panel outputs root (default %s)' % panel.PANEL_DIR)
        p.add_argument('--workdir', help='aligner scratch (default $VNTR_WORK/realign)')
        p.add_argument('--no-runtime', action='store_true', help='do not write the runtime table')
        p.add_argument('--panel-engine', default='native', choices=['auto', 'vg', 'native'],
                       help='msa_graph engine for the full-panel graph (vg construct -M is ~50x slower on '
                            'hundreds of rows; default native)')
    p = sub.add_parser('run', help='a batch of regions x methods')
    p.add_argument('--methods', default=','.join(BASE_METHODS))
    p.add_argument('--regions', default='all')
    p.add_argument('--stratum')
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--budget-mb', type=float, default=DEFAULT_BUDGET_MB,
                   help='predicted MB of all running jobs together (default %d)' % DEFAULT_BUDGET_MB)
    p.add_argument('--force', action='store_true')
    p.add_argument('--retry-failed', action='store_true')
    p.add_argument('--name', help='output method name for a variant (one method)')
    common(p)
    p = sub.add_parser('one', help='one job in this process')
    p.add_argument('method')
    p.add_argument('region')
    p.add_argument('--base', help='the POA method a variant is built from')
    p.add_argument('--result', help='also write the runtime row here (JSON)')
    common(p)
    p = sub.add_parser('predict', help='predicted memory per region and method')
    p.add_argument('--regions', default='all')
    p.add_argument('--mem-mb', type=float, default=DEFAULT_MEM_MB)
    a = ap.parse_args(argv)

    if a.cmd == 'predict':
        w = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
        w.writerow(['region_id', 'stratum', 'n_distinct', 'max_len', 'len2', 'total_bp', 'abpoa_mb', 'spoa_mb',
                    'abpoa_fits', 'spoa_fits'])
        for rid in pick_regions(a.regions):
            r = region_info(rid)
            pa, ps = predict_mb('abpoa', r['lens']), predict_mb('spoa', r['lens'])
            w.writerow([rid, r['stratum'], r['n_distinct'], r['max_len'], r['len2'], r['total_bp'], pa, ps,
                        int(pa <= a.mem_mb), int(ps <= a.mem_mb)])
        return 0
    if a.cmd == 'one':
        base = a.base or a.method
        if base not in BASE_METHODS:
            raise SystemExit('poa_panel: base method must be one of %s' % ', '.join(BASE_METHODS))
        row = run_job(a.method, a.region, base=base, order=a.order, timeout=a.timeout, mem_mb=a.mem_mb,
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
