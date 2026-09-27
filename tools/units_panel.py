#!/usr/bin/env python3
"""units_panel.py -- the full-panel arm of M4, the repeat-unit-aware realigner (unit_aware).

For every region this aligns the deduplicated FULL panel (work/panel/union/<id>.fa, written by
`tools/panel.py union`: every distinct sequence of hprc.fa.gz and hap32.fa, CHM13's first) with
realign_units.align_units_fasta -- the same method, parameters, unit/fallback decision, guard and
fallback (mafft L-INS-i, or FFT-NS-i above 100 distinct sequences) as the hap32 arm -- keeps the
full MSA and its full-panel graph, and projects the MSA onto the 34 hap32 rows:

    tools/panel.py union            work/panel/union/<id>.fa + <id>.map.tsv (made here when missing)
    realign_units.align_units_fasta the MSA-only mode on the union (rows u0001..; the repeat model
                                    comes from region.json and CHM13's hap32.fa record)
    -> work/panel/unit_aware/<id>.msa.fa.gz  the full MSA, rows named by distinct id (the map TSV
                                    lists who carries each)
    -> work/panel/unit_aware/<id>.gfa        its graph, one P line per hprc.fa.gz record (panel.py
                                    panel-graph: distinct rows expanded back to member names;
                                    every path's spelling re-checked here)
    -> candidates/unit_aware__all/<id>.msa.fa    the hap32 rows (panel.py project; all-gap columns
                                    dropped; every row asserted to spell hap32.fa)
    -> candidates/unit_aware__all/<id>.gfa       msa_graph.py (default engine, as the hap32 arm), one
                                    P line per hap32.fa name; re-read and asserted to spell hap32.fa
    -> candidates/unit_aware__all/<id>.realign.json  how it was made (also for timeouts, memouts and
                                    failures, so a re-run skips them unless --force/--retry-failed)
    -> results/realign_runtime.all.units.tsv  one row per region (region_id, method, panel=all,
                                    n_distinct, seconds, status, mode, ...); the latest run wins

Usage (from the repository root):

    python3 tools/units_panel.py run [--regions ID,...|all] [--stratum S] [--jobs 2] [--threads 3]
        [--timeout 1800] [--mem-mb 12000] [--force | --retry-failed]
        [--cand-root DIR --panel-root DIR --no-runtime]      (trial runs go to scratch)
    python3 tools/units_panel.py one REGION [--threads 3] [--timeout 1800] ...   # one job, this process
    python3 tools/units_panel.py table [--log BATCH.log]   # rewrite the runtime table from the JSONs

Scheduling: every region is a separate process (its own session), at most --jobs at a time,
regions smallest first (total distinct bp), so one slow region does not starve the rest. The
parent samples the summed RSS of the job's whole process tree (mafft runs in sessions of its
own) and kills the tree when the job's wall clock passes --timeout (status timeout) or the RSS
passes --mem-mb (status memout); the region then gets a realign.json and a runtime row with that
status. A second job starts only when the OS reports enough available memory for it
(predicted MB + 1.5 GB); when nothing runs, the next job starts regardless.
Inside a job, every mafft call (flank pieces, the fallback, the adequacy guard's fallback) is
clamped to what is left of the region's budget less a margin for the graphs (min(300 s, 15%)),
and to --mem-mb less this process's peak RSS; so a guard fallback that cannot finish in time or
memory stops on its own and the unit MSA is kept (guard fallback_status timeout / memout),
instead of the region being killed from outside. The method's own caps stay: unit-level
refinement and the base-level polish stop after 300 s each.
"""
import argparse
import collections
import csv
import datetime
import fcntl
import gzip
import json
import os
import resource
import shutil
import signal
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
import realign_units  # noqa: E402

METHOD = realign_units.METHOD            # unit_aware
RUNTIME_TSV = os.path.join(config.RESULTS_DIR, 'realign_runtime.all.units.tsv')
RUNTIME_COLS = ['region_id', 'method', 'panel', 'n_distinct', 'seconds', 'status', 'mode', 'fallback_method',
                'reason', 'stratum', 'span_bp', 'n_panel_records', 'n_hap32_only', 'max_len', 'len2', 'total_bp',
                'predicted_mb', 'peak_rss_mb', 'mem_cap_mb', 'timeout_s', 'threads', 'align_s', 'panel_graph_s',
                'project_s', 'graph_s', 'symbols', 'unit_strings_distinct', 'units_per_seq_max', 'guide',
                'refine_timed_out', 'polish_timed_out', 'guard_ratio', 'guard_fallback_method', 'guard_fallback_status',
                'guard_fallback_ratio', 'guard_decision',
                'phase_one_start', 'columns_full', 'nodes_full', 'edges_full', 'columns_proj', 'nodes_proj',
                'nodes_per_kb_proj', 'frac_nodes_1bp_proj', 'masked_runs', 'masked_bp', 'host_load', 'finished',
                'note']
DEFAULT_TIMEOUT = 1800
DEFAULT_MEM_MB = 12000
DEFAULT_THREADS = 3
T0 = time.time()


def log(msg):
    sys.stderr.write('[units_panel %7.1fs] %s\n' % (time.time() - T0, msg))
    sys.stderr.flush()


# ---------------------------------------------------------------- paths and sizes

def union_files(rid):
    return (os.path.join(panel.UNION_DIR, rid + '.fa'), os.path.join(panel.UNION_DIR, rid + '.map.tsv'))


def ensure_union(rd, rid):
    fa, mp = union_files(rid)
    if not (os.path.exists(fa) and os.path.exists(mp)):
        panel.union_to_files(rd, fa, mp)
    return fa, mp


def out_paths(rid, cand_root=None, panel_root=None, method=METHOD):
    cd = os.path.join(cand_root or config.CANDIDATES_DIR, method + '__all')
    pd = os.path.join(panel_root or panel.PANEL_DIR, method)
    return {'cand_dir': cd, 'panel_dir': pd,
            'msa': os.path.join(cd, rid + '.msa.fa'), 'gfa': os.path.join(cd, rid + '.gfa'),
            'json': os.path.join(cd, rid + '.realign.json'),
            'full_msa': os.path.join(pd, rid + '.msa.fa.gz'), 'full_gfa': os.path.join(pd, rid + '.gfa')}


def previous_status(rid, cand_root=None, panel_root=None):
    o = out_paths(rid, cand_root, panel_root)
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
    rd = realign.region_dir(rid)
    fa, mp = ensure_union(rd, rid)
    lens, n_only, n_rec = [], 0, 0
    for r in panel.read_map(mp):
        lens.append(int(r['length']))
        n_only += int(r['weight']) == 0
        n_rec += int(r['weight'])
    stratum, span = realign.region_meta(rd)
    s = sorted(lens, reverse=True)
    return {'region_id': rid, 'rd': rd, 'stratum': stratum, 'span_bp': span, 'n_distinct': len(lens),
            'n_hap32_only': n_only, 'n_panel_records': n_rec, 'max_len': s[0], 'len2': s[1] if len(s) > 1 else 0,
            'total_bp': sum(lens), 'predicted_mb': predict_mb(len(lens), s[0])}


def predict_mb(n_distinct, max_len, cap=DEFAULT_MEM_MB):
    """Rough peak RSS (MB) of one job's process tree, for scheduling only. Measured on the 149
    full-panel runs: the unit path peaks at 290-610 B per (distinct sequence x longest sequence)
    cell (the base-level polish copies the n x columns matrix per row; L012184 6.6 GB,
    L015347 7.5 GB), so 300 MB + 450 B per cell, capped at the RSS cap. mafft fallbacks are not
    predicted: L-INS-i on 55 x 11 kb took 7.0 GB (L015743), and the guard's FFT-NS-i on 432 x
    10 kb 9.7 GB (L014160)."""
    return round(min(cap, 300 + 450.0 * n_distinct * max_len / 2 ** 20), 0)


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


# ---------------------------------------------------------------- checks

def check_paths(gfa, fasta_recs):
    """Re-read gfa; assert every named sequence has a path spelling it and there are no others."""
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


# ---------------------------------------------------------------- one job

def summary_row(rid, info_r, info, times, status, seconds, threads, timeout, mem_mb, peak=None, note=None):
    al = info.get('align') or {}
    pg = info.get('panel_graph') or {}
    g = info.get('graph') or {}
    pr = info.get('projection') or {}
    u = al.get('units') or {}
    um = u.get('unit_msa') or {}
    gd = al.get('guard') or {}
    pol = al.get('polish') or {}
    ph = al.get('phase_check') or {}
    fb = al.get('fallback') or {}
    fbm = fb.get('method') or gd.get('fallback_method') if al.get('mode') == 'fallback' else ''
    ups = u.get('units_per_seq') or []
    if note is None:
        note = info.get('message') or ''
    return {'region_id': rid, 'method': METHOD, 'panel': 'all', 'n_distinct': info_r['n_distinct'],
            'seconds': seconds, 'status': status, 'mode': al.get('mode') or '', 'fallback_method': fbm or '',
            'reason': (al.get('reason') or '')[:160].replace('\t', ' ').replace('\n', ' '),
            'stratum': info_r['stratum'], 'span_bp': info_r['span_bp'], 'n_panel_records': info_r['n_panel_records'],
            'n_hap32_only': info_r['n_hap32_only'], 'max_len': info_r['max_len'], 'len2': info_r['len2'],
            'total_bp': info_r['total_bp'], 'predicted_mb': info_r['predicted_mb'],
            'peak_rss_mb': '' if peak is None else round(peak, 0), 'mem_cap_mb': mem_mb, 'timeout_s': timeout,
            'threads': threads, 'align_s': times.get('align_s', ''), 'panel_graph_s': times.get('panel_graph_s', ''),
            'project_s': times.get('project_s', ''), 'graph_s': times.get('graph_s', ''),
            'symbols': u.get('symbols', ''), 'unit_strings_distinct': u.get('unit_strings_distinct', ''),
            'units_per_seq_max': ups[-1] if ups else '', 'guide': um.get('guide', ''),
            'refine_timed_out': um.get('refine_timed_out', ''), 'polish_timed_out': pol.get('timed_out', ''),
            'guard_ratio': gd.get('unit_sp_over_opt', ''), 'guard_fallback_method': gd.get('fallback_method', ''),
            'guard_fallback_status': gd.get('fallback_status', ''),
            'guard_fallback_ratio': gd.get('fallback_sp_over_opt', ''),
            'guard_decision': gd.get('decision', ''),
            'phase_one_start': ph.get('units_in_one_start_columns', ''),
            'columns_full': al.get('columns', ''), 'nodes_full': pg.get('nodes', ''), 'edges_full': pg.get('edges', ''),
            'columns_proj': pr.get('columns_out', ''), 'nodes_proj': g.get('nodes', ''),
            'nodes_per_kb_proj': g.get('nodes_per_kb', ''), 'frac_nodes_1bp_proj': g.get('frac_nodes_1bp', ''),
            'masked_runs': al.get('masked_runs', ''), 'masked_bp': al.get('masked_bp', ''),
            'host_load': round(os.getloadavg()[0], 1),
            'finished': datetime.datetime.now().isoformat(timespec='seconds'),
            'note': str(note)[:200].replace('\t', ' ').replace('\n', ' ')}


def write_json(path, info):
    tmp = path + '.tmp%d' % os.getpid()
    with open(tmp, 'w') as f:
        json.dump(info, f, indent=1)
    os.replace(tmp, path)


def run_job(rid, threads=DEFAULT_THREADS, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB, cand_root=None,
            panel_root=None, workroot=None, record_runtime=True, panel_engine='native', params=None):
    """Align region rid's full-panel union with unit_aware; write every output; return the runtime row."""
    info_r = region_info(rid)
    rd = info_r['rd']
    fa, mp = union_files(rid)
    hap_fa = os.path.join(rd, 'hap32.fa')
    rjp = os.path.join(rd, 'region.json')
    o = out_paths(rid, cand_root, panel_root)
    os.makedirs(o['cand_dir'], exist_ok=True)
    os.makedirs(o['panel_dir'], exist_ok=True)
    workroot = workroot or os.path.join(config.WORK_DIR, 'realign_units')
    os.makedirs(workroot, exist_ok=True)
    for k in ('msa', 'gfa', 'full_msa', 'full_gfa'):      # never leave an older result beside a new attempt
        if os.path.exists(o[k]):
            os.remove(o[k])
    P = params or realign_units.Params()
    P.timeout = int(timeout)
    P.mem_mb = int(mem_mb)
    t0 = time.time()
    # Every mafft call of the method (flank pieces, the fallback, the guard's fallback) goes
    # through realign.align_fasta; clamp its timeout to what is left of this region's budget,
    # less a margin for building the graphs. A guard fallback that cannot finish in time then
    # times out on its own and the unit MSA is kept (the guard records fallback_status
    # 'timeout'), instead of the whole region being killed from outside. Likewise for memory.
    deadline = t0 + timeout - min(300.0, 0.15 * timeout)
    orig_align_fasta = realign.align_fasta

    def clamped_align_fasta(*args, **kw):
        left = deadline - time.time()
        kw['timeout'] = max(30, int(min(kw.get('timeout') or timeout, left)))
        # the parent caps this job's whole process tree at mem_mb: leave room for this process
        # (its peak so far) so an oversized mafft is stopped here, gracefully, not the region
        own = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (2.0 ** 20 if sys.platform == 'darwin' else 1024.0)
        kw['mem_mb'] = max(500, int(min(kw.get('mem_mb') or mem_mb, mem_mb - own - 500)))
        return orig_align_fasta(*args, **kw)
    realign.align_fasta = clamped_align_fasta
    info = collections.OrderedDict([
        ('method', METHOD), ('panel', 'all'), ('region_id', rid), ('stratum', info_r['stratum']),
        ('span_bp', info_r['span_bp']), ('status', None), ('union_fa', realign._rel(fa)),
        ('union_map', realign._rel(mp)), ('n_distinct', info_r['n_distinct']),
        ('n_hap32_only', info_r['n_hap32_only']), ('n_panel_records', info_r['n_panel_records']),
        ('max_len', info_r['max_len']), ('len2', info_r['len2']), ('total_bp', info_r['total_bp']),
        ('timeout_s', timeout), ('mafft_deadline_s', round(deadline - t0, 1)), ('mem_cap_mb', mem_mb),
        ('threads', threads)])
    times = {}
    wd = tempfile.mkdtemp(prefix='units_panel.%s.' % rid, dir=workroot)
    try:
        full_msa = os.path.join(wd, 'full.msa.fa')
        al = realign_units.align_units_fasta(fa, rjp, full_msa, threads=threads, params=P, workdir=wd)
        al.pop('params', None)
        info['params'] = P.as_dict()
        info['align'] = al
        times['align_s'] = al.get('align_s')
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
                st['paths_verified'] = check_paths(tmp_gfa, msa_graph.read_fasta(hap_fa))
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
        realign.align_fasta = orig_align_fasta
        shutil.rmtree(wd, ignore_errors=True)
        for p in os.listdir(o['cand_dir']):
            if p.startswith(rid + '.') and '.tmp' in p:
                os.remove(os.path.join(o['cand_dir'], p))
    info.update(times)
    info['seconds'] = round(time.time() - t0, 2)
    info['finished'] = datetime.datetime.now().isoformat(timespec='seconds')
    info['host_load'] = round(os.getloadavg()[0], 1)
    info['units_panel_version'] = 1
    write_json(o['json'], info)
    row = summary_row(rid, info_r, info, times, info['status'], info['seconds'], threads, timeout, mem_mb)
    if record_runtime:
        update_runtime(row)
    return row


def record_killed(rid, status, seconds, peak, a, why):
    """The job was killed from outside (timeout / memout / crash): write its realign.json and row."""
    info_r = region_info(rid)
    o = out_paths(rid, a.cand_root, a.panel_root)
    os.makedirs(o['cand_dir'], exist_ok=True)
    for k in ('msa', 'gfa', 'full_msa', 'full_gfa'):
        if os.path.exists(o[k]):
            os.remove(o[k])
    for d in (o['cand_dir'], o['panel_dir']):
        if os.path.isdir(d):
            for p in os.listdir(d):
                if p.startswith(rid + '.') and '.tmp' in p:
                    os.remove(os.path.join(d, p))
    info = collections.OrderedDict([
        ('method', METHOD), ('panel', 'all'), ('region_id', rid), ('stratum', info_r['stratum']),
        ('span_bp', info_r['span_bp']), ('status', status), ('message', why),
        ('n_distinct', info_r['n_distinct']), ('n_hap32_only', info_r['n_hap32_only']),
        ('n_panel_records', info_r['n_panel_records']), ('max_len', info_r['max_len']),
        ('total_bp', info_r['total_bp']), ('timeout_s', a.timeout), ('mem_cap_mb', a.mem_mb),
        ('threads', a.threads), ('peak_rss_mb', round(peak, 0)), ('seconds', round(seconds, 1)),
        ('finished', datetime.datetime.now().isoformat(timespec='seconds')), ('units_panel_version', 1)])
    write_json(o['json'], info)
    row = summary_row(rid, info_r, info, {}, status, round(seconds, 1), a.threads, a.timeout, a.mem_mb,
                      peak=peak, note=why)
    if not a.no_runtime:
        update_runtime(row)
    return row


# ---------------------------------------------------------------- process-tree control

def _ps_table():
    try:
        out = subprocess.run(['ps', '-A', '-o', 'pid=,ppid=,pgid=,rss='], capture_output=True, text=True,
                             timeout=20).stdout
    except (OSError, subprocess.SubprocessError):
        return []
    rows = []
    for line in out.splitlines():
        f = line.split()
        if len(f) == 4:
            try:
                rows.append(tuple(int(x) for x in f))
            except ValueError:
                pass
    return rows


def tree_of(pid, table=None):
    """[(pid, pgid, rss_kb)] of pid and all its descendants (mafft runs in sessions of its own)."""
    table = table if table is not None else _ps_table()
    kids = collections.defaultdict(list)
    info = {}
    for p, pp, g, rss in table:
        kids[pp].append(p)
        info[p] = (g, rss)
    out, todo = [], [pid]
    while todo:
        p = todo.pop()
        if p in info:
            out.append((p, info[p][0], info[p][1]))
        todo.extend(kids.get(p, ()))
    return out


def kill_tree(proc):
    """SIGTERM then SIGKILL every process group and process of proc's tree."""
    for sig in (signal.SIGTERM, signal.SIGKILL):
        tree = tree_of(proc.pid)
        if not tree and proc.poll() is not None:
            return
        for p, g, _ in tree:
            for fn, x in ((os.killpg, g), (os.kill, p)):
                try:
                    fn(x, sig)
                except (ProcessLookupError, PermissionError, OSError):
                    pass
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            pass
        time.sleep(1)
        if not tree_of(proc.pid) and proc.poll() is not None:
            return


# ---------------------------------------------------------------- scheduler

def pick_regions(arg, stratum=None):
    rows = realign.region_list()
    if arg and arg != 'all':
        want = [x.strip() for x in arg.split(',') if x.strip()]
        known = set(r['region_id'] for r in rows)
        missing = [x for x in want if x not in known]
        if missing:
            raise SystemExit('units_panel: no region package for %s' % ', '.join(missing))
        rows = [r for r in rows if r['region_id'] in set(want)]
    if stratum:
        keep = set(stratum.split(','))
        rows = [r for r in rows if r.get('stratum') in keep]
    return [r['region_id'] for r in rows]


def job_argv(rid, a, result):
    cmd = [sys.executable, os.path.abspath(__file__), 'one', rid, '--timeout', str(a.timeout),
           '--mem-mb', str(a.mem_mb), '--threads', str(a.threads), '--result', result]
    for opt in ('cand_root', 'panel_root', 'workdir'):
        v = getattr(a, opt)
        if v:
            cmd += ['--' + opt.replace('_', '-'), v]
    if a.no_runtime:
        cmd.append('--no-runtime')
    return cmd


def set_peak(rid, peak, a):
    """Record the parent-sampled peak RSS of a finished job in its realign.json."""
    o = out_paths(rid, a.cand_root, a.panel_root)
    try:
        with open(o['json']) as f:
            info = json.load(f, object_pairs_hook=collections.OrderedDict)
    except (OSError, ValueError):
        return
    info['peak_rss_mb'] = peak
    write_json(o['json'], info)


def rebuild_table(a, log_path=None):
    """Rewrite the runtime table from every candidates/unit_aware__all/<id>.realign.json (one
    schema for all rows); peak RSS comes from the JSON, else from a batch log's 'rss N MB'."""
    import re
    peaks = {}
    if log_path:
        for line in open(log_path):
            m = re.match(r'\[units_panel\s+[\d.]+s\] (\S+)\s.*\brss\s+([\d.]+) MB', line)
            if m:
                peaks[m.group(1)] = float(m.group(2))
    rows = []
    for rid in pick_regions('all'):
        o = out_paths(rid, a.cand_root, a.panel_root)
        try:
            with open(o['json']) as f:
                info = json.load(f)
        except (OSError, ValueError):
            continue
        times = dict((k, info[k]) for k in ('align_s', 'panel_graph_s', 'project_s', 'graph_s') if k in info)
        peak = info.get('peak_rss_mb', peaks.get(rid))
        info_r = region_info(rid)
        row = summary_row(rid, info_r, info, times, info.get('status'), info.get('seconds'), info.get('threads'),
                          info.get('timeout_s'), info.get('mem_cap_mb'), peak=peak)
        row['finished'] = info.get('finished', '')
        row['host_load'] = info.get('host_load', '')
        rows.append(row)
    path = RUNTIME_TSV
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + '.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        tmp = path + '.tmp%d' % os.getpid()
        with open(tmp, 'w', newline='') as f:
            w = csv.DictWriter(f, RUNTIME_COLS, delimiter='\t', extrasaction='ignore', lineterminator='\n')
            w.writeheader()
            for r in sorted(rows, key=lambda r: (r['method'], r['region_id'])):
                w.writerow(dict((c, r.get(c, '')) for c in RUNTIME_COLS))
        os.replace(tmp, path)
        fcntl.flock(lk, fcntl.LOCK_UN)
    return len(rows)


def available_mb():
    try:
        import realign_poa
        return realign_poa.available_mb()
    except Exception:  # noqa: BLE001
        return None


def run_batch(a):
    rids = pick_regions(a.regions, a.stratum)
    infos = dict((rid, region_info(rid)) for rid in rids)
    rids.sort(key=lambda r: (infos[r]['total_bp'], infos[r]['max_len'], r))
    pending = []
    for rid in rids:
        prev = previous_status(rid, a.cand_root, a.panel_root)
        if prev and not a.force and not (a.retry_failed and prev != 'ok'):
            continue
        pending.append(rid)
    log('%d regions to run (of %d); %d at a time, %d threads, timeout %d s, cap %d MB'
        % (len(pending), len(rids), a.jobs, a.threads, a.timeout, a.mem_mb))
    logdir = os.path.join(a.workdir or os.path.join(config.WORK_DIR, 'realign_units'), 'units_panel_logs')
    os.makedirs(logdir, exist_ok=True)
    running = {}
    counts = collections.Counter()
    last_wait = 0
    while pending or running:
        for p in list(running):
            j = running[p]
            if p.poll() is not None:
                running.pop(p)
                try:
                    with open(j['result']) as f:
                        row = json.load(f)
                    row['peak_rss_mb'] = round(j['peak'], 0)      # sampled here, over the whole tree
                    set_peak(j['rid'], row['peak_rss_mb'], a)
                    if not a.no_runtime:
                        update_runtime(row)
                except (OSError, ValueError):
                    row = record_killed(j['rid'], 'crash', time.time() - j['t0'], j['peak'], a,
                                        'no result; exit %s; log %s' % (p.returncode, j['log']))
                counts[row.get('status')] += 1
                log('%-9s %-14s %-9s %-8s %7ss align %7ss  rss %6.0f MB  n %4s  nodes full %-6s proj %-6s %s' % (
                    j['rid'], infos[j['rid']]['stratum'], row.get('status'), row.get('mode', ''), row.get('seconds', ''),
                    row.get('align_s', ''), j['peak'], row.get('n_distinct', ''), row.get('nodes_full', ''),
                    row.get('nodes_proj', ''), (row.get('reason') or row.get('note') or '')[:80]))
                continue
            now = time.time()
            if now - j['last_ps'] >= 3:
                j['last_ps'] = now
                rss = sum(r for _, _, r in tree_of(p.pid)) / 1024.0
                j['peak'] = max(j['peak'], rss)
                why = None
                if now - j['t0'] > a.timeout:
                    why = ('timeout', 'killed after %.0f s (cap %d s)' % (now - j['t0'], a.timeout))
                elif rss > a.mem_mb:
                    why = ('memout', 'killed at %.0f MB RSS (cap %d MB)' % (rss, a.mem_mb))
                if why:
                    kill_tree(p)
                    running.pop(p)
                    row = record_killed(j['rid'], why[0], time.time() - j['t0'], j['peak'], a, why[1])
                    counts[why[0]] += 1
                    log('%-9s %-14s %-9s %s' % (j['rid'], infos[j['rid']]['stratum'], why[0], why[1]))
        started = False
        if pending and len(running) < a.jobs:
            rid = pending[0]
            need = infos[rid]['predicted_mb']
            avail = available_mb()
            if not running or avail is None or avail >= need + 1500:
                pending.pop(0)
                result = os.path.join(logdir, '%s.result.json' % rid)
                lg = os.path.join(logdir, '%s.log' % rid)
                if os.path.exists(result):
                    os.remove(result)
                with open(lg, 'w') as lf:
                    p = subprocess.Popen(job_argv(rid, a, result), stdout=lf, stderr=subprocess.STDOUT,
                                         start_new_session=True, cwd=config.REPO)
                running[p] = {'rid': rid, 't0': time.time(), 'last_ps': 0, 'peak': 0.0, 'result': result, 'log': lg}
                started = True
            elif time.time() - last_wait > 300:
                log('waiting for memory: next %s needs ~%.0f MB, available %.0f MB' % (rid, need + 1500, avail))
                last_wait = time.time()
        if not started:
            time.sleep(1.0)
    log('done: %s' % ', '.join('%s %d' % kv for kv in sorted(counts.items(), key=lambda kv: str(kv[0]))))
    return 0


# ---------------------------------------------------------------- CLI

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    def common(p):
        p.add_argument('--timeout', type=int, default=DEFAULT_TIMEOUT, help='wall clock per region (s)')
        p.add_argument('--mem-mb', type=int, default=DEFAULT_MEM_MB, help='RSS cap per region (MB)')
        p.add_argument('--threads', type=int, default=DEFAULT_THREADS, help='mafft threads per region')
        p.add_argument('--cand-root', help='candidates root (default %s)' % config.CANDIDATES_DIR)
        p.add_argument('--panel-root', help='full-panel outputs root (default %s)' % panel.PANEL_DIR)
        p.add_argument('--workdir', help='scratch (default $VNTR_WORK/realign_units)')
        p.add_argument('--no-runtime', action='store_true', help='do not write the runtime table')
    p = sub.add_parser('run', help='a batch of regions, smallest first')
    p.add_argument('--regions', default='all')
    p.add_argument('--stratum')
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--force', action='store_true')
    p.add_argument('--retry-failed', action='store_true')
    common(p)
    p = sub.add_parser('one', help='one region in this process')
    p.add_argument('region')
    p.add_argument('--result', help='also write the runtime row here (JSON)')
    common(p)
    p = sub.add_parser('table', help='rewrite the runtime table from the realign.json files')
    p.add_argument('--log', help='a batch log to take peak RSS from where the JSON has none')
    p.add_argument('--cand-root')
    p.add_argument('--panel-root')
    a = ap.parse_args(argv)
    if a.cmd == 'table':
        n = rebuild_table(a, a.log)
        log('%d rows -> %s' % (n, RUNTIME_TSV))
        return 0
    if a.cmd == 'one':
        rid = os.path.basename(os.path.normpath(a.region))
        row = run_job(rid, threads=a.threads, timeout=a.timeout, mem_mb=a.mem_mb, cand_root=a.cand_root,
                      panel_root=a.panel_root, workroot=a.workdir, record_runtime=not a.no_runtime)
        if a.result:
            with open(a.result + '.tmp', 'w') as f:
                json.dump(row, f)
            os.replace(a.result + '.tmp', a.result)
        log('%s: %s (%s) in %.1f s %s' % (rid, row['status'], row.get('mode'), row['seconds'], row.get('note', '')))
        return 0 if row['status'] == 'ok' else 2
    return run_batch(a)


if __name__ == '__main__':
    sys.exit(main())
