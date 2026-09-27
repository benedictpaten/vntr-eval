#!/usr/bin/env python3
"""poa_panel_mc.py -- the full-panel arm of poa_abpoa_mc (tools/realign_poa_mc.py).

tools/poa_panel.py runs poa_abpoa / poa_spoa on the full panel, but it cannot take poa_abpoa_mc
without changing its own behaviour: it registers its methods through realign_poa's variant table,
aligns the deduplicated union, and predicts memory from abPOA's default band. This is the thin
driver for the one extra method. It reuses poa_panel's union files and path checks and panel.py's
projection and panel graph, so the outputs have the same layout:

    work/panel/union/<id>.fa + .map.tsv     tools/panel.py union (made here when missing)
    realign_poa_mc.align_fasta              Cactus's abPOA call and BAR's 10 kb end rule, N masking and
                                            row checks as in the hap32 arm (see that module)
    -> <panel-root>/poa_abpoa_mc/<id>.msa.fa.gz   the full MSA: one row per aligned record, named by
                                            its hprc.fa.gz name (CHM13/GRCh38 and hap32-only rows by
                                            their hap32.fa name), in input order
    -> <panel-root>/poa_abpoa_mc/<id>.gfa   the full-panel graph (panel.py panel-graph), one P line per
                                            hprc.fa.gz record; every path's spelling re-read and checked
    -> <cand-root>/poa_abpoa_mc__all/<id>.msa.fa    the hap32 rows (panel.py project; all-gap columns dropped)
    -> <cand-root>/poa_abpoa_mc__all/<id>.gfa       msa_graph.py graph of those; every path re-read and checked
    -> <cand-root>/poa_abpoa_mc__all/<id>.realign.json   how it was made (also for memouts, timeouts,
                                            skips and failures, so a re-run skips them)
    -> results/mcpoa_runtime.all.tsv        one row per region; the latest run wins

What abPOA aligns (one difference from poa_panel, which aligns the distinct sequences once each):
every hprc.fa.gz record, duplicates included -- Cactus aligns every haplotype as a row, and both
the -p guide tree (minimizer Jaccard over all rows) and BAR's column scores count copies -- plus
each hap32.fa sequence that no full-panel record carries (a sampled recombinant; needed for the
projection), once. Rows are grouped by distinct sequence in the union's order (CHM13 first, then
decreasing weight; hprc.fa.gz order within a group); abPOA then takes them longest first and -p
re-orders them by its guide tree. The panel graph and the projection keep the first row of each
distinct sequence (panel.py's rule); `dup_split` counts distinct sequences whose copies abPOA
aligned differently (their other rows are in the gzipped full MSA, not in either graph).

Usage (from the repository root):

    python3 tools/poa_panel_mc.py run [--regions ID,..|all] [--stratum S] [--jobs 2] [--timeout 1800]
        [--mem-mb 10000] [--budget-mb 20000] [--force | --retry-failed]
        [--cand-root DIR] [--panel-root DIR] [--no-runtime]
    python3 tools/poa_panel_mc.py one REGION [same options]      # one job in this process
    python3 tools/poa_panel_mc.py plan [--regions ..]            # records, rule, predicted MB per region

--timeout caps the summed abPOA wall clock of a region (both end calls together: the second call
gets what the first left). --mem-mb stops a region whose predicted peak is above it (status memout,
"predicted" in the note, not run) and kills an abPOA call whose sampled RSS passes it (memout).

Scheduling: jobs are separate processes, at most --jobs at a time, smallest region first (the bp
abPOA is given: the sum over its calls of the strings' lengths). A job starts only when the
predicted memory of the running jobs plus its own (1.25 x prediction, at most --mem-mb) stays
within --budget-mb and, for a job predicted above 4 GB, memory_pressure reports its prediction
plus 2 GB free; a smaller job may overtake one that waits. A job that has waited --max-wait
seconds with nothing else of ours running is recorded as skipped_memory.

Memory model (MB) per abPOA call, realign_poa_mc.predict_call_mb with the full-panel constants
below: 60 + TOUCHED x (20 B x (L1 + 4) x (NODES_PER_BP x L1 + 2), rounded up to a power of two),
L1 the longest string of the call (at most 10 kb under BAR's rule). abPOA allocates 5 int32 DP
matrices of graph nodes x query length and rounds the allocation up to a power of two; hundreds of
rows grow a graph of more nodes per bp than hap32's 34 (calibration below). The prediction is
conservative on purpose; memory_pressure, not the budget, is the real gate for big jobs.
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
import poa_panel  # noqa: E402
import realign  # noqa: E402
import realign_poa_mc as mc  # noqa: E402

METHOD = mc.METHOD
ALL_METHOD = METHOD + '__all'
RUNTIME_TSV = os.path.join(config.RESULTS_DIR, 'mcpoa_runtime.all.tsv')
RUNTIME_COLS = ['region_id', 'method', 'panel', 'status', 'stratum', 'span_bp', 'n_records', 'n_distinct',
                'n_hap32_only', 'max_interior', 'len2', 'aligned_bp', 'bar_rule', 'rows_ge_10kb', 'rows_gt_20kb',
                'unaligned_bp', 'unaligned_rows', 'abpoa_calls', 'progressive', 'predicted_mb', 'peak_rss_mb',
                'sampled_rss_mb', 'waited_s', 'free_mb_at_start', 'mem_cap_mb', 'timeout_s', 'seconds', 'align_s',
                'abpoa_s', 'abpoa_nodes_max', 'columns_full', 'dup_split', 'panel_graph_s', 'project_s', 'graph_s',
                'nodes_full', 'edges_full', 'columns_proj', 'nodes_proj', 'nodes_per_kb_proj', 'frac_nodes_1bp_proj',
                'masked_runs', 'masked_bp', 'finished', 'note']
DEFAULT_TIMEOUT = 1800
DEFAULT_MEM_MB = 10000
DEFAULT_BUDGET_MB = 20000
WAIT_ABOVE_MB = 4000
# The MSA-based pipeline (realign.align_fasta, panel.project/panel_graph) holds the full MSA as Python
# strings; a region whose unaligned middles would make it larger than this many cells is not run.
MAX_MSA_CELLS = 2e9
# Full-panel calibration, /usr/bin/time -l on 7 calls of 407-459 rows (L003826, L014297, L000034 single;
# L005990 and L004145 both ends, 10 kb strings): abPOA nodes per bp of the longest string 1.28-2.19, and
# (peak - 60 MB) over the power-of-two allocation 0.97-1.09 -- abPOA's band touches nearly all of it at
# b=1000, f=0.1. A 10 kb call sits at a step: below 21.4k nodes the allocation is 4 GB (measured peaks
# 3.9-4.4 GB), above it 8 GB. L004145's right end ended at 21.9k nodes (2.19/bp) but its last
# reallocation stayed under the step (3.9 GB). 2.2 nodes/bp puts every 10 kb call on the 8 GB step
# (predicted 9.1 GB, under the 10 GB cap), so the memory_pressure gate always asks for that much.
# The full run (2026-09-26, 184 calls; results/mcpoa_runtime.all.tsv) then measured 1.0-4.1 nodes/bp
# (median 1.2; the 8 GB step holds up to 4.29) and peaks of 0-1.6x the allocation (median 0.82): measured
# over predicted was median 0.40, max 1.23 (L003828, a 4.6 GB prediction; L013531 1.11). No 10 kb call
# went above 0.78x its prediction (max 7.1 GB). A larger TOUCHED would push every 10 kb call over the
# 10 GB cap, so the constants were left as run; the sampled-RSS kill at --mem-mb is the real bound.
PANEL_NODES_PER_BP = 2.2
PANEL_TOUCHED = 1.1
T0 = time.time()


def log(msg):
    sys.stderr.write('[poa_panel_mc %7.1fs] %s\n' % (time.time() - T0, msg))
    sys.stderr.flush()


def predict_call_mb(lens):
    return mc.predict_call_mb(lens, nodes_per_bp=PANEL_NODES_PER_BP, touched=PANEL_TOUCHED)


def out_paths(rid, cand_root=None, panel_root=None):
    cd = os.path.join(cand_root or config.CANDIDATES_DIR, ALL_METHOD)
    pd = os.path.join(panel_root or panel.PANEL_DIR, METHOD)
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


# ---------------------------------------------------------------- inputs and plan

def panel_records(rd, map_tsv):
    """[(name, sequence, distinct id)]: every hprc.fa.gz record, then (in its union slot) each hap32
    sequence no full-panel record carries, once; grouped by distinct sequence in union order."""
    hprc = [(n, s.upper()) for n, _, s in panel.read_fasta(os.path.join(rd, 'hprc.fa.gz'))]
    hap = [(n, s.upper()) for n, s in msa_graph.read_fasta(os.path.join(rd, 'hap32.fa'))]
    mp = panel.read_map(map_tsv)
    uid_of = {}
    for r in mp:
        for m in r['members']:
            uid_of[m] = r['id']
    by_uid = collections.defaultdict(list)
    for n, s in hprc:
        by_uid[uid_of[n]].append((n, s))
    hap_seq = dict(hap)
    out = []
    for r in mp:
        if int(r['weight']) > 0:
            recs = by_uid[r['id']]
            if len(recs) != int(r['weight']):
                raise RuntimeError('%s: %s has %d hprc.fa.gz records, map says %s' % (rd, r['id'], len(recs), r['weight']))
        else:
            n = r['members'][0]
            recs = [(n, hap_seq[n])]
        for n, s in recs:
            if len(s) != int(r['length']):
                raise RuntimeError('%s: %s length %d, map says %s' % (rd, n, len(s), r['length']))
            out.append((n, s, r['id']))
    return out


def plan_region(rid):
    rd = realign.region_dir(rid)
    fa, mp = poa_panel.ensure_union(rd, rid)
    recs = panel_records(rd, mp)
    L, R = mc.region_flanks(rd)
    lens = sorted((len(realign.mask_runs(s)[0]) - L - R for _, s, _ in recs), reverse=True)
    if lens[0] < mc.BANDING_LIMIT:
        call_lens, rule = [lens], 'single'
    else:
        pre = [min(x, mc.BANDING_LIMIT) for x in lens]
        call_lens, rule = [pre, pre], 'ends'
    stratum, span = realign.region_meta(rd)
    maprows = panel.read_map(mp)
    return {'region_id': rid, 'rd': rd, 'union_fa': fa, 'union_map': mp, 'stratum': stratum, 'span_bp': span,
            'flanks': [L, R], 'n_records': len(recs), 'n_distinct': len(maprows),
            'n_hap32_only': sum(1 for r in maprows if int(r['weight']) == 0),
            'max_interior': lens[0], 'len2': lens[1] if len(lens) > 1 else 0, 'bar_rule': rule,
            'rows_ge_10kb': sum(1 for x in lens if x >= mc.BANDING_LIMIT),
            'rows_gt_20kb': sum(1 for x in lens if x > 2 * mc.BANDING_LIMIT),
            'unaligned_bp_planned': sum(max(0, x - 2 * mc.BANDING_LIMIT) for x in lens),
            'aligned_bp': sum(sum(c) for c in call_lens),
            'predicted_mb': max(predict_call_mb(c) for c in call_lens),
            # columns ~ 1.5 x the longest string per call (measured 1.0-2.0), plus one column per unaligned bp
            'msa_cells_predicted': float(len(recs)) * (1.5 * sum(max(c) for c in call_lens)
                                                       + sum(max(0, x - 2 * mc.BANDING_LIMIT) for x in lens))}


# ---------------------------------------------------------------- runtime table

def update_runtime(row, path=RUNTIME_TSV):
    """Replace the region's row of the runtime table, under a file lock."""
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

def dup_split(full_msa, recs):
    """Distinct sequences whose copies abPOA aligned to different rows."""
    rows = dict(msa_graph.read_msa(full_msa))
    seen = collections.defaultdict(set)
    for n, _, u in recs:
        seen[u].add(rows[n])
    return sum(1 for v in seen.values() if len(v) > 1)


def run_job(rid, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB, cand_root=None, panel_root=None, workroot=None,
            record_runtime=True, wait_memory=True, max_wait=3600, runtime_tsv=RUNTIME_TSV):
    """Align region rid's full panel with poa_abpoa_mc; write every output; return the runtime row."""
    t0 = time.time()
    p = plan_region(rid)
    rd, mp = p['rd'], p['union_map']
    hap_fa = os.path.join(rd, 'hap32.fa')
    o = out_paths(rid, cand_root, panel_root)
    os.makedirs(o['cand_dir'], exist_ok=True)
    os.makedirs(o['panel_dir'], exist_ok=True)
    workroot = workroot or os.path.join(config.WORK_DIR, 'mcpoa', 'tmp')
    os.makedirs(workroot, exist_ok=True)
    for k in ('msa', 'gfa', 'full_msa', 'full_gfa'):     # never leave an older result beside a new attempt
        if os.path.exists(o[k]):
            os.remove(o[k])
    info = collections.OrderedDict([
        ('method', ALL_METHOD), ('base_method', METHOD), ('panel', 'all'), ('region_id', rid),
        ('stratum', p['stratum']), ('span_bp', p['span_bp']), ('status', None), ('description', mc.DESCRIPTION),
        ('input', 'every hprc.fa.gz record (duplicates included) + hap32-only sequences once'),
        ('union_fa', realign._rel(p['union_fa'])), ('union_map', realign._rel(mp))])
    for k in ('flanks', 'n_records', 'n_distinct', 'n_hap32_only', 'max_interior', 'len2', 'bar_rule', 'rows_ge_10kb',
              'rows_gt_20kb', 'unaligned_bp_planned', 'aligned_bp', 'predicted_mb', 'msa_cells_predicted'):
        info[k] = p[k]
    info.update(mem_cap_mb=mem_mb, timeout_s=timeout, predictor={'nodes_per_bp': PANEL_NODES_PER_BP,
                                                                 'touched': PANEL_TOUCHED})
    times, waited, free = {}, 0.0, None
    wd = tempfile.mkdtemp(prefix='poa_panel_mc.%s.' % rid, dir=workroot)
    try:
        if p['msa_cells_predicted'] > MAX_MSA_CELLS:
            info['status'] = 'skipped'
            info['message'] = ('not run: BAR would leave %d bp unaligned over %d rows, a full MSA of ~%.1e cells '
                               '(cap %.0e) that the MSA-based pipeline cannot hold' % (
                                   p['unaligned_bp_planned'], p['rows_gt_20kb'], p['msa_cells_predicted'],
                                   MAX_MSA_CELLS))
        elif mem_mb and p['predicted_mb'] > mem_mb:
            info['status'] = 'memout'
            info['message'] = 'not run: predicted %.0f MB > cap %.0f MB (longest string %d, rule %s)' % (
                p['predicted_mb'], mem_mb, min(p['max_interior'], mc.BANDING_LIMIT), p['bar_rule'])
        else:
            if wait_memory and p['predicted_mb'] > WAIT_ABOVE_MB:
                waited, free = mc.wait_for_memory(p['predicted_mb'] + 2000, max_wait=max_wait, log=log)
                if free is not None and free < p['predicted_mb'] + 2000:
                    info['status'] = 'skipped_memory'
                    info['message'] = 'not run: memory_pressure reports %.0f MB free after %.0f s, need %.0f MB' % (
                        free, waited, p['predicted_mb'] + 2000)
            if info['status'] is None:
                recs = panel_records(rd, mp)
                in_fa = os.path.join(wd, 'in.fa')
                with open(in_fa, 'w') as f:
                    for n, s, _ in recs:
                        f.write('>%s\n%s\n' % (n, s))
                full_msa = os.path.join(wd, 'full.msa.fa')
                al = mc.align_fasta(in_fa, full_msa, tuple(p['flanks']), timeout=timeout, mem_mb=mem_mb, workdir=wd,
                                    region_timeout=timeout)
                al.pop('input', None)
                al.pop('msa', None)
                info['align'] = al
                info['status'] = al['status']
                times['align_s'] = al.get('align_s')
                if al.get('message'):
                    info['message'] = al['message']
                if al['status'] == 'ok':
                    try:
                        info['dup_split'] = dup_split(full_msa, recs)
                        tmp = o['full_msa'] + '.tmp%d' % os.getpid()
                        with open(full_msa, 'rb') as fi, gzip.open(tmp, 'wb', compresslevel=6) as fo:
                            shutil.copyfileobj(fi, fo)
                        os.replace(tmp, o['full_msa'])
                        tg = time.time()
                        pg = panel.panel_graph(full_msa, mp, rd, o['full_gfa'], engine='native')
                        hprc = [(n, s) for n, _, s in panel.read_fasta(os.path.join(rd, 'hprc.fa.gz'))]
                        pg['paths_verified'] = poa_panel.check_paths(o['full_gfa'], hprc)
                        pg['gfa'] = realign._rel(o['full_gfa'])
                        pg['msa'] = realign._rel(o['full_msa'])
                        info['panel_graph'] = pg
                        times['panel_graph_s'] = round(time.time() - tg, 2)
                        tp = time.time()
                        tmp_msa = o['msa'] + '.tmp%d' % os.getpid()
                        info['projection'] = panel.project(full_msa, mp, hap_fa, tmp_msa)
                        times['project_s'] = round(time.time() - tp, 2)
                        th = time.time()
                        tmp_gfa = o['gfa'] + '.tmp%d' % os.getpid()
                        st = msa_graph.msa_to_gfa(tmp_msa, hap_fa, tmp_gfa, workdir=wd)
                        st['paths_verified'] = poa_panel.check_paths(tmp_gfa, msa_graph.read_fasta(hap_fa))
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
        for x in os.listdir(o['cand_dir']):
            if x.startswith(rid + '.') and '.tmp' in x:
                os.remove(os.path.join(o['cand_dir'], x))
    info.update(times)
    info['memory_wait_s'] = waited
    info['memory_free_mb_at_start'] = free
    info['seconds'] = round(time.time() - t0, 2)
    info['finished'] = datetime.datetime.now().isoformat(timespec='seconds')
    info['host_load'] = round(os.getloadavg()[0], 1)
    tmpj = o['json'] + '.tmp%d' % os.getpid()
    with open(tmpj, 'w') as f:
        json.dump(info, f, indent=1)
    os.replace(tmpj, o['json'])
    al = info.get('align') or {}
    aln = al.get('aligner') or {}
    calls = aln.get('calls') or []
    pg = info.get('panel_graph') or {}
    g = info.get('graph') or {}
    pr = info.get('projection') or {}
    row = {'region_id': rid, 'method': ALL_METHOD, 'panel': 'all', 'status': info['status'],
           'stratum': p['stratum'], 'span_bp': p['span_bp'], 'n_records': p['n_records'], 'n_distinct': p['n_distinct'],
           'n_hap32_only': p['n_hap32_only'], 'max_interior': p['max_interior'], 'len2': p['len2'],
           'aligned_bp': p['aligned_bp'], 'bar_rule': p['bar_rule'], 'rows_ge_10kb': p['rows_ge_10kb'],
           'rows_gt_20kb': p['rows_gt_20kb'], 'unaligned_bp': aln.get('unaligned_bp', ''),
           'unaligned_rows': aln.get('unaligned_rows', ''), 'abpoa_calls': aln.get('abpoa_calls', len(calls) or ''),
           'progressive': aln.get('progressive', ''), 'predicted_mb': p['predicted_mb'],
           'peak_rss_mb': aln.get('peak_rss_mb', ''), 'sampled_rss_mb': aln.get('sampled_rss_mb', ''),
           'waited_s': waited, 'free_mb_at_start': '' if free is None else round(free),
           'mem_cap_mb': mem_mb, 'timeout_s': timeout, 'seconds': info['seconds'], 'align_s': times.get('align_s', ''),
           'abpoa_s': aln.get('abpoa_s', ''),
           'abpoa_nodes_max': max([c.get('nodes') or 0 for c in calls] or [0]) or '',
           'columns_full': al.get('columns', ''), 'dup_split': info.get('dup_split', ''),
           'panel_graph_s': times.get('panel_graph_s', ''), 'project_s': times.get('project_s', ''),
           'graph_s': times.get('graph_s', ''), 'nodes_full': pg.get('nodes', ''), 'edges_full': pg.get('edges', ''),
           'columns_proj': pr.get('columns_out', ''), 'nodes_proj': g.get('nodes', ''),
           'nodes_per_kb_proj': g.get('nodes_per_kb', ''), 'frac_nodes_1bp_proj': g.get('frac_nodes_1bp', ''),
           'masked_runs': al.get('masked_runs', ''), 'masked_bp': al.get('masked_bp', ''),
           'finished': info['finished'],
           'note': (info.get('message') or '')[:200].replace('\t', ' ').replace('\n', ' ')}
    if record_runtime:
        update_runtime(row, runtime_tsv)
    return row


# ---------------------------------------------------------------- scheduler

def job_argv(rid, a):
    cmd = [sys.executable, os.path.abspath(__file__), 'one', rid, '--timeout', str(a.timeout), '--mem-mb',
           str(a.mem_mb), '--no-wait-memory']
    for opt in ('cand_root', 'panel_root', 'workdir'):
        v = getattr(a, opt)
        if v:
            cmd += ['--' + opt.replace('_', '-'), v]
    if a.no_runtime:
        cmd.append('--no-runtime')
    return cmd


def run_batch(a):
    rids = poa_panel.pick_regions(a.regions, a.stratum)
    plans = {rid: plan_region(rid) for rid in rids}
    rids.sort(key=lambda r: (plans[r]['aligned_bp'], plans[r]['max_interior'], r))
    jobs = []
    for rid in rids:
        prev = previous_status(rid, a.cand_root, a.panel_root)
        if prev and not a.force and not (a.retry_failed and prev != 'ok'):
            continue
        pred = plans[rid]['predicted_mb']
        jobs.append({'rid': rid, 'pred': pred, 'reserve': min(1.25 * pred, a.mem_mb)})
    log('%d jobs of %d regions; budget %d MB, cap %d MB, %d at a time, region abPOA cap %d s'
        % (len(jobs), len(rids), a.budget_mb, a.mem_mb, a.jobs, a.timeout))
    logdir = os.path.join(a.workdir or os.path.join(config.WORK_DIR, 'mcpoa', 'tmp'), 'poa_panel_mc_logs')
    os.makedirs(logdir, exist_ok=True)
    running, counts, pending = {}, collections.Counter(), list(jobs)
    last_wait, wait_since = 0, {}
    while pending or running:
        for pr in list(running):
            if pr.poll() is not None:
                j = running.pop(pr)
                try:
                    row = json.load(open(j['result']))
                except (OSError, ValueError):
                    row = {'status': 'crash', 'note': 'no result; exit %s; log %s' % (pr.returncode, j['log'])}
                    if not a.no_runtime:
                        update_runtime({'region_id': j['rid'], 'method': ALL_METHOD, 'panel': 'all', 'status': 'crash',
                                        'seconds': round(time.time() - j['t0'], 1), 'note': row['note']})
                counts[row.get('status')] += 1
                log('%-9s %-14s %7ss abpoa %7ss peak %7s MB (pred %s) rule %-6s unaligned %-7s nodes full %-6s '
                    'proj %-6s dup_split %s %s' % (
                        j['rid'], row.get('status'), row.get('seconds', ''), row.get('abpoa_s', ''),
                        row.get('peak_rss_mb', ''), row.get('predicted_mb', ''), row.get('bar_rule', ''),
                        row.get('unaligned_bp', ''), row.get('nodes_full', ''), row.get('nodes_proj', ''),
                        row.get('dup_split', ''), (row.get('note') or '')[:120]))
            elif time.time() - running[pr]['t0'] > a.timeout + 1800:
                log('killing %s: job wall clock over %d s' % (running[pr]['rid'], a.timeout + 1800))
                realign._kill_group(pr)
        started = False
        if pending and len(running) < a.jobs:
            reserved = sum(j['reserve'] for j in running.values())
            free = None
            for k, j in enumerate(pending):
                fits = reserved + j['reserve'] <= a.budget_mb or j['pred'] > a.mem_mb
                if fits and WAIT_ABOVE_MB < j['pred'] <= a.mem_mb:
                    if free is None:
                        free = mc.memory_free_mb()
                    fits = free is None or free >= j['pred'] + 2000
                if not fits:
                    wait_since.setdefault(j['rid'], time.time())
                    if not running and time.time() - wait_since[j['rid']] > a.max_wait:
                        pending.pop(k)
                        note = 'not run: waited %.0f s for memory (memory_pressure free %s MB, need %.0f)' % (
                            time.time() - wait_since[j['rid']], free, j['pred'] + 2000)
                        log('%s skipped_memory: %s' % (j['rid'], note))
                        counts['skipped_memory'] += 1
                        if not a.no_runtime:
                            update_runtime({'region_id': j['rid'], 'method': ALL_METHOD, 'panel': 'all',
                                            'status': 'skipped_memory', 'predicted_mb': j['pred'], 'note': note})
                        break
                    continue
                pending.pop(k)
                j['result'] = os.path.join(logdir, '%s.result.json' % j['rid'])
                j['log'] = os.path.join(logdir, '%s.log' % j['rid'])
                if os.path.exists(j['result']):
                    os.remove(j['result'])
                cmd = job_argv(j['rid'], a) + ['--result', j['result']]
                j['t0'] = time.time()
                with open(j['log'], 'w') as lf:
                    pr = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT, start_new_session=True,
                                          cwd=config.REPO)
                running[pr] = j
                started = True
                break
            if not started and pending and time.time() - last_wait > 300:
                log('waiting for memory: %d pending, next %s needs %.0f MB, reserved %.0f, memory_pressure free %s'
                    % (len(pending), pending[0]['rid'], pending[0]['reserve'], reserved, free))
                last_wait = time.time()
        if not started:
            time.sleep(1.0 if running else 10)
    log('done: %s' % ', '.join('%s %d' % kv for kv in sorted(counts.items(), key=lambda kv: str(kv[0]))))
    return 0


# ---------------------------------------------------------------- CLI

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    def common(p):
        p.add_argument('--timeout', type=float, default=DEFAULT_TIMEOUT,
                       help='summed abPOA wall clock per region (s, default %d)' % DEFAULT_TIMEOUT)
        p.add_argument('--mem-mb', type=float, default=DEFAULT_MEM_MB,
                       help='predicted or sampled abPOA RSS cap (MB, default %d)' % DEFAULT_MEM_MB)
        p.add_argument('--cand-root', help='candidates root (default %s)' % config.CANDIDATES_DIR)
        p.add_argument('--panel-root', help='full-panel outputs root (default %s)' % panel.PANEL_DIR)
        p.add_argument('--workdir', help='abPOA scratch (default $VNTR_WORK/mcpoa/tmp)')
        p.add_argument('--no-runtime', action='store_true', help='do not write results/mcpoa_runtime.all.tsv')
        p.add_argument('--max-wait', type=float, default=3600, help='seconds a job may wait for memory')
    p = sub.add_parser('run', help='a batch of regions')
    p.add_argument('--regions', default='all')
    p.add_argument('--stratum')
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--budget-mb', type=float, default=DEFAULT_BUDGET_MB)
    p.add_argument('--force', action='store_true')
    p.add_argument('--retry-failed', action='store_true')
    common(p)
    p = sub.add_parser('one', help='one region in this process')
    p.add_argument('region')
    p.add_argument('--result', help='also write the runtime row here (JSON)')
    p.add_argument('--no-wait-memory', dest='wait_memory', action='store_false')
    common(p)
    p = sub.add_parser('plan', help='records, BAR rule and predicted memory per region')
    p.add_argument('--regions', default='all')
    p.add_argument('--stratum')
    a = ap.parse_args(argv)

    if a.cmd == 'plan':
        cols = ['region_id', 'stratum', 'n_records', 'n_distinct', 'n_hap32_only', 'max_interior', 'len2', 'bar_rule',
                'rows_ge_10kb', 'rows_gt_20kb', 'unaligned_bp_planned', 'aligned_bp', 'predicted_mb']
        w = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
        w.writerow(cols)
        for rid in poa_panel.pick_regions(a.regions, a.stratum):
            p = plan_region(rid)
            w.writerow([p[c] for c in cols])
        return 0
    if a.cmd == 'one':
        row = run_job(a.region, timeout=a.timeout, mem_mb=a.mem_mb, cand_root=a.cand_root, panel_root=a.panel_root,
                      workroot=a.workdir, record_runtime=not a.no_runtime, wait_memory=a.wait_memory,
                      max_wait=a.max_wait)
        if a.result:
            with open(a.result + '.tmp', 'w') as f:
                json.dump(row, f)
            os.replace(a.result + '.tmp', a.result)
        log('%s: %s in %.1f s %s' % (a.region, row['status'], row['seconds'], row.get('note', '')))
        return 0 if row['status'] == 'ok' else (4 if row['status'] in ('timeout', 'memout') else 2)
    return run_batch(a)


if __name__ == '__main__':
    sys.exit(main())
