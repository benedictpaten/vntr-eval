#!/usr/bin/env python3
"""iterate.py -- Stage 4h: iterate on the full-panel abPOA arm over a fixed chr20 test set.

A VARIANT is a way to build a hap32 candidate graph for a region: MC itself, an existing Stage 4
candidate set (poa_abpoa, poa_abpoa__all), an abPOA flag string / input order run on the region's
full HG002-free panel (tools/panel.py union) and projected onto the 34 hap32 rows, or a named
structural variant (seeded incremental alignment, post-alignment node merging). For every region of
the test set (work/iterate/testset.tsv) the harness

  build   -> work/iterate/candidates/<variant>/<id>.{gfa,msa.fa,json}   (+ full MSA in panel/<variant>/)
  stage0  -> work/iterate/stage0/<variant>/<id>.json   (evaluate.py --skip truth: kmer_frac_extra,
             all_cost_over_opt, nodes/kb)
  call    -> local vg call on the hybrid graph with re-mapped reads, 3 replicates (call_local.py:
             hybrid200k, hybrid50k, hybrid200kids), under work/iterate/w/stage3/
  score   -> work/iterate/score/<variant>[@arm]/<id>.json (score_haplotypes.score_region, no truvari;
             the score is sensitivity.ed_suppress_nested, as tools/stage3_rule.py)
  clean   -> drops the hybrid graphs, indexes and GAFs of scored calls (VCFs and build.json kept)
  table   -> work/iterate/results.tsv (one row per region x variant) and a per-variant summary
             (median kmer_frac_extra and cost/opt, summed median ED, replicate-rule better/worse
             against mc, poa_abpoa and poa_abpoa__all)

Replicate rule (tools/stage3_rule.py): a variant is BETTER than a baseline at a region only if every
one of its replicate EDs is below every one of the baseline's, WORSE only if every one is above;
otherwise a tie. Both need >= 2 replicates.

Usage (from the repository root):

    python3 tools/iterate.py list
    python3 tools/iterate.py build  fp_G [--jobs 4]              # a registered variant
    python3 tools/iterate.py build  myvar --flags '-X 6 -O 6,26' [--order longest|shortest|freq]
    python3 tools/iterate.py stage0 fp_G,fp_R [--jobs 4]
    python3 tools/iterate.py call   fp_G [--jobs 4]
    python3 tools/iterate.py score  fp_G ; python3 tools/iterate.py clean fp_G
    python3 tools/iterate.py all    fp_G                          # build, stage0, call, score, clean, table
    python3 tools/iterate.py table  [--variants a,b,...]

Environment: VNTR_ITERATE (default work/iterate), VNTR_REGIONS (default work/stage4/regions),
VNTR_STAGE2_WORK (the directory holding fraglen.json; default work/iterate/stage2).
"""
import argparse
import collections
import concurrent.futures
import csv
import glob
import gzip
import heapq
import json
import os
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

TOOLS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(TOOLS)
os.environ.setdefault('VNTR_REGIONS', os.path.join(REPO, 'work', 'stage4', 'regions'))
sys.path.insert(0, TOOLS)
import config  # noqa: E402
import msa_graph  # noqa: E402
import panel  # noqa: E402
import poa_panel  # noqa: E402
import realign  # noqa: E402
import realign_poa  # noqa: E402

IT = os.environ.get('VNTR_ITERATE') or os.path.join(REPO, 'work', 'iterate')
REGIONS = config.REGIONS_DIR
TESTSET = os.path.join(IT, 'testset.tsv')
CAND = os.path.join(IT, 'candidates')
FULLMSA = os.path.join(IT, 'panel')
STAGE0 = os.path.join(IT, 'stage0')
SCORE = os.path.join(IT, 'score')
CALLWORK = os.path.join(IT, 'w')                       # VNTR_WORK of call_local.py
STAGE3 = os.path.join(CALLWORK, 'stage3')
STAGE2 = os.environ.get('VNTR_STAGE2_WORK') or os.path.join(IT, 'stage2')
RESULTS = os.path.join(IT, 'results.tsv')
SUMMARY = os.path.join(IT, 'summary.tsv')
TMP = os.path.join(IT, 'tmp')
S4 = os.path.join(REPO, 'work', 'stage4')
VG = os.path.join(REPO, 'work', 'bin', 'vg-2a6a228a5')
ABPOA = realign_poa.ABPOA
ARMS = ['hybrid50k', 'hybrid200kids']                  # replicates besides the Stage 3 call (hybrid200k)
REPS = [''] + ['@' + a for a in ARMS]
BASELINES = ['mc', 'poa_abpoa', 'poa_abpoa__all']
TIMEOUT = 900
MEM_MB = 12000

# the full-panel memory model (poa_panel), not the hap32 arm's len1^2 refusal
realign_poa.predict_mb = poa_panel._aligner_predict

# ------------------------------------------------------------------ variants

VARIANTS = collections.OrderedDict()


def _reg(name, desc, **kw):
    kw['desc'] = desc
    VARIANTS[name] = kw


_reg('mc', 'Minigraph-Cactus (regions/<id>/mc.gfa)', kind='mc')
_reg('poa_abpoa', 'abPOA on the 34 hap32 rows (Stage 4 candidates)', kind='link',
     src=os.path.join(S4, 'candidates', 'poa_abpoa'), full=None)
_reg('poa_abpoa__all', 'abPOA defaults on the full panel, projected (Stage 4c candidates)', kind='link',
     src=os.path.join(S4, 'candidates', 'poa_abpoa__all'), full=os.path.join(S4, 'panel', 'poa_abpoa'))
_reg('fp_default', 'full panel, abPOA defaults, longest first (rebuild of poa_abpoa__all)', kind='abpoa', flags=[])
_reg('fp_G', 'full panel, -G (log path score)', kind='abpoa', flags=['-G'])
_reg('fp_R', 'full panel, -R (indels right-most)', kind='abpoa', flags=['-R'])
_reg('fp_J', 'full panel, -J (indels at the end)', kind='abpoa', flags=['-J'])
_reg('fp_RJ', 'full panel, -R -J', kind='abpoa', flags=['-R', '-J'])
_reg('fp_nb', 'full panel, unbanded (-b -1)', kind='abpoa', flags=['-b', '-1'])
_reg('fp_b100', 'full panel, wide band (-b 100 -f 0.1)', kind='abpoa', flags=['-b', '100', '-f', '0.1'])
_reg('fp_X6', 'full panel, mismatch 6 (-X 6)', kind='abpoa', flags=['-X', '6'])
_reg('fp_O6', 'full panel, convex gaps -O 6,26 -E 2,1', kind='abpoa', flags=['-O', '6,26', '-E', '2,1'])
_reg('fp_aff', 'full panel, affine gaps -O 6,0 -E 2,0', kind='abpoa', flags=['-O', '6,0', '-E', '2,0'])
_reg('fp_short', 'full panel, defaults, shortest first', kind='abpoa', flags=[], order='shortest')
_reg('fp_freq', 'full panel, defaults, most frequent distinct allele first (panel carriers)', kind='abpoa',
     flags=[], order='freq')
_reg('s1_seed10', 'seed: 10 most frequent distinct alleles + longest + shortest by abPOA, the rest by abpoa -i',
     kind='abpoa', flags=[], seed_k=10)
_reg('s1_seed20', 'seed: 20 most frequent distinct alleles + longest + shortest by abPOA, the rest by abpoa -i',
     kind='abpoa', flags=[], seed_k=20)
_reg('s2_modn_pre', 'poa_abpoa__all full-panel graph, vg mod -U 10 (normalize), then restricted to hap32',
     kind='merge', base='poa_abpoa__all', where='pre', how='modn')
_reg('s2_modn_post', 'poa_abpoa__all projected graph, vg mod -U 10 (normalize)',
     kind='merge', base='poa_abpoa__all', where='post', how='modn')
_reg('s2_sib_pre', 'poa_abpoa__all full-panel graph, merge sibling nodes with identical sequence, predecessors '
     'and successors, then restricted to hap32', kind='merge', base='poa_abpoa__all', where='pre', how='sib')
_reg('s2_zip_pre', 'poa_abpoa__all full-panel graph, zip: merge nodes with identical sequence and identical '
     'predecessors or identical successors, to a fixed point, then restricted to hap32',
     kind='merge', base='poa_abpoa__all', where='pre', how='zip')
_reg('s2_zip_post', 'poa_abpoa__all projected graph, zip merge (as s2_zip_pre) on the hap32 graph',
     kind='merge', base='poa_abpoa__all', where='post', how='zip')


def variant(name):
    if name in VARIANTS:
        return VARIANTS[name]
    p = os.path.join(CAND, name, 'variant.json')
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    raise SystemExit('unknown variant %s (python3 tools/iterate.py list)' % name)


def log(*a):
    sys.stderr.write('[%s] %s\n' % (time.strftime('%H:%M:%S'), ' '.join(str(x) for x in a)))
    sys.stderr.flush()


def testset():
    rows = []
    with open(TESTSET) as f:
        for r in csv.DictReader((x for x in f if not x.startswith('#')), delimiter='\t'):
            rows.append(r)
    return rows


def region_ids(spec=None):
    ids = [r['region_id'] for r in testset()]
    if spec:
        want = spec.split(',')
        bad = [x for x in want if x not in ids]
        if bad:
            log('note: regions outside the test set: %s' % ','.join(bad))
        return want
    return ids


def cand_gfa(v, rid):
    if v == 'mc':
        return os.path.join(REGIONS, rid, 'mc.gfa')
    return os.path.join(CAND, v, rid + '.gfa')


def full_msa_of(v, rid):
    """The full-panel MSA (rows named by distinct id) a variant was projected from, or None."""
    spec = variant(v)
    if spec['kind'] == 'link':
        return os.path.join(spec['full'], rid + '.msa.fa.gz') if spec.get('full') else None
    if spec['kind'] == 'abpoa':
        return os.path.join(FULLMSA, v, rid + '.msa.fa.gz')
    return None


# ------------------------------------------------------------------ build: abPOA variants

def _ordered(recs, order, weight):
    k = lambda n: int(n[1:])   # noqa: E731  (realign.align_fasta names the aligner input s<k>)
    if order == 'shortest':
        key = lambda r: (len(r[1]), k(r[0]))   # noqa: E731
    elif order == 'freq':
        key = lambda r: (-weight[k(r[0])], -len(r[1]), k(r[0]))   # noqa: E731
    else:
        key = lambda r: (-len(r[1]), k(r[0]))   # noqa: E731
    return sorted(recs, key=key)


def _abpoa(flags, recs, out_fa, wd, timeout, mem_mb, incr=None):
    os.makedirs(wd, exist_ok=True)
    fa = os.path.join(wd, 'in.fa')
    realign._write_fa(fa, recs)
    cmd = [ABPOA, '-m', '0', '-r', '1'] + list(flags) + (['-i', incr] if incr else [])
    return realign_poa._run('abpoa', cmd, fa, out_fa, wd, timeout, mem_mb, 'given', {'flags': list(flags)})


def variant_align(in_fa, out_fa, threads=1, workdir=None, timeout=TIMEOUT, mem_mb=MEM_MB, flags=(),
                  order='longest', weight=None, seed_k=None, **_):
    """realign.py plugin: abpoa -m 0 -r 1 FLAGS on the distinct masked sequences in ORDER; with seed_k,
    first the seed (the seed_k heaviest distinct alleles by panel carriers, plus the longest and the
    shortest; longest first), then every other sequence (longest first) added with abpoa -i."""
    recs = msa_graph.read_fasta(in_fa)
    if not seed_k:
        return _abpoa(flags, _ordered(recs, order, weight), out_fa, os.path.join(workdir, 'one'), timeout, mem_mb)
    byw = _ordered(recs, 'freq', weight)
    seed = {n for n, _ in byw[:seed_k]}
    bylen = _ordered(recs, 'longest', weight)
    seed |= {bylen[0][0], bylen[-1][0]}
    A = [r for r in bylen if r[0] in seed]
    B = [r for r in bylen if r[0] not in seed]
    wa = os.path.join(workdir, 'seed')
    amsa = os.path.join(wa, 'seed.msa.fa')
    r1 = _abpoa(flags, A, amsa, wa, timeout, mem_mb)
    if r1['status'] != 'ok' or not B:
        if r1['status'] == 'ok':
            shutil.copy(amsa, out_fa)
        return r1
    r2 = _abpoa(flags, B, out_fa, os.path.join(workdir, 'rest'), timeout, mem_mb, incr=amsa)
    r2['stage1'] = {k: r1.get(k) for k in ('status', 'seconds', 'peak_rss_mb')}
    r2['seconds'] = round((r1.get('seconds') or 0) + (r2.get('seconds') or 0), 2)
    r2['peak_rss_mb'] = max(r1.get('peak_rss_mb') or 0, r2.get('peak_rss_mb') or 0)
    r2['n_seed'], r2['n_added'] = len(A), len(B)
    return r2


def finish_graph(rid, proj_msa, out_gfa, wd):
    hap_fa = os.path.join(REGIONS, rid, 'hap32.fa')
    st = msa_graph.msa_to_gfa(proj_msa, hap_fa, out_gfa + '.tmp', workdir=wd, engine='native')
    poa_panel.check_paths(out_gfa + '.tmp', msa_graph.read_fasta(hap_fa))
    os.replace(out_gfa + '.tmp', out_gfa)
    return {k: st.get(k) for k in ('nodes', 'edges', 'columns')}


def build_abpoa(v, spec, rid, timeout, mem_mb):
    out_dir = os.path.join(CAND, v)
    rd = os.path.join(REGIONS, rid)
    fa, mp = poa_panel.ensure_union(rd, rid)
    rows = panel.read_map(mp)
    weight = [int(r['n_panel_samples']) for r in rows]
    lens = sorted((int(r['length']) for r in rows), reverse=True)
    info = {'variant': v, 'region_id': rid, 'n_distinct': len(rows), 'flags': spec.get('flags', []),
            'order': spec.get('order', 'longest'), 'seed_k': spec.get('seed_k'),
            'predicted_mb': poa_panel.predict_mb('abpoa', lens)}
    wd = tempfile.mkdtemp(prefix='it.%s.%s.' % (v, rid), dir=TMP)
    try:
        if info['predicted_mb'] > mem_mb:
            info.update(status='memout', message='predicted %.0f MB' % info['predicted_mb'])
            return info
        m = realign._as_method(v, {'align': variant_align, 'tool': 'abpoa', 'description': spec['desc'],
                                   'params': {'flags': spec.get('flags', []), 'order': spec.get('order', 'longest'),
                                              'weight': weight, 'seed_k': spec.get('seed_k')}})
        full = os.path.join(wd, 'full.msa.fa')
        al = realign.align_fasta(m, fa, full, threads=1, timeout=timeout, mem_mb=mem_mb, workdir=wd, dedup=True)
        a = al.get('aligner') or {}
        info['align'] = {'status': al['status'], 'message': al.get('message'), 'seconds': a.get('seconds'),
                         'peak_rss_mb': a.get('peak_rss_mb'), 'command': a.get('command'),
                         'stage1': a.get('stage1'), 'n_seed': a.get('n_seed')}
        if al['status'] != 'ok':
            info.update(status=al['status'], message=al.get('message'))
            return info
        os.makedirs(os.path.join(FULLMSA, v), exist_ok=True)
        with open(full, 'rb') as fi, gzip.open(os.path.join(FULLMSA, v, rid + '.msa.fa.gz'), 'wb') as fo:
            shutil.copyfileobj(fi, fo)
        proj = os.path.join(wd, 'proj.msa.fa')
        info['projection'] = panel.project(full, mp, os.path.join(rd, 'hap32.fa'), proj)
        info['graph'] = finish_graph(rid, proj, os.path.join(out_dir, rid + '.gfa'), wd)
        shutil.copy(proj, os.path.join(out_dir, rid + '.msa.fa'))
        info['status'] = 'ok'
        return info
    finally:
        shutil.rmtree(wd, ignore_errors=True)


# ------------------------------------------------------------------ build: node-merging variants

def read_gfa(path):
    """(seqs {id: seq}, paths {name: [(id, '+'|'-')]}) from S and P lines."""
    seqs, paths = {}, collections.OrderedDict()
    with open(path) as f:
        for line in f:
            if line.startswith('S\t'):
                x = line.rstrip('\n').split('\t')
                seqs[x[1]] = x[2].upper()
            elif line.startswith('P\t'):
                x = line.rstrip('\n').split('\t')
                paths[x[1]] = [(s[:-1], s[-1]) for s in x[2].split(',') if s]
    return seqs, paths


def spell(seqs, steps):
    return ''.join(seqs[n] if o == '+' else panel.revcomp(seqs[n]) for n, o in steps)


def unchop(seqs, paths):
    """Merge every non-branching run (a -> b where a's only successor is b, b's only predecessor is a,
    no path ends at a and none starts at b), as msa_graph.py's column blocks are; paths are rewritten."""
    succ, pred = collections.defaultdict(set), collections.defaultdict(set)
    for steps in paths.values():
        for (a, _), (b, _) in zip(steps, steps[1:]):
            succ[a].add(b)
            pred[b].add(a)
        if steps:
            pred[steps[0][0]].add('^')
            succ[steps[-1][0]].add('$')
    nxt = {a: next(iter(s)) for a, s in succ.items() if len(s) == 1 and '$' not in s}
    nxt = {a: b for a, b in nxt.items() if pred[b] == {a} and b != a}
    heads = [n for n in nxt if n not in set(nxt.values())]
    chain_of, seqs2 = {}, dict(seqs)
    for h in heads:
        ch = [h]
        while ch[-1] in nxt:
            ch.append(nxt[ch[-1]])
        chain_of[h] = ch
        seqs2[h] = ''.join(seqs[x] for x in ch)
    inner = {x for ch in chain_of.values() for x in ch[1:]}
    out = collections.OrderedDict()
    for name, steps in paths.items():
        out[name] = [(n, o) for n, o in steps if n not in inner]
    return seqs2, out


def write_gfa(path, seqs, paths):
    """Unchop, keep only the nodes and edges the paths use, renumber 1.. in a topological order of those
    edges (call_local.py lays out span node IDs in ID order). Raises on a cycle."""
    seqs, paths = unchop(seqs, paths)
    used = collections.OrderedDict()
    edges = set()
    for steps in paths.values():
        for n, o in steps:
            if o != '+':
                raise ValueError('reverse step %s%s' % (n, o))
            used.setdefault(n, None)
        for (a, _), (b, _) in zip(steps, steps[1:]):
            edges.add((a, b))
    succ, indeg = collections.defaultdict(list), {n: 0 for n in used}
    for a, b in edges:
        succ[a].append(b)
        indeg[b] += 1
    first = {}
    for steps in paths.values():            # tie-break: first appearance along the paths
        for i, (n, _) in enumerate(steps):
            first[n] = min(first.get(n, 1 << 60), i * 100000 + len(first))
    heap = [(first[n], n) for n, d in indeg.items() if d == 0]
    heapq.heapify(heap)
    order = []
    while heap:
        _, n = heapq.heappop(heap)
        order.append(n)
        for m in succ[n]:
            indeg[m] -= 1
            if indeg[m] == 0:
                heapq.heappush(heap, (first[m], m))
    if len(order) != len(used):
        raise ValueError('cycle: %d of %d nodes ordered' % (len(order), len(used)))
    new = {n: str(i + 1) for i, n in enumerate(order)}
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        f.write('H\tVN:Z:1.0\n')
        for n in order:
            f.write('S\t%s\t%s\n' % (new[n], seqs[n]))
        for a, b in sorted(edges, key=lambda e: (int(new[e[0]]), int(new[e[1]]))):
            f.write('L\t%s\t+\t%s\t+\t0M\n' % (new[a], new[b]))
        for name, steps in paths.items():
            f.write('P\t%s\t%s\t*\n' % (name, ','.join(new[n] + '+' for n, _ in steps)))
    os.replace(tmp, path)
    return {'nodes': len(order), 'edges': len(edges)}


def merge_sib(seqs, paths, mode):
    """Merge forward nodes with identical sequence and identical predecessor AND successor sets
    (mode 'sib'), or identical predecessor OR identical successor sets (mode 'zip'), repeated until
    nothing merges. Merging two nodes with the same predecessors (or successors) cannot close a cycle
    in a DAG; paths are relabelled, so they spell the same sequences."""
    paths = {k: list(v) for k, v in paths.items()}
    rounds = merged = 0
    while True:
        rounds += 1
        pred, succ = collections.defaultdict(set), collections.defaultdict(set)
        for steps in paths.values():
            for (a, _), (b, _) in zip(steps, steps[1:]):
                succ[a].add(b)
                pred[b].add(a)
            if steps:
                pred[steps[0][0]].add('^')
                succ[steps[-1][0]].add('$')
        nodes = {n for steps in paths.values() for n, _ in steps}
        keys = [lambda n: (seqs[n], frozenset(pred[n]), frozenset(succ[n]))] if mode == 'sib' else \
            [lambda n: (seqs[n], frozenset(pred[n])), lambda n: (seqs[n], frozenset(succ[n]))]
        remap = {}
        for key in keys:
            groups = collections.defaultdict(list)
            for n in sorted(nodes, key=lambda x: int(x)):
                if n not in remap:
                    groups[key(n)].append(n)
            for g in groups.values():
                for n in g[1:]:
                    remap[n] = g[0]
            if remap:
                break                        # re-derive neighbours before the other key
        if not remap:
            return paths, {'rounds': rounds, 'merged': merged}
        merged += len(remap)
        paths = {k: [(remap.get(n, n), o) for n, o in v] for k, v in paths.items()}


def merge_modn(seqs, paths, wd):
    """vg mod -U 10 (iterated normalize) on a GFA whose paths are renamed p<i>; back to (seqs, paths)."""
    names = list(paths)
    gin, gvg, gout = (os.path.join(wd, x) for x in ('in.gfa', 'mod.vg', 'mod.gfa'))
    ren = collections.OrderedDict(('p%d' % i, paths[n]) for i, n in enumerate(names))
    nodes = sorted({n for s in ren.values() for n, _ in s}, key=int)
    edges = sorted({(a, b) for s in ren.values() for (a, _), (b, _) in zip(s, s[1:])}, key=lambda e: (int(e[0]), int(e[1])))
    with open(gin, 'w') as f:
        f.write('H\tVN:Z:1.0\n')
        for n in nodes:
            f.write('S\t%s\t%s\n' % (n, seqs[n]))
        for a, b in edges:
            f.write('L\t%s\t+\t%s\t+\t0M\n' % (a, b))
        for k, s in ren.items():
            f.write('P\t%s\t%s\t*\n' % (k, ','.join(n + o for n, o in s)))
    env = config.tool_env()
    with open(gvg, 'wb') as fo:
        subprocess.run([VG, 'mod', '-U', '10', gin], stdout=fo, check=True, env=env, stderr=subprocess.DEVNULL)
    with open(gout, 'w') as fo:
        subprocess.run([VG, 'convert', '-f', gvg], stdout=fo, check=True, env=env, stderr=subprocess.DEVNULL)
    s2, p2 = read_gfa(gout)
    back = collections.OrderedDict((names[int(k[1:])], p2[k]) for k in ren)
    return s2, back, {}


def build_merge(v, spec, rid):
    rd = os.path.join(REGIONS, rid)
    hap = msa_graph.read_fasta(os.path.join(rd, 'hap32.fa'))
    info = {'variant': v, 'region_id': rid, 'base': spec['base'], 'where': spec['where'], 'how': spec['how']}
    wd = tempfile.mkdtemp(prefix='it.%s.%s.' % (v, rid), dir=TMP)
    try:
        if spec['where'] == 'post':
            src = cand_gfa(spec['base'], rid)
            if not os.path.exists(src):
                info.update(status='no_base', message=src)
                return info
            seqs, paths = read_gfa(src)
            upaths = paths
        else:
            fmsa = full_msa_of(spec['base'], rid)
            if not fmsa or not os.path.exists(fmsa):
                info.update(status='no_base', message=str(fmsa))
                return info
            fa, mp = poa_panel.union_files(rid)
            msa = os.path.join(wd, 'full.msa.fa')
            with gzip.open(fmsa, 'rb') as fi, open(msa, 'wb') as fo:
                shutil.copyfileobj(fi, fo)
            full_gfa = os.path.join(wd, 'full.gfa')
            msa_graph.msa_to_gfa(msa, fa, full_gfa, workdir=wd, engine='native', ignore_extra=True)
            seqs, upaths = read_gfa(full_gfa)
            info['full_nodes_before'] = len(seqs)
        t0 = time.time()
        if spec['how'] == 'modn':
            seqs2, up2, st = merge_modn(seqs, upaths, wd)
        else:
            up2, st = merge_sib(seqs, upaths, spec['how'])
            seqs2 = seqs
        info['merge'] = dict(st, seconds=round(time.time() - t0, 2))
        if spec['where'] == 'post':
            out_paths = up2
        else:
            uid = {}
            for r in panel.read_map(mp):
                for m in r['members']:
                    uid[m] = r['id']
            out_paths = collections.OrderedDict((n, up2[uid[n]]) for n, _ in hap)
            info['full_nodes_after'] = len({x for s in up2.values() for x, _ in s})
        bad = [n for n, s in hap if spell(seqs2, out_paths[n]) != s.upper()]
        if bad:
            raise ValueError('%d paths do not spell hap32.fa after merging: %s' % (len(bad), bad[:3]))
        gfa = os.path.join(CAND, v, rid + '.gfa')
        info['graph'] = write_gfa(gfa, seqs2, out_paths)
        poa_panel.check_paths(gfa, hap)
        info['status'] = 'ok'
        return info
    finally:
        shutil.rmtree(wd, ignore_errors=True)


def build_one(v, rid, timeout=TIMEOUT, mem_mb=MEM_MB, force=False):
    spec = variant(v)
    out_dir = os.path.join(CAND, v)
    js = os.path.join(out_dir, rid + '.json')
    if os.path.exists(js) and not force:
        with open(js) as f:
            return json.load(f)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(TMP, exist_ok=True)
    t0 = time.time()
    try:
        if spec['kind'] == 'mc':
            return {'variant': v, 'region_id': rid, 'status': 'ok'}
        if spec['kind'] == 'link':
            src = os.path.join(spec['src'], rid + '.gfa')
            if os.path.exists(src):
                shutil.copy(src, os.path.join(out_dir, rid + '.gfa'))
                info = {'variant': v, 'region_id': rid, 'status': 'ok', 'src': src}
                rj = os.path.join(spec['src'], rid + '.realign.json')
                if os.path.exists(rj):
                    a = (json.load(open(rj)).get('align') or {}).get('aligner') or {}
                    info['align'] = {'seconds': a.get('seconds'), 'peak_rss_mb': a.get('peak_rss_mb')}
            else:
                info = {'variant': v, 'region_id': rid, 'status': 'missing', 'src': src}
        elif spec['kind'] == 'abpoa':
            info = build_abpoa(v, spec, rid, timeout, mem_mb)
        elif spec['kind'] == 'merge':
            info = build_merge(v, spec, rid)
        else:
            raise ValueError('unknown kind %s' % spec['kind'])
    except Exception as e:  # noqa: BLE001  (one region must not end a batch)
        info = {'variant': v, 'region_id': rid, 'status': 'error', 'message': repr(e)[:500]}
    info['seconds'] = round(time.time() - t0, 1)
    with open(js, 'w') as f:
        json.dump(info, f, indent=1, default=str)
    return info


def _pool(fn, jobs, items):
    with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(fn, *it): it for it in items}
        for fu in concurrent.futures.as_completed(futs):
            yield futs[fu], fu.result()


def cmd_build(variants, ids, jobs, force, timeout, mem_mb):
    for v in variants:
        spec = variant(v)
        os.makedirs(os.path.join(CAND, v), exist_ok=True)
        with open(os.path.join(CAND, v, 'variant.json'), 'w') as f:
            json.dump(dict(spec, name=v), f, indent=1)
        if spec['kind'] == 'merge':
            cmd_build([spec['base']], ids, jobs, False, timeout, mem_mb)
        # biggest first so the long jobs start early
        size = {r: os.path.getsize(poa_panel.union_files(r)[0]) if os.path.exists(poa_panel.union_files(r)[0]) else 0
                for r in ids}
        items = [(v, r, timeout, mem_mb, force) for r in sorted(ids, key=lambda r: -size[r])]
        for (_, rid, *_), info in _pool(build_one, jobs, items):
            if info.get('status') != 'ok' or force:
                log('build', v, rid, info.get('status'), info.get('seconds'), info.get('message', '') or '')


# ------------------------------------------------------------------ stage 0

def stage0_one(v, rid, force=False):
    out = os.path.join(STAGE0, v, rid + '.json')
    if os.path.exists(out) and not force:
        return 'cached'
    gfa = cand_gfa(v, rid)
    if not os.path.exists(gfa):
        return 'no_graph'
    os.makedirs(os.path.dirname(out), exist_ok=True)
    lg = out.replace('.json', '.log')
    with open(lg, 'w') as f:
        p = subprocess.run([sys.executable, os.path.join(TOOLS, 'evaluate.py'), os.path.join(REGIONS, rid), gfa,
                            '--name', 'it_' + v, '--out', out, '--skip', 'truth', '--threads', '2'],
                           stdout=f, stderr=subprocess.STDOUT, env=dict(os.environ))
    return 'ok' if p.returncode == 0 else 'error rc=%d' % p.returncode


def cmd_stage0(variants, ids, jobs, force):
    items = [(v, r, force) for v in variants for r in ids]
    for (v, rid, _), st in _pool(stage0_one, jobs, items):
        if st not in ('ok', 'cached'):
            log('stage0', v, rid, st)


def stage0_metrics(v, rid):
    p = os.path.join(STAGE0, v, rid + '.json')
    if not os.path.exists(p):
        return {}
    d = json.load(open(p))
    return {'nodes': (d.get('size') or {}).get('nodes'), 'nodes_per_kb': (d.get('size') or {}).get('nodes_per_kb'),
            'kmer_frac_extra': (d.get('redundancy') or {}).get('kmer_frac_extra'),
            'cost_over_opt': (d.get('alignment') or {}).get('all_cost_over_opt')}


# ------------------------------------------------------------------ local calling and scoring

def call_env():
    env = dict(os.environ, VNTR_REGIONS=REGIONS, VNTR_CANDIDATES=CAND, VNTR_WORK=CALLWORK, VNTR_STAGE2_WORK=STAGE2)
    env.setdefault('TMPDIR', TMP)
    return env


def vcf_path(v, rid, rep):
    if rep == '':
        return os.path.join(STAGE3, 'calls', v, rid + '.vcf.gz')
    return os.path.join(STAGE3, 'diag', v + rep, rid + '.vcf.gz')


def cmd_call(variants, ids, jobs, threads):
    if not os.path.exists(os.path.join(STAGE2, 'fraglen.json')):
        raise SystemExit('no fraglen.json in %s (set VNTR_STAGE2_WORK)' % STAGE2)
    os.makedirs(TMP, exist_ok=True)
    for v in variants:
        have = [r for r in ids if os.path.exists(cand_gfa(v, r))]
        todo = [r for r in have if not all(os.path.exists(vcf_path(v, r, rep)) for rep in REPS)]
        if not todo:
            continue
        log('call', v, len(todo), 'regions x 3 replicates')
        lg = os.path.join(IT, 'logs', 'call.%s.log' % v)
        os.makedirs(os.path.dirname(lg), exist_ok=True)
        with open(lg, 'a') as f:
            subprocess.run([sys.executable, os.path.join(TOOLS, 'call_local.py'), 'run', ','.join(todo), '--graphs', v,
                            '--arms', ','.join(ARMS), '--arm-graphs', v, '--jobs', str(jobs), '--threads',
                            str(threads), '--compact'], stdout=f, stderr=subprocess.STDOUT, env=call_env(), cwd=REPO)
        miss = [(r, rep) for r in todo for rep in REPS if not os.path.exists(vcf_path(v, r, rep))]
        if miss:
            log('call', v, 'missing', len(miss), 'calls, e.g.', miss[:4], '(see %s)' % lg)


def score_one(v, rid, rep, force=False):
    import score_haplotypes as sh
    out = os.path.join(SCORE, v + rep, rid + '.json')
    if os.path.exists(out) and not force:
        return 'cached'
    vcf = vcf_path(v, rid, rep)
    if not os.path.exists(vcf):
        return 'no_vcf'
    os.makedirs(os.path.dirname(out), exist_ok=True)
    try:
        res = sh.score_region(os.path.join(REGIONS, rid), vcf, label=v + rep, truvari=False)
    except Exception as e:  # noqa: BLE001
        res = {'region_id': rid, 'status': 'error: %r' % (e,)}
    with open(out + '.tmp', 'w') as f:
        json.dump(res, f, indent=1)
    os.replace(out + '.tmp', out)
    return res.get('status')


def cmd_score(variants, ids, jobs, force):
    items = [(v, r, rep, force) for v in variants for r in ids for rep in REPS]
    for (v, rid, rep, _), st in _pool(score_one, jobs, items):
        if st not in ('ok', 'cached', 'no_vcf'):
            log('score', v + rep, rid, st)


def cmd_clean(variants, ids):
    """Drop what a finished call no longer needs: the hybrid graph, its GBZ, the GAF and the between-anchors
    GFA of every replicate whose score exists (build.json and the logs stay; the shared read fetch stays)."""
    freed = 0
    for v in variants:
        for rid in ids:
            if not all(os.path.exists(os.path.join(SCORE, v + rep, rid + '.json')) for rep in REPS):
                continue
            for d in glob.glob(os.path.join(STAGE3, 'hybrid*', v, rid)):
                for x in ('hybrid.gfa', 'between.gfa', 'call.gbz', 'reads.gaf', 'reads.gaf.gz'):
                    p = os.path.join(d, x)
                    if os.path.exists(p):
                        freed += os.path.getsize(p)
                        os.remove(p)
                shutil.rmtree(os.path.join(d, 'tmp'), ignore_errors=True)
                for p in glob.glob(os.path.join(d, 'idx*')):
                    freed += os.path.getsize(p)
                    os.remove(p)
    log('clean: freed %.1f MB' % (freed / 2 ** 20))


def rep_scores(v, rid):
    out = []
    for rep in REPS:
        p = os.path.join(SCORE, v + rep, rid + '.json')
        if os.path.exists(p):
            d = json.load(open(p))
            if d.get('status') in (None, 'ok'):
                out.append(d.get('sensitivity', {}).get('ed_suppress_nested', d['ed']))
    return out


def rule(a, b):
    if len(a) < 2 or len(b) < 2:
        return ''
    if max(a) < min(b):
        return 'better'
    if min(a) > max(b):
        return 'worse'
    return 'tie'


# ------------------------------------------------------------------ table

def known_variants():
    vs = list(BASELINES)
    for d in sorted(glob.glob(os.path.join(CAND, '*'))):
        n = os.path.basename(d)
        if n not in vs and os.path.isdir(d):
            vs.append(n)
    for d in sorted(glob.glob(os.path.join(STAGE0, '*'))):
        n = os.path.basename(d)
        if n not in vs:
            vs.append(n)
    return vs


def _med(xs):
    xs = [x for x in xs if x is not None]
    return statistics.median(xs) if xs else None


def cmd_table(variants=None, quiet=False):
    ts = testset()
    group = {r['region_id']: r['group'] for r in ts}
    ids = [r['region_id'] for r in ts]
    variants = variants or known_variants()
    cols = ['variant', 'region_id', 'group', 'build_status', 'align_s', 'peak_rss_mb', 'nodes', 'nodes_per_kb',
            'kmer_frac_extra', 'cost_over_opt', 'ed_reps', 'ed_median', 'vs_mc', 'vs_poa_abpoa', 'vs_poa_abpoa__all']
    base_scores = {b: {r: rep_scores(b, r) for r in ids} for b in BASELINES}
    rows = []
    for v in variants:
        for rid in ids:
            bj = os.path.join(CAND, v, rid + '.json')
            b = json.load(open(bj)) if os.path.exists(bj) else ({'status': 'ok'} if v == 'mc' else {})
            m = stage0_metrics(v, rid)
            sc = rep_scores(v, rid)
            al = b.get('align') or {}
            row = {'variant': v, 'region_id': rid, 'group': group[rid], 'build_status': b.get('status', ''),
                   'align_s': al.get('seconds'), 'peak_rss_mb': al.get('peak_rss_mb'),
                   'nodes': m.get('nodes'), 'nodes_per_kb': m.get('nodes_per_kb'),
                   'kmer_frac_extra': m.get('kmer_frac_extra'), 'cost_over_opt': m.get('cost_over_opt'),
                   'ed_reps': ','.join(str(x) for x in sc), 'ed_median': _med(sc) if sc else None}
            for bname in BASELINES:
                row['vs_' + bname] = '' if bname == v else rule(sc, base_scores[bname][rid])
            rows.append(row)
    with open(RESULTS, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join('' if r[c] is None else str(r[c]) for c in cols) + '\n')
    scols = ['variant', 'built', 'med_kmer_frac_extra', 'med_cost_over_opt', 'med_nodes_per_kb', 'max_align_s',
             'called', 'sum_ed_median', 'sum_ed_median_mc', 'sum_ed_median_poa', 'bw_vs_mc', 'bw_vs_poa_abpoa',
             'bw_vs_poa_abpoa__all']
    summ = []
    for v in variants:
        rs = [r for r in rows if r['variant'] == v]
        called = [r for r in rs if r['ed_median'] is not None]
        cid = {r['region_id'] for r in called}
        s = {'variant': v, 'built': sum(1 for r in rs if r['build_status'] == 'ok'),
             'med_kmer_frac_extra': _med([r['kmer_frac_extra'] for r in rs]),
             'med_cost_over_opt': _med([r['cost_over_opt'] for r in rs]),
             'med_nodes_per_kb': _med([r['nodes_per_kb'] for r in rs]),
             'max_align_s': max([r['align_s'] for r in rs if r['align_s'] is not None] or [None], key=lambda x: x or 0),
             'called': len(called), 'sum_ed_median': sum(r['ed_median'] for r in called) if called else None}
        for bname, key in (('mc', 'sum_ed_median_mc'), ('poa_abpoa', 'sum_ed_median_poa')):
            bs = [_med(base_scores[bname][r]) for r in cid if base_scores[bname][r]]
            s[key] = sum(bs) if called and len(bs) == len(cid) else None
        for bname in BASELINES:
            vals = [r['vs_' + bname] for r in called]
            s['bw_vs_' + bname] = '' if bname == v else '%d/%d' % (vals.count('better'), vals.count('worse'))
        summ.append(s)
    with open(SUMMARY, 'w') as f:
        f.write('\t'.join(scols) + '\n')
        for s in summ:
            f.write('\t'.join('' if s[c] is None else (('%.4g' % s[c]) if isinstance(s[c], float) else str(s[c]))
                              for c in scols) + '\n')
    if not quiet:
        fmt = lambda x: '-' if x is None else (('%.4f' % x) if isinstance(x, float) else str(x))   # noqa: E731
        print('%-16s %5s %8s %8s %7s %7s %6s %8s %8s %8s %7s %7s %7s' % (
            'variant', 'built', 'kmerx', 'cost/opt', 'nod/kb', 'maxs', 'called', 'sumED', 'ED(mc)', 'ED(poa)',
            'b/w mc', 'b/w poa', 'b/w all'))
        for s in summ:
            print('%-16s %5s %8s %8s %7s %7s %6s %8s %8s %8s %7s %7s %7s' % (
                s['variant'], s['built'], fmt(s['med_kmer_frac_extra']), fmt(s['med_cost_over_opt']),
                fmt(s['med_nodes_per_kb']), fmt(s['max_align_s']), s['called'], fmt(s['sum_ed_median']),
                fmt(s['sum_ed_median_mc']), fmt(s['sum_ed_median_poa']), s['bw_vs_mc'] or '-',
                s['bw_vs_poa_abpoa'] or '-', s['bw_vs_poa_abpoa__all'] or '-'))
    return rows, summ


# ------------------------------------------------------------------ main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('list')
    for c in ('build', 'stage0', 'call', 'score', 'clean', 'all'):
        p = sub.add_parser(c)
        p.add_argument('variants', help='comma list')
        p.add_argument('--regions', default=None, help='comma list (default: the test set)')
        p.add_argument('--jobs', type=int, default=4)
        p.add_argument('--threads', type=int, default=2, help='threads per local call')
        p.add_argument('--force', action='store_true')
        p.add_argument('--flags', default=None, help='define VARIANT (one) as abPOA with these flags')
        p.add_argument('--order', default='longest', choices=['longest', 'shortest', 'freq'])
        p.add_argument('--seed-k', type=int, default=None)
        p.add_argument('--timeout', type=float, default=TIMEOUT)
        p.add_argument('--mem-mb', type=float, default=MEM_MB)
    p = sub.add_parser('table')
    p.add_argument('--variants', default=None)
    a = ap.parse_args(argv)
    if a.cmd == 'list':
        for n, s in VARIANTS.items():
            print('%-16s %-6s %s' % (n, s['kind'], s['desc']))
        for d in sorted(glob.glob(os.path.join(CAND, '*', 'variant.json'))):
            n = os.path.basename(os.path.dirname(d))
            if n not in VARIANTS:
                s = json.load(open(d))
                print('%-16s %-6s %s' % (n, s['kind'], s['desc']))
        return
    if a.cmd == 'table':
        cmd_table(a.variants.split(',') if a.variants else None)
        return
    variants = a.variants.split(',')
    if a.flags is not None or a.seed_k:
        if len(variants) != 1 or variants[0] in VARIANTS:
            raise SystemExit('--flags/--seed-k define one new variant name')
        v = variants[0]
        flags = shlex.split(a.flags or '')
        VARIANTS[v] = {'kind': 'abpoa', 'flags': flags, 'order': a.order, 'seed_k': a.seed_k,
                       'desc': 'full panel, abPOA %s, order %s%s' % (' '.join(flags) or 'defaults', a.order,
                                                                   ', seed %d' % a.seed_k if a.seed_k else '')}
    ids = region_ids(a.regions)
    if a.cmd in ('build', 'all'):
        cmd_build(variants, ids, a.jobs, a.force, a.timeout, a.mem_mb)
    if a.cmd in ('stage0', 'all'):
        cmd_stage0(variants, ids, a.jobs, a.force)
    if a.cmd in ('call', 'all'):
        cmd_call(variants, ids, a.jobs, a.threads)
    if a.cmd in ('score', 'all'):
        cmd_score(variants, ids, a.jobs, a.force)
    if a.cmd in ('clean', 'all'):
        cmd_clean(variants, ids)
    if a.cmd == 'all':
        cmd_table()


if __name__ == '__main__':
    main()
