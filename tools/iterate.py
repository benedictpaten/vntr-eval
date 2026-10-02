#!/usr/bin/env python3
"""iterate.py -- Stage 4h: iterate on the full-panel abPOA arm over a fixed chr20 test set.

A VARIANT is a way to build a hap32 candidate graph for a region: MC itself, an existing Stage 4
candidate set (poa_abpoa, poa_abpoa__all), an abPOA flag string / input order run on the region's
full HG002-free panel (tools/panel.py union) and projected onto the 34 hap32 rows, or a named
structural variant (seeded incremental alignment, post-alignment node merging, coarsening, gap
normalisation of the full-panel MSA, a centre-star or profile-aligner MSA of the full panel,
nearest-neighbour threading of the full panel on a k-mer minimum spanning tree, or an abPOA backbone of
k-medoid representatives with the rest threaded on it). For
every region of the test set (work/iterate/testset.tsv) the harness

  build   -> work/iterate/candidates/<variant>/<id>.{gfa,msa.fa,json}   (+ full MSA in panel/<variant>/)
  stage0  -> work/iterate/stage0/<variant>/<id>.json   (evaluate.py --skip truth: kmer_frac_extra,
             all_cost_over_opt, nodes/kb) and <id>.snarls.json (top-level snarls = sites, all snarls)
  call    -> local vg call on the hybrid graph with re-mapped reads, 3 replicates (call_local.py:
             hybrid200k, hybrid50k, hybrid200kids), under work/iterate/w/stage3/
  score   -> work/iterate/score/<variant>[@arm]/<id>.json (score_haplotypes.score_region, no truvari;
             the score is sensitivity.ed_suppress_nested, as tools/stage3_rule.py)
  clean   -> drops the hybrid graphs, indexes and GAFs of scored calls (VCFs and build.json kept)
  table   -> work/iterate/results.tsv (one row per region x variant) and a per-variant summary
             (median and mean kmer_frac_extra, cost/opt and sites, summed median ED, replicate-rule better/worse
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
SNARL_VG = os.environ.get('VNTR_SNARL_VG') or os.path.join(REPO, 'work', 'bin', 'vg-91d38c802')
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


# Coarsened graphs: fewer, larger sites. Runs of >= k alignment columns where every hap32 row has the
# same base become shared anchor nodes; between two anchors each distinct row sequence is one allele node.
for _base in ('poa_abpoa', 'poa_abpoa__all', 'fp_G'):
    for _k in (16, 32, 64):
        _reg('cz%d_%s' % (_k, _base), 'coarsened %s: anchors = runs of >= %d identical columns, whole alleles '
             'between them' % (_base, _k), kind='coarsen', base=_base, k=_k)


# Gap normalisation (4m): the base's full-panel MSA with each row's gap runs shifted left (or right)
# against the column consensus to a fixed point (tools/gap_norm.py), then projected and rebuilt as
# the base was. 'merge' first joins a row's consecutive gap runs when moving the bases between them
# costs at most that many edits per other row (linear sum of pairs). where='proj' normalises the
# hap32 projection instead (against the hap32 rows' consensus): sample-dependent, a reference only.
for _base in ('poa_abpoa__all', 'fp_G'):
    for _d in ('left', 'right'):
        _reg('gn%s_%s' % (_d[0], _base), '%s full-panel MSA, gaps %s-normalised to the column consensus, '
             'projected' % (_base, _d), kind='gapnorm', base=_base, where='full', direction=_d)
_reg('gml_poa_abpoa__all', 'poa_abpoa__all full-panel MSA, gap runs merged (1 edit per other row), then '
     'left-normalised, projected', kind='gapnorm', base='poa_abpoa__all', where='full', direction='left', merge=1.0)
_reg('gjl_poa_abpoa__all', 'poa_abpoa__all full-panel MSA, gap runs left-normalised jointly (a run shared by '
     'several rows moves only if it can in all of them), projected', kind='gapnorm', base='poa_abpoa__all',
     where='full', direction='left', joint=True)
_reg('gnpl_poa_abpoa__all', 'poa_abpoa__all projection, gaps left-normalised to the hap32 consensus',
     kind='gapnorm', base='poa_abpoa__all', where='proj', direction='left')


# A different alignment objective for the full panel (4n). 'star': every distinct allele is aligned
# pairwise to one centre (abPOA on the two sequences, unbanded, its default two-piece affine scores),
# each pairwise alignment's indels are left-normalised, and the alignments are merged on the centre's
# columns; the insertions that fall between the same two centre bases are aligned to each other with
# abPOA. 'profile': a progressive aligner that scores a sequence or profile against a profile, whose
# position-specific gap costs make indels stack in shared columns. Both are projected onto hap32 as
# the abPOA variants are.
_reg('st_cons', 'full panel, centre-star to the abPOA consensus of the full panel (abpoa -r 0, longest '
     'first, as poa_abpoa__all), projected', kind='star', centre='consensus')
_reg('st_chm13', 'full panel, centre-star to CHM13\'s sequence, projected', kind='star', centre='ref')
_reg('st_long', 'full panel, centre-star to the longest distinct allele, projected', kind='star', centre='longest')
_reg('st_medoid', 'full panel, centre-star to the medoid allele (least panel-weighted k-mer distance to the rest), '
     'projected', kind='star', centre='medoid')
_reg('st_maj', 'full panel, centre-star to the majority consensus of the st_chm13 full-panel MSA (columns carried '
     'by at least half the panel weight), projected', kind='star', centre='majority', source='st_chm13', maj=0.5)
_reg('st_maj_med', 'full panel, centre-star to the majority consensus of the st_medoid full-panel MSA, so that no '
     'step depends on CHM13, projected', kind='star', centre='majority', source='st_medoid', maj=0.5)
_reg('pf_famsa', 'full panel, FAMSA defaults, projected', kind='profile', tool='famsa', args=[])
_reg('pf_famsa_go2', 'full panel, FAMSA with twice the gap-open cost (-go -29700), projected', kind='profile',
     tool='famsa', args=['-go', '-29700'])
# FAMSA gap-cost sweep (4t): its default gap costs are tuned for protein, and on whole chr20 its projection
# cuts single large deletions into many pieces. Refinement off keeps the higher gap-open runs inside the cap.
_reg('pf_famsa_r0', 'full panel, FAMSA with refinement off, projected', kind='profile', tool='famsa',
     args=['-refine_mode', 'off'])
_reg('pf_famsa_go2r0', 'full panel, FAMSA with twice the gap-open cost and refinement off, projected',
     kind='profile', tool='famsa', args=['-go', '-29700', '-refine_mode', 'off'])
_reg('pf_famsa_go4r0', 'full panel, FAMSA with four times the gap-open cost and refinement off, projected',
     kind='profile', tool='famsa', args=['-go', '-59400', '-refine_mode', 'off'])
_reg('pf_famsa_upgma', 'full panel, FAMSA with a UPGMA guide tree, projected', kind='profile', tool='famsa',
     args=['-gt', 'upgma'])
_reg('gnl_pf_famsa', 'pf_famsa full-panel MSA, gaps left-normalised per row, projected', kind='gapnorm',
     base='pf_famsa', where='full', direction='left')
_reg('gjl_pf_famsa', 'pf_famsa full-panel MSA, gap runs left-normalised jointly, projected', kind='gapnorm',
     base='pf_famsa', where='full', direction='left', joint=True)
_reg('pf_famsa_h32', 'FAMSA defaults on the hap32 rows alone (sample-dependent control)', kind='profile',
     tool='famsa', args=[], subset='hap32')
_reg('pf_famsa_ge2', 'full panel, FAMSA with half the gap-extension cost (-ge -625), projected', kind='profile',
     tool='famsa', args=['-ge', '-625'])
_reg('pf_famsa_ge4', 'full panel, FAMSA with a quarter of the gap-extension cost (-ge -312), projected',
     kind='profile', tool='famsa', args=['-ge', '-312'])
_reg('pf_famsa_ge8', 'full panel, FAMSA with an eighth of the gap-extension cost (-ge -156), projected',
     kind='profile', tool='famsa', args=['-ge', '-156'])
_reg('pf_kalign', 'full panel, Kalign 3 --type dna (defaults), projected', kind='profile', tool='kalign', args=[])
_reg('pf_muscle', 'full panel, MUSCLE5 -super5, projected', kind='profile', tool='muscle', args=[])
_reg('pf_mafft', 'full panel, mafft FFT-NS-2 (--retree 2 --maxiterate 0), projected', kind='profile',
     tool='mafft', args=['--retree', '2', '--maxiterate', '0'], threads=2)

# Nearest-neighbour threading (4q): each distinct allele is placed against ONE similar allele already in
# the MSA instead of a graph of all of them. The alleles join a minimum spanning tree over multiset
# k-mer distance in Prim's order from CHM13; each is aligned pairwise to its parent (abPOA, unbanded,
# indels left-normalised) and inherits the parent's columns, with its insertions placed in the slot's
# existing gap columns where they fit (so recurrent insertions stack) and new columns otherwise. mst3
# aligns each allele to its 3 nearest already-aligned alleles and keeps the lowest-cost alignment.
_reg('mst', 'full panel, nearest-neighbour threading on a k-mer MST from CHM13 (pairwise to the tree parent), '
     'projected', kind='mst', nn=1)
_reg('mst3', 'full panel, nearest-neighbour threading, best of the 3 nearest aligned alleles, projected',
     kind='mst', nn=3)

# Backbone plus threading (4r): K k-medoid representatives of the full panel's distinct sequences (multiset
# k-mer distance as mst; CHM13 and GRCh38 fixed as medoids; no use of HG002 or the hap32 rows) are aligned
# with abPOA defaults (as poa_abpoa), and every other sequence is threaded onto its nearest already-placed
# member of its cluster with mst's pairwise column-inheritance merge.
# bbt32m and bbt64m thread every member on its cluster's medoid instead.
for _k in (32, 64):
    _reg('bbt%d' % _k, 'full panel, abPOA backbone of %d k-medoid representatives plus nearest-neighbour '
         'threading within clusters, projected' % _k, kind='bbt', K=_k)
    _reg('bbt%dm' % _k, 'full panel, abPOA backbone of %d k-medoid representatives, every other allele threaded '
         'on its medoid, projected' % _k, kind='bbt', K=_k, attach='medoid')

# The repeat-unit-aware aligner (realign_units, 4o). ua32 is its hap32 arm (Stage 4 candidates);
# ua_all and ua_all_poa run it on the full panel (tools/units_panel.py run --panel-root
# work/stage4/panel [--fallback-engine abpoa]) and project the MSA onto hap32. They differ only in
# the aligner of the non-unit parts (flank pieces, regions without a usable motif, the guard's
# fallback): mafft, or abPOA.
_reg('ua32', 'unit_aware on the 34 hap32 rows (Stage 4 candidates)', kind='link',
     src=os.path.join(S4, 'candidates', 'unit_aware'), full=None)
_reg('ua_all', 'unit_aware on the full panel, mafft for non-unit parts, projected', kind='link',
     src=os.path.join(S4, 'candidates', 'unit_aware__all'), full=os.path.join(S4, 'panel', 'unit_aware'))
_reg('ua_all_poa', 'unit_aware on the full panel, abPOA for non-unit parts, projected', kind='link',
     src=os.path.join(S4, 'candidates', 'unit_aware_poa__all'), full=os.path.join(S4, 'panel', 'unit_aware_poa'))


# Linkage strength: the same graphs called with a stronger Li-Stephens linkage model (vg call
# --linkage-weight, default 2). lw2_mc is mc laid out as a candidate (like the lw*_mc arms) with the
# default weight, the baseline for them.
MC_LINK = os.path.join(IT, 'mc_link')
for _w in (2, 4, 8):
    _reg('lw%d_mc' % _w, 'mc graph as a candidate, --linkage-weight %d' % _w, kind='link', src=MC_LINK, full=None,
         extra=['--linkage-weight', str(_w)])
for _base in ('poa_abpoa', 'poa_abpoa__all'):
    for _w in (4, 8):
        _reg('lw%d_%s' % (_w, _base), '%s, --linkage-weight %d' % (_base, _w), kind='link',
             src=VARIANTS[_base]['src'], full=None, extra=['--linkage-weight', str(_w)])

# Re-baseline on a newer vg (set VNTR_STAGE3_VG to the pinned binary when calling these): the same
# graphs as mc, lw2_mc, poa_abpoa and poa_abpoa__all, called again.
_reg('rb_mc_link', 'mc as a candidate (new node IDs), re-baseline binary', kind='link', src=MC_LINK, full=None)
for _base in ('poa_abpoa', 'poa_abpoa__all'):
    _reg('rb_' + _base, '%s, re-baseline binary' % _base, kind='link', src=VARIANTS[_base]['src'], full=None)


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
    if spec['kind'] in ('abpoa', 'star', 'profile', 'mst', 'bbt') or (spec['kind'] == 'gapnorm' and spec['where'] == 'full'):
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
        if spec['kind'] == 'abpoa' and info['predicted_mb'] > mem_mb:
            info.update(status='memout', message='predicted %.0f MB' % info['predicted_mb'])
            return info
        if spec['kind'] == 'star':
            sp = {'centre': spec['centre']}
            if spec['centre'] == 'medoid':
                sp['weight'] = [int(r['weight']) for r in rows]
            elif spec['centre'] == 'majority':
                sp.update(source=full_msa_of(spec['source'], rid), source_map=mp, maj=spec.get('maj', 0.5))
            m = realign._as_method(v, {'align': star_align, 'tool': 'abpoa', 'description': spec['desc'],
                                       'params': sp})
        elif spec['kind'] == 'mst':
            m = realign._as_method(v, {'align': mst_align, 'tool': 'abpoa', 'description': spec['desc'],
                                       'params': {'nn': spec.get('nn', 1)}})
        elif spec['kind'] == 'bbt':
            grch = [i for i, r in enumerate(rows) if int(r['weight']) > 0 and
                    any(x.startswith('GRCh38#') for x in r['members'])]
            m = realign._as_method(v, {'align': bbt_align, 'tool': 'abpoa', 'description': spec['desc'],
                                       'params': {'K': spec['K'], 'attach': spec.get('attach', 'cluster'),
                                                  'grch38': grch[0] if grch else None,
                                                  'panel': [i for i, r in enumerate(rows) if int(r['weight']) > 0]}})
        elif spec['kind'] == 'profile':
            m = realign._as_method(v, {'align': profile_align, 'tool': spec['tool'], 'description': spec['desc'],
                                       'params': {'tool': spec['tool'], 'args': spec.get('args', [])}})
        else:
            m = realign._as_method(v, {'align': variant_align, 'tool': 'abpoa', 'description': spec['desc'],
                                       'params': {'flags': spec.get('flags', []), 'order': spec.get('order', 'longest'),
                                                  'weight': weight, 'seed_k': spec.get('seed_k')}})
        if spec.get('subset') == 'hap32':
            # A sample-dependent control: align only the distinct sequences hap32 uses.
            names = set(n for n, _ in msa_graph.read_fasta(os.path.join(rd, 'hap32.fa')))
            keep = set(r['id'] for r in rows if names & set(r['members']))
            sub = os.path.join(wd, 'hap32.union.fa')
            msa_graph.write_msa([(n, q) for n, q in msa_graph.read_fasta(fa) if n in keep], sub)
            fa = sub
        full = os.path.join(wd, 'full.msa.fa')
        al = realign.align_fasta(m, fa, full, threads=spec.get('threads', 1), timeout=timeout, mem_mb=mem_mb,
                                 workdir=wd, dedup=True)
        a = al.get('aligner') or {}
        info['align'] = {'status': al['status'], 'message': al.get('message'), 'seconds': a.get('seconds'),
                         'peak_rss_mb': a.get('peak_rss_mb'), 'command': a.get('command'),
                         'stage1': a.get('stage1'), 'n_seed': a.get('n_seed'),
                         'tool_version': a.get('tool_version'), 'star': a.get('star'), 'mst': a.get('mst'),
                         'bbt': a.get('bbt')}
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
        elif spec['kind'] in ('abpoa', 'star', 'profile', 'mst', 'bbt'):
            info = build_abpoa(v, spec, rid, timeout, mem_mb)
        elif spec['kind'] == 'merge':
            info = build_merge(v, spec, rid)
        elif spec['kind'] == 'coarsen':
            info = build_coarsen(v, spec, rid)
        elif spec['kind'] == 'gapnorm':
            info = build_gapnorm(v, spec, rid)
        else:
            raise ValueError('unknown kind %s' % spec['kind'])
    except Exception as e:  # noqa: BLE001  (one region must not end a batch)
        info = {'variant': v, 'region_id': rid, 'status': 'error', 'message': repr(e)[:500]}
    info['seconds'] = round(time.time() - t0, 1)
    with open(js, 'w') as f:
        json.dump(info, f, indent=1, default=str)
    return info


def proj_msa_of(v, rid):
    """The hap32 MSA a variant's graph was induced from."""
    spec = variant(v)
    if spec['kind'] == 'link':
        return os.path.join(spec['src'], rid + '.msa.fa')
    return os.path.join(CAND, v, rid + '.msa.fa')


def build_coarsen(v, spec, rid):
    recs = msa_graph.read_fasta(proj_msa_of(spec['base'], rid))
    recs = list(recs.items()) if isinstance(recs, dict) else list(recs)
    names = [n for n, _ in recs]
    rows = [q.upper() for _, q in recs]
    ncol = len(rows[0])
    cons = [len({r[c] for r in rows}) == 1 and rows[0][c] not in '-.' for c in range(ncol)]
    # anchor runs of >= k conserved columns (always keep the first and last run: the region anchors)
    runs, c = [], 0
    while c < ncol:
        if cons[c]:
            d = c
            while d < ncol and cons[d]:
                d += 1
            runs.append((c, d))
            c = d
        else:
            c += 1
    keep = [r for i, r in enumerate(runs) if r[1] - r[0] >= spec['k'] or i in (0, len(runs) - 1)]
    seqs, paths, nid = {}, collections.OrderedDict((n, []) for n in names), [0]
    def node(sq):
        nid[0] += 1
        seqs[str(nid[0])] = sq
        return str(nid[0])
    prev = 0
    for (a, b) in keep + [(ncol, ncol)]:
        alle = {}
        for n, r in zip(names, rows):
            seg = r[prev:a].replace('-', '').replace('.', '')
            if seg:
                if seg not in alle:
                    alle[seg] = node(seg)
                paths[n].append((alle[seg], '+'))
        if a < b:
            an = node(rows[0][a:b])
            for n in names:
                paths[n].append((an, '+'))
        prev = b
    out = os.path.join(CAND, v, rid + '.gfa')
    write_gfa(out + '.tmp', seqs, paths)
    poa_panel.check_paths(out + '.tmp', msa_graph.read_fasta(os.path.join(REGIONS, rid, 'hap32.fa')))
    os.replace(out + '.tmp', out)
    return {'variant': v, 'region_id': rid, 'status': 'ok', 'anchors': len(keep), 'nodes': len(seqs)}


def build_gapnorm(v, spec, rid):
    """Normalise the gap placement of the base's full-panel MSA (where='full') or of its hap32
    projection (where='proj') with gap_norm.py, then project and build the graph as finish_graph does
    for every abPOA variant (column graph, path check)."""
    import gap_norm
    rd = os.path.join(REGIONS, rid)
    info = {'variant': v, 'region_id': rid, 'base': spec['base'], 'where': spec['where'],
            'direction': spec['direction'], 'merge': spec.get('merge'), 'joint': bool(spec.get('joint'))}
    if spec['where'] == 'full':
        src = full_msa_of(spec['base'], rid)
    else:
        src = proj_msa_of(spec['base'], rid)
    if not src or not os.path.exists(src):
        info.update(status='no_base', message=str(src))
        return info
    named = [(n, r) for n, r in msa_graph.read_msa(src) if n not in msa_graph.CONSENSUS_NAMES]
    rows = [r for _, r in named]
    t0 = time.time()
    if spec.get('merge') is not None:
        rows, info['merge_stats'] = gap_norm.merge_runs(rows, spec['merge'])
    norm = gap_norm.normalise_joint if spec.get('joint') else gap_norm.normalise
    rows, info['norm'] = norm(rows, spec['direction'])
    info['norm']['seconds'] = round(time.time() - t0, 2)
    wd = tempfile.mkdtemp(prefix='it.%s.%s.' % (v, rid), dir=TMP)
    try:
        out = os.path.join(wd, 'norm.msa.fa')
        msa_graph.write_msa(list(zip([n for n, _ in named], rows)), out)
        if spec['where'] == 'full':
            os.makedirs(os.path.join(FULLMSA, v), exist_ok=True)
            with open(out, 'rb') as fi, gzip.open(os.path.join(FULLMSA, v, rid + '.msa.fa.gz'), 'wb') as fo:
                shutil.copyfileobj(fi, fo)
            proj = os.path.join(wd, 'proj.msa.fa')
            _, mp = poa_panel.union_files(rid)
            info['projection'] = panel.project(out, mp, os.path.join(rd, 'hap32.fa'), proj)
        else:
            proj = out
        info['graph'] = finish_graph(rid, proj, os.path.join(CAND, v, rid + '.gfa'), wd)
        shutil.copy(proj, os.path.join(CAND, v, rid + '.msa.fa'))
        info['status'] = 'ok'
        return info
    finally:
        shutil.rmtree(wd, ignore_errors=True)


# ------------------------------------------------------------------ build: centre-star and profile MSAs (4n)

class Capped(Exception):
    """A per-region time or memory cap was hit (status 'timeout' or 'memout')."""

    def __init__(self, status, msg):
        Exception.__init__(self, msg)
        self.status = status


_VERSIONS = {}


def tool_version(tool):
    """'<tool> <first version number it prints>', cached per process."""
    if tool not in _VERSIONS:
        if tool == 'mafft':
            _VERSIONS[tool] = 'mafft ' + realign.mafft_version()
        else:
            import re
            cmd = {'famsa': ['famsa'], 'kalign': ['kalign', '--version'], 'muscle': ['muscle', '-version'],
                   'abpoa': [ABPOA, '-v']}[tool]
            try:
                p = subprocess.run(cmd, capture_output=True, text=True, timeout=60, env=config.tool_env())
                m = re.search(r'\d+\.\d+(\.\d+)?', p.stdout + p.stderr)
                _VERSIONS[tool] = '%s %s' % (tool, m.group(0) if m else 'unknown')
            except (OSError, subprocess.SubprocessError):
                _VERSIONS[tool] = tool + ' unknown'
    return _VERSIONS[tool]


def _abpoa_rows(recs, wd, tag, flags, deadline, mem_mb, consensus=False):
    """abpoa -m 0 -r 1 FLAGS on recs in the given order -> {name: row} (with consensus: -r 0, the
    consensus sequence). Raises Capped when the region's deadline or the memory cap is reached."""
    left = deadline - time.time()
    if left <= 0:
        raise Capped('timeout', 'time cap reached before %s' % tag)
    fa, out = os.path.join(wd, tag + '.fa'), os.path.join(wd, tag + '.out.fa')
    realign._write_fa(fa, recs)
    r = realign.run_proc([ABPOA, '-m', '0', '-r', '0' if consensus else '1'] + list(flags) + [fa], out, None,
                         timeout=left, mem_mb=mem_mb, env=config.tool_env())
    if r['status'] in ('timeout', 'memout'):
        raise Capped(r['status'], 'abpoa %s: %s after %.0f s' % (tag, r['status'], r['seconds']))
    if r['status'] != 'ok':
        raise RuntimeError('abpoa %s: exit %s' % (tag, r.get('returncode')))
    return collections.OrderedDict(msa_graph.read_msa(out)), r.get('peak_rss_mb') or 0.0


def left_normalise_pair(a, b):
    """Left-normalise the indels of a two-row alignment without changing its cost. A run of gaps in one
    row over columns [s, e) moves to [s-1, e-1) when column s-1 aligns two bases and the other row
    has the same base at s-1 and e-1, so the base that moves keeps its match or mismatch. Repeats to a
    fixed point (runs that meet merge, which can only lower an affine cost). Returns the two rows."""
    a, b = list(a), list(b)
    n = len(a)
    moved = True
    while moved:
        moved = False
        for g, o in ((a, b), (b, a)):
            c = 0
            while c < n:
                if g[c] != '-' or o[c] == '-':
                    c += 1
                    continue
                e = c
                while e < n and g[e] == '-' and o[e] != '-':
                    e += 1
                end, s = e, c
                while s > 0 and g[s - 1] != '-' and o[s - 1] != '-' and o[s - 1] == o[e - 1]:
                    g[e - 1], g[s - 1] = g[s - 1], '-'
                    s, e = s - 1, e - 1
                    moved = True
                c = end
    return ''.join(a), ''.join(b)


def star_merge(centre, pairs, wd, deadline, mem_mb, stats):
    """Merge pairwise alignments to the centre into one MSA. pairs: {name: (centre row, sequence row)}.
    The columns are the centre's bases, and between each two of them (and at both ends) an insertion
    slot; the distinct insertions of one slot are aligned to each other with abPOA (longest first).
    Returns {name: MSA row}."""
    n = len(centre)
    cols, ins = {}, {}
    for k, (ca, sa) in pairs.items():
        col, slot, p = ['-'] * n, collections.defaultdict(list), 0
        for x, y in zip(ca, sa):
            if x == '-':
                if y != '-':
                    slot[p].append(y)
            else:
                col[p] = y
                p += 1
        if p != n:
            raise ValueError('pairwise centre row of %s spells %d of %d bases' % (k, p, n))
        cols[k], ins[k] = col, {q: ''.join(v) for q, v in slot.items()}
    byslot = collections.defaultdict(set)
    for d in ins.values():
        for q, s in d.items():
            byslot[q].add(s)
    slot_rows = {}
    for q, strs in byslot.items():
        if len(strs) == 1:
            s = next(iter(strs))
            slot_rows[q] = (len(s), {s: s})
            continue
        recs = sorted(strs, key=lambda x: (-len(x), x))
        got, pk = _abpoa_rows([('i%d' % i, s) for i, s in enumerate(recs)], wd, 'slot', [], deadline, mem_mb)
        al = {s: got['i%d' % i] for i, s in enumerate(recs)}
        slot_rows[q] = (len(al[recs[0]]), al)
        stats['slots_aligned'] += 1
        stats['slot_seqs_aligned'] += len(recs)
        stats['peak_rss_mb'] = max(stats['peak_rss_mb'], pk)
    stats['slots'] = len(slot_rows)
    stats['insert_columns'] = sum(w for w, _ in slot_rows.values())
    out = {}
    for k in pairs:
        parts = []
        for q in range(n + 1):
            if q in slot_rows:
                w, al = slot_rows[q]
                s = ins[k].get(q)
                parts.append(al[s] if s else '-' * w)
            if q < n:
                parts.append(cols[k][q])
        out[k] = ''.join(parts)
    return out


def star_medoid(recs, weight=None, k=None):
    """The sequence of recs (named s<i>) minimising the sum over all others of w_j (1 - multiset k-mer
    Jaccard), the mst distance; w_j is weight[j] (panel haplotypes carrying s<j>), or 1. Ties go to the
    lower s-number."""
    toks = [kmer_tokens(s, k or MST_K) for _, s in recs]
    w = [(weight[int(n[1:])] if weight else 1) for n, _ in recs]
    best = None
    for i in range(len(recs)):
        d = 0.0
        for j in range(len(recs)):
            if i != j:
                inter = len(toks[i] & toks[j])
                uni = len(toks[i]) + len(toks[j]) - inter
                d += w[j] * (1.0 - (inter / uni if uni else 1.0))
        if best is None or d < best[0]:
            best = (d, i)
    return recs[best[1]][1]


def majority_consensus(msa_gz, map_tsv, maj=0.5):
    """Majority consensus of a full-panel MSA (rows named by distinct id): a column is kept when the
    alleles with a base there carry at least MAJ of the panel weight (map 'weight', the panel haplotypes
    with that allele), and spells its heaviest base. So an indel carried by less than MAJ of the panel
    is left out of the centre."""
    w = {r['id']: max(int(r['weight']), 1) for r in panel.read_map(map_tsv)}
    rows = [(n, r) for n, r in msa_graph.read_msa(msa_gz) if n not in msa_graph.CONSENSUS_NAMES]
    tot = float(sum(w.get(n, 1) for n, _ in rows))
    out = []
    for c in range(len(rows[0][1])):
        cnt = collections.Counter()
        for n, r in rows:
            if r[c] in 'ACGT':
                cnt[r[c]] += w.get(n, 1)
        if cnt and sum(cnt.values()) >= maj * tot:
            out.append(max(sorted(cnt), key=lambda b: cnt[b]))
    return ''.join(out)


def star_align(in_fa, out_fa, threads=1, workdir=None, timeout=TIMEOUT, mem_mb=MEM_MB, centre='consensus',
               weight=None, source=None, source_map=None, maj=0.5, **_):
    """realign.py plugin: centre-star MSA of the distinct masked sequences (named s<k>). The centre is
    abPOA's consensus of all of them (abpoa -r 0, longest first, as poa_abpoa__all is built), with
    centre='ref' s0, which is CHM13's sequence (panel.union lists it first), with centre='longest'
    the longest sequence, with centre='medoid' the sequence of least weighted k-mer distance to all the
    others (star_medoid), or with centre='majority' the majority consensus of the full-panel MSA at
    SOURCE (majority_consensus). Each sequence is aligned
    to the centre with abpoa -m 0 -b -1 (global, unbanded, default scores), its indels left-normalised
    (left_normalise_pair), and the alignments merged (star_merge). The whole region shares one time cap;
    each abPOA run has the memory cap."""
    t0 = time.time()
    deadline = t0 + timeout
    recs = msa_graph.read_fasta(in_fa)
    wd = tempfile.mkdtemp(prefix='star.', dir=workdir)
    stats = collections.Counter(peak_rss_mb=0.0)
    res = {'command': ['abpoa', '-m', '0', '-b', '-1', '(pairwise to the centre)'], 'tool_version': tool_version('abpoa')}
    cseq = ''
    try:
        if centre == 'ref':
            if recs[0][0] != 's0':
                raise ValueError('centre=ref needs s0 (CHM13) first, got %s' % recs[0][0])
            cseq = recs[0][1]
        elif centre == 'longest':
            cseq = sorted(recs, key=lambda r: (-len(r[1]), int(r[0][1:])))[0][1]
        elif centre == 'medoid':
            cseq = star_medoid(recs, weight)
        elif centre == 'majority':
            cseq = majority_consensus(source, source_map, maj)
        else:
            byl = sorted(recs, key=lambda r: (-len(r[1]), int(r[0][1:])))
            got, pk = _abpoa_rows(byl, wd, 'cons', [], deadline, mem_mb, consensus=True)
            cseq = next(iter(got.values()))
            stats['peak_rss_mb'] = max(stats['peak_rss_mb'], pk)
        if not cseq or set(cseq) - set('ACGT'):
            raise ValueError('centre is empty or not ACGT')
        pairs = collections.OrderedDict()
        for name, s in recs:
            if s == cseq:
                pairs[name] = (cseq, s)
                stats['identical_to_centre'] += 1
                continue
            got, pk = _abpoa_rows([('c', cseq), ('q', s)], wd, 'pair', ['-b', '-1'], deadline, mem_mb)
            stats['peak_rss_mb'] = max(stats['peak_rss_mb'], pk)
            ca, sa = got['c'], got['q']
            keep = [i for i in range(len(ca)) if ca[i] != '-' or sa[i] != '-']
            ca, sa = ''.join(ca[i] for i in keep), ''.join(sa[i] for i in keep)
            if ca.replace('-', '') != cseq or sa.replace('-', '') != s:
                raise ValueError('abpoa pairwise rows do not spell their inputs (%s)' % name)
            ca, sa = left_normalise_pair(ca, sa)
            pairs[name] = (ca, sa)
        stats['pairs'] = len(pairs)
        rows = star_merge(cseq, pairs, wd, deadline, mem_mb, stats)
        msa_graph.write_msa([(n, rows[n]) for n, _ in recs], out_fa)
        res['status'] = 'ok'
    except Capped as e:
        res.update(status=e.status, message=str(e))
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    res['seconds'] = round(time.time() - t0, 2)
    res['peak_rss_mb'] = stats.pop('peak_rss_mb')
    res['star'] = dict(stats, centre=centre, centre_len=len(cseq))
    return res


def profile_align(in_fa, out_fa, threads=1, workdir=None, timeout=TIMEOUT, mem_mb=MEM_MB, tool='kalign', args=(),
                  **_):
    """realign.py plugin: a progressive profile aligner on the distinct masked sequences, under the
    region's time and memory cap (realign.run_proc). mafft goes through realign.mafft_align. Rows are
    rewritten upper-case with '-' gaps."""
    if tool == 'mafft':
        r = realign.mafft_align(in_fa, out_fa, threads=threads, workdir=workdir, timeout=timeout, mem_mb=mem_mb,
                                args=args)
        r['tool_version'] = tool_version('mafft')
        return r
    exe = shutil.which(tool, path=config.tool_env()['PATH']) or tool
    part = os.path.join(workdir, 'profile.out.fa')
    if tool == 'famsa':
        cmd = [exe, '-t', str(threads)] + list(args) + [in_fa, part]
    elif tool == 'kalign':
        cmd = [exe, '-i', in_fa, '-o', part, '--type', 'dna', '-n', str(threads)] + list(args)
    elif tool == 'muscle':
        cmd = [exe, '-super5', in_fa, '-output', part, '-threads', str(threads)] + list(args)
    else:
        raise ValueError('unknown profile tool %s' % tool)
    lg = os.path.join(workdir, tool + '.log')
    r = realign.run_proc(cmd, lg + '.out', lg, timeout=timeout, mem_mb=mem_mb, env=config.tool_env(), cwd=workdir)
    r['command'] = [tool] + [x if x not in (in_fa, part) else ('<in.fa>' if x == in_fa else '<out.fa>')
                             for x in cmd[1:]]
    r['tool_version'] = tool_version(tool)
    if r['status'] == 'ok' and not (os.path.exists(part) and os.path.getsize(part)):
        r['status'] = 'error'
    if r['status'] == 'ok':
        msa_graph.write_msa(msa_graph.read_msa(part), out_fa)
    else:
        r['message'] = realign._tail(lg)
    return r


# ------------------------------------------------------------------ build: nearest-neighbour threading (4q)

MST_K = 15
ABPOA_COST = {'mismatch': 4, 'open1': 4, 'ext1': 2, 'open2': 24, 'ext2': 1}   # abPOA's default scores
PLACE = {'match': 2.0, 'mismatch': -3.0, 'new': -1.0, 'skip': -2.5}           # place_insert scores
PLACE_MAX_CELLS = 4000000
PAIR_MAX_CELLS = 600000000     # a 22 kb x 22 kb unbanded pair fits the 12 GB cap; 30 kb x 30 kb does not
_BIT = {'A': 1, 'C': 2, 'G': 4, 'T': 8}


def kmer_tokens(seq, k=MST_K):
    """The k-mer multiset of seq as a set: the j-th occurrence of a k-mer is its own token, so the Jaccard
    index of two token sets is the multiset Jaccard (sum of min counts over sum of max counts), which sees
    copy-number differences that a k-mer set would not."""
    seen, out = collections.Counter(), set()
    for i in range(len(seq) - k + 1):
        w = seq[i:i + k]
        seen[w] += 1
        out.add(hash((w, seen[w])))
    return out


def mst_order(recs, k=MST_K):
    """Prim's minimum spanning tree over the distinct sequences, distance 1 - multiset k-mer Jaccard,
    grown from recs[0] (CHM13, s0). Returns (order, parent, dist): order is the order the sequences join
    the tree, each joining at its nearest tree member parent[i]; dist(i, j) is the distance function."""
    toks = [kmer_tokens(s, k) for _, s in recs]
    n = len(recs)
    cache = {}

    def dist(i, j):
        if i == j:
            return 0.0
        key = (i, j) if i < j else (j, i)
        if key not in cache:
            a, b = toks[i], toks[j]
            inter = len(a & b)
            uni = len(a) + len(b) - inter
            cache[key] = 1.0 - (inter / uni if uni else (1.0 if recs[i][1] == recs[j][1] else 0.0))
        return cache[key]

    inside = [False] * n
    best = [float('inf')] * n
    parent = [None] * n
    order = []
    cur = 0
    inside[0] = True
    order.append(0)
    for _ in range(n - 1):
        for j in range(n):
            if not inside[j]:
                d = dist(cur, j)
                if d < best[j]:
                    best[j], parent[j] = d, cur
        cur = min((j for j in range(n) if not inside[j]), key=lambda j: (best[j], j))
        inside[cur] = True
        order.append(cur)
    return order, parent, dist


def pair_cost(a, b):
    """abPOA's default cost of a two-row alignment (no match bonus): 4 per mismatch and min(4 + 2L, 24 + L)
    per gap run of length L in either row."""
    c, n = 0, len(a)
    i = 0
    while i < n:
        x, y = a[i], b[i]
        if x != '-' and y != '-':
            c += ABPOA_COST['mismatch'] if x != y else 0
            i += 1
            continue
        g = 0 if x == '-' else 1
        e = i
        while e < n and (a[e] if g == 0 else b[e]) == '-' and (b[e] if g == 0 else a[e]) != '-':
            e += 1
        L = e - i
        c += min(ABPOA_COST['open1'] + ABPOA_COST['ext1'] * L, ABPOA_COST['open2'] + ABPOA_COST['ext2'] * L)
        i = e
    return c


def pair_align(pseq, qseq, wd, deadline, mem_mb, stats):
    """Global pairwise alignment of qseq to pseq (abpoa -m 0 -b -1, default scores), both-gap columns
    dropped and indels left-normalised. A pair over PAIR_MAX_CELLS cells, or one whose unbanded run hits the
    memory cap, is aligned with abPOA's default adaptive band instead. Returns (parent row, query row,
    cost)."""
    try:
        if len(pseq) * len(qseq) > PAIR_MAX_CELLS:
            raise Capped('memout', 'predicted')
        got, pk = _abpoa_rows([('p', pseq), ('q', qseq)], wd, 'pair', ['-b', '-1'], deadline, mem_mb)
    except Capped as e:
        if e.status != 'memout':
            raise
        got, pk = _abpoa_rows([('p', pseq), ('q', qseq)], wd, 'pair', [], deadline, mem_mb)
        stats['banded'] += 1
    stats['peak_rss_mb'] = max(stats['peak_rss_mb'], pk)
    stats['pairwise'] += 1
    pa, qa = got['p'], got['q']
    keep = [i for i in range(len(pa)) if pa[i] != '-' or qa[i] != '-']
    pa, qa = ''.join(pa[i] for i in keep), ''.join(qa[i] for i in keep)
    if pa.replace('-', '') != pseq or qa.replace('-', '') != qseq:
        raise ValueError('abpoa pairwise rows do not spell their inputs')
    pa, qa = left_normalise_pair(pa, qa)
    return pa, qa, pair_cost(pa, qa)


def place_insert(ins, masks, sc=PLACE):
    """Place an inserted string into a run of existing columns (masks: the bases each holds, as _BIT
    masks; the parent has a gap in all of them) plus new columns. Global DP: a base into an existing
    column scores match or mismatch, a base into a new column 'new', and skipping an existing column is
    free before the first and after the last placed base and costs 'skip' between them; ties put the
    insertion as far left as possible. Returns the slot as a list of (existing column index or None for a
    new column, query base index or None)."""
    m, g = len(ins), len(masks)
    H = [[0.0] * (g + 1) for _ in range(m + 1)]
    P = [bytearray(g + 1) for _ in range(m + 1)]    # 0 diagonal, 1 new column, 2 skip
    for j in range(1, g + 1):
        P[0][j] = 2
    M, X, N, S = sc['match'], sc['mismatch'], sc['new'], sc['skip']
    for i in range(1, m + 1):
        bit = _BIT.get(ins[i - 1], 0)
        prev, cur, pr = H[i - 1], H[i], P[i]
        cur[0] = prev[0] + N
        pr[0] = 1
        last = i == m
        sk = 0.0 if last else S
        for j in range(1, g + 1):
            d = prev[j - 1] + (M if masks[j - 1] & bit else X)
            u = prev[j] + N
            lft = cur[j - 1] + sk
            if last and lft >= d and lft >= u:
                cur[j], pr[j] = lft, 2
            elif d >= u and d >= lft:
                cur[j], pr[j] = d, 0
            elif u >= lft:
                cur[j], pr[j] = u, 1
            else:
                cur[j], pr[j] = lft, 2
    out, i, j = [], m, g
    while i > 0 or j > 0:
        p = P[i][j]
        if p == 0:
            out.append((j - 1, i - 1))
            i, j = i - 1, j - 1
        elif p == 1:
            out.append((None, i - 1))
            i -= 1
        else:
            out.append((j - 1, None))
            j -= 1
    out.reverse()
    return out


class ThreadedMSA(object):
    """An MSA grown one row at a time, each new row placed against one parent row already in it. Columns
    are ids in self.order; rowcols[name] lists the column of each of the row's bases; mask[c] holds the
    bases column c has received."""

    def __init__(self, name, seq):
        self.order = list(range(len(seq)))
        self.mask = [_BIT.get(b, 0) for b in seq]
        self.rowcols = {name: list(range(len(seq)))}

    @classmethod
    def from_msa(cls, rows):
        """Start from an existing MSA (rows: {name: aligned row}) instead of one row; all-gap columns are
        dropped."""
        rows = collections.OrderedDict(rows)
        width = len(next(iter(rows.values())))
        keep = [c for c in range(width) if any(r[c] != '-' for r in rows.values())]
        self = cls.__new__(cls)
        self.order = list(range(len(keep)))
        self.mask = [0] * len(keep)
        self.rowcols = {}
        for name, r in rows.items():
            cols = []
            for j, c in enumerate(keep):
                b = r[c]
                if b != '-':
                    cols.append(j)
                    self.mask[j] |= _BIT.get(b, 0)
            self.rowcols[name] = cols
        return self

    def _new(self, base):
        self.mask.append(_BIT.get(base, 0))
        return len(self.mask) - 1

    def add(self, name, qseq, parent, pa, qa, stats):
        """Add qseq, aligned to row PARENT as (pa, qa): a base aligned to a parent base takes that base's
        column; a deletion is gaps in the parent's columns; the bases inserted between two parent bases go
        into that slot's existing columns (the parent has gaps there) and new ones, by place_insert."""
        pc = self.rowcols[parent]
        n = len(pc)
        pos = {c: i for i, c in enumerate(self.order)}
        qcols = []
        slots = {}                     # parent slot p (before parent base p) -> (first qcols index, string)
        p, ins = 0, []
        for x, y in zip(pa, qa):
            if x == '-':
                ins.append(y)
                continue
            if ins:
                slots[p] = (len(qcols), ''.join(ins))
                qcols.extend([None] * len(ins))
                ins = []
            if y != '-':
                qcols.append(pc[p])
            p += 1
        if ins:
            slots[p] = (len(qcols), ''.join(ins))
            qcols.extend([None] * len(ins))
        if p != n:
            raise ValueError('parent row spells %d of %d bases' % (p, n))
        segs = []                      # (start, end) of self.order replaced by a new column list
        for sp in sorted(slots):
            q0, s = slots[sp]
            lo = pos[pc[sp - 1]] + 1 if sp > 0 else 0
            hi = pos[pc[sp]] if sp < n else len(self.order)
            old = self.order[lo:hi]
            stats['slots'] += 1
            if not old or len(s) * len(old) > PLACE_MAX_CELLS:
                if old:
                    stats['place_capped'] += 1
                ops = [(None, t) for t in range(len(s))] + [(j, None) for j in range(len(old))]
            else:
                ops = place_insert(s, [self.mask[c] for c in old])
            seg = []
            for j, t in ops:
                if j is None:
                    c = self._new(s[t])
                    stats['new_columns'] += 1
                else:
                    c = old[j]
                    if t is not None:
                        stats['reused_columns'] += 1
                seg.append(c)
                if t is not None:
                    qcols[q0 + t] = c
            segs.append((lo, hi, seg))
        if segs:
            out, last = [], 0
            for lo, hi, seg in segs:
                out.extend(self.order[last:lo])
                out.extend(seg)
                last = hi
            out.extend(self.order[last:])
            self.order = out
        for c, b in zip(qcols, qseq):
            self.mask[c] |= _BIT.get(b, 0)
        self.rowcols[name] = qcols

    def write(self, recs, out_fa):
        """Write the rows of recs (name, sequence) as an MSA; returns the column count."""
        idx = {c: i for i, c in enumerate(self.order)}
        ncol = len(self.order)
        out = []
        for name, seq in recs:
            row = ['-'] * ncol
            last = -1
            for c, b in zip(self.rowcols[name], seq):
                k = idx[c]
                if k <= last:
                    raise ValueError('row %s: columns out of order' % name)
                row[k] = b
                last = k
            out.append((name, ''.join(row)))
        msa_graph.write_msa(out, out_fa)
        return ncol


def mst_align(in_fa, out_fa, threads=1, workdir=None, timeout=TIMEOUT, mem_mb=MEM_MB, nn=1, k=MST_K, **_):
    """realign.py plugin: nearest-neighbour threading of the distinct masked sequences (named s<k>, s0 =
    CHM13). The sequences join a minimum spanning tree (mst_order) in Prim's order from CHM13; each is
    aligned pairwise (pair_align) to its tree parent, or with nn > 1 to each of its nn nearest
    already-aligned sequences by k-mer distance keeping the lowest-cost alignment, and added to the MSA on
    that parent's columns (ThreadedMSA.add). One time cap for the region; each abPOA run has the memory
    cap."""
    t0 = time.time()
    deadline = t0 + timeout
    recs = msa_graph.read_fasta(in_fa)
    wd = tempfile.mkdtemp(prefix='mst.', dir=workdir)
    stats = collections.Counter(peak_rss_mb=0.0)
    res = {'command': ['abpoa', '-m', '0', '-b', '-1', '(pairwise to the nearest aligned allele)'],
           'tool_version': tool_version('abpoa')}
    try:
        if recs[0][0] != 's0':
            raise ValueError('mst needs s0 (CHM13) first, got %s' % recs[0][0])
        order, parent, dist = mst_order(recs, k)
        stats['tree_seconds'] = round(time.time() - t0, 2)
        msa = ThreadedMSA(recs[0][0], recs[0][1])
        done = [order[0]]
        for x in order[1:]:
            if nn > 1:
                cand = sorted(done, key=lambda j: (dist(x, j), j))[:nn]
            else:
                cand = [parent[x]]
            best = None
            for j in cand:
                pa, qa, c = pair_align(recs[j][1], recs[x][1], wd, deadline, mem_mb, stats)
                if best is None or c < best[3]:
                    best = (j, pa, qa, c)
            if best[0] != parent[x]:
                stats['parent_not_tree'] += 1
            stats['cost'] += best[3]
            msa.add(recs[x][0], recs[x][1], recs[best[0]][0], best[1], best[2], stats)
            done.append(x)
        stats['columns'] = msa.write(recs, out_fa)
        res['status'] = 'ok'
    except Capped as e:
        res.update(status=e.status, message=str(e))
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    res['seconds'] = round(time.time() - t0, 2)
    res['peak_rss_mb'] = stats.pop('peak_rss_mb')
    res['mst'] = dict(stats, nn=nn, k=k)
    return res


# ------------------------------------------------------------------ build: backbone plus threading (4r)

BBT_MAX_SWEEPS = 50


def kmedoids(n, dist, K, fixed=(), max_sweeps=BBT_MAX_SWEEPS):
    """K medoids of points 0..n-1 under dist. The fixed medoids come first and never move; farthest-first
    traversal fills the rest up to K, then FasterPAM's eager swaps (Schubert and Rousseeuw 2021) move the
    free medoids while a swap lowers the summed distance of the points to their nearest medoid, until a
    pass over every non-medoid finds none (at most max_sweeps passes). Ties go to the lower index. Returns
    (medoids, cluster, swaps): cluster[p] is the index into medoids of p's nearest medoid (ties to the
    earlier medoid)."""
    D = [[dist(i, j) for j in range(n)] for i in range(n)]
    med = list(dict.fromkeys(fixed)) or [0]
    nfix = len(list(dict.fromkeys(fixed)))
    if K >= n:
        med += [i for i in range(n) if i not in med]
    near = [min(D[m][p] for m in med) for p in range(n)]
    while len(med) < min(K, n):
        ism = set(med)
        c = max((p for p in range(n) if p not in ism), key=lambda p: (near[p], -p))
        med.append(c)
        near = [min(a, b) for a, b in zip(near, D[c])]
    inf = float('inf')

    def nearest_two():
        nn, dn, ds = [0] * n, [0.0] * n, [0.0] * n
        for p in range(n):
            b1, b2, i1 = inf, inf, 0
            for i, m in enumerate(med):
                d = D[m][p]
                if d < b1:
                    b1, b2, i1 = d, b1, i
                elif d < b2:
                    b2 = d
            nn[p], dn[p], ds[p] = i1, b1, b2
        rem = [0.0] * len(med)
        for p in range(n):
            rem[nn[p]] += ds[p] - dn[p]
        return nn, dn, ds, rem

    nn, dn, ds, rem = nearest_two()
    swaps = 0
    free = list(range(nfix, len(med)))
    if free and K < n:
        since, xc, passes = 0, 0, 0
        while since < n and passes < max_sweeps * n:
            passes += 1
            if xc not in set(med):
                delta = list(rem)
                acc = 0.0
                Dc = D[xc]
                for p in range(n):
                    d, a = Dc[p], dn[p]
                    if d < a:
                        acc += d - a
                        delta[nn[p]] += a - ds[p]
                    elif d < ds[p]:
                        delta[nn[p]] += d - ds[p]
                i = min(free, key=lambda j: (delta[j], j))
                if delta[i] + acc < -1e-12:
                    med[i] = xc
                    swaps += 1
                    nn, dn, ds, rem = nearest_two()
                    since = 0
                else:
                    since += 1
            else:
                since += 1
            xc = (xc + 1) % n
    return med, nn, swaps


def _prim_from(root, members, dist):
    """Prim's tree over members grown from root: [(member, parent)] in joining order (root excluded)."""
    best = {m: (dist(root, m), root) for m in members if m != root}
    out = []
    while best:
        x = min(best, key=lambda m: (best[m][0], m))
        out.append((x, best.pop(x)[1]))
        for m in best:
            d = dist(x, m)
            if d < best[m][0]:
                best[m] = (d, x)
    return out


def bbt_plan(recs, panel, grch38, K, k=MST_K, attach='cluster'):
    """The threading plan of bbt_align, from the sequences and full-panel membership alone. The distinct
    full-panel sequences (indices panel into recs; recs[0] is CHM13's) are ordered CHM13, then GRCh38's
    when it is distinct, then by length and sequence, and split into K clusters by kmedoids on 1 - multiset
    k-mer Jaccard (the mst distance), with CHM13 and GRCh38 fixed as medoids. Returns (reps, steps, info):
    reps the medoids; steps [(x, parent)] in threading order: each panel sequence on its nearest
    already-placed member of its cluster (Prim's tree grown from the medoid; with attach='medoid', on the
    medoid itself, nearest first), then each sequence outside the panel (hap32-only recombinant paths) on
    the nearest panel member of its nearest medoid's cluster (or on that medoid), so that no panel row's
    placement depends on them."""
    toks = [kmer_tokens(s, k) for _, s in recs]
    cache = {}

    def dist(i, j):
        if i == j:
            return 0.0
        key = (i, j) if i < j else (j, i)
        if key not in cache:
            a, b = toks[i], toks[j]
            inter = len(a & b)
            uni = len(a) + len(b) - inter
            cache[key] = 1.0 - (inter / uni if uni else (1.0 if recs[i][1] == recs[j][1] else 0.0))
        return cache[key]

    have = set(range(len(recs)))
    pan = [x for x in dict.fromkeys(panel) if x in have]
    fixed = [0] + ([grch38] if grch38 is not None and grch38 in pan and grch38 != 0 else [])
    pan = fixed + sorted((x for x in pan if x not in fixed), key=lambda x: (len(recs[x][1]), recs[x][1]))
    med, nn, swaps = kmedoids(len(pan), lambda a, b: dist(pan[a], pan[b]), K, fixed=range(len(fixed)))
    reps = [pan[m] for m in med]
    clusters = collections.OrderedDict((r, [r]) for r in reps)
    for p, i in enumerate(nn):
        if pan[p] not in clusters:
            clusters[reps[i]].append(pan[p])
    steps = []
    for r, mem in clusters.items():
        if attach == 'medoid':
            steps.extend((x, r) for x in sorted((m for m in mem if m != r), key=lambda m: (dist(r, m), m)))
        else:
            steps.extend(_prim_from(r, mem, dist))
    inpan = set(pan)
    extra = [x for x in range(len(recs)) if x not in inpan]
    for x in extra:
        r = reps[min(range(len(reps)), key=lambda i: (dist(x, reps[i]), i))]
        steps.append((x, r if attach == 'medoid' else min(clusters[r], key=lambda m: (dist(x, m), m))))
    sizes = sorted(len(v) for v in clusters.values())
    info = {'attach': attach, 'K': len(reps), 'n_panel': len(pan), 'n_outside_panel': len(extra), 'swaps': swaps,
            'fixed': [recs[x][0] for x in fixed], 'cluster_max': sizes[-1],
            'cluster_median': statistics.median(sizes), 'singletons': sum(1 for s in sizes if s == 1)}
    return reps, steps, info


def bbt_align(in_fa, out_fa, threads=1, workdir=None, timeout=TIMEOUT, mem_mb=MEM_MB, K=32, panel=None,
              grch38=None, k=MST_K, attach='cluster', **_):
    """realign.py plugin: backbone plus threading of the distinct masked sequences (named s<i>, s0 =
    CHM13; panel and grch38 are such indices i, from the region's panel map). bbt_plan picks K k-medoid
    representatives of the full-panel sequences; abPOA aligns them (abpoa -m 0 -r 1, defaults, longest
    first, as poa_abpoa) into the backbone, and every other sequence is aligned pairwise (pair_align) to
    its plan parent and added on that parent's columns (ThreadedMSA.add: insertions reuse the slot's gap
    columns where they fit, so recurrent insertions stack). One time cap for the region; each abPOA run
    has the memory cap."""
    t0 = time.time()
    deadline = t0 + timeout
    recs = msa_graph.read_fasta(in_fa)
    idx = {int(n[1:]): i for i, (n, _) in enumerate(recs)}
    wd = tempfile.mkdtemp(prefix='bbt.', dir=workdir)
    stats = collections.Counter(peak_rss_mb=0.0)
    res = {'command': ['abpoa', '-m', '0', '-r', '1', '(backbone)', '+', 'abpoa', '-m', '0', '-b', '-1',
                       '(pairwise to the plan parent)'], 'tool_version': tool_version('abpoa')}
    plan = {}
    try:
        if recs[0][0] != 's0':
            raise ValueError('bbt needs s0 (CHM13) first, got %s' % recs[0][0])
        pan = [idx[x] for x in (panel if panel is not None else sorted(idx)) if x in idx]
        reps, steps, plan = bbt_plan(recs, pan, idx.get(grch38) if grch38 is not None else None, K, k, attach)
        stats['plan_seconds'] = round(time.time() - t0, 2)
        if len(reps) == 1:
            msa = ThreadedMSA(recs[reps[0]][0], recs[reps[0]][1])
        else:
            byl = sorted((recs[r] for r in reps), key=lambda r: (-len(r[1]), int(r[0][1:])))
            t1 = time.time()
            rows, pk = _abpoa_rows(byl, wd, 'backbone', [], deadline, mem_mb)
            stats['backbone_seconds'] = round(time.time() - t1, 2)
            stats['peak_rss_mb'] = max(stats['peak_rss_mb'], pk)
            for n, s in byl:
                if rows[n].replace('-', '') != s:
                    raise ValueError('abpoa backbone row %s does not spell its input' % n)
            msa = ThreadedMSA.from_msa(rows)
            stats['backbone_columns'] = len(msa.order)
        for x, p in steps:
            pa, qa, c = pair_align(recs[p][1], recs[x][1], wd, deadline, mem_mb, stats)
            stats['cost'] += c
            msa.add(recs[x][0], recs[x][1], recs[p][0], pa, qa, stats)
        stats['columns'] = msa.write(recs, out_fa)
        for (n, s), (n2, row) in zip(recs, msa_graph.read_msa(out_fa)):
            if n != n2 or row.replace('-', '') != s:
                raise ValueError('bbt MSA row %s does not spell its input' % n)
        res['status'] = 'ok'
    except Capped as e:
        res.update(status=e.status, message=str(e))
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    res['seconds'] = round(time.time() - t0, 2)
    res['peak_rss_mb'] = stats.pop('peak_rss_mb')
    res['bbt'] = dict(stats, **plan)
    return res


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
        if spec['kind'] in ('merge', 'coarsen', 'gapnorm'):
            cmd_build([spec['base']], ids, jobs, False, timeout, mem_mb)
        # biggest first so the long jobs start early
        size = {r: os.path.getsize(poa_panel.union_files(r)[0]) if os.path.exists(poa_panel.union_files(r)[0]) else 0
                for r in ids}
        items = [(v, r, timeout, mem_mb, force) for r in sorted(ids, key=lambda r: -size[r])]
        for (_, rid, *_), info in _pool(build_one, jobs, items):
            if info.get('status') != 'ok' or force:
                log('build', v, rid, info.get('status'), info.get('seconds'), info.get('message', '') or '')


# ------------------------------------------------------------------ stage 0

def snarl_counts(gfa):
    """{'snarls': all non-trivial snarls, 'sites': the top-level ones} of a GFA (vg snarls | vg view -R)."""
    p1 = subprocess.Popen([SNARL_VG, 'snarls', gfa], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.run([SNARL_VG, 'view', '-R', '-'], stdin=p1.stdout, capture_output=True, text=True)
    p1.stdout.close()
    if p1.wait() != 0 or p2.returncode != 0:
        return {'status': 'error'}
    sn = [json.loads(x) for x in p2.stdout.splitlines() if x.strip()]
    return {'snarls': len(sn), 'sites': sum(1 for s in sn if 'parent' not in s), 'vg': os.path.basename(SNARL_VG)}


def snarls_one(v, rid, gfa, force=False):
    """Site count of a variant's graph: `vg snarls` (non-trivial snarls) on the GFA, read with `vg view -R`;
    'snarls' counts all of them and 'sites' the top-level ones (no parent), which are vg call's sites.
    Cached in stage0/<variant>/<id>.snarls.json."""
    out = os.path.join(STAGE0, v, rid + '.snarls.json')
    if os.path.exists(out) and not force:
        return json.load(open(out))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    res = snarl_counts(gfa)
    if res.get('status') == 'error':
        return res
    with open(out, 'w') as f:
        json.dump(res, f)
    return res


def stage0_one(v, rid, force=False):
    out = os.path.join(STAGE0, v, rid + '.json')
    gfa = cand_gfa(v, rid)
    if not os.path.exists(gfa):
        return 'no_graph'
    snarls_one(v, rid, gfa, force)
    if os.path.exists(out) and not force:
        return 'cached'
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
    sp = os.path.join(STAGE0, v, rid + '.snarls.json')
    sn = json.load(open(sp)) if os.path.exists(sp) else {}
    return {'nodes': (d.get('size') or {}).get('nodes'), 'nodes_per_kb': (d.get('size') or {}).get('nodes_per_kb'),
            'kmer_frac_extra': (d.get('redundancy') or {}).get('kmer_frac_extra'),
            'cost_over_opt': (d.get('alignment') or {}).get('all_cost_over_opt'), 'sites': sn.get('sites')}


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
                            str(threads), '--compact'], stdout=f, stderr=subprocess.STDOUT,
                           env=dict(call_env(), VG_CALL_EXTRA=' '.join(variant(v).get('extra', []))), cwd=REPO)
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


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return statistics.mean(xs) if xs else None


def cmd_table(variants=None, quiet=False):
    ts = testset()
    group = {r['region_id']: r['group'] for r in ts}
    ids = [r['region_id'] for r in ts]
    variants = variants or known_variants()
    cols = ['variant', 'region_id', 'group', 'build_status', 'align_s', 'peak_rss_mb', 'nodes', 'nodes_per_kb',
            'kmer_frac_extra', 'cost_over_opt', 'sites', 'ed_reps', 'ed_median', 'vs_mc', 'vs_poa_abpoa', 'vs_poa_abpoa__all']
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
                   'sites': m.get('sites'), 'ed_reps': ','.join(str(x) for x in sc), 'ed_median': _med(sc) if sc else None}
            for bname in BASELINES:
                row['vs_' + bname] = '' if bname == v else rule(sc, base_scores[bname][rid])
            rows.append(row)
    with open(RESULTS, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join('' if r[c] is None else str(r[c]) for c in cols) + '\n')
    scols = ['variant', 'built', 'med_kmer_frac_extra', 'med_cost_over_opt', 'med_nodes_per_kb', 'max_align_s',
             'called', 'sum_ed_median', 'sum_ed_median_mc', 'sum_ed_median_poa', 'bw_vs_mc', 'bw_vs_poa_abpoa',
             'bw_vs_poa_abpoa__all', 'mean_kmer_frac_extra', 'mean_cost_over_opt', 'med_sites', 'mean_sites']
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
             'called': len(called), 'sum_ed_median': sum(r['ed_median'] for r in called) if called else None,
             'mean_kmer_frac_extra': _mean([r['kmer_frac_extra'] for r in rs]),
             'mean_cost_over_opt': _mean([r['cost_over_opt'] for r in rs]),
             'med_sites': _med([r['sites'] for r in rs]), 'mean_sites': _mean([r['sites'] for r in rs])}
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
        print('%-16s %5s %8s %8s %6s %7s %7s %6s %8s %8s %8s %7s %7s %7s' % (
            'variant', 'built', 'kmerx', 'cost/opt', 'sites', 'nod/kb', 'maxs', 'called', 'sumED', 'ED(mc)', 'ED(poa)',
            'b/w mc', 'b/w poa', 'b/w all'))
        for s in summ:
            print('%-16s %5s %8s %8s %6s %7s %7s %6s %8s %8s %8s %7s %7s %7s' % (
                s['variant'], s['built'], fmt(s['med_kmer_frac_extra']), fmt(s['med_cost_over_opt']),
                '-' if s['med_sites'] is None else '%g' % s['med_sites'],
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
