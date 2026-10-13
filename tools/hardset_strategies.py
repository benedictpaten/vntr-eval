#!/usr/bin/env python3
"""hardset_strategies.py -- realignment strategies on the hard set (tools/hardset.py), scored by
pgrealign qc's pair and reference F1 against the same pairwise optima qc uses.

    python3 tools/hardset_strategies.py optima --jobs 3            # cache each region's optimal pairs
    python3 tools/hardset_strategies.py run star famsa --jobs 3    # align with strategies (resumes)
    python3 tools/hardset_strategies.py score                      # -> results/hardset_chr20.tsv, .md

`run` runs one (strategy, region) per child process, records its wall time and the peak resident memory
of its whole process group (sampled every second), and kills it past --timeout seconds or --max-gb GB;
such a region is recorded as 'timeout' or 'memory' and has no alignment. Alignments are written to
work/hardset/<contig>/msa/<strategy>/<id>.fa.gz (one row per allele, '-' for gaps, allele order).

Strategies (STRATEGIES below):
  star          the current method: pgrealign's medoid star with abPOA (star.star_msa, its SCORES)
  star_longest  the same star, centred on the longest allele instead of the medoid
  star_gapnorm  the star's alignment with every row's gaps shifted left against the column consensus
                (tools/gap_norm.py), as the optima are left-normalised
  abpoa_mc      progressive abPOA (-p guide tree) with Minigraph-Cactus BAR's settings: HOXD70-like matrix,
                gaps -O 400,1200 -E 30,1, band -b 1000 -f 0.1, -k 15 -w 5 -n 500, longest first
                (tools/realign_poa_mc.py; BAR's 10 kb end rule is not applied)
  abpoa_prog    progressive abPOA with pgrealign's SCORES and the same band
  famsa         FAMSA 2 defaults
  mafft         mafft FFT-NS-2 (--retree 2 --maxiterate 0)
  cluster_par   k-mer clusters (greedy, sketch distance < CLUSTER_D to the cluster's heaviest allele), a
                medoid star inside each, joined as parallel alternatives: the alleles' common exact
                prefix and suffix are shared columns, each cluster's middles get their own columns
  cluster_merge the same clusters and stars, joined into one alignment by mafft --merge (each cluster's
                alignment kept fixed)
  select        not run: per region, whichever of SELECT has the best pair F1 on 20 other allele pairs
                (selection_pairs, `optima` caches their optima too), scored on qc's pairs; its time is
                the sum of the candidates'
  units         the repeat-unit-aware aligner (tools/realign_units.py, M4): alleles cut into catalogue-motif
                units, unit strings aligned progressively, expanded and polished; it falls back to mafft
                when the motif does not fit or its alignment costs > 1.2x the pairwise optima
  star_centrolign     the medoid star with each allele's alignment to the centre from centrolign
                      (jeizenga/centrolign, two sequences in, a CIGAR out; work/tools/centrolign or
                      $VNTR_CENTROLIGN) instead of abPOA, left-normalised and merged exactly as `star`
                      (abPOA still aligns the insertion slots); pairs run $VNTR_CENTROLIGN_JOBS (3) at a time
  star_centrolign_c0  the same with centrolign's constraint_method 0: its default (3) calls a segment
                      with too little anchor score unalignable and leaves it as a deletion plus an insertion

Scoring (`score`): qc.compared_pairs picks the same pairs as `pgrealign qc` (20 random allele pairs, and
CHM13's allele against the 20 heaviest), whose optimal alignments `optima` computed once with
star.align_pair under SCORES (banded above qc's max_cells; a pair whose alignment would exceed
--max-gb is skipped, as qc does). Complexity of an alignment's column graph: `columns`; `sequence`, the
bases of the graph (qc's definition); `nodes` and `edges` once unary chains are merged. `cost_ratio`:
the summed cost of the induced allele pairs over that of their optima (pair_cost), which, unlike F1, does
not penalise a co-optimal gap placement.
"""
import argparse
import collections
import csv
import gzip
import json
import os
import random
import re
import shutil
import signal
import statistics
import subprocess
import sys
import tempfile
import threading
import time
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from hardset import PGREALIGN_DIR, CATALOG  # noqa: E402,F401  (hardset puts pgrealign on sys.path)
from pgrealign import qc, star  # noqa: E402
from pgrealign.realign import MAX_CELLS  # noqa: E402

CONTIG = 'chr20'
ABPOA = config._tool('VNTR_ABPOA', 'abpoa')
FAMSA = config._tool('VNTR_FAMSA', 'famsa')
CENTROLIGN = config._tool('VNTR_CENTROLIGN', 'centrolign',
                          os.path.join(config.WORK_DIR, 'tools', 'centrolign', 'build', 'centrolign'))
CENTROLIGN_PAIR_JOBS = int(os.environ.get('VNTR_CENTROLIGN_JOBS', 3))
CIGAR_RE = re.compile(r'(\d+)([=XMID])')
CLUSTER_D = 0.5
MC_MATRIX = [91, -114, -61, -123, -100, -114, 100, -125, -61, -100, -61, -125, 100, -114, -100,
             -123, -61, -114, 91, -100, -100, -100, -100, -100, 100]
SELECT = ['star', 'star_longest', 'star_gapnorm']     # `select` chooses among these, first on ties
STRATEGIES = ['star', 'star_longest', 'star_gapnorm', 'abpoa_mc', 'abpoa_prog', 'famsa', 'mafft',
              'cluster_par', 'cluster_merge', 'units', 'star_centrolign', 'star_centrolign_c0']


def paths(contig):
    root = os.path.join(config.WORK_DIR, 'hardset', contig)
    return {'root': root, 'fasta': os.path.join(root, 'fasta'), 'optima': os.path.join(root, 'optima'),
            'optima_sel': os.path.join(root, 'optima_sel'),
            'msa': os.path.join(root, 'msa'), 'log': os.path.join(root, 'runs.tsv'),
            'list': os.path.join(config.REPO, 'hardset', contig + '.tsv')}


def read_list(contig):
    with open(paths(contig)['list']) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def read_alleles(path):
    """[(id, seq, weight)] from a hardset FASTA."""
    out, name, w, chunks = [], None, 1, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ''.join(chunks), w))
                x = line[1:].split()
                name, chunks = x[0], []
                w = int(next((t.split('=')[1] for t in x[1:] if t.startswith('weight=')), 1))
            else:
                chunks.append(line.strip().upper())
    if name is not None:
        out.append((name, ''.join(chunks), w))
    return out


def write_fasta(path, recs):
    with open(path, 'w') as f:
        for n, s in recs:
            f.write('>%s\n%s\n' % (n, s))


def run_tool(cmd, out_path=None, cwd=None):
    """Run cmd, its stdout to out_path, in the directory of its last argument (the input, in the work
    directory), since abPOA -p and FAMSA leave files such as GuideTree behind."""
    cwd = cwd or os.path.dirname(os.path.abspath(cmd[-1]))
    with open(out_path, 'w') if out_path else open(os.devnull, 'w') as o:
        p = subprocess.run(cmd, stdout=o, stderr=subprocess.PIPE, cwd=cwd, env=config.tool_env())
    if p.returncode:
        raise RuntimeError('%s exited %d: %s' % (cmd[0], p.returncode, p.stderr.decode(errors='replace')[-300:]))


def check_rows(rows, alleles):
    width = {len(r) for r in rows.values()}
    if len(width) != 1:
        raise ValueError('rows of unequal length')
    for a, s, _ in alleles:
        if rows[a].replace('-', '') != s:
            raise ValueError('row %s does not spell its allele' % a)
    return rows


# ----------------------------------------------------------------------------- strategies

def s_star(alleles, wd, ref, centre=None):
    with star.Abpoa(ABPOA, star.SCORES, workdir=wd) as align:
        return star.star_msa(alleles, align=align, max_cells=MAX_CELLS, centre=centre).rows


def s_star_longest(alleles, wd, ref):
    return s_star(alleles, wd, ref, centre=max(alleles, key=lambda a: (len(a[1]), a[2]))[0])


def s_star_gapnorm(alleles, wd, ref, region=None):
    """The star's alignment (from `run star` when there), every row's gap runs shifted left against the
    column consensus to a fixed point (tools/gap_norm.py)."""
    import gap_norm
    path = os.path.join(paths(region['contig'])['msa'], 'star', region['id'] + '.fa.gz')
    rows = read_msa(path) if os.path.exists(path) else s_star(alleles, wd, ref)
    names = [a for a, _, _ in alleles]
    out, _ = gap_norm.normalise([rows[a] for a in names], 'left')
    return dict(zip(names, out))


def abpoa_msa(alleles, wd, flags, longest_first=False):
    order = sorted(alleles, key=lambda a: -len(a[1])) if longest_first else list(alleles)
    fa, out = os.path.join(wd, 'in.fa'), os.path.join(wd, 'out.fa')
    write_fasta(fa, [(a, s) for a, s, _ in order if s])
    run_tool([ABPOA, '-m', '0', '-r', '1'] + flags + [fa], out)
    rows = {k: v.upper() for k, v in star.read_fasta_rows(out).items()}
    width = len(next(iter(rows.values())))
    for a, s, _ in alleles:
        rows.setdefault(a, '-' * width)
    return rows


def s_abpoa_mc(alleles, wd, ref):
    mtx = os.path.join(wd, 'mc.mtx')
    with open(mtx, 'w') as f:
        f.write('\tA\tC\tG\tT\tN\n')
        for i in range(5):
            f.write('ACGTN'[i] + ''.join('\t%d' % MC_MATRIX[i * 5 + j] for j in range(5)) + '\n')
    return abpoa_msa(alleles, wd, ['-O', '400,1200', '-E', '30,1', '-b', '1000', '-f', '0.1', '-t', mtx, '-p',
                                   '-k', '15', '-w', '5', '-n', '500'], longest_first=True)


def s_abpoa_prog(alleles, wd, ref):
    return abpoa_msa(alleles, wd, star.SCORES.flags() + ['-b', '1000', '-f', '0.1', '-p'])


def msa_tool(alleles, wd, cmd):
    fa, out = os.path.join(wd, 'in.fa'), os.path.join(wd, 'out.fa')
    write_fasta(fa, [(a, s) for a, s, _ in alleles if s])
    run_tool(cmd(fa, out))
    rows = {k: v.upper() for k, v in star.read_fasta_rows(out).items()}
    width = len(next(iter(rows.values())))
    for a, s, _ in alleles:
        rows.setdefault(a, '-' * width)
    return rows


def s_famsa(alleles, wd, ref):
    return msa_tool(alleles, wd, lambda fa, out: [FAMSA, '-t', '2', fa, out])


def s_mafft(alleles, wd, ref):
    fa, out = os.path.join(wd, 'in.fa'), os.path.join(wd, 'out.fa')
    write_fasta(fa, [(a, s) for a, s, _ in alleles if s])
    run_tool([config.MAFFT, '--retree', '2', '--maxiterate', '0', '--nuc', '--thread', '2', '--quiet', fa], out)
    rows = {k: v.upper() for k, v in star.read_fasta_rows(out).items()}
    width = len(next(iter(rows.values())))
    for a, s, _ in alleles:
        rows.setdefault(a, '-' * width)
    return rows


def clusters(alleles, d=CLUSTER_D):
    """Greedy k-mer clusters: alleles by weight, heaviest first, each joins the cluster whose founder is
    nearest by sketch distance when that is below d, or founds a new cluster. [[allele index]]."""
    sk = [star.sketch(star.kmer_counts(s)) for _, s, _ in alleles]
    founders, members = [], []
    for i in sorted(range(len(alleles)), key=lambda i: (-alleles[i][2], i)):
        best = min(((star.sketch_distance(sk[i], sk[f]), c) for c, f in enumerate(founders)), default=(1.0, -1))
        if best[0] < d:
            members[best[1]].append(i)
        else:
            founders.append(i)
            members.append([i])
    return members


def cluster_stars(alleles, wd, seqs=None):
    """A medoid star of each cluster, over `seqs` (default the alleles' own sequences)."""
    seqs = seqs or {a: s for a, s, _ in alleles}
    out = []
    with star.Abpoa(ABPOA, star.SCORES, workdir=wd) as align:
        for c in clusters(alleles):
            sub = [(alleles[i][0], seqs[alleles[i][0]], alleles[i][2]) for i in c]
            if all(not s for _, s, _ in sub):
                out.append({a: '' for a, _, _ in sub})
            elif len(sub) == 1:
                out.append({sub[0][0]: sub[0][1]})
            else:
                out.append(star.star_msa(sub, align=align, max_cells=MAX_CELLS).rows)
    return out


def s_cluster_par(alleles, wd, ref):
    ss = [s for _, s, _ in alleles]
    short = min(len(s) for s in ss)
    p = 0
    while p < short and all(s[p] == ss[0][p] for s in ss):
        p += 1
    q = 0
    while q < short - p and all(s[-1 - q] == ss[0][-1 - q] for s in ss):
        q += 1
    mids = {a: s[p:len(s) - q] for a, s, _ in alleles}
    blocks = cluster_stars(alleles, wd, mids)
    widths = [len(next(iter(b.values()))) for b in blocks]
    rows = {}
    for k, b in enumerate(blocks):
        for a, r in b.items():
            rows[a] = ss[0][:p] + ''.join(r if j == k else '-' * widths[j] for j in range(len(blocks))) + \
                (ss[0][len(ss[0]) - q:] if q else '')
    return rows


def s_cluster_merge(alleles, wd, ref):
    blocks = cluster_stars(alleles, wd)
    if len(blocks) == 1:
        return blocks[0]
    recs, table, n = [], [], 0
    for b in blocks:
        idx = []
        for a, r in b.items():
            if r.replace('-', ''):
                recs.append((a, r))
                n += 1
                idx.append(n)
        if len(idx) > 1:
            table.append(' '.join(map(str, idx)))
    fa, tab, out = os.path.join(wd, 'in.fa'), os.path.join(wd, 'sub.table'), os.path.join(wd, 'out.fa')
    write_fasta(fa, recs)
    with open(tab, 'w') as f:
        f.write('\n'.join(table) + '\n')
    run_tool([config.MAFFT, '--merge', tab, '--nuc', '--thread', '2', '--quiet', fa], out)
    rows = {k: v.upper() for k, v in star.read_fasta_rows(out).items()}
    width = len(next(iter(rows.values())))
    for a, s, _ in alleles:
        rows.setdefault(a, '-' * width)
    return rows


def s_units(alleles, wd, ref, region=None):
    import realign_units
    rj = {'region_id': region['id'], 'stratum': 'hotspot_vntr', 'span_start': int(region['ref_start']),
          'core_start': int(region['core_start']), 'core_end': int(region['core_end']),
          'period': region['period'], 'motif': region['motif'],
          'tr_annotation': {'period': region['period'], 'tr_class': region['classes'].split(',')[-1]}}
    name = {a: ('CHM13#' + a if a == ref else a) for a, _, _ in alleles}
    fa, out = os.path.join(wd, 'in.fa'), os.path.join(wd, 'out.fa')
    write_fasta(fa, [(name[a], s) for a, s, _ in alleles])
    info = realign_units.align_units_fasta(fa, rj, out, threads=2, workdir=wd)
    if info.get('status') != 'ok':
        raise RuntimeError('units: %s' % (info.get('message') or info.get('reason')))
    got = {k: v.upper() for k, v in star.read_fasta_rows(out).items()}
    sys.stderr.write('units mode=%s reason=%s\n' % (info.get('mode'), info.get('reason')))
    return {a: got[name[a]] for a, _, _ in alleles}


def cigar_rows(a, b, cigar):
    """The two rows a CIGAR spells, a being the sequence that D consumes (centrolign's first input)."""
    ra, rb, i, j = [], [], 0, 0
    for n, op in CIGAR_RE.findall(cigar):
        n = int(n)
        if op in '=XM':
            ra.append(a[i:i + n])
            rb.append(b[j:j + n])
            i, j = i + n, j + n
        elif op == 'D':
            ra.append(a[i:i + n])
            rb.append('-' * n)
            i += n
        else:
            ra.append('-' * n)
            rb.append(b[j:j + n])
            j += n
    if i != len(a) or j != len(b):
        raise ValueError('CIGAR spells %d/%d of %d/%d bases' % (i, j, len(a), len(b)))
    return ''.join(ra), ''.join(rb)


def centrolign_pair(centre, seq, wd, conf=None):
    """centrolign's alignment of seq to centre (two sequences in, a CIGAR out), as (centre row, row),
    left-normalised as star.align_pair does. `conf` is a config file (-C) holding every parameter."""
    if seq == centre:
        return centre, seq
    if not centre or not seq:
        return centre + '-' * len(seq), '-' * len(centre) + seq
    fa = os.path.join(wd, 'pair.fa')
    write_fasta(fa, [('centre', centre), ('seq', seq)])
    if conf:
        cfg = os.path.join(wd, 'pair.yaml')
        with open(conf) as f, open(cfg, 'w') as o:
            o.write(f.read().replace('fasta_name: ""', 'fasta_name: "%s"' % fa))
        cmd = [CENTROLIGN, '-C', cfg]
    else:
        cmd = [CENTROLIGN, '-v', '0', fa]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=wd)
    if p.returncode:
        raise RuntimeError('centrolign exited %d: %s' % (p.returncode, p.stderr[-300:]))
    ca, sa = cigar_rows(centre, seq, p.stdout.strip())
    return star.left_normalise_pair(ca, sa)


def s_star_centrolign(alleles, wd, ref, conf=None):
    """The medoid star, each allele's alignment to the centre from centrolign instead of abPOA; abPOA
    still aligns the insertion slots, and star_msa merges exactly as in `star`. The pairs run
    CENTROLIGN_PAIR_JOBS at a time (centrolign is ~40x abPOA's time), and their summed CPU time goes to
    stderr as cpu_seconds."""
    import concurrent.futures
    import resource
    seq = {a: s if set(s) <= star._ALIGNABLE else ''.join(c if c in star._ALIGNABLE else 'N' for c in s)
           for a, s, _ in alleles}
    names = [a for a, _, _ in alleles]
    cid = names[star.medoid([seq[a] for a in names], [w for _, _, w in alleles])]

    def one(i, s):
        d = os.path.join(wd, 'p%d' % i)
        os.makedirs(d, exist_ok=True)
        return s, centrolign_pair(seq[cid], s, d, conf)

    with concurrent.futures.ThreadPoolExecutor(CENTROLIGN_PAIR_JOBS) as ex:
        done = dict(ex.map(lambda x: one(*x), enumerate(seq[a] for a in names if a != cid)))
    orig = star.align_pair
    star.align_pair = lambda centre, s, align, banded=False: (centre, s) if s == centre else done[s]
    try:
        rows = s_star(alleles, wd, ref, centre=cid)
    finally:
        star.align_pair = orig
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    sys.stderr.write('cpu_seconds=%.1f centre=%s\n' % (ru.ru_utime + ru.ru_stime, cid))
    return rows


def s_star_centrolign_c0(alleles, wd, ref):
    """star_centrolign with constraint_method 0: every segment counts as alignable, so no stretch is
    left as an unaligned deletion plus insertion for want of anchor score."""
    conf = os.path.join(wd, 'c0.yaml')
    text = subprocess.run([CENTROLIGN, '-G'], capture_output=True, text=True).stdout
    if 'constraint_method: 3' not in text:
        raise RuntimeError('centrolign -G has no constraint_method: 3 to change')
    with open(conf, 'w') as f:
        f.write(text.replace('constraint_method: 3', 'constraint_method: 0').replace(
            'logging_level: 2', 'logging_level: 0'))
    return s_star_centrolign(alleles, wd, ref, conf)


FUNCS = {'star': s_star, 'star_longest': s_star_longest, 'abpoa_mc': s_abpoa_mc, 'abpoa_prog': s_abpoa_prog,
         'famsa': s_famsa, 'mafft': s_mafft, 'cluster_par': s_cluster_par, 'cluster_merge': s_cluster_merge,
         'units': s_units, 'star_gapnorm': s_star_gapnorm, 'star_centrolign': s_star_centrolign,
         'star_centrolign_c0': s_star_centrolign_c0}


def core_coords(region):
    """Catalogue start and end of the region's longest target (whose motif hardset.py lists), for the
    units aligner."""
    want = set(region['targets'].split(','))
    best = (0, None, None)
    with gzip.open(CATALOG, 'rt') as f:
        for line in f:
            if line.startswith(region['contig'] + '\t'):
                x = line.split('\t')
                if x[3] in want:
                    best = max(best, (int(x[2]) - int(x[1]), int(x[1]), int(x[2])))
    return best[1], best[2]


def cmd_align(a):
    """Child process: one strategy on one region (or `optima`)."""
    P = paths(a.contig)
    region = next(r for r in read_list(a.contig) if r['id'] == a.region)
    alleles = read_alleles(os.path.join(P['fasta'], a.region + '.fa'))
    os.makedirs(os.path.join(config.WORK_DIR, 'tmp'), exist_ok=True)
    wd = tempfile.mkdtemp(prefix='hs.', dir=os.path.join(config.WORK_DIR, 'tmp'))
    try:
        if a.strategy in ('optima', 'optima_sel'):
            return optima(region, alleles, a.out, wd, a.max_gb, select=a.strategy == 'optima_sel')
        kw = {}
        if a.strategy == 'units':
            region['core_start'], region['core_end'] = core_coords(region)
            kw['region'] = region
        elif a.strategy == 'star_gapnorm':
            kw['region'] = region
        rows = check_rows(FUNCS[a.strategy](alleles, wd, region['ref_allele'], **kw), alleles)
        rows = star.drop_gap_columns({k: rows[k] for k, _, _ in alleles})
        with gzip.open(a.out + '.tmp', 'wt') as f:
            for k, _, _ in alleles:
                f.write('>%s\n%s\n' % (k, rows[k]))
        os.replace(a.out + '.tmp', a.out)
    finally:
        shutil.rmtree(wd, ignore_errors=True)


def selection_pairs(pkg, rows, n=20):
    """Up to n random allele pairs (CHM13's allele included) that qc.compared_pairs does not compare,
    drawn with a seed of their own: what `select` chooses a region's alignment on, so the choice is
    not scored on the pairs that made it."""
    used = {(a, b) for _, a, b in qc.compared_pairs(pkg, rows)}
    ids = [a['id'] for a in pkg['alleles']]
    cand = [(ids[i], ids[j]) for i in range(len(ids)) for j in range(i + 1, len(ids))
            if (ids[i], ids[j]) not in used and (ids[j], ids[i]) not in used]
    rng = random.Random(zlib.crc32((pkg['id'] + '/select').encode()))
    return [('pair', a, b) for a, b in sorted(rng.sample(cand, min(n, len(cand))))]


def optima(region, alleles, out, wd, max_gb, select=False):
    """The optimal alignments of the pairs qc compares (or, with `select`, of selection_pairs), as JSON
    lines (kind, a, b, row a, row b)."""
    pkg = {'id': region['id'], 'ref': region['ref_allele'],
           'alleles': [{'id': a, 'seq': s, 'weight': w} for a, s, w in alleles]}
    seqs = {a: s for a, s, _ in alleles}
    rows = dict.fromkeys(seqs)
    limit = max_gb * 2 ** 30
    with star.Abpoa(ABPOA, star.SCORES, workdir=wd) as align, gzip.open(out + '.tmp', 'wt') as f:
        for kind, x, y in (selection_pairs if select else qc.compared_pairs)(pkg, rows):
            n, m = len(seqs[x]), len(seqs[y])
            rec = {'kind': kind, 'a': x, 'b': y}
            if qc.pair_bytes(n, m, MAX_CELLS) > limit:
                rec['skipped'] = True
            else:
                rec['banded'] = n * m > MAX_CELLS
                rec['rows'] = star.align_pair(seqs[x], seqs[y], align, banded=rec['banded'])
            f.write(json.dumps(rec) + '\n')
    os.replace(out + '.tmp', out)


# ----------------------------------------------------------------------------- the job runner

def group_rss_kb(pgid):
    p = subprocess.run(['ps', '-A', '-o', 'pgid=,rss='], capture_output=True, text=True)
    return sum(int(r) for g, r in (line.split() for line in p.stdout.splitlines() if line.strip())
               if int(g) == pgid)


def run_job(contig, strategy, region, out, timeout, max_gb):
    cmd = [sys.executable, os.path.abspath(__file__), 'align', strategy, region, out, '--contig', contig,
           '--max-gb', str(max_gb)]
    errp = out + '.err'
    t0, peak, status = time.time(), 0, 'ok'
    with open(errp, 'w') as err:
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=err, start_new_session=True,
                             env=config.tool_env())
        while p.poll() is None:
            time.sleep(1)
            rss = group_rss_kb(p.pid)
            peak = max(peak, rss)
            if time.time() - t0 > timeout or rss > max_gb * 2 ** 20:
                status = 'timeout' if time.time() - t0 > timeout else 'memory'
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
    if status == 'ok' and p.returncode:
        status = 'error'
    return {'strategy': strategy, 'id': region, 'status': status, 'seconds': round(time.time() - t0, 1),
            'peak_gb': round(peak / 2 ** 20, 2)}


def cmd_run(a):
    P = paths(a.contig)
    regions = read_list(a.contig)
    if a.regions:
        regions = [r for r in regions if r['id'] in set(a.regions)]
    regions.sort(key=lambda r: int(r['alleles']) * int(r['longest']))
    done = set()
    if os.path.exists(P['log']):
        with open(P['log']) as f:
            done = {(r['strategy'], r['id']) for r in csv.DictReader(f, delimiter='\t')}
    jobs = []                        # smallest region first, every strategy on it
    for r in regions:
        for s in a.strategies:
            d = P.get(s) or os.path.join(P['msa'], s)
            os.makedirs(d, exist_ok=True)
            if (s, r['id']) not in done or a.force:
                jobs.append((s, r['id'], os.path.join(d, r['id'] + ('.jsonl.gz' if s in P else '.fa.gz'))))
    lock = threading.Lock()
    new = not os.path.exists(P['log'])
    log = open(P['log'], 'a')
    if new:
        log.write('strategy\tid\tstatus\tseconds\tpeak_gb\n')
        log.flush()
    todo = list(jobs)

    def worker():
        while True:
            with lock:
                if not todo:
                    return
                s, rid, out = todo.pop(0)
            res = run_job(a.contig, s, rid, out, a.timeout, a.max_gb)
            with lock:
                log.write('\t'.join(str(res[k]) for k in ('strategy', 'id', 'status', 'seconds', 'peak_gb')) + '\n')
                log.flush()
                print('%s %s %s %ss %sGB' % (s, rid, res['status'], res['seconds'], res['peak_gb']), flush=True)

    threads = [threading.Thread(target=worker) for _ in range(a.jobs)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    log.close()


# ----------------------------------------------------------------------------- scoring

def complexity(rows):
    """columns, sequence (distinct bases per column, summed), and nodes and edges of the column graph
    once unary chains are merged."""
    names = list(rows)
    width = len(rows[names[0]])
    seq = qc.sequence_bp(rows)
    edges = set()
    for r in rows.values():
        prev = None
        for c, x in enumerate(r):
            if x != '-':
                v = (c, x)
                if prev is not None:
                    edges.add((prev, v))
                prev = v
    outd, ind = collections.Counter(), collections.Counter()
    for u, v in edges:
        outd[u] += 1
        ind[v] += 1
    merged = sum(1 for u, v in edges if outd[u] == 1 and ind[v] == 1)
    return {'columns': width, 'sequence': seq, 'nodes': seq - merged, 'edges': len(edges) - merged}


def read_msa(path):
    out, name = {}, None
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].strip()
            else:
                out[name] = line.strip()
    return out


def pair_cost(ra, rb, sc=star.SCORES):
    """The cost of a two-row alignment under abPOA's scores, as M (n + m) - 2 score: 0 for a match column,
    2 (M + X) for a mismatch, M per gapped base plus twice the gap run's cost. Minimising it maximises
    the score, so it does not depend on where co-optimal gaps are put, as the counts of aligned pairs do."""
    cost, run, side = 0, 0, None

    def gap(n):
        return 2 * min(sc.gap_open[0] + n * sc.gap_ext[0], sc.gap_open[1] + n * sc.gap_ext[1]) + sc.match * n

    for x, y in zip(ra, rb):
        if x == '-' and y == '-':
            continue
        g = 'a' if x == '-' else 'b' if y == '-' else None
        if g != side and run:
            cost += gap(run)
            run = 0
        side = g
        if g:
            run += 1
        elif x != y:
            cost += 2 * (sc.match + sc.mismatch)
    if run:
        cost += gap(run)
    return cost


def score_region(rows, opt):
    s = dict.fromkeys(['pair_tp', 'pair_opt', 'pair_msa', 'ref_tp', 'ref_opt', 'ref_msa', 'skipped',
                       'pair_cost', 'pair_cost_opt'], 0)
    for rec in opt:
        if rec.get('skipped'):
            s['skipped'] += 1
            continue
        o = qc.aligned_pairs(*rec['rows'])
        m = qc.aligned_pairs(*qc.induced(rows[rec['a']], rows[rec['b']]))
        k = rec['kind']
        s[k + '_tp'] += len(o & m)
        s[k + '_opt'] += len(o)
        s[k + '_msa'] += len(m)
        if k == 'pair':
            s['pair_cost'] += pair_cost(*qc.induced(rows[rec['a']], rows[rec['b']]))
            s['pair_cost_opt'] += pair_cost(*rec['rows'])
    return s


def f1(tp, opt, msa):
    return 2 * tp / (opt + msa) if opt + msa else None


def cmd_score(a):
    P = paths(a.contig)
    regions = read_list(a.contig)
    runs = {}
    with open(P['log']) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            runs[(r['strategy'], r['id'])] = r
    strategies = [s for s in STRATEGIES if any(k[0] == s for k in runs)]
    out = []
    for reg in regions:
        rid = reg['id']
        with gzip.open(os.path.join(P['optima'], rid + '.jsonl.gz'), 'rt') as f:
            opt = [json.loads(line) for line in f]
        with gzip.open(os.path.join(P['optima_sel'], rid + '.jsonl.gz'), 'rt') as f:
            sel = [json.loads(line) for line in f]
        first = len(out)
        for s in strategies:
            run = runs.get((s, rid))
            row = {'strategy': s, 'id': rid, 'status': run['status'] if run else 'not_run',
                   'seconds': run['seconds'] if run else '', 'peak_gb': run['peak_gb'] if run else ''}
            path = os.path.join(P['msa'], s, rid + '.fa.gz')
            if row['status'] == 'ok' and os.path.exists(path):
                rows = read_msa(path)
                sc = score_region(rows, opt)
                ss = score_region(rows, sel)
                v = f1(ss['pair_tp'], ss['pair_opt'], ss['pair_msa'])
                row['sel_f1'] = round(v, 4) if v is not None else ''      # '': no pair left to choose on
                row.update(sc)
                row.update(complexity(rows))
                row['pair_f1'] = round(f1(sc['pair_tp'], sc['pair_opt'], sc['pair_msa']), 4)
                row['ref_f1'] = round(f1(sc['ref_tp'], sc['ref_opt'], sc['ref_msa']), 4) if sc['ref_opt'] else ''
                row['cost_ratio'] = round(sc['pair_cost'] / sc['pair_cost_opt'], 4) if sc['pair_cost_opt'] else ''
            out.append(row)
            print(rid, s, row.get('pair_f1'), file=sys.stderr)
        cands = [r for r in out[first:] if r['strategy'] in SELECT and 'sel_f1' in r]
        if cands:
            best = max(cands, key=lambda r: (r['sel_f1'] if r['sel_f1'] != '' else -1, -SELECT.index(r['strategy'])))
            out.append(dict(best, strategy='select', chosen=best['strategy'],
                            seconds=round(sum(float(r['seconds']) for r in out[first:] if r['strategy'] in SELECT), 1),
                            peak_gb=max(float(r['peak_gb']) for r in out[first:] if r['strategy'] in SELECT)))
    fields = ['strategy', 'id', 'status', 'seconds', 'peak_gb', 'pair_f1', 'ref_f1', 'cost_ratio', 'sel_f1', 'chosen', 'pair_tp',
              'pair_opt', 'pair_msa', 'ref_tp', 'ref_opt', 'ref_msa', 'pair_cost', 'pair_cost_opt', 'skipped',
              'columns', 'sequence', 'nodes', 'edges']
    path = a.out or os.path.join(config.RESULTS_DIR, 'hardset_%s.tsv' % a.contig)
    with open(path, 'w') as f:
        f.write('\t'.join(fields) + '\n')
        for r in out:
            f.write('\t'.join(str(r.get(k, '')) for k in fields) + '\n')
    print('wrote', path, file=sys.stderr)
    summary(out, regions, path[:-4] + '.summary.tsv')


def quartile(xs):
    xs = sorted(xs)
    return xs[(len(xs) - 1) // 4] if xs else None


def summary(out, regions, path, eps=0.005):
    """One row per strategy: regions aligned, median and lower-quartile pair F1 and ref F1 (a region it
    did not align counts at the original graph's F1, `pair_f1_graph`, which is what would stay), regions
    better and worse than `star` by more than eps, summed run time, largest peak memory, and median
    graph nodes and sequence relative to star's."""
    graph = {r['id']: (float(r['pair_f1_graph']), float(r['ref_f1_graph'] or 'nan')) for r in regions}
    by = collections.defaultdict(dict)
    for r in out:
        by[r['strategy']][r['id']] = r
    base = by.get('star', {})
    fields = ['strategy', 'aligned', 'pair_f1_median', 'pair_f1_q1', 'ref_f1_median', 'better', 'worse',
              'cost_ratio_median', 'seconds', 'max_peak_gb', 'nodes_vs_star', 'sequence_vs_star']
    rows = []
    for s, res in [('graph', None)] + list(by.items()):
        if res is None:
            pf = [graph[i][0] for i in graph]
            rf = [graph[i][1] for i in graph]
            rows.append({'strategy': 'original graph', 'aligned': len(graph), 'pair_f1_median': statistics.median(pf),
                         'pair_f1_q1': quartile(pf), 'ref_f1_median': statistics.median(rf)})
            continue
        ok = {i: r for i, r in res.items() if 'pair_f1' in r}
        pf = [ok[i]['pair_f1'] if i in ok else graph[i][0] for i in graph]
        rf = [ok[i]['ref_f1'] if i in ok and ok[i]['ref_f1'] != '' else graph[i][1] for i in graph]
        both = [i for i in ok if i in base and 'pair_f1' in base[i]]
        rows.append({
            'strategy': s, 'aligned': len(ok), 'pair_f1_median': statistics.median(pf), 'pair_f1_q1': quartile(pf),
            'ref_f1_median': statistics.median(rf),
            'cost_ratio_median': statistics.median(ok[i]['cost_ratio'] for i in ok if ok[i]['cost_ratio'] != ''),
            'better': sum(ok[i]['pair_f1'] > base[i]['pair_f1'] + eps for i in both),
            'worse': sum(ok[i]['pair_f1'] < base[i]['pair_f1'] - eps for i in both),
            'seconds': round(sum(float(r['seconds'] or 0) for r in res.values())),
            'max_peak_gb': max((float(r['peak_gb'] or 0) for r in res.values()), default=0),
            'nodes_vs_star': statistics.median(ok[i]['nodes'] / base[i]['nodes'] for i in both) if both else None,
            'sequence_vs_star': statistics.median(ok[i]['sequence'] / base[i]['sequence'] for i in both) if both else None})
    with open(path, 'w') as f:
        f.write('\t'.join(fields) + '\n')
        for r in rows:
            f.write('\t'.join('' if r.get(k) is None else ('%.4f' % r[k] if isinstance(r[k], float) else str(r[k]))
                              for k in fields) + '\n')
    print('wrote', path, file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('optima')
    p.add_argument('regions', nargs='*')
    p = sub.add_parser('run')
    p.add_argument('strategies', nargs='+', choices=STRATEGIES)
    p.add_argument('--regions', nargs='*')
    p = sub.add_parser('align')
    p.add_argument('strategy', choices=STRATEGIES + ['optima', 'optima_sel'])
    p.add_argument('region')
    p.add_argument('out')
    p = sub.add_parser('score')
    p.add_argument('-o', '--out')
    for p in sub.choices.values():
        p.add_argument('--contig', default=CONTIG)
        p.add_argument('--jobs', type=int, default=3)
        p.add_argument('--timeout', type=float, default=1800)
        p.add_argument('--max-gb', type=float, default=9.0)
        p.add_argument('--force', action='store_true')
    a = ap.parse_args(argv)
    if a.cmd == 'optima':
        a.strategies, a.regions = ['optima', 'optima_sel'], a.regions or None
        cmd_run(a)
    elif a.cmd == 'run':
        cmd_run(a)
    elif a.cmd == 'align':
        cmd_align(a)
    else:
        cmd_score(a)


if __name__ == '__main__':
    main()
