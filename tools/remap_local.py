#!/usr/bin/env python3
"""remap_local.py -- Stage 2: re-map a region's short reads to local graphs and compare them.

For one region package (regions/<id>) and a set of local graphs of it (the Minigraph-Cactus
baseline mc.gfa, realigned candidates candidates/<method>/<id>.gfa, and a "truth graph" built
from HG002's two haplotypes), this tool

  1. fetches the reads the genome-wide mapping placed in or next to the region, with their
     sequences and base qualities (GAF-Base: the read is rebuilt from the path it aligned to and
     its cs difference string; bq:Z holds the qualities), keeping mates together;
  2. extends every local graph by the same CHM13 flank on both sides (default 1 kb, one node
     each, joined to every path), so reads that run past the anchors and mates in the flanks
     map the same way whatever the graph between the anchors is;
  3. indexes each graph (vg autoindex -w sr-giraffe) and maps the reads with vg giraffe
     (paired where both mates were recovered, with one fixed fragment-length distribution for
     every graph; single-end otherwise);
  4. scores each mapping: fraction mapped, MAPQ distribution, core/flank depth, haplotype
     compatibility and spread, best-pair explanation, and against the truth: which reads fit
     HG002's haplotypes, whether those reads are placed confidently, and (where both truth
     haplotypes are exact copies of hap32 paths) whether they are placed at the right position.

Subcommands (run from the repository root):

    python3 tools/remap_local.py fetch   L014297 [L009656 ...]      # reads -> <work>/<id>/reads*
    python3 tools/remap_local.py fraglen [ids...]                   # pooled fragment length -> <work>/fraglen.json
    python3 tools/remap_local.py map     L014297 --methods mc,unit_aware,truth [--threads 2]
    python3 tools/remap_local.py analyse L014297 [--methods ...]    # -> <work>/<id>/stage2.json
    python3 tools/remap_local.py run     L014297 L009656 ... [--jobs 2 --threads 2]   # all of the above
    python3 tools/remap_local.py table   [ids...] --tsv results/stage2_remap.tsv      # one row per region x graph

--work defaults to $VNTR_WORK/remap (git-ignored). Methods default to DEFAULT_METHODS; a method
whose candidate GFA does not exist for a region is skipped and recorded. `truth` is the graph of
HG002's two haplotypes (mafft G-INS-i of truth.fa, then msa_graph.py) -- evaluation only: it
gives the ceiling a perfect graph for this sample would reach, and the reference placements.
`orig` (analysis only) is the genome-wide giraffe mapping to the whole hap32 graph as stored in
the GAF-Base, scored the same way with mc.gfa's paths.

Definitions (also in results/stage2_remap.md):

  core       each sequence's copy of the census core interval: the CHM13 flank just outside the
             core (up to 300 bp each side) is aligned semi-globally to the sequence, and the core is
             what lies between the two matches. Per graph, every node base gets the fraction of the
             paths through it that have that base inside their core; read depth sums the aligned
             graph bases weighted by that fraction. Depth = aligned bases / mean HG002 core
             (flank) length, so it is in x coverage.
  graph-core read   a mapped read with >= 10 aligned bases in the core of that graph.
  truth-core read   a read whose truth-graph alignment has >= 10 bases in HG002's core (a read set
             that does not depend on the graph being scored).
  compatible a read's aligned node walk is a contiguous sub-walk of some hap32 path (either strand).
  truth-consistent  the truth-graph alignment has at most 1 high-quality edit: substitutions at
             base quality >= 20, internal indel events, and terminal soft clips holding >= 5
             bases of quality >= 20 each count 1.
  placed correctly  the graph alignment lies on the hap32 path that is a copy of the truth
             haplotype the read comes from, within 5 bp of the truth-graph position (carried
             through a global alignment when the copy is near, edit <= 10, not exact) and on the
             same strand. Evaluated for truth-consistent reads with truth-graph MAPQ >= 30 whose
             truth haplotype(s) have such a copy.
"""
import argparse
import collections
import concurrent.futures
import gzip
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import time

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402

REPO = config.REPO
REGIONS = config.REGIONS_DIR
CANDIDATES = config.CANDIDATES_DIR
RESULTS = config.RESULTS_DIR
WORK = os.path.join(config.WORK_DIR, 'remap')

DEFAULT_METHODS = ['mc', 'mafft_linsi', 'unit_aware', 'poa_spoa', 'poa_abpoa',
                   'unit_aware__all', 'mafft_linsi__all', 'poa_spoa__all', 'poa_abpoa__all',
                   'poa_abpoa_mc__all', 'truth']
FLANK = 1000
CORE_READ_BP = 10
HQ = 20
TRUTH_MAPQ_UNIQUE = 30
POS_TOL = 5
NEAR_EDIT = 10
FLANK_PROBE = 300

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CS_RE = re.compile(r'(:\d+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+|=[A-Za-z]+)')
STEP_RE = re.compile(r'([<>])([^<>]+)')


def rc(s):
    return s.translate(COMP)[::-1]


def log(*a):
    print('[remap %s]' % time.strftime('%H:%M:%S'), *a, file=sys.stderr, flush=True)


def run(cmd, stdout=None, stderr=None, cwd=None, check=True):
    env = config.tool_env()
    so = open(stdout, 'w') if isinstance(stdout, str) else stdout
    se = open(stderr, 'w') if isinstance(stderr, str) else stderr
    try:
        p = subprocess.run([str(c) for c in cmd], stdout=so if so else subprocess.PIPE,
                           stderr=se if se else subprocess.PIPE, cwd=cwd, env=env, text=True)
    finally:
        if isinstance(stdout, str):
            so.close()
        if isinstance(stderr, str):
            se.close()
    if check and p.returncode != 0:
        msg = p.stderr[-2000:] if isinstance(p.stderr, str) else ''
        raise RuntimeError('command failed (%d): %s\n%s' % (p.returncode, ' '.join(map(str, cmd)), msg))
    return p


# ------------------------------------------------------------------ region / files

def region_dir(x):
    if os.path.isdir(x):
        return os.path.abspath(x)
    return os.path.join(REGIONS, x)


def load_region(x):
    rd = region_dir(x)
    with open(os.path.join(rd, 'region.json')) as f:
        d = json.load(f)
    d['dir'] = rd
    return d


def open_any(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_fasta(path):
    out, name, buf = [], None, []
    with open_any(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ''.join(buf).upper()))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        out.append((name, ''.join(buf).upper()))
    return out


def parse_gfa(path):
    """nodes {id: seq}, edges [(a, ao, b, bo)], paths [(name, [(id, '+'|'-')])] from S/L/P/W lines."""
    nodes, edges, paths = {}, [], []
    with open_any(path) as f:
        for line in f:
            if not line or line[0] not in 'SLPW':
                continue
            x = line.rstrip('\n').split('\t')
            if x[0] == 'S':
                nodes[x[1]] = x[2].upper()
            elif x[0] == 'L':
                edges.append((x[1], x[2], x[3], x[4]))
            elif x[0] == 'P':
                steps = [(s[:-1], s[-1]) for s in x[2].split(',') if s]
                paths.append((x[1], steps))
            elif x[0] == 'W':
                name = '%s#%s#%s' % (x[1], x[2], x[3])
                steps = [(n, '+' if o == '>' else '-') for o, n in STEP_RE.findall(x[6])]
                paths.append((name, steps))
    return nodes, edges, paths


def walk_seq(nodes, steps):
    return ''.join(nodes[n] if o == '+' else rc(nodes[n]) for n, o in steps)


def flip_steps(steps):
    return [(n, '-' if o == '+' else '+') for n, o in reversed(steps)]


def ref_flanks(reg, flank):
    """CHM13 sequence of `flank` bp left of span_start and right of span_end."""
    fa = config.data_paths(reg['contig'])['ref_fa']
    a, b = reg['span_start'], reg['span_end']
    lo = max(1, a - flank)
    p = run([config.SAMTOOLS, 'faidx', fa, '%s:%d-%d' % (reg['contig'], lo, a - 1),
             '%s:%d-%d' % (reg['contig'], b + 1, b + flank)])
    seqs, cur = [], None
    for line in p.stdout.splitlines():
        if line.startswith('>'):
            cur = []
            seqs.append(cur)
        else:
            cur.append(line.strip())
    left, right = (''.join(s).upper() for s in seqs)
    return left, right


# ------------------------------------------------------------------ 1. fetch reads

def parse_gaf_line(line):
    x = line.rstrip('\n').split('\t')
    if len(x) < 12:
        return None
    tags = {}
    for t in x[12:]:
        p = t.split(':', 2)
        if len(p) == 3:
            tags[p[0]] = p[2]
    path = x[5]
    steps = [] if path == '*' else [(n, '+' if o == '>' else '-') for o, n in STEP_RE.findall(path)]
    mate = 1 if 'fn' in tags else 2 if 'fp' in tags else 0
    return {'name': x[0], 'qlen': int(x[1]), 'qs': int(x[2]) if x[2] != '*' else 0,
            'qe': int(x[3]) if x[3] != '*' else 0, 'steps': steps,
            'plen': int(x[6]) if x[6] != '*' else 0, 'ps': int(x[7]) if x[7] != '*' else 0,
            'pe': int(x[8]) if x[8] != '*' else 0, 'mapq': int(x[11]) if x[11] != '*' else 0,
            'AS': int(tags['AS']) if 'AS' in tags else None, 'cs': tags.get('cs', ''),
            'bq': tags.get('bq'), 'mate': mate, 'raw_path': path, 'tags': tags}


def read_gaf(path):
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith('@') or not line.strip():
                continue
            r = parse_gaf_line(line)
            if r:
                out.append(r)
    return out


def read_from_cs(pseq, cs):
    """Aligned read bases, in the orientation of the path, from the path substring and cs."""
    out, i = [], 0
    for o in CS_RE.findall(cs):
        c = o[0]
        if c == ':':
            n = int(o[1:])
            out.append(pseq[i:i + n])
            i += n
        elif c == '=':
            out.append(o[1:].upper())
            i += len(o) - 1
        elif c == '*':
            out.append(o[2].upper())
            i += 1
        elif c == '+':
            out.append(o[1:].upper())
        elif c == '-':
            i += len(o) - 1
    if i != len(pseq):
        return None
    return ''.join(out)


def gbz_query(args, out_gfa):
    cmd = [config.GBZ_BASE, 'query'] + args + [config.data_paths()['gbz_db']]
    p = run(cmd, stdout=out_gfa)
    return p


def fetch(reg, wd, flank=FLANK):
    """Reads overlapping the local graph's extent: the anchor-to-anchor subgraph (--between) plus
    the CHM13 window span +- flank (--interval, context 100). Writes reads_1.fq, reads_2.fq
    (pairs), reads_se.fq, reads.tsv and reads.orig.gaf (the genome-wide alignments)."""
    os.makedirs(wd, exist_ok=True)
    gaf_db = config.data_paths()['gaf_db']
    c = reg['contig']
    a0, b0 = reg['span_start'] - 1, reg['span_end']          # 0-based half-open span
    L, R = reg['anchor_left'], reg['anchor_right']
    t0 = time.time()
    qs = []
    q1 = os.path.join(wd, 'q_between')
    gbz_query(['--between', '%s:%s' % (L, R), '--limit', '2000000', '--gaf-base', gaf_db,
               '--gaf-output', q1 + '.gaf', '--alignments', 'overlapping'], q1 + '.gfa')
    qs.append(q1)
    q2 = os.path.join(wd, 'q_window')
    gbz_query(['--sample', 'CHM13', '--contig', c, '--interval', '%d..%d' % (max(0, a0 - flank), b0 + flank),
               '--context', '100', '--gaf-base', gaf_db, '--gaf-output', q2 + '.gaf',
               '--alignments', 'overlapping'], q2 + '.gfa')
    qs.append(q2)
    q3 = os.path.join(wd, 'q_nodes')
    gbz_query(['--sample', 'CHM13', '--contig', c, '--interval',
               '%d..%d' % (max(0, a0 - flank - 2000), b0 + flank + 2000), '--context', '2000'], q3 + '.gfa')
    nodes = {}
    for q in (q1, q2, q3):
        n, _, _ = parse_gfa(q + '.gfa')
        nodes.update(n)
    alns, dup_diff = {}, 0
    for q in (q1, q2):
        for r in read_gaf(q + '.gaf'):
            k = (r['name'], r['mate'])
            if k in alns:
                if alns[k]['raw_path'] != r['raw_path'] or alns[k]['ps'] != r['ps']:
                    dup_diff += 1
                continue
            alns[k] = r
    missing = sorted({n for r in alns.values() for n, _ in r['steps'] if n not in nodes})
    if missing:
        q4 = os.path.join(wd, 'q_missing')
        args = []
        for n in missing:
            args += ['-n', n]
        gbz_query(args + ['--context', '0'], q4 + '.gfa')
        n, _, _ = parse_gfa(q4 + '.gfa')
        nodes.update(n)
    stats = collections.Counter()
    recs = {}
    for k, r in alns.items():
        stats['alignments'] += 1
        if not r['steps']:
            stats['no_path'] += 1
            continue
        if any(n not in nodes for n, _ in r['steps']):
            stats['missing_node'] += 1
            continue
        pseq = walk_seq(nodes, r['steps'])
        if len(pseq) != r['plen']:
            stats['path_length_mismatch'] += 1
            continue
        s = read_from_cs(pseq[r['ps']:r['pe']], r['cs'])
        if s is None:
            stats['cs_mismatch'] += 1
            continue
        if r['qs'] != 0 or r['qe'] != r['qlen'] or len(s) != r['qlen']:
            stats['clipped_or_short'] += 1
            continue
        q = r['bq'] if r['bq'] and len(r['bq']) == len(s) else None
        if q is None:
            stats['no_quality'] += 1
            q = 'I' * len(s)
        recs[k] = (s, q)
    names = collections.defaultdict(dict)
    for (nm, mate), v in recs.items():
        names[nm][mate] = v
    pairs = sorted(nm for nm, d in names.items() if 1 in d and 2 in d)
    singles = sorted((nm, m) for nm, d in names.items() for m in d if not (1 in d and 2 in d))
    rows = []
    with open(os.path.join(wd, 'reads_1.fq'), 'w') as f1, open(os.path.join(wd, 'reads_2.fq'), 'w') as f2:
        for i, nm in enumerate(pairs):
            ln = 'p%d' % i
            for m, fh in ((1, f1), (2, f2)):
                s, q = names[nm][m]
                fh.write('@%s\n%s\n+\n%s\n' % (ln, s, q))
                r = alns[(nm, m)]
                rows.append((ln, m, nm, len(s), r['mapq'], r['raw_path'], r['ps'], r['pe'], r['plen'], r['cs'],
                             r['tags'].get('pd', '')))
    with open(os.path.join(wd, 'reads_se.fq'), 'w') as fs:
        for i, (nm, m) in enumerate(singles):
            ln = 's%d' % i
            s, q = names[nm][m]
            fs.write('@%s\n%s\n+\n%s\n' % (ln, s, q))
            r = alns[(nm, m)]
            rows.append((ln, m, nm, len(s), r['mapq'], r['raw_path'], r['ps'], r['pe'], r['plen'], r['cs'],
                         r['tags'].get('pd', '')))
    with open(os.path.join(wd, 'reads.tsv'), 'w') as f:
        f.write('local\tmate\tname\tlen\torig_mapq\torig_path\torig_ps\torig_pe\torig_plen\torig_cs\torig_pd\n')
        for row in rows:
            f.write('\t'.join(map(str, row)) + '\n')
    # CHM13 offsets of the nodes (for the fragment-length estimate)
    chm = {}
    with open(q2 + '.gfa') as f:
        for line in f:
            if line.startswith('W\tCHM13'):
                x = line.rstrip('\n').split('\t')
                off = int(x[4])
                for o, n in STEP_RE.findall(x[6]):
                    chm[n] = (off, off + len(nodes[n]), o)
                    off += len(nodes[n])
    frags = []
    for nm in pairs:
        iv = []
        ok = True
        for m in (1, 2):
            r = alns[(nm, m)]
            if not all(n in chm for n, _ in r['steps']):
                ok = False
                break
            lo = min(chm[n][0] for n, _ in r['steps'])
            hi = max(chm[n][1] for n, _ in r['steps'])
            if hi - lo != r['plen']:          # not a colinear CHM13 walk
                ok = False
                break
            fwd = all((o == '+') == (chm[n][2] == '>') for n, o in r['steps'])
            rev = all((o == '-') == (chm[n][2] == '>') for n, o in r['steps'])
            if fwd:
                iv.append((lo + r['ps'], lo + r['pe']))
            elif rev:
                iv.append((hi - r['pe'], hi - r['ps']))
            else:
                ok = False
                break
        if ok:
            frags.append(max(iv[0][1], iv[1][1]) - min(iv[0][0], iv[1][0]))
    info = {'region_id': reg['region_id'], 'flank': flank, 'alignments': stats['alignments'],
            'reads_rebuilt': len(recs), 'pairs': len(pairs), 'singles': len(singles),
            'dropped': {k: v for k, v in stats.items() if k != 'alignments'},
            'duplicate_alignments_differing': dup_diff, 'missing_nodes_fetched': len(missing),
            'fragment_lengths_chm13_colinear': frags, 'seconds': round(time.time() - t0, 1),
            'queries': {'between': '%s:%s' % (L, R),
                        'window_0based': [max(0, a0 - flank), b0 + flank]}}
    with open(os.path.join(wd, 'reads.json'), 'w') as f:
        json.dump(info, f)
    shutil.copy(q1 + '.gaf', os.path.join(wd, 'reads.orig.between.gaf'))
    for q in (q1, q2):
        os.replace(q + '.gaf', q + '.gaf.done')
    return info


def fraglen(ids, work):
    fr = []
    for rid in ids:
        p = os.path.join(work, rid, 'reads.json')
        if os.path.exists(p):
            with open(p) as f:
                fr += json.load(f)['fragment_lengths_chm13_colinear']
    fr = [x for x in fr if 100 <= x <= 2000]
    fr.sort()
    med = statistics.median(fr)
    q1, q3 = fr[len(fr) // 4], fr[3 * len(fr) // 4]
    mad = statistics.median(abs(x - med) for x in fr)
    # trimmed mean/sd within 3 MAD-sigmas
    sig = 1.4826 * mad
    kept = [x for x in fr if abs(x - med) <= 4 * sig]
    d = {'n_pairs': len(fr), 'median': med, 'iqr': [q1, q3], 'mad_sigma': sig,
         'mean': statistics.mean(kept), 'stdev': statistics.pstdev(kept), 'n_kept': len(kept),
         'source': 'pairs with both mates on colinear CHM13 walks in the regions %s' % ','.join(ids)}
    with open(os.path.join(work, 'fraglen.json'), 'w') as f:
        json.dump(d, f, indent=1)
    return d


# ------------------------------------------------------------------ 2. local graphs

def source_gfa(reg, method):
    rid = reg['region_id']
    if method == 'mc':
        return os.path.join(reg['dir'], 'mc.gfa')
    if method == 'truth':
        return None
    return os.path.join(CANDIDATES, method, rid + '.gfa')


def truth_gfa(reg, wd):
    """Graph of HG002's two haplotypes: mafft G-INS-i of truth.fa, induced with msa_graph."""
    out = os.path.join(wd, 'truth.src.gfa')
    if os.path.exists(out):
        return out
    import realign
    import msa_graph
    tfa = os.path.join(reg['dir'], 'truth.fa')
    msa = os.path.join(wd, 'truth.msa.fa')
    os.makedirs(os.path.join(wd, 'truth_mafft'), exist_ok=True)
    # G-INS-i's pairwise DP is quadratic in memory: above 15 kb use FFT-NS-2 (two sequences)
    method = 'mafft_ginsi' if max(len(x) for _, x in read_fasta(tfa)) <= 15000 else 'mafft_fftns2'
    info = realign.align_fasta(method, tfa, msa, threads=2, timeout=3600, mem_mb=10000,
                               workdir=os.path.join(wd, 'truth_mafft'))
    if info['status'] != 'ok':
        raise RuntimeError('truth alignment failed: %s' % info)
    msa_graph.msa_to_gfa(msa, tfa, out + '.tmp', engine='native')
    os.replace(out + '.tmp', out)
    return out


def build_local(src, want, flanks, out_gfa, prefix='h'):
    """Copy of `src` with a CHM13 flank node on each side joined to every path. `want` is the
    [(name, seq)] every path must spell (between the flanks). Paths are renamed prefix+index
    (vg reads '#' as PanSN) and written as haplotype W lines (sample prefix+index, haplotype 0,
    contig 'region'); returns {local W name 'h00#0#region': original name}."""
    nodes, edges, paths = parse_gfa(src)
    byname = {n: s for n, s in paths}
    ids = [int(n) for n in nodes]
    lf, rf = str(max(ids) + 1), str(max(ids) + 2)
    left, right = flanks
    names = {}
    lines_p = []
    new_edges = set((a, ao, b, bo) for a, ao, b, bo in edges)
    for i, (nm, seq) in enumerate(want):
        if nm not in byname:
            raise RuntimeError('%s: no path named %s' % (src, nm))
        st = byname[nm]
        if walk_seq(nodes, st) != seq:
            raise RuntimeError('%s: path %s does not spell its sequence' % (src, nm))
        first, last = st[0], st[-1]
        new_edges.add((lf, '+', first[0], first[1]))
        new_edges.add((last[0], last[1], rf, '+'))
        # W lines, so vg autoindex builds the GBWT from these haplotypes (with only P lines it
        # makes a greedy path cover instead, and giraffe's extension is not haplotype-aware)
        ln = '%s%02d' % (prefix, i)
        names[ln + '#0#region'] = nm
        walk = '>' + lf + ''.join(('>' if o == '+' else '<') + n for n, o in st) + '>' + rf
        lines_p.append('W\t%s\t0\tregion\t0\t%d\t%s\n' % (ln, len(left) + len(seq) + len(right), walk))
    with open(out_gfa + '.tmp', 'w') as f:
        f.write('H\tVN:Z:1.0\n')
        f.write('S\t%s\t%s\n' % (lf, left))
        for n, s in nodes.items():
            f.write('S\t%s\t%s\n' % (n, s))
        f.write('S\t%s\t%s\n' % (rf, right))
        for a, ao, b, bo in sorted(new_edges):
            f.write('L\t%s\t%s\t%s\t%s\t0M\n' % (a, ao, b, bo))
        f.writelines(lines_p)
    # re-check the written file
    n2, _, p2 = parse_gfa(out_gfa + '.tmp')
    want_d = dict(want)
    for ln, st in p2:
        if walk_seq(n2, st) != left + want_d[names[ln]] + right:
            raise RuntimeError('written path %s does not spell flank+sequence+flank' % ln)
    os.replace(out_gfa + '.tmp', out_gfa)
    return {'names': names, 'left_flank_node': lf, 'right_flank_node': rf}


# ------------------------------------------------------------------ 3. index + map

def index_files(prefix):
    fs = {'gbz': prefix + '.giraffe.gbz', 'dist': prefix + '.dist'}
    for mn in (prefix + '.shortread.withzip.min', prefix + '.min'):
        if os.path.exists(mn):
            fs['min'] = mn
            break
    zp = prefix + '.shortread.zipcodes'
    if os.path.exists(zp):
        fs['zip'] = zp
    return fs


def map_graph(reg, method, wd, frag, threads=2, force=False):
    gd = os.path.join(wd, method)
    done = os.path.join(gd, 'map.json')
    if os.path.exists(done) and not force:
        with open(done) as f:
            return json.load(f)
    os.makedirs(gd, exist_ok=True)
    t0 = time.time()
    flanks = ref_flanks(reg, FLANK)
    if method == 'truth':
        src = truth_gfa(reg, wd)
        want = read_fasta(os.path.join(reg['dir'], 'truth.fa'))
        pre = 't'
    else:
        src = source_gfa(reg, method)
        if not os.path.exists(src):
            return {'method': method, 'status': 'no_graph'}
        want = read_fasta(os.path.join(reg['dir'], 'hap32.fa'))
        pre = 'h'
    local = os.path.join(gd, 'local.gfa')
    meta = build_local(src, want, flanks, local, prefix=pre)
    prefix = os.path.join(gd, 'idx')
    tmp = os.path.join(gd, 'tmp')
    os.makedirs(tmp, exist_ok=True)
    t1 = time.time()
    for x in os.listdir(gd):
        if x.startswith('idx.'):
            os.remove(os.path.join(gd, x))
    run([config.VG, 'autoindex', '-n', '-w', 'sr-giraffe', '-g', local, '-p', prefix, '-t', threads,
         '-T', tmp, '-M', '4G'], stdout=os.path.join(gd, 'autoindex.out'), stderr=os.path.join(gd, 'autoindex.log'))
    fs = index_files(prefix)
    t2 = time.time()
    base = [config.VG, 'giraffe', '-Z', fs['gbz'], '-d', fs['dist'], '-m', fs['min']]
    if 'zip' in fs:
        base += ['-z', fs['zip']]
    base += ['-o', 'gaf', '--named-coordinates', '-t', threads]
    out = os.path.join(gd, 'reads.gaf')
    r1, r2, se = (os.path.join(wd, x) for x in ('reads_1.fq', 'reads_2.fq', 'reads_se.fq'))
    with open(out + '.tmp', 'w') as fo:
        if os.path.getsize(r1):
            run(base + ['-f', r1, '-f', r2, '--fragment-mean', '%.1f' % frag['mean'],
                        '--fragment-stdev', '%.1f' % frag['stdev']], stdout=fo,
                stderr=os.path.join(gd, 'giraffe.paired.log'))
        if os.path.getsize(se):
            run(base + ['-f', se], stdout=fo, stderr=os.path.join(gd, 'giraffe.single.log'))
    os.replace(out + '.tmp', out)
    t3 = time.time()
    nodes, _, _ = parse_gfa(local)
    info = {'method': method, 'status': 'ok', 'source_gfa': os.path.relpath(src, REPO) if src.startswith(REPO) else src,
            'local_gfa_nodes': len(nodes), 'local_gfa_bp': sum(len(s) for s in nodes.values()),
            'names': meta['names'], 'left_flank_node': meta['left_flank_node'],
            'right_flank_node': meta['right_flank_node'], 'fragment': frag,
            'seconds': {'build': round(t1 - t0, 1), 'index': round(t2 - t1, 1), 'map': round(t3 - t2, 1)}}
    with open(done, 'w') as f:
        json.dump(info, f)
    for x in os.listdir(tmp):
        p = os.path.join(tmp, x)
        if os.path.isdir(p):
            shutil.rmtree(p, ignore_errors=True)
        else:
            os.remove(p)
    return info


# ------------------------------------------------------------------ 4. analysis helpers

def flank_probe_end(pattern, text, expect):
    """Semi-global: `pattern` aligned end to end, free start and end in `text` (Myers' bit-vector
    approximate matching, unit costs). Returns the text position just after the best match (ties:
    the one closest to `expect`) and its cost."""
    m = len(pattern)
    if m == 0:
        return 0, 0
    peq = {}
    for i, ch in enumerate(pattern):
        peq[ch] = peq.get(ch, 0) | (1 << i)
    mask = (1 << m) - 1
    high = 1 << (m - 1)
    pv, mv, score = mask, 0, m
    best, where = m, [0]
    for j, ch in enumerate(text, 1):
        eq = peq.get(ch, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (~(xh | pv) & mask)
        mh = pv & xh
        if ph & high:
            score += 1
        elif mh & high:
            score -= 1
        ph = (ph << 1) & mask          # row 0 is all zeros: free start in the text
        mh = (mh << 1) & mask
        pv = mh | (~(xv | ph) & mask)
        mv = ph & xv
        if score < best:
            best, where = score, [j]
        elif score == best:
            where.append(j)
    j = min(where, key=lambda x: abs(x - expect))
    return j, best


def core_interval(seq, left_probe, right_probe, exp_left, exp_right_from_end):
    """[start, end) of the core in `seq`: after the match of the CHM13 probe left of the core,
    before the match of the probe right of it."""
    s, c1 = flank_probe_end(left_probe, seq, exp_left)
    j, c2 = flank_probe_end(right_probe[::-1], seq[::-1], exp_right_from_end)
    e = len(seq) - j
    if e < s:
        e = s
    return s, e, c1, c2


def region_cores(reg, seqs):
    """{name: (start, end, cost_left, cost_right)} for the anchor-to-anchor sequences."""
    chm = dict(read_fasta(os.path.join(reg['dir'], 'hap32.fa')))
    ref = next(s for n, s in chm.items() if n.startswith('CHM13'))
    cs0 = reg['core_start'] - reg['span_start']
    ce0 = reg['core_end'] - reg['span_start'] + 1
    lp = ref[max(0, cs0 - FLANK_PROBE):cs0]
    rp = ref[ce0:ce0 + FLANK_PROBE]
    out = {}
    for nm, s in seqs:
        out[nm] = core_interval(s, lp, rp, cs0, len(ref) - ce0)
    return out


def cs_stats(cs, bq):
    """High-quality edit count and details from a cs string and the read's qualities (read
    orientation = the alignment's path orientation for giraffe GAF)."""
    ops = CS_RE.findall(cs)
    q = [ord(c) - 33 for c in bq] if bq else None
    qpos = 0
    hq = 0
    subs = subs_hq = ins_ev = del_ev = 0
    clip = [0, 0]
    for k, o in enumerate(ops):
        c = o[0]
        if c == ':':
            qpos += int(o[1:])
        elif c == '=':
            qpos += len(o) - 1
        elif c == '*':
            subs += 1
            if q is None or q[qpos] >= HQ:
                subs_hq += 1
                hq += 1
            qpos += 1
        elif c == '+':
            n = len(o) - 1
            terminal = k == 0 or k == len(ops) - 1
            if terminal:
                nhq = sum(1 for x in (q[qpos:qpos + n] if q else [HQ] * n) if x >= HQ)
                clip[0 if k == 0 else 1] = n
                if nhq >= 5:
                    hq += 1
            else:
                ins_ev += 1
                hq += 1
            qpos += n
        elif c == '-':
            del_ev += 1
            hq += 1
    return {'hq_edits': hq, 'subs': subs, 'subs_hq': subs_hq, 'ins': ins_ev, 'dels': del_ev,
            'clip5': clip[0], 'clip3': clip[1]}


class PathIndex:
    """Positions of every node occurrence on every path; sub-walk queries."""

    def __init__(self, nodes, paths):
        self.nodes = nodes
        self.paths = paths                       # [(name, steps)]
        self.occ = collections.defaultdict(list)  # node -> [(path idx, step idx)]
        self.offs = []
        for pi, (nm, st) in enumerate(paths):
            off, o = [], 0
            for si, (n, ori) in enumerate(st):
                self.occ[n].append((pi, si))
                off.append(o)
                o += len(nodes[n])
            off.append(o)
            self.offs.append(off)

    def placements(self, steps, ps, pe):
        """[(path idx, fwd start, fwd end, strand)] where the walk is a contiguous sub-walk.
        strand '+' = the walk runs along the path; the interval is on the path's forward strand."""
        out = []
        if not steps:
            return out
        k = len(steps)
        wl = sum(len(self.nodes[n]) for n, _ in steps)
        n0, o0 = steps[0]
        for pi, si in self.occ.get(n0, ()):
            st = self.paths[pi][1]
            if st[si][1] == o0 and st[si:si + k] == steps:
                a = self.offs[pi][si] + ps
                out.append((pi, a, a + (pe - ps), '+'))
        rs = flip_steps(steps)
        n1, o1 = rs[0]
        for pi, si in self.occ.get(n1, ()):
            st = self.paths[pi][1]
            if st[si][1] == o1 and st[si:si + k] == rs:
                base = self.offs[pi][si]
                out.append((pi, base + (wl - pe), base + (wl - ps), '-'))
        return out


def node_core_fraction(nodes, paths, cores):
    """{node: prefix sums of the per-base core fraction}: the fraction of the paths through a node
    that have that base inside their core. `cores` = {path name: (s, e)} in path coordinates."""
    diff = {n: [0] * (len(s) + 1) for n, s in nodes.items()}
    tw = collections.Counter()
    for nm, st in paths:
        s0, e0 = cores[nm]
        o = 0
        for n, ori in st:
            L = len(nodes[n])
            tw[n] += 1
            a, b = max(s0, o), min(e0, o + L)
            if a < b:
                if ori == '+':
                    x, y = a - o, b - o
                else:
                    x, y = L - (b - o), L - (a - o)
                diff[n][x] += 1
                diff[n][y] -= 1
            o += L
    pref = {}
    for n, s in nodes.items():
        t = tw[n]
        acc, run_, out = 0.0, 0, [0.0]
        for i in range(len(s)):
            run_ += diff[n][i]
            acc += (run_ / t) if t else 0.0
            out.append(acc)
        pref[n] = out
    return pref


def aligned_core_bases(steps, ps, pe, nodes, pref):
    """(core, total) graph bases covered by [ps, pe) of the walk, core weighted by fraction."""
    core = tot = 0.0
    o = 0
    for n, ori in steps:
        L = len(nodes[n])
        a, b = max(ps, o), min(pe, o + L)
        if a < b:
            tot += b - a
            if n in pref:
                p = pref[n]
                if ori == '+':
                    core += p[b - o] - p[a - o]
                else:
                    core += p[L - (a - o)] - p[L - (b - o)]
        o += L
    return core, tot


def mapq_bins(mq):
    b = collections.Counter()
    for q in mq:
        b['0' if q == 0 else '1-4' if q < 5 else '5-29' if q < 30 else '30-59' if q < 60 else '60'] += 1
    n = len(mq) or 1
    return {k: round(b[k] / n, 4) for k in ('0', '1-4', '5-29', '30-59', '60')}


def frac(a, b):
    return round(a / b, 4) if b else None


def load_reads_tsv(wd):
    rows = {}
    with open(os.path.join(wd, 'reads.tsv')) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for line in f:
            x = dict(zip(hdr, line.rstrip('\n').split('\t')))
            rows[(x['local'], int(x['mate']) if x['local'].startswith('p') else 0)] = x
    return rows


def local_alignments(gaf):
    """{(local name, mate)}: alignment. Giraffe keeps the FASTQ name; mates carry fn/fp."""
    out = {}
    for r in read_gaf(gaf):
        mate = r['mate'] if r['name'].startswith('p') else 0
        out[(r['name'], mate)] = r
    return out


# ------------------------------------------------------------------ 5. analyse one region

def analyse(reg, wd, methods):
    rid = reg['region_id']
    reads = load_reads_tsv(wd)
    keys = sorted(reads)
    # reads the genome-wide mapping placed on the anchor-to-anchor subgraph (the --between query);
    # the rest come from the flank window and include edge reads that barely touch the local graph
    between = set()
    bg = os.path.join(wd, 'reads.orig.between.gaf')
    if os.path.exists(bg):
        between = {(r['name'], r['mate']) for r in read_gaf(bg)}
    span_keys = [k for k in keys if (reads[k]['name'], int(reads[k]['mate'])) in between]
    # sequences and qualities for cs_stats (in read orientation)
    fq = {}
    for fn, mate in (('reads_1.fq', 1), ('reads_2.fq', 2), ('reads_se.fq', 0)):
        with open(os.path.join(wd, fn)) as f:
            while True:
                h = f.readline()
                if not h:
                    break
                s = f.readline().strip()
                f.readline()
                q = f.readline().strip()
                fq[(h[1:].strip(), mate)] = (s, q)
    hap = read_fasta(os.path.join(reg['dir'], 'hap32.fa'))
    truth = read_fasta(os.path.join(reg['dir'], 'truth.fa'))
    hapd = dict(hap)
    cores_h = region_cores(reg, hap)
    cores_t = region_cores(reg, truth)
    tcore_len = statistics.mean(e - s for s, e, _, _ in (cores_t[n] for n, _ in truth))
    tlen = statistics.mean(len(s) for _, s in truth)
    tflank_len = tlen - tcore_len + 2 * FLANK
    # distinct haplotype sequences and the truth-closest ones
    import region as regionlib
    distinct = {}
    for nm, s in hap:
        distinct.setdefault(s, []).append(nm)
    dseqs = list(distinct)
    hap_group = {nm: gi for gi, s in enumerate(dseqs) for nm in distinct[s]}
    closest = {}
    for tn, ts in truth:
        ds = [(regionlib.edit_distance(ts, s), gi) for gi, s in enumerate(dseqs)]
        d, gi = min(ds)
        closest[tn] = {'group': gi, 'edit': d, 'paths': distinct[dseqs[gi]]}
    # truth haplotypes with an exact (or near: edit <= NEAR_EDIT) hap32 copy; for a near copy the
    # truth coordinates are carried to the copy's through a global alignment of the two
    exact = {}
    for tn, ts in truth:
        v = closest[tn]
        if v['edit'] > NEAR_EDIT:
            continue
        if v['edit'] == 0:
            exact[tn] = (v['paths'], None, len(ts), len(ts), None)
            continue
        pseq = dseqs[v['group']]
        _, cig = regionlib.align_cigar(pseq, ts)
        if cig is None:
            continue
        m, i, j = [0] * (len(ts) + 1), 0, 0
        e = [0] * (len(ts) + 1)            # 1 where the truth differs from the copy
        for n, op in re.findall(r'(\d+)([=XID])', cig):
            for _ in range(int(n)):
                if op in '=X':
                    m[j] = i
                    e[j] = 1 if op == 'X' else 0
                    i += 1
                    j += 1
                elif op == 'I':
                    m[j] = i
                    e[j] = 1
                    j += 1
                else:
                    e[j] = 1
                    i += 1
        m[len(ts)] = i
        epref = [0]
        for x in e:
            epref.append(epref[-1] + x)
        exact[tn] = (v['paths'], m, len(ts), len(pseq), epref)

    def expected_pos(tn, ta):
        """truth-graph path offset -> offset on the hap32 copy's local path"""
        _, m, lt, lp, _ = exact[tn]
        if m is None or ta < FLANK:
            return ta
        if ta > FLANK + lt:
            return ta - lt + lp
        return FLANK + m[ta - FLANK]

    def copy_identical(tn, ta, tb):
        """True when the truth interval [ta, tb) (local truth-path offsets) is spelled identically
        by the hap32 copy, so the read's right place on the copy is defined."""
        _, m, lt, lp, epref = exact[tn]
        if m is None:
            return True
        lo, hi = max(0, ta - FLANK), min(lt, tb - FLANK)
        if hi < lo:
            return True
        return epref[min(lt + 1, hi + 1)] - epref[lo] == 0

    res = {'region_id': rid, 'stratum': reg['stratum'], 'n_reads': len(keys),
           'n_hap32_paths': len(hap), 'n_distinct_haplotypes': len(dseqs),
           'truth_core_len_mean': round(tcore_len, 1), 'truth_flank_len_mean': round(tflank_len, 1),
           'truth_closest': {tn: {'edit': v['edit'], 'paths': v['paths']} for tn, v in closest.items()},
           'cores_hap32': {n: cores_h[n][:2] for n in cores_h}, 'cores_truth': {n: cores_t[n][:2] for n in cores_t},
           'core_probe_costs_max': max(max(v[2], v[3]) for v in list(cores_h.values()) + list(cores_t.values())),
           'graphs': {}}

    # ---- truth graph: reference placements and truth-consistency
    tinfo = json.load(open(os.path.join(wd, 'truth', 'map.json')))
    tn_nodes, _, tn_paths = parse_gfa(os.path.join(wd, 'truth', 'local.gfa'))
    t_orig = [(tinfo['names'][ln], st) for ln, st in tn_paths]
    tpi = PathIndex(tn_nodes, t_orig)
    tcores = {nm: (FLANK + cores_t[nm][0], FLANK + cores_t[nm][1]) for nm, _ in truth}
    tpref = node_core_fraction(tn_nodes, t_orig, tcores)
    talns = local_alignments(os.path.join(wd, 'truth', 'reads.gaf'))
    tread = {}
    for k in keys:
        r = talns.get(k)
        s, q = fq[k]
        d = {'mapped': bool(r and r['steps']), 'mapq': r['mapq'] if r else 0}
        if d['mapped']:
            cst = cs_stats(r['cs'], q)
            d.update(cst)
            d['AS'] = r['AS']
            core, tot = aligned_core_bases(r['steps'], r['ps'], r['pe'], tn_nodes, tpref)
            d['core_bp'] = core
            d['placements'] = [(t_orig[pi][0], a, b, sd) for pi, a, b, sd in tpi.placements(r['steps'], r['ps'], r['pe'])]
            d['consistent'] = cst['hq_edits'] <= 1
        else:
            d['consistent'] = False
            d['core_bp'] = 0.0
        tread[k] = d
    truth_core = {k for k in keys if tread[k]['mapped'] and tread[k]['core_bp'] >= CORE_READ_BP}
    consistent = {k for k in keys if tread[k]['consistent']}
    evaluable = set()
    for k in consistent:
        d = tread[k]
        if d['mapq'] >= TRUTH_MAPQ_UNIQUE and d['placements'] and all(
                p[0] in exact and copy_identical(p[0], p[1], p[2]) for p in d['placements']):
            evaluable.add(k)
    res['truth_reads'] = {'truth_core_reads': len(truth_core), 'truth_consistent': len(consistent),
                          'truth_consistent_frac': frac(len(consistent), len(keys)),
                          'truth_consistent_core': len(consistent & truth_core),
                          'truth_consistent_core_frac': frac(len(consistent & truth_core), len(truth_core)),
                          'placement_evaluable': len(evaluable),
                          'placement_evaluable_core': len(evaluable & truth_core),
                          'truth_copies_used': {tn: {'paths': v[0], 'edit': closest[tn]['edit']}
                                                for tn, v in exact.items()},
                          'hq_edit_hist_all': dict(collections.Counter(min(tread[k].get('hq_edits', 99), 5) for k in keys if tread[k]['mapped'])),
                          'hq_edit_hist_core': dict(collections.Counter(min(tread[k].get('hq_edits', 99), 5) for k in truth_core))}

    # ---- every graph (and the original genome-wide mapping)
    for method in methods + ['orig']:
        if method == 'orig':
            gd = os.path.join(wd, 'mc')
            if not os.path.exists(os.path.join(gd, 'map.json')):
                continue
        else:
            gd = os.path.join(wd, method)
            mp = os.path.join(gd, 'map.json')
            if not os.path.exists(mp):
                src = source_gfa(reg, method)
                res['graphs'][method] = {'status': 'no_graph' if src and not os.path.exists(src) else 'not_run'}
                continue
        info = json.load(open(os.path.join(gd, 'map.json')))
        if info.get('status') != 'ok':
            res['graphs'][method] = {'status': info.get('status')}
            continue
        nodes, _, lpaths = parse_gfa(os.path.join(gd, 'local.gfa'))
        names = info['names']
        opaths = [(names[ln], st) for ln, st in lpaths]
        if method == 'truth':
            pcores = tcores
        else:
            pcores = {nm: (FLANK + cores_h[nm][0], FLANK + cores_h[nm][1]) for nm, _ in opaths}
        pref = node_core_fraction(nodes, opaths, pcores)
        pi_ = PathIndex(nodes, opaths)
        flank_nodes = {info['left_flank_node'], info['right_flank_node']}
        inner = set(nodes) - flank_nodes
        if method == 'orig':
            alns = {}
            for k in keys:
                x = reads[k]
                steps = [] if x['orig_path'] == '*' else [(n, '+' if o == '>' else '-') for o, n in STEP_RE.findall(x['orig_path'])]
                alns[k] = {'steps': steps, 'ps': int(x['orig_ps']), 'pe': int(x['orig_pe']), 'mapq': int(x['orig_mapq']),
                           'cs': x['orig_cs'], 'AS': None, 'plen': int(x['orig_plen'])}
        else:
            alns = local_alignments(os.path.join(gd, 'reads.gaf'))
        per = {}
        core_bp = flank_bp = core_bp5 = flank_bp5 = core_eff = flank_eff = 0.0
        for k in keys:
            r = alns.get(k)
            d = {'mapped': bool(r and r['steps']), 'mapq': (r['mapq'] if r and r['steps'] else None)}
            if d['mapped']:
                steps, ps, pe = r['steps'], r['ps'], r['pe']
                if method == 'orig':
                    # trim steps outside the anchor-to-anchor graph; their bases count as flank.
                    # (The flank node ids of the local graph may be real hap32 ids outside the
                    # span, so membership is tested against the inner nodes only.)
                    lo, hi = 0, len(steps)
                    while lo < hi and steps[lo][0] not in inner:
                        lo += 1
                    while hi > lo and steps[hi - 1][0] not in inner:
                        hi -= 1
                    if lo >= hi or any(steps[i][0] not in inner for i in range(lo, hi)):
                        # entirely outside (or re-leaves): all flank
                        d.update({'core_bp': 0.0, 'aln_bp': pe - ps, 'placements': [], 'outside': True,
                                  'trimmed': True})
                        flank_bp += pe - ps
                        if d['mapq'] >= 5:
                            flank_bp5 += pe - ps
                        flank_eff += (pe - ps) * (1 - 10 ** (-d['mapq'] / 10))
                        per[k] = d
                        continue
                    pre = 0
                    for i in range(lo):
                        pre += _orig_len(steps[i][0], r)
                    inner_len = sum(len(nodes[steps[i][0]]) for i in range(lo, hi))
                    ps2, pe2 = max(0, ps - pre), min(inner_len, pe - pre)
                    outside_bp = (pe - ps) - max(0, pe2 - ps2)
                    d['trimmed'] = lo > 0 or hi < len(steps)
                    steps, ps, pe = steps[lo:hi], ps2, pe2
                else:
                    outside_bp = 0
                core, tot = aligned_core_bases(steps, ps, pe, nodes, pref)
                fl = tot - core + outside_bp
                d['core_bp'] = core
                d['aln_bp'] = tot + outside_bp
                core_bp += core
                flank_bp += fl
                if d['mapq'] >= 5:
                    core_bp5 += core
                    flank_bp5 += fl
                w = 1 - 10 ** (-d['mapq'] / 10)
                core_eff += core * w
                flank_eff += fl * w
                d['span'] = any(n not in flank_nodes for n, _ in steps)
                d['placements'] = pi_.placements(steps, ps, pe) if steps else []
                if method != 'orig':
                    s, q = fq[k]
                    d.update(cs_stats(r['cs'], q))
                    d['AS'] = r['AS']
            per[k] = d
        mapped = [k for k in keys if per[k]['mapped']]
        sk_m = [k for k in span_keys if per[k]['mapped']]
        gcore = [k for k in mapped if per[k]['core_bp'] >= CORE_READ_BP]
        span = [k for k in mapped if per[k].get('span')]
        g = {'status': 'ok', 'nodes': len(nodes) - 2,
             'reads': len(keys), 'mapped_frac': frac(len(mapped), len(keys)),
             'mapq_lt5_frac': frac(sum(1 for k in mapped if per[k]['mapq'] < 5), len(mapped)),
             'mapq0_frac': frac(sum(1 for k in mapped if per[k]['mapq'] == 0), len(mapped)),
             'mapq_mean': round(statistics.mean(per[k]['mapq'] for k in mapped), 2) if mapped else None,
             'mapq_bins': mapq_bins([per[k]['mapq'] for k in mapped]),
             'span_input_reads': len(span_keys),
             'span_input_mapped_frac': frac(len(sk_m), len(span_keys)),
             'span_input_mapq_lt5_frac': frac(sum(1 for k in sk_m if per[k]['mapq'] < 5), len(sk_m)),
             'graph_core_reads': len(gcore),
             'core_mapq_lt5_frac': frac(sum(1 for k in gcore if per[k]['mapq'] < 5), len(gcore)),
             'core_mapq_bins': mapq_bins([per[k]['mapq'] for k in gcore]),
             'core_depth': round(core_bp / tcore_len, 3) if tcore_len else None,
             'flank_depth': round(flank_bp / tflank_len, 3),
             'core_depth_mq5': round(core_bp5 / tcore_len, 3) if tcore_len else None,
             'flank_depth_mq5': round(flank_bp5 / tflank_len, 3),
             'core_depth_eff': round(core_eff / tcore_len, 3) if tcore_len else None,
             'flank_depth_eff': round(flank_eff / tflank_len, 3)}
        if method == 'orig':
            # genome-wide alignments that leave the window count their outside bases as flank, so the
            # orig flank depth is inflated; only the core depths compare with the local graphs
            for x in ('flank_depth', 'flank_depth_mq5', 'flank_depth_eff'):
                g[x] = None
        g['core_flank_ratio'] = round(g['core_depth'] / g['flank_depth'], 3) if g['flank_depth'] else None
        g['core_flank_ratio_mq5'] = round(g['core_depth_mq5'] / g['flank_depth_mq5'], 3) if g['flank_depth_mq5'] else None
        g['core_flank_ratio_eff'] = round(g['core_depth_eff'] / g['flank_depth_eff'], 3) if g['flank_depth_eff'] else None
        # truth-core reads (graph-independent set)
        tc = sorted(truth_core)
        tcm = [k for k in tc if per[k]['mapped']]
        g['truth_core_mapped_frac'] = frac(len(tcm), len(tc))
        g['truth_core_mapq_lt5_frac'] = frac(sum(1 for k in tcm if per[k]['mapq'] < 5), len(tcm))
        if method != 'truth':
            # haplotype compatibility on this graph's hap32 paths (by sequence group)
            def groups(k):
                return {hap_group[opaths[p[0]][0]] for p in per[k]['placements']}
            for label, sel in (('span', span), ('core', gcore), ('truth_core', [k for k in tc if per[k]['mapped'] and per[k].get('span')])):
                comp = [k for k in sel if per[k]['placements']]
                g['%s_reads' % label] = len(sel)
                g['%s_compatible_frac' % label] = frac(len(comp), len(sel))
                if label == 'span':
                    continue
                gsets = {k: groups(k) for k in comp}
                used = collections.Counter(x for k in comp for x in gsets[k])
                g['%s_haps_touched' % label] = len(used)
                g['%s_haps_per_read_mean' % label] = round(statistics.mean(len(v) for v in gsets.values()), 2) if comp else None
                # greedy cover of the compatible reads
                left = set(comp)
                ncov = 0
                target = 0.95 * len(comp)
                covered = 0
                while left and covered < target:
                    best = max(range(len(dseqs)), key=lambda gi: sum(1 for k in left if gi in gsets[k]))
                    got = {k for k in left if best in gsets[k]}
                    if not got:
                        break
                    left -= got
                    covered += len(got)
                    ncov += 1
                g['%s_haps_cover95' % label] = ncov
                # best pair
                bits = [0] * len(dseqs)
                for i, k in enumerate(comp):
                    for x in gsets[k]:
                        bits[x] |= 1 << i
                best, bp = 0, None
                for a in range(len(dseqs)):
                    for b in range(a, len(dseqs)):
                        c = (bits[a] | bits[b]).bit_count()
                        if c > best:
                            best, bp = c, (a, b)
                g['%s_best_pair_frac' % label] = frac(best, len(sel))
                g['%s_best_pair' % label] = sorted({distinct[dseqs[x]][0] for x in bp}) if bp else None
                ta, tb = (closest[tn]['group'] for tn, _ in truth)
                g['%s_truth_pair_frac' % label] = frac((bits[ta] | bits[tb]).bit_count(), len(sel))
                g['%s_best_pair_is_truth_closest' % label] = (bool(bp) and set(bp) == {ta, tb})
        # truth-consistent reads: confidence, representation, placement
        tcs = sorted(consistent)
        tcsm = [k for k in tcs if per[k]['mapped']]
        g['truth_consistent_mapped_frac'] = frac(len(tcsm), len(tcs))
        g['truth_consistent_mapq_ge5_frac'] = frac(sum(1 for k in tcsm if per[k]['mapq'] >= 5), len(tcs))
        tcc = sorted(consistent & truth_core)
        g['truth_consistent_core_mapq_ge5_frac'] = frac(sum(1 for k in tcc if per[k]['mapped'] and per[k]['mapq'] >= 5), len(tcc))
        if method != 'orig':
            g['truth_consistent_graph_as_good_frac'] = frac(
                sum(1 for k in tcsm if per[k]['hq_edits'] <= tread[k]['hq_edits']), len(tcs))
            g['truth_consistent_core_graph_as_good_frac'] = frac(
                sum(1 for k in tcc if per[k]['mapped'] and per[k]['hq_edits'] <= tread[k]['hq_edits']), len(tcc))
        name_to_pi = {nm: i for i, (nm, _) in enumerate(opaths)}
        if method != 'truth' and exact:

            def correct(k):
                d = per[k]
                if not d['mapped'] or not d['placements']:
                    return False
                return any(pj == name_to_pi[pn] and sd == tsd and abs(a - expected_pos(tn, ta_)) <= POS_TOL
                           for tn, ta_, tb_, tsd in tread[k]['placements'] for pn in exact[tn][0]
                           for pj, a, b, sd in d['placements'])
            # orig: reads whose genome-wide alignment leaves the anchor span are not evaluable
            # (their flank lies on real graph nodes, not on the local CHM13 flank)
            ev = sorted(k for k in evaluable if not per[k].get('trimmed'))
            for label, sel in (('', ev), ('core_', [k for k in ev if k in truth_core])):
                ok = sum(1 for k in sel if correct(k))
                conf = [k for k in sel if per[k]['mapped'] and per[k]['mapq'] >= 5]
                conf_ok = sum(1 for k in conf if correct(k))
                g[label + 'placement_evaluable'] = len(sel)
                g[label + 'placed_correctly_frac'] = frac(ok, len(sel))
                g[label + 'placed_confident_frac'] = frac(len(conf), len(sel))
                g[label + 'placed_confident_correct_frac'] = frac(conf_ok, len(sel))
                g[label + 'confident_wrong_frac'] = frac(len(conf) - conf_ok, len(conf))
        if method != 'orig':
            g['seconds'] = info.get('seconds')
            g['source_gfa'] = info.get('source_gfa')
        res['graphs'][method] = g
        # per-read table for inspection
        pdir = os.path.join(wd, method if method != 'orig' else 'mc')
        with open(os.path.join(pdir, 'per_read.%s.tsv' % method), 'w') as f:
            f.write('read\tmate\tmapped\tmapq\tcore_bp\tspan\tn_paths\thaps\thq_edits\ttruth_hq_edits\t'
                    'truth_mapq\ttruth_core\ttruth_consistent\tevaluable\tcorrect\tplacements\ttruth_placements\n')
            for k in keys:
                d = per[k]
                t = tread[k]
                pl = d.get('placements') or []
                ok = ''
                if method != 'truth' and exact and k in evaluable and not d.get('trimmed'):
                    ok = int(correct(k))
                f.write('\t'.join(map(str, [
                    k[0], k[1], int(d['mapped']), '' if d['mapq'] is None else d['mapq'],
                    round(d.get('core_bp', 0.0), 1), int(bool(d.get('span'))), len(pl),
                    ','.join(sorted({str(hap_group.get(opaths[p[0]][0], opaths[p[0]][0])) for p in pl})),
                    d.get('hq_edits', ''), t.get('hq_edits', ''), t['mapq'], int(k in truth_core),
                    int(k in consistent), int(k in evaluable), ok,
                    ';'.join('%s:%d%s' % (opaths[p[0]][0], p[1], p[3]) for p in pl[:4]),
                    ';'.join('%s:%d%s' % (p[0], p[1], p[3]) for p in (t.get('placements') or [])[:2])])) + '\n')
    with open(os.path.join(wd, 'stage2.json'), 'w') as f:
        json.dump(res, f, indent=1)
    return res


_ORIG_NODE_LEN = {}


def _orig_len(n, r):
    return _ORIG_NODE_LEN[n]


def load_orig_node_lengths(wd):
    for q in ('q_between.gfa', 'q_window.gfa', 'q_nodes.gfa', 'q_missing.gfa'):
        p = os.path.join(wd, q)
        if os.path.exists(p):
            n, _, _ = parse_gfa(p)
            for k, v in n.items():
                _ORIG_NODE_LEN[k] = len(v)


# ------------------------------------------------------------------ driver

def region_work(work, rid):
    return os.path.join(work, rid)


def do_region(rid, methods, work, threads, force=False, steps=('fetch', 'map', 'analyse'), frag=None):
    reg = load_region(rid)
    wd = region_work(work, rid)
    out = {'region_id': rid}
    if 'fetch' in steps and (force or not os.path.exists(os.path.join(wd, 'reads.json'))):
        out['fetch'] = fetch(reg, wd)
    if 'map' in steps:
        if frag is None:
            frag = json.load(open(os.path.join(work, 'fraglen.json')))
        ms = list(methods)
        if 'truth' not in ms:
            ms.append('truth')
        for m in ms:
            try:
                info = map_graph(reg, m, wd, frag, threads=threads, force=force)
                log(rid, m, info.get('status'), info.get('seconds', ''))
            except Exception as e:  # noqa: BLE001
                log(rid, m, 'FAILED', e)
                os.makedirs(os.path.join(wd, m), exist_ok=True)
                with open(os.path.join(wd, m, 'map.json'), 'w') as f:
                    json.dump({'method': m, 'status': 'failed', 'error': str(e)[:2000]}, f)
    if 'analyse' in steps:
        load_orig_node_lengths(wd)
        out['analyse'] = analyse(reg, wd, [m for m in methods if m != 'truth'] + ['truth'])
    return out


TABLE_COLS = ['region_id', 'stratum', 'graph', 'status', 'nodes', 'stage0_cost_over_opt', 'reads', 'mapped_frac',
              'span_input_reads', 'span_input_mapped_frac', 'span_input_mapq_lt5_frac', 'mapq_mean',
              'mapq_lt5_frac', 'mapq0_frac', 'graph_core_reads', 'core_mapq_lt5_frac',
              'truth_core_mapq_lt5_frac', 'core_depth', 'flank_depth', 'core_flank_ratio',
              'core_depth_mq5', 'flank_depth_mq5', 'core_flank_ratio_mq5', 'core_flank_ratio_eff',
              'span_compatible_frac', 'core_compatible_frac', 'core_haps_touched', 'core_haps_per_read_mean',
              'core_haps_cover95', 'core_best_pair_frac', 'core_truth_pair_frac', 'core_best_pair_is_truth_closest',
              'truth_core_compatible_frac', 'truth_core_best_pair_frac', 'truth_core_truth_pair_frac',
              'truth_consistent_mapq_ge5_frac', 'truth_consistent_core_mapq_ge5_frac',
              'truth_consistent_graph_as_good_frac', 'truth_consistent_core_graph_as_good_frac',
              'placement_evaluable', 'placed_correctly_frac', 'placed_confident_correct_frac', 'confident_wrong_frac',
              'core_placement_evaluable', 'core_placed_correctly_frac', 'core_placed_confident_correct_frac',
              'core_confident_wrong_frac',
              'n_truth_core_reads', 'n_truth_consistent', 'truth_consistent_core_frac', 'n_distinct_haplotypes',
              'truth_closest_edit']


def table(ids, work, out_tsv, stage0_dirs=()):
    rows = []
    for rid in ids:
        p = os.path.join(work, rid, 'stage2.json')
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        for gname, g in d['graphs'].items():
            row = dict(g)
            row.update({'region_id': rid, 'stratum': d['stratum'], 'graph': gname,
                        'n_truth_core_reads': d['truth_reads']['truth_core_reads'],
                        'n_truth_consistent': d['truth_reads']['truth_consistent'],
                        'truth_consistent_core_frac': d['truth_reads']['truth_consistent_core_frac'],
                        'n_distinct_haplotypes': d['n_distinct_haplotypes'],
                        'truth_closest_edit': '/'.join(str(v['edit']) for v in d['truth_closest'].values()),
                        'stage0_cost_over_opt': stage0_cost(gname, rid, stage0_dirs) if g.get('status') == 'ok' else None})
            rows.append(row)
    with open(out_tsv, 'w') as f:
        f.write('\t'.join(TABLE_COLS) + '\n')
        for r in rows:
            f.write('\t'.join('' if r.get(c) is None else str(r.get(c)) for c in TABLE_COLS) + '\n')
    return rows


SUMMARY_METRICS = [
    ('stage0_cost_over_opt', 'Stage 0: all-pairs cost/opt of the graph (for reference)', 'low'),
    ('truth_core_mapq_lt5_frac', 'MAPQ<5, truth-core reads', 'low'),
    ('truth_consistent_core_mapq_ge5_frac', 'truth-consistent core reads at MAPQ>=5', 'high'),
    ('core_placed_confident_correct_frac', 'uniquely-placeable core reads placed right at MAPQ>=5', 'high'),
    ('core_confident_wrong_frac', 'core reads at MAPQ>=5 placed wrong', 'low'),
    ('core_placed_correctly_frac', 'uniquely-placeable core reads placed right (any MAPQ)', 'high'),
    ('core_compatible_frac', 'core reads on a single haplotype path', 'high'),
    ('core_haps_cover95', 'haplotypes to cover 95% of compatible core reads', 'low'),
    ('core_best_pair_frac', 'core reads explained by the best path pair', 'high'),
    ('core_truth_pair_frac', 'core reads explained by the truth-closest pair', 'high'),
    ('span_input_mapq_lt5_frac', 'MAPQ<5, all reads the genome-wide mapping put on the anchor span', 'low'),
    ('core_flank_ratio', 'core/flank depth, all reads', 'high'),
    ('core_flank_ratio_mq5', 'core/flank depth, MAPQ>=5', 'high'),
    ('core_depth_mq5', 'core depth from MAPQ>=5 reads (x)', 'high'),
    ('flank_depth_mq5', 'flank depth from MAPQ>=5 reads (x)', 'high'),
    ('truth_consistent_core_graph_as_good_frac', 'truth-consistent core reads aligned as well as to truth', 'high'),
    ('span_input_mapped_frac', 'reads mapped, of those the genome-wide mapping put on the anchor span', 'high'),
    ('mapped_frac', 'reads mapped, all fetched reads (includes window-edge reads)', 'high'),
]
SUMMARY_GRAPHS = ['orig', 'mc', 'mafft_linsi', 'unit_aware', 'poa_spoa', 'poa_abpoa', 'unit_aware__all',
                  'mafft_linsi__all', 'poa_spoa__all', 'poa_abpoa__all', 'poa_abpoa_mc__all', 'truth']
GROUPS = [('VNTR hotspots', ('hotspot_vntr',)), ('non-VNTR hotspot', ('hotspot_other',)),
          ('matched VNTR controls', ('control_vntr_matched',)), ('non-repeat SV controls', ('control_nontr_sv',))]


def _fmt(x, nd=3):
    if x is None or x == '':
        return '-'
    if isinstance(x, bool):
        return 'yes' if x else 'no'
    if isinstance(x, float):
        return ('%.' + str(nd) + 'f') % x
    return str(x)


def stage0_cost(method, rid, extra_dirs=()):
    """Stage 0 all-pairs cost/opt of a graph from evaluate.py output (results/<method>/<id>.json,
    or the same layout under one of extra_dirs)."""
    if method in ('orig', 'truth'):
        method = 'mc' if method == 'orig' else None
    if not method:
        return None
    for d in (RESULTS,) + tuple(extra_dirs):
        p = os.path.join(d, method, rid + '.json')
        if os.path.exists(p):
            try:
                with open(p) as f:
                    return json.load(f)['alignment'].get('all_cost_over_opt')
            except (ValueError, KeyError):
                return None
    return None


def summary(ids, work, out_md=None, stage0_dirs=()):
    """Markdown tables: Stage 0 cost/opt for reference, then one table per metric (regions x
    graphs) with per-group medians and paired comparisons with mc."""
    data = {}
    for rid in ids:
        p = os.path.join(work, rid, 'stage2.json')
        if os.path.exists(p):
            data[rid] = json.load(open(p))
    ids = [r for r in ids if r in data]
    graphs = [g for g in SUMMARY_GRAPHS if any(g in data[r]['graphs'] and data[r]['graphs'][g].get('status') == 'ok'
                                               for r in ids)]
    L = []
    for rid in ids:
        for gname in graphs:
            g = data[rid]['graphs'].get(gname)
            if g and g.get('status') == 'ok':
                g['stage0_cost_over_opt'] = stage0_cost(gname, rid, stage0_dirs)
    for key, label, better in SUMMARY_METRICS:
        L.append('\n#### %s (`%s`, %s is better)\n' % (label, key, better))
        L.append('| region | stratum | ' + ' | '.join(graphs) + ' |')
        L.append('|---|---|' + '---|' * len(graphs))
        for rid in ids:
            d = data[rid]
            row = [rid, d['stratum'].replace('control_', 'c_').replace('hotspot_', 'hs_')]
            for gname in graphs:
                g = d['graphs'].get(gname, {})
                row.append(_fmt(g.get(key)) if g.get('status') == 'ok' else '')
            L.append('| ' + ' | '.join(row) + ' |')
        # paired against mc, by group
        L.append('')
        L.append('| group: median (paired diff vs mc; better/worse) | ' + ' | '.join(graphs) + ' |')
        L.append('|---|' + '---|' * len(graphs))
        for gl, strata in GROUPS:
            rs = [r for r in ids if data[r]['stratum'] in strata]
            if not rs:
                continue
            cells = []
            for gname in graphs:
                vals, diffs = [], []
                for r in rs:
                    g = data[r]['graphs'].get(gname, {})
                    m = data[r]['graphs'].get('mc', {})
                    v, v0 = g.get(key), m.get(key)
                    if g.get('status') == 'ok' and isinstance(v, (int, float)) and not isinstance(v, bool):
                        vals.append(v)
                        if isinstance(v0, (int, float)):
                            diffs.append(v - v0)
                if not vals:
                    cells.append('')
                    continue
                med = statistics.median(vals)
                if gname == 'mc' or not diffs:
                    cells.append('%s (n=%d)' % (_fmt(float(med)), len(vals)))
                    continue
                sgn = 1 if better == 'high' else -1
                nb = sum(1 for x in diffs if sgn * x > 1e-9)
                nw = sum(1 for x in diffs if sgn * x < -1e-9)
                cells.append('%s (%+.3f; %d/%d)' % (_fmt(float(med)), statistics.median(diffs), nb, nw))
            L.append('| %s (%d) | ' % (gl, len(rs)) + ' | '.join(cells) + ' |')
    text = '\n'.join(L) + '\n'
    if out_md:
        with open(out_md, 'w') as f:
            f.write(text)
    return text


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['fetch', 'fraglen', 'map', 'analyse', 'run', 'table', 'summary'])
    ap.add_argument('ids', nargs='*')
    ap.add_argument('--methods', default=','.join(DEFAULT_METHODS))
    ap.add_argument('--work', default=WORK)
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--tsv', default=os.path.join(RESULTS, 'stage2_remap.tsv'))
    ap.add_argument('--md', default=None, help='summary: also write the markdown tables here')
    ap.add_argument('--stage0-dir', action='append', default=[],
                    help='summary: extra directory laid out like results/ holding evaluate.py JSONs')
    a = ap.parse_args(argv)
    methods = [m for m in a.methods.split(',') if m]
    os.makedirs(a.work, exist_ok=True)
    ids = a.ids
    if a.cmd == 'fraglen':
        print(json.dumps(fraglen(ids, a.work), indent=1))
        return 0
    if a.cmd == 'summary':
        if not ids:
            ids = sorted(x for x in os.listdir(a.work) if os.path.exists(os.path.join(a.work, x, 'stage2.json')))
        print(summary(ids, a.work, a.md, tuple(a.stage0_dir)))
        return 0
    if a.cmd == 'table':
        if not ids:
            ids = sorted(x for x in os.listdir(a.work) if os.path.exists(os.path.join(a.work, x, 'stage2.json')))
        rows = table(ids, a.work, a.tsv, tuple(a.stage0_dir))
        print('%d rows -> %s' % (len(rows), a.tsv))
        return 0
    steps = {'fetch': ('fetch',), 'map': ('map',), 'analyse': ('analyse',),
             'run': ('fetch', 'map', 'analyse')}[a.cmd]
    if a.cmd == 'run':
        # fetch everything first, then the pooled fragment length, then map and analyse
        for rid in ids:
            do_region(rid, methods, a.work, a.threads, a.force, steps=('fetch',))
        if a.force or not os.path.exists(os.path.join(a.work, 'fraglen.json')):
            fraglen(ids, a.work)
        steps = ('map', 'analyse')
    if a.jobs > 1 and len(ids) > 1:
        with concurrent.futures.ProcessPoolExecutor(a.jobs) as ex:
            futs = {ex.submit(do_region, rid, methods, a.work, a.threads, a.force, steps): rid for rid in ids}
            for fu in concurrent.futures.as_completed(futs):
                rid = futs[fu]
                try:
                    fu.result()
                    log(rid, 'done')
                except Exception as e:  # noqa: BLE001
                    log(rid, 'FAILED', repr(e))
    else:
        for rid in ids:
            do_region(rid, methods, a.work, a.threads, a.force, steps)
            log(rid, 'done')
    return 0


if __name__ == '__main__':
    sys.exit(main())
