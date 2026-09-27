#!/usr/bin/env python3
"""evaluate.py -- Stage 0 (graph only) and Stage 1 (existing reads) metrics of one graph of one region.

    python3 tools/evaluate.py <region_dir> <graph.gfa> [--name METHOD] [--reads]
                              [--out results/<METHOD>/<region_id>.json]
                              [--max-pairs 600] [--affine-budget 5e10] [--no-truvari] [--no-refine]
                              [--threads 4] [--seed 1]
    python3 tools/evaluate.py <region_dir> <panel_graph.gfa> --panel <panel.fa[.gz]> [--max-pairs 300]

<region_dir> is a region package (regions/<region_id>/: region.json, hap32.fa, truth.fa, ...).
<graph.gfa> is the graph to judge: GFA 1.0 with one P (or W) line per hap32.fa sequence, with the
same names, spelling exactly the same sequences. The script first asserts that; a graph that
changes a sequence is invalid (exit status 3, and the JSON says so). The graph must be acyclic on
the handles between the paths' first and last handles (it is asserted; a cyclic graph gets the
metrics that do not need a DAG and NA for the others).

Every metric is defined in tools/METRICS.md. In short:
  size         nodes, edges, bp, nodes per CHM13 kb, node lengths, top-level bubbles
  alignment    graph-implied pairwise cost / optimal pairwise cost over path pairs (unit and
               gap-affine), excess edits per kb, unaligned homology per kb
  inflation    each path expressed against CHM13 through the graph: SV-sized pieces, indel bp
  redundancy   21-mer graph positions beyond the most copies any one path has
  truth        HG002 haplotypes vs any graph path (exact sequence-to-DAG DP, banded by the
               distance to the closest panel path) and vs the closest panel path; the truth written
               the way the graph decomposes it, scored by truvari: raw, after truvari refine on its
               candidate regions (as the investigation), and whole-span phab (-w)
  reads        (--reads) placement redundancy of the region's existing short reads

--panel FASTA judges a graph whose paths are an arbitrary panel instead of hap32.fa, e.g. a
full-panel graph (work/panel/<method>/<id>.gfa, one path per record of regions/<id>/hprc.fa.gz;
see tools/panel.py). Every FASTA record must be a path of the same name spelling it exactly
(exit 3 otherwise); the CHM13 record is the reference. Only size, alignment, inflation and
redundancy are computed: size, bubbles and k-mer redundancy on the distinct walks (they do not
depend on multiplicity), the alignment metrics on seeded samples of --max-pairs (default 300)
uniform haplotype pairs (all_*) and CHM13 x --max-pairs sampled haplotypes (ref_*), inflation
over every haplotype. Output defaults to results/panel/<METHOD>/<region_id>.json.

Method-independent optima (pairwise and truth-vs-panel distances) are cached under
config.WORK_DIR/cache/evaluate/<region_id>/ so that comparing many graphs of a region pays for
them once. truvari and bcftools work files go to config.WORK_DIR/evaluate/<METHOD>/<region_id>/.
"""
import argparse
import bisect
import collections
import concurrent.futures
import ctypes
import gzip
import hashlib
import json
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
import time

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402

REPO = os.path.dirname(TOOLS)
LIB_PATH = os.path.join(TOOLS, 'bin', 'libvntreval.so')
TRUVARI = os.environ.get('VNTR_TRUVARI', os.path.join(config.EVAL_DIR, 'work/truvari-venv/bin/truvari'))
AFFINE = (4, 6, 2)          # mismatch, gap open, gap extend: a gap of k bases costs 6 + 2k
K = 21                       # k-mer length for the redundancy metric
SV = 50                      # SV size threshold (bp)
READ_MIN = 30                # minimum aligned graph sequence of a read to test its placements
VERSION = 2

_RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_RC)[::-1]


def log(msg):
    print('[evaluate %6.1fs] %s' % (time.time() - T0, msg), file=sys.stderr, flush=True)


T0 = time.time()


# ----------------------------------------------------------------------------- C helpers

class CLib:
    """ctypes wrapper of tools/bin/libvntreval.so (built from tools/c by make on first use)."""

    def __init__(self):
        src = os.path.join(TOOLS, 'c', 'vntreval.c')
        if not os.path.exists(LIB_PATH) or os.path.getmtime(LIB_PATH) < os.path.getmtime(src):
            subprocess.run(['make', '-s', '-C', os.path.join(TOOLS, 'c')], check=True,
                           env=dict(os.environ, CC=config.CC))
        lib = ctypes.CDLL(LIB_PATH)
        lib.ed_unit.restype = ctypes.c_long
        lib.ed_unit.argtypes = [ctypes.c_char_p, ctypes.c_long, ctypes.c_char_p, ctypes.c_long]
        lib.ed_affine.restype = ctypes.c_long
        lib.ed_affine.argtypes = [ctypes.c_char_p, ctypes.c_long, ctypes.c_char_p, ctypes.c_long,
                                  ctypes.c_int, ctypes.c_int, ctypes.c_int]
        lib.dag_align.restype = ctypes.c_long
        lib.dag_align.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_long),
                                  ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                  ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int),
                                  ctypes.c_char_p, ctypes.c_long, ctypes.c_long,
                                  ctypes.c_double, ctypes.POINTER(ctypes.c_int),
                                  ctypes.POINTER(ctypes.c_int)]
        self.lib = lib

    def ed(self, a, b):
        """exact global unit-cost edit distance of two bytes objects"""
        if not a:
            return len(b)
        if not b:
            return len(a)
        if a == b:
            return 0
        # the pattern is the shorter one (fewer 64-bit blocks)
        if len(a) > len(b):
            a, b = b, a
        return self.lib.ed_unit(a, len(a), b, len(b))

    def affine(self, a, b):
        """exact global gap-affine cost (mismatch 4, gap 6 + 2k)"""
        if a == b:
            return 0
        r = self.lib.ed_affine(a, len(a), b, len(b), *AFFINE)
        if r < 0:
            raise MemoryError('ed_affine')
        return r


CL = None


def clib():
    global CL
    if CL is None:
        CL = CLib()
    return CL


# ----------------------------------------------------------------------------- inputs

def read_fasta(path):
    """ordered list of (name, SEQUENCE); name = first word of the header"""
    out = []
    if not os.path.exists(path):
        return out
    op = gzip.open if path.endswith('.gz') else open
    name, buf = None, []
    with op(path, 'rt') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ''.join(buf).upper()))
                name, buf = line[1:].split()[0] if line[1:].strip() else '', []
            elif line:
                buf.append(line.strip())
    if name is not None:
        out.append((name, ''.join(buf).upper()))
    return out


class Graph:
    """A GFA 1.0 graph: segments, links and paths (P or W lines).
    Nodes are indexed 0..N-1; a handle is 2*index + (1 if reverse)."""

    def __init__(self, path):
        self.names = []          # segment name per index
        self.seqs = []           # forward sequence per index (upper case)
        idx = {}
        links = []
        praw = []
        op = gzip.open if path.endswith('.gz') else open
        with op(path, 'rt') as f:
            for line in f:
                if not line or line[0] not in 'SLPW':
                    continue
                x = line.rstrip('\n\r').split('\t')
                if x[0] == 'S':
                    idx[x[1]] = len(self.names)
                    self.names.append(x[1])
                    self.seqs.append(x[2].upper())
                elif x[0] == 'L':
                    links.append((x[1], x[2], x[3], x[4]))
                elif x[0] == 'P':
                    steps = [(s[:-1], s[-1]) for s in x[2].split(',') if s]
                    praw.append((x[1], steps, None))
                elif x[0] == 'W':
                    steps = [(n, '+' if o == '>' else '-') for o, n in re.findall(r'([<>])([^<>]+)', x[6])]
                    praw.append(('%s#%s#%s' % (x[1], x[2], x[3]), steps, (x[4], x[5])))
        self.idx = idx
        self.n = len(self.names)
        self.links = set()
        self.bad_links = 0
        for a, ao, b, bo in links:
            if a not in idx or b not in idx:
                self.bad_links += 1
                continue
            ha = 2 * idx[a] + (ao == '-')
            hb = 2 * idx[b] + (bo == '-')
            self.links.add(self.canon_edge(ha, hb))
        self.paths = collections.OrderedDict()   # name -> list of handles
        self.path_range = {}
        self.aliases = {}                        # path key -> names it may be listed under
        self.bad_steps = []
        for name, steps, rng in praw:
            hs = []
            for n, o in steps:
                if n not in idx:
                    self.bad_steps.append(name)
                    hs = None
                    break
                hs.append(2 * idx[n] + (o == '-'))
            if hs is None:
                continue
            key = name
            k = 1
            while key in self.paths:          # W lines of one sample/hap/contig: disambiguate
                k += 1
                key = '%s#%d' % (name, k)
            self.paths[key] = hs
            al = {key, name}
            if rng:
                self.path_range[key] = rng
                # vg writes a PanSN path 'S#H#C[start]' (a subrange) as a W line with that start, and
                # 'S#H#C#0' (phase block 0) as a W line 'S H C' (the block number is dropped)
                for b in (name, name + '#0'):
                    al.update((b, '%s[%s]' % (b, rng[0]), '%s[%s-%s]' % (b, rng[0], rng[1])))
            self.aliases[key] = al
        self._oseq = {}

    @staticmethod
    def canon_edge(ha, hb):
        """an edge ha->hb equals (hb^1)->(ha^1); keep the smaller representation"""
        e1, e2 = (ha, hb), (hb ^ 1, ha ^ 1)
        return min(e1, e2)

    def hseq(self, h):
        s = self._oseq.get(h)
        if s is None:
            s = self.seqs[h >> 1]
            if h & 1:
                s = revcomp(s)
            self._oseq[h] = s
        return s

    def hlen(self, h):
        return len(self.seqs[h >> 1])

    def spell(self, hs):
        return ''.join(self.hseq(h) for h in hs)

    def hname(self, h):
        return ('<' if h & 1 else '>') + self.names[h >> 1]


SUBRANGE_RE = re.compile(r'\[\d+(?:-\d+)?\]$')


def match_paths(g, names):
    """Map each hap32 name to a graph path. In order: the exact name; a W line's aliases
    (sample#hap#contig, with '[start]' or '[start-end]' appended, which is how vg writes a PanSN
    name carrying a subrange such as GRCh38#0#chr4[58921381]); the same with any trailing
    '[...]' subrange removed from both names; a unique graph path whose name extends the hap32
    name after a '#', '[' or ':' delimiter. A W line also answers to 'sample#hap#contig#0' (vg drops
    phase block 0 when it writes a W line). Each graph path may be claimed by one name only."""
    out, missing = {}, []
    keys = list(g.paths)
    alias = collections.defaultdict(set)
    stripped = collections.defaultdict(set)
    for k in keys:
        for a in g.aliases.get(k, {k}):
            alias[a].add(k)
            stripped[SUBRANGE_RE.sub('', a)].add(k)
    used = set()
    for nm in names:
        hit = None
        for cand in (alias.get(nm, set()), stripped.get(SUBRANGE_RE.sub('', nm), set()),
                     set(k for k in keys if k.startswith(nm) and k[len(nm):len(nm) + 1] in ('#', '[', ':'))):
            cand = cand - used
            if len(cand) == 1:
                hit = next(iter(cand))
                break
            if len(cand) > 1:
                break           # ambiguous: do not guess
        if hit is None:
            missing.append(nm)
        else:
            out[nm] = hit
            used.add(hit)
    return out, missing


# ----------------------------------------------------------------------------- cache

class Cache:
    """method-independent optima of one region, keyed by sequence hashes (JSON on disk)"""

    def __init__(self, region_id, enabled=True):
        self.path = os.path.join(config.WORK_DIR, 'cache', 'evaluate', region_id, 'optima.json')
        self.d = {}
        self.dirty = False
        self.enabled = enabled
        if enabled and os.path.exists(self.path):
            try:
                self.d = json.load(open(self.path))
            except (ValueError, OSError):
                self.d = {}

    @staticmethod
    def key(kind, a, b):
        ha = hashlib.sha1(a).hexdigest()[:20]
        hb = hashlib.sha1(b).hexdigest()[:20]
        if hb < ha:
            ha, hb = hb, ha
        return '%s:%s:%s' % (kind, ha, hb)

    def get(self, k):
        return self.d.get(k)

    def put(self, k, v):
        self.d[k] = v
        self.dirty = True

    def save(self):
        if not (self.enabled and self.dirty):
            return
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        tmp = self.path + '.%d.tmp' % os.getpid()
        with open(tmp, 'w') as f:
            json.dump(self.d, f)
        os.replace(tmp, self.path)


# ----------------------------------------------------------------------------- DAG

def build_dag(g, handle_paths):
    """Handles between the paths' first and last handles, from the L lines plus the paths' own
    steps. Returns dict with topological order (None when cyclic) and predecessor lists."""
    succ = collections.defaultdict(set)
    pred = collections.defaultdict(set)
    missing_edges = set()
    edges = set(g.links)
    for hs in handle_paths:
        for a, b in zip(hs, hs[1:]):
            e = Graph.canon_edge(a, b)
            if e not in edges:
                missing_edges.add(e)
    for a, b in edges | missing_edges:
        for x, y in ((a, b), (b ^ 1, a ^ 1)):
            succ[x].add(y)
            pred[y].add(x)
    sources = set(hs[0] for hs in handle_paths)
    sinks = set(hs[-1] for hs in handle_paths)

    def reach(starts, nbr):
        seen = set(starts)
        st = list(starts)
        while st:
            v = st.pop()
            for w in nbr[v]:
                if w not in seen:
                    seen.add(w)
                    st.append(w)
        return seen
    keep = reach(sources, succ) & reach(sinks, pred)
    indeg = {h: 0 for h in keep}
    for h in keep:
        for w in succ[h]:
            if w in keep:
                indeg[w] += 1
    q = collections.deque(sorted(h for h in keep if indeg[h] == 0))
    order = []
    while q:
        v = q.popleft()
        order.append(v)
        for w in sorted(succ[v]):
            if w in keep:
                indeg[w] -= 1
                if indeg[w] == 0:
                    q.append(w)
    acyclic = len(order) == len(keep)
    return {'succ': succ, 'pred': pred, 'keep': keep, 'order': order if acyclic else None,
            'acyclic': acyclic, 'sources': sources, 'sinks': sinks,
            'missing_edges': len(missing_edges),
            'both_orientations': sum(1 for h in keep if (h ^ 1) in keep) // 2}


def dag_align(g, dag, query, max_bytes, prio=None, ub=-1):
    """min unit-cost edit distance of `query` to any source->sink path; returns (d, handles, err).
    prio: handle -> tie-break priority of the traceback (higher wins; default 0).
    ub: an upper bound on the distance (e.g. to the closest panel path); cells above it are pruned,
    which keeps the result exact and makes the DP banded. -1 = no bound."""
    order = dag['order']
    pos = {h: i for i, h in enumerate(order)}
    H = len(order)
    seqparts, off = [], [0]
    pstart, padj = [0], []
    src = bytearray(H)
    snk = bytearray(H)
    for i, h in enumerate(order):
        s = g.hseq(h).encode()
        seqparts.append(s)
        off.append(off[-1] + len(s))
        for p in dag['pred'][h]:
            if p in pos:
                padj.append(pos[p])
        pstart.append(len(padj))
        src[i] = h in dag['sources']
        snk[i] = h in dag['sinks']
    seq = b''.join(seqparts)
    pr = (ctypes.c_int * max(1, H))(*[(prio or {}).get(h, 0) for h in order])
    Off = (ctypes.c_long * len(off))(*off)
    Ps = (ctypes.c_int * len(pstart))(*pstart)
    Pa = (ctypes.c_int * max(1, len(padj)))(*padj)
    out = (ctypes.c_int * (H + 1))()
    npth = ctypes.c_int(0)
    q = query.encode()
    d = clib().lib.dag_align(H, seq, Off, Ps, Pa, bytes(src), bytes(snk), pr, q, len(q), int(ub),
                             float(max_bytes), out, ctypes.byref(npth))
    if d < 0:
        return None, None, {-2: 'dp_too_large', -3: 'distance_ge_65535', -4: 'traceback_error',
                            -5: 'no_source_to_sink', -6: 'upper_bound_violated'}.get(d, 'error%d' % d)
    return d, [order[out[i]] for i in range(npth.value)], None


def bubbles(g, dag, handle_paths, ref_hs):
    """Top-level bubbles: a handle is a cut when no edge jumps over it in topological order (then
    every source->sink path passes through it; a virtual source/sink joins the paths' ends).
    Consecutive cuts that are not adjacent bound a bubble. Alleles = distinct traversals by the
    paths between the two cuts."""
    order = dag['order']
    if order is None:
        return None
    pos = {h: i + 1 for i, h in enumerate(order)}      # 0 = virtual source, n+1 = virtual sink
    n = len(order) + 2
    diff = [0] * (n + 1)

    def add(i, j):
        if j - i > 1:
            diff[i + 1] += 1
            diff[j] -= 1
    for h in order:
        for w in dag['succ'][h]:
            if w in pos:
                add(pos[h], pos[w])
    for h in dag['sources']:
        add(0, pos[h])
    for h in dag['sinks']:
        add(pos[h], n - 1)
    cuts = []
    c = 0
    for i in range(n):
        c += diff[i]
        if c == 0:
            cuts.append(i)
    cut_set = set(cuts)
    inv = [None] + list(order) + [None]
    rpos = {}
    o = 0
    for h in ref_hs:
        rpos[h] = (o, o + g.hlen(h))
        o += g.hlen(h)
    alle = collections.defaultdict(collections.Counter)
    for hs in handle_paths:
        idxs = [(k, pos[h]) for k, h in enumerate(hs) if pos.get(h) in cut_set]
        full = [(-1, 0)] + idxs + [(len(hs), n - 1)]
        for (ka, pa), (kb, pb) in zip(full, full[1:]):
            if pb - pa > 1:
                alle[(pa, pb)][tuple(hs[ka + 1:kb])] += 1
    res = []
    for a, b in zip(cuts, cuts[1:]):
        if b - a <= 1:
            continue
        als = alle.get((a, b), collections.Counter())
        lens = sorted(set(sum(g.hlen(h) for h in al) for al in als))
        ha, hb = inv[a], inv[b]
        rs = rpos[ha][1] if ha in rpos else (0 if a == 0 else None)
        re_ = rpos[hb][0] if hb in rpos else (o if b == n - 1 else None)
        rbp = re_ - rs if rs is not None and re_ is not None else None
        res.append({'inner_handles': b - a - 1, 'n_alleles': len(als), 'ref_bp': rbp,
                    'allele_len_min': lens[0] if lens else None, 'allele_len_max': lens[-1] if lens else None,
                    'max_len_diff_vs_ref': max((abs(x - rbp) for x in lens), default=0) if rbp is not None else None,
                    'size': max(rbp or 0, lens[-1] if lens else 0)})
    return res


# ----------------------------------------------------------------------------- metrics

def q_stats(xs, pre=''):
    xs = [x for x in xs if x is not None]
    if not xs:
        return {pre + 'n': 0}
    s = sorted(xs)

    def qq(p):
        k = (len(s) - 1) * p
        lo = int(k)
        hi = min(lo + 1, len(s) - 1)
        return s[lo] + (s[hi] - s[lo]) * (k - lo)
    return {pre + 'n': len(s), pre + 'median': round(qq(0.5), 4), pre + 'mean': round(sum(s) / len(s), 4),
            pre + 'q1': round(qq(0.25), 4), pre + 'q3': round(qq(0.75), 4), pre + 'max': round(s[-1], 4)}


def size_metrics(g, paths, ref, span_bp, dag):
    used = set()
    edges = set()
    for hs in paths.values():
        used.update(h >> 1 for h in hs)
        for a, b in zip(hs, hs[1:]):
            edges.add(Graph.canon_edge(a, b))
    lens = sorted(g.hlen(2 * i) for i in used)
    ref_nodes = set(h >> 1 for h in paths[ref])
    nonref = [i for i in used if i not in ref_nodes]
    nonref_bp = sum(g.hlen(2 * i) for i in nonref)
    walks = set(tuple(hs) for hs in paths.values())
    seqs = set(g.spell(hs) for hs in paths.values())
    out = {
        'nodes': len(used), 'edges': len(edges), 'node_bp': sum(lens),
        'gfa_nodes': g.n, 'gfa_edges': len(g.links), 'gfa_bp': sum(len(s) for s in g.seqs),
        'span_bp': span_bp,
        'nodes_per_kb': round(1000.0 * len(used) / span_bp, 2),
        'edges_per_kb': round(1000.0 * len(edges) / span_bp, 2),
        'node_len_median': statistics.median(lens) if lens else None,
        'node_len_mean': round(sum(lens) / len(lens), 3) if lens else None,
        'node_frac_1bp': round(sum(1 for x in lens if x == 1) / len(lens), 4) if lens else None,
        'node_frac_le4bp': round(sum(1 for x in lens if x <= 4) / len(lens), 4) if lens else None,
        'nonref_nodes': len(nonref), 'nonref_bp': nonref_bp,
        'nonref_bp_per_kb': round(1000.0 * nonref_bp / span_bp, 2),
        'distinct_walks': len(walks), 'distinct_seqs': len(seqs),
        'same_seq_different_walk': len(walks) - len(seqs),
        'acyclic': dag['acyclic'], 'handles_on_paths_dag': len(dag['keep']),
        'nodes_both_orientations': dag['both_orientations'],
        'path_edges_missing_from_L': dag['missing_edges'],
    }
    bl = bubbles(g, dag, list(paths.values()), paths[ref])
    if bl is not None:
        out['bubbles'] = len(bl)
        out['bubbles_per_kb'] = round(1000.0 * len(bl) / span_bp, 3)
        out['bubbles_sv'] = sum(1 for b in bl if (b['max_len_diff_vs_ref'] or 0) >= SV)
        out['bubbles_ge5_alleles'] = sum(1 for b in bl if b['n_alleles'] >= 5)
        out['max_alleles_any_bubble'] = max((b['n_alleles'] for b in bl), default=0)
        if bl:
            big = max(bl, key=lambda b: (b['size'], b['n_alleles']))
            out['largest_bubble_size'] = big['size']
            out['largest_bubble_alleles'] = big['n_alleles']
            out['largest_bubble_ref_bp'] = big['ref_bp']
            out['largest_bubble_inner_handles'] = big['inner_handles']
            out['largest_bubble_frac_of_span'] = round((big['ref_bp'] or 0) / span_bp, 4)
        else:
            out.update({'largest_bubble_size': 0, 'largest_bubble_alleles': 0, 'largest_bubble_ref_bp': 0,
                        'largest_bubble_inner_handles': 0, 'largest_bubble_frac_of_span': 0})
    return out


def segments(a, b, ca, cb, pos_b):
    """Segments of walks a, b between consecutive shared handles: list of (a0, a1, b0, b1) base
    intervals (half-open) that are not on a shared handle. ca/cb: cumulative base offsets."""
    segs = []
    pa = pb = -1
    last = -1
    for i, h in enumerate(a):
        j = pos_b.get(h)
        if j is None or j <= last:
            continue
        if i - pa > 1 or j - pb > 1:
            segs.append((ca[pa + 1], ca[i], cb[pb + 1], cb[j]))
        pa, pb, last = i, j, j
    if len(a) - pa > 1 or len(b) - pb > 1:
        segs.append((ca[pa + 1], ca[len(a)], cb[pb + 1], cb[len(b)]))
    return segs


def pair_metrics(g, names, paths, ref, cache, max_pairs, affine_budget, seed, threads, panel_sample=False):
    """graph-implied vs optimal pairwise alignment cost, unit and gap-affine.

    panel_sample (--panel, a panel of hundreds of paths): the `all_*` metrics are aggregated over a
    seeded uniform sample of at most max_pairs path pairs, and the `ref_*` metrics over CHM13
    paired with a seeded sample of at most max_pairs other paths; both samples are scored in one
    pass but aggregated apart, so neither is biased by the other."""
    C = clib()
    N = len(names)
    walks = [paths[n] for n in names]
    if panel_sample:                     # a panel repeats walks: spell and index each one once
        memo = {}
        for hs in walks:
            k = tuple(hs)
            if k not in memo:
                c = [0]
                for h in hs:
                    c.append(c[-1] + g.hlen(h))
                memo[k] = (g.spell(hs).encode(), c, {h: i for i, h in enumerate(hs)})
        seqs = [memo[tuple(hs)][0] for hs in walks]
        cum = [memo[tuple(hs)][1] for hs in walks]
        posd = [memo[tuple(hs)][2] for hs in walks]
    else:
        seqs = [g.spell(paths[n]).encode() for n in names]
        cum = []
        for hs in walks:
            c = [0]
            for h in hs:
                c.append(c[-1] + g.hlen(h))
            cum.append(c)
        posd = [{h: i for i, h in enumerate(hs)} for hs in walks]
    r = names.index(ref)
    all_set = ref_set = None
    if panel_sample:
        rng = random.Random(seed)
        tot = N * (N - 1) // 2
        k = min(max_pairs, tot)
        all_set = set()
        while len(all_set) < k:         # uniform over unordered pairs, without listing them all
            i, j = rng.randrange(N), rng.randrange(N)
            if i != j:
                all_set.add((min(i, j), max(i, j)))
        others = [x for x in range(N) if x != r]
        ref_set = set((min(r, x), max(r, x)) for x in rng.sample(others, min(max_pairs, len(others))))
        pairs = sorted(all_set | ref_set)
        allp = range(tot)
        sampled = k < tot
    else:
        allp = [(i, j) for i in range(N) for j in range(i + 1, N)]
        if len(allp) <= max_pairs:
            pairs = allp
            sampled = False
        else:
            rng = random.Random(seed)
            refp = [p for p in allp if r in p]
            rest = [p for p in allp if r not in p]
            pairs = sorted(refp + rng.sample(rest, max(0, max_pairs - len(refp))))
            sampled = True

    # optimal costs (method independent; cached)
    def opt(kind, i, j):
        a, b = seqs[i], seqs[j]
        if a == b:
            return 0
        k = Cache.key(kind, a, b)
        v = cache.get(k)
        if v is None:
            v = C.ed(a, b) if kind == 'u' else C.affine(a, b)
            cache.put(k, v)
        return v

    # which pairs get the affine computation: decided from the sequences alone (never from what
    # the cache holds), so every graph of a region is scored on the same pairs
    need = set()
    for i, j in pairs:
        if seqs[i] != seqs[j]:
            need.add(Cache.key('a', seqs[i], seqs[j]) + ':%d' % (len(seqs[i]) * len(seqs[j])))
    cells = sum(int(k.rsplit(':', 1)[1]) for k in need)
    ref_pairs = [p for p in pairs if r in p]
    if cells <= affine_budget:
        aff_pairs = list(pairs)
        aff_mode = 'all_pairs'
    else:
        rng = random.Random(seed + 1)
        rest = [p for p in pairs if r not in p]
        # add as many non-reference pairs as fit the budget
        ref_cells = sum(len(seqs[i]) * len(seqs[j]) for i, j in ref_pairs)
        per = (sum(len(seqs[i]) * len(seqs[j]) for i, j in rest) / len(rest)) if rest else 1
        kmax = int(max(0, affine_budget - ref_cells) / max(per, 1))
        aff_pairs = sorted(ref_pairs + rng.sample(rest, min(len(rest), kmax)))
        aff_mode = 'ref_pairs_plus_%d_sampled' % min(len(rest), kmax)
    aff_set = set(aff_pairs)

    # compute the optima in parallel (ctypes releases the GIL)
    jobs = []
    for i, j in pairs:
        if seqs[i] != seqs[j]:
            jobs.append(('u', i, j))
            if (i, j) in aff_set:
                jobs.append(('a', i, j))
    jobs = [jb for jb in jobs if cache.get(Cache.key(jb[0], seqs[jb[1]], seqs[jb[2]])) is None]
    # largest first for load balance; deduplicate identical sequence pairs
    seen = set()
    uj = []
    for jb in sorted(jobs, key=lambda x: -len(seqs[x[1]]) * len(seqs[x[2]])):
        k = Cache.key(jb[0], seqs[jb[1]], seqs[jb[2]])
        if k not in seen:
            seen.add(k)
            uj.append(jb)
    if uj:
        t = time.time()
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
            futs = {ex.submit(C.ed if kind == 'u' else C.affine, seqs[i], seqs[j]): (kind, i, j)
                    for kind, i, j in uj}
            for f in concurrent.futures.as_completed(futs):
                kind, i, j = futs[f]
                cache.put(Cache.key(kind, seqs[i], seqs[j]), f.result())
        log('optimal pairwise costs: %d computed in %.1fs' % (len(uj), time.time() - t))

    # graph-implied costs: segments between shared handles, each segment pair aligned optimally
    t = time.time()
    pair_keys = {}
    need_u, need_a = set(), set()
    for i, j in pairs:
        if walks[i] == walks[j]:
            continue
        keys = []
        for a0, a1, b0, b1 in segments(walks[i], walks[j], cum[i], cum[j], posd[j]):
            sa, sb = seqs[i][a0:a1], seqs[j][b0:b1]
            key = (sa, sb) if sa <= sb else (sb, sa)
            keys.append(key)
            need_u.add(key)
            if (i, j) in aff_set:
                need_a.add(key)
        pair_keys[(i, j)] = keys
    segu, sega = {}, {}
    big = []
    for key in need_u:
        if len(key[0]) * len(key[1]) > 200000:
            big.append(('u', key))
        else:
            segu[key] = C.ed(*key)
    for key in need_a:
        if len(key[0]) * len(key[1]) > 20000:
            big.append(('a', key))
        else:
            sega[key] = C.affine(*key)
    if big:
        big.sort(key=lambda x: -len(x[1][0]) * len(x[1][1]))
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
            futs = {ex.submit(C.ed if kind == 'u' else C.affine, *key): (kind, key) for kind, key in big}
            for f in concurrent.futures.as_completed(futs):
                kind, key = futs[f]
                (segu if kind == 'u' else sega)[key] = f.result()
    rows = []
    for i, j in pairs:
        same_seq = seqs[i] == seqs[j]
        L = (len(seqs[i]) + len(seqs[j])) / 2.0
        if walks[i] == walks[j]:
            rows.append({'i': i, 'j': j, 'cost': 0, 'dopt': 0, 'U': 0, 'segs': 0, 'same_seq': True, 'same_walk': True,
                         'L': L,
                         'acost': 0 if (i, j) in aff_set else None, 'aopt': 0 if (i, j) in aff_set else None})
            continue
        keys = pair_keys[(i, j)]
        cost = sum(segu[k] for k in keys)
        U = sum(max(len(k[0]), len(k[1])) - segu[k] for k in keys)
        acost = sum(sega[k] for k in keys) if (i, j) in aff_set else None
        rows.append({'i': i, 'j': j, 'cost': cost, 'dopt': opt('u', i, j), 'U': U, 'segs': len(keys),
                     'same_seq': same_seq, 'same_walk': False, 'L': L,
                     'acost': acost, 'aopt': opt('a', i, j) if acost is not None else None})
    log('graph-implied costs over %d pairs (%d unit, %d affine segment pairs) in %.1fs' % (
        len(pairs), len(need_u), len(need_a), time.time() - t))

    def agg(sub, tag):
        """ratio of sums (as topology/topo.py), per-pair ratios, per-kb excess"""
        o = {}
        cost = sum(x['cost'] for x in sub)
        dopt = sum(x['dopt'] for x in sub)
        Ld = sum(x['L'] for x in sub if not x['same_walk'])   # as topology/topo.py: pairs of distinct walks
        o['pairs'] = len(sub)
        o['pairs_distinct_seq'] = sum(1 for x in sub if not x['same_seq'])
        o['pairs_same_seq_diff_walk'] = sum(1 for x in sub if x['same_seq'] and not x['same_walk'])
        o['cost_over_opt'] = round(cost / dopt, 4) if dopt else (1.0 if cost == 0 else None)
        o['excess_per_kb'] = round(1000.0 * (cost - dopt) / Ld, 2) if Ld else 0.0
        o['opt_per_kb'] = round(1000.0 * dopt / Ld, 2) if Ld else 0.0
        o['U_per_kb'] = round(1000.0 * sum(x['U'] for x in sub) / Ld, 2) if Ld else 0.0
        o['segments_per_kb'] = round(1000.0 * sum(x['segs'] for x in sub) / Ld, 2) if Ld else 0.0
        rat = [x['cost'] / x['dopt'] for x in sub if x['dopt'] > 0]
        o.update(q_stats(rat, 'pair_ratio_'))
        dw = [x for x in sub if not x['same_walk']]
        o['frac_pairs_optimal'] = round(sum(1 for x in dw if x['cost'] == x['dopt']) / len(dw), 4) if dw else None
        return {tag + '_' + k: v for k, v in o.items()}

    def agg_aff(sub, tag):
        sub = [x for x in sub if x['acost'] is not None]
        o = {}
        cost = sum(x['acost'] for x in sub)
        opt_ = sum(x['aopt'] for x in sub)
        Ld = sum(x['L'] for x in sub if not x['same_walk'])   # as topology/topo.py: pairs of distinct walks
        o['pairs'] = len(sub)
        o['cost_over_opt'] = round(cost / opt_, 4) if opt_ else (1.0 if cost == 0 else None)
        o['excess_per_kb'] = round(1000.0 * (cost - opt_) / Ld, 2) if Ld else 0.0
        o['opt_per_kb'] = round(1000.0 * opt_ / Ld, 2) if Ld else 0.0
        rat = [x['acost'] / x['aopt'] for x in sub if x['aopt'] > 0]
        o.update(q_stats(rat, 'pair_ratio_'))
        return {tag + '_' + k: v for k, v in o.items()}

    res = {'n_paths': N, 'pairs_total': len(allp), 'pairs_sampled': sampled, 'affine_pairs': aff_mode}
    if panel_sample:
        res['pairs_mode'] = 'panel: %d uniform pairs (all_*), CHM13 x %d sampled paths (ref_*), seed %d' % (
            len(all_set), len(ref_set), seed)
        all_rows = [x for x in rows if (x['i'], x['j']) in all_set]
        ref_rows = [x for x in rows if (x['i'], x['j']) in ref_set]
    else:
        all_rows = rows
        ref_rows = [x for x in rows if r in (x['i'], x['j'])]
    res.update(agg(all_rows, 'all'))
    res.update(agg(ref_rows, 'ref'))
    res.update(agg_aff(all_rows, 'affine_all'))
    res.update(agg_aff(ref_rows, 'affine_ref'))
    if panel_sample:
        # sampling uncertainty of the ratio of sums: 95% bootstrap interval over the sampled pairs
        # (two graphs of one panel FASTA are scored on the same pairs, so their difference is paired)
        brng = random.Random(seed + 2)
        for sub, tag, ck, ok in ((all_rows, 'all', 'cost', 'dopt'), (ref_rows, 'ref', 'cost', 'dopt'),
                                 (all_rows, 'affine_all', 'acost', 'aopt'), (ref_rows, 'affine_ref', 'acost', 'aopt')):
            xs = [(x[ck], x[ok]) for x in sub if x[ck] is not None]
            if not xs:
                continue
            bs = []
            for _ in range(1000):
                c = o = 0
                for _k in range(len(xs)):
                    a, b = xs[brng.randrange(len(xs))]
                    c += a
                    o += b
                if o:
                    bs.append(c / o)
            bs.sort()
            if bs:
                res[tag + '_cost_over_opt_ci95'] = [round(bs[int(0.025 * (len(bs) - 1))], 4),
                                                    round(bs[int(0.975 * (len(bs) - 1))], 4)]
    dref = {}
    for x in rows:
        if r in (x['i'], x['j']):
            dref[x['j'] if x['i'] == r else x['i']] = x['dopt']
    return res, dref


def ref_coords(g, ref_hs):
    idx, off = {}, []
    o = 0
    for i, h in enumerate(ref_hs):
        if h not in idx:
            idx[h] = i
        off.append(o)
        o += g.hlen(h)
    off.append(o)
    return idx, off


def deconstruct(g, hs, ref_hs, refidx, refoff, refseq, left_base):
    """Records (offset0, REF, ALT) of walk `hs` against the CHM13 walk, one per stretch between
    consecutive shared CHM13 handles (visited in increasing CHM13 order, greedily), padded with
    the preceding reference base, common suffix then prefix trimmed (oracle_truvari.deconstruct).
    offset0 = 0-based offset of the padding base in the CHM13 walk (-1 = the base before it,
    `left_base`)."""
    anchors = [(-1, -1)]
    last = -1
    for wi, h in enumerate(hs):
        ri = refidx.get(h)
        if ri is not None and ri > last:
            anchors.append((wi, ri))
            last = ri
    anchors.append((len(hs), len(ref_hs)))
    recs = []
    for (wa, ra), (wb, rb) in zip(anchors, anchors[1:]):
        wseg = hs[wa + 1:wb]
        rseg = ref_hs[ra + 1:rb]
        if wseg == rseg:
            continue
        r0 = refoff[ra + 1]           # first base after the left anchor (0 when no anchor)
        r1 = refoff[rb]               # first base of the right anchor (len when none)
        pad = refseq[r0 - 1] if r0 > 0 else left_base
        ref = pad + refseq[r0:r1]
        alt = pad + g.spell(wseg)
        p0 = r0 - 1
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0] and ref[1] == alt[1]:
            ref, alt, p0 = ref[1:], alt[1:], p0 + 1
        if ref != alt:
            recs.append((p0, ref, alt))
    return recs


def inflation_metrics(g, names, paths, ref, dref, refseq):
    C = clib()
    ref_hs = paths[ref]
    refidx, refoff = ref_coords(g, ref_hs)
    rows = []
    edcache = {}
    walkmemo = {}                 # a panel repeats walks; each is deconstructed once
    for k, n in enumerate(names):
        if n == ref:
            continue
        hs = paths[n]
        wk = tuple(hs)
        if wk not in walkmemo:
            recs = deconstruct(g, hs, ref_hs, refidx, refoff, refseq, 'N')
            anch = 0
            for _, a, b in recs:
                key = (a, b)
                d = edcache.get(key)
                if d is None:
                    d = C.ed(a.encode(), b.encode())
                    edcache[key] = d
                anch += d
            L = sum(g.hlen(h) for h in hs)
            walkmemo[wk] = {'records': len(recs), 'sv': sum(1 for _, a, b in recs if abs(len(a) - len(b)) >= SV),
                            'indel_sum': sum(abs(len(a) - len(b)) for _, a, b in recs),
                            'net': abs(L - len(refseq)), 'anchored': anch, 'L': (L + len(refseq)) / 2.0}
        rows.append(dict(walkmemo[wk], free=dref.get(k)))
    out = {}
    out.update(q_stats([x['sv'] for x in rows], 'sv_pieces_'))
    out.update(q_stats([x['records'] for x in rows], 'records_'))
    out['paths_with_sv_piece'] = sum(1 for x in rows if x['sv'] > 0)
    out['sv_pieces_total'] = sum(x['sv'] for x in rows)
    isum = sum(x['indel_sum'] for x in rows)
    net = sum(x['net'] for x in rows)
    Lt = sum(x['L'] for x in rows)
    out['indel_bp_graph'] = isum
    out['indel_bp_net'] = net
    out['indel_bp_ratio'] = round(isum / net, 4) if net else (1.0 if isum == 0 else None)
    out['indel_bp_excess_per_kb'] = round(1000.0 * (isum - net) / Lt, 2) if Lt else 0.0
    nz = [x for x in rows if x['free']]
    fa = sum(x['free'] for x in nz)
    aa = sum(x['anchored'] for x in nz)
    out['anchored_over_free'] = round(aa / fa, 4) if fa else 1.0
    out.update(q_stats([x['anchored'] / x['free'] for x in nz], 'anchored_over_free_path_'))
    return out


def kmer_redundancy(g, paths, k=K):
    """P(K) = distinct graph positions (node, forward offset, orientation of the k-mer's first
    base) of k-mer K over all paths; C(K) = the most times one path contains K.
    extra = sum(P - C). (viewer/view_region.py kmer_redundancy)"""
    pos = collections.defaultdict(set)
    maxocc = collections.Counter()
    for hs in paths.values():
        seq = g.spell(hs)
        bp = []
        for h in hs:
            L = g.hlen(h)
            nid = h >> 1
            if h & 1:
                bp.extend((nid, L - 1 - j, 1) for j in range(L))
            else:
                bp.extend((nid, j, 0) for j in range(L))
        occ = collections.Counter()
        for x in range(len(seq) - k + 1):
            km = seq[x:x + k]
            if 'N' in km:
                continue
            pos[km].add(bp[x])
            occ[km] += 1
        for km, c in occ.items():
            if c > maxocc[km]:
                maxocc[km] = c
    P = sum(len(v) for v in pos.values())
    extra = sum(len(v) - maxocc[km] for km, v in pos.items())
    multi = sum(1 for km, v in pos.items() if len(v) > maxocc[km])
    return {'k': k, 'kmer_positions': P, 'kmer_extra_positions': extra,
            'kmer_frac_extra': round(extra / P, 5) if P else 0.0,
            'kmers': len(pos), 'kmers_with_extra': multi}


# ----------------------------------------------------------------------------- truth

def fetch_ref(contig, a1, b1):
    fa = config.data_paths(contig)['ref_fa']
    if not os.path.exists(fa):
        return None
    p = subprocess.run([config.SAMTOOLS, 'faidx', fa, '%s:%d-%d' % (contig, a1, b1)], capture_output=True,
                       text=True, env=config.tool_env())
    if p.returncode != 0:
        return None
    return ''.join(p.stdout.split('\n')[1:]).upper()


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, env=config.tool_env(), **kw)


def write_vcf(recs, contig, ref_fa, od, tag):
    raw = os.path.join(od, tag + '.raw.vcf')
    fai = {l.split('\t')[0]: l.split('\t')[1] for l in open(ref_fa + '.fai')}
    with open(raw, 'w') as f:
        f.write('##fileformat=VCFv4.2\n##contig=<ID=%s,length=%s>\n' % (contig, fai[contig]))
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n')
        for (p, ref, alt), hs in sorted(recs.items()):
            gt = '%d|%d' % (1 if 1 in hs else 0, 1 if 2 in hs else 0)
            f.write('%s\t%d\t.\t%s\t%s\t60\tPASS\t.\tGT\t%s\n' % (contig, p, ref, alt, gt))
    vcf = os.path.join(od, tag + '.vcf.gz')
    norm = sh([config.BCFTOOLS, 'norm', '-m-any', '-f', ref_fa, raw]).stdout
    subprocess.run([config.BCFTOOLS, 'sort', '-Oz', '-o', vcf], input=norm, text=True, capture_output=True,
                   check=True, env=config.tool_env())
    sh([config.BCFTOOLS, 'index', '-t', '-f', vcf])
    return vcf


def f1(tpb, fn, tpc, fp):
    if None in (tpb, fn, tpc, fp):
        return None
    rec = tpb / (tpb + fn) if tpb + fn else None
    prec = tpc / (tpc + fp) if tpc + fp else None
    if rec is None and prec is None:
        return None
    if not rec or not prec:
        return 0.0
    return round(2 * rec * prec / (rec + prec), 4)


def truvari_score(recs, contig, a1, b1, od, refine, refine_timeout):
    """truvari bench (pipeline parameters) of the records against the scored truth in the span
    intersected with the SV benchmark; then truvari refine -u (phab, POA) on its default candidate
    regions (refined_*: as the investigation), and truvari refine -u -w over the whole region
    (phab_*: counted in the harmonised representation)."""
    p = config.data_paths(contig)
    ref_fa = p['ref_fa']
    truth_norm = os.path.join(config.EVAL_DIR, 'work/wgs-mm095/score/%s.truth.norm.vcf.gz' % contig)
    bench_bed = os.path.join(config.EVAL_DIR, 'work/wgs-mm095', contig, 'truth.%s.stvar.bed' % contig)
    if not os.path.exists(bench_bed):
        bench_bed = p['stvar_bed']
    for x in (ref_fa, truth_norm, bench_bed, TRUVARI):
        if not os.path.exists(x):
            return {'status': 'missing_input:%s' % x}
    if os.path.exists(od):
        shutil.rmtree(od)
    os.makedirs(od)
    iv = []
    with open(bench_bed) as f:
        for line in f:
            x = line.split('\t')
            if x[0] != contig:
                continue
            a, b = max(int(x[1]), a1 - 1), min(int(x[2]), b1)
            if a < b:
                iv.append((a, b))
    res = {'bench_bp': sum(b - a for a, b in iv)}
    if not iv:
        res['status'] = 'span_outside_sv_benchmark'
        return res
    bed = os.path.join(od, 'region.bed')
    with open(bed, 'w') as f:
        for a, b in iv:
            f.write('%s\t%d\t%d\n' % (contig, a, b))
    truth = os.path.join(od, 'truth.vcf.gz')
    sh([config.BCFTOOLS, 'view', '-r', '%s:%d-%d' % (contig, max(1, a1 - 2000), b1 + 2000), '-Oz', '-o', truth,
        truth_norm])
    sh([config.BCFTOOLS, 'index', '-t', '-f', truth])
    vcf = write_vcf(recs, contig, ref_fa, od, 'graph_truth')
    res['records'] = len(recs)
    res['records_sv'] = sum(1 for (p_, a, b) in recs if abs(len(a) - len(b)) >= SV)
    bd = os.path.join(od, 'bench')
    r = subprocess.run([TRUVARI, 'bench', '-b', truth, '-c', vcf, '-f', ref_fa, '-o', bd,
                        '--includebed', bed, '--sizemin', '50', '--sizefilt', '50',
                        '--bSample', 'HG002', '--cSample', 'HG002', '--pick', 'ac'],
                       capture_output=True, text=True, env=config.tool_env())
    sj = os.path.join(bd, 'summary.json')
    if r.returncode != 0 or not os.path.exists(sj):
        res['status'] = 'bench_failed: ' + r.stderr[-300:]
        return res
    s = json.load(open(sj))
    for k in ('TP-base', 'FN', 'TP-comp', 'FP'):
        res['raw_' + k] = s.get(k)
    res['raw_f1'] = f1(s.get('TP-base'), s.get('FN'), s.get('TP-comp'), s.get('FP'))
    res['status'] = 'ok'
    if refine:
        # a pristine copy of the bench directory for the whole-span phab scoring below
        pb = bd + '_phab'
        shutil.copytree(bd, pb)
        t = time.time()
        vs = os.path.join(bd, 'refine.variant_summary.json')
        try:
            subprocess.run([TRUVARI, 'refine', '-u', '-t', '1', '.'], cwd=bd, capture_output=True, text=True,
                           timeout=refine_timeout, env=config.tool_env())
        except subprocess.TimeoutExpired:
            res['refine_status'] = 'timeout'
        if os.path.exists(vs):
            s2 = json.load(open(vs))
            for k in ('TP-base', 'FN', 'TP-comp', 'FP'):
                res['refined_' + k] = s2.get(k)
            res['refined_f1'] = f1(s2.get('TP-base'), s2.get('FN'), s2.get('TP-comp'), s2.get('FP'))
            res['refine_status'] = 'ok'
        elif 'refine_status' not in res:
            rl = os.path.join(bd, 'refine.log.txt')
            if os.path.exists(rl) and 'No regions to be refined' in open(rl).read():
                # no FP/FN region to harmonise: the refined result is the raw one
                for k in ('TP-base', 'FN', 'TP-comp', 'FP', 'f1'):
                    res['refined_' + k] = res['raw_' + k]
                res['refine_status'] = 'nothing_to_refine'
            else:
                res['refine_status'] = 'failed'
        res['refine_seconds'] = round(time.time() - t, 1)
        # whole-span harmonisation counted in phab's representation (-w): the region is the whole
        # span (intersected with the SV benchmark), so the score depends only on the haplotype
        # sequences the graph writes, not on how the graph decomposes them
        t = time.time()
        vs = os.path.join(pb, 'refine.variant_summary.json')
        try:
            subprocess.run([TRUVARI, 'refine', '-u', '-w', '-t', '1', '--regions', os.path.abspath(bed), '.'],
                           cwd=pb, capture_output=True, text=True, timeout=refine_timeout, env=config.tool_env())
        except subprocess.TimeoutExpired:
            res['phab_status'] = 'timeout'
        if os.path.exists(vs):
            s3 = json.load(open(vs))
            for k in ('TP-base', 'FN', 'TP-comp', 'FP'):
                res['phab_' + k] = s3.get(k)
            res['phab_f1'] = f1(s3.get('TP-base'), s3.get('FN'), s3.get('TP-comp'), s3.get('FP'))
            res['phab_status'] = 'ok'
        elif 'phab_status' not in res:
            rl = os.path.join(pb, 'refine.log.txt')
            if os.path.exists(rl) and 'No regions to be refined' in open(rl).read():
                for k in ('TP-base', 'FN', 'TP-comp', 'FP', 'f1'):     # no FP/FN: phab = raw
                    res['phab_' + k] = res['raw_' + k]
                res['phab_status'] = 'nothing_to_refine'
            else:
                res['phab_status'] = 'failed'
        res['phab_seconds'] = round(time.time() - t, 1)
    return res


def truth_metrics(g, dag, names, paths, ref, truth, region, cache, work, args):
    C = clib()
    out = {}
    ref_hs = paths[ref]
    refidx, refoff = ref_coords(g, ref_hs)
    refseq = g.spell(ref_hs)
    pseqs = [g.spell(paths[n]).encode() for n in names]
    # traceback tie-break among equally good graph paths: CHM13 handles first, then handles on more
    # panel paths (so the truth is written as close to the reference as the graph allows, and the
    # choice does not depend on the order of the GFA's lines)
    prio = collections.Counter()
    for n in names:
        for h in set(paths[n]):
            prio[h] += 1
    for h in ref_hs:
        prio[h] += 1 << 20
    tpaths = {}
    for h, (tn, ts) in enumerate(truth[:2], start=1):
        tb = ts.encode()
        # closest panel path (method independent)
        k = Cache.key('tp', tb, hashlib.sha1(b''.join(sorted(set(pseqs)))).hexdigest().encode())
        v = cache.get(k)
        if v is None:
            best = None
            for nm, s in zip(names, pseqs):
                d = C.ed(tb, s)
                if best is None or d < best[0]:
                    best = (d, nm)
            v = list(best)
            cache.put(k, v)
        out['h%d_name' % h] = tn
        out['h%d_len' % h] = len(ts)
        out['h%d_d_panel' % h] = v[0]
        out['h%d_closest_panel_path' % h] = v[1]
        out['h%d_d_ref' % h] = C.ed(tb, refseq.encode())
        if dag['order'] is None:
            out['h%d_d_graph' % h] = None
            out['h%d_note' % h] = 'graph_cyclic'
            continue
        t = time.time()
        # the closest panel path is a source->sink path, so its distance bounds the DP (exact pruning);
        # try tighter bounds first (narrower bands), raising the bound while it proves too small
        ub = min(v[0], 64)
        while True:
            d, hp, err = dag_align(g, dag, ts, args.dp_max_bytes, prio, ub=ub)
            if err != 'upper_bound_violated' or ub >= v[0]:
                break
            ub = min(v[0], ub * 4)
        out['h%d_dp_bound' % h] = ub
        out['h%d_d_graph' % h] = d
        out['h%d_dp_seconds' % h] = round(time.time() - t, 2)
        if err:
            out['h%d_note' % h] = err
            continue
        # a valid graph path spells a sequence at distance d from the truth: check it
        ps = g.spell(hp)
        chk = C.ed(tb, ps.encode())
        if chk != d:
            out['h%d_note' % h] = 'traceback_distance_mismatch:%d' % chk
            continue
        out['h%d_graph_path_is_panel_path' % h] = int(any(paths[n] == hp for n in names))
        tpaths[h] = hp
    for h in (1, 2):
        dp, dg = out.get('h%d_d_panel' % h), out.get('h%d_d_graph' % h)
        if dp is not None and dg is not None:
            out['h%d_panel_minus_graph' % h] = dp - dg
    out['d_graph_sum'] = sum(out.get('h%d_d_graph' % h) or 0 for h in (1, 2)) \
        if all(out.get('h%d_d_graph' % h) is not None for h in (1, 2)) else None
    out['d_panel_sum'] = sum(out.get('h%d_d_panel' % h) or 0 for h in (1, 2))
    # the truth written by the graph
    if args.no_truvari:
        out['truvari'] = {'status': 'skipped'}
        return out
    if len(tpaths) < 2:
        out['truvari'] = {'status': 'no_truth_graph_paths'}
        return out
    contig = region['contig']
    a1, b1 = int(region['span_start']), int(region['span_end'])
    if len(refseq) != b1 - a1 + 1:
        out['truvari'] = {'status': 'chm13_path_length_%d_ne_span_%d' % (len(refseq), b1 - a1 + 1)}
        return out
    if not os.path.exists(config.data_paths(contig)['ref_fa']):
        # the truvari scoring needs the CHM13 FASTA, the normalised truth VCF and the SV benchmark
        # BED of the vg-call-eval repository; the region package alone does not carry them
        out['truvari'] = {'status': 'needs_big_data: no CHM13 FASTA at %s (set VNTR_EVAL_DIR, or pass '
                                    '--no-truvari)' % config.data_paths(contig)['ref_fa']}
        log('truvari scoring skipped: it needs the big data (CHM13 FASTA, truth VCF, benchmark BED '
            'under VNTR_EVAL_DIR=%s)' % config.EVAL_DIR)
        return out
    ctx = fetch_ref(contig, a1 - 1, b1)
    if ctx is None or ctx[1:] != refseq:
        out['truvari'] = {'status': 'chm13_path_differs_from_reference_span'}
        return out
    recs = {}
    for h, hp in tpaths.items():
        for p0, a, b in deconstruct(g, hp, ref_hs, refidx, refoff, refseq, ctx[0]):
            recs.setdefault((a1 + p0, a, b), set()).add(h)   # 1-based POS = a1 + p0
    t = time.time()
    out['truvari'] = truvari_score(recs, contig, a1, b1, os.path.join(work, 'truvari'), not args.no_refine,
                                   args.refine_timeout)
    out['truvari']['seconds'] = round(time.time() - t, 1)
    return out


# ----------------------------------------------------------------------------- reads (Stage 1)

CS_RE = re.compile(r'(:\d+|\*[a-zA-Z][a-zA-Z]|[+][a-zA-Z]+|-[a-zA-Z]+|=[a-zA-Z]+)')
STEP_RE = re.compile(r'([<>])(\d+)')


def fetch_reads(region, region_dir, mc):
    """The region's short reads as a GAF against the hap32 graph: <region_dir>/reads.gaf(.gz) if
    the package has it, else queried from GAF-Base between the anchors (cached)."""
    for nm in ('reads.gaf', 'reads.gaf.gz'):
        p = os.path.join(region_dir, nm)
        if os.path.exists(p):
            return p, 'package'
    cdir = os.path.join(config.WORK_DIR, 'cache', 'evaluate', region['region_id'])
    p = os.path.join(cdir, 'reads.gaf')
    if os.path.exists(p):
        return p, 'cache'
    dp = config.data_paths(region['contig'])
    if not (os.path.exists(dp['gaf_db']) and os.path.exists(dp['gbz_db'])):
        return None, 'no_gaf_base'
    os.makedirs(cdir, exist_ok=True)

    def tok(h):
        h = h.strip()
        if h[0] in '<>':
            return h[1:] + ('+' if h[0] == '>' else '-')
        return h
    btw = '%s:%s' % (tok(region['anchor_left']), tok(region['anchor_right']))
    tmp = p + '.%d.tmp' % os.getpid()
    subprocess.run([config.GBZ_BASE, 'query', dp['gbz_db'], '--between', btw, '--limit', '1000000',
                    '--gaf-base', dp['gaf_db'], '--gaf-output', tmp, '--alignments', 'overlapping'],
                   check=True, capture_output=True, text=True, env=config.tool_env())
    os.replace(tmp, p)
    return p, 'gaf_base'


def parse_gaf(path):
    op = gzip.open if path.endswith('.gz') else open
    out = []
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('@') or not line.strip():
                continue
            x = line.rstrip('\n').split('\t')
            if len(x) < 12:
                continue
            tags = {}
            for t in x[12:]:
                pp = t.split(':', 2)
                if len(pp) == 3:
                    tags[pp[0]] = pp[2]
            out.append({'name': x[0], 'qlen': int(x[1]), 'qs': int(x[2]), 'qe': int(x[3]),
                        'steps': [(n, o) for o, n in STEP_RE.findall(x[5])], 'plen': int(x[6]),
                        'ps': int(x[7]), 'pe': int(x[8]), 'mapq': int(x[11]), 'cs': tags.get('cs', '')})
    return out


def core_nodes(mc, mc_paths, ref, region):
    """Nodes of the baseline graph that belong to the core (reads_analysis.Region._core_nodes):
    CHM13 nodes overlapping the core, plus non-CHM13 nodes that some path places in the core
    when each run of non-CHM13 nodes is spread linearly between the CHM13 positions around it."""
    a1 = int(region['span_start'])
    a, b = int(region['core_start']) - a1, int(region['core_end']) - a1 + 1   # 0-based offsets in span
    ref_hs = mc_paths[ref]
    rpos = {}
    o = 0
    for h in ref_hs:
        rpos[h >> 1] = (o, o + mc.hlen(h))
        o += mc.hlen(h)
    span_len = o
    core = set(n for n, (s, e) in rpos.items() if s < b and e > a)

    def place(run, P, Q):
        lo, hi = min(P, Q), max(P, Q)
        tot = sum(mc.hlen(2 * n) for n in run) or 1
        off = 0
        for n in run:
            L = mc.hlen(2 * n)
            p = lo + (off + L / 2.0) / tot * (hi - lo)
            off += L
            if a <= p <= b:
                core.add(n)
    for nm, hs in mc_paths.items():
        if nm == ref:
            continue
        prev_end = 0
        pending = []
        for h in hs:
            n = h >> 1
            if n in rpos:
                if pending:
                    place(pending, prev_end, rpos[n][0])
                pending = []
                prev_end = rpos[n][1]
            else:
                pending.append(n)
        if pending:
            place(pending, prev_end, span_len)
    return core


def read_metrics(g, names, paths, region, region_dir, args):
    """Placement redundancy of the region's existing reads (reads/reads_analysis.py placements)."""
    mcp = os.path.join(region_dir, 'mc.gfa')
    if not os.path.exists(mcp):
        return {'status': 'no_mc_gfa'}
    mc = Graph(mcp)
    mmap, miss = match_paths(mc, names)
    if miss:
        return {'status': 'mc_gfa_missing_paths'}
    mc_paths = {n: mc.paths[mmap[n]] for n in names}
    ref = [n for n in names if n.startswith('CHM13')][0]
    gaf, src = fetch_reads(region, region_dir, mc)
    if gaf is None:
        log('--reads skipped: it needs the big data (the hap32 GBZ-Base and read GAF-Base under '
            'VNTR_EVAL_DIR=%s, or a cached work/cache/evaluate/<id>/reads.gaf)' % config.EVAL_DIR)
        return {'status': 'needs_big_data: ' + src}
    reads = parse_gaf(gaf)
    core = core_nodes(mc, mc_paths, ref, region)
    # the aligned graph sequence of each read that lies wholly inside the baseline subgraph and
    # places >= 1 base on a core node
    rseqs = []
    n_outside = 0
    for r in reads:
        hs = []
        ok = True
        for n, o in r['steps']:
            i = mc.idx.get(n)
            if i is None:
                ok = False
                break
            hs.append(2 * i + (o == '<'))
        if not ok or not hs:
            n_outside += 1
            continue
        cum = 0
        touch = False
        for h in hs:
            L = mc.hlen(h)
            if cum + L > r['ps'] and cum < r['pe'] and (h >> 1) in core:
                touch = True
            cum += L
        if not touch:
            continue
        pseq = mc.spell(hs)[r['ps']:r['pe']]
        if len(pseq) < READ_MIN:
            continue
        rseqs.append((pseq, r['mapq']))
    # placements in the graph under test
    pinfo = []
    for nm in names:
        hs = paths[nm]
        seq = g.spell(hs)
        cum = [0]
        for h in hs:
            cum.append(cum[-1] + g.hlen(h))
        pinfo.append((seq, cum, hs))
    memo = {}

    def placements(seq):
        v = memo.get(seq)
        if v is not None:
            return v
        rc = revcomp(seq)
        L = len(seq)
        places = set()
        max_occ = 0
        nwalk = 0
        for ps, cum, hs in pinfo:
            occ = 0
            for pat, strand in ((seq, 0), (rc, 1)) if rc != seq else ((seq, 0),):
                st = ps.find(pat)
                while st != -1:
                    occ += 1
                    i0 = bisect.bisect_right(cum, st) - 1
                    i1 = bisect.bisect_right(cum, st + L - 1) - 1
                    # a graph placement: its handle path plus the start offset in the first handle
                    np_ = tuple(hs[i0:i1 + 1])
                    s_off = st - cum[i0]
                    e_off = (st + L) - cum[i1]
                    # canonical over the two orientations of the same placement
                    rev = tuple(h ^ 1 for h in reversed(np_))
                    last_len = cum[i1 + 1] - cum[i1]
                    first_len = cum[i0 + 1] - cum[i0]
                    fw = (np_, s_off, e_off)
                    bw = (rev, last_len - e_off, first_len - s_off)
                    places.add(min(fw, bw))
                    st = ps.find(pat, st + 1)
            max_occ = max(max_occ, occ)
            nwalk += occ > 0
        v = (len(places), max_occ, nwalk)
        memo[seq] = v
        return v
    rows = []
    for pseq, mq in rseqs:
        npl, mo, nw = placements(pseq)
        rows.append((npl, mo, mq))

    def split(sub, tag):
        m = len(sub)
        o = {'n': m}
        if not m:
            return {tag + '_' + k: v for k, v in o.items()}
        o['absent'] = round(sum(1 for x in sub if x[0] == 0) / m, 4)
        o['unique'] = round(sum(1 for x in sub if x[0] == 1) / m, 4)
        o['tandem'] = round(sum(1 for x in sub if x[1] >= 2) / m, 4)
        o['redundant_only'] = round(sum(1 for x in sub if x[0] > 1 and x[1] == 1) / m, 4)
        o['redundant'] = round(sum(1 for x in sub if x[0] > max(1, x[1])) / m, 4)
        pr = [x for x in sub if x[0] > 0]
        o['ratio_mean'] = round(sum(x[0] / max(1, x[1]) for x in pr) / len(pr), 4) if pr else None
        o['ratio_median'] = statistics.median([x[0] / max(1, x[1]) for x in pr]) if pr else None
        o['excess_placements_mean'] = round(sum(x[0] - x[1] for x in pr) / len(pr), 4) if pr else None
        return {tag + '_' + k: v for k, v in o.items()}
    out = {'status': 'ok', 'source': src, 'gaf_records': len(reads), 'outside_subgraph': n_outside,
           'core_reads_tested': len(rows)}
    out.update(split(rows, 'all'))
    out.update(split([x for x in rows if x[2] < 5], 'mqlt5'))
    out.update(split([x for x in rows if x[2] >= 30], 'mqge30'))
    return out


# ----------------------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('region_dir')
    ap.add_argument('graph')
    ap.add_argument('--name', default=None, help='method name (default: the graph file\'s directory name)')
    ap.add_argument('--out', default=None, help='output JSON (default results/<METHOD>/<region_id>.json)')
    ap.add_argument('--reads', action='store_true', help='Stage 1: placement redundancy of the existing reads')
    ap.add_argument('--max-pairs', type=int, default=None,
                    help='evaluate all path pairs up to this many, else a seeded sample (default 600); '
                         'with --panel, the size of each seeded sample (default 300)')
    ap.add_argument('--panel', default=None, metavar='FASTA',
                    help='judge a graph whose paths are this panel (e.g. work/panel/mc/<id>.gfa with '
                         'regions/<id>/hprc.fa.gz) instead of hap32.fa: one path per FASTA record, '
                         'same name, spelling it exactly; the CHM13 record is the reference. Size, '
                         'alignment (seeded pair samples), inflation and redundancy only (no truth, '
                         'no reads); output defaults to results/panel/<METHOD>/<region_id>.json')
    ap.add_argument('--affine-budget', type=float, default=5e10,
                    help='affine optima: if the distinct path pairs need more than this many full DP '
                         'cells (len a * len b summed), score the CHM13 pairs plus a seeded sample of '
                         'the others that fits (default 5e10); the choice depends only on the '
                         'sequences, so it is the same for every graph of a region')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--dp-max-bytes', type=float, default=3e9,
                    help='memory cap of the truth-vs-graph DP (default 3e9)')
    ap.add_argument('--no-truvari', action='store_true', help='skip the truvari scoring of the truth')
    ap.add_argument('--no-refine', action='store_true', help='skip truvari refine')
    ap.add_argument('--refine-timeout', type=int, default=900)
    ap.add_argument('--no-cache', action='store_true', help='do not read or write the optima cache')
    ap.add_argument('--skip', default='', help='comma list of sections to skip: alignment,inflation,'
                    'redundancy,truth')
    args = ap.parse_args(argv)
    skip = set(x for x in args.skip.split(',') if x)

    rd = os.path.abspath(args.region_dir)
    region = json.load(open(os.path.join(rd, 'region.json')))
    rid = region.get('region_id') or os.path.basename(rd)
    region['region_id'] = rid
    method = args.name or os.path.basename(os.path.dirname(os.path.abspath(args.graph)))
    if os.path.basename(args.graph) == 'mc.gfa' and not args.name:
        method = 'mc'
    if args.panel:
        return main_panel(args, rd, region, rid, method, skip)
    if args.max_pairs is None:
        args.max_pairs = 600
    outp = args.out or os.path.join(config.RESULTS_DIR, method, rid + '.json')
    work = os.path.join(config.WORK_DIR, 'evaluate', method, rid)
    os.makedirs(work, exist_ok=True)
    res = collections.OrderedDict()
    res['region_id'] = rid
    res['method'] = method
    gp = os.path.abspath(args.graph)
    res['graph'] = os.path.relpath(gp, REPO) if gp.startswith(REPO + os.sep) else gp
    res['evaluate_version'] = VERSION
    for k in ('stratum', 'contig', 'core_start', 'core_end', 'span_start', 'span_end', 'period', 'motif',
              'copies', 'vg_fp', 'vg_fn', 'pg_fp', 'pg_fn', 'in_benchmark'):
        if k in region:
            res[k] = region[k]

    hap = read_fasta(os.path.join(rd, 'hap32.fa'))
    names = [n for n, _ in hap]
    if len(set(names)) != len(names):
        raise SystemExit('duplicate names in hap32.fa')
    refs = [n for n in names if n.startswith('CHM13')]
    if len(refs) != 1:
        raise SystemExit('need exactly one CHM13 sequence in hap32.fa, found %s' % refs)
    ref = refs[0]
    truth = read_fasta(os.path.join(rd, 'truth.fa'))

    g = Graph(args.graph)
    mp, missing = match_paths(g, names)
    mismatch = []
    for n, s in hap:
        if n in mp and g.spell(g.paths[mp[n]]) != s:
            mismatch.append(n)
    extra = [k for k in g.paths if k not in set(mp.values())]
    valid = not missing and not mismatch and not g.bad_steps
    res['validity'] = {'valid': valid, 'paths_expected': len(names), 'paths_missing': missing,
                       'paths_sequence_mismatch': mismatch, 'extra_paths': len(extra),
                       'paths_with_unknown_nodes': g.bad_steps, 'links_to_unknown_nodes': g.bad_links}
    if not valid:
        res['status'] = 'invalid'
        write_out(outp, res)
        log('INVALID graph: missing %d, sequence mismatch %d, unknown nodes %d' %
            (len(missing), len(mismatch), len(g.bad_steps)))
        return 3
    paths = collections.OrderedDict((n, g.paths[mp[n]]) for n in names)
    dag = build_dag(g, list(paths.values()))
    res['validity']['acyclic'] = dag['acyclic']
    span_bp = sum(g.hlen(h) for h in paths[ref])
    log('%s %s: %d paths, %d nodes, span %d bp, acyclic=%s' % (rid, method, len(names), g.n, span_bp,
                                                               dag['acyclic']))
    cache = Cache(rid, enabled=not args.no_cache)
    t = time.time()
    res['size'] = size_metrics(g, paths, ref, span_bp, dag)
    res.setdefault('timing', {})['size'] = round(time.time() - t, 2)
    dref = {}
    if 'alignment' not in skip:
        t = time.time()
        res['alignment'], dref = pair_metrics(g, names, paths, ref, cache, args.max_pairs, args.affine_budget,
                                              args.seed, args.threads)
        res['timing']['alignment'] = round(time.time() - t, 2)
        cache.save()
    if 'inflation' not in skip:
        t = time.time()
        if not dref:
            C = clib()
            rs = g.spell(paths[ref]).encode()
            dref = {k: C.ed(rs, g.spell(paths[n]).encode()) for k, n in enumerate(names) if n != ref}
        res['inflation'] = inflation_metrics(g, names, paths, ref, dref, g.spell(paths[ref]))
        res['timing']['inflation'] = round(time.time() - t, 2)
    if 'redundancy' not in skip:
        t = time.time()
        red = kmer_redundancy(g, paths)
        red['kmer_extra_per_kb'] = round(1000.0 * red['kmer_extra_positions'] / span_bp, 2)
        if 'alignment' in res:
            red['U_per_kb'] = res['alignment'].get('all_U_per_kb')
        res['redundancy'] = red
        res['timing']['redundancy'] = round(time.time() - t, 2)
    if 'truth' not in skip:
        t = time.time()
        if len(truth) >= 2:
            res['truth'] = truth_metrics(g, dag, names, paths, ref, truth, region, cache, work, args)
        else:
            res['truth'] = {'status': 'no_truth_fa'}
        res['timing']['truth'] = round(time.time() - t, 2)
        cache.save()
    if args.reads:
        t = time.time()
        res['reads'] = read_metrics(g, names, paths, region, rd, args)
        res['timing']['reads'] = round(time.time() - t, 2)
    res['status'] = 'ok'
    res['timing']['total'] = round(time.time() - T0, 2)
    write_out(outp, res)
    a = res.get('alignment', {})
    tr = res.get('truth', {})
    log('done: nodes/kb %s, cost/opt %s (ref %s, affine %s), excess/kb %s, SV pieces median %s, '
        'k-mer extra %s, d_graph %s/%s, truvari raw %s refined %s phab %s -> %s' % (
            res['size']['nodes_per_kb'], a.get('all_cost_over_opt'), a.get('ref_cost_over_opt'),
            a.get('affine_all_cost_over_opt'), a.get('all_excess_per_kb'),
            res.get('inflation', {}).get('sv_pieces_median'), res.get('redundancy', {}).get('kmer_frac_extra'),
            tr.get('h1_d_graph'), tr.get('h2_d_graph'), tr.get('truvari', {}).get('raw_f1'),
            tr.get('truvari', {}).get('refined_f1'), tr.get('truvari', {}).get('phab_f1'), outp))
    return 0


PANEL_EXCLUDED = ('HG002', 'HG003', 'HG004', 'NA24385', 'NA24149', 'NA24143', 'gref_CHM13')


def main_panel(args, rd, region, rid, method, skip):
    """--panel FASTA: Stage 0 graph metrics of a graph whose paths are an arbitrary panel (e.g. the
    full Minigraph-Cactus graph over the region, one path per full-panel haplotype). Checks that
    every FASTA record is a path spelling it; no hap32, truth or read checks. Metrics that do not
    depend on how often a walk occurs (size, bubbles, k-mer redundancy) are computed on the
    distinct walks; pair metrics sample haplotype pairs uniformly (so frequent walks weigh more,
    as duplicated hap32 paths do), inflation medians are over haplotypes."""
    if args.max_pairs is None:
        args.max_pairs = 300
    outp = args.out or os.path.join(config.RESULTS_DIR, 'panel', method, rid + '.json')
    res = collections.OrderedDict()
    res['region_id'] = rid
    res['method'] = method
    gp = os.path.abspath(args.graph)
    res['graph'] = os.path.relpath(gp, REPO) if gp.startswith(REPO + os.sep) else gp
    res['evaluate_version'] = VERSION
    for k in ('stratum', 'contig', 'core_start', 'core_end', 'span_start', 'span_end', 'period', 'motif',
              'copies', 'vg_fp', 'vg_fn', 'pg_fp', 'pg_fn', 'in_benchmark'):
        if k in region:
            res[k] = region[k]
    pan = read_fasta(args.panel)
    names = [n for n, _ in pan]
    if not names:
        raise SystemExit('empty panel FASTA %s' % args.panel)
    if len(set(names)) != len(names):
        raise SystemExit('duplicate names in %s' % args.panel)
    bad = [n for n in names if n.split('#', 1)[0] in PANEL_EXCLUDED]
    if bad:
        raise SystemExit('the panel holds excluded haplotypes (HG002 family or gref_CHM13): %s' % bad[:5])
    refs = [n for n in names if n.startswith('CHM13')]
    if len(refs) != 1:
        raise SystemExit('need exactly one CHM13 record in %s, found %s' % (args.panel, refs[:5]))
    ref = refs[0]
    pp = os.path.abspath(args.panel)
    g = Graph(args.graph)
    mp, missing = match_paths(g, names)
    spelled = {}
    mismatch = []
    for n, s in pan:
        if n in mp:
            k = tuple(g.paths[mp[n]])
            if k not in spelled:
                spelled[k] = g.spell(g.paths[mp[n]])
            if spelled[k] != s:
                mismatch.append(n)
    extra = [k for k in g.paths if k not in set(mp.values())]
    valid = not missing and not mismatch and not g.bad_steps
    res['panel'] = {'fasta': os.path.relpath(pp, REPO) if pp.startswith(REPO + os.sep) else pp,
                    'records': len(names), 'reference': ref,
                    'distinct_seqs': len(set(s for _, s in pan))}
    res['validity'] = {'valid': valid, 'paths_expected': len(names), 'paths_missing': missing[:20],
                       'n_paths_missing': len(missing), 'paths_sequence_mismatch': mismatch[:20],
                       'n_paths_sequence_mismatch': len(mismatch), 'extra_paths': len(extra),
                       'paths_with_unknown_nodes': g.bad_steps[:20], 'links_to_unknown_nodes': g.bad_links}
    if not valid:
        res['status'] = 'invalid'
        write_out(outp, res)
        log('INVALID graph: missing %d, sequence mismatch %d, unknown nodes %d' %
            (len(missing), len(mismatch), len(g.bad_steps)))
        return 3
    paths = collections.OrderedDict((n, g.paths[mp[n]]) for n in names)
    distinct = collections.OrderedDict([(ref, paths[ref])])
    seen = {tuple(paths[ref])}
    for n, hs in paths.items():
        if tuple(hs) not in seen:
            seen.add(tuple(hs))
            distinct[n] = hs
    res['panel']['distinct_walks'] = len(distinct)
    dag = build_dag(g, list(distinct.values()))
    res['validity']['acyclic'] = dag['acyclic']
    span_bp = sum(g.hlen(h) for h in paths[ref])
    log('%s %s (panel): %d paths, %d distinct walks, %d nodes, span %d bp, acyclic=%s' % (
        rid, method, len(names), len(distinct), g.n, span_bp, dag['acyclic']))
    cache = Cache(rid, enabled=not args.no_cache)
    t = time.time()
    res['size'] = size_metrics(g, distinct, ref, span_bp, dag)
    res['size']['n_paths'] = len(names)
    res.setdefault('timing', {})['size'] = round(time.time() - t, 2)
    dref = {}
    if 'alignment' not in skip:
        t = time.time()
        res['alignment'], dref = pair_metrics(g, names, paths, ref, cache, args.max_pairs, args.affine_budget,
                                              args.seed, args.threads, panel_sample=True)
        res['timing']['alignment'] = round(time.time() - t, 2)
        cache.save()
    if 'inflation' not in skip:
        t = time.time()
        # anchored/free uses the CHM13 optima of the sampled ref_* pairs only
        res['inflation'] = inflation_metrics(g, names, paths, ref, dref, g.spell(paths[ref]))
        res['inflation']['anchored_over_free_paths'] = len(dref)
        res['timing']['inflation'] = round(time.time() - t, 2)
    if 'redundancy' not in skip:
        t = time.time()
        red = kmer_redundancy(g, distinct)
        red['kmer_extra_per_kb'] = round(1000.0 * red['kmer_extra_positions'] / span_bp, 2)
        if 'alignment' in res:
            red['U_per_kb'] = res['alignment'].get('all_U_per_kb')
        res['redundancy'] = red
        res['timing']['redundancy'] = round(time.time() - t, 2)
    res['status'] = 'ok'
    res['timing']['total'] = round(time.time() - T0, 2)
    write_out(outp, res)
    a = res.get('alignment', {})
    log('done (panel): nodes/kb %s, cost/opt %s (ref %s, affine %s), excess/kb %s, U/kb %s, SV pieces median %s, '
        'k-mer extra %s -> %s' % (res['size']['nodes_per_kb'], a.get('all_cost_over_opt'), a.get('ref_cost_over_opt'),
                                  a.get('affine_all_cost_over_opt'), a.get('all_excess_per_kb'), a.get('all_U_per_kb'),
                                  res.get('inflation', {}).get('sv_pieces_median'),
                                  res.get('redundancy', {}).get('kmer_frac_extra'), outp))
    return 0


def write_out(path, res):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + '.%d.tmp' % os.getpid()
    with open(tmp, 'w') as f:
        json.dump(res, f, indent=1)
    os.replace(tmp, path)


if __name__ == '__main__':
    sys.exit(main())
