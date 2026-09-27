#!/usr/bin/env python3
"""viewer.py -- one self-contained HTML page for a VNTR region: its alignment, its graphs
(the Minigraph-Cactus baseline and any candidate graphs), calls and reads.

    python3 tools/viewer.py TARGET [--candidate GFA [--candidate-label NAME]] ... [--out FILE.html]
                            [--no-data] [--mafft auto|linsi|einsi|ginsi|fftns] [options]

TARGET is one of
    L012184                     a region id: regions/L012184/ when that package exists,
                                otherwise looked up in regions/regions.tsv or the census
    regions/L012184             a region package directory (needs none of the big data)
    chr4:191498564-191501883    CHM13 coordinates, 1-based inclusive (needs the big data)
    chr4 191498564 191501883

--candidate FILE.gfa (repeatable) adds a graph panel for a candidate graph of the same region:
a GFA with one P (or W) line per haplotype sequence of the region, named as in hap32.fa and
spelling exactly those sequences (tools/msa_graph.py writes such files). When FILE.msa.fa exists
next to it (candidates/<method>/L012184.gfa + L012184.msa.fa) the alignment panel can switch to
that MSA; the truth and called rows it lacks are added with mafft --add, which keeps its columns.
The panel is named by --candidate-label, else by the candidate's directory (the method).

Data sources, in order of preference:
  * the region package (region.json, hap32.fa, truth.fa, mc.gfa, calls.tsv);
  * the big data (config.py), when present and --no-data is not given: region.py is run with the
    package's anchors into the work directory, for reads, called haplotypes and truth records;
  * tools/evaluate.py results (results/mc/<id>.json, results/<method>/<id>.json), shown next to
    the viewer's own graph metrics when they exist.
Without the big data (or with --no-data) the page has no reads panel and no called haplotypes;
everything else works from the package alone.

The page: header (locus, TR annotation, one-line summary); alignment (an MSA of CHM13, the
end-to-end panel haplotypes, the HG002 truth haplotypes and each caller's genotype applied to
CHM13: mafft by default, or a candidate's own MSA; colour by differences or by the nodes of any
graph); graphs (a metrics table with one column per graph, a Bandage image per graph coloured by
path sharing, and one node track per graph that follows the MSA columns); calls (vg, PanGenie,
truth records with truvari status); reads (MAPQ, reads per walk, a pileup). Light and dark.
Default output: <work>/viewer/<name>.html. Python standard library only; needs mafft for the
MSA and Bandage for the graph images (both optional: --no-bandage).
"""
import argparse
import base64
import collections
import hashlib
import html
import json
import os
import random
import re
import subprocess
import sys
import time
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
import region as vr  # noqa: E402

GAP_RUN = re.compile(r'-{4,}')
SHARE_BINS = [(1, 1), (2, 3), (4, 8), (9, 16), (17, 10 ** 9)]
SHARE_COLS = ['#d7301f', '#fc8d59', '#e3b448', '#74a9cf', '#0570b0']
REF_GREY = '#9a9a9a'
CLASS_ORDER = {'spanning': 0, 'enters_L': 1, 'exits_R': 2, 'internal': 3}


class ViewError(RuntimeError):
    pass


def log(msg, quiet=False):
    if not quiet:
        sys.stderr.write('[viewer] %s\n' % msg)


# ---------------------------------------------------------------- readers

def read_fasta(path):
    out, name, buf = [], None, []
    if not path or not os.path.exists(path):
        return out
    opener = open
    if path.endswith('.gz'):
        import gzip
        opener = gzip.open
    with opener(path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ''.join(buf)))
                name, buf = line[1:], []
            elif line:
                buf.append(line)
    if name is not None:
        out.append((name, ''.join(buf)))
    return out


def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path) as f:
        lines = f.read().splitlines()
    lines = [ln for ln in lines if ln and not ln.startswith('##')]
    if not lines:
        return []
    head = lines[0].lstrip('#').split('\t')
    return [dict(zip(head, ln.split('\t'))) for ln in lines[1:]]


def read_json(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def num(x, default=None):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return default
    if v != v:
        return default
    return int(v) if v.is_integer() and abs(v) < 2 ** 53 else v


def first_word(name):
    return name.split()[0] if name else name


# ---------------------------------------------------------------- graphs

STEP_W = re.compile(r'([<>])([^<>\s]+)')
STEP_P = re.compile(r'([^,\s]+?)([+-])(?:,|$)')


def flip(h):
    return (h[0], '<' if h[1] == '>' else '>')


def rev_steps(steps):
    return tuple(flip(h) for h in reversed(steps))


def short_name(name):
    """recombination#12#chr4#0 -> r12, CHM13#0#chr4 -> CHM13, GRCh38#0#chr4[113545] -> GRCh38."""
    p = name.split('#')
    if len(p) >= 2:
        s, h = p[0], p[1]
        if s == 'recombination':
            return 'r' + h
        return s if h in ('0', '') else s + '#' + h
    return name


class Graph:
    """A GFA with paths, reduced to distinct walks in CHM13 orientation.

    walks: list of dicts (steps: tuple of (node, '>'|'<'), weight, names, seq, cls, is_ref, id);
    the CHM13 walk comes first, then spanning walks by weight, then fragments."""

    def __init__(self, text, name, kind, src, known_names=None, subpath_names=None):
        self.name, self.kind, self.src = name, kind, src
        self.seqs, self.edges, self.warnings = {}, [], []
        raw = []   # (name, steps, weight, cls_hint, hid)
        toolkit = False
        for line in text.splitlines():
            if not line:
                continue
            f = line.split('\t')
            t = f[0]
            if t == 'S':
                self.seqs[f[1]] = f[2].upper()
            elif t == 'L':
                self.edges.append((f[1], f[2], f[3], f[4]))
            elif t == 'P':
                steps = tuple((n, '>' if o == '+' else '<') for n, o in STEP_P.findall(f[2]))
                raw.append([f[1], steps, 1, None, None])
            elif t == 'W':
                tags = dict((x[:2], x[5:]) for x in f[7:] if len(x) > 5)
                steps = tuple((n, o) for o, n in STEP_W.findall(f[6]))
                wt = int(tags.get('WT', 1))
                if re.match(r'^h\d+$', f[1]) and f[3] in CLASS_ORDER:
                    toolkit = True
                    raw.append([f[1], steps, wt, f[3], f[1]])
                elif f[1] == 'CHM13' and 'WT' in tags and f[3] not in CLASS_ORDER:
                    raw.append(['CHM13', steps, wt, 'spanning', 'CHM13'])
                else:
                    nm = '%s#%s#%s' % (f[1], f[2], f[3])
                    if f[4] not in ('*', '0', ''):
                        nm += '[%s]' % f[4]
                    if known_names is not None and nm not in known_names and nm + '#0' in known_names:
                        nm += '#0'
                    raw.append([nm, steps, wt, None, None])
        self.toolkit = toolkit
        if not raw:
            raise ViewError('%s: no P or W lines' % src)
        missing = set(n for r in raw for n, _ in r[1] if n not in self.seqs)
        if missing:
            raise ViewError('%s: %d path steps use undefined segments' % (src, len(missing)))
        self.node_len = dict((n, len(s)) for n, s in self.seqs.items())
        refs = [r for r in raw if r[0].split('#')[0] == 'CHM13']
        if not refs:
            raise ViewError('%s: no CHM13 path' % src)
        ref = refs[0]
        self.ref_orient = {}
        for n, o in ref[1]:
            self.ref_orient.setdefault(n, o)
        L0, R0 = ref[1][0], ref[1][-1]
        groups = collections.OrderedDict()
        for nm, steps, wt, cls, hid in raw:
            sc = sum(self.node_len[n] * (1 if o == self.ref_orient.get(n, '>') else -1) for n, o in steps)
            if sc < 0:
                steps = rev_steps(steps)
            if cls is None:
                a, b = steps[0], steps[-1]
                if a == L0 and b == R0:
                    cls = 'spanning'
                elif a == L0:
                    cls = 'enters_L'
                elif b == R0:
                    cls = 'exits_R'
                else:
                    cls = 'internal'
            g = groups.setdefault(steps, {'steps': steps, 'weight': 0, 'names': [], 'cls': cls, 'hids': []})
            g['weight'] += wt
            if hid is None:
                g['names'].append(nm)
            else:
                g['hids'].append(hid)
        ref_steps = ref[1] if sum(self.node_len[n] * (1 if o == self.ref_orient.get(n, '>') else -1)
                                  for n, o in ref[1]) >= 0 else rev_steps(ref[1])
        walks = list(groups.values())
        for w in walks:
            w['is_ref'] = 1 if w['steps'] == ref_steps else 0
            w['seq'] = self.walk_seq(w['steps'])
            if toolkit and subpath_names:
                for hid in w['hids']:
                    w['names'].extend(subpath_names.get(hid, []))
            w['names'].sort(key=lambda x: (not x.startswith('CHM13'), not x.startswith('GRCh38'),
                                           [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', x)]))
        walks.sort(key=lambda w: (-w['is_ref'], CLASS_ORDER.get(w['cls'], 4), -w['weight'], -len(w['seq']),
                                  w['names'][:1], w['hids'][:1]))
        for w in walks:
            if toolkit and w['hids']:
                w['id'] = w['hids'][0] if not w['is_ref'] else 'CHM13'
            else:
                w['id'] = short_name(w['names'][0]) if w['names'] else '?'
        seen = collections.Counter(w['id'] for w in walks)
        dup = collections.Counter()
        for w in walks:
            if seen[w['id']] > 1:
                dup[w['id']] += 1
                w['id'] = '%s.%d' % (w['id'], dup[w['id']])
        self.walks = walks
        self.ref_index = next(i for i, w in enumerate(walks) if w['is_ref'])
        self.ref_nodes = set(n for n, _ in walks[self.ref_index]['steps'])
        self.n_paths = len(raw)
        self.weight_total = sum(w['weight'] for w in walks)

    def walk_seq(self, steps):
        return ''.join(self.seqs[n] if o == '>' else vr.revcomp(self.seqs[n]) for n, o in steps)


def load_graph(path, name, kind, known_names=None, subpath_names=None):
    with open(path) as f:
        text = f.read()
    return Graph(text, name, kind, os.path.abspath(path), known_names, subpath_names)


# ---------------------------------------------------------------- graph metrics

def pair_metrics(wa, wb, cache):
    """Graph-implied alignment of two walks against the optimum: the nodes the two walks share
    (in the same orientation and order) are the graph's alignment anchors. Returns
    (cost_graph, unaligned_homology, segments, segments >= 50 bp) or None when the shared
    nodes are not in the same order in both walks."""
    a, b = wa['steps'], wb['steps']
    sb = set(b)
    common = [s for s in a if s in sb]
    ib = {}
    for i, s in enumerate(b):
        ib.setdefault(s, i)
    ia = {}
    for i, s in enumerate(a):
        ia.setdefault(s, i)
    last = -1
    for s in common:
        if ib[s] <= last:
            return None
        last = ib[s]
    if len(set(common)) != len(common):
        return None
    oa, ob = wa['off'], wb['off']
    sa_full, sb_full = wa['seq'], wb['seq']
    cost = U = segs = big = 0
    pa = pb = -1
    for s in common + [None]:
        na, nb = (len(a), len(b)) if s is None else (ia[s], ib[s])
        if na - pa > 1 or nb - pb > 1:
            x = sa_full[oa[pa + 1]:oa[na]]
            y = sb_full[ob[pb + 1]:ob[nb]]
            key = (x, y) if x <= y else (y, x)
            d = cache.get(key)
            if d is None:
                d = vr.edit_distance(x, y)
                cache[key] = d
            cost += d
            U += max(len(x), len(y)) - d
            segs += 1
            if max(len(x), len(y)) >= 50:
                big += 1
        pa, pb = na, nb
    return cost, U, segs, big


def graph_metrics(G, span_len, dopt_cache, max_pairs, seed_key):
    """Metrics for one graph, all from its spanning distinct walks (CHM13 included)."""
    span = [w for w in G.walks if w['cls'] == 'spanning']
    for w in G.walks:
        o, acc = [0], 0
        for n, _ in w['steps']:
            acc += G.node_len[n]
            o.append(acc)
        w['off'] = o
    node_use = collections.Counter()
    for w in G.walks:
        for n in set(x for x, _ in w['steps']):
            node_use[n] += w['weight']
    total_bp = sum(G.node_len.values())
    ref = G.walks[G.ref_index]
    m = collections.OrderedDict()
    m['nodes'] = len(G.seqs)
    m['edges'] = len(G.edges)
    m['node_bp'] = total_bp
    m['nodes_per_kb'] = round(1000.0 * len(G.seqs) / span_len, 1) if span_len else None
    m['mean_node_bp'] = round(total_bp / len(G.seqs), 2) if G.seqs else None
    m['ref_steps'] = len(ref['steps'])
    m['ref_bp_per_node'] = round(len(ref['seq']) / max(1, len(ref['steps'])), 2)
    m['paths'] = G.n_paths
    m['distinct_walks'] = len(G.walks)
    m['spanning_distinct_walks'] = len(span)
    m['spanning_weight'] = sum(w['weight'] for w in span)
    m['weight_total'] = G.weight_total
    m['nonref_nodes'] = sum(1 for n in G.seqs if n not in G.ref_nodes)
    m['nonref_bp'] = sum(G.node_len[n] for n in G.seqs if n not in G.ref_nodes)
    m['private_bp'] = sum(G.node_len[n] for n in G.seqs if node_use[n] <= 1)
    m['directed_cycle'] = bool(vr.graph_has_directed_cycle(
        [n for n in G.seqs], [(a, ao, b, bo) for a, ao, b, bo in G.edges]))
    m['reversing_edges'] = sum(1 for a, ao, b, bo in G.edges if ao != bo)
    revisit = sum(1 for w in G.walks if len(set(n for n, _ in w['steps'])) < len(w['steps']))
    m['walks_revisiting_a_node'] = revisit
    m['kmer'] = kmer_redundancy(G)
    # graph-implied alignment cost vs the optimum
    cache = {}

    def dopt(x, y):
        key = (x['seq'], y['seq']) if x['seq'] <= y['seq'] else (y['seq'], x['seq'])
        d = dopt_cache.get(key)
        if d is None:
            d = vr.edit_distance(x['seq'], y['seq'])
            dopt_cache[key] = d
        return d
    tot = collections.Counter()
    for w in span:
        if w is ref:
            continue
        pm = pair_metrics(ref, w, cache)
        if pm is None:
            tot['unordered'] += 1
            continue
        c, U, segs, big = pm
        d = dopt(ref, w)
        wt = w['weight']
        L = (len(ref['seq']) + len(w['seq'])) / 2.0
        tot['cost'] += wt * c
        tot['dopt'] += wt * d
        tot['U'] += wt * U
        tot['L'] += wt * L
        tot['big'] += wt * big
        tot['w'] += wt
    if tot['w']:
        m['ref_cost_over_opt'] = round(tot['cost'] / tot['dopt'], 3) if tot['dopt'] else 1.0
        m['ref_M_per_kb'] = round(1000.0 * (tot['cost'] - tot['dopt']) / tot['L'], 1)
        m['ref_sv_segments_per_walk'] = round(tot['big'] / tot['w'], 2)
    pairs = [(i, j) for i in range(len(span)) for j in range(i + 1, len(span))]
    sampled = False
    if len(pairs) > max_pairs:
        rng = random.Random(zlib.crc32(seed_key.encode()))
        pairs = rng.sample(pairs, max_pairs)
        sampled = True
    tp = collections.Counter()
    for i, j in pairs:
        pm = pair_metrics(span[i], span[j], cache)
        if pm is None:
            tp['unordered'] += 1
            continue
        c, U, segs, big = pm
        d = dopt(span[i], span[j])
        wt = span[i]['weight'] * span[j]['weight']
        L = (len(span[i]['seq']) + len(span[j]['seq'])) / 2.0
        tp['cost'] += wt * c
        tp['dopt'] += wt * d
        tp['U'] += wt * U
        tp['L'] += wt * L
        tp['w'] += wt
        tp['pairs'] += 1
    if tp['w']:
        m['cost_over_opt'] = round(tp['cost'] / tp['dopt'], 3) if tp['dopt'] else 1.0
        m['M_per_kb'] = round(1000.0 * (tp['cost'] - tp['dopt']) / tp['L'], 1)
        m['U_per_kb'] = round(1000.0 * tp['U'] / tp['L'], 1)
        m['dopt_per_kb'] = round(1000.0 * tp['dopt'] / tp['L'], 1)
    m['pairs'] = tp['pairs']
    m['pairs_sampled'] = sampled
    m['pairs_unordered'] = tp['unordered'] + tot['unordered']
    return m, node_use


def kmer_redundancy(G, k=21):
    """MSA-independent measure of sequence the graph keeps in parallel copies. Every k-mer
    spelled by a walk has a graph position (node, offset, orientation of its first base).
    P(K) = distinct positions of K over all walks; C(K) = most copies of K in any one walk.
    A graph that merges shared copies needs only C(K) positions: extra = sum(P - C) / sum(P).
    Reordered repeat copies also add to it, so read it as an upper bound on unmerged sequence."""
    pos = collections.defaultdict(set)
    maxocc = collections.Counter()
    for w in G.walks:
        seq = w['seq']
        bp = []
        for n, o in w['steps']:
            L = G.node_len[n]
            bp.extend((n, j, o) for j in range(L))
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
    return {'k': k, 'positions': P, 'extra_positions': extra, 'frac_extra': round(extra / P, 4) if P else 0.0,
            'kmers': len(pos), 'kmers_with_extra': multi}


# ---------------------------------------------------------------- Bandage

def share_colour(w):
    for (lo, hi), c in zip(SHARE_BINS, SHARE_COLS):
        if lo <= w <= hi:
            return c
    return SHARE_COLS[-1]


def bandage_png(G, node_use, workdir, quiet, height=1400):
    if not os.path.exists(config.BANDAGE):
        return None, 'Bandage not found'
    lines = ['H\tVN:Z:1.0']
    for n in sorted(G.seqs, key=lambda x: (len(x), x)):
        col = REF_GREY if n in G.ref_nodes else share_colour(node_use.get(n, 0))
        lines.append('S\t%s\t%s\tCL:z:%s\tDP:f:%d' % (n, G.seqs[n], col, max(1, node_use.get(n, 1))))
    for a, ao, b, bo in G.edges:
        lines.append('L\t%s\t%s\t%s\t%s\t0M' % (a, ao, b, bo))
    text = '\n'.join(lines) + '\n'
    total = sum(len(s) for s in G.seqs.values())
    nodewidth = max(5.0, min(200.0, total / 150.0))
    key = hashlib.sha1((text + str(nodewidth) + str(height)).encode()).hexdigest()[:16]
    gfa = os.path.join(workdir, 'bandage.%s.gfa' % key)
    png = os.path.join(workdir, 'bandage.%s.png' % key)
    if not os.path.exists(png):
        with open(gfa, 'w') as f:
            f.write(text)
        t0 = time.time()
        cmd = [config.BANDAGE, 'image', gfa, png, '--height', str(height), '--colour', 'custom',
               '--nodewidth', '%.1f' % nodewidth, '--depwidth', '0.5', '--edgewidth', '%.1f' % max(1.5, nodewidth / 15)]
        try:
            subprocess.run(cmd, capture_output=True, text=True, timeout=900, check=True)
        except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            return None, 'Bandage failed: %s' % str(e)[:200]
        log('Bandage %s done in %.1f s (%d nodes)' % (G.name, time.time() - t0, len(G.seqs)), quiet)
    if not os.path.exists(png):
        return None, 'Bandage produced no image'
    with open(png, 'rb') as f:
        return f.read(), None


# ---------------------------------------------------------------- MSA

def cigar_ref_to_query(cigar):
    """Map each reference offset to the query offset aligned at or after it (list, len ref+1)."""
    m, qi = [], 0
    for n, op in re.findall(r'(\d+)([=XIDM])', cigar):
        n = int(n)
        if op in '=XM':
            for _ in range(n):
                m.append(qi)
                qi += 1
        elif op == 'I':
            qi += n
        elif op == 'D':
            m.extend([qi] * n)
    m.append(qi)
    return m


def mafft_args(mode):
    return {'linsi': ['--localpair', '--maxiterate', '1000'],
            'einsi': ['--genafpair', '--maxiterate', '1000'],
            'ginsi': ['--globalpair', '--maxiterate', '1000'],
            'fftns': ['--retree', '2', '--maxiterate', '0']}.get(mode, ['--auto'])


def run_mafft(seqs, workdir, mode, threads, quiet):
    """seqs: unique non-empty sequences. Returns their aligned strings (upper case), cached."""
    blob = ('%s\n' % mode + '\n'.join(seqs)).encode()
    key = hashlib.sha1(blob).hexdigest()[:16]
    outp = os.path.join(workdir, 'msa.%s.fa' % key)
    if not os.path.exists(outp):
        inp = os.path.join(workdir, 'msa.%s.in.fa' % key)
        with open(inp, 'w') as f:
            for i, s in enumerate(seqs):
                f.write('>s%d\n%s\n' % (i, s))
        cmd = [config.MAFFT] + mafft_args(mode) + ['--thread', str(threads), '--quiet', inp]
        t0 = time.time()
        log('mafft %s on %d sequences (%d bp)' % (mode, len(seqs), sum(map(len, seqs))), quiet)
        with open(outp + '.tmp', 'w') as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=True, env=config.tool_env())
        os.replace(outp + '.tmp', outp)
        log('mafft done in %.1f s' % (time.time() - t0), quiet)
    aln = dict((first_word(n), s.upper()) for n, s in read_fasta(outp))
    out = [aln['s%d' % i] for i in range(len(seqs))]
    for s, a in zip(seqs, out):
        if a.replace('-', '') != s.upper():
            raise ViewError('mafft output does not match its input')
    return out, key


def strip_gap_columns(strs):
    if not strs:
        return strs
    keep = [c for c in range(len(strs[0])) if any(s[c] != '-' for s in strs)]
    return [''.join(s[c] for c in keep) for s in strs]


def mafft_add(existing, new, workdir, threads, quiet):
    """Add new sequences to a fixed alignment (mafft --add). Returns (existing', new') aligned
    strings, or (None, reason) if the existing alignment did not survive unchanged."""
    blob = ('add\n' + '\n'.join(existing) + '\n#\n' + '\n'.join(new)).encode()
    key = hashlib.sha1(blob).hexdigest()[:16]
    outp = os.path.join(workdir, 'add.%s.fa' % key)
    if not os.path.exists(outp):
        ex = os.path.join(workdir, 'add.%s.existing.fa' % key)
        nw = os.path.join(workdir, 'add.%s.new.fa' % key)
        with open(ex, 'w') as f:
            for i, s in enumerate(existing):
                f.write('>e%d\n%s\n' % (i, s))
        with open(nw, 'w') as f:
            for i, s in enumerate(new):
                f.write('>n%d\n%s\n' % (i, s))
        t0 = time.time()
        log('mafft --add of %d sequences to a %d-row MSA' % (len(new), len(existing)), quiet)
        with open(outp + '.tmp', 'w') as out:
            subprocess.run([config.MAFFT, '--add', nw, '--thread', str(threads), '--quiet', ex],
                           stdout=out, stderr=subprocess.PIPE, check=True, env=config.tool_env())
        os.replace(outp + '.tmp', outp)
        log('mafft --add done in %.1f s' % (time.time() - t0), quiet)
    aln = dict((first_word(n), s.upper()) for n, s in read_fasta(outp))
    e2 = [aln.get('e%d' % i) for i in range(len(existing))]
    n2 = [aln.get('n%d' % i) for i in range(len(new))]
    if any(x is None for x in e2 + n2):
        return None, 'mafft --add lost sequences'
    for s, a in zip(new, n2):
        if a.replace('-', '') != s.upper():
            return None, 'mafft --add changed an added sequence'
    if strip_gap_columns(e2) != strip_gap_columns([s.upper() for s in existing]):
        return None, 'mafft --add changed the existing alignment'
    return e2, n2


# ---------------------------------------------------------------- reads

def parse_gaf(path, max_reads):
    reads = []
    if not path or not os.path.exists(path):
        return reads, 0
    n_total = 0
    with open(path) as f:
        for line in f:
            if line.startswith('@') or not line.strip():
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 12:
                continue
            n_total += 1
            if len(reads) >= max_reads:
                continue
            steps = [(n, o) for o, n in STEP_W.findall(p[5])]
            if not steps:
                continue
            tags = dict((t[:2], t[5:]) for t in p[12:] if len(t) > 5 and t[:2] in ('AS', 'pd'))
            reads.append({'name': p[0], 'qlen': int(p[1]), 'qs': int(p[2]), 'qe': int(p[3]), 'strand': p[4],
                          'steps': steps, 'plen': int(p[6]), 'ps': int(p[7]), 'pe': int(p[8]),
                          'matches': int(p[9]), 'blen': int(p[10]), 'mapq': int(p[11]),
                          'AS': num(tags.get('AS')), 'pd': tags.get('pd')})
    return reads, n_total


def place_reads(reads, G, draw_walks):
    """For every read: the distinct walks of G that contain its in-subgraph node path
    contiguously (either orientation), and a placement (walk, base start, base end, partial)
    on one of draw_walks (the walks that have an MSA row)."""
    node_len = G.node_len
    walks = G.walks
    index, offs = [], []
    for w in walks:
        d = collections.defaultdict(list)
        for i, h in enumerate(w['steps']):
            d[h].append(i)
        index.append(d)
        o, acc = [], 0
        for h in w['steps']:
            o.append(acc)
            acc += node_len[h[0]]
        o.append(acc)
        offs.append(o)
    draw = set(draw_walks)
    out = []
    for rd in reads:
        st = rd['steps']
        runs, cur = [], []
        for k, h in enumerate(st):
            if h[0] in node_len:
                cur.append(k)
            elif cur:
                runs.append(cur)
                cur = []
        if cur:
            runs.append(cur)
        rec = {'status': 'outside', 'compat': [], 'place': None, 'partial': 0, 'ori': 1, 'inside_all': 0}
        if runs:
            best = max(runs, key=lambda r: sum(node_len[st[k][0]] for k in r))
            i0, i1 = best[0], best[-1] + 1
            run = tuple(st[i0:i1])
            rlen = sum(node_len[h[0]] for h in run)
            lead_out, trail_out = i0 > 0, i1 < len(st)
            if not lead_out and not trail_out:
                s, e = rd['ps'], rd['pe']
            elif lead_out and not trail_out:
                lead = rd['plen'] - rlen
                s, e = max(0, rd['ps'] - lead), rd['pe'] - lead
            elif trail_out and not lead_out:
                s, e = rd['ps'], rlen
            else:
                s, e = 0, rlen
            s, e = max(0, min(s, rlen)), max(0, min(e, rlen))
            orients = [(run, s, e, 1), (rev_steps(run), rlen - e, rlen - s, -1)]
            rec['inside_all'] = int(len(runs) == 1 and not lead_out and not trail_out)
            full = []
            for wi, w in enumerate(walks):
                hit = None
                for seq, ss, ee, ori in orients:
                    for i in index[wi].get(seq[0], ()):
                        if w['steps'][i:i + len(seq)] == seq:
                            hit = (i, ss, ee, ori)
                            break
                    if hit:
                        break
                if hit:
                    full.append((wi, hit))
            rec['compat'] = [wi for wi, _ in full]
            on_draw = [(wi, h) for wi, h in full if wi in draw]
            if on_draw:
                rec['status'] = 'walk'
                wi, (i, ss, ee, ori) = max(on_draw, key=lambda t: (walks[t[0]]['weight'], -t[0]))
                rec['place'] = (wi, offs[wi][i] + ss, offs[wi][i] + ee)
                rec['ori'] = ori
            else:
                rec['status'] = 'fragment' if full else 'novel'
                bestp = None
                for wi in draw_walks:
                    w = walks[wi]
                    for seq, ss, ee, ori in orients:
                        L = len(seq)
                        for j in range(L):
                            for i in index[wi].get(seq[j], ()):
                                k = 0
                                while j + k < L and i + k < len(w['steps']) and w['steps'][i + k] == seq[j + k]:
                                    k += 1
                                bp = offs[wi][i + k] - offs[wi][i]
                                if bestp is None or bp > bestp[0]:
                                    bestp = (bp, wi, i, j, k, ss, ee, L, ori, seq)
                if bestp:
                    bp, wi, i, j, k, ss, ee, L, ori, seq = bestp
                    b0 = offs[wi][i] + (ss if j == 0 else 0)
                    b1 = offs[wi][i + k]
                    if j + k == L:
                        b1 = offs[wi][i] + (ee - sum(node_len[h[0]] for h in seq[:j]))
                    if b1 > b0:
                        rec['place'] = (wi, b0, b1)
                        rec['partial'] = 1
                        rec['ori'] = ori
        out.append(rec)
    return out


# ---------------------------------------------------------------- census / annotation

def awk_overlap(path, contig, a0, b0, ccol, scol, ecol):
    """Rows of a headered TSV on <contig> overlapping [a0,b0) (0-based half-open)."""
    if not path or not os.path.exists(path):
        return [], []
    with open(path) as f:
        head = f.readline().rstrip('\n').split('\t')
    prog = 'NR>1 && $%d==c && $%d<b && $%d>a' % (ccol + 1, scol + 1, ecol + 1)
    try:
        out = subprocess.run(['awk', '-F\t', '-v', 'c=' + contig, '-v', 'a=%d' % a0, '-v', 'b=%d' % b0, prog, path],
                             capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return head, []
    return head, [dict(zip(head, ln.split('\t'))) for ln in out.splitlines() if ln]


def census_lookup(contig, req0, req1, span0, span1):
    res = {'locus': None, 'loci': [], 'strata': [], 'regions': []}
    lp = config.census_file('loci.tsv')
    _, loci = awk_overlap(lp, contig, span0, span1, 1, 2, 3)
    for r in loci:
        r['_ov'] = max(0, min(int(r['end']), req1) - max(int(r['start']), req0))
    loci.sort(key=lambda r: -r['_ov'])
    res['loci'] = loci
    if loci and loci[0]['_ov'] > 0:
        res['locus'] = loci[0]
    for r in read_tsv(config.census_file('strata.tsv')):
        if r.get('contig') != contig:
            continue
        s, e = int(r['start']), int(r['end'])
        if (res['locus'] and r.get('locus_id') == res['locus']['locus_id']) or (s < req1 and e > req0):
            res['strata'].append(r)
    _, regs = awk_overlap(config.census_file('vntr_regions.tsv'), contig, span0, span1, 1, 2, 3)
    _, cat = awk_overlap(config.census_file('catalogue.tsv'), contig, span0, span1, 1, 2, 3)
    seen = set(r['region_id'] for r in regs)
    regs += [r for r in cat if r['region_id'] not in seen]
    regs.sort(key=lambda r: int(r['start']))
    res['regions'] = regs
    return res


def resolve_target(t):
    """Returns ('package', dir) or ('coords', contig, start1, end1, region_id or None)."""
    if os.path.isdir(t) and os.path.exists(os.path.join(t, 'region.json')):
        return ('package', os.path.abspath(t))
    m = re.match(r'^([^:\s]+):([\d,]+)-([\d,]+)$', t)
    if m:
        return ('coords', m.group(1), int(m.group(2).replace(',', '')), int(m.group(3).replace(',', '')), None)
    pdir = os.path.join(config.REGIONS_DIR, t)
    if os.path.exists(os.path.join(pdir, 'region.json')):
        return ('package', pdir)
    for r in read_tsv(os.path.join(config.REGIONS_DIR, 'regions.tsv')):
        if r.get('region_id') == t:
            return ('coords', r['contig'], int(r['core_start']), int(r['core_end']), t)
    if re.match(r'^L\d+$', t):
        for fn in ('strata.tsv', 'loci.tsv'):
            for r in read_tsv(config.census_file(fn)):
                if r.get('locus_id') == t:
                    # census tables are 0-based half-open
                    return ('coords', r['contig'], int(r['start']) + 1, int(r['end']), t)
        raise ViewError('region %s: no package in %s and not in the census (%s)' % (t, config.REGIONS_DIR,
                                                                                    config.CENSUS_DIR))
    if re.match(r'^TR\d+$', t):
        for fn in ('vntr_regions.tsv', 'catalogue.tsv'):
            path = config.census_file(fn)
            if not path:
                continue
            out = subprocess.run(['awk', '-F\t', '-v', 'id=' + t, '$1==id{print $2"\t"$3"\t"$4; exit}', path],
                                 capture_output=True, text=True).stdout.split()
            if out:
                return ('coords', out[0], int(out[1]) + 1, int(out[2]), t)
        raise ViewError('region %s not in the census catalogue' % t)
    raise ViewError('give a region id, a region package directory, chr:start-end or <contig> <start> <end>')


# ---------------------------------------------------------------- package / toolkit sources

class Sources:
    """Everything the page needs, from a region package and/or a region.py output directory."""

    def __init__(self):
        self.pkg = None          # package dir
        self.rj = {}             # region.json
        self.tk = None           # toolkit dir
        self.S = {}              # toolkit summary.json
        self.warnings = []


def run_toolkit(contig, start, end, outdir, anchors, args, quiet):
    if args.rerun or not os.path.exists(os.path.join(outdir, 'summary.json')):
        log('running region.py into %s' % outdir, quiet)
        vr.analyse(contig, start, end, outdir, pad=args.pad, reads=not args.no_reads, names=args.names,
                   verbose=not quiet, anchors=anchors)
    return outdir


def gather(args):
    src = Sources()
    quiet = args.quiet
    tgt = resolve_target(args.target) if args.start is None else \
        ('coords', args.target, args.start, args.end, None)
    use_data = config.have_graph_data() and not args.no_data
    if tgt[0] == 'package':
        src.pkg = tgt[1]
        src.rj = read_json(os.path.join(src.pkg, 'region.json'))
        for k in ('contig', 'span_start', 'span_end'):
            if k not in src.rj:
                raise ViewError('%s/region.json lacks %s' % (src.pkg, k))
        rid = src.rj.get('region_id') or os.path.basename(src.pkg.rstrip('/'))
        name = rid
        contig = src.rj['contig']
        core = [int(src.rj.get('core_start') or src.rj['span_start']), int(src.rj.get('core_end') or src.rj['span_end'])]
        if args.dir:
            src.tk = os.path.abspath(args.dir)
        elif use_data:
            anchors = (src.rj.get('anchor_left'), src.rj.get('anchor_right'))
            if not all(anchors):
                anchors = None
                src.warnings.append('region.json has no anchors; region.py chose its own, which may '
                                    'give a different span')
            try:
                src.tk = run_toolkit(contig, core[0], core[1], os.path.join(args.workroot, name + '.region'),
                                     anchors, args, quiet)
            except vr.ToolError as e:
                src.warnings.append('region.py failed, showing the package only: %s' % str(e)[:300])
                src.tk = None
    else:
        _, contig, s1, e1, rid = tgt
        name = rid or '%s_%d_%d' % (contig, s1, e1)
        core = [s1, e1]
        if args.dir:
            src.tk = os.path.abspath(args.dir)
        else:
            if not use_data:
                raise ViewError('%s needs the big data (%s) or a region package' % (args.target, config.EVAL_DIR))
            src.tk = run_toolkit(contig, s1, e1, os.path.join(args.workroot, name + '.region'), None, args, quiet)
    if src.tk:
        if not os.path.exists(os.path.join(src.tk, 'summary.json')):
            raise ViewError('%s has no summary.json' % src.tk)
        src.S = read_json(os.path.join(src.tk, 'summary.json'))
        if not src.pkg:
            core = list(src.S['requested_interval_1based'])
    src.name, src.contig, src.core, src.region_id = name, contig, core, rid
    return src


# ---------------------------------------------------------------- main build

def build(args):
    t_start = time.time()
    quiet = args.quiet
    os.makedirs(args.workroot, exist_ok=True)
    src = gather(args)
    name, contig = src.name, src.contig
    req1, req2 = src.core
    vdir = os.path.join(args.workroot, name + '.view')
    os.makedirs(vdir, exist_ok=True)
    warnings = list(src.warnings)
    S = src.S
    tk = src.tk

    # ---- span, reference, baseline graph
    hap32 = read_fasta(os.path.join(src.pkg, 'hap32.fa')) if src.pkg else []
    known = set(first_word(n) for n, _ in hap32) if hap32 else None
    hap_seq_by_name = dict((first_word(n), s.upper()) for n, s in hap32)
    sub_names = collections.defaultdict(list)
    tk_graph = None
    if tk:
        for sp in read_tsv(os.path.join(tk, 'subpaths.tsv')):
            if sp.get('path_name') not in (None, '.', ''):
                sub_names[sp['distinct_id']].append(sp['path_name'])
        tk_graph = load_graph(os.path.join(tk, 'subgraph.gfa'), 'toolkit', 'baseline', None, sub_names)
    if src.pkg and os.path.exists(os.path.join(src.pkg, 'mc.gfa')):
        base = load_graph(os.path.join(src.pkg, 'mc.gfa'), 'Minigraph-Cactus (hap32)', 'baseline', known)
        A1, B1 = int(src.rj['span_start']), int(src.rj['span_end'])
    elif tk_graph:
        base = tk_graph
        base.name = 'Minigraph-Cactus (hap32)'
        A1, B1 = S['analysed_span_1based']
    else:
        raise ViewError('no baseline graph: the package has no mc.gfa and region.py was not run')
    refseq = base.walks[base.ref_index]['seq']
    if B1 - A1 + 1 != len(refseq):
        warnings.append('CHM13 walk is %d bp but the span %s:%d-%d is %d bp' % (len(refseq), contig, A1, B1, B1 - A1 + 1))
    if tk and src.pkg:
        tA, tB = S.get('analysed_span_1based', [None, None])
        if [tA, tB] != [A1, B1]:
            warnings.append('region.py span %s:%s-%s differs from the package span %d-%d; its reads, '
                            'calls and called haplotypes are not used' % (contig, tA, tB, A1, B1))
            tk, S, tk_graph = None, {}, None
    if hap32:
        bad = 0
        for w in base.walks:
            for nm in w['names']:
                if nm in hap_seq_by_name and hap_seq_by_name[nm] != w['seq']:
                    bad += 1
        if bad:
            warnings.append('%d mc.gfa paths do not spell their hap32.fa sequence' % bad)

    # ---- truth, called haplotypes
    truth = []
    if src.pkg and os.path.exists(os.path.join(src.pkg, 'truth.fa')):
        truth = [s.upper() for _, s in read_fasta(os.path.join(src.pkg, 'truth.fa'))][:2]
    elif tk:
        truth = [s.upper() for f in ('truth.hap1.fa', 'truth.hap2.fa')
                 for _, s in read_fasta(os.path.join(tk, f))[:1]]
    called = []
    if tk:
        for lab, fn in (('vg', 'vg.called.fa'), ('PanGenie', 'pangenie.called.fa')):
            for k, (_, s) in enumerate(read_fasta(os.path.join(tk, fn))):
                called.append((lab, k + 1, s.upper()))

    # ---- per-walk distances (baseline), closest walk to each truth haplotype
    span_walks = [i for i, w in enumerate(base.walks) if w['cls'] == 'spanning']
    wmeta = {}
    for wi in span_walks:
        w = base.walks[wi]
        s = w['seq']
        e = {'edit_ref': vr.edit_distance(s, refseq)}
        e['id_ref'] = round(vr.identity(e['edit_ref'], len(s), len(refseq)), 4)
        for h, t in enumerate(truth):
            d = vr.edit_distance(s, t)
            e['edit_h%d' % (h + 1)] = d
            e['id_h%d' % (h + 1)] = round(vr.identity(d, len(s), len(t)), 4)
        e['nonref_nodes'] = len(set(n for n, _ in w['steps'] if n not in base.ref_nodes))
        e['nonref_bp'] = sum(base.node_len[n] for n, _ in w['steps'] if n not in base.ref_nodes)
        wmeta[wi] = e
    closest = {}
    for h, t in enumerate(truth):
        k = 'edit_h%d' % (h + 1)
        ranked = sorted(span_walks, key=lambda wi: (wmeta[wi][k], abs(len(base.walks[wi]['seq']) - len(t))))
        nonref = [wi for wi in ranked if not base.walks[wi]['is_ref']]
        closest['h%d' % (h + 1)] = {
            'walk': base.walks[ranked[0]]['id'], 'edit': wmeta[ranked[0]][k], 'identity': wmeta[ranked[0]]['id_h%d' % (h + 1)],
            'walk_nonref': base.walks[nonref[0]]['id'] if nonref else None,
            'edit_nonref': wmeta[nonref[0]][k] if nonref else None,
            'edit_ref': vr.edit_distance(refseq, t), 'len': len(t),
            'exact_weight': sum(base.walks[wi]['weight'] for wi in ranked if wmeta[wi][k] == 0)}
    by_seq = collections.defaultdict(list)
    for wi in span_walks:
        by_seq[base.walks[wi]['seq']].append(base.walks[wi]['id'])
    same_seq_groups = [g for g in by_seq.values() if len(g) > 1]
    same_seq_of = dict((hid, [x for x in g if x != hid]) for g in same_seq_groups for hid in g)

    # ---- rows and sequence ids
    sid_of, seqs = {}, []

    def sid(s):
        if s not in sid_of:
            sid_of[s] = len(seqs)
            seqs.append(s)
        return sid_of[s]
    rows = []
    rw = base.walks[base.ref_index]
    rows.append({'kind': 'ref', 'label': 'CHM13', 'sid': sid(refseq), 'walk': base.ref_index, 'weight': rw['weight']})
    for h, t in enumerate(truth):
        c = closest.get('h%d' % (h + 1), {})
        rows.append({'kind': 'truth', 'label': 'HG002 h%d' % (h + 1), 'sid': sid(t), 'closest': c.get('walk'),
                     'closest_edit': c.get('edit'), 'closest_id': c.get('identity'),
                     'closest_nonref': c.get('walk_nonref'), 'closest_nonref_edit': c.get('edit_nonref'),
                     'edit_ref': c.get('edit_ref')})
    chv = S.get('called_haplotypes_vs_truth', {}) if S else {}
    for lab, slot, s in called:
        cv = chv.get('vg' if lab == 'vg' else 'pangenie', {})
        m = re.search(r'slot%d-(h\d)' % slot, cv.get('pairing', ''))
        r = {'kind': 'called', 'label': '%s slot%d' % (lab, slot), 'sid': sid(s), 'caller': lab}
        if m:
            r.update({'paired_truth': m.group(1), 'edit_truth': cv.get('d_' + m.group(1)),
                      'phase_reliable': cv.get('phase_reliable')})
        rows.append(r)
    comp_ids = collections.defaultdict(list)
    for k, c in closest.items():
        comp_ids[c['walk']].append(k)
    for wi in span_walks:
        w = base.walks[wi]
        if w['is_ref']:
            rows[0].update(wmeta[wi])
            rows[0]['names'] = w['names'][:60]
            continue
        r = {'kind': 'panel', 'label': w['id'], 'sid': sid(w['seq']), 'walk': wi, 'weight': w['weight'],
             'names': w['names'][:60], 'closest_to': comp_ids.get(w['id'], []),
             'same_seq': same_seq_of.get(w['id'], []), 'steps': len(w['steps'])}
        r.update(wmeta[wi])
        rows.append(r)
    for r in rows:
        r['len'] = len(seqs[r['sid']])

    # ---- MSAs
    msas = []
    msas.append(build_mafft_msa(rows, seqs, refseq, base, A1, req1, req2, args, vdir, warnings, S, tk))
    # ---- candidate graphs
    graphs = [base]
    for k, cpath in enumerate(args.candidate or []):
        label = (args.candidate_label[k] if args.candidate_label and k < len(args.candidate_label)
                 else os.path.basename(os.path.dirname(os.path.abspath(cpath))) or 'candidate')
        try:
            G = load_graph(cpath, label, 'candidate', set(hap_seq_by_name) if hap_seq_by_name else None)
        except (ViewError, OSError) as e:
            warnings.append('candidate %s not loaded: %s' % (cpath, e))
            continue
        # every candidate path must spell a panel sequence
        panel_seqs = set(seqs[r['sid']] for r in rows if r['kind'] in ('ref', 'panel'))
        nbad = sum(1 for w in G.walks if w['seq'] not in panel_seqs)
        if hap_seq_by_name:
            wrong = [nm for w in G.walks for nm in w['names']
                     if nm in hap_seq_by_name and hap_seq_by_name[nm] != w['seq']]
            unknown = [nm for w in G.walks for nm in w['names'] if nm not in hap_seq_by_name]
            if wrong:
                G.warnings.append('%d paths do not spell their hap32.fa sequence (%s)' % (len(wrong), ', '.join(wrong[:4])))
            if unknown:
                G.warnings.append('%d path names are not in hap32.fa (%s)' % (len(unknown), ', '.join(unknown[:4])))
            missing = set(hap_seq_by_name) - set(nm for w in G.walks for nm in w['names'])
            if missing:
                G.warnings.append('%d hap32.fa sequences have no path (%s)' % (len(missing), ', '.join(sorted(missing)[:4])))
        if nbad:
            G.warnings.append('%d distinct walks spell no panel sequence; they have no MSA row' % nbad)
        if G.walks[G.ref_index]['seq'] != refseq:
            G.warnings.append('the CHM13 path does not spell the CHM13 span')
        graphs.append(G)
        stem = os.path.splitext(os.path.abspath(cpath))[0]
        mpath = stem + '.msa.fa'
        if os.path.exists(mpath) and not args.no_candidate_msa:
            try:
                msas.append(build_file_msa(mpath, label, rows, seqs, A1, req1, req2, args, vdir, warnings))
            except (ViewError, subprocess.CalledProcessError) as e:
                warnings.append('MSA %s not used: %s' % (mpath, e))

    # ---- graph metrics, node tables, Bandage
    dopt_cache = {}
    graphs_out = []
    span_len = len(refseq)
    for gi, G in enumerate(graphs):
        t0 = time.time()
        met, node_use = graph_metrics(G, span_len, dopt_cache, args.max_pairs, name + G.name)
        log('metrics for %s in %.1f s' % (G.name, time.time() - t0), quiet)
        nidx = {}
        for w in G.walks:
            for n, _ in w['steps']:
                if n not in nidx:
                    nidx[n] = len(nidx)
        for n in G.seqs:
            if n not in nidx:
                nidx[n] = len(nidx)
        nodes_list = sorted(nidx, key=lambda n: nidx[n])
        nw = collections.Counter()
        for w in G.walks:
            for n in set(x for x, _ in w['steps']):
                nw[n] += 1
        # which walk of this graph gives each row its nodes (ref and panel rows only)
        by_seq_g = {}
        for wi, w in enumerate(G.walks):
            if w['cls'] != 'spanning':
                continue
            s = w['seq']
            if s not in by_seq_g or G.walks[by_seq_g[s]]['weight'] < w['weight']:
                by_seq_g[s] = wi
        row_walk = []
        for r in rows:
            if r['kind'] not in ('ref', 'panel'):
                row_walk.append(-1)
            elif gi == 0:
                row_walk.append(r['walk'])
            else:
                row_walk.append(by_seq_g.get(seqs[r['sid']], -1))
        png, note = (None, 'Bandage disabled') if args.no_bandage else bandage_png(G, node_use, vdir, quiet)
        if note and not args.no_bandage:
            warnings.append('%s: %s' % (G.name, note))
        graphs_out.append({
            'name': G.name, 'kind': G.kind, 'src': G.src, 'metrics': met, 'warnings': G.warnings,
            'nodes': {'id': nodes_list, 'len': [G.node_len[n] for n in nodes_list],
                      'w': [node_use[n] for n in nodes_list], 'nw': [nw[n] for n in nodes_list],
                      'ref': [1 if n in G.ref_nodes else 0 for n in nodes_list]},
            'walks': [{'id': w['id'], 'weight': w['weight'], 'cls': w['cls'], 'len': len(w['seq']),
                       'is_ref': w['is_ref'], 'sid': sid_of.get(w['seq'], -1), 'names': w['names'][:60],
                       'runs': encode_runs(w['steps'], nidx)} for w in G.walks],
            'row_walk': row_walk, 'png': png, 'bandage_note': note,
            'evaluate': evaluate_results(src.region_id, 'mc' if gi == 0 else method_of(G.src), G.src)})

    # ---- reads (placed on the baseline walks)
    reads, n_reads_total, placed = [], 0, []
    if tk and not args.no_reads:
        reads, n_reads_total = parse_gaf(os.path.join(tk, 'reads.gaf'), args.max_reads)
        if n_reads_total > len(reads):
            warnings.append('%d of %d read alignments kept (--max-reads)' % (len(reads), n_reads_total))
        t0 = time.time()
        draw = [i for i, w in enumerate(base.walks) if w['cls'] == 'spanning']
        rn = set(n for rd in reads for n, _ in rd['steps'] if tk_graph and n in tk_graph.node_len)
        overlap = (sum(1 for n in rn if n in base.node_len) / len(rn)) if rn else 1.0
        if overlap >= 0.99 or not tk_graph:
            placed = place_reads(reads, base, draw)
        else:
            # the package graph has other node ids: place on region.py's graph, map walks by sequence
            warnings.append('mc.gfa node ids differ from the read alignments (%.0f%% shared); reads were placed on '
                            "region.py's subgraph and mapped to baseline walks by sequence" % (100 * overlap))
            tdraw = [i for i, w in enumerate(tk_graph.walks) if w['cls'] == 'spanning']
            placed = place_reads(reads, tk_graph, tdraw)
            bseq = {}
            for i in draw:
                bseq.setdefault(base.walks[i]['seq'], i)
            tmap = dict((i, bseq.get(w['seq'], -1)) for i, w in enumerate(tk_graph.walks))
            for p in placed:
                p['compat'] = [tmap[x] for x in p['compat'] if tmap.get(x, -1) >= 0]
                if p['place']:
                    wi = tmap.get(p['place'][0], -1)
                    p['place'] = (wi, p['place'][1], p['place'][2]) if wi >= 0 else None
        log('placed %d reads in %.1f s' % (len(reads), time.time() - t0), quiet)
    status_n = collections.Counter()
    compat_n, unique_n, compat_mq5, placed_n = (collections.Counter() for _ in range(4))
    reads_out = []
    for rd, pl in zip(reads, placed):
        status_n[pl['status']] += 1
        for wi in pl['compat']:
            compat_n[wi] += 1
            if rd['mapq'] >= 5:
                compat_mq5[wi] += 1
        if len(pl['compat']) == 1:
            unique_n[pl['compat'][0]] += 1
        place = pl.get('place')
        if place and not pl['partial']:
            placed_n[place[0]] += 1
        strand_ref = '+' if (rd['strand'] == '+') == (pl.get('ori', 1) == 1) else '-'
        first, last = rd['steps'][0], rd['steps'][-1]
        reads_out.append([
            rd['name'], rd['mapq'], strand_ref, rd['qlen'],
            round(rd['matches'] / rd['blen'], 3) if rd['blen'] else 0,
            len(rd['steps']), '%s%s' % (first[1], first[0]), '%s%s' % (last[1], last[0]),
            pl.get('inside_all', 0), pl['status'], pl['compat'],
            place[0] if place else -1, place[1] if place else -1, place[2] if place else -1,
            pl['partial'], rd['pd'] or ''])
    path_counts = [{'walk': wi, 'compat': compat_n[wi], 'unique': unique_n[wi], 'compat_mq5': compat_mq5[wi],
                    'placed': placed_n[wi]} for wi in range(len(base.walks))]
    for r in rows:
        if 'walk' in r:
            r['reads'] = compat_n[r['walk']]
            r['reads_unique'] = unique_n[r['walk']]

    # ---- annotation: census, region.json
    census = census_lookup(contig, req1 - 1, req2, A1 - 1, B1)
    locus = census['locus']
    meta = region_meta(src, census, locus)
    trs = []
    main_region = None
    if locus and locus.get('main_region_id'):
        main_region = next((r for r in census['regions'] if r['region_id'] == locus['main_region_id']), None)
    vntr_regs = [r for r in census['regions'] if r.get('tr_class') == 'VNTR']
    if main_region is None and vntr_regs:
        ov = [r for r in vntr_regs if int(r['start']) < req2 and int(r['end']) > req1 - 1] or vntr_regs
        main_region = max(ov, key=lambda r: int(r['length']))
    for r in census['regions']:
        trs.append({'id': r['region_id'], 's': int(r['start']) + 1, 'e': int(r['end']), 'cls': r.get('tr_class'),
                    'period': num(r.get('period')), 'motif': r.get('motif', ''), 'copies': num(r.get('copies')),
                    'gc': num(r.get('region_gc')), 'src': r.get('src', ''),
                    'vg_fp': num(r.get('vg_fp')), 'vg_fn': num(r.get('vg_fn')),
                    'main': int(main_region is not None and r['region_id'] == main_region['region_id'])})
    if not trs and meta.get('period'):
        trs.append({'id': meta.get('region_id') or 'core', 's': req1, 'e': req2, 'cls': meta.get('tr_class') or 'VNTR',
                    'period': num(meta.get('period')), 'motif': meta.get('motif') or '', 'copies': num(meta.get('copies')),
                    'gc': None, 'src': 'region.json', 'vg_fp': num(meta.get('vg_fp')), 'vg_fn': num(meta.get('vg_fn')),
                    'main': 1})

    # ---- calls & truth records
    callrows, truthrows = load_calls(tk, src.pkg)
    n = collections.Counter()
    for c in callrows:
        if c['called'] and c['truvari'] in ('FP', 'TP'):
            n[('vg' if c['caller'] == 'vg' else 'pg') + '_' + c['truvari'].lower()] += 1
    for t in truthrows:
        if t['source'] == 'stvar':
            if t['vg'] == 'FN':
                n['vg_fn'] += 1
            if t['pg'] == 'FN':
                n['pg_fn'] += 1
            if t['sv50']:
                n['truth_sv'] += 1
    counts = {k: n[k] for k in ('vg_fp', 'vg_tp', 'vg_fn', 'pg_fp', 'pg_tp', 'pg_fn', 'truth_sv')}
    counts['source'] = 'calls in the span' if (callrows or truthrows) else 'none'

    # ---- MSA strings: runs of >= 4 gaps written as -<n>; (halves the page; decoded by the page)
    for m in msas:
        m['strs'] = [GAP_RUN.sub(lambda x: '-%d;' % len(x.group()), t) for t in m['strs']]
        m['enc'] = 'gaprun'

    # ---- Bandage images into data URIs, within the page budget
    budget = args.max_page_mb * 1e6 * 0.8 - 400000 - sum(len(s) for m in msas for s in m['strs']) - 60 * len(reads_out)
    total_png = sum(len(g['png'] or b'') for g in graphs_out)
    if total_png * 4 / 3 > budget and not args.no_bandage:
        for gi, g in enumerate(graphs_out):
            if g['png'] and len(g['png']) * 4 / 3 > budget / len(graphs_out):
                data, _ = bandage_png(graphs[gi], graph_metrics_node_use(graphs[gi]), vdir, quiet, height=800)
                if data and len(data) < len(g['png']):
                    g['png'] = data
    for g in graphs_out:
        g['bandage'] = ('data:image/png;base64,' + base64.b64encode(g['png']).decode()) if g['png'] else None
        del g['png']

    rs = (S.get('reads') or read_json(os.path.join(tk, 'reads_summary.json'))) if tk else {}
    D = {
        'contig': contig, 'req': [req1, req2], 'span': [A1, B1], 'name': name, 'region_id': src.region_id,
        'meta': meta, 'package_dir': src.pkg, 'toolkit_dir': tk, 'view_dir': vdir,
        'boundaries': (S.get('boundaries') if S else None) or {
            'L': {'node': src.rj.get('anchor_left')}, 'R': {'node': src.rj.get('anchor_right')}},
        'tandem_scan': S.get('reference_tandem_scan') if S else vr.tandem_scan(refseq),
        'called_vs_truth': chv or None, 'closest': closest, 'reads_summary': rs,
        'rows': rows, 'msas': msas, 'graphs': graphs_out, 'same_seq_groups': same_seq_groups,
        'trs': trs, 'calls': callrows, 'truth': truthrows,
        'reads': reads_out if tk else None, 'n_reads_total': n_reads_total, 'read_status': dict(status_n),
        'path_counts': path_counts, 'counts': counts,
        'census': {'locus': strip_private(locus), 'strata': census['strata'],
                   'main_region': main_region and main_region['region_id']},
        'share_bins': SHARE_BINS[:-1] + [[17, None]], 'share_cols': SHARE_COLS, 'ref_grey': REF_GREY,
        'warnings': (S.get('warnings') or [] if S else []) + warnings,
        'generated': time.strftime('%Y-%m-%d %H:%M'), 'build_seconds': None,
    }
    D['build_seconds'] = round(time.time() - t_start, 1)
    page = render_html(D)
    out = args.out or os.path.join(args.workroot, name + '.html')
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, 'w') as f:
        f.write(page)
    size = os.path.getsize(out)
    if size > args.max_page_mb * 1e6:
        log('WARNING: %s is %.1f MB (> %g MB); lower --max-cols or --max-reads, or use --no-bandage'
            % (out, size / 1e6, args.max_page_mb), False)
    summary_line = {
        'html': out, 'bytes': size, 'span': '%s:%d-%d' % (contig, A1, B1), 'rows': len(rows),
        'msas': [m['name'] for m in msas], 'graphs': [
            {'name': g['name'], 'nodes': g['metrics']['nodes'], 'cost_over_opt': g['metrics'].get('cost_over_opt'),
             'kmer_extra': g['metrics']['kmer']['frac_extra']} for g in graphs_out],
        'reads': len(reads_out), 'read_status': dict(status_n), 'seconds': D['build_seconds'],
        'warnings': len(D['warnings'])}
    print(json.dumps(summary_line))
    return summary_line


EVAL_KEYS = [('size', 'nodes'), ('size', 'nodes_per_kb'), ('size', 'node_len_mean'), ('size', 'node_frac_1bp'),
             ('size', 'bubbles'), ('size', 'bubbles_sv'), ('size', 'max_alleles_any_bubble'),
             ('alignment', 'all_cost_over_opt'), ('alignment', 'ref_cost_over_opt'),
             ('alignment', 'affine_all_cost_over_opt'), ('alignment', 'all_excess_per_kb'),
             ('alignment', 'all_U_per_kb'), ('alignment', 'pairs_sampled'),
             ('inflation', 'sv_pieces_median'), ('inflation', 'indel_bp_ratio'),
             ('redundancy', 'kmer_frac_extra'),
             ('truth', 'h1_d_graph'), ('truth', 'h2_d_graph'), ('truth', 'd_graph_sum'), ('truth', 'd_panel_sum'),
             ('truth.truvari', 'raw_f1'), ('truth.truvari', 'refined_f1')]


def method_of(gfa_path):
    """Method name of a candidate GFA: its directory (candidates/<method>/<region_id>.gfa)."""
    return os.path.basename(os.path.dirname(os.path.abspath(gfa_path)))


def evaluate_results(region_id, method, gfa_path):
    """Key numbers from tools/evaluate.py's results/<method>/<region_id>.json, if it exists."""
    if not region_id:
        return None
    p = os.path.join(config.RESULTS_DIR, method, region_id + '.json')
    if not os.path.exists(p):
        return None
    try:
        r = read_json(p)
    except ValueError:
        return None
    out = {'file': p, 'status': r.get('status'), 'valid': (r.get('validity') or {}).get('valid'),
           'graph': r.get('graph'),
           'graph_matches': bool(r.get('graph')) and os.path.basename(r['graph']) == os.path.basename(gfa_path),
           'm': {}}
    for sec, k in EVAL_KEYS:
        d = r
        for part in sec.split('.'):
            d = d.get(part) if isinstance(d, dict) else None
        if isinstance(d, dict) and k in d:
            out['m'][sec.split('.')[0] + '.' + k] = d[k]
    return out


def graph_metrics_node_use(G):
    use = collections.Counter()
    for w in G.walks:
        for n in set(x for x, _ in w['steps']):
            use[n] += w['weight']
    return use


def encode_runs(steps, nidx):
    """Walk as runs of consecutive node indices: [signed index+1, run length, index step]."""
    out = []
    k = 0
    while k < len(steps):
        n, o = steps[k]
        x = nidx[n]
        sgn = 1 if o == '>' else -1
        d, c = 0, 1
        if k + 1 < len(steps) and steps[k + 1][1] == o and abs(nidx[steps[k + 1][0]] - x) == 1:
            d = nidx[steps[k + 1][0]] - x
            while k + c < len(steps) and steps[k + c][1] == o and nidx[steps[k + c][0]] == x + d * c:
                c += 1
        out.extend([sgn * (x + 1), c, d])
        k += c
    return out


def crop_msa(msa_rows, rows, seqs, refrow_aln, A1, req1, req2, max_cols):
    """msa_rows: per row aligned string (or None), each covering its sequence from base0.
    Crops to max_cols columns centred on the core interval. Returns the MSA dict parts."""
    ncols = len(refrow_aln)
    base0 = [r.get('_b0', 0) for r in rows]
    crop = None
    c_lo, c_hi = 0, ncols
    if ncols > max_cols:
        b2c = [c for c, ch in enumerate(refrow_aln) if ch != '-']
        o1 = max(0, min(len(b2c) - 1, req1 - A1 - base0[0]))
        o2 = max(0, min(len(b2c) - 1, req2 - A1 - base0[0]))
        cmid = (b2c[o1] + b2c[o2]) // 2
        c_lo = max(0, min(cmid - max_cols // 2, ncols - max_cols))
        c_hi = c_lo + max_cols
        crop = [c_lo, c_hi]
    strs, sidx, a, b0 = [], {}, [], []
    for i, s in enumerate(msa_rows):
        if s is None:
            a.append(-1)
            b0.append(0)
            continue
        b0.append(base0[i] + sum(1 for ch in s[:c_lo] if ch != '-'))
        t = s[c_lo:c_hi]
        if t not in sidx:
            sidx[t] = len(strs)
            strs.append(t)
        a.append(sidx[t])
    return {'strs': strs, 'a': a, 'b0': b0, 'ncols': c_hi - c_lo, 'ncols_full': ncols, 'crop': crop}


def build_mafft_msa(rows, seqs, refseq, base, A1, req1, req2, args, vdir, warnings, S, tk):
    """mafft over every row's sequence; long spans are first cut to the CHM13 window around the
    core (projected through each row's alignment to CHM13)."""
    rseq = [seqs[r['sid']] for r in rows]
    for r in rows:
        r['_b0'] = 0
    precrop = None
    if max(len(s) for s in rseq) > args.max_seq:
        c = (req1 + req2) // 2 - A1
        half = args.max_seq // 2 - 1
        w0 = max(0, min(req1 - A1 - 1000, c - half))
        w1 = min(len(refseq), max(req2 - A1 + 1000, c + half))
        if w1 - w0 > args.max_seq:
            w0, w1 = max(0, c - half), min(len(refseq), c + half)
        precrop = [w0, w1]
        warnings.append('span is longer than --max-seq %d; every sequence was cut to the part aligned to CHM13 '
                        '%s:%d-%d before the mafft MSA' % (args.max_seq, base.walks[0]['names'][0].split('#')[-1]
                                                          if base.walks[0]['names'] else '', A1 + w0, A1 + w1 - 1))
        cig = {}
        if tk:
            for h in read_tsv(os.path.join(tk, 'haplotypes.tsv')):
                cig[h.get('cigar')] = None
        for i, r in enumerate(rows):
            s = rseq[i]
            if r['kind'] == 'ref':
                rseq[i], r['_b0'] = refseq[w0:w1], w0
                continue
            if len(refseq) * len(s) > 4e8:
                warnings.append('row %s dropped from the mafft MSA: too long to project' % r['label'])
                rseq[i] = ''
                continue
            cg = vr.align_cigar(refseq, s)[1]
            m = cigar_ref_to_query(cg)
            q0, q1 = m[w0], m[w1]
            rseq[i], r['_b0'] = s[q0:q1], q0
    uniq, uidx = [], {}
    for s in rseq:
        if s and s not in uidx:
            uidx[s] = len(uniq)
            uniq.append(s)
    if not os.path.exists(config.MAFFT):
        raise ViewError('mafft not found (%s)' % config.MAFFT)
    key = None
    if len(uniq) == 1:
        aln = uniq[:]
    else:
        aln, key = run_mafft(uniq, vdir, args.mafft, args.threads, args.quiet)
    msa_rows = [aln[uidx[s]] if s else None for s in rseq]
    m = crop_msa(msa_rows, rows, seqs, msa_rows[0], A1, req1, req2, args.max_cols)
    for r in rows:
        r.pop('_b0', None)
    m.update({'name': 'mafft %s (all rows)' % args.mafft, 'kind': 'mafft', 'precrop': precrop, 'key': key,
              'note': 'mafft %s of every row on the page (CHM13, panel walks%s)'
                      % (' '.join(mafft_args(args.mafft)),
                         ''.join(', ' + k for k in ('truth', 'called') if any(r['kind'] == k for r in rows)))})
    if m['crop']:
        warnings.append('%s has %d columns; showing columns %d-%d (--max-cols %d)' % (
            m['name'], m['ncols_full'], m['crop'][0] + 1, m['crop'][1], args.max_cols))
    return m


def build_file_msa(path, label, rows, seqs, A1, req1, req2, args, vdir, warnings):
    """A candidate's own MSA (FASTA, rows named like hap32.fa). Rows are matched to the page rows
    by sequence; rows it lacks (truth, calls) are added with mafft --add, which keeps the
    candidate's alignment fixed."""
    recs = read_fasta(path)
    if not recs:
        raise ViewError('empty MSA')
    L = set(len(s) for _, s in recs)
    if len(L) != 1:
        raise ViewError('rows of different lengths')
    by_seq = {}
    unmatched = []
    known = set(seqs[r['sid']] for r in rows)
    for nm, s in recs:
        s = s.upper().replace('.', '-')
        u = s.replace('-', '')
        if u not in known:
            unmatched.append(first_word(nm))
        by_seq.setdefault(u, s)
    note = []
    if unmatched:
        note.append('%d MSA rows spell no page sequence (%s)' % (len(unmatched), ', '.join(unmatched[:4])))
    need = [seqs[r['sid']] for r in rows if seqs[r['sid']] not in by_seq]
    need = list(collections.OrderedDict.fromkeys(need))
    added = 0
    if need and not args.no_add:
        existing = list(collections.OrderedDict.fromkeys(by_seq.values()))
        e2, n2 = mafft_add(existing, need, vdir, args.threads, args.quiet)
        if e2 is None:
            note.append('truth and called rows could not be added (%s)' % n2)
        else:
            remap = dict(zip(existing, e2))
            by_seq = dict((u, remap[s]) for u, s in by_seq.items())
            for u, s in zip(need, n2):
                by_seq[u] = s
            added = len(need)
            note.append('%d rows it lacks (truth, calls) were added with mafft --add, which keeps its '
                        'alignment fixed' % added)
    msa_rows = [by_seq.get(seqs[r['sid']]) for r in rows]
    for r in rows:
        r['_b0'] = 0
    if msa_rows[0] is None:
        raise ViewError('no row spells the CHM13 sequence')
    ncols = len(msa_rows[0])
    msa_rows = [s if s is None or len(s) == ncols else None for s in msa_rows]
    m = crop_msa(msa_rows, rows, seqs, msa_rows[0], A1, req1, req2, args.max_cols)
    for r in rows:
        r.pop('_b0', None)
    m.update({'name': '%s MSA' % label, 'kind': 'file', 'precrop': None, 'key': None, 'src': os.path.abspath(path),
              'note': '; '.join(['from %s' % os.path.basename(path)] + note)})
    for x in note:
        if 'could not' in x or 'spell no' in x:
            warnings.append('%s: %s' % (m['name'], x))
    return m


def region_meta(src, census, locus):
    rj = src.rj or {}
    meta = {'region_id': src.region_id, 'source': 'region.json' if rj else ('census' if locus else None)}
    for k in ('stratum', 'period', 'motif', 'copies', 'tr_class', 'notes', 'vg_fp', 'vg_fn', 'pg_fp', 'pg_fn',
              'in_benchmark', 'n_hap32', 'n_hprc', 'anchor_left', 'anchor_right', 'matched_to'):
        if rj.get(k) not in (None, ''):
            meta[k] = rj[k]
    if locus:
        for k, k2 in (('tr_class', 'tr_class'), ('period', 'period'), ('motif', 'motif'), ('copies', 'copies'),
                      ('gc', 'region_gc'), ('region_len', 'region_length'), ('main_region_id', 'main_region_id'),
                      ('error_status', 'error_status'), ('vg_fp', 'vg_fp'), ('vg_fn', 'vg_fn'), ('pg_fp', 'pg_fp'),
                      ('pg_fn', 'pg_fn'), ('fp_pctseq', 'vg_fp_median_pctseq'), ('locus_id', 'locus_id')):
            if locus.get(k2) not in (None, '') and k not in meta:
                meta[k] = locus[k2]
    strata = sorted(set(s['stratum'] + (' (matched to %s)' % s['matched_to'] if s.get('matched_to') else '')
                        for s in census['strata']))
    if strata and 'stratum' not in meta:
        meta['stratum'] = ', '.join(strata)
    return meta


def load_calls(tk, pkg):
    """vg / PanGenie call rows and truth records: from region.py output when present, else from the
    package's calls.tsv (the same columns; truth records as caller=truth or source=stvar/smvar)."""
    calls, truth = [], []
    if tk and os.path.exists(os.path.join(tk, 'calls.tsv')):
        crow, trow = read_tsv(os.path.join(tk, 'calls.tsv')), read_tsv(os.path.join(tk, 'truth_records.tsv'))
    elif pkg:
        # package calls.tsv (tools/package_regions.py): one table, 'source' = vg | pangenie |
        # truth_stvar | truth_smvar; truth rows carry max_allele_len_diff in len_diff and the
        # alt allele lengths in alt_len
        allr = read_tsv(os.path.join(pkg, 'calls.tsv'))
        crow, trow = [], []
        for r in allr:
            s = (r.get('source') or r.get('caller') or '').lower()
            if s in ('vg', 'pangenie'):
                r = dict(r, caller=s)
                crow.append(r)
            elif s.startswith('truth') or s in ('stvar', 'smvar'):
                r = dict(r, source=s.replace('truth_', ''))
                r.setdefault('alt_lens', r.get('alt_len', ''))
                trow.append(r)
    else:
        return calls, truth
    for c in crow:
        caller = 'pangenie' if c.get('caller', '').lower() == 'pangenie' else c.get('caller')
        ld = num(c.get('len_diff'))
        if ld is None and num(c.get('alt_len')) is not None and num(c.get('ref_len')) is not None:
            ld = num(c.get('alt_len')) - num(c.get('ref_len'))
        calls.append({'caller': caller, 'pos': num(c.get('pos')), 'end': num(c.get('end')), 'id': c.get('id', ''),
                      'ref_len': num(c.get('ref_len')), 'allele': num(c.get('allele')), 'alt_len': num(c.get('alt_len')),
                      'len_diff': ld, 'called': num(c.get('called'), 1), 'gt': c.get('gt', ''),
                      'qual': num(c.get('qual'), c.get('qual')), 'filter': c.get('filter', ''),
                      'in_interval': num(c.get('in_interval')), 'truvari': c.get('truvari', '-') or '-',
                      'aardvark': c.get('aardvark', '-') or '-', 'pctseq': num(c.get('PctSeqSimilarity')),
                      'pctsize': num(c.get('PctSizeSimilarity')), 'sizediff': num(c.get('SizeDiff')),
                      'startdist': num(c.get('StartDistance')), 'edge': num(c.get('overlaps_span_edge'))})
    for t in trow:
        dl = num(t.get('max_allele_len_diff'), num(t.get('len_diff')))
        truth.append({'source': t.get('source') or 'stvar', 'pos': num(t.get('pos')), 'end': num(t.get('end')),
                      'ref_len': num(t.get('ref_len')), 'alt_lens': t.get('alt_lens', t.get('alt_len', '')),
                      'dlen': dl, 'sv50': num(t.get('is_sv50'), 1 if abs(dl or 0) >= 50 else 0),
                      'svtype': t.get('svtype', ''), 'gt': t.get('gt', ''), 'h1': t.get('applied_h1', ''),
                      'h2': t.get('applied_h2', ''), 'in_req': num(t.get('in_requested_interval')),
                      'bench': num(t.get('in_sv_benchmark')), 'vg': t.get('vg_truvari', '-') or '-',
                      'pg': t.get('pg_truvari', '-') or '-', 'trf_period': num(t.get('TRFperiod')),
                      'trf_motif': t.get('TRFrepeat', '.'), 'trf_copies': num(t.get('TRFcopies')),
                      'rm': t.get('RM_clsfam', '.'), 'lcr': num(t.get('LCR'))})
    return calls, truth


def strip_private(d):
    if not d:
        return d
    return dict((k, v) for k, v in d.items() if not k.startswith('_'))


# ---------------------------------------------------------------- HTML

def render_html(D):
    rid = D.get('region_id')
    title = '%s:%s-%s' % (D['contig'], format(D['req'][0], ','), format(D['req'][1], ','))
    data = json.dumps(D, separators=(',', ':')).replace('</', '<\\/')
    page = TEMPLATE.replace('__TITLE__', html.escape(('%s %s' % (rid, title)) if rid else title))
    page = page.replace('__CSS__', CSS).replace('__JS__', JS).replace('__DATA__', data)
    return page


TEMPLATE = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>__CSS__</style>
</head>
<body>
<header id="hdr"></header>
<nav id="toc"><a href="#p-aln">Alignment</a><a href="#p-graph">Graphs</a><a href="#p-calls">Calls</a><a href="#p-reads">Reads</a>
<button id="themebtn" type="button" title="cycle theme: auto / light / dark">theme: auto</button></nav>
<main>
<section id="p-aln" class="panel">
  <h2>Alignment <span class="sub" id="aln-sub"></span></h2>
  <div class="controls" id="aln-controls"></div>
  <div class="legend" id="aln-legend"></div>
  <canvas id="overview" class="overview"></canvas>
  <div class="track" id="aln-track"><div class="labels" id="aln-labels"></div>
    <div class="scroller msa-scroll" id="aln-scroll"><div class="spacer"></div><canvas id="aln-canvas"></canvas></div></div>
  <details class="more"><summary>Distinct walks of the baseline graph</summary><div id="walk-table"></div></details>
</section>
<section id="p-graph" class="panel">
  <h2>Graphs <span class="sub" id="graph-sub"></span></h2>
  <div id="metrics-table"></div>
  <div id="graph-warn"></div>
  <div class="bgrid" id="bandage-grid"></div>
  <div class="legend" id="bandage-legend"></div>
  <h3>Node tracks</h3>
  <div class="controls" id="node-controls"></div>
  <div class="legend" id="node-legend"></div>
  <div id="node-blocks"></div>
</section>
<section id="p-calls" class="panel">
  <h2>Calls <span class="sub" id="calls-sub"></span></h2>
  <div id="calls-body">
  <div id="called-dist"></div>
  <h3>vg call alleles</h3><div class="controls" id="vg-controls"></div><div id="vg-table"></div>
  <h3>PanGenie called alleles</h3><div class="controls" id="pg-controls"></div><div id="pg-table"></div>
  <h3>Truth records (T2T-Q100 v1.1)</h3><div class="controls" id="truth-controls"></div><div id="truth-table"></div>
  </div>
</section>
<section id="p-reads" class="panel">
  <h2>Reads <span class="sub" id="reads-sub"></span></h2>
  <div id="reads-body">
  <div class="readsgrid"><div id="mapq-hist"></div><div id="reads-facts"></div></div>
  <h3>Reads per walk</h3>
  <div id="path-table"></div>
  <h3>Pileup on the walk each read follows</h3>
  <div class="controls" id="pile-controls"></div>
  <div class="legend" id="pile-legend"></div>
  <div class="track" id="pile-track"><div class="labels" id="pile-labels"></div>
    <div class="scroller msa-scroll" id="pile-scroll"><div class="spacer"></div><canvas id="pile-canvas"></canvas></div></div>
  <h3>Read alignments</h3>
  <div class="controls" id="rt-controls"></div>
  <div id="reads-table"></div>
  </div>
</section>
<footer id="ftr"></footer>
</main>
<div id="tip" class="tip" hidden></div>
<script type="application/json" id="data">__DATA__</script>
<script>__JS__</script>
</body>
</html>
'''

CSS = r'''
:root{--bg:#ffffff;--panel:#f6f7f9;--fg:#1d2127;--muted:#5f6873;--border:#d9dde3;--grid:#e6e9ee;
--match:#dde1e6;--matchtext:#8a929c;--mismatch:#e8590c;--ins:#7048e8;--delline:#9aa2ab;--sel:#1c7ed6;
--tp:#2f9e44;--fp:#e03131;--fn:#e67700;--na:#868e96;--mq0:#e03131;--mq1:#f08c00;--mq5:#c99a06;--mq30:#4dabf7;--mq60:#1864ab;
--ref1:#c3c8cf;--ref2:#a4abb4;--vntr:#ae3ec9;--str:#15aabf;--oth:#868e96;--tick:#495057;--hl:#fab005;
--truth:#2b8a3e;--called:#9c36b5;--panelc:#495057;--refc:#1d2127;--link:#1971c2;--warnbg:#fff4e6;--warnfg:#8f4a00;--absent:#f1f3f5;}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#121417;--panel:#1b1e23;--fg:#e6e8eb;--muted:#9aa3ad;--border:#30353c;--grid:#262a30;
--match:#3b4149;--matchtext:#7d8590;--mismatch:#ff8a3d;--ins:#9775fa;--delline:#6b737c;--sel:#4dabf7;
--tp:#51cf66;--fp:#ff6b6b;--fn:#ffa94d;--na:#868e96;--mq0:#ff6b6b;--mq1:#ffa94d;--mq5:#ffd43b;--mq30:#74c0fc;--mq60:#339af0;
--ref1:#4a5058;--ref2:#5f666f;--vntr:#da77f2;--str:#3bc9db;--oth:#868e96;--tick:#ced4da;--hl:#ffd43b;
--truth:#69db7c;--called:#e599f7;--panelc:#adb5bd;--refc:#f1f3f5;--link:#74c0fc;--warnbg:#3b2a12;--warnfg:#ffd8a8;--absent:#1e2126;}}
:root[data-theme="dark"]{--bg:#121417;--panel:#1b1e23;--fg:#e6e8eb;--muted:#9aa3ad;--border:#30353c;--grid:#262a30;
--match:#3b4149;--matchtext:#7d8590;--mismatch:#ff8a3d;--ins:#9775fa;--delline:#6b737c;--sel:#4dabf7;
--tp:#51cf66;--fp:#ff6b6b;--fn:#ffa94d;--na:#868e96;--mq0:#ff6b6b;--mq1:#ffa94d;--mq5:#ffd43b;--mq30:#74c0fc;--mq60:#339af0;
--ref1:#4a5058;--ref2:#5f666f;--vntr:#da77f2;--str:#3bc9db;--oth:#868e96;--tick:#ced4da;--hl:#ffd43b;
--truth:#69db7c;--called:#e599f7;--panelc:#adb5bd;--refc:#f1f3f5;--link:#74c0fc;--warnbg:#3b2a12;--warnfg:#ffd8a8;--absent:#1e2126;}
*{box-sizing:border-box}
html,body{margin:0;overflow-wrap:anywhere;background:var(--bg);color:var(--fg);font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
header,main,nav{padding:0 16px;max-width:100%}
header{padding-top:14px}
h1{font-size:20px;margin:0 0 4px}
h2{font-size:17px;margin:0 0 8px}
h3{font-size:14px;margin:14px 0 6px}
.sub{font-weight:400;color:var(--muted);font-size:13px}
.panel{border:1px solid var(--border);border-radius:8px;background:var(--bg);padding:12px 14px;margin:12px 0}
nav#toc{position:sticky;top:0;z-index:5;background:var(--bg);border-bottom:1px solid var(--border);display:flex;gap:16px;align-items:center;padding-top:6px;padding-bottom:6px}
nav#toc a{color:var(--link);text-decoration:none}
nav#toc button{margin-left:auto}
button,select,input{font:inherit;font-size:13px;color:var(--fg);background:var(--panel);border:1px solid var(--border);border-radius:5px;padding:2px 7px}
button{cursor:pointer}
.mono,code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:12px}
.line{margin:2px 0}
.kv{display:inline-block;margin-right:14px;white-space:nowrap}
.kv b{font-weight:600}
.muted{color:var(--muted)}
.warn{background:var(--warnbg);color:var(--warnfg);border-radius:6px;padding:6px 10px;margin:6px 0;font-size:13px}
.controls{display:flex;flex-wrap:wrap;gap:6px 14px;align-items:center;margin:6px 0;font-size:13px}
.controls label{display:inline-flex;gap:4px;align-items:center;white-space:nowrap}
.legend{font-size:12px;color:var(--muted);margin:4px 0 6px;line-height:1.8}
.sw{display:inline-block;width:12px;height:10px;border-radius:2px;vertical-align:-1px;margin-right:4px;border:1px solid var(--border)}
.track{display:flex;align-items:flex-start;border:1px solid var(--border);border-radius:6px;overflow:hidden;background:var(--bg)}
.labels{flex:0 0 auto;width:200px;border-right:1px solid var(--border);background:var(--panel);font-size:12px;overflow:hidden}
.labels .lab{display:flex;align-items:center;gap:4px;padding:0 4px;white-space:nowrap;overflow:hidden;cursor:pointer;border-bottom:1px solid transparent}
.labels .lab:hover{background:var(--grid)}
.labels .lab.sel{outline:2px solid var(--sel);outline-offset:-2px}
.labels .lab .x{margin-left:auto;color:var(--muted);padding:0 3px}
.labels .lab .x:hover{color:var(--fp)}
.labels .hd{color:var(--muted);padding:0 4px;white-space:nowrap;overflow:hidden;font-size:11px}
.badge{font-size:10px;border-radius:3px;padding:0 3px;border:1px solid var(--border);color:var(--muted)}
.k-ref{color:var(--refc);font-weight:600}.k-truth{color:var(--truth);font-weight:600}.k-called{color:var(--called)}.k-panel{color:var(--panelc)}
.scroller{flex:1 1 auto;overflow-x:auto;overflow-y:hidden;position:relative;min-width:0}
.scroller .spacer{height:1px}
.scroller canvas{position:sticky;left:0;display:block;margin-top:-1px}
canvas.overview{width:100%;height:60px;display:block;border:1px solid var(--border);border-radius:6px;margin:4px 0 6px;cursor:pointer}
.tip{position:fixed;z-index:20;max-width:460px;background:var(--panel);color:var(--fg);border:1px solid var(--border);border-radius:6px;padding:6px 8px;font-size:12px;pointer-events:none;box-shadow:0 2px 10px rgba(0,0,0,.18);white-space:normal;word-break:break-word}
table{border-collapse:collapse;font-size:12px;width:100%}
th,td{border-bottom:1px solid var(--grid);padding:2px 6px;text-align:left;white-space:nowrap}
th{position:sticky;top:0;background:var(--panel);cursor:pointer;user-select:none;font-weight:600}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
.tablewrap{max-height:460px;overflow:auto;border:1px solid var(--border);border-radius:6px}
.s-TP{color:var(--tp);font-weight:600}.s-FP{color:var(--fp);font-weight:600}.s-FN{color:var(--fn);font-weight:600}
table.mt{width:auto;min-width:50%}
table.mt th{cursor:default}
table.mt td.grp{color:var(--muted);font-size:11px;padding-top:8px;border-bottom:none}
table.mt td.best{font-weight:700}
table.mt th{white-space:normal;max-width:180px;overflow-wrap:normal;word-break:normal}
table.mt td:first-child{white-space:normal;min-width:150px}
.readsgrid{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:14px;align-items:start}
.bgrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:12px;align-items:start;margin-top:10px}
.bgrid figure.full{grid-column:1/-1;overflow:auto;max-height:85vh}
.bgrid figure img{width:100%;height:auto;background:#fff;border:1px solid var(--border);border-radius:6px;cursor:zoom-in;display:block}
.bgrid figure.full img{width:auto;max-width:none;cursor:zoom-out}
.btitle{font-size:13px;margin-bottom:4px}
@media (max-width:820px){.readsgrid{grid-template-columns:1fr}.labels{width:130px}.kv,.controls label{white-space:normal}code{word-break:break-all}nav#toc{gap:10px;flex-wrap:wrap}}
figure{margin:0}
figcaption{font-size:12px;color:var(--muted);margin-top:4px}
.status{font-size:12px;color:var(--muted);min-height:18px;margin-top:4px}
details.more{margin-top:8px}
details.more summary{cursor:pointer;color:var(--link)}
footer{font-size:12px;color:var(--muted);padding:6px 0 30px}
.facts div{margin:2px 0}
.ntblock{margin-top:10px}
.ntblock h4{font-size:13px;margin:8px 0 4px;font-weight:600}
svg text{fill:var(--fg)}
'''

JS = r'''
'use strict';
const D = JSON.parse(document.getElementById('data').textContent);
const $ = (id) => document.getElementById(id);
const esc = (s) => String(s === null || s === undefined ? '' : s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
const fmt = (x) => (x === null || x === undefined || x === '') ? '–' : (typeof x === 'number' ? (Number.isInteger(x) ? x.toLocaleString('en-US') : (Math.abs(x) < 1 ? x.toFixed(3) : x.toFixed(Math.abs(x) < 10 ? 2 : 1))) : (typeof x === 'boolean' ? (x ? 'yes' : 'no') : x));
const pct = (x) => (x === null || x === undefined || x === '') ? '–' : (100 * x).toFixed(1) + '%';
let COL = {};
function css(n) { return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }
function isDark() { const t = document.documentElement.getAttribute('data-theme'); if (t) return t === 'dark'; return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches; }
function loadColors() {
  for (const k of ['bg','panel','fg','muted','border','grid','match','matchtext','mismatch','ins','delline','sel','tp','fp','fn','na',
                   'mq0','mq1','mq5','mq30','mq60','ref1','ref2','vntr','str','oth','tick','hl','truth','called','panelc','refc','absent'])
    COL[k] = css('--' + k);
  COL.dark = isDark();
}
const tip = $('tip');
function showTip(ev, htmlText) { tip.innerHTML = htmlText; tip.hidden = false;
  const w = tip.offsetWidth, h = tip.offsetHeight; let x = ev.clientX + 14, y = ev.clientY + 14;
  if (x + w > window.innerWidth - 8) x = ev.clientX - w - 14; if (y + h > window.innerHeight - 8) y = ev.clientY - h - 14;
  tip.style.left = Math.max(4, x) + 'px'; tip.style.top = Math.max(4, y) + 'px'; }
function hideTip() { tip.hidden = true; }
function sw(c) { return '<span class="sw" style="background:' + c + '"></span>'; }

// ------------------------------------------------------------------ data prep
for (const m of D.msas) if (m.enc === 'gaprun') m.strs = m.strs.map(s => s.replace(/-(\d+);/g, (x, n) => '-'.repeat(+n)));
const ROWS = D.rows, NR = ROWS.length, A1 = D.span[0];
function expandRuns(r) { const out = []; for (let i = 0; i < r.length; i += 3) { const s = r[i], n = r[i + 1], d = r[i + 2];
  const sg = s > 0 ? 1 : -1, x = Math.abs(s) - 1; for (let k = 0; k < n; k++) out.push(sg * (x + d * k + 1)); } return out; }
const GR = D.graphs.map((g, gi) => { const N = g.nodes;
  const walks = g.walks.map(w => { const steps = expandRuns(w.runs); const off = new Int32Array(steps.length + 1); let o = 0;
    for (let i = 0; i < steps.length; i++) { off[i] = o; o += N.len[Math.abs(steps[i]) - 1]; } off[steps.length] = o;
    return Object.assign({}, w, {steps, off}); });
  return {gi, g, N, walks, sel: -1, visits: new Map()}; });
const rowOfBaseWalk = new Map();
ROWS.forEach((r, i) => { if (r.walk !== undefined && !rowOfBaseWalk.has(r.walk)) rowOfBaseWalk.set(r.walk, i); });
const MAIN = D.trs.find(r => r.main);
const PERIOD = MAIN && MAIN.period ? Math.round(MAIN.period) : null;

// ------------------------------------------------------------------ the MSA on show
let MI = 0, M = null, NC = 0, SEQ = null;
let col2ref = null, refGap = null, refCol = [], refFirst = A1;
const colNodeCache = new Map(), splitCache = new Map(), b2cCache = new Map(), tickCache = new Map();
function rowStr(r) { const a = M.a[r]; return a >= 0 ? SEQ[a] : null; }
function setMsa(i) {
  MI = i; M = D.msas[i]; NC = M.ncols; SEQ = M.strs;
  colNodeCache.clear(); splitCache.clear(); b2cCache.clear(); tickCache.clear();
  col2ref = new Int32Array(NC); refGap = new Uint8Array(NC); refCol = [];
  const rs = rowStr(0) || ''; let k = M.b0[0];
  for (let c = 0; c < NC; c++) { if (rs[c] && rs[c] !== '-') { refCol.push(c); k++; } else refGap[c] = 1; col2ref[c] = A1 + Math.max(0, k - 1); }
  refFirst = A1 + M.b0[0];
}
function colOfRef(p) { const i = p - refFirst; if (i < 0) return -1; if (i >= refCol.length) return NC; return refCol[i]; }
// per row of the MSA: graph node index of each column (-1 at gaps), for rows with a walk in graph gi
function colNode(gi, r) {
  const key = gi * 1000003 + r; if (colNodeCache.has(key)) return colNodeCache.get(key);
  const G = GR[gi], wi = G.g.row_walk[r], s = rowStr(r); let out = null;
  if (wi !== undefined && wi >= 0 && s) {
    const w = G.walks[wi], L = G.N.len, st = w.steps; out = new Int32Array(NC).fill(-1);
    let si = 0, b = M.b0[r];
    while (si < st.length && b >= L[Math.abs(st[si]) - 1]) { b -= L[Math.abs(st[si]) - 1]; si++; }
    let rem = si < st.length ? L[Math.abs(st[si]) - 1] - b : 0;
    for (let c = 0; c < NC; c++) { if (s[c] === '-') continue;
      while (rem === 0 && si < st.length - 1) { si++; rem = L[Math.abs(st[si]) - 1]; }
      if (si >= st.length || rem === 0) break;
      out[c] = Math.abs(st[si]) - 1; rem--; }
  }
  colNodeCache.set(key, out); return out;
}
// graph/MSA split for graph gi: same base in a column, different node
function splitOf(gi) {
  if (splitCache.has(gi)) return splitCache.get(gi);
  const wr = [], cn = [], ss = [], cells = ROWS.map(() => null);
  for (let i = 0; i < NR; i++) { const x = colNode(gi, i); if (x) { wr.push(i); cn.push(x); ss.push(rowStr(i)); cells[i] = new Float32Array(NC); } }
  const rowSame = new Float64Array(NR), rowDiff = new Float64Array(NR);
  const wb = new Map(), wbn = new Map(), sq = new Map(), sqn = new Map(); let sameTot = 0, diffTot = 0;
  for (let c = 0; c < NC; c++) { wb.clear(); wbn.clear(); sq.clear(); sqn.clear();
    for (let k = 0; k < wr.length; k++) { const ch = ss[k][c]; if (ch === '-') continue; const w = ROWS[wr[k]].weight || 1, key = ch + cn[k][c];
      wb.set(ch, (wb.get(ch) || 0) + w); wbn.set(key, (wbn.get(key) || 0) + w); sq.set(ch, (sq.get(ch) || 0) + w * w); sqn.set(key, (sqn.get(key) || 0) + w * w); }
    if (!wb.size) continue;
    let same = 0, samen = 0; for (const [x, v] of wb) same += (v * v - sq.get(x)) / 2; for (const [x, v] of wbn) samen += (v * v - sqn.get(x)) / 2;
    sameTot += same; diffTot += same - samen;
    for (let k = 0; k < wr.length; k++) { const ch = ss[k][c]; if (ch === '-') continue; const i = wr[k], w = ROWS[i].weight || 1, key = ch + cn[k][c];
      const sm = wb.get(ch) - w; cells[i][c] = sm > 0 ? (wb.get(ch) - wbn.get(key)) / sm : 0; rowSame[i] += sm; rowDiff[i] += wb.get(ch) - wbn.get(key); } }
  const res = {cells, frac: sameTot ? diffTot / sameTot : 0, rowSame, rowDiff}; splitCache.set(gi, res); return res;
}
function rowSplit(gi, r) { const s = splitOf(gi); return s.rowSame[r] ? s.rowDiff[r] / s.rowSame[r] : (colNode(gi, r) ? 0 : null); }
function nodeHue(i) { return (i * 137.50776) % 360; }
function nodeColor(G, i) { if (i < 0) return null; if (G.N.ref[i]) return (i % 2) ? COL.ref1 : COL.ref2;
  return 'hsl(' + nodeHue(i).toFixed(0) + ',' + (COL.dark ? '62%' : '68%') + ',' + (COL.dark ? '60%' : '50%') + ')'; }
function shareColor(w) { const b = D.share_bins, c = D.share_cols; for (let k = 0; k < b.length; k++) { if (w >= b[k][0] && (b[k][1] === null || w <= b[k][1])) return c[k]; } return c[c.length - 1]; }
function splitColor(v) { if (v <= 0) return COL.match; const t = Math.min(1, v); return 'rgba(224,49,49,' + (0.25 + 0.75 * t).toFixed(2) + ')'; }
function shareLegend() { return D.share_bins.map((b, k) => sw(D.share_cols[k]) + (b[1] === null ? b[0] + '+' : (b[0] === b[1] ? b[0] : b[0] + '-' + b[1]))).join(' '); }

// ------------------------------------------------------------------ header
function header() {
  const m = D.meta || {}, loc = D.census.locus, c = D.counts, rs = D.reads_summary || {}, bm = D.graphs[0].metrics, B = D.boundaries || {};
  let h = '<h1>' + (D.region_id ? esc(D.region_id) + ' <span class="sub">' : '<span>') + esc(D.contig) + ':' + fmt(D.req[0]) + '-' + fmt(D.req[1]) + (m.stratum ? ' · ' + esc(m.stratum) : '') + '</span></h1>';
  h += '<div class="line muted">span ' + esc(D.contig) + ':' + fmt(D.span[0]) + '-' + fmt(D.span[1]) + ' (' + fmt(D.span[1] - D.span[0] + 1) + ' bp of CHM13, anchor node to anchor node: ' +
       esc((B.L || {}).node || m.anchor_left || '?') + ' … ' + esc((B.R || {}).node || m.anchor_right || '?') + ')' + (D.package_dir ? ' · package <code>' + esc(D.package_dir) + '</code>' : '') + '</div>';
  let tr = '';
  if (m.tr_class) tr += '<span class="kv"><b>' + esc(m.tr_class) + '</b></span>';
  if (m.period) tr += '<span class="kv">period <b>' + esc(m.period) + '</b></span>';
  if (m.copies) tr += '<span class="kv">copies <b>' + esc(m.copies) + '</b></span>';
  if (m.region_len) tr += '<span class="kv">TR ' + esc(m.main_region_id || '') + ' <b>' + fmt(+m.region_len) + ' bp</b></span>';
  if (m.gc) tr += '<span class="kv">GC ' + esc(m.gc) + '</span>';
  if (m.motif) tr += '<span class="kv">motif <code title="' + esc(m.motif) + '">' + esc(m.motif.length > 48 ? m.motif.slice(0, 48) + '…' : m.motif) + '</code></span>';
  if (!tr && MAIN) tr += '<span class="kv"><b>' + esc(MAIN.cls) + '</b> ' + esc(MAIN.id) + ' period ' + esc(MAIN.period) + '</span>';
  if (!tr) tr = '<span class="kv muted">no tandem-repeat annotation</span>';
  if (m.vg_fp !== undefined || m.pg_fp !== undefined) tr += '<span class="kv">census: vg FP <b class="s-FP">' + esc(m.vg_fp) + '</b> FN <b class="s-FN">' + esc(m.vg_fn) + '</b> · PanGenie FP <b class="s-FP">' + esc(m.pg_fp) + '</b> FN <b class="s-FN">' + esc(m.pg_fn) + '</b></span>';
  if (m.error_status) tr += '<span class="kv">' + esc(m.error_status) + '</span>';
  h += '<div class="line">' + tr + '</div>';
  if (m.notes) h += '<div class="line muted">' + esc(Array.isArray(m.notes) ? m.notes.join('; ') : m.notes) + '</div>';
  const ts = D.tandem_scan || {};
  if (ts.vntr) h += '<div class="line muted">truth-independent scan of the CHM13 span: VNTR period ' + ts.vntr.period + ' over ' + fmt(ts.vntr.covered_bp) + ' bp (' + pct(ts.vntr.covered_frac) + ')</div>';
  let one = '<span class="kv">graphs <b>' + D.graphs.length + '</b> (' + D.graphs.map(g => esc(g.name)).join(', ') + ')</span>' +
    '<span class="kv">baseline nodes <b>' + fmt(bm.nodes) + '</b> (' + fmt(bm.nodes_per_kb) + '/kb)</span>' +
    '<span class="kv">end-to-end walks <b>' + fmt(bm.spanning_distinct_walks) + '</b> distinct, ' + fmt(bm.spanning_weight) + ' paths</span>';
  if (D.reads) one += '<span class="kv">reads <b>' + fmt(rs.n_alignments) + '</b>' + (rs.expected_from_truth_haplotypes ? ' (' + fmt(rs.expected_from_truth_haplotypes) + ' expected)' : '') + '</span><span class="kv">MAPQ&lt;5 <b>' + pct(rs.frac_mapq_lt5) + '</b></span>';
  if (c.source !== 'none') one += '<span class="kv">calls in span: vg FP <b class="s-FP">' + c.vg_fp + '</b> FN <b class="s-FN">' + c.vg_fn + '</b> · PanGenie FP <b class="s-FP">' + c.pg_fp + '</b> FN <b class="s-FN">' + c.pg_fn + '</b> (truth SVs ' + c.truth_sv + ')</span>';
  h += '<div class="line" style="margin-top:6px">' + one + '</div>';
  let l3 = '';
  for (const k of ['h1', 'h2']) { const x = (D.closest || {})[k]; if (x) l3 += '<span class="kv">closest panel walk to HG002 ' + k + ': <b>' + esc(x.walk) + '</b> edit ' + fmt(x.edit) + ' (' + pct(x.identity) + ')' + (x.exact_weight ? ', exact in ' + x.exact_weight + ' paths' : '') + '</span>'; }
  const chv = D.called_vs_truth;
  if (chv && chv.vg) l3 += '<span class="kv">edit to truth (h1+h2): CHM13 ' + fmt((chv.reference_as_call || {}).total) + ' · vg <b>' + fmt(chv.vg.total) + '</b> · PanGenie ' + fmt((chv.pangenie || {}).total) + '</span>';
  if (D.same_seq_groups.length) l3 += '<span class="kv">identical sequence on different baseline walks: <b>' + D.same_seq_groups.map(g => g.join('=')).join(', ') + '</b></span>';
  h += '<div class="line">' + l3 + '</div>';
  for (const w of D.warnings) h += '<div class="warn">' + esc(w) + '</div>';
  $('hdr').innerHTML = h;
  $('ftr').innerHTML = 'Built ' + esc(D.generated) + ' in ' + D.build_seconds + ' s by tools/viewer.py' + (D.package_dir ? ' from the region package <code>' + esc(D.package_dir) + '</code>' : '') +
    (D.toolkit_dir ? '; reads, calls and called haplotypes from region.py output <code>' + esc(D.toolkit_dir) + '</code>' : '; no region.py output (no reads, no called haplotypes)') +
    '; intermediates in <code>' + esc(D.view_dir) + '</code>. Coordinates are 1-based CHM13. Identity = 1 - unit-cost edit distance / longer length.';
}

// ------------------------------------------------------------------ shared MSA view state
const VIEW = {px: 10, rowH: 14, mode: 'diff', g: 0, cmp: 0, dots: true, ticks: 'ruler', hidden: new Set()};
const NODE = {mode: 'msa', color: 'node', rowH: 12, frag: false, zoom: 'fit'};
const scrollers = [];
let drawPending = false;
function requestDraw() { if (drawPending) return; drawPending = true; requestAnimationFrame(() => { drawPending = false; drawAllMsa(); }); }
function registerScroller(el, draw, active) { const s = {el, draw, active: active || (() => true)}; scrollers.push(s);
  el.addEventListener('scroll', () => { if (!s.active()) { draw(); return; }
    for (const o of scrollers) if (o !== s && o.active() && Math.abs(o.el.scrollLeft - el.scrollLeft) > 0.5) o.el.scrollLeft = el.scrollLeft;
    requestDraw(); }, {passive: true});
  let dragX = null, startLeft = 0;
  el.addEventListener('mousedown', e => { dragX = e.clientX; startLeft = el.scrollLeft; el.style.cursor = 'grabbing'; e.preventDefault(); });
  window.addEventListener('mouseup', () => { if (dragX !== null) { dragX = null; el.style.cursor = ''; } });
  window.addEventListener('mousemove', e => { if (dragX !== null) el.scrollLeft = startLeft - (e.clientX - dragX); });
}
function setSpacers() { for (const s of scrollers) if (s.active()) s.el.querySelector('.spacer').style.width = Math.max(1, Math.ceil(NC * VIEW.px)) + 'px'; }
function viewCol0(el) { return el.scrollLeft / VIEW.px; }
function drawAllMsa() { drawAlignment(); drawOverviewBox(); if (NODE.mode === 'msa') drawAllNodes(); drawPile(); }
function setZoom(px, centerCol) {
  const el = $('aln-scroll'); const w = el.clientWidth;
  const cc = centerCol !== undefined ? centerCol : (el.scrollLeft + w / 2) / VIEW.px;
  VIEW.px = px; setSpacers();
  const left = Math.max(0, cc * px - w / 2);
  for (const s of scrollers) if (s.active()) s.el.scrollLeft = left;
  drawAllMsa(); drawOverview();
}
function fitPx() { return Math.max(0.02, ($('aln-scroll').clientWidth - 2) / NC); }
function prepCanvas(cv, w, h) { const dpr = window.devicePixelRatio || 1; cv.style.width = w + 'px'; cv.style.height = h + 'px';
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
  const ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, w, h); return ctx; }

// ------------------------------------------------------------------ annotations (CHM13 coordinates -> columns)
let ANN = [];
function annItem(group, s, e, colorKey, label, tipHtml, thin) {
  let c0 = colOfRef(s), c1 = colOfRef(Math.max(s, e)); if (c1 < 0 || c0 >= NC) return; c0 = Math.max(0, c0); c1 = Math.min(NC - 1, c1);
  ANN.push({group, c0, c1, colorKey, label, tip: tipHtml, thin: !!thin}); }
function buildAnnotations() {
  ANN = [];
  for (const r of D.trs) annItem('TR', r.s, r.e, r.cls === 'VNTR' ? 'vntr' : (r.cls === 'STR' ? 'str' : 'oth'), (r.cls || '') + ' p' + r.period,
    '<b>' + esc(r.id) + '</b> ' + esc(r.cls) + ' ' + fmt(r.s) + '-' + fmt(r.e) + ' (' + fmt(r.e - r.s + 1) + ' bp)<br>period ' + esc(r.period) + ', copies ' + esc(r.copies) + (r.gc !== null && r.gc !== undefined ? ', GC ' + esc(r.gc) : '') + '<br>motif <code>' + esc((r.motif || '').slice(0, 120)) + '</code><br>source ' + esc(r.src) + (r.vg_fp !== null && r.vg_fp !== undefined ? '<br>census counts: vg FP ' + r.vg_fp + ', FN ' + r.vg_fn : ''));
  for (const t of D.truth) { if (t.source !== 'stvar') continue; if (!t.sv50 && Math.abs(t.dlen || 0) < 10) continue;
    const st = t.vg === 'TP' ? 'tp' : (t.vg === 'FN' ? 'fn' : 'na');
    annItem('truth SV', t.pos, t.end, st, (t.svtype || '') + ' ' + t.dlen, '<b>truth ' + esc(t.svtype || 'indel') + '</b> ' + fmt(t.pos) + '-' + fmt(t.end) + ' Δlen ' + fmt(t.dlen) + ' GT ' + esc(t.gt) + '<br>vg ' + esc(t.vg) + ' · PanGenie ' + esc(t.pg) + ' · benchmark ' + esc(t.bench) + '<br>applied h1 ' + esc(t.h1) + ', h2 ' + esc(t.h2), !t.sv50); }
  for (const c of D.calls) { if (!c.called) continue; const sv = Math.abs(c.len_diff || 0) >= 50; if (!sv && Math.abs(c.len_diff || 0) < 10) continue;
    const st = c.truvari === 'TP' ? 'tp' : (c.truvari === 'FP' ? 'fp' : 'na'); const g = c.caller === 'vg' ? 'vg calls' : 'PanGenie calls';
    annItem(g, c.pos, c.end, st, 'Δ' + c.len_diff, '<b>' + esc(c.caller) + '</b> ' + fmt(c.pos) + '-' + fmt(c.end) + ' Δlen ' + fmt(c.len_diff) + ' GT ' + esc(c.gt) + '<br>truvari <b>' + esc(c.truvari) + '</b>' + (c.pctseq !== null && c.pctseq !== undefined ? ' · PctSeqSim ' + fmt(c.pctseq) + ' · PctSizeSim ' + fmt(c.pctsize) : '') + '<br>id ' + esc(c.id), !sv); }
  const groups = ['TR', 'truth SV', 'vg calls', 'PanGenie calls'];
  ANN.lanes = {};
  for (const g of groups) { const items = ANN.filter(a => a.group === g).sort((a, b) => a.c0 - b.c0); const ends = [];
    for (const a of items) { let l = ends.findIndex(e => e < a.c0 - 2); if (l < 0) { l = ends.length; ends.push(-1); } ends[l] = a.c1; a.lane = l; }
    ANN.lanes[g] = Math.max(1, Math.min(ends.length, 8)); }
  ANN.groups = groups;
}
const LANE_H = 7, RULER_H = 20;
function annLayout() { let y = RULER_H; const lay = []; for (const g of ANN.groups) { const h = ANN.lanes[g] * LANE_H + 4; lay.push({g, y, h}); y += h; } return {lay, h: y + 4}; }

// ------------------------------------------------------------------ alignment panel
let rowOrder = ROWS.map((_, i) => i);
function visibleRows() { return rowOrder.filter(i => !VIEW.hidden.has(i)); }
function rowLabelHtml(i) { const r = ROWS[i]; let s = '<span class="k-' + r.kind + '">' + esc(r.label) + '</span>';
  if (r.weight !== undefined) s += ' <span class="badge" title="paths with this walk in the baseline graph">×' + r.weight + '</span>';
  if (r.closest_to && r.closest_to.length) s += ' <span class="badge" title="closest panel walk to truth ' + r.closest_to.join(',') + '">★' + r.closest_to.join(',') + '</span>';
  if (r.same_seq && r.same_seq.length) s += ' <span class="badge" title="same sequence as ' + esc(r.same_seq.join(',')) + ' but a different baseline node path">=' + esc(r.same_seq.join(',')) + '</span>';
  if (M && M.a[i] < 0) s += ' <span class="badge" title="this row is not in the MSA on show">absent</span>';
  s += ' <span class="muted">' + fmt(r.len) + '</span>';
  return s; }
function rowTip(i) { const r = ROWS[i]; let s = '<b>' + esc(r.label) + '</b> (' + esc(r.kind) + ')<br>length ' + fmt(r.len) + ' bp';
  if (r.walk !== undefined) { s += ', ' + fmt(r.steps) + ' baseline node steps, ' + r.weight + ' paths<br>non-reference nodes ' + fmt(r.nonref_nodes) + ' (' + fmt(r.nonref_bp) + ' bp)' +
    '<br>edit to CHM13 ' + fmt(r.edit_ref) + ' · to HG002 h1 ' + fmt(r.edit_h1) + ' (' + pct(r.id_h1) + ') · h2 ' + fmt(r.edit_h2) + ' (' + pct(r.id_h2) + ')';
    s += '<br>graph/MSA split (' + esc(D.graphs[VIEW.g].name) + '): ' + pct(rowSplit(VIEW.g, i));
    if (r.reads !== undefined && D.reads) s += '<br>reads following this walk ' + r.reads + ' (only this walk: ' + r.reads_unique + ')';
    if (r.names && r.names.length) s += '<br>paths: ' + esc(r.names.join(', ')); }
  if (r.kind === 'truth') s += '<br>edit to CHM13 ' + fmt(r.edit_ref) + '<br>closest panel walk ' + esc(r.closest) + ' (edit ' + fmt(r.closest_edit) + ', ' + pct(r.closest_id) + ')<br>closest non-reference walk ' + esc(r.closest_nonref) + ' (edit ' + fmt(r.closest_nonref_edit) + ')<br><i>evaluation only</i>';
  if (r.kind === 'called') s += '<br>' + esc(r.label.split(' ')[0]) + ' genotype applied to CHM13<br>paired with truth ' + esc(r.paired_truth) + ', edit ' + fmt(r.edit_truth) + (r.phase_reliable === false ? ' (phase unreliable)' : '');
  if (M.a[i] < 0) s += '<br><b>not in the MSA on show</b>';
  s += '<br><span class="muted">click: compare every row against this one</span>';
  return s; }
function buildLabels() {
  const el = $('aln-labels'); const {lay, h} = annLayout(); let s = '<div class="hd" style="height:' + RULER_H + 'px;line-height:' + RULER_H + 'px">CHM13 ' + esc(D.contig) + (PERIOD ? ' · period ' + PERIOD : '') + '</div>';
  for (const L of lay) s += '<div class="hd" style="height:' + L.h + 'px;line-height:' + Math.min(L.h, 14) + 'px">' + esc(L.g) + '</div>';
  s += '<div style="height:' + (h - lay.reduce((a, L) => a + L.h, RULER_H)) + 'px"></div>';
  for (const i of visibleRows()) s += '<div class="lab' + (VIEW.cmp === i ? ' sel' : '') + '" data-i="' + i + '" style="height:' + VIEW.rowH + 'px;line-height:' + VIEW.rowH + 'px;font-size:' + Math.min(12, VIEW.rowH - 2) + 'px">' + rowLabelHtml(i) + '<span class="x" title="hide row">×</span></div>';
  el.innerHTML = s;
  el.querySelectorAll('.lab').forEach(d => { const i = +d.dataset.i;
    d.addEventListener('mousemove', e => showTip(e, rowTip(i))); d.addEventListener('mouseleave', hideTip);
    d.querySelector('.x').addEventListener('click', e => { e.stopPropagation(); VIEW.hidden.add(i); refreshRows(); });
    d.addEventListener('click', () => { VIEW.cmp = i; $('cmpsel').value = i; if (VIEW.mode !== 'diff') { VIEW.mode = 'diff'; $('modesel').value = 'diff'; legendAln(); } refreshRows(); }); });
}
function refreshRows() { buildLabels(); buildAllNodeLabels(); drawAllMsa(); drawOverview(); updateHiddenNote(); if (NODE.mode !== 'msa') drawAllNodes(); }
function updateHiddenNote() { const n = VIEW.hidden.size; $('hiddennote').textContent = n ? n + ' hidden' : ''; $('showall').style.display = n ? '' : 'none'; }
function cmpString() { return rowStr(VIEW.cmp) || rowStr(0); }
// 0 empty, 1 match, 2 mismatch, 3 deletion (row gap), 4 insertion (row base, comparison gap)
function cellState(s, c, cmpS) { if (!s) return 0; const a = s[c], b = cmpS[c]; if (a === '-') return b === '-' ? 0 : 3; if (b === '-') return 4; return a === b ? 1 : 2; }
const PRI = [0, 1, 4, 2, 3];
function stateColor(st) { return st === 1 ? COL.match : st === 2 ? COL.mismatch : st === 4 ? COL.ins : null; }
function rowColoring(r) { if (VIEW.mode === 'diff') return null; const cn = colNode(VIEW.g, r); if (!cn) return null;
  return {cn, G: GR[VIEW.g], sc: VIEW.mode === 'split' ? splitOf(VIEW.g).cells[r] : null}; }
function cellColor(rc, c, st) {
  if (!rc) return stateColor(st); const n = rc.cn[c]; if (n < 0) return null;
  if (VIEW.mode === 'node') return nodeColor(rc.G, n);
  if (VIEW.mode === 'share') return shareColor(rc.G.N.w[n]);
  if (VIEW.mode === 'split') return splitColor(rc.sc[c]);
  return stateColor(st); }
function drawAlignment() {
  const el = $('aln-scroll'), cv = $('aln-canvas'); const W = el.clientWidth; const {lay, h: annH} = annLayout(); const vr = visibleRows();
  const H = annH + vr.length * VIEW.rowH + 2; const ctx = prepCanvas(cv, W, H); const px = VIEW.px, col0 = viewCol0(el);
  const cStart = Math.max(0, Math.floor(col0)), cEnd = Math.min(NC, Math.ceil(col0 + W / px) + 1);
  const X = (c) => (c - col0) * px;
  ctx.fillStyle = COL.muted; ctx.font = '10px ui-monospace,Menlo,monospace'; ctx.textBaseline = 'top';
  const steps = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000];
  const tickEvery = steps.find(s => s * px >= 70) || 100000;
  ctx.strokeStyle = COL.border; ctx.beginPath(); ctx.moveTo(0, RULER_H - 0.5); ctx.lineTo(W, RULER_H - 0.5); ctx.stroke();
  const p0 = col2ref[cStart] || A1, p1 = col2ref[Math.max(0, cEnd - 1)] || A1;
  for (let p = Math.ceil(p0 / tickEvery) * tickEvery; p <= p1; p += tickEvery) { const c = colOfRef(p); if (c < 0 || c >= NC) continue; const x = X(c) + px / 2;
    ctx.fillStyle = COL.muted; ctx.fillRect(x, RULER_H - 6, 1, 6); ctx.fillText(p.toLocaleString('en-US'), x + 2, 1); }
  { const a = colOfRef(D.req[0]), b = colOfRef(D.req[1]); ctx.fillStyle = COL.sel; ctx.globalAlpha = 0.5; ctx.fillRect(X(Math.max(0, a)), RULER_H - 3, (Math.min(NC - 1, b) - Math.max(0, a) + 1) * px, 3); ctx.globalAlpha = 1; }
  if (PERIOD && MAIN && VIEW.ticks !== 'off') { ctx.fillStyle = COL.vntr; const tStep = PERIOD * px < 3 ? Math.ceil(3 / (PERIOD * px)) : 1;
    for (let p = MAIN.s, k = 0; p <= MAIN.e; p += PERIOD, k++) { if (k % tStep) continue; const c = colOfRef(p); if (c < cStart - 1 || c > cEnd) continue; ctx.fillRect(X(c), RULER_H - 11, 1, 5); } }
  for (const L of lay) { ctx.fillStyle = COL.grid; ctx.fillRect(0, L.y + L.h - 1, W, 1);
    for (const a of ANN) { if (a.group !== L.g || a.lane >= 8) continue; if (a.c1 < cStart - 1 || a.c0 > cEnd) continue;
      const y = L.y + 2 + a.lane * LANE_H; ctx.fillStyle = COL[a.colorKey] || COL.na; const x0 = X(a.c0), w = Math.max(3, (a.c1 - a.c0 + 1) * px);
      ctx.globalAlpha = a.thin ? 0.55 : 1; ctx.fillRect(x0, y + (a.thin ? 2 : 0), w, a.thin ? LANE_H - 4 : LANE_H - 1); ctx.globalAlpha = 1; } }
  const cmpS = cmpString(); const rh = VIEW.rowH; const text = px >= 7 && rh >= 10;
  ctx.font = Math.min(12, Math.floor(Math.min(px * 1.25, rh - 1))) + 'px ui-monospace,SFMono-Regular,Menlo,monospace'; ctx.textBaseline = 'middle'; ctx.textAlign = 'center';
  vr.forEach((r, k) => { const y = annH + k * rh; const s = rowStr(r); const rc = rowColoring(r);
    if (k % 2) { ctx.fillStyle = COL.grid; ctx.globalAlpha = 0.35; ctx.fillRect(0, y, W, rh); ctx.globalAlpha = 1; }
    if (!s) { ctx.fillStyle = COL.absent; ctx.fillRect(0, y + 1, W, rh - 2); ctx.fillStyle = COL.muted; ctx.textAlign = 'start'; ctx.font = '10px sans-serif'; ctx.fillText('not in this MSA', 4, y + rh / 2); ctx.textAlign = 'center';
      ctx.font = Math.min(12, Math.floor(Math.min(px * 1.25, rh - 1))) + 'px ui-monospace,SFMono-Regular,Menlo,monospace'; return; }
    if (px >= 1) {
      let runC = null, runX = 0, runW = 0;
      const flush = () => { if (runC) { ctx.fillStyle = runC; ctx.fillRect(runX, y + 1, runW, rh - 2); } runC = null; };
      for (let c = cStart; c < cEnd; c++) { const st = cellState(s, c, cmpS); const col = cellColor(rc, c, st); const x = X(c);
        if (col === runC && col) { runW = x + px - runX; continue; } flush(); if (col) { runC = col; runX = x; runW = px; }
        if (st === 3) { ctx.fillStyle = COL.delline; ctx.fillRect(x, y + rh / 2 - 0.5, px, 1); } }
      flush();
      if (text) { for (let c = cStart; c < cEnd; c++) { const ch = s[c]; if (ch === '-') continue; const st = cellState(s, c, cmpS);
          if (st === 1 && VIEW.mode === 'diff') { ctx.fillStyle = COL.matchtext; ctx.fillText(VIEW.dots && r !== VIEW.cmp ? '·' : ch, X(c) + px / 2, y + rh / 2 + 0.5); }
          else { ctx.fillStyle = (VIEW.mode === 'diff' && (st === 2 || st === 4)) ? '#ffffff' : COL.fg; ctx.fillText(ch, X(c) + px / 2, y + rh / 2 + 0.5); } } }
    } else {
      for (let x = 0; x < W; x++) { const c0 = Math.floor(col0 + x / px), c1 = Math.min(NC, Math.floor(col0 + (x + 1) / px)); if (c0 >= NC) break;
        let best = 0, bc = c0; for (let c = c0; c < Math.max(c0 + 1, c1); c++) { const st = cellState(s, c, cmpS); if (PRI[st] > PRI[best]) { best = st; bc = c; } }
        const col = cellColor(rc, bc, best); if (col) { ctx.fillStyle = col; ctx.fillRect(x, y + 1, 1, rh - 2); }
        else if (best === 3) { ctx.fillStyle = COL.delline; ctx.fillRect(x, y + rh / 2 - 0.5, 1, 1); } }
    }
    if (VIEW.ticks === 'rows' && PERIOD && MAIN) drawRowTicks(ctx, r, y, X, cStart, cEnd);
    if (r === VIEW.cmp) { ctx.strokeStyle = COL.sel; ctx.lineWidth = 1; ctx.strokeRect(0.5, y + 0.5, W - 1, rh - 1); }
  });
  ctx.textAlign = 'start';
}
function drawRowTicks(ctx, r, y, X, cStart, cEnd) {
  const a = Math.max(0, colOfRef(MAIN.s)), b = Math.min(NC - 1, colOfRef(MAIN.e)); if (b < a) return;
  let cols = tickCache.get(r); if (!cols) { cols = []; const s = rowStr(r) || ''; let k = 0; for (let c = a; c <= b; c++) { if (s[c] === '-') continue; if (k % PERIOD === 0) cols.push(c); k++; } tickCache.set(r, cols); }
  ctx.fillStyle = COL.tick; const tStep = PERIOD * VIEW.px < 3 ? Math.ceil(3 / (PERIOD * VIEW.px)) : 1;
  cols.forEach((c, k) => { if (k % tStep || c < cStart - 1 || c > cEnd) return; ctx.fillRect(X(c), y, 1, 4); });
}
function alnHover(e) {
  const el = $('aln-scroll'), rect = $('aln-canvas').getBoundingClientRect(); const x = e.clientX - rect.left, y = e.clientY - rect.top;
  const c = Math.floor(viewCol0(el) + x / VIEW.px); if (c < 0 || c >= NC) { hideTip(); return; }
  const {lay, h: annH} = annLayout(); let s = 'column ' + (c + 1) + ' · CHM13 ' + D.contig + ':' + col2ref[c].toLocaleString('en-US') + (refGap[c] ? ' (insertion relative to CHM13)' : '');
  if (y >= annH) { const k = Math.floor((y - annH) / VIEW.rowH); const vr = visibleRows(); if (k >= vr.length) { hideTip(); return; } const r = vr[k]; const rs = rowStr(r); const ch = rs ? rs[c] : null;
    s = '<b>' + esc(ROWS[r].label) + '</b> ' + (ch === null ? 'not in this MSA' : ch === '-' ? 'gap' : 'base <b>' + ch + '</b>') + '<br>' + s;
    GR.forEach((G, gi) => { const cn = colNode(gi, r); if (cn && cn[c] >= 0) { const n = cn[c];
      s += '<br>' + esc(G.g.name) + ': node ' + esc(G.N.id[n]) + ' (' + G.N.len[n] + ' bp, ' + (G.N.ref[n] ? 'on CHM13' : 'off CHM13') + ', ' + G.N.w[n] + ' paths)' + (gi === VIEW.g ? '; same base on another node here ' + pct(splitOf(gi).cells[r][c]) : ''); } });
  } else if (y >= RULER_H) { const L = lay.find(L => y >= L.y && y < L.y + L.h); if (L) { const lane = Math.floor((y - L.y - 2) / LANE_H);
      const hits = ANN.filter(a => a.group === L.g && a.lane === lane && c >= a.c0 - 2 / VIEW.px - 1 && c <= a.c1 + 2 / VIEW.px + 1); if (hits.length) s = hits.map(a => a.tip).join('<hr>') + '<br><span class="muted">' + s + '</span>'; } }
  showTip(e, s);
}
function msaSub() {
  $('aln-sub').innerHTML = esc(M.name) + ': ' + ROWS.length + ' rows, ' + NC.toLocaleString('en-US') + ' columns' + (M.crop ? ' (of ' + M.ncols_full.toLocaleString('en-US') + ')' : '') + ' · ' + esc(M.note || '') + '; hover for bases and nodes, click a label to compare against it, drag to pan';
}
function alnControls() {
  const opts = ROWS.map((r, i) => '<option value="' + i + '">' + esc(r.label) + '</option>').join('');
  const mopts = D.msas.map((m, i) => '<option value="' + i + '">' + esc(m.name) + '</option>').join('');
  const gopts = D.graphs.map((g, i) => '<option value="' + i + '">' + esc(g.name) + '</option>').join('');
  $('aln-controls').innerHTML =
    '<label>alignment <select id="msasel">' + mopts + '</select></label>' +
    '<label>zoom <select id="zoomsel"><option value="12">bases 12px</option><option value="9">bases 9px</option><option value="4">4 px/col</option><option value="2">2 px/col</option><option value="1">1 px/col</option><option value="0.5">0.5 px/col</option><option value="fit">fit whole MSA</option></select></label>' +
    '<button id="zin" type="button" title="zoom in">+</button><button id="zout" type="button" title="zoom out">−</button>' +
    '<label>colour <select id="modesel"><option value="diff">differences vs comparison row</option><option value="node">graph node (same node = same colour)</option><option value="share">node sharing (paths through node)</option><option value="split">graph/MSA split (same base, other node)</option></select></label>' +
    '<label>graph <select id="gsel">' + gopts + '</select></label>' +
    '<label>compare to <select id="cmpsel">' + opts + '</select></label>' +
    '<label><input type="checkbox" id="dots" checked> matches as dots</label>' +
    '<label>period ticks <select id="ticksel"><option value="ruler">ruler</option><option value="rows">every row</option><option value="off">off</option></select></label>' +
    '<label>rows <select id="rowh"><option value="14">normal</option><option value="8">compact</option><option value="20">tall</option></select></label>' +
    '<label>sort <select id="sortsel"><option value="default">truth & calls first, walks by weight</option><option value="len">by length</option><option value="h1">by edit to HG002 h1</option><option value="h2">by edit to HG002 h2</option><option value="split">by graph/MSA split</option><option value="msa">by MSA similarity</option></select></label>' +
    '<label><input type="checkbox" id="g-called" checked> calls</label><label><input type="checkbox" id="g-truth" checked> truth</label><label><input type="checkbox" id="g-panel" checked> panel walks</label>' +
    '<span id="hiddennote" class="muted"></span><button id="showall" type="button" style="display:none">show all rows</button>' +
    '<label>go to CHM13 <input id="gopos" size="11" placeholder="position"></label>';
  $('zoomsel').value = '12';
  $('msasel').addEventListener('change', e => switchMsa(+e.target.value));
  $('zoomsel').addEventListener('change', e => { const v = e.target.value; setZoom(v === 'fit' ? fitPx() : +v); });
  $('zin').addEventListener('click', () => setZoom(Math.min(16, VIEW.px * 2)));
  $('zout').addEventListener('click', () => setZoom(Math.max(fitPx(), VIEW.px / 2)));
  $('modesel').addEventListener('change', e => { VIEW.mode = e.target.value; legendAln(); drawAlignment(); drawOverview(); });
  $('gsel').addEventListener('change', e => { VIEW.g = +e.target.value; legendAln(); drawAlignment(); drawOverview(); });
  $('cmpsel').addEventListener('change', e => { VIEW.cmp = +e.target.value; refreshRows(); });
  $('dots').addEventListener('change', e => { VIEW.dots = e.target.checked; drawAlignment(); });
  $('ticksel').addEventListener('change', e => { VIEW.ticks = e.target.value; drawAlignment(); });
  $('rowh').addEventListener('change', e => { VIEW.rowH = +e.target.value; refreshRows(); });
  $('sortsel').addEventListener('change', e => { sortRows(e.target.value); refreshRows(); });
  for (const [id, kinds] of [['g-called', ['called']], ['g-truth', ['truth']], ['g-panel', ['panel']]])
    $(id).addEventListener('change', e => { ROWS.forEach((r, i) => { if (kinds.includes(r.kind)) { if (e.target.checked) VIEW.hidden.delete(i); else VIEW.hidden.add(i); } }); refreshRows(); });
  $('showall').addEventListener('click', () => { VIEW.hidden.clear(); ['g-called', 'g-truth', 'g-panel'].forEach(id => $(id).checked = true); refreshRows(); });
  $('gopos').addEventListener('keydown', e => { if (e.key !== 'Enter') return; const p = parseInt(e.target.value.replace(/[, ]/g, ''), 10); if (!p) return; const c = colOfRef(p); if (c < 0 || c >= NC) { e.target.style.borderColor = 'var(--fp)'; return; } e.target.style.borderColor = ''; setZoom(Math.max(VIEW.px, 9), c); });
}
function switchMsa(i) {
  const el = $('aln-scroll'); const cc = Math.max(0, Math.min(NC - 1, Math.floor((el.scrollLeft + el.clientWidth / 2) / VIEW.px))); const refPos = col2ref[cc];
  const wasFit = Math.abs(VIEW.px - fitPx()) < 1e-6;
  setMsa(i); buildAnnotations(); placeReads(); layoutPile(); msaSub(); metricsTable();
  buildLabels(); buildAllNodeLabels(); setZoom(wasFit ? fitPx() : VIEW.px, colOfRef(refPos)); drawAllNodes();
}
function legendAln() {
  let s; const gname = esc(D.graphs[VIEW.g].name);
  if (VIEW.mode === 'diff') s = sw(COL.match) + 'same as comparison row' + ' &nbsp;' + sw(COL.mismatch) + 'mismatch' + ' &nbsp;' + sw(COL.ins) + 'base where comparison row has a gap' + ' &nbsp;<span class="sw" style="background:linear-gradient(transparent 45%,' + COL.delline + ' 45%,' + COL.delline + ' 60%,transparent 60%)"></span>gap where comparison row has a base';
  else if (VIEW.mode === 'node') s = gname + ': ' + sw(COL.ref1) + sw(COL.ref2) + 'node on the CHM13 path (alternating) &nbsp;' + sw('hsl(30,68%,50%)') + sw('hsl(200,68%,50%)') + 'other node, colour hashed from the node; same colour in one column = same node. Rows with no path in this graph (truth, calls) stay in difference colours.';
  else if (VIEW.mode === 'share') s = gname + ': ' + shareLegend() + ' paths through the node';
  else s = gname + ': ' + sw(COL.match) + 'every row with this base in this column is on the same node &nbsp;' + sw('rgba(224,49,49,0.5)') + sw('rgba(224,49,49,1)') + 'share of same-base rows (by paths) on a different node: unmerged copies, or an indel the MSA and the graph place at different copies of a repeat.';
  $('aln-legend').innerHTML = s + (PERIOD ? ' &nbsp;' + sw(COL.vntr) + 'period-' + PERIOD + ' ticks from ' + esc(MAIN.id) : '');
}
function sortRows(how) {
  const base = ROWS.map((_, i) => i);
  const kindRank = {ref: 0, truth: 1, called: 2, panel: 3};
  if (how === 'msa') { rowOrder = msaOrder(); return; }
  const key = {
    default: i => [kindRank[ROWS[i].kind], 0],
    len: i => [ROWS[i].kind === 'ref' ? -1 : 0, -(ROWS[i].len || 0)],
    h1: i => [kindRank[ROWS[i].kind] < 3 ? kindRank[ROWS[i].kind] - 5 : 0, ROWS[i].edit_h1 === undefined || ROWS[i].edit_h1 === null ? 1e12 : ROWS[i].edit_h1],
    h2: i => [kindRank[ROWS[i].kind] < 3 ? kindRank[ROWS[i].kind] - 5 : 0, ROWS[i].edit_h2 === undefined || ROWS[i].edit_h2 === null ? 1e12 : ROWS[i].edit_h2],
    split: i => [kindRank[ROWS[i].kind] < 3 ? kindRank[ROWS[i].kind] - 5 : 0, -(rowSplit(VIEW.g, i) || 0)],
  }[how];
  rowOrder = base.slice().sort((a, b) => { const ka = key(a), kb = key(b); return ka[0] - kb[0] || ka[1] - kb[1] || a - b; });
}
function msaOrder() {
  const dist = (i, j) => { const a = rowStr(i), b = rowStr(j); if (!a || !b) return 1e12; if (M.a[i] === M.a[j]) return 0; let d = 0; for (let c = 0; c < NC; c++) if (a[c] !== b[c]) d++; return d; };
  const left = new Set(ROWS.map((_, i) => i)); left.delete(0); const out = [0]; let cur = 0;
  while (left.size) { let best = -1, bd = Infinity; for (const j of left) { const d = dist(cur, j); if (d < bd) { bd = d; best = j; } } out.push(best); left.delete(best); cur = best; }
  return out; }
let ovImg = null;
function drawOverview() {
  const cv = $('overview'); const W = cv.clientWidth || 800; const vr = visibleRows(); const rh = Math.max(1, Math.min(3, Math.floor(46 / Math.max(1, vr.length)))); const H = 12 + vr.length * rh + 2;
  cv.style.height = H + 'px'; const ctx = prepCanvas(cv, W, H); const cmpS = cmpString(); const per = NC / W;
  for (const a of ANN) { if (a.thin || (a.group !== 'truth SV' && a.group !== 'vg calls')) continue; ctx.fillStyle = COL[a.colorKey]; const y = a.group === 'truth SV' ? 1 : 6; ctx.fillRect(a.c0 / per, y, Math.max(1, (a.c1 - a.c0 + 1) / per), 4); }
  vr.forEach((r, k) => { const y = 12 + k * rh; const s = rowStr(r); if (!s) return; const rc = rowColoring(r);
    for (let x = 0; x < W; x++) { const c0 = Math.floor(x * per), c1 = Math.max(c0 + 1, Math.floor((x + 1) * per)); let best = 0, bc = c0;
      for (let c = c0; c < c1 && c < NC; c++) { const st = cellState(s, c, cmpS); if (PRI[st] > PRI[best]) { best = st; bc = c; } }
      const col = cellColor(rc, bc, best); if (col) { ctx.fillStyle = col; ctx.fillRect(x, y, 1, rh); } } });
  ovImg = ctx.getImageData(0, 0, cv.width, cv.height); drawOverviewBox();
}
function drawOverviewBox() { const cv = $('overview'); if (!ovImg) return; const ctx = cv.getContext('2d'); const dpr = window.devicePixelRatio || 1; ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.putImageData(ovImg, 0, 0); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const W = cv.clientWidth, H = cv.clientHeight; const el = $('aln-scroll'); const c0 = viewCol0(el), cw = el.clientWidth / VIEW.px;
  ctx.strokeStyle = COL.sel; ctx.lineWidth = 2; ctx.strokeRect(Math.max(1, c0 / NC * W), 1, Math.max(3, Math.min(W, cw / NC * W)), H - 2);
  const a = colOfRef(D.req[0]), b = colOfRef(D.req[1]); ctx.fillStyle = COL.sel; ctx.globalAlpha = 0.35; ctx.fillRect(Math.max(0, a) / NC * W, H - 3, Math.max(2, (Math.min(NC - 1, b) - Math.max(0, a)) / NC * W), 3); ctx.globalAlpha = 1; }
function overviewNav(e) { const cv = $('overview'), rect = cv.getBoundingClientRect(); const c = (e.clientX - rect.left) / rect.width * NC; const el = $('aln-scroll'); el.scrollLeft = Math.max(0, c * VIEW.px - el.clientWidth / 2); }

// ------------------------------------------------------------------ generic table
function renderTable(el, rows, cols, opts) {
  opts = opts || {}; let sortK = opts.sort || null, sortDir = opts.dir || 1; const max = opts.max || 2000;
  function draw() {
    let rs = rows.slice(); if (sortK !== null) { const col = cols.find(c => c[0] === sortK); rs.sort((a, b) => { let x = a[sortK], y = b[sortK];
      if (x === null || x === undefined || x === '' || x === '.') return 1; if (y === null || y === undefined || y === '' || y === '.') return -1;
      if (col && col[2] === 'n') return sortDir * (+x - +y); return sortDir * String(x).localeCompare(String(y)); }); }
    let h = '<div class="tablewrap"><table><thead><tr>' + cols.map(c => '<th class="' + (c[2] || '') + '" data-k="' + c[0] + '">' + esc(c[1]) + (sortK === c[0] ? (sortDir > 0 ? ' ▲' : ' ▼') : '') + '</th>').join('') + '</tr></thead><tbody>';
    for (const r of rs.slice(0, max)) { h += '<tr>' + cols.map(c => { const v = r[c[0]]; const f = c[3] ? c[3](v, r) : (typeof v === 'number' ? fmt(v) : esc(v === undefined || v === null ? '' : v));
      const cls = (c[2] || '') + (/^(TP|FP|FN)$/.test(String(v)) ? ' s-' + v : ''); return '<td class="' + cls + '">' + f + '</td>'; }).join('') + '</tr>'; }
    h += '</tbody></table></div>'; if (rs.length > max) h += '<div class="muted">showing ' + max + ' of ' + rs.length + ' rows</div>'; if (!rs.length) h += '<div class="muted">no rows</div>';
    el.innerHTML = h; el.querySelectorAll('th').forEach(th => th.addEventListener('click', () => { const k = th.dataset.k; if (sortK === k) sortDir = -sortDir; else { sortK = k; sortDir = 1; } draw(); }));
  }
  draw(); return {update: (r) => { rows = r; draw(); }};
}
function walkTable() {
  const B = GR[0]; const rows = D.graphs[0].walks.map((w, i) => { const ri = rowOfBaseWalk.get(i); const r = ri !== undefined ? ROWS[ri] : {}; const pc = (D.path_counts || [])[i] || {};
    return {id: w.id, cls: w.cls, weight: w.weight, len: w.len, steps: B.walks[i].steps.length, nonref_bp: r.nonref_bp, edit_ref: r.edit_ref, edit_h1: r.edit_h1, id_h1: r.id_h1, edit_h2: r.edit_h2, id_h2: r.id_h2,
      reads: D.reads ? pc.compat : null, unique: D.reads ? pc.unique : null, same: (r.same_seq || []).join(','), closest: (r.closest_to || []).join(','), names: (w.names || []).slice(0, 8).join(', ')}; });
  renderTable($('walk-table'), rows, [
    ['id', 'walk'], ['cls', 'class'], ['weight', 'paths', 'n'], ['len', 'bp', 'n'], ['steps', 'node steps', 'n'], ['nonref_bp', 'non-ref bp', 'n'], ['edit_ref', 'edit CHM13', 'n'],
    ['edit_h1', 'edit h1', 'n'], ['id_h1', 'ident h1', 'n', pct], ['edit_h2', 'edit h2', 'n'], ['id_h2', 'ident h2', 'n', pct],
    ['reads', 'reads following', 'n'], ['unique', 'reads only this', 'n'], ['closest', 'closest to truth'], ['same', 'same sequence as'], ['names', 'path names']], {max: 500});
}

// ------------------------------------------------------------------ graphs: metrics, Bandage, node tracks
const MET = [
  ['grp', 'graph (viewer)'],
  ['nodes', 'nodes', 'viewer', m => m.nodes, 'low'],
  ['nodes_per_kb', 'nodes per CHM13 kb', 'viewer', m => m.nodes_per_kb, 'low'],
  ['mean_node_bp', 'mean node bp', 'viewer', m => m.mean_node_bp, 'high'],
  ['node_bp', 'bp in nodes', 'viewer', m => m.node_bp, 'low'],
  ['private_bp', 'bp on single-path nodes', 'viewer', m => m.private_bp, 'low'],
  ['walks', 'distinct end-to-end walks', 'viewer', m => m.spanning_distinct_walks, null],
  ['cyc', 'walks revisiting a node / directed cycle', 'viewer', m => fmt(m.walks_revisiting_a_node) + ' / ' + (m.directed_cycle ? 'yes' : 'no'), null],
  ['cost', 'graph-implied cost / optimum, path pairs', 'viewer', m => m.cost_over_opt, 'low', 'For each pair of end-to-end walks, the nodes they share (same order and orientation) fix an alignment; its unit edit cost over the optimal pairwise edit distance, summed over pairs weighted by paths. Findings: hotspots 1.56, matched controls 1.09.'],
  ['exc', 'excess edits per kb', 'viewer', m => m.M_per_kb, 'low'],
  ['unal', 'unaligned homology per kb', 'viewer', m => m.U_per_kb, 'low', 'Bases between shared nodes that the optimum would match but the graph leaves unaligned.'],
  ['rcost', 'cost / optimum vs CHM13', 'viewer', m => m.ref_cost_over_opt, 'low'],
  ['svseg', 'SV-sized pieces per walk vs CHM13', 'viewer', m => m.ref_sv_segments_per_walk, 'low', 'Segments of >= 50 bp between nodes shared with CHM13, per walk.'],
  ['kmer', 'k-mer redundancy', 'viewer', m => m.kmer ? m.kmer.frac_extra : null, 'low', 'Share of 21-mer graph positions beyond the most copies any one walk has (upper bound on unmerged sequence).'],
  ['split', 'graph/MSA split vs the MSA on show', 'js', gi => splitOf(gi).frac, 'low', 'Over the columns of the MSA on show: pairs of rows with the same base on different nodes of this graph, as a share of all same-base pairs.'],
  ['grp', 'tools/evaluate.py (results/&lt;method&gt;/&lt;region&gt;.json)'],
  ['e_npk', 'nodes per kb', 'eval', 'size.nodes_per_kb', 'low'],
  ['e_1bp', 'share of 1-bp nodes', 'eval', 'size.node_frac_1bp', 'low'],
  ['e_bub', 'max alleles in a bubble', 'eval', 'size.max_alleles_any_bubble', 'low'],
  ['e_cost', 'cost / optimum (unit, all pairs)', 'eval', 'alignment.all_cost_over_opt', 'low'],
  ['e_rcost', 'cost / optimum vs CHM13', 'eval', 'alignment.ref_cost_over_opt', 'low'],
  ['e_aff', 'cost / optimum (affine)', 'eval', 'alignment.affine_all_cost_over_opt', 'low'],
  ['e_exc', 'excess edits per kb', 'eval', 'alignment.all_excess_per_kb', 'low'],
  ['e_sv', 'SV pieces per path (median)', 'eval', 'inflation.sv_pieces_median', 'low'],
  ['e_kmer', 'k-mer redundancy', 'eval', 'redundancy.kmer_frac_extra', 'low'],
  ['e_dg', 'truth edit to best graph path (h1+h2)', 'eval', 'truth.d_graph_sum', 'low'],
  ['e_dp', 'truth edit to closest panel path (h1+h2)', 'eval', 'truth.d_panel_sum', 'low'],
  ['e_raw', 'truth-in-graph truvari F1, raw', 'eval', 'truth.raw_f1', 'high'],
  ['e_ref', 'truth-in-graph truvari F1, refined', 'eval', 'truth.refined_f1', 'high'],
];
function metricValue(def, gi) { const g = D.graphs[gi];
  if (def[2] === 'viewer') return def[3](g.metrics);
  if (def[2] === 'js') return def[3](gi);
  if (def[2] === 'eval') return g.evaluate && g.evaluate.m ? g.evaluate.m[def[3]] : undefined;
  return undefined; }
function metricsTable() {
  const anyEval = D.graphs.some(g => g.evaluate);
  let h = '<div style="overflow-x:auto"><table class="mt"><thead><tr><th>metric</th>' + D.graphs.map(g => '<th class="n" title="' + esc(g.src) + '">' + esc(g.name) + '</th>').join('') + '</tr></thead><tbody>';
  let skipGroup = false;
  for (const def of MET) {
    if (def[0] === 'grp') { skipGroup = def[1].indexOf('evaluate') >= 0 && !anyEval; if (!skipGroup) h += '<tr><td class="grp" colspan="' + (D.graphs.length + 1) + '">' + def[1] + '</td></tr>'; continue; }
    if (skipGroup) continue;
    const vals = D.graphs.map((g, gi) => metricValue(def, gi));
    const nums = vals.filter(v => typeof v === 'number');
    const best = def[4] && nums.length > 1 ? (def[4] === 'low' ? Math.min(...nums) : Math.max(...nums)) : null;
    h += '<tr><td' + (def[5] ? ' title="' + esc(def[5]) + '"' : '') + '>' + esc(def[1]) + (def[5] ? ' <span class="muted">ⓘ</span>' : '') + '</td>' + vals.map(v => '<td class="n' + (best !== null && v === best ? ' best' : '') + '">' + (def[0] === 'split' || def[0] === 'kmer' || def[0] === 'e_kmer' || def[0] === 'e_1bp' ? pct(v) : fmt(v)) + '</td>').join('') + '</tr>';
  }
  if (anyEval) h += '<tr><td class="grp" colspan="' + (D.graphs.length + 1) + '">' + D.graphs.map(g => esc(g.name) + ': ' + (g.evaluate ? esc(g.evaluate.file) + (g.evaluate.graph_matches ? '' : ' <b>(made for another GFA: ' + esc(g.evaluate.graph) + ')</b>') + (g.evaluate.valid === false ? ' <b>INVALID</b>' : '') : 'no evaluate.py result')).join('<br>') + '</td></tr>';
  h += '</tbody></table></div><div class="legend">Bold = best of the graphs shown. Viewer metrics are computed by the viewer from the end-to-end paths (up to ' + '--max-pairs pairs); tools/evaluate.py is the defined yardstick (tools/METRICS.md) and is shown when its results file exists.</div>';
  $('metrics-table').innerHTML = h;
}
function graphLine(gi) { const m = D.graphs[gi].metrics;
  return fmt(m.nodes) + ' nodes (' + fmt(m.nodes_per_kb) + '/kb, mean ' + fmt(m.mean_node_bp) + ' bp) · cost/opt ' + fmt(m.cost_over_opt) + ' · vs CHM13 ' + fmt(m.ref_cost_over_opt) + ' · SV pieces/walk ' + fmt(m.ref_sv_segments_per_walk) + ' · k-mer extra ' + pct(m.kmer ? m.kmer.frac_extra : null); }
function graphsSection() {
  $('graph-sub').textContent = D.graphs.length + ' graph' + (D.graphs.length > 1 ? 's' : '') + ': ' + D.graphs.map(g => g.name).join(' · ');
  metricsTable();
  let w = ''; D.graphs.forEach(g => { for (const x of (g.warnings || [])) w += '<div class="warn">' + esc(g.name) + ': ' + esc(x) + '</div>'; }); $('graph-warn').innerHTML = w;
  let h = '';
  D.graphs.forEach((g, gi) => { h += '<figure id="bf-' + gi + '"><div class="btitle"><b>' + esc(g.name) + '</b> <span class="muted">' + esc(g.kind) + ' · ' + graphLine(gi) + '</span></div>' +
    (g.bandage ? '<img alt="Bandage rendering of ' + esc(g.name) + '" src="' + g.bandage + '">' : '<div class="muted">No Bandage image' + (g.bandage_note ? ': ' + esc(g.bandage_note) : '') + '.</div>') +
    '<figcaption><code>' + esc(g.src) + '</code></figcaption></figure>'; });
  $('bandage-grid').innerHTML = h;
  D.graphs.forEach((g, gi) => { const f = $('bf-' + gi); const img = f.querySelector('img'); if (img) img.addEventListener('click', () => f.classList.toggle('full')); });
  $('bandage-legend').innerHTML = 'Bandage renderings (force-directed layout, not to scale with the MSA; click to toggle full size). Node colour: ' + sw(D.ref_grey) + 'on the CHM13 path; others by paths through the node: ' + shareLegend() + '. Node width grows with that count.';
  $('node-controls').innerHTML = '<label>x axis <select id="nmode"><option value="msa">MSA columns (synced with the alignment)</option><option value="bp">path bp (unaligned)</option><option value="steps">node steps (1 step = 1 block)</option></select></label>' +
    '<label>colour <select id="ncolor"><option value="node">node identity</option><option value="share">node sharing</option><option value="split">graph/MSA split</option></select></label>' +
    '<label id="nzoomlab" style="display:none">zoom <select id="nzoom"><option value="fit">fit</option><option value="1">1 px</option><option value="3">3 px</option><option value="8">8 px</option><option value="16">16 px</option></select></label>' +
    '<label><input type="checkbox" id="nfrag"> include fragment walks (non-MSA modes)</label><span class="muted">one track per graph; rows follow the alignment\'s rows; click a block to outline every visit to that node</span>';
  $('nmode').addEventListener('change', e => { NODE.mode = e.target.value; $('nzoomlab').style.display = NODE.mode === 'msa' ? 'none' : ''; buildAllNodeLabels(); drawAllNodes(); });
  $('ncolor').addEventListener('change', e => { NODE.color = e.target.value; legendNodes(); drawAllNodes(); });
  $('nzoom').addEventListener('change', e => { NODE.zoom = e.target.value; drawAllNodes(); });
  $('nfrag').addEventListener('change', e => { NODE.frag = e.target.checked; buildAllNodeLabels(); drawAllNodes(); });
  let b = '';
  D.graphs.forEach((g, gi) => { b += '<div class="ntblock"><h4>' + esc(g.name) + ' <span class="sub">' + graphLine(gi) + '</span></h4><div class="track"><div class="labels" id="nl-' + gi + '"></div><div class="scroller" id="ns-' + gi + '"><div class="spacer"></div><canvas id="nc-' + gi + '"></canvas></div></div><div class="status" id="nst-' + gi + '"></div></div>'; });
  $('node-blocks').innerHTML = b;
  D.graphs.forEach((g, gi) => {
    registerScroller($('ns-' + gi), () => drawNodes(gi), () => NODE.mode === 'msa');
    const cv = $('nc-' + gi);
    cv.addEventListener('mousemove', e => { const h = nodeAt(gi, e); if (h) showTip(e, nodeTip(gi, h)); else hideTip(); }); cv.addEventListener('mouseleave', hideTip);
    cv.addEventListener('click', e => { const h = nodeAt(gi, e); const G = GR[gi]; G.sel = h ? (G.sel === h.n ? -1 : h.n) : -1;
      $('nst-' + gi).innerHTML = G.sel >= 0 ? 'outlined: node ' + esc(G.N.id[G.sel]) + ' (' + G.N.len[G.sel] + ' bp), used by ' + G.N.nw[G.sel] + ' walks / ' + G.N.w[G.sel] + ' paths' : ''; drawNodes(gi); });
  });
  legendNodes();
}
function legendNodes() { let s;
  if (NODE.color === 'node') s = sw(COL.ref1) + sw(COL.ref2) + 'nodes on the CHM13 path (alternating) &nbsp;' + sw('hsl(30,68%,50%)') + sw('hsl(200,68%,50%)') + 'other nodes, hashed colour: the same node has the same colour in every row of its graph (colours are not comparable between graphs). A thick outline marks a node the walk visits more than once.';
  else if (NODE.color === 'share') s = shareLegend() + ' paths through the node';
  else s = sw(COL.match) + 'no split &nbsp;' + sw('rgba(224,49,49,1)') + 'the MSA on show aligns this base to the same base in other rows that sit on a different node (MSA-column mode only)';
  $('node-legend').innerHTML = s; }
function nodeRows(gi) { const G = GR[gi], rw = G.g.row_walk;
  if (NODE.mode === 'msa') return visibleRows().filter(i => colNode(gi, i)).map(i => ({walk: rw[i], row: i}));
  const used = new Set(), out = [], hidden = new Set();
  for (const i of VIEW.hidden) if (rw[i] >= 0) hidden.add(rw[i]);
  for (const i of visibleRows()) { const wi = rw[i]; if (wi >= 0 && !used.has(wi)) { used.add(wi); out.push({walk: wi, row: i}); } }
  G.walks.forEach((w, wi) => { if (!used.has(wi) && !hidden.has(wi) && (w.cls === 'spanning' || NODE.frag)) { used.add(wi); out.push({walk: wi, row: -1}); } });
  return out; }
function buildAllNodeLabels() { for (let gi = 0; gi < GR.length; gi++) buildNodeLabels(gi); }
function buildNodeLabels(gi) { const el = $('nl-' + gi); if (!el) return; const G = GR[gi];
  let s = '<div class="hd" style="height:16px">' + (NODE.mode === 'msa' ? 'MSA columns' : NODE.mode === 'bp' ? 'path bp' : 'node steps') + '</div>';
  for (const o of nodeRows(gi)) { const w = G.walks[o.walk]; const lab = o.row >= 0 ? ROWS[o.row].label : w.id;
    s += '<div class="lab" style="height:' + NODE.rowH + 'px;line-height:' + NODE.rowH + 'px;font-size:11px"><span class="k-' + (w.is_ref ? 'ref' : 'panel') + '">' + esc(lab) + '</span> <span class="badge">×' + w.weight + '</span> <span class="muted">' + (w.cls === 'spanning' ? '' : esc(w.cls) + ' ') + fmt(w.steps.length) + ' steps</span></div>'; }
  el.innerHTML = s; }
function nodeScale(gi, rows) { if (NODE.mode === 'msa') return VIEW.px; const el = $('ns-' + gi); const G = GR[gi];
  const maxLen = Math.max(1, ...rows.map(o => NODE.mode === 'bp' ? G.walks[o.walk].off[G.walks[o.walk].steps.length] : G.walks[o.walk].steps.length));
  return NODE.zoom === 'fit' ? Math.max(0.01, (el.clientWidth - 2) / maxLen) : +NODE.zoom; }
function drawAllNodes() { for (let gi = 0; gi < GR.length; gi++) drawNodes(gi); }
function blockColor(G, n, row, c) { if (NODE.color === 'share') return shareColor(G.N.w[n]); if (NODE.color === 'split') return row >= 0 && NODE.mode === 'msa' ? splitColor(splitOf(G.gi).cells[row][c]) : COL.match; return nodeColor(G, n); }
function drawNodes(gi) {
  const G = GR[gi]; const el = $('ns-' + gi), cv = $('nc-' + gi); if (!el) return; const W = el.clientWidth; const rows = nodeRows(gi); const rh = NODE.rowH; const H = 16 + rows.length * rh + 2;
  const px = nodeScale(gi, rows); const spacer = el.querySelector('.spacer');
  if (NODE.mode === 'msa') { spacer.style.width = Math.ceil(NC * VIEW.px) + 'px'; if (Math.abs(el.scrollLeft - $('aln-scroll').scrollLeft) > 0.5) el.scrollLeft = $('aln-scroll').scrollLeft; }
  else { const maxLen = Math.max(1, ...rows.map(o => NODE.mode === 'bp' ? G.walks[o.walk].off[G.walks[o.walk].steps.length] : G.walks[o.walk].steps.length)); spacer.style.width = Math.ceil(maxLen * px) + 'px'; }
  const ctx = prepCanvas(cv, W, H); const x0 = el.scrollLeft / px;
  ctx.fillStyle = COL.muted; ctx.font = '10px ui-monospace,Menlo,monospace'; ctx.textBaseline = 'top';
  ctx.fillText(NODE.mode === 'msa' ? 'column ' + Math.floor(x0 + 1) + ' · CHM13 ' + (col2ref[Math.min(NC - 1, Math.floor(x0))] || '').toLocaleString('en-US') : (NODE.mode === 'bp' ? 'bp ' : 'step ') + Math.floor(x0 + 1), 2, 2);
  rows.forEach((o, k) => { const y = 16 + k * rh; const w = G.walks[o.walk];
    if (k % 2) { ctx.fillStyle = COL.grid; ctx.globalAlpha = 0.35; ctx.fillRect(0, y, W, rh); ctx.globalAlpha = 1; }
    let vis = G.visits.get(o.walk); if (!vis) { vis = new Map(); for (const s of w.steps) { const n = Math.abs(s) - 1; vis.set(n, (vis.get(n) || 0) + 1); } G.visits.set(o.walk, vis); }
    if (NODE.mode === 'msa') { const cn = colNode(gi, o.row); if (!cn) return; const c0 = Math.max(0, Math.floor(x0)), c1 = Math.min(NC, Math.ceil(x0 + W / px) + 1);
      let run = -2, rs = 0; const flush = (cE) => { if (run >= 0) { const xa = (rs - x0) * px, xb = (cE - x0) * px; ctx.fillStyle = blockColor(G, run, o.row, rs); ctx.fillRect(xa, y + 1, Math.max(0.6, xb - xa - (px >= 3 ? 1 : 0)), rh - 2);
          if (run === G.sel || vis.get(run) > 1) { ctx.strokeStyle = run === G.sel ? COL.hl : COL.fg; ctx.lineWidth = run === G.sel ? 2 : 1.5; ctx.strokeRect(xa + 0.5, y + 1.5, Math.max(1, xb - xa - 1), rh - 3); } } };
      for (let c = c0; c < c1; c++) { const n = cn[c]; if (n !== run) { flush(c); run = n; rs = c; } } flush(c1);
    } else { const bp = NODE.mode === 'bp'; let i0 = 0; if (bp) { let lo = 0, hi = w.steps.length; while (lo < hi) { const m = (lo + hi) >> 1; if (w.off[m + 1] <= x0) lo = m + 1; else hi = m; } i0 = lo; } else i0 = Math.max(0, Math.floor(x0));
      for (let i = i0; i < w.steps.length; i++) { const a = bp ? w.off[i] : i, b = bp ? w.off[i + 1] : i + 1; const xa = (a - x0) * px, xb = (b - x0) * px; if (xa > W) break; const n = Math.abs(w.steps[i]) - 1;
        ctx.fillStyle = blockColor(G, n, -1, 0); ctx.fillRect(xa, y + 1, Math.max(0.6, xb - xa - (xb - xa >= 3 ? 1 : 0)), rh - 2);
        if (w.steps[i] < 0 && xb - xa >= 6) { ctx.fillStyle = COL.bg; ctx.fillRect(xa + 1, y + rh / 2 - 0.5, xb - xa - 2, 1); }
        if (n === G.sel || vis.get(n) > 1) { ctx.strokeStyle = n === G.sel ? COL.hl : COL.fg; ctx.lineWidth = n === G.sel ? 2 : 1.5; ctx.strokeRect(xa + 0.5, y + 1.5, Math.max(1, xb - xa - 1), rh - 3); } } }
  });
}
function nodeAt(gi, e) { const G = GR[gi]; const el = $('ns-' + gi), rect = $('nc-' + gi).getBoundingClientRect(); const x = e.clientX - rect.left, y = e.clientY - rect.top; const rows = nodeRows(gi); const k = Math.floor((y - 16) / NODE.rowH);
  if (k < 0 || k >= rows.length) return null; const o = rows[k]; const w = G.walks[o.walk]; const px = nodeScale(gi, rows); const u = el.scrollLeft / px + x / px;
  if (NODE.mode === 'msa') { const c = Math.floor(u); if (c < 0 || c >= NC) return null; const n = colNode(gi, o.row)[c]; if (n < 0) return null; let step = -1; for (let i = 0; i < w.steps.length; i++) if (Math.abs(w.steps[i]) - 1 === n) { step = i; break; } return {n, o, step, c}; }
  let i; if (NODE.mode === 'bp') { i = -1; for (let s = 0; s < w.steps.length; s++) if (w.off[s + 1] > u) { i = s; break; } } else i = Math.floor(u);
  if (i < 0 || i >= w.steps.length) return null; return {n: Math.abs(w.steps[i]) - 1, o, step: i, rev: w.steps[i] < 0}; }
function nodeTip(gi, h) { const G = GR[gi]; const w = G.walks[h.o.walk]; return '<b>' + esc(G.g.name) + ' node ' + esc(G.N.id[h.n]) + '</b> (' + G.N.len[h.n] + ' bp, ' + (G.N.ref[h.n] ? 'on the CHM13 path' : 'off CHM13') + ')<br>' + esc(h.o.row >= 0 ? ROWS[h.o.row].label : w.id) + ' step ' + (h.step + 1) + (h.rev ? ' (reverse)' : '') + '<br>used by ' + G.N.nw[h.n] + ' distinct walks, ' + G.N.w[h.n] + ' paths' + (h.c !== undefined ? '<br>MSA column ' + (h.c + 1) + ', CHM13 ' + col2ref[h.c].toLocaleString('en-US') + '<br>same base on another node: ' + pct(splitOf(gi).cells[h.o.row][h.c]) : ''); }

// ------------------------------------------------------------------ calls panel
function callsPanel() {
  const chv = D.called_vs_truth || {}; const c = D.counts;
  if (!D.calls.length && !D.truth.length) { $('calls-sub').textContent = 'none'; $('calls-body').innerHTML = '<div class="muted">No calls: the package has no calls.tsv and region.py was not run.</div>'; return; }
  $('calls-sub').textContent = 'span; vg: ' + c.vg_tp + ' TP, ' + c.vg_fp + ' FP alleles, ' + c.vg_fn + ' FN truth SVs · PanGenie: ' + c.pg_tp + ' TP, ' + c.pg_fp + ' FP, ' + c.pg_fn + ' FN';
  let h = '';
  if (chv.vg || chv.pangenie) { h = '<div style="overflow-x:auto"><table style="width:auto"><thead><tr><th>called haplotypes vs truth (whole span)</th><th class="n">edit to h1</th><th class="n">edit to h2</th><th class="n">total</th><th>pairing</th><th>note</th></tr></thead><tbody>';
    for (const [k, lab] of [['reference_as_call', 'CHM13 (call nothing)'], ['vg', 'vg'], ['pangenie', 'PanGenie']]) { const v = chv[k]; if (!v) continue;
      const worse = k !== 'reference_as_call' && chv.reference_as_call && v.total > chv.reference_as_call.total;
      h += '<tr><td>' + lab + '</td><td class="n">' + fmt(v.d_h1) + '</td><td class="n">' + fmt(v.d_h2) + '</td><td class="n"><b' + (worse ? ' class="s-FP"' : '') + '>' + fmt(v.total) + '</b></td><td>' + esc(v.pairing) + '</td><td>' + (worse ? 'worse than calling nothing' : '') + (v.phase_reliable === false ? ' phase unreliable' : '') + '</td></tr>'; }
    h += '</tbody></table></div>'; }
  $('called-dist').innerHTML = h;
  const callCols = [['pos', 'pos', 'n'], ['end', 'end', 'n'], ['ref_len', 'ref', 'n'], ['alt_len', 'alt', 'n'], ['len_diff', 'Δlen', 'n'], ['allele', 'allele', 'n'], ['gt', 'GT'], ['qual', 'QUAL', 'n'], ['filter', 'FILTER'],
    ['truvari', 'truvari'], ['pctseq', 'PctSeqSim', 'n'], ['pctsize', 'PctSizeSim', 'n'], ['sizediff', 'SizeDiff', 'n'], ['startdist', 'StartDist', 'n'], ['aardvark', 'aardvark'], ['id', 'id', '', v => '<span class="mono" title="' + esc(v) + '">' + esc(String(v).length > 34 ? String(v).slice(0, 34) + '…' : v) + '</span>']];
  for (const [caller, tid, cid] of [['vg', 'vg-table', 'vg-controls'], ['pangenie', 'pg-table', 'pg-controls']]) {
    const all = D.calls.filter(r => r.caller === caller);
    $(cid).innerHTML = '<label><input type="checkbox" class="f-sv" checked> |Δlen| ≥ 50 only</label><label><input type="checkbox" class="f-called" checked> called alleles only</label><label><input type="checkbox" class="f-err"> truvari FP only</label><span class="muted">' + all.length + ' alleles in span</span>';
    const tb = renderTable($(tid), [], callCols, {sort: 'pos'});
    const upd = () => { const sv = $(cid).querySelector('.f-sv').checked, cl = $(cid).querySelector('.f-called').checked, er = $(cid).querySelector('.f-err').checked;
      tb.update(all.filter(r => (!sv || Math.abs(r.len_diff || 0) >= 50) && (!cl || r.called) && (!er || r.truvari === 'FP'))); };
    $(cid).querySelectorAll('input').forEach(i => i.addEventListener('change', upd)); upd(); }
  const tcols = [['source', 'file'], ['pos', 'pos', 'n'], ['end', 'end', 'n'], ['ref_len', 'ref', 'n'], ['alt_lens', 'alt lens'], ['dlen', 'max Δlen', 'n'], ['svtype', 'type'], ['gt', 'GT'], ['h1', 'applied h1'], ['h2', 'applied h2'],
    ['vg', 'vg'], ['pg', 'PanGenie'], ['bench', 'in bench', 'n'], ['trf_period', 'TRF period', 'n'], ['rm', 'RM']];
  $('truth-controls').innerHTML = '<label><input type="checkbox" id="t-sv" checked> SV ≥ 50 bp (stvar) only</label><label><input type="checkbox" id="t-err"> vg or PanGenie FN only</label><span class="muted">' + D.truth.length + ' truth records in span</span>';
  const tt = renderTable($('truth-table'), [], tcols, {sort: 'pos'});
  const tu = () => { const sv = $('t-sv').checked, er = $('t-err').checked; tt.update(D.truth.filter(t => (!sv || (t.source === 'stvar' && t.sv50)) && (!er || t.vg === 'FN' || t.pg === 'FN'))); };
  $('t-sv').addEventListener('change', tu); $('t-err').addEventListener('change', tu); tu();
}

// ------------------------------------------------------------------ reads panel
// read tuple: 0 name,1 mapq,2 strand,3 qlen,4 ident,5 nsteps,6 first,7 last,8 inside,9 status,10 compat,11 walk,12 b0,13 b1,14 partial,15 pd
const RD = D.reads || [];
let PLC = [];
function baseCols(r) { if (b2cCache.has(r)) return b2cCache.get(r); const s = rowStr(r); let o = null;
  if (s) { const cols = []; for (let c = 0; c < s.length; c++) if (s[c] !== '-') cols.push(c); o = {b0: M.b0[r], cols}; } b2cCache.set(r, o); return o; }
function placeReads() { PLC = RD.map(rd => { if (rd[11] < 0) return null; const row = rowOfBaseWalk.get(rd[11]); if (row === undefined) return null; const bc = baseCols(row); if (!bc || !bc.cols.length) return null;
  const i0 = rd[12] - bc.b0, i1 = rd[13] - 1 - bc.b0; if (i1 < 0 || i0 >= bc.cols.length) return null;
  const c0 = bc.cols[Math.max(0, i0)], c1 = bc.cols[Math.min(bc.cols.length - 1, i1)]; return c1 >= c0 ? [c0, c1] : null; }); }
function mapqCol(q) { return q === 0 ? COL.mq0 : q < 5 ? COL.mq1 : q < 30 ? COL.mq5 : q < 60 ? COL.mq30 : COL.mq60; }
function readsPanel() {
  if (!D.reads) { $('reads-sub').textContent = 'none'; $('reads-body').innerHTML = '<div class="muted">No reads: the big data (GAF-Base) is not available here, or --no-reads was given. The region package holds sequences only.</div>'; return; }
  const rs = D.reads_summary || {}; const st = D.read_status || {};
  $('reads-sub').textContent = fmt(rs.n_alignments) + ' alignments (' + fmt(rs.n_distinct_read_names) + ' names) overlapping the subgraph, from GAF-Base; placed on the baseline graph';
  const hist = rs.mapq_histogram || {}; const ks = ['0', '1-4', '5-9', '10-19', '20-29', '30-59', '60', '61-255'].filter(k => k in hist); const mx = Math.max(1, ...ks.map(k => hist[k]));
  const cols = {'0': COL.mq0, '1-4': COL.mq1, '5-9': COL.mq5, '10-19': COL.mq5, '20-29': COL.mq5, '30-59': COL.mq30, '60': COL.mq60, '61-255': COL.mq60};
  const bw = 44, H = 130; let svg = '<svg viewBox="0 0 ' + (ks.length * bw + 10) + ' ' + (H + 34) + '" width="100%" style="max-width:' + (ks.length * bw + 10) + 'px" role="img" aria-label="MAPQ histogram">';
  ks.forEach((k, i) => { const v = hist[k], h = v / mx * H; svg += '<rect x="' + (i * bw + 6) + '" y="' + (H - h + 12) + '" width="' + (bw - 10) + '" height="' + h + '" fill="' + cols[k] + '" rx="2"></rect><text x="' + (i * bw + 6 + (bw - 10) / 2) + '" y="' + (H - h + 9) + '" font-size="10" text-anchor="middle">' + v + '</text><text x="' + (i * bw + 6 + (bw - 10) / 2) + '" y="' + (H + 26) + '" font-size="10" text-anchor="middle">' + k + '</text>'; });
  svg += '</svg>'; $('mapq-hist').innerHTML = '<div class="muted" style="font-size:12px">MAPQ histogram</div>' + svg;
  let f = '<div class="facts">';
  f += '<div>observed ' + fmt(rs.n_alignments) + ' vs expected ' + fmt(rs.expected_from_truth_haplotypes) + ' from the HG002 haplotype lengths (' + fmt(rs.observed_over_expected_truth) + '×) and ' + fmt(rs.expected_from_reference_span) + ' from CHM13 (' + fmt(rs.observed_over_expected_reference) + '×), at 30x</div>';
  f += '<div>MAPQ 0: <b>' + pct(rs.frac_mapq0) + '</b> · MAPQ&lt;5: <b>' + pct(rs.frac_mapq_lt5) + '</b> · median MAPQ ' + fmt(rs.median_mapq) + '</div>';
  f += '<div>block identity median ' + fmt(rs.median_block_identity) + ', &lt;0.9 in ' + pct(rs['frac_block_identity_lt_0.9']) + '</div>';
  f += '<div>node path followed by: an end-to-end walk <b>' + fmt(st.walk || 0) + '</b> · only a fragment walk ' + fmt(st.fragment || 0) + ' · <b>no walk</b> (novel combination of edges) ' + fmt(st.novel || 0) + (st.outside ? ' · no subgraph node ' + st.outside : '') + '</div>';
  if (D.n_reads_total > RD.length) f += '<div class="warn">only ' + RD.length + ' of ' + D.n_reads_total + ' alignments loaded</div>';
  f += '</div>'; $('reads-facts').innerHTML = f;
  const bw0 = D.graphs[0].walks;
  const prow = (D.path_counts || []).map(p => { const w = bw0[p.walk]; return Object.assign({}, p, {id: w.id, cls: w.cls, weight: w.weight, len: w.len, per_kb: w.len ? +(p.compat / w.len * 1000).toFixed(1) : 0, row: rowOfBaseWalk.has(p.walk) ? 'yes' : 'no'}); });
  renderTable($('path-table'), prow, [['id', 'walk'], ['cls', 'class'], ['weight', 'paths', 'n'], ['len', 'bp', 'n'], ['row', 'in MSA'], ['compat', 'reads following', 'n'], ['compat_mq5', 'of which MAPQ≥5', 'n'], ['unique', 'following only this walk', 'n'], ['placed', 'drawn on this walk', 'n'], ['per_kb', 'following per kb', 'n']], {sort: 'compat', dir: -1, max: 400});
  $('pile-controls').innerHTML = '<label>colour <select id="pcol"><option value="mapq">MAPQ</option><option value="ident">block identity</option><option value="strand">strand vs CHM13</option><option value="status">path status</option></select></label>' +
    '<label>show <select id="pfilt"><option value="all">all reads</option><option value="mq5">MAPQ ≥ 5</option><option value="mq0">MAPQ &lt; 5</option></select></label>' +
    '<label>layout <select id="pgroup"><option value="packed">packed</option><option value="walk">grouped by walk drawn on</option></select></label>' +
    '<label><input type="checkbox" id="phl" checked> dim reads not following the comparison row\'s walk</label><span class="muted" id="pnote"></span>';
  ['pcol', 'pfilt', 'pgroup', 'phl'].forEach(id => $(id).addEventListener('change', () => { layoutPile(); drawPile(); legendPile(); }));
  legendPile(); layoutPile();
  $('rt-controls').innerHTML = '<label>show <select id="rtf"><option value="all">all</option><option value="novel">path followed by no walk</option><option value="mq0">MAPQ &lt; 5</option><option value="partial">placed partially</option></select></label>';
  const rows = RD.map((r, i) => ({i, name: r[0], mapq: r[1], strand: r[2], qlen: r[3], ident: r[4], nsteps: r[5], range: r[6] + ' … ' + r[7], inside: r[8] ? 'yes' : 'no', status: r[9],
    ncompat: r[10].length, compat: r[10].slice(0, 6).map(w => bw0[w].id).join(',') + (r[10].length > 6 ? ',…' : ''), walk: r[11] >= 0 ? bw0[r[11]].id : '', partial: r[14] ? 'partial' : ''}));
  const tb = renderTable($('reads-table'), rows, [['name', 'read', '', v => '<span class="mono">' + esc(v) + '</span>'], ['mapq', 'MAPQ', 'n'], ['strand', 'strand'], ['qlen', 'len', 'n'], ['ident', 'block ident', 'n'], ['nsteps', 'path nodes', 'n'], ['range', 'path node range', '', v => '<span class="mono">' + esc(v) + '</span>'],
    ['inside', 'path inside subgraph'], ['status', 'path status'], ['ncompat', '# walks followed', 'n'], ['compat', 'walks followed'], ['walk', 'drawn on'], ['partial', 'placement']], {max: 3000});
  $('rtf').addEventListener('change', e => { const v = e.target.value; tb.update(rows.filter(r => v === 'all' || (v === 'novel' && r.status === 'novel') || (v === 'mq0' && r.mapq < 5) || (v === 'partial' && r.partial))); });
}
function legendPile() { if (!D.reads) return; const v = $('pcol').value; let s;
  if (v === 'mapq') s = sw(COL.mq0) + '0 ' + sw(COL.mq1) + '1-4 ' + sw(COL.mq5) + '5-29 ' + sw(COL.mq30) + '30-59 ' + sw(COL.mq60) + '60';
  else if (v === 'ident') s = sw(COL.fp) + '&lt;0.9 ' + sw(COL.fn) + '0.9-0.98 ' + sw(COL.tp) + '≥0.98 block identity';
  else if (v === 'strand') s = sw(COL.mq30) + 'forward vs CHM13 ' + sw(COL.called) + 'reverse';
  else s = sw(COL.tp) + 'follows an end-to-end walk ' + sw(COL.fn) + 'only a fragment walk (drawn partially) ' + sw(COL.fp) + 'no walk follows its path (drawn on the best partial match)';
  $('pile-legend').innerHTML = s + ' &nbsp; Hatched = partial placement. Reads are placed on the row of the highest-weight end-to-end baseline walk whose node path contains theirs, then projected onto the MSA on show.'; }
let PILE = {items: [], groups: [], H: 0};
const PLH = 4;
function readColor(r) { const v = $('pcol').value; if (v === 'mapq') return mapqCol(r[1]); if (v === 'ident') return r[4] < 0.9 ? COL.fp : r[4] < 0.98 ? COL.fn : COL.tp;
  if (v === 'strand') return r[2] === '+' ? COL.mq30 : COL.called; return r[9] === 'walk' ? COL.tp : r[9] === 'fragment' ? COL.fn : COL.fp; }
function layoutPile() {
  if (!D.reads) return;
  const f = $('pfilt').value; const grp = $('pgroup').value;
  const items = RD.map((r, i) => i).filter(i => PLC[i] && (f === 'all' || (f === 'mq5' ? RD[i][1] >= 5 : RD[i][1] < 5)));
  items.sort((a, b) => PLC[a][0] - PLC[b][0] || PLC[b][1] - PLC[a][1]);
  const groups = []; const byG = new Map();
  for (const i of items) { const g = grp === 'walk' ? RD[i][11] : -1; if (!byG.has(g)) { byG.set(g, []); groups.push(g); } byG.get(g).push(i); }
  if (grp === 'walk') groups.sort((a, b) => (rowOfBaseWalk.get(a) || 0) - (rowOfBaseWalk.get(b) || 0));
  const MAXL = 400; let y = 16; const lay = []; let overflow = 0; const out = [];
  for (const g of groups) { const ends = []; const top = y + (grp === 'walk' ? 12 : 0); let nl = 0;
    for (const i of byG.get(g)) { const c0 = PLC[i][0]; let l = -1; for (let k = 0; k < ends.length; k++) if (ends[k] < c0 - 1) { l = k; break; }
      if (l < 0) { if (ends.length >= MAXL) { overflow++; continue; } l = ends.length; ends.push(0); } ends[l] = PLC[i][1]; out.push({i, y: top + l * PLH}); nl = Math.max(nl, l + 1); }
    lay.push({g, y, n: byG.get(g).length}); y = top + nl * PLH + 6; }
  PILE = {items: out, groups: lay, H: y + 4, grp};
  $('pnote').textContent = out.length + ' reads drawn' + (overflow ? ', ' + overflow + ' not drawn (lane limit)' : '') + '; ' + (RD.length - items.length) + ' not placeable or filtered';
  const el = $('pile-labels'); let s = '<div class="hd" style="height:16px">pileup</div>';
  const bw0 = D.graphs[0].walks;
  if (grp === 'walk') { let prevY = 16; for (const L of lay) { s += '<div style="height:' + (L.y - prevY) + 'px"></div><div class="hd" style="height:12px">' + esc(bw0[L.g] ? bw0[L.g].id : '?') + ' (' + L.n + ')</div>'; prevY = L.y + 12; } }
  el.innerHTML = s;
}
function drawPile() {
  if (!D.reads) return;
  const el = $('pile-scroll'), cv = $('pile-canvas'); const W = el.clientWidth; const ctx = prepCanvas(cv, W, Math.max(40, PILE.H)); const px = VIEW.px, col0 = viewCol0(el);
  ctx.fillStyle = COL.muted; ctx.font = '10px ui-monospace,Menlo,monospace'; ctx.textBaseline = 'top'; ctx.fillText('column ' + Math.floor(col0 + 1) + ' · CHM13 ' + (col2ref[Math.min(NC - 1, Math.floor(col0))] || '').toLocaleString('en-US'), 2, 2);
  const selWalk = ROWS[VIEW.cmp] && ROWS[VIEW.cmp].walk !== undefined ? ROWS[VIEW.cmp].walk : -1; const dim = $('phl').checked && selWalk >= 0 && VIEW.cmp !== 0;
  if (PILE.grp === 'walk') for (const L of PILE.groups) { ctx.fillStyle = COL.grid; ctx.fillRect(0, L.y + 11, W, 1); }
  for (const it of PILE.items) { const r = RD[it.i], p = PLC[it.i]; const x0 = (p[0] - col0) * px, x1 = (p[1] + 1 - col0) * px; if (x1 < 0 || x0 > W) continue;
    ctx.globalAlpha = dim && !r[10].includes(selWalk) ? 0.15 : 1; ctx.fillStyle = readColor(r); ctx.fillRect(x0, it.y, Math.max(1, x1 - x0), PLH - 1);
    if (r[14]) { ctx.fillStyle = COL.bg; for (let x = Math.ceil(x0 / 6) * 6; x < x1; x += 6) ctx.fillRect(x, it.y, 2, PLH - 1); } }
  ctx.globalAlpha = 1;
}
function pileHover(e) { const el = $('pile-scroll'), rect = $('pile-canvas').getBoundingClientRect(); const x = e.clientX - rect.left, y = e.clientY - rect.top; const c = viewCol0(el) + x / VIEW.px;
  const hit = PILE.items.find(it => y >= it.y && y < it.y + PLH && c >= PLC[it.i][0] && c <= PLC[it.i][1] + 1); if (!hit) { hideTip(); return; } const r = RD[hit.i], p = PLC[hit.i]; const bw0 = D.graphs[0].walks;
  showTip(e, '<b class="mono">' + esc(r[0]) + '</b><br>MAPQ ' + r[1] + ' · strand ' + r[2] + ' · ' + r[3] + ' bp · block identity ' + r[4] + '<br>path ' + r[5] + ' nodes ' + esc(r[6]) + ' … ' + esc(r[7]) + (r[8] ? '' : ' (leaves the subgraph)') +
    '<br>status ' + esc(r[9]) + (r[14] ? ' (partial placement)' : '') + '<br>follows ' + r[10].length + ' walks: ' + esc(r[10].slice(0, 12).map(w => bw0[w].id).join(', ')) + (r[10].length > 12 ? ', …' : '') + '<br>drawn on ' + esc(bw0[r[11]] ? bw0[r[11]].id : '') + ', MSA columns ' + (p[0] + 1) + '-' + (p[1] + 1) + ', CHM13 ≈' + col2ref[p[0]].toLocaleString('en-US') + '-' + col2ref[p[1]].toLocaleString('en-US')); }

// ------------------------------------------------------------------ init
function init() {
  loadColors(); setMsa(0); header();
  buildAnnotations(); alnControls(); legendAln(); sortRows('default'); msaSub();
  registerScroller($('aln-scroll'), drawAlignment);
  if (D.reads) registerScroller($('pile-scroll'), drawPile);
  graphsSection(); callsPanel(); placeReads(); readsPanel(); walkTable();
  buildLabels(); buildAllNodeLabels(); updateHiddenNote();
  setSpacers();
  const c = colOfRef(D.req[0]); setZoom(12, Math.max(0, c) + $('aln-scroll').clientWidth / 24);
  $('aln-canvas').addEventListener('mousemove', alnHover); $('aln-canvas').addEventListener('mouseleave', hideTip);
  if (D.reads) { $('pile-canvas').addEventListener('mousemove', pileHover); $('pile-canvas').addEventListener('mouseleave', hideTip); }
  const ov = $('overview'); let ovDown = false; ov.addEventListener('mousedown', e => { ovDown = true; overviewNav(e); }); window.addEventListener('mouseup', () => ovDown = false); ov.addEventListener('mousemove', e => { if (ovDown) overviewNav(e); });
  drawOverview(); drawAllNodes();
  let rT = null; window.addEventListener('resize', () => { clearTimeout(rT); rT = setTimeout(() => { drawAllMsa(); drawOverview(); drawAllNodes(); }, 120); });
  const recolor = () => { loadColors(); legendAln(); legendNodes(); legendPile(); drawAllMsa(); drawOverview(); drawAllNodes(); };
  if (window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', recolor);
  const tb = $('themebtn'); const themes = ['auto', 'light', 'dark']; let ti = 0;
  tb.addEventListener('click', () => { ti = (ti + 1) % 3; if (themes[ti] === 'auto') document.documentElement.removeAttribute('data-theme'); else document.documentElement.setAttribute('data-theme', themes[ti]); tb.textContent = 'theme: ' + themes[ti]; recolor(); });
}
init();
'''


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('target', help='region id (L012184), region package directory, chr:start-end, or a contig '
                                   'followed by start and end')
    ap.add_argument('start', type=int, nargs='?', help='1-based inclusive (with a contig)')
    ap.add_argument('end', type=int, nargs='?', help='1-based inclusive (with a contig)')
    ap.add_argument('--candidate', action='append', default=[], metavar='GFA',
                    help='candidate graph of the region (repeatable); FILE.msa.fa next to it is used as its MSA')
    ap.add_argument('--candidate-label', action='append', default=[], metavar='NAME',
                    help='name of the corresponding --candidate (default: its directory name)')
    ap.add_argument('--out', help='output HTML (default: <work>/viewer/<name>.html)')
    ap.add_argument('--workroot', default=os.path.join(config.WORK_DIR, 'viewer'),
                    help='viewer intermediates and region.py output (default: %(default)s)')
    ap.add_argument('--dir', help='existing region.py output directory to use (never written)')
    ap.add_argument('--no-data', action='store_true',
                    help='do not use the big data even if present (package-only page, as a collaborator sees it)')
    ap.add_argument('--pad', type=int, default=0, help='region.py --pad when region.py is run on coordinates')
    ap.add_argument('--names', action='store_true', help='region.py --names when region.py is run (slow)')
    ap.add_argument('--no-reads', action='store_true', help='no reads (region.py --no-reads; no reads panel)')
    ap.add_argument('--rerun', action='store_true', help='re-run region.py even if its output exists')
    ap.add_argument('--mafft', choices=['auto', 'linsi', 'einsi', 'ginsi', 'fftns'], default='auto',
                    help='mafft strategy of the all-rows MSA (default --auto)')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--max-cols', type=int, default=40000, help='show at most this many MSA columns')
    ap.add_argument('--max-seq', type=int, default=40000,
                    help='cut sequences longer than this to the CHM13 window around the core before the mafft MSA')
    ap.add_argument('--max-reads', type=int, default=20000)
    ap.add_argument('--max-pairs', type=int, default=600,
                    help='distinct-walk pairs for the viewer\'s graph metrics (default 600, i.e. all of them for hap32)')
    ap.add_argument('--no-bandage', action='store_true')
    ap.add_argument('--no-candidate-msa', action='store_true', help='ignore <candidate>.msa.fa files')
    ap.add_argument('--no-add', action='store_true',
                    help='do not add truth and called rows to a candidate MSA with mafft --add')
    ap.add_argument('--max-page-mb', type=float, default=5.0, help='page size to aim for (default 5)')
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args(argv)
    if (args.start is None) != (args.end is None):
        ap.error('give both start and end')
    if args.start is not None and args.end < args.start:
        ap.error('end < start')
    if args.candidate_label and len(args.candidate_label) > len(args.candidate):
        ap.error('more --candidate-label than --candidate')
    try:
        build(args)
    except (ViewError, vr.ToolError) as e:
        sys.stderr.write('viewer: ERROR: %s\n' % e)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
