#!/usr/bin/env python3
"""msa_graph.py -- the shared MSA -> graph builder for every realigner.

    python3 tools/msa_graph.py MSA.fa regions/<id>/hap32.fa OUT.gfa
                               [--engine auto|vg|native] [--merge-blocks] [--block-max 50]
                               [--allow-missing] [--ignore-extra] [--json STATS.json] [--quiet]

    from msa_graph import msa_to_gfa
    stats = msa_to_gfa('L012184.msa.fa', 'regions/L012184/hap32.fa', 'L012184.gfa')

Input: aligned rows in FASTA ('-' or '.' gaps, any case, sequences may wrap), or PIR
(">P1;name", a description line, the sequence ending in '*'); abPOA -r 1/-r 2 and spoa -r 1
output are plain FASTA. A row is matched to hap32.fa by the first word of its header. Every
row's ungapped sequence must equal its hap32.fa sequence (case-insensitive): a mismatch is an
error, as are duplicate rows and hap32.fa sequences with no row (unless --allow-missing).
Rows whose name is not in hap32.fa are an error unless --ignore-extra; abPOA's
"Consensus_sequence" row is always dropped (with a note).

Output: GFA 1.0 with S lines (upper case), L lines (all forward) and one P line per row,
named as in hap32.fa and in hap32.fa order, each spelling exactly its hap32.fa sequence
(re-checked from the written file). Node ids are 1..n in column order, so every edge goes
from a lower to a higher id and the graph is acyclic.

Graph induction. Column induction: in each MSA column, the rows with the same base share a
1-bp node; each row's path joins its consecutive bases. The graph is then compacted: a node
whose only successor has it as only predecessor is merged with it, unless a path ends or
starts between them (this is `vg mod -u`, unchop). --engine vg builds the column graph with
`vg construct -M` (rows renamed first, because vg reads '#' in a name as PanSN and rewrites
the name) followed by `vg mod -u`; --engine native does the same in Python; auto uses vg when
it is installed and checks that the native graph is identical (node sequences of every path).
Both then pass through the same compaction and renumbering.

--merge-blocks (native only) removes the 1-bp mesh that column induction leaves where many
SNVs sit close together: the MSA is cut into blocks, maximal runs of columns in which the set
of rows with a base (not a gap) does not change; a block of at most --block-max columns gets
one node per distinct row string instead of one per distinct base per column. The default
50 keeps every merged bubble below SV size. Node counts are reported both ways, whichever is
written.

N: an N is a base like any other, never merged with A/C/G/T, so N runs become nodes of their
own and every path still spells its original sequence; the stats report how many N bp (and
nodes) the graph has. All-gap columns are dropped.
"""
import argparse
import collections
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
try:
    import config  # noqa: E402
    VG = config.VG
    TOOL_ENV = config.tool_env()
except ImportError:  # usable as a single file
    VG = shutil.which('vg') or 'vg'
    TOOL_ENV = dict(os.environ)

GAP_CHARS = '-.'
CONSENSUS_NAMES = ('Consensus_sequence',)


class MsaGraphError(RuntimeError):
    pass


# ---------------------------------------------------------------- input

def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_fasta(path):
    """[(first word of header, sequence)] in file order."""
    out, name, buf = [], None, []
    with _open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ''.join(buf)))
                h = line[1:].split()
                name, buf = (h[0] if h else ''), []
            elif line:
                buf.append(line)
    if name is not None:
        out.append((name, ''.join(buf)))
    return out


def read_msa(path):
    """Aligned rows [(name, row)] from FASTA or PIR. Rows are upper-cased, '.' becomes '-'."""
    with _open(path) as f:
        text = f.read()
    lines = text.splitlines()
    pir = any(ln.startswith('>') and len(ln) > 4 and ln[3] == ';' for ln in lines[:5])
    rows = []
    if pir:
        i = 0
        while i < len(lines):
            ln = lines[i].strip()
            if not ln.startswith('>'):
                i += 1
                continue
            name = ln[4:].split()[0] if ln[4:].split() else ''
            i += 2                        # skip the description line
            buf = []
            while i < len(lines) and not lines[i].startswith('>'):
                buf.append(lines[i].strip())
                i += 1
            s = ''.join(buf)
            if s.endswith('*'):
                s = s[:-1]
            rows.append((name, s))
    else:
        name, buf = None, []
        for ln in lines:
            ln = ln.strip()
            if ln.startswith('>'):
                if name is not None:
                    rows.append((name, ''.join(buf)))
                h = ln[1:].split()
                name, buf = (h[0] if h else ''), []
            elif ln:
                buf.append(ln)
        if name is not None:
            rows.append((name, ''.join(buf)))
    out = []
    for n, s in rows:
        s = s.upper().replace(' ', '')
        for g in GAP_CHARS[1:]:
            s = s.replace(g, '-')
        out.append((n, s))
    return out


def write_msa(rows, path, width=0):
    """Write [(name, aligned row)] as FASTA ('-' gaps); width 0 = one line per row."""
    with open(path, 'w') as f:
        for n, s in rows:
            f.write('>%s\n' % n)
            if width:
                for i in range(0, len(s), width):
                    f.write(s[i:i + width] + '\n')
            else:
                f.write(s + '\n')


def match_rows(msa, hap, allow_missing=False, ignore_extra=False):
    """Pair MSA rows with hap32.fa records. Returns (names, rows, notes): names in hap32.fa
    order, rows the matching aligned strings."""
    notes = []
    want = collections.OrderedDict()
    for n, s in hap:
        if n in want:
            raise MsaGraphError('duplicate name in hap32.fa: %s' % n)
        want[n] = s.upper()
    got = collections.OrderedDict()
    extra = []
    for n, s in msa:
        if n in got:
            raise MsaGraphError('duplicate MSA row: %s' % n)
        if n not in want:
            if n in CONSENSUS_NAMES:
                notes.append('dropped MSA row %s' % n)
                continue
            extra.append(n)
            continue
        got[n] = s
    if extra:
        if not ignore_extra:
            raise MsaGraphError('%d MSA rows are not in hap32.fa: %s' % (len(extra), ', '.join(extra[:5])))
        notes.append('ignored %d MSA rows not in hap32.fa (%s)' % (len(extra), ', '.join(extra[:5])))
    missing = [n for n in want if n not in got]
    if missing:
        if not allow_missing:
            raise MsaGraphError('%d hap32.fa sequences have no MSA row: %s' % (len(missing), ', '.join(missing[:5])))
        notes.append('%d hap32.fa sequences have no MSA row (%s)' % (len(missing), ', '.join(missing[:5])))
    names = [n for n in want if n in got]
    if not names:
        raise MsaGraphError('no MSA row matches hap32.fa')
    lens = set(len(got[n]) for n in names)
    if len(lens) != 1:
        raise MsaGraphError('MSA rows have different lengths: %s' % sorted(lens)[:6])
    bad = []
    for n in names:
        if got[n].replace('-', '') != want[n]:
            bad.append(n)
    if bad:
        n = bad[0]
        u = got[n].replace('-', '')
        i = next((k for k in range(min(len(u), len(want[n]))) if u[k] != want[n][k]), min(len(u), len(want[n])))
        raise MsaGraphError('%d MSA rows do not spell their hap32.fa sequence (%s: first difference at base %d, '
                            'row %d bp vs hap32 %d bp)' % (len(bad), n, i + 1, len(u), len(want[n])))
    rows = [got[n] for n in names]
    # drop all-gap columns
    ncol = len(rows[0])
    keep = [c for c in range(ncol) if any(r[c] != '-' for r in rows)]
    if len(keep) < ncol:
        notes.append('dropped %d all-gap columns' % (ncol - len(keep)))
        rows = [''.join(r[c] for c in keep) for r in rows]
    return names, rows, notes


# ---------------------------------------------------------------- graphs

class G:
    """seq: {id: str}; edges: set of (a, b) forward; paths: [list of ids] (row order)."""

    def __init__(self):
        self.seq = {}
        self.edges = set()
        self.paths = []

    def spell(self, i):
        return ''.join(self.seq[n] for n in self.paths[i])


def columns_graph(rows, lo=0, hi=None, g=None, prev=None):
    """Column induction over columns [lo, hi) of rows, extending g (per-row last node in prev)."""
    if g is None:
        g = G()
        g.paths = [[] for _ in rows]
    if prev is None:
        prev = [None] * len(rows)
    hi = len(rows[0]) if hi is None else hi
    nid = max(g.seq) + 1 if g.seq else 1
    seq, edges, paths = g.seq, g.edges, g.paths
    for c in range(lo, hi):
        col = {}
        for r, row in enumerate(rows):
            b = row[c]
            if b == '-':
                continue
            n = col.get(b)
            if n is None:
                n = col[b] = nid
                seq[nid] = b
                nid += 1
            p = prev[r]
            if p is not None:
                edges.add((p, n))
            prev[r] = n
            paths[r].append(n)
    return g, prev


def blocks_of(rows):
    """Maximal runs of columns with the same set of non-gap rows: [(lo, hi, rowset)]."""
    out = []
    ncol = len(rows[0])
    cur, lo = None, 0
    for c in range(ncol):
        s = frozenset(r for r, row in enumerate(rows) if row[c] != '-')
        if s != cur:
            if cur is not None:
                out.append((lo, c, cur))
            cur, lo = s, c
    if cur is not None:
        out.append((lo, ncol, cur))
    return out


def blocks_graph(rows, block_max):
    """Column induction, except that blocks of <= block_max columns get one node per distinct
    row string."""
    g = G()
    g.paths = [[] for _ in rows]
    prev = [None] * len(rows)
    merged = 0
    for lo, hi, rs in blocks_of(rows):
        if hi - lo > block_max or len(rs) < 2:
            columns_graph(rows, lo, hi, g, prev)
            continue
        merged += 1
        nid = max(g.seq) + 1 if g.seq else 1
        by = {}
        for r in sorted(rs):
            s = rows[r][lo:hi]
            n = by.get(s)
            if n is None:
                n = by[s] = nid
                g.seq[nid] = s
                nid += 1
            if prev[r] is not None:
                g.edges.add((prev[r], n))
            prev[r] = n
            g.paths[r].append(n)
    return g, merged


def unchop(g):
    """Merge u->v when v is u's only successor, u is v's only predecessor, and no path ends
    at u or starts at v. Node ids of the result are 1..n in the order of the input ids."""
    out_n = collections.defaultdict(set)
    in_n = collections.defaultdict(set)
    for a, b in g.edges:
        out_n[a].add(b)
        in_n[b].add(a)
    ends = set(p[-1] for p in g.paths if p)
    starts = set(p[0] for p in g.paths if p)

    def joins(u):
        """the node u merges into its successor with, or None"""
        if len(out_n[u]) != 1 or u in ends:
            return None
        v = next(iter(out_n[u]))
        if len(in_n[v]) != 1 or v in starts or v == u:
            return None
        return v
    nxt = {}
    has_prev = set()
    for u in g.seq:
        v = joins(u)
        if v is not None:
            nxt[u] = v
            has_prev.add(v)
    order = sorted(g.seq)
    new_id, tail = {}, {}
    ng = G()
    k = 0
    for u in order:
        if u in has_prev:
            continue
        k += 1
        parts = [g.seq[u]]
        new_id[u] = k
        x = u
        while x in nxt:
            x = nxt[x]
            parts.append(g.seq[x])
            new_id[x] = k
        tail[k] = x
        ng.seq[k] = ''.join(parts)
    for a, b in g.edges:
        if a in nxt and nxt[a] == b:
            continue
        ng.edges.add((new_id[a], new_id[b]))
    for p in g.paths:
        q = []
        for n in p:
            m = new_id[n]
            if not q or q[-1] != m or n not in has_prev:
                q.append(m)
        ng.paths.append(q)
    # renumber topologically (ids along edges must increase): Kahn by smallest id
    indeg = collections.Counter()
    succ = collections.defaultdict(list)
    for a, b in ng.edges:
        indeg[b] += 1
        succ[a].append(b)
    import heapq
    heap = [n for n in ng.seq if indeg[n] == 0]
    heapq.heapify(heap)
    topo = []
    while heap:
        n = heapq.heappop(heap)
        topo.append(n)
        for m in succ[n]:
            indeg[m] -= 1
            if indeg[m] == 0:
                heapq.heappush(heap, m)
    if len(topo) != len(ng.seq):
        raise MsaGraphError('the induced graph has a cycle (internal error)')
    ren = dict((n, i + 1) for i, n in enumerate(topo))
    fg = G()
    fg.seq = dict((ren[n], s) for n, s in ng.seq.items())
    fg.edges = set((ren[a], ren[b]) for a, b in ng.edges)
    fg.paths = [[ren[n] for n in p] for p in ng.paths]
    return fg


def vg_graph(rows, workdir, quiet=True):
    """Column graph from `vg construct -M` + `vg mod -u` on rows renamed s0..sN."""
    if not (os.path.exists(VG) or shutil.which(VG)):
        raise MsaGraphError('vg not found (%s)' % VG)
    tmp = tempfile.mkdtemp(prefix='msa_graph.', dir=workdir)
    try:
        fa = os.path.join(tmp, 'msa.fa')
        write_msa([('s%d' % i, r) for i, r in enumerate(rows)], fa)
        vgf = os.path.join(tmp, 'c.vg')
        ugf = os.path.join(tmp, 'u.vg')
        with open(vgf, 'wb') as o:
            p = subprocess.run([VG, 'construct', '-M', fa, '-m', '1000000'], stdout=o, stderr=subprocess.PIPE,
                               env=TOOL_ENV)
        if p.returncode:
            raise MsaGraphError('vg construct -M failed: %s' % p.stderr.decode(errors='replace')[:400])
        with open(ugf, 'wb') as o:
            p = subprocess.run([VG, 'mod', '-u', vgf], stdout=o, stderr=subprocess.PIPE, env=TOOL_ENV)
        if p.returncode:
            raise MsaGraphError('vg mod -u failed: %s' % p.stderr.decode(errors='replace')[:400])
        p = subprocess.run([VG, 'view', ugf], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=TOOL_ENV)
        if p.returncode:
            raise MsaGraphError('vg view failed: %s' % p.stderr.decode(errors='replace')[:400])
        text = p.stdout.decode()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    seq, edges, praw = {}, set(), {}
    for line in text.splitlines():
        f = line.split('\t')
        if f[0] == 'S':
            seq[int(f[1])] = f[2].upper()
        elif f[0] == 'L':
            if f[2] != '+' or f[4] != '+':
                raise MsaGraphError('vg made a reversing edge %s' % line)
            edges.add((int(f[1]), int(f[3])))
        elif f[0] == 'P':
            st = [x for x in f[2].split(',') if x]
            if any(x[-1] != '+' for x in st):
                raise MsaGraphError('vg path %s visits a node in reverse' % f[1])
            praw[f[1]] = [int(x[:-1]) for x in st]
        elif f[0] == 'W':
            raise MsaGraphError('vg wrote W lines despite the renamed rows')
    g = G()
    g.seq, g.edges = seq, edges
    g.paths = []
    for i in range(len(rows)):
        if 's%d' % i not in praw:
            raise MsaGraphError('vg lost path s%d' % i)
        g.paths.append(praw['s%d' % i])
    return g


def canonical(g):
    """Engine-independent form: per path, the list of its node sequences."""
    return [tuple(g.seq[n] for n in p) for p in g.paths]


# ---------------------------------------------------------------- output

def write_gfa(g, names, path):
    tmp = path + '.tmp%d' % os.getpid()
    with open(tmp, 'w') as f:
        f.write('H\tVN:Z:1.0\n')
        for n in sorted(g.seq):
            f.write('S\t%d\t%s\n' % (n, g.seq[n]))
        for a, b in sorted(g.edges):
            f.write('L\t%d\t+\t%d\t+\t0M\n' % (a, b))
        for nm, p in zip(names, g.paths):
            f.write('P\t%s\t%s\t*\n' % (nm, ','.join('%d+' % n for n in p)))
    os.replace(tmp, path)


def check_gfa(path, want):
    """Re-read a GFA and check every P line against want {name: seq}; returns problems."""
    seq, paths = {}, {}
    with open(path) as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            if x[0] == 'S':
                seq[x[1]] = x[2]
            elif x[0] == 'P':
                paths[x[1]] = [s for s in x[2].split(',') if s]
    probs = []
    for nm, s in want.items():
        if nm not in paths:
            probs.append('%s: no P line' % nm)
            continue
        sp = ''.join(seq[st[:-1]] for st in paths[nm])
        if sp != s:
            probs.append('%s: P line spells %d bp, expected %d' % (nm, len(sp), len(s)))
    return probs


def graph_stats(g, ref_len=None):
    lens = [len(s) for s in g.seq.values()]
    nbp = sum(s.count('N') for s in g.seq.values())
    return {'nodes': len(lens), 'edges': len(g.edges), 'node_bp': sum(lens),
            'nodes_per_kb': round(1000.0 * len(lens) / ref_len, 1) if ref_len else None,
            'mean_node_bp': round(sum(lens) / len(lens), 2) if lens else None,
            'nodes_1bp': sum(1 for x in lens if x == 1),
            'frac_nodes_1bp': round(sum(1 for x in lens if x == 1) / len(lens), 4) if lens else None,
            'max_node_bp': max(lens) if lens else 0,
            'N_bp_in_nodes': nbp, 'N_nodes': sum(1 for s in g.seq.values() if 'N' in s)}


# ---------------------------------------------------------------- API

def msa_to_gfa(msa_fasta, hap32_fasta, out_gfa, engine='auto', merge_blocks=False, block_max=50,
               allow_missing=False, ignore_extra=False, workdir=None, quiet=True):
    """Build out_gfa from msa_fasta (rows named as in hap32_fasta). Returns a stats dict;
    raises MsaGraphError when the MSA does not reproduce hap32.fa or the graph fails its check."""
    t0 = time.time()
    msa = read_msa(msa_fasta)
    hap = [(n, s.upper()) for n, s in read_fasta(hap32_fasta)]
    lower = sum(1 for n, s in read_fasta(hap32_fasta) if s != s.upper())
    names, rows, notes = match_rows(msa, hap, allow_missing, ignore_extra)
    if lower:
        notes.append('%d hap32.fa sequences have lower-case bases; the graph is upper case' % lower)
    want = dict((n, s) for n, s in hap if n in set(names))
    ref = next((n for n in names if n.split('#')[0] == 'CHM13'), None)
    ref_len = len(want[ref]) if ref else None
    stats = collections.OrderedDict()
    stats['msa'] = os.path.abspath(msa_fasta)
    stats['hap32'] = os.path.abspath(hap32_fasta)
    stats['gfa'] = os.path.abspath(out_gfa)
    stats['rows'] = len(names)
    stats['columns'] = len(rows[0])
    stats['ref_len'] = ref_len
    stats['gap_frac'] = round(sum(r.count('-') for r in rows) / float(len(rows) * len(rows[0])), 4)
    stats['N_bp_in_rows'] = sum(r.count('N') for r in rows)

    native = unchop(columns_graph(rows)[0])
    stats['columns_mode'] = graph_stats(native, ref_len)
    bg, nmerged = blocks_graph(rows, block_max)
    bg = unchop(bg)
    stats['blocks_mode'] = graph_stats(bg, ref_len)
    stats['blocks_mode']['block_max'] = block_max
    stats['blocks_mode']['blocks_merged'] = nmerged
    stats['blocks_mode']['blocks_total'] = len(blocks_of(rows))

    use_vg = engine == 'vg' or (engine == 'auto' and not merge_blocks and
                                (os.path.exists(VG) or shutil.which(VG)))
    if merge_blocks and engine == 'vg':
        raise MsaGraphError('--merge-blocks needs --engine native (or auto)')
    chosen = bg if merge_blocks else native
    stats['engine'] = 'native'
    if use_vg:
        try:
            vg_dir = workdir or os.path.dirname(os.path.abspath(out_gfa)) or '.'
            os.makedirs(vg_dir, exist_ok=True)   # the output directory may not exist yet
            raw = vg_graph(rows, vg_dir, quiet)
            stats['vg_raw_nodes'] = len(raw.seq)
            vgg = unchop(raw)
            stats['vg_unchop_changed'] = len(vgg.seq) != len(raw.seq)
            same = canonical(vgg) == canonical(native)
            stats['vg_matches_native'] = same
            if not same:
                notes.append('vg construct -M graph differs from the native column graph (%d vs %d nodes); '
                             'the vg graph is written' % (len(vgg.seq), len(native.seq)))
            chosen = vgg
            stats['engine'] = 'vg construct -M + vg mod -u'
        except MsaGraphError as e:
            if engine == 'vg':
                raise
            notes.append('vg engine failed (%s); native graph written' % str(e)[:200])
    stats['mode'] = 'blocks' if merge_blocks else 'columns'
    stats.update(graph_stats(chosen, ref_len))
    for i, n in enumerate(names):
        if chosen.spell(i) != want[n]:
            raise MsaGraphError('internal: path %s does not spell its sequence' % n)
    if any(a >= b for a, b in chosen.edges):
        raise MsaGraphError('internal: an edge goes backwards')
    os.makedirs(os.path.dirname(os.path.abspath(out_gfa)) or '.', exist_ok=True)
    write_gfa(chosen, names, out_gfa)
    probs = check_gfa(out_gfa, want)
    if probs:
        raise MsaGraphError('written GFA fails its path check: ' + '; '.join(probs[:5]))
    stats['paths_checked'] = len(want)
    stats['notes'] = notes
    stats['seconds'] = round(time.time() - t0, 2)
    return stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('msa', help='aligned FASTA / PIR, rows named as in hap32.fa')
    ap.add_argument('hap32', help='regions/<id>/hap32.fa')
    ap.add_argument('out', help='output GFA')
    ap.add_argument('--engine', choices=['auto', 'vg', 'native'], default='auto')
    ap.add_argument('--merge-blocks', action='store_true',
                    help='merge short blocks (constant set of non-gap rows) into whole-string alleles')
    ap.add_argument('--block-max', type=int, default=50, help='longest block to merge (default 50)')
    ap.add_argument('--allow-missing', action='store_true', help='allow hap32.fa sequences without a row')
    ap.add_argument('--ignore-extra', action='store_true', help='drop rows whose name is not in hap32.fa')
    ap.add_argument('--json', help='also write the stats here')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args(argv)
    try:
        st = msa_to_gfa(a.msa, a.hap32, a.out, engine=a.engine, merge_blocks=a.merge_blocks,
                        block_max=a.block_max, allow_missing=a.allow_missing, ignore_extra=a.ignore_extra,
                        quiet=a.quiet)
    except MsaGraphError as e:
        sys.stderr.write('msa_graph: ERROR: %s\n' % e)
        return 2
    if a.json:
        with open(a.json, 'w') as f:
            json.dump(st, f, indent=1)
    if not a.quiet:
        sys.stderr.write('msa_graph: %d rows x %d columns -> %d nodes (%s, %s mode; columns %d, blocks<=%d %d), '
                         '%d N bp in nodes, %.1f s\n' % (
                             st['rows'], st['columns'], st['nodes'], st['engine'], st['mode'],
                             st['columns_mode']['nodes'], a.block_max, st['blocks_mode']['nodes'],
                             st['N_bp_in_nodes'], st['seconds']))
    print(json.dumps(st))
    return 0


if __name__ == '__main__':
    sys.exit(main())
