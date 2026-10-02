#!/usr/bin/env python3
"""gap_norm.py -- Stage 4m: gap placement in the full-panel abPOA MSA.

Two things:

  normalise(rows, direction)  shift every row's gap runs left (or right) against the column consensus,
                              as VCF left-alignment against a reference, iterated to a fixed point.
                              The consensus of a column is its most common base, gaps not counted, so
                              insertion columns have a consensus too (that is what lets an insertion
                              move: the other rows' gap runs slide past it). A gap run [g, e) moves one
                              column left when the row's base at g-1 equals the consensus of column e-1;
                              that base then sits at e-1. Shifts only go one way, so the loop ends. The
                              consensus is recomputed between passes. All-gap columns are dropped. Row
                              sequences never change (asserted).

  diagnose(rid)               where the projected full-panel graph's extra repeated k-mers come from,
                              against abPOA on hap32 alone. Positions of a 21-mer in a column-induced
                              graph are the start columns of its occurrences (evaluate.kmer_redundancy
                              counts graph positions; in a column graph a position is a column). An
                              occurrence pair that the hap32-only MSA puts in one column but the
                              projected full MSA puts in two is SPLIT. For each split pair the two rows'
                              induced pairwise alignments are compared between the nearest base pairs
                              both MSAs align (so both align the same two subsequences there):
                                equal    same edit cost: the full panel chose an equally good alignment
                                         that places the indel elsewhere (placement inconsistency)
                                f_worse  the full-panel alignment of these two rows costs more there
                                f_better it costs less (hap32-only is the worse one)
                              Merging the split positions joined by 'equal' pairs gives the excess that
                              placement alone explains.

    python3 tools/gap_norm.py diagnose TR773368,TR770009 [--json out.json]
    python3 tools/gap_norm.py kmerx some.msa.fa
"""
import argparse
import bisect
import collections
import json
import os
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import msa_graph  # noqa: E402

K = 21


# ------------------------------------------------------------------ normalisation

def consensus(rows):
    """Per column, the most common non-gap base (ties: alphabetical), or None for an all-gap column."""
    out = []
    for col in zip(*rows):
        c = collections.Counter(col)
        c.pop('-', None)
        out.append(min(c, key=lambda b: (-c[b], b)) if c else None)
    return out


def _shift_left(row, cons):
    """Shift row's (a list) gap runs left in place while the moved base matches the consensus of the
    column it moves into. Returns the number of single-column moves."""
    n, moved, c = len(row), 0, 0
    while c < n:
        if row[c] != '-':
            c += 1
            continue
        g = e = c
        while e < n and row[e] == '-':
            e += 1
        old_e = e
        while g > 0:
            b = row[g - 1]
            if b == '-':                      # ran into the previous gap run: they are one run now
                while g > 0 and row[g - 1] == '-':
                    g -= 1
                continue
            if b != cons[e - 1]:
                break
            row[e - 1], row[g - 1] = b, '-'
            g -= 1
            e -= 1
            moved += 1
        c = old_e
    return moved


def normalise(rows, direction='left', cons_rows=None, max_passes=100000):
    """(normalised rows, stats). rows: aligned strings. cons_rows: indices of the rows the consensus is
    taken over (default all). Columns that end up all-gap are dropped."""
    if direction not in ('left', 'right'):
        raise ValueError(direction)
    rev = direction == 'right'
    work = [list(r[::-1] if rev else r) for r in rows]
    passes = moves = 0
    while passes < max_passes:
        cons = consensus([work[i] for i in cons_rows] if cons_rows is not None else work)
        m = sum(_shift_left(r, cons) for r in work)
        passes += 1
        moves += m
        if not m:
            break
    else:
        raise RuntimeError('no fixed point after %d passes' % max_passes)
    out = [''.join(r[::-1] if rev else r) for r in work]
    keep = [c for c in range(len(out[0])) if any(r[c] != '-' for r in out)] if out else []
    if out and len(keep) < len(out[0]):
        out = [''.join(r[c] for c in keep) for r in out]
    for a, b in zip(rows, out):
        if a.replace('-', '') != b.replace('-', ''):
            raise AssertionError('normalisation changed a row sequence')
    return out, {'passes': passes, 'moves': moves, 'columns_in': len(rows[0]) if rows else 0,
                 'columns_out': len(out[0]) if out else 0}


def normalise_joint(rows, direction='left', max_passes=100000):
    """As normalise, but gap runs that several rows share exactly (same start and end column) move
    together: a shared run moves one column only if the move is allowed in every row that has it, and
    stops where any of those rows would merge it with its previous run. Per-row normalisation scatters
    a shared run to as many places as its rows have variant bases next to it; this keeps it whole."""
    if direction not in ('left', 'right'):
        raise ValueError(direction)
    rev = direction == 'right'
    work = [list(r[::-1] if rev else r) for r in rows]
    passes = moves = 0
    while passes < max_passes:
        cons = consensus(work)
        groups = collections.defaultdict(list)
        for i, r in enumerate(work):
            for run in gap_runs(r):
                groups[run].append(i)
        m = 0
        for (g, e), members in sorted(groups.items()):
            while g > 0:
                want = cons[e - 1]
                if want is None or any(work[i][g - 1] != want for i in members):
                    break                       # a gap at g-1 (merge) never equals a base, so stops too
                for i in members:
                    work[i][e - 1], work[i][g - 1] = work[i][g - 1], '-'
                g -= 1
                e -= 1
                m += 1
        passes += 1
        moves += m
        if not m:
            break
    else:
        raise RuntimeError('no fixed point after %d passes' % max_passes)
    out = [''.join(r[::-1] if rev else r) for r in work]
    for a, b in zip(rows, out):
        if a.replace('-', '') != b.replace('-', ''):
            raise AssertionError('normalisation changed a row sequence')
    out = drop_gap_cols(out)
    return out, {'passes': passes, 'moves': moves, 'columns_in': len(rows[0]) if rows else 0,
                 'columns_out': len(out[0]) if out else 0}


def gap_runs(row):
    """[(start, end)] of the maximal gap runs of an aligned row."""
    out, c, n = [], 0, len(row)
    while c < n:
        if row[c] == '-':
            d = c
            while d < n and row[d] == '-':
                d += 1
            out.append((c, d))
            c = d
        else:
            c += 1
    return out


def merge_runs(rows, lam=1.0, max_passes=1000):
    """(rows, stats): consolidate each row's fragmented indels. For two consecutive gap runs G1 = [a, b),
    G2 = [c, d) of a row and the bases S = [b, c) between them, S may move left (to [a, a+|S|)) or right
    (to [d-|S|, d)), which joins G1 and G2 into one run. The change in the row's linear sum-of-pairs cost
    against the other rows is
        delta = sum over the columns a..d-1 of (after - before),
    where a base x at column j costs N - n_j(x) and a gap costs N - n_j('-') (n_j = the other rows'
    counts at column j, N = the number of other rows). The cheaper move is made when delta <= lam * N,
    i.e. when it costs at most lam edits per other row (the price of the gap open it saves). One left-to-
    right sweep per row per pass, repeated until no row changes; every merge removes a gap run."""
    work = [list(r) for r in rows]
    R = len(work)
    if R < 2:
        return [''.join(r) for r in work], {'passes': 0, 'merges': 0}
    N = R - 1
    counts = [collections.Counter(col) for col in zip(*work)]
    passes = merges = 0
    while passes < max_passes:
        passes += 1
        m_pass = 0
        for row in work:
            def own(j, ch):                      # other rows' count of ch at column j
                return counts[j][ch] - (1 if row[j] == ch else 0)

            def cost_seg(chars, lo):             # chars placed at columns lo.. (gaps included)
                s = 0
                for t, ch in enumerate(chars):
                    s += N - own(lo + t, ch)
                return s
            runs = gap_runs(row)
            i = 0
            while i + 1 < len(runs):
                (a, b), (c, d) = runs[i], runs[i + 1]
                seg = row[b:c]
                m = len(seg)
                before = cost_seg(row[a:d], a)
                left = seg + ['-'] * (d - a - m)
                right = ['-'] * (d - a - m) + seg
                cl, cr = cost_seg(left, a), cost_seg(right, a)
                best, new = (cl, left) if cl <= cr else (cr, right)
                if best - before <= lam * N:
                    for t, ch in enumerate(new):
                        j = a + t
                        if row[j] != ch:
                            counts[j][row[j]] -= 1
                            counts[j][ch] += 1
                            row[j] = ch
                    m_pass += 1
                    runs[i:i + 2] = [(a + m, d) if new is left else (a, d - m)]
                    i = max(i - 1, 0)
                    continue
                i += 1
        merges += m_pass
        if not m_pass:
            break
    out = [''.join(r) for r in work]
    for a, b in zip(rows, out):
        if a.replace('-', '') != b.replace('-', ''):
            raise AssertionError('merge changed a row sequence')
    return drop_gap_cols(out), {'passes': passes, 'merges': merges}


# ------------------------------------------------------------------ k-mer redundancy of a column graph

def kmer_excess(rows, k=K):
    """evaluate.kmer_redundancy on the column-induced graph of rows: positions are start columns."""
    pos = collections.defaultdict(set)
    maxocc = collections.Counter()
    for row in rows:
        cols = [c for c, ch in enumerate(row) if ch != '-']
        seq = ''.join(row[c] for c in cols)
        occ = collections.Counter()
        for x in range(len(seq) - k + 1):
            km = seq[x:x + k]
            if 'N' in km:
                continue
            pos[km].add(cols[x])
            occ[km] += 1
        for km, c in occ.items():
            if c > maxocc[km]:
                maxocc[km] = c
    P = sum(len(v) for v in pos.values())
    extra = sum(len(v) - maxocc[km] for km, v in pos.items())
    return {'positions': P, 'extra': extra, 'kmer_frac_extra': round(extra / P, 5) if P else 0.0}


# ------------------------------------------------------------------ diagnosis

class _Rows:
    def __init__(self, rows):
        self.rows = rows
        self.col = [[c for c, ch in enumerate(r) if ch != '-'] for r in rows]
        self.idx = []
        for r, cols in zip(rows, self.col):
            a = [-1] * len(r)
            for i, c in enumerate(cols):
                a[c] = i
            self.idx.append(a)

    def partner(self, a, b):
        """For each base i of row a, the base of row b in the same column (or -1)."""
        ib = self.idx[b]
        return [ib[c] for c in self.col[a]]

    def cost(self, a, b, i0, i1):
        """Edit cost of the induced a/b alignment strictly between a's bases i0 and i1 (-1 and len are the
        row ends), where both bounding bases are aligned to each other (or are the ends)."""
        ra, rb = self.rows[a], self.rows[b]
        c0 = self.col[a][i0] if i0 >= 0 else -1
        c1 = self.col[a][i1] if i1 < len(self.col[a]) else len(ra)
        n = 0
        for c in range(c0 + 1, c1):
            x, y = ra[c], rb[c]
            if x == '-' and y == '-':
                continue
            if x != y:
                n += 1
        return n


def diagnose_rows(names, f_rows, h_rows, k=K, examples=0):
    """f_rows: projected full-panel MSA; h_rows: hap32-only MSA; same names and sequences, same order."""
    import region   # region.edit_distance (C)
    seqs = [r.replace('-', '') for r in f_rows]
    if seqs != [r.replace('-', '') for r in h_rows]:
        raise ValueError('row sequences differ')
    F, H = _Rows(f_rows), _Rows(h_rows)
    # distinct rows only: identical sequences share an MSA row in both
    first = {}
    uniq = [i for i, s in enumerate(seqs) if first.setdefault(s, i) == i]
    inst = collections.defaultdict(list)            # kmer -> [(row, offset)]
    for i in uniq:
        s = seqs[i]
        for x in range(len(s) - k + 1):
            km = s[x:x + k]
            if 'N' not in km:
                inst[km].append((i, x))
    pair_cache = {}

    def pair(a, b):
        key = (a, b)
        if key not in pair_cache:
            pf, ph = F.partner(a, b), H.partner(a, b)
            comm = [i for i, (u, v) in enumerate(zip(pf, ph)) if u == v and u >= 0]
            pair_cache[key] = (pf, ph, comm)
        return pair_cache[key]

    def classify(a, x1, b, x2):
        pf, ph, comm = pair(a, b)
        La = len(seqs[a])
        j = bisect.bisect_left(comm, x1)
        i0 = comm[j - 1] if j > 0 else -1
        j = bisect.bisect_left(comm, x1 + k)
        i1 = comm[j] if j < len(comm) else La
        j0 = ph[i0] if i0 >= 0 else -1
        j1 = ph[i1] if i1 < La else len(seqs[b])
        cf, ch = F.cost(a, b, i0, i1), H.cost(a, b, i0, i1)
        opt = region.edit_distance(seqs[a][i0 + 1:i1], seqs[b][j0 + 1:j1])
        cls = 'equal' if cf == ch else ('f_worse' if cf > ch else 'f_better')
        return cls, cf, ch, opt, i1 - i0 - 1

    totals = collections.Counter()
    link_n = collections.Counter()
    link_cost = collections.defaultdict(lambda: [0, 0, 0])      # class -> [cost_f, cost_h, opt] summed
    seg_len = collections.defaultdict(list)
    shown = []
    for km, occ in inst.items():
        fpos = {F.col[i][x] for i, x in occ}
        hpos = {H.col[i][x] for i, x in occ}
        totals['P_f'] += len(fpos)
        totals['P_h'] += len(hpos)
        if len(fpos) <= 1:
            continue
        parent = {p: p for p in fpos}                          # union-find on f positions, 'equal' links
        parent_all = {p: p for p in fpos}                      # any link (all G_h co-locations)

        def find(par, p):
            while par[p] != p:
                par[p] = par[par[p]]
                p = par[p]
            return p
        byh = collections.defaultdict(dict)
        for i, x in occ:
            byh[H.col[i][x]].setdefault(F.col[i][x], (i, x))
        for hc, fp in byh.items():
            if len(fp) < 2:
                continue
            reps = sorted(fp.items())
            p0, (a, x1) = reps[0]
            for p1, (b, x2) in reps[1:]:
                cls, cf, ch, opt, L = classify(a, x1, b, x2)
                link_n[cls] += 1
                lc = link_cost[cls]
                lc[0] += cf
                lc[1] += ch
                lc[2] += opt
                seg_len[cls].append(L)
                if cls == 'equal':
                    parent[find(parent, p1)] = find(parent, p0)
                parent_all[find(parent_all, p1)] = find(parent_all, p0)
                if examples and len(shown) < examples and cls == 'equal' and ch == opt:
                    shown.append({'kmer': km, 'row_a': a, 'row_b': b, 'f_cols': [p0, p1], 'h_col': hc,
                                  'cost_f': cf, 'cost_h': ch, 'opt': opt})
        totals['merged_equal'] += len(fpos) - len({find(parent, p) for p in fpos})
        totals['merged_any'] += len(fpos) - len({find(parent_all, p) for p in fpos})
    ef, eh = kmer_excess(f_rows, k), kmer_excess(h_rows, k)
    gap = ef['extra'] - eh['extra']
    med = lambda v: sorted(v)[len(v) // 2] if v else None   # noqa: E731
    return {
        'rows': len(names), 'distinct_rows': len(uniq), 'kmerx_f': ef['kmer_frac_extra'],
        'kmerx_h': eh['kmer_frac_extra'], 'extra_f': ef['extra'], 'extra_h': eh['extra'], 'excess': gap,
        'split_links': dict(link_n), 'link_cost_f_h_opt': {c: v for c, v in link_cost.items()},
        'segment_len_median': {c: med(v) for c, v in seg_len.items()},
        'merged_by_equal': totals['merged_equal'], 'merged_by_any': totals['merged_any'],
        'explained_by_placement': round(totals['merged_equal'] / gap, 3) if gap > 0 else None,
        'explained_by_any_link': round(totals['merged_any'] / gap, 3) if gap > 0 else None,
        'examples': shown,
    }


def _ordered_like(names, rows_named):
    d = dict(rows_named)
    return [d[n] for n in names]


def drop_gap_cols(rows):
    keep = [c for c in range(len(rows[0])) if any(r[c] != '-' for r in rows)]
    return [''.join(r[c] for c in keep) for r in rows]


def diagnose(rid, f_msa, h_msa, full_msa=None, map_tsv=None, examples=3):
    """Diagnosis of one region; with the full MSA and map, also the kmerx of the projection after
    normalising the full MSA (left/right) and after normalising the projection itself."""
    f = msa_graph.read_msa(f_msa)
    names = [n for n, _ in f]
    f_rows = [r for _, r in f]
    h_rows = drop_gap_cols(_ordered_like(names, msa_graph.read_msa(h_msa)))
    out = {'region_id': rid}
    out.update(diagnose_rows(names, f_rows, h_rows, examples=examples))
    for d in ('left', 'right'):
        pr, st = normalise(f_rows, d)
        out['kmerx_proj_norm_' + d] = kmer_excess(pr)['kmer_frac_extra']
    if full_msa and map_tsv:
        import panel
        full = [(n, r) for n, r in msa_graph.read_msa(full_msa) if n not in msa_graph.CONSENSUS_NAMES]
        uid_of = {}
        for r in panel.read_map(map_tsv):
            uid_of[r['id']] = r['id']
            for m in r['members']:
                uid_of[m] = r['id']
        fidx = {}
        for i, (n, _) in enumerate(full):
            fidx.setdefault(uid_of[n], i)
        sel = [fidx[uid_of[n]] for n in names]
        for d in ('left', 'right'):
            nr, st = normalise([r for _, r in full], d)
            pr = drop_gap_cols([nr[i] for i in sel])
            out['kmerx_full_norm_' + d] = kmer_excess(pr)['kmer_frac_extra']
            out['full_norm_' + d + '_stats'] = st
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    d = sub.add_parser('diagnose')
    d.add_argument('regions')
    d.add_argument('--f-dir', default=os.path.join(os.path.dirname(TOOLS), 'work', 'stage4', 'candidates',
                                                   'poa_abpoa__all'))
    d.add_argument('--h-dir', default=os.path.join(os.path.dirname(TOOLS), 'work', 'stage4', 'candidates',
                                                   'poa_abpoa'))
    d.add_argument('--full-dir', default=os.path.join(os.path.dirname(TOOLS), 'work', 'stage4', 'panel', 'poa_abpoa'))
    d.add_argument('--examples', type=int, default=3)
    d.add_argument('--json')
    km = sub.add_parser('kmerx')
    km.add_argument('msa')
    a = ap.parse_args(argv)
    if a.cmd == 'kmerx':
        print(json.dumps(kmer_excess([r for _, r in msa_graph.read_msa(a.msa)])))
        return
    import poa_panel
    res = []
    for rid in a.regions.split(','):
        _, mp = poa_panel.union_files(rid)
        r = diagnose(rid, os.path.join(a.f_dir, rid + '.msa.fa'), os.path.join(a.h_dir, rid + '.msa.fa'),
                     os.path.join(a.full_dir, rid + '.msa.fa.gz'), mp, a.examples)
        res.append(r)
        print(json.dumps(r))
        sys.stdout.flush()
    if a.json:
        with open(a.json, 'w') as f:
            json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
