#!/usr/bin/env python3
"""Compare pgrealign qc tables made with --graph: per-region F1 of the realigned alignment against the
original graph's ("before"), and between realign versions on the regions they share.

usage: qc_compare.py LABEL=qc.tsv [LABEL=qc.tsv ...]
Each table's graph_* columns are the graph extract read; tables measured against different scores have
different optima, so compare only tables made with the same --scores.
"""
import csv
import statistics
import sys

BINS = [(0, 1000, '<1 kb'), (1000, 10000, '1-10 kb'), (10000, 50000, '10-50 kb'), (50000, 10 ** 12, '>=50 kb')]
MARGIN = 0.01


def f1(tp, opt, msa):
    return 2 * tp / (opt + msa) if opt + msa else None


def load(path):
    out = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        r = {k: (v if k == 'id' else int(v)) for k, v in r.items()}
        for kind in ('pair', 'ref'):
            r[kind] = f1(r[kind + '_tp'], r[kind + '_opt'], r[kind + '_msa'])
            if 'graph_sequence' in r:
                r['g_' + kind] = f1(r['graph_%s_tp' % kind], r[kind + '_opt'], r['graph_%s_msa' % kind])
        out[r['id']] = r
    return out


def sums(rows, kind, graph=False):
    g = 'graph_' if graph else ''
    tp = sum(r[g + kind + '_tp'] for r in rows)
    opt = sum(r[kind + '_opt'] for r in rows)
    msa = sum(r[g + kind + '_msa'] for r in rows)
    return tp / msa if msa else float('nan'), tp / opt if opt else float('nan'), f1(tp, opt, msa) or float('nan')


def before_after(label, t):
    rows = list(t.values())
    print('\n## %s: %d regions, %d banded pairs' % (label, len(rows), sum(r['banded'] for r in rows)))
    print('\n| set | pairs | P after | R after | F1 after | P before | R before | F1 before | growth after | growth before |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    for lo, hi, name in [(0, 10 ** 12, 'all')] + BINS:
        sub = [r for r in rows if lo <= r['longest'] < hi]
        if not sub:
            continue
        L = sum(r['longest'] for r in sub)
        for kind in ('pair', 'ref'):
            a, b = sums(sub, kind), sums(sub, kind, True)
            print('| %s (%d) | %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.3f | %.3f |' % (
                name, len(sub), kind, *a, *b, sum(r['sequence'] for r in sub) / L,
                sum(r['graph_sequence'] for r in sub) / L))
    print('\nPer region (F1 of each region with pairs of that kind):')
    print('\n| kind | regions | mean after | mean before | after < 0.9 | before < 0.9 | after worse by >%.2f | after better by >%.2f |' % (MARGIN, MARGIN))
    print('|---|---|---|---|---|---|---|---|')
    for kind in ('pair', 'ref'):
        sub = [r for r in rows if r[kind] is not None and r['g_' + kind] is not None]
        if not sub:
            continue
        print('| %s | %d | %.4f | %.4f | %d | %d | %d | %d |' % (
            kind, len(sub), statistics.mean(r[kind] for r in sub), statistics.mean(r['g_' + kind] for r in sub),
            sum(r[kind] < 0.9 for r in sub), sum(r['g_' + kind] < 0.9 for r in sub),
            sum(r[kind] < r['g_' + kind] - MARGIN for r in sub), sum(r[kind] > r['g_' + kind] + MARGIN for r in sub)))
    worse = sorted((r for r in rows if r['ref'] is not None and r['g_ref'] is not None
                    and r['ref'] < r['g_ref'] - MARGIN), key=lambda r: r['ref'] - r['g_ref'])[:10]
    if worse:
        print('\nRegions whose reference F1 fell most:')
        print('\n| region | alleles | longest | ref F1 after | before | pair F1 after | before |')
        print('|---|---|---|---|---|---|---|')
        for r in worse:
            fmt = lambda x: '-' if x is None else '%.3f' % x
            print('| %s | %d | %d | %s | %s | %s | %s |' % (r['id'], r['alleles'], r['longest'], fmt(r['ref']),
                                                          fmt(r['g_ref']), fmt(r['pair']), fmt(r['g_pair'])))


def versus(la, ta, lb, tb):
    common = sorted(set(ta) & set(tb))
    print('\n## %s vs %s on %d shared regions (%s only %d, %s only %d)' % (
        la, lb, len(common), la, len(set(ta) - set(tb)), lb, len(set(tb) - set(ta))))
    print('\n| kind | F1 %s | F1 %s | regions %s better by >%.2f | %s better |' % (la, lb, la, MARGIN, lb))
    print('|---|---|---|---|---|')
    for kind in ('pair', 'ref'):
        A, B = [ta[i] for i in common], [tb[i] for i in common]
        both = [(x[kind], y[kind]) for x, y in zip(A, B) if x[kind] is not None and y[kind] is not None]
        print('| %s | %.4f | %.4f | %d | %d |' % (kind, sums(A, kind)[2], sums(B, kind)[2],
                                                   sum(x > y + MARGIN for x, y in both), sum(y > x + MARGIN for x, y in both)))


def main(argv):
    tables = []
    for a in argv:
        label, path = a.split('=', 1)
        tables.append((label, load(path)))
    for label, t in tables:
        if t and 'graph_sequence' in next(iter(t.values())):
            before_after(label, t)
    for i in range(len(tables)):
        for j in range(i + 1, len(tables)):
            versus(*tables[i], *tables[j])


if __name__ == '__main__':
    main(sys.argv[1:])
