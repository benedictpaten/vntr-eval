#!/usr/bin/env python3
"""Fragmentation of the hap32 rows in each variant's projected MSA: gap runs summed over distinct rows,
distinct gap runs (start, end column), and top-level sites. Usage: frag.py v1,v2,... r1,r2,..."""
import os, sys
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
sys.path.insert(0, R + '/tools')
import msa_graph, iterate

def runs(row):
    out, i, n = [], 0, len(row)
    while i < n:
        if row[i] == '-':
            j = i
            while j < n and row[j] == '-':
                j += 1
            out.append((i, j))
            i = j
        else:
            i += 1
    return out

def one(v, rid):
    f = '%s/work/iterate/candidates/%s/%s.msa.fa' % (R, v, rid)
    if not os.path.exists(f):
        f = '%s/work/stage4/candidates/%s/%s.msa.fa' % (R, v, rid)
    if not os.path.exists(f):
        return None
    rows = {}
    for n, r in msa_graph.read_msa(f):
        if n not in msa_graph.CONSENSUS_NAMES:
            rows.setdefault(r.replace('-', ''), r)
    # gap runs measured between bases only (leading and trailing gaps are ignored)
    tot, distinct = 0, set()
    for r in rows.values():
        for a, b in runs(r):
            if a > 0 and b < len(r):
                tot += 1
                distinct.add((a, b))
    g = f[:-len('.msa.fa')] + '.gfa'
    s = iterate.snarl_counts(g).get('sites') if os.path.exists(g) else None
    return tot, len(distinct), s

def main():
    vs, rids = sys.argv[1].split(','), sys.argv[2].split(',')
    print('%-10s' % 'region' + ''.join('%22s' % v[:21] for v in vs))
    T = {v: [0, 0, 0] for v in vs}
    for rid in rids:
        line = '%-10s' % rid
        for v in vs:
            x = one(v, rid)
            line += '%22s' % ('-' if x is None else '%d/%d/%s' % x)
            if x:
                T[v] = [a + (b or 0) for a, b in zip(T[v], x)]
        print(line)
    print('%-10s' % 'total' + ''.join('%22s' % ('%d/%d/%d' % tuple(T[v])) for v in vs))
    print('cells: gap runs over distinct hap32 rows / distinct gap runs / top-level sites')


if __name__ == '__main__':
    main()
