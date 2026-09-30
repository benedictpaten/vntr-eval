#!/usr/bin/env python3
"""Paired Stage 0 comparison of every variant against poa_abpoa__all on the test set (work/iterate/results.tsv)."""
import csv, statistics as st, sys
rows = list(csv.DictReader(open('/Users/benedictpaten/PycharmProjects/vntr-eval/work/iterate/results.tsv'), delimiter='\t'))
d = {}
for r in rows:
    d.setdefault(r['variant'], {})[r['region_id']] = r
base, hp = d['poa_abpoa__all'], d['poa_abpoa']
ids = list(base)
f = lambda v, r, k: float(d[v][r][k]) if d[v][r][k] else None
want = sys.argv[1].split(',') if len(sys.argv) > 1 else list(d)
print('%-14s %7s %7s %5s %5s %7s %5s %5s %6s  %s' % ('variant', 'kmerx', 'd_med', 'b', 'w', 'cost', 'b', 'w', 'gap%', 'mean kmerx by group hurt/helped/hotspot/mc_right'))
for v in want:
    if v not in d or f(v, ids[0], 'kmer_frac_extra') is None:
        continue
    dk = [f(v, r, 'kmer_frac_extra') - f('poa_abpoa__all', r, 'kmer_frac_extra') for r in ids]
    dc = [f(v, r, 'cost_over_opt') - f('poa_abpoa__all', r, 'cost_over_opt') for r in ids]
    mv = st.mean(f(v, r, 'kmer_frac_extra') for r in ids)
    ma = st.mean(f('poa_abpoa__all', r, 'kmer_frac_extra') for r in ids)
    mp = st.mean(f('poa_abpoa', r, 'kmer_frac_extra') for r in ids)
    grp = []
    for g in ('hurt', 'helped', 'hotspot', 'mc_right'):
        grp.append('%.3f' % st.mean(f(v, r, 'kmer_frac_extra') for r in ids if base[r]['group'] == g))
    print('%-14s %7.4f %7.4f %5d %5d %7.4f %5d %5d %5.0f%%  %s' % (v, mv, st.median(dk), sum(x < -0.005 for x in dk), sum(x > 0.005 for x in dk),
          st.mean(f(v, r, 'cost_over_opt') for r in ids), sum(x < -0.005 for x in dc), sum(x > 0.005 for x in dc), 100 * (ma - mv) / (ma - mp), '/'.join(grp)))
