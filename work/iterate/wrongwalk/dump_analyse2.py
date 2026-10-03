#!/usr/bin/env python3
"""dump_analyse.py with the offset validated: the offset is the one value for which every vg allele's T_h plus it is
a hap32 allele length (so vg's alleles are the panel's). Then vg's best genotype against the best genotype of the
right lengths, split into read and depth terms, and the read term again under flat weights."""
import collections, csv, sys, math
import numpy as np
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
sys.path.insert(0, R + '/tools'); import msa_graph
cls = {r['region']: r for r in csv.DictReader(open('classify.tsv'), delimiter='\t')}
S = collections.defaultdict(lambda: {'reads': [], 'gt': {}})
for l in open('dump.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[1] == '#lengths': S[(f[0], f[2])]['T'] = [int(x) for x in f[3:]]
    elif f[1] == '#gt':
        a, b = (int(x) for x in f[3].split('/')); S[(f[0], f[2])]['gt'][(a, b)] = tuple(float(x) for x in f[4:7])
    elif f[1] != '#depth': S[(f[0], f[1])]['reads'].append([float(f[3])] + [float(x) for x in f[5:]])
C = collections.Counter(); rows = []
for (rid, site), s in S.items():
    r = cls[rid]
    if int(r['n_wrong']) == 0 or 'T' not in s or not s['reads']:
        continue
    H = {len(x) for _, x in msa_graph.read_fasta(R + '/work/stage4/regions/%s/hap32.fa' % rid)}
    T = s['T']
    offs = [o for o in {h - T[0] for h in H} if all(t + o in H for t in T)]
    if len(offs) != 1:
        C['offset ambiguous or none (%d)' % len(offs)] += 1; continue
    off = offs[0]; g = s['gt']
    right = [int(x) for x in r['ceil_len'].split('/')]
    cand = [[i for i, t in enumerate(T) if t + off == x] for x in right]
    if not cand[0] or not cand[1]:
        C['right length missing'] += 1; rows.append((rid, 'missing', right, sorted(set(t + off for t in T)))); continue
    best = max(g, key=lambda k: g[k][0])
    rg = max(((min(i, j), max(i, j)) for i in cand[0] for j in cand[1]), key=lambda k: g[k][0])
    if rg == best:
        C['vg best has the right lengths (error is sequence or later)'] += 1; continue
    M = np.array(s['reads']); e = M[:, 0]; rel = M[:, 1:]
    flat = lambda k: float(np.log((1 - e) * (0.5 * rel[:, k[0]] + 0.5 * rel[:, k[1]]) + e).sum())
    d_read, d_depth, d_flat = g[best][1] - g[rg][1], g[best][2] - g[rg][2], flat(best) - flat(rg)
    # per read: does vg's rel favour the right genotype's alleles or the called's (best fit within each genotype)
    fr = np.maximum(rel[:, rg[0]], rel[:, rg[1]]); fc = np.maximum(rel[:, best[0]], rel[:, best[1]])
    k = ('read term favours called' if d_read > 0 else 'read term favours right') + \
        (' | flat weights too' if d_flat > 0 else ' | flat weights favour right')
    C[k] += 1
    rows.append((rid, int(r['excess']), round(d_read, 1), round(d_depth, 1), round(d_flat, 1), int((fr > fc).sum()), int((fc > fr).sum()), len(e)))
for k, v in C.most_common(): print('%-70s %d' % (k, v))
print('\nregion excess d_read d_depth d_read_flat reads_favouring_right reads_favouring_called n_reads')
for x in sorted([x for x in rows if isinstance(x[1], int)], key=lambda x: -x[1])[:25]: print(x)
print([x for x in rows if x[1] == 'missing'])
