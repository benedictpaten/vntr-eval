#!/usr/bin/env python3
"""Offline, from dump.tsv: at every covered repeat site, the best genotype under vg's length-weighted mixture
(vg's own choice) and under flat weights (with vg's depth term), and whether each gets both allele lengths right
(the ceiling pair's span lengths; alleles mapped to span lengths by the offset the called pair fixes)."""
import collections, csv, math
import numpy as np
cls = {r['region']: r for r in csv.DictReader(open('classify.tsv'), delimiter='\t')}
sites = collections.defaultdict(lambda: {'reads': [], 'gt': {}})
for l in open('dump.tsv'):
    f = l.rstrip('\n').split('\t')
    rid, tag = f[0], f[1]
    if tag == '#lengths':
        sites[(rid, f[2])]['T'] = [int(x) for x in f[3:]]
    elif tag == '#gt':
        a, b = (int(x) for x in f[3].split('/'))
        sites[(rid, f[2])]['gt'][(a, b)] = (float(f[4]), float(f[5]), float(f[6]))
    elif tag == '#depth':
        pass
    else:
        sites[(rid, tag)]['reads'].append([float(f[3])] + [float(x) for x in f[5:]])
C = collections.Counter(); ex = collections.Counter()
for (rid, site), s in sites.items():
    r = cls[rid]
    if not s['reads'] or 'T' not in s:
        continue
    M = np.array(s['reads']); e = M[:, 0]; rel = M[:, 1:]
    n = rel.shape[1]; T = s['T']; g = s['gt']
    vgbest = max(g, key=lambda k: g[k][0])
    best, bk = -1e300, None
    for i in range(n):
        for j in range(i, n):
            v = np.log((1 - e) * (0.5 * rel[:, i] + 0.5 * rel[:, j]) + e).sum() + g[(i, j)][2]
            if v > best:
                best, bk = v, (i, j)
    called = sorted(int(x) for x in r['called_len'].split('/'))
    right = sorted(int(x) for x in r['ceil_len'].split('/'))
    offs = {called[0] - T[vgbest[0]]} if sorted([called[0] - T[vgbest[0]], called[1] - T[vgbest[1]]]) and \
        called[0] - T[vgbest[0]] == called[1] - T[vgbest[1]] else set()
    if not offs:
        for p in ((0, 1), (1, 0)):
            if called[p[0]] - T[vgbest[0]] == called[p[1]] - T[vgbest[1]]:
                offs = {called[p[0]] - T[vgbest[0]]}
    if not offs:
        C['offset?'] += 1; continue
    off = offs.pop()
    vl = sorted(T[k] + off for k in vgbest); fl = sorted(T[k] + off for k in bk)
    grp = 'right walks' if int(r['n_wrong']) == 0 else 'wrong walks'
    C[(grp, 'vg lengths right', vl == right)] += 1
    C[(grp, 'flat lengths right', fl == right)] += 1
    C[(grp, 'flat changes the genotype', bk != vgbest)] += 1
for k in sorted(C, key=str):
    print(k, C[k])
