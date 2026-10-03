#!/usr/bin/env python3
"""For each covered repeat site in dump.tsv (vg-diag --dump-likelihoods): vg's best genotype against the genotype
of the right panel alleles (the ceiling pair), each split into read term and depth term.

Alleles are identified by length: the span length of an allele is T_h plus a constant offset (the parts of the
span outside the site's interior), fixed by the called pair's span lengths.

    python3 work/iterate/wrongwalk/dump_analyse.py
"""
import collections, csv, math
W = '/Users/benedictpaten/PycharmProjects/vntr-eval/work/iterate/wrongwalk'
cls = {r['region']: r for r in csv.DictReader(open(W + '/classify.tsv'), delimiter='\t')}
L, DEP, GT, NR = {}, {}, collections.defaultdict(dict), collections.Counter()
EVID = collections.defaultdict(list)
for l in open(W + '/dump.tsv'):
    f = l.rstrip('\n').split('\t')
    rid = f[0]
    if f[1] == '#lengths':
        L.setdefault(rid, {})[f[2]] = [int(x) for x in f[3:]]
    elif f[1] == '#depth':
        DEP.setdefault(rid, {})[f[2]] = [float(x) for x in f[3:]]
    elif f[1] == '#gt':
        a, b = (int(x) for x in f[3].split('/'))
        GT[(rid, f[2])][(a, b)] = (float(f[4]), float(f[5]), float(f[6]))
    else:
        NR[(rid, f[1])] += 1
        EVID[(rid, f[1])].append((float(f[3]), [float(x) for x in f[5:]]))
out = collections.Counter()
rows = []
for (rid, site), g in GT.items():
    r = cls[rid]
    if int(r['n_wrong']) == 0:
        continue
    T = L[rid][site]
    best = max(g, key=lambda k: g[k][0])
    called = [int(x) for x in r['called_len'].split('/')]
    right = [int(x) for x in r['ceil_len'].split('/')]
    # offset from the best genotype against the called span lengths (either pairing)
    offs = {called[0] - T[best[0]], called[1] - T[best[1]]} if called[0] - T[best[0]] == called[1] - T[best[1]] else \
           {called[0] - T[best[1]]} if called[0] - T[best[1]] == called[1] - T[best[0]] else set()
    if len(offs) != 1:
        out['offset not determined'] += 1
        rows.append((rid, 'offset?', T, best, called, right)); continue
    off = offs.pop()
    cand = [[i for i, t in enumerate(T) if t + off == x] for x in right]
    if not cand[0] or not cand[1]:
        out['a right allele length is not among vg\'s alleles'] += 1
        rows.append((rid, 'right allele missing (len %s; vg lengths %s)' % (right, sorted(set(t + off for t in T))), int(r['excess'])))
        continue
    # among the alleles of the right lengths, the best-scoring genotype (ties of length are rare)
    rg = max(((min(i, j), max(i, j)) for i in cand[0] for j in cand[1]), key=lambda k: g[k][0])
    b, t = g[best], g[rg]
    if rg == best:
        out['right genotype is vg\'s best (the error is length-ambiguity or later)'] += 1
        continue
    d_read, d_depth = b[1] - t[1], b[2] - t[2]
    key = 'read term favours the called' if d_read > 0 else 'read term favours the right'
    key += '; depth favours the called' if d_depth > 0 else '; depth favours the right'
    out[key] += 1
    rows.append((rid, key, int(r['excess']), round(b[0] - t[0], 1), round(d_read, 1), round(d_depth, 1), NR[(rid, site)], best, rg))
for k, v in out.most_common():
    print('%-80s %d' % (k, v))
print()
for x in sorted(rows, key=lambda x: -x[2] if isinstance(x[2], int) else 0)[:30]:
    print(x)
