#!/usr/bin/env python3
"""Called-haplotype edit distance on the covered regions for each --repeat-sites arm, split by what
rep_site got wrong there (classify.tsv), with paired bootstraps against rep_site and the prototype.

    python3 work/iterate/wrongwalk/compare_arms.py rep_d1 rep_raw rep_d1raw [vgreads ...]
An arm named vgreads / vgreads_nodepth is read from rescore.tsv.
"""
import csv, json, os, random, sys
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
W = R + '/work/iterate/wrongwalk'
rows = {r['region']: r for r in csv.DictReader(open(W + '/classify.tsv'), delimiter='\t')}
cmp = json.load(open(R + '/work/iterate/atomic/compare_chr20.json'))
resc = {r['region']: r for r in csv.DictReader(open(W + '/rescore.tsv'), delimiter='\t')} if os.path.exists(W + '/rescore.tsv') else {}
arms = sys.argv[1:]


def ed(arm, rid):
    if arm in ('vgreads', 'vgreads_nodepth'):
        r = resc.get(rid)
        return int(r[arm]) if r else None
    p = R + '/work/iterate/hapscore/%s/%s.json' % (arm, rid)
    if not os.path.exists(p):
        return None
    d = json.load(open(p))
    return d['ed'] if d.get('status') == 'ok' else None


def group(r):
    if int(r['n_wrong']) == 0:
        return 'right walks'
    cl = sum(int(x) for x in r['called_len'].split('/'))
    tl = sum(int(x) for x in r['ceil_len'].split('/'))
    return 'called shorter' if cl - tl < -10 else ('called longer' if cl - tl > 10 else 'same total length')


vals = {}
for rid, r in rows.items():
    v = {'ceiling': int(r['ceil']), 'atomic': int(r['atomic']), 'rep_site': int(r['ed'])}
    ok = True
    for a in arms:
        x = ed(a, rid)
        if x is None:
            ok = False
            break
        v[a] = x
    if ok:
        vals[rid] = (group(r), v)
cols = ['ceiling', 'atomic', 'rep_site'] + arms
print('covered regions scored by every arm: %d' % len(vals))
print('  %-18s %4s ' % ('rep_site got', 'n') + ' '.join('%10s' % c for c in cols))
for g in ['right walks', 'same total length', 'called shorter', 'called longer', None]:
    ks = [k for k, (gg, _) in vals.items() if g is None or gg == g]
    print('  %-18s %4d ' % (g or 'all', len(ks)) + ' '.join('%10d' % sum(vals[k][1][c] for k in ks) for c in cols))
random.seed(1)
print('\npaired bootstrap over regions (negative: first is better)')
for a in arms:
    for b in ('rep_site', 'atomic'):
        d = [v[a] - v[b] for _, v in vals.values()]
        bs = sorted(sum(random.choice(d) for _ in d) for _ in range(2000))
        print('  %-16s - %-9s %+7d [%+d, %+d]  better/worse %d/%d' % (a, b, sum(d), bs[50], bs[1949], sum(x < 0 for x in d), sum(x > 0 for x in d)))
