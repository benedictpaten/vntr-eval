#!/usr/bin/env python3
"""Paired per-region bootstrap of SV errors, realigned vs original full chr20 graph (test_chain.sh arms).
Regions: the realigned regions (patch.json spans +- 100 bp); a call outside every region is 'outside'.
Raw: truvari bench fp/fn VCFs. Refined: truvari refine's BD (FP in refine.comp, FN in refine.base)."""
import bisect, collections, gzip, json, random, sys
T = '/Users/benedictpaten/PycharmProjects/vntr-eval/work/full/chr20/test'
spans = sorted((v['span'][0] - 100, v['span'][1] + 100, r) for r, v in json.load(open(T + '/patch.json'))['regions'].items())
starts = [s[0] for s in spans]

def region(pos):
    i = bisect.bisect_right(starts, pos) - 1
    while i >= 0 and spans[i][0] <= pos:
        if pos <= spans[i][1]:
            return spans[i][2]
        i -= 1
    return 'outside'

def raw(path):
    c = collections.Counter()
    for l in gzip.open(path, 'rt'):
        if l[0] != '#':
            c[region(int(l.split('\t', 2)[1]))] += 1
    return c

def refined(path, want):
    c = collections.Counter()
    for l in gzip.open(path, 'rt'):
        if l[0] == '#':
            continue
        x = l.rstrip('\n').split('\t')
        if x[9].split(':')[x[8].split(':').index('BD')] == want:
            c[region(int(x[1]))] += 1
    return c

random.seed(1)
for rd in ('short', 'ont'):
    E = {}
    for arm in ('original', 'realigned'):
        d = '%s/score/%s_%s/score/chr20.truvari/' % (T, arm, rd)
        E[arm] = {'raw': (raw(d + 'fp.vcf.gz'), raw(d + 'fn.vcf.gz')),
                  'refined': (refined(d + 'refine.comp.vcf.gz', 'FP'), refined(d + 'refine.base.vcf.gz', 'FN'))}
    for kind in ('raw', 'refined'):
        keys = set()
        for arm in E:
            for c in E[arm][kind]:
                keys |= set(c)
        keys.discard('outside')
        diff = {k: (E['realigned'][kind][0][k] + E['realigned'][kind][1][k]) - (E['original'][kind][0][k] + E['original'][kind][1][k]) for k in keys}
        d = list(diff.values())
        bs = sorted(sum(random.choice(d) for _ in d) for _ in range(2000))
        out_o = E['original'][kind][0]['outside'] + E['original'][kind][1]['outside']
        out_r = E['realigned'][kind][0]['outside'] + E['realigned'][kind][1]['outside']
        top = sorted(diff.items(), key=lambda kv: kv[1])
        tot = sum(d)
        ten = sum(v for _, v in top[:10])
        print('%-5s %-7s inside regions: realigned - original FP+FN %+d [%+d, %+d], better/worse %d/%d (regions with errors %d); outside %d -> %d; the 10 most-improved regions hold %+d' % (
            rd, kind, tot, bs[50], bs[1949], sum(v < 0 for v in d), sum(v > 0 for v in d), len(d), out_o, out_r, ten))
        print('       most improved: %s' % ', '.join('%s %+d' % kv for kv in top[:6]))
        print('       most worsened: %s' % ', '.join('%s %+d' % kv for kv in top[::-1][:6]))
