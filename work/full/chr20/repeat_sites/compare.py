#!/usr/bin/env python3
"""Repeat-sites arms against base on the latest chr20 graph (run.sh): F1, the per-region paired bootstrap of SV FP+FN
inside the realigned regions (patch.json spans +- 100 bp), raw and refined, and called-haplotype edit distance and
off-panel strands (called_arms.tsv). usage: compare.py"""
import bisect, collections, csv, gzip, json, random
O = '/Users/benedictpaten/PycharmProjects/vntr-eval/work/full/chr20/repeat_sites'
T = '/Users/benedictpaten/PycharmProjects/vntr-eval/work/full/chr20/test'
ARMS = ['base', 'rep', 'desc', 'link']
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


print('| arm | ALL F1 | SNV F1 | indel F1 | SV F1 raw (FP / FN) | SV F1 refined (FP / FN) |')
print('|---|---|---|---|---|---|')
E = {}
for arm in ARMS:
    log = open('%s/score.%s.log' % (O, arm)).read()
    i = log.index('{', log.index('score: {'))
    s = json.loads(log[i:log.rindex('}') + 1])
    d = '%s/score/rs_%s/score/chr20.truvari/' % (O, arm)
    r = json.load(open(d + 'refine.variant_summary.json'))
    print('| %s | %.4f | %.4f | %.4f | %.4f (%d / %d) | %.4f (%d / %d) |' % (
        arm, s['all_f1'], s['snv_f1'], s['indel_f1'], s['sv_f1'], s['sv_fp'], s['sv_fn'], r['f1'], r['FP'], r['FN']))
    E[arm] = {'raw': (raw(d + 'fp.vcf.gz'), raw(d + 'fn.vcf.gz')),
              'refined': (refined(d + 'refine.comp.vcf.gz', 'FP'), refined(d + 'refine.base.vcf.gz', 'FN'))}
random.seed(1)
print()
for arm in ARMS[1:]:
    for kind in ('raw', 'refined'):
        keys = set()
        for a in ('base', arm):
            for c in E[a][kind]:
                keys |= set(c)
        keys.discard('outside')
        diff = {k: (E[arm][kind][0][k] + E[arm][kind][1][k]) - (E['base'][kind][0][k] + E['base'][kind][1][k]) for k in keys}
        d = list(diff.values())
        bs = sorted(sum(random.choice(d) for _ in d) for _ in range(2000))
        oa = E['base'][kind][0]['outside'] + E['base'][kind][1]['outside']
        ob = E[arm][kind][0]['outside'] + E[arm][kind][1]['outside']
        print('%-5s %-8s SV FP+FN in regions, arm - base: %+d [%+d, %+d]; regions better/worse %d/%d; outside %d -> %d' % (
            arm, kind, sum(d), bs[50], bs[1949], sum(v < 0 for v in d), sum(v > 0 for v in d), oa, ob))
# called haplotypes: summed edit distance to the paired truth and off-panel strands, regions inside the truth BED
rows = [r for r in csv.DictReader(open(O + '/called_arms.tsv'), delimiter='\t') if r['in_bed'] == '1' and r['skipped'] == '0']
by = collections.defaultdict(dict)
for r in rows:
    by[(r['region'], r['strand'])][r['arm']] = r
common = [k for k, v in by.items() if all(a in v and v[a]['ed_truth'].isdigit() for a in ARMS)]
print('\ncalled haplotypes, %d strands in %d regions inside the truth BED, scored in every arm:' % (
    len(common), len({k[0] for k in common})))
print('| arm | edit distance to truth | off-panel strands |')
print('|---|---|---|')
for arm in ARMS:
    ed = sum(int(by[k][arm]['ed_truth']) for k in common)
    off = sum(by[k][arm]['in_panel'] == '0' for k in common)
    print('| %s | %d | %d |' % (arm, ed, off))
for arm in ARMS[1:]:
    regs = collections.defaultdict(int)
    for k in common:
        regs[k[0]] += int(by[k][arm]['ed_truth']) - int(by[k]['base']['ed_truth'])
    d = list(regs.values())
    bs = sorted(sum(random.choice(d) for _ in d) for _ in range(2000))
    print('%-5s edit distance, arm - base: %+d [%+d, %+d]; regions better/worse %d/%d' % (
        arm, sum(d), bs[50], bs[1949], sum(v < 0 for v in d), sum(v > 0 for v in d)))
