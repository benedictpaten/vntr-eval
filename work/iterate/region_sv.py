#!/usr/bin/env python3
"""Per-region chr20 SV FP/FN (inside patched_full spans +-100 bp) for production, unpatched,
patched_full, patched_all, plus full-panel abPOA cost; -> work/iterate/region_sv.tsv"""
import gzip, json, os, subprocess, collections
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
S = R + '/work/stage4/chr20/score_full/%s/score/chr20.truvari/%s.vcf.gz'
P = os.path.expanduser('~/PycharmProjects/vg-call-eval/work/wgs-mm095/score/chr20.truvari/%s.vcf.gz')
patch = json.load(open(R + '/work/stage4/chr20/patched_full/patch.json'))['regions']
spans = {r: v['span'] for r, v in patch.items() if isinstance(v, dict) and 'span' in v}
iv = sorted((s[0] - 100, s[1] + 100, r) for r, s in spans.items())
def count(path):
    c = collections.Counter()
    for line in gzip.open(path, 'rt'):
        if line[0] == '#': continue
        x = line.split('\t', 3); pos = int(x[1])
        for a, b, r in iv:
            if a <= pos <= b: c[r] += 1
    return c
arms = {'prod': lambda k: P % k}
for a in ('unpatched', 'patched_full', 'patched_all'):
    arms[a] = (lambda a: lambda k: S % (a, k))(a)
cnt = {(a, k): count(f(k)) for a, f in arms.items() for k in ('fp', 'fn', 'tp-base')}
out = open(R + '/work/iterate/region_sv.tsv', 'w')
cols = ['region_id', 'span_bp', 'n_distinct', 'full_align_s', 'full_rss_mb', 'hap32_distinct', 'period']
for a in arms: cols += [a + '_fp', a + '_fn']
cols += ['prod_tp', 'd_full', 'd_all']
out.write('\t'.join(cols) + '\n')
for r in sorted(spans):
    j = R + '/work/stage4/candidates/poa_abpoa__all/%s.realign.json' % r
    d = json.load(open(j)) if os.path.exists(j) else {}
    rj = json.load(open(R + '/work/stage4/regions/%s/region.json' % r))
    al = (d.get('align') or {}).get('aligner') or {}
    row = [r, spans[r][1] - spans[r][0] + 1, d.get('n_distinct'), al.get('seconds'), al.get('peak_rss_mb'),
           (rj.get('hap32') or {}).get('n_distinct_sequences'), rj.get('period')]
    e = {}
    for a in arms:
        e[a] = cnt[(a, 'fp')][r] + cnt[(a, 'fn')][r]
        row += [cnt[(a, 'fp')][r], cnt[(a, 'fn')][r]]
    row += [cnt[('prod', 'tp-base')][r], e['patched_full'] - e['unpatched'], e['patched_all'] - e['unpatched']]
    out.write('\t'.join(map(str, row)) + '\n')
out.close()
