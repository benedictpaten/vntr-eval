#!/usr/bin/env python3
"""Depth ratio (DR = effective reads / expected reads of the written genotype) and read count (DP) of each
covered repeat site, from an arm's VCF, against what the true pair's lengths would expect.

DR_truth ~ DR_called * sum(len_called + Lbar - 1) / sum(len_truth_best + Lbar - 1), with allele lengths over the
region span (the site's interior differs from the span only by flank shared by every allele).

    python3 work/iterate/wrongwalk/depth_ratio.py [arm]   (default rep_site)
"""
import csv, json, math, statistics, subprocess, sys, collections
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
ARM = sys.argv[1] if len(sys.argv) > 1 else 'rep_site'
VCF = {'rep_site': 'chr20.rep', 'rep_d1': 'chr20.rep_d1', 'rep_raw': 'chr20.rep_raw', 'rep_d1raw': 'chr20.rep_d1raw'}[ARM]
V = R + '/work/stage4/chr20/patched_st_medoid/%s.vcf.gz' % VCF
LBAR = 150.0
cls = json.load(open(R + '/work/iterate/wrongwalk/classify_%s.json' % ARM)) if False else None
rows = {r['region']: r for r in csv.DictReader(open(R + '/work/iterate/wrongwalk/classify.tsv'), delimiter='\t')}
out = collections.defaultdict(list)
for rid, r in rows.items():
    j = json.load(open(R + '/work/stage4/regions/%s/region.json' % rid))
    a, b = int(j['span_start']), int(j['span_end'])
    q = subprocess.run(['bcftools', 'query', '-r', 'chr20:%d-%d' % (a, b), '-f', '%POS\t%REF\t[%DP\t%DR]\n', V],
                       capture_output=True, text=True).stdout
    best = None
    for l in q.splitlines():
        p, ref, dp, dr = l.split('\t')
        if dp in ('.', '') or dr in ('.', ''):
            continue
        dp, dr = int(dp), float(dr)
        if best is None or dp > best[0]:
            best = (dp, dr)
    if best is None:
        continue
    dp, dr = best
    cl = [int(x) for x in r['called_len'].split('/')]
    tl = [int(x) for x in r['ceil_len'].split('/')]
    mu_ratio = sum(x + LBAR - 1 for x in cl) / sum(x + LBAR - 1 for x in tl)
    wrong = int(r['n_wrong']) > 0
    if not wrong:
        g = 'right walks'
    else:
        d = sum(cl) - sum(tl)
        g = 'wrong, called shorter' if d < -10 else ('wrong, called longer' if d > 10 else 'wrong, |total dlen| <= 10')
    out[g].append((dr, dr * mu_ratio, dp, int(r['excess']), rid))
print('arm', ARM, '(DR of the site record with the most reads; DR_truth estimated with Lbar = %d)' % LBAR)
for g in ['right walks', 'wrong, |total dlen| <= 10', 'wrong, called shorter', 'wrong, called longer']:
    xs = out[g]
    if not xs:
        continue
    print('  %-26s %3d sites  DR median %.3f [IQR %.3f-%.3f]  DR_truth median %.3f  DP median %d  excess %d' % (
        g, len(xs), statistics.median(x[0] for x in xs), *statistics.quantiles([x[0] for x in xs], n=4)[::2],
        statistics.median(x[1] for x in xs), statistics.median(x[2] for x in xs), sum(x[3] for x in xs)))
for g in ['wrong, called shorter', 'wrong, called longer']:
    xs = sorted(out[g], key=lambda x: -x[3])[:8]
    print('  %s, largest:' % g, '  '.join('%s DR %.2f->truth %.2f (DP %d, excess %d)' % (x[4], x[0], x[1], x[2], x[3]) for x in xs))
