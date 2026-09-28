#!/usr/bin/env python3
"""stage4_report.py -- tables for results/stage4_chr20.md from patch_contig.py's score directories.

    python3 tools/stage4_report.py --score work/stage4/chr20/score --patch work/stage4/chr20/patched/patch.json \
        --detection work/stage4/detection.tsv --out results/stage4_chr20

-> <out>.tsv (one row per arm: chr20 SV and small-variant metrics, SV FP/FN inside/outside the patched
regions, aardvark-labelled small-variant FP/FN records inside/outside, summed haplotype ED over the
patched regions split by production SV FP > 0 / = 0) and <out>_regions.tsv (per patched region: its
production SV FP/FN and each arm's haplotype ED)."""
import argparse
import csv
import gzip
import json
import os

ARMS = ('production', 'unpatched', 'mcnull', 'patched')


def inside(ivs, pos, end):
    return any(s <= end and pos <= e for s, e in ivs)


def bd_counts(path, label, ivs):
    """(snv_in, snv_out, indel_in, indel_out) of aardvark records whose BD is `label`."""
    c = [0, 0, 0, 0]
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            x = line.rstrip('\n').split('\t')
            fmt = x[8].split(':')
            if 'BD' not in fmt or x[9].split(':')[fmt.index('BD')] != label:
                continue
            pos, ref = int(x[1]), x[3]
            snv = len(ref) == 1 and all(len(a) == 1 for a in x[4].split(','))
            k = (0 if snv else 2) + (0 if inside(ivs, pos, pos + len(ref) - 1) else 1)
            c[k] += 1
    return c


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--score', required=True)
    p.add_argument('--patch', required=True)
    p.add_argument('--detection', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--pad', type=int, default=100)
    a = p.parse_args()
    pj = json.load(open(a.patch))
    ids = sorted(pj['regions'])
    ivs = [(g['span'][0] - a.pad, g['span'][1] + a.pad) for g in pj['regions'].values()]
    det = {r['region_id']: r for r in csv.DictReader(open(a.detection), delimiter='\t')}
    rows, per = [], {i: {'region_id': i, 'prod_sv_fp': int(det[i]['vg_sv_fp']), 'prod_sv_fn': int(det[i]['vg_sv_fn']),
                         'span_bp': pj['regions'][i]['span'][1] - pj['regions'][i]['span'][0] + 1} for i in ids}
    for arm in ARMS:
        d = os.path.join(a.score, arm)
        r = json.load(open(os.path.join(d, 'summary.json')))
        ad = os.path.join(d, 'score', 'chr20.aardvark')
        q = bd_counts(os.path.join(ad, 'query.vcf.gz'), 'FP', ivs)
        t = bd_counts(os.path.join(ad, 'truth.vcf.gz'), 'FN', ivs)
        r.update(small_fp_snv_in=q[0], small_fp_snv_out=q[1], small_fp_indel_in=q[2], small_fp_indel_out=q[3],
                 small_fn_snv_in=t[0], small_fn_snv_out=t[1], small_fn_indel_in=t[2], small_fn_indel_out=t[3])
        eds = {}
        for i in ids:
            h = json.load(open(os.path.join(d, 'haps', i + '.json')))
            assert h['status'] == 'ok', (arm, i, h['status'])
            eds[i] = h['ed']
            per[i]['ed_' + arm] = h['ed']
            per[i]['ed_ref'] = h['ed_ref']
        fp = [i for i in ids if per[i]['prod_sv_fp'] > 0]
        nofp = [i for i in ids if per[i]['prod_sv_fp'] == 0]
        r.update(hap_ed_all=sum(eds.values()), hap_ed_fp_regions=sum(eds[i] for i in fp),
                 hap_ed_nofp_regions=sum(eds[i] for i in nofp), n_fp_regions=len(fp), n_nofp_regions=len(nofp),
                 hap_ed_ref_all=sum(per[i]['ed_ref'] for i in ids))
        rows.append(r)
    keys = list(rows[0])
    with open(a.out + '.tsv', 'w') as f:
        f.write('\t'.join(keys) + '\n')
        for r in rows:
            f.write('\t'.join(('%.4f' % r[k]) if isinstance(r[k], float) else str(r[k]) for k in keys) + '\n')
    cols = ['region_id', 'span_bp', 'prod_sv_fp', 'prod_sv_fn', 'ed_ref'] + ['ed_' + x for x in ARMS]
    with open(a.out + '_regions.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for i in ids:
            f.write('\t'.join(str(per[i][c]) for c in cols) + '\n')
    for r in rows:
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}))


if __name__ == '__main__':
    main()
