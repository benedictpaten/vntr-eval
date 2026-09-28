#!/usr/bin/env python3
"""stage4_detect.py -- Stage 4 step D: detection quality against production vg SV errors on one contig.

    python3 tools/stage4_detect.py --scan work/stage4/scan.tsv --hotspots work/stage4/hotspots.tsv \
        --truvari ~/PycharmProjects/vg-call-eval/work/wgs-mm095/score/chr20.truvari \
        --out work/stage4/detection.tsv [--json work/stage4/detection.json]

Each truvari FP (comp) and FN (base) record [POS-1, POS-1+len(REF)) is assigned to every scanned TR
region whose padded interval (scan.tsv pad_start..pad_end) it overlaps. Reported for three nested
sets -- every merged TR region, the prefilter (max length diff >= 50 bp), the detected regions --
the regions, the distinct FP records they hold (share of all TR FPs), and the regions with >= 1 FP.
Standard library only (bcftools for reading the VCFs).
"""
import argparse
import bisect
import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402


def records(vcf):
    out = subprocess.run([config.BCFTOOLS, 'query', '-f', '%CHROM\t%POS\t%REF\n', vcf],
                         capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        c, p, r = line.split('\t')
        yield c, int(p) - 1, int(p) - 1 + len(r)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scan', required=True)
    ap.add_argument('--hotspots', required=True)
    ap.add_argument('--truvari', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--json')
    a = ap.parse_args(argv)
    regs = list(csv.DictReader(open(a.scan), delimiter='\t'))
    contig = regs[0]['contig']
    regs.sort(key=lambda r: int(r['pad_start']))
    starts = [int(r['pad_start']) for r in regs]
    maxlen = max(int(r['pad_end']) - int(r['pad_start']) for r in regs)
    det = {r['region_id']: r for r in csv.DictReader(open(a.hotspots), delimiter='\t')}
    hits = {r['region_id']: {'fp': set(), 'fn': set()} for r in regs}
    n_total = {}
    for kind in ('fp', 'fn'):
        recs = [x for x in records(os.path.join(a.truvari, kind + '.vcf.gz')) if x[0] == contig]
        n_total[kind] = len(recs)
        for k, (_, s, e) in enumerate(recs):
            i = bisect.bisect_left(starts, s - maxlen)
            while i < len(regs) and starts[i] < e:
                if int(regs[i]['pad_end']) > s:
                    hits[regs[i]['region_id']][kind].add(k)
                i += 1
    sets = {'all_tr': [r['region_id'] for r in regs],
            'prefilter': [r['region_id'] for r in regs if r['prefilter'] == '1'],
            'packaged_stage0': [k for k, d in det.items() if d['status'] == 'ok'],
            'detected': [k for k, d in det.items() if d['detected'] == '1']}
    allfp = set().union(*(hits[r]['fp'] for r in sets['all_tr']))
    allfn = set().union(*(hits[r]['fn'] for r in sets['all_tr']))
    rep = {'contig': contig, 'truvari_fp_total': n_total['fp'], 'truvari_fn_total': n_total['fn'],
           'tr_fp': len(allfp), 'tr_fn': len(allfn), 'sets': {}}
    for name, ids in sets.items():
        fp = set().union(*(hits[r]['fp'] for r in ids)) if ids else set()
        fn = set().union(*(hits[r]['fn'] for r in ids)) if ids else set()
        nfp = sum(1 for r in ids if hits[r]['fp'])
        rep['sets'][name] = {'regions': len(ids), 'fp_records': len(fp),
                             'share_of_tr_fp': round(len(fp) / max(1, len(allfp)), 4),
                             'fn_records': len(fn), 'share_of_tr_fn': round(len(fn) / max(1, len(allfn)), 4),
                             'regions_with_fp': nfp, 'share_regions_with_fp': round(nfp / max(1, len(ids)), 4)}
    with open(a.out, 'w') as fo:
        fo.write('region_id\ttr_class\tprefilter\tstage0\tdetected\tvg_sv_fp\tvg_sv_fn\n')
        for r in regs:
            if r['prefilter'] != '1':
                continue
            d = det.get(r['region_id'], {})
            fo.write('\t'.join([r['region_id'], r['tr_class'], r['prefilter'], d.get('status', 'missing'),
                                d.get('detected', '0'), str(len(hits[r['region_id']]['fp'])),
                                str(len(hits[r['region_id']]['fn']))]) + '\n')
    if a.json:
        json.dump(rep, open(a.json, 'w'), indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == '__main__':
    main()
