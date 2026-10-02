#!/usr/bin/env python3
"""Called haplotypes of every chr20 repeat region, from each arm's whole-contig VCF (vg 91d38c802), scored
against HG002's truth with tools/score_haplotypes.py (no truvari), with the two haplotypes kept as FASTA.

    python3 work/iterate/hapscore_chr20.py [--jobs 4]   -> work/iterate/hapscore/<arm>/<id>.{json,fa}
"""
import argparse, concurrent.futures, os, subprocess, sys
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
C = R + '/work/stage4/chr20'
ARMS = {'mc': C + '/rb/unpatched/chr20.vcf.gz', 'abpoa_h32': C + '/rb/patched_all/chr20.vcf.gz',
        'st_chm13': C + '/patched_st_chm13/chr20.vcf.gz', 'st_medoid': C + '/patched_st_medoid/chr20.vcf.gz',
        'st_medoid_chains': C + '/patched_st_medoid/chr20.chains.vcf.gz',
        'rep_site': C + '/patched_st_medoid/chr20.rep.vcf.gz', 'rep_desc': C + '/patched_st_medoid/chr20.desc.vcf.gz',
        'rep_link': C + '/patched_st_medoid/chr20.link.vcf.gz'}
OUT = R + '/work/iterate/hapscore'


def one(args):
    arm, rid = args
    od = os.path.join(OUT, arm)
    js = os.path.join(od, rid + '.json')
    if os.path.exists(js):
        return arm, rid, 'cached'
    p = subprocess.run([sys.executable, R + '/tools/score_haplotypes.py', 'score', R + '/work/stage4/regions/' + rid,
                        ARMS[arm], '--label', arm, '--no-truvari', '--no-refine', '--out', js,
                        '--haplotypes', os.path.join(od, rid + '.fa')], capture_output=True, text=True)
    return arm, rid, 'ok' if p.returncode == 0 else p.stderr[-200:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=4)
    a = ap.parse_args()
    ids = open(C + '/patch_all_ids.txt').read().split()
    for arm in ARMS:
        if not os.path.exists(ARMS[arm]):
            continue
        os.makedirs(os.path.join(OUT, arm), exist_ok=True)
    bad = 0
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
        for arm, rid, st in ex.map(one, [(arm, r) for arm in ARMS if os.path.exists(ARMS[arm]) for r in ids]):
            if st not in ('ok', 'cached'):
                bad += 1
                print(arm, rid, st.replace('\n', ' '))
    print('done, failures', bad)


if __name__ == '__main__':
    main()
