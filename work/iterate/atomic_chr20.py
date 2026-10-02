#!/usr/bin/env python3
"""Genotype every chr20 repeat region as one site whose alleles are the distinct hap32 sequences
(tools/panel_genotype.py, 4j's prototype), with the reads the production mapping placed within 5 kb
of the region (4j used a 200 kb window, which only the flank depth estimate needs 2 kb of).

    python3 work/iterate/atomic_chr20.py [--jobs 4] [--lam 3]   -> work/iterate/atomic/atomic_chr20.tsv
PANELGT_DEPTH_W sets the depth weight (default 1, as 4j).
"""
import argparse, concurrent.futures, os, sys
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
sys.path.insert(0, R + '/tools')
os.environ.setdefault('VNTR_REGIONS', R + '/work/stage4/regions')
import call_local, panel_genotype as pg  # noqa: E402

FLANK = 5000
OUT = R + '/work/iterate/atomic'


def region_reads(rid):
    wd = call_local.fetch_window_reads(rid, FLANK)
    wd = wd if isinstance(wd, str) and os.path.isdir(wd) else os.path.join(
        call_local.STAGE3, 'hybrid%dk' % (FLANK // 1000), 'reads', rid)
    reads = []
    for f in ('reads_1.fq', 'reads_2.fq', 'reads_se.fq'):
        reads += pg.read_fastq(os.path.join(wd, f))
    return reads


def one(args):
    rid, lam = args
    pg.region_reads = region_reads
    try:
        return pg.genotype(rid, lam, OUT, R + '/work/panelgt/semig')
    except Exception as e:  # noqa: BLE001
        return {'region_id': rid, 'error': repr(e)[:200]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--lam', type=float, default=3.0)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    ids = open(R + '/work/stage4/chr20/patch_all_ids.txt').read().split()
    keys = ['region_id', 'haps', 'reads_kept', 'kappa', 'n_flank', 'called_len', 'truth_len', 'panel_call_ed',
            'ceiling_ed', 'error']
    with open(os.path.join(OUT, 'atomic_chr20_dw%s.tsv' % os.environ.get('PANELGT_DEPTH_W', '1')), 'w') as f, \
            concurrent.futures.ProcessPoolExecutor(a.jobs) as ex:
        f.write('\t'.join(keys) + '\n')
        for x in ex.map(one, [(r, a.lam) for r in ids]):
            f.write('\t'.join(str(x.get(k, '')) for k in keys) + '\n')
            f.flush()
    print('done')


if __name__ == '__main__':
    main()
