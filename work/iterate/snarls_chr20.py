#!/usr/bin/env python3
"""4q: top-level and all snarls of every chr20 candidate graph (work/stage4/chr20/patch_all_ids.txt) for each
method, to compare the region-graph site count with the VCF site count of the whole-chr20 call.

    python3 work/iterate/snarls_chr20.py mc,poa_abpoa,poa_abpoa__all,unit_aware_poa__all,mst__all,mst3__all
-> work/iterate/snarls_chr20.json  {method: {id: [snarls, sites]}} (cached per method and region)
"""
import concurrent.futures
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'tools'))
import iterate  # noqa: E402

OUT = os.path.join(REPO, 'work', 'iterate', 'snarls_chr20.json')


def gfa_of(m, rid):
    if m == 'mc':
        return os.path.join(REPO, 'work', 'stage4', 'regions', rid, 'mc.gfa')
    return os.path.join(REPO, 'work', 'stage4', 'candidates', m, rid + '.gfa')


def one(m, rid):
    g = gfa_of(m, rid)
    if not os.path.exists(g):
        return m, rid, None
    r = iterate.snarl_counts(g)
    return m, rid, None if r.get('status') == 'error' else [r['snarls'], r['sites']]


def main():
    ids = open(os.path.join(REPO, 'work', 'stage4', 'chr20', 'patch_all_ids.txt')).read().split()
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    todo = [(m, r) for m in sys.argv[1].split(',') for r in ids if not res.get(m, {}).get(r)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(os.environ.get('JOBS', '3'))) as ex:
        for m, rid, c in ex.map(lambda a: one(*a), todo):
            res.setdefault(m, {})[rid] = c
    with open(OUT, 'w') as f:
        json.dump(res, f)
    for m in sys.argv[1].split(','):
        v = [x for x in res[m].values() if x]
        print('%-22s regions %3d  snarls %6d  sites %6d' % (m, len(v), sum(x[0] for x in v), sum(x[1] for x in v)))


if __name__ == '__main__':
    main()
