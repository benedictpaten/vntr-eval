import csv, collections, concurrent.futures, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frag import runs  # noqa
from frag import one
ARMS = {'st_long': 'st_long__all', 'h32f': 'pf_famsa_h32__all', 'all': 'poa_abpoa', 'full': 'poa_abpoa__all', 'ua': 'unit_aware_poa__all', 'mst': 'mst__all',
        'mst3': 'mst3__all', 'bbt64m': 'bbt64m__all', 'pf_famsa': 'pf_famsa__all'}
def gr(a):
    arm, m, rid = a
    import frag
    f = '/Users/benedictpaten/PycharmProjects/vntr-eval/work/stage4/candidates/%s/%s.msa.fa' % (m, rid)
    if not os.path.exists(f): f = '/Users/benedictpaten/PycharmProjects/vntr-eval/work/iterate/candidates/%s/%s.msa.fa' % (m, rid)
    if not os.path.exists(f): return arm, rid, None
    rows = {}
    for n, r in frag.msa_graph.read_msa(f):
        if n not in frag.msa_graph.CONSENSUS_NAMES: rows.setdefault(r.replace('-', ''), r)
    t, d = 0, set()
    for r in rows.values():
        for x, y in runs(r):
            if x > 0 and y < len(r): t += 1; d.add((x, y))
    return arm, rid, (t, len(d))
if __name__ == '__main__':
    rows = list(csv.DictReader(open('region_metrics.tsv'), delimiter='\t'))
    rids = sorted({r['region'] for r in rows})
    with concurrent.futures.ProcessPoolExecutor(6) as ex:
        res = list(ex.map(gr, [(a, m, r) for a, m in ARMS.items() for r in rids], chunksize=8))
    with open('gapruns_chr20.tsv', 'w') as o:
        o.write('arm\tregion\truns\tdistinct\n')
        for a, r, v in res:
            if v: o.write('%s\t%s\t%d\t%d\n' % (a, r, v[0], v[1]))
    print('done')
