#!/usr/bin/env python3
"""For every chr20 region, which of pair_recall.py's distinct hap32 sequences (by index) are the
references (CHM13, GRCh38) and which is the medoid star's centre, so that the pairwise consistency
can be scored without the pairs a star reproduces by construction (a pair with the star's own
centre). -> work/iterate/pair_centres.json  {region: {"ref": [...], "chm13": [...], "medoid": [...]}}"""
import collections, concurrent.futures, csv, json, os, sys
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
sys.path.insert(0, R + '/tools')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import msa_graph, panel, iterate, poa_panel, pair_recall as pr  # noqa: E402


def info(rid):
    rows, names = {}, collections.defaultdict(list)
    for n, r in msa_graph.read_msa('%s/poa_abpoa/%s.msa.fa' % (pr.CAND, rid)):
        if n in msa_graph.CONSENSUS_NAMES:
            continue
        s = r.replace('-', '')
        rows.setdefault(s, r)
        names[s].append(n)
    seqs = sorted(rows, key=len)        # pair_recall.py's order
    ref = [i for i, s in enumerate(seqs) if any(n.startswith(('CHM13', 'GRCh38')) for n in names[s])]
    chm = [i for i, s in enumerate(seqs) if any(n.startswith('CHM13') for n in names[s])]
    fa, mp = poa_panel.union_files(rid)
    recs = msa_graph.read_fasta(fa)
    w = {r['id']: int(r['weight']) for r in panel.read_map(mp)}
    med = iterate.star_medoid([('s%d' % k, s) for k, (n, s) in enumerate(recs)], [w.get(n, 1) for n, _ in recs])
    return rid, {'ref': ref, 'chm13': chm, 'medoid': [i for i, s in enumerate(seqs) if s == med]}


def main():
    ids = sorted({r['region'] for r in csv.DictReader(open(R + '/work/iterate/pair_recall_mc.tsv'), delimiter='\t')})
    with concurrent.futures.ProcessPoolExecutor(6) as ex:
        out = dict(ex.map(info, ids, chunksize=4))
    json.dump(out, open(R + '/work/iterate/pair_centres.json', 'w'))
    print('regions', len(out), 'with a CHM13 row', sum(1 for v in out.values() if v['chm13']),
          'medoid among the hap32 rows', sum(1 for v in out.values() if v['medoid']))


if __name__ == '__main__':
    main()
