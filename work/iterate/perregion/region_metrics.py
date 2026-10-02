#!/usr/bin/env python3
"""Per-region refined SV FP/FN (inside patched_all spans +-100 bp) for every whole-chr20 arm, beside the
region graph's kmerx (extra k-mer positions of the projected MSA), node count and top-level site count.
-> region_metrics.tsv; prints arm totals and per-region correlations of the difference from abPOA hap32."""
import collections, concurrent.futures, gzip, json, os, sys
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
sys.path.insert(0, R + '/tools')
import msa_graph, gap_norm
D = R + '/work/stage4/chr20'
ARMS = {'unpatched': None, 'all': 'poa_abpoa', 'full': 'poa_abpoa__all', 'ua': 'unit_aware_poa__all',
        'mst': 'mst__all', 'mst3': 'mst3__all', 'bbt64m': 'bbt64m__all', 'pf_famsa': 'pf_famsa__all', 'pf_famsa_h32': 'pf_famsa_h32__all', 'st_long': 'st_long__all', 'st_chm13': 'st_chm13__all'}
patch = json.load(open(D + '/patched_all/patch.json'))['regions']
spans = {r: v['span'] for r, v in patch.items() if isinstance(v, dict) and 'span' in v}
iv = sorted((s[0] - 100, s[1] + 100, r) for r, s in spans.items())

def count(path, want):
    c = collections.Counter()
    for line in gzip.open(path, 'rt'):
        if line[0] == '#':
            continue
        x = line.rstrip('\n').split('\t')
        if x[9].split(':')[x[8].split(':').index('BD')] != want:
            continue
        pos = int(x[1])
        for a, b, r in iv:
            if a <= pos <= b:
                c[r] += 1
                break
    return c

def kmx(args):
    arm, m, rid = args
    f = '%s/work/stage4/candidates/%s/%s.msa.fa' % (R, m, rid)
    if not os.path.exists(f):
        return arm, rid, None
    rows = [r for n, r in msa_graph.read_msa(f) if n not in msa_graph.CONSENSUS_NAMES]
    e = gap_norm.kmer_excess(rows)
    g = '%s/work/stage4/candidates/%s/%s.gfa' % (R, m, rid)
    nodes = sum(1 for l in open(g) if l[0] == 'S') if os.path.exists(g) else None
    return arm, rid, (e['positions'], e['extra'], nodes)

def main():
    sv = {}
    for arm in ARMS:
        t = '%s/score_rb/rb_%s/score/chr20.truvari/' % (D, 'unpatched' if arm == 'unpatched' else 'patched_' + arm)
        sv[arm] = (count(t + 'refine.comp.vcf.gz', 'FP'), count(t + 'refine.base.vcf.gz', 'FN'))
    sn = json.load(open(R + '/work/iterate/snarls_chr20.json'))
    jobs = [(a, m, r) for a, m in ARMS.items() if m for r in spans]
    km = {}
    with concurrent.futures.ProcessPoolExecutor(6) as ex:
        for a, r, v in ex.map(kmx, jobs, chunksize=8):
            km[(a, r)] = v
    out = open(os.path.dirname(os.path.abspath(__file__)) + '/region_metrics.tsv', 'w')
    out.write('region\tarm\tfp\tfn\tpositions\textra\tnodes\tsites\n')
    for r in sorted(spans):
        for a, m in ARMS.items():
            k = km.get((a, r)) or (None, None, None)
            s = (sn.get(m or 'mc', {}).get(r) or [None, None])[1]
            out.write('\t'.join(str(x) for x in (r, a, sv[a][0][r], sv[a][1][r], k[0], k[1], k[2], s)) + '\n')
    out.close()
    print('wrote region_metrics.tsv')


if __name__ == '__main__':
    main()
