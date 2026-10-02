#!/usr/bin/env python3
"""The medoid of a region's alleles estimated from an alignment instead of from k-mers.

For an MSA, the weighted sum of one row's induced distances to every other row, where the distance
of two rows is the number of columns in which they differ (a mismatch, or a base against a gap;
gap against gap is not counted), is sum over columns of (W - w_c(x)): x is the row's character in
that column, w_c(x) the panel weight of the rows that have x there, W the total weight. So the
medoid is one pass over the MSA. Weights are the union map's `weight` (panel haplotypes carrying
the allele; 0 for a sequence only the sampled haplotypes carry).

    python3 work/iterate/aln_medoid.py SOURCE_VARIANT   -> compares with the k-mer medoid per region
"""
import collections, concurrent.futures, json, os, sys
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
sys.path.insert(0, R + '/tools')
import msa_graph, panel, poa_panel, iterate  # noqa: E402


def aln_medoid(rows, weight):
    """rows: {name: aligned row}; weight: {name: w}. Returns (name, scores)."""
    names = sorted(rows)
    width = len(rows[names[0]])
    total = float(sum(weight.get(n, 0) for n in names))
    score = dict.fromkeys(names, 0.0)
    for c in range(width):
        w = collections.Counter()
        for n in names:
            w[rows[n][c]] += weight.get(n, 0)
        for n in names:
            score[n] += total - w[rows[n][c]]
    best = min(names, key=lambda n: (score[n], n))
    return best, score


def region(args):
    rid, source = args
    f = iterate.full_msa_of(source, rid)
    if not f or not os.path.exists(f):
        return rid, None
    rows = {n: r for n, r in msa_graph.read_msa(f) if n not in msa_graph.CONSENSUS_NAMES}
    fa, mp = poa_panel.union_files(rid)
    w = {r['id']: int(r['weight']) for r in panel.read_map(mp)}
    best, score = aln_medoid(rows, w)
    # the k-mer medoid, as star_medoid picks it
    recs = msa_graph.read_fasta(fa)
    km = iterate.star_medoid([('s%d' % k, s) for k, (n, s) in enumerate(recs)], [w.get(n, 0) for n, _ in recs])
    kmid = next(n for n, s in recs if s == km)
    seq = {n: r.replace('-', '') for n, r in rows.items()}
    return rid, {'aln': best, 'kmer': kmid, 'same': best == kmid or seq.get(best) == seq.get(kmid),
                 'len_aln': len(seq[best]), 'len_kmer': len(seq.get(kmid, '')), 'n': len(rows),
                 'w_aln': w.get(best, 0), 'w_kmer': w.get(kmid, 0),
                 'rel_gap': (score.get(kmid, 0) - score[best]) / score[best] if score[best] else 0.0}


def main():
    source = sys.argv[1] if len(sys.argv) > 1 else 'st_medoid'
    ids = open(R + '/work/stage4/chr20/patch_all_ids.txt').read().split()
    with concurrent.futures.ProcessPoolExecutor(6) as ex:
        res = {rid: v for rid, v in ex.map(region, [(r, source) for r in ids], chunksize=4) if v}
    json.dump(res, open(R + '/work/iterate/aln_medoid_%s.json' % source, 'w'))
    diff = [v for v in res.values() if not v['same']]
    print('regions', len(res), 'alignment medoid differs from the k-mer medoid in', len(diff))
    if diff:
        import statistics
        print('  of those: median length |aln - kmer| %d bp; k-mer medoid scores %.1f%% worse (median) by the alignment distance'
              % (statistics.median(abs(v['len_aln'] - v['len_kmer']) for v in diff), 100 * statistics.median(v['rel_gap'] for v in diff)))
        print('  panel weight of the chosen allele, median: alignment %d, k-mer %d' %
              (statistics.median(v['w_aln'] for v in diff), statistics.median(v['w_kmer'] for v in diff)))


if __name__ == '__main__':
    main()
