#!/usr/bin/env python3
"""Is an alignment's pairwise consistency biased towards alleles that resemble CHM13? Measured on pairs
of full-panel alleles that have nothing to do with HG002.

Per region, pairs of distinct panel alleles are drawn at random (seeded by the region), each allele
with probability proportional to the panel haplotypes that carry it, leaving out CHM13's and GRCh38's
alleles and every star's centre (the medoid and the longest allele). Each pair is aligned with abPOA
(unbanded, default scores) and its aligned position pairs are compared with those of each full-panel
MSA, as pair_recall.py does for the hap32 rows. Each allele's distance from CHM13's allele is 1 minus
the multiset 15-mer Jaccard index (iterate.kmer_tokens).

    python3 work/iterate/pair_bias.py [--pairs 40] [--jobs 6]   -> work/iterate/pair_bias.tsv
"""
import argparse, concurrent.futures, os, random, sys, tempfile, time
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
sys.path.insert(0, R + '/tools')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import msa_graph, panel, poa_panel, iterate  # noqa: E402
from pair_recall import pairs  # noqa: E402

ARMS = {'abpoa_full': 'poa_abpoa__all', 'famsa_full': 'pf_famsa', 'mst3': 'mst3', 'bbt64m': 'bbt64m',
        'st_long': 'st_long', 'st_chm13': 'st_chm13', 'st_maj': 'st_maj', 'st_medoid': 'st_medoid',
        'st_cons': 'st_cons', 'st_maj_med': 'st_maj_med'}


def jaccard_distance(a, b):
    inter = len(a & b)
    uni = len(a) + len(b) - inter
    return 1.0 - (inter / uni if uni else 1.0)


USE = None


def region(args):
    rid, npairs, use = args
    global USE
    USE = use
    fa, mp = poa_panel.union_files(rid)
    rows_map = panel.read_map(mp)
    weight = {r['id']: max(int(r['weight']), 1) for r in rows_map}
    refs = {r['id'] for r in rows_map if any(m.startswith(('CHM13#', 'GRCh38#')) for m in r['members'])}
    chm = next((r['id'] for r in rows_map if any(m.startswith('CHM13#') for m in r['members'])), None)
    msas = {}
    for arm, v in ARMS.items():
        if USE and arm not in USE:
            continue
        f = iterate.full_msa_of(v, rid)
        if f and os.path.exists(f):
            msas[arm] = {n: r for n, r in msa_graph.read_msa(f) if n not in msa_graph.CONSENSUS_NAMES}
    # The alleles and their sequences come from st_chm13's MSA whichever arms are scored, so that a
    # run on a subset of arms draws the same pairs.
    f = iterate.full_msa_of('st_chm13', rid)
    if chm is None or not f or not os.path.exists(f):
        return []
    seq = {n: r.replace('-', '') for n, r in msa_graph.read_msa(f) if n not in msa_graph.CONSENSUS_NAMES}
    # every star's centre is left out: the medoid (as star_medoid picks it) and the longest allele
    ids = sorted(seq)
    med = iterate.star_medoid([('s%d' % k, seq[n]) for k, n in enumerate(ids)], [weight.get(n, 1) for n in ids])
    longest = max(ids, key=lambda n: (len(seq[n]), -int(n[1:])))
    skip = refs | {n for n in ids if seq[n] == med} | {longest}
    pool = [n for n in ids if n not in skip]
    if len(pool) < 2:
        return []
    rng = random.Random(rid)
    w = [weight.get(n, 1) for n in pool]
    chosen = set()
    tries = 0
    while len(chosen) < npairs and tries < npairs * 20:
        tries += 1
        a, b = rng.choices(pool, weights=w, k=2)
        if a != b:
            chosen.add((min(a, b), max(a, b)))
    tok = {n: iterate.kmer_tokens(seq[n]) for n in {x for p in chosen for x in p} | {chm}}
    wd = tempfile.mkdtemp(prefix='bias.', dir=R + '/work/tmp')
    out = []
    for a, b in sorted(chosen):
        sa, sb = seq[a], seq[b]
        try:
            if len(sa) * len(sb) > iterate.PAIR_MAX_CELLS:
                raise iterate.Capped('memout', 'predicted')
            got, _ = iterate._abpoa_rows([('p', sa), ('q', sb)], wd, 'pair', ['-b', '-1'], time.time() + 3600,
                                         iterate.MEM_MB)
        except iterate.Capped:
            continue
        opt = pairs(got['p'], got['q'])
        da, db = jaccard_distance(tok[a], tok[chm]), jaccard_distance(tok[b], tok[chm])
        for arm, rows in msas.items():
            ra, rb = rows.get(a), rows.get(b)
            if ra is None or rb is None or ra.replace('-', '') != sa or rb.replace('-', '') != sb:
                continue
            m = pairs(ra, rb)
            out.append((rid, arm, a, b, len(sa), len(sb), '%.4f' % da, '%.4f' % db, len(opt), len(opt & m), len(m)))
    try:
        os.rmdir(wd)
    except OSError:
        pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pairs', type=int, default=40)
    ap.add_argument('--jobs', type=int, default=6)
    ap.add_argument('--out', default=R + '/work/iterate/pair_bias.tsv')
    ap.add_argument('--arms', help='comma list; default every arm')
    a = ap.parse_args()
    use = set(a.arms.split(',')) if a.arms else None
    ids = open(R + '/work/stage4/chr20/patch_all_ids.txt').read().split()
    with open(a.out, 'w') as o, concurrent.futures.ProcessPoolExecutor(a.jobs) as ex:
        o.write('region\tarm\ta\tb\tlen_a\tlen_b\tdist_a\tdist_b\topt_pairs\tshared\tmsa_pairs\n')
        for lines in ex.map(region, [(r, a.pairs, use) for r in ids]):
            for l in lines:
                o.write('\t'.join(map(str, l)) + '\n')
    print('done')


if __name__ == '__main__':
    main()
