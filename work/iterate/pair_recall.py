#!/usr/bin/env python3
"""Fraction of the aligned position pairs of each optimal pairwise alignment that an MSA also aligns.

For every region and every pair of distinct hap32 sequences, the pair is aligned with abPOA on the two
sequences (unbanded, default scores, as iterate.pair_align but without left normalisation). Its aligned
position pairs (i, j) -- columns where both have a base -- are compared with the pairs each arm's projected
MSA puts in one column. No re-alignment between gaps.

    python3 work/iterate/pair_recall.py [--regions @ids.txt|R1,R2] [--jobs 4]
-> work/iterate/pair_recall.tsv (one line per region, arm and pair) and a per-arm summary.
"""
import argparse, collections, concurrent.futures, gzip, json, os, sys, tempfile, time
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
sys.path.insert(0, R + '/tools')
import msa_graph, iterate  # noqa: E402

ARMS = {'abpoa_h32': 'poa_abpoa', 'abpoa_full': 'poa_abpoa__all', 'famsa_h32': 'pf_famsa_h32__all',
        'famsa_full': 'pf_famsa__all', 'mst': 'mst__all', 'st_long': 'st_long__all', 'st_chm13': 'st_chm13__all',
        'ua_full': 'unit_aware_poa__all', 'bbt64m': 'bbt64m__all', 'mst3': 'mst3__all', 'st_medoid': 'st_medoid__all',
        'st_maj': 'st_maj__all', 'st_cons': 'st_cons__all', 'st_maj_med': 'st_maj_med__all', 'st_amed': 'st_amed__all'}
CAND = R + '/work/stage4/candidates'
CACHE = R + '/work/iterate/pairrecall_cache'   # optimal pairwise alignments, as partner arrays, per region
USE = list(ARMS)


def rows_by_seq(path):
    """{ungapped sequence: its first MSA row}."""
    out = {}
    for n, r in msa_graph.read_msa(path):
        if n not in msa_graph.CONSENSUS_NAMES:
            out.setdefault(r.replace('-', ''), r)
    return out


def pairs(ra, rb):
    """Aligned position pairs (i, j) of two rows of one alignment."""
    out, i, j = set(), 0, 0
    for x, y in zip(ra, rb):
        if x != '-' and y != '-':
            out.add((i, j))
        i += x != '-'
        j += y != '-'
    return out


def region(rid, arms=None):
    arms = arms or USE
    msas = {a: rows_by_seq(f) for a, m in ARMS.items() if a in arms or a == 'abpoa_h32'
            for f in ['%s/%s/%s.msa.fa' % (CAND, m, rid)] if os.path.exists(f)}
    seqs = sorted(msas['abpoa_h32'], key=len)
    cf = '%s/%s.json.gz' % (CACHE, rid)
    cache = json.load(gzip.open(cf, 'rt')) if os.path.exists(cf) else {}
    fresh = False
    wd = tempfile.mkdtemp(prefix='pr.', dir=R + '/work/tmp')
    lines = []
    for x in range(len(seqs)):
        for y in range(x + 1, len(seqs)):
            a, b = seqs[x], seqs[y]
            key = '%d,%d' % (x, y)
            if key in cache:
                band, part = cache[key]
                opt = {(i, j) for i, j in enumerate(part) if j >= 0}
                for arm, rows in msas.items():
                    if arm in arms and a in rows and b in rows:
                        m = pairs(rows[a], rows[b])
                        lines.append((rid, arm, x, y, len(a), len(b), band, len(opt), len(opt & m), len(m)))
                continue
            try:
                if len(a) * len(b) > iterate.PAIR_MAX_CELLS:
                    raise iterate.Capped('memout', 'predicted')
                got, _ = iterate._abpoa_rows([('p', a), ('q', b)], wd, 'pair', ['-b', '-1'], time.time() + 3600,
                                             iterate.MEM_MB)
                band = 0
            except iterate.Capped:
                got, _ = iterate._abpoa_rows([('p', a), ('q', b)], wd, 'pair', [], time.time() + 3600, iterate.MEM_MB)
                band = 1
            opt = pairs(got['p'], got['q'])
            part = [-1] * len(a)
            for i, j in opt:
                part[i] = j
            cache[key] = [band, part]
            fresh = True
            for arm, rows in msas.items():
                if arm in arms and a in rows and b in rows:
                    m = pairs(rows[a], rows[b])
                    lines.append((rid, arm, x, y, len(a), len(b), band, len(opt), len(opt & m), len(m)))
    os.rmdir(wd) if not os.listdir(wd) else None
    if fresh:
        os.makedirs(CACHE, exist_ok=True)
        with gzip.open(cf + '.tmp', 'wt') as f:
            json.dump(cache, f, separators=(',', ':'))
        os.replace(cf + '.tmp', cf)
    return lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--regions', default='@' + R + '/work/stage4/chr20/patch_all_ids.txt')
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--out', default=R + '/work/iterate/pair_recall.tsv')
    ap.add_argument('--arms', default=','.join(ARMS))
    a = ap.parse_args()
    arms = a.arms.split(',')
    ids = open(a.regions[1:]).read().split() if a.regions.startswith('@') else a.regions.split(',')
    ids = [r for r in ids if os.path.exists('%s/poa_abpoa/%s.msa.fa' % (CAND, r))]
    tot = collections.defaultdict(lambda: [0, 0, 0, 0])
    with open(a.out, 'w') as o, concurrent.futures.ProcessPoolExecutor(a.jobs) as ex:
        o.write('region\tarm\ta\tb\tlen_a\tlen_b\tbanded\topt_pairs\tshared\tmsa_pairs\n')
        for lines in ex.map(region, ids, [arms] * len(ids)):
            for l in lines:
                o.write('\t'.join(map(str, l)) + '\n')
                t = tot[l[1]]
                t[0] += 1; t[1] += l[7]; t[2] += l[8]; t[3] += l[9]
    print('arm          pairs   opt pairs contained (recall)   MSA pairs on the optimal alignment (precision)')
    for arm in arms:
        if arm in tot:
            n, op, sh, mp = tot[arm]
            print('%-11s %7d   %.4f   %.4f' % (arm, n, sh / op, sh / mp))


if __name__ == '__main__':
    main()
