#!/usr/bin/env python3
"""Optimal pairwise pairs contained in the MC graph (the unpatched arm), as pair_recall.py measures them
for the MSAs. Two bases are aligned in the graph when they occupy the same position (segment, offset
in the segment's forward orientation); a path that visits a segment twice gives a base more than one
partner. Uses the optimal alignments cached by pair_recall.py, in its order of sequences.

    python3 work/iterate/pair_recall_mc.py [--jobs 6]   -> work/iterate/pair_recall_mc.tsv
"""
import argparse, collections, concurrent.futures, gzip, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pair_recall as pr  # noqa: E402

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def gfa_paths(path):
    seg, paths = {}, {}
    for line in open(path):
        if line[0] == 'S':
            t = line.rstrip('\n').split('\t')
            seg[t[1]] = t[2]
        elif line[0] == 'P':
            t = line.rstrip('\n').split('\t')
            paths[t[1]] = [(s[:-1], s[-1]) for s in t[2].split(',')]
    out = {}
    for name, steps in paths.items():
        seq, pos = [], []
        for sid, o in steps:
            s = seg[sid]
            if o == '+':
                seq.append(s)
                pos.extend((sid, k) for k in range(len(s)))
            else:
                seq.append(s[::-1].translate(COMP))
                pos.extend((sid, len(s) - 1 - k) for k in range(len(s)))
        out[name] = (''.join(seq), pos)
    return out


def region(rid):
    base, _ = pr.rows_by_seq('%s/poa_abpoa/%s.msa.fa' % (pr.CAND, rid)), None
    seqs = sorted(base, key=len)
    cf = '%s/%s.json.gz' % (pr.CACHE, rid)
    if not os.path.exists(cf):
        return [], {'no_cache': 1}
    cache = json.load(gzip.open(cf, 'rt'))
    g = gfa_paths('%s/work/stage4/regions/%s/mc.gfa' % (pr.R, rid))
    by_seq = {}
    for name, (s, pos) in g.items():
        by_seq.setdefault(s, pos)
    stats = collections.Counter(seqs=len(seqs), seqs_in_graph=sum(1 for s in seqs if s in by_seq))
    lines = []
    for x in range(len(seqs)):
        for y in range(x + 1, len(seqs)):
            a, b = seqs[x], seqs[y]
            key = '%d,%d' % (x, y)
            if a not in by_seq or b not in by_seq or key not in cache:
                stats['pairs_skipped'] += 1
                continue
            band, part = cache[key]
            opt = {(i, j) for i, j in enumerate(part) if j >= 0}
            where = collections.defaultdict(list)
            for j, p in enumerate(by_seq[b]):
                where[p].append(j)
            m = {(i, j) for i, p in enumerate(by_seq[a]) for j in where.get(p, ())}
            lines.append((rid, 'mc', x, y, len(a), len(b), band, len(opt), len(opt & m), len(m)))
    return lines, stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=6)
    ap.add_argument('--out', default=pr.R + '/work/iterate/pair_recall_mc.tsv')
    a = ap.parse_args()
    ids = open(pr.R + '/work/stage4/chr20/patch_all_ids.txt').read().split()
    ids = [r for r in ids if os.path.exists('%s/poa_abpoa/%s.msa.fa' % (pr.CAND, r))]
    tot = collections.Counter()
    with open(a.out, 'w') as o, concurrent.futures.ProcessPoolExecutor(a.jobs) as ex:
        o.write('region\tarm\ta\tb\tlen_a\tlen_b\tbanded\topt_pairs\tshared\tmsa_pairs\n')
        for lines, stats in ex.map(region, ids):
            tot.update(stats)
            for l in lines:
                o.write('\t'.join(map(str, l)) + '\n')
    print(dict(tot))


if __name__ == '__main__':
    main()
