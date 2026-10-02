#!/usr/bin/env python3
"""Data for the alignment viewer: per region, the distinct hap32 sequences, each arm's projected MSA rows
for them, the optimal pairwise alignment of every pair (as partner arrays: for each base of sequence i,
the base of sequence j it is aligned to, or -1), and the full-panel MSA of two arms with the hap32 rows
marked. -> one gzip-compressed JSON (base64) per run, embedded by the viewer page.

    python3 work/iterate/viz_alignments.py OUT.b64 RID[,RID...] [--full RID,...]
"""
import base64, gzip, json, os, sys, tempfile, time
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
sys.path.insert(0, R + '/tools')
import msa_graph, iterate  # noqa: E402

ARMS = [('abpoa_h32', 'poa_abpoa', 'abPOA on the hap32 rows (sample-dependent)'),
        ('abpoa_full', 'poa_abpoa__all', 'abPOA on the full panel, projected'),
        ('st_chm13', 'st_chm13__all', 'centre-star to CHM13 on the full panel, projected'),
        ('famsa_full', 'pf_famsa__all', 'FAMSA on the full panel, projected'),
        ('st_maj', 'st_maj__all', 'centre-star to the majority consensus of the full panel, projected'),
        ('st_medoid', 'st_medoid__all', 'centre-star to the medoid allele of the full panel, projected'),
        ('st_long', 'st_long__all', 'centre-star to the longest allele of the full panel, projected'),
        ('mst3', 'mst3__all', 'nearest-neighbour threading on the full panel, projected')]
FULL = [('abpoa_full', 'poa_abpoa__all'), ('st_chm13', 'st_chm13')]
CAND = R + '/work/stage4/candidates'


def rows_by_seq(path):
    out, names = {}, {}
    for n, r in msa_graph.read_msa(path):
        if n not in msa_graph.CONSENSUS_NAMES:
            s = r.replace('-', '')
            out.setdefault(s, r)
            names.setdefault(s, []).append(n)
    return out, names


def partner(ra, rb):
    """For each base of row ra, the index of the base of rb in its column, or -1."""
    out, j = [], 0
    for x, y in zip(ra, rb):
        if x != '-':
            out.append(j if y != '-' else -1)
        j += y != '-'
    return out


def region(rid, with_full):
    base, names = rows_by_seq('%s/poa_abpoa/%s.msa.fa' % (CAND, rid))
    seqs = sorted(base, key=lambda s: (len(s), s))
    d = {'rid': rid, 'names': [names[s] for s in seqs], 'lens': [len(s) for s in seqs], 'arms': {}, 'desc': {}}
    for arm, m, desc in ARMS:
        f = '%s/%s/%s.msa.fa' % (CAND, m, rid)
        if os.path.exists(f):
            rows, _ = rows_by_seq(f)
            d['arms'][arm] = [rows[s] for s in seqs]
            d['desc'][arm] = desc
    wd = tempfile.mkdtemp(prefix='viz.', dir=R + '/work/tmp')
    opt = {}
    # pair_recall.py's cache, in the same order of sequences (by length, then MSA order)
    cf = '%s/work/iterate/pairrecall_cache/%s.json.gz' % (R, rid)
    cache = json.load(gzip.open(cf, 'rt')) if os.path.exists(cf) else {}
    order = sorted(base, key=len)
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            a, b = order.index(seqs[i]), order.index(seqs[j])
            hit = cache.get('%d,%d' % (min(a, b), max(a, b)))
            if hit is not None:
                part = hit[1]
                if a > b:   # cached for (b, a): invert
                    inv = [-1] * len(seqs[i])
                    for k, v in enumerate(part):
                        if v >= 0:
                            inv[v] = k
                    part = inv
                opt['%d,%d' % (i, j)] = part
                continue
            got, _ = iterate._abpoa_rows([('p', seqs[i]), ('q', seqs[j])], wd, 'pair', ['-b', '-1'],
                                         time.time() + 3600, iterate.MEM_MB)
            opt['%d,%d' % (i, j)] = partner(got['p'], got['q'])
    d['opt'] = opt
    if with_full:
        idx = {s: k for k, s in enumerate(seqs)}
        d['full'] = {}
        for arm, v in FULL:
            rows = [r for n, r in msa_graph.read_msa(iterate.full_msa_of(v, rid)) if n not in msa_graph.CONSENSUS_NAMES]
            hap = [idx.get(r.replace('-', ''), -1) for r in rows]
            order = sorted(range(len(rows)), key=lambda k: len(rows[k].replace('-', '')))
            d['full'][arm] = {'rows': [rows[k] for k in order], 'hap': [hap[k] for k in order]}
    return d


def main():
    out, rids = sys.argv[1], sys.argv[2].split(',')
    full = set(sys.argv[sys.argv.index('--full') + 1].split(',')) if '--full' in sys.argv else set()
    data = [region(r, r in full) for r in rids]
    raw = json.dumps(data, separators=(',', ':')).encode()
    b = base64.b64encode(gzip.compress(raw, 9)).decode()
    open(out, 'w').write(b)
    print('regions', len(data), 'json %.1f MB, embedded %.1f MB' % (len(raw) / 1e6, len(b) / 1e6))


if __name__ == '__main__':
    main()
