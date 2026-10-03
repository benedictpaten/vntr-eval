#!/usr/bin/env python3
"""Covered chr20 repeat regions: does vg pick the wrong panel walk because of the reads it has, or because of
how it scores them?

vg scores a read against an allele by pairing the node visits the mapper gave it with the allele's visits, so a
read the mapper put on nodes an allele does not visit fits that allele badly, whatever its bases. The
whole-allele prototype (tools/panel_genotype.py) aligns each read's bases to each allele instead. Here the
prototype's model is applied to vg's own reads: the alignments in the medoid-star GAF that visit a node inside
the region, their sequence rebuilt from the path and the cs difference string.

Per region it writes the called pair under that model, and for every read:
  - the alleles it fits best by sequence (lowest semi-global edit distance),
  - the alleles whose node walk contains every node the read visits inside the region (its placement).

    python3 work/iterate/wrongwalk/rescore.py [--jobs 6]
        -> work/iterate/wrongwalk/rescore.tsv, rescore_reads.tsv, and a summary on stdout
"""
import argparse, collections, concurrent.futures, csv, gzip, json, math, os, re, subprocess, sys
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
P = R + '/work/stage4/chr20/patched_st_medoid'
W = R + '/work/iterate/wrongwalk'
SEMIG = R + '/work/panelgt/semig'
sys.path.insert(0, R + '/tools')
import msa_graph                                   # noqa: E402
from score_haplotypes import edit_distance as ED   # noqa: E402

FL, RL, LAM, MAX_ERR = 2000, 150.0, 3.0, 0.08     # the prototype's constants
COMP = str.maketrans('ACGTN', 'TGCAN')
TOK = re.compile(r'([<>])(\d+)')


def rc(s):
    return s.translate(COMP)[::-1]


def load_regions():
    rows = {r['region']: r for r in csv.DictReader(open(W + '/classify.tsv'), delimiter='\t')}
    span = {}
    for rid in rows:
        j = json.load(open(R + '/work/stage4/regions/%s/region.json' % rid))
        span[rid] = (int(j['span_start']), int(j['span_end']))   # 1-based inclusive
    return rows, span


def graph_pass(span):
    """CHM13 node positions, every haplotype's walk between each region's two anchor nodes (the CHM13 nodes
    just outside the span), and the CHM13 flank nodes (2 kb each side) used for the read-start rate."""
    seglen, walks = {}, []
    with open(P + '/graph.gfa') as f:
        for l in f:
            if l[0] == 'S':
                a = l.split('\t', 3)
                seglen[int(a[1])] = len(a[2].rstrip('\n'))
            elif l[0] == 'W':
                a = l.rstrip('\n').split('\t')
                walks.append(('%s#%s#%s#%s' % (a[1], a[2], a[3], a[4]), a[6]))
    chm = next(w for n, w in walks if n.startswith('CHM13#'))
    toks = [int(x) for _, x in TOK.findall(chm)]
    pos, p = {}, 0
    starts = []
    for n in toks:
        starts.append((p, n))
        pos.setdefault(n, p)
        p += seglen[n]
    import bisect
    sp = [s for s, _ in starts]

    def node_at(x):          # 0-based position -> node
        return starts[bisect.bisect_right(sp, x) - 1][1]
    anchors, flank = {}, {}
    for rid, (a, b) in span.items():
        anchors[rid] = (node_at(a - 2), node_at(b))
        for lo, hi in ((max(0, a - 1 - FL), a - 1), (b, b + FL)):
            k = bisect.bisect_right(sp, lo) - 1
            while k < len(starts) and starts[k][0] < hi:
                flank.setdefault(starts[k][1], set()).add(rid)
                k += 1
    for rid, (l, r) in anchors.items():       # anchors and interior are not flank
        flank.get(l, set()).discard(rid)
        flank.get(r, set()).discard(rid)
    want = collections.defaultdict(list)
    for rid, (l, r) in anchors.items():
        want[l].append((rid, 'L'))
        want[r].append((rid, 'R'))
    allele_walks = collections.defaultdict(list)       # rid -> [(walk name, [(orient, node), ...])]
    maxlen, dropped = {}, collections.Counter()
    for rid in span:
        hl = max(len(s) for _, s in msa_graph.read_fasta(R + '/work/stage4/regions/%s/hap32.fa' % rid))
        l, r = anchors[rid]
        maxlen[rid] = 2 * hl + seglen[l] + seglen[r] + 1000
    for name, w in walks:
        t = TOK.findall(w)
        hits = collections.defaultdict(dict)
        for i, (o, x) in enumerate(t):
            n = int(x)
            if n in want:
                for rid, side in want[n]:
                    hits[rid].setdefault(side, []).append(i)
        for rid, h in hits.items():
            if 'L' not in h or 'R' not in h or len(h['L']) != 1 or len(h['R']) != 1:
                continue
            i, j = h['L'][0], h['R'][0]
            bp = sum(seglen[int(x)] for _, x in t[min(i, j):max(i, j) + 1])
            if bp > maxlen[rid]:               # a walk that wanders off between the anchors: not an allele
                dropped[rid] += 1
                continue
            if i < j:
                allele_walks[rid].append((name, [(o, int(x)) for o, x in t[i:j + 1]]))
            else:                              # the walk crosses the region backwards: flip it
                seg = t[j:i + 1][::-1]
                allele_walks[rid].append((name, [('<' if o == '>' else '>', int(x)) for o, x in seg]))
    del walks
    print('walks dropped as too long: %d, in %d regions' % (sum(dropped.values()), len(dropped)), flush=True)
    return seglen, pos, anchors, flank, allele_walks


def read_pass(node_region, flank, pos):
    reads = collections.defaultdict(list)
    nflank = collections.Counter()
    with gzip.open(P + '/reads.gaf.gz', 'rt') as f:
        for l in f:
            if l[0] == '@':
                continue
            a = l.split('\t')
            nodes = [int(x) for _, x in TOK.findall(a[5])]
            rs = {node_region[n] for n in nodes if n in node_region}
            if rs:
                cs = next((x[5:] for x in a[12:] if x.startswith('cs:Z:')), None)
                for rid in rs:
                    reads[rid].append((a[0], int(a[11]), a[4], a[5], int(a[7]), int(a[8]), cs.rstrip('\n') if cs else None))
            else:
                # a read of the read-start rate: every reference node it visits is in one region's flanks
                # (non-reference nodes, a SNP's ALT say, are allowed)
                fl = None
                for n in nodes:
                    if n not in pos:
                        continue
                    s = flank.get(n)
                    if not s:
                        fl = None
                        break
                    fl = s if fl is None else fl & s
                    if not fl:
                        break
                for rid in fl or ():
                    nflank[rid] += 1
    return reads, nflank


CSOP = re.compile(r'(:\d+|\*[a-zA-Z][a-zA-Z]|\+[a-zA-Z]+|-[a-zA-Z]+)')


def rebuild(path, pstart, pend, cs, seq):
    """The read's aligned bases in path orientation, from the path's sequence and the cs string."""
    ps = ''.join(seq[int(x)] if o == '>' else rc(seq[int(x)]) for o, x in TOK.findall(path))[pstart:pend]
    out, i = [], 0
    for op in CSOP.findall(cs):
        if op[0] == ':':
            k = int(op[1:]); out.append(ps[i:i + k]); i += k
        elif op[0] == '*':
            out.append(op[2].upper()); i += 1
        elif op[0] == '+':
            out.append(op[1:].upper())
        else:
            i += len(op) - 1
    return ''.join(out)


def genotype(D, L, kappa, depth_w=1.0):
    n = len(L)
    N = len(D)
    best, bp = -math.inf, (0, 0)
    for i in range(n):
        for j in range(i, n):
            ll = 0.0
            for d in D:
                x, y = -LAM * d[i], -LAM * d[j]
                m = max(x, y)
                ll += m + math.log(0.5 * math.exp(x - m) + 0.5 * math.exp(y - m))
            if depth_w > 0 and kappa > 0:
                mu = kappa * (L[i] + L[j] + 2 * (RL - 1))
                ll += depth_w * (N * math.log(mu) - mu - math.lgamma(N + 1))
            if ll > best:
                best, bp = ll, (i, j)
    return bp


def one(job):
    rid, row, rreads, kappa, walks, seq = job
    d = R + '/work/stage4/regions/' + rid
    truth = [s.upper() for _, s in msa_graph.read_fasta(d + '/truth.fa')]
    hap = [(n, s.upper()) for n, s in msa_graph.read_fasta(d + '/hap32.fa')]
    alle = sorted({s for _, s in hap})
    idx = {s: i for i, s in enumerate(alle)}
    # each walk's span sequence: spell anchor..anchor and cut the anchors' outside parts by the CHM13 offsets
    (lo, ro), wa = row['_cut'], collections.defaultdict(set)
    nwalk_match = 0
    for name, w in walks:
        s = ''.join(seq[n] if o == '>' else rc(seq[n]) for o, n in w)
        s = s[lo:len(s) - ro] if ro else s[lo:]
        if s in idx:
            nwalk_match += 1
            wa[idx[s]].add(frozenset(n for _, n in w[1:-1]))
    seqs, places = [], []
    for name, mq, strand, path, ps, pe, cs in rreads:
        if not cs:
            continue
        s = rebuild(path, ps, pe, cs, seq)
        inside = frozenset(int(x) for _, x in TOK.findall(path)) & row['_interior']
        seqs.append(s)
        places.append((name, mq, inside))
    wd = W + '/rescore/' + rid
    os.makedirs(wd, exist_ok=True)
    with open(wd + '/haps.fa', 'w') as f:
        for i, h in enumerate(alle):
            f.write('>h%d\n%s\n' % (i, h))
    with open(wd + '/reads.fa', 'w') as f:
        for i, s in enumerate(seqs):
            f.write('>r%d\n%s\n' % (i, s))
    out = subprocess.run([SEMIG, wd + '/haps.fa', wd + '/reads.fa'], capture_output=True, text=True, check=True).stdout
    D, keep = [], []
    for line, s, pl in zip(out.splitlines(), seqs, places):
        dd = [int(x) for x in line.split('\t')[1:]]
        if min(dd) <= MAX_ERR * len(s):
            D.append(dd); keep.append(pl)
    L = [len(h) for h in alle]
    res = {'region': rid}
    for lab, w in (('vgreads', 1.0), ('vgreads_nodepth', 0.0)):
        i, j = genotype(D, L, kappa, w) if D else (0, 0)
        res[lab] = min(ED(alle[i], truth[0]) + ED(alle[j], truth[1]), ED(alle[i], truth[1]) + ED(alle[j], truth[0]))
        res[lab + '_len'] = '%d/%d' % (L[i], L[j])
    # the ceiling pair and vg's called pair, as allele indices
    dt = [[ED(a, t) for t in truth] for a in alle]
    ceil = min((dt[i][0] + dt[j][1], i, j) for i in range(len(alle)) for j in range(len(alle)))
    called = [s.upper() for _, s in msa_graph.read_fasta(R + '/work/iterate/hapscore/rep_site/%s.fa' % rid)]
    ci = [idx.get(c, -1) for c in called]
    right = {ceil[1], ceil[2]}
    wrong_called = {c for c in ci if c >= 0 and c not in right}
    # per read: best-by-sequence alleles, and whether its placement fits a right / a wrongly called allele
    rows_out = []
    for dd, (name, mq, inside) in zip(D, keep):
        m = min(dd)
        seq_best = {k for k, x in enumerate(dd) if x == m}
        fits = {k for k, ws in wa.items() if any(inside <= w for w in ws)} if inside else set(wa)
        rows_out.append((rid, name, mq, len(inside), int(bool(seq_best & right)), int(bool(seq_best & wrong_called)),
                         int(bool(fits & right)), int(bool(fits & wrong_called)), min(dd[k] for k in right) - m,
                         min((dd[k] for k in wrong_called), default=-1) - m if wrong_called else -1))
    res.update({'reads': len(D), 'kappa': round(kappa, 4), 'alleles': len(alle), 'walks': len(walks),
                'walks_matching_hap32': nwalk_match, 'alleles_with_walk': len(wa), 'wrong_called': len(wrong_called)})
    return res, rows_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=6)
    a = ap.parse_args()
    rows, span = load_regions()
    print('graph pass', flush=True)
    seglen, pos, anchors, flank, aw = graph_pass(span)
    node_region, interior = {}, {}
    for rid, ws in aw.items():
        l, r = anchors[rid]
        nodes = {n for _, w in ws for _, n in w} - {l, r}
        interior[rid] = frozenset(nodes)
        for n in nodes:
            node_region[n] = rid
    print('regions with walks %d of %d; interior nodes %d' % (len(aw), len(rows), len(node_region)), flush=True)
    print('read pass', flush=True)
    reads, nflank = read_pass(node_region, flank, pos)
    print('reads kept %d' % sum(len(v) for v in reads.values()), flush=True)
    need = set(node_region) | {n for l, r in anchors.values() for n in (l, r)}
    for v in reads.values():
        for r in v:
            need.update(int(x) for _, x in TOK.findall(r[3]))
    seq = {}
    with open(P + '/graph.gfa') as f:
        for l in f:
            if l[0] == 'S':
                x = l.split('\t', 3)
                n = int(x[1])
                if n in need:
                    seq[n] = x[2].rstrip('\n').upper()
    jobs = []
    for rid, row in rows.items():
        if rid not in aw:
            continue
        a1, b1 = span[rid]
        l, r = anchors[rid]
        row['_cut'] = ((a1 - 1) - pos[l], (pos[r] + seglen[r]) - b1)
        row['_interior'] = interior[rid]
        kappa = nflank[rid] / (4.0 * (FL - RL + 1))
        rn = {n for _, w in aw[rid] for _, n in w}
        for r in reads.get(rid, []):
            rn.update(int(x) for _, x in TOK.findall(r[3]))
        sub = {n: seq[n] for n in rn if n in seq}
        jobs.append((rid, row, reads.get(rid, []), kappa, aw[rid], sub))
    print('scoring %d regions' % len(jobs), flush=True)
    res_all, reads_all = [], []
    with concurrent.futures.ProcessPoolExecutor(a.jobs) as ex:
        for res, ro in ex.map(one, jobs, chunksize=4):
            res_all.append(res); reads_all += ro
    keys = list(res_all[0])
    with open(W + '/rescore.tsv', 'w') as f:
        f.write('\t'.join(keys) + '\n')
        for r in res_all:
            f.write('\t'.join(str(r[k]) for k in keys) + '\n')
    with open(W + '/rescore_reads.tsv', 'w') as f:
        f.write('region\tread\tmapq\tinside_nodes\tseq_fits_right\tseq_fits_wrong\tplaced_fits_right\tplaced_fits_wrong\tright_minus_best\twrong_minus_best\n')
        for r in reads_all:
            f.write('\t'.join(str(x) for x in r) + '\n')
    print('done')


if __name__ == '__main__':
    main()
