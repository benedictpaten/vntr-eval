#!/usr/bin/env python3
"""Per realigned region (test/patch.json) and per called strand of fix4.star2.{short,ont}: is the strand's
sequence across the region spelled by a panel walk of star2/sampled.gbz (panel_star2.tsv.gz), and how far
is it from HG002's truth haplotypes compared with the best panel walk?

Haplotypes are built with vntr-eval tools/score_haplotypes.py's rule (Prepared/build), over the span
s-1..e+1 (the interior plus the last base of L and the first base of R, so an insertion at the interior's
edge is inside), in the VCF's own phase (one PS for the contig; strand 0 = GT slot 1). Truth: the
T2T-Q100 v1.1 dipcall VCF (truth.chr20.stvar.vcf.gz, all sizes, phased), built the same way.
Output called.tsv: one row per region x arm x strand."""
import bisect, collections, gzip, json, sys
sys.path.insert(0, '/Users/benedictpaten/PycharmProjects/vntr-eval/tools')
import score_haplotypes as S

C = '/Users/benedictpaten/PycharmProjects/vntr-eval/work/full/chr20'
D = C + '/investigate/nested_dup'
W = '/Users/benedictpaten/PycharmProjects/vg-call-eval/work/wgs/chr20'
OUT = '/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-2026-07-vg-call-refactor/e061d8cf-7996-49a4-986b-3686374ef250/scratchpad/linkage'
REF = ''.join(l.strip() for l in open(W + '/chr20.fa') if l[0] != '>').upper()
A = {r: tuple(v['span']) for r, v in json.load(open(C + '/test/patch.json'))['regions'].items()}
regs = sorted(A, key=lambda r: A[r])
starts = [A[r][0] - 1 for r in regs]
ends = [A[r][1] + 1 for r in regs]
bed = [tuple(map(int, l.split()[1:3])) for l in open(W + '/truth.chr20.stvar.bed')]
bed.sort()
bs = [b[0] for b in bed]


def in_bed(a1, b1):          # 1-based inclusive a1..b1 inside one 0-based half-open interval
    i = bisect.bisect_right(bs, a1 - 1) - 1
    return i >= 0 and bed[i][0] <= a1 - 1 and b1 <= bed[i][1]


def load(vcf):
    by = collections.defaultdict(list)
    n = 0
    for line in gzip.open(vcf, 'rt'):
        if line[0] == '#':
            continue
        r = S.parse_record(line, 0, n)
        n += 1
        lo = bisect.bisect_left(ends, r['pos'])
        hi = bisect.bisect_right(starts, r['end'])
        for k in range(lo, hi):
            by[regs[k]].append(r)
    return by


panel = collections.defaultdict(lambda: collections.defaultdict(list))
for l in gzip.open('/Users/benedictpaten/PycharmProjects/vntr-eval/work/full/chr20/repeat_sites/panel_banded.tsv.gz', 'rt'):
    r, name, st, o, s = l.rstrip('\n').split('\t')
    a, b = A[r]
    panel[r][REF[a - 2] + s.upper() + REF[b]].append(name)
print('panel loaded', file=sys.stderr, flush=True)

truth = load(W + '/truth.chr20.stvar.vcf.gz')
arms = {a.split('=')[0]: load(a.split('=')[1]) for a in sys.argv[1:]}
print('vcfs loaded', file=sys.stderr, flush=True)


def build(recs, a1, b1, refseq):
    p = S.Prepared(recs, a1, b1, refseq, 2)
    seqs, skipped, applied = p.build([False] * p.n_blocks)
    return seqs, p, skipped


out = open('/Users/benedictpaten/PycharmProjects/vntr-eval/work/full/chr20/repeat_sites/called_arms.tsv', 'w')
cols = ['region', 'arm', 'strand', 's', 'e', 'ref_len', 'n_panel_walks', 'n_panel_seqs', 'in_bed', 't_cross', 't_in_panel',
        'sites', 'nonref_sites', 'cross', 'skipped', 'under_nonref_parent', 'len', 'in_panel', 'panel_haps',
        'paired_truth', 'ed_truth', 'best_panel_ed_truth', 'nearest_panel_ed', 'ceiling_pair', 'call_pair', 'ref_pair']
out.write('\t'.join(cols) + '\n')
for k, r in enumerate(regs):
    s, e = A[r]
    a1, b1 = s - 1, e + 1
    refseq = REF[a1 - 1:b1]
    P = panel[r]
    (t1, t2), tp, _ = build(truth.get(r, []), a1, b1, refseq)
    tcross = tp.flags.get('records_crossing_span_edge', 0)
    tin = (t1 in P) + (t2 in P)
    ib = in_bed(a1, b1)
    ed = {}

    def E(x, y):
        key = (x, y) if x <= y else (y, x)
        v = ed.get(key)
        if v is None:
            v = ed[key] = S.edit_distance(x, y)
        return v
    calls = {}
    for arm, by in arms.items():
        recs = by.get(r, [])
        (c1, c2), p, sk = build(recs, a1, b1, refseq)
        calls[arm] = (c1, c2, p, sk, recs)
    need_panel = any(c not in P for c1, c2, *_ in calls.values() for c in (c1, c2))
    pseqs = list(P)
    if need_panel:
        pt = [(E(q, t1), E(q, t2)) for q in pseqs]
        ceil = min(min(pt[i][0] + pt[j][1], pt[i][1] + pt[j][0]) for i in range(len(pseqs)) for j in range(len(pseqs)))
        bp = (min(x[0] for x in pt), min(x[1] for x in pt))
    else:
        ceil, bp = '.', ('.', '.')
    refpair = '.'
    for arm, (c1, c2, p, sk, recs) in calls.items():
        sites = {x['id'].split('_')[0] for x in recs}
        nonref = {x['id'].split('_')[0] for x in recs if any(x['gt'])}
        a = (E(c1, t1), E(c2, t2))
        b = (E(c1, t2), E(c2, t1))
        pair = (0, 1) if sum(a) <= sum(b) else (1, 0)
        cp = min(sum(a), sum(b))
        for h, c in enumerate((c1, c2)):
            tj = pair[h]
            inp = c in P
            near = 0 if inp else min(E(c, q) for q in pseqs)
            row = [r, arm, h, s, e, b1 - a1 + 1, sum(len(v) for v in P.values()), len(P), int(ib), tcross, tin,
                   len(sites), len(nonref), p.flags.get('records_crossing_span_edge', 0), sk[h],
                   p.last_under_nonref_parent[h], len(c), int(inp), ','.join(P[c]) if inp else '.',
                   tj + 1, E(c, (t1, t2)[tj]), bp[tj], near, ceil, cp, refpair]
            out.write('\t'.join(map(str, row)) + '\n')
    if k % 1000 == 0:
        print(k, r, file=sys.stderr, flush=True)
out.close()
