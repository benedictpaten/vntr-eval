#!/usr/bin/env python3
"""scan_tr_regions.py -- Stage 4 step A: one pass over a contig's hap32 GBZ, per TR region.

    python3 tools/scan_tr_regions.py --catalogue CAT.tsv --contig chr20 --gbz chr20.gbz \
        --out scan.tsv [--merge 100] [--pad 100] [--min-diff 50] [--bed prefilter.bed]

Catalogue rows (region_id chrom start end length tr_class period motif ...; start 0-based,
end exclusive) of one contig are merged when the gap between them is <= --merge bp (the merged
region keeps the first member's id, and the class, period and motif of its longest member).
Each merged region is padded by --pad bp and bracketed by the CHM13 nodes holding the first
and the last padded base. One streaming pass over `vg convert -f` (S lines, then W lines) gives,
per region and per path (the hap32 haplotype fragments, CHM13 and GRCh38) that visits both
bracket nodes in order (forward L..R, or reverse -R..-L; the first L, then the first R after
it), the bp strictly between the brackets; CHM13's own length is the same quantity on CHM13.
Written per region: n paths spanning, distinct walks (node sequence between the brackets),
min/max length, max |length - CHM13 length|. The prefilter (--min-diff) is max diff >= 50 bp.
Standard library only.
"""
import argparse
import bisect
import collections
import os
import subprocess
import sys
from array import array

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402


def load_catalogue(path, contig):
    rows = []
    with open(path) as f:
        h = f.readline().rstrip('\n').split('\t')
        for line in f:
            x = line.rstrip('\n').split('\t')
            if x[1] != contig:
                continue
            r = dict(zip(h, x))
            rows.append((int(r['start']), int(r['end']), r))
    rows.sort(key=lambda t: (t[0], t[1]))
    return rows


def merge(rows, gap):
    out = []
    for s, e, r in rows:
        if out and s - out[-1]['end'] <= gap:
            m = out[-1]
            m['end'] = max(m['end'], e)
            m['members'].append(r)
        else:
            out.append({'start': s, 'end': e, 'members': [r]})
    for m in out:
        big = max(m['members'], key=lambda r: int(r['end']) - int(r['start']))
        m['region_id'] = m['members'][0]['region_id']
        m['tr_class'], m['period'], m['motif'] = big['tr_class'], big['period'], big['motif']
    return out


def walk_ids(walk):
    return list(map(int, walk.replace('>', ' ').replace('<', ' -').split()))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalogue', required=True)
    ap.add_argument('--contig', required=True)
    ap.add_argument('--gbz', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--bed', help='write the prefiltered regions as BED (contig start0 end id period motif)')
    ap.add_argument('--merge', type=int, default=100)
    ap.add_argument('--pad', type=int, default=100)
    ap.add_argument('--min-diff', type=int, default=50)
    ap.add_argument('--ref-sample', default='CHM13')
    a = ap.parse_args(argv)

    regs = merge(load_catalogue(a.catalogue, a.contig), a.merge)
    print('%d catalogue regions merged into %d' % (sum(len(m['members']) for m in regs), len(regs)),
          file=sys.stderr)
    p = subprocess.Popen([config.VG, 'convert', '-f', a.gbz], stdout=subprocess.PIPE, text=True,
                         bufsize=1 << 22, env=config.tool_env())
    ids_s, lens_s = array('q'), array('i')
    walks = []                       # (path name, ids) for every W line; CHM13 first resolved later
    nodelen = None
    for line in p.stdout:
        t = line[0]
        if t == 'S':
            f = line.split('\t', 3)
            ids_s.append(int(f[1]))
            lens_s.append(len(f[2].rstrip('\n')))
        elif t == 'W':
            if nodelen is None:
                off = min(ids_s)
                nodelen = array('i', [0]) * (max(ids_s) - off + 1)
                for i, ln in zip(ids_s, lens_s):
                    nodelen[i - off] = ln
                del ids_s, lens_s
            f = line.rstrip('\n').split('\t')
            name = '%s#%s#%s[%s]' % (f[1], f[2], f[3], f[4])
            walks.append((name, f[1] == a.ref_sample, f[6]))
    if p.wait() != 0:
        raise SystemExit('vg convert failed')
    ref = [w for w in walks if w[1]]
    if len(ref) != 1:
        raise SystemExit('expected one %s walk, found %d' % (a.ref_sample, len(ref)))

    # CHM13 node starts -> bracket nodes
    rids = walk_ids(ref[0][2])
    if any(i < 0 for i in rids):
        raise SystemExit('reverse step on the reference walk')
    rstart, pos = [], 0
    for i in rids:
        rstart.append(pos)
        pos += nodelen[i - off]
    clen = pos
    bnodes = set()
    for m in regs:
        ps, pe = max(0, m['start'] - a.pad), min(clen, m['end'] + a.pad)
        li = bisect.bisect_right(rstart, ps) - 1
        ri = bisect.bisect_right(rstart, pe - 1) - 1
        m.update(pad_start=ps, pad_end=pe, L=rids[li], R=rids[ri],
                 chm13_len=rstart[ri] - (rstart[li] + nodelen[rids[li] - off]), lens=[], walks=set())
        if li != ri:
            bnodes.add(rids[li])
            bnodes.add(rids[ri])
    by_L = collections.defaultdict(list)
    for k, m in enumerate(regs):
        if m['L'] != m['R']:
            by_L[m['L']].append(k)
    del rids, rstart

    for name, _, walk in walks:
        ids = walk_ids(walk)
        occ = collections.defaultdict(list)       # signed node -> [step index]
        for k, i in enumerate(ids):
            if (i if i > 0 else -i) in bnodes:
                occ[i].append(k)
        if not occ:
            continue
        cum = [0]
        acc = 0
        for i in ids:
            acc += nodelen[(i if i > 0 else -i) - off]
            cum.append(acc)                       # cum[k] = start of step k, cum[k+1] = its end
        for Lnode in {abs(s) for s in occ} & set(by_L):
            for k in by_L[Lnode]:
                m = regs[k]
                L, R = m['L'], m['R']
                got = None
                if occ.get(L) and occ.get(R):
                    iL = occ[L][0]
                    j = bisect.bisect_right(occ[R], iL)
                    if j < len(occ[R]):
                        iR = occ[R][j]
                        got = (cum[iR] - cum[iL + 1], tuple(ids[iL + 1:iR]))
                if got is None and occ.get(-R) and occ.get(-L):
                    iR = occ[-R][0]
                    j = bisect.bisect_right(occ[-L], iR)
                    if j < len(occ[-L]):
                        iL = occ[-L][j]
                        got = (cum[iL] - cum[iR + 1], tuple(-x for x in reversed(ids[iR + 1:iL])))
                if got:
                    m['lens'].append(got[0])
                    m['walks'].add(hash(got[1]))
        del ids, cum, occ

    n_pass = collections.Counter()
    n_all = collections.Counter()
    with open(a.out, 'w') as fo:
        fo.write('\t'.join(['region_id', 'contig', 'start', 'end', 'tr_class', 'period', 'motif',
                            'n_members', 'members', 'pad_start', 'pad_end', 'L', 'R', 'chm13_len',
                            'n_paths', 'n_walks', 'min_len', 'max_len', 'max_diff', 'prefilter']) + '\n')
        bed = open(a.bed, 'w') if a.bed else None
        for m in regs:
            lens = m['lens']
            md = max((abs(x - m['chm13_len']) for x in lens), default=0)
            ok = md >= a.min_diff
            n_all[m['tr_class']] += 1
            n_pass[m['tr_class']] += ok
            fo.write('\t'.join(map(str, [
                m['region_id'], a.contig, m['start'], m['end'], m['tr_class'], m['period'], m['motif'],
                len(m['members']), ','.join(r['region_id'] for r in m['members']), m['pad_start'],
                m['pad_end'], m['L'], m['R'], m['chm13_len'], len(lens), len(m['walks']),
                min(lens, default=''), max(lens, default=''), md, int(ok)])) + '\n')
            if bed and ok:
                bed.write('\t'.join(map(str, [a.contig, m['start'], m['end'], m['region_id'],
                                              m['period'], m['motif']])) + '\n')
        if bed:
            bed.close()
    for c in sorted(n_all):
        print('%s\t%d merged regions\t%d pass (max diff >= %d)' % (c, n_all[c], n_pass[c], a.min_diff),
              file=sys.stderr)
    print('total\t%d\t%d' % (sum(n_all.values()), sum(n_pass.values())), file=sys.stderr)


if __name__ == '__main__':
    main()
