#!/usr/bin/env python3
"""snarl_anchors.py -- truth-free region anchors from the snarl decomposition of a contig graph.

    python3 tools/snarl_anchors.py build --contig chr20 [--gfa orig.gfa] [--cache DIR]
    python3 tools/snarl_anchors.py query --contig chr20 START1 END1 [--pad 200] [--max-span 250000]

The anchors of a padded interval [a0, b0) (0-based, half-open on CHM13) are two boundary nodes of
ONE chain of the snarl tree (`vg snarls -T -P CHM13` on the contig GBZ, trivial snarls included,
so consecutive snarls of a chain share a boundary node): L = the chain's last boundary node that
ends at or before a0, R = its first boundary node that starts at or after b0, over every chain
(top-level or nested) that has both; the pair with the smallest CHM13 span wins. So the anchors are
the boundaries of the smallest snarl, or run of consecutive snarls of one chain, that encloses the
padded interval. Any two boundary nodes of one chain separate the graph, so everything between
them is exactly that run of snarls, whoever's paths are in it. Only boundary nodes that CHM13
visits once, in forward orientation, are used. No truth record is consulted.

`build` computes the snarls once (cached as <cache>/<contig>.snarls.json.gz) and writes
<cache>/<contig>.anchors.json.gz: per chain, the CHM13-forward boundary nodes with their CHM13
coordinates. Chains are the connected components of "shares a boundary node" over all snarls
(a nested chain's boundary nodes lie inside its parent snarl, so they never touch the parent's
chain).
"""
import argparse
import bisect
import collections
import gzip
import json
import os
import re
import subprocess
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402

DEFAULT_CACHE = os.path.join(config.WORK_DIR, 'stage4', 'snarls') if hasattr(config, 'WORK_DIR') else '.'
STEP_RE = re.compile(r'([<>])(\d+)')


class SnarlAnchorError(RuntimeError):
    pass


def cache_paths(contig, cache):
    return (os.path.join(cache, '%s.snarls.json.gz' % contig),
            os.path.join(cache, '%s.anchors.json.gz' % contig))


def compute_snarls(contig, out):
    """vg snarls -T -P CHM13 on the contig GBZ -> JSON lines (gzipped)."""
    gbz = config.data_paths(contig)['contig_gbz']
    pb = out + '.pb.tmp'
    with open(pb, 'wb') as f:
        subprocess.run([config.VG, 'snarls', '-T', '-P', 'CHM13', '-t', '2', gbz], stdout=f, check=True)
    with gzip.open(out + '.tmp', 'wb') as f:
        p = subprocess.Popen([config.VG, 'view', '-Rj', pb], stdout=subprocess.PIPE)
        for line in p.stdout:
            f.write(line)
        if p.wait():
            raise SnarlAnchorError('vg view -Rj failed')
    os.replace(out + '.tmp', out)
    os.remove(pb)


def chm13_walk(gfa):
    """[(node, orient '+'/'-', start0, end0)] along the CHM13 W line of a GFA."""
    lens, walk = {}, None
    with open(gfa) as f:
        for line in f:
            if line[0] == 'S':
                x = line.split('\t', 3)
                lens[int(x[1])] = len(x[2].rstrip('\n'))
            elif line[0] == 'W' and line.split('\t', 2)[1] == 'CHM13':
                if walk is not None:
                    raise SnarlAnchorError('more than one CHM13 W line in %s' % gfa)
                x = line.rstrip('\n').split('\t')
                walk = (int(x[4]), STEP_RE.findall(x[6]))
    if walk is None:
        raise SnarlAnchorError('no CHM13 W line in %s' % gfa)
    pos, out = walk[0], []
    for o, n in walk[1]:
        n = int(n)
        out.append((n, '+' if o == '>' else '-', pos, pos + lens[n]))
        pos += lens[n]
    return out


def build(contig, gfa, cache):
    os.makedirs(cache, exist_ok=True)
    sj, aj = cache_paths(contig, cache)
    if not os.path.exists(sj):
        compute_snarls(contig, sj)
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    n_snarls = 0
    with gzip.open(sj, 'rt') as f:
        for line in f:
            s = json.loads(line)
            a, b = int(s['start']['node_id']), int(s['end']['node_id'])
            n_snarls += 1
            for x in (a, b):
                parent.setdefault(x, x)
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    walk = chm13_walk(gfa)
    visits = collections.Counter(n for n, _, _, _ in walk)
    chains = collections.defaultdict(list)
    for n, o, st, en in walk:
        if n in parent and o == '+' and visits[n] == 1:
            chains[find(n)].append((st, en, n))
    out = {'contig': contig, 'snarls': os.path.abspath(sj), 'gfa': os.path.abspath(gfa),
           'n_snarls': n_snarls, 'n_boundary_nodes': len(parent),
           'chains': [sorted(v) for v in chains.values() if len(v) >= 2]}
    with gzip.open(aj + '.tmp', 'wt') as f:
        json.dump(out, f)
    os.replace(aj + '.tmp', aj)
    return aj


class Anchors:
    """Query object over <cache>/<contig>.anchors.json.gz."""

    def __init__(self, contig, cache=DEFAULT_CACHE):
        aj = cache_paths(contig, cache)[1]
        if not os.path.exists(aj):
            raise SnarlAnchorError('no snarl anchor cache %s (run snarl_anchors.py build)' % aj)
        d = json.load(gzip.open(aj, 'rt'))
        self.contig = contig
        pts = []
        for ci, ch in enumerate(d['chains']):
            for st, en, n in ch:
                pts.append((st, en, n, ci))
        pts.sort()
        self.pts = pts
        self.starts = [p[0] for p in pts]
        self.ends_sorted = sorted((p[1], k) for k, p in enumerate(pts))
        self.ends = [e for e, _ in self.ends_sorted]

    def query(self, a0, b0, max_span=250000):
        """(L, R, info): signed node ids (always +) of the tightest chain pair around [a0, b0)."""
        # nearest left boundary per chain: largest end <= a0, within max_span of b0
        left = {}
        i = bisect.bisect_right(self.ends, a0)
        lo_lim = b0 - max_span
        while i > 0:
            i -= 1
            e, k = self.ends_sorted[i]
            st, en, n, ci = self.pts[k]
            if st < lo_lim:
                if e < lo_lim:
                    break
                continue
            if ci not in left or en > left[ci][1]:
                left[ci] = (st, en, n)
        right = {}
        j = bisect.bisect_left(self.starts, b0)
        hi_lim = a0 + max_span
        while j < len(self.pts):
            st, en, n, ci = self.pts[j]
            if st > hi_lim:
                break
            if en <= hi_lim and (ci not in right or st < right[ci][0]):
                right[ci] = (st, en, n)
            j += 1
        best = None
        for ci in set(left) & set(right):
            L, R = left[ci], right[ci]
            span = R[1] - L[0]
            if span <= max_span and (best is None or span < best[0]):
                best = (span, L, R, ci)
        if best is None:
            # report the tightest pair ignoring max_span, for the skip reason
            wide = self._unbounded(a0, b0)
            raise SnarlAnchorError('no snarl chain encloses %d-%d within --max-span %d%s' % (
                a0, b0, max_span, '' if wide is None else ' (tightest enclosing span %d bp, L=%d R=%d)' % wide))
        span, L, R, ci = best
        return L[2], R[2], {'span_bp': span, 'L': list(L), 'R': list(R), 'chain': ci,
                            'rule': 'snarl chain boundaries (vg snarls -T -P CHM13), truth-free'}

    def _unbounded(self, a0, b0):
        left, right = {}, {}
        for st, en, n, ci in self.pts:
            if en <= a0:
                if ci not in left or en > left[ci][1]:
                    left[ci] = (st, en, n)
            elif st >= b0 and (ci not in right or st < right[ci][0]):
                right[ci] = (st, en, n)
        spans = [(right[c][1] - left[c][0], left[c][0], right[c][1]) for c in set(left) & set(right)]
        return min(spans) if spans else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('build')
    b.add_argument('--contig', required=True)
    b.add_argument('--gfa', required=True, help='GFA dump of the contig GBZ (vg convert -f)')
    b.add_argument('--cache', default=DEFAULT_CACHE)
    q = sub.add_parser('query')
    q.add_argument('--contig', required=True)
    q.add_argument('start1', type=int)
    q.add_argument('end1', type=int)
    q.add_argument('--pad', type=int, default=200)
    q.add_argument('--max-span', type=int, default=250000)
    q.add_argument('--cache', default=DEFAULT_CACHE)
    a = ap.parse_args(argv)
    if a.cmd == 'build':
        print(build(a.contig, a.gfa, a.cache))
    else:
        L, R, info = Anchors(a.contig, a.cache).query(max(0, a.start1 - 1 - a.pad), a.end1 + a.pad, a.max_span)
        print('%d+\t%d+\t%s' % (L, R, json.dumps(info)))


if __name__ == '__main__':
    main()
