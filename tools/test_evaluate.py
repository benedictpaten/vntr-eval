#!/usr/bin/env python3
"""test_evaluate.py -- offline tests of the exact DPs and metric helpers behind evaluate.py.

    python3 tools/test_evaluate.py

Checks the C routines against plain-Python dynamic programs on random inputs: unit edit distance,
gap-affine cost (including the band-doubling exactness bound), and the query-vs-DAG distance with
its traceback (against brute force over every source->sink path of small random DAGs, with and
without the upper-bound pruning); the path
name matching of vg-written W lines and the traceback's tie-break. Then runs
evaluate.py end to end on a synthetic region whose graph is known (a SNP bubble and a
misaligned repeat) and checks the headline numbers.
"""
import itertools
import json
import os
import random
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evaluate as E  # noqa: E402


def ed_py(a, b):
    prev = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        cur = [i] + [0] * len(b)
        for j in range(1, len(b) + 1):
            cur[j] = min(prev[j - 1] + (a[i - 1] != b[j - 1]), prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return prev[-1]


def affine_py(a, b, X=4, GO=6, GE=2):
    n, m = len(a), len(b)
    INF = 1 << 40
    if n == 0 and m == 0:
        return 0
    H = [0] + [GO + GE * j for j in range(1, m + 1)]
    F = [INF] * (m + 1)
    for i in range(1, n + 1):
        Hn = [GO + GE * i] + [0] * m
        Fn = [INF] * (m + 1)
        Ev = INF
        for j in range(1, m + 1):
            Fn[j] = min(F[j] + GE, H[j] + GO + GE)
            Ev = min(Ev + GE, Hn[j - 1] + GO + GE)
            Hn[j] = min(H[j - 1] + (0 if a[i - 1] == b[j - 1] else X), Fn[j], Ev)
        H, F = Hn, Fn
    return H[m]


def rand_seq(rng, n, alpha='ACGT'):
    return ''.join(rng.choice(alpha) for _ in range(n))


def mutate(rng, s, rate):
    out = []
    for c in s:
        r = rng.random()
        if r < rate / 3:
            continue
        if r < 2 * rate / 3:
            out.append(rng.choice('ACGT'))
        elif r < rate:
            out.append(c + rng.choice('ACGT'))
        else:
            out.append(c)
    return ''.join(out)


def test_ed_affine():
    C = E.clib()
    rng = random.Random(7)
    for t in range(400):
        n = rng.randint(0, 90)
        a = rand_seq(rng, n, 'ACGTN' if t % 5 == 0 else 'ACGT')
        b = mutate(rng, a, rng.choice([0.0, 0.05, 0.3, 0.9])) if t % 3 else rand_seq(rng, rng.randint(0, 150))
        if t % 7 == 0:
            b = b + rand_seq(rng, rng.randint(50, 200))     # long length difference
        assert C.ed(a.encode(), b.encode()) == ed_py(a, b), (a, b)
        assert C.affine(a.encode(), b.encode()) == affine_py(a, b), (a, b, C.affine(a.encode(), b.encode()),
                                                                     affine_py(a, b))
    # long sequences that force band doubling: repeats with a shifted expansion
    for t in range(6):
        u = rand_seq(rng, rng.randint(5, 40))
        a = u * rng.randint(10, 30)
        b = mutate(rng, u * rng.randint(3, 30), 0.1)
        assert C.ed(a.encode(), b.encode()) == ed_py(a, b)
        assert C.affine(a.encode(), b.encode()) == affine_py(a, b)
    print('ed / affine: ok')


def write_gfa(path, nodes, edges, paths):
    with open(path, 'w') as f:
        f.write('H\tVN:Z:1.0\n')
        for i, s in nodes.items():
            f.write('S\t%s\t%s\n' % (i, s))
        for a, ao, b, bo in edges:
            f.write('L\t%s\t%s\t%s\t%s\t0M\n' % (a, ao, b, bo))
        for n, steps in paths.items():
            f.write('P\t%s\t%s\t*\n' % (n, ','.join(steps)))


def all_paths(succ, src, snk, cur, acc, out):
    acc.append(cur)
    if cur in snk:
        out.append(list(acc))
    for w in succ.get(cur, ()):
        all_paths(succ, src, snk, w, acc, out)
    acc.pop()


def test_dag():
    rng = random.Random(11)
    with tempfile.TemporaryDirectory() as td:
        for t in range(120):
            k = rng.randint(2, 9)
            nodes = {str(i + 1): rand_seq(rng, rng.randint(1, 6)) for i in range(k)}
            edges = []
            for i in range(1, k + 1):
                for j in range(i + 1, k + 1):
                    if j == i + 1 or rng.random() < 0.3:
                        edges.append((str(i), '+', str(j), '+'))
            # a few reverse-strand traversals: node j entered on its reverse strand
            if t % 4 == 0 and k >= 3:
                j = rng.randint(2, k - 1)
                edges = [e for e in edges if e[2] != str(j) and e[0] != str(j)]
                edges.append((str(j - 1), '+', str(j), '-'))
                edges.append((str(j), '-', str(j + 1), '+'))
            path = ['1+']
            # a path along the chain (possibly through the reversed node)
            succ = {}
            for a, ao, b, bo in edges:
                succ.setdefault((a, ao), []).append((b, bo))
            cur = ('1', '+')
            while cur[0] != str(k):
                cur = rng.choice(succ[cur])
                path.append(cur[0] + cur[1])
            gfa = os.path.join(td, 'g.gfa')
            write_gfa(gfa, nodes, edges, {'p': path})
            g = E.Graph(gfa)
            hs = g.paths['p']
            dag = E.build_dag(g, [hs])
            assert dag['acyclic']
            # brute force over every source->sink handle path
            hsucc = {}
            for h in dag['order']:
                hsucc[h] = [w for w in dag['succ'][h] if w in dag['keep']]
            cand = []
            for s in dag['sources']:
                all_paths(hsucc, s, dag['sinks'], s, [], cand)
            q = mutate(rng, g.spell(hs), 0.3) if t % 2 else rand_seq(rng, rng.randint(0, 25))
            if not q:
                q = 'A'
            best = min(ed_py(g.spell(p), q) for p in cand)
            d, hp, err = E.dag_align(g, dag, q, 1e9)
            assert err is None, err
            assert d == best, (d, best)
            assert ed_py(g.spell(hp), q) == d
            assert hp[0] in dag['sources'] and hp[-1] in dag['sinks']
            for a, b in zip(hp, hp[1:]):
                assert b in dag['succ'][a]
            # the banded DP (cells above an upper bound pruned) is exact for any valid bound
            for ub in (best, best + rng.randint(0, 5), 2 * best + 3):
                d2, hp2, err2 = E.dag_align(g, dag, q, 1e9, None, ub)
                assert err2 is None and d2 == best, (ub, d2, best, err2)
                assert ed_py(g.spell(hp2), q) == best
            if best > 0:
                d3, hp3, err3 = E.dag_align(g, dag, q, 1e9, None, best - 1)
                assert d3 is None and err3 == 'upper_bound_violated', (d3, err3)
    # larger bubbly DAGs: banded (ub = distance to a panel path) == unbanded, and the traceback's
    # sequence is at that distance
    with tempfile.TemporaryDirectory() as td:
        for t in range(25):
            nodes, edges, walks = {}, [], []
            nid = 1
            nodes[str(nid)] = rand_seq(rng, 20)
            prev = [str(nid)]
            nw = 6
            walks = [['1+'] for _ in range(nw)]
            for b in range(rng.randint(3, 12)):
                alle = []
                base = rand_seq(rng, rng.randint(0, 30))
                for a in range(rng.randint(1, 4)):
                    nid += 1
                    nodes[str(nid)] = mutate(rng, base, 0.3) or 'A'
                    alle.append(str(nid))
                nid += 1
                nodes[str(nid)] = rand_seq(rng, rng.randint(1, 8))
                join = str(nid)
                for a in alle:
                    for p_ in prev:
                        edges.append((p_, '+', a, '+'))
                    edges.append((a, '+', join, '+'))
                for w in walks:
                    w.extend([rng.choice(alle) + '+', join + '+'])
                prev = [join]
            gfa = os.path.join(td, 'b.gfa')
            write_gfa(gfa, nodes, edges, {'w%d' % i: w for i, w in enumerate(walks)})
            g = E.Graph(gfa)
            hps = list(g.paths.values())
            dag = E.build_dag(g, hps)
            q = mutate(rng, g.spell(rng.choice(hps)), rng.choice([0.02, 0.1, 0.4]))
            ub = min(ed_py(g.spell(h), q) for h in hps)
            d1, hp1, e1 = E.dag_align(g, dag, q, 1e9)
            d2, hp2, e2 = E.dag_align(g, dag, q, 1e9, None, ub)
            assert e1 is None and e2 is None and d1 == d2 <= ub, (d1, d2, ub, e1, e2)
            assert ed_py(g.spell(hp1), q) == d1 and ed_py(g.spell(hp2), q) == d2
    print('dag_align: ok')


def test_end_to_end():
    """A synthetic region: flank + SNP bubble + a 3-unit repeat whose extra unit the 'bad' graph
    places at a different end in different paths."""
    rng = random.Random(3)
    fl, fr = rand_seq(rng, 60), rand_seq(rng, 60)
    unit = 'ACGTTGCAAT'
    ref = fl + 'A' + unit * 3 + fr
    h1 = fl + 'G' + unit * 3 + fr
    h2 = fl + 'A' + unit * 4 + fr
    h3 = fl + 'A' + unit * 4 + fr
    with tempfile.TemporaryDirectory() as td:
        rd = os.path.join(td, 'R1')
        os.makedirs(rd)
        json.dump({'region_id': 'R1', 'stratum': 'test', 'contig': 'chrT', 'core_start': 62, 'core_end': 91,
                   'span_start': 1, 'span_end': len(ref), 'anchor_left': '>1', 'anchor_right': '>9'},
                  open(os.path.join(rd, 'region.json'), 'w'))
        with open(os.path.join(rd, 'hap32.fa'), 'w') as f:
            for n, s in (('CHM13#0#chrT', ref), ('h1', h1), ('h2', h2), ('h3', h3)):
                f.write('>%s\n%s\n' % (n, s))
        with open(os.path.join(rd, 'truth.fa'), 'w') as f:
            f.write('>t1\n%s\n>t2\n%s\n' % (h1, fl + 'G' + unit * 4 + fr))
        # good graph: extra unit inserted at the end of the array in both h2 and h3 (one node)
        nodes = {'1': fl, '2': 'A', '3': 'G', '4': unit * 3, '5': unit, '9': fr}
        edges = [('1', '+', '2', '+'), ('1', '+', '3', '+'), ('2', '+', '4', '+'), ('3', '+', '4', '+'),
                 ('4', '+', '5', '+'), ('4', '+', '9', '+'), ('5', '+', '9', '+')]
        paths = {'CHM13#0#chrT': ['1+', '2+', '4+', '9+'], 'h1': ['1+', '3+', '4+', '9+'],
                 'h2': ['1+', '2+', '4+', '5+', '9+'], 'h3': ['1+', '2+', '4+', '5+', '9+']}
        good = os.path.join(td, 'good.gfa')
        write_gfa(good, nodes, edges, paths)
        # misaligned graph: h3 puts its extra unit at the start of the array, h2 at the end
        nodes_b = dict(nodes)
        nodes_b['7'] = unit
        edges_b = edges + [('2', '+', '7', '+'), ('7', '+', '4', '+')]
        paths_b = dict(paths)
        paths_b['h3'] = ['1+', '2+', '7+', '4+', '9+']
        bad = os.path.join(td, 'bad.gfa')
        write_gfa(bad, nodes_b, edges_b, paths_b)
        # under-aligned graph: h3 carries its 4 units on one parallel node (unaligned homology)
        nodes_p = dict(nodes)
        nodes_p['6'] = unit * 4
        edges_p = edges + [('2', '+', '6', '+'), ('6', '+', '9', '+')]
        paths_p = dict(paths)
        paths_p['h3'] = ['1+', '2+', '6+', '9+']
        par = os.path.join(td, 'par.gfa')
        write_gfa(par, nodes_p, edges_p, paths_p)
        # invalid graph: h1 spelled differently
        paths_i = dict(paths)
        paths_i['h1'] = ['1+', '2+', '4+', '9+']
        inv = os.path.join(td, 'inv.gfa')
        write_gfa(inv, nodes, edges, paths_i)
        env = dict(os.environ, VNTR_WORK=os.path.join(td, 'work'))
        out = {}
        for tag, gfa in (('good', good), ('bad', bad), ('par', par), ('inv', inv)):
            o = os.path.join(td, tag + '.json')
            p = subprocess.run([sys.executable, os.path.join(HERE, 'evaluate.py'), rd, gfa, '--name', tag,
                                '--out', o, '--no-truvari'], capture_output=True, text=True, env=env)
            out[tag] = (p.returncode, json.load(open(o)))
        rc, gd = out['good']
        assert rc == 0, gd
        # vs CHM13 every path is optimally aligned; one non-reference pair (G+3 units vs A+4 units)
        # has a cheaper pairwise optimum than any single alignment of all four (10 vs 11 edits)
        assert gd['alignment']['ref_cost_over_opt'] == 1.0, gd['alignment']
        assert gd['alignment']['all_pair_ratio_max'] == 1.1
        # a repeat-unit insertion leaves a floor: 'last unit + flank' k-mers sit on the CHM13 array
        # node and on the inserted unit's node (the metric's documented upper-bound caveat)
        assert gd['redundancy']['kmer_extra_positions'] <= 2 * len(unit)
        assert gd['truth']['h1_d_graph'] == 0 and gd['truth']['h2_d_graph'] == 0   # t2 = G + 4 units: a graph path
        assert gd['truth']['h2_d_panel'] == 1
        rc, bd = out['bad']
        assert rc == 0
        # h2 vs h3: identical sequence, but the graph stacks the extra unit at opposite ends: 20 edits
        assert bd['alignment']['all_pairs_same_seq_diff_walk'] == 1
        assert bd['alignment']['all_cost_over_opt'] > gd['alignment']['all_cost_over_opt']
        assert bd['size']['same_seq_different_walk'] == 1
        rc, pd = out['par']
        assert rc == 0
        # a parallel copy bracketed by shared nodes costs no extra edits but is unaligned homology
        assert pd['alignment']['all_U_per_kb'] > 0 and gd['alignment']['all_U_per_kb'] == 0
        assert pd['redundancy']['kmer_extra_positions'] > gd['redundancy']['kmer_extra_positions'] + 20
        rc, iv = out['inv']
        assert rc == 3 and iv['status'] == 'invalid' and iv['validity']['paths_sequence_mismatch'] == ['h1']
    print('end to end: ok')


def test_names_and_tiebreak():
    """W-line names as vg writes them, and the CHM13-first traceback tie-break"""
    with tempfile.TemporaryDirectory() as td:
        gfa = os.path.join(td, 'w.gfa')
        with open(gfa, 'w') as f:
            f.write('H\tVN:Z:1.0\nS\t1\tACGT\nS\t2\tGG\nS\t3\tGG\nS\t4\tTTAA\n')
            for a, b in (('1', '2'), ('1', '3'), ('2', '4'), ('3', '4')):
                f.write('L\t%s\t+\t%s\t+\t0M\n' % (a, b))
            # listed so that the non-reference allele comes first in every predecessor list
            f.write('W\tsampleA\t1\tchr1\t0\t10\t>1>3>4\n')
            f.write('W\tGRCh38\t0\tchr1\t5000\t5010\t>1>3>4\n')
            f.write('W\tCHM13\t0\tchr1\t0\t10\t>1>2>4\n')
            f.write('P\tsampleB#2#chr1#1\t1+,3+,4+\t*\n')
        g = E.Graph(gfa)
        names = ['CHM13#0#chr1', 'GRCh38#0#chr1[5000]', 'sampleA#1#chr1#0', 'sampleB#2#chr1#1']
        mp, miss = E.match_paths(g, names)
        assert not miss, miss
        assert mp['GRCh38#0#chr1[5000]'] == 'GRCh38#0#chr1' and mp['sampleA#1#chr1#0'] == 'sampleA#1#chr1'
        paths = {n: g.paths[mp[n]] for n in names}
        dag = E.build_dag(g, list(paths.values()))
        ref = paths['CHM13#0#chr1']
        d, hp, err = E.dag_align(g, dag, 'ACGTGGTTAA', 1e9, {h: (1 << 20) for h in ref})
        assert d == 0 and hp == ref, (d, hp, ref)
        d, hp, err = E.dag_align(g, dag, 'ACGTGGTTAA', 1e9, {h: (1 << 20) for h in paths['sampleA#1#chr1#0']})
        assert d == 0 and hp == paths['sampleA#1#chr1#0']
    print('names / tie-break: ok')


if __name__ == '__main__':
    test_ed_affine()
    test_dag()
    test_names_and_tiebreak()
    test_end_to_end()
    print('all tests passed')
