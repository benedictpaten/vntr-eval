#!/usr/bin/env python3
"""Offline unit tests for tools/region.py (no graph or VCF access).
    python3 tools/test_region.py"""
import random
import re
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import region as vr


def dp(a, b):
    prev = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        cur = [i] + [0] * len(b)
        for j in range(1, len(b) + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
        prev = cur
    return prev[-1]


def cigar_check(ref, hap, cig):
    i = j = cost = 0
    for n, op in re.findall(r'(\d+)([=XID])', cig):
        n = int(n)
        if op == '=':
            assert ref[i:i + n] == hap[j:j + n]
            i += n; j += n
        elif op == 'X':
            assert all(a != b for a, b in zip(ref[i:i + n], hap[j:j + n]))
            i += n; j += n; cost += n
        elif op == 'D':
            i += n; cost += n
        else:
            j += n; cost += n
    assert i == len(ref) and j == len(hap), (i, j, len(ref), len(hap))
    return cost


def test_edit_and_cigar():
    random.seed(7)
    for _ in range(2000):
        a = ''.join(random.choice('ACGT'[:random.randint(1, 4)]) for _ in range(random.randint(0, 30)))
        b = ''.join(random.choice('ACGT') for _ in range(random.randint(0, 30)))
        d = dp(a, b)
        assert vr.edit_distance(a, b) == d
        dist, cig = vr.align_cigar(a, b)
        assert dist == d
        assert cigar_check(a, b, cig) == d
    # budget cut-off still gives the exact distance
    d, cig = vr.align_cigar('ACGT' * 10, 'ACG' * 10, max_cells=10)
    assert cig is None and d == dp('ACGT' * 10, 'ACG' * 10)


def test_walks():
    gfa = vr.Gfa('S\t1\tAC\nS\t2\tGGT\nS\t3\tA\nL\t1\t+\t2\t-\t0M\n'
                 'W\tx\t1\tc\t0\t6\t>1<2>3\tWT:i:3\n')
    w = gfa.walks[0]
    assert w['steps'] == (1, -2, 3) and w['weight'] == 3
    assert gfa.walk_seq(w['steps']) == 'AC' + 'ACC' + 'A'
    assert vr.reverse_steps((1, -2, 3)) == (-3, 2, -1)
    assert gfa.walk_seq(vr.reverse_steps(w['steps'])) == vr.revcomp('ACACCA')
    reg = vr.Region(gfa, (1, -2, 3), 0, 1, 3)
    st, fl = reg.normalise((-3, 2, -1))   # a reverse traversal is flipped back
    assert st == (1, -2, 3) and fl
    assert reg.classify((1, -2, 3)) == 'spanning'
    assert reg.classify((1, -2)) == 'enters_L'
    assert reg.classify((-2, 3)) == 'exits_R'
    m = reg.walk_metrics((1, 2, 1, 3))
    assert m['max_node_visits'] == 2 and m['n_nodes_revisited'] == 1


def test_cycle_detection():
    assert not vr.graph_has_directed_cycle([1, 2, 3], [(1, '+', 2, '+'), (2, '+', 3, '+')])
    assert vr.graph_has_directed_cycle([1, 2], [(1, '+', 2, '+'), (2, '+', 1, '+')])
    assert vr.graph_has_directed_cycle([1], [(1, '+', 1, '+')])


def test_apply():
    ref = 'AAAACCCCGGGG'
    recs = [{'pos': 3, 'ref': 'AAC', 'alt': 'A'}, {'pos': 10, 'ref': 'G', 'alt': 'GTT'}]
    assert vr.apply_python(ref, 1, recs) == 'AAACCCGTTGGG'[:0] + 'AA' + 'A' + 'CCCG' + 'GTT' + 'GG'


def test_bp_union():
    assert vr.bp_union([(0, 10), (5, 20), (30, 40)], 0, 100) == 30
    assert vr.bp_union([(0, 10)], 5, 8) == 3


def test_runs_in_walk():
    # subgraph = {10 (L), 11, 12, 13 (R)}; outside nodes 1, 2, 9, 99
    S = {10, 11, 12, 13}
    tok = vr.boundary_token_re(10, 13)
    # full traversal, then a second pass through L that returns out via L (cycle through boundary)
    w = '>1>10>11>13>2>10>12<10>9'
    runs = vr.runs_in_walk(w, S, tok)
    assert [r[1] for r in runs] == [(10, 11, 13), (10, 12, -10)], runs
    # reverse traversal and a path fragment that starts inside and ends inside
    assert [r[1] for r in vr.runs_in_walk('<13<12<10>1', S, tok)] == [(-13, -12, -10)]
    runs = vr.runs_in_walk('>11>12', S, tok)
    assert runs == [(0, (11, 12), True)]
    # node ids that merely share a prefix with a boundary are not boundaries
    assert vr.runs_in_walk('>101>1300>2', S, tok) == []
    # a run reaching the end of the path is flagged
    assert vr.runs_in_walk('>2>10>11', S, tok)[0][2] is True


def test_tandem_scan():
    random.seed(3)
    rnd = lambda n: ''.join(random.choice('ACGT') for _ in range(n))
    motif = rnd(23)
    seq = rnd(500) + motif * 20 + rnd(500)
    r = vr.tandem_scan(seq)
    assert r['vntr']['period'] == 23, r
    assert 440 <= r['vntr']['covered_bp'] <= 520, r
    r = vr.tandem_scan(rnd(2000))
    assert r['vntr'] is None or r['vntr']['covered_bp'] < 100, r
    r = vr.tandem_scan(rnd(300) + 'CA' * 60 + rnd(300))
    assert r['str']['period'] == 2, r


def test_c_edit_distance_matches_python():
    random.seed(11)
    for _ in range(300):
        a = ''.join(random.choice('ACGTN') for _ in range(random.randint(0, 200)))
        b = ''.join(random.choice('ACGTN') for _ in range(random.randint(0, 200)))
        assert vr.edit_distance(a, b) == vr.edit_distance_py(a, b)
    a = 'ACGT' * 300
    b = 'ACGA' * 290
    assert vr.edit_distance(a, b) == vr.edit_distance_py(a, b)


def test_parse_anchor():
    assert vr.parse_anchor('123+') == 123 and vr.parse_anchor('123-') == -123
    assert vr.parse_anchor('>123') == 123 and vr.parse_anchor('<123') == -123
    try:
        vr.parse_anchor('12x')
    except vr.ToolError:
        pass
    else:
        raise AssertionError('bad anchor accepted')


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn()
            print('ok', name)
