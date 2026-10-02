#!/usr/bin/env python3
"""Tests for the nearest-neighbour threading aligners in tools/iterate.py (kinds mst, mst3, bbt).

    python3 tools/test_iterate_mst.py
"""
import collections
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import iterate as it  # noqa: E402
import msa_graph  # noqa: E402


def masks(s):
    return [it._BIT[b] for b in s]


class TestPlaceInsert(unittest.TestCase):
    def test_fits_into_existing_columns_at_the_left(self):
        self.assertEqual(it.place_insert('ACG', masks('ACGTTT')),
                         [(0, 0), (1, 1), (2, 2), (3, None), (4, None), (5, None)])

    def test_overflow_goes_to_new_columns(self):
        self.assertEqual(it.place_insert('ACGAAA', masks('ACG')),
                         [(0, 0), (1, 1), (2, 2), (None, 3), (None, 4), (None, 5)])

    def test_ties_go_left(self):
        self.assertEqual(it.place_insert('A', masks('AA')), [(0, 0), (1, None)])

    def test_matching_columns_are_found_past_unrelated_ones(self):
        ops = it.place_insert('GGG', masks('TTGGG'))
        self.assertEqual([j for j, t in ops if t is not None], [2, 3, 4])


class TestThreadedMSA(unittest.TestCase):
    def test_recurrent_insertions_stack_and_deletions_inherit(self):
        st = collections.Counter()
        m = it.ThreadedMSA('r', 'AAACCCGGG')
        m.add('a', 'AAACCCTTTGGG', 'r', 'AAACCC---GGG', 'AAACCCTTTGGG', st)
        m.add('b', 'AAACCCTTGGG', 'r', 'AAACCC--GGG', 'AAACCCTTGGG', st)      # reuses a's TT columns
        m.add('c', 'AAAGGG', 'a', 'AAACCCTTTGGG', 'AAA------GGG', st)
        m.add('d', 'AAACCCGGG', 'c', 'AAA---GGG', 'AAACCCGGG', st)            # CCC back in r's columns
        recs = [('r', 'AAACCCGGG'), ('a', 'AAACCCTTTGGG'), ('b', 'AAACCCTTGGG'), ('c', 'AAAGGG'),
                ('d', 'AAACCCGGG')]
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, 'm.fa')
            self.assertEqual(m.write(recs, f), 12)
            rows = dict(msa_graph.read_msa(f))
        self.assertEqual(rows, {'r': 'AAACCC---GGG', 'a': 'AAACCCTTTGGG', 'b': 'AAACCCTT-GGG',
                                'c': 'AAA------GGG', 'd': 'AAACCC---GGG'})
        self.assertEqual(st['new_columns'], 3)


class TestDistance(unittest.TestCase):
    def test_multiset_jaccard_sees_copy_number(self):
        unit = 'ACGTTGCAAGT'
        recs = [('s0', unit * 10), ('s1', unit * 6), ('s2', unit * 9)]
        order, parent, dist = it.mst_order(recs, k=5)
        self.assertEqual(order, [0, 2, 1])
        self.assertEqual(parent[1], 2)
        self.assertLess(dist(0, 2), dist(0, 1))

    def test_pair_cost_is_abpoa_convex(self):
        self.assertEqual(it.pair_cost('ACGT', 'AGGT'), 4)
        self.assertEqual(it.pair_cost('ACGTAC', 'A----C'), 12)            # min(4 + 2*4, 24 + 4)
        self.assertEqual(it.pair_cost('A' * 30, 'A' + '-' * 28 + 'A'), 52)   # min(4 + 56, 24 + 28)


class TestBackbone(unittest.TestCase):
    def test_kmedoids_fixed_medoid_stays_and_finds_the_other_cluster(self):
        pts = [0.0, 0.1, 0.2, 10.0, 10.1, 10.3, 10.35, 13.0]
        med, cl, _ = it.kmedoids(len(pts), lambda i, j: abs(pts[i] - pts[j]), 2, fixed=[2])
        self.assertEqual(med[0], 2)
        self.assertEqual(med[1], 5)                     # the far cluster's medoid, not the outlier it started at
        self.assertEqual(cl, [0, 0, 0, 1, 1, 1, 1, 1])

    def test_kmedoids_k_at_least_n_takes_everything(self):
        med, cl, _ = it.kmedoids(3, lambda i, j: abs(i - j), 5, fixed=[1])
        self.assertEqual(sorted(med), [0, 1, 2])
        self.assertEqual(med[0], 1)

    def test_from_msa_then_add_reuses_backbone_columns(self):
        st = collections.Counter()
        m = it.ThreadedMSA.from_msa(collections.OrderedDict([('a', 'AC--GT'), ('b', 'ACTTGT'), ('c', '------')]))
        m.add('d', 'ACTGT', 'a', 'AC-GT', 'ACTGT', st)   # insertion relative to a lands in b's T column
        recs = [('a', 'ACGT'), ('b', 'ACTTGT'), ('d', 'ACTGT')]
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, 'm.fa')
            self.assertEqual(m.write(recs, f), 6)
            rows = dict(msa_graph.read_msa(f))
        self.assertEqual(rows, {'a': 'AC--GT', 'b': 'ACTTGT', 'd': 'ACT-GT'})
        self.assertEqual(st['new_columns'], 0)

    def test_plan_keeps_panel_rows_independent_of_outside_rows(self):
        unit = 'ACGTTGCAAGT'
        recs = [('s0', unit * 10), ('s1', unit * 4), ('s2', unit * 9), ('s3', unit * 5), ('s4', unit * 4 + 'A')]
        reps, steps, info = it.bbt_plan(recs, [0, 1, 2, 3], None, 2, k=5)
        self.assertEqual(reps[0], 0)
        self.assertEqual(info['n_outside_panel'], 1)
        self.assertEqual(steps[-1][0], 4)                # the outside row is threaded last
        placed = set(reps)
        for x, p in steps:
            self.assertIn(p, placed)
            placed.add(x)
        self.assertEqual(placed, set(range(5)))


if __name__ == '__main__':
    unittest.main()
