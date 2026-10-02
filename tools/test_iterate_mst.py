#!/usr/bin/env python3
"""Tests for the nearest-neighbour threading aligner in tools/iterate.py (kinds mst, mst3).

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


if __name__ == '__main__':
    unittest.main()
