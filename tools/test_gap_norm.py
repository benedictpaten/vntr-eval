#!/usr/bin/env python3
"""Tests for tools/gap_norm.py (gap normalisation of an MSA).

    python3 tools/test_gap_norm.py
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gap_norm  # noqa: E402


class TestNormalise(unittest.TestCase):
    def test_deletion_in_a_repeat_moves_to_the_left_end(self):
        rows = ['ACGTACGTAA', 'ACGT----AA']
        out, st = gap_norm.normalise(rows, 'left')
        self.assertEqual(out, ['ACGTACGTAA', '----ACGTAA'])
        self.assertEqual(st['moves'], 4)

    def test_right_slides_to_the_right_end(self):
        out, _ = gap_norm.normalise(['ACGTACGTAA', '----ACGTAA'], 'right')
        self.assertEqual(out, ['ACGTACGTAA', 'ACGTA----A'])   # one past the repeat: the next A matches too

    def test_a_mismatching_base_stops_the_shift(self):
        out, _ = gap_norm.normalise(['ACGTACGTAA', 'ACGT----AA', 'ACCT----AA'], 'left')
        self.assertEqual(out, ['ACGTACGTAA', '----ACGTAA', 'ACC----TAA'])

    def test_joint_keeps_a_shared_run_whole(self):
        out, _ = gap_norm.normalise_joint(['ACGTACGTAA', 'ACGT----AA', 'ACCT----AA'], 'left')
        self.assertEqual(out, ['ACGTACGTAA', 'ACG----TAA', 'ACC----TAA'])

    def test_sequences_kept_and_all_gap_columns_dropped(self):
        rows = ['AC-GT-ACGT', 'A--GTTACGT', 'ACG-T-AC-T']
        for fn in (gap_norm.normalise, gap_norm.normalise_joint):
            for d in ('left', 'right'):
                out, _ = fn(rows, d)
                self.assertEqual([r.replace('-', '') for r in out], [r.replace('-', '') for r in rows])
                self.assertTrue(all(any(r[c] != '-' for r in out) for c in range(len(out[0]))))


class TestKmers(unittest.TestCase):
    def test_a_shifted_copy_is_an_extra_position(self):
        k = 4
        same = gap_norm.kmer_excess(['ACGTTGCA--', 'ACGTTGCA--'], k)
        self.assertEqual(same['extra'], 0)
        shifted = gap_norm.kmer_excess(['ACGTTGCA--', '--ACGTTGCA'], k)
        self.assertEqual(shifted['extra'], 5)          # each of the 5 4-mers sits at two columns


if __name__ == '__main__':
    unittest.main()
