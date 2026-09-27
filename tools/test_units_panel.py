#!/usr/bin/env python3
"""Tests of units_panel.py (the full-panel arm of unit_aware); about 10 s, needs mafft.

    python3 tools/test_units_panel.py

1. One real region end to end into a temporary root: the full MSA spells the union, the
   full-panel graph has one path per hprc.fa.gz record spelling it, the projected MSA and graph
   spell hap32.fa, and the realign.json / runtime row agree.
2. The process-tree killer takes down a child in a session of its own (as mafft runs).
"""
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402
import msa_graph  # noqa: E402
import panel  # noqa: E402
import units_panel  # noqa: E402

REGION = 'L000506'          # a non-TR control: fallback mode, ~2 s


@unittest.skipUnless(os.path.exists(os.path.join(config.REGIONS_DIR, REGION, 'hprc.fa.gz')) and
                     (shutil.which(config.MAFFT) or os.path.exists(config.MAFFT)), 'needs the region package and mafft')
class TestOneRegion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.makedirs(config.WORK_DIR, exist_ok=True)
        cls.td = tempfile.mkdtemp(prefix='test_units_panel.', dir=config.WORK_DIR)
        cls.row = units_panel.run_job(REGION, threads=1, timeout=600, cand_root=os.path.join(cls.td, 'cand'),
                                      panel_root=os.path.join(cls.td, 'panel'), workroot=os.path.join(cls.td, 'work'),
                                      record_runtime=False)
        cls.o = units_panel.out_paths(REGION, os.path.join(cls.td, 'cand'), os.path.join(cls.td, 'panel'))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.td, ignore_errors=True)

    def test_status_and_files(self):
        self.assertEqual(self.row['status'], 'ok', self.row)
        for k in ('msa', 'gfa', 'json', 'full_msa', 'full_gfa'):
            self.assertTrue(os.path.exists(self.o[k]), k)
        with open(self.o['json']) as f:
            info = json.load(f)
        self.assertEqual(info['status'], 'ok')
        self.assertEqual(info['panel'], 'all')
        self.assertEqual(self.row['panel'], 'all')
        self.assertEqual(int(self.row['n_distinct']), info['n_distinct'])

    def test_full_msa_spells_union(self):
        fa, mp = units_panel.union_files(REGION)
        want = dict((n, s.upper()) for n, _, s in panel.read_fasta(fa))
        tmp = os.path.join(self.td, 'full.msa.fa')
        with gzip.open(self.o['full_msa'], 'rb') as fi, open(tmp, 'wb') as fo:
            shutil.copyfileobj(fi, fo)
        rows = msa_graph.read_msa(tmp)
        self.assertEqual(sorted(n for n, _ in rows), sorted(want))
        self.assertEqual(len(set(len(r) for _, r in rows)), 1)
        for n, r in rows:
            self.assertEqual(r.replace('-', ''), want[n], n)

    def test_panel_graph_paths(self):
        hprc = [(n, s) for n, _, s in panel.read_fasta(os.path.join(config.REGIONS_DIR, REGION, 'hprc.fa.gz'))]
        self.assertEqual(units_panel.check_paths(self.o['full_gfa'], hprc), len(hprc))

    def test_projection_spells_hap32(self):
        hap = msa_graph.read_fasta(os.path.join(config.REGIONS_DIR, REGION, 'hap32.fa'))
        rows = msa_graph.read_msa(self.o['msa'])
        self.assertEqual([n for n, _ in rows], [n for n, _ in hap])
        for (n, r), (_, s) in zip(rows, hap):
            self.assertEqual(r.replace('-', ''), s.upper())
        L = len(rows[0][1])
        self.assertTrue(all(any(r[c] != '-' for _, r in rows) for c in range(L)), 'all-gap column kept')
        self.assertEqual(units_panel.check_paths(self.o['gfa'], hap), len(hap))

    def test_previous_status_skips(self):
        self.assertEqual(units_panel.previous_status(REGION, os.path.join(self.td, 'cand'),
                                                     os.path.join(self.td, 'panel')), 'ok')


class TestKillTree(unittest.TestCase):
    def test_kills_grandchild_in_own_session(self):
        code = ('import subprocess, time, sys\n'
                'p = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"], start_new_session=True)\n'
                'print(p.pid, flush=True)\n'
                'time.sleep(120)\n')
        p = subprocess.Popen([sys.executable, '-c', code], stdout=subprocess.PIPE, start_new_session=True, text=True)
        gpid = int(p.stdout.readline())
        tree = units_panel.tree_of(p.pid)
        self.assertIn(gpid, [x for x, _, _ in tree])
        units_panel.kill_tree(p)
        time.sleep(0.5)
        self.assertIsNotNone(p.poll())
        alive = True
        try:
            os.kill(gpid, 0)
        except ProcessLookupError:
            alive = False
        self.assertFalse(alive, 'grandchild in its own session survived')


if __name__ == '__main__':
    unittest.main(verbosity=2)
