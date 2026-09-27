#!/usr/bin/env python3
"""Tests for tools/poa_panel.py (the full-panel arm of the POA realigners).

    python3 tools/test_poa_panel.py

Runs both POA methods end to end on the smallest packaged region (L000506: 27 distinct full-panel
sequences) into a temporary directory, plus the memory model, the refusal path, a variant order
and the runtime table. Needs regions/L000506, work/panel/union (made on the fly if missing),
abpoa and spoa.
"""
import csv
import json
import os
import shutil
import sys
import tempfile
import unittest

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402
import msa_graph  # noqa: E402
import panel  # noqa: E402
import poa_panel  # noqa: E402

RID = 'L000506'
HAVE = os.path.exists(os.path.join(config.REGIONS_DIR, RID, 'hap32.fa'))


class TestPredict(unittest.TestCase):
    def test_spoa_is_quadratic_in_the_longest(self):
        self.assertAlmostEqual(poa_panel.predict_mb('spoa', [1000, 10]), 50 + 80e6 / 2 ** 20, places=0)
        # only the longest matters
        self.assertEqual(poa_panel.predict_mb('spoa', [8000, 10]), poa_panel.predict_mb('spoa', [8000, 7999]))

    def test_abpoa_uses_the_two_longest_and_widens_cells(self):
        a = poa_panel.predict_mb('abpoa', [12000, 6000, 10])
        self.assertAlmostEqual(a, 50 + 11 * 12000 * 6000 / 2 ** 20, places=0)
        b = poa_panel.predict_mb('abpoa', [20000, 6000])
        self.assertAlmostEqual(b, 50 + 22 * 20000 * 6000 / 2 ** 20, places=0)
        # a lone outlier costs far less than the hap32 arm's max_len^2 rule
        self.assertLess(poa_panel.predict_mb('abpoa', [40000, 5000]), 22 * 40000 ** 2 / 2 ** 20 / 5)


@unittest.skipUnless(HAVE, 'no regions/%s package' % RID)
class TestRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.td = tempfile.mkdtemp(prefix='test_poa_panel.')
        cls.cand, cls.pan, cls.work = (os.path.join(cls.td, x) for x in ('cand', 'panel', 'work'))
        cls.rows = {}
        for m in poa_panel.BASE_METHODS:
            cls.rows[m] = poa_panel.run_job(m, RID, cand_root=cls.cand, panel_root=cls.pan, workroot=cls.work,
                                            record_runtime=False)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.td, ignore_errors=True)

    def test_outputs_and_paths(self):
        rd = os.path.join(config.REGIONS_DIR, RID)
        hap = msa_graph.read_fasta(os.path.join(rd, 'hap32.fa'))
        hprc = [(n, s) for n, _, s in panel.read_fasta(os.path.join(rd, 'hprc.fa.gz'))]
        umap = panel.read_map(os.path.join(panel.UNION_DIR, RID + '.map.tsv'))
        for m, row in self.rows.items():
            self.assertEqual(row['status'], 'ok', row)
            self.assertEqual(row['panel'], 'all')
            self.assertEqual(row['n_distinct'], len(umap))
            o = poa_panel.out_paths(m, RID, self.cand, self.pan)
            for k in ('msa', 'gfa', 'json', 'full_msa', 'full_gfa'):
                self.assertTrue(os.path.exists(o[k]), (m, k))
            # the candidate graph spells hap32.fa, the full-panel graph spells hprc.fa.gz
            self.assertEqual(poa_panel.check_paths(o['gfa'], hap), len(hap))
            self.assertEqual(poa_panel.check_paths(o['full_gfa'], hprc), len(hprc))
            # projected MSA: hap32 names in hap32 order, no all-gap column
            pm = msa_graph.read_msa(o['msa'])
            self.assertEqual([n for n, _ in pm], [n for n, _ in hap])
            self.assertTrue(all(any(r[c] != '-' for _, r in pm) for c in range(len(pm[0][1]))))
            # full MSA: one row per distinct union sequence, each spelling it
            fm = dict(msa_graph.read_msa(o['full_msa']))
            useq = dict((n, s) for n, _, s in panel.read_fasta(os.path.join(panel.UNION_DIR, RID + '.fa')))
            self.assertEqual(set(fm), set(useq))
            self.assertTrue(all(fm[n].replace('-', '') == useq[n] for n in useq))
            j = json.load(open(o['json']))
            self.assertEqual(j['status'], 'ok')
            self.assertEqual(j['graph']['paths_verified'], len(hap))
            self.assertEqual(j['order'], 'longest')
            self.assertEqual(poa_panel.previous_status(m, RID, self.cand, self.pan), 'ok')

    def test_refused_by_prediction(self):
        row = poa_panel.run_job('poa_spoa', RID, mem_mb=10, cand_root=self.cand + '_x', panel_root=self.pan + '_x',
                                workroot=self.work, record_runtime=False)
        self.assertEqual(row['status'], 'memout')
        self.assertIn('predicted', row['note'])
        o = poa_panel.out_paths('poa_spoa', RID, self.cand + '_x', self.pan + '_x')
        self.assertFalse(os.path.exists(o['gfa']))
        self.assertEqual(json.load(open(o['json']))['status'], 'memout')

    def test_variant_order(self):
        row = poa_panel.run_job('poa_abpoa_random1', RID, base='poa_abpoa', order='random:1',
                                cand_root=self.cand, panel_root=self.pan, workroot=self.work, record_runtime=False)
        self.assertEqual(row['status'], 'ok')
        self.assertEqual(row['order'], 'random:1')

    def test_runtime_table(self):
        path = os.path.join(self.td, 'rt.tsv')
        poa_panel.update_runtime({'region_id': 'A', 'method': 'm', 'status': 'error'}, path)
        poa_panel.update_runtime({'region_id': 'A', 'method': 'm', 'status': 'ok'}, path)
        poa_panel.update_runtime({'region_id': 'B', 'method': 'm', 'status': 'ok'}, path)
        rows = list(csv.DictReader(open(path), delimiter='\t'))
        self.assertEqual([(r['region_id'], r['status']) for r in rows], [('A', 'ok'), ('B', 'ok')])
        self.assertEqual(list(rows[0])[:6], ['region_id', 'method', 'panel', 'n_distinct', 'seconds', 'status'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
