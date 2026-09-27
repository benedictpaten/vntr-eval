#!/usr/bin/env python3
"""Tests for tools/mafft_panel.py (the full-panel arm of the mafft realigners).

    python3 tools/test_mafft_panel.py

Runs mafft FFT-NS-2 and L-INS-i end to end on the smallest packaged region (L014899: 7 distinct
full-panel sequences) into a temporary directory, plus the cost and memory models, the skip rule,
a not-run record and the runtime table. Needs regions/L014899, work/panel/union (made on the fly
if missing) and mafft.
"""
import argparse
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
import mafft_panel  # noqa: E402

RID = 'L014899'
HAVE = os.path.exists(os.path.join(config.REGIONS_DIR, RID, 'hap32.fa'))


def fake_info(lens):
    t = sum(lens)
    return {'mlens': lens, 'pair_cells': (t * t - sum(x * x for x in lens)) // 2}


class TestModels(unittest.TestCase):
    def test_pair_cells(self):
        self.assertEqual(fake_info([3, 5, 7])['pair_cells'], 3 * 5 + 3 * 7 + 5 * 7)

    def test_pairwise_time_scales_with_pair_cells_over_threads(self):
        i = fake_info([1000] * 10)
        self.assertAlmostEqual(mafft_panel.predict_s('mafft_linsi', i, 1),
                               mafft_panel.C_PAIRWISE * 45 * 1e6, places=6)
        self.assertAlmostEqual(mafft_panel.predict_s('mafft_linsi', i, 3) * 3,
                               mafft_panel.predict_s('mafft_linsi', i, 1), places=6)

    def test_pairwise_drops_to_one_thread_when_the_dp_does_not_fit(self):
        a = argparse.Namespace(threads_pairwise=3, threads_progressive=2, mem_mb=12000)
        t, pm = mafft_panel.plan_threads('mafft_einsi', fake_info([30000, 30000, 100]), a)
        self.assertEqual(t, 1)
        self.assertAlmostEqual(pm, 400 + 20.0 * 30000 * 30000 / 2 ** 20, places=3)
        t, _ = mafft_panel.plan_threads('mafft_einsi', fake_info([3000, 3000]), a)
        self.assertEqual(t, 3)
        t, _ = mafft_panel.plan_threads('mafft_fftnsi', fake_info([30000, 30000]), a)
        self.assertEqual(t, 2)

    def test_skip_rule(self):
        a = argparse.Namespace(skip_after=2, skip_factor=3.0, timeout=100.0, no_skip={'mafft_fftns2'})
        job = {'method': 'mafft_linsi', 'cost': 50, 'pred_s': 10.0}
        self.assertIsNone(mafft_panel.skip_reason(job, {'mafft_linsi': [(10, 'A')]}, a))
        self.assertIn('not attempted', mafft_panel.skip_reason(job, {'mafft_linsi': [(10, 'A'), (40, 'B')]}, a))
        # a timeout on a larger region does not count, nor one of another method
        self.assertIsNone(mafft_panel.skip_reason(job, {'mafft_linsi': [(10, 'A'), (60, 'B')],
                                                        'mafft_einsi': [(1, 'C'), (2, 'D')]}, a))
        # the progressive iterative method is skipped the same way (its own cost measure)
        prog = {'method': 'mafft_fftnsi', 'cost': 50, 'pred_s': 10.0}
        self.assertIn('N x L1 x L2', mafft_panel.skip_reason(prog, {'mafft_fftnsi': [(1, 'A'), (2, 'B')]}, a))
        prog['pred_s'] = 301.0
        self.assertIn('predicted', mafft_panel.skip_reason(prog, {}, a))
        # a no-skip method is always attempted
        fft = {'method': 'mafft_fftns2', 'cost': 50, 'pred_s': 1e6}
        self.assertIsNone(mafft_panel.skip_reason(fft, {'mafft_fftns2': [(1, 'A'), (2, 'B')]}, a))

    def test_cost_measure(self):
        i = fake_info([10, 8, 3])
        self.assertEqual(mafft_panel.cost_of('mafft_ginsi', i), i['pair_cells'])
        self.assertEqual(mafft_panel.cost_of('mafft_fftnsi', i), 3 * 10 * 8)

    def test_priority(self):
        self.assertEqual(mafft_panel.parse_priority('a=0.5, b=2'), {'a': 0.5, 'b': 2.0})


@unittest.skipUnless(HAVE, 'no regions/%s package' % RID)
class TestRun(unittest.TestCase):
    METHODS = ('mafft_fftns2', 'mafft_linsi')

    @classmethod
    def setUpClass(cls):
        cls.td = tempfile.mkdtemp(prefix='test_mafft_panel.')
        cls.cand, cls.pan, cls.work = (os.path.join(cls.td, x) for x in ('cand', 'panel', 'work'))
        cls.rows = {}
        for m in cls.METHODS:
            cls.rows[m] = mafft_panel.run_job(m, RID, threads=1, cand_root=cls.cand, panel_root=cls.pan,
                                              workroot=cls.work, record_runtime=False)

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
            o = mafft_panel.out_paths(m, RID, self.cand, self.pan)
            for k in ('msa', 'gfa', 'json', 'full_msa', 'full_gfa'):
                self.assertTrue(os.path.exists(o[k]), (m, k))
            self.assertEqual(mafft_panel.check_paths(o['gfa'], hap), len(hap))
            self.assertEqual(mafft_panel.check_paths(o['full_gfa'], hprc), len(hprc))
            pm = msa_graph.read_msa(o['msa'])
            self.assertEqual([n for n, _ in pm], [n for n, _ in hap])
            self.assertTrue(all(any(r[c] != '-' for _, r in pm) for c in range(len(pm[0][1]))))
            fm = dict(msa_graph.read_msa(o['full_msa']))
            useq = dict((n, s) for n, _, s in panel.read_fasta(os.path.join(panel.UNION_DIR, RID + '.fa')))
            self.assertEqual(set(fm), set(useq))
            self.assertTrue(all(fm[n].replace('-', '') == useq[n] for n in useq))
            with open(o['json']) as f:
                j = json.load(f)
            self.assertEqual(j['status'], 'ok')
            self.assertEqual(j['graph']['paths_verified'], len(hap))
            self.assertEqual(j['panel_graph']['paths_verified'], len(hprc))
            self.assertIn('mafft', ' '.join(j['align']['aligner']['command']))
            self.assertEqual(mafft_panel.previous_status(m, RID, self.cand, self.pan), 'ok')

    def test_not_run_record(self):
        a = argparse.Namespace(cand_root=self.cand + '_x', panel_root=self.pan + '_x', mem_mb=1, timeout=1)
        row = mafft_panel.record_not_run('mafft_einsi', RID, 'skipped', 'not attempted: test', 5.0, 10, 3, a,
                                         record_runtime=False)
        self.assertEqual(row['status'], 'skipped')
        o = mafft_panel.out_paths('mafft_einsi', RID, a.cand_root, a.panel_root)
        self.assertFalse(os.path.exists(o['gfa']))
        with open(o['json']) as f:
            self.assertEqual(json.load(f)['status'], 'skipped')
        self.assertEqual(mafft_panel.previous_status('mafft_einsi', RID, a.cand_root, a.panel_root), 'skipped')

    def test_runtime_table(self):
        path = os.path.join(self.td, 'rt.tsv')
        mafft_panel.update_runtime({'region_id': 'A', 'method': 'mafft_linsi', 'status': 'error'}, path)
        mafft_panel.update_runtime({'region_id': 'A', 'method': 'mafft_linsi', 'status': 'ok'}, path)
        mafft_panel.update_runtime({'region_id': 'B', 'method': 'mafft_fftns2', 'status': 'ok'}, path)
        with open(path) as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual([(r['method'], r['region_id'], r['status']) for r in rows],
                         [('mafft_fftns2', 'B', 'ok'), ('mafft_linsi', 'A', 'ok')])
        self.assertEqual(list(rows[0])[:6], ['region_id', 'method', 'panel', 'n_distinct', 'seconds', 'status'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
