#!/usr/bin/env python3
"""Tests for tools/poa_panel_mc.py (the full-panel arm of poa_abpoa_mc).

    python3 tools/test_poa_panel_mc.py

1. The full-panel memory model and the plan: under BAR's 10 kb rule no call is longer than 10 kb,
   so no region can be predicted above the prediction for a 10 kb call.
2. What abPOA is given: every hprc.fa.gz record once (duplicates included) plus each hap32-only
   sequence once, grouped by distinct sequence in union order.
3. The region cap on the summed abPOA wall clock (realign_poa_mc.bar_align, region_timeout): the
   second end call gets what the first left, and nothing left is a timeout.
4. End to end on the smallest full panel (L014899: 457 records, 7 distinct) into a temporary
   directory: every output exists, the graphs spell hap32.fa and hprc.fa.gz, the full MSA has one
   row per record, and the two refusal paths (predicted memory, predicted MSA size) write their
   records without running abPOA. Needs regions/L014899, the union
   files (made on the fly) and Cactus's abPOA.
"""
import json
import os
import shutil
import sys
import tempfile
import time
import unittest

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402
import msa_graph  # noqa: E402
import panel  # noqa: E402
import poa_panel  # noqa: E402
import poa_panel_mc as ppm  # noqa: E402
import realign_poa_mc as mc  # noqa: E402

RID = 'L014899'
HAVE = os.path.exists(os.path.join(config.REGIONS_DIR, RID, 'hap32.fa'))
HAVE_ABPOA = os.path.exists(mc.ABPOA_MC)


class TestPredict(unittest.TestCase):
    def test_full_panel_constants_and_the_10kb_ceiling(self):
        lens = [mc.BANDING_LIMIT] * 400
        top = ppm.predict_call_mb(lens)
        alloc = 20.0 * (mc.BANDING_LIMIT + 4) * (ppm.PANEL_NODES_PER_BP * mc.BANDING_LIMIT + 2)
        p2 = 1 << int(alloc - 1).bit_length()
        self.assertAlmostEqual(top, 60 + ppm.PANEL_TOUCHED * p2 / 2 ** 20, places=0)
        # the hap32 constants are untouched by the full-panel ones
        self.assertEqual(mc.predict_call_mb([3000, 2000]), mc.predict_call_mb([3000, 2000], mc.NODES_PER_BP,
                                                                                mc.TOUCHED))
        self.assertLessEqual(top, ppm.DEFAULT_MEM_MB)

    @unittest.skipUnless(HAVE, 'no region packages')
    def test_plan_is_bounded_by_the_rule(self):
        top = ppm.predict_call_mb([mc.BANDING_LIMIT] * 2)
        for rid in ('L014899', 'L012184', 'L012272'):
            if not os.path.exists(os.path.join(config.REGIONS_DIR, rid, 'hap32.fa')):
                continue
            p = ppm.plan_region(rid)
            self.assertLessEqual(p['predicted_mb'], top)
            self.assertEqual(p['bar_rule'], 'single' if p['max_interior'] < mc.BANDING_LIMIT else 'ends')
            if p['bar_rule'] == 'ends':
                self.assertEqual(p['aligned_bp'], 2 * sum(min(mc.BANDING_LIMIT, x) for x in self._lens(rid, p)))

    @staticmethod
    def _lens(rid, p):
        import realign
        rd = os.path.join(config.REGIONS_DIR, rid)
        recs = ppm.panel_records(rd, p['union_map'])
        L, R = p['flanks']
        return [len(realign.mask_runs(s)[0]) - L - R for _, s, _ in recs]


@unittest.skipUnless(HAVE, 'no regions/%s package' % RID)
class TestRecords(unittest.TestCase):
    def test_every_record_once_and_hap32_only_once(self):
        rd = os.path.join(config.REGIONS_DIR, RID)
        fa, mp = poa_panel.ensure_union(rd, RID)
        recs = ppm.panel_records(rd, mp)
        hprc = [(n, s.upper()) for n, _, s in panel.read_fasta(os.path.join(rd, 'hprc.fa.gz'))]
        umap = panel.read_map(mp)
        names = [n for n, _, _ in recs]
        self.assertEqual(len(names), len(set(names)))
        self.assertEqual(len(recs), sum(max(1, int(r['weight'])) for r in umap))
        self.assertEqual(set(n for n, _ in hprc) - set(names), set())
        seq = dict((n, s) for n, s, _ in recs)
        for n, s in hprc:
            self.assertEqual(seq[n], s)
        # grouped by distinct id, in union order
        order = [u for _, _, u in recs]
        firsts = [u for k, u in enumerate(order) if k == 0 or order[k - 1] != u]
        self.assertEqual(firsts, [r['id'] for r in umap])
        # every hap32 sequence is represented
        hs = set(s for _, s, _ in recs)
        for n, s in msa_graph.read_fasta(os.path.join(rd, 'hap32.fa')):
            self.assertIn(s.upper(), hs)


class TestRegionCap(unittest.TestCase):
    def test_second_call_gets_what_is_left(self):
        got = []

        def slow(strings, tag=None, calls=None, timeout=None, **_):
            got.append((tag, timeout))
            time.sleep(0.4)
            w = max(len(x) for x in strings)
            return [x + '-' * (w - len(x)) for x in strings]      # a valid (left-flush) MSA
        L, R = 'AC', 'G'
        seqs = [L + 'A' * 30 + R, L + 'A' * 25 + R]
        old = mc.BANDING_LIMIT
        mc.BANDING_LIMIT = 10                    # force the two-end path
        try:
            with self.assertRaises(mc.AlignerFailed) as cm:
                mc.bar_align(seqs, (2, 1), poa=slow, timeout=100, region_timeout=0.3)
            self.assertEqual(cm.exception.status, 'timeout')
            self.assertEqual([t for t, _ in got], ['left'])
            got[:] = []
            rows, info = mc.bar_align(seqs, (2, 1), poa=slow, timeout=100, region_timeout=5)
            self.assertEqual([r.replace('-', '') for r in rows], seqs)
            self.assertEqual(info['region_timeout'], 5)
            self.assertEqual([t for t, _ in got], ['left', 'right'])
            self.assertAlmostEqual(got[0][1], 5, delta=0.1)
            self.assertLess(got[1][1], 4.8)
        finally:
            mc.BANDING_LIMIT = old


@unittest.skipUnless(HAVE and HAVE_ABPOA, 'no regions/%s package or no Cactus abPOA' % RID)
class TestRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.td = tempfile.mkdtemp(prefix='test_poa_panel_mc.')
        cls.cand, cls.pan, cls.work = (os.path.join(cls.td, x) for x in ('cand', 'panel', 'work'))
        cls.row = ppm.run_job(RID, cand_root=cls.cand, panel_root=cls.pan, workroot=cls.work, record_runtime=False)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.td, ignore_errors=True)

    def test_outputs_and_paths(self):
        rd = os.path.join(config.REGIONS_DIR, RID)
        hap = msa_graph.read_fasta(os.path.join(rd, 'hap32.fa'))
        hprc = [(n, s) for n, _, s in panel.read_fasta(os.path.join(rd, 'hprc.fa.gz'))]
        row = self.row
        self.assertEqual(row['status'], 'ok', row)
        self.assertEqual(row['method'], 'poa_abpoa_mc__all')
        self.assertEqual(row['bar_rule'], 'single')
        self.assertEqual(row['abpoa_calls'], 1)
        o = ppm.out_paths(RID, self.cand, self.pan)
        for k in ('msa', 'gfa', 'json', 'full_msa', 'full_gfa'):
            self.assertTrue(os.path.exists(o[k]), k)
        self.assertEqual(poa_panel.check_paths(o['gfa'], hap), len(hap))
        self.assertEqual(poa_panel.check_paths(o['full_gfa'], hprc), len(hprc))
        pm = msa_graph.read_msa(o['msa'])
        self.assertEqual([n for n, _ in pm], [n for n, _ in hap])
        self.assertTrue(all(any(r[c] != '-' for _, r in pm) for c in range(len(pm[0][1]))))
        fm = msa_graph.read_msa(o['full_msa'])
        self.assertEqual(len(fm), row['n_records'])
        want = dict((n, s) for n, s, _ in ppm.panel_records(rd, ppm.plan_region(RID)['union_map']))
        for n, r in fm:
            self.assertEqual(r.replace('-', ''), want[n])
        with open(o['json']) as f:
            j = json.load(f)
        self.assertEqual(j['status'], 'ok')
        self.assertIsInstance(j['dup_split'], int)
        c = j['align']['aligner']['calls'][0]
        self.assertEqual(c['n'], row['n_records'])
        self.assertIn('-p', c['command'])
        self.assertGreater(c['nodes'], 0)
        self.assertEqual(ppm.previous_status(RID, self.cand, self.pan), 'ok')

    def test_refused_by_prediction(self):
        cand = os.path.join(self.td, 'cand_refused')
        row = ppm.run_job(RID, mem_mb=1, cand_root=cand, panel_root=os.path.join(self.td, 'panel_refused'),
                          workroot=self.work, record_runtime=False)
        self.assertEqual(row['status'], 'memout')
        self.assertIn('predicted', row['note'])
        o = ppm.out_paths(RID, cand, os.path.join(self.td, 'panel_refused'))
        self.assertTrue(os.path.exists(o['json']))
        self.assertFalse(os.path.exists(o['gfa']))

    def test_skipped_when_the_msa_would_be_too_large(self):
        cand, pan = os.path.join(self.td, 'cand_big'), os.path.join(self.td, 'panel_big')
        old = ppm.MAX_MSA_CELLS
        ppm.MAX_MSA_CELLS = 10
        try:
            row = ppm.run_job(RID, cand_root=cand, panel_root=pan, workroot=self.work, record_runtime=False)
        finally:
            ppm.MAX_MSA_CELLS = old
        self.assertEqual(row['status'], 'skipped')
        self.assertIn('cannot hold', row['note'])
        self.assertEqual(ppm.previous_status(RID, cand, pan), 'skipped')


if __name__ == '__main__':
    unittest.main(verbosity=1)
