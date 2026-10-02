#!/usr/bin/env python3
"""Tests for tools/realign_units.py (M4, repeat-unit-aware realignment) and tools/unit_dp.c.

    python3 tools/test_realign_units.py            # a few seconds; mafft tests are skipped without mafft

Covers: the flank/cycle/flank DP (segmentation, partial units, free ends, bypass, local search),
the repeat model (period, consensus, rotation), the exact sum-of-pairs profile aligner against
brute force, the affine SP scorer, the progressive aligner and its refinement, the unit-column
expansion, the base-level polish, an end-to-end synthetic region (rows spell their input, units
stacked in phase, the graph builds), the fallback for a non-TR region, and mafft's text mode.
"""
import itertools
import json
import os
import random
import shutil
import sys
import tempfile
import unittest

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import realign_units as ru  # noqa: E402
import msa_graph  # noqa: E402
import config  # noqa: E402

HAVE_MAFFT = bool(shutil.which(config.MAFFT) or os.path.exists(config.MAFFT))


def brute_unit_cost(a, b, sc):
    """Linear-gap unit-level DP by plain Python (reference)."""
    n, m = len(a), len(b)
    D = [[0.0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        D[i][0] = D[i - 1][0] + sc.cost(a[i - 1], sc.G)
    for j in range(1, m + 1):
        D[0][j] = D[0][j - 1] + sc.cost(sc.G, b[j - 1])
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            D[i][j] = min(D[i - 1][j - 1] + sc.cost(a[i - 1], b[j - 1]),
                          D[i - 1][j] + sc.cost(a[i - 1], sc.G), D[i][j - 1] + sc.cost(sc.G, b[j - 1]))
    return D[n][m]


def pair_affine(ra, rb, sc, go):
    c, opens, st = 0.0, 0, 0
    for x, y in zip(ra, rb):
        if x < 0 and y < 0:
            continue
        if x < 0:
            c += sc.cost(sc.G, y)
            if st != 1:
                opens += 1
                st = 1
        elif y < 0:
            c += sc.cost(x, sc.G)
            if st != 2:
                opens += 1
                st = 2
        else:
            c += sc.cost(x, y)
            st = 0
    return c + go * opens


def mutate(s, rng, rate=0.05):
    out = []
    for ch in s:
        r = rng.random()
        if r < rate / 3:
            continue
        if r < 2 * rate / 3:
            out.append(rng.choice('ACGT'))
        elif r < rate:
            out.append(ch)
            out.append(rng.choice('ACGT'))
        else:
            out.append(ch)
    return ''.join(out)


class TestDP(unittest.TestCase):
    P = ru.Params()
    fl = 'GATTACAGATTACACCGGTTACGTAGGCTAGCTAAGG'
    fr = 'TTTTGGGGCCCCAAAATGCAGTCAGTCCAGTAGCTGA'
    u = 'ACGTTGCAT'

    def test_segmentation(self):
        arr = self.u * 3 + 'ACGTGCAT' + self.u + 'ACGTTGGCAT' + self.u[:5]
        seq = self.fl + arr + self.fr
        s = ru.segment(seq, self.fl, self.u, self.fr, self.P, False, False)
        self.assertEqual((s.a, s.b), (len(self.fl), len(self.fl) + len(arr)))
        self.assertEqual(''.join(x.seq for x in s.units), arr)
        self.assertEqual(len(s.units), 7)
        self.assertEqual(s.units[-1].seq, self.u[:5])            # partial last unit
        self.assertEqual(s.units[3].seq, 'ACGTGCAT')              # unit with a deletion
        for x in s.units:                                          # slots non-decreasing inside a unit
            self.assertEqual(list(x.slots), sorted(x.slots))

    def test_free_ends_and_bypass(self):
        s = ru.segment('AAAAAAAAAA' + self.fl + self.u * 2 + self.fr + 'CCCCCC', self.fl, self.u, self.fr,
                       self.P, True, True)
        self.assertEqual(s.left, 'AAAAAAAAAA' + self.fl)
        self.assertEqual([x.seq for x in s.units], [self.u, self.u])
        s = ru.segment(self.fl + self.fr, self.fl, self.u, self.fr, self.P, False, False)
        self.assertEqual((s.a, s.b, len(s.units)), (len(self.fl), len(self.fl), 0))

    def test_local_free_entry(self):
        seq = 'CCCCCCTTTT' + self.u[3:] + self.u * 4 + 'GAGAGAGAGA'
        s = ru.segment(seq, '', self.u, '', self.P, True, True, free_entry=True)
        self.assertEqual(s.a, 10)
        self.assertEqual(s.units[0].events[0][:2], ('M', 3))       # entered mid-unit
        self.assertEqual([x.seq for x in s.units[1:]], [self.u] * 4)

    def test_model_rotation(self):
        rng = random.Random(3)
        motif = 'GGTGACCTAGGACAT'
        flank_l = ''.join(rng.choice('ACGT') for _ in range(300))
        flank_r = ''.join(rng.choice('ACGT') for _ in range(300))
        arr = ''.join(mutate(motif, rng, 0.03) for _ in range(40))
        chm = flank_l + arr + flank_r
        rj = {'stratum': 'hotspot_vntr', 'span_start': 1001, 'core_start': 1001 + 300, 'core_end': 1000 + 300 + len(arr),
              'period': 15, 'motif': motif[5:] + motif[:5], 'tr_annotation': {}}
        m, why = ru.derive_model(chm, rj, self.P)
        self.assertIsNone(why)
        self.assertIsNone(m.reject)
        self.assertEqual(m.period, 15)
        self.assertIn(m.motif, [motif[k:] + motif[:k] for k in range(15)])
        self.assertLessEqual(abs(m.a - 300), 3)
        self.assertGreater(m.identity, 0.9)
        # CHM13's own array starts at unit position 0
        s = ru.segment(chm, chm[max(0, m.a - 300):m.a], m.motif, chm[m.b:m.b + 300], self.P, m.a > 300,
                       len(chm) - m.b > 300)
        self.assertEqual(s.units[0].events[0][1], 0)


class TestProfileAligner(unittest.TestCase):
    def setUp(self):
        rng = random.Random(7)
        base = 'ACGTAGGCTA'
        self.reps = list(dict.fromkeys([mutate(base, rng, 0.3) or 'A' for _ in range(12)]))
        self.sc = ru.SymCost(self.reps)
        self.rng = rng

    def rand_string(self, n):
        return tuple(self.rng.randrange(len(self.reps)) for _ in range(n))

    def test_pairwise_optimal(self):
        for _ in range(20):
            a, b = self.rand_string(self.rng.randint(1, 9)), self.rand_string(self.rng.randint(1, 9))
            A, LA = ru.array('i', a), len(a)
            B, LB = ru.array('i', b), len(b)
            merged, L, cost = ru.align_flat(A, 1, LA, [1.0], B, 1, LB, [1.0], self.sc)
            self.assertAlmostEqual(cost, brute_unit_cost(a, b, self.sc))
            rows = [list(merged[i * L:(i + 1) * L]) for i in range(2)]
            self.assertEqual(tuple(x for x in rows[0] if x >= 0), a)
            self.assertEqual(tuple(x for x in rows[1] if x >= 0), b)
            self.assertAlmostEqual(ru.sp_total(merged, 2, L, [1, 1], self.sc), cost)

    def test_profile_optimal_brute_force(self):
        """Merging a 2-row profile with a row is optimal over all column-preserving merges."""
        for _ in range(6):
            a, b = self.rand_string(4), self.rand_string(3)
            m2, L2, _ = ru.align_flat(ru.array('i', a), 1, 4, [1.0], ru.array('i', b), 1, 3, [1.0], self.sc)
            c = self.rand_string(self.rng.randint(1, 4))
            merged, L, cost = ru.align_flat(m2, 2, L2, [1.0, 1.0], ru.array('i', c), 1, len(c), [1.0], self.sc)
            prof = [list(m2[i * L2:(i + 1) * L2]) for i in range(2)]
            best = None
            for Ltot in range(max(L2, len(c)), L2 + len(c) + 1):
                for cols in itertools.combinations(range(Ltot), L2):      # where the profile columns go
                    free = [k for k in range(Ltot) if k not in cols]
                    if len(free) > len(c):
                        continue
                    for rpos in itertools.combinations(range(Ltot), len(c)):
                        if not set(free) <= set(rpos):
                            continue
                        rows = [[-1] * Ltot for _ in range(3)]
                        for k, col in enumerate(cols):
                            rows[0][col], rows[1][col] = prof[0][k], prof[1][k]
                        for k, col in enumerate(rpos):
                            rows[2][col] = c[k]
                        v = sum(pair_affine(rows[x], rows[2], self.sc, 0.0) for x in range(2))
                        best = v if best is None else min(best, v)
            self.assertAlmostEqual(cost, best)

    def test_sp_scorers(self):
        for go in (0.0, 1.5):
            rows = []
            for _ in range(5):
                rows.append([self.rng.randrange(len(self.reps)) if self.rng.random() < 0.7 else -1 for _ in range(15)])
            flat = ru.array('i', [x for r in rows for x in r])
            w = [1.0, 2.0, 1.0, 1.0, 3.0]
            ref = sum(w[a] * w[b] * pair_affine(rows[a], rows[b], self.sc, go)
                      for a, b in itertools.combinations(range(5), 2))
            self.assertAlmostEqual(ru.sp_total(flat, 5, 15, w, self.sc, go), ref)
            if go == 0.0:
                self.assertAlmostEqual(ru.sp_cost_py(rows, w, self.sc), ref)
            cr = sum(w[a] * w[b] * pair_affine(rows[a], rows[b], self.sc, go) for a in (0, 1) for b in (2, 3, 4))
            self.assertAlmostEqual(ru.cross_total(flat, 15, [0, 1], [2, 3, 4], w, self.sc, go), cr)

    def test_progressive_and_refinement(self):
        P = ru.Params()
        strings = [self.rand_string(self.rng.randint(5, 14)) for _ in range(9)]
        strings = list(dict.fromkeys(strings))
        w = [1.0] * len(strings)
        rows, st = ru.prog_unit_msa(strings, w, self.sc, P)
        for r, s in zip(rows, strings):
            self.assertEqual(tuple(x for x in r if x >= 0), s)
        self.assertLessEqual(st['sp_final'], st['sp_progressive'] + 1e-6)
        lower = sum(brute_unit_cost(a, b, self.sc) for a, b in itertools.combinations(strings, 2))
        self.assertGreaterEqual(st['sp_final'], lower - 1e-6)     # never below the pairwise optima
        P2 = ru.Params(unit_go=3.0)
        rows2, st2 = ru.prog_unit_msa(strings, w, self.sc, P2)
        for r, s in zip(rows2, strings):
            self.assertEqual(tuple(x for x in r if x >= 0), s)


class TestExpansionAndPolish(unittest.TestCase):
    def test_expand_unit_column(self):
        P = ru.Params()
        motif = 'ACGTTGCAT'
        toks = ['ACGTTGCAT', 'ACGTGCAT', 'ACGTTGGCAT', 'TTGCAT', 'ACGTTGGGCAT']
        slots = {}
        for t in toks:
            s = ru.segment(t, '', motif, '', P, False, False)      # global, as the region DP is
            self.assertEqual(len(s.units), 1)
            self.assertEqual(s.units[0].seq, t)
            slots[t] = tuple(s.units[0].slots)
        out = ru.expand_unit_column(slots, len(motif))
        self.assertEqual(len(set(len(v) for v in out.values())), 1)
        for t, v in out.items():
            self.assertEqual(v.replace('-', ''), t)

    def test_star_msa(self):
        al = ru.star_msa(['GGA', 'GGAT', 'GA', 'GGA'])
        self.assertEqual(len(set(len(v) for v in al.values())), 1)
        for s, v in al.items():
            self.assertEqual(v.replace('-', ''), s)

    def test_polish_fixes_boundary_shift(self):
        # a SNP next to a unit boundary written as a deletion plus an insertion
        rows = ['GTTCTCC-TTTCATC', 'GTTCTC-ATTTCATC', 'GTTCTCC-TTTCATC']
        P = ru.Params()
        new, st = ru.polish_msa(rows, P, 10)
        for a, b in zip(rows, new):
            self.assertEqual(a.replace('-', ''), b.replace('-', ''))
        self.assertLess(st['sp_after'], st['sp_before'])
        self.assertEqual(st['sp_after'], 2.0)     # one mismatch against each of the two other rows


class TestEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='test_realign_units.')
        rng = random.Random(11)
        motif = 'GGAGGGTGAGGTGTGTA'
        variants = [motif, mutate(motif, rng, 0.12), mutate(motif, rng, 0.12), motif[:10] + 'C' + motif[11:]]
        left = ''.join(rng.choice('ACGT') for _ in range(250))
        right = ''.join(rng.choice('ACGT') for _ in range(250))
        base = [rng.randrange(len(variants)) for _ in range(30)]
        chm_arr = ''.join(variants[k] for k in base)
        recs = [('CHM13#0#chrT', left + chm_arr + right)]
        for h in range(10):
            units = list(base)
            if h % 3 == 1:
                del units[5:5 + h]
            if h % 3 == 2:
                units[12:12] = units[3:3 + h]
            arr = ''.join(mutate(variants[k], rng, 0.01) for k in units)
            l2 = left if h != 4 else left[:100] + 'A' + left[101:]
            recs.append(('sample#%d#chrT' % h, l2 + arr + right))
        s = recs[3][1]
        recs[3] = (recs[3][0], s[:400] + 'N' * 25 + s[400:])       # an assembly gap inside the array
        cls.recs = recs
        cls.fa = os.path.join(cls.tmp, 'in.fa')
        with open(cls.fa, 'w') as f:
            for n, s in recs:
                f.write('>%s\n%s\n' % (n, s))
        cls.rj = {'region_id': 'TEST', 'stratum': 'hotspot_vntr', 'span_start': 1, 'core_start': 251,
                  'core_end': 250 + len(chm_arr), 'period': 17, 'motif': motif, 'tr_annotation': {}}
        cls.rjp = os.path.join(cls.tmp, 'region.json')
        json.dump(cls.rj, open(cls.rjp, 'w'))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def run_one(self, **kw):
        out = os.path.join(self.tmp, 'out%d.msa.fa' % len(kw))
        info = ru.align_units_fasta(self.fa, self.rjp, out, threads=1, params=ru.Params(**kw), workdir=self.tmp)
        return out, info

    @unittest.skipUnless(HAVE_MAFFT, 'mafft not installed (flank pieces use realign.py)')
    def test_unit_path(self):
        out, info = self.run_one()
        self.assertEqual(info['status'], 'ok', info.get('reason'))
        self.assertEqual(info['mode'], 'unit')
        self.assertEqual(info['model']['period'], 17)
        rows = msa_graph.read_msa(out)
        self.assertEqual([n for n, _ in rows], [n for n, _ in self.recs])
        for (n, r), (_, s) in zip(rows, self.recs):
            self.assertEqual(r.replace('-', ''), s.upper())
        self.assertEqual(info['masked_runs'], 1)
        self.assertEqual(info['phase_check']['units_in_one_start_columns'], 1.0)
        hap = os.path.join(self.tmp, 'hap.fa')
        with open(hap, 'w') as f:
            for n, s in self.recs:
                f.write('>%s\n%s\n' % (n, s))
        st = msa_graph.msa_to_gfa(out, hap, os.path.join(self.tmp, 'g.gfa'), engine='native')
        self.assertEqual(st['paths_checked'], len(self.recs))
        self.assertGreaterEqual(st['N_bp_in_nodes'], 25)

    @unittest.skipUnless(HAVE_MAFFT, 'mafft not installed')
    def test_mafft_text_mode(self):
        out, info = self.run_one(unit_aligner='mafft', polish_rounds=0)
        self.assertEqual(info['status'], 'ok', info.get('reason'))
        self.assertEqual(info['mode'], 'unit', info.get('reason'))
        self.assertEqual(info['units']['unit_msa']['mode'], 'linsi')

    @unittest.skipUnless(HAVE_MAFFT, 'mafft not installed')
    def test_fallback_non_tr(self):
        rj = dict(self.rj, stratum='control_nontr_sv')
        p = os.path.join(self.tmp, 'nontr.json')
        json.dump(rj, open(p, 'w'))
        out = os.path.join(self.tmp, 'fb.msa.fa')
        info = ru.align_units_fasta(self.fa, p, out, threads=1, workdir=self.tmp)
        self.assertEqual(info['mode'], 'fallback')
        self.assertEqual(info['status'], 'ok')
        self.assertIn('not a tandem repeat', info['reason'])
        for (n, r), (_, s) in zip(msa_graph.read_msa(out), self.recs):
            self.assertEqual(r.replace('-', ''), s.upper())

    def test_fallback_engine_abpoa(self):
        P = ru.Params(fallback_engine='abpoa')
        self.assertEqual((P.flank_method, P.fallback, P.fallback_large), ('poa_abpoa',) * 3)
        self.assertEqual(ru.Params(fallback_engine='abpoa', fallback='mafft_linsi').fallback, 'mafft_linsi')
        self.assertEqual(ru.Params().fallback, ru.DEFAULT_FALLBACK)
        self.assertRaises(ValueError, ru.Params, fallback_engine='spoa')
        import realign_poa
        if not shutil.which(realign_poa.ABPOA):
            self.skipTest('abpoa not installed')
        rj = dict(self.rj, stratum='control_nontr_sv')
        p = os.path.join(self.tmp, 'nontr_poa.json')
        json.dump(rj, open(p, 'w'))
        out = os.path.join(self.tmp, 'fb_poa.msa.fa')
        info = ru.align_units_fasta(self.fa, p, out, threads=1, params=P, workdir=self.tmp)
        self.assertEqual((info['status'], info['mode'], info['fallback']['method']), ('ok', 'fallback', 'poa_abpoa'))
        for (n, r), (_, s) in zip(msa_graph.read_msa(out), self.recs):
            self.assertEqual(r.replace('-', ''), s.upper())


if __name__ == '__main__':
    unittest.main()
