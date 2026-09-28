#!/usr/bin/env python3
"""Offline tests for score_haplotypes.py (no bcftools, truvari or big data needed).

    python3 tools/test_score_haplotypes.py
"""
import itertools
import json
import os
import random
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score_haplotypes as sc  # noqa: E402

REF = 'ACGTACGTTTGCAGGCATCCGATTACAGGATCCATGACCTGAAGTCCATGCATGCAAC'   # 58 bp
A1 = 1001
B1 = A1 + len(REF) - 1


def ed(a, b):
    return sc.edit_distance(a, b)


def apply(recs):
    """Reference implementation for the tests: apply (pos, ref, alt) records, rightmost first."""
    s = REF
    for pos, ref, alt in sorted(recs, key=lambda r: -r[0]):
        i = pos - A1
        assert s[i:i + len(ref)] == ref, (pos, ref, s[i:i + len(ref)])
        s = s[:i] + alt + s[i + len(ref):]
    return s


def vcf_text(records, contig='chr1', sample='HG002', header_contigs=('chr1', 'chrX')):
    """records: (pos, ref, alts, gt[, filter[, ps[, id[, info]]]])."""
    lines = ['##fileformat=VCFv4.2'] + ['##contig=<ID=%s>' % c for c in header_contigs] + [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample]
    for r in records:
        pos, ref, alts, gt = r[:4]
        filt = r[4] if len(r) > 4 and r[4] else 'PASS'
        ps = r[5] if len(r) > 5 else None
        rid = r[6] if len(r) > 6 else '.'          # vg snarl ID, e.g. '>1>20'
        info = r[7] if len(r) > 7 else '.'         # e.g. 'AT=>1>2>20,>1>20'
        fmt, val = ('GT:PS', '%s:%s' % (gt, ps)) if ps is not None else ('GT', gt)
        lines.append('\t'.join([contig, str(pos), rid, ref, ','.join(alts), '60', filt, info, fmt, val]))
    return '\n'.join(lines) + '\n'


class Region:
    """A synthetic region package in a temporary directory."""

    def __init__(self, t1, t2, contig='chr1', a1=A1, ref=REF, rid='T000001', stratum='hotspot_vntr'):
        self.dir = tempfile.mkdtemp(prefix='sh_test_')
        self.contig = contig
        rj = {'region_id': rid, 'stratum': stratum, 'contig': contig, 'span_start': a1,
              'span_end': a1 + len(ref) - 1, 'core_start': a1 + 5, 'core_end': a1 + len(ref) - 6}
        with open(os.path.join(self.dir, 'region.json'), 'w') as f:
            json.dump(rj, f)
        with open(os.path.join(self.dir, 'hap32.fa'), 'w') as f:
            f.write('>CHM13#0#%s len=%d\n%s\n>other#1#%s#0\n%s\n' % (contig, len(ref), ref, contig, ref))
        with open(os.path.join(self.dir, 'truth.fa'), 'w') as f:
            f.write('>HG002#1 x\n%s\n>HG002#2 y\n%s\n' % (t1, t2))

    def vcf(self, records, **kw):
        p = os.path.join(self.dir, 'calls.vcf')
        with open(p, 'w') as f:
            f.write(vcf_text(records, contig=kw.pop('contig', self.contig), **kw))
        return p

    def score(self, records, **kw):
        vkw = {k: kw.pop(k) for k in ('sample', 'header_contigs', 'contig') if k in kw}
        return sc.score_region(self.dir, self.vcf(records, **vkw), label='test', truvari=False, **kw)

    def close(self):
        shutil.rmtree(self.dir)


# a small truth: a SNV (hom), a 5 bp deletion (het, h1), a 4 bp insertion (het, h2)
SNV = (1003, 'G', 'T')
DEL = (1010, 'TGCAGG', 'T')
INS = (1028, 'G', 'GAAAA')
T1 = apply([SNV, DEL])
T2 = apply([SNV, INS])
TRUTH_RECS = [(1003, 'G', ['T'], '1|1'), (1010, 'TGCAGG', ['T'], '1|0'), (1028, 'G', ['GAAAA'], '0|1')]


class TestTrimAndEdit(unittest.TestCase):
    def test_trim(self):
        self.assertEqual(sc.trim_allele(10, 'ACGT', 'A'), (11, 14, ''))
        self.assertEqual(sc.trim_allele(10, 'A', 'ACC'), (11, 11, 'CC'))
        self.assertEqual(sc.trim_allele(10, 'G', 'T'), (10, 11, 'T'))
        self.assertEqual(sc.trim_allele(10, 'CAT', 'CGT'), (11, 12, 'G'))

    def test_inside_is_raw(self):
        r = sc.parse_record('chr1\t1010\t.\tTGCAGG\tT\t.\tPASS\t.\tGT\t0|1', 0, 0)
        why, info, crossing = sc.allele_edit(r, 1, A1, B1, REF)
        self.assertEqual((why, crossing), ('ok', False))
        self.assertEqual({k: info[k] for k in ('raw', 'trim', 'apply', 'tseq')},
                         {'raw': (1010, 1015), 'trim': (1011, 1016), 'apply': (1011, 1016, ''), 'tseq': ''})
        self.assertEqual(sc.apply_edits(REF, A1, [info['apply']]), apply([(1010, 'TGCAGG', 'T')]))

    def test_left_crossing_deletion_is_clipped(self):
        # deletion from 996 (outside) to 1005: removes the span's first 5 bases
        r = {'pos': 995, 'ref': 'N' * 6 + REF[:5], 'alts': ['N'], 'idx': 0}
        r['ref'] = 'NNNNNN' + REF[:5]
        why, e, crossing = sc.allele_edit(r, 1, A1, B1, REF)
        self.assertEqual(why, 'ok')
        self.assertTrue(crossing)
        self.assertEqual(e['apply'], (A1, A1 + 5, ''))
        self.assertEqual(sc.apply_edits(REF, A1, [e['apply']]), REF[5:])

    def test_right_crossing_substitution_is_clipped(self):
        # REF = the span's last 2 bases + 2 outside, ALT = 2 bases: the 2 REF bases outside the span
        # consume the 2 trailing ALT bases, so inside the span the record is a 2 bp deletion
        r = {'pos': B1 - 1, 'ref': REF[-2:] + 'AA', 'alts': ['GG'], 'idx': 0}
        why, e, crossing = sc.allele_edit(r, 1, A1, B1, REF)
        self.assertEqual((why, crossing), ('ok', True))
        self.assertEqual(sc.apply_edits(REF, A1, [e['apply']]), REF[:-2])
        r['alts'] = ['GGTTCC']                     # 6 bp ALT: the inner 4 stay
        why, e, crossing = sc.allele_edit(r, 1, A1, B1, REF)
        self.assertEqual(sc.apply_edits(REF, A1, [e['apply']]), REF[:-2] + 'GGTT')

    def test_edge_insertion_dropped(self):
        # a boundary-crossing record that trims to an insertion after the span's last base
        r = {'pos': B1, 'ref': REF[-1] + 'A', 'alts': [REF[-1] + 'TTTA'], 'idx': 0}
        self.assertEqual(sc.allele_edit(r, 1, A1, B1, REF)[0], 'edge_insertion')
        # a record wholly inside the span by REF is applied as written, even an insertion after B1
        r = {'pos': B1, 'ref': REF[-1], 'alts': [REF[-1] + 'TTT'], 'idx': 0}
        why, e, crossing = sc.allele_edit(r, 1, A1, B1, REF)
        self.assertEqual((why, crossing), ('ok', False))
        self.assertEqual(sc.apply_edits(REF, A1, [e['apply']]), REF + 'TTT')

    def test_padding_only_overlap_is_outside(self):
        # a deletion whose REF touches the span only through its padding base
        r = {'pos': B1, 'ref': REF[-1] + 'ACGT', 'alts': [REF[-1]], 'idx': 0}
        self.assertEqual(sc.allele_edit(r, 1, A1, B1, REF)[0], 'outside_after_trim')

    def test_ref_mismatch_and_symbolic(self):
        r = {'pos': 1003, 'ref': 'A', 'alts': ['T', '<DEL>', '*'], 'idx': 0}
        self.assertEqual(sc.allele_edit(r, 1, A1, B1, REF)[0], 'ref_mismatch')
        r['ref'] = 'G'
        self.assertEqual(sc.allele_edit(r, 2, A1, B1, REF)[0], 'symbolic')
        self.assertEqual(sc.allele_edit(r, 3, A1, B1, REF)[0], 'star')


class TestPrimary(unittest.TestCase):
    def setUp(self):
        self.reg = Region(T1, T2)

    def tearDown(self):
        self.reg.close()

    def test_truth_gives_zero(self):
        r = self.reg.score(TRUTH_RECS)
        self.assertEqual(r['status'], 'ok')
        self.assertEqual(r['ed'], 0)
        self.assertTrue(r['exact'])
        self.assertEqual(r['ed_ref'], ed(REF, T1) + ed(REF, T2))
        self.assertEqual(r['gain'], 1.0)
        self.assertEqual(r['called_len'], [len(T1), len(T2)])
        self.assertEqual(r['phase']['method'], 'fully_phased')

    def test_empty_gives_ed_ref(self):
        r = self.reg.score([])
        self.assertEqual(r['ed'], r['ed_ref'])
        self.assertEqual(r['gain'], 0.0)

    def test_pairing_is_symmetric(self):
        swapped = [(p, rf, a, g[::-1]) for p, rf, a, g in TRUTH_RECS]
        r = self.reg.score(swapped)
        self.assertEqual(r['ed'], 0)
        self.assertEqual(r['pairing'], 'slot1-h2,slot2-h1')

    def test_dropped_record(self):
        r = self.reg.score(TRUTH_RECS[:2])
        self.assertEqual(r['ed'], 4)
        self.assertEqual(r['ed_h2'], 4)

    def test_phase_switch_is_detected(self):
        recs = [TRUTH_RECS[0], TRUTH_RECS[1], (1028, 'G', ['GAAAA'], '1|0')]
        r = self.reg.score(recs)
        self.assertGreater(r['ed'], 0)
        self.assertEqual(r['phase']['n_blocks'], 1)

    def test_unphased_record_is_rephased(self):
        recs = [TRUTH_RECS[0], TRUTH_RECS[1], (1028, 'G', ['GAAAA'], '1/0')]
        r = self.reg.score(recs)
        self.assertEqual(r['ed'], 0)
        self.assertEqual(r['phase']['method'], 'enumerated')
        self.assertGreater(r['phase']['ed_given_phase'], 0)
        self.assertIn('unphased_het_records', r['flags'])

    def test_ps_blocks_are_rephased(self):
        recs = [(1003, 'G', ['T'], '1|1', None, 5), (1010, 'TGCAGG', ['T'], '1|0', None, 5),
                (1028, 'G', ['GAAAA'], '1|0', None, 9)]
        r = self.reg.score(recs)
        self.assertEqual(r['ed'], 0)
        self.assertEqual(r['phase']['n_phase_sets'], 2)
        self.assertEqual(r['phase']['n_blocks'], 2)
        # the same records in one PS: the switch is real and is charged
        recs[2] = (1028, 'G', ['GAAAA'], '1|0', None, 5)
        self.assertGreater(self.reg.score(recs)['ed'], 0)

    def test_greedy_branch(self):
        recs = [(1003, 'G', ['T'], '1|1', None, 5), (1010, 'TGCAGG', ['T'], '1|0', None, 5),
                (1028, 'G', ['GAAAA'], '1|0', None, 9)]
        r = self.reg.score(recs, enum_cap=0)
        self.assertEqual(r['phase']['method'], 'greedy')
        self.assertTrue(r['phase']['enum_cap_hit'])
        self.assertIn('phase_greedy', r['flags'])
        self.assertEqual(r['ed'], 0)

    def test_missing_and_filtered_are_reference(self):
        recs = [TRUTH_RECS[0], (1010, 'TGCAGG', ['T'], '1|.'), (1028, 'G', ['GAAAA'], '0|1', 'lowconf')]
        r = self.reg.score(recs)
        self.assertEqual(r['records']['missing_alleles_as_ref'], 1)
        self.assertEqual(r['records']['filtered_records_as_ref'], 1)
        self.assertEqual(r['ed'], 4)             # the filtered insertion is missing from h2
        r2 = self.reg.score(recs, apply_filtered=True)
        self.assertEqual(r2['ed'], 0)
        self.assertEqual(r2['records']['filtered_records_applied'], 1)

    def test_multiallelic(self):
        t1 = apply([(1003, 'G', 'T')])
        t2 = apply([(1003, 'G', 'C')])
        reg = Region(t1, t2)
        try:
            r = reg.score([(1003, 'G', ['T', 'C'], '1|2')])
            self.assertEqual(r['ed'], 0)
            r = reg.score([(1003, 'G', ['T', 'C'], '1|1')])
            self.assertEqual(r['ed'], 1)
        finally:
            reg.close()

    def test_nested_overlap_rule(self):
        """A parent record (long REF) wins in the slot where it is non-reference; its child fills in
        the slot where the parent is reference; the overlapped child allele is counted."""
        parent_alt = 'TGAAGG'                          # parent's own ALT for REF TGCAGG (C->A at 1012)
        child = (1012, 'C', ['G'], '1|1')
        parent = (1010, 'TGCAGG', [parent_alt], '1|0')
        t1 = apply([(1010, 'TGCAGG', parent_alt)])
        t2 = apply([(1012, 'C', 'G')])
        reg = Region(t1, t2)
        try:
            r = reg.score([child, parent])           # input order must not matter
            self.assertEqual(r['ed'], 0)
            self.assertEqual(r['records']['overlap_skipped_alleles'], [1, 0])
            self.assertEqual(r['records']['applied_alleles'], [1, 1])
        finally:
            reg.close()

    def test_crossing_record(self):
        # truth h1 starts 3 bp into the span (a deletion that begins before the span)
        t1 = REF[3:]
        reg = Region(t1, REF)
        try:
            r = reg.score([(997, 'NNNN' + REF[:3], ['N'], '1|0')])
            self.assertEqual(r['ed'], 0)
            self.assertEqual(r['records']['records_crossing_span_edge'], 1)
            self.assertIn('records_crossing_span_edge', r['flags'])
        finally:
            reg.close()

    def test_contig_check_and_sample(self):
        r = sc.score_region(self.reg.dir, self.reg.vcf(TRUTH_RECS, header_contigs=('chr2',)), truvari=False)
        self.assertTrue(r['status'].startswith('error: contig chr1 not in the VCF header'))
        r = self.reg.score(TRUTH_RECS, sample='SAMPLE1')
        self.assertEqual(r['ed'], 0)
        self.assertIn('sample_name_SAMPLE1_used_as_HG002', r['flags'])

    def test_shifted_coordinates_are_an_error(self):
        # every record one base to the right: REF no longer matches CHM13, which must not score as
        # a plausible ED (the whole VCF silently reduced to reference)
        shifted = [(p + 1, rf, a, g) for p, rf, a, g in TRUTH_RECS]
        r = self.reg.score(shifted)
        self.assertTrue(r['status'].startswith('error: REF disagrees with CHM13'), r['status'])
        self.assertEqual(r['records']['ref_mismatch_alleles_as_ref'], 2)   # the hom SNV counts once

    def test_records_on_other_contigs_are_ignored(self):
        p = os.path.join(self.reg.dir, 'two.vcf')
        body = vcf_text(TRUTH_RECS)
        other = vcf_text([(1003, 'G', ['A'], '1|1')], contig='chrX').splitlines()[-1]
        with open(p, 'w') as f:
            f.write(body + other + '\n')
        self.assertEqual(sc.score_region(self.reg.dir, p, truvari=False)['ed'], 0)


class TestOverlapModes(unittest.TestCase):
    """The two ways a truth-style VCF packs one haplotype into overlapping records, found in the
    T2T-Q100 truth itself; the trimmed rule must apply both, the raw rule drops one."""

    def test_conflict_rule(self):
        def e(pos, ref, alt):
            why, info, _ = sc.allele_edit({'pos': pos, 'ref': ref, 'alts': [alt]}, 1, A1, B1, REF)
            self.assertEqual(why, 'ok')
            return info

        def c(x, y, mode='trimmed'):
            return sc.conflict(x, y, mode), sc.conflict(y, x, mode)

        R = lambda p, n=1: REF[p - A1:p - A1 + n]   # noqa: E731
        ins = e(1010, R(1010), R(1010) + 'AA')                   # insertion before 1011
        self.assertEqual(c(ins, e(1010, R(1010), R(1010) + 'AA')), (True, True))       # duplicate
        self.assertEqual(c(ins, e(1010, R(1010), R(1010) + 'CC')), (True, True))       # same point
        # an insertion before 1011 and a SNV at 1011: REF spans share no base, both apply
        self.assertEqual(c(ins, e(1011, R(1011), 'A')), (False, False))
        # a deletion 1011..1015 whose REF (padding 1010) overlaps the insertion's: the insertion
        # sits on the deletion's left boundary and the deletion does not spell it
        self.assertEqual(c(ins, e(1010, R(1010, 6), R(1010))), (False, False))
        # the same insertion is spelled by a change that starts with it: redundant
        self.assertEqual(c(ins, e(1010, R(1010, 3), R(1010) + 'AA' + 'TT')), (True, True))
        # a padding-base SNV and a deletion from the next base: the changes do not overlap
        snv = e(1010, R(1010), 'A' if R(1010) != 'A' else 'C')
        self.assertEqual(c(snv, e(1010, R(1010, 6), R(1010))), (False, False))
        self.assertEqual(c(snv, e(1010, R(1010, 6), R(1010)), 'raw'), (True, True))
        # a long multi-allelic REF whose allele changes only its start, and a SNV inside it
        packed = e(1010, R(1010, 16), R(1010) + 'A' + R(1011, 15))
        self.assertEqual(c(packed, e(1020, R(1020), 'A' if R(1020) != 'A' else 'C')), (False, False))
        # overlapping changes
        self.assertEqual(c(e(1010, R(1010, 6), R(1010)), e(1012, R(1012), 'A' if R(1012) != 'A' else 'C')),
                         (True, True))

    def test_same_event_in_a_repeat(self):
        ref = 'GGT' + 'CA' * 6 + 'TTG'
        a1 = 1
        reg = {'pos': 3, 'ref': 'TCA', 'alts': ['T']}           # delete the first CA
        shifted = {'pos': 7, 'ref': 'ACA', 'alts': ['A']}       # the same deletion two units right
        _, x, _ = sc.allele_edit(reg, 1, a1, len(ref), ref)
        _, y, _ = sc.allele_edit(shifted, 1, a1, len(ref), ref)
        self.assertFalse(sc._same_event(x, y))                  # REF spans do not overlap: not tested
        shifted = {'pos': 5, 'ref': 'ACA', 'alts': ['A']}       # one unit right, REF overlaps
        _, y, _ = sc.allele_edit(shifted, 1, a1, len(ref), ref)
        self.assertTrue(sc._same_event(x, y))
        self.assertTrue(sc.conflict(x, y, 'trimmed'))

    def test_multiallelic_packing(self):
        # h1: 1 bp insertion after 1011 + SNV at 1020; h2: deletion 1011..1025. Written as the smvar
        # VCF writes it: one record with REF 1010..1025 (alleles: the deletion, and the insertion
        # spelled over the whole REF), plus the SNV with '*' for the deleted haplotype.
        ref_long = REF[9:25]
        ins_long = REF[9] + 'A' + REF[10:25]
        t1 = apply([(1010, REF[9], REF[9] + 'A'), (1020, REF[19], 'A' if REF[19] != 'A' else 'G')])
        snv = 'A' if REF[19] != 'A' else 'G'
        t2 = apply([(1010, ref_long, REF[9])])
        reg = Region(t1, t2)
        try:
            recs = [(1010, ref_long, [REF[9], ins_long], '2|1'), (1020, REF[19], [snv, '*'], '1|2')]
            r = reg.score(recs)
            self.assertEqual(r['ed'], 0)
            self.assertEqual(r['records']['star_alleles'], 1)
            self.assertNotIn('symbolic_alleles_as_ref', r['flags'])
            r = reg.score(recs, overlap='raw')
            self.assertEqual(r['ed'], 1)                       # the SNV is dropped under the raw rule
        finally:
            reg.close()

    def test_snv_on_a_deletions_padding_base(self):
        snv = 'A' if REF[9] != 'A' else 'C'
        t1 = apply([(1010, 'TGCAGG', snv)])
        reg = Region(t1, REF)
        try:
            recs = [(1010, 'T', [snv], '1|0'), (1010, 'TGCAGG', ['T'], '1|0')]
            self.assertEqual(reg.score(recs)['ed'], 0)
            self.assertEqual(reg.score(recs, overlap='raw')['ed'], 1)     # longest first: the SNV is dropped
        finally:
            reg.close()

    def test_duplicate_records_apply_once(self):
        reg = Region(T1, T2)
        try:
            recs = TRUTH_RECS + [(1010, 'TGCAGG', ['T', '*'], '1|2'), (1028, 'G', ['GAAAA'], '0|1')]
            for mode in ('trimmed', 'raw'):
                r = reg.score(recs, overlap=mode)
                self.assertEqual(r['ed'], 0, mode)
                self.assertEqual(r['records']['overlap_skipped_alleles'], [1, 1], mode)
        finally:
            reg.close()


class TestHaploid(unittest.TestCase):
    def test_chrx_outside_par(self):
        self.assertEqual(sc.truth_ploidy('chrX', 3252648, 3259266), 1)
        self.assertEqual(sc.truth_ploidy('chrX', 1000, 2000), 2)          # PAR1
        self.assertEqual(sc.truth_ploidy('chr7', 1000, 2000), 2)
        t1 = apply([SNV, DEL])
        reg = Region(t1, REF, contig='chrX', a1=5000001)
        try:
            recs = [(5000003, 'G', ['T'], '1'), (5000010, 'TGCAGG', ['T'], '1')]
            r = reg.score(recs, contig='chrX')
            self.assertEqual(r['truth_ploidy'], 1)
            self.assertEqual(r['ed'], 0)
            self.assertEqual(r['ed_ref'], 2 * ed(REF, t1))
            r = reg.score(recs[:1], contig='chrX')
            self.assertEqual(r['ed'], 2 * 5)                                # twice the haploid distance
            self.assertAlmostEqual(r['ed_per_kb'], 1000.0 * 5 / len(t1), places=3)
        finally:
            reg.close()


class TestPhaseSearchIsExact(unittest.TestCase):
    def test_random_against_brute_force(self):
        rnd = random.Random(7)
        ref = ''.join(rnd.choice('ACGT') for _ in range(300))
        a1 = 1
        for trial in range(25):
            t = [list(ref), list(ref)]
            recs = []
            positions = sorted(rnd.sample(range(5, 290, 3), 8))
            for p in positions:
                alt = rnd.choice([b for b in 'ACGT' if b != ref[p]])
                g = rnd.choice([(1, 0), (0, 1), (1, 1)])
                recs.append({'idx': len(recs), 'chrom': 'c', 'pos': p + 1, 'end': p + 1, 'id': '.', 'ref': ref[p],
                             'alts': [alt], 'filter': 'PASS', 'gt': list(g), 'phased': rnd.random() < 0.5,
                             'ps': rnd.choice(['1', '2', None])})
                for h in (0, 1):
                    if g[h] and rnd.random() < 0.85:       # truth: mostly the call, with some noise
                        t[h][p] = alt
            truth = (''.join(t[0]), ''.join(t[1]))
            prep = sc.Prepared(recs, a1, len(ref), ref, 2)
            cache = sc.EdCache(truth)
            got = sc.best_phase(prep, cache, enum_cap=12)
            brute = min(cache.pair(*prep.build(list(f))[0])[0]
                        for f in itertools.product([False, True], repeat=prep.n_blocks)) if prep.n_blocks \
                else cache.pair(*prep.build([])[0])[0]
            self.assertEqual(got['ed'], brute, trial)
            greedy = sc.best_phase(prep, cache, enum_cap=0)
            self.assertGreaterEqual(greedy['ed'], brute)
            self.assertLessEqual(greedy['ed'], got['ed_given_phase'])


class TestPerturb(unittest.TestCase):
    def test_kinds(self):
        region = {'a1': A1, 'b1': B1}
        big = (1020, REF[19:19 + 60 if 19 + 60 < len(REF) else len(REF) - 1], None)
        lines = vcf_text(TRUTH_RECS + [(1040, REF[39:52], [REF[39]], '1|1')]).splitlines()
        header = [x for x in lines if x.startswith('#')]
        body = [x for x in lines if not x.startswith('#')]
        self.assertIsNone(sc.perturb(region, body, header, 'drop_sv'))   # no SV-sized record here
        new, desc, exp = sc.perturb(region, body, header, 'swap_phase')
        self.assertEqual(exp, 'gt0')
        self.assertEqual(len(new), len(body))
        self.assertNotEqual(new, body)
        new, desc, exp = sc.perturb(region, body, header, 'unphase_swapped')
        self.assertEqual(exp, 'eq0')
        self.assertTrue(any('/' in x.split('\t')[9] for x in new))
        new, desc, exp = sc.perturb(region, body, header, 'ps_split_swapped')
        self.assertTrue(all(':' in x.split('\t')[9] for x in new))
        new, desc, exp = sc.perturb(region, body, header, 'flip_gt')
        self.assertEqual(exp, 'gt0')
        del big

    def test_perturbations_score_as_expected(self):
        reg = Region(T1, T2)
        try:
            lines = vcf_text(TRUTH_RECS).splitlines()
            header = [x for x in lines if x.startswith('#')]
            body = [x for x in lines if not x.startswith('#')]
            region = sc.load_region(reg.dir)
            for kind in ('swap_phase', 'flip_gt', 'unphase_swapped', 'ps_split_swapped'):
                new, desc, exp = sc.perturb(region, body, header, kind)
                p = os.path.join(reg.dir, kind + '.vcf')
                with open(p, 'w') as f:
                    f.write('\n'.join(header + new) + '\n')
                r = sc.score_region(reg.dir, p, truvari=False)
                if exp == 'gt0':
                    self.assertGreater(r['ed'], 0, kind)
                else:
                    self.assertEqual(r['ed'], 0, kind)
        finally:
            reg.close()


def at(*paths):
    """INFO/AT from node-number lists: at([1, 2, 20], [1, 20]) -> 'AT=>1>2>20,>1>20'."""
    return 'AT=' + ','.join(''.join('>%d' % n for n in p) for p in paths)


class TestNesting(unittest.TestCase):
    """vg nests a child snarl's record inside its parent's, and trims the shared prefix of the
    parent's alleles, so in a tandem repeat a parent's POS can lie to the right of its child's."""

    # parent snarl >1>20 deletes 1011..1030 (REF written from 1010); its child >5>8, written at
    # 1009 (left of the parent's POS), is a SNV at 1011, inside the parent's deletion
    PARENT = (1010, REF[9:30], [REF[9]], '1|1', None, None, '>1>20', at(range(1, 21), [1, 20]))
    CHILD = (1009, REF[8:11], [REF[8:10] + 'C'], '1|1', None, None, '>5>8', at([5, 6, 8], [5, 7, 8]))

    def test_ancestors_from_at(self):
        recs = [sc.parse_record(vcf_text([r]).splitlines()[-1], 0, i) for i, r in enumerate([self.PARENT, self.CHILD])]
        self.assertEqual(sc.vg_ancestors(recs), [set(), {0}])
        ordered, anc, st = sc.nesting_order(recs)
        self.assertEqual([r['idx'] for r in ordered], [0, 1])          # parent first despite POS
        self.assertEqual(st['descendant_before_ancestor_by_pos'], 1)
        # without vg IDs nothing is nested and the order is the plain (POS, -end, idx) order
        plain = [dict(r, bnd=None, interior=None) for r in recs]
        self.assertEqual([r['idx'] for r in sc.nesting_order(plain)[0]], [1, 0])
        self.assertEqual(sc.vg_ancestors(plain), [set(), set()])

    def test_parent_wins_whatever_its_pos(self):
        """The check that fired on scorer version 1: the child (smaller POS) was applied first and
        the parent's deletion skipped as a conflict."""
        t = apply([(1010, REF[9:30], REF[9])])
        reg = Region(t, t)
        try:
            r = reg.score([self.CHILD, self.PARENT])
            self.assertEqual(r['ed'], 0)
            self.assertEqual(r['called_len'], [len(t), len(t)])
            self.assertEqual(r['records']['overlap_skipped_alleles'], [1, 1])
            self.assertIn('descendant_before_ancestor_by_pos', r['flags'])
            # the same records without vg IDs (no nesting known): the plain order, child first
            r0 = reg.score([c[:6] for c in (self.CHILD, self.PARENT)])
            self.assertGreater(r0['ed'], 0)
        finally:
            reg.close()

    def test_child_applies_where_parent_is_reference(self):
        """Parent 1|0: on slot 2 the parent is reference and its child applies."""
        parent = self.PARENT[:3] + ('1|0',) + self.PARENT[4:]
        child = self.CHILD[:3] + ('.|1',) + self.CHILD[4:]
        t1 = apply([(1010, REF[9:30], REF[9])])
        t2 = apply([(1009, REF[8:11], REF[8:10] + 'C')])
        reg = Region(t1, t2)
        try:
            r = reg.score([child, parent])
            self.assertEqual(r['ed'], 0)
            self.assertEqual(r['records']['applied_under_nonref_parent'], [0, 0])
        finally:
            reg.close()

    def test_disagreeing_child_and_the_suppress_sensitivity(self):
        """A child insertion at the edge of the parent's change that the parent's ALT does not spell:
        the rule applies it (counted in applied_under_nonref_parent); the suppress_nested
        sensitivity leaves the parent's allele alone."""
        parent = (1010, REF[9:20], [REF[9]], '1|1', None, None, '>1>20', at(range(1, 21), [1, 20]))
        # insertion of AAA before base 1021 (right after the deleted stretch), in child >15>18
        child = (1020, REF[19], [REF[19] + 'AAA'], '1|1', None, None, '>15>18', at([15, 16, 18], [15, 17, 18]))
        t = apply([(1010, REF[9:20], REF[9])])
        reg = Region(t, t)
        try:
            r = reg.score([parent, child])
            self.assertEqual(r['records']['applied_under_nonref_parent'], [1, 1])
            self.assertEqual(r['ed'], 6)
            self.assertEqual(r['sensitivity']['ed_suppress_nested'], 0)
            self.assertIn('applied_under_nonref_parent', r['flags'])
        finally:
            reg.close()

    def test_cycle_is_released(self):
        a = (1010, REF[9:12], [REF[9]], '1|1', None, None, '>1>4', at([1, 2, 3, 5, 6, 4]))
        b = (1020, REF[19:22], [REF[19]], '1|1', None, None, '>5>6', at([5, 1, 4, 6]))
        recs = [sc.parse_record(vcf_text([r]).splitlines()[-1], 0, i) for i, r in enumerate([a, b])]
        ordered, _, st = sc.nesting_order(recs)
        self.assertEqual(len(ordered), 2)
        self.assertEqual(st['nesting_cycle_released'], 1)


class TestSummarise(unittest.TestCase):
    def test_paired_counts(self):
        d = tempfile.mkdtemp(prefix='sh_sum_')
        try:
            rows = {r['region_id']: r for r in sc.regions_table()}
            hs = [r for r in rows if rows[r]['stratum'] == 'hotspot_vntr'][:3]
            ctl = [r for r in rows if rows[r]['stratum'] == 'control_vntr_matched'][:2]
            vals = {'mc': {hs[0]: 100, hs[1]: 50, hs[2]: 0, ctl[0]: 5, ctl[1]: 0},
                    'new': {hs[0]: 80, hs[1]: 55, hs[2]: 0, ctl[0]: 5, ctl[1]: 20}}
            for lab, v in vals.items():
                os.makedirs(os.path.join(d, lab))
                for rid, e in v.items():
                    with open(os.path.join(d, lab, rid + '.json'), 'w') as f:
                        json.dump({'region_id': rid, 'stratum': rows[rid]['stratum'], 'status': 'ok', 'ed': e,
                                   'ed_ref': 200, 'ed_per_kb': e / 2.0, 'gain': 1 - e / 200.0, 'exact': e == 0,
                                   'truth_len': [1000, 1000], 'phase': {'n_blocks': 1, 'enum_cap_hit': False},
                                   'records': {}, 'flags': [], 'truvari': {}}, f)
            table, per_region = sc.summarise(['mc', 'new'], 'mc', results_dir=d, regions=','.join(hs + ctl),
                                             out=os.path.join(d, 'summary'))
            t = [x for x in table if x['label'] == 'new' and x['stratum'] == 'hotspot_vntr'][0]
            self.assertEqual((t['n'], t['better_gt0'], t['worse_gt0'], t['equal']), (3, 1, 1, 1))
            self.assertEqual((t['better_gt10'], t['worse_gt10']), (1, 0))
            self.assertEqual((t['ed_sum'], t['base_ed_sum'], t['exact']), (135, 150, 1))
            c = [x for x in table if x['label'] == 'new' and x['stratum'] == 'control_vntr_matched'][0]
            self.assertEqual((c['worse_gt0'], c['worse_gt10']), (1, 1))
            with open(os.path.join(d, 'summary.md')) as f:
                md = f.read()
            self.assertIn('controls worse: control_vntr_matched', md)
            self.assertTrue(os.path.exists(os.path.join(d, 'summary.tsv')))
        finally:
            shutil.rmtree(d)


if __name__ == '__main__':
    unittest.main(verbosity=1)
