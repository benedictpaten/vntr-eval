#!/usr/bin/env python3
"""score_haplotypes.py -- Stage 3 scorer: how close are the haplotypes a VCF writes over a region's
span to HG002's truth haplotypes?

    python3 tools/score_haplotypes.py score REGION_DIR VCF [--label L] [--out JSON] [--sample HG002]
                                      [--no-truvari] [--no-refine] [--enum-cap 12] [--labels-dir DIR]
    python3 tools/score_haplotypes.py batch --label L --vcf PATTERN [--regions pilot|all|ID,ID|FILE]
                                      [--jobs 4] [--out-dir results/stage3/L] [--force] [...]
    python3 tools/score_haplotypes.py summarise --labels mc,unit_aware__all,... [--baseline mc]
                                      [--regions pilot] [--out results/stage3_summary]
    python3 tools/score_haplotypes.py validate [--jobs 4]      # the scorer's own checks (see below)

PATTERN is a path with {id}, {contig} and {stratum} placeholders, for example
    work/stage3/call/mc/{id}.vcf.gz
    ~/PycharmProjects/vg-call-eval/work/wgs-mm095/{contig}/{contig}.vcf.gz   (production, genome-wide)

THE SCORING RULE (fixed before any Stage 3 call was scored)

Per region, over the anchor-to-anchor span [span_start, span_end] of region.json:

1. Called haplotypes. Every VCF record whose REF span overlaps the span is applied to the CHM13 span
   sequence (hap32.fa's CHM13 record), once per GT slot.
   - Phase: records with a phased GT share a phase block per PS value (phased records without PS
     share one block, as the VCF specification says); every unphased heterozygous record is a block
     of its own. The relative orientation of blocks is unknown, so it is chosen to minimise ED
     (below): exact enumeration of all 2^(k-1) assignments of k blocks while k-1 <= 12
     (--enum-cap); beyond that a greedy search (start from the VCF's own orientation, flip one block
     at a time while ED falls, repeat until a pass changes nothing), flagged 'phase_greedy'. The ED
     under the VCF's own orientation is reported as phase.ed_given_phase.
   - Records crossing a span boundary are included, clipped to the span, and flagged: the allele is
     trimmed to its minimal form (common suffix, then common prefix), and the bases of REF outside
     the span are cut together with the same number of leading/trailing ALT bases. A crossing
     record that trims to an insertion exactly on a span edge is outside the span and is dropped
     (flagged). A record wholly inside the span by its REF is applied as written.
   - Missing alleles ('.'), records with a FILTER other than PASS or '.', symbolic/breakend alleles
     and records whose REF disagrees with CHM13 are treated as reference, and counted. A '*' allele
     (the base is deleted by an upstream record) applies nothing and is counted, not flagged.
   - Overlapping records within one GT slot: records are taken in (POS, longest REF first, input)
     order, except that a vg record is never taken before a record it is nested in (nesting_order:
     C is nested in P when both of C's snarl boundary nodes, from vg's '>start>end' ID, are
     interior nodes of P's INFO/AT traversals), and an allele that conflicts with an allele already
     applied to that slot is skipped and counted (overlap_skipped_alleles; a reference allele
     occupies nothing). The nesting clause is how "a nested vg parent record wins" (below) is
     implemented: vg trims the shared prefix of a parent's alleles, so in a tandem repeat the
     parent's POS can lie to the right of its own child's, and scorer version 1 (plain POS order)
     then applied the child first and skipped the parent's allele as a conflict (MC calls at
     L011138, L012272, L001909). Records without vg IDs (the truth VCFs) have no nesting, so for
     them the order is the plain one. Two alleles whose REF
     spans share no base never conflict, as in any VCF without overlapping records. Alleles whose
     REF spans overlap conflict when the bases they change (trimmed intervals) overlap, when an
     insertion lies strictly inside the other's change or two insertions share a point, when they
     are the same event written twice, or when one is an insertion at the boundary of the other's
     change that the other's ALT already spells there. So a nested vg parent record wins where its
     change covers a child, and a child applies where the parent's allele is reference there.
     Decided on the truth VCF alone (validate): this rule gives ED 0 at all 149 regions for both the
     stvar truth and the redundant stvar+smvar union; the raw rule (any shared REF base conflicts:
     bcftools consensus, region.py) drops a SNV on a deletion's padding base (L001909) and 18
     regions of the union's multi-allelic packing; a purely trimmed rule double-applies 24 vg child
     insertions that repeat their parent's boundary insertion, and drops an insertion followed by a
     SNV at the next base (8 regions of the truth). The ED under the raw rule is kept in every JSON
     as sensitivity.ed_raw_overlap; the summary reports both.
     Where a child allele applies on a slot where a record it is nested in is also non-reference
     (the parent's trimmed change does not cover it: parent and child disagree, or trimming in a
     repeat hides that they write one event), the count is records.applied_under_nonref_parent,
     and sensitivity.ed_suppress_nested is the ED when every such allele is left out (the parent's
     allele alone spells the haplotype through its children). Neither changes the rule's ED.
   - Ploidy: HG002 is male, so outside the PARs of chrX/chrY the truth is one haplotype
     (truth.fa's HG002#2 is then the reference placeholder). There both sides are compared as a
     homozygous pair: truth (t1, t1), and a haploid GT fills both slots, so every ED is twice the
     haploid distance and gain and ED/kb are unchanged. A haploid GT elsewhere fills both slots too
     (flagged).
2. PRIMARY: ED = min(ed(c1,t1)+ed(c2,t2), ed(c1,t2)+ed(c2,t1)), exact unit-cost edit distance
   (tools/fastedit.c via region.edit_distance, Myers bit-vector). Reported per region, per kb of
   truth (ED over the summed truth length), 'exact' (ED == 0), ED_ref = the same for CHM13 on both
   haplotypes, and gain = 1 - ED/ED_ref (None when ED_ref == 0: 6 regions whose truth is CHM13).
3. SECONDARY: truvari bench with the genome-wide pipeline's parameters (bcftools norm -m-any, sort;
   --sizemin 50 --sizefilt 50 --pick ac, sample HG002) of the VCF's records against the normalised
   truth stvar (work/wgs-mm095/score/<contig>.truth.norm.vcf.gz), restricted to the span
   intersected with the SV benchmark; then truvari refine -u on its own candidate regions
   ('refined') and truvari refine -u -w over the whole span ('phab'). F1 uses TP-base for recall
   and TP-comp for precision. Raw F1 is a representation score only. The VCF goes to truvari as it
   is (FILTERed records included, as in the genome-wide scoring). With --labels-dir, the genome-wide
   run's own truvari labels are also counted over the span and the core.

Output: one JSON per region (results/stage3/<label>/<id>.json in batch mode).

The validate command runs the checks that must be able to fail, over all 149 regions, into
results/stage3_scorer_validation.tsv: the stvar truth VCF (what truth.fa was built from) and the
stvar+smvar union must give ED 0; perturbed stvar VCFs (the largest SV record dropped, one
heterozygous record's phase swapped, one GT flipped consistently) must give ED > 0; the swap with the
phase made unknown (unphased, or split into two PS blocks with one block swapped, also through the
greedy branch) must give ED 0 again; an empty VCF must give ED == ED_ref.
"""
import argparse
import collections
import concurrent.futures
import csv
import glob
import gzip
import hashlib
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import time

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402

SCORER_VERSION = 2   # 2: nesting_order (a vg parent is applied before its children whatever their POS)
SV = 50
TRUVARI = os.environ.get('VNTR_TRUVARI', os.path.join(config.EVAL_DIR, 'work/truvari-venv/bin/truvari'))
WGS_RUN = os.path.join(config.EVAL_DIR, 'work/wgs-mm095')
GENOMEWIDE_VCF = os.path.join(WGS_RUN, '{contig}', '{contig}.vcf.gz')
GENOMEWIDE_LABELS = os.path.join(WGS_RUN, 'score', '{contig}.truvari')
STAGE3_WORK = os.path.join(config.WORK_DIR, 'stage3')
STAGE3_RESULTS = os.path.join(config.RESULTS_DIR, 'stage3')
ENUM_CAP = 12
REF_MISMATCH_ERROR = 0.2   # a region errors when this share of its alleles has a REF that is not CHM13's
OVERLAP_MODE = 'trimmed'   # see conflict(); decided on the truth VCF alone (module docstring)
GREEDY_MAX_PASSES = 20
PASS_FILTERS = ('PASS', '.', '')

PILOT = ['L009656', 'L014297', 'L015415', 'L005990', 'L012184', 'L002013', 'L011138', 'L012272', 'L001909',
         'L016870',
         'L016124', 'L007687', 'L007172', 'L011430', 'L007065', 'L004021', 'L014111',
         'L000506', 'L006556']
STRATA_ORDER = ['hotspot_vntr', 'hotspot_other', 'control_vntr_matched', 'control_vntr_correct',
                'control_nontr_sv']
CONTROL_STRATA = ('control_vntr_matched', 'control_vntr_correct', 'control_nontr_sv')
HOTSPOT_STRATA = ('hotspot_vntr', 'hotspot_other')
# CHM13 chrX PARs (0-based half-open), as vg-call-eval/scripts/wgs/chrX.par.bed: read off the
# T2T-Q100 truth, which is diploid in exactly these two blocks and haploid elsewhere on chrX.
PARS = {'chrX': [(0, 2394370), (153926003, 154259566)]}
HAPLOID_CONTIGS = ('chrX', 'chrY')


class ScoreError(RuntimeError):
    pass


# ----------------------------------------------------------------------------- edit distance

def _load_edit_distance():
    try:
        from region import edit_distance as ed   # C (tools/fastedit.c) when a compiler exists
        return ed
    except Exception:  # noqa: BLE001  (region.py missing or broken: exact Python fallback)
        def ed(a, b):
            if not a:
                return len(b)
            if not b:
                return len(a)
            if len(a) > len(b):
                a, b = b, a
            m = len(a)
            mask = (1 << m) - 1
            high = 1 << (m - 1)
            peq = {}
            for i, c in enumerate(a):
                peq[c] = peq.get(c, 0) | (1 << i)
            vp, vn, score = mask, 0, m
            for c in b:
                eq = peq.get(c, 0)
                xv = eq | vn
                xh = (((eq & vp) + vp) ^ vp) | eq
                ph = vn | (~(xh | vp) & mask)
                mh = vp & xh
                if ph & high:
                    score += 1
                elif mh & high:
                    score -= 1
                ph = ((ph << 1) | 1) & mask
                mh = (mh << 1) & mask
                vp = mh | (~(xv | ph) & mask)
                vn = ph & xv
            return score
        return ed


edit_distance = _load_edit_distance()


class EdCache:
    """ed(seq, truth[i]) memoised on the sequence's digest (phase enumeration rebuilds the same
    haplotypes many times)."""

    def __init__(self, truths):
        self.truths = truths
        self.memo = {}
        self.computed = 0

    def ed(self, seq, i):
        k = (hashlib.sha1(seq.encode()).digest(), i)
        v = self.memo.get(k)
        if v is None:
            v = edit_distance(seq, self.truths[i])
            self.memo[k] = v
            self.computed += 1
        return v

    def pair(self, c1, c2):
        """(ED, pairing, d_h1, d_h2): the better pairing of called slots to truth haplotypes."""
        a = (self.ed(c1, 0), self.ed(c2, 1))
        b = (self.ed(c1, 1), self.ed(c2, 0))
        if sum(a) <= sum(b):
            return sum(a), 'slot1-h1,slot2-h2', a[0], a[1]
        return sum(b), 'slot1-h2,slot2-h1', b[1], b[0]


# ----------------------------------------------------------------------------- inputs

def read_fasta(path):
    recs, name, buf = [], None, []
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    recs.append((name, ''.join(buf).upper()))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(buf).upper()))
    return recs


def load_region(region_dir):
    region_dir = os.path.abspath(region_dir)
    with open(os.path.join(region_dir, 'region.json')) as f:
        rj = json.load(f)
    hap = read_fasta(os.path.join(region_dir, 'hap32.fa'))
    chm = [s for n, s in hap if n.startswith('CHM13')]
    if not chm:
        raise ScoreError('no CHM13 record in hap32.fa')
    refseq = chm[0]
    a1, b1 = int(rj['span_start']), int(rj['span_end'])
    if len(refseq) != b1 - a1 + 1:
        raise ScoreError('CHM13 record is %d bp, span is %d bp' % (len(refseq), b1 - a1 + 1))
    truth = dict(read_fasta(os.path.join(region_dir, 'truth.fa')))
    if 'HG002#1' not in truth or 'HG002#2' not in truth:
        raise ScoreError('truth.fa lacks HG002#1/HG002#2')
    return {'dir': region_dir, 'id': rj['region_id'], 'json': rj, 'contig': rj['contig'], 'a1': a1, 'b1': b1,
            'ref': refseq, 'truth': (truth['HG002#1'], truth['HG002#2']), 'stratum': rj.get('stratum'),
            'core': (rj.get('core_start'), rj.get('core_end'))}


def truth_ploidy(contig, a1, b1):
    """1 where HG002 (male) carries one copy: chrX/chrY outside the PARs; else 2."""
    if contig not in HAPLOID_CONTIGS:
        return 2
    for s, e in PARS.get(contig, []):
        if a1 - 1 < e and s < b1:
            return 2
    return 1


def _parse_info_end(info):
    for kv in info.split(';'):
        if kv.startswith('END='):
            try:
                return int(kv[4:])
            except ValueError:
                return None
    return None


VG_ID_RE = re.compile(r'^([<>])(\d+)([<>])(\d+)$')
NODE_RE = re.compile(r'[<>](\d+)')


def parse_vg_snarl(rec_id, info):
    """vg's snarl identity of a record: (boundary nodes, interior nodes of its INFO/AT traversals).
    vg call names a record by its snarl's two boundary nodes ('>start>end') and writes each allele's
    traversal in INFO/AT (the boundary nodes included; a difference block (SB) lists only its own
    interior stretch). Returns (None, None) for a record that is not a vg snarl record (any VCF
    without that ID form, e.g. the truth VCFs)."""
    m = VG_ID_RE.match(rec_id or '')
    if not m:
        return None, None
    bnd = (m.group(2), m.group(4))
    at = None
    for kv in (info or '').split(';'):
        if kv.startswith('AT='):
            at = kv[3:]
            break
    if at is None:
        return bnd, None
    nodes = set(NODE_RE.findall(at))
    nodes.discard(bnd[0])
    nodes.discard(bnd[1])
    return bnd, nodes


def vg_ancestors(recs):
    """For each record, the indices (into recs) of the records it is nested in: a vg record C is
    nested in a vg record P of a different snarl when both of C's boundary nodes are interior nodes
    of P's traversals (INFO/AT). A graph identity, not a coordinate or an ID order: vg trims the
    shared prefix of a parent's alleles, so in a tandem repeat a parent's POS can lie to the right
    of its own child's, and node IDs need not be topologically ordered. Records that are not vg
    snarl records have no ancestors."""
    holders = collections.defaultdict(set)
    for j, r in enumerate(recs):
        if r.get('interior'):
            for n in r['interior']:
                holders[n].add(j)
    out = []
    for i, r in enumerate(recs):
        b = r.get('bnd')
        if not b or b[0] not in holders or b[1] not in holders:
            out.append(set())
            continue
        out.append({j for j in holders[b[0]] & holders[b[1]] if j != i and recs[j].get('bnd') != b})
    return out


def nesting_order(recs):
    """The records in (POS, longest REF first, input order) order, except that no record precedes a
    record it is nested in (vg_ancestors): the smallest such key among the records whose ancestors
    have all been taken comes next. For records without nesting this is exactly the plain order.
    Returns (ordered records, ancestor sets keyed by input idx, stats)."""
    import heapq
    anc = vg_ancestors(recs)
    key = [(r['pos'], -r['end'], r['idx']) for r in recs]
    waiting = [len(a) for a in anc]
    kids = collections.defaultdict(list)
    for i, a in enumerate(anc):
        for j in a:
            kids[j].append(i)
    heap = [(key[i], i) for i in range(len(recs)) if not waiting[i]]
    heapq.heapify(heap)
    taken, order = [False] * len(recs), []
    stats = collections.Counter()
    while len(order) < len(recs):
        if not heap:        # a nesting cycle (not expected from a snarl tree): release the smallest key
            stats['nesting_cycle_released'] += 1
            i = min((key[i], i) for i in range(len(recs)) if not taken[i])[1]
        else:
            i = heapq.heappop(heap)[1]
            if taken[i]:
                continue
        taken[i] = True
        order.append(i)
        for c in kids[i]:
            waiting[c] -= 1
            if waiting[c] == 0 and not taken[c]:
                heapq.heappush(heap, (key[c], c))
    for i, a in enumerate(anc):
        if a:
            stats['nested_records'] += 1
        for j in a:
            if key[i] < key[j] and recs[i]['pos'] <= recs[j]['end'] and recs[j]['pos'] <= recs[i]['end']:
                stats['descendant_before_ancestor_by_pos'] += 1
    if [key[i] for i in order] != sorted(key):
        stats['order_changed_by_nesting'] = 1
    return [recs[i] for i in order], {recs[i]['idx']: {recs[j]['idx'] for j in a} for i, a in enumerate(anc)}, stats


def parse_record(line, sample_index, idx):
    f = line.rstrip('\n').split('\t')
    pos, ref = int(f[1]), f[3].upper()
    alts = [a.upper() for a in f[4].split(',')] if f[4] not in ('.', '') else []
    end = pos + len(ref) - 1
    if any(a.startswith('<') for a in alts):
        e = _parse_info_end(f[7])
        if e is not None:
            end = e
    bnd, interior = parse_vg_snarl(f[2], f[7] if len(f) > 7 else '')
    gt, ps, phased = [], None, True
    if len(f) > 9:
        keys = f[8].split(':')
        vals = f[9 + sample_index].split(':')
        d = dict(zip(keys, vals))
        g = d.get('GT', '.')
        phased = ('/' not in g)
        gt = [None if x in ('.', '') else int(x) for x in re.split(r'[|/]', g)]
        p = d.get('PS')
        ps = None if p in (None, '.', '') else p
    return {'idx': idx, 'chrom': f[0], 'pos': pos, 'end': end, 'id': f[2], 'ref': ref, 'alts': alts,
            'filter': f[6], 'gt': gt, 'phased': phased, 'ps': ps, 'bnd': bnd, 'interior': interior}


def _indexed(vcf):
    return vcf.endswith('.gz') and (os.path.exists(vcf + '.tbi') or os.path.exists(vcf + '.csi'))


def _have_bcftools():
    return bool(os.path.exists(config.BCFTOOLS) or shutil.which(config.BCFTOOLS, path=config.TOOL_PATH))


def vcf_header_lines(vcf):
    if _indexed(vcf) or vcf.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    out = []
    with opener(vcf, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                break
            out.append(line.rstrip('\n'))
    return out


def vcf_samples_contigs(header):
    samples, contigs = [], []
    for h in header:
        if h.startswith('#CHROM'):
            samples = h.split('\t')[9:]
        elif h.startswith('##contig=<'):
            m = re.search(r'ID=([^,>]+)', h)
            if m:
                contigs.append(m.group(1))
    return samples, contigs


def read_vcf_records(vcf, contig, a1, b1, sample='HG002'):
    """Records of `contig` whose REF span (INFO/END for symbolic alleles) overlaps a1..b1 (1-based).
    Returns (records, sample_used, notes). bcftools -r on an indexed VCF, else a Python scan."""
    header = vcf_header_lines(vcf)
    samples, contigs = vcf_samples_contigs(header)
    notes = []
    if contigs and contig not in contigs:
        raise ScoreError('contig %s not in the VCF header (%s...)' % (contig, ','.join(contigs[:5])))
    if not samples:
        si, used = None, None
    elif sample in samples:
        si, used = samples.index(sample), sample
    elif len(samples) == 1:
        si, used = 0, samples[0]
        notes.append('sample_name_%s_used_as_%s' % (samples[0], sample))
    else:
        raise ScoreError('sample %s not among %s' % (sample, samples))
    recs = []
    if _indexed(vcf) and _have_bcftools():
        p = subprocess.run([config.BCFTOOLS, 'view', '-H', '-r', '%s:%d-%d' % (contig, max(1, a1), b1), vcf],
                           capture_output=True, text=True, env=config.tool_env())
        if p.returncode != 0:
            raise ScoreError('bcftools view failed: %s' % p.stderr[-300:])
        lines = p.stdout.splitlines()
    else:
        opener = gzip.open if vcf.endswith('.gz') else open
        lines = []
        with opener(vcf, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                c, rest = line.split('\t', 1)
                if c != contig:
                    continue
                lines.append(line)
    for i, line in enumerate(lines):
        r = parse_record(line, si if si is not None else 0, i)
        if r['chrom'] != contig or r['pos'] > b1 or r['end'] < a1:
            continue
        if si is None:
            r['gt'] = []
        recs.append(r)
    return recs, used, notes


# ----------------------------------------------------------------------------- applying records

def trim_allele(pos, ref, alt):
    """Minimal representation: common suffix, then common prefix removed. Returns (s, e, alt) with
    ref[s..e) (1-based, half-open) replaced by alt; s == e is an insertion before base s."""
    r, a = ref, alt
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    k = 0
    while k < len(r) and k < len(a) and r[k] == a[k]:
        k += 1
    return pos + k, pos + len(r), a[k:]


def allele_edit(rec, allele, a1, b1, refseq):
    """How allele `allele` (>= 1) of rec changes the span. Returns
    ('ok', {'apply': (s, e, alt), 'raw': (pos, end), 'trim': (ts, te)}, crossing), where `apply`
    replaces 1-based half-open [s, e) of the span by alt (clipped to the span), `raw` is the
    record's REF span (inclusive) and `trim` the half-open interval the allele actually changes
    (ts == te: an insertion before base ts); or (reason, None, False) when the allele cannot be
    applied (it is then treated as reference)."""
    if allele > len(rec['alts']):
        return 'allele_index_out_of_range', None, False
    alt = rec['alts'][allele - 1]
    if alt == '*':
        return 'star', None, False               # spanned by an upstream deletion: nothing to apply
    if alt in ('.', '') or alt.startswith('<') or '[' in alt or ']' in alt or not re.fullmatch(r'[ACGTN]+', alt):
        return 'symbolic', None, False
    ref = rec['ref']
    pos, end = rec['pos'], rec['pos'] + len(ref) - 1
    # REF must agree with CHM13 where it lies inside the span
    lo, hi = max(pos, a1), min(end, b1)
    if lo <= hi and ref[lo - pos:hi - pos + 1] != refseq[lo - a1:hi - a1 + 1]:
        return 'ref_mismatch', None, False
    s, e, a = trim_allele(pos, ref, alt)
    info = {'raw': (pos, end), 'trim': (s, e), 'tseq': a, 'rec': (pos, ref, alt)}
    if pos >= a1 and end <= b1:
        info['apply'] = (s, e, a)                # same result as the raw record
        return 'ok', info, False
    if s == e:                                    # insertion before base s
        if a1 < s <= b1:
            info['apply'] = (s, e, a)
            return 'ok', info, False
        return 'edge_insertion', None, False
    if e - 1 < a1 or s > b1:                      # every changed base lies outside the span
        return 'outside_after_trim', None, False
    crossing = s < a1 or e - 1 > b1
    cs, ce = s, e
    if cs < a1:
        n = a1 - cs
        a = a[min(n, len(a)):]
        cs = a1
    if ce - 1 > b1:
        m = ce - 1 - b1
        a = a[:len(a) - min(m, len(a))]
        ce = b1 + 1
    info['apply'] = (cs, ce, a)
    return 'ok', info, crossing


def apply_edits(refseq, a1, edits):
    out, cur = [], a1
    for s, e, alt in sorted(edits, key=lambda x: (x[0], x[1])):
        if s < cur:
            raise ScoreError('overlapping edits %s' % ((s, e, alt),))
        out.append(refseq[cur - a1:s - a1])
        out.append(alt)
        cur = e
    out.append(refseq[cur - a1:])
    return ''.join(out)


def _same_event(x, y):
    """Whether two alleles whose REF spans overlap make the same change to the reference (one event
    written twice, e.g. an insertion or deletion placed at two positions of a repeat)."""
    (p1, r1, a1_), (p2, r2, a2_) = x['rec'], y['rec']
    lo, hi = min(p1, p2), max(p1 + len(r1), p2 + len(r2))
    win = [None] * (hi - lo)
    for p, r in ((p1, r1), (p2, r2)):
        for i, c in enumerate(r):
            win[p - lo + i] = c
    if None in win:
        return False
    w = ''.join(win)
    return w[:p1 - lo] + a1_ + w[p1 - lo + len(r1):] == w[:p2 - lo] + a2_ + w[p2 - lo + len(r2):]


def conflict(x, y, mode):
    """Whether two alleles cannot both be applied to one haplotype.

    raw: their REF spans share a base (bcftools consensus's rule, and region.py's).
    trimmed (the rule): alleles whose REF spans share no base never conflict, as in any VCF without
    overlapping records. Alleles whose REF spans overlap conflict when
      - the bases they change (trimmed intervals) overlap, an insertion lies strictly inside the
        other's change, or two insertions share a point;
      - they are the same event written twice (_same_event);
      - one is an insertion at the boundary of the other's change that the other's ALT already
        spells there (vg's nested child records repeat their parent's boundary insertion so).
    The raw rule drops real variants where a truth VCF packs one haplotype into overlapping records
    (smvar's multi-allelic records; a SNV on a deletion's padding base); the purely trimmed rule
    double-applies vg's nested child insertions and drops an insertion followed by a SNV at the next
    base (see the module docstring)."""
    (rs1, re1), (rs2, re2) = x['raw'], y['raw']
    raw_overlap = max(rs1, rs2) <= min(re1, re2)
    if mode == 'raw':
        return raw_overlap
    if not raw_overlap:
        return False
    (s1, e1), (s2, e2) = x['trim'], y['trim']
    i1, i2 = s1 == e1, s2 == e2
    if i1 and i2:
        if s1 == s2:
            return True
    elif i1 or i2:
        (p, ins), (s, e, alt) = ((s1, x['tseq']), (s2, e2, y['tseq'])) if i1 else ((s2, y['tseq']), (s1, e1, x['tseq']))
        if s < p < e:
            return True
        if (p == e and alt.endswith(ins)) or (p == s and alt.startswith(ins)):
            return True
    elif max(s1, s2) < min(e1, e2):
        return True
    return _same_event(x, y)


class Prepared:
    """The records of one region reduced to what the haplotype builder needs, with every flag.

    Records are taken in (POS, longest REF first, input order) order, except that a vg record never
    precedes a record it is nested in (nesting_order); within one GT slot an allele that conflicts
    (see conflict()) with an allele already applied to that slot is skipped. Before scorer version
    2 the order was plain (POS, longest REF, input): where vg's prefix trimming put a child's POS
    left of its parent's, the child was applied first and the parent's allele skipped as a conflict.

    Also kept for each item: its GT's non-reference slots and its ancestors among the items, so
    that build() can count (and, with suppress_nested, remove) alleles applied on a slot where a
    record they are nested in is non-reference too -- a parent and a child that the conflict rule
    lets both apply, i.e. where they disagree or where trimming hides that they repeat one event."""

    def __init__(self, recs, a1, b1, refseq, ploidy, apply_filtered=False, overlap=None):
        self.a1, self.b1, self.refseq = a1, b1, refseq
        self.overlap = overlap or OVERLAP_MODE
        self.flags = collections.Counter()
        self.items = []          # ((edit id slot 1, edit id slot 2), block or None)
        self.item_nonref = []    # (slot 1 non-ref, slot 2 non-ref) per item, from the GT
        self.item_anc = []       # item indices of the item's ancestors (all earlier in self.items)
        self.edits = []          # allele edits (allele_edit info dicts), shared by both slots of a hom
        blocks = collections.OrderedDict()
        ps_seen = set()
        ordered, anc_by_idx, nstats = nesting_order(recs)
        self.flags.update(nstats)
        item_of_idx = {}
        for r in ordered:
            self.flags['records_in_span'] += 1
            gt = list(r['gt'])
            if not gt:
                self.flags['no_gt'] += 1
                continue
            if len(gt) == 1:
                self.flags['haploid_gt'] += 1
                if ploidy == 2:
                    self.flags['haploid_gt_in_diploid_region'] += 1
                gt = gt * 2
            elif len(gt) > 2:
                self.flags['polyploid_gt_truncated'] += 1
                gt = gt[:2]
            elif ploidy == 1 and gt[0] != gt[1]:
                self.flags['diploid_het_gt_in_haploid_region'] += 1
            nonref = [a for a in gt if a]
            if nonref:
                self.flags['nonref_records'] += 1
            missing = sum(1 for a in gt if a is None)
            if missing:
                self.flags['missing_alleles_as_ref'] += missing
            gt = [0 if a is None else a for a in gt]
            if r['filter'] not in PASS_FILTERS and any(gt):
                if apply_filtered:
                    self.flags['filtered_records_applied'] += 1
                else:
                    self.flags['filtered_records_as_ref'] += 1
                    continue
            ids = []
            crossing_any = False
            by_allele = {}
            for a in gt:
                if a == 0:
                    ids.append(None)
                    continue
                if a in by_allele:
                    ids.append(by_allele[a])
                    continue
                why, info, crossing = allele_edit(r, a, a1, b1, refseq)
                if why == 'star':
                    self.flags['star_alleles'] += 1
                    by_allele[a] = None
                    ids.append(None)
                elif why != 'ok':
                    self.flags[why + '_alleles_as_ref'] += 1
                    by_allele[a] = None
                    ids.append(None)
                else:
                    info['item'] = len(self.items)
                    self.edits.append(info)
                    by_allele[a] = len(self.edits) - 1
                    ids.append(by_allele[a])
                    crossing_any = crossing_any or crossing
            if crossing_any:
                self.flags['records_crossing_span_edge'] += 1
            if ids[0] is None and ids[1] is None:
                continue
            block = None
            if ids[0] != ids[1]:                 # orientation matters
                if r['phased']:
                    block = ('PS', r['ps'])
                    ps_seen.add(r['ps'])
                else:
                    block = ('UNPHASED', r['idx'])
                    self.flags['unphased_het_records'] += 1
                blocks.setdefault(block, len(blocks))
            item_of_idx[r['idx']] = len(self.items)
            self.items.append(((ids[0], ids[1]), None if block is None else blocks[block]))
            self.item_nonref.append((bool(gt[0]), bool(gt[1])))
            self.item_anc.append(sorted(item_of_idx[j] for j in anc_by_idx.get(r['idx'], ()) if j in item_of_idx))
        self.n_blocks = len(blocks)
        self.block_keys = list(blocks)
        self.n_phase_sets = len(ps_seen)
        # the conflict graph between alleles of different records (orientation-independent)
        self.conf = [set() for _ in self.edits]
        # every conflict needs overlapping REF spans (inclusive), so sweep on those
        order = sorted(range(len(self.edits)), key=lambda i: self.edits[i]['raw'][0])
        for x, i in enumerate(order):
            ei = self.edits[i]
            hi = ei['raw'][1]
            for j in order[x + 1:]:
                ej = self.edits[j]
                if ej['raw'][0] > hi:
                    break
                if ei['item'] != ej['item'] and conflict(ei, ej, self.overlap):
                    self.conf[i].add(j)
                    self.conf[j].add(i)
        self.n_conflicting_pairs = sum(len(c) for c in self.conf) // 2

    def build(self, flips, suppress_nested=False):
        """Both slot sequences for a block orientation (flips[b] True swaps block b's slots).
        suppress_nested (a sensitivity, not the rule): an allele is not applied on a slot where any
        record it is nested in has a non-reference allele, i.e. the parent's allele alone spells
        the haplotype through its children."""
        applied = (set(), set())
        nonref_items = (set(), set())
        edits = ([], [])
        skipped = [0, 0]
        under_nonref_parent = [0, 0]
        for k, ((x0, x1), b) in enumerate(self.items):
            nr = self.item_nonref[k]
            if b is not None and flips[b]:
                x0, x1 = x1, x0
                nr = (nr[1], nr[0])
            anc = self.item_anc[k]
            for h, x in ((0, x0), (1, x1)):
                if nr[h]:
                    nonref_items[h].add(k)
                if x is None:
                    continue
                nested_under = bool(anc) and not nonref_items[h].isdisjoint(anc)
                if nested_under and suppress_nested:
                    skipped[h] += 1
                    continue
                if not self.conf[x].isdisjoint(applied[h]):
                    skipped[h] += 1
                    continue
                applied[h].add(x)
                edits[h].append(self.edits[x]['apply'])
                under_nonref_parent[h] += nested_under
        seqs = (apply_edits(self.refseq, self.a1, edits[0]), apply_edits(self.refseq, self.a1, edits[1]))
        self.last_under_nonref_parent = under_nonref_parent
        return seqs, skipped, (len(edits[0]), len(edits[1]))


def best_phase(prep, cache, enum_cap=ENUM_CAP, suppress_nested=False):
    """Choose the block orientation minimising ED. Returns a dict with the chosen sequences."""
    k = prep.n_blocks
    evaluated = 0

    def score(flips):
        nonlocal evaluated
        evaluated += 1
        seqs, skipped, applied = prep.build(flips, suppress_nested=suppress_nested)
        ed, pairing, d1, d2 = cache.pair(seqs[0], seqs[1])
        return ed, pairing, d1, d2, seqs, skipped, applied

    given = [False] * k
    res_given = score(given)
    best, best_flips = res_given, given
    if k <= 1:
        method = 'fully_phased' if k == 1 else 'no_heterozygous_records'
    elif k - 1 <= enum_cap:
        method = 'enumerated'
        for mask in range(1, 1 << (k - 1)):
            flips = [False] + [bool(mask >> (i - 1) & 1) for i in range(1, k)]
            r = score(flips)
            if r[0] < best[0]:
                best, best_flips = r, flips
    else:
        method = 'greedy'
        cur = list(given)
        for _ in range(GREEDY_MAX_PASSES):
            improved = False
            for b in range(1, k):
                trial = list(cur)
                trial[b] = not trial[b]
                r = score(trial)
                if r[0] < best[0]:
                    best, cur, improved = r, trial, True
            if not improved:
                break
        best_flips = cur
    ed, pairing, d1, d2, seqs, skipped, applied = best
    prep.build(best_flips, suppress_nested=suppress_nested)      # the chosen orientation's counts
    return {'ed': ed, 'pairing': pairing, 'd_h1': d1, 'd_h2': d2, 'seqs': seqs, 'skipped': skipped,
            'applied': applied, 'under_nonref_parent': list(prep.last_under_nonref_parent),
            'method': method, 'n_blocks': k, 'evaluated': evaluated,
            'ed_given_phase': res_given[0],
            'flipped_blocks': [i for i, f in enumerate(best_flips) if f],
            'enum_cap_hit': method == 'greedy'}


def score_haplotypes(region, recs, enum_cap=ENUM_CAP, apply_filtered=False, overlap=None):
    """The primary score of one region for already-read records."""
    a1, b1, ref = region['a1'], region['b1'], region['ref']
    ploidy = truth_ploidy(region['contig'], a1, b1)
    t1, t2 = region['truth']
    notes = []
    if ploidy == 1:
        if t2 != ref:
            notes.append('haploid_region_but_truth_h2_is_not_the_reference')
        t2 = t1
    cache = EdCache((t1, t2))
    prep = Prepared(recs, a1, b1, ref, ploidy, apply_filtered=apply_filtered, overlap=overlap)
    ph = best_phase(prep, cache, enum_cap)
    ed_ref = cache.ed(ref, 0) + cache.ed(ref, 1)
    sens = {}
    if prep.overlap != 'raw':
        # the same score under the raw overlap rule (bcftools consensus / region.py): how much the
        # rule for overlapping records moves this region
        pr = Prepared(recs, a1, b1, ref, ploidy, apply_filtered=apply_filtered, overlap='raw')
        sens['ed_raw_overlap'] = best_phase(pr, cache, enum_cap)['ed']
    if prep.flags.get('nested_records'):
        # a parent's non-reference allele suppresses every allele nested in it on that slot: how
        # much the rule's handling of parent/child disagreement (a child applies where the
        # parent's trimmed change does not cover it) moves this region
        sens['ed_suppress_nested'] = best_phase(prep, cache, enum_cap, suppress_nested=True)['ed']
    else:
        sens['ed_suppress_nested'] = ph['ed']
    tlen = len(t1) + len(t2)
    flags = dict(prep.flags)
    flags['overlap_skipped_alleles'] = list(ph['skipped'])
    flags['applied_alleles'] = list(ph['applied'])
    flags['applied_under_nonref_parent'] = list(ph['under_nonref_parent'])
    warn = []
    for k in ('records_crossing_span_edge', 'edge_insertion_alleles_as_ref', 'missing_alleles_as_ref',
              'filtered_records_as_ref', 'filtered_records_applied', 'symbolic_alleles_as_ref',
              'ref_mismatch_alleles_as_ref', 'allele_index_out_of_range_alleles_as_ref',
              'haploid_gt_in_diploid_region', 'diploid_het_gt_in_haploid_region', 'polyploid_gt_truncated',
              'unphased_het_records', 'nesting_cycle_released'):
        if flags.get(k):
            warn.append(k)
    if sum(ph['skipped']):
        warn.append('overlap_skipped_alleles')
    if flags.get('descendant_before_ancestor_by_pos'):
        warn.append('descendant_before_ancestor_by_pos')
    if sum(ph['under_nonref_parent']):
        warn.append('applied_under_nonref_parent')
    if ph['n_blocks'] > 1:
        warn.append('phase_optimised_%s' % ph['method'])
    if ph['enum_cap_hit']:
        warn.append('phase_greedy')
    warn.extend(notes)
    c1, c2 = ph['seqs']
    return {
        'truth_ploidy': ploidy,
        'ref_len': len(ref),
        'truth_len': [len(t1), len(t2)],
        'called_len': [len(c1), len(c2)],
        'ed': ph['ed'],
        'ed_per_kb': round(1000.0 * ph['ed'] / tlen, 4) if tlen else None,
        'exact': ph['ed'] == 0,
        'ed_ref': ed_ref,
        'ed_ref_per_kb': round(1000.0 * ed_ref / tlen, 4) if tlen else None,
        'gain': round(1.0 - ph['ed'] / ed_ref, 5) if ed_ref else None,
        'pairing': ph['pairing'], 'ed_h1': ph['d_h1'], 'ed_h2': ph['d_h2'],
        'phase': {'method': ph['method'], 'n_blocks': ph['n_blocks'], 'n_phase_sets': prep.n_phase_sets,
                  'assignments_evaluated': ph['evaluated'], 'enum_cap': enum_cap,
                  'enum_cap_hit': ph['enum_cap_hit'], 'ed_given_phase': ph['ed_given_phase'],
                  'flipped_blocks': ph['flipped_blocks'][:50]},
        'records': flags,
        'overlap_mode': prep.overlap,
        'conflicting_allele_pairs': prep.n_conflicting_pairs,
        'flags': warn,
        'sensitivity': sens,
        'edit_distances_computed': cache.computed,
        '_seqs': (c1, c2),
    }


# ----------------------------------------------------------------------------- truvari

def run(cmd, **kw):
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, env=config.tool_env(), **kw)


def sh(cmd, **kw):
    p = run(cmd, **kw)
    if p.returncode != 0:
        raise ScoreError('%s failed: %s' % (os.path.basename(str(cmd[0])) + ' ' + str(cmd[1]), p.stderr[-400:]))
    return p


def f1(tpb, fn, tpc, fp):
    if None in (tpb, fn, tpc, fp):
        return None
    rec = tpb / (tpb + fn) if tpb + fn else None
    prec = tpc / (tpc + fp) if tpc + fp else None
    if rec is None and prec is None:
        return None
    if not rec or not prec:
        return 0.0
    return round(2 * rec * prec / (rec + prec), 4)


def bench_bed_path(contig):
    p = os.path.join(WGS_RUN, contig, 'truth.%s.stvar.bed' % contig)
    return p if os.path.exists(p) else config.data_paths(contig)['stvar_bed']


def truvari_score(region, vcf, sample_used, od, refine=True, refine_timeout=900, window=2000):
    """Secondary: truvari bench / refine / whole-span phab of the VCF's records against the truth
    stvar over the span intersected with the SV benchmark (evaluate.py's recipe)."""
    contig, a1, b1 = region['contig'], region['a1'], region['b1']
    ref_fa = config.data_paths(contig)['ref_fa']
    truth_norm = os.path.join(WGS_RUN, 'score', '%s.truth.norm.vcf.gz' % contig)
    bench_bed = bench_bed_path(contig)
    for x in (ref_fa, truth_norm, bench_bed, TRUVARI):
        if not os.path.exists(x):
            return {'status': 'missing_input:%s' % x}
    if os.path.exists(od):
        shutil.rmtree(od)
    os.makedirs(od)
    iv = []
    with open(bench_bed) as f:
        for line in f:
            x = line.split('\t')
            if x[0] != contig:
                continue
            a, b = max(int(x[1]), a1 - 1), min(int(x[2]), b1)
            if a < b:
                iv.append((a, b))
    res = {'bench_bp': sum(b - a for a, b in iv), 'span_bp': b1 - a1 + 1}
    if not iv:
        res['status'] = 'span_outside_sv_benchmark'
        return res
    bed = os.path.join(od, 'region.bed')
    with open(bed, 'w') as f:
        for a, b in iv:
            f.write('%s\t%d\t%d\n' % (contig, a, b))
    win = '%s:%d-%d' % (contig, max(1, a1 - window), b1 + window)
    truth = os.path.join(od, 'truth.vcf.gz')
    sh([config.BCFTOOLS, 'view', '-r', win, '-Oz', '-o', truth, truth_norm])
    sh([config.BCFTOOLS, 'index', '-t', '-f', truth])
    # the calls: window, one sample named HG002, split multiallelics and left-align, sort
    src = vcf
    if not _indexed(vcf):
        src = os.path.join(od, 'calls.input.vcf.gz')
        sh([config.BCFTOOLS, 'sort', '-Oz', '-o', src, vcf])
        sh([config.BCFTOOLS, 'index', '-t', '-f', src])
    sel = os.path.join(od, 'calls.window.vcf.gz')
    cmd = [config.BCFTOOLS, 'view', '-r', win, '-Oz', '-o', sel]
    if sample_used:
        cmd[2:2] = ['-s', sample_used]
    sh(cmd + [src])
    if sample_used and sample_used != 'HG002':
        names = os.path.join(od, 'sample.txt')
        with open(names, 'w') as f:
            f.write('HG002\n')
        rh = os.path.join(od, 'calls.renamed.vcf.gz')
        sh([config.BCFTOOLS, 'reheader', '-s', names, '-o', rh, sel])
        sel = rh
    unsorted = os.path.join(od, 'calls.norm.unsorted.vcf.gz')
    calls = os.path.join(od, 'calls.vcf.gz')
    sh([config.BCFTOOLS, 'norm', '-m-any', '-f', ref_fa, '-Oz', '-o', unsorted, sel])
    sh([config.BCFTOOLS, 'sort', '-Oz', '-o', calls, unsorted])
    sh([config.BCFTOOLS, 'index', '-t', '-f', calls])
    os.remove(unsorted)
    # SV-sized called alleles inside the span (a count of the records, for context)
    p = sh([config.BCFTOOLS, 'query', '-r', '%s:%d-%d' % (contig, a1, b1), '-f', '%POS\t%REF\t%ALT\t[%GT]\n', calls])
    n_sv = 0
    for line in p.stdout.splitlines():
        pos, rf, alt, gt = line.split('\t')
        called = [x for x in re.split(r'[|/]', gt) if x not in ('.', '0')]
        if called and abs(len(alt) - len(rf)) >= SV and not alt.startswith('<'):
            n_sv += 1
    res['call_sv_records_in_span'] = n_sv
    bd = os.path.join(od, 'bench')
    r = run([TRUVARI, 'bench', '-b', truth, '-c', calls, '-f', ref_fa, '-o', bd,
             '--includebed', bed, '--sizemin', '50', '--sizefilt', '50',
             '--bSample', 'HG002', '--cSample', 'HG002', '--pick', 'ac'])
    sj = os.path.join(bd, 'summary.json')
    if r.returncode != 0 or not os.path.exists(sj):
        res['status'] = 'bench_failed: ' + r.stderr[-300:]
        return res
    s = json.load(open(sj))
    for k in ('TP-base', 'FN', 'TP-comp', 'FP'):
        res['raw_' + k] = s.get(k)
    res['raw_f1'] = f1(s.get('TP-base'), s.get('FN'), s.get('TP-comp'), s.get('FP'))
    res['status'] = 'ok'
    if not refine:
        return res
    pb = bd + '_phab'
    shutil.copytree(bd, pb)
    for tag, d, extra in (('refined', bd, []), ('phab', pb, ['-w', '--regions', os.path.abspath(bed)])):
        t = time.time()
        vs = os.path.join(d, 'refine.variant_summary.json')
        st = None
        try:
            subprocess.run([TRUVARI, 'refine', '-u'] + extra + ['-t', '1', '.'], cwd=d, capture_output=True,
                           text=True, timeout=refine_timeout, env=config.tool_env())
        except subprocess.TimeoutExpired:
            st = 'timeout'
        if os.path.exists(vs):
            s2 = json.load(open(vs))
            for k in ('TP-base', 'FN', 'TP-comp', 'FP'):
                res['%s_%s' % (tag, k)] = s2.get(k)
            res[tag + '_f1'] = f1(s2.get('TP-base'), s2.get('FN'), s2.get('TP-comp'), s2.get('FP'))
            st = 'ok'
        elif st is None:
            rl = os.path.join(d, 'refine.log.txt')
            if os.path.exists(rl) and 'No regions to be refined' in open(rl).read():
                for k in ('TP-base', 'FN', 'TP-comp', 'FP', 'f1'):     # nothing to harmonise: = raw
                    res['%s_%s' % (tag, k)] = res['raw_' + k]
                st = 'nothing_to_refine'
            else:
                st = 'failed'
        res[tag + '_status'] = st
        res[tag + '_seconds'] = round(time.time() - t, 1)
    return res


def label_counts(region, labels_dir):
    """Genome-wide truvari labels (tp-base/fn/tp-comp/fp) of records overlapping the span and the
    core, for comparison with regions.tsv vg_fp / vg_fn (counted in the core)."""
    contig, a1, b1 = region['contig'], region['a1'], region['b1']
    c1, c2 = region['core']
    out = {'dir': labels_dir}
    for name in ('tp-base', 'fn', 'tp-comp', 'fp'):
        path = os.path.join(labels_dir, name + '.vcf.gz')
        if not os.path.exists(path):
            out['status'] = 'missing %s' % path
            return out
        lo = min(a1, c1 or a1)
        hi = max(b1, c2 or b1)
        p = sh([config.BCFTOOLS, 'query', '-r', '%s:%d-%d' % (contig, lo, hi), '-f', '%POS\t%REF\n', path])
        span = core = core_pos = 0
        for line in p.stdout.splitlines():
            pos, rf = line.split('\t')
            s, e = int(pos), int(pos) + len(rf) - 1
            if s <= b1 and e >= a1:
                span += 1
            if c1 and s <= c2 and e >= c1:
                core += 1
            if c1 and c1 <= s <= c2:
                core_pos += 1
        out['span_' + name] = span
        out['core_' + name] = core
        out['corepos_' + name] = core_pos
    out['status'] = 'ok'
    return out


# ----------------------------------------------------------------------------- one region

def score_region(region_dir, vcf, label=None, sample='HG002', truvari=True, refine=True, enum_cap=ENUM_CAP,
                 labels_dir=None, work=None, apply_filtered=False, write_haplotypes=None, refine_timeout=900,
                 overlap=None):
    t0 = time.time()
    region = load_region(region_dir)
    out = {'region_id': region['id'], 'stratum': region['stratum'], 'contig': region['contig'],
           'span_start': region['a1'], 'span_end': region['b1'], 'core_start': region['core'][0],
           'core_end': region['core'][1], 'label': label, 'vcf': os.path.abspath(vcf),
           'scorer_version': SCORER_VERSION, 'apply_filtered': apply_filtered}
    try:
        recs, used, notes = read_vcf_records(vcf, region['contig'], region['a1'], region['b1'], sample)
    except ScoreError as e:
        out['status'] = 'error: %s' % e
        return out
    out['sample'] = used
    prim = score_haplotypes(region, recs, enum_cap=enum_cap, apply_filtered=apply_filtered, overlap=overlap)
    seqs = prim.pop('_seqs')
    bad = prim['records'].get('ref_mismatch_alleles_as_ref', 0)
    if bad:
        ok_alleles = sum(prim['records']['applied_alleles']) + sum(prim['records']['overlap_skipped_alleles'])
        if bad >= REF_MISMATCH_ERROR * (bad + ok_alleles):
            out.update(prim)
            out['status'] = ('error: REF disagrees with CHM13 for %d of %d alleles (wrong contig offset or '
                             'coordinates?)' % (bad, bad + ok_alleles))
            return out
    out.update(prim)
    out['flags'] = notes + out['flags']
    out['status'] = 'ok'
    if write_haplotypes:
        os.makedirs(os.path.dirname(os.path.abspath(write_haplotypes)), exist_ok=True)
        with open(write_haplotypes, 'w') as f:
            for i, s in enumerate(seqs, start=1):
                f.write('>%s#called%d %s:%d-%d len=%d\n%s\n' % (region['id'], i, region['contig'], region['a1'],
                                                              region['b1'], len(s), s))
    if truvari:
        work = work or os.path.join(STAGE3_WORK, 'score', label or 'unlabelled', region['id'])
        try:
            out['truvari'] = truvari_score(region, os.path.abspath(vcf), used, os.path.join(work, 'truvari'),
                                           refine=refine, refine_timeout=refine_timeout)
        except ScoreError as e:
            out['truvari'] = {'status': 'error: %s' % e}
    else:
        out['truvari'] = {'status': 'skipped'}
    if labels_dir:
        try:
            out['genomewide_labels'] = label_counts(region, labels_dir.format(contig=region['contig']))
        except ScoreError as e:
            out['genomewide_labels'] = {'status': 'error: %s' % e}
    out['seconds'] = round(time.time() - t0, 2)
    return out


# ----------------------------------------------------------------------------- batch

def regions_table(path=None):
    path = path or os.path.join(config.REGIONS_DIR, 'regions.tsv')
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def select_regions(spec):
    rows = regions_table()
    ids = [r['region_id'] for r in rows]
    if spec in (None, 'all'):
        return ids
    if spec == 'pilot':
        return list(PILOT)
    if os.path.exists(spec):
        with open(spec) as f:
            want = [x.split()[0] for x in f if x.strip() and not x.startswith('#')]
    else:
        want = [x for x in spec.split(',') if x]
    bad = [x for x in want if x not in ids]
    if bad:
        raise SystemExit('unknown regions: %s' % ','.join(bad))
    return want


def _batch_one(kw):
    rid, out_path = kw.pop('rid'), kw.pop('out_path')
    try:
        res = score_region(**kw)
    except Exception as e:  # noqa: BLE001  (one region must not end a batch)
        res = {'region_id': rid, 'label': kw.get('label'), 'vcf': kw.get('vcf'), 'status': 'error: %r' % (e,)}
    tmp = out_path + '.tmp%d' % os.getpid()
    with open(tmp, 'w') as f:
        json.dump(res, f, indent=1)
    os.replace(tmp, out_path)
    return rid, res.get('status'), res.get('ed'), res.get('ed_ref')


def batch(label, pattern, regions, out_dir=None, jobs=4, force=False, truvari=True, refine=True,
          enum_cap=ENUM_CAP, labels_dir=None, sample='HG002', apply_filtered=False, quiet=False,
          missing_ok=False, overlap=None):
    out_dir = out_dir or os.path.join(STAGE3_RESULTS, label)
    os.makedirs(out_dir, exist_ok=True)
    rows = {r['region_id']: r for r in regions_table()}
    tasks = []
    missing = []
    for rid in regions:
        r = rows[rid]
        vcf = os.path.expanduser(pattern.format(id=rid, contig=r['contig'], stratum=r['stratum']))
        out_path = os.path.join(out_dir, rid + '.json')
        if not force and os.path.exists(out_path):
            try:
                if json.load(open(out_path)).get('status') == 'ok':
                    continue
            except ValueError:
                pass
        if not os.path.exists(vcf):
            missing.append((rid, vcf))
            if not missing_ok:
                with open(out_path, 'w') as f:
                    json.dump({'region_id': rid, 'stratum': r['stratum'], 'label': label, 'vcf': vcf,
                               'status': 'error: no VCF'}, f, indent=1)
            continue
        tasks.append({'rid': rid, 'out_path': out_path, 'region_dir': os.path.join(config.REGIONS_DIR, rid),
                      'vcf': vcf, 'label': label, 'sample': sample, 'truvari': truvari, 'refine': refine,
                      'enum_cap': enum_cap, 'labels_dir': labels_dir, 'apply_filtered': apply_filtered,
                      'overlap': overlap,
                      'work': os.path.join(STAGE3_WORK, 'score', label, rid)})
    done = []
    if jobs <= 1:
        for t in tasks:
            done.append(_batch_one(t))
            if not quiet:
                print('%s %s ED=%s ED_ref=%s' % done[-1][:4], flush=True)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as ex:
            for r in ex.map(_batch_one, tasks):
                done.append(r)
                if not quiet:
                    print('%s %s ED=%s ED_ref=%s' % r[:4], flush=True)
    for rid, vcf in missing:
        print('%s: no VCF at %s' % (rid, vcf), file=sys.stderr)
    return done, missing


# ----------------------------------------------------------------------------- summarise

def result_status(results_dir, label, keep):
    """{status class: [region ids]} of a label's JSONs over the kept regions (missing = no JSON)."""
    st = collections.defaultdict(list)
    seen = set()
    for p in glob.glob(os.path.join(results_dir, label, '*.json')):
        try:
            with open(p) as f:
                d = json.load(f)
        except ValueError:
            continue
        rid = d.get('region_id')
        if rid not in keep:
            continue
        seen.add(rid)
        st['ok' if d.get('status') == 'ok' else 'error'].append(rid)
    st['missing'] = sorted(set(keep) - seen)
    return st


def load_results(results_dir, label):
    out = {}
    for p in glob.glob(os.path.join(results_dir, label, '*.json')):
        try:
            with open(p) as f:
                d = json.load(f)
        except ValueError:
            continue
        if d.get('status') == 'ok':
            out[d['region_id']] = d
    return out


def _med(xs):
    xs = [x for x in xs if x is not None]
    return float(statistics.median(xs)) if xs else None


def _fmt(x, nd=1):
    if x is None:
        return '-'
    if isinstance(x, float):
        return ('%.' + str(nd) + 'f') % x
    return str(x)


def load_noise(path):
    """{graph: {region: noise}} from replicate_noise.py's TSV (the truth-free replicate spread of each
    graph's own call: max over its replicate arms of the pair ED to the Stage 3 call)."""
    out = collections.defaultdict(dict)
    if not path:
        return out
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row.get('noise') not in (None, '', 'None'):
                out[row['graph']][row['region_id']] = int(row['noise'])
    return out


def summarise(labels, baseline, results_dir=None, regions=None, out=None, title=None, noise=None):
    """Per-stratum paired tables. noise (replicate_noise.py's TSV) adds a second, noise-aware count
    beside the rule's: a region counts as better (worse) only when the label's ED is below (above)
    the baseline's by more than the two graphs' replicate noise summed -- by the triangle
    inequality, then no measured replicate of either graph could close the gap. The rule's own
    counts are unchanged."""
    results_dir = results_dir or STAGE3_RESULTS
    nz = load_noise(noise)
    rows = {r['region_id']: r for r in regions_table()}
    keep = set(select_regions(regions)) if regions else set(rows)
    data = {lab: {k: v for k, v in load_results(results_dir, lab).items() if k in keep} for lab in labels}
    if baseline not in data:
        data[baseline] = {k: v for k, v in load_results(results_dir, baseline).items() if k in keep}
    base = data[baseline]
    status = {lab: result_status(results_dir, lab, keep) for lab in set(labels) | {baseline}}
    table = []
    strata = [s for s in STRATA_ORDER if any(rows[r]['stratum'] == s for r in keep)]
    for stratum in strata:
        for lab in labels:
            d = data[lab]
            ids = sorted(r for r in d if rows[r]['stratum'] == stratum and r in base)
            if not ids:
                table.append({'stratum': stratum, 'label': lab, 'n': 0})
                continue
            ed = [d[r]['ed'] for r in ids]
            eb = [base[r]['ed'] for r in ids]
            er = [d[r]['ed_ref'] for r in ids]
            diff = [a - b for a, b in zip(ed, eb)]
            eraw = [d[r].get('sensitivity', {}).get('ed_raw_overlap', d[r]['ed']) for r in ids]
            ebraw = [base[r].get('sensitivity', {}).get('ed_raw_overlap', base[r]['ed']) for r in ids]
            draw = [a - b for a, b in zip(eraw, ebraw)]
            esup = [d[r].get('sensitivity', {}).get('ed_suppress_nested', d[r]['ed']) for r in ids]
            ebsup = [base[r].get('sensitivity', {}).get('ed_suppress_nested', base[r]['ed']) for r in ids]
            dsup = [a - b for a, b in zip(esup, ebsup)]
            if noise:
                margin = [None if (r not in nz.get(lab, {}) or r not in nz.get(baseline, {}))
                          else (0 if lab == baseline else nz[lab][r] + nz[baseline][r]) for r in ids]
                have = [(x, m) for x, m in zip(diff, margin) if m is not None]
                nstat = {'noise_n': len(have), 'noise_margin_sum': sum(m for _, m in have),
                         'better_noise': sum(1 for x, m in have if x < -m),
                         'worse_noise': sum(1 for x, m in have if x > m),
                         'noise_missing': ','.join(r for r, m in zip(ids, margin) if m is None)}
            else:
                nstat = {}
            tr = [d[r].get('truvari', {}) for r in ids]

            def pooled(prefix):
                keys = ('TP-base', 'FN', 'TP-comp', 'FP')
                vals = [[t.get('%s_%s' % (prefix, k)) for k in keys] for t in tr]
                vals = [v for v in vals if None not in v]
                if not vals:
                    return None, 0
                s = [sum(v[i] for v in vals) for i in range(4)]
                return f1(*s), len(vals)
            rf, rn = pooled('raw')
            ff, fn_ = pooled('refined')
            pf, pn = pooled('phab')
            fp_sum = sum(t.get('raw_FP') or 0 for t in tr)
            fn_sum = sum(t.get('raw_FN') or 0 for t in tr)
            flags = collections.Counter(f for r in ids for f in d[r].get('flags', []))
            table.append({
                'stratum': stratum, 'label': lab, 'baseline': baseline, 'n': len(ids),
                'ed_median': _med(ed), 'ed_sum': sum(ed), 'base_ed_median': _med(eb), 'base_ed_sum': sum(eb),
                'ed_per_kb_median': _med([d[r]['ed_per_kb'] for r in ids]),
                'ed_ref_sum': sum(er),
                'gain_pooled': round(1 - sum(ed) / sum(er), 4) if sum(er) else None,
                'gain_median': _med([d[r]['gain'] for r in ids]),
                'exact': sum(1 for r in ids if d[r]['exact']),
                'base_exact': sum(1 for r in ids if base[r]['exact']),
                'diff_median': _med(diff), 'diff_sum': sum(diff),
                'better_gt0': sum(1 for x in diff if x < 0), 'worse_gt0': sum(1 for x in diff if x > 0),
                'equal': sum(1 for x in diff if x == 0),
                'better_gt10': sum(1 for x in diff if x < -10), 'worse_gt10': sum(1 for x in diff if x > 10),
                'raw_f1_pooled': rf, 'refined_f1_pooled': ff, 'phab_f1_pooled': pf, 'truvari_n': rn,
                'raw_fp_sum': fp_sum, 'raw_fn_sum': fn_sum,
                'rawrule_ed_sum': sum(eraw), 'rawrule_base_ed_sum': sum(ebraw),
                'rawrule_better_gt0': sum(1 for x in draw if x < 0), 'rawrule_worse_gt0': sum(1 for x in draw if x > 0),
                'rule_changed_regions': sum(1 for r, a in zip(ids, eraw) if a != d[r]['ed']),
                'suppress_ed_sum': sum(esup), 'suppress_base_ed_sum': sum(ebsup),
                'suppress_better_gt0': sum(1 for x in dsup if x < 0), 'suppress_worse_gt0': sum(1 for x in dsup if x > 0),
                'suppress_changed_regions': sum(1 for r, a in zip(ids, esup) if a != d[r]['ed']),
                'applied_under_nonref_parent': sum(sum(d[r]['records'].get('applied_under_nonref_parent', [0, 0]))
                                                   for r in ids),
                'descendant_before_ancestor': sum(1 for r in ids
                                                  if d[r]['records'].get('descendant_before_ancestor_by_pos')),
                'phase_optimised': sum(1 for r in ids if d[r]['phase']['n_blocks'] > 1),
                'phase_greedy': sum(1 for r in ids if d[r]['phase']['enum_cap_hit']),
                'crossing_edge': sum(1 for r in ids if d[r]['records'].get('records_crossing_span_edge')),
                'missing_or_filtered': sum(1 for r in ids if d[r]['records'].get('missing_alleles_as_ref')
                                           or d[r]['records'].get('filtered_records_as_ref')),
                'regions': ','.join(ids),
                'top_flags': ';'.join('%s:%d' % kv for kv in flags.most_common(6)),
            })
            table[-1].update(nstat)
    # per-region table
    per_region = []
    for rid in sorted(keep, key=lambda r: (STRATA_ORDER.index(rows[r]['stratum']), r)):
        row = {'region_id': rid, 'stratum': rows[rid]['stratum']}
        anyv = False
        for lab in labels:
            d = data[lab].get(rid)
            row[lab] = d['ed'] if d else None
            anyv = anyv or d is not None
        refd = next((data[lab][rid] for lab in labels if rid in data[lab]), None)
        row['ed_ref'] = refd['ed_ref'] if refd else None
        row['truth_bp'] = sum(refd['truth_len']) if refd else None
        if anyv:
            per_region.append(row)
    if out:
        cols = ['stratum', 'label', 'baseline', 'n', 'ed_median', 'base_ed_median', 'ed_sum', 'base_ed_sum',
                'diff_median', 'diff_sum', 'better_gt0', 'worse_gt0', 'equal', 'better_gt10', 'worse_gt10',
                'exact', 'base_exact', 'ed_per_kb_median', 'ed_ref_sum', 'gain_pooled', 'gain_median',
                'raw_f1_pooled', 'refined_f1_pooled', 'phab_f1_pooled', 'truvari_n', 'raw_fp_sum', 'raw_fn_sum',
                'rawrule_ed_sum', 'rawrule_base_ed_sum', 'rawrule_better_gt0', 'rawrule_worse_gt0', 'rule_changed_regions',
                'suppress_ed_sum', 'suppress_base_ed_sum', 'suppress_better_gt0', 'suppress_worse_gt0',
                'suppress_changed_regions', 'applied_under_nonref_parent', 'descendant_before_ancestor',
                'phase_optimised', 'phase_greedy', 'crossing_edge', 'missing_or_filtered', 'top_flags', 'regions']
        if noise:
            cols += ['noise_n', 'noise_margin_sum', 'better_noise', 'worse_noise', 'noise_missing']
        with open(out + '.tsv', 'w') as f:
            f.write('\t'.join(cols) + '\n')
            for t in table:
                f.write('\t'.join(_fmt(t.get(c), 4) for c in cols) + '\n')
        with open(out + '.md', 'w') as f:
            f.write(summary_markdown(table, per_region, labels, baseline, results_dir, regions, title, status,
                                     noise=noise, nz=nz))
    return table, per_region


def summary_markdown(table, per_region, labels, baseline, results_dir, regions, title=None, status=None,
                     noise=None, nz=None):
    L = ['# %s' % (title or 'Stage 3: called haplotypes against HG002, paired per stratum'), '',
         'Written by `tools/score_haplotypes.py summarise` from `%s/<label>/<id>.json`; the scoring rule is in '
         'the docstring of `tools/score_haplotypes.py`. Regions: %s. '
         'Baseline: `%s`.' % (os.path.relpath(results_dir, config.REPO), regions or 'all', baseline), '',
         'ED is the summed unit edit distance between the called and the true haplotype pair over the '
         'anchor-to-anchor span, best pairing, 0 = both haplotypes exact. Differences are label minus '
         'baseline on the regions both have, so negative is better. "better/worse >0" count regions whose ED '
         'fell/rose at all, ">10" by more than 10. Gain = 1 - ED/ED_ref, where ED_ref calls CHM13 on both '
         'haplotypes; "pooled" sums ED and ED_ref over the stratum. F1s are pooled over the stratum (TP/FP/FN '
         'summed) and are representation scores only.', '']
    if status:
        for lab in labels + ([baseline] if baseline not in labels else []):
            st = status[lab]
            L.append('- `%s`: %d scored%s%s' % (
                lab, len(st['ok']),
                ('; errors: %s' % ', '.join(sorted(st['error']))) if st['error'] else '',
                ('; no result: %s' % ', '.join(st['missing'])) if st['missing'] else ''))
        L.append('')
    L += ['| stratum | label | n | median ED (base) | sum ED (base) | median diff | better / worse >0 | '
         'better / worse >10 | exact (base) | median ED/kb | gain pooled / median | raw / refined / phab F1 | '
         'FP / FN (raw) |', '|' + '---|' * 13]
    for t in table:
        if not t.get('n'):
            L.append('| %s | %s | 0 | | | | | | | | | | |' % (t['stratum'], t['label']))
            continue
        L.append('| %s | %s | %d | %s (%s) | %d (%d) | %s | %d / %d | %d / %d | %d (%d) | %s | %s / %s | %s / %s / %s | '
                 '%d / %d |' % (
                     t['stratum'], t['label'], t['n'], _fmt(t['ed_median']), _fmt(t['base_ed_median']), t['ed_sum'],
                     t['base_ed_sum'], _fmt(t['diff_median']), t['better_gt0'], t['worse_gt0'], t['better_gt10'],
                     t['worse_gt10'], t['exact'], t['base_exact'], _fmt(t['ed_per_kb_median'], 2),
                     _fmt(t['gain_pooled'], 3), _fmt(t['gain_median'], 3), _fmt(t['raw_f1_pooled'], 3),
                     _fmt(t['refined_f1_pooled'], 3), _fmt(t['phab_f1_pooled'], 3), t['raw_fp_sum'], t['raw_fn_sum']))
    L += ['', '## Sensitivity to the overlap rule', '',
          'The same tables under the raw overlap rule (an allele is skipped when its REF span shares a base with '
          'an allele already applied to that haplotype; bcftools consensus\'s and region.py\'s rule), which the '
          'scorer does not use because it drops real variants from the truth VCF itself '
          '(results/stage3_scorer_validation.tsv). "changed" counts regions whose ED the rule moves.', '',
          '| stratum | label | n | sum ED (base), raw rule | better / worse >0, raw rule | changed |',
          '|---|---|---|---|---|---|']
    for t in table:
        if t.get('n'):
            L.append('| %s | %s | %d | %d (%d) | %d / %d | %d |' % (
                t['stratum'], t['label'], t['n'], t['rawrule_ed_sum'], t['rawrule_base_ed_sum'],
                t['rawrule_better_gt0'], t['rawrule_worse_gt0'], t['rule_changed_regions']))
    L += ['', '## Sensitivity to the nesting clause', '',
          'The same tables when a parent\'s non-reference allele suppresses every allele nested in it on that slot '
          '(sensitivity.ed_suppress_nested), against the rule, under which a child applies where the parent\'s '
          'trimmed change does not cover it. "under parent" counts the child alleles the rule applies on a slot '
          'where a record they are nested in is non-reference (records.applied_under_nonref_parent); "reordered" '
          'counts regions where a vg child\'s POS lies left of its parent\'s, which scorer version 1 got wrong.', '',
          '| stratum | label | n | sum ED (base), suppress | better / worse >0, suppress | changed | under parent | reordered |',
          '|---|---|---|---|---|---|---|---|']
    for t in table:
        if t.get('n'):
            L.append('| %s | %s | %d | %d (%d) | %d / %d | %d | %d | %d |' % (
                t['stratum'], t['label'], t['n'], t['suppress_ed_sum'], t['suppress_base_ed_sum'],
                t['suppress_better_gt0'], t['suppress_worse_gt0'], t['suppress_changed_regions'],
                t['applied_under_nonref_parent'], t['descendant_before_ancestor']))
    if noise:
        L += ['', '## Noise-aware counts', '',
              'From `%s` (tools/replicate_noise.py): each graph\'s noise at a region is the largest pair ED between '
              'its Stage 3 call and one of its replicate calls (flank 50 kb; span node IDs re-laid; span nodes on mc\'s '
              'native IDs), computed without the truth. A region counts as better (worse) here only when the label\'s ED is below (above) '
              'the baseline\'s by more than the two graphs\' noise summed. The rule\'s counts above are unchanged. '
              '"n" counts the regions where both graphs have a noise measurement.' % os.path.relpath(noise, config.REPO),
              '', '| stratum | label | n | summed margin | better / worse beyond noise | better / worse >0 (rule) | no noise |',
              '|---|---|---|---|---|---|---|']
        for t in table:
            if t.get('n') and t['label'] != baseline:
                L.append('| %s | %s | %d | %d | %d / %d | %d / %d | %s |' % (
                    t['stratum'], t['label'], t.get('noise_n', 0), t.get('noise_margin_sum', 0),
                    t.get('better_noise', 0), t.get('worse_noise', 0), t['better_gt0'], t['worse_gt0'],
                    t.get('noise_missing') or '-'))
    L += ['', '## Mechanical reading of the decision rule', '',
          'Per label: hotspot ED falls = summed ED over the hotspot VNTRs is lower than the baseline\'s and more '
          'regions are better than worse (at >0); a control stratum is worse = its summed ED is higher or more '
          'regions are worse than better (at >0). This is only the arithmetic; the decision is made on the '
          'tables above.', '']
    for lab in labels:
        if lab == baseline:
            continue
        hs = [t for t in table if t['label'] == lab and t['stratum'] == 'hotspot_vntr' and t.get('n')]
        ctl = [t for t in table if t['label'] == lab and t['stratum'] in CONTROL_STRATA and t.get('n')]
        if not hs:
            L.append('- `%s`: no paired hotspot VNTRs.' % lab)
            continue
        h = hs[0]
        falls = h['ed_sum'] < h['base_ed_sum'] and h['better_gt0'] > h['worse_gt0']
        worse = [t['stratum'] for t in ctl if t['ed_sum'] > t['base_ed_sum'] or t['worse_gt0'] > t['better_gt0']]
        L.append('- `%s`: hotspot VNTR ED %s (sum %d vs %d, better/worse %d/%d); controls worse: %s.' % (
            lab, 'falls' if falls else 'does not fall', h['ed_sum'], h['base_ed_sum'], h['better_gt0'],
            h['worse_gt0'], ', '.join(worse) if worse else 'none'))
        if noise and h.get('noise_n'):
            worse_n = [t['stratum'] for t in ctl if t.get('worse_noise')]
            L.append('  - beyond noise: hotspot VNTRs better/worse %d/%d of %d; controls worse beyond noise: %s.' % (
                h['better_noise'], h['worse_noise'], h['noise_n'], ', '.join(worse_n) if worse_n else 'none'))
    L += ['', '## Per region (ED)', '']
    if noise:
        L += ['Each graph\'s replicate noise is in brackets.', '']
    L += ['| region | stratum | ' + ' | '.join(labels) + ' | ED_ref | truth bp |',
          '|' + '---|' * (len(labels) + 4)]
    nz = nz or {}

    def cell(lab, rid, v):
        if not noise or v is None:
            return _fmt(v)
        n = nz.get(lab, {}).get(rid)
        return '%s [%s]' % (_fmt(v), '-' if n is None else n)
    for r in per_region:
        L.append('| %s | %s | %s | %s | %s |' % (r['region_id'], r['stratum'],
                                                  ' | '.join(cell(lab, r['region_id'], r[lab]) for lab in labels),
                                                  _fmt(r['ed_ref']), _fmt(r['truth_bp'])))
    return '\n'.join(L) + '\n'


# ----------------------------------------------------------------------------- validation

def _write_vcf(path, header, lines):
    with open(path, 'w') as f:
        f.write('\n'.join(header) + '\n')
        for x in lines:
            f.write(x.rstrip('\n') + '\n')


def build_truth_vcf(regions, od, pad=5000):
    """The truth over the spans +- pad: (stvar, merged). stvar is the VCF truth.fa was built from
    (region.py applied stvar first and smvar only where it neither duplicated nor overlapped an
    applied stvar record: 1 smvar record in all 149 spans). merged is stvar plus the smvar records
    that are not exact duplicates (same POS, REF, ALT): the same haplotypes written redundantly,
    with smvar's multi-allelic packing, a harder input for the overlap rule."""
    rows = {r['region_id']: r for r in regions_table()}
    os.makedirs(od, exist_ok=True)
    bed = os.path.join(od, 'spans.bed')
    iv = sorted((rows[r]['contig'], max(0, int(rows[r]['span_start']) - 1 - pad), int(rows[r]['span_end']) + pad)
                for r in regions)
    with open(bed, 'w') as f:
        for c, a, b in iv:
            f.write('%s\t%d\t%d\n' % (c, a, b))
    p = config.data_paths()
    parts = []
    for src in ('stvar', 'smvar'):
        o = os.path.join(od, 'truth.%s.vcf.gz' % src)
        sh([config.BCFTOOLS, 'view', '-R', bed, '-Oz', '-o', o, p[src]])
        sh([config.BCFTOOLS, 'index', '-t', '-f', o])
        parts.append(o)
    cat = os.path.join(od, 'truth.concat.vcf.gz')
    sh([config.BCFTOOLS, 'concat', '-a', '-Oz', '-o', cat] + parts)
    srt = os.path.join(od, 'truth.sorted.vcf.gz')
    sh([config.BCFTOOLS, 'sort', '-Oz', '-o', srt, cat])
    out = os.path.join(od, 'truth.merged.vcf.gz')
    sh([config.BCFTOOLS, 'norm', '-d', 'exact', '-Oz', '-o', out, srt])
    sh([config.BCFTOOLS, 'index', '-t', '-f', out])
    for x in (cat, srt):
        os.remove(x)
    return parts[0], out


def perturb(region, recs_lines, header, kind):
    """A perturbed copy of the truth records of one region, or None when the region offers no
    record to perturb. Returns (lines, description, expectation) with expectation 'gt0' or 'eq0'."""
    a1, b1 = region['a1'], region['b1']
    parsed = [(ln, parse_record(ln, 0, i)) for i, ln in enumerate(recs_lines)]

    def ok(r):
        return r['filter'] in PASS_FILTERS and r['pos'] >= a1 and r['end'] <= b1 and \
            all(re.fullmatch(r'[ACGTN]+', a) for a in r['alts'])

    def eff(r):
        g = [0 if a is None else a for a in r['gt']]
        return g * 2 if len(g) == 1 else g

    def set_gt(line, gt_str, ps=None):
        f = line.rstrip('\n').split('\t')
        keys = f[8].split(':')
        vals = f[9].split(':')
        vals[keys.index('GT')] = gt_str
        if ps is not None:
            if 'PS' in keys:
                vals[keys.index('PS')] = str(ps)
            else:
                keys.append('PS')
                vals.append(str(ps))
        f[8], f[9] = ':'.join(keys), ':'.join(vals)
        return '\t'.join(f)

    def is_het(r):
        g = eff(r)
        return len(r['gt']) == 2 and g[0] != g[1]

    def applied_nonref(r):
        return ok(r) and any(eff(r))

    if kind == 'drop_sv':
        cand = [(ln, r) for ln, r in parsed if applied_nonref(r) and
                max(abs(len(a) - len(r['ref'])) for a in r['alts']) >= SV]
        if not cand:
            return None
        ln0, r0 = max(cand, key=lambda x: max(abs(len(a) - len(x[1]['ref'])) for a in x[1]['alts']))
        return [ln for ln, r in parsed if ln is not ln0], 'dropped SV record at %d (%s)' % (
            r0['pos'], '/'.join(str(len(a) - len(r0['ref'])) for a in r0['alts'])), 'gt0'
    hets = [(ln, r) for ln, r in parsed if applied_nonref(r) and is_het(r)]
    if kind in ('swap_phase', 'unphase_swapped'):
        # needs two het records with distinct allele pairs, else a swap is undetectable by design
        if len(hets) < 2:
            return None
        ln0, r0 = hets[len(hets) // 2]
        g = r0['gt']
        sw = '%s|%s' % tuple('.' if a is None else str(a) for a in (g[1], g[0]))
        if kind == 'unphase_swapped':
            sw = sw.replace('|', '/')
        new = [set_gt(ln, sw) if ln is ln0 else ln for ln, r in parsed]
        return new, '%s record at %d: %s -> %s' % ('swapped' if kind == 'swap_phase' else 'swapped and unphased',
                                                   r0['pos'], '|'.join('.' if a is None else str(a) for a in g),
                                                   sw), ('gt0' if kind == 'swap_phase' else 'eq0')
    if kind == 'ps_split_swapped':
        if len(hets) < 2:
            return None
        cut = hets[len(hets) // 2][1]['pos']
        new = []
        for ln, r in parsed:
            if len(r['gt']) == 2 and '|' in ln.split('\t')[9].split(':')[0]:
                if r['pos'] >= cut:
                    g = r['gt']
                    new.append(set_gt(ln, '%s|%s' % tuple('.' if a is None else str(a) for a in (g[1], g[0])), ps=2))
                else:
                    new.append(set_gt(ln, ln.split('\t')[9].split(':')[0], ps=1))
            else:
                new.append(ln)
        return new, 'records from %d in PS 2 with every GT swapped, the rest PS 1' % cut, 'eq0'
    if kind == 'flip_gt':
        # a change that leaves the VCF consistent: a hom-alt record made heterozygous, or else a
        # heterozygous record made hom-alt where the other haplotype has no overlapping allele
        # (otherwise the new allele conflicts with one already there, is skipped, and nothing changes)
        def size(r):
            return max(abs(len(a) - len(r['ref'])) for a in r['alts']), len(r['ref'])
        homs = [(ln, r) for ln, r in parsed if applied_nonref(r) and len(r['gt']) == 2 and eff(r)[0] == eff(r)[1]]
        if homs:
            ln0, r0 = max(homs, key=lambda x: size(x[1]))
            g = eff(r0)
            new_gt = '%d|0' % g[0]
        else:
            cand = []
            for ln, r in hets:
                g = eff(r)
                empty = 0 if not g[0] else 1
                clash = any(q is not r and len(q['gt']) == 2 and eff(q)[empty] and q['pos'] <= r['end'] and
                            r['pos'] <= q['end'] for _, q in parsed)
                if not clash:
                    cand.append((ln, r))
            if not cand:
                return None
            ln0, r0 = max(cand, key=lambda x: size(x[1]))
            g = eff(r0)
            a = g[0] or g[1]
            new_gt = '%d|%d' % (a, a)
        new = [set_gt(ln, new_gt) if ln is ln0 else ln for ln, r in parsed]
        return new, 'GT at %d: %s -> %s' % (r0['pos'], '|'.join(str(x) for x in eff(r0)), new_gt), 'gt0'
    raise ValueError(kind)


EXPECT = {'truth_stvar': 'eq0', 'truth_stvar_raw': 'eq0', 'truth_stvar_apply_filtered': 'eq0',
          'truth_merged': 'eq0', 'truth_merged_raw': 'eq0', 'empty_vcf': 'eq_ref',
          'drop_sv': 'gt0', 'swap_phase': 'gt0', 'flip_gt': 'gt0',
          'unphase_swapped': 'eq0', 'ps_split_swapped': 'eq0', 'ps_split_swapped_greedy': 'eq0'}


def validate(jobs=4, regions=None, out=None, truvari=True):
    regions = regions or select_regions('all')
    od = os.path.join(STAGE3_WORK, 'validation')
    os.makedirs(od, exist_ok=True)
    rows = {r['region_id']: r for r in regions_table()}
    t0 = time.time()
    stvar_vcf, merged_vcf = build_truth_vcf(regions, os.path.join(od, 'truth'))
    header = vcf_header_lines(stvar_vcf)
    res_dir = os.path.join(od, 'results')
    report = collections.OrderedDict()
    descr = collections.defaultdict(dict)
    # (1) the truth VCF itself: the rule as fixed (trimmed overlap, FILTERed records as reference),
    # the raw overlap rule for comparison, FILTERed records applied, and the redundant merged VCF
    for lab, vcf, kw in (('truth_stvar', stvar_vcf, {'truvari': truvari}),
                         ('truth_stvar_raw', stvar_vcf, {'overlap': 'raw'}),
                         ('truth_stvar_apply_filtered', stvar_vcf, {'apply_filtered': True}),
                         ('truth_merged', merged_vcf, {}),
                         ('truth_merged_raw', merged_vcf, {'overlap': 'raw'})):
        kw.setdefault('truvari', False)
        batch(lab, vcf, regions, out_dir=os.path.join(res_dir, lab), jobs=jobs, force=True, refine=False,
              quiet=True, **kw)
        report[lab] = load_results(res_dir, lab)
    # (3) an empty VCF
    empty = os.path.join(od, 'empty.vcf')
    contigs = sorted(set(rows[r]['contig'] for r in regions))
    hdr = ['##fileformat=VCFv4.2'] + ['##contig=<ID=%s>' % c for c in contigs] + \
          ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
           '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002']
    _write_vcf(empty, hdr, [])
    batch('empty_vcf', empty, regions, out_dir=os.path.join(res_dir, 'empty_vcf'), jobs=jobs, force=True,
          truvari=False, quiet=True)
    report['empty_vcf'] = load_results(res_dir, 'empty_vcf')
    # (2) perturbations of the stvar truth, one VCF per region and kind
    for kind in ('drop_sv', 'swap_phase', 'flip_gt', 'unphase_swapped', 'ps_split_swapped'):
        kd = os.path.join(od, 'perturbed', kind)
        os.makedirs(kd, exist_ok=True)
        have = []
        for rid in regions:
            region = load_region(os.path.join(config.REGIONS_DIR, rid))
            p = run([config.BCFTOOLS, 'view', '-H', '-r', '%s:%d-%d' % (region['contig'], region['a1'], region['b1']),
                     stvar_vcf])
            lines = [x for x in p.stdout.splitlines() if x]
            pr = perturb(region, lines, header, kind)
            path = os.path.join(kd, rid + '.vcf')
            if pr is None:
                if os.path.exists(path):
                    os.remove(path)
                continue
            _write_vcf(path, header, pr[0])
            descr[kind][rid] = pr[1]
            have.append(rid)
        batch(kind, os.path.join(kd, '{id}.vcf'), have, out_dir=os.path.join(res_dir, kind), jobs=jobs,
              force=True, truvari=False, quiet=True)
        report[kind] = load_results(res_dir, kind)
    # the phase search forced into its greedy branch on the PS-split case
    batch('ps_split_swapped_greedy', os.path.join(od, 'perturbed', 'ps_split_swapped', '{id}.vcf'),
          list(descr['ps_split_swapped']), out_dir=os.path.join(res_dir, 'ps_split_swapped_greedy'), jobs=jobs,
          force=True, truvari=False, enum_cap=0, quiet=True)
    report['ps_split_swapped_greedy'] = load_results(res_dir, 'ps_split_swapped_greedy')
    descr['ps_split_swapped_greedy'] = descr['ps_split_swapped']
    # verdicts
    lines_tsv = []
    summary = collections.OrderedDict()
    for lab, d in report.items():
        exp = EXPECT[lab]
        for rid in regions:
            if lab in descr and rid not in descr[lab]:
                continue
            r = d.get(rid)
            if r is None:
                verdict, ed, edr = 'MISSING', None, None
            else:
                ed, edr = r['ed'], r['ed_ref']
                ok = {'eq0': ed == 0, 'eq_ref': ed == edr, 'gt0': ed > 0}[exp]
                verdict = 'pass' if ok else 'FAIL'
            tv = (r or {}).get('truvari', {})
            rec = (r or {}).get('records', {})
            lines_tsv.append([lab, exp, rid, rows[rid]['stratum'], verdict, ed, edr,
                              ';'.join((r or {}).get('flags', [])), rec.get('filtered_records_as_ref', 0),
                              rec.get('missing_alleles_as_ref', 0), descr.get(lab, {}).get(rid, ''),
                              tv.get('raw_f1') if tv.get('status') == 'ok' else ''])
            summary.setdefault(lab, collections.Counter())[verdict] += 1
    out = out or os.path.join(config.RESULTS_DIR, 'stage3_scorer_validation')
    with open(out + '.tsv', 'w') as f:
        f.write('check\texpect\tregion_id\tstratum\tverdict\ted\ted_ref\tflags\tfiltered_as_ref\t'
                'missing_alleles\tperturbation\ttruvari_raw_f1\n')
        for x in lines_tsv:
            f.write('\t'.join('' if v is None else str(v) for v in x) + '\n')
    return summary, lines_tsv, time.time() - t0


# ----------------------------------------------------------------------------- CLI

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('score', help='one region, one VCF')
    s.add_argument('region_dir')
    s.add_argument('vcf')
    s.add_argument('--label')
    s.add_argument('--out', help='JSON path (default: stdout)')
    s.add_argument('--haplotypes', help='also write the two called haplotypes as FASTA')
    b = sub.add_parser('batch', help='many regions, one label')
    b.add_argument('--label', required=True)
    b.add_argument('--vcf', required=True, help='path pattern with {id}, {contig}, {stratum}')
    b.add_argument('--regions', default='all', help="'all', 'pilot', a comma list or a file of ids")
    b.add_argument('--out-dir', help='default results/stage3/<label>')
    b.add_argument('--jobs', type=int, default=4)
    b.add_argument('--force', action='store_true')
    b.add_argument('--missing-ok', action='store_true', help='skip regions without a VCF silently')
    for p in (s, b):
        p.add_argument('--sample', default='HG002')
        p.add_argument('--no-truvari', action='store_true')
        p.add_argument('--no-refine', action='store_true')
        p.add_argument('--refine-timeout', type=int, default=900)
        p.add_argument('--enum-cap', type=int, default=ENUM_CAP)
        p.add_argument('--apply-filtered', action='store_true',
                       help='apply FILTERed records (the rule treats them as reference; for checks only)')
        p.add_argument('--overlap', choices=('trimmed', 'raw'), default=None,
                       help='overlap rule within a GT slot (default %s; raw is for checks only)' % OVERLAP_MODE)
        p.add_argument('--labels-dir', help='genome-wide truvari label directory pattern with {contig}; '
                                            "'genomewide' = %s" % GENOMEWIDE_LABELS)
    m = sub.add_parser('summarise', help='per-stratum paired tables')
    m.add_argument('--labels', required=True)
    m.add_argument('--baseline', default='mc')
    m.add_argument('--results', default=STAGE3_RESULTS)
    m.add_argument('--regions', default=None)
    m.add_argument('--out', default=os.path.join(config.RESULTS_DIR, 'stage3_summary'))
    m.add_argument('--title')
    m.add_argument('--noise', help="replicate_noise.py's TSV: add noise-aware counts beside the rule's")
    v = sub.add_parser('validate', help="the scorer's own checks")
    v.add_argument('--jobs', type=int, default=4)
    v.add_argument('--regions', default='all')
    v.add_argument('--no-truvari', action='store_true')
    args = ap.parse_args(argv)

    if args.cmd in ('score', 'batch'):
        labels_dir = GENOMEWIDE_LABELS if args.labels_dir == 'genomewide' else args.labels_dir
    if args.cmd == 'score':
        res = score_region(args.region_dir, args.vcf, label=args.label, sample=args.sample,
                           truvari=not args.no_truvari, refine=not args.no_refine, enum_cap=args.enum_cap,
                           labels_dir=labels_dir, apply_filtered=args.apply_filtered,
                           write_haplotypes=args.haplotypes, refine_timeout=args.refine_timeout,
                           overlap=args.overlap)
        txt = json.dumps(res, indent=1)
        if args.out:
            with open(args.out, 'w') as f:
                f.write(txt + '\n')
        else:
            print(txt)
    elif args.cmd == 'batch':
        vcf = GENOMEWIDE_VCF if args.vcf == 'genomewide' else args.vcf
        done, missing = batch(args.label, vcf, select_regions(args.regions), out_dir=args.out_dir, jobs=args.jobs,
                              force=args.force, truvari=not args.no_truvari, refine=not args.no_refine,
                              enum_cap=args.enum_cap, labels_dir=labels_dir, sample=args.sample,
                              apply_filtered=args.apply_filtered, missing_ok=args.missing_ok, overlap=args.overlap)
        bad = [d for d in done if d[1] != 'ok']
        print('%d scored, %d errors, %d without a VCF' % (len(done), len(bad), len(missing)))
    elif args.cmd == 'summarise':
        labels = [x for x in args.labels.split(',') if x]
        table, _ = summarise(labels, args.baseline, results_dir=args.results, regions=args.regions, out=args.out,
                             title=args.title, noise=args.noise)
        print(open(args.out + '.md').read())
    elif args.cmd == 'validate':
        summary, _, secs = validate(jobs=args.jobs, regions=select_regions(args.regions),
                                    truvari=not args.no_truvari)
        for lab, c in summary.items():
            print('%-28s %s' % (lab, dict(c)))
        print('%.0f s' % secs)


if __name__ == '__main__':
    main()
