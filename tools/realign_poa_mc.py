#!/usr/bin/env python3
"""realign_poa_mc.py -- poa_abpoa_mc: abPOA run the way Minigraph-Cactus runs it in its BAR phase.

One method, `poa_abpoa_mc`: every setting is Cactus's (pangenome mode), not abPOA's CLI default, and
long regions are cut the way BAR cuts them. It is a realigner like `poa_abpoa` (tools/realign_poa.py)
and shares realign.py's pipeline (N masking, output checks, msa_graph.py graph, candidate files); it
differs in the aligner call and in two details listed under "Around the aligner".

Region mode: align regions/<id>/hap32.fa, build the graph, write
candidates/poa_abpoa_mc/<id>.{msa.fa,gfa,realign.json} and one row per region in
results/mcpoa_runtime.tsv (not realign_runtime.tsv):

    python3 tools/realign_poa_mc.py regions/L012184 [--out CANDIDATES_ROOT]
    python3 tools/realign_poa_mc.py all --jobs 2 [--stratum hotspot_vntr] [--force | --retry-failed]

MSA-only mode, for any FASTA whose records all run anchor to anchor (e.g. a region's full panel);
--region (or --flanks L,R) gives the anchor lengths, without them the flanks are 0 bp:

    python3 tools/realign_poa_mc.py --input seqs.fa[.gz] --msa-out out.msa.fa --region L014297 [--json info.json]
        [--project-to regions/<id>/hap32.fa --projected-out hap32rows.msa.fa]

    import sys; sys.path.insert(0, 'tools')
    import realign_poa_mc
    info = realign_poa_mc.align_fasta('seqs.fa', 'out.msa.fa', flanks=(1, 1))   # info['status']

What Cactus does, and where it is written (Cactus checkout ~/CLionProjects/cactus, commit c07e5b4):

  * Scores and band (src/cactus/cactus_progressive_config.xml, <bar><poa>, read by
    abpoaParamaters_constructFromCactusParams in bar/impl/poaBarAligner.c): global mode; adaptive band
    b=1000, f=0.1; a 5x5 ACGTN substitution matrix (a HOXD70 variant: transitions -61, N -100 against
    every base, N/N 100), not match/mismatch; convex gaps O=400,1200 E=30,1; minimizer seeding off;
    progressive (guide-tree) order on, with k=15, w=5, min_w=500; progressive is dropped above 5000
    rows (and when 1 - shortest/longest > 1.0, which never happens).
  * The command line: Cactus's own debug dump (dump_abpoa_input) prints
    `abpoa in.fa -O 400,1200 -E 30,1 -b 1000 -f 0.1 -t matrix.mtx -r 1 -m 0 -p`, the matrix written
    tab-separated under the header "\\tA\\tC\\tG\\tT\\tN". It omits -k/-w/-n because seeding is off, but
    abPOA's -p guide tree sketches minimizers with abpt->k and abpt->w whether or not seeding is on
    (abpoa_build_guide_tree_partition in submodules/abPOA/src/abpoa_seed.c), and the library call sets
    k=15, w=5 where the CLI defaults are 19 and 10. So `-k 15 -w 5 -n 500` are passed here too.
  * The binary: Cactus's vendored abPOA 1.5.4 (submodules/abPOA/bin/abpoa), not the 1.5.7 poa_abpoa uses.
  * Input order: longest first (get_end_sequences sorts caps by adjacency length, descending; ties
    here keep input order, where Cactus's qsort leaves them unspecified). With -p abPOA then
    re-orders by its guide tree, so the sort mostly breaks ties.
  * Length rule (cactus_align.py --pangenome sets bar bandingLimit = 10000; make_flower_alignment_poa):
    the region is treated as one BAR flower whose two ends are the anchor nodes, with one adjacency
    string per haplotype, the sequence strictly between the anchor nodes.
      - every string < 10 kb: one abPOA call on the whole strings (the "dominant end" branch; its
        strand is the left anchor's here, where Cactus's is whichever end its iterator meets first);
      - otherwise each end is aligned separately (make_consistent_partial_order_alignments): at the
        left anchor the first min(len, 10 kb) bases of each string, at the right anchor the first
        min(len, 10 kb) bases of each reverse-complemented string. A string shorter than 20 kb is in
        both MSAs over an overlap of 2 x min(len, 10 kb) - len bases, and trim() cuts it once, at the
        point that keeps the largest sum of column scores (column score = bases in the column - 1;
        trim, trim_msa_suffix, sum_column_scores and make_column_scores are ported line for line,
        including the order rows are trimmed in and the tie-breaks). The middle of a string > 20 kb
        (len - 20 kb) is in no MSA: it stays unaligned, one private node per haplotype.
    The anchors' bases are added back as shared columns at both ends. The per-region record says
    whether the rule fired, how many rows were >= 10 kb and > 20 kb and how many bp stayed unaligned.
    The sliding window (partialOrderAlignmentWindow, 10 kb) never engages, since no string passed to
    abPOA is longer than bandingLimit.
    CAVEAT: this is BAR alone on the whole region. In the real pipeline the CAF phase first pinches
    blocks from the minigraph mappings inside the region, so BAR sees many smaller flowers, and the
    10 kb rule rarely leaves anything unaligned there.

Around the aligner (realign.align_fasta, the shared pipeline): N runs are cut out before alignment
and put back as insertion columns of their own (as for every realigner; the length rule is applied to
the masked strings, which is what abPOA sees); each output row is checked against its input. Two
differences from poa_abpoa:
  * identical sequences are NOT collapsed (dedup off): Cactus aligns every haplotype as a row, and
    -p's guide tree (minimizer Jaccard over all rows) and BAR's column scores both count duplicates;
  * the memory predictor is this module's (predict_call_mb): abPOA allocates 5 int32 DP matrices of
    (graph nodes) x (query length), rounded up to a power of two, and Cactus's scores force int32
    cells at any length (max_mat=100). Measured peaks are in results/mcpoa_runtime.tsv.

align_fasta(..., region_timeout=S) caps the summed wall clock of a region's abPOA calls (the second end
call gets what the first left); the full-panel driver uses it for its per-region cap. Each call's record
(info['aligner']['calls']) carries abPOA's graph size (`nodes`: distinct bases per MSA column, summed).

Region mode runs regions smallest first, --jobs at a time (threads; each abPOA call is a separate
single-threaded process under --timeout seconds and a sampled-RSS cap of --mem-mb, and every call is
wrapped in /usr/bin/time -l for its exact peak RSS). A region predicted over --mem-mb is not run
(status memout, "predicted" in the note); one predicted over 4 GB first checks memory_pressure and
waits for that much free memory.

Not in realign.py's METHODS table on purpose: the BAR emulation needs each region's anchor lengths,
which realign.py's plugin interface does not pass.
"""
import argparse
import collections
import concurrent.futures
import csv
import datetime
import fcntl
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402
import msa_graph  # noqa: E402

METHOD = 'poa_abpoa_mc'
CACTUS_DIR = os.path.expanduser(os.environ.get('VNTR_CACTUS', '~/CLionProjects/cactus'))
ABPOA_MC = os.path.expanduser(os.environ.get('VNTR_ABPOA_MC') or
                              os.path.join(CACTUS_DIR, 'submodules/abPOA/bin/abpoa'))
TIME = '/usr/bin/time'

# <bar><poa> of src/cactus/cactus_progressive_config.xml (Cactus c07e5b4); pangenome mode sets
# partialOrderAlignmentDisableSeeding=1 (already the default) and bar bandingLimit=10000.
MC_POA = collections.OrderedDict([
    ('align_mode', 0),                  # ABPOA_GLOBAL_MODE ("only global works")
    ('band_constant', 1000),            # partialOrderAlignmentBandConstant  -> abpt->wb
    ('band_fraction', 0.1),             # partialOrderAlignmentBandFraction  -> abpt->wf
    ('gap_open1', 400), ('gap_ext1', 30),
    ('gap_open2', 1200), ('gap_ext2', 1),
    ('disable_seeding', 1),
    ('k', 15), ('w', 5), ('min_w', 500),
    ('progressive', 1),
    ('progressive_max_rows', 5000),
    ('progressive_max_length_diff', 1.0),
    ('window', 10000),                  # partialOrderAlignmentWindow (never engaged, see docstring)
    ('mask_filter', -1),                # partialOrderAlignmentMaskFilter (off)
])
# partialOrderAlignmentSubMatrix, row-major over A C G T N
SUB_MATRIX = [91, -114, -61, -123, -100,
              -114, 100, -125, -61, -100,
              -61, -125, 100, -114, -100,
              -123, -61, -114, 91, -100,
              -100, -100, -100, -100, 100]
BANDING_LIMIT = 10000                   # cactus_align.py: --pangenome -> maxLen 10000 -> bar bandingLimit

DEFAULT_TIMEOUT = 1800
DEFAULT_MEM_MB = 10000
WAIT_ABOVE_MB = 4000                    # check memory_pressure before a region predicted above this
RUNTIME_TSV = os.path.join(config.RESULTS_DIR, 'mcpoa_runtime.tsv')
RUNTIME_COLS = ['method', 'region_id', 'stratum', 'span_bp', 'n_seqs', 'n_distinct', 'max_len', 'max_interior',
                'total_bp', 'status', 'bar_rule', 'rows_ge_10kb', 'rows_gt_20kb', 'unaligned_bp',
                'abpoa_calls', 'progressive', 'align_s', 'abpoa_s', 'graph_s', 'total_s', 'peak_rss_mb',
                'sampled_rss_mb', 'predicted_mb', 'waited_s', 'timeout_s', 'mem_cap_mb', 'columns', 'gap_frac',
                'nodes', 'nodes_per_kb', 'frac_nodes_1bp', 'masked_runs', 'masked_bp', 'finished', 'note']
DESCRIPTION = ('abPOA 1.5.4 as Minigraph-Cactus BAR runs it (pangenome): global, -b 1000 -f 0.1, Cactus 5x5 '
               'matrix, -O 400,1200 -E 30,1, -p (k15 w5), longest first, 10 kb end rule, no dedup')

_COMP = str.maketrans('ACGTNacgtn-', 'TGCANtgcan-')


class AlignerFailed(Exception):
    def __init__(self, status, message, record=None):
        Exception.__init__(self, message)
        self.status, self.message, self.record = status, message, record


def revcomp(s):
    return s.translate(_COMP)[::-1]


def tool_version():
    try:
        p = subprocess.run([ABPOA_MC, '-v'], capture_output=True, text=True, timeout=30)
        return 'abpoa %s (%s)' % ((p.stdout + p.stderr).strip().splitlines()[0], ABPOA_MC)
    except (OSError, subprocess.SubprocessError, IndexError):
        return 'abpoa unknown (%s)' % ABPOA_MC


# ---------------------------------------------------------------- the abPOA call

def write_matrix(path):
    """The substitution matrix file, exactly as Cactus's dump_abpoa_input writes it."""
    with open(path, 'w') as f:
        f.write('\tA\tC\tG\tT\tN\n')
        for i in range(5):
            f.write('ACGTN'[i] + ''.join('\t%d' % SUB_MATRIX[i * 5 + j] for j in range(5)) + '\n')


def abpoa_command(in_fa, matrix, progressive=True):
    """Cactus's dump_abpoa_input command line, plus the k/w/min_w the library sets (see docstring)."""
    p = MC_POA
    cmd = [ABPOA_MC, in_fa, '-O', '%d,%d' % (p['gap_open1'], p['gap_open2']),
           '-E', '%d,%d' % (p['gap_ext1'], p['gap_ext2']), '-b', str(p['band_constant']),
           '-f', '%g' % p['band_fraction'], '-t', matrix, '-r', '1', '-m', str(p['align_mode'])]
    if not p['disable_seeding']:
        cmd.append('-S')
    cmd += ['-k', str(p['k']), '-w', str(p['w']), '-n', str(p['min_w'])]
    if progressive:
        cmd.append('-p')
    return cmd


def use_progressive(lens):
    """Cactus's test in msa_make_partial_order_alignment (lens in input order, longest first)."""
    p = MC_POA
    if not p['progressive']:
        return False
    if len(lens) > p['progressive_max_rows']:
        return False
    if lens and lens[0] > 0 and 1.0 - float(lens[-1]) / lens[0] > p['progressive_max_length_diff']:
        return False
    return True


def predict_call_mb(lens, nodes_per_bp=None, touched=None):
    """Predicted peak RSS (MB) of one abPOA call on strings of these lengths: 5 int32 matrices of
    (graph nodes) x (longest query), the allocation rounded up to a power of two, of which the band
    touches a fraction; the graph is taken as nodes_per_bp x the longest string. Calibrated on
    measured peaks (results/mcpoa_runtime.tsv), meant to refuse what cannot fit and to flag big jobs.
    The defaults are the hap32 calibration; the full-panel driver (tools/poa_panel_mc.py) passes its
    own."""
    if len(lens) < 2:
        return 0.0
    npb = NODES_PER_BP if nodes_per_bp is None else nodes_per_bp
    tch = TOUCHED if touched is None else touched
    l1 = max(lens)
    alloc = 20.0 * (l1 + 4) * (npb * l1 + 2)
    p2 = 1 << max(0, int(alloc - 1).bit_length())
    return round(60 + tch * p2 / 2 ** 20, 1)


# Calibration, /usr/bin/time -l on the 151 abPOA calls of the hap32 run with strings >= 1.5 kb
# (32-34 rows): MSA columns were 1.00-1.42 x the longest string (median 1.02; the graph has a few more
# nodes than columns), and the peak was 0.16-1.49 x (median 0.61) the power-of-two allocation computed
# with 1.5 nodes per bp; TOUCHED is 1.1 x that maximum. (The hap32 run used 1.4, fitted on a 9-call
# pilot; it under-predicted L000289 and L015904 by 4% and 2%, both far below the 4 GB flag, so no
# decision changed. Max measured peak 5.5 GB, L010328.) A full panel (hundreds of rows) grows a larger
# graph per bp: tools/poa_panel_mc.py passes its own constants (PANEL_NODES_PER_BP, PANEL_TOUCHED).
NODES_PER_BP = 1.5
TOUCHED = 1.65


def _peak_from_time(path):
    try:
        for line in open(path):
            if 'maximum resident set size' in line:
                return round(int(line.split()[0]) / 2.0 ** 20, 1)   # bytes on macOS
    except (OSError, ValueError):
        pass
    return None


def graph_nodes(rows):
    """abPOA's graph size behind an MSA it wrote: one node per distinct base in each column (the
    columns of -r 1 output are abPOA's aligned-node groups, one node per base in a group)."""
    return sum(len(set(col)) - (1 if '-' in col else 0) for col in zip(*rows))


def poa_msa(strings, workdir, tag, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB, calls=None):
    """One window of msa_make_partial_order_alignment: align `strings` (in the order given; the
    caller sorts them longest first) with one abPOA call; returns the aligned rows, same order.
    One string: returned as is. An empty string is aligned as one 'N' and masked out afterwards
    (Cactus's hack for abPOA, which cannot take empty sequences). Appends a record to `calls`."""
    import realign                              # late import (realign imports tools/realign_*.py)
    n = len(strings)
    if n == 0:
        return []
    if n == 1:
        return [strings[0]]
    if max(len(s) for s in strings) > MC_POA['window']:
        raise ValueError('string longer than the POA window (%d); BAR would slide a window' % MC_POA['window'])
    empty = [len(s) == 0 for s in strings]
    inp = ['N' if e else s for s, e in zip(strings, empty)]
    lens = [len(s) for s in inp]
    prog = use_progressive(lens)
    d = tempfile.mkdtemp(prefix='abpoa_mc.%s.' % tag, dir=workdir)
    rec = collections.OrderedDict([('tag', tag), ('n', n), ('max_len', max(lens)), ('total_bp', sum(lens)),
                                   ('empty', sum(empty)), ('progressive', prog),
                                   ('predicted_mb', predict_call_mb(lens))])
    try:
        in_fa = os.path.join(d, 'in.fa')
        with open(in_fa, 'w') as f:
            for i, s in enumerate(inp):
                f.write('>r%d\n%s\n' % (i, s))
        mat = os.path.join(d, 'matrix.mtx')
        write_matrix(mat)
        cmd = abpoa_command(in_fa, mat, prog)
        rec['command'] = ['abpoa'] + [('<in.fa>' if c == in_fa else '<matrix.mtx>' if c == mat else c)
                                      for c in cmd[1:]]
        tfile = os.path.join(d, 'time.txt')
        out = os.path.join(d, 'out.msa')
        log = os.path.join(d, 'abpoa.log')
        r = realign.run_proc([TIME, '-l', '-o', tfile] + cmd, out, log, timeout=timeout, mem_mb=mem_mb,
                             env=config.tool_env(), cwd=d)
        rec['status'] = r['status']
        rec['seconds'] = r['seconds']
        rec['sampled_rss_mb'] = r['peak_rss_mb']
        rec['peak_rss_mb'] = _peak_from_time(tfile)
        if calls is not None:
            calls.append(rec)
        if r['status'] == 'ok' and os.path.getsize(out) == 0:
            rec['status'] = 'error'
        if rec['status'] != 'ok':
            raise AlignerFailed(rec['status'], '%s after %.0f s (peak sampled RSS %.0f MB)%s' % (
                rec['status'], r['seconds'], r['peak_rss_mb'] or 0,
                ('; log: ' + realign._tail(log, 3)) if realign._tail(log, 3) else ''), rec)
        got = dict(msa_graph.read_msa(out))
        rows = []
        for i, s in enumerate(inp):
            row = got.get('r%d' % i)
            if row is None:
                raise AlignerFailed('error', 'abPOA output has no row r%d' % i, rec)
            if row.replace('-', '') != s:
                raise AlignerFailed('error', 'abPOA row r%d does not spell its input' % i, rec)
            if empty[i]:
                row = '-' * len(row)
            rows.append(row)
        if len(set(len(x) for x in rows)) != 1:
            raise AlignerFailed('error', 'abPOA rows differ in length', rec)
        rec['columns'] = len(rows[0])
        rec['nodes'] = graph_nodes(rows)
        return rows
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ---------------------------------------------------------------- BAR's trim (poaBarAligner.c, ported)
#
# An MSA here is a list of bytearrays (b'-' = gap); column scores are ints (Cactus uses floats, but
# every value is a whole number far below 2^24, so the arithmetic is identical).

GAP = ord('-')


def make_column_scores(msa):
    """max(number of bases in the column - 1, 0) for every column."""
    if not msa:
        return []
    ncol = len(msa[0])
    cs = [0] * ncol
    for row in msa:
        for i in range(ncol):
            if row[i] != GAP:
                cs[i] += 1
    return [c - 1 if c >= 1 else 0 for c in cs]


def sum_column_scores(row, msa, cs):
    """Cumulative column scores over the columns in which `row` has a base (one entry per base)."""
    out, cu = [], 0
    r = msa[row]
    for i in range(len(r)):
        if r[i] != GAP:
            cu += cs[i]
            out.append(cu)
    return out


def trim_msa_suffix(msa, cs, row, suffix_start):
    """Remove the bases of `row` from its suffix_start-th base on, updating the column scores."""
    k = 0
    r = msa[row]
    for i in range(len(r)):
        if r[i] != GAP:
            if k >= suffix_start:
                r[i] = GAP
                cs[i] = cs[i] - 1 if cs[i] > 1 else 0
            k += 1


def trim(row1, msa1, cs1, row2, msa2, cs2, overlap):
    """Make two MSAs consistent for one sequence that is a prefix of row1 (forward) and of row2
    (reverse complement) with `overlap` bases in both; keeps the cut with the largest summed column
    score, ties resolved as poaBarAligner.c does. Returns the cut point (bases of the overlap kept in
    msa1)."""
    if overlap == 0:
        return None
    assert overlap > 0
    seq_len1 = sum(1 for c in msa1[row1] if c != GAP)
    seq_len2 = sum(1 for c in msa2[row2] if c != GAP)
    assert overlap <= seq_len1 and overlap <= seq_len2
    cu1 = sum_column_scores(row1, msa1, cs1)
    cu2 = sum_column_scores(row2, msa2, cs2)
    max_cut_score = cu2[seq_len2 - 1]
    if overlap < seq_len1:
        max_cut_score += cu1[seq_len1 - overlap - 1]
    cut = 0
    for i in range(overlap - 1):
        s = cu1[seq_len1 - overlap + i] + cu2[seq_len2 - i - 2]
        if s > max_cut_score:
            cut = i + 1
            max_cut_score = s
    f = cu1[seq_len1 - 1]
    if overlap < seq_len2:
        f += cu2[seq_len2 - overlap - 1]
    if f > max_cut_score:
        max_cut_score = f
        cut = overlap
    trim_msa_suffix(msa1, cs1, row1, seq_len1 - overlap + cut)
    trim_msa_suffix(msa2, cs2, row2, seq_len2 - cut)
    return cut


# ---------------------------------------------------------------- one region as one BAR flower

def bar_align(seqs, flanks=(0, 0), workdir=None, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB,
              poa=None, region_timeout=None):
    """Align anchor-to-anchor sequences as BAR aligns one flower whose two ends are the anchors.
    flanks = (bases of the left anchor, bases of the right anchor), shared by every sequence.
    Returns (rows, info): rows = aligned strings in input order, each spelling its sequence.
    `poa` replaces poa_msa (for tests). `region_timeout` caps the summed wall clock of the abPOA
    calls (each call gets min(timeout, what is left); status timeout when nothing is left)."""
    poa_one = poa or poa_msa
    t_start = time.time()

    def poa(strings, tag, calls=None, timeout=timeout, **kw):
        if region_timeout:
            left = region_timeout - (time.time() - t_start)
            if left <= 0:
                raise AlignerFailed('timeout', 'region cap of %.0f s used up before the %s call'
                                    % (region_timeout, tag))
            timeout = min(timeout, left) if timeout else left
        return poa_one(strings, tag=tag, calls=calls, timeout=timeout, **kw)
    L, R = flanks
    n = len(seqs)
    first = seqs[0]
    for s in seqs:
        if len(s) < L + R or s[:L] != first[:L] or (R and s[len(s) - R:] != first[len(first) - R:]):
            raise ValueError('sequences do not share %d + %d bp of anchor flank' % (L, R))
    inner = [s[L:len(s) - R] for s in seqs]
    lens = [len(x) for x in inner]
    order = sorted(range(n), key=lambda i: (-lens[i], i))      # longest first; ties in input order
    calls = []
    info = collections.OrderedDict([('flanks', [L, R]), ('banding_limit', BANDING_LIMIT),
                                    ('max_interior', max(lens) if lens else 0), ('region_timeout', region_timeout)])
    wk = dict(workdir=workdir, timeout=timeout, mem_mb=mem_mb, calls=calls)
    if not lens or max(lens) < BANDING_LIMIT:
        info['bar_rule'] = 'single'
        got = poa([inner[i] for i in order], tag='all', **wk)
        rows_in = [None] * n
        for k, i in enumerate(order):
            rows_in[i] = got[k]
        info.update(rows_ge_10kb=0, rows_gt_20kb=0, unaligned_bp=0, unaligned_rows=0)
    else:
        info['bar_rule'] = 'ends'
        pre = [min(x, BANDING_LIMIT) for x in lens]
        ov = [max(0, 2 * p - x) for p, x in zip(pre, lens)]
        left = poa([inner[i][:pre[i]] for i in order], tag='left', **wk)
        right = poa([revcomp(inner[i])[:pre[i]] for i in order], tag='right', **wk)
        mL = [bytearray(r, 'ascii') for r in left]
        mR = [bytearray(r, 'ascii') for r in right]
        csL = make_column_scores(mL)
        csR = make_column_scores(mR)
        cuts = {}
        for k, i in enumerate(order):                           # end 0 (left) row order, as Cactus
            if ov[i] > 0:
                cuts[i] = trim(k, mL, csL, k, mR, csR, ov[i])
        nL = len(mL[0]) if mL else 0
        nR = len(mR[0]) if mR else 0
        keptL, keptR, mid = [0] * n, [0] * n, [''] * n
        rowL, rowR = [None] * n, [None] * n
        for k, i in enumerate(order):
            a = mL[k].decode()
            b = revcomp(mR[k].decode())                         # back to the forward strand
            keptL[i] = len(a) - a.count('-')
            keptR[i] = len(b) - b.count('-')
            mid[i] = inner[i][keptL[i]:lens[i] - keptR[i]]
            if keptL[i] + len(mid[i]) + keptR[i] != lens[i] or a.replace('-', '') != inner[i][:keptL[i]] \
                    or b.replace('-', '') != inner[i][lens[i] - keptR[i]:]:
                raise ValueError('internal: BAR pieces of row %d do not spell it' % i)
            rowL[i], rowR[i] = a, b
        # left MSA columns, then each row's unaligned middle in columns of its own, then the right MSA
        tot_mid = sum(len(x) for x in mid)
        rows_in = []
        off = 0
        for i in range(n):
            m = '-' * off + mid[i] + '-' * (tot_mid - off - len(mid[i]))
            off += len(mid[i])
            rows_in.append((rowL[i] if nL else '') + m + (rowR[i] if nR else ''))
        info.update(rows_ge_10kb=sum(1 for x in lens if x >= BANDING_LIMIT),
                    rows_gt_20kb=sum(1 for x in lens if x > 2 * BANDING_LIMIT),
                    unaligned_bp=tot_mid, unaligned_rows=sum(1 for x in mid if x),
                    overlap_rows=len(cuts), cut_at_left_end=sum(1 for c in cuts.values() if c == 0),
                    cut_at_right_end=sum(1 for i, c in cuts.items() if c == ov[i]))
    rows = [first[:L] + r + first[len(first) - R:] for r in rows_in] if R else [first[:L] + r for r in rows_in]
    for i in range(n):
        if rows[i].replace('-', '') != seqs[i]:
            raise ValueError('internal: row %d does not spell its sequence' % i)
    if len(set(len(r) for r in rows)) > 1:
        raise ValueError('internal: rows differ in length')
    info['abpoa_calls'] = len(calls)
    info['calls'] = calls
    info['progressive'] = ','.join(str(int(c['progressive'])) for c in calls)
    return rows, info


# ---------------------------------------------------------------- realign.py plugin-style entry points

def mcpoa_align(in_fa, out_fa, threads=1, workdir=None, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB,
                flanks=(0, 0), region_timeout=None, **_):
    """realign.Method-compatible aligner: in_fa holds the masked sequences s0..sK (anchor to anchor);
    writes the MSA to out_fa and returns the status dict."""
    recs = msa_graph.read_fasta(in_fa)
    names = [n for n, _ in recs]
    seqs = [s.upper() for _, s in recs]
    t0 = time.time()
    res = collections.OrderedDict([('tool_version', tool_version()), ('order', 'longest (then -p guide tree)'),
                                   ('command', None)])
    try:
        rows, info = bar_align(seqs, tuple(flanks), workdir=workdir, timeout=timeout, mem_mb=mem_mb,
                               region_timeout=region_timeout)
    except AlignerFailed as e:
        res.update(status=e.status, message=e.message, seconds=round(time.time() - t0, 2))
        if e.record:
            res['failed_call'] = e.record
            res['peak_rss_mb'] = e.record.get('peak_rss_mb') or e.record.get('sampled_rss_mb')
        return res
    res.update(info)
    res['command'] = info['calls'][0]['command'] if info['calls'] else None
    res['status'] = 'ok'
    res['seconds'] = round(time.time() - t0, 2)
    res['abpoa_s'] = round(sum(c['seconds'] for c in info['calls']), 2)
    peaks = [c['peak_rss_mb'] or 0 for c in info['calls']]
    res['peak_rss_mb'] = max(peaks) if peaks else 0.0
    res['sampled_rss_mb'] = max([c['sampled_rss_mb'] or 0 for c in info['calls']] or [0.0])
    res['predicted_mb'] = max([c['predicted_mb'] for c in info['calls']] or [0.0])
    msa_graph.write_msa(list(zip(names, rows)), out_fa)
    return res


def make_method(flanks, region_timeout=None):
    import realign
    params = {'flanks': list(flanks)}
    if region_timeout:
        params['region_timeout'] = region_timeout
    return realign.Method(METHOD, mcpoa_align, params, DESCRIPTION, 'abpoa')


def align_fasta(input_fa, msa_out, flanks=(0, 0), timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB, workdir=None,
                region_timeout=None):
    """MSA-only mode through realign's shared pipeline (N masking, row checks), dedup off.
    region_timeout: cap on the summed abPOA wall clock (both end calls), see bar_align."""
    import realign
    wd = workdir or os.path.join(config.WORK_DIR, 'mcpoa', 'tmp')
    os.makedirs(wd, exist_ok=True)
    return realign.align_fasta(make_method(flanks, region_timeout), input_fa, msa_out, threads=1, timeout=timeout,
                               mem_mb=mem_mb, workdir=wd, dedup=False)


# ---------------------------------------------------------------- regions

def region_flanks(rd):
    """(len(anchor_left_seq), len(anchor_right_seq)) from region.json."""
    j = json.load(open(os.path.join(rd, 'region.json')))
    return len(j['anchor_left_seq']), len(j['anchor_right_seq'])


def plan_region(seqs, flanks):
    """What BAR will align in this region, without running it: the abPOA calls' string lengths
    (after N masking, as the pipeline does) and the predicted peak."""
    import realign
    L, R = flanks
    lens = [len(realign.mask_runs(s.upper())[0]) - L - R for s in seqs]
    if max(lens) < BANDING_LIMIT:
        call_lens = [lens]
    else:
        pre = [min(x, BANDING_LIMIT) for x in lens]
        call_lens = [pre, pre]
    return {'max_interior': max(lens), 'predicted_mb': max(predict_call_mb(c) for c in call_lens),
            'bar_rule': 'single' if len(call_lens) == 1 else 'ends'}


def memory_free_mb():
    """Free memory per memory_pressure (its 'free percentage' x physical memory); None if unknown."""
    try:
        out = subprocess.run(['memory_pressure'], capture_output=True, text=True, timeout=60).stdout
        tot = int(subprocess.run(['sysctl', '-n', 'hw.memsize'], capture_output=True, text=True,
                                 timeout=20).stdout.strip())
        m = re.search(r'free percentage:\s*(\d+)%', out)
        return tot / 2.0 ** 20 * int(m.group(1)) / 100.0 if m else None
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


def wait_for_memory(need_mb, max_wait=3600, poll=30, log=None):
    t0 = time.time()
    while True:
        free = memory_free_mb()
        if free is None or free >= need_mb:
            return round(time.time() - t0, 1), free
        if time.time() - t0 > max_wait:
            return round(time.time() - t0, 1), free
        if log:
            log('waiting for %.0f MB free (memory_pressure: %.0f MB)' % (need_mb, free))
        time.sleep(poll)



def update_runtime(row, path=RUNTIME_TSV):
    """Replace the (method, region_id) row of the runtime table, under a file lock."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + '.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        rows = []
        if os.path.exists(path):
            with open(path) as f:
                rows = list(csv.DictReader(f, delimiter='\t'))
        rows = [r for r in rows if (r.get('method'), r.get('region_id')) != (row['method'], row['region_id'])]
        rows.append(row)
        rows.sort(key=lambda r: (r['method'], r['region_id']))
        tmp = path + '.tmp%d' % os.getpid()
        with open(tmp, 'w', newline='') as f:
            w = csv.DictWriter(f, RUNTIME_COLS, delimiter='\t', extrasaction='ignore', lineterminator='\n')
            w.writeheader()
            for r in rows:
                w.writerow(dict((c, r.get(c, '')) for c in RUNTIME_COLS))
        os.replace(tmp, path)
        fcntl.flock(lk, fcntl.LOCK_UN)


def run_region(rd, cand_root=None, workroot=None, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB, force=False,
               retry_failed=False, runtime_tsv=RUNTIME_TSV, wait_memory=True, merge_blocks=False, log=None):
    """Align regions/<id>/hap32.fa with poa_abpoa_mc, build the graph, write
    <cand_root>/poa_abpoa_mc/<id>.{msa.fa,gfa,realign.json}; returns the runtime row."""
    import realign
    rid = os.path.basename(os.path.normpath(rd))
    o = realign.outputs(METHOD, rid, cand_root)
    prev = realign.previous_status(METHOD, rid, cand_root)
    if prev and not force and not (retry_failed and prev != 'ok'):
        return {'method': METHOD, 'region_id': rid, 'status': 'skipped', 'note': 'kept %s' % prev}
    os.makedirs(o['dir'], exist_ok=True)
    workroot = workroot or os.path.join(config.WORK_DIR, 'mcpoa', 'tmp')
    os.makedirs(workroot, exist_ok=True)
    stratum, span = realign.region_meta(rd)
    hap = os.path.join(rd, 'hap32.fa')
    recs = msa_graph.read_fasta(hap)
    flanks = region_flanks(rd)
    plan = plan_region([s for _, s in recs], flanks)
    t0 = time.time()
    for k in ('msa', 'gfa', 'frag'):
        if os.path.exists(o[k]):
            os.remove(o[k])
    info = None
    waited, free = 0.0, None
    if mem_mb and plan['predicted_mb'] > mem_mb:
        info = collections.OrderedDict([('method', METHOD), ('status', 'memout'), ('message',
                'not run: predicted %.0f MB > cap %.0f MB' % (plan['predicted_mb'], mem_mb))])
    else:
        if wait_memory and plan['predicted_mb'] > WAIT_ABOVE_MB:
            waited, free = wait_for_memory(plan['predicted_mb'] + 2000, log=log)
            if free is not None and free < plan['predicted_mb'] + 2000:
                info = collections.OrderedDict([('method', METHOD), ('status', 'skipped_memory'), ('message',
                        'not run: memory_pressure reports %.0f MB free after %.0f s, need %.0f MB' % (
                            free, waited, plan['predicted_mb'] + 2000))])
        if info is None:
            tmp_msa = o['msa'] + '.tmp%d' % os.getpid()
            info = realign.align_fasta(make_method(flanks), hap, tmp_msa, threads=1, timeout=timeout, mem_mb=mem_mb,
                                       workdir=workroot, dedup=False)
    info['region_id'] = rid
    info['stratum'] = stratum
    info['span_bp'] = span
    info['description'] = DESCRIPTION
    info['flanks'] = list(flanks)
    info['plan'] = plan
    info['memory_wait_s'] = waited
    info['memory_free_mb_at_start'] = free
    info['n_distinct_seqs'] = len(set(s.upper() for _, s in recs))
    info['msa'] = realign._rel(o['msa']) if info['status'] == 'ok' else None
    graph_s = None
    if info['status'] == 'ok':
        tg = time.time()
        tmp_gfa = o['gfa'] + '.tmp%d' % os.getpid()
        try:
            st = msa_graph.msa_to_gfa(tmp_msa, hap, tmp_gfa, merge_blocks=merge_blocks, workdir=workroot)
            os.replace(tmp_msa, o['msa'])
            os.replace(tmp_gfa, o['gfa'])
            st['msa'] = realign._rel(o['msa'])
            st['gfa'] = realign._rel(o['gfa'])
            st['hap32'] = realign._rel(hap)
            info['graph'] = st
            info['gfa'] = st['gfa']
        except msa_graph.MsaGraphError as e:
            info['status'] = 'graph_error'
            info['message'] = str(e)
        graph_s = round(time.time() - tg, 2)
        for p in (tmp_msa, tmp_gfa):
            if os.path.exists(p):
                os.remove(p)
    info['graph_s'] = graph_s
    info['total_s'] = round(time.time() - t0, 2)
    info['finished'] = datetime.datetime.now().isoformat(timespec='seconds')
    info['host_load'] = round(os.getloadavg()[0], 1)
    with open(o['json'] + '.tmp', 'w') as f:
        json.dump(info, f, indent=1)
    os.replace(o['json'] + '.tmp', o['json'])
    g = info.get('graph') or {}
    al = info.get('aligner') or {}
    row = {'method': METHOD, 'region_id': rid, 'stratum': stratum, 'span_bp': span,
           'n_seqs': info.get('n_seqs', len(recs)), 'n_distinct': info['n_distinct_seqs'],
           'max_len': info.get('max_len', max(len(s) for _, s in recs)),
           'max_interior': al.get('max_interior', plan['max_interior']), 'total_bp': info.get('total_bp'),
           'status': info['status'], 'bar_rule': al.get('bar_rule', plan['bar_rule']),
           'rows_ge_10kb': al.get('rows_ge_10kb', ''), 'rows_gt_20kb': al.get('rows_gt_20kb', ''),
           'unaligned_bp': al.get('unaligned_bp', ''), 'abpoa_calls': al.get('abpoa_calls', ''),
           'progressive': al.get('progressive', ''), 'align_s': info.get('align_s'), 'abpoa_s': al.get('abpoa_s', ''),
           'graph_s': graph_s, 'total_s': info['total_s'], 'peak_rss_mb': al.get('peak_rss_mb', ''),
           'sampled_rss_mb': al.get('sampled_rss_mb', ''), 'predicted_mb': plan['predicted_mb'],
           'waited_s': waited, 'timeout_s': timeout, 'mem_cap_mb': mem_mb, 'columns': info.get('columns', ''),
           'gap_frac': info.get('gap_frac', ''), 'nodes': g.get('nodes', ''), 'nodes_per_kb': g.get('nodes_per_kb', ''),
           'frac_nodes_1bp': g.get('frac_nodes_1bp', ''), 'masked_runs': info.get('masked_runs', ''),
           'masked_bp': info.get('masked_bp', ''), 'finished': info['finished'],
           'note': (info.get('message') or '')[:200].replace('\t', ' ').replace('\n', ' ')}
    if runtime_tsv:
        update_runtime(row, runtime_tsv)
    return row


# ---------------------------------------------------------------- CLI

def main(argv=None):
    import realign
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('regions', nargs='*', help="region dirs or ids; 'all' = every packaged region")
    ap.add_argument('--method', default=METHOD, choices=[METHOD], help='the one method (for CLI parity)')
    ap.add_argument('--out', help='candidates root (default %s)' % config.CANDIDATES_DIR)
    ap.add_argument('--input', help='MSA-only mode: align this FASTA ...')
    ap.add_argument('--msa-out', help='... and write only the MSA here')
    ap.add_argument('--json', help='MSA-only mode: also write the info dict here')
    ap.add_argument('--region', help='MSA-only mode: take the anchor flanks from this region (dir or id)')
    ap.add_argument('--flanks', help='MSA-only mode: anchor flank lengths L,R (default 0,0)')
    ap.add_argument('--project-to', help='MSA-only mode: also project the MSA onto this hap32.fa ...')
    ap.add_argument('--projected-out', help='... and write the projected MSA here')
    ap.add_argument('--timeout', type=float, default=DEFAULT_TIMEOUT, help='seconds per abPOA call (default 1800)')
    ap.add_argument('--mem-mb', type=float, default=DEFAULT_MEM_MB,
                    help='predicted or sampled RSS above this stops the region (default %d)' % DEFAULT_MEM_MB)
    ap.add_argument('--jobs', type=int, default=1, help='regions at once (default 1; at most 2)')
    ap.add_argument('--stratum', help='comma list of strata (with all)')
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--retry-failed', action='store_true')
    ap.add_argument('--runtime-tsv', default=RUNTIME_TSV, help='runtime table (default results/mcpoa_runtime.tsv)')
    ap.add_argument('--no-runtime', action='store_true')
    ap.add_argument('--merge-blocks', action='store_true')
    ap.add_argument('--workdir', help='scratch for abPOA runs (default $VNTR_WORK/mcpoa/tmp)')
    ap.add_argument('--no-wait-memory', dest='wait_memory', action='store_false')
    a = ap.parse_args(argv)
    wd = a.workdir or os.path.join(config.WORK_DIR, 'mcpoa', 'tmp')
    os.makedirs(wd, exist_ok=True)

    if a.input or a.msa_out:
        if not (a.input and a.msa_out) or a.regions:
            ap.error('MSA-only mode takes --input FASTA and --msa-out FILE and no regions')
        if bool(a.project_to) != bool(a.projected_out):
            ap.error('--project-to and --projected-out go together')
        if a.region and a.flanks:
            ap.error('--region or --flanks, not both')
        flanks = (region_flanks(realign.region_dir(a.region)) if a.region else
                  tuple(int(x) for x in a.flanks.split(',')) if a.flanks else (0, 0))
        info = align_fasta(a.input, a.msa_out, flanks, timeout=a.timeout, mem_mb=a.mem_mb, workdir=wd)
        info['flanks'] = list(flanks)
        if a.json:
            with open(a.json, 'w') as f:
                json.dump(info, f, indent=1)
        al = info.get('aligner') or {}
        sys.stderr.write('realign_poa_mc: %s: %s, %s rows, %s columns, rule %s, unaligned %s bp, %.1f s%s\n' % (
            a.input, info['status'], info.get('n_seqs'), info.get('columns'), al.get('bar_rule'),
            al.get('unaligned_bp'), info.get('align_s', 0), (' -- ' + info['message']) if info.get('message') else ''))
        if info['status'] == 'ok' and a.project_to:
            import realign_poa
            try:
                pr = realign_poa.project_to_hap32(a.msa_out, a.project_to, a.projected_out)
            except ValueError as e:
                sys.stderr.write('realign_poa_mc: projection failed: %s\n' % e)
                return 2
            sys.stderr.write('realign_poa_mc: projected %d hap32 rows, %d -> %d columns: %s\n' % (
                pr['rows'], pr['columns_in'], pr['columns_out'], a.projected_out))
        return 0 if info['status'] == 'ok' else (4 if info['status'] in ('timeout', 'memout') else 2)

    if not a.regions:
        ap.error("give region dirs or ids, or 'all'")
    if a.regions == ['all']:
        rows = realign.region_list()
        if a.stratum:
            keep = set(a.stratum.split(','))
            rows = [r for r in rows if r.get('stratum') in keep]
        rdirs = [os.path.join(config.REGIONS_DIR, r['region_id']) for r in rows]
    else:
        rdirs = [realign.region_dir(x) for x in a.regions]

    def say(msg):
        sys.stderr.write('realign_poa_mc: %s\n' % msg)
        sys.stderr.flush()

    def job(rd):
        try:
            return run_region(rd, a.out, wd, a.timeout, a.mem_mb, a.force, a.retry_failed,
                              None if a.no_runtime else a.runtime_tsv, a.wait_memory, a.merge_blocks,
                              log=lambda m: say('%s %s' % (os.path.basename(rd), m)))
        except Exception as e:                  # keep the batch going
            row = {'method': METHOD, 'region_id': os.path.basename(rd), 'status': 'crash',
                   'note': '%s: %s' % (type(e).__name__, e)}
            if not a.no_runtime:
                update_runtime(row, a.runtime_tsv)
            return row

    t0 = time.time()
    counts = collections.Counter()
    worst = 0
    ex = concurrent.futures.ThreadPoolExecutor(max_workers=max(1, a.jobs)) if a.jobs > 1 else None
    it = ex.map(job, rdirs) if ex else map(job, rdirs)
    for row in it:
        counts[row['status']] += 1
        if row['status'] not in ('ok', 'skipped'):
            worst = 2
        say('%-9s %-8s %-6s align %6ss abpoa %6ss graph %6ss peak %6s MB nodes %-6s unaligned %-6s %s' % (
            row['region_id'], row['status'], row.get('bar_rule', ''), row.get('align_s', ''), row.get('abpoa_s', ''),
            row.get('graph_s', ''), row.get('peak_rss_mb', ''), row.get('nodes', ''), row.get('unaligned_bp', ''),
            row.get('note', '')))
    if ex:
        ex.shutdown()
    say('done in %.0f s: %s' % (time.time() - t0, ', '.join('%s %d' % kv for kv in sorted(counts.items()))))
    return worst


if __name__ == '__main__':
    sys.exit(main())
