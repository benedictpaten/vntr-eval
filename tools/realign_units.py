#!/usr/bin/env python3
"""realign_units.py -- M4: repeat-unit-aware two-level realignment of a region's haplotypes.

Region mode (default): align regions/<id>/hap32.fa and build the candidate graph; writes

    candidates/unit_aware/<id>.msa.fa        aligned rows ('-' gaps), hap32.fa names and order
    candidates/unit_aware/<id>.gfa           GFA 1.0 from tools/msa_graph.py (one P line per row)
    candidates/unit_aware/<id>.realign.json  how it was made: mode (unit / fallback) and why, the
                                             repeat model, symbols, timings, graph stats
    results/realign_runtime.tsv              one row per (method, region) (shared with realign.py)

    python3 tools/realign_units.py L012184                 # one region (id or regions/<id>)
    python3 tools/realign_units.py all --jobs 2 --threads 2
    python3 tools/realign_units.py all --stratum hotspot_vntr --force

MSA-only mode, for an ARBITRARY input FASTA (e.g. the deduplicated full HPRC panel of a region);
the motif / period and the reference array come from --region-json (and the CHM13 record of the
input, or of the region's hap32.fa when the input has none); writes only the MSA:

    python3 tools/realign_units.py --input seqs.fa[.gz] --region-json regions/<id>/region.json \\
                                   --msa-out out.msa.fa [--info-out info.json] [--threads 2] \\
                                   [--fallback mafft_linsi]

    import sys; sys.path.insert(0, 'tools')
    from realign_units import align_units_fasta
    info = align_units_fasta('seqs.fa', 'regions/L012184/region.json', 'out.msa.fa', threads=2)
    info['status'], info['mode']      # 'ok' / 'error';  'unit' or 'fallback' (and info['reason'])

The algorithm (details, parameters and measurements in tools/UNIT_ALIGN.md):
  1. Repeat model from CHM13 + region.json: the period is chosen among the annotated period, the
     motif length, the CHM13 tandem-scan periods and their divisors by lag identity on the CHM13
     core; the unit is refined to the consensus of a wrap-around alignment of the CHM13 core to it,
     and rotated so that the CHM13 array starts at unit position 0. No period, a period < 2, fewer
     than 3 units in CHM13, unit identity < 0.6 or an array covering < 50% of the core: fall back.
  2. Every (distinct, N-masked) sequence is aligned to [left flank] -> (unit cycle)* -> [right flank]
     (the CHM13 sequence next to the array, up to --flank-k bp each side; gap-affine DP in C,
     tools/unit_dp.c). This finds the array in every sequence and cuts it into units at the
     cycle's wrap, so all sequences are phased by the same unit; partial units at the ends and
     units with indels are allowed, and every base knows its unit position.
  3. Distinct unit strings are the symbols (identical units share a symbol; near-identical ones
     stay distinct; their similarity is the edit distance between the unit strings).
  4. The unit strings are aligned by an exact sum-of-pairs progressive aligner (C): substituting
     unit a for b costs ed(a, b), gapping a unit costs its length -- so the unit-level objective is
     the base-level edit cost -- with a UPGMA guide tree, then refinement (every guide-tree split
     and every row realigned against the rest, kept when the cost falls). Whole units are gapped
     and units stay in phase. mafft --textmatrix on the symbols is kept as an option
     (--unit-aligner mafft; --anysymbol cannot be used: it scores every non-amino-acid character
     as unknown); it was measured and is worse (UNIT_ALIGN.md).
  5. Each unit column is expanded to bases through the unit (every unit variant is already
     aligned to it by the DP; insertions between the same two unit positions are star-aligned).
     The left and right flank pieces are aligned with mafft (realign.py's pipeline, L-INS-i).
  6. Base-level polish: every row is realigned against the others within a band around its
     path (same exact SP cost), kept when its pairs get cheaper; this removes the unit-boundary
     shifts a single SNP can cause in a degenerate array.
  7. Adequacy guard: when the MSA costs more than --guard-ratio (1.2) times the sum of pairwise
     optimal edit distances, the fallback is computed too and the MSA with the lower cost is kept.
  8. N runs are cut out before step 2 and put back as their own columns (as realign.py does);
     the rows are checked to spell their input; the graph is built by msa_graph.msa_to_gfa.

Fallback (non-TR strata, unusable motif, the guard, or any failure of the unit path): realign.py's
mafft_linsi pipeline on the whole sequences (--fallback picks another realign.py method).
--fallback-engine abpoa uses abPOA (poa_abpoa) for the fallback, the guard's fallback and the flank
pieces instead of mafft.
"""
import argparse
import collections
import concurrent.futures
import ctypes
import datetime
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from array import array

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402
import msa_graph  # noqa: E402

METHOD = 'unit_aware'
DEFAULT_FALLBACK = 'mafft_linsi'
VERSION = 1

# The aligner for everything that is not a unit array: the flank pieces, the fallback (no usable
# motif, or the unit path failed) and the adequacy guard's fallback. Each engine sets the three
# realign.py methods; an explicit flank_method / fallback / fallback_large overrides its choice.
FALLBACK_ENGINES = collections.OrderedDict([
    ('mafft', {'flank_method': 'auto', 'fallback': DEFAULT_FALLBACK, 'fallback_large': 'mafft_fftnsi'}),
    ('abpoa', {'flank_method': 'poa_abpoa', 'fallback': 'poa_abpoa', 'fallback_large': 'poa_abpoa'}),
])

# mafft --text: every byte but newline, CR, space, '-', '<', '=', '>' (checked with maffttext2hex)
_BAD = set([0x0a, 0x0d, 0x20, 0x2d, 0x3c, 0x3d, 0x3e])
SYMBOLS = bytes(b for b in range(1, 256) if b not in _BAD)       # 248
MAX_SYMBOLS = len(SYMBOLS)


class Params(object):
    """Every tunable of the method (defaults are the shipped configuration)."""

    def __init__(self, **kw):
        self.ma, self.mi, self.go, self.ge, self.nsc = 2, 4, 6, 1, -1   # DP scores (match, mismatch, open, extend)
        self.flank_k = 300          # CHM13 flank bases each side of the array used to anchor it
        self.min_period = 2
        self.max_period = 1000
        self.min_units = 3          # CHM13 units (>= half a unit long) needed
        self.min_identity = 0.60    # CHM13 array identity to the consensus unit
        self.min_core_cov = 0.50    # share of the annotated core covered by the CHM13 array
        self.period_tol = 0.02      # smallest period whose lag identity is within this of the best
        self.unit_aligner = 'prog'  # prog: exact sum-of-pairs progressive + refinement; mafft: --textmatrix
        self.refine_rounds = 3      # prog: leave-one-out refinement rounds (stops early without gain)
        self.unit_go = 0.0          # prog: gap-open cost per pair of rows and unit-level gap run (base-edit units)
        self.polish_rounds = 2      # base-level leave-one-out polish of the final MSA (0 = off)
        self.polish_go = 0.5        # gap-open cost per pair in the polish (keeps gaps contiguous on ties)
        self.polish_band = 100      # band half-width: max(polish_band, 3 x period) bases around the row's path
        self.polish_max_cells = 30000000   # band cells per row realignment (24 bytes each)
        self.polish_max_s = 300
        self.tree_exact_max = 150   # prog: exact unit-level pair costs for the guide tree up to this many strings
        self.refine = 'tree'        # prog: tree (every guide-tree split, then leave-one-out) | loo | none
        self.refine_tree_max = 120  # prog: tree splits only up to this many strings (leave-one-out beyond)
        self.refine_max_s = 300     # prog: stop refining after this many seconds
        self.unit_mafft = 'auto'    # linsi | einsi | ginsi | nwnsi | nwns2 | auto
        self.unit_op = None         # mafft --op for the unit alignment (None = mafft default)
        self.unit_ep = None         # mafft --ep
        self.matrix = 'norm'        # symbol score: norm = 10*(1 - 2*ed/max len); sp = la+lb-ed
        self.max_symbols = MAX_SYMBOLS
        self.flank_method = 'auto'  # realign.py method for the flank pieces (auto: linsi, fftns2 if long)
        self.fallback = DEFAULT_FALLBACK
        self.fallback_large = 'mafft_fftnsi'   # fallback when there are more than large_n distinct sequences
        self.large_n = 100
        self.guard_ratio = 1.2      # MSA SP / sum of pairwise optima above this: also run the fallback, keep the lower SP
        self.guard_max_pairs = 3000 # pairs used for that ratio (seeded sample beyond)
        self.guard_max_cells = 3e10 # bit-parallel edit-distance budget for the sample (sum of len a x len b)
        self.timeout = 900
        self.mem_mb = 10000
        self.fallback_engine = kw.pop('fallback_engine', None) or 'mafft'
        if self.fallback_engine not in FALLBACK_ENGINES:
            raise ValueError('unknown fallback engine %s (%s)' % (self.fallback_engine, ', '.join(FALLBACK_ENGINES)))
        for k, v in FALLBACK_ENGINES[self.fallback_engine].items():
            setattr(self, k, v)
        for k, v in kw.items():
            if not hasattr(self, k):
                raise TypeError('unknown parameter %s' % k)
            setattr(self, k, v)

    def as_dict(self):
        return dict((k, v) for k, v in self.__dict__.items())


# ---------------------------------------------------------------- the C dynamic program

_LIB = None


def lib():
    """tools/bin/libunitdp.so, (re)built from tools/unit_dp.c when missing or stale."""
    global _LIB
    if _LIB is not None:
        return _LIB
    src = os.path.join(TOOLS, 'unit_dp.c')
    so = os.path.join(TOOLS, 'bin', 'libunitdp.so')
    if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
        os.makedirs(os.path.dirname(so), exist_ok=True)
        tmp = so + '.tmp%d' % os.getpid()
        subprocess.run([config.CC, '-O3', '-shared', '-fPIC', '-o', tmp, src], check=True, env=config.tool_env())
        os.replace(tmp, so)
    L = ctypes.CDLL(so)
    i32, cp = ctypes.c_int32, ctypes.c_char_p
    L.unit_align.restype = i32
    L.unit_align.argtypes = [cp, i32, cp, i32, cp, i32, cp, i32, i32, i32, i32, i32, i32, i32,
                             ctypes.POINTER(ctypes.c_uint8), ctypes.POINTER(i32), i32, ctypes.POINTER(i32)]
    L.edit_distance.restype = i32
    L.edit_distance.argtypes = [cp, i32, cp, i32]
    L.edit_matrix.restype = None
    L.edit_matrix.argtypes = [cp, ctypes.POINTER(i32), i32, ctypes.POINTER(i32)]
    L.edit_cross.restype = None
    L.edit_cross.argtypes = [cp, ctypes.POINTER(i32), i32, cp, ctypes.POINTER(i32), i32, ctypes.POINTER(i32)]
    pi, pd = ctypes.POINTER(i32), ctypes.POINTER(ctypes.c_double)
    L.profile_align.restype = ctypes.c_double
    L.profile_align.argtypes = [i32, pi, pi, pd, ctypes.c_double, i32, pi, pi, pd, ctypes.c_double, pd, i32,
                                ctypes.POINTER(ctypes.c_uint8), pi]
    L.profile_align_affine.restype = ctypes.c_double
    L.profile_align_affine.argtypes = [i32, pi, pi, pd, ctypes.c_double, i32, pi, pi, pd, ctypes.c_double, pd, i32,
                                       ctypes.c_double, ctypes.POINTER(ctypes.c_uint8), pi]
    L.sp_affine.restype = ctypes.c_double
    L.sp_affine.argtypes = [pi, i32, i32, pd, pd, i32, ctypes.c_double]
    L.unit_pair_costs.restype = None
    L.unit_pair_costs.argtypes = [pi, pi, i32, pd, i32, pd]
    L.kmer_dists.restype = None
    L.kmer_dists.argtypes = [pi, pi, i32, i32, pd]
    L.profile_align_flat.restype = i32
    L.profile_align_flat.argtypes = [pi, i32, i32, pd, pi, i32, i32, pd, pd, i32, ctypes.c_double,
                                     ctypes.POINTER(ctypes.c_uint8), pd]
    L.merge_flat.restype = None
    L.merge_flat.argtypes = [pi, i32, i32, pi, i32, i32, ctypes.POINTER(ctypes.c_uint8), i32, pi]
    L.drop_gap_cols.restype = i32
    L.drop_gap_cols.argtypes = [pi, i32, i32, pi]
    L.row_to_profile_banded.restype = i32
    L.row_to_profile_banded.argtypes = [pi, i32, i32, pd, pi, i32, pi, pi, pd, i32, ctypes.c_double,
                                        ctypes.POINTER(ctypes.c_uint8), pd]
    L.band_extract.restype = i32
    L.band_extract.argtypes = [pi, i32, i32, i32, i32, pi, pi, pi, pi, pi]
    L.ed_bitpar.restype = ctypes.c_int64
    L.ed_bitpar.argtypes = [cp, ctypes.c_int64, cp, ctypes.c_int64]
    L.induced_cost.restype = ctypes.c_int64
    L.induced_cost.argtypes = [cp, cp, ctypes.c_int64]
    L.cross_cost.restype = ctypes.c_double
    L.cross_cost.argtypes = [pi, i32, pi, i32, pi, i32, pd, pd, i32, ctypes.c_double]
    _LIB = L
    return L


def dp_align(seq, fl, motif, fr, P, free_prefix, free_suffix, free_entry=False):
    """Align seq to FL -> motif* -> FR; returns (score, ops bytes, nodes list)."""
    L = lib()
    n = len(seq)
    flags = (1 if free_prefix else 0) | (2 if free_suffix else 0) | (4 if free_entry else 0)
    cap = n + len(fl) + len(fr) + 4 * len(motif) + 1024
    for _ in range(3):
        op = (ctypes.c_uint8 * cap)()
        nd = (ctypes.c_int32 * cap)()
        sc = ctypes.c_int32()
        r = L.unit_align(seq.encode(), n, fl.encode(), len(fl), motif.encode(), len(motif), fr.encode(), len(fr),
                         P.ma, P.mi, P.go, P.ge, P.nsc, flags, op, nd, cap, ctypes.byref(sc))
        if r == -2 ** 31:
            raise MemoryError('unit_align: allocation failed (%d bp x %d states)' % (n, len(fl) + len(motif) + len(fr)))
        if r < 0:
            cap = -r + 16
            continue
        return sc.value, bytes(op[:r]), list(nd[:r])
    raise RuntimeError('unit_align: path buffer')


def _packed(strs):
    buf = b''.join(s.encode() for s in strs)
    off = (ctypes.c_int32 * (len(strs) + 1))()
    k = 0
    for i, s in enumerate(strs):
        off[i] = k
        k += len(s)
    off[len(strs)] = k
    return buf, off


def edit_matrix(strs):
    n = len(strs)
    buf, off = _packed(strs)
    out = (ctypes.c_int32 * (n * n))()
    lib().edit_matrix(buf, off, n, out)
    return [list(out[i * n:(i + 1) * n]) for i in range(n)]


def edit_cross(qs, rs):
    qb, qo = _packed(qs)
    rb, ro = _packed(rs)
    out = (ctypes.c_int32 * (len(qs) * len(rs)))()
    lib().edit_cross(qb, qo, len(qs), rb, ro, len(rs), out)
    nr = len(rs)
    return [list(out[i * nr:(i + 1) * nr]) for i in range(len(qs))]


# ---------------------------------------------------------------- path parsing

class Unit(object):
    __slots__ = ('bases', 'slots', 'events')

    def __init__(self):
        self.bases, self.slots, self.events = [], [], []   # slot: 2j = unit position j, 2j+1 = after j

    @property
    def seq(self):
        return ''.join(self.bases)


class Seg(object):
    """One sequence cut into left piece, units, right piece."""

    def __init__(self, seq, a, b, units, score):
        self.seq, self.a, self.b, self.units, self.score = seq, a, b, units, score

    @property
    def left(self):
        return self.seq[:self.a]

    @property
    def right(self):
        return self.seq[self.b:]


def parse_path(seq, ops, nodes, k1, p, score):
    i, a, b, units, cur = 0, None, None, [], None
    c1 = k1 + p
    for o, v in zip(ops, nodes):
        if o == 0 or o == 4:
            i += 1
        elif o == 1 or o == 2:
            if k1 <= v < c1:
                j = v - k1
                cur.bases.append(seq[i])
                cur.slots.append(2 * j + (1 if o == 2 else 0))
                cur.events.append(('M' if o == 1 else 'I', j, seq[i]))
            i += 1
        elif o == 3:
            if k1 <= v < c1:
                cur.events.append(('D', v - k1, ''))
        elif o == 5 or o == 6:
            if o == 6:
                a = i
            cur = Unit()
            units.append(cur)
        elif o == 7:
            b = i
            cur = None
        elif o == 8:
            a = b = i
    if a is None:
        a = b = i
    if b is None:
        b = i
    units = [u for u in units if u.bases]
    if i != len(seq):
        raise RuntimeError('unit path consumed %d of %d bases' % (i, len(seq)))
    if sum(len(u.bases) for u in units) != b - a:
        raise RuntimeError('unit path: array bases do not add up')
    return Seg(seq, a, b, units, score)


def segment(seq, fl, motif, fr, P, free_prefix, free_suffix, free_entry=False):
    sc, ops, nodes = dp_align(seq, fl, motif, fr, P, free_prefix, free_suffix, free_entry)
    return parse_path(seq, ops, nodes, len(fl), len(motif), sc)


# ---------------------------------------------------------------- the repeat model

def lag_identity(s, p):
    n = len(s) - p
    if n <= 0:
        return 0.0
    return sum(1 for i in range(n) if s[i] == s[i + p]) / float(n)


def _divisors(p):
    return [d for d in range(2, p + 1) if p % d == 0]


def _period_candidates(rj, core_len):
    c = set()
    ta = rj.get('tr_annotation') or {}
    motif = (rj.get('motif') or '').upper()
    for v in (rj.get('period'), ta.get('period'),
              len(motif) if motif and set(motif) <= set('ACGT') and len(motif) < 200 else None,
              ((ta.get('reference_tandem_scan') or {}).get('vntr') or {}).get('period'),
              ((ta.get('reference_tandem_scan') or {}).get('str') or {}).get('period')):
        try:
            v = int(round(float(v)))
        except (TypeError, ValueError):
            continue
        if v >= 2:
            c.update(_divisors(v))
    return sorted(x for x in c if 2 <= x <= core_len // 3)


def _best_window(s, c0, c1, p):
    """Start of the p-window in s[c0:c1] that agrees best with its neighbours at +-p, +-2p."""
    best, bx = -1, c0
    step = max(1, p // 4)
    for x in range(c0, max(c0 + 1, c1 - p + 1), step):
        w = s[x:x + p]
        agree = 0
        for k in (-2, -1, 1, 2):
            y = x + k * p
            if y < c0 or y + p > c1:
                continue
            u = s[y:y + p]
            agree += sum(1 for q in range(p) if w[q] == u[q])
        if agree > best:
            best, bx = agree, x
    return s[bx:bx + p]


def consensus_unit(units, p):
    """Majority unit of a wrap-around alignment: per position the majority base (or a deletion),
    plus an insertion after a position carried by more than half the units covering it."""
    base = [collections.Counter() for _ in range(p)]
    dels = [0] * p
    ins = [collections.Counter() for _ in range(p)]
    cover = [0] * p
    for u in units:
        seen = set()
        cur_ins = collections.defaultdict(list)
        for kind, j, ch in u.events:
            if kind == 'M':
                base[j][ch] += 1
                seen.add(j)
            elif kind == 'D':
                dels[j] += 1
                seen.add(j)
            else:
                cur_ins[j].append(ch)
        for j in seen:
            cover[j] += 1
        for j, chs in cur_ins.items():
            ins[j][''.join(chs)] += 1
    out = []
    for j in range(p):
        nm = sum(base[j].values())
        if nm and nm >= dels[j]:
            out.append(base[j].most_common(1)[0][0])
        ni = sum(ins[j].values())
        if cover[j] and ni * 2 > cover[j]:
            out.append(ins[j].most_common(1)[0][0])
    return ''.join(out)


def array_stats(seg, p):
    n_m = n_match = n_i = n_d = 0
    for u in seg.units:
        for kind, j, ch in u.events:
            if kind == 'M':
                n_m += 1
            elif kind == 'I':
                n_i += 1
            else:
                n_d += 1
    return n_m, n_i, n_d


class Model(object):
    pass


def derive_model(chm, rj, P):
    """Repeat model from the CHM13 span sequence and region.json; returns (model, None) or
    (None, reason). chm starts at span_start."""
    stratum = rj.get('stratum', '')
    ta = rj.get('tr_annotation') or {}
    try:
        span_start = int(rj['span_start'])
        c0 = int(rj['core_start']) - span_start
        c1 = int(rj['core_end']) - span_start + 1
    except (KeyError, TypeError, ValueError):
        return None, 'region.json lacks span/core coordinates'
    if stratum == 'control_nontr_sv' or ta.get('tr_class') == 'nonTR':
        return None, 'not a tandem repeat (%s)' % (stratum or ta.get('tr_class'))
    c0, c1 = max(0, c0), min(len(chm), c1)
    core = chm[c0:c1]
    if len(core) < 3 * P.min_period:
        return None, 'core too short (%d bp)' % len(core)
    cands = _period_candidates(rj, len(core))
    cands = [p for p in cands if P.min_period <= p <= P.max_period]
    if not cands:
        return None, 'no usable period (annotated period %s, motif %r)' % (rj.get('period'), (rj.get('motif') or '')[:20])
    ident = dict((p, lag_identity(core, p)) for p in cands)
    best = max(ident.values())
    p = min(q for q in cands if ident[q] >= best - P.period_tol)
    m = Model()
    m.period_candidates = dict((str(q), round(ident[q], 4)) for q in cands)
    m.period_lag_identity = round(ident[p], 4)
    motif = _best_window(chm, c0, c1, p)
    margin = max(2 * p, 50)
    w0, w1 = max(0, c0 - margin), min(len(chm), c1 + margin)
    win = chm[w0:w1]
    history = []
    seg = None
    for it in range(8):
        seg = segment(win, '', motif, '', P, True, True, free_entry=True)
        new = consensus_unit(seg.units, len(motif))
        history.append(motif)
        if new == motif or len(new) < P.min_period or new in history:
            break
        motif = new
    # rotate so the array starts at unit position 0
    first = seg.units[0] if seg.units else None
    j0 = 0
    if first is not None:
        j0 = next((j for kind, j, ch in first.events if kind in ('M', 'D')), 0)
    motif = motif[j0:] + motif[:j0]
    seg = segment(win, '', motif, '', P, True, True, free_entry=True)
    a, b = seg.a + w0, seg.b + w0
    n_m, n_i, n_d = array_stats(seg, len(motif))
    matches = 0
    for u in seg.units:
        for kind, j, ch in u.events:
            if kind == 'M' and motif[j] == ch:
                matches += 1
    alen = n_m + n_i + n_d
    m.period = len(motif)
    m.period_chosen = p
    m.motif = motif
    m.motif_annotated = rj.get('motif')
    m.rotation = j0
    m.consensus_iterations = len(history)
    m.a, m.b = a, b
    m.chm_len = len(chm)
    m.identity = round(matches / float(alen), 4) if alen else 0.0
    m.units_chm13 = sum(1 for u in seg.units if len(u.bases) * 2 >= len(motif))
    ov = max(0, min(b, c1) - max(a, c0))
    m.core_cov = round(ov / float(c1 - c0), 4) if c1 > c0 else 0.0
    m.array_start_chm13 = span_start + a        # 1-based CHM13 coordinate of the first array base
    m.array_end_chm13 = span_start + b - 1
    m.core = (c0, c1)
    reasons = []
    if m.units_chm13 < P.min_units:
        reasons.append('%d CHM13 units < %d' % (m.units_chm13, P.min_units))
    if m.identity < P.min_identity:
        reasons.append('CHM13 unit identity %.2f < %.2f' % (m.identity, P.min_identity))
    if m.core_cov < P.min_core_cov:
        reasons.append('array covers %.2f of the core < %.2f' % (m.core_cov, P.min_core_cov))
    if m.period < P.min_period:
        reasons.append('period %d < %d' % (m.period, P.min_period))
    m.reject = '; '.join(reasons) or None
    return m, None


def model_dict(m):
    return collections.OrderedDict([
        ('period', m.period), ('period_chosen', m.period_chosen), ('motif', m.motif),
        ('motif_annotated', m.motif_annotated), ('rotation', m.rotation),
        ('consensus_iterations', m.consensus_iterations),
        ('period_lag_identity', m.period_lag_identity), ('period_candidates', m.period_candidates),
        ('chm13_array_span_offsets', [m.a, m.b]), ('chm13_array', [m.array_start_chm13, m.array_end_chm13]),
        ('chm13_array_bp', m.b - m.a), ('chm13_units', m.units_chm13), ('chm13_identity', m.identity),
        ('core_cov', m.core_cov), ('reject', m.reject)])


# ---------------------------------------------------------------- mafft on symbols

_UNIT_MODES = {
    'linsi': ['--localpair', '--maxiterate', '1000'],
    'einsi': ['--genafpair', '--maxiterate', '1000'],
    'ginsi': ['--globalpair', '--maxiterate', '1000'],
    'nwnsi': ['--retree', '2', '--maxiterate', '1000'],
    'nwns2': ['--retree', '2', '--maxiterate', '0'],
}


def pick_unit_mode(mode, n, maxlen):
    if mode != 'auto':
        return mode
    # all-pairs modes cost ~ n^2/2 * L^2 cells; keep them under ~2e10
    if n * (n - 1) / 2.0 * maxlen * maxlen <= 2e10 and n <= 200:
        return 'linsi'
    if n <= 600:
        return 'nwnsi'
    return 'nwns2'


def score_matrix(reps, formula):
    d = edit_matrix(reps)
    n = len(reps)
    S = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            la, lb = len(reps[i]), len(reps[j])
            if formula == 'sp':
                S[i][j] = la + lb - d[i][j]
            else:
                S[i][j] = int(round(10.0 * (1.0 - 2.0 * d[i][j] / float(max(la, lb, 1)))))
    return S


def mafft_symbols(strings, reps, P, threads, workdir):
    """Align symbol strings (bytes) with mafft --textmatrix; returns aligned bytes rows."""
    wd = tempfile.mkdtemp(prefix='units.', dir=workdir)
    try:
        used = sorted(set(b for s in strings for b in s))
        S = score_matrix(reps, P.matrix)
        mfile = os.path.join(wd, 'matrix.txt')
        with open(mfile, 'w') as f:
            for x in used:
                for y in used:
                    f.write('0x%02x 0x%02x %d\n' % (x, y, S[SYMBOLS.index(x)][SYMBOLS.index(y)]))
        infile = os.path.join(wd, 'in.txt')
        with open(infile, 'wb') as f:
            for k, s in enumerate(strings):
                f.write(b'>s%d\n' % k + s + b'\n')
        mode = pick_unit_mode(P.unit_mafft, len(strings), max(len(s) for s in strings))
        args = list(_UNIT_MODES[mode])
        if P.unit_op is not None:
            args += ['--op', str(P.unit_op)]
        if P.unit_ep is not None:
            args += ['--ep', str(P.unit_ep)]
        cmd = [config.MAFFT, '--textmatrix', mfile] + args + ['--thread', str(threads), '--threadit', '0', infile]
        env = dict(config.tool_env(), TMPDIR=wd, LC_ALL='C', LANG='C')
        env.pop('LC_CTYPE', None)
        t = time.time()
        with open(os.path.join(wd, 'out.txt'), 'wb') as o, open(os.path.join(wd, 'log.txt'), 'wb') as e:
            p = subprocess.run(cmd, stdout=o, stderr=e, env=env, cwd=wd, timeout=P.timeout)
        secs = time.time() - t
        if p.returncode:
            with open(os.path.join(wd, 'log.txt'), 'rb') as f:
                raise RuntimeError('mafft --textmatrix failed: %s' % f.read()[-300:])
        with open(os.path.join(wd, 'out.txt'), 'rb') as f:
            data = f.read()
        rows = {}
        name, buf = None, []
        for line in data.split(b'\n'):
            if line.startswith(b'>'):
                if name is not None:
                    rows[name] = b''.join(buf)
                name, buf = line[1:].split()[0].decode(), []
            elif line:
                buf.append(line)
        if name is not None:
            rows[name] = b''.join(buf)
        out = []
        for k, s in enumerate(strings):
            r = rows.get('s%d' % k)
            if r is None or r.replace(b'-', b'') != s:
                raise RuntimeError('mafft --textmatrix output row s%d missing or wrong (LC_ALL?)' % k)
            out.append(r)
        if len(set(len(r) for r in out)) != 1:
            raise RuntimeError('mafft --textmatrix rows differ in length')
        return out, {'mode': mode, 'args': args, 'seconds': round(secs, 2)}
    finally:
        shutil.rmtree(wd, ignore_errors=True)



# ---------------------------------------------------------------- exact sum-of-pairs unit aligner

class SymCost(object):
    """Cost matrix over unit symbols 0..S-1 plus the gap G = S: C[a][b] = edit distance of the unit
    strings, C[a][G] = C[G][a] = unit length, C[G][G] = 0 (the base-level edit cost of the pair)."""

    def __init__(self, reps):
        S = len(reps)
        self.S, self.G, self.S1 = S, S, S + 1
        d = edit_matrix(reps) if S else []
        flat = [0.0] * (self.S1 * self.S1)
        for a in range(S):
            row = d[a]
            base = a * self.S1
            for b in range(S):
                flat[base + b] = float(row[b])
            flat[base + S] = float(len(reps[a]))
            flat[S * self.S1 + a] = float(len(reps[a]))
        self.flat = flat
        self.c = (ctypes.c_double * len(flat))(*flat)

    def cost(self, a, b):
        return self.flat[a * self.S1 + b]


def _cint(a):
    """ctypes view of an array('i') (no copy)."""
    if len(a) == 0:
        return (ctypes.c_int32 * 1)()
    return (ctypes.c_int32 * len(a)).from_buffer(a)


def _cdbl(xs):
    return (ctypes.c_double * max(1, len(xs)))(*[float(x) for x in xs])


def _zeros(n):
    return array('i', bytes(4 * n))


def take_rows(R, L, idx):
    """Rows idx of a flat MSA (array('i'), L columns) as a new flat MSA."""
    out = array('i')
    for i in idx:
        out.extend(R[i * L:(i + 1) * L])
    return out


def drop_gap_cols(R, n, L):
    out = _zeros(n * L)
    k = lib().drop_gap_cols(_cint(R), n, L, _cint(out))
    return out[:n * k], k


def align_flat(A, nA, LA, wA, B, nB, LB, wB, sc, go=0.0):
    """Merge flat MSAs A and B optimally under the SP cost (linear gaps, plus the approximate
    gap-open cost go > 0). Returns (merged flat MSA with A's rows then B's, columns, DP cost)."""
    ops = (ctypes.c_uint8 * (LA + LB + 1))()
    cost = ctypes.c_double()
    n = lib().profile_align_flat(_cint(A), nA, LA, _cdbl(wA), _cint(B), nB, LB, _cdbl(wB), sc.c, sc.S1,
                                 float(go), ops, ctypes.byref(cost))
    if n < 0:
        raise MemoryError('profile_align_flat: allocation failed (%d x %d columns)' % (LA, LB))
    out = _zeros((nA + nB) * n)
    lib().merge_flat(_cint(A), nA, LA, _cint(B), nB, LB, ops, n, _cint(out))
    return out, n, cost.value


def sp_total(R, n, L, weights, sc, go=0.0):
    """Exact weighted sum-of-pairs cost of a flat MSA (go per gap run and pair)."""
    return lib().sp_affine(_cint(R), n, L, _cdbl(weights), sc.c, sc.S1, float(go))


def cross_total(R, L, ia, ib, weights, sc, go=0.0):
    """Exact cost summed over pairs (a in ia, b in ib) of a flat MSA; weights by row of R."""
    a = (ctypes.c_int32 * max(1, len(ia)))(*ia)
    b = (ctypes.c_int32 * max(1, len(ib)))(*ib)
    return lib().cross_cost(_cint(R), L, a, len(ia), b, len(ib), _cdbl(weights), sc.c, sc.S1, float(go))


def sp_cost_py(rows, weights, sc):
    """Linear-gap weighted sum of pairs of an MSA given as lists of rows (reference for tests)."""
    if not rows:
        return 0.0
    tot = 0.0
    G = sc.G
    for c in range(len(rows[0])):
        cnt = collections.defaultdict(float)
        for r, w in zip(rows, weights):
            x = r[c]
            cnt[G if x < 0 else x] += w
        items = list(cnt.items())
        for x in range(len(items)):
            a, na = items[x]
            for y in range(x + 1, len(items)):
                b, nb = items[y]
                tot += na * nb * sc.cost(a, b)
    return tot


def _flat_seqs(strings):
    off = (ctypes.c_int32 * (len(strings) + 1))()
    buf = []
    for i, s in enumerate(strings):
        off[i] = len(buf)
        buf.extend(s)
    off[len(strings)] = len(buf)
    return (ctypes.c_int32 * max(1, len(buf)))(*buf), off


def guide_distances(strings, sc, exact_max):
    n = len(strings)
    seqs, off = _flat_seqs(strings)
    out = (ctypes.c_double * (n * n))()
    maxl = max(len(x) for x in strings)
    if n <= exact_max and n * (n - 1) / 2.0 * maxl * maxl <= 3e9:
        lib().unit_pair_costs(seqs, off, n, sc.c, sc.S1, out)
        kind = 'exact'
    else:
        lib().kmer_dists(seqs, off, n, 3, out)
        kind = 'kmer3'
    return [list(out[i * n:(i + 1) * n]) for i in range(n)], kind


def upgma_order(D, weights):
    """UPGMA merges [(a, b, new)] over clusters 0..n-1 (new ids from n on); lazy heap."""
    import heapq
    n = len(D)
    size = dict((i, float(weights[i])) for i in range(n))
    dist = dict((i, {}) for i in range(n))
    heap = []
    for i in range(n):
        for j in range(i + 1, n):
            dist[i][j] = dist[j][i] = D[i][j]
            heap.append((D[i][j], i, j))
    heapq.heapify(heap)
    active = set(range(n))
    merges, nxt = [], n
    while len(active) > 1:
        d, i, j = heapq.heappop(heap)
        if i not in active or j not in active:
            continue
        k = nxt
        nxt += 1
        active.discard(i)
        active.discard(j)
        dist[k] = {}
        for m in active:
            v = (size[i] * dist[i][m] + size[j] * dist[j][m]) / (size[i] + size[j])
            dist[k][m] = v
            dist[m][k] = v
            heapq.heappush(heap, (v, min(k, m), max(k, m)))
        size[k] = size[i] + size[j]
        for m in active:
            del dist[m][i]
            del dist[m][j]
        del dist[i], dist[j]
        active.add(k)
        merges.append((i, j, k))
    return merges


def prog_unit_msa(strings, weights, sc, P):
    """Progressive alignment of symbol strings (UPGMA guide tree) under the exact SP cost, then
    refinement: every guide-tree split (and every single row) is cut out, both sides are
    realigned, and the result is kept when the cost of the pairs across the split falls (the
    pairs within each side are unchanged). Returns (rows in input order as lists, stats)."""
    n = len(strings)
    t0 = time.time()
    st = collections.OrderedDict()
    if n == 1:
        return [list(strings[0])], {'note': 'one string'}
    go = P.unit_go
    D, kind = guide_distances(strings, sc, P.tree_exact_max)
    st['guide'] = kind
    merges = upgma_order(D, weights)
    prof = dict((i, ([i], array('i', strings[i]), len(strings[i]))) for i in range(n))
    for a, b, k in merges:
        ia, Ra, La = prof.pop(a)
        ib, Rb, Lb = prof.pop(b)
        R, L, _ = align_flat(Ra, len(ia), La, [weights[x] for x in ia], Rb, len(ib), Lb,
                             [weights[x] for x in ib], sc, go)
        prof[k] = (ia + ib, R, L)
    (order, R, L), = prof.values()
    pos = dict((x, r) for r, x in enumerate(order))
    msa = take_rows(R, L, [pos[i] for i in range(n)])
    st['progressive_s'] = round(time.time() - t0, 2)
    cost = sp_total(msa, n, L, weights, sc, go)
    st['sp_progressive'] = cost
    t1 = time.time()
    rounds = accepted = 0
    leaves = dict((i, [i]) for i in range(n))
    for a, b, k in merges:
        leaves[k] = leaves[a] + leaves[b]
    splits = []
    if P.refine == 'tree' and n <= P.refine_tree_max:
        splits = [leaves[k] for a, b, k in reversed(merges[:-1])]
        splits += [leaves[a] for a, b, k in merges[-1:]]
        splits = [g for g in splits if 1 < len(g) < n]
    if P.refine in ('tree', 'loo'):
        splits += [[i] for i in range(n)]
    st['refine'] = P.refine
    st['refine_splits'] = len(splits)
    timed_out = False
    for rnd in range(P.refine_rounds if splits else 0):
        improved = False
        for grp in splits:
            if time.time() - t1 > P.refine_max_s:
                timed_out = True
                break
            gs = set(grp)
            others = [i for i in range(n) if i not in gs]
            old = cross_total(msa, L, grp, others, weights, sc, go)
            A, LA = drop_gap_cols(take_rows(msa, L, grp), len(grp), L)
            B, LB = drop_gap_cols(take_rows(msa, L, others), len(others), L)
            wo, wg = [weights[i] for i in others], [weights[i] for i in grp]
            new, Ln, _ = align_flat(B, len(others), LB, wo, A, len(grp), LA, wg, sc, go)
            no = len(others)
            newc = cross_total(new, Ln, list(range(no, n)), list(range(no)), wo + wg, sc, go)
            if newc < old - 1e-6:
                p2 = dict((x, r) for r, x in enumerate(others + list(grp)))
                msa = take_rows(new, Ln, [p2[i] for i in range(n)])
                L = Ln
                cost += newc - old
                improved = True
                accepted += 1
        rounds += 1
        if not improved or timed_out:
            break
    st['refine_rounds'] = rounds
    st['refine_accepted'] = accepted
    st['refine_timed_out'] = timed_out
    st['sp_final'] = sp_total(msa, n, L, weights, sc, go)
    st['refine_s'] = round(time.time() - t1, 2)
    rows = [list(msa[i * L:(i + 1) * L]) for i in range(n)]
    return rows, st


def unit_level_msa(segs, P, threads, workdir):
    """Symbols, then the symbol-string MSA. Returns (aligned rows per seg as lists of unit index
    or None per unit column, stats)."""
    counts = collections.Counter(u.seq for s in segs for u in s.units)
    distinct = sorted(counts, key=lambda t: (-counts[t], len(t), t))
    st = collections.OrderedDict()
    st['unit_tokens'] = sum(counts.values())
    st['unit_variants'] = len(distinct)
    maxs = min(P.max_symbols, MAX_SYMBOLS) if P.unit_aligner == 'mafft' else len(distinct)
    rep = {}
    if len(distinct) <= maxs:
        kept = distinct
        for t in distinct:
            rep[t] = t
        st['symbols'] = len(kept)
        st['variants_hashed'] = 0
    else:
        kept = distinct[:maxs]
        rest = distinct[maxs:]
        for t in kept:
            rep[t] = t
        d = edit_cross(rest, kept)
        eds = []
        for i, t in enumerate(rest):
            row = d[i]
            j = min(range(len(kept)), key=lambda k: (row[k], k))
            rep[t] = kept[j]
            eds.append(row[j] / float(max(len(t), len(kept[j]), 1)))
        st['symbols'] = len(kept)
        st['variants_hashed'] = len(rest)
        st['tokens_hashed'] = sum(counts[t] for t in rest)
        st['hashed_median_rel_ed'] = round(sorted(eds)[len(eds) // 2], 3)
        st['hashed_max_rel_ed'] = round(max(eds), 3)
    if P.unit_aligner == 'prog':
        sid = dict((t, i) for i, t in enumerate(kept))
        strings = [tuple(sid[u.seq] for u in s.units) for s in segs]
    else:
        sym = dict((t, SYMBOLS[i]) for i, t in enumerate(kept))
        strings = [bytes(sym[rep[u.seq]] for u in s.units) for s in segs]
    uniq = collections.OrderedDict()
    for k, s in enumerate(strings):
        if s:
            uniq.setdefault(s, []).append(k)
    st['unit_strings_distinct'] = len(uniq)
    st['units_per_seq'] = [min(len(s.units) for s in segs), sorted(len(s.units) for s in segs)[len(segs) // 2],
                           max(len(s.units) for s in segs)]
    if not uniq:
        return [[] for _ in segs], st
    keys = list(uniq)
    st['unit_aligner'] = P.unit_aligner
    if len(keys) == 1:
        aligned = [keys[0]]
        st['unit_msa'] = {'mode': 'none', 'note': 'one distinct unit string'}
    elif P.unit_aligner == 'prog':
        sc = SymCost(kept)
        rows, info = prog_unit_msa(keys, [len(uniq[k]) for k in keys], sc, P)
        aligned = rows
        st['unit_msa'] = info
    else:
        aligned, info = mafft_symbols(keys, kept, P, threads, workdir)
        st['unit_msa'] = info
    ncol = len(aligned[0])
    st['unit_columns'] = ncol
    arow = {}
    for key, row in zip(keys, aligned):
        arow[key] = row
    out = []
    for k, s in enumerate(strings):
        if not s:
            out.append([None] * ncol)
            continue
        row = arow[s]
        gap = -1 if P.unit_aligner == 'prog' else 0x2d
        cols, ui = [], 0
        for ch in row:
            if ch == gap:
                cols.append(None)
            else:
                cols.append(ui)
                ui += 1
        out.append(cols)
    return out, st


# ---------------------------------------------------------------- expansion to bases

def nw_pair(a, b):
    """Unit-cost global alignment; returns aligned (a', b')."""
    n, m = len(a), len(b)
    D = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        D[i][0] = i
    for j in range(m + 1):
        D[0][j] = j
    for i in range(1, n + 1):
        ai, Di, Dp = a[i - 1], D[i], D[i - 1]
        for j in range(1, m + 1):
            x = Dp[j - 1] + (ai != b[j - 1])
            y = Dp[j] + 1
            z = Di[j - 1] + 1
            Di[j] = x if x <= y and x <= z else (y if y <= z else z)
    i, j, ra, rb = n, m, [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and D[i][j] == D[i - 1][j - 1] + (a[i - 1] != b[j - 1]):
            ra.append(a[i - 1]); rb.append(b[j - 1]); i -= 1; j -= 1
        elif i > 0 and D[i][j] == D[i - 1][j] + 1:
            ra.append(a[i - 1]); rb.append('-'); i -= 1
        else:
            ra.append('-'); rb.append(b[j - 1]); j -= 1
    return ''.join(reversed(ra)), ''.join(reversed(rb))


def star_msa(strs, weight=None):
    """MSA of a few short strings: identical lengths are stacked; otherwise a centre-star
    alignment around the most frequent (then longest) string. Returns {string: aligned}."""
    strs = list(dict.fromkeys(strs))
    if len(strs) == 1:
        return {strs[0]: strs[0]}
    if len(set(len(s) for s in strs)) == 1:
        return dict((s, s) for s in strs)
    w = weight or {}
    centre = max(strs, key=lambda s: (w.get(s, 0), len(s), s))
    # insertions relative to the centre, per centre gap position (0..len)
    pairs = {}
    gapw = [0] * (len(centre) + 1)
    for s in strs:
        if s == centre:
            continue
        ac, as_ = nw_pair(centre, s)
        pos, cur, ins = 0, 0, collections.defaultdict(str)
        for x, y in zip(ac, as_):
            if x == '-':
                ins[pos] += y
            else:
                pos += 1
        for k, v in ins.items():
            gapw[k] = max(gapw[k], len(v))
        pairs[s] = (ac, as_)
    out = {}
    for s in strs:
        if s == centre:
            parts = []
            for i, ch in enumerate(centre):
                parts.append('-' * gapw[i])
                parts.append(ch)
            parts.append('-' * gapw[len(centre)])
            out[s] = ''.join(parts)
            continue
        ac, as_ = pairs[s]
        ins = collections.defaultdict(str)
        col = {}
        pos = 0
        for x, y in zip(ac, as_):
            if x == '-':
                ins[pos] += y
            else:
                col[pos] = y
                pos += 1
        parts = []
        for i in range(len(centre)):
            parts.append(ins[i].ljust(gapw[i], '-'))
            parts.append(col[i])
        parts.append(ins[len(centre)].ljust(gapw[len(centre)], '-'))
        out[s] = ''.join(parts)
    return out


def expand_unit_column(toks, p):
    """toks: {unit string: slots}. Base-level block of one unit column, through the unit's
    positions; returns {unit string: aligned string} (all-gap columns removed)."""
    match = {}
    ins = {}
    for t, slots in toks.items():
        mm, ii = {}, collections.defaultdict(list)
        for ch, sl in zip(t, slots):
            j, isins = divmod(sl, 2)
            if isins:
                ii[j].append(ch)
            else:
                mm[j] = ch
        match[t] = mm
        ins[t] = dict((j, ''.join(v)) for j, v in ii.items())
    ins_al = {}
    for j in range(p):
        present = [ins[t][j] for t in toks if j in ins[t]]
        if present:
            wt = collections.Counter(present)
            al = star_msa(present, wt)
            width = len(next(iter(al.values())))
            ins_al[j] = (al, width)
    out = {}
    for t in toks:
        parts = []
        for j in range(p):
            parts.append(match[t].get(j, '-'))
            if j in ins_al:
                al, width = ins_al[j]
                s = ins[t].get(j)
                parts.append(al[s] if s is not None else '-' * width)
        out[t] = ''.join(parts)
    # drop columns that are gaps in every token of this unit column
    L = len(next(iter(out.values())))
    keep = [c for c in range(L) if any(v[c] != '-' for v in out.values())]
    if len(keep) < L:
        out = dict((t, ''.join(v[c] for c in keep)) for t, v in out.items())
    return out


def expand(segs, colmaps, p):
    """Base rows of the array part of every seg from the unit-level MSA."""
    pat = collections.defaultdict(collections.Counter)
    for s in segs:
        for u in s.units:
            pat[u.seq][tuple(u.slots)] += 1
    canon = dict((t, c.most_common(1)[0][0]) for t, c in pat.items())
    ncol = len(colmaps[0]) if colmaps else 0
    parts = [[] for _ in segs]
    widths = []
    for c in range(ncol):
        toks = {}
        for k, s in enumerate(segs):
            ui = colmaps[k][c]
            if ui is not None:
                t = s.units[ui].seq
                toks[t] = canon[t]
        if not toks:
            continue
        block = expand_unit_column(toks, p)
        w = len(next(iter(block.values())))
        widths.append(w)
        for k, s in enumerate(segs):
            ui = colmaps[k][c]
            parts[k].append(block[s.units[ui].seq] if ui is not None else '-' * w)
    return [''.join(x) for x in parts], widths



# ---------------------------------------------------------------- base-level polish

_BASE_CODE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


class BaseCost(object):
    """Unit costs over A, C, G, T, other (4) and the gap (5)."""

    def __init__(self):
        self.S, self.G, self.S1 = 5, 5, 6
        flat = [0.0] * 36
        for a in range(6):
            for b in range(6):
                flat[a * 6 + b] = 0.0 if a == b else 1.0
        self.flat = flat
        self.c = (ctypes.c_double * 36)(*flat)


def polish_msa(rows, P, band):
    """Leave-one-out base-level realignment of each row against the others (banded around its
    current path, exact SP cost with a small gap-open term); a new row placement is kept only
    when the cost of its pairs with the other rows falls. rows: aligned strings (no N expected;
    any non-ACGT is one symbol). Returns (rows, stats)."""
    t0 = time.time()
    st = collections.OrderedDict([('rounds', 0), ('accepted', 0)])
    n = len(rows)
    if n < 2 or P.polish_rounds <= 0:
        return rows, st
    L = len(rows[0])
    msa = array('i', [(-1 if ch == '-' else _BASE_CODE.get(ch, 4)) for r in rows for ch in r])
    sc = BaseCost()
    w = [1.0] * n
    go = P.polish_go
    st['sp_before'] = sp_total(msa, n, L, w, sc, 0.0)
    st['band'] = band
    timed_out = False
    for rnd in range(P.polish_rounds):
        improved = False
        for x in range(n):
            if time.time() - t0 > P.polish_max_s:
                timed_out = True
                break
            rest = _zeros((n - 1) * L)
            row = _zeros(L)
            m = ctypes.c_int32()
            lo = _zeros(L + 1)
            hi = _zeros(L + 1)
            Lr = lib().band_extract(_cint(msa), n, L, x, band, _cint(rest), _cint(row), ctypes.byref(m),
                                    _cint(lo), _cint(hi))
            m = m.value
            cells = sum(max(0, hi[i] - lo[i] + 1) for i in range(Lr + 1))
            if cells > P.polish_max_cells or m == 0:
                st['rows_skipped'] = st.get('rows_skipped', 0) + 1
                continue
            rest = rest[:(n - 1) * Lr]
            row = row[:m]
            others = [i for i in range(n) if i != x]
            old = cross_total(msa, L, [x], others, w, sc, go)
            ops = (ctypes.c_uint8 * (Lr + m + 1))()
            cost = ctypes.c_double()
            k = lib().row_to_profile_banded(_cint(rest), n - 1, Lr, _cdbl(w[:n - 1]), _cint(row), m, _cint(lo),
                                            _cint(hi), sc.c, sc.S1, float(go), ops, ctypes.byref(cost))
            if k < 0:
                continue
            new = _zeros(n * k)
            lib().merge_flat(_cint(rest), n - 1, Lr, _cint(row), 1, m, ops, k, _cint(new))
            newc = cross_total(new, k, [n - 1], list(range(n - 1)), w, sc, go)
            if newc < old - 1e-9:
                p2 = dict((q, r) for r, q in enumerate(others + [x]))
                msa = take_rows(new, k, [p2[i] for i in range(n)])
                L = k
                improved = True
                st['accepted'] += 1
        st['rounds'] += 1
        if not improved or timed_out:
            break
    st['timed_out'] = timed_out
    st['sp_after'] = sp_total(msa, n, L, w, sc, 0.0)
    st['seconds'] = round(time.time() - t0, 2)
    inv = 'ACGT'
    out = []
    for i, r in enumerate(rows):
        orig = iter(ch for ch in r if ch != '-')
        s = []
        for v in msa[i * L:(i + 1) * L]:
            if v < 0:
                s.append('-')
            else:
                ch = next(orig)
                s.append(ch)
        out.append(''.join(s))
    return out, st

# ---------------------------------------------------------------- flank pieces

def align_pieces(pieces, P, threads, workdir, label):
    """MSA rows of arbitrary pieces (in order) via realign.py's pipeline; identical and empty
    pieces are handled without mafft."""
    import realign
    distinct = list(dict.fromkeys(x for x in pieces if x))
    info = collections.OrderedDict([('pieces_distinct', len(distinct)),
                                    ('max_len', max([len(x) for x in pieces] or [0]))])
    if not distinct:
        return [''] * len(pieces), info
    if len(distinct) == 1:
        L = len(distinct[0])
        info['method'] = 'none'
        return [x if x else '-' * L for x in pieces], info
    method = P.flank_method
    if method == 'auto':
        if info['max_len'] > 20000:
            method = 'mafft_fftns2'
        elif len(distinct) > 100:
            method = 'mafft_fftnsi'
        else:
            method = 'mafft_linsi'
    wd = tempfile.mkdtemp(prefix='flank.%s.' % label, dir=workdir)
    try:
        fin = os.path.join(wd, 'in.fa')
        fout = os.path.join(wd, 'out.fa')
        with open(fin, 'w') as f:
            for k, x in enumerate(distinct):
                f.write('>p%d\n%s\n' % (k, x))
        r = realign.align_fasta(method, fin, fout, threads=threads, timeout=P.timeout, mem_mb=P.mem_mb, workdir=wd)
        if r.get('status') != 'ok' and method != 'mafft_fftns2':
            info['first_attempt'] = {'method': method, 'status': r.get('status'), 'message': r.get('message')}
            method = 'mafft_fftns2'
            r = realign.align_fasta(method, fin, fout, threads=threads, timeout=P.timeout, mem_mb=P.mem_mb, workdir=wd)
        if r.get('status') != 'ok':
            raise RuntimeError('flank alignment (%s) failed: %s %s' % (label, r.get('status'), r.get('message')))
        got = dict(msa_graph.read_msa(fout))
        info['method'] = method
        info['seconds'] = r.get('align_s')
        row = dict((x, got['p%d' % k]) for k, x in enumerate(distinct))
        L = len(next(iter(row.values())))
        return [row[x] if x else '-' * L for x in pieces], info
    finally:
        shutil.rmtree(wd, ignore_errors=True)



# ---------------------------------------------------------------- adequacy guard

def guard_pairs(seqs, P):
    """Pairs of distinct sequences for the SP / optimum ratio: all of them, or a seeded sample
    that fits guard_max_pairs and guard_max_cells."""
    import random
    n = len(seqs)
    allp = [(a, b) for a in range(n) for b in range(a + 1, n)]
    rng = random.Random(1)
    if len(allp) > P.guard_max_pairs:
        allp = rng.sample(allp, P.guard_max_pairs)
    out, cells = [], 0.0
    for a, b in allp:
        c = len(seqs[a]) * float(len(seqs[b])) / 64.0
        if out and cells + c > P.guard_max_cells / 64.0:
            break
        out.append((a, b))
        cells += c
    return out


def msa_sp_over_opt(rows_by_seq, pairs, opt=None):
    """(sum over pairs of the MSA-induced unit cost, sum of the optimal edit distances, ratio).
    rows_by_seq: [(sequence, aligned row)] of distinct sequences."""
    L = lib()
    ind = 0
    need_opt = opt is None
    opt = opt or 0
    for a, b in pairs:
        ra, rb = rows_by_seq[a][1].encode(), rows_by_seq[b][1].encode()
        ind += L.induced_cost(ra, rb, len(ra))
        if need_opt:
            sa, sb = rows_by_seq[a][0].encode(), rows_by_seq[b][0].encode()
            opt += L.ed_bitpar(sa, len(sa), sb, len(sb))
    return ind, opt, (ind / float(opt) if opt else 1.0)

# ---------------------------------------------------------------- the whole method

def fallback_method(P, n_distinct):
    """realign.py method for the fallback: P.fallback, or P.fallback_large above P.large_n sequences."""
    return P.fallback_large if (P.fallback_large and n_distinct > P.large_n) else P.fallback


def find_chm13(recs):
    for n, s in recs:
        if n.split('#')[0] == 'CHM13':
            return n, s
    return None, None


def unit_align_records(recs, rj, chm, P, threads, workdir):
    """recs: [(name, seq)] (unique names). Returns (rows [(name, aligned)], info). Raises on failure
    of the unit path (the caller falls back)."""
    import realign
    t0 = time.time()
    info = collections.OrderedDict()
    model, why = derive_model(chm.upper(), rj, P)
    if model is None:
        return None, {'reason': why}
    info['model'] = model_dict(model)
    if model.reject:
        return None, {'reason': 'motif unusable: ' + model.reject, 'model': info['model']}
    info['model_s'] = round(time.time() - t0, 2)
    seqs = [s.upper() for _, s in recs]
    idx = collections.OrderedDict()
    for s in seqs:
        idx.setdefault(s, len(idx))
    distinct = list(idx)
    row_of = [idx[s] for s in seqs]
    masked, runs = [], {}
    for k, s in enumerate(distinct):
        ms, rs = realign.mask_runs(s, 1)
        masked.append(ms)
        runs[k] = rs
    info['n_seqs'] = len(recs)
    info['n_distinct'] = len(distinct)
    info['masked_runs'] = sum(len(v) for v in runs.values())
    info['masked_bp'] = sum(len(x) for v in runs.values() for _, x in v)
    chmU = chm.upper()
    a, b = model.a, model.b
    K = P.flank_k
    fl = chmU[max(0, a - K):a]
    fr = chmU[b:b + K]
    free_p, free_s = a > K, len(chmU) - b > K
    t1 = time.time()
    # Diagnostic: CHM13's own segmentation normally reproduces the model's array. It can run past
    # it when the flank carries more repeat-like sequence than the model's search window saw (a
    # second array beyond a gap, or an AT-rich flank under a 2-bp unit): the cycle then absorbs it.
    # Stricter scores were tried for those regions and gave worse alignments, so this is only flagged.
    tol = max(2 * model.period, 20)
    chm_masked = realign.mask_runs(chmU, 1)[0]
    cs = segment(chm_masked, fl, model.motif, fr, P, free_p, free_s)
    info['chm13_dp_array'] = [cs.a, cs.b]
    info['array_beyond_model'] = bool(abs(cs.a - a) > tol or abs(cs.b - b) > tol)
    SP = P
    segs = []
    for s in masked:
        segs.append(segment(s, fl, model.motif, fr, SP, free_p, free_s))
    info['segment_s'] = round(time.time() - t1, 2)
    info['flank_anchor'] = {'left_bp': len(fl), 'right_bp': len(fr), 'free_prefix': free_p, 'free_suffix': free_s}
    arr = sorted(s.b - s.a for s in segs)
    info['array_bp'] = [arr[0], arr[len(arr) // 2], arr[-1]]
    info['no_array'] = sum(1 for s in segs if not s.units)
    # CHM13's own segmentation must put its array where the model says
    cn, _ = find_chm13(recs)
    if cn is not None:
        cs = segs[row_of[[n for n, _ in recs].index(cn)]]
        info['chm13_array_offsets'] = [cs.a, cs.b]
    t2 = time.time()
    colmaps, ust = unit_level_msa(segs, P, threads, workdir)
    info['units'] = ust
    info['unit_msa_s'] = round(time.time() - t2, 2)
    t3 = time.time()
    arows, widths = expand(segs, colmaps, model.period)
    info['expand_s'] = round(time.time() - t3, 2)
    info['array_columns'] = len(arows[0]) if arows else 0
    t4 = time.time()
    lrows, linfo = align_pieces([s.left for s in segs], P, threads, workdir, 'left')
    rrows, rinfo = align_pieces([s.right for s in segs], P, threads, workdir, 'right')
    info['flank_s'] = round(time.time() - t4, 2)
    info['left_flank'] = linfo
    info['right_flank'] = rinfo
    rows = [lrows[k] + arows[k] + rrows[k] for k in range(len(segs))]
    for k, s in enumerate(masked):
        if rows[k].replace('-', '') != s:
            raise RuntimeError('internal: unit MSA row %d does not spell its (masked) sequence' % k)
    if P.polish_rounds > 0 and len(rows) > 1:
        band = min(max(P.polish_band, 3 * model.period), 2000)
        rows, pst = polish_msa(rows, P, band)
        info['polish'] = pst
        for k, s in enumerate(masked):
            if rows[k].replace('-', '') != s:
                raise RuntimeError('internal: polished row %d does not spell its (masked) sequence' % k)
    rows = realign.reinsert_runs(rows, runs)
    if rows:
        keep = [c for c in range(len(rows[0])) if any(r[c] != '-' for r in rows)]
        if len(keep) < len(rows[0]):
            rows = [''.join(r[c] for c in keep) for r in rows]
    for k, s in enumerate(distinct):
        if rows[k].replace('-', '') != s:
            raise RuntimeError('internal: row %d does not spell its sequence after N reinsertion' % k)
    out = [(recs[i][0], rows[row_of[i]]) for i in range(len(recs))]
    info['columns'] = len(rows[0]) if rows else 0
    info['phase_check'] = phase_check(segs, colmaps)
    info['unit_s'] = round(time.time() - t0, 2)
    return out, info


def phase_check(segs, colmaps):
    """Share of unit columns whose units all start at the same unit position (in phase), and of
    all units that sit in such columns."""
    if not colmaps or not colmaps[0]:
        return None
    ncol = len(colmaps[0])
    same = tot = units_in = units_all = 0
    for c in range(ncol):
        starts = []
        for k, s in enumerate(segs):
            ui = colmaps[k][c]
            if ui is not None:
                u = s.units[ui]
                starts.append(next((j for kind, j, ch in u.events if kind in ('M', 'D')), -1))
        if not starts:
            continue
        tot += 1
        units_all += len(starts)
        if len(set(starts)) == 1:
            same += 1
            units_in += len(starts)
    return {'unit_columns_one_start': same, 'unit_columns': tot,
            'units_in_one_start_columns': round(units_in / float(units_all), 4) if units_all else None}


def align_units_fasta(input_fa, region_json, msa_out, threads=2, params=None, workdir=None, chm13=None):
    """MSA-only entry point: align every record of input_fa (unit-aware, or the fallback), write
    msa_out (only when status is 'ok'); returns the info dict."""
    import realign
    P = params or Params()
    t0 = time.time()
    if isinstance(region_json, str):
        with open(region_json) as f:
            rj = json.load(f)
    else:
        rj = dict(region_json)
    recs = msa_graph.read_fasta(input_fa)
    info = collections.OrderedDict([('method', METHOD), ('version', VERSION),
                                    ('input', os.path.abspath(input_fa)), ('msa', os.path.abspath(msa_out)),
                                    ('region_id', rj.get('region_id')), ('stratum', rj.get('stratum')),
                                    ('status', None), ('mode', None), ('reason', None)])
    info['params'] = P.as_dict()
    info['n_records'] = len(recs)
    names = [n for n, _ in recs]
    if not recs or len(set(names)) != len(names):
        info.update(status='error', message='empty input or duplicate names')
        return info
    wroot = workdir or os.path.join(config.WORK_DIR, 'realign_units')
    os.makedirs(wroot, exist_ok=True)
    wd = tempfile.mkdtemp(prefix='ua.', dir=wroot)
    try:
        chm = chm13
        cname, cseq = find_chm13(recs)
        if chm is None:
            chm = cseq
        if chm is None and isinstance(region_json, str):
            hap = os.path.join(os.path.dirname(os.path.abspath(region_json)), 'hap32.fa')
            if os.path.exists(hap):
                chm = find_chm13(msa_graph.read_fasta(hap))[1]
        info['chm13_source'] = 'input' if cseq is not None and chm13 is None else ('given' if chm13 else 'hap32.fa')
        rows, uinfo = None, {}
        if chm is None:
            uinfo = {'reason': 'no CHM13 sequence (input or hap32.fa next to region.json)'}
        else:
            try:
                rows, uinfo = unit_align_records(recs, rj, chm, P, threads, wd)
            except Exception as e:     # any failure of the unit path -> fallback, recorded
                rows = None
                uinfo = dict(uinfo or {}, reason='unit path failed: %s: %s' % (type(e).__name__, str(e)[:300]))
        info.update(uinfo)
        if rows is not None and P.guard_ratio and P.guard_ratio > 0:
            # adequacy guard: does the unit MSA cost much more than the pairwise optima?
            tg = time.time()
            byseq = collections.OrderedDict()
            for n, r in rows:
                byseq.setdefault(r.replace('-', ''), r)
            rs = list(byseq.items())
            pairs = guard_pairs([x for x, _ in rs], P)
            ind, opt, ratio = msa_sp_over_opt(rs, pairs)
            g = collections.OrderedDict([('pairs', len(pairs)), ('unit_sp', ind), ('opt_sum', opt),
                                         ('unit_sp_over_opt', round(ratio, 4)), ('threshold', P.guard_ratio)])
            if ratio > P.guard_ratio:
                fb_out = os.path.join(wd, 'guard_fallback.msa.fa')
                fbm = fallback_method(P, len(rs))
                fb = realign.align_fasta(fbm, input_fa, fb_out, threads=threads, timeout=P.timeout,
                                         mem_mb=P.mem_mb, workdir=wd)
                g['fallback_method'] = fbm
                g['fallback_status'] = fb.get('status')
                if fb.get('status') == 'ok':
                    fb_rows = dict(msa_graph.read_msa(fb_out))
                    name_of = {}
                    for n, r in rows:
                        name_of.setdefault(r.replace('-', ''), n)
                    fb_byseq = [(x, fb_rows[name_of[x]]) for x, _ in rs]
                    find, _, fratio = msa_sp_over_opt(fb_byseq, pairs, opt)
                    g['fallback_sp'] = find
                    g['fallback_sp_over_opt'] = round(fratio, 4)
                    if find < ind:
                        g['decision'] = 'fallback'
                        info['guard'] = g
                        info['unit_rejected'] = {'reason': 'guard', 'units': info.get('units'),
                                                 'model': info.get('model')}
                        rows = None
                        uinfo = {'reason': 'guard: unit MSA costs %.3fx the pairwise optima (> %.2f) and %s is '
                                           'lower (%.3fx)' % (ratio, P.guard_ratio, fbm, fratio)}
                        info['reason'] = uinfo['reason']
                        info['guard_fallback_msa'] = fb_out
                        info['fallback_info'] = fb
                    else:
                        g['decision'] = 'unit'
                else:
                    g['decision'] = 'unit'
            else:
                g['decision'] = 'unit'
            g['seconds'] = round(time.time() - tg, 2)
            info['guard'] = g
        if rows is not None:
            tmp = msa_out + '.tmp%d' % os.getpid()
            os.makedirs(os.path.dirname(os.path.abspath(msa_out)) or '.', exist_ok=True)
            msa_graph.write_msa(rows, tmp)
            os.replace(tmp, msa_out)
            info['status'] = 'ok'
            info['mode'] = 'unit'
            info['reason'] = None
            info['gap_frac'] = round(sum(r.count('-') for _, r in rows) / float(max(1, len(rows) * len(rows[0][1]))), 4)
        elif info.get('guard_fallback_msa'):
            info['mode'] = 'fallback'
            fb = info.pop('fallback_info')
            os.makedirs(os.path.dirname(os.path.abspath(msa_out)) or '.', exist_ok=True)
            shutil.move(info.pop('guard_fallback_msa'), msa_out)
            info['fallback'] = collections.OrderedDict((k, fb.get(k)) for k in (
                'method', 'status', 'message', 'n_distinct', 'align_s', 'masked_runs', 'masked_bp', 'columns'))
            info['status'] = 'ok'
            info['columns'] = fb.get('columns')
            info['gap_frac'] = fb.get('gap_frac')
        else:
            info['mode'] = 'fallback'
            fb = realign.align_fasta(fallback_method(P, len(set(s.upper() for _, s in recs))), input_fa, msa_out,
                                     threads=threads, timeout=P.timeout, mem_mb=P.mem_mb, workdir=wd)
            info['fallback'] = collections.OrderedDict((k, fb.get(k)) for k in (
                'method', 'status', 'message', 'n_distinct', 'align_s', 'masked_runs', 'masked_bp', 'columns'))
            info['status'] = fb.get('status')
            info['columns'] = fb.get('columns')
            info['gap_frac'] = fb.get('gap_frac')
            if fb.get('status') != 'ok':
                info['message'] = fb.get('message')
        return info
    finally:
        info['align_s'] = round(time.time() - t0, 2)
        shutil.rmtree(wd, ignore_errors=True)


# ---------------------------------------------------------------- region mode

def outputs(rid, method=METHOD):
    d = os.path.join(config.CANDIDATES_DIR, method)
    return {'dir': d, 'msa': os.path.join(d, rid + '.msa.fa'), 'gfa': os.path.join(d, rid + '.gfa'),
            'json': os.path.join(d, rid + '.realign.json')}


def realign_region(rdir, threads=2, params=None, force=False, method=METHOD, record_runtime=True):
    import realign
    P = params or Params()
    rid = os.path.basename(os.path.normpath(rdir))
    o = outputs(rid, method)
    if not force and os.path.exists(o['json']) and os.path.exists(o['gfa']) and os.path.exists(o['msa']):
        try:
            if json.load(open(o['json'])).get('status') == 'ok':
                return {'method': method, 'region_id': rid, 'status': 'skipped', 'note': 'kept ok'}
        except ValueError:
            pass
    os.makedirs(o['dir'], exist_ok=True)
    for k in ('msa', 'gfa'):
        if os.path.exists(o[k]):
            os.remove(o[k])
    t0 = time.time()
    hap = os.path.join(rdir, 'hap32.fa')
    rjp = os.path.join(rdir, 'region.json')
    tmp_msa = o['msa'] + '.tmp%d' % os.getpid()
    wroot = os.path.join(config.WORK_DIR, 'realign_units')
    info = align_units_fasta(hap, rjp, tmp_msa, threads=threads, params=P, workdir=wroot)
    info['method'] = method
    info['msa'] = os.path.relpath(o['msa'], config.REPO) if info['status'] == 'ok' else None
    graph_s = None
    if info['status'] == 'ok':
        tg = time.time()
        tmp_gfa = o['gfa'] + '.tmp%d' % os.getpid()
        try:
            st = msa_graph.msa_to_gfa(tmp_msa, hap, tmp_gfa, workdir=wroot)
            os.replace(tmp_msa, o['msa'])
            os.replace(tmp_gfa, o['gfa'])
            st['msa'] = os.path.relpath(o['msa'], config.REPO)
            st['gfa'] = os.path.relpath(o['gfa'], config.REPO)
            st['hap32'] = os.path.relpath(hap, config.REPO)
            info['graph'] = st
            info['gfa'] = st['gfa']
        except msa_graph.MsaGraphError as e:
            info['status'] = 'graph_error'
            info['message'] = str(e)
        graph_s = round(time.time() - tg, 2)
        for pth in (tmp_msa, tmp_gfa):
            if os.path.exists(pth):
                os.remove(pth)
    info['graph_s'] = graph_s
    info['total_s'] = round(time.time() - t0, 2)
    info['threads'] = threads
    info['finished'] = datetime.datetime.now().isoformat(timespec='seconds')
    info['host_load'] = round(os.getloadavg()[0], 1)
    with open(o['json'] + '.tmp', 'w') as f:
        json.dump(info, f, indent=1)
    os.replace(o['json'] + '.tmp', o['json'])
    g = info.get('graph') or {}
    try:
        rj = json.load(open(rjp))
        span = int(rj['span_end']) - int(rj['span_start']) + 1
    except (OSError, ValueError, KeyError):
        span = ''
    note = info.get('mode') or ''
    if info.get('reason'):
        note += ': ' + info['reason']
    if info.get('message'):
        note += ' | ' + str(info['message'])
    u = info.get('units') or {}
    row = {'method': method, 'region_id': rid, 'stratum': info.get('stratum'), 'span_bp': span,
           'n_seqs': info.get('n_records'),
           'n_distinct': info.get('n_distinct') or (info.get('fallback') or {}).get('n_distinct'),
           'status': info['status'], 'align_s': info.get('align_s'), 'graph_s': graph_s,
           'total_s': info['total_s'], 'threads': threads, 'timeout_s': P.timeout,
           'columns': info.get('columns', ''), 'gap_frac': info.get('gap_frac', ''), 'nodes': g.get('nodes', ''),
           'nodes_per_kb': g.get('nodes_per_kb', ''), 'frac_nodes_1bp': g.get('frac_nodes_1bp', ''),
           'masked_runs': info.get('masked_runs', ''), 'masked_bp': info.get('masked_bp', ''),
           'N_nodes': g.get('N_nodes', ''), 'finished': info['finished'],
           'note': (note + (' | symbols %s/%s' % (u.get('symbols'), u.get('unit_variants')) if u else ''))[:200]
           .replace('\t', ' ').replace('\n', ' ')}
    if record_runtime:
        try:
            realign._lock_update_runtime(row)
        except Exception as e:  # the table is a convenience; never lose the candidate over it
            sys.stderr.write('realign_units: runtime table not updated: %s\n' % e)
    row['mode'] = info.get('mode')
    row['reason'] = info.get('reason')
    row['symbols'] = u.get('symbols')
    row['unit_variants'] = u.get('unit_variants')
    return row


def _job(args):
    rd, kw = args
    try:
        return realign_region(rd, **kw)
    except Exception as e:
        return {'method': kw.get('method', METHOD), 'region_id': os.path.basename(rd), 'status': 'crash',
                'note': '%s: %s' % (type(e).__name__, e)}


def params_from_args(a):
    kw = {}
    if a.dp_scores:
        ma, mi, go, ge = [int(x) for x in a.dp_scores.split(',')]
        kw.update(ma=ma, mi=mi, go=go, ge=ge)
    for k in ('fallback_engine', 'fallback_large', 'large_n', 'flank_k', 'unit_aligner', 'refine_rounds', 'polish_rounds', 'refine', 'unit_go', 'unit_mafft', 'unit_op', 'unit_ep', 'matrix', 'max_symbols', 'fallback', 'timeout',
              'mem_mb', 'flank_method', 'guard_ratio', 'min_units', 'min_identity', 'min_core_cov', 'period_tol'):
        v = getattr(a, k, None)
        if v is not None:
            kw[k] = v
    return Params(**kw)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('regions', nargs='*', help="region ids or dirs; 'all' = every packaged region")
    ap.add_argument('--input', help='MSA-only mode: align this FASTA ...')
    ap.add_argument('--region-json', help='... with the motif / reference array of this region.json ...')
    ap.add_argument('--msa-out', help='... and write only the MSA here')
    ap.add_argument('--info-out', help='MSA-only mode: also write the info JSON here')
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--stratum', help='comma list of strata (region mode)')
    ap.add_argument('--method-name', default=METHOD, help='candidates/<name>/ (default unit_aware)')
    ap.add_argument('--no-runtime', action='store_true', help='do not update results/realign_runtime.tsv')
    ap.add_argument('--dp-scores', help='match,mismatch,open,extend for the unit DP (default 2,4,6,1)')
    ap.add_argument('--flank-k', type=int)
    ap.add_argument('--unit-aligner', choices=['prog', 'mafft'], help='unit-level aligner (default prog)')
    ap.add_argument('--refine-rounds', type=int)
    ap.add_argument('--polish-rounds', type=int, help='base-level polish rounds (default 2, 0 = off)')
    ap.add_argument('--refine', choices=['tree', 'loo', 'none'], help='unit-level refinement (default tree)')
    ap.add_argument('--unit-go', type=float, help='unit-level gap-open cost (prog)')
    ap.add_argument('--unit-mafft', choices=sorted(_UNIT_MODES) + ['auto'])
    ap.add_argument('--unit-op', type=float)
    ap.add_argument('--unit-ep', type=float)
    ap.add_argument('--matrix', choices=['norm', 'sp'])
    ap.add_argument('--max-symbols', type=int)
    ap.add_argument('--fallback-engine', choices=list(FALLBACK_ENGINES),
                    help='aligner for the flank pieces, the fallback and the guard fallback (default mafft); '
                         '--flank-method, --fallback and --fallback-large override it')
    ap.add_argument('--flank-method', help='realign.py method for flank pieces (default auto)')
    ap.add_argument('--fallback', help='realign.py method for the fallback (default mafft_linsi)')
    ap.add_argument('--fallback-large', help='fallback above --large-n distinct sequences (default mafft_fftnsi)')
    ap.add_argument('--large-n', type=int, help='distinct-sequence count above which --fallback-large is used (100)')
    ap.add_argument('--guard-ratio', type=float, help='adequacy guard threshold (default 1.2; 0 = off)')
    ap.add_argument('--timeout', type=int)
    ap.add_argument('--mem-mb', type=int)
    ap.add_argument('--min-units', type=int)
    ap.add_argument('--min-identity', type=float)
    ap.add_argument('--min-core-cov', type=float)
    ap.add_argument('--period-tol', type=float)
    a = ap.parse_args(argv)
    P = params_from_args(a)
    if a.input:
        if not (a.msa_out and a.region_json):
            ap.error('--input needs --msa-out and --region-json')
        info = align_units_fasta(a.input, a.region_json, a.msa_out, threads=a.threads, params=P)
        if a.info_out:
            with open(a.info_out, 'w') as f:
                json.dump(info, f, indent=1)
        sys.stderr.write('realign_units: %s (%s%s) %s s\n' % (info['status'], info.get('mode'),
                                                              ': ' + info['reason'] if info.get('reason') else '',
                                                              info.get('align_s')))
        return 0 if info['status'] == 'ok' else 1
    import realign
    if not a.regions:
        ap.error('give region ids, all, or --input')
    if a.regions == ['all']:
        rows = realign.region_list()
        if a.stratum:
            keep = set(a.stratum.split(','))
            rows = [r for r in rows if r.get('stratum') in keep]
        dirs = [os.path.join(config.REGIONS_DIR, r['region_id']) for r in rows]
    else:
        dirs = [realign.region_dir(x) for x in a.regions]
    kw = dict(threads=a.threads, params=P, force=a.force, method=a.method_name, record_runtime=not a.no_runtime)
    jobs = [(d, kw) for d in dirs]
    t0 = time.time()
    n_ok = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, a.jobs)) as ex:
        for r in ex.map(_job, jobs):
            n_ok += r.get('status') in ('ok', 'skipped')
            sys.stderr.write('%s\t%s\t%s\t%s\t%ss\t%s\n' % (r.get('region_id'), r.get('status'), r.get('mode', ''),
                                                        'sym %s/%s' % (r.get('symbols'), r.get('unit_variants'))
                                                        if r.get('symbols') is not None else '',
                                                        r.get('total_s', ''), (r.get('reason') or r.get('note') or '')[:120]))
    sys.stderr.write('realign_units: %d/%d ok in %.0f s\n' % (n_ok, len(jobs), time.time() - t0))
    return 0 if n_ok == len(jobs) else 1


if __name__ == '__main__':
    sys.exit(main())
