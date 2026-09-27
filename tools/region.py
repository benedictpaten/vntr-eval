#!/usr/bin/env python3
"""region.py -- everything about one CHM13 interval of the hap32 graph.

    python3 tools/region.py <contig> <start> <end> <outdir> [--pad N] [--no-reads]
                            [--names] [--anchors L+:R+] [--mode between|interval]
                            [--context N] [--msa] [--bandage]

<start>/<end> are 1-based inclusive (VCF convention). tools/region.md lists every output
file and column. Data and tool locations come from config.py.

Standard library only; shells out to gbz-base, bcftools, samtools, bgzip and
(optionally) vg, mafft and Bandage.

The subgraph is NOT the gbz-base `--interval --context 0` subgraph. That subgraph contains
only the nodes of the CHM13 walk, so every haplotype that deviates from CHM13 leaves it and
is chopped into fragments. Instead (mode `between`, the default) two CHM13 anchor nodes are
picked flanking the padded interval: every haplotype walks through each exactly once and no
truth variant touches them. Everything between them is extracted with
`gbz-base query --between`. The analysed span (anchor to anchor, both anchors included)
therefore usually extends beyond the requested interval; both are reported. `--anchors`
fixes the two anchor nodes instead of choosing them, so that the output lines up with a
region package (regions/<id>/region.json: anchor_left, anchor_right).
"""

import argparse
import collections
import json
import os
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402

# ----------------------------------------------------------------------------
# Locations (all from config.py)

CACHE_DIR = config.CACHE_DIR
GBZ_BASE = config.GBZ_BASE
BCFTOOLS = config.BCFTOOLS
SAMTOOLS = config.SAMTOOLS
BGZIP = config.BGZIP
MAFFT = config.MAFFT
BANDAGE = config.BANDAGE
VG = config.VG
data_paths = config.data_paths


# ----------------------------------------------------------------------------
# Small utilities

class ToolError(RuntimeError):
    pass


class Log:
    def __init__(self, path=None, verbose=True):
        self.fh = open(path, 'w') if path else None
        self.verbose = verbose
        self.t0 = time.time()

    def write(self, msg):
        if self.fh:
            self.fh.write(msg if msg.endswith('\n') else msg + '\n')
            self.fh.flush()

    def info(self, msg):
        line = '[%6.1fs] %s' % (time.time() - self.t0, msg)
        self.write(line)
        if self.verbose:
            print(line, file=sys.stderr, flush=True)


NULL_LOG = Log(None, verbose=False)


def run(cmd, log=NULL_LOG, check=True, text=True, stdin_data=None, stdout_path=None):
    """Run a command; return CompletedProcess.  stderr is kept in the log."""
    log.write('$ ' + ' '.join(shlex.quote(str(c)) for c in cmd))
    env = config.tool_env()
    if stdout_path:
        with open(stdout_path, 'w' if text else 'wb') as out:
            p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=text,
                               input=stdin_data, env=env)
    else:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=text,
                           input=stdin_data, env=env)
    err = p.stderr if text else p.stderr.decode(errors='replace')
    if err:
        log.write('  stderr: ' + err.strip().replace('\n', '\n  stderr: '))
    if check and p.returncode != 0:
        raise ToolError('command failed (%d): %s\n%s' % (p.returncode, ' '.join(map(str, cmd)), err))
    return p


_RC = str.maketrans('ACGTNacgtnRYKMSWBDHVrykmswbdhv', 'TGCANtgcanYRMKSWVHDByrmkswvhdb')


def revcomp(s):
    return s.translate(_RC)[::-1]


def rle(ops):
    out = []
    prev, n = None, 0
    for o in ops:
        if o == prev:
            n += 1
        else:
            if prev is not None:
                out.append('%d%s' % (n, prev))
            prev, n = o, 1
    if prev is not None:
        out.append('%d%s' % (n, prev))
    return ''.join(out)


def write_tsv(path, header, rows):
    with open(path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for r in rows:
            f.write('\t'.join('.' if r.get(h) is None or r.get(h) == '' else str(r.get(h))
                              for h in header) + '\n')


def write_fasta(path, records, width=0):
    with open(path, 'w') as f:
        for name, seq in records:
            f.write('>' + name + '\n')
            if width:
                for i in range(0, len(seq), width):
                    f.write(seq[i:i + width] + '\n')
            else:
                f.write(seq + '\n')


def fmt_float(x, nd=4):
    return None if x is None else round(float(x), nd)


# ----------------------------------------------------------------------------
# Edit distance: Myers (1999) / Hyyro bit-vector algorithm on Python ints.
# Global unit-cost Levenshtein distance, exact.  ~0.05 s for 10 kb x 10 kb,
# ~1 s for 50 kb x 50 kb.

def _peq(a):
    return {c: int(''.join('1' if x == c else '0' for x in reversed(a)), 2) for c in set(a)}


def _load_fastedit():
    """ctypes handle on fastedit.c, compiled into the cache on first use; None if unavailable."""
    if os.environ.get('VNTR_NO_C'):
        return None
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fastedit.c')
    try:
        import ctypes
        import hashlib
        with open(src, 'rb') as f:
            key = hashlib.sha1(f.read()).hexdigest()[:12]
        lib = os.path.join(CACHE_DIR, 'libfastedit.%s.%s' % (key, 'dylib' if sys.platform == 'darwin' else 'so'))
        if not os.path.exists(lib):
            os.makedirs(CACHE_DIR, exist_ok=True)
            tmp = lib + '.tmp%d' % os.getpid()
            subprocess.run([config.CC, '-O3', '-shared', '-fPIC', '-o', tmp, src], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=config.tool_env())
            os.replace(tmp, lib)
        h = ctypes.CDLL(lib)
        h.edit_distance.restype = ctypes.c_long
        h.edit_distance.argtypes = [ctypes.c_char_p, ctypes.c_long, ctypes.c_char_p, ctypes.c_long]
        return h.edit_distance
    except (OSError, subprocess.CalledProcessError, AttributeError):
        return None


_FAST = None
_FAST_TRIED = False


def edit_distance(a, b):
    """Exact global unit-cost Levenshtein distance (C when available, else Python)."""
    global _FAST, _FAST_TRIED
    if not _FAST_TRIED:
        _FAST_TRIED = True
        _FAST = _load_fastedit()
    if _FAST is not None and a and b:
        if len(a) > len(b):
            a, b = b, a
        ab, bb = a.encode('latin-1', 'replace'), b.encode('latin-1', 'replace')
        d = _FAST(ab, len(ab), bb, len(bb))
        if d >= 0:
            return d
    return edit_distance_py(a, b)


def edit_distance_py(a, b):
    if not a:
        return len(b)
    if not b:
        return len(a)
    if len(a) > len(b):
        a, b = b, a
    m = len(a)
    mask = (1 << m) - 1
    high = 1 << (m - 1)
    peq = _peq(a)
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


def align_cigar(ref, hap, max_cells=4e8):
    """Global unit-cost alignment of hap against ref.  Returns (distance, cigar)
    with ops = X I D (I consumes hap only, D consumes ref only).  cigar is None
    when len(ref)*len(hap) > max_cells (the distance is still exact)."""
    m, n = len(ref), len(hap)
    if m == 0:
        return n, ('%dI' % n if n else '')
    if n == 0:
        return m, '%dD' % m
    if m * n > max_cells:
        return edit_distance(ref, hap), None
    mask = (1 << m) - 1
    high = 1 << (m - 1)
    peq = _peq(ref)
    vp, vn, score = mask, 0, m
    vps, vns = [vp], [vn]
    for c in hap:
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
        vps.append(vp)
        vns.append(vn)

    def cell(i, j):  # D[i][j]
        mi = (1 << i) - 1
        return j + (vps[j] & mi).bit_count() - (vns[j] & mi).bit_count()

    i, j, d = m, n, score
    ops = []
    while i > 0 and j > 0:
        dd = cell(i - 1, j - 1)
        same = ref[i - 1] == hap[j - 1]
        if dd + (0 if same else 1) == d:
            ops.append('=' if same else 'X')
            i -= 1
            j -= 1
            d = dd
            continue
        v = ((vps[j] >> (i - 1)) & 1) - ((vns[j] >> (i - 1)) & 1)  # D[i][j]-D[i-1][j]
        if v == 1:
            ops.append('D')
            i -= 1
            d -= 1
            continue
        ops.append('I')
        j -= 1
        d -= 1
    ops.extend('D' * i)
    ops.extend('I' * j)
    ops.reverse()
    return score, rle(ops)


def identity(dist, la, lb):
    m = max(la, lb)
    return 1.0 if m == 0 else 1.0 - dist / m


# ----------------------------------------------------------------------------
# GFA parsing.  A walk is a tuple of signed node ids (+id forward, -id reverse).

STEP_RE = re.compile(r'([<>])(\d+)')


def parse_steps(walk):
    return tuple(int(n) if o == '>' else -int(n) for o, n in STEP_RE.findall(walk))


def steps_to_walk(steps):
    return ''.join(('>' if s > 0 else '<') + str(abs(s)) for s in steps)


def reverse_steps(steps):
    return tuple(-s for s in reversed(steps))


class Gfa:
    def __init__(self, text):
        self.seqs = {}
        self.edges = []
        self.walks = []
        self.header = []
        for line in text.splitlines():
            if not line:
                continue
            f = line.split('\t')
            t = f[0]
            if t == 'S':
                self.seqs[int(f[1])] = f[2].upper()
            elif t == 'L':
                self.edges.append((int(f[1]), f[2], int(f[3]), f[4]))
            elif t == 'W':
                w = {'sample': f[1], 'hap': f[2], 'contig': f[3],
                     'start': int(f[4]) if f[4] != '*' else None,
                     'end': int(f[5]) if f[5] != '*' else None,
                     'steps': parse_steps(f[6]), 'weight': 1, 'cigar': None}
                for tag in f[7:]:
                    if tag.startswith('WT:i:'):
                        w['weight'] = int(tag[5:])
                    elif tag.startswith('CG:Z:'):
                        w['cigar'] = tag[5:]
                self.walks.append(w)
            elif t == 'H':
                self.header.append(line)

    def seqlen(self, node):
        return len(self.seqs[abs(node)])

    def walk_seq(self, steps):
        return ''.join(self.seqs[s] if s > 0 else revcomp(self.seqs[-s]) for s in steps)

    def walk_len(self, steps):
        return sum(len(self.seqs[abs(s)]) for s in steps)


def gbz_query(args, log):
    p = run([GBZ_BASE, 'query'] + args, log=log)
    return p.stdout, p.stderr


# ----------------------------------------------------------------------------
# VCF helpers

def parse_info(s):
    d = {}
    if s in ('.', ''):
        return d
    for kv in s.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def vcf_records(vcf, contig, a, b, log=NULL_LOG):
    """Records overlapping contig:a-b (1-based inclusive; htslib uses the REF span)."""
    if not os.path.exists(vcf):
        return None
    out = run([BCFTOOLS, 'view', '-H', '-r', '%s:%d-%d' % (contig, max(1, a), b), vcf], log=log).stdout
    recs = []
    for line in out.splitlines():
        f = line.split('\t')
        r = {'chrom': f[0], 'pos': int(f[1]), 'id': f[2], 'ref': f[3].upper(),
             'alts': [x.upper() for x in f[4].split(',')], 'qual': f[5], 'filter': f[6],
             'info': parse_info(f[7]), 'line': line}
        r['end'] = r['pos'] + len(r['ref']) - 1
        if 'END' in r['info'] and any(x.startswith('<') for x in r['alts']):
            try:
                r['end'] = int(r['info']['END'])
            except ValueError:
                pass
        fmt = f[8].split(':') if len(f) > 8 else []
        smp = f[9].split(':') if len(f) > 9 else []
        r['fmt'] = dict(zip(fmt, smp))
        gt = r['fmt'].get('GT', '.')
        r['gt'] = gt
        r['phased'] = '|' in gt
        r['gt_alleles'] = [None if x == '.' else int(x) for x in re.split(r'[|/]', gt)] if gt else []
        recs.append(r)
    return recs


def vcf_header(vcf, log=NULL_LOG):
    return run([BCFTOOLS, 'view', '-h', vcf], log=log).stdout


def normalise_records(recs, header_vcf, ref_fa, workdir, tag, log=NULL_LOG):
    """Left-align every ALT allele the way the scoring pipeline did (bcftools
    norm -m-any -f ref).  Alleles are split here, one line per ALT stamped
    ID=rec<i>_<allele>, so the mapping back survives any re-ordering by norm.
    Returns {record_index: {allele_index: (pos, REF, ALT)}}."""
    out = collections.defaultdict(dict)
    if not recs:
        return out
    hdr = [h for h in vcf_header(header_vcf, log).rstrip('\n').split('\n')]
    hdr[-1] = '\t'.join(hdr[-1].split('\t')[:8])  # drop FORMAT/sample columns
    body = []
    for i, r in enumerate(recs):
        f = r['line'].split('\t')
        for ai, alt in enumerate(r['alts'], start=1):
            if alt == '*' or alt.startswith('<'):
                continue
            body.append('\t'.join([f[0], f[1], 'rec%d_%d' % (i, ai), f[3], alt, '.', 'PASS', '.']))
    if not body:
        return out
    tmp = os.path.join(workdir, tag + '.pre_norm.vcf')
    with open(tmp, 'w') as fh:
        fh.write('\n'.join(hdr) + '\n' + '\n'.join(body) + '\n')
    p = run([BCFTOOLS, 'norm', '-m-any', '-f', ref_fa, tmp], log=log, check=False)
    if p.returncode != 0:
        log.info('WARNING: bcftools norm failed for %s; joins use raw alleles' % tag)
        for i, r in enumerate(recs):
            for ai, a in enumerate(r['alts'], start=1):
                out[i][ai] = (r['pos'], r['ref'], a)
        return out
    for line in p.stdout.splitlines():
        if line.startswith('#'):
            continue
        f = line.split('\t')
        m = re.match(r'rec(\d+)_(\d+)$', f[2])
        if m:
            out[int(m.group(1))][int(m.group(2))] = (int(f[1]), f[3].upper(), f[4].upper())
    os.remove(tmp)
    return out


def status_table(vcf_dir, names, contig, a, b, log=NULL_LOG):
    """{(pos, REF, ALT): (status, info)} from truvari/aardvark output files."""
    tab = {}
    for fname, status in names:
        path = os.path.join(vcf_dir, fname)
        recs = vcf_records(path, contig, a, b, log)
        if recs is None:
            continue
        for r in recs:
            for alt in r['alts']:
                st = status
                if status == 'BD':  # aardvark: per-record FORMAT/BD
                    st = r['fmt'].get('BD', '?')
                tab[(r['pos'], r['ref'], alt)] = (st, r['info'])
    return tab


def load_bed(path, contig):
    iv = []
    if not os.path.exists(path):
        return iv
    with open(path) as f:
        for line in f:
            if line.startswith(contig + '\t'):
                x = line.split('\t')
                iv.append((int(x[1]), int(x[2])))
    iv.sort()
    return iv


def bed_coverage(iv, a0, b0):
    """Fraction of 0-based half-open [a0,b0) covered by sorted intervals iv."""
    if b0 <= a0:
        return 0.0
    cov = 0
    for s, e in iv:
        if e <= a0:
            continue
        if s >= b0:
            break
        cov += min(e, b0) - max(s, a0)
    return cov / (b0 - a0)


def rmsk_overlaps(rmsk, contig, a0, b0, log=NULL_LOG):
    """RepeatMasker elements overlapping [a0,b0), from a per-contig cache."""
    if not os.path.exists(rmsk):
        return []
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache = os.path.join(CACHE_DIR, 'rmsk.%s.bed' % contig)
    if not os.path.exists(cache):
        tmp = cache + '.tmp%d' % os.getpid()
        with open(rmsk) as f, open(tmp, 'w') as out:
            for line in f:
                if line.startswith(contig + '\t'):
                    out.write(line)
        os.replace(tmp, cache)
    res = []
    with open(cache) as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            s, e = int(x[1]), int(x[2])
            if s < b0 and e > a0:
                name, cls = x[3], x[6] if len(x) > 6 else ''
                m = re.match(r'^\(([ACGTN]+)\)n$', name)
                res.append({'start0': s, 'end': e, 'name': name, 'class': cls,
                            'family': x[7] if len(x) > 7 else '',
                            'period': len(m.group(1)) if m else None})
    return res


# ----------------------------------------------------------------------------
# Reference probe and boundary choice

def probe_reference(db, contig, lo0, hi0, log):
    """--interval lo0..hi0 --context 0 --distinct.  Returns (ref_nodes, cov,
    maxvisit, H, gfa) where ref_nodes = [(signed_step, start0, end0)] along the
    CHM13 walk and cov[node] = number of haplotype visits on that node."""
    text, _ = gbz_query([db, '--sample', 'CHM13', '--contig', contig,
                         '--interval', '%d..%d' % (lo0, hi0), '--context', '0', '--distinct'], log)
    g = Gfa(text)
    ref = [w for w in g.walks if w['sample'] == 'CHM13']
    if not ref:
        raise ToolError('reference walk not found in probe %s:%d..%d' % (contig, lo0, hi0))
    ref = ref[0]
    pos = ref['start']
    nodes = []
    for s in ref['steps']:
        L = g.seqlen(s)
        nodes.append((s, pos, pos + L))
        pos += L
    cov = collections.Counter()
    maxvisit = collections.Counter()
    for w in g.walks:
        c = collections.Counter(abs(s) for s in w['steps'])
        for n, k in c.items():
            cov[n] += k * w['weight']
            if k > maxvisit[n]:
                maxvisit[n] = k
    H = max(cov.values()) if cov else 0
    return nodes, cov, maxvisit, H, g


def truth_spans(paths, contig, a1, b1, log):
    """1-based closed spans of truth records (either file) that change either
    haplotype, for keeping boundaries clear of truth variants."""
    spans = []
    for key in ('stvar', 'smvar'):
        recs = vcf_records(paths[key], contig, a1, b1, log) or []
        for r in recs:
            alleles = [x for x in r['gt_alleles'] if x]
            alts = [r['alts'][x - 1] for x in alleles if x - 1 < len(r['alts'])]
            if not any(a != '*' for a in alts):
                continue
            end = r['end']
            if any(len(a) > len(r['ref']) for a in alts):
                end += 1  # an insertion after the last REF base touches the next base
            spans.append((r['pos'], end))
    spans.sort()
    return spans


def span_hits(spans, s1, e1):
    for a, b in spans:
        if a > e1:
            break
        if b >= s1:
            return True
    return False


def choose_boundaries(db, paths, contig, contig_len, a0, b0, log, window=20000, max_window=320000,
                      avoid_truth=True, skip=(set(), set())):
    """Pick L (last good node ending <= a0) and R (first good node starting >= b0).
    Good: haplotype coverage == max coverage in the probe, each haplotype
    subpath visits it at most once, reference visits it once, and (optionally)
    no truth record touches it.  Returns dict or raises ToolError."""
    w = window
    while True:
        lo0 = max(0, a0 - w)
        hi0 = min(contig_len, b0 + w)
        nodes, cov, maxvisit, H, g = probe_reference(db, contig, lo0, hi0, log)
        refcount = collections.Counter(abs(s) for s, _, _ in nodes)
        spans = truth_spans(paths, contig, lo0 + 1, hi0, log) if avoid_truth else []
        good = []
        for idx, (s, st, en) in enumerate(nodes):
            n = abs(s)
            if cov[n] != H or maxvisit[n] > 1 or refcount[n] != 1:
                continue
            if spans and span_hits(spans, st + 1, en):
                continue
            good.append(idx)
        Lc = [i for i in good if nodes[i][2] <= a0 and abs(nodes[i][0]) not in skip[0]]
        Rc = [i for i in good if nodes[i][1] >= b0 and abs(nodes[i][0]) not in skip[1]]
        edge = []
        # at a contig end there may be no anchor outside the interval: use the end node
        if not Lc and lo0 == 0 and nodes:
            Lc = [0]
            edge.append('L at contig start (coverage %d of %d)' % (cov[abs(nodes[0][0])], H))
        if not Rc and hi0 == contig_len and nodes:
            Rc = [len(nodes) - 1]
            edge.append('R at contig end (coverage %d of %d)' % (cov[abs(nodes[-1][0])], H))
        if Lc and Rc:
            li, ri = Lc[-1], Rc[0]
            return {'L': nodes[li], 'R': nodes[ri], 'ref_nodes': nodes[li:ri + 1], 'H': H,
                    'probe': (lo0, hi0), 'probe_gfa': g, 'cov': cov, 'edge': edge}
        if (lo0 == 0 and hi0 == contig_len) or w >= max_window:
            raise ToolError('no boundary nodes found within %d bp (L %s, R %s)' %
                            (w, 'ok' if Lc else 'missing', 'ok' if Rc else 'missing'))
        w *= 2
        log.info('  widening boundary search window to %d bp' % w)


def parse_anchor(s):
    """'123+', '123-', '>123' or '<123' -> signed node id."""
    s = s.strip()
    m = re.match(r'^([<>])(\d+)$', s) or re.match(r'^(\d+)([+-])$', s)
    if not m:
        raise ToolError('cannot parse anchor %r (use 123+ or >123)' % s)
    a, b = m.groups()
    if a in '<>':
        return int(b) if a == '>' else -int(b)
    return int(a) if b == '+' else -int(a)


def fixed_boundaries(db, contig, contig_len, a0, b0, anchors, log, window=20000, max_window=320000):
    """Use the given anchor nodes (signed ids) as L and R.  The probe window around the
    interval is widened until both lie on the CHM13 walk.  Returns the same dict as
    choose_boundaries, with 'notes' for anything that breaks the anchor rule."""
    Lid, Rid = abs(anchors[0]), abs(anchors[1])
    w = window
    while True:
        lo0 = max(0, a0 - w)
        hi0 = min(contig_len, b0 + w)
        nodes, cov, maxvisit, H, g = probe_reference(db, contig, lo0, hi0, log)
        idx = {}
        for i, (s, _, _) in enumerate(nodes):
            idx.setdefault(abs(s), []).append(i)
        if Lid in idx and Rid in idx:
            li, ri = idx[Lid][0], idx[Rid][-1]
            if ri <= li:
                raise ToolError('anchor %s is not left of %s on CHM13' % (anchors[0], anchors[1]))
            notes = []
            for name, want, i in (('L', anchors[0], li), ('R', anchors[1], ri)):
                s = nodes[i][0]
                if (s > 0) != (want > 0):
                    notes.append('anchor %s: CHM13 walks it as %s' % (name, steps_to_walk([s])))
                if cov[abs(s)] != H:
                    notes.append('anchor %s: haplotype coverage %d of %d' % (name, cov[abs(s)], H))
                if maxvisit[abs(s)] > 1 or len(idx[abs(s)]) > 1:
                    notes.append('anchor %s: visited more than once' % name)
            return {'L': nodes[li], 'R': nodes[ri], 'ref_nodes': nodes[li:ri + 1], 'H': H,
                    'probe': (lo0, hi0), 'probe_gfa': g, 'cov': cov, 'edge': [], 'notes': notes}
        if (lo0 == 0 and hi0 == contig_len) or w >= max_window:
            raise ToolError('anchors %s / %s not on the CHM13 walk within %d bp of the interval' %
                            (steps_to_walk([anchors[0]]), steps_to_walk([anchors[1]]), w))
        w *= 2


# ----------------------------------------------------------------------------
# Walk analysis

class Region:
    """Everything known about the extracted subgraph."""

    def __init__(self, gfa, ref_steps, ref_start0, Lstep, Rstep):
        self.g = gfa
        self.ref_steps = tuple(ref_steps)
        self.ref_start0 = ref_start0
        self.L = Lstep
        self.R = Rstep
        self.ref_orient = {}
        for s in ref_steps:
            self.ref_orient.setdefault(abs(s), 1 if s > 0 else -1)
        self.ref_nodes = set(abs(s) for s in ref_steps)
        self.nodes = set(gfa.seqs)

    def orient_score(self, steps):
        sc = 0
        for s in steps:
            n = abs(s)
            want = self.ref_orient.get(n, 1)
            sc += len(self.g.seqs[n]) * (1 if (s > 0) == (want > 0) else -1)
        return sc

    def normalise(self, steps):
        """Return (steps in reference orientation, flipped?)."""
        sc = self.orient_score(steps)
        rev = reverse_steps(steps)
        if sc < 0 or (sc == 0 and rev < steps):
            return rev, True
        return steps, False

    def boundary_label(self, s):
        if s == self.L:
            return 'L'
        if s == -self.L:
            return 'L-'
        if s == self.R:
            return 'R'
        if s == -self.R:
            return 'R-'
        return 'inner'

    def classify(self, steps):
        a, b = self.boundary_label(steps[0]), self.boundary_label(steps[-1])
        if a == 'L' and b == 'R':
            return 'spanning'
        if a == 'L':
            return 'enters_L'
        if b == 'R':
            return 'exits_R'
        return 'internal'

    def walk_metrics(self, steps):
        g = self.g
        c = collections.Counter(abs(s) for s in steps)
        nonref = [n for n in c if n not in self.ref_nodes]
        rev_steps = sum(1 for s in steps if (s > 0) != (self.ref_orient.get(abs(s), 1) > 0))
        return {
            'length_bp': g.walk_len(steps),
            'n_steps': len(steps),
            'n_distinct_nodes': len(c),
            'n_nonref_nodes': len(nonref),
            'nonref_bp': sum(len(g.seqs[abs(s)]) for s in steps if abs(s) not in self.ref_nodes),
            'n_reverse_steps': rev_steps,
            'max_node_visits': max(c.values()) if c else 0,
            'n_nodes_revisited': sum(1 for v in c.values() if v > 1),
            'starts_at': self.boundary_label(steps[0]),
            'ends_at': self.boundary_label(steps[-1]),
            'first_node': steps_to_walk(steps[:1]),
            'last_node': steps_to_walk(steps[-1:]),
        }


def graph_has_directed_cycle(nodes, edges):
    """Cycle in the bidirected graph viewed as a directed graph on handles."""
    succ = collections.defaultdict(list)
    for a, ao, b, bo in edges:
        ha, hb = (a, ao == '-'), (b, bo == '-')
        succ[ha].append(hb)
        succ[(b, bo != '-')].append((a, ao != '-'))
    color = {}
    for n in nodes:
        for start in ((n, False), (n, True)):
            if start in color:
                continue
            stack = [(start, iter(succ[start]))]
            color[start] = 1
            while stack:
                v, it = stack[-1]
                nxt = next(it, None)
                if nxt is None:
                    color[v] = 2
                    stack.pop()
                    continue
                cs = color.get(nxt, 0)
                if cs == 1:
                    return True
                if cs == 0:
                    color[nxt] = 1
                    stack.append((nxt, iter(succ[nxt])))
    return False


# ----------------------------------------------------------------------------
# Named haplotype runs from `vg paths -A` (opt-in: ~40 s on chr20)

_STEP_AT = re.compile(r'([<>])(\d+)')


def runs_in_walk(walk, node_set, tok_re):
    """Maximal runs of consecutive steps of `walk` (a GAF path string) whose
    nodes are all in node_set.  Runs can only begin at the walk start or at a
    boundary token matched by tok_re (entry into a between-subgraph must pass a
    boundary node).  Returns [(char_offset, steps, reaches_walk_end)]."""
    walk = walk.rstrip('\n')
    starts = []
    m0 = _STEP_AT.match(walk)
    if m0 and int(m0.group(2)) in node_set:
        starts.append(0)
    for m in tok_re.finditer(walk):
        a = m.start()
        if a == 0:
            continue
        j = a - 1
        while j >= 0 and walk[j].isdigit():
            j -= 1
        prev = int(walk[j + 1:a]) if j + 1 < a else None
        if prev is None or prev not in node_set:
            starts.append(a)
    out = []
    for a in starts:
        steps = []
        pos = a
        while True:
            m = _STEP_AT.match(walk, pos)
            if not m:
                break
            n = int(m.group(2))
            if n not in node_set:
                break
            steps.append(n if m.group(1) == '>' else -n)
            pos = m.end()
        if steps:
            out.append((a, tuple(steps), pos >= len(walk)))
    return out


def boundary_token_re(L, R):
    return re.compile(r'[<>](%d|%d)(?![0-9])' % (abs(L), abs(R)))


def named_runs(contig_gbz, node_set, L, R, log):
    """Emulate gbz-base subpath extraction on every named path of the contig
    GBZ (streamed from `vg paths -A`).  Returns ([(path_name, offset, steps,
    reaches_path_end)], n_paths)."""
    tok = boundary_token_re(L, R)
    cmd = [VG, 'paths', '-x', contig_gbz, '-A']
    log.write('$ ' + ' '.join(cmd))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                         env=config.tool_env())
    runs = []
    n_paths = 0
    for line in p.stdout:
        if line.startswith('@'):
            continue
        f = line.split('\t', 6)
        if len(f) < 6:
            continue
        n_paths += 1
        for off, steps, at_end in runs_in_walk(f[5], node_set, tok):
            runs.append((f[0], off, steps, at_end))
    p.wait()
    err = p.stderr.read()
    if err:
        log.write('  stderr: ' + err.strip()[:2000])
    if p.returncode != 0:
        raise ToolError('vg paths -A failed: %s' % err[:500])
    return runs, n_paths


def hap_key(name):
    parts = name.split('#')
    if len(parts) >= 2:
        return parts[0] + '#' + parts[1]
    return name


# ----------------------------------------------------------------------------
# Truth haplotypes

def choose_truth(paths, contig, A1, B1, log):
    """Select, per haplotype, the truth records to apply over [A1,B1] (1-based).
    stvar first; smvar only where it does not duplicate/overlap an applied stvar
    record on that haplotype.  Returns (rows, applied{1:[..],2:[..]})."""
    rows = []
    applied = {1: [], 2: []}
    occupied = {1: [], 2: []}  # (start, end) of applied records per hap

    def overlaps(h, s, e):
        for a, b, src_ in occupied[h]:
            if a <= e and s <= b:
                return src_
        return None

    for src in ('stvar', 'smvar'):
        recs = vcf_records(paths[src], contig, A1, B1, log) or []
        recs.sort(key=lambda r: (r['pos'], r['end']))
        for r in recs:
            row = {'source': src, 'chrom': r['chrom'], 'pos': r['pos'], 'end': r['end'],
                   'ref_len': len(r['ref']), 'alt_lens': ','.join(str(len(a)) for a in r['alts']),
                   'svtype': r['info'].get('SVTYPE'), 'svlen': r['info'].get('SVLEN'),
                   'gt': r['gt'], 'TRF': 1 if r['info'].get('TRF') else 0,
                   'TRFperiod': r['info'].get('TRFperiod'), 'TRFrepeat': r['info'].get('TRFrepeat'),
                   'TRFstart': r['info'].get('TRFstart'), 'TRFend': r['info'].get('TRFend'),
                   'TRFcopies': r['info'].get('TRFcopies'), 'TRFdiff': r['info'].get('TRFdiff'),
                   'LCR': r['info'].get('LCR'), 'REMAP': r['info'].get('REMAP'),
                   'RM_clsfam': r['info'].get('RM_clsfam'),
                   'within_span': 1 if (r['pos'] >= A1 and r['end'] <= B1) else 0,
                   '_rec': r}
            real = [a for a in r['alts'] if a != '*' and not a.startswith('<')]
            row['max_allele_len_diff'] = max([abs(len(a) - len(r['ref'])) for a in real] or [0])
            row['is_sv50'] = 1 if row['max_allele_len_diff'] >= 50 else 0
            for h in (1, 2):
                col = 'applied_h%d' % h
                al = r['gt_alleles'][h - 1] if len(r['gt_alleles']) >= h else None
                if al is None:
                    row[col] = 'no:missing_gt'
                    continue
                if al == 0:
                    row[col] = 'ref'
                    continue
                alt = r['alts'][al - 1]
                if alt == '*':
                    row[col] = 'no:star_allele'
                    continue
                if alt.startswith('<') or '[' in alt or ']' in alt:
                    row[col] = 'no:symbolic'
                    continue
                if not row['within_span']:
                    row[col] = 'no:crosses_span_edge'
                    continue
                key = (r['pos'], r['ref'], alt)
                if src == 'smvar' and any((x['pos'], x['ref'], x['alt']) == key for x in applied[h]):
                    row[col] = 'no:duplicate_of_stvar'
                    continue
                hit = overlaps(h, r['pos'], r['end'])
                if hit:
                    row[col] = 'no:overlaps_applied_' + hit
                    continue
                applied[h].append({'pos': r['pos'], 'ref': r['ref'], 'alt': alt, 'src': src})
                occupied[h].append((r['pos'], r['end'], src))
                row[col] = 'yes'
            rows.append(row)
    return rows, applied


def apply_python(refseq, A1, recs):
    s = refseq
    for r in sorted(recs, key=lambda x: -x['pos']):
        i = r['pos'] - A1
        assert s[i:i + len(r['ref'])].upper() == r['ref'], ('REF mismatch', r)
        s = s[:i] + r['alt'] + s[i + len(r['ref']):]
    return s


def build_truth(paths, contig, A1, B1, refseq, outdir, workdir, log):
    rows, applied = choose_truth(paths, contig, A1, B1, log)
    # One resolved VCF, GT per haplotype, then bcftools consensus -H 1 / -H 2.
    merged = collections.OrderedDict()
    for h in (1, 2):
        for r in applied[h]:
            k = (r['pos'], r['ref'], r['alt'])
            merged.setdefault(k, set()).add(h)
    keys = sorted(merged)
    vcf = os.path.join(workdir, 'truth.resolved.vcf')
    ctg_len = None
    fai = paths['ref_fa'] + '.fai'
    if os.path.exists(fai):
        for line in open(fai):
            x = line.split('\t')
            if x[0] == contig:
                ctg_len = int(x[1])
    with open(vcf, 'w') as f:
        f.write('##fileformat=VCFv4.2\n##contig=<ID=%s%s>\n' %
                (contig, (',length=%d' % ctg_len) if ctg_len else ''))
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n')
        for k in keys:
            hs = merged[k]
            gt = '%d|%d' % (1 if 1 in hs else 0, 1 if 2 in hs else 0)
            f.write('%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s\n' % (contig, k[0], k[1], k[2], gt))
    region_fa = os.path.join(workdir, 'ref.region.fa')
    with open(region_fa, 'w') as f:
        f.write('>%s:%d-%d\n%s\n' % (contig, A1, B1, refseq))
    seqs = {}
    notes = []
    if keys:
        run([BGZIP, '-f', vcf], log=log)
        run([BCFTOOLS, 'index', '-f', vcf + '.gz'], log=log)
    for h in (1, 2):
        py = apply_python(refseq, A1, applied[h])
        if keys:
            out = run([BCFTOOLS, 'consensus', '-f', region_fa, '-H', str(h), vcf + '.gz'], log=log).stdout
            bc = ''.join(out.split('\n', 1)[1].split()).upper()
        else:
            bc = refseq
        if bc != py:
            notes.append('hap%d: bcftools consensus and python application differ (%d vs %d bp); '
                         'bcftools used' % (h, len(bc), len(py)))
        seqs[h] = bc
        write_fasta(os.path.join(outdir, 'truth.hap%d.fa' % h),
                    [('truth_h%d %s:%d-%d applied=%d' % (h, contig, A1, B1, len(applied[h])), bc)], 80)
    return rows, applied, seqs, notes


# ----------------------------------------------------------------------------
# Reads

def summarise_reads(gaf_path, region, span_len, truth_lens, depth, log):
    mapq = []
    strands = collections.Counter()
    strands_mq = collections.Counter()
    spandiff = []
    inside = crosses = 0
    nonref_touch = 0
    qlens = []
    idents = []
    names = set()
    n = 0
    with open(gaf_path) as f:
        for line in f:
            if line.startswith('@') or not line.strip():
                continue
            x = line.rstrip('\n').split('\t')
            if len(x) < 12:
                continue
            n += 1
            names.add(x[0])
            try:
                qlens.append(int(x[1]))
                mq = int(x[11])
            except ValueError:
                mq = 255
            mapq.append(mq)
            steps = parse_steps(x[5])
            ids = set(abs(s) for s in steps)
            if ids and ids <= region.nodes:
                inside += 1
            else:
                crosses += 1
            if any(i in region.nodes and i not in region.ref_nodes for i in ids):
                nonref_touch += 1
            known = [s for s in steps if abs(s) in region.g.seqs]
            sc = region.orient_score(known) if known else 0
            strand = '+' if sc >= 0 else '-'
            strands[strand] += 1
            strands_mq[('mapq<5' if mq < 5 else 'mapq>=5') + ' ' + strand] += 1
            try:
                if int(x[10]) > 0:
                    idents.append(int(x[9]) / int(x[10]))
                qspan = int(x[3]) - int(x[2])
                pspan = int(x[8]) - int(x[7])
                spandiff.append(qspan - pspan)
            except ValueError:
                pass
    bins = [(0, 0), (1, 4), (5, 9), (10, 19), (20, 29), (30, 59), (60, 60), (61, 255)]
    hist = collections.OrderedDict()
    for a, b in bins:
        hist['%d' % a if a == b else '%d-%d' % (a, b)] = sum(1 for q in mapq if a <= q <= b)
    rl = statistics.median(qlens) if qlens else 150
    exp_ref = depth * (span_len + rl - 1) / rl
    exp_truth = None
    if truth_lens:
        exp_truth = sum((depth / 2.0) * (lh + rl - 1) / rl for lh in truth_lens)
    return {
        'n_alignments': n,
        'n_distinct_read_names': len(names),
        'mapq_histogram': hist,
        'frac_mapq_lt5': fmt_float(sum(1 for q in mapq if q < 5) / n) if n else None,
        'frac_mapq0': fmt_float(sum(1 for q in mapq if q == 0) / n) if n else None,
        'median_mapq': statistics.median(mapq) if mapq else None,
        'strand_vs_reference': dict(strands),
        'median_read_length': rl,
        'strand_by_mapq': dict(sorted(strands_mq.items())),
        'median_block_identity': fmt_float(statistics.median(idents)) if idents else None,
        'frac_block_identity_lt_0.9': fmt_float(sum(1 for i in idents if i < 0.9) / len(idents)) if idents else None,
        'frac_block_identity_lt_0.5': fmt_float(sum(1 for i in idents if i < 0.5) / len(idents)) if idents else None,
        'frac_query_minus_path_span_ge20': fmt_float(sum(1 for d in spandiff if d >= 20) / len(spandiff))
        if spandiff else None,
        'frac_path_minus_query_span_ge20': fmt_float(sum(1 for d in spandiff if d <= -20) / len(spandiff))
        if spandiff else None,
        'n_path_inside_subgraph': inside,
        'n_path_leaves_subgraph': crosses,
        'frac_path_inside_subgraph': fmt_float(inside / n) if n else None,
        'n_touching_nonref_nodes': nonref_touch,
        'depth_assumed': depth,
        'expected_formula': 'expected = depth * (L + r - 1) / r, r = median read length; '
                            'from reference: L = CHM13 span; from truth: sum over the two HG002 '
                            'haplotypes of (depth/2) * (L_h + r - 1) / r',
        'expected_from_reference_span': fmt_float(exp_ref, 1),
        'expected_from_truth_haplotypes': fmt_float(exp_truth, 1),
        'observed_over_expected_truth': fmt_float(n / exp_truth, 3) if exp_truth else None,
        'observed_over_expected_reference': fmt_float(n / exp_ref, 3) if exp_ref else None,
    }


# ----------------------------------------------------------------------------
# Calls

CALL_COLS = ['caller', 'rec', 'chrom', 'pos', 'end', 'id', 'ref_len', 'allele', 'alt_len',
             'len_diff', 'called', 'gt', 'qual', 'filter', 'in_interval', 'norm_pos',
             'norm_ref_len', 'norm_alt_len', 'truvari', 'aardvark', 'PctSeqSimilarity',
             'PctSizeSimilarity', 'SizeDiff', 'StartDistance', 'MatchId', 'overlaps_span_edge']


def calls_rows(caller, vcf, truvari_dir, aardvark_dir, contig, A1, B1, I1, I2, paths, workdir,
               log, already_normalised=False):
    recs = vcf_records(vcf, contig, A1, B1, log)
    if recs is None:
        return [], 'missing %s' % vcf
    margin = 1000
    tv = status_table(truvari_dir, [('tp-comp.vcf.gz', 'TP'), ('fp.vcf.gz', 'FP')],
                      contig, A1 - margin, B1 + margin, log)
    ad = status_table(aardvark_dir, [('query.vcf.gz', 'BD')], contig, A1 - margin, B1 + margin, log)
    if already_normalised:
        norm = {i: {ai: (r['pos'], r['ref'], a) for ai, a in enumerate(r['alts'], start=1)}
                for i, r in enumerate(recs)}
    else:
        norm = normalise_records(recs, vcf, paths['ref_fa'], workdir, caller, log)
    rows = []
    for i, r in enumerate(recs):
        called = set(x for x in r['gt_alleles'] if x)
        for ai, alt in enumerate(r['alts'], start=1):
            nk = norm.get(i, {}).get(ai)
            raw = (r['pos'], r['ref'], alt)
            tvs = tv.get(nk) if nk else None
            if tvs is None:
                tvs = tv.get(raw)
            ads = ad.get(raw) or (ad.get(nk) if nk else None)
            info = tvs[1] if tvs else {}
            rows.append({
                'caller': caller, 'rec': i, 'chrom': r['chrom'], 'pos': r['pos'], 'end': r['end'],
                'id': r['id'], 'ref_len': len(r['ref']), 'allele': ai, 'alt_len': len(alt),
                'len_diff': len(alt) - len(r['ref']), 'called': 1 if ai in called else 0,
                'gt': r['gt'], 'qual': r['qual'], 'filter': r['filter'],
                'in_interval': 1 if (r['pos'] <= I2 and r['end'] >= I1) else 0,
                'norm_pos': nk[0] if nk else None, 'norm_ref_len': len(nk[1]) if nk else None,
                'norm_alt_len': len(nk[2]) if nk else None,
                'truvari': tvs[0] if tvs else '-', 'aardvark': ads[0] if ads else '-',
                'PctSeqSimilarity': info.get('PctSeqSimilarity'),
                'PctSizeSimilarity': info.get('PctSizeSimilarity'),
                'SizeDiff': info.get('SizeDiff'), 'StartDistance': info.get('StartDistance'),
                'MatchId': info.get('MatchId'),
                'overlaps_span_edge': 1 if (r['pos'] < A1 or r['end'] > B1) else 0,
            })
    return rows, None


def called_haplotypes(recs, A1, B1, refseq):
    """Apply a caller's genotypes over [A1,B1] to CHM13, one sequence per GT slot.
    Records crossing the span edges, '*'/symbolic alleles and records overlapping
    an already-applied record in the same slot are skipped and counted."""
    out = {'skipped_edge': 0, 'skipped_overlap': 0, 'skipped_other': 0, 'missing_allele': 0,
           'unphased_het': 0, 'phase_sets': set(), 'n_applied': [0, 0]}
    slots = [[], []]
    occ = [[], []]
    # longest first at equal POS, so a nested parent record wins over its children
    for r in sorted(recs, key=lambda r: (r['pos'], -r['end'])):
        gts = r['gt_alleles']
        if len(gts) == 1:
            gts = gts * 2
        if len(gts) < 2:
            continue
        nonref = [x for x in gts if x]
        if not nonref:
            continue
        if r['pos'] < A1 or r['end'] > B1:
            out['skipped_edge'] += 1
            continue
        if not r['phased'] and gts[0] != gts[1]:
            out['unphased_het'] += 1
        if 'PS' in r['fmt'] and r['fmt']['PS'] not in ('.', ''):
            out['phase_sets'].add(r['fmt']['PS'])
        for h in (0, 1):
            a = gts[h]
            if a is None:
                out['missing_allele'] += 1
                continue
            if a == 0:
                continue
            alt = r['alts'][a - 1]
            if alt == '*' or alt.startswith('<'):
                out['skipped_other'] += 1
                continue
            if any(s <= r['end'] and r['pos'] <= e for s, e in occ[h]):
                out['skipped_overlap'] += 1
                continue
            slots[h].append({'pos': r['pos'], 'ref': r['ref'], 'alt': alt})
            occ[h].append((r['pos'], r['end']))
            out['n_applied'][h] += 1
    seqs = []
    for h in (0, 1):
        try:
            seqs.append(apply_python(refseq, A1, slots[h]))
        except AssertionError:
            seqs.append(None)
    out['phase_sets'] = len(out['phase_sets'])
    return seqs, out


def pair_distance(c1, c2, t1, t2):
    """Best pairing of called slots to truth haplotypes by summed edit distance."""
    a = (edit_distance(c1, t1), edit_distance(c2, t2))
    b = (edit_distance(c1, t2), edit_distance(c2, t1))
    if sum(a) <= sum(b):
        return {'pairing': 'slot1-h1,slot2-h2', 'd_h1': a[0], 'd_h2': a[1], 'total': sum(a)}
    return {'pairing': 'slot1-h2,slot2-h1', 'd_h1': b[1], 'd_h2': b[0], 'total': sum(b)}


def truth_truvari_status(rows, paths, contig, A1, B1, workdir, log):
    """Add vg/PanGenie truvari base status (TP/FN) to stvar truth rows."""
    st = [row for row in rows if row['source'] == 'stvar']
    recs = [row['_rec'] for row in st]
    norm = normalise_records(recs, paths['stvar'], paths['ref_fa'], workdir, 'truth', log)
    margin = 1000
    tabs = {}
    for caller, d in (('vg', paths['vg_truvari']), ('pg', paths['pg_truvari'])):
        tabs[caller] = status_table(d, [('tp-base.vcf.gz', 'TP'), ('fn.vcf.gz', 'FN')],
                                    contig, A1 - margin, B1 + margin, log)
    for i, row in enumerate(st):
        keys = list(norm.get(i, {}).values())
        for caller in ('vg', 'pg'):
            sts = []
            for k in keys:
                if k in tabs[caller]:
                    sts.append(tabs[caller][k][0])
            row[caller + '_truvari'] = ','.join(sts) if sts else '-'
    for row in rows:
        if row['source'] != 'stvar':
            row.setdefault('vg_truvari', '-')
            row.setdefault('pg_truvari', '-')


# ----------------------------------------------------------------------------
# Main driver

HAP_COLS = ['hap_id', 'weight', 'is_reference', 'class', 'length_bp', 'len_minus_ref',
            'n_steps', 'n_distinct_nodes', 'n_nonref_nodes', 'nonref_bp', 'n_reverse_steps',
            'max_node_visits', 'n_nodes_revisited', 'flipped', 'starts_at', 'ends_at',
            'edit_to_ref', 'identity_to_ref', 'edit_to_truth_h1', 'identity_to_truth_h1',
            'edit_to_truth_h2', 'identity_to_truth_h2', 'first_node', 'last_node',
            'cigar_gbz', 'cigar']

SUBPATH_COLS = ['subpath_id', 'distinct_id', 'class', 'length_bp', 'n_steps',
                'max_node_visits', 'n_nodes_revisited', 'flipped', 'starts_at', 'ends_at',
                'haplotype', 'path_name']

NAMED_COLS = ['haplotype', 'path_fragments', 'n_subpaths', 'n_spanning', 'classes',
              'total_bp', 'spanning_bp', 'has_cycle', 'max_node_visits',
              'n_nodes_revisited', 'distinct_ids', 'fragment_ends_inside']


def contig_length(fa, contig):
    fai = fa + '.fai'
    if os.path.exists(fai):
        for line in open(fai):
            x = line.split('\t')
            if x[0] == contig:
                return int(x[1])
    raise ToolError('contig %s not in %s' % (contig, fai))


def fetch_ref(fa, contig, a1, b1, log):
    out = run([SAMTOOLS, 'faidx', fa, '%s:%d-%d' % (contig, a1, b1)], log=log).stdout
    return ''.join(out.split('\n', 1)[1].split()).upper()


def analyse(contig, start1, end1, outdir, pad=0, reads=True, names=False, mode='between',
            context=0, depth=30.0, msa=False, bandage=False, window=20000, max_span=250000,
            limit=300000, max_align_cells=4e8, verbose=True, anchors=None):
    """Analyse contig:start1-end1 (1-based inclusive) into outdir; returns the summary dict.
    anchors: optional (L, R) as signed node ids or strings ('123+', '>123') to use as the
    boundary nodes instead of choosing them (mode 'between' only)."""
    if end1 < start1:
        raise ToolError('end < start')
    os.makedirs(outdir, exist_ok=True)
    workdir = tempfile.mkdtemp(prefix='work.', dir=outdir)
    log = Log(os.path.join(outdir, 'run.log'), verbose)
    paths = data_paths(contig)
    for k in ('gbz_db', 'ref_fa'):
        if not os.path.exists(paths[k]):
            raise ToolError('missing %s' % paths[k])
    summary = collections.OrderedDict()
    summary['contig'] = contig
    summary['requested_interval_1based'] = [start1, end1]
    summary['pad'] = pad
    summary['mode'] = mode
    if anchors:
        anchors = tuple(parse_anchor(a) if isinstance(a, str) else int(a) for a in anchors)
        if mode != 'between':
            raise ToolError('--anchors needs --mode between')
        summary['fixed_anchors'] = [steps_to_walk([a]) for a in anchors]
    summary['tools'] = {'gbz-base': GBZ_BASE, 'bcftools': BCFTOOLS, 'samtools': SAMTOOLS, 'vg': VG}
    warnings = []
    clen = contig_length(paths['ref_fa'], contig)
    a0 = max(0, start1 - 1 - pad)          # padded interval, 0-based half-open
    b0 = min(clen, end1 + pad)
    summary['padded_interval_1based'] = [a0 + 1, b0]

    # ---- 1. subgraph extraction -------------------------------------------
    log.info('%s:%d-%d  mode=%s' % (contig, start1, end1, mode))
    if mode == 'between':
        skip = (set(), set())
        for attempt in range(4):
            if anchors:
                if attempt:
                    raise ToolError('the fixed anchors do not separate the graph')
                bd = fixed_boundaries(paths['gbz_db'], contig, clen, a0, b0, anchors, log, window=window)
            else:
                bd = choose_boundaries(paths['gbz_db'], paths, contig, clen, a0, b0, log, window=window,
                                       skip=skip)
            Ls, Lst, Len_ = bd['L']
            Rs, Rst, Ren = bd['R']
            span = Ren - Lst
            if span > max_span:
                raise ToolError('boundary span %d bp exceeds --max-span %d (L=%d R=%d); '
                                'try --mode interval' % (span, max_span, Lst, Ren))
            btw = '%d%s:%d%s' % (abs(Ls), '+' if Ls > 0 else '-', abs(Rs), '+' if Rs > 0 else '-')
            log.info('boundaries L=%s [%d,%d) R=%s [%d,%d)  span %d bp  H=%d' %
                     (steps_to_walk([Ls]), Lst, Len_, steps_to_walk([Rs]), Rst, Ren, span, bd['H']))
            try:
                text, err = gbz_query([paths['gbz_db'], '--between', btw, '--limit', str(limit),
                                       '--distinct'], log)
            except ToolError as e:
                if 'Found more than' in str(e):
                    log.info('  between query leaked past the boundaries; moving outward')
                    skip[0].add(abs(Ls))
                    skip[1].add(abs(Rs))
                    continue
                raise
            g = Gfa(text)
            ref_steps = [s for s, _, _ in bd['ref_nodes']]
            # every CHM13 node inside the subgraph must lie between L and R
            inside_ref = set(abs(s) for s in ref_steps)
            probe_ref = set(abs(s) for s, _, _ in probe_reference_nodes(bd))
            outside = [n for n in g.seqs if n in probe_ref and n not in inside_ref]
            if outside:
                log.info('  %d reference nodes outside [L,R] were pulled in; moving outward' % len(outside))
                skip[0].add(abs(Ls))
                skip[1].add(abs(Rs))
                continue
            break
        else:
            raise ToolError('could not find separating boundary nodes')
        region = Region(g, ref_steps, Lst, Ls, Rs)
        query_args = [paths['gbz_db'], '--between', btw, '--limit', str(limit)]
        for e in bd.get('edge', []):
            warnings.append('boundary fallback: ' + e)
        for e in bd.get('notes', []):
            warnings.append('fixed anchors: ' + e)
        summary['boundaries'] = {
            'L': {'node': steps_to_walk([Ls]), 'start0': Lst, 'end0': Len_},
            'R': {'node': steps_to_walk([Rs]), 'start0': Rst, 'end0': Ren},
            'haplotypes_at_boundaries_H': bd['H'],
            'probe_window_0based': list(bd['probe']),
            'rule': 'CHM13 nodes with haplotype coverage == max coverage in the probe window, '
                    'visited at most once by every haplotype subpath and once by CHM13, and not '
                    'touched by any truth record; L = last such node ending at or before the '
                    'padded start, R = first such node starting at or after the padded end',
        }
        if anchors:
            summary['boundaries']['rule'] = 'fixed by --anchors'
        with open(os.path.join(outdir, 'subgraph.query.gfa'), 'w') as f:
            f.write(text)
    else:
        text, err = gbz_query([paths['gbz_db'], '--sample', 'CHM13', '--contig', contig,
                               '--interval', '%d..%d' % (a0, b0), '--context', str(context),
                               '--distinct', '--cigar'], log)
        g = Gfa(text)
        ref = [w for w in g.walks if w['sample'] == 'CHM13']
        if not ref:
            raise ToolError('no reference walk in interval query')
        ref = ref[0]
        ref_steps = list(ref['steps'])
        Ls, Rs = ref_steps[0], ref_steps[-1]
        Lst = ref['start']
        Ren = ref['end']
        region = Region(g, ref_steps, Lst, Ls, Rs)
        query_args = [paths['gbz_db'], '--sample', 'CHM13', '--contig', contig,
                      '--interval', '%d..%d' % (a0, b0), '--context', str(context)]
        summary['boundaries'] = {'L': {'node': steps_to_walk([Ls]), 'start0': Lst},
                                 'R': {'node': steps_to_walk([Rs]), 'end0': Ren},
                                 'rule': 'interval mode: first/last CHM13 nodes of the query'}
        with open(os.path.join(outdir, 'subgraph.query.gfa'), 'w') as f:
            f.write(text)
        warnings.append('interval mode: haplotypes that leave the query subgraph are cut into '
                        'fragments (with --context 0 the subgraph is only the CHM13 nodes)')
    A1, B1 = Lst + 1, Ren
    span_len = B1 - A1 + 1
    summary['analysed_span_1based'] = [A1, B1]
    summary['analysed_span_bp'] = span_len
    refseq = fetch_ref(paths['ref_fa'], contig, A1, B1, log)
    ref_walk_seq = g.walk_seq(tuple(ref_steps))
    if ref_walk_seq != refseq:
        warnings.append('CHM13 walk sequence differs from the FASTA over the span '
                        '(%d vs %d bp)' % (len(ref_walk_seq), len(refseq)))

    # ---- 2. distinct haplotypes --------------------------------------------
    log.info('analysing %d nodes, %d distinct walks' % (len(g.seqs), len(g.walks)))
    ref_norm, _ = region.normalise(tuple(ref_steps))
    distinct = []
    for w in g.walks:
        steps, flipped = region.normalise(w['steps'])
        d = {'steps': steps, 'weight': w['weight'], 'flipped': 1 if flipped else 0,
             'cigar_gbz': w['cigar'] if not flipped else (w['cigar'] and 'reversed:' + w['cigar'])}
        d.update(region.walk_metrics(steps))
        d['class'] = region.classify(steps)
        d['is_reference'] = 1 if steps == ref_norm else 0
        distinct.append(d)
    # order: reference first, then spanning by weight, then fragments
    order = {'spanning': 0, 'enters_L': 1, 'exits_R': 2, 'internal': 3}
    distinct.sort(key=lambda d: (-d['is_reference'], order[d['class']], -d['weight'], -d['length_bp']))
    for i, d in enumerate(distinct):
        d['hap_id'] = 'h%d' % i
    ref_len = len(refseq)
    seq_of = {}
    for d in distinct:
        s = g.walk_seq(d['steps'])
        seq_of[d['hap_id']] = s
        if d['class'] == 'spanning':
            d['len_minus_ref'] = d['length_bp'] - ref_len
    n_ref_walk = sum(1 for d in distinct if d['is_reference'])
    if n_ref_walk != 1:
        warnings.append('reference walk matched %d distinct walks' % n_ref_walk)

    # ---- 3. truth ------------------------------------------------------------
    log.info('reconstructing truth haplotypes over %s:%d-%d' % (contig, A1, B1))
    truth_rows, applied, tseqs, tnotes = build_truth(paths, contig, A1, B1, refseq, outdir, workdir, log)
    warnings.extend(tnotes)
    try:
        truth_truvari_status(truth_rows, paths, contig, A1, B1, workdir, log)
    except ToolError as e:
        warnings.append('truth truvari join failed: %s' % str(e)[:200])
    sv_bed = load_bed(paths['stvar_bed'], contig)
    sm_bed = load_bed(paths['smvar_bed'], contig)
    for row in truth_rows:
        row['in_sv_benchmark'] = 1 if bed_coverage(sv_bed, row['pos'] - 1, row['end']) == 1.0 else 0
        row['in_requested_interval'] = 1 if (row['pos'] <= end1 and row['end'] >= start1) else 0
    tcols = ['source', 'chrom', 'pos', 'end', 'ref_len', 'alt_lens', 'max_allele_len_diff', 'is_sv50',
             'svtype', 'svlen', 'gt', 'applied_h1', 'applied_h2', 'within_span',
             'in_requested_interval', 'in_sv_benchmark', 'vg_truvari', 'pg_truvari', 'TRF',
             'TRFperiod', 'TRFrepeat', 'TRFstart', 'TRFend', 'TRFcopies', 'TRFdiff', 'LCR', 'REMAP',
             'RM_clsfam']
    write_tsv(os.path.join(outdir, 'truth_records.tsv'), tcols, truth_rows)
    th1, th2 = tseqs[1], tseqs[2]
    summary['truth'] = collections.OrderedDict([
        ('hap1_len', len(th1)), ('hap2_len', len(th2)), ('reference_len', ref_len),
        ('n_records_in_span', len(truth_rows)),
        ('n_applied_h1', len(applied[1])), ('n_applied_h2', len(applied[2])),
        ('n_sv50_stvar_records', sum(1 for r in truth_rows if r['is_sv50'] and r['source'] == 'stvar')),
        ('n_sv50_stvar_TRF_period_ge7', sum(1 for r in truth_rows if r['is_sv50'] and r['source'] == 'stvar'
                                            and r['TRFperiod'] and int(r['TRFperiod']) >= 7)),
        ('n_sv50_stvar_in_requested_interval', sum(1 for r in truth_rows if r['is_sv50'] and r['source'] == 'stvar'
                                                   and r['in_requested_interval'])),
        ('requested_interval_frac_in_sv_benchmark',
         fmt_float(bed_coverage(sv_bed, start1 - 1, end1))),
        ('span_frac_in_sv_benchmark', fmt_float(bed_coverage(sv_bed, A1 - 1, B1))),
        ('span_frac_in_smvar_benchmark', fmt_float(bed_coverage(sm_bed, A1 - 1, B1))),
        ('edit_h1_vs_ref', edit_distance(th1, refseq)), ('edit_h2_vs_ref', edit_distance(th2, refseq)),
        ('edit_h1_vs_h2', edit_distance(th1, th2)),
    ])
    if summary['truth']['requested_interval_frac_in_sv_benchmark'] < 1.0:
        warnings.append('requested interval is %.0f%% inside the SV benchmark BED; truth outside it '
                        'is unreliable and truvari does not score calls there' %
                        (100 * summary['truth']['requested_interval_frac_in_sv_benchmark']))

    # ---- 4. distances: every spanning distinct walk vs ref and truth -----------
    log.info('edit distances for %d spanning walks' % sum(1 for d in distinct if d['class'] == 'spanning'))
    for d in distinct:
        if d['class'] != 'spanning':
            continue
        s = seq_of[d['hap_id']]
        dist, cig = align_cigar(refseq, s, max_align_cells)
        d['edit_to_ref'] = dist
        d['identity_to_ref'] = fmt_float(identity(dist, len(s), ref_len))
        d['cigar'] = cig if cig is not None else 'NA:too_long'
        for h, t in ((1, th1), (2, th2)):
            e = edit_distance(s, t)
            d['edit_to_truth_h%d' % h] = e
            d['identity_to_truth_h%d' % h] = fmt_float(identity(e, len(s), len(t)))
    write_tsv(os.path.join(outdir, 'haplotypes.tsv'), HAP_COLS, distinct)
    write_fasta(os.path.join(outdir, 'haplotypes.fa'),
                [('%s weight=%d class=%s len=%d%s' % (d['hap_id'], d['weight'], d['class'],
                                                      d['length_bp'], ' reference=CHM13' if d['is_reference'] else ''),
                  seq_of[d['hap_id']]) for d in distinct])

    # closest panel haplotype per truth haplotype
    spanning = [d for d in distinct if d['class'] == 'spanning']
    closest = collections.OrderedDict()
    for h, t in ((1, th1), (2, th2)):
        ranked = sorted(spanning, key=lambda d: (d['edit_to_truth_h%d' % h], abs(d['length_bp'] - len(t))))
        entry = collections.OrderedDict()
        entry['truth_len'] = len(t)
        entry['method'] = 'exact global unit-cost edit distance (Myers bit-vector), all spanning distinct walks'
        if ranked:
            b = ranked[0]
            entry['best'] = {'hap_id': b['hap_id'], 'weight': b['weight'], 'is_reference': b['is_reference'],
                             'edit_distance': b['edit_to_truth_h%d' % h],
                             'length_diff_panel_minus_truth': b['length_bp'] - len(t),
                             'identity': b['identity_to_truth_h%d' % h]}
            nr = [d for d in ranked if not d['is_reference']]
            if nr:
                b2 = nr[0]
                entry['best_non_reference'] = {'hap_id': b2['hap_id'], 'weight': b2['weight'],
                                               'edit_distance': b2['edit_to_truth_h%d' % h],
                                               'length_diff_panel_minus_truth': b2['length_bp'] - len(t),
                                               'identity': b2['identity_to_truth_h%d' % h]}
            entry['reference_edit_distance'] = edit_distance(refseq, t)
            entry['n_exact_matches_weight'] = sum(d['weight'] for d in ranked if d['edit_to_truth_h%d' % h] == 0)
            entry['top5'] = [(d['hap_id'], d['weight'], d['edit_to_truth_h%d' % h], d['length_bp'] - len(t))
                             for d in ranked[:5]]
        else:
            entry['best'] = None
        closest['truth_h%d' % h] = entry
    with open(os.path.join(outdir, 'closest.json'), 'w') as f:
        json.dump(closest, f, indent=1)
    write_tsv(os.path.join(outdir, 'closest.tsv'),
              ['hap_id', 'weight', 'is_reference', 'length_bp', 'edit_to_ref', 'edit_to_truth_h1',
               'identity_to_truth_h1', 'edit_to_truth_h2', 'identity_to_truth_h2'], spanning)

    # ---- 5. individual (non-distinct) subpaths, plus reads -------------------
    log.info('individual haplotype subpaths%s' % (' + reads' if reads else ''))
    qargs = list(query_args)
    gaf_path = os.path.join(outdir, 'reads.gaf')
    if reads:
        if not os.path.exists(paths['gaf_db']):
            warnings.append('GAF-base missing; no reads')
            reads = False
        else:
            qargs += ['--gaf-base', paths['gaf_db'], '--gaf-output', gaf_path, '--alignments', 'overlapping']
    text2, err2 = gbz_query(qargs, log)
    with open(os.path.join(outdir, 'subpaths.gfa'), 'w') as f:
        f.write(text2)
    g2 = Gfa(text2)
    dist_index = {d['steps']: d['hap_id'] for d in distinct}
    subpaths = []
    for i, w in enumerate(g2.walks):
        steps, flipped = region.normalise(w['steps'])
        m = region.walk_metrics(steps)
        subpaths.append({'subpath_id': 's%d' % i, 'distinct_id': dist_index.get(steps, '?'),
                         'class': region.classify(steps), 'length_bp': m['length_bp'],
                         'n_steps': m['n_steps'], 'max_node_visits': m['max_node_visits'],
                         'n_nodes_revisited': m['n_nodes_revisited'], 'flipped': 1 if flipped else 0,
                         'starts_at': m['starts_at'], 'ends_at': m['ends_at'], '_steps': steps})
    n_unmatched = sum(1 for s in subpaths if s['distinct_id'] == '?')
    if n_unmatched:
        warnings.append('%d individual subpaths match no distinct walk' % n_unmatched)

    # ---- 6. optional: names from vg paths -A ---------------------------------
    named_summary = None
    if names:
        if mode != 'between':
            warnings.append('--names needs --mode between; skipped')
        elif not os.path.exists(paths['contig_gbz']):
            warnings.append('--names: %s missing; skipped' % paths['contig_gbz'])
        else:
            log.info('naming haplotypes with vg paths -A on %s (slow)' % paths['contig_gbz'])
            runs, n_paths = named_runs(paths['contig_gbz'], region.nodes, Ls, Rs, log)
            nruns = []
            for name, off, steps, at_end in runs:
                st, fl = region.normalise(steps)
                nruns.append((name, off, st, fl, at_end))
            anon = collections.Counter(s['_steps'] for s in subpaths)
            namedc = collections.Counter(r[2] for r in nruns)
            match = anon == namedc
            if not match:
                warnings.append('named runs (%d) != gbz-base subpaths (%d) as multisets' %
                                (sum(namedc.values()), sum(anon.values())))
            # attach names to subpaths (assignment within identical walks is arbitrary)
            pool = collections.defaultdict(list)
            for r in nruns:
                pool[r[2]].append(r)
            for s in subpaths:
                lst = pool.get(s['_steps'])
                if lst:
                    r = lst.pop(0)
                    s['path_name'] = r[0]
                    s['haplotype'] = hap_key(r[0])
            byhap = collections.OrderedDict()
            for name, off, steps, fl, at_end in sorted(nruns, key=lambda r: (hap_key(r[0]), r[0], r[1])):
                hk = hap_key(name)
                e = byhap.setdefault(hk, {'haplotype': hk, 'path_fragments': set(), 'n_subpaths': 0,
                                          'n_spanning': 0, 'classes': collections.Counter(),
                                          'total_bp': 0, 'spanning_bp': [], 'visits': collections.Counter(),
                                          'max_node_visits': 0, 'distinct_ids': [],
                                          'fragment_ends_inside': 0})
                e['path_fragments'].add(name)
                e['n_subpaths'] += 1
                cls = region.classify(steps)
                e['classes'][cls] += 1
                L_ = g.walk_len(steps)
                e['total_bp'] += L_
                if cls == 'spanning':
                    e['n_spanning'] += 1
                    e['spanning_bp'].append(L_)
                c = collections.Counter(abs(x) for x in steps)
                e['max_node_visits'] = max(e['max_node_visits'], max(c.values()))
                e['visits'].update(c)
                e['distinct_ids'].append(dist_index.get(steps, '?'))
                # a named run can only start/end on an inner node where its path fragment starts/ends
                e['fragment_ends_inside'] += (region.boundary_label(steps[0]) == 'inner') + \
                    (region.boundary_label(steps[-1]) == 'inner')
            named_rows = []
            for hk, e in byhap.items():
                named_rows.append({
                    'haplotype': hk, 'path_fragments': ','.join(sorted(e['path_fragments'])),
                    'n_subpaths': e['n_subpaths'], 'n_spanning': e['n_spanning'],
                    'classes': ','.join('%s:%d' % kv for kv in sorted(e['classes'].items())),
                    'total_bp': e['total_bp'], 'spanning_bp': ','.join(map(str, e['spanning_bp'])),
                    'has_cycle': 1 if e['max_node_visits'] > 1 else 0,
                    'max_node_visits': e['max_node_visits'],
                    'n_nodes_revisited': sum(1 for v in e['visits'].values() if v > 1),
                    'distinct_ids': ','.join(e['distinct_ids']),
                    'fragment_ends_inside': e['fragment_ends_inside']})
            write_tsv(os.path.join(outdir, 'haplotype_names.tsv'), NAMED_COLS, named_rows)
            named_summary = collections.OrderedDict([
                ('n_paths_scanned', n_paths), ('n_named_runs', len(nruns)),
                ('named_runs_equal_gbz_subpaths', match),
                ('n_haplotypes_touching', len(named_rows)),
                ('n_haplotypes_single_spanning_traversal',
                 sum(1 for r in named_rows if r['n_subpaths'] == 1 and r['n_spanning'] == 1)),
                ('n_haplotypes_fragmented', sum(1 for r in named_rows if r['n_spanning'] == 0)),
                ('n_haplotypes_with_cycle', sum(1 for r in named_rows if r['has_cycle'])),
                ('n_haplotypes_multiple_subpaths', sum(1 for r in named_rows if r['n_subpaths'] > 1)),
            ])
    write_tsv(os.path.join(outdir, 'subpaths.tsv'), SUBPATH_COLS, subpaths)

    # ---- 7. topology -----------------------------------------------------------
    refc = collections.Counter(abs(s) for s in ref_steps)
    tot_w = sum(d['weight'] for d in distinct)
    span_w = sum(d['weight'] for d in distinct if d['class'] == 'spanning')
    cyc_d = [d for d in distinct if d['max_node_visits'] > 1]
    node_visit_any = collections.Counter()
    for s in subpaths:
        c = collections.Counter(abs(x) for x in s['_steps'])
        for n, k in c.items():
            if k > node_visit_any[n]:
                node_visit_any[n] = k
    topo = collections.OrderedDict([
        ('mode', mode), ('query', ' '.join(['gbz-base', 'query'] + query_args)),
        ('analysed_span_1based', [A1, B1]), ('reference_path_length_bp', ref_len),
        ('nodes', len(g.seqs)), ('edges', len(g.edges)),
        ('total_node_bp', sum(len(s) for s in g.seqs.values())),
        ('nodes_per_kb_reference', fmt_float(len(g.seqs) / (ref_len / 1000.0), 2) if ref_len else None),
        ('reference_steps', len(ref_steps)), ('reference_distinct_nodes', len(refc)),
        ('reference_nodes_visited_more_than_once', sum(1 for v in refc.values() if v > 1)),
        ('nonreference_nodes', len(region.nodes - region.ref_nodes)),
        ('nonreference_node_bp', sum(len(g.seqs[n]) for n in region.nodes - region.ref_nodes)),
        ('distinct_haplotype_walks', len(distinct)),
        ('haplotype_walks_total_weight', tot_w),
        ('individual_subpaths', len(subpaths)),
        ('spanning_distinct_walks', len(spanning)),
        ('spanning_walks_weight', span_w),
        ('fragment_walks_weight', tot_w - span_w),
        ('fragment_distinct_walks', len(distinct) - len(spanning)),
        ('haplotypes_at_boundaries_H', summary['boundaries'].get('haplotypes_at_boundaries_H')),
        ('cycles', collections.OrderedDict([
            ('distinct_walks_revisiting_a_node', len(cyc_d)),
            ('weight_of_walks_revisiting_a_node', sum(d['weight'] for d in cyc_d)),
            ('max_visits_of_one_node_by_one_walk', max([d['max_node_visits'] for d in distinct] or [0])),
            ('nodes_revisited_by_some_walk', sum(1 for v in node_visit_any.values() if v > 1)),
            ('self_loop_edges', sum(1 for a, ao, b, bo in g.edges if a == b)),
            ('reversing_edges', sum(1 for a, ao, b, bo in g.edges if ao != bo)),
            ('graph_has_directed_cycle', graph_has_directed_cycle(g.seqs, g.edges)),
        ])),
    ])
    if named_summary:
        topo['named'] = named_summary
    with open(os.path.join(outdir, 'topology.json'), 'w') as f:
        json.dump(topo, f, indent=1)
    summary['topology'] = topo
    frag = collections.OrderedDict([
        ('haplotypes_at_boundaries', summary['boundaries'].get('haplotypes_at_boundaries_H')),
        ('walks_spanning_end_to_end', span_w),
        ('walks_fragmentary', tot_w - span_w),
        ('distinct_spanning', len(spanning)),
        ('distinct_fragment', len(distinct) - len(spanning)),
    ])
    summary['fragmentation'] = frag
    summary['closest'] = {k: v.get('best') for k, v in closest.items()}
    summary['closest_non_reference'] = {k: v.get('best_non_reference') for k, v in closest.items()}

    # ---- 8. reads -------------------------------------------------------------
    if reads and os.path.exists(gaf_path):
        log.info('summarising reads')
        rs = summarise_reads(gaf_path, region, ref_len, [len(th1), len(th2)], depth, log)
        with open(os.path.join(outdir, 'reads_summary.json'), 'w') as f:
            json.dump(rs, f, indent=1)
        summary['reads'] = rs

    # ---- 9. calls -------------------------------------------------------------
    log.info('calls')
    vrows, vnote = calls_rows('vg', paths['vg_vcf'], paths['vg_truvari'], paths['vg_aardvark'],
                              contig, A1, B1, start1, end1, paths, workdir, log)
    prows, pnote = calls_rows('pangenie', paths['pg_norm'], paths['pg_truvari'], paths['pg_aardvark'],
                              contig, A1, B1, start1, end1, paths, workdir, log, already_normalised=True)
    for n in (vnote, pnote):
        if n:
            warnings.append(n)
    prows = [r for r in prows if r['called']]  # PanGenie genotypes every panel allele
    # representation-independent check: rebuild each caller's haplotypes over the span
    callhap = collections.OrderedDict()
    base = pair_distance(refseq, refseq, th1, th2)
    callhap['reference_as_call'] = base
    for caller, vcf in (('vg', paths['vg_vcf']), ('pangenie', paths['pg_norm'])):
        recs = vcf_records(vcf, contig, A1, B1, log)
        if recs is None:
            continue
        (c1, c2), info = called_haplotypes(recs, A1, B1, refseq)
        if c1 is None or c2 is None:
            callhap[caller] = {'error': 'REF mismatch while applying calls'}
            continue
        e = pair_distance(c1, c2, th1, th2)
        e.update({'slot1_len': len(c1), 'slot2_len': len(c2)})
        e.update(info)
        e['phase_reliable'] = (info['unphased_het'] <= 1 and info['phase_sets'] <= 1)
        callhap[caller] = e
        write_fasta(os.path.join(outdir, '%s.called.fa' % caller),
                    [('%s_slot1 %s:%d-%d' % (caller, contig, A1, B1), c1),
                     ('%s_slot2 %s:%d-%d' % (caller, contig, A1, B1), c2)], 80)
    summary['called_haplotypes_vs_truth'] = callhap
    write_tsv(os.path.join(outdir, 'calls.tsv'), CALL_COLS, vrows + prows)

    def ctab(rows):
        c = collections.Counter()
        for r in rows:
            if not r['called']:
                continue
            big = abs(r['len_diff']) >= 50
            c['sv50_' + r['truvari'] if big else 'small_' + r['aardvark']] += 1
        return dict(sorted(c.items()))
    summary['calls'] = {'vg': ctab(vrows), 'pangenie': ctab(prows)}
    tv = collections.Counter()
    for r in truth_rows:
        if r['source'] == 'stvar' and r['is_sv50']:
            tv['vg_' + r.get('vg_truvari', '-')] += 1
            tv['pg_' + r.get('pg_truvari', '-')] += 1
    summary['truth_sv50_status'] = dict(sorted(tv.items()))

    # ---- 10. repeat annotation ---------------------------------------------------
    rm = rmsk_overlaps(paths['rmsk'], contig, A1 - 1, B1, log)
    tr = [e for e in rm if e['class'] in ('Simple_repeat', 'Low_complexity', 'Satellite')]
    summary['repeatmasker'] = {
        'n_elements': len(rm),
        'tandem_like': [{'start1': e['start0'] + 1, 'end': e['end'], 'name': e['name'],
                         'class': e['class'], 'period': e['period']} for e in tr][:50],
        'bp_simple_repeat_period_ge7': bp_union([(e['start0'], e['end']) for e in rm
                                                 if e['class'] == 'Simple_repeat' and e['period']
                                                 and e['period'] >= 7], A1 - 1, B1),
        'bp_simple_repeat_period_1_6': bp_union([(e['start0'], e['end']) for e in rm
                                                 if e['class'] == 'Simple_repeat' and e['period']
                                                 and e['period'] < 7], A1 - 1, B1),
    }

    # annotated GFA: reference walk labelled, distinct walks named h<i>
    write_annotated_gfa(os.path.join(outdir, 'subgraph.gfa'), g, distinct, contig, A1, B1)

    summary['reference_tandem_scan'] = tandem_scan(refseq)
    summary['truth_h1_tandem_scan'] = tandem_scan(th1)

    # ---- 11. optional extras --------------------------------------------------------
    if msa:
        seqs = [(d['hap_id'] + '_w%d' % d['weight'], seq_of[d['hap_id']]) for d in spanning]
        seqs += [('truth_h1', th1), ('truth_h2', th2)]
        tot = sum(len(s) for _, s in seqs)
        if tot > 2_000_000:
            warnings.append('--msa skipped: %d bp total' % tot)
        else:
            inp = os.path.join(workdir, 'msa.in.fa')
            write_fasta(inp, seqs)
            log.info('mafft on %d sequences (%d bp)' % (len(seqs), tot))
            run([MAFFT, '--auto', '--thread', '4', '--quiet', inp], log=log,
                stdout_path=os.path.join(outdir, 'msa.fa'))
    if bandage:
        try:
            run([BANDAGE, 'image', os.path.join(outdir, 'subgraph.gfa'),
                 os.path.join(outdir, 'subgraph.png')], log=log)
        except (ToolError, OSError) as e:
            warnings.append('Bandage failed: %s' % str(e)[:200])

    summary['warnings'] = warnings
    with open(os.path.join(outdir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=1)
    shutil.rmtree(workdir, ignore_errors=True)
    log.info('done: %s' % outdir)
    return summary


def probe_reference_nodes(bd):
    """All CHM13 nodes of the probe (list of (step, start0, end0))."""
    g = bd['probe_gfa']
    ref = [w for w in g.walks if w['sample'] == 'CHM13'][0]
    pos = ref['start']
    out = []
    for s in ref['steps']:
        L = g.seqlen(s)
        out.append((s, pos, pos + L))
        pos += L
    return out


def tandem_scan(seq, max_period=200, max_len=50000, min_ident=0.8):
    """Truth-independent tandem-repeat scan of one sequence (RepeatMasker misses
    many VNTRs).  For each period p, a window of w = max(2p, 40) bp starting at i
    is repetitive when >= min_ident of its bases equal the base p downstream;
    covered bp = matching bases (and their copies p downstream) inside such windows.  Returns the best STR
    (p 1-6) and VNTR (p >= 7) period by covered bp, replacing the best period by
    its smallest divisor that covers >= 85% as much (so multiples do not win).
    Sequences longer than max_len are scanned on their central max_len bp."""
    import itertools
    s = seq.upper()
    if len(s) > max_len:
        a = (len(s) - max_len) // 2
        s = s[a:a + max_len]
    L = len(s)
    res = {'scanned_bp': L, 'min_identity': min_ident}
    cov = {}
    for p in range(1, min(max_period, L // 3) + 1):
        w = max(2 * p, 40)
        if w + p > L:
            break
        match = [1 if a == b else 0 for a, b in zip(s, s[p:])]
        n = len(match)
        pre = [0] + list(itertools.accumulate(match))
        need = min_ident * w
        diff = [0] * (n + 1)
        hit = False
        for i in range(0, n - w + 1):
            if pre[i + w] - pre[i] >= need:
                diff[i] += 1
                diff[i + w] -= 1
                hit = True
        if not hit:
            cov[p] = 0
            continue
        covered = bytearray(L)
        run_ = 0
        for j in range(n):
            run_ += diff[j]
            if run_ > 0 and match[j]:   # a matching base inside a repetitive window,
                covered[j] = 1          # and its copy one period downstream
                covered[j + p] = 1
        cov[p] = sum(covered)
    for label, lo, hi in (('str', 1, 6), ('vntr', 7, max_period)):
        cand = {p: c for p, c in cov.items() if lo <= p <= hi}
        best = max(cand.values()) if cand else 0
        if best == 0:
            res[label] = None
            continue
        top = min(q for q, c in cand.items() if c == best)
        p = min(q for q in cand if top % q == 0 and cand[q] >= 0.85 * best)
        res[label] = {'period': p, 'covered_bp': cand[p], 'covered_frac': round(cand[p] / L, 4)}
    return res


def bp_union(iv, a0, b0):
    iv = sorted((max(a0, s), min(b0, e)) for s, e in iv if e > a0 and s < b0)
    tot, cur_s, cur_e = 0, None, None
    for s, e in iv:
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        tot += cur_e - cur_s
    return tot


def write_annotated_gfa(path, g, distinct, contig, A1, B1):
    with open(path, 'w') as f:
        f.write('H\tVN:Z:1.1\tRS:Z:CHM13\n')
        for n in sorted(g.seqs):
            f.write('S\t%d\t%s\n' % (n, g.seqs[n]))
        for a, ao, b, bo in g.edges:
            f.write('L\t%d\t%s\t%d\t%s\t0M\n' % (a, ao, b, bo))
        for d in distinct:
            if d['is_reference']:
                f.write('W\tCHM13\t0\t%s\t%d\t%d\t%s\tWT:i:%d\n' %
                        (contig, A1 - 1, B1, steps_to_walk(d['steps']), d['weight']))
        for d in distinct:
            if d['is_reference']:
                continue
            tags = 'WT:i:%d' % d['weight']
            if d.get('cigar') and not str(d['cigar']).startswith('NA'):
                tags += '\tCG:Z:' + d['cigar']
            f.write('W\t%s\t%d\t%s\t0\t%d\t%s\t%s\n' %
                    (d['hap_id'], 0, d['class'], d['length_bp'], steps_to_walk(d['steps']), tags))


# ----------------------------------------------------------------------------

def print_summary(s):
    b = s.get('boundaries', {})
    t = s.get('topology', {})
    print('%s requested %d-%d, analysed %d-%d (%d bp, mode %s)' % (
        s['contig'], s['requested_interval_1based'][0], s['requested_interval_1based'][1],
        s['analysed_span_1based'][0], s['analysed_span_1based'][1], s['analysed_span_bp'], s['mode']))
    print('  graph: %d nodes, %d edges, %.1f nodes/kb; %d distinct walks (total weight %d), '
          '%d spanning (weight %d)' % (t['nodes'], t['edges'], t['nodes_per_kb_reference'] or 0,
                                       t['distinct_haplotype_walks'], t['haplotype_walks_total_weight'],
                                       t['spanning_distinct_walks'], t['spanning_walks_weight']))
    c = t['cycles']
    print('  cycles: %d distinct walks revisit a node (max %d visits); directed cycle in graph: %s' % (
        c['distinct_walks_revisiting_a_node'], c['max_visits_of_one_node_by_one_walk'],
        c['graph_has_directed_cycle']))
    tr = s['truth']
    print('  truth: hap1 %d bp, hap2 %d bp, CHM13 %d bp; %d records (%d SV>=50); SV-benchmark coverage '
          'of requested interval %.2f' % (tr['hap1_len'], tr['hap2_len'], tr['reference_len'],
                                         tr['n_records_in_span'], tr['n_sv50_stvar_records'],
                                         tr['requested_interval_frac_in_sv_benchmark']))
    for k, v in s['closest'].items():
        if v:
            print('  closest to %s: %s (weight %d) edit %d, len diff %+d, identity %.4f' % (
                k, v['hap_id'], v['weight'], v['edit_distance'], v['length_diff_panel_minus_truth'],
                v['identity']))
    if 'reads' in s:
        r = s['reads']
        print('  reads: %d (expected %.0f from truth haps, %.0f from CHM13 span); MAPQ<5 %.3f; '
              'inside subgraph %.3f' % (r['n_alignments'], r['expected_from_truth_haplotypes'] or 0,
                                        r['expected_from_reference_span'] or 0,
                                        r['frac_mapq_lt5'] or 0, r['frac_path_inside_subgraph'] or 0))
    print('  calls: vg %s | pangenie %s' % (s['calls']['vg'], s['calls']['pangenie']))
    ch = s.get('called_haplotypes_vs_truth', {})
    parts = []
    for k, v in ch.items():
        if 'total' in v:
            parts.append('%s %d+%d=%d%s' % (k, v['d_h1'], v['d_h2'], v['total'],
                                           '' if v.get('phase_reliable', True) else ' (phase unreliable)'))
    if parts:
        print('  called haplotypes, edit distance to truth h1+h2: ' + '; '.join(parts))
    if 'named' in t:
        n = t['named']
        print('  named: %d haplotypes touch the span; %d single spanning traversal, %d fragmented, '
              '%d with a cycle; runs==gbz subpaths: %s' % (
                  n['n_haplotypes_touching'], n['n_haplotypes_single_spanning_traversal'],
                  n['n_haplotypes_fragmented'], n['n_haplotypes_with_cycle'],
                  n['named_runs_equal_gbz_subpaths']))
    print('  truth SV>=50 status: %s' % s['truth_sv50_status'])
    ts = s.get('reference_tandem_scan', {})
    if ts:
        f = lambda x: ('p%d over %d bp (%.0f%%)' % (x['period'], x['covered_bp'], 100 * x['covered_frac'])) if x else 'none'
        print('  CHM13 tandem scan: VNTR %s; STR %s' % (f(ts.get('vntr')), f(ts.get('str'))))
    for w in s['warnings']:
        print('  WARNING: ' + w)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('contig')
    ap.add_argument('start', type=int, help='1-based inclusive')
    ap.add_argument('end', type=int, help='1-based inclusive')
    ap.add_argument('outdir')
    ap.add_argument('--pad', type=int, default=0, help='bp added to both sides before choosing boundaries')
    ap.add_argument('--no-reads', action='store_true', help='skip the GAF-base read query')
    ap.add_argument('--names', action='store_true',
                    help='name haplotype subpaths via `vg paths -A` on the contig GBZ (~40 s on chr20)')
    ap.add_argument('--anchors', help='fix the boundary nodes, e.g. 114943418+:114943990+ (mode between)')
    ap.add_argument('--mode', choices=['between', 'interval'], default='between')
    ap.add_argument('--context', type=int, default=0, help='--context for --mode interval')
    ap.add_argument('--depth', type=float, default=30.0, help='sequencing depth for expected read count')
    ap.add_argument('--window', type=int, default=20000, help='initial boundary search window (bp)')
    ap.add_argument('--max-span', type=int, default=250000, help='refuse spans longer than this')
    ap.add_argument('--limit', type=int, default=300000, help='gbz-base --between node limit')
    ap.add_argument('--max-align-cells', type=float, default=4e8,
                    help='skip CIGAR traceback when len(ref)*len(hap) exceeds this')
    ap.add_argument('--msa', action='store_true', help='mafft MSA of spanning walks + truth haps -> msa.fa')
    ap.add_argument('--bandage', action='store_true', help='render subgraph.png with Bandage')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args(argv)
    try:
        s = analyse(a.contig, a.start, a.end, a.outdir, pad=a.pad, reads=not a.no_reads, names=a.names,
                    mode=a.mode, context=a.context, depth=a.depth, msa=a.msa, bandage=a.bandage,
                    window=a.window, max_span=a.max_span, limit=a.limit,
                    max_align_cells=a.max_align_cells, verbose=not a.quiet,
                    anchors=a.anchors.split(':') if a.anchors else None)
    except ToolError as e:
        print('ERROR: %s' % e, file=sys.stderr)
        return 2
    print_summary(s)
    return 0


if __name__ == '__main__':
    sys.exit(main())
