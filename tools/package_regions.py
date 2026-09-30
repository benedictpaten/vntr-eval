#!/usr/bin/env python3
"""package_regions.py -- build the region package (regions/<region_id>/...) for a locus list.

    python3 tools/package_regions.py build  [--loci FILE] [--out DIR] [--work DIR] [--pad 200]
                                            [--max-span 250000] [--jobs 3] [--only ID,ID] [--force]
                                            [--anchor-mode truth|snarl] [--snarl-cache DIR]
    python3 tools/package_regions.py hprc-fetch [--loci FILE] [--work DIR] [--only ID,ID]
    python3 tools/package_regions.py hprc   [--out DIR] [--work DIR] [--jobs 1] [--only ID,ID]
                                            [--force]
    python3 tools/package_regions.py index  [--out DIR]
    python3 tools/package_regions.py prune-hprc [--out DIR] [--max-mb 150] [--keep-strata hotspot_vntr]
    python3 tools/package_regions.py bed    [--vntr-regions FILE] [--out DIR] [--min-len 1000]
    python3 tools/package_regions.py validate [--out DIR]

`build --anchor-mode snarl` (truth-free; tools/snarl_anchors.py) takes the anchors from the snarl
decomposition instead: the boundaries of the smallest snarl, or run of consecutive snarls of one
chain, enclosing the padded interval (`snarl_anchors.py build` caches it once per contig). A locus
whose enclosing span exceeds --max-span is skipped (no fallback), and loci whose spans overlap are
merged into one region named after the first (region.json 'merged_loci'), with the anchors of
the union. Everything after anchor choice is the same.

`build` (default --anchor-mode truth) runs, per locus: anchor choice and extraction with tools/region.py (anchors are CHM13
nodes that every hap32 path visits once and no truth record touches; the interval is padded
by --pad bp on each side first, so the anchor-to-anchor span holds the whole locus plus at
least --pad bp of flank; when that rule finds no anchors, or only ones more than --max-span
apart, the nearest best-covered nodes are used and the haplotypes that miss them become
fragments), then one `vg paths -A` pass per contig to name every haplotype subpath (the named
runs must equal gbz-base's anonymous subpaths as multisets), then writes hap32.fa,
hap32.fragments.fa, mc.gfa, truth.fa, calls.tsv and region.json. Every mc.gfa path is
re-spelled from its S lines and asserted equal to its hap32.fa sequence.

`hprc` adds hprc.fa.gz: the same anchor-to-anchor span for every sample haplotype of the full
HPRC v2.1 graph (config.data_paths()['hprc_gbz']; HG002 is not in it), plus CHM13 and GRCh38.
The full graph has the hap32 node ids but the hap32 anchors do not separate it, so each region
is an `--interval` query with context, cut at the anchors, with the reference-sense paths
subtracted exactly (see Stage D below). Each query loads the 5.7 GB GBZ (~45-65 s, ~11 GB
RSS), so run one at a time; `hprc-fetch` runs the queries next to `build`, as soon as a
locus has its anchors, and `hprc` reuses them. Spans over HPRC_MAX_SPAN are skipped.

`index` (re)writes regions/regions.tsv from the region.json files; `bed` writes
regions/vntr_ge1kb.bed from the census VNTR table; `validate` re-reads every package and
checks the invariants (counts, files, path spelling, flanks).

--loci accepts the census strata.tsv (locus_id contig start(0-based) end stratum ...), a
regions.tsv written by this script (region_id contig core_start core_end stratum ...), or a
BED file (contig start0 end name [period motif vg_fp vg_fn]; stratum "target").

Data and tool locations come from config.py. Nothing is written outside --out and --work.
"""

import argparse
import collections
import datetime
import gzip
import json
import multiprocessing
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
import traceback

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402
import region as vr  # noqa: E402

PILOT = {
    'L009656': 'pilot: alleles present but misaligned',
    'L014297': 'pilot: alleles present but misaligned',
    'L015415': 'pilot: alleles present but misaligned (both HG002 haplotypes exactly in the panel)',
    'L005990': 'pilot: alleles present but misaligned',
    'L012184': 'pilot: under-aligned (parallel duplicated path)',
    'L011138': 'pilot: under-aligned',
    'L012272': 'pilot: under-aligned',
    'L002013': 'pilot: alleles absent from the panel',
    'L016870': 'pilot: negative control',
}

STRATUM_ROLE = {
    'hotspot_vntr': 'test',
    'hotspot_other': 'test',
    'control_vntr_matched': 'control (must not regress)',
    'control_vntr_correct': 'control (must not regress)',
    'control_nontr_sv': 'control (must not regress)',
    'target': 'target (truth-free list)',
}

INDEX_COLS = ['region_id', 'stratum', 'contig', 'core_start', 'core_end', 'span_start', 'span_end',
              'anchor_left', 'anchor_right', 'period', 'motif', 'copies', 'n_hap32', 'n_hprc',
              'vg_fp', 'vg_fn', 'pg_fp', 'pg_fn', 'in_benchmark', 'notes']

STEP = re.compile(r'([<>])(\d+)')
MAX_SPAN = 250000        # anchor-to-anchor limit before falling back to nearer anchors
HPRC_MAX_SPAN = 100000   # longer spans are not queried in the full graph
_RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_RC)[::-1]


def now():
    return datetime.datetime.now().isoformat(timespec='seconds')


def default_work():
    return os.path.join(config.WORK_DIR, 'package')


# ----------------------------------------------------------------------------
# Locus lists

def _num(x, typ=float):
    try:
        if x in (None, '', '.', 'NA', 'nan'):
            return None
        v = typ(x)
        return v
    except ValueError:
        return None


def load_loci(path):
    """Return [locus dict] with region_id, contig, core_start, core_end (1-based inclusive),
    stratum and whatever annotation the file carries (kept under 'census')."""
    loci = []
    if path.endswith('.bed') or path.endswith('.bed.gz'):
        op = gzip.open if path.endswith('.gz') else open
        with op(path, 'rt') as f:
            for line in f:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                x = line.rstrip('\n').split('\t')
                d = {'region_id': x[3] if len(x) > 3 else '%s_%s_%s' % (x[0], x[1], x[2]),
                     'contig': x[0], 'core_start': int(x[1]) + 1, 'core_end': int(x[2]),
                     'stratum': 'target', 'census': {}}
                for i, k in ((4, 'period'), (5, 'motif'), (6, 'vg_fp'), (7, 'vg_fn')):
                    if len(x) > i:
                        d['census'][k] = x[i]
                loci.append(d)
        return loci
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            if not line.strip():
                continue
            r = dict(zip(header, line.rstrip('\n').split('\t')))
            if 'locus_id' in r:            # census strata.tsv / loci.tsv (start is 0-based)
                d = {'region_id': r['locus_id'], 'contig': r['contig'],
                     'core_start': int(r['start']) + 1, 'core_end': int(r['end']),
                     'stratum': r.get('stratum', 'target')}
            elif 'core_start' in r:        # regions.tsv
                d = {'region_id': r['region_id'], 'contig': r['contig'],
                     'core_start': int(r['core_start']), 'core_end': int(r['core_end']),
                     'stratum': r.get('stratum', 'target')}
            else:
                raise SystemExit('%s: need locus_id/start/end or region_id/core_start/core_end' % path)
            d['census'] = r
            loci.append(d)
    return loci


# ----------------------------------------------------------------------------
# Stage A: anchors + extraction (tools/region.py)

def fallback_anchors(contig, a0, b0, log, window=20000, max_window=80000):
    """Anchors when region.py's rule finds none within reach, or only ones more than
    MAX_SPAN apart. region.py wants coverage == the probe's maximum; here each side takes,
    within the probe window, the node with the highest coverage among those that keep the
    rest of the rule (each haplotype subpath visits it at most once, CHM13 once, no truth
    record touches it), the nearest one on ties. Haplotypes that miss an anchor become
    fragments. Returns (L, R, note) as signed node ids."""
    paths = config.data_paths(contig)
    clen = vr.contig_length(paths['ref_fa'], contig)
    w = window
    while True:
        lo0, hi0 = max(0, a0 - w), min(clen, b0 + w)
        nodes, cov, maxvisit, H, g = vr.probe_reference(paths['gbz_db'], contig, lo0, hi0, log)
        refc = collections.Counter(abs(s) for s, _, _ in nodes)
        spans = vr.truth_spans(paths, contig, lo0 + 1, hi0, log)
        elig = [(i, s, st, en) for i, (s, st, en) in enumerate(nodes)
                if maxvisit[abs(s)] <= 1 and refc[abs(s)] == 1 and not vr.span_hits(spans, st + 1, en)]
        left = [e for e in elig if e[3] <= a0]
        right = [e for e in elig if e[2] >= b0]
        if left and right:
            cl = max(cov[abs(e[1])] for e in left)
            cr = max(cov[abs(e[1])] for e in right)
            L = [e for e in left if cov[abs(e[1])] == cl][-1]
            R = [e for e in right if cov[abs(e[1])] == cr][0]
            note = ('fallback anchors: haplotype coverage L %d, R %d of %d (nearest best-covered '
                    'nodes within %d bp of the padded locus)' % (cl, cr, H, w))
            return L[1], R[1], note
        if (lo0 == 0 and hi0 == clen) or w >= max_window:
            raise vr.ToolError('fallback: no anchor within %d bp (L %s, R %s)' %
                               (w, 'ok' if left else 'missing', 'ok' if right else 'missing'))
        w *= 2


def stage_extract(args):
    locus, work, pad, force, max_span = args
    rid = locus['region_id']
    tdir = os.path.join(work, rid, 'toolkit')
    sj = os.path.join(tdir, 'summary.json')
    t0 = time.time()
    if not force and os.path.exists(sj):
        return rid, 'cached', None, 0.0
    if os.path.exists(tdir):
        shutil.rmtree(tdir)
    kw = dict(pad=pad, reads=False, names=False, verbose=False)
    note = None
    try:
        if locus.get('anchors'):
            # snarl mode: the anchors are fixed; no fallback
            note = locus['anchor_note']
            vr.analyse(locus['contig'], locus['core_start'], locus['core_end'], tdir,
                       anchors=locus['anchors'], max_span=max_span, **kw)
            s = json.load(open(sj))
            s.setdefault('warnings', []).append(note)
            s['package_anchor_note'] = note
            json.dump(s, open(sj, 'w'), indent=1)
            return rid, 'ok', note, time.time() - t0
        try:
            vr.analyse(locus['contig'], locus['core_start'], locus['core_end'], tdir,
                       max_span=max_span, **kw)
        except vr.ToolError as e:
            # the default rule wants every haplotype on both anchors; one haplotype that
            # bypasses CHM13 for a long stretch pushes the anchors far out (L000601: 321 kb)
            # or out of reach. Then take the nearest well-covered nodes instead.
            if 'no boundary nodes' not in str(e) and 'exceeds --max-span' not in str(e):
                raise
            first = str(e).split('\n')[0][:200]
            log = vr.Log(os.path.join(work, rid, 'fallback.log'), verbose=False)
            a0 = max(0, locus['core_start'] - 1 - pad)
            b0 = locus['core_end'] + pad
            L, R, note = fallback_anchors(locus['contig'], a0, b0, log)
            note += ' [default rule: %s]' % first
            shutil.rmtree(tdir, ignore_errors=True)
            vr.analyse(locus['contig'], locus['core_start'], locus['core_end'], tdir,
                       anchors=(L, R), max_span=max_span, **kw)
        if note:
            s = json.load(open(sj))
            s.setdefault('warnings', []).append(note)
            s['package_anchor_note'] = note
            json.dump(s, open(sj, 'w'), indent=1)
        return rid, 'ok', note, time.time() - t0
    except Exception as e:  # noqa: BLE001 -- report every failure per locus
        with open(os.path.join(work, rid, 'error.txt'), 'w') as f:
            f.write(traceback.format_exc())
        return rid, 'error', '%s: %s' % (type(e).__name__, str(e)[:300]), time.time() - t0


# ----------------------------------------------------------------------------
# Stage B: names from `vg paths -A` (one pass per contig)

def scan_contig(contig, infos, log_path):
    """infos: [{'id', 'nodes' (set), 'L', 'R'}]. Streams `vg paths -A` over the contig GBZ
    and returns ({id: [(path_name, offset, steps, reaches_path_end)]}, n_paths). A run
    begins at the path start or where the path enters a locus through one of its anchors
    (a between-subgraph can only be entered through an anchor), and extends while the nodes
    stay in the locus node set -- the same rule as region.runs_in_walk."""
    gbz = config.data_paths(contig)['contig_gbz']
    by_node = collections.defaultdict(list)
    for i, inf in enumerate(infos):
        by_node[abs(inf['L'])].append(i)
        by_node[abs(inf['R'])].append(i)
    tok = re.compile(r'[<>](%s)(?![0-9])' % '|'.join(str(n) for n in sorted(by_node)))
    runs = {inf['id']: [] for inf in infos}
    cmd = [config.VG, 'paths', '-x', gbz, '-A']
    t0 = time.time()
    with open(log_path, 'a') as lg:
        lg.write('$ %s\n' % ' '.join(cmd))
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=lg, text=True, env=config.tool_env(),
                             bufsize=1 << 20)
        n_paths = 0
        for line in p.stdout:
            if line.startswith('@'):
                continue
            f = line.split('\t', 6)
            if len(f) < 6:
                continue
            n_paths += 1
            name, walk = f[0], f[5].rstrip('\n')
            starts = collections.defaultdict(set)
            m0 = STEP.match(walk)
            if m0:
                n0 = int(m0.group(2))
                for i, inf in enumerate(infos):
                    if n0 in inf['nodes']:
                        starts[i].add(0)
            for m in tok.finditer(walk):
                a = m.start()
                n = int(m.group(1))
                for i in by_node[n]:
                    if a == 0:
                        starts[i].add(0)
                        continue
                    j = a - 1
                    while j >= 0 and walk[j].isdigit():
                        j -= 1
                    prev = int(walk[j + 1:a])
                    if prev not in infos[i]['nodes']:
                        starts[i].add(a)
            for i, offs in starts.items():
                ns = infos[i]['nodes']
                for a in sorted(offs):
                    steps = []
                    pos = a
                    while True:
                        m = STEP.match(walk, pos)
                        if not m:
                            break
                        n = int(m.group(2))
                        if n not in ns:
                            break
                        steps.append(n if m.group(1) == '>' else -n)
                        pos = m.end()
                    if steps:
                        runs[infos[i]['id']].append((name, a, tuple(steps), pos >= len(walk)))
        p.wait()
        lg.write('  %s: %d paths scanned in %.1f s, exit %d\n' % (contig, n_paths, time.time() - t0,
                                                                   p.returncode))
    if p.returncode != 0:
        raise vr.ToolError('vg paths -A failed on %s' % gbz)
    return runs, n_paths


def stage_names(args):
    contig, infos, work = args
    t0 = time.time()
    try:
        runs, n_paths = scan_contig(contig, infos, os.path.join(work, 'names.%s.log' % contig))
    except Exception as e:  # noqa: BLE001
        return contig, None, '%s: %s' % (type(e).__name__, e), time.time() - t0
    for inf in infos:
        out = os.path.join(work, inf['id'], 'named_runs.json')
        with open(out, 'w') as f:
            json.dump({'contig': contig, 'n_paths_scanned': n_paths,
                       'runs': [[n, off, list(st), end] for n, off, st, end in runs[inf['id']]]}, f)
    return contig, n_paths, None, time.time() - t0


# ----------------------------------------------------------------------------
# Stage C: write the package

def read_fasta(path):
    recs = []
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        name, desc, buf = None, '', []
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    recs.append((name, desc, ''.join(buf)))
                h = line[1:].split(None, 1)
                name, desc, buf = h[0], (h[1] if len(h) > 1 else ''), []
            else:
                buf.append(line.strip())
        if name is not None:
            recs.append((name, desc, ''.join(buf)))
    return recs


def write_fasta(path, recs):
    """recs: (name, description, sequence); one line per sequence."""
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'wt') as f:
        for name, desc, seq in recs:
            f.write('>%s%s\n%s\n' % (name, (' ' + desc) if desc else '', seq))


def name_key(name):
    parts = name.split('#')
    rank = 0 if parts[0] == 'CHM13' else 1 if parts[0] == 'GRCh38' else 2
    nums = [int(x) if x.isdigit() else x for x in re.split(r'(\d+)', name)]
    return (rank, parts[0], [str(n).zfill(12) if isinstance(n, int) else n for n in nums])


def canon_edge(a, b):
    """Edge between consecutive signed steps a -> b, canonical over both strands."""
    e1 = (abs(a), '+' if a > 0 else '-', abs(b), '+' if b > 0 else '-')
    e2 = (abs(b), '-' if b > 0 else '+', abs(a), '-' if a > 0 else '+')
    return min(e1, e2)


def spell(seqs, steps):
    return ''.join(seqs[s] if s > 0 else revcomp(seqs[-s]) for s in steps)


def gfa_spell_check(path, want):
    """Re-read a GFA 1.0 file; spell every P line from S lines; compare with want{name: seq}.
    Returns (n_checked, [problems])."""
    seqs, paths = {}, {}
    with open(path) as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            if x[0] == 'S':
                seqs[x[1]] = x[2]
            elif x[0] == 'P':
                steps = [(s[:-1], s[-1]) for s in x[2].split(',')]
                paths[x[1]] = steps
    probs = []
    for name, seq in want.items():
        if name not in paths:
            probs.append('missing P line %s' % name)
            continue
        got = ''.join(seqs[n] if o == '+' else revcomp(seqs[n]) for n, o in paths[name])
        if got != seq:
            probs.append('P line %s spells %d bp, hap32.fa has %d bp' % (name, len(got), len(seq)))
    extra = set(paths) - set(want)
    if extra:
        probs.append('P lines without a sequence: %s' % ','.join(sorted(extra)[:5]))
    return len(want), probs


def signed_to_anchor(s):
    return '%d%s' % (abs(s), '+' if s > 0 else '-')


def nrun_count(seq, minlen=10):
    return len(re.findall('N{%d,}' % minlen, seq.upper()))


def bed_frac(bed, a0, b0):
    return vr.bed_coverage(bed, a0, b0)


def write_calls(tdir, out):
    """calls.tsv: vg and PanGenie alleles and truth records over the span, with truvari status."""
    cols = ['source', 'chrom', 'pos', 'end', 'id', 'ref_len', 'alt_len', 'len_diff', 'gt', 'called',
            'qual', 'filter', 'truvari', 'aardvark', 'vg_truvari', 'pg_truvari', 'applied_h1',
            'applied_h2', 'in_sv_benchmark', 'svtype', 'TRFperiod', 'MatchId', 'overlaps_span_edge']
    rows = []

    def read(path):
        if not os.path.exists(path):
            return []
        with open(path) as f:
            h = f.readline().rstrip('\n').split('\t')
            return [dict(zip(h, l.rstrip('\n').split('\t'))) for l in f if l.strip()]
    for r in read(os.path.join(tdir, 'calls.tsv')):
        rows.append({'source': 'vg' if r['caller'] == 'vg' else 'pangenie', 'chrom': r['chrom'],
                     'pos': r['pos'], 'end': r['end'], 'id': r['id'], 'ref_len': r['ref_len'],
                     'alt_len': r['alt_len'], 'len_diff': r['len_diff'], 'gt': r['gt'],
                     'called': r['called'], 'qual': r['qual'], 'filter': r['filter'],
                     'truvari': r['truvari'], 'aardvark': r['aardvark'], 'MatchId': r.get('MatchId', '.'),
                     'overlaps_span_edge': r['overlaps_span_edge']})
    for r in read(os.path.join(tdir, 'truth_records.tsv')):
        rows.append({'source': 'truth_' + r['source'], 'chrom': r['chrom'], 'pos': r['pos'],
                     'end': r['end'], 'id': '.', 'ref_len': r['ref_len'], 'alt_len': r['alt_lens'],
                     'len_diff': r['max_allele_len_diff'], 'gt': r['gt'], 'called': '.',
                     'vg_truvari': r.get('vg_truvari', '-'), 'pg_truvari': r.get('pg_truvari', '-'),
                     'applied_h1': r['applied_h1'], 'applied_h2': r['applied_h2'],
                     'in_sv_benchmark': r['in_sv_benchmark'], 'svtype': r['svtype'],
                     'TRFperiod': r['TRFperiod'], 'overlaps_span_edge': '0' if r['within_span'] == '1' else '1'})
    order = {'truth_stvar': 0, 'truth_smvar': 1, 'vg': 2, 'pangenie': 3}
    rows.sort(key=lambda r: (int(r['pos']), order.get(r['source'], 9)))
    with open(out, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join('.' if r.get(c) in (None, '') else str(r.get(c)) for c in cols) + '\n')
    n = collections.Counter()
    for r in rows:
        if r['source'] in ('vg', 'pangenie') and r['called'] == '1':
            try:
                big = abs(int(r['len_diff'])) >= 50
            except ValueError:
                big = False
            if big:
                n['%s_sv50_%s' % (r['source'], r['truvari'])] += 1
    return len(rows), dict(sorted(n.items()))


def census_fields(locus):
    c = locus.get('census', {})
    get = lambda k, t=float: _num(c.get(k), t)  # noqa: E731
    return {
        'period': get('period', float), 'motif': c.get('motif') or None, 'copies': get('copies', float),
        'vg_fp': get('vg_fp', int), 'vg_fn': get('vg_fn', int), 'pg_fp': get('pg_fp', int),
        'pg_fn': get('pg_fn', int),
    }


def stage_write(args):
    locus, work, out, pad = args
    rid = locus['region_id']
    tdir = os.path.join(work, rid, 'toolkit')
    rdir = os.path.join(out, rid)
    try:
        s = json.load(open(os.path.join(tdir, 'summary.json')))
        nr = json.load(open(os.path.join(work, rid, 'named_runs.json')))
        g = vr.Gfa(open(os.path.join(tdir, 'subgraph.query.gfa')).read())
        seqs = g.seqs
        L = vr.parse_steps(s['boundaries']['L']['node'])[0]
        R = vr.parse_steps(s['boundaries']['R']['node'])[0]
        A1, B1 = s['analysed_span_1based']
        contig = s['contig']
        notes = []
        # the CHM13 walk, from the named CHM13 run (it is in CHM13 orientation by definition)
        ref_runs = [r for r in nr['runs'] if r[0].split('#')[0] == 'CHM13']
        if len(ref_runs) != 1:
            raise vr.ToolError('%d CHM13 runs in the span' % len(ref_runs))
        ref_steps = tuple(ref_runs[0][2])
        region = vr.Region(g, ref_steps, A1 - 1, L, R)
        # gbz-base subpaths (anonymous) vs named runs: must be equal as multisets
        g2 = vr.Gfa(open(os.path.join(tdir, 'subpaths.gfa')).read())
        anon = collections.Counter(region.normalise(w['steps'])[0] for w in g2.walks)
        named = []
        for name, off, steps, at_end in nr['runs']:
            st, flipped = region.normalise(tuple(steps))
            named.append({'name': name, 'steps': st, 'flipped': flipped, 'cls': region.classify(st)})
        namedc = collections.Counter(r['steps'] for r in named)
        if anon != namedc:
            raise vr.ToolError('named runs (%d) != gbz-base subpaths (%d) as multisets' %
                               (sum(namedc.values()), sum(anon.values())))
        if ref_steps[0] != L or ref_steps[-1] != R:
            raise vr.ToolError('CHM13 run does not go from L to R')
        refseq = spell(seqs, ref_steps)
        Lseq, Rseq = spell(seqs, (L,)), spell(seqs, (R,))
        # ---- hap32.fa (spanning) and hap32.fragments.fa -------------------
        span_recs, frag_recs, used = [], [], collections.Counter()
        anchor_only = []      # a path that ends on an anchor, or starts on one, and no more
        for r in sorted(named, key=lambda r: (name_key(r['name']), r['cls'])):
            seq = spell(seqs, r['steps'])
            if r['cls'] != 'spanning' and all(abs(x) in (abs(L), abs(R)) for x in r['steps']):
                anchor_only.append('%s(%s)' % (r['name'], vr.steps_to_walk(r['steps'])))
                continue
            if r['cls'] == 'spanning':
                nm = r['name']
                used[nm] += 1
                if used[nm] > 1:
                    nm = '%s:%d' % (nm, used[nm])
                    notes.append('path %s spans twice' % r['name'])
                assert seq.startswith(Lseq) and seq.endswith(Rseq)
                span_recs.append((nm, r, seq))
            else:
                nm = '%s:%s' % (r['name'], r['cls'])
                used[nm] += 1
                if used[nm] > 1:
                    nm = '%s:%d' % (nm, used[nm])
                frag_recs.append((nm, r, seq))
        if not span_recs:
            raise vr.ToolError('no spanning paths')
        os.makedirs(rdir, exist_ok=True)
        write_fasta(os.path.join(rdir, 'hap32.fa'),
                    [(nm, 'len=%d' % len(sq), sq) for nm, r, sq in span_recs])
        write_fasta(os.path.join(rdir, 'hap32.fragments.fa'),
                    [(nm, 'len=%d first=%s last=%s' % (len(sq), vr.steps_to_walk(r['steps'][:1]),
                                                       vr.steps_to_walk(r['steps'][-1:])), sq)
                     for nm, r, sq in frag_recs])
        # ---- mc.gfa: the subgraph the spanning paths induce -------------------
        nodes, edges = set(), set()
        for nm, r, sq in span_recs:
            st = r['steps']
            nodes.update(abs(x) for x in st)
            for a, b in zip(st, st[1:]):
                edges.add(canon_edge(a, b))
        gpath = os.path.join(rdir, 'mc.gfa')
        with open(gpath, 'w') as f:
            f.write('H\tVN:Z:1.0\n')
            for n in sorted(nodes):
                f.write('S\t%d\t%s\n' % (n, seqs[n]))
            for a, ao, b, bo in sorted(edges):
                f.write('L\t%d\t%s\t%d\t%s\t0M\n' % (a, ao, b, bo))
            for nm, r, sq in span_recs:
                f.write('P\t%s\t%s\t*\n' % (nm, ','.join('%d%s' % (abs(x), '+' if x > 0 else '-')
                                                          for x in r['steps'])))
        n_chk, probs = gfa_spell_check(gpath, {nm: sq for nm, r, sq in span_recs})
        if probs:
            raise AssertionError('mc.gfa path check failed: ' + '; '.join(probs[:5]))
        # ---- truth.fa -------------------------------------------------------------
        th = {}
        for h in (1, 2):
            recs = read_fasta(os.path.join(tdir, 'truth.hap%d.fa' % h))
            th[h] = recs[0][2].upper()
        truth_ok = all(t.startswith(Lseq) and t.endswith(Rseq) for t in th.values())
        if not truth_ok:
            notes.append('truth haplotypes do not share the anchor sequence')
        write_fasta(os.path.join(rdir, 'truth.fa'),
                    [('HG002#%d' % h, 'T2T-Q100 v1.1 GT slot %d, %s:%d-%d, len=%d; EVALUATION ONLY'
                      % (h, contig, A1, B1, len(th[h])), th[h]) for h in (1, 2)])
        # ---- calls.tsv --------------------------------------------------------------
        n_calls, call_counts = write_calls(tdir, os.path.join(rdir, 'calls.tsv'))
        # ---- region.json ------------------------------------------------------------
        paths = config.data_paths(contig)
        sv_bed = vr.load_bed(paths['stvar_bed'], contig)
        core_frac = bed_frac(sv_bed, locus['core_start'] - 1, locus['core_end'])
        span_frac = bed_frac(sv_bed, A1 - 1, B1)
        cf = census_fields(locus)
        lens = [len(sq) for _, _, sq in span_recs]
        distinct = len(set(sq for _, _, sq in span_recs))
        nlen = [len(seqs[n]) for n in nodes]
        topo = s.get('topology', {})
        n_with_N = sum(1 for _, _, sq in span_recs if nrun_count(sq))
        grch38 = [nm for nm, _, _ in span_recs if nm.startswith('GRCh38')]
        if locus['region_id'] in PILOT:
            notes.insert(0, PILOT[locus['region_id']])
        if s.get('package_anchor_note'):
            notes.append(s['package_anchor_note'])
        nonspan_haps = sorted(set(vr.hap_key(r['name']) for r in named) -
                              set(vr.hap_key(r['name']) for _, r, _ in span_recs), key=name_key)
        if nonspan_haps:
            notes.append('%d haplotype(s) touch the span without spanning it: %s' %
                         (len(nonspan_haps), ','.join(nonspan_haps)))
        if frag_recs:
            notes.append('%d non-spanning path pieces (hap32.fragments.fa)' % len(frag_recs))
        if anchor_only:
            notes.append('%d path piece(s) are an anchor node alone (a path ends or starts on it; '
                         'left out of hap32.fragments.fa)' % len(anchor_only))
        if n_with_N:
            notes.append('%d spanning sequences carry N runs >= 10 bp' % n_with_N)
        if not grch38:
            notes.append('GRCh38 does not span')
        H = s['boundaries'].get('haplotypes_at_boundaries_H')
        if H is not None and H < 34:
            notes.append('anchor coverage %d of 34 haplotype paths' % H)
        for w in s.get('warnings', []):
            if 'interval is' in w and 'SV benchmark' in w:
                notes.append('core %.0f%% inside the SV benchmark' % (100 * core_frac))
        c = locus.get('census', {})
        rj = collections.OrderedDict()
        rj['region_id'] = rid
        rj['stratum'] = locus['stratum']
        rj['anchor_mode'] = locus.get('anchor_mode', 'truth')
        if locus.get('merged_loci'):
            rj['merged_loci'] = locus['merged_loci']
        rj['role'] = STRATUM_ROLE.get(locus['stratum'], locus['stratum'])
        rj['contig'] = contig
        rj['core_start'] = locus['core_start']
        rj['core_end'] = locus['core_end']
        rj['span_start'] = A1
        rj['span_end'] = B1
        rj['anchor_left'] = signed_to_anchor(L)
        rj['anchor_right'] = signed_to_anchor(R)
        rj['period'] = cf['period']
        rj['motif'] = cf['motif']
        rj['copies'] = cf['copies']
        rj['n_hap32'] = len(span_recs)
        rj['n_hprc'] = 0
        rj['vg_fp'] = cf['vg_fp']
        rj['vg_fn'] = cf['vg_fn']
        rj['pg_fp'] = cf['pg_fp']
        rj['pg_fn'] = cf['pg_fn']
        rj['in_benchmark'] = 1 if core_frac >= 0.999999 else 0
        rj['notes'] = notes
        rj['coordinates'] = 'CHM13v2.0, 1-based inclusive; span = first base of the left anchor node to the last base of the right anchor node'
        rj['pad_bp'] = pad
        rj['span_bp'] = B1 - A1 + 1
        rj['core_bp'] = locus['core_end'] - locus['core_start'] + 1
        rj['flank_left_bp'] = locus['core_start'] - A1
        rj['flank_right_bp'] = B1 - locus['core_end']
        rj['anchor_left_seq'] = Lseq
        rj['anchor_right_seq'] = Rseq
        rj['sv_benchmark_frac_core'] = round(core_frac, 4)
        rj['sv_benchmark_frac_span'] = round(span_frac, 4)
        rj['tr_annotation'] = {
            'census_region_id': c.get('region_id'), 'tr_class': c.get('tr_class'),
            'period': cf['period'], 'motif': cf['motif'], 'copies': cf['copies'],
            'region_gc': _num(c.get('region_gc')), 'dist_telomere': _num(c.get('dist_telomere'), int),
            'matched_to': c.get('matched_to') or None, 'match': c.get('match') or None,
            'reference_tandem_scan': s.get('reference_tandem_scan'),
            'truth_h1_tandem_scan': s.get('truth_h1_tandem_scan'),
            'repeatmasker_tandem_like': (s.get('repeatmasker') or {}).get('tandem_like', [])[:10],
        }
        rj['census_counts'] = {k: _num(c.get(k), int) for k in
                               ('vg_fp', 'vg_fn', 'vg_tp_base', 'vg_tp_comp', 'pg_fp', 'pg_fn',
                                'truth_sv', 'truth_indel10') if k in c}
        rj['hap32'] = {
            'n_spanning_paths': len(span_recs), 'n_distinct_sequences': distinct,
            'n_fragments': len(frag_recs), 'anchor_only_pieces': anchor_only,
            'haplotypes_not_spanning': nonspan_haps,
            'len_min': min(lens), 'len_median': statistics.median(lens),
            'len_max': max(lens), 'chm13_len': len(refseq), 'grch38_spans': bool(grch38),
            'n_sequences_with_N_runs': n_with_N,
            'anchor_coverage_H': H,
        }
        rj['baseline_mc'] = {
            'graph': 'Minigraph-Cactus HPRC v2.1, hap32 sampled panel (32 haplotypes + CHM13 + GRCh38)',
            'nodes': len(nodes), 'edges': len(edges), 'node_bp': sum(nlen),
            'nodes_per_kb_chm13': round(len(nodes) / (len(refseq) / 1000.0), 2),
            'mean_node_len': round(sum(nlen) / len(nlen), 2), 'median_node_len': statistics.median(nlen),
            'between_subgraph_nodes': topo.get('nodes'), 'between_subgraph_edges': topo.get('edges'),
            'fragment_only_nodes': len(set(seqs) - nodes),
            'directed_cycle': (topo.get('cycles') or {}).get('graph_has_directed_cycle'),
            'census_graph_features': {k: _num(c.get(k)) for k in c if k.startswith('g_')},
        }
        rj['truth'] = s.get('truth')
        rj['closest_hap32_to_truth'] = s.get('closest')
        rj['called_haplotypes_vs_truth'] = s.get('called_haplotypes_vs_truth')
        rj['calls_in_span'] = {'n_rows': n_calls, 'called_sv50': call_counts,
                               'truth_sv50_status': s.get('truth_sv50_status')}
        rj['reads_census'] = {k: _num(c.get(k)) for k in ('reads_per_kb', 'mapq0_frac', 'mapq_lt5_frac')
                              if k in c}
        rj['hprc'] = {'status': 'not run'}
        rj['provenance'] = {
            'built': now(), 'script': 'tools/package_regions.py', 'region_py': 'tools/region.py',
            'graph_db': paths['gbz_db'], 'contig_gbz_for_names': paths['contig_gbz'],
            'truth': [paths['stvar'], paths['smvar']], 'vg_calls': paths['vg_vcf'],
            'pangenie_calls': paths['pg_norm'],
            'names_check': 'named runs from vg paths -A == gbz-base --between subpaths (multiset)',
            'n_paths_scanned': nr.get('n_paths_scanned'),
            'toolkit_warnings': s.get('warnings', []),
        }
        rj['files'] = file_table(rdir)
        with open(os.path.join(rdir, 'region.json'), 'w') as f:
            json.dump(rj, f, indent=1)
        return rid, 'ok', None
    except Exception as e:  # noqa: BLE001
        with open(os.path.join(work, rid, 'write_error.txt'), 'w') as f:
            f.write(traceback.format_exc())
        return rid, 'error', '%s: %s' % (type(e).__name__, str(e)[:300])


def file_table(rdir):
    out = collections.OrderedDict()
    for fn in ('hap32.fa', 'hap32.fragments.fa', 'hprc.fa.gz', 'truth.fa', 'mc.gfa', 'calls.tsv'):
        p = os.path.join(rdir, fn)
        if os.path.exists(p):
            e = {'bytes': os.path.getsize(p)}
            if fn.endswith(('.fa', '.fa.gz')):
                e['n_sequences'] = sum(1 for _ in read_fasta(p))
            out[fn] = e
    return out


# ----------------------------------------------------------------------------
# Stage D: full HPRC panel
#
# The full graph (config.data_paths()['hprc_gbz']) has the hap32 node ids, but the hap32 anchors
# do not separate it (`--between` leaks: other haplotypes add edges out of the span), so the span
# is taken from an `--interval` query on CHM13 with enough context that a haplotype's allele
# stays inside the subgraph, and each walk is cut at the two anchor nodes. gbz-base reports walks
# without names, and the query returns REFERENCE-sense paths too: CHM13, GRCh38, and the gref
# cover (gref_CHM13#0#<contig>, a copy of the CHM13 path, plus fragments that by construction
# never contain a reference node -- vg src/gref.hpp -- so never an anchor). Those are subtracted
# exactly: two copies of the CHM13 segment and one per spanning GRCh38 path of hap32.fa. What
# remains are the panel's sample haplotypes (456 in this graph; HG002 is not among them).

HPRC_SAMPLE_HAPLOTYPES = 456


def hprc_params(contig, L, R, a0, b0):
    span = b0 - a0
    return {'contig': contig, 'L': L, 'R': R, 'a0': a0, 'b0': b0,
            'context': max(10000, span)}


def hprc_params_from_summary(tdir):
    s = json.load(open(os.path.join(tdir, 'summary.json')))
    L = vr.parse_steps(s['boundaries']['L']['node'])[0]
    R = vr.parse_steps(s['boundaries']['R']['node'])[0]
    return hprc_params(s['contig'], L, R, s['boundaries']['L']['start0'], s['boundaries']['R']['end0'])


def hprc_fetch(rid, work, prm, force=False):
    """Run (or reuse) the full-graph interval query; returns (gfa_gz_path, meta)."""
    wdir = os.path.join(work, rid)
    os.makedirs(wdir, exist_ok=True)
    gz = os.path.join(wdir, 'hprc.interval.gfa.gz')
    mp = os.path.join(wdir, 'hprc.interval.json')
    if not force and os.path.exists(gz) and os.path.exists(mp):
        meta = json.load(open(mp))
        if all(meta.get(k) == prm[k] for k in prm):
            return gz, meta
    gbz = config.data_paths(prm['contig'])['hprc_gbz']
    cmd = [config.GBZ_BASE, 'query', gbz, '--sample', 'CHM13', '--contig', prm['contig'],
           '--interval', '%d..%d' % (prm['a0'], prm['b0']), '--context', str(prm['context']),
           '--distinct']
    t0 = time.time()
    tmp = gz + '.tmp%d' % os.getpid()      # two fetchers may race on one locus
    with open(os.path.join(wdir, 'hprc.log'), 'w') as fe:
        fe.write('$ %s\n' % ' '.join(cmd))
        fe.flush()
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=fe, env=config.tool_env())
        with gzip.open(tmp, 'wb', compresslevel=3) as fo:
            shutil.copyfileobj(p.stdout, fo, 1 << 20)
        p.wait()
    meta = dict(prm, cmd=' '.join(cmd), seconds=round(time.time() - t0, 1), returncode=p.returncode)
    if p.returncode != 0:
        meta['error'] = open(os.path.join(wdir, 'hprc.log')).read()[-400:].strip()
    os.replace(tmp, gz)
    with open(mp + '.tmp%d' % os.getpid(), 'w') as f:
        json.dump(meta, f, indent=1)
    os.replace(mp + '.tmp%d' % os.getpid(), mp)
    return gz, meta


def cut_at_anchors(g, L, R):
    """Cut every walk of the interval subgraph at the anchors. Returns ({steps: weight} of the
    L..R segments in reference orientation, Counter of the other anchor-touching walks)."""
    spans, other = collections.Counter(), collections.Counter()
    aL, aR = abs(L), abs(R)
    for w in g.walks:
        st, wt = w['steps'], w['weight']
        iL = [i for i, s in enumerate(st) if abs(s) == aL]
        iR = [i for i, s in enumerate(st) if abs(s) == aR]
        if not iL and not iR:
            continue
        if len(iL) == 1 and len(iR) == 1:
            i, j = iL[0], iR[0]
            if st[i] == L and st[j] == R and i < j:
                spans[tuple(st[i:j + 1])] += wt
            elif st[i] == -L and st[j] == -R and j < i:
                spans[vr.reverse_steps(st[j:i + 1])] += wt
            else:
                other['anchors_out_of_order'] += wt
        elif len(iL) > 1 or len(iR) > 1:
            other['anchor_visited_twice'] += wt
        elif iL:
            other['left_anchor_only'] += wt
        else:
            other['right_anchor_only'] += wt
    return spans, other


def hprc_one(args):
    rid, out, work, force = args
    rdir = os.path.join(out, rid)
    rjp = os.path.join(rdir, 'region.json')
    t0 = time.time()
    try:
        rj = json.load(open(rjp))
        fa = os.path.join(rdir, 'hprc.fa.gz')
        if not force and os.path.exists(fa) and rj.get('hprc', {}).get('status') == 'ok':
            return rid, 'cached', None, 0.0
        L = vr.parse_anchor(rj['anchor_left'])
        R = vr.parse_anchor(rj['anchor_right'])
        if rj['span_end'] - rj['span_start'] + 1 > HPRC_MAX_SPAN:
            info = {'status': 'skipped', 'reason': 'span %d bp exceeds %d bp: the full-graph query '
                    'would carry ~456 haplotypes of it' % (rj['span_end'] - rj['span_start'] + 1,
                                                           HPRC_MAX_SPAN)}
            return finish_hprc(rjp, rj, rdir, info, None, t0)
        prm = hprc_params(rj['contig'], L, R, rj['span_start'] - 1, rj['span_end'])
        gz, meta = hprc_fetch(rid, work, prm, force=False)
        method = ('gbz-base query --sample CHM13 --contig %s --interval %d..%d --context %d --distinct '
                  'on %s; walks cut at the anchors' % (prm['contig'], prm['a0'], prm['b0'],
                                                       prm['context'],
                                                       os.path.basename(config.data_paths()['hprc_gbz'])))
        if meta.get('returncode') != 0:
            info = {'status': 'failed', 'reason': 'gbz-base failed: ' + meta.get('error', '')[-200:],
                    'method': method}
            return finish_hprc(rjp, rj, rdir, info, None, t0)
        with gzip.open(gz, 'rt') as f:
            g = vr.Gfa(f.read())
        hap32 = read_fasta(os.path.join(rdir, 'hap32.fa'))
        chm13 = [(nm, sq) for nm, _, sq in hap32 if nm.startswith('CHM13#')]
        grch38 = [(nm, sq) for nm, _, sq in hap32 if nm.startswith('GRCh38#')]
        for node, want in ((L, rj['anchor_left_seq']), (R, rj['anchor_right_seq'])):
            if abs(node) not in g.seqs:
                info = {'status': 'failed', 'reason': 'anchor %d absent from the full graph' % abs(node),
                        'method': method}
                return finish_hprc(rjp, rj, rdir, info, None, t0)
            got = g.seqs[abs(node)] if node > 0 else revcomp(g.seqs[abs(node)])
            if got != want:
                info = {'status': 'failed', 'reason': 'anchor %d has another sequence in the full graph'
                        % abs(node), 'method': method}
                return finish_hprc(rjp, rj, rdir, info, None, t0)
        spans, other = cut_at_anchors(g, L, R)
        groups = collections.Counter()
        for st, wt in spans.items():
            groups[g.walk_seq(st)] += wt
        n_all = sum(groups.values())
        chm13_weight = groups[chm13[0][1]]
        # subtract the REFERENCE-sense paths: CHM13 and its gref copy, and GRCh38
        sub = []
        c_seq = chm13[0][1]
        if groups[c_seq] < 2:
            info = {'status': 'failed', 'method': method,
                    'reason': 'the CHM13 segment has weight %d in the full graph (expected >= 2: CHM13 '
                              'and gref_CHM13)' % groups[c_seq]}
            return finish_hprc(rjp, rj, rdir, info, None, t0)
        groups[c_seq] -= 2
        sub.append('CHM13 x2 (CHM13 and its gref_CHM13 copy)')
        for nm, sq in grch38:
            if groups[sq] < 1:
                info = {'status': 'failed', 'method': method,
                        'reason': '%s from hap32.fa has no copy among the full-graph segments' % nm}
                return finish_hprc(rjp, rj, rdir, info, None, t0)
            groups[sq] -= 1
            sub.append(nm)
        groups = +groups
        n_span = sum(groups.values())
        # more segments than haplotypes is possible: an assembly can carry two contigs through
        # the span (a false duplication), or one contig can pass it twice
        excess = max(0, n_span - HPRC_SAMPLE_HAPLOTYPES)
        by32 = collections.defaultdict(list)
        for nm, _, sq in hap32:
            by32[sq].append(nm)
        order = sorted(groups.items(), key=lambda kv: (-kv[1], len(kv[0]), kv[0]))
        recs = [(nm, 'reference path (also in hap32.fa) len=%d' % len(sq), sq) for nm, sq in chm13 + grch38]
        k = 0
        for gi, (sq, wt) in enumerate(order, start=1):
            tag = ''
            if sq == c_seq:
                tag = ' same_as=CHM13'
            elif sq in by32:
                tag = ' same_as_hap32=%s' % ','.join(by32[sq][:3]) + (',...' if len(by32[sq]) > 3 else '')
            for _ in range(wt):
                k += 1
                recs.append(('hprc#%d' % k, 'group=%d group_size=%d len=%d%s' % (gi, wt, len(sq), tag), sq))
        tmp = fa + '.tmp.gz'
        write_fasta(tmp, recs)
        os.replace(tmp, fa)
        lens = [len(s) for s in groups]
        info = {'status': 'ok', 'method': method,
                'graph': 'HPRC v2.1 Minigraph-Cactus, CHM13-based eval build with gref cover: 456 sample '
                         'haplotypes (HG002 not among them) + CHM13 + GRCh38 + gref_CHM13',
                'n_sample_haplotypes_spanning': n_span,
                'n_distinct_sample_sequences': len(groups),
                'n_distinct_in_hap32': sum(1 for s in groups if s in by32),
                'sample_haplotypes_in_hap32_sequences': sum(w for s, w in groups.items() if s in by32),
                'n_reference_records': len(chm13) + len(grch38),
                'subtracted_reference_segments': sub,
                'segments_before_subtraction': n_all,
                'segments_beyond_456_haplotypes': excess,
                'chm13_segment_weight_before_subtraction': chm13_weight,
                'anchor_touching_walks_not_spanning': dict(other),
                'len_min': min(lens) if lens else None,
                'len_max': max(lens) if lens else None,
                'subgraph_nodes': len(g.seqs), 'query_seconds': meta.get('seconds'),
                'naming': 'hprc#k: gbz-base reports walks without sample names; k is an ordinal over '
                          'the sample haplotypes, grouped by identical sequence (group, by decreasing '
                          'group size), so names are stable for this graph and anchor pair but carry '
                          'no sample identity. CHM13 and GRCh38 keep their hap32.fa names.'}
        return finish_hprc(rjp, rj, rdir, info, n_span, t0)
    except Exception as e:  # noqa: BLE001
        with open(os.path.join(work, rid, 'hprc_error.txt'), 'w') as f:
            f.write(traceback.format_exc())
        return rid, 'error', '%s: %s' % (type(e).__name__, str(e)[:300]), time.time() - t0


def finish_hprc(rjp, rj, rdir, info, n_span, t0):
    info['seconds'] = round(time.time() - t0, 1)
    rj = json.load(open(rjp))     # re-read: another stage may have touched it
    rj['hprc'] = info
    rj['n_hprc'] = n_span or 0
    rj['notes'] = [n for n in rj['notes'] if not n.startswith('hprc:')]
    if info['status'] != 'ok':
        rj['notes'].append('hprc: ' + info['reason'])
        fa = os.path.join(rdir, 'hprc.fa.gz')
        if os.path.exists(fa):
            os.remove(fa)
    else:
        if info.get('segments_beyond_456_haplotypes'):
            rj['notes'].append('hprc: %d more spanning segments than the 456 sample haplotypes (an '
                               'assembly carries two contigs through the span, or one passes it twice)'
                               % info['segments_beyond_456_haplotypes'])
        nn = sum(info['anchor_touching_walks_not_spanning'].values())
        if nn:
            rj['notes'].append('hprc: %d full-panel walks touch an anchor without spanning (includes '
                               'reference paths)' % nn)
    rj['files'] = file_table(rdir)
    tmp = rjp + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(rj, f, indent=1)
    os.replace(tmp, rjp)
    return rj['region_id'], info['status'], info.get('reason'), time.time() - t0


def hprc_priority(rid, stratum):
    return (0 if rid in PILOT else 1 if stratum == 'hotspot_vntr' else
            2 if stratum == 'control_vntr_matched' else 3, rid)


def cmd_hprc_fetch(a):
    """Query the full graph for every locus of --loci as soon as stage A has written its
    anchors (so that this can run next to `build`), pilots and hotspots first."""
    loci = load_loci(a.loci)
    if a.only:
        keep = set(a.only.split(','))
        loci = [l for l in loci if l['region_id'] in keep]
    if not a.keep_order:
        loci.sort(key=lambda l: hprc_priority(l['region_id'], l['stratum']))
    print('[%s] hprc-fetch: %d loci' % (now(), len(loci)), flush=True)
    for l in loci:
        rid = l['region_id']
        tdir = os.path.join(a.work, rid, 'toolkit')
        waited = 0
        while not os.path.exists(os.path.join(tdir, 'summary.json')):
            if os.path.exists(os.path.join(a.work, rid, 'error.txt')) or waited > a.wait:
                break
            time.sleep(10)
            waited += 10
        if not os.path.exists(os.path.join(tdir, 'summary.json')):
            print('  %-10s skipped (no anchors)' % rid, flush=True)
            continue
        prm = hprc_params_from_summary(tdir)
        if prm['b0'] - prm['a0'] > HPRC_MAX_SPAN:
            print('  %-10s skipped (span %d bp)' % (rid, prm['b0'] - prm['a0']), flush=True)
            continue
        gz, meta = hprc_fetch(rid, a.work, prm, force=a.force)
        print('  %-10s %s rc=%s %ss %.1f MB' % (rid, prm['contig'], meta.get('returncode'),
                                              meta.get('seconds'), os.path.getsize(gz) / 1e6), flush=True)
    print('[%s] hprc-fetch done' % now(), flush=True)


# ----------------------------------------------------------------------------
# Index, BED, pruning, validation

def region_dirs(out):
    return sorted(d for d in os.listdir(out) if os.path.exists(os.path.join(out, d, 'region.json')))


def fmt(v):
    if v is None:
        return '.'
    if isinstance(v, float):
        return ('%.1f' % v) if v != int(v) else str(int(v))
    if isinstance(v, list):
        return '; '.join(str(x) for x in v) if v else '.'
    return str(v)


def write_index(out, order=None):
    rows = [json.load(open(os.path.join(out, d, 'region.json'))) for d in region_dirs(out)]
    rank = {s: i for i, s in enumerate(['hotspot_vntr', 'control_vntr_matched', 'control_vntr_correct',
                                        'hotspot_other', 'control_nontr_sv', 'target'])}
    pos = {rid: i for i, rid in enumerate(order or [])}
    rows.sort(key=lambda r: (rank.get(r['stratum'], 9), pos.get(r['region_id'], 0), r['region_id']))
    with open(os.path.join(out, 'regions.tsv'), 'w') as f:
        f.write('\t'.join(INDEX_COLS) + '\n')
        for r in rows:
            f.write('\t'.join(fmt(r.get(c)).replace('\t', ' ') for c in INDEX_COLS) + '\n')
    return len(rows)


def write_bed(vntr_regions, out, min_len):
    n = 0
    bp = 0
    path = os.path.join(out, 'vntr_ge%dkb.bed' % (min_len // 1000) if min_len % 1000 == 0
                        else 'vntr_ge%dbp.bed' % min_len)
    rows = []
    with open(vntr_regions) as f:
        h = f.readline().rstrip('\n').split('\t')
        for line in f:
            r = dict(zip(h, line.rstrip('\n').split('\t')))
            if r.get('in_bench', '1') not in ('1', '1.0', 'True'):
                continue
            if r.get('tr_class', 'VNTR') != 'VNTR':
                continue
            s, e = int(r['start']), int(r['end'])
            if e - s < min_len:
                continue
            rows.append((r['chrom'], s, e, r['region_id'], r['period'], r['motif'], r['vg_fp'], r['vg_fn']))
    order = {('chr%d' % i): i for i in range(1, 23)}
    order.update({'chrX': 23, 'chrY': 24, 'chrM': 25})
    rows.sort(key=lambda x: (order.get(x[0], 99), x[1]))
    with open(path, 'w') as f:
        f.write('#contig\tstart\tend\tregion_id\tperiod\tmotif\tvg_fp\tvg_fn\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')
            n += 1
            bp += r[2] - r[1]
    return path, n, bp


def prune_hprc(out, max_mb, keep_strata):
    dirs = region_dirs(out)
    sizes = {d: os.path.getsize(os.path.join(out, d, 'hprc.fa.gz'))
             for d in dirs if os.path.exists(os.path.join(out, d, 'hprc.fa.gz'))}
    total = sum(sizes.values())
    removed = []
    if total > max_mb * 1e6:
        for d in dirs:
            rj = json.load(open(os.path.join(out, d, 'region.json')))
            if d in sizes and rj['stratum'] not in keep_strata and d not in PILOT:
                os.remove(os.path.join(out, d, 'hprc.fa.gz'))
                rj['hprc']['file_removed'] = 'hprc.fa.gz dropped to keep the package under %d MB' % max_mb
                rj['files'] = file_table(os.path.join(out, d))
                json.dump(rj, open(os.path.join(out, d, 'region.json'), 'w'), indent=1)
                removed.append(d)
    return total, removed


def validate(out):
    rep = collections.OrderedDict()
    dirs = region_dirs(out)
    per_stratum = collections.Counter()
    missing, bad_paths, bad_flanks, bad_names = [], [], [], []
    total_bytes = collections.Counter()
    n_seq = collections.Counter()
    for d in dirs:
        rdir = os.path.join(out, d)
        rj = json.load(open(os.path.join(rdir, 'region.json')))
        per_stratum[rj['stratum']] += 1
        for fn in ('hap32.fa', 'hap32.fragments.fa', 'mc.gfa', 'truth.fa', 'calls.tsv', 'region.json'):
            if not os.path.exists(os.path.join(rdir, fn)):
                missing.append('%s/%s' % (d, fn))
        for fn in os.listdir(rdir):
            total_bytes[fn] += os.path.getsize(os.path.join(rdir, fn))
        h = read_fasta(os.path.join(rdir, 'hap32.fa'))
        names = [x[0] for x in h]
        if len(names) != len(set(names)) or len(h) != rj['n_hap32']:
            bad_names.append(d)
        n_seq['hap32'] += len(h)
        want = {x[0]: x[2] for x in h}
        _, probs = gfa_spell_check(os.path.join(rdir, 'mc.gfa'), want)
        if probs:
            bad_paths.append((d, probs[:3]))
        Ls, Rs = rj['anchor_left_seq'], rj['anchor_right_seq']
        t = read_fasta(os.path.join(rdir, 'truth.fa'))
        allseq = [x[2] for x in h] + [x[2] for x in t]
        hp = os.path.join(rdir, 'hprc.fa.gz')
        if os.path.exists(hp):
            hh = read_fasta(hp)
            n_seq['hprc'] += len(hh)
            allseq += [x[2] for x in hh]
            nref = sum(1 for x in hh if not x[0].startswith('hprc#'))
            if len(hh) - nref != rj['n_hprc'] or len(set(x[0] for x in hh)) != len(hh):
                bad_names.append(d + ':hprc_count')
        if not all(s.startswith(Ls) and s.endswith(Rs) for s in allseq):
            bad_flanks.append(d)
        chm = [x[2] for x in h if x[0].startswith('CHM13#')]
        if len(chm) != 1 or len(chm[0]) != rj['span_end'] - rj['span_start'] + 1:
            bad_flanks.append(d + ':chm13_len')
    rep['regions'] = len(dirs)
    rep['per_stratum'] = dict(per_stratum)
    rep['missing_files'] = missing
    rep['mc_gfa_path_spelling_failures'] = bad_paths
    rep['flank_or_length_failures'] = bad_flanks
    rep['name_or_count_failures'] = bad_names
    rep['sequences'] = dict(n_seq)
    rep['bytes_by_file'] = dict(total_bytes)
    rep['total_mb'] = round(sum(total_bytes.values()) / 1e6, 2)
    return rep


# ----------------------------------------------------------------------------
# Drivers

def run_pool(fn, jobs, items):
    if jobs <= 1:
        for it in items:
            yield fn(it)
        return
    with multiprocessing.Pool(jobs) as pool:
        for res in pool.imap_unordered(fn, items):
            yield res


def snarl_assign(loci, pad, max_span, cache):
    """Snarl-mode anchors for every locus; overlapping spans are merged (repeatedly, since the
    union's anchors can reach further). Returns (loci to build, {region_id: skip reason})."""
    import snarl_anchors as sa
    by_contig = collections.defaultdict(list)
    for l in loci:
        by_contig[l['contig']].append(l)
    out, skipped = [], {}
    for contig, ls in by_contig.items():
        A = sa.Anchors(contig, cache)
        clen = vr.contig_length(config.data_paths(contig)['ref_fa'], contig)
        # at a contig end there is no boundary node beyond the interval: clamp the padded
        # interval to the first / last boundary node on CHM13
        first_end = min(p[1] for p in A.pts)
        last_start = max(p[0] for p in A.pts)

        def anchors(a0, b0):
            return A.query(max(a0, first_end), min(b0, last_start), max_span)

        groups = []
        for l in ls:
            a0, b0 = max(0, l['core_start'] - 1 - pad), min(clen, l['core_end'] + pad)
            try:
                L, R, info = anchors(a0, b0)
            except sa.SnarlAnchorError as e:
                skipped[l['region_id']] = str(e)
                continue
            groups.append({'loci': [l], 'a0': a0, 'b0': b0, 'L': L, 'R': R, 'info': info})
        while True:
            groups.sort(key=lambda g: g['info']['L'][0])
            merged, changed = [], False
            for g in groups:
                if merged and g['info']['L'][0] < merged[-1]['info']['R'][1]:
                    m = merged[-1]
                    m['loci'] += g['loci']
                    m['a0'], m['b0'] = min(m['a0'], g['a0']), max(m['b0'], g['b0'])
                    changed = True
                    try:
                        m['L'], m['R'], m['info'] = anchors(m['a0'], m['b0'])
                    except sa.SnarlAnchorError as e:
                        m['error'] = str(e)
                else:
                    merged.append(g)
            groups = merged
            if not changed:
                break
        for g in groups:
            ids = [l['region_id'] for l in g['loci']]
            if g.get('error'):
                for rid in ids:
                    skipped[rid] = 'merged group %s: %s' % (','.join(ids), g['error'])
                continue
            l = dict(g['loci'][0])
            if len(ids) > 1:
                l['core_start'] = min(x['core_start'] for x in g['loci'])
                l['core_end'] = max(x['core_end'] for x in g['loci'])
                l['merged_loci'] = ids
            l['anchors'] = (g['L'], g['R'])
            l['anchor_mode'] = 'snarl'
            l['anchor_note'] = 'snarl anchors: %d+ / %d+, chain span %d bp (%s)%s' % (
                g['L'], g['R'], g['info']['span_bp'], g['info']['rule'],
                ('; merged overlapping loci ' + ','.join(ids)) if len(ids) > 1 else '')
            out.append(l)
    return out, skipped


def cmd_build(a):
    loci = load_loci(a.loci)
    if a.only:
        keep = set(a.only.split(','))
        loci = [l for l in loci if l['region_id'] in keep]
    snarl_skipped = {}
    if a.anchor_mode == 'snarl':
        n0 = len(loci)
        loci, snarl_skipped = snarl_assign(loci, a.pad, a.max_span, a.snarl_cache)
        print('[%s] snarl anchors: %d loci -> %d regions (%d merged away, %d skipped)' % (
            now(), n0, len(loci), n0 - len(loci) - len(snarl_skipped), len(snarl_skipped)), flush=True)
        for k, v in snarl_skipped.items():
            print('  SKIP %s: %s' % (k, v), flush=True)
    os.makedirs(a.work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    for l in loci:
        os.makedirs(os.path.join(a.work, l['region_id']), exist_ok=True)
    t0 = time.time()
    status = {k: ('skipped', v) for k, v in snarl_skipped.items()}
    print('[%s] stage A: extract %d loci (jobs %d, pad %d)' % (now(), len(loci), a.jobs, a.pad), flush=True)
    for rid, st, msg, sec in run_pool(stage_extract, a.jobs, [(l, a.work, a.pad, a.force, a.max_span) for l in loci]):
        status[rid] = (st, msg)
        print('  %-10s %-6s %5.1fs %s' % (rid, st, sec, msg or ''), flush=True)
    ok = [l for l in loci if status[l['region_id']][0] in ('ok', 'cached')]
    # stage B: names, one vg paths -A pass per contig
    by_contig = collections.defaultdict(list)
    for l in ok:
        rid = l['region_id']
        if not a.force and os.path.exists(os.path.join(a.work, rid, 'named_runs.json')) and \
                status[rid][0] == 'cached':
            continue
        tdir = os.path.join(a.work, rid, 'toolkit')
        s = json.load(open(os.path.join(tdir, 'summary.json')))
        g = vr.Gfa(open(os.path.join(tdir, 'subgraph.query.gfa')).read())
        by_contig[s['contig']].append({'id': rid, 'nodes': set(g.seqs),
                                       'L': vr.parse_steps(s['boundaries']['L']['node'])[0],
                                       'R': vr.parse_steps(s['boundaries']['R']['node'])[0]})
    print('[%s] stage B: name paths on %d contigs' % (now(), len(by_contig)), flush=True)
    # largest contigs first so the pool finishes together
    items = sorted(((c, infos, a.work) for c, infos in by_contig.items()),
                   key=lambda x: -os.path.getsize(config.data_paths(x[0])['contig_gbz']))
    for contig, n_paths, err, sec in run_pool(stage_names, a.jobs, items):
        print('  %-6s %s %6.1fs %s' % (contig, n_paths, sec, err or ''), flush=True)
        if err:
            for inf in by_contig[contig]:
                status[inf['id']] = ('error', 'names: ' + err)
    ok = [l for l in ok if status[l['region_id']][0] in ('ok', 'cached')]
    print('[%s] stage C: write %d packages' % (now(), len(ok)), flush=True)
    for rid, st, msg in run_pool(stage_write, a.jobs, [(l, a.work, a.out, a.pad) for l in ok]):
        if st != 'ok':
            status[rid] = (st, 'write: ' + (msg or ''))
            print('  %-10s %s %s' % (rid, st, msg), flush=True)
    n = write_index(a.out, [l['region_id'] for l in load_loci(a.loci)])
    failed = {k: v for k, v in status.items() if v[0] == 'error'}
    with open(os.path.join(a.work, 'build_status.json'), 'w') as f:
        json.dump({'finished': now(), 'seconds': round(time.time() - t0, 1), 'status': status}, f, indent=1)
    print('[%s] done in %.0f s: %d packages indexed, %d failed' % (now(), time.time() - t0, n, len(failed)))
    for k, v in failed.items():
        print('  FAILED %s: %s' % (k, v[1]))


def cmd_hprc(a):
    dirs = region_dirs(a.out)
    if a.only:
        keep = set(a.only.split(','))
        dirs = [d for d in dirs if d in keep]
    # hotspots and pilots first, so that a partial run has the loci that matter most
    dirs.sort(key=lambda d: hprc_priority(d, json.load(open(os.path.join(a.out, d, 'region.json')))['stratum']))
    print('[%s] hprc: %d regions, jobs %d' % (now(), len(dirs), a.jobs), flush=True)
    t0 = time.time()
    for rid, st, msg, sec in run_pool(hprc_one, a.jobs, [(d, a.out, a.work, a.force) for d in dirs]):
        print('  %-10s %-7s %5.1fs %s' % (rid, st, sec, msg or ''), flush=True)
    write_index(a.out)
    print('[%s] hprc done in %.0f s' % (now(), time.time() - t0), flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('build')
    b.add_argument('--loci', default=config.census_file('strata.tsv'))
    b.add_argument('--out', default=config.REGIONS_DIR)
    b.add_argument('--work', default=default_work())
    b.add_argument('--pad', type=int, default=200)
    b.add_argument('--max-span', type=int, default=MAX_SPAN,
                   help='longest anchor-to-anchor span before falling back to nearer anchors')
    b.add_argument('--jobs', type=int, default=3)
    b.add_argument('--anchor-mode', choices=('truth', 'snarl'), default='truth',
                   help='truth: CHM13 nodes every hap32 path visits once and no truth record touches '
                        '(with a nearest-best-covered fallback); snarl: truth-free snarl chain boundaries')
    b.add_argument('--snarl-cache', default=os.path.join(config.WORK_DIR, 'stage4', 'snarls'),
                   help='snarl_anchors.py build output directory (snarl mode)')
    b.add_argument('--only')
    b.add_argument('--force', action='store_true')
    h = sub.add_parser('hprc')
    h.add_argument('--out', default=config.REGIONS_DIR)
    h.add_argument('--work', default=default_work())
    h.add_argument('--jobs', type=int, default=1)
    h.add_argument('--only')
    h.add_argument('--force', action='store_true')
    hf = sub.add_parser('hprc-fetch')
    hf.add_argument('--loci', default=config.census_file('strata.tsv'))
    hf.add_argument('--work', default=default_work())
    hf.add_argument('--only')
    hf.add_argument('--wait', type=int, default=7200, help='seconds to wait for a locus\'s anchors')
    hf.add_argument('--force', action='store_true')
    hf.add_argument('--keep-order', action='store_true',
                    help='fetch in the order of --loci (default: pilots, hotspots, matched controls, rest)')
    i = sub.add_parser('index')
    i.add_argument('--out', default=config.REGIONS_DIR)
    i.add_argument('--loci')
    p = sub.add_parser('prune-hprc')
    p.add_argument('--out', default=config.REGIONS_DIR)
    p.add_argument('--max-mb', type=float, default=150)
    p.add_argument('--keep-strata', default='hotspot_vntr')
    e = sub.add_parser('bed')
    e.add_argument('--vntr-regions', default=config.census_file('vntr_regions.tsv'))
    e.add_argument('--out', default=config.REGIONS_DIR)
    e.add_argument('--min-len', type=int, default=1000)
    v = sub.add_parser('validate')
    v.add_argument('--out', default=config.REGIONS_DIR)
    a = ap.parse_args(argv)
    if a.cmd == 'build':
        if not a.loci:
            raise SystemExit('no locus list: pass --loci (census strata.tsv, a regions.tsv, or a BED)')
        cmd_build(a)
    elif a.cmd == 'hprc':
        cmd_hprc(a)
    elif a.cmd == 'hprc-fetch':
        cmd_hprc_fetch(a)
    elif a.cmd == 'index':
        n = write_index(a.out, [l['region_id'] for l in load_loci(a.loci)] if a.loci else None)
        print('%d regions indexed in %s' % (n, os.path.join(a.out, 'regions.tsv')))
    elif a.cmd == 'prune-hprc':
        total, removed = prune_hprc(a.out, a.max_mb, set(a.keep_strata.split(',')))
        print('hprc.fa.gz total %.1f MB; removed %d' % (total / 1e6, len(removed)))
    elif a.cmd == 'bed':
        if not a.vntr_regions:
            raise SystemExit('pass --vntr-regions (the census vntr_regions.tsv)')
        path, n, bp = write_bed(a.vntr_regions, a.out, a.min_len)
        print('%s: %d regions, %.2f Mb' % (path, n, bp / 1e6))
    elif a.cmd == 'validate':
        print(json.dumps(validate(a.out), indent=1))


if __name__ == '__main__':
    main()
