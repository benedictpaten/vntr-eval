#!/usr/bin/env python3
"""realign_poa.py -- partial-order alignment (abPOA, spoa) realigners, plugged into realign.py.

Region mode: align regions/<id>/hap32.fa with a POA tool, build the graph with msa_graph.py and
write candidates/<METHOD>/<id>.msa.fa, <id>.gfa and <id>.realign.json, plus one row per region in
results/realign_runtime.tsv (all through realign.py's shared pipeline, the same as the mafft
methods):

    python3 tools/realign_poa.py regions/L012184 --method poa_abpoa [--out candidates]
    python3 tools/realign_poa.py L012184 L005990 --method poa_spoa
    python3 tools/realign_poa.py all --method poa_abpoa --jobs 2 [--stratum hotspot_vntr] [--force]

MSA-only mode, for an ARBITRARY input FASTA (e.g. the deduplicated full HPRC panel); writes only
the MSA (rows named by the first word of each header, in input order, '-' gaps, every row spelling
its input sequence; identical sequences are aligned once and their rows copied):

    python3 tools/realign_poa.py --method poa_abpoa --input seqs.fa[.gz] --msa-out out.msa.fa [--json info.json]
        [--project-to regions/<id>/hap32.fa --projected-out hap32rows.msa.fa]

--project-to keeps, for each hap32.fa record, a row with the identical sequence (full-panel
records are named hprc#k, not by hap32 path, so rows are matched by sequence), renames it, drops
all-gap columns and writes --projected-out, ready for
`python3 tools/msa_graph.py hap32rows.msa.fa regions/<id>/hap32.fa out.gfa`.

    import sys; sys.path.insert(0, 'tools')
    import realign_poa, realign
    info = realign.align_fasta('poa_spoa', 'seqs.fa', 'out.msa.fa', timeout=900, mem_mb=12000)

Because this module defines METHODS, realign.py knows the two methods too
(`python3 tools/realign.py poa_abpoa all`, `realign.py --list`).

Methods (both always GLOBAL: every sequence runs anchor to anchor, so a local or semi-global
mode would leave the shared flanks unaligned; see tools/POA.md for the pilot that chose them):

  poa_abpoa   abpoa -m 0 -r 1     (default scores 2/4, convex gaps O=4,24 E=2,1, adaptive band
              b=10 f=0.01), sequences added longest first (what abpoa -L does; the sort is done
              here, ties in input order)
  poa_spoa    spoa -l 1 -r 1      (default scores 5/-4, convex gaps g=-8 e=-6 q=-10 c=-4;
              unbanded), sequences added longest first

Variants for pilots (--config, --order; the output method name then gets a suffix unless --name):
  --config (abPOA):  default (no extra flags), p (-p), S (-S), Sp (-S -p), nb (-b -1, unbanded),
                     or --tool-args '...' for anything else
  --order:           given (input order; CHM13 first in hap32.fa), longest (length descending),
                     shortest, guide (nearest-first by a k-mer count distance: start at the medoid,
                     then always add the sequence closest to one already added, i.e. Prim's order
                     on a minimum spanning tree), random:SEED

What the shared pipeline (realign.align_fasta) does around the POA call: identical sequences are
aligned once; runs of non-ACGT characters (N) are cut out before alignment and put back as
insertion columns of their own after the preceding base, so msa_graph.py makes each run an N node
(the same for every realigner); the aligner runs under a wall-clock timeout (default 900 s) and a
cap on its RSS (killed at --mem-mb); its output is checked row by row.

Memory: spoa is unbanded and preallocates its DP matrices for 4 x max_len graph nodes x max_len
(80 B x max_len^2 measured); abPOA without -S allocates (graph nodes) x (sequence length), 18-29 B
x max_len^2 measured. Before running, the aligner's memory is predicted from the longest sequence
(MEM_C below; POA.md has the measurements); a region predicted over --mem-mb is not run and gets
status 'memout' with 'predicted' in its message. Regions that pass the prediction are still
killed when their sampled RSS passes --mem-mb.
"""
import argparse
import collections
import concurrent.futures
import json
import os
import random
import shutil
import subprocess
import sys
import time

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402
import msa_graph  # noqa: E402

DEFAULT_TIMEOUT = 900
DEFAULT_MEM_MB = 12000


def _tool(env, name):
    v = os.environ.get(env)
    if v:
        return os.path.expanduser(v)
    return shutil.which(name, path=config.TOOL_PATH) or name


ABPOA = _tool('VNTR_ABPOA', 'abpoa')
SPOA = _tool('VNTR_SPOA', 'spoa')

ABPOA_CONFIGS = collections.OrderedDict([
    ('default', []),
    ('p', ['-p']),
    ('S', ['-S']),
    ('Sp', ['-S', '-p']),
    ('nb', ['-b', '-1']),
])

# Predicted peak RSS in bytes = C * max_len^2 (+ 50 MB), measured with /usr/bin/time -l (POA.md).
# spoa preallocates its DP matrices for 4 x max_len graph nodes x max_len columns, five int32
# matrices under convex gaps: 80 B x max_len^2, observed 80-81 at four regions of 3-12 kb.
# abPOA (no -S) allocates (graph nodes) x (length): 18-29 B x max_len^2 observed at 6-20 kb,
# lower with -p; -S is irregular (up to 90). The abPOA figure is only used to refuse regions
# that cannot fit; the sampled-RSS kill is the real guard.
MEM_C = {'spoa': 80.0, 'abpoa': 22.0, 'abpoa_long': 24.0, 'abpoa_S': 24.0}

# the chosen configurations (POA.md explains the pilot)
DEFAULTS = {
    'poa_abpoa': {'config': 'default', 'order': 'longest'},
    'poa_spoa': {'config': None, 'order': 'longest'},
}

_VERSIONS = {}


def tool_version(tool):
    if tool not in _VERSIONS:
        exe = ABPOA if tool == 'abpoa' else SPOA
        flag = '-v' if tool == 'abpoa' else '--version'
        try:
            p = subprocess.run([exe, flag], capture_output=True, text=True, timeout=30)
            _VERSIONS[tool] = '%s %s' % (tool, (p.stdout + p.stderr).strip().splitlines()[0])
        except (OSError, subprocess.SubprocessError, IndexError):
            _VERSIONS[tool] = tool + ' unknown'
    return _VERSIONS[tool]


# ---------------------------------------------------------------- input order

_CODE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def _kmer_profile(s, k=16, scale=1):
    """Counts of the k-mers of s (2-bit codes; k-mers with a non-ACGT base are skipped), keeping
    only those whose mixed code is divisible by scale (FracMinHash; deterministic)."""
    c = collections.Counter()
    mask = (1 << (2 * k)) - 1
    x, n = 0, 0
    for ch in s:
        v = _CODE.get(ch)
        if v is None:
            x, n = 0, 0
            continue
        x = ((x << 2) | v) & mask
        n += 1
        if n >= k and (scale == 1 or ((x * 0x9E3779B97F4A7C15) >> 17) % scale == 0):
            c[x] += 1
    return c


def kmer_distance(a, b):
    """1 - weighted Jaccard (sum of min counts / sum of max counts); copy number counts."""
    if not a and not b:
        return 0.0
    num = den = 0
    for h in set(a) | set(b):
        x, y = a.get(h, 0), b.get(h, 0)
        num += min(x, y)
        den += max(x, y)
    return 1.0 - num / float(den) if den else 1.0


def guide_order(seqs, k=16):
    """Nearest-first order: the medoid first, then repeatedly the sequence closest to any one
    already chosen (Prim's order on the k-mer distance). Ties: input order."""
    n = len(seqs)
    if n <= 2:
        return list(range(n))
    scale = 1 if n <= 64 else 8
    prof = [_kmer_profile(s, k, scale) for s in seqs]
    d = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d[i][j] = d[j][i] = kmer_distance(prof[i], prof[j])
    start = min(range(n), key=lambda i: (sum(d[i]), i))
    order, best = [start], list(d[start])
    left = set(range(n)) - {start}
    while left:
        nxt = min(left, key=lambda j: (best[j], j))
        order.append(nxt)
        left.discard(nxt)
        for j in left:
            if d[nxt][j] < best[j]:
                best[j] = d[nxt][j]
    return order


def order_indices(seqs, order):
    n = len(seqs)
    if order in (None, '', 'given'):
        return list(range(n))
    if order == 'longest':
        return sorted(range(n), key=lambda i: (-len(seqs[i]), i))
    if order == 'shortest':
        return sorted(range(n), key=lambda i: (len(seqs[i]), i))
    if order == 'guide':
        return guide_order(seqs)
    if order.startswith('random:'):
        idx = list(range(n))
        random.Random(int(order.split(':', 1)[1])).shuffle(idx)
        return idx
    raise ValueError('unknown order %r (given, longest, shortest, guide, random:SEED)' % order)


# ---------------------------------------------------------------- aligners (realign.py plugin API)

def _reordered_input(in_fa, order, workdir):
    recs = msa_graph.read_fasta(in_fa)
    idx = order_indices([s for _, s in recs], order)
    path = os.path.join(workdir, 'ordered.fa')
    with open(path, 'w') as f:
        for i in idx:
            f.write('>%s\n%s\n' % recs[i])
    return path, recs, [recs[i][0] for i in idx]


def predict_mb(tool, lens, args=()):
    m = max(lens) if lens else 0
    if tool == 'spoa':
        c = MEM_C['spoa']
    elif '-S' in args:
        c = MEM_C['abpoa_S']
    else:
        c = MEM_C['abpoa_long'] if m > 16000 else MEM_C['abpoa']
    return round(50 + c * m * m / 2 ** 20, 1)


def available_mb():
    """Memory the OS could hand out now (free + inactive + speculative + purgeable pages), from
    vm_stat on macOS or /proc/meminfo elsewhere; None when unknown."""
    try:
        if os.path.exists('/proc/meminfo'):
            for line in open('/proc/meminfo'):
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1]) / 1024.0
            return None
        out = subprocess.run(['vm_stat'], capture_output=True, text=True, timeout=20).stdout
        page = 4096
        tot = 0
        for line in out.splitlines():
            if 'page size of' in line:
                page = int(line.split('page size of')[1].split()[0])
            for key in ('Pages free', 'Pages inactive', 'Pages speculative', 'Pages purgeable'):
                if line.startswith(key + ':'):
                    tot += int(line.split(':')[1].strip().rstrip('.'))
        return tot * page / 2.0 ** 20
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


def wait_for_memory(need_mb, max_wait=1800, poll=15, log=None):
    """Block until available_mb() >= need_mb (or max_wait seconds pass). Returns seconds waited."""
    t0 = time.time()
    while True:
        av = available_mb()
        if av is None or av >= need_mb or time.time() - t0 > max_wait:
            return round(time.time() - t0, 1)
        if log:
            log('waiting for %.0f MB (available %.0f MB)' % (need_mb, av))
        time.sleep(poll)


def _run(tool, cmd, in_fa, out_fa, workdir, timeout, mem_mb, order, extra):
    import realign                      # late import: realign imports this module as a plugin
    ordered, recs, names = _reordered_input(in_fa, order, workdir)
    pred = predict_mb(tool, [len(s) for _, s in recs], cmd)
    base = {'tool_version': tool_version(tool), 'order': order or 'given', 'input_order': names[:500],
            'predicted_mb': pred, 'n_aligned': len(recs), 'command': [tool] + cmd[1:] + ['<in.fa>']}
    base.update(extra)
    if mem_mb and pred > mem_mb:
        base.update(status='memout', returncode=None, seconds=0.0, peak_rss_mb=0.0,
                    message='not run: predicted %.0f MB > cap %.0f MB (max_len %d)' % (
                        pred, mem_mb, max(len(s) for _, s in recs)))
        return base
    log = os.path.join(workdir, tool + '.log')
    part = out_fa + '.part'
    r = realign.run_proc(cmd + [ordered], part, log, timeout=timeout, mem_mb=mem_mb,
                         env=config.tool_env(), cwd=workdir)
    base.update(r)
    if r['status'] == 'ok' and os.path.getsize(part) == 0:
        base['status'] = 'error'
    if base['status'] == 'ok':
        os.replace(part, out_fa)
    else:
        base['message'] = '%s after %.0f s (peak sampled RSS %.0f MB)%s' % (
            base['status'], base.get('seconds') or 0, base.get('peak_rss_mb') or 0,
            ('; log: ' + realign._tail(log, 3)) if realign._tail(log, 3) else '')
        if os.path.exists(part):
            os.remove(part)
    return base


def abpoa_align(in_fa, out_fa, threads=1, workdir=None, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB,
                config_name='default', order='longest', tool_args=None, **_):
    """abpoa -m 0 -r 1 [config flags] on in_fa reordered by `order`; abPOA is single-threaded."""
    flags = list(tool_args) if tool_args is not None else list(ABPOA_CONFIGS[config_name])
    cmd = [ABPOA, '-m', '0', '-r', '1'] + flags
    return _run('abpoa', cmd, in_fa, out_fa, workdir, timeout, mem_mb, order,
                {'config': config_name if tool_args is None else 'custom', 'flags': flags})


def spoa_align(in_fa, out_fa, threads=1, workdir=None, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB,
               order='longest', tool_args=None, **_):
    """spoa -l 1 -r 1 (global, default scores) on in_fa reordered by `order`; single-threaded."""
    flags = list(tool_args or [])
    cmd = [SPOA, '-l', '1', '-r', '1'] + flags
    return _run('spoa', cmd, in_fa, out_fa, workdir, timeout, mem_mb, order, {'flags': flags})


def _spec(tool, config_name=None, order=None, tool_args=None):
    if tool == 'abpoa':
        params = {'config_name': config_name or 'default', 'order': order or 'longest'}
        if tool_args is not None:
            params['tool_args'] = list(tool_args)
        flags = tool_args if tool_args is not None else ABPOA_CONFIGS[params['config_name']]
        desc = 'abPOA global (-m 0) %s, order %s' % (' '.join(flags) or 'default flags', params['order'])
        return {'align': abpoa_align, 'params': params, 'description': desc, 'tool': 'abpoa'}
    params = {'order': order or 'longest'}
    if tool_args:
        params['tool_args'] = list(tool_args)
    desc = 'spoa global (-l 1) %sorder %s' % ((' '.join(tool_args) + ', ') if tool_args else '', params['order'])
    return {'align': spoa_align, 'params': params, 'description': desc, 'tool': 'spoa'}


METHODS = collections.OrderedDict([
    ('poa_abpoa', _spec('abpoa', DEFAULTS['poa_abpoa']['config'], DEFAULTS['poa_abpoa']['order'])),
    ('poa_spoa', _spec('spoa', None, DEFAULTS['poa_spoa']['order'])),
])


def variant(method, config_name=None, order=None, tool_args=None, name=None):
    """(name, spec) for a method with a non-default configuration or order; the name gets a
    suffix (e.g. poa_abpoa_Sp_guide) unless given."""
    tool = 'abpoa' if method == 'poa_abpoa' else 'spoa' if method == 'poa_spoa' else None
    if tool is None:
        raise SystemExit('realign_poa: unknown method %s (poa_abpoa, poa_spoa)' % method)
    d = DEFAULTS[method]
    cfg = config_name or d['config']
    od = order or d['order']
    if tool == 'spoa' and config_name:
        raise SystemExit('realign_poa: --config applies to poa_abpoa only')
    if tool == 'abpoa' and cfg not in ABPOA_CONFIGS and tool_args is None:
        raise SystemExit('realign_poa: unknown abPOA config %s (%s)' % (cfg, ', '.join(ABPOA_CONFIGS)))
    spec = _spec(tool, cfg, od, tool_args)
    if name is None:
        if cfg == d['config'] and od == d['order'] and tool_args is None:
            name = method
        else:
            parts = [method]
            if tool == 'abpoa' and (cfg != d['config'] or tool_args is not None):
                parts.append('custom' if tool_args is not None else cfg)
            if od != d['order'] or len(parts) > 1:
                parts.append(od.replace(':', ''))
            name = '_'.join(parts)
    return name, spec


# ---------------------------------------------------------------- projection (full panel -> hap32)

def project_to_hap32(msa_fa, hap32_fa, out_fa, allow_missing=False):
    """Keep, for every hap32.fa record, an MSA row that spells the same sequence (matched by
    sequence, since full-panel records are not named like hap32 paths; case-insensitive), named
    as in hap32.fa and in its order; drop all-gap columns; write out_fa. Returns a dict with the
    counts; raises ValueError when a hap32 sequence has no identical row (unless allow_missing)."""
    msa = msa_graph.read_msa(msa_fa)
    by_seq = {}
    for n, r in msa:
        by_seq.setdefault(r.replace('-', ''), r)
    hap = msa_graph.read_fasta(hap32_fa)
    rows, missing = [], []
    for n, sq in hap:
        r = by_seq.get(sq.upper())
        if r is None:
            missing.append(n)
        else:
            rows.append((n, r))
    if missing and not allow_missing:
        raise ValueError('%d hap32 sequences have no identical MSA row: %s' % (len(missing), ', '.join(missing[:5])))
    ncol = len(rows[0][1]) if rows else 0
    keep = [c for c in range(ncol) if any(r[c] != '-' for _, r in rows)]
    if len(keep) < ncol:
        rows = [(n, ''.join(r[c] for c in keep)) for n, r in rows]
    msa_graph.write_msa(rows, out_fa)
    return {'rows': len(rows), 'missing': missing, 'columns_in': ncol, 'columns_out': len(keep),
            'msa_rows': len(msa)}


# ---------------------------------------------------------------- fragments (abpoa -i)

# Fragment classes (package_regions.py names, plus the anchor-only pieces realign.py derives):
# a prefix of the region is aligned in extension mode from the left, a suffix in extension mode
# on the reversed MSA and sequence (so from the right), an internal piece in local mode.
FRAG_PASSES = (('prefix', ('enters_L', 'anchor_L'), 2, False),
               ('suffix', ('exits_R', 'anchor_R'), 2, True),
               ('internal', ('internal',), 1, False))


def _drop_gap_columns(rows):
    if not rows:
        return rows
    keep = [c for c in range(len(rows[0])) if any(r[c] != '-' for r in rows)]
    if len(keep) == len(rows[0]):
        return rows
    return [''.join(r[c] for c in keep) for r in rows]


def add_fragments_inc(rows, frags, workdir, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB):
    """Add fragment sequences to an existing MSA with `abpoa -i`, one pass per class (FRAG_PASSES),
    each pass against the MSA of everything added so far (spanning rows and earlier fragments).
    rows: [(name, aligned row)]; frags: [(name, class, sequence)]. Returns (rows, info): the input
    rows followed by one row per fragment, all in the same columns (all-gap columns dropped), or
    (None, info) when a pass fails. The input rows keep their alignment exactly (only inserted
    all-gap columns); every row is checked to spell its sequence."""
    import realign                      # late import: realign imports this module as a plugin
    names = [n for n, _ in rows]
    cur = [r.upper().replace('.', '-') for _, r in rows]
    orig = list(cur)
    info = collections.OrderedDict([('method', 'abpoa -i per class: prefix -m 2, suffix -m 2 on the '
                                     'reversed MSA, internal -m 1'), ('passes', [])])
    known = set(c for _, cls, _, _ in FRAG_PASSES for c in cls)
    other = [n for n, c, _ in frags if c not in known]
    if other:
        info.update(status='error', message='unknown fragment class: %s' % ', '.join(other[:5]))
        return None, info
    for tag, classes, mode, rev in FRAG_PASSES:
        todo = [(n, s.upper()) for n, c, s in frags if c in classes and s]
        if not todo:
            continue
        f = (lambda x: x[::-1]) if rev else (lambda x: x)
        drow = collections.OrderedDict()
        for r in cur:
            drow.setdefault(r, len(drow))
        dseq = collections.OrderedDict()
        for _, s in todo:
            dseq.setdefault(s, len(dseq))
        wd = os.path.join(workdir, 'frag_' + tag)
        os.makedirs(wd, exist_ok=True)
        mf, sf, of = (os.path.join(wd, x) for x in ('msa.fa', 'in.fa', 'out.fa'))
        with open(mf, 'w') as fo:
            for r, k in drow.items():
                fo.write('>m%d\n%s\n' % (k, f(r)))
        with open(sf, 'w') as fo:
            for s, k in dseq.items():
                fo.write('>f%d\n%s\n' % (k, f(s)))
        cmd = [ABPOA, '-m', str(mode), '-r', '1', '-i', mf, sf]
        p = {'pass': tag, 'mode': mode, 'reversed': rev, 'n_fragments': len(todo), 'n_distinct': len(dseq),
             'msa_rows_in': len(drow), 'columns_in': len(cur[0])}
        r = realign.run_proc(cmd, of, os.path.join(wd, 'abpoa.log'), timeout=timeout, mem_mb=mem_mb,
                             env=config.tool_env(), cwd=wd)
        p.update(status=r['status'], seconds=r.get('seconds'), peak_rss_mb=r.get('peak_rss_mb'))
        info['passes'].append(p)
        if r['status'] != 'ok' or not os.path.getsize(of):
            info.update(status='error' if r['status'] == 'ok' else r['status'],
                        message='abpoa -i %s pass: %s' % (tag, r.get('status')))
            return None, info
        got = {n: f(row.upper().replace('.', '-')) for n, row in msa_graph.read_msa(of)}
        # abpoa re-lays gaps of the rows it reads (always on the reversed pass), so the fragments
        # are carried back into the columns of `cur`: a fragment base goes to the old column of
        # the rows that share its new column when they all came from one old column (and it keeps
        # the fragment's column order); otherwise it gets an inserted column of its own, shared
        # with the other fragments' bases of the same new column in the same gap.
        oldset = None
        for r0, k in drow.items():
            row = got.get('m%d' % k)
            if row is None or row.replace('-', '') != r0.replace('-', ''):
                info.update(status='error', message='abpoa -i %s pass changed MSA row m%d' % (tag, k))
                return None, info
            if oldset is None:
                oldset = [set() for _ in range(len(row))]
            for a, b in zip((c for c, x in enumerate(r0) if x != '-'), (c for c, x in enumerate(row) if x != '-')):
                oldset[b].add(a)
        placed, ins_keys = [], set()
        unanchored = 0
        for n, s in todo:
            row = got.get('f%d' % dseq[s])
            if row is None or row.replace('-', '') != s or len(row) != len(oldset):
                info.update(status='error', message='abpoa -i %s pass: fragment %s does not spell its sequence'
                            % (tag, n))
                return None, info
            prev, pl = -1, []
            for c, x in enumerate(row):
                if x == '-':
                    continue
                o = oldset[c]
                if len(o) == 1 and min(o) > prev:
                    prev = min(o)
                    pl.append((prev, None, x))
                else:
                    pl.append((prev, c, x))
                    ins_keys.add((prev, c))
                    unanchored += 1
            placed.append(pl)
        ncol0 = len(cur[0])
        by_prev = collections.defaultdict(list)
        for j, c in ins_keys:
            by_prev[j].append(c)
        order, colpos = [], {}
        for j in range(-1, ncol0):
            if j >= 0:
                colpos[(j, None)] = len(order)
                order.append((j, None))
            for c in sorted(by_prev.get(j, ())):
                colpos[(j, c)] = len(order)
                order.append((j, c))
        ncol = len(order)
        oldpos = [colpos[(j, None)] for j in range(ncol0)]
        new_cur = []
        for r0 in cur:
            buf = ['-'] * ncol
            for j, x in enumerate(r0):
                if x != '-':
                    buf[oldpos[j]] = x
            new_cur.append(''.join(buf))
        fr = []
        for pl in placed:
            buf = ['-'] * ncol
            for j, c, x in pl:
                buf[colpos[(j, None) if c is None else (j, c)]] = x
            fr.append(''.join(buf))
        for (n, s), row in zip(todo, fr):
            assert row.replace('-', '') == s, n
        cur = new_cur + fr
        names += [n for n, _ in todo]
        p['columns_out'] = ncol
        p['relaid_columns'] = sum(1 for o in oldset if len(o) > 1)
        p['fragment_bases_unanchored'] = unanchored
    cur = _drop_gap_columns(cur)
    proj = _drop_gap_columns(cur[:len(orig)])
    info['input_alignment_preserved'] = proj == orig
    info['columns'] = len(cur[0]) if cur else 0
    info['n_fragment_rows'] = len(cur) - len(orig)
    info['status'] = 'ok'
    return list(zip(names, cur)), info


# ---------------------------------------------------------------- CLI

def main(argv=None):
    import realign
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('regions', nargs='*', help="region dirs or ids; 'all' = every packaged region")
    ap.add_argument('--method', required=True, choices=list(METHODS))
    ap.add_argument('--out', help='candidates root (default %s)' % config.CANDIDATES_DIR)
    ap.add_argument('--config', dest='config_name', help='abPOA flag set: %s' % ', '.join(ABPOA_CONFIGS))
    ap.add_argument('--order', help='given, longest, shortest, guide, random:SEED')
    ap.add_argument('--tool-args', help="extra flags for the tool, one string (replaces --config's flags)")
    ap.add_argument('--name', help='method name for a variant (default: the method plus a suffix)')
    ap.add_argument('--input', help='MSA-only mode: align this FASTA ...')
    ap.add_argument('--msa-out', help='... and write only the MSA here')
    ap.add_argument('--json', help='MSA-only mode: also write the info dict here')
    ap.add_argument('--project-to', help='MSA-only mode: also project the MSA onto this hap32.fa (rows matched '
                                         'by sequence) ...')
    ap.add_argument('--projected-out', help='... and write the projected MSA here (hap32 names and order)')
    ap.add_argument('--timeout', type=float, default=DEFAULT_TIMEOUT, help='seconds per alignment (default 900)')
    ap.add_argument('--mem-mb', type=float, default=DEFAULT_MEM_MB,
                    help='memory cap: predicted or sampled RSS above it stops the region (default %d)'
                         % DEFAULT_MEM_MB)
    ap.add_argument('--jobs', type=int, default=1, help='regions at once (default 1; at most 2 advised)')
    ap.add_argument('--stratum', help='comma list of strata (with all)')
    ap.add_argument('--force', action='store_true', help='re-run regions that already have a result')
    ap.add_argument('--retry-failed', action='store_true', help='re-run only earlier failures and timeouts')
    ap.add_argument('--no-runtime', action='store_true', help='do not write results/realign_runtime.tsv')
    ap.add_argument('--no-dedup', action='store_true', help='align identical sequences separately')
    ap.add_argument('--merge-blocks', action='store_true', help='pass --merge-blocks to msa_graph (not the default)')
    ap.add_argument('--workdir', help='scratch for aligner runs (default $VNTR_WORK/realign)')
    ap.add_argument('--no-wait-memory', dest='wait_memory', action='store_false',
                    help='region mode: do not wait for free memory before a region predicted to need > 2 GB')
    a = ap.parse_args(argv)

    targs = a.tool_args.split() if a.tool_args is not None else None
    name, spec = variant(a.method, a.config_name, a.order, targs, a.name)
    ms = realign.all_methods()
    ms[name] = realign._as_method(name, spec)       # register the variant for realign's pipeline
    wd = a.workdir or os.path.join(config.WORK_DIR, 'realign')
    os.makedirs(wd, exist_ok=True)

    if a.input or a.msa_out:
        if not (a.input and a.msa_out) or a.regions:
            ap.error('MSA-only mode takes --input FASTA and --msa-out FILE and no regions')
        if bool(a.project_to) != bool(a.projected_out):
            ap.error('--project-to and --projected-out go together')
        info = realign.align_fasta(name, a.input, a.msa_out, threads=1, timeout=a.timeout, mem_mb=a.mem_mb,
                                   workdir=wd, dedup=not a.no_dedup)
        if a.json:
            with open(a.json, 'w') as f:
                json.dump(info, f, indent=1)
        sys.stderr.write('realign_poa: %s %s: %s, %s rows (%s distinct), %s columns, %.1f s%s\n' % (
            name, a.input, info['status'], info.get('n_seqs'), info.get('n_distinct'), info.get('columns'),
            info.get('align_s', 0), (' -- ' + info['message']) if info.get('message') else ''))
        if info['status'] == 'ok' and a.project_to:
            try:
                pr = project_to_hap32(a.msa_out, a.project_to, a.projected_out)
            except ValueError as e:
                sys.stderr.write('realign_poa: projection failed: %s\n' % e)
                return 2
            sys.stderr.write('realign_poa: projected %d hap32 rows, %d -> %d columns: %s\n' % (
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
    kw = dict(threads=1, timeout=a.timeout, mem_mb=a.mem_mb, force=a.force, retry_failed=a.retry_failed,
              dedup=not a.no_dedup, cand_root=a.out, workroot=wd, record_runtime=not a.no_runtime,
              merge_blocks=a.merge_blocks)
    t0 = time.time()
    counts = collections.Counter()
    worst = 0

    tool = spec['tool']
    flags = [ABPOA] + ABPOA_CONFIGS.get(spec['params'].get('config_name') or '', []) if tool == 'abpoa' else [SPOA]
    if spec['params'].get('tool_args'):
        flags = flags + list(spec['params']['tool_args'])

    def job(rd):
        realign.all_methods()[name] = realign._as_method(name, spec)
        # a big region waits until the machine has the memory it is predicted to need
        if a.wait_memory and (a.force or a.retry_failed or not realign.previous_status(name, os.path.basename(rd), a.out)):
            try:
                lens = [len(s) for _, s in msa_graph.read_fasta(os.path.join(rd, 'hap32.fa'))]
                need = predict_mb(tool, lens, flags)
                if 2000 < need <= a.mem_mb:
                    w = wait_for_memory(need + 1500, log=lambda m: sys.stderr.write(
                        'realign_poa: %s %s\n' % (os.path.basename(rd), m)))
                    if w:
                        sys.stderr.write('realign_poa: %s waited %.0f s for memory\n' % (os.path.basename(rd), w))
            except OSError:
                pass
        try:
            return realign.realign_region(name, rd, **kw)
        except Exception as e:  # keep the batch going
            return {'method': name, 'region_id': os.path.basename(rd), 'status': 'crash',
                    'note': '%s: %s' % (type(e).__name__, e)}

    if a.jobs <= 1:
        it = map(job, rdirs)
        ex = None
    else:
        ex = concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs)   # the work is in subprocesses
        it = ex.map(job, rdirs)
    for row in it:
        counts[row['status']] += 1
        if row['status'] not in ('ok', 'skipped'):
            worst = 2
        sys.stderr.write('realign_poa: %-18s %-9s %-8s align %6ss graph %6ss nodes %-6s %s\n' % (
            name, row['region_id'], row['status'], row.get('align_s', ''), row.get('graph_s', ''),
            row.get('nodes', ''), row.get('note', '')))
        sys.stderr.flush()
    if ex:
        ex.shutdown()
    sys.stderr.write('realign_poa: %s done in %.0f s: %s\n' % (
        name, time.time() - t0, ', '.join('%s %d' % kv for kv in sorted(counts.items()))))
    return worst


if __name__ == '__main__':
    sys.exit(main())
