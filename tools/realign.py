#!/usr/bin/env python3
"""realign.py -- re-align a region's haplotype sequences and build a candidate graph from the MSA.

Region mode (the default): align regions/<id>/hap32.fa with METHOD, then build the graph with
msa_graph.py, and write

    candidates/<METHOD>/<id>.msa.fa        aligned rows, '-' gaps, one line per row, hap32.fa names
                                           and order; every row spells its hap32.fa sequence
    candidates/<METHOD>/<id>.gfa           GFA 1.0, one P line per hap32.fa sequence (msa_graph.py)
    candidates/<METHOD>/<id>.realign.json  how it was made: command, versions, N masking, timings,
                                           graph stats; written for failures and timeouts too
    results/realign_runtime.tsv            one row per (method, region), the latest run wins

    python3 tools/realign.py mafft_linsi L012184                    # one region (id or regions/<id>)
    python3 tools/realign.py mafft_fftns2,mafft_linsi all --jobs 2 --threads 2
    python3 tools/realign.py mafft_einsi all --stratum hotspot_vntr --timeout 900
    python3 tools/realign.py --list                                  # the method table

A region whose candidate already exists (or whose last attempt timed out) is skipped unless
--force; --retry-failed re-runs only failures and timeouts, --retry memout[,timeout,...] only
those statuses (e.g. memouts again with --threads 1 and a larger --mem-mb).

MSA-only mode, for an ARBITRARY input FASTA (e.g. the deduplicated full HPRC panel): write only
the MSA, no graph, no candidates/ files, no runtime row:

    python3 tools/realign.py mafft_linsi --input seqs.fa[.gz] --msa-out out.msa.fa \\
                             [--threads 2] [--timeout 900] [--mem-mb 10000] [--json info.json]

    import sys; sys.path.insert(0, 'tools')
    from realign import align_fasta
    info = align_fasta('mafft_linsi', 'seqs.fa', 'out.msa.fa', threads=2, timeout=900)
    info['status']      # 'ok', 'timeout', 'memout' or 'error' (out.msa.fa exists only when 'ok')

The CLI exits 0 when the MSA was written, 4 on a timeout or memout, 2 on an error.

Input rows are named by the first word of their header (names must be unique); the MSA has one
row per input record, in input order, named the same, upper case, '-' gaps, one line per row, and
every row's ungapped sequence equals its input sequence (checked). Projecting such an MSA onto a
subset of rows (keep the rows, drop all-gap columns) is msa_graph.match_rows(..., ignore_extra=True).

What every method gets (the shared pipeline, align_fasta):
  1. Identical input sequences are aligned once (--no-dedup to turn off) and their rows copied.
  2. N masking: every maximal run of non-ACGT characters (in practice only N; --mask-min-run sets
     the shortest run masked, default 1) is cut out before alignment, so an assembly gap is never
     aligned against real sequence. After alignment the run is put back as columns of its own,
     right after the column of the base that precedes it, in which only that row has bases (the
     N-run haplotype reads through the columns of its remaining sequence, as a deletion, and the
     run is an insertion). So the MSA still spells the original sequence, and msa_graph turns each
     run into one N node of its own: an N-gap haplotype adds one N node and no other node. evaluate.py
     requires paths to spell the Ns (N matches only N), so path gaps are not an option. The count
     of runs and N bp is in the .realign.json and the runtime table.
  3. The aligner sees the distinct masked sequences, upper case, named s0..sK (no name mangling),
     under a wall-clock timeout (--timeout, default 900 s) and a memory cap on its whole process
     group (--mem-mb, default 10000; sampled RSS), and its output is checked row by row.
     Measured on the 149 regions: at --threads 2, L-/E-/G-INS-i exceed 10 GB once the longest
     sequence passes ~20 kb (E-INS-i from ~12 kb); with --threads 1 the progressive step needs
     about half, and every such region up to 33.6 kb fitted in 16 GB (peak 14.5 GB). 55 kb
     (L012272, 42 kb of it flank) does not fit.

Methods (the table is METHODS below; --list prints it):
  mafft_fftns2  mafft --retree 2 --maxiterate 0                  FFT-NS-2, progressive (the default)
  mafft_linsi   mafft --localpair --maxiterate 1000               L-INS-i
  mafft_einsi   mafft --genafpair --maxiterate 1000 --ep 0        E-INS-i
  mafft_ginsi   mafft --globalpair --maxiterate 1000              G-INS-i
  mafft_fftnsi  mafft --retree 2 --maxiterate 1000                FFT-NS-i (for the full-panel arm)
Every mafft call adds --nuc --thread T --threadit 0 (iterative refinement single-threaded, which
mafft documents as needed for reproducible output) and runs with TMPDIR in a private work dir.

Extending the table: any module tools/realign_*.py that defines METHODS = {name: spec} is picked
up (tools/realign_poa.py adds poa_*). A spec is a dict

    {'align': func, 'params': {...}, 'description': '...', 'tool': 'abpoa'}

(or a Method, or a bare function), where func(in_fa, out_fa, threads=, workdir=, timeout=,
mem_mb=, **params) aligns the sequences of in_fa (distinct, upper-case ACGT, named s0..sK) into
out_fa (aligned FASTA or PIR, rows named as in in_fa, any gap character and case) and returns a
dict: {'status': 'ok'|'timeout'|'memout'|'error', 'command': [...], 'message': '...', ...};
raising an exception counts as 'error'. run_proc() (timeout + memory cap on a process group) is
there to be reused. The shared pipeline above (dedup, N masking, checks, graph, outputs) applies
to every method, built in or not.

Fragments (optional, --fragments): after a successful alignment, hap32.fragments.fa (pieces of
haplotypes that do not span the region) is added with mafft --addfragments --keeplength and
written to candidates/<METHOD>/<id>.fragments.msa.fa (fragment rows only, in the columns of
<id>.msa.fa; --keeplength deletes fragment bases that would need a new column, so these rows do
not always spell their fragment). The graph never includes fragments.
"""
import argparse
import collections
import concurrent.futures
import csv
import datetime
import fcntl
import glob
import gzip
import importlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402
import msa_graph  # noqa: E402

RUNTIME_TSV = os.path.join(config.RESULTS_DIR, 'realign_runtime.tsv')
RUNTIME_COLS = ['method', 'region_id', 'stratum', 'span_bp', 'n_seqs', 'n_distinct', 'max_len', 'total_bp',
                'status', 'align_s', 'graph_s', 'total_s', 'peak_rss_mb', 'threads', 'timeout_s',
                'columns', 'gap_frac', 'nodes', 'nodes_per_kb', 'frac_nodes_1bp', 'masked_runs',
                'masked_bp', 'N_nodes', 'finished', 'note']
DEFAULT_TIMEOUT = 900
DEFAULT_MEM_MB = 10000
ACGT = re.compile(r'[^ACGT]+')


# ---------------------------------------------------------------- process control

def _group_rss_mb(pgid):
    """Summed RSS (MB) of the processes in process group pgid (sampled with ps)."""
    try:
        out = subprocess.run(['ps', '-A', '-o', 'pgid=,rss='], capture_output=True, text=True,
                             timeout=20).stdout
    except (OSError, subprocess.SubprocessError):
        return 0.0
    kb = 0
    for line in out.splitlines():
        f = line.split()
        if len(f) == 2 and f[0] == str(pgid):
            try:
                kb += int(f[1])
            except ValueError:
                pass
    return kb / 1024.0


def _kill_group(p):
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(p.pid, sig)
        except (ProcessLookupError, PermissionError):
            return
        try:
            p.wait(timeout=5)
            return
        except subprocess.TimeoutExpired:
            continue


def run_proc(cmd, stdout_path=None, stderr_path=None, timeout=None, mem_mb=None, env=None, cwd=None,
             poll=1.0):
    """Run cmd in its own process group. stdout/stderr go to files (or /dev/null). The group is
    killed when it runs longer than timeout seconds or its summed RSS exceeds mem_mb.
    Returns {'status': ok|timeout|memout|error, 'returncode', 'seconds', 'peak_rss_mb'}."""
    out = open(stdout_path, 'wb') if stdout_path else subprocess.DEVNULL
    err = open(stderr_path, 'wb') if stderr_path else subprocess.DEVNULL
    t0 = time.time()
    peak = 0.0
    status = None
    try:
        p = subprocess.Popen(cmd, stdout=out, stderr=err, env=env, cwd=cwd, start_new_session=True)
        wait = min(0.2, poll)
        while True:
            try:
                p.wait(timeout=wait)
                break
            except subprocess.TimeoutExpired:
                pass
            wait = poll
            peak = max(peak, _group_rss_mb(p.pid))
            if timeout and time.time() - t0 > timeout:
                status = 'timeout'
            elif mem_mb and peak > mem_mb:
                status = 'memout'
            if status:
                _kill_group(p)
                break
    finally:
        if stdout_path:
            out.close()
        if stderr_path:
            err.close()
    rc = p.returncode
    if status is None:
        status = 'ok' if rc == 0 else 'error'
    return {'status': status, 'returncode': rc, 'seconds': round(time.time() - t0, 2),
            'peak_rss_mb': round(peak, 1)}


def _tail(path, n=6):
    try:
        with open(path, errors='replace') as f:
            lines = [ln.rstrip() for ln in f if ln.strip()]
        return ' | '.join(lines[-n:])[-600:]
    except OSError:
        return ''


# ---------------------------------------------------------------- mafft

_MAFFT_VERSION = None


def mafft_version():
    global _MAFFT_VERSION
    if _MAFFT_VERSION is None:
        try:
            p = subprocess.run([config.MAFFT, '--version'], capture_output=True, text=True,
                               env=config.tool_env(), timeout=60)
            _MAFFT_VERSION = (p.stdout + p.stderr).strip().splitlines()[0]
        except (OSError, subprocess.SubprocessError, IndexError):
            _MAFFT_VERSION = 'unknown'
    return _MAFFT_VERSION


def mafft_align(in_fa, out_fa, threads=2, workdir=None, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB,
                args=(), threadit=0, **_):
    """mafft ARGS --nuc --thread T --threadit 0 in_fa > out_fa, with TMPDIR private to the call."""
    wd = tempfile.mkdtemp(prefix='mafft.', dir=workdir)
    cmd = [config.MAFFT] + list(args) + ['--nuc', '--thread', str(threads), '--threadit', str(threadit), in_fa]
    env = dict(config.tool_env(), TMPDIR=wd)
    log = os.path.join(os.path.dirname(os.path.abspath(out_fa)), os.path.basename(out_fa) + '.mafft.log')
    tmp_out = out_fa + '.part'
    r = run_proc(cmd, tmp_out, log, timeout=timeout, mem_mb=mem_mb, env=env, cwd=wd)
    shutil.rmtree(wd, ignore_errors=True)
    r['command'] = ['mafft'] + cmd[1:-1] + ['<in.fa>']
    r['tool_version'] = mafft_version()
    if r['status'] == 'ok' and os.path.getsize(tmp_out) == 0:
        r['status'] = 'error'
    if r['status'] == 'ok':
        os.replace(tmp_out, out_fa)
    else:
        r['message'] = _tail(log)
        if os.path.exists(tmp_out):
            os.remove(tmp_out)
    if os.path.exists(log) and r['status'] == 'ok':
        os.remove(log)
    return r


def mafft_add_fragments(msa_fa, frag_fa, out_fa, threads=2, workdir=None, timeout=DEFAULT_TIMEOUT,
                        mem_mb=DEFAULT_MEM_MB):
    """mafft --addfragments frag_fa --keeplength msa_fa > out_fa."""
    wd = tempfile.mkdtemp(prefix='mafftadd.', dir=workdir)
    cmd = [config.MAFFT, '--addfragments', frag_fa, '--keeplength', '--nuc', '--thread', str(threads),
           '--threadit', '0', msa_fa]
    log = out_fa + '.mafft.log'
    r = run_proc(cmd, out_fa, log, timeout=timeout, mem_mb=mem_mb, env=dict(config.tool_env(), TMPDIR=wd),
                 cwd=wd)
    shutil.rmtree(wd, ignore_errors=True)
    if r['status'] != 'ok':
        r['message'] = _tail(log)
    elif os.path.exists(log):
        os.remove(log)
    r['command'] = ['mafft'] + cmd[1:2] + ['<fragments.fa>'] + cmd[3:-1] + ['<msa.fa>']
    return r


class Method(object):
    """A realigner: align(in_fa, out_fa, threads=, workdir=, timeout=, mem_mb=, **params) -> dict."""

    def __init__(self, name, align, params=None, description='', tool=''):
        self.name, self.align, self.params = name, align, dict(params or {})
        self.description, self.tool = description, tool

    def run(self, in_fa, out_fa, threads, workdir, timeout, mem_mb):
        try:
            r = self.align(in_fa, out_fa, threads=threads, workdir=workdir, timeout=timeout, mem_mb=mem_mb,
                           **self.params)
        except Exception as e:  # a plugin's failure is recorded, not fatal
            return {'status': 'error', 'message': '%s: %s' % (type(e).__name__, e)}
        r = dict(r or {})
        r.setdefault('status', 'ok' if os.path.exists(out_fa) else 'error')
        return r


def _mafft(name, args, desc):
    return Method(name, mafft_align, {'args': list(args)}, desc, 'mafft')


METHODS = collections.OrderedDict((m.name, m) for m in [
    _mafft('mafft_fftns2', ['--retree', '2', '--maxiterate', '0'], 'mafft FFT-NS-2: progressive, two guide trees'),
    _mafft('mafft_linsi', ['--localpair', '--maxiterate', '1000'],
           'mafft L-INS-i: local pairwise consistency + iterative refinement'),
    _mafft('mafft_einsi', ['--genafpair', '--maxiterate', '1000', '--ep', '0'],
           'mafft E-INS-i: generalized affine pairwise (unalignable stretches allowed) + refinement'),
    _mafft('mafft_ginsi', ['--globalpair', '--maxiterate', '1000'],
           'mafft G-INS-i: global pairwise consistency + iterative refinement'),
    _mafft('mafft_fftnsi', ['--retree', '2', '--maxiterate', '1000'],
           'mafft FFT-NS-i: FFT-NS-2 + iterative refinement (for the full-panel arm)'),
])


def _as_method(name, spec):
    if hasattr(spec, 'align') and callable(getattr(spec, 'align')):
        return Method(name, spec.align, getattr(spec, 'params', {}), getattr(spec, 'description', ''),
                      getattr(spec, 'tool', ''))
    if isinstance(spec, dict) and callable(spec.get('align')):
        return Method(name, spec['align'], spec.get('params'), spec.get('description', ''), spec.get('tool', ''))
    if callable(spec):
        return Method(name, spec, {}, getattr(spec, '__doc__', '') or '', '')
    raise ValueError('method %s: not a Method, dict with align, or function' % name)


_ALL = None
PLUGIN_ERRORS = {}


def all_methods():
    """The built-in table plus METHODS from every tools/realign_*.py module."""
    global _ALL
    if _ALL is None:
        _ALL = collections.OrderedDict(METHODS)
        for path in sorted(glob.glob(os.path.join(TOOLS, 'realign_*.py'))):
            mod = os.path.basename(path)[:-3]
            try:
                m = importlib.import_module(mod)
                for name, spec in (getattr(m, 'METHODS', None) or {}).items():
                    _ALL[name] = _as_method(name, spec)
            except Exception as e:
                PLUGIN_ERRORS[mod] = '%s: %s' % (type(e).__name__, e)
    return _ALL


def get_method(name):
    ms = all_methods()
    if name not in ms:
        extra = (' (plugin import errors: %s)' % PLUGIN_ERRORS) if PLUGIN_ERRORS else ''
        raise KeyError('unknown method %s; known: %s%s' % (name, ', '.join(ms), extra))
    return ms[name]


# ---------------------------------------------------------------- the shared pipeline

def read_fasta(path):
    """[(first word of header, sequence)]; sequences may wrap; .gz is read."""
    return msa_graph.read_fasta(path)


def mask_runs(seq, min_run=1):
    """Cut maximal non-ACGT runs of >= min_run characters out of seq.
    Returns (masked, runs) with runs = [(bases of masked before the run, run string)]."""
    parts, runs, last, kept = [], [], 0, 0
    for m in ACGT.finditer(seq):
        if m.end() - m.start() < min_run:
            continue
        parts.append(seq[last:m.start()])
        kept += m.start() - last
        runs.append((kept, m.group()))
        last = m.end()
    parts.append(seq[last:])
    return ''.join(parts), runs


def reinsert_runs(rows, runs_by_row):
    """rows: aligned strings; runs_by_row: {row index: [(pos, run)]} from mask_runs. Each run is
    inserted as new columns, after the column of the row's base pos-1 (at the start when pos 0),
    in which only that row has bases. Returns new rows."""
    if not any(runs_by_row.values()):
        return list(rows)
    ncol = len(rows[0]) if rows else 0
    ins = collections.defaultdict(list)       # column index -> [(row, run order, run)]
    for r, runs in sorted(runs_by_row.items()):
        if not runs:
            continue
        cols = [c for c, ch in enumerate(rows[r]) if ch != '-']
        for k, (pos, run) in enumerate(runs):
            at = cols[pos - 1] + 1 if pos > 0 else (cols[0] if cols else 0)
            ins[at].append((r, k, run))
    out = [[] for _ in rows]
    for c in range(ncol + 1):
        for r, _k, run in sorted(ins.get(c, [])):
            for i in range(len(rows)):
                out[i].append(run if i == r else '-' * len(run))
        if c < ncol:
            for i, row in enumerate(rows):
                out[i].append(row[c])
    return [''.join(x) for x in out]


def _write_fa(path, recs):
    with open(path, 'w') as f:
        for n, s in recs:
            f.write('>%s\n%s\n' % (n, s))


def align_fasta(method, input_fa, msa_out, threads=2, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB,
                workdir=None, dedup=True, mask_min_run=1, keep_work=False):
    """Align every record of input_fa with METHOD and write the MSA to msa_out (only when the status
    is 'ok'). Returns an info dict (status, counts, masking, timings, the aligner's command)."""
    t0 = time.time()
    m = get_method(method) if isinstance(method, str) else method
    recs = read_fasta(input_fa)
    info = collections.OrderedDict([('method', m.name), ('input', os.path.abspath(input_fa)),
                                    ('msa', os.path.abspath(msa_out)), ('status', None)])
    names = [n for n, _ in recs]
    if not recs:
        info.update(status='error', message='no sequences in %s' % input_fa)
        return info
    dup = [n for n, c in collections.Counter(names).items() if c > 1]
    if dup:
        info.update(status='error', message='duplicate names in %s: %s' % (input_fa, ', '.join(dup[:5])))
        return info
    seqs = [s.upper() for _, s in recs]
    # distinct sequences (first occurrence order)
    if dedup:
        idx = collections.OrderedDict()
        for s in seqs:
            idx.setdefault(s, len(idx))
        distinct = list(idx)
        row_of = [idx[s] for s in seqs]
    else:
        distinct, row_of = list(seqs), list(range(len(seqs)))
    masked, runs = [], {}
    for k, s in enumerate(distinct):
        ms, rs = mask_runs(s, mask_min_run)
        masked.append(ms)
        runs[k] = rs
    info['n_seqs'] = len(recs)
    info['n_distinct'] = len(distinct)
    info['max_len'] = max(len(s) for s in seqs)
    info['total_bp'] = sum(len(s) for s in seqs)
    info['masked_runs'] = sum(len(runs[row_of[i]]) for i in range(len(recs)))
    info['masked_bp'] = sum(sum(len(x) for _, x in runs[row_of[i]]) for i in range(len(recs)))
    info['masked_runs_distinct'] = sum(len(v) for v in runs.values())
    info['masked_detail'] = [{'row': names[row_of.index(k)], 'offset_masked': p, 'len': len(x),
                              'chars': ''.join(sorted(set(x)))}
                             for k in sorted(runs) for p, x in runs[k]][:50]
    info['threads'] = threads
    info['timeout_s'] = timeout
    info['mem_mb'] = mem_mb
    wd = tempfile.mkdtemp(prefix='realign.%s.' % m.name, dir=workdir)
    info['workdir'] = wd if keep_work else None
    try:
        # the aligner's input: distinct, non-empty masked sequences
        todo = [k for k, s in enumerate(masked) if s]
        aligned = {}
        if len(todo) == 1:
            aligned[todo[0]] = masked[todo[0]]
            info['aligner'] = {'status': 'ok', 'note': 'one distinct sequence; no aligner run', 'seconds': 0.0}
        elif todo:
            in_fa = os.path.join(wd, 'in.fa')
            out_fa = os.path.join(wd, 'out.fa')
            _write_fa(in_fa, [('s%d' % k, masked[k]) for k in todo])
            r = m.run(in_fa, out_fa, threads, wd, timeout, mem_mb)
            info['aligner'] = r
            if r.get('status') != 'ok':
                info['status'] = r.get('status') or 'error'
                info['message'] = r.get('message', '')
                return info
            got = dict(msa_graph.read_msa(out_fa))
            for k in todo:
                row = got.get('s%d' % k)
                if row is None:
                    info.update(status='error', message='aligner output has no row s%d' % k)
                    return info
                if row.replace('-', '') != masked[k]:
                    info.update(status='error', message='aligner row s%d does not spell its input '
                                '(%d vs %d bp)' % (k, len(row.replace('-', '')), len(masked[k])))
                    return info
                aligned[k] = row
            lens = set(len(v) for v in aligned.values())
            if len(lens) != 1:
                info.update(status='error', message='aligner rows have different lengths %s' % sorted(lens)[:5])
                return info
        ncol = len(next(iter(aligned.values()))) if aligned else 0
        rows = [aligned.get(k, '-' * ncol) for k in range(len(distinct))]
        rows = reinsert_runs(rows, runs)
        # drop all-gap columns
        if rows:
            keep = [c for c in range(len(rows[0])) if any(r[c] != '-' for r in rows)]
            if len(keep) < len(rows[0]):
                rows = [''.join(r[c] for c in keep) for r in rows]
        for k, s in enumerate(distinct):
            if rows[k].replace('-', '') != s:
                info.update(status='error', message='internal: row for %s does not spell its sequence'
                            % names[row_of.index(k)])
                return info
        out_rows = [(names[i], rows[row_of[i]]) for i in range(len(recs))]
        tmp = msa_out + '.tmp%d' % os.getpid()
        os.makedirs(os.path.dirname(os.path.abspath(msa_out)) or '.', exist_ok=True)
        msa_graph.write_msa(out_rows, tmp)
        os.replace(tmp, msa_out)
        info['status'] = 'ok'
        info['columns'] = len(rows[0]) if rows else 0
        info['gap_frac'] = round(sum(r.count('-') for _, r in out_rows) / float(max(1, len(out_rows) * info['columns'])), 4)
        return info
    finally:
        info['align_s'] = round(time.time() - t0, 2)
        if not keep_work:
            shutil.rmtree(wd, ignore_errors=True)


# ---------------------------------------------------------------- region mode

def region_list():
    """Region rows from regions/regions.tsv (or every regions/*/region.json), smallest span first."""
    tsv = os.path.join(config.REGIONS_DIR, 'regions.tsv')
    out = []
    if os.path.exists(tsv):
        with open(tsv) as f:
            out = list(csv.DictReader(f, delimiter='\t'))
    else:
        for p in sorted(glob.glob(os.path.join(config.REGIONS_DIR, '*', 'region.json'))):
            out.append(json.load(open(p)))
    out = [r for r in out if os.path.exists(os.path.join(config.REGIONS_DIR, r['region_id'], 'hap32.fa'))]
    for r in out:
        try:
            r['_bp'] = int(r['span_end']) - int(r['span_start']) + 1
        except (KeyError, ValueError, TypeError):
            r['_bp'] = 0
    return sorted(out, key=lambda r: (r['_bp'], r['region_id']))


def region_dir(x):
    if os.path.isdir(x):
        return os.path.abspath(x)
    d = os.path.join(config.REGIONS_DIR, x)
    if os.path.isdir(d):
        return d
    raise SystemExit('realign: no region directory %s' % x)


def region_meta(rd):
    p = os.path.join(rd, 'region.json')
    try:
        j = json.load(open(p))
    except (OSError, ValueError):
        j = {}
    span = None
    try:
        span = int(j['span_end']) - int(j['span_start']) + 1
    except (KeyError, ValueError, TypeError):
        pass
    return j.get('stratum', ''), span


def _lock_update_runtime(row):
    """Replace the (method, region_id) row of results/realign_runtime.tsv, under a file lock."""
    os.makedirs(os.path.dirname(RUNTIME_TSV), exist_ok=True)
    lock = RUNTIME_TSV + '.lock'
    with open(lock, 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        rows = []
        if os.path.exists(RUNTIME_TSV):
            with open(RUNTIME_TSV) as f:
                rows = list(csv.DictReader(f, delimiter='\t'))
        rows = [r for r in rows if (r.get('method'), r.get('region_id')) != (row['method'], row['region_id'])]
        rows.append(row)
        order = dict((m, i) for i, m in enumerate(all_methods()))
        rows.sort(key=lambda r: (order.get(r['method'], 999), r['method'], r['region_id']))
        tmp = RUNTIME_TSV + '.tmp%d' % os.getpid()
        with open(tmp, 'w', newline='') as f:
            w = csv.DictWriter(f, RUNTIME_COLS, delimiter='\t', extrasaction='ignore', lineterminator='\n')
            w.writeheader()
            for r in rows:
                w.writerow(dict((c, r.get(c, '')) for c in RUNTIME_COLS))
        os.replace(tmp, RUNTIME_TSV)
        fcntl.flock(lk, fcntl.LOCK_UN)


def _rel(p):
    """p relative to the repository when it is inside it, else absolute."""
    p = os.path.abspath(p)
    return os.path.relpath(p, config.REPO) if p.startswith(config.REPO + os.sep) else p


def outputs(method, rid, cand_root=None):
    d = os.path.join(cand_root or config.CANDIDATES_DIR, method)
    return {'dir': d, 'msa': os.path.join(d, rid + '.msa.fa'), 'gfa': os.path.join(d, rid + '.gfa'),
            'json': os.path.join(d, rid + '.realign.json'), 'frag': os.path.join(d, rid + '.fragments.msa.fa')}


def previous_status(method, rid, cand_root=None):
    o = outputs(method, rid, cand_root)
    try:
        j = json.load(open(o['json']))
    except (OSError, ValueError):
        return None
    st = j.get('status')
    if st == 'ok' and not (os.path.exists(o['gfa']) and os.path.exists(o['msa'])):
        return None
    return st


def realign_region(method, rdir, threads=2, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB, force=False,
                   retry_failed=False, fragments=False, dedup=True, mask_min_run=1, cand_root=None,
                   workroot=None, record_runtime=True, merge_blocks=False):
    """Align regions/<id>/hap32.fa, build the graph, write candidates/<method>/<id>.*; returns the
    runtime row (status 'skipped' when an earlier result is kept)."""
    rid = os.path.basename(os.path.normpath(rdir))
    m = get_method(method)
    o = outputs(m.name, rid, cand_root)
    prev = previous_status(m.name, rid, cand_root)
    if prev and not force:
        if isinstance(retry_failed, str):
            retry_failed = set(retry_failed.split(','))
        if isinstance(retry_failed, (set, frozenset, list, tuple)):
            again = prev in set(retry_failed)
        else:
            again = bool(retry_failed) and prev != 'ok'
        if not again:
            return {'method': m.name, 'region_id': rid, 'status': 'skipped', 'note': 'kept %s' % prev}
    os.makedirs(o['dir'], exist_ok=True)
    workroot = workroot or os.path.join(config.WORK_DIR, 'realign')
    os.makedirs(workroot, exist_ok=True)
    stratum, span = region_meta(rdir)
    hap = os.path.join(rdir, 'hap32.fa')
    t0 = time.time()
    for k in ('msa', 'gfa', 'frag'):          # never leave an older result beside a new attempt
        if os.path.exists(o[k]):
            os.remove(o[k])
    tmp_msa = o['msa'] + '.tmp%d' % os.getpid()
    info = align_fasta(m, hap, tmp_msa, threads=threads, timeout=timeout, mem_mb=mem_mb, workdir=workroot,
                       dedup=dedup, mask_min_run=mask_min_run)
    info['region_id'] = rid
    info['stratum'] = stratum
    info['span_bp'] = span
    info['description'] = m.description
    info['msa'] = _rel(o['msa']) if info['status'] == 'ok' else None
    graph_s = None
    if info['status'] == 'ok':
        tg = time.time()
        tmp_gfa = o['gfa'] + '.tmp%d' % os.getpid()
        try:
            st = msa_graph.msa_to_gfa(tmp_msa, hap, tmp_gfa, merge_blocks=merge_blocks, workdir=workroot)
            os.replace(tmp_msa, o['msa'])
            os.replace(tmp_gfa, o['gfa'])
            st['msa'] = _rel(o['msa'])
            st['gfa'] = _rel(o['gfa'])
            st['hap32'] = _rel(hap)
            info['graph'] = st
            info['gfa'] = st['gfa']
        except msa_graph.MsaGraphError as e:
            info['status'] = 'graph_error'
            info['message'] = str(e)
        graph_s = round(time.time() - tg, 2)
        for p in (tmp_msa, tmp_gfa):
            if os.path.exists(p):
                os.remove(p)
    if info['status'] == 'ok' and fragments:
        frag = os.path.join(rdir, 'hap32.fragments.fa')
        if os.path.exists(frag) and os.path.getsize(frag) > 0:
            info['fragments'] = add_fragments(o['msa'], frag, o['frag'], threads, workroot, timeout, mem_mb)
        else:
            info['fragments'] = {'status': 'none', 'note': 'no hap32.fragments.fa records'}
    info['graph_s'] = graph_s
    info['total_s'] = round(time.time() - t0, 2)
    info['finished'] = datetime.datetime.now().isoformat(timespec='seconds')
    info['host_load'] = round(os.getloadavg()[0], 1)
    info['realign_version'] = 1
    with open(o['json'] + '.tmp', 'w') as f:
        json.dump(info, f, indent=1)
    os.replace(o['json'] + '.tmp', o['json'])
    g = info.get('graph') or {}
    al = info.get('aligner') or {}
    row = {'method': m.name, 'region_id': rid, 'stratum': stratum, 'span_bp': span,
           'n_seqs': info.get('n_seqs'), 'n_distinct': info.get('n_distinct'), 'max_len': info.get('max_len'),
           'total_bp': info.get('total_bp'), 'status': info['status'], 'align_s': info.get('align_s'),
           'graph_s': graph_s, 'total_s': info['total_s'], 'peak_rss_mb': al.get('peak_rss_mb', ''),
           'threads': threads, 'timeout_s': timeout, 'columns': info.get('columns', ''),
           'gap_frac': info.get('gap_frac', ''), 'nodes': g.get('nodes', ''), 'nodes_per_kb': g.get('nodes_per_kb', ''),
           'frac_nodes_1bp': g.get('frac_nodes_1bp', ''), 'masked_runs': info.get('masked_runs'),
           'masked_bp': info.get('masked_bp'), 'N_nodes': g.get('N_nodes', ''), 'finished': info['finished'],
           'note': (info.get('message') or '')[:200].replace('\t', ' ').replace('\n', ' ')}
    if record_runtime:
        _lock_update_runtime(row)
    return row


def add_fragments(msa_fa, frag_fa, out_fa, threads, workroot, timeout, mem_mb):
    """Fragments added to msa_fa with --addfragments --keeplength; writes only the fragment rows."""
    wd = tempfile.mkdtemp(prefix='frag.', dir=workroot)
    try:
        frags = read_fasta(frag_fa)
        msa = msa_graph.read_msa(msa_fa)
        ren_m = [('m%d' % i, r) for i, (_, r) in enumerate(msa)]
        clean, dropped = [], 0
        for i, (n, s) in enumerate(frags):
            ms, _ = mask_runs(s.upper(), 1)
            if ms:
                clean.append(('f%d' % i, ms))
            else:
                dropped += 1
        if not clean:
            return {'status': 'none', 'note': 'all fragments empty after N masking'}
        mf = os.path.join(wd, 'msa.fa')
        ff = os.path.join(wd, 'frag.fa')
        of = os.path.join(wd, 'out.fa')
        msa_graph.write_msa(ren_m, mf)
        _write_fa(ff, clean)
        r = mafft_add_fragments(mf, ff, of, threads, wd, timeout, mem_mb)
        if r['status'] != 'ok':
            return r
        got = dict(msa_graph.read_msa(of))
        ncol = len(msa[0][1])
        rows, spelled = [], 0
        for i, (n, s) in enumerate(frags):
            row = got.get('f%d' % i)
            if row is None:
                continue
            if len(row) != ncol:
                return dict(r, status='error', message='fragment row length %d != %d columns' % (len(row), ncol))
            spelled += row.replace('-', '') == mask_runs(s.upper(), 1)[0]
            rows.append((n, row))
        for i, (_, row) in enumerate(ren_m):
            if got.get('m%d' % i) != row:
                return dict(r, status='error', message='--keeplength changed an existing row')
        msa_graph.write_msa(rows, out_fa)
        r.update(n_fragments=len(frags), n_written=len(rows), n_spell_exactly=spelled, n_empty_after_mask=dropped,
                 out=_rel(out_fa))
        return r
    finally:
        shutil.rmtree(wd, ignore_errors=True)


def _job(args):
    method, rd, kw = args
    try:
        return realign_region(method, rd, **kw)
    except Exception as e:  # keep the batch going
        return {'method': method, 'region_id': os.path.basename(rd), 'status': 'crash',
                'note': '%s: %s' % (type(e).__name__, e)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('method', nargs='?', help='method name, or a comma list (region mode)')
    ap.add_argument('regions', nargs='*', help="region ids or dirs; 'all' = every packaged region")
    ap.add_argument('--list', action='store_true', help='print the method table and exit')
    ap.add_argument('--input', help='MSA-only mode: align this FASTA (any records) ...')
    ap.add_argument('--msa-out', help='... and write only the MSA here')
    ap.add_argument('--json', help='MSA-only mode: also write the info dict here')
    ap.add_argument('--threads', type=int, default=2, help='aligner threads (default 2)')
    ap.add_argument('--jobs', type=int, default=1, help='regions aligned at once (default 1; at most 2 advised)')
    ap.add_argument('--timeout', type=float, default=DEFAULT_TIMEOUT, help='seconds per alignment (default 900)')
    ap.add_argument('--mem-mb', type=float, default=DEFAULT_MEM_MB,
                    help='kill the aligner when its process group exceeds this RSS (default 10000)')
    ap.add_argument('--stratum', help='comma list of strata (with all)')
    ap.add_argument('--force', action='store_true', help='re-run regions that already have a result')
    ap.add_argument('--retry-failed', action='store_true', help='re-run regions whose last attempt failed or timed out')
    ap.add_argument('--retry', help='re-run regions whose last status is in this comma list (e.g. memout)')
    ap.add_argument('--fragments', action='store_true', help='also write <id>.fragments.msa.fa (mafft --addfragments)')
    ap.add_argument('--no-dedup', action='store_true', help='align identical sequences separately')
    ap.add_argument('--mask-min-run', type=int, default=1, help='shortest non-ACGT run cut out before alignment')
    ap.add_argument('--merge-blocks', action='store_true', help='pass --merge-blocks to msa_graph (not the default)')
    ap.add_argument('--candidates', help='candidates root (default %s)' % config.CANDIDATES_DIR)
    ap.add_argument('--workdir', help='scratch for aligner runs (default $VNTR_WORK/realign)')
    ap.add_argument('--keep-work', action='store_true', help='MSA-only mode: keep the aligner work dir')
    a = ap.parse_args(argv)

    if a.list:
        for n, m in all_methods().items():
            p = m.params.get('args')
            print('%-16s %-8s %s%s' % (n, m.tool, m.description, ('  [%s]' % ' '.join(p)) if p else ''))
        for mod, e in PLUGIN_ERRORS.items():
            print('# %s failed to import: %s' % (mod, e))
        return 0
    if not a.method:
        ap.error('a method is required (see --list)')

    if a.input or a.msa_out:
        if not (a.input and a.msa_out) or a.regions:
            ap.error('MSA-only mode takes one method, --input FASTA and --msa-out FILE, and no regions')
        wd = a.workdir or os.path.join(config.WORK_DIR, 'realign')
        os.makedirs(wd, exist_ok=True)
        info = align_fasta(a.method, a.input, a.msa_out, threads=a.threads, timeout=a.timeout, mem_mb=a.mem_mb,
                           workdir=wd, dedup=not a.no_dedup, mask_min_run=a.mask_min_run, keep_work=a.keep_work)
        if a.json:
            with open(a.json, 'w') as f:
                json.dump(info, f, indent=1)
        sys.stderr.write('realign: %s %s: %s, %s rows (%s distinct), %s columns, %.1f s%s\n' % (
            a.method, a.input, info['status'], info.get('n_seqs'), info.get('n_distinct'), info.get('columns'),
            info.get('align_s', 0), (' -- ' + info['message']) if info.get('message') else ''))
        return 0 if info['status'] == 'ok' else (4 if info['status'] in ('timeout', 'memout') else 2)

    methods = a.method.split(',')
    for mname in methods:
        get_method(mname)
    if not a.regions:
        ap.error("give region ids, or 'all'")
    if a.regions == ['all']:
        rows = region_list()
        if a.stratum:
            keep = set(a.stratum.split(','))
            rows = [r for r in rows if r.get('stratum') in keep]
        rdirs = [os.path.join(config.REGIONS_DIR, r['region_id']) for r in rows]
    else:
        rdirs = [region_dir(x) for x in a.regions]
    kw = dict(threads=a.threads, timeout=a.timeout, mem_mb=a.mem_mb, force=a.force,
              retry_failed=set(a.retry.split(',')) if a.retry else a.retry_failed,
              fragments=a.fragments, dedup=not a.no_dedup, mask_min_run=a.mask_min_run, cand_root=a.candidates,
              workroot=a.workdir, merge_blocks=a.merge_blocks)
    counts = collections.Counter()
    worst = 0
    for mname in methods:
        jobs = [(mname, rd, kw) for rd in rdirs]
        t0 = time.time()
        if a.jobs <= 1:
            it = map(_job, jobs)
            ex = None
        else:
            ex = concurrent.futures.ProcessPoolExecutor(max_workers=a.jobs)
            it = (f.result() for f in concurrent.futures.as_completed([ex.submit(_job, j) for j in jobs]))
        for row in it:
            counts[(mname, row['status'])] += 1
            if row['status'] not in ('ok', 'skipped'):
                worst = 2
            sys.stderr.write('realign: %-14s %-9s %-11s align %6ss graph %6ss nodes %-6s %s\n' % (
                mname, row['region_id'], row['status'], row.get('align_s', ''), row.get('graph_s', ''),
                row.get('nodes', ''), row.get('note', '')))
            sys.stderr.flush()
        if ex:
            ex.shutdown()
        sys.stderr.write('realign: %s done in %.0f s: %s\n' % (
            mname, time.time() - t0, ', '.join('%s %d' % (s, c) for (mm, s), c in sorted(counts.items()) if mm == mname)))
    return worst


if __name__ == '__main__':
    sys.exit(main())
