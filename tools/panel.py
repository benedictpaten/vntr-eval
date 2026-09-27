#!/usr/bin/env python3
"""panel.py -- the full-panel arm: align every haplotype of the eval graph, not only hap32.

The baseline mc.gfa is not an alignment of the 34 hap32 sequences: Minigraph-Cactus aligned the
whole HPRC v2.1 panel and hap32 keeps 32 sampled rows of it (plus CHM13 and GRCh38). A fair
replacement for MC therefore aligns the full panel and is then projected onto hap32. This tool
provides the pieces:

    union     regions/<id>  ->  distinct sequences (u0001..) + a map TSV of who carries each
    project   a full-panel MSA  ->  an MSA of hap32.fa (hap32 rows, all-gap columns dropped)
    panel-graph  a full-panel MSA  ->  its graph with one path per full-panel haplotype
              (work/panel/<method>/<id>.gfa, for evaluate.py --panel)
    mc-graph  the full Minigraph-Cactus graph between the region's anchors, one path per
              spanning full-panel haplotype  ->  work/panel/mc/<id>.gfa
    stats     per-region panel statistics  ->  results/panel_stats.tsv
    premise   is mc.gfa the full graph restricted to the hap32 walks?  (implied pairwise
              alignments of hap32 paths in both graphs)
    gref-check  one `vg paths -R -A` pass over the eval GBZ: the gref_CHM13 cover never spans an
              anchor, its base copy equals CHM13, GRCh38 spans where hap32.fa says
    check-hprc  every regions/<id>/hprc.fa.gz loads and matches region.json
    import-queries / fetch  the cached full-graph interval queries (gbz-base) this tool reads

Usage (from the repository root):

    python3 tools/panel.py union regions/L012184 work/panel/union/L012184.fa [--map work/panel/union/L012184.map.tsv]
    python3 tools/panel.py project FULL.msa.fa work/panel/union/L012184.map.tsv regions/L012184/hap32.fa \\
                                   candidates/<method>__all/L012184.msa.fa
    python3 tools/panel.py panel-graph FULL.msa.fa work/panel/union/L012184.map.tsv regions/L012184 \\
                                   work/panel/<method>/L012184.gfa [--merge-blocks]
    python3 tools/panel.py union-all [--regions ID,...]          # every region -> work/panel/union/
    python3 tools/panel.py mc-graph [--regions ID,...]           # -> work/panel/mc/<id>.gfa
    python3 tools/panel.py stats [--out results/panel_stats.tsv]
    python3 tools/panel.py premise [--regions ID,...] [--out results/panel_premise.tsv]
    python3 tools/panel.py gref-check [--regions ...]            # ~15 GB RSS for one vg load
    python3 tools/panel.py check-hprc

Full-panel records (regions/<id>/hprc.fa.gz, written by package_regions.py) are the sample
haplotypes of hprc-v2.1-mc-chm13-eval.gref.gbz that pass through both anchor nodes, anonymous
(`hprc#k`: gbz-base reports walks without sample names), preceded by CHM13 and GRCh38 under
their hap32.fa names. HG002 and its parents are absent from that graph (checked from its path
metadata, see gref-check); gref_CHM13, the reference-cover chimera, is not a haplotype and is
excluded (its base copy is subtracted, its fragments never reach an anchor).

The node ids of the hap32 graph are those of the eval graph (vg haplotypes keeps them), so the
full graph's span is located by the same anchor node ids; `premise` checks every mc.gfa node
against the full graph's sequence for that id.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402

REPO = os.path.dirname(TOOLS)
PANEL_DIR = os.path.join(config.WORK_DIR, 'panel')
QUERY_DIR = os.path.join(PANEL_DIR, 'query')
UNION_DIR = os.path.join(PANEL_DIR, 'union')
MC_DIR = os.path.join(PANEL_DIR, 'mc')
HPRC_SAMPLE_HAPLOTYPES = 456          # haplotypes of the 228 HPRC samples in the eval graph
EXCLUDED_SAMPLES = ('HG002', 'HG003', 'HG004', 'NA24385', 'NA24149', 'NA24143')
NON_HAPLOTYPES = ('gref_CHM13',)      # the reference cover: never a panel member
STEP_RE = re.compile(r'([<>])(\d+)')
_RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
T0 = time.time()


class PanelError(RuntimeError):
    pass


def log(msg):
    print('[panel %6.1fs] %s' % (time.time() - T0, msg), file=sys.stderr, flush=True)


def revcomp(s):
    return s.translate(_RC)[::-1]


def _open(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def read_fasta(path):
    """[(name, header after the name, SEQUENCE)] in file order (sequence upper-cased)."""
    out, name, desc, buf = [], None, '', []
    with _open(path) as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if name is not None:
                    out.append((name, desc, ''.join(buf).upper()))
                h = line[1:].split(None, 1)
                name, desc, buf = (h[0] if h else ''), (h[1] if len(h) > 1 else ''), []
            elif line:
                buf.append(line.strip())
    if name is not None:
        out.append((name, desc, ''.join(buf).upper()))
    return out


def write_fasta(path, recs):
    """recs: [(name, description, sequence)]; one line per sequence; .gz by extension."""
    d = os.path.dirname(os.path.abspath(path))
    os.makedirs(d, exist_ok=True)
    tmp = path + '.tmp%d' % os.getpid() + ('.gz' if path.endswith('.gz') else '')
    with _open(tmp, 'wt') as f:
        for n, desc, s in recs:
            f.write('>%s%s\n%s\n' % (n, (' ' + desc) if desc else '', s))
    os.replace(tmp, path)


def sample_of(name):
    return name.split('#', 1)[0]


def assert_no_excluded(names, where):
    bad = [n for n in names if sample_of(n) in EXCLUDED_SAMPLES + NON_HAPLOTYPES]
    if bad:
        raise PanelError('%s holds excluded haplotypes (HG002 family or gref_CHM13): %s'
                         % (where, ', '.join(bad[:5])))


def region_dir(x):
    """A region id or a region directory -> (region_dir, region_id)."""
    if os.path.isdir(x):
        d = os.path.abspath(x)
    else:
        d = os.path.join(config.REGIONS_DIR, x)
        if not os.path.isdir(d):
            raise PanelError('no region package %s' % x)
    return d, os.path.basename(d.rstrip('/'))


def all_region_ids():
    p = os.path.join(config.REGIONS_DIR, 'regions.tsv')
    with open(p) as f:
        return [r['region_id'] for r in csv.DictReader(f, delimiter='\t')]


def pick_regions(arg):
    if not arg:
        return all_region_ids()
    out = []
    for x in arg.split(','):
        x = x.strip()
        if x.endswith('.tsv') and os.path.exists(x):
            with open(x) as f:
                out.extend(r['region_id'] for r in csv.DictReader(f, delimiter='\t'))
        elif x:
            out.append(x)
    return out


def load_region(x):
    rd, rid = region_dir(x)
    rj = json.load(open(os.path.join(rd, 'region.json')))
    hap = read_fasta(os.path.join(rd, 'hap32.fa'))
    hp = os.path.join(rd, 'hprc.fa.gz')
    hprc = read_fasta(hp) if os.path.exists(hp) else []
    return rd, rid, rj, hap, hprc


def is_reference_name(n):
    return n.startswith('CHM13#') or n.startswith('GRCh38#')


# ============================================================================ union / project

def union(rdir):
    """Deduplicated sequences of hap32.fa (every spanning hap32 path, so a sampled recombinant
    whose junction falls inside the region is covered) and hprc.fa.gz (every full-panel
    haplotype spanning the anchors, CHM13 and GRCh38 included).

    Returns (records, map_rows): records [(uid, description, seq)]; map_rows dicts with id,
    length, weight (full-panel haplotypes with this sequence, CHM13 and GRCh38 counted once
    each), n_hap32, n_panel_samples, n_bp_N and members (hap32 names, then full-panel names).
    Order: CHM13's sequence first (u0001), then decreasing weight, then n_hap32, length and
    sequence -- deterministic for given inputs."""
    rd, rid, rj, hap, hprc = load_region(rdir)
    if not hap:
        raise PanelError('%s: empty hap32.fa' % rid)
    names = [n for n, _, _ in hap] + [n for n, _, _ in hprc]
    assert_no_excluded(names, rid)
    chm = [s for n, _, s in hap if n.startswith('CHM13#')]
    if len(chm) != 1:
        raise PanelError('%s: need one CHM13 record in hap32.fa' % rid)
    groups = collections.OrderedDict()

    def g(s):
        if s not in groups:
            groups[s] = {'hap32': [], 'panel': [], 'samples': 0}
        return groups[s]
    for n, _, s in hap:
        g(s)['hap32'].append(n)
    seen = set()
    for n, _, s in hprc:
        if n in seen:
            raise PanelError('%s: duplicate name %s in hprc.fa.gz' % (rid, n))
        seen.add(n)
        e = g(s)
        e['panel'].append(n)
        if not is_reference_name(n):
            e['samples'] += 1
    order = sorted(groups.items(), key=lambda kv: (kv[0] != chm[0], -len(kv[1]['panel']),
                                                  -len(kv[1]['hap32']), len(kv[0]), kv[0]))
    w = max(4, len(str(len(order))))
    recs, rows = [], []
    for i, (s, e) in enumerate(order, start=1):
        uid = 'u%0*d' % (w, i)
        members = e['hap32'] + [n for n in e['panel'] if n not in e['hap32']]
        nN = s.count('N')
        recs.append((uid, 'len=%d weight=%d n_hap32=%d' % (len(s), len(e['panel']), len(e['hap32'])), s))
        rows.append({'id': uid, 'length': len(s), 'weight': len(e['panel']), 'n_hap32': len(e['hap32']),
                     'n_panel_samples': e['samples'], 'n_bp_N': nN, 'members': ','.join(members)})
    return recs, rows


MAP_COLS = ['id', 'length', 'weight', 'n_hap32', 'n_panel_samples', 'n_bp_N', 'members']


def write_map(path, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + '.tmp%d' % os.getpid()
    with open(tmp, 'w') as f:
        f.write('\t'.join(MAP_COLS) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in MAP_COLS) + '\n')
    os.replace(tmp, path)


def read_map(path):
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    for r in rows:
        r['members'] = r['members'].split(',') if r['members'] else []
    return rows


def union_to_files(rdir, out_fa, map_tsv=None):
    recs, rows = union(rdir)
    if map_tsv is None:
        base = out_fa[:-3] if out_fa.endswith('.gz') else out_fa
        base = base[:-3] if base.endswith('.fa') else base
        map_tsv = base + '.map.tsv'
    write_fasta(out_fa, recs)
    write_map(map_tsv, rows)
    return {'distinct': len(recs), 'fasta': out_fa, 'map': map_tsv,
            'hap32_only': sum(1 for r in rows if r['weight'] == 0),
            'total_bp': sum(r['length'] for r in rows), 'max_bp': max(r['length'] for r in rows)}


def project(full_msa, map_tsv, hap32_fa, out_msa):
    """Restrict a full-panel MSA to hap32.fa: one row per hap32 name (the row of the distinct
    sequence it belongs to), all-gap columns dropped; asserts that every ungapped row equals its
    hap32.fa sequence. MSA rows may be named by distinct id (u0001) or by any member name; a
    consensus row (abPOA) is ignored. Returns a stats dict."""
    import msa_graph
    rows = msa_graph.read_msa(full_msa)
    rows = [(n, s) for n, s in rows if n not in msa_graph.CONSENSUS_NAMES]
    if not rows:
        raise PanelError('%s: no rows' % full_msa)
    lens = set(len(s) for _, s in rows)
    if len(lens) != 1:
        raise PanelError('%s: rows have different lengths %s' % (full_msa, sorted(lens)[:5]))
    mp = read_map(map_tsv)
    uid_of = {}
    for r in mp:
        uid_of[r['id']] = r['id']
        for m in r['members']:
            uid_of[m] = r['id']
    row_of = {}
    unknown = []
    for n, s in rows:
        u = uid_of.get(n)
        if u is None:
            unknown.append(n)
            continue
        row_of.setdefault(u, s)          # first row of a distinct sequence wins
    if unknown:
        raise PanelError('%d MSA rows are not in the map: %s' % (len(unknown), ', '.join(unknown[:5])))
    hap = read_fasta(hap32_fa)
    out = []
    missing = []
    for n, _, s in hap:
        u = uid_of.get(n)
        if u is None or u not in row_of:
            missing.append(n)
            continue
        out.append((n, row_of[u]))
    if missing:
        raise PanelError('%d hap32.fa sequences have no row in %s: %s' % (len(missing), full_msa,
                                                                          ', '.join(missing[:5])))
    ncol = len(out[0][1])
    keep = [c for c in range(ncol) if any(r[c] != '-' for _, r in out)]
    if len(keep) < ncol:
        out = [(n, ''.join(r[c] for c in keep)) for n, r in out]
    want = {n: s for n, _, s in hap}
    bad = [n for n, r in out if r.replace('-', '') != want[n]]
    if bad:
        raise PanelError('%d projected rows do not spell their hap32.fa sequence: %s' % (len(bad), ', '.join(bad[:5])))
    os.makedirs(os.path.dirname(os.path.abspath(out_msa)), exist_ok=True)
    tmp = out_msa + '.tmp%d' % os.getpid()
    msa_graph.write_msa(out, tmp)
    os.replace(tmp, out_msa)
    return {'rows_in': len(rows), 'distinct_rows_used': len(set(uid_of[n] for n, _, _ in hap)),
            'rows_out': len(out), 'columns_in': ncol, 'columns_out': len(keep),
            'all_gap_columns_dropped': ncol - len(keep)}


def panel_graph(full_msa, map_tsv, region, out_gfa, merge_blocks=False, engine='auto', block_max=50):
    """The full-panel graph of a full-panel MSA: the graph is induced from the distinct rows with
    tools/msa_graph.py (same induction and compaction as the hap32 arm), then written with one
    P line per regions/<id>/hprc.fa.gz record (name kept) following its distinct sequence's
    path, and every path's spelling is checked against hprc.fa.gz. The map must be the one
    `union` writes for this region (it is recomputed and compared). This is the graph that
    `evaluate.py --panel regions/<id>/hprc.fa.gz` judges (work/panel/<method>/<id>.gfa)."""
    import tempfile
    import msa_graph
    rd, rid, rj, hap, hprc = load_region(region)
    recs, rows = union(rd)
    mp = read_map(map_tsv)
    if [(r['id'], int(r['length']), r['members']) for r in mp] != \
            [(r['id'], r['length'], r['members'].split(',')) for r in rows]:
        raise PanelError('%s does not match the union of %s (re-run `panel.py union`)' % (map_tsv, rid))
    uid_of = {}
    for r in mp:
        uid_of[r['id']] = r['id']
        for m in r['members']:
            uid_of[m] = r['id']
    msa = [(n, s) for n, s in msa_graph.read_msa(full_msa) if n not in msa_graph.CONSENSUS_NAMES]
    got = collections.OrderedDict()
    for n, s in msa:
        u = uid_of.get(n)
        if u is None:
            raise PanelError('MSA row %s is not in the map' % n)
        got.setdefault(u, s)
    need = [r['id'] for r in rows if r['weight'] > 0]            # ids that carry full-panel records
    missing = [u for u in need if u not in got]
    if missing:
        raise PanelError('%d distinct full-panel sequences have no MSA row: %s' % (len(missing), ', '.join(missing[:5])))
    seq_of = {u: s for u, _, s in recs}
    os.makedirs(os.path.dirname(os.path.abspath(out_gfa)) or '.', exist_ok=True)   # may not exist yet
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(out_gfa)) or '.') as td:
        ufa, umsa, ugfa = (os.path.join(td, x) for x in ('u.fa', 'u.msa.fa', 'u.gfa'))
        write_fasta(ufa, [(u, '', seq_of[u]) for u in got])
        msa_graph.write_msa(list(got.items()), umsa)
        stats = msa_graph.msa_to_gfa(umsa, ufa, ugfa, engine=engine, merge_blocks=merge_blocks,
                                     block_max=block_max, workdir=td)
        lines, upath = [], {}
        with open(ugfa) as f:
            for line in f:
                if line.startswith('P\t'):
                    x = line.rstrip('\n').split('\t')
                    upath[x[1]] = x[2]
                else:
                    lines.append(line)
    tmp = out_gfa + '.tmp%d' % os.getpid()
    with open(tmp, 'w') as f:
        f.writelines(lines)
        for n, _, _ in hprc:
            f.write('P\t%s\t%s\t*\n' % (n, upath[uid_of[n]]))
    seqs, _, paths = read_gfa_paths(tmp)
    want = {n: s for n, _, s in hprc}
    bad = [n for n, w in paths.items()
           if ''.join(seqs[x] if x > 0 else revcomp(seqs[-x]) for x in w) != want.get(n)]
    if bad or len(paths) != len(want):
        os.remove(tmp)
        raise PanelError('%s: panel graph fails its spelling check (%d bad, %d/%d paths)'
                         % (rid, len(bad), len(paths), len(want)))
    os.replace(tmp, out_gfa)
    return {'region_id': rid, 'gfa': out_gfa, 'panel_paths': len(paths), 'distinct_rows': len(got),
            'nodes': stats['nodes'], 'edges': stats.get('edges'), 'columns': stats['columns'],
            'mode': stats['mode'], 'engine': stats['engine']}


# ============================================================================ full-graph query

def find_query(rid, query_dir=None):
    """The cached full-graph interval query of a region (gzipped GFA from gbz-base --distinct)."""
    cands = []
    if query_dir:
        cands += [os.path.join(query_dir, rid + '.gfa.gz'),
                  os.path.join(query_dir, rid, 'hprc.interval.gfa.gz')]
    cands += [os.path.join(QUERY_DIR, rid + '.gfa.gz'),
              os.path.join(config.WORK_DIR, 'package', rid, 'hprc.interval.gfa.gz')]
    for c in cands:
        if os.path.exists(c) and os.path.getsize(c) > 0:
            return c
    return None


def query_params(rj):
    """The interval query package_regions.py ran (and that hprc.fa.gz was cut from)."""
    a0, b0 = rj['span_start'] - 1, rj['span_end']
    return {'contig': rj['contig'], 'a0': a0, 'b0': b0, 'context': max(10000, b0 - a0)}


def query_cmd(rj):
    p = query_params(rj)
    gbz = config.data_paths(p['contig'])['hprc_gbz']
    return [config.GBZ_BASE, 'query', gbz, '--sample', 'CHM13', '--contig', p['contig'],
            '--interval', '%d..%d' % (p['a0'], p['b0']), '--context', str(p['context']), '--distinct']


def fetch_query(rid, force=False):
    """Run the gbz-base interval query on the eval GBZ (~45-65 s, ~5-11 GB RSS)."""
    rd, rid = region_dir(rid)
    rj = json.load(open(os.path.join(rd, 'region.json')))
    out = os.path.join(QUERY_DIR, rid + '.gfa.gz')
    if os.path.exists(out) and not force:
        return out
    os.makedirs(QUERY_DIR, exist_ok=True)
    cmd = query_cmd(rj)
    tmp = out + '.tmp%d.gz' % os.getpid()
    t = time.time()
    with open(os.path.join(QUERY_DIR, rid + '.log'), 'w') as fe:
        fe.write('$ %s\n' % ' '.join(cmd))
        fe.flush()
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=fe, env=config.tool_env())
        with gzip.open(tmp, 'wb', compresslevel=3) as fo:
            shutil.copyfileobj(p.stdout, fo, 1 << 20)
        p.wait()
    if p.returncode != 0:
        os.remove(tmp)
        raise PanelError('%s: gbz-base failed (see %s.log)' % (rid, os.path.join(QUERY_DIR, rid)))
    os.replace(tmp, out)
    log('%s: queried in %.0fs' % (rid, time.time() - t))
    return out


def import_queries(src, rids):
    """Copy package_regions.py's cached queries (<src>/<id>/hprc.interval.gfa.gz) into
    work/panel/query/, after checking that each was run with this region's parameters."""
    os.makedirs(QUERY_DIR, exist_ok=True)
    n = 0
    for rid in rids:
        rd, rid = region_dir(rid)
        rj = json.load(open(os.path.join(rd, 'region.json')))
        gz = os.path.join(src, rid, 'hprc.interval.gfa.gz')
        mj = os.path.join(src, rid, 'hprc.interval.json')
        if not (os.path.exists(gz) and os.path.exists(mj)):
            log('%s: no cached query in %s' % (rid, src))
            continue
        meta = json.load(open(mj))
        want = query_params(rj)
        if any(meta.get(k) != v for k, v in want.items()) or meta.get('returncode') != 0:
            log('%s: cached query parameters differ (%s), skipped' % (rid, {k: meta.get(k) for k in want}))
            continue
        dst = os.path.join(QUERY_DIR, rid + '.gfa.gz')
        if not os.path.exists(dst) or os.path.getsize(dst) != os.path.getsize(gz):
            shutil.copyfile(gz, dst + '.tmp')
            os.replace(dst + '.tmp', dst)
        with open(os.path.join(QUERY_DIR, rid + '.json'), 'w') as f:
            json.dump(dict(meta, imported_from=gz), f, indent=1)
        n += 1
    return n


def parse_steps(w):
    return tuple(int(n) if o == '>' else -int(n) for o, n in STEP_RE.findall(w))


def steps_str(steps):
    return ''.join(('>' if s > 0 else '<') + str(abs(s)) for s in steps)


def rev_steps(steps):
    return tuple(-s for s in reversed(steps))


def parse_anchor(s):
    s = str(s).strip()
    m = re.match(r'^([<>])(\d+)$', s) or re.match(r'^(\d+)([+-])$', s)
    if not m:
        raise PanelError('cannot parse anchor %r' % s)
    a, b = m.groups()
    if a in '<>':
        return int(b) if a == '>' else -int(b)
    return int(a) if b == '+' else -int(a)


def canon_edge(a, b):
    """signed-id edge a->b equals (-b)->(-a); keep the smaller representation"""
    return min((a, b), (-b, -a))


class QueryGraph:
    """A gbz-base GFA: S (int id -> sequence), L (canonical signed edges), W (walks + WT weight)."""

    def __init__(self, path):
        self.seqs, self.edges, self.walks = {}, set(), []
        with _open(path) as f:
            for line in f:
                t = line[0]
                if t == 'S':
                    x = line.rstrip('\n').split('\t')
                    self.seqs[int(x[1])] = x[2].upper()
                elif t == 'L':
                    x = line.rstrip('\n').split('\t')
                    a = int(x[1]) * (1 if x[2] == '+' else -1)
                    b = int(x[3]) * (1 if x[4] == '+' else -1)
                    self.edges.add(canon_edge(a, b))
                elif t == 'W':
                    x = line.rstrip('\n').split('\t')
                    wt = 1
                    for tag in x[7:]:
                        if tag.startswith('WT:i:'):
                            wt = int(tag[5:])
                    self.walks.append({'sample': x[1], 'hap': x[2], 'contig': x[3],
                                       'steps': parse_steps(x[6]), 'weight': wt})
        self._rc = {}

    def nseq(self, s):
        if s > 0:
            return self.seqs[s]
        r = self._rc.get(s)
        if r is None:
            r = self._rc[s] = revcomp(self.seqs[-s])
        return r

    def spell(self, steps):
        return ''.join(self.nseq(s) for s in steps)


def cut_at_anchors(walks, L, R):
    """Cut every walk at the anchors, exactly as package_regions.cut_at_anchors (which wrote
    hprc.fa.gz). Returns (Counter {L..R steps in reference orientation: weight}, Counter of
    the other anchor-touching walks by kind)."""
    spans, other = collections.Counter(), collections.Counter()
    aL, aR = abs(L), abs(R)
    for w in walks:
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
                spans[rev_steps(st[j:i + 1])] += wt
            else:
                other['anchors_out_of_order'] += wt
        elif len(iL) > 1 or len(iR) > 1:
            other['anchor_visited_twice'] += wt
        elif iL:
            other['left_anchor_only'] += wt
        else:
            other['right_anchor_only'] += wt
    return spans, other


def read_gfa_paths(path):
    """(seqs {int id: seq}, edges set, paths OrderedDict name -> steps) of a GFA with P lines."""
    seqs, edges, paths = {}, set(), collections.OrderedDict()
    with _open(path) as f:
        for line in f:
            t = line[0]
            if t == 'S':
                x = line.rstrip('\n').split('\t')
                seqs[int(x[1])] = x[2].upper()
            elif t == 'L':
                x = line.rstrip('\n').split('\t')
                edges.add(canon_edge(int(x[1]) * (1 if x[2] == '+' else -1),
                                     int(x[3]) * (1 if x[4] == '+' else -1)))
            elif t == 'P':
                x = line.rstrip('\n').split('\t')
                paths[x[1]] = tuple(int(s[:-1]) * (1 if s[-1] == '+' else -1) for s in x[2].split(',') if s)
    return seqs, edges, paths


class Panel:
    """Everything about one region's full panel, reconstructed from the cached interval query
    and checked against the package: the spanning walks, the reference subtraction, and the
    assignment of walks to the anonymous hprc#k records."""

    def __init__(self, rid, query_dir=None, need_graph=True):
        self.rd, self.rid, self.rj, self.hap, self.hprc = load_region(rid)
        rj = self.rj
        if rj.get('hprc', {}).get('status') != 'ok' or not self.hprc:
            raise PanelError('%s: no full-panel records (hprc status %s)' % (self.rid, rj.get('hprc', {}).get('status')))
        self.L, self.R = parse_anchor(rj['anchor_left']), parse_anchor(rj['anchor_right'])
        mc = os.path.join(self.rd, 'mc.gfa')
        self.mc_seqs, self.mc_edges, self.mc_paths = read_gfa_paths(mc)
        qp = find_query(self.rid, query_dir)
        if qp is None:
            raise PanelError('%s: no cached full-graph query; run `panel.py import-queries` or `panel.py fetch`'
                             % self.rid)
        self.query_path = qp
        self.q = QueryGraph(qp)
        self._build()

    def _build(self):
        q, rid = self.q, self.rid
        for node, want in ((self.L, self.rj['anchor_left_seq']), (self.R, self.rj['anchor_right_seq'])):
            if abs(node) not in q.seqs or q.nseq(node) != want:
                raise PanelError('%s: anchor %d absent or different in the full graph' % (rid, node))
        spans, other = cut_at_anchors(q.walks, self.L, self.R)
        self.other = other
        self.n_segments = sum(spans.values())
        ref = [w for w in q.walks if w['sample'] == 'CHM13']
        if len(ref) != 1:
            raise PanelError('%s: %d CHM13 walks in the query' % (rid, len(ref)))
        rs, _ = cut_at_anchors(ref, self.L, self.R)
        if len(rs) != 1:
            raise PanelError('%s: the CHM13 walk does not span the anchors once' % rid)
        self.chm13_walk = next(iter(rs))
        chm_name = [n for n, _, _ in self.hap if n.startswith('CHM13#')][0]
        if self.mc_paths.get(chm_name) != self.chm13_walk:
            raise PanelError('%s: CHM13 walk differs between mc.gfa and the full graph' % rid)
        self.chm13_walk_weight = spans[self.chm13_walk]
        if spans[self.chm13_walk] < 2:
            raise PanelError('%s: the CHM13 walk has weight %d (< 2: CHM13 and gref_CHM13)' % (rid, spans[self.chm13_walk]))
        spans[self.chm13_walk] -= 2           # CHM13 itself and its gref_CHM13 copy
        self.ref_walks = [(chm_name, self.chm13_walk)]
        self.grch38 = []
        for n, _, _ in self.hap:
            if n.startswith('GRCh38#'):
                w = self.mc_paths.get(n)
                if w is None or spans[w] < 1:
                    raise PanelError('%s: the GRCh38 walk of %s is not among the full-graph walks' % (rid, n))
                spans[w] -= 1
                self.ref_walks.append((n, w))
                self.grch38.append(n)
        spans = +spans
        self.sample_walks = spans            # {steps: weight} of the sample haplotypes
        # assign walks to the hprc.fa.gz records: group by sequence, check the multiset
        by_seq = collections.defaultdict(list)
        for w, wt in spans.items():
            by_seq[q.spell(w)].append((w, wt))
        want = collections.Counter(s for n, _, s in self.hprc if not is_reference_name(n))
        got = collections.Counter({s: sum(wt for _, wt in v) for s, v in by_seq.items()})
        if want != got:
            d1 = sum((want - got).values())
            d2 = sum((got - want).values())
            raise PanelError('%s: reconstructed sample haplotypes differ from hprc.fa.gz (%d missing, %d extra)'
                             % (rid, d1, d2))
        queues = {}
        for s, v in by_seq.items():
            v.sort(key=lambda x: (-x[1], x[0]))
            queues[s] = [w for w, wt in v for _ in range(wt)]
        self.walk_of = collections.OrderedDict()
        for n, _, s in self.hprc:
            if is_reference_name(n):
                w = dict(self.ref_walks).get(n)
                if w is None:
                    raise PanelError('%s: reference record %s of hprc.fa.gz is not in hap32.fa' % (rid, n))
                self.walk_of[n] = w
            else:
                self.walk_of[n] = queues[s].pop(0)
        for n, _, s in self.hprc:
            if q.spell(self.walk_of[n]) != s:
                raise PanelError('%s: walk of %s does not spell its sequence' % (rid, n))

    # ------------------------------------------------------------------ outputs

    def write_gfa(self, out):
        """The full MC graph between the anchors: nodes of the spanning walks (samples, CHM13,
        GRCh38), the query's edges among them, one P line per hprc.fa.gz record (spelling it)."""
        q = self.q
        nodes = set()
        for w in self.walk_of.values():
            nodes.update(abs(s) for s in w)
        edges = set()
        for w in set(self.walk_of.values()):
            for a, b in zip(w, w[1:]):
                edges.add(canon_edge(a, b))
        for a, b in q.edges:
            if abs(a) in nodes and abs(b) in nodes:
                edges.add((a, b))
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        tmp = out + '.tmp%d' % os.getpid() + ('.gz' if out.endswith('.gz') else '')
        with _open(tmp, 'wt') as f:
            f.write('H\tVN:Z:1.0\n')
            for n in sorted(nodes):
                f.write('S\t%d\t%s\n' % (n, q.seqs[n]))
            for a, b in sorted(edges, key=lambda e: (abs(e[0]), abs(e[1]), e)):
                f.write('L\t%d\t%s\t%d\t%s\t0M\n' % (abs(a), '+' if a > 0 else '-', abs(b), '+' if b > 0 else '-'))
            for n, w in self.walk_of.items():
                f.write('P\t%s\t%s\t*\n' % (n, ','.join('%d%s' % (abs(s), '+' if s > 0 else '-') for s in w)))
        # re-read and check every path's spelling
        seqs, _, paths = read_gfa_paths(tmp)
        want = {n: s for n, _, s in self.hprc}
        comp = {}
        bad = []
        for n, w in paths.items():
            s = ''.join(seqs[x] if x > 0 else comp.setdefault(x, revcomp(seqs[-x])) for x in w)
            if s != want.get(n):
                bad.append(n)
        if bad or len(paths) != len(want):
            os.remove(tmp)
            raise PanelError('%s: written graph fails the spelling check (%d bad, %d/%d paths)'
                             % (self.rid, len(bad), len(paths), len(want)))
        os.replace(tmp, out)
        lens = [len(q.seqs[n]) for n in nodes]
        return {'region_id': self.rid, 'gfa': out, 'paths': len(paths), 'nodes': len(nodes),
                'edges': len(edges), 'node_bp': sum(lens),
                'distinct_walks': len(set(self.walk_of.values())),
                'distinct_seqs': len(set(want.values())),
                'mean_node_bp': round(sum(lens) / len(lens), 2),
                'nodes_per_kb': round(1000.0 * len(nodes) / len(q.spell(self.chm13_walk)), 1),
                'bytes': os.path.getsize(out)}


# ============================================================================ walk analysis

KIND_PRIORITY = {'sample': 0, 'partial': 1, 'CHM13': 2, 'GRCh38': 3}


def classify_hap32_walk(w, cands, pos_index):
    """How a hap32 walk relates to the full-panel walks. cands: [(kind, walk)] with kind 'sample'
    (a full-panel haplotype's anchor-to-anchor walk), 'CHM13' / 'GRCh38' (the reference walks,
    which are not panel members), 'partial' (a full-panel walk that does not run anchor to anchor,
    either orientation). Returns (label, piece kinds, novel steps):
      verbatim           the walk of some full-panel haplotype
      copy_of_CHM13 / copy_of_GRCh38   the reference walk, which no full-panel haplotype takes
      junction           two pieces of full-panel haplotype walks (a recombination inside the span)
      mosaic_<k>         k > 2 such pieces
      mixed(<kinds>)     pieces that need a reference or a non-spanning walk
    Pieces come from greedy longest extension (optimal for the number of pieces when covering a
    string by substrings of a set); a tie prefers sample > partial > reference."""
    for kind, cw in cands:
        if cw == w and kind == 'sample':
            return 'verbatim', ['sample'], 0
    for kind, cw in cands:
        if cw == w and kind in ('CHM13', 'GRCh38'):
            return 'copy_of_' + kind, [kind], 0
    i, n, novel, kinds = 0, len(w), 0, []
    while i < n:
        best, bk = 0, None
        for (k, p) in pos_index.get(w[i], ()):
            kind, pw = cands[k]
            L = 0
            while i + L < n and p + L < len(pw) and pw[p + L] == w[i + L]:
                L += 1
            if L > best or (L == best and L and KIND_PRIORITY[kind] < KIND_PRIORITY[bk]):
                best, bk = L, kind
        if best == 0:
            novel += 1
            kinds.append('novel')
            i += 1
        else:
            kinds.append(bk)
            i += best
    if all(k == 'sample' for k in kinds):
        return ('junction' if len(kinds) == 2 else 'mosaic_%d' % len(kinds)), kinds, novel
    return 'mixed(%s)' % '+'.join(kinds), kinds, novel


def panel_index(cands):
    idx = collections.defaultdict(list)
    for k, (_, w) in enumerate(cands):
        for p, s in enumerate(w):
            idx[s].append((k, p))
    return idx


def walk_candidates(P):
    """[(kind, walk)] for classify_hap32_walk: the sample walks, the reference walks, and every
    query walk that does not run anchor to anchor (both orientations)."""
    cands = [('sample', w) for w in P.sample_walks]
    cands += [('CHM13' if n.startswith('CHM13') else 'GRCh38', w) for n, w in P.ref_walks]
    for x in P.q.walks:
        sp, _ = cut_at_anchors([x], P.L, P.R)
        if not sp:
            cands.append(('partial', x['steps']))
            cands.append(('partial', rev_steps(x['steps'])))
    return cands


# ============================================================================ stats

STATS_COLS = [
    'region_id', 'stratum', 'contig', 'span_bp', 'chm13_bp',
    'n_hap32', 'n_hap32_sampled', 'n_hap32_distinct',
    'hprc_spanning', 'hprc_not_spanning', 'hprc_beyond_456',
    'walks_left_anchor_only', 'walks_right_anchor_only', 'walks_anchor_other', 'reference_removed',
    'chm13_walk_weight_before', 'gref_fragments_in_span', 'gref_fragments_with_anchor',
    'gref_copy_equals_chm13', 'grch38_spanning_in_graph',
    'panel_distinct_seqs', 'panel_distinct_walks', 'union_distinct_seqs',
    'hap32_sampled_verbatim_seq', 'hap32_sampled_verbatim_walk', 'hap32_not_verbatim',
    'hap32_not_verbatim_kinds', 'hap32_not_verbatim_why',
    'distinct_bp_max', 'distinct_bp_min', 'distinct_bp_total', 'distinct_N_bp', 'distinct_seqs_with_N',
    'panel_haplotypes_with_N', 'hap32_sampled_top_seq_share',
]


def region_stats(rid, query_dir=None, gref=None):
    P = Panel(rid, query_dir)
    rj, q = P.rj, P.q
    sampled = [(n, s) for n, _, s in P.hap if not is_reference_name(n)]
    panel_seqs = collections.Counter(s for n, _, s in P.hprc if not is_reference_name(n))
    walk_lookup = set(P.sample_walks)
    cands = idx = None
    verb_seq = verb_walk = 0
    why = []
    kinds = collections.Counter()
    for n, s in sampled:
        w = P.mc_paths[n]
        in_seq = s in panel_seqs
        in_walk = w in walk_lookup
        verb_seq += in_seq
        verb_walk += in_walk
        if not in_walk:
            if idx is None:
                cands = walk_candidates(P)
                idx = panel_index(cands)
            label, pk, novel = classify_hap32_walk(w, cands, idx)
            if in_seq:
                label = 'same_seq_other_walk'
            kinds[label.split('(')[0]] += 1
            why.append('%s:%s%s' % (n, label, '(%d novel steps)' % novel if novel else ''))
    union_seqs = set(s for _, _, s in P.hap) | set(s for _, _, s in P.hprc)
    dl = [len(s) for s in union_seqs]
    nN = [s.count('N') for s in union_seqs]
    other = P.other
    n_span = sum(P.sample_walks.values())
    top = panel_seqs.most_common(1)[0][1] if panel_seqs else 0
    g = (gref or {}).get(P.rid, {})
    row = {
        'region_id': P.rid, 'stratum': rj['stratum'], 'contig': rj['contig'],
        'span_bp': rj['span_end'] - rj['span_start'] + 1,
        'chm13_bp': len(q.spell(P.chm13_walk)),
        'n_hap32': len(P.hap), 'n_hap32_sampled': len(sampled),
        'n_hap32_distinct': len(set(s for _, _, s in P.hap)),
        'hprc_spanning': n_span,
        'hprc_not_spanning': max(0, HPRC_SAMPLE_HAPLOTYPES - n_span),
        'hprc_beyond_456': max(0, n_span - HPRC_SAMPLE_HAPLOTYPES),
        'walks_left_anchor_only': other.get('left_anchor_only', 0),
        'walks_right_anchor_only': other.get('right_anchor_only', 0),
        'walks_anchor_other': other.get('anchor_visited_twice', 0) + other.get('anchors_out_of_order', 0),
        'reference_removed': ';'.join(['CHM13', 'gref_CHM13'] + P.grch38),
        'chm13_walk_weight_before': P.chm13_walk_weight,
        'gref_fragments_in_span': g.get('fragments_in_span', ''),
        'gref_fragments_with_anchor': g.get('fragments_with_anchor_or_ref_node', ''),
        'gref_copy_equals_chm13': g.get('gref_copy_equals_chm13', ''),
        'grch38_spanning_in_graph': g.get('grch38_spanning', ''),
        'panel_distinct_seqs': len(panel_seqs),
        'panel_distinct_walks': len(P.sample_walks),
        'union_distinct_seqs': len(union_seqs),
        'hap32_sampled_verbatim_seq': verb_seq,
        'hap32_sampled_verbatim_walk': verb_walk,
        'hap32_not_verbatim': len(sampled) - verb_seq,
        'hap32_not_verbatim_kinds': ';'.join('%s=%d' % kv for kv in sorted(kinds.items())) or '.',
        'hap32_not_verbatim_why': ';'.join(why) if why else '.',
        'distinct_bp_max': max(dl), 'distinct_bp_min': min(dl), 'distinct_bp_total': sum(dl),
        'distinct_N_bp': sum(nN), 'distinct_seqs_with_N': sum(1 for x in nN if x),
        'panel_haplotypes_with_N': sum(c for s, c in panel_seqs.items() if 'N' in s),
        'hap32_sampled_top_seq_share': round(top / max(1, n_span), 3),
    }
    return row


# ============================================================================ premise check

def implied_matches(wa, wb, lens):
    """The base pairs two walks align through shared node visits, as evaluate.segments() pairs
    them (greedy, increasing in both walks): [(offset in a, offset in b, node length)]."""
    posb = {}
    for j, h in enumerate(wb):
        posb.setdefault(h, j)
    ca = [0]
    for h in wa:
        ca.append(ca[-1] + lens[abs(h)])
    cb = [0]
    for h in wb:
        cb.append(cb[-1] + lens[abs(h)])
    out = []
    last = -1
    for i, h in enumerate(wa):
        j = posb.get(h)
        if j is None or j <= last:
            continue
        out.append((ca[i], cb[j], lens[abs(h)]))
        last = j
    return out, ca, cb


def implied_cost(sa, sb, wa, wb, lens, C, memo):
    """graph-implied unit cost: optimal alignment of every stretch between shared node visits"""
    m, ca, cb = implied_matches(wa, wb, lens)
    cost = 0
    pa = pb = 0
    for (a0, b0, L) in m + [(len(sa), len(sb), 0)]:
        x, y = sa[pa:a0], sb[pb:b0]
        if x or y:
            key = (x, y) if x <= y else (y, x)
            d = memo.get(key)
            if d is None:
                d = memo[key] = C.ed(key[0].encode(), key[1].encode())
            cost += d
        pa, pb = a0 + L, b0 + L
    return cost, m


def premise_region(rid, query_dir=None, max_pairs=None):
    """Is mc.gfa the full graph restricted to the hap32 paths? Node ids and sequences, edges,
    each hap32 walk among the full-panel walks, and the implied pairwise alignment of every pair
    of hap32 paths in mc.gfa against that of their full-panel counterparts in the full graph."""
    import evaluate
    C = evaluate.clib()
    P = Panel(rid, query_dir)
    q = P.q
    node_missing = sum(1 for n in P.mc_seqs if n not in q.seqs)
    node_diff = sum(1 for n, s in P.mc_seqs.items() if n in q.seqs and q.seqs[n] != s)
    edge_missing = sum(1 for e in P.mc_edges if e not in q.edges)
    # the counterpart of each hap32 path in the full panel: the same walk when a full-panel
    # haplotype (or reference path) takes it, else the most frequent walk spelling its sequence
    ref_walks = dict(P.ref_walks)
    by_seq = collections.defaultdict(list)
    for w, wt in P.sample_walks.items():
        by_seq[q.spell(w)].append((wt, w))
    names = [n for n, _, _ in P.hap]
    mcw = {n: P.mc_paths[n] for n in names}
    counterpart, kind = {}, {}
    for n, _, s in P.hap:
        w = mcw[n]
        if n in ref_walks:
            counterpart[n], kind[n] = ref_walks[n], 'reference'
        elif w in P.sample_walks:
            counterpart[n], kind[n] = w, 'same_walk'
        elif by_seq.get(s):
            counterpart[n], kind[n] = max(by_seq[s])[1], 'same_seq_other_walk'
        else:
            counterpart[n], kind[n] = None, 'absent'
    lens = dict((k, len(v)) for k, v in q.seqs.items())
    lens.update((k, len(v)) for k, v in P.mc_seqs.items())
    seqs = {n: s for n, _, s in P.hap}
    memo = {}
    rows = []
    pairs = [(i, j) for i in range(len(names)) for j in range(i + 1, len(names))]
    if max_pairs and len(pairs) > max_pairs:
        import random
        pairs = sorted(random.Random(1).sample(pairs, max_pairs))
    for i, j in pairs:
        a, b = names[i], names[j]
        if counterpart[a] is None or counterpart[b] is None:
            rows.append({'pair': (a, b), 'status': 'no_counterpart'})
            continue
        c1, m1 = implied_cost(seqs[a], seqs[b], mcw[a], mcw[b], lens, C, memo)
        if counterpart[a] == mcw[a] and counterpart[b] == mcw[b]:
            rows.append({'pair': (a, b), 'status': 'identical_walks', 'cost_mc': c1, 'cost_full': c1,
                         'same_matches': True})
            continue
        c2, m2 = implied_cost(seqs[a], seqs[b], counterpart[a], counterpart[b], lens, C, memo)
        # compare as aligned base pairs (node boundaries may differ between two walks)
        def bases(m):
            return set((x + k, y + k) for x, y, L in m for k in range(L))
        rows.append({'pair': (a, b), 'status': 'different_walks', 'cost_mc': c1, 'cost_full': c2,
                     'same_matches': bases(m1) == bases(m2)})
    done = [r for r in rows if r['status'] != 'no_counterpart']
    kinds = collections.Counter(kind.values())
    res = {
        'region_id': P.rid, 'stratum': P.rj['stratum'],
        'mc_nodes': len(P.mc_seqs), 'mc_nodes_missing_in_full': node_missing,
        'mc_nodes_seq_differs': node_diff, 'mc_edges': len(P.mc_edges), 'mc_edges_missing_in_full': edge_missing,
        'hap32_paths': len(names),
        'hap32_reference': kinds.get('reference', 0), 'hap32_same_walk': kinds.get('same_walk', 0),
        'hap32_same_seq_other_walk': kinds.get('same_seq_other_walk', 0), 'hap32_absent': kinds.get('absent', 0),
        'pairs': len(rows), 'pairs_compared': len(done),
        'pairs_identical_walks': sum(1 for r in done if r['status'] == 'identical_walks'),
        'pairs_same_alignment': sum(1 for r in done if r['same_matches']),
        'pairs_same_cost': sum(1 for r in done if r['cost_mc'] == r['cost_full']),
        'cost_mc_sum': sum(r['cost_mc'] for r in done),
        'cost_full_sum': sum(r['cost_full'] for r in done),
        'not_same_walk': ';'.join('%s:%s' % (n, k) for n, k in kind.items() if k in ('same_seq_other_walk', 'absent')) or '.',
    }
    diff = [r for r in done if r['cost_mc'] != r['cost_full']]
    res['pairs_cost_differs'] = len(diff)
    res['example_differences'] = ';'.join('%s|%s:%d/%d' % (r['pair'][0], r['pair'][1], r['cost_mc'], r['cost_full'])
                                          for r in diff[:3]) or '.'
    return res


# ============================================================================ gref check

def gref_check(rids, query_dir=None, out_json=None, vg=None, gbz=None):
    """Stream `vg paths -x <eval GBZ> -R -A` (every REFERENCE-sense path: CHM13, GRCh38 and
    the gref_CHM13 cover) once and check, per region: (1) the gref_CHM13 base copy of the
    contig takes exactly the CHM13 walk between the anchors (and the whole contig's walks are
    identical); (2) no gref fragment (`gref_CHM13#0#<contig>_<N>_alt`) contains an anchor or any
    CHM13 node of the span, while counting the fragments that enter the span's other nodes;
    (3) which GRCh38 paths run anchor to anchor, and whether their walk is the one hap32.fa's
    GRCh38 record takes in mc.gfa. That is the accounting package_regions.py subtracted by weight."""
    vg = vg or config.VG
    gbz = gbz or config.data_paths()['hprc_gbz']
    info = {}
    anchor_regions = collections.defaultdict(list)     # abs node id -> [(rid, 'L'/'R')]
    span_nodes = collections.defaultdict(list)         # abs node id -> [rid] (non-CHM13 span nodes)
    ref_nodes = collections.defaultdict(list)          # abs node id -> [rid] (CHM13 span nodes)
    by_contig = collections.defaultdict(list)
    for rid in rids:
        try:
            P = Panel(rid, query_dir)
        except PanelError as e:
            log('skip %s: %s' % (rid, e))
            continue
        chm = set(abs(s) for s in P.chm13_walk)
        allspan = set()
        for w in P.sample_walks:
            allspan.update(abs(s) for s in w)
        for _, w in P.ref_walks:
            allspan.update(abs(s) for s in w)
        for n in chm:
            ref_nodes[n].append(P.rid)
        for n in allspan - chm:
            span_nodes[n].append(P.rid)
        anchor_regions[abs(P.L)].append((P.rid, 'L'))
        anchor_regions[abs(P.R)].append((P.rid, 'R'))
        info[P.rid] = {'contig': P.rj['contig'], 'L': P.L, 'R': P.R, 'chm13_walk': P.chm13_walk,
                       'grch38_hap32': {n: w for n, w in P.ref_walks if n.startswith('GRCh38#')},
                       'span_nonref_nodes': len(allspan - chm), 'span_chm13_nodes': len(chm)}
        by_contig[P.rj['contig']].append(P.rid)
        del P
    log('gref-check: %d regions, %d CHM13 span nodes, %d other span nodes' % (len(info), len(ref_nodes), len(span_nodes)))
    anchor_re = re.compile(r'[<>](%s)(?=[<>]|$)' % '|'.join(str(n) for n in sorted(anchor_regions)))
    res = {rid: {'fragments_in_span': 0, 'fragments_with_anchor_or_ref_node': 0, 'fragment_names_with_anchor': [],
                 'gref_copy_equals_chm13': None, 'chm13_path_walk_ok': None, 'grch38_spanning': 0,
                 'grch38_spanning_names': [], 'grch38_walk_matches_hap32': []} for rid in info}
    contig_hash = {}
    counts = collections.Counter()
    cmd = [vg, 'paths', '-x', gbz, '-R', '-A']
    log('$ ' + ' '.join(cmd))
    t = time.time()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=config.tool_env(),
                         bufsize=1 << 22)

    def segment(path, rid):
        """the L..R walk of a reference path (searched in either orientation)"""
        d = info[rid]
        L, R = d['L'], d['R']
        ml = re.search(r'[<>]%d(?=[<>]|$)' % abs(L), path)
        mr = re.search(r'[<>]%d(?=[<>]|$)' % abs(R), path)
        if not ml or not mr:
            return None
        a, b = sorted((ml.start(), mr.start()))
        e = re.match(r'[<>]\d+', path[b:]).end() + b
        st = parse_steps(path[a:e])
        if st and st[0] == -R:
            st = rev_steps(st)
        return st

    for raw in p.stdout:
        line = raw.decode('ascii', 'replace') if isinstance(raw, bytes) else raw
        x = line.split('\t', 6)
        if len(x) < 6:
            continue
        name, path = x[0], x[5]
        sample = sample_of(name)
        counts[sample] += 1
        if sample == 'gref_CHM13' and name.endswith('_alt'):
            counts['fragments'] += 1
            ids = set(int(n) for n in re.findall(r'\d+', path))
            hit = set()
            for n in ids:
                for rid in span_nodes.get(n, ()):
                    hit.add(rid)
            for rid in hit:
                res[rid]['fragments_in_span'] += 1
            bad = set()
            for n in ids:
                for rid in ref_nodes.get(n, ()):
                    bad.add(rid)
                for rid, _ in anchor_regions.get(n, ()):
                    bad.add(rid)
            for rid in bad:
                res[rid]['fragments_with_anchor_or_ref_node'] += 1
                res[rid]['fragment_names_with_anchor'].append(name)
            continue
        contig = name.split('#')[2] if name.count('#') >= 2 else name
        if sample in ('CHM13', 'gref_CHM13'):
            contig_hash[(sample, contig)] = (hashlib.sha1(path.encode()).hexdigest(), len(path))
            for rid in by_contig.get(contig, ()):
                st = segment(path, rid)
                key = 'gref_copy_walk' if sample == 'gref_CHM13' else 'chm13_path_walk'
                info[rid][key] = st
            continue
        if sample == 'GRCh38':
            hit = collections.defaultdict(set)
            for m in anchor_re.finditer(path):
                for rid, side in anchor_regions[int(m.group(1))]:
                    hit[rid].add(side)
            for rid, sides in hit.items():
                if sides == {'L', 'R'}:
                    st = segment(path, rid)
                    res[rid]['grch38_spanning'] += 1
                    res[rid]['grch38_spanning_names'].append(name)
                    res[rid]['grch38_walk_matches_hap32'].append(
                        any(st == w for w in info[rid]['grch38_hap32'].values()))
    err = p.stderr.read().decode('ascii', 'replace')
    p.wait()
    if p.returncode != 0:
        raise PanelError('vg paths failed (%d): %s' % (p.returncode, err[-400:]))
    log('vg paths -R -A streamed in %.0fs: %s' % (time.time() - t, dict(counts)))
    for rid, d in info.items():
        r = res[rid]
        c = d['contig']
        r['chm13_path_walk_ok'] = d.get('chm13_path_walk') == d['chm13_walk']
        r['gref_copy_equals_chm13'] = (d.get('gref_copy_walk') == d['chm13_walk']
                                       and contig_hash.get(('CHM13', c)) == contig_hash.get(('gref_CHM13', c)))
        r['gref_copy_walk_equals_chm13_in_span'] = d.get('gref_copy_walk') == d['chm13_walk']
        r['grch38_hap32'] = sorted(d['grch38_hap32'])
        r['span_nonref_nodes'] = d['span_nonref_nodes']
        r['span_chm13_nodes'] = d['span_chm13_nodes']
    whole = {c: contig_hash.get(('CHM13', c)) == contig_hash.get(('gref_CHM13', c))
             for (s, c) in contig_hash if s == 'CHM13'}
    out = {'graph': gbz, 'command': ' '.join(cmd), 'seconds': round(time.time() - t, 1),
           'paths_by_sample': dict(counts), 'whole_contig_gref_copy_equals_chm13': whole, 'regions': res}
    if out_json:
        os.makedirs(os.path.dirname(os.path.abspath(out_json)), exist_ok=True)
        with open(out_json, 'w') as f:
            json.dump(out, f, indent=1)
    return out


# ============================================================================ check-hprc

def check_hprc(rid):
    """hprc.fa.gz loads, its counts match region.json, names are clean, flanks are shared."""
    rd, rid, rj, hap, hprc = load_region(rid)
    h = rj.get('hprc', {})
    probs = []
    if h.get('status') != 'ok':
        return {'region_id': rid, 'ok': False, 'problems': 'hprc status %s' % h.get('status')}
    names = [n for n, _, _ in hprc]
    if len(set(names)) != len(names):
        probs.append('duplicate names')
    try:
        assert_no_excluded(names, rid)
    except PanelError as e:
        probs.append(str(e))
    samples = [(n, s) for n, _, s in hprc if not is_reference_name(n)]
    refs = [n for n, _, _ in hprc if is_reference_name(n)]
    hap_refs = [n for n, _, _ in hap if is_reference_name(n)]
    if sorted(refs) != sorted(hap_refs):
        probs.append('reference records %s vs hap32 %s' % (refs, hap_refs))
    if len(samples) != rj.get('n_hprc') or len(samples) != h.get('n_sample_haplotypes_spanning'):
        probs.append('samples %d vs n_hprc %s / %s' % (len(samples), rj.get('n_hprc'), h.get('n_sample_haplotypes_spanning')))
    nd = len(set(s for _, s in samples))
    if nd != h.get('n_distinct_sample_sequences'):
        probs.append('distinct %d vs %s' % (nd, h.get('n_distinct_sample_sequences')))
    s32 = set(s for _, _, s in hap)
    if sum(1 for s in set(s for _, s in samples) if s in s32) != h.get('n_distinct_in_hap32'):
        probs.append('distinct in hap32 differs')
    la, ra = rj['anchor_left_seq'], rj['anchor_right_seq']
    nf = sum(1 for _, _, s in hprc if not (s.startswith(la) and s.endswith(ra)))
    if nf:
        probs.append('%d records without the anchor flanks' % nf)
    hap_seq = {n: s for n, _, s in hap}
    for n, _, s in hprc:
        if is_reference_name(n) and hap_seq.get(n) != s:
            probs.append('reference record %s differs from hap32.fa' % n)
    grp = collections.Counter()
    size = {}
    for n, d, s in hprc:
        m = re.search(r'group=(\d+) group_size=(\d+)', d)
        if m and not is_reference_name(n):
            grp[m.group(1)] += 1
            size[m.group(1)] = int(m.group(2))
    if any(grp[g] != size[g] for g in grp):
        probs.append('group sizes do not match their records')
    return {'region_id': rid, 'ok': not probs, 'samples': len(samples), 'distinct': nd,
            'records': len(hprc), 'problems': '; '.join(probs) or '.'}


# ============================================================================ CLI

def write_tsv(path, rows, cols):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + '.tmp%d' % os.getpid()
    with open(tmp, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join('' if r.get(c) is None else str(r.get(c)) for c in cols) + '\n')
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('union', help='distinct sequences of hap32.fa + hprc.fa.gz, with a map TSV')
    s.add_argument('region')
    s.add_argument('out_fa')
    s.add_argument('--map', default=None, help='map TSV (default: OUT without .fa + .map.tsv)')
    s = sub.add_parser('union-all', help='union for every region into work/panel/union/')
    s.add_argument('--regions', default='')
    s.add_argument('--out-dir', default=UNION_DIR)
    s = sub.add_parser('project', help='restrict a full-panel MSA to hap32.fa')
    s.add_argument('full_msa')
    s.add_argument('map_tsv')
    s.add_argument('hap32_fa')
    s.add_argument('out_msa')
    s = sub.add_parser('panel-graph', help='full-panel MSA -> graph with one path per hprc.fa.gz record')
    s.add_argument('full_msa')
    s.add_argument('map_tsv')
    s.add_argument('region')
    s.add_argument('out_gfa', help='by convention work/panel/<method>/<id>.gfa')
    s.add_argument('--merge-blocks', action='store_true')
    s.add_argument('--block-max', type=int, default=50)
    s.add_argument('--engine', default='auto', choices=['auto', 'vg', 'native'])
    s = sub.add_parser('mc-graph', help='full MC graph between the anchors -> work/panel/mc/<id>.gfa')
    s.add_argument('--regions', default='')
    s.add_argument('--out-dir', default=MC_DIR)
    s.add_argument('--query-dir', default=None)
    s.add_argument('--force', action='store_true')
    s = sub.add_parser('stats', help='per-region panel statistics -> results/panel_stats.tsv')
    s.add_argument('--regions', default='')
    s.add_argument('--out', default=os.path.join(config.RESULTS_DIR, 'panel_stats.tsv'))
    s.add_argument('--query-dir', default=None)
    s.add_argument('--gref-json', default=os.path.join(PANEL_DIR, 'gref_check.json'))
    s = sub.add_parser('premise', help='mc.gfa vs the full graph restricted to the hap32 walks')
    s.add_argument('--regions', default='')
    s.add_argument('--out', default=os.path.join(config.RESULTS_DIR, 'panel_premise.tsv'))
    s.add_argument('--query-dir', default=None)
    s.add_argument('--max-pairs', type=int, default=0)
    s = sub.add_parser('gref-check', help='one vg paths -R -A pass over the eval GBZ')
    s.add_argument('--regions', default='')
    s.add_argument('--query-dir', default=None)
    s.add_argument('--out', default=os.path.join(PANEL_DIR, 'gref_check.json'))
    s = sub.add_parser('check-hprc', help='hprc.fa.gz of every region loads and matches region.json')
    s.add_argument('--regions', default='')
    s = sub.add_parser('import-queries', help='copy package_regions.py query caches into work/panel/query')
    s.add_argument('src', help='package work dir holding <id>/hprc.interval.gfa.gz')
    s.add_argument('--regions', default='')
    s = sub.add_parser('fetch', help='run the full-graph interval query for regions lacking one')
    s.add_argument('--regions', default='')
    s.add_argument('--force', action='store_true')
    a = ap.parse_args(argv)

    if a.cmd == 'union':
        print(json.dumps(union_to_files(a.region, a.out_fa, a.map), indent=1))
    elif a.cmd == 'union-all':
        tot = collections.Counter()
        for rid in pick_regions(a.regions):
            r = union_to_files(rid, os.path.join(a.out_dir, rid + '.fa'))
            tot['regions'] += 1
            tot['distinct'] += r['distinct']
            tot['hap32_only'] += r['hap32_only']
        print(json.dumps(dict(tot)))
    elif a.cmd == 'project':
        print(json.dumps(project(a.full_msa, a.map_tsv, a.hap32_fa, a.out_msa), indent=1))
    elif a.cmd == 'panel-graph':
        print(json.dumps(panel_graph(a.full_msa, a.map_tsv, a.region, a.out_gfa, a.merge_blocks, a.engine,
                                     a.block_max), indent=1))
    elif a.cmd == 'mc-graph':
        rows = []
        for rid in pick_regions(a.regions):
            out = os.path.join(a.out_dir, rid + '.gfa')
            if os.path.exists(out) and not a.force:
                continue
            t = time.time()
            try:
                P = Panel(rid, a.query_dir)
                r = P.write_gfa(out)
            except PanelError as e:
                log('FAILED %s: %s' % (rid, e))
                continue
            r['seconds'] = round(time.time() - t, 1)
            log('%s: %d paths, %d nodes, %d distinct walks, %.1f MB, %.0fs' % (
                rid, r['paths'], r['nodes'], r['distinct_walks'], r['bytes'] / 1e6, r['seconds']))
            rows.append(r)
        if rows:
            cols = ['region_id', 'paths', 'nodes', 'edges', 'node_bp', 'distinct_walks', 'distinct_seqs',
                    'mean_node_bp', 'nodes_per_kb', 'bytes', 'seconds']
            sp = os.path.join(a.out_dir, 'summary.tsv')
            old = []
            if os.path.exists(sp):
                with open(sp) as f:
                    old = [r for r in csv.DictReader(f, delimiter='\t')
                           if r['region_id'] not in set(x['region_id'] for x in rows)]
            write_tsv(sp, sorted(old + rows, key=lambda r: r['region_id']), cols)
    elif a.cmd == 'stats':
        gref = None
        if a.gref_json and os.path.exists(a.gref_json):
            gref = json.load(open(a.gref_json))['regions']
        rows = []
        for rid in pick_regions(a.regions):
            try:
                rows.append(region_stats(rid, a.query_dir, gref))
            except PanelError as e:
                log('FAILED %s: %s' % (rid, e))
        write_tsv(a.out, rows, STATS_COLS)
        log('%d regions -> %s' % (len(rows), a.out))
    elif a.cmd == 'premise':
        rows = []
        for rid in pick_regions(a.regions):
            t = time.time()
            r = premise_region(rid, a.query_dir, a.max_pairs or None)
            r['seconds'] = round(time.time() - t, 1)
            log('%s: nodes missing %d, seq differs %d, edges missing %d; hap32 same walk %d/%d; pairs same '
                'alignment %d/%d, same cost %d' % (rid, r['mc_nodes_missing_in_full'], r['mc_nodes_seq_differs'],
                                                   r['mc_edges_missing_in_full'], r['hap32_same_walk'] + r['hap32_reference'],
                                                   r['hap32_paths'], r['pairs_same_alignment'], r['pairs_compared'],
                                                   r['pairs_same_cost']))
            rows.append(r)
        cols = list(rows[0].keys()) if rows else []
        write_tsv(a.out, rows, cols)
    elif a.cmd == 'gref-check':
        out = gref_check(pick_regions(a.regions), a.query_dir, a.out)
        R = out['regions']
        summ = {'regions': len(R),
                'gref_copy_equals_chm13': sum(1 for r in R.values() if r['gref_copy_equals_chm13']),
                'chm13_path_walk_ok': sum(1 for r in R.values() if r['chm13_path_walk_ok']),
                'regions_with_fragment_on_anchor_or_ref_node': sum(1 for r in R.values() if r['fragments_with_anchor_or_ref_node']),
                'fragments_in_span_total': sum(r['fragments_in_span'] for r in R.values()),
                'regions_with_fragments_in_span': sum(1 for r in R.values() if r['fragments_in_span']),
                'grch38_matches': sum(1 for r in R.values() if r['grch38_spanning'] == len(r['grch38_hap32'])
                                      and all(r['grch38_walk_matches_hap32'])),
                'whole_contigs_equal': sum(out['whole_contig_gref_copy_equals_chm13'].values()),
                'contigs': len(out['whole_contig_gref_copy_equals_chm13'])}
        print(json.dumps(summ, indent=1))
    elif a.cmd == 'check-hprc':
        rows = [check_hprc(rid) for rid in pick_regions(a.regions)]
        bad = [r for r in rows if not r['ok']]
        for r in bad:
            print('%s\t%s' % (r['region_id'], r['problems']))
        print(json.dumps({'regions': len(rows), 'ok': len(rows) - len(bad),
                          'sample_haplotypes': sum(r.get('samples', 0) for r in rows),
                          'distinct_sample_sequences': sum(r.get('distinct', 0) for r in rows)}))
        return 1 if bad else 0
    elif a.cmd == 'import-queries':
        print(json.dumps({'imported': import_queries(a.src, pick_regions(a.regions))}))
    elif a.cmd == 'fetch':
        for rid in pick_regions(a.regions):
            if find_query(rid) and not a.force:
                continue
            fetch_query(rid, force=a.force)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except PanelError as e:
        print('panel.py: error: %s' % e, file=sys.stderr)
        sys.exit(2)
