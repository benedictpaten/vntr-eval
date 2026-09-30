#!/usr/bin/env python3
"""patch_contig.py -- Stage 4: splice realigned hotspot graphs into a whole contig's hap32 graph,
re-map the contig's reads to it and to the unpatched graph, call both, score both.

The contig-wide generalisation of call_local.py's hybrid (build_hybrid, fetch_window_reads,
map_hybrid, vg_call_command), for many regions at once:

  gfa    GFA dump of the contig's hap32 GBZ (vg convert -f) -> <out>/orig.gfa
  patch  every listed region's candidate graph (candidates/<method>/<id>.gfa) replaces the hap32
         graph between its anchors, anchors included. Every path of the GBZ is kept: where it
         crosses a region anchor to anchor its walk there becomes the candidate path of the same
         name (flipped where the GBZ stores the path reversed); W-line fields 1-5 stay as they are,
         so CHM13 and GRCh38 stay REFERENCE paths exactly as in the source GBZ (header RS tag kept;
         the NM/SG graph-name tags dropped, since the graph changes). With --fragments, a path
         that starts or ends inside a region (or on its right / left anchor) is threaded too: its
         run from its first to its last step inside the region becomes the candidate path of the
         fragment's name (<path>:enters_L|exits_R|internal|anchor_L|anchor_R, from
         realign.py --fragments-poa), which must spell the run. A region is left UNPATCHED
         (and recorded) when any path fragment starts or ends inside it (without --fragments, or
         with no matching candidate path), visits it without crossing
         anchor to anchor, crosses it under a name the candidate lacks, or when a candidate path is
         not threaded exactly once; overlapping regions keep the first by span start. Node IDs: the
         region's candidate segments (>1,024 bp chopped, as gbwtgraph would) are laid, in the
         candidate's topological order, evenly over the free IDs of [anchor_left, anchor_right]
         (the region's own IDs, freed), so vg call's 4,096-ID depth-rate windows mix a region and
         its flanks as in the source graph; where the range has too few free IDs the region's IDs
         follow the contig's largest ID instead (ids_interleaved false, recorded). `--ids none`
         writes the unpatched graph through the same code (the unpatched arm).
         Then `vg gbwt -G --gbz-format` and ASSERT (verify): same path names and metadata, every
         path spells the same sequence as in the source GBZ (md5 per path over `vg paths -F`), the
         GBZ's node count/length are the GFA's and its CHM13 path has the GFA's node IDs (no
         segment translation).
  reads  the reads the production mapping placed on the contig, rebuilt from GAF-Base exactly as
         remap_local.fetch rebuilds them (path sequence + cs, base qualities, clipped alignments
         dropped, mates kept together): gbz-base --alignments overlapping over CHM13 chunks of the
         contig (--context 10 kb), then every contig node no chunk covered queried by ID, one copy
         per (read name, mate). -> reads_1.fq.gz, reads_2.fq.gz, reads_se.fq.gz, reads.json
  map    vg autoindex -n -w sr-giraffe on the GFA, giraffe paired (Stage 2's fragment-length
         distribution) + single-end, GAF --named-coordinates, bgzipped; index files deleted.
  call   the pinned vg with the production flags (call_local.vg_call_command), --gaf-reads.
  score  vg-call-eval's short-read harness (scripts/wgs/bench_wgs.score_contig: aardvark for small
         variants, truvari with the pipeline's parameters for SVs), then truvari FP/FN inside vs
         outside the patched regions, and score_haplotypes.py batch over the patched regions.

    python3 tools/patch_contig.py gfa   --out work/stage4/chr20
    python3 tools/patch_contig.py patch --out work/stage4/chr20/patched --ids FILE --method poa_abpoa [--fragments]
    python3 tools/patch_contig.py patch --out work/stage4/chr20/unpatched --ids none
    python3 tools/patch_contig.py reads --out work/stage4/chr20/reads --jobs 4
    python3 tools/patch_contig.py map   --gfa .../patched/graph.gfa --reads .../reads --out .../patched --threads 8
    python3 tools/patch_contig.py call  --gbz .../patched/graph.gbz --gaf .../patched/reads.gaf.gz --out .../patched/chr20.vcf.gz
    python3 tools/patch_contig.py score --vcf V --label L --out work/stage4/chr20/score --patch .../patched/patch.json
Regions and candidates come from $VNTR_REGIONS / $VNTR_CANDIDATES (config.py).
"""
import argparse
import collections
import gzip
import hashlib
import json
import multiprocessing
import os
import re
import subprocess
import sys
import time

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402
import call_local as cl  # noqa: E402
import remap_local as rl  # noqa: E402

VG = cl.VG_PINNED
STEP_RE = re.compile(r'([<>])([^<>]+)')
MAX_NODE = cl.MAX_NODE
FRAGLEN = os.path.join(cl.STAGE2_WORK, 'fraglen.json')


def log(*a):
    print('[%s] %s' % (time.strftime('%H:%M:%S'), ' '.join(str(x) for x in a)), file=sys.stderr, flush=True)


def run(cmd, **kw):
    return cl.run([str(c) for c in cmd], **kw)


def rc(s):
    return cl.rc(s)


def flip(steps):
    return [('<' if o == '>' else '>', n) for o, n in reversed(steps)]


def contig_gbz(contig):
    return config.data_paths(contig)['contig_gbz']


# ------------------------------------------------------------------ gfa

def cmd_gfa(a):
    os.makedirs(a.out, exist_ok=True)
    out = os.path.join(a.out, 'orig.gfa')
    if not os.path.exists(out):
        run([VG, 'convert', '-f', contig_gbz(a.contig)], stdout=out + '.tmp')
        os.replace(out + '.tmp', out)
    return out


# ------------------------------------------------------------------ patch

def load_regions(ids, method):
    regs = []
    for rid in ids:
        r = rl.load_region(rid)
        r['candidate'] = os.path.join(config.CANDIDATES_DIR, method, rid + '.gfa')
        regs.append(r)
    return regs


def between(reg, od):
    """{node: seq} of the hap32 graph between the anchors, anchors included (gbz-base --between)."""
    p = os.path.join(od, 'between', reg['region_id'] + '.gfa')
    if not os.path.exists(p):
        os.makedirs(os.path.dirname(p), exist_ok=True)
        run([config.GBZ_BASE, 'query', '--between', '%s:%s' % (reg['anchor_left'], reg['anchor_right']),
             '--limit', '2000000', config.data_paths(reg['contig'])['gbz_db']], stdout=p + '.tmp')
        os.replace(p + '.tmp', p)
    return rl.parse_gfa(p)[0]


def classify_fragment(seg, AL, AR, at_start, at_end):
    """A run of a path inside a region that does not cross anchor to anchor, as a threadable
    fragment: (class, reversed) with class enters_L / anchor_L (a prefix of the region, the path
    ends inside or on the left anchor), exits_R / anchor_R (a suffix, the path starts inside or on
    the right anchor) or internal (the whole path is inside; reversed is None, fixed later by
    sequence), named as in hap32.fragments.fa / realign.region_fragments. None when the run is
    not at an end of its path (a visit that leaves the way it came)."""
    L, R = ('>', AL), ('>', AR)
    Lr, Rr = ('<', AL), ('<', AR)
    if seg[0] == L and at_end:
        return ('anchor_L' if len(seg) == 1 else 'enters_L'), False
    if seg[-1] == Lr and at_start:
        return ('anchor_L' if len(seg) == 1 else 'enters_L'), True
    if seg[-1] == R and at_start:
        return ('anchor_R' if len(seg) == 1 else 'exits_R'), False
    if seg[0] == Rr and at_end:
        return ('anchor_R' if len(seg) == 1 else 'exits_R'), True
    if at_start and at_end and not {L, Lr, R, Rr} & set(seg):
        return 'internal', None
    return None


def cmd_patch(a):
    os.makedirs(a.out, exist_ok=True)
    src = a.gfa or os.path.join(os.path.dirname(a.out.rstrip('/')), 'orig.gfa')
    ids = [] if a.ids == 'none' else [x.split()[0] for x in open(a.ids) if x.strip() and not x.startswith('#')]
    regs = sorted(load_regions(ids, a.method), key=lambda r: (r['span_start'], r['span_end']))
    status = collections.OrderedDict()      # rid -> reason ('patched' at the end)
    inside, node2reg, info = {}, {}, {}
    last_end = -1
    for r in regs:
        rid = r['region_id']
        if r['anchor_left'][-1] != '+' or r['anchor_right'][-1] != '+':
            status[rid] = 'anchor_on_reverse_strand'
            continue
        if not os.path.exists(r['candidate']):
            status[rid] = 'no_candidate'
            continue
        nodes = between(r, a.out)
        if r['span_start'] <= last_end or any(n in node2reg for n in nodes):
            status[rid] = 'overlaps_earlier_region'
            continue
        last_end = r['span_end']
        inside[rid] = nodes
        for n in nodes:
            node2reg[n] = rid
        info[rid] = r
    # pass 1: node IDs, and every path's runs of steps inside a region
    log('patch: pass 1 over', src)
    all_ids = set()
    runs = collections.defaultdict(list)     # W-line index -> [(i, j, rid, rev)]
    spans = collections.defaultdict(collections.Counter)   # rid -> name -> spanning runs
    bad = collections.defaultdict(collections.Counter)
    fcount = collections.defaultdict(collections.Counter)
    frag_pending = collections.defaultdict(list)   # rid -> [(W index, i, j, name, seq as stored, rev|None)]
    run_seq = {}
    wi = -1
    with open(src) as f:
        for line in f:
            if line[0] == 'S':
                all_ids.add(int(line.split('\t', 2)[1]))
            elif line[0] == 'W':
                wi += 1
                x = line.rstrip('\n').split('\t')
                name = '%s#%s#%s' % (x[1], x[2], x[3])
                if x[1] == 'recombination':
                    name += '#%s' % x[4]
                elif x[1] == 'GRCh38':
                    name += '[%s]' % x[4]
                steps = STEP_RE.findall(x[6])
                idx = [i for i, (_, n) in enumerate(steps) if n in node2reg]
                k = 0
                while k < len(idx):
                    rid = node2reg[steps[idx[k]][1]]
                    m = k
                    while m + 1 < len(idx) and idx[m + 1] == idx[m] + 1 and node2reg[steps[idx[m + 1]][1]] == rid:
                        m += 1
                    i, j = idx[k], idx[m] + 1
                    AL, AR = info[rid]['anchor_left'][:-1], info[rid]['anchor_right'][:-1]
                    seg = steps[i:j]
                    if seg[0] == ('>', AL) and seg[-1] == ('>', AR):
                        rev = False
                    elif seg[0] == ('<', AR) and seg[-1] == ('<', AL):
                        rev = True
                    else:
                        fr = classify_fragment(seg, AL, AR, i == 0, j == len(steps)) if a.fragments else None
                        if fr is None:
                            why = ('fragment_starts_inside' if i == 0 else 'fragment_ends_inside' if j == len(steps)
                                   else 'visit_not_anchor_to_anchor')
                            bad[rid][why] += 1
                        else:
                            cls, frev = fr
                            nd = inside[rid]
                            sq = ''.join(nd[n] if o == '>' else rc(nd[n]) for o, n in seg)
                            base = '%s:%s' % (name, cls)
                            fcount[rid][base] += 1
                            fname = base if fcount[rid][base] == 1 else '%s:%d' % (base, fcount[rid][base])
                            frag_pending[rid].append((wi, i, j, fname, sq, frev))
                        k = m + 1
                        continue
                    fw = flip(seg) if rev else seg
                    nd = inside[rid]
                    run_seq[(rid, name, spans[rid][name])] = ''.join(nd[n] if o == '>' else rc(nd[n]) for o, n in fw)
                    spans[rid][name] += 1
                    runs[wi].append((i, j, rid, rev, name))
                    k = m + 1
    # candidates: every path spells hap32.fa and the path's own run; every hap32 name threaded once
    cand, nfrag = {}, {}
    for rid in list(inside):
        if bad[rid]:
            status[rid] = 'unpatched:' + ','.join('%s=%d' % kv for kv in sorted(bad[rid].items()))
            continue
        snodes, sedges, spaths = rl.parse_gfa(info[rid]['candidate'])
        hap = dict(cl.read_fasta(os.path.join(info[rid]['dir'], 'hap32.fa')))
        sp = {nm: [('>' if o == '+' else '<', n) for n, o in st] for nm, st in spaths}
        why = None
        for nm, st in sp.items():
            if nm in hap and ''.join(snodes[n] if o == '>' else rc(snodes[n]) for o, n in st) != hap[nm]:
                raise SystemExit('%s: candidate path %s does not spell hap32.fa' % (rid, nm))
        for nm, c in spans[rid].items():
            if nm not in sp:
                why = 'crossing_path_not_in_candidate:%s' % nm
            elif c != 1:
                why = 'path_crosses_%d_times:%s' % (c, nm)
            elif run_seq[(rid, nm, 0)] != hap.get(nm):
                raise SystemExit('%s: %s crosses with a sequence other than hap32.fa' % (rid, nm))
        for nm in hap:
            if spans[rid][nm] != 1:
                why = why or 'hap32_path_not_crossing:%s' % nm
        # fragments: each is threaded along the candidate path of its name, in the orientation
        # that spells the run (fixed by the anchor it touches; either for an internal piece)
        fruns = []
        for wi_, i, j, fname, sq, frev in frag_pending[rid]:
            if fname not in sp:
                why = why or 'fragment_not_in_candidate:%s' % fname
                continue
            cs = ''.join(snodes[n] if o == '>' else rc(snodes[n]) for o, n in sp[fname])
            ok = [r for r in ((False, True) if frev is None else (frev,)) if cs == (rc(sq) if r else sq)]
            if not ok:
                why = why or 'fragment_sequence_differs:%s' % fname
                continue
            fruns.append((wi_, i, j, ok[0], fname))
        if why:
            status[rid] = 'unpatched:' + why
            continue
        for wi_, i, j, rev, fname in fruns:
            runs[wi_].append((i, j, rid, rev, fname))
        nfrag[rid] = len(fruns)
        cand[rid] = (snodes, sedges, sp)
    removed = set(n for rid in cand for n in inside[rid])
    taken = {i for i in all_ids if str(i) not in removed}
    top = max(all_ids) + 1
    chop, new_nodes, new_edges, reginfo = {}, [], set(), {}
    for rid in sorted(cand, key=lambda r: info[r]['span_start']):
        snodes, sedges, sp = cand[rid]
        AL, AR = int(info[rid]['anchor_left'][:-1]), int(info[rid]['anchor_right'][:-1])
        lo, hi = min(AL, AR), max(AL, AR)
        order = sorted(snodes, key=int)
        pieces = {n: [snodes[n][i:i + MAX_NODE] for i in range(0, len(snodes[n]), MAX_NODE)] for n in order}
        n_new = sum(len(v) for v in pieces.values())
        free = [i for i in range(lo, hi + 1) if i not in taken]
        if len(free) >= n_new:
            step = len(free) / float(n_new)
            slots = [free[min(len(free) - 1, int(k * step))] for k in range(n_new)]
            inter = True
        else:
            slots = list(range(top, top + n_new))
            top += n_new
            inter = False
        assert len(set(slots)) == n_new and not (set(slots) & taken)
        taken.update(slots)
        ch = {}
        k = 0
        for n in order:
            ch[n] = []
            for s in pieces[n]:
                new_nodes.append((slots[k], s))
                ch[n].append(str(slots[k]))
                k += 1
        for u in ch.values():
            for p, q in zip(u, u[1:]):
                new_edges.add((p, '+', q, '+'))
        for x, xo, y, yo in sedges:
            new_edges.add((ch[x][-1] if xo == '+' else ch[x][0], xo, ch[y][0] if yo == '+' else ch[y][-1], yo))
        chop[rid] = {nm: [(o, z) for o, n in st for z in (ch[n] if o == '>' else ch[n][::-1])] for nm, st in sp.items()}
        status[rid] = 'patched'
        reginfo[rid] = {'old_nodes': len(inside[rid]), 'new_nodes': n_new, 'candidate_segments': len(snodes),
                        'id_range': [lo, hi], 'free_ids': len(free), 'ids_interleaved': inter,
                        'paths_crossing': sum(spans[rid].values()), 'fragments_threaded': nfrag.get(rid, 0),
                        'span': [info[rid]['span_start'], info[rid]['span_end']]}
    # pass 2: write
    out_gfa = os.path.join(a.out, 'graph.gfa')
    log('patch: writing', out_gfa, '(%d regions patched)' % len(cand))
    wi = -1
    n_s = n_len = 0
    with open(src) as f, open(out_gfa + '.tmp', 'w') as fo:
        for line in f:
            t = line[0]
            if t == 'H':
                if 'RS:Z:' in line:
                    fo.write(line)
            elif t == 'S':
                n = line.split('\t', 2)[1]
                if n not in removed:
                    fo.write(line)
                    n_s += 1
                    n_len += len(line.rstrip('\n').split('\t')[2])
            elif t == 'L':
                x = line.split('\t', 4)
                if x[1] not in removed and x[3] not in removed:
                    fo.write(line)
            elif t == 'W':
                wi += 1
                rr = sorted(r for r in runs.get(wi, []) if r[2] in cand)
                if not rr:
                    fo.write(line)
                    continue
                x = line.rstrip('\n').split('\t')
                steps = STEP_RE.findall(x[6])
                out, prev = [], 0
                for i, j, rid, rev, name in rr:
                    out += steps[prev:i]
                    rep = chop[rid][name]
                    rep = flip(rep) if rev else rep
                    for p, q in ((steps[i - 1] if i > 0 else None, rep[0]), (rep[-1], steps[j] if j < len(steps) else None)):
                        if p and q:
                            new_edges.add((p[1], '+' if p[0] == '>' else '-', q[1], '+' if q[0] == '>' else '-'))
                    for p, q in zip(rep, rep[1:]):
                        new_edges.add((p[1], '+' if p[0] == '>' else '-', q[1], '+' if q[0] == '>' else '-'))
                    out += rep
                    prev = j
                out += steps[prev:]
                x[6] = ''.join(o + n for o, n in out)
                fo.write('\t'.join(x) + '\n')
        for nid, s in new_nodes:
            fo.write('S\t%d\t%s\n' % (nid, s))
            n_s += 1
            n_len += len(s)
        for e in sorted(new_edges):
            fo.write('L\t%s\t%s\t%s\t%s\t0M\n' % e)
    os.replace(out_gfa + '.tmp', out_gfa)
    gbz = os.path.join(a.out, 'graph.gbz')
    log('patch: vg gbwt -G')
    run([VG, 'gbwt', '-G', out_gfa, '--gbz-format', '-g', gbz], stderr=os.path.join(a.out, 'gbwt.log'))
    res = {'source_gfa': src, 'method': a.method, 'regions_listed': len(ids), 'patched': len(cand),
           'status': status, 'regions': reginfo, 'nodes': n_s, 'length': n_len,
           'removed_nodes': len(removed), 'added_nodes': len(new_nodes)}
    res['verify'] = verify(contig_gbz(a.contig), gbz, out_gfa, n_s, n_len)
    with open(os.path.join(a.out, 'patch.json'), 'w') as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(a.out, 'patch.tsv'), 'w') as f:
        f.write('region_id\tstatus\tspan_start\tspan_end\told_nodes\tnew_nodes\tids_interleaved\n')
        for rid, st in status.items():
            g = reginfo.get(rid, {})
            r = next(x for x in regs if x['region_id'] == rid)
            f.write('%s\t%s\t%d\t%d\t%s\t%s\t%s\n' % (rid, st, r['span_start'], r['span_end'], g.get('old_nodes', ''),
                                                    g.get('new_nodes', ''), g.get('ids_interleaved', '')))
    log('patch: done', json.dumps({k: res[k] for k in ('patched', 'nodes', 'removed_nodes', 'added_nodes')}))


def path_digests(gbz):
    """{path name: (md5 of its sequence, length)} over `vg paths -F`."""
    p = subprocess.Popen([VG, 'paths', '-x', gbz, '-F'], stdout=subprocess.PIPE, text=True, env=config.tool_env())
    out, name, h, n = {}, None, None, 0
    for line in p.stdout:
        if line.startswith('>'):
            if name is not None:
                out[name] = (h.hexdigest(), n)
            name, h, n = line[1:].split()[0], hashlib.md5(), 0
        else:
            s = line.strip().upper()
            h.update(s.encode())
            n += len(s)
    if name is not None:
        out[name] = (h.hexdigest(), n)
    if p.wait():
        raise SystemExit('vg paths -F failed on %s' % gbz)
    return out


def verify(orig, new, gfa, n_s, n_len):
    log('verify: path metadata')
    meta = [run([VG, 'paths', '-x', g, '-M']).stdout.splitlines() for g in (orig, new)]
    if sorted(meta[0]) != sorted(meta[1]):
        raise SystemExit('ASSERT: path metadata differ between %s and %s' % (orig, new))
    log('verify: path sequences')
    d0, d1 = path_digests(orig), path_digests(new)
    if d0 != d1:
        diff = [k for k in set(d0) | set(d1) if d0.get(k) != d1.get(k)]
        raise SystemExit('ASSERT: %d paths spell different sequences, e.g. %s' % (len(diff), diff[:3]))
    st = dict(l.split('\t')[:2] for l in run([VG, 'stats', '-z', '-l', new]).stdout.splitlines() if '\t' in l)
    if int(st['nodes']) != n_s or int(st['length']) != n_len:
        raise SystemExit('ASSERT: GBZ has %s nodes/%s bp, the GFA %d/%d' % (st['nodes'], st['length'], n_s, n_len))
    chm = run([VG, 'paths', '-x', new, '-A', '-Q', 'CHM13']).stdout.split('\n')
    chm = [l for l in chm if l and not l.startswith('@')]
    want = None
    with open(gfa) as f:
        for line in f:
            if line.startswith('W\tCHM13\t'):
                want = line.rstrip('\n').split('\t')[6]
    if len(chm) != 1 or chm[0].split('\t')[5] != want:
        raise SystemExit('ASSERT: the GBZ CHM13 path does not use the GFA node IDs (translation?)')
    log('verify: OK (%d paths identical)' % len(d0))
    return {'paths': len(d0), 'metadata_identical': True, 'sequences_identical': True,
            'gbz_nodes': int(st['nodes']), 'gbz_length': int(st['length']), 'chm13_node_ids_as_gfa': True}


# ------------------------------------------------------------------ reads

NODES = {}


def rebuild_gaf(gaf, tsv):
    """remap_local.fetch's rebuild, per alignment -> TSV name, mate, seq, qual. Returns counts."""
    st = collections.Counter()
    with open(gaf) as f, open(tsv, 'w') as fo:
        for line in f:
            if line.startswith('@'):
                continue
            r = rl.parse_gaf_line(line)
            if not r:
                continue
            st['alignments'] += 1
            if not r['steps']:
                st['no_path'] += 1
                continue
            if any(n not in NODES for n, _ in r['steps']):
                st['missing_node'] += 1
                continue
            pseq = rl.walk_seq(NODES, r['steps'])
            if len(pseq) != r['plen']:
                st['path_length_mismatch'] += 1
                continue
            s = rl.read_from_cs(pseq[r['ps']:r['pe']], r['cs'])
            if s is None:
                st['cs_mismatch'] += 1
                continue
            if r['qs'] != 0 or r['qe'] != r['qlen'] or len(s) != r['qlen']:
                st['clipped_or_short'] += 1
                continue
            q = r['bq'] if r['bq'] and len(r['bq']) == len(s) else None
            if q is None:
                st['no_quality'] += 1
                q = 'I' * len(s)
            fo.write('%s\t%d\t%s\t%s\n' % (r['name'], r['mate'], s, q))
            st['rebuilt'] += 1
    return st


def _chunk(job):
    kind, arg, wd, contig = job
    base = os.path.join(wd, 'chunks', '%s_%s' % (kind, arg[0]))
    if os.path.exists(base + '.json'):
        with open(base + '.json') as f:
            d = json.load(f)
        return collections.Counter(d['counts']), set(d['covered'])
    dp = config.data_paths(contig)
    q = [config.GBZ_BASE, 'query']
    if kind == 'iv':
        q += ['--sample', 'CHM13', '--contig', contig, '--interval', '%d..%d' % arg, '--context', '10000']
    else:
        for n in arg[1]:
            q += ['-n', n]
        q += ['--context', '0']
    q += ['--gaf-base', dp['gaf_db'], '--gaf-output', base + '.gaf', '--alignments', 'overlapping', dp['gbz_db']]
    run(q, stdout=base + '.gfa')
    covered = []
    with open(base + '.gfa') as f:
        for line in f:
            if line[0] == 'S':
                covered.append(line.split('\t', 2)[1])
    st = rebuild_gaf(base + '.gaf', base + '.tsv')
    os.remove(base + '.gaf')
    os.remove(base + '.gfa')
    with open(base + '.json', 'w') as f:
        json.dump({'counts': st, 'covered': covered if kind == 'iv' else []}, f)
    return st, set(covered)


def cmd_reads(a):
    wd = a.out
    os.makedirs(os.path.join(wd, 'chunks'), exist_ok=True)
    gfa = a.gfa or os.path.join(os.path.dirname(wd.rstrip('/')), 'orig.gfa')
    log('reads: loading node sequences from', gfa)
    with open(gfa) as f:
        for line in f:
            if line[0] == 'S':
                x = line.rstrip('\n').split('\t')
                NODES[x[1]] = x[2].upper()
    clen = cl.contig_length(a.contig)
    ivs = [(s, min(clen, s + a.chunk)) for s in range(0, clen, a.chunk)]
    ctx = multiprocessing.get_context('fork')
    total = collections.Counter()
    covered = set()
    with ctx.Pool(a.jobs) as pool:
        for st, cov in pool.imap_unordered(_chunk, [('iv', iv, wd, a.contig) for iv in ivs]):
            total.update(st)
            covered |= cov
        miss = sorted(n for n in NODES if n not in covered)
        log('reads: %d interval chunks, %d alignments; %d contig nodes outside every chunk' % (len(ivs), total['alignments'], len(miss)))
        batches = [(k, miss[k:k + 3000]) for k in range(0, len(miss), 3000)]
        for st, _ in pool.imap_unordered(_chunk, [('nodes', b, wd, a.contig) for b in batches]):
            total.update(st)
    log('reads: sorting')
    tsvs = sorted(os.path.join(wd, 'chunks', x) for x in os.listdir(os.path.join(wd, 'chunks')) if x.endswith('.tsv'))
    srt = os.path.join(wd, 'sorted.tsv')
    tmp = os.path.join(wd, 'sorttmp')
    os.makedirs(tmp, exist_ok=True)
    env = dict(os.environ, LC_ALL='C')
    subprocess.run(['sort', '-t', '\t', '-k1,1', '-k2,2n', '-S', '3G', '--parallel=%d' % a.jobs, '-T', tmp, '-o', srt] + tsvs,
                   check=True, env=env)
    log('reads: pairing')
    outs = {k: subprocess.Popen('bgzip -@ 2 -c > %s' % os.path.join(wd, 'reads_%s.fq.gz' % k), shell=True,
                                stdin=subprocess.PIPE, text=True) for k in ('1', '2', 'se')}
    cnt = collections.Counter()

    def flush(nm, d):
        if '1' in d and '2' in d:
            ln = 'p%d' % cnt['pairs']
            cnt['pairs'] += 1
            for m in ('1', '2'):
                outs[m].stdin.write('@%s\n%s\n+\n%s\n' % (ln, d[m][0], d[m][1]))
        else:
            for m in sorted(d):
                outs['se'].stdin.write('@s%d\n%s\n+\n%s\n' % (cnt['singles'], d[m][0], d[m][1]))
                cnt['singles'] += 1

    cur, d = None, {}
    with open(srt) as f:
        for line in f:
            nm, m, s, q = line.rstrip('\n').split('\t')
            if nm != cur:
                if cur is not None:
                    flush(cur, d)
                cur, d = nm, {}
            if m in d:
                cnt['duplicate_copies'] += 1
                if d[m] != (s, q):
                    cnt['duplicates_differing'] += 1
                continue
            d[m] = (s, q)
        if cur is not None:
            flush(cur, d)
    for p in outs.values():
        p.stdin.close()
        p.wait()
    os.remove(srt)
    for t in tsvs:
        os.remove(t)
    info = {'contig': a.contig, 'chunk_bp': a.chunk, 'context_bp': 10000, 'counts': dict(total), **cnt,
            'reads_rebuilt_unique': 2 * cnt['pairs'] + cnt['singles']}
    with open(os.path.join(wd, 'reads.json'), 'w') as f:
        json.dump(info, f, indent=1)
    log('reads: done', json.dumps(info))


# ------------------------------------------------------------------ map

def cmd_map(a):
    od = a.out
    out = os.path.join(od, 'reads.gaf.gz')
    if os.path.exists(out):
        return out
    with open(FRAGLEN) as f:
        frag = json.load(f)
    prefix = os.path.join(od, 'idx')
    tmp = os.path.join(od, 'tmp')
    os.makedirs(tmp, exist_ok=True)
    t = str(a.threads)
    log('map: autoindex', a.gfa)
    run([VG, 'autoindex', '-n', '-w', 'sr-giraffe', '-g', a.gfa, '-p', prefix, '-t', t, '-T', tmp, '-M', a.mem],
        stdout=os.path.join(od, 'autoindex.out'), stderr=os.path.join(od, 'autoindex.log'))
    mn = prefix + '.shortread.withzip.min'
    if not os.path.exists(mn):
        mn = prefix + '.min'
    base = [VG, 'giraffe', '-Z', prefix + '.giraffe.gbz', '-d', prefix + '.dist', '-m', mn]
    if os.path.exists(prefix + '.shortread.zipcodes'):
        base += ['-z', prefix + '.shortread.zipcodes']
    base += ['-o', 'gaf', '--named-coordinates', '-t', t]
    r1, r2, se = (os.path.join(a.reads, 'reads_%s.fq.gz' % k) for k in ('1', '2', 'se'))
    parts = []
    for tag, args in (('paired', ['-f', r1, '-f', r2, '--fragment-mean', '%.1f' % frag['mean'],
                                  '--fragment-stdev', '%.1f' % frag['stdev']]), ('single', ['-f', se])):
        part = os.path.join(od, 'reads.%s.gaf.gz' % tag)
        log('map: giraffe', tag)
        with open(part + '.tmp', 'wb') as fo, open(os.path.join(od, 'giraffe.%s.log' % tag), 'w') as fe:
            g = subprocess.Popen([str(c) for c in base + args], stdout=subprocess.PIPE, stderr=fe, env=config.tool_env())
            z = subprocess.Popen(['bgzip', '-@', '2', '-c'], stdin=g.stdout, stdout=fo)
            g.stdout.close()
            if g.wait() or z.wait():
                raise SystemExit('giraffe %s failed' % tag)
        os.replace(part + '.tmp', part)
        parts.append(part)
    # every GAF node is a GFA segment (segment names are node IDs: nothing is chopped by the indexer)
    seg = set()
    with open(a.gfa) as f:
        for line in f:
            if line[0] == 'S':
                seg.add(line.split('\t', 2)[1])
    n = unmapped = 0
    for part in parts:
        with gzip.open(part, 'rt') as f:
            for line in f:
                if line[0] == '@':
                    continue
                x = line.split('\t', 7)
                n += 1
                if x[5] == '*':
                    unmapped += 1
                elif any(v not in seg for _, v in STEP_RE.findall(x[5])):
                    raise SystemExit('GAF names a node the GFA does not have: %s' % x[5][:80])
    with open(out + '.tmp', 'wb') as fo:
        for part in parts:
            with open(part, 'rb') as f:
                while True:
                    b = f.read(1 << 24)
                    if not b:
                        break
                    fo.write(b)
    os.replace(out + '.tmp', out)
    for part in parts:
        os.remove(part)
    for x in os.listdir(od):
        if x.startswith('idx.') and not x.endswith('.log'):
            os.remove(os.path.join(od, x))
    subprocess.run(['rm', '-rf', tmp])
    with open(os.path.join(od, 'map.json'), 'w') as f:
        json.dump({'records': n, 'unmapped': unmapped, 'fraglen': [frag['mean'], frag['stdev']]}, f)
    log('map: done', n, 'records', unmapped, 'unmapped')
    return out


# ------------------------------------------------------------------ call

def cmd_call(a):
    out = a.out
    raw = out[:-3] if out.endswith('.gz') else out
    cmd = cl.vg_call_command({'contig': a.contig}, a.gbz, a.gaf, raw[:-4] + '.mosaic.tsv', a.threads)
    with open(raw[:-4] + '.cmd.json', 'w') as f:
        json.dump({'cmd': cmd, 'production_source': cl.PRODUCTION_SOURCE}, f, indent=1)
    log('call:', ' '.join(cmd))
    run(['/usr/bin/time', '-l'] + cmd, stdout=raw, stderr=raw[:-4] + '.log')
    run(['bgzip', '-f', raw])
    run(['tabix', '-f', '-p', 'vcf', raw + '.gz'])
    log('call: done', raw + '.gz')


# ------------------------------------------------------------------ score

def overlaps(ivs, pos, end):
    return any(s <= end and pos <= e for s, e in ivs)


def vcf_intervals(path):
    out = []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            x = line.split('\t', 5)
            out.append((int(x[1]), int(x[1]) + len(x[3]) - 1))
    return out


def cmd_score(a):
    from pathlib import Path
    ev = config.EVAL_DIR
    wd = os.path.join(a.out, a.label)
    d = os.path.join(wd, a.contig)
    os.makedirs(d, exist_ok=True)
    prod = os.path.join(ev, 'work/wgs-mm095', a.contig)
    for x in ('%s.fa', '%s.fa.fai', 'truth.%s.smvar.bed', 'truth.%s.smvar.vcf.gz', 'truth.%s.smvar.vcf.gz.tbi',
              'truth.%s.stvar.bed', 'truth.%s.stvar.vcf.gz', 'truth.%s.stvar.vcf.gz.tbi'):
        p = os.path.join(d, x % a.contig)
        if not os.path.lexists(p):
            os.symlink(os.path.join(prod, x % a.contig), p)
    for ext in ('', '.tbi'):
        p = os.path.join(d, '%s.vcf.gz%s' % (a.contig, ext))
        if os.path.lexists(p):
            os.remove(p)
        os.symlink(os.path.abspath(a.vcf + ext), p)
    sys.path.insert(0, os.path.join(ev, 'scripts/wgs'))
    sys.path.insert(0, os.path.join(ev, 'scripts'))
    import bench_wgs
    import bench_metrics as bm
    res = os.path.join(wd, 'score.json')
    if not os.path.exists(res) or os.path.getmtime(res) < os.path.getmtime(a.vcf):
        log('score: bench_wgs.score_contig', a.label)
        r = bench_wgs.score_contig(Path(wd), a.contig, 'HG002', a.threads, cl.TRUVARI)
        with open(res, 'w') as f:
            json.dump(r, f, indent=1)
    with open(res) as f:
        r = json.load(f)
    rows = {'label': a.label}
    sv = bm.sv([r], [a.contig])
    rows.update(sv_precision=sv.precision, sv_recall=sv.recall, sv_f1=sv.f1, sv_fp=sv.query_fp, sv_fn=sv.truth_fn,
                sv_tp_base=sv.truth_tp, sv_tp_comp=sv.query_tp)
    # JointIndel, as bench_wgs reports indels: aardvark's plain Indel row is query-only
    for vt, lab in (('ALL', 'all'), ('Snv', 'snv'), ('JointIndel', 'indel')):
        c = bm.small([r], [a.contig], vt)
        rows.update({'%s_f1' % lab: c.f1, '%s_precision' % lab: c.precision,
                     '%s_recall' % lab: c.recall, '%s_fp' % lab: c.query_fp, '%s_fn' % lab: c.truth_fn})
    # SV FP/FN inside vs outside the patched regions (anchor-to-anchor span +- pad)
    with open(a.patch) as f:
        pj = json.load(f)
    ivs = [(g['span'][0] - a.pad, g['span'][1] + a.pad) for g in pj['regions'].values()]
    tdir = os.path.join(wd, 'score', '%s.truvari' % a.contig)
    for kind, fn in (('fp', 'fp.vcf.gz'), ('fn', 'fn.vcf.gz')):
        iv = vcf_intervals(os.path.join(tdir, fn))
        k = sum(1 for p, e in iv if overlaps(ivs, p, e))
        rows['sv_%s_in_patched' % kind] = k
        rows['sv_%s_outside' % kind] = len(iv) - k
    # haplotype ED over the patched regions
    ids = sorted(pj['regions'])
    idf = os.path.join(wd, 'patched_ids.txt')
    with open(idf, 'w') as f:
        f.write('\n'.join(ids) + '\n')
    hd = os.path.join(wd, 'haps')
    if not a.no_haps:
        run([sys.executable, os.path.join(TOOLS, 'score_haplotypes.py'), 'batch', '--label', a.label, '--vcf',
             os.path.abspath(a.vcf), '--regions', idf, '--jobs', str(a.threads), '--no-truvari', '--out-dir', hd]
            + (['--force'] if a.force else []), stdout=os.path.join(wd, 'haps.log'), stderr=subprocess.STDOUT)
    with open(os.path.join(wd, 'summary.json'), 'w') as f:
        json.dump(rows, f, indent=1)
    log('score:', json.dumps(rows))


# ------------------------------------------------------------------ main

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='cmd', required=True)
    for name in ('gfa', 'patch', 'reads', 'map', 'call', 'score'):
        s = sp.add_parser(name)
        s.add_argument('--contig', default='chr20')
        s.add_argument('--out', required=True)
        if name in ('patch', 'reads'):
            s.add_argument('--gfa', help='the contig GFA dump (default: <out>/../orig.gfa)')
        if name == 'patch':
            s.add_argument('--ids', required=True, help="file of region IDs, or 'none' (unpatched)")
            s.add_argument('--method', default='poa_abpoa')
            s.add_argument('--fragments', action='store_true',
                           help='thread paths that start or end inside a region along the candidate\'s '
                                'fragment paths (realign.py --fragments-poa) instead of leaving the region unpatched')
        if name == 'reads':
            s.add_argument('--jobs', type=int, default=4)
            s.add_argument('--chunk', type=int, default=2000000)
        if name == 'map':
            s.add_argument('--gfa', required=True)
            s.add_argument('--reads', required=True)
            s.add_argument('--threads', type=int, default=8)
            s.add_argument('--mem', default='8G')
        if name == 'call':
            s.add_argument('--gbz', required=True)
            s.add_argument('--gaf', required=True)
            s.add_argument('--threads', type=int, default=8)
        if name == 'score':
            s.add_argument('--vcf', required=True)
            s.add_argument('--label', required=True)
            s.add_argument('--patch', required=True, help="the patched arm's patch.json")
            s.add_argument('--pad', type=int, default=100)
            s.add_argument('--threads', type=int, default=4)
            s.add_argument('--no-haps', action='store_true')
            s.add_argument('--force', action='store_true')
    a = p.parse_args(argv)
    {'gfa': cmd_gfa, 'patch': cmd_patch, 'reads': cmd_reads, 'map': cmd_map, 'call': cmd_call,
     'score': cmd_score}[a.cmd](a)


if __name__ == '__main__':
    main()
