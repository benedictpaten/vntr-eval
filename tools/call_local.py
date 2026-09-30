#!/usr/bin/env python3
"""call_local.py -- Stage 3: `vg call` on a region's local graph, with the reads re-mapped to it.

The Stage 3 call (FINAL_ARM = hybrid200k; `call ID --graph G` with no --arm) of one region on one
graph G (mc = regions/<id>/mc.gfa, a candidate = candidates/<G>/<id>.gfa, or truth):

  1. HYBRID GRAPH (build_hybrid). The genome-wide hap32 graph of the contig
     (work/wgs-mm095/<contig>/<contig>.gbz) over CHM13 [span_start-1-200 kb, span_end+200 kb]
     (gbz-base --interval, --context 5000, plus everything between the anchors), with the span --
     the anchors and everything between them -- replaced by G's graph. Every panel path is the GBZ's
     own path cut to the window (named_pieces: one `vg paths -A` pass per contig), and its span walk
     is G's path of the same hap32 name. So CHM13 and GRCh38 are represented exactly as genome-wide
     (GBWT tag reference_samples = "CHM13 GRCh38"; CHM13#0#<contig> and GRCh38#0#<contig>[offset]
     REFERENCE paths; the 32 sampled haplotypes HAPLOTYPE paths recombination#<k>#<contig>#<block>,
     fragments of one haplotype one panel member), the panel is the same 34 haplotypes, enumeration
     uses the same 32, and the flanking sites give the linkage layer its genome-wide context.
     Node IDs: flank nodes and mc's span nodes keep their genome-wide IDs; a candidate's span nodes
     are spread over the genome-wide span's ID range, so vg call's 4,096-ID depth-rate windows
     (allele_likelihood.cpp local_read_stats, RATE_WINDOW) mix flank and span much as they do
     genome-wide -- except where a candidate has more span nodes than the range has free IDs
     (build.json ids_interleaved false; 8 graph/region pairs of the pilot, all controls or the truth
     graph): its span IDs then follow the window's largest ID, in windows of mostly span. The call
     depends on this numbering through the depth term (null_decomp.py): in the pilot, re-laying the
     span IDs moved mc's hotspot haplotypes by up to 2,819 edits and mc_unchop's by up to 3,080, the
     three candidates' by at most 159 (results/stage3_noise.md). The truth graph's span is mafft(CHM13, HG002's truth haplotypes); the truth
     haplotypes (sample 'truthhap') carry CHM13's flank walks.
  2. READS (fetch_window_reads, map_hybrid). The reads the genome-wide mapping placed in the window
     (remap_local.fetch with flank 200 kb: rebuilt from GAF-Base with sequence and base qualities,
     mates together), mapped to the hybrid graph exactly as Stage 2 mapped its local graphs
     (vg autoindex -n -w sr-giraffe; giraffe paired with Stage 2's pooled fragment-length
     distribution, single-end without a mate; GAF).
  3. CALL with the pinned vg (work/bin/vg-2a6a228a5, `pin`) and EXACTLY the flags of the genome-wide
     short-read production run (PRODUCTION below; re-running that command on chr21 with the pinned
     binary reproduces work/wgs-mm095/chr21/chr21.vcf byte for byte), replacing only the graph and
     the read source (--gaf-reads <local GAF> for --gaf-base/--gbz-base). Ploidy as call_wgs.sh:
     -d 1 on chrY and on chrX outside the PARs.
  4. FINISH (finish_vcf). Shift to CHM13 coordinates (CHROM = contig, POS + the offset of the local
     CHM13 path's start, which is checked to spell CHM13 there), check every REF allele against the
     CHM13 FASTA, bgzip in vg's own record order (not bcftools sort: see finish_vcf) and index
     -> work/stage3/calls/<graph>/<id>.vcf.gz (+ .json with the command, .log, .raw.vcf).

Diagnostic arms (`call --arm A`, mc unless stated; -> work/stage3/diag/<graph>@<arm>/<id>.vcf.gz):
  stage2          the Stage 2 local graph (1 kb CHM13 flank, any graph) with the Stage 2 GAF
  stage2orig      the Stage 2 local mc graph with the genome-wide alignments projected into it
  gwsub<F>k       the genome-wide hap32 graph over span +- F kb with the genome-wide alignments
                  contained in it: calling locally with production's own alignments
  hybrid<F>k      the hybrid graph with re-mapped reads (any graph); hybrid200k is the Stage 3 call
  hybrid<F>korig  the hybrid mc graph with the genome-wide alignments (it keeps every genome-wide ID)
  hybrid<F>kids   the hybrid graph with the span's node IDs re-laid (half a spread step on; mc's
                  too), reads re-mapped (any graph)
  hybrid<F>kmcids the hybrid graph with the span's nodes on mc's native span IDs (any graph but mc)
  ...nodepth      any of these with --depth-term 0

Replicate arms (`run pilot --graphs '' --arms hybrid50k,hybrid200kids,hybrid200kmcids --arm-graphs G,...
--compact`): hybrid50k, hybrid200kids and hybrid200kmcids are neutral changes of the Stage 3 call -- the window size, the node
numbering -- run for EVERY graph, so that each graph's own noise can be measured without the truth
(replicate_noise.py). MC-only nulls cannot stand in for it: the candidates move where MC does not
(controls) and hardly move where MC does (node numbering).

Null graphs (`run pilot --graphs mc_relabel,mc_unchop`; -> work/stage3/diag/<null>@hybrid200k/<id>.vcf.gz):
MC's own alignment written the way a candidate is (NULL_GRAPHS, unchop_gfa), re-mapped and called like
any graph, to measure how far the Stage 3 call moves when the alignment does not change at all.

The GATE (`gate`) compares the arms of mc with the genome-wide production calls over each
region's anchor-to-anchor span: the two called haplotype sequences built by the Stage 3 scoring
rule (span_haplotypes; sequences, never records keyed by POS), a record-level decomposition keyed
by vg's snarl ID (same / phase only / genotype / one side only, and the DR ratio, i.e. the
depth-rate estimate), and truvari FP/FN from one span-restricted truvari run per VCF (the
genome-wide pipeline's parameters) plus the genome-wide truvari labels in the span.

Subcommands (from the repository root; VNTR_STAGE2_WORK or --stage2-work = remap_local.py's work
directory, which holds the Stage 2 local graphs, GAFs and fraglen.json):

    python3 tools/call_local.py pin                                  # copy + re-sign the pinned vg
    python3 tools/call_local.py call L014297 --graph mc [--arm A] [--threads 2]
    python3 tools/call_local.py run  pilot --graphs mc,unit_aware__all,mafft_linsi,unit_aware,truth \
                                     [--arms gwsub50k,hybrid200korig,...] --jobs 3 --threads 2
    python3 tools/call_local.py haps L014297 work/stage3/calls/mc/L014297.vcf.gz   # haplotypes + ED
    python3 tools/call_local.py gate pilot --tsv results/stage3_gate.tsv [--arms ...]
    python3 tools/call_local.py validate-scorer pilot                # truth VCF -> ED 0, perturbed -> ED > 0
    python3 tools/call_local.py refinish work/stage3/calls/mc/*.vcf.gz   # rebuild .vcf.gz from .raw.vcf

Everything written goes to $VNTR_WORK/stage3 (default work/stage3) and, for `gate`, --tsv/--json.
"""
import argparse
import collections
import concurrent.futures
import gzip
import itertools
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

REPO = config.REPO
REGIONS = config.REGIONS_DIR
STAGE3 = os.path.join(config.WORK_DIR, 'stage3')
STAGE2_WORK = os.environ.get('VNTR_STAGE2_WORK') or os.path.join(config.WORK_DIR, 'remap')

# ------------------------------------------------------------------ the caller and its flags

VG_SOURCE = os.path.join(config.HOME, 'CLionProjects/vg/bin/vg')      # branch snarl-tree-order
VG_COMMIT = '2a6a228a5'
VG_PINNED = os.environ.get('VNTR_STAGE3_VG') or os.path.join(REPO, 'work', 'bin', 'vg-' + VG_COMMIT)

# The genome-wide short-read production run whose calls are
# ~/PycharmProjects/vg-call-eval/work/wgs-mm095/<contig>/<contig>.vcf.gz.
#
# Found in ~/PycharmProjects/vg-call-eval/scripts/wgs/call_wgs.sh, function call_one_bed, launched per
# contig by scripts/wgs/schedule_wgs.py (--work work/wgs-mm095, --threads 5 default, --extra
# "--mismap-max 0.95": the run directory is named for it, and docs/tier2-parameters.md "Whole genome"
# describes work/wgs-mm095 as the `--mismap-max 0.95` arm on binary 8acbb43a2; 0cab3fbd4 later made
# 0.95 the default, so the pinned 2a6a228a5 would behave the same without the flag). The command, with
# REF_SAMPLE=CHM13, SAMPLE=HG002, READS_DB=work/reads.hap32.gaf.db, GRAPH_DB=work/graph.hap32.gbz.db:
#
#   vg call work/wgs-mm095/<C>/<C>.gbz -p CHM13#0#<C> -s HG002 -d <ploidy> -t 5 --progress \
#       --mismap-max 0.95 --read-likelihood --phased --mosaic-out <C>.mosaic.tsv \
#       [--ploidy-bed scripts/wgs/chrX.par.bed   (chrX only, with -d 1; chrY -d 1)] \
#       [-r <C>.snarls.pb   (only if prep made one; none exists for wgs-mm095, so snarls computed)] \
#       --gaf-base work/reads.hap32.gaf.db --gbz-base work/graph.hap32.gbz.db
#
# Checked, not assumed (2026-09-27): this command for chr21, run with the pinned binary at -t 4,
# reproduced work/wgs-mm095/chr21/chr21.vcf byte for byte (113,609 records, header included;
# 1,154 s, 6.2 GB peak). So the flags are the production flags and 2a6a228a5 calls short reads as
# the production binary did.
PRODUCTION_SOURCE = ('~/PycharmProjects/vg-call-eval/scripts/wgs/call_wgs.sh (call_one_bed), run by '
                     'scripts/wgs/schedule_wgs.py --work work/wgs-mm095 --extra "--mismap-max 0.95"')
PRODUCTION_FLAGS = ['-s', 'HG002', '--progress', '--mismap-max', '0.95', '--read-likelihood', '--phased']
PRODUCTION_READS = ['--gaf-base', 'work/reads.hap32.gaf.db', '--gbz-base', 'work/graph.hap32.gbz.db']
PRODUCTION_THREADS = 5
# chrX PARs as scripts/wgs/call_wgs.sh draws them (CHM13 chrX, 1-based): PAR1 <= 2394370, PAR2 >= 153926003
PAR = {'chrX': [(1, 2394370), (153926003, 10 ** 10)]}

MAX_NODE = 1024          # GBZ node length limit (gbwtgraph chops longer segments)
FLANK = 1000             # remap_local.py's CHM13 flank

PILOT = {
    'hotspot_vntr': ['L009656', 'L014297', 'L015415', 'L005990', 'L012184', 'L002013', 'L011138',
                     'L012272', 'L001909'],
    'hotspot_other': ['L016870'],
    'control_vntr_matched': ['L016124', 'L007687', 'L007172', 'L011430', 'L007065', 'L004021', 'L014111'],
    'control_nontr_sv': ['L000506', 'L006556'],
}
PILOT_IDS = [x for v in PILOT.values() for x in v]
STAGE3_GRAPHS = ['mc', 'unit_aware__all', 'mafft_linsi', 'unit_aware', 'truth']
TRUVARI = os.path.join(config.EVAL_DIR, 'work/truvari-venv/bin/truvari')

STEP_RE = re.compile(r'([<>])([^<>]+)')
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(s):
    return s.translate(COMP)[::-1]


def log(*a):
    print('[call_local %s]' % time.strftime('%H:%M:%S'), *a, file=sys.stderr, flush=True)


class CallError(RuntimeError):
    pass


def run(cmd, stdout=None, stderr=None, check=True, cwd=None, env=None):
    so = open(stdout, 'w') if isinstance(stdout, str) else stdout
    se = open(stderr, 'w') if isinstance(stderr, str) else stderr
    try:
        p = subprocess.run([str(c) for c in cmd], stdout=so if so else subprocess.PIPE,
                           stderr=se if se else subprocess.PIPE, text=True, cwd=cwd,
                           env=env or config.tool_env())
    finally:
        if isinstance(stdout, str):
            so.close()
        if isinstance(stderr, str):
            se.close()
    if check and p.returncode != 0:
        msg = p.stderr[-3000:] if isinstance(p.stderr, str) else ''
        if isinstance(stderr, str) and os.path.exists(stderr):
            with open(stderr) as f:
                msg = f.read()[-3000:]
        raise CallError('command failed (%d): %s\n%s' % (p.returncode, ' '.join(map(str, cmd)), msg))
    return p


# ------------------------------------------------------------------ regions

def region_ids(spec):
    if spec in ('pilot', 'all19'):
        return list(PILOT_IDS)
    out = []
    for x in spec.split(','):
        x = x.strip()
        if x in PILOT:
            out += PILOT[x]
        elif x:
            out.append(x)
    return out


def load_region(rid):
    with open(os.path.join(REGIONS, rid, 'region.json')) as f:
        d = json.load(f)
    d['dir'] = os.path.join(REGIONS, rid)
    return d


def ploidy_of(reg):
    """HG002 is male: chrY haploid, chrX haploid outside the PARs (call_wgs.sh's -d 1 + PAR BED)."""
    c = reg['contig']
    if c == 'chrY':
        return 1
    if c == 'chrX':
        for s, e in PAR['chrX']:
            if reg['span_start'] <= e and s <= reg['span_end']:
                return 2
        return 1
    return 2


def read_fasta(path):
    out, name, buf = [], None, []
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ''.join(buf).upper()))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        out.append((name, ''.join(buf).upper()))
    return out


def fetch_ref(contig, start0, end0):
    """CHM13 bases [start0, end0) (0-based half-open)."""
    fa = config.data_paths(contig)['ref_fa']
    p = run([config.SAMTOOLS, 'faidx', fa, '%s:%d-%d' % (contig, start0 + 1, end0)])
    return ''.join(p.stdout.split('\n', 1)[1].split()).upper()


def contig_length(contig):
    fai = config.data_paths(contig)['ref_fa'] + '.fai'
    with open(fai) as f:
        for line in f:
            x = line.split('\t')
            if x[0] == contig:
                return int(x[1])
    raise CallError('no %s in %s' % (contig, fai))


# ------------------------------------------------------------------ 1. the calling graph

def parse_local_gfa(path):
    nodes, edges, walks = collections.OrderedDict(), [], []
    with open(path) as f:
        for line in f:
            if not line or line[0] not in 'SLW':
                continue
            x = line.rstrip('\n').split('\t')
            if x[0] == 'S':
                nodes[x[1]] = x[2].upper()
            elif x[0] == 'L':
                edges.append((x[1], x[2], x[3], x[4]))
            else:
                walks.append({'sample': x[1], 'hap': x[2], 'seqid': x[3], 'start': int(x[4]),
                              'end': int(x[5]), 'steps': STEP_RE.findall(x[6])})
    return nodes, edges, walks


def walk_name_fields(orig, contig):
    """(sample, hap, seqid, start) of a W line for a hap32 path name, following the genome-wide GBZ:
    CHM13#0#chr6 -> CHM13 0 chr6 0; GRCh38#0#chr6[167641449] -> GRCh38 0 chr6 167641449;
    recombination#12#chr6#1 -> recombination 12 chr6 1 (gbwtgraph reads a haplotype W line's start
    as the phase block / fragment number). Anything else (the truth graph's HG002#1) keeps its
    sample and haplotype on the region's contig."""
    p = orig.split('#')
    if p[0] == 'CHM13':
        return 'CHM13', '0', contig, 0
    if p[0] == 'GRCh38':
        m = re.match(r'(.+)\[(\d+)\]$', p[2])
        return 'GRCh38', p[1], (m.group(1) if m else p[2]), (int(m.group(2)) if m else 0)
    if p[0] == 'recombination' and len(p) == 4:
        return 'recombination', p[1], p[2], int(p[3])
    return p[0], (p[1] if len(p) > 1 else '0'), contig, 0


def graph_dir(graph, rid):
    return os.path.join(STAGE3, 'graphs', graph, rid)


def prep(rid, graph, stage2=None, force=False):
    """Build <stage3>/graphs/<graph>/<id>/{call.gfa, call.gbz, reads.gaf, prep.json}."""
    stage2 = stage2 or STAGE2_WORK
    reg = load_region(rid)
    od = graph_dir(graph, rid)
    done = os.path.join(od, 'prep.json')
    if os.path.exists(done) and not force:
        with open(done) as f:
            return json.load(f)
    sd = os.path.join(stage2, rid, graph)
    src = os.path.join(sd, 'local.gfa')
    if graph == 'truth':
        raise CallError('the Stage 2 truth graph has no CHM13 path, so vg call cannot use CHM13 as its '
                        'reference there; calling it needs a CHM13 + HG002 graph and a re-map')
    if not os.path.exists(src) or not os.path.exists(os.path.join(sd, 'reads.gaf')):
        raise CallError('%s: no Stage 2 local graph/GAF for %s (%s)' % (rid, graph, sd))
    with open(os.path.join(sd, 'map.json')) as f:
        m = json.load(f)
    names = m['names']
    os.makedirs(od, exist_ok=True)
    contig = reg['contig']
    nodes, edges, walks = parse_local_gfa(src)
    # chop long segments with explicit new IDs
    top = max(int(n) for n in nodes)
    chop = {}
    new_nodes = collections.OrderedDict()
    for n, s in nodes.items():
        if len(s) <= MAX_NODE:
            new_nodes[n] = s
            continue
        pieces = []
        for i in range(0, len(s), MAX_NODE):
            top += 1
            new_nodes[str(top)] = s[i:i + MAX_NODE]
            pieces.append(str(top))
        chop[n] = pieces

    def first(n, o):
        pc = chop.get(n, [n])
        return pc[0] if o == '+' else pc[-1]

    def last(n, o):
        pc = chop.get(n, [n])
        return pc[-1] if o == '+' else pc[0]

    new_edges = set()
    for a, ao, b, bo in edges:
        new_edges.add((last(a, ao), ao, first(b, bo), bo))
    for n, pc in chop.items():
        for u, v in zip(pc, pc[1:]):
            new_edges.add((u, '+', v, '+'))

    def expand(steps):
        out = []
        for o, n in steps:
            pc = chop.get(n, [n])
            out += [(o, x) for x in (pc if o == '>' else pc[::-1])]
        return out

    ref_samples = set()
    wl = []
    seqs = {}
    chm13 = None
    for w in walks:
        local_name = '%s#%s#%s' % (w['sample'], w['hap'], w['seqid'])
        orig = names[local_name]
        sample, hap, seqid, start = walk_name_fields(orig, contig)
        if sample in ('CHM13', 'GRCh38'):
            ref_samples.add(sample)
        steps = expand(w['steps'])
        seq = ''.join(new_nodes[n] if o == '>' else rc(new_nodes[n]) for o, n in steps)
        old = ''.join(nodes[n] if o == '>' else rc(nodes[n]) for o, n in w['steps'])
        if seq != old:
            raise CallError('%s %s: chopped walk does not spell the Stage 2 walk' % (rid, orig))
        key = (sample, hap, seqid, start)
        if key in seqs:
            raise CallError('%s: two paths map to one W name %s' % (rid, key))
        seqs[key] = orig
        wl.append('W\t%s\t%s\t%s\t%d\t%d\t%s\n' % (sample, hap, seqid, start, start + len(seq),
                                                    ''.join(o + n for o, n in steps)))
        if sample == 'CHM13':
            chm13 = {'name': orig, 'seq': seq, 'steps': steps}
    if chm13 is None:
        raise CallError('%s %s: no CHM13 path' % (rid, graph))
    gfa = os.path.join(od, 'call.gfa')
    with open(gfa + '.tmp', 'w') as f:
        f.write('H\tVN:Z:1.0\tRS:Z:%s\n' % ' '.join(s for s in ('CHM13', 'GRCh38') if s in ref_samples))
        for n, s in new_nodes.items():
            f.write('S\t%s\t%s\n' % (n, s))
        for a, ao, b, bo in sorted(new_edges):
            f.write('L\t%s\t%s\t%s\t%s\t0M\n' % (a, ao, b, bo))
        f.writelines(wl)
    os.replace(gfa + '.tmp', gfa)
    gbz = os.path.join(od, 'call.gbz')
    run([VG_PINNED, 'gbwt', '-G', gfa, '--gbz-format', '-g', gbz], stderr=os.path.join(od, 'gbwt.log'))
    # checks: node IDs and sequences kept, reference paths as in the genome-wide GBZ
    conv = run([VG_PINNED, 'convert', '-f', gbz]).stdout
    got = {}
    for line in conv.splitlines():
        if line.startswith('S\t'):
            x = line.split('\t')
            got[x[1]] = x[2]
    if got != dict(new_nodes):
        raise CallError('%s %s: GBZ nodes differ from call.gfa (IDs were renumbered?)' % (rid, graph))
    tags = run([VG_PINNED, 'gbwt', '-Z', gbz, '--tags']).stdout
    meta = run([VG_PINNED, 'paths', '-x', gbz, '-M', '-S', 'CHM13']).stdout
    if ('CHM13#0#%s\tREFERENCE' % contig) not in meta:
        raise CallError('%s %s: CHM13#0#%s is not a REFERENCE path:\n%s' % (rid, graph, contig, meta))
    ref_tag = [l.split('\t', 1)[1] for l in tags.splitlines() if l.startswith('reference_samples')]
    # offset of the local CHM13 path on the contig
    lf = m['left_flank_node']
    left_len = len(nodes[lf])
    offset0 = reg['span_start'] - 1 - left_len
    ref = fetch_ref(contig, offset0, offset0 + len(chm13['seq']))
    if ref != chm13['seq']:
        raise CallError('%s %s: local CHM13 path does not spell CHM13 at offset %d' % (rid, graph, offset0))
    if chm13['steps'][0][1] not in (lf,) + tuple(chop.get(lf, [])):
        raise CallError('%s %s: CHM13 path does not start on the left flank node' % (rid, graph))
    # translate the GAF
    gaf_out = os.path.join(od, 'reads.gaf')
    stats = translate_gaf(os.path.join(sd, 'reads.gaf'), gaf_out, new_nodes, chop)
    info = {'region_id': rid, 'graph': graph, 'contig': contig, 'stage2_dir': sd,
            'source_gfa': m.get('source_gfa'), 'nodes': len(new_nodes), 'segments': len(nodes),
            'segments_chopped': len(chop), 'edges': len(new_edges), 'paths': len(wl),
            'reference_samples_tag': ref_tag, 'chm13_path_len': len(chm13['seq']),
            'left_flank_node': lf, 'left_flank_len': left_len,
            'right_flank_len': len(nodes[m['right_flank_node']]),
            'offset0': offset0, 'gaf': stats, 'vg': VG_PINNED}
    with open(done, 'w') as f:
        json.dump(info, f, indent=1)
    return info


def translate_gaf(src, dst, nodes, chop):
    """Segment-space GAF -> node-space GAF (chopped segments expanded; offsets unchanged)."""
    st = collections.Counter()
    with open(src) as f, open(dst + '.tmp', 'w') as fo:
        for line in f:
            if line.startswith('@') or not line.strip():
                fo.write(line)
                continue
            x = line.rstrip('\n').split('\t')
            st['alignments'] += 1
            if x[5] == '*':
                st['unmapped'] += 1
                fo.write(line)
                continue
            steps = STEP_RE.findall(x[5])
            out = []
            for o, n in steps:
                pc = chop.get(n, [n])
                if len(pc) > 1:
                    st['steps_chopped'] += 1
                out += [(o, p) for p in (pc if o == '>' else pc[::-1])]
            plen = sum(len(nodes[n]) for _, n in out)
            if plen != int(x[6]):
                raise CallError('GAF path length mismatch after translation: %s' % x[0])
            x[5] = ''.join(o + n for o, n in out)
            fo.write('\t'.join(x) + '\n')
    os.replace(dst + '.tmp', dst)
    return dict(st)


# ------------------------------------------------------------------ the genome-wide graph around the span (mc@gwsub<F>k)

def window_query(rid, flank, context=5000):
    """gbz-base query of the hap32 graph over CHM13 [span_start-1-flank, span_end+flank) with
    --context, and the GAF-Base alignments contained in it -> <graphs>/gwsub<F>k/<id>/{q.gfa,reads.gaf}."""
    reg = load_region(rid)
    od = graph_dir('gwsub%dk' % (flank // 1000), rid)
    os.makedirs(od, exist_ok=True)
    contig = reg['contig']
    dp = config.data_paths(contig)
    lo = max(0, reg['span_start'] - 1 - flank)
    hi = min(contig_length(contig), reg['span_end'] + flank)
    qgfa, qgaf = os.path.join(od, 'q.gfa'), os.path.join(od, 'q.gaf')
    if not (os.path.exists(qgfa) and os.path.exists(qgaf)):
        run([config.GBZ_BASE, 'query', '--sample', 'CHM13', '--contig', contig, '--interval', '%d..%d' % (lo, hi),
             '--context', str(context), '--gaf-base', dp['gaf_db'], '--gaf-output', qgaf + '.tmp',
             '--alignments', 'contained', dp['gbz_db']], stdout=qgfa + '.tmp', stderr=os.path.join(od, 'query.log'))
        os.replace(qgaf + '.tmp', qgaf)
        os.replace(qgfa + '.tmp', qgfa)
    # The window's nodes: the interval query's plus everything between the anchors. The --context
    # of an interval query can fall short inside long parallel paths (L012272, L011138), which would
    # cut haplotypes in the middle of the span.
    wgfa = os.path.join(od, 'window.gfa')
    if not os.path.exists(wgfa) or not os.path.exists(os.path.join(od, 'reads.gaf')):
        bn = os.path.join(od, 'between.gfa')
        bg = os.path.join(od, 'between.gaf')
        if not (os.path.exists(bn) and os.path.exists(bg)):
            run([config.GBZ_BASE, 'query', '--between', '%s:%s' % (reg['anchor_left'], reg['anchor_right']),
                 '--limit', '2000000', '--gaf-base', dp['gaf_db'], '--gaf-output', bg + '.tmp',
                 '--alignments', 'overlapping', dp['gbz_db']], stdout=bn + '.tmp')
            os.replace(bg + '.tmp', bg)
            os.replace(bn + '.tmp', bn)
        qn, qe, _ = parse_local_gfa(qgfa)
        bnodes, bedges, _ = parse_local_gfa(bn)
        added = [n for n in bnodes if n not in qn]
        # the genome-wide alignments: contained in the interval subgraph, plus those overlapping the
        # between-anchors subgraph that lie wholly in the window; one copy of each (name, mate, path, start)
        wn = set(qn) | set(bnodes)
        seen = set()
        with open(os.path.join(od, 'reads.gaf.tmp'), 'w') as fo:
            for src in (qgaf, bg):
                with open(src) as f:
                    for line in f:
                        if line.startswith('@'):
                            continue
                        x = line.rstrip('\n').split('\t')
                        mate = 'fn' if '\tfn:Z:' in line else 'fp' if '\tfp:Z:' in line else ''
                        key = (x[0], mate, x[5], x[7])
                        if key in seen:
                            continue
                        if x[5] != '*' and any(n not in wn for _, n in STEP_RE.findall(x[5])):
                            continue
                        seen.add(key)
                        fo.write(line)
        os.replace(os.path.join(od, 'reads.gaf.tmp'), os.path.join(od, 'reads.gaf'))
        with open(wgfa + '.tmp', 'w') as f:
            f.write('H\tVN:Z:1.0\tXA:i:%d\n' % len(added))
            for n, x in list(qn.items()) + [(n, bnodes[n]) for n in added]:
                f.write('S\t%s\t%s\n' % (n, x))
            for e in sorted(set(qe) | set(bedges)):
                f.write('L\t%s\t%s\t%s\t%s\t0M\n' % e)
        os.replace(wgfa + '.tmp', wgfa)
    return od, lo, hi


def window_added_nodes(od):
    with open(os.path.join(od, 'window.gfa')) as f:
        return int(f.readline().rstrip('\n').split('XA:i:')[1])


def named_pieces(rids, flanks, context=5000):
    """The hap32 GBZ's own haplotype and GRCh38 paths, cut to each region's window subgraph (maximal
    runs of steps on its nodes), with their names -- gbz-base anonymises haplotypes, and a haplotype
    that enters a window twice (a GBWT fragment boundary, or a stray visit) must stay ONE panel member,
    as it is genome-wide. One streaming pass of `vg paths -A` per contig. -> <od>/pieces.json"""
    if isinstance(flanks, int):
        flanks = [flanks]
    todo = collections.defaultdict(list)
    for rid in rids:
        for flank in flanks:
            od, _, _ = window_query(rid, flank, context)
            if not os.path.exists(os.path.join(od, 'pieces.json')):
                todo[load_region(rid)['contig']].append((od, od))
    for contig, rs in todo.items():
        log('named_pieces: streaming %s paths for %d windows' % (contig, len(rs)))
        nsets = {}
        for rid, od in rs:
            nodes, _, _ = parse_local_gfa(os.path.join(od, 'window.gfa'))
            nsets[rid] = set(nodes)
        pieces = {rid: [] for rid, _ in rs}
        gbz = config.data_paths(contig)['contig_gbz']
        for sel in (['-S', 'recombination'], ['-Q', 'GRCh38#0#%s' % contig]):
            proc = subprocess.Popen([VG_PINNED, 'paths', '-x', gbz, '-A'] + sel, stdout=subprocess.PIPE,
                                    stderr=subprocess.DEVNULL, text=True, env=config.tool_env())
            for line in proc.stdout:
                if line.startswith('@'):
                    continue
                x = line.split('\t', 6)
                steps = STEP_RE.findall(x[5])
                for rid, _ in rs:
                    ns = nsets[rid]
                    idx = [i for i, st in enumerate(steps) if st[1] in ns]
                    cur = []
                    for i in idx:
                        if cur and i != cur[-1] + 1:
                            pieces[rid].append((x[0], [steps[j] for j in cur]))
                            cur = []
                        cur.append(i)
                    if cur:
                        pieces[rid].append((x[0], [steps[j] for j in cur]))
            if proc.wait() != 0:
                raise CallError('vg paths -A %s failed on %s' % (' '.join(sel), gbz))
        for rid, od in rs:
            with open(os.path.join(od, 'pieces.json'), 'w') as f:
                json.dump(pieces[rid], f)


def prep_gwsub(rid, flank, context=5000, force=False):
    """The genome-wide hap32 graph over CHM13 [span_start - flank, span_end + flank] (gbz-base
    --interval with --context), with its node IDs; CHM13 its REFERENCE path CHM13#0#<contig> (start 0,
    shifted afterwards); GRCh38 and the 32 sampled haplotypes as the GBZ's own paths cut to the window
    (named_pieces), GRCh38 a REFERENCE sample; the genome-wide alignments contained in it (GAF-Base)."""
    reg = load_region(rid)
    graph = 'gwsub%dk' % (flank // 1000)
    od, lo, hi = window_query(rid, flank, context)
    done = os.path.join(od, 'prep.json')
    if os.path.exists(done) and not force:
        with open(done) as f:
            return json.load(f)
    named_pieces([rid], flank, context)
    contig = reg['contig']
    _, _, walks = parse_local_gfa(os.path.join(od, 'q.gfa'))
    nodes, edges, _ = parse_local_gfa(os.path.join(od, 'window.gfa'))
    with open(os.path.join(od, 'pieces.json')) as f:
        pieces = json.load(f)
    chm13 = offset0 = None
    wl = []
    for w in walks:
        if w['sample'] == 'CHM13':
            steps = [tuple(s) for s in w['steps']]
            offset0 = w['start']
            chm13 = ''.join(nodes[n] if o == '>' else rc(nodes[n]) for o, n in steps)
            wl.append('W\tCHM13\t0\t%s\t0\t%d\t%s\n' % (contig, len(chm13), ''.join(o + n for o, n in steps)))
    if chm13 is None:
        raise CallError('%s: no CHM13 walk in the window query' % rid)
    ref = fetch_ref(contig, offset0, offset0 + len(chm13))
    if ref != chm13:
        raise CallError('%s: window CHM13 walk does not spell CHM13 at %d' % (rid, offset0))
    per_path = collections.Counter()
    samples = collections.Counter()
    haps = set()
    for name, steps in pieces:
        seq_len = sum(len(nodes[n]) for _, n in steps)
        i = per_path[name]
        per_path[name] += 1
        p = name.split('#')
        if p[0] == 'GRCh38':
            m = re.match(r'(.+)\[(\d+)\]$', p[2])
            start = (int(m.group(2)) if m else 0) + i
            wl.append('W\tGRCh38\t0\t%s\t%d\t%d\t%s\n' % (contig, start, start + seq_len,
                                                          ''.join(o + n for o, n in steps)))
            samples['GRCh38'] += 1
            haps.add(('GRCh38', '0'))
        else:
            start = int(p[3]) * 10 ** 6 + i
            wl.append('W\trecombination\t%s\t%s\t%d\t%d\t%s\n' % (p[1], contig, start, start + seq_len,
                                                                 ''.join(o + n for o, n in steps)))
            samples['recombination'] += 1
            haps.add(('recombination', p[1]))
    # the named pieces must be gbz-base's (anonymous) walks, up to orientation: gbz-base writes every
    # walk in CHM13's direction, the GBZ stores some haplotypes reversed, and the named pieces keep
    # the stored orientation. A walk left over is another visit of CHM13 to the window (a nearby
    # duplication), which gbz-base writes as an anonymous walk and which this graph leaves out.
    def canon(st):
        st = tuple(tuple(x) for x in st)
        rv = tuple(('<' if o == '>' else '>', n) for o, n in reversed(st))
        return min(st, rv)
    anon = collections.Counter(canon(w['steps']) for w in walks if w['sample'] != 'CHM13')
    named = collections.Counter(canon(st) for _, st in pieces)
    reversed_pieces = sum(1 for _, st in pieces
                          if tuple(tuple(x) for x in st) != canon(st))
    gfa = os.path.join(od, 'call.gfa')
    with open(gfa, 'w') as f:
        f.write('H\tVN:Z:1.0\tRS:Z:CHM13%s\n' % (' GRCh38' if samples['GRCh38'] else ''))
        for n, s in nodes.items():
            f.write('S\t%s\t%s\n' % (n, s))
        for a, ao, b, bo in edges:
            f.write('L\t%s\t%s\t%s\t%s\t0M\n' % (a, ao, b, bo))
        f.writelines(wl)
    gbz = os.path.join(od, 'call.gbz')
    run([VG_PINNED, 'gbwt', '-G', gfa, '--gbz-format', '-g', gbz], stderr=os.path.join(od, 'gbwt.log'))
    conv = run([VG_PINNED, 'convert', '-f', gbz]).stdout
    got = {l.split('\t')[1]: l.split('\t')[2] for l in conv.splitlines() if l.startswith('S\t')}
    if got != dict(nodes):
        raise CallError('%s: GBZ nodes differ from the window GFA' % rid)
    with open(os.path.join(od, 'reads.gaf')) as f:
        nreads = sum(1 for l in f if not l.startswith('@'))
    info = {'region_id': rid, 'graph': graph, 'contig': contig, 'window0': [lo, hi], 'context': context,
            'nodes': len(nodes), 'panel_haplotypes': len(haps) + 1, 'gbwt_paths': len(wl),
            'pieces': dict(samples), 'paths_in_several_pieces': sum(1 for v in per_path.values() if v > 1),
            'window_nodes_added_between_anchors': window_added_nodes(od),
            'pieces_equal_gbz_base_walks': anon == named,
            'gbz_base_walks_unnamed': sum((anon - named).values()),
            'named_pieces_not_in_gbz_base': sum((named - anon).values()),
            'pieces_stored_reversed': reversed_pieces,
            'chm13_path_len': len(chm13), 'offset0': offset0, 'reads': nreads}
    with open(done, 'w') as f:
        json.dump(info, f, indent=1)
    return info


# ------------------------------------------------------------------ hybrid graphs: the genome-wide flank around a graph's span

TRUTH_SAMPLE = 'truthhap'    # HG002's truth haplotypes in the truth graph (not 'HG002', the sample called)


def truth_span_gfa(reg):
    """The truth graph's span: CHM13 and HG002's truth haplotypes (one outside the PARs of chrX/Y),
    aligned with mafft (G-INS-i up to 15 kb, else FFT-NS-2, as Stage 2's truth graph) and induced
    with msa_graph.py; paths CHM13#0#<contig>, truthhap#1, truthhap#2. The Stage 2 truth graph
    lacks CHM13, so vg call could not use CHM13 as its reference there."""
    import realign
    import msa_graph
    od = os.path.join(STAGE3, 'truth_span', reg['region_id'])
    out = os.path.join(od, 'truth.gfa')
    if os.path.exists(out):
        return out
    os.makedirs(od, exist_ok=True)
    hap = read_fasta(os.path.join(reg['dir'], 'hap32.fa'))
    chm = [(n, x) for n, x in hap if n.startswith('CHM13')][0]
    t = truth_haps(reg)
    fa = os.path.join(od, 'seqs.fa')
    with open(fa, 'w') as f:
        f.write('>%s\n%s\n' % chm)
        for i, x in enumerate(t, 1):
            f.write('>%s#%d\n%s\n' % (TRUTH_SAMPLE, i, x))
    method = 'mafft_ginsi' if max(len(chm[1]), max(len(x) for x in t)) <= 15000 else 'mafft_fftns2'
    msa = os.path.join(od, 'seqs.msa.fa')
    os.makedirs(os.path.join(od, 'mafft'), exist_ok=True)
    info = realign.align_fasta(method, fa, msa, threads=2, timeout=3600, mem_mb=10000,
                               workdir=os.path.join(od, 'mafft'))
    if info['status'] != 'ok':
        raise CallError('%s: truth span alignment failed: %s' % (reg['region_id'], info['status']))
    msa_graph.msa_to_gfa(msa, fa, out + '.tmp', engine='native')
    os.replace(out + '.tmp', out)
    return out


def span_gfa(reg, graph):
    if graph in ('mc', 'mc_relabel'):
        return os.path.join(reg['dir'], 'mc.gfa')
    if graph == 'mc_unchop':
        return unchop_gfa(reg)
    if graph == 'truth':
        return truth_span_gfa(reg)
    return os.path.join(config.CANDIDATES_DIR, graph, reg['region_id'] + '.gfa')


# NULL GRAPHS: MC's own alignment, written the way a candidate is written, to measure how far the
# Stage 3 call moves when the alignment does not change at all. Both go through build_hybrid's
# candidate branch (span node IDs spread over the free IDs of the genome-wide span's range,
# flank-to-span edges from the path junctions) and are re-mapped and called like any graph.
#   mc_relabel  mc.gfa itself: only the node IDs (and hence the index, the depth-rate windows and
#               any ID-ordered tie-breaks) differ from the mc call
#   mc_unchop   mc.gfa with every non-branching run of nodes merged (unchop_gfa): the same graph
#               and paths with different node boundaries, as a coarser MSA-induced graph has
# Their calls go to <stage3>/diag/<graph>@hybrid200k/, not calls/.
NULL_GRAPHS = ('mc_relabel', 'mc_unchop')


def unchop_gfa(reg):
    """mc.gfa with every non-branching run of nodes merged into one node (vg mod -u's rule): a->b is
    merged when a's right side has one edge, to b's left side, which has no other, and no path starts
    or ends at that junction. Merged node ID = the run's first ID; paths are the same walks and are
    checked to spell the same sequences. -> <stage3>/null/mc_unchop/<id>.gfa (P lines)"""
    import remap_local as rl
    od = os.path.join(STAGE3, 'null', 'mc_unchop')
    out = os.path.join(od, reg['region_id'] + '.gfa')
    if os.path.exists(out):
        return out
    os.makedirs(od, exist_ok=True)
    nodes, edges, paths = rl.parse_gfa(os.path.join(reg['dir'], 'mc.gfa'))

    def exit_side(n, o):
        return (n, 'R' if o == '+' else 'L')

    def entry_side(n, o):
        return (n, 'L' if o == '+' else 'R')
    adj = collections.defaultdict(set)
    for a, ao, b, bo in edges:
        s, t = exit_side(a, ao), entry_side(b, bo)
        adj[s].add(t)
        adj[t].add(s)
    termini = set()
    for _, st in paths:
        termini.add(entry_side(*st[0]))
        termini.add(exit_side(*st[-1]))
    nxt, prv = {}, {}
    for a in nodes:
        s = (a, 'R')
        if len(adj[s]) != 1 or s in termini:
            continue
        (b, side), = adj[s]
        if side != 'L' or b == a or len(adj[(b, 'L')]) != 1 or (b, 'L') in termini:
            continue
        nxt[a] = b
        prv[b] = a
    chain_of = {}      # node -> (chain id, is head, is tail)
    new_nodes = collections.OrderedDict()
    for a in sorted(nodes, key=int):
        if a in prv or a in chain_of:
            continue
        run_ = [a]
        while run_[-1] in nxt:
            run_.append(nxt[run_[-1]])
        cid = run_[0]
        new_nodes[cid] = ''.join(nodes[x] for x in run_)
        for i, x in enumerate(run_):
            chain_of[x] = (cid, i == 0, i == len(run_) - 1)
    for a in nodes:                     # a cycle of merges has no head: leave its nodes unmerged
        if a not in chain_of:
            new_nodes[a] = nodes[a]
            chain_of[a] = (a, True, True)
            nxt.pop(a, None)
    new_edges = set()
    for a, ao, b, bo in edges:
        if nxt.get(a) == b and ao == '+' and bo == '+':
            continue
        if nxt.get(b) == a and ao == '-' and bo == '-':
            continue
        new_edges.add((chain_of[a][0], ao, chain_of[b][0], bo))
    new_paths = []
    for name, st in paths:
        w = [(chain_of[n][0], o) for n, o in st if (chain_of[n][1] if o == '+' else chain_of[n][2])]
        if rl.walk_seq(new_nodes, w) != rl.walk_seq(nodes, st):
            raise CallError('%s: unchopped path %s does not spell mc.gfa\'s' % (reg['region_id'], name))
        new_paths.append((name, w))
    with open(out + '.tmp', 'w') as f:
        f.write('H\tVN:Z:1.0\n')
        for n, s in new_nodes.items():
            f.write('S\t%s\t%s\n' % (n, s))
        for a, ao, b, bo in sorted(new_edges):
            f.write('L\t%s\t%s\t%s\t%s\t0M\n' % (a, ao, b, bo))
        for name, w in new_paths:
            f.write('P\t%s\t%s\t*\n' % (name, ','.join(n + o for n, o in w)))
    os.replace(out + '.tmp', out)
    return out


def between_nodes(reg, od):
    """Node set of the hap32 graph between the two anchors, anchors included (gbz-base --between)."""
    p = os.path.join(od, 'between.gfa')
    if not os.path.exists(p):
        L, R = reg['anchor_left'], reg['anchor_right']
        run([config.GBZ_BASE, 'query', '--between', '%s:%s' % (L, R), '--limit', '2000000',
             config.data_paths(reg['contig'])['gbz_db']], stdout=p + '.tmp')
        os.replace(p + '.tmp', p)
    nodes, _, _ = parse_local_gfa(p)
    return set(nodes)


def hybrid_dir(graph, rid, flank, tag=''):
    return os.path.join(STAGE3, 'hybrid%dk%s' % (flank // 1000, tag), graph, rid)


# vg call's depth-rate windows are blocks of RATE_WINDOW node IDs (allele_likelihood.cpp,
# local_read_stats). The replicate arm hybrid<F>kids re-lays the span's node IDs with this offset
# (a fraction of the spread step, or of a window where the IDs cannot interleave): a neutral change
# of the graph's numbering that moves which reads and bases share a window.
RATE_WINDOW = 4096
TAG_ID_OFFSET = {'': 0.0, 'ids': 0.5, 'mcids': 0.0}
ID_TAGS = ('ids', 'mcids')
HYBRID_ARM_RE = re.compile(r'hybrid(\d+)k(orig|ids|mcids)?$')


def parse_hybrid_setup(setup):
    """'hybrid200k' -> (200000, ''), 'hybrid50kids' -> (50000, 'ids'), 'hybrid200korig' -> (200000, 'orig'),
    'hybrid200kmcids' -> (200000, 'mcids')."""
    m = HYBRID_ARM_RE.match(setup)
    if not m:
        raise CallError('not a hybrid setup: %s' % setup)
    return int(m.group(1)) * 1000, m.group(2) or ''


def build_hybrid(rid, graph, flank, context=5000, force=False, tag=''):
    """The genome-wide hap32 graph over the window [span_start-1-flank, span_end+flank) with everything
    between the anchors (anchors included) replaced by `graph`'s span graph (mc.gfa or a candidate).

    Every panel path keeps its genome-wide flank walk (named_pieces) and its span walk is the graph's
    path of the same hap32 name, so the panel, the flank variation, the linkage context and CHM13's
    and GRCh38's reference status are the genome-wide run's. Pieces that enter the span but are not
    hap32 spanning paths (fragments ending inside the span) keep only their flank parts.
    Node IDs: flank nodes keep their genome-wide IDs, and so do mc's span nodes; a candidate's span
    nodes (topological order) are spread over the free IDs of the genome-wide span's range
    [anchor_left, anchor_right] in the same order, so vg call's 4,096-ID depth-rate windows mix flank
    and span much as they do genome-wide. Where the candidate has more span nodes than there are
    free IDs in that range (build.json ids_interleaved false: 8 graph/region pairs of the pilot,
    all controls or the truth graph), its span IDs follow the window's largest ID instead, and the
    depth rate there comes from windows that hold mostly span. Segments > 1,024 bp are chopped first.
    tag 'ids' (the replicate arm hybrid<F>kids) spreads the span IDs -- mc's too -- with an offset of
    TAG_ID_OFFSET['ids'] of a step (of a window, where they cannot interleave): the same graph with
    different numbering. tag 'mcids' (hybrid<F>kmcids, not for mc, whose IDs these are) puts the
    graph's span nodes, in order, on mc's own native span IDs (subsampled evenly), so the span has
    mc's ID density: the layout under which mc's calls are made. A graph with more span nodes than
    mc has cannot be laid out so (CallError).
    -> <stage3>/hybrid<F>k<tag>/<graph>/<id>/{hybrid.gfa, call.gbz, build.json}"""
    import remap_local as rl
    reg = load_region(rid)
    od = hybrid_dir(graph, rid, flank, tag)
    id_offset = TAG_ID_OFFSET[tag]
    done = os.path.join(od, 'build.json')
    if os.path.exists(done) and not force:
        with open(done) as f:
            return json.load(f)
    os.makedirs(od, exist_ok=True)
    contig = reg['contig']
    wd, lo, hi = window_query(rid, flank, context)
    named_pieces([rid], flank, context)
    _, _, wwalks = parse_local_gfa(os.path.join(wd, 'q.gfa'))
    wnodes, wedges, _ = parse_local_gfa(os.path.join(wd, 'window.gfa'))
    with open(os.path.join(wd, 'pieces.json')) as f:
        pieces = [(n, [tuple(s) for s in st]) for n, st in json.load(f)]
    inside = between_nodes(reg, od)
    if reg['anchor_left'][-1] != '+' or reg['anchor_right'][-1] != '+':
        raise CallError('%s: anchors on the reverse strand are not handled' % rid)
    AL, AR = reg['anchor_left'][:-1], reg['anchor_right'][:-1]
    snodes, sedges, spaths = rl.parse_gfa(span_gfa(reg, graph))
    hap = dict(read_fasta(os.path.join(reg['dir'], 'hap32.fa')))
    if graph == 'truth':
        # the span panel is CHM13 and HG002's truth haplotypes; the 32 sampled haplotypes and GRCh38
        # keep only their flank walks
        hap = {n: x for n, x in hap.items() if n.startswith('CHM13')}
        for i, x in enumerate(truth_haps(reg), 1):
            hap['%s#%d' % (TRUTH_SAMPLE, i)] = x
    spath = {}
    for name, steps in spaths:
        if name in hap:
            if rl.walk_seq(snodes, steps) != hap[name]:
                raise CallError('%s %s: span path %s does not spell hap32.fa' % (rid, graph, name))
            spath[name] = [('>' if o == '+' else '<', n) for n, o in steps]
    missing = [n for n in hap if n not in spath]
    if missing:
        raise CallError('%s %s: span graph lacks hap32 paths %s' % (rid, graph, missing[:3]))
    flank_nodes = collections.OrderedDict((n, s) for n, s in wnodes.items() if n not in inside)
    order = sorted(snodes, key=int)
    seg_pieces = {n: [snodes[n][i:i + MAX_NODE] for i in range(0, len(snodes[n]), MAX_NODE)] for n in order}
    n_new = sum(len(v) for v in seg_pieces.values())
    new_nodes = collections.OrderedDict()
    chop = {}
    ids_interleaved = True
    if graph == 'mc' and id_offset == 0 and all(len(v) == 1 for v in seg_pieces.values()):
        if any(n in flank_nodes for n in order):
            raise CallError('%s: mc span node IDs clash with flank nodes' % rid)
        for n in order:
            new_nodes[n] = snodes[n]
            chop[n] = [n]
    else:
        taken = {int(n) for n in flank_nodes}
        free = [i for i in range(int(AL), int(AR) + 1) if i not in taken]
        if tag == 'mcids':
            if graph == 'mc':
                raise CallError('tag mcids is mc\'s own layout; the mc call is the native one')
            mc_ids = sorted(int(n) for n in rl.parse_gfa(os.path.join(reg['dir'], 'mc.gfa'))[0])
            if n_new > len(mc_ids):
                raise CallError('%s %s: %d span nodes, more than mc\'s %d native span IDs' % (rid, graph, n_new, len(mc_ids)))
            if any(i in taken for i in mc_ids):
                raise CallError('%s: mc span node IDs clash with flank nodes' % rid)
            slots = [mc_ids[int(k * len(mc_ids) / float(n_new))] for k in range(n_new)]
        elif len(free) >= n_new:
            step = len(free) / float(n_new)
            slots = [free[min(len(free) - 1, int(k * step + id_offset * step))] for k in range(n_new)]
        else:
            ids_interleaved = False
            top = max(int(n) for n in wnodes) + 1 + int(id_offset * RATE_WINDOW)
            slots = list(range(top, top + n_new))
        if len(set(slots)) != n_new:
            raise CallError('%s %s: span node ID slots collide' % (rid, graph))
        k = 0
        for n in order:
            chop[n] = []
            for s in seg_pieces[n]:
                nid = str(slots[k])
                k += 1
                new_nodes[nid] = s
                chop[n].append(nid)

    def expand(steps):
        out = []
        for o, n in steps:
            pc = chop[n]
            out += [(o, x) for x in (pc if o == '>' else pc[::-1])]
        return out

    edges = set()
    for a, ao, b, bo in wedges:
        # flank edges; for mc (genome-wide IDs kept) also the window's flank-to-span edges
        if (a in flank_nodes and b in flank_nodes) or \
                (graph == 'mc' and chop.get(a) == [a] and b in flank_nodes) or \
                (graph == 'mc' and a in flank_nodes and chop.get(b) == [b]):
            edges.add((a, ao, b, bo))
    for a, ao, b, bo in sedges:
        la = chop[a][-1] if ao == '+' else chop[a][0]
        fb = chop[b][0] if bo == '+' else chop[b][-1]
        edges.add((la, ao, fb, bo))
    for n, pc in chop.items():
        for u, v in zip(pc, pc[1:]):
            edges.add((u, '+', v, '+'))
    all_nodes = collections.OrderedDict(list(flank_nodes.items()) + list(new_nodes.items()))

    def seq_of(steps):
        return ''.join(all_nodes[n] if o == '>' else rc(all_nodes[n]) for o, n in steps)

    def flip(steps):
        return [('<' if o == '>' else '>', n) for o, n in reversed(steps)]

    def split_piece(steps):
        """('span', (left, right)) for a walk crossing the span anchor to anchor in CHM13's direction
        (callers flip a reversed walk first), ('flank', [walk]) for one outside it, else
        ('truncated', [its parts outside the span])."""
        ins = [n in inside for _, n in steps]
        if not any(ins):
            return 'flank', [steps]
        i = ins.index(True)
        j = len(ins) - ins[::-1].index(True)
        if all(ins[i:j]) and steps[i] == ('>', AL) and steps[j - 1] == ('>', AR):
            return 'span', (steps[:i], steps[j:])
        parts, cur = [], []
        for st, x in zip(steps, ins):
            if x:
                if cur:
                    parts.append(cur)
                cur = []
            else:
                cur.append(st)
        if cur:
            parts.append(cur)
        return 'truncated', parts

    wl = []
    stats = collections.Counter()
    used_span = collections.Counter()
    chm = [w for w in wwalks if w['sample'] == 'CHM13']
    if len(chm) != 1:
        raise CallError('%s: %d CHM13 walks in the window' % (rid, len(chm)))
    offset0 = chm[0]['start']
    kind, parts = split_piece([tuple(s) for s in chm[0]['steps']])
    chm_name = [n for n in hap if n.startswith('CHM13')][0]
    if kind != 'span':
        raise CallError('%s: CHM13 does not cross the span anchor to anchor' % rid)
    chm_parts = parts
    chm_steps = parts[0] + expand(spath[chm_name]) + parts[1]
    chm_seq = seq_of(chm_steps)
    if fetch_ref(contig, offset0, offset0 + len(chm_seq)) != chm_seq:
        raise CallError('%s %s: hybrid CHM13 path does not spell CHM13 at %d' % (rid, graph, offset0))
    wl.append('W\tCHM13\t0\t%s\t0\t%d\t%s\n' % (contig, len(chm_seq), ''.join(o + n for o, n in chm_steps)))
    used_span[chm_name] += 1
    per_path = collections.Counter()
    grch38 = False
    for name, steps in pieces:
        rev = ('<', AL) in steps or ('<', AR) in steps
        kind, parts = split_piece(flip(steps) if rev else steps)
        if rev:
            stats['stored_reversed'] += 1
        if kind == 'span' and name in spath:
            w_ = parts[0] + expand(spath[name]) + parts[1]
            outs = [flip(w_) if rev else w_]      # keep the orientation the GBZ stores
            used_span[name] += 1
            stats['spanning'] += 1
        elif kind == 'span':
            outs = [(flip(p) if rev else p) for p in parts if p]
            stats['crossing_not_in_hap32'] += 1
        elif kind == 'truncated':
            outs = [(flip(p) if rev else p) for p in parts if p]
            stats['truncated_at_span'] += 1
        else:
            outs = parts
            stats['flank_only'] += 1
        p = name.split('#')
        for st in outs:
            i = per_path[name]
            per_path[name] += 1
            ln = len(seq_of(st))
            if p[0] == 'GRCh38':
                m = re.match(r'(.+)\[(\d+)\]$', p[2])
                start = (int(m.group(2)) if m else 0) + i
                wl.append('W\tGRCh38\t0\t%s\t%d\t%d\t%s\n' % (contig, start, start + ln, ''.join(o + n for o, n in st)))
                grch38 = True
            else:
                start = int(p[3]) * 10 ** 6 + i
                wl.append('W\trecombination\t%s\t%s\t%d\t%d\t%s\n' % (p[1], contig, start, start + ln,
                                                                     ''.join(o + n for o, n in st)))
    for name in sorted(hap):
        if name.startswith(TRUTH_SAMPLE + '#'):
            # HG002's truth haplotypes carry CHM13's flank walks, so that they are panel members
            # continuous across the span's edges as every other panel member is (without flanks the
            # linkage layer could only cross the span edge on CHM13, the one continuous path)
            st = chm_parts[0] + expand(spath[name]) + chm_parts[1]
            wl.append('W\t%s\t%s\t%s\t0\t%d\t%s\n' % (TRUTH_SAMPLE, name.split('#')[1], contig, len(seq_of(st)),
                                                          ''.join(o + n for o, n in st)))
            used_span[name] += 1
            stats['truth_span_only'] += 1
    not_used = [n for n in hap if used_span[n] != 1]
    if not_used:
        raise CallError('%s %s: hap32 paths not threaded exactly once through the span: %s'
                        % (rid, graph, [(n, used_span[n]) for n in not_used][:5]))
    for line in wl:
        st = STEP_RE.findall(line.rstrip('\n').split('\t')[6])
        for (o1, a), (o2, b) in zip(st, st[1:]):
            e = (a, '+' if o1 == '>' else '-', b, '+' if o2 == '>' else '-')
            er = (b, '-' if o2 == '>' else '+', a, '-' if o1 == '>' else '+')
            if e not in edges and er not in edges:
                edges.add(e)
                stats['junction_edges'] += 1
    gfa = os.path.join(od, 'hybrid.gfa')
    with open(gfa + '.tmp', 'w') as f:
        f.write('H\tVN:Z:1.0\tRS:Z:CHM13%s\n' % (' GRCh38' if grch38 else ''))
        for n, s in all_nodes.items():
            f.write('S\t%s\t%s\n' % (n, s))
        for a, ao, b, bo in sorted(edges):
            f.write('L\t%s\t%s\t%s\t%s\t0M\n' % (a, ao, b, bo))
        f.writelines(wl)
    os.replace(gfa + '.tmp', gfa)
    gbz = os.path.join(od, 'call.gbz')
    run([VG_PINNED, 'gbwt', '-G', gfa, '--gbz-format', '-g', gbz], stderr=os.path.join(od, 'gbwt.log'))
    conv = run([VG_PINNED, 'convert', '-f', gbz]).stdout
    got = {l.split('\t')[1]: l.split('\t')[2] for l in conv.splitlines() if l.startswith('S\t')}
    if got != dict(all_nodes):
        raise CallError('%s %s: hybrid GBZ nodes differ from hybrid.gfa' % (rid, graph))
    info = {'region_id': rid, 'graph': graph, 'flank': flank, 'contig': contig, 'window0': [lo, hi],
            'offset0': offset0, 'chm13_path_len': len(chm_seq), 'flank_nodes': len(flank_nodes),
            'span_nodes': len(new_nodes), 'span_segments': len(snodes), 'ids_interleaved': ids_interleaved,
            'id_tag': tag, 'id_offset': id_offset, 'native_ids': graph == 'mc' and id_offset == 0,
            'paths': len(wl), 'pieces': dict(stats), 'grch38': grch38, 'source_gfa': span_gfa(reg, graph)}
    with open(done, 'w') as f:
        json.dump(info, f, indent=1)
    return info


def hybrid_orig_gaf(rid, flank):
    """The genome-wide alignments of the window (prep_gwsub's reads.gaf) that lie wholly on the
    hybrid mc graph's nodes -- the hybrid mc graph keeps every genome-wide node ID, so they need no
    translation. Alignments through span nodes mc.gfa lacks (nodes only non-spanning fragments
    use) are dropped and counted. -> <hybrid dir>/reads.orig.gaf"""
    od = hybrid_dir('mc', rid, flank)
    out = os.path.join(od, 'reads.orig.gaf')
    if os.path.exists(out):
        return out
    prep_gwsub(rid, flank)
    nodes, _, _ = parse_local_gfa(os.path.join(od, 'hybrid.gfa'))
    kept = dropped = 0
    with open(os.path.join(graph_dir('gwsub%dk' % (flank // 1000), rid), 'reads.gaf')) as f, open(out + '.tmp', 'w') as fo:
        for line in f:
            x = line.split('\t', 7)
            if len(x) > 5 and x[5] != '*' and any(n not in nodes for _, n in STEP_RE.findall(x[5])):
                dropped += 1
                continue
            kept += 1
            fo.write(line)
    os.replace(out + '.tmp', out)
    with open(os.path.join(od, 'reads.orig.json'), 'w') as f:
        json.dump({'kept': kept, 'dropped_off_graph': dropped}, f)
    return out


def fetch_window_reads(rid, flank):
    """The reads the genome-wide mapping placed in the hybrid window, rebuilt with sequence and base
    qualities, mates kept together: remap_local.fetch with flank = F (the anchor-to-anchor subgraph
    plus the CHM13 window span +- F, context 100). -> <stage3>/hybrid<F>k/reads/<id>/"""
    import remap_local as rl
    reg = rl.load_region(rid)
    wd = os.path.join(STAGE3, 'hybrid%dk' % (flank // 1000), 'reads', rid)
    if os.path.exists(os.path.join(wd, 'reads.json')):
        with open(os.path.join(wd, 'reads.json')) as f:
            return json.load(f)
    return rl.fetch(reg, wd, flank=flank)


def map_hybrid(rid, graph, flank, threads=2, stage2=None, force=False, tag=''):
    """giraffe the window's reads to the hybrid graph as Stage 2 mapped its local graphs: vg autoindex
    -n -w sr-giraffe; giraffe paired with Stage 2's pooled fragment-length distribution, single-end
    for reads without their mate; GAF with --named-coordinates (segment names are the node IDs here,
    since nothing is chopped by the indexer: checked)."""
    od = hybrid_dir(graph, rid, flank, tag)
    out = os.path.join(od, 'reads.gaf')
    # kept bgzipped (vg call --gaf-reads reads it; byte-identical VCF, checked at L015415): the
    # uncompressed GAFs of a 149-region run need more disk than the machine has
    if os.path.exists(out + '.gz') and not force:
        return out + '.gz'
    if os.path.exists(out) and not force:
        run(['bgzip', '-f', out])
        return out + '.gz'
    build_hybrid(rid, graph, flank, tag=tag)
    rd = os.path.join(STAGE3, 'hybrid%dk' % (flank // 1000), 'reads', rid)
    fetch_window_reads(rid, flank)
    with open(os.path.join(stage2 or STAGE2_WORK, 'fraglen.json')) as f:
        frag = json.load(f)
    prefix = os.path.join(od, 'idx')
    tmp = os.path.join(od, 'tmp')
    os.makedirs(tmp, exist_ok=True)
    run([VG_PINNED, 'autoindex', '-n', '-w', 'sr-giraffe', '-g', os.path.join(od, 'hybrid.gfa'), '-p', prefix,
         '-t', threads, '-T', tmp, '-M', '4G'], stdout=os.path.join(od, 'autoindex.out'),
        stderr=os.path.join(od, 'autoindex.log'))
    mn = prefix + '.shortread.withzip.min'
    if not os.path.exists(mn):
        mn = prefix + '.min'
    base = [VG_PINNED, 'giraffe', '-Z', prefix + '.giraffe.gbz', '-d', prefix + '.dist', '-m', mn]
    if os.path.exists(prefix + '.shortread.zipcodes'):
        base += ['-z', prefix + '.shortread.zipcodes']
    base += ['-o', 'gaf', '--named-coordinates', '-t', threads]
    r1, r2, se = (os.path.join(rd, x) for x in ('reads_1.fq', 'reads_2.fq', 'reads_se.fq'))
    with open(out + '.tmp', 'w') as fo:
        if os.path.getsize(r1):
            run(base + ['-f', r1, '-f', r2, '--fragment-mean', '%.1f' % frag['mean'],
                        '--fragment-stdev', '%.1f' % frag['stdev']], stdout=fo,
                stderr=os.path.join(od, 'giraffe.paired.log'))
        if os.path.getsize(se):
            run(base + ['-f', se], stdout=fo, stderr=os.path.join(od, 'giraffe.single.log'))
    nodes, _, _ = parse_local_gfa(os.path.join(od, 'hybrid.gfa'))
    with open(out + '.tmp') as f:
        for line in f:
            if line.startswith('@'):
                continue
            x = line.split('\t', 7)
            if x[5] != '*' and any(n not in nodes for _, n in STEP_RE.findall(x[5])):
                raise CallError('%s %s: GAF names a node the hybrid graph does not have' % (rid, graph))
    os.replace(out + '.tmp', out)
    shutil.rmtree(tmp, ignore_errors=True)
    for x in os.listdir(od):
        if x.startswith('idx.') and not x.endswith('.log'):
            os.remove(os.path.join(od, x))
    run(['bgzip', '-f', out])
    return out + '.gz'


# ------------------------------------------------------------------ genome-wide alignments, projected (mc@orig)

def project_orig(rid, stage2=None):
    """The genome-wide alignments of the Stage 2 reads, moved into the local mc graph's node space.
    Inside the anchors mc.gfa keeps the hap32 node IDs, so steps there are unchanged; a run of steps
    outside the anchors (at the start or end of the read's path) must be a contiguous piece of the
    CHM13 walk next to the span and is replaced by the flank node, with the path offsets shifted.
    Reads whose outside part leaves the CHM13 walk, or runs past the 1 kb flank, are dropped."""
    stage2 = stage2 or STAGE2_WORK
    reg = load_region(rid)
    gd = graph_dir('mc', rid)
    dst = os.path.join(gd, 'reads.orig.gaf')
    if os.path.exists(dst) and os.path.exists(os.path.join(gd, 'reads.orig.json')):
        with open(os.path.join(gd, 'reads.orig.json')) as f:
            return json.load(f)
    s2 = os.path.join(stage2, rid)
    with open(os.path.join(s2, 'mc', 'map.json')) as f:
        m = json.load(f)
    lf, rf = m['left_flank_node'], m['right_flank_node']
    nodes, _, _ = parse_local_gfa(os.path.join(gd, 'call.gfa'))
    inside = set(nodes) - {lf, rf}
    L, R = len(nodes[lf]), len(nodes[rf])
    # CHM13 coordinates of genome-wide nodes around the span (0-based start, end, orientation)
    chm = {}
    for q in ('q_window.gfa', 'q_nodes.gfa'):
        p = os.path.join(s2, q)
        if not os.path.exists(p):
            continue
        qn = {}
        with open(p) as f:
            for line in f:
                if line.startswith('S\t'):
                    x = line.rstrip('\n').split('\t')
                    qn[x[1]] = len(x[2])
        with open(p) as f:
            for line in f:
                if line.startswith('W\tCHM13\t'):
                    x = line.rstrip('\n').split('\t')
                    off = int(x[4])
                    for o, n in STEP_RE.findall(x[6]):
                        chm.setdefault(n, (off, off + qn[n], o))
                        off += qn[n]
    a0, b0 = reg['span_start'] - 1, reg['span_end']        # span, 0-based half-open
    # reads.tsv: local name, mate, original alignment; fastq for qualities
    quals = {}
    for fq, mate in (('reads_1.fq', 1), ('reads_2.fq', 2), ('reads_se.fq', 0)):
        with open(os.path.join(s2, fq)) as f:
            while True:
                h = f.readline()
                if not h:
                    break
                s = f.readline().strip()
                f.readline()
                q = f.readline().strip()
                quals[(h[1:].strip(), mate)] = (s, q)
    st = collections.Counter()
    out = []
    with open(os.path.join(s2, 'reads.tsv')) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for line in f:
            r = dict(zip(hdr, line.rstrip('\n').split('\t')))
            st['reads'] += 1
            ln, mate = r['local'], int(r['mate'])
            key = (ln, mate if ln.startswith('p') else 0)
            seq, qual = quals[key]
            steps = [(o, n) for o, n in STEP_RE.findall(r['orig_path'])]
            ps, pe, plen = int(r['orig_ps']), int(r['orig_pe']), int(r['orig_plen'])
            ins = [n in inside for _, n in steps]
            if all(ins):
                new, nps, npe, nplen = steps, ps, pe, plen
            else:
                # leading and trailing outside runs; an outside step between inside ones is not allowed
                i = 0
                while i < len(steps) and not ins[i]:
                    i += 1
                j = len(steps)
                while j > i and not ins[j - 1]:
                    j -= 1
                if any(not x for x in ins[i:j]):
                    st['drop_outside_inside_span'] += 1
                    continue
                lead, mid, trail = steps[:i], steps[i:j], steps[j:]
                ok = True
                new = list(mid)
                nps, npe = ps, pe
                shift = 0
                for part, where in ((lead, 'lead'), (trail, 'trail')):
                    if not part:
                        continue
                    if any(n not in chm for _, n in part):
                        ok = False
                        st['drop_off_chm13_walk'] += 1
                        break
                    fwd = all((o == '>') == (chm[n][2] == '>') for o, n in part)
                    rev = all((o == '<') == (chm[n][2] == '>') for o, n in part)
                    ivs = sorted((chm[n][0], chm[n][1]) for _, n in part)
                    lo, hi = ivs[0][0], ivs[-1][1]
                    if sum(e - s for s, e in ivs) != hi - lo or not (fwd or rev):
                        ok = False
                        st['drop_not_colinear'] += 1
                        break
                    plen_part = hi - lo
                    # which flank: left flank ends at a0, right flank starts at b0
                    if hi == a0 and plen_part <= L:
                        flank, flen = lf, L
                    elif lo == b0 and plen_part <= R:
                        flank, flen = rf, R
                    elif not mid and lo >= a0 - L and hi <= a0:
                        flank, flen = lf, L
                    elif not mid and lo >= b0 and hi <= b0 + R:
                        flank, flen = rf, R
                    else:
                        ok = False
                        st['drop_beyond_flank'] += 1
                        break
                    # orientation of the flank step in the path
                    o = '>' if fwd else '<'
                    # bases of the flank node before the part, in path orientation
                    if flank == lf:
                        before = (lo - (a0 - L)) if o == '>' else ((a0) - hi)
                    else:
                        before = (lo - b0) if o == '>' else ((b0 + R) - hi)
                    after = flen - plen_part - before
                    if where == 'lead':
                        new = [(o, flank)] + new
                        shift += before
                    else:
                        new = new + [(o, flank)]
                    if not mid:
                        # the whole read is in one flank: path = that node only
                        new = [(o, flank)]
                        shift = before
                        break
                    _ = after
                if not ok:
                    continue
                nps, npe = ps + shift, pe + shift
                nplen = sum(len(nodes[n]) for _, n in new)
                st['projected_with_flank'] += 1
            # the aligned read must still be spelled (sequence from path + cs)
            pseq = ''.join(nodes[n] if o == '>' else rc(nodes[n]) for o, n in new)
            if nplen != len(pseq) or not (0 <= nps <= npe <= nplen):
                st['drop_bad_offsets'] += 1
                continue
            rebuilt = read_from_cs(pseq[nps:npe], r['orig_cs'])
            if rebuilt != seq:
                st['drop_cs_mismatch'] += 1
                continue
            tags = ['AS:i:0', 'bq:Z:' + qual, 'cs:Z:' + r['orig_cs']]
            if mate == 1 and ln.startswith('p'):
                tags.append('fn:Z:' + ln)
            elif mate == 2 and ln.startswith('p'):
                tags.append('fp:Z:' + ln)
            out.append('\t'.join([ln, str(len(seq)), '0', str(len(seq)), '+',
                                  ''.join(o + n for o, n in new), str(nplen), str(nps), str(npe),
                                  str(len(seq)), str(len(seq)), r['orig_mapq']] + tags) + '\n')
            st['kept'] += 1
    with open(dst + '.tmp', 'w') as f:
        f.writelines(out)
    os.replace(dst + '.tmp', dst)
    with open(os.path.join(gd, 'reads.orig.json'), 'w') as f:
        json.dump(dict(st), f, indent=1)
    return dict(st)


CS_RE = re.compile(r'(:\d+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+|=[A-Za-z]+)')


def read_from_cs(pseq, cs):
    out, i = [], 0
    for o in CS_RE.findall(cs):
        c = o[0]
        if c == ':':
            n = int(o[1:])
            out.append(pseq[i:i + n])
            i += n
        elif c == '=':
            out.append(o[1:].upper())
            i += len(o) - 1
        elif c == '*':
            out.append(o[2].upper())
            i += 1
        elif c == '+':
            out.append(o[1:].upper())
        elif c == '-':
            i += len(o) - 1
    if i != len(pseq):
        return None
    return ''.join(out)


# ------------------------------------------------------------------ 2. call, shift, check

def calls_path(graph, rid, arm=None):
    if arm:
        return os.path.join(STAGE3, 'diag', '%s@%s' % (graph, arm), rid + '.vcf.gz')
    if graph in NULL_GRAPHS:
        return os.path.join(STAGE3, 'diag', '%s@%s' % (graph, FINAL_ARM), rid + '.vcf.gz')
    return os.path.join(STAGE3, 'calls', graph, rid + '.vcf.gz')


def vg_call_command(reg, gbz, gaf, mosaic, threads, extra=()):
    ploidy = ploidy_of(reg)
    cmd = [VG_PINNED, 'call', gbz, '-p', 'CHM13#0#%s' % reg['contig'], '-d', str(ploidy),
           '-t', str(threads)] + PRODUCTION_FLAGS + ['--mosaic-out', mosaic] + list(extra) + \
          os.environ.get('VG_CALL_EXTRA', '').split() + ['--gaf-reads', gaf]
    return cmd


# The setup `call` uses when no arm is named, chosen by the gate (see the module docstring):
# the genome-wide flank of 200 kb on each side around the graph's span, reads re-mapped to it.
FINAL_ARM = 'hybrid200k'


def call(rid, graph, arm=None, threads=1, stage2=None, force=False, compact=False):
    """Call one region on one graph. arm None writes the Stage 3 call (FINAL_ARM) to
    <stage3>/calls/<graph>/<id>.vcf.gz; a named arm writes <stage3>/diag/<graph>@<arm>/<id>.vcf.gz:
      stage2            the Stage 2 local graph (1 kb CHM13 flank) with the Stage 2 re-mapped reads
      stage2nodepth     the same with --depth-term 0
      stage2orig        the Stage 2 local mc graph with the genome-wide alignments projected into it
      gwsub<F>k         the genome-wide hap32 graph over span +- F kb with the genome-wide alignments
      hybrid<F>k        the genome-wide flank (+- F kb) around the graph's span, reads re-mapped to it
      hybrid<F>kids     the same with the span's node IDs re-laid (build_hybrid tag 'ids'), re-mapped
      hybrid<F>kmcids   the same with the span's nodes on mc's native span IDs (tag 'mcids'; not mc)
      ...nodepth        any of these with --depth-term 0
    hybrid50k, hybrid200kids and hybrid200kmcids are the REPLICATE arms: neutral changes (window
    size; node numbering) of the Stage 3 call, run for every graph to measure each graph's own noise
    (replicate_noise.py). compact (arms only): after the call, gzip the arm's reads.gaf and delete
    its hybrid.gfa (both regenerable; call.gbz and build.json are kept)."""
    reg = load_region(rid)
    out = calls_path(graph, rid, arm)
    info_path = out.replace('.vcf.gz', '.json')
    if os.path.exists(out) and os.path.exists(info_path) and not force:
        with open(info_path) as f:
            return json.load(f)
    setup = arm or FINAL_ARM
    extra = []
    if setup.endswith('nodepth'):
        extra = ['--depth-term', '0']
        setup = setup[:-len('nodepth')]
    if setup.startswith('hybrid'):
        F, kind = parse_hybrid_setup(setup)
        tag = kind if kind in ID_TAGS else ''
        pj = build_hybrid(rid, graph, F, tag=tag)
        gd = hybrid_dir(graph, rid, F, tag)
        if kind == 'orig':
            if graph != 'mc':
                raise CallError('arm %s exists for mc only (genome-wide node IDs)' % arm)
            gaf = hybrid_orig_gaf(rid, F)
        else:
            gaf = map_hybrid(rid, graph, F, threads=max(1, threads), stage2=stage2, tag=tag)
    elif setup.startswith('gwsub'):
        if graph != 'mc':
            raise CallError('arm %s exists for mc only' % arm)
        F = int(re.match(r'gwsub(\d+)k$', setup).group(1)) * 1000
        pj = prep_gwsub(rid, F)
        gd = graph_dir('gwsub%dk' % (F // 1000), rid)
        gaf = os.path.join(gd, 'reads.gaf')
    elif setup in ('stage2', 'stage2orig'):
        pj = prep(rid, graph, stage2)
        gd = graph_dir(graph, rid)
        gaf = os.path.join(gd, 'reads.gaf')
        if setup == 'stage2orig':
            if graph != 'mc':
                raise CallError('arm stage2orig exists for mc only (its node IDs are the genome-wide ones)')
            project_orig(rid, stage2)
            gaf = os.path.join(gd, 'reads.orig.gaf')
    else:
        raise CallError('unknown arm %s' % arm)
    od = os.path.dirname(out)
    os.makedirs(od, exist_ok=True)
    raw = out.replace('.vcf.gz', '.raw.vcf')
    errlog = out.replace('.vcf.gz', '.log')
    mosaic = out.replace('.vcf.gz', '.mosaic.tsv')
    cmd = vg_call_command(reg, os.path.join(gd, 'call.gbz'), gaf, mosaic, threads, extra)
    t0 = time.time()
    run(cmd, stdout=raw, stderr=errlog)
    secs = time.time() - t0
    chk = finish_vcf(raw, out, reg, pj)
    with open(errlog) as f:
        err = f.read()
    loaded = re.search(r'Loaded (\d+) reads \((\d+) filtered out\)', err)
    enum = re.search(r'Enumerating alleles from the (\d+) GBZ panel haplotypes', err)
    link = re.search(r'Linkage: (\d+) panel haplotypes over (\d+) GBWT sequences', err)
    info = {'region_id': rid, 'graph': graph, 'arm': arm, 'setup': arm or FINAL_ARM, 'vcf': out, 'command': cmd,
            'seconds': round(secs, 1), 'shift': chk, 'offset0': pj['offset0'],
            'reads_loaded': int(loaded.group(1)) if loaded else None,
            'reads_filtered': int(loaded.group(2)) if loaded else None,
            'enumeration_haplotypes': int(enum.group(1)) if enum else None,
            'linkage_haplotypes': int(link.group(1)) if link else None,
            'production': {'source': PRODUCTION_SOURCE, 'flags': PRODUCTION_FLAGS,
                           'reads_replaced': PRODUCTION_READS, 'threads': PRODUCTION_THREADS}}
    if compact and arm and setup.startswith('hybrid') and not setup.endswith('orig'):
        if os.path.exists(gaf) and not gaf.endswith('.gz'):
            run(['gzip', '-f', '-1', gaf])
            info['compacted'] = [os.path.basename(gaf) + '.gz']
        hg = os.path.join(gd, 'hybrid.gfa')
        if os.path.exists(hg):
            os.remove(hg)
            info['compacted'] = info.get('compacted', []) + ['hybrid.gfa removed']
    with open(info_path, 'w') as f:
        json.dump(info, f, indent=1)
    return info


def finish_vcf(raw, out, reg, pj):
    """Shift, check, bgzip and index. vg's own record order is KEPT (bgzip, not bcftools sort): two
    records can share POS and REF length -- a nested parent and its child both written at one
    position -- and a haplotype builder that applies records in file order then resolves their
    overlap by that order, as it does on the genome-wide VCF, which is vg's order too. bcftools sort
    reorders such ties (by ALT), which alone changed L012184's called haplotypes by 25 edits."""
    shifted = out[:-3]
    chk = shift_vcf(raw, shifted, reg, pj)
    last, unsorted = -1, 0
    with open(shifted) as f:
        for line in f:
            if not line.startswith('#'):
                pos = int(line.split('\t', 2)[1])
                unsorted += pos < last
                last = pos
    if unsorted:
        chk['resorted'] = unsorted
        run([config.BCFTOOLS, 'sort', '-Oz', '-o', out, shifted])
    else:
        run([config.BGZIP, '-f', shifted])
    run([config.BCFTOOLS, 'index', '-f', '-t', out])
    if os.path.exists(shifted):
        os.remove(shifted)
    return chk


def refinish(vcf_gz):
    """Rebuild <id>.vcf.gz from the kept <id>.raw.vcf with the current finish_vcf (no vg run)."""
    info_path = vcf_gz.replace('.vcf.gz', '.json')
    with open(info_path) as f:
        info = json.load(f)
    reg = load_region(info['region_id'])
    setup = (info.get('setup') or info.get('arm') or 'stage2').replace('nodepth', '')
    if info.get('arm') in (None, 'orig', 'nodepth') and not info.get('setup'):
        setup = 'stage2'
    if setup.startswith('hybrid'):
        F, kind = parse_hybrid_setup(setup)
        pj = build_hybrid(info['region_id'], info['graph'], F, tag=kind if kind in ID_TAGS else '')
    elif setup.startswith('gwsub'):
        pj = prep_gwsub(info['region_id'], int(re.match(r'gwsub(\d+)k', setup).group(1)) * 1000)
    else:
        pj = prep(info['region_id'], info['graph'])
    info['shift'] = finish_vcf(vcf_gz.replace('.vcf.gz', '.raw.vcf'), vcf_gz, reg, pj)
    with open(info_path, 'w') as f:
        json.dump(info, f, indent=1)
    return info['shift']


def shift_vcf(raw, dst, reg, pj):
    """Local path coordinates -> CHM13: CHROM = contig, POS += offset0; REF checked against CHM13."""
    contig = reg['contig']
    off = pj['offset0']
    clen = contig_length(contig)
    lo, hi = off, off + pj['chm13_path_len']
    ref = fetch_ref(contig, lo, hi)
    n = bad = 0
    chroms = collections.Counter()
    with open(raw) as f, open(dst, 'w') as fo:
        for line in f:
            if line.startswith('##contig='):
                fo.write('##contig=<ID=%s,length=%d>\n' % (contig, clen))
                continue
            if line.startswith('#'):
                fo.write(line)
                continue
            x = line.rstrip('\n').split('\t')
            chroms[x[0]] += 1
            pos = int(x[1]) + off
            r = x[3].upper()
            i = pos - 1 - lo
            if i < 0 or ref[i:i + len(r)] != r:
                bad += 1
            x[0] = contig
            x[1] = str(pos)
            if 'END=' in x[7]:
                x[7] = re.sub(r'(^|;)END=(\d+)', lambda mm: '%sEND=%d' % (mm.group(1), int(mm.group(2)) + off), x[7])
            fo.write('\t'.join(x) + '\n')
            n += 1
    if bad:
        raise CallError('%s: %d of %d REF alleles disagree with CHM13 after the shift' % (reg['region_id'], bad, n))
    return {'records': n, 'ref_checked': n, 'ref_mismatch': 0, 'raw_chrom': dict(chroms)}


# ------------------------------------------------------------------ 3. haplotypes over the span (the scoring rule)

_ED = None


def edit_distance(a, b):
    """Exact unit-cost edit distance (tools/fastedit.c via region.py's loader)."""
    global _ED
    if _ED is None:
        import region as _region
        _ED = _region.edit_distance
    if a == b:
        return 0
    return _ED(a, b)


def parse_vcf_span(vcf, contig, a1, b1, sample='HG002'):
    """Records whose REF span overlaps [a1, b1] (1-based), with GT, PS, FILTER."""
    p = run([config.BCFTOOLS, 'view', '-r', '%s:%d-%d' % (contig, a1, b1), vcf])
    recs = []
    si = None
    for line in p.stdout.splitlines():
        if line.startswith('##'):
            continue
        x = line.split('\t')
        if line.startswith('#'):
            si = x.index(sample) if sample in x else 9
            continue
        fmt = x[8].split(':')
        val = x[si].split(':')
        d = dict(zip(fmt, val))
        gt = d.get('GT', '.')
        al = [None if g == '.' else int(g) for g in re.split(r'[|/]', gt)]
        ref = x[3].upper()
        bnd, interior = _scorer().parse_vg_snarl(x[2], x[7])
        recs.append({'pos': int(x[1]), 'end': int(x[1]) + len(ref) - 1, 'id': x[2], 'ref': ref,
                     'alts': [a.upper() for a in x[4].split(',')], 'filter': x[6],
                     'gt': al, 'phased': '|' in gt, 'ps': d.get('PS'), 'raw_gt': gt,
                     'idx': len(recs), 'bnd': bnd, 'interior': interior})
    return recs


def _scorer():
    import score_haplotypes
    return score_haplotypes


def _trim(pos, ref, alt):
    r, a = ref, alt
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    k = 0
    while k < len(r) and k < len(a) and r[k] == a[k]:
        k += 1
    return pos + k, pos + len(r), a[k:]         # [s, e) 1-based half-open replaced by a


def _allele_edit(rec, a, a1, b1, refseq, flags):
    alt = rec['alts'][a - 1] if a - 1 < len(rec['alts']) else None
    if alt is None or alt in ('*', '.') or alt.startswith('<') or '[' in alt or ']' in alt \
            or not re.fullmatch(r'[ACGTN]+', alt):
        flags['symbolic_as_ref'] += 1
        return None
    pos, ref = rec['pos'], rec['ref']
    end = pos + len(ref) - 1
    lo, hi = max(pos, a1), min(end, b1)
    if lo <= hi and ref[lo - pos:hi - pos + 1] != refseq[lo - a1:hi - a1 + 1]:
        flags['ref_mismatch_as_ref'] += 1
        return None
    if pos >= a1 and end <= b1:
        return (pos, end + 1, alt)
    flags['crossing_span_edge'] += 1
    s, e, t = _trim(pos, ref, alt)
    if s == e:
        if a1 < s <= b1:
            return (s, e, t)
        flags['edge_insertion_dropped'] += 1
        return None
    if e - 1 < a1 or s > b1:
        flags['outside_after_trim'] += 1
        return None
    if s < a1:
        n = a1 - s
        t = t[min(n, len(t)):]
        s = a1
    if e - 1 > b1:
        k = e - 1 - b1
        t = t[:len(t) - min(k, len(t))]
        e = b1 + 1
    return (s, e, t)


def _apply(refseq, a1, edits):
    out, cur = [], a1
    for s, e, t in sorted(edits):
        out.append(refseq[cur - a1:s - a1])
        out.append(t)
        cur = e
    out.append(refseq[cur - a1:])
    return ''.join(out)


def prepare_haps(recs, a1, b1, refseq, ploidy):
    """Items (raw_start, raw_end, (edit_slot0, edit_slot1), unit) and flags, in the order they are
    applied: (POS, longest REF, file order), except that a vg record never precedes a record it is
    nested in (score_haplotypes.nesting_order; vg's prefix trimming can put a child's POS left of
    its parent's, and the child then won the overlap and the parent's allele was dropped)."""
    flags = collections.Counter()
    items, units = [], collections.OrderedDict()
    ordered, _, nstats = _scorer().nesting_order(recs)
    flags.update(nstats)
    for r in ordered:
        flags['records'] += 1
        gt = list(r['gt'])
        if len(gt) == 1:
            flags['haploid_gt'] += 1
            gt = gt * 2
        gt = gt[:2]
        if any(g is None for g in gt):
            flags['missing_allele_as_ref'] += sum(1 for g in gt if g is None)
            gt = [0 if g is None else g for g in gt]
        if not any(gt):
            continue
        if r['filter'] not in ('PASS', '.'):
            flags['filtered_as_ref'] += 1
            continue
        ed = [None if g == 0 else _allele_edit(r, g, a1, b1, refseq, flags) for g in gt]
        if ed[0] is None and ed[1] is None:
            continue
        unit = None
        if ed[0] != ed[1]:
            key = ('PS', r['ps']) if r['phased'] else ('rec', r['idx'])
            if not r['phased']:
                flags['unphased_het'] += 1
            unit = units.setdefault(key, len(units))
        items.append((r['pos'], r['end'], tuple(ed), unit))
    return items, len(units), flags


def build_pair(items, flips, refseq, a1):
    """The raw overlap rule in the items' order: an allele is skipped when its REF span shares a
    base with an allele already applied to that slot (checked against every applied span, so the
    rule does not depend on the items being sorted by POS)."""
    edits = ([], [])
    spans = ([], [])
    skipped = 0
    for s, e, (x0, x1), u in items:
        if u is not None and flips[u]:
            x0, x1 = x1, x0
        for h, xx in ((0, x0), (1, x1)):
            if xx is None:
                continue
            if any(s <= pe and ps_ <= e for ps_, pe in spans[h]):
                skipped += 1
                continue
            spans[h].append((s, e))
            edits[h].append(xx)
    return (_apply(refseq, a1, edits[0]), _apply(refseq, a1, edits[1])), skipped


def pair_ed(c, t):
    if len(t) == 1:
        return edit_distance(c[0], t[0]) + edit_distance(c[1], t[0])
    return min(edit_distance(c[0], t[0]) + edit_distance(c[1], t[1]),
               edit_distance(c[0], t[1]) + edit_distance(c[1], t[0]))


def truth_haps(reg):
    t = dict(read_fasta(os.path.join(reg['dir'], 'truth.fa')))
    if ploidy_of(reg) == 1:
        return (t['HG002#1'],)
    return (t['HG002#1'], t['HG002#2'])


def span_haplotypes(vcf, reg, enum_cap=12, target=None):
    """The called pair over the span, phase units oriented to minimise ED to `target` (default:
    the truth). Returns dict with seqs, ED, flags."""
    a1, b1 = reg['span_start'], reg['span_end']
    refseq = fetch_ref(reg['contig'], a1 - 1, b1)
    recs = parse_vcf_span(vcf, reg['contig'], a1, b1)
    ploidy = ploidy_of(reg)
    items, k, flags = prepare_haps(recs, a1, b1, refseq, ploidy)
    tgt = target or truth_haps(reg)
    cache = {}

    def score(fl):
        key = tuple(fl)
        if key not in cache:
            seqs, sk = build_pair(items, fl, refseq, a1)
            cache[key] = (pair_ed(seqs, tgt), seqs, sk)
        return cache[key]

    flips = [False] * k
    greedy = False
    if k >= 2 and k - 1 <= enum_cap:
        best = None
        for combo in itertools.product([False, True], repeat=k - 1):
            fl = [False] + list(combo)
            sc = score(fl)
            if best is None or sc[0] < best[0][0]:
                best = (sc, fl)
        flips = best[1]
    elif k - 1 > enum_cap:
        greedy = True
        cur = score(flips)[0]
        changed = True
        while changed:
            changed = False
            for u in range(1, k):
                fl = list(flips)
                fl[u] = not fl[u]
                v = score(fl)[0]
                if v < cur:
                    flips, cur, changed = fl, v, True
    ed, seqs, skipped = score(flips)
    written = score([False] * k)[1]
    if greedy:
        flags['phase_greedy'] += 1
    flags['overlap_skipped'] += skipped
    ed_ref = pair_ed((refseq, refseq), tgt)
    return {'seqs': seqs, 'written_seqs': written, 'ed': ed, 'ed_ref': ed_ref, 'units': k,
            'records': len(recs), 'flags': dict(flags), 'ploidy': ploidy,
            'truth_len': sum(len(t) for t in tgt) * (2 if len(tgt) == 1 else 1)}


# ------------------------------------------------------------------ 4. truvari over the span

def truvari_span(vcf, reg, od):
    """truvari bench (the genome-wide pipeline's parameters) on the records overlapping the span."""
    contig, a1, b1 = reg['contig'], reg['span_start'], reg['span_end']
    os.makedirs(od, exist_ok=True)
    ref = config.data_paths(contig)['ref_fa']
    reg_s = '%s:%d-%d' % (contig, a1, b1)
    comp_raw = os.path.join(od, 'comp.raw.vcf.gz')
    run([config.BCFTOOLS, 'view', '-r', reg_s, '-Oz', '-o', comp_raw, vcf])
    samp = os.path.join(od, 'sample.txt')
    with open(samp, 'w') as f:
        f.write('HG002\n')
    comp_rh = os.path.join(od, 'comp.rh.vcf.gz')
    run([config.BCFTOOLS, 'reheader', '-s', samp, '-o', comp_rh, comp_raw])
    tmp = os.path.join(od, 'comp.split.vcf.gz')
    run([config.BCFTOOLS, 'norm', '-m-any', '-f', ref, '-Oz', '-o', tmp, comp_rh])
    comp = os.path.join(od, 'comp.norm.vcf.gz')
    run([config.BCFTOOLS, 'sort', '-Oz', '-o', comp, tmp])
    run([config.BCFTOOLS, 'index', '-f', '-t', comp])
    tn = os.path.join(config.EVAL_DIR, 'work/wgs-mm095/score', contig + '.truth.norm.vcf.gz')
    base = os.path.join(od, 'base.vcf.gz')
    run([config.BCFTOOLS, 'view', '-r', reg_s, '-Oz', '-o', base, tn])
    run([config.BCFTOOLS, 'index', '-f', '-t', base])
    bed = os.path.join(config.EVAL_DIR, 'work/wgs-mm095', contig, 'truth.%s.stvar.bed' % contig)
    tv = os.path.join(od, 'tv')
    if os.path.exists(tv):
        shutil.rmtree(tv)
    run([TRUVARI, 'bench', '-b', base, '-c', comp, '-f', ref, '-o', tv, '--includebed', bed,
         '--sizemin', '50', '--sizefilt', '50', '--bSample', 'HG002', '--cSample', 'HG002',
         '--pick', 'ac'])
    with open(os.path.join(tv, 'summary.json')) as f:
        s = json.load(f)
    return {k: s.get(k) for k in ('TP-base', 'TP-comp', 'FP', 'FN', 'precision', 'recall', 'f1')}


def gw_labels(reg):
    """Genome-wide truvari labels (work/wgs-mm095/score/<contig>.truvari) overlapping the span."""
    d = config.data_paths(reg['contig'])['vg_truvari']
    out = {}
    for k, fn in (('TP-base', 'tp-base'), ('FN', 'fn'), ('TP-comp', 'tp-comp'), ('FP', 'fp')):
        p = run([config.BCFTOOLS, 'view', '-H', '-r', '%s:%d-%d' % (reg['contig'], reg['span_start'],
                 reg['span_end']), os.path.join(d, fn + '.vcf.gz')])
        out[k] = sum(1 for l in p.stdout.splitlines() if l.strip())
    return out


# ------------------------------------------------------------------ 6. the gate

def gw_vcf(reg):
    return config.data_paths(reg['contig'])['vg_vcf']


def snarl_records(vcf, reg):
    """{(snarl ID, SB block): (called allele sequences in GT order, DR, DP)} over the span. The ID
    column is vg's snarl name, its two boundary nodes, which is a graph identity (not a POS) for
    every graph that keeps the genome-wide node IDs: mc in every setup here."""
    p = run([config.BCFTOOLS, 'view', '-H', '-r', '%s:%d-%d' % (reg['contig'], reg['span_start'],
             reg['span_end']), vcf])
    out = {}
    for line in p.stdout.splitlines():
        x = line.split('\t')
        f = dict(zip(x[8].split(':'), x[9].split(':')))
        sb = re.search(r'(?:^|;)SB=([^;]+)', x[7])
        alle = [x[3]] + x[4].split(',')
        gt = [None if g == '.' else alle[int(g)] for g in re.split(r'[|/]', f.get('GT', '.'))]
        out[(x[2], sb.group(1) if sb else '')] = (tuple(gt), f.get('DR'), f.get('DP'), f.get('GQ'))
    return out


def record_diff(vcf_a, vcf_b, reg):
    """How two mc-family call sets differ record by record: same genotype, same genotype in the
    other phase orientation, different genotype, or a snarl/block present in one only; and the
    median ratio of DR (observed/expected reads, so the depth-rate estimate) where both have it."""
    A, B = snarl_records(vcf_a, reg), snarl_records(vcf_b, reg)
    c = collections.Counter()
    ratios = []
    for k in set(A) | set(B):
        if k not in A:
            c['only_b'] += 1
            continue
        if k not in B:
            c['only_a'] += 1
            continue
        ga, gb = A[k][0], B[k][0]
        if ga == gb:
            c['same'] += 1
        elif sorted(ga, key=str) == sorted(gb, key=str):
            c['phase_only'] += 1
        else:
            c['genotype'] += 1
        try:
            ratios.append(float(A[k][1]) / float(B[k][1]))
        except (TypeError, ValueError, ZeroDivisionError):
            pass
    ratios.sort()
    c = dict(c)
    c['dr_ratio_median'] = round(ratios[len(ratios) // 2], 3) if ratios else None
    c['dr_ratio_eq1_frac'] = round(sum(1 for r in ratios if abs(r - 1) < 1e-3) / len(ratios), 3) if ratios else None
    return c


GATE_ARMS = ['stage2', 'stage2nodepth', 'stage2orig', 'gwsub10k', 'gwsub50k', 'gwsub200k', 'hybrid200korig',
             'hybrid50k', 'hybrid200k', 'hybrid200knodepth']
GATE_TRUVARI_ARMS = ['stage2', 'hybrid200korig', 'hybrid200k']


def gate_one(rid, arms=None, truvari_arms=None):
    arms = arms or GATE_ARMS
    truvari_arms = GATE_TRUVARI_ARMS if truvari_arms is None else truvari_arms
    reg = load_region(rid)
    gw = gw_vcf(reg)
    row = {'region_id': rid, 'stratum': reg['stratum'], 'contig': reg['contig'],
           'span': '%s:%d-%d' % (reg['contig'], reg['span_start'], reg['span_end']), 'span_bp': reg['span_bp'],
           'ploidy': ploidy_of(reg)}
    hg = span_haplotypes(gw, reg)
    row['gw'] = {'records': hg['records'], 'ed_truth': hg['ed'], 'ed_ref': hg['ed_ref'], 'units': hg['units'],
                 'flags': hg['flags']}
    row['ed_ref'] = hg['ed_ref']
    row['truth_len'] = hg['truth_len']
    tvd = os.path.join(STAGE3, 'gate', rid)
    row['gw']['truvari_span'] = truvari_span(gw, reg, os.path.join(tvd, 'gw'))
    row['gw']['labels'] = gw_labels(reg)
    row['arms'] = {}
    for arm in arms:
        p = calls_path('mc', rid, arm)
        if not os.path.exists(p):
            continue
        ha = span_haplotypes(p, reg)
        same = sorted(ha['seqs']) == sorted(hg['seqs'])
        d = {'records': ha['records'], 'identical': same,
             'ed_vs_gw': 0 if same else pair_ed(ha['seqs'], tuple(hg['seqs'])),
             'ed_truth': ha['ed'], 'units': ha['units'], 'flags': ha['flags'],
             'records_vs_gw': record_diff(p, gw, reg)}
        if arm in truvari_arms:
            d['truvari_span'] = truvari_span(p, reg, os.path.join(tvd, arm))
        with open(p.replace('.vcf.gz', '.json')) as f:
            ci = json.load(f)
        d['reads_loaded'] = ci.get('reads_loaded')
        d['enumeration_haplotypes'] = ci.get('enumeration_haplotypes')
        d['linkage_haplotypes'] = ci.get('linkage_haplotypes')
        row['arms'][arm] = d
    return row


def gate_tables(rows, arms):
    """TSV rows and a markdown table of the gate."""
    tsv = []
    hdr = ['region_id', 'stratum', 'span_bp', 'ed_ref', 'gw_ed_truth', 'gw_FP', 'gw_FN', 'gw_labels_FP', 'gw_labels_FN']
    for a in arms:
        hdr += ['%s_identical' % a, '%s_ed_vs_gw' % a, '%s_ed_truth' % a, '%s_rec_same' % a, '%s_rec_phase_only' % a,
                '%s_rec_genotype' % a, '%s_rec_only_local' % a, '%s_rec_only_gw' % a, '%s_dr_ratio' % a,
                '%s_FP' % a, '%s_FN' % a]
    for r in rows:
        g = r['gw']
        line = [r['region_id'], r['stratum'], r['span_bp'], r['ed_ref'], g['ed_truth'],
                g['truvari_span']['FP'], g['truvari_span']['FN'], g['labels']['FP'], g['labels']['FN']]
        for a in arms:
            d = r['arms'].get(a)
            if not d:
                line += [''] * 11
                continue
            rv = d['records_vs_gw']
            tv = d.get('truvari_span') or {}
            line += [int(d['identical']), d['ed_vs_gw'], d['ed_truth'], rv.get('same', 0), rv.get('phase_only', 0),
                     rv.get('genotype', 0), rv.get('only_a', 0), rv.get('only_b', 0), rv.get('dr_ratio_median'),
                     tv.get('FP', ''), tv.get('FN', '')]
        tsv.append(line)
    return hdr, tsv


def gate_markdown(rows, arms, labels=None):
    """One markdown table: per region the genome-wide call's ED and, per arm, whether the called
    haplotypes are identical to it (= or the edit distance between the two pairs) and the arm's ED."""
    labels = labels or {}
    head = ['region', 'stratum', 'ED_ref', 'genome-wide ED'] + [labels.get(a, a) for a in arms]
    out = ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)]
    for r in rows:
        cells = [r['region_id'], r['stratum'], str(r['ed_ref']), str(r['gw']['ed_truth'])]
        for a in arms:
            d = r['arms'].get(a)
            if not d:
                cells.append('-')
                continue
            rv = d['records_vs_gw']
            kind = ''
            if not d['identical']:
                if not rv.get('genotype') and not rv.get('only_a') and not rv.get('only_b'):
                    kind = ' phase'
                else:
                    kind = ' gt%d' % (rv.get('genotype', 0) + rv.get('only_a', 0) + rv.get('only_b', 0))
            cells.append('=' if d['identical'] else 'd%d%s (ED %d)' % (d['ed_vs_gw'], kind, d['ed_truth']))
        out.append('| ' + ' | '.join(cells) + ' |')
    return '\n'.join(out)


# ------------------------------------------------------------------ 7. scorer validationdef gate_markdown(rows, arms, labels=None):
    """One markdown table: per region the genome-wide call's ED and, per arm, whether the called
    haplotypes are identical to it (= or the edit distance between the two pairs) and the arm's ED."""
    labels = labels or {}
    head = ['region', 'stratum', 'ED_ref', 'genome-wide ED'] + [labels.get(a, a) for a in arms]
    out = ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)]
    for r in rows:
        cells = [r['region_id'], r['stratum'], str(r['ed_ref']), str(r['gw']['ed_truth'])]
        for a in arms:
            d = r['arms'].get(a)
            if not d:
                cells.append('-')
                continue
            rv = d['records_vs_gw']
            kind = ''
            if not d['identical']:
                if not rv.get('genotype') and not rv.get('only_a') and not rv.get('only_b'):
                    kind = ' phase'
                else:
                    kind = ' gt%d' % (rv.get('genotype', 0) + rv.get('only_a', 0) + rv.get('only_b', 0))
            cells.append('=' if d['identical'] else 'd%d%s (ED %d)' % (d['ed_vs_gw'], kind, d['ed_truth']))
        out.append('| ' + ' | '.join(cells) + ' |')
    return '\n'.join(out)


# ------------------------------------------------------------------ 7. scorer validation

def validate_scorer(rid):
    """The haplotype builder must reproduce truth.fa from the truth records (ED 0), and must see a
    perturbation (ED > 0). Uses region.py's own truth-record choice, written as a phased VCF."""
    import region as _region
    reg = load_region(rid)
    a1, b1, c = reg['span_start'], reg['span_end'], reg['contig']
    paths = config.data_paths(c)
    rows, applied = _region.choose_truth(paths, c, a1, b1, _region.NULL_LOG)
    wd = os.path.join(STAGE3, 'validate', rid)
    os.makedirs(wd, exist_ok=True)
    merged = collections.OrderedDict()
    for h in (1, 2):
        for r in applied[h]:
            merged.setdefault((r['pos'], r['ref'], r['alt']), set()).add(h)
    ploidy = ploidy_of(reg)

    def write(path, keys, gtf):
        with open(path[:-3], 'w') as f:
            f.write('##fileformat=VCFv4.2\n##contig=<ID=%s,length=%d>\n' % (c, contig_length(c)))
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n')
            for i, k in enumerate(keys):
                gt, ps = gtf(i, k)
                f.write('%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT:PS\t%s:%s\n' % (c, k[0], k[1], k[2], gt, ps))
        run([config.BCFTOOLS, 'sort', '-Oz', '-o', path, path[:-3]])
        run([config.BCFTOOLS, 'index', '-f', '-t', path])
        return path

    keys = sorted(merged)

    def gt_of(k):
        hs = merged[k]
        if ploidy == 1:
            return '1' if 1 in hs else '0'
        return '%d|%d' % (1 if 1 in hs else 0, 1 if 2 in hs else 0)

    out = {'region_id': rid, 'truth_records': len(keys)}
    v = write(os.path.join(wd, 'truth.vcf.gz'), keys, lambda i, k: (gt_of(k), 1))
    out['truth_ed'] = span_haplotypes(v, reg)['ed']
    # unphased: every het its own unit; the search must find ED 0 again (when it is exact)
    v = write(os.path.join(wd, 'truth_unphased.vcf.gz'), keys, lambda i, k: (gt_of(k).replace('|', '/'), '.'))
    out['truth_unphased_ed'] = span_haplotypes(v, reg)['ed']
    # perturbations: drop the largest record; swap one het's phase within the block
    if keys:
        big = max(keys, key=lambda k: abs(len(k[1]) - len(k[2])))
        v = write(os.path.join(wd, 'drop_largest.vcf.gz'), [k for k in keys if k != big],
                  lambda i, k: (gt_of(k), 1))
        out['drop_largest_ed'] = span_haplotypes(v, reg)['ed']
        hets = [k for k in keys if len(merged[k]) == 1]
        if hets and ploidy == 2 and len(keys) > 1:
            h0 = hets[len(hets) // 2]
            v = write(os.path.join(wd, 'swap_one_het.vcf.gz'), keys,
                      lambda i, k: ((gt_of(k)[::-1] if k == h0 else gt_of(k)), 1))
            out['swap_one_het_ed'] = span_haplotypes(v, reg)['ed']
    v = write(os.path.join(wd, 'empty.vcf.gz'), [], lambda i, k: ('0|0', 1))
    h = span_haplotypes(v, reg)
    out['empty_ed'], out['ed_ref'] = h['ed'], h['ed_ref']
    return out


# ------------------------------------------------------------------ CLI

def pin():
    os.makedirs(os.path.dirname(VG_PINNED), exist_ok=True)
    if not os.path.exists(VG_PINNED):
        shutil.copy2(VG_SOURCE, VG_PINNED)
        run(['codesign', '-s', '-', '-f', VG_PINNED])
    v = run([VG_PINNED, 'version']).stdout.splitlines()[0]
    if VG_COMMIT not in v:
        raise CallError('pinned vg is %s, not %s' % (v, VG_COMMIT))
    return v


def _run_one(kw):
    try:
        return call(**kw)
    except Exception as e:  # noqa: BLE001
        return {'region_id': kw['rid'], 'graph': kw['graph'], 'arm': kw.get('arm'), 'error': str(e)[:2000]}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--stage2-work', default=None)
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('pin')
    p = sub.add_parser('prep')
    p.add_argument('ids')
    p.add_argument('--graph', default='mc')
    p.add_argument('--force', action='store_true')
    p = sub.add_parser('call')
    p.add_argument('ids')
    p.add_argument('--graph', default='mc')
    p.add_argument('--arm', default=None)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--force', action='store_true')
    p = sub.add_parser('run')
    p.add_argument('ids')
    p.add_argument('--graphs', default=','.join(STAGE3_GRAPHS))
    p.add_argument('--arms', default='', help='comma list of diagnostic arms (hybrid50k,hybrid200kids,...)')
    p.add_argument('--arm-graphs', default='mc', help='the graphs the arms are run for (default mc)')
    p.add_argument('--compact', action='store_true', help="gzip each arm's reads.gaf and drop its hybrid.gfa")
    p.add_argument('--jobs', type=int, default=3)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--force', action='store_true')
    p = sub.add_parser('refinish')
    p.add_argument('vcfs', nargs='+')
    p = sub.add_parser('haps')
    p.add_argument('id')
    p.add_argument('vcf')
    p = sub.add_parser('gate')
    p.add_argument('ids', nargs='?', default='pilot')
    p.add_argument('--arms', default='')
    p.add_argument('--json', default=os.path.join(STAGE3, 'gate', 'gate.json'))
    p.add_argument('--tsv', default=None)
    p.add_argument('--md', default=None, help='write the haplotype-identity table (markdown) here')
    p.add_argument('--jobs', type=int, default=3)
    p = sub.add_parser('validate-scorer')
    p.add_argument('ids', nargs='?', default='pilot')
    a = ap.parse_args(argv)
    s2 = a.stage2_work
    if a.cmd == 'pin':
        print(pin())
    elif a.cmd == 'prep':
        for rid in region_ids(a.ids):
            print(json.dumps(prep(rid, a.graph, s2, a.force)))
    elif a.cmd == 'call':
        for rid in region_ids(a.ids):
            print(json.dumps(call(rid, a.graph, a.arm, a.threads, s2, a.force)))
    elif a.cmd == 'run':
        ids = region_ids(a.ids)
        jobs = []
        for rid in ids:
            for g in [x for x in a.graphs.split(',') if x]:
                jobs.append({'rid': rid, 'graph': g, 'arm': None, 'threads': a.threads, 'stage2': s2,
                             'force': a.force})
            for arm in [x for x in a.arms.split(',') if x]:
                for g in [x for x in a.arm_graphs.split(',') if x]:
                    jobs.append({'rid': rid, 'graph': g, 'arm': arm, 'threads': a.threads, 'stage2': s2,
                                 'force': a.force, 'compact': a.compact})
        # a graph without a span graph for a region (mafft_linsi has no L012272) is skipped
        jobs = [j for j in jobs if os.path.exists(span_gfa(load_region(j['rid']), j['graph']))
                or j['graph'] in ('truth', 'mc_unchop')]
        # shared inputs first, serially: window queries, named panel pieces, window reads, Stage 2 preps
        setups = {(j['arm'] or FINAL_ARM).replace('nodepth', '') for j in jobs}
        flanks = sorted({int(re.match(r'(?:gwsub|hybrid)(\d+)k', x).group(1)) * 1000
                         for x in setups if x.startswith(('gwsub', 'hybrid'))})
        if flanks:
            named_pieces(ids, flanks)
        for x in setups:
            if x.startswith('hybrid') and not x.endswith('orig'):
                F = parse_hybrid_setup(x)[0]
                for rid in ids:
                    try:
                        fetch_window_reads(rid, F)
                    except Exception as e:                      # one region must not stop the batch
                        log('fetch failed %s %s: %s' % (rid, x, str(e).splitlines()[-1][:300]))
        for j in jobs:
            x = (j['arm'] or FINAL_ARM).replace('nodepth', '')
            if os.path.exists(calls_path(j['graph'], j['rid'], j['arm'])) and not a.force:
                continue
            try:
                if x.startswith('stage2'):
                    prep(j['rid'], j['graph'], s2, False)
                elif x.startswith('gwsub'):
                    prep_gwsub(j['rid'], int(re.match(r'gwsub(\d+)k$', x).group(1)) * 1000)
                elif x.startswith('hybrid'):
                    F, kind = parse_hybrid_setup(x)
                    build_hybrid(j['rid'], j['graph'], F, tag=kind if kind in ID_TAGS else '')
            except Exception as e:  # noqa: BLE001
                log('prep failed', j['rid'], j['graph'], x, str(e)[:300])
        # an arm and its nodepth twin share one hybrid mapping: map before calling in parallel
        mapped = set()
        maps = []
        for j in jobs:
            x = (j['arm'] or FINAL_ARM).replace('nodepth', '')
            if os.path.exists(calls_path(j['graph'], j['rid'], j['arm'])) and not a.force:
                continue
            if x.startswith('hybrid') and not x.endswith('orig') and (j['rid'], j['graph'], x) not in mapped:
                mapped.add((j['rid'], j['graph'], x))
                F, kind = parse_hybrid_setup(x)
                maps.append((j['rid'], j['graph'], F, kind if kind in ID_TAGS else ''))

        def _map(m):
            try:
                map_hybrid(m[0], m[1], m[2], threads=max(1, a.threads), stage2=s2, tag=m[3])
            except Exception as e:  # noqa: BLE001
                log('map failed', m, str(e)[:300])
        with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
            list(ex.map(_map, maps))
        with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
            for r in ex.map(_run_one, jobs):
                if 'error' in r:
                    log('FAILED', r['region_id'], r['graph'], r.get('arm'), r['error'][:500])
                else:
                    log('ok', r['region_id'], r['graph'], r.get('arm') or FINAL_ARM, '%.1fs' % r['seconds'],
                        'records', r['shift']['records'], 'reads', r['reads_loaded'],
                        'enum', r['enumeration_haplotypes'], 'link', r['linkage_haplotypes'])
    elif a.cmd == 'refinish':
        for v in a.vcfs:
            print(v, json.dumps(refinish(v)))
    elif a.cmd == 'haps':
        reg = load_region(a.id)
        h = span_haplotypes(a.vcf, reg)
        h.pop('seqs')
        h.pop('written_seqs')
        print(json.dumps(h, indent=1))
    elif a.cmd == 'gate':
        ids = region_ids(a.ids)
        arms = a.arms.split(',') if a.arms else GATE_ARMS
        with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
            rows = list(ex.map(lambda rid: gate_one(rid, arms), ids))
        os.makedirs(os.path.dirname(a.json), exist_ok=True)
        with open(a.json, 'w') as f:
            json.dump(rows, f, indent=1)
        hdr, tsv = gate_tables(rows, arms)
        if a.tsv:
            with open(a.tsv, 'w') as f:
                f.write('\t'.join(hdr) + '\n')
                for line in tsv:
                    f.write('\t'.join(str(x) for x in line) + '\n')
        if a.md:
            with open(a.md, 'w') as f:
                f.write(gate_markdown(rows, arms) + '\n')
        for r in rows:
            print(r['region_id'], r['stratum'][:14], 'gwED', r['gw']['ed_truth'],
                  ' '.join('%s:%s(d%s,ED%s)' % (k, 'SAME' if v['identical'] else 'diff', v['ed_vs_gw'], v['ed_truth'])
                           for k, v in r['arms'].items()))
    elif a.cmd == 'validate-scorer':
        out = [validate_scorer(rid) for rid in region_ids(a.ids)]
        with open(os.path.join(STAGE3, 'validate', 'validate.json'), 'w') as f:
            json.dump(out, f, indent=1)
        for r in out:
            print(json.dumps(r))


if __name__ == '__main__':
    main()
