# VNTR region package

149 regions of the human pangenome (of 150 selected; one could not be cut out, see *Not packaged*), each
given as the haplotype sequences that the graph holds there, for multiple-alignment and
graph-construction experiments. Each region comes with two panels of haplotypes: the 34 that the
genotyping graph holds (`hap32.fa`), and the full HPRC panel that Minigraph-Cactus aligned
(`hprc.fa.gz`). This page says what the regions are, what is in each file, and what we would like back.
The results so far are in [../docs/results.md](../docs/results.md).

## Why these regions

In the HPRC v2.1 Minigraph-Cactus graph, a small set of long VNTRs (about 1% of VNTR regions) holds
half of the structural-variant false positives that short-read genotyping (`vg call`) makes. At these
loci the graph does not tangle; it misaligns the panel. The haplotypes' repeat arrays are stacked out of
phase in one giant bubble of ~3.6 bp nodes, so the pairwise alignment that the graph implies between two
haplotypes costs 1.56x the optimal pairwise alignment, against 1.09x at VNTRs of the same length where
genotyping works. The first hypothesis we are testing is that re-aligning each region's haplotypes with a
suitable multiple aligner, and building the region's graph from that alignment, makes the graph much
cleaner. The background is in [../docs/findings.md](../docs/findings.md) and the plan in
[../docs/plan.md](../docs/plan.md).

Every region is cut between two **anchor nodes**: CHM13 nodes of the graph that every panel path present
there passes through, none of them more than once, that no HG002 truth variant touches, and that lie at
least 200 bp outside the locus (the nearest such nodes; median flank 222 bp, but up to 42 kb where
variation or truth records extend far). So every sequence of a region starts with the same left-anchor
sequence and ends with the same right-anchor sequence: the flanks are shared, and a global (end-to-end)
alignment is the right model. Do not use local alignment (spoa's default mode, for example): it can
leave the shared flanks unaligned.

## Two panels, and what the baseline is

- **The full panel** (`hprc.fa.gz`): every haplotype of the HPRC v2.1 Minigraph-Cactus graph (its
  HG002-free evaluation build, below) that runs through the region from anchor to anchor, plus CHM13
  and GRCh38. That is up to 456
  haplotypes of 228 HPRC samples (median 452 per region). Minigraph-Cactus aligned all of them.
- **hap32** (`hap32.fa`): CHM13, GRCh38 and 32 haplotypes sampled from that graph by
  `vg haplotypes`. This is the graph `vg call` genotyped HG002 with; the sampling is personalised
  to HG002's reads, which is the standard vg pipeline. Within these regions, 4,620 of the 4,734
  sampled sequences are exact copies of one full-panel haplotype. The other 114 (in 32 regions)
  join two haplotypes inside the span, or copy a reference path.
- **The baseline `mc.gfa` is Minigraph-Cactus's alignment of the full panel, restricted to the
  hap32 paths.** It is not an alignment of the 34 hap32 sequences. Every node and edge of every
  `mc.gfa` is in the full graph with the same sequence, and every pair of hap32 paths implies the
  same alignment in both (`results/panel_premise.tsv`).
- **So there are two ways to replace it,** and we would like both (see *What we would like back*):
  - aligning `hap32.fa` alone: an easier problem (at the median hotspot, 30 distinct sequences
    against 374 in the full panel), and what a post-hoc fix of a sampled graph would do;
  - aligning the full panel and then keeping the hap32 rows: the like-for-like replacement for
    Minigraph-Cactus.

**Both panels are free of HG002.** Both come from the HPRC v2.1 *evaluation* graph
(`hprc-v2.1-mc-chm13-eval.gref.gbz`):
- Its samples are 228 HPRC individuals, CHM13 and GRCh38 (checked from the GBWT sample list).
- HG002 and its parents HG003 and HG004 are absent, and so are their aliases NA24385, NA24149 and
  NA24143.
- The graph's `gref_CHM13` path is excluded from both panels. It is vg's reference-cover path, a
  chimera stitched from other haplotypes, not a haplotype Minigraph-Cactus aligned. We checked that
  none of its pieces contains an anchor node in any region.

**Do not add HG002's, HG003's or HG004's sequences from any other source** (HPRC or T2T
assemblies, the non-eval HPRC graph, read data). A graph that contains the evaluated individual's
own haplotypes has an inherent advantage when genotyping that individual, and could not be scored
fairly.

## Strata: what each region is for

| stratum | n | role | what it is |
|---|---|---|---|
| `hotspot_vntr` | 40 | **test** | VNTRs (period >= 7) with the most `vg call` SV false positives (>= 12 here, <= 3 per contig). These are the regions to fix. |
| `control_vntr_matched` | 40 | **must not regress** | VNTRs with no false positive or false negative, matched to the hotspots on length and period bin. A new method should leave them as simple as they are. |
| `control_vntr_correct` | 25 | **must not regress** | VNTRs that carry a truth SV that `vg call` genotypes correctly, with no false positive. |
| `hotspot_other` | 24 | **test** (secondary) | Non-VNTR loci with >= 3 false positives: mostly short tandem repeats (period < 7), a few low-complexity, satellite and non-repeat loci. |
| `control_nontr_sv` | 20 | **must not regress** | Non-repeat loci with a truth SV that `vg call` genotypes correctly. The graph is already right here. |

Nine regions are the pilot set named in the findings (their `notes` field starts with `pilot:`):
alleles present but misaligned (L009656, L014297, L015415, L005990); under-aligned, i.e. homologous
sequence that the graph never aligned (L012184, L011138, L012272); alleles absent from the panel
(L002013); and a negative control (L016870).

All loci lie inside the HG002 T2T-Q100 v1.1 SV benchmark on chr1-22 and chrX, and were at most 20 kb
long when selected.

## Files

```
regions/
  README.md            this page
  regions.tsv          one row per region (the index)
  vntr_ge1kb.bed       the truth-free target list for scaling up (coordinates only, no sequences)
  <region_id>/
    region.json        everything known about the region
    hap32.fa           the input: every hap32 panel path that spans the region, anchor to anchor
    hap32.fragments.fa panel path pieces that enter the region but do not span it (may be empty)
    hprc.fa.gz         the same span for the haplotypes of the full HPRC v2.1 panel (where available)
    mc.gfa             the baseline: the Minigraph-Cactus graph of the region, with one path per hap32.fa sequence
    truth.fa           HG002's two haplotypes over the span -- EVALUATION ONLY, see below
    calls.tsv          vg, PanGenie and truth records over the span, for context
```

The directory is 92 MB: `mc.gfa` 56 MB, `hap32.fa` 20 MB (5,026 sequences), `hprc.fa.gz` 9.6 MB
(65,199 sequences), the rest 6 MB. Spans run from 658 bp to 55 kb (median 3.0 kb; hotspots median
5.7 kb).

`region_id` is the census locus id (`L` + six digits), or, for 16 matched controls that have no call
or truth record and so no locus, the census VNTR region id (`TR` + six digits).

**Coordinates** are CHM13v2.0. In `regions.tsv` and `region.json` they are 1-based and inclusive. The
*core* (`core_start..core_end`) is the VNTR or locus interval. The *span* (`span_start..span_end`)
runs from the first base of the left anchor node to the last base of the right anchor node, and it is
exactly the CHM13 sequence in `hap32.fa`. `vntr_ge1kb.bed` is a BED file (0-based, half-open).

### regions.tsv

Tab-separated, one header line. `.` means not applicable or unknown.

| column | meaning |
|---|---|
| `region_id` | census locus id; also the directory name |
| `stratum` | one of the five strata above |
| `contig` | CHM13 contig |
| `core_start`, `core_end` | the VNTR / locus interval, 1-based inclusive |
| `span_start`, `span_end` | the anchor-to-anchor interval, 1-based inclusive |
| `anchor_left`, `anchor_right` | anchor node ids in the Minigraph-Cactus graph, with the orientation CHM13 walks them (`123+`) |
| `period`, `motif`, `copies` | the repeat annotation of the census (tandem-repeat finder); `.` for non-VNTR loci |
| `n_hap32` | number of sequences in `hap32.fa` (34 when every panel path spans) |
| `n_hprc` | number of HPRC sample haplotypes in `hprc.fa.gz` (not counting CHM13 and GRCh38), 0 if none |
| `vg_fp`, `vg_fn` | `vg call` SV (>= 50 bp) false positives and false negatives in the core (truvari) |
| `pg_fp`, `pg_fn` | the same for PanGenie (which genotyped the full HPRC panel, not hap32) |
| `in_benchmark` | 1 if the core lies wholly inside the SV benchmark BED |
| `notes` | `;`-separated remarks: pilot role, haplotypes that do not span, N runs, anchor fallbacks, full-panel status |

### hap32.fa -- the input

One record per panel path that runs through the whole region, one line per sequence, upper case. The
panel ("hap32") is the graph that `vg call` used: CHM13, GRCh38 and 32 haplotypes sampled from the
HPRC v2.1 eval graph (`vg haplotypes`). The sampler may switch between real haplotypes, hence the name
`recombination`; within these regions almost every sampled sequence is a copy of one full-panel
haplotype (see *Two panels*). Every sequence starts with the left anchor's sequence and ends with the right
anchor's sequence (both are in `region.json`: `anchor_left_seq`, `anchor_right_seq`). Exactly one
record is CHM13, and it is the first.

Record names are the graph's path names, unchanged, and are unique within a file:

- `CHM13#0#chr4` -- the reference;
- `GRCh38#0#chr4[58921381]` -- GRCh38; the bracket is the GRCh38 path fragment's start offset, part of
  the name;
- `recombination#12#chr4#0` -- sampled haplotype 12 (the last field is the path fragment number).

The description after the name is `len=<bp>`. A handful of regions contain sequences with runs of `N`
(assembly gaps carried into the graph as alleles); `region.json` counts them
(`hap32.n_sequences_with_N_runs`). They are given as they are. Identical sequences are kept as separate
records: the number of haplotypes carrying an allele matters for genotyping.

### hap32.fragments.fa

Pieces of panel paths that enter the region but do not run through it (a haplotype that has a gap in
the region, or that leaves it through a translocation-like edge). Names are
`<path name>:<class>` with class `enters_L` (starts at the left anchor and stops inside), `exits_R`
(starts inside and leaves through the right anchor) or `internal`; the description gives the first and
last node. They are not part of the required output; you may use them (e.g. `mafft --addfragments`) or
ignore them. Pieces that are only an anchor node are left out and counted in `region.json`.

### hprc.fa.gz

The same anchor-to-anchor span for the haplotypes of the full HPRC v2.1 eval graph (456 sample
haplotypes of 228 samples; HG002, HG003 and HG004 are not among them, and `gref_CHM13` is left out).
This is the panel Minigraph-Cactus aligned, and the input for the full-panel arm. The first records are
CHM13 and GRCh38 with their `hap32.fa` names; the sample haplotypes follow as `hprc#1`, `hprc#2`, ...
The graph query that produces them does not report sample names, so `hprc#k` carries no sample identity:
records with identical sequence are grouped (`group=<g> group_size=<n>` in the description, groups by
decreasing size) and `same_as_hap32=` flags sequences also present in `hap32.fa`. Only haplotypes that
pass through both anchor nodes are included, so that the flanks stay identical; a full-panel haplotype
that carries a variant on an anchor node, has a contig end inside the span, or leaves it is left out,
and `region.json` (`hprc.anchor_touching_walks_not_spanning`) counts those. It is present for all 149
regions: a median of 452 sample haplotypes per region (fewer on chrX, where the male samples carry one
copy: median 340), 267-457. One region (L012272) has 457, one more than there are haplotypes: an
assembly there carries two contigs through the span. The candidate graph you return still needs only
the hap32 paths (see *What we would like back*). `tools/panel.py union` turns `hap32.fa` plus
`hprc.fa.gz` into one FASTA of distinct sequences to align, with a map back to every record.

### mc.gfa -- the baseline

GFA 1.0: `S` lines with the Minigraph-Cactus node ids and sequences, `L` lines for every edge the paths
use, and one `P` line per `hap32.fa` record, with the same name, spelling exactly its sequence (this is
asserted when the package is built). It is the graph those paths induce: nodes and edges that only
fragments use are left out. This is "M0", what every method is compared with. It is the full-panel
Minigraph-Cactus alignment restricted to the hap32 paths (see *Two panels*), so its node cuts include
those made by haplotypes that are not in `hap32.fa`. Its node ids are large and meaningless to you.

### truth.fa -- held out

`HG002#1` and `HG002#2`: the two haplotypes of HG002 over the same span, built by applying the phased
T2T-Q100 v1.1 truth (small and structural variants) to CHM13. **Please do not use truth.fa, the
`truth_*` rows of calls.tsv or the truth-derived fields of region.json (`truth`,
`closest_hap32_to_truth`, `called_haplotypes_vs_truth`, and the like) to build or tune graphs.** They
are there so that the evaluation can run anywhere; a method that has seen them cannot be scored
fairly. HG002 and its parents are in neither panel, so it is a fair test of the panel's alleles; keep
it that way (see *Two panels*).

### calls.tsv

Context, not input. One row per record overlapping the span, sorted by position: `source` is
`truth_stvar`, `truth_smvar` (the T2T-Q100 structural and small-variant VCFs), `vg` or `pangenie`; then
`pos`, `end`, `ref_len`, `alt_len`, `len_diff`, `gt`. For calls, `truvari` is the status of an SV-sized
call (`TP`/`FP`) and `aardvark` that of a small one (`-` where a record was not scored). For truth rows,
`vg_truvari` / `pg_truvari` say whether each caller found the truth SV (`TP`/`FN`), `applied_h1` /
`applied_h2` whether it went into truth.fa's haplotypes, and `in_sv_benchmark`, `svtype` and
`TRFperiod` annotate it. `overlaps_span_edge` is 1 for a record that reaches past an anchor.

### region.json

The `regions.tsv` fields plus: `role`, the anchor sequences, flank lengths, benchmark coverage;
`tr_annotation` (census period, motif, copies, GC, telomere distance, tandem scans of CHM13 and a
repeatmasker summary); `hap32` (number of spanning paths, distinct sequences, fragments, lengths,
haplotypes that do not span); `baseline_mc` (nodes, edges, bp, nodes per CHM13 kb, node lengths of
`mc.gfa`); `hprc` (how `hprc.fa.gz` was made, or why it is missing); `calls_in_span`; truth-derived
fields (held out, see above); `files` (sizes and record counts); `provenance`.

### vntr_ge1kb.bed

All VNTRs of period >= 7 that lie in the SV benchmark and are at least 1 kb long in CHM13: 3,111
regions, 6.13 Mb, which hold 46% of all `vg call` SV false positives genome-wide. The selection uses
length only, no truth. Columns: contig, start, end (0-based half-open), census region id, period,
motif, vg FP, vg FN (the last two for reference). It has no sequences: it is the list to scale a method
up to once it works on the 150 regions (`tools/package_regions.py build --loci regions/vntr_ge1kb.bed`
packages any subset of it the same way).

## Not packaged

**L000601** (`hotspot_other`, chr1:103,777,433-103,780,251, 24 vg false positives) sits at the end of a
~94 kb tandemly duplicated segment in which HG002 carries 94 kb deletions of different copies (truth
records from 103.52 to 103.78 Mb). Every CHM13 node for at least 170 kb to its left is touched by one
of them, so no anchor pair that keeps the truth haplotypes' flanks identical exists within 320 kb: it is
a copy-number locus at the 100 kb scale, not a region that fits this format.

In 17 of the 149 regions fewer than 34 panel paths span (`n_hap32` < 34): a sampled haplotype or GRCh38
has a gap there or does not reach the anchors. `notes` names them; their pieces are in
`hap32.fragments.fa`.

## What we would like back

For each region and each method, a graph with one path per `hap32.fa` record, from each of the two
arms:

```
candidates/<method>/<region_id>.gfa            arm 1: hap32.fa aligned on its own
candidates/<method>/<region_id>.msa.fa         (alignment-based methods: the aligned rows, '-' for gaps)
candidates/<method>__all/<region_id>.gfa       arm 2: the full panel aligned, then restricted to the hap32 rows
candidates/<method>__all/<region_id>.msa.fa
```

- **Same names, same sequences.** Every `hap32.fa` record must be a path of the same name that spells
  its sequence exactly (N runs included, so put them in nodes of their own). A graph that changes a
  sequence is rejected. Extra paths (a consensus, say) are allowed and ignored.
- The graph should be a DAG on the handles the paths use (a cycle is reported and costs some metrics).
- Node ids and node lengths are up to you. So is whether the graph is compacted; we compact nothing.
- If you build the graph from an MSA with vg (`vg construct -M msa.fa -F fasta`, then `vg mod -u` or
  `vg unchop` and `vg convert -f`), vg may write the PanSN-style names (`sample#hap#contig`) as `W` lines
  and the GRCh38 bracket as the W line's start. That is fine: the evaluation maps those back to the
  names. `tools/msa_graph.py` builds the same graph without vg and checks every path.
- Please say, per method, the exact command lines, the input order (partial-order aligners depend on
  it) and anything you dropped or masked.

**Arm 2, step by step** (everything here runs offline, from this repository alone):

```bash
python3 tools/panel.py union regions/L014297 work/union/L014297.fa     # distinct sequences of hap32 + full panel
#   -> work/union/L014297.fa (rows u0001.., CHM13 first) and work/union/L014297.map.tsv (which records each row stands for)
<your aligner> work/union/L014297.fa > L014297.full.msa.fa             # rows keep their u0001.. names
python3 tools/panel.py project L014297.full.msa.fa work/union/L014297.map.tsv regions/L014297/hap32.fa \
    candidates/<method>__all/L014297.msa.fa                            # hap32 rows, all-gap columns dropped
python3 tools/msa_graph.py candidates/<method>__all/L014297.msa.fa regions/L014297/hap32.fa \
    candidates/<method>__all/L014297.gfa
```

If you also send the full-panel MSA (gzipped), we can build and score the full-panel graph itself:
`tools/panel.py panel-graph FULL.msa.fa MAP.tsv regions/<id> OUT.gfa`, then
`tools/evaluate.py regions/<id> OUT.gfa --panel regions/<id>/hprc.fa.gz`.

## How it will be scored

The stages are in [../docs/plan.md](../docs/plan.md) ("How a candidate graph is judged"). In short:

- **Stage 0, graph only**:
  - nodes per kb and node lengths;
  - the pairwise alignment cost that the graph implies between two haplotypes, over their optimal
    pairwise cost (unit and affine), and excess edits per kb;
  - SV-sized pieces per haplotype;
  - k-mer redundancy (the same 21-mer on different nodes), which predicted read placement best in
    our local re-mapping;
  - whether each HG002 haplotype can be spelled by a path, and at what edit cost;
  - "the truth written by the graph", scored with truvari raw, after `truvari refine`, and as one
    whole-span phab region.

  Target: the hotspots move to the matched controls' distribution (cost/optimum <= 1.1) and the
  controls do not change.
- **Stage 1**: placement redundancy of the region's existing short reads in the graph.
- **Stages 2-4**: local re-mapping, local genotyping with `vg call`, and splicing the rebuilt regions
  back into a whole-contig graph (done by us).

Stage 0 is one command per region and graph ([../tools/METRICS.md](../tools/METRICS.md) defines every
number):

```bash
cd vntr-eval                                                                 # the repository root
python3 tools/evaluate.py regions/L014297 regions/L014297/mc.gfa             # the baseline -> results/mc/L014297.json
python3 tools/evaluate.py regions/L014297 candidates/<method>/L014297.gfa    # -> results/<method>/L014297.json
python3 tools/evaluate.py regions/L014297 candidates/<method>__all/L014297.gfa
```

**What it needs.** Evaluation needs `regions/` and `tools/` of this repository, Python 3.9 or later
(standard library only; tested with 3.9.6 and 3.14), `make` and a C compiler. The first run builds
`tools/bin/libvntreval.so` itself. It needs neither vg nor mafft. Two parts need the big data, which is not in this repository:

| part | needs | without it |
|---|---|---|
| truvari scoring of the truth written by the graph (`truth.truvari`) | the CHM13 FASTA, the normalised truth VCF and the SV benchmark BED of the evaluation repository (`VNTR_EVAL_DIR`), and truvari | reported as `needs_big_data`; pass `--no-truvari` to skip it quietly |
| Stage 1 (`--reads`) | the hap32 graph and read databases (GBZ-Base, GAF-Base) | reported as `needs_big_data` |

Everything else is computed from the package, including the truth-to-graph edit distances
(`truth.h1_d_graph`, `h2_d_graph`). We tested this from a clean shell on a copy of `regions/` and
`tools/` only, with `HOME` pointing nowhere, `PATH=/usr/bin:/bin` (macOS's Python 3.9.6) and no
evaluation repository. `evaluate.py` on `mc.gfa`, on a candidate and with `--panel` gave Stage 0
numbers identical to ours, and the arm-2 steps above reproduced our `unit_aware__all` graph.
We run the truvari scoring and Stage 1 on our side.

`tools/viewer.py regions/<id> --no-data --candidate candidates/<method>/<id>.gfa` draws a region's
alignment and graphs as one HTML page from the package alone. It needs mafft for the all-rows MSA,
and Bandage for the graph images (or `--no-bandage`). `python3 tools/config.py` prints where every
tool looks for data and binaries.

## Rebuilding the package

`tools/package_regions.py` builds it from the graph databases of the evaluation repository (paths in
`tools/config.py`). **This needs the big data** (the hap32 and full HPRC graphs, the truth, the calls,
`gbz-base` and `vg`), so it runs only on our side; a collaborator never needs to rebuild:

```bash
python3 tools/package_regions.py build --loci <census strata.tsv | regions/regions.tsv | a BED> --work W   # packages
python3 tools/package_regions.py hprc-fetch --loci <same> --work W   # full-graph queries (~1 min, ~11 GB each)
python3 tools/package_regions.py hprc --work W                       # adds hprc.fa.gz from the cached queries
python3 tools/package_regions.py index         # rewrites regions.tsv
python3 tools/package_regions.py validate      # re-checks every package
python3 tools/package_regions.py bed --vntr-regions <census vntr_regions.tsv>                   # vntr_ge1kb.bed
```
