# Stage 2: local re-map of the existing reads (pilot, 19 regions)

The question: do HG002's short reads map better to a realigned local graph than to the
Minigraph-Cactus (MC) local graph?

We re-mapped the reads that the genome-wide giraffe run placed in each region. For each region we
mapped them to MC's graph, to the realigned candidates from both arms, and to a graph of HG002's
own two haplotypes. That last graph is for evaluation only; call it the "truth graph".

## Answer

**Yes, but only somewhat, and it depends on the method. The headroom is limited because much of
the ambiguity is irreducible.**

### mafft L-INS-i (hap32 arm) is the only method better than MC at every hotspot, without losing at the matched controls

Hotspot VNTRs, n = 9. Numbers are medians; paired differences are against MC.

| Measure | MC | L-INS-i (better / worse vs MC) | Truth graph |
|---|---|---|---|
| Reads in HG002's core with MAPQ < 5 | 0.817 | 0.671 (8/0; paired −0.075) | 0.518 |
| Reads that fit HG002's haplotypes and land in the core at MAPQ ≥ 5 | 0.172 | 0.327 (8/0; +0.078) | 0.495 |
| Uniquely placeable core reads placed at the right position at MAPQ ≥ 5 (6 evaluable hotspots) | 0.384 | 0.483 (6/0; +0.196) | – |
| Core depth from MAPQ ≥ 5 reads | 2.8x | 5.2x (8/0) | 13.4x |

- Confident placements are not more often wrong. Among core reads at MAPQ ≥ 5, the share
  placed wrong is 0.131 for L-INS-i against 0.132 for MC (paired −0.033; 4 regions better, 2 worse).
- The gap closed is modest. L-INS-i closes about a quarter of the gap between MC and the truth
  graph on the MAPQ < 5 fraction (paired −0.075 of −0.272).
- Most of the MAPQ loss is inherent. Even the truth graph, which holds only HG002's two alleles,
  leaves 52% of core reads at MAPQ < 5. That is tandem ambiguity, which no realignment removes.
  This is the local counterpart of the 56% "tandem" read class in [findings.md](../docs/findings.md).
- At the 7 matched VNTR controls L-INS-i is neutral: median MAPQ < 5 is 0.053 against 0.039 for
  MC (paired −0.007; 4 better, 2 worse).
- **L-INS-i fails at one non-repeat control, L006556.** MAPQ < 5 rises from 0.007 to 0.255
  (all reads on the span: 0.01 → 0.20). This is the Stage 0 regression already known there
  (cost/optimum 1.19 against 1.00).

### The full-panel arm helps no more than MC

The full-panel arm is the like-for-like replacement for MC.

- **unit_aware__all** is 4 better / 5 worse than MC on MAPQ < 5 at hotspots (median 0.751), and
  worse at 6 of 7 matched controls.
- **poa_spoa__all and poa_abpoa__all** are worse than MC at every hotspot where they exist
  (0/4 and 0/5). At controls they do real damage: L007172 goes from 0.32 to 0.78 with spoa__all.
- **poa_abpoa_mc__all** (abPOA with Cactus's settings) is mixed: 5 better / 3 worse at hotspots
  and 4/2 at controls.
- **mafft_linsi__all** exists at only one hotspot (L015415), where it equals L-INS-i.

So on read placement, as at Stage 0, aligning the full panel and projecting onto hap32 is worse
than aligning hap32 alone.

### Stage 0 cost/optimum does not predict read placement; k-mer redundancy does

**The Stage 0 headline metric does not predict read placement.** unit_aware has the lowest
median cost/optimum at the hotspots (1.023 against 1.335 for MC), yet it barely changes MAPQ
there (6 better / 3 worse) and is worse at 6 of 7 matched controls.

What predicts placement is **k-mer redundancy**, the same 21-mer sitting on different nodes.
Across 142 pairs of region and candidate graph, each compared with MC for that region:

| Change against MC (Spearman with the change in core MAPQ < 5) | All 142 pairs | Hotspots | Controls |
|---|---|---|---|
| k-mer extra fraction | 0.75 | 0.54 | 0.81 |
| nodes per kb | 0.65 | 0.44 | 0.63 |
| cost/optimum | 0.26 | −0.08 | 0.18 |

Examples:

- At L009656, unit_aware has cost/optimum 1.02 but k-mer extra 0.53, the same as MC's 0.54, and
  it maps worse than MC. poa_abpoa_mc__all has cost/optimum 1.29 but the lowest k-mer extra (0.40),
  and maps best.
- The control regressions all raise k-mer extra: L007687 goes from 0.22 to 0.28 with abPOA and
  to 0.37 with unit_aware__all; L006556 goes from 0.01 to 0.23 with L-INS-i.

A realigner that merges identical sequence onto shared nodes is what helps the reads. An
equal-cost alignment that leaves identical sequence on separate nodes does not. **Stage 0 should
rank candidates on k-mer redundancy alongside cost/optimum.**

### One gross fault is fixed outright: L012272

L012272's 55 kb span holds a 5 kb VNTR and 42 kb of anchor flank.

- **The genome-wide mapping puts 92% of the span's reads at MAPQ < 5.** The local MC graph
  reproduces this (91%), so the ambiguity is inside MC's own subgraph. The flank MAPQ ≥ 5 depth is
  4.0x.
- **unit_aware's realignment of the flanks fixes it.** Span MAPQ < 5 drops to 8% (7% with
  unit_aware__all; 5% in the truth graph), and the flank MAPQ ≥ 5 depth rises to 30.4x.
- k-mer extra falls from 0.52 to 0.12. This is the "parallel duplicated paths" fault from
  [plan.md](../docs/plan.md), over tens of kb.
- Only unit_aware, unit_aware__all, mafft_fftns2 and poa_abpoa_mc have graphs there, and only the
  first two were mapped.

### Read-set limits that no graph changes

- **The core depth deficit belongs to the read set.** Raw core/flank depth is the same in every
  local graph: median 0.62 at hotspots, 0.25 at L012184, 0.98 at controls. The reads missing from
  hotspot cores were placed elsewhere or left unmapped by the genome-wide run, so they are not in
  the GAF-Base window. A local re-map cannot bring them back.
- **Haplotype compatibility does not tell the graphs apart.** 98% of core reads lie on a single
  hap32 path in every graph at hotspots, MC included. What does differ is which paths:
  - The truth-closest pair of hap32 paths explains 38% of the core reads in MC and 49% in L-INS-i
    (7 better / 1 worse).
  - The best pair explains 52% in MC and 64% in L-INS-i.
  - 95% of compatible core reads need 13 haplotypes in MC and 10 in L-INS-i.
  - At the VNTR hotspots the truth-closest pair is never the best pair, in any graph (it is at 23
    of 105 control cells). Reads are spread over many paths because tandem ambiguity puts them on
    other haplotypes at shifted positions.

## Method

### Regions

The 7 pilot loci named in findings.md, plus L011138 and L012272 (both named there) and L001909:

- 9 VNTR hotspots: L009656, L014297, L015415, L005990, L012184, L002013, L011138, L012272, L001909.
- 1 non-VNTR hotspot: L016870.
- 7 matched VNTR controls, the census matches of the pilot hotspots: L016124, L007687, L007172,
  L011430, L007065, L004021, L014111.
- 2 non-repeat SV controls: L000506, L006556.

### Reads

For each region the reads come from the GAF-Base: all alignments overlapping the anchor-to-anchor
subgraph (`--between`), plus the CHM13 window span ± 1 kb (`--interval`, context 100).

- **Rebuilding the reads.** GAF-Base carries no read sequence, but each read can be rebuilt from
  the path it aligned to and its `cs` string. Every genome-wide alignment covers the whole read
  (qs = 0, qe = qlen). Base qualities come from `bq:Z`.
- **Counts.** 36,395 reads were rebuilt and none dropped. Each rebuilt read was checked against
  its path length and read length.
- **Mates.** Mates share one read name and are told apart by the `fn`/`fp` tags. 95.4% of reads
  came with their mate (17,362 pairs) and were mapped paired. The 1,671 without a mate were mapped
  single-end.
- **Fragment length.** One fragment-length distribution is used for every graph. It was estimated
  from 5,795 pairs whose mates both lie on colinear CHM13 walks: mean 402, sd 166 (median 372).

### Graphs

Every graph got the same 1 kb of CHM13 flank added on each side, as one node joined to every path.
Reads that run past the anchors, and mates in the flanks, therefore map identically whatever lies
between the anchors.

**Paths must be written as haplotype W lines.** With P lines, `vg autoindex` finds no haplotypes
and builds a greedy path cover. Giraffe's extension is then not haplotype-aware: at L014297, core
reads on a single haplotype fell to 0.70, against 0.98 genome-wide. With W lines the local MC
mapping reproduces the genome-wide one.

The graphs mapped:

- **MC:** mc.gfa.
- **hap32 arm:** mafft_linsi, unit_aware, poa_spoa and poa_abpoa. These are the best by Stage 0
  cost/optimum; both POAs were cheap enough to map.
- **Full-panel arm (`__all`):** unit_aware__all (the best of that arm at Stage 0),
  poa_abpoa_mc__all (the second with near-full coverage), poa_spoa__all, poa_abpoa__all, and
  mafft_linsi__all where it exists.
- **Truth graph:** mafft G-INS-i of truth.fa, or FFT-NS-2 above 15 kb, induced with msa_graph.py.

That makes 180 mappings in all: `vg autoindex -w sr-giraffe`, then `vg giraffe` with the default
preset, 2 threads, and GAF in segment coordinates.

**The truth graph is an MSA-induced graph too, so it is not a strict ceiling.** Where mafft
misaligns HG002's two alleles, it inherits the misalignment. This happens at L006556 (MAPQ < 5
0.12 against 0.01 for MC) and L016870 (0.21 against 0.04 for L-INS-i). Elsewhere it is the best
graph or within 0.01 of the best.

### Validation: the local mapping reproduces the genome-wide one

`orig` is the genome-wide giraffe alignment of the same reads to the whole hap32 graph, scored with
mc.gfa's paths.

| Measure, hotspots and VNTR controls | Genome-wide vs local MC |
|---|---|
| Core MAPQ < 5, median at hotspots | 0.795 vs 0.817 (paired difference 0.003) |
| Largest per-region difference | 0.04 (L002013, L005990, L016870) |
| Core reads on one haplotype, hotspots | 0.983 vs 0.983 |

So at VNTRs the MAPQ loss is local, as findings.md reported (97% of MAPQ-0 alternatives local).

The two non-repeat controls are the exception: genome-wide MAPQ < 5 is 0.061 and 0.083, against
0.020 and 0.007 locally. That is the genome-wide competition a local map cannot see.

### Definitions

All are in the tool's docstring.

- **Core.** The census core interval, located in each haplotype and in HG002's haplotypes by
  aligning the 300 bp of CHM13 on either side semi-globally. Per graph, each node base gets the
  fraction of its paths that have that base inside their core.
- **Truth-core read.** A read with ≥ 10 bases in HG002's core in its truth-graph alignment. This
  read set is the same whichever graph is scored.
- **Truth-consistent.** The truth-graph alignment has at most 1 high-quality edit: a substitution
  at base quality ≥ 20, an internal indel, or a soft clip with ≥ 5 bases of quality ≥ 20.
- **Placed right.** On the hap32 path that copies the read's truth haplotype, within 5 bp and on
  the same strand. Evaluated for truth-consistent reads with truth-graph MAPQ ≥ 30 whose truth
  haplotype has an exact copy, or a near copy (edit ≤ 10) that spells the read's interval exactly.
- **Depth.** Aligned bases over the mean HG002 core (or flank) length.

## Caveats

- **Local mapping loses genome-wide competition.** Real alternatives elsewhere in the genome are
  absent. This matters little at the VNTRs (the table above), but it shows at the non-repeat
  controls. Reads that the genome-wide run placed elsewhere or left unmapped are not in the read
  set at all.
- **This is a pilot.**
  - 9 VNTR hotspots and 7 matched controls, so the better/worse counts are small.
  - Some methods lack graphs at the largest hotspots: spoa at L011138, L012272 and L001909; the
    `__all` POAs at 3-5 hotspots; mafft_linsi__all at all but L015415.
  - Per-group medians are over the regions where each graph exists; paired counts are against MC
    in the same regions.
- **Placement is evaluable at 6 of 9 hotspots.** L002013, L012272 and L001909 have no hap32 path
  within 10 edits of a truth haplotype. At L005990 and L012184 only 24 core reads are evaluable.
- **MAPQ is giraffe's.** A different mapper or preset could rank the graphs differently.
- **Some metrics are indirect.** Haplotype compatibility and the best pair are judged on each
  read's primary alignment only. Depth uses graph-base core fractions; where a graph's flank is
  itself ambiguous (L012272), the MAPQ ≥ 5 core/flank ratio is confounded, so the absolute depths
  are reported too.
- **Genome-wide (`orig`) flank depths are not comparable** (bases outside the window count as
  flank), so they are left blank.
- **Stage 0 numbers here** (cost/optimum, k-mer extra, nodes per kb) were computed with
  `evaluate.py --skip inflation,truth` into scratch. results/ holds only mc and poa_abpoa_mc so
  far; the Evaluate phase will fill the rest. The TSV's `stage0_cost_over_opt` column used the
  scratch copies.
- **Resources.** The machine was loaded (load 12-27 from other agents); runs used 2 jobs × 2
  threads. Mapping took 2-208 s per graph (37 min in total) and indexing under 1 s per graph.
  Work files (321 MB) are in scratch.

## Reproduce

```
W=<work dir>            # default $VNTR_WORK/remap (git-ignored)
python3 tools/remap_local.py run L009656 L014297 ... --work $W --jobs 2 --threads 2   # fetch, fraglen, map, analyse
python3 tools/remap_local.py table   <ids> --work $W --tsv results/stage2_remap.tsv [--stage0-dir DIR]
python3 tools/remap_local.py summary <ids> --work $W [--stage0-dir DIR] [--md tables.md]
```

- `--methods` (default: mc, the 4 hap32 candidates, 5 `__all` candidates, truth) picks the
  graphs.
- Per region, the work directory holds:
  - `reads_{1,2,se}.fq`, `reads.tsv` (with the genome-wide alignment of each read) and
    `reads.json`;
  - per graph: `<method>/local.gfa`, the indexes, `reads.gaf` and `per_read.<method>.tsv`;
  - `stage2.json` with every metric.
- results/stage2_remap.tsv has one row per region × graph (228 rows, including `orig` and
  `no_graph` rows). Its columns are named as in the tables below.

## Appendix: every metric, per region and per group

- **Cells** in the region rows are the graph's value.
- **Group rows** give the median, then (paired median difference against MC; regions better / worse).
- **Strata:** hs_vntr = hotspot VNTR, hs_other = non-VNTR hotspot, c_vntr_matched = matched
  control, c_nontr_sv = non-repeat SV control.

#### Stage 0: all-pairs cost/opt of the graph (for reference) (`stage0_cost_over_opt`, low is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 1.335 | 1.335 | 1.099 | 1.023 | 1.050 | 1.062 | 1.031 |  | 1.114 | 1.127 | 1.286 | - |
| L014297 | hs_vntr | 1.303 | 1.303 | 1.051 | 1.022 | 1.040 | 1.067 | 1.056 |  | 1.108 | 1.133 | 1.452 | - |
| L015415 | hs_vntr | 1.184 | 1.184 | 1.049 | 1.006 | 1.016 | 1.011 | 1.015 | 1.022 | 1.024 | 1.020 | 1.059 | - |
| L005990 | hs_vntr | 1.211 | 1.211 | 1.063 | 1.014 | 1.027 | 1.044 | 1.042 |  | 1.220 | 1.184 | 1.245 | - |
| L012184 | hs_vntr | 1.599 | 1.599 | 1.143 | 1.136 | 1.081 | 1.091 | 1.183 |  |  |  | 1.640 | - |
| L002013 | hs_vntr | 1.406 | 1.406 | 1.100 | 1.056 | 1.077 | 1.098 | 1.080 |  |  | 1.179 | 1.470 | - |
| L011138 | hs_vntr | 1.339 | 1.339 | 1.041 | 1.016 |  |  | 1.024 |  |  |  | 1.105 | - |
| L012272 | hs_vntr | 1.254 | 1.254 |  | 1.048 |  |  | 1.081 |  |  |  |  | - |
| L001909 | hs_vntr | 1.959 | 1.959 | 1.159 | 1.107 |  | 1.138 | 1.115 |  |  |  | 1.143 | - |
| L016870 | hs_other | 1.030 | 1.030 | 1.023 | 1.021 | 1.009 | 1.017 | 1.024 |  |  |  | 1.070 | - |
| L016124 | c_vntr_matched | 1.585 | 1.585 | 1.054 | 1.048 | 1.019 | 1.045 | 1.051 | 1.054 | 1.090 | 1.254 | 1.868 | - |
| L007687 | c_vntr_matched | 1.112 | 1.112 | 1.034 | 1.027 | 1.036 | 1.047 | 1.049 | 1.045 | 1.073 | 1.101 | 1.056 | - |
| L007172 | c_vntr_matched | 1.103 | 1.103 | 1.031 | 1.026 | 1.138 | 1.094 | 1.045 |  | 1.431 | 1.230 | 1.128 | - |
| L011430 | c_vntr_matched | 1.361 | 1.361 | 1.203 | 1.141 | 1.179 | 1.268 | 1.158 | 1.214 | 1.268 | 1.383 | 1.437 | - |
| L007065 | c_vntr_matched | 1.000 | 1.000 | 1.002 | 1.009 | 1.000 | 1.000 | 1.009 | 1.028 | 1.022 | 1.000 | 1.000 | - |
| L004021 | c_vntr_matched | 1.022 | 1.022 | 1.022 | 1.021 | 1.021 | 1.021 | 1.022 |  | 1.104 | 1.149 | 1.021 | - |
| L014111 | c_vntr_matched | 1.112 | 1.112 | 1.153 | 1.019 | 1.055 | 1.073 | 1.125 | 1.153 | 1.076 | 1.161 | 1.112 | - |
| L000506 | c_nontr_sv | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | - |
| L006556 | c_nontr_sv | 1.003 | 1.003 | 1.189 | 1.189 | 1.003 | 1.003 | 1.183 | 1.183 | 1.003 | 1.003 | 1.003 | - |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 1.335 (+0.000; 0/0) | 1.335 (n=9) | 1.081 (-0.275; 8/0) | 1.023 (-0.311; 9/0) | 1.045 (-0.274; 6/0) | 1.067 (-0.273; 7/0) | 1.056 (-0.304; 9/0) | 1.022 (-0.162; 1/0) | 1.111 (-0.177; 3/1) | 1.133 (-0.170; 5/0) | 1.266 (-0.007; 4/4) |  |
| non-VNTR hotspot (1) | 1.030 (+0.000; 0/0) | 1.030 (n=1) | 1.023 (-0.007; 1/0) | 1.021 (-0.009; 1/0) | 1.009 (-0.021; 1/0) | 1.017 (-0.012; 1/0) | 1.024 (-0.006; 1/0) |  |  |  | 1.070 (+0.040; 0/1) |  |
| matched VNTR controls (7) | 1.112 (+0.000; 0/0) | 1.112 (n=7) | 1.034 (-0.072; 4/2) | 1.026 (-0.084; 6/1) | 1.036 (-0.057; 5/1) | 1.047 (-0.039; 6/0) | 1.049 (-0.058; 4/3) | 1.054 (-0.066; 3/2) | 1.090 (-0.036; 4/3) | 1.161 (+0.022; 2/4) | 1.112 (+0.000; 2/3) |  |
| non-repeat SV controls (2) | 1.002 (+0.000; 0/0) | 1.002 (n=2) | 1.095 (+0.093; 0/1) | 1.095 (+0.093; 0/1) | 1.002 (+0.000; 0/0) | 1.002 (+0.000; 0/0) | 1.091 (+0.090; 0/1) | 1.091 (+0.090; 0/1) | 1.002 (+0.000; 0/0) | 1.002 (+0.000; 0/0) | 1.002 (+0.000; 0/0) |  |

#### MAPQ<5, truth-core reads (`truth_core_mapq_lt5_frac`, low is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 0.592 | 0.597 | 0.545 | 0.680 | 0.617 | 0.542 | 0.642 |  | 0.666 | 0.631 | 0.482 | 0.134 |
| L014297 | hs_vntr | 0.748 | 0.752 | 0.564 | 0.745 | 0.691 | 0.715 | 0.775 |  | 0.852 | 0.856 | 0.825 | 0.315 |
| L015415 | hs_vntr | 0.710 | 0.724 | 0.669 | 0.740 | 0.716 | 0.716 | 0.733 | 0.668 | 0.801 | 0.780 | 0.771 | 0.607 |
| L005990 | hs_vntr | 0.569 | 0.531 | 0.499 | 0.475 | 0.496 | 0.482 | 0.554 |  | 0.593 | 0.558 | 0.539 | 0.297 |
| L012184 | hs_vntr | 0.819 | 0.817 | 0.673 | 0.795 | 0.842 | 0.814 | 0.717 |  |  |  | 0.686 | 0.652 |
| L002013 | hs_vntr | 0.795 | 0.834 | 0.739 | 0.787 | 0.781 | 0.733 | 0.803 |  |  | 0.836 | 0.774 | 0.562 |
| L011138 | hs_vntr | 0.952 | 0.954 | 0.821 | 0.824 |  |  | 0.844 |  |  |  | 0.892 | 0.518 |
| L012272 | hs_vntr | 0.899 | 0.907 |  | 0.789 |  |  | 0.751 |  |  |  |  | 0.450 |
| L001909 | hs_vntr | 0.868 | 0.854 | 0.824 | 0.882 |  | 0.822 | 0.909 |  |  |  | 0.822 | 0.633 |
| L016870 | hs_other | 0.200 | 0.240 | 0.038 | 0.133 | 0.057 | 0.057 | 0.164 |  |  |  | 0.057 | 0.209 |
| L016124 | c_vntr_matched | 0.033 | 0.034 | 0.043 | 0.043 | 0.043 | 0.043 | 0.040 | 0.043 | 0.045 | 0.097 | 0.054 | 0.043 |
| L007687 | c_vntr_matched | 0.304 | 0.311 | 0.304 | 0.330 | 0.333 | 0.734 | 0.739 | 0.321 | 0.599 | 0.736 | 0.306 | 0.228 |
| L007172 | c_vntr_matched | 0.315 | 0.324 | 0.315 | 0.413 | 0.621 | 0.610 | 0.538 |  | 0.784 | 0.590 | 0.422 | 0.046 |
| L011430 | c_vntr_matched | 0.231 | 0.220 | 0.186 | 0.334 | 0.236 | 0.184 | 0.318 | 0.190 | 0.310 | 0.262 | 0.190 | 0.080 |
| L007065 | c_vntr_matched | 0.034 | 0.034 | 0.053 | 0.045 | 0.032 | 0.032 | 0.045 | 0.100 | 0.321 | 0.032 | 0.032 | 0.032 |
| L004021 | c_vntr_matched | 0.034 | 0.039 | 0.014 | 0.014 | 0.014 | 0.014 | 0.033 |  | 0.267 | 0.222 | 0.007 | 0.009 |
| L014111 | c_vntr_matched | 0.009 | 0.006 | 0.006 | 0.025 | 0.006 | 0.025 | 0.048 | 0.006 | 0.006 | 0.025 | 0.006 | 0.006 |
| L000506 | c_nontr_sv | 0.061 | 0.020 | 0.020 | 0.020 | 0.020 | 0.020 | 0.020 | 0.020 | 0.020 | 0.020 | 0.020 | 0.000 |
| L006556 | c_nontr_sv | 0.083 | 0.007 | 0.255 | 0.255 | 0.007 | 0.007 | 0.255 | 0.255 | 0.007 | 0.007 | 0.007 | 0.124 |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 0.795 (-0.003; 6/3) | 0.817 (n=9) | 0.671 (-0.075; 8/0) | 0.787 (-0.022; 6/3) | 0.704 (-0.021; 4/2) | 0.716 (-0.037; 7/0) | 0.751 (+0.009; 4/5) | 0.668 (-0.055; 1/0) | 0.734 (+0.074; 0/4) | 0.780 (+0.034; 0/5) | 0.773 (-0.046; 5/3) | 0.518 (-0.272; 9/0) |
| non-VNTR hotspot (1) | 0.200 (-0.040; 1/0) | 0.240 (n=1) | 0.038 (-0.202; 1/0) | 0.133 (-0.107; 1/0) | 0.057 (-0.183; 1/0) | 0.057 (-0.183; 1/0) | 0.164 (-0.077; 1/0) |  |  |  | 0.057 (-0.183; 1/0) | 0.209 (-0.031; 1/0) |
| matched VNTR controls (7) | 0.034 (-0.001; 4/2) | 0.039 (n=7) | 0.053 (-0.007; 4/2) | 0.045 (+0.019; 1/6) | 0.043 (+0.009; 2/4) | 0.043 (+0.009; 3/4) | 0.048 (+0.042; 1/6) | 0.100 (+0.009; 1/3) | 0.310 (+0.228; 0/6) | 0.222 (+0.062; 1/6) | 0.054 (-0.003; 4/2) | 0.043 (-0.030; 5/1) |
| non-repeat SV controls (2) | 0.072 (+0.058; 0/2) | 0.014 (n=2) | 0.138 (+0.124; 0/1) | 0.138 (+0.124; 0/1) | 0.014 (+0.000; 0/0) | 0.014 (+0.000; 0/0) | 0.138 (+0.124; 0/1) | 0.138 (+0.124; 0/1) | 0.014 (+0.000; 0/0) | 0.014 (+0.000; 0/0) | 0.014 (+0.000; 0/0) | 0.062 (+0.048; 1/1) |

#### truth-consistent core reads at MAPQ>=5 (`truth_consistent_core_mapq_ge5_frac`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 0.423 | 0.423 | 0.460 | 0.342 | 0.406 | 0.463 | 0.371 |  | 0.342 | 0.377 | 0.542 | 0.893 |
| L014297 | hs_vntr | 0.283 | 0.279 | 0.429 | 0.258 | 0.317 | 0.287 | 0.233 |  | 0.133 | 0.129 | 0.171 | 0.696 |
| L015415 | hs_vntr | 0.282 | 0.268 | 0.328 | 0.257 | 0.282 | 0.280 | 0.254 | 0.328 | 0.192 | 0.223 | 0.223 | 0.384 |
| L005990 | hs_vntr | 0.441 | 0.478 | 0.512 | 0.535 | 0.520 | 0.529 | 0.456 |  | 0.420 | 0.461 | 0.471 | 0.706 |
| L012184 | hs_vntr | 0.173 | 0.169 | 0.326 | 0.198 | 0.151 | 0.178 | 0.286 |  |  |  | 0.311 | 0.342 |
| L002013 | hs_vntr | 0.205 | 0.172 | 0.269 | 0.216 | 0.227 | 0.269 | 0.202 |  |  | 0.169 | 0.229 | 0.468 |
| L011138 | hs_vntr | 0.044 | 0.043 | 0.182 | 0.168 |  |  | 0.150 |  |  |  | 0.110 | 0.495 |
| L012272 | hs_vntr | 0.083 | 0.075 |  | 0.198 |  |  | 0.250 |  |  |  |  | 0.580 |
| L001909 | hs_vntr | 0.129 | 0.138 | 0.161 | 0.107 |  | 0.164 | 0.081 |  |  |  | 0.163 | 0.378 |
| L016870 | hs_other | 0.796 | 0.748 | 0.952 | 0.864 | 0.942 | 0.942 | 0.825 |  |  |  | 0.942 | 0.786 |
| L016124 | c_vntr_matched | 0.983 | 0.981 | 0.975 | 0.975 | 0.975 | 0.975 | 0.978 | 0.975 | 0.971 | 0.916 | 0.959 | 0.975 |
| L007687 | c_vntr_matched | 0.702 | 0.694 | 0.704 | 0.677 | 0.674 | 0.263 | 0.261 | 0.687 | 0.404 | 0.263 | 0.702 | 0.777 |
| L007172 | c_vntr_matched | 0.699 | 0.687 | 0.714 | 0.607 | 0.408 | 0.423 | 0.498 |  | 0.225 | 0.436 | 0.607 | 0.963 |
| L011430 | c_vntr_matched | 0.801 | 0.804 | 0.823 | 0.703 | 0.787 | 0.831 | 0.719 | 0.826 | 0.733 | 0.771 | 0.826 | 0.940 |
| L007065 | c_vntr_matched | 0.974 | 0.975 | 0.953 | 0.961 | 0.975 | 0.975 | 0.961 | 0.908 | 0.685 | 0.975 | 0.975 | 0.976 |
| L004021 | c_vntr_matched | 0.970 | 0.962 | 0.987 | 0.987 | 0.987 | 0.987 | 0.968 |  | 0.742 | 0.783 | 0.994 | 0.992 |
| L014111 | c_vntr_matched | 0.994 | 0.995 | 0.997 | 0.976 | 0.997 | 0.977 | 0.953 | 0.997 | 0.997 | 0.977 | 0.997 | 0.997 |
| L000506 | c_nontr_sv | 0.939 | 0.980 | 0.980 | 0.980 | 0.980 | 0.980 | 0.980 | 0.980 | 0.980 | 0.980 | 0.980 | 1.000 |
| L006556 | c_nontr_sv | 0.910 | 0.992 | 0.853 | 0.853 | 0.992 | 0.992 | 0.853 | 0.853 | 0.992 | 0.992 | 0.992 | 0.885 |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 0.205 (+0.004; 6/2) | 0.172 (n=9) | 0.327 (+0.078; 8/0) | 0.216 (+0.029; 5/4) | 0.300 (+0.026; 4/2) | 0.280 (+0.026; 7/0) | 0.250 (-0.014; 4/5) | 0.328 (+0.059; 1/0) | 0.267 (-0.079; 0/4) | 0.223 (-0.045; 0/5) | 0.226 (+0.041; 5/3) | 0.495 (+0.296; 9/0) |
| non-VNTR hotspot (1) | 0.796 (+0.048; 1/0) | 0.748 (n=1) | 0.952 (+0.204; 1/0) | 0.864 (+0.116; 1/0) | 0.942 (+0.194; 1/0) | 0.942 (+0.194; 1/0) | 0.825 (+0.078; 1/0) |  |  |  | 0.942 (+0.194; 1/0) | 0.786 (+0.039; 1/0) |
| matched VNTR controls (7) | 0.970 (+0.002; 4/3) | 0.962 (n=7) | 0.953 (+0.010; 5/2) | 0.961 (-0.018; 1/6) | 0.975 (-0.006; 2/4) | 0.975 (-0.006; 2/4) | 0.953 (-0.042; 1/6) | 0.908 (-0.006; 2/3) | 0.733 (-0.220; 1/6) | 0.783 (-0.065; 0/6) | 0.959 (+0.002; 4/2) | 0.975 (+0.030; 6/1) |
| non-repeat SV controls (2) | 0.924 (-0.061; 0/2) | 0.986 (n=2) | 0.916 (-0.070; 0/1) | 0.916 (-0.070; 0/1) | 0.986 (+0.000; 0/0) | 0.986 (+0.000; 0/0) | 0.916 (-0.070; 0/1) | 0.916 (-0.070; 0/1) | 0.986 (+0.000; 0/0) | 0.986 (+0.000; 0/0) | 0.986 (+0.000; 0/0) | 0.943 (-0.043; 1/1) |

#### uniquely-placeable core reads placed right at MAPQ>=5 (`core_placed_confident_correct_frac`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 0.399 | 0.403 | 0.465 | 0.368 | 0.424 | 0.469 | 0.387 |  | 0.358 | 0.428 | 0.588 | - |
| L014297 | hs_vntr | 0.418 | 0.392 | 0.443 | 0.342 | 0.430 | 0.367 | 0.279 |  | 0.089 | 0.101 | 0.114 | - |
| L015415 | hs_vntr | 0.648 | 0.648 | 0.915 | 0.803 | 0.803 | 0.789 | 0.845 | 0.915 | 0.549 | 0.549 | 0.592 | - |
| L005990 | hs_vntr | 0.375 | 0.375 | 0.500 | 0.583 | 0.375 | 0.500 | 0.292 |  | 0.292 | 0.375 | 0.500 | - |
| L012184 | hs_vntr | 0.125 | 0.125 | 0.750 | 0.625 | 0.542 | 0.542 | 0.667 |  |  |  | 0.583 | - |
| L002013 | hs_vntr | - | - | - | - | - | - | - |  |  | - | - | - |
| L011138 | hs_vntr | 0.089 | 0.089 | 0.446 | 0.287 |  |  | 0.268 |  |  |  | 0.261 | - |
| L012272 | hs_vntr | - | - |  | - |  |  | - |  |  |  |  | - |
| L001909 | hs_vntr | - | - | - | - |  | - | - |  |  |  | - | - |
| L016870 | hs_other | 0.667 | 0.625 | 0.931 | 0.778 | 0.903 | 0.931 | 0.764 |  |  |  | 0.917 | - |
| L016124 | c_vntr_matched | 0.986 | 0.986 | 0.990 | 0.988 | 0.988 | 0.988 | 0.986 | 0.990 | 0.984 | 0.922 | 0.960 | - |
| L007687 | c_vntr_matched | 0.864 | 0.864 | 0.911 | 0.834 | 0.834 | 0.396 | 0.421 | 0.843 | 0.596 | 0.404 | 0.860 | - |
| L007172 | c_vntr_matched | 0.723 | 0.718 | 0.754 | 0.607 | 0.429 | 0.446 | 0.532 |  | 0.218 | 0.459 | 0.652 | - |
| L011430 | c_vntr_matched | 0.865 | 0.869 | 0.892 | 0.788 | 0.838 | 0.873 | 0.803 | 0.892 | 0.819 | 0.857 | 0.873 | - |
| L007065 | c_vntr_matched | 0.986 | 0.989 | 0.968 | 0.973 | 0.988 | 0.988 | 0.971 | 0.913 | 0.701 | 0.988 | 0.988 | - |
| L004021 | c_vntr_matched | 0.972 | 0.959 | 0.979 | 0.980 | 0.980 | 0.980 | 0.966 |  | 0.749 | 0.784 | 0.994 | - |
| L014111 | c_vntr_matched | 0.995 | 0.995 | 0.992 | 0.969 | 0.992 | 0.972 | 0.944 | 0.992 | 0.993 | 0.972 | 0.992 | - |
| L000506 | c_nontr_sv | 0.835 | 0.856 | 0.856 | 0.856 | 0.856 | 0.856 | 0.856 | 0.856 | 0.856 | 0.856 | 0.856 | - |
| L006556 | c_nontr_sv | 0.882 | 0.987 | 0.962 | 0.962 | 0.987 | 0.987 | 0.962 | 0.962 | 0.987 | 0.987 | 0.987 | - |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 0.387 (+0.000; 1/1) | 0.384 (n=6) | 0.483 (+0.196; 6/0) | 0.476 (+0.176; 4/2) | 0.430 (+0.038; 4/0) | 0.500 (+0.125; 4/1) | 0.339 (+0.081; 3/3) | 0.915 (+0.268; 1/0) | 0.325 (-0.091; 0/4) | 0.401 (-0.049; 1/2) | 0.542 (+0.148; 4/2) |  |
| non-VNTR hotspot (1) | 0.667 (+0.042; 1/0) | 0.625 (n=1) | 0.931 (+0.306; 1/0) | 0.778 (+0.153; 1/0) | 0.903 (+0.278; 1/0) | 0.931 (+0.306; 1/0) | 0.764 (+0.139; 1/0) |  |  |  | 0.917 (+0.292; 1/0) |  |
| matched VNTR controls (7) | 0.972 (+0.000; 2/2) | 0.959 (n=7) | 0.968 (+0.020; 5/2) | 0.969 (-0.026; 2/5) | 0.980 (-0.003; 2/5) | 0.972 (-0.001; 3/4) | 0.944 (-0.051; 1/5) | 0.913 (-0.003; 2/3) | 0.749 (-0.210; 0/7) | 0.857 (-0.064; 0/7) | 0.960 (-0.003; 2/5) |  |
| non-repeat SV controls (2) | 0.858 (-0.063; 0/2) | 0.921 (n=2) | 0.909 (-0.013; 0/1) | 0.909 (-0.013; 0/1) | 0.921 (+0.000; 0/0) | 0.921 (+0.000; 0/0) | 0.909 (-0.013; 0/1) | 0.909 (-0.013; 0/1) | 0.921 (+0.000; 0/0) | 0.921 (+0.000; 0/0) | 0.921 (+0.000; 0/0) |  |

#### core reads at MAPQ>=5 placed wrong (`core_confident_wrong_frac`, low is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 0.191 | 0.179 | 0.119 | 0.158 | 0.129 | 0.123 | 0.174 |  | 0.156 | 0.093 | 0.088 | - |
| L014297 | hs_vntr | 0.083 | 0.139 | 0.146 | 0.100 | 0.128 | 0.147 | 0.083 |  | 0.364 | 0.333 | 0.438 | - |
| L015415 | hs_vntr | 0.061 | 0.061 | 0.030 | 0.034 | 0.034 | 0.051 | 0.016 | 0.030 | 0.071 | 0.114 | 0.067 | - |
| L005990 | hs_vntr | 0.000 | 0.000 | 0.250 | 0.125 | 0.438 | 0.200 | 0.000 |  | 0.000 | 0.182 | 0.077 | - |
| L012184 | hs_vntr | 0.250 | 0.250 | 0.143 | 0.062 | 0.071 | 0.133 | 0.111 |  |  |  | 0.176 | - |
| L002013 | hs_vntr | - | - | - | - | - | - | - |  |  | - | - | - |
| L011138 | hs_vntr | 0.176 | 0.125 | 0.091 | 0.210 |  |  | 0.250 |  |  |  | 0.068 | - |
| L012272 | hs_vntr | - | - |  | - |  |  | - |  |  |  |  | - |
| L001909 | hs_vntr | - | - | - | - |  | - | - |  |  |  | - | - |
| L016870 | hs_other | 0.094 | 0.100 | 0.043 | 0.067 | 0.044 | 0.029 | 0.035 |  |  |  | 0.057 | - |
| L016124 | c_vntr_matched | 0.012 | 0.012 | 0.010 | 0.012 | 0.012 | 0.012 | 0.014 | 0.010 | 0.012 | 0.017 | 0.014 | - |
| L007687 | c_vntr_matched | 0.024 | 0.024 | 0.023 | 0.049 | 0.062 | 0.070 | 0.010 | 0.066 | 0.021 | 0.031 | 0.060 | - |
| L007172 | c_vntr_matched | 0.021 | 0.019 | 0.024 | 0.050 | 0.045 | 0.044 | 0.053 |  | 0.103 | 0.038 | 0.017 | - |
| L011430 | c_vntr_matched | 0.055 | 0.055 | 0.049 | 0.042 | 0.057 | 0.062 | 0.037 | 0.049 | 0.045 | 0.047 | 0.054 | - |
| L007065 | c_vntr_matched | 0.012 | 0.011 | 0.008 | 0.012 | 0.012 | 0.012 | 0.014 | 0.018 | 0.017 | 0.012 | 0.012 | - |
| L004021 | c_vntr_matched | 0.010 | 0.009 | 0.016 | 0.015 | 0.015 | 0.015 | 0.009 |  | 0.014 | 0.011 | 0.006 | - |
| L014111 | c_vntr_matched | 0.003 | 0.003 | 0.008 | 0.008 | 0.008 | 0.008 | 0.012 | 0.008 | 0.007 | 0.008 | 0.008 | - |
| L000506 | c_nontr_sv | 0.110 | 0.126 | 0.126 | 0.126 | 0.126 | 0.126 | 0.126 | 0.126 | 0.126 | 0.126 | 0.126 | - |
| L006556 | c_nontr_sv | 0.015 | 0.013 | 0.026 | 0.026 | 0.013 | 0.013 | 0.026 | 0.026 | 0.013 | 0.013 | 0.013 | - |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 0.130 (+0.000; 1/2) | 0.132 (n=6) | 0.131 (-0.033; 4/2) | 0.113 (-0.024; 4/2) | 0.128 (-0.027; 4/1) | 0.133 (-0.010; 3/2) | 0.097 (-0.025; 4/1) | 0.030 (-0.031; 1/0) | 0.113 (+0.005; 1/2) | 0.148 (+0.117; 1/3) | 0.082 (-0.026; 3/3) |  |
| non-VNTR hotspot (1) | 0.094 (-0.006; 1/0) | 0.100 (n=1) | 0.043 (-0.057; 1/0) | 0.067 (-0.033; 1/0) | 0.044 (-0.056; 1/0) | 0.029 (-0.071; 1/0) | 0.035 (-0.065; 1/0) |  |  |  | 0.057 (-0.043; 1/0) |  |
| matched VNTR controls (7) | 0.012 (+0.000; 0/4) | 0.012 (n=7) | 0.016 (-0.001; 4/3) | 0.015 (+0.005; 1/5) | 0.015 (+0.005; 0/6) | 0.015 (+0.006; 0/6) | 0.014 (+0.002; 3/4) | 0.018 (+0.005; 2/3) | 0.017 (+0.003; 2/5) | 0.017 (+0.005; 1/6) | 0.014 (+0.001; 3/4) |  |
| non-repeat SV controls (2) | 0.062 (-0.007; 1/1) | 0.069 (n=2) | 0.076 (+0.006; 0/1) | 0.076 (+0.006; 0/1) | 0.069 (+0.000; 0/0) | 0.069 (+0.000; 0/0) | 0.076 (+0.006; 0/1) | 0.076 (+0.006; 0/1) | 0.069 (+0.000; 0/0) | 0.069 (+0.000; 0/0) | 0.069 (+0.000; 0/0) |  |

#### uniquely-placeable core reads placed right (any MAPQ) (`core_placed_correctly_frac`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 0.598 | 0.604 | 0.629 | 0.525 | 0.604 | 0.632 | 0.560 |  | 0.519 | 0.663 | 0.852 | - |
| L014297 | hs_vntr | 0.570 | 0.557 | 0.722 | 0.519 | 0.582 | 0.595 | 0.494 |  | 0.468 | 0.304 | 0.380 | - |
| L015415 | hs_vntr | 0.803 | 0.831 | 0.944 | 0.901 | 0.915 | 0.915 | 0.915 | 0.915 | 0.803 | 0.704 | 0.761 | - |
| L005990 | hs_vntr | 0.458 | 0.500 | 0.667 | 0.708 | 0.500 | 0.667 | 0.833 |  | 0.667 | 0.667 | 0.750 | - |
| L012184 | hs_vntr | 0.375 | 0.417 | 0.792 | 0.625 | 0.667 | 0.750 | 0.750 |  |  |  | 0.708 | - |
| L002013 | hs_vntr | - | - | - | - | - | - | - |  |  | - | - | - |
| L011138 | hs_vntr | 0.452 | 0.478 | 0.612 | 0.490 |  |  | 0.427 |  |  |  | 0.484 | - |
| L012272 | hs_vntr | - | - |  | - |  |  | - |  |  |  |  | - |
| L001909 | hs_vntr | - | - | - | - |  | - | - |  |  |  | - | - |
| L016870 | hs_other | 0.917 | 0.917 | 0.958 | 0.917 | 0.931 | 0.944 | 0.958 |  |  |  | 0.931 | - |
| L016124 | c_vntr_matched | 0.988 | 0.988 | 0.990 | 0.988 | 0.988 | 0.988 | 0.986 | 0.990 | 0.988 | 0.984 | 0.986 | - |
| L007687 | c_vntr_matched | 0.953 | 0.957 | 0.962 | 0.932 | 0.928 | 0.706 | 0.672 | 0.919 | 0.847 | 0.728 | 0.906 | - |
| L007172 | c_vntr_matched | 0.916 | 0.902 | 0.918 | 0.818 | 0.804 | 0.825 | 0.766 |  | 0.682 | 0.752 | 0.861 | - |
| L011430 | c_vntr_matched | 0.942 | 0.938 | 0.946 | 0.903 | 0.930 | 0.930 | 0.911 | 0.950 | 0.927 | 0.938 | 0.927 | - |
| L007065 | c_vntr_matched | 0.988 | 0.989 | 0.988 | 0.988 | 0.988 | 0.988 | 0.986 | 0.959 | 0.883 | 0.988 | 0.988 | - |
| L004021 | c_vntr_matched | 0.989 | 0.991 | 0.983 | 0.984 | 0.984 | 0.984 | 0.991 |  | 0.871 | 0.925 | 0.994 | - |
| L014111 | c_vntr_matched | 0.997 | 0.995 | 0.992 | 0.989 | 0.992 | 0.989 | 0.972 | 0.992 | 0.993 | 0.989 | 0.992 | - |
| L000506 | c_nontr_sv | 0.876 | 0.876 | 0.866 | 0.866 | 0.866 | 0.866 | 0.866 | 0.866 | 0.866 | 0.866 | 0.866 | - |
| L006556 | c_nontr_sv | 0.987 | 0.987 | 0.975 | 0.975 | 0.987 | 0.987 | 0.962 | 0.962 | 0.987 | 0.987 | 0.987 | - |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 0.514 (-0.027; 1/5) | 0.528 (n=6) | 0.694 (+0.149; 6/0) | 0.575 (+0.042; 4/2) | 0.604 (+0.025; 3/0) | 0.667 (+0.085; 5/0) | 0.655 (+0.020; 3/3) | 0.915 (+0.085; 1/0) | 0.593 (-0.057; 1/3) | 0.665 (-0.034; 2/2) | 0.729 (+0.127; 4/2) |  |
| non-VNTR hotspot (1) | 0.917 (+0.000; 0/0) | 0.917 (n=1) | 0.958 (+0.042; 1/0) | 0.917 (+0.000; 0/0) | 0.931 (+0.014; 1/0) | 0.944 (+0.028; 1/0) | 0.958 (+0.042; 1/0) |  |  |  | 0.931 (+0.014; 1/0) |  |
| matched VNTR controls (7) | 0.988 (+0.000; 3/3) | 0.988 (n=7) | 0.983 (+0.002; 4/3) | 0.984 (-0.007; 0/6) | 0.984 (-0.007; 0/6) | 0.984 (-0.007; 0/6) | 0.972 (-0.023; 0/6) | 0.959 (-0.003; 2/3) | 0.883 (-0.106; 0/6) | 0.938 (-0.007; 0/6) | 0.986 (-0.003; 1/6) |  |
| non-repeat SV controls (2) | 0.932 (-0.000; 0/1) | 0.932 (n=2) | 0.920 (-0.011; 0/2) | 0.920 (-0.011; 0/2) | 0.927 (-0.005; 0/1) | 0.927 (-0.005; 0/1) | 0.914 (-0.018; 0/2) | 0.914 (-0.018; 0/2) | 0.927 (-0.005; 0/1) | 0.927 (-0.005; 0/1) | 0.927 (-0.005; 0/1) |  |

#### core reads on a single haplotype path (`core_compatible_frac`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 0.988 | 0.983 | 0.969 | 0.976 | 0.967 | 0.979 | 0.967 |  | 0.971 | 0.962 | 0.976 | - |
| L014297 | hs_vntr | 0.983 | 0.983 | 0.990 | 0.993 | 0.983 | 0.980 | 0.980 |  | 0.980 | 0.990 | 0.987 | - |
| L015415 | hs_vntr | 0.988 | 0.988 | 0.998 | 0.985 | 0.995 | 0.990 | 0.988 | 0.990 | 0.990 | 0.988 | 0.988 | - |
| L005990 | hs_vntr | 0.964 | 0.962 | 0.979 | 0.965 | 0.976 | 0.974 | 0.949 |  | 0.942 | 0.953 | 0.983 | - |
| L012184 | hs_vntr | 0.942 | 0.945 | 0.937 | 0.945 | 0.916 | 0.905 | 0.932 |  |  |  | 0.940 | - |
| L002013 | hs_vntr | 0.989 | 0.986 | 0.983 | 0.984 | 0.983 | 0.980 | 0.984 |  |  | 0.986 | 0.987 | - |
| L011138 | hs_vntr | 0.990 | 0.991 | 0.995 | 0.997 |  |  | 0.996 |  |  |  | 0.996 | - |
| L012272 | hs_vntr | 0.970 | 0.972 |  | 0.966 |  |  | 0.966 |  |  |  |  | - |
| L001909 | hs_vntr | 0.976 | 0.979 | 0.986 | 0.996 |  | 0.995 | 0.998 |  |  |  | 0.989 | - |
| L016870 | hs_other | 0.990 | 0.989 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |  |  |  | 1.000 | - |
| L016124 | c_vntr_matched | 0.997 | 0.997 | 0.994 | 0.993 | 0.993 | 0.994 | 0.993 | 0.994 | 0.993 | 0.993 | 0.999 | - |
| L007687 | c_vntr_matched | 0.998 | 0.998 | 1.000 | 1.000 | 0.998 | 0.998 | 0.998 | 0.998 | 0.998 | 0.998 | 1.000 | - |
| L007172 | c_vntr_matched | 0.989 | 0.990 | 0.989 | 0.979 | 0.982 | 0.987 | 0.984 |  | 0.980 | 0.989 | 0.987 | - |
| L011430 | c_vntr_matched | 0.972 | 0.974 | 0.970 | 0.974 | 0.976 | 0.978 | 0.968 | 0.972 | 0.962 | 0.966 | 0.974 | - |
| L007065 | c_vntr_matched | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 | 0.999 | 1.000 | 1.000 | - |
| L004021 | c_vntr_matched | 0.998 | 0.999 | 0.999 | 0.999 | 0.999 | 0.999 | 0.997 |  | 0.993 | 0.993 | 0.998 | - |
| L014111 | c_vntr_matched | 0.999 | 1.000 | 1.000 | 0.991 | 1.000 | 0.994 | 0.993 | 1.000 | 0.999 | 0.991 | 1.000 | - |
| L000506 | c_nontr_sv | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | - |
| L006556 | c_nontr_sv | 1.000 | 1.000 | 0.972 | 0.972 | 1.000 | 1.000 | 0.972 | 0.972 | 1.000 | 1.000 | 1.000 | - |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 0.983 (+0.000; 4/4) | 0.983 (n=9) | 0.984 (+0.005; 5/3) | 0.984 (+0.000; 5/4) | 0.979 (-0.002; 2/3) | 0.980 (-0.003; 3/4) | 0.980 (-0.003; 2/6) | 0.990 (+0.002; 1/0) | 0.975 (-0.008; 1/3) | 0.986 (-0.000; 1/4) | 0.987 (+0.002; 5/2) |  |
| non-VNTR hotspot (1) | 0.990 (+0.000; 1/0) | 0.989 (n=1) | 1.000 (+0.011; 1/0) | 1.000 (+0.011; 1/0) | 1.000 (+0.011; 1/0) | 1.000 (+0.011; 1/0) | 1.000 (+0.011; 1/0) |  |  |  | 1.000 (+0.011; 1/0) |  |
| matched VNTR controls (7) | 0.998 (-0.001; 0/4) | 0.998 (n=7) | 0.999 (+0.000; 1/3) | 0.993 (+0.000; 1/3) | 0.998 (+0.000; 1/2) | 0.994 (+0.000; 1/3) | 0.993 (-0.004; 0/5) | 0.998 (-0.001; 0/3) | 0.993 (-0.004; 0/6) | 0.993 (-0.004; 0/5) | 0.999 (+0.000; 2/2) |  |
| non-repeat SV controls (2) | 1.000 (+0.000; 0/0) | 1.000 (n=2) | 0.986 (-0.014; 0/1) | 0.986 (-0.014; 0/1) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 0.986 (-0.014; 0/1) | 0.986 (-0.014; 0/1) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) |  |

#### haplotypes to cover 95% of compatible core reads (`core_haps_cover95`, low is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 7 | 8 | 9 | 9 | 8 | 9 | 9 |  | 10 | 9 | 7 | - |
| L014297 | hs_vntr | 9 | 9 | 9 | 7 | 8 | 8 | 8 |  | 10 | 10 | 10 | - |
| L015415 | hs_vntr | 7 | 7 | 6 | 7 | 6 | 6 | 8 | 6 | 8 | 8 | 7 | - |
| L005990 | hs_vntr | 13 | 13 | 11 | 11 | 12 | 12 | 12 |  | 16 | 13 | 12 | - |
| L012184 | hs_vntr | 17 | 18 | 15 | 17 | 17 | 16 | 18 |  |  |  | 17 | - |
| L002013 | hs_vntr | 19 | 19 | 19 | 18 | 19 | 19 | 20 |  |  | 20 | 19 | - |
| L011138 | hs_vntr | 11 | 11 | 7 | 8 |  |  | 12 |  |  |  | 12 | - |
| L012272 | hs_vntr | 18 | 18 |  | 16 |  |  | 18 |  |  |  |  | - |
| L001909 | hs_vntr | 28 | 28 | 27 | 28 |  | 28 | 28 |  |  |  | 28 | - |
| L016870 | hs_other | 4 | 4 | 2 | 3 | 3 | 2 | 3 |  |  |  | 3 | - |
| L016124 | c_vntr_matched | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | - |
| L007687 | c_vntr_matched | 2 | 3 | 2 | 3 | 3 | 3 | 3 | 2 | 4 | 3 | 3 | - |
| L007172 | c_vntr_matched | 3 | 3 | 3 | 4 | 5 | 4 | 5 |  | 5 | 5 | 4 | - |
| L011430 | c_vntr_matched | 5 | 4 | 5 | 5 | 5 | 5 | 5 | 4 | 5 | 5 | 5 | - |
| L007065 | c_vntr_matched | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 3 | 2 | 2 | - |
| L004021 | c_vntr_matched | 2 | 2 | 2 | 2 | 2 | 2 | 2 |  | 3 | 3 | 2 | - |
| L014111 | c_vntr_matched | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | - |
| L000506 | c_nontr_sv | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | - |
| L006556 | c_nontr_sv | 2 | 2 | 3 | 3 | 2 | 2 | 3 | 3 | 2 | 2 | 2 | - |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 13.000 (+0.000; 2/0) | 13.000 (n=9) | 10.000 (-1.000; 5/1) | 11.000 (-1.000; 6/1) | 10.000 (-1.000; 4/0) | 12.000 (-1.000; 4/1) | 12.000 (+0.000; 2/4) | 6.000 (-1.000; 1/0) | 10.000 (+1.500; 0/4) | 10.000 (+1.000; 0/4) | 12.000 (+0.000; 3/2) |  |
| non-VNTR hotspot (1) | 4.000 (+0.000; 0/0) | 4.000 (n=1) | 2.000 (-2.000; 1/0) | 3.000 (-1.000; 1/0) | 3.000 (-1.000; 1/0) | 2.000 (-2.000; 1/0) | 3.000 (-1.000; 1/0) |  |  |  | 3.000 (-1.000; 1/0) |  |
| matched VNTR controls (7) | 2.000 (+0.000; 1/1) | 2.000 (n=7) | 2.000 (+0.000; 1/1) | 2.000 (+0.000; 0/2) | 2.000 (+0.000; 0/2) | 2.000 (+0.000; 0/2) | 2.000 (+0.000; 0/2) | 2.000 (+0.000; 1/0) | 3.000 (+1.000; 0/5) | 3.000 (+0.000; 0/3) | 2.000 (+0.000; 0/2) |  |
| non-repeat SV controls (2) | 3.000 (+0.000; 0/0) | 3.000 (n=2) | 3.500 (+0.500; 0/1) | 3.500 (+0.500; 0/1) | 3.000 (+0.000; 0/0) | 3.000 (+0.000; 0/0) | 3.500 (+0.500; 0/1) | 3.500 (+0.500; 0/1) | 3.000 (+0.000; 0/0) | 3.000 (+0.000; 0/0) | 3.000 (+0.000; 0/0) |  |

#### core reads explained by the best path pair (`core_best_pair_frac`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 0.694 | 0.691 | 0.667 | 0.632 | 0.606 | 0.694 | 0.628 |  | 0.576 | 0.598 | 0.717 | - |
| L014297 | hs_vntr | 0.571 | 0.567 | 0.705 | 0.607 | 0.627 | 0.614 | 0.527 |  | 0.399 | 0.379 | 0.444 | - |
| L015415 | hs_vntr | 0.581 | 0.584 | 0.658 | 0.593 | 0.643 | 0.666 | 0.602 | 0.634 | 0.553 | 0.529 | 0.563 | - |
| L005990 | hs_vntr | 0.576 | 0.579 | 0.612 | 0.582 | 0.562 | 0.568 | 0.521 |  | 0.456 | 0.517 | 0.602 | - |
| L012184 | hs_vntr | 0.447 | 0.443 | 0.561 | 0.419 | 0.387 | 0.439 | 0.367 |  |  |  | 0.407 | - |
| L002013 | hs_vntr | 0.480 | 0.515 | 0.480 | 0.451 | 0.451 | 0.483 | 0.469 |  |  | 0.398 | 0.442 | - |
| L011138 | hs_vntr | 0.472 | 0.447 | 0.717 | 0.617 |  |  | 0.466 |  |  |  | 0.507 | - |
| L012272 | hs_vntr | 0.300 | 0.283 |  | 0.389 |  |  | 0.360 |  |  |  |  | - |
| L001909 | hs_vntr | 0.230 | 0.223 | 0.227 | 0.206 |  | 0.230 | 0.222 |  |  |  | 0.230 | - |
| L016870 | hs_other | 0.895 | 0.904 | 0.979 | 0.926 | 0.947 | 0.979 | 0.936 |  |  |  | 0.926 | - |
| L016124 | c_vntr_matched | 0.996 | 0.994 | 0.989 | 0.986 | 0.986 | 0.987 | 0.986 | 0.989 | 0.987 | 0.986 | 0.992 | - |
| L007687 | c_vntr_matched | 0.948 | 0.941 | 0.967 | 0.948 | 0.926 | 0.931 | 0.850 | 0.953 | 0.907 | 0.919 | 0.936 | - |
| L007172 | c_vntr_matched | 0.902 | 0.897 | 0.899 | 0.838 | 0.812 | 0.868 | 0.809 |  | 0.668 | 0.837 | 0.868 | - |
| L011430 | c_vntr_matched | 0.817 | 0.823 | 0.841 | 0.801 | 0.819 | 0.839 | 0.793 | 0.855 | 0.786 | 0.806 | 0.839 | - |
| L007065 | c_vntr_matched | 0.984 | 0.984 | 0.986 | 0.986 | 0.986 | 0.986 | 0.984 | 0.974 | 0.947 | 0.986 | 0.986 | - |
| L004021 | c_vntr_matched | 0.986 | 0.988 | 0.979 | 0.979 | 0.979 | 0.979 | 0.984 |  | 0.883 | 0.909 | 0.987 | - |
| L014111 | c_vntr_matched | 0.999 | 0.999 | 0.997 | 0.990 | 0.997 | 0.991 | 0.973 | 0.997 | 0.996 | 0.988 | 0.997 | - |
| L000506 | c_nontr_sv | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | - |
| L006556 | c_nontr_sv | 0.991 | 0.992 | 0.923 | 0.923 | 0.992 | 0.992 | 0.897 | 0.897 | 0.992 | 0.992 | 0.992 | - |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 0.480 (+0.003; 6/3) | 0.515 (n=9) | 0.635 (+0.054; 6/2) | 0.582 (+0.003; 5/4) | 0.584 (-0.036; 2/4) | 0.568 (+0.003; 4/3) | 0.469 (-0.040; 3/6) | 0.634 (+0.050; 1/0) | 0.504 (-0.119; 0/4) | 0.517 (-0.093; 0/5) | 0.476 (-0.007; 4/4) |  |
| non-VNTR hotspot (1) | 0.895 (-0.010; 0/1) | 0.904 (n=1) | 0.979 (+0.074; 1/0) | 0.926 (+0.022; 1/0) | 0.947 (+0.043; 1/0) | 0.979 (+0.075; 1/0) | 0.936 (+0.032; 1/0) |  |  |  | 0.926 (+0.022; 1/0) |  |
| matched VNTR controls (7) | 0.984 (+0.000; 3/2) | 0.984 (n=7) | 0.979 (+0.001; 4/3) | 0.979 (-0.009; 2/5) | 0.979 (-0.008; 1/6) | 0.979 (-0.008; 2/5) | 0.973 (-0.025; 0/6) | 0.974 (-0.002; 2/3) | 0.907 (-0.037; 0/7) | 0.919 (-0.018; 1/6) | 0.986 (-0.002; 2/5) |  |
| non-repeat SV controls (2) | 0.939 (-0.000; 0/1) | 0.940 (n=2) | 0.905 (-0.034; 0/1) | 0.905 (-0.034; 0/1) | 0.940 (+0.000; 0/0) | 0.940 (+0.000; 0/0) | 0.892 (-0.048; 0/1) | 0.892 (-0.048; 0/1) | 0.940 (+0.000; 0/0) | 0.940 (+0.000; 0/0) | 0.940 (+0.000; 0/0) |  |

#### core reads explained by the truth-closest pair (`core_truth_pair_frac`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 0.487 | 0.484 | 0.481 | 0.442 | 0.466 | 0.546 | 0.451 |  | 0.409 | 0.508 | 0.658 | - |
| L014297 | hs_vntr | 0.376 | 0.379 | 0.530 | 0.440 | 0.436 | 0.436 | 0.376 |  | 0.269 | 0.205 | 0.253 | - |
| L015415 | hs_vntr | 0.453 | 0.450 | 0.612 | 0.522 | 0.584 | 0.598 | 0.476 | 0.527 | 0.406 | 0.434 | 0.456 | - |
| L005990 | hs_vntr | 0.441 | 0.447 | 0.504 | 0.490 | 0.472 | 0.476 | 0.420 |  | 0.391 | 0.406 | 0.523 | - |
| L012184 | hs_vntr | 0.382 | 0.382 | 0.495 | 0.355 | 0.367 | 0.399 | 0.321 |  |  |  | 0.353 | - |
| L002013 | hs_vntr | 0.150 | 0.152 | 0.186 | 0.181 | 0.165 | 0.173 | 0.172 |  |  | 0.112 | 0.186 | - |
| L011138 | hs_vntr | 0.324 | 0.312 | 0.482 | 0.395 |  |  | 0.312 |  |  |  | 0.396 | - |
| L012272 | hs_vntr | 0.247 | 0.235 |  | 0.309 |  |  | 0.294 |  |  |  |  | - |
| L001909 | hs_vntr | 0.138 | 0.145 | 0.152 | 0.139 |  | 0.161 | 0.151 |  |  |  | 0.164 | - |
| L016870 | hs_other | 0.895 | 0.894 | 0.947 | 0.916 | 0.916 | 0.937 | 0.936 |  |  |  | 0.884 | - |
| L016124 | c_vntr_matched | 0.992 | 0.989 | 0.986 | 0.983 | 0.983 | 0.984 | 0.983 | 0.986 | 0.984 | 0.983 | 0.989 | - |
| L007687 | c_vntr_matched | 0.931 | 0.931 | 0.955 | 0.903 | 0.900 | 0.696 | 0.644 | 0.922 | 0.803 | 0.677 | 0.912 | - |
| L007172 | c_vntr_matched | 0.855 | 0.848 | 0.868 | 0.786 | 0.770 | 0.779 | 0.688 |  | 0.650 | 0.693 | 0.789 | - |
| L011430 | c_vntr_matched | 0.817 | 0.823 | 0.841 | 0.787 | 0.813 | 0.833 | 0.773 | 0.851 | 0.786 | 0.806 | 0.839 | - |
| L007065 | c_vntr_matched | 0.983 | 0.984 | 0.984 | 0.984 | 0.984 | 0.984 | 0.983 | 0.950 | 0.870 | 0.984 | 0.984 | - |
| L004021 | c_vntr_matched | 0.985 | 0.986 | 0.977 | 0.978 | 0.978 | 0.978 | 0.982 |  | 0.857 | 0.908 | 0.986 | - |
| L014111 | c_vntr_matched | 0.999 | 0.999 | 0.997 | 0.990 | 0.997 | 0.991 | 0.973 | 0.997 | 0.996 | 0.988 | 0.997 | - |
| L000506 | c_nontr_sv | 0.878 | 0.878 | 0.867 | 0.867 | 0.867 | 0.867 | 0.867 | 0.867 | 0.867 | 0.867 | 0.867 | - |
| L006556 | c_nontr_sv | 0.991 | 0.992 | 0.923 | 0.923 | 0.992 | 0.992 | 0.897 | 0.897 | 0.992 | 0.992 | 0.992 | - |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 0.376 (+0.000; 4/4) | 0.379 (n=9) | 0.489 (+0.085; 7/1) | 0.395 (+0.042; 6/3) | 0.451 (+0.019; 4/2) | 0.436 (+0.029; 7/0) | 0.321 (-0.001; 4/5) | 0.527 (+0.077; 1/0) | 0.398 (-0.066; 0/4) | 0.406 (-0.040; 1/4) | 0.374 (+0.027; 6/2) |  |
| non-VNTR hotspot (1) | 0.895 (+0.001; 1/0) | 0.894 (n=1) | 0.947 (+0.053; 1/0) | 0.916 (+0.022; 1/0) | 0.916 (+0.022; 1/0) | 0.937 (+0.043; 1/0) | 0.936 (+0.043; 1/0) |  |  |  | 0.884 (-0.009; 0/1) |  |
| matched VNTR controls (7) | 0.983 (+0.000; 2/3) | 0.984 (n=7) | 0.977 (+0.000; 3/3) | 0.978 (-0.009; 0/6) | 0.978 (-0.009; 0/6) | 0.978 (-0.008; 1/5) | 0.973 (-0.025; 0/7) | 0.950 (-0.003; 1/4) | 0.857 (-0.114; 0/7) | 0.908 (-0.018; 0/6) | 0.984 (+0.000; 1/3) |  |
| non-repeat SV controls (2) | 0.934 (-0.000; 0/1) | 0.935 (n=2) | 0.895 (-0.040; 0/2) | 0.895 (-0.040; 0/2) | 0.930 (-0.005; 0/1) | 0.930 (-0.005; 0/1) | 0.882 (-0.053; 0/2) | 0.882 (-0.053; 0/2) | 0.930 (-0.005; 0/1) | 0.930 (-0.005; 0/1) | 0.930 (-0.005; 0/1) |  |

#### MAPQ<5, all reads the genome-wide mapping put on the anchor span (`span_input_mapq_lt5_frac`, low is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 0.505 | 0.507 | 0.474 | 0.543 | 0.497 | 0.432 | 0.516 |  | 0.534 | 0.507 | 0.382 | 0.106 |
| L014297 | hs_vntr | 0.640 | 0.646 | 0.505 | 0.662 | 0.595 | 0.608 | 0.686 |  | 0.762 | 0.765 | 0.732 | 0.254 |
| L015415 | hs_vntr | 0.603 | 0.613 | 0.566 | 0.628 | 0.608 | 0.610 | 0.620 | 0.563 | 0.679 | 0.655 | 0.654 | 0.551 |
| L005990 | hs_vntr | 0.241 | 0.218 | 0.207 | 0.196 | 0.205 | 0.200 | 0.228 |  | 0.245 | 0.231 | 0.223 | 0.123 |
| L012184 | hs_vntr | 0.704 | 0.705 | 0.515 | 0.600 | 0.649 | 0.628 | 0.544 |  |  |  | 0.525 | 0.493 |
| L002013 | hs_vntr | 0.714 | 0.748 | 0.658 | 0.701 | 0.695 | 0.652 | 0.714 |  |  | 0.751 | 0.689 | 0.507 |
| L011138 | hs_vntr | 0.937 | 0.938 | 0.803 | 0.805 |  |  | 0.824 |  |  |  | 0.875 | 0.512 |
| L012272 | hs_vntr | 0.915 | 0.905 |  | 0.076 |  |  | 0.072 |  |  |  |  | 0.047 |
| L001909 | hs_vntr | 0.593 | 0.571 | 0.543 | 0.640 |  | 0.547 | 0.603 |  |  |  | 0.547 | 0.427 |
| L016870 | hs_other | 0.155 | 0.169 | 0.025 | 0.087 | 0.037 | 0.037 | 0.106 |  |  |  | 0.037 | 0.168 |
| L016124 | c_vntr_matched | 0.029 | 0.029 | 0.035 | 0.035 | 0.035 | 0.035 | 0.033 | 0.035 | 0.038 | 0.076 | 0.044 | 0.035 |
| L007687 | c_vntr_matched | 0.262 | 0.270 | 0.260 | 0.281 | 0.283 | 0.631 | 0.635 | 0.274 | 0.509 | 0.633 | 0.264 | 0.192 |
| L007172 | c_vntr_matched | 0.264 | 0.271 | 0.260 | 0.341 | 0.512 | 0.504 | 0.443 |  | 0.672 | 0.486 | 0.357 | 0.038 |
| L011430 | c_vntr_matched | 0.190 | 0.179 | 0.150 | 0.265 | 0.190 | 0.147 | 0.253 | 0.154 | 0.248 | 0.208 | 0.154 | 0.066 |
| L007065 | c_vntr_matched | 0.030 | 0.030 | 0.046 | 0.039 | 0.028 | 0.028 | 0.039 | 0.088 | 0.284 | 0.028 | 0.028 | 0.028 |
| L004021 | c_vntr_matched | 0.052 | 0.037 | 0.013 | 0.013 | 0.013 | 0.013 | 0.031 |  | 0.250 | 0.207 | 0.006 | 0.009 |
| L014111 | c_vntr_matched | 0.014 | 0.005 | 0.005 | 0.023 | 0.005 | 0.023 | 0.044 | 0.005 | 0.005 | 0.023 | 0.005 | 0.005 |
| L000506 | c_nontr_sv | 0.031 | 0.010 | 0.010 | 0.010 | 0.010 | 0.010 | 0.010 | 0.010 | 0.010 | 0.010 | 0.010 | 0.000 |
| L006556 | c_nontr_sv | 0.064 | 0.005 | 0.201 | 0.201 | 0.005 | 0.005 | 0.201 | 0.201 | 0.005 | 0.005 | 0.005 | 0.169 |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 0.640 (-0.001; 6/3) | 0.646 (n=9) | 0.529 (-0.069; 8/0) | 0.628 (-0.022; 5/4) | 0.601 (-0.032; 6/0) | 0.608 (-0.038; 7/0) | 0.603 (+0.007; 4/5) | 0.563 (-0.050; 1/0) | 0.607 (+0.047; 0/4) | 0.655 (+0.012; 1/4) | 0.601 (-0.042; 5/3) | 0.427 (-0.242; 9/0) |
| non-VNTR hotspot (1) | 0.155 (-0.014; 1/0) | 0.169 (n=1) | 0.025 (-0.144; 1/0) | 0.087 (-0.082; 1/0) | 0.037 (-0.132; 1/0) | 0.037 (-0.132; 1/0) | 0.106 (-0.063; 1/0) |  |  |  | 0.037 (-0.132; 1/0) | 0.168 (-0.001; 1/0) |
| matched VNTR controls (7) | 0.052 (+0.000; 2/3) | 0.037 (n=7) | 0.046 (-0.010; 5/2) | 0.039 (+0.012; 1/6) | 0.035 (+0.006; 3/4) | 0.035 (+0.006; 3/4) | 0.044 (+0.038; 1/6) | 0.088 (+0.004; 2/3) | 0.250 (+0.213; 1/6) | 0.207 (+0.047; 1/6) | 0.044 (-0.002; 5/2) | 0.035 (-0.028; 6/1) |
| non-repeat SV controls (2) | 0.047 (+0.039; 0/2) | 0.008 (n=2) | 0.106 (+0.098; 0/1) | 0.106 (+0.098; 0/1) | 0.008 (+0.000; 0/0) | 0.008 (+0.000; 0/0) | 0.106 (+0.098; 0/1) | 0.106 (+0.098; 0/1) | 0.008 (+0.000; 0/0) | 0.008 (+0.000; 0/0) | 0.008 (+0.000; 0/0) | 0.085 (+0.077; 1/1) |

#### core/flank depth, all reads (`core_flank_ratio`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | - | 0.625 | 0.634 | 0.616 | 0.623 | 0.624 | 0.619 |  | 0.618 | 0.625 | 0.626 | 0.622 |
| L014297 | hs_vntr | - | 0.618 | 0.619 | 0.618 | 0.618 | 0.619 | 0.619 |  | 0.617 | 0.617 | 0.618 | 0.616 |
| L015415 | hs_vntr | - | 0.628 | 0.627 | 0.630 | 0.625 | 0.628 | 0.628 | 0.629 | 0.628 | 0.626 | 0.629 | 0.631 |
| L005990 | hs_vntr | - | 0.528 | 0.526 | 0.529 | 0.528 | 0.526 | 0.529 |  | 0.526 | 0.528 | 0.526 | 0.521 |
| L012184 | hs_vntr | - | 0.248 | 0.284 | 0.274 | 0.285 | 0.285 | 0.282 |  |  |  | 0.278 | 0.279 |
| L002013 | hs_vntr | - | 0.591 | 0.587 | 0.588 | 0.586 | 0.587 | 0.585 |  |  | 0.585 | 0.584 | 0.596 |
| L011138 | hs_vntr | - | 1.250 | 1.258 | 1.245 |  |  | 1.237 |  |  |  | 1.251 | 1.268 |
| L012272 | hs_vntr | - | 0.809 |  | 0.809 |  |  | 0.810 |  |  |  |  | 0.807 |
| L001909 | hs_vntr | - | 0.535 | 0.521 | 0.604 |  | 0.602 | 0.601 |  |  |  | 0.577 | 0.666 |
| L016870 | hs_other | - | 0.182 | 0.183 | 0.183 | 0.183 | 0.183 | 0.183 |  |  |  | 0.183 | 0.192 |
| L016124 | c_vntr_matched | - | 1.327 | 1.328 | 1.328 | 1.328 | 1.328 | 1.328 | 1.328 | 1.328 | 1.328 | 1.328 | 1.328 |
| L007687 | c_vntr_matched | - | 0.891 | 0.891 | 0.891 | 0.891 | 0.891 | 0.891 | 0.891 | 0.891 | 0.891 | 0.891 | 0.891 |
| L007172 | c_vntr_matched | - | 0.804 | 0.804 | 0.806 | 0.804 | 0.805 | 0.804 |  | 0.805 | 0.807 | 0.806 | 0.803 |
| L011430 | c_vntr_matched | - | 0.689 | 0.689 | 0.685 | 0.689 | 0.688 | 0.686 | 0.691 | 0.691 | 0.692 | 0.688 | 0.685 |
| L007065 | c_vntr_matched | - | 0.975 | 0.976 | 0.976 | 0.976 | 0.976 | 0.976 | 0.976 | 0.975 | 0.976 | 0.976 | 0.975 |
| L004021 | c_vntr_matched | - | 1.468 | 1.468 | 1.468 | 1.468 | 1.468 | 1.468 |  | 1.468 | 1.467 | 1.468 | 1.468 |
| L014111 | c_vntr_matched | - | 0.981 | 0.981 | 0.983 | 0.982 | 0.981 | 0.982 | 0.981 | 0.982 | 0.981 | 0.981 | 0.982 |
| L000506 | c_nontr_sv | - | 1.004 | 1.004 | 1.004 | 1.004 | 1.004 | 1.004 | 1.004 | 1.004 | 1.004 | 1.004 | 1.004 |
| L006556 | c_nontr_sv | - | 1.292 | 1.217 | 1.217 | 1.292 | 1.292 | 1.220 | 1.220 | 1.292 | 1.292 | 1.292 | 1.241 |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) |  | 0.618 (n=9) | 0.603 (+0.000; 4/4) | 0.616 (+0.000; 4/3) | 0.602 (-0.001; 1/3) | 0.602 (+0.000; 3/3) | 0.619 (+0.001; 5/3) | 0.629 (+0.001; 1/0) | 0.617 (-0.002; 0/3) | 0.617 (-0.001; 0/3) | 0.601 (+0.001; 5/2) | 0.622 (+0.003; 5/4) |
| non-VNTR hotspot (1) |  | 0.182 (n=1) | 0.183 (+0.001; 1/0) | 0.183 (+0.001; 1/0) | 0.183 (+0.001; 1/0) | 0.183 (+0.001; 1/0) | 0.183 (+0.001; 1/0) |  |  |  | 0.183 (+0.001; 1/0) | 0.192 (+0.010; 1/0) |
| matched VNTR controls (7) |  | 0.975 (n=7) | 0.976 (+0.000; 2/0) | 0.976 (+0.001; 4/1) | 0.976 (+0.000; 3/0) | 0.976 (+0.000; 3/1) | 0.976 (+0.000; 3/1) | 0.976 (+0.001; 3/0) | 0.975 (+0.001; 4/0) | 0.976 (+0.001; 4/1) | 0.976 (+0.000; 3/1) | 0.975 (+0.000; 2/2) |
| non-repeat SV controls (2) |  | 1.148 (n=2) | 1.111 (-0.037; 0/1) | 1.111 (-0.037; 0/1) | 1.148 (+0.000; 0/0) | 1.148 (+0.000; 0/0) | 1.112 (-0.036; 0/1) | 1.112 (-0.036; 0/1) | 1.148 (+0.000; 0/0) | 1.148 (+0.000; 0/0) | 1.148 (+0.000; 0/0) | 1.123 (-0.025; 0/1) |

#### core/flank depth, MAPQ>=5 (`core_flank_ratio_mq5`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | - | 0.282 | 0.317 | 0.202 | 0.238 | 0.285 | 0.224 |  | 0.205 | 0.232 | 0.335 | 0.566 |
| L014297 | hs_vntr | - | 0.149 | 0.323 | 0.186 | 0.197 | 0.172 | 0.163 |  | 0.105 | 0.100 | 0.129 | 0.452 |
| L015415 | hs_vntr | - | 0.178 | 0.237 | 0.177 | 0.186 | 0.194 | 0.189 | 0.227 | 0.124 | 0.132 | 0.143 | 0.350 |
| L005990 | hs_vntr | - | 0.259 | 0.291 | 0.298 | 0.288 | 0.297 | 0.248 |  | 0.222 | 0.243 | 0.255 | 0.408 |
| L012184 | hs_vntr | - | 0.058 | 0.096 | 0.058 | 0.046 | 0.056 | 0.080 |  |  |  | 0.086 | 0.103 |
| L002013 | hs_vntr | - | 0.086 | 0.145 | 0.117 | 0.118 | 0.146 | 0.109 |  |  | 0.081 | 0.123 | 0.294 |
| L011138 | hs_vntr | - | 0.052 | 0.240 | 0.216 |  |  | 0.174 |  |  |  | 0.138 | 0.756 |
| L012272 | hs_vntr | - | 0.551 |  | 0.166 |  |  | 0.202 |  |  |  |  | 0.502 |
| L001909 | hs_vntr | - | 0.093 | 0.101 | 0.090 |  | 0.112 | 0.062 |  |  |  | 0.109 | 0.255 |
| L016870 | hs_other | - | 0.136 | 0.178 | 0.162 | 0.173 | 0.175 | 0.154 |  |  |  | 0.176 | 0.176 |
| L016124 | c_vntr_matched | - | 1.309 | 1.297 | 1.298 | 1.297 | 1.297 | 1.301 | 1.297 | 1.293 | 1.220 | 1.278 | 1.297 |
| L007687 | c_vntr_matched | - | 0.622 | 0.625 | 0.608 | 0.602 | 0.217 | 0.212 | 0.611 | 0.339 | 0.213 | 0.625 | 0.688 |
| L007172 | c_vntr_matched | - | 0.555 | 0.562 | 0.464 | 0.293 | 0.305 | 0.372 |  | 0.169 | 0.319 | 0.478 | 0.777 |
| L011430 | c_vntr_matched | - | 0.555 | 0.588 | 0.475 | 0.544 | 0.585 | 0.489 | 0.589 | 0.497 | 0.538 | 0.579 | 0.652 |
| L007065 | c_vntr_matched | - | 0.950 | 0.934 | 0.938 | 0.952 | 0.952 | 0.938 | 0.883 | 0.664 | 0.952 | 0.952 | 0.953 |
| L004021 | c_vntr_matched | - | 1.413 | 1.450 | 1.450 | 1.450 | 1.450 | 1.422 |  | 1.087 | 1.148 | 1.461 | 1.457 |
| L014111 | c_vntr_matched | - | 0.981 | 0.981 | 0.965 | 0.982 | 0.963 | 0.942 | 0.981 | 0.982 | 0.963 | 0.981 | 0.982 |
| L000506 | c_nontr_sv | - | 0.976 | 0.976 | 0.976 | 0.976 | 0.976 | 0.976 | 0.976 | 0.976 | 0.976 | 0.976 | 1.004 |
| L006556 | c_nontr_sv | - | 1.287 | 1.127 | 1.127 | 1.287 | 1.287 | 1.136 | 1.136 | 1.287 | 1.287 | 1.287 | 1.188 |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) |  | 0.149 (n=9) | 0.238 (+0.049; 8/0) | 0.177 (+0.000; 4/4) | 0.192 (+0.018; 4/2) | 0.172 (+0.019; 6/1) | 0.174 (+0.011; 5/4) | 0.227 (+0.049; 1/0) | 0.164 (-0.049; 0/4) | 0.132 (-0.046; 0/5) | 0.134 (+0.022; 5/3) | 0.408 (+0.172; 8/1) |
| non-VNTR hotspot (1) |  | 0.136 (n=1) | 0.178 (+0.042; 1/0) | 0.162 (+0.026; 1/0) | 0.173 (+0.037; 1/0) | 0.175 (+0.039; 1/0) | 0.154 (+0.018; 1/0) |  |  |  | 0.176 (+0.040; 1/0) | 0.176 (+0.040; 1/0) |
| matched VNTR controls (7) |  | 0.950 (n=7) | 0.934 (+0.003; 4/2) | 0.938 (-0.014; 1/6) | 0.952 (-0.011; 3/4) | 0.952 (-0.012; 3/4) | 0.938 (-0.039; 1/6) | 0.883 (-0.011; 1/3) | 0.664 (-0.283; 1/6) | 0.952 (-0.089; 1/6) | 0.952 (+0.002; 4/2) | 0.953 (+0.044; 6/1) |
| non-repeat SV controls (2) |  | 1.131 (n=2) | 1.051 (-0.080; 0/1) | 1.051 (-0.080; 0/1) | 1.131 (+0.000; 0/0) | 1.131 (+0.000; 0/0) | 1.056 (-0.076; 0/1) | 1.056 (-0.076; 0/1) | 1.131 (+0.000; 0/0) | 1.131 (+0.000; 0/0) | 1.131 (+0.000; 0/0) | 1.096 (-0.035; 1/1) |

#### core depth from MAPQ>=5 reads (x) (`core_depth_mq5`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 9.390 | 9.164 | 10.124 | 6.929 | 8.240 | 9.956 | 7.702 |  | 7.037 | 7.994 | 11.787 | 20.150 |
| L014297 | hs_vntr | 4.192 | 4.158 | 8.792 | 4.983 | 5.500 | 4.840 | 4.351 |  | 2.743 | 2.620 | 3.407 | 13.389 |
| L015415 | hs_vntr | 4.371 | 4.175 | 5.556 | 4.152 | 4.378 | 4.530 | 4.442 | 5.434 | 2.892 | 3.133 | 3.388 | 7.263 |
| L005990 | hs_vntr | 8.245 | 8.659 | 9.743 | 9.966 | 9.644 | 9.945 | 8.310 |  | 7.425 | 8.120 | 8.530 | 13.655 |
| L012184 | hs_vntr | 1.391 | 1.390 | 2.737 | 1.652 | 1.281 | 1.572 | 2.280 |  |  |  | 2.441 | 2.930 |
| L002013 | hs_vntr | 3.533 | 2.758 | 4.745 | 3.847 | 3.871 | 4.785 | 3.568 |  |  | 2.617 | 4.038 | 9.381 |
| L011138 | hs_vntr | 1.158 | 0.993 | 4.871 | 4.483 |  |  | 3.706 |  |  |  | 2.751 | 14.218 |
| L012272 | hs_vntr | 2.399 | 2.212 |  | 5.051 |  |  | 6.152 |  |  |  |  | 15.496 |
| L001909 | hs_vntr | 1.818 | 2.429 | 2.696 | 2.026 |  | 2.936 | 1.601 |  |  |  | 2.879 | 6.701 |
| L016870 | hs_other | 2.184 | 2.163 | 2.920 | 2.645 | 2.839 | 2.882 | 2.516 |  |  |  | 2.879 | 2.647 |
| L016124 | c_vntr_matched | 39.112 | 39.163 | 38.804 | 38.817 | 38.804 | 38.804 | 38.927 | 38.804 | 38.689 | 36.479 | 38.224 | 38.804 |
| L007687 | c_vntr_matched | 17.448 | 17.283 | 17.409 | 16.937 | 16.760 | 5.947 | 5.825 | 17.034 | 9.450 | 5.843 | 17.387 | 19.409 |
| L007172 | c_vntr_matched | 19.677 | 19.483 | 19.834 | 16.373 | 10.310 | 10.702 | 13.139 |  | 5.656 | 11.258 | 16.615 | 27.442 |
| L011430 | c_vntr_matched | 17.965 | 18.039 | 19.116 | 15.433 | 17.665 | 19.022 | 15.871 | 19.134 | 16.143 | 17.479 | 18.821 | 21.167 |
| L007065 | c_vntr_matched | 27.869 | 27.855 | 27.324 | 27.517 | 27.925 | 27.925 | 27.517 | 25.905 | 19.403 | 27.925 | 27.925 | 27.933 |
| L004021 | c_vntr_matched | 33.809 | 33.527 | 34.401 | 34.401 | 34.401 | 34.401 | 33.734 |  | 25.667 | 27.177 | 34.649 | 34.573 |
| L014111 | c_vntr_matched | 24.820 | 25.015 | 25.023 | 24.583 | 25.036 | 24.567 | 24.027 | 25.023 | 25.036 | 24.553 | 25.023 | 25.036 |
| L000506 | c_nontr_sv | 25.816 | 27.430 | 27.430 | 27.430 | 27.430 | 27.430 | 27.430 | 27.430 | 27.430 | 27.430 | 27.430 | 28.222 |
| L006556 | c_nontr_sv | 32.164 | 36.642 | 29.757 | 29.757 | 36.642 | 36.642 | 29.951 | 29.951 | 36.642 | 36.642 | 36.642 | 31.826 |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 3.533 (+0.165; 7/2) | 2.758 (n=9) | 5.213 (+1.364; 8/0) | 4.483 (+0.825; 6/3) | 4.939 (+0.594; 4/2) | 4.785 (+0.682; 7/0) | 4.351 (+0.267; 6/3) | 5.434 (+1.259; 1/0) | 4.965 (-1.349; 0/4) | 3.133 (-1.042; 0/5) | 3.397 (+0.751; 5/3) | 13.389 (+6.623; 9/0) |
| non-VNTR hotspot (1) | 2.184 (+0.021; 1/0) | 2.163 (n=1) | 2.920 (+0.757; 1/0) | 2.645 (+0.482; 1/0) | 2.839 (+0.676; 1/0) | 2.882 (+0.719; 1/0) | 2.516 (+0.353; 1/0) |  |  |  | 2.879 (+0.716; 1/0) | 2.647 (+0.484; 1/0) |
| matched VNTR controls (7) | 24.820 (+0.014; 4/3) | 25.015 (n=7) | 25.023 (+0.126; 5/2) | 24.583 (-0.346; 1/6) | 25.036 (-0.359; 3/4) | 24.567 (-0.359; 3/4) | 24.027 (-0.988; 1/6) | 25.023 (-0.249; 2/3) | 19.403 (-7.833; 1/6) | 24.553 (-2.684; 1/6) | 25.023 (+0.070; 5/2) | 27.442 (+1.046; 6/1) |
| non-repeat SV controls (2) | 28.990 (-3.046; 0/2) | 32.036 (n=2) | 28.593 (-3.443; 0/1) | 28.593 (-3.443; 0/1) | 32.036 (+0.000; 0/0) | 32.036 (+0.000; 0/0) | 28.691 (-3.346; 0/1) | 28.691 (-3.346; 0/1) | 32.036 (+0.000; 0/0) | 32.036 (+0.000; 0/0) | 32.036 (+0.000; 0/0) | 30.024 (-2.012; 1/1) |

#### flank depth from MAPQ>=5 reads (x) (`flank_depth_mq5`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | - | 32.502 | 31.983 | 34.321 | 34.574 | 34.911 | 34.317 |  | 34.327 | 34.454 | 35.186 | 35.581 |
| L014297 | hs_vntr | - | 27.936 | 27.186 | 26.848 | 27.913 | 28.205 | 26.722 |  | 26.138 | 26.117 | 26.505 | 29.626 |
| L015415 | hs_vntr | - | 23.479 | 23.465 | 23.490 | 23.500 | 23.351 | 23.460 | 23.890 | 23.409 | 23.820 | 23.668 | 20.769 |
| L005990 | hs_vntr | - | 33.441 | 33.441 | 33.441 | 33.441 | 33.441 | 33.441 |  | 33.389 | 33.436 | 33.441 | 33.440 |
| L012184 | hs_vntr | - | 23.895 | 28.423 | 28.327 | 28.143 | 28.142 | 28.381 |  |  |  | 28.238 | 28.467 |
| L002013 | hs_vntr | - | 32.151 | 32.826 | 32.775 | 32.815 | 32.863 | 32.703 |  |  | 32.287 | 32.798 | 31.882 |
| L011138 | hs_vntr | - | 18.991 | 20.267 | 20.767 |  |  | 21.316 |  |  |  | 19.990 | 18.802 |
| L012272 | hs_vntr | - | 4.014 |  | 30.425 |  |  | 30.498 |  |  |  |  | 30.854 |
| L001909 | hs_vntr | - | 26.061 | 26.655 | 22.635 |  | 26.194 | 25.828 |  |  |  | 26.477 | 26.257 |
| L016870 | hs_other | - | 15.944 | 16.423 | 16.354 | 16.423 | 16.423 | 16.354 |  |  |  | 16.370 | 15.025 |
| L016124 | c_vntr_matched | - | 29.911 | 29.911 | 29.911 | 29.911 | 29.911 | 29.911 | 29.911 | 29.911 | 29.911 | 29.911 | 29.911 |
| L007687 | c_vntr_matched | - | 27.800 | 27.861 | 27.874 | 27.861 | 27.465 | 27.460 | 27.874 | 27.861 | 27.447 | 27.800 | 28.195 |
| L007172 | c_vntr_matched | - | 35.086 | 35.311 | 35.265 | 35.177 | 35.119 | 35.310 |  | 33.398 | 35.252 | 34.776 | 35.310 |
| L011430 | c_vntr_matched | - | 32.489 | 32.489 | 32.489 | 32.489 | 32.489 | 32.489 | 32.489 | 32.489 | 32.489 | 32.489 | 32.489 |
| L007065 | c_vntr_matched | - | 29.322 | 29.251 | 29.322 | 29.322 | 29.322 | 29.322 | 29.322 | 29.228 | 29.322 | 29.322 | 29.322 |
| L004021 | c_vntr_matched | - | 23.724 | 23.724 | 23.724 | 23.724 | 23.724 | 23.724 |  | 23.614 | 23.671 | 23.724 | 23.725 |
| L014111 | c_vntr_matched | - | 25.508 | 25.518 | 25.487 | 25.489 | 25.507 | 25.516 | 25.518 | 25.489 | 25.507 | 25.508 | 25.489 |
| L000506 | c_nontr_sv | - | 28.115 | 28.115 | 28.115 | 28.115 | 28.115 | 28.115 | 28.115 | 28.115 | 28.115 | 28.115 | 28.117 |
| L006556 | c_nontr_sv | - | 28.469 | 26.404 | 26.404 | 28.469 | 28.469 | 26.376 | 26.376 | 28.469 | 28.469 | 28.469 | 26.794 |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) |  | 26.061 (n=9) | 27.804 (+0.297; 4/3) | 28.327 (+0.624; 6/2) | 30.479 (+0.342; 4/1) | 28.205 (+0.269; 5/1) | 28.381 (+0.552; 5/3) | 23.890 (+0.411; 1/0) | 29.764 (-0.061; 1/3) | 32.287 (+0.136; 3/2) | 27.371 (+0.531; 6/1) | 29.626 (+0.196; 5/4) |
| non-VNTR hotspot (1) |  | 15.944 (n=1) | 16.423 (+0.479; 1/0) | 16.354 (+0.410; 1/0) | 16.423 (+0.479; 1/0) | 16.423 (+0.479; 1/0) | 16.354 (+0.410; 1/0) |  |  |  | 16.370 (+0.426; 1/0) | 15.025 (-0.919; 0/1) |
| matched VNTR controls (7) |  | 29.322 (n=7) | 29.251 (+0.000; 3/1) | 29.322 (+0.000; 2/1) | 29.322 (+0.000; 2/1) | 29.322 (+0.000; 1/2) | 29.322 (+0.000; 2/1) | 29.322 (+0.000; 2/0) | 29.228 (-0.019; 1/4) | 29.322 (+0.000; 1/3) | 29.322 (+0.000; 0/1) | 29.322 (+0.000; 3/1) |
| non-repeat SV controls (2) |  | 28.292 (n=2) | 27.259 (-1.033; 0/1) | 27.259 (-1.033; 0/1) | 28.292 (+0.000; 0/0) | 28.292 (+0.000; 0/0) | 27.245 (-1.046; 0/1) | 27.245 (-1.046; 0/1) | 28.292 (+0.000; 0/0) | 28.292 (+0.000; 0/0) | 28.292 (+0.000; 0/0) | 27.456 (-0.836; 1/1) |

#### truth-consistent core reads aligned as well as to truth (`truth_consistent_core_graph_as_good_frac`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | - | 0.954 | 0.958 | 0.941 | 0.956 | 0.961 | 0.945 |  | 0.945 | 0.952 | 0.956 | 1.000 |
| L014297 | hs_vntr | - | 0.942 | 0.925 | 0.929 | 0.933 | 0.946 | 0.925 |  | 0.925 | 0.946 | 0.938 | 1.000 |
| L015415 | hs_vntr | - | 0.989 | 0.992 | 0.989 | 0.983 | 0.986 | 0.989 | 0.977 | 0.992 | 0.986 | 0.986 | 1.000 |
| L005990 | hs_vntr | - | 0.955 | 0.964 | 0.951 | 0.953 | 0.966 | 0.932 |  | 0.951 | 0.955 | 0.962 | 1.000 |
| L012184 | hs_vntr | - | 0.903 | 0.944 | 0.892 | 0.954 | 0.960 | 0.958 |  |  |  | 0.946 | 1.000 |
| L002013 | hs_vntr | - | 0.916 | 0.915 | 0.911 | 0.920 | 0.922 | 0.911 |  |  | 0.915 | 0.911 | 1.000 |
| L011138 | hs_vntr | - | 0.938 | 0.953 | 0.958 |  |  | 0.950 |  |  |  | 0.967 | 1.000 |
| L012272 | hs_vntr | - | 0.974 |  | 0.984 |  |  | 0.979 |  |  |  |  | 1.000 |
| L001909 | hs_vntr | - | 0.778 | 0.792 | 0.844 |  | 0.846 | 0.842 |  |  |  | 0.850 | 1.000 |
| L016870 | hs_other | - | 0.990 | 0.990 | 1.000 | 1.000 | 1.000 | 0.990 |  |  |  | 1.000 | 1.000 |
| L016124 | c_vntr_matched | - | 1.000 | 1.000 | 0.998 | 0.998 | 1.000 | 1.000 | 1.000 | 0.998 | 0.997 | 1.000 | 1.000 |
| L007687 | c_vntr_matched | - | 0.998 | 1.000 | 0.998 | 1.000 | 0.995 | 0.995 | 0.995 | 0.998 | 0.995 | 1.000 | 1.000 |
| L007172 | c_vntr_matched | - | 0.996 | 0.998 | 0.987 | 0.993 | 0.993 | 0.989 |  | 0.993 | 0.993 | 0.993 | 1.000 |
| L011430 | c_vntr_matched | - | 0.981 | 0.976 | 0.978 | 0.984 | 0.984 | 0.978 | 0.978 | 0.967 | 0.973 | 0.978 | 1.000 |
| L007065 | c_vntr_matched | - | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.996 | 1.000 | 1.000 | 1.000 |
| L004021 | c_vntr_matched | - | 0.999 | 0.999 | 0.999 | 0.999 | 0.999 | 0.998 |  | 0.999 | 0.997 | 0.998 | 1.000 |
| L014111 | c_vntr_matched | - | 0.998 | 1.000 | 0.995 | 1.000 | 0.997 | 0.997 | 1.000 | 1.000 | 0.998 | 1.000 | 1.000 |
| L000506 | c_nontr_sv | - | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| L006556 | c_nontr_sv | - | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) |  | 0.942 (n=9) | 0.949 (+0.006; 6/2) | 0.941 (-0.004; 3/5) | 0.953 (-0.000; 3/3) | 0.960 (+0.007; 6/1) | 0.945 (+0.000; 4/4) | 0.977 (-0.011; 0/1) | 0.948 (-0.007; 1/3) | 0.952 (-0.002; 1/3) | 0.951 (+0.004; 5/3) | 1.000 (+0.058; 9/0) |
| non-VNTR hotspot (1) |  | 0.990 (n=1) | 0.990 (+0.000; 0/0) | 1.000 (+0.010; 1/0) | 1.000 (+0.010; 1/0) | 1.000 (+0.010; 1/0) | 0.990 (+0.000; 0/0) |  |  |  | 1.000 (+0.010; 1/0) | 1.000 (+0.010; 1/0) |
| matched VNTR controls (7) |  | 0.998 (n=7) | 1.000 (+0.001; 4/1) | 0.998 (-0.002; 1/4) | 0.999 (+0.001; 4/2) | 0.997 (+0.000; 2/3) | 0.997 (-0.002; 0/5) | 1.000 (+0.000; 1/2) | 0.998 (-0.002; 1/4) | 0.997 (-0.003; 0/5) | 1.000 (+0.000; 2/3) | 1.000 (+0.002; 5/0) |
| non-repeat SV controls (2) |  | 1.000 (n=2) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) |

#### reads mapped, of those the genome-wide mapping put on the anchor span (`span_input_mapped_frac`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 1.000 | 0.995 | 0.998 | 0.981 | 0.993 | 0.997 | 0.983 |  | 0.983 | 0.996 | 0.997 | 0.995 |
| L014297 | hs_vntr | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |  | 1.000 | 1.000 | 0.997 | 1.000 |
| L015415 | hs_vntr | 1.000 | 0.996 | 0.996 | 0.994 | 0.992 | 0.994 | 0.994 | 0.994 | 0.992 | 0.992 | 0.994 | 0.994 |
| L005990 | hs_vntr | 1.000 | 0.997 | 0.997 | 0.996 | 0.996 | 0.996 | 0.994 |  | 0.994 | 0.996 | 0.996 | 0.992 |
| L012184 | hs_vntr | 1.000 | 0.995 | 0.984 | 0.953 | 0.990 | 0.995 | 0.983 |  |  |  | 0.972 | 0.968 |
| L002013 | hs_vntr | 1.000 | 0.999 | 0.989 | 0.987 | 0.987 | 0.992 | 0.983 |  |  | 0.983 | 0.984 | 1.000 |
| L011138 | hs_vntr | 1.000 | 0.994 | 0.998 | 0.992 |  |  | 0.985 |  |  |  | 0.997 | 1.000 |
| L012272 | hs_vntr | 1.000 | 1.000 |  | 1.000 |  |  | 1.000 |  |  |  |  | 1.000 |
| L001909 | hs_vntr | 1.000 | 0.980 | 0.973 | 0.978 |  | 0.981 | 0.974 |  |  |  | 0.988 | 0.976 |
| L016870 | hs_other | 1.000 | 0.994 | 0.994 | 1.000 | 1.000 | 1.000 | 0.994 |  |  |  | 1.000 | 1.000 |
| L016124 | c_vntr_matched | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| L007687 | c_vntr_matched | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| L007172 | c_vntr_matched | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |  | 1.000 | 1.000 | 1.000 | 1.000 |
| L011430 | c_vntr_matched | 1.000 | 0.994 | 0.994 | 0.992 | 0.994 | 0.994 | 0.992 | 0.994 | 0.994 | 0.994 | 0.994 | 0.992 |
| L007065 | c_vntr_matched | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| L004021 | c_vntr_matched | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |  | 1.000 | 1.000 | 1.000 | 1.000 |
| L014111 | c_vntr_matched | 1.000 | 0.999 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| L000506 | c_nontr_sv | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| L006556 | c_nontr_sv | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 1.000 (+0.004; 8/0) | 0.996 (n=9) | 0.997 (+0.000; 2/3) | 0.992 (-0.002; 0/8) | 0.992 (-0.003; 0/5) | 0.995 (-0.001; 2/4) | 0.985 (-0.006; 0/8) | 0.994 (-0.002; 0/1) | 0.993 (-0.003; 0/3) | 0.996 (-0.001; 1/3) | 0.995 (-0.001; 3/5) | 0.995 (-0.000; 2/5) |
| non-VNTR hotspot (1) | 1.000 (+0.006; 1/0) | 0.994 (n=1) | 0.994 (+0.000; 0/0) | 1.000 (+0.006; 1/0) | 1.000 (+0.006; 1/0) | 1.000 (+0.006; 1/0) | 0.994 (+0.000; 0/0) |  |  |  | 1.000 (+0.006; 1/0) | 1.000 (+0.006; 1/0) |
| matched VNTR controls (7) | 1.000 (+0.000; 2/0) | 1.000 (n=7) | 1.000 (+0.000; 1/0) | 1.000 (+0.000; 1/1) | 1.000 (+0.000; 1/0) | 1.000 (+0.000; 1/0) | 1.000 (+0.000; 1/1) | 1.000 (+0.000; 1/0) | 1.000 (+0.000; 1/0) | 1.000 (+0.000; 1/0) | 1.000 (+0.000; 1/0) | 1.000 (+0.000; 1/1) |
| non-repeat SV controls (2) | 1.000 (+0.000; 0/0) | 1.000 (n=2) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) | 1.000 (+0.000; 0/0) |

#### reads mapped, all fetched reads (includes window-edge reads) (`mapped_frac`, high is better)

| region | stratum | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hs_vntr | 1.000 | 0.922 | 0.923 | 0.912 | 0.920 | 0.923 | 0.914 |  | 0.914 | 0.922 | 0.923 | 0.922 |
| L014297 | hs_vntr | 1.000 | 0.952 | 0.952 | 0.952 | 0.952 | 0.952 | 0.952 |  | 0.952 | 0.952 | 0.951 | 0.952 |
| L015415 | hs_vntr | 1.000 | 0.937 | 0.937 | 0.936 | 0.935 | 0.936 | 0.936 | 0.936 | 0.935 | 0.935 | 0.936 | 0.936 |
| L005990 | hs_vntr | 1.000 | 0.955 | 0.955 | 0.954 | 0.954 | 0.954 | 0.952 |  | 0.952 | 0.954 | 0.954 | 0.951 |
| L012184 | hs_vntr | 1.000 | 0.968 | 0.960 | 0.938 | 0.964 | 0.968 | 0.959 |  |  |  | 0.951 | 0.949 |
| L002013 | hs_vntr | 1.000 | 0.946 | 0.940 | 0.939 | 0.939 | 0.942 | 0.936 |  |  | 0.936 | 0.937 | 0.946 |
| L011138 | hs_vntr | 1.000 | 0.972 | 0.976 | 0.970 |  |  | 0.964 |  |  |  | 0.974 | 0.977 |
| L012272 | hs_vntr | 1.000 | 0.990 |  | 0.990 |  |  | 0.990 |  |  |  |  | 0.990 |
| L001909 | hs_vntr | 1.000 | 0.942 | 0.937 | 0.940 |  | 0.943 | 0.937 |  |  |  | 0.949 | 0.939 |
| L016870 | hs_other | 1.000 | 0.916 | 0.916 | 0.918 | 0.918 | 0.918 | 0.916 |  |  |  | 0.918 | 0.918 |
| L016124 | c_vntr_matched | 1.000 | 0.947 | 0.947 | 0.947 | 0.947 | 0.947 | 0.947 | 0.947 | 0.947 | 0.947 | 0.947 | 0.947 |
| L007687 | c_vntr_matched | 1.000 | 0.934 | 0.934 | 0.934 | 0.934 | 0.934 | 0.934 | 0.934 | 0.934 | 0.934 | 0.934 | 0.934 |
| L007172 | c_vntr_matched | 1.000 | 0.932 | 0.932 | 0.932 | 0.932 | 0.932 | 0.932 |  | 0.932 | 0.932 | 0.932 | 0.932 |
| L011430 | c_vntr_matched | 1.000 | 0.948 | 0.948 | 0.947 | 0.948 | 0.948 | 0.947 | 0.948 | 0.948 | 0.948 | 0.948 | 0.947 |
| L007065 | c_vntr_matched | 1.000 | 0.957 | 0.957 | 0.957 | 0.957 | 0.957 | 0.957 | 0.957 | 0.957 | 0.957 | 0.957 | 0.957 |
| L004021 | c_vntr_matched | 1.000 | 0.971 | 0.971 | 0.971 | 0.971 | 0.971 | 0.971 |  | 0.971 | 0.971 | 0.971 | 0.971 |
| L014111 | c_vntr_matched | 1.000 | 0.956 | 0.958 | 0.958 | 0.957 | 0.958 | 0.958 | 0.958 | 0.957 | 0.958 | 0.957 | 0.957 |
| L000506 | c_nontr_sv | 1.000 | 0.905 | 0.905 | 0.905 | 0.905 | 0.905 | 0.905 | 0.905 | 0.905 | 0.905 | 0.905 | 0.905 |
| L006556 | c_nontr_sv | 1.000 | 0.904 | 0.904 | 0.904 | 0.904 | 0.904 | 0.904 | 0.904 | 0.904 | 0.904 | 0.904 | 0.904 |

| group: median (paired diff vs mc; better/worse) | orig | mc | mafft_linsi | unit_aware | poa_spoa | poa_abpoa | unit_aware__all | mafft_linsi__all | poa_spoa__all | poa_abpoa__all | poa_abpoa_mc__all | truth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| VNTR hotspots (9) | 1.000 (+0.048; 9/0) | 0.952 (n=9) | 0.946 (+0.000; 2/3) | 0.940 (-0.002; 0/8) | 0.945 (-0.002; 0/5) | 0.943 (-0.001; 2/4) | 0.952 (-0.005; 0/8) | 0.936 (-0.001; 0/1) | 0.943 (-0.002; 0/3) | 0.936 (-0.001; 1/3) | 0.950 (-0.001; 3/5) | 0.949 (-0.000; 2/5) |
| non-VNTR hotspot (1) | 1.000 (+0.084; 1/0) | 0.916 (n=1) | 0.916 (+0.000; 0/0) | 0.918 (+0.003; 1/0) | 0.918 (+0.003; 1/0) | 0.918 (+0.003; 1/0) | 0.916 (+0.000; 0/0) |  |  |  | 0.918 (+0.003; 1/0) | 0.918 (+0.003; 1/0) |
| matched VNTR controls (7) | 1.000 (+0.052; 7/0) | 0.948 (n=7) | 0.948 (+0.000; 1/0) | 0.947 (+0.000; 1/1) | 0.948 (+0.000; 1/0) | 0.948 (+0.000; 1/0) | 0.947 (+0.000; 1/1) | 0.948 (+0.000; 1/0) | 0.948 (+0.000; 1/0) | 0.948 (+0.000; 1/0) | 0.948 (+0.000; 1/0) | 0.947 (+0.000; 1/1) |
| non-repeat SV controls (2) | 1.000 (+0.095; 2/0) | 0.905 (n=2) | 0.905 (+0.000; 0/0) | 0.905 (+0.000; 0/0) | 0.905 (+0.000; 0/0) | 0.905 (+0.000; 0/0) | 0.905 (+0.000; 0/0) | 0.905 (+0.000; 0/0) | 0.905 (+0.000; 0/0) | 0.905 (+0.000; 0/0) | 0.905 (+0.000; 0/0) | 0.905 (+0.000; 0/0) |
