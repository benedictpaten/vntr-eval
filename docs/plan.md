# Plan: cleaner VNTR alignments in the pangenome graph

Why: ~1% of VNTR regions hold half of short-read `vg call`'s SV false positives, and in them the
Minigraph-Cactus graph (as sampled into the hap32 panel) misaligns the panel haplotypes: the pairwise
alignment the graph implies costs 1.56x the optimum, against 1.09x at matched error-free VNTRs
([findings.md](findings.md)). The hypothesis to test first is the simplest one: **re-aligning the
haplotype sequences of each region with an appropriate multiple aligner, and rebuilding the region's
graph from that alignment, makes it much cleaner** -- and a cleaner graph lets reads be placed and
alleles be genotyped.

## What "gross construction faults" means

Most hotspots are aligned but badly: the repeat units of different haplotypes are stacked out of
phase. A minority show something cruder -- homologous sequence that was never aligned at all, or
assembly artefacts carried into the graph:

1. **Parallel duplicated paths.** A group of haplotypes runs through its own copy of a stretch that
   the other haplotypes share with CHM13. At chr4:191.50 Mb, 12 of 30 sampled haplotypes bypass
   1,254 bp of CHM13 through a parallel path 96.7% identical to it, including 542 bp of unique
   flank that is identical in CHM13, both HG002 haplotypes and other panel haplotypes. Reads on
   either copy of that flank are 76-85% MAPQ < 5 (7-17% on the flank that is not duplicated).
2. **N-gap alleles.** Assembly gaps carried in as alleles: at the same locus one haplotype carries
   a 3.96 kb segment containing 3,146 bp of N.
3. **Paths that end inside the array.** Haplotype fragments stopping inside the VNTR (contig ends):
   8/40 hotspots against 1/40 matched controls.
4. **A lone divergent haplotype shredded into many SV-sized bubbles.** L000601 has 24 false
   positives and no truth SV.

A good multiple alignment fixes (1) and (4) by construction, since it aligns all homologous
sequence. (2) needs the Ns masked or the gap haplotypes left out of the alignment; (3) needs
fragments added to an alignment of the spanning haplotypes (`mafft --addfragments`) or left out.

## The alignment methods

The input to every method is the same: for each region, the sequences of all panel haplotypes that
span it from a left anchor node to a right anchor node (anchors are CHM13 nodes every path visits
once, outside any variation), so all sequences share identical flanks. Each method produces a
multiple alignment; the region's graph is then induced from it (`vg construct -M msa.fa`, then
compacted), with every haplotype as a path, and compared with the Minigraph-Cactus subgraph.

- **M0, baseline**: the Minigraph-Cactus graph as sampled (hap32), read as an alignment.
- **M1, mafft FFT-NS-2**: fast progressive alignment; the cheap baseline MSA.
- **M2, mafft L-INS-i**: iterative refinement with consistency from local pairwise alignments. The
  most accurate mafft mode for tens of sequences of a few kb; assumes one alignable domain.
- **M3, mafft E-INS-i** (generalized affine pairwise, with unalignable regions allowed): made for
  sequences with several conserved blocks separated by long indels, which is what arrays of
  different lengths look like. Also G-INS-i (global) as a check.
- **M4, repeat-unit-aware two-level alignment** (built here): split each array into repeat units
  using the region's motif, cluster the unit variants into symbols, align the unit strings with
  `mafft --anysymbol` (so whole units are gapped, and units stay in phase), then expand each unit
  column to bases by aligning the unit variants within it; flanks aligned normally.
- **M5, partial-order alignment**: abPOA 1.5.7 (adaptive banded, convex gaps; optionally a progressive
  guide tree `-p` and minimizer seeding `-S` for long regions) and spoa 4.1.5 (unbanded), both in
  global mode (spoa's default is local, which would leave the shared flanks unaligned). POA aligns
  each sequence to the graph of those before it, so it never has to force a repeat expansion into one
  place the way a progressive profile alignment can, but the result depends on input order: a pilot
  of ~10 regions picks one configuration and one ordering per tool, and measures how much order
  matters.

**Two input panels.** The baseline is not an alignment of the 34 hap32 sequences: Minigraph-Cactus
aligned the full HPRC v2.1 panel, and hap32 is 32
haplotypes sampled from that graph (`vg haplotypes`, personalised to HG002) plus CHM13 and GRCh38. So
`mc.gfa` is "align everything, keep the sampled rows". Within these regions every sampled sequence
checked is an exact copy of one real haplotype (180/180 at six hotspots), and a hotspot's full panel
holds 280-440 distinct sequences against about 30 in hap32.

Both panels come from the **eval** graph, `hprc-v2.1-mc-chm13-eval.gref.gbz`: 228 HPRC individuals
plus CHM13 and GRCh38, with HG002 and its parents HG003 and HG004 absent (checked from the GBWT's
sample list). A graph that contained HG002's own haplotypes would have an inherent advantage when
genotyping HG002. The `gref_CHM13` path is excluded too: it is the reference-cover path, a chimera
of other haplotypes, not a haplotype Minigraph-Cactus aligned.

Each method is therefore run twice:

- **m(hap32)**: align the 34 hap32 sequences. The easier problem; what a post-hoc fix of the
  sampled graph would do.
- **m(all) projected to hap32** (`<method>__all`): align the deduplicated full panel, keep the hap32
  rows, drop all-gap columns, build the graph. The like-for-like replacement for MC.

The difference between the two measures what the extra haplotypes do to the alignment of the
sampled ones; the full-panel graphs are also compared with MC's full-panel subgraph directly. The
all-pairs mafft modes (L-/E-/G-INS-i) will not fit hundreds of multi-kb sequences everywhere, so
the full-panel arm adds FFT-NS-i and records every region a method could not finish; comparisons are
made only on regions where both graphs exist.

For N-gap haplotypes: replace N runs by gaps in the alignment (they become path gaps, not nodes).
For fragments: add with `--addfragments` after the spanning set is aligned.

Caveats kept in view: a repeat MSA is one of many equal-cost optima, and different optima change how
the truth VCF's records line up with graph bubbles; arrays with differently ordered unit types cannot
be stacked acyclically; inducing a graph from columns can recreate a mesh of 1 bp nodes, so the
graph is compacted and, if needed, small bubbles are merged.

## How a candidate graph is judged

The yardstick matters: raw truvari F1 cannot credit a better representation here (the exact truth,
written the way the current graph breaks it up, scores 0.27). So:

- **Stage 0, graph only** (no reads): nodes per kb and node lengths; the graph-implied pairwise
  alignment cost over the optimal pairwise edit distance (unit and affine); SV-sized pieces per
  haplotype; k-mer redundancy (the same sequence on different nodes); whether each HG002 haplotype
  can be spelled by a path, and at what edit cost; "the truth written by the graph" scored with
  truvari, raw and after `truvari refine`. Target: hotspots move to the matched controls'
  distribution (cost/optimum <= 1.1), and the controls do not change.
- **Stage 1, existing reads**: placement redundancy of the region's current reads -- how many node
  paths their sequence has in the graph, against how many copies one haplotype carries.
- **Stage 2, local re-map**: map the same reads to the old and new local graphs (giraffe); MAPQ < 5
  fraction, core/flank depth, how much of the reads the best haplotype pair explains.
- **Stage 3, local genotyping**: `vg call` on the local graphs; score by the called haplotypes' edit
  distance to HG002 and by refined truvari F1; all controls must hold.
- **Stage 4**: splice the rebuilt regions into one contig, re-index, re-map and re-call it.

Test set: the 150 stratified loci from the census (40 FP-dense VNTR hotspots, 40 FP-free VNTRs
matched on length and period, 25 VNTRs vg calls correctly, 25 non-VNTR hotspots, 20 non-repeat SV
controls), each packaged as sequences in [../regions](../regions); plus the truth-free target list
of all benchmarked VNTRs >= 1 kb (3,111 regions) as coordinates, for scaling up.

Longer-term options, once the local result is known: genotype each problem VNTR as one site of whole
alleles (caller side); do the re-alignment upstream in Minigraph-Cactus (mask VNTRs before Cactus,
re-insert them as one bubble aligned as here); a genome-wide detector for the gross faults above.
