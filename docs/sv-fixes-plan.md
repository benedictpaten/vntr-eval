# Plan: the SV errors left on the realigned chr20 graph

**Status: planned, not started (2026-10-03).** The measurements behind it are in
[results/full_graph_chr20.md](../results/full_graph_chr20.md). The analysis scripts are in
`work/full/chr20/investigate/residual_sv/`.

## What the errors are

After realignment the HG002 chr20 calls have 562 short-read and 533 ONT raw SV errors (FP + FN).
- **94-97% are inside realigned tandem-repeat regions.**
- **They are concentrated.** 23 loci hold half of the short-read and ONT errors combined, and the 3% of
  regions with more than 100 alleles hold 67-71%.
- **The 270 regions with an error fall into three classes.** These come from comparing called and truth
  haplotypes by edit distance:
  - **Representation.** The called haplotype is within 10 bp of the truth, but its records are cut
    differently. This is 46% of ONT errors and 16% of short-read errors.
  - **Not fixable by realignment.** This is about half of the real errors:
    - HG002's allele is more than 1% from every panel allele (14 regions, 27-30%);
    - or one truth haplotype was not sampled (14-25%).
  - **Fixable in the graph or the caller.**
    - On short reads, 154 errors are in regions where the sampled graph holds both truth haplotypes.
    - 111 of those, in 35 regions, are calls that recombine haplotypes inside one region.
- **A haplotype-level check confirms the ONT gain but not the short-read gain.** The nesting-aware
  `score_haplotypes` gives:
  - ONT: +0.068 refined SV F1 [+0.013, +0.125];
  - short reads: +0.019 [−0.030, +0.064].

## The fixes, in order

Each fix needs a pass/fail check. The checks are listed after the fixes.

1. **Write SV records from the called haplotypes inside realigned regions.**
   - What: align each called haplotype to CHM13 with fixed affine scores, left-normalise, and write the
     records from that alignment. It is a post-processing step on the vg call VCF, so no vg change is
     needed.
   - Targets: the representation class.
2. **Make the star consistent with the reference.**
   - What: centre on CHM13 wherever the medoid star's induced CHM13-versus-allele alignment has gap runs
     of 50 bp or more that abPOA's pairwise alignment lacks. Alternatively, after merging, move
     insertion slots so each induced pairwise alignment matches the pairwise one.
   - Targets: cancelling INS/DEL pairs (33 short-read and 50 ONT FPs), and TR_chr20_18259720, whose
     ONT call is exact but whose star gives 706 bp of spurious gaps.
3. **Drop nested duplicates.**
   - What: drop a child record when the overlapping parent record on the same haplotype already carries
     an SV of the same sign and size within 20%.
   - Targets: 19-20 FPs.
4. **Genotype each realigned region as one site.**
   - What: the alleles are whole sampled-haplotype traversals between the region's anchors. Children are
     chosen only within the chosen parent traversal.
   - Targets: short-read recombinant calls (111 errors). It overlaps the parked work in
     [repeat-sites.md](repeat-sites.md).
5. **Fix the centre in expanded regions.** Either:
   - remove the medoid's 64-comparator cap, or sample comparators stratified by length (the cap is
     active in 898 chr20 regions); or
   - centre on the longest allele when the star expands more than 3x the longest allele. On the 40 most
     expanded regions this gives 1.595 → 0.278 Mb, at pair F1 −0.021.

   Targets: TR_chr20_64970081, whose sampled graph is 37.4 kb against 9.7 kb and whose call is a 616 bp
   allele where the truth is 4,669 bp.
6. **Treat novel-allele regions separately.**
   - What: report the 14 novel-allele regions on their own, since a per-record benchmark is fragile there.
     Fixing them needs allele discovery (local ONT assembly or augmentation), not graph realignment.
7. **Correct the scoring.**
   - What: report truvari refine's harmonised counts (`--write-phab`, or the `out_*` columns of
     `refine.regions.txt`), not `refine.variant_summary.json`.
   - Why: regions refine judges to hold no SV keep their original labels, which is about 116 errors per
     arm. This is not yet checked independently.

**Not to do:**
- cap growth by falling back to the original subgraph, since the expanded regions carry most of the
  gain;
- slot pooling;
- progressive POA centres;
- a higher convex gap open.

## How each fix is judged

- Run the same chr20 test chain on both graphs, using `test_arm.sh`.
- Gate on the per-region paired bootstrap and on the haplotype-level gain from `score_haplotypes`, not
  on raw truvari F1.
- Fixes 1, 3 and 7 change only records, so the haplotype-level numbers must not move.
- Fixes 2, 4 and 5 change genotypes, so check them on chr6 before adopting them.
