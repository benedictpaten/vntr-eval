# Stage 4: patch chr20's hotspots with abPOA graphs and re-call the whole contig

**Result.** Splicing the abPOA graphs of 49 detected tandem-repeat hotspots into the chr20 hap32 graph
cuts truvari SV false positives inside the patched regions from 176-179 to 124 (-30%). False negatives
do not change, and the SV FPs outside the regions stay within the null's range. chr20 SV F1 rises from 0.5325
(unpatched, re-mapped) and 0.5364 (renumbering null) to **0.5542**. Production's F1 is 0.5365.
Small-variant F1 does not move: SNV 0.9851 in every arm, and indel 0.9270 against 0.9273-0.9275, a
difference inside the null's own shift.

Summed called-haplotype edit distance (ED) over the 49 regions does **not** improve. It is 37,302 patched,
34,426 in the renumbering null and 31,990 unpatched. The median region improves by 7 edits, and 28 of 49
regions improve against the null. A few large regressions outweigh them, and so do the 17 regions that the
unpatched call already gets almost right. The increase over the null (+2,876) is about the size of the null's
own shift from unpatched (+2,436). So the ED result is noise-level. It does not reproduce Stage 3's 22-40% ED cut.

## Method

1. **Detection** (results/stage4_detect.md). The rule kmer_frac_extra >= 0.32 was fitted on the 145
   non-chr20 regions and applied to 680 prefiltered chr20 TR regions. It detected 57 regions, which hold
   47% of chr20's TR SV FPs. abPOA (poa_abpoa) realigned 50 of them. The other 7 were refused as predicted
   memouts: the D20S16 satellite and six ALR packages.
2. **Patch** (tools/patch_contig.py `patch`). The source is chr20.gbz, dumped to GFA with `vg convert -f`.
   Each region's inside set comes from `gbz-base --between` (anchors included), and is replaced by the
   candidate graph. Segments longer than 1,024 bp are chopped first.
   - **Paths.** Every one of the 137 GBZ paths is kept (CHM13, 9 GRCh38 fragments, 127 recombination
     fragments). Where a path crosses a region anchor to anchor, its walk there becomes the candidate
     path of the same name. It is flipped where the GBZ stores the path reversed. W-line fields 1-5 are
     unchanged and the RS tag is kept, so CHM13 and GRCh38 remain reference samples.
   - **Unpatched region.** A region is left out if a path fragment starts or ends inside it. That
     excluded one region: TR773254, a 6.7 kb VNTR with 13 production SV FPs, where 10 hap32
     fragment ends fall inside. So 49 regions are patched. Together they remove 74,781 nodes and add
     32,995.
   - **Node IDs.** Each region's new node IDs are spread over the free IDs of [anchor_left, anchor_right],
     which is the region's own range. All 49 regions interleave.
   - **Assertions (all pass).** Checked with `vg paths -F`: every path spells the same sequence as in
     chr20.gbz, compared by md5 per path over all 137 paths. The path metadata are identical. The
     GBZ's node count and length equal the GFA's. The GBZ's CHM13 path carries the GFA's node IDs, so
     there is no segment translation.
   - **Unpatched arm.** It is written by the same code with no regions. It also passes the identity
     check against chr20.gbz.
3. **Reads** (`reads`). These are the reads production placed on chr20, rebuilt from GAF-Base with
   remap_local.fetch's logic (path sequence + cs, base qualities, mates together).
   - Retrieval used 34 CHM13 chunks of 2 Mb with 10 kb context. The 256,302 nodes that no chunk
     covered were then queried by ID.
   - Keeping one copy per (name, mate) gives 13,278,235 reads: 6,600,310 pairs and 77,615 singles.
     GAF-Base holds 13,278,000-13,280,000 alignments in chr20's handle range, so the rebuild is complete.
     None were dropped.
4. **Map** (`map`). This follows call_local.map_hybrid: `vg autoindex -n -w sr-giraffe` on each arm's
   GFA, then giraffe paired with Stage 2's fragment-length distribution (401.7 ± 166.4) plus single-end
   reads, writing GAF with --named-coordinates, bgzipped. All arms map the same reads, and every GAF
   node is a GFA segment (asserted). Unmapped reads: 2,442 unpatched and 2,398 patched.
5. **Call** (`call`). This uses the pinned vg 2a6a228a5 with call_local.vg_call_command:
   `-p CHM13#0#chr20 -d 2 -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased --mosaic-out ... --gaf-reads`,
   run at -t 8. Each call took about 380 s. Peak memory was 14-18 GB.
6. **Score** (`score`, tools/stage4_report.py).
   - **Harness.** vg-call-eval's bench_wgs.score_contig, as the pipeline runs it: aardvark for small
     variants, and truvari bench with --sizemin/--sizefilt 50 and --pick ac after multiallelic
     splitting. There is no refine step, because the harness does not run one. F1 uses both sides' TP
     (bench_metrics). Indels are aardvark's JointIndel row.
   - **Inside vs outside.** Truvari FP and FN records overlapping a patched region's anchor-to-anchor
     span ± 100 bp are counted as inside.
   - **Haplotype ED.** tools/score_haplotypes.py batch over the 49 regions (--no-truvari). The 34
     regions with a production SV FP count as hotspots; the other 15 are detected regions with no FP.
7. **Renumbering null** (`mcnull`). The same 49 regions are patched with their own MC graph (mc.gfa),
   then re-mapped, called and scored. The graph is identical to the source and only the region node IDs
   are re-laid by the patch machinery. This arm measures what numbering alone does through vg call's
   4,096-ID depth-rate windows, and through giraffe.

## Results (chr20, whole contig)

| | production | unpatched | renumbering null | **patched** |
|---|---|---|---|---|
| SV precision | 0.5017 | 0.4937 | 0.5006 | **0.5324** |
| SV recall | 0.5765 | 0.5778 | 0.5778 | 0.5778 |
| SV F1 | 0.5365 | 0.5325 | 0.5364 | **0.5542** |
| SV FP / FN | 428 / 324 | 445 / 323 | 432 / 323 | **382** / 323 |
| SV FP in patched / outside | 181 / 247 | 179 / 266 | 176 / 256 | **124** / 258 |
| SV FN in patched / outside | 134 / 190 | 131 / 192 | 133 / 190 | 133 / 190 |
| SNV F1 | 0.9851 | 0.9851 | 0.9851 | 0.9851 |
| Indel F1 (JointIndel) | 0.9288 | 0.9275 | 0.9273 | 0.9270 |
| Indel FP / FN | 1635 / 1355 | 1667 / 1378 | 1687 / 1368 | 1692 / 1378 |
| ALL small F1 | 0.9725 | 0.9721 | 0.9721 | 0.9720 |
| hap ED, 49 patched regions | 42,134 | 31,990 | 34,426 | 37,302 |
| hap ED, 34 hotspots (prod FP > 0) | 24,658 | 21,901 | 24,228 | 26,782 |
| hap ED, 15 detected with no FP | 17,476 | 10,089 | 10,198 | 10,520 |

All 49 regions together have an ED of 67,711 when no variant is called at all (the reference).

**Per-region ED, patched vs the null.** The median change is -7 edits: 28 regions are better, 20 worse
and 1 the same.
- **Large gains:** TR763614 (-4,004, 10 production FPs), TR772945 (-807) and TR759446 (-584).
- **Large regressions:** TR756034 (+3,373; 14.4 kb, 23 FPs) and TR772976 (+3,110). At TR772976 the patched
  call writes a nested record with a 3,745 bp ALT over a 528 bp REF. That makes both haplotypes about
  4.6 kb long, against truth's 2.8 kb, so the ED ends up worse than the reference's.
- **Regions MC already gets right:** 17 regions have unpatched ED <= 50, summing to 275. Patching them
  adds about 2,900 edits in total, for example TR768271 (2 to 1,587) and TR772531 (2 to 884).

**What the null shows.**
- Renumbering alone moves summed ED by +2,436 (at most 1,058 in one region), total SV FP by -13
  (only -3 inside the regions) and indel FP by +20.
- So the SV FP drop inside the patched regions (-52 against the null) is well outside the numbering
  noise. The ED change and the indel change are not.

## Caveats

- **Node IDs and depth windows.** vg call measures its read rate over 4,096-node-ID windows, so node
  numbering changes calls. Patched IDs are spread over each region's original ID range, and all 49
  regions interleave. The renumbering null measures what is left, and ED moves by thousands of edits
  from numbering alone.
- **Contig-local vs genome-wide mapping.** Both arms re-map production's chr20 reads to chr20 only.
  Reads cannot move to or from other contigs, and the minimizer index is chr20's. Re-mapping by itself
  moves the result away from production: summed ED over these regions falls from 42,134 to 31,990, and
  SV FP rises by 17, 19 of them outside the regions. Compare patched with unpatched or the null, not
  with production.
- **hap32-only patch.** Only the 34 hap32 paths are threaded. Truth is used only for scoring, and for
  choosing anchors as in Stages 0-3. One hotspot (TR773254) could not be patched because hap32
  fragment ends fall inside it, and 7 detected regions had no candidate.
- **Candidate anchors are not shared nodes.** abPOA can fold an anchor base into the next node, so
  some candidate paths do not visit the candidate's first node. An example is TR755508's
  recombination#3#chr20#1. Sequences are still identical (asserted), but the region's edge is then not
  a single-node snarl boundary for every haplotype.
- **Single replicate.** Patched and unpatched have one call each. ED differences below about 2.5k in
  total, or about 1k in one region, are within the measured numbering noise.
- **Different intervals.** Inside/outside uses the anchor-to-anchor span ± 100 bp. Detection's truth
  overlap used the scan interval ± 100 bp. Production's 181 FPs in the patched regions equal
  detection's 194 for the 50 realigned regions minus TR773254's 13.

## Files

- tools/patch_contig.py (gfa, patch, reads, map, call, score) and tools/stage4_report.py
- results/stage4_chr20.tsv (per arm) and results/stage4_chr20_regions.tsv (per region: production FP/FN and each arm's ED)
- work/stage4/chr20/:
  - {unpatched,patched,mcnull}/{graph.gfa,graph.gbz,patch.json,patch.tsv,reads.gaf.gz,chr20.vcf.gz}
  - reads/ (FASTQ.gz, reads.json)
  - score/<arm>/ (harness output, haps/, summary.json)
  - run_rest.sh and run_mcnull.sh
