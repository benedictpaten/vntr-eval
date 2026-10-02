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

## Stage 4b: realign every eligible repeat (no detection rule)

All 669 packaged eligible TR regions (>= 50 bp haplotype length variation) realigned with abPOA
defaults on hap32: 636 got a graph (33 refused on predicted memory), 624 patched (12 refused because a
panel-haplotype fragment ends inside). Same reads, mapping and vg call as the other arms
(`work/stage4/chr20/run_all.sh`; scores in `work/stage4/chr20/score_all/`).

| chr20 | production | unpatched | patched, 49 detected | patched, 624 eligible |
|---|---|---|---|---|
| SV precision / recall | 0.5017 / 0.5765 | 0.4937 / 0.5778 | 0.5324 / 0.5778 | **0.5908** / 0.5765 |
| SV F1 | 0.5365 | 0.5325 | 0.5542 | **0.5835** |
| SV FP / FN | 428 / 324 | 445 / 323 | 382 / 323 | **302** / 324 |
| SV FP inside / outside the 624 regions | 385 / 43 | 390 / 55 | - | 246 / 56 |
| SNV F1 / indel F1 | 0.9851 / 0.9288 | 0.9851 / 0.9275 | 0.9851 / 0.9270 | 0.9853 / 0.9282 |

Eligible regions left unpatched hold 30 production SV FPs: 12 patch-refused (fragment ends, 19 FP),
33 memory-refused (5 FP), 11 not packageable (6 FP).

## Stage 4c/4d: full-panel abPOA and the motif-guided aligner, every eligible repeat

Same 669 eligible regions, reads, mapping, vg call and scoring as Stage 4b (`work/stage4/chr20/run_full.sh`,
`run_units.sh`; scores in `work/stage4/chr20/score_{all,full,units}/`).
- Full panel: the HG002-free eval panel fetched per region (660 of 669; spans > 100 kb skipped), all
  distinct sequences aligned with abPOA defaults, projected onto the 34 hap32 rows: 633 aligned (36 over
  the memory cap), 621 patched.
- Motif-guided (repeat-unit-aware) on hap32, fallback abPOA where the motif is unusable (Stage 3 showed
  the mafft fallback breaks simple SVs): 659 candidates, 629 patched.

| chr20 | production | unpatched | abPOA hap32 (624) | abPOA full panel -> hap32 (621) | motif-guided hap32 (629) |
|---|---|---|---|---|---|
| SV precision / recall | 0.5017 / 0.5765 | 0.4937 / 0.5778 | **0.5908** / 0.5765 | 0.4675 / 0.4928 | 0.4762 / 0.5556 |
| SV F1 | 0.5365 | 0.5325 | **0.5835** | 0.4798 | 0.5128 |
| SV FP / FN | 428 / 324 | 445 / 323 | **302** / 324 | 426 / 388 | 462 / 340 |
| SV FP inside / outside patched | 385 / 43 | 390 / 55 | 246 / 56 | 358 / 68 | 412 / 50 |
| SV FN inside / outside patched | 278 / 46 | 275 / 48 | 277 / 47 | 328 / 60 | 305 / 35 |
| SNV F1 / indel F1 | 0.9851 / 0.9288 | 0.9851 / 0.9275 | 0.9853 / 0.9282 | 0.9836 / 0.9181 | 0.9853 / 0.9266 |

(Inside/outside use each arm's own patched set.) Read mapping is the same in every arm (first 3M GAF
records: unmapped 0.02%, MAPQ < 5 6.64-6.69%) and only 7-14 regions per arm fell back to appended node
IDs, so the differences come from how the patched graphs represent alleles, not from mapping.

- **abPOA on hap32 is the only arm that helps**: FP -143, FN flat, small variants flat or slightly up.
- **Aligning the full panel and projecting hurts calling**: FN +65, SNV F1 -0.0015, indel F1 -0.0094,
  FP only -19. The ~460-way alignment is the harder problem (Stage 0 saw deletions cut into unit-sized
  pieces when aligned inside the full panel), and here that costs recall and small-variant accuracy.
  Where the small-variant errors fall (inside vs outside patched regions) was not checked.
- **The motif-guided aligner, applied to every eligible repeat, is net negative** (FP +17, FN +17),
  consistent with Stage 3's record-level truvari: its haplotypes are close to the truth but written as
  records the truth VCF does not use, and most eligible regions are not hotspots.

## Stage 4e: why the full-panel abPOA arm hurts, and three variants

**Result.** When hap32's rows are aligned inside the 460-way POA, they come out aligned worse to each
other. The projected graphs roughly double the repeated-k-mer fraction (`kmer_frac_extra`, the Stage 0
metric that tracked read placement): 0.212 against 0.104 for abPOA on hap32. On chr20 this costs nothing
outside the patched regions. Inside them it adds small-variant false positives. None of the variants that
let the panel shape the hap32 alignment beats `poa_abpoa__all`. The one variant that keeps hap32's own
alignment turns out to be `poa_abpoa` itself once projected, so no new chr20 arm was run.

**Where chr20's small-variant losses fall** (aardvark records, SNVs and indels < 50 bp; "inside" means
within ± 100 bp of the 621 patched_full spans; the patched_all spans give the same counts):

| chr20 | SNV FP in / out | indel FP in / out | SNV FN in / out | indel FN in / out |
|---|---|---|---|---|
| unpatched | 101 / 264 | 172 / 1424 | 73 / 1773 | 89 / 1289 |
| abPOA hap32 (patched_all) | 96 / 262 | 177 / 1415 | 60 / 1764 | 95 / 1280 |
| abPOA full panel (patched_full) | **329** / 261 | **611** / 1423 | 69 / 1765 | 122 / 1281 |

- Every loss is a false positive inside the patched regions. Outside them the counts do not change.
- The net +667 FPs are spread over 71 regions that got worse (21 got better). The top 10 regions hold 56%.
- This fits with the arm's SV FN +65: in the patched regions the caller picks a wrong haplotype.

**Stage 0 on the 149 packaged regions** (138 paired; the other 11 are predicted memouts for every
full-panel arm). Each cell is the median `all_cost_over_opt` / `kmer_frac_extra`. "b/w" counts regions
better / worse than `poa_abpoa__all` on cost/opt by more than 0.005.

| stratum (n) | mc | abPOA hap32 | abPOA full (`__all`) | v1: full, `-p` | v3: `__all` polished |
|---|---|---|---|---|---|
| hotspot VNTR (30) | 1.582 / 0.583 | 1.083 / 0.414 | 1.130 / 0.557 | 1.210 / 0.546 (5/23) | 1.281 / 0.562 (4/25) |
| hotspot other (23) | 1.653 / 0.440 | 1.164 / 0.367 | 1.244 / 0.515 | 1.277 / 0.498 (8/15) | 1.313 / 0.466 (1/21) |
| matched control (40) | 1.102 / 0.053 | 1.024 / 0.038 | 1.084 / 0.062 | 1.078 / 0.062 (13/11) | 1.042 / 0.054 (17/4) |
| correct control (25) | 1.015 / 0.041 | 1.006 / 0.031 | 1.077 / 0.136 | 1.033 / 0.082 (10/6) | 1.072 / 0.076 (4/3) |
| all (138) | 1.176 / 0.143 | 1.037 / 0.104 | 1.100 / 0.212 | 1.103 / 0.206 (36/55) | 1.124 / 0.161 (26/53) |

- **v1, abPOA `-p` (guide tree) on the full panel.** It helps the controls but is worse at the hotspots
  (VNTR hotspots: 5 better, 23 worse). Overall it is no better than `__all`.
- **v3, polishing.** Starting from the `__all` projection, each distinct hap32 row is removed and
  re-added with `abpoa -i` onto the other rows' MSA, one round. This helps matched controls but is
  worse at the hotspots, because a cheap alignment to a graph path is not a good pairwise alignment.
- **v2, hap32 first, then the other panel sequences added with `abpoa -i`.** abPOA keeps the hap32
  rows' alignment exactly: in all 21 regions tested, the projection has the same column graph as the
  stage-1 MSA. So the projection is `poa_abpoa`'s graph. In 19 of 21 regions it is identical. In the
  other 2, equal-length sequences enter in union order instead of hap32.fa order (confirmed by
  re-running stage 1 in that order). On a 25-region sample it scores 1.065 / 0.109, the same as
  `poa_abpoa`, against `__all`'s 1.129 / 0.270. For a patch that threads only hap32 paths, v2 therefore
  *is* patched_all (SV F1 0.5835). Running it on chr20 would only measure tie-order noise, so it was
  not run, and there is no patched_fullv arm.
- **Why.** In a POA each sequence is aligned to the graph built so far. In `__all`, hap32 rows enter
  among hundreds of other alleles. A row can then match its repeat units through other sequences'
  nodes, so identical k-mers in two hap32 rows end up on different nodes. The damage grows with panel
  diversity: across regions, the extra `kmer_frac_extra` of `__all` over hap32 correlates with the
  number of distinct panel sequences (Spearman 0.45, n = 118). The median extra is 0.015 in the half
  with ≤ 163 distinct sequences and 0.110 in the other half.

Files: `work/fullpanel/` holds fp.py (variants), batch.py, sig.py (column-graph identity), table.py,
inout.py, and `candidates/<variant>/` and `results/<variant>/` (Stage 0 JSON, no truth section).

## Stage 4f: truth-free snarl anchors and threaded fragments (abPOA on hap32, every eligible repeat)

Same 680 prefiltered TR loci, reads, giraffe mapping, pinned vg 2a6a228a5 call and scoring as Stage 4b
(`work/stage4/chr20/run_v2.sh`; packages in `work/stage4/regions_v2`, candidates in `work/stage4/candidates_v2`,
arm in `work/stage4/chr20/patched_v2`, scores in `work/stage4/chr20/score_v2`). Two things change:
- **Anchors** (`package_regions.py build --anchor-mode snarl`, tools/snarl_anchors.py). The anchors are the
  boundaries of the smallest snarl, or run of consecutive snarls of one chain, that encloses the padded
  (+-200 bp) locus. The snarls come from `vg snarls -T -P CHM13` on chr20.gbz (2.25 M snarls, trivial ones
  included; cached in `work/stage4/snarls`). Only boundary nodes that CHM13 visits once, forward, are used.
  No truth record is consulted. A locus whose enclosing span exceeds 250 kb is skipped. Loci with
  overlapping spans are merged into one region, named after the first locus (`merged_loci` in region.json).
  At a contig end the interval is clamped to the first or last boundary node.
- **Fragments** (`realign.py --fragments-poa`, `patch_contig.py patch --fragments`). A hap32 path piece that
  starts or ends inside a region is added to the region's finished abPOA MSA with `abpoa -i`. This covers
  hap32.fragments.fa plus the pieces that are only an anchor node. Prefixes use extension mode, suffixes
  use extension mode on the reversed MSA, and internal pieces use local mode. The fragment rows are
  carried back into the original columns, so the spanning rows' alignment is unchanged (tested). The
  candidate graph is built from all rows. The patch then threads each fragment's run, from its first to
  its last step inside the region, along the candidate path of the same name. All 137 paths still spell
  their chr20.gbz sequence (md5, asserted).

**Coverage**

| chr20 | 4b: truth-rule anchors | 4f: snarl anchors + fragments |
|---|---|---|
| packaged | 669 of 680 loci (11 failed) | 657 regions = 677 loci (14 merged regions hold 34 loci; 3 skipped) |
| abPOA graph | 636 (33 over the memory cap) | 642 regions = 652 loci (15 regions = 25 loci over the cap) |
| patched | 624 (9 overlap an earlier region, 3 have a fragment end inside) | **642 of 642** (none refused) |
| fragments threaded | - | 14 paths in 4 regions (TR773254: 10, TR756022: 2, TR755542, TR763163) |
| production SV FP in eligible loci left unpatched | 30 | 5 (all in TR768070, over the cap) |

- **Skipped (span > 250 kb, 0 production FP):** TR762632 is the 3.18 Mb centromeric array. TR763173 needs 954 kb and TR772623 needs 1.06 Mb.
- **Over the cap:** 15 regions, each with a haplotype of 22-238 kb.
- **Anchor agreement:** of the 649 regions both builds share, 603 get the same anchors as the truth rule,
  36 get a narrower span and 10 a wider one. The total anchor-to-anchor span is 1.83 Mb against 3.80 Mb,
  because the old fallback anchors reached far.

**Calls**

| chr20 | production | unpatched | 4b: abPOA hap32 (624) | **4f: snarl anchors + fragments (642)** |
|---|---|---|---|---|
| SV precision / recall | 0.5017 / 0.5765 | 0.4937 / 0.5778 | 0.5908 / 0.5765 | 0.5898 / 0.5647 |
| SV F1 | 0.5365 | 0.5325 | **0.5835** | 0.5770 |
| SV FP / FN | 428 / 324 | 445 / 323 | 302 / 324 | 297 / 333 |
| SV FP inside / outside own patched regions | 385 / 43 | 390 / 55 | 246 / 56 | 263 / 34 |
| SV FN inside / outside own patched regions | 278 / 46 | 275 / 48 | 277 / 47 | 309 / 24 |
| SV FP / FN inside the 642 4f regions | 404 / 302 | 411 / 299 | 268 / 301 | 263 / 309 |
| SNV F1 / indel F1 | 0.9851 / 0.9288 | 0.9851 / 0.9275 | 0.9853 / 0.9282 | 0.9852 / 0.9283 |

(The production and unpatched inside/outside rows use 4b's 624 regions.)

- **One region explains the whole F1 drop.** TR773254 is a 7.1 kb VNTR with 13 production FPs, where 10 of
  the 34 hap32 haplotypes break inside the repeat. It could not be patched in 4b. Patched with its 20
  threaded pieces, its FP/FN goes from 14/14 to 15/23. That is the whole FN increase (+9 of +9) and
  -0.0065 SV F1.
- **Everywhere else the arm matches 4b:** FP -5, FN +0, which is inside the renumbering null's +-13 FP.
  TR756022, the other threaded region with calls, goes from 1/1 to 0/0.
- **Small variants do not move.**
- **What this does not settle:** it cannot tell whether TR773254's recall loss comes from the threaded
  fragments or from the region's abPOA graph. Without the fragments the region cannot be patched, so
  there is no control.

So the snarl anchors do what the truth rule did without consulting truth: 603/649 identical anchors, and
every realigned region patched. Threading fragments in the one fragment-heavy hotspot costs recall there
and gains nothing elsewhere. Caveats: single replicate, and truvari refine was not run on this arm.

## Stage 4h: full-panel test-set iteration

**Result.** On a fixed 22-region chr20 test set, no variant of the full-panel arm closes the gap to abPOA on
hap32. abPOA's `-G` (path-score heuristic), combined with a higher mismatch penalty, removes up to about
half of the extra repeated k-mers. But the called haplotypes barely improve: abPOA on hap32 stays
reproducibly better at 10-14 regions and worse at 1-2 against every full-panel variant. `-G` alone is the
best full-panel variant for calling, with 8 regions better and 5 worse than `poa_abpoa__all`. Merging
parallel nodes after alignment changes nothing, so the duplication is not identical parallel nodes.

**Test set** (`work/iterate/testset.tsv`; per-region chr20 SV FP+FN inside the patched_full span ± 100 bp,
from `work/iterate/region_sv.py`, which reproduces the 98/111 (+44) and 115/60 (-131) counts):
- **8 regions where the full-panel arm hurt most**: +38 errors in total (the hap32 arm: -11).
- **6 where it helped most**: -27 (the hap32 arm: -19).
- **4 hotspots**: 10-12 production SV FPs each.
- **4 where MC is already right**: 0 production FP/FN and 2 matched truth SVs each, with 41-186 distinct
  panel alleles.
- **On chr20 these 22 regions hold**: unpatched 132 SV errors, patched_full 136, patched_all 91.
- **Every region's full-panel abPOA finishes in under 20 s.** The hotspot TR756034 (106 s, 5.5 GB) was left out.

**Harness** (`tools/iterate.py`). Each variant goes through the same steps:
1. **Build.** The region's full-panel union goes through the same `realign.align_fasta` masking, projection
   (`panel.project`) and `msa_graph.py` path checks as `poa_panel.py`.
2. **Stage 0.** `evaluate.py --skip truth`.
3. **Call.** `call_local.py` hybrid200k, hybrid50k and hybrid200kids, i.e. 3 replicates, with reads re-mapped.
4. **Score.** `score_haplotypes.score_region` without truvari; the score is `ed_suppress_nested`.
5. **Compare.** The replicate rule of `stage3_rule.py`.

`fp_default` rebuilds `poa_abpoa__all`'s full MSA byte for byte (same Stage 0), so the pipeline reproduces
Stage 4c. The fragment-length setting is Stage 2's documented mean 402 / sd 166; the original
`fraglen.json` was lost, and a re-fetch gives 399 / 166. At TR772413 the hybrid50k replicate fails for
every graph: a 1-bp flank node on no walk is dropped from the GBZ. That region therefore has 2 replicates.

**Baselines on the test set** (medians over the 22 regions; ED = summed median called-haplotype edit
distance; b/w = regions reproducibly better/worse under the replicate rule):

| graph | kmer_frac_extra | cost/opt | nodes/kb | summed ED | exact | b/w vs mc | b/w vs poa_abpoa | b/w vs poa_abpoa__all |
|---|---|---|---|---|---|---|---|---|
| mc | 0.259 | 1.408 | 215 | 6,840 | 4 | - | 5/9 | 9/5 |
| abPOA hap32 (`poa_abpoa`) | 0.175 | 1.057 | 104 | 10,222 | 5 | 9/5 | - | 11/2 |
| abPOA full panel (`poa_abpoa__all`) | 0.286 | 1.117 | 136 | 7,741 | 2 | 5/9 | 2/11 | - |

The harness reproduces the direction of the chr20 result:
- **By region counts, abPOA on hap32 beats the full-panel arm 11/2 and MC 9/5.** The full-panel arm
  loses to MC 5/9.
- **In the "hurt" group, the full-panel arm is worse than MC at 5 of 8 regions** and better at 1.
- **The "helped" group does not reproduce against MC** (1 better, 3 worse). Its chr20 gains were 3-7
  truvari records per region.
- **Summed ED is dominated by single blow-ups**, so the region counts are the verdict. abPOA on hap32
  has 3,532 at TR772976 and 3,106 at TR772999 (one replicate there is 515).

**Sweep** (full panel, projected; 22 regions each; ranked by median `kmer_frac_extra`). "gap" is the
share of the mean `kmer_frac_extra` difference between `poa_abpoa__all` (0.320) and `poa_abpoa`
(0.237) that the variant removes. "b/w" counts regions whose `kmer_frac_extra` moved by more than 0.005.

| variant | abPOA flags / structure | kmerx median / mean | gap | kmerx b/w vs `__all` | cost/opt median / mean | nodes/kb |
|---|---|---|---|---|---|---|
| fp_GX8O6 | -G -X 8 -O 6,26 | 0.234 / 0.283 | 45% | 15/1 | 1.104 / 1.224 | 111 |
| fp_GX6 | -G -X 6 | 0.237 / 0.288 | 38% | 15/2 | 1.112 / 1.221 | 96 |
| fp_GX10 | -G -X 10 | 0.239 / 0.282 | 46% | 14/2 | 1.084 / 1.224 | 109 |
| fp_GX8 | -G -X 8 | 0.240 / 0.288 | 39% | 14/3 | 1.084 / 1.224 | 109 |
| fp_GX8nb | -G -X 8 -b -1 | 0.242 / 0.274 | **55%** | 17/1 | 1.110 / 1.183 | 100 |
| fp_X8 | -X 8 | 0.256 / 0.301 | 22% | 12/3 | 1.105 / 1.225 | 105 |
| fp_GX5 | -G -X 5 | 0.256 / 0.288 | 39% | 15/2 | 1.103 / 1.215 | 105 |
| fp_GX6nb | -G -X 6 -b -1 | 0.259 / 0.276 | 53% | 14/2 | 1.094 / 1.173 | 95 |
| fp_GX6O6 | -G -X 6 -O 6,26 -E 2,1 | 0.273 / 0.296 | 29% | 13/3 | 1.113 / 1.213 | 102 |
| fp_nb, fp_b100 | -b -1; -b 100 -f 0.1 (identical graphs) | 0.279 / 0.314 | 6% | 6/2 | 1.123 / 1.212 | 134 |
| fp_G | -G | 0.280 / 0.305 | 17% | 9/5 | 1.107 / 1.211 | 114 |
| fp_Gnb | -G -b -1 | 0.285 / 0.295 | 30% | 12/3 | 1.107 / 1.218 | 109 |
| s2 merges (5) | identical siblings; zip (same seq + same preds or succs); `vg mod -U 10`; before or after projection | 0.286 / 0.319-0.320 | 0-1% | 0/0 | 1.117 / 1.203 | 134-136 |
| fp_J | -J | 0.286 / 0.319 | 0% | 1/0 | 1.117 / 1.202 | 136 |
| `poa_abpoa__all` | defaults | 0.286 / 0.320 | 0% | - | 1.117 / 1.203 | 136 |
| fp_GX6M1 | -G -M 1 -X 3 | 0.287 / 0.302 | 21% | 10/3 | 1.101 / 1.170 | 111 |
| fp_O6 | -O 6,26 -E 2,1 | 0.291 / 0.323 | -4% | 4/8 | 1.121 / 1.221 | 127 |
| fp_freq | most frequent distinct allele first | 0.294 / 0.347 | -33% | 10/12 | 1.318 / 1.430 | 123 |
| fp_X6 | -X 6 | 0.297 / 0.306 | 16% | 11/4 | 1.119 / 1.227 | 118 |
| s1_seed20 | 20 most frequent + longest + shortest, then `abpoa -i` | 0.315 / 0.324 | -6% | 8/9 | 1.163 / 1.385 | 109 |
| fp_aff | -O 6,0 -E 2,0 (affine) | 0.319 / 0.323 | -4% | 7/8 | 1.103 / 1.174 | 135 |
| fp_R, fp_RJ | -R; -R -J (identical) | 0.328 / 0.342 | -27% | 2/14 | 1.121 / 1.191 | 133 |
| s1_seed10 | 10 most frequent + longest + shortest, then `abpoa -i` | 0.334 / 0.343 | -28% | 8/9 | 1.378 / 1.504 | 139 |
| fp_short | shortest first | 0.357 / 0.360 | -48% | 6/13 | 1.105 / 1.200 | 179 |

What moves the duplication:
- **`-G` combined with a mismatch penalty of 6 or more is the only lever that removes much of it.**
  Alone, each does little (`-G` 17%, `-X 8` 22%), but together they remove 38-46%, and 53-55% unbanded.
  Unbanded without `-G` does little (6%).
- **Gap placement does nothing or harms:** `-J` changes nothing and `-R` adds duplication. Stronger gap
  penalties and affine gaps also do not help.
- **Input orders other than longest-first hurt.** Frequency-first and the seeded variants also raise
  cost/opt (1.3-1.5 mean), because long alleles are then aligned globally onto graphs of short ones. A
  sample-independent seed does not reproduce v2's effect: the hap32 rows still enter among hundreds of
  other alleles.
- **Node merging finds almost nothing to merge.** Identical-sibling merge: 0 nodes in all 22 regions.
  Zip: 517 nodes in the full graphs, 213 after projection. `vg mod -U 10`: -1.5% nodes. None changes
  `kmer_frac_extra` by more than 0.0005, and every path still spells its sequence. The duplicated units
  sit at different offsets on different branches, not on parallel identical nodes.

**Local calling of the best variants** (3 replicates each):

| graph | summed ED | exact | b/w vs mc | b/w vs poa_abpoa | b/w vs poa_abpoa__all |
|---|---|---|---|---|---|
| mc | 6,840 | 4 | - | 5/9 | 9/5 |
| poa_abpoa | 10,222 | 5 | 9/5 | - | 11/2 |
| poa_abpoa__all | 7,741 | 2 | 5/9 | 2/11 | - |
| fp_G (-G) | 7,406 | 4 | 6/7 | 2/10 | **8/5** |
| fp_GX8 (-G -X 8) | 7,766 | 4 | 4/8 | 2/11 | 7/4 |
| fp_GX8nb (-G -X 8 -b -1) | 7,685 | 4 | 6/6 | 2/10 | 5/3 |
| fp_GX6 (-G -X 6) | 9,228 | 4 | 6/9 | 1/12 | 6/6 |
| fp_GX6nb (-G -X 6 -b -1) | 11,600 | 3 | 5/9 | 1/14 | 4/8 |

- **No full-panel variant approaches abPOA on hap32**: at best 2 regions better against 10 worse.
- **`-G` is a modest gain over `poa_abpoa__all`** (8/5, summed ED -335, exact regions 2 -> 4). It is not
  better than MC (6/7).
- **Less duplication does not buy better calls.** fp_GX6 and fp_GX6nb have the least duplication and call
  no better than `-G` alone. Their summed ED is inflated by single regions: fp_GX6 has 3,149 at TR755910,
  where every other graph has 604.
- **Stage 0 barely ranks the graphs.** Within a region, the Spearman correlation of `kmer_frac_extra`
  with the median ED across the 8 called graphs has median 0.09; for cost/opt it is 0.17.
- **The replicates disagree in 56 of 176 region x graph cells**, so a single call would have been
  misleading.

So the full-panel arm's calling loss is not the repeated k-mers that abPOA's knobs remove. If the full
panel is to be used, `-G` is the setting to use. Otherwise hap32-only alignment remains the arm to patch with.

Files: `tools/iterate.py`; `work/iterate/`:
- `testset.tsv`, `region_sv.py` / `region_sv.tsv` (per-region chr20 SV errors, all patched regions)
- `results.tsv` (region x variant), `summary.tsv`, `paired.py` (paired Stage 0 vs `__all`)
- `candidates/<variant>/`, `panel/<variant>/` (full MSAs), `stage0/`, `score/`
- `w/stage3/` (calls; hybrid graphs, indexes and GAFs cleaned after scoring)
- `stage2/fraglen.json`

## Stage 4i: coarsened sites, and where the called-haplotype error comes from

Same 22-region test set and 3-replicate local vg call as Stage 4h (`tools/iterate.py`; results in
`work/iterate/results.tsv`, `coarsen_table.txt`).

**Coarsening does not help.** Across the Stage 4h variants, within-region call error tracked the number of
top-level snarls (more sites, worse calls in 18 of 20 regions; median Spearman +0.37), so we tried fewer,
larger sites: `coarsen` keeps runs of >= k columns where every hap32 row agrees as anchors and makes each
distinct row sequence between anchors one allele (k = 16/32/64, applied to abPOA hap32, full panel, full
panel -G). Every coarsened graph calls far worse (summed ED 23,400-28,700 vs mc 6,840 and abPOA hap32
10,222; loses most regions to abPOA hap32), as expected once whole alleles repeat shared sequence
(repeated-k-mer fraction 0.75-0.86): reads can no longer be placed.

**Where the error is.** Called haplotypes compared with the truth (best pairing) and with the 34 hap32 panel
sequences, final-arm calls, 22 regions (44 haplotypes per graph):

| graph | called haplotypes identical to a panel haplotype | ED from those | ED from off-panel (stitched) haplotypes | VCF records |
|---|---|---|---|---|
| mc | 55% | 907 | 5,940 | 862 |
| abPOA hap32 | 52% | 798 | 9,435 | 1,084 |
| abPOA full panel | 45% | 849 | 7,006 | 1,646 |
| full panel -G | 45% | 884 | 6,539 | 1,480 |

In every graph 87-92% of the error comes from called haplotypes that no panel haplotype carries: the caller
combines alleles of different panel haplotypes across adjacent sites. The full-panel graphs have more sites
(more records) and produce such haplotypes slightly more often, which is how they lose calls. The ceiling
for calling only panel haplotypes is far lower: each truth haplotype's nearest panel haplotype is 973 edits
away in total over the 22 regions (exact in 7 regions), against 6,840-10,222 for the calls.

**Implication.** At these regions the limiting step is how vg call builds alleles across the sites of one
repeat, not the graph alignment. Genotyping each repeat as one site whose alleles are the panel haplotypes
(scored with reads aligned to the merged, fine-grained graph) could cut the error several-fold; the graph
realignment choice matters much less.

## Stage 4j: prototype genotyping of a repeat as one site of panel-haplotype alleles

`tools/panel_genotype.py` (+ `tools/c/semig.c`), offline, same 22-region test set. Alleles = the distinct hap32
panel sequences over the region (anchor to anchor). Reads = the reads the genome-wide mapping placed in the
region's window that share a 21-mer with the region; each gets its semi-global edit distance d to every
allele (both strands) and is dropped if > 8% of its length from all of them. A diploid pair (h1, h2) scores
sum_r log(0.5 e^(-3 d_r1) + 0.5 e^(-3 d_r2)) plus w x a Poisson depth term: the number of kept reads against
kappa (L1 + L2 + 2(R - 1)), with kappa (read starts per bp per haplotype) estimated from reads lying wholly
in 2 kb of CHM13 flank on each side (0.09-0.10, as 30x predicts). No truth is used except for scoring.
Single run (vg call numbers are 3-replicate medians).

| method | summed ED to truth | regions better / worse vs this row's baseline |
|---|---|---|
| ceiling: nearest panel haplotype to each truth haplotype | 973 | – |
| panel genotyper, read content only (w = 0) | 9,458 | vs mc 7/8 |
| panel genotyper, depth weight w = 0.3 | 5,107 | vs mc 9/7 |
| **panel genotyper, depth weight w = 1** | **4,127** | vs mc 9/8; vs abPOA hap32 10/7; vs abPOA full 13/6 |
| vg call on mc | 6,840 | – |
| vg call on abPOA hap32 | 10,222 | – |
| vg call on abPOA full panel | 7,741 | – |

- Without a depth term, content alone cannot tell copy number (reads inside an array fit every allele), and
  the prototype is no better than vg call. With it, summed error falls 40% below vg call on mc and 60% below
  vg call on the abPOA hap32 graph, reaching the ceiling in 7 of 22 regions.
- Region counts are only even against mc (9 better, 8 worse): the gain is concentrated in regions where vg
  call builds badly wrong alleles (TR772999: 467 vs 1,621), while the prototype still makes a few large
  length errors (TR772976: called 2,828/1,283 against truth 2,810/2,801). The depth weight was chosen from two
  values on this same set.
- So whole-allele genotyping of repeats is promising but not yet a clear win. Next: use mate pairs (insert
  size spans arrays and measures length), base qualities and mapping quality in the read likelihood; test on
  held-out regions; then, if it holds, genotype chains of snarls inside repeats jointly in vg call itself.

## 4k. Is the linkage prior too weak? (test set, `--linkage-weight`)

Question: 87-92% of the called-haplotype error in 4i came from called haplotypes that no panel haplotype
carries (recombinant paths), yet vg call has a Li-Stephens linkage HMM that should penalise switching
between panel haplotypes. Test: raise `--linkage-weight` (default 2) and re-call the 22-region test set,
3 replicates, pinned vg 2a6a228a5 (`tools/iterate.py`, kind `link`; the flag reaches vg through
`VG_CALL_EXTRA` in `tools/call_local.py`, checked in each call's command line).

The `lwN_*` rows build the graph through the candidate path (GFA -> GBZ, new node IDs), so `lw2_mc` is mc
with only its node layout changed: 7,779 against 6,840 for mc with its native IDs. Node layout alone moves
summed ED by ~14% here (the depth-rate window was keyed on node ID in this binary; see the PR's positional
window). Rows are therefore compared at the same weight and the same layout path.

| graph | summed ED, w=2 | w=4 | w=8 | w=8 vs w=2, better / worse |
|---|---|---|---|---|
| mc (candidate layout) | 7,779 | 7,294 | 6,508 | 5/3 |
| abPOA hap32 | 10,222 | 7,014 | 5,933 | 7/5 |
| abPOA full panel, projected | 7,741 | 8,059 | 5,883 | 6/5 |

| at w=8 | better / worse |
|---|---|
| abPOA hap32 vs mc | 6/4 (ED -9%) |
| abPOA full panel vs mc | 3/9 |
| abPOA full panel vs abPOA hap32 | 2/11 |

Called haplotypes found in the panel (44 per graph) and the error split between in-panel and off-panel calls:

| graph | in panel, w=2 -> w=8 | ED in / off panel, w=2 | w=8 |
|---|---|---|---|
| mc (candidate layout) | 24 -> 32 | 907 / 6,879 | 1,944 / 4,604 |
| abPOA hap32 | 23 -> 31 | 798 / 9,435 | 1,144 / 4,936 |
| abPOA full panel | 20 -> 26 | 849 / 7,006 | 1,432 / 4,358 |

- Yes, the prior is too weak for these regions. Quadrupling the weight cuts the recombinant calls
  by a third and the summed error by 16-42% in every graph. The gain is concentrated in a few regions (most
  region counts are ties under the replicate rule), mostly the hotspots: abPOA hap32 hotspot ED 4,440 -> 1,510.
- The weight helps the aligned graphs most. That fits the idea that finer, more even snarls give the
  HMM more junctions at which to switch haplotypes.
- More linkage does not rescue the full-panel projection: it has the lowest summed ED but loses region-wise
  to both mc and abPOA hap32. Its duplicated repeat units (4g) still cost whole regions.
- Still off-panel: 12-18 of 44 called haplotypes, holding 75-80% of the remaining error. w=8 is the top of
  this sweep, not a fitted optimum. Whether it costs anything outside repeats is a whole-contig question; see 4l.

## 4l. Linkage weight 8 on the whole of chr20

The four Stage 4 chr20 arms re-called at `--linkage-weight 8` with their existing graphs and read mappings
(`work/iterate/lw_chr20.sh`, pinned vg 2a6a228a5, one draw each). Haplotype ED is summed over the 624
patched repeat regions scored in 4b (median-free, single draw; better/worse = regions with lower/higher ED).

| arm | weight | SV F1 | SV FP / FN | SNV F1 | indel F1 | ALL F1 | haplotype ED, 624 regions |
|---|---|---|---|---|---|---|---|
| unpatched (mc) | 2 | 0.5325 | 445 / 323 | 0.9851 | 0.9275 | 0.9721 | 88,099 |
| unpatched (mc) | 8 | 0.5335 | 438 / 324 | 0.9834 | 0.9182 | 0.9688 | 81,015 |
| mc re-laid (mcnull) | 2 | 0.5317 | 438 / 326 | 0.9851 | 0.9273 | 0.9721 | 84,912 |
| mc re-laid (mcnull) | 8 | 0.5417 | 428 / 319 | 0.9834 | 0.9178 | 0.9687 | 82,235 |
| abPOA hap32 patched | 2 | 0.5835 | 302 / 324 | 0.9853 | 0.9282 | 0.9724 | 93,529 |
| abPOA hap32 patched | 8 | 0.5861 | 304 / 319 | 0.9833 | 0.9174 | 0.9686 | 81,378 |
| full panel patched | 2 | 0.4798 | 426 / 388 | 0.9836 | 0.9181 | 0.9687 | – |
| full panel patched | 8 | 0.4775 | 428 / 389 | 0.9820 | 0.9091 | 0.9655 | – |

- The test-set result generalises to the repeats: at weight 8 the haplotype ED over the 624 patched
  regions falls 8% on mc and 13% on the patched graph (95/56 and 92/62 regions better/worse).
- SV F1 barely moves (+0.001 to +0.010), and the patching gain is unchanged (+0.053 at both weights).
- Small variants pay for it genome-wide: SNV F1 -0.0017, indel F1 -0.009, ALL F1 -0.0033 in every arm
  (unpatched: +108 SNV FP, +139 SNV FN, +175 indel FP, +201 indel FN). A global weight of 8 is a bad trade.
- So the prior is too weak inside repeats and about right elsewhere. The lever to test is a
  repeat-conditional linkage weight, not a new default.

## 4m. Gap placement in the full-panel MSA

Question: the hap32 projection of the full-panel abPOA MSA (`poa_abpoa__all`) has more repeated k-mers than
abPOA run on the hap32 rows alone (`poa_abpoa`): kmerx (`kmer_frac_extra`) median 0.286 against 0.175. Is
this inconsistent gap placement, meaning two rows with the same unit content whose indels sit in different
columns? If so, normalising each row's gaps against the column consensus should remove it. Tools:
`tools/gap_norm.py` (with tests in `tools/test_gap_norm.py`) and `iterate.py` kind `gapnorm`. Stage 0 was
run on the 22-region test set.

**Diagnosis** (`gap_norm.py diagnose`; `work/iterate/gapnorm/diagnose22.json`). In a column graph, a 21-mer's
positions are the columns where its occurrences start. This reproduces Stage 0's kmerx exactly. A pair of
occurrences is *split* when abPOA on hap32 puts them in one column and the projection puts them in two. For
each split pair, the two rows' induced pairwise alignments are compared between the nearest base pairs
that both MSAs align. Over that stretch both MSAs align the same two subsequences, so the costs compare
directly.

| | TR773368 | TR770009 | TR765271 | TR755896 | all 22 |
|---|---|---|---|---|---|
| extra k-mer positions, projection minus hap32-only | 136 | 675 | 450 | 1,358 | 47,063 |
| explained by equal-cost placement (split positions merged by equal-cost pairs) | 0 | 0 | 113 | 0 | ≤ 2,080 (≤ 4.4%) |
| split pairs: projection's alignment worse / equal / better | 534 / 0 / 0 | 1,084 / 0 / 0 | 263 / 113 / 224 | 3,516 / 0 / 0 | 88% / 1% / 10% |
| gap runs over distinct hap32 rows, projection vs hap32-only | 8 vs 2 | 191 vs 36 | 336 vs 169 | 56 vs 21 | 33,878 vs 20,171 |

- **Gap placement is not the cause.** Equally good alignments that place an indel elsewhere explain at
  most 4.4% of the excess, and 0% in 3 of the 4 regions examined by hand. In 88% of split pairs, the
  projection aligns the two rows worse than hap32-only does over the stretch between shared anchors (summed
  cost +16%).
- **The projection fragments indels.** Each hap32 row's indels are cut into 1.7 times as many gap runs,
  for the same total gap length. Example: in TR773368 the 816 bp allele's 419 bp deletion becomes six
  pieces, against one in the hap32-only MSA.
- **The fragments are not noise.** 99.7% of the hap32 rows' gap runs in the full MSA have exactly the same
  start and end as a run of some non-hap32 panel allele. abPOA aligns each sequence to its best path through
  a graph of up to 460 alleles. Each unit of a short allele therefore matches whichever panel allele fits it
  best, and its deletion breaks at those alleles' junctions. The full MSA is self-consistent. The projection
  keeps the junctions but drops most of the alleles that justified them.
- **How the extra k-mers arise.** 86% of the extra positions are k-mers that span a gap in their row at
  that position (junction k-mers). The rest are unit copies placed against differently varied units.
- **The panel supports these junctions.** A sum-of-pairs merge test on the full MSA confirms this.
  `gap_norm.merge_runs` moves the bases between two of a row's gap runs so as to join the runs. It accepts
  a move costing up to λ edits per other row (λ = 0-2). On the 4 regions it moves kmerx by at most 0.015,
  and on the 22 regions `gml_poa_abpoa__all` (λ = 1, then left normalisation) does no better than left
  normalisation alone.

**Variants** (Stage 0, 22 regions; b/w = regions where the metric fell / rose by more than 0.005; lower is
better for both metrics):

| variant | MSA processing | kmerx median / mean | kmerx b/w vs `__all` | cost/opt median / mean | cost/opt b/w vs `__all` |
|---|---|---|---|---|---|
| mc | - | 0.259 / 0.307 | 11/7 | 1.408 / 1.507 | 4/18 |
| `poa_abpoa` | abPOA on hap32 alone | 0.175 / 0.237 | 20/1 | 1.057 / 1.092 | 16/2 |
| `poa_abpoa__all` | none | 0.286 / 0.320 | - | 1.117 / 1.203 | - |
| gnl_poa_abpoa__all | full MSA, per-row left normalisation | 0.308 / 0.320 | 12/3 | 1.183 / 1.259 | 3/15 |
| gnr_poa_abpoa__all | full MSA, per-row right normalisation | 0.338 / 0.356 | 2/18 | 1.151 / 1.247 | 5/14 |
| gnl_fp_G | `-G` full MSA, left (vs fp_G: kmerx 11/3, cost 2/17) | 0.289 / 0.305 | 13/6 | 1.144 / 1.261 | 7/14 |
| gnr_fp_G | `-G` full MSA, right (vs fp_G: 1/19, 6/13) | 0.308 / 0.341 | 6/16 | 1.130 / 1.247 | 5/12 |
| gml_poa_abpoa__all | full MSA, merge (λ = 1), then left | 0.308 / 0.322 | 8/4 | 1.175 / 1.249 | 4/13 |
| gjl_poa_abpoa__all | full MSA, joint left: a run shared by several rows moves only if it can in all | 0.282 / 0.318 | 5/1 | 1.127 / 1.213 | 2/5 |
| gnpl_poa_abpoa__all | the projection, left, against the hap32 consensus (sample-dependent; reference only) | 0.266 / 0.303 | 20/0 | 1.150 / 1.245 | 4/14 |

- **Per-row normalisation lowers kmerx in most regions and raises cost/opt in most.**
  - It scatters gap runs that rows shared, because each row stops at its own first variant base. Over the
    22 regions, distinct gap runs rise from 4,394 to 6,380, though total runs fall from 33,878 to 28,892.
  - Where that happens, kmerx rises too. In TR772945 it goes from 0.276 to 0.426.
  - Right normalisation is worse on both metrics.
- **Joint normalisation avoids the scattering and is nearly inert:** kmerx −0.004, cost/opt +0.010. That is
  the most a placement fix can do, consistent with the ≤ 4.4% ceiling above.
- **Even on the projection itself, normalisation does not close the gap.** The sample-dependent version
  reaches 0.266, against 0.175 for hap32-only, and still raises cost/opt.
- **No variant beats `poa_abpoa__all` on kmerx without raising cost/opt**, so none was called.

So the full-panel projection's duplication is not a gap-placement artefact. Its junctions are real
junctions of the full panel, which a 34-row subset does not need. No normalisation that leaves the panel's
alignment intact can remove them. Two directions remain:
- **A different alignment objective for the full panel.** Score each row against the column profile (sum of
  pairs, affine gaps, iterative refinement) instead of against abPOA's best path, so that an indel is cut
  only where many alleles support the cut.
- **Accept the sample dependence.** Re-align, on the hap32 rows, each projected window whose junctions no
  two hap32 rows share.

Files: `tools/gap_norm.py`, `tools/test_gap_norm.py`; `work/iterate/gapnorm/diagnose{4,22}.json`; candidates
and full MSAs under `work/iterate/candidates/<variant>/` and `work/iterate/panel/<variant>/`.

## 4n. A different alignment objective: centre-star and profile aligners

Question: 4m traced the full-panel projection's duplication to abPOA cutting each allele's indels at
junctions that other panel alleles support. Does an objective that stacks indels in shared columns give a
sample-independent full-panel MSA whose hap32 projection is as good as aligning hap32 alone? Gate for
calling: median kmerx at least 0.03 below `poa_abpoa__all` (≤ 0.256) and median cost/opt not above it
(≤ 1.117).

Variants (`tools/iterate.py`, kinds `star` and `profile`). Each goes through the same masking, projection,
path check and Stage 0 as every full-panel variant, with a cap of 900 s and 12 GB per region:
- **Centre-star** (`st_cons`, `st_chm13`, `st_long`). Each distinct allele is aligned to one centre by
  abPOA on the two sequences (`-m 0 -b -1`, default two-piece affine scores). Its indels are then
  left-normalised where that leaves the cost unchanged. The pairwise alignments are merged on the centre's
  columns, and the insertions that fall between the same two centre bases are aligned to each other with
  abPOA. The centre is abPOA's consensus of the full panel (`abpoa -r 0`, longest first: the graph
  `poa_abpoa__all` comes from), CHM13's sequence, or the longest allele.
- **Profile aligners**: FAMSA 2.5.2 defaults (`pf_famsa`; it scores with a protein matrix, having no DNA
  mode), FAMSA with twice the gap-open cost (`pf_famsa_go2`), Kalign 3.5.1 `--type dna` (`pf_kalign`), mafft
  v7.526 FFT-NS-2 (`pf_mafft`) and MUSCLE 5.3 `-super5` (`pf_muscle`). Kalign came from Homebrew core,
  FAMSA and MUSCLE from the brewsci/bio Homebrew tap. No full-panel mafft MSAs of these regions were on
  disk: `results/realign_runtime.all.mafft.tsv` covers the L-numbered region set.

Stage 0, 22 regions. b/w = regions where the metric fell / rose by more than 0.005 against
`poa_abpoa__all`; lower is better for both metrics.

| variant | kmerx median / mean | kmerx b/w | cost/opt median / mean | cost/opt b/w | nodes/kb | max align s | gate |
|---|---|---|---|---|---|---|---|
| mc | 0.259 / 0.307 | 11/7 | 1.408 / 1.507 | 4/18 | 215 | - | - |
| `poa_abpoa` (hap32 alone) | 0.175 / 0.237 | 20/1 | 1.057 / 1.092 | 16/2 | 104 | - | sample-dependent |
| `poa_abpoa__all` | 0.286 / 0.320 | - | 1.117 / 1.203 | - | 136 | 17 | - |
| st_chm13 | 0.192 / 0.276 | 16/3 | 1.259 / 1.347 | 8/14 | 79 | 27 | fails cost/opt |
| st_cons | 0.231 / 0.295 | 15/4 | 1.153 / 1.248 | 11/10 | 102 | 40 | fails cost/opt |
| st_long | 0.234 / 0.288 | 13/6 | 1.215 / 1.427 | 6/14 | 88 | 41 | fails cost/opt |
| pf_famsa | 0.287 / 0.304 | 11/6 | **1.066 / 1.104** | 18/1 | 115 | 108 | fails kmerx |
| pf_mafft | 0.345 / 0.331 | 11/9 | 1.262 / 1.493 | 6/15 | 154 | 344 | fails both |
| pf_kalign | 0.403 / 0.407 | 1/21 | 1.342 / 1.366 | 2/18 | 387 | 257 | fails both |
| pf_muscle | the 8 largest regions hit 900 s (in UCLUST clustering); only the smallest, TR773368, was run to the end (162 s); stopped | | | | | | fails |
| pf_famsa_go2 | 2 of 11 regions finished (3 over 12 GB, 4 over 900 s, 2 exited with an error); stopped | | | | | | fails |

- **Stacking indels on one centre removes much of the duplication.** Star cuts median kmerx by 0.05-0.09,
  and with CHM13 as the centre reaches 0.192, close to hap32 alone. The largest gains are where 4m's
  fragmentation was worst:
  - TR770009: 0.297 -> 0.111 (`st_cons`).
  - TR755896: 0.343 -> 0.177.
  - TR773368: 0.096 -> 0.000, the same graph as hap32 alone.

  This confirms 4m: the fragmentation comes from aligning each allele to its best path through hundreds of
  others.
- **Every star raises cost/opt.** Two alleles meet only through the centre, so their induced alignment can
  cost the sum of their distances to it.
  - Insertions in different slots of the centre are never aligned to each other, and the same goes for
    deletions of different centre bases.
  - Worst case, diverse repeats of equal length: in TR764119, the induced pairwise edit cost summed over
    the hap32 pairs is 50,892, against 24,968 for `poa_abpoa__all` (cost/opt 2.67 vs 1.29).
  - The same slot effect adds duplication in some regions: TR755910 kmerx 0.109 -> 0.400.
  - The consensus is the most central centre and costs least. The longest allele almost removes
    insertion slots (8,193 insertion columns over the 22 regions, against 244,348 for the consensus), but
    not the cost.
- **FAMSA is the best full-panel aligner on cost/opt.** It reaches 1.066 against 1.057 for hap32 alone,
  and is better than `poa_abpoa__all` in 18 regions and worse in 1, but its median kmerx does not move.
  - mafft FFT-NS-2 and Kalign are worse than abPOA on both metrics.
  - MUSCLE5 is too slow for 400 or more alleles of a few kb.
  - Raising FAMSA's gap-open cost so that gaps stack exceeds the time and memory caps.
- **No variant passes the gate, so none was called.** The two objectives trade off: a single centre
  stacks indels but loses pairwise optimality, and progressive profile alignment keeps pairwise optimality
  but still fragments.

Next, if the full panel is kept: refine the star by leave-one-out realignment of each allele to the
profile of the rest, which targets exactly the cross-slot cost. Otherwise use 4m's sample-dependent
fallback.

Files: `tools/iterate.py` (`star_align`, `left_normalise_pair`, `star_merge`, `profile_align`); candidates,
full MSAs and Stage 0 under `work/iterate/{candidates,panel,stage0}/<variant>/`.

## 4o. The repeat-unit-aware aligner on the full panel

Question: does `unit_aware` (`tools/realign_units.py`), which aligns whole repeat units, give a
sample-independent full-panel MSA whose hap32 projection is as good as aligning hap32 alone? Same gate as 4n:
median kmerx ≤ 0.256 and median cost/opt ≤ 1.117.

Variants (`tools/iterate.py`, kind `link`):
- `ua32`: the existing hap32 arm, `work/stage4/candidates/unit_aware`.
- `ua_all`: the full panel aligned by `tools/units_panel.py run --panel-root work/stage4/panel`, with its
  existing caps (1,800 s, 12 GB per region), then projected.
- `ua_all_poa`: the same, with the new `--fallback-engine abpoa`. abPOA (`poa_abpoa`) then aligns
  everything that is not a unit array: the flank pieces, regions without a usable motif, and the adequacy
  guard's fallback. The switch is `realign_units.FALLBACK_ENGINES`. Its output goes under the method name
  `unit_aware_poa`.

Every test region's `region.json` carries a motif. 21 regions take the unit path. TR773234 (period 2) falls
back in both variants because its CHM13 array covers 0.37 of the core. The guard ran at TR761129 and
TR764119, and both times the fallback cost more, so the unit MSA was kept.

Stage 0, 22 regions; b/w against `poa_abpoa__all` as in 4n.

| variant | kmerx median / mean | kmerx b/w | cost/opt median / mean | cost/opt b/w | nodes/kb | max align s / peak RSS | gate |
|---|---|---|---|---|---|---|---|
| mc | 0.259 / 0.307 | 11/7 | 1.408 / 1.507 | 4/18 | 215 | - | - |
| `poa_abpoa` (hap32 alone) | 0.175 / 0.237 | 20/1 | 1.057 / 1.092 | 16/2 | 104 | - | sample-dependent |
| `poa_abpoa__all` | 0.286 / 0.320 | - | 1.117 / 1.203 | - | 136 | 17 s | - |
| st_chm13 | 0.192 / 0.276 | 16/3 | 1.259 / 1.347 | 8/14 | 79 | 27 s | fails cost/opt |
| pf_famsa | 0.287 / 0.304 | 11/6 | 1.066 / 1.104 | 18/1 | 115 | 108 s | fails kmerx |
| ua32 (hap32 alone) | 0.231 / 0.266 | 18/3 | 1.041 / 1.062 | 21/0 | 115 | - | sample-dependent |
| **ua_all** | **0.257** / 0.294 | 17/3 | **1.061** / 1.084 | 19/0 | 116 | 274 s / 7.3 GB | fails kmerx by 0.001 |
| ua_all_poa | 0.258 / 0.295 | 16/3 | 1.061 / 1.091 | 18/0 | 116 | 101 s / 1.5 GB | fails kmerx |

- **`ua_all` is the first full-panel variant that beats `poa_abpoa__all` on both metrics in most
  regions.**
  - Its cost/opt is within 0.004 of hap32 alone, like FAMSA's.
  - Its median kmerx (0.2572) misses the gate (0.256) by 0.001, so it was not called. It is also above the
    0.24 that would have justified calling a near miss.
- **Projecting costs the unit-aware aligner much less than it costs abPOA.**
  - From hap32 alone to the full-panel projection, `unit_aware` gains 0.026 kmerx and 0.019 cost/opt.
    abPOA gains 0.111 and 0.060.
  - Most of the remaining kmerx gap to `poa_abpoa` (0.175) is already there in `ua32` (0.231). The unit
    aligner duplicates more k-mers than abPOA even on the hap32 rows alone.
- **The fallback engine hardly matters on repeats.**
  - Apart from TR773234, the two engines differ only in the flanks. Per-region kmerx is within 0.004 and
    cost/opt within 0.012.
  - At TR773234, mafft FFT-NS-i (0.036 / 1.127) beats abPOA, which reproduces `poa_abpoa__all` exactly
    (0.072 / 1.266). That is the whole difference in the means.
  - abPOA is faster and lighter. At TR761129 the guard's mafft fallback took 274 s and 7.3 GB against
    48 s for abPOA.
  - The mafft deletion scattering that hurt `unit_aware` on whole chr20 came from non-repeat fallback
    regions, and this test set has only one.

Next: `unit_aware`'s own kmerx on hap32 is the lever, not the panel. Find which unit columns hold the extra
k-mers in `ua32` against `poa_abpoa`. Or call `ua_all` anyway, since it misses the gate by 0.001 and its
cost/opt is the best of any sample-independent variant.

Files: `tools/units_panel.py` and `tools/realign_units.py` (`--fallback-engine`), `tools/iterate.py`
(`ua32`, `ua_all`, `ua_all_poa`); full MSAs in `work/stage4/panel/unit_aware{,_poa}/`, projections in
`work/stage4/candidates/unit_aware{,_poa}__all/`; runtime rows in `results/realign_runtime.all.units.tsv`.

## 4p. Re-baseline on vg 91d38c802, and the unit-aware full panel on whole chr20

All arms called with the current PR head (vg 91d38c802: positional depth-rate window, nested kappa from
the region's ploidy, depth-term normaliser, code-review fixes), on the existing graphs and mappings.

| chr20 arm | sites (VCF IDs) | SV F1 raw | SV FP / FN raw | SV F1 refined | refined FP / FN | indel F1 | ALL F1 |
|---|---|---|---|---|---|---|---|
| unpatched (mc) | 112,156 | 0.5323 | 440 / 325 | 0.6397 | 347 / 243 | 0.9274 | 0.9721 |
| abPOA on hap32 | 114,217 | **0.5838** | 305 / 321 | **0.6654** | 251 / 252 | 0.9281 | 0.9724 |
| abPOA on the full panel, projected | 117,298 | 0.4769 | 425 / 391 | 0.6139 | 322 / 280 | 0.9180 | 0.9687 |
| unit-aware on the full panel (abPOA fallback), projected | 116,231 | 0.4756 | 438 / 390 | 0.5799 | 367 / 295 | 0.9227 | 0.9706 |

Refined = `truvari refine -u -a mafft`, insensitive to how an allele is split into records.

- The Stage 4 conclusions hold on the current caller: patching with abPOA on hap32 gains +0.051 raw and
  +0.026 refined SV F1; projecting a full-panel alignment loses, raw and refined.
- The unit-aware full panel (4o) has hap32-level stage-0 cost/opt, yet calls worst after refinement:
  stage-0 metrics do not predict calling here.
- `--no-atomize-blocks` does not rescue the full-panel arms (SV FN 381 and 368 against 316 for
  abPOA-hap32), so their losses are not records split by block emission.
- Haplotype ED summed over the 624 patched regions does not track SV F1 either (unpatched 83,037,
  abPOA-hap32 90,618, unit-aware full 85,670, abPOA full 95,414; region counts mixed); a few regions
  dominate the sums.
- What does separate them: sites. The full-panel projections add 4,075-5,142 sites over the unpatched
  graph against 2,061 for abPOA-hap32. Each indel the full-panel alignment cuts at a breakpoint supported
  by an absent allele (4m) becomes a cut point, and so a site boundary, in the projection; more sites
  give the genotyper more junctions at which to stitch recombinant alleles (4i: 87-92% of the error).
- Test set on the new binary (summed called-haplotype ED): mc as a candidate 7,061 (6,840 native on the
  old binary; the renumbering effect shrank from +939 to +221), abPOA hap32 10,388, abPOA full 8,288,
  unit-aware full 8,118 (mafft fallback 8,274).

## 4q. Nearest-neighbour threading: each allele placed against one similar allele

Question: 4m and 4p trace the full-panel projection's loss to POA threading each allele through many
existing allele paths, so that its indels are cut at breakpoints supported by other alleles, and those
breakpoints become extra sites in the projection. Does an MSA in which each distinct allele is placed
relative to ONE similar allele remove the extra sites, and does that recover hap32-level calling? Gate for
calling (site count is now the main gate): a median site count clearly below `poa_abpoa__all` and near
`poa_abpoa`, with kmerx not worse than `poa_abpoa__all`.

Method (`tools/iterate.py`, kind `mst`; `mst_order`, `pair_align`, `place_insert`, `ThreadedMSA`,
`mst_align`; tests in `tools/test_iterate_mst.py`):
- **Distance.** 1 minus the multiset 15-mer Jaccard index: the j-th copy of a k-mer is its own token, so
  copy-number differences count. Exact, over all pairs (at most 33 s per region on chr20).
- **Order.** Prim's minimum spanning tree grown from CHM13's sequence (s0). Each allele joins at its
  nearest tree member.
- **Alignment.** Each allele, in tree order, is aligned pairwise to its parent: abPOA on the two sequences
  (`-m 0 -b -1`, default convex gaps), with indels left-normalised (`left_normalise_pair`, as 4n). A pair
  above 6x10^8 cells, or one that hits the 12 GB cap unbanded, uses abPOA's adaptive band instead. This
  applied only to the three chr20 regions with 30-37 kb alleles, which also have no `poa_abpoa__all`
  candidate.
- **Merge.** The allele inherits the parent's columns.
  - A base aligned to a parent base, match or mismatch, takes that base's column.
  - A deletion is gaps in the parent's columns.
  - The bases inserted between two parent bases go into that slot's existing columns, where the parent
    has gaps, and into new columns. `place_insert` decides by a small global DP: +2 for a base a column
    already holds, -3 for a mismatching one, -1 for a new column, -2.5 for skipping a column inside the
    insertion, free leading and trailing skips, ties to the left. So a recurrent insertion stacks in the
    columns of the first, and an allele that restores bases its parent deleted reuses the deleted bases'
    columns.
- **`mst3`.** Each allele, in the same order, is aligned to its 3 nearest already-aligned alleles by k-mer
  distance. The alignment with the lowest abPOA cost (4 per mismatch, min(4 + 2L, 24 + L) per gap) is kept.
  On the test set 34% of alleles took a parent other than their tree parent.
- Both variants then go through the same projection, `write_gfa`/path check and Stage 0 as every
  full-panel variant (900 s and 12 GB per region).

**Site count.** Stage 0 now reports the top-level snarls of each graph: `vg snarls` on the GFA, read with
`vg view -R`, using vg 91d38c802 (`iterate.snarl_counts`, cached as `stage0/<variant>/<id>.snarls.json`).
It reproduces `work/iterate/snarl_stats.json` exactly for mc, `poa_abpoa` and `poa_abpoa__all`.

Stage 0, 22 regions. b/w = regions where the metric fell / rose against `poa_abpoa__all` (by more than
0.005 for kmerx and cost/opt, by any amount for sites); lower is better for all three.

| variant | kmerx median / mean | kmerx b/w | cost/opt median / mean | cost/opt b/w | sites median / mean (total) | sites b/w | nodes/kb | max align s | gate |
|---|---|---|---|---|---|---|---|---|---|
| mc | 0.259 / 0.307 | 11/7 | 1.408 / 1.507 | 4/18 | 21 / 25.3 (557) | 20/0 | 215 | - | - |
| `poa_abpoa` (hap32 alone) | 0.175 / 0.237 | 20/1 | 1.057 / 1.092 | 16/2 | 30 / 64.1 (1,411) | 18/3 | 104 | - | sample-dependent |
| `poa_abpoa__all` | 0.286 / 0.320 | - | 1.117 / 1.203 | - | 47.5 / 74.9 (1,648) | - | 136 | 17 | - |
| st_chm13 (4n) | 0.192 / 0.276 | 16/3 | 1.259 / 1.347 | 8/14 | 14 / 26.3 (579) | 22/0 | 79 | 27 | passes |
| pf_famsa (4n) | 0.287 / 0.304 | 11/6 | 1.066 / 1.104 | 18/1 | 27.5 / 55.7 (1,226) | 20/1 | 115 | 108 | fails kmerx by 0.001 |
| ua_all_poa (4o) | 0.258 / 0.295 | 16/3 | 1.061 / 1.091 | 18/0 | 29 / 61.3 (1,349) | 16/3 | 116 | 101 | passes |
| **mst** | 0.282 / 0.285 | 17/4 | 1.238 / 1.459 | 4/15 | **23.5 / 34.4** (756) | 20/2 | 97 | 40 | passes |
| **mst3** | 0.245 / 0.271 | 17/3 | 1.260 / 1.370 | 4/16 | **18 / 35.7** (786) | 21/1 | 103 | 96 | passes |

- **Threading on one parent removes the extra sites.**
  - Median sites fall from 47.5 to 23.5 (`mst`) and 18 (`mst3`), below hap32 alone (30). Total sites over
    the 22 regions fall to about half of hap32 alone.
  - On whole chr20 the same holds: summed over the 636 region graphs (`work/iterate/snarls_chr20.py`),
    `mst` has 11,424 top-level sites, against 14,212 for hap32 alone, 18,334 for `poa_abpoa__all` (633
    regions), 16,289 for unit-aware full and 9,385 for mc. All snarls, nested ones included: 25,866,
    against 26,557, 34,099, 32,282 and 23,865.
- **kmerx improves in most regions:** 17 better, 4 worse for `mst`. `mst3` also lowers the median, to
  0.245.
- **cost/opt is worse, as for the star (4n).** 1.24-1.26 against 1.117. Two alleles in different
  subtrees meet only through their tree path, so their induced alignment can cost the sum of the steps.
  Re-choosing the parent among the 3 nearest (`mst3`) does not change this.
- **The test-set site count misranked unit-aware full.** On the test set its sites are at hap32 level
  (median 29, total 1,349 against 1,411). On whole chr20 it has 2,077 more top-level sites than hap32
  alone, and calling there added twice as many VCF sites (4p). The chr20 region-graph totals rank all four
  earlier arms in the same order as their VCF site counts.

**Test set, called** (vg 91d38c802, 3 replicates; summed median called-haplotype ED, lower is better; b/w
by the replicate rule against the re-baseline arms):

| variant | summed ED | b/w vs mc | b/w vs `poa_abpoa` | b/w vs `poa_abpoa__all` |
|---|---|---|---|---|
| mc (`rb_mc_link`) | 7,061 | - | 5/10 | 10/7 |
| `poa_abpoa` (`rb_poa_abpoa`) | 10,388 | 10/5 | - | 10/3 |
| `poa_abpoa__all` (`rb_poa_abpoa__all`) | 8,288 | 7/10 | 3/10 | - |
| ua_all_poa | 8,118 | 9/7 | 4/8 | 9/6 |
| st_chm13 | 10,778 | 9/7 | 6/6 | 10/7 |
| **mst** | **7,750** | 9/8 | 5/8 | 10/6 |
| mst3 | 8,065 | 10/5 | 4/8 | 9/4 |

`mst` has the lowest summed ED of any full-panel arm, but as 4p found, test-set ED does not rank the
whole-chr20 SV F1. hap32 alone has the worst ED of the baselines and the best SV F1.

**Whole chr20** (all 636 eligible repeats aligned with `mst`, written to `work/stage4/candidates/mst__all/`;
`patch_contig.py patch/map/call/score` as 4p, on vg 91d38c802; `truvari refine -u -a mafft` as
`refine4.sh`; scripts `work/iterate/mst_chr20_align.sh`, `mst_chr20.sh`):

| chr20 arm | sites (VCF IDs) | SV F1 raw | SV FP / FN raw | SV F1 refined | refined FP / FN | indel F1 | ALL F1 |
|---|---|---|---|---|---|---|---|
| unpatched (mc) | 112,156 | 0.5323 | 440 / 325 | 0.6397 | 347 / 243 | 0.9274 | 0.9721 |
| abPOA on hap32 | 114,217 | **0.5838** | 305 / 321 | **0.6654** | 251 / 252 | 0.9281 | 0.9724 |
| abPOA on the full panel, projected | 117,298 | 0.4769 | 425 / 391 | 0.6139 | 322 / 280 | 0.9180 | 0.9687 |
| unit-aware on the full panel, projected | 116,231 | 0.4756 | 438 / 390 | 0.5799 | 367 / 295 | 0.9227 | 0.9706 |
| **mst on the full panel, projected** | **113,616** | 0.5446 | 447 / 310 | 0.6438 | 358 / **235** | 0.9248 | 0.9714 |
| mst3 on the full panel, projected | 113,812 | 0.5449 | 443 / 310 | 0.6349 | 369 / **235** | 0.9277 | 0.9723 |

`mst3`'s three largest regions (TR762710, TR762711, TR762722; 18-22 kb alleles) went over the 900 s cap
and use the `mst` graph (`work/stage4/candidates/mst3__all/fallback_mst.txt`; `work/iterate/mst3_chain.sh`).

- **Threading recovers most of the full-panel loss, but not hap32-level calling.** Refined SV F1 is
  0.644 for `mst`, against 0.614 for abPOA on the full panel and 0.580 for unit-aware. That is above
  unpatched (0.640) but 0.022 below abPOA on hap32 (0.665). `mst3` is lower refined (0.635) and equal raw.
  It has the best indel and ALL F1 of any full-panel arm (0.9277 / 0.9723, against 0.9281 / 0.9724 for
  hap32).
- **Fewer sites did not buy fewer false positives.**
  - Both `mst` arms call fewer sites than hap32 patching (113,616 and 113,812, against 114,217), and their
    region graphs have 20% fewer top-level sites.
  - Their gain is in recall. Raw SV TP-base is 455, the highest of any arm (444 for hap32). Refined FN is
    235, the lowest (252 for hap32, 243 unpatched).
  - Their loss is in precision. Raw SV FP inside the patched spans is 392 (`mst`) and 388 (`mst3`),
    against 251 for hap32 and 385 unpatched.
  - The FP excess over hap32 is spread out. `mst` has more FPs in 85 regions and fewer in 36, and the
    10 worst regions hold 75 of the 141. Two of those 10 are banded-fallback regions: TR768271 (+14) and
    TR773472 (+8).
- **So the site count is not sufficient on its own.** The extra junctions of the abPOA projection fit its
  lost recall: its refined FN is 280, and `mst`, with fewer sites, gets 235. They do not account for hap32's
  precision. The likely cause of the remaining FPs is the cost/opt that threading gives up (1.24 against
  1.06 for hap32 alone). Alleles in different subtrees are aligned only through their tree path, so a
  real SV allele can be spelled by columns that a better alignment would share with its neighbours.

Recommendation: keep abPOA on hap32 as the patch method. `mst` is the best sample-independent full-panel
arm so far, with refined SV F1 0.644 against 0.614. One-parent threading fixes the junction and recall
problem, so the full panel's loss is not inherent to using the full panel. To close the remaining 0.02, gate
on cost/opt as well as sites. Keep the tree threading for column placement, then refine pairwise
optimality without re-introducing many-allele junctions. Two ways to try: realign each allele to the
profile of its tree neighbourhood (leave-one-out, 4n's next step), or thread on the hap32 rows first and
add the rest of the panel as children of their nearest hap32 allele.

Files: `tools/iterate.py` (kinds `mst`, `mst3`; `snarl_counts`, `snarls_one`; the summary now carries
mean kmerx and cost/opt and median and mean sites), `tools/test_iterate_mst.py`,
`work/iterate/snarls_chr20.py` (and its `snarls_chr20.json`), `work/iterate/mst_chr20.sh`; candidates in
`work/iterate/candidates/mst{,3}/` and `work/stage4/candidates/mst{,3}__all/`, full MSAs in
`work/iterate/panel/mst{,3}/`, chr20 arms in `work/stage4/chr20/patched_mst{,3}/` and
`score_rb/rb_patched_mst{,3}/`. The three `mst` regions over 6x10^8 cells were rebuilt with the banded
fallback (`iterate.py build mst --regions TR768271,TR768619,TR773472 --force`) before patching.

## 4r. Backbone plus threading: an abPOA backbone of k-medoid representatives

Question: `mst` (4q) has the best recall of any full-panel arm but loses precision. Its cost/opt is 1.238,
against 1.117 for `poa_abpoa__all` and 1.057 for hap32 alone, because alleles in different branches of the
tree meet only through the tree. Does a two-tier alignment fix the cross-branch alignment without bringing back the
many-allele junctions? Tier 1 is abPOA on a small, diverse set of representatives, where abPOA does well.
Tier 2 threads the rest onto that backbone. Gate for calling: median sites near `mst` (≤ ~24) and median
cost/opt clearly below 1.238 (aim ≤ 1.15).

Method (`tools/iterate.py`, kind `bbt`; `kmedoids`, `bbt_plan`, `bbt_align`, `ThreadedMSA.from_msa`; tests
in `tools/test_iterate_mst.py`). The method is sample-independent: it uses no HG002 data and does not
consult the hap32 rows.
- **Clusters.** The region's distinct full-panel sequences (panel weight > 0) are split into K groups by
  k-medoids on `mst`'s distance (1 minus the multiset 15-mer Jaccard index). CHM13's sequence is always a
  fixed medoid, and so is GRCh38's when it differs. Farthest-first traversal seeds the free medoids, and
  FasterPAM's eager swaps then refine them. Ties are broken by sequence length and then the sequence
  itself, never by hap32 membership.
- **Backbone.** abPOA aligns the K medoids with `abpoa -m 0 -r 1` defaults, longest first, as `poa_abpoa`.
- **Threading.** Every other panel sequence is aligned pairwise to a parent (`pair_align`, as `mst`). It is
  added on the parent's columns with `mst`'s column-inheritance merge (`ThreadedMSA.add`). Insertions go into
  the slot's existing gap columns where they fit, so recurrent insertions stack, including insertions
  against backbone columns that another medoid fills. Two parent rules:
  - `bbt32`, `bbt64`: the nearest already-placed member of the cluster (Prim's tree grown from the medoid).
  - `bbt32m`, `bbt64m`: the cluster's medoid itself.

  The hap32-only recombinant paths (panel weight 0; one test region has any) go last, on the nearest
  medoid's cluster, so that no panel row's placement depends on them.
- Then the same projection, `write_gfa`/path check and Stage 0 as every full-panel variant (900 s and 12 GB
  per region). The longest build took 38 s. Clustering took at most 12 s and the backbone at most 2.4 s.

Stage 0, 22 regions. b/w = regions where the metric fell / rose against `mst` (by more than 0.005 for kmerx and
cost/opt, by any amount for sites); lower is better for all three.

| variant | kmerx median / mean | kmerx b/w | cost/opt median / mean | cost/opt b/w | sites median / mean (total) | sites b/w | nodes/kb | gate |
|---|---|---|---|---|---|---|---|---|
| `poa_abpoa` (hap32 alone) | 0.175 / 0.237 | 18/2 | 1.057 / 1.092 | 18/1 | 30 / 64.1 (1,411) | 6/14 | 104 | sample-dependent |
| `poa_abpoa__all` | 0.286 / 0.320 | 4/17 | 1.117 / 1.203 | 15/4 | 47.5 / 74.9 (1,648) | 2/20 | 136 | - |
| `mst` | 0.282 / 0.285 | - | 1.238 / 1.459 | - | 23.5 / 34.4 (756) | - | 97 | - |
| bbt32 | 0.269 / 0.300 | 8/11 | 1.195 / 1.255 | 13/7 | 29.5 / 54.0 (1,188) | 3/18 | 118 | fails both |
| bbt64 | 0.283 / 0.312 | 6/14 | 1.158 / 1.247 | 15/4 | 30.5 / 58.4 (1,285) | 2/18 | 129 | fails both (cost/opt by 0.008) |
| bbt32m | 0.280 / 0.296 | 8/12 | 1.195 / 1.236 | 14/5 | 29.5 / 56.0 (1,232) | 3/18 | 119 | fails both |
| **bbt64m** | 0.286 / 0.308 | 7/14 | **1.147** / 1.234 | 16/3 | 32.5 / 61.5 (1,352) | 2/19 | 127 | fails sites |

- **The backbone improves cost/opt, and K is a dial between `mst` and full-panel abPOA.**
  - Median cost/opt falls from 1.238 to 1.195 at K = 32 and to 1.147-1.158 at K = 64. That is better than
    `mst` in 13-16 regions and worse in 3-7.
  - Threading on the medoid is slightly better than the in-cluster tree at K = 64. Each member's alignment
    to its medoid is then one optimal pairwise step, not a composition of several.
  - Means barely move (1.23-1.26 against 1.459 for `mst`), because a few regions dominate them. In TR761129,
    `mst` is at 3.79 and bbt at 2.20-2.45, against 1.34 for hap32 alone.
- **The sites come back with the backbone.**
  - Median sites rise from 23.5 to 29.5-32.5, at the level of hap32 alone (30), and rise in 18-19 of the 22
    regions.
  - Total sites rise from 756 to 1,188-1,352, between hap32 alone (1,411) and `mst`.
  - This is 4m's mechanism, at the scale of K alleles instead of 460. The POA cuts each representative's
    indels at junctions that other representatives support, and the projection keeps them.
  - Larger K raises sites and lowers cost/opt. No K reaches both of `mst`'s sites and `poa_abpoa__all`'s
    cost/opt.
- **Where the remaining cost/opt is** (`work/iterate/bbt_decomp.py`; graph-implied unit cost of every hap32
  pair, through a new opt-in per-pair dump in `evaluate.py`, `VNTR_EVAL_PAIRS_OUT`).
  - Pairs whose two alleles share a cluster carry under 1% of the optimal cost, so the gate is decided by
    pairs across clusters.
  - Pooled over the test set, cross-cluster pairs reach cost/opt 1.20-1.23 for the four bbt arms, against
    1.41 for `mst`, 1.17 for full-panel abPOA and 1.09 for abPOA on hap32 alone.
  - So the backbone recovers about three quarters of `mst`'s cross-branch loss. The rest is composition: a
    member reaches another cluster only through its own medoid's backbone row, and a pairwise alignment
    composed with a backbone alignment is not an optimal alignment.
- **No variant passes, so none was called.**
  - `bbt64m` meets the cost/opt aim (1.147) but fails the sites gate: its median of 32.5 is at the level of
    hap32 alone, not of `mst`.
  - `bbt64` misses on both counts.

Recommendation: keep abPOA on hap32 as the patch method, and `mst` as the best sample-independent arm
(refined SV F1 0.644). Over this family, sites and cost/opt trade against each other through a single
quantity: how many alleles a sequence is aligned against at once. `mst` uses one, bbt uses K, full-panel
abPOA uses all of them. To test 4q's precision hypothesis directly, call `bbt64m` anyway. Its median
cost/opt (1.147) is well below `mst`'s (1.238) and near `poa_abpoa__all`'s (1.117), and its sites are at the
level of hap32 alone, the best-calling arm. If its refined SV F1 lands near hap32's 0.665, sites at that
level are not the problem, and the gate should be on cost/opt. Otherwise the next lever is to cut the composition loss without a graph of many alleles:
realign each member to the backbone rows of its own and its nearest other cluster (two to three alleles),
keeping `place_insert`'s column reuse.

Files: `tools/iterate.py` (kind `bbt`: `kmedoids`, `_prim_from`, `bbt_plan`, `bbt_align`,
`ThreadedMSA.from_msa`; variants `bbt32`, `bbt64`, `bbt32m`, `bbt64m`), `tools/test_iterate_mst.py`
(`TestBackbone`), `tools/evaluate.py` (`VNTR_EVAL_PAIRS_OUT`), `work/iterate/bbt_decomp.py` and
`bbt_decomp.json`; candidates in `work/iterate/candidates/bbt*/`, full MSAs in `work/iterate/panel/bbt*/`,
Stage 0 in `work/iterate/stage0/bbt*/`.
