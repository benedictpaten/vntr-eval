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
