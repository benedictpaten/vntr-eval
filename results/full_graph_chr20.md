# The full chr20 graph with every tandem repeat realigned, tested end to end

2026-10-03. The pipeline is `pgrealign` (docs/pipeline-design.md). The drivers are work/full/overnight.sh,
work/full/chr20/redo_realigned.sh, test_chain.sh and test_arm.sh. The per-region bootstrap is
work/full/chr20/test/region_boot.py.

## What was built

**Input.** chr20 of the HG002-free eval graph, `hprc-v2.1-mc-chm13-eval.gref.gbz`.
- 458 haplotypes: 456 HPRC, plus CHM13 and GRCh38.
- The gref_CHM13 cover was removed.
- 4,618,817 nodes, 99.15 Mb.

**The pipeline, step by step:**

| step | result | time |
|---|---|---|
| catalogue (RepeatMasker + trfind) | 18,098 chr20 targets | (genome-wide, 26 s) |
| regions (anchors from the full graph's snarls) | 14,511 regions; 21 too long (> 100 kb) | 15 s |
| extract (one streaming pass) | all 2,173 walks; 0 complex; 248 internal-fragment candidates | 8.6 min |
| realign (medoid star, every region whose alleles differ) | 14,492 ok, 19 too big | 36 min at 8 jobs (4.1 CPU-h) |
| replace (dense ids) | 845,062 nodes out, 1,011,456 in; 6.6 M runs rewritten | 20 min, 13.2 GB |
| verify | every one of 2,173 paths spells its original sequence (28.7 Gbp); metadata identical | 12 min |

The realigned graph has 4,785,211 nodes and 101.54 Mb.

**It holds 2.4% more sequence, and 18% more inside the replaced regions.** The growth sits in the most
polymorphic regions, because the star aligns every allele only to the centre:

| distinct alleles | regions | change in sequence |
|---|---|---|
| ≤ 5 | 628 | +30% |
| 6-20 | 9,500 | +0.2% |
| 21-100 | 3,908 | +4.7% |
| > 100 | 456 | +112% |

## How it was tested

The same chain ran on both graphs, with one binary (vg 91d38c802) throughout.

1. **Haplotype sampling.**
   - The indexes were built with `vg autoindex -w sampling`.
   - Sampling was `vg haplotypes --num-haplotypes 32 --include-reference --set-reference CHM13
     --set-reference GRCh38`, without diploid sampling. That is how the production hap32 graph appears
     to have been made: 32 recombinants plus CHM13 and GRCh38.
   - The k-mer counts are HG002's 29-mers from the chr20 short reads below: jellyfish, converted by
     tools/c/jf2kff. On vg's own test data, jf2kff samples a byte-identical graph to KMC's counts.
   - In both sampled graphs the CHM13 path is byte-identical to chr20.fa.
2. **Short reads.** The 13,278,235 HG002 reads that production placed on chr20, mapped with giraffe as in
   stage 4: paired, fragment length 401.7 ± 166.4, then the singles.
3. **ONT reads.** The 85,373 reads (2.86 Gbp, 43x) that the genome-wide ONT build placed on chr20. They
   were rebuilt as FASTQ from their alignments to E821 with `vg convert -F` against the whole E821 graph:
   219 of them lie on parts of chr20 the chr20 subgraph lacks. They were mapped with giraffe `-b r10`.
4. **Calling.** Through a per-graph GAF-Base. Short reads used the production flags (`--mismap-max 0.95
   --read-likelihood --phased`). ONT used `--preset ont`.
5. **Scoring.** bench_wgs, against T2T-Q100 v1.1: aardvark for small variants, truvari for SVs, and
   truvari refine.

## Results (chr20, HG002)

| | ALL | SNV | indel | raw SV F1 (FP / FN) | refined SV F1 (FP / FN) |
|---|---|---|---|---|---|
| original, short reads | 0.9724 | 0.9853 | 0.9277 | 0.5435 (411 / 323) | 0.6342 (328 / 261) |
| **realigned, short reads** | 0.9723 | 0.9855 | 0.9262 | **0.6144** (247 / 315) | **0.6741** (202 / 273) |
| original, ONT | 0.9712 | 0.9876 | 0.9133 | 0.5475 (469 / 294) | 0.6771 (346 / 201) |
| **realigned, ONT** | 0.9715 | 0.9879 | 0.9128 | **0.6452** (255 / 278) | **0.7611** (171 / 188) |

**Paired per-region bootstrap.** These are SV FP+FN inside the realigned regions (span ± 100 bp),
realigned minus original:

| | change [95% CI] | regions better / worse | outside the regions |
|---|---|---|---|
| short, raw | −170 [−247, −98] | 138 / 70 | 26 → 24 |
| short, refined | −112 [−186, −41] | 106 / 68 | 21 → 19 |
| ONT, raw | −228 [−307, −156] | 122 / 51 | 24 → 22 |
| ONT, refined | −188 [−278, −106] | 83 / 40 | 5 → 5 |

- **Realigning every tandem repeat in the full graph improves SV calling for both read types.** How much
  of that is real at the level of haplotypes is below.
  - Short reads: refined SV F1 +0.040, raw +0.071. False positives fall by 40% (411 → 247 raw).
  - ONT: refined SV F1 +0.084, raw +0.098.
- **The gain is broad.** The ten most-improved regions hold 41-57% of it, and calls outside the regions
  do not change.
- **Small variants are unchanged.** The short-read indel F1 cost (0.9277 → 0.9262) is examined below.
- **The consistent losers** are TR_chr20_64970081 (+11 short, +9 ONT; 419 alleles, its sequence grew
  12.6 kb → 145 kb) and TR_chr20_18259720 (+8, +7).
- **Comparison with stage 4.** The hap32-projected medoid star (4x) gave +0.060 refined over its own
  re-mapped control. That was a different chain: the 34 production haplotypes, no re-sampling.

**Cost.** Calls take as long as on the original graph: short 296 s against 349 s, ONT 794 s against
857 s, at 8 threads. The comparable timings are below.

### What the SV gain is made of

A record-level SV error can be pure representation: the called haplotype matches the truth haplotype,
but its records are cut differently. To separate the two, each error locus was scored on the
haplotypes themselves, using the nesting-aware `score_haplotypes`. A locus counts as representation
when the called haplotype is within 10 bp edit distance of the truth.
- 514 loci were bootstrapped.
- 4 regions inside one 474 kb enclosing record were excluded.

| refined SV F1 gain | all errors [95% CI] | representation removed [95% CI] |
|---|---|---|
| short reads | +0.040 [−0.002, +0.081] | +0.019 [−0.030, +0.064] |
| ONT | +0.084 [+0.033, +0.135] | +0.068 [+0.013, +0.125] |

- **The ONT gain is real at the level of haplotypes.**
  - ONT haplotypes come within 10 bp of the truth in 36 loci on the realigned graph against 17 on the
    original (p = 0.013).
  - The gain holds at every threshold from 0 to 50 bp: +0.064 to +0.087.
- **About half the short-read gain is representation.**
  - The realigned graph's short-read haplotypes are no closer to the truth: 95 loci better, 93 worse.
  - The remainder cannot be told from zero.
  - The per-region bootstrap above (−112 [−186, −41]) counts records, so it includes the
    representation gain.
- **What is left.**
  - 94-97% of the remaining SV errors are inside realigned regions, and 23 loci hold half of them.
  - About half (45-52%) cannot be fixed by realignment: HG002's allele is novel to the panel, or one of
    its haplotypes was not sampled.
  - **TR_chr20_64970081 is a real regression.** The call is 616 bp where the truth is 4,669 bp.
  - **TR_chr20_18259720 is pure representation.** The realigned ONT haplotypes are exactly HG002's, but
    the records score 36 FP and 4 FN.

### The short-read indel cost

- **About half of it is counting.**
  - The realigned short call writes 429 fewer TP indel records for the same truth, and that alone costs
    0.0008 F1.
  - Base-pair indel F1 rises on both read types: short 0.8508 → 0.8662, ONT 0.8987 → 0.9148.
- **The rest is 62 regions (0.4%) where the star forced a non-homologous allele into shared columns.**
  - There the short call adds 125 indel and 85 SNV FPs. Everywhere else FPs fall (−108 indel, −63 SNV).
  - TR_chr20_62032077 alone adds 101 indel and 60 SNV FPs.
  - Scoring those 62 regions as on the original graph gives indel F1 0.9291, above the original's
    0.9277.
  - ONT is unaffected there (−21 FP).
- **Candidate fix: a forced-homology guard in the star.**
  - After each allele is aligned to the centre, scan 40-column windows for more than 30% mismatches.
  - In a flagged window, put the allele's bases in a private insertion slot instead of sharing columns,
    so the stretch becomes one bubble rather than many single-column ones.

### Why the graph grew 2.4%

The replaced regions grew from 13.11 Mb to 15.50 Mb. Their floor, the longest allele of each region, is
12.20 Mb.

| bp above the floor | |
|---|---|
| cross-slot scatter | 2,957,359 (89.6%) |
| substitution nodes | 223,624 |
| insertions left unmerged within a slot | 75,914 |
| backbone slack | 44,166 |
| **realigned graph** | **3,301,063** |
| original graph | 911,225 |

- **Cross-slot scatter.**
  - The star puts each allele's copy-number expansion in its own insertion slot, and never aligns two
    slots to each other. Within a slot, abPOA merges almost perfectly.
  - Left normalisation cannot gather the slots, because the repeat units are not exact copies.
  - In TR_chr20_1893502, alleles' largest blocks sit in 187 different slots.
- **The medoid's 64-comparator cap moves the centre to shorter alleles** (11 of 14 regions checked).
  TR_chr20_64970081's centre is 1,520 bp instead of 3,448. Removing the cap cuts the expansion only
  partly.
- **Expansion goes with the SV gain, not against it.** Regions whose graph more than doubled carry −85
  of the −112 short-read and −114 of the −188 ONT refined FP+FN change. A growth cap that falls back to
  the original subgraph would hand back most of the gain, so it is rejected.
- **The HG002-sampled graph grew only 69 kb (+0.09%).**
- **TR_chr20_66205689 is the q-arm telomere.** 410 of its 415 walks end inside it, and 142 kb of its
  163 kb is used only by fragments aligned to the centre.
- **Candidate fix: centre the star on the longest allele when the medoid star expands.**
  - On the 40 most-expanded regions this cuts sequence from 1.595 Mb to 0.278 Mb.
  - Against optimal pairwise alignment, precision is unchanged and pair F1 falls 0.021.
  - Re-centring the 66 regions whose graph more than doubled would remove 1.9-2.1 of the 2.39 Mb.
  - Whether it keeps the SV gain needs one more end-to-end run.

## Deliverables: ONT anchors, mosaics and phasing

`work/full/chr20/deliver.sh` makes the collaborator bundle in `work/full/chr20/collaborator_HG002_chr20/`.
Its `README.md` has the commands. The bundle holds:
- the ONT call with `--anchors-out --mosaic-out`;
- the short-read call with `--mosaic-out`.

`--anchors-out` turns on the off-reference descent, which changes the ONT phasing, so the ONT call was
re-run and re-scored. The short-read call is the test chain's.

| | ALL | SNV | indel | raw SV F1 | refined SV F1 | switch error (switches / pairs) |
|---|---|---|---|---|---|---|
| original, ONT + anchors | 0.9712 | 0.9875 | 0.9132 | 0.5482 | 0.6778 | 0.561% (335 / 59,676) |
| **realigned, ONT + anchors** | 0.9716 | 0.9879 | 0.9129 | **0.6452** | **0.7611** | 0.590% (358 / 60,695) |
| original, short | 0.9724 | 0.9853 | 0.9277 | 0.5435 | 0.6342 | 3.03% (1,789 / 59,112) |
| realigned, short | 0.9723 | 0.9855 | 0.9262 | 0.6144 | 0.6741 | 3.20% (1,919 / 59,951) |

- Switch error is from whatshap compare against the T2T-Q100 v1.1 small variants. Every call is one
  chromosome-length phase block.
- **Phasing is slightly worse on the realigned graph:** +23 switch-error pairs on ONT and +130 on short
  reads, over about 1,000 more het pairs. Decomposed with `scripts/switch_bed.py` (vg-call-eval), that is
  31 → 38 true switches and 152 → 160 flips on ONT, and 1,171 → 1,249 switches and 309 → 335 flips on
  short reads.
- **The ONT difference is churn, not a shift.**
  - 19 of the original's 31 switches are gone and 25 new ones appear.
  - 8 of those 44 changed switches are in TR_chr20_26343418. It was too long to realign and is the same
    in both graphs, so they come from sampling and mapping.
  - 97-100% of ONT switches lie within 1 kb of a tandem-repeat region, against 63% of phased het sites.
  - Most come in pairs inside one region (22 of 31, 28 of 38): the repeat's own hets phased against its
    flanks, a local error rather than a long-range break.
- **Short-read switches are no more concentrated in the repeats on the realigned graph:** 81% within
  1 kb of a region, against 80% on the original.
- `work/full/chr20/investigate/switches/where.py` does this.
- **Anchors.** Realigned: 591,327 anchors from 16,620,105 read placements, with all 18,332,729 pins
  verified. Original: 589,036.

## What went wrong on the way, and the fixes

1. **Sparse node ids slowed calling about 8x.** The first build gave new nodes ids above the genome-wide
   maximum.
   - vg call fetches reads in 16,384-id windows. A site straddling two windows is fetched on its own,
     uncached, and with `--gbz-base` set to a GBZ each fetch reloads the graph.
   - The short-read call took 54 min, and the ONT call was stopped after more than 2 h.
   - `replace --id-mode dense`, now the default, puts each region's new nodes right after its left
     anchor.
2. **The ONT read rebuild against the E821 chr20 subgraph crashed** on 219 reads. It now uses the whole
   E821 graph.

## Caveats

- **Only chr20, and contig-local mapping.** The numbers are comparable between the two graphs, not with
  the genome-wide production figures.
- **Each graph samples its own 32 haplotypes,** so a difference mixes graph quality with what was
  sampled.
- **Dense ids renumber the contig.** Contigs need a global renumbering before they are merged into one
  graph.
- **chr6, the held-out contig, has not been run.** It waits for disk: its patched GFA is about 22 GB.
