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

**Cost.** Calling takes about as long on the realigned graph as on the original. See "Run time and
memory" below.

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
- **Switch errors are of three sizes.**
  - whatshap's switch-error count is every adjacent het pair whose relative phase is wrong. Its
    switch/flip decomposition splits that into **flips** (one het phased wrong, two errors each) and
    **true switches**.
  - True switches are of two kinds:
    - two that bound a **short mis-phased stretch** of 2-20 hets;
    - **long-range switches**, after which the phase stays wrong until the next long-range switch.
  - On ONT the split is clean: no short stretch has more than 13 hets, and no phase segment has fewer
    than 29.
  - The decomposition is asserted equal to whatshap's own `all_switchflips` for all four calls.

| | switch errors (whatshap) | flips | true switches | in short stretches (stretches) | long-range | mis-phased bp |
|---|---|---|---|---|---|---|
| original, ONT | 335 (0.561%) | 152 | 31 | 24 (15) | 7 | 25.8 Mb |
| realigned, ONT | 358 (0.590%) | 160 | 38 | 33 (20) | **5** | 19.1 Mb |
| original, short | 1,789 (3.03%) | 309 | 1,171 | 796 (519) | 375 | 31.8 Mb |
| realigned, short | 1,919 (3.20%) | 335 | 1,249 | 878 (588) | 371 | 29.5 Mb |

- **Flips are mostly 1 bp indels and are not a repeat effect.**
  - On ONT, 79% (original) and 81% (realigned) of flipped hets are 1 bp indels, though 1 bp indels are
    only 7% of the assessed hets.
  - Flips lie within 1 kb of a tandem-repeat region at the background rate (61-65% against 63% of hets).
  - Short-read flips are mostly SNVs (55-62%), with 1 bp indels 3x enriched.
- **Short mis-phased stretches are a repeat's own hets phased against its flanks.**
  - On ONT they are small: a median of 3 hets over 588 bp (original) or 762 bp (realigned).
  - 14 of 15 and 20 of 20 sit within 1 kb of a tandem-repeat region.
  - The realigned graph has 5 more.
- **Long-range switches fall from 7 to 5 on ONT.**
  - Each sits in a gap between assessed hets, of 350 bp to 1.26 Mb. Five of the original's seven gaps
    exceed 80 kb, and two of the realigned graph's five.
  - Two of the original's, and one of the realigned graph's, are in regions too long to realign, which
    are the same in both graphs.
  - Only the switch at 36,099,955 is shared.
  - Mis-phased bp (the smaller of the two alternating halves) falls 25.8 → 19.1 Mb. With five to seven
    switches that mostly reflects where they land.
- **Short-read long-range switches are unchanged (375 → 371).** The extra short-read errors are short
  stretches (+69) and flips (+26). Short-read mis-phased bp is 45-48%, near its 50% ceiling, so it says
  nothing.
- **In sum, the realigned graph adds local phasing errors inside repeats and removes two long-range ONT
  switches.**
- Two scripts in `work/full/chr20/investigate/switches/` do this:
  - `classes.py` reads whatshap `--longest-block-tsv`, which is 0-based, so positions are +1;
  - `where.py` reads the switch BED.
- **Anchors.** Realigned: 591,327 anchors from 16,620,105 read placements, with all 18,332,729 pins
  verified. Original: 589,036.

## Run time and memory

### Building the realigned graph

pgrealign on the full chr20 graph (458 haplotypes), 8 threads:

| step | wall | CPU | peak RSS |
|---|---|---|---|
| prepare: vg chunk, drop gref, snarls, reference index | 199 s | 243 s | 14.7 GB |
| extract: one streaming pass | 517 s | 2,905 s | 2.6 GB |
| realign: 14,511 regions, medoid star | 2,137 s | 12,702 s | 8.5 GB |
| replace (dense ids) + vg gbwt -G | 1,203 s | 4,349 s | 13.2 GB |

### Preparing each arm

Per arm, 8 threads; original and realigned graph are within 5% of each other:

| step | wall | peak RSS |
|---|---|---|
| sampling index (vg autoindex -w sampling) | 1,115-1,157 s | 4.8 GB |
| vg haplotypes, 32 haplotypes | 175 s | 4.1 GB |
| giraffe indexes | 75 s | 12.3 GB |
| giraffe, short reads | 530 s | 2.2 GB |
| giraffe, ONT (-b r10) | 871-919 s | 20.9 GB |

### The calls

- **Setup.** `-t 5`, with the graph served from a GBZ-Base database (`gbz-base construct`), as in
  the tier-2 docs.
- **The earlier timings were not comparable.** The test chain passed the GBZ itself as `--gbz-base`,
  at 8 threads, so every read fetch reloaded the graph.
- **Genotypes do not depend on the setup.** Every call below was checked to give the same genotypes
  as the test chain or the deliverable.

Each cell is wall / CPU / peak RSS:

| call | graph | vg 91d38c802 (the deliverables) | vg 8b993a339 (the fix below) |
|---|---|---|---|
| short reads, with mosaic | original | 242 s / 898 s / 8.5 GB | 191 s / 695 s / 10.6 GB |
| | realigned | 256 s / 869 s / 8.5 GB | 218 s / 683 s / 8.6 GB |
| ONT, `--preset ont`, with mosaic | original | 604 s / 2,383 s / 9.1 GB | 484 s / 1,581 s / 10.2 GB |
| | realigned | 750 s / 2,477 s / 7.9 GB | 542 s / 1,564 s / 8.9 GB |
| ONT with `--anchors-out` | original | 1,408 s / 3,468 s / 12.7 GB | 1,030 s / 2,265 s / 11.4 GB |
| | realigned | 1,617 s / 3,624 s / 11.1 GB | 1,064 s / 2,169 s / 12.1 GB |

**Against the earlier reference runs** (vg-call-eval):

- **Short reads match the tier-2 run.** The tier-2 `readlik` chr20 call (hap32 graph, 0cab3fbd4)
  took 172 s / 716 s / 7.5 GB. Without the mosaic, 8b993a339 on the original arm takes 173 s / 667 s
  / 11.3 GB.
- **The genome-wide ONT build's chr20 took 243 s.** That run used the 18-haplotype E821 graph. Its
  own build, 2a6a228a5, takes 484 s on this 34-haplotype graph with the same reads. So the factor of
  two comes from the graph, not from vg.
- **`--anchors-out` costs about 2.1x the ONT call.** It turns on descent into chains the reference
  does not cross: 94 k child calls against 58 k. The anchors run keeps only 2.2 of 5 threads busy,
  so part of it is serial. That has not been looked at.

### The regression the timings exposed, and its fix

The calls at 91d38c802 were ~25% slower than at 2a6a228a5, a week earlier. That is the build the
genome-wide ONT run used.
- **Short reads, no mosaic, original arm:** 164 s / 680 s CPU, against 206 s / 846 s at 283808454.
- **The cause was the depth-rate window on reference coordinates** (vg 6421bb2bc). It counted read
  starts over three 16 kb reference buckets per site by fetching their reads through the read source.
  - Buckets near a fetch-window edge loaded neighbouring 16,384-id windows into the 2-entry per-thread
    cache, evicting the window the sites were using.
  - Buckets straddling a window boundary ran an extra gbz-base subprocess each.
  - gbz-base subprocesses rose from 576 to 939, and reads fetched from 15.7 M to 22.8 M.
  - In a 60 s profile, threads waited longer on gbz-base for the buckets than for all the sites
    together.
- **The other 49 commits cost nothing measurable.**
- **The fix (vg 8b993a339, branch `depth-rate-tallies` in `~/CLionProjects/vg-speed`, not yet on the
  PR).**
  - The read source counts read starts itself, from tallies it makes the first time it fetches each
    window.
  - Its cache is shared by all threads and holds 4 windows per thread.
  - Every output (VCF, mosaic, anchors) is byte-identical, and TAP 18_vg_call passes 453/453.
- **Speed is back.**

  | call (original arm) | 2a6a228a5 | 283808454 | 8b993a339 |
  |---|---|---|---|
  | short reads | 164 s | 206 s | 173 s |
  | ONT | 484 s | 581 s | 484 s |
  | ONT + anchors | 1,282 s | 1,393 s | 1,030 s |

- **It costs ~1-2.7 GB more peak memory, from the larger cache.** A 2-window cache saves that memory
  but takes 189 s on short reads.
- **Fetching is still 70-80% of worker time, as it was at 2a6a228a5.** That is waiting on the
  gbz-base child process, then parsing its GAF text back into alignments. Prefetching the next
  window would hide the wait.

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
