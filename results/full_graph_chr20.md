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

- **Realigning every tandem repeat in the full graph improves SV calling for both read types, and
  clearly beyond noise.**
  - Short reads: refined SV F1 +0.040, raw +0.071. False positives fall by 40% (411 → 247 raw).
  - ONT: refined SV F1 +0.084, raw +0.098.
- **The gain is broad.** The ten most-improved regions hold 41-57% of it, and calls outside the regions
  do not change.
- **Small variants are unchanged, except a small indel cost on short reads** (0.9277 → 0.9262). Its
  regions have not been looked at yet.
- **The consistent losers** are TR_chr20_64970081 (+11 short, +9 ONT; 419 alleles, its sequence grew
  12.6 kb → 145 kb) and TR_chr20_18259720 (+8, +7).
- **Comparison with stage 4.** The hap32-projected medoid star (4x) gave +0.060 refined over its own
  re-mapped control. That was a different chain: the 34 production haplotypes, no re-sampling.

**Cost.** Calls take as long as on the original graph: short 296 s against 349 s, ONT 794 s against
857 s, at 8 threads.

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
