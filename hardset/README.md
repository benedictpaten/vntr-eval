# The hard tandem-repeat set

The chr20 tandem-repeat regions whose realigned alignment agrees worst with the pairwise optima. It is a
benchmark for realignment strategies.

## The set

`chr20.tsv` lists the 45 regions of the latest chr20 graph whose pair F1 is below 0.8, in
`work/full/chr20/qc/latest.tsv` (`pgrealign qc` with `--graph`). Pair F1 is 2 pair_tp / (pair_opt +
pair_msa), over qc's 20 random allele pairs per region. Across all 14,511 regions the median is 1.0.

- **Classes:** 27 VNTR, 7 STR+VNTR, 8 STR, and one each of LC, LC+VNTR and SAT.
- **Size:** median 220 distinct alleles (5-452), and the longest allele has a median of 4.3 kb
  (1.3-45.8 kb).
- **Phase flips:** 8 regions hold 9 of the 10 local ONT mis-phased stretches that remain on chr20
  (`flip_stretches`). The stretches come from vg-call-eval `work/offpanel/switches/switches_chr20_ont_1c8f89b7a.tsv`.

Columns of `chr20.tsv`:

| column | meaning |
|---|---|
| id, contig, ref_start, ref_end, span, L, R, targets, classes | as in `chr20.regions.tsv.gz` (CHM13 coordinates, boundary node ids) |
| period, motif | the catalogue period and motif of the region's longest target repeat |
| alleles, longest | distinct alleles, and the longest one's length |
| ref_allele | the allele CHM13 spells |
| pair_f1_graph, ref_f1_graph | qc's pair F1 and reference F1 of the original (Minigraph-Cactus) graph |
| pair_f1_realigned, ref_f1_realigned | the same for the current realignment (`msas.banded.jsonl.gz`) |
| realign_seconds | the current realignment's run time for the region, fragments included |
| flip_stretches | local ONT mis-phased stretches (pairs of switches) overlapping ref_start..ref_end |

Each region's alleles are in `work/hardset/chr20/fasta/<id>.fa`, which is git-ignored. These are all
the distinct panel alleles, as `pgrealign extract` wrote them, in the region's forward orientation. Each
header is `>a<k> weight=<w>`, where w is the number of panel haplotypes that spell the allele.

Regenerate both with `python3 tools/hardset.py` (`--max-f1` sets the threshold).

## Baseline

| | pair F1, median | pair F1, lower quartile | ref F1, median | run time |
|---|---|---|---|---|
| original graph | 0.647 | 0.509 | 0.688 | - |
| current realignment (medoid star) | 0.703 | 0.628 | 0.742 | 1,961 s in all, at most 601 s per region |

The realignment is better than the graph in 33 regions and worse in 8 (by more than 0.005).

## Strategies

`tools/hardset_strategies.py` runs each strategy on all 45 regions. It scores each one against the same
optimal alignments of the same pairs that qc uses. The cached optima reproduce qc's `pair_opt` and
`ref_opt` exactly, and the rerun star reproduces qc's pair F1 exactly. The strategies are described in
the script's docstring. Each (strategy, region) ran in its own process, three at a time, with limits of
9 GB resident and 1,800 s.

Results are in `results/hardset_chr20.tsv` (per region) and `results/hardset_chr20.summary.tsv`:

| strategy | aligned | pair F1 median | lower quartile | better / worse than star | cost ratio | time (s) | peak GB | nodes / star | sequence / star |
|---|---|---|---|---|---|---|---|---|---|
| original graph | 45 | 0.647 | 0.509 | | | | | | |
| **star** (current) | 45 | 0.703 | 0.628 | - | 1.34 | 1,773 | 6.6 | 1 | 1 |
| star_longest | 44 | 0.699 | 0.580 | 26 / 18 | 1.72 | 1,807 | 10.4 | 1.06 | 0.32 |
| **star_gapnorm** | 45 | **0.717** | **0.642** | 26 / 11 | 2.08 | star + 282 | 1.5 | 0.98 | 0.98 |
| **select** | 45 | **0.720** | 0.633 | **24 / 2** | 1.50 | 3,862 | 10.4 | 1.00 | 0.99 |
| abpoa_mc (progressive, MC BAR settings) | 42 | 0.621 | 0.487 | 3 / 36 | 1.42 | 2,382 | 11.0 | 0.91 | 0.32 |
| abpoa_prog (progressive, pgrealign scores) | 41 | 0.621 | 0.501 | 2 / 36 | 1.41 | 2,327 | 9.7 | 0.89 | 0.31 |
| famsa | 44 | 0.667 | 0.506 | 10 / 33 | 1.37 | 6,396 | 8.2 | 1.27 | 0.36 |
| mafft FFT-NS-2 | 43 | 0.613 | 0.486 | 10 / 33 | 1.50 | 8,502 | 8.3 | 1.33 | 0.76 |
| cluster_par (clusters as parallel blocks) | 45 | 0.598 | 0.478 | 4 / 29 | 2.21 | 1,525 | 8.1 | 1.02 | 1.02 |
| cluster_merge (clusters joined by mafft --merge) | 44 | 0.700 | 0.566 | 12 / 16 | 1.36 | 3,648 | 9.5 | 1.03 | 0.87 |
| units (motif-aware, tools/realign_units.py) | 44 | 0.630 | 0.467 | 6 / 35 | 1.43 | 3,207 | 5.5 | 1.23 | 0.34 |

How to read the table:
- **aligned** counts the regions a strategy finished. Wherever a strategy did not finish, the medians use
  the original graph's F1, since that graph is what would stay. Failures were concentrated in the four
  regions over 27 kb:
  - abPOA ran out of memory 6 times and timed out once;
  - FAMSA and mafft timed out 3 times;
  - star_longest and cluster_merge ran out of memory once each;
  - units failed once, when its mafft fallback timed out.
- **better / worse** is the region count against star, by more than 0.005 pair F1.
- **cost ratio** is the median, over regions, of the summed cost of the induced pairs under pgrealign's
  scores, divided by the cost of their optima. Unlike F1, it does not penalise a gap placed elsewhere at
  equal cost.
- **nodes / star** and **sequence / star** are medians of per-region ratios, measured on the
  alignment's column graph after unary chains are merged. Sequence is qc's measure: the distinct bases
  per column, summed.
- **select** chooses, per region, among star, star_longest and star_gapnorm. It picks the best pair F1 on
  20 other allele pairs, so the choice is not scored on the pairs that made it. It chose star 18 times,
  star_longest 15 times and star_gapnorm 12 times.

At the 8 phase-flip regions, select scores 0.22, 0.40, 0.61, 0.63, 0.65, 0.72, 0.71 and 0.79. Against
star, it moves them by at most +0.05. So no strategy repairs these regions.

## What it means

- **No aligner rescues the hard set.** Taking each region's best of all ten strategies gives a median of
  only 0.750.
- **The alternatives lose on pair F1.** The progressive aligners score 0.05-0.09 below star in median
  F1: abPOA with either setting, FAMSA, mafft and the unit-aware aligner. Joining clusters in parallel is
  the worst, at 0.598. The progressive aligners do build much smaller graphs: about a third of star's
  sequence. That is the star's known growth in many-allele regions.
- **star_gapnorm gains by matching the optima's gap convention.** Its F1 rises because qc's optima are
  left-normalised. But its induced pairs cost 2.08x their optima, where star's cost 1.34x, because
  shifting gaps against the column consensus splits them.
- **Recommendation for the graph test:**
  - Use **select**. It gains +0.018 median pair F1 over star and is worse in only 2 of 45 regions. It
    keeps star's graph size, and it costs about twice star's run time.
  - **star_gapnorm** is the cheap second choice: +0.014 at +16% run time. Use it only if the
    graph-level test shows that its extra gap cost does no harm.
  - If graph size matters more than pair F1, **star_longest** is the option. It keeps star's median F1
    to within 0.004 with a third of the sequence, but its lower quartile falls 0.05.

Each region scores only 20 pairs, so a difference of 0.01 in one region is noise. The medians and the
better/worse counts are the signal.
