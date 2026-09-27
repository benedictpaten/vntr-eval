# Results: does re-aligning the VNTR regions make the graph cleaner?

The question was whether simply re-aligning each region's haplotypes with a suitable multiple
aligner, and building the region's graph from that alignment, makes the hotspot VNTRs much cleaner
([plan.md](plan.md)). Background: in ~1% of VNTR regions, which hold half of short-read `vg call`'s
SV false positives, the Minigraph-Cactus (MC) graph implies pairwise alignments that cost 1.56x the
optimum, against 1.09x at matched VNTRs where genotyping works ([findings.md](findings.md)).

Every number below comes from a file in `results/`. Unless a line says otherwise, the metric is
cost/opt: the edit cost of the pairwise alignment the graph implies between two haplotypes, divided
by their optimal pairwise edit distance, taken as the median over the regions of a stratum. The
target was cost/opt <= 1.1 at hotspots, with no change at the controls. All metrics are defined in
[../tools/METRICS.md](../tools/METRICS.md).

## Summary

1. **On the graph itself, yes.** Re-aligning the 34 hap32 sequences takes hotspot VNTRs from 1.566
   (MC) to 1.048 with the repeat-unit-aware aligner, to 1.07-1.08 with the partial-order aligners,
   and to 1.11 with mafft L-INS-i. Hotspots at or below 1.1 go from 2/40 to 33/40. Homologous
   sequence left on parallel nodes almost vanishes: 40.5 bp/kb down to 0.36 (unit-aware), and at
   most 4.6 for any method.
   Source: `results/stage01_tables.md` (A).
2. **The controls mostly hold.** Failures: FFT-NS-2 regresses 10/40 matched controls by more than
   0.05. L-INS-i and G-INS-i, and the unit-aware method through its L-INS-i fallback, damage 3/20
   non-repeat SVs; the worst is L006556, 1.00 -> 1.19.
3. **The fair comparison is the full-panel arm.** MC aligned all ~450 HPRC haplotypes, and
   `mc.gfa` keeps only the 34 sampled rows. The like-for-like replacement is to align the full
   panel and project it onto hap32. Aligning hap32 alone is an easier problem.
4. **The full-panel arm keeps most of the gain at hotspots, but not at the controls.** Median
   cost/opt at hotspot VNTRs, MC / aligned on hap32 / full panel projected: unit-aware 1.566 /
   1.048 / 1.059 (98% of the gain kept, n=40); abPOA 1.582 / 1.083 / 1.130 (91%, n=30); spoa
   1.599 / 1.050 / 1.115 (88%, n=20). Matched controls keep only 23-66% of their small gain.
   Projection is worse than hap32 alone in 34/40 hotspots for unit-aware.
5. **MC is just as far from optimal on its own full panel** (1.513 at hotspot VNTRs). Realigned
   full-panel graphs reach 1.054 (unit-aware), better than full MC in 40/40 hotspots.
6. **Reads gain only modestly.** In a 19-region local re-map, only mafft L-INS-i (hap32 arm) places
   reads better at every hotspot where it ran (8/8). Core MAPQ<5 falls from 0.817 to 0.671;
   HG002's own two-haplotype graph gets 0.518, so tandem ambiguity is the floor. The full-panel
   arms place reads no better than MC. Read placement tracks k-mer redundancy (Spearman 0.75), not
   cost/opt (0.26).
7. **Representation shifts.** HG002's haplotypes stay spellable. But when HG002's truth is written
   through a realigned graph, raw truvari F1 at matched controls falls from 0.958 to 0.30-0.81,
   while whole-span phab stays at 1.000: the graphs cut the same sequences into different records
   than the truth VCF does.
8. **Where it fails:**
   - very long regions and outlier alleles (a 115 kb allele at L012272), which exhaust memory or
     time for spoa, abPOA and the all-pairs mafft modes;
   - differently ordered or compound arrays: 10/64 hotspots stay above 1.1 with every method on
     hap32, 21/64 on the full panel;
   - non-repeat SVs, for L-INS-i;
   - the all-pairs mafft modes, which could not run the full panel at hotspots at all
     (4/1/1 regions).
9. **For Minigraph-Cactus.** abPOA run with Cactus's own settings on each region alone already
   beats MC on the full panel (1.515 -> 1.269, better in 36/39 hotspot VNTRs). So most of the damage
   comes before the POA step (probably the pinches from minigraph mappings; not tested). A better
   aligner for these regions helps further: 1.111 with default abPOA, 1.054 unit-aware.
10. **Next:**
    - Stage 3 (local `vg call`), scored by haplotype edit distance and whole-span phab F1, not raw
      F1.
    - Rank candidates by k-mer redundancy as well as cost/opt.
    - Fix the unit-aware method's scattered unit gaps and its non-repeat fallback.
    - Get the one method that helped reads (L-INS-i) through the full panel.

## What was run

- **Regions:** 149 of the 150 stratified loci, packaged in `regions/`. L000601 has no anchor pair
  within 320 kb (see `regions/README.md`). They are 40 hotspot VNTRs, 24 other hotspots, 40
  matched VNTR controls, 25 VNTRs `vg call` gets right, and 20 non-repeat SV controls.
- **Methods:** eight realigners plus the MC baseline:
  - mafft FFT-NS-2, L-INS-i, E-INS-i and G-INS-i;
  - abPOA and spoa (global mode, longest sequence first);
  - abPOA with every Minigraph-Cactus BAR setting, including its 10 kb rule (`poa_abpoa_mc`);
  - the repeat-unit-aware aligner (`unit_aware`, [../tools/UNIT_ALIGN.md](../tools/UNIT_ALIGN.md));
  - on the full panel only, mafft FFT-NS-i as a cheaper iterative mode.

  Each graph is induced from the MSA by `tools/msa_graph.py`. Every path of every graph was checked
  to spell its sequence: 3,488 graphs scored, 0 invalid (`results/stage01.md`).
- **Two arms per method:**
  - `candidates/<m>/` aligns the 34 hap32 sequences.
  - `candidates/<m>__all/` aligns the distinct sequences of the full panel (median 452 haplotypes
    per region), keeps the hap32 rows, drops all-gap columns and builds the graph.

  The full-panel graphs themselves are in `work/panel/<m>/` (git-ignored) and were scored on
  full-panel pairs in `results/full/`.
- **Stages run:**
  - Stage 0 (graph only) and Stage 1 (placement redundancy of the existing reads) on every graph.
  - Stage 2 (local giraffe re-map) on 19 regions: 9 hotspot VNTRs, 1 other hotspot, 7 matched
    controls and 2 non-repeat controls (`results/stage2_remap.md`).
  - Stages 3-4 (local `vg call`, whole-contig splice) have not been run.

## 1. Aligning the 34 hap32 sequences

Median cost/opt per stratum, with the number of regions at or below 1.1 in brackets. n is given
where a method did not finish every region. The medians are over each method's own regions; paired
comparisons are in `results/stage01_tables.md`. Source: `results/<method>/*.json`, recomputed for
this page.

| method | hotspot VNTR (40) | hotspot other (24) | matched control (40) | correct control (25) | non-repeat SV (20) |
|---|---|---|---|---|---|
| MC (mc.gfa) | 1.566 (2) | 1.625 (3) | 1.102 (19) | 1.015 (17) | 1.000 (20) |
| mafft FFT-NS-2 | 1.327 (4) | 1.313 (6) | 1.100 (20) | 1.057 (15) | 1.000 (20) |
| mafft L-INS-i | 1.111 (18, n=39) | 1.164 (8) | 1.028 (33) | 1.015 (22) | 1.000 (19) |
| mafft E-INS-i | 1.182 (9, n=37) | 1.221 (6) | 1.030 (33) | 1.007 (22) | 1.000 (20) |
| mafft G-INS-i | 1.110 (15, n=39) | 1.163 (8) | 1.027 (34) | 1.034 (21) | 1.000 (19) |
| abPOA | 1.083 (25, n=37) | 1.150 (9) | 1.024 (33) | 1.006 (22) | 1.000 (20) |
| spoa | 1.067 (20, n=27) | 1.125 (10) | 1.016 (34, n=39) | 1.005 (24) | 1.000 (20) |
| abPOA, Cactus settings | 1.170 (11) | 1.284 (6) | 1.031 (31) | 1.007 (22) | 1.000 (20) |
| repeat-unit-aware | **1.048 (33)** | **1.053 (20)** | 1.015 (37) | 1.015 (25) | 1.000 (19) |

What else changes at hotspot VNTRs (medians; `results/stage01_tables.md` A):

- **Parallel homologous sequence** (bp/kb of sequence left on separate nodes when other haplotypes
  hold the same sequence) falls from 40.5 to 0.36 (unit-aware), 0.8-0.9 (POA) and 1.4-1.6 (L-/G-INS-i).
  This was the gross fault at L012184 and L012272. FFT-NS-2 leaves 4.6.
- **Gap-affine cost/opt** falls from 1.475 to 1.09-1.15 for the good methods. Here unit-aware
  (1.093), abPOA (1.100) and L-INS-i (1.110) are level.
- **SV-sized pieces per haplotype** fall only from 9.2 to 4.1 (L-INS-i) and 6.7 (unit-aware). In
  the median hotspot every non-CHM13 haplotype differs from CHM13 by >= 50 bp in length, so about 1
  piece per haplotype is the floor.
- **k-mer redundancy** (the same 21-mer on different nodes) falls only from 0.60 to 0.40-0.52.
  unit-aware is the worst of the good methods here (0.518), because it puts identical unit
  variants into different columns.
- **Nodes per kb** halves (764 to 265-381). It cannot be compared directly with MC, which keeps
  node cuts made by the other ~420 haplotypes: 48.7 nodes/kb even at non-repeat SVs, against 17-23
  for any MSA graph.

**Which method.** On the 135 regions where all eight methods finished, unit-aware is best or within
0.005 of the best on cost/opt at 24/27 hotspot VNTRs and 21/24 other hotspots
(`results/stage01_tables.md` C). It is not best on the other measures:

- L-INS-i and G-INS-i give the fewest SV pieces (3.7 against 5.4) and match it on gap-affine cost
  (1.115 against 1.117).
- E-INS-i has the lowest k-mer and read redundancy.
- unit-aware's objective charges whole-unit gaps linearly, so where a deletion goes is free and it
  scatters them.

**The controls.** Paired against MC, with regressions counted above 0.05 (`results/stage01.md` b):

- FFT-NS-2 regresses 10/40 matched and 7/25 correct controls; the worst is L000709, 1.00 -> 2.97.
- L-INS-i, G-INS-i and unit-aware (whose fallback for non-repeat regions is L-INS-i) regress
  3/20 non-repeat SVs: L006556 goes 1.00 -> 1.19 (affine 1.37), plus L002292 and L003883.
  E-INS-i and the POAs reproduce MC exactly there.
- abPOA and spoa regress one short-period correct control, L011224 (a 10 bp CA repeat): 1.027 ->
  1.137 and 1.115.
- Under gap-affine cost, unit-aware is worse than MC at 11/40 matched and 11/25 correct controls.

## 2. The panel question

**Why there is a question.** `mc.gfa` is not an alignment of the 34 hap32 sequences.
Minigraph-Cactus aligned the whole HPRC v2.1 panel. hap32 is 32 haplotypes sampled from that
graph for HG002 (`vg haplotypes`), plus CHM13 and GRCh38. So `mc.gfa` means: align ~450
haplotypes, then keep 34 rows.

This was checked directly, not assumed (`results/panel_premise.tsv`):
- All 345,486 nodes and 406,612 edges of the 149 `mc.gfa` files are in the full graph, with
  identical sequences.
- All 79,212 hap32 pairs that have a full-panel counterpart imply the same alignment and the same
  cost in both graphs.
- 4,620 of 4,734 sampled paths are verbatim full-panel walks. The other 114 are recombinant
  junctions (88) or reference copies (24) (`results/panel_stats.tsv`).

**Which comparison is fair.** Aligning the 34 hap32 sequences alone (m(hap32)) is an easier
problem than MC solved: 30 distinct sequences at the median hotspot, against 374 in the full panel
(`results/panel_stats.tsv`). It is also what a post-hoc fix of each sampled graph would do. The
like-for-like replacement for MC is to align the same full panel and project it onto hap32
(m__all).

Both panels come from the eval graph (`hprc-v2.1-mc-chm13-eval.gref.gbz`), which has no HG002,
HG003 or HG004. `gref_CHM13`, the reference-cover path stitched from other haplotypes, is
excluded. It was checked on all 149 regions: no gref fragment contains an anchor node, and the
reference copies were subtracted exactly (`results/panel_stats.tsv`).

**MC on its own panel.** The sampling does not explain MC's poor score. At hotspot VNTRs full MC
scores 1.513 on full-panel pairs against 1.566 for mc.gfa on hap32 (1.603 with a second sampling
seed). Excess edits per kb are 242 against 182, and k-mer redundancy 0.78 against 0.60. Matched
controls are worse on the full panel: 1.199 against 1.102. Source: `results/full/mc/`,
`results/stage01_tables.md` D.

**The paired comparison, projected graphs scored on hap32.** The columns are MC / m(hap32) / m__all.
"Kept" is the share of m(hap32)'s gain over MC that m__all keeps. "b / w" counts regions where m__all
is better / worse than m(hap32) by more than 0.005. Source: `results/<m>/`, `results/<m>__all/`,
recomputed for this page; it matches `results/stage01_tables.md` D.

| method | hotspot VNTR | kept | b / w | matched control | kept | b / w |
|---|---|---|---|---|---|---|
| unit-aware | 1.566 / 1.048 / 1.059 (n=40) | 98% | 2 / 34 | 1.102 / 1.015 / 1.045 (n=40) | 66% | 0 / 23 |
| abPOA | 1.582 / 1.083 / 1.130 (n=30) | 91% | 1 / 22 | 1.102 / 1.024 / 1.084 (n=40) | 23% | 1 / 22 |
| spoa | 1.599 / 1.050 / 1.115 (n=20) | 88% | 1 / 19 | 1.103 / 1.016 / 1.069 (n=37) | 40% | 0 / 25 |
| abPOA, Cactus settings | 1.567 / 1.170 / 1.275 (n=39) | 73% | 6 / 32 | 1.102 / 1.031 / 1.084 (n=40) | 26% | 2 / 19 |
| mafft FFT-NS-2 | 1.598 / 1.288 / 1.420 (n=31) | 57% | 5 / 26 | 1.102 / 1.100 / 1.206 (n=40) | - | 7 / 23 |
| mafft L-INS-i | 1.307 / 1.042 / 1.035 (**n=4**) | 102% | 2 / 1 | 1.112 / 1.022 / 1.054 (n=27) | 64% | 2 / 12 |

Against MC, unit-aware's projection is better at 40/40 hotspot VNTRs and 24/24 other hotspots.
abPOA's is better at 30/30 of the hotspot VNTRs it finished.

What changes when the full panel is aligned:

- **Gap placements are inherited.** The full-panel columns are placed to serve ~450 sequences.
  A hap32 row that carries one clean deletion in its own alignment inherits a deletion cut into
  pieces:
  - At L014297 (a 17 bp unit), spoa's indel runs per haplotype go from 28.5 (mean 37 bp) to 45.4
    (mean 24 bp), deleting one unit at a time (`tools/POA.md`).
  - At L005990 the unit-aware projection splits the shorter alleles' single deletion into many
    short gaps, with 2,648 nodes against 1,612.
- **Nodes per kb at hotspot VNTRs** rise from 381 to 454 (unit-aware, n=40) and 274 to 436 (abPOA,
  n=30).
- **The controls pay the most.** They start near the optimum, so any inherited fragmentation shows.
- **Truth spellability goes the other way.** The summed edit distance from HG002's haplotypes to
  the closest graph path at hotspot VNTRs is 815 for MC, 1,071 for unit-aware on hap32 and 670 for
  its projection. The full panel's extra cuts leave more recombinant routes (`results/stage01.md` e).
- **Coverage is the practical limit:**
  - The all-pairs mafft modes finished the full panel at only 4 (L-INS-i), 1 (E-INS-i) and 1
    (G-INS-i) hotspot VNTRs, so they cannot be judged there.
  - abPOA misses the 10 hotspot VNTRs with the longest alleles, and spoa misses 20.
  - Only unit-aware (149/149) and Cactus-settings abPOA (148/149) cover almost everything
    (`results/stage01.md`, coverage table).

**Full-panel graphs scored on full-panel pairs** (300 sampled pairs per region, the same pairs for
every graph; `results/full/`). At hotspot VNTRs:
- unit-aware goes from full MC's 1.513 to 1.054: better in 40/40, and 32/40 at or below 1.1.
- spoa reaches 1.094 (n=20), abPOA 1.111 (n=30), Cactus-settings abPOA 1.269 (n=39) and FFT-NS-2
  1.361 (n=31).
- A second sampling seed moved full MC's hotspot median to 1.603, but no comparison changed sign
  and no count moved by more than 4 regions (`results/stage01.md` d).

**So:** the fair comparison is m__all against mc. On it, the good realigners still solve most of
the hotspot problem on the graph metrics: unit-aware keeps 98% of its gain, abPOA 91%, spoa 88%.
The 34-sequence results overstate what an upstream change would achieve, by a little at hotspots
and by more at the controls.

## 3. Reads

**Stage 1: placement redundancy of the existing reads.** Pooled over MAPQ<5 core reads at hotspot
VNTRs, each method paired with MC on its own regions (`results/stage01.md` f):
- The "redundant only" class (the same read sequence on different nodes, where no haplotype
  carries it twice) is 29-31% for MC. It falls to 14.1% (E-INS-i), 17.7% (L-INS-i), 18.6%
  (Cactus-settings abPOA), 20.9% (abPOA) and 23.3% (unit-aware).
- The projected arms keep less of the fall: unit-aware__all 24.7%, Cactus-settings abPOA__all
  24.0%. abPOA__all and spoa__all remove none of it.
- Tandem ambiguity stays at 56% in every graph, as the investigation predicted.

**Stage 2: local re-map with giraffe** (19 regions; `results/stage2_remap.md`,
`results/stage2_remap.tsv`). Medians at the 9 hotspot VNTRs:

| measure | MC | L-INS-i (hap32 arm) | HG002's own two haplotypes |
|---|---|---|---|
| core reads at MAPQ<5 | 0.817 | 0.671 (better at 8/8) | 0.518 |
| truth-consistent core reads at MAPQ>=5 | 0.172 | 0.327 | 0.495 |
| core depth from MAPQ>=5 reads | 2.8x | 5.2x | 13.4x |

- **Local mapping reproduces the genome-wide one** at VNTRs: 0.795 genome-wide against 0.817
  locally, a paired difference of 0.003.
- **L-INS-i closes about a quarter of the gap** to HG002's own graph, without making confident
  placements more often wrong (0.131 against 0.132).
- **The other methods help less or not at all:**
  - unit-aware (the Stage 0 winner) is better at 6 and worse at 3 hotspots, and worse at 6/7
    matched controls.
  - The full-panel arms are no better than MC: unit-aware__all is better at 4 hotspots and worse
    at 5; spoa__all and abPOA__all are worse at every hotspot where they exist.
  - L-INS-i fails at the non-repeat control L006556 (MAPQ<5 0.007 -> 0.255), its known Stage 0
    regression.
- **What predicts placement is k-mer redundancy, not cost/opt.** Over 142 region-by-graph changes
  against MC, Spearman with the change in core MAPQ<5 is 0.75 for k-mer redundancy and 0.26 for
  cost/opt (-0.08 at hotspots).
- **One gross fault is fixed outright.** At L012272 (55 kb span, 42 kb of it flank) 91% of the
  span's reads are MAPQ<5 on MC's graph. On unit-aware's graph it is 8% (7% for unit-aware__all),
  and flank depth from confident reads rises from 4.0x to 30.4x.
- **The floor:** even HG002's own graph leaves 52% of hotspot core reads at MAPQ<5 (tandem
  ambiguity), and the reads missing from hotspot cores are missing from the read set, not from the
  graph.

## 4. Truth in the graph, and how it is written

- **Spellability:**
  - Every HG002 haplotype that is an exact panel copy stays exactly spellable in every graph (17
    at hotspot VNTRs).
  - MC spells 22 more exactly only by recombining through its mesh of 1 bp nodes. The realigned
    graphs keep 10-17 of those 22, and unit-aware__all keeps 19.
  - At the "alleles absent" pilot L002013, HG002's haplotypes are 594 and 578 edits from the
    closest hap32 haplotype, yet MC's graph spells both exactly and unit-aware's within 2 edits.
    "Absent from the panel" does not mean "absent from the graph".
  - Source: `results/stage01.md` e, `results/mc/L002013.json`.
- **Truvari at hotspot VNTRs.** Refined F1 of the truth written through the graph rises from 0.867
  (MC) to 0.92-0.98. Raw F1 rises only a little (0.26 -> 0.30-0.38 for mafft and POA) and not at
  all for unit-aware (0.260).
- **Truvari at controls: a red flag for Stage 3.** Raw F1 of the truth written through the graph
  at matched controls falls from 0.958 (MC) to 0.771 (L-INS-i), 0.648 (abPOA), 0.405 (spoa),
  0.323 (unit-aware) and 0.304 (unit-aware__all). Refined F1 falls to 0.74-0.96. Whole-span phab
  stays at 1.000, and the haplotypes are spelled exactly.
  - The cause: the truth VCF writes some repeat alleles as compensating SV pairs. At L001894 it has
    a 76 bp deletion plus a 72 bp insertion 300 bp apart. MC's graph writes the same records. The
    realigned graphs write the same haplotype with no SV-sized record, so truvari counts FNs.
  - unit-aware goes the other way, writing 55 SV records at matched controls where the truth has
    34.
  - Source: `results/stage01.md` b.

## 5. Where re-alignment fails

| failure | where | methods |
|---|---|---|
| memory or time on long regions and outlier alleles | L012272 (55 kb span; longest full-panel allele 115 kb), L015347, L001909, L007040, L011138; on the full panel L012184 (one 39.6 kb allele against CHM13's 4.8 kb) | spoa (14 hap32 / 26 full-panel regions not run), abPOA (3 / 11), L-/E-/G-INS-i (1-3 on hap32; 72-84 on the full panel, mostly skipped for time), FFT-NS-2 full panel (9) |
| arrays whose unit types are ordered differently, or with several motifs | L011463 (best 1.241), L009658 (1.175), L004145 (1.142), L000034 (1.127): 10/64 hotspots stay above 1.1 with every method on hap32, 21/64 projected | all |
| non-repeat SVs | L006556, L002292, L003883 | L-INS-i, G-INS-i, unit-aware's fallback |
| short-period controls | L011224 | abPOA, spoa |
| parallel copies recreated | L012184 (FFT-NS-2 1.86, worse than MC), L002143 (FFT-NS-2 9.5) | FFT-NS-2 |
| BAR's 10 kb rule leaves the middle unaligned | L012272 (hap32), L012184 full panel (140 kb unaligned) | Cactus-settings abPOA, by design |
| representation shift in truvari | every VNTR control stratum | every method (worst: unit-aware, spoa) |

Not tested: fragments (paths ending inside the array) are in no graph, and every metric uses
spanning paths. N-gap haplotypes work: the one in hap32 (3,146 bp of N at L012184) became a single
N node in every graph. Sources: `results/stage01.md` c, `results/realign_runtime*.tsv`,
`results/mcpoa_runtime*.tsv`.

## 6. What this means for Stages 3-4

- **Score Stage 3 on haplotypes, not records.** Raw truvari F1 at the controls will fall for
  representation alone (section 4). Judge local `vg call` by the called haplotypes' edit distance
  to HG002, whole-span phab F1 and refined F1. Treat raw F1 as a representation score.
- **Candidates for Stage 3:**
  - **mc**, the baseline.
  - **unit_aware__all**, the fair replacement that covers every region with near-full Stage 0 gain.
  - **mafft_linsi** (hap32 arm), the only graph that improved read placement at every hotspot
    (though not at the non-repeat control L006556).
  - **unit_aware** (hap32 arm).
  - **HG002's own graph**, as the ceiling.

  Run it on the 19 Stage 2 regions first: `tools/remap_local.py` already rebuilds their reads,
  estimates fragment lengths and indexes each local graph.
- **Expect modest gains from alignment alone.** MAPQ improves by about a quarter of what is
  possible, the floor is tandem ambiguity, and SV pieces per haplotype stay at 4-7 where about 1
  is ideal. The caller-side change in [findings.md](findings.md) (strategy 2: genotype each problem
  VNTR as one site of whole alleles) is what addresses the rest, and it composes with a cleaner
  graph.
- **Stage 4 should splice full-panel graphs, not hap32 graphs.** hap32 is re-sampled from the full
  graph for each sample, so the production path is:
  1. splice the realigned full-panel region graphs (`work/panel/<m>/`) into the eval graph
     between their anchor nodes;
  2. re-run haplotype sampling, giraffe and `vg call` on one contig.

  Splicing hap32-arm graphs into hap32 only tests the easier problem. The anchor nodes are
  shared by all graphs, so the splice points exist; new node ids mean re-indexing.

## 7. What this means for doing it upstream in Minigraph-Cactus

- **Most of the damage seems to come before the POA step.** abPOA run with every Cactus BAR
  setting on each whole region already beats MC: full panel 1.515 -> 1.269 at hotspot VNTRs,
  better in 36/39. In the real pipeline, CAF first pinches blocks from the minigraph mappings
  inside the region, and BAR only fills the gaps between them (`tools/POA.md`). The likely source
  of MC's misaligned arrays is therefore those pinches, not BAR's scoring. This is inferred, not
  tested.
- **The upstream change this suggests:**
  1. Do not pinch inside catalogue VNTRs (for example the >= 1 kb set, 3,111 benchmarked regions,
     6.1 Mb), so that each VNTR becomes one bubble.
  2. Align that bubble's alleles with a better aligner than BAR's settings: unit-aware (1.054 on
     the full panel) or abPOA with default scoring (1.111), both better than Cactus's settings
     (1.269). Handle alleles over BAR's 10 kb banding limit, which left 140 kb unaligned at
     L012184.
- **Cost is affordable.** unit-aware aligned the full panel of all 149 regions in 12,497 s of
  aligner time: median 13 s, maximum 1,722 s, peak 9.7 GB (`results/realign_runtime.all.units.tsv`).
  POA memory is set by single outlier alleles: 11 regions for abPOA and 26 for spoa were refused on
  predicted memory above 12 GB.
- **Minimising edit cost is not enough for reads.** The full-panel arms did not improve read
  placement in Stage 2, and the measure that tracks placement is k-mer redundancy. An upstream
  aligner should also put identical sequence on shared nodes; L-INS-i and E-INS-i do this best on
  hap32. It is not known whether they keep that advantage on the full panel, because they could not
  run it at hotspots.

## 8. What to do next

1. **Stage 3 on the Stage 2 pilot**, with the candidates and scoring in section 6. This decides
   whether the Stage 0 gains reach genotypes.
2. **Rank candidates on k-mer redundancy alongside cost/opt**, and add it to the target: hotspots
   should approach the controls' 0.05, not only cost/opt 1.1.
3. **Improve the unit-aware method** where it loses:
   - put identical unit variants in shared columns (k-mer redundancy 0.52 against L-INS-i's
     0.42);
   - use an affine whole-unit gap only where it helps (tried globally: it helps VNTRs and hurts
     STRs);
   - replace the L-INS-i fallback for non-repeat regions with E-INS-i or abPOA, which reproduce MC
     exactly there;
   - model several motifs, for the 10 hotspots nothing fixes.
4. **Get a read-friendly aligner through the full panel:** L-INS-i or E-INS-i on the distinct
   full-panel sequences of the hotspots. This needs more memory or time than this machine gave
   them: on hap32 alone, scaling from the other runs put L012272 at about 30 GB or more.
   Alternatively, seed them from the unit-aware alignment.
5. **Reduce the projection loss** (untried): weight sequences by haplotype count in the full-panel
   alignment, or re-polish the hap32 rows after projecting.
6. **Test the Cactus hypothesis directly:** run `cactus-pangenome` on one contig with CAF pinching
   disabled inside the >= 1 kb VNTRs, and compare the resulting subgraphs with `evaluate.py --panel`.
7. **Scale** to `regions/vntr_ge1kb.bed` (3,111 regions) with the method that survives Stage 3.

## Runtime and memory

Per method and panel, from the runtime tables every driver wrote (`results/runtime_summary.md`,
`tools/runtime_summary.py`). Time is the aligner's wall clock on a shared, often heavily loaded
10-core Mac, with 1-3 threads depending on the method, so it is not CPU time. Peak memory is exact
for Cactus-settings abPOA (`/usr/bin/time -l`) and sampled about once a second for the others, which
makes those lower bounds; unit-aware on hap32 was not memory-sampled at all. "Refused" regions were
not run because a predicted time or memory was over the cap (hap32: 900 s per mafft call, 12 GB;
full panel: 1,800 s per region, 10-12 GB), so they are the largest regions and every statistic
below is censored from above.

| method | hap32: finished | median (hotspot) / max time | max peak | full panel: finished | median (hotspot) / max time | max peak | largest refused, predicted |
|---|---|---|---|---|---|---|---|
| mafft FFT-NS-2 | 149 | 0.6 s (5 s) / 188 s | 8.6 GB | 140 | 5 s (302 s) / 3,054 s | 8.4 GB | 8 GB, 15 h |
| mafft FFT-NS-i | - | - | - | 105 | 6 s (204 s) / 1,563 s | 6.5 GB | 8 GB, 30 h |
| mafft L-INS-i | 148 | 5 s (47 s) / 737 s | 12.2 GB | 77 | 31 s (463 s) / 1,029 s | 5.6 GB | 244 GB, 515 h |
| mafft E-INS-i | 146 | 3 s (49 s) / 397 s | 11.9 GB | 65 | 19 s (190 s) / 812 s | 7.9 GB | 244 GB, 515 h |
| mafft G-INS-i | 148 | 4 s (46 s) / 847 s | 14.2 GB | 65 | 23 s (128 s) / 321 s | 5.6 GB | 244 GB, 515 h |
| abPOA | 146 | 0.1 s (1 s) / 10 s | 7.3 GB | 138 | 1 s (25 s) / 269 s | 8.7 GB | 268 GB |
| spoa | 135 | 0.9 s (5 s) / 22 s | 9.5 GB | 123 | 9 s (91 s) / 269 s | 10.2 GB | 993 GB |
| abPOA, Cactus settings | 149 | 1.1 s (4 s) / 18 s | 5.4 GB (exact) | 148 | 21 s (124 s) / 416 s | 6.9 GB (exact) | MSA too large to build (L012272) |
| repeat-unit-aware | 149 | 2.2 s (6 s) / 421 s | not sampled | 149 | 13 s (132 s) / 1,722 s | 9.4 GB | - |

- **On hap32 everything is cheap.** Every method's median is under 5 s. Only the all-pairs mafft
  modes need minutes (up to 14 min) and more than 10 GB, at the few largest regions; L012272
  killed all three at 10-11 GB, and E-INS-i also timed out at L001909 and L015347.
- **On the full panel, coverage is set by memory and all-pairs time, not typical runtime.**
  - Only unit-aware finishes every region: 3.5 h of aligner time in total, at most 29 min and
    9.4 GB for one region.
  - Cactus-settings abPOA finishes 148 because BAR's 10 kb rule caps every string it aligns.
  - Default abPOA and spoa are fast where they run, but their memory grows with the square of the
    longest allele, so single outlier alleles put 11 and 26 regions out of reach; spoa's prediction
    for L012272 is ~1 TB.
  - The all-pairs mafft modes are predicted to need up to 515 h and 244 GB at the largest region.
- **To compare methods fairly, use the shared-region tables** in `results/runtime_summary.md`
  (135 regions on hap32, 104 on the full panel).
  - On hap32: abPOA 0.1 s, FFT-NS-2 0.5 s, spoa 0.9 s, Cactus-settings abPOA 1.1 s, unit-aware 2.0 s
    and L-/E-/G-INS-i 3-4 s median per region.
  - On the full panel: abPOA 0.5 s, FFT-NS-2 2.5 s, FFT-NS-i 5.9 s, spoa 6.2 s, unit-aware 6.6 s
    and Cactus-settings abPOA 14.4 s.
- Graph building, projection and evaluation are not included.

## Before/after pages

Six self-contained viewer pages are in [../results/pages/](../results/pages/). Each shows MC, the
best hap32-arm candidate and the best full-panel-arm candidate. "Best" means the lowest cost/opt
in `results/<method>/<id>.json`. Each page has an MSA panel that switches between the three
alignments, Bandage images, node tracks, calls and reads. Every page is 0.9-2.9 MB, loads no
external resource, and rendered with all three graph images and no console errors.

| page | stratum | MC | best hap32 arm | best full-panel arm (projected) |
|---|---|---|---|---|
| [L012184](../results/pages/L012184.html) | hotspot, parallel duplicated path, N-gap haplotype | 1.599 | spoa 1.081 | unit-aware 1.183 |
| [L005990](../results/pages/L005990.html) | hotspot, block of units near the array start | 1.211 | unit-aware 1.014 | unit-aware 1.042 |
| [L014297](../results/pages/L014297.html) | hotspot, 17 bp unit, 153 copies | 1.303 | unit-aware 1.022 | unit-aware 1.057 |
| [L002013](../results/pages/L002013.html) | hotspot, HG002 alleles absent from the panel | 1.406 | unit-aware 1.056 | unit-aware 1.080 |
| [L007172](../results/pages/L007172.html) | matched control (the projected POAs damage its reads) | 1.103 | unit-aware 1.026 | unit-aware 1.045 |
| [L001894](../results/pages/L001894.html) | matched control (the truvari representation shift) | 1.273 | unit-aware 1.003 | unit-aware 1.060 |

Rebuild one with:

```
python3 tools/viewer.py L014297 \
  --candidate candidates/unit_aware/L014297.gfa --candidate-label "unit_aware (hap32 arm)" \
  --candidate candidates/unit_aware__all/L014297.gfa --candidate-label "unit_aware (full panel, projected)" \
  --out results/pages/L014297.html
```

## Sources

- **Stage 0-1:**
  - `results/stage01.md` (narrative), `results/stage01_tables.md` and `results/stage01_pairs.tsv`
    (`tools/stage01.py`);
  - `results/summary.{md,tsv}` and `results/per_region.tsv` (`tools/summarise.py`);
  - per graph, `results/<method>/<id>.json`, and for the full-panel graphs
    `results/full/<method>/<id>.json` (`tools/evaluate.py`).
- **Stage 2:** `results/stage2_remap.md` and `results/stage2_remap.tsv` (`tools/remap_local.py`).
- **Panel:** `results/panel_stats.tsv`, `results/panel_premise.tsv` (`tools/panel.py`).
- **Runtimes and failures:**
  - `results/realign_runtime.tsv` (hap32 arm);
  - `results/realign_runtime.all.{mafft,poa,units}.tsv` (full panel);
  - `results/mcpoa_runtime{,.all}.tsv` (Cactus-settings abPOA).
- **Methods:** [../tools/README.md](../tools/README.md), [../tools/UNIT_ALIGN.md](../tools/UNIT_ALIGN.md),
  [../tools/POA.md](../tools/POA.md), [../tools/METRICS.md](../tools/METRICS.md).

The key numbers here were recomputed from the per-region JSONs and TSVs. They agree with the
files listed above, including these cases:
- the three-way comparison mc / m / m__all for unit-aware, abPOA, spoa, Cactus-settings abPOA,
  FFT-NS-2 and L-INS-i;
- full MC against each full-panel graph;
- truth sums and pooled truvari F1;
- the pooled read classes;
- the Stage 2 medians.

Stage 2 counts a region better or worse on any difference. The Stage 0-1 tables use a 0.005
threshold. Other summaries of these runs that quote slightly different counts counted any
difference: for example, the projected unit-aware graph is worse than the hap32 one in 38/40
hotspot VNTRs at threshold 0, and in 34/40 above 0.005.
