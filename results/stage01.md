# Stage 0-1: realigned VNTR graphs against Minigraph-Cactus

Every candidate graph in `candidates/` was scored with `tools/evaluate.py`, both arms: `<method>`
(aligned on the 34 hap32 sequences) and `<method>__all` (the full HPRC panel aligned, then projected
onto hap32). The baseline is `mc` (`regions/<id>/mc.gfa`). All of them were also scored with the
existing reads (`--reads`). The full-panel graphs were scored on the full panel (`--panel`).

- **What was scored:**
  - 2,329 hap32-scored graphs: mc on 149 regions plus 17 candidate sets, in `results/<method>/`;
  - 1,159 full-panel graphs: full MC plus 9 methods, in `results/full/<method>/`.
- **Validity:** all 3,488 are valid. Every path spells its sequence, and no evaluation failed.
- **Tables:** `results/summary.{md,tsv}` and `results/full/summary.{md,tsv}` (tools/summarise.py), and
  [stage01_tables.md](stage01_tables.md) with its long form `stage01_pairs.tsv` (tools/stage01.py).
  Every number below comes from these.
- **Pairing:** every comparison is paired over the regions where all the graphs it names exist, and
  n is given each time.
- **Metrics:** defined in [../tools/METRICS.md](../tools/METRICS.md).
  - cost/opt is the graph-implied pairwise edit cost over the optimum.
  - "Parallel homologous bp/kb" is `U_per_kb`.
  - SV pieces/path is the per-region mean.

## Short answers

**(a) Hotspots against controls.** On cost/opt, yes: every method except FFT-NS-2 moves hotspot
VNTRs most of the way to the matched controls.
- **Cost/opt on hap32:** the median goes from 1.566 to 1.048 (repeat-unit-aware), 1.067 (spoa, n=27),
  1.083 (abPOA, n=37) and 1.11 (L-INS-i / G-INS-i). That is 98-112% of the distance to mc's matched
  controls (1.102).
- **Target:** hotspots at cost/opt <= 1.1 go from 2/40 to 33/40 with unit_aware.
- **Parallel homologous sequence:** almost all of it goes (40.5 to 0.4-4.6 bp/kb).
- **What does not close:**
  - SV pieces per haplotype: 9.2 to 4.1-8.4 (10-60% closed). The ideal here is about 1 per haplotype.
  - k-mer redundancy: 0.60 to 0.41-0.55 (10-36% closed).
  - Mean node length: 3.6 to 4.3-9.1 bp, against 11.6 at the controls.

  These metrics also depend on how diverse the hotspot haplotypes really are.
- **Nodes/kb:** it halves, but it cannot be compared with mc. MC's subgraph keeps node cuts from the
  full panel: 48.7 nodes/kb even at non-repeat SV controls, against 17-23 for every MSA graph.

**(b) Controls.** On cost/opt the controls hold or improve, with two exceptions.
- **FFT-NS-2:** it regresses a quarter of the VNTR controls by more than 0.05 (10/40 matched, 7/25
  correct), for example L000709 1.00 to 2.97.
- **L-INS-i and G-INS-i:** they damage 3 of 20 non-repeat SV controls, for example L006556 1.00 to
  1.19. unit_aware inherits this, because L-INS-i is its fallback.
- **Truth written by the graph (red flag):** raw truvari F1 at the VNTR controls falls for every
  method.
  - Matched controls: 0.958 to 0.32-0.81.
  - Correct controls: 0.985 to 0.53-0.91.
  - Refined F1 recovers only partly (0.75-0.96 at matched controls). Whole-span phab stays at 1.00,
    so the sequences are right; the change is in how they are cut into records.
  - Cause: MC's graph writes HG002's alleles the way the truth VCF does, including compensating
    deletion/insertion pairs inside repeats. The realigned graphs write them more parsimoniously as
    small edits, or as scattered sub-SV pieces.

**(c) Which method.**
- **By cost/opt (unit edit cost):** repeat-unit-aware (unit_aware) is best. It is best or within
  0.005 of the best at 24/27 hotspot VNTRs and 21/24 other hotspots, on the 135 regions where all
  eight methods exist.
- **By the other metrics it is not best.**
  - mafft L-INS-i/G-INS-i match or beat unit_aware on gap-affine cost, SV pieces, k-mer redundancy and
    read redundancy.
  - E-INS-i and abPOA-with-MC-settings beat it on SV pieces and redundancy, but not on affine cost.
  - unit_aware's linear unit-gap objective scatters whole-unit gaps.
- **POA:** spoa and abPOA have the lowest unit cost after unit_aware on hap32, but fall behind L-INS-i
  under affine scoring. Their projected arms regress the most controls and remove none of the read
  redundancy.
- **FFT-NS-2 fails outright.** It leaves unaligned copies and is catastrophic at L002143 (9.5).
- **Where methods fail:**
  - very long regions and outlier alleles (spoa 14, abPOA 3, L-/E-/G-INS-i 1-3 missing on hap32;
    far more on the full panel);
  - mixed-motif and differently ordered arrays: at 10 of 64 hotspots no method reaches 1.1 (L011463,
    L009658, L004145, L000034);
  - non-repeat SVs for L-INS-i;
  - BAR's 10 kb rule for abPOA-MC (L012272).
- **N-gap haplotypes (one, L012184)** are handled by every method: the 3,146 bp N run becomes one N
  node.
- **Fragments** are in no graph, so these metrics cannot judge them.

**(d) The panel question.** Aligning the full panel and projecting keeps most, but not all, of the
gain, and costs more at the controls.
- **Hotspot VNTRs, cost/opt, mc / m(hap32) / m__all:**
  - unit_aware: 1.566 / 1.048 / 1.059 (n=40; 98% of the gain kept; m__all worse than m(hap32) in
    34/40);
  - abPOA: 1.582 / 1.083 / 1.130 (n=30; 91%);
  - spoa: 1.599 / 1.050 / 1.115 (n=20; 88%);
  - abPOA-MC: 1.567 / 1.170 / 1.275 (n=39; 73%);
  - FFT-NS-2: 1.598 / 1.288 / 1.420 (n=31; 57%).
- **Controls:** the projection costs the most here. Matched controls go from 1.02-1.03 on hap32 to
  1.05-1.08 when projected (worse in 19-25 of 37-40 for POA and unit_aware).
- **The all-pairs mafft modes** finished the full panel on only 4 (L-INS-i) and 1 (E-/G-INS-i)
  hotspot VNTRs, so they cannot be judged on it. On those few regions they hold their gain.
- **The full MC graph on the full panel** is as far from optimal as mc.gfa is on hap32:
  - hotspot VNTRs 1.513 against 1.566 (paired: 18 better, 21 worse; 1.603 with a second sampling
    seed);
  - matched controls 1.199 against 1.102;
  - excess edits/kb 242 against 182; parallel homologous sequence 55 against 41 bp/kb; k-mer extra
    0.78 against 0.60.
- **Realigned full-panel graphs** get close to optimal on the full panel. At hotspot VNTRs:
  - unit_aware 1.054, better than full MC in 40/40, with 32/40 at <= 1.1;
  - spoa 1.094 (n=20), abPOA 1.111 (n=30), abPOA-MC 1.269 (n=39), FFT-NS-2 1.361 (n=31).

**(e) Truth in the graph.**
- **Exact spelling:** every HG002 haplotype that is an exact panel copy stays exactly spellable in
  every graph (17/17 at hotspot VNTRs).
  - mc also spells 22 more exactly by recombining pieces of its 1-bp mesh. The realigned graphs
    keep 10-17 of those 22 (unit_aware 17, L-INS-i 11), and unit_aware__all keeps 19.
  - Summed edit distance at hotspot VNTRs rises from 815 (mc) to 826-1,845 in the hap32 arm, but
    falls to 670 for unit_aware__all. All are 8-23x below the closest single panel path
    (15,596).
- **Truvari F1 at hotspot VNTRs:**
  - Refined F1 rises for every method: 0.867 to 0.92-0.98.
  - Raw F1 rises modestly for the mafft and POA modes (0.26 to 0.30-0.38) and not at all for
    unit_aware (0.260).
  - phab stays at 0.99-1.00.

**(f) Stage 1: reads.** The redundant class of MAPQ<5 core reads (same sequence on different nodes,
no haplotype holding it twice) falls in the hap32 arm.
- **Pooled at hotspot VNTRs:**
  - 29.7% to 14.1% (E-INS-i, n=37), 17.7% (L-/G-INS-i), 18.6% (abPOA-MC), 20.9% (abPOA),
    22.7% (FFT-NS-2), 23.3% (unit_aware);
  - spoa: 35.4% to 28.7% (n=27).
- **Unique placements** rise from 13% to 20-30%.
- **Tandem ambiguity** stays at 56% by construction.
- **Projected arms** keep less of the fall: unit_aware__all 24.7%, abPOA-MC__all 24.0%. abPOA__all
  (33.4% from 34.3%) and spoa__all (42.3% from 41.5%) remove none.

## Coverage: which graphs exist

| method | family | hap32 arm | `__all` arm (hotspot VNTRs) | full-panel graphs | why regions are missing |
|---|---|---|---|---|---|
| mc | Minigraph-Cactus | 149 | - | 149 | - |
| `mafft_fftns2` | progressive MSA | 149 | 140 (31) | 140 | all: 8 not attempted after timeouts, 1 timeout (all hotspot VNTRs) |
| `mafft_fftnsi` | progressive MSA | - | 105 (9) | 105 | all: 43 skipped (predicted too slow or batch budget), 1 timeout |
| `mafft_linsi` | progressive MSA | 148 | 77 (4) | 77 | hap32: L012272 memout; all: 70 skipped, 1 timeout, 1 memout |
| `mafft_einsi` | progressive MSA | 146 | 65 (1) | 65 | hap32: 2 timeouts, 1 memout; all: 83 skipped, 1 memout |
| `mafft_ginsi` | progressive MSA | 148 | 65 (1) | 65 | hap32: L012272 memout; all: 83 skipped, 1 memout |
| `poa_abpoa` | POA | 146 | 138 (30) | 138 | predicted memory > 12 GB (hotspot VNTRs with long alleles) |
| `poa_spoa` | POA | 135 | 123 (20) | 123 | predicted memory > 12 GB (13 + 1 on hap32; 26 on the full panel) |
| `poa_abpoa_mc` | POA, Cactus settings | 149 | 148 (39) | 148 | all: L012272 skipped (every haplotype > 20 kb) |
| `unit_aware` | repeat-unit-aware | 149 | 149 (40) | 149 | - |

Missing regions are not random. They are the longest hotspot VNTRs: L012272, L015347, L001909,
L007040, L007009, L011138 and others. A method with few `__all` graphs is therefore compared on
small and medium regions only, and every table says so through its n.

## (a) Do hotspot VNTRs come to look like the matched controls?

Hap32 arm, hotspot VNTRs. Each cell is the method's median. In brackets is the share of the distance
from mc's hotspot median to mc's matched-control median that it closes. Both medians are taken on the
same hotspots, so methods with fewer regions are compared fairly with mc.

| method | n | cost/opt | excess edits/kb | affine cost/opt | SV pieces/path | k-mer extra | parallel hom. bp/kb | nodes/kb | mean node bp | hotspots <= 1.1 | own matched controls, cost/opt |
|---|---|---|---|---|---|---|---|---|---|---|---|
| mc | 40 | 1.566 | 182 | 1.475 | 9.21 | 0.601 | 40.5 | 764 | 3.6 | 2/40 | 1.102 |
| `mafft_fftns2` | 40 | 1.327 (52%) | 97.0 (48%) | 1.290 (42%) | 8.36 (10%) | 0.547 (10%) | 4.64 (89%) | 372 | 7.0 | 4/40 | 1.100 |
| `mafft_linsi` | 39 | 1.111 (98%) | 36.3 (83%) | 1.110 (83%) | 4.13 (60%) | 0.419 (33%) | 1.58 (96%) | 355 | 6.0 | 18/39 | 1.028 |
| `mafft_einsi` | 37 | 1.182 (83%) | 58.8 (68%) | 1.149 (74%) | 7.16 (23%) | 0.405 (36%) | 2.94 (92%) | 265 | 9.1 | 9/37 | 1.030 |
| `mafft_ginsi` | 39 | 1.110 (98%) | 39.6 (81%) | 1.112 (83%) | 4.64 (54%) | 0.418 (33%) | 1.41 (97%) | 352 | 6.1 | 15/39 | 1.027 |
| `poa_abpoa` | 37 | 1.083 (104%) | 26.0 (88%) | 1.100 (87%) | 5.44 (45%) | 0.479 (22%) | 0.78 (98%) | 300 | 5.5 | 25/37 | 1.024 |
| `poa_spoa` | 27 | 1.067 (107%) | 16.8 (93%) | 1.135 (82%) | 5.15 (43%) | 0.515 (12%) | 0.87 (97%) | 345 | 4.5 | 20/27 | 1.016 |
| `poa_abpoa_mc` | 40 | 1.170 (85%) | 49.0 (75%) | 1.148 (75%) | 4.70 (53%) | 0.423 (32%) | 4.56 (89%) | 299 | 6.5 | 11/40 | 1.031 |
| `unit_aware` | 40 | 1.048 (112%) | 15.2 (94%) | 1.093 (87%) | 6.74 (29%) | 0.518 (15%) | 0.36 (99%) | 381 | 4.3 | 33/40 | 1.015 |
| *mc, matched controls* | 40 | 1.102 | 4.4 | 1.037 | 0.67 | 0.053 | 0.15 | 101 | 11.6 | 19/40 | |

- **Hotspot_other (24)** behaves the same way. Cost/opt goes from 1.625 to 1.053 (unit_aware, 20/24 at
  <= 1.1), 1.125 (spoa), 1.150 (abPOA) and 1.164 (L-INS-i).
- **Measured against each method's own controls,** the hotspot/control gap is much smaller than under mc
  on cost/opt: 1.048 against 1.015 for unit_aware, 1.111 against 1.028 for L-INS-i. The residual gap is
  still large on excess edits/kb (15 against 1.1), SV pieces (6.7 against 0.47) and k-mer extra (0.52
  against 0.07).
- **SV pieces have a floor near 1 per haplotype.** In the median hotspot VNTR every non-CHM13 haplotype
  differs from CHM13 by >= 50 bp in length (median fraction 1.00), against 0.53 at matched controls. So
  an ideal graph writes about 1 SV piece per haplotype at hotspots; the best method still writes 4.1.
- **The full-panel arm projected (m__all)** closes less of the gap:
  - unit_aware__all 1.059 (109%, n=40, 31/40 at <= 1.1);
  - spoa__all 1.115 (97%, n=20); abPOA__all 1.130 (94%, n=30);
  - abPOA-MC__all 1.275 (63%, n=39); FFT-NS-2__all 1.420 (36%, n=31).

## (b) Are the controls left unchanged?

Paired against mc. b / w = regions better / worse than mc by more than 0.005 cost/opt; "> 0.05" counts
regressions larger than 0.05.

| method | matched control (n=40): median, b / w, > 0.05 | correct control (n=25) | non-TR SV control (n=20) | largest control regressions |
|---|---|---|---|---|
| mc | 1.102 | 1.015 | 1.000 | |
| `mafft_fftns2` | 1.100, 14 / 15, **10** | 1.057, 9 / 8, **7** | 1.000, 0 / 1, 0 | [L000709](../regions/L000709/) +1.97, [L008887](../regions/L008887/) +1.62, [L011121](../regions/L011121/) +1.29 |
| `mafft_linsi` | 1.028, 22 / 3, 0 | 1.015, 11 / 4, 2 | 1.000, 0 / 4, **3** | [L006556](../regions/L006556/) +0.19, [L002292](../regions/L002292/) +0.09, [L003883](../regions/L003883/) +0.06, [L007033](../regions/L007033/) +0.06 |
| `mafft_einsi` | 1.030, 22 / 1, 0 | 1.007, 10 / 1, 0 | 1.000, 0 / 1, 0 | none above 0.02 |
| `mafft_ginsi` | 1.027, 23 / 2, 0 | 1.034, 11 / 6, 3 | 1.000, 0 / 4, **3** | [L006556](../regions/L006556/) +0.19, [L003119](../regions/L003119/) +0.12 |
| `poa_abpoa` | 1.024, 24 / 0, 0 | 1.006, 12 / 4, 1 | 1.000, 0 / 1, 0 | [L011224](../regions/L011224/) +0.11 |
| `poa_spoa` | 1.016, 23 / 3, 0 (n=39) | 1.005, 13 / 1, 1 | 1.000, 0 / 0, 0 | [L011224](../regions/L011224/) +0.09 |
| `poa_abpoa_mc` | 1.031, 20 / 2, 0 | 1.007, 9 / 2, 0 | 1.000, 0 / 1, 0 | none above 0.03 |
| `unit_aware` | 1.015, 24 / 3, 0 | 1.015, 13 / 4, 1 | 1.000, 0 / 4, **3** | L006556, L002292, L003883 (its L-INS-i fallback); [L015743](../regions/L015743/) +0.05 |

**Projected arm (m__all):** the controls fare worse.
- abPOA__all regresses 8/40 matched and 7/25 correct controls by more than 0.05. spoa__all regresses
  5 and 8.
- FFT-NS-2__all regresses 16 and 6. At [TR687489](../regions/TR687489/) cost/opt reaches 31.8: its
  panel is nearly identical (optimum 2.2 edits/kb), so this is 68 excess edits/kb.
- unit_aware__all regresses 3 and 2. The worst is [TR101140](../regions/TR101140/), +0.43, where its
  fallback becomes FFT-NS-i above 100 distinct sequences.

**Gap-affine scoring** exposes a cost that unit cost hides. unit_aware is worse than mc under affine
scoring at 11/40 matched and 11/25 correct controls; abPOA at 4 and 8; spoa at 9 and 6. L-INS-i and
E-INS-i are worse at 0-5. unit_aware's median affine cost/opt at correct controls is 1.024, against
1.012 for mc and 1.002-1.005 for the mafft INS-i modes.

**Red flag: how the truth is written at the controls.** HG002's haplotypes stay exactly spellable at
every control (70/70 matched, 43-44/44 correct). But the truth written through the realigned graphs
no longer matches the truth VCF's records.

| pooled truvari F1, raw / refined | mc | FFT-NS-2 | L-INS-i | E-INS-i | G-INS-i | abPOA | spoa | abPOA-MC | unit_aware | unit_aware__all |
|---|---|---|---|---|---|---|---|---|---|---|
| matched controls (40) | 0.958 / 0.958 | 0.584 / 0.832 | 0.771 / 0.954 | 0.732 / 0.929 | 0.722 / 0.942 | 0.648 / 0.881 | 0.405 / 0.746 | 0.811 / 0.958 | 0.323 / 0.903 | 0.304 / 0.744 |
| correct controls (25) | 0.985 / 0.985 | 0.735 / 0.989 | 0.839 / 0.971 | 0.822 / 0.985 | 0.855 / 0.971 | 0.812 / 0.941 | 0.675 / 0.924 | 0.913 / 0.985 | 0.527 / 0.944 | 0.428 / 0.874 |

- Whole-span phab F1 is 1.000 for every method except FFT-NS-2 (0.967), so the sequences are right.
- The truth VCF writes some repeat alleles as compensating SV pairs. At
  [L001894](../regions/L001894/) it has a 76 bp deletion and a 72 bp insertion 300 bp apart, and a
  60 bp deletion next to a 48 bp insertion. mc's graph reproduces those records. L-INS-i, E-INS-i,
  G-INS-i, abPOA, spoa and unit_aware write the same haplotype, which is an exact panel path, with
  no SV-sized record (abPOA-MC writes one), so truvari counts 3 FN. Refine does not recover them.
- unit_aware goes the other way: it writes 55 SV records for the 34 truth SVs at matched controls,
  14 of them TP (mc: 37, 34 TP).
- A genotyper on these graphs would inherit the representation. This is a benchmarking and
  representation issue rather than an alignment error, but Stage 3's raw F1 will show it.

## (c) Which method is best, and where does each fail?

**Hap32 arm, on the 135 regions where mc and all eight methods have a graph** (hotspot VNTRs n=27,
because spoa is missing 13 long ones). Cells are medians, except "truth edits", which is summed over
regions.

| hotspot VNTR (27) | mc | FFT-NS-2 | L-INS-i | E-INS-i | G-INS-i | abPOA | spoa | abPOA-MC | unit_aware |
|---|---|---|---|---|---|---|---|---|---|
| cost/opt (regions <= 1.1) | 1.599 (1) | 1.288 (3) | 1.114 (11) | 1.186 (7) | 1.126 (9) | 1.083 (17) | 1.067 (20) | 1.176 (8) | **1.053 (22)** |
| affine cost/opt | 1.571 | 1.246 | **1.115** | 1.173 | **1.115** | 1.119 | 1.135 | 1.173 | 1.117 |
| SV pieces/path | 8.49 | 6.12 | **3.67** | 5.09 | 4.61 | 4.00 | 5.15 | 4.46 | 5.36 |
| k-mer extra | 0.579 | 0.537 | 0.416 | **0.393** | 0.415 | 0.448 | 0.515 | 0.413 | 0.526 |
| parallel hom. bp/kb | 24.0 | 4.59 | 1.58 | 2.81 | 1.29 | 1.02 | 0.87 | 4.33 | **0.62** |
| MAPQ<5 reads with excess placements | 0.820 | 0.746 | 0.547 | **0.505** | 0.547 | 0.661 | 0.709 | 0.620 | 0.695 |
| truth edits, sum | **370** | 441 | 626 | 584 | 651 | 971 | 712 | 913 | 641 |
| best or within 0.005 on cost/opt | - | 0 | 0 | 1 | 0 | 1 | 6 | 1 | **24** |

| hotspot other (24) | mc | FFT-NS-2 | L-INS-i | E-INS-i | G-INS-i | abPOA | spoa | abPOA-MC | unit_aware |
|---|---|---|---|---|---|---|---|---|---|
| cost/opt (regions <= 1.1) | 1.625 (3) | 1.313 (6) | 1.164 (8) | 1.221 (6) | 1.163 (8) | 1.150 (9) | 1.125 (10) | 1.284 (6) | **1.053 (20)** |
| affine cost/opt | 1.598 | 1.258 | **1.138** | 1.162 | **1.138** | 1.197 | 1.230 | 1.176 | 1.218 |
| SV pieces/path | 4.70 | 3.35 | 2.44 | 2.76 | 2.59 | 3.00 | 2.94 | 2.80 | **1.56** |
| truth edits, sum | **240** | 285 | 406 | 471 | 357 | 793 | 869 | 658 | 542 |
| best or within 0.005 on cost/opt | - | 0 | 0 | 0 | 0 | 1 | 5 | 0 | **21** |

- **Repeat-unit-aware (unit_aware)** is the only method that reaches the target at most hotspots: 33/40
  hotspot VNTRs and 20/24 other hotspots at <= 1.1.
  - It also leaves the least homologous sequence on parallel nodes.
  - Its objective is linear in whole-unit gaps. Where a deletion goes is then free, so it scatters unit
    gaps. The gap-affine ratio is no better than L-INS-i's.
  - At hotspot VNTRs it writes more SV pieces (5.4 against 3.7) and more k-mer and read redundancy
    (identical units on different columns).
  - At controls it writes the truth as scattered pieces (section b).
- **Progressive MSA (mafft):**
  - **L-INS-i and G-INS-i** give the best gap-affine alignments and the fewest SV pieces, and remove
    the most read redundancy along with E-INS-i. They are 0.06-0.07 behind unit_aware on unit cost at
    hotspots, and they regress 3 non-repeat SV controls.
  - **E-INS-i** leaves more unaligned sequence (2.9 bp/kb) and more SV pieces, but has the lowest k-mer
    and read redundancy.
  - **FFT-NS-2** is not usable. It recreates the parallel-path fault (L012184 1.86, worse than mc) and
    stacks units out of register (L002143 9.51, median pair 16x the optimum). It regresses a quarter of
    the VNTR controls by more than 0.05.
- **Partial-order alignment:**
  - **spoa and abPOA** have the lowest unit cost after unit_aware at hotspots (1.067, 1.083) and leave
    little parallel sequence. But under affine scoring they fall behind L-INS-i at hotspot_other (1.23
    and 1.20 against 1.14).
  - With abPOA-MC they lose the most truth spellability (hotspot edits 971, 712 and 913 against mc's
    370), and both fail by memory on the longest alleles.
  - **abPOA with Minigraph-Cactus's own settings (poa_abpoa_mc)** is between mc and default abPOA
    (1.170). Cactus's scoring is not what makes MC's hotspots bad, but it does not fix them.

**Where each method fails:**

| failure mode | affected regions | methods that fail there |
|---|---|---|
| very long regions and outlier alleles (memory or time) | [L012272](../regions/L012272/) (55 kb span, 115 kb longest full-panel allele), [L015347](../regions/L015347/), [L001909](../regions/L001909/), [L007040](../regions/L007040/), [L011138](../regions/L011138/), [L012184](../regions/L012184/) (full panel: one 39.6 kb allele) | spoa (14 hap32, 26 full-panel regions missing), abPOA (3 / 11), all-pairs mafft modes (1-3 / 72-84), FFT-NS-2 full panel (9) |
| BAR 10 kb rule leaves the middle unaligned | L012272 (hap32 arm U = 607 bp/kb, cost/opt 1.467 against mc 1.254), L012184 full panel (1.640 projected) | `poa_abpoa_mc` only (by design) |
| compound arrays and differently ordered units | [L011463](../regions/L011463/) (60 bp unit, 234 unit symbols: best 1.241), [L009658](../regions/L009658/) (1.175), [L004145](../regions/L004145/) (8 bp, 15.6 kb; unit_aware fell back: 1.228, best abPOA 1.142), [L000034](../regions/L000034/) (1.127) | all: 10 of 64 hotspots stay above 1.1 with every method |
| non-repeat SVs | [L006556](../regions/L006556/) (1.00 to 1.19, affine 1.37, 153 against 16 nodes/kb), L002292, L003883 | L-INS-i, G-INS-i, unit_aware (fallback); E-INS-i and POA are exact |
| short-period correct controls | [L011224](../regions/L011224/) (10 bp CA repeat: 1.027 to 1.137 / 1.115) | abPOA, spoa |
| N-gap haplotype | L012184 (`recombination#29`, 3,146 bp of N) | none: every method makes it one N node |
| fragments (paths ending inside the array) | 7/40 hotspot VNTRs and 1/40 matched controls have `hap32.fragments.fa` records | untested: no graph includes fragments, and every metric uses spanning paths |

Cost/opt by span is flat for unit_aware (hotspots of both strata: 1.053 / 1.048 / 1.048 for spans
< 5 kb, 5-10 kb and > 10 kb; n = 39 / 16 / 9). abPOA goes 1.091 / 1.079 / 1.117 but loses 3 of the 9
longest; spoa has only 1 of the 9.

## (d) The panel question: m(hap32) vs m__all vs mc

**Projected graphs, scored on hap32.** Paired over the regions where mc, m(hap32) and m__all all
exist.
- **Kept** = (mc - m__all) / (mc - m(hap32)) on the medians, shown only where m(hap32) gained at
  least 0.02.
- **all vs hap32** = regions where m__all is better / worse than m(hap32) by more than 0.005 (p is
  the Wilcoxon test).

| method | stratum (n) | cost/opt mc / m / m__all | kept | all vs hap32 b / w (p) | m__all vs mc b / w | <= 1.1 mc / m / m__all | affine mc / m / m__all |
|---|---|---|---|---|---|---|---|
| `unit_aware` | hotspot VNTR (40) | 1.566 / 1.048 / 1.059 | 98% | 2 / 34 (2e-7) | 40 / 0 | 2 / 33 / 31 | 1.475 / 1.093 / 1.138 |
| | hotspot other (24) | 1.625 / 1.053 / 1.097 | 92% | 0 / 22 (1e-7) | 24 / 0 | 3 / 20 / 12 | 1.598 / 1.218 / 1.312 |
| | matched control (40) | 1.102 / 1.015 / 1.045 | 66% | 0 / 23 (2e-8) | 22 / 8 | 19 / 37 / 31 | 1.037 / 1.020 / 1.072 |
| `poa_abpoa` | hotspot VNTR (30) | 1.582 / 1.083 / 1.130 | 91% | 1 / 22 (1e-6) | 30 / 0 | 1 / 19 / 11 | 1.498 / 1.117 / 1.238 |
| | hotspot other (23) | 1.653 / 1.164 / 1.244 | 84% | 1 / 20 | 21 / 2 | 2 / 8 / 3 | 1.603 / 1.207 / 1.476 |
| | matched control (40) | 1.102 / 1.024 / 1.084 | 23% | 1 / 22 | 15 / 12 | 19 / 33 / 20 | 1.037 / 1.006 / 1.067 |
| `poa_spoa` | hotspot VNTR (20) | 1.599 / 1.050 / 1.115 | 88% | 1 / 19 (5e-4) | 19 / 1 | 1 / 14 / 8 | 1.544 / 1.119 / 1.299 |
| | hotspot other (22) | 1.658 / 1.127 / 1.239 | 79% | 0 / 22 | 20 / 2 | 2 / 8 / 3 | 1.604 / 1.234 / 1.529 |
| | matched control (37) | 1.103 / 1.016 / 1.069 | 40% | 0 / 25 | 19 / 10 | 17 / 32 / 26 | 1.038 / 1.018 / 1.099 |
| `poa_abpoa_mc` | hotspot VNTR (39) | 1.567 / 1.170 / 1.275 | 73% | 6 / 32 (2e-7) | 29 / 10 | 2 / 11 / 6 | 1.475 / 1.146 / 1.317 |
| | hotspot other (24) | 1.625 / 1.284 / 1.356 | 79% | 5 / 19 | 18 / 5 | 3 / 6 / 3 | 1.598 / 1.176 / 1.343 |
| | matched control (40) | 1.102 / 1.031 / 1.084 | 26% | 2 / 19 | 13 / 11 | 19 / 31 / 20 | 1.037 / 1.008 / 1.038 |
| `mafft_fftns2` | hotspot VNTR (31) | 1.598 / 1.288 / 1.420 | 57% | 5 / 26 (1e-3) | 22 / 9 | 1 / 3 / 0 | 1.521 / 1.249 / 1.369 |
| | hotspot other (24) | 1.625 / 1.313 / 1.360 | 85% | 1 / 21 | 19 / 3 | 3 / 6 / 1 | 1.598 / 1.258 / 1.300 |
| | matched control (40) | 1.102 / 1.100 / 1.206 | - | 7 / 23 | 12 / 21 | 19 / 20 / 15 | 1.037 / 1.057 / 1.151 |
| `mafft_linsi` | hotspot VNTR (**4**) | 1.307 / 1.042 / 1.035 | 102% | 2 / 1 | 4 / 0 | 1 / 3 / 3 | 1.324 / 1.048 / 1.040 |
| | hotspot other (6) | 1.518 / 1.158 / 1.159 | 100% | 2 / 4 | 6 / 0 | 1 / 2 / 1 | 1.489 / 1.136 / 1.153 |
| | matched control (27) | 1.112 / 1.022 / 1.054 | 64% | 2 / 12 (3e-3) | 15 / 3 | 12 / 23 / 20 | 1.036 / 1.004 / 1.022 |
| `mafft_einsi` | hotspot VNTR (**1**) / other (4) / matched (22) | 1.431 / 1.050 / 1.045; 1.518 / 1.210 / 1.330; 1.106 / 1.016 / 1.018 | | 1 / 0; 0 / 4; 1 / 3 | | | |
| `mafft_ginsi` | hotspot VNTR (**1**) / other (4) / matched (22) | 1.431 / 1.035 / 1.033; 1.518 / 1.158 / 1.163; 1.106 / 1.014 / 1.055 | | 0 / 0; 0 / 2; 0 / 8 | | | |

- **Only unit_aware (40) and abPOA-MC (39) cover nearly all hotspot VNTRs in the projected arm.**
  - FFT-NS-2 covers 31. abPOA and spoa miss the 10 and 20 with the longest alleles.
  - L-/E-/G-INS-i cover 4, 1 and 1. Their "gain kept" of about 100% is on small, easy regions (mc
    median 1.307 at those 4, against 1.566 overall), so it says nothing about hotspots in general.
  - FFT-NS-i has no hap32 arm. Against mc only: hotspot VNTRs 1.335 to 1.173 (n=9, 7 better / 2
    worse); matched controls 1.102 to 1.111 (n=34, 12 / 14).
- **The answer, per method:**
  - unit_aware keeps almost all of its gain at hotspots.
  - POA keeps most of it (84-91% abPOA, 79-88% spoa), with a larger loss under affine scoring
    (abPOA 1.117 to 1.238).
  - abPOA-MC and FFT-NS-2 lose a quarter to a half.
  - No method improves on its hap32 arm at hotspots. The projection is worse in 19-34 of 20-40
    regions for every method that covers them.
- **The controls pay the most** (matched controls, kept 23-66%). The full panel's columns were placed
  for about 450 sequences, so the hap32 rows inherit gaps cut into pieces.

**The full panel itself** (`results/full/`, 300 sampled pairs per region, the same pairs for every
graph of a region).

| stratum (n) | cost/opt, mc.gfa on hap32 / full MC on the full panel | excess edits/kb | parallel hom. bp/kb | k-mer extra | nodes/kb | <= 1.1 |
|---|---|---|---|---|---|---|
| hotspot VNTR (40) | 1.566 / 1.513 (18 better, 21 worse) | 182 / 242 | 40.5 / 55.0 | 0.601 / 0.784 | 764 / 1,370 | 2 / 2 |
| matched control (40) | 1.102 / 1.199 (14 / 18) | 4.4 / 8.9 | 0.15 / 0.40 | 0.053 / 0.159 | 101 / 138 | 19 / 15 |
| correct control (25) | 1.015 / 1.056 (7 / 15) | 2.0 / 6.9 | 0.28 / 0.28 | 0.041 / 0.141 | 86 / 112 | 17 / 15 |
| hotspot other (24) | 1.625 / 1.427 (11 / 12) | 96 / 107 | 22.7 / 23.6 | 0.435 / 0.654 | 492 / 875 | 3 / 1 |
| non-TR SV control (20) | 1.000 / 1.000 | 0.005 / 0.025 | 0 / 0 | 0.007 / 0.006 | 49 / 56 | 20 / 20 |

So the full MC graph is as far from optimal on the whole panel as mc.gfa is on hap32. At hotspot
VNTRs the ratio is similar (1.51 against 1.57), the absolute excess is larger, and the controls are
worse. Hap32's view of MC is not a sampling artefact. (mc.gfa is exactly the full graph restricted
to the hap32 walks; see tools/README.md, "premise".)

**Realigned full-panel graphs against the full MC graph** (paired; b / w by more than 0.005):

| method | hotspot VNTR: n, full MC -> m, b / w, <= 1.1 | matched control | correct control | hotspot other |
|---|---|---|---|---|
| `unit_aware` | 40: 1.513 -> **1.054**, 40 / 0, 2 -> 32 | 40: 1.199 -> 1.054, 28 / 8 | 25: 1.056 -> 1.023, 18 / 3 | 24: 1.427 -> 1.063, 24 / 0 |
| `poa_spoa` | 20: 1.677 -> 1.094, 20 / 0, 1 -> 11 | 37: 1.218 -> 1.090, 22 / 6 | 24: 1.055 -> 1.116, 11 / 9 | 22: 1.456 -> 1.183, 22 / 0 |
| `poa_abpoa` | 30: 1.624 -> 1.111, 29 / 1, 1 -> 13 | 40: 1.199 -> 1.114, 20 / 11 | 25: 1.056 -> 1.109, 9 / 9 | 23: 1.449 -> 1.205, 22 / 1 |
| `poa_abpoa_mc` | 39: 1.515 -> 1.269, 36 / 3, 2 -> 8 | 40: 1.199 -> 1.149, 15 / 13 | 25: 1.056 -> 1.045, 11 / 4 | 24: 1.427 -> 1.324, 19 / 4 |
| `mafft_fftns2` | 31: 1.610 -> 1.361, 28 / 3, 1 -> 1 | 40: 1.199 -> **1.309**, 13 / 22 | 25: 1.056 -> 1.118, 6 / 12 | 24: 1.427 -> 1.275, 20 / 3 |
| `mafft_fftnsi` | 9: 1.592 -> 1.289, 9 / 0 | 34: 1.210 -> 1.243, 14 / 13 | 25: 1.056 -> 1.103, 7 / 11 | 17: 1.463 -> 1.251, 16 / 1 |
| `mafft_linsi` | **4**: 1.515 -> 1.040, 4 / 0 | 27: 1.196 -> 1.083, 17 / 5 | 20: 1.070 -> 1.043, 10 / 4 | 6: 1.286 -> 1.134, 6 / 0 |
| `mafft_einsi`, `mafft_ginsi` | **1**: 1.592 -> 1.055 / 1.045 | 22: 1.149 -> 1.068 / 1.059 | 18: 1.070 -> 1.037 / 1.048 | 4: 1.207 -> 1.192 / 1.088 |

- **Hotspot VNTRs:** on the full panel the realigned graphs close more than all of the gap to the
  matched controls. unit_aware reaches 1.054, and full MC's matched controls are 1.199.
- **Parallel homologous sequence at hotspot VNTRs** goes from 36-55 bp/kb (full MC, on each method's
  regions) to 0.6 (unit_aware), 2.7-2.9 (POA) and 12-13 (FFT-NS-2, abPOA-MC).
- **The other structural metrics do not close:**
  - SV pieces per haplotype: 8.6 to 7.3 for unit_aware, 3.3 for spoa;
  - k-mer extra: about 0.7 for everyone;
  - nodes/kb: 1,370 to 1,449 for unit_aware; the POA graphs fall to about 1,070.
- **At controls:** unit_aware and the INS-i modes improve on full MC. abPOA and spoa are mixed. FFT-NS-2
  is worse at the matched and correct controls, and FFT-NS-i is no better.

**Seed check.** The full-panel metrics sample 300 pairs. I re-scored full MC, unit_aware,
abPOA-MC and abPOA on hotspot VNTRs and matched controls with `--seed 2`: 310 graphs, written to
scratch, not to `results/`.

| method | stratum (n) | seed 1: full MC -> m, b / w, <= 1.1 | seed 2: full MC -> m, b / w, <= 1.1 | same sign vs MC |
|---|---|---|---|---|
| `unit_aware` | hotspot VNTR (40) | 1.513 -> 1.054, 40 / 0, 32 | 1.603 -> 1.054, 40 / 0, 32 | 40 / 40 |
| | matched control (40) | 1.199 -> 1.054, 28 / 8, 31 | 1.154 -> 1.056, 29 / 6, 30 | 39 / 40 |
| `poa_abpoa_mc` | hotspot VNTR (39) | 1.515 -> 1.269, 36 / 3, 8 | 1.615 -> 1.320, 37 / 2, 9 | 38 / 39 |
| | matched control (40) | 1.199 -> 1.149, 15 / 13, 16 | 1.154 -> 1.135, 16 / 13, 18 | 37 / 40 |
| `poa_abpoa` | hotspot VNTR (30) | 1.624 -> 1.111, 29 / 1, 13 | 1.650 -> 1.115, 29 / 1, 13 | 30 / 30 |
| | matched control (40) | 1.199 -> 1.114, 20 / 11, 17 | 1.154 -> 1.108, 24 / 10, 19 | 39 / 40 |

- **Full MC's own value moves:** the median per-region change between seeds is 0.036 at hotspots, up
  to 0.47. The realigned graphs move by 0.003-0.019.
- **The comparisons hold:** every conclusion stands, and the better / worse and <= 1.1 counts move by
  at most 4 regions.
- **So the full-panel MC medians carry about ±0.05 of sampling noise at hotspots** (1.51 with seed 1,
  1.60 with seed 2), and the "as far from optimal as on hap32" statement holds under both seeds.

## (e) Truth in the graph

Hotspot VNTRs, 80 HG002 haplotypes (40 regions); each method is paired with mc on its own regions.

| | mc | FFT-NS-2 | L-INS-i | E-INS-i | G-INS-i | abPOA | spoa | abPOA-MC | unit_aware | unit_aware__all | abPOA-MC__all |
|---|---|---|---|---|---|---|---|---|---|---|---|
| regions | 40 | 40 | 39 | 37 | 39 | 37 | 27 | 40 | 40 | 40 | 39 |
| haplotypes spelled exactly (mc's exact ones kept) | 39 | 36 (32/39) | 29 (28/39) | 28 (27/39) | 30 (29/39) | 29 (28/38) | 27 (26/31) | 28 (28/39) | 38 (34/39) | 43 (36/39) | 35 (31/39) |
| of mc's 22 recombinant exact spellings, kept | 22 | 15 | 11 | 10 | 12 | 12 | 13 of 18 | 11 | 17 | 19 | 14 |
| summed edits, mc -> method | 815 | 826 | 812 -> 1,042 | 550 -> 845 | 812 -> 1,080 | 604 -> 1,523 | 370 -> 712 | 815 -> 1,845 | 815 -> 1,071 | 815 -> **670** | 812 -> 981 |
| raw F1, mc -> method | 0.262 | 0.319 | 0.259 -> 0.337 | 0.270 -> 0.382 | 0.259 -> 0.324 | 0.273 -> 0.327 | 0.301 -> 0.347 | 0.299 | 0.260 | 0.265 | 0.259 -> 0.284 |
| refined F1, mc -> method | 0.867 | 0.951 | 0.864 -> 0.964 | 0.906 -> 0.964 | 0.864 -> 0.918 | 0.909 -> 0.983 | 0.894 -> 0.978 | 0.919 | 0.941 | 0.946 | 0.864 -> 0.932 |
| phab F1 | 0.989 | 0.997 | 1.000 | 0.998 | 1.000 | 1.000 | 0.997 | 0.999 | 0.996 | 1.000 | 1.000 |

- **Exact panel copies** (17 of the 80 haplotypes) stay exactly spellable in every graph, as they
  must.
- **Recombinant spellings:** mc spells 22 more exactly only by recombining through its mesh of 1-bp
  nodes. Column-induced graphs allow fewer recombinant paths and lose 5-12 of those 22.
  - So summed edits rise in the hap32 arm, driven by a few loci: [L005990](../regions/L005990/) 21 to
    75-116, [L001909](../regions/L001909/) 55 to 114-280, [L002013](../regions/L002013/) 0 to 24
    (L-INS-i).
  - Every graph stays 8-23x closer to the truth than the closest single panel path (15,596
    edits).
- **Projection:** the projected arm spells the truth better than the hap32 arm. unit_aware__all has
  670 edits against mc's 815, and FFT-NS-2__all 236 against 404 (n=31). The full panel's gap
  placements leave more recombinant routes.
- **Refined F1** rises at hotspots for every method (0.92-0.98 against 0.86-0.91).
- **Raw F1** rises modestly for the mafft and POA modes (by 0.04-0.11) and not for unit_aware (0.260
  against 0.262).
  - At hotspot_other it falls for spoa (0.31 to 0.19) and unit_aware (0.20).
  - unit_aware writes more SV records than the truth has: 995 against 558 at hotspot VNTRs, against
    mc's 743.
  - L-INS-i writes fewer: 418 for 543.

## (f) Stage 1: placement redundancy of the existing reads

MAPQ<5 core reads, pooled over reads, mc -> method, paired regions.

| method | hotspot VNTR: n regions / reads | unique | redundant only | any excess placement | regions redundant-only fell / rose | matched control: redundant only |
|---|---|---|---|---|---|---|
| `mafft_fftns2` | 40 / 23,465 | 13.4 -> 20.4% | 29.7 -> 22.7% | 82.5 -> 75.5% | 31 / 9 | 19.5 -> 16.8% |
| `mafft_linsi` | 39 / 22,726 | 13.7 -> 24.9% | 28.9 -> 17.7% | 82.2 -> 67.1% | 37 / 2 | 19.5 -> 11.3% |
| `mafft_einsi` | 37 / 19,262 | 13.4 -> 30.4% | 31.0 -> **14.1%** | 82.6 -> 58.7% | 36 / 1 | 19.5 -> 10.0% |
| `mafft_ginsi` | 39 / 22,726 | 13.7 -> 25.0% | 28.9 -> 17.6% | 82.2 -> 67.0% | 37 / 2 | 19.5 -> 11.3% |
| `poa_abpoa` | 37 / 18,601 | 14.5 -> 23.2% | 29.7 -> 20.9% | 81.1 -> 70.5% | 35 / 2 | 19.5 -> 12.6% |
| `poa_spoa` | 27 / 10,962 | 15.5 -> 22.2% | 35.4 -> 28.7% | 79.8 -> 72.9% | 24 / 3 | 19.8 -> 15.1% |
| `poa_abpoa_mc` | 40 / 23,465 | 13.4 -> 24.5% | 29.7 -> 18.6% | 82.5 -> 67.8% | 39 / 1 | 19.5 -> 11.9% |
| `unit_aware` | 40 / 23,465 | 13.4 -> 19.9% | 29.7 -> 23.3% | 82.5 -> 75.5% | 33 / 6 | 19.5 -> 15.7% |
| `unit_aware__all` | 40 / 23,465 | 13.4 -> 18.4% | 29.7 -> 24.7% | 82.5 -> 77.6% | 29 / 10 | 19.5 -> 16.4% |
| `poa_abpoa_mc__all` | 39 / 22,726 | 13.7 -> 18.6% | 28.9 -> 24.0% | 82.2 -> 76.8% | 33 / 5 | 19.5 -> 15.1% |
| `mafft_fftns2__all` | 31 / 12,493 | 14.5 -> 20.5% | 33.5 -> 27.6% | 80.9 -> 75.1% | 22 / 9 | 19.5 -> 19.3% |
| `poa_abpoa__all` | 30 / 11,834 | 14.6 -> 15.5% | 34.3 -> 33.4% | 80.8 -> 81.4% | 16 / 12 | 19.5 -> 17.4% |
| `poa_spoa__all` | 20 / 7,548 | 14.9 -> 14.1% | 41.5 -> 42.3% | 80.8 -> 83.8% | 8 / 11 | 19.0 -> 19.3% |

- **Tandem ambiguity** is 56% of hotspot MAPQ<5 reads in every graph. It depends only on the
  haplotype sequences, so no alignment can remove it, as the investigation predicted.
- **Most of the removable class goes:** the best graph (E-INS-i) removes about half the redundant
  class, and L-INS-i, G-INS-i and abPOA-MC remove about 40%.
- **unit_aware removes the least of the hap32-arm methods apart from spoa** (22% of the class, against
  24% for FFT-NS-2 and 39-55% for the INS-i modes). It stacks units in phase but puts identical unit
  variants in different columns, so the same read sequence sits on several nodes.
- **The projected POA arms remove none.**
- **Matched controls follow the same pattern:** 19.5% redundant only, down to 10-17%.
- **Why a placement ratio can rise:** it is computed with start/end offsets. A compacted graph with
  long nodes holding several tandem copies can raise it, as FFT-NS-2 does (3.17 to 3.39) while
  redundant-only falls.

## Per-locus observations

1. **[L012184](../regions/L012184/)** (chr4, parallel duplicated path plus an N-gap haplotype).
   - mc has cost/opt 1.599 with 41.2 bp/kb of parallel homologous sequence.
   - Every realigner cuts parallel homologous sequence to 0.4-5.8 bp/kb. abPOA and spoa give 1.091
     and 1.081, unit_aware 1.136. FFT-NS-2 is still worse than mc (1.860).
   - unit_aware spells both HG002 haplotypes exactly (mc: 6 edits).
   - The 3,146 bp N run is one N node in every graph.
   - On the full panel, abPOA and spoa cannot run: one allele is 39.6 kb against CHM13's 4.8 kb.
     abPOA-MC's projection is worse than mc (1.640; BAR left 140 kb unaligned). unit_aware__all
     gives 1.183 projected and 1.055 on the full panel, where full MC scores 1.468.
2. **[L002143](../regions/L002143/)** (chr11, 69 bp unit, 101 copies).
   - FFT-NS-2 scores cost/opt 9.51. The median pair costs 16x its optimum, with 594 nodes/kb, so the
     units are stacked out of register. FFT-NS-2__all scores 5.01.
   - Every other method scores 1.08-1.38 (unit_aware 1.083), against mc's 2.404.
3. **[L011138](../regions/L011138/) and [L012272](../regions/L012272/)** (under-aligned hotspots;
   mc has 340 and 457 bp/kb of parallel homologous sequence).
   - unit_aware removes it (0.04 and 0.01 bp/kb) and gives 1.339 to 1.016 and 1.254 to 1.048.
   - L012272 (55 kb span, 42 kb of it anchor flank):
     - Only FFT-NS-2 (1.325), abPOA-MC and unit_aware finish in the hap32 arm.
     - abPOA-MC scores 1.467, worse than mc, with 607 bp/kb unaligned by the 10 kb rule.
     - Only unit_aware__all runs on the full panel (1.081).
4. **[L005990](../regions/L005990/)** (chr17, a block of about 38 units near the start of the array).
   - Cost/opt improves from 1.211 to 1.014 (unit_aware) and 1.063 (L-INS-i).
   - HG002's edit distance to the graph rises from 21 to 109 and 75, and spoa reaches 116.
   - unit_aware's affine ratio is 1.175 against L-INS-i's 1.063: its unit gaps are scattered.
   - The projections recover part of the truth distance (FFT-NS-2__all 9, abPOA__all 23,
     unit_aware__all 71).
5. **[L014297](../regions/L014297/)** (chr6, 17 bp unit, 153 CHM13 copies).
   - Every hap32 arm except FFT-NS-2 (1.288) scores 1.02-1.09 (mc 1.303).
   - FFT-NS-2__all scores 1.931 with 31.4 bp/kb of parallel sequence: units aligned only to
     haplotypes the projection removed.
   - unit_aware__all scores 1.056.
6. **[L001894](../regions/L001894/)** (matched control).
   - Both HG002 haplotypes are exact panel paths in every graph, and cost/opt improves from 1.273 to
     1.003-1.059.
   - Yet the truth written by six of the realigned graphs has no SV-sized record (abPOA-MC: one;
     FFT-NS-2: four), where the truth VCF has a 76 bp deletion plus a 72 bp insertion, and a 60 bp
     deletion plus a 48 bp insertion.
   - Raw F1 goes from 1.0 to 0.0, and refine does not recover it. This is the representation shift
     behind the control F1 drop in (b).
7. **[L006556](../regions/L006556/)** (non-repeat SV control).
   - L-INS-i and G-INS-i, and unit_aware through its fallback, turn mc's exact alignment (1.003) into
     1.189. Affine is 1.365, with 153 nodes/kb against 16.
   - E-INS-i and abPOA reproduce mc exactly. A non-repeat SV needs a global or generalized-affine
     aligner, not a local-pair one.
8. **[L011224](../regions/L011224/)** (correct control, 10 bp CA repeat).
   - abPOA and spoa regress from 1.027 to 1.137 and 1.115 (affine 1.29 and 1.35).
   - unit_aware's adequacy check switched it to L-INS-i (1.036).
9. **[L000034](../regions/L000034/), [L011463](../regions/L011463/), [L009658](../regions/L009658/)**
   (arrays whose unit types are ordered differently between haplotypes, or which have several
   motifs; UNIT_ALIGN.md names L000034).
   - No method reaches 1.1: best 1.127, 1.241 and 1.175, all unit_aware, against mc's 2.73, 2.17 and
     2.95.
   - UNIT_ALIGN.md shows that L000034's differently ordered units cannot be stacked acyclically. I
     did not check the other two for the same cause.
10. **[L001909](../regions/L001909/)** (20 kb; 22 of its 32 hap32 sequences are recombinant junctions
    of two full-panel haplotypes).
    - Cost/opt improves from 1.959 to 1.107 (unit_aware) and 1.159 (L-INS-i).
    - HG002's edit distance rises from 55 to 114-280.
    - Only unit_aware and abPOA-MC run it on the full panel (1.115 and 1.143 projected; 1.039 and
      1.088 on the full panel, against full MC's 1.841).
11. **[L016870](../regions/L016870/)** (chrX negative control from the plan).
    - Every hap32 arm scores 1.009-1.032 against mc's 1.030.
    - abPOA-MC__all scores 1.070 (BAR left 25.5 kb unaligned).

## Caveats

- **Unit cost and affine cost disagree on which method is best.** Cost/opt uses unit edit cost,
  which is blind to where a gap goes. unit_aware's win shrinks or reverses under gap-affine scoring,
  SV pieces, k-mer and read redundancy, and raw truvari F1. The plan's target (cost/opt <= 1.1) is met
  best by unit_aware, while L-INS-i/G-INS-i and abPOA-MC give tidier gaps.
- **Many equal-cost optima.** A repeat MSA is one of many alignments with the same cost. The graph
  metrics judge the one each method returned.
- **Full-panel scores are sampled** (300 pairs, seed 1). The same pairs are used for every graph of a
  region, so method-vs-MC differences are paired. Absolute values at hotspots move between seeds: full
  MC's hotspot median is 1.513 with seed 1 and 1.603 with seed 2, while the realigned graphs move by
  less than 0.02. The seed check in (d) shows that the comparisons hold.
- **Biased subsets.**
  - The all-pairs mafft `__all` arms cover only small and medium regions (4 and 1 hotspot VNTRs).
  - spoa's hap32 arm misses 13 long hotspot VNTRs; abPOA misses 3 in its hap32 arm and 10 in `__all`.
  - The common-subset table (c) has 27 hotspot VNTRs, not 40.
- **nodes/kb and node length cannot be compared across mc and MSA graphs.** mc keeps cuts made by
  variation in the full panel (48.7 nodes/kb at non-repeat SV controls, against 17-23 for every MSA
  graph).
- **Raw truvari F1 at controls is a representation score.** Its fall at the VNTR controls comes with
  exact spelling and phab F1 1.0 (section b).
- **Read pooling** weights each region by its MAPQ<5 read count.
  - mc's hotspot numbers here are 13.4% unique, 56.0% tandem and 29.7% redundant only (n=23,465).
    The metrics agent reported 12.7 / 57.6 / 28.8 (n=22,838).
  - The per-region values are identical to that run (all 38 summarised metrics match on all 149
    regions), so the difference is in how that figure was pooled.
- **Existing results.**
  - `results/mc/` was recomputed with the current evaluate.py (identical per-region metrics).
  - `results/poa_abpoa_mc*` were kept from an earlier run with the same code.
  - `results/panel/{mc,poa_abpoa_mc}` (the panel agents' output) are identical to
    `results/full/{mc,poa_abpoa_mc}` in every metric section, and are superseded by `results/full/`.
- **Stage 1 only.** Stage 1 tests where the existing reads' sequences would land. It is not a
  re-mapping; that is Stage 2.
- **Code changes.**
  - `tools/evaluate_all.py` gained `--panel-subdir` (default `panel`, unchanged behaviour) so that
    full-panel results could go to `results/full/`.
  - `tools/stage01.py` is new, and `tools/README.md` has a section for it.
  - evaluate.py and summarise.py are unchanged.

## Commands

```bash
cd ~/PycharmProjects/vntr-eval
python3 tools/evaluate_all.py mc --reads --force --jobs 3 --threads 1                 # baseline
python3 tools/evaluate_all.py $(ls candidates) --reads --jobs 3 --threads 1            # both arms
python3 tools/evaluate_all.py mc $(ls work/panel | grep -v -e query -e union -e mc) \
    --panel --panel-subdir full --jobs 1 --threads 1                                    # full-panel graphs
python3 tools/summarise.py --methods mc,<all 17> --metrics headline --md results/summary.md
python3 tools/summarise.py --methods mc,<all 17> --metrics all --tsv results/summary.tsv --per-region results/per_region.tsv
python3 tools/summarise.py --results results/full --methods mc,<9 methods> --md results/full/summary.md --tsv results/full/summary.tsv
python3 tools/stage01.py                                  # -> results/stage01_tables.md, results/stage01_pairs.tsv
```
