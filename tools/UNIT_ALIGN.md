# M4: repeat-unit-aware two-level realignment (`tools/realign_units.py`)

M4 of [docs/plan.md](../docs/plan.md). It re-aligns a region's haplotype sequences so that the
repeat units of different haplotypes stack in phase, then builds the region's graph from that
alignment with [msa_graph.py](msa_graph.py), exactly as the mafft and POA realigners do.
The C code is [unit_dp.c](unit_dp.c), built on demand into `tools/bin/libunitdp.so`. The tests are
[test_realign_units.py](test_realign_units.py).

## Result on the 149 packaged regions (hap32)

`python3 tools/realign_units.py all --jobs 2 --threads 2` wrote 149 candidates, every one valid,
in 795 s wall-clock (median 3.2 s per region, 1,358 s summed).
- **Unit mode: 121 regions.**
- **Fallback: 28 regions.**
  - 20 non-TR controls;
  - 6 TR regions whose motif does not fit, where the CHM13 array covers < 50% of the core:
    L000709, TR101140, L003119, L000698, L015743, L013074. Most have periods of 225-1990 bp with
    2-8 copies.
  - 2 regions switched by the guard: L011224 (unit MSA 1.65x the pairwise optima, L-INS-i 1.06x)
    and L004145 (1.247x vs 1.221x, a near tie).

**Stage-0 scores.** evaluate.py was run without reads or truvari, into a scratch directory, against
`results/mc` (MC = the Minigraph-Cactus hap32 baseline). Medians, M4 / MC:

| stratum | n (unit) | cost/opt | affine | CHM13 pairs | unaligned homology bp/kb | SV pieces/path | nodes/kb | cost/opt <= 1.1 | better / worse than MC (> 0.005) |
|---|---|---|---|---|---|---|---|---|---|
| hotspot_vntr | 40 (39) | 1.048 / 1.566 | 1.093 / 1.475 | 1.040 / 1.670 | 0.35 / 40.5 | 5 / 8 | 381 / 764 | 33 / 2 | 40 / 0 |
| control_vntr_matched | 40 (38) | 1.015 / 1.102 | 1.020 / 1.037 | 1.012 / 1.069 | 0 / 0.15 | 0 / 1 | 40 / 101 | 37 / 19 | 24 / 3 |
| control_vntr_correct | 25 (20) | 1.015 / 1.015 | 1.024 / 1.012 | 1.011 / 1.012 | 0 / 0.28 | 1 / 1 | 26 / 86 | 25 / 17 | 13 / 4 |
| hotspot_other | 24 (24) | 1.053 / 1.625 | 1.218 / 1.598 | 1.050 / 1.565 | 0.60 / 22.7 | 0.5 / 5 | 398 / 492 | 20 / 3 | 24 / 0 |
| control_nontr_sv | 20 (0) | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 0 / 0 | 1 / 1 | 23 / 49 | 19 / 20 | 0 / 4 |

**Regressions.** Every unit-mode region where M4 is worse than MC is a control. The largest is
L010113 (1.038 vs 1.000); the others are at most +0.009, except L000381 (1.026 vs 1.017). The four
non-TR controls that regress (L006556 1.19, L002292 1.09, L003883 1.06, L005688 1.02) are plain
mafft L-INS-i, i.e. M2's own behaviour on a large SV, not the unit method.

**HG002 distance to the graph** (`h1_d_graph + h2_d_graph`, summed) gets worse at the hotspots:
1,071 vs 815 at hotspot VNTRs, and 542 vs 240 at hotspot_other. Per region it is mixed (22 worse,
23 better, 103 equal in the go = 0 test run) and driven by a few loci: L007844, L007081 and L001909
worse, L015347 better. MC's 1-bp mesh allows more recombinant walks; a column-induced graph
recombines only where rows share a column.

**Units and symbols.** Median periods and symbols per unit-mode region:

| stratum | period | symbols (distinct unit variants) | distinct unit strings |
|---|---|---|---|
| hotspot_vntr | 29 | 51 (10-297) | 29 |
| control_vntr_matched | 26 | 76 (23-374) | 8 |
| control_vntr_correct | 39 | 48 (12-97) | 8 |
| hotspot_other | 4 | 23 (7-148) | 27 |

Five regions have more than mafft's 248 symbols (L007376 374, L000381 335, L015347 297,
L014160 291, L010350 267). The shipped aligner has no alphabet limit.

**Phase.** In every unit-mode region, every unit column holds units that start at one unit
position (`phase_check` = 1.0). The spot checks in viewer.py's alignment panel (L012184, L005990,
L000034) show the inserted units stacked in shared columns across haplotypes, at the same phase.

**Full HPRC panel, MSA-only mode:**

| region | records | distinct sequences | unit strings | time |
|---|---|---|---|---|
| L014297 | 457 | 299 | 298 | 43 s |
| L012184 | 447 | 364 | 356 | 171 s |

L012184 was measured before the polish cap was raised to 3e7 cells; with the earlier cap its polish
skipped every row.

**How the method got here** (10 pilot hotspots, evaluate.py all-pairs cost/optimum):

| region | MC | mafft on units (best setting) | exact-SP progressive (leave-one-out) | + tree refinement + polish | mafft L-INS-i (M2) | abPOA |
|---|---|---|---|---|---|---|
| L012184 | 1.599 | 1.230 | 1.134 | 1.136 | 1.143 | 1.091 |
| L005990 | 1.211 | 1.096 | 1.020 | 1.014 | 1.063 | 1.044 |
| L009656 | 1.335 | 1.125 | 1.035 | 1.023 | 1.099 | 1.062 |
| L014297 | 1.303 | 1.062 | 1.022 | 1.022 | 1.051 | 1.067 |
| L015415 | 1.184 | 1.016 | 1.006 | 1.006 | 1.049 | 1.011 |
| L002013 | 1.406 | 1.224 | 1.080 | 1.056 | 1.100 | 1.098 |
| L000034 | 2.726 | 1.352 | 1.147 | 1.127 | 1.328 | 1.331 |
| L009803 | 1.329 | 1.171 | 1.064 | 1.052 | 1.108 | 1.087 |
| L014210 | 2.286 | 1.246 | 1.071 | 1.059 | 1.122 | 1.100 |
| L015331 | 1.810 | 1.266 | 1.064 | 1.060 | 1.152 | 1.083 |
| median | 1.371 | 1.197 | 1.064 | 1.054 | 1.104 | 1.085 |

"mafft on units, best setting" is L-INS-i with the SP-form matrix and `--ep 2`. The mafft_linsi and
abPOA columns are those agents' candidates, scored the same way.

**On all 149 regions,** before the polish and the guard, controls regressed badly: TR499276 went
from 1.00 to 3.95, L011224 from 1.03 to 2.46, L010113 from 1.00 to 1.63, TR499275 from 1.00 to
1.59. The polish brought these to 1.00, 1.51 (then the guard switched it to L-INS-i), 1.04 and 1.00.

**The unit-level gap-open cost, `--unit-go`** (all regions, go = 5 against the shipped 0):
- at hotspot VNTRs it improves the affine ratio (1.069 vs 1.092) and mean SV pieces per path
  (5.8 vs 7.0), at +0.004 cost/optimum;
- at STR hotspots it is worse: cost/optimum 1.088 vs 1.053, 13 vs 20 regions <= 1.1, median SV
  pieces 3.5 vs 0.5;
- overall 125 vs 134 regions are at cost/optimum <= 1.1.

So 0 ships, and `--unit-go 5` is the setting for hotspot VNTRs when SV pieces matter more.

**Array extent.** In three regions (L007376, TR743754, TR831034; flag `array_beyond_model`),
CHM13's own flank-cycle-flank segmentation runs past the model's array, because the flank carries
more repeat-like sequence than the model's search window saw.
- At L007376 and TR743754 the repeat genuinely continues.
- At TR831034 a 2-bp unit absorbs an AT-rich flank.
- Stricter DP scores (gap extend -3, and TRF-like 2/-7/-7) were tried on the three. They fixed
  TR831034's extent but gave a worse alignment (SP ratio 1.170 vs 1.131; MC 1.705), so the
  default scores stay and the case is only flagged.


## How to run it

```bash
# region mode: regions/<id>/hap32.fa -> candidates/unit_aware/<id>.{msa.fa,gfa,realign.json}
python3 tools/realign_units.py L012184                       # one region (id or regions/<id>)
python3 tools/realign_units.py all --jobs 2 --threads 2      # every packaged region; existing ok results kept
python3 tools/realign_units.py all --stratum hotspot_vntr --force

# MSA-only mode, for ANY input FASTA (e.g. a region's deduplicated full HPRC panel); writes only the MSA.
# The motif / period and the reference array come from --region-json, CHM13 from the input's CHM13#...
# record or, when the input has none, from the hap32.fa next to that region.json.
python3 tools/realign_units.py --input regions/L012184/hprc.fa.gz --region-json regions/L012184/region.json \
        --msa-out /path/L012184.all.msa.fa [--info-out /path/L012184.all.json] [--threads 2] [--fallback mafft_linsi]
```

```python
import sys; sys.path.insert(0, 'tools')
from realign_units import align_units_fasta, Params
info = align_units_fasta('seqs.fa', 'regions/L012184/region.json', 'out.msa.fa', threads=2,
                         params=Params(unit_go=0.0))          # every tunable is a Params field
info['status'], info['mode'], info['reason']                  # 'ok'; 'unit' or 'fallback'; why it fell back
```

The MSA has one row per input record, in input order, named by the first word of each header, with
`-` for gaps. Rows are upper case and each spells its input sequence exactly; this is checked before
the file is written. Projecting a full-panel MSA onto the hap32 rows means keeping those rows and
dropping all-gap columns: `msa_graph.match_rows(..., ignore_extra=True)`, or
`msa_graph.msa_to_gfa(msa, hap32.fa, out.gfa, ignore_extra=True)`.

Outputs in region mode:

| file | content |
|---|---|
| `candidates/unit_aware/<id>.msa.fa` | the MSA (hap32.fa names and order) |
| `candidates/unit_aware/<id>.gfa` | the graph from `msa_graph.msa_to_gfa` (vg construct -M + vg mod -u, one P line per row) |
| `candidates/unit_aware/<id>.realign.json` | everything below: mode and reason, repeat model, segmentation, symbols, unit MSA, polish, guard, timings, graph stats |
| `results/realign_runtime.tsv` | one row per (method, region), shared with realign.py (file lock) |

## The algorithm

### 1. Repeat model (CHM13 + region.json, independent of the panel)

- **Period.** The candidates are the annotated period, the annotated motif's length (when it is
  plain ACGT and not truncated), the periods of region.json's CHM13 tandem scan, and every divisor
  (>= 2) of those. The CHM13 core is scored by lag identity (the share of positions i with
  s[i] = s[i+p]), which is insensitive to sparse indels. The smallest candidate within 0.02 of the
  best is chosen.
  - The annotated motif is not used as is. In 23 of the 129 TR regions its length differs from the
    period, and above 203 bp it is truncated.
  - The truth-haplotype tandem scan in region.json is never read.
- **Unit.** It starts from the p-window of the core that agrees best with its neighbours at
  ±p and ±2p. It is then refined to the consensus of a local wrap-around alignment of the CHM13
  core (± max(2p, 50) bp) to it: the majority base per position, a deletion where most units
  delete, an insertion where most units insert. This is iterated to a fixed point (1-3 rounds).
- **Rotation.** The unit is rotated so that CHM13's array starts at unit position 0. The array is
  the local alignment's extent, and it defines the left and right CHM13 flanks.
- **Applicability.** The model is rejected, and the region falls back, when:
  - the region is non-TR, or has no period >= 2;
  - CHM13 has fewer than 3 units (each at least half a unit long);
  - unit identity is below 0.60;
  - the array covers under 50% of the annotated core.

### 2. Segmentation: flank -> (unit)* -> flank DP (C, `unit_align`)

Each distinct, N-masked sequence is aligned to a model: a free prefix, the last `flank_k` (300) bp of
the CHM13 left flank, the unit as a cycle traversed any number of times, the first 300 bp of the
right flank, and a free suffix. The prefix and suffix are free only when the flank is longer than
300 bp; otherwise the alignment is global, since every sequence starts and ends with the shared
anchors.

- **Scores.** Gap-affine, maximised: match +2, mismatch -4, gap open -6 (first base included),
  extend -1.
- **Cycle.** It is entered at unit position 0; a unit that starts later begins with deletions. It
  may be left after any position, so the last unit may be partial. It may also be bypassed, for a
  haplotype with no array. A deletion chain around the ring takes two passes per row.
- **Traceback.** It cuts the array at every wrap, so the units of all sequences are phased by the
  same unit. Each base keeps its unit position: a match at position j, or an insertion after j.
- **Cost.** n x (flanks + period) cells, 0.1-10 s per region for 34 sequences. For the full panel of
  L012184 (364 distinct sequences) it was 10 s.

Most of the time the unit representation is lossless. For the L012184 pair that is 454 edits apart,
the optimal unit-level pairwise alignment costs 449 against 447 for the base edit distance of the
two arrays.

### 3. Symbols

Every distinct unit string is its own symbol. Identical units share a symbol; near-identical ones
stay distinct, and their similarity enters through the cost matrix below, so nothing is merged. A
partial unit is simply a shorter string.

Only the optional mafft path (`--unit-aligner mafft`) is limited to 248 symbols (below). There, the
247 most frequent variants keep their own symbol and each rarer one is mapped to its nearest kept
variant by edit distance, for the unit-level alignment only. The expansion always uses each unit's
own sequence.

### 4. The unit-level MSA: an exact sum-of-pairs progressive aligner (C)

The objective is the base-level edit cost, written on units:
- substituting unit a for unit b costs ed(a, b);
- gapping a unit costs its length;
- gaps are linear, so the sum of pairs of a profile merge is exact.

The steps:
- **Guide tree.** UPGMA on the exact unit-level pair costs. Above 150 strings, or when that would
  exceed 3e9 cells, it uses hashed symbol-3-mer L1 distances instead.
- **Progressive merges.** Each merge is an optimal profile-profile DP under that cost
  (`profile_align`).
- **Refinement.** Every guide-tree split, then every single row, is realigned against the rest. The
  result is kept when the cost of the pairs across the split falls; pairs within each side are
  unchanged. Up to 3 rounds; the tree splits only apply up to 120 strings.
- **Gap open.** `--unit-go G` adds a gap-open cost per pair and run of gapped units, as a
  Gotoh-on-profiles approximation in the DP and exactly in the acceptance test. It trades
  cost/optimum for fewer SV-sized pieces (measured below).

**Why not `mafft --anysymbol`, and why not mafft at all.** With `--anysymbol`, mafft scores every
character that is not an amino-acid letter as unknown; it runs `replaceu` before aligning and
`restoreu` after. So it cannot tell unit variants apart.
- **The text route.** mafft's way to align arbitrary symbols is text mode with a user score matrix
  (`--textmatrix FILE`, lines `0xAA 0xBB score`). It accepts 248 byte values: everything except
  newline, CR, space, `-`, `<`, `=` and `>` (checked with `maffttext2hex`). On a toy case,
  `--anysymbol` and plain `--text` misplaced a gap that `--textmatrix` got right.
- **The locale trap.** Text mode silently outputs only the first header line when `LC_CTYPE` is
  UTF-8. Python sets `LC_CTYPE=C.UTF-8` for its subprocesses (PEP 538), so this happens whenever
  mafft is run from Python. The mafft path therefore runs with `LC_ALL=C`.
- **Measured.** On 10 pilot hotspots, mafft L-INS-i / E-INS-i / G-INS-i on the symbols, with two
  matrix forms and a sweep of `--op` / `--ep`, reached a mean cost/optimum of 1.18 at best. The
  exact-SP aligner reached 1.064 before refinement, polish and guard (table below).
- **Why.** mafft normalises a user matrix (it subtracts the average), so near-identical unit
  variants end up scoring close to zero and gaps look cheap. Under edit cost, gapping a unit costs a
  whole period.
- **Kept as an option.** mafft stays available as `--unit-aligner mafft`, with `--unit-mafft`,
  `--unit-op`, `--unit-ep` and `--matrix`.

### 5. Expansion to bases

**Unit columns.** A unit column is expanded through the unit: every unit variant is already aligned
to the unit by the DP, so it occupies the match columns 0..p-1. Insertions between the same two unit
positions are star-aligned: stacked when they have the same length, otherwise aligned around the
most frequent one. Columns that are gaps in every unit of that column are dropped. Each distinct
unit string uses one placement (the most common), so identical units in one column give identical
rows.

**Flanks.** The left and right flank pieces are aligned with realign.py's pipeline (dedup, N
masking, a checked MSA):
- mafft L-INS-i by default;
- FFT-NS-i above 100 distinct pieces;
- FFT-NS-2 above 20 kb.

Decomposition on two pilots:

| region | array: sum of pairwise optima | unit-level SP | after expansion | flanks |
|---|---|---|---|---|
| L012184 | 845,964 | 893,817 (+5.7%) | 897,351 (+0.4%) | 0 excess |
| L000034 | 306,572 | 352,416 (+15%) | 356,121 (+1.1%) | 0 excess |

So the expansion is nearly lossless, and what remains is the unit-level MSA against the pairwise
optima. Part of that is inherent, since one MSA cannot realise every pairwise optimum at once.

### 6. Base-level polish

In a degenerate array, one SNP next to a unit boundary can move the boundary. For example,
`GTTCTCC|TTTCATC` against `GTTCTC|ATTTCATC` turns one mismatch into an insertion plus a deletion. At
low unit identity such shifts can cascade.

The polish fixes this. Each row is realigned against all the others (`row_to_profile_banded`):
- the cost is the same exact SP cost at base level, plus a gap-open term of 0.5 per pair so that ties
  keep gaps contiguous;
- the DP is restricted to a band of max(100, 3 x period) bases around the row's current path, capped
  at 3e7 cells;
- a change is kept only when the row's pairs get cheaper;
- 2 rounds.

At TR499276, 5 near-identical haplotypes whose MC graph is perfect, the unit MSA cost 6.1x the
pairwise optima before the polish and 1.00x after.

### 7. Adequacy guard

After the polish, the MSA's sum-of-pairs unit cost is divided by the sum of the pairs' optimal edit
distances (bit-parallel Myers/Hyyrö in C). All pairs of distinct sequences are used, or a seeded
sample of at most 3,000 pairs within a 3e10-cell budget. This truth-free ratio tracks evaluate.py's
graph cost/optimum closely.

When it exceeds `--guard-ratio` (1.2), the fallback MSA is computed too, and the MSA with the lower
cost over the same pairs is kept. The info JSON records both values and the decision.

### 8. N runs, fallback

**N runs.** Every maximal non-ACGT run is cut out before segmentation and put back afterwards as
columns of its own. This is realign.py's `mask_runs` / `reinsert_runs`, so an N-gap haplotype adds
one N node and no misalignment.

**Fallback.** realign.py's `mafft_linsi` pipeline (`--fallback` picks another realign.py method) is
used for:
- non-TR regions;
- a rejected model;
- the guard;
- any exception in the unit path, recorded as `reason`.

## Parameters (`Params`, CLI flags)

| parameter | default | flag |
|---|---|---|
| DP scores match, mismatch, open, extend | 2, 4, 6, 1 | `--dp-scores 2,4,6,1` |
| flank anchor length | 300 bp | `--flank-k` |
| unit-level aligner | `prog` | `--unit-aligner prog\|mafft` |
| refinement | tree splits + single rows, 3 rounds | `--refine tree\|loo\|none`, `--refine-rounds` |
| unit-level gap open | 0 (5 = fewer SV pieces at VNTR hotspots, worse at STRs) | `--unit-go` |
| polish rounds / gap open / band | 2 / 0.5 / max(100, 3p) | `--polish-rounds` |
| guard ratio | 1.2 | `--guard-ratio` (0 = off) |
| model acceptance | >= 3 units, identity >= 0.60, core cover >= 0.50 | `--min-units`, `--min-identity`, `--min-core-cov` |
| fallback | `mafft_linsi` | `--fallback` |

## The info JSON (`<id>.realign.json`, `--info-out`)

- `status`, `mode` (`unit` / `fallback`), `reason`;
- `model`: period, motif (rotated), rotation, lag identities of the candidates, CHM13 array
  coordinates, unit count, identity, core cover;
- `array_bp` (min / median / max over sequences), `no_array`, `masked_runs` / `masked_bp`;
- `units`:
  - `unit_tokens`, `unit_variants` (distinct unit strings), `symbols`, `variants_hashed`
    (mafft path only);
  - `unit_strings_distinct`, `units_per_seq`, `unit_columns`;
  - `unit_msa`: guide tree, SP before and after refinement, accepted refinements, seconds;
- `phase_check`: the share of unit columns (and of units) whose units all start at the same unit
  position;
- `polish` (SP before and after, rows accepted), `guard` (ratio, fallback ratio, decision);
- `left_flank` / `right_flank` (method, seconds), `columns`, `gap_frac`;
- `graph` (msa_graph stats), timings (`segment_s`, `unit_msa_s`, `flank_s`, `align_s`, `graph_s`,
  `total_s`).

## Limitations

- **One unit per region.** A compound array made of two motifs is modelled by the chosen period,
  and the other motif's stretch becomes degenerate units. The polish and the guard contain the
  damage, but an arrangement of several cycles (flank, cycle A, spacer, cycle B, flank) would model
  it properly.
- **Unit order.** Arrays with differently ordered unit types cannot be stacked acyclically. The
  exact-SP MSA chooses the cheapest compromise, which is up to 15% above the pairwise optima at
  L000034.
- **The objective is linear.** Scattered unit-level gaps give more SV-sized pieces against CHM13
  than a gap-affine objective would; `--unit-go` trades cost/optimum for fewer pieces.
- **Full panel.** Exact guide-tree costs and tree-split refinement switch off above 150 and 120
  strings; leave-one-out refinement and the polish remain, with a 300 s cap on each.
