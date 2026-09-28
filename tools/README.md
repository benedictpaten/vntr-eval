# tools

Python 3 (standard library only) plus a little C. Data and tool locations come from
[config.py](config.py): defaults point at `~/PycharmProjects/vg-call-eval` (read only) and at this
repository; every location has an environment override (`VNTR_EVAL_DIR`, `VNTR_REGIONS`,
`VNTR_CANDIDATES`, `VNTR_RESULTS`, `VNTR_CENSUS`, `VNTR_WORK`, `VNTR_VG`, `VNTR_MAFFT`, ...).
`python3 tools/config.py` prints what is in use. Caches and intermediates go to `work/` (untracked).

| tool | what it does | needs the big data |
|---|---|---|
| [msa_graph.py](#msa_graphpy) | aligned rows -> GFA with one path per `hap32.fa` sequence (the builder every realigner uses) | no |
| [viewer.py](#viewerpy) | one self-contained HTML page per region: MSA, baseline and candidate graphs, calls, reads | no (reads and called haplotypes only with it) |
| [region.py](region.md) | everything about one CHM13 interval of the hap32 graph: anchored subgraph, haplotypes, truth haplotypes, calls, reads | yes |
| [package_regions.py](../regions/README.md) | builds the region packages in `regions/` | yes |
| [evaluate.py](METRICS.md) | Stage 0/1 metrics of one graph of one region -> `results/<method>/<id>.json` | no (Stage 1 reads: yes) |
| [realign.py](#realignpy) | re-align a region's `hap32.fa` (or any FASTA) with mafft FFT-NS-2 / L-INS-i / E-INS-i / G-INS-i (and plugin methods: `realign_poa.py`) and build `candidates/<method>/<id>.gfa` | no |
| [panel.py](#panelpy-the-full-panel-arm) | the full-panel arm: union of all haplotypes, projection of a full-panel MSA onto hap32, the full MC graph over each span, panel statistics and checks | `union`/`project`/`panel-graph`: no; the rest reads cached full-graph queries |
| [poa_panel.py](POA.md#full-panel-arm-poa_abpoa__all-poa_spoa__all) | the POA methods on the full panel: union -> abPOA/spoa MSA -> full-panel graph (`work/panel/<m>/`) and projected hap32 graph (`candidates/<m>__all/`), runtimes in `results/realign_runtime.all.poa.tsv` | no (reads the union files `panel.py` wrote) |
| [mafft_panel.py](#mafft_panelpy-the-mafft-methods-on-the-full-panel) | the mafft methods (FFT-NS-2, FFT-NS-i, L-/E-/G-INS-i) on the full panel: union -> MSA -> full-panel graph (`work/panel/<m>/`) and projected hap32 graph (`candidates/<m>__all/`), runtimes and every timeout or skip in `results/realign_runtime.all.mafft.tsv` | no (reads the union files `panel.py` wrote) |
| [stage01.py](#stage01py-the-stage-0-1-comparison) | the Stage 0-1 comparison of every realigner against mc, both arms and the full-panel graphs, paired per stratum -> `results/stage01_tables.md`, `results/stage01_pairs.tsv` (narrative: `results/stage01.md`) | no (reads `results/`) |
| [evaluate_all.py](METRICS.md) / [summarise.py](METRICS.md) | run evaluate.py over every region of one or more methods (`--panel` for full-panel graphs); per-stratum comparison tables, paired against a baseline | no (Stage 1 reads: yes) |
| [realign_units.py](UNIT_ALIGN.md) | the repeat-unit-aware realigner (`unit_aware`; C core `unit_dp.c`): cut arrays into units, align unit strings with an exact edit-cost progressive aligner, expand, polish; falls back to mafft L-INS-i where the motif does not fit | no |
| [units_panel.py](UNIT_ALIGN.md) | `unit_aware` on the full panel: union -> MSA -> full-panel graph (`work/panel/unit_aware/`) and projected hap32 graph (`candidates/unit_aware__all/`), runtimes in `results/realign_runtime.all.units.tsv` | no (reads the union files `panel.py` wrote) |
| [realign_poa_mc.py](POA.md#poa_abpoa_mc-abpoa-run-the-way-minigraph-cactus-runs-it) / [poa_panel_mc.py](POA.md#poa_abpoa_mc-abpoa-run-the-way-minigraph-cactus-runs-it) | abPOA with Minigraph-Cactus's BAR settings and 10 kb rule (`poa_abpoa_mc`), hap32 arm and full-panel arm; runtimes in `results/mcpoa_runtime{,.all}.tsv` | no |
| [remap_local.py](../results/stage2_remap.md) | Stage 2: rebuild a region's reads from GAF-Base, map them with giraffe to each local graph (MC, candidates, HG002's own graph), compare placement | yes |
| [call_local.py](../results/stage3_gate.md) | Stage 3 caller: the graph's span inside the genome-wide hap32 graph +- 200 kb (same panel, CHM13/GRCh38 as reference samples), reads re-mapped with giraffe, `vg call` with the production flags (pinned vg 2a6a228a5), shifted to CHM13 -> `work/stage3/calls/<graph>/<id>.vcf.gz`; `gate` compares local MC with production (diagnostic arms decompose every difference); null graphs `mc_relabel` (MC with candidate-style node IDs) and `mc_unchop` (MC unchopped) -> `work/stage3/diag/<null>@hybrid200k/`; replicate arms `hybrid50k` (50 kb flank) and `hybrid200kids` (span node IDs re-laid) for any graph (`run --arms ... --arm-graphs ...`) -> `work/stage3/diag/<graph>@<arm>/` | yes |
| [replicate_noise.py](../results/stage3_noise.md) | Stage 3 noise per graph, without the truth: the pair edit distance between each graph's Stage 3 call and its replicate calls -> `results/stage3_noise.{tsv,md}`, read by `score_haplotypes.py summarise --noise` | no (reads the calls) |
| [null_decomp.py](../results/stage3_pilot.md) | Stage 3: where a null graph's call difference comes from (alignments vs caller, depth term); `rescore-*` recomputes its EDs from the kept VCFs | yes (vg runs) |
| [score_haplotypes.py](#score_haplotypespy-the-stage-3-scorer) | Stage 3 scorer: the haplotypes a VCF writes over a region's span against HG002's truth (edit distance, paired per stratum), plus truvari bench/refine/phab; batch, summary and the scorer's own validation | truvari part and `validate`: yes |

Tests (offline, seconds each): every
`tools/test_*.py`, for example `python3 tools/test_msa_graph.py`, `test_evaluate.py`, `test_panel.py`,
`test_realign.py`, `test_realign_poa.py`, `test_realign_poa_mc.py`, `test_realign_units.py`,
`test_poa_panel.py`, `test_poa_panel_mc.py`, `test_mafft_panel.py`, `test_units_panel.py`, `test_region.py`,
`test_score_haplotypes.py`.

## msa_graph.py

```bash
python3 tools/msa_graph.py candidates/<method>/L012184.msa.fa regions/L012184/hap32.fa candidates/<method>/L012184.gfa
python3 tools/msa_graph.py MSA.fa hap32.fa OUT.gfa --merge-blocks [--block-max 50]   # remove 1-bp mesh
```

```python
import sys; sys.path.insert(0, 'tools')
from msa_graph import msa_to_gfa, read_msa, write_msa, MsaGraphError
stats = msa_to_gfa(msa_fa, hap32_fa, out_gfa)          # raises MsaGraphError on any mismatch
```

- **Input**: FASTA with `-` or `.` gaps (any case, wrapped or not) or PIR; abPOA `-r 1`/`-r 2` and
  spoa `-r 1` output work as is (abPOA's `Consensus_sequence` row is dropped). Rows are matched to
  `hap32.fa` by the first word of the header. Every row's ungapped sequence must equal its `hap32.fa`
  sequence (case-insensitive); a mismatch, a duplicate, a missing sequence (unless `--allow-missing`)
  or an unknown row (unless `--ignore-extra`) is an error.
- **Output**: GFA 1.0, upper-case S lines, forward L lines, one P line per row in `hap32.fa` order
  and with its name, re-read and checked against `hap32.fa` after writing. Node ids run 1..n in
  topological order.
- **Induction**: column induction (rows with the same base in a column share a node), then
  compaction (unchop: merge u->v when that is the only edge out of u and into v and no path ends
  or starts between them). With vg installed, `--engine auto` (default) builds the graph with
  `vg construct -M` + `vg mod -u` on rows renamed `s0..sN` (vg rewrites `#` names as PanSN W lines)
  and checks it equals the native build; `--engine native` needs no vg.
- **`--merge-blocks`**: cut the MSA into blocks, maximal runs of columns with the same set of
  non-gap rows; a block of at most `--block-max` columns (default 50, below SV size) gets one node
  per distinct row string instead of one per base. The stats report node counts both ways.
- **N**: an N is a base like any other (never merged with A/C/G/T), so N runs become their own
  nodes and paths keep spelling their sequence; `N_bp_in_nodes` counts them.
- **Stats** (printed as JSON, `--json FILE`): rows, columns, gap fraction, nodes, edges, bp in
  nodes, nodes per CHM13 kb, mean node length, share of 1-bp nodes, N bp, the same for
  `columns_mode` and `blocks_mode`, the engine, and whether vg and the native build agree.

On L012184 (chr4, 32 spanning sequences), `mafft --auto` gives 17,261 columns -> 3,239 nodes
(67% 1 bp), 1,855 with `--merge-blocks`; the Minigraph-Cactus subgraph has 9,292. On L005990
(chr17): 14,289 columns -> 838 nodes, 619 merged, against 3,484. vg and the native build agreed
on both, and on 25 random MSAs in the test.

## viewer.py

```bash
python3 tools/viewer.py L012184                                   # regions/L012184 (+ big data if present)
python3 tools/viewer.py L012184 --candidate candidates/mafft_linsi/L012184.gfa \
                                --candidate candidates/abpoa/L012184.gfa        # one panel per graph
python3 tools/viewer.py regions/L012184 --no-data                 # what a collaborator without the data sees
python3 tools/viewer.py chr4:191498564-191501883                  # coordinates, via region.py (big data)
```

The page (default `work/viewer/<id>.html`, self-contained, light and dark, ~1-4 MB) has:

- **Alignment**: an MSA of CHM13, the panel sequences, the HG002 truth haplotypes (evaluation only)
  and vg's and PanGenie's genotypes applied to CHM13. `alignment` switches between the mafft MSA of
  all rows and each candidate's own MSA (`<candidate stem>.msa.fa`, with the truth and called rows
  added by `mafft --add`, which keeps its columns). `colour` shows differences against any row, or
  the nodes of any graph: node identity, node sharing, or the graph/MSA split (same base in a
  column on different nodes, i.e. unmerged or misaligned copies).
- **Graphs**: a metrics table with one column per graph (the viewer's own numbers, plus
  `tools/evaluate.py`'s from `results/mc/<id>.json` and `results/<method>/<id>.json` when they
  exist); a Bandage image per graph coloured by path sharing; one node track per graph that
  scrolls with the alignment.
- **Calls** (vg, PanGenie, truth records with truvari status) and **Reads** (MAPQ, reads per
  walk, a pileup on the walk each read follows; only with the big data).

A candidate is a GFA with one P (or W) line per `hap32.fa` sequence, same names, same sequences;
the viewer warns about paths that do not match. `--no-bandage`, `--max-reads`, `--max-cols` and
`--max-page-mb` keep a page small; `--mafft linsi|einsi|ginsi|fftns` changes the all-rows MSA.

## realign.py

```bash
python3 tools/realign.py --list                                         # the method table (+ plugins)
python3 tools/realign.py mafft_linsi L012184                            # -> candidates/mafft_linsi/L012184.{msa.fa,gfa,realign.json}
python3 tools/realign.py mafft_fftns2,mafft_linsi all --jobs 2 --threads 2 [--stratum hotspot_vntr] [--timeout 900] [--mem-mb 10000]
python3 tools/realign.py mafft_einsi all --retry memout --threads 1 --mem-mb 16000   # memouts again, one thread, more memory
python3 tools/realign.py mafft_linsi --input ANY.fa[.gz] --msa-out OUT.msa.fa [--json info.json]   # MSA only (full-panel arm)
```

```python
import sys; sys.path.insert(0, 'tools')
from realign import align_fasta, realign_region
info = align_fasta('mafft_linsi', 'seqs.fa', 'out.msa.fa', threads=2, timeout=900)   # info['status']: ok|timeout|memout|error
```

- **Methods**: `mafft_fftns2` (`--retree 2 --maxiterate 0`), `mafft_linsi` (`--localpair --maxiterate 1000`),
  `mafft_einsi` (`--genafpair --maxiterate 1000 --ep 0`), `mafft_ginsi` (`--globalpair --maxiterate 1000`),
  `mafft_fftnsi` (`--retree 2 --maxiterate 1000`, for the full-panel arm); every call adds `--nuc --thread T
  --threadit 0` (single-threaded refinement, for reproducible output). Any `tools/realign_*.py` defining
  `METHODS = {name: {'align': func, 'params': {...}, 'description': ..., 'tool': ...}}` adds methods
  (`realign_poa.py`: `poa_abpoa`, `poa_spoa`); `func(in_fa, out_fa, threads=, workdir=, timeout=, mem_mb=, **params)`
  returns `{'status': ..., 'command': [...]}`, and `run_proc()` gives it the timeout and memory cap.
- **Shared pipeline** (every method): identical sequences aligned once; every run of non-ACGT characters
  (N) cut out before alignment and put back as columns of its own right after the preceding base, so the
  MSA still spells the sequence and `msa_graph.py` makes the run one N node (paths must spell their Ns;
  see evaluate.py); the aligner sees the distinct masked sequences as `s0..sK`, under `--timeout` (default
  900 s) and `--mem-mb` (default 10000, summed RSS of its process group, sampled each second); each output
  row is checked against its input.
- **Outputs** (region mode): `<id>.msa.fa` (hap32 names and order, `-` gaps, one line per row), `<id>.gfa`
  (`msa_graph.msa_to_gfa`, column induction + unchop, vg construct -M cross-check), `<id>.realign.json`
  (status, command, mafft version, N masking, graph stats, timings; also written for a timeout or memout,
  which then has no MSA or GFA), and one row per (method, region) in `results/realign_runtime.tsv`
  (latest run wins). Existing results are skipped unless `--force`, `--retry-failed` or `--retry STATUS`.
  `--fragments` also writes `<id>.fragments.msa.fa` (`mafft --addfragments --keeplength`; not in the graph).

## panel.py: the full-panel arm

`mc.gfa` is not an alignment of the 34 hap32 sequences: Minigraph-Cactus aligned the whole HPRC
v2.1 panel and hap32 keeps 32 sampled rows of it plus CHM13 and GRCh38. A like-for-like
replacement for MC therefore aligns the full panel and is projected onto hap32. For a method *m*
the comparison is m(hap32) (`candidates/<m>/`) against m(all) projected onto hap32
(`candidates/<m>__all/`) against `mc`; the full-panel graphs themselves are compared with MC's
full-panel graph under `evaluate.py --panel`.

The full panel of a region is `regions/<id>/hprc.fa.gz`: CHM13 and GRCh38 under their hap32 names,
then every sample haplotype of the eval graph (`hprc-v2.1-mc-chm13-eval.gref.gbz`; 456 haplotypes
of 228 HPRC samples, no HG002, HG003 or HG004) that passes through both anchor nodes, as anonymous
`hprc#k` records grouped by identical sequence. `gref_CHM13` is not a panel member.

```bash
# for a method: distinct sequences of the region, align them, then project and/or build the panel graph
python3 tools/panel.py union regions/L012184 work/panel/union/L012184.fa     # + work/panel/union/L012184.map.tsv
<aligner> work/panel/union/L012184.fa > FULL.msa.fa                           # rows named u0001..
python3 tools/panel.py project FULL.msa.fa work/panel/union/L012184.map.tsv regions/L012184/hap32.fa \
                               candidates/<m>__all/L012184.msa.fa
python3 tools/msa_graph.py candidates/<m>__all/L012184.msa.fa regions/L012184/hap32.fa candidates/<m>__all/L012184.gfa
python3 tools/panel.py panel-graph FULL.msa.fa work/panel/union/L012184.map.tsv regions/L012184 \
                               work/panel/<m>/L012184.gfa [--merge-blocks]
gzip -c FULL.msa.fa > work/panel/<m>/L012184.msa.fa.gz

# judge them
python3 tools/evaluate.py regions/L012184 candidates/<m>__all/L012184.gfa --reads                  # hap32 arm, as any candidate
python3 tools/evaluate.py regions/L012184 work/panel/<m>/L012184.gfa --panel regions/L012184/hprc.fa.gz   # full panel
python3 tools/evaluate_all.py mc <m> --panel --jobs 1 --threads 3                                  # all regions
python3 tools/summarise.py --results results/panel --baseline mc --methods mc,<m>                   # full-panel tables

# set-up and checks (already run; outputs listed below)
python3 tools/panel.py union-all                     # every region -> work/panel/union/<id>.fa + .map.tsv
python3 tools/panel.py import-queries <package work dir>   # or: python3 tools/panel.py fetch (gbz-base, ~60 s each)
python3 tools/panel.py mc-graph                      # -> work/panel/mc/<id>.gfa (+ summary.tsv)
python3 tools/panel.py stats                         # -> results/panel_stats.tsv
python3 tools/panel.py premise                       # -> results/panel_premise.tsv
python3 tools/panel.py gref-check                    # -> work/panel/gref_check.json (one vg load, ~12 GB, ~6 min)
python3 tools/panel.py check-hprc
python3 tools/evaluate_all.py mc --panel --jobs 1 --threads 3   # -> results/panel/mc/<id>.json (MC's full-panel baseline)
```

**union** writes the distinct sequences of `hap32.fa` and `hprc.fa.gz` together (every hap32
sequence is included, so a sampled recombinant whose junction falls inside the span is covered),
CHM13's sequence first as `u0001`, then by decreasing weight. Header `>u0001 len=L weight=W
n_hap32=H`. The map TSV has one row per distinct sequence: `id`, `length`, `weight` (full-panel
records with this sequence, CHM13 and GRCh38 counted once each; 0 for a hap32-only sequence),
`n_hap32`, `n_panel_samples` (weight without the reference records), `n_bp_N`, and `members`
(hap32 names, then full-panel names). Names of HG002, HG003, HG004 (or their NA aliases) or
gref_CHM13 are an error.

**project** keeps one row per `hap32.fa` name (the row of its distinct sequence; MSA rows may be
named by distinct id or by any member name, abPOA's consensus row is ignored), drops the columns
that are all gap among those rows, and asserts that every ungapped row equals its `hap32.fa`
sequence. The output is an MSA of `hap32.fa` for `msa_graph.py`, like any hap32-arm MSA.

**panel-graph** induces the graph from the distinct rows with `msa_graph.py` (same column
induction, compaction and `--merge-blocks`), then writes one P line per `hprc.fa.gz` record,
following its distinct sequence's path, and checks every spelling. The map must be the one
`union` wrote for the region (it is recomputed and compared).

**mc-graph** writes `work/panel/mc/<id>.gfa`, the full Minigraph-Cactus graph between the region's
anchors: the nodes and edges of the eval graph that the spanning walks use, with one P line per
`hprc.fa.gz` record spelling it exactly (asserted). The eval graph keeps the hap32 graph's node
ids (vg haplotypes does not renumber), so the span is located by the same anchor nodes; the walks
come from the cached interval query package_regions.py ran (`work/panel/query/<id>.gfa.gz`), cut at
the anchors exactly as it did, minus the reference walks (the CHM13 walk twice: CHM13 and its
gref_CHM13 copy; each GRCh38 path of `hap32.fa` once). The reconstruction must reproduce
`hprc.fa.gz`'s multiset of sample sequences, or the region fails. Within a group of identical
sequences, records take the walks in order of decreasing walk weight (two walks spelling one
sequence occur once in 149 regions). 548 MB for 149 regions.

**stats** (`results/panel_stats.tsv`, one row per region):

| column | meaning |
|---|---|
| `span_bp`, `chm13_bp` | anchor-to-anchor CHM13 length |
| `n_hap32`, `n_hap32_sampled`, `n_hap32_distinct` | spanning hap32 paths, those that are sampled (not CHM13/GRCh38), distinct sequences |
| `hprc_spanning`, `hprc_not_spanning`, `hprc_beyond_456` | sample haplotypes running anchor to anchor; 456 minus that; segments above 456 (an assembly carrying two contigs through the span) |
| `walks_left_anchor_only`, `walks_right_anchor_only`, `walks_anchor_other` | query walks that touch an anchor without spanning (walks, not haplotypes: a haplotype broken inside the span gives one of each) |
| `reference_removed` | the reference paths subtracted: CHM13, gref_CHM13 and the GRCh38 path(s) |
| `chm13_walk_weight_before` | haplotypes on the CHM13 walk before subtraction (2 = no sample carries it) |
| `gref_fragments_in_span`, `gref_fragments_with_anchor`, `gref_copy_equals_chm13`, `grch38_spanning_in_graph` | from `gref-check` (below) |
| `panel_distinct_seqs`, `panel_distinct_walks`, `union_distinct_seqs` | distinct sample sequences and walks; distinct sequences of the union |
| `hap32_sampled_verbatim_seq`, `hap32_sampled_verbatim_walk`, `hap32_not_verbatim` | sampled hap32 paths whose sequence / walk is that of some full-panel haplotype |
| `hap32_not_verbatim_kinds`, `hap32_not_verbatim_why` | per path: `junction` (two pieces of full-panel walks: a recombination inside the span), `mosaic_<k>`, `copy_of_CHM13` / `copy_of_GRCh38` (the reference walk, which no full-panel haplotype takes), `mixed(...)` (needs a walk that does not span, or a reference) |
| `distinct_bp_max`, `distinct_bp_min`, `distinct_bp_total` | lengths of the union's distinct sequences |
| `distinct_N_bp`, `distinct_seqs_with_N`, `panel_haplotypes_with_N` | N content |
| `hap32_sampled_top_seq_share` | share of the spanning haplotypes carried by the most common sequence |

**premise** tests that `mc.gfa` is the full graph restricted to the hap32 walks: every `mc.gfa` node
id has the same sequence in the full graph, every edge exists there, each hap32 path's walk is a
full-panel walk (or a reference walk), and for every pair of hap32 paths the alignment implied by
shared node visits (as `evaluate.py` pairs them) and its unit cost are the same in `mc.gfa` and
between the corresponding full-panel walks. On all 149 regions: 345,486 nodes and 406,612 edges,
none missing or different; all 79,212 pairs of hap32 paths with a full-panel counterpart imply the
same alignment. The 114 sampled paths without one (2.4%, in 32 regions: 88 junctions, 24 reference
copies, 2 mixed) are walks of the full graph too, so the implied alignment is the full graph's
for them as well. Differences between an `mc` hap32 score and an `mc` full-panel score therefore
come from which haplotype pairs are scored, not from the graph.

**gref-check** streams `vg paths -x <eval GBZ> -R -A` (every reference-sense path) once and checks,
for every region, that the gref_CHM13 base copy takes exactly the CHM13 walk (and the whole
contig's walks are byte-identical: 25/25), that no gref fragment contains an anchor or any CHM13
node of the span (149/149 clean; 3,372 fragments enter spans off the reference, in 132 regions),
and that the GRCh38 paths running anchor to anchor are the hap32.fa ones with the same walk
(149/149). That is exactly the subtraction package_regions.py made by weight.

**evaluate.py --panel FASTA** judges a graph whose paths are a panel (see its docstring and
[METRICS.md](METRICS.md#the-full-panel-arm---panel)): path-spells-sequence is asserted, hap32,
truth and read checks are skipped, alignment metrics use seeded samples of `--max-pairs` (300)
uniform haplotype pairs and CHM13 x 300 haplotypes. Pairs are drawn by record index, so two graphs
of one `hprc.fa.gz` are scored on the same pairs (a paired comparison); across seeds the
all-pairs ratio of sums is noisy at hotspots (L012184: 1.43, 1.47, 1.72 for seeds 1-3; CHM13 pairs
1.28-1.30), and `all_cost_over_opt_ci95` (bootstrap over the sampled pairs) understates that,
because a few very long haplotypes dominate the sums. Use `--max-pairs 1000` for tighter numbers.

## mafft_panel.py: the mafft methods on the full panel

```bash
python3 tools/mafft_panel.py run --jobs 2 --timeout 1800 --mem-mb 12000 --budget-mb 18000   # every region x method
python3 tools/mafft_panel.py one mafft_linsi L015415 [--threads 3 --timeout 1800]           # one job
python3 tools/mafft_panel.py status                                                          # counts per method, status, stratum
python3 tools/mafft_panel.py predict --methods mafft_linsi                                   # the cost and memory model
python3 tools/evaluate.py regions/L015415 candidates/mafft_linsi__all/L015415.gfa --reads            # projected graph, hap32 arm
python3 tools/evaluate.py regions/L015415 work/panel/mafft_linsi/L015415.gfa --panel regions/L015415/hprc.fa.gz   # full-panel graph
```

Per region and method: `panel.py union` -> `realign.align_fasta` (the same method definition, dedup,
N masking and checks as the hap32 arm) -> `work/panel/<m>/<id>.msa.fa.gz` (full MSA, rows `u0001..`)
and `work/panel/<m>/<id>.gfa` (one P line per `hprc.fa.gz` record, spelling checked) ->
`panel.py project` -> `candidates/<m>__all/<id>.msa.fa` and `.gfa` (every path re-read and
asserted to spell `hap32.fa`) plus `<id>.realign.json` (also for timeouts and skips).

The all-pairs modes cannot run on hundreds of multi-kb sequences everywhere, so jobs start
smallest first (rank = predicted seconds x a per-method weight) under a 30 min cap, and a method
that has timed out on two smaller regions (by its cost measure: pair cells for L-/E-/G-INS-i,
N x L1 x L2 for FFT-NS) is not attempted on larger ones; the batch also had a wall budget. Every
region x method has a row in `results/realign_runtime.all.mafft.tsv` (first columns `region_id,
method, panel, n_distinct, seconds, status`; status `ok`, `timeout`, `memout` or `skipped`, the
reason in `note`). Compare methods only on regions where both graphs exist. The model, the
thresholds and the options are in the module docstring.

## stage01.py: the Stage 0-1 comparison

```bash
python3 tools/evaluate_all.py mc $(ls candidates) --reads --jobs 3 --threads 1        # results/<method>/ (both arms)
python3 tools/evaluate_all.py mc <methods> --panel --panel-subdir full --jobs 1       # results/full/<method>/
python3 tools/stage01.py [--md results/stage01_tables.md] [--tsv results/stage01_pairs.tsv] [--per-region FILE]
```

It reads `results/<method>/` (mc, every `<m>` and `<m>__all`) and `results/full/<method>/`, and writes
tables for the questions in [../results/stage01.md](../results/stage01.md):
- coverage and why regions are missing (from the runtime tables);
- A: hotspot VNTRs against matched controls, with the share of mc's hotspot-control gap each method
  closes on its own hotspots;
- B: controls against mc (regressions > 0.01 / 0.05, affine);
- C: all methods of an arm on their common regions, with best-per-region counts;
- D: m(hap32) against m__all against mc, then full MC against mc.gfa, then realigned full-panel graphs
  against full MC;
- E: truth spellability and pooled truvari F1 (raw / refined / phab);
- F: pooled MAPQ<5 read classes.

Every comparison is paired over the regions where all the graphs it names exist, with n.
`stage01_pairs.tsv` holds the same numbers in long form (one row per table, method, stratum and metric).

## score_haplotypes.py: the Stage 3 scorer

```bash
python3 tools/score_haplotypes.py score regions/L014297 calls.vcf.gz --label mc --out L014297.json [--haplotypes called.fa]
python3 tools/score_haplotypes.py batch --label mc --vcf 'work/stage3/call/mc/{id}.vcf.gz' --regions pilot --jobs 4
python3 tools/score_haplotypes.py batch --label genomewide --vcf genomewide --labels-dir genomewide --regions all   # production
python3 tools/score_haplotypes.py summarise --labels mc,unit_aware__all,mafft_linsi,unit_aware,truth,genomewide \
                                            --baseline mc --regions pilot [--noise results/stage3_noise.tsv]
                                                                                   # -> results/stage3_summary.{tsv,md}
python3 tools/score_haplotypes.py validate --jobs 4                                # -> results/stage3_scorer_validation.tsv
python3 tools/test_score_haplotypes.py
```

- **Input**: a region package and a VCF in CHM13 coordinates (the region's contig name, e.g.
  `chr6`; sample `HG002`, or a single-sample VCF under any name). bgzipped and tabix-indexed, or
  plain text. `--vcf` in batch mode is a pattern with `{id}`, `{contig}`, `{stratum}`; `genomewide`
  means `vg-call-eval/work/wgs-mm095/{contig}/{contig}.vcf.gz`.
- **The rule** is stated in full in the module docstring. In short: every record overlapping the
  span is applied to CHM13 once per GT slot; unknown phase between blocks (PS sets, unphased hets)
  is chosen to minimise ED, exactly up to 2^12 assignments, greedily beyond (flagged); records
  crossing the span edge are clipped (flagged); missing/FILTERed/symbolic alleles count as
  reference (flagged); overlapping records in one slot are resolved in (POS, longest first)
  order by what each allele actually changes, except that a vg record is never applied before
  a record it is nested in (nesting from vg's `>start>end` snarl ID and INFO/AT; scorer
  version 2 -- version 1 let a child whose trimmed POS lay left of its parent's win, and dropped
  the parent's allele). ED = the better pairing of summed unit edit distances to HG002's two
  haplotypes; chrX/chrY outside the PARs compare a homozygous pair.
- **Output** (`results/stage3/<label>/<id>.json`): `ed`, `ed_per_kb`, `exact`, `ed_ref` (CHM13 on
  both haplotypes), `gain`, `called_len`, `truth_len`, the pairing, `phase` (method, blocks,
  assignments tried, ED under the VCF's own phase), `records` (every counter), `flags`,
  `sensitivity.ed_raw_overlap` (ED under the raw overlap rule), `sensitivity.ed_suppress_nested`
  (ED when a parent's non-reference allele suppresses everything nested in it on that slot;
  `records.applied_under_nonref_parent` counts the child alleles the rule applies there), `truvari` (raw, refined, phab: TP,
  FP, FN, F1) and, with `--labels-dir`, the genome-wide truvari labels over the span and the core.
  A region errors (no ED) when the contig is not in the VCF header or when at least 20% of its
  alleles have a REF that is not CHM13's (wrong coordinates).
- **Validation** (`validate`, all 149 regions): the stvar truth VCF and the stvar+smvar union give
  ED 0 everywhere; dropped SVs, a swapped phase and a flipped GT give ED > 0; the same swap with the
  phase unknown gives ED 0 again through the phase search (exact and greedy); an empty VCF gives
  ED_ref. The raw overlap rule is kept as a check and fails where the docstring says it does.
  The truth VCFs carry no vg nesting, so the nesting clause is checked by the unit tests
  (`TestNesting`: a parent written right of its child must win; version 1 scores that case ED 40).
- **Noise-aware counts** (`summarise --noise`): beside the rule's better/worse counts, a region
  counts only when the paired difference exceeds both graphs' replicate noise summed
  (`replicate_noise.py`, truth-free). The rule's own counts are unchanged.

