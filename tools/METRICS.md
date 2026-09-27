# Metrics of a region graph (Stage 0 and Stage 1)

`tools/evaluate.py` scores one graph of one region; `tools/summarise.py` compares graphs across
regions. This page defines every number they report. The definitions follow the investigation
that motivated the work ([../docs/findings.md](../docs/findings.md)). On the same subgraphs,
`evaluate.py` reproduces its per-locus values (see *Validation* at the end).

```bash
make -C tools/c                    # builds tools/bin/libvntreval.so (evaluate.py also runs make if needed)
python3 tools/test_evaluate.py     # checks the C dynamic programs against brute force, then runs a synthetic end-to-end test
python3 tools/evaluate.py regions/L012184 regions/L012184/mc.gfa --reads                 # -> results/mc/L012184.json
python3 tools/evaluate.py regions/L012184 candidates/mafft_linsi/L012184.gfa --reads     # -> results/mafft_linsi/L012184.json
python3 tools/evaluate_all.py mc mafft_linsi --reads --jobs 2 --threads 2               # every region; keeps existing results unless --force
python3 tools/summarise.py --md results/summary.md --tsv results/summary.tsv --per-region results/per_region.tsv
```

## Inputs, validity and conventions

- **Region package**: `regions/<id>/`, containing `region.json`, `hap32.fa`, `truth.fa`, `mc.gfa` and
  `calls.tsv` (the schema is in [../docs/plan.md](../docs/plan.md)). `hap32.fa` holds every hap32 path that
  spans the region. Each is written anchor to anchor, so all of them share identical flanks. Exactly one
  of them is CHM13, and its sequence is the CHM13 span `span_start..span_end`.
- **Graph**: GFA 1.0 with S, L, and P or W lines.
  - **Name matching.** Every `hap32.fa` name must match exactly one graph path. The matcher tries, in order:
    1. the exact name;
    2. a W line written as `sample#hap#contig`, `...[start]` or `...[start-end]`. vg writes a PanSN name
       carrying a subrange this way (for example `GRCh38#0#chr4[58921381]` becomes a W line with start
       58921381);
    3. the name with a trailing `[...]` removed;
    4. a unique graph path whose name extends the name after `#`, `[` or `:`.
  - **Validity.** A graph is **valid** when the following all hold:
    - every name matches;
    - every matched path spells its `hap32.fa` sequence exactly (upper-cased);
    - no path steps onto an unknown node.

    An invalid graph gets `status: invalid`, the list of offending paths and exit status 3. No other
    metric is computed for it.
  - **Extra paths** (for example a consensus path) are allowed, counted, and otherwise ignored.
- **Handles**. A node visited forward and a node visited in reverse are different handles. The DAG is
  built on the handles that are both reachable from a path's first handle and able to reach a path's last
  handle. It uses the L lines plus every consecutive step pair of the paths; a step pair that has no L
  line is counted as `path_edges_missing_from_L`. The DAG is topologically sorted.
  - `validity.acyclic` / `size.acyclic` is false when that sort fails.
  - A cyclic graph still gets size, alignment, inflation and redundancy. It gets NA for bubbles and for the
    truth-versus-graph DP, since both need a DAG.
- **Characters**. Sequences are upper-cased. A, C, G and T are symbols, and every other byte (N, IUPAC)
  is one extra symbol N. Two bases match iff they are the same symbol, so N matches N and nothing else.
- **Paths vs walks**. A *path* is one `hap32.fa` sequence, and the panel has 34 paths when every
  haplotype spans. A *walk* is a path's handle sequence, and two paths can share a walk. Statistics are
  over paths, so a walk carried by k haplotypes counts k times. This is the investigation's "weighted
  by walk weight".
- **Span**. `span_bp` is the length of the CHM13 path, so it is the same for every graph of a region.
  Every "per kb" of span divides by it.
- **Pairs of distinct walks** are the denominator of the per-kb pair metrics. Two paths on the same walk
  are a trivial pair: cost 0, optimum 0, and not counted in any denominator.

## size

| field | definition |
|---|---|
| `nodes`, `edges`, `node_bp` | nodes visited by at least one hap32 path; distinct edges (canonical over orientation) used by consecutive path steps; total bp of those nodes |
| `gfa_nodes`, `gfa_edges`, `gfa_bp` | the same over the whole GFA (unused nodes included) |
| `nodes_per_kb`, `edges_per_kb` | `nodes` (`edges`) / `span_bp` x 1000 |
| `node_len_median`, `node_len_mean` | over the used nodes, one value per node |
| `node_frac_1bp`, `node_frac_le4bp` | fraction of used nodes of length 1 (<= 4) |
| `nonref_nodes`, `nonref_bp`, `nonref_bp_per_kb` | used nodes not on the CHM13 path, and their bp (per kb of span) |
| `distinct_walks`, `distinct_seqs`, `same_seq_different_walk` | distinct handle sequences and distinct spelled sequences among the paths, and their difference. A difference > 0 means identical haplotypes routed through different nodes |
| `bubbles`, `bubbles_per_kb` | top-level bubbles (below) |
| `bubbles_sv` | bubbles whose longest or shortest allele differs from the CHM13 allele by >= 50 bp |
| `bubbles_ge5_alleles`, `max_alleles_any_bubble` | as named |
| `largest_bubble_*` | the bubble with the largest *size*, where size = max(CHM13 bp inside, longest allele). Reported: `_alleles` (distinct traversals by the paths), `_ref_bp`, `_inner_handles`, `_frac_of_span` (CHM13 bp inside / span) |

**Top-level bubbles** (topology/topo.py):
- In the topological order, with a virtual source before every first handle and a virtual sink after
  every last handle, a handle is a **cut** iff no edge jumps over it. Then every source-to-sink path
  passes through it. This does not depend on which topological order is used.
- Two consecutive cuts that are not adjacent bound a bubble.
- The bubble's **alleles** are the distinct handle sequences the paths take between the two cuts. The
  empty traversal (a deletion edge) counts as an allele.
- The investigation checked that these bubbles are exactly vg's top-level snarls: 2,343 of 2,343 vg
  snarl IDs with both ends on cuts join consecutive cuts.

## alignment: the alignment the graph implies vs the optimum

For a pair of paths A, B:
- **Shared handles.** The handles they share, taken in A's order, are the graph's implied alignment
  anchors. In an acyclic graph they occur in the same order in B.
- **Segments.** Between consecutive shared handles, and before the first and after the last, A and B
  each have a segment of bases that are not on a shared handle (either segment may be empty).
- **Implied cost.** `cost` is the sum over segments of the optimal global edit distance between A's and
  B's segment. Bases on shared handles are free matches. This is the cheapest alignment of A and B that
  stacks every shared node on itself, so it is exactly "stack the bases that sit on the same node
  positions and align the rest as well as possible".
- **Optimum.** `dopt` is the unconstrained optimal global edit distance of A and B. Always
  `cost >= dopt`.
- **Unaligned homology.** `U` is the sum over segments of `max(|a_seg|, |b_seg|) - ed(a_seg, b_seg)`: bases
  that could be matched but sit on different nodes.
- **Pairs.** Every unordered pair of paths (561 for 34 paths). Above `--max-pairs` (default 600) the
  pairs are a seeded random sample that always keeps every CHM13 pair.

| field | definition |
|---|---|
| `all_cost_over_opt` | sum(cost) / sum(dopt) over all pairs (the investigation's locus value; hotspot median 1.56, matched controls 1.09, non-TR SV 1.00) |
| `all_excess_per_kb` | 1000 x sum(cost - dopt) / sum(L), where L = (len A + len B) / 2 and the sum runs over pairs of distinct walks (the investigation's M per kb; 223 vs 3.9) |
| `all_opt_per_kb` | 1000 x sum(dopt) / sum(L): panel diversity, which is the same for every graph of a region |
| `all_U_per_kb` | 1000 x sum(U) / sum(L): homologous sequence left on parallel nodes |
| `all_segments_per_kb` | segments per kb of pair length |
| `all_pair_ratio_{median,mean,q1,q3,max}` | per-pair `cost / dopt`, over pairs with dopt > 0 |
| `all_frac_pairs_optimal` | fraction of distinct-walk pairs with cost == dopt |
| `all_pairs_same_seq_diff_walk` | pairs with the same sequence but different walks. Their cost is > 0 and their dopt is 0, so they add to the excess but not to the per-pair ratios |
| `ref_*` | the same over only the pairs that contain CHM13, which is each path against the reference |
| `affine_all_*`, `affine_ref_*` | the same under gap-affine scoring (below) |
| `affine_pairs` | `all_pairs`, or `ref_pairs_plus_<k>_sampled` when the budget forced a sample |

**Unit cost**: global Levenshtein distance, computed exactly with the Myers/Hyyro bit-vector algorithm
(`ed_unit` in `tools/c/vntreval.c`).

**Gap-affine** (topology/affine.py):
- **Costs.** A mismatch costs 4 and a gap of k bases costs 6 + 2k; implied cost and optimum both use them.
- **DP.** Exact Gotoh (`ed_affine`), banded around the diagonals 0..(len B - len A) with a band that doubles
  until the banded optimum is provably global:
  - an alignment that leaves a corridor of half-width w has at least 2(w+1) + |len A - len B| gap bases,
    so it costs at least 6 + 2(2(w+1) + |dlen|);
  - once the banded optimum is at or below that bound, it is the global optimum.
- **Budget.** Pairs get the affine optimum while the distinct pairs need at most `--affine-budget` full
  DP cells (sum of len A x len B, default 5e10). Above that budget, the CHM13 pairs are scored plus a
  seeded sample of the others that fits.
  - The choice depends only on the sequences, so every graph of a region is scored on the same pairs.
  - The investigation reported the CHM13 pairs: `affine_ref_cost_over_opt`, hotspots 1.64 vs 1.03.

The optima do not depend on the method. They are cached per region in
`$VNTR_WORK/cache/evaluate/<id>/optima.json`, keyed by sequence hashes, so a second graph of a region pays
only for its segments.

## inflation: each path written against CHM13 through the graph

This is `representability/inflation.py` and `oracle_truvari.deconstruct`.
- **Anchors.** Walk each non-CHM13 path. The CHM13 handles it visits, taken greedily in increasing CHM13
  order, are its anchors.
- **Records.** Each stretch between consecutive anchors whose handles differ from the CHM13 handles
  between the same anchors becomes one record (REF = the CHM13 bases, ALT = the path's bases). A record is
  padded with the preceding reference base and trimmed of a common suffix, then of a common prefix. This
  is the representation a graph genotyper (vg deconstruct / vg call) writes.

| field | definition |
|---|---|
| `sv_pieces_{median,mean,q1,q3,max}` | per path, the records with abs(len REF - len ALT) >= 50: "SV-sized pieces per haplotype" (hotspots ~9 vs 1 at controls in the investigation, whose median was over distinct walks) |
| `paths_with_sv_piece`, `sv_pieces_total` | as named |
| `records_*` | records per path, any size |
| `indel_bp_graph` | sum over paths and records of abs(len REF - len ALT): the indel bp the graph's decomposition uses |
| `indel_bp_net` | sum over paths of abs(len path - len CHM13): the least indel bp any alignment needs |
| `indel_bp_ratio` | `indel_bp_graph / indel_bp_net` (1 = no compensating insertion/deletion pairs) |
| `indel_bp_excess_per_kb` | 1000 x (graph - net) / sum over paths of (len path + len CHM13)/2 |
| `anchored_over_free` | sum over paths of the edit distance through the records (sum of ed(REF, ALT)), divided by the sum of the free edit distance path vs CHM13, over paths with free > 0 (the investigation's walk inflation; 1.64 vs 1.02 medians of the per-walk values) |
| `anchored_over_free_path_*` | the per-path ratios |

## redundancy

| field | definition |
|---|---|
| `kmer_positions` | over all 21-mers K spelled by the paths (without N), the sum over K of P(K), the number of distinct graph positions of K. A position is the (node, forward offset, strand) of K's first base |
| `kmer_extra_positions` | sum over K of P(K) - C(K), where C(K) is the most times any single path contains K. A graph that merges every shared copy needs only C(K) positions |
| `kmer_frac_extra` | extra / positions |
| `kmer_extra_per_kb` | extra per kb of span |
| `kmers`, `kmers_with_extra` | distinct 21-mers, and those with P > C |
| `U_per_kb` | = `alignment.all_U_per_kb`: homologous bp left on parallel nodes |

`kmer_*` is the viewer's measure (viewer/view_region.py; the viewer keyed a reverse-strand base by its
offset along the walk, here it is keyed by its forward offset) and should be read as an upper bound on unmerged
sequence. Reordered repeat copies also add to it. A unit inserted at the end of an array leaves up to
about one unit of k-mers that cross from the last reference unit into the flank and also exist
on the inserted unit's node.

## truth: can the graph spell HG002, and how does it write it

`truth.fa` holds HG002 hap1 and hap2 (T2T-Q100 v1.1) over the same anchor-to-anchor span. It is for
evaluation only.

| field | definition |
|---|---|
| `h{1,2}_d_graph` | minimum unit edit distance from the truth haplotype to the sequence of any source-to-sink path through the graph. Sources are the paths' first handles and sinks their last handles, and the path is any path through the DAG, recombinations included. The DP is exact (`dag_align` in `vntreval.c`): per handle, the column of the edit DP is propagated through the handle's bases, and a handle's entry column is the element-wise minimum of its predecessors' exit columns. A traceback recovers the handle path, and the script checks that the path's sequence is at that distance. Among equally good paths (common where identical sequence sits on different nodes), the traceback prefers CHM13 handles, then handles on more panel paths. The truth is therefore written as close to the reference as the graph allows, and the choice does not depend on the order of the GFA's lines. 0 means the graph can spell the haplotype |
| `h{1,2}_graph_path_is_panel_path` | the best path found is one of the panel paths |
| `h{1,2}_d_panel`, `h{1,2}_closest_panel_path` | minimum edit distance to any hap32 path, which is the same for every graph |
| `h{1,2}_d_ref` | edit distance to CHM13 |
| `d_graph_sum`, `d_panel_sum`, `h*_panel_minus_graph` | sums and the gain from recombination |
| `truvari.*` | the truth written by the graph (below) |

**The truth written by the graph** (representability/oracle_paths.py + oracle_truvari.py):
- **Decomposition.** Each truth haplotype's best graph path is decomposed against CHM13 exactly as in
  *inflation*. Its sequence is the truth itself when d_graph = 0.
- **VCF.** The records are written as one VCF, with hap1's records `1|0`, hap2's `0|1`, and shared records
  `1|1`. It is run through `bcftools norm -m-any -f` and sorted.
- **Scoring.** `truvari bench` is run with the pipeline's parameters (`--sizemin 50 --sizefilt 50 --pick ac`,
  samples HG002) against `work/wgs-mm095/score/<contig>.truth.norm.vcf.gz`. The region is the span
  intersected with the SV benchmark BED.
- **Refine, twice.**
  1. `truvari refine -u -t 1` (phab/POA harmonisation) on truvari's default candidate regions (the
     FP/FN clusters plus 100 bp), as the investigation ran it. Without `-w`, a region's original records
     become TP only when the harmonised region has no FP/FN left, otherwise they keep their bench labels,
     so the credit is all-or-nothing per region. It is also sensitive to where a region's boundary cuts
     the records: the same exact truth sequence, written two ways, went from 5/0/9/0 to 4/1/4/5 at L014297.
  2. `truvari refine -u -w -t 1 --regions <span ∩ SV benchmark>`, on a pristine copy of the bench
     directory. It harmonises the whole region and counts in phab's representation, so it depends only
     on the two haplotype sequences the graph writes and not on how the graph decomposes them: the
     exact truth scores 1.0, and a near-exact path loses only its differences. This is the
     representation-free score (`phab_*`).

| field | definition |
|---|---|
| `records`, `records_sv` | records written (any size / >= 50 bp) |
| `raw_TP-base, raw_FN, raw_TP-comp, raw_FP, raw_f1` | truvari bench; F1 = 2PR/(P+R), with recall TP-base/(TP-base+FN) and precision TP-comp/(TP-comp+FP) |
| `refined_*` | from `refine.variant_summary.json` of refine run 1. `refine_status` is `ok`, `nothing_to_refine` (no FP/FN, so the raw counts are copied), `timeout` (after `--refine-timeout`, default 900 s) or `failed` |
| `phab_*` | the same from refine run 2 (whole span, phab representation); `phab_status` as above |
| `status` | `ok`, `span_outside_sv_benchmark`, `no_truth_graph_paths` or a missing input |

**How to read it.**
- The exact truth, written the way the hap32 graph breaks it up, scores pooled raw F1 0.27 at hotspots,
  and 0.87-0.89 after refine run 1.
- A better graph should raise raw F1, which is the score a genotyper that picked the right haplotypes
  would get with this graph's representation.
- `phab_*` checks only that the right sequences can be written; `h*_d_graph` measures the same thing in
  edits.
- `summarise.py` also pools TP/FP/FN over a stratum.
- truvari refine's own memory peaks at about 5 GB on the largest (19 kb) regions.

## reads (`--reads`, Stage 1): placement redundancy of the existing short reads

This is `reads/reads_analysis.py placements()`.
- **Reads.** The region's reads are `<region>/reads.gaf(.gz)` if present. Otherwise they come from
  GAF-Base: `gbz-base query --between anchor_left:anchor_right --alignments overlapping`, cached in
  `$VNTR_WORK/cache/evaluate/<id>/reads.gaf`.
- **Filter.** A read is tested when all of its alignment path lies on nodes of the baseline `mc.gfa` and
  it places at least one aligned base on a **core** node. Core nodes are CHM13 nodes overlapping
  `core_start..core_end`, plus non-CHM13 nodes that some path places in the core when each run of
  non-CHM13 nodes is spread linearly between the CHM13 positions around it.
- **Sequence.** The read's aligned graph sequence (the baseline path's sequence `[ps, pe)`, at least 30 bp)
  is the query. Every graph is therefore tested with the same sequences.
- **Search.** In the graph under test, the query and its reverse complement are searched exactly in
  every hap32 path. Each occurrence becomes a **placement**: its handle sequence plus the start offset in
  the first handle and the end offset in the last, made canonical over the two orientations.
- **Counts.** `places` is the number of distinct placements. `max_occ` is the most occurrences in any
  single path.
- **Offsets.** The investigation keyed a placement by its node path alone. Adding the offsets makes
  `places >= max_occ` in an acyclic graph (the ideal ratio is 1). It also stops a compacted graph, whose
  long nodes hold several tandem copies, from counting those copies as one placement.

| field | definition (fractions of tested reads) |
|---|---|
| `*_unique` | places == 1 |
| `*_tandem` | max_occ >= 2 (some haplotype holds it twice: ambiguity no graph can remove) |
| `*_redundant_only` | places > 1 and max_occ == 1: the same sequence on different nodes in different haplotypes, which a better-aligned graph removes (investigation: 31% of hotspot MAPQ<5 reads) |
| `*_redundant` | places > max(1, max_occ): any excess placement |
| `*_absent` | places == 0 (the read's baseline sequence is spelled by no panel path, e.g. it crosses a recombination) |
| `*_ratio_mean`, `*_ratio_median` | places / max(1, max_occ), over reads with places > 0 |
| `*_excess_placements_mean` | places - max_occ |

The prefix `all_` covers every tested read, `mqlt5_` MAPQ < 5 and `mqge30_` MAPQ >= 30. The MAPQ is the
one the read had in the hap32 mapping. `core_reads_tested`, `outside_subgraph` and `gaf_records` are
counts.

## The full-panel arm (`--panel`)

`evaluate.py <region_dir> <graph> --panel <panel.fa[.gz]>` judges a graph whose paths are an
arbitrary panel, typically a full-panel graph `work/panel/<method>/<id>.gfa` against
`regions/<id>/hprc.fa.gz` (built by [panel.py](README.md#panelpy-the-full-panel-arm)).

- **Validity**: every FASTA record must be a path of the same name (the same name matching as
  for hap32) spelling it exactly; otherwise exit 3. The CHM13 record is the reference. Names of
  HG002, HG003, HG004 (or their NA aliases) and gref_CHM13 are refused.
- **Skipped**: the hap32, truth and reads sections.
- **size**, bubbles and **redundancy** are computed on the distinct walks (their values do not
  depend on how many haplotypes share a walk); `size.n_paths` gives the panel size.
- **alignment**: `all_*` over a seeded uniform sample of `--max-pairs` (default 300) unordered
  haplotype pairs, `ref_*` over CHM13 paired with a seeded sample of `--max-pairs` other
  haplotypes; scored in one pass, aggregated apart (`pairs_mode` says so). Frequent walks weigh
  more, as duplicated hap32 paths do. Pairs are drawn by record index, so two graphs of one panel
  FASTA are scored on the same pairs. `*_cost_over_opt_ci95` is a 1000-fold bootstrap over the
  sampled pairs; it understates the between-seed spread at hotspots, where a few very long
  haplotypes dominate the sums.
- **inflation** over every haplotype; `anchored_over_free` uses the CHM13 optima of the `ref_*`
  sample only (`anchored_over_free_paths`).
- Output defaults to `results/panel/<METHOD>/<region_id>.json`; `evaluate_all.py <methods> --panel`
  runs it over every region.

## summarise.py

Per stratum and metric (METRICS list in the script, each with its direction), it reports:
- per method: n, median [q1-q3] and mean;
- against the baseline method (`mc`) on the regions both have: the median and IQR of the per-region
  difference, how many regions got better / worse / stayed equal, and a two-sided Wilcoxon signed-rank p;
- truvari F1 pooled over the stratum's regions (sum of TP/FN/FP), raw and refined;
- in the headline table, the number of regions at the plan's target, `all_cost_over_opt <= 1.1`.

The Wilcoxon test drops zero differences and averages tied ranks. It uses the exact permutation
distribution of the observed ranks up to 60 pairs, and a tie-corrected normal approximation above that.
Only valid graphs (`status: ok`) are used; invalid ones are counted per method.

## Validation

**Against the investigation, on its own subgraphs.** `evaluate.py` was run on the hap32 subgraph of
all 150 stratified loci, rebuilt from the investigation's `topology/regions/*/subgraph.gfa` with the same
anchor-to-anchor spans. Per locus, the investigation's rounding is the tolerance.

| investigation table | columns | loci that agree |
|---|---|---|
| `topology/topology_loci.tsv` | `all_pair_cost_over_opt`, `all_pair_M_per_kb`, `all_pair_U_per_kb`, `all_pair_dopt_per_kb`, `ref_pair_cost_over_opt`, `ref_pair_M_per_kb`, bubbles, SV bubbles, largest-bubble alleles | 150 / 150 each |
| `topology/affine.tsv` | `affine_graph_over_opt`, `affine_excess_per_kb` (as `affine_ref_*`) | 89 / 89 |
| `representability/inflation_loci.tsv` | `walks_inflation` (as `anchored_over_free`) | 149 / 149 |
| `representability/oracle_paths.tsv` | truth-path edit distance, h1 and h2 | 149 / 149 |
| `representability/oracle_paths.tsv`, `refine.tsv` | truth-by-graph truvari TP-base | raw 145 / 149, refined 145 / 149 |

- **Why truvari differs at 4 loci.** The DP returns one of several equally good paths through the graph
  (identical sequence on different nodes). The investigation took the first one in GFA order; this tool
  takes the CHM13-first one.
- **nodes_per_kb** differs at 5 loci, because the investigation also counted nodes that only fragments
  use.

Per-stratum medians of `all_cost_over_opt` [IQR], against the investigation's table:

| subgraphs | hotspot VNTR | matched control | correct control | hotspot other | non-TR SV |
|---|---|---|---|---|---|
| investigation subgraphs | 1.56 [1.34-1.81] | 1.09 [1.01-1.32] | 1.01 [1.00-1.13] | 1.66 | 1.00 |
| region packages (`results/mc`) | 1.57 [1.34-1.81] | 1.10 [1.01-1.32] | 1.01 [1.00-1.13] | 1.63 | 1.00 |

On the packages, hotspot VNTRs also show:
- excess edits per kb 182 against 4.4 at matched controls;
- SV pieces per path 8 against 1;
- k-mer extra fraction 0.60 against 0.05.

Pooled truth-by-graph truvari F1 at hotspot VNTRs is 0.262 raw, 0.867 refined and 0.989 phab (the
investigation reported 0.271 raw and 0.892 refined). Pooled MAPQ<5 read classes at hotspot VNTRs are
unique 12.7%, tandem 57.6% and redundant 28.8%, against the investigation's 13.2%, 56.1% and 30.7%. That
comparison is on a different read window, and placements here include offsets.

**Where definitions differ from the investigation's tables:**
- `nodes_per_kb` counts only nodes that a spanning path uses. The investigation also counted nodes used
  only by fragments.
- Per-path medians weight each walk by its multiplicity. `inflation.py`'s medians were over distinct
  walks.
- Read placements include the start and end offsets (see *reads*).
- The truth-path tie-break prefers CHM13 (see *truth*).

**Cost.**
- A typical region takes 3-15 s. The largest take 1-3 min: 19 kb spans of ~40k nodes (L015347), and the
  55 kb L012272.
- Pairwise optima are cached per region, so a second graph of the region skips them.
- `evaluate.py` itself stays under ~1.5 GB. The truth DP is banded by the distance to the closest panel
  path, which is exact (cells above a valid upper bound cannot be on an optimal alignment).
- truvari refine peaks at ~5 GB on the largest regions.
