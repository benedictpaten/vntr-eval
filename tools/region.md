# tools/region.py: what it extracts and every output

Reference for `tools/region.py` (a port of the investigation's `vntr_region.py`). Data and tool
locations come from `tools/config.py`; the exact edit distance uses `tools/fastedit.c`, compiled on
first use into `config.CACHE_DIR` (pure Python when no C compiler is available, same results).

## Why the subgraph is not `--interval --context 0`

`gbz-base query --interval A..B --context 0` returns **only the nodes of the
CHM13 walk**. On chr20:1902001-1904000, all 1,205 nodes are CHM13 nodes. Every
haplotype that deviates from CHM13 therefore leaves the subgraph and re-enters
it, and shows up as hundreds of 2-27 bp fragments (1,242 distinct walks there).
That looks like a tangle, but it is an extraction artifact.

In the default `--mode between`, the tool picks two CHM13 anchor nodes, L and R,
that flank the padded interval, and extracts everything between them with
`gbz-base query --between L:R`. An anchor must:

* carry the maximum haplotype coverage seen in a probe window;
* be visited at most once by every haplotype subpath, and exactly once by CHM13;
* not be touched by any truth record (stvar or smvar), so that the truth and
  panel haplotypes share identical ends.

L is the last such node ending at or before the padded start; R is the first
such node starting at or after the padded end. The probe window starts at 20 kb
and doubles, up to 320 kb. Spans over `--max-span` (250 kb) are refused.

The subgraph is then bounded by L and R. Every haplotype that walks through L
and R appears as one complete walk. Walks only fragment where a path fragment
genuinely starts or ends inside the span, or where a walk cycles back out
through a boundary.

The analysed span (L start to R end) is usually larger than the requested
interval. Both are reported. `--mode interval [--context N]` reproduces the old
extraction, with a warning.

The graph holds 34 haplotype paths: 32 sampled `recombination#k` haplotypes
(themselves split into fragments `recombination#k#chrN#f`), plus the CHM13 and
GRCh38 reference paths. `H` is the coverage used for the anchors, normally 34.
gbz-base W lines never carry path names (they are written as `unknown`, because
GBWT records have no path identity). Names come only from `--names`, which
streams `vg paths -A` over the per-contig GBZ and re-derives the same subpaths.
The run checks that the two multisets are equal (`named_runs_equal_gbz_subpaths`).

## Usage

```
python3 tools/region.py <contig> <start> <end> <outdir> [--pad N] [--no-reads] [--names]
        [--anchors L+:R+]
        [--mode between|interval] [--context N] [--depth 30] [--msa] [--bandage]
        [--window 20000] [--max-span 250000] [--limit 300000] [--max-align-cells 4e8] [--quiet]
```

`<start>` and `<end>` are 1-based and inclusive, as in a VCF. Internally the
padded interval is `[start-1-pad, end+pad)`, 0-based half-open, which is
gbz-base's `--interval` convention.

| option | effect |
|---|---|
| `--pad N` | widen the interval by N bp on each side before choosing anchors |
| `--no-reads` | skip the GAF-base read query |
| `--names` | name every subpath via `vg paths -A` (slow; `--mode between` only) |
| `--anchors L+:R+` | use these two CHM13 nodes as L and R instead of choosing them (e.g. a region package's `anchor_left`/`anchor_right`); anything that breaks the anchor rule is reported as a warning |
| `--msa` | `mafft --auto` alignment of all spanning distinct walks plus both truth haplotypes, written to `msa.fa` |
| `--bandage` | render `subgraph.gfa` to `subgraph.png` |
| `--max-align-cells` | skip the CIGAR traceback when len(ref) × len(hap) exceeds this; distances stay exact |

From Python:

```python
import sys
sys.path.insert(0, 'tools')
import region as vr

s = vr.analyse('chr20', 1902001, 1904000, 'out', reads=False)   # returns summary dict
vr.edit_distance(a, b)        # exact Levenshtein, Myers bit-vector (10 kb x 10 kb in 0.05 s)
vr.align_cigar(ref, hap)      # (distance, CIGAR with = X I D)
vr.tandem_scan(seq)           # truth-independent STR/VNTR period scan
vr.Gfa(open('subgraph.gfa').read())
```

Printed summary: requested and analysed span; graph size; cycles; truth
haplotype lengths; the closest panel haplotype to each truth haplotype; reads
observed vs expected; call status; each caller's haplotypes scored against
truth; the CHM13 tandem scan; warnings.

### Output files

| file | content |
|---|---|
| `summary.json` | all headline numbers (keys below) and `warnings` |
| `run.log` | every command run, with its stderr |
| `subgraph.gfa` | annotated distinct-haplotype GFA (S, L and W lines; see below) |
| `subgraph.query.gfa` | raw gbz-base output of the `--distinct` query |
| `subpaths.gfa` | raw gbz-base output without `--distinct`: one W line per individual haplotype subpath |
| `haplotypes.tsv`, `haplotypes.fa` | one row / record per distinct walk |
| `subpaths.tsv` | one row per individual (non-distinct) subpath |
| `haplotype_names.tsv` | `--names` only: one row per haplotype (sample#hap) |
| `topology.json` | graph and walk statistics |
| `truth.hap1.fa`, `truth.hap2.fa` | HG002 sequence of the analysed span, per haplotype |
| `truth_records.tsv` | every truth record overlapping the span, and what was done with it |
| `closest.json`, `closest.tsv` | closest panel haplotype to each truth haplotype |
| `vg.called.fa`, `pangenie.called.fa` | each caller's genotype applied to CHM13 over the span (two slots) |
| `calls.tsv` | vg and PanGenie alleles in the span, with truvari and aardvark status |
| `reads.gaf`, `reads_summary.json` | read alignments overlapping the subgraph, and their summary |
| `msa.fa`, `subgraph.png` | with `--msa` / `--bandage` |

`subgraph.gfa` W lines:

* The reference walk is written as
  `W CHM13 0 <contig> <span start0> <span end> <walk> WT:i:<weight>`.
* Every other distinct walk is written as
  `W h<i> 0 <class> 0 <length> <walk> WT:i:<weight> CG:Z:<cigar>`. The CG tag
  is the base-level alignment against the CHM13 span, and is present for
  spanning walks only.
* `h<i>` matches `hap_id` in `haplotypes.tsv`.

### haplotypes.tsv (and haplotypes.fa)

One row per distinct walk, as `gbz-base --distinct` reports them. Rows are
ordered reference first, then spanning walks by weight and length, then fragments.
Walks are put into CHM13 orientation, so sequences read in reference orientation.

| column | meaning |
|---|---|
| `hap_id` | `h0`, `h1`, ...; `h0` is the CHM13 walk when found |
| `weight` | number of the 34 haplotype paths with exactly this walk |
| `is_reference` | 1 if identical to the CHM13 walk (weight-1 other paths share it) |
| `class` | `spanning` = starts at L and ends at R in reference orientation; `enters_L` = starts at L, ends elsewhere; `exits_R` = ends at R, starts elsewhere; `internal` = neither |
| `length_bp` | sum of node lengths |
| `len_minus_ref` | `length_bp` minus the CHM13 span length (spanning walks only) |
| `n_steps` | node visits |
| `n_distinct_nodes` | distinct nodes visited |
| `n_nonref_nodes` | distinct nodes not on the CHM13 walk |
| `nonref_bp` | bp of steps on such nodes |
| `n_reverse_steps` | steps that traverse a node opposite to CHM13's orientation (reverse, for non-CHM13 nodes) |
| `max_node_visits` | most visits of one node by this walk (>1 means the walk is a cycle) |
| `n_nodes_revisited` | nodes this walk visits more than once |
| `flipped` | 1 if gbz-base reported the walk reversed and it was flipped to reference orientation. Orientation is chosen by bp-weighted agreement with CHM13 node orientations; gbz-base keeps only a "canonical" orientation of each walk. |
| `starts_at`, `ends_at` | boundary label of the first / last step: `L`, `L-`, `R`, `R-` (reversed), `inner` |
| `edit_to_ref`, `identity_to_ref` | unit-cost global edit distance to the CHM13 span; identity = 1 - d / max(len) |
| `edit_to_truth_h1`, `identity_to_truth_h1`, `edit_to_truth_h2`, `identity_to_truth_h2` | the same against each truth haplotype (spanning walks only) |
| `first_node`, `last_node` | first / last step, e.g. `>114943418` |
| `cigar_gbz` | `--mode interval` only: gbz-base's node-level CIGAR (`reversed:` prefix if flipped) |
| `cigar` | base-level unit-cost global alignment of this walk (query) against the CHM13 span (reference): `=` match, `X` mismatch, `I` bases only in the walk, `D` bases only in CHM13. Spanning walks only; `NA:too_long` past `--max-align-cells`. |

FASTA headers read `>h<i> weight=<w> class=<class> len=<bp> [reference=CHM13]`.

### subpaths.tsv

One row per individual haplotype subpath from the non-distinct query. The
weights in `haplotypes.tsv` sum to this row count.

| column | meaning |
|---|---|
| `subpath_id` | `s0`, `s1`, ... |
| `distinct_id` | `hap_id` of the identical distinct walk |
| `class`, `length_bp`, `n_steps`, `max_node_visits`, `n_nodes_revisited`, `flipped`, `starts_at`, `ends_at` | as in `haplotypes.tsv` |
| `haplotype` | `--names` only: `sample#hap`, e.g. `recombination#5`, `CHM13#0`, `GRCh38#0` (assignment among identical walks is arbitrary) |
| `path_name` | `--names` only: full GBZ path name, including the fragment index |

### haplotype_names.tsv (`--names`)

| column | meaning |
|---|---|
| `haplotype` | `sample#hap` |
| `path_fragments` | GBZ path fragments that have subpaths in the span |
| `n_subpaths` | subpaths of this haplotype in the span (1 for a clean traversal) |
| `n_spanning` | how many of them are L-to-R traversals |
| `classes` | `class:count` list |
| `total_bp`, `spanning_bp` | bp over all subpaths / length of each spanning traversal |
| `has_cycle`, `max_node_visits` | whether any node is visited more than once by one subpath, and the maximum visits |
| `n_nodes_revisited` | nodes visited more than once, counted over all this haplotype's subpaths |
| `distinct_ids` | `hap_id` of each subpath |
| `fragment_ends_inside` | subpath ends on an inner node, i.e. a path fragment starts or ends inside the span |

### topology.json

| key | meaning |
|---|---|
| `mode`, `query`, `analysed_span_1based` | how the subgraph was extracted |
| `nodes`, `edges`, `total_node_bp` | subgraph size |
| `nodes_per_kb_reference` | nodes / CHM13 span kb |
| `reference_path_length_bp`, `reference_steps`, `reference_distinct_nodes`, `reference_nodes_visited_more_than_once` | the CHM13 walk |
| `nonreference_nodes`, `nonreference_node_bp` | nodes not on the CHM13 walk |
| `distinct_haplotype_walks`, `haplotype_walks_total_weight`, `individual_subpaths` | walk counts |
| `spanning_distinct_walks`, `spanning_walks_weight`, `fragment_walks_weight`, `fragment_distinct_walks` | fragmentation |
| `haplotypes_at_boundaries_H` | haplotype coverage of the anchors |
| `cycles.distinct_walks_revisiting_a_node`, `cycles.weight_of_walks_revisiting_a_node`, `cycles.max_visits_of_one_node_by_one_walk`, `cycles.nodes_revisited_by_some_walk` | walk-level cycles |
| `cycles.self_loop_edges`, `cycles.reversing_edges` | edges whose two ends are on the same node / in opposite orientations |
| `cycles.graph_has_directed_cycle` | whether the subgraph, viewed as a directed graph on handles, has any cycle |
| `named` | `--names` only: `n_paths_scanned`, `n_named_runs`, `named_runs_equal_gbz_subpaths`, `n_haplotypes_touching`, `n_haplotypes_single_spanning_traversal`, `n_haplotypes_fragmented` (no spanning traversal), `n_haplotypes_with_cycle`, `n_haplotypes_multiple_subpaths` |

### Truth haplotypes and truth_records.tsv

The truth haplotypes are rebuilt over the analysed span (not the requested
interval), so they have the same ends as the spanning panel walks. The steps:

1. Collect records from `CHM13v2.0_HG2-T2TQ100-V1.1_stvar.vcf.gz` first, then
   `..._smvar.vcf.gz`.
2. For each haplotype, apply the allele named by GT slot 1 or 2.
3. Skip `*` (spanning-deletion) alleles, symbolic alleles, and records crossing
   the span edge.
4. Skip smvar records that duplicate or overlap an applied stvar record on that
   haplotype. In SV regions the two files repeat each other: stvar splits
   smvar's multi-allelic `1|2` records into `1|0` / `0|1` records.
5. Skip any later record that overlaps an applied record on the same haplotype.
6. Write the resolved records to one VCF, with GT per haplotype, and run
   `bcftools consensus -H 1` and `-H 2` against the span FASTA.
7. As a cross-check, apply the same records in Python. Any mismatch goes into
   `warnings`; the result has always agreed so far.

| column | meaning |
|---|---|
| `source` | `stvar` or `smvar` |
| `chrom`, `pos`, `end` | `end` = POS + len(REF) - 1 |
| `ref_len`, `alt_lens` | allele lengths |
| `max_allele_len_diff` | max \|len(ALT) - len(REF)\| over real alleles (excludes `*`) |
| `is_sv50` | `max_allele_len_diff` >= 50 |
| `svtype`, `svlen`, `gt` | from the record |
| `applied_h1`, `applied_h2` | `yes`, `ref`, or `no:<reason>`, where the reason is one of `star_allele`, `symbolic`, `crosses_span_edge`, `duplicate_of_stvar`, `overlaps_applied_stvar`, `overlaps_applied_smvar`, `missing_gt` |
| `within_span` | record lies entirely inside the analysed span |
| `in_requested_interval` | overlaps the requested interval |
| `in_sv_benchmark` | REF span fully inside the SV benchmark BED; truth outside the BED is unreliable, and truvari does not score it |
| `vg_truvari`, `pg_truvari` | truvari base status (`TP` = tp-base, `FN` = fn) for vg / PanGenie, joined after the scoring pipeline's `bcftools norm -m-any -f ref`. `-` = not scored (<50 bp, outside the BED, or an smvar row). |
| `TRF`, `TRFperiod`, `TRFrepeat`, `TRFstart`, `TRFend`, `TRFcopies`, `TRFdiff`, `LCR`, `REMAP`, `RM_clsfam` | truth INFO annotations. `LCR` is an entropy score, so high means complex sequence; period >= 7 is a VNTR, 1-6 an STR. |

`summary.json:truth` also records the haplotype lengths, the counts of stvar
SVs >= 50 bp (total, TRF period >= 7, in the requested interval), the SV and
small-variant benchmark coverage of the interval and of the span, and the edit
distances h1-ref, h2-ref and h1-h2.

### closest.json / closest.tsv

For each truth haplotype, the spanning distinct walks are ranked by exact edit
distance, with ties broken by length difference. All spanning walks are scored;
no approximation is used.

* `closest.json` per truth haplotype:
  * `truth_len`, `method`;
  * `best` (`hap_id`, `weight`, `is_reference`, `edit_distance`,
    `length_diff_panel_minus_truth`, `identity`);
  * `best_non_reference`;
  * `reference_edit_distance`;
  * `n_exact_matches_weight` (panel paths identical to truth);
  * `top5` as (hap_id, weight, edit, length difference).
* `closest.tsv` has one row per spanning distinct walk, with its distances to
  CHM13 and to both truth haplotypes.

### calls.tsv

One row per ALT allele. vg records come from
`work/wgs-mm095/<contig>/<contig>.vcf.gz`, overlapping the span. PanGenie rows
come from `work/pangenie/score/<contig>.norm.vcf.gz` and include called alleles
only, because PanGenie genotypes every panel allele.

To join vg alleles to the truvari output, each allele is left-aligned exactly as
the scoring pipeline did (`bcftools norm -m-any -f ref`, split per allele here so
the mapping is exact). The allele is then looked up by (POS, REF, ALT) in
`tp-comp.vcf.gz` / `fp.vcf.gz`. Small variants are looked up in aardvark's
`query.vcf.gz` (FORMAT/BD).

| column | meaning |
|---|---|
| `caller` | `vg` or `pangenie` |
| `rec` | record index within this span |
| `chrom`, `pos`, `end`, `id` | `id` is the snarl for vg |
| `ref_len`, `allele`, `alt_len`, `len_diff` | allele index and lengths |
| `called` | the allele appears in GT |
| `gt`, `qual`, `filter` | from the record |
| `in_interval` | overlaps the requested interval |
| `norm_pos`, `norm_ref_len`, `norm_alt_len` | after left-alignment |
| `truvari` | `TP`, `FP`, or `-` = not scored (<50 bp, outside the benchmark BED, or monoallelic) |
| `aardvark` | BD value (`TP` / `FP` ...), or `-` = not scored |
| `PctSeqSimilarity`, `PctSizeSimilarity`, `SizeDiff`, `StartDistance`, `MatchId` | truvari INFO of the matched record; for FPs these describe the nearest truth SV |
| `overlaps_span_edge` | the record extends beyond the analysed span |

`summary.json:calls` counts called alleles as `sv50_TP` / `sv50_FP` / `sv50_-`
(\|len_diff\| >= 50, truvari) and `small_<aardvark status>`.
`summary.json:truth_sv50_status` counts stvar SVs >= 50 bp by vg / PanGenie
truvari status.

### Called haplotypes vs truth (summary.json:called_haplotypes_vs_truth)

This is a representation-independent score. Each caller's genotypes are applied
to CHM13 over the span, one sequence per GT slot, and written to
`<caller>.called.fa`. Records are applied longest first at equal POS, so a
nested vg parent record wins over its `.|1` children. Records crossing the span
edge and records overlapping in the same slot are skipped and counted.

The two slots are paired with the truth haplotypes by the smaller summed edit
distance:

* `d_h1`, `d_h2`, `total`, `pairing`;
* `slot1_len`, `slot2_len`, `n_applied`;
* skip counters: `skipped_edge`, `skipped_overlap`, `skipped_other`,
  `missing_allele`;
* `unphased_het`, `phase_sets`, `phase_reliable`.

`reference_as_call` is the baseline of calling nothing.

PanGenie's genotypes are unphased. With more than one het, the slot assignment
is arbitrary, and the result is flagged `phase_reliable: false`. A caller whose
`total` exceeds `reference_as_call` made the region worse than doing nothing.

### reads.gaf / reads_summary.json

`reads.gaf` is `gbz-base query --between ... --gaf-base work/reads.hap32.gaf.db
--alignments overlapping`: full alignments overlapping any subgraph node. The
backend is GAF-Base.

| key | meaning |
|---|---|
| `n_alignments`, `n_distinct_read_names` | mates share a name |
| `mapq_histogram` | bins 0, 1-4, 5-9, 10-19, 20-29, 30-59, 60, 61-255 |
| `frac_mapq_lt5`, `frac_mapq0`, `median_mapq` | MAPQ summaries |
| `strand_vs_reference`, `strand_by_mapq` | read orientation relative to CHM13 (bp-weighted agreement of path node orientations), and the same split by MAPQ < 5 / >= 5 |
| `median_read_length` | median query length |
| `median_block_identity`, `frac_block_identity_lt_0.9`, `frac_block_identity_lt_0.5` | GAF matches / block length |
| `frac_query_minus_path_span_ge20`, `frac_path_minus_query_span_ge20` | read bases not placed on the path (insertions) / path bases not covered by the read (deletions), >= 20 bp |
| `n_path_inside_subgraph`, `n_path_leaves_subgraph`, `frac_path_inside_subgraph` | whether every node of the read's path is in the subgraph |
| `n_touching_nonref_nodes` | reads touching a node that is not on the CHM13 walk |
| `expected_from_reference_span` | `depth × (L + r - 1) / r`, with L = CHM13 span and r = median read length; depth is `--depth`, default 30 |
| `expected_from_truth_haplotypes` | the sum over the two HG002 haplotypes of `(depth / 2) × (L_h + r - 1) / r`. Use this one in VNTRs, where the HG002 length can be far from CHM13's. |
| `observed_over_expected_truth`, `observed_over_expected_reference` | observed / expected |

### Other summary.json keys

| key | meaning |
|---|---|
| `boundaries` | L and R nodes with coordinates, `H`, the probe window, the rule |
| `fragmentation` | spanning vs fragmentary walk weight |
| `closest`, `closest_non_reference` | copies of `closest.json:best` |
| `repeatmasker` | overlapping RepeatMasker elements, including Simple_repeat bp with period >= 7 and 1-6. See the caveat below. |
| `reference_tandem_scan`, `truth_h1_tandem_scan` | `tandem_scan()` of the CHM13 span and of truth hap1 |

`tandem_scan()` works as follows:

1. For each period p (1-200), a window of max(2p, 40) bp is repetitive when at
   least 80% of its bases equal the base p downstream.
2. Covered bp = matching bases, plus their copies, inside such windows.
3. The best STR (p 1-6) and VNTR (p >= 7) period by covered bp is reported,
   replaced by its smallest divisor that covers at least 85% as much.
4. Sequences over 50 kb are scanned on their central 50 kb.

It detects the presence of a VNTR reliably. The reported period can differ from
TRF's for compound repeats.

## Conventions and caveats

* **Coordinates.** Everything the user types or reads is 1-based inclusive. BED
  columns and `start0` are 0-based.
* **Identity** = 1 - edit_distance / max(len_a, len_b), using unit-cost
  Levenshtein distance. It is not a gap-affine alignment identity.
* **Truth outside the SV benchmark BED is unreliable.** `in_sv_benchmark` and
  `summary.json:truth.*_frac_in_sv_benchmark` say where you are.
* **Named mode costs one pass over the contig's paths** (~40 s and ~670 MB of
  streamed GAF for chr20); budget accordingly in batches.
* **CIGAR traceback memory** is about 2 × len(hap) × len(ref) bits: 100 MB at
  20 kb × 20 kb. Past `--max-align-cells`, only the distance is reported.
* **Robustness.** Reverse-oriented walks are flipped into CHM13 orientation.
  Fragments are classified rather than dropped. Empty results (no calls, no
  truth records, one walk) produce empty tables. At contig ends the edge node is
  used as the anchor, with a warning. If the `--between` query leaks past the
  anchors, which would mean they do not separate the graph, the next anchors
  outward are tried.

