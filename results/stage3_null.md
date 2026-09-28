# Stage 3 null graphs: MC's own alignment written as a candidate, paired against the Stage 3 MC call (pilot, scorer version 2)

Written by `tools/score_haplotypes.py summarise` from `results/stage3/<label>/<id>.json`; the scoring rule is in the docstring of `tools/score_haplotypes.py`. Regions: pilot. Baseline: `mc`.

ED is the summed unit edit distance between the called and the true haplotype pair over the anchor-to-anchor span, best pairing, 0 = both haplotypes exact. Differences are label minus baseline on the regions both have, so negative is better. "better/worse >0" count regions whose ED fell/rose at all, ">10" by more than 10. Gain = 1 - ED/ED_ref, where ED_ref calls CHM13 on both haplotypes; "pooled" sums ED and ED_ref over the stratum. F1s are pooled over the stratum (TP/FP/FN summed) and are representation scores only.

- `mc`: 19 scored
- `mc_relabel`: 19 scored
- `mc_unchop`: 19 scored
- `genomewide`: 19 scored
- `truth`: 19 scored

| stratum | label | n | median ED (base) | sum ED (base) | median diff | better / worse >0 | better / worse >10 | exact (base) | median ED/kb | gain pooled / median | raw / refined / phab F1 | FP / FN (raw) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hotspot_vntr | mc | 9 | 2166.0 (2166.0) | 27140 (27140) | 0.0 | 0 / 0 | 0 / 0 | 0 (0) | 126.11 | 0.301 / 0.659 | 0.237 / 0.237 / 0.576 | 160 / 86 |
| hotspot_vntr | mc_relabel | 9 | 1355.0 (2166.0) | 25137 (27140) | 0.0 | 1 / 4 | 1 / 3 | 0 (0) | 102.60 | 0.353 / 0.737 | 0.249 / 0.249 / 0.546 | 152 / 85 |
| hotspot_vntr | mc_unchop | 9 | 1977.0 (2166.0) | 23089 (27140) | -584.0 | 5 / 3 | 5 / 2 | 0 (0) | 102.60 | 0.406 / 0.736 | 0.263 / 0.263 / 0.554 | 148 / 83 |
| hotspot_vntr | genomewide | 9 | 1977.0 (2166.0) | 25849 (27140) | 41.0 | 3 / 5 | 3 / 5 | 0 (0) | 102.60 | 0.335 / 0.736 | 0.236 / 0.236 / 0.502 | 157 / 86 |
| hotspot_vntr | truth | 9 | 1.0 (2166.0) | 941 (27140) | -2166.0 | 9 / 0 | 9 / 0 | 4 (0) | 0.03 | 0.976 / 1.000 | 0.309 / 0.629 / 0.933 | 61 / 89 |
| hotspot_other | mc | 1 | 0.0 (0.0) | 0 (0) | 0.0 | 0 / 0 | 0 / 0 | 1 (1) | 0.00 | 1.000 / 1.000 | 0.000 / 0.000 / 0.800 | 13 / 2 |
| hotspot_other | mc_relabel | 1 | 0.0 (0.0) | 0 (0) | 0.0 | 0 / 0 | 0 / 0 | 1 (1) | 0.00 | 1.000 / 1.000 | 0.000 / 0.000 / 0.800 | 13 / 2 |
| hotspot_other | mc_unchop | 1 | 0.0 (0.0) | 0 (0) | 0.0 | 0 / 0 | 0 / 0 | 1 (1) | 0.00 | 1.000 / 1.000 | 0.000 / 0.000 / 0.800 | 13 / 2 |
| hotspot_other | genomewide | 1 | 2.0 (0.0) | 2 (0) | 2.0 | 0 / 1 | 0 / 0 | 0 (1) | 0.65 | 1.000 / 1.000 | 0.000 / 0.000 / 0.800 | 13 / 2 |
| hotspot_other | truth | 1 | 0.0 (0.0) | 0 (0) | 0.0 | 0 / 0 | 0 / 0 | 1 (1) | 0.00 | 1.000 / 1.000 | 0.667 / 1.000 / 1.000 | 0 / 1 |
| control_vntr_matched | mc | 7 | 0.0 (0.0) | 2 (2) | 0.0 | 0 / 0 | 0 / 0 | 6 (6) | 0.00 | 0.999 / 1.000 | 1.000 / 1.000 / 1.000 | 0 / 0 |
| control_vntr_matched | mc_relabel | 7 | 0.0 (0.0) | 2 (2) | 0.0 | 0 / 0 | 0 / 0 | 6 (6) | 0.00 | 0.999 / 1.000 | 1.000 / 1.000 / 1.000 | 0 / 0 |
| control_vntr_matched | mc_unchop | 7 | 0.0 (0.0) | 2 (2) | 0.0 | 0 / 0 | 0 / 0 | 6 (6) | 0.00 | 0.999 / 1.000 | 1.000 / 1.000 / 1.000 | 0 / 0 |
| control_vntr_matched | genomewide | 7 | 0.0 (0.0) | 6 (2) | 0.0 | 0 / 1 | 0 / 0 | 6 (6) | 0.00 | 0.996 / 1.000 | 1.000 / 1.000 / 1.000 | 0 / 0 |
| control_vntr_matched | truth | 7 | 0.0 (0.0) | 0 (2) | 0.0 | 1 / 0 | 0 / 0 | 7 (6) | 0.00 | 1.000 / 1.000 | 0.842 / 1.000 / 1.000 | 1 / 2 |
| control_nontr_sv | mc | 2 | 0.5 (0.5) | 1 (1) | 0.0 | 0 / 0 | 0 / 0 | 1 (1) | 0.26 | 0.998 / 0.998 | 1.000 / 1.000 / 1.000 | 0 / 0 |
| control_nontr_sv | mc_relabel | 2 | 0.5 (0.5) | 1 (1) | 0.0 | 0 / 0 | 0 / 0 | 1 (1) | 0.26 | 0.998 / 0.998 | 1.000 / 1.000 / 1.000 | 0 / 0 |
| control_nontr_sv | mc_unchop | 2 | 0.5 (0.5) | 1 (1) | 0.0 | 0 / 0 | 0 / 0 | 1 (1) | 0.26 | 0.998 / 0.998 | 1.000 / 1.000 / 1.000 | 0 / 0 |
| control_nontr_sv | genomewide | 2 | 0.0 (0.5) | 0 (1) | -0.5 | 1 / 0 | 0 / 0 | 2 (1) | 0.00 | 1.000 / 1.000 | 1.000 / 1.000 / 1.000 | 0 / 0 |
| control_nontr_sv | truth | 2 | 155.0 (0.5) | 310 (1) | 154.5 | 1 / 1 | 0 / 1 | 1 (1) | 100.19 | 0.512 / 0.495 | 0.667 / 0.667 / 0.667 | 0 / 1 |

## Sensitivity to the overlap rule

The same tables under the raw overlap rule (an allele is skipped when its REF span shares a base with an allele already applied to that haplotype; bcftools consensus's and region.py's rule), which the scorer does not use because it drops real variants from the truth VCF itself (results/stage3_scorer_validation.tsv). "changed" counts regions whose ED the rule moves.

| stratum | label | n | sum ED (base), raw rule | better / worse >0, raw rule | changed |
|---|---|---|---|---|---|
| hotspot_vntr | mc | 9 | 27109 (27109) | 0 / 0 | 5 |
| hotspot_vntr | mc_relabel | 9 | 25113 (27109) | 1 / 4 | 5 |
| hotspot_vntr | mc_unchop | 9 | 23075 (27109) | 5 / 2 | 3 |
| hotspot_vntr | genomewide | 9 | 25844 (27109) | 2 / 5 | 5 |
| hotspot_vntr | truth | 9 | 904 (27109) | 9 / 0 | 2 |
| hotspot_other | mc | 1 | 0 (0) | 0 / 0 | 0 |
| hotspot_other | mc_relabel | 1 | 0 (0) | 0 / 0 | 0 |
| hotspot_other | mc_unchop | 1 | 0 (0) | 0 / 0 | 0 |
| hotspot_other | genomewide | 1 | 2 (0) | 0 / 1 | 0 |
| hotspot_other | truth | 1 | 0 (0) | 0 / 0 | 0 |
| control_vntr_matched | mc | 7 | 2 (2) | 0 / 0 | 0 |
| control_vntr_matched | mc_relabel | 7 | 2 (2) | 0 / 0 | 0 |
| control_vntr_matched | mc_unchop | 7 | 2 (2) | 0 / 0 | 0 |
| control_vntr_matched | genomewide | 7 | 6 (2) | 0 / 1 | 0 |
| control_vntr_matched | truth | 7 | 0 (2) | 1 / 0 | 0 |
| control_nontr_sv | mc | 2 | 1 (1) | 0 / 0 | 0 |
| control_nontr_sv | mc_relabel | 2 | 1 (1) | 0 / 0 | 0 |
| control_nontr_sv | mc_unchop | 2 | 1 (1) | 0 / 0 | 0 |
| control_nontr_sv | genomewide | 2 | 0 (1) | 1 / 0 | 0 |
| control_nontr_sv | truth | 2 | 310 (1) | 1 / 1 | 0 |

## Sensitivity to the nesting clause

The same tables when a parent's non-reference allele suppresses every allele nested in it on that slot (sensitivity.ed_suppress_nested), against the rule, under which a child applies where the parent's trimmed change does not cover it. "under parent" counts the child alleles the rule applies on a slot where a record they are nested in is non-reference (records.applied_under_nonref_parent); "reordered" counts regions where a vg child's POS lies left of its parent's, which scorer version 1 got wrong.

| stratum | label | n | sum ED (base), suppress | better / worse >0, suppress | changed | under parent | reordered |
|---|---|---|---|---|---|---|---|
| hotspot_vntr | mc | 9 | 25999 (25999) | 0 / 0 | 7 | 23 | 4 |
| hotspot_vntr | mc_relabel | 9 | 24800 (25999) | 1 / 4 | 6 | 24 | 4 |
| hotspot_vntr | mc_unchop | 9 | 21477 (25999) | 5 / 2 | 6 | 18 | 4 |
| hotspot_vntr | genomewide | 9 | 24683 (25999) | 2 / 5 | 7 | 26 | 4 |
| hotspot_vntr | truth | 9 | 904 (25999) | 9 / 0 | 2 | 2 | 0 |
| hotspot_other | mc | 1 | 0 (0) | 0 / 0 | 0 | 0 | 0 |
| hotspot_other | mc_relabel | 1 | 0 (0) | 0 / 0 | 0 | 0 | 0 |
| hotspot_other | mc_unchop | 1 | 0 (0) | 0 / 0 | 0 | 0 | 0 |
| hotspot_other | genomewide | 1 | 2 (0) | 0 / 1 | 0 | 0 | 0 |
| hotspot_other | truth | 1 | 0 (0) | 0 / 0 | 0 | 0 | 0 |
| control_vntr_matched | mc | 7 | 2 (2) | 0 / 0 | 0 | 0 | 0 |
| control_vntr_matched | mc_relabel | 7 | 2 (2) | 0 / 0 | 0 | 0 | 0 |
| control_vntr_matched | mc_unchop | 7 | 2 (2) | 0 / 0 | 0 | 0 | 0 |
| control_vntr_matched | genomewide | 7 | 6 (2) | 0 / 1 | 0 | 0 | 0 |
| control_vntr_matched | truth | 7 | 0 (2) | 1 / 0 | 0 | 0 | 0 |
| control_nontr_sv | mc | 2 | 1 (1) | 0 / 0 | 0 | 0 | 0 |
| control_nontr_sv | mc_relabel | 2 | 1 (1) | 0 / 0 | 0 | 0 | 0 |
| control_nontr_sv | mc_unchop | 2 | 1 (1) | 0 / 0 | 0 | 0 | 0 |
| control_nontr_sv | genomewide | 2 | 0 (1) | 1 / 0 | 0 | 0 | 0 |
| control_nontr_sv | truth | 2 | 310 (1) | 1 / 1 | 0 | 0 | 0 |

## Noise-aware counts

From `results/stage3_noise.tsv` (tools/replicate_noise.py): each graph's noise at a region is the largest pair ED between its Stage 3 call and one of its replicate calls (flank 50 kb; span node IDs re-laid; span nodes on mc's native IDs), computed without the truth. A region counts as better (worse) here only when the label's ED is below (above) the baseline's by more than the two graphs' noise summed. The rule's counts above are unchanged. "n" counts the regions where both graphs have a noise measurement.

| stratum | label | n | summed margin | better / worse beyond noise | better / worse >0 (rule) | no noise |
|---|---|---|---|---|---|---|
| hotspot_vntr | mc_relabel | 9 | 13683 | 0 / 0 | 1 / 4 | - |
| hotspot_vntr | mc_unchop | 9 | 19958 | 2 / 0 | 5 / 3 | - |
| hotspot_vntr | genomewide | 0 | 0 | 0 / 0 | 3 / 5 | L001909,L002013,L005990,L009656,L011138,L012184,L012272,L014297,L015415 |
| hotspot_vntr | truth | 9 | 9761 | 7 / 0 | 9 / 0 | - |
| hotspot_other | mc_relabel | 1 | 0 | 0 / 0 | 0 / 0 | - |
| hotspot_other | mc_unchop | 1 | 0 | 0 / 0 | 0 / 0 | - |
| hotspot_other | genomewide | 0 | 0 | 0 / 0 | 0 / 1 | L016870 |
| hotspot_other | truth | 1 | 0 | 0 / 0 | 0 / 0 | - |
| control_vntr_matched | mc_relabel | 7 | 0 | 0 / 0 | 0 / 0 | - |
| control_vntr_matched | mc_unchop | 7 | 0 | 0 / 0 | 0 / 0 | - |
| control_vntr_matched | genomewide | 0 | 0 | 0 / 0 | 0 / 1 | L004021,L007065,L007172,L007687,L011430,L014111,L016124 |
| control_vntr_matched | truth | 7 | 0 | 1 / 0 | 1 / 0 | - |
| control_nontr_sv | mc_relabel | 2 | 0 | 0 / 0 | 0 / 0 | - |
| control_nontr_sv | mc_unchop | 2 | 0 | 0 / 0 | 0 / 0 | - |
| control_nontr_sv | genomewide | 0 | 0 | 0 / 0 | 1 / 0 | L000506,L006556 |
| control_nontr_sv | truth | 2 | 0 | 1 / 1 | 1 / 1 | - |

## Mechanical reading of the decision rule

Per label: hotspot ED falls = summed ED over the hotspot VNTRs is lower than the baseline's and more regions are better than worse (at >0); a control stratum is worse = its summed ED is higher or more regions are worse than better (at >0). This is only the arithmetic; the decision is made on the tables above.

- `mc_relabel`: hotspot VNTR ED does not fall (sum 25137 vs 27140, better/worse 1/4); controls worse: none.
  - beyond noise: hotspot VNTRs better/worse 0/0 of 9; controls worse beyond noise: none.
- `mc_unchop`: hotspot VNTR ED falls (sum 23089 vs 27140, better/worse 5/3); controls worse: none.
  - beyond noise: hotspot VNTRs better/worse 2/0 of 9; controls worse beyond noise: none.
- `genomewide`: hotspot VNTR ED does not fall (sum 25849 vs 27140, better/worse 3/5); controls worse: control_vntr_matched.
- `truth`: hotspot VNTR ED falls (sum 941 vs 27140, better/worse 9/0); controls worse: control_nontr_sv.
  - beyond noise: hotspot VNTRs better/worse 7/0 of 9; controls worse beyond noise: control_nontr_sv.

## Per region (ED)

Each graph's replicate noise is in brackets.

| region | stratum | mc | mc_relabel | mc_unchop | genomewide | truth | ED_ref | truth bp |
|---|---|---|---|---|---|---|---|---|
| L001909 | hotspot_vntr | 7469 [4237] | 7503 [1376] | 6612 [1881] | 6373 [-] | 1 [0] | 7652 | 33326 |
| L002013 | hotspot_vntr | 5031 [54] | 5084 [11] | 5096 [3080] | 5079 [-] | 36 [0] | 2458 | 8354 |
| L005990 | hotspot_vntr | 6630 [0] | 6630 [0] | 4518 [0] | 6671 [-] | 0 [0] | 846 | 13712 |
| L009656 | hotspot_vntr | 2166 [1003] | 2166 [1003] | 3504 [808] | 2241 [-] | 0 [0] | 673 | 9598 |
| L011138 | hotspot_vntr | 2561 [2823] | 422 [1858] | 1977 [4189] | 1977 [-] | 88 [0] | 7502 | 20308 |
| L012184 | hotspot_vntr | 1030 [514] | 1030 [534] | 236 [752] | 1188 [-] | 763 [0] | 8317 | 17818 |
| L012272 | hotspot_vntr | 1311 [52] | 1355 [138] | 199 [43] | 1407 [-] | 53 [1034] | 7022 | 104194 |
| L014297 | hotspot_vntr | 252 [44] | 257 [36] | 257 [0] | 223 [-] | 0 [0] | 1754 | 4282 |
| L015415 | hotspot_vntr | 690 [0] | 690 [0] | 690 [478] | 690 [-] | 0 [0] | 2621 | 6725 |
| L016870 | hotspot_other | 0 [0] | 0 [0] | 0 [0] | 2 [-] | 0 [0] | 10162 | 3080 |
| L004021 | control_vntr_matched | 2 [0] | 2 [0] | 2 [0] | 6 [-] | 0 [0] | 222 | 15325 |
| L007065 | control_vntr_matched | 0 [0] | 0 [0] | 0 [0] | 0 [-] | 0 [0] | 73 | 8575 |
| L007172 | control_vntr_matched | 0 [0] | 0 [0] | 0 [0] | 0 [-] | 0 [0] | 338 | 6356 |
| L007687 | control_vntr_matched | 0 [0] | 0 [0] | 0 [0] | 0 [-] | 0 [0] | 175 | 5460 |
| L011430 | control_vntr_matched | 0 [0] | 0 [0] | 0 [0] | 0 [-] | 0 [0] | 276 | 6303 |
| L014111 | control_vntr_matched | 0 [0] | 0 [0] | 0 [0] | 0 [-] | 0 [0] | 303 | 8375 |
| L016124 | control_vntr_matched | 0 [0] | 0 [0] | 0 [0] | 0 [-] | 0 [0] | 299 | 7230 |
| L000506 | control_nontr_sv | 1 [0] | 1 [0] | 1 [0] | 0 [-] | 0 [0] | 328 | 1930 |
| L006556 | control_nontr_sv | 0 [0] | 0 [0] | 0 [0] | 0 [-] | 310 [0] | 307 | 1547 |
