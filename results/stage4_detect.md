# Stage 4, part 1: detect and realign (chr20)

Every rule below was fixed before any chr20 truth was read. The files are in `work/stage4/`, and the tools are `tools/scan_tr_regions.py`, `tools/stage4_rule.py` and `tools/stage4_detect.py`, plus the existing `package_regions.py`, `evaluate(_all).py` and `realign.py poa_abpoa`, run with `VNTR_REGIONS/CANDIDATES/RESULTS/WORK` pointing into `work/stage4/`.

## A. Scan (one pass over chr20.gbz, 20 s, 1.6 GB RSS)

The chr20 TR catalogue has 18,061 rows. Merging rows less than 100 bp apart gives 17,594 regions. A merged region takes its class from its longest member. Each region is padded by 100 bp and bracketed by the CHM13 nodes that hold its padded ends. 419 regions sit inside a single CHM13 node and so have no variation. The prefilter keeps regions whose max |hap length - CHM13 length| is at least 50 bp (`scan.tsv`, `prefilter.bed`).

| class | merged regions | prefilter |
|---|---|---|
| VNTR | 3,347 | 428 |
| STR | 12,167 | 186 |
| SAT | 136 | 41 |
| LC | 1,944 | 25 |
| total | 17,594 | **680** |

## B. Package (hap32 only) and Stage 0 of mc.gfa

- **Packaged:** 669 of the 680 (`work/stage4/regions/`), and `validate` passes. 11 failed: 8 had "fixed anchors do not separate the graph", 2 exceeded the 250 kb max span (a 3.2 Mb ALR region and a 0.95 Mb STR), and 1 had no separating boundary nodes.
- **Duplicate spans:** there are 657 distinct anchor spans, because 7 spans are shared by 2-4 packages each. TR762947/948/949/951, for example, are byte-identical.
- **Stage 0 (`--skip truth`):** done for all 669. Three distinct pericentromeric ALR regions (185-226 kb spans, TR762945, TR762947, TR763165) had only the size and redundancy metrics computed. The full all-pairs alignment was still running after more than 20 minutes. The rule does not use alignment, and abPOA cannot run on these regions anyway. The copies TR762948/949/951 reuse TR762947's JSON.

## C. Detection rule (fitted on regions/, chr20 excluded)

The training set is 145 regions: 61 hotspots (hotspot_vntr + hotspot_other) and 84 controls. The brief said 146, but `regions.tsv` holds 149 regions and 4 of them are on chr20.

| feature (results/mc) | AUC | best threshold | BA |
|---|---|---|---|
| redundancy.kmer_frac_extra | **0.990** | 0.3199 | **0.959** |
| redundancy.kmer_extra_per_kb | 0.989 | 406 | 0.948 |
| alignment.all_cost_over_opt | 0.874 | 1.1765 | 0.818 |
| best two-feature: cost/opt >= 3.19 OR kfe >= 0.32 | | | 0.967 (under the +0.02 needed to replace a one-feature rule) |

**The rule is kmer_frac_extra >= 0.3199.** On the training set it catches 56/61 hotspots (35/37 VNTR, 21/24 other) and flags 0/84 controls. The margin is thin: the highest control is 0.310 and the lowest detected hotspot is 0.330.

Applied unchanged to chr20, it detects **57 of the 680** (`work/stage4/hotspots.tsv`): 35 VNTR, 13 STR, 8 SAT and 1 LC. Six of the 57 are the giant pericentromeric ALR packages, 4 of them identical copies. All three chr20 hotspots of the original set (L009000, L009448, L009398) are detected.

## D. Detection against chr20 truth (production wgs-mm095, truvari SV FP/FN overlapping the padded region)

chr20 has 428 SV FPs in total, and 422 of them (98.6%) lie in a TR region.

| set | regions | TR SV FPs held | regions with >= 1 FP | TR SV FNs held |
|---|---|---|---|---|
| all merged TR | 17,594 | 422 (100%) | 161 (0.9%) | 319 |
| prefilter (>= 50 bp) | 680 | 409 (**96.9%**) | 153 (**22.5%**) | 313 (98.1%) |
| packaged + Stage 0 | 669 | 403 (95.5%) | 150 (22.4%) | 306 |
| **detected** | **57** | **199 (47.2%)** | **36 (63.2%)**; 36/51 (70.6%) without the 6 pericentromeric ALR | 159 (49.8%) |

- **AUC on chr20:** within the prefilter, kmer_frac_extra separates the regions with at least one FP from the rest with AUC 0.777. On training it was 0.99, where the controls were chosen to be clean.
- **Missed regions:** 21 undetected regions have at least 3 FPs each. Most are VNTRs of about 1.5-2 kb with kfe between 0.16 and 0.32.
- **FPs per class:** VNTR holds 354 FP-region hits in the prefilter and 179 in the detected set. STR has 39 and 12, SAT 14 and 8, LC 2 and 0.

## E. Realignment (poa_abpoa: abpoa -m 0 -r 1, longest first; default 900 s / 10 GB caps)

- **Realigned:** 50 of the 57 detected regions (87.7%), and they hold 194 of the detected set's 199 FPs. msa_graph asserted that every path spells its sequence, and evaluate found all 50 valid. The outputs are `work/stage4/candidates/poa_abpoa/<id>.{msa.fa,gfa,realign.json}`.
- **Failures:** 7 regions, all memouts predicted before the run. One is TR768070 (D20S16 satellite, 22 kb, predicted 11.1 GB, 5 FPs). The other six are the 185-226 kb ALR packages (0 FPs, predicted about 1 TB). These stay unpatched.
- **Stage 0 of the candidates** (`work/stage4/results/poa_abpoa/`, paired in `candidates_stage0.tsv`):

| metric | median MC | median poa_abpoa |
|---|---|---|
| all_cost_over_opt | 1.488 | **1.099** |
| kmer_frac_extra | 0.415 | 0.329 |
| nodes | 738 | 366 |

  cost/opt improved at 50/50 regions and got worse at none.

## Caveats

- **Anchor choice reads the truth:** package_regions.py picks anchors that no HG002 truth record touches, as it did in Stages 0-3. This shapes the region spans but none of the detection decisions.
- **Detection runs on the packaged span:** Stage 0 metrics come from the anchor-to-anchor span (pad 200 bp), while the truth overlap uses the scan interval (pad 100 bp). The packaged spans of neighbouring regions can overlap, and 7 spans are exact duplicates. Patching must de-duplicate by span.
- **The rule transfers imperfectly:** it was trained on hotspots against matched clean controls. On chr20's broader prefilter population, 37% of detected regions have no SV FP (29% without the ALR), and about half of the TR SV FPs lie outside the detected set.
- **Pericentromeric ALR regions:** the 3 distinct giant ALR regions (6 packages) have 0 bp inside the stvar benchmark, so their 0 FPs are uninformative.
- **Depth-window caveat for the next step:** vg call measures read rate over 4,096-node-ID windows, so any patched graph must keep node IDs spread over each region's original ID range.
