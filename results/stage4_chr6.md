# Stage 4e: held-out chr6, the chr20 winner unchanged

chr6 was never used to fit anything in this project. The pipeline is the one that won on chr20, with no
parameter changed: every tandem repeat whose hap32 haplotype lengths differ by >= 50 bp is realigned with
abPOA defaults on the 34 hap32 sequences, spliced into the chr6 hap32 graph, the reads production placed
on chr6 are re-mapped to the patched and the unpatched graph alike, and both are called with the pinned
PR-branch vg and the production flags (`work/stage4_chr6/run_chr6.sh`, `run_chr6b.sh`; scores in
`work/stage4_chr6/chr6/score/`). One difference from chr20: vg call reads the alignments from a per-arm
GAF-Base (`--gaf-reads` would need ~40 GB RAM on chr6); both chr6 arms use it.

Coverage: 46,699 merged TR regions, 1,120 eligible (671 VNTR, 398 STR, rest SAT/LC), 1,119 packaged,
1,102 realigned, 1,080 patched. 34.0M reads rebuilt and re-mapped per arm.

| chr6 | production | unpatched | patched (1,080 regions) |
|---|---|---|---|
| SV precision / recall | 0.5670 / 0.6063 | 0.5734 / 0.6063 | **0.6753 / 0.6128** |
| SV F1 | 0.5860 | 0.5894 | **0.6425** |
| SV FP / FN | 717 / 609 | 697 / 609 | **453 / 599** |
| SV FP inside / outside patched | 655 / 62 | 642 / 55 | 395 / 58 |
| SV FN inside / outside patched | 557 / 52 | 553 / 56 | 544 / 55 |
| SNV F1 / indel F1 | 0.9879 / 0.9396 | 0.9879 / 0.9390 | 0.9878 / 0.9386 |

The chr20 result holds out: SV F1 +0.053 (chr20 +0.051), false positives -244 (-35%; chr20 -32%),
false negatives slightly down (-10), SV calls outside the patched regions unchanged, small-variant F1
within 0.0004. Same caveats as chr20: hap32-only patch, contig-local re-mapping, the node-ID depth-rate
window in vg call (patched IDs spread over each region's original range), anchors chosen away from truth
records by package_regions.py.
