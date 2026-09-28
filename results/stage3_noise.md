# Stage 3 replicate noise: how far each graph's own call moves under neutral changes (pilot)

Written by `tools/replicate_noise.py` (the method is in its docstring) into `results/stage3_noise.tsv`. For each graph and region, d is the pair edit distance between the span haplotypes of the Stage 3 call (hybrid 200 kb, call_local.py) and of a replicate call of the same graph: `flank50k` (50 kb genome-wide flank, reads re-fetched and re-mapped), `ids` (the Stage 3 hybrid graph with its span node IDs re-laid by half a spread step; for mc, its native IDs spread) and `mcids` (the span nodes on mc's native span IDs; not for mc or mc_relabel), each re-mapped. Noise = the largest. No truth is used, so candidate graphs are measured blind.

| graph | stratum | n | summed noise | median | max (region) | regions noise > 0 / > 10 | summed d flank50k / ids / mcids |
|---|---|---|---|---|---|---|---|
| mc | hotspot_vntr | 9 | 8727 | 54.0 | 4237 (L001909) | 7 / 7 | 8684 (9) / 5767 (9) / - |
| mc | hotspot_other | 1 | 0 | 0.0 | 0 (L016870) | 0 / 0 | 0 (1) / 0 (1) / - |
| mc | control_vntr_matched | 7 | 0 | 0.0 | 0 (L016124) | 0 / 0 | 0 (7) / 0 (7) / - |
| mc | control_nontr_sv | 2 | 0 | 0.0 | 0 (L000506) | 0 / 0 | 0 (2) / 0 (2) / - |
| unit_aware__all | hotspot_vntr | 9 | 2291 | 209.0 | 815 (L015415) | 9 / 9 | 2205 (9) / 0 (9) / 199 (8) |
| unit_aware__all | hotspot_other | 1 | 1238 | 1238.0 | 1238 (L016870) | 1 / 1 | 1238 (1) / 0 (1) / 0 (1) |
| unit_aware__all | control_vntr_matched | 7 | 76 | 0.0 | 76 (L011430) | 1 / 1 | 76 (7) / 0 (7) / 0 (4) |
| unit_aware__all | control_nontr_sv | 2 | 0 | 0.0 | 0 (L000506) | 0 / 0 | 0 (2) / 0 (2) / 0 (1) |
| mafft_linsi | hotspot_vntr | 8 | 936 | 2.5 | 803 (L005990) | 6 / 2 | 936 (8) / 0 (8) / 56 (8) |
| mafft_linsi | hotspot_other | 1 | 0 | 0.0 | 0 (L016870) | 0 / 0 | 0 (1) / 0 (1) / 0 (1) |
| mafft_linsi | control_vntr_matched | 7 | 2 | 0.0 | 2 (L004021) | 1 / 0 | 2 (7) / 0 (7) / 0 (7) |
| mafft_linsi | control_nontr_sv | 2 | 310 | 155.0 | 310 (L006556) | 1 / 1 | 310 (2) / 0 (2) / 0 (1) |
| unit_aware | hotspot_vntr | 9 | 2580 | 184.0 | 932 (L002013) | 7 / 7 | 2580 (9) / 0 (9) / 253 (9) |
| unit_aware | hotspot_other | 1 | 0 | 0.0 | 0 (L016870) | 0 / 0 | 0 (1) / 0 (1) / 0 (1) |
| unit_aware | control_vntr_matched | 7 | 242 | 0.0 | 183 (L007172) | 3 / 2 | 242 (7) / 0 (7) / 0 (5) |
| unit_aware | control_nontr_sv | 2 | 310 | 155.0 | 310 (L006556) | 1 / 1 | 310 (2) / 0 (2) / 0 (1) |
| truth | hotspot_vntr | 9 | 1034 | 0.0 | 1034 (L012272) | 1 / 1 | 1034 (9) / 0 (9) / 0 (9) |
| truth | hotspot_other | 1 | 0 | 0.0 | 0 (L016870) | 0 / 0 | 0 (1) / 0 (1) / 0 (1) |
| truth | control_vntr_matched | 7 | 0 | 0.0 | 0 (L016124) | 0 / 0 | 0 (7) / 0 (7) / 0 (7) |
| truth | control_nontr_sv | 2 | 0 | 0.0 | 0 (L000506) | 0 / 0 | 0 (2) / 0 (2) / 0 (1) |
| mc_relabel | hotspot_vntr | 9 | 4956 | 138.0 | 1858 (L011138) | 7 / 7 | 4956 (9) / - / - |
| mc_relabel | hotspot_other | 1 | 0 | 0.0 | 0 (L016870) | 0 / 0 | 0 (1) / - / - |
| mc_relabel | control_vntr_matched | 7 | 0 | 0.0 | 0 (L016124) | 0 / 0 | 0 (7) / - / - |
| mc_relabel | control_nontr_sv | 2 | 0 | 0.0 | 0 (L000506) | 0 / 0 | 0 (2) / - / - |
| mc_unchop | hotspot_vntr | 9 | 11231 | 752.0 | 4189 (L011138) | 7 / 7 | 8202 (9) / 0 (9) / 3081 (9) |
| mc_unchop | hotspot_other | 1 | 0 | 0.0 | 0 (L016870) | 0 / 0 | 0 (1) / 0 (1) / 0 (1) |
| mc_unchop | control_vntr_matched | 7 | 0 | 0.0 | 0 (L016124) | 0 / 0 | 0 (7) / 0 (7) / 0 (7) |
| mc_unchop | control_nontr_sv | 2 | 0 | 0.0 | 0 (L000506) | 0 / 0 | 0 (2) / 0 (2) / 0 (2) |

## Per region

d flank50k / d ids / d mcids per graph; "-" where the replicate does not exist (summed columns above: sum (regions)).

| region | stratum | mc | unit_aware__all | mafft_linsi | unit_aware | truth | mc_relabel | mc_unchop |
|---|---|---|---|---|---|---|---|---|
| L009656 | hotspot_vntr | 1003 / 0 / - | 53 / 0 / - | 3 / 0 / 0 | 26 / 0 / 0 | 0 / 0 / 0 | 1003 / - / - | 808 / 0 / 0 |
| L014297 | hotspot_vntr | 44 / 33 / - | 17 / 0 / 17 | 2 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 36 / - / - | 0 / 0 / 0 |
| L015415 | hotspot_vntr | 0 / 0 / - | 815 / 0 / 0 | 2 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / - / - | 478 / 0 / 0 |
| L005990 | hotspot_vntr | 0 / 0 / - | 256 / 0 / 0 | 803 / 0 / 0 | 412 / 0 / 0 | 0 / 0 / 0 | 0 / - / - | 0 / 0 / 0 |
| L012184 | hotspot_vntr | 514 / 0 / - | 209 / 0 / 0 | 0 / 0 / 0 | 758 / 0 / 0 | 0 / 0 / 0 | 534 / - / - | 752 / 0 / 0 |
| L002013 | hotspot_vntr | 11 / 54 / - | 31 / 0 / 117 | 118 / 0 / 56 | 932 / 0 / 159 | 0 / 0 / 0 | 11 / - / - | 51 / 0 / 3080 |
| L011138 | hotspot_vntr | 2823 / 2817 / - | 484 / 0 / 57 | 0 / 0 / 0 | 184 / 0 / 93 | 0 / 0 / 0 | 1858 / - / - | 4189 / 0 / 0 |
| L012272 | hotspot_vntr | 52 / 44 / - | 245 / 0 / 8 | - | 29 / 0 / 1 | 1034 / 0 / 0 | 138 / - / - | 43 / 0 / 1 |
| L001909 | hotspot_vntr | 4237 / 2819 / - | 95 / 0 / 0 | 8 / 0 / 0 | 239 / 0 / 0 | 0 / 0 / 0 | 1376 / - / - | 1881 / 0 / 0 |
| L016870 | hotspot_other | 0 / 0 / - | 1238 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / - / - | 0 / 0 / 0 |
| L016124 | control_vntr_matched | 0 / 0 / - | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / - / - | 0 / 0 / 0 |
| L007687 | control_vntr_matched | 0 / 0 / - | 0 / 0 / - | 0 / 0 / 0 | 57 / 0 / 0 | 0 / 0 / 0 | 0 / - / - | 0 / 0 / 0 |
| L007172 | control_vntr_matched | 0 / 0 / - | 0 / 0 / - | 0 / 0 / 0 | 183 / 0 / - | 0 / 0 / 0 | 0 / - / - | 0 / 0 / 0 |
| L011430 | control_vntr_matched | 0 / 0 / - | 76 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / - / - | 0 / 0 / 0 |
| L007065 | control_vntr_matched | 0 / 0 / - | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / - / - | 0 / 0 / 0 |
| L004021 | control_vntr_matched | 0 / 0 / - | 0 / 0 / 0 | 2 / 0 / 0 | 2 / 0 / 0 | 0 / 0 / 0 | 0 / - / - | 0 / 0 / 0 |
| L014111 | control_vntr_matched | 0 / 0 / - | 0 / 0 / - | 0 / 0 / 0 | 0 / 0 / - | 0 / 0 / 0 | 0 / - / - | 0 / 0 / 0 |
| L000506 | control_nontr_sv | 0 / 0 / - | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / - / - | 0 / 0 / 0 |
| L006556 | control_nontr_sv | 0 / 0 / - | 0 / 0 / - | 310 / 0 / - | 310 / 0 / - | 0 / 0 / - | 0 / - / - | 0 / 0 / 0 |

## Replicate ED to the truth (graphs whose Stage 3 scores are already known only)

| region | graph | ED Stage 3 call | ED flank50k | ED ids | ED mcids |
|---|---|---|---|---|---|
| L009656 | mc | 2166 | 2521 | 2166 | - |
| L014297 | mc | 252 | 226 | 257 | - |
| L015415 | mc | 690 | 690 | 690 | - |
| L005990 | mc | 6630 | 6630 | 6630 | - |
| L012184 | mc | 1030 | 588 | 1030 | - |
| L002013 | mc | 5031 | 5032 | 5084 | - |
| L011138 | mc | 2561 | 1617 | 422 | - |
| L012272 | mc | 1311 | 1353 | 1355 | - |
| L001909 | mc | 7469 | 3402 | 7503 | - |
| L004021 | mc | 2 | 2 | 2 | - |
| L000506 | mc | 1 | 1 | 1 | - |
| L012184 | truth | 763 | 763 | 763 | 763 |
| L002013 | truth | 36 | 36 | 36 | 36 |
| L011138 | truth | 88 | 88 | 88 | 88 |
| L012272 | truth | 53 | 1087 | 53 | 53 |
| L001909 | truth | 1 | 1 | 1 | 1 |
| L006556 | truth | 310 | 310 | 310 | - |
| L009656 | mc_relabel | 2166 | 2521 | - | - |
| L014297 | mc_relabel | 257 | 228 | - | - |
| L015415 | mc_relabel | 690 | 690 | - | - |
| L005990 | mc_relabel | 6630 | 6630 | - | - |
| L012184 | mc_relabel | 1030 | 569 | - | - |
| L002013 | mc_relabel | 5084 | 5084 | - | - |
| L011138 | mc_relabel | 422 | 1617 | - | - |
| L012272 | mc_relabel | 1355 | 1482 | - | - |
| L001909 | mc_relabel | 7503 | 7522 | - | - |
| L004021 | mc_relabel | 2 | 2 | - | - |
| L000506 | mc_relabel | 1 | 1 | - | - |
| L009656 | mc_unchop | 3504 | 3269 | 3504 | 3504 |
| L014297 | mc_unchop | 257 | 257 | 257 | 257 |
| L015415 | mc_unchop | 690 | 212 | 690 | 690 |
| L005990 | mc_unchop | 4518 | 4518 | 4518 | 4518 |
| L012184 | mc_unchop | 236 | 798 | 236 | 236 |
| L002013 | mc_unchop | 5096 | 5060 | 5096 | 5030 |
| L011138 | mc_unchop | 1977 | 2405 | 1977 | 1977 |
| L012272 | mc_unchop | 199 | 168 | 199 | 198 |
| L001909 | mc_unchop | 6612 | 5479 | 6612 | 6612 |
| L004021 | mc_unchop | 2 | 2 | 2 | 2 |
| L000506 | mc_unchop | 1 | 1 | 1 | 1 |

Rows where every ED is 0 are left out.
