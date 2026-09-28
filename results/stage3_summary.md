# Stage 3 baseline: today's genome-wide production calls (wgs-mm095) scored over each span (scorer version 2)

Written by `tools/score_haplotypes.py summarise` from `results/stage3/<label>/<id>.json`; the scoring rule is in the docstring of `tools/score_haplotypes.py`. Regions: all. Baseline: `genomewide`.

ED is the summed unit edit distance between the called and the true haplotype pair over the anchor-to-anchor span, best pairing, 0 = both haplotypes exact. Differences are label minus baseline on the regions both have, so negative is better. "better/worse >0" count regions whose ED fell/rose at all, ">10" by more than 10. Gain = 1 - ED/ED_ref, where ED_ref calls CHM13 on both haplotypes; "pooled" sums ED and ED_ref over the stratum. F1s are pooled over the stratum (TP/FP/FN summed) and are representation scores only.

- `genomewide`: 149 scored

| stratum | label | n | median ED (base) | sum ED (base) | median diff | better / worse >0 | better / worse >10 | exact (base) | median ED/kb | gain pooled / median | raw / refined / phab F1 | FP / FN (raw) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hotspot_vntr | genomewide | 40 | 1729.0 (1729.0) | 124168 (124168) | 0.0 | 0 / 0 | 0 / 0 | 0 (0) | 186.46 | 0.363 / 0.409 | 0.245 / 0.248 / 0.557 | 642 / 389 |
| hotspot_other | genomewide | 24 | 644.0 (644.0) | 28376 (28376) | 0.0 | 0 / 0 | 0 / 0 | 0 (0) | 126.53 | 0.390 / 0.416 | 0.229 / 0.229 / 0.572 | 229 / 66 |
| control_vntr_matched | genomewide | 40 | 0.0 (0.0) | 254 (254) | 0.0 | 0 / 0 | 0 / 0 | 28 (28) | 0.00 | 0.988 / 1.000 | 0.986 / 0.986 / 0.986 | 1 / 0 |
| control_vntr_correct | genomewide | 25 | 0.0 (0.0) | 1356 (1356) | 0.0 | 0 / 0 | 0 / 0 | 15 (15) | 0.00 | 0.961 / 1.000 | 1.000 / 1.000 / 1.000 | 0 / 0 |
| control_nontr_sv | genomewide | 20 | 0.0 (0.0) | 2 (2) | 0.0 | 0 / 0 | 0 / 0 | 19 (19) | 0.00 | 1.000 / 1.000 | 1.000 / 1.000 / 1.000 | 0 / 0 |

## Sensitivity to the overlap rule

The same tables under the raw overlap rule (an allele is skipped when its REF span shares a base with an allele already applied to that haplotype; bcftools consensus's and region.py's rule), which the scorer does not use because it drops real variants from the truth VCF itself (results/stage3_scorer_validation.tsv). "changed" counts regions whose ED the rule moves.

| stratum | label | n | sum ED (base), raw rule | better / worse >0, raw rule | changed |
|---|---|---|---|---|---|
| hotspot_vntr | genomewide | 40 | 124459 (124459) | 0 / 0 | 20 |
| hotspot_other | genomewide | 24 | 28450 (28450) | 0 / 0 | 13 |
| control_vntr_matched | genomewide | 40 | 254 (254) | 0 / 0 | 0 |
| control_vntr_correct | genomewide | 25 | 1356 (1356) | 0 / 0 | 0 |
| control_nontr_sv | genomewide | 20 | 2 (2) | 0 / 0 | 0 |

## Sensitivity to the nesting clause

The same tables when a parent's non-reference allele suppresses every allele nested in it on that slot (sensitivity.ed_suppress_nested), against the rule, under which a child applies where the parent's trimmed change does not cover it. "under parent" counts the child alleles the rule applies on a slot where a record they are nested in is non-reference (records.applied_under_nonref_parent); "reordered" counts regions where a vg child's POS lies left of its parent's, which scorer version 1 got wrong.

| stratum | label | n | sum ED (base), suppress | better / worse >0, suppress | changed | under parent | reordered |
|---|---|---|---|---|---|---|---|
| hotspot_vntr | genomewide | 40 | 123529 (123529) | 0 / 0 | 24 | 87 | 11 |
| hotspot_other | genomewide | 24 | 28384 (28384) | 0 / 0 | 14 | 35 | 5 |
| control_vntr_matched | genomewide | 40 | 173 (173) | 0 / 0 | 4 | 15 | 0 |
| control_vntr_correct | genomewide | 25 | 1338 (1338) | 0 / 0 | 1 | 1 | 0 |
| control_nontr_sv | genomewide | 20 | 2 (2) | 0 / 0 | 0 | 0 | 0 |

## Mechanical reading of the decision rule

Per label: hotspot ED falls = summed ED over the hotspot VNTRs is lower than the baseline's and more regions are better than worse (at >0); a control stratum is worse = its summed ED is higher or more regions are worse than better (at >0). This is only the arithmetic; the decision is made on the tables above.


## Per region (ED)

| region | stratum | genomewide | ED_ref | truth bp |
|---|---|---|---|---|
| L000034 | hotspot_vntr | 692 | 1166 | 7369 |
| L000289 | hotspot_vntr | 6089 | 3498 | 8244 |
| L001084 | hotspot_vntr | 2716 | 2530 | 15532 |
| L001264 | hotspot_vntr | 10386 | 2755 | 11293 |
| L001909 | hotspot_vntr | 6373 | 7652 | 33326 |
| L002013 | hotspot_vntr | 5079 | 2458 | 8354 |
| L002096 | hotspot_vntr | 731 | 4219 | 9951 |
| L002143 | hotspot_vntr | 1708 | 865 | 14465 |
| L002182 | hotspot_vntr | 100 | 1748 | 3930 |
| L004145 | hotspot_vntr | 1750 | 4581 | 27206 |
| L005411 | hotspot_vntr | 2241 | 2444 | 10413 |
| L005481 | hotspot_vntr | 1584 | 647 | 5588 |
| L005990 | hotspot_vntr | 6671 | 846 | 13712 |
| L006952 | hotspot_vntr | 353 | 1957 | 2895 |
| L007009 | hotspot_vntr | 2305 | 13050 | 18182 |
| L007040 | hotspot_vntr | 1440 | 31963 | 6877 |
| L007081 | hotspot_vntr | 1641 | 4677 | 7663 |
| L008216 | hotspot_vntr | 271 | 460 | 5804 |
| L008864 | hotspot_vntr | 7142 | 3694 | 7474 |
| L008976 | hotspot_vntr | 11670 | 10187 | 21309 |
| L009000 | hotspot_vntr | 2060 | 19625 | 9197 |
| L009398 | hotspot_vntr | 868 | 1211 | 2695 |
| L009448 | hotspot_vntr | 5062 | 3218 | 13657 |
| L009656 | hotspot_vntr | 2241 | 673 | 9598 |
| L009658 | hotspot_vntr | 887 | 8714 | 13005 |
| L009803 | hotspot_vntr | 1443 | 4745 | 5511 |
| L010328 | hotspot_vntr | 5276 | 2504 | 8328 |
| L011138 | hotspot_vntr | 1977 | 7502 | 20308 |
| L011463 | hotspot_vntr | 1410 | 2411 | 9747 |
| L012184 | hotspot_vntr | 1188 | 8317 | 17818 |
| L012272 | hotspot_vntr | 1407 | 7022 | 104194 |
| L012294 | hotspot_vntr | 10715 | 1476 | 6717 |
| L013244 | hotspot_vntr | 4119 | 1558 | 10235 |
| L014160 | hotspot_vntr | 622 | 2603 | 13956 |
| L014210 | hotspot_vntr | 2270 | 2443 | 12494 |
| L014297 | hotspot_vntr | 223 | 1754 | 4282 |
| L014358 | hotspot_vntr | 670 | 1231 | 4069 |
| L015331 | hotspot_vntr | 764 | 1547 | 8058 |
| L015347 | hotspot_vntr | 9334 | 12274 | 30329 |
| L015415 | hotspot_vntr | 690 | 2621 | 6725 |
| L002940 | hotspot_other | 701 | 1063 | 3641 |
| L003695 | hotspot_other | 59 | 430 | 4037 |
| L003828 | hotspot_other | 1677 | 7178 | 3596 |
| L004848 | hotspot_other | 864 | 306 | 4063 |
| L004989 | hotspot_other | 105 | 1132 | 2468 |
| L005064 | hotspot_other | 798 | 296 | 6829 |
| L005261 | hotspot_other | 710 | 1272 | 2968 |
| L005618 | hotspot_other | 1180 | 2474 | 4751 |
| L007161 | hotspot_other | 25 | 717 | 6521 |
| L007184 | hotspot_other | 256 | 489 | 2101 |
| L007480 | hotspot_other | 1360 | 259 | 3632 |
| L007798 | hotspot_other | 9528 | 2691 | 8295 |
| L007844 | hotspot_other | 2795 | 3063 | 6721 |
| L008877 | hotspot_other | 587 | 4243 | 9583 |
| L010060 | hotspot_other | 423 | 187 | 4796 |
| L012256 | hotspot_other | 2604 | 2048 | 5910 |
| L012489 | hotspot_other | 37 | 384 | 5838 |
| L013064 | hotspot_other | 312 | 512 | 4382 |
| L013531 | hotspot_other | 446 | 1836 | 3399 |
| L015126 | hotspot_other | 1126 | 632 | 5672 |
| L015708 | hotspot_other | 66 | 68 | 5817 |
| L015847 | hotspot_other | 584 | 3302 | 5374 |
| L015904 | hotspot_other | 2131 | 1764 | 6382 |
| L016870 | hotspot_other | 2 | 10162 | 3080 |
| L000308 | control_vntr_matched | 0 | 2123 | 6119 |
| L000381 | control_vntr_matched | 0 | 548 | 14160 |
| L000709 | control_vntr_matched | 4 | 238 | 12287 |
| L001894 | control_vntr_matched | 0 | 447 | 12942 |
| L001899 | control_vntr_matched | 0 | 183 | 3735 |
| L002884 | control_vntr_matched | 2 | 679 | 8222 |
| L003826 | control_vntr_matched | 110 | 121 | 6992 |
| L003872 | control_vntr_matched | 2 | 290 | 5775 |
| L004021 | control_vntr_matched | 6 | 222 | 15325 |
| L007065 | control_vntr_matched | 0 | 73 | 8575 |
| L007172 | control_vntr_matched | 0 | 338 | 6356 |
| L007376 | control_vntr_matched | 1 | 682 | 38058 |
| L007687 | control_vntr_matched | 0 | 175 | 5460 |
| L007717 | control_vntr_matched | 2 | 702 | 8266 |
| L009770 | control_vntr_matched | 0 | 2108 | 7268 |
| L010350 | control_vntr_matched | 0 | 225 | 14971 |
| L011227 | control_vntr_matched | 0 | 170 | 7056 |
| L011430 | control_vntr_matched | 0 | 276 | 6303 |
| L012313 | control_vntr_matched | 1 | 212 | 6404 |
| L014111 | control_vntr_matched | 0 | 303 | 8375 |
| L014181 | control_vntr_matched | 0 | 400 | 7966 |
| L015223 | control_vntr_matched | 0 | 2682 | 8140 |
| L016124 | control_vntr_matched | 0 | 299 | 7230 |
| L017325 | control_vntr_matched | 6 | 6654 | 7036 |
| TR101140 | control_vntr_matched | 0 | 0 | 14800 |
| TR368009 | control_vntr_matched | 0 | 7 | 21928 |
| TR423938 | control_vntr_matched | 0 | 0 | 4512 |
| TR424743 | control_vntr_matched | 0 | 108 | 7232 |
| TR446002 | control_vntr_matched | 1 | 15 | 4620 |
| TR455654 | control_vntr_matched | 0 | 7 | 4807 |
| TR460913 | control_vntr_matched | 0 | 0 | 2626 |
| TR499275 | control_vntr_matched | 0 | 57 | 8648 |
| TR499276 | control_vntr_matched | 0 | 10 | 8966 |
| TR573259 | control_vntr_matched | 0 | 0 | 5838 |
| TR668316 | control_vntr_matched | 0 | 0 | 6126 |
| TR687433 | control_vntr_matched | 1 | 1 | 6704 |
| TR687489 | control_vntr_matched | 0 | 0 | 5462 |
| TR724750 | control_vntr_matched | 0 | 1 | 5956 |
| TR743754 | control_vntr_matched | 0 | 7 | 5846 |
| TR831034 | control_vntr_matched | 118 | 200 | 13876 |
| L000407 | control_vntr_correct | 0 | 95 | 5536 |
| L000698 | control_vntr_correct | 0 | 1834 | 8678 |
| L001491 | control_vntr_correct | 3 | 914 | 3336 |
| L001922 | control_vntr_correct | 54 | 1423 | 5699 |
| L003119 | control_vntr_correct | 1102 | 2204 | 3614 |
| L003672 | control_vntr_correct | 0 | 246 | 5066 |
| L007033 | control_vntr_correct | 0 | 3956 | 1116 |
| L008887 | control_vntr_correct | 0 | 723 | 5898 |
| L008902 | control_vntr_correct | 0 | 294 | 5808 |
| L009771 | control_vntr_correct | 0 | 230 | 2568 |
| L009885 | control_vntr_correct | 69 | 364 | 8022 |
| L010113 | control_vntr_correct | 0 | 107 | 5944 |
| L010508 | control_vntr_correct | 0 | 75 | 7126 |
| L011121 | control_vntr_correct | 0 | 692 | 5788 |
| L011224 | control_vntr_correct | 1 | 535 | 5672 |
| L013000 | control_vntr_correct | 11 | 1602 | 7024 |
| L013074 | control_vntr_correct | 0 | 9862 | 3826 |
| L014324 | control_vntr_correct | 0 | 916 | 5934 |
| L014426 | control_vntr_correct | 1 | 552 | 3225 |
| L015285 | control_vntr_correct | 0 | 205 | 5162 |
| L015743 | control_vntr_correct | 2 | 5098 | 4486 |
| L016062 | control_vntr_correct | 0 | 61 | 5201 |
| L016523 | control_vntr_correct | 0 | 142 | 5579 |
| L016830 | control_vntr_correct | 1 | 2221 | 7418 |
| L016884 | control_vntr_correct | 112 | 170 | 4024 |
| L000506 | control_nontr_sv | 0 | 328 | 1930 |
| L001527 | control_nontr_sv | 0 | 297 | 1977 |
| L002292 | control_nontr_sv | 0 | 1479 | 2893 |
| L003171 | control_nontr_sv | 0 | 668 | 1464 |
| L003883 | control_nontr_sv | 0 | 2886 | 4542 |
| L004620 | control_nontr_sv | 0 | 313 | 1632 |
| L005688 | control_nontr_sv | 0 | 580 | 1930 |
| L006556 | control_nontr_sv | 0 | 307 | 1547 |
| L007192 | control_nontr_sv | 0 | 240 | 1610 |
| L008631 | control_nontr_sv | 0 | 336 | 1741 |
| L009283 | control_nontr_sv | 0 | 135 | 1435 |
| L009600 | control_nontr_sv | 0 | 210 | 2234 |
| L010452 | control_nontr_sv | 0 | 122 | 1724 |
| L011931 | control_nontr_sv | 0 | 837 | 2765 |
| L012530 | control_nontr_sv | 0 | 123 | 1442 |
| L013399 | control_nontr_sv | 0 | 1562 | 2876 |
| L014899 | control_nontr_sv | 0 | 51 | 1950 |
| L016040 | control_nontr_sv | 0 | 65 | 1497 |
| L016274 | control_nontr_sv | 0 | 528 | 2012 |
| L017097 | control_nontr_sv | 2 | 1714 | 3238 |
