# Candidate graphs by stratum

Results: `/Users/benedictpaten/PycharmProjects/vntr-eval/results`. Baseline: `mc`. Methods: `mc`, `mafft_fftns2`, `mafft_fftns2__all`, `mafft_fftnsi__all`, `mafft_linsi`, `mafft_linsi__all`, `mafft_einsi`, `mafft_einsi__all`, `mafft_ginsi`, `mafft_ginsi__all`, `poa_abpoa`, `poa_abpoa__all`, `poa_spoa`, `poa_spoa__all`, `poa_abpoa_mc`, `poa_abpoa_mc__all`, `unit_aware`, `unit_aware__all`.

| method | valid regions | invalid |
|---|---|---|
| mc | 149 | 0 |
| mafft_fftns2 | 149 | 0 |
| mafft_fftns2__all | 140 | 0 |
| mafft_fftnsi__all | 105 | 0 |
| mafft_linsi | 148 | 0 |
| mafft_linsi__all | 77 | 0 |
| mafft_einsi | 146 | 0 |
| mafft_einsi__all | 65 | 0 |
| mafft_ginsi | 148 | 0 |
| mafft_ginsi__all | 65 | 0 |
| poa_abpoa | 146 | 0 |
| poa_abpoa__all | 138 | 0 |
| poa_spoa | 135 | 0 |
| poa_spoa__all | 123 | 0 |
| poa_abpoa_mc | 149 | 0 |
| poa_abpoa_mc__all | 148 | 0 |
| unit_aware | 149 | 0 |
| unit_aware__all | 149 | 0 |

Cells: median [q1-q3] over regions. Paired columns compare each method with `mc` on the regions both have: median difference (method - baseline), regions better/worse/equal (by the metric's direction; '.' = no preferred direction), Wilcoxon signed-rank two-sided p.

## hotspot_vntr (40 regions)

| method | n | nodes / CHM13 kb | cost/opt, all pairs | excess edits / kb | affine cost/opt, all | unaligned homology bp / kb | SV pieces / path (median) | indel bp / net length change | k-mer extra positions (frac) | truth edits to best graph path (h1+h2) | truth-by-graph truvari F1, raw | truth-by-graph truvari F1, refined | reads: redundant (cross-walk) only | reads: placements / max per path | regions with cost/opt <= 1.1 | pooled truvari F1 raw / refined |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc | 40 | 764 | 1.566 | 182 | 1.475 | 40.5 | 8 | 1.685 | 0.601 | 3.5 | 0.276 | 1 | 0.237 | 2.389 | 2/40 | 0.262 / 0.867 |
| mafft_fftns2 | 40 | 372 | 1.327 | 97.0 | 1.29 | 4.64 | 8 | 1.267 | 0.547 | 4 | 0.296 | 1 | 0.188 | 2.556 | 4/40 | 0.319 / 0.951 |
| mafft_fftns2__all | 31 | 391 | 1.42 | 120 | 1.369 | 9.19 | 9 | 1.376 | 0.598 | 3 | 0.333 | 1 | 0.223 | 2.699 | 0/31 | 0.333 / 0.951 |
| mafft_fftnsi__all | 9 | 272 | 1.173 | 63.6 | 1.176 | 5.17 | 3 | 1.274 | 0.554 | 1 | 0.364 | 1 | 0.235 | 3.397 | 2/9 | 0.389 / 1 |
| mafft_linsi | 39 | 355 | 1.111 | 36.3 | 1.11 | 1.58 | 4 | 1.087 | 0.419 | 9 | 0.308 | 1 | 0.155 | 1.728 | 18/39 | 0.337 / 0.964 |
| mafft_linsi__all | 4 | 162 | 1.035 | 13.3 | 1.04 | 2.595 | 1 | 1.061 | 0.314 | 4 | 0.567 | 1 | 0.155 | 1.827 | 3/4 | 0.526 / 1 |
| mafft_einsi | 37 | 265 | 1.182 | 58.8 | 1.149 | 2.94 | 6 | 1.139 | 0.405 | 11 | 0.4 | 1 | 0.11 | 1.679 | 9/37 | 0.382 / 0.964 |
| mafft_einsi__all | 1 | 68.7 | 1.045 | 16.1 | 1.044 | 3.89 | 1.5 | 1.023 | 0.189 | 6 | 0.571 | 1 | 0.272 | 2.061 | 1/1 | 0.571 / 1 |
| mafft_ginsi | 39 | 352 | 1.11 | 39.6 | 1.112 | 1.41 | 5 | 1.091 | 0.418 | 9 | 0.333 | 1 | 0.148 | 1.724 | 15/39 | 0.324 / 0.918 |
| mafft_ginsi__all | 1 | 100.0 | 1.033 | 11.6 | 1.031 | 2.04 | 1 | 1.023 | 0.28 | 9 | 0.8 | 1 | 0.106 | 2.03 | 1/1 | 0.8 / 1 |
| poa_abpoa | 37 | 300 | 1.083 | 26.0 | 1.1 | 0.78 | 5 | 1.04 | 0.479 | 12 | 0.32 | 1 | 0.164 | 2.019 | 25/37 | 0.327 / 0.983 |
| poa_abpoa__all | 30 | 436 | 1.13 | 43.4 | 1.238 | 1.995 | 4.5 | 1.08 | 0.557 | 2 | 0.265 | 1 | 0.288 | 2.655 | 11/30 | 0.3 / 0.777 |
| poa_spoa | 27 | 345 | 1.067 | 16.8 | 1.135 | 0.87 | 4 | 1.034 | 0.515 | 6 | 0.316 | 1 | 0.239 | 1.981 | 20/27 | 0.347 / 0.978 |
| poa_spoa__all | 20 | 386 | 1.115 | 29.1 | 1.299 | 2.13 | 2 | 1.075 | 0.577 | 1 | 0.227 | 1 | 0.397 | 2.966 | 8/20 | 0.241 / 0.881 |
| poa_abpoa_mc | 40 | 299 | 1.17 | 49.0 | 1.148 | 4.56 | 5 | 1.166 | 0.423 | 18.5 | 0.357 | 1 | 0.15 | 1.832 | 11/40 | 0.299 / 0.919 |
| poa_abpoa_mc__all | 39 | 360 | 1.275 | 87.0 | 1.317 | 11.6 | 7 | 1.296 | 0.549 | 6 | 0.263 | 1 | 0.226 | 2.329 | 6/39 | 0.284 / 0.932 |
| unit_aware | 40 | 381 | 1.048 | 15.2 | 1.093 | 0.355 | 5 | 1.02 | 0.518 | 3 | 0.238 | 1 | 0.211 | 2.085 | 33/40 | 0.26 / 0.941 |
| unit_aware__all | 40 | 454 | 1.059 | 19.4 | 1.138 | 0.6 | 5 | 1.026 | 0.542 | 2 | 0.21 | 1 | 0.221 | 2.49 | 31/40 | 0.265 / 0.946 |

| metric | method | n | median [q1-q3] | mean | paired n | median diff [q1-q3] | better / worse / equal | p |
|---|---|---|---|---|---|---|---|---|
| nodes / CHM13 kb | mc | 40 | 764 [446-1465] | 1033 | - | - | - | - |
| nodes / CHM13 kb | mafft_fftns2 | 40 | 372 [198-672] | 469 | 40 | -407 [-773--228] | 39 / 1 / 0 | 2e-10 |
| nodes / CHM13 kb | mafft_fftns2__all | 31 | 391 [226-840] | 553 | 31 | -334 [-550--145] | 30 / 1 / 0 | 1.8e-08 |
| nodes / CHM13 kb | mafft_fftnsi__all | 9 | 272 [160-295] | 311 | 9 | -171 [-237--62] | 9 / 0 / 0 | 0.0039 |
| nodes / CHM13 kb | mafft_linsi | 39 | 355 [159-652] | 442 | 39 | -410 [-858--290] | 39 / 0 / 0 | 3.6e-12 |
| nodes / CHM13 kb | mafft_linsi__all | 4 | 162 [134-329] | 301 | 4 | -190 [-237--141] | 4 / 0 / 0 | 0.12 |
| nodes / CHM13 kb | mafft_einsi | 37 | 265 [113-624] | 366 | 37 | -418 [-853--345] | 37 / 0 / 0 | 1.5e-11 |
| nodes / CHM13 kb | mafft_einsi__all | 1 | 68.7 [68.7-68.7] | 68.7 | 1 | -95.8 [-95.8--95.8] | 1 / 0 / 0 | 1 |
| nodes / CHM13 kb | mafft_ginsi | 39 | 352 [157-652] | 437 | 39 | -453 [-859--290] | 39 / 0 / 0 | 3.6e-12 |
| nodes / CHM13 kb | mafft_ginsi__all | 1 | 100.0 [100.0-100.0] | 100.0 | 1 | -64.6 [-64.6--64.6] | 1 / 0 / 0 | 1 |
| nodes / CHM13 kb | poa_abpoa | 37 | 300 [160-559] | 391 | 37 | -424 [-760--334] | 37 / 0 / 0 | 1.5e-11 |
| nodes / CHM13 kb | poa_abpoa__all | 30 | 436 [290-746] | 604 | 30 | -234 [-418--132] | 26 / 4 / 0 | 6.9e-07 |
| nodes / CHM13 kb | poa_spoa | 27 | 345 [188-629] | 452 | 27 | -364 [-665--254] | 27 / 0 / 0 | 1.5e-08 |
| nodes / CHM13 kb | poa_spoa__all | 20 | 386 [305-651] | 540 | 20 | -168 [-297--37.6] | 16 / 4 / 0 | 0.00039 |
| nodes / CHM13 kb | poa_abpoa_mc | 40 | 299 [139-488] | 389 | 40 | -419 [-844--330] | 40 / 0 / 0 | 1.8e-12 |
| nodes / CHM13 kb | poa_abpoa_mc__all | 39 | 360 [193-749] | 507 | 39 | -390 [-748--203] | 39 / 0 / 0 | 3.6e-12 |
| nodes / CHM13 kb | unit_aware | 40 | 381 [218-647] | 500 | 40 | -379 [-699--222] | 39 / 1 / 0 | 5.5e-12 |
| nodes / CHM13 kb | unit_aware__all | 40 | 454 [245-763] | 582 | 40 | -346 [-582--137] | 37 / 3 / 0 | 6.7e-10 |
| cost/opt, all pairs | mc | 40 | 1.566 [1.338-1.81] | 1.755 | - | - | - | - |
| cost/opt, all pairs | mafft_fftns2 | 40 | 1.327 [1.183-1.469] | 1.57 | 40 | -0.2 [-0.374--0.037] | 33 / 7 / 0 | 1.1e-05 |
| cost/opt, all pairs | mafft_fftns2__all | 31 | 1.42 [1.27-1.74] | 1.628 | 31 | -0.146 [-0.36-0.027] | 22 / 9 / 0 | 0.0044 |
| cost/opt, all pairs | mafft_fftnsi__all | 9 | 1.173 [1.151-1.546] | 1.474 | 9 | -0.052 [-0.28--0.02] | 7 / 2 / 0 | 0.16 |
| cost/opt, all pairs | mafft_linsi | 39 | 1.111 [1.06-1.156] | 1.149 | 39 | -0.401 [-0.641--0.296] | 39 / 0 / 0 | 3.6e-12 |
| cost/opt, all pairs | mafft_linsi__all | 4 | 1.035 [1.03-1.054] | 1.049 | 4 | -0.28 [-0.832--0.133] | 4 / 0 / 0 | 0.12 |
| cost/opt, all pairs | mafft_einsi | 37 | 1.182 [1.112-1.252] | 1.207 | 37 | -0.36 [-0.56--0.243] | 36 / 1 / 0 | 1.5e-10 |
| cost/opt, all pairs | mafft_einsi__all | 1 | 1.045 [1.045-1.045] | 1.045 | 1 | -0.386 [-0.386--0.386] | 1 / 0 / 0 | 1 |
| cost/opt, all pairs | mafft_ginsi | 39 | 1.11 [1.063-1.167] | 1.151 | 39 | -0.405 [-0.627--0.295] | 39 / 0 / 0 | 3.6e-12 |
| cost/opt, all pairs | mafft_ginsi__all | 1 | 1.033 [1.033-1.033] | 1.033 | 1 | -0.398 [-0.398--0.398] | 1 / 0 / 0 | 1 |
| cost/opt, all pairs | poa_abpoa | 37 | 1.083 [1.044-1.111] | 1.121 | 37 | -0.454 [-0.699--0.306] | 37 / 0 / 0 | 1.5e-11 |
| cost/opt, all pairs | poa_abpoa__all | 30 | 1.13 [1.083-1.332] | 1.232 | 30 | -0.363 [-0.691--0.218] | 30 / 0 / 0 | 1.9e-09 |
| cost/opt, all pairs | poa_spoa | 27 | 1.067 [1.04-1.098] | 1.099 | 27 | -0.518 [-0.78--0.298] | 27 / 0 / 0 | 1.5e-08 |
| cost/opt, all pairs | poa_spoa__all | 20 | 1.115 [1.078-1.34] | 1.179 | 20 | -0.415 [-0.914--0.224] | 19 / 1 / 0 | 3.8e-06 |
| cost/opt, all pairs | poa_abpoa_mc | 40 | 1.17 [1.082-1.215] | 1.2 | 40 | -0.391 [-0.603--0.253] | 39 / 1 / 0 | 6e-11 |
| cost/opt, all pairs | poa_abpoa_mc__all | 39 | 1.275 [1.141-1.571] | 1.464 | 39 | -0.186 [-0.415--0.004] | 29 / 10 / 0 | 4.4e-05 |
| cost/opt, all pairs | unit_aware | 40 | 1.048 [1.021-1.06] | 1.062 | 40 | -0.464 [-0.716--0.321] | 40 / 0 / 0 | 1.8e-12 |
| cost/opt, all pairs | unit_aware__all | 40 | 1.059 [1.033-1.098] | 1.079 | 40 | -0.454 [-0.711--0.312] | 40 / 0 / 0 | 1.8e-12 |
| excess edits / kb | mc | 40 | 182 [113-299] | 202 | - | - | - | - |
| excess edits / kb | mafft_fftns2 | 40 | 97.0 [73.0-164] | 122 | 40 | -73.8 [-156--18.0] | 33 / 7 / 0 | 8.2e-06 |
| excess edits / kb | mafft_fftns2__all | 31 | 120 [96.0-189] | 144 | 31 | -57.3 [-111-6.88] | 22 / 9 / 0 | 0.0034 |
| excess edits / kb | mafft_fftnsi__all | 9 | 63.6 [53.9-91.9] | 94.4 | 9 | -11.0 [-58.0--6.01] | 7 / 2 / 0 | 0.13 |
| excess edits / kb | mafft_linsi | 39 | 36.3 [24.5-51.3] | 42.6 | 39 | -160 [-239--77.9] | 39 / 0 / 0 | 3.6e-12 |
| excess edits / kb | mafft_linsi__all | 4 | 13.3 [12.0-16.8] | 15.5 | 4 | -116 [-181--75.2] | 4 / 0 / 0 | 0.12 |
| excess edits / kb | mafft_einsi | 37 | 58.8 [31.6-93.4] | 64.0 | 37 | -122 [-175--67.9] | 36 / 1 / 0 | 1.5e-10 |
| excess edits / kb | mafft_einsi__all | 1 | 16.1 [16.1-16.1] | 16.1 | 1 | -137 [-137--137] | 1 / 0 / 0 | 1 |
| excess edits / kb | mafft_ginsi | 39 | 39.6 [24.7-55.8] | 43.2 | 39 | -158 [-238--78.0] | 39 / 0 / 0 | 3.6e-12 |
| excess edits / kb | mafft_ginsi__all | 1 | 11.6 [11.6-11.6] | 11.6 | 1 | -142 [-142--142] | 1 / 0 / 0 | 1 |
| excess edits / kb | poa_abpoa | 37 | 26.0 [15.4-41.3] | 31.3 | 37 | -163 [-246--77.7] | 37 / 0 / 0 | 1.5e-11 |
| excess edits / kb | poa_abpoa__all | 30 | 43.4 [25.9-75.5] | 51.8 | 30 | -124 [-220--47.3] | 30 / 0 / 0 | 1.9e-09 |
| excess edits / kb | poa_spoa | 27 | 16.8 [10.4-33.0] | 24.0 | 27 | -155 [-263--86.8] | 27 / 0 / 0 | 1.5e-08 |
| excess edits / kb | poa_spoa__all | 20 | 29.1 [21.0-44.6] | 33.6 | 20 | -127 [-219--52.9] | 19 / 1 / 0 | 3.8e-06 |
| excess edits / kb | poa_abpoa_mc | 40 | 49.0 [24.8-90.2] | 56.8 | 40 | -139 [-205--62.2] | 39 / 1 / 0 | 5.5e-12 |
| excess edits / kb | poa_abpoa_mc__all | 39 | 87.0 [56.3-155] | 112 | 39 | -74.0 [-184--3.725] | 29 / 10 / 0 | 3.3e-06 |
| excess edits / kb | unit_aware | 40 | 15.2 [6.447-28.6] | 17.9 | 40 | -175 [-267--96.4] | 40 / 0 / 0 | 1.8e-12 |
| excess edits / kb | unit_aware__all | 40 | 19.4 [11.4-33.8] | 22.7 | 40 | -166 [-260--92.8] | 40 / 0 / 0 | 1.8e-12 |
| affine cost/opt, all | mc | 40 | 1.475 [1.347-1.688] | 1.656 | - | - | - | - |
| affine cost/opt, all | mafft_fftns2 | 40 | 1.29 [1.172-1.389] | 1.486 | 40 | -0.172 [-0.352--0.063] | 33 / 7 / 0 | 8.2e-06 |
| affine cost/opt, all | mafft_fftns2__all | 31 | 1.369 [1.282-1.68] | 1.561 | 31 | -0.156 [-0.309-0.011] | 22 / 9 / 0 | 0.0029 |
| affine cost/opt, all | mafft_fftnsi__all | 9 | 1.176 [1.164-1.616] | 1.476 | 9 | -0.06 [-0.293--0.014] | 7 / 2 / 0 | 0.16 |
| affine cost/opt, all | mafft_linsi | 39 | 1.11 [1.07-1.141] | 1.118 | 39 | -0.396 [-0.602--0.289] | 39 / 0 / 0 | 3.6e-12 |
| affine cost/opt, all | mafft_linsi__all | 4 | 1.04 [1.032-1.058] | 1.05 | 4 | -0.292 [-0.79--0.141] | 4 / 0 / 0 | 0.12 |
| affine cost/opt, all | mafft_einsi | 37 | 1.149 [1.077-1.225] | 1.157 | 37 | -0.348 [-0.536--0.257] | 36 / 1 / 0 | 2.9e-11 |
| affine cost/opt, all | mafft_einsi__all | 1 | 1.044 [1.044-1.044] | 1.044 | 1 | -0.414 [-0.414--0.414] | 1 / 0 / 0 | 1 |
| affine cost/opt, all | mafft_ginsi | 39 | 1.112 [1.07-1.152] | 1.119 | 39 | -0.393 [-0.588--0.289] | 39 / 0 / 0 | 3.6e-12 |
| affine cost/opt, all | mafft_ginsi__all | 1 | 1.031 [1.031-1.031] | 1.031 | 1 | -0.426 [-0.426--0.426] | 1 / 0 / 0 | 1 |
| affine cost/opt, all | poa_abpoa | 37 | 1.1 [1.069-1.135] | 1.124 | 37 | -0.408 [-0.659--0.259] | 37 / 0 / 0 | 1.5e-11 |
| affine cost/opt, all | poa_abpoa__all | 30 | 1.238 [1.146-1.49] | 1.347 | 30 | -0.218 [-0.517--0.067] | 26 / 4 / 0 | 4.4e-05 |
| affine cost/opt, all | poa_spoa | 27 | 1.135 [1.106-1.192] | 1.157 | 27 | -0.417 [-0.664--0.235] | 27 / 0 / 0 | 1.5e-08 |
| affine cost/opt, all | poa_spoa__all | 20 | 1.299 [1.183-1.498] | 1.346 | 20 | -0.23 [-0.651--0.049] | 16 / 4 / 0 | 0.0014 |
| affine cost/opt, all | poa_abpoa_mc | 40 | 1.148 [1.098-1.201] | 1.159 | 40 | -0.37 [-0.558--0.258] | 39 / 1 / 0 | 9.1e-12 |
| affine cost/opt, all | poa_abpoa_mc__all | 39 | 1.317 [1.148-1.492] | 1.417 | 39 | -0.169 [-0.383--0.022] | 31 / 8 / 0 | 5.5e-05 |
| affine cost/opt, all | unit_aware | 40 | 1.093 [1.065-1.141] | 1.114 | 40 | -0.387 [-0.572--0.259] | 40 / 0 / 0 | 1.8e-12 |
| affine cost/opt, all | unit_aware__all | 40 | 1.138 [1.092-1.179] | 1.169 | 40 | -0.355 [-0.524--0.175] | 38 / 2 / 0 | 1.6e-10 |
| unaligned homology bp / kb | mc | 40 | 40.5 [18.5-87.3] | 73.8 | - | - | - | - |
| unaligned homology bp / kb | mafft_fftns2 | 40 | 4.64 [2.388-7.685] | 5.808 | 40 | -33.9 [-82.8--12.5] | 37 / 3 / 0 | 3.1e-10 |
| unaligned homology bp / kb | mafft_fftns2__all | 31 | 9.19 [5.9-21.3] | 13.5 | 31 | -18.0 [-35.1--1.41] | 24 / 7 / 0 | 0.00029 |
| unaligned homology bp / kb | mafft_fftnsi__all | 9 | 5.17 [1.72-6.8] | 6.21 | 9 | -22.0 [-27.9--6.41] | 7 / 2 / 0 | 0.074 |
| unaligned homology bp / kb | mafft_linsi | 39 | 1.58 [0.735-2.165] | 1.637 | 39 | -39.0 [-74.8--14.9] | 39 / 0 / 0 | 3.6e-12 |
| unaligned homology bp / kb | mafft_linsi__all | 4 | 2.595 [1.782-3.38] | 2.567 | 4 | -16.3 [-23.5--7.03] | 3 / 1 / 0 | 0.25 |
| unaligned homology bp / kb | mafft_einsi | 37 | 2.94 [1.32-6.52] | 4.412 | 37 | -27.7 [-68.2--10.5] | 37 / 0 / 0 | 1.5e-11 |
| unaligned homology bp / kb | mafft_einsi__all | 1 | 3.89 [3.89-3.89] | 3.89 | 1 | -7.69 [-7.69--7.69] | 1 / 0 / 0 | 1 |
| unaligned homology bp / kb | mafft_ginsi | 39 | 1.41 [0.65-2.15] | 1.604 | 39 | -39.0 [-74.4--14.8] | 39 / 0 / 0 | 3.6e-12 |
| unaligned homology bp / kb | mafft_ginsi__all | 1 | 2.04 [2.04-2.04] | 2.04 | 1 | -9.54 [-9.54--9.54] | 1 / 0 / 0 | 1 |
| unaligned homology bp / kb | poa_abpoa | 37 | 0.78 [0.52-1.46] | 1.38 | 37 | -35.9 [-58.1--16.5] | 37 / 0 / 0 | 1.5e-11 |
| unaligned homology bp / kb | poa_abpoa__all | 30 | 1.995 [1.205-4.3] | 3.431 | 30 | -23.6 [-44.7--9.325] | 30 / 0 / 0 | 1.9e-09 |
| unaligned homology bp / kb | poa_spoa | 27 | 0.87 [0.54-2.51] | 1.891 | 27 | -22.3 [-43.2--10.6] | 27 / 0 / 0 | 1.5e-08 |
| unaligned homology bp / kb | poa_spoa__all | 20 | 2.13 [1.538-3.918] | 2.881 | 20 | -19.4 [-35.2--9.305] | 20 / 0 / 0 | 1.9e-06 |
| unaligned homology bp / kb | poa_abpoa_mc | 40 | 4.56 [0.92-12.6] | 23.3 | 40 | -30.5 [-64.9--8.547] | 35 / 5 / 0 | 5.3e-07 |
| unaligned homology bp / kb | poa_abpoa_mc__all | 39 | 11.6 [4.78-22.1] | 17.0 | 39 | -22.1 [-53.5--6.315] | 35 / 4 / 0 | 1.3e-07 |
| unaligned homology bp / kb | unit_aware | 40 | 0.355 [0.065-0.985] | 0.693 | 40 | -40.2 [-86.6--17.7] | 40 / 0 / 0 | 1.8e-12 |
| unaligned homology bp / kb | unit_aware__all | 40 | 0.6 [0.158-1.24] | 0.809 | 40 | -40.1 [-86.6--17.1] | 40 / 0 / 0 | 1.8e-12 |
| SV pieces / path (median) | mc | 40 | 8 [6-11] | 8.9 | - | - | - | - |
| SV pieces / path (median) | mafft_fftns2 | 40 | 8 [4-12] | 8.95 | 40 | -1 [-5-4] | 26 / 12 / 2 | 0.54 |
| SV pieces / path (median) | mafft_fftns2__all | 31 | 9 [5-11.5] | 8.629 | 31 | -1 [-4-2.5] | 16 / 12 / 3 | 0.46 |
| SV pieces / path (median) | mafft_fftnsi__all | 9 | 3 [2-5] | 4.333 | 9 | -4 [-5--2] | 8 / 1 / 0 | 0.027 |
| SV pieces / path (median) | mafft_linsi | 39 | 4 [2-5] | 4.615 | 39 | -5 [-6--2] | 35 / 2 / 2 | 2.4e-09 |
| SV pieces / path (median) | mafft_linsi__all | 4 | 1 [1-1.25] | 1.25 | 4 | -6.5 [-7--6] | 4 / 0 / 0 | 0.12 |
| SV pieces / path (median) | mafft_einsi | 37 | 6 [3-8] | 6.811 | 37 | -3 [-5-1] | 25 / 10 / 2 | 0.013 |
| SV pieces / path (median) | mafft_einsi__all | 1 | 1.5 [1.5-1.5] | 1.5 | 1 | -6.5 [-6.5--6.5] | 1 / 0 / 0 | 1 |
| SV pieces / path (median) | mafft_ginsi | 39 | 5 [2-6.5] | 5.167 | 39 | -4 [-6--2] | 34 / 1 / 4 | 2.1e-07 |
| SV pieces / path (median) | mafft_ginsi__all | 1 | 1 [1-1] | 1 | 1 | -7 [-7--7] | 1 / 0 / 0 | 1 |
| SV pieces / path (median) | poa_abpoa | 37 | 5 [3-7] | 6.311 | 37 | -4 [-6.5-0] | 27 / 7 / 3 | 0.00079 |
| SV pieces / path (median) | poa_abpoa__all | 30 | 4.5 [2.25-8] | 6.233 | 30 | -4 [-6--0.25] | 22 / 7 / 1 | 0.006 |
| SV pieces / path (median) | poa_spoa | 27 | 4 [1.75-6] | 4.722 | 27 | -4 [-7--2.5] | 24 / 2 / 1 | 4.6e-06 |
| SV pieces / path (median) | poa_spoa__all | 20 | 2 [1.375-5.25] | 3.725 | 20 | -6 [-7--2] | 18 / 2 / 0 | 9.2e-05 |
| SV pieces / path (median) | poa_abpoa_mc | 40 | 5 [3-7] | 5.1 | 40 | -4 [-5.25--1] | 33 / 4 / 3 | 3.4e-08 |
| SV pieces / path (median) | poa_abpoa_mc__all | 39 | 7 [4.5-9] | 7.282 | 39 | -2 [-5-0] | 27 / 9 / 3 | 0.016 |
| SV pieces / path (median) | unit_aware | 40 | 5 [2-16] | 10.3 | 40 | -1.5 [-5-2.25] | 22 / 15 / 3 | 0.38 |
| SV pieces / path (median) | unit_aware__all | 40 | 5 [1-12.5] | 9.85 | 40 | -2.5 [-7-1] | 26 / 12 / 2 | 0.24 |
| indel bp / net length change | mc | 40 | 1.685 [1.398-1.974] | 2.095 | - | - | - | - |
| indel bp / net length change | mafft_fftns2 | 40 | 1.267 [1.098-1.503] | 1.857 | 40 | -0.304 [-0.628--0.035] | 32 / 8 / 0 | 0.00099 |
| indel bp / net length change | mafft_fftns2__all | 31 | 1.376 [1.202-1.655] | 1.926 | 31 | -0.325 [-0.489--0.022] | 24 / 7 / 0 | 0.00076 |
| indel bp / net length change | mafft_fftnsi__all | 9 | 1.274 [1.247-1.97] | 1.592 | 9 | -0.379 [-1.877-0.006] | 6 / 3 / 0 | 0.074 |
| indel bp / net length change | mafft_linsi | 39 | 1.087 [1.024-1.139] | 1.143 | 39 | -0.608 [-0.823--0.384] | 39 / 0 / 0 | 3.6e-12 |
| indel bp / net length change | mafft_linsi__all | 4 | 1.061 [1.018-1.18] | 1.137 | 4 | -0.437 [-1.647--0.189] | 4 / 0 / 0 | 0.12 |
| indel bp / net length change | mafft_einsi | 37 | 1.139 [1.058-1.294] | 1.219 | 37 | -0.552 [-0.705--0.226] | 35 / 2 / 0 | 1e-10 |
| indel bp / net length change | mafft_einsi__all | 1 | 1.023 [1.023-1.023] | 1.023 | 1 | -0.645 [-0.645--0.645] | 1 / 0 / 0 | 1 |
| indel bp / net length change | mafft_ginsi | 39 | 1.091 [1.026-1.149] | 1.154 | 39 | -0.601 [-0.804--0.404] | 39 / 0 / 0 | 3.6e-12 |
| indel bp / net length change | mafft_ginsi__all | 1 | 1.023 [1.023-1.023] | 1.023 | 1 | -0.645 [-0.645--0.645] | 1 / 0 / 0 | 1 |
| indel bp / net length change | poa_abpoa | 37 | 1.04 [1.006-1.106] | 1.12 | 37 | -0.623 [-0.938--0.412] | 37 / 0 / 0 | 1.5e-11 |
| indel bp / net length change | poa_abpoa__all | 30 | 1.08 [1.016-1.26] | 1.239 | 30 | -0.596 [-0.926--0.375] | 30 / 0 / 0 | 1.9e-09 |
| indel bp / net length change | poa_spoa | 27 | 1.034 [1.005-1.104] | 1.069 | 27 | -0.664 [-1.1--0.564] | 27 / 0 / 0 | 1.5e-08 |
| indel bp / net length change | poa_spoa__all | 20 | 1.075 [1.006-1.307] | 1.158 | 20 | -0.655 [-1.861--0.503] | 20 / 0 / 0 | 1.9e-06 |
| indel bp / net length change | poa_abpoa_mc | 40 | 1.166 [1.034-1.285] | 1.232 | 40 | -0.504 [-0.761--0.25] | 39 / 1 / 0 | 5.5e-12 |
| indel bp / net length change | poa_abpoa_mc__all | 39 | 1.296 [1.045-1.702] | 1.659 | 39 | -0.208 [-0.536--0.061] | 35 / 4 / 0 | 1e-06 |
| indel bp / net length change | unit_aware | 40 | 1.02 [1.004-1.052] | 1.054 | 40 | -0.654 [-0.937--0.373] | 40 / 0 / 0 | 1.8e-12 |
| indel bp / net length change | unit_aware__all | 40 | 1.026 [1.006-1.088] | 1.072 | 40 | -0.631 [-0.923--0.398] | 40 / 0 / 0 | 1.8e-12 |
| k-mer extra positions (frac) | mc | 40 | 0.601 [0.431-0.709] | 0.578 | - | - | - | - |
| k-mer extra positions (frac) | mafft_fftns2 | 40 | 0.547 [0.416-0.659] | 0.53 | 40 | -0.043 [-0.103-0.006] | 29 / 11 / 0 | 0.0039 |
| k-mer extra positions (frac) | mafft_fftns2__all | 31 | 0.598 [0.482-0.715] | 0.573 | 31 | 0.005 [-0.055-0.075] | 15 / 16 / 0 | 0.44 |
| k-mer extra positions (frac) | mafft_fftnsi__all | 9 | 0.554 [0.418-0.583] | 0.527 | 9 | 0.025 [-0.053-0.125] | 4 / 5 / 0 | 0.5 |
| k-mer extra positions (frac) | mafft_linsi | 39 | 0.419 [0.305-0.577] | 0.453 | 39 | -0.123 [-0.169--0.074] | 39 / 0 / 0 | 3.6e-12 |
| k-mer extra positions (frac) | mafft_linsi__all | 4 | 0.314 [0.279-0.453] | 0.417 | 4 | -0.111 [-0.136--0.069] | 3 / 1 / 0 | 0.25 |
| k-mer extra positions (frac) | mafft_einsi | 37 | 0.405 [0.277-0.555] | 0.424 | 37 | -0.143 [-0.198--0.072] | 36 / 1 / 0 | 7.3e-11 |
| k-mer extra positions (frac) | mafft_einsi__all | 1 | 0.189 [0.189-0.189] | 0.189 | 1 | -0.22 [-0.22--0.22] | 1 / 0 / 0 | 1 |
| k-mer extra positions (frac) | mafft_ginsi | 39 | 0.418 [0.302-0.575] | 0.452 | 39 | -0.123 [-0.171--0.071] | 39 / 0 / 0 | 3.6e-12 |
| k-mer extra positions (frac) | mafft_ginsi__all | 1 | 0.28 [0.28-0.28] | 0.28 | 1 | -0.129 [-0.129--0.129] | 1 / 0 / 0 | 1 |
| k-mer extra positions (frac) | poa_abpoa | 37 | 0.479 [0.288-0.54] | 0.449 | 37 | -0.13 [-0.172--0.079] | 36 / 1 / 0 | 3.6e-10 |
| k-mer extra positions (frac) | poa_abpoa__all | 30 | 0.557 [0.443-0.681] | 0.563 | 30 | 0.012 [-0.017-0.068] | 14 / 16 / 0 | 0.42 |
| k-mer extra positions (frac) | poa_spoa | 27 | 0.515 [0.355-0.566] | 0.491 | 27 | -0.069 [-0.116--0.034] | 25 / 2 / 0 | 6.7e-06 |
| k-mer extra positions (frac) | poa_spoa__all | 20 | 0.577 [0.501-0.676] | 0.573 | 20 | 0.056 [0.006-0.098] | 4 / 16 / 0 | 0.044 |
| k-mer extra positions (frac) | poa_abpoa_mc | 40 | 0.423 [0.272-0.591] | 0.455 | 40 | -0.138 [-0.174--0.057] | 39 / 1 / 0 | 1.6e-08 |
| k-mer extra positions (frac) | poa_abpoa_mc__all | 39 | 0.549 [0.388-0.642] | 0.525 | 39 | -0.04 [-0.124-0.014] | 26 / 13 / 0 | 0.00044 |
| k-mer extra positions (frac) | unit_aware | 40 | 0.518 [0.351-0.592] | 0.479 | 40 | -0.078 [-0.155--0.032] | 37 / 3 / 0 | 9.2e-08 |
| k-mer extra positions (frac) | unit_aware__all | 40 | 0.542 [0.406-0.629] | 0.526 | 40 | -0.031 [-0.138--0.002] | 31 / 9 / 0 | 0.0051 |
| truth edits to best graph path (h1+h2) | mc | 40 | 3.5 [0.75-24] | 20.4 | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2 | 40 | 4 [0-19] | 20.6 | 40 | 0 [-4.25-4.5] | 16 / 12 / 12 | 1 |
| truth edits to best graph path (h1+h2) | mafft_fftns2__all | 31 | 3 [0-7.5] | 7.613 | 31 | 0 [-5.5-0] | 13 / 7 / 11 | 0.068 |
| truth edits to best graph path (h1+h2) | mafft_fftnsi__all | 9 | 1 [0-7] | 3.444 | 9 | 0 [-4-0] | 3 / 1 / 5 | 0.5 |
| truth edits to best graph path (h1+h2) | mafft_linsi | 39 | 9 [1.5-31] | 26.7 | 39 | 1 [0-11] | 8 / 20 / 11 | 0.065 |
| truth edits to best graph path (h1+h2) | mafft_linsi__all | 4 | 4 [0.75-7.5] | 4.25 | 4 | 3.5 [0-7.5] | 0 / 2 / 2 | 0.5 |
| truth edits to best graph path (h1+h2) | mafft_einsi | 37 | 11 [1-30] | 22.8 | 37 | 1 [0-10] | 7 / 19 / 11 | 0.044 |
| truth edits to best graph path (h1+h2) | mafft_einsi__all | 1 | 6 [6-6] | 6 | 1 | 6 [6-6] | 0 / 1 / 0 | 1 |
| truth edits to best graph path (h1+h2) | mafft_ginsi | 39 | 9 [1.5-32] | 27.7 | 39 | 2 [0-10] | 6 / 22 / 11 | 0.031 |
| truth edits to best graph path (h1+h2) | mafft_ginsi__all | 1 | 9 [9-9] | 9 | 1 | 9 [9-9] | 0 / 1 / 0 | 1 |
| truth edits to best graph path (h1+h2) | poa_abpoa | 37 | 12 [1-35] | 41.2 | 37 | 1 [0-26] | 7 / 19 / 11 | 0.0036 |
| truth edits to best graph path (h1+h2) | poa_abpoa__all | 30 | 2 [0-5] | 14.5 | 30 | 0 [-6.75-0] | 12 / 6 / 12 | 0.11 |
| truth edits to best graph path (h1+h2) | poa_spoa | 27 | 6 [1-15] | 26.4 | 27 | 0 [0-3.5] | 6 / 13 / 8 | 0.24 |
| truth edits to best graph path (h1+h2) | poa_spoa__all | 20 | 1 [0-3.25] | 5.6 | 20 | 0 [-8.25-0] | 9 / 3 / 8 | 0.072 |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc | 40 | 18.5 [1.75-38] | 46.1 | 40 | 5 [0-21.8] | 3 / 25 / 12 | 1.3e-05 |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc__all | 39 | 6 [1-17] | 25.2 | 39 | 0 [-2.5-3.5] | 15 / 13 / 11 | 0.93 |
| truth edits to best graph path (h1+h2) | unit_aware | 40 | 3 [0-9.75] | 26.8 | 40 | 0 [-6-2.25] | 16 / 12 / 12 | 0.79 |
| truth edits to best graph path (h1+h2) | unit_aware__all | 40 | 2 [0-6.25] | 16.8 | 40 | 0 [-6.25-0] | 18 / 7 / 15 | 0.063 |
| truth-by-graph truvari F1, raw | mc | 40 | 0.276 [0.157-0.364] | 0.264 | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftns2 | 40 | 0.296 [0.19-0.467] | 0.318 | 40 | 0.082 [-0.059-0.177] | 25 / 14 / 1 | 0.11 |
| truth-by-graph truvari F1, raw | mafft_fftns2__all | 31 | 0.333 [0.274-0.424] | 0.353 | 31 | 0.067 [-0.043-0.155] | 21 / 10 / 0 | 0.076 |
| truth-by-graph truvari F1, raw | mafft_fftnsi__all | 9 | 0.364 [0.286-0.4] | 0.408 | 9 | 0.059 [-0.048-0.165] | 6 / 3 / 0 | 0.25 |
| truth-by-graph truvari F1, raw | mafft_linsi | 39 | 0.308 [0.17-0.539] | 0.369 | 39 | 0.097 [-0.042-0.277] | 26 / 13 / 0 | 0.029 |
| truth-by-graph truvari F1, raw | mafft_linsi__all | 4 | 0.567 [0.329-0.85] | 0.612 | 4 | 0.259 [-0.015-0.54] | 3 / 1 / 0 | 0.38 |
| truth-by-graph truvari F1, raw | mafft_einsi | 37 | 0.4 [0.301-0.545] | 0.406 | 37 | 0.172 [0.032-0.238] | 29 / 7 / 1 | 1.2e-05 |
| truth-by-graph truvari F1, raw | mafft_einsi__all | 1 | 0.571 [0.571-0.571] | 0.571 | 1 | 0.264 [0.264-0.264] | 1 / 0 / 0 | 1 |
| truth-by-graph truvari F1, raw | mafft_ginsi | 39 | 0.333 [0.16-0.517] | 0.367 | 39 | 0.077 [-0.084-0.262] | 26 / 13 / 0 | 0.051 |
| truth-by-graph truvari F1, raw | mafft_ginsi__all | 1 | 0.8 [0.8-0.8] | 0.8 | 1 | 0.492 [0.492-0.492] | 1 / 0 / 0 | 1 |
| truth-by-graph truvari F1, raw | poa_abpoa | 37 | 0.32 [0.222-0.5] | 0.347 | 37 | 0.049 [-0.088-0.192] | 23 / 14 / 0 | 0.12 |
| truth-by-graph truvari F1, raw | poa_abpoa__all | 30 | 0.265 [0.2-0.41] | 0.27 | 30 | 0 [-0.087-0.103] | 14 / 14 / 2 | 0.97 |
| truth-by-graph truvari F1, raw | poa_spoa | 27 | 0.316 [0.186-0.493] | 0.35 | 27 | 0.039 [-0.051-0.148] | 16 / 10 / 1 | 0.27 |
| truth-by-graph truvari F1, raw | poa_spoa__all | 20 | 0.227 [0.116-0.331] | 0.218 | 20 | -0.081 [-0.235-0.072] | 6 / 12 / 2 | 0.14 |
| truth-by-graph truvari F1, raw | poa_abpoa_mc | 40 | 0.357 [0.161-0.451] | 0.341 | 40 | 0.063 [-0.111-0.195] | 26 / 13 / 1 | 0.09 |
| truth-by-graph truvari F1, raw | poa_abpoa_mc__all | 39 | 0.263 [0.178-0.406] | 0.282 | 39 | 0.005 [-0.08-0.148] | 20 / 18 / 1 | 0.42 |
| truth-by-graph truvari F1, raw | unit_aware | 40 | 0.238 [0.163-0.319] | 0.251 | 40 | -0.012 [-0.128-0.076] | 17 / 21 / 2 | 0.48 |
| truth-by-graph truvari F1, raw | unit_aware__all | 40 | 0.21 [0.117-0.375] | 0.242 | 40 | -0.023 [-0.14-0.099] | 17 / 21 / 2 | 0.48 |
| truth-by-graph truvari F1, refined | mc | 40 | 1 [1-1] | 0.89 | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftns2 | 40 | 1 [1-1] | 0.95 | 40 | 0 [0-0] | 6 / 4 / 30 | 0.32 |
| truth-by-graph truvari F1, refined | mafft_fftns2__all | 31 | 1 [1-1] | 0.952 | 31 | 0 [0-0] | 4 / 2 / 25 | 0.69 |
| truth-by-graph truvari F1, refined | mafft_fftnsi__all | 9 | 1 [1-1] | 1 | 9 | 0 [0-0] | 1 / 0 / 8 | 1 |
| truth-by-graph truvari F1, refined | mafft_linsi | 39 | 1 [1-1] | 0.957 | 39 | 0 [0-0] | 6 / 3 / 30 | 0.36 |
| truth-by-graph truvari F1, refined | mafft_linsi__all | 4 | 1 [1-1] | 1 | 4 | 0 [0-0] | 0 / 0 / 4 | - |
| truth-by-graph truvari F1, refined | mafft_einsi | 37 | 1 [1-1] | 0.948 | 37 | 0 [0-0] | 5 / 3 / 29 | 0.64 |
| truth-by-graph truvari F1, refined | mafft_einsi__all | 1 | 1 [1-1] | 1 | 1 | 0 [0-0] | 0 / 0 / 1 | - |
| truth-by-graph truvari F1, refined | mafft_ginsi | 39 | 1 [1-1] | 0.938 | 39 | 0 [0-0] | 6 / 4 / 29 | 0.49 |
| truth-by-graph truvari F1, refined | mafft_ginsi__all | 1 | 1 [1-1] | 1 | 1 | 0 [0-0] | 0 / 0 / 1 | - |
| truth-by-graph truvari F1, refined | poa_abpoa | 37 | 1 [1-1] | 0.961 | 37 | 0 [0-0] | 5 / 2 / 30 | 0.3 |
| truth-by-graph truvari F1, refined | poa_abpoa__all | 30 | 1 [0.561-1] | 0.802 | 30 | 0 [-0.027-0] | 3 / 8 / 19 | 0.1 |
| truth-by-graph truvari F1, refined | poa_spoa | 27 | 1 [1-1] | 0.942 | 27 | 0 [0-0] | 4 / 3 / 20 | 0.58 |
| truth-by-graph truvari F1, refined | poa_spoa__all | 20 | 1 [0.947-1] | 0.814 | 20 | 0 [-0.053-0] | 1 / 5 / 14 | 0.16 |
| truth-by-graph truvari F1, refined | poa_abpoa_mc | 40 | 1 [1-1] | 0.926 | 40 | 0 [0-0] | 6 / 4 / 30 | 0.56 |
| truth-by-graph truvari F1, refined | poa_abpoa_mc__all | 39 | 1 [1-1] | 0.928 | 39 | 0 [0-0] | 6 / 4 / 29 | 0.56 |
| truth-by-graph truvari F1, refined | unit_aware | 40 | 1 [1-1] | 0.837 | 40 | 0 [0-0] | 5 / 8 / 27 | 0.26 |
| truth-by-graph truvari F1, refined | unit_aware__all | 40 | 1 [1-1] | 0.894 | 40 | 0 [0-0] | 6 / 4 / 30 | 0.92 |
| reads: redundant (cross-walk) only | mc | 40 | 0.237 [0.183-0.339] | 0.264 | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftns2 | 40 | 0.188 [0.142-0.303] | 0.226 | 40 | -0.032 [-0.078-0.011] | 27 / 13 / 0 | 0.018 |
| reads: redundant (cross-walk) only | mafft_fftns2__all | 31 | 0.223 [0.144-0.368] | 0.254 | 31 | -0.025 [-0.053-0.02] | 21 / 10 / 0 | 0.13 |
| reads: redundant (cross-walk) only | mafft_fftnsi__all | 9 | 0.235 [0.166-0.261] | 0.251 | 9 | 0.015 [-0.056-0.034] | 4 / 5 / 0 | 0.82 |
| reads: redundant (cross-walk) only | mafft_linsi | 39 | 0.155 [0.098-0.195] | 0.158 | 39 | -0.081 [-0.16--0.049] | 37 / 2 / 0 | 3.9e-09 |
| reads: redundant (cross-walk) only | mafft_linsi__all | 4 | 0.155 [0.099-0.206] | 0.151 | 4 | -0.104 [-0.19--0.058] | 4 / 0 / 0 | 0.12 |
| reads: redundant (cross-walk) only | mafft_einsi | 37 | 0.11 [0.084-0.156] | 0.127 | 37 | -0.139 [-0.182--0.078] | 36 / 1 / 0 | 1.3e-09 |
| reads: redundant (cross-walk) only | mafft_einsi__all | 1 | 0.272 [0.272-0.272] | 0.272 | 1 | -0.183 [-0.183--0.183] | 1 / 0 / 0 | 1 |
| reads: redundant (cross-walk) only | mafft_ginsi | 39 | 0.148 [0.098-0.194] | 0.155 | 39 | -0.085 [-0.164--0.045] | 36 / 3 / 0 | 4.6e-09 |
| reads: redundant (cross-walk) only | mafft_ginsi__all | 1 | 0.106 [0.106-0.106] | 0.106 | 1 | -0.35 [-0.35--0.35] | 1 / 0 / 0 | 1 |
| reads: redundant (cross-walk) only | poa_abpoa | 37 | 0.164 [0.122-0.245] | 0.186 | 37 | -0.059 [-0.089--0.037] | 34 / 3 / 0 | 4.6e-08 |
| reads: redundant (cross-walk) only | poa_abpoa__all | 30 | 0.288 [0.17-0.409] | 0.296 | 30 | 0.027 [-0.014-0.056] | 9 / 21 / 0 | 0.033 |
| reads: redundant (cross-walk) only | poa_spoa | 27 | 0.239 [0.137-0.282] | 0.229 | 27 | -0.032 [-0.066--0.016] | 23 / 3 / 1 | 7.3e-06 |
| reads: redundant (cross-walk) only | poa_spoa__all | 20 | 0.397 [0.246-0.435] | 0.368 | 20 | 0.054 [0.035-0.1] | 4 / 16 / 0 | 0.00085 |
| reads: redundant (cross-walk) only | poa_abpoa_mc | 40 | 0.15 [0.112-0.23] | 0.176 | 40 | -0.068 [-0.116--0.023] | 38 / 2 / 0 | 3.1e-10 |
| reads: redundant (cross-walk) only | poa_abpoa_mc__all | 39 | 0.226 [0.144-0.301] | 0.237 | 39 | -0.016 [-0.047-0.007] | 24 / 14 / 1 | 0.03 |
| reads: redundant (cross-walk) only | unit_aware | 40 | 0.211 [0.118-0.26] | 0.211 | 40 | -0.046 [-0.082--0.011] | 32 / 8 / 0 | 4.2e-05 |
| reads: redundant (cross-walk) only | unit_aware__all | 40 | 0.221 [0.131-0.325] | 0.242 | 40 | -0.021 [-0.061-0.009] | 25 / 13 / 2 | 0.049 |
| reads: placements / max per path | mc | 40 | 2.389 [2.128-3.074] | 2.66 | - | - | - | - |
| reads: placements / max per path | mafft_fftns2 | 40 | 2.556 [1.865-3.402] | 2.804 | 40 | -0.136 [-0.425-0.357] | 25 / 15 / 0 | 0.68 |
| reads: placements / max per path | mafft_fftns2__all | 31 | 2.699 [2.054-3.787] | 3.01 | 31 | 0.076 [-0.19-0.853] | 14 / 17 / 0 | 0.073 |
| reads: placements / max per path | mafft_fftnsi__all | 9 | 3.397 [2.087-3.41] | 2.984 | 9 | 0.246 [-0.134-1.256] | 4 / 5 / 0 | 0.25 |
| reads: placements / max per path | mafft_linsi | 39 | 1.728 [1.501-2.396] | 2.033 | 39 | -0.531 [-0.837--0.282] | 38 / 1 / 0 | 7.3e-12 |
| reads: placements / max per path | mafft_linsi__all | 4 | 1.827 [1.725-1.94] | 1.838 | 4 | -0.373 [-0.647--0.254] | 4 / 0 / 0 | 0.12 |
| reads: placements / max per path | mafft_einsi | 37 | 1.679 [1.384-2.155] | 1.805 | 37 | -0.733 [-0.976--0.482] | 36 / 1 / 0 | 4.4e-11 |
| reads: placements / max per path | mafft_einsi__all | 1 | 2.061 [2.061-2.061] | 2.061 | 1 | -1.103 [-1.103--1.103] | 1 / 0 / 0 | 1 |
| reads: placements / max per path | mafft_ginsi | 39 | 1.724 [1.501-2.414] | 2.022 | 39 | -0.543 [-0.841--0.325] | 38 / 1 / 0 | 7.3e-12 |
| reads: placements / max per path | mafft_ginsi__all | 1 | 2.03 [2.03-2.03] | 2.03 | 1 | -1.134 [-1.134--1.134] | 1 / 0 / 0 | 1 |
| reads: placements / max per path | poa_abpoa | 37 | 2.019 [1.582-2.423] | 2.147 | 37 | -0.417 [-0.626--0.218] | 36 / 1 / 0 | 1e-10 |
| reads: placements / max per path | poa_abpoa__all | 30 | 2.655 [2.211-3.308] | 2.827 | 30 | 0.242 [0.082-0.597] | 6 / 24 / 0 | 0.0016 |
| reads: placements / max per path | poa_spoa | 27 | 1.981 [1.848-2.454] | 2.145 | 27 | -0.24 [-0.347--0.13] | 23 / 4 / 0 | 6.7e-06 |
| reads: placements / max per path | poa_spoa__all | 20 | 2.966 [2.493-3.264] | 3.064 | 20 | 0.794 [0.33-0.945] | 3 / 17 / 0 | 0.00017 |
| reads: placements / max per path | poa_abpoa_mc | 40 | 1.832 [1.572-2.479] | 2.35 | 40 | -0.508 [-0.786--0.347] | 37 / 3 / 0 | 2.4e-08 |
| reads: placements / max per path | poa_abpoa_mc__all | 39 | 2.329 [1.912-3] | 2.531 | 39 | -0.127 [-0.241-0.073] | 26 / 13 / 0 | 0.025 |
| reads: placements / max per path | unit_aware | 40 | 2.085 [1.798-2.925] | 2.398 | 40 | -0.238 [-0.398--0.059] | 31 / 8 / 1 | 5.3e-06 |
| reads: placements / max per path | unit_aware__all | 40 | 2.49 [1.98-3.121] | 2.671 | 40 | -0.003 [-0.277-0.258] | 20 / 20 / 0 | 0.9 |

| method | regions scored (raw / refined / phab) | raw TP-base / FN / TP-comp / FP | raw F1 | refined TP-base / FN / TP-comp / FP | refined F1 | phab TP-base / FN / TP-comp / FP | phab F1 |
|---|---|---|---|---|---|---|---|
| mc | 40 / 40 / 40 | 173 / 385 / 168 / 575 | 0.262 | 496 / 62 / 628 / 115 | 0.867 | 374 / 6 / 369 / 2 | 0.989 |
| mafft_fftns2 | 40 / 40 / 40 | 215 / 343 / 212 / 567 | 0.319 | 540 / 18 / 728 / 51 | 0.951 | 373 / 1 / 371 / 1 | 0.997 |
| mafft_fftns2__all | 31 / 31 / 31 | 153 / 205 / 151 / 402 | 0.333 | 342 / 16 / 523 / 30 | 0.951 | 230 / 1 / 231 / 1 | 0.996 |
| mafft_fftnsi__all | 9 / 9 / 9 | 23 / 16 / 20 / 49 | 0.389 | 39 / 0 / 69 / 0 | 1 | 27 / 0 / 25 / 0 | 1 |
| mafft_linsi | 39 / 39 / 39 | 164 / 379 / 159 / 259 | 0.337 | 517 / 26 / 408 / 10 | 0.964 | 357 / 0 / 355 / 0 | 1 |
| mafft_linsi__all | 4 / 4 / 4 | 10 / 11 / 10 / 7 | 0.526 | 21 / 0 / 17 / 0 | 1 | 17 / 0 / 17 / 0 | 1 |
| mafft_einsi | 37 / 37 / 37 | 190 / 285 / 189 / 328 | 0.382 | 452 / 23 / 505 / 12 | 0.964 | 327 / 0 / 325 / 1 | 0.998 |
| mafft_einsi__all | 1 / 1 / 1 | 2 / 1 / 2 / 2 | 0.571 | 3 / 0 / 4 / 0 | 1 | 3 / 0 / 3 / 0 | 1 |
| mafft_ginsi | 39 / 39 / 39 | 165 / 378 / 162 / 304 | 0.324 | 497 / 46 / 429 / 37 | 0.918 | 364 / 0 / 365 / 0 | 1 |
| mafft_ginsi__all | 1 / 1 / 1 | 2 / 1 / 2 / 0 | 0.8 | 3 / 0 / 2 / 0 | 1 | 4 / 0 / 4 / 0 | 1 |
| poa_abpoa | 37 / 37 / 37 | 155 / 326 / 150 / 303 | 0.327 | 474 / 7 / 444 / 9 | 0.983 | 328 / 0 / 331 / 0 | 1 |
| poa_abpoa__all | 30 / 30 / 30 | 115 / 224 / 111 / 303 | 0.3 | 274 / 65 / 310 / 104 | 0.777 | 222 / 3 / 225 / 2 | 0.989 |
| poa_spoa | 27 / 27 / 27 | 100 / 173 / 98 / 200 | 0.347 | 263 / 10 / 296 / 2 | 0.978 | 196 / 1 / 196 / 0 | 0.997 |
| poa_spoa__all | 20 / 20 / 20 | 38 / 122 / 39 / 120 | 0.241 | 139 / 21 / 142 / 17 | 0.881 | 109 / 0 / 110 / 0 | 1 |
| poa_abpoa_mc | 40 / 40 / 40 | 146 / 412 / 143 / 267 | 0.299 | 503 / 55 / 384 / 26 | 0.919 | 377 / 1 / 375 / 0 | 0.999 |
| poa_abpoa_mc__all | 39 / 39 / 39 | 158 / 385 / 152 / 397 | 0.284 | 496 / 47 / 522 / 27 | 0.932 | 367 / 0 / 369 / 0 | 1 |
| unit_aware | 40 / 40 / 40 | 204 / 354 / 201 / 794 | 0.26 | 514 / 44 / 956 / 39 | 0.941 | 387 / 2 / 384 / 1 | 0.996 |
| unit_aware__all | 40 / 40 / 40 | 203 / 355 / 200 / 761 | 0.265 | 533 / 25 / 901 / 60 | 0.946 | 365 / 0 / 371 / 0 | 1 |

## control_vntr_matched (40 regions)

| method | n | nodes / CHM13 kb | cost/opt, all pairs | excess edits / kb | affine cost/opt, all | unaligned homology bp / kb | SV pieces / path (median) | indel bp / net length change | k-mer extra positions (frac) | truth edits to best graph path (h1+h2) | truth-by-graph truvari F1, raw | truth-by-graph truvari F1, refined | reads: redundant (cross-walk) only | reads: placements / max per path | regions with cost/opt <= 1.1 | pooled truvari F1 raw / refined |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc | 40 | 101 | 1.102 | 4.405 | 1.037 | 0.145 | 1 | 1.113 | 0.053 | 0 | 1 | 1 | 0.054 | 1.071 | 19/40 | 0.958 / 0.958 |
| mafft_fftns2 | 40 | 38.3 | 1.1 | 8.47 | 1.057 | 0.02 | 1 | 1.106 | 0.052 | 0 | 0.833 | 1 | 0.069 | 1.094 | 20/40 | 0.584 / 0.832 |
| mafft_fftns2__all | 40 | 39.2 | 1.206 | 15.1 | 1.151 | 0.1 | 1 | 1.202 | 0.059 | 0 | 0.571 | 1 | 0.08 | 1.104 | 15/40 | 0.39 / 0.635 |
| mafft_fftnsi__all | 34 | 37.7 | 1.111 | 8.685 | 1.053 | 0.06 | 1 | 1.12 | 0.056 | 0 | 0.667 | 1 | 0.077 | 1.101 | 16/34 | 0.525 / 0.879 |
| mafft_linsi | 40 | 28.5 | 1.028 | 1.79 | 1.007 | 0 | 1 | 1.02 | 0.051 | 0 | 1 | 1 | 0.042 | 1.057 | 33/40 | 0.771 / 0.954 |
| mafft_linsi__all | 27 | 36.6 | 1.054 | 1.84 | 1.022 | 0 | 0 | 1.018 | 0.065 | 0 | 1 | 1 | 0.057 | 1.067 | 20/27 | 0.731 / 0.958 |
| mafft_einsi | 40 | 27.7 | 1.03 | 1.29 | 1.004 | 0 | 1 | 1.014 | 0.032 | 0 | 1 | 1 | 0.031 | 1.043 | 33/40 | 0.732 / 0.929 |
| mafft_einsi__all | 22 | 18.2 | 1.018 | 0.99 | 1.003 | 0 | 0 | 1.008 | 0.015 | 0 | 1 | 1 | 0.03 | 1.036 | 17/22 | 0.732 / 1 |
| mafft_ginsi | 40 | 28.5 | 1.027 | 1.79 | 1.008 | 0 | 1 | 1.02 | 0.051 | 0 | 1 | 1 | 0.042 | 1.061 | 34/40 | 0.722 / 0.942 |
| mafft_ginsi__all | 22 | 34.1 | 1.055 | 2.18 | 1.021 | 0 | 0 | 1.017 | 0.062 | 0 | 1 | 1 | 0.055 | 1.068 | 17/22 | 0.636 / 1 |
| poa_abpoa | 40 | 31.9 | 1.024 | 1.06 | 1.006 | 0 | 0.5 | 1.016 | 0.038 | 0 | 1 | 1 | 0.035 | 1.042 | 33/40 | 0.648 / 0.881 |
| poa_abpoa__all | 40 | 54.0 | 1.084 | 6.115 | 1.067 | 0.075 | 0 | 1.108 | 0.062 | 0 | 0.417 | 1 | 0.058 | 1.092 | 20/40 | 0.4 / 0.832 |
| poa_spoa | 39 | 44.5 | 1.016 | 1.21 | 1.018 | 0 | 1 | 1.016 | 0.058 | 0 | 0.667 | 1 | 0.057 | 1.075 | 34/39 | 0.405 / 0.746 |
| poa_spoa__all | 37 | 78.1 | 1.069 | 4.9 | 1.099 | 0.17 | 0 | 1.046 | 0.146 | 0 | 0 | 1 | 0.157 | 1.238 | 26/37 | 0.296 / 0.696 |
| poa_abpoa_mc | 40 | 31.7 | 1.031 | 1.65 | 1.008 | 0 | 1 | 1.031 | 0.032 | 0 | 1 | 1 | 0.035 | 1.044 | 31/40 | 0.811 / 0.958 |
| poa_abpoa_mc__all | 40 | 29.2 | 1.084 | 5.735 | 1.038 | 0.03 | 1 | 1.13 | 0.045 | 0 | 1 | 1 | 0.039 | 1.047 | 20/40 | 0.73 / 0.893 |
| unit_aware | 40 | 40.3 | 1.015 | 1.09 | 1.02 | 0 | 0 | 1.009 | 0.071 | 0 | 0.383 | 1 | 0.061 | 1.082 | 37/40 | 0.323 / 0.903 |
| unit_aware__all | 40 | 45.0 | 1.045 | 1.825 | 1.072 | 0.01 | 0 | 1.022 | 0.092 | 0 | 0.45 | 1 | 0.065 | 1.093 | 31/40 | 0.304 / 0.744 |

| metric | method | n | median [q1-q3] | mean | paired n | median diff [q1-q3] | better / worse / equal | p |
|---|---|---|---|---|---|---|---|---|
| nodes / CHM13 kb | mc | 40 | 101 [51.2-154] | 125 | - | - | - | - |
| nodes / CHM13 kb | mafft_fftns2 | 40 | 38.3 [16.3-55.1] | 54.5 | 40 | -58.8 [-87.4--30.3] | 39 / 1 / 0 | 3.6e-12 |
| nodes / CHM13 kb | mafft_fftns2__all | 40 | 39.2 [16.8-87.6] | 65.5 | 40 | -54.0 [-86.7--24.3] | 38 / 2 / 0 | 3.5e-08 |
| nodes / CHM13 kb | mafft_fftnsi__all | 34 | 37.7 [13.7-77.1] | 55.0 | 34 | -52.5 [-76.3--22.6] | 32 / 2 / 0 | 2.9e-09 |
| nodes / CHM13 kb | mafft_linsi | 40 | 28.5 [14.6-70.9] | 63.9 | 40 | -45.0 [-94.8--27.4] | 39 / 1 / 0 | 7.8e-09 |
| nodes / CHM13 kb | mafft_linsi__all | 27 | 36.6 [15.0-70.7] | 45.0 | 27 | -35.4 [-93.6--21.6] | 25 / 2 / 0 | 2e-06 |
| nodes / CHM13 kb | mafft_einsi | 40 | 27.7 [11.6-65.3] | 53.1 | 40 | -64.9 [-95.1--29.9] | 40 / 0 / 0 | 1.8e-12 |
| nodes / CHM13 kb | mafft_einsi__all | 22 | 18.2 [8.982-59.8] | 32.3 | 22 | -34.4 [-83.7--24.5] | 22 / 0 / 0 | 4.8e-07 |
| nodes / CHM13 kb | mafft_ginsi | 40 | 28.5 [15.6-71.1] | 65.4 | 40 | -42.6 [-94.8--24.9] | 39 / 1 / 0 | 1.2e-08 |
| nodes / CHM13 kb | mafft_ginsi__all | 22 | 34.1 [12.9-67.7] | 44.0 | 22 | -30.3 [-81.7--18.1] | 20 / 2 / 0 | 2.1e-05 |
| nodes / CHM13 kb | poa_abpoa | 40 | 31.9 [13.0-71.1] | 62.1 | 40 | -57.5 [-88.9--27.0] | 38 / 2 / 0 | 4.3e-09 |
| nodes / CHM13 kb | poa_abpoa__all | 40 | 54.0 [19.9-128] | 87.1 | 40 | -39.2 [-76.0--16.6] | 34 / 6 / 0 | 7.7e-05 |
| nodes / CHM13 kb | poa_spoa | 39 | 44.5 [16.4-93.0] | 70.5 | 39 | -40.4 [-90.9--22.1] | 37 / 2 / 0 | 2.6e-07 |
| nodes / CHM13 kb | poa_spoa__all | 37 | 78.1 [27.4-153] | 104 | 37 | -22.5 [-51.3--2.84] | 28 / 9 / 0 | 0.0013 |
| nodes / CHM13 kb | poa_abpoa_mc | 40 | 31.7 [14.0-61.3] | 48.4 | 40 | -65.6 [-99.8--28.8] | 40 / 0 / 0 | 1.8e-12 |
| nodes / CHM13 kb | poa_abpoa_mc__all | 40 | 29.2 [15.0-57.9] | 52.4 | 40 | -60.8 [-90.9--30.5] | 40 / 0 / 0 | 1.8e-12 |
| nodes / CHM13 kb | unit_aware | 40 | 40.3 [16.9-109] | 75.3 | 40 | -45.7 [-85.2--20.5] | 33 / 7 / 0 | 1.9e-06 |
| nodes / CHM13 kb | unit_aware__all | 40 | 45.0 [20.6-118] | 88.4 | 40 | -31.9 [-78.1--14.9] | 31 / 9 / 0 | 0.00013 |
| cost/opt, all pairs | mc | 40 | 1.102 [1.005-1.316] | 1.284 | - | - | - | - |
| cost/opt, all pairs | mafft_fftns2 | 40 | 1.1 [1.022-1.388] | 1.296 | 40 | 0 [-0.098-0.034] | 17 / 17 / 6 | 0.79 |
| cost/opt, all pairs | mafft_fftns2__all | 40 | 1.206 [1.028-1.809] | 2.864 | 40 | 0.007 [-0.059-0.439] | 16 / 22 / 2 | 0.079 |
| cost/opt, all pairs | mafft_fftnsi__all | 34 | 1.111 [1.024-1.439] | 2.219 | 34 | 0 [-0.103-0.02] | 15 / 16 / 3 | 0.87 |
| cost/opt, all pairs | mafft_linsi | 40 | 1.028 [1.002-1.066] | 1.075 | 40 | -0.018 [-0.254-0] | 24 / 9 / 7 | 2.5e-05 |
| cost/opt, all pairs | mafft_linsi__all | 27 | 1.054 [1.007-1.105] | 1.096 | 27 | -0.02 [-0.17-0] | 16 / 6 / 5 | 0.0042 |
| cost/opt, all pairs | mafft_einsi | 40 | 1.03 [1.002-1.061] | 1.067 | 40 | -0.019 [-0.254-0] | 24 / 6 / 10 | 3.9e-06 |
| cost/opt, all pairs | mafft_einsi__all | 22 | 1.018 [1-1.08] | 1.097 | 22 | -0.02 [-0.244-0] | 14 / 2 / 6 | 0.00015 |
| cost/opt, all pairs | mafft_ginsi | 40 | 1.027 [1.003-1.061] | 1.073 | 40 | -0.019 [-0.256-0] | 25 / 8 / 7 | 6.4e-06 |
| cost/opt, all pairs | mafft_ginsi__all | 22 | 1.055 [1.004-1.093] | 1.098 | 22 | -0.019 [-0.172-0] | 14 / 4 / 4 | 0.0056 |
| cost/opt, all pairs | poa_abpoa | 40 | 1.024 [1.002-1.067] | 1.077 | 40 | -0.019 [-0.15-0] | 29 / 3 / 8 | 6.4e-08 |
| cost/opt, all pairs | poa_abpoa__all | 40 | 1.084 [1.006-1.286] | 1.262 | 40 | 0 [-0.122-0.031] | 16 / 17 / 7 | 0.47 |
| cost/opt, all pairs | poa_spoa | 39 | 1.016 [1.003-1.047] | 1.047 | 39 | -0.057 [-0.273-0] | 27 / 6 / 6 | 2.6e-06 |
| cost/opt, all pairs | poa_spoa__all | 37 | 1.069 [1.015-1.124] | 1.103 | 37 | -0.012 [-0.206-0.011] | 20 / 11 / 6 | 0.0081 |
| cost/opt, all pairs | poa_abpoa_mc | 40 | 1.031 [1.005-1.067] | 1.093 | 40 | -0.005 [-0.145-0] | 25 / 6 / 9 | 1.5e-05 |
| cost/opt, all pairs | poa_abpoa_mc__all | 40 | 1.084 [1.007-1.476] | 1.304 | 40 | 0 [-0.053-0.031] | 15 / 14 / 11 | 0.86 |
| cost/opt, all pairs | unit_aware | 40 | 1.015 [1.003-1.033] | 1.033 | 40 | -0.065 [-0.286-0] | 28 / 7 / 5 | 1.1e-06 |
| cost/opt, all pairs | unit_aware__all | 40 | 1.045 [1.015-1.083] | 1.101 | 40 | -0.021 [-0.206-0.002] | 23 / 13 / 4 | 0.002 |
| excess edits / kb | mc | 40 | 4.405 [0.21-28.6] | 24.6 | - | - | - | - |
| excess edits / kb | mafft_fftns2 | 40 | 8.47 [0.49-26.3] | 22.4 | 40 | 0 [-3.895-2.983] | 17 / 17 / 6 | 0.95 |
| excess edits / kb | mafft_fftns2__all | 40 | 15.1 [0.902-63.4] | 55.0 | 40 | 0.085 [-2.958-42.5] | 15 / 22 / 3 | 0.1 |
| excess edits / kb | mafft_fftnsi__all | 34 | 8.685 [0.512-27.3] | 27.1 | 34 | 0 [-5.237-2.212] | 15 / 15 / 4 | 1 |
| excess edits / kb | mafft_linsi | 40 | 1.79 [0.087-6.345] | 9.855 | 40 | -1.055 [-13.4-0] | 24 / 9 / 7 | 7.5e-05 |
| excess edits / kb | mafft_linsi__all | 27 | 1.84 [0.44-9.485] | 12.3 | 27 | -0.96 [-6.66-0] | 17 / 6 / 4 | 0.0031 |
| excess edits / kb | mafft_einsi | 40 | 1.29 [0.05-5.442] | 8.732 | 40 | -1.01 [-10.2-0] | 24 / 6 / 10 | 4.4e-06 |
| excess edits / kb | mafft_einsi__all | 22 | 0.99 [0.003-6.792] | 12.5 | 22 | -0.965 [-5.982-0] | 14 / 2 / 6 | 0.00015 |
| excess edits / kb | mafft_ginsi | 40 | 1.79 [0.087-5.335] | 9.787 | 40 | -1.15 [-13.5-0] | 25 / 8 / 7 | 2.5e-05 |
| excess edits / kb | mafft_ginsi__all | 22 | 2.18 [0.158-8.938] | 14.0 | 22 | -0.905 [-5.975-0] | 14 / 4 / 4 | 0.0047 |
| excess edits / kb | poa_abpoa | 40 | 1.06 [0.05-5.125] | 9.863 | 40 | -1.58 [-10.1-0] | 29 / 3 / 8 | 8.2e-08 |
| excess edits / kb | poa_abpoa__all | 40 | 6.115 [0.318-27.7] | 19.2 | 40 | 0 [-3.767-2.185] | 16 / 17 / 7 | 0.71 |
| excess edits / kb | poa_spoa | 39 | 1.21 [0.055-4.04] | 3.989 | 39 | -2.5 [-24.9-0] | 27 / 6 / 6 | 9.3e-06 |
| excess edits / kb | poa_spoa__all | 37 | 4.9 [0.3-12.2] | 9.285 | 37 | -0.29 [-21.8-0.15] | 20 / 11 / 6 | 0.01 |
| excess edits / kb | poa_abpoa_mc | 40 | 1.65 [0.105-8.607] | 10.9 | 40 | -0.305 [-6.16-0] | 25 / 6 / 9 | 5.6e-05 |
| excess edits / kb | poa_abpoa_mc__all | 40 | 5.735 [0.235-30.8] | 25.2 | 40 | 0 [-3.91-2.305] | 16 / 14 / 10 | 1 |
| excess edits / kb | unit_aware | 40 | 1.09 [0.057-2.963] | 3.042 | 40 | -3.21 [-21.4-0] | 28 / 7 / 5 | 2.4e-06 |
| excess edits / kb | unit_aware__all | 40 | 1.825 [0.472-7.433] | 9.044 | 40 | -0.495 [-15.6-0.045] | 23 / 13 / 4 | 0.0022 |
| affine cost/opt, all | mc | 40 | 1.037 [1.004-1.231] | 1.204 | - | - | - | - |
| affine cost/opt, all | mafft_fftns2 | 40 | 1.057 [1.003-1.274] | 1.233 | 40 | 0 [-0.072-0.062] | 19 / 16 / 5 | 0.85 |
| affine cost/opt, all | mafft_fftns2__all | 40 | 1.151 [1.026-1.69] | 2.398 | 40 | 0.004 [-0.037-0.395] | 16 / 22 / 2 | 0.071 |
| affine cost/opt, all | mafft_fftnsi__all | 34 | 1.053 [1.013-1.411] | 1.735 | 34 | 0 [-0.076-0.035] | 15 / 16 / 3 | 0.9 |
| affine cost/opt, all | mafft_linsi | 40 | 1.007 [1.001-1.045] | 1.039 | 40 | -0.021 [-0.197-0] | 28 / 6 / 6 | 7.2e-06 |
| affine cost/opt, all | mafft_linsi__all | 27 | 1.022 [1.004-1.07] | 1.076 | 27 | -0.013 [-0.142-0] | 18 / 5 / 4 | 0.0031 |
| affine cost/opt, all | mafft_einsi | 40 | 1.004 [1-1.035] | 1.03 | 40 | -0.021 [-0.197-0] | 29 / 2 / 9 | 9.3e-09 |
| affine cost/opt, all | mafft_einsi__all | 22 | 1.003 [1-1.023] | 1.048 | 22 | -0.015 [-0.191-0] | 14 / 2 / 6 | 0.00021 |
| affine cost/opt, all | mafft_ginsi | 40 | 1.008 [1.001-1.041] | 1.04 | 40 | -0.021 [-0.197-0] | 27 / 7 / 6 | 2.9e-05 |
| affine cost/opt, all | mafft_ginsi__all | 22 | 1.021 [1.001-1.079] | 1.083 | 22 | -0.011 [-0.191-0] | 13 / 5 / 4 | 0.038 |
| affine cost/opt, all | poa_abpoa | 40 | 1.006 [1-1.051] | 1.065 | 40 | -0.015 [-0.126-0] | 26 / 5 / 9 | 0.00041 |
| affine cost/opt, all | poa_abpoa__all | 40 | 1.067 [1.008-1.366] | 1.358 | 40 | 0 [-0.032-0.098] | 14 / 19 / 7 | 0.47 |
| affine cost/opt, all | poa_spoa | 39 | 1.018 [1.001-1.063] | 1.058 | 39 | -0.029 [-0.147-0] | 24 / 9 / 6 | 0.00048 |
| affine cost/opt, all | poa_spoa__all | 37 | 1.099 [1.03-1.255] | 1.184 | 37 | 0.026 [-0.023-0.112] | 12 / 20 / 5 | 0.43 |
| affine cost/opt, all | poa_abpoa_mc | 40 | 1.008 [1.001-1.038] | 1.046 | 40 | -0.017 [-0.142-0] | 29 / 2 / 9 | 1.2e-06 |
| affine cost/opt, all | poa_abpoa_mc__all | 40 | 1.038 [1.004-1.36] | 1.232 | 40 | 0 [-0.041-0.043] | 17 / 14 / 9 | 0.99 |
| affine cost/opt, all | unit_aware | 40 | 1.02 [1.002-1.071] | 1.061 | 40 | -0.01 [-0.16-0.006] | 23 / 13 / 4 | 0.009 |
| affine cost/opt, all | unit_aware__all | 40 | 1.072 [1.017-1.158] | 1.149 | 40 | 0 [-0.13-0.033] | 16 / 20 / 4 | 0.86 |
| unaligned homology bp / kb | mc | 40 | 0.145 [0-3.473] | 8.982 | - | - | - | - |
| unaligned homology bp / kb | mafft_fftns2 | 40 | 0.02 [0-0.258] | 1.486 | 40 | -0.035 [-2.827-0] | 22 / 4 / 14 | 3.1e-05 |
| unaligned homology bp / kb | mafft_fftns2__all | 40 | 0.1 [0-0.772] | 1.836 | 40 | 0 [-1.595-0.01] | 18 / 13 / 9 | 0.035 |
| unaligned homology bp / kb | mafft_fftnsi__all | 34 | 0.06 [0-0.78] | 2.216 | 34 | 0 [-0.935-0.338] | 16 / 10 / 8 | 0.48 |
| unaligned homology bp / kb | mafft_linsi | 40 | 0 [0-0.12] | 0.219 | 40 | -0.01 [-3.473-0] | 20 / 3 / 17 | 6e-05 |
| unaligned homology bp / kb | mafft_linsi__all | 27 | 0 [0-0.415] | 0.301 | 27 | 0 [-0.685-0.01] | 10 / 7 / 10 | 0.15 |
| unaligned homology bp / kb | mafft_einsi | 40 | 0 [0-0.542] | 0.804 | 40 | -0.01 [-3.473-0] | 20 / 4 / 16 | 0.00028 |
| unaligned homology bp / kb | mafft_einsi__all | 22 | 0 [0-0.145] | 0.638 | 22 | 0 [-0.317-0] | 7 / 5 / 10 | 0.34 |
| unaligned homology bp / kb | mafft_ginsi | 40 | 0 [0-0.212] | 0.235 | 40 | -0.01 [-3.473-0] | 20 / 3 / 17 | 0.00013 |
| unaligned homology bp / kb | mafft_ginsi__all | 22 | 0 [0-0.103] | 0.2 | 22 | 0 [-0.287-0] | 7 / 5 / 10 | 0.18 |
| unaligned homology bp / kb | poa_abpoa | 40 | 0 [0-0.06] | 0.09 | 40 | -0.085 [-3.427-0] | 23 / 2 / 15 | 2e-06 |
| unaligned homology bp / kb | poa_abpoa__all | 40 | 0.075 [0-0.492] | 0.596 | 40 | 0 [-3.337-0] | 19 / 8 / 13 | 0.0015 |
| unaligned homology bp / kb | poa_spoa | 39 | 0 [0-0.02] | 0.106 | 39 | -0.15 [-3.32-0] | 23 / 0 / 16 | 2.4e-07 |
| unaligned homology bp / kb | poa_spoa__all | 37 | 0.17 [0-0.6] | 0.568 | 37 | 0 [-2.25-0.02] | 15 / 13 / 9 | 0.044 |
| unaligned homology bp / kb | poa_abpoa_mc | 40 | 0 [0-0.118] | 0.446 | 40 | -0.05 [-2.857-0] | 21 / 6 / 13 | 5.5e-05 |
| unaligned homology bp / kb | poa_abpoa_mc__all | 40 | 0.03 [0-1.065] | 3.287 | 40 | 0 [-1.27-0.035] | 14 / 14 / 12 | 0.25 |
| unaligned homology bp / kb | unit_aware | 40 | 0 [0-0.022] | 0.154 | 40 | -0.01 [-3.46-0] | 20 / 5 / 15 | 0.00014 |
| unaligned homology bp / kb | unit_aware__all | 40 | 0.01 [0-0.215] | 0.266 | 40 | -0.01 [-3.388-0] | 21 / 7 / 12 | 0.0011 |
| SV pieces / path (median) | mc | 40 | 1 [0-1] | 0.7 | - | - | - | - |
| SV pieces / path (median) | mafft_fftns2 | 40 | 1 [0-1.25] | 0.975 | 40 | 0 [0-0] | 1 / 8 / 31 | 0.074 |
| SV pieces / path (median) | mafft_fftns2__all | 40 | 1 [0-1.25] | 1.35 | 40 | 0 [0-1] | 3 / 12 / 25 | 0.006 |
| SV pieces / path (median) | mafft_fftnsi__all | 34 | 1 [0-1] | 0.912 | 34 | 0 [0-0] | 4 / 7 / 23 | 0.51 |
| SV pieces / path (median) | mafft_linsi | 40 | 1 [0-1] | 0.725 | 40 | 0 [0-0] | 3 / 3 / 34 | 0.97 |
| SV pieces / path (median) | mafft_linsi__all | 27 | 0 [0-1] | 0.667 | 27 | 0 [0-0] | 2 / 1 / 24 | 1 |
| SV pieces / path (median) | mafft_einsi | 40 | 1 [0-1] | 0.75 | 40 | 0 [0-0] | 3 / 3 / 34 | 0.72 |
| SV pieces / path (median) | mafft_einsi__all | 22 | 0 [0-1] | 0.682 | 22 | 0 [0-0] | 1 / 1 / 20 | 1 |
| SV pieces / path (median) | mafft_ginsi | 40 | 1 [0-1] | 0.738 | 40 | 0 [0-0] | 3 / 3 / 34 | 0.78 |
| SV pieces / path (median) | mafft_ginsi__all | 22 | 0 [0-1] | 0.727 | 22 | 0 [0-0] | 1 / 1 / 20 | 1 |
| SV pieces / path (median) | poa_abpoa | 40 | 0.5 [0-1] | 0.75 | 40 | 0 [0-0] | 4 / 3 / 33 | 0.8 |
| SV pieces / path (median) | poa_abpoa__all | 40 | 0 [0-1] | 0.875 | 40 | 0 [0-0] | 5 / 7 / 28 | 0.72 |
| SV pieces / path (median) | poa_spoa | 39 | 1 [0-1] | 1.205 | 39 | 0 [0-0] | 2 / 7 / 30 | 0.12 |
| SV pieces / path (median) | poa_spoa__all | 37 | 0 [0-1] | 1.095 | 37 | 0 [0-0] | 5 / 7 / 25 | 0.37 |
| SV pieces / path (median) | poa_abpoa_mc | 40 | 1 [0-1] | 0.85 | 40 | 0 [0-0] | 3 / 4 / 33 | 0.67 |
| SV pieces / path (median) | poa_abpoa_mc__all | 40 | 1 [0-1] | 0.7 | 40 | 0 [0-0] | 3 / 4 / 33 | 1 |
| SV pieces / path (median) | unit_aware | 40 | 0 [0-1] | 1.275 | 40 | 0 [0-0] | 4 / 8 / 28 | 0.45 |
| SV pieces / path (median) | unit_aware__all | 40 | 0 [0-1] | 1.25 | 40 | 0 [0-0] | 5 / 7 / 28 | 0.43 |
| indel bp / net length change | mc | 40 | 1.113 [1-1.927] | 1.858 | - | - | - | - |
| indel bp / net length change | mafft_fftns2 | 40 | 1.106 [1.008-1.725] | 1.942 | 40 | 0 [-0.073-0.01] | 17 / 12 / 11 | 0.28 |
| indel bp / net length change | mafft_fftns2__all | 40 | 1.202 [1.01-3.047] | 4.996 | 40 | 0.003 [-0.01-0.644] | 12 / 21 / 7 | 0.083 |
| indel bp / net length change | mafft_fftnsi__all | 34 | 1.12 [1.002-1.849] | 3.6 | 34 | 0 [-0.096-0.009] | 12 / 11 / 11 | 0.66 |
| indel bp / net length change | mafft_linsi | 40 | 1.02 [1-1.154] | 1.377 | 40 | -0.02 [-0.177-0] | 24 / 2 / 14 | 1.5e-07 |
| indel bp / net length change | mafft_linsi__all | 27 | 1.018 [1-1.119] | 1.484 | 27 | -0.013 [-0.205-0] | 15 / 2 / 10 | 0.00015 |
| indel bp / net length change | mafft_einsi | 40 | 1.014 [1-1.113] | 1.37 | 40 | -0.02 [-0.179-0] | 24 / 0 / 16 | 1.2e-07 |
| indel bp / net length change | mafft_einsi__all | 22 | 1.008 [1-1.097] | 1.325 | 22 | -0.013 [-0.231-0] | 12 / 0 / 10 | 0.00049 |
| indel bp / net length change | mafft_ginsi | 40 | 1.02 [1-1.154] | 1.36 | 40 | -0.03 [-0.267-0] | 25 / 2 / 13 | 7.5e-08 |
| indel bp / net length change | mafft_ginsi__all | 22 | 1.017 [1-1.097] | 1.253 | 22 | -0.037 [-0.322-0] | 13 / 2 / 7 | 0.00043 |
| indel bp / net length change | poa_abpoa | 40 | 1.016 [1-1.122] | 1.343 | 40 | -0.011 [-0.26-0] | 25 / 0 / 15 | 6e-08 |
| indel bp / net length change | poa_abpoa__all | 40 | 1.108 [1-1.694] | 1.559 | 40 | 0 [-0.185-0] | 17 / 9 / 14 | 0.059 |
| indel bp / net length change | poa_spoa | 39 | 1.016 [1-1.119] | 1.364 | 39 | -0.028 [-0.368-0] | 26 / 1 / 12 | 1.5e-07 |
| indel bp / net length change | poa_spoa__all | 37 | 1.046 [1-1.222] | 1.479 | 37 | -0.002 [-0.469-0] | 20 / 6 / 11 | 0.00036 |
| indel bp / net length change | poa_abpoa_mc | 40 | 1.031 [1-1.17] | 1.483 | 40 | -0.007 [-0.114-0] | 21 / 1 / 18 | 1.6e-05 |
| indel bp / net length change | poa_abpoa_mc__all | 40 | 1.13 [1-1.867] | 1.746 | 40 | 0 [-0.032-0.004] | 15 / 10 / 15 | 0.63 |
| indel bp / net length change | unit_aware | 40 | 1.009 [1-1.162] | 1.264 | 40 | -0.034 [-0.285-0] | 26 / 1 / 13 | 3e-08 |
| indel bp / net length change | unit_aware__all | 40 | 1.022 [1-1.245] | 1.421 | 40 | -0.021 [-0.181-0] | 24 / 4 / 12 | 0.00015 |
| k-mer extra positions (frac) | mc | 40 | 0.053 [0.009-0.157] | 0.091 | - | - | - | - |
| k-mer extra positions (frac) | mafft_fftns2 | 40 | 0.052 [0.011-0.145] | 0.092 | 40 | -0.001 [-0.019-0.034] | 22 / 17 / 1 | 1 |
| k-mer extra positions (frac) | mafft_fftns2__all | 40 | 0.059 [0.033-0.223] | 0.132 | 40 | 0.017 [-0.002-0.064] | 13 / 25 / 2 | 0.0035 |
| k-mer extra positions (frac) | mafft_fftnsi__all | 34 | 0.056 [0.032-0.213] | 0.115 | 34 | 0.016 [-0.001-0.042] | 9 / 22 / 3 | 0.011 |
| k-mer extra positions (frac) | mafft_linsi | 40 | 0.051 [0.006-0.108] | 0.072 | 40 | -0.004 [-0.033-0] | 29 / 10 / 1 | 0.012 |
| k-mer extra positions (frac) | mafft_linsi__all | 27 | 0.065 [0.009-0.111] | 0.076 | 27 | -0.002 [-0.022-0.006] | 18 / 8 / 1 | 0.35 |
| k-mer extra positions (frac) | mafft_einsi | 40 | 0.032 [0.004-0.103] | 0.064 | 40 | -0.005 [-0.038--0.002] | 33 / 6 / 1 | 0.00071 |
| k-mer extra positions (frac) | mafft_einsi__all | 22 | 0.015 [0.001-0.078] | 0.046 | 22 | -0.004 [-0.048--0.001] | 17 / 4 / 1 | 0.024 |
| k-mer extra positions (frac) | mafft_ginsi | 40 | 0.051 [0.006-0.109] | 0.074 | 40 | -0.004 [-0.033-0.001] | 29 / 10 / 1 | 0.016 |
| k-mer extra positions (frac) | mafft_ginsi__all | 22 | 0.062 [0.005-0.109] | 0.071 | 22 | -0.003 [-0.038-0.008] | 15 / 6 / 1 | 0.36 |
| k-mer extra positions (frac) | poa_abpoa | 40 | 0.038 [0.004-0.123] | 0.08 | 40 | -0.005 [-0.025--0.001] | 31 / 8 / 1 | 0.032 |
| k-mer extra positions (frac) | poa_abpoa__all | 40 | 0.062 [0.011-0.242] | 0.131 | 40 | 0.003 [-0.003-0.09] | 18 / 21 / 1 | 0.044 |
| k-mer extra positions (frac) | poa_spoa | 39 | 0.058 [0.008-0.15] | 0.098 | 39 | -0.001 [-0.015-0.037] | 21 / 17 / 1 | 0.85 |
| k-mer extra positions (frac) | poa_spoa__all | 37 | 0.146 [0.019-0.269] | 0.161 | 37 | 0.082 [-0.001-0.137] | 10 / 26 / 1 | 5.6e-05 |
| k-mer extra positions (frac) | poa_abpoa_mc | 40 | 0.032 [0.005-0.094] | 0.06 | 40 | -0.008 [-0.036--0.003] | 36 / 3 / 1 | 4.4e-06 |
| k-mer extra positions (frac) | poa_abpoa_mc__all | 40 | 0.045 [0.008-0.165] | 0.089 | 40 | -0.002 [-0.021-0.009] | 26 / 13 / 1 | 0.43 |
| k-mer extra positions (frac) | unit_aware | 40 | 0.071 [0.008-0.146] | 0.103 | 40 | -0 [-0.006-0.058] | 20 / 19 / 1 | 0.26 |
| k-mer extra positions (frac) | unit_aware__all | 40 | 0.092 [0.013-0.197] | 0.129 | 40 | 0.007 [-0.001-0.094] | 14 / 25 / 1 | 0.0049 |
| truth edits to best graph path (h1+h2) | mc | 40 | 0 [0-0.25] | 0.425 | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2 | 40 | 0 [0-0] | 0.45 | 40 | 0 [0-0] | 1 / 1 / 38 | 1 |
| truth edits to best graph path (h1+h2) | mafft_fftns2__all | 40 | 0 [0-0] | 0.35 | 40 | 0 [0-0] | 2 / 0 / 38 | 0.5 |
| truth edits to best graph path (h1+h2) | mafft_fftnsi__all | 34 | 0 [0-0] | 0.324 | 34 | 0 [0-0] | 2 / 0 / 32 | 0.5 |
| truth edits to best graph path (h1+h2) | mafft_linsi | 40 | 0 [0-0.25] | 0.6 | 40 | 0 [0-0] | 0 / 2 / 38 | 0.5 |
| truth edits to best graph path (h1+h2) | mafft_linsi__all | 27 | 0 [0-0] | 0.259 | 27 | 0 [0-0] | 0 / 0 / 27 | - |
| truth edits to best graph path (h1+h2) | mafft_einsi | 40 | 0 [0-0.25] | 0.55 | 40 | 0 [0-0] | 0 / 2 / 38 | 0.5 |
| truth edits to best graph path (h1+h2) | mafft_einsi__all | 22 | 0 [0-0] | 0.318 | 22 | 0 [0-0] | 0 / 0 / 22 | - |
| truth edits to best graph path (h1+h2) | mafft_ginsi | 40 | 0 [0-0.25] | 0.625 | 40 | 0 [0-0] | 0 / 2 / 38 | 0.5 |
| truth edits to best graph path (h1+h2) | mafft_ginsi__all | 22 | 0 [0-0] | 0.318 | 22 | 0 [0-0] | 0 / 0 / 22 | - |
| truth edits to best graph path (h1+h2) | poa_abpoa | 40 | 0 [0-0.25] | 0.4 | 40 | 0 [0-0] | 1 / 0 / 39 | 1 |
| truth edits to best graph path (h1+h2) | poa_abpoa__all | 40 | 0 [0-0] | 0.375 | 40 | 0 [0-0] | 1 / 0 / 39 | 1 |
| truth edits to best graph path (h1+h2) | poa_spoa | 39 | 0 [0-0] | 0.436 | 39 | 0 [0-0] | 0 / 1 / 38 | 1 |
| truth edits to best graph path (h1+h2) | poa_spoa__all | 37 | 0 [0-0] | 0.378 | 37 | 0 [0-0] | 1 / 0 / 36 | 1 |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc | 40 | 0 [0-0.25] | 0.475 | 40 | 0 [0-0] | 0 / 1 / 39 | 1 |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc__all | 40 | 0 [0-0] | 0.375 | 40 | 0 [0-0] | 2 / 0 / 38 | 0.5 |
| truth edits to best graph path (h1+h2) | unit_aware | 40 | 0 [0-0.25] | 0.375 | 40 | 0 [0-0] | 2 / 0 / 38 | 0.5 |
| truth edits to best graph path (h1+h2) | unit_aware__all | 40 | 0 [0-0] | 0.375 | 40 | 0 [0-0] | 1 / 0 / 39 | 1 |
| truth-by-graph truvari F1, raw | mc | 24 | 1 [1-1] | 0.975 | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftns2 | 26 | 0.833 [0.259-1] | 0.611 | 24 | 0 [-0.629-0] | 0 / 11 / 13 | 0.00098 |
| truth-by-graph truvari F1, raw | mafft_fftns2__all | 29 | 0.571 [0-1] | 0.536 | 24 | -0.256 [-0.562-0] | 0 / 13 / 11 | 0.00024 |
| truth-by-graph truvari F1, raw | mafft_fftnsi__all | 23 | 0.667 [0.25-1] | 0.621 | 22 | -0.057 [-0.688-0] | 0 / 11 / 11 | 0.00098 |
| truth-by-graph truvari F1, raw | mafft_linsi | 24 | 1 [1-1] | 0.848 | 24 | 0 [0-0] | 1 / 4 / 19 | 0.12 |
| truth-by-graph truvari F1, raw | mafft_linsi__all | 15 | 1 [0.619-1] | 0.828 | 15 | 0 [-0.381-0] | 0 / 5 / 10 | 0.062 |
| truth-by-graph truvari F1, raw | mafft_einsi | 25 | 1 [0.667-1] | 0.774 | 24 | 0 [0-0] | 1 / 5 / 18 | 0.062 |
| truth-by-graph truvari F1, raw | mafft_einsi__all | 11 | 1 [0.786-1] | 0.842 | 11 | 0 [-0.214-0] | 0 / 3 / 8 | 0.25 |
| truth-by-graph truvari F1, raw | mafft_ginsi | 25 | 1 [0.667-1] | 0.774 | 24 | 0 [0-0] | 1 / 5 / 18 | 0.062 |
| truth-by-graph truvari F1, raw | mafft_ginsi__all | 11 | 1 [0.486-1] | 0.717 | 11 | 0 [-0.514-0] | 0 / 5 / 6 | 0.062 |
| truth-by-graph truvari F1, raw | poa_abpoa | 25 | 1 [0.4-1] | 0.719 | 24 | 0 [-0.407-0] | 0 / 7 / 17 | 0.016 |
| truth-by-graph truvari F1, raw | poa_abpoa__all | 24 | 0.417 [0-1] | 0.458 | 24 | -0.5 [-1-0] | 0 / 15 / 9 | 6.1e-05 |
| truth-by-graph truvari F1, raw | poa_spoa | 24 | 0.667 [0-1] | 0.532 | 23 | -0.333 [-1-0] | 0 / 12 / 11 | 0.00049 |
| truth-by-graph truvari F1, raw | poa_spoa__all | 24 | 0 [0-1] | 0.359 | 23 | -1 [-1-0] | 0 / 16 / 7 | 3.1e-05 |
| truth-by-graph truvari F1, raw | poa_abpoa_mc | 24 | 1 [1-1] | 0.906 | 24 | 0 [0-0] | 1 / 3 / 20 | 0.25 |
| truth-by-graph truvari F1, raw | poa_abpoa_mc__all | 26 | 1 [0.571-1] | 0.755 | 24 | 0 [-0.1-0] | 0 / 6 / 18 | 0.031 |
| truth-by-graph truvari F1, raw | unit_aware | 24 | 0.383 [0-1] | 0.476 | 24 | -0.45 [-1-0] | 0 / 14 / 10 | 0.00012 |
| truth-by-graph truvari F1, raw | unit_aware__all | 24 | 0.45 [0-1] | 0.482 | 24 | -0.45 [-1-0] | 0 / 14 / 10 | 0.00012 |
| truth-by-graph truvari F1, refined | mc | 24 | 1 [1-1] | 0.975 | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftns2 | 26 | 1 [0.625-1] | 0.796 | 24 | 0 [0-0] | 0 / 5 / 19 | 0.062 |
| truth-by-graph truvari F1, refined | mafft_fftns2__all | 29 | 1 [0.429-1] | 0.729 | 24 | 0 [0-0] | 0 / 5 / 19 | 0.062 |
| truth-by-graph truvari F1, refined | mafft_fftnsi__all | 23 | 1 [1-1] | 0.889 | 22 | 0 [0-0] | 0 / 3 / 19 | 0.25 |
| truth-by-graph truvari F1, refined | mafft_linsi | 24 | 1 [1-1] | 0.958 | 24 | 0 [0-0] | 1 / 1 / 22 | 1 |
| truth-by-graph truvari F1, refined | mafft_linsi__all | 15 | 1 [1-1] | 0.967 | 15 | 0 [0-0] | 0 / 1 / 14 | 1 |
| truth-by-graph truvari F1, refined | mafft_einsi | 25 | 1 [1-1] | 0.907 | 24 | 0 [0-0] | 1 / 1 / 22 | 1 |
| truth-by-graph truvari F1, refined | mafft_einsi__all | 11 | 1 [1-1] | 1 | 11 | 0 [0-0] | 0 / 0 / 11 | - |
| truth-by-graph truvari F1, refined | mafft_ginsi | 25 | 1 [1-1] | 0.92 | 24 | 0 [0-0] | 1 / 1 / 22 | 1 |
| truth-by-graph truvari F1, refined | mafft_ginsi__all | 11 | 1 [1-1] | 1 | 11 | 0 [0-0] | 0 / 0 / 11 | - |
| truth-by-graph truvari F1, refined | poa_abpoa | 25 | 1 [1-1] | 0.84 | 24 | 0 [0-0] | 0 / 3 / 21 | 0.25 |
| truth-by-graph truvari F1, refined | poa_abpoa__all | 24 | 1 [0.875-1] | 0.771 | 24 | 0 [-0.1-0] | 0 / 6 / 18 | 0.031 |
| truth-by-graph truvari F1, refined | poa_spoa | 24 | 1 [1-1] | 0.792 | 23 | 0 [0-0] | 0 / 4 / 19 | 0.12 |
| truth-by-graph truvari F1, refined | poa_spoa__all | 24 | 1 [0-1] | 0.632 | 23 | 0 [-0.75-0] | 0 / 9 / 14 | 0.0039 |
| truth-by-graph truvari F1, refined | poa_abpoa_mc | 24 | 1 [1-1] | 0.965 | 24 | 0 [0-0] | 1 / 1 / 22 | 1 |
| truth-by-graph truvari F1, refined | poa_abpoa_mc__all | 26 | 1 [1-1] | 0.865 | 24 | 0 [0-0] | 0 / 2 / 22 | 0.5 |
| truth-by-graph truvari F1, refined | unit_aware | 24 | 1 [1-1] | 0.833 | 24 | 0 [0-0] | 0 / 4 / 20 | 0.12 |
| truth-by-graph truvari F1, refined | unit_aware__all | 24 | 1 [0.375-1] | 0.729 | 24 | 0 [-0.425-0] | 0 / 7 / 17 | 0.016 |
| reads: redundant (cross-walk) only | mc | 40 | 0.054 [0.012-0.131] | 0.084 | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftns2 | 40 | 0.069 [0.021-0.149] | 0.112 | 40 | 0.007 [-0.006-0.051] | 15 / 21 / 4 | 0.044 |
| reads: redundant (cross-walk) only | mafft_fftns2__all | 40 | 0.08 [0.027-0.183] | 0.146 | 40 | 0.015 [-0.003-0.091] | 12 / 26 / 2 | 0.0028 |
| reads: redundant (cross-walk) only | mafft_fftnsi__all | 34 | 0.077 [0.025-0.182] | 0.127 | 34 | 0.016 [-0.007-0.081] | 10 / 22 / 2 | 0.018 |
| reads: redundant (cross-walk) only | mafft_linsi | 40 | 0.042 [0.01-0.098] | 0.066 | 40 | 0 [-0.017-0.006] | 18 / 15 / 7 | 0.27 |
| reads: redundant (cross-walk) only | mafft_linsi__all | 27 | 0.057 [0.015-0.089] | 0.07 | 27 | 0 [-0.018-0.009] | 13 / 10 / 4 | 0.43 |
| reads: redundant (cross-walk) only | mafft_einsi | 40 | 0.031 [0.011-0.082] | 0.062 | 40 | -0.002 [-0.027-0.004] | 22 / 13 / 5 | 0.062 |
| reads: redundant (cross-walk) only | mafft_einsi__all | 22 | 0.03 [0.008-0.064] | 0.047 | 22 | 0 [-0.056-0.006] | 10 / 7 / 5 | 0.21 |
| reads: redundant (cross-walk) only | mafft_ginsi | 40 | 0.042 [0.01-0.098] | 0.067 | 40 | 0 [-0.017-0.006] | 18 / 15 / 7 | 0.29 |
| reads: redundant (cross-walk) only | mafft_ginsi__all | 22 | 0.055 [0.011-0.078] | 0.067 | 22 | 0 [-0.013-0.014] | 8 / 10 / 4 | 0.94 |
| reads: redundant (cross-walk) only | poa_abpoa | 40 | 0.035 [0.011-0.124] | 0.088 | 40 | -0.002 [-0.014-0.007] | 21 / 15 / 4 | 0.49 |
| reads: redundant (cross-walk) only | poa_abpoa__all | 40 | 0.058 [0.015-0.187] | 0.138 | 40 | 0.011 [-0.004-0.066] | 11 / 25 / 4 | 0.0082 |
| reads: redundant (cross-walk) only | poa_spoa | 39 | 0.057 [0.013-0.167] | 0.109 | 39 | 0 [-0.007-0.025] | 15 / 19 / 5 | 0.37 |
| reads: redundant (cross-walk) only | poa_spoa__all | 37 | 0.157 [0.022-0.272] | 0.181 | 37 | 0.049 [0.001-0.21] | 6 / 28 / 3 | 3.8e-05 |
| reads: redundant (cross-walk) only | poa_abpoa_mc | 40 | 0.035 [0.011-0.12] | 0.073 | 40 | -0.002 [-0.023-0.005] | 21 / 14 / 5 | 0.19 |
| reads: redundant (cross-walk) only | poa_abpoa_mc__all | 40 | 0.039 [0.012-0.126] | 0.087 | 40 | 0.002 [-0.01-0.017] | 16 / 20 / 4 | 0.56 |
| reads: redundant (cross-walk) only | unit_aware | 40 | 0.061 [0.016-0.134] | 0.102 | 40 | 0.001 [-0.007-0.017] | 15 / 20 / 5 | 0.34 |
| reads: redundant (cross-walk) only | unit_aware__all | 40 | 0.065 [0.02-0.203] | 0.132 | 40 | 0.008 [-0.002-0.065] | 12 / 25 / 3 | 0.012 |
| reads: placements / max per path | mc | 40 | 1.071 [1.012-1.168] | 1.16 | - | - | - | - |
| reads: placements / max per path | mafft_fftns2 | 40 | 1.094 [1.024-1.225] | 1.223 | 40 | 0.007 [-0.007-0.096] | 15 / 22 / 3 | 0.024 |
| reads: placements / max per path | mafft_fftns2__all | 40 | 1.104 [1.029-1.317] | 1.315 | 40 | 0.022 [-0.005-0.216] | 13 / 25 / 2 | 0.0015 |
| reads: placements / max per path | mafft_fftnsi__all | 34 | 1.101 [1.03-1.283] | 1.272 | 34 | 0.019 [-0.009-0.172] | 12 / 20 / 2 | 0.011 |
| reads: placements / max per path | mafft_linsi | 40 | 1.057 [1.014-1.115] | 1.113 | 40 | -0 [-0.023-0.005] | 20 / 15 / 5 | 0.17 |
| reads: placements / max per path | mafft_linsi__all | 27 | 1.067 [1.017-1.166] | 1.101 | 27 | -0.002 [-0.025-0.009] | 14 / 10 / 3 | 0.39 |
| reads: placements / max per path | mafft_einsi | 40 | 1.043 [1.011-1.115] | 1.102 | 40 | -0.003 [-0.045-0.004] | 23 / 13 / 4 | 0.043 |
| reads: placements / max per path | mafft_einsi__all | 22 | 1.036 [1.01-1.086] | 1.064 | 22 | 0 [-0.062-0.006] | 10 / 8 / 4 | 0.18 |
| reads: placements / max per path | mafft_ginsi | 40 | 1.061 [1.014-1.115] | 1.113 | 40 | -0 [-0.023-0.005] | 20 / 15 / 5 | 0.19 |
| reads: placements / max per path | mafft_ginsi__all | 22 | 1.068 [1.011-1.118] | 1.094 | 22 | 0 [-0.028-0.014] | 9 / 10 / 3 | 0.78 |
| reads: placements / max per path | poa_abpoa | 40 | 1.042 [1.013-1.17] | 1.14 | 40 | -0.002 [-0.02-0.008] | 21 / 16 / 3 | 0.44 |
| reads: placements / max per path | poa_abpoa__all | 40 | 1.092 [1.015-1.318] | 1.24 | 40 | 0.012 [-0.004-0.13] | 12 / 26 / 2 | 0.0083 |
| reads: placements / max per path | poa_spoa | 39 | 1.075 [1.014-1.255] | 1.186 | 39 | 0 [-0.008-0.027] | 16 / 19 / 4 | 0.46 |
| reads: placements / max per path | poa_spoa__all | 37 | 1.238 [1.031-1.45] | 1.317 | 37 | 0.117 [0.001-0.268] | 5 / 30 / 2 | 2.2e-05 |
| reads: placements / max per path | poa_abpoa_mc | 40 | 1.044 [1.012-1.141] | 1.111 | 40 | -0.004 [-0.046-0.005] | 23 / 13 / 4 | 0.08 |
| reads: placements / max per path | poa_abpoa_mc__all | 40 | 1.047 [1.012-1.215] | 1.167 | 40 | 0.001 [-0.012-0.029] | 17 / 20 / 3 | 0.41 |
| reads: placements / max per path | unit_aware | 40 | 1.082 [1.016-1.191] | 1.167 | 40 | 0.003 [-0.011-0.017] | 15 / 22 / 3 | 0.46 |
| reads: placements / max per path | unit_aware__all | 40 | 1.093 [1.022-1.326] | 1.24 | 40 | 0.009 [0-0.14] | 9 / 29 / 2 | 0.0056 |

| method | regions scored (raw / refined / phab) | raw TP-base / FN / TP-comp / FP | raw F1 | refined TP-base / FN / TP-comp / FP | refined F1 | phab TP-base / FN / TP-comp / FP | phab F1 |
|---|---|---|---|---|---|---|---|
| mc | 40 / 40 / 40 | 34 / 0 / 34 / 3 | 0.958 | 34 / 0 / 34 / 3 | 0.958 | 33 / 0 / 33 / 0 | 1 |
| mafft_fftns2 | 40 / 40 / 40 | 26 / 8 / 26 / 29 | 0.584 | 32 / 2 / 41 / 14 | 0.832 | 29 / 0 / 29 / 2 | 0.967 |
| mafft_fftns2__all | 40 / 40 / 40 | 24 / 10 / 24 / 65 | 0.39 | 29 / 5 / 45 / 44 | 0.635 | 26 / 0 / 26 / 2 | 0.963 |
| mafft_fftnsi__all | 34 / 34 / 34 | 21 / 11 / 21 / 27 | 0.525 | 29 / 3 / 41 / 7 | 0.879 | 24 / 0 / 24 / 2 | 0.96 |
| mafft_linsi | 40 / 40 / 40 | 27 / 7 / 27 / 9 | 0.771 | 31 / 3 / 36 / 0 | 0.954 | 29 / 0 / 29 / 0 | 1 |
| mafft_linsi__all | 27 / 27 / 27 | 19 / 6 / 19 / 8 | 0.731 | 23 / 2 / 27 / 0 | 0.958 | 21 / 0 / 21 / 0 | 1 |
| mafft_einsi | 40 / 40 / 40 | 26 / 8 / 26 / 11 | 0.732 | 31 / 3 / 35 / 2 | 0.929 | 29 / 0 / 29 / 0 | 1 |
| mafft_einsi__all | 22 / 22 / 22 | 15 / 4 / 15 / 7 | 0.732 | 19 / 0 / 22 / 0 | 1 | 18 / 0 / 18 / 0 | 1 |
| mafft_ginsi | 40 / 40 / 40 | 26 / 8 / 26 / 12 | 0.722 | 31 / 3 / 37 / 1 | 0.942 | 29 / 0 / 29 / 0 | 1 |
| mafft_ginsi__all | 22 / 22 / 22 | 14 / 5 / 14 / 11 | 0.636 | 19 / 0 / 25 / 0 | 1 | 18 / 0 / 18 / 0 | 1 |
| poa_abpoa | 40 / 40 / 40 | 23 / 11 / 23 / 14 | 0.648 | 28 / 6 / 35 / 2 | 0.881 | 27 / 0 / 27 / 0 | 1 |
| poa_abpoa__all | 40 / 40 / 40 | 16 / 18 / 16 / 30 | 0.4 | 26 / 8 / 42 / 4 | 0.832 | 27 / 0 / 27 / 0 | 1 |
| poa_spoa | 39 / 39 / 39 | 17 / 16 / 17 / 34 | 0.405 | 27 / 6 / 35 / 16 | 0.746 | 27 / 0 / 27 / 0 | 1 |
| poa_spoa__all | 37 / 37 / 37 | 12 / 21 / 12 / 36 | 0.296 | 24 / 9 / 32 / 16 | 0.696 | 25 / 0 / 25 / 0 | 1 |
| poa_abpoa_mc | 40 / 40 / 40 | 30 / 4 / 30 / 10 | 0.811 | 32 / 2 / 39 / 1 | 0.958 | 30 / 0 / 30 / 0 | 1 |
| poa_abpoa_mc__all | 40 / 40 / 40 | 27 / 7 / 27 / 13 | 0.73 | 31 / 3 / 35 / 5 | 0.893 | 29 / 0 / 29 / 2 | 0.967 |
| unit_aware | 40 / 40 / 40 | 15 / 19 / 14 / 41 | 0.323 | 28 / 6 / 55 / 0 | 0.903 | 26 / 0 / 26 / 0 | 1 |
| unit_aware__all | 40 / 40 / 40 | 14 / 20 / 14 / 44 | 0.304 | 26 / 8 / 42 / 16 | 0.744 | 26 / 0 / 26 / 0 | 1 |

## control_vntr_correct (25 regions)

| method | n | nodes / CHM13 kb | cost/opt, all pairs | excess edits / kb | affine cost/opt, all | unaligned homology bp / kb | SV pieces / path (median) | indel bp / net length change | k-mer extra positions (frac) | truth edits to best graph path (h1+h2) | truth-by-graph truvari F1, raw | truth-by-graph truvari F1, refined | reads: redundant (cross-walk) only | reads: placements / max per path | regions with cost/opt <= 1.1 | pooled truvari F1 raw / refined |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc | 25 | 85.8 | 1.015 | 1.95 | 1.012 | 0.28 | 1 | 1.047 | 0.041 | 0 | 1 | 1 | 0.048 | 1.05 | 17/25 | 0.985 / 0.985 |
| mafft_fftns2 | 25 | 18.9 | 1.057 | 2.89 | 1.013 | 0 | 1 | 1.069 | 0.029 | 0 | 1 | 1 | 0.066 | 1.077 | 15/25 | 0.735 / 0.989 |
| mafft_fftns2__all | 25 | 19.8 | 1.042 | 3.4 | 1.02 | 0.01 | 1 | 1.081 | 0.038 | 0 | 1 | 1 | 0.095 | 1.105 | 15/25 | 0.61 / 0.974 |
| mafft_fftnsi__all | 25 | 19.2 | 1.036 | 3.17 | 1.01 | 0.02 | 1 | 1.063 | 0.039 | 0 | 1 | 1 | 0.065 | 1.081 | 18/25 | 0.601 / 0.985 |
| mafft_linsi | 25 | 18.6 | 1.015 | 1.67 | 1.005 | 0 | 1 | 1.009 | 0.026 | 0 | 1 | 1 | 0.037 | 1.046 | 22/25 | 0.839 / 0.971 |
| mafft_linsi__all | 20 | 15.4 | 1.034 | 1.77 | 1.017 | 0.055 | 1 | 1.008 | 0.022 | 0 | 1 | 1 | 0.025 | 1.034 | 18/20 | 0.863 / 1 |
| mafft_einsi | 25 | 15.3 | 1.007 | 1.54 | 1.002 | 0 | 1 | 1.009 | 0.021 | 0 | 1 | 1 | 0.031 | 1.037 | 22/25 | 0.822 / 0.985 |
| mafft_einsi__all | 18 | 14.2 | 1.006 | 0.53 | 1.001 | 0 | 1 | 1.006 | 0.018 | 0 | 1 | 1 | 0.023 | 1.025 | 16/18 | 0.877 / 1 |
| mafft_ginsi | 25 | 30.2 | 1.034 | 2.19 | 1.003 | 0.01 | 1 | 1.026 | 0.026 | 0 | 1 | 1 | 0.05 | 1.051 | 21/25 | 0.855 / 0.971 |
| mafft_ginsi__all | 18 | 19.6 | 1.034 | 1.77 | 1.014 | 0.045 | 1 | 1.009 | 0.028 | 0 | 1 | 1 | 0.023 | 1.028 | 15/18 | 0.842 / 1 |
| poa_abpoa | 25 | 18.6 | 1.006 | 1.54 | 1.007 | 0 | 1 | 1.009 | 0.031 | 0 | 1 | 1 | 0.051 | 1.06 | 22/25 | 0.812 / 0.941 |
| poa_abpoa__all | 25 | 67.9 | 1.077 | 9.33 | 1.141 | 0.07 | 1 | 1.075 | 0.136 | 0 | 0.4 | 1 | 0.169 | 1.199 | 15/25 | 0.463 / 0.893 |
| poa_spoa | 25 | 18.3 | 1.005 | 0.93 | 1.008 | 0 | 1 | 1.007 | 0.026 | 0 | 1 | 1 | 0.033 | 1.05 | 24/25 | 0.675 / 0.924 |
| poa_spoa__all | 24 | 45.7 | 1.066 | 3.32 | 1.162 | 0.04 | 1 | 1.027 | 0.095 | 0 | 0 | 1 | 0.096 | 1.122 | 15/24 | 0.308 / 0.889 |
| poa_abpoa_mc | 25 | 15.3 | 1.007 | 1.54 | 1.002 | 0 | 1 | 1.016 | 0.019 | 0 | 1 | 1 | 0.025 | 1.041 | 22/25 | 0.913 / 0.985 |
| poa_abpoa_mc__all | 25 | 21.7 | 1.024 | 1.7 | 1.011 | 0.02 | 1 | 1.04 | 0.025 | 0 | 1 | 1 | 0.05 | 1.055 | 21/25 | 0.707 / 0.959 |
| unit_aware | 25 | 25.8 | 1.015 | 1.12 | 1.024 | 0 | 1 | 1.009 | 0.033 | 0 | 1 | 1 | 0.057 | 1.071 | 25/25 | 0.527 / 0.944 |
| unit_aware__all | 25 | 33.3 | 1.017 | 1.62 | 1.046 | 0 | 1 | 1.009 | 0.044 | 0 | 0.25 | 1 | 0.074 | 1.081 | 22/25 | 0.428 / 0.874 |

| metric | method | n | median [q1-q3] | mean | paired n | median diff [q1-q3] | better / worse / equal | p |
|---|---|---|---|---|---|---|---|---|
| nodes / CHM13 kb | mc | 25 | 85.8 [52.7-136] | 104 | - | - | - | - |
| nodes / CHM13 kb | mafft_fftns2 | 25 | 18.9 [14.7-58.3] | 48.3 | 25 | -43.8 [-76.0--22.9] | 22 / 3 / 0 | 4.5e-05 |
| nodes / CHM13 kb | mafft_fftns2__all | 25 | 19.8 [14.2-52.6] | 48.5 | 25 | -48.6 [-78.7--26.2] | 23 / 2 / 0 | 0.0001 |
| nodes / CHM13 kb | mafft_fftnsi__all | 25 | 19.2 [14.2-52.2] | 46.9 | 25 | -49.8 [-77.2--25.9] | 23 / 2 / 0 | 3.2e-05 |
| nodes / CHM13 kb | mafft_linsi | 25 | 18.6 [13.9-53.6] | 39.2 | 25 | -49.8 [-77.2--26.2] | 23 / 2 / 0 | 2e-06 |
| nodes / CHM13 kb | mafft_linsi__all | 20 | 15.4 [12.8-54.5] | 33.5 | 20 | -42.8 [-76.9--25.1] | 19 / 1 / 0 | 3.6e-05 |
| nodes / CHM13 kb | mafft_einsi | 25 | 15.3 [12.6-29.7] | 30.2 | 25 | -49.8 [-89.2--26.7] | 25 / 0 / 0 | 6e-08 |
| nodes / CHM13 kb | mafft_einsi__all | 18 | 14.2 [10.3-21.2] | 20.8 | 18 | -46.0 [-77.1--26.3] | 18 / 0 / 0 | 7.6e-06 |
| nodes / CHM13 kb | mafft_ginsi | 25 | 30.2 [14.1-58.4] | 43.0 | 25 | -49.8 [-77.2--26.2] | 23 / 2 / 0 | 1e-05 |
| nodes / CHM13 kb | mafft_ginsi__all | 18 | 19.6 [13.7-58.4] | 35.5 | 18 | -42.8 [-74.9--23.0] | 16 / 2 / 0 | 0.00067 |
| nodes / CHM13 kb | poa_abpoa | 25 | 18.6 [12.6-77.1] | 50.5 | 25 | -34.1 [-73.9--19.2] | 22 / 3 / 0 | 5.4e-05 |
| nodes / CHM13 kb | poa_abpoa__all | 25 | 67.9 [18.9-119] | 80.4 | 25 | -21.7 [-36.3--2.49] | 19 / 6 / 0 | 0.08 |
| nodes / CHM13 kb | poa_spoa | 25 | 18.3 [12.2-80.5] | 51.5 | 25 | -44.0 [-76.8--12.1] | 20 / 5 / 0 | 0.0012 |
| nodes / CHM13 kb | poa_spoa__all | 24 | 45.7 [21.3-124] | 78.3 | 24 | -15.4 [-38.1-17.8] | 16 / 8 / 0 | 0.14 |
| nodes / CHM13 kb | poa_abpoa_mc | 25 | 15.3 [12.6-29.5] | 29.5 | 25 | -52.9 [-83.7--30.8] | 25 / 0 / 0 | 6e-08 |
| nodes / CHM13 kb | poa_abpoa_mc__all | 25 | 21.7 [13.7-40.9] | 35.5 | 25 | -52.3 [-76.8--26.4] | 25 / 0 / 0 | 6e-08 |
| nodes / CHM13 kb | unit_aware | 25 | 25.8 [14.3-74.2] | 54.7 | 25 | -50.2 [-76.0--12.1] | 19 / 6 / 0 | 0.0081 |
| nodes / CHM13 kb | unit_aware__all | 25 | 33.3 [16.2-55.9] | 51.9 | 25 | -49.4 [-71.7--19.2] | 21 / 4 / 0 | 0.0015 |
| cost/opt, all pairs | mc | 25 | 1.015 [1.002-1.133] | 1.155 | - | - | - | - |
| cost/opt, all pairs | mafft_fftns2 | 25 | 1.057 [1.004-1.196] | 1.232 | 25 | 0 [-0.015-0.108] | 10 / 10 / 5 | 0.76 |
| cost/opt, all pairs | mafft_fftns2__all | 25 | 1.042 [1.007-1.246] | 1.641 | 25 | 0.003 [-0.001-0.029] | 7 / 13 / 5 | 0.058 |
| cost/opt, all pairs | mafft_fftnsi__all | 25 | 1.036 [1.007-1.193] | 1.552 | 25 | 0 [-0.002-0.022] | 8 / 12 / 5 | 0.37 |
| cost/opt, all pairs | mafft_linsi | 25 | 1.015 [1.003-1.049] | 1.047 | 25 | 0 [-0.029-0.001] | 12 / 8 / 5 | 0.07 |
| cost/opt, all pairs | mafft_linsi__all | 20 | 1.034 [1.002-1.051] | 1.051 | 20 | 0 [-0.038-0.008] | 6 / 9 / 5 | 0.52 |
| cost/opt, all pairs | mafft_einsi | 25 | 1.007 [1.002-1.045] | 1.043 | 25 | -0 [-0.016-0] | 14 / 3 / 8 | 0.0053 |
| cost/opt, all pairs | mafft_einsi__all | 18 | 1.006 [1.001-1.042] | 1.047 | 18 | 0 [-0.073-0.001] | 6 / 5 / 7 | 0.21 |
| cost/opt, all pairs | mafft_ginsi | 25 | 1.034 [1.004-1.053] | 1.052 | 25 | 0 [-0.029-0.002] | 12 / 8 / 5 | 0.16 |
| cost/opt, all pairs | mafft_ginsi__all | 18 | 1.034 [1.002-1.056] | 1.059 | 18 | 0 [-0.076-0.011] | 6 / 7 / 5 | 0.45 |
| cost/opt, all pairs | poa_abpoa | 25 | 1.006 [1.002-1.049] | 1.032 | 25 | -0 [-0.131-0] | 13 / 6 / 6 | 0.026 |
| cost/opt, all pairs | poa_abpoa__all | 25 | 1.077 [1.003-1.157] | 1.192 | 25 | 0 [-0.054-0.057] | 12 / 10 / 3 | 1 |
| cost/opt, all pairs | poa_spoa | 25 | 1.005 [1.002-1.026] | 1.02 | 25 | -0.009 [-0.131-0] | 16 / 4 / 5 | 0.0027 |
| cost/opt, all pairs | poa_spoa__all | 24 | 1.066 [1.016-1.126] | 1.09 | 24 | -0 [-0.117-0.082] | 12 / 9 / 3 | 0.47 |
| cost/opt, all pairs | poa_abpoa_mc | 25 | 1.007 [1.002-1.045] | 1.047 | 25 | -0 [-0.014-0] | 14 / 4 / 7 | 0.016 |
| cost/opt, all pairs | poa_abpoa_mc__all | 25 | 1.024 [1.003-1.06] | 1.097 | 25 | 0 [-0.001-0.008] | 9 / 9 / 7 | 0.74 |
| cost/opt, all pairs | unit_aware | 25 | 1.015 [1.003-1.043] | 1.025 | 25 | -0.008 [-0.11-0.002] | 14 / 9 / 2 | 0.025 |
| cost/opt, all pairs | unit_aware__all | 25 | 1.017 [1.004-1.055] | 1.047 | 25 | -0 [-0.048-0.004] | 13 / 11 / 1 | 0.28 |
| excess edits / kb | mc | 25 | 1.95 [0.43-9.61] | 9.954 | - | - | - | - |
| excess edits / kb | mafft_fftns2 | 25 | 2.89 [0.44-20.6] | 13.3 | 25 | 0 [-1.03-2.86] | 10 / 10 / 5 | 0.71 |
| excess edits / kb | mafft_fftns2__all | 25 | 3.4 [1.22-23.8] | 26.2 | 25 | 0.61 [-0.16-4.82] | 7 / 13 / 5 | 0.036 |
| excess edits / kb | mafft_fftnsi__all | 25 | 3.17 [0.5-18.2] | 20.2 | 25 | 0 [-0.25-2.6] | 8 / 12 / 5 | 0.14 |
| excess edits / kb | mafft_linsi | 25 | 1.67 [0.27-8.81] | 7.862 | 25 | 0 [-1.64-0.18] | 12 / 8 / 5 | 0.2 |
| excess edits / kb | mafft_linsi__all | 20 | 1.77 [0.075-9.195] | 4.921 | 20 | 0 [-2.202-0.558] | 6 / 9 / 5 | 0.8 |
| excess edits / kb | mafft_einsi | 25 | 1.54 [0.09-3.34] | 4.284 | 25 | -0.03 [-1.64-0] | 14 / 3 / 8 | 0.0053 |
| excess edits / kb | mafft_einsi__all | 18 | 0.53 [0.045-3.192] | 4.099 | 18 | 0 [-1.685-0.143] | 6 / 5 / 7 | 0.17 |
| excess edits / kb | mafft_ginsi | 25 | 2.19 [0.28-9.96] | 8.756 | 25 | 0 [-1.64-0.45] | 12 / 8 / 5 | 0.55 |
| excess edits / kb | mafft_ginsi__all | 18 | 1.77 [0.045-11.1] | 5.868 | 18 | 0 [-3.327-0.793] | 6 / 7 / 5 | 0.84 |
| excess edits / kb | poa_abpoa | 25 | 1.54 [0.09-3.07] | 3.864 | 25 | -0.04 [-5.47-0] | 13 / 6 / 6 | 0.065 |
| excess edits / kb | poa_abpoa__all | 25 | 9.33 [1-22.1] | 12.5 | 25 | 0 [-3.14-6.51] | 11 / 10 / 4 | 0.68 |
| excess edits / kb | poa_spoa | 25 | 0.93 [0.09-1.67] | 2.295 | 25 | -0.84 [-7.94-0] | 16 / 4 / 5 | 0.0023 |
| excess edits / kb | poa_spoa__all | 24 | 3.32 [0.877-11.7] | 8.955 | 24 | -0.015 [-7.48-4.088] | 12 / 9 / 3 | 0.73 |
| excess edits / kb | poa_abpoa_mc | 25 | 1.54 [0.23-3.32] | 4.842 | 25 | -0.02 [-1.64-0] | 14 / 4 / 7 | 0.01 |
| excess edits / kb | poa_abpoa_mc__all | 25 | 1.7 [0.32-8.33] | 10.3 | 25 | 0 [-0.07-0.45] | 8 / 9 / 8 | 0.85 |
| excess edits / kb | unit_aware | 25 | 1.12 [0.3-3.32] | 2.489 | 25 | -0.79 [-7.94-0.27] | 14 / 9 / 2 | 0.015 |
| excess edits / kb | unit_aware__all | 25 | 1.62 [0.33-3.8] | 4.016 | 25 | -0.24 [-1.94-0.34] | 13 / 11 / 1 | 0.14 |
| affine cost/opt, all | mc | 25 | 1.012 [1.001-1.052] | 1.123 | - | - | - | - |
| affine cost/opt, all | mafft_fftns2 | 25 | 1.013 [1.002-1.192] | 1.21 | 25 | 0 [-0.021-0.184] | 10 / 12 / 3 | 0.5 |
| affine cost/opt, all | mafft_fftns2__all | 25 | 1.02 [1.004-1.188] | 1.612 | 25 | 0.004 [-0.006-0.068] | 8 / 14 / 3 | 0.063 |
| affine cost/opt, all | mafft_fftnsi__all | 25 | 1.01 [1.003-1.088] | 1.525 | 25 | 0.001 [-0.008-0.035] | 9 / 13 / 3 | 0.41 |
| affine cost/opt, all | mafft_linsi | 25 | 1.005 [1-1.035] | 1.025 | 25 | -0.002 [-0.026-0.002] | 13 / 8 / 4 | 0.1 |
| affine cost/opt, all | mafft_linsi__all | 20 | 1.017 [1-1.054] | 1.03 | 20 | 0 [-0.041-0.023] | 7 / 9 / 4 | 0.71 |
| affine cost/opt, all | mafft_einsi | 25 | 1.002 [1-1.007] | 1.012 | 25 | -0.003 [-0.027-0] | 17 / 2 / 6 | 0.00034 |
| affine cost/opt, all | mafft_einsi__all | 18 | 1.001 [1-1.021] | 1.013 | 18 | -0.001 [-0.072-0] | 9 / 4 / 5 | 0.065 |
| affine cost/opt, all | mafft_ginsi | 25 | 1.003 [1-1.042] | 1.037 | 25 | -0 [-0.027-0.006] | 13 / 8 / 4 | 0.32 |
| affine cost/opt, all | mafft_ginsi__all | 18 | 1.014 [1-1.068] | 1.044 | 18 | 0 [-0.072-0.04] | 7 / 7 / 4 | 0.67 |
| affine cost/opt, all | poa_abpoa | 25 | 1.007 [1-1.061] | 1.051 | 25 | 0 [-0.025-0.009] | 11 / 9 / 5 | 0.5 |
| affine cost/opt, all | poa_abpoa__all | 25 | 1.141 [1.009-1.324] | 1.36 | 25 | 0.02 [-0-0.163] | 7 / 17 / 1 | 0.013 |
| affine cost/opt, all | poa_spoa | 25 | 1.008 [1-1.111] | 1.056 | 25 | 0 [-0.021-0.004] | 11 / 9 / 5 | 0.5 |
| affine cost/opt, all | poa_spoa__all | 24 | 1.162 [1.013-1.322] | 1.204 | 24 | 0.012 [-0.004-0.219] | 6 / 16 / 2 | 0.054 |
| affine cost/opt, all | poa_abpoa_mc | 25 | 1.002 [1-1.016] | 1.016 | 25 | -0.002 [-0.025-0] | 16 / 4 / 5 | 0.0056 |
| affine cost/opt, all | poa_abpoa_mc__all | 25 | 1.011 [1.002-1.049] | 1.069 | 25 | 0 [-0.023-0.013] | 9 / 11 / 5 | 0.7 |
| affine cost/opt, all | unit_aware | 25 | 1.024 [1.005-1.124] | 1.074 | 25 | 0.004 [-0.024-0.019] | 9 / 15 / 1 | 0.88 |
| affine cost/opt, all | unit_aware__all | 25 | 1.046 [1.009-1.106] | 1.074 | 25 | 0.008 [-0.005-0.046] | 7 / 17 / 1 | 0.26 |
| unaligned homology bp / kb | mc | 25 | 0.28 [0-0.65] | 1.96 | - | - | - | - |
| unaligned homology bp / kb | mafft_fftns2 | 25 | 0 [0-0.03] | 0.168 | 25 | -0.22 [-0.5-0] | 15 / 1 / 9 | 0.001 |
| unaligned homology bp / kb | mafft_fftns2__all | 25 | 0.01 [0-0.61] | 2.486 | 25 | 0 [-0.4-0.49] | 9 / 8 / 8 | 0.93 |
| unaligned homology bp / kb | mafft_fftnsi__all | 25 | 0.02 [0-0.35] | 1.016 | 25 | -0.03 [-0.4-0] | 13 / 5 / 7 | 0.12 |
| unaligned homology bp / kb | mafft_linsi | 25 | 0 [0-0.2] | 0.452 | 25 | -0.03 [-0.45-0] | 13 / 6 / 6 | 0.018 |
| unaligned homology bp / kb | mafft_linsi__all | 20 | 0.055 [0-0.44] | 0.261 | 20 | 0 [-0.29-0.072] | 8 / 6 / 6 | 0.59 |
| unaligned homology bp / kb | mafft_einsi | 25 | 0 [0-0.13] | 0.238 | 25 | 0 [-0.45-0] | 12 / 5 / 8 | 0.0045 |
| unaligned homology bp / kb | mafft_einsi__all | 18 | 0 [0-0.133] | 2.081 | 18 | 0 [-0.35-0] | 7 / 4 / 7 | 0.29 |
| unaligned homology bp / kb | mafft_ginsi | 25 | 0.01 [0-0.5] | 0.618 | 25 | -0.03 [-0.45-0] | 13 / 6 / 6 | 0.11 |
| unaligned homology bp / kb | mafft_ginsi__all | 18 | 0.045 [0-0.432] | 0.243 | 18 | 0 [-0.24-0.033] | 6 / 6 / 6 | 0.81 |
| unaligned homology bp / kb | poa_abpoa | 25 | 0 [0-0] | 0.095 | 25 | -0.28 [-0.65-0] | 15 / 2 / 8 | 7.6e-05 |
| unaligned homology bp / kb | poa_abpoa__all | 25 | 0.07 [0-0.63] | 0.428 | 25 | 0 [-0.5-0.06] | 11 / 9 / 5 | 0.28 |
| unaligned homology bp / kb | poa_spoa | 25 | 0 [0-0.01] | 0.072 | 25 | -0.28 [-0.65-0] | 15 / 4 / 6 | 0.00028 |
| unaligned homology bp / kb | poa_spoa__all | 24 | 0.04 [0-0.388] | 0.311 | 24 | 0 [-0.537-0.013] | 11 / 7 / 6 | 0.048 |
| unaligned homology bp / kb | poa_abpoa_mc | 25 | 0 [0-0] | 0.267 | 25 | -0.12 [-0.47-0] | 13 / 3 / 9 | 0.0063 |
| unaligned homology bp / kb | poa_abpoa_mc__all | 25 | 0.02 [0-0.11] | 0.421 | 25 | -0.14 [-0.65-0] | 13 / 5 / 7 | 0.058 |
| unaligned homology bp / kb | unit_aware | 25 | 0 [0-0.15] | 0.168 | 25 | -0.23 [-0.47-0] | 16 / 3 / 6 | 0.0032 |
| unaligned homology bp / kb | unit_aware__all | 25 | 0 [0-0.08] | 0.1 | 25 | -0.28 [-0.65-0] | 15 / 4 / 6 | 0.0032 |
| SV pieces / path (median) | mc | 25 | 1 [1-1] | 0.88 | - | - | - | - |
| SV pieces / path (median) | mafft_fftns2 | 25 | 1 [1-1] | 1.12 | 25 | 0 [0-0] | 0 / 3 / 22 | 0.25 |
| SV pieces / path (median) | mafft_fftns2__all | 25 | 1 [0-2] | 1.4 | 25 | 0 [0-1] | 3 / 8 / 14 | 0.066 |
| SV pieces / path (median) | mafft_fftnsi__all | 25 | 1 [0-2] | 1.24 | 25 | 0 [0-0] | 3 / 6 / 16 | 0.29 |
| SV pieces / path (median) | mafft_linsi | 25 | 1 [1-1] | 0.96 | 25 | 0 [0-0] | 0 / 2 / 23 | 0.5 |
| SV pieces / path (median) | mafft_linsi__all | 20 | 1 [0.75-1] | 0.9 | 20 | 0 [0-0] | 0 / 1 / 19 | 1 |
| SV pieces / path (median) | mafft_einsi | 25 | 1 [1-1] | 1 | 25 | 0 [0-0] | 0 / 3 / 22 | 0.25 |
| SV pieces / path (median) | mafft_einsi__all | 18 | 1 [0.25-1] | 0.833 | 18 | 0 [0-0] | 0 / 1 / 17 | 1 |
| SV pieces / path (median) | mafft_ginsi | 25 | 1 [1-1] | 0.96 | 25 | 0 [0-0] | 0 / 2 / 23 | 0.5 |
| SV pieces / path (median) | mafft_ginsi__all | 18 | 1 [0.25-1] | 0.778 | 18 | 0 [0-0] | 0 / 0 / 18 | - |
| SV pieces / path (median) | poa_abpoa | 25 | 1 [1-1] | 0.92 | 25 | 0 [0-0] | 2 / 2 / 21 | 1 |
| SV pieces / path (median) | poa_abpoa__all | 25 | 1 [1-2] | 1.28 | 25 | 0 [0-0] | 0 / 6 / 19 | 0.031 |
| SV pieces / path (median) | poa_spoa | 25 | 1 [1-1] | 0.96 | 25 | 0 [0-0] | 2 / 2 / 21 | 1 |
| SV pieces / path (median) | poa_spoa__all | 24 | 1 [0.75-2] | 1.5 | 24 | 0 [0-1] | 4 / 7 / 13 | 0.13 |
| SV pieces / path (median) | poa_abpoa_mc | 25 | 1 [1-1] | 0.96 | 25 | 0 [0-0] | 0 / 1 / 24 | 1 |
| SV pieces / path (median) | poa_abpoa_mc__all | 25 | 1 [1-1] | 0.96 | 25 | 0 [0-0] | 2 / 3 / 20 | 0.75 |
| SV pieces / path (median) | unit_aware | 25 | 1 [0-1] | 1.04 | 25 | 0 [0-0] | 2 / 3 / 20 | 0.5 |
| SV pieces / path (median) | unit_aware__all | 25 | 1 [0-2] | 1.2 | 25 | 0 [0-1] | 3 / 7 / 15 | 0.13 |
| indel bp / net length change | mc | 25 | 1.047 [1.001-1.19] | 1.23 | - | - | - | - |
| indel bp / net length change | mafft_fftns2 | 25 | 1.069 [1.001-1.246] | 1.219 | 25 | 0 [-0.023-0] | 11 / 6 / 8 | 0.64 |
| indel bp / net length change | mafft_fftns2__all | 25 | 1.081 [1.007-1.398] | 1.811 | 25 | 0 [-0.001-0.02] | 7 / 9 / 9 | 0.46 |
| indel bp / net length change | mafft_fftnsi__all | 25 | 1.063 [1.001-1.156] | 1.728 | 25 | 0 [-0.01-0] | 9 / 6 / 10 | 0.76 |
| indel bp / net length change | mafft_linsi | 25 | 1.009 [1.002-1.125] | 1.153 | 25 | -0 [-0.043-0] | 13 / 4 / 8 | 0.0026 |
| indel bp / net length change | mafft_linsi__all | 20 | 1.008 [1-1.137] | 1.176 | 20 | -0 [-0.012-0] | 10 / 2 / 8 | 0.02 |
| indel bp / net length change | mafft_einsi | 25 | 1.009 [1.001-1.125] | 1.157 | 25 | -0.001 [-0.02-0] | 14 / 2 / 9 | 0.00043 |
| indel bp / net length change | mafft_einsi__all | 18 | 1.006 [1-1.117] | 1.19 | 18 | -0 [-0.016-0] | 9 / 0 / 9 | 0.0039 |
| indel bp / net length change | mafft_ginsi | 25 | 1.026 [1.002-1.125] | 1.154 | 25 | 0 [-0.043-0] | 12 / 5 / 8 | 0.0093 |
| indel bp / net length change | mafft_ginsi__all | 18 | 1.009 [1-1.117] | 1.186 | 18 | 0 [-0.016-0] | 8 / 2 / 8 | 0.064 |
| indel bp / net length change | poa_abpoa | 25 | 1.009 [1.001-1.114] | 1.095 | 25 | -0 [-0.044-0] | 14 / 4 / 7 | 0.03 |
| indel bp / net length change | poa_abpoa__all | 25 | 1.075 [1.005-1.144] | 1.151 | 25 | 0 [-0.081-0.026] | 11 / 8 / 6 | 0.66 |
| indel bp / net length change | poa_spoa | 25 | 1.007 [1-1.093] | 1.077 | 25 | -0 [-0.14-0] | 14 / 4 / 7 | 0.0038 |
| indel bp / net length change | poa_spoa__all | 24 | 1.027 [1.001-1.128] | 1.107 | 24 | -0 [-0.062-0.005] | 12 / 6 / 6 | 0.19 |
| indel bp / net length change | poa_abpoa_mc | 25 | 1.016 [1-1.125] | 1.159 | 25 | 0 [-0.018-0] | 12 / 2 / 11 | 0.0013 |
| indel bp / net length change | poa_abpoa_mc__all | 25 | 1.04 [1.001-1.125] | 1.169 | 25 | 0 [-0.01-0] | 9 / 5 / 11 | 0.11 |
| indel bp / net length change | unit_aware | 25 | 1.009 [1.001-1.069] | 1.075 | 25 | 0 [-0.14-0] | 12 / 4 / 9 | 0.0042 |
| indel bp / net length change | unit_aware__all | 25 | 1.009 [1.002-1.083] | 1.145 | 25 | 0 [-0.018-0] | 12 / 6 / 7 | 0.021 |
| k-mer extra positions (frac) | mc | 25 | 0.041 [0.014-0.07] | 0.069 | - | - | - | - |
| k-mer extra positions (frac) | mafft_fftns2 | 25 | 0.029 [0.019-0.067] | 0.078 | 25 | -0.003 [-0.015-0.015] | 15 / 9 / 1 | 0.77 |
| k-mer extra positions (frac) | mafft_fftns2__all | 25 | 0.038 [0.015-0.116] | 0.09 | 25 | 0 [-0.009-0.035] | 12 / 11 / 2 | 0.56 |
| k-mer extra positions (frac) | mafft_fftnsi__all | 25 | 0.039 [0.015-0.092] | 0.084 | 25 | 0 [-0.015-0.026] | 12 / 11 / 2 | 0.99 |
| k-mer extra positions (frac) | mafft_linsi | 25 | 0.026 [0.011-0.053] | 0.043 | 25 | -0.009 [-0.037-0.006] | 15 / 8 / 2 | 0.092 |
| k-mer extra positions (frac) | mafft_linsi__all | 20 | 0.022 [0.007-0.046] | 0.042 | 20 | -0.004 [-0.015-0.006] | 11 / 7 / 2 | 0.39 |
| k-mer extra positions (frac) | mafft_einsi | 25 | 0.021 [0.007-0.032] | 0.028 | 25 | -0.012 [-0.043--0.005] | 20 / 3 / 2 | 7.9e-06 |
| k-mer extra positions (frac) | mafft_einsi__all | 18 | 0.018 [0.004-0.027] | 0.021 | 18 | -0.009 [-0.025--0.001] | 14 / 2 / 2 | 0.011 |
| k-mer extra positions (frac) | mafft_ginsi | 25 | 0.026 [0.014-0.059] | 0.05 | 25 | -0.009 [-0.036-0.006] | 15 / 8 / 2 | 0.18 |
| k-mer extra positions (frac) | mafft_ginsi__all | 18 | 0.028 [0.006-0.059] | 0.047 | 18 | -0.001 [-0.016-0.009] | 9 / 7 / 2 | 0.67 |
| k-mer extra positions (frac) | poa_abpoa | 25 | 0.031 [0.005-0.104] | 0.063 | 25 | -0.001 [-0.03-0.007] | 13 / 9 / 3 | 0.48 |
| k-mer extra positions (frac) | poa_abpoa__all | 25 | 0.136 [0.026-0.219] | 0.139 | 25 | 0.05 [0-0.154] | 5 / 18 / 2 | 0.00025 |
| k-mer extra positions (frac) | poa_spoa | 25 | 0.026 [0.005-0.106] | 0.069 | 25 | 0 [-0.017-0.012] | 12 / 10 / 3 | 0.82 |
| k-mer extra positions (frac) | poa_spoa__all | 24 | 0.095 [0.043-0.211] | 0.122 | 24 | 0.037 [0.003-0.077] | 5 / 18 / 1 | 0.00085 |
| k-mer extra positions (frac) | poa_abpoa_mc | 25 | 0.019 [0.005-0.028] | 0.027 | 25 | -0.013 [-0.043--0.005] | 20 / 2 / 3 | 1.2e-05 |
| k-mer extra positions (frac) | poa_abpoa_mc__all | 25 | 0.025 [0.009-0.073] | 0.049 | 25 | -0.006 [-0.03-0.001] | 15 / 8 / 2 | 0.11 |
| k-mer extra positions (frac) | unit_aware | 25 | 0.033 [0.014-0.109] | 0.073 | 25 | 0 [-0.015-0.043] | 11 / 12 / 2 | 0.62 |
| k-mer extra positions (frac) | unit_aware__all | 25 | 0.044 [0.021-0.091] | 0.08 | 25 | 0.001 [-0.01-0.036] | 11 / 13 / 1 | 0.49 |
| truth edits to best graph path (h1+h2) | mc | 25 | 0 [0-0] | 0.28 | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2 | 25 | 0 [0-0] | 0.28 | 25 | 0 [0-0] | 0 / 0 / 25 | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2__all | 25 | 0 [0-0] | 0.28 | 25 | 0 [0-0] | 0 / 0 / 25 | - |
| truth edits to best graph path (h1+h2) | mafft_fftnsi__all | 25 | 0 [0-0] | 0.28 | 25 | 0 [0-0] | 0 / 0 / 25 | - |
| truth edits to best graph path (h1+h2) | mafft_linsi | 25 | 0 [0-0] | 0.28 | 25 | 0 [0-0] | 1 / 1 / 23 | 1 |
| truth edits to best graph path (h1+h2) | mafft_linsi__all | 20 | 0 [0-0] | 0.25 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | mafft_einsi | 25 | 0 [0-0] | 0.28 | 25 | 0 [0-0] | 0 / 0 / 25 | - |
| truth edits to best graph path (h1+h2) | mafft_einsi__all | 18 | 0 [0-0] | 0.222 | 18 | 0 [0-0] | 0 / 0 / 18 | - |
| truth edits to best graph path (h1+h2) | mafft_ginsi | 25 | 0 [0-0] | 0.32 | 25 | 0 [0-0] | 1 / 1 / 23 | 1 |
| truth edits to best graph path (h1+h2) | mafft_ginsi__all | 18 | 0 [0-0] | 0.222 | 18 | 0 [0-0] | 0 / 0 / 18 | - |
| truth edits to best graph path (h1+h2) | poa_abpoa | 25 | 0 [0-1] | 0.36 | 25 | 0 [0-0] | 0 / 1 / 24 | 1 |
| truth edits to best graph path (h1+h2) | poa_abpoa__all | 25 | 0 [0-0] | 0.24 | 25 | 0 [0-0] | 1 / 0 / 24 | 1 |
| truth edits to best graph path (h1+h2) | poa_spoa | 25 | 0 [0-0] | 0.28 | 25 | 0 [0-0] | 0 / 0 / 25 | - |
| truth edits to best graph path (h1+h2) | poa_spoa__all | 24 | 0 [0-0.25] | 0.292 | 24 | 0 [0-0] | 0 / 0 / 24 | - |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc | 25 | 0 [0-1] | 0.92 | 25 | 0 [0-0] | 0 / 1 / 24 | 1 |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc__all | 25 | 0 [0-0] | 0.28 | 25 | 0 [0-0] | 0 / 0 / 25 | - |
| truth edits to best graph path (h1+h2) | unit_aware | 25 | 0 [0-1] | 0.36 | 25 | 0 [0-0] | 0 / 1 / 24 | 1 |
| truth edits to best graph path (h1+h2) | unit_aware__all | 25 | 0 [0-0] | 0.28 | 25 | 0 [0-0] | 0 / 0 / 25 | - |
| truth-by-graph truvari F1, raw | mc | 25 | 1 [1-1] | 0.987 | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftns2 | 25 | 1 [0.8-1] | 0.817 | 25 | 0 [0-0] | 0 / 6 / 19 | 0.031 |
| truth-by-graph truvari F1, raw | mafft_fftns2__all | 25 | 1 [0-1] | 0.661 | 25 | 0 [-0.818-0] | 0 / 10 / 15 | 0.002 |
| truth-by-graph truvari F1, raw | mafft_fftnsi__all | 25 | 1 [0-1] | 0.638 | 25 | 0 [-0.778-0] | 0 / 11 / 14 | 0.00098 |
| truth-by-graph truvari F1, raw | mafft_linsi | 25 | 1 [1-1] | 0.853 | 25 | 0 [0-0] | 0 / 4 / 21 | 0.12 |
| truth-by-graph truvari F1, raw | mafft_linsi__all | 20 | 1 [1-1] | 0.903 | 20 | 0 [0-0] | 0 / 3 / 17 | 0.25 |
| truth-by-graph truvari F1, raw | mafft_einsi | 25 | 1 [1-1] | 0.816 | 25 | 0 [0-0] | 0 / 5 / 20 | 0.062 |
| truth-by-graph truvari F1, raw | mafft_einsi__all | 18 | 1 [1-1] | 0.911 | 18 | 0 [0-0] | 0 / 2 / 16 | 0.5 |
| truth-by-graph truvari F1, raw | mafft_ginsi | 25 | 1 [1-1] | 0.893 | 25 | 0 [0-0] | 0 / 3 / 22 | 0.25 |
| truth-by-graph truvari F1, raw | mafft_ginsi__all | 18 | 1 [1-1] | 0.911 | 18 | 0 [0-0] | 0 / 2 / 16 | 0.5 |
| truth-by-graph truvari F1, raw | poa_abpoa | 25 | 1 [1-1] | 0.854 | 25 | 0 [0-0] | 0 / 5 / 20 | 0.062 |
| truth-by-graph truvari F1, raw | poa_abpoa__all | 25 | 0.4 [0-1] | 0.512 | 25 | -0.6 [-1-0] | 1 / 14 / 10 | 0.00018 |
| truth-by-graph truvari F1, raw | poa_spoa | 25 | 1 [0.5-1] | 0.703 | 25 | 0 [-0.5-0] | 0 / 10 / 15 | 0.002 |
| truth-by-graph truvari F1, raw | poa_spoa__all | 24 | 0 [0-0.75] | 0.333 | 24 | -1 [-1--0.25] | 0 / 18 / 6 | 7.6e-06 |
| truth-by-graph truvari F1, raw | poa_abpoa_mc | 25 | 1 [1-1] | 0.931 | 25 | 0 [0-0] | 0 / 2 / 23 | 0.5 |
| truth-by-graph truvari F1, raw | poa_abpoa_mc__all | 25 | 1 [0.5-1] | 0.727 | 25 | 0 [-0.5-0] | 0 / 8 / 17 | 0.0078 |
| truth-by-graph truvari F1, raw | unit_aware | 25 | 1 [0-1] | 0.623 | 25 | 0 [-1-0] | 0 / 11 / 14 | 0.00098 |
| truth-by-graph truvari F1, raw | unit_aware__all | 25 | 0.25 [0-1] | 0.469 | 25 | -0.667 [-1-0] | 0 / 15 / 10 | 6.1e-05 |
| truth-by-graph truvari F1, refined | mc | 25 | 1 [1-1] | 0.987 | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftns2 | 25 | 1 [1-1] | 0.987 | 25 | 0 [0-0] | 0 / 0 / 25 | - |
| truth-by-graph truvari F1, refined | mafft_fftns2__all | 25 | 1 [1-1] | 0.92 | 25 | 0 [0-0] | 0 / 2 / 23 | 0.5 |
| truth-by-graph truvari F1, refined | mafft_fftnsi__all | 25 | 1 [1-1] | 0.92 | 25 | 0 [0-0] | 0 / 2 / 23 | 0.5 |
| truth-by-graph truvari F1, refined | mafft_linsi | 25 | 1 [1-1] | 0.96 | 25 | 0 [0-0] | 0 / 1 / 24 | 1 |
| truth-by-graph truvari F1, refined | mafft_linsi__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | mafft_einsi | 25 | 1 [1-1] | 0.96 | 25 | 0 [0-0] | 0 / 1 / 24 | 1 |
| truth-by-graph truvari F1, refined | mafft_einsi__all | 18 | 1 [1-1] | 1 | 18 | 0 [0-0] | 0 / 0 / 18 | - |
| truth-by-graph truvari F1, refined | mafft_ginsi | 25 | 1 [1-1] | 0.96 | 25 | 0 [0-0] | 0 / 1 / 24 | 1 |
| truth-by-graph truvari F1, refined | mafft_ginsi__all | 18 | 1 [1-1] | 1 | 18 | 0 [0-0] | 0 / 0 / 18 | - |
| truth-by-graph truvari F1, refined | poa_abpoa | 25 | 1 [1-1] | 0.92 | 25 | 0 [0-0] | 0 / 2 / 23 | 0.5 |
| truth-by-graph truvari F1, refined | poa_abpoa__all | 25 | 1 [1-1] | 0.84 | 25 | 0 [0-0] | 1 / 4 / 20 | 0.12 |
| truth-by-graph truvari F1, refined | poa_spoa | 25 | 1 [1-1] | 0.88 | 25 | 0 [0-0] | 0 / 3 / 22 | 0.25 |
| truth-by-graph truvari F1, refined | poa_spoa__all | 24 | 1 [1-1] | 0.792 | 24 | 0 [0-0] | 0 / 5 / 19 | 0.062 |
| truth-by-graph truvari F1, refined | poa_abpoa_mc | 25 | 1 [1-1] | 0.96 | 25 | 0 [0-0] | 0 / 1 / 24 | 1 |
| truth-by-graph truvari F1, refined | poa_abpoa_mc__all | 25 | 1 [1-1] | 0.92 | 25 | 0 [0-0] | 0 / 2 / 23 | 0.5 |
| truth-by-graph truvari F1, refined | unit_aware | 25 | 1 [1-1] | 0.84 | 25 | 0 [0-0] | 0 / 4 / 21 | 0.12 |
| truth-by-graph truvari F1, refined | unit_aware__all | 25 | 1 [1-1] | 0.8 | 25 | 0 [0-0] | 0 / 5 / 20 | 0.062 |
| reads: redundant (cross-walk) only | mc | 25 | 0.048 [0.02-0.123] | 0.099 | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftns2 | 25 | 0.066 [0.033-0.129] | 0.122 | 25 | 0 [-0.039-0.069] | 10 / 11 / 4 | 0.59 |
| reads: redundant (cross-walk) only | mafft_fftns2__all | 25 | 0.095 [0.045-0.201] | 0.149 | 25 | 0.035 [-0.003-0.078] | 8 / 16 / 1 | 0.019 |
| reads: redundant (cross-walk) only | mafft_fftnsi__all | 25 | 0.065 [0.032-0.191] | 0.135 | 25 | 0.02 [-0.002-0.066] | 7 / 17 / 1 | 0.053 |
| reads: redundant (cross-walk) only | mafft_linsi | 25 | 0.037 [0.012-0.103] | 0.072 | 25 | 0 [-0.039-0.03] | 11 / 9 / 5 | 0.55 |
| reads: redundant (cross-walk) only | mafft_linsi__all | 20 | 0.025 [0.012-0.099] | 0.069 | 20 | 0 [-0.007-0.033] | 8 / 8 / 4 | 0.94 |
| reads: redundant (cross-walk) only | mafft_einsi | 25 | 0.031 [0.012-0.055] | 0.06 | 25 | -0.002 [-0.06-0.003] | 13 / 8 / 4 | 0.1 |
| reads: redundant (cross-walk) only | mafft_einsi__all | 18 | 0.023 [0.011-0.054] | 0.039 | 18 | -0.003 [-0.031-0] | 10 / 4 / 4 | 0.12 |
| reads: redundant (cross-walk) only | mafft_ginsi | 25 | 0.05 [0.012-0.109] | 0.077 | 25 | 0 [-0.039-0.03] | 11 / 9 / 5 | 0.65 |
| reads: redundant (cross-walk) only | mafft_ginsi__all | 18 | 0.023 [0.011-0.106] | 0.073 | 18 | 0 [-0.01-0.034] | 8 / 7 / 3 | 0.89 |
| reads: redundant (cross-walk) only | poa_abpoa | 25 | 0.051 [0.017-0.098] | 0.098 | 25 | 0 [-0.036-0.031] | 10 / 11 / 4 | 0.93 |
| reads: redundant (cross-walk) only | poa_abpoa__all | 25 | 0.169 [0.05-0.287] | 0.194 | 25 | 0.059 [0.009-0.179] | 3 / 19 / 3 | 3.3e-05 |
| reads: redundant (cross-walk) only | poa_spoa | 25 | 0.033 [0.012-0.113] | 0.103 | 25 | 0 [-0.029-0.011] | 11 / 8 / 6 | 0.83 |
| reads: redundant (cross-walk) only | poa_spoa__all | 24 | 0.096 [0.054-0.291] | 0.181 | 24 | 0.057 [0.014-0.109] | 0 / 20 / 4 | 1.9e-06 |
| reads: redundant (cross-walk) only | poa_abpoa_mc | 25 | 0.025 [0.012-0.048] | 0.05 | 25 | -0.002 [-0.057-0] | 13 / 6 / 6 | 0.019 |
| reads: redundant (cross-walk) only | poa_abpoa_mc__all | 25 | 0.05 [0.021-0.112] | 0.088 | 25 | 0 [-0.028-0.027] | 10 / 11 / 4 | 0.97 |
| reads: redundant (cross-walk) only | unit_aware | 25 | 0.057 [0.012-0.103] | 0.107 | 25 | 0 [-0.023-0.03] | 8 / 12 / 5 | 0.65 |
| reads: redundant (cross-walk) only | unit_aware__all | 25 | 0.074 [0.024-0.136] | 0.114 | 25 | 0.008 [0-0.035] | 5 / 18 / 2 | 0.096 |
| reads: placements / max per path | mc | 25 | 1.05 [1.024-1.138] | 1.17 | - | - | - | - |
| reads: placements / max per path | mafft_fftns2 | 25 | 1.077 [1.033-1.143] | 1.22 | 25 | 0 [-0.039-0.072] | 12 / 11 / 2 | 0.71 |
| reads: placements / max per path | mafft_fftns2__all | 25 | 1.105 [1.051-1.249] | 1.333 | 25 | 0.035 [-0.002-0.108] | 7 / 17 / 1 | 0.0065 |
| reads: placements / max per path | mafft_fftnsi__all | 25 | 1.081 [1.038-1.249] | 1.274 | 25 | 0.024 [0-0.088] | 6 / 18 / 1 | 0.011 |
| reads: placements / max per path | mafft_linsi | 25 | 1.046 [1.016-1.117] | 1.108 | 25 | -0.002 [-0.041-0.012] | 13 / 9 / 3 | 0.43 |
| reads: placements / max per path | mafft_linsi__all | 20 | 1.034 [1.018-1.103] | 1.087 | 20 | 0 [-0.012-0.013] | 9 / 8 / 3 | 0.75 |
| reads: placements / max per path | mafft_einsi | 25 | 1.037 [1.016-1.059] | 1.096 | 25 | -0.007 [-0.091-0] | 16 / 6 / 3 | 0.021 |
| reads: placements / max per path | mafft_einsi__all | 18 | 1.025 [1.012-1.079] | 1.056 | 18 | -0.002 [-0.033-0.004] | 10 / 5 / 3 | 0.23 |
| reads: placements / max per path | mafft_ginsi | 25 | 1.051 [1.016-1.124] | 1.116 | 25 | -0.002 [-0.039-0.026] | 13 / 9 / 3 | 0.57 |
| reads: placements / max per path | mafft_ginsi__all | 18 | 1.028 [1.015-1.117] | 1.091 | 18 | 0 [-0.016-0.028] | 8 / 7 / 3 | 0.8 |
| reads: placements / max per path | poa_abpoa | 25 | 1.06 [1.024-1.145] | 1.152 | 25 | 0 [-0.039-0.06] | 12 / 10 / 3 | 0.87 |
| reads: placements / max per path | poa_abpoa__all | 25 | 1.199 [1.06-1.388] | 1.322 | 25 | 0.088 [0.012-0.207] | 4 / 19 / 2 | 3.3e-05 |
| reads: placements / max per path | poa_spoa | 25 | 1.05 [1.022-1.202] | 1.176 | 25 | 0 [-0.039-0.011] | 12 / 9 / 4 | 0.97 |
| reads: placements / max per path | poa_spoa__all | 24 | 1.122 [1.059-1.331] | 1.269 | 24 | 0.067 [0.011-0.129] | 2 / 19 / 3 | 0.00013 |
| reads: placements / max per path | poa_abpoa_mc | 25 | 1.041 [1.016-1.06] | 1.084 | 25 | -0.008 [-0.096-0] | 16 / 6 / 3 | 0.023 |
| reads: placements / max per path | poa_abpoa_mc__all | 25 | 1.055 [1.026-1.177] | 1.148 | 25 | 0 [-0.017-0.02] | 11 / 11 / 3 | 0.92 |
| reads: placements / max per path | unit_aware | 25 | 1.071 [1.024-1.206] | 1.158 | 25 | 0.004 [-0.01-0.022] | 8 / 13 / 4 | 0.61 |
| reads: placements / max per path | unit_aware__all | 25 | 1.081 [1.027-1.19] | 1.177 | 25 | 0.005 [0-0.033] | 6 / 17 / 2 | 0.3 |

| method | regions scored (raw / refined / phab) | raw TP-base / FN / TP-comp / FP | raw F1 | refined TP-base / FN / TP-comp / FP | refined F1 | phab TP-base / FN / TP-comp / FP | phab F1 |
|---|---|---|---|---|---|---|---|
| mc | 25 / 25 / 25 | 34 / 0 / 33 / 1 | 0.985 | 34 / 0 / 33 / 1 | 0.985 | 34 / 0 / 33 / 0 | 1 |
| mafft_fftns2 | 25 / 25 / 25 | 30 / 4 / 29 / 17 | 0.735 | 34 / 0 / 45 / 1 | 0.989 | 34 / 0 / 33 / 0 | 1 |
| mafft_fftns2__all | 25 / 25 / 25 | 25 / 9 / 24 / 22 | 0.61 | 33 / 1 / 45 / 1 | 0.974 | 32 / 0 / 31 / 0 | 1 |
| mafft_fftnsi__all | 25 / 25 / 25 | 24 / 10 / 23 / 21 | 0.601 | 33 / 1 / 44 / 0 | 0.985 | 32 / 0 / 31 / 0 | 1 |
| mafft_linsi | 25 / 25 / 25 | 30 / 4 / 28 / 7 | 0.839 | 33 / 1 / 34 / 1 | 0.971 | 34 / 0 / 32 / 0 | 1 |
| mafft_linsi__all | 20 / 20 / 20 | 26 / 2 / 25 / 6 | 0.863 | 28 / 0 / 31 / 0 | 1 | 27 / 0 / 26 / 0 | 1 |
| mafft_einsi | 25 / 25 / 25 | 29 / 5 / 27 / 7 | 0.822 | 33 / 1 / 34 / 0 | 0.985 | 34 / 0 / 32 / 0 | 1 |
| mafft_einsi__all | 18 / 18 / 18 | 22 / 2 / 21 / 4 | 0.877 | 24 / 0 / 25 / 0 | 1 | 23 / 0 / 22 / 0 | 1 |
| mafft_ginsi | 25 / 25 / 25 | 31 / 3 / 29 / 7 | 0.855 | 33 / 1 / 35 / 1 | 0.971 | 34 / 0 / 32 / 0 | 1 |
| mafft_ginsi__all | 18 / 18 / 18 | 22 / 2 / 21 / 6 | 0.842 | 24 / 0 / 27 / 0 | 1 | 23 / 0 / 22 / 0 | 1 |
| poa_abpoa | 25 / 25 / 25 | 28 / 6 / 28 / 7 | 0.812 | 31 / 3 / 34 / 1 | 0.941 | 32 / 0 / 32 / 0 | 1 |
| poa_abpoa__all | 25 / 25 / 25 | 19 / 15 / 19 / 29 | 0.463 | 29 / 5 / 45 / 3 | 0.893 | 32 / 0 / 32 / 0 | 1 |
| poa_spoa | 25 / 25 / 25 | 24 / 10 / 22 / 12 | 0.675 | 30 / 4 / 33 / 1 | 0.924 | 32 / 0 / 30 / 0 | 1 |
| poa_spoa__all | 24 / 24 / 24 | 12 / 21 / 12 / 33 | 0.308 | 28 / 5 / 42 / 3 | 0.889 | 30 / 0 / 30 / 0 | 1 |
| poa_abpoa_mc | 25 / 25 / 25 | 32 / 2 / 31 / 4 | 0.913 | 33 / 1 / 35 / 0 | 0.985 | 33 / 0 / 32 / 0 | 1 |
| poa_abpoa_mc__all | 25 / 25 / 25 | 26 / 8 / 25 / 13 | 0.707 | 33 / 1 / 36 / 2 | 0.959 | 33 / 0 / 31 / 0 | 1 |
| unit_aware | 25 / 25 / 25 | 22 / 12 / 20 / 25 | 0.527 | 31 / 3 / 44 / 1 | 0.944 | 34 / 0 / 32 / 0 | 1 |
| unit_aware__all | 25 / 25 / 25 | 19 / 15 / 17 / 32 | 0.428 | 32 / 2 / 40 / 9 | 0.874 | 33 / 0 / 31 / 0 | 1 |

## hotspot_other (24 regions)

| method | n | nodes / CHM13 kb | cost/opt, all pairs | excess edits / kb | affine cost/opt, all | unaligned homology bp / kb | SV pieces / path (median) | indel bp / net length change | k-mer extra positions (frac) | truth edits to best graph path (h1+h2) | truth-by-graph truvari F1, raw | truth-by-graph truvari F1, refined | reads: redundant (cross-walk) only | reads: placements / max per path | regions with cost/opt <= 1.1 | pooled truvari F1 raw / refined |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc | 24 | 492 | 1.625 | 96.3 | 1.598 | 22.7 | 5 | 1.548 | 0.435 | 0 | 0.182 | 1 | 0.295 | 2.108 | 3/24 | 0.309 / 0.974 |
| mafft_fftns2 | 24 | 261 | 1.313 | 54.6 | 1.258 | 2.945 | 3 | 1.135 | 0.363 | 2.5 | 0.297 | 1 | 0.251 | 1.961 | 6/24 | 0.317 / 0.954 |
| mafft_fftns2__all | 24 | 308 | 1.36 | 72.0 | 1.3 | 7.31 | 3 | 1.255 | 0.418 | 1.5 | 0.25 | 1 | 0.263 | 2.305 | 1/24 | 0.311 / 0.964 |
| mafft_fftnsi__all | 17 | 244 | 1.296 | 56.3 | 1.268 | 4.44 | 2 | 1.239 | 0.369 | 1 | 0.182 | 1 | 0.29 | 2.092 | 2/17 | 0.223 / 0.927 |
| mafft_linsi | 24 | 293 | 1.164 | 28.0 | 1.138 | 1.325 | 2 | 1.079 | 0.299 | 3.5 | 0.4 | 1 | 0.146 | 1.544 | 8/24 | 0.401 / 0.971 |
| mafft_linsi__all | 6 | 292 | 1.159 | 29.8 | 1.153 | 4.555 | 2 | 1.106 | 0.312 | 0.5 | 0.2 | 1 | 0.129 | 1.495 | 1/6 | 0.326 / 1 |
| mafft_einsi | 24 | 247 | 1.221 | 37.9 | 1.162 | 3.88 | 2.5 | 1.141 | 0.299 | 5.5 | 0.422 | 1 | 0.171 | 1.445 | 6/24 | 0.423 / 0.915 |
| mafft_einsi__all | 4 | 176 | 1.33 | 55.6 | 1.278 | 8.235 | 2.5 | 1.271 | 0.396 | 3.5 | 0.5 | 1 | 0.152 | 1.501 | 1/4 | 0.408 / 1 |
| mafft_ginsi | 24 | 308 | 1.163 | 28.4 | 1.138 | 1.13 | 2.5 | 1.069 | 0.304 | 2.5 | 0.4 | 1 | 0.141 | 1.542 | 8/24 | 0.402 / 0.982 |
| mafft_ginsi__all | 4 | 258 | 1.163 | 26.7 | 1.151 | 5.395 | 2 | 1.043 | 0.349 | 1.5 | 0.2 | 1 | 0.136 | 1.475 | 1/4 | 0.308 / 1 |
| poa_abpoa | 24 | 270 | 1.15 | 26.6 | 1.197 | 1.245 | 3 | 1.058 | 0.352 | 2 | 0.3 | 1 | 0.225 | 1.851 | 9/24 | 0.329 / 0.846 |
| poa_abpoa__all | 23 | 474 | 1.244 | 47.0 | 1.476 | 3.53 | 2 | 1.101 | 0.515 | 0 | 0.2 | 1 | 0.331 | 2.34 | 3/23 | 0.244 / 0.803 |
| poa_spoa | 24 | 316 | 1.125 | 16.8 | 1.23 | 0.94 | 2 | 1.031 | 0.387 | 1.5 | 0.167 | 1 | 0.275 | 1.974 | 10/24 | 0.193 / 0.827 |
| poa_spoa__all | 22 | 557 | 1.239 | 44.2 | 1.529 | 3.855 | 1 | 1.082 | 0.57 | 0 | 0 | 1 | 0.348 | 2.486 | 3/22 | 0.161 / 0.939 |
| poa_abpoa_mc | 24 | 214 | 1.284 | 39.5 | 1.176 | 2.09 | 3 | 1.199 | 0.295 | 5 | 0.4 | 1 | 0.195 | 1.623 | 6/24 | 0.431 / 0.915 |
| poa_abpoa_mc__all | 24 | 315 | 1.356 | 79.5 | 1.343 | 7.85 | 3 | 1.281 | 0.432 | 0 | 0.37 | 1 | 0.268 | 2.185 | 3/24 | 0.383 / 0.89 |
| unit_aware | 24 | 398 | 1.053 | 11.3 | 1.218 | 0.605 | 0.5 | 1.024 | 0.433 | 0.5 | 0.133 | 1 | 0.261 | 2.024 | 20/24 | 0.203 / 0.895 |
| unit_aware__all | 24 | 502 | 1.097 | 16.3 | 1.312 | 1.27 | 0 | 1.03 | 0.475 | 0 | 0 | 1 | 0.278 | 2.077 | 12/24 | 0.182 / 0.905 |

| metric | method | n | median [q1-q3] | mean | paired n | median diff [q1-q3] | better / worse / equal | p |
|---|---|---|---|---|---|---|---|---|
| nodes / CHM13 kb | mc | 24 | 492 [337-872] | 691 | - | - | - | - |
| nodes / CHM13 kb | mafft_fftns2 | 24 | 261 [130-626] | 415 | 24 | -206 [-378--83.8] | 24 / 0 / 0 | 1.2e-07 |
| nodes / CHM13 kb | mafft_fftns2__all | 24 | 308 [183-699] | 495 | 24 | -173 [-298--57.7] | 22 / 2 / 0 | 6.4e-05 |
| nodes / CHM13 kb | mafft_fftnsi__all | 17 | 244 [187-477] | 423 | 17 | -135 [-206--63.2] | 15 / 2 / 0 | 0.0011 |
| nodes / CHM13 kb | mafft_linsi | 24 | 293 [178-502] | 397 | 24 | -223 [-400--128] | 23 / 1 / 0 | 3.6e-07 |
| nodes / CHM13 kb | mafft_linsi__all | 6 | 292 [163-376] | 279 | 6 | -140 [-293--92.6] | 6 / 0 / 0 | 0.031 |
| nodes / CHM13 kb | mafft_einsi | 24 | 247 [147-425] | 347 | 24 | -281 [-471--142] | 23 / 1 / 0 | 2.4e-07 |
| nodes / CHM13 kb | mafft_einsi__all | 4 | 176 [71.6-335] | 231 | 4 | -150 [-252--131] | 4 / 0 / 0 | 0.12 |
| nodes / CHM13 kb | mafft_ginsi | 24 | 308 [173-490] | 397 | 24 | -226 [-412--123] | 23 / 1 / 0 | 3.6e-07 |
| nodes / CHM13 kb | mafft_ginsi__all | 4 | 258 [98.0-435] | 274 | 4 | -108 [-229--67.0] | 3 / 1 / 0 | 0.25 |
| nodes / CHM13 kb | poa_abpoa | 24 | 270 [206-506] | 368 | 24 | -247 [-471--140] | 22 / 2 / 0 | 6.6e-06 |
| nodes / CHM13 kb | poa_abpoa__all | 23 | 474 [370-748] | 613 | 23 | -78.4 [-132--9.55] | 17 / 6 / 0 | 0.016 |
| nodes / CHM13 kb | poa_spoa | 24 | 316 [236-526] | 414 | 24 | -230 [-438--87.8] | 22 / 2 / 0 | 3e-06 |
| nodes / CHM13 kb | poa_spoa__all | 22 | 557 [394-744] | 698 | 22 | 30.0 [-38.6-114] | 10 / 12 / 0 | 0.23 |
| nodes / CHM13 kb | poa_abpoa_mc | 24 | 214 [141-478] | 359 | 24 | -284 [-498--166] | 24 / 0 / 0 | 1.2e-07 |
| nodes / CHM13 kb | poa_abpoa_mc__all | 24 | 315 [200-696] | 522 | 24 | -150 [-193--96.2] | 21 / 3 / 0 | 6.6e-06 |
| nodes / CHM13 kb | unit_aware | 24 | 398 [284-641] | 514 | 24 | -146 [-290--48.9] | 20 / 4 / 0 | 0.00028 |
| nodes / CHM13 kb | unit_aware__all | 24 | 502 [392-794] | 630 | 24 | -50.5 [-119-41.0] | 14 / 10 / 0 | 0.17 |
| cost/opt, all pairs | mc | 24 | 1.625 [1.268-1.838] | 1.688 | - | - | - | - |
| cost/opt, all pairs | mafft_fftns2 | 24 | 1.313 [1.147-1.417] | 1.312 | 24 | -0.284 [-0.516--0.159] | 22 / 2 / 0 | 2.3e-06 |
| cost/opt, all pairs | mafft_fftns2__all | 24 | 1.36 [1.207-1.576] | 1.449 | 24 | -0.103 [-0.347--0.046] | 20 / 4 / 0 | 5.3e-05 |
| cost/opt, all pairs | mafft_fftnsi__all | 17 | 1.296 [1.227-1.52] | 1.437 | 17 | -0.17 [-0.348--0.096] | 16 / 1 / 0 | 7.6e-05 |
| cost/opt, all pairs | mafft_linsi | 24 | 1.164 [1.09-1.282] | 1.196 | 24 | -0.37 [-0.611--0.191] | 24 / 0 / 0 | 1.2e-07 |
| cost/opt, all pairs | mafft_linsi__all | 6 | 1.159 [1.118-1.278] | 1.206 | 6 | -0.294 [-0.448--0.098] | 6 / 0 / 0 | 0.031 |
| cost/opt, all pairs | mafft_einsi | 24 | 1.221 [1.1-1.401] | 1.263 | 24 | -0.278 [-0.523--0.182] | 24 / 0 / 0 | 1.2e-07 |
| cost/opt, all pairs | mafft_einsi__all | 4 | 1.33 [1.2-1.435] | 1.304 | 4 | -0.175 [-0.208--0.109] | 3 / 1 / 0 | 0.25 |
| cost/opt, all pairs | mafft_ginsi | 24 | 1.163 [1.09-1.252] | 1.199 | 24 | -0.369 [-0.61--0.191] | 24 / 0 / 0 | 1.2e-07 |
| cost/opt, all pairs | mafft_ginsi__all | 4 | 1.163 [1.113-1.22] | 1.169 | 4 | -0.294 [-0.376--0.196] | 4 / 0 / 0 | 0.12 |
| cost/opt, all pairs | poa_abpoa | 24 | 1.15 [1.073-1.207] | 1.154 | 24 | -0.438 [-0.634--0.198] | 24 / 0 / 0 | 1.2e-07 |
| cost/opt, all pairs | poa_abpoa__all | 23 | 1.244 [1.158-1.31] | 1.271 | 23 | -0.342 [-0.537--0.139] | 21 / 2 / 0 | 3.3e-06 |
| cost/opt, all pairs | poa_spoa | 24 | 1.125 [1.052-1.165] | 1.119 | 24 | -0.48 [-0.712--0.229] | 24 / 0 / 0 | 1.2e-07 |
| cost/opt, all pairs | poa_spoa__all | 22 | 1.239 [1.172-1.299] | 1.259 | 22 | -0.393 [-0.552--0.168] | 20 / 2 / 0 | 4.8e-06 |
| cost/opt, all pairs | poa_abpoa_mc | 24 | 1.284 [1.105-1.399] | 1.274 | 24 | -0.258 [-0.44--0.185] | 22 / 2 / 0 | 2.3e-06 |
| cost/opt, all pairs | poa_abpoa_mc__all | 24 | 1.356 [1.16-1.562] | 1.418 | 24 | -0.132 [-0.283--0.022] | 18 / 6 / 0 | 0.0016 |
| cost/opt, all pairs | unit_aware | 24 | 1.053 [1.037-1.092] | 1.064 | 24 | -0.533 [-0.758--0.243] | 24 / 0 / 0 | 1.2e-07 |
| cost/opt, all pairs | unit_aware__all | 24 | 1.097 [1.055-1.126] | 1.092 | 24 | -0.473 [-0.733--0.229] | 24 / 0 / 0 | 1.2e-07 |
| excess edits / kb | mc | 24 | 96.3 [55.0-177] | 127 | - | - | - | - |
| excess edits / kb | mafft_fftns2 | 24 | 54.6 [37.3-79.0] | 58.3 | 24 | -48.3 [-105--25.2] | 22 / 2 / 0 | 3.9e-06 |
| excess edits / kb | mafft_fftns2__all | 24 | 72.0 [50.9-109] | 78.6 | 24 | -26.0 [-71.6--3.045] | 20 / 4 / 0 | 0.00028 |
| excess edits / kb | mafft_fftnsi__all | 17 | 56.3 [40.0-82.9] | 60.2 | 17 | -31.9 [-53.5--13.2] | 16 / 1 / 0 | 0.00029 |
| excess edits / kb | mafft_linsi | 24 | 28.0 [22.5-49.7] | 38.5 | 24 | -66.1 [-129--27.2] | 24 / 0 / 0 | 1.2e-07 |
| excess edits / kb | mafft_linsi__all | 6 | 29.8 [26.0-54.5] | 39.2 | 6 | -57.0 [-69.6--25.1] | 6 / 0 / 0 | 0.031 |
| excess edits / kb | mafft_einsi | 24 | 37.9 [25.5-68.1] | 49.7 | 24 | -60.5 [-108--27.4] | 24 / 0 / 0 | 1.2e-07 |
| excess edits / kb | mafft_einsi__all | 4 | 55.6 [42.5-71.4] | 58.4 | 4 | -28.6 [-33.3--20.3] | 3 / 1 / 0 | 0.25 |
| excess edits / kb | mafft_ginsi | 24 | 28.4 [20.6-49.2] | 39.3 | 24 | -66.1 [-115--27.5] | 24 / 0 / 0 | 1.2e-07 |
| excess edits / kb | mafft_ginsi__all | 4 | 26.7 [23.4-36.0] | 32.7 | 4 | -57.1 [-68.3--39.4] | 4 / 0 / 0 | 0.12 |
| excess edits / kb | poa_abpoa | 24 | 26.6 [18.0-42.2] | 31.5 | 24 | -69.3 [-144--49.7] | 24 / 0 / 0 | 1.2e-07 |
| excess edits / kb | poa_abpoa__all | 23 | 47.0 [34.5-63.5] | 50.4 | 23 | -54.4 [-131--26.2] | 21 / 2 / 0 | 7.9e-06 |
| excess edits / kb | poa_spoa | 24 | 16.8 [13.2-31.0] | 22.8 | 24 | -80.5 [-150--46.3] | 24 / 0 / 0 | 1.2e-07 |
| excess edits / kb | poa_spoa__all | 22 | 44.2 [33.0-58.4] | 46.3 | 22 | -52.9 [-127--37.9] | 20 / 2 / 0 | 4.8e-06 |
| excess edits / kb | poa_abpoa_mc | 24 | 39.5 [26.5-79.3] | 53.7 | 24 | -52.9 [-85.5--27.1] | 22 / 2 / 0 | 1e-05 |
| excess edits / kb | poa_abpoa_mc__all | 24 | 79.5 [49.8-92.7] | 86.8 | 24 | -41.0 [-62.9--6.597] | 18 / 6 / 0 | 0.0025 |
| excess edits / kb | unit_aware | 24 | 11.3 [7.862-17.1] | 13.9 | 24 | -85.6 [-164--51.9] | 24 / 0 / 0 | 1.2e-07 |
| excess edits / kb | unit_aware__all | 24 | 16.3 [11.8-23.2] | 18.9 | 24 | -77.9 [-159--49.8] | 24 / 0 / 0 | 1.2e-07 |
| affine cost/opt, all | mc | 24 | 1.598 [1.259-1.727] | 1.584 | - | - | - | - |
| affine cost/opt, all | mafft_fftns2 | 24 | 1.258 [1.097-1.329] | 1.247 | 24 | -0.292 [-0.474--0.173] | 23 / 1 / 0 | 8.3e-07 |
| affine cost/opt, all | mafft_fftns2__all | 24 | 1.3 [1.186-1.491] | 1.39 | 24 | -0.127 [-0.249--0.035] | 21 / 3 / 0 | 7.6e-05 |
| affine cost/opt, all | mafft_fftnsi__all | 17 | 1.268 [1.203-1.442] | 1.369 | 17 | -0.164 [-0.269--0.124] | 16 / 1 / 0 | 7.6e-05 |
| affine cost/opt, all | mafft_linsi | 24 | 1.138 [1.064-1.205] | 1.139 | 24 | -0.383 [-0.595--0.201] | 24 / 0 / 0 | 1.2e-07 |
| affine cost/opt, all | mafft_linsi__all | 6 | 1.153 [1.123-1.241] | 1.182 | 6 | -0.297 [-0.461--0.146] | 6 / 0 / 0 | 0.031 |
| affine cost/opt, all | mafft_einsi | 24 | 1.162 [1.064-1.277] | 1.182 | 24 | -0.327 [-0.521--0.196] | 24 / 0 / 0 | 1.2e-07 |
| affine cost/opt, all | mafft_einsi__all | 4 | 1.278 [1.18-1.369] | 1.271 | 4 | -0.208 [-0.252--0.143] | 4 / 0 / 0 | 0.12 |
| affine cost/opt, all | mafft_ginsi | 24 | 1.138 [1.066-1.203] | 1.142 | 24 | -0.386 [-0.572--0.199] | 24 / 0 / 0 | 1.2e-07 |
| affine cost/opt, all | mafft_ginsi__all | 4 | 1.151 [1.104-1.201] | 1.155 | 4 | -0.297 [-0.38--0.219] | 4 / 0 / 0 | 0.12 |
| affine cost/opt, all | poa_abpoa | 24 | 1.197 [1.097-1.289] | 1.206 | 24 | -0.265 [-0.485--0.096] | 23 / 1 / 0 | 3.6e-07 |
| affine cost/opt, all | poa_abpoa__all | 23 | 1.476 [1.269-1.577] | 1.455 | 23 | -0.083 [-0.252-0.063] | 14 / 9 / 0 | 0.11 |
| affine cost/opt, all | poa_spoa | 24 | 1.23 [1.114-1.323] | 1.22 | 24 | -0.28 [-0.549--0.086] | 23 / 1 / 0 | 3e-06 |
| affine cost/opt, all | poa_spoa__all | 22 | 1.529 [1.4-1.687] | 1.551 | 22 | 0.017 [-0.192-0.107] | 10 / 12 / 0 | 0.59 |
| affine cost/opt, all | poa_abpoa_mc | 24 | 1.176 [1.074-1.345] | 1.211 | 24 | -0.278 [-0.445--0.185] | 23 / 1 / 0 | 1.7e-06 |
| affine cost/opt, all | poa_abpoa_mc__all | 24 | 1.343 [1.177-1.507] | 1.397 | 24 | -0.102 [-0.249-0.005] | 17 / 7 / 0 | 0.012 |
| affine cost/opt, all | unit_aware | 24 | 1.218 [1.118-1.299] | 1.218 | 24 | -0.253 [-0.512--0.1] | 22 / 2 / 0 | 8.3e-07 |
| affine cost/opt, all | unit_aware__all | 24 | 1.312 [1.145-1.464] | 1.307 | 24 | -0.233 [-0.395--0.008] | 18 / 6 / 0 | 0.00057 |
| unaligned homology bp / kb | mc | 24 | 22.7 [9.555-37.4] | 30.3 | - | - | - | - |
| unaligned homology bp / kb | mafft_fftns2 | 24 | 2.945 [1.695-4.205] | 3.677 | 24 | -17.9 [-35.2--8.113] | 22 / 2 / 0 | 1.2e-06 |
| unaligned homology bp / kb | mafft_fftns2__all | 24 | 7.31 [3.228-12.2] | 8.892 | 24 | -11.6 [-30.2--3.445] | 20 / 4 / 0 | 1.6e-05 |
| unaligned homology bp / kb | mafft_fftnsi__all | 17 | 4.44 [1.87-6.24] | 5.506 | 17 | -12.1 [-18.8--3.12] | 15 / 2 / 0 | 0.00021 |
| unaligned homology bp / kb | mafft_linsi | 24 | 1.325 [0.792-2.518] | 1.728 | 24 | -20.5 [-36.6--8.233] | 24 / 0 / 0 | 1.2e-07 |
| unaligned homology bp / kb | mafft_linsi__all | 6 | 4.555 [1.978-6.705] | 4.25 | 6 | -16.1 [-39.9--3.938] | 6 / 0 / 0 | 0.031 |
| unaligned homology bp / kb | mafft_einsi | 24 | 3.88 [1.635-7.255] | 4.954 | 24 | -18.4 [-34.9--6.85] | 21 / 3 / 0 | 1.7e-06 |
| unaligned homology bp / kb | mafft_einsi__all | 4 | 8.235 [3.655-13.9] | 9.315 | 4 | -2.8 [-8.49-2.145] | 2 / 2 / 0 | 0.62 |
| unaligned homology bp / kb | mafft_ginsi | 24 | 1.13 [0.627-2.06] | 1.808 | 24 | -21.0 [-36.6--8.143] | 24 / 0 / 0 | 1.2e-07 |
| unaligned homology bp / kb | mafft_ginsi__all | 4 | 5.395 [2.065-9.178] | 5.848 | 4 | -3.94 [-11.8-0.82] | 2 / 2 / 0 | 0.62 |
| unaligned homology bp / kb | poa_abpoa | 24 | 1.245 [0.44-2.188] | 1.51 | 24 | -21.1 [-35.7--8.898] | 22 / 2 / 0 | 6e-07 |
| unaligned homology bp / kb | poa_abpoa__all | 23 | 3.53 [2.575-5.09] | 4.407 | 23 | -17.7 [-34.0--6.785] | 20 / 3 / 0 | 6e-06 |
| unaligned homology bp / kb | poa_spoa | 24 | 0.94 [0.462-2.312] | 1.564 | 24 | -21.1 [-35.7--8.985] | 23 / 1 / 0 | 3.6e-07 |
| unaligned homology bp / kb | poa_spoa__all | 22 | 3.855 [2.272-5.877] | 4.582 | 22 | -16.8 [-34.6--5.967] | 19 / 3 / 0 | 1.2e-05 |
| unaligned homology bp / kb | poa_abpoa_mc | 24 | 2.09 [1.385-5.035] | 3.969 | 24 | -17.3 [-35.1--6.955] | 23 / 1 / 0 | 2.4e-07 |
| unaligned homology bp / kb | poa_abpoa_mc__all | 24 | 7.85 [3.45-17.6] | 14.6 | 24 | -8.14 [-23.3-2.18] | 16 / 8 / 0 | 0.0098 |
| unaligned homology bp / kb | unit_aware | 24 | 0.605 [0.302-1.107] | 1.066 | 24 | -21.6 [-35.7--9.025] | 24 / 0 / 0 | 1.2e-07 |
| unaligned homology bp / kb | unit_aware__all | 24 | 1.27 [0.53-1.677] | 1.541 | 24 | -20.8 [-35.3--9.045] | 21 / 3 / 0 | 1.7e-06 |
| SV pieces / path (median) | mc | 24 | 5 [3-7] | 5.917 | - | - | - | - |
| SV pieces / path (median) | mafft_fftns2 | 24 | 3 [2-5] | 3.917 | 24 | -1.5 [-3.25-0] | 16 / 2 / 6 | 0.0014 |
| SV pieces / path (median) | mafft_fftns2__all | 24 | 3 [1.75-7] | 3.875 | 24 | -2 [-3-0] | 16 / 2 / 6 | 0.00056 |
| SV pieces / path (median) | mafft_fftnsi__all | 17 | 2 [2-3] | 2.765 | 17 | -2 [-3-0] | 12 / 1 / 4 | 0.0012 |
| SV pieces / path (median) | mafft_linsi | 24 | 2 [2-3] | 2.417 | 24 | -3 [-4--1.75] | 21 / 1 / 2 | 6.7e-06 |
| SV pieces / path (median) | mafft_linsi__all | 6 | 2 [1.25-2.75] | 2 | 6 | -1 [-2.5--1] | 5 / 0 / 1 | 0.062 |
| SV pieces / path (median) | mafft_einsi | 24 | 2.5 [1.75-4] | 2.875 | 24 | -2.5 [-4--1] | 20 / 1 / 3 | 5.7e-06 |
| SV pieces / path (median) | mafft_einsi__all | 4 | 2.5 [1.75-3] | 2.25 | 4 | -1.5 [-2.25--0.75] | 3 / 0 / 1 | 0.25 |
| SV pieces / path (median) | mafft_ginsi | 24 | 2.5 [1.75-3.25] | 2.708 | 24 | -3 [-3.25--1] | 21 / 1 / 2 | 7.6e-06 |
| SV pieces / path (median) | mafft_ginsi__all | 4 | 2 [1.75-2] | 1.75 | 4 | -2 [-2.25--1.75] | 4 / 0 / 0 | 0.12 |
| SV pieces / path (median) | poa_abpoa | 24 | 3 [2-5] | 3.292 | 24 | -3 [-3.25-0] | 17 / 4 / 3 | 8.5e-05 |
| SV pieces / path (median) | poa_abpoa__all | 23 | 2 [0.5-4] | 3.087 | 23 | -3 [-4--1] | 19 / 3 / 1 | 0.00063 |
| SV pieces / path (median) | poa_spoa | 24 | 2 [1-5.25] | 3.833 | 24 | -2 [-3-0] | 17 / 2 / 5 | 0.0013 |
| SV pieces / path (median) | poa_spoa__all | 22 | 1 [0-2.75] | 2 | 22 | -3.5 [-4.75--2.25] | 19 / 2 / 1 | 1.1e-05 |
| SV pieces / path (median) | poa_abpoa_mc | 24 | 3 [2-4.25] | 3.25 | 24 | -2 [-3.25--1] | 19 / 1 / 4 | 4e-05 |
| SV pieces / path (median) | poa_abpoa_mc__all | 24 | 3 [2-5] | 3.708 | 24 | -2 [-3--1] | 19 / 4 / 1 | 0.00086 |
| SV pieces / path (median) | unit_aware | 24 | 0.5 [0-3] | 2.917 | 24 | -3 [-4.25--1] | 20 / 3 / 1 | 0.00054 |
| SV pieces / path (median) | unit_aware__all | 24 | 0 [0-2] | 3.25 | 24 | -3 [-4.25--1] | 19 / 2 / 3 | 0.0051 |
| indel bp / net length change | mc | 24 | 1.548 [1.247-2.39] | 2.044 | - | - | - | - |
| indel bp / net length change | mafft_fftns2 | 24 | 1.135 [1.083-1.403] | 1.385 | 24 | -0.341 [-0.784--0.083] | 22 / 2 / 0 | 8.3e-06 |
| indel bp / net length change | mafft_fftns2__all | 24 | 1.255 [1.143-1.554] | 1.483 | 24 | -0.228 [-0.617--0.019] | 18 / 6 / 0 | 0.00037 |
| indel bp / net length change | mafft_fftnsi__all | 17 | 1.239 [1.141-1.639] | 1.543 | 17 | -0.449 [-0.693--0.045] | 15 / 2 / 0 | 0.00021 |
| indel bp / net length change | mafft_linsi | 24 | 1.079 [1.046-1.19] | 1.192 | 24 | -0.456 [-0.756--0.171] | 24 / 0 / 0 | 1.2e-07 |
| indel bp / net length change | mafft_linsi__all | 6 | 1.106 [1.027-1.379] | 1.199 | 6 | -0.391 [-0.824--0.057] | 6 / 0 / 0 | 0.031 |
| indel bp / net length change | mafft_einsi | 24 | 1.141 [1.103-1.37] | 1.308 | 24 | -0.352 [-0.81--0.062] | 22 / 2 / 0 | 3.9e-06 |
| indel bp / net length change | mafft_einsi__all | 4 | 1.271 [1.192-1.437] | 1.358 | 4 | -0.213 [-0.417-0.002] | 3 / 1 / 0 | 0.38 |
| indel bp / net length change | mafft_ginsi | 24 | 1.069 [1.044-1.193] | 1.186 | 24 | -0.459 [-0.787--0.175] | 24 / 0 / 0 | 1.2e-07 |
| indel bp / net length change | mafft_ginsi__all | 4 | 1.043 [1.009-1.176] | 1.143 | 4 | -0.391 [-0.699--0.11] | 4 / 0 / 0 | 0.12 |
| indel bp / net length change | poa_abpoa | 24 | 1.058 [1.03-1.176] | 1.13 | 24 | -0.484 [-1.082--0.162] | 24 / 0 / 0 | 1.2e-07 |
| indel bp / net length change | poa_abpoa__all | 23 | 1.101 [1.041-1.319] | 1.2 | 23 | -0.49 [-0.987--0.152] | 21 / 2 / 0 | 1.2e-06 |
| indel bp / net length change | poa_spoa | 24 | 1.031 [1.014-1.13] | 1.098 | 24 | -0.523 [-1.188--0.171] | 24 / 0 / 0 | 1.2e-07 |
| indel bp / net length change | poa_spoa__all | 22 | 1.082 [1.037-1.262] | 1.202 | 22 | -0.55 [-1.02--0.153] | 22 / 0 / 0 | 4.8e-07 |
| indel bp / net length change | poa_abpoa_mc | 24 | 1.199 [1.036-1.34] | 1.285 | 24 | -0.356 [-0.684--0.059] | 23 / 1 / 0 | 2.3e-06 |
| indel bp / net length change | poa_abpoa_mc__all | 24 | 1.281 [1.063-1.574] | 1.391 | 24 | -0.266 [-0.567--0.057] | 21 / 3 / 0 | 1.3e-05 |
| indel bp / net length change | unit_aware | 24 | 1.024 [1.008-1.083] | 1.062 | 24 | -0.529 [-1.246--0.189] | 24 / 0 / 0 | 1.2e-07 |
| indel bp / net length change | unit_aware__all | 24 | 1.03 [1.013-1.082] | 1.076 | 24 | -0.517 [-1.204--0.179] | 24 / 0 / 0 | 1.2e-07 |
| k-mer extra positions (frac) | mc | 24 | 0.435 [0.37-0.593] | 0.467 | - | - | - | - |
| k-mer extra positions (frac) | mafft_fftns2 | 24 | 0.363 [0.262-0.516] | 0.41 | 24 | -0.06 [-0.102--0.013] | 20 / 4 / 0 | 0.0014 |
| k-mer extra positions (frac) | mafft_fftns2__all | 24 | 0.418 [0.324-0.597] | 0.471 | 24 | 0.013 [-0.058-0.061] | 10 / 14 / 0 | 0.64 |
| k-mer extra positions (frac) | mafft_fftnsi__all | 17 | 0.369 [0.264-0.532] | 0.421 | 17 | 0.012 [-0.058-0.045] | 7 / 10 / 0 | 0.96 |
| k-mer extra positions (frac) | mafft_linsi | 24 | 0.299 [0.217-0.444] | 0.357 | 24 | -0.101 [-0.161--0.05] | 23 / 1 / 0 | 3.6e-07 |
| k-mer extra positions (frac) | mafft_linsi__all | 6 | 0.312 [0.273-0.416] | 0.32 | 6 | -0.114 [-0.135--0.056] | 6 / 0 / 0 | 0.031 |
| k-mer extra positions (frac) | mafft_einsi | 24 | 0.299 [0.242-0.437] | 0.356 | 24 | -0.12 [-0.14--0.058] | 23 / 1 / 0 | 3.6e-07 |
| k-mer extra positions (frac) | mafft_einsi__all | 4 | 0.396 [0.314-0.478] | 0.396 | 4 | -0.045 [-0.092--0.003] | 3 / 1 / 0 | 0.38 |
| k-mer extra positions (frac) | mafft_ginsi | 24 | 0.304 [0.223-0.43] | 0.356 | 24 | -0.104 [-0.161--0.047] | 22 / 2 / 0 | 6e-07 |
| k-mer extra positions (frac) | mafft_ginsi__all | 4 | 0.349 [0.213-0.456] | 0.32 | 4 | -0.139 [-0.174--0.09] | 4 / 0 / 0 | 0.12 |
| k-mer extra positions (frac) | poa_abpoa | 24 | 0.352 [0.268-0.51] | 0.393 | 24 | -0.073 [-0.125--0.02] | 20 / 4 / 0 | 0.0016 |
| k-mer extra positions (frac) | poa_abpoa__all | 23 | 0.515 [0.41-0.678] | 0.52 | 23 | 0.056 [0.021-0.115] | 5 / 18 / 0 | 0.011 |
| k-mer extra positions (frac) | poa_spoa | 24 | 0.387 [0.29-0.584] | 0.425 | 24 | -0.026 [-0.106-0.013] | 16 / 8 / 0 | 0.029 |
| k-mer extra positions (frac) | poa_spoa__all | 22 | 0.57 [0.44-0.719] | 0.555 | 22 | 0.103 [0.048-0.142] | 3 / 19 / 0 | 0.0017 |
| k-mer extra positions (frac) | poa_abpoa_mc | 24 | 0.295 [0.212-0.468] | 0.35 | 24 | -0.094 [-0.165--0.055] | 24 / 0 / 0 | 1.2e-07 |
| k-mer extra positions (frac) | poa_abpoa_mc__all | 24 | 0.432 [0.37-0.646] | 0.472 | 24 | 0.011 [-0.013-0.042] | 9 / 15 / 0 | 0.28 |
| k-mer extra positions (frac) | unit_aware | 24 | 0.433 [0.321-0.593] | 0.452 | 24 | -0.01 [-0.083-0.052] | 13 / 11 / 0 | 0.47 |
| k-mer extra positions (frac) | unit_aware__all | 24 | 0.475 [0.361-0.643] | 0.496 | 24 | 0.026 [-0.036-0.101] | 9 / 15 / 0 | 0.21 |
| truth edits to best graph path (h1+h2) | mc | 24 | 0 [0-5] | 10 | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2 | 24 | 2.5 [0-14.2] | 11.9 | 24 | 0 [0-6.25] | 3 / 11 / 10 | 0.081 |
| truth edits to best graph path (h1+h2) | mafft_fftns2__all | 24 | 1.5 [0-9] | 6.75 | 24 | 0 [0-2] | 4 / 9 / 11 | 0.26 |
| truth edits to best graph path (h1+h2) | mafft_fftnsi__all | 17 | 1 [0-7] | 7.706 | 17 | 0 [0-3] | 2 / 7 / 8 | 0.18 |
| truth edits to best graph path (h1+h2) | mafft_linsi | 24 | 3.5 [0-22.2] | 16.9 | 24 | 1.5 [0-11.8] | 1 / 13 / 10 | 0.00085 |
| truth edits to best graph path (h1+h2) | mafft_linsi__all | 6 | 0.5 [0-6.25] | 6 | 6 | 0 [0-0] | 1 / 1 / 4 | 1 |
| truth edits to best graph path (h1+h2) | mafft_einsi | 24 | 5.5 [0-23] | 19.6 | 24 | 3.5 [0-13.5] | 1 / 15 / 8 | 9.2e-05 |
| truth edits to best graph path (h1+h2) | mafft_einsi__all | 4 | 3.5 [1.5-9.75] | 7.75 | 4 | 1 [-1.5-2.75] | 1 / 2 / 1 | 1 |
| truth edits to best graph path (h1+h2) | mafft_ginsi | 24 | 2.5 [0-15.5] | 14.9 | 24 | 1 [0-7] | 1 / 14 / 9 | 0.0012 |
| truth edits to best graph path (h1+h2) | mafft_ginsi__all | 4 | 1.5 [0-9] | 7.5 | 4 | 0 [-0.75-0.75] | 1 / 1 / 2 | 1 |
| truth edits to best graph path (h1+h2) | poa_abpoa | 24 | 2 [0-16.2] | 33.0 | 24 | 2 [0-7] | 2 / 13 / 9 | 0.0024 |
| truth edits to best graph path (h1+h2) | poa_abpoa__all | 23 | 0 [0-2] | 24.9 | 23 | 0 [-3-0] | 8 / 4 / 11 | 0.41 |
| truth edits to best graph path (h1+h2) | poa_spoa | 24 | 1.5 [0-12.2] | 36.2 | 24 | 0 [0-4.75] | 3 / 10 / 11 | 0.033 |
| truth edits to best graph path (h1+h2) | poa_spoa__all | 22 | 0 [0-1.75] | 23.7 | 22 | 0 [-1.75-0] | 8 / 4 / 10 | 0.48 |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc | 24 | 5 [0-16] | 27.4 | 24 | 4 [0-9.25] | 0 / 17 / 7 | 1.5e-05 |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc__all | 24 | 0 [0-2.25] | 17.1 | 24 | 0 [-1-0.25] | 7 / 6 / 11 | 0.85 |
| truth edits to best graph path (h1+h2) | unit_aware | 24 | 0.5 [0-6] | 22.6 | 24 | 0 [0-2.25] | 5 / 9 / 10 | 0.64 |
| truth edits to best graph path (h1+h2) | unit_aware__all | 24 | 0 [0-1.25] | 17.2 | 24 | 0 [-4-0] | 8 / 4 / 12 | 0.41 |
| truth-by-graph truvari F1, raw | mc | 24 | 0.182 [0.153-0.508] | 0.301 | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftns2 | 24 | 0.297 [0.065-0.427] | 0.29 | 24 | 0 [-0.144-0.185] | 10 / 10 / 4 | 0.98 |
| truth-by-graph truvari F1, raw | mafft_fftns2__all | 24 | 0.25 [0.075-0.378] | 0.271 | 24 | 0 [-0.164-0.153] | 11 / 10 / 3 | 0.84 |
| truth-by-graph truvari F1, raw | mafft_fftnsi__all | 17 | 0.182 [0-0.286] | 0.21 | 17 | 0 [-0.226-0.068] | 6 / 8 / 3 | 0.22 |
| truth-by-graph truvari F1, raw | mafft_linsi | 23 | 0.4 [0.202-0.603] | 0.425 | 23 | 0.04 [-0.067-0.373] | 15 / 6 / 2 | 0.13 |
| truth-by-graph truvari F1, raw | mafft_linsi__all | 6 | 0.2 [0-0.433] | 0.23 | 6 | 0.109 [0-0.252] | 3 / 1 / 2 | 0.88 |
| truth-by-graph truvari F1, raw | mafft_einsi | 24 | 0.422 [0.222-0.586] | 0.438 | 24 | 0.104 [-0.016-0.283] | 15 / 6 / 3 | 0.058 |
| truth-by-graph truvari F1, raw | mafft_einsi__all | 4 | 0.5 [0.375-0.518] | 0.393 | 4 | 0.159 [-0.042-0.336] | 2 / 1 / 1 | 0.5 |
| truth-by-graph truvari F1, raw | mafft_ginsi | 23 | 0.4 [0.202-0.667] | 0.46 | 23 | 0.119 [-0.124-0.373] | 14 / 9 / 0 | 0.085 |
| truth-by-graph truvari F1, raw | mafft_ginsi__all | 4 | 0.2 [0-0.467] | 0.267 | 4 | 0.109 [-0.167-0.285] | 2 / 1 / 1 | 1 |
| truth-by-graph truvari F1, raw | poa_abpoa | 23 | 0.3 [0.174-0.445] | 0.321 | 23 | 0.008 [-0.175-0.218] | 12 / 9 / 2 | 0.75 |
| truth-by-graph truvari F1, raw | poa_abpoa__all | 23 | 0.2 [0-0.369] | 0.213 | 23 | -0.081 [-0.198-0] | 4 / 14 / 5 | 0.033 |
| truth-by-graph truvari F1, raw | poa_spoa | 23 | 0.167 [0-0.249] | 0.159 | 23 | -0.167 [-0.245-0] | 5 / 15 / 3 | 0.013 |
| truth-by-graph truvari F1, raw | poa_spoa__all | 22 | 0 [0-0.241] | 0.116 | 22 | -0.182 [-0.413--0.02] | 3 / 16 / 3 | 0.0044 |
| truth-by-graph truvari F1, raw | poa_abpoa_mc | 23 | 0.4 [0.309-0.58] | 0.435 | 23 | 0.144 [0-0.309] | 16 / 5 / 2 | 0.026 |
| truth-by-graph truvari F1, raw | poa_abpoa_mc__all | 23 | 0.37 [0.17-0.5] | 0.347 | 23 | 0 [-0.137-0.191] | 11 / 10 / 2 | 0.65 |
| truth-by-graph truvari F1, raw | unit_aware | 23 | 0.133 [0-0.309] | 0.187 | 23 | -0.133 [-0.306-0.075] | 7 / 14 / 2 | 0.066 |
| truth-by-graph truvari F1, raw | unit_aware__all | 23 | 0 [0-0.265] | 0.152 | 23 | -0.137 [-0.246-0] | 4 / 16 / 3 | 0.0041 |
| truth-by-graph truvari F1, refined | mc | 24 | 1 [1-1] | 0.924 | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftns2 | 24 | 1 [1-1] | 0.916 | 24 | 0 [0-0] | 0 / 1 / 23 | 1 |
| truth-by-graph truvari F1, refined | mafft_fftns2__all | 24 | 1 [1-1] | 0.896 | 24 | 0 [0-0] | 2 / 2 / 20 | 0.88 |
| truth-by-graph truvari F1, refined | mafft_fftnsi__all | 17 | 1 [1-1] | 0.817 | 17 | 0 [0-0] | 1 / 3 / 13 | 0.5 |
| truth-by-graph truvari F1, refined | mafft_linsi | 23 | 1 [1-1] | 0.957 | 23 | 0 [0-0] | 2 / 1 / 20 | 1 |
| truth-by-graph truvari F1, refined | mafft_linsi__all | 6 | 1 [1-1] | 1 | 6 | 0 [0-0] | 0 / 0 / 6 | - |
| truth-by-graph truvari F1, refined | mafft_einsi | 24 | 1 [1-1] | 0.859 | 24 | 0 [0-0] | 1 / 3 / 20 | 0.38 |
| truth-by-graph truvari F1, refined | mafft_einsi__all | 4 | 1 [1-1] | 1 | 4 | 0 [0-0] | 0 / 0 / 4 | - |
| truth-by-graph truvari F1, refined | mafft_ginsi | 23 | 1 [1-1] | 0.978 | 23 | 0 [0-0] | 2 / 1 / 20 | 1 |
| truth-by-graph truvari F1, refined | mafft_ginsi__all | 4 | 1 [1-1] | 1 | 4 | 0 [0-0] | 0 / 0 / 4 | - |
| truth-by-graph truvari F1, refined | poa_abpoa | 23 | 1 [1-1] | 0.858 | 23 | 0 [0-0] | 1 / 5 / 17 | 0.062 |
| truth-by-graph truvari F1, refined | poa_abpoa__all | 23 | 1 [0.388-1] | 0.715 | 23 | 0 [-0.3-0] | 0 / 7 / 16 | 0.016 |
| truth-by-graph truvari F1, refined | poa_spoa | 23 | 1 [0.3-1] | 0.701 | 23 | 0 [-0.55-0] | 1 / 8 / 14 | 0.0078 |
| truth-by-graph truvari F1, refined | poa_spoa__all | 22 | 1 [0-1] | 0.614 | 22 | 0 [-0.875-0] | 0 / 8 / 14 | 0.0078 |
| truth-by-graph truvari F1, refined | poa_abpoa_mc | 23 | 1 [1-1] | 0.932 | 23 | 0 [0-0] | 2 / 2 / 19 | 0.62 |
| truth-by-graph truvari F1, refined | poa_abpoa_mc__all | 23 | 1 [1-1] | 0.853 | 23 | 0 [0-0] | 1 / 4 / 18 | 0.19 |
| truth-by-graph truvari F1, refined | unit_aware | 23 | 1 [0.2-1] | 0.713 | 23 | 0 [-0.55-0] | 0 / 7 / 16 | 0.016 |
| truth-by-graph truvari F1, refined | unit_aware__all | 23 | 1 [0.083-1] | 0.667 | 23 | 0 [-0.667-0] | 1 / 9 / 13 | 0.0059 |
| reads: redundant (cross-walk) only | mc | 24 | 0.295 [0.212-0.357] | 0.289 | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftns2 | 24 | 0.251 [0.128-0.305] | 0.229 | 24 | -0.056 [-0.091-0.003] | 17 / 6 / 1 | 0.0089 |
| reads: redundant (cross-walk) only | mafft_fftns2__all | 24 | 0.263 [0.15-0.354] | 0.258 | 24 | -0.019 [-0.109-0.026] | 16 / 8 / 0 | 0.14 |
| reads: redundant (cross-walk) only | mafft_fftnsi__all | 17 | 0.29 [0.157-0.355] | 0.265 | 17 | -0.019 [-0.065-0.009] | 12 / 5 / 0 | 0.16 |
| reads: redundant (cross-walk) only | mafft_linsi | 24 | 0.146 [0.098-0.236] | 0.168 | 24 | -0.128 [-0.174--0.043] | 23 / 1 / 0 | 3e-06 |
| reads: redundant (cross-walk) only | mafft_linsi__all | 6 | 0.129 [0.099-0.148] | 0.126 | 6 | -0.142 [-0.202--0.054] | 6 / 0 / 0 | 0.031 |
| reads: redundant (cross-walk) only | mafft_einsi | 24 | 0.171 [0.111-0.224] | 0.171 | 24 | -0.109 [-0.164--0.061] | 22 / 2 / 0 | 6e-07 |
| reads: redundant (cross-walk) only | mafft_einsi__all | 4 | 0.152 [0.128-0.186] | 0.161 | 4 | -0.134 [-0.199--0.062] | 4 / 0 / 0 | 0.12 |
| reads: redundant (cross-walk) only | mafft_ginsi | 24 | 0.141 [0.098-0.236] | 0.166 | 24 | -0.13 [-0.158--0.067] | 22 / 2 / 0 | 2.3e-06 |
| reads: redundant (cross-walk) only | mafft_ginsi__all | 4 | 0.136 [0.108-0.159] | 0.131 | 4 | -0.184 [-0.21--0.131] | 4 / 0 / 0 | 0.12 |
| reads: redundant (cross-walk) only | poa_abpoa | 24 | 0.225 [0.154-0.303] | 0.238 | 24 | -0.046 [-0.103-0.004] | 17 / 7 / 0 | 0.013 |
| reads: redundant (cross-walk) only | poa_abpoa__all | 23 | 0.331 [0.294-0.422] | 0.352 | 23 | 0.052 [-0.019-0.079] | 8 / 15 / 0 | 0.01 |
| reads: redundant (cross-walk) only | poa_spoa | 24 | 0.275 [0.211-0.337] | 0.262 | 24 | -0.015 [-0.075-0.008] | 16 / 7 / 1 | 0.098 |
| reads: redundant (cross-walk) only | poa_spoa__all | 22 | 0.348 [0.291-0.449] | 0.356 | 22 | 0.052 [-0.017-0.116] | 7 / 15 / 0 | 0.054 |
| reads: redundant (cross-walk) only | poa_abpoa_mc | 24 | 0.195 [0.124-0.24] | 0.187 | 24 | -0.102 [-0.141--0.038] | 23 / 1 / 0 | 6e-07 |
| reads: redundant (cross-walk) only | poa_abpoa_mc__all | 24 | 0.268 [0.187-0.341] | 0.258 | 24 | -0.019 [-0.045-0.004] | 17 / 7 / 0 | 0.023 |
| reads: redundant (cross-walk) only | unit_aware | 24 | 0.261 [0.212-0.33] | 0.265 | 24 | -0.019 [-0.064-0.018] | 14 / 9 / 1 | 0.17 |
| reads: redundant (cross-walk) only | unit_aware__all | 24 | 0.278 [0.23-0.348] | 0.284 | 24 | -0.004 [-0.066-0.037] | 12 / 11 / 1 | 0.48 |
| reads: placements / max per path | mc | 24 | 2.108 [1.76-2.379] | 2.07 | - | - | - | - |
| reads: placements / max per path | mafft_fftns2 | 24 | 1.961 [1.49-2.271] | 1.946 | 24 | -0.174 [-0.463-0.17] | 16 / 8 / 0 | 0.14 |
| reads: placements / max per path | mafft_fftns2__all | 24 | 2.305 [1.619-2.604] | 2.262 | 24 | 0.059 [-0.213-0.53] | 12 / 12 / 0 | 0.24 |
| reads: placements / max per path | mafft_fftnsi__all | 17 | 2.092 [1.494-2.815] | 2.15 | 17 | 0.146 [-0.116-0.537] | 7 / 10 / 0 | 0.17 |
| reads: placements / max per path | mafft_linsi | 24 | 1.544 [1.317-1.739] | 1.538 | 24 | -0.54 [-0.736--0.208] | 23 / 1 / 0 | 8.3e-07 |
| reads: placements / max per path | mafft_linsi__all | 6 | 1.495 [1.295-1.713] | 1.505 | 6 | -0.394 [-0.747--0.142] | 6 / 0 / 0 | 0.031 |
| reads: placements / max per path | mafft_einsi | 24 | 1.445 [1.31-1.761] | 1.513 | 24 | -0.501 [-0.816--0.171] | 24 / 0 / 0 | 1.2e-07 |
| reads: placements / max per path | mafft_einsi__all | 4 | 1.501 [1.469-1.56] | 1.528 | 4 | -0.597 [-0.931--0.198] | 3 / 1 / 0 | 0.25 |
| reads: placements / max per path | mafft_ginsi | 24 | 1.542 [1.326-1.724] | 1.537 | 24 | -0.558 [-0.78--0.2] | 22 / 2 / 0 | 6e-07 |
| reads: placements / max per path | mafft_ginsi__all | 4 | 1.475 [1.289-1.672] | 1.486 | 4 | -0.669 [-0.859--0.385] | 4 / 0 / 0 | 0.12 |
| reads: placements / max per path | poa_abpoa | 24 | 1.851 [1.613-1.973] | 1.792 | 24 | -0.247 [-0.543--0.067] | 20 / 4 / 0 | 0.00043 |
| reads: placements / max per path | poa_abpoa__all | 23 | 2.34 [2.109-2.75] | 2.485 | 23 | 0.385 [0.022-0.688] | 6 / 17 / 0 | 0.00055 |
| reads: placements / max per path | poa_spoa | 24 | 1.974 [1.705-2.12] | 1.915 | 24 | -0.107 [-0.36-0.079] | 15 / 9 / 0 | 0.056 |
| reads: placements / max per path | poa_spoa__all | 22 | 2.486 [2.286-2.978] | 2.661 | 22 | 0.45 [0.21-1.004] | 2 / 20 / 0 | 0.00014 |
| reads: placements / max per path | poa_abpoa_mc | 24 | 1.623 [1.354-1.764] | 1.596 | 24 | -0.419 [-0.65--0.229] | 24 / 0 / 0 | 1.2e-07 |
| reads: placements / max per path | poa_abpoa_mc__all | 24 | 2.185 [1.787-2.244] | 2.037 | 24 | -0.049 [-0.11-0.105] | 15 / 9 / 0 | 0.64 |
| reads: placements / max per path | unit_aware | 24 | 2.024 [1.69-2.128] | 1.935 | 24 | -0.075 [-0.281-0.093] | 15 / 9 / 0 | 0.11 |
| reads: placements / max per path | unit_aware__all | 24 | 2.077 [1.87-2.404] | 2.144 | 24 | 0.038 [-0.157-0.31] | 11 / 13 / 0 | 0.36 |

| method | regions scored (raw / refined / phab) | raw TP-base / FN / TP-comp / FP | raw F1 | refined TP-base / FN / TP-comp / FP | refined F1 | phab TP-base / FN / TP-comp / FP | phab F1 |
|---|---|---|---|---|---|---|---|
| mc | 24 / 24 / 24 | 52 / 57 / 53 / 179 | 0.309 | 108 / 1 / 222 / 10 | 0.974 | 84 / 0 / 86 / 0 | 1 |
| mafft_fftns2 | 24 / 24 / 24 | 45 / 64 / 45 / 130 | 0.317 | 108 / 1 / 161 / 14 | 0.954 | 85 / 0 / 87 / 0 | 1 |
| mafft_fftns2__all | 24 / 24 / 24 | 47 / 62 / 46 / 143 | 0.311 | 107 / 2 / 179 / 10 | 0.964 | 94 / 0 / 93 / 1 | 0.995 |
| mafft_fftnsi__all | 17 / 17 / 17 | 18 / 43 / 17 / 78 | 0.223 | 56 / 5 / 89 / 6 | 0.927 | 42 / 0 / 42 / 0 | 1 |
| mafft_linsi | 24 / 24 / 24 | 42 / 67 / 43 / 60 | 0.401 | 107 / 2 / 99 / 4 | 0.971 | 87 / 0 / 86 / 0 | 1 |
| mafft_linsi__all | 6 / 6 / 6 | 7 / 11 / 7 / 18 | 0.326 | 18 / 0 / 25 / 0 | 1 | 12 / 0 / 12 / 0 | 1 |
| mafft_einsi | 24 / 24 / 24 | 52 / 57 / 51 / 83 | 0.423 | 101 / 8 / 121 / 13 | 0.915 | 79 / 1 / 78 / 0 | 0.994 |
| mafft_einsi__all | 4 / 4 / 4 | 6 / 3 / 5 / 12 | 0.408 | 9 / 0 / 17 / 0 | 1 | 9 / 0 / 9 / 0 | 1 |
| mafft_ginsi | 24 / 24 / 24 | 46 / 63 / 45 / 72 | 0.402 | 106 / 3 / 116 / 1 | 0.982 | 84 / 0 / 83 / 0 | 1 |
| mafft_ginsi__all | 4 / 4 / 4 | 4 / 5 / 4 / 13 | 0.308 | 9 / 0 / 17 / 0 | 1 | 8 / 0 / 9 / 0 | 1 |
| poa_abpoa | 24 / 24 / 24 | 39 / 70 / 39 / 89 | 0.329 | 90 / 19 / 111 / 17 | 0.846 | 73 / 1 / 72 / 2 | 0.98 |
| poa_abpoa__all | 23 / 23 / 23 | 30 / 77 / 31 / 113 | 0.244 | 80 / 27 / 125 / 19 | 0.803 | 70 / 1 / 70 / 1 | 0.986 |
| poa_spoa | 24 / 24 / 24 | 26 / 83 / 27 / 140 | 0.193 | 84 / 25 / 149 / 18 | 0.827 | 85 / 2 / 87 / 0 | 0.988 |
| poa_spoa__all | 22 / 22 / 22 | 16 / 84 / 17 / 88 | 0.161 | 90 / 10 / 103 / 2 | 0.939 | 77 / 0 / 79 / 0 | 1 |
| poa_abpoa_mc | 24 / 24 / 24 | 50 / 59 / 50 / 73 | 0.431 | 102 / 7 / 110 / 13 | 0.915 | 83 / 0 / 84 / 1 | 0.994 |
| poa_abpoa_mc__all | 24 / 24 / 24 | 51 / 58 / 51 / 106 | 0.383 | 99 / 10 / 137 / 20 | 0.89 | 78 / 1 / 79 / 0 | 0.994 |
| unit_aware | 24 / 24 / 24 | 28 / 81 / 28 / 139 | 0.203 | 93 / 16 / 157 / 10 | 0.895 | 81 / 0 / 81 / 0 | 1 |
| unit_aware__all | 24 / 24 / 24 | 25 / 84 / 25 / 140 | 0.182 | 94 / 15 / 157 / 8 | 0.905 | 81 / 0 / 82 / 0 | 1 |

## control_nontr_sv (20 regions)

| method | n | nodes / CHM13 kb | cost/opt, all pairs | excess edits / kb | affine cost/opt, all | unaligned homology bp / kb | SV pieces / path (median) | indel bp / net length change | k-mer extra positions (frac) | truth edits to best graph path (h1+h2) | truth-by-graph truvari F1, raw | truth-by-graph truvari F1, refined | reads: redundant (cross-walk) only | reads: placements / max per path | regions with cost/opt <= 1.1 | pooled truvari F1 raw / refined |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc | 20 | 48.7 | 1 | 0.005 | 1 | 0 | 1 | 1 | 0.007 | 0 | 1 | 1 | 0.023 | 1.047 | 20/20 | 1 / 1 |
| mafft_fftns2 | 20 | 17.3 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.01 | 1.01 | 20/20 | 1 / 1 |
| mafft_fftns2__all | 20 | 17.3 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.003 | 1.003 | 20/20 | 1 / 1 |
| mafft_fftnsi__all | 20 | 17.3 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.003 | 1.003 | 20/20 | 1 / 1 |
| mafft_linsi | 20 | 22.9 | 1 | 0.01 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.026 | 1.047 | 19/20 | 0.976 / 1 |
| mafft_linsi__all | 20 | 22.4 | 1 | 0.01 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.026 | 1.047 | 19/20 | 0.976 / 1 |
| mafft_einsi | 20 | 17.7 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.01 | 1.01 | 20/20 | 1 / 1 |
| mafft_einsi__all | 20 | 17.7 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.01 | 1.01 | 20/20 | 1 / 1 |
| mafft_ginsi | 20 | 22.9 | 1 | 0.01 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.026 | 1.047 | 19/20 | 0.976 / 1 |
| mafft_ginsi__all | 20 | 22.4 | 1 | 0.01 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.026 | 1.047 | 19/20 | 0.976 / 1 |
| poa_abpoa | 20 | 17.7 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.01 | 1.01 | 20/20 | 1 / 1 |
| poa_abpoa__all | 20 | 17.7 | 1 | 0.005 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.01 | 1.01 | 20/20 | 1 / 1 |
| poa_spoa | 20 | 17.7 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.01 | 1.01 | 20/20 | 1 / 1 |
| poa_spoa__all | 20 | 17.7 | 1 | 0.005 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.01 | 1.01 | 20/20 | 1 / 1 |
| poa_abpoa_mc | 20 | 17.7 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.01 | 1.01 | 20/20 | 1 / 1 |
| poa_abpoa_mc__all | 20 | 17.7 | 1 | 0.01 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.01 | 1.01 | 20/20 | 1 / 1 |
| unit_aware | 20 | 22.9 | 1 | 0.01 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.026 | 1.047 | 19/20 | 0.976 / 1 |
| unit_aware__all | 20 | 22.4 | 1 | 0.01 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 0.026 | 1.047 | 19/20 | 0.976 / 1 |

| metric | method | n | median [q1-q3] | mean | paired n | median diff [q1-q3] | better / worse / equal | p |
|---|---|---|---|---|---|---|---|---|
| nodes / CHM13 kb | mc | 20 | 48.7 [33.7-59.5] | 55.9 | - | - | - | - |
| nodes / CHM13 kb | mafft_fftns2 | 20 | 17.3 [9.015-28.1] | 27.4 | 20 | -25.8 [-31.6--22.6] | 20 / 0 / 0 | 1.9e-06 |
| nodes / CHM13 kb | mafft_fftns2__all | 20 | 17.3 [9.015-29.3] | 27.6 | 20 | -25.8 [-31.6--21.3] | 20 / 0 / 0 | 1.9e-06 |
| nodes / CHM13 kb | mafft_fftnsi__all | 20 | 17.3 [9.015-29.2] | 27.5 | 20 | -25.8 [-31.6--22.6] | 20 / 0 / 0 | 1.9e-06 |
| nodes / CHM13 kb | mafft_linsi | 20 | 22.9 [10.7-43.9] | 54.2 | 20 | -23.8 [-28.6--14.7] | 17 / 3 / 0 | 0.074 |
| nodes / CHM13 kb | mafft_linsi__all | 20 | 22.4 [10.7-43.9] | 54.2 | 20 | -23.8 [-28.6--14.7] | 17 / 3 / 0 | 0.074 |
| nodes / CHM13 kb | mafft_einsi | 20 | 17.7 [9.015-28.5] | 27.7 | 20 | -25.8 [-31.4--22.2] | 20 / 0 / 0 | 1.9e-06 |
| nodes / CHM13 kb | mafft_einsi__all | 20 | 17.7 [9.015-28.5] | 27.7 | 20 | -25.8 [-31.4--21.8] | 20 / 0 / 0 | 1.9e-06 |
| nodes / CHM13 kb | mafft_ginsi | 20 | 22.9 [10.7-43.9] | 54.1 | 20 | -23.8 [-28.6--14.7] | 17 / 3 / 0 | 0.074 |
| nodes / CHM13 kb | mafft_ginsi__all | 20 | 22.4 [10.7-43.9] | 54.1 | 20 | -23.8 [-28.6--14.7] | 17 / 3 / 0 | 0.074 |
| nodes / CHM13 kb | poa_abpoa | 20 | 17.7 [9.015-28.5] | 27.7 | 20 | -25.8 [-31.4--22.2] | 20 / 0 / 0 | 1.9e-06 |
| nodes / CHM13 kb | poa_abpoa__all | 20 | 17.7 [9.015-29.6] | 27.8 | 20 | -25.8 [-31.4--22.2] | 20 / 0 / 0 | 1.9e-06 |
| nodes / CHM13 kb | poa_spoa | 20 | 17.7 [9.015-28.5] | 27.4 | 20 | -26.5 [-31.4--22.2] | 20 / 0 / 0 | 1.9e-06 |
| nodes / CHM13 kb | poa_spoa__all | 20 | 17.7 [9.015-29.6] | 27.6 | 20 | -26.5 [-31.4--22.2] | 20 / 0 / 0 | 1.9e-06 |
| nodes / CHM13 kb | poa_abpoa_mc | 20 | 17.7 [9.015-28.5] | 27.7 | 20 | -25.8 [-31.4--22.2] | 20 / 0 / 0 | 1.9e-06 |
| nodes / CHM13 kb | poa_abpoa_mc__all | 20 | 17.7 [9.015-30.0] | 27.9 | 20 | -25.8 [-31.4--22.2] | 20 / 0 / 0 | 1.9e-06 |
| nodes / CHM13 kb | unit_aware | 20 | 22.9 [10.7-43.9] | 54.2 | 20 | -23.8 [-28.6--14.7] | 17 / 3 / 0 | 0.074 |
| nodes / CHM13 kb | unit_aware__all | 20 | 22.4 [10.7-43.9] | 54.2 | 20 | -23.8 [-28.6--14.7] | 17 / 3 / 0 | 0.074 |
| cost/opt, all pairs | mc | 20 | 1 [1-1.002] | 1.002 | - | - | - | - |
| cost/opt, all pairs | mafft_fftns2 | 20 | 1 [1-1.002] | 1.002 | 20 | 0 [0-0] | 2 / 1 / 17 | 1 |
| cost/opt, all pairs | mafft_fftns2__all | 20 | 1 [1-1.003] | 1.002 | 20 | 0 [0-0] | 2 / 2 / 16 | 0.62 |
| cost/opt, all pairs | mafft_fftnsi__all | 20 | 1 [1-1.002] | 1.002 | 20 | 0 [0-0] | 2 / 1 / 17 | 1 |
| cost/opt, all pairs | mafft_linsi | 20 | 1 [1-1.004] | 1.019 | 20 | 0 [0-0] | 2 / 5 / 13 | 0.094 |
| cost/opt, all pairs | mafft_linsi__all | 20 | 1 [1-1.004] | 1.019 | 20 | 0 [0-0] | 2 / 5 / 13 | 0.094 |
| cost/opt, all pairs | mafft_einsi | 20 | 1 [1-1.002] | 1.002 | 20 | 0 [0-0] | 2 / 1 / 17 | 1 |
| cost/opt, all pairs | mafft_einsi__all | 20 | 1 [1-1.002] | 1.002 | 20 | 0 [0-0] | 2 / 1 / 17 | 1 |
| cost/opt, all pairs | mafft_ginsi | 20 | 1 [1-1.004] | 1.019 | 20 | 0 [0-0] | 2 / 5 / 13 | 0.094 |
| cost/opt, all pairs | mafft_ginsi__all | 20 | 1 [1-1.004] | 1.019 | 20 | 0 [0-0] | 2 / 5 / 13 | 0.094 |
| cost/opt, all pairs | poa_abpoa | 20 | 1 [1-1.002] | 1.002 | 20 | 0 [0-0] | 2 / 2 / 16 | 0.75 |
| cost/opt, all pairs | poa_abpoa__all | 20 | 1 [1-1.002] | 1.002 | 20 | 0 [0-0] | 1 / 2 / 17 | 0.75 |
| cost/opt, all pairs | poa_spoa | 20 | 1 [1-1.002] | 1.002 | 20 | 0 [0-0] | 2 / 0 / 18 | 0.5 |
| cost/opt, all pairs | poa_spoa__all | 20 | 1 [1-1.002] | 1.002 | 20 | 0 [0-0] | 1 / 1 / 18 | 1 |
| cost/opt, all pairs | poa_abpoa_mc | 20 | 1 [1-1.002] | 1.002 | 20 | 0 [0-0] | 2 / 2 / 16 | 0.75 |
| cost/opt, all pairs | poa_abpoa_mc__all | 20 | 1 [1-1.003] | 1.002 | 20 | 0 [0-0] | 1 / 3 / 16 | 0.38 |
| cost/opt, all pairs | unit_aware | 20 | 1 [1-1.004] | 1.019 | 20 | 0 [0-0] | 2 / 5 / 13 | 0.094 |
| cost/opt, all pairs | unit_aware__all | 20 | 1 [1-1.004] | 1.019 | 20 | 0 [0-0] | 2 / 5 / 13 | 0.094 |
| excess edits / kb | mc | 20 | 0.005 [0-0.6] | 0.442 | - | - | - | - |
| excess edits / kb | mafft_fftns2 | 20 | 0 [0-0.6] | 0.513 | 20 | 0 [0-0] | 2 / 1 / 17 | 1 |
| excess edits / kb | mafft_fftns2__all | 20 | 0 [0-0.6] | 0.536 | 20 | 0 [0-0] | 2 / 2 / 16 | 0.62 |
| excess edits / kb | mafft_fftnsi__all | 20 | 0 [0-0.6] | 0.513 | 20 | 0 [0-0] | 2 / 1 / 17 | 1 |
| excess edits / kb | mafft_linsi | 20 | 0.01 [0-1.125] | 10.7 | 20 | 0 [0-0.015] | 2 / 5 / 13 | 0.078 |
| excess edits / kb | mafft_linsi__all | 20 | 0.01 [0-1.125] | 10.5 | 20 | 0 [0-0.015] | 2 / 5 / 13 | 0.078 |
| excess edits / kb | mafft_einsi | 20 | 0 [0-0.6] | 0.513 | 20 | 0 [0-0] | 2 / 1 / 17 | 1 |
| excess edits / kb | mafft_einsi__all | 20 | 0 [0-0.6] | 0.513 | 20 | 0 [0-0] | 2 / 1 / 17 | 1 |
| excess edits / kb | mafft_ginsi | 20 | 0.01 [0-1.125] | 10.7 | 20 | 0 [0-0.015] | 2 / 5 / 13 | 0.078 |
| excess edits / kb | mafft_ginsi__all | 20 | 0.01 [0-1.125] | 10.5 | 20 | 0 [0-0.015] | 2 / 5 / 13 | 0.078 |
| excess edits / kb | poa_abpoa | 20 | 0 [0-0.6] | 0.516 | 20 | 0 [0-0] | 2 / 2 / 16 | 0.62 |
| excess edits / kb | poa_abpoa__all | 20 | 0.005 [0-0.6] | 0.517 | 20 | 0 [0-0] | 1 / 2 / 17 | 0.5 |
| excess edits / kb | poa_spoa | 20 | 0 [0-0.6] | 0.44 | 20 | 0 [0-0] | 2 / 0 / 18 | 0.5 |
| excess edits / kb | poa_spoa__all | 20 | 0.005 [0-0.6] | 0.444 | 20 | 0 [0-0] | 1 / 1 / 18 | 1 |
| excess edits / kb | poa_abpoa_mc | 20 | 0 [0-0.6] | 0.516 | 20 | 0 [0-0] | 2 / 2 / 16 | 0.62 |
| excess edits / kb | poa_abpoa_mc__all | 20 | 0.01 [0-0.795] | 0.554 | 20 | 0 [0-0] | 1 / 3 / 16 | 0.25 |
| excess edits / kb | unit_aware | 20 | 0.01 [0-1.125] | 10.7 | 20 | 0 [0-0.015] | 2 / 5 / 13 | 0.078 |
| excess edits / kb | unit_aware__all | 20 | 0.01 [0-1.125] | 10.5 | 20 | 0 [0-0.015] | 2 / 5 / 13 | 0.078 |
| affine cost/opt, all | mc | 20 | 1 [1-1] | 1 | - | - | - | - |
| affine cost/opt, all | mafft_fftns2 | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 3 / 0 / 17 | 0.25 |
| affine cost/opt, all | mafft_fftns2__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 3 / 1 / 16 | 0.88 |
| affine cost/opt, all | mafft_fftnsi__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 3 / 0 / 17 | 0.25 |
| affine cost/opt, all | mafft_linsi | 20 | 1 [1-1] | 1.033 | 20 | 0 [0-0] | 3 / 4 / 13 | 0.41 |
| affine cost/opt, all | mafft_linsi__all | 20 | 1 [1-1] | 1.032 | 20 | 0 [0-0] | 3 / 4 / 13 | 0.41 |
| affine cost/opt, all | mafft_einsi | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 3 / 0 / 17 | 0.25 |
| affine cost/opt, all | mafft_einsi__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 3 / 0 / 17 | 0.25 |
| affine cost/opt, all | mafft_ginsi | 20 | 1 [1-1] | 1.033 | 20 | 0 [0-0] | 3 / 4 / 13 | 0.41 |
| affine cost/opt, all | mafft_ginsi__all | 20 | 1 [1-1] | 1.032 | 20 | 0 [0-0] | 3 / 4 / 13 | 0.41 |
| affine cost/opt, all | poa_abpoa | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 3 / 1 / 16 | 0.38 |
| affine cost/opt, all | poa_abpoa__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 2 / 1 / 17 | 0.75 |
| affine cost/opt, all | poa_spoa | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 2 / 1 / 17 | 1 |
| affine cost/opt, all | poa_spoa__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 1 / 2 / 17 | 0.75 |
| affine cost/opt, all | poa_abpoa_mc | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 3 / 1 / 16 | 0.38 |
| affine cost/opt, all | poa_abpoa_mc__all | 20 | 1 [1-1] | 1.001 | 20 | 0 [0-0] | 2 / 2 / 16 | 1 |
| affine cost/opt, all | unit_aware | 20 | 1 [1-1] | 1.033 | 20 | 0 [0-0] | 3 / 4 / 13 | 0.41 |
| affine cost/opt, all | unit_aware__all | 20 | 1 [1-1] | 1.032 | 20 | 0 [0-0] | 3 / 4 / 13 | 0.41 |
| unaligned homology bp / kb | mc | 20 | 0 [0-0] | 0.367 | - | - | - | - |
| unaligned homology bp / kb | mafft_fftns2 | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 1 / 0 / 19 | 1 |
| unaligned homology bp / kb | mafft_fftns2__all | 20 | 0 [0-0] | 0.037 | 20 | 0 [0-0] | 1 / 1 / 18 | 1 |
| unaligned homology bp / kb | mafft_fftnsi__all | 20 | 0 [0-0] | 0.037 | 20 | 0 [0-0] | 1 / 1 / 18 | 1 |
| unaligned homology bp / kb | mafft_linsi | 20 | 0 [0-0] | 0.554 | 20 | 0 [0-0] | 1 / 3 / 16 | 0.88 |
| unaligned homology bp / kb | mafft_linsi__all | 20 | 0 [0-0] | 0.574 | 20 | 0 [0-0] | 1 / 4 / 15 | 0.62 |
| unaligned homology bp / kb | mafft_einsi | 20 | 0 [0-0] | 0.001 | 20 | 0 [0-0] | 1 / 1 / 18 | 1 |
| unaligned homology bp / kb | mafft_einsi__all | 20 | 0 [0-0] | 0.021 | 20 | 0 [0-0] | 1 / 2 / 17 | 1 |
| unaligned homology bp / kb | mafft_ginsi | 20 | 0 [0-0] | 0.646 | 20 | 0 [0-0] | 1 / 3 / 16 | 0.88 |
| unaligned homology bp / kb | mafft_ginsi__all | 20 | 0 [0-0] | 0.553 | 20 | 0 [0-0] | 1 / 4 / 15 | 0.62 |
| unaligned homology bp / kb | poa_abpoa | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 1 / 0 / 19 | 1 |
| unaligned homology bp / kb | poa_abpoa__all | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 1 / 0 / 19 | 1 |
| unaligned homology bp / kb | poa_spoa | 20 | 0 [0-0] | 0.002 | 20 | 0 [0-0] | 1 / 0 / 19 | 1 |
| unaligned homology bp / kb | poa_spoa__all | 20 | 0 [0-0] | 0.002 | 20 | 0 [0-0] | 1 / 0 / 19 | 1 |
| unaligned homology bp / kb | poa_abpoa_mc | 20 | 0 [0-0] | 0.001 | 20 | 0 [0-0] | 1 / 1 / 18 | 1 |
| unaligned homology bp / kb | poa_abpoa_mc__all | 20 | 0 [0-0] | 0.018 | 20 | 0 [0-0] | 1 / 1 / 18 | 1 |
| unaligned homology bp / kb | unit_aware | 20 | 0 [0-0] | 0.554 | 20 | 0 [0-0] | 1 / 3 / 16 | 0.88 |
| unaligned homology bp / kb | unit_aware__all | 20 | 0 [0-0] | 0.574 | 20 | 0 [0-0] | 1 / 4 / 15 | 0.62 |
| SV pieces / path (median) | mc | 20 | 1 [0-1] | 0.6 | - | - | - | - |
| SV pieces / path (median) | mafft_fftns2 | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_fftns2__all | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_fftnsi__all | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_linsi | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_linsi__all | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_einsi | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_einsi__all | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_ginsi | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_ginsi__all | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | poa_abpoa | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | poa_abpoa__all | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | poa_spoa | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | poa_spoa__all | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | poa_abpoa_mc | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | poa_abpoa_mc__all | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | unit_aware | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | unit_aware__all | 20 | 1 [0-1] | 0.6 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | mc | 20 | 1 [1-1] | 1.003 | - | - | - | - |
| indel bp / net length change | mafft_fftns2 | 20 | 1 [1-1] | 1.003 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | mafft_fftns2__all | 20 | 1 [1-1] | 1.003 | 20 | 0 [0-0] | 0 / 1 / 19 | 1 |
| indel bp / net length change | mafft_fftnsi__all | 20 | 1 [1-1] | 1.003 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | mafft_linsi | 20 | 1 [1-1.001] | 1.006 | 20 | 0 [0-0] | 0 / 3 / 17 | 0.25 |
| indel bp / net length change | mafft_linsi__all | 20 | 1 [1-1.001] | 1.006 | 20 | 0 [0-0] | 0 / 3 / 17 | 0.25 |
| indel bp / net length change | mafft_einsi | 20 | 1 [1-1] | 1.003 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | mafft_einsi__all | 20 | 1 [1-1] | 1.003 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | mafft_ginsi | 20 | 1 [1-1.001] | 1.006 | 20 | 0 [0-0] | 0 / 3 / 17 | 0.25 |
| indel bp / net length change | mafft_ginsi__all | 20 | 1 [1-1.001] | 1.006 | 20 | 0 [0-0] | 0 / 3 / 17 | 0.25 |
| indel bp / net length change | poa_abpoa | 20 | 1 [1-1] | 1.003 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | poa_abpoa__all | 20 | 1 [1-1] | 1.003 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | poa_spoa | 20 | 1 [1-1] | 1.003 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | poa_spoa__all | 20 | 1 [1-1] | 1.003 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | poa_abpoa_mc | 20 | 1 [1-1] | 1.003 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | poa_abpoa_mc__all | 20 | 1 [1-1] | 1.003 | 20 | 0 [0-0] | 0 / 1 / 19 | 1 |
| indel bp / net length change | unit_aware | 20 | 1 [1-1.001] | 1.006 | 20 | 0 [0-0] | 0 / 3 / 17 | 0.25 |
| indel bp / net length change | unit_aware__all | 20 | 1 [1-1.001] | 1.006 | 20 | 0 [0-0] | 0 / 3 / 17 | 0.25 |
| k-mer extra positions (frac) | mc | 20 | 0.007 [0-0.009] | 0.006 | - | - | - | - |
| k-mer extra positions (frac) | mafft_fftns2 | 20 | 0 [0-0.007] | 0.003 | 20 | -0 [-0.007-0] | 10 / 3 / 7 | 0.027 |
| k-mer extra positions (frac) | mafft_fftns2__all | 20 | 0 [0-0.008] | 0.004 | 20 | 0 [-0.007-0] | 9 / 2 / 9 | 0.12 |
| k-mer extra positions (frac) | mafft_fftnsi__all | 20 | 0 [0-0.008] | 0.003 | 20 | 0 [-0.007-0] | 9 / 2 / 9 | 0.067 |
| k-mer extra positions (frac) | mafft_linsi | 20 | 0 [0-0] | 0.022 | 20 | -0.004 [-0.009-0] | 11 / 4 / 5 | 0.36 |
| k-mer extra positions (frac) | mafft_linsi__all | 20 | 0 [0-0] | 0.022 | 20 | -0.004 [-0.009-0] | 11 / 3 / 6 | 0.43 |
| k-mer extra positions (frac) | mafft_einsi | 20 | 0 [0-0] | 0.001 | 20 | -0.004 [-0.008-0] | 13 / 0 / 7 | 0.00024 |
| k-mer extra positions (frac) | mafft_einsi__all | 20 | 0 [0-0] | 0.001 | 20 | -0.004 [-0.008-0] | 13 / 0 / 7 | 0.00024 |
| k-mer extra positions (frac) | mafft_ginsi | 20 | 0 [0-0] | 0.022 | 20 | -0.004 [-0.009-0] | 11 / 4 / 5 | 0.36 |
| k-mer extra positions (frac) | mafft_ginsi__all | 20 | 0 [0-0] | 0.022 | 20 | -0.004 [-0.009-0] | 11 / 3 / 6 | 0.43 |
| k-mer extra positions (frac) | poa_abpoa | 20 | 0 [0-0] | 0.001 | 20 | -0.006 [-0.009-0] | 13 / 1 / 6 | 0.00037 |
| k-mer extra positions (frac) | poa_abpoa__all | 20 | 0 [0-0] | 0.001 | 20 | -0.006 [-0.009-0] | 13 / 0 / 7 | 0.00024 |
| k-mer extra positions (frac) | poa_spoa | 20 | 0 [0-0] | 0.001 | 20 | -0.004 [-0.008-0] | 13 / 1 / 6 | 0.00037 |
| k-mer extra positions (frac) | poa_spoa__all | 20 | 0 [0-0] | 0.001 | 20 | -0.006 [-0.009-0] | 13 / 0 / 7 | 0.00024 |
| k-mer extra positions (frac) | poa_abpoa_mc | 20 | 0 [0-0] | 0.001 | 20 | -0.006 [-0.009-0] | 13 / 0 / 7 | 0.00024 |
| k-mer extra positions (frac) | poa_abpoa_mc__all | 20 | 0 [0-0] | 0.001 | 20 | -0.006 [-0.009-0] | 13 / 1 / 6 | 0.0012 |
| k-mer extra positions (frac) | unit_aware | 20 | 0 [0-0] | 0.022 | 20 | -0.004 [-0.009-0] | 11 / 4 / 5 | 0.36 |
| k-mer extra positions (frac) | unit_aware__all | 20 | 0 [0-0] | 0.022 | 20 | -0.004 [-0.009-0] | 11 / 3 / 6 | 0.43 |
| truth edits to best graph path (h1+h2) | mc | 20 | 0 [0-0] | 0 | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2 | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2__all | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | mafft_fftnsi__all | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | mafft_linsi | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | mafft_linsi__all | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | mafft_einsi | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | mafft_einsi__all | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | mafft_ginsi | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | mafft_ginsi__all | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | poa_abpoa | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | poa_abpoa__all | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | poa_spoa | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | poa_spoa__all | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc__all | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | unit_aware | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth edits to best graph path (h1+h2) | unit_aware__all | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | mc | 20 | 1 [1-1] | 1 | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftns2 | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | mafft_fftns2__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | mafft_fftnsi__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | mafft_linsi | 20 | 1 [1-1] | 0.983 | 20 | 0 [0-0] | 0 / 1 / 19 | 1 |
| truth-by-graph truvari F1, raw | mafft_linsi__all | 20 | 1 [1-1] | 0.983 | 20 | 0 [0-0] | 0 / 1 / 19 | 1 |
| truth-by-graph truvari F1, raw | mafft_einsi | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | mafft_einsi__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | mafft_ginsi | 20 | 1 [1-1] | 0.983 | 20 | 0 [0-0] | 0 / 1 / 19 | 1 |
| truth-by-graph truvari F1, raw | mafft_ginsi__all | 20 | 1 [1-1] | 0.983 | 20 | 0 [0-0] | 0 / 1 / 19 | 1 |
| truth-by-graph truvari F1, raw | poa_abpoa | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | poa_abpoa__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | poa_spoa | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | poa_spoa__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | poa_abpoa_mc | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | poa_abpoa_mc__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, raw | unit_aware | 20 | 1 [1-1] | 0.983 | 20 | 0 [0-0] | 0 / 1 / 19 | 1 |
| truth-by-graph truvari F1, raw | unit_aware__all | 20 | 1 [1-1] | 0.983 | 20 | 0 [0-0] | 0 / 1 / 19 | 1 |
| truth-by-graph truvari F1, refined | mc | 20 | 1 [1-1] | 1 | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftns2 | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | mafft_fftns2__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | mafft_fftnsi__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | mafft_linsi | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | mafft_linsi__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | mafft_einsi | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | mafft_einsi__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | mafft_ginsi | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | mafft_ginsi__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | poa_abpoa | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | poa_abpoa__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | poa_spoa | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | poa_spoa__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | poa_abpoa_mc | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | poa_abpoa_mc__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | unit_aware | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| truth-by-graph truvari F1, refined | unit_aware__all | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| reads: redundant (cross-walk) only | mc | 20 | 0.023 [0-0.06] | 0.037 | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftns2 | 20 | 0.01 [0-0.044] | 0.024 | 20 | 0 [0-0] | 4 / 2 / 14 | 0.16 |
| reads: redundant (cross-walk) only | mafft_fftns2__all | 20 | 0.003 [0-0.041] | 0.022 | 20 | 0 [-0.01-0] | 6 / 1 / 13 | 0.047 |
| reads: redundant (cross-walk) only | mafft_fftnsi__all | 20 | 0.003 [0-0.041] | 0.022 | 20 | 0 [-0.01-0] | 6 / 1 / 13 | 0.047 |
| reads: redundant (cross-walk) only | mafft_linsi | 20 | 0.026 [0-0.066] | 0.034 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |
| reads: redundant (cross-walk) only | mafft_linsi__all | 20 | 0.026 [0-0.066] | 0.034 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |
| reads: redundant (cross-walk) only | mafft_einsi | 20 | 0.01 [0-0.041] | 0.025 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: redundant (cross-walk) only | mafft_einsi__all | 20 | 0.01 [0-0.051] | 0.027 | 20 | 0 [0-0] | 3 / 4 / 13 | 0.81 |
| reads: redundant (cross-walk) only | mafft_ginsi | 20 | 0.026 [0-0.066] | 0.034 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |
| reads: redundant (cross-walk) only | mafft_ginsi__all | 20 | 0.026 [0-0.066] | 0.034 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |
| reads: redundant (cross-walk) only | poa_abpoa | 20 | 0.01 [0-0.046] | 0.025 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: redundant (cross-walk) only | poa_abpoa__all | 20 | 0.01 [0-0.046] | 0.025 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: redundant (cross-walk) only | poa_spoa | 20 | 0.01 [0-0.041] | 0.025 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: redundant (cross-walk) only | poa_spoa__all | 20 | 0.01 [0-0.046] | 0.025 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: redundant (cross-walk) only | poa_abpoa_mc | 20 | 0.01 [0-0.046] | 0.025 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: redundant (cross-walk) only | poa_abpoa_mc__all | 20 | 0.01 [0-0.046] | 0.025 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: redundant (cross-walk) only | unit_aware | 20 | 0.026 [0-0.066] | 0.034 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |
| reads: redundant (cross-walk) only | unit_aware__all | 20 | 0.026 [0-0.066] | 0.034 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |
| reads: placements / max per path | mc | 20 | 1.047 [1-1.107] | 1.066 | - | - | - | - |
| reads: placements / max per path | mafft_fftns2 | 20 | 1.01 [1-1.058] | 1.053 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.22 |
| reads: placements / max per path | mafft_fftns2__all | 20 | 1.003 [1-1.053] | 1.052 | 20 | 0 [-0.01-0] | 6 / 2 / 12 | 0.078 |
| reads: placements / max per path | mafft_fftnsi__all | 20 | 1.003 [1-1.053] | 1.052 | 20 | 0 [-0.01-0] | 6 / 2 / 12 | 0.078 |
| reads: placements / max per path | mafft_linsi | 20 | 1.047 [1-1.084] | 1.064 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |
| reads: placements / max per path | mafft_linsi__all | 20 | 1.047 [1-1.084] | 1.064 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |
| reads: placements / max per path | mafft_einsi | 20 | 1.01 [1-1.065] | 1.055 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: placements / max per path | mafft_einsi__all | 20 | 1.01 [1-1.076] | 1.058 | 20 | 0 [0-0] | 3 / 4 / 13 | 0.81 |
| reads: placements / max per path | mafft_ginsi | 20 | 1.047 [1-1.084] | 1.064 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |
| reads: placements / max per path | mafft_ginsi__all | 20 | 1.047 [1-1.084] | 1.064 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |
| reads: placements / max per path | poa_abpoa | 20 | 1.01 [1-1.065] | 1.055 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: placements / max per path | poa_abpoa__all | 20 | 1.01 [1-1.065] | 1.055 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: placements / max per path | poa_spoa | 20 | 1.01 [1-1.065] | 1.055 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: placements / max per path | poa_spoa__all | 20 | 1.01 [1-1.065] | 1.055 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: placements / max per path | poa_abpoa_mc | 20 | 1.01 [1-1.065] | 1.055 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: placements / max per path | poa_abpoa_mc__all | 20 | 1.01 [1-1.065] | 1.055 | 20 | 0 [0-0] | 4 / 3 / 13 | 0.38 |
| reads: placements / max per path | unit_aware | 20 | 1.047 [1-1.084] | 1.064 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |
| reads: placements / max per path | unit_aware__all | 20 | 1.047 [1-1.084] | 1.064 | 20 | 0 [0-0.002] | 4 / 5 / 11 | 0.82 |

| method | regions scored (raw / refined / phab) | raw TP-base / FN / TP-comp / FP | raw F1 | refined TP-base / FN / TP-comp / FP | refined F1 | phab TP-base / FN / TP-comp / FP | phab F1 |
|---|---|---|---|---|---|---|---|
| mc | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| mafft_fftns2 | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| mafft_fftns2__all | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| mafft_fftnsi__all | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| mafft_linsi | 20 / 20 / 20 | 20 / 0 / 20 / 1 | 0.976 | 20 / 0 / 21 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| mafft_linsi__all | 20 / 20 / 20 | 20 / 0 / 20 / 1 | 0.976 | 20 / 0 / 21 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| mafft_einsi | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| mafft_einsi__all | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| mafft_ginsi | 20 / 20 / 20 | 20 / 0 / 20 / 1 | 0.976 | 20 / 0 / 21 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| mafft_ginsi__all | 20 / 20 / 20 | 20 / 0 / 20 / 1 | 0.976 | 20 / 0 / 21 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| poa_abpoa | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| poa_abpoa__all | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| poa_spoa | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| poa_spoa__all | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| poa_abpoa_mc | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| poa_abpoa_mc__all | 20 / 20 / 20 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| unit_aware | 20 / 20 / 20 | 20 / 0 / 20 / 1 | 0.976 | 20 / 0 / 21 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |
| unit_aware__all | 20 / 20 / 20 | 20 / 0 / 20 / 1 | 0.976 | 20 / 0 / 21 / 0 | 1 | 20 / 0 / 20 / 0 | 1 |

