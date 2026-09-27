# Candidate graphs by stratum

Results: `results/full`. Baseline: `mc`. Methods: `mc`, `mafft_fftns2`, `mafft_fftnsi`, `mafft_linsi`, `mafft_einsi`, `mafft_ginsi`, `poa_abpoa`, `poa_spoa`, `poa_abpoa_mc`, `unit_aware`.

| method | valid regions | invalid |
|---|---|---|
| mc | 149 | 0 |
| mafft_fftns2 | 140 | 0 |
| mafft_fftnsi | 105 | 0 |
| mafft_linsi | 77 | 0 |
| mafft_einsi | 65 | 0 |
| mafft_ginsi | 65 | 0 |
| poa_abpoa | 138 | 0 |
| poa_spoa | 123 | 0 |
| poa_abpoa_mc | 148 | 0 |
| unit_aware | 149 | 0 |

Cells: median [q1-q3] over regions. Paired columns compare each method with `mc` on the regions both have: median difference (method - baseline), regions better/worse/equal (by the metric's direction; '.' = no preferred direction), Wilcoxon signed-rank two-sided p.

## hotspot_vntr (40 regions)

| method | n | nodes / CHM13 kb | cost/opt, all pairs | excess edits / kb | affine cost/opt, all | unaligned homology bp / kb | SV pieces / path (median) | indel bp / net length change | k-mer extra positions (frac) | truth edits to best graph path (h1+h2) | truth-by-graph truvari F1, raw | truth-by-graph truvari F1, refined | reads: redundant (cross-walk) only | reads: placements / max per path | regions with cost/opt <= 1.1 | pooled truvari F1 raw / refined |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc | 40 | 1370 | 1.513 | 242 | 1.473 | 55.0 | 8 | 1.811 | 0.784 | - | - | - | - | - | 2/40 | - / - |
| mafft_fftns2 | 31 | 1683 | 1.361 | 139 | 1.307 | 12.2 | 7 | 1.371 | 0.815 | - | - | - | - | - | 1/31 | - / - |
| mafft_fftnsi | 9 | 1130 | 1.289 | 91.9 | 1.256 | 6.35 | 3 | 1.406 | 0.71 | - | - | - | - | - | 1/9 | - / - |
| mafft_linsi | 4 | 723 | 1.04 | 19.1 | 1.044 | 0.86 | 1 | 1.034 | 0.545 | - | - | - | - | - | 3/4 | - / - |
| mafft_einsi | 1 | 208 | 1.055 | 16.6 | 1.054 | 6.29 | 1 | 1.02 | 0.486 | - | - | - | - | - | 1/1 | - / - |
| mafft_ginsi | 1 | 258 | 1.045 | 13.4 | 1.047 | 0.64 | 1 | 1.02 | 0.358 | - | - | - | - | - | 1/1 | - / - |
| poa_abpoa | 30 | 1067 | 1.111 | 43.1 | 1.226 | 2.91 | 4 | 1.08 | 0.681 | - | - | - | - | - | 13/30 | - / - |
| poa_spoa | 20 | 1080 | 1.094 | 29.6 | 1.245 | 2.685 | 1.5 | 1.064 | 0.689 | - | - | - | - | - | 11/20 | - / - |
| poa_abpoa_mc | 39 | 1013 | 1.269 | 104 | 1.289 | 13.2 | 7 | 1.278 | 0.706 | - | - | - | - | - | 8/39 | - / - |
| unit_aware | 40 | 1449 | 1.054 | 21.1 | 1.127 | 0.58 | 4 | 1.037 | 0.708 | - | - | - | - | - | 32/40 | - / - |

| metric | method | n | median [q1-q3] | mean | paired n | median diff [q1-q3] | better / worse / equal | p |
|---|---|---|---|---|---|---|---|---|
| nodes / CHM13 kb | mc | 40 | 1370 [801-3635] | 2283 | - | - | - | - |
| nodes / CHM13 kb | mafft_fftns2 | 31 | 1683 [1058-3871] | 3001 | 31 | 659 [242-995] | 2 / 29 / 0 | 5.1e-08 |
| nodes / CHM13 kb | mafft_fftnsi | 9 | 1130 [980-1416] | 1222 | 9 | 316 [111-433] | 0 / 9 / 0 | 0.0039 |
| nodes / CHM13 kb | mafft_linsi | 4 | 723 [526-970] | 772 | 4 | -23.6 [-118-58.7] | 2 / 2 / 0 | 0.88 |
| nodes / CHM13 kb | mafft_einsi | 1 | 208 [208-208] | 208 | 1 | -45.6 [-45.6--45.6] | 1 / 0 / 0 | 1 |
| nodes / CHM13 kb | mafft_ginsi | 1 | 258 [258-258] | 258 | 1 | 5.13 [5.13-5.13] | 0 / 1 / 0 | 1 |
| nodes / CHM13 kb | poa_abpoa | 30 | 1067 [808-1768] | 1510 | 30 | -101 [-590-70.5] | 19 / 11 / 0 | 0.017 |
| nodes / CHM13 kb | poa_spoa | 20 | 1080 [800-1375] | 1379 | 20 | 62.8 [-248-196] | 9 / 11 / 0 | 0.96 |
| nodes / CHM13 kb | poa_abpoa_mc | 39 | 1013 [745-1824] | 1774 | 39 | -231 [-806-16.8] | 27 / 12 / 0 | 0.00012 |
| nodes / CHM13 kb | unit_aware | 40 | 1449 [985-2537] | 2236 | 40 | -15.0 [-356-478] | 20 / 20 / 0 | 0.87 |
| cost/opt, all pairs | mc | 40 | 1.513 [1.444-1.855] | 1.733 | - | - | - | - |
| cost/opt, all pairs | mafft_fftns2 | 31 | 1.361 [1.26-1.519] | 1.514 | 31 | -0.239 [-0.425--0.138] | 28 / 3 / 0 | 4e-06 |
| cost/opt, all pairs | mafft_fftnsi | 9 | 1.289 [1.173-1.357] | 1.429 | 9 | -0.4 [-0.633--0.266] | 9 / 0 / 0 | 0.0039 |
| cost/opt, all pairs | mafft_linsi | 4 | 1.04 [1.034-1.098] | 1.093 | 4 | -0.476 [-1.063--0.318] | 4 / 0 / 0 | 0.12 |
| cost/opt, all pairs | mafft_einsi | 1 | 1.055 [1.055-1.055] | 1.055 | 1 | -0.537 [-0.537--0.537] | 1 / 0 / 0 | 1 |
| cost/opt, all pairs | mafft_ginsi | 1 | 1.045 [1.045-1.045] | 1.045 | 1 | -0.547 [-0.547--0.547] | 1 / 0 / 0 | 1 |
| cost/opt, all pairs | poa_abpoa | 30 | 1.111 [1.064-1.294] | 1.215 | 30 | -0.459 [-0.728--0.344] | 29 / 1 / 0 | 9.3e-09 |
| cost/opt, all pairs | poa_spoa | 20 | 1.094 [1.053-1.283] | 1.167 | 20 | -0.594 [-0.813--0.39] | 20 / 0 / 0 | 1.9e-06 |
| cost/opt, all pairs | poa_abpoa_mc | 39 | 1.269 [1.12-1.515] | 1.389 | 39 | -0.276 [-0.44--0.109] | 36 / 3 / 0 | 2e-08 |
| cost/opt, all pairs | unit_aware | 40 | 1.054 [1.036-1.078] | 1.072 | 40 | -0.465 [-0.794--0.404] | 40 / 0 / 0 | 1.8e-12 |
| excess edits / kb | mc | 40 | 242 [169-328] | 239 | - | - | - | - |
| excess edits / kb | mafft_fftns2 | 31 | 139 [100-186] | 143 | 31 | -87.4 [-150--41.0] | 28 / 3 / 0 | 1e-06 |
| excess edits / kb | mafft_fftnsi | 9 | 91.9 [40.2-93.8] | 88.9 | 9 | -121 [-158--72.0] | 9 / 0 / 0 | 0.0039 |
| excess edits / kb | mafft_linsi | 4 | 19.1 [17.0-23.6] | 21.5 | 4 | -192 [-251--132] | 4 / 0 / 0 | 0.12 |
| excess edits / kb | mafft_einsi | 1 | 16.6 [16.6-16.6] | 16.6 | 1 | -162 [-162--162] | 1 / 0 / 0 | 1 |
| excess edits / kb | mafft_ginsi | 1 | 13.4 [13.4-13.4] | 13.4 | 1 | -165 [-165--165] | 1 / 0 / 0 | 1 |
| excess edits / kb | poa_abpoa | 30 | 43.1 [29.8-75.2] | 51.5 | 30 | -207 [-279--130] | 29 / 1 / 0 | 9.3e-09 |
| excess edits / kb | poa_spoa | 20 | 29.6 [22.8-47.8] | 36.0 | 20 | -193 [-289--128] | 20 / 0 / 0 | 1.9e-06 |
| excess edits / kb | poa_abpoa_mc | 39 | 104 [62.9-153] | 112 | 39 | -140 [-199--52.8] | 36 / 3 / 0 | 5e-10 |
| excess edits / kb | unit_aware | 40 | 21.1 [14.4-32.9] | 23.4 | 40 | -227 [-300--151] | 40 / 0 / 0 | 1.8e-12 |
| affine cost/opt, all | mc | 40 | 1.473 [1.396-1.774] | 1.644 | - | - | - | - |
| affine cost/opt, all | mafft_fftns2 | 31 | 1.307 [1.257-1.463] | 1.468 | 31 | -0.198 [-0.335--0.109] | 29 / 2 / 0 | 2.6e-06 |
| affine cost/opt, all | mafft_fftnsi | 9 | 1.256 [1.185-1.358] | 1.393 | 9 | -0.407 [-0.494--0.249] | 9 / 0 / 0 | 0.0039 |
| affine cost/opt, all | mafft_linsi | 4 | 1.044 [1.042-1.085] | 1.083 | 4 | -0.48 [-0.973--0.317] | 4 / 0 / 0 | 0.12 |
| affine cost/opt, all | mafft_einsi | 1 | 1.054 [1.054-1.054] | 1.054 | 1 | -0.564 [-0.564--0.564] | 1 / 0 / 0 | 1 |
| affine cost/opt, all | mafft_ginsi | 1 | 1.047 [1.047-1.047] | 1.047 | 1 | -0.572 [-0.572--0.572] | 1 / 0 / 0 | 1 |
| affine cost/opt, all | poa_abpoa | 30 | 1.226 [1.127-1.383] | 1.282 | 30 | -0.328 [-0.548--0.142] | 29 / 1 / 0 | 1e-07 |
| affine cost/opt, all | poa_spoa | 20 | 1.245 [1.127-1.37] | 1.3 | 20 | -0.405 [-0.56--0.258] | 19 / 1 / 0 | 9.5e-06 |
| affine cost/opt, all | poa_abpoa_mc | 39 | 1.289 [1.142-1.43] | 1.351 | 39 | -0.222 [-0.367--0.11] | 36 / 3 / 0 | 6.1e-08 |
| affine cost/opt, all | unit_aware | 40 | 1.127 [1.082-1.189] | 1.163 | 40 | -0.365 [-0.587--0.2] | 36 / 4 / 0 | 9.7e-10 |
| unaligned homology bp / kb | mc | 40 | 55.0 [31.8-109] | 84.4 | - | - | - | - |
| unaligned homology bp / kb | mafft_fftns2 | 31 | 12.2 [7.72-19.4] | 13.6 | 31 | -25.4 [-54.2--13.5] | 29 / 2 / 0 | 4.7e-09 |
| unaligned homology bp / kb | mafft_fftnsi | 9 | 6.35 [3.15-8.53] | 6.32 | 9 | -25.6 [-32.9--15.6] | 9 / 0 / 0 | 0.0039 |
| unaligned homology bp / kb | mafft_linsi | 4 | 0.86 [0.59-1.45] | 1.18 | 4 | -20.2 [-27.7--16.3] | 4 / 0 / 0 | 0.12 |
| unaligned homology bp / kb | mafft_einsi | 1 | 6.29 [6.29-6.29] | 6.29 | 1 | -7.52 [-7.52--7.52] | 1 / 0 / 0 | 1 |
| unaligned homology bp / kb | mafft_ginsi | 1 | 0.64 [0.64-0.64] | 0.64 | 1 | -13.2 [-13.2--13.2] | 1 / 0 / 0 | 1 |
| unaligned homology bp / kb | poa_abpoa | 30 | 2.91 [1.58-4.412] | 3.904 | 30 | -38.7 [-59.5--23.9] | 30 / 0 / 0 | 1.9e-09 |
| unaligned homology bp / kb | poa_spoa | 20 | 2.685 [1.573-4.838] | 3.366 | 20 | -32.5 [-60.8--19.6] | 20 / 0 / 0 | 1.9e-06 |
| unaligned homology bp / kb | poa_abpoa_mc | 39 | 13.2 [6.47-24.0] | 19.0 | 39 | -31.1 [-76.0--15.6] | 37 / 2 / 0 | 2.8e-08 |
| unaligned homology bp / kb | unit_aware | 40 | 0.58 [0.12-1.18] | 0.924 | 40 | -52.4 [-106--30.1] | 40 / 0 / 0 | 1.8e-12 |
| SV pieces / path (median) | mc | 40 | 8 [7-11] | 8.9 | - | - | - | - |
| SV pieces / path (median) | mafft_fftns2 | 31 | 7 [4.5-10.5] | 8.194 | 31 | -0.5 [-3-1.5] | 16 / 11 / 4 | 0.47 |
| SV pieces / path (median) | mafft_fftnsi | 9 | 3 [2-6] | 4.111 | 9 | -2 [-4--0.5] | 7 / 0 / 2 | 0.016 |
| SV pieces / path (median) | mafft_linsi | 4 | 1 [1-1.25] | 1.25 | 4 | -5.5 [-6--4.875] | 4 / 0 / 0 | 0.12 |
| SV pieces / path (median) | mafft_einsi | 1 | 1 [1-1] | 1 | 1 | -6 [-6--6] | 1 / 0 / 0 | 1 |
| SV pieces / path (median) | mafft_ginsi | 1 | 1 [1-1] | 1 | 1 | -6 [-6--6] | 1 / 0 / 0 | 1 |
| SV pieces / path (median) | poa_abpoa | 30 | 4 [2-7.75] | 5.833 | 30 | -4 [-5.75--0.25] | 22 / 5 / 3 | 0.0013 |
| SV pieces / path (median) | poa_spoa | 20 | 1.5 [1-3.5] | 3 | 20 | -5.5 [-8.25--1.75] | 17 / 2 / 1 | 0.00027 |
| SV pieces / path (median) | poa_abpoa_mc | 39 | 7 [4.5-7] | 6.987 | 39 | -2 [-5-0] | 28 / 7 / 4 | 0.0027 |
| SV pieces / path (median) | unit_aware | 40 | 4 [1.75-11.1] | 9.287 | 40 | -2 [-6.25-2.5] | 23 / 15 / 2 | 0.22 |
| indel bp / net length change | mc | 40 | 1.811 [1.432-2.061] | 2.177 | - | - | - | - |
| indel bp / net length change | mafft_fftns2 | 31 | 1.371 [1.236-1.808] | 1.861 | 31 | -0.352 [-0.663-0.008] | 22 / 9 / 0 | 0.0007 |
| indel bp / net length change | mafft_fftnsi | 9 | 1.406 [1.314-1.959] | 1.718 | 9 | -0.704 [-1.97--0.099] | 8 / 1 / 0 | 0.02 |
| indel bp / net length change | mafft_linsi | 4 | 1.034 [1.016-1.169] | 1.152 | 4 | -0.815 [-1.895--0.433] | 4 / 0 / 0 | 0.12 |
| indel bp / net length change | mafft_einsi | 1 | 1.02 [1.02-1.02] | 1.02 | 1 | -1.076 [-1.076--1.076] | 1 / 0 / 0 | 1 |
| indel bp / net length change | mafft_ginsi | 1 | 1.02 [1.02-1.02] | 1.02 | 1 | -1.075 [-1.075--1.075] | 1 / 0 / 0 | 1 |
| indel bp / net length change | poa_abpoa | 30 | 1.08 [1.019-1.184] | 1.231 | 30 | -0.751 [-1.003--0.404] | 30 / 0 / 0 | 1.9e-09 |
| indel bp / net length change | poa_spoa | 20 | 1.064 [1.006-1.175] | 1.184 | 20 | -0.909 [-1.37--0.565] | 20 / 0 / 0 | 1.9e-06 |
| indel bp / net length change | poa_abpoa_mc | 39 | 1.278 [1.064-1.674] | 1.531 | 39 | -0.366 [-0.64--0.166] | 38 / 1 / 0 | 1.8e-11 |
| indel bp / net length change | unit_aware | 40 | 1.037 [1.006-1.08] | 1.082 | 40 | -0.781 [-1.015--0.421] | 40 / 0 / 0 | 1.8e-12 |
| k-mer extra positions (frac) | mc | 40 | 0.784 [0.643-0.826] | 0.719 | - | - | - | - |
| k-mer extra positions (frac) | mafft_fftns2 | 31 | 0.815 [0.723-0.865] | 0.778 | 31 | 0.06 [0.025-0.108] | 4 / 27 / 0 | 5e-07 |
| k-mer extra positions (frac) | mafft_fftnsi | 9 | 0.71 [0.686-0.772] | 0.712 | 9 | 0.056 [0.018-0.112] | 0 / 9 / 0 | 0.0039 |
| k-mer extra positions (frac) | mafft_linsi | 4 | 0.545 [0.469-0.648] | 0.572 | 4 | -0.041 [-0.095-0.008] | 2 / 2 / 0 | 0.62 |
| k-mer extra positions (frac) | mafft_einsi | 1 | 0.486 [0.486-0.486] | 0.486 | 1 | 0.137 [0.137-0.137] | 0 / 1 / 0 | 1 |
| k-mer extra positions (frac) | mafft_ginsi | 1 | 0.358 [0.358-0.358] | 0.358 | 1 | 0.01 [0.01-0.01] | 0 / 1 / 0 | 1 |
| k-mer extra positions (frac) | poa_abpoa | 30 | 0.681 [0.582-0.803] | 0.669 | 30 | -0.031 [-0.089-0.016] | 22 / 8 / 0 | 0.045 |
| k-mer extra positions (frac) | poa_spoa | 20 | 0.689 [0.629-0.815] | 0.678 | 20 | 0.011 [-0.069-0.045] | 9 / 11 / 0 | 0.99 |
| k-mer extra positions (frac) | poa_abpoa_mc | 39 | 0.706 [0.559-0.792] | 0.669 | 39 | -0.05 [-0.099--0.014] | 34 / 5 / 0 | 5.6e-07 |
| k-mer extra positions (frac) | unit_aware | 40 | 0.708 [0.626-0.799] | 0.701 | 40 | -0.048 [-0.088-0.019] | 26 / 14 / 0 | 0.048 |
| truth edits to best graph path (h1+h2) | mc | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_linsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_einsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_abpoa | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_spoa | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | unit_aware | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_linsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_einsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_abpoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_spoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | unit_aware | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_linsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_einsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_abpoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_spoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | unit_aware | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mc | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftns2 | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftnsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_linsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_einsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_ginsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_abpoa | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_spoa | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | unit_aware | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mc | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_fftns2 | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_fftnsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_linsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_einsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_ginsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_abpoa | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_spoa | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| reads: placements / max per path | unit_aware | 0 | - | - | - | - | - | - |

| method | regions scored (raw / refined / phab) | raw TP-base / FN / TP-comp / FP | raw F1 | refined TP-base / FN / TP-comp / FP | refined F1 | phab TP-base / FN / TP-comp / FP | phab F1 |
|---|---|---|---|---|---|---|---|
| mc | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_fftns2 | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_fftnsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_linsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_einsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_ginsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_abpoa | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_spoa | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_abpoa_mc | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| unit_aware | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |

## control_vntr_matched (40 regions)

| method | n | nodes / CHM13 kb | cost/opt, all pairs | excess edits / kb | affine cost/opt, all | unaligned homology bp / kb | SV pieces / path (median) | indel bp / net length change | k-mer extra positions (frac) | truth edits to best graph path (h1+h2) | truth-by-graph truvari F1, raw | truth-by-graph truvari F1, refined | reads: redundant (cross-walk) only | reads: placements / max per path | regions with cost/opt <= 1.1 | pooled truvari F1 raw / refined |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc | 40 | 138 | 1.199 | 8.875 | 1.143 | 0.395 | 0 | 1.343 | 0.159 | - | - | - | - | - | 15/40 | - / - |
| mafft_fftns2 | 40 | 143 | 1.309 | 29.0 | 1.234 | 0.29 | 1 | 1.523 | 0.248 | - | - | - | - | - | 10/40 | - / - |
| mafft_fftnsi | 34 | 124 | 1.243 | 11.5 | 1.16 | 0.685 | 0.5 | 1.368 | 0.212 | - | - | - | - | - | 12/34 | - / - |
| mafft_linsi | 27 | 156 | 1.083 | 3.55 | 1.057 | 0.02 | 0 | 1.058 | 0.173 | - | - | - | - | - | 16/27 | - / - |
| mafft_einsi | 22 | 80.3 | 1.068 | 2.205 | 1.012 | 0.02 | 0 | 1.054 | 0.057 | - | - | - | - | - | 14/22 | - / - |
| mafft_ginsi | 22 | 110 | 1.059 | 2.01 | 1.026 | 0.02 | 0 | 1.035 | 0.098 | - | - | - | - | - | 16/22 | - / - |
| poa_abpoa | 40 | 152 | 1.114 | 7.38 | 1.125 | 0.1 | 0 | 1.115 | 0.186 | - | - | - | - | - | 17/40 | - / - |
| poa_spoa | 37 | 202 | 1.09 | 5.03 | 1.157 | 0.21 | 0 | 1.1 | 0.26 | - | - | - | - | - | 21/37 | - / - |
| poa_abpoa_mc | 40 | 107 | 1.149 | 8.665 | 1.075 | 0.105 | 0 | 1.233 | 0.119 | - | - | - | - | - | 16/40 | - / - |
| unit_aware | 40 | 187 | 1.054 | 3.245 | 1.069 | 0.015 | 0 | 1.034 | 0.185 | - | - | - | - | - | 31/40 | - / - |

| metric | method | n | median [q1-q3] | mean | paired n | median diff [q1-q3] | better / worse / equal | p |
|---|---|---|---|---|---|---|---|---|
| nodes / CHM13 kb | mc | 40 | 138 [64.3-208] | 184 | - | - | - | - |
| nodes / CHM13 kb | mafft_fftns2 | 40 | 143 [79.3-356] | 272 | 40 | 12.9 [-4.693-76.3] | 13 / 26 / 1 | 0.0021 |
| nodes / CHM13 kb | mafft_fftnsi | 34 | 124 [71.4-285] | 223 | 34 | 13.7 [-6.743-61.8] | 12 / 22 / 0 | 0.017 |
| nodes / CHM13 kb | mafft_linsi | 27 | 156 [71.0-198] | 178 | 27 | 2.5 [-9.545-42.3] | 13 / 14 / 0 | 0.22 |
| nodes / CHM13 kb | mafft_einsi | 22 | 80.3 [49.0-117] | 116 | 22 | -10.3 [-57.8--0.872] | 19 / 3 / 0 | 0.00018 |
| nodes / CHM13 kb | mafft_ginsi | 22 | 110 [57.6-192] | 157 | 22 | 0.345 [-10.2-29.0] | 11 / 11 / 0 | 0.59 |
| nodes / CHM13 kb | poa_abpoa | 40 | 152 [91.7-271] | 241 | 40 | 26.7 [-1.22-111] | 14 / 26 / 0 | 0.00071 |
| nodes / CHM13 kb | poa_spoa | 37 | 202 [90.7-380] | 272 | 37 | 44.9 [6.76-135] | 7 / 30 / 0 | 7.2e-08 |
| nodes / CHM13 kb | poa_abpoa_mc | 40 | 107 [73.6-178] | 184 | 40 | -2.395 [-36.8-8.858] | 25 / 15 / 0 | 0.28 |
| nodes / CHM13 kb | unit_aware | 40 | 187 [86.1-299] | 276 | 40 | 20.8 [-0.33-126] | 11 / 29 / 0 | 2.4e-05 |
| cost/opt, all pairs | mc | 40 | 1.199 [1.031-1.448] | 1.28 | - | - | - | - |
| cost/opt, all pairs | mafft_fftns2 | 40 | 1.309 [1.121-2.073] | 1.732 | 40 | 0.052 [-0.041-0.796] | 14 / 24 / 2 | 0.012 |
| cost/opt, all pairs | mafft_fftnsi | 34 | 1.243 [1.08-1.534] | 1.42 | 34 | 0 [-0.096-0.187] | 15 / 16 / 3 | 0.44 |
| cost/opt, all pairs | mafft_linsi | 27 | 1.083 [1.021-1.147] | 1.115 | 27 | -0.051 [-0.269-0] | 19 / 5 / 3 | 0.00015 |
| cost/opt, all pairs | mafft_einsi | 22 | 1.068 [1.004-1.19] | 1.14 | 22 | -0.024 [-0.126-0] | 15 / 4 / 3 | 0.0011 |
| cost/opt, all pairs | mafft_ginsi | 22 | 1.059 [1.016-1.101] | 1.081 | 22 | -0.043 [-0.269-0] | 15 / 4 / 3 | 0.00042 |
| cost/opt, all pairs | poa_abpoa | 40 | 1.114 [1.031-1.25] | 1.24 | 40 | -0.01 [-0.147-0.006] | 24 / 12 / 4 | 0.032 |
| cost/opt, all pairs | poa_spoa | 37 | 1.09 [1.034-1.159] | 1.14 | 37 | -0.06 [-0.284-0] | 26 / 8 / 3 | 0.00043 |
| cost/opt, all pairs | poa_abpoa_mc | 40 | 1.149 [1.015-1.411] | 1.274 | 40 | 0 [-0.052-0.074] | 18 / 17 / 5 | 0.85 |
| cost/opt, all pairs | unit_aware | 40 | 1.054 [1.019-1.085] | 1.105 | 40 | -0.109 [-0.306-0] | 29 / 9 / 2 | 5e-05 |
| excess edits / kb | mc | 40 | 8.875 [1.56-35.0] | 31.1 | - | - | - | - |
| excess edits / kb | mafft_fftns2 | 40 | 29.0 [7.812-73.6] | 52.1 | 40 | 0.66 [-2.78-32.7] | 14 / 24 / 2 | 0.082 |
| excess edits / kb | mafft_fftnsi | 34 | 11.5 [2.237-53.2] | 35.4 | 34 | 0 [-7.05-7.125] | 15 / 16 / 3 | 0.96 |
| excess edits / kb | mafft_linsi | 27 | 3.55 [0.425-16.1] | 12.3 | 27 | -3.48 [-13.5-0] | 19 / 5 / 3 | 6.4e-05 |
| excess edits / kb | mafft_einsi | 22 | 2.205 [0.068-32.0] | 19.1 | 22 | -1.345 [-6.01-0] | 14 / 2 / 6 | 0.0021 |
| excess edits / kb | mafft_ginsi | 22 | 2.01 [0.228-16.5] | 10.6 | 22 | -2.065 [-14.8-0] | 15 / 4 / 3 | 0.00034 |
| excess edits / kb | poa_abpoa | 40 | 7.38 [0.873-26.0] | 19.1 | 40 | -0.565 [-9.643-0.072] | 24 / 12 / 4 | 0.015 |
| excess edits / kb | poa_spoa | 37 | 5.03 [1.71-16.8] | 13.7 | 37 | -1.65 [-20.9-0] | 26 / 8 / 3 | 0.00058 |
| excess edits / kb | poa_abpoa_mc | 40 | 8.665 [0.825-39.4] | 28.1 | 40 | 0 [-5.247-1.873] | 18 / 17 / 5 | 0.83 |
| excess edits / kb | unit_aware | 40 | 3.245 [0.483-9.582] | 7.551 | 40 | -6.345 [-28.0-0] | 29 / 9 / 2 | 5e-06 |
| affine cost/opt, all | mc | 40 | 1.143 [1.012-1.299] | 1.195 | - | - | - | - |
| affine cost/opt, all | mafft_fftns2 | 40 | 1.234 [1.067-1.881] | 1.649 | 40 | 0.031 [-0.04-0.702] | 15 / 22 / 3 | 0.013 |
| affine cost/opt, all | mafft_fftnsi | 34 | 1.16 [1.051-1.392] | 1.354 | 34 | 0 [-0.072-0.17] | 16 / 16 / 2 | 0.36 |
| affine cost/opt, all | mafft_linsi | 27 | 1.057 [1.013-1.11] | 1.074 | 27 | -0.056 [-0.125--0.001] | 20 / 4 / 3 | 0.00013 |
| affine cost/opt, all | mafft_einsi | 22 | 1.012 [1.001-1.145] | 1.082 | 22 | -0.028 [-0.114--0.001] | 16 / 3 / 3 | 0.00042 |
| affine cost/opt, all | mafft_ginsi | 22 | 1.026 [1.012-1.087] | 1.054 | 22 | -0.046 [-0.15-0] | 15 / 4 / 3 | 0.0024 |
| affine cost/opt, all | poa_abpoa | 40 | 1.125 [1.043-1.255] | 1.304 | 40 | -0.001 [-0.069-0.06] | 20 / 16 / 4 | 0.92 |
| affine cost/opt, all | poa_spoa | 37 | 1.157 [1.074-1.242] | 1.196 | 37 | 0.001 [-0.084-0.059] | 15 / 19 / 3 | 0.97 |
| affine cost/opt, all | poa_abpoa_mc | 40 | 1.075 [1.007-1.292] | 1.198 | 40 | -0 [-0.03-0.032] | 20 / 15 / 5 | 0.95 |
| affine cost/opt, all | unit_aware | 40 | 1.069 [1.034-1.132] | 1.147 | 40 | -0.012 [-0.13-0.039] | 22 / 15 / 3 | 0.2 |
| unaligned homology bp / kb | mc | 40 | 0.395 [0-6.388] | 13.5 | - | - | - | - |
| unaligned homology bp / kb | mafft_fftns2 | 40 | 0.29 [0.025-2.013] | 3.026 | 40 | -0.215 [-5.805-0.045] | 23 / 13 / 4 | 0.0026 |
| unaligned homology bp / kb | mafft_fftnsi | 34 | 0.685 [0.095-2.81] | 2.348 | 34 | -0.19 [-2.725-0.142] | 19 / 12 / 3 | 0.1 |
| unaligned homology bp / kb | mafft_linsi | 27 | 0.02 [0-0.5] | 0.511 | 27 | -0.15 [-2.455-0.01] | 14 / 8 / 5 | 0.014 |
| unaligned homology bp / kb | mafft_einsi | 22 | 0.02 [0-0.883] | 1.251 | 22 | 0 [-2.155-0.033] | 9 / 7 / 6 | 0.093 |
| unaligned homology bp / kb | mafft_ginsi | 22 | 0.02 [0-0.747] | 0.547 | 22 | 0 [-2.165-0.007] | 10 / 6 / 6 | 0.058 |
| unaligned homology bp / kb | poa_abpoa | 40 | 0.1 [0-0.37] | 0.528 | 40 | -0.33 [-5.808-0] | 25 / 7 / 8 | 1.8e-05 |
| unaligned homology bp / kb | poa_spoa | 37 | 0.21 [0.04-0.8] | 0.592 | 37 | -0.26 [-5.69-0] | 23 / 6 / 8 | 0.00013 |
| unaligned homology bp / kb | poa_abpoa_mc | 40 | 0.105 [0-2.025] | 5.794 | 40 | -0.095 [-2.63-0.01] | 21 / 12 / 7 | 0.0092 |
| unaligned homology bp / kb | unit_aware | 40 | 0.015 [0-0.323] | 0.372 | 40 | -0.195 [-6.248-0] | 24 / 8 / 8 | 0.00019 |
| SV pieces / path (median) | mc | 40 | 0 [0-1] | 0.9 | - | - | - | - |
| SV pieces / path (median) | mafft_fftns2 | 40 | 1 [0-3] | 1.675 | 40 | 0 [0-1] | 2 / 15 / 23 | 0.0042 |
| SV pieces / path (median) | mafft_fftnsi | 34 | 0.5 [0-2] | 1.265 | 34 | 0 [0-0] | 3 / 8 / 23 | 0.23 |
| SV pieces / path (median) | mafft_linsi | 27 | 0 [0-1] | 0.704 | 27 | 0 [0-0] | 2 / 2 / 23 | 0.62 |
| SV pieces / path (median) | mafft_einsi | 22 | 0 [0-1] | 0.591 | 22 | 0 [0-0] | 2 / 1 / 19 | 0.75 |
| SV pieces / path (median) | mafft_ginsi | 22 | 0 [0-1] | 0.591 | 22 | 0 [0-0] | 1 / 0 / 21 | 1 |
| SV pieces / path (median) | poa_abpoa | 40 | 0 [0-1] | 0.925 | 40 | 0 [0-0] | 7 / 5 / 28 | 0.89 |
| SV pieces / path (median) | poa_spoa | 37 | 0 [0-1] | 1 | 37 | 0 [0-0] | 5 / 4 / 28 | 0.97 |
| SV pieces / path (median) | poa_abpoa_mc | 40 | 0 [0-1] | 0.925 | 40 | 0 [0-0] | 2 / 5 / 33 | 0.81 |
| SV pieces / path (median) | unit_aware | 40 | 0 [0-1] | 1.075 | 40 | 0 [0-0] | 6 / 5 / 29 | 1 |
| indel bp / net length change | mc | 40 | 1.343 [1.055-1.7] | 1.736 | - | - | - | - |
| indel bp / net length change | mafft_fftns2 | 40 | 1.523 [1.152-2.688] | 2.524 | 40 | 0.096 [-0.04-0.955] | 13 / 24 / 3 | 0.034 |
| indel bp / net length change | mafft_fftnsi | 34 | 1.368 [1.11-2.069] | 1.943 | 34 | 0 [-0.169-0.176] | 15 / 16 / 3 | 0.72 |
| indel bp / net length change | mafft_linsi | 27 | 1.058 [1.016-1.378] | 1.388 | 27 | -0.123 [-0.345-0] | 17 / 4 / 6 | 0.00029 |
| indel bp / net length change | mafft_einsi | 22 | 1.054 [1.017-1.409] | 1.33 | 22 | -0.04 [-0.219-0] | 13 / 2 / 7 | 0.00085 |
| indel bp / net length change | mafft_ginsi | 22 | 1.035 [1.012-1.077] | 1.215 | 22 | -0.189 [-0.469-0] | 14 / 3 / 5 | 0.00021 |
| indel bp / net length change | poa_abpoa | 40 | 1.115 [1.022-1.554] | 1.487 | 40 | -0.034 [-0.289-0] | 26 / 7 / 7 | 0.0011 |
| indel bp / net length change | poa_spoa | 37 | 1.1 [1.027-1.502] | 1.443 | 37 | -0.092 [-0.536-0] | 27 / 4 / 6 | 1.5e-05 |
| indel bp / net length change | poa_abpoa_mc | 40 | 1.233 [1.032-1.737] | 1.655 | 40 | -0.001 [-0.094-0.013] | 21 / 11 / 8 | 0.29 |
| indel bp / net length change | unit_aware | 40 | 1.034 [1.008-1.435] | 1.348 | 40 | -0.154 [-0.549--0.003] | 31 / 4 / 5 | 9.7e-06 |
| k-mer extra positions (frac) | mc | 40 | 0.159 [0.036-0.29] | 0.193 | - | - | - | - |
| k-mer extra positions (frac) | mafft_fftns2 | 40 | 0.248 [0.12-0.476] | 0.303 | 40 | 0.077 [0.001-0.18] | 9 / 30 / 1 | 8.2e-06 |
| k-mer extra positions (frac) | mafft_fftnsi | 34 | 0.212 [0.074-0.492] | 0.275 | 34 | 0.05 [-0.003-0.154] | 10 / 23 / 1 | 0.00095 |
| k-mer extra positions (frac) | mafft_linsi | 27 | 0.173 [0.038-0.28] | 0.19 | 27 | -0.002 [-0.042-0.013] | 17 / 9 / 1 | 0.34 |
| k-mer extra positions (frac) | mafft_einsi | 22 | 0.057 [0.017-0.27] | 0.153 | 22 | -0.008 [-0.081--0] | 16 / 5 / 1 | 0.046 |
| k-mer extra positions (frac) | mafft_ginsi | 22 | 0.098 [0.032-0.278] | 0.168 | 22 | -0.003 [-0.043-0.008] | 14 / 7 / 1 | 0.24 |
| k-mer extra positions (frac) | poa_abpoa | 40 | 0.186 [0.059-0.317] | 0.215 | 40 | 0.006 [-0.01-0.053] | 17 / 22 / 1 | 0.1 |
| k-mer extra positions (frac) | poa_spoa | 37 | 0.26 [0.071-0.392] | 0.255 | 37 | 0.041 [0.004-0.096] | 7 / 29 / 1 | 9.3e-07 |
| k-mer extra positions (frac) | poa_abpoa_mc | 40 | 0.119 [0.037-0.202] | 0.149 | 40 | -0.031 [-0.088--0.002] | 33 / 6 / 1 | 3.3e-06 |
| k-mer extra positions (frac) | unit_aware | 40 | 0.185 [0.069-0.403] | 0.243 | 40 | 0.03 [-0.006-0.102] | 16 / 23 / 1 | 0.001 |
| truth edits to best graph path (h1+h2) | mc | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_linsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_einsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_abpoa | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_spoa | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | unit_aware | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_linsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_einsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_abpoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_spoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | unit_aware | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_linsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_einsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_abpoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_spoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | unit_aware | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mc | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftns2 | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftnsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_linsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_einsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_ginsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_abpoa | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_spoa | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | unit_aware | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mc | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_fftns2 | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_fftnsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_linsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_einsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_ginsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_abpoa | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_spoa | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| reads: placements / max per path | unit_aware | 0 | - | - | - | - | - | - |

| method | regions scored (raw / refined / phab) | raw TP-base / FN / TP-comp / FP | raw F1 | refined TP-base / FN / TP-comp / FP | refined F1 | phab TP-base / FN / TP-comp / FP | phab F1 |
|---|---|---|---|---|---|---|---|
| mc | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_fftns2 | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_fftnsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_linsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_einsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_ginsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_abpoa | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_spoa | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_abpoa_mc | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| unit_aware | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |

## control_vntr_correct (25 regions)

| method | n | nodes / CHM13 kb | cost/opt, all pairs | excess edits / kb | affine cost/opt, all | unaligned homology bp / kb | SV pieces / path (median) | indel bp / net length change | k-mer extra positions (frac) | truth edits to best graph path (h1+h2) | truth-by-graph truvari F1, raw | truth-by-graph truvari F1, refined | reads: redundant (cross-walk) only | reads: placements / max per path | regions with cost/opt <= 1.1 | pooled truvari F1 raw / refined |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc | 25 | 112 | 1.056 | 6.91 | 1.044 | 0.28 | 1 | 1.171 | 0.141 | - | - | - | - | - | 15/25 | - / - |
| mafft_fftns2 | 25 | 138 | 1.118 | 8.82 | 1.094 | 0.36 | 1 | 1.205 | 0.209 | - | - | - | - | - | 11/25 | - / - |
| mafft_fftnsi | 25 | 137 | 1.103 | 8.54 | 1.088 | 0.2 | 1 | 1.174 | 0.207 | - | - | - | - | - | 12/25 | - / - |
| mafft_linsi | 20 | 101 | 1.043 | 2.4 | 1.032 | 0.025 | 1 | 1.028 | 0.11 | - | - | - | - | - | 18/20 | - / - |
| mafft_einsi | 18 | 81.9 | 1.037 | 1.995 | 1.024 | 0.03 | 1 | 1.035 | 0.105 | - | - | - | - | - | 14/18 | - / - |
| mafft_ginsi | 18 | 101 | 1.048 | 2.405 | 1.038 | 0.025 | 1 | 1.027 | 0.13 | - | - | - | - | - | 17/18 | - / - |
| poa_abpoa | 25 | 202 | 1.109 | 9.11 | 1.209 | 0.18 | 1 | 1.061 | 0.232 | - | - | - | - | - | 12/25 | - / - |
| poa_spoa | 24 | 212 | 1.116 | 6.93 | 1.262 | 0.22 | 1 | 1.079 | 0.253 | - | - | - | - | - | 11/24 | - / - |
| poa_abpoa_mc | 25 | 109 | 1.045 | 5.89 | 1.037 | 0.1 | 1 | 1.079 | 0.139 | - | - | - | - | - | 15/25 | - / - |
| unit_aware | 25 | 188 | 1.023 | 1.83 | 1.056 | 0.01 | 1 | 1.026 | 0.228 | - | - | - | - | - | 23/25 | - / - |

| metric | method | n | median [q1-q3] | mean | paired n | median diff [q1-q3] | better / worse / equal | p |
|---|---|---|---|---|---|---|---|---|
| nodes / CHM13 kb | mc | 25 | 112 [69.0-192] | 163 | - | - | - | - |
| nodes / CHM13 kb | mafft_fftns2 | 25 | 138 [76.4-293] | 231 | 25 | 11.5 [-1.42-56.8] | 8 / 16 / 1 | 0.0039 |
| nodes / CHM13 kb | mafft_fftnsi | 25 | 137 [72.9-293] | 220 | 25 | 15.9 [-2.5-28.4] | 10 / 14 / 1 | 0.0072 |
| nodes / CHM13 kb | mafft_linsi | 20 | 101 [72.6-194] | 153 | 20 | 2.93 [-1.323-25.7] | 7 / 13 / 0 | 0.11 |
| nodes / CHM13 kb | mafft_einsi | 18 | 81.9 [48.3-131] | 120 | 18 | -4.685 [-17.9--1.557] | 14 / 3 / 1 | 0.0093 |
| nodes / CHM13 kb | mafft_ginsi | 18 | 101 [82.1-189] | 161 | 18 | 7.26 [-0.218-54.3] | 5 / 13 / 0 | 0.067 |
| nodes / CHM13 kb | poa_abpoa | 25 | 202 [68.2-286] | 220 | 25 | 46.5 [-0.86-104] | 8 / 17 / 0 | 0.0042 |
| nodes / CHM13 kb | poa_spoa | 24 | 212 [98.5-310] | 212 | 24 | 46.1 [12.3-123] | 5 / 19 / 0 | 0.00011 |
| nodes / CHM13 kb | poa_abpoa_mc | 25 | 109 [62.2-213] | 159 | 25 | -3.44 [-8.08-2.55] | 16 / 9 / 0 | 0.26 |
| nodes / CHM13 kb | unit_aware | 25 | 188 [93.8-286] | 225 | 25 | 28.6 [9.86-115] | 4 / 21 / 0 | 0.0001 |
| cost/opt, all pairs | mc | 25 | 1.056 [1.022-1.335] | 1.213 | - | - | - | - |
| cost/opt, all pairs | mafft_fftns2 | 25 | 1.118 [1.026-1.316] | 1.387 | 25 | 0.002 [-0.005-0.122] | 7 / 14 / 4 | 0.27 |
| cost/opt, all pairs | mafft_fftnsi | 25 | 1.103 [1.026-1.285] | 1.323 | 25 | 0.002 [-0.032-0.079] | 8 / 13 / 4 | 0.56 |
| cost/opt, all pairs | mafft_linsi | 20 | 1.043 [1.024-1.057] | 1.049 | 20 | -0.023 [-0.146-0.001] | 11 / 7 / 2 | 0.018 |
| cost/opt, all pairs | mafft_einsi | 18 | 1.037 [1.022-1.065] | 1.055 | 18 | -0.028 [-0.136-0] | 11 / 3 / 4 | 0.0031 |
| cost/opt, all pairs | mafft_ginsi | 18 | 1.048 [1.033-1.058] | 1.047 | 18 | -0.023 [-0.144-0.005] | 10 / 6 / 2 | 0.034 |
| cost/opt, all pairs | poa_abpoa | 25 | 1.109 [1.022-1.213] | 1.222 | 25 | 0 [-0.047-0.068] | 10 / 12 / 3 | 0.77 |
| cost/opt, all pairs | poa_spoa | 24 | 1.116 [1.018-1.193] | 1.121 | 24 | -0.001 [-0.1-0.054] | 13 / 10 / 1 | 0.5 |
| cost/opt, all pairs | poa_abpoa_mc | 25 | 1.045 [1.012-1.19] | 1.128 | 25 | -0.002 [-0.08-0] | 15 / 7 / 3 | 0.079 |
| cost/opt, all pairs | unit_aware | 25 | 1.023 [1.01-1.042] | 1.037 | 25 | -0.05 [-0.165-0] | 18 / 7 / 0 | 0.00029 |
| excess edits / kb | mc | 25 | 6.91 [1.13-16.0] | 24.2 | - | - | - | - |
| excess edits / kb | mafft_fftns2 | 25 | 8.82 [2.82-46.4] | 29.8 | 25 | 0.11 [-0.79-8.3] | 7 / 14 / 4 | 0.3 |
| excess edits / kb | mafft_fftnsi | 25 | 8.54 [1.96-27.4] | 24.9 | 25 | 0.04 [-4.55-6.7] | 8 / 13 / 4 | 0.54 |
| excess edits / kb | mafft_linsi | 20 | 2.4 [0.963-6.673] | 5.501 | 20 | -0.71 [-12.4-0.12] | 11 / 7 / 2 | 0.038 |
| excess edits / kb | mafft_einsi | 18 | 1.995 [0.54-8.623] | 6.405 | 18 | -1.645 [-9.473-0] | 11 / 3 / 4 | 0.0085 |
| excess edits / kb | mafft_ginsi | 18 | 2.405 [1.14-6.08] | 4.917 | 18 | -0.71 [-10.1-0.32] | 10 / 6 / 2 | 0.074 |
| excess edits / kb | poa_abpoa | 25 | 9.11 [0.8-25.9] | 14.8 | 25 | 0 [-6.67-4.93] | 10 / 11 / 4 | 0.86 |
| excess edits / kb | poa_spoa | 24 | 6.93 [0.802-14.1] | 10.6 | 24 | -0.08 [-7.995-3.45] | 13 / 10 / 1 | 0.41 |
| excess edits / kb | poa_abpoa_mc | 25 | 5.89 [0.83-12.9] | 13.8 | 25 | -0.08 [-5.53-0.01] | 15 / 7 / 3 | 0.068 |
| excess edits / kb | unit_aware | 25 | 1.83 [0.63-3.97] | 4.33 | 25 | -4.05 [-15.4-0] | 18 / 6 / 1 | 0.00049 |
| affine cost/opt, all | mc | 25 | 1.044 [1.024-1.184] | 1.182 | - | - | - | - |
| affine cost/opt, all | mafft_fftns2 | 25 | 1.094 [1.024-1.305] | 1.366 | 25 | 0.005 [-0.023-0.105] | 9 / 14 / 2 | 0.3 |
| affine cost/opt, all | mafft_fftnsi | 25 | 1.088 [1.024-1.213] | 1.303 | 25 | 0.001 [-0.023-0.084] | 10 / 13 / 2 | 0.56 |
| affine cost/opt, all | mafft_linsi | 20 | 1.032 [1.019-1.066] | 1.038 | 20 | -0.02 [-0.105-0.002] | 13 / 6 / 1 | 0.023 |
| affine cost/opt, all | mafft_einsi | 18 | 1.024 [1.013-1.043] | 1.037 | 18 | -0.021 [-0.131-0] | 12 / 4 / 2 | 0.0027 |
| affine cost/opt, all | mafft_ginsi | 18 | 1.038 [1.027-1.07] | 1.043 | 18 | -0.017 [-0.113-0.011] | 11 / 6 / 1 | 0.098 |
| affine cost/opt, all | poa_abpoa | 25 | 1.209 [1.021-1.383] | 1.345 | 25 | 0 [-0.019-0.245] | 10 / 13 / 2 | 0.25 |
| affine cost/opt, all | poa_spoa | 24 | 1.262 [1.03-1.344] | 1.23 | 24 | 0.011 [-0.016-0.232] | 7 / 16 / 1 | 0.18 |
| affine cost/opt, all | poa_abpoa_mc | 25 | 1.037 [1.013-1.171] | 1.103 | 25 | -0.01 [-0.059-0] | 15 / 8 / 2 | 0.08 |
| affine cost/opt, all | unit_aware | 25 | 1.056 [1.02-1.088] | 1.07 | 25 | -0.02 [-0.096-0.007] | 15 / 10 / 0 | 0.11 |
| unaligned homology bp / kb | mc | 25 | 0.28 [0.03-3.17] | 10.3 | - | - | - | - |
| unaligned homology bp / kb | mafft_fftns2 | 25 | 0.36 [0.06-3.23] | 2.466 | 25 | 0 [-0.84-0.32] | 10 / 11 / 4 | 0.54 |
| unaligned homology bp / kb | mafft_fftnsi | 25 | 0.2 [0.03-2.22] | 1.388 | 25 | 0 [-1.46-0.3] | 12 / 10 / 3 | 0.24 |
| unaligned homology bp / kb | mafft_linsi | 20 | 0.025 [0-0.267] | 0.618 | 20 | -0.035 [-1.842-0.003] | 12 / 5 / 3 | 0.041 |
| unaligned homology bp / kb | mafft_einsi | 18 | 0.03 [0-0.475] | 3.092 | 18 | -0.01 [-0.843-0] | 9 / 4 / 5 | 0.24 |
| unaligned homology bp / kb | mafft_ginsi | 18 | 0.025 [0-0.128] | 0.197 | 18 | -0.035 [-2.712-0] | 11 / 4 / 3 | 0.024 |
| unaligned homology bp / kb | poa_abpoa | 25 | 0.18 [0.03-0.69] | 0.464 | 25 | -0.06 [-3.03-0.07] | 14 / 8 / 3 | 0.05 |
| unaligned homology bp / kb | poa_spoa | 24 | 0.22 [0.007-0.633] | 0.364 | 24 | -0.06 [-1.757-0.17] | 13 / 8 / 3 | 0.081 |
| unaligned homology bp / kb | poa_abpoa_mc | 25 | 0.1 [0.01-0.37] | 1.246 | 25 | -0.05 [-2.02-0] | 16 / 5 / 4 | 0.034 |
| unaligned homology bp / kb | unit_aware | 25 | 0.01 [0.01-0.13] | 0.222 | 25 | -0.27 [-3.16-0] | 18 / 5 / 2 | 0.0011 |
| SV pieces / path (median) | mc | 25 | 1 [1-1] | 0.96 | - | - | - | - |
| SV pieces / path (median) | mafft_fftns2 | 25 | 1 [1-2] | 1.4 | 25 | 0 [0-0] | 3 / 6 / 16 | 0.29 |
| SV pieces / path (median) | mafft_fftnsi | 25 | 1 [1-2] | 1.36 | 25 | 0 [0-0] | 3 / 6 / 16 | 0.33 |
| SV pieces / path (median) | mafft_linsi | 20 | 1 [0-1] | 0.85 | 20 | 0 [0-0] | 2 / 2 / 16 | 1 |
| SV pieces / path (median) | mafft_einsi | 18 | 1 [0-1] | 0.833 | 18 | 0 [0-0] | 1 / 2 / 15 | 1 |
| SV pieces / path (median) | mafft_ginsi | 18 | 1 [0-1] | 0.778 | 18 | 0 [0-0] | 1 / 1 / 16 | 1 |
| SV pieces / path (median) | poa_abpoa | 25 | 1 [1-1] | 1.04 | 25 | 0 [0-0] | 1 / 3 / 21 | 0.75 |
| SV pieces / path (median) | poa_spoa | 24 | 1 [0-1] | 1.083 | 24 | 0 [0-0] | 4 / 4 / 16 | 0.63 |
| SV pieces / path (median) | poa_abpoa_mc | 25 | 1 [1-1] | 1 | 25 | 0 [0-0] | 2 / 3 / 20 | 1 |
| SV pieces / path (median) | unit_aware | 25 | 1 [0-2] | 0.92 | 25 | 0 [0-0] | 6 / 6 / 13 | 1 |
| indel bp / net length change | mc | 25 | 1.171 [1.041-1.466] | 1.321 | - | - | - | - |
| indel bp / net length change | mafft_fftns2 | 25 | 1.205 [1.04-1.302] | 1.547 | 25 | 0 [-0.019-0.113] | 11 / 11 / 3 | 0.68 |
| indel bp / net length change | mafft_fftnsi | 25 | 1.174 [1.039-1.302] | 1.484 | 25 | 0 [-0.019-0.085] | 11 / 11 / 3 | 0.9 |
| indel bp / net length change | mafft_linsi | 20 | 1.028 [1.014-1.135] | 1.089 | 20 | -0.057 [-0.174-0] | 14 / 4 / 2 | 0.0013 |
| indel bp / net length change | mafft_einsi | 18 | 1.035 [1.017-1.153] | 1.097 | 18 | -0.057 [-0.13-0] | 12 / 3 / 3 | 0.00085 |
| indel bp / net length change | mafft_ginsi | 18 | 1.027 [1.017-1.09] | 1.078 | 18 | -0.055 [-0.143-0] | 12 / 4 / 2 | 0.0034 |
| indel bp / net length change | poa_abpoa | 25 | 1.061 [1.022-1.301] | 1.189 | 25 | -0.006 [-0.148-0.006] | 14 / 8 / 3 | 0.18 |
| indel bp / net length change | poa_spoa | 24 | 1.079 [1.014-1.191] | 1.127 | 24 | -0.018 [-0.169-0.002] | 16 / 7 / 1 | 0.0091 |
| indel bp / net length change | poa_abpoa_mc | 25 | 1.079 [1.021-1.26] | 1.171 | 25 | -0.001 [-0.133-0] | 15 / 5 / 5 | 0.024 |
| indel bp / net length change | unit_aware | 25 | 1.026 [1.006-1.124] | 1.085 | 25 | -0.078 [-0.24--0.003] | 21 / 2 / 2 | 2.1e-05 |
| k-mer extra positions (frac) | mc | 25 | 0.141 [0.064-0.272] | 0.194 | - | - | - | - |
| k-mer extra positions (frac) | mafft_fftns2 | 25 | 0.209 [0.113-0.486] | 0.28 | 25 | 0.025 [-0.005-0.127] | 7 / 16 / 2 | 0.0039 |
| k-mer extra positions (frac) | mafft_fftnsi | 25 | 0.207 [0.104-0.469] | 0.263 | 25 | 0.015 [-0.005-0.066] | 7 / 16 / 2 | 0.015 |
| k-mer extra positions (frac) | mafft_linsi | 20 | 0.11 [0.049-0.175] | 0.144 | 20 | -0.009 [-0.019-0.009] | 13 / 7 / 0 | 0.28 |
| k-mer extra positions (frac) | mafft_einsi | 18 | 0.105 [0.04-0.171] | 0.145 | 18 | -0.016 [-0.043-0.003] | 12 / 5 / 1 | 0.13 |
| k-mer extra positions (frac) | mafft_ginsi | 18 | 0.13 [0.068-0.177] | 0.155 | 18 | -0.009 [-0.017-0.018] | 12 / 6 / 0 | 0.52 |
| k-mer extra positions (frac) | poa_abpoa | 25 | 0.232 [0.062-0.388] | 0.236 | 25 | 0 [-0.01-0.138] | 12 / 11 / 2 | 0.22 |
| k-mer extra positions (frac) | poa_spoa | 24 | 0.253 [0.152-0.388] | 0.246 | 24 | 0.04 [0.009-0.122] | 6 / 18 / 0 | 0.00096 |
| k-mer extra positions (frac) | poa_abpoa_mc | 25 | 0.139 [0.042-0.208] | 0.152 | 25 | -0.022 [-0.083--0.007] | 19 / 5 / 1 | 0.0014 |
| k-mer extra positions (frac) | unit_aware | 25 | 0.228 [0.117-0.356] | 0.237 | 25 | 0.029 [-0.012-0.106] | 8 / 17 / 0 | 0.02 |
| truth edits to best graph path (h1+h2) | mc | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_linsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_einsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_abpoa | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_spoa | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | unit_aware | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_linsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_einsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_abpoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_spoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | unit_aware | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_linsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_einsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_abpoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_spoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | unit_aware | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mc | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftns2 | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftnsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_linsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_einsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_ginsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_abpoa | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_spoa | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | unit_aware | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mc | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_fftns2 | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_fftnsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_linsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_einsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_ginsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_abpoa | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_spoa | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| reads: placements / max per path | unit_aware | 0 | - | - | - | - | - | - |

| method | regions scored (raw / refined / phab) | raw TP-base / FN / TP-comp / FP | raw F1 | refined TP-base / FN / TP-comp / FP | refined F1 | phab TP-base / FN / TP-comp / FP | phab F1 |
|---|---|---|---|---|---|---|---|
| mc | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_fftns2 | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_fftnsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_linsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_einsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_ginsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_abpoa | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_spoa | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_abpoa_mc | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| unit_aware | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |

## hotspot_other (24 regions)

| method | n | nodes / CHM13 kb | cost/opt, all pairs | excess edits / kb | affine cost/opt, all | unaligned homology bp / kb | SV pieces / path (median) | indel bp / net length change | k-mer extra positions (frac) | truth edits to best graph path (h1+h2) | truth-by-graph truvari F1, raw | truth-by-graph truvari F1, refined | reads: redundant (cross-walk) only | reads: placements / max per path | regions with cost/opt <= 1.1 | pooled truvari F1 raw / refined |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc | 24 | 875 | 1.427 | 107 | 1.415 | 23.6 | 4.5 | 1.496 | 0.654 | - | - | - | - | - | 1/24 | - / - |
| mafft_fftns2 | 24 | 1457 | 1.275 | 74.8 | 1.272 | 9.25 | 2.5 | 1.246 | 0.736 | - | - | - | - | - | 2/24 | - / - |
| mafft_fftnsi | 17 | 1109 | 1.251 | 60.1 | 1.247 | 5.08 | 2 | 1.272 | 0.702 | - | - | - | - | - | 2/17 | - / - |
| mafft_linsi | 6 | 1063 | 1.134 | 33.6 | 1.128 | 4.75 | 1 | 1.068 | 0.585 | - | - | - | - | - | 2/6 | - / - |
| mafft_einsi | 4 | 804 | 1.192 | 56.4 | 1.169 | 10.9 | 2.5 | 1.144 | 0.691 | - | - | - | - | - | 0/4 | - / - |
| mafft_ginsi | 4 | 852 | 1.088 | 23.3 | 1.077 | 3.71 | 1 | 1.028 | 0.571 | - | - | - | - | - | 2/4 | - / - |
| poa_abpoa | 23 | 884 | 1.205 | 46.9 | 1.382 | 4.37 | 2 | 1.074 | 0.607 | - | - | - | - | - | 4/23 | - / - |
| poa_spoa | 22 | 1005 | 1.183 | 45.3 | 1.47 | 4.995 | 1 | 1.055 | 0.678 | - | - | - | - | - | 5/22 | - / - |
| poa_abpoa_mc | 24 | 974 | 1.324 | 77.6 | 1.32 | 9.23 | 3 | 1.163 | 0.627 | - | - | - | - | - | 4/24 | - / - |
| unit_aware | 24 | 1462 | 1.063 | 16.2 | 1.273 | 1.205 | 0 | 1.027 | 0.708 | - | - | - | - | - | 16/24 | - / - |

| metric | method | n | median [q1-q3] | mean | paired n | median diff [q1-q3] | better / worse / equal | p |
|---|---|---|---|---|---|---|---|---|
| nodes / CHM13 kb | mc | 24 | 875 [503-1820] | 1294 | - | - | - | - |
| nodes / CHM13 kb | mafft_fftns2 | 24 | 1457 [846-3063] | 2203 | 24 | 665 [207-1431] | 0 / 24 / 0 | 1.2e-07 |
| nodes / CHM13 kb | mafft_fftnsi | 17 | 1109 [723-1556] | 1436 | 17 | 343 [168-642] | 0 / 17 / 0 | 1.5e-05 |
| nodes / CHM13 kb | mafft_linsi | 6 | 1063 [692-1229] | 1123 | 6 | 215 [59.9-255] | 1 / 5 / 0 | 0.094 |
| nodes / CHM13 kb | mafft_einsi | 4 | 804 [487-1426] | 1109 | 4 | 135 [35.0-274] | 0 / 4 / 0 | 0.12 |
| nodes / CHM13 kb | mafft_ginsi | 4 | 852 [500-1474] | 1122 | 4 | 136 [12.8-310] | 0 / 4 / 0 | 0.12 |
| nodes / CHM13 kb | poa_abpoa | 23 | 884 [694-1585] | 1172 | 23 | 3.03 [-104-137] | 11 / 12 / 0 | 1 |
| nodes / CHM13 kb | poa_spoa | 22 | 1005 [789-1595] | 1247 | 22 | 126 [10.6-201] | 5 / 17 / 0 | 0.033 |
| nodes / CHM13 kb | poa_abpoa_mc | 24 | 974 [646-2018] | 1389 | 24 | 17.8 [-59.5-92.9] | 10 / 14 / 0 | 0.53 |
| nodes / CHM13 kb | unit_aware | 24 | 1462 [1146-2521] | 1910 | 24 | 616 [303-768] | 1 / 23 / 0 | 2.4e-07 |
| cost/opt, all pairs | mc | 24 | 1.427 [1.258-1.758] | 1.565 | - | - | - | - |
| cost/opt, all pairs | mafft_fftns2 | 24 | 1.275 [1.177-1.549] | 1.398 | 24 | -0.128 [-0.245--0.065] | 20 / 4 / 0 | 0.00021 |
| cost/opt, all pairs | mafft_fftnsi | 17 | 1.251 [1.175-1.571] | 1.399 | 17 | -0.21 [-0.354--0.114] | 16 / 1 / 0 | 0.0017 |
| cost/opt, all pairs | mafft_linsi | 6 | 1.134 [1.078-1.167] | 1.138 | 6 | -0.156 [-0.222--0.124] | 6 / 0 / 0 | 0.031 |
| cost/opt, all pairs | mafft_einsi | 4 | 1.192 [1.131-1.325] | 1.264 | 4 | -0.007 [-0.036-0.014] | 2 / 2 / 0 | 0.62 |
| cost/opt, all pairs | mafft_ginsi | 4 | 1.088 [1.063-1.146] | 1.122 | 4 | -0.13 [-0.173--0.114] | 4 / 0 / 0 | 0.12 |
| cost/opt, all pairs | poa_abpoa | 23 | 1.205 [1.136-1.299] | 1.245 | 23 | -0.267 [-0.488--0.177] | 22 / 1 / 0 | 4.8e-07 |
| cost/opt, all pairs | poa_spoa | 22 | 1.183 [1.109-1.28] | 1.229 | 22 | -0.286 [-0.516--0.225] | 22 / 0 / 0 | 4.8e-07 |
| cost/opt, all pairs | poa_abpoa_mc | 24 | 1.324 [1.158-1.602] | 1.41 | 24 | -0.123 [-0.255--0.041] | 20 / 4 / 0 | 0.0031 |
| cost/opt, all pairs | unit_aware | 24 | 1.063 [1.042-1.118] | 1.079 | 24 | -0.377 [-0.636--0.228] | 24 / 0 / 0 | 1.2e-07 |
| excess edits / kb | mc | 24 | 107 [86.2-164] | 128 | - | - | - | - |
| excess edits / kb | mafft_fftns2 | 24 | 74.8 [48.9-105] | 83.4 | 24 | -41.0 [-72.5--18.7] | 20 / 4 / 0 | 0.00021 |
| excess edits / kb | mafft_fftnsi | 17 | 60.1 [44.1-71.1] | 64.1 | 17 | -46.7 [-77.3--22.8] | 16 / 1 / 0 | 0.0011 |
| excess edits / kb | mafft_linsi | 6 | 33.6 [24.1-42.2] | 36.4 | 6 | -44.5 [-51.3--34.2] | 6 / 0 / 0 | 0.031 |
| excess edits / kb | mafft_einsi | 4 | 56.4 [43.0-69.8] | 56.4 | 4 | -1.45 [-10.3-2.715] | 2 / 2 / 0 | 0.88 |
| excess edits / kb | mafft_ginsi | 4 | 23.3 [17.7-31.0] | 25.4 | 4 | -37.6 [-47.3--27.4] | 4 / 0 / 0 | 0.12 |
| excess edits / kb | poa_abpoa | 23 | 46.9 [36.0-71.2] | 52.1 | 23 | -61.9 [-106--46.7] | 22 / 1 / 0 | 4.8e-07 |
| excess edits / kb | poa_spoa | 22 | 45.3 [26.4-56.5] | 46.4 | 22 | -61.6 [-102--51.6] | 22 / 0 / 0 | 4.8e-07 |
| excess edits / kb | poa_abpoa_mc | 24 | 77.6 [54.6-121] | 92.5 | 24 | -30.9 [-59.0--13.0] | 20 / 4 / 0 | 0.00065 |
| excess edits / kb | unit_aware | 24 | 16.2 [10.8-21.6] | 17.9 | 24 | -90.5 [-141--70.5] | 24 / 0 / 0 | 1.2e-07 |
| affine cost/opt, all | mc | 24 | 1.415 [1.29-1.658] | 1.485 | - | - | - | - |
| affine cost/opt, all | mafft_fftns2 | 24 | 1.272 [1.178-1.451] | 1.35 | 24 | -0.121 [-0.184--0.07] | 21 / 3 / 0 | 0.00028 |
| affine cost/opt, all | mafft_fftnsi | 17 | 1.247 [1.18-1.386] | 1.334 | 17 | -0.183 [-0.288--0.104] | 16 / 1 / 0 | 0.0021 |
| affine cost/opt, all | mafft_linsi | 6 | 1.128 [1.071-1.17] | 1.135 | 6 | -0.169 [-0.204--0.16] | 6 / 0 / 0 | 0.031 |
| affine cost/opt, all | mafft_einsi | 4 | 1.169 [1.124-1.283] | 1.238 | 4 | -0.057 [-0.072--0.049] | 4 / 0 / 0 | 0.12 |
| affine cost/opt, all | mafft_ginsi | 4 | 1.077 [1.063-1.127] | 1.113 | 4 | -0.168 [-0.211--0.146] | 4 / 0 / 0 | 0.12 |
| affine cost/opt, all | poa_abpoa | 23 | 1.382 [1.282-1.539] | 1.429 | 23 | -0.076 [-0.167-0.074] | 16 / 7 / 0 | 0.16 |
| affine cost/opt, all | poa_spoa | 22 | 1.47 [1.366-1.616] | 1.512 | 22 | -0.002 [-0.162-0.163] | 11 / 11 / 0 | 0.9 |
| affine cost/opt, all | poa_abpoa_mc | 24 | 1.32 [1.215-1.444] | 1.379 | 24 | -0.08 [-0.206--0.021] | 19 / 5 / 0 | 0.015 |
| affine cost/opt, all | unit_aware | 24 | 1.273 [1.138-1.345] | 1.28 | 24 | -0.158 [-0.285--0.081] | 21 / 3 / 0 | 1.3e-05 |
| unaligned homology bp / kb | mc | 24 | 23.6 [17.6-42.9] | 42.3 | - | - | - | - |
| unaligned homology bp / kb | mafft_fftns2 | 24 | 9.25 [6.405-11.8] | 9.415 | 24 | -15.3 [-30.1--8.282] | 22 / 2 / 0 | 1.2e-06 |
| unaligned homology bp / kb | mafft_fftnsi | 17 | 5.08 [3.87-6.56] | 5.648 | 17 | -18.6 [-24.2--6.75] | 16 / 1 / 0 | 4.6e-05 |
| unaligned homology bp / kb | mafft_linsi | 6 | 4.75 [2.958-5.35] | 4.283 | 6 | -14.6 [-38.8--6.418] | 6 / 0 / 0 | 0.031 |
| unaligned homology bp / kb | mafft_einsi | 4 | 10.9 [7.093-14.1] | 10.3 | 4 | -2.635 [-6.902-1.498] | 2 / 2 / 0 | 0.62 |
| unaligned homology bp / kb | mafft_ginsi | 4 | 3.71 [1.993-5.362] | 3.645 | 4 | -7.78 [-12.7--4.477] | 4 / 0 / 0 | 0.12 |
| unaligned homology bp / kb | poa_abpoa | 23 | 4.37 [2.51-5.165] | 4.539 | 23 | -20.3 [-36.9--12.9] | 23 / 0 / 0 | 2.4e-07 |
| unaligned homology bp / kb | poa_spoa | 22 | 4.995 [2.25-6.112] | 4.912 | 22 | -19.9 [-38.6--11.7] | 21 / 1 / 0 | 9.5e-07 |
| unaligned homology bp / kb | poa_abpoa_mc | 24 | 9.23 [4.857-20.7] | 13.3 | 24 | -12.9 [-22.6--4.932] | 21 / 2 / 1 | 1e-05 |
| unaligned homology bp / kb | unit_aware | 24 | 1.205 [0.725-2.107] | 1.698 | 24 | -22.6 [-41.6--16.0] | 23 / 1 / 0 | 2.4e-07 |
| SV pieces / path (median) | mc | 24 | 4.5 [3-6.25] | 5.5 | - | - | - | - |
| SV pieces / path (median) | mafft_fftns2 | 24 | 2.5 [1-7.25] | 3.958 | 24 | -1 [-4-0] | 15 / 4 / 5 | 0.015 |
| SV pieces / path (median) | mafft_fftnsi | 17 | 2 [1-3] | 2.824 | 17 | -1 [-3-0] | 12 / 3 / 2 | 0.055 |
| SV pieces / path (median) | mafft_linsi | 6 | 1 [1-1.75] | 2 | 6 | -0.5 [-1.75-0] | 3 / 1 / 2 | 0.62 |
| SV pieces / path (median) | mafft_einsi | 4 | 2.5 [1.75-3] | 2.25 | 4 | 0 [-0.75-0.25] | 1 / 1 / 2 | 1 |
| SV pieces / path (median) | mafft_ginsi | 4 | 1 [1-1.25] | 1.25 | 4 | -0.5 [-2-0] | 2 / 0 / 2 | 0.5 |
| SV pieces / path (median) | poa_abpoa | 23 | 2 [1-5] | 3.304 | 23 | -2 [-3.5-0] | 16 / 5 / 2 | 0.029 |
| SV pieces / path (median) | poa_spoa | 22 | 1 [0-2] | 1.682 | 22 | -3 [-4.75--1.25] | 18 / 2 / 2 | 4.6e-05 |
| SV pieces / path (median) | poa_abpoa_mc | 24 | 3 [3-5] | 3.792 | 24 | -0.5 [-2.25-0] | 12 / 4 / 8 | 0.015 |
| SV pieces / path (median) | unit_aware | 24 | 0 [0-3.75] | 3.792 | 24 | -2.5 [-4.25-0] | 15 / 5 / 4 | 0.087 |
| indel bp / net length change | mc | 24 | 1.496 [1.21-2.022] | 1.83 | - | - | - | - |
| indel bp / net length change | mafft_fftns2 | 24 | 1.246 [1.14-1.535] | 1.462 | 24 | -0.269 [-0.602--0.025] | 18 / 6 / 0 | 0.00018 |
| indel bp / net length change | mafft_fftnsi | 17 | 1.272 [1.09-1.834] | 1.536 | 17 | -0.378 [-0.814--0.043] | 14 / 3 / 0 | 0.00066 |
| indel bp / net length change | mafft_linsi | 6 | 1.068 [1.033-1.196] | 1.149 | 6 | -0.249 [-0.372--0.104] | 6 / 0 / 0 | 0.031 |
| indel bp / net length change | mafft_einsi | 4 | 1.144 [1.117-1.386] | 1.358 | 4 | -0.001 [-0.123-0.093] | 2 / 2 / 0 | 0.88 |
| indel bp / net length change | mafft_ginsi | 4 | 1.028 [1.025-1.206] | 1.203 | 4 | -0.113 [-0.212--0.085] | 4 / 0 / 0 | 0.12 |
| indel bp / net length change | poa_abpoa | 23 | 1.074 [1.039-1.326] | 1.247 | 23 | -0.43 [-0.796--0.156] | 23 / 0 / 0 | 2.4e-07 |
| indel bp / net length change | poa_spoa | 22 | 1.055 [1.039-1.306] | 1.229 | 22 | -0.448 [-0.909--0.151] | 22 / 0 / 0 | 4.8e-07 |
| indel bp / net length change | poa_abpoa_mc | 24 | 1.163 [1.108-1.592] | 1.432 | 24 | -0.203 [-0.443--0.063] | 23 / 1 / 0 | 3.6e-07 |
| indel bp / net length change | unit_aware | 24 | 1.027 [1.013-1.125] | 1.103 | 24 | -0.451 [-0.951--0.176] | 24 / 0 / 0 | 1.2e-07 |
| k-mer extra positions (frac) | mc | 24 | 0.654 [0.538-0.75] | 0.629 | - | - | - | - |
| k-mer extra positions (frac) | mafft_fftns2 | 24 | 0.736 [0.603-0.776] | 0.707 | 24 | 0.07 [0.036-0.149] | 2 / 22 / 0 | 0.00028 |
| k-mer extra positions (frac) | mafft_fftnsi | 17 | 0.702 [0.594-0.736] | 0.662 | 17 | 0.062 [0.042-0.115] | 2 / 15 / 0 | 0.00011 |
| k-mer extra positions (frac) | mafft_linsi | 6 | 0.585 [0.55-0.602] | 0.562 | 6 | -0.011 [-0.046-0.008] | 3 / 3 / 0 | 0.56 |
| k-mer extra positions (frac) | mafft_einsi | 4 | 0.691 [0.644-0.712] | 0.666 | 4 | 0.096 [0.053-0.121] | 1 / 3 / 0 | 0.25 |
| k-mer extra positions (frac) | mafft_ginsi | 4 | 0.571 [0.495-0.625] | 0.549 | 4 | -0.024 [-0.069-0.007] | 2 / 2 / 0 | 0.62 |
| k-mer extra positions (frac) | poa_abpoa | 23 | 0.607 [0.557-0.759] | 0.623 | 23 | 0.013 [-0.007-0.037] | 8 / 15 / 0 | 0.18 |
| k-mer extra positions (frac) | poa_spoa | 22 | 0.678 [0.596-0.81] | 0.659 | 22 | 0.047 [0.026-0.061] | 2 / 20 / 0 | 0.0042 |
| k-mer extra positions (frac) | poa_abpoa_mc | 24 | 0.627 [0.483-0.72] | 0.598 | 24 | -0.001 [-0.062-0.015] | 12 / 12 / 0 | 0.16 |
| k-mer extra positions (frac) | unit_aware | 24 | 0.708 [0.6-0.811] | 0.682 | 24 | 0.059 [0.008-0.116] | 5 / 19 / 0 | 0.023 |
| truth edits to best graph path (h1+h2) | mc | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_linsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_einsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_abpoa | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_spoa | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | unit_aware | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_linsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_einsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_abpoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_spoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | unit_aware | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_linsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_einsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_abpoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_spoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | unit_aware | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mc | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftns2 | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftnsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_linsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_einsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_ginsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_abpoa | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_spoa | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | unit_aware | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mc | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_fftns2 | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_fftnsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_linsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_einsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_ginsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_abpoa | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_spoa | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| reads: placements / max per path | unit_aware | 0 | - | - | - | - | - | - |

| method | regions scored (raw / refined / phab) | raw TP-base / FN / TP-comp / FP | raw F1 | refined TP-base / FN / TP-comp / FP | refined F1 | phab TP-base / FN / TP-comp / FP | phab F1 |
|---|---|---|---|---|---|---|---|
| mc | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_fftns2 | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_fftnsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_linsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_einsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_ginsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_abpoa | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_spoa | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_abpoa_mc | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| unit_aware | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |

## control_nontr_sv (20 regions)

| method | n | nodes / CHM13 kb | cost/opt, all pairs | excess edits / kb | affine cost/opt, all | unaligned homology bp / kb | SV pieces / path (median) | indel bp / net length change | k-mer extra positions (frac) | truth edits to best graph path (h1+h2) | truth-by-graph truvari F1, raw | truth-by-graph truvari F1, refined | reads: redundant (cross-walk) only | reads: placements / max per path | regions with cost/opt <= 1.1 | pooled truvari F1 raw / refined |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc | 20 | 56.4 | 1 | 0.025 | 1 | 0 | 0 | 1 | 0.006 | - | - | - | - | - | 20/20 | - / - |
| mafft_fftns2 | 20 | 48.8 | 1 | 0.035 | 1 | 0 | 0 | 1 | 0 | - | - | - | - | - | 20/20 | - / - |
| mafft_fftnsi | 20 | 48.8 | 1 | 0.025 | 1 | 0 | 0 | 1 | 0 | - | - | - | - | - | 20/20 | - / - |
| mafft_linsi | 20 | 53.2 | 1 | 0.035 | 1 | 0 | 0 | 1 | 0 | - | - | - | - | - | 19/20 | - / - |
| mafft_einsi | 20 | 49.3 | 1 | 0.025 | 1 | 0 | 0 | 1 | 0 | - | - | - | - | - | 20/20 | - / - |
| mafft_ginsi | 20 | 53.2 | 1 | 0.035 | 1 | 0 | 0 | 1 | 0 | - | - | - | - | - | 19/20 | - / - |
| poa_abpoa | 20 | 49.3 | 1 | 0.025 | 1 | 0 | 0 | 1 | 0 | - | - | - | - | - | 20/20 | - / - |
| poa_spoa | 20 | 49.3 | 1 | 0.025 | 1 | 0 | 0 | 1 | 0 | - | - | - | - | - | 20/20 | - / - |
| poa_abpoa_mc | 20 | 49.3 | 1 | 0.035 | 1 | 0 | 0 | 1 | 0 | - | - | - | - | - | 20/20 | - / - |
| unit_aware | 20 | 53.2 | 1 | 0.035 | 1 | 0 | 0 | 1 | 0 | - | - | - | - | - | 19/20 | - / - |

| metric | method | n | median [q1-q3] | mean | paired n | median diff [q1-q3] | better / worse / equal | p |
|---|---|---|---|---|---|---|---|---|
| nodes / CHM13 kb | mc | 20 | 56.4 [45.2-73.5] | 67.4 | - | - | - | - |
| nodes / CHM13 kb | mafft_fftns2 | 20 | 48.8 [40.6-63.2] | 62.5 | 20 | -3.675 [-5.588--1.942] | 18 / 0 / 2 | 7.6e-06 |
| nodes / CHM13 kb | mafft_fftnsi | 20 | 48.8 [40.6-63.2] | 62.4 | 20 | -3.795 [-5.588--1.942] | 18 / 0 / 2 | 7.6e-06 |
| nodes / CHM13 kb | mafft_linsi | 20 | 53.2 [44.2-90.0] | 88.3 | 20 | -2.56 [-5.435--0.773] | 15 / 4 / 1 | 0.11 |
| nodes / CHM13 kb | mafft_einsi | 20 | 49.3 [40.9-63.9] | 62.8 | 20 | -3.675 [-5.21--2.412] | 17 / 1 / 2 | 1.5e-05 |
| nodes / CHM13 kb | mafft_ginsi | 20 | 53.2 [44.2-90.0] | 88.2 | 20 | -2.56 [-5.435--0.773] | 15 / 4 / 1 | 0.11 |
| nodes / CHM13 kb | poa_abpoa | 20 | 49.3 [40.6-64.2] | 62.8 | 20 | -3.675 [-5.435--1.707] | 18 / 1 / 1 | 7.6e-06 |
| nodes / CHM13 kb | poa_spoa | 20 | 49.3 [40.6-65.2] | 62.7 | 20 | -3.925 [-5.435--1.942] | 18 / 2 / 0 | 3.6e-05 |
| nodes / CHM13 kb | poa_abpoa_mc | 20 | 49.3 [40.6-64.3] | 62.8 | 20 | -3.775 [-5.435--1.942] | 18 / 1 / 1 | 7.6e-06 |
| nodes / CHM13 kb | unit_aware | 20 | 53.2 [44.2-90.0] | 88.3 | 20 | -2.56 [-5.435--0.773] | 15 / 4 / 1 | 0.11 |
| cost/opt, all pairs | mc | 20 | 1 [1-1.001] | 1.002 | - | - | - | - |
| cost/opt, all pairs | mafft_fftns2 | 20 | 1 [1-1.002] | 1.003 | 20 | 0 [0-0] | 2 / 3 / 15 | 0.44 |
| cost/opt, all pairs | mafft_fftnsi | 20 | 1 [1-1.001] | 1.002 | 20 | 0 [0-0] | 2 / 2 / 16 | 0.88 |
| cost/opt, all pairs | mafft_linsi | 20 | 1 [1-1.004] | 1.018 | 20 | 0 [0-0] | 2 / 5 / 13 | 0.11 |
| cost/opt, all pairs | mafft_einsi | 20 | 1 [1-1.001] | 1.002 | 20 | 0 [0-0] | 2 / 3 / 15 | 0.44 |
| cost/opt, all pairs | mafft_ginsi | 20 | 1 [1-1.004] | 1.018 | 20 | 0 [0-0] | 2 / 5 / 13 | 0.11 |
| cost/opt, all pairs | poa_abpoa | 20 | 1 [1-1.001] | 1.002 | 20 | 0 [0-0] | 1 / 3 / 16 | 0.25 |
| cost/opt, all pairs | poa_spoa | 20 | 1 [1-1.001] | 1.002 | 20 | 0 [0-0] | 1 / 2 / 17 | 0.5 |
| cost/opt, all pairs | poa_abpoa_mc | 20 | 1 [1-1.002] | 1.003 | 20 | 0 [0-0] | 1 / 3 / 16 | 0.25 |
| cost/opt, all pairs | unit_aware | 20 | 1 [1-1.004] | 1.018 | 20 | 0 [0-0] | 2 / 5 / 13 | 0.11 |
| excess edits / kb | mc | 20 | 0.025 [0.01-0.26] | 0.221 | - | - | - | - |
| excess edits / kb | mafft_fftns2 | 20 | 0.035 [0.01-0.36] | 0.323 | 20 | 0 [0-0] | 2 / 3 / 15 | 0.31 |
| excess edits / kb | mafft_fftnsi | 20 | 0.025 [0.007-0.267] | 0.264 | 20 | 0 [0-0] | 2 / 2 / 16 | 0.62 |
| excess edits / kb | mafft_linsi | 20 | 0.035 [0.007-0.46] | 3.696 | 20 | 0 [0-0.007] | 2 / 5 / 13 | 0.078 |
| excess edits / kb | mafft_einsi | 20 | 0.025 [0.007-0.267] | 0.27 | 20 | 0 [0-0] | 2 / 3 / 15 | 0.31 |
| excess edits / kb | mafft_ginsi | 20 | 0.035 [0.007-0.46] | 3.696 | 20 | 0 [0-0.007] | 2 / 5 / 13 | 0.078 |
| excess edits / kb | poa_abpoa | 20 | 0.025 [0.01-0.267] | 0.271 | 20 | 0 [0-0] | 1 / 3 / 16 | 0.25 |
| excess edits / kb | poa_spoa | 20 | 0.025 [0.01-0.275] | 0.229 | 20 | 0 [0-0] | 1 / 2 / 17 | 0.5 |
| excess edits / kb | poa_abpoa_mc | 20 | 0.035 [0.01-0.345] | 0.313 | 20 | 0 [0-0] | 1 / 3 / 16 | 0.25 |
| excess edits / kb | unit_aware | 20 | 0.035 [0.007-0.46] | 3.696 | 20 | 0 [0-0.007] | 2 / 5 / 13 | 0.078 |
| affine cost/opt, all | mc | 20 | 1 [1-1] | 1 | - | - | - | - |
| affine cost/opt, all | mafft_fftns2 | 20 | 1 [1-1] | 1.002 | 20 | 0 [0-0] | 3 / 2 / 15 | 1 |
| affine cost/opt, all | mafft_fftnsi | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 3 / 1 / 16 | 0.38 |
| affine cost/opt, all | mafft_linsi | 20 | 1 [1-1] | 1.032 | 20 | 0 [0-0] | 3 / 4 / 13 | 0.38 |
| affine cost/opt, all | mafft_einsi | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 3 / 2 / 15 | 0.62 |
| affine cost/opt, all | mafft_ginsi | 20 | 1 [1-1] | 1.032 | 20 | 0 [0-0] | 3 / 4 / 13 | 0.38 |
| affine cost/opt, all | poa_abpoa | 20 | 1 [1-1] | 1 | 20 | 0 [0-0] | 2 / 2 / 16 | 1 |
| affine cost/opt, all | poa_spoa | 20 | 1 [1-1] | 1.001 | 20 | 0 [0-0] | 1 / 3 / 16 | 0.25 |
| affine cost/opt, all | poa_abpoa_mc | 20 | 1 [1-1] | 1.004 | 20 | 0 [0-0] | 2 / 2 / 16 | 0.88 |
| affine cost/opt, all | unit_aware | 20 | 1 [1-1] | 1.032 | 20 | 0 [0-0] | 3 / 4 / 13 | 0.38 |
| unaligned homology bp / kb | mc | 20 | 0 [0-0] | 0.216 | - | - | - | - |
| unaligned homology bp / kb | mafft_fftns2 | 20 | 0 [0-0] | 0.046 | 20 | 0 [0-0] | 1 / 2 / 17 | 1 |
| unaligned homology bp / kb | mafft_fftnsi | 20 | 0 [0-0] | 0.046 | 20 | 0 [0-0] | 1 / 2 / 17 | 1 |
| unaligned homology bp / kb | mafft_linsi | 20 | 0 [0-0.013] | 0.199 | 20 | 0 [0-0.003] | 1 / 5 / 14 | 0.44 |
| unaligned homology bp / kb | mafft_einsi | 20 | 0 [0-0.003] | 0.037 | 20 | 0 [0-0] | 1 / 4 / 15 | 0.56 |
| unaligned homology bp / kb | mafft_ginsi | 20 | 0 [0-0.013] | 0.189 | 20 | 0 [0-0.003] | 1 / 5 / 14 | 0.44 |
| unaligned homology bp / kb | poa_abpoa | 20 | 0 [0-0] | 0 | 20 | 0 [0-0] | 2 / 0 / 18 | 0.5 |
| unaligned homology bp / kb | poa_spoa | 20 | 0 [0-0] | 0.01 | 20 | 0 [0-0] | 1 / 1 / 18 | 1 |
| unaligned homology bp / kb | poa_abpoa_mc | 20 | 0 [0-0] | 0.024 | 20 | 0 [0-0] | 1 / 1 / 18 | 1 |
| unaligned homology bp / kb | unit_aware | 20 | 0 [0-0.013] | 0.199 | 20 | 0 [0-0.003] | 1 / 5 / 14 | 0.44 |
| SV pieces / path (median) | mc | 20 | 0 [0-1] | 0.4 | - | - | - | - |
| SV pieces / path (median) | mafft_fftns2 | 20 | 0 [0-1] | 0.4 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_fftnsi | 20 | 0 [0-1] | 0.4 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_linsi | 20 | 0 [0-1] | 0.4 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_einsi | 20 | 0 [0-1] | 0.4 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | mafft_ginsi | 20 | 0 [0-1] | 0.4 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | poa_abpoa | 20 | 0 [0-1] | 0.4 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | poa_spoa | 20 | 0 [0-1] | 0.4 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | poa_abpoa_mc | 20 | 0 [0-1] | 0.4 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| SV pieces / path (median) | unit_aware | 20 | 0 [0-1] | 0.4 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | mc | 20 | 1 [1-1] | 1.004 | - | - | - | - |
| indel bp / net length change | mafft_fftns2 | 20 | 1 [1-1] | 1.005 | 20 | 0 [0-0] | 0 / 1 / 19 | 1 |
| indel bp / net length change | mafft_fftnsi | 20 | 1 [1-1] | 1.004 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | mafft_linsi | 20 | 1 [1-1.006] | 1.007 | 20 | 0 [0-0] | 0 / 3 / 17 | 0.25 |
| indel bp / net length change | mafft_einsi | 20 | 1 [1-1] | 1.004 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | mafft_ginsi | 20 | 1 [1-1.006] | 1.007 | 20 | 0 [0-0] | 0 / 3 / 17 | 0.25 |
| indel bp / net length change | poa_abpoa | 20 | 1 [1-1] | 1.004 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | poa_spoa | 20 | 1 [1-1] | 1.004 | 20 | 0 [0-0] | 0 / 0 / 20 | - |
| indel bp / net length change | poa_abpoa_mc | 20 | 1 [1-1] | 1.005 | 20 | 0 [0-0] | 0 / 1 / 19 | 1 |
| indel bp / net length change | unit_aware | 20 | 1 [1-1.006] | 1.007 | 20 | 0 [0-0] | 0 / 3 / 17 | 0.25 |
| k-mer extra positions (frac) | mc | 20 | 0.006 [0.001-0.009] | 0.006 | - | - | - | - |
| k-mer extra positions (frac) | mafft_fftns2 | 20 | 0 [0-0.007] | 0.004 | 20 | -0 [-0.006-0] | 10 / 3 / 7 | 0.27 |
| k-mer extra positions (frac) | mafft_fftnsi | 20 | 0 [0-0.007] | 0.004 | 20 | -0 [-0.006-0] | 10 / 3 / 7 | 0.17 |
| k-mer extra positions (frac) | mafft_linsi | 20 | 0 [0-0] | 0.019 | 20 | -0.005 [-0.008-0] | 12 / 4 / 4 | 0.56 |
| k-mer extra positions (frac) | mafft_einsi | 20 | 0 [0-0] | 0.001 | 20 | -0.004 [-0.008-0] | 14 / 0 / 6 | 0.00012 |
| k-mer extra positions (frac) | mafft_ginsi | 20 | 0 [0-0] | 0.019 | 20 | -0.005 [-0.008-0] | 12 / 4 / 4 | 0.56 |
| k-mer extra positions (frac) | poa_abpoa | 20 | 0 [0-0] | 0 | 20 | -0.005 [-0.008-0] | 14 / 0 / 6 | 0.00012 |
| k-mer extra positions (frac) | poa_spoa | 20 | 0 [0-0] | 0 | 20 | -0.005 [-0.008-0] | 14 / 0 / 6 | 0.00012 |
| k-mer extra positions (frac) | poa_abpoa_mc | 20 | 0 [0-0] | 0.002 | 20 | -0.005 [-0.008-0] | 14 / 1 / 5 | 0.0084 |
| k-mer extra positions (frac) | unit_aware | 20 | 0 [0-0] | 0.019 | 20 | -0.005 [-0.008-0] | 12 / 4 / 4 | 0.56 |
| truth edits to best graph path (h1+h2) | mc | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_linsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_einsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_abpoa | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_spoa | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth edits to best graph path (h1+h2) | unit_aware | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_linsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_einsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_abpoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_spoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, raw | unit_aware | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftns2 | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_fftnsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_linsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_einsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | mafft_ginsi | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_abpoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_spoa | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| truth-by-graph truvari F1, refined | unit_aware | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mc | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftns2 | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_fftnsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_linsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_einsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | mafft_ginsi | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_abpoa | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_spoa | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| reads: redundant (cross-walk) only | unit_aware | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mc | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_fftns2 | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_fftnsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_linsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_einsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | mafft_ginsi | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_abpoa | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_spoa | 0 | - | - | - | - | - | - |
| reads: placements / max per path | poa_abpoa_mc | 0 | - | - | - | - | - | - |
| reads: placements / max per path | unit_aware | 0 | - | - | - | - | - | - |

| method | regions scored (raw / refined / phab) | raw TP-base / FN / TP-comp / FP | raw F1 | refined TP-base / FN / TP-comp / FP | refined F1 | phab TP-base / FN / TP-comp / FP | phab F1 |
|---|---|---|---|---|---|---|---|
| mc | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_fftns2 | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_fftnsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_linsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_einsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| mafft_ginsi | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_abpoa | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_spoa | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| poa_abpoa_mc | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |
| unit_aware | 0 / 0 / 0 | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - | 0 / 0 / 0 / 0 | - |

