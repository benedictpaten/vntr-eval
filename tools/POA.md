# Partial-order alignment realigners: `poa_abpoa` and `poa_spoa`

`tools/realign_poa.py` aligns a region's haplotype sequences with a partial-order aligner (abPOA 1.5.7 or
spoa 4.1.5), and `tools/msa_graph.py` builds the graph from the resulting MSA. It plugs into
`tools/realign.py`, so the two POA methods share everything but the aligner call with the mafft methods:
identical sequences are aligned once, N runs are handled the same way, and outputs, checks and the runtime
table are the same.

## The two methods

| method | command | input order | memory cap |
|---|---|---|---|
| `poa_abpoa` | `abpoa -m 0 -r 1`: global; default scores (match 2, mismatch 4); convex gaps `-O 4,24 -E 2,1`; adaptive band `-b 10 -f 0.01` | longest sequence first (what `abpoa -L` does; sorted in Python, ties in input order) | 12 GB |
| `poa_spoa` | `spoa -l 1 -r 1`: global; default scores (5 / -4); convex gaps `-g -8 -e -6 -q -10 -c -4`; unbanded | longest sequence first | 12 GB |

Both tools always run in global mode. Every sequence runs from anchor to anchor, and local or semi-global
alignment would leave the shared flanks unaligned. spoa's default, `-l 0`, is local, so it is never used.
Both tools are single-threaded. The pilot below chose the configurations and the input order.

## How to run

```bash
cd ~/PycharmProjects/vntr-eval
# Region mode: writes candidates/<method>/<id>.msa.fa, <id>.gfa, <id>.realign.json,
# and one row per region in results/realign_runtime.tsv
python3 tools/realign_poa.py regions/L012184 --method poa_abpoa [--out candidates]
python3 tools/realign_poa.py all --method poa_spoa [--jobs 1] [--stratum hotspot_vntr] [--force | --retry-failed]
python3 tools/realign.py poa_abpoa all            # the same, through realign.py's own CLI

# MSA-only mode, for ANY FASTA (e.g. the full HPRC panel). It writes only the MSA: one row per input
# record, named by the first word of its header, in input order, with '-' gaps, and every row spells its
# input sequence. Identical sequences are aligned once and their rows copied. The exit status is 0 when
# ok, 4 on timeout or memout, and 2 on error.
python3 tools/realign_poa.py --method poa_abpoa --input regions/L014297/hprc.fa.gz --msa-out all.msa.fa [--json info.json]
# Optionally project onto hap32 in the same call. For each hap32.fa record it keeps the row with the
# identical sequence (full-panel records are named hprc#k, so rows are matched by sequence, not name),
# renames it and drops all-gap columns:
python3 tools/realign_poa.py --method poa_abpoa --input regions/L014297/hprc.fa.gz --msa-out all.msa.fa \
    --project-to regions/L014297/hap32.fa --projected-out hap32rows.msa.fa
mkdir -p candidates/poa_abpoa__all      # msa_graph.py's CLI needs the output directory to exist (see Notes)
python3 tools/msa_graph.py hap32rows.msa.fa regions/L014297/hap32.fa candidates/poa_abpoa__all/L014297.gfa

# Pilot variants (the method name gets a suffix unless you pass --name):
python3 tools/realign_poa.py L012184 --method poa_abpoa --config p --order given --name abpoa_p --out /scratch/cand --no-runtime
python3 tools/realign_poa.py L012184 --method poa_spoa --order random:2 --out /scratch/cand --no-runtime
python3 tools/test_realign_poa.py                 # order functions and a synthetic run of both tools
```

In Python:

```python
import sys; sys.path.insert(0, 'tools')
import realign, realign_poa               # importing realign_poa registers poa_abpoa and poa_spoa
info = realign.align_fasta('poa_abpoa', 'seqs.fa', 'out.msa.fa', timeout=900, mem_mb=12000)
info['status']                            # 'ok', 'timeout', 'memout' or 'error'; out.msa.fa is written only when 'ok'
name, spec = realign_poa.variant('poa_spoa', order='guide')   # a variant as a realign method spec
realign_poa.project_to_hap32('all.msa.fa', 'regions/L014297/hap32.fa', 'hap32rows.msa.fa')
```

Options: `--config` takes an abPOA flag set: `default` (none), `p` (`-p`), `S` (`-S`), `Sp` (`-S -p`) or
`nb` (`-b -1`, unbanded). `--tool-args '...'` passes any other flags. `--order` takes `given` (input order;
CHM13 first in hap32.fa), `longest`, `shortest`, `guide` or `random:SEED`. `guide` is nearest-first on a
k-mer distance (1 minus the weighted Jaccard of 16-mer counts, so copy number counts): it starts at the
medoid, then keeps adding the sequence closest to any already added (Prim's order). The remaining
options are `--timeout` (default 900 s), `--mem-mb` (default 12000), `--no-dedup`, `--merge-blocks` (not
the default), `--workdir` (default `$VNTR_WORK/realign`) and `--no-wait-memory`. The tools can be
overridden with `VNTR_ABPOA` and `VNTR_SPOA`.

## What happens around the aligner (shared with the mafft methods)

- **Duplicates:** identical sequences are aligned once and their rows copied. POA alignment ignores
  edge weights, so this changes nothing except the run time.
- **N runs:** runs of non-ACGT characters are cut out before alignment. Each run is then put back as
  columns of its own right after the preceding base, where only that row has bases, and msa_graph.py turns
  it into one N node.
  - The hap32 panel has one N run: L012184 `recombination#29#chr4#0`, 3,146 bp. In both POA graphs it is
    one N node of 3,146 bp.
  - This differs from docs/plan.md, which says "N runs become path gaps". evaluate.py rejects a path that
    does not spell its sequence, so N runs have to stay in the graph as nodes.
- **Checks:** the aligner's output is checked row by row (every row must spell its input), and the graph
  must spell every hap32 sequence. Every region that ran passed both checks.
- **Memory:** before a run, peak memory is predicted from the longest sequence, and a region predicted
  over the cap is not run (status `memout`, with "predicted" in its message). The prediction comes from
  `/usr/bin/time -l` measurements on the pilot:
  - spoa preallocates DP matrices for 4 x max_len graph nodes x max_len positions (five int32 matrices
    under convex gaps). That is 80 bytes x max_len², measured at 80-81 on four regions of 3-12 kb. On
    hap32 the graph never outgrew the preallocation. With 12 GB, spoa fits sequences up to ~12.2 kb.
  - abPOA without `-S` allocates graph nodes x length: 18-29 B x max_len² measured at 6-20 kb (L001909,
    20 kb: 9.0 GB default, 9.4 GB with `-p`). The predictor uses 22, or 24 above 16 kb. It is only
    meant to refuse regions that cannot fit.

  A running aligner is also killed if its RSS passes the cap. RSS is sampled once a second, so a fast
  allocator can overshoot: abPOA `-S -p` at L001909 reached 16.3 GB before a 12 GB kill. In region mode,
  a region predicted to need more than 2 GB first waits until the machine has that much free
  (`--no-wait-memory` turns this off).

## Pilot

**Design.** The pilot has 11 regions:
- the 7 hotspot VNTRs: L012184, L005990 and L002013 (named in the task), L009656, L014297 and L015415
  (the pilot loci from findings.md), and L001909 (20 kb, the long-region test);
- the 4 matched controls of four of them: L011430 (matched to L012184), L007172 (to L005990), L007065 (to
  L002013) and L016124 (to L014297).

Each variant was built with `realign_poa.py` and scored with
`evaluate.py <region_dir> <gfa> --skip truth` (Stage 0 graph metrics, no truth, no reads). The criterion
is `all_cost_over_opt`: the cost of the pairwise alignment the graph implies, over the optimal pairwise
edit distance, summed over all path pairs (lower is better, 1 is optimal). The `mc` column is the
baseline from results/mc. The mafft FFT-NS-2 and L-INS-i columns score the mafft agent's candidates the
same way; L-INS-i had not finished L001909 at the time.

### abPOA: flags and input order (cost/opt)

A `*` marks a hotspot. Variant names are `<flags>/<order>`; `given` means CHM13 first.

| region | mc | FFT-NS-2 | L-INS-i | **default/longest** | default/given | default/guide | -p/given | -p/longest | -S/given | -S -p/given | -b -1/given |
|---|---|---|---|---|---|---|---|---|---|---|---|
| L012184* | 1.599 | 1.860 | 1.143 | **1.091** | 1.238 | 1.167 | 1.154 | 1.154 | 3.289 | 2.682 | 1.251 |
| L005990* | 1.211 | 1.098 | 1.063 | **1.044** | 1.152 | 1.047 | 1.048 | 1.048 | 1.451 | 1.246 | 1.101 |
| L002013* | 1.406 | 1.143 | 1.100 | **1.098** | 1.230 | 1.177 | 1.185 | 1.167 | 1.894 | 1.457 | 1.178 |
| L009656* | 1.335 | 1.075 | 1.099 | **1.062** | 1.085 | 1.073 | 1.064 | 1.064 | 1.884 | 2.428 | 1.067 |
| L014297* | 1.303 | 1.288 | 1.051 | **1.067** | 1.097 | 1.068 | 1.068 | 1.068 | 2.410 | 1.503 | 1.278 |
| L015415* | 1.184 | 1.154 | 1.049 | **1.011** | 1.021 | 1.012 | 1.012 | 1.012 | 1.375 | 1.562 | 1.048 |
| L001909* | 1.959 | 1.388 | - | **1.138** | 1.403 | 1.172 | 1.141 | 1.141 | 2.520 | memout | 1.502 |
| L011430 | 1.361 | 1.175 | 1.203 | **1.268** | 1.214 | 1.478 | 1.377 | 1.377 | 2.016 | 1.247 | 1.214 |
| L007172 | 1.103 | 1.076 | 1.031 | **1.094** | 1.073 | 1.227 | 1.070 | 1.070 | 1.227 | 1.236 | 1.050 |
| L007065 | 1.000 | 1.000 | 1.002 | **1.000** | 1.000 | 1.000 | 1.000 | 1.000 | 1.732 | 1.004 | 1.000 |
| L016124 | 1.585 | 1.713 | 1.054 | **1.045** | 1.045 | 1.045 | 1.045 | 1.045 | 1.820 | 3.724 | 1.045 |
| median, hotspots | 1.335 | 1.154 | 1.081 (6) | **1.067** | 1.152 | 1.073 | 1.068 | 1.068 | 1.894 | 1.533 (6) | 1.178 |
| mean, hotspots | | | | **1.073** | 1.175 | 1.102 | 1.096 | 1.093 | 2.118 | 1.813 (6) | 1.204 |
| median, controls | 1.232 | 1.126 | 1.043 | **1.069** | 1.059 | 1.136 | 1.058 | 1.058 | 1.776 | 1.242 | 1.047 |
| mean, all 11 | | | | **1.083** | 1.142 | 1.133 | 1.106 | 1.104 | 1.965 | 1.809 (10) | 1.158 |
| better / worse than default/longest | | | | | 2 / 7 | 0 / 7 | 1 / 6 | 1 / 6 | 0 / 11 | 1 / 9 | 2 / 7 |

**Choice: default flags, longest first.** It has the lowest hotspot mean (1.073) and the lowest mean over
all 11 regions (1.083), and it is at least as good as `-p` at all 7 hotspots. `-p`, abPOA's guide tree,
is second: its hotspot median ties (1.068), but at the matched control L011430 it places CHM13 badly
(CHM13-pair cost/opt 1.73). Longest-first lost to `-p` at one region, the control L007172 (1.094 against
1.070).

What the flags do:
- **`-S` (minimizer seeding and partition) fails in these repeats.** Its anchors chain wrongly, and the
  MSA has 92,688 columns at L012184, against 13,051 with default flags (both counts include the 3,146
  N columns), and 152,769 at L001909, against 20,692: haplotypes are left side by side instead of
  aligned. `-S -p` is no better, and at L001909 it passed the 12 GB cap and was killed (RSS 16.3 GB at
  the kill).
- **Unbanded (`-b -1`) is no better than the adaptive band** and is 3x slower, so the band is not what
  limits abPOA here.

### spoa: input order (cost/opt)

| region | mc | L-INS-i | given | **longest** | guide | abPOA default/longest |
|---|---|---|---|---|---|---|
| L012184* | 1.599 | 1.143 | 1.196 | **1.081** | 1.132 | 1.091 |
| L005990* | 1.211 | 1.063 | 1.063 | **1.027** | 1.024 | 1.044 |
| L002013* | 1.406 | 1.100 | 1.222 | **1.077** | 1.119 | 1.098 |
| L009656* | 1.335 | 1.099 | 1.069 | **1.050** | 1.049 | 1.062 |
| L014297* | 1.303 | 1.051 | 2.111 | **1.040** | 1.065 | 1.067 |
| L015415* | 1.184 | 1.049 | 1.035 | **1.016** | 1.014 | 1.011 |
| L001909* | 1.959 | - | memout | **memout** | memout | 1.138 |
| L011430 | 1.361 | 1.203 | 1.162 | **1.179** | 1.156 | 1.268 |
| L007172 | 1.103 | 1.031 | 1.112 | **1.138** | 1.096 | 1.094 |
| L007065 | 1.000 | 1.002 | 1.000 | **1.000** | 1.000 | 1.000 |
| L016124 | 1.585 | 1.054 | 1.019 | **1.019** | 1.022 | 1.045 |
| median, hotspots (6) | 1.319 | 1.081 | 1.132 | **1.045** | 1.057 | 1.064 |
| mean, hotspots (6) | | | 1.282 | **1.049** | 1.067 | |
| median, controls | 1.232 | 1.043 | 1.066 | **1.079** | 1.059 | 1.069 |
| mean, all 10 | | | 1.199 | **1.063** | 1.068 | |

**Choice: longest first.** It has the lowest hotspot median and mean. Guide order is close (5 regions
better, 4 worse) and slightly better on the controls. The given order (CHM13 first) is worst: at L014297
it scores 2.11, worse than MC's 1.30. L001909 (20 kb) is predicted to need 31 GB, so it was not run.

### Does the input order matter?

Three hotspots were aligned with each tool from three random input orders (seeds 1-3), next to the given,
longest-first and guide orders:

| region | tool | given | longest | guide | random 1 | random 2 | random 3 | spread over random | spread over all 6 |
|---|---|---|---|---|---|---|---|---|---|
| L012184 | abPOA (default flags) | 1.238 | **1.091** | 1.167 | 1.197 | 1.255 | 1.147 | 0.108 | 0.164 |
| L012184 | spoa | 1.196 | **1.081** | 1.132 | 1.484 | 1.256 | 1.176 | 0.309 | 0.404 |
| L012184 | abPOA -p | 1.154 | 1.154 | | 1.154 | 1.154 | 1.154 | 0 | 0 |
| L005990 | abPOA (default flags) | 1.152 | **1.044** | 1.047 | 1.394 | 1.218 | 1.204 | 0.191 | 0.351 |
| L005990 | spoa | 1.063 | 1.027 | **1.024** | 1.419 | 1.280 | 1.093 | 0.325 | 0.394 |
| L005990 | abPOA -p | 1.048 | 1.048 | | 1.048 | 1.048 | 1.048 | 0 | 0 |
| L002013 | abPOA (default flags) | 1.230 | **1.098** | 1.177 | 1.235 | 1.273 | 1.363 | 0.129 | 0.265 |
| L002013 | spoa | 1.222 | **1.077** | 1.119 | 1.209 | 1.277 | 1.434 | 0.225 | 0.357 |
| L002013 | abPOA -p | 1.185 | 1.167 | | 1.185 | 1.185 | 1.185 | 0 | 0.017 |

- **Order matters a great deal for plain POA.** Between random orders the cost ratio moves by 0.11-0.19
  with abPOA and 0.23-0.33 with spoa. That is as large as the gap between MC and the best alignment.
- **A bad order can be worse than MC.** spoa from random order 1 scores 1.419 at L005990, where MC
  scores 1.211.
- **Longest-first is the best of the six orders** in 5 of 6 region-tool pairs; guide order is 0.003
  better in the sixth.
- **With `-p` the order hardly matters.** The three random orders gave byte-identical MSAs. Only the
  longest-first order changed L002013 (1.167 against 1.185).

## Full runs (hap32, all 149 packaged regions)

Both chosen configurations were run on every packaged region with
`realign_poa.py all --method <m> --jobs 1`, one region at a time per method, and the two methods ran
at the same time.

| | poa_abpoa | poa_spoa |
|---|---|---|
| built (valid, acyclic graph) | **146 / 149** | **135 / 149** |
| memory (predicted over 12 GB, not run) | 3 | 14 |
| timeouts (900 s) | 0 | 0 |
| aligner time: total, median, max | 103 s, 0.11 s, 10.4 s (L001084) | 323 s, 0.92 s, 21.6 s (L005990) |
| aligner time, median at hotspot VNTRs | 1.2 s | 5.2 s |

- **Not run, abPOA (3):** L011138 (27.7 kb, predicted 17.7 GB), L015347 (33.6 kb, 25.8 GB), L012272
  (55.4 kb, 70.2 GB).
- **Not run, spoa (14):** all 13 hotspot VNTRs with a sequence over 12.2 kb, plus the matched control
  L007376 (19.4 kb). By predicted memory: L012294 12.2 GB (just over the cap), L001264 15.5, L009000 15.9,
  L004145 18.8, L007009 28.3, L007376 28.6, L007040 28.8, L013244 29.0, L008976 29.9, L001909 31.3,
  L001084 38.9, L011138 58.8, L015347 86.0 and L012272 234 GB.
- **Retries.** These were not retried: the machine is shared, and at the end of the run swap was 12.6 of
  14.3 GB used with 13.6 GB free. On an idle machine, spoa L012294 and abPOA L011138 should fit:
  `realign_poa.py L012294 --method poa_spoa --mem-mb 13000 --force` and
  `realign_poa.py L011138 --method poa_abpoa --mem-mb 19000 --force`. The rest cannot fit in 32 GB.

Every region that ran passed both checks: its MSA rows spell their sequences, and its graph spells every
hap32 path (re-checked by evaluate.py; `vg construct -M` and the native builder gave identical graphs
every time). Each attempt has a row in `results/realign_runtime.tsv`, including the memouts with their
reason in `note`. That file is the table realign.py maintains under a file lock (columns `method`,
`region_id`, `status`, `align_s`, `graph_s`, `total_s`, `peak_rss_mb` and others), so no separate
`.poa.tsv` was needed.

**Stage 0 on all regions** (`evaluate.py --skip truth`, no reads; the numbers are in the scratch dir,
not in results/). The table is limited to the 135 regions where abPOA, spoa and the mafft agent's
L-INS-i all have a graph, so it covers 27 of the 40 hotspot VNTRs. Cells are median [IQR], with the
count of regions at cost/opt <= 1.1:

| stratum | n | mc | poa_abpoa | poa_spoa | mafft_linsi |
|---|---|---|---|---|---|
| hotspot_vntr | 27 | 1.599 [1.335-2.167], 1 | 1.083 [1.046-1.189], 17 | **1.067 [1.040-1.105], 20** | 1.114 [1.079-1.183], 11 |
| control_vntr_matched | 39 | 1.103 [1.006-1.357], 18 | 1.028 [1.002-1.073], 32 | **1.016 [1.002-1.049], 34** | 1.031 [1.002-1.066], 32 |
| control_vntr_correct | 25 | 1.015, 17 | 1.006, 22 | **1.005, 24** | 1.015, 22 |
| hotspot_other | 24 | 1.625 [1.269-1.878], 3 | 1.150 [1.074-1.229], 9 | **1.125 [1.054-1.171], 10** | 1.164 [1.090-1.297], 8 |
| control_nontr_sv | 20 | 1.000, 20 | 1.000, 20 | 1.000, 20 | 1.000, 19 |

On their own coverage, abPOA takes the 37 hotspot VNTRs it built from a median of 1.567 [1.352-1.810] to
**1.083 [1.044-1.111]**: 25 of 37 at or under 1.1, where MC has 2 of 37.

The other metrics agree, with one exception: **gap-affine cost/opt favours L-INS-i.**
- Hotspot VNTRs: L-INS-i 1.115, abPOA 1.119, spoa 1.135 (MC 1.571).
- hotspot_other: L-INS-i 1.138, abPOA 1.197, spoa 1.230.
- L-INS-i also leaves less k-mer redundancy (0.42 against 0.45 for abPOA and 0.52 for spoa at hotspot
  VNTRs).

All three methods:
- cut unaligned homology from 24 to 0.9-1.6 bp/kb;
- cut nodes/kb from 682 to 293-345;
- cut SV pieces per path from 8 to 4-5 at hotspot VNTRs.

**Controls.** Regions worse than MC by more than 0.01 in cost/opt:

| method | region | stratum | mc | candidate |
|---|---|---|---|---|
| abPOA | L011224 | control_vntr_correct | 1.027 | 1.137 |
| abPOA | L001491 | control_vntr_correct | 1.013 | 1.050 |
| abPOA | L013000 | control_vntr_correct | 1.005 | 1.034 |
| spoa | L011224 | control_vntr_correct | 1.027 | 1.115 |
| spoa | L007172 | control_vntr_matched | 1.103 | 1.138 |

L011224 is a 10 bp CA-rich repeat (`GCACACACAT`, ~211 copies) with alleles of 2.6-3.1 kb, and both POA
tools misplace it. Every other control change is under 0.01, and most controls improve.

To get the full metrics (truth, truvari, reads) into results/, run
`python3 tools/evaluate_all.py poa_abpoa poa_spoa --reads --jobs 2 --threads 2`. That has not been done.

## Notes and surprises

- **Input order matters as much as the choice of aligner.** At the three hotspots, a random order spreads
  the cost ratio by up to 0.33, and one spoa order was worse than MC. Adding the longest sequence first
  was best or within 0.003 of best everywhere. A likely reason: once the longest array is in the graph,
  every later haplotype aligns to it with deletions, instead of each longer haplotype having to insert
  its extra units somewhere. So the order is part of the method.
  - **For the full panel, check that rule first.** Longest-first there starts from an outlier: at
    L012184 the longest full-panel haplotype is 39.6 kb, against a median of 6.5 kb and a hap32 maximum
    of 9.7 kb. 45 of 149 regions have a full-panel maximum over 1.5 times the hap32 maximum.
- **`-S` is unusable on VNTRs.** Minimizer anchors in a tandem array chain wrongly, and whole haplotypes
  end up side by side (92,688 columns against 13,051 at L012184).
- **Equally good alignments disagree.** Take the residue pairs each MSA aligns between CHM13 and every
  other haplotype. At four of six pilot hotspots, abPOA, spoa and L-INS-i share only 18-31% of them
  (Jaccard; L012184, L002013, L014297, L015415), though their cost ratios are within 0.07 of each other.
  At L005990 and L009656 they share 61-77%. FFT-NS-2 against L-INS-i is no closer (20-81%). This is the
  "many equal-cost optima" caveat in plan.md, measured: where the units are stacked is arbitrary, so how
  truth records line up with bubbles will differ between methods whose Stage-0 scores are the same.
- **Spot check at L005990 (chr17), from coarse per-row gap maps of the two MSAs.** abPOA aligns every
  array from its left end and puts the extra units of the longer haplotypes at the 3' end. L-INS-i
  instead leaves CHM13 gapped for about 1 kb near the start of the array, where the other rows have
  bases. That resembles MC's layout, which inserts about 38 units near the start in every non-CHM13
  haplotype. The two score nearly the same (1.044 and 1.063). CHM13 is the shortest sequence there, so
  longest-first adds it last.
- **POA wins on unit edit cost and L-INS-i on affine cost.** POA spreads a length difference over more
  indels. At L014297, CHM13 against each haplotype comes out as 5 SV-sized pieces with abPOA or spoa,
  and 1 with L-INS-i. Compensating indels are gone in all three: graph indel bp over net indel bp is
  1.00-1.20, against 1.22-4.95 in MC at the pilot hotspots.
- **spoa's memory is fixed by its preallocation**, 80 B x max_len², whatever the graph does. It measured
  the same with 299 distinct full-panel sequences at L014297: 5.1 GB predicted, 5.1 GB measured. So the
  longest sequence alone decides whether spoa can run.
  - On the full panel, abPOA fits under 12 GB (by prediction) in 135 of 149 regions and spoa in 123.
  - abPOA's graph, and so its memory, also grows with the number of haplotypes, so expect abPOA to need
    more than predicted there.
- **Full-panel smoke test (L014297 only; not a result).** The full panel ran in MSA-only mode: 457
  records, 299 distinct; abPOA took 13 s, spoa 67 s. The MSA was projected to the 34 hap32 rows by
  sequence identity and built into a graph. Cost/opt came out 1.132 (abPOA) and 1.109 (spoa), against
  1.067 and 1.040 from aligning hap32 alone.
  - hprc.fa.gz names records `hprc#k`, not by hap32 name, so `msa_graph.py --ignore-extra` alone does
    not project. It reports that all 32 recombination rows have no MSA row. Rows have to be matched to
    hap32.fa by sequence, which `--project-to` / `project_to_hap32` now does. At L014297 all 34 hap32
    sequences had an identical full-panel row, and the MSA went from 8,418 to 3,020 columns.
- **Bug in msa_graph.py (not fixed here; it is the tools agent's file).** Its CLI fails with
  FileNotFoundError when the output GFA's directory does not exist. `vg_graph` makes its temporary
  directory there before anything creates it. Region mode is not affected, because realign.py creates
  the directory first.

## Full-panel arm (`poa_abpoa__all`, `poa_spoa__all`)

[poa_panel.py](poa_panel.py) runs the two chosen configurations, unchanged (same flags, longest first,
same N masking and output checks), on the deduplicated full panel of each region
(`work/panel/union/<id>.fa`: every distinct sequence of `hprc.fa.gz` and `hap32.fa`; HG002 and its
parents are not in the eval graph, gref_CHM13 is excluded). Per region and method it writes:

| output | what |
|---|---|
| `work/panel/<method>/<id>.msa.fa.gz` | the full MSA, one row per distinct sequence (`u0001..`, see the union's `.map.tsv`) |
| `work/panel/<method>/<id>.gfa` | its graph (`panel.py panel-graph`), one P line per `hprc.fa.gz` record; every path re-read and checked against `hprc.fa.gz` |
| `candidates/<method>__all/<id>.msa.fa` | the MSA restricted to the hap32 rows (`panel.py project`), all-gap columns dropped |
| `candidates/<method>__all/<id>.gfa` | `msa_graph.py` graph of that; every path re-read and checked against `hap32.fa` |
| `candidates/<method>__all/<id>.realign.json` | how it was made (also for regions that were not run) |
| `results/realign_runtime.all.poa.tsv` | one row per (method, region): `region_id, method, panel, n_distinct, seconds, status`, then sizes, predicted and sampled memory, per-step times, graph sizes, notes |

```bash
python3 tools/poa_panel.py run --methods poa_abpoa,poa_spoa --jobs 2 --timeout 1800 --mem-mb 12000 --budget-mb 14000
python3 tools/poa_panel.py one poa_spoa L014297                 # one job in this process
python3 tools/poa_panel.py run --methods poa_abpoa --order random:1 --regions L014297 \
    --cand-root SCRATCH/cand --panel-root SCRATCH/panel --no-runtime     # a variant, kept out of candidates/
python3 tools/poa_panel.py predict                              # predicted memory per region
python3 tools/test_poa_panel.py
python3 tools/evaluate_all.py mc poa_abpoa poa_spoa --panel     # the full-panel graphs vs full MC
```

Regions run smallest first (total distinct bp), at most two jobs at once, each a separate process (the
graph building is Python). A job starts only when the predicted memory of all running jobs stays within
`--budget-mb` and the OS reports enough available. The full-panel graph is built with msa_graph's native
engine: `vg construct -M` takes 45-75 s on 300-400 rows against 1-2 s. At L014297 and L005990 the two
engines gave the same graph up to node numbering (every path cut into the same nodes).

**Memory model.** spoa: 80 B x len1² (its preallocation for 4 x len1 graph nodes; the full-panel
graphs stay near 1.2 x len1 nodes, and no job went above its prediction). abPOA: 11 B x len1 x len2
(int16 cells) or 22 B when len1 > 13 kb (int32), where len1 >= len2 are the two longest sequences. With
longest first, the second sequence is aligned to a graph of at least len1 nodes. The hap32 arm's
22-24 B x len1² rule would refuse every region whose longest haplotype is a lone outlier. Measured over
predicted for abPOA: median 0.78, maximum 1.40 (L004145: 7.4 GB against 5.3 GB predicted). Inside a job
this replaces `realign_poa.predict_mb`. A region predicted over the 12 GB cap is not run.

### Result of the run (2026-09-26)

| | poa_abpoa | poa_spoa |
|---|---|---|
| done (graphs valid, every path checked) | 138 / 149 | 123 / 149 |
| not run: predicted over 12 GB | 11 | 26 |
| timeouts (1800 s), errors | 0, 0 | 0, 0 |
| aligner time: total, median, max | 2,010 s, 1.0 s, 269 s (L009658) | 3,724 s, 9.0 s, 269 s (L014160) |
| aligner time over the hap32 arm, median per region | 10x | 10.5x |
| peak sampled RSS | 8.9 GB (L002013) | 10.5 GB (L005990) |

Not run, abPOA: L000289, L001084, L001909, L007009, L007040, L008976, L011138, L012184, L012272, L015347
(hotspot VNTRs) and L016870 (hotspot_other). spoa: those 11, plus L001264, L002013, L004145, L007081,
L009000, L009448, L009658, L010328, L012294, L013244 (hotspot VNTRs), L007798 (hotspot_other), L007376,
TR368009, TR424743 (matched controls) and L001922 (correct control). In every case one or two outlier
haplotypes set the memory: at L012184 the two longest are 39.6 and 34.3 kb, against 4.8 kb for CHM13. At
L000289 (12.7 GB predicted), L011138 (14.5) and L016870 (14.4), abPOA might fit on an idle machine with
a 20 GB cap. Given the 0.44-1.40 spread of measured over predicted, that is not certain.

`tools/poa_panel.py run --methods poa_abpoa --regions L000289,L011138,L016870 --mem-mb 20000 --budget-mb 20000 --retry-failed --jobs 1`

### Aligning the full panel and projecting is worse than aligning hap32 alone

Stage 0 (`evaluate.py --skip truth`, no reads), paired over the regions where MC, m(hap32) and m__all
all exist. Medians; the MC column is over the abPOA set (the spoa set differs only by its missing
regions). "worse" counts regions where m__all's cost/opt is more than 0.005 above the other:

| stratum (n abPOA / spoa) | MC | abPOA hap32 | abPOA all | spoa hap32 | spoa all | all worse than hap32 (abPOA, spoa) | all worse than MC |
|---|---|---|---|---|---|---|---|
| hotspot VNTR (30 / 20), cost/opt | 1.58 | 1.08 | 1.13 | 1.05 | 1.11 | 22/30, 19/20 | 0/30, 1/20 |
| hotspot VNTR, nodes/kb | 647 | 274 | 436 | 251 | 386 | | |
| hotspot VNTR, regions at cost/opt <= 1.1 | 1 | 19 | 11 | 14 | 8 | | |
| hotspot other (23 / 22), cost/opt | 1.65 | 1.16 | 1.24 | 1.13 | 1.24 | 20/23, 22/22 | 2/23, 2/22 |
| matched control (40 / 37), cost/opt | 1.10 | 1.02 | 1.08 | 1.02 | 1.07 | 22/40, 25/37 | 12/40, 10/37 |
| matched control, regions at cost/opt <= 1.1 | 19 | 33 | 20 | 32 | 26 | | |
| correct control (25 / 24), cost/opt | 1.01 | 1.01 | 1.08 | 1.00 | 1.07 | 16/25, 16/24 | 9/25, 8/24 |
| non-TR SV (20 / 20) | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0, 0 | 0-1 |

- Paired median change from hap32 to all, hotspot VNTRs:
  - cost/opt: +0.052 (abPOA) and +0.063 (spoa);
  - affine cost/opt: +0.17 and +0.15;
  - nodes/kb: +171 and +151;
  - k-mer extra: +0.13 and +0.12.
  - SV-sized pieces per haplotype hold or fall.
- **Against MC, the full-panel arm still wins at hotspots:** 30/30 (abPOA) and 19/20 (spoa) regions are
  better. At the controls it is about MC's level, and worse than MC in about a quarter to a third of
  regions.
- **The intersection is biased:** it leaves out the 10 (abPOA) and 20 (spoa) hotspot VNTRs with the
  longest outlier alleles.
- **What changes, seen by eye at L014297 (spoa):**
  - CHM13 carries 153 copies of a 17 bp unit (3,012 bp). The hap32 rows are mostly 1.4-2.3 kb; the
    full-panel union is longer (median 3.5 kb, maximum 8.1 kb).
  - Aligned alone, each hap32 row's deletion against CHM13 is a few long blocks.
  - Aligned inside the full panel, starting from the 8.1 kb allele, the same deletion is cut into
    unit-sized pieces: 17 bp deleted, one unit aligned, 17 bp deleted, and so on.
  - Per distinct haplotype, the indel runs against CHM13 go from 28.5 (mean 37 bp) to 45.4 (mean 24 bp).
    The same holds at L005990, L009656 and L015415 for both tools (1.2-2.8x as many runs, shorter).
  - Only 25-27% of the CHM13-to-haplotype residue pairs are the same in the two MSAs.
  - Each extra run costs a gap open in every pair that crosses it, hence the higher affine and unit
    cost. The pieces are below SV size, hence the unchanged SV-piece count.
- **Full-panel graphs judged on full-panel pairs** (`evaluate.py --panel`, the same 300 seeded pairs as
  MC's), cost/opt for MC / abPOA / spoa:
  - hotspots: L005990 1.48 / 1.11 / 1.11; L009656 4.34 / 1.33 / 1.43; L014297 1.83 / 1.08 / 1.07;
    L015415 1.44 / 1.02 / 1.02; L002013 1.44 / 1.12 / -;
  - matched controls: L011430 1.26 / 1.28 / 1.26; L007172 1.20 / 1.47 / 1.45; L007065 1.00 / 1.01 / 1.02;
    L016124 1.40 / 1.22 / 1.07.

### Input order at full-panel size

Three random orders against longest first, at three hotspots that both tools could run. The hap32 rows
reproduce the hap32 table above at L005990. All-pairs cost/opt:

| region | tool | scored on | longest | random 1 | random 2 | random 3 | spread over random |
|---|---|---|---|---|---|---|---|
| L014297 | abPOA | hap32 alone | 1.067 | 1.278 | 1.065 | 1.065 | 0.21 |
| | | all, projected | **1.133** | 1.185 | 1.750 | 1.421 | 0.57 |
| | | all, full graph | **1.080** | 1.422 | 1.357 | 1.459 | 0.10 |
| L014297 | spoa | hap32 alone | 1.040 | 1.120 | 1.039 | 1.041 | 0.08 |
| | | all, projected | **1.108** | 1.800 | 2.583 | 1.639 | 0.94 |
| | | all, full graph | **1.069** | 1.515 | 1.543 | 1.545 | 0.03 |
| L009656 | abPOA | hap32 alone | 1.062 | 1.059 | 1.077 | 1.060 | 0.02 |
| | | all, projected | 1.127 | **1.094** | 1.152 | 1.100 | 0.06 |
| | | all, full graph | 1.330 | **1.292** | 1.376 | 1.450 | 0.16 |
| L009656 | spoa | hap32 alone | 1.050 | 1.054 | 1.114 | 1.092 | 0.06 |
| | | all, projected | 1.114 | 1.122 | 1.196 | **1.073** | 0.12 |
| | | all, full graph | 1.432 | **1.310** | 1.687 | 1.523 | 0.38 |
| L005990 | abPOA | hap32 alone | 1.044 | 1.394 | 1.218 | 1.204 | 0.19 |
| | | all, projected | **1.184** | 1.657 | 1.441 | 1.539 | 0.22 |
| | | all, full graph | **1.112** | 1.368 | 1.288 | 1.325 | 0.08 |
| L005990 | spoa | hap32 alone | 1.027 | 1.419 | 1.280 | 1.093 | 0.33 |
| | | all, projected | **1.220** | 1.642 | 1.617 | 1.346 | 0.30 |
| | | all, full graph | **1.105** | 1.397 | 1.361 | 1.345 | 0.05 |

- **Order matters at least as much with 250-370 distinct sequences as with 34.** On the projected graph
  the spread over random orders is 0.06-0.94, against 0.02-0.33 on hap32. One spoa order scores 2.58 at
  L014297, where MC scores 1.30.
- **Longest first is still the best order in 4 of 6 region-tool pairs.** L009656 is the exception for
  both tools: there it loses to one or two random orders, by 0.03-0.04 projected and 0.04-0.12 on the
  full graph. So the rule holds on the full panel, even though it starts from outlier alleles.
- These numbers are Stage 0 only, at three regions. The Evaluate phase has the full comparison.

## `poa_abpoa_mc`: abPOA run the way Minigraph-Cactus runs it

One more aligner: abPOA with every setting Minigraph-Cactus uses in its BAR phase (pangenome mode), in
place of abPOA's CLI defaults. It is [realign_poa_mc.py](realign_poa_mc.py) (hap32 arm, tests in
`test_realign_poa_mc.py`) and [poa_panel_mc.py](poa_panel_mc.py) (full-panel arm, tests in
`test_poa_panel_mc.py`), not an entry in realign.py's method table: BAR's length rule needs each region's
anchor lengths, which that plugin interface does not pass.

### MC's settings against abPOA's defaults

Cactus checkout `~/CLionProjects/cactus`, commit c07e5b4. `config.xml` is
`src/cactus/cactus_progressive_config.xml`; `poaBarAligner.c` is `bar/impl/poaBarAligner.c`.

| setting | `poa_abpoa` (abPOA CLI defaults) | `poa_abpoa_mc` (Cactus) | source |
|---|---|---|---|
| binary | abPOA 1.5.7 | abPOA 1.5.4, Cactus's vendored copy | `submodules/abPOA/src/abpoa.c:19` |
| mode | global (`-m 0`) | global | `poaBarAligner.c:33` |
| substitution scores | match 2, mismatch 4 | 5x5 ACGTN matrix: match 91 (A, T) or 100 (C, G), transitions -61, transversions -114 to -125, N -100 (N/N 100) | `config.xml:308`, `poaBarAligner.c:60-77` |
| convex gaps | `-O 4,24 -E 2,1` | `-O 400,1200 -E 30,1` | `config.xml:309-312`, `poaBarAligner.c:40-43` |
| adaptive band | `-b 10 -f 0.01` | `-b 1000 -f 0.1` | `config.xml:306-307`, `poaBarAligner.c:36-37` |
| minimizer seeding | off | off | `config.xml:315`, `poaBarAligner.c:46` |
| k, w, min_w | 19, 10, 500 (unused) | 15, 5, 500, and used: `-p` sketches minimizers with them even with seeding off, so `-k 15 -w 5 -n 500` are passed | `config.xml:316-318`, `poaBarAligner.c:48-50`, `abPOA src/abpoa_seed.c:691-709` |
| order | longest first (sorted in Python) | longest first, then `-p` (guide tree); `-p` dropped above 5,000 rows | `config.xml:319-321`, `poaBarAligner.c:1054-1080, 566-571` |
| identical sequences | aligned once | every haplotype is a row | (BAR aligns every cap) |
| long sequences | aligned whole | BAR's 10 kb rule (below) | `src/cactus/setup/cactus_align.py:336-341`, `poaBarAligner.c:1117-1130` |

The command Cactus itself prints for a BAR call (`dump_abpoa_input`, `poaBarAligner.c:171-227`) is
`abpoa in.fa -O 400,1200 -E 30,1 -b 1000 -f 0.1 -t matrix.mtx -r 1 -m 0 -p`, with the matrix written
tab-separated under the header `\tA\tC\tG\tT\tN`. `realign_poa_mc.py` runs exactly that, plus
`-k 15 -w 5 -n 500`.

### The length rule, emulated

`cactus-pangenome` sets BAR's `bandingLimit` to 10 kb. Each region is treated as one BAR flower. Its two
ends are the anchor nodes, and each haplotype's adjacency string is the sequence strictly between them.

- **Every string under 10 kb:** one abPOA call on the whole strings.
- **Otherwise, each end is aligned separately** (`make_consistent_partial_order_alignments`,
  `poaBarAligner.c:751`):
  - At the left end, the first min(len, 10 kb) bases of each string.
  - At the right end, the same prefix of each reverse-complemented string.
  - A string shorter than 20 kb is in both MSAs. `trim()` cuts it once, at the point with the largest
    summed column score. `trim`, `trim_msa_suffix` and `make_column_scores` (`poaBarAligner.c:323-431`)
    are ported line for line.
  - The middle of a string longer than 20 kb is in neither MSA. It stays unaligned: one private node
    per haplotype.
- **The anchor bases** are added back as shared columns at both ends.

**Caveat:** this is BAR alone on the whole region. In the real pipeline, the CAF phase first pinches
blocks from the minigraph mappings inside the region. So BAR sees smaller flowers, and the 10 kb rule
rarely leaves anything unaligned.

Where the rule fired:
- **hap32:** 21 of 149 regions; a middle stayed unaligned in 5.
- **Full panel:** 37 of 149 regions. One or two outlier haplotypes are enough to trigger it.
  - A middle stayed unaligned in 16 of the 36 regions that were built. From L007798 (44 bp, one row) to
    L015347 (699 kb over 119 of 373 rows).
  - L012272 was not run: all 459 of its haplotypes are over 20 kb. The 15.1 Mb of unaligned middles
    would make a full MSA of ~7e9 cells, which the MSA-based pipeline cannot hold. It is recorded as
    `skipped`.

### Full-panel run

`poa_panel_mc.py` follows `poa_panel.py`: union, align, full MSA and full-panel graph in
`work/panel/poa_abpoa_mc/`, projection onto hap32 in `candidates/poa_abpoa_mc__all/`. It differs in two
ways:
- **What abPOA is given:** every `hprc.fa.gz` record, duplicates included, plus each hap32 sequence no
  full-panel record carries, once each (about 450 rows per region).
- **Which row the projection uses:** the first row of each distinct sequence. abPOA aligned the copies of
  a sequence differently for only 5 sequences in 4 regions.

**Outcome:**
- **148 of 149 regions built.** No timeouts (30 min region cap) and no memouts.
- **This includes 10 of the 11 regions `poa_abpoa__all` could not run:** L000289, L001084, L001909,
  L007009, L007040, L008976, L011138, L012184, L015347, L016870. Only L012272 is missing.

**Cost:**
- **abPOA time:** 8,039 s in total, about 18x `poa_abpoa__all` per region (median).
  - Single-call regions: median 15.5 s.
  - 10 kb-rule regions: median 166 s, maximum 401 s (L012184).
- **Batch wall clock:** 75 min at 2 jobs.
- **Peak RSS:** at most 7.1 GB (L016870). 21 regions went above 4 GB.

**Memory model.** abPOA allocates 5 int32 DP matrices of graph nodes x query length, rounded up to a power
of two. With Cactus's wide band most of it is resident: per call, the peak was a median 0.82x that
allocation (range 0-1.6x).
- **Graph size:** over calls on strings of at least 1.5 kb, abPOA graphs came out at 1.0-4.1 nodes per bp
  of the longest string (median 1.2).
- **What the predictor assumes:** 2.2 nodes per bp. A 10 kb call allocates 4 GB below 21.4k nodes and
  8 GB up to 42.9k, so this puts every 10 kb call on the 8 GB step (predicted 9.1 GB). The largest graph
  seen, 41.2k nodes (L016870), was still on that step.
- **Accuracy:** measured peaks were 0.40x the prediction (median) and at most 1.23x.
  - Two single-call regions went over their prediction: L003828 by 23% and L013531 by 11%. Both were
    far below the 10 GB cap.
  - Peak RSS itself varies between runs: L014297 peaked at 2.3 GB in one run and 3.1 GB in another,
    with byte-identical output.
- **Memory checks:** `memory_pressure` was checked before each of the 49 jobs predicted above 4 GB, and
  no job had to wait.

**Determinism:** six regions were run twice and gave byte-identical full MSAs, full graphs and projected
MSAs.

### Comparison

Each table is paired over the regions where all three graphs exist; n is given per stratum. Cells:
- **cost/opt:** median `all_cost_over_opt`, with (k) = regions at or below 1.1.
- **nodes/kb:** median.
- **b/w:** regions where `poa_abpoa_mc` is better / worse by more than 0.005 in cost/opt.

The mc columns come from `results/mc` and `results/panel/mc`. The poa_abpoa columns come from
`work/mcpoa/eval_ref/` (evaluated with the same options, because `results/poa_abpoa*` did not exist
yet).

**hap32 panel** (m(hap32); 146 regions):

| stratum (n) | cost/opt mc | poa_abpoa | poa_abpoa_mc | nodes/kb mc | poa_abpoa | poa_abpoa_mc | b/w vs mc | b/w vs poa_abpoa |
|---|---|---|---|---|---|---|---|---|
| hotspot VNTR (37) | 1.567 (2) | 1.083 (25) | 1.171 (10) | 725 | 300 | 301 | 37/0 | 1/32 |
| hotspot other (24) | 1.625 (3) | 1.150 (9) | 1.284 (6) | 492 | 270 | 214 | 22/1 | 3/19 |
| matched control (40) | 1.102 (19) | 1.024 (33) | 1.031 (31) | 101 | 32 | 32 | 20/2 | 5/16 |
| correct control (25) | 1.015 (17) | 1.006 (22) | 1.007 (22) | 86 | 19 | 15 | 9/2 | 4/6 |
| non-TR SV (20) | 1.000 (20) | 1.000 (20) | 1.000 (20) | 49 | 18 | 18 | 0/1 | 0/0 |

**Full panel, projected onto hap32** (m__all; 138 regions):

| stratum (n) | cost/opt mc | poa_abpoa | poa_abpoa_mc | nodes/kb mc | poa_abpoa | poa_abpoa_mc | b/w vs mc | b/w vs poa_abpoa |
|---|---|---|---|---|---|---|---|---|
| hotspot VNTR (30) | 1.582 (1) | 1.130 (11) | 1.359 (4) | 647 | 436 | 311 | 21/9 | 0/30 |
| hotspot other (23) | 1.653 (2) | 1.244 (3) | 1.370 (2) | 481 | 474 | 329 | 18/4 | 5/17 |
| matched control (40) | 1.102 (19) | 1.084 (20) | 1.084 (20) | 101 | 54 | 29 | 13/11 | 10/14 |
| correct control (25) | 1.015 (17) | 1.077 (15) | 1.024 (21) | 86 | 68 | 22 | 6/7 | 11/6 |
| non-TR SV (20) | 1.000 (20) | 1.000 (20) | 1.000 (20) | 49 | 18 | 18 | 0/1 | 0/0 |

**Full-panel graphs on full-panel pairs** (`evaluate.py --panel`; 138 regions):

| stratum (n) | cost/opt mc | poa_abpoa | poa_abpoa_mc | nodes/kb mc | poa_abpoa | poa_abpoa_mc | b/w vs mc | b/w vs poa_abpoa |
|---|---|---|---|---|---|---|---|---|
| hotspot VNTR (30) | 1.624 (1) | 1.111 (13) | 1.350 (3) | 1157 | 1067 | 982 | 27/3 | 1/29 |
| hotspot other (23) | 1.449 (0) | 1.205 (4) | 1.327 (3) | 779 | 884 | 959 | 18/4 | 1/21 |
| matched control (40) | 1.199 (15) | 1.114 (17) | 1.149 (16) | 138 | 152 | 107 | 15/13 | 11/19 |
| correct control (25) | 1.056 (15) | 1.109 (12) | 1.045 (15) | 112 | 202 | 109 | 11/4 | 11/7 |
| non-TR SV (20) | 1.000 (20) | 1.000 (20) | 1.000 (20) | 56 | 49 | 49 | 0/2 | 0/1 |

- **At hotspots, `poa_abpoa_mc` sits between MC and `poa_abpoa` in every arm, and the full panel widens
  the gap to `poa_abpoa`.** Paired median change from m(hap32) to m__all at hotspot VNTRs:
  - `poa_abpoa_mc`: +0.094 cost/opt (32 of 39 regions worse);
  - `poa_abpoa`: +0.052 (22 of 30 worse).
  - Projected, `poa_abpoa_mc__all` is still better than MC at 21 of 30 hotspot VNTRs, but worse at 9.
    On the full-panel graphs it is better at 27 of 30.
- **The 10 kb rule is not the main cause.** On the projected full panel, hotspot VNTRs aligned in one
  call score 1.321 (n=15) and those split at the ends 1.376 (n=15), against 1.127 and 1.156 for
  `poa_abpoa`.
- **Controls fare better under MC's settings on the full panel.** Correct controls: 1.024 against 1.077
  for `poa_abpoa__all`, at a third of its nodes/kb. In the projected table `poa_abpoa_mc__all` has the
  fewest nodes/kb in every stratum except non-TR SV, where it ties.
- **Unaligned homology rises instead.** At the 30 hotspot VNTRs of the projected table, U/kb goes from 4.6
  (hap32) to 11.1 (full panel), against 0.8 to 2.0 for `poa_abpoa` (MC: 26.0). That is homologous
  sequence left on parallel nodes.
- **The 10 regions only `poa_abpoa_mc` could build on the full panel** go from a median of 1.507 (MC) to
  1.134, projected. There are two exceptions, both with sequence left unaligned by the rule:
  - L012184: 1.599 to 1.640, with 140 kb unaligned;
  - L016870: 1.030 to 1.070, with 25.5 kb unaligned.

### Fragmentation check (full panel)

For each distinct non-CHM13 haplotype in the projected MSA: runs of insertion or deletion against CHM13
(insertion and deletion runs counted separately), then the mean over haplotypes. Mean run length in bp is
in brackets. The script is `work/mcpoa/panel/scripts/indel_runs.py split`.

| region | spoa hap32 -> all | poa_abpoa hap32 -> all | poa_abpoa_mc hap32 -> all |
|---|---|---|---|
| L014297 | 28.5 (37) -> 45.4 (24) | 9.9 (108) -> 21.4 (50) | 7.6 (140) -> 11.1 (129) |
| L005990 | 108.5 (16) -> 189.6 (10) | 75.5 (24) -> 157.3 (12) | 50.7 (51) -> 47.8 (77) |
| L009656 | 43.3 (18) -> 121.8 (7.5) | 31.0 (26) -> 76.0 (11) | 24.0 (33) -> 29.9 (36) |
| L015415 | 15.3 (99) -> 18.2 (83) | 4.9 (307) -> 14.5 (104) | 3.5 (497) -> 8.0 (192) |

- The 28.5 -> 45.4 quoted above at L014297 is spoa's. Default abPOA goes 9.9 -> 21.4 there.
- Under MC's settings the full panel does not cut deletions into unit-sized pieces:
  - runs change by 0.9-2.3x, where default abPOA's change by 2.1-3.0x;
  - runs stay longer than default abPOA's: 36-192 bp on the full panel, against 11-104 bp.
- So `poa_abpoa_mc`'s larger full-panel loss is not the fragmentation described above. The U/kb rise is
  where it shows.

```bash
python3 tools/poa_panel_mc.py plan                                   # rows, rule, predicted MB per region
python3 tools/poa_panel_mc.py run --jobs 2 --timeout 1800 --mem-mb 10000 --budget-mb 20000 \
    [--cand-root DIR]                                                # this run built in work/mcpoa/candidates, then moved
python3 tools/poa_panel_mc.py one L014297                            # one region in this process
python3 tools/test_poa_panel_mc.py                                   # ~4 s
python3 tools/evaluate_all.py poa_abpoa_mc__all --reads --jobs 2 --threads 2
python3 tools/evaluate_all.py poa_abpoa_mc --panel --jobs 2 --threads 2
python3 work/mcpoa/panel/scripts/summary_panel.py                    # the three tables
```

The runtime table is `results/mcpoa_runtime.all.tsv`: one row per region, with the rule, unaligned bp,
predicted and measured memory, per-step times, graph sizes, and the reason for any region that was not
run.
