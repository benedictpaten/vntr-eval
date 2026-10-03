# vntr-eval

Making the pangenome graph align VNTRs better.

In the HPRC v2.1 Minigraph-Cactus graph, about 1% of VNTR regions hold half of short-read `vg call`'s
SV false positives. There the graph misaligns the panel's repeat arrays: the pairwise alignments it
implies cost 1.56x the optimum, against 1.09x at matched VNTRs where genotyping works
([docs/findings.md](docs/findings.md)).

This repository tests whether re-aligning each such region's haplotypes with a suitable multiple
aligner, and rebuilding the region's graph from the alignment, fixes that. It holds:
- the test regions as self-contained sequence packages;
- the realigners and the evaluation;
- every candidate graph and its scores.

## Status

All stages are done on chr20, and the held-out test on chr6 is done.
- Stages 3-4 patch every tandem repeat whose haplotype lengths differ into the whole-contig graph, then
  re-call. They are in [results/stage4_chr20.md](results/stage4_chr20.md) and
  [results/stage4_chr6.md](results/stage4_chr6.md).
- The chosen protocol aligns the full panel to its medoid with a star alignment. Section 4x has it.
- Calling each repeat as one site in `vg call` is parked; see [docs/repeat-sites.md](docs/repeat-sites.md).
- The full chr20 graph with every tandem repeat realigned is built and tested end to end; see
  [results/full_graph_chr20.md](results/full_graph_chr20.md). The plan for its remaining SV errors is
  [docs/sv-fixes-plan.md](docs/sv-fixes-plan.md).

The Stage 0-2 results are in [docs/results.md](docs/results.md). In short:

- **The graph gets much cleaner.** On the 34 sampled haplotypes, hotspot VNTRs go from cost/opt
  1.566 to 1.048 with a repeat-unit-aware aligner, and 33/40 reach the 1.1 target (MC: 2/40).
- **The fair comparison keeps most of that.** Aligning the full ~450-haplotype panel and projecting
  onto the sampled haplotypes, which is how the MC baseline was made, gives 1.059.
- **The controls mostly hold.**
- **Read placement improves only modestly,** and only for some methods: mafft L-INS-i best. It
  follows k-mer redundancy, not cost/opt.
- **The failures:** very long regions, compound arrays and non-repeat SVs.

## Layout

```
docs/
  findings.md       the investigation that motivated this work
  plan.md           the plan: methods, the two panels, how a candidate is judged
  results.md        what was found (Stages 0-2), what it means, what to do next
regions/            149 region packages, one directory per region; README.md explains every file
  regions.tsv       the index (stratum, coordinates, anchors, repeat annotation, vg/PanGenie errors)
  <id>/             hap32.fa, hprc.fa.gz, mc.gfa, truth.fa (evaluation only), calls.tsv, region.json
  vntr_ge1kb.bed    the 3,111 benchmarked VNTRs >= 1 kb, for scaling up
candidates/<method>/<id>.{gfa,msa.fa,realign.json}         a method aligned on the 34 hap32 sequences
candidates/<method>__all/<id>.{gfa,msa.fa,realign.json}    the same method on the full panel, projected onto hap32
results/
  <method>/<id>.json        Stage 0-1 metrics of each graph (tools/evaluate.py); mc/ is the baseline
  full/<method>/<id>.json   the full-panel graphs scored on full-panel pairs
  stage01.md, stage01_tables.md, summary.md, ...   Stage 0-1 comparison
  stage2_remap.md, stage2_remap.tsv                 Stage 2 local re-map (19 regions)
  panel_stats.tsv, panel_premise.tsv                the two panels, and the check that mc.gfa is the full graph restricted to hap32
  realign_runtime*.tsv, mcpoa_runtime*.tsv          runtimes, and every region a method could not finish
  pages/<id>.html                                   before/after viewer pages for six loci
tools/              Python (standard library) and a little C; README.md lists every tool
work/               caches, full-panel MSAs and graphs, intermediates (git-ignored; set VNTR_WORK to move it)
```

## Two panels

- **hap32** (`regions/<id>/hap32.fa`): the 34 haplotypes of the graph `vg call` used. These are
  CHM13, GRCh38 and 32 haplotypes sampled for HG002 from the HPRC v2.1 graph.
- **The full panel** (`regions/<id>/hprc.fa.gz`): the up to 456 HPRC haplotypes that
  Minigraph-Cactus aligned, plus CHM13 and GRCh38.
- **The baseline `mc.gfa`** is MC's full-panel alignment restricted to the hap32 paths. So each
  method is run twice: on hap32 (`candidates/<m>/`) and on the full panel, projected
  (`candidates/<m>__all/`).
- **Both panels come from the HG002-free evaluation graph.** HG002, HG003 and HG004 are absent,
  and the `gref_CHM13` reference-cover path is excluded. Do not add their sequences from anywhere
  else. See [regions/README.md](regions/README.md).

## Look at a region

```bash
python3 tools/viewer.py L014297                                  # -> work/viewer/L014297.html
python3 tools/viewer.py L014297 --candidate candidates/unit_aware/L014297.gfa \
                                --candidate candidates/unit_aware__all/L014297.gfa --out /tmp/L014297.html
python3 tools/viewer.py regions/L014297 --no-data                # from the package alone
```

The page is one self-contained HTML file. It has:
- an MSA that can switch to each candidate's own alignment, coloured by differences or by graph
  node;
- a metrics table, a Bandage image and a node track per graph;
- calls with truvari status;
- a read pileup (only with the big data).

It needs mafft, and Bandage for the images (`--no-bandage` otherwise). The six pages in
`results/pages/` are ready to open.

## Evaluate a candidate

A candidate is a GFA with one path per `hap32.fa` record, same names, spelling the same sequences
(`regions/README.md`, "What we would like back"). From the repository root:

```bash
make -C tools/c                                                              # optional; evaluate.py builds it when needed
python3 tools/evaluate.py regions/L014297 regions/L014297/mc.gfa             # baseline -> results/mc/L014297.json
python3 tools/evaluate.py regions/L014297 candidates/<method>/L014297.gfa    # -> results/<method>/L014297.json
python3 tools/evaluate_all.py <method> --jobs 2 --threads 2                  # every region of a method
python3 tools/summarise.py --methods mc,<method> --md results/summary.md     # per-stratum table, paired against mc
```

- **Requirements:** Python 3.9 or later (standard library only), `make` and a C compiler; no vg or
  mafft.
- **What needs the big data:** the truvari part of the truth scoring, and Stage 1 (`--reads`).
  They need the evaluation repository (`VNTR_EVAL_DIR`) and report `needs_big_data` without it.
  Everything else runs from `regions/` and `tools/` alone.
- **Where the metrics are defined:** [tools/METRICS.md](tools/METRICS.md).

To build a candidate from an alignment:

```bash
python3 tools/msa_graph.py MY.msa.fa regions/<id>/hap32.fa candidates/<method>/<id>.gfa
```

For the full-panel arm, see `tools/panel.py union` and `project` in
[regions/README.md](regions/README.md). To re-run the realigners used here, see
`tools/realign.py --list` and [tools/README.md](tools/README.md).

## Data and configuration

`tools/config.py` holds every data and tool location, with environment overrides
(`python3 tools/config.py` prints them). The defaults point at the read-only evaluation repository
`~/PycharmProjects/vg-call-eval`, which holds the graphs, reads, truth and calls, and at this
repository for everything written.

These tools need that big data:
- `region.py` and `package_regions.py` (which build the packages);
- `remap_local.py` (Stage 2);
- `panel.py fetch`, `mc-graph` and `gref-check`;
- `evaluate.py --reads` and its truvari scoring.

`region.py` also needs Python 3.10 or later. Caches go to `work/` (git-ignored).
