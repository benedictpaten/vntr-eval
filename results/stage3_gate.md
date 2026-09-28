# Stage 3 gate: does local `vg call` reproduce the production calls?

> **Revised after verification (scorer version 2).** Both haplotype builders -- this gate's and
> `tools/score_haplotypes.py` -- applied vg records in POS order, so where vg's prefix trimming put a
> nested child's POS left of its parent's, the child won and the parent's allele was dropped. Both
> now apply a record before anything nested in it. The gate was re-run (`results/stage3_gate.tsv`,
> `work/stage3/gate/gate.json`): **no identity verdict changed** (the "=" / "phase only" columns and
> every count below are the same), but the ED values at L011138, L012272 and L001909 did, and the
> ED table and the numbers quoted from it below are the re-run's. The pilot report
> ([stage3_pilot.md](stage3_pilot.md)) has the consequences.

The gate calls the local MC graph for the 19 pilot regions and compares the result with the
genome-wide production calls (`work/wgs-mm095`) over each region's anchor-to-anchor span. It
compares two things:

- the two called haplotype sequences, built from each VCF by the Stage 3 scoring rule;
- truvari FP/FN.

Tool: `tools/call_local.py`. Per-region numbers: `results/stage3_gate.tsv`, and
`work/stage3/gate/gate.json`, which also holds the flags and the record-level decomposition.

## Verdict

**The calling itself passes.** The local graph has to be the genome-wide graph around the span, not
the Stage 2 graph. Given production's own alignments and a 200 kb genome-wide flank, the local call
reproduces production's haplotypes in 17 of 19 regions:

- The other two, L014297 and L012184, have every genotype, DP, DR and GQ identical. They differ
  only in the orientation of the phase, from one record onward.
- Truvari FP/FN are identical in 19 of 19 regions.

**Reads re-mapped to the local graph do not reproduce production.** Stage 3 has to re-map the reads
to every graph. With re-mapped reads, the MC call is identical to production in 8 of 19 regions.
L014297 is a ninth that differs in phase only, as above.

| | identical to production | differs |
|---|---|---|
| controls | 7 of 9 | L004021 (phase only; ED 2 against production's 6), L000506 (ED 1 against 0) |
| non-VNTR hotspot | – | L016870 (ED 0 against production's 2) |
| VNTR hotspots | 1 of 9 | 7 in 7 to 111 records; L014297 in phase only |

- **The cause of every re-map difference is the re-map.** The same graph with production's
  alignments is identical to production, except for the two phase-only cases.
- **The re-map moves 5-33% of the span reads at hotspot VNTRs,** in placement or in MAPQ.
- **A hotspot call is decided by near-ties,** so this is enough to change it.

**`gate_passed` is false under the stated criterion** (identical in nearly all regions). What this
means for Stage 3 is in the last section.

## The caller and its flags

The caller is pinned at `work/bin/vg-2a6a228a5`: a copy of `~/CLionProjects/vg/bin/vg`
(snarl-tree-order, 2a6a228a5), re-signed with `codesign -s - -f`.

**Where the flags come from.** The production command is in
`~/PycharmProjects/vg-call-eval/scripts/wgs/call_wgs.sh` (function `call_one_bed`). It was launched
by `scripts/wgs/schedule_wgs.py --work work/wgs-mm095 --extra "--mismap-max 0.95"`, whose
`--threads` defaults to 5. There was no `-r` snarls file, because none exists for wgs-mm095.

```
vg call <C>.gbz -p CHM13#0#<C> -s HG002 -d <ploidy> -t 5 --progress --mismap-max 0.95 \
    --read-likelihood --phased --mosaic-out <C>.mosaic.tsv [--ploidy-bed chrX.par.bed on chrX] \
    --gaf-base work/reads.hap32.gaf.db --gbz-base work/graph.hap32.gbz.db
```

**The flags and the binary were checked, not assumed.** Re-running this command on chr21 with the
pinned binary (`-t 4`) reproduced `work/wgs-mm095/chr21/chr21.vcf` byte for byte: 113,609 records,
header included.

**The local call keeps these flags.** It replaces only the graph and the read source
(`--gaf-reads <local GAF>`). Ploidy follows `call_wgs.sh`: `-d 1` on chrX outside the PARs, which
applies to L016870.

- Every call enumerates from 32 panel haplotypes and links over 34, as production does.
- Every REF allele matches CHM13 after the coordinate shift, in all 304 calls: 94 Stage 3 calls and
  210 diagnostic calls. `shift_vcf` refuses a VCF with any mismatch.
- The local calls are deterministic. Reversing the read order or using `-t 4` instead of `-t 1`
  changed no genotype in 6 regions tested.

## Why the Stage 2 local graph fails, and the fix

The Stage 2 local graph is `mc.gfa` plus 1 kb of CHM13 flank, called with the Stage 2 GAF. It is
identical to production in only 7 of 19 regions. Two causes are systematic.

### 1. The depth term's read rate

**How vg measures the rate.** It is in `allele_likelihood.cpp`, `local_read_stats`, and is the same
for `--gaf-reads` and GAF-Base (`RATE_WINDOW = 4096`).

- The rate is the sum over reads of (1 − e_r), counted by the node-ID window of each read's first
  node, divided by the bp of every graph node in that window.
- The windows are fixed blocks of 4,096 node IDs.
- The same window statistics also set the mean read length used by the mixture weights.

**Why the Stage 2 graph gets it wrong.** `mc.gfa` keeps the genome-wide node IDs, so its windows
are production's windows. But they contain only the span, and none of the flank nodes and reads
that production's windows hold. The flank nodes get IDs max+1 and max+2, which puts the left flank
in the last window.

**Measured.** Compare DR (observed/expected reads) with the same reads (arm `stage2orig`) against
production. The median ratio is:

- 3.02 at L009656 and 2.11 at L015415;
- 1.11-1.73 at the controls;
- 0.17 at L004021.

**The term moves hotspot calls.** Switching it off changes the called haplotypes in 5 of 19
regions on the Stage 2 graph. On the final setup it changes them in 7 of 19, all hotspots; for
example, L009656 goes from ED 2,166 to 4,333 and L015415 from 690 to 161.

### 2. Linkage context

Production's panel linkage uses sites outside the span, and the Stage 2 graph has none.

**The case that shows it.** At control L007687, production calls a het insertion and a het SNV
that the reads alone call hom and ref. Both have negative GQN: the panel overrode the reads. Both
are right (production ED 0). The Stage 2 graph calls what the reads say (ED 171).

**A genome-wide flank fixes it.** With ±10 kb of genome-wide flank (arm `gwsub10k`), L007687 is
identical to production.

### The fix: a hybrid graph

The local graph is the genome-wide hap32 graph over CHM13 span ± 200 kb, with the span replaced by
the graph under test (`build_hybrid`).

**The graph:**

- The window is taken with `gbz-base --interval --context 5000`, plus everything between the
  anchors. The context alone missed 185 and 1,748 span nodes at L011138 and L012272.
- The panel paths are the GBZ's own paths cut to the window. They are streamed once per contig
  with `vg paths -A`, because gbz-base anonymises haplotypes.
  - Fragments of one haplotype stay one panel member.
  - CHM13 and GRCh38 are REFERENCE paths under the `reference_samples` tag "CHM13 GRCh38", as in
    the genome-wide GBZ.
  - These pieces are gbz-base's walks up to orientation. The one exception is a second CHM13
    visit, which is left out, in 3 windows.

**Node IDs.** Flank nodes keep their genome-wide IDs, so does MC's span, and a candidate's span
nodes are spread over the same ID range. So the depth-rate windows mix flank and span as they do
genome-wide.

**The reads** are those the genome-wide mapping placed in the window, fetched by
`remap_local.fetch` with a 200 kb flank. They are mapped with Stage 2's exact giraffe protocol.

**Checks.** The window arms confirm each step, using production's alignments:

- The genome-wide graph over ±200 kb (`gwsub200k`) is identical to production in 17 of 19 regions.
- The hybrid MC graph (`hybrid200korig`) is identical to `gwsub200k` in 19 of 19.
- DR equals production's at every shared record in 18 of 19. At L004021 the median ratio is 0.92,
  with the calls identical anyway.

**Why 200 kb and not 50 kb.** The window size matters at hotspots:

- `gwsub50k` is identical in 16 of 19. It leaves a one-record phase flip at L009656 and 20 records
  at L001909 (phase and GQ-0 ties).
- `gwsub200k` fixes both of those, but its decode flips L012184's phase.
- Mapping cost is small: all 94 hybrid maps and calls took 14 min at 6 threads.

## The gate, per region

Legend for the table:

- "=" means the called haplotype pair is identical to production's.
- "dN" is the edit distance between the two pairs.
- "phase only" means every record has the same genotype, only oriented differently.
- "k rec" counts the records whose genotype differs or that are present in one VCF only, keyed by
  vg's snarl ID.
- FP/FN come from truvari over the span, with the genome-wide pipeline's parameters. For production,
  the span-restricted run equals the genome-wide labels in all 19 regions.

| region | stratum | ED_ref | production ED | Stage 2 graph, re-mapped | gw graph ±200 kb, gw alignments | hybrid ±200 kb, gw alignments | **hybrid ±200 kb, re-mapped (Stage 3 call)** | hybrid ±50 kb, re-mapped | FP/FN production | FP/FN Stage 3 call |
|---|---|---|---|---|---|---|---|---|---|---|
| L009656 | hotspot VNTR | 673 | 2241 | d1692, 26 rec | = | = | d87, 7 rec | d981, 8 rec | 12/2 | 12/2 |
| L014297 | hotspot VNTR | 1754 | 257 | d578, 27 rec | d242, phase only | d242, phase only | d242, phase only | d236, 27 rec | 14/1 | 14/1 |
| L015415 | hotspot VNTR | 2621 | 690 | d478, 2 rec | = | = | = | = | 12/8 | 12/8 |
| L005990 | hotspot VNTR | 846 | 6677 | d917, 17 rec | = | = | d895, 15 rec | d895, 15 rec | 23/1 | 23/1 |
| L012184 | hotspot VNTR | 8317 | 1160 | d1173, 48 rec | d2225, phase only | d2225, phase only | d434, 25 rec | d920, 23 rec | 21/8 | 22/7 |
| L002013 | hotspot VNTR | 2458 | 5063 | d3547, 74 rec | = | = | d98, 25 rec | d109, 25 rec | 20/11 | 20/11 |
| L011138 | hotspot VNTR | 7502 | 5258 | d6802, 107 rec | = | = | d1310, 34 rec | d7362, 113 rec | 17/18 | 19/19 |
| L012272 | hotspot VNTR | 7022 | 5693 | d241, 25 rec | = | = | d257, 26 rec | d242, 26 rec | 12/12 | 11/12 |
| L001909 | hotspot VNTR | 7652 | 5348 | d8687, 184 rec | = | = | d1632, 111 rec | d5684, 99 rec | 26/25 | 27/25 |
| L016870 | hotspot other | 10162 | 2 | = | = | = | d2, 1 rec | d2, 1 rec | 13/2 | 13/2 |
| L016124 | control VNTR | 299 | 0 | = | = | = | = | = | 0/0 | 0/0 |
| L007687 | control VNTR | 175 | 0 | d171, 2 rec | = | = | = | = | 0/0 | 0/0 |
| L007172 | control VNTR | 338 | 0 | = | = | = | = | = | 0/0 | 0/0 |
| L011430 | control VNTR | 276 | 0 | = | = | = | = | = | 0/0 | 0/0 |
| L007065 | control VNTR | 73 | 0 | = | = | = | = | = | 0/0 | 0/0 |
| L004021 | control VNTR | 222 | 6 | d6, phase only | = | = | d8, phase only | d8, phase only | 0/0 | 0/0 |
| L014111 | control VNTR | 303 | 0 | = | = | = | = | = | 0/0 | 0/0 |
| L000506 | control non-repeat SV | 328 | 0 | d1, 1 rec | = | = | d1, 1 rec | d1, 1 rec | 0/0 | 0/0 |
| L006556 | control non-repeat SV | 307 | 0 | = | = | = | = | = | 0/0 | 0/0 |

ED to the truth under each setup (this gate's builder, which uses the raw overlap rule; the
scorer's rule differs by up to a few tens of edits at some hotspots, e.g. the Stage 3 call at
L012184 is 1,030 under the scorer's rule):

| region | production | hybrid ±200 kb, gw alignments | **Stage 3 call** (hybrid ±200 kb, re-mapped) | hybrid ±50 kb, re-mapped | Stage 3 call without the depth term | Stage 2 graph |
|---|---|---|---|---|---|---|
| L009656 | 2241 | 2241 | 2166 | 2521 | 4333 | 3845 |
| L014297 | 257 | 257 | 257 | 251 | 253 | 833 |
| L015415 | 690 | 690 | 690 | 690 | 161 | 212 |
| L005990 | 6677 | 6677 | 6636 | 6636 | 6636 | 6666 |
| L012184 | 1160 | 1343 | 1002 | 563 | 1002 | 447 |
| L002013 | 5063 | 5063 | 5015 | 5016 | 5036 | 1704 |
| L011138 | 1977 | 1977 | 2561 | 1617 | 14788 | 14625 |
| L012272 | 1407 | 1407 | 1311 | 1353 | 1355 | 1353 |
| L001909 | 6372 | 6372 | 7471 | 3401 | 6166 | 1451 |
| L016870 | 2 | 2 | 0 | 0 | 0 | 2 |
| L004021 | 6 | 6 | 2 | 2 | 2 | 0 |
| L000506 | 0 | 0 | 1 | 1 | 1 | 1 |
| the other 7 controls | 0 | 0 | 0 | 0 | 0 | 0 (L007687: 171) |

## Every difference of the Stage 3 call, explained

**Phase only, from the decode.** At L014297 and L012184 the local chain is shorter than
production's whole-contig chain. It orients the phase differently from one record on: every
genotype, DP, DR and GQ is identical. This persists with production's alignments at ±50 and
±200 kb. At L012184 it appears only at ±200 kb.

**The re-map, at the other 10.** L009656, L005990, L012184, L002013, L011138, L012272, L001909,
L016870, L004021 and L000506 all reproduce production when the same graph is given production's
alignments (`hybrid200korig` is "=" in each). So the difference is the local re-map.

- **Placement and MAPQ.** The local re-map keeps the placement and MAPQ of 67-95% of span reads at
  hotspot VNTRs (L001909 0.67, L011138 0.68, L012184 0.70, L009656 0.74, L015415 0.82,
  L005990 0.86, L012272 0.87, L002013 0.87, L014297 0.95), against 81-98% at the controls and
  82-98% in the flanks (`work/stage3/gate/remap_diff_200k.txt`).
- **Why a local map cannot match.** It has no genome-wide competitors, and giraffe's minimizer
  filtering depends on the whole graph. The same reads mapped to the ±50 kb and ±200 kb hybrids
  of one region differ in placement for 649 of 2,083 span reads at L011138, and in MAPQ alone for
  82 more.
- **The direction on the truth is mixed:**
  - closer to the truth at L009656 (−75), L005990 (−41), L012184 (−341 against the same graph with
    production's alignments), L002013 (−48), L012272 (−96), L016870 (−2) and L004021 (−4);
  - further at L001909 (+1,099), L011138 (+584) and L000506 (+1).

**The Stage 2 graph's differences:**

- L007687 is linkage context (above).
- L004021 is phase only.
- L000506 is the re-map: `stage2orig` is identical to production.
- The hotspots add the depth-rate error to the re-map.

## What this means for Stage 3

1. **Compare the graphs against the local MC call, never against production.**
   - The Stage 3 MC call (`work/stage3/calls/mc/`) goes through the same fetch, re-map, flank and
     flags as every candidate.
   - Production is scored separately, as the task asks.
2. **The decision rule needs a noise floor at hotspots.** Two neutral changes to the MC call alone
   move hotspot ED by hundreds to thousands of edits:
   - re-mapping instead of production's alignments: |ΔED| up to 1,099 (L001909), 584 (L011138);
   - a 50 kb instead of a 200 kb flank: |ΔED| up to 4,070 (L001909), 944 (L011138),
     439 (L012184).
   - The MC controls do not move under either.

   "Any decrease counts" will therefore count noise at hotspots. **Correction:** this report went
   on to call the control check sound because MC's controls are stable. That does not follow: the
   candidates' controls move under the same neutral changes (a 50 kb flank moves mafft_linsi and
   unit_aware by 310 edits at L006556, where every MC call is exact). Each graph's own noise is
   now measured, without the truth, in [stage3_noise.md](stage3_noise.md).
3. **Tied records depend on file order.** vg can write a nested parent and its child at one POS
   with equal REF length. A builder that applies records in (POS, longest REF, file order), as the
   scoring rule does, resolves the tie by file order.
   - `bcftools sort` reorders such ties by ALT. That alone changed L012184's haplotypes by 25 edits.
   - `call_local.py` therefore keeps vg's record order, as production's VCFs do.
   - Since scorer version 2 a nested parent is applied first whatever the file order: the
     bcftools-sorted MC VCFs of L012184, L011138, L001909, L002013 and L009656 now score exactly
     as vg's order does (version 1: 1,005 against 1,030 at L012184).
4. **The scorer and this gate agree.** This gate's haplotype builder (`span_haplotypes`) is
   independent of `tools/score_haplotypes.py`.
   - Under the scorer's `raw` overlap rule they agree on ED in 209 of 209 VCFs (every arm of the
     re-run gate and production). Before version 2 they agreed too, because they shared the same
     POS-order flaw; they now share only the nesting order (`score_haplotypes.nesting_order`).
   - The scorer's default `trimmed` rule differs from `raw` at 5 hotspots of the Stage 3 MC call,
     by 2-28 edits.
   - Validation (`validate-scorer`), in all 19 regions: the truth VCF gives ED 0; the truth made
     unphased gives ED 0 again (the phase search finds it); dropping the largest record gives
     ED > 0; swapping one het's phase gives ED > 0 where there is a het; an empty VCF gives ED_ref.
5. **The truth graph is a good ceiling in 15 of 19 regions.** Its span is mafft(CHM13, HG002#1,
   HG002#2), and HG002's haplotypes carry CHM13's flank walks. It returns HG002 exactly (ED 0) in
   15 of 19 regions. The exceptions:
   - L012184 (ED 763), L011138 (88) and L012272 (53), long hotspot VNTRs;
   - L006556 (310, above its ED_ref of 307). There the 3-way mafft alignment scatters HG002's
     305 bp deletion into single-base columns: the misalignment Stage 2 reported for its own truth
     graph.

   Without CHM13's flank walks the truth haplotypes existed only inside the span, and L014297 came
   out at ED 810, not 0.
6. **Calls exist for all five graphs,** as the Stage 3 setup, in `work/stage3/calls/<graph>/<id>.vcf.gz`:
   - mc, unit_aware__all, unit_aware and truth: 19 regions each;
   - mafft_linsi: 18 (it has no L012272 graph).

   They are not scored here, and they carry the same re-map noise.

## Reproduce

```
export VNTR_STAGE2_WORK=<remap_local.py work dir>      # Stage 2 graphs, GAFs, fraglen.json
python3 tools/call_local.py pin
python3 tools/call_local.py run pilot --graphs mc,unit_aware__all,mafft_linsi,unit_aware,truth \
    --arms stage2,stage2nodepth,stage2orig,gwsub10k,gwsub50k,gwsub200k,hybrid200korig,hybrid50k,hybrid200knodepth \
    --jobs 3 --threads 2
python3 tools/call_local.py gate pilot --tsv results/stage3_gate.tsv --md work/stage3/gate/gate_table.md
python3 tools/call_local.py validate-scorer pilot
```

- **Run order.** The panel pieces are streamed once per contig; the first `run` over all 19
  regions takes about 20 min for that step. The final setup then takes 14 min for 94 maps and calls.
- **Disk.** `work/stage3` is 5.5 GB.
- **Re-running production** (on chr21) is the command in `call_local.py`'s `PRODUCTION_*` notes,
  with the pinned binary.
