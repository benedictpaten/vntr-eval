# Stage 3 pilot: stopped at the gate (revised after verification)

This version replaces the first pilot report. Verification found two blocking problems:

1. The scorer dropped parent alleles.
2. The noise remedy the report proposed measured MC's noise only.

Both are addressed in the tools below, and every affected call and score was re-run. Every ED
here comes from scorer version 2. The earlier version's numbers are quoted only as "before".

## Verdict

**The graphs were still not compared.** The gate still fails: re-run with the fixed builder, no
identity verdict of [stage3_gate.md](stage3_gate.md) changed.

**The three candidates are still unscored against the truth.**
- unit_aware__all, mafft_linsi and unit_aware have no truth ED anywhere, for their Stage 3 calls
  or their replicates.
- For them I computed only distances between two calls of the same graph, which use no truth.
- So the decision rule can still be settled before anyone sees their scores.

**What the two fixes changed:**

- **Issue 1 (scorer).** The builder now applies a vg record before anything nested in it.
  - Summed hotspot VNTR ED of the local MC call falls from 31,577 to 27,140.
  - The null spread is now −7.4% (relabel) and −14.9% (unchop) of summed ED, not −11.7% and
    −10.1%.
  - Under the fixed rule read mechanically, mc_unchop, a graph with MC's own alignment, now
    passes as a win.
  - Production is now closer to the truth than the local MC call, not further.
- **Issue 2 (noise).** Every graph now has replicate calls of its own (three per graph, two for
  mc and one for mc_relabel), and each graph's noise is measured without the truth.
  - The candidates' controls move under a neutral 50 kb flank where every MC call is unchanged
    and exact.
  - So the fixed "no control stratum worse" check is not sound for the candidates.
  - The tools now report counts that account for noise, next to the rule's own counts. They
    cannot make the fixed rule sound without changing it; that choice is yours (last section).

## Issue 1: the haplotype builder dropped parent alleles

**What was wrong.** `score_haplotypes.Prepared` applied records in (POS, longest REF, file order)
order, and the first allele applied to a slot won every conflict. vg trims the shared prefix of a
nested parent's alleles, so in a tandem repeat the parent's POS can lie to the right of its own
child's. The child was then applied first and the parent's allele skipped. The gate's builder
(`call_local.span_haplotypes`, used by `null_decomp.py`) had the same order. No check could fire:
the truth VCFs carry no vg nesting.

**The fix, without changing the rule.** The docstring's rule already says "a nested vg parent
record wins where its change covers a child". The code now does that:

- `nesting_order` keeps the same (POS, longest REF, file order) order, except that a record is
  never applied before a record it is nested in.
- Nesting is a graph identity, read from vg's `>start>end` snarl ID and INFO/AT: C is nested in P
  when both of C's boundary nodes are interior nodes of P's traversals. It is not a coordinate or
  an ID order.
- The conflict rule is unchanged. Records without vg IDs keep exactly the old order.
- The gate's builder uses the same order, and its raw overlap rule no longer assumes POS order.
- `null_decomp.py` now builds haplotypes with the scorer's own rule (`rescore-*` recomputes its
  EDs from the kept VCFs).

**Checks:**

- **The unit test fires on the old scorer.** `TestNesting` writes a parent to the right of its
  child. Scorer version 1 scores it ED 40; version 2 scores it 0.
- **Only reordered regions changed.** Every region whose ED changed has a child written left of a
  parent it overlaps (`descendant_before_ancestor_by_pos`):
  - 3 of 19 for mc and for mc_relabel, 4 of 19 for mc_unchop;
  - 12 of 149 for production;
  - no truth-graph call, and no control anywhere.
- **Truvari is unchanged** in every JSON. The scorer validation (149 regions,
  `stage3_scorer_validation.tsv`) is identical in every ED and verdict.
- **The file-order dependence is gone.** The bcftools-sorted MC VCFs of 5 hotspots now score
  exactly as vg's order does. Version 1 moved L012184 by 25 edits this way.
- **The two builders agree.** The re-run gate's builder and the scorer's raw rule agree in 209
  of 209 VCFs.
- **An independent builder agrees.** The version 2 figures equal every figure the verifier
  reported from its independent order-fixed builder (27,140; −2,003; −4,051; −1,291; and the
  per-region values at L011138, L012272 and L001909).

### What changed in the numbers (hotspot VNTRs; ED to HG002, before → now)

| region | ED_ref | mc | mc_relabel | mc_unchop | production | truth graph |
|---|---|---|---|---|---|---|
| L011138 | 7,502 | 6,207 → **2,561** | 3,095 → **422** | 5,610 → **1,977** | 5,258 → **1,977** | 88 |
| L012272 | 7,022 | 5,597 → **1,311** | 5,641 → **1,355** | 4,485 → **199** | 5,693 → **1,407** | 53 |
| L001909 | 7,652 | 3,974 → **7,469** | 3,279 → **7,503** | 3,989 → **6,612** | 5,351 → **6,373** | 1 |
| L002013 | 2,458 | 5,031 | 5,084 | 5,102 → **5,096** | 5,079 | 36 |
| L009656 | 673 | 2,166 | 2,166 | 3,504 | 2,241 | 0 |
| L014297 | 1,754 | 252 | 257 | 257 | 223 | 0 |
| L015415 | 2,621 | 690 | 690 | 690 | 690 | 0 |
| L005990 | 846 | 6,630 | 6,630 | 4,518 | 6,671 | 0 |
| L012184 | 8,317 | 1,030 | 1,030 | 236 | 1,188 | 763 |
| **sum** | 38,845 | 31,577 → **27,140** | 27,872 → **25,137** | 28,391 → **23,089** | 32,394 → **25,849** | 941 |

At L011138 the parent at chr4:1427914 (REF 11,683 bp, GT 2|1) now applies on both slots. Called
lengths go from 11,862 / 14,260 to 7,925 / 14,769; the truth is 13,717 / 6,591.

**Paired against the local MC call (label − mc on the 9 hotspots):**

| | Δ summed ED, before → now | better / worse at > 0, before → now | at > 10, before → now |
|---|---|---|---|
| mc_relabel | −3,705 (−11.7%) → **−2,003 (−7.4%)** | 2 / 3 → **1 / 4** | 2 / 2 → **1 / 3** |
| mc_unchop | −3,186 (−10.1%) → **−4,051 (−14.9%)** | 4 / 4 → **5 / 3** | 4 / 3 → **5 / 2** |
| production | +817 → **−1,291** | 2 / 6 → **3 / 5** | 2 / 6 → **3 / 5** |

Controls and L016870: no ED changed for any label. Full tables are in
[stage3_null.md](stage3_null.md).

### What changed in the conclusions

1. **A null now passes the fixed rule.** Read mechanically ("summed hotspot ED falls, more
   regions better than worse, no control stratum worse"), mc_unchop passes: 23,089 against
   27,140, 5 better / 3 worse, no control worse. Yet it carries MC's own alignment. mc_relabel
   does not pass (1 / 4).
2. **The resolution limit is larger.** Nulls move summed hotspot ED by 7-15%, not "about 10%".
   MC's own replicates (Issue 2) put its summed hotspot ED anywhere from 22,059 to 27,140:
   - 27,140 for the Stage 3 call;
   - 25,137 with its IDs spread (hybrid200kids);
   - 22,059 with a 50 kb flank (hybrid50k).
3. **Production is now closer to the truth than the local MC call.** By sum: 25,849 against
   27,140. By region count it is still 3 better / 5 worse. The first report's reason for pairing
   candidates with the local MC call ("the local re-map helps") no longer holds. The pairing
   stands for the other reason: the local MC call goes through the same fetch, re-map, flank and
   flags as every candidate, and production does not.
4. **The decomposition numbers moved** (`work/stage3/null/decomp/*/decomp.json`, rescored):
   - **L011138.** The ID effect is −2,139 (2,561 → 422, was −3,112). The unchopped graph given
     MC's own alignments scores 422, the relabel null's value; its own re-placed reads then add
     +1,555 (was +2,515).
   - **L001909.** The relabel difference changes sign, from −695 to +34.
   - **Unchopped graph with MC's alignments.** It now reproduces MC's haplotypes at L005990 and
     L002013. It was L005990 only, because version 1 mis-built L002013.
   - **Unchanged:** with `--depth-term 0` the relabelled and MC graphs give identical haplotypes
     in 5 of 5 regions, and the alignment and VCF-identity checks, which do not use the builder.
   - **What version 1 was hiding.** The `--depth-term 0` MC call at L011138 scores 14,788; version 1
     had scored it 6,574 by dropping its parent's 27 kb allele.
5. **Production over all 149 regions** ([stage3_summary.md](stage3_summary.md)):
   - hotspot VNTRs 130,597 → 124,168;
   - other hotspots 28,271 → 28,376;
   - controls unchanged.
6. **The truth graph is unchanged:** exact in 13 of 19. The misses are L012184 (763), L011138
   (88), L012272 (53), L002013 (36), L001909 (1) and L006556 (310). The first report left
   L002013 and L001909 out of that list.

### The parent/child disagreement: settled by the rule, measured by a sensitivity

After the fix, the parent/child question is what remains of the verifier's case: a child allele
applied on a slot where its parent is also non-reference. The rule applies such a child wherever
the parent's trimmed change does not cover it. That is what the docstring says, so it is kept.
Two new outputs show how much it matters:

- **The count.** `records.applied_under_nonref_parent` counts these child alleles over the 9
  hotspots: 18-26 per MC-derived label and production, 2 for the truth graph.
- **The sensitivity.** `sensitivity.ed_suppress_nested` is the ED when a parent's non-reference
  allele suppresses everything nested in it.

Summed hotspot VNTR ED, rule → suppress:

| mc | mc_relabel | mc_unchop | production | truth graph |
|---|---|---|---|---|
| 27,140 → 25,999 | 25,137 → 24,800 | 23,089 → 21,477 | 25,849 → 24,683 | 941 → 904 |

- **L011138 is the largest case.** The child is a 1,472 bp insertion at chr4:1439597, at the end
  of the parent's trimmed change. The parent's ALT does not contain it. Suppressing it scores
  1,767, against the rule's 2,561.
- **The known double apply is one of these.** The truth graph's double-applied repeat insertion
  at L002013 (first report, observation 2) falls in this class. Under suppress the truth graph is
  exact at L002013 and L001909, so in 15 of 19 regions.
- **The choice matters little for the paired counts.** Against mc, mc_unchop is 5 / 2 under
  suppress and 5 / 3 under the rule; mc_relabel is 1 / 4 either way; production is 2 / 5 and
  3 / 5.

Whether to keep the rule or switch the primary to suppress is a rule change, so it is yours to
make, before the candidates are scored.

## Issue 2: every graph's own noise

**What was wrong.** The first report called the control criterion sound because both MC nulls
left every control unchanged. That shows only that MC is stable. The noise remedy it proposed
(count a candidate as better or worse only outside an ensemble of MC nulls) never measures the
candidates' own noise, and puts the band in the wrong places.

**The fix in the tools:**

- **Replicate arms for any graph** (`call_local.py run --arms ... --arm-graphs ...`), each
  re-mapped and called with the pinned vg and the production flags:
  - `hybrid50k`: 50 kb genome-wide flank, reads re-fetched;
  - `hybrid200kids`: the span's node IDs re-laid by half a spread step (for mc, its native IDs
    spread);
  - `hybrid200kmcids`: the span's nodes placed on mc's own native span IDs, i.e. at mc's ID
    density. Not run for mc, whose layout this is.
- **Every graph gets a numbering change of the kind that moves mc**, native density against
  spread.
- **Checks on the replicates.** The re-laid graph is the same graph up to node IDs (checked at
  L014297 for unit_aware__all: 534 span IDs moved, no node, edge or path different). Re-mapping
  gives identical alignments after ID translation (79,871 of 79,871 reads).
- **The replicate calls.** 310 new calls. All 461 Stage 3, null and replicate calls carry the
  same command: pinned vg, production flags, `--gaf-reads`.
  - 10 `mcids` replicates are impossible: the graph has more span nodes than mc has IDs
    (unit_aware__all at L009656, L007687, L007172, L014111, L006556; unit_aware at L007172,
    L014111, L006556; mafft_linsi and truth at L006556).
  - mafft_linsi has no L012272 graph.
- **`tools/replicate_noise.py`.** For each graph and region it computes the pair edit distance
  between the Stage 3 call's span haplotypes and each replicate's, with the scorer's rule and
  phase oriented to the Stage 3 call, never to the truth. Noise is the largest of these
  distances.
  - Truth EDs of replicates are computed only for graphs already scored: mc, the nulls and the
    truth graph.
  - Output: [stage3_noise.md](stage3_noise.md), `stage3_noise.tsv`.
- **`score_haplotypes.py summarise --noise`** adds a count next to the rule's that accounts for
  noise. A region counts as better (worse) only when the paired difference exceeds both graphs'
  noise summed. By the triangle inequality, no measured replicate of either graph could then
  close the gap.
  - The rule's own counts are untouched.
  - This is an extra reading, not a new rule.

**Summed noise.** The per-stratum sum of each graph's noise, with the largest region:

| graph | hotspot VNTRs (n) | L016870 | matched VNTR controls (7) | non-repeat SV controls (2) |
|---|---|---|---|---|
| mc | 8,727 (9): max 4,237 at L001909 | 0 | 0 | 0 |
| unit_aware__all | 2,291 (9): max 815 at L015415 | 1,238 | 76 (L011430) | 0 |
| mafft_linsi | 936 (8): max 803 at L005990 | 0 | 2 (L004021) | 310 (L006556) |
| unit_aware | 2,580 (9): max 932 at L002013 | 0 | 242 (L007172 183, L007687 57, L004021 2) | 310 (L006556) |
| truth graph | 1,034 (9): L012272 only | 0 | 0 | 0 |
| mc_relabel (flank only) | 4,956 (9) | 0 | 0 | 0 |
| mc_unchop | 11,231 (9): max 4,189 at L011138 | 0 | 0 | 0 |

### What this shows

1. **The control criterion is not sound for the candidates.**
   - At L006556, L011430, L007172 and L007687 every MC call is exact (ED 0 and noise 0), but the
     candidate's two calls differ.
   - So at least one of the candidate's calls is worse than MC there, whatever the truth:
     unit_aware__all at L011430; mafft_linsi at L006556; unit_aware at L006556, L007172 and
     L007687.
   - Whether each candidate passes the fixed control check can therefore depend on an arbitrary
     flank size. The verifier's finding holds, reproduced here and extended to unit_aware.
2. **The node-ID dependence is mostly MC's.**
   - Re-laying IDs moves mc's hotspot haplotypes by up to 2,819 edits (L001909, L011138). It moves
     mc_unchop's by 3,080 (L002013, mc's density).
   - It moves the candidates' by at most 159: 0 everywhere under the half-step, and 1-159 under
     mc's density at L002013, L011138, L012272 and L014297.
   - So mc_relabel is not a like-for-like null for the candidates, and mc's native-ID value at a
     hotspot is a fragile draw (L011138: 2,561 native, 422 spread, 1,617 at 50 kb).
3. **The flank is the candidates' main noise.** At hotspots it is 936 (mafft_linsi) to 2,580
   (unit_aware) summed. MC's flank noise is 8,684.
   - At L015415 and L016870 mc's noise is 0 and unit_aware__all's is 815 and 1,238. An MC-only
     band would have counted that as a real effect.
4. **Reading the MC-derived labels with noise taken into account.** Hotspot VNTRs better / worse
   beyond noise:

   | label | beyond noise | the rule's count (> 0) |
   |---|---|---|
   | mc_relabel | 0 / 0 | 1 / 4 |
   | mc_unchop | 2 / 0 | 5 / 3 |
   | truth graph | 7 / 0 | 9 / 0 |

   - **The relabel null comes out as no effect,** as it should.
   - **The unchop null keeps two reproducible gains:** L005990 −2,112 (both graphs' noise 0) and
     L012272 −1,112 (noise 52 and 43). Merely unchopping MC's graph moves reads (51-639 span reads
     re-placed) and improves those two calls in every replicate. That is a real effect of node
     boundaries, not noise, and a realigned candidate would get it too.
   - **Production** has no replicates, so it has no noise-aware count.
5. **The margins before any candidate is scored.** Summed over the hotspots, the noise margin
   against mc is 11,018 (unit_aware__all), 9,611 (mafft_linsi, 8 regions) and 11,307
   (unit_aware), against mc's summed ED of 27,140. Most of it is at L001909 and L011138 (mc's
   noise 4,237 and 2,823). A candidate can register a noise-aware win or loss only where its
   difference from mc exceeds that region's margin:
   - 0-81 at L015415 and L012272 for mafft_linsi and unit_aware;
   - 2,823-4,476 at L011138 and L001909 for all three.

**Can Issue 2 be fixed in the tools? Only in part.** The tools now measure every graph's noise,
blind, and report counts that account for it. The fixed control criterion ("no control stratum
worse", any difference counts) stays unsound for the candidates. Only changing the decision rule
can fix that.

## What still holds from the first report

- **The gate failure and the reason for stopping.** The re-run gate has the same identity counts:
  - hybrid200k 8 of 19 identical to production;
  - hybrid200korig and gwsub200k 17 of 19;
  - Stage 2 graph 7 of 19.
- **The relabel mechanism.** vg call's depth-term rate windows are keyed on node IDs:
  - identical alignments after ID translation;
  - each graph reproduces its own VCF byte for byte from the other's reads, in 5 of 5 regions;
  - identical calls with `--depth-term 0`, in 5 of 5 regions.
  This is a caller property worth fixing on the vg branch.
- **MC is stable at the controls,** under every builder and every replicate. The candidates are
  not.

## Smaller corrections

- **mafft_linsi has 18 Stage 3 calls, not 19.** There is no L012272 graph, so its hotspot VNTR
  pairing is n = 8, and L012272 is one of the three regions the scorer bug had distorted.
- **The truth-graph miss list** is L012184, L011138, L012272, L002013, L001909 and L006556 (see
  Issue 1).
- **`build_hybrid`'s ID fallback.** Where a graph has more span nodes than there are free IDs in
  the span's range (`ids_interleaved` false; 8 graph/region pairs, all controls or the truth
  graph), its span IDs follow the window's largest ID. The docstring's claim that windows always
  "mix flank and span as they do genome-wide" was false and is corrected. The verifier found the
  calls there unchanged by `--depth-term 0`; I did not re-test that.
- **Crossing records and phase optimisation were never exercised.**
  - No scored label has a record crossing a span edge (0 of 95 pilot JSONs, 0 of 149 production).
  - Every scored call is a single phase block, so phase optimisation never ran.
  - The verifier reports one crossing record in a candidate call (L014111).
- **`--depth-term 0` is not a cleaner arm.** The verifier measured that switching it off moves
  hotspot haplotypes for every graph by thousands of edits summed (mc 7,515, mafft_linsi 10,367).
  It is a different experiment. I did not re-measure this.

## Decisions for you, before any candidate is scored

1. **Nesting clause.** Keep the rule (a child applies outside its parent's trimmed change), or
   make `ed_suppress_nested` primary. At the hotspots it moves summed ED by 1-7% per label and
   changes no paired count by more than one region.
2. **Hotspot criterion.** Options:
   - Keep "any decrease counts". Its mechanical reading passes mc_unchop.
   - Count only differences beyond both graphs' replicate noise. This can be registered now,
     because the noise is measured blind. It reads mc_relabel as no effect and mc_unchop as 2 / 0.
   - Require every replicate of the candidate to beat every replicate of mc. This is tighter; it
     needs the candidates' replicate EDs, so it is scored once the rule is fixed.
3. **Control criterion.** "No control worse" in the Stage 3 call alone is unsound for the
   candidates. Options:
   - Apply it beyond noise. This is lenient where the candidate itself is noisy: L006556 carries
     310 of noise for mafft_linsi and unit_aware.
   - Apply it to every replicate: a candidate fails a control if any of its calls is worse than
     every mc call. This is strict. Because every mc call is exact at L011430, L006556, L007172
     and L007687, the flank noise alone already fails each candidate at one control or more under
     it: unit_aware__all at L011430, mafft_linsi at L006556, unit_aware at three. That follows from
     truth-free distances and mc's known scores, without any candidate ED.
4. **A second baseline.** Also pair candidates against mc_unchop. Unchopping alone gives
   reproducible hotspot changes, so a candidate's win over mc is only partly attributable to
   realignment.
5. **Or reduce the noise at its source.** Key vg call's depth-rate windows on something other
   than node IDs (a caller change, with production equivalence to re-check), or scale to the
   40 hotspots.

## Reproduce

```
export VNTR_STAGE2_WORK=<remap_local.py work dir>
# replicate calls, every graph (the Stage 3 calls in work/stage3/calls/ and diag/<null>@hybrid200k/ already exist)
python3 tools/call_local.py run pilot --graphs '' --arms hybrid50k,hybrid200kids \
    --arm-graphs mc,unit_aware__all,mafft_linsi,unit_aware,truth,mc_unchop --compact --jobs 3 --threads 2
python3 tools/call_local.py run pilot --graphs '' --arms hybrid50k --arm-graphs mc_relabel --compact --jobs 3 --threads 2
python3 tools/call_local.py run pilot --graphs '' --arms hybrid200kmcids \
    --arm-graphs unit_aware__all,mafft_linsi,unit_aware,truth,mc_unchop --compact --jobs 3 --threads 2
# scores (scorer version 2)
for L in mc truth; do python3 tools/score_haplotypes.py batch --label $L --vcf "work/stage3/calls/$L/{id}.vcf.gz" --regions pilot --jobs 5 --force; done
for L in mc_relabel mc_unchop; do python3 tools/score_haplotypes.py batch --label $L --vcf "work/stage3/diag/$L@hybrid200k/{id}.vcf.gz" --regions pilot --jobs 5 --force; done
python3 tools/score_haplotypes.py batch --label genomewide --vcf genomewide --labels-dir genomewide --regions all --jobs 5 --force
python3 tools/replicate_noise.py --regions pilot --out results/stage3_noise --jobs 5      # truth-free for candidates
python3 tools/score_haplotypes.py summarise --labels mc,mc_relabel,mc_unchop,genomewide,truth --baseline mc \
    --regions pilot --noise results/stage3_noise.tsv --out results/stage3_null
python3 tools/score_haplotypes.py summarise --labels genomewide --baseline genomewide --regions all --out results/stage3_summary
python3 tools/null_decomp.py rescore-relabel; python3 tools/null_decomp.py rescore-unchop
python3 tools/call_local.py gate pilot --tsv results/stage3_gate.tsv --md work/stage3/gate/gate_table.md
python3 tools/score_haplotypes.py validate --jobs 3
python3 tools/test_score_haplotypes.py
```

- **Runtime.** The 310 replicate calls took about 55 min at 3 jobs × 2 threads: 36 min for the
  50 kb and ID arms, 18 min for the mc-density arm. Scoring, noise, gate and validation took a
  few minutes each.
- **Disk.** `work/stage3` is 9.6 GB. The replicate arms' reads.gaf are gzipped and their
  hybrid.gfa removed (`--compact`); both can be regenerated.
- **To score the candidates once the rule is settled:**
  `score_haplotypes.py batch --label <graph> --vcf 'work/stage3/calls/<graph>/{id}.vcf.gz' --regions pilot`,
  then `summarise --baseline mc --noise results/stage3_noise.tsv`. For a replicate-based rule,
  also run `replicate_noise.py --truth-labels mc,...,<graph>`.
