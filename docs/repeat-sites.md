# Repeat sites in vg call (parked)

**Status:** parked on 2026-10-02. The code is on branch
[`repeat-sites`](https://github.com/benedictpaten/vg/tree/repeat-sites) of benedictpaten/vg (head
`7d8c48d30`, based on the PR #4990 head `283808454`). It is not part of the PR. TAP `18_vg_call.t`
passes, 468/468. The measurements are in
[results/stage4_chr20.md](../results/stage4_chr20.md), sections 4y, 4z and 4aa.

**Rebased 2026-10-05** onto the PR head `0575c3825`, as branch
[`repeat-sites-rebased`](https://github.com/benedictpaten/vg/tree/repeat-sites-rebased) (head `60d166b5f`). Only
the TAP plan count conflicted. TAP `18_vg_call.t` passes, 476/476, and without the options the binary's records and
mosaic are identical to `0575c3825`'s. It is still parked and not part of the PR. See
[Rebased, on the latest chr20 graph](#rebased-on-the-latest-chr20-graph-2026-10-05).

## What it is

`vg call --repeat-sites BED` genotypes each BED region as one site.
- **The site.** A fake snarl spanning the region's top-level snarls, which must be consecutive links of
  one top-level chain.
- **The alleles.** The panel haplotypes' whole walks across the region.
- **Reads, scoring and output.** Reads are mapped to the fine-grained graph as before. Scoring and the
  linkage model are unchanged. The output is one record per difference, through the block-record code.

Two variants test alternatives:
- `--repeat-descent` descends into the snarls nested in the chosen walks.
- `--repeat-linkage` keeps per-snarl sites but lets the linkage model allow no switch inside a
  region.

## What was found (chr20 HG002, short reads, medoid-star graph)

**Repeat sites work.**
- Refined SV F1 goes from 0.6957 to 0.7204. Per-region bootstrap of refined FP+FN: −49 [−93, −10].
- Small variants are unchanged.
- Haplotypes that no panel haplotype spells fall from 182 to 29.

**Neither variant beats plain repeat sites.**
- Descent gives part of the gain back (0.7165): it re-opens stitching inside the chosen walks.
- Linkage alone is not enough (0.7069): each snarl is still genotyped on its own.

**The whole-allele prototype does not beat vg call** once six regions the scorer cannot read are left
out. Those six lie inside a 474 kb record that encloses them (see 4z). Without them, repeat sites score
45,545 edits, the prototype 52,190 and per-snarl calling 55,270.

**Where repeat sites apply (490 regions), no stitching is left.** What remains is the wrong panel walk:
30,352 edits above the best panel pair, almost all copy number. Section 4aa breaks it down:

| cause | edits |
|---|---|
| the true allele is not in the panel; vg picks the walk the reads fit best | 12,616 |
| TR755508 (telomere): a right allele is not among vg's alleles | 9,379 |
| real model failures (9 sites) | 3,405 |
| other | 4,952 |

**The real model failures come from MAPQ.**
- Reads from the expanded allele map ambiguously within the repeat. They get low MAPQ (median
  mismapping probability 0.79), so each counts at most ln(1/e_r).
- High-MAPQ reads that a homozygote fits on both haplotypes cost the heterozygote ln 2 each.
- Read placement, the linkage model and the depth term were each ruled out.

**No global option fixes it.** Each fixes some sites and breaks others. The best single option is
`--depth-count-raw`, which gives refined SV F1 0.7304.

## Rebased, on the latest chr20 graph (2026-10-05)

The latest chr20 graph has every tandem-repeat region realigned with pgrealign: 14,511 regions. HG002 short reads
were mapped to its 32 sampled haplotypes and called with the PR head's flags. Each arm reuses the same reads and
the same GAF-Base. The BED is the realigned regions, from left anchor start to right anchor end, with touching
regions merged: 14,147 intervals. Repeat sites took 13,420 of them, covering 340,951 top-level snarls. The other 727
hold no top-level snarl on the reference.

| arm | ALL F1 | indel F1 | SV F1 raw | SV F1 refined | called-haplotype edit distance | off-panel strands | call time |
|---|---|---|---|---|---|---|---|
| per-snarl (PR head) | 0.9738 | 0.9302 | 0.6403 | 0.7066 | 53,805 | 389 | 271 s |
| `--repeat-sites` | 0.9730 | 0.9286 | 0.6412 | 0.7336 | 42,297 | 58 | 360 s |
| `--repeat-descent` | 0.9730 | 0.9284 | 0.6399 | 0.7036 | 49,413 | 140 | 484 s |
| `--repeat-linkage` | 0.9736 | 0.9296 | 0.6421 | 0.7086 | 49,737 | 264 | 272 s |

Edit distance is summed over the 27,971 called strands of 13,990 regions inside the truth BED. Each arm is compared
with per-snarl calling by a paired per-region bootstrap:
- **Repeat sites cut haplotype error.** Edit distance falls by 11,508 [−21,094, −3,607].
- **The refined SV gain cannot be told from zero.** SV FP+FN inside the regions changes by −19 [−73, +26]; raw SV
  F1 is flat.
- **They cost small variants.** ALL F1 falls 0.0008 and indel F1 0.0016. Whether that is the block records'
  representation or real genotype changes has not been measured.
- **`--repeat-linkage` gets a third of the edit-distance gain**: −4,068 [−7,918, −983], with no SV or small-variant
  change. Forbidding switches is not what makes repeat sites work. What does is one joint genotype over whole panel
  walks, with each read scored once against whole walks and a depth term over the whole region.
- **ONT was not run.** Its recombinant calls are mostly correct (HG002's allele is often off-panel), so repeat sites,
  which can only return panel walks, are not expected to help it.

Drivers: `work/full/chr20/repeat_sites/run.sh` (calls and scores), `compare.py` (tables and bootstraps).

## Open when this is picked up again

1. **Do not count ambiguity between copies inside a repeat site as mismapping.** A read whose low MAPQ
   comes only from alternative placements within the site belongs to the site either way. This targets
   the 9 real failures and is untested.
2. **Repeat sites inside a larger snarl.** 71 of the 78 regions that are not sites lie in one 1.06 Mb
   top-level snarl at chr20:64.1-65.2 Mb. A nested site could recover at most 4,848 edits.
3. **Alleles the panel lacks** cost the most (12,616). Only a richer panel, or alleles built from the
   reads, can fix that.
4. **Held-out test.** Everything above is chr20. chr6 has not been run with repeat sites.

## How to reproduce

- **The chr20 arms.** `work/iterate/repeat_eval.sh` and `repeat_eval2.sh` run them.
- **The wrong-walk analysis.** It is in `work/iterate/wrongwalk/`.
- **The diagnostic dump.** It needs the patch `wrongwalk/diag_dump.patch` applied to the branch. The
  patch makes `--dump-likelihoods` also write each allele's length and every diploid genotype's read
  and depth terms. A binary built with it writes a byte-identical VCF.
