# Repeat sites in vg call (parked)

**Status:** parked on 2026-10-02. The code is on branch
[`repeat-sites`](https://github.com/benedictpaten/vg/tree/repeat-sites) of benedictpaten/vg (head
`7d8c48d30`, based on the PR #4990 head `283808454`). It is not part of the PR. TAP `18_vg_call.t`
passes, 468/468. The measurements are in
[results/stage4_chr20.md](../results/stage4_chr20.md), sections 4y, 4z and 4aa.

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
