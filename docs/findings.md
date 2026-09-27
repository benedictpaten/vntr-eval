# Short-read SV errors in VNTRs: where they are, what the graph does there, and what to change

HG002, ~30x Illumina, the hap32 graph (32 sampled HPRC v2.1 haplotypes + CHM13 + GRCh38), `vg call`
at the PR default (`work/wgs-mm095`), T2T-Q100 v1.1. 2026-09-25. Every number below comes from the
files named in *Sources*; the workflow run is `wf_679112d3-b8a`.

## 1. Only a small subset of VNTRs produces the false positives

The denominator is every VNTR region in the SV benchmark, not only those where the truth has an SV:
137,197 regions of period >= 7 (28.1 Mb). RepeatMasker could not supply them on its own -- its
Simple_repeat track has no motif over 30 bp and misses 62% of the truth's TRF VNTR spans -- so the
catalogue adds a truth-independent tandem-repeat finder (periods 7-2,000 bp; it overlaps 95.1% of the
truth-TRF VNTR spans).

- VNTRs are **1.03% of benchmarked sequence but hold 81% of vg's SV false positives** (10,174 of
  12,553), 71% of its false negatives and 82% of PanGenie's false positives. Outside tandem repeats
  both callers score SV F1 about 0.93; inside VNTRs vg scores 0.455 and PanGenie 0.462.
- **2.8% of VNTR regions have any vg FP, 0.9% have three or more, 0.07% have ten or more.** Half the
  VNTR FPs sit in 796 regions (0.58%, 1.96 Mb) and 80% in 2,068 (1.5%).
- The 1,231 regions with >= 3 FPs cover 2.5 Mb -- 0.09% of the benchmarked genome -- and hold
  **49.9% of all genome-wide vg SV FPs** and 33% of FNs.
- FNs and PanGenie's FPs concentrate in the same regions (0.52-0.55% of regions hold half of each;
  Spearman 0.79 between vg and PanGenie FP counts; 88.7% of vg's >= 3-FP regions also have a
  PanGenie FP).
- **The subset is predictable without the truth.** It is long (median 1.4 kb against 74 bp for
  FP-free VNTRs), high-copy (35 against 2.6 copies), GC-rich (0.60 against 0.43) and subtelomeric
  (59% within 5 Mb of a telomere, against 15%). Length alone separates it with AUC 0.97: 0.03% of
  VNTRs under 50 bp have an FP, 43% at 1-2 kb, 85% at >= 5 kb. 93% of the subset carries a truth SV.
  Graph complexity adds signal at fixed length (AUC 0.75-0.85).
- A truth-free target set for production: benchmarked VNTRs >= 1 kb, 3,111 regions (6.1 Mb), hold
  46% of all vg SV FPs and 40% of FNs (>= 500 bp: 16,650 regions, 64% and 53%).

A separate 22% of VNTR FNs sit in regions with no FP. They are small (median 75 bp) and 62% have a
20-49 bp vg indel within 100 bp: near-threshold size misses, not a graph problem.

## 2. The graph does not tangle there -- it misaligns the panel

Compared on 40 FP-dense VNTR hotspots against 40 FP-free VNTRs matched on length and period:

- **No cycles.** Zero node revisits by any haplotype or by CHM13 in 3,148 regions and in all 150
  stratified loci; zero directed cycles. (The first test's "1,242 fragments" at chr20:1.90 Mb was an
  artefact of querying with `--context 0`, which cuts haplotypes wherever they leave the reference.)
- **Hotspots are giant multi-allelic bubbles of tiny nodes.** 848 against 113 nodes per kb (mean node
  3.6 bp against 11), 30 against 7 distinct haplotype walks, 26 against 3 alleles in the largest
  bubble. 32/40 hotspots are one giant bubble, 8/40 a chain of fragments; 32/40 controls are simple.
- **The alignment the graph implies is far from optimal.** For every pair of haplotypes, the edit
  cost allowed by the nodes they share is **1.56x the pairwise optimum at hotspots, 1.09x at matched
  controls and 1.00x at non-repeat SVs** (223 against 3.9 excess edits per kb; AUC 0.94; 1.64 against
  1.03 under affine scoring). Each haplotype becomes about 9 SV-sized pieces against 1. An L-INS-i
  MSA stacks 93% of pairwise-identical bases on shared columns where the graph stacks 80% (chr10).
- **What the misalignment looks like** (deep dives): repeat expansions placed at different ends of
  the array in different haplotypes, repaid by compensating deletions (chr10: a haplotype 9 bp longer
  than CHM13 and 384 edits from it is encoded as 1,430 bp inserted plus 1,421 deleted; one 78 bp unit
  carried by 31 haplotypes is spelled by three different node paths); a whole family of haplotypes
  routed through a parallel 1,254 bp copy that is 96.7% identical to the CHM13 sequence it bypasses,
  including 542 bp of unique flank (chr4); a block of ~38 units inserted near the start of the array
  in every non-CHM13 haplotype and deleted again 3' (chr17). Also N-gap haplotypes and paths ending
  inside the array.
- **FPs sit in the big bubbles.** In the hotspot spans, the 78 bubbles with >= 21 alleles are 4% of
  bubbles but hold 58% of vg's FPs; bubbles under 50 bp hold none.

## 3. Representation, not only genotyping, caps the score

- **Is HG002's allele in the graph?** At hotspots, 57% of HG002 haplotypes are within 3% of some panel
  haplotype (20% exact); 27% have none. At matched controls, 100% are within 1%. Both alleles are
  present at 19/40 hotspots, one absent at 8, both absent at 13. Checking six absent-allele loci
  against the full 459-haplotype HPRC v2.1 panel brought only 1 of 11 missing haplotypes within 3%.
- **Scoring cannot credit a better answer as things stand.** The exact truth, written the way this
  graph decomposes it, scores truvari F1 0.271 at the hotspots (vg 0.246) and reproduces 93% of vg's
  FP+FN. After `truvari refine` the exact truth scores 0.892, the best panel walk 0.487, PanGenie 0.323,
  vg 0.248. **Any graph change here has to be judged on haplotype edit distance and refined F1**, with
  raw bench F1 secondary.
- **vg also genotypes these sites badly.** It builds recombinant, over-long alleles from dozens of
  nested child snarls: its called haplotype matches the closest panel walk in 4 of 80 cases and is
  further from the truth than CHM13 in 30 of 80. Summed edit distance to truth over the hotspots:
  calling nothing 194,796; vg 130,860; PanGenie 124,891; the best panel walk 15,586. At chr4 vg's
  called allele is 360 units, longer than any panel haplotype (335); at chr17 it is 201 units against
  a truth of 69. Its depth-term rate is ~6x too low in these dense graphs because it divides by graph
  bp rather than haplotype bp, which favours expansions.

## 4. Reads reach the regions but cannot be placed

At hotspot cores against matched controls: depth 0.62x the flank (0.92x); 74% of reads MAPQ < 5 (19%);
so the effective depth under vg's MAPQ weighting is 6.6x against 24x. 97% of the ambiguity is local
(MAPQ-0 reads rarely have a copy elsewhere in CHM13). Of the MAPQ < 5 reads, **31% are redundant --
the same sequence on different nodes in different panel haplotypes, which a better-aligned graph would
remove -- and 56% are tandem ambiguity, which it would not.** 44% of hotspot core reads end in a
>= 20 bp low-quality tail (13% at controls): sequencing degradation, not missing alleles.

## 5. Strategies, ranked

1. **Re-align the alleles inside each multi-allelic VNTR bubble and splice the bubble back** between
   anchor nodes that every path visits once (the toolkit finds them). Align the spanning haplotype
   sequences with a unit-aware MSA or POA (mafft L-INS-i takes ~25 s for 36 x 7 kb), induce a graph
   from the MSA, re-thread each path as its row. Fixes the misalignment (the representation FPs and the
   redundant read class), not absent alleles or tandem ambiguity. Risks: repeat MSAs have many equal
   optima; differently ordered arrays cannot be stacked acyclically; column induction can recreate
   1-bp mesh; new node IDs mean re-indexing.
2. **Genotype each problematic VNTR as one site of whole alleles** (caller-side companion): panel
   alleles, optionally +-k units, scored as whole alleles with a per-haplotype-bp depth term and no
   nested-child recombination. Offline whole-allele scoring picked the right allele pair at chr4 and
   ranked the truth first at chr10; the best panel pair would take refined F1 from 0.25 to 0.49.
3. **First repair gross construction defects**, found by a genome-wide detector: parallel non-reference
   runs >= 90% identical to the reference they bypass, N-gap haplotypes, paths ending inside a VNTR.
   10/40 hotspots are under-aligned this way against 1/40 controls, and it is the cheapest read-side win.
4. **Do it upstream in Minigraph-Cactus** once a pilot works: mask catalogue VNTRs (e.g. >= 1 kb)
   before Cactus and re-insert per-haplotype alleles as one bubble aligned as in (1).
5. Adding haplotypes or changing sampling: low priority (see section 3). Untested: whether hap32's
   recombinant sampling puts switches inside VNTRs and creates chimeric alleles.
6. Acyclicity work: unnecessary (no cycles). Collapsing VNTRs into cycles: not recommended.

**A pilot needs no `vg call` until its last stage.** On ~15 hotspot loci plus all VNTR controls, the
graph-only metrics (cost/optimum, excess edits per kb, pieces per haplotype), "truth written by the
graph" raw and refined F1, and read-placement redundancy on the existing reads can all be computed
before and after re-aligning a bubble; targets from the controls are cost/optimum <= 1.1 and
inflation <= 1.1. Then a local re-map of the same reads, then `vg call` on the local graph scored by
haplotype edit distance and refined F1, with the VNTR and non-TR controls as must-not-regress.
Suggested loci: alleles present but misaligned -- L009656 chr21:41,564,799-41,568,906, L014297
chr6:171,505,695-171,508,299, L015415 chr8:770,452-772,838 (both HG002 haplotypes exactly in the
panel), L005990 chr17; under-aligned -- L012184 chr4:191,498,564-191,501,883, L011138, L012272;
alleles absent -- L002013 chr10:133,306,611-133,310,510; negative control -- L016870 chrX.

## 6. A correction to an existing page

`docs/pangenie-comparison.md` says PanGenie and `vg call` ran on the same hap32 graph. They did not:
PanGenie genotyped the HPRC v2.1 Minigraph-Cactus leave-one-out VCF (the file name says so, and its
variant IDs use different node IDs from hap32 at the same position -- e.g. chr1:5,546 is on node
663079 in PanGenie's IDs and 661901 in hap32). So PanGenie had the full panel minus HG002, not 34
haplotypes. The comparison is same reads and same truth, but not the same graph.

## 7. How to look at a region

```bash
export PATH=~/.local/bin:$PATH
V=/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-2026-07-vg-call-refactor/e061d8cf-7996-49a4-986b-3686374ef250/scratchpad/vntr
# a self-contained HTML page for any locus: MSA of CHM13, panel walks, truth and calls; Bandage graph
# and per-walk node track; calls with truvari status; read pileup by walk
python3 $V/viewer/view_region.py L009656                     # by census locus id
python3 $V/viewer/view_region.py chr16:864311-870232 --mafft linsi
# the data behind it (1-based inclusive)
python3 $V/lib/vntr_region.py chr21 41564800 41568906 OUTDIR --msa --bandage
```

Pages already built: `$V/viewer/chr4_191498564_191501883.html`, `chr17_369877_372692.html`,
`chr10_133306611_133310510.html` (hotspots); `chr19_4017784_4020263.html`, `chr4_39782252_39784838.html`,
`chr18_80168214_80171889.html` (their matched controls); `chr1_80566080_80566608.html` (non-TR control);
`chr20_1902001_1904000.html`. In the node track, colour by node identity in MSA-column mode: identical
columns sitting on different nodes are the misalignment.

## Sources

Census: `$V/census/census.md`, `vntr_regions.tsv`, `loci.tsv`, `strata.tsv`. Topology:
`$V/topology/topology.md`. Truth in the graph: `$V/representability/representability.md`. Reads:
`$V/reads/reads.md`. Deep dives: `$V/deep1/deep.md` (chr17), `$V/deep2/deep.md` (chr4),
`$V/deep3/deep.md` (chr10). Figures: `$V/deep2/msa_vs_graph.png` (chr4), `$V/deep1/alignment_tracks.png`
(chr17), `$V/deep3/haplotype_paths.png` (chr10).
