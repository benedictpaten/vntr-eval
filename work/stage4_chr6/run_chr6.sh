#!/bin/bash
# Stage 4e, held-out chr6: the chr20 winner unchanged -- every eligible TR (>= 50 bp haplotype length
# variation) realigned with abPOA defaults on hap32, patched, re-mapped, called, scored -- vs unpatched.
# vg call reads alignments from a per-arm GAF-Base (--gaf-reads would need ~40 GB RAM on chr6).
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
source work/stage4_chr6/env.sh
S=work/stage4_chr6; D=$S/chr6; E=/Users/benedictpaten/PycharmProjects/vg-call-eval
VG=work/bin/vg-2a6a228a5
CAT=/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-2026-07-vg-call-refactor/e061d8cf-7996-49a4-986b-3686374ef250/scratchpad/vntr/census/catalogue.tsv
step() { echo "[$(date +%T)] $*"; }
step scan
python3 tools/scan_tr_regions.py --catalogue $CAT --contig chr6 --gbz $E/work/wgs-mm095/chr6/chr6.gbz --out $S/scan.tsv --bed $S/prefilter.bed > $S/scan.log 2>&1 || { echo SCAN_FAILED; exit 1; }
step "package $(wc -l < $S/prefilter.bed) eligible"
python3 tools/package_regions.py build --loci $S/prefilter.bed --out $S/regions --work $S/pkgwork --jobs 6 > $S/package.log 2>&1
IDS=$(for d in $S/regions/*/; do [ -s "$d/hap32.fa" ] && basename "$d"; done)
step "realign $(echo $IDS | wc -w)"
python3 tools/realign.py poa_abpoa $IDS --jobs 3 --threads 1 > $S/realign.log 2>&1
for i in $IDS; do [ -s $S/candidates/poa_abpoa/$i.gfa ] && echo $i; done > $S/patch_ids.txt
step "graphs ($(wc -l < $S/patch_ids.txt) candidates)"
python3 tools/patch_contig.py gfa --contig chr6 --out $D > $S/gfa.log 2>&1 || { echo GFA_FAILED; exit 1; }
python3 tools/patch_contig.py patch --contig chr6 --out $D/unpatched --ids none > $S/patch_unpatched.log 2>&1 || { echo PATCH_FAILED; exit 1; }
python3 tools/patch_contig.py patch --contig chr6 --out $D/patched --ids $S/patch_ids.txt --method poa_abpoa > $S/patch_patched.log 2>&1 || { echo PATCH_FAILED; exit 1; }
step reads
python3 tools/patch_contig.py reads --contig chr6 --out $D/reads --jobs 4 > $S/reads.log 2>&1 || { echo READS_FAILED; exit 1; }
for arm in unpatched patched; do
  step "map $arm"
  python3 tools/patch_contig.py map --contig chr6 --gfa $D/$arm/graph.gfa --reads $D/reads --out $D/$arm --threads 8 > $S/map_$arm.log 2>&1 || { echo MAP_FAILED; exit 1; }
  step "dbs $arm"
  rm -f $D/$arm/fifo; mkfifo $D/$arm/fifo
  gbz-base construct $D/$arm/graph.gbz -o $D/$arm/graph.db --overwrite > $S/gbzdb_$arm.log 2>&1
  gaf-base sort $D/$arm/reads.gaf.gz -o $D/$arm/fifo -p > $S/sort_$arm.log 2>&1 &
  gaf-base construct $D/$arm/fifo -r $D/$arm/graph.gbz -o $D/$arm/reads.db --overwrite > $S/gafdb_$arm.log 2>&1 || { echo GAFDB_FAILED; exit 1; }
  wait; rm -f $D/$arm/fifo
  step "call $arm"
  $VG call $D/$arm/graph.gbz -p CHM13#0#chr6 -d 2 -t 8 -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased \
    --gaf-base $D/$arm/reads.db --gbz-base $D/$arm/graph.db > $D/$arm/chr6.vcf 2> $S/call_$arm.log || { echo CALL_FAILED; exit 1; }
  bgzip -f $D/$arm/chr6.vcf && tabix -f -p vcf $D/$arm/chr6.vcf.gz
done
step score
for x in "patched:$D/patched/chr6.vcf.gz" "unpatched:$D/unpatched/chr6.vcf.gz" "production:$E/work/wgs-mm095/chr6/chr6.vcf.gz"; do
  python3 tools/patch_contig.py score --contig chr6 --vcf ${x#*:} --label ${x%%:*} --out $D/score --patch $D/patched/patch.json --threads 2 > $S/score_${x%%:*}.log 2>&1 &
done
wait
step RUN_CHR6_DONE
