#!/bin/bash
# Continuation of run_chr6.sh from the per-arm databases (GAF header stripped: its graph name is the
# giraffe index graph's, which gaf-base construct refuses to match against graph.gbz).
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
source work/stage4_chr6/env.sh
S=work/stage4_chr6; D=$S/chr6; E=/Users/benedictpaten/PycharmProjects/vg-call-eval
VG=work/bin/vg-2a6a228a5
step() { echo "[$(date +%T)] $*"; }
for arm in unpatched patched; do
  if [ ! -s $D/$arm/reads.gaf.gz ]; then
    step "map $arm"
    python3 tools/patch_contig.py map --contig chr6 --gfa $D/$arm/graph.gfa --reads $D/reads --out $D/$arm --threads 8 > $S/map_$arm.log 2>&1 || { echo MAP_FAILED; exit 1; }
  fi
  step "dbs $arm"
  gzip -dc $D/$arm/reads.gaf.gz | grep -v '^@' | bgzip -@ 4 > $D/$arm/reads.nohdr.gaf.gz
  [ -s $D/$arm/graph.db ] || gbz-base construct $D/$arm/graph.gbz -o $D/$arm/graph.db --overwrite > $S/gbzdb_$arm.log 2>&1
  rm -f $D/$arm/fifo; mkfifo $D/$arm/fifo
  gaf-base sort $D/$arm/reads.nohdr.gaf.gz -o $D/$arm/fifo -p > $S/sort_$arm.log 2>&1 &
  gaf-base construct $D/$arm/fifo -r $D/$arm/graph.gbz -o $D/$arm/reads.db --overwrite > $S/gafdb_$arm.log 2>&1 || { echo GAFDB_FAILED; exit 1; }
  wait; rm -f $D/$arm/fifo $D/$arm/reads.nohdr.gaf.gz
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
