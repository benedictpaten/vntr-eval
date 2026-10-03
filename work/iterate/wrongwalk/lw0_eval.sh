#!/bin/bash
# --repeat-sites with the linkage model all but off (--linkage-weight 1e-6; --phased refuses 0), so the written
# genotype is the direct pass's; and with only its allele-frequency prior off (--linkage-prior 0).
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
until grep -q DEPTH_EVAL_DONE work/iterate/logs/depth_eval.out 2>/dev/null; do sleep 30; done
source work/stage4/chr20/env.sh
export TMPDIR=$PWD/work/tmp
D=work/stage4/chr20; P=$D/patched_st_medoid; VG=work/bin/vg-58faa2261
for arm in rep_lwtiny rep_lp0; do
  case $arm in rep_lwtiny) X="--linkage-weight 0.000001";; rep_lp0) X="--linkage-prior 0";; esac
  echo "[$(date +%T)] call $arm ($X)"
  $VG call $P/graph.gbz -p CHM13#0#chr20 -d 2 -t 8 -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased \
      --repeat-sites $D/repeat_sites_st_medoid.bed $X --gaf-reads $P/reads.gaf.gz > $P/chr20.$arm.vcf 2> $P/chr20.$arm.err
  rc=$?; echo "call exit $rc"
  [ $rc = 0 ] && bgzip -f $P/chr20.$arm.vcf && tabix -f -p vcf $P/chr20.$arm.vcf.gz
done
python3 work/iterate/hapscore_chr20.py --jobs 6
echo "[$(date +%T)] LINKAGE_ARMS_DONE"
