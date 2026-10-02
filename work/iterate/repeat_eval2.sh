#!/bin/bash
# After repeat_eval.sh: --repeat-descent and --repeat-linkage on the same graph, reads and BED; then every arm's
# per-region bootstrap and called-haplotype edit distance.
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
until grep -q REPEAT_EVAL_DONE work/iterate/logs/repeat_eval.out 2>/dev/null; do sleep 30; done
source work/stage4/chr20/env.sh
export TMPDIR=$PWD/work/tmp
D=work/stage4/chr20; P=$D/patched_st_medoid; VG=work/bin/vg-58faa2261
for arm in desc link; do
  if [ $arm = desc ]; then X="--repeat-descent"; else X="--repeat-linkage"; fi
  echo "[$(date +%T)] call $arm"
  $VG call $P/graph.gbz -p CHM13#0#chr20 -d 2 -t 8 -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased \
      --repeat-sites $D/repeat_sites_st_medoid.bed $X --gaf-reads $P/reads.gaf.gz > $P/chr20.$arm.vcf 2> $P/chr20.$arm.err
  echo "call exit $?"
  bgzip -f $P/chr20.$arm.vcf && tabix -f -p vcf $P/chr20.$arm.vcf.gz
  grep "repeat " $P/chr20.$arm.err
  python3 tools/patch_contig.py score --vcf $P/chr20.$arm.vcf.gz --label rb_patched_st_medoid_$arm --out $D/score_rb \
      --patch $D/patched_all/patch.json > $D/score_st_medoid_$arm.log 2>&1 || echo SCORE_FAIL
  ( cd $D && PATH=$HOME/.local/bin:/opt/homebrew/bin:$PATH /Users/benedictpaten/PycharmProjects/vg-call-eval/work/truvari-venv/bin/truvari refine -u -a mafft -t 6 \
      -r score_rb/rb_patched_st_medoid_$arm/score/chr20.truvari/candidate.refine.bed -f score_rb/rb_patched_st_medoid_$arm/chr20/chr20.fa \
      score_rb/rb_patched_st_medoid_$arm/score/chr20.truvari > score_rb/rb_patched_st_medoid_$arm/refine.log 2>&1 ) || echo REFINE_FAIL
done
echo "[$(date +%T)] scoring haplotypes"
python3 work/iterate/hapscore_chr20.py --jobs 6
cd work/iterate/perregion && python3 region_metrics.py > /dev/null && cd ../../..
echo "[$(date +%T)] REPEAT_EVAL2_DONE"
