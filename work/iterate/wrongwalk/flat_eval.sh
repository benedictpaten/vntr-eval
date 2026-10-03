#!/bin/bash
# The wrong walk at covered repeat sites is the length-weighted mixture: --repeat-sites with --flat-mixture,
# alone and with --depth-count-raw.
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
source work/stage4/chr20/env.sh
export TMPDIR=$PWD/work/tmp
D=work/stage4/chr20; P=$D/patched_st_medoid; VG=work/bin/vg-58faa2261
for arm in rep_flat rep_flatraw; do
  case $arm in rep_flat) X="--flat-mixture";; rep_flatraw) X="--flat-mixture --depth-count-raw";; esac
  echo "[$(date +%T)] call $arm ($X)"
  $VG call $P/graph.gbz -p CHM13#0#chr20 -d 2 -t 8 -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased \
      --repeat-sites $D/repeat_sites_st_medoid.bed $X --gaf-reads $P/reads.gaf.gz > $P/chr20.$arm.vcf 2> $P/chr20.$arm.err
  echo "call exit $?"
  bgzip -f $P/chr20.$arm.vcf && tabix -f -p vcf $P/chr20.$arm.vcf.gz
  python3 tools/patch_contig.py score --vcf $P/chr20.$arm.vcf.gz --label rb_patched_st_medoid_$arm --out $D/score_rb \
      --patch $D/patched_all/patch.json > $D/score_st_medoid_$arm.log 2>&1 || echo SCORE_FAIL
  ( cd $D && PATH=$HOME/.local/bin:/opt/homebrew/bin:$PATH /Users/benedictpaten/PycharmProjects/vg-call-eval/work/truvari-venv/bin/truvari refine -u -a mafft -t 6 \
      -r score_rb/rb_patched_st_medoid_$arm/score/chr20.truvari/candidate.refine.bed -f score_rb/rb_patched_st_medoid_$arm/chr20/chr20.fa \
      score_rb/rb_patched_st_medoid_$arm/score/chr20.truvari > score_rb/rb_patched_st_medoid_$arm/refine.log 2>&1 ) || echo REFINE_FAIL
done
echo "[$(date +%T)] scoring haplotypes"
python3 work/iterate/hapscore_chr20.py --jobs 6
echo "[$(date +%T)] FLAT_EVAL_DONE"
