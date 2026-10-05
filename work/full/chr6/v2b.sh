#!/bin/bash
# Continuation of v2.sh from the databases (GAF header stripped, as stage4_chr6/run_chr6b.sh), then call and score.
set -u -o pipefail
R=/Users/benedictpaten/PycharmProjects/vntr-eval; cd $R
export PATH=$HOME/.local/bin:/opt/homebrew/bin:$PATH TMPDIR=$R/work/tmp
V=work/full/chr6/v2; A=$V/patched; D=work/stage4_chr6/chr6; E=/Users/benedictpaten/PycharmProjects/vg-call-eval
VG=work/bin/vg-nestfix4
step() { echo "[$(date +%T)] $*"; }
step dbs
gzip -dc $A/reads.gaf.gz | grep -v '^@' | bgzip -@ 4 > $A/reads.nohdr.gaf.gz
rm -f $A/fifo; mkfifo $A/fifo
gaf-base sort $A/reads.nohdr.gaf.gz -o $A/fifo -p > $V/sort.log 2>&1 &
gaf-base construct $A/fifo -r $A/graph.gbz -o $A/reads.db --overwrite > $V/gafdb.log 2>&1 || { echo GAFDB_FAILED; exit 1; }
wait; rm -f $A/fifo $A/reads.nohdr.gaf.gz
step call
$VG call $A/graph.gbz -p CHM13#0#chr6 -d 2 -t 8 -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased \
  --gaf-base $A/reads.db --gbz-base $A/graph.db > $A/chr6.vcf 2> $V/call.log || { echo CALL_FAILED; exit 1; }
bgzip -f $A/chr6.vcf && tabix -f -p vcf $A/chr6.vcf.gz
step score
python3 tools/patch_contig.py score --contig chr6 --vcf $A/chr6.vcf.gz --label pgrealign_v2 --out $D/score --patch $A/patch.json --threads 4 > $V/score.log 2>&1 || { echo SCORE_FAILED; exit 1; }
for L in pgrealign_v2 unpatched; do
  [ -s $D/score/$L/score/chr6.truvari/refine.variant_summary.json ] && continue
  step "refine $L"
  ( cd $D/score && $E/work/truvari-venv/bin/truvari refine -u -a mafft -t 4 -r $L/score/chr6.truvari/candidate.refine.bed -f $L/chr6/chr6.fa $L/score/chr6.truvari > $L/refine.log 2>&1 ) || echo "$L REFINE_FAILED"
done
echo V2B_DONE
