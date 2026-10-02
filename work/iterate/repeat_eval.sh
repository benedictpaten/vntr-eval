#!/bin/bash
# --repeat-sites on chr20: the medoid-star patched graph called with each patched repeat region as one site,
# scored like the other arms. Waits for the identity gate (same binary, no option) to finish first.
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
G=/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-2026-07-vg-call-refactor/e061d8cf-7996-49a4-986b-3686374ef250/scratchpad/gate
until [ -f $G/rep.done ]; do sleep 30; done
echo "[$(date +%T)] gate: $(cat $G/rep.done)"
/usr/bin/cmp $G/new2.vcf $G/rep.vcf && echo "GATE VCF IDENTICAL"
/usr/bin/cmp $G/new2.anchors.tsv $G/rep.anchors.tsv && echo "GATE ANCHORS IDENTICAL"
/usr/bin/cmp $G/new2.mosaic.tsv $G/rep.mosaic.tsv && echo "GATE MOSAIC IDENTICAL"
source work/stage4/chr20/env.sh
export TMPDIR=$PWD/work/tmp
D=work/stage4/chr20; P=$D/patched_st_medoid; VG=work/bin/vg-2a287c328
echo "[$(date +%T)] call"
$VG call $P/graph.gbz -p CHM13#0#chr20 -d 2 -t 8 -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased \
    --repeat-sites $D/repeat_sites_st_medoid.bed --gaf-reads $P/reads.gaf.gz > $P/chr20.rep.vcf 2> $P/chr20.rep.err
echo "call exit $?"
bgzip -f $P/chr20.rep.vcf && tabix -f -p vcf $P/chr20.rep.vcf.gz
grep "repeat sites" $P/chr20.rep.err
echo "[$(date +%T)] score"
python3 tools/patch_contig.py score --vcf $P/chr20.rep.vcf.gz --label rb_patched_st_medoid_rep --out $D/score_rb \
    --patch $D/patched_all/patch.json > $D/score_st_medoid_rep.log 2>&1 || echo SCORE_FAIL
( cd $D && PATH=$HOME/.local/bin:/opt/homebrew/bin:$PATH /Users/benedictpaten/PycharmProjects/vg-call-eval/work/truvari-venv/bin/truvari refine -u -a mafft -t 6 \
    -r score_rb/rb_patched_st_medoid_rep/score/chr20.truvari/candidate.refine.bed -f score_rb/rb_patched_st_medoid_rep/chr20/chr20.fa \
    score_rb/rb_patched_st_medoid_rep/score/chr20.truvari > score_rb/rb_patched_st_medoid_rep/refine.log 2>&1 ) || echo REFINE_FAIL
echo "[$(date +%T)] REPEAT_EVAL_DONE"
