#!/bin/bash
# The parked repeat-sites branch, rebased onto the PR head 0575c3825 (vg 60d166b5f), on the latest chr20 graph with every
# region realigned (test/banded: same 32 sampled haplotypes, short-read GAF-Base and calling flags as its baseline call).
# Arms: base (no option; must match the baseline's records), --repeat-sites, + --repeat-descent, --repeat-linkage.
set -u -o pipefail
R=/Users/benedictpaten/PycharmProjects/vntr-eval; E=/Users/benedictpaten/PycharmProjects/vg-call-eval; cd $R
export PATH=$HOME/.local/bin:/opt/homebrew/bin:$PATH TMPDIR=$R/work/tmp
VG=/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-2026-07-vg-call-refactor/e061d8cf-7996-49a4-986b-3686374ef250/scratchpad/vg-rs/bin/vg
C=work/full/chr20; D=$C/test/banded; O=$C/repeat_sites; BED=$O/realigned_all.bed
step() { echo "[$(date +%T)] $*"; }
for arm in base rep desc link; do
  case $arm in base) X=();; rep) X=(--repeat-sites $BED);; desc) X=(--repeat-sites $BED --repeat-descent);; link) X=(--repeat-sites $BED --repeat-linkage);; esac
  step "call $arm ${X[*]:-}"
  if [ ! -s $O/$arm.vcf.gz ]; then
    /usr/bin/time -l $VG call $D/sampled.gbz -p CHM13#0#chr20 -d 2 -t 4 -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased \
       ${X[@]+"${X[@]}"} --mosaic-out $O/$arm.mosaic.tsv --gaf-base $D/short.gaf.db --gbz-base $D/sampled.gbz > $O/$arm.vcf 2> $O/$arm.err || { echo "$arm CALL_FAILED"; tail -3 $O/$arm.err; exit 1; }
    bgzip -f $O/$arm.vcf && tabix -f -p vcf $O/$arm.vcf.gz
  fi
  grep -E ' real|footprint' $O/$arm.err | tr -s ' ' | tr '\n' ' '; echo; grep -i "repeat" $O/$arm.err | grep -v "^\[vg call\] options" | head -4 | cut -c1-200
  if [ $arm = base ]; then
    cmp -s <(gzip -dc $O/base.vcf.gz | grep -v '^#') <(gzip -dc $D/short.vcf.gz | grep -v '^#') && echo "GATE: base records identical to the baseline call (vg 0575c3825)" || { echo "GATE_FAILED: base differs from the baseline"; exit 1; }
    cmp -s $O/base.mosaic.tsv $D/short.mosaic.tsv && echo "GATE: mosaic identical" || echo "GATE: mosaic differs"
  fi
  L=rs_$arm
  python3 tools/patch_contig.py score --contig chr20 --vcf $O/$arm.vcf.gz --label $L --out $O/score --patch $C/test/patch.json --threads 4 --no-haps > $O/score.$arm.log 2>&1 || echo "$arm SCORE_FAILED"
  ( cd $O/score && $E/work/truvari-venv/bin/truvari refine -u -a mafft -t 4 -r $L/score/chr20.truvari/candidate.refine.bed -f $L/chr20/chr20.fa $L/score/chr20.truvari > $L/refine.log 2>&1 ) || echo "$arm REFINE_FAILED"
done
step "panel walks and called haplotypes"
[ -s $O/panel_banded.tsv.gz ] || python3 $O/panel_banded.py > $O/panel_banded.log 2>&1 || echo PANEL_FAILED
python3 $O/called_arms.py base=$O/base.vcf.gz rep=$O/rep.vcf.gz desc=$O/desc.vcf.gz link=$O/link.vcf.gz > $O/called_arms.log 2>&1 || echo CALLED_FAILED
echo RS_DONE
