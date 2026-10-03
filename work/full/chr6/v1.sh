#!/bin/bash
# V1, held-out chr6: the medoid-star protocol on the full panel through pgrealign, applied to the chr6 hap32
# graph, re-mapped, called and scored against the earlier chr6 arms (stage4_chr6: unpatched, abPOA-hap32).
#   union of the full graph's and the hap32 graph's packages -> realign (medoid star, select any) -> project onto
#   hap32 -> replace in the hap32 chr6 GBZ -> map (giraffe, as stage4_chr6) -> GAF-Base -> vg call -> score.
set -u
R=/Users/benedictpaten/PycharmProjects/vntr-eval
cd $R
export PYTHONPATH=$R PATH=$HOME/CLionProjects/vg/bin:$PATH
source work/stage4_chr6/env.sh
export VNTR_STAGE2_WORK=$R/work/iterate/stage2 TMPDIR=$R/work/tmp   # env.sh points both into a cleared scratchpad
F=work/full/chr6; V=$F/v1; S=work/stage4_chr6; D=$S/chr6; E=/Users/benedictpaten/PycharmProjects/vg-call-eval
VG=work/bin/vg-2a6a228a5
mkdir -p $V
step() { echo "[$(date +%T)] $*"; }
until [ -s $F/extract.v2/summary.json ]; do sleep 30; done
step union
python3 -m pgrealign.cli union --base $F/extract.v2/packages.jsonl.gz --add $F/hap32/extract/packages.jsonl.gz -o $V/union.jsonl.gz > $V/union.json || exit 1
step realign
/usr/bin/time -l python3 -m pgrealign.cli realign --packages $V/union.jsonl.gz -o $V/msas.union.jsonl.gz --select any -j 8 > $V/realign.json 2> $V/realign.err || { echo REALIGN_FAILED; exit 1; }
cat $V/realign.json
step project
python3 -m pgrealign.cli project --union $V/union.jsonl.gz --msas $V/msas.union.jsonl.gz -o $V/msas.hap32.jsonl.gz > $V/project.json || exit 1
step replace
mkdir -p $V/patched
/usr/bin/time -l python3 -m pgrealign.cli replace --graph $E/work/wgs/chr6/chr6.gbz --extract $F/hap32/extract --msas $V/msas.hap32.jsonl.gz \
  -o $V/patched/graph.gfa --id-start 300000000 --gbz $V/patched/graph.gbz --keep-gfa > $V/replace.json 2> $V/replace.err || { echo REPLACE_FAILED; exit 1; }
cat $V/replace.json
python3 - <<PY
import gzip, json
regs = {}
for l in gzip.open('$V/msas.hap32.jsonl.gz', 'rt'):
    m = json.loads(l); regs[m['id']] = None
spans = {}
for l in gzip.open('$F/chr6.regions.tsv.gz', 'rt'):
    if l.startswith('#'): continue
    x = l.rstrip('\n').split('\t')
    if x[0] in regs: spans[x[0]] = {'span': [int(x[4]) + 1, int(x[5])]}
json.dump({'method': 'pgrealign medoid star (full panel, projected)', 'regions': spans}, open('$V/patched/patch.json', 'w'))
print('patch.json regions', len(spans))
PY
step map
python3 tools/patch_contig.py map --contig chr6 --gfa $V/patched/graph.gfa --reads $D/reads --out $V/patched --threads 8 > $V/map.log 2>&1 || { echo MAP_FAILED; exit 1; }
step dbs
A=$V/patched
rm -f $A/fifo; mkfifo $A/fifo
gbz-base construct $A/graph.gbz -o $A/graph.db --overwrite > $V/gbzdb.log 2>&1
gaf-base sort $A/reads.gaf.gz -o $A/fifo -p > $V/sort.log 2>&1 &
gaf-base construct $A/fifo -r $A/graph.gbz -o $A/reads.db --overwrite > $V/gafdb.log 2>&1 || { echo GAFDB_FAILED; exit 1; }
wait; rm -f $A/fifo
step call
$VG call $A/graph.gbz -p CHM13#0#chr6 -d 2 -t 8 -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased \
  --gaf-base $A/reads.db --gbz-base $A/graph.db > $A/chr6.vcf 2> $V/call.log || { echo CALL_FAILED; exit 1; }
bgzip -f $A/chr6.vcf && tabix -f -p vcf $A/chr6.vcf.gz
step score
python3 tools/patch_contig.py score --contig chr6 --vcf $A/chr6.vcf.gz --label pgrealign_v1 --out $D/score --patch $A/patch.json --threads 4 > $V/score.log 2>&1
for x in patched unpatched; do
  python3 tools/patch_contig.py score --contig chr6 --vcf $D/$x/chr6.vcf.gz --label ${x}_vs_v1 --out $D/score --patch $A/patch.json --threads 2 > $V/score_$x.log 2>&1 &
done
wait
step V1_DONE
