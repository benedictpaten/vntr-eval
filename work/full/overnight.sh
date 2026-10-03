#!/bin/bash
# chr20 first (2026-10-02): the full ~456-haplotype chr20 graph realigned end to end and verified. chr6 (realign,
# replace, V1) waits until disk is freed.
set -u
R=/Users/benedictpaten/PycharmProjects/vntr-eval
cd $R
export PYTHONPATH=$R PATH=$HOME/CLionProjects/vg/bin:$PATH TMPDIR=$R/work/tmp
mkdir -p $TMPDIR
G=/Users/benedictpaten/PycharmProjects/vg-call-eval/data/hprc-v2.1-mc-chm13-eval.gref.gbz
PG="python3 -m pgrealign.cli"
step() { echo "[$(date +%T)] $*"; }
T() { /usr/bin/time -l "$@"; }

step "genome-wide node id range"
vg stats -r $G > work/full/eval.idrange.txt 2>/dev/null; cat work/full/eval.idrange.txt
MAXID=$(awk -F'[:\t]' '{print $NF}' work/full/eval.idrange.txt)   # 'node-id-range<TAB>1:224038840'
ID20=$((MAXID + 1)); ID6=$((MAXID + 50000001))
echo "new ids: chr20 from $ID20, chr6 from $ID6"

# ---------------------------------------------------------------- (1) chr20, full graph
C=work/full/chr20; mkdir -p $C
(
  step "chr20 prepare"
  T $PG prepare --gbz $G --contig chr20 -o $C -t 8 > $C/prepare.json 2> $C/prepare.err || { echo PREPARE_FAILED; exit 1; }
  step "chr20 regions"
  $PG regions --targets catalogs/chm13v2.0_tr.bed.gz --contig chr20 --snarls $C/chr20.snarls.pb --ref $C/chr20.ref.tsv.gz -o $C/chr20.regions.tsv.gz > $C/regions.json 2> $C/regions.err || { echo REGIONS_FAILED; exit 1; }
  cat $C/regions.json
  step "chr20 extract"
  T $PG extract --graph $C/chr20.gbz --regions $C/chr20.regions.tsv.gz -o $C/extract -t 8 > $C/extract.json 2> $C/extract.err || { echo EXTRACT_FAILED; exit 1; }
  cat $C/extract.json
  step "chr20 realign"
  T $PG realign --packages $C/extract/packages.jsonl.gz -o $C/msas.jsonl.gz --select any -j 8 > $C/realign.json 2> $C/realign.err || { echo REALIGN_FAILED; exit 1; }
  cat $C/realign.json
  step "chr20 replace"
  T $PG replace --graph $C/chr20.gbz --extract $C/extract --msas $C/msas.jsonl.gz -o $C/chr20.realigned.gfa \
     --id-start $ID20 --gbz $C/chr20.realigned.gbz -t 8 > $C/replace.json 2> $C/replace.err || { echo REPLACE_FAILED; exit 1; }
  cat $C/replace.json
  step "chr20 verify"
  T $PG verify --before $C/chr20.gbz --after $C/chr20.realigned.gbz -o $C/verify.json > /dev/null 2> $C/verify.err || echo VERIFY_FAILED
  python3 -c "import json; d=json.load(open('$C/verify.json')); print({k: d[k] for k in ('ok', 'paths', 'bases', 'metadata_identical', 'paths_differ')}); print(d['before']); print(d['after'])"
  step CHR20_DONE
)

step OVERNIGHT_DONE
