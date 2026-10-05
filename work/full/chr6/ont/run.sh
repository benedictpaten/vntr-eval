#!/bin/bash
# chr6 held out, ONT: before (original hap32 chr6 graph) and latest (realigned, v2), both mapped and called with the latest vg.
# Steps as chr20's test_arm.sh, at 4 threads (chr20's ONT giraffe at 8 threads peaked at a 31 GB footprint).
set -u -o pipefail
R=/Users/benedictpaten/PycharmProjects/vntr-eval; E=/Users/benedictpaten/PycharmProjects/vg-call-eval; cd $R
export PATH=$HOME/.local/bin:/opt/homebrew/bin:$PATH TMPDIR=$R/work/tmp
VG=$R/work/bin/vg-nestfix4; O=work/full/chr6/ont; SC=work/stage4_chr6/chr6/score; T=4
step() { echo "[$(date +%T)] $*"; }
peak() { grep -E ' real|footprint' "$1" | tr -s ' ' | tr '\n' ' '; echo; }
step "ONT FASTQ from the genome-wide build's chr6 alignments"
[ -s $O/ont.chr6.fq.gz ] || { gzip -dc $E/work/chr6.ont.gaf.gz | $VG convert -F /dev/stdin $E/data/E821-16-sampled.gbz | $VG view -X - | gzip -1 > $O/ont.chr6.fq.gz.tmp && mv $O/ont.chr6.fq.gz.tmp $O/ont.chr6.fq.gz; } || { echo FASTQ_FAILED; exit 1; }
echo "reads: $(gzip -dc $O/ont.chr6.fq.gz | awk 'NR%4==2{n++; b+=length($0)} END{print n, b}')"
for arm in before latest; do
  if [ $arm = before ]; then G=work/stage4_chr6/chr6/unpatched/graph.gbz; DB=work/stage4_chr6/chr6/unpatched/graph.db
  else G=work/full/chr6/v2/patched/graph.gbz; DB=work/full/chr6/v2/patched/graph.db; fi
  D=$O/$arm; mkdir -p $D
  if [ ! -s $D/ont.gaf.gz ]; then
    step "$arm: lr-giraffe indexes"
    /usr/bin/time -l $VG autoindex -n -w lr-giraffe -G $G -p $D/map -t $T -T $TMPDIR -M 8G > $D/autoindex.log 2> $D/autoindex.err || { echo "$arm AUTOINDEX_FAILED"; exit 1; }
    peak $D/autoindex.err
    step "$arm: giraffe ONT (-b r10)"
    /usr/bin/time -l $VG giraffe -b r10 -Z $G -d $D/map.dist -m $D/map.longread.withzip.min -z $D/map.longread.zipcodes -o gaf -t $T \
       -f $O/ont.chr6.fq.gz 2> $D/giraffe.err | grep -v '^@' | bgzip -@ 2 > $D/ont.gaf.gz.tmp && mv $D/ont.gaf.gz.tmp $D/ont.gaf.gz || { echo "$arm GIRAFFE_FAILED"; exit 1; }
    peak $D/giraffe.err
    rm -f $D/map.dist $D/map.*.min $D/map.*.zipcodes
  fi
  gzip -dc $D/ont.gaf.gz | awk -F'\t' '{n++; if($6=="*") u++} END{print "ONT alignments", n, "unmapped", u+0}'
  step "$arm: GAF-Base"
  rm -f $D/fifo; mkfifo $D/fifo
  gaf-base sort $D/ont.gaf.gz --preset long -o $D/fifo -p > $D/sort.log 2>&1 &
  gaf-base construct $D/fifo -r $G -o $D/ont.gaf.db --overwrite > $D/gafdb.log 2>&1 || { echo "$arm GAFDB_FAILED"; exit 1; }
  wait; rm -f $D/fifo
  step "$arm: call --preset ont"
  /usr/bin/time -l $VG call $G -p CHM13#0#chr6 -d 2 -t $T -s HG002 --progress --preset ont --read-likelihood --phased \
     --gaf-base $D/ont.gaf.db --gbz-base $DB > $D/ont.vcf 2> $D/call.err || { echo "$arm CALL_FAILED"; exit 1; }
  peak $D/call.err
  bgzip -f $D/ont.vcf && tabix -f -p vcf $D/ont.vcf.gz
  step "$arm: score"
  L=ont_$arm
  python3 tools/patch_contig.py score --contig chr6 --vcf $D/ont.vcf.gz --label $L --out $SC --patch work/full/chr6/v2/patched/patch.json --threads 4 --no-haps > $D/score.log 2>&1 || { echo "$arm SCORE_FAILED"; exit 1; }
  ( cd $SC && $E/work/truvari-venv/bin/truvari refine -u -a mafft -t 4 -r $L/score/chr6.truvari/candidate.refine.bed -f $L/chr6/chr6.fa $L/score/chr6.truvari > $L/refine.log 2>&1 ) || echo "$arm REFINE_FAILED"
done
echo ONT_DONE
