#!/bin/bash
# chr20 test of the realigned full graph (2026-10-02): for the original and the realigned ~456-haplotype chr20
# graphs alike, sample 32 HG002 haplotypes (as the production hap32 graph was sampled), map HG002 short reads and
# HG002 ONT reads (only those the genome-wide builds placed on chr20), call, and score against T2T-Q100 v1.1.
#
#   shared:  HG002 29-mer counts from the chr20 short reads (jellyfish -> tools/c/jf2kff, no KMC; checked to
#            sample identically to KMC's KFF on vg's own test data); ONT FASTQ rebuilt from the genome-wide
#            build's chr20 alignments (vg convert -F with the E821 chr20 node sequences)
#   per arm: vg autoindex -w sampling -> vg haplotypes --num-haplotypes 32 --include-reference (no diploid
#            sampling: production hap32 has 32 recombinants + CHM13 + GRCh38) -> CHM13 path == chr20.fa check ->
#            vg autoindex -w sr-giraffe -w lr-giraffe on the sampled GBZ -> giraffe (short: paired + single,
#            fragment 401.7 +- 166.4 as stage 4; ONT: -b r10) -> GAF-Base -> vg call (production flags; ONT
#            --preset ont) -> score (bench_wgs via patch_contig.py score) + truvari refine
# One binary throughout: vg 91d38c802 (positional depth-rate window, so the realigned graph's new node ids do
# not move the depth term).
set -u -o pipefail
R=/Users/benedictpaten/PycharmProjects/vntr-eval
E=/Users/benedictpaten/PycharmProjects/vg-call-eval
cd $R
VG=$R/work/bin/vg-91d38c802
export PATH=$HOME/.local/bin:/opt/homebrew/bin:$PATH TMPDIR=$R/work/tmp
C=work/full/chr20; T=$C/test; mkdir -p $T $TMPDIR
THREADS=8
step() { echo "[$(date +%T)] $*"; }
TM() { /usr/bin/time -l "$@"; }
peak() { grep -E 'real|maximum resident' "$1" | tr -s ' ' | tr '\n' ' '; echo; }

# One graph ("arm") of the chr20 test: test_arm.sh NAME GRAPH. Needs $T/HG002.chr20.kff, $T/ont.chr20.fq.gz and
# $T/patch.json from test_chain.sh.
arm() {
  A=$1; G=$2; D=$T/$A; mkdir -p $D
  step "$A: sampling index (vg autoindex -w sampling)"
  [ -s $D/full.hapl ] || TM $VG autoindex -n -w sampling -G $G -p $D/full -t $THREADS -T $TMPDIR -M 16G > $D/sampling.log 2> $D/sampling.err || { echo "$A SAMPLING_INDEX_FAILED"; return 1; }
  peak $D/sampling.err
  step "$A: sample 32 HG002 haplotypes"
  [ -s $D/sampled.gbz ] || TM $VG haplotypes -v 2 -t $THREADS -i $D/full.hapl -k $T/HG002.chr20.kff --include-reference \
      --set-reference CHM13 --set-reference GRCh38 --num-haplotypes 32 -g $D/sampled.gbz $G > $D/haplotypes.log 2> $D/haplotypes.err || { echo "$A SAMPLING_FAILED"; return 1; }
  peak $D/haplotypes.err; grep -i -E 'coverage|haplotypes' $D/haplotypes.err | head -5
  $VG paths -M -x $D/sampled.gbz | awk 'NR>1{print $2, $3}' | sort | uniq -c
  $VG paths -F -x $D/sampled.gbz -Q CHM13#0#chr20 | sed '1s/.*/>chr20/' > $D/chm13.fa
  if /usr/bin/cmp -s <(grep -v '>' $D/chm13.fa | tr -d '\n') <(grep -v '>' $E/work/wgs/chr20/chr20.fa | tr -d '\n'); then echo "$A: CHM13#0#chr20 identical to chr20.fa"; else echo "$A: CHM13 PATH DIFFERS"; return 1; fi
  rm -f $D/chm13.fa
  step "$A: giraffe indexes"
  [ -s $D/map.dist ] || TM $VG autoindex -n -w sr-giraffe -w lr-giraffe -G $D/sampled.gbz -p $D/map -t $THREADS -T $TMPDIR -M 12G > $D/autoindex.log 2> $D/autoindex.err || { echo "$A AUTOINDEX_FAILED"; return 1; }
  peak $D/autoindex.err
  # ---- short reads
  step "$A: giraffe short reads"
  if [ ! -s $D/short.gaf.gz ]; then
    TM $VG giraffe -Z $D/sampled.gbz -d $D/map.dist -m $D/map.shortread.withzip.min -z $D/map.shortread.zipcodes -o gaf -t $THREADS \
       -f work/stage4/chr20/reads/reads_1.fq.gz -f work/stage4/chr20/reads/reads_2.fq.gz --fragment-mean 401.7 --fragment-stdev 166.4 2> $D/giraffe.paired.err | bgzip -@ 2 > $D/short.paired.gaf.gz
    peak $D/giraffe.paired.err
    TM $VG giraffe -Z $D/sampled.gbz -d $D/map.dist -m $D/map.shortread.withzip.min -z $D/map.shortread.zipcodes -o gaf -t $THREADS \
       -f work/stage4/chr20/reads/reads_se.fq.gz 2> $D/giraffe.single.err | bgzip -@ 2 > $D/short.single.gaf.gz
    (gzip -dc $D/short.paired.gaf.gz; gzip -dc $D/short.single.gaf.gz) | grep -v '^@' | bgzip -@ 4 > $D/short.gaf.gz && rm -f $D/short.paired.gaf.gz $D/short.single.gaf.gz
  fi
  gzip -dc $D/short.gaf.gz | awk -F'\t' '{n++; if($6=="*") u++} END{print "short alignments", n, "unmapped", u+0}'
  step "$A: short GAF-Base and call"
  rm -f $D/fifo; mkfifo $D/fifo
  gaf-base sort $D/short.gaf.gz -o $D/fifo -p > $D/short.sort.log 2>&1 &
  TM gaf-base construct $D/fifo -r $D/sampled.gbz -o $D/short.gaf.db --overwrite > $D/short.gafdb.log 2> $D/short.gafdb.err || { echo "$A SHORT_GAFDB_FAILED"; return 1; }
  wait; rm -f $D/fifo
  TM $VG call $D/sampled.gbz -p CHM13#0#chr20 -d 2 -t $THREADS -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased \
     --mosaic-out $D/short.mosaic.tsv --gaf-base $D/short.gaf.db --gbz-base $D/sampled.gbz > $D/short.vcf 2> $D/short.call.err || { echo "$A SHORT_CALL_FAILED"; return 1; }
  peak $D/short.call.err
  bgzip -f $D/short.vcf && tabix -f -p vcf $D/short.vcf.gz
  # ---- ONT
  step "$A: giraffe ONT (-b r10)"
  if [ ! -s $D/ont.gaf.gz ]; then
    TM $VG giraffe -b r10 -Z $D/sampled.gbz -d $D/map.dist -m $D/map.longread.withzip.min -z $D/map.longread.zipcodes -o gaf -t $THREADS \
       -f $T/ont.chr20.fq.gz 2> $D/giraffe.ont.err | grep -v '^@' | bgzip -@ 2 > $D/ont.gaf.gz
    peak $D/giraffe.ont.err
  fi
  gzip -dc $D/ont.gaf.gz | awk -F'\t' '{n++; if($6=="*") u++} END{print "ONT alignments", n, "unmapped", u+0}'
  step "$A: ONT GAF-Base and call"
  rm -f $D/fifo; mkfifo $D/fifo
  gaf-base sort $D/ont.gaf.gz --preset long -o $D/fifo -p > $D/ont.sort.log 2>&1 &
  TM gaf-base construct $D/fifo -r $D/sampled.gbz -o $D/ont.gaf.db --overwrite > $D/ont.gafdb.log 2> $D/ont.gafdb.err || { echo "$A ONT_GAFDB_FAILED"; return 1; }
  wait; rm -f $D/fifo
  TM $VG call $D/sampled.gbz -p CHM13#0#chr20 -d 2 -t $THREADS -s HG002 --progress --preset ont --read-likelihood --phased \
     --mosaic-out $D/ont.mosaic.tsv --gaf-base $D/ont.gaf.db --gbz-base $D/sampled.gbz > $D/ont.vcf 2> $D/ont.call.err || { echo "$A ONT_CALL_FAILED"; return 1; }
  peak $D/ont.call.err
  bgzip -f $D/ont.vcf && tabix -f -p vcf $D/ont.vcf.gz
  # ---- score
  step "$A: score"
  for rd in short ont; do
    L=${A}_${rd}
    python3 tools/patch_contig.py score --contig chr20 --vcf $D/$rd.vcf.gz --label $L --out $T/score --patch $T/patch.json --threads 4 --no-haps > $D/score.$rd.log 2>&1 || echo "$A $rd SCORE_FAILED"
    ( cd $T/score && $E/work/truvari-venv/bin/truvari refine -u -a mafft -t 4 -r $L/score/chr20.truvari/candidate.refine.bed \
        -f $L/chr20/chr20.fa $L/score/chr20.truvari > $L/refine.log 2>&1 ) || echo "$A $rd REFINE_FAILED"
  done
  # the giraffe and sampling indexes are rebuilt in minutes: drop them to keep disk
  rm -f $D/map.dist $D/map.*.min $D/map.*.zipcodes
  step "$A done"
}


[ -s $T/ont.chr20.fq.gz ] && [ -s $T/HG002.chr20.kff ] && [ -s $T/patch.json ] || { echo "shared inputs missing"; exit 1; }
arm "$1" "$2"
