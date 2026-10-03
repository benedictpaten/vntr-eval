#!/bin/bash
# --repeat-sites with --dump-likelihoods (vg-diag: the dump also writes each allele's length and every diploid
# genotype's read and depth terms), filtered through a FIFO to the covered repeat sites. The VCF must equal rep's.
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
source work/stage4/chr20/env.sh
export TMPDIR=$PWD/work/tmp
D=work/stage4/chr20; P=$D/patched_st_medoid; VG=work/bin/vg-diag; O=work/iterate/wrongwalk
rm -f $O/dump.fifo; mkfifo $O/dump.fifo
awk -F'\t' 'NR==FNR {keep[$1]=$2; next}
  { s = ($1 ~ /^#(lengths|depth|gt)$/) ? $2 : $1; if (s ~ /^#/) next;
    split(s, p, "_"); k = p[1] "_" p[2]; gsub(/[+-]/, "", k);
    if (k in keep) print keep[k] "\t" $0 }' $O/dump_sites.tsv $O/dump.fifo > $O/dump.tsv &
echo "[$(date +%T)] call"
$VG call $P/graph.gbz -p CHM13#0#chr20 -d 2 -t 8 -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased \
    --repeat-sites $D/repeat_sites_st_medoid.bed --dump-likelihoods $O/dump.fifo --gaf-reads $P/reads.gaf.gz > $P/chr20.diag.vcf 2> $P/chr20.diag.err
echo "call exit $?"
wait
rm -f $O/dump.fifo
cmp <(bgzip -dc $P/chr20.rep.vcf.gz | grep -v '^##') <(grep -v '^##' $P/chr20.diag.vcf) && echo "DIAG VCF IDENTICAL TO rep (body)"
wc -l $O/dump.tsv
echo "[$(date +%T)] DUMP_DONE"
