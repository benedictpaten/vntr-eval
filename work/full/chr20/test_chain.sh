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

# ------------------------------------------------------------------ shared inputs (start at once)
if [ ! -s $T/HG002.chr20.kff ]; then
  step "k-mer counts (jellyfish, chr20 short reads)"
  TM jellyfish count -m 29 -C -L 2 -s 600M -t 6 -o $TMPDIR/HG002.chr20.jf \
     <(gzip -dc work/stage4/chr20/reads/reads_1.fq.gz work/stage4/chr20/reads/reads_2.fq.gz work/stage4/chr20/reads/reads_se.fq.gz) 2> $T/jellyfish.err
  peak $T/jellyfish.err
  jellyfish dump -c -L 2 $TMPDIR/HG002.chr20.jf | work/bin/jf2kff $T/HG002.chr20.kff.tmp 29 255 2> $T/jf2kff.err && mv $T/HG002.chr20.kff.tmp $T/HG002.chr20.kff
  cat $T/jf2kff.err; work/bin/jf2kff --check $T/HG002.chr20.kff
  rm -f $TMPDIR/HG002.chr20.jf     # the converter's own intermediate
fi
if [ ! -s $T/ont.chr20.fq.gz ]; then
  step "ONT FASTQ (genome-wide build's chr20 reads, rebuilt from their E821 alignments)"
  gzip -dc $E/work/chr20.ont.gaf.gz | $VG convert -F /dev/stdin $E/data/E821-16-sampled.gbz | $VG view -X - | gzip -1 > $T/ont.chr20.fq.gz.tmp \
    && mv $T/ont.chr20.fq.gz.tmp $T/ont.chr20.fq.gz
  gzip -dc $T/ont.chr20.fq.gz | awk 'NR%4==2{n++; b+=length($0)} END{print "ONT reads", n, "bases", b}'
fi

# ------------------------------------------------------------------ wait for the graph build
until grep -q OVERNIGHT_DONE work/full/overnight.log 2>/dev/null; do sleep 60; done
python3 -c "import json,sys; d=json.load(open('$C/verify.json')); sys.exit(0 if d['ok'] else 1)" || { echo "VERIFY NOT OK: not testing the realigned graph"; exit 1; }

# inside/outside spans for the scorer: the realigned regions
python3 - <<PY
import gzip, json
ok = {json.loads(l)['id'] for l in gzip.open('$C/msas.jsonl.gz', 'rt') if json.loads(l).get('status') == 'ok'}
spans = {}
for l in gzip.open('$C/chr20.regions.tsv.gz', 'rt'):
    if l.startswith('#'): continue
    x = l.rstrip('\n').split('\t')
    if x[0] in ok: spans[x[0]] = {'span': [int(x[4]) + 1, int(x[5])]}
json.dump({'method': 'pgrealign medoid star, full graph', 'regions': spans}, open('$T/patch.json', 'w'))
print('realigned regions for inside/outside counts:', len(spans))
PY

bash $C/test_arm.sh original $C/chr20.gbz
bash $C/test_arm.sh realigned $C/chr20.realigned.gbz
step "summary"
python3 - <<PY
import json, os
T = '$T/score'
for arm in ('original', 'realigned'):
    for rd in ('short', 'ont'):
        L = '%s_%s' % (arm, rd)
        log = open('$T/%s/score.%s.log' % (arm, rd)).read() if os.path.exists('$T/%s/score.%s.log' % (arm, rd)) else ''
        try:
            s = json.loads(log[log.index('{'):log.rindex('}') + 1])
        except ValueError:
            print(L, 'no score'); continue
        rf = '%s/%s/score/chr20.truvari/refine.variant_summary.json' % (T, L)
        r = json.load(open(rf)) if os.path.exists(rf) else {}
        print('%-18s ALL %.4f SNV %.4f indel %.4f | SV raw %.4f (FP %d FN %d) refined %s (FP %s FN %s) | SV FP in/out %s/%s FN in/out %s/%s' % (
            L, s['all_f1'], s['snv_f1'], s['indel_f1'], s['sv_f1'], s['sv_fp'], s['sv_fn'],
            ('%.4f' % r['f1']) if 'f1' in r else '?', r.get('FP', '?'), r.get('FN', '?'),
            s.get('sv_fp_in_patched'), s.get('sv_fp_outside'), s.get('sv_fn_in_patched'), s.get('sv_fn_outside')))
PY
step TEST_CHAIN_DONE
