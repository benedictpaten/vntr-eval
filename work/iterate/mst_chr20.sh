#!/bin/bash
# 4q: whole-chr20 patch with a nearest-neighbour threading arm (candidates in work/stage4/candidates/<v>__all/),
# mapped, called on vg 91d38c802, scored like patched_all, then truvari refine as refine4.sh.
# Usage: mst_chr20.sh mst [mst3 ...]
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
source work/stage4/chr20/env.sh
export VNTR_STAGE2_WORK=$PWD/work/iterate/stage2
export VNTR_STAGE3_VG=$PWD/work/bin/vg-91d38c802
export TMPDIR=$PWD/work/tmp; mkdir -p $TMPDIR
D=work/stage4/chr20
for v in "$@"; do
  P=$D/patched_$v
  echo "[$(date +%T)] $v patch"
  python3 tools/patch_contig.py patch --out $P --ids $D/patch_all_ids.txt --method ${v}__all > $D/patch_$v.log 2>&1 || { echo "PATCH_FAIL $v"; continue; }
  echo "[$(date +%T)] $v map"
  python3 tools/patch_contig.py map --gfa $P/graph.gfa --reads $D/reads --out $P --threads 8 > $D/map_$v.log 2>&1 || { echo "MAP_FAIL $v"; continue; }
  echo "[$(date +%T)] $v call"
  python3 tools/patch_contig.py call --gbz $P/graph.gbz --gaf $P/reads.gaf.gz --out $P/chr20.vcf.gz --threads 8 > $D/call_$v.log 2>&1 || { echo "CALL_FAIL $v"; continue; }
  echo "[$(date +%T)] $v score"
  python3 tools/patch_contig.py score --vcf $P/chr20.vcf.gz --label rb_patched_$v --out $D/score_rb --patch $D/patched_all/patch.json > $D/score_$v.log 2>&1 || { echo "SCORE_FAIL $v"; continue; }
  echo "[$(date +%T)] $v refine"
  ( cd $D && PATH=$HOME/.local/bin:/opt/homebrew/bin:$PATH /Users/benedictpaten/PycharmProjects/vg-call-eval/work/truvari-venv/bin/truvari refine -u -a mafft -t 6 \
      -r score_rb/rb_patched_$v/score/chr20.truvari/candidate.refine.bed -f score_rb/rb_patched_$v/chr20/chr20.fa \
      score_rb/rb_patched_$v/score/chr20.truvari > score_rb/rb_patched_$v/refine.log 2>&1 ) || echo "REFINE_FAIL $v"
  echo "[$(date +%T)] $v DONE"
done
echo "[$(date +%T)] MST_CHR20_DONE"
