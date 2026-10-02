#!/bin/bash
# 4q: nearest-neighbour threading on every eligible chr20 repeat (work/stage4/chr20/patch_all_ids.txt); the hap32
# projections are copied to work/stage4/candidates/<v>__all/ for patch_contig.py.
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
source work/stage4/chr20/env.sh
export TMPDIR=$PWD/work/tmp; mkdir -p $TMPDIR
D=work/stage4/chr20
for v in "$@"; do
  echo "[$(date +%T)] align $v on $(wc -l < $D/patch_all_ids.txt) regions"
  python3 tools/iterate.py build $v --regions $(paste -sd, $D/patch_all_ids.txt) --jobs 3 > work/iterate/logs/$v.chr20_build.log 2>&1 || echo "ALIGN_FAIL $v"
  mkdir -p work/stage4/candidates/${v}__all
  for r in $(cat $D/patch_all_ids.txt); do
    for x in gfa msa.fa json; do
      [ -f work/iterate/candidates/$v/$r.$x ] && cp work/iterate/candidates/$v/$r.$x work/stage4/candidates/${v}__all/
    done
  done
  echo "[$(date +%T)] $v: $(ls work/stage4/candidates/${v}__all/*.gfa | wc -l) graphs"
done
echo "[$(date +%T)] ALIGN_DONE"
