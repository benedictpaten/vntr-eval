#!/bin/bash
# 4q: test-set calls of the nearest-neighbour threading variants (and st_chm13, which passes the site gate) on vg 91d38c802.
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
source work/stage4/chr20/env.sh
export VNTR_STAGE2_WORK=$PWD/work/iterate/stage2
export TMPDIR=$PWD/work/tmp; mkdir -p $TMPDIR
export VNTR_STAGE3_VG=$PWD/work/bin/vg-91d38c802
for v in mst mst3 st_chm13; do
  echo "[$(date +%T)] test set $v"
  python3 tools/iterate.py all $v --jobs 4 > work/iterate/logs/$v.all.log 2>&1 || echo "FAIL $v"
done
echo "[$(date +%T)] MST_TEST_DONE"
