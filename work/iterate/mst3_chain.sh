#!/bin/bash
# 4q: once mst3's chr20 alignment is copied into work/stage4/candidates/mst3__all, give every region where mst3
# did not finish under the 900 s cap the mst graph (listed in mst3__all/fallback_mst.txt), then run the
# whole-chr20 pipeline.
cd /Users/benedictpaten/PycharmProjects/vntr-eval
until grep -q "\] mst3: .* graphs" work/iterate/mst_chr20_align.log; do sleep 30; done
C=work/stage4/candidates
python3 -c "
import json,glob
print('\n'.join(json.load(open(p))['region_id'] for p in glob.glob('work/iterate/candidates/mst3/TR*.json') if json.load(open(p))['status']!='ok'))" > $C/mst3__all/fallback_mst.txt
for r in $(cat $C/mst3__all/fallback_mst.txt); do
  grep -qx $r work/stage4/chr20/patch_all_ids.txt || continue
  for x in gfa msa.fa; do cp $C/mst__all/$r.$x $C/mst3__all/; done
done
echo "[$(date +%T)] mst3: $(ls $C/mst3__all/*.gfa | wc -l) graphs, $(wc -l < $C/mst3__all/fallback_mst.txt) from mst"
work/iterate/mst_chr20.sh mst3
echo "[$(date +%T)] CHAIN_DONE"
