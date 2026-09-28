#!/bin/bash
# Stage 4b chr20: realign EVERY eligible (packaged) TR region with abPOA, patch, re-map, call, score.
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
source work/stage4/chr20/env.sh
D=work/stage4/chr20
IDS=$(for d in work/stage4/regions/*/; do [ -s "$d/hap32.fa" ] && basename "$d"; done)
echo "[$(date +%T)] realign $(echo $IDS | wc -w) regions"
python3 tools/realign.py poa_abpoa $IDS --jobs 3 --threads 1 > $D/realign_all.log 2>&1
for i in $IDS; do [ -s work/stage4/candidates/poa_abpoa/$i.gfa ] && echo $i; done > $D/patch_all_ids.txt
echo "[$(date +%T)] patch $(wc -l < $D/patch_all_ids.txt) candidates"
python3 tools/patch_contig.py patch --out $D/patched_all --ids $D/patch_all_ids.txt --method poa_abpoa > $D/patch_all.log 2>&1 || { echo PATCH_FAILED; exit 1; }
echo "[$(date +%T)] map"
python3 tools/patch_contig.py map --gfa $D/patched_all/graph.gfa --reads $D/reads --out $D/patched_all --threads 8 > $D/map_all.log 2>&1 || { echo MAP_FAILED; exit 1; }
echo "[$(date +%T)] call"
python3 tools/patch_contig.py call --gbz $D/patched_all/graph.gbz --gaf $D/patched_all/reads.gaf.gz --out $D/patched_all/chr20.vcf.gz --threads 8 > $D/call_all.log 2>&1 || { echo CALL_FAILED; exit 1; }
echo "[$(date +%T)] score"
for x in "patched_all:$D/patched_all/chr20.vcf.gz" "unpatched:$D/unpatched/chr20.vcf.gz" "production:/Users/benedictpaten/PycharmProjects/vg-call-eval/work/wgs-mm095/chr20/chr20.vcf.gz"; do
  python3 tools/patch_contig.py score --vcf ${x#*:} --label ${x%%:*} --out $D/score_all --patch $D/patched_all/patch.json --threads 2 > $D/score_all_${x%%:*}.log 2>&1 &
done
wait
echo "[$(date +%T)] RUN_ALL_DONE"
