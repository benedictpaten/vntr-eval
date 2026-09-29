#!/bin/bash
# Stage 4d chr20: repeat-unit-aware aligner on hap32 (abPOA fallback) for every eligible TR region,
# realigned now; patch/map/call/score only after the full-panel arm (run_full.sh) has finished.
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
source work/stage4/chr20/env.sh
D=work/stage4/chr20
IDS=$(cat $D/full_ids.txt)
echo "[$(date +%T)] unit-aware realign $(echo $IDS | wc -w) regions"
python3 tools/realign_units.py $IDS --jobs 2 --threads 2 --fallback poa_abpoa --no-runtime > $D/units_realign.log 2>&1
for i in $IDS; do [ -s work/stage4/candidates/unit_aware/$i.gfa ] && echo $i; done > $D/patch_units_ids.txt
echo "[$(date +%T)] $(wc -l < $D/patch_units_ids.txt) candidates; waiting for the full-panel arm"
until grep -q -E "RUN_FULL_DONE|FAILED" $D/run_full.log; do sleep 60; done
echo "[$(date +%T)] patch"
python3 tools/patch_contig.py patch --out $D/patched_units --ids $D/patch_units_ids.txt --method unit_aware > $D/patch_units.log 2>&1 || { echo PATCH_FAILED; exit 1; }
echo "[$(date +%T)] map"
python3 tools/patch_contig.py map --gfa $D/patched_units/graph.gfa --reads $D/reads --out $D/patched_units --threads 8 > $D/map_units.log 2>&1 || { echo MAP_FAILED; exit 1; }
echo "[$(date +%T)] call"
python3 tools/patch_contig.py call --gbz $D/patched_units/graph.gbz --gaf $D/patched_units/reads.gaf.gz --out $D/patched_units/chr20.vcf.gz --threads 8 > $D/call_units.log 2>&1 || { echo CALL_FAILED; exit 1; }
echo "[$(date +%T)] score"
for x in "patched_units:$D/patched_units/chr20.vcf.gz" "patched_all:$D/patched_all/chr20.vcf.gz" "unpatched:$D/unpatched/chr20.vcf.gz"; do
  python3 tools/patch_contig.py score --vcf ${x#*:} --label ${x%%:*} --out $D/score_units --patch $D/patched_units/patch.json --threads 2 > $D/score_units_${x%%:*}.log 2>&1 &
done
wait
echo "[$(date +%T)] RUN_UNITS_DONE"
