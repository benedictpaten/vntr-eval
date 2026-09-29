#!/bin/bash
# Stage 4c chr20: full HG002-free panel per eligible TR region -> abPOA -> project to hap32 -> patch -> re-map -> call -> score.
set -u
cd /Users/benedictpaten/PycharmProjects/vntr-eval
source work/stage4/chr20/env.sh
D=work/stage4/chr20
IDS=$(for d in work/stage4/regions/*/; do [ -s "$d/hap32.fa" ] && basename "$d"; done)
echo "$IDS" > $D/full_ids.txt
awk "NR%2==1" $D/full_ids.txt > $D/full_half.aa; awk "NR%2==0" $D/full_ids.txt > $D/full_half.ab
echo "[$(date +%T)] hprc-fetch $(wc -l < $D/full_ids.txt) regions, 2 fetchers"
for h in $D/full_half.aa $D/full_half.ab; do
  python3 tools/package_regions.py hprc-fetch --loci work/stage4/prefilter.bed --work work/stage4/pkgwork --only "$(paste -sd, $h)" --keep-order > $h.log 2>&1 &
done
wait
echo "[$(date +%T)] hprc files"
python3 tools/package_regions.py hprc --out work/stage4/regions --work work/stage4/pkgwork --jobs 4 > $D/hprc_files.log 2>&1
echo "[$(date +%T)] full-panel abPOA ($(ls work/stage4/regions/*/hprc.fa.gz 2>/dev/null | wc -l) regions with a panel)"
python3 tools/poa_panel.py run --methods poa_abpoa --regions "$(paste -sd, $D/full_ids.txt)" --jobs 2 --mem-mb 12000 --budget-mb 20000 \
  --cand-root work/stage4/candidates --panel-root work/stage4/panel --no-runtime > $D/full_align.log 2>&1
for i in $IDS; do [ -s work/stage4/candidates/poa_abpoa__all/$i.gfa ] && echo $i; done > $D/patch_full_ids.txt
echo "[$(date +%T)] patch $(wc -l < $D/patch_full_ids.txt) candidates"
python3 tools/patch_contig.py patch --out $D/patched_full --ids $D/patch_full_ids.txt --method poa_abpoa__all > $D/patch_full.log 2>&1 || { echo PATCH_FAILED; exit 1; }
echo "[$(date +%T)] map"
python3 tools/patch_contig.py map --gfa $D/patched_full/graph.gfa --reads $D/reads --out $D/patched_full --threads 8 > $D/map_full.log 2>&1 || { echo MAP_FAILED; exit 1; }
echo "[$(date +%T)] call"
python3 tools/patch_contig.py call --gbz $D/patched_full/graph.gbz --gaf $D/patched_full/reads.gaf.gz --out $D/patched_full/chr20.vcf.gz --threads 8 > $D/call_full.log 2>&1 || { echo CALL_FAILED; exit 1; }
echo "[$(date +%T)] score"
for x in "patched_full:$D/patched_full/chr20.vcf.gz" "patched_all:$D/patched_all/chr20.vcf.gz" "unpatched:$D/unpatched/chr20.vcf.gz"; do
  python3 tools/patch_contig.py score --vcf ${x#*:} --label ${x%%:*} --out $D/score_full --patch $D/patched_full/patch.json --threads 2 > $D/score_full_${x%%:*}.log 2>&1 &
done
wait
echo "[$(date +%T)] RUN_FULL_DONE"
