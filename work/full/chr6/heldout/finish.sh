#!/bin/bash
# chr6 held out: finish the revised realign, then qc it against the original graph, under a 12 GB budget.
set -u -o pipefail
cd /Users/benedictpaten/PycharmProjects/vntr-eval/work/full/chr6/heldout
export PYTHONPATH=/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-2026-07-vg-call-refactor/e061d8cf-7996-49a4-986b-3686374ef250/scratchpad/pgrealign-run PATH=$HOME/.local/bin:/opt/homebrew/bin:$PATH TMPDIR=$HOME/PycharmProjects/vntr-eval/work/tmp
echo "[$(date +%T)] realign"
[ -s msas.new.jsonl.gz ] || /usr/bin/time -l python3 -m pgrealign.cli realign --packages ../extract.v2/packages.jsonl.gz -o msas.new.jsonl.gz -j 4 --max-memory 12 > realign.new.json 2> realign.new.err || { echo REALIGN_FAILED; tail -5 realign.new.err; exit 1; }
cat realign.new.json; grep -E ' real|footprint' realign.new.err
echo "[$(date +%T)] qc"
/usr/bin/time -l python3 -m pgrealign.cli qc --extract ../extract.v2 --msas msas.new.jsonl.gz --graph ../chr6.nogref.gbz --vg $HOME/CLionProjects/vg/bin/vg -j 4 --max-memory 12 -o qc.new.tsv > qc.new.json 2> qc.new.err || { echo QC_FAILED; tail -5 qc.new.err; exit 1; }
grep -E ' real|footprint' qc.new.err
python3 ../../qc_compare.py chr6=qc.new.tsv
echo FINISH_DONE
