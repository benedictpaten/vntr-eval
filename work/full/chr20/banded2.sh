#!/bin/bash
# chr20 with the 19 regions only the banded fallback realigns patched in (msas.banded = star2 + those 19),
# through the chr20 test chain, calls with the latest vg. One step at a time; giraffe at 4 threads.
set -u -o pipefail
R=/Users/benedictpaten/PycharmProjects/vntr-eval; cd $R
C=work/full/chr20; VG=$R/work/bin/vg-f7e130b16; N=banded
export PYTHONPATH=/Users/benedictpaten/PycharmProjects/pgrealign PATH=$HOME/.local/bin:/opt/homebrew/bin:$PATH TMPDIR=$R/work/tmp
PG="python3 -P -m pgrealign.cli"
step() { echo "[$(date +%T)] $*"; }
step "replace (dense ids)"
[ -s $C/chr20.$N.gbz ] || /usr/bin/time -l $PG replace --graph $C/chr20.gbz --extract $C/extract --msas $C/msas.$N.jsonl.gz \
   -o $C/chr20.$N.gfa --id-mode dense --gbz $C/chr20.$N.gbz -t 4 --vg $VG > $C/replace.$N.json 2> $C/replace.$N.err || { echo REPLACE_FAILED; tail -3 $C/replace.$N.err; exit 1; }
grep -E ' real|footprint' $C/replace.$N.err | tr -s ' ' | tr '\n' ' '; echo; df -h $R | tail -1
step verify
$PG verify --before $C/chr20.gbz --after $C/chr20.$N.gbz --vg $VG -o $C/verify.$N.json > /dev/null 2> $C/verify.$N.err || { echo VERIFY_FAILED; exit 1; }
python3 -c "import json; d=json.load(open('$C/verify.$N.json')); print({k: d[k] for k in ('ok', 'paths', 'bases', 'metadata_identical')}); assert d['ok']" || { echo VERIFY_FAILED; exit 1; }
step "test arm"
bash $C/test_arm_latest.sh $N $C/chr20.$N.gbz || { echo ARM_FAILED; exit 1; }
echo BANDED_DONE
