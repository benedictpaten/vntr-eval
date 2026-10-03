#!/bin/bash
# Redo the realigned chr20 graph with dense node ids (replace --id-mode dense), verify it, then run its arm of the
# test (test_arm.sh) and print both arms' scores. The first build put 182,106 new nodes at ids >= 224,038,841, far
# from their regions; vg call's id-windowed read fetches then straddled windows (short call 54 min vs 6 min).
set -u -o pipefail
R=/Users/benedictpaten/PycharmProjects/vntr-eval
cd $R
export PYTHONPATH=$R PATH=$HOME/CLionProjects/vg/bin:$PATH TMPDIR=$R/work/tmp
C=work/full/chr20; PG="python3 -m pgrealign.cli"
step() { echo "[$(date +%T)] $*"; }
step "replace (dense ids)"
/usr/bin/time -l $PG replace --graph $C/chr20.gbz --extract $C/extract --msas $C/msas.jsonl.gz -o $C/chr20.realigned.gfa \
   --id-mode dense --gbz $C/chr20.realigned.gbz -t 8 > $C/replace.dense.json 2> $C/replace.dense.err || { echo REPLACE_FAILED; exit 1; }
cat $C/replace.dense.json; grep -E 'real|maximum' $C/replace.dense.err
step verify
$PG verify --before $C/chr20.gbz --after $C/chr20.realigned.gbz -o $C/verify.json > /dev/null 2> $C/verify.err || { echo VERIFY_FAILED; exit 1; }
python3 -c "import json; d=json.load(open('$C/verify.json')); print({k: d[k] for k in ('ok', 'paths', 'bases', 'metadata_identical')}); print(d['before']); print(d['after'])"
vg stats -r $C/chr20.realigned.gbz
step "test arm: realigned"
bash $C/test_arm.sh realigned $C/chr20.realigned.gbz
step summary
python3 - <<PY
import json, os
T = '$C/test/score'
for arm in ('original', 'realigned'):
    for rd in ('short', 'ont'):
        L = '%s_%s' % (arm, rd)
        p = '$C/test/%s/score.%s.log' % (arm, rd)
        log = open(p).read() if os.path.exists(p) else ''
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
step REDO_DONE
