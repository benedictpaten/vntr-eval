#!/bin/bash
# The adopted alignment (star2) with the banded fallback: star2's 14,492 ok regions plus the 19 regions it left
# out (18 too_big, 1 timeout) as realigned by pgrealign's banded fallback (investigate/coverage/excluded), tested
# through the same chr20 vg call arm as star2.
set -u -o pipefail
C=/Users/benedictpaten/PycharmProjects/vntr-eval/work/full/chr20; E=$C/investigate/coverage/excluded
while [ ! -s $E/msas.jsonl.gz ]; do kill -0 30468 2>/dev/null || { [ -s $E/msas.jsonl.gz ] || { echo MISSING; exit 1; }; }; sleep 60; done
python3 - <<PY
import gzip, json
new = {}
for l in gzip.open('$E/msas.jsonl.gz', 'rt'):
    m = json.loads(l); new[m['id']] = l
assert all(json.loads(l)['status'] == 'ok' for l in new.values()), 'an excluded region failed'
n = 0
with gzip.open('$C/msas.banded.jsonl.gz.tmp', 'wt') as out:
    for l in gzip.open('$C/msas.star2.jsonl.gz', 'rt'):
        m = json.loads(l)
        if m['id'] in new:
            assert m['status'] != 'ok', m['id']
            out.write(new.pop(m['id'])); n += 1
        else:
            out.write(l)
assert not new, sorted(new)
print('msas.banded: replaced', n)
PY
mv $C/msas.banded.jsonl.gz.tmp $C/msas.banded.jsonl.gz
bash $C/realign_arm.sh banded
