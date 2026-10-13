#!/usr/bin/env python3
"""centrolign_paircost.py -- how far centrolign's pairwise alignments are from the optima, on qc's pairs.

    python3 tools/centrolign_paircost.py TR_chr20_1782084 TR_chr20_8073301 ...

For each hard-set region, aligns each of qc's 20 allele pairs (cached by `hardset_strategies.py optima`)
with centrolign, and prints the summed pair_cost over that of the optima. There are three settings:
the defaults (def), constraint_method 0 (c0), and c0 with its stitching scores set to 10x pgrealign's
(c0s).
"""
import gzip
import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hardset_strategies as H  # noqa: E402


def configs(wd):
    text = subprocess.run([H.CENTROLIGN, '-G'], capture_output=True, text=True).stdout
    c0 = text.replace('constraint_method: 3', 'constraint_method: 0').replace('logging_level: 2', 'logging_level: 0')
    c0s = c0.replace('stitch_mismatch: 80', 'stitch_mismatch: 40').replace(
        'stitch_gap_open: 60,800,2500', 'stitch_gap_open: 40,240,240').replace(
        'stitch_gap_extend: 30,5,1', 'stitch_gap_extend: 20,10,10')
    out = {'def': None}
    for k, t in (('c0', c0), ('c0s', c0s)):
        out[k] = os.path.join(wd, k + '.yaml')
        with open(out[k], 'w') as f:
            f.write(t)
    return out


def main(regions):
    P = H.paths(H.CONTIG)
    wd = tempfile.mkdtemp(prefix='clpc.', dir=os.path.join(H.config.WORK_DIR, 'tmp'))
    conf = configs(wd)
    for rid in regions:
        seqs = {a: s for a, s, _ in H.read_alleles(os.path.join(P['fasta'], rid + '.fa'))}
        cost = dict.fromkeys(['opt'] + list(conf), 0)
        with gzip.open(os.path.join(P['optima'], rid + '.jsonl.gz'), 'rt') as f:
            for line in f:
                r = json.loads(line)
                if r['kind'] != 'pair' or r.get('skipped'):
                    continue
                cost['opt'] += H.pair_cost(*r['rows'])
                for k, c in conf.items():
                    cost[k] += H.pair_cost(*H.centrolign_pair(seqs[r['a']], seqs[r['b']], wd, c))
        print(rid, ' '.join('%s %.3f' % (k, cost[k] / cost['opt']) for k in conf), flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
