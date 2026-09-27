#!/usr/bin/env python3
"""Tests of realign_poa.py: input orders, and a synthetic end-to-end run of both tools through
realign.py's pipeline (dedup, N masking, row checks) and msa_graph.py.

    python3 tools/test_realign_poa.py
"""
import os
import random
import shutil
import sys
import tempfile

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import msa_graph  # noqa: E402
import realign  # noqa: E402
import realign_poa as rp  # noqa: E402


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print('ok  ', msg)


def synthetic(seed=7):
    """A VNTR-like set: shared flanks, arrays of 8-20 copies of a 21-bp unit with point
    variants, one duplicate, one N run."""
    rnd = random.Random(seed)
    left = ''.join(rnd.choice('ACGT') for _ in range(120))
    right = ''.join(rnd.choice('ACGT') for _ in range(120))
    unit = 'GGAGGGTGAGGTGTGTACCAT'
    recs = []
    for i in range(9):
        arr = []
        for _ in range(rnd.randint(8, 20)):
            u = list(unit)
            if rnd.random() < 0.3:
                u[rnd.randrange(len(u))] = rnd.choice('ACGT')
            arr.append(''.join(u))
        recs.append(('hap#%d' % i, left + ''.join(arr) + right))
    recs.insert(0, ('CHM13#0#chrT', recs[0][1][:-5] + 'ACGTA'))
    recs.append(('dup#1', recs[3][1]))
    s = recs[5][1]
    recs[5] = (recs[5][0], s[:150] + 'N' * 30 + s[150:])
    return recs


def main():
    seqs = ['ACGT' * 10, 'ACGT' * 30, 'ACGA' * 20, 'ACGT' * 11]
    check(rp.order_indices(seqs, 'given') == [0, 1, 2, 3], 'given order')
    check(rp.order_indices(seqs, 'longest') == [1, 2, 3, 0], 'longest-first order (stable)')
    check(rp.order_indices(seqs, 'shortest') == [0, 3, 2, 1], 'shortest-first order')
    r1, r2 = rp.order_indices(seqs, 'random:1'), rp.order_indices(seqs, 'random:1')
    check(r1 == r2 and sorted(r1) == [0, 1, 2, 3], 'random:SEED is a reproducible permutation')
    g = rp.guide_order(seqs)
    check(sorted(g) == [0, 1, 2, 3], 'guide order is a permutation')
    # the two nearly identical sequences (0 and 3) must be adjacent in Prim's order
    check(abs(g.index(0) - g.index(3)) == 1 or g[0] in (0, 3), 'guide order keeps near-identical sequences together')
    a, b = rp._kmer_profile('ACGT' * 10), rp._kmer_profile('ACGT' * 20)
    check(0 < rp.kmer_distance(a, b) < 1, 'k-mer distance counts copy number')
    check(rp.kmer_distance(a, a) == 0.0, 'k-mer distance of a sequence to itself is 0')
    check(rp.predict_mb('spoa', [10000]) > rp.predict_mb('abpoa', [10000], ['-p']), 'spoa predicted above abPOA')

    tmp = tempfile.mkdtemp(prefix='test_realign_poa.')
    try:
        recs = synthetic()
        fa = os.path.join(tmp, 'in.fa')
        with open(fa, 'w') as f:
            for n, s in recs:
                f.write('>%s desc\n%s\n' % (n, s))
        for method in ('poa_abpoa', 'poa_spoa'):
            for order in (None, 'guide', 'random:3'):
                name, spec = rp.variant(method, None, order)
                realign.all_methods()[name] = realign._as_method(name, spec)
                out = os.path.join(tmp, name + '.msa.fa')
                info = realign.align_fasta(name, fa, out, timeout=120, mem_mb=4000, workdir=tmp)
                check(info['status'] == 'ok', '%s: status ok (%s)' % (name, info.get('message', '')))
                check(info['n_distinct'] == len(recs) - 1 and info['masked_bp'] == 30,
                      '%s: duplicate aligned once, N run masked' % name)
                rows = msa_graph.read_msa(out)
                check([n for n, _ in rows] == [n for n, _ in recs], '%s: rows in input order, named by first word' % name)
                check(all(r.replace('-', '') == s for (_, r), (_, s) in zip(rows, recs)), '%s: rows spell the input' % name)
                check(rows[3][1] == rows[-1][1], '%s: duplicate rows identical' % name)
                gfa = os.path.join(tmp, name + '.gfa')
                st = msa_graph.msa_to_gfa(out, fa, gfa, workdir=tmp)
                check(st['paths_checked'] == len(recs) and st['N_nodes'] == 1,
                      '%s: graph paths spell the input; the N run is one node' % name)
        # projection onto a subset named differently (full panel -> hap32), matched by sequence
        sub = os.path.join(tmp, 'sub.fa')
        with open(sub, 'w') as f:
            for n, s in [recs[0], recs[2], recs[5]]:
                f.write('>renamed_%s\n%s\n' % (n, s))
        proj = os.path.join(tmp, 'proj.msa.fa')
        pr = rp.project_to_hap32(os.path.join(tmp, 'poa_abpoa.msa.fa'), sub, proj)
        prow = msa_graph.read_msa(proj)
        check(pr['rows'] == 3 and [n for n, _ in prow] == ['renamed_' + recs[i][0] for i in (0, 2, 5)],
              'projection keeps the subset by sequence, with the subset names')
        check(all(any(r[c] != '-' for _, r in prow) for c in range(len(prow[0][1]))), 'projection drops all-gap columns')
        st = msa_graph.msa_to_gfa(proj, sub, os.path.join(tmp, 'proj.gfa'), workdir=tmp)
        check(st['paths_checked'] == 3, 'projected MSA builds a graph')
        with open(sub, 'a') as f:
            f.write('>not_there\nACGTACGT\n')
        try:
            rp.project_to_hap32(os.path.join(tmp, 'poa_abpoa.msa.fa'), sub, proj)
            check(False, 'projection refuses a sequence with no identical row')
        except ValueError:
            check(True, 'projection refuses a sequence with no identical row')
        # predicted-memory refusal
        name, spec = rp.variant('poa_spoa', None, None)
        realign.all_methods()[name] = realign._as_method(name, spec)
        info = realign.align_fasta(name, fa, os.path.join(tmp, 'x.msa.fa'), timeout=120, mem_mb=1, workdir=tmp)
        check(info['status'] == 'memout' and 'predicted' in info.get('message', ''), 'predicted memory over the cap: not run')
        # variant names
        check(rp.variant('poa_abpoa', 'p', 'given')[0] == ('poa_abpoa' if rp.DEFAULTS['poa_abpoa']['config'] == 'p'
                                                          else 'poa_abpoa_p_given'), 'variant naming')
        check(rp.variant('poa_spoa', None, 'random:2')[0] == 'poa_spoa_random2', 'spoa variant naming')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print('all tests passed')


if __name__ == '__main__':
    main()
