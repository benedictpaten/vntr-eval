#!/usr/bin/env python3
"""Offline tests for tools/msa_graph.py.
    python3 tools/test_msa_graph.py            (vg is used for the cross-check when installed)"""
import os
import random
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import msa_graph as mg  # noqa: E402

TMP = tempfile.mkdtemp(prefix='test_msa_graph.')


def fa(path, recs):
    with open(path, 'w') as f:
        for n, s in recs:
            f.write('>%s\n%s\n' % (n, s))
    return path


def random_msa(rng, nrows, ncols, n_run=False):
    base = [rng.choice('ACGT') for _ in range(ncols)]
    rows = []
    for r in range(nrows):
        row = []
        for c in range(ncols):
            x = rng.random()
            row.append('-' if x < 0.08 else rng.choice('ACGT') if x < 0.15 else base[c])
        rows.append(row)
    if n_run and ncols >= 30:
        for c in range(5, 25):
            rows[1][c] = 'N'
    names = ['CHM13#0#chr1'] + ['recombination#%d#chr1#0' % r for r in range(1, nrows - 1)] + ['GRCh38#0#chr1[12345]']
    return [(names[r], ''.join(rows[r])) for r in range(nrows)]


def gfa_paths(path):
    seq, paths, edges = {}, {}, []
    for line in open(path):
        x = line.rstrip('\n').split('\t')
        if x[0] == 'S':
            seq[x[1]] = x[2]
        elif x[0] == 'L':
            edges.append((int(x[1]), int(x[3])))
        elif x[0] == 'P':
            paths[x[1]] = [s[:-1] for s in x[2].split(',')]
    return seq, paths, edges


def test_random_msas_spell_and_match_vg():
    rng = random.Random(5)
    have_vg = os.path.exists(mg.VG) or shutil.which(mg.VG)
    for k in range(25):
        msa = random_msa(rng, rng.randint(2, 12), rng.randint(30, 300) if k % 3 == 0 else rng.randint(5, 300),
                         n_run=(k % 3 == 0))
        if any(not s.replace('-', '') for _, s in msa):
            continue
        hap = [(n, s.replace('-', '')) for n, s in msa]
        mp = fa(os.path.join(TMP, 'm%d.fa' % k), msa)
        hp = fa(os.path.join(TMP, 'h%d.fa' % k), hap)
        out = os.path.join(TMP, 'g%d.gfa' % k)
        st = mg.msa_to_gfa(mp, hp, out, engine='auto')
        seq, paths, edges = gfa_paths(out)
        for n, s in hap:
            assert ''.join(seq[i] for i in paths[n]) == s
        assert all(a < b for a, b in edges)
        if have_vg:
            assert st['vg_matches_native'], st
        # blocks mode also spells everything and never has more nodes
        st2 = mg.msa_to_gfa(mp, hp, out + '.b', merge_blocks=True, block_max=10 ** 6)
        seq, paths, _ = gfa_paths(out + '.b')
        for n, s in hap:
            assert ''.join(seq[i] for i in paths[n]) == s
        assert st2['nodes'] <= st['columns_mode']['nodes'] or st2['blocks_mode']['blocks_merged'] == 0
        if k % 3 == 0:
            assert st['N_bp_in_nodes'] >= 20


def test_unchop_keeps_path_ends():
    # a path that ends where another continues must not be merged away
    msa = [('CHM13#0#c', 'ACGTACGT'), ('x', 'ACGT----'), ('y', '----ACGT')]
    hap = [(n, s.replace('-', '')) for n, s in msa]
    out = os.path.join(TMP, 'ends.gfa')
    mg.msa_to_gfa(fa(os.path.join(TMP, 'ends.fa'), msa), fa(os.path.join(TMP, 'ends.h.fa'), hap), out,
                  engine='native')
    seq, paths, _ = gfa_paths(out)
    assert len(seq) == 2 and paths['CHM13#0#c'] == ['1', '2'] and paths['x'] == ['1'] and paths['y'] == ['2']


def test_formats_and_errors():
    msa = [('CHM13#0#c', 'ACG.TA'), ('b', 'acgcta')]
    hap = [('CHM13#0#c', 'ACGTA'), ('b', 'ACGCTA')]
    hp = fa(os.path.join(TMP, 'fh.fa'), hap)
    out = os.path.join(TMP, 'f.gfa')
    mg.msa_to_gfa(fa(os.path.join(TMP, 'f.fa'), msa), hp, out)
    # PIR
    with open(os.path.join(TMP, 'f.pir'), 'w') as f:
        f.write('>P1;CHM13#0#c\nsequence\nACG-TA*\n>P1;b\nsequence\nACGCTA*\n')
    mg.msa_to_gfa(os.path.join(TMP, 'f.pir'), hp, out)
    # abPOA -r 2: consensus row is dropped
    mg.msa_to_gfa(fa(os.path.join(TMP, 'c.fa'), msa + [('Consensus_sequence', 'ACGCTA')]), hp, out)
    for bad, why in (([('CHM13#0#c', 'ACG-TA'), ('b', 'ACGCTT')], 'spell'),
                     ([('CHM13#0#c', 'ACG-TA')], 'no MSA row'),
                     ([('CHM13#0#c', 'ACG-TA'), ('b', 'ACGCTA'), ('z', 'ACGCTA')], 'not in hap32'),
                     ([('CHM13#0#c', 'ACG-TA'), ('b', 'ACGCTA-')], 'lengths')):
        try:
            mg.msa_to_gfa(fa(os.path.join(TMP, 'bad.fa'), bad), hp, out)
        except mg.MsaGraphError as e:
            assert why in str(e), (why, e)
        else:
            raise AssertionError('accepted a bad MSA (%s)' % why)


if __name__ == '__main__':
    try:
        for name, fn in list(globals().items()):
            if name.startswith('test_'):
                fn()
                print('ok', name)
    finally:
        shutil.rmtree(TMP, ignore_errors=True)
