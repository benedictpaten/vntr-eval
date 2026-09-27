#!/usr/bin/env python3
"""Offline tests of tools/panel.py on a synthetic region (a few seconds, no big data).

    python3 tools/test_panel.py

The synthetic full-graph query has the shapes the real ones have: a CHM13 walk shared by
CHM13, its gref_CHM13 copy and two samples (weight 4); a GRCh38 walk that samples also take;
a sample walking the region in reverse; a walk that touches only the left anchor; and flanks
outside the anchors. union/project are tested for their assertions, Panel for the reference
subtraction, the walk-to-record assignment and the written graph, stats and premise for their
classifications."""
import gzip
import json
import os
import sys
import tempfile

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import panel  # noqa: E402

L_SEQ, R_SEQ = 'ACGT', 'GGTT'
SEQ = {1: L_SEQ, 2: 'A', 3: 'C', 4: R_SEQ, 5: 'TTTT', 6: 'CCCC', 7: 'G'}
CHM, G38, S3 = L_SEQ + 'A' + R_SEQ, L_SEQ + 'C' + R_SEQ, L_SEQ + 'G' + R_SEQ


def write(path, text):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'wt') as f:
        f.write(text)


def make_region(root, rid='T0001', hap_extra=(), hprc_extra=()):
    rd = os.path.join(root, 'regions', rid)
    os.makedirs(rd)
    hap = [('CHM13#0#chrT', CHM), ('GRCh38#0#chrT[100]', G38), ('recombination#1#chrT#0', G38),
           ('recombination#2#chrT#0', CHM)] + list(hap_extra)
    write(os.path.join(rd, 'hap32.fa'), ''.join('>%s len=%d\n%s\n' % (n, len(s), s) for n, s in hap))
    hprc = [('CHM13#0#chrT', 'reference', CHM), ('GRCh38#0#chrT[100]', 'reference', G38)]
    hprc += [('hprc#%d' % k, 'group=1 group_size=4', G38) for k in (1, 2, 3, 4)]
    hprc += [('hprc#%d' % k, 'group=2 group_size=2 same_as=CHM13', CHM) for k in (5, 6)]
    hprc += list(hprc_extra)
    write(os.path.join(rd, 'hprc.fa.gz'), ''.join('>%s %s\n%s\n' % r for r in hprc))
    mc = ['H\tVN:Z:1.0'] + ['S\t%d\t%s' % (i, SEQ[i]) for i in (1, 2, 3, 4)]
    mc += ['L\t1\t+\t2\t+\t0M', 'L\t1\t+\t3\t+\t0M', 'L\t2\t+\t4\t+\t0M', 'L\t3\t+\t4\t+\t0M']
    mc += ['P\tCHM13#0#chrT\t1+,2+,4+\t*', 'P\tGRCh38#0#chrT[100]\t1+,3+,4+\t*',
           'P\trecombination#1#chrT#0\t1+,3+,4+\t*', 'P\trecombination#2#chrT#0\t1+,2+,4+\t*']
    write(os.path.join(rd, 'mc.gfa'), '\n'.join(mc) + '\n')
    rj = {'region_id': rid, 'stratum': 'test', 'contig': 'chrT', 'span_start': 101, 'span_end': 109,
          'anchor_left': '1+', 'anchor_right': '4+', 'anchor_left_seq': L_SEQ, 'anchor_right_seq': R_SEQ,
          'n_hprc': 6, 'hprc': {'status': 'ok', 'n_sample_haplotypes_spanning': 6,
                                'n_distinct_sample_sequences': 2, 'n_distinct_in_hap32': 2}}
    json.dump(rj, open(os.path.join(rd, 'region.json'), 'w'))
    qd = os.path.join(root, 'query')
    os.makedirs(qd, exist_ok=True)
    q = ['H\tVN:Z:1.1\tRS:Z:CHM13'] + ['S\t%d\t%s' % (i, s) for i, s in SEQ.items()]
    q += ['L\t5\t+\t1\t+\t0M', 'L\t1\t+\t2\t+\t0M', 'L\t1\t+\t3\t+\t0M', 'L\t1\t+\t7\t+\t0M',
          'L\t2\t+\t4\t+\t0M', 'L\t3\t+\t4\t+\t0M', 'L\t7\t+\t4\t+\t0M', 'L\t4\t+\t6\t+\t0M']
    q += ['W\tCHM13\t0\tchrT\t96\t113\t>5>1>2>4>6\tWT:i:4',        # CHM13 + gref copy + 2 samples
          'W\tunknown\t1\tchrT\t0\t13\t>1>3>4\tWT:i:1',              # GRCh38 (clipped window)
          'W\tunknown\t2\tchrT\t0\t17\t>5>1>3>4>6\tWT:i:3',          # 3 samples
          'W\tunknown\t3\tchrT\t0\t17\t<6<4<3<1<5\tWT:i:1',          # 1 sample, reverse strand
          'W\tunknown\t4\tchrT\t0\t9\t>5>1>2\tWT:i:1']               # left anchor only
    write(os.path.join(qd, rid + '.gfa.gz'), '\n'.join(q) + '\n')
    return rd, qd


def expect_error(fn, *a, **k):
    try:
        fn(*a, **k)
    except panel.PanelError as e:
        return str(e)
    raise AssertionError('expected PanelError from %s' % fn.__name__)


def main():
    with tempfile.TemporaryDirectory() as root:
        rd, qd = make_region(root)
        # union: CHM13 first, members, weights
        recs, rows = panel.union(rd)
        assert [r[2] for r in recs] == [CHM, G38], recs
        assert rows[0]['weight'] == 3 and rows[0]['n_hap32'] == 2, rows[0]
        assert rows[1]['weight'] == 5 and rows[1]['n_panel_samples'] == 4, rows[1]
        assert rows[0]['members'].split(',')[:2] == ['CHM13#0#chrT', 'recombination#2#chrT#0']
        fa, mp = os.path.join(root, 'u.fa'), os.path.join(root, 'u.map.tsv')
        panel.union_to_files(rd, fa, mp)
        # project: a full MSA with an extra column that is all-gap among hap32 rows
        msa = os.path.join(root, 'full.msa.fa')
        write(msa, '>u0001\nACGTA-GGTT\n>u0002\nACGT-CGGTT\n>Consensus_sequence\nACGTACGGTT\n')
        out = os.path.join(root, 'proj.msa.fa')
        st = panel.project(msa, mp, os.path.join(rd, 'hap32.fa'), out)
        assert st['rows_out'] == 4 and st['columns_out'] == 10, st
        rows_out = dict((n, s) for n, _, s in panel.read_fasta(out))
        assert rows_out['recombination#1#chrT#0'] == 'ACGT-CGGTT'
        write(msa, '>u0001\nACGTA-GGTT-\n>u0002\nACGT-CGGTTT\n')       # u0002 spells a longer sequence
        assert 'do not spell' in expect_error(panel.project, msa, mp, os.path.join(rd, 'hap32.fa'), out)
        write(msa, '>u0001\nACGTAGGTT\n')                                # u0002 missing
        assert 'no row' in expect_error(panel.project, msa, mp, os.path.join(rd, 'hap32.fa'), out)
        # all-gap columns among the hap32 rows are dropped
        write(msa, '>u0001\nACGTA--GGTT\n>u0002\nACGT-C-GGTT\n>u0003\nACGT--GGGTT\n')
        write(mp + '.3', open(mp).read() + 'u0003\t9\t1\t0\t1\t0\thprc#9\n')
        st = panel.project(msa, mp + '.3', os.path.join(rd, 'hap32.fa'), out)
        assert st['all_gap_columns_dropped'] == 1, st
        # panel-graph: one path per hprc.fa.gz record, rows may be named by a member name
        write(msa, '>CHM13#0#chrT\nACGTA-GGTT\n>u0002\nACGT-CGGTT\n')
        pg = os.path.join(root, 'panel.gfa')
        st = panel.panel_graph(msa, mp, rd, pg, engine='native')
        assert st['panel_paths'] == 8 and st['nodes'] == 4, st
        _, _, pp = panel.read_gfa_paths(pg)
        assert list(pp)[:3] == ['CHM13#0#chrT', 'GRCh38#0#chrT[100]', 'hprc#1'] and pp['hprc#1'] == pp['hprc#4']
        write(mp + '.bad', open(mp).read().replace('hprc#5', 'hprc#7'))
        assert 'does not match' in expect_error(panel.panel_graph, msa, mp + '.bad', rd, pg, engine='native')

        # Panel: reference subtraction and walk assignment
        P = panel.Panel(rd, qd)
        assert P.chm13_walk == (1, 2, 4) and P.chm13_walk_weight == 4
        assert dict(P.sample_walks) == {(1, 2, 4): 2, (1, 3, 4): 4}, P.sample_walks
        assert P.other == {'left_anchor_only': 1}, P.other
        assert P.walk_of['hprc#1'] == (1, 3, 4) and P.walk_of['hprc#6'] == (1, 2, 4)
        g = os.path.join(root, 'mc_full.gfa')
        r = P.write_gfa(g)
        assert r['paths'] == 8 and r['nodes'] == 4 and r['distinct_walks'] == 2, r
        seqs, edges, paths = panel.read_gfa_paths(g)
        assert 5 not in seqs and 7 not in seqs, 'flank or unused node leaked into the span graph'
        assert (1, 7) not in edges
        # stats and premise
        s = panel.region_stats(rd, qd)
        assert s['hprc_spanning'] == 6 and s['hprc_not_spanning'] == 450 and s['hap32_not_verbatim'] == 0, s
        assert s['reference_removed'] == 'CHM13;gref_CHM13;GRCh38#0#chrT[100]'
        pr = panel.premise_region(rd, qd)
        assert pr['mc_nodes_missing_in_full'] == 0 and pr['pairs_same_cost'] == pr['pairs_compared'] == 6, pr

        # a hap32 path that copies the GRCh38 walk nobody else takes, and a junction
        root2 = os.path.join(root, 'b')
        os.makedirs(root2)
        rd2, qd2 = make_region(root2, 'T0002')
        qp = os.path.join(qd2, 'T0002.gfa.gz')
        txt = gzip.open(qp, 'rt').read().replace('>5>1>3>4>6\tWT:i:3', '>5>1>7>4>6\tWT:i:3')
        txt = txt.replace('<6<4<3<1<5\tWT:i:1', '<6<4<7<1<5\tWT:i:1')
        write(qp, txt)
        hp = os.path.join(rd2, 'hprc.fa.gz')
        write(hp, gzip.open(hp, 'rt').read().replace(G38 + '\n', S3 + '\n').replace(
            '>GRCh38#0#chrT[100] reference\n' + S3, '>GRCh38#0#chrT[100] reference\n' + G38))
        s2 = panel.region_stats(rd2, qd2)
        assert s2['hap32_not_verbatim_kinds'] == 'copy_of_GRCh38=1', s2['hap32_not_verbatim_kinds']
        # excluded samples
        root3 = os.path.join(root, 'c')
        os.makedirs(root3)
        rd3, _ = make_region(root3, 'T0003', hprc_extra=[('HG002#1#chrT', 'x', CHM)])
        assert 'excluded' in expect_error(panel.union, rd3)
        print('union / project: ok')
        print('panel reconstruction, graph, stats, premise: ok')
        print('all tests passed')


if __name__ == '__main__':
    main()
