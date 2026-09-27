#!/usr/bin/env python3
"""Offline tests for tools/realign.py (a few seconds; mafft is used when installed).
    python3 tools/test_realign.py"""
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import realign as ra  # noqa: E402
import msa_graph as mg  # noqa: E402

TMP = tempfile.mkdtemp(prefix='test_realign.')
HAVE_MAFFT = bool(shutil.which(ra.config.MAFFT) or os.path.exists(ra.config.MAFFT))


def fa(path, recs):
    with open(path, 'w') as f:
        for n, s in recs:
            f.write('>%s some description\n' % n)
            for i in range(0, len(s), 60):          # wrapped, to exercise the reader
                f.write(s[i:i + 60] + '\n')
    return path


def test_mask_runs():
    assert ra.mask_runs('ACGT') == ('ACGT', [])
    assert ra.mask_runs('ACNNNGT') == ('ACGT', [(2, 'NNN')])
    assert ra.mask_runs('NNACGTRN') == ('ACGT', [(0, 'NN'), (4, 'RN')])
    assert ra.mask_runs('ACNGTNNNA', min_run=2) == ('ACNGTA', [(5, 'NNN')])
    assert ra.mask_runs('NNNN') == ('', [(0, 'NNNN')])


def test_reinsert_runs_spells_and_is_exclusive():
    rng = random.Random(3)
    for _ in range(200):
        nrow = rng.randint(1, 6)
        ncol = rng.randint(1, 40)
        rows = [''.join(rng.choice('ACGT--') for _ in range(ncol)) for _ in range(nrow)]
        runs, want = {}, []
        for r, row in enumerate(rows):
            m = row.replace('-', '')
            rs = []
            # random runs at random positions of the ungapped row
            pos = sorted(set(rng.randint(0, len(m)) for _ in range(rng.randint(0, 3))))
            orig, last = [], 0
            for p in pos:
                run = 'N' * rng.randint(1, 5)
                orig.append(m[last:p] + run)
                rs.append((p, run))
                last = p
            orig.append(m[last:])
            runs[r] = rs
            want.append(''.join(orig))
        out = ra.reinsert_runs(rows, runs)
        assert len(set(len(x) for x in out)) == 1
        for r in range(nrow):
            assert out[r].replace('-', '') == want[r], (rows, runs, out)
        # every inserted column has bases in one row only, and they are Ns
        orig_cols = sum(1 for c in range(ncol))
        n_new = len(out[0]) - orig_cols
        assert n_new == sum(len(x) for rs in runs.values() for _, x in rs)
        for c in range(len(out[0])):
            col = [x[c] for x in out]
            if 'N' in col:
                assert sum(ch != '-' for ch in col) == 1


def test_run_proc_timeout_and_ok():
    r = ra.run_proc(['sleep', '30'], timeout=1.5, poll=0.3)
    assert r['status'] == 'timeout' and r['seconds'] < 10, r
    r = ra.run_proc(['true'])
    assert r['status'] == 'ok' and r['returncode'] == 0
    r = ra.run_proc(['false'])
    assert r['status'] == 'error'
    out = os.path.join(TMP, 'echo.txt')
    r = ra.run_proc(['sh', '-c', 'echo hello'], stdout_path=out)
    assert open(out).read() == 'hello\n'


def fake_align(in_fa, out_fa, threads=1, workdir=None, timeout=None, mem_mb=None, pad=0, **_):
    """Test aligner: right-pads every sequence with gaps to the longest (+pad all-gap columns)."""
    recs = mg.read_fasta(in_fa)
    L = max(len(s) for _, s in recs) + pad
    with open(out_fa, 'w') as f:
        for n, s in recs:
            f.write('>%s\n%s\n' % (n, (s + '-' * (L - len(s))).lower()))
    return {'status': 'ok', 'command': ['fake']}


def bad_align(in_fa, out_fa, **_):
    raise RuntimeError('boom')


def test_pipeline_with_fake_method():
    m = ra.Method('fake', fake_align, {'pad': 2})
    recs = [('CHM13#0#chr1', 'ACGTACGTAAAC'), ('a', 'ACGTNNNNACGTAAAC'), ('b', 'ACGTACGTAAAC'),
            ('c', 'acgtac'), ('d', 'NNNN')]
    inp = fa(os.path.join(TMP, 'in.fa'), recs)
    out = os.path.join(TMP, 'out.msa.fa')
    info = ra.align_fasta(m, inp, out, workdir=TMP)
    assert info['status'] == 'ok', info
    assert info['n_seqs'] == 5 and info['n_distinct'] == 4       # CHM13 and b are identical
    assert info['masked_runs'] == 2 and info['masked_bp'] == 8
    rows = mg.read_msa(out)
    assert [n for n, _ in rows] == [n for n, _ in recs]
    for (n, s), (_, r) in zip(recs, rows):
        assert r.replace('-', '') == s.upper()
    assert rows[0][1] == rows[2][1]
    # the graph builder accepts it, and the N run is one node of its own
    hap = fa(os.path.join(TMP, 'hap.fa'), recs)
    st = mg.msa_to_gfa(out, hap, os.path.join(TMP, 'out.gfa'), engine='native')
    assert st['N_nodes'] == 2 and st['N_bp_in_nodes'] == 8, st
    # a failing plugin is an error status, not an exception
    info = ra.align_fasta(ra.Method('bad', bad_align), inp, os.path.join(TMP, 'bad.fa'), workdir=TMP)
    assert info['status'] == 'error' and 'boom' in info['message']
    assert not os.path.exists(os.path.join(TMP, 'bad.fa'))
    # duplicate names are refused
    fa(os.path.join(TMP, 'dup.fa'), [('x', 'ACGT'), ('x', 'ACGA')])
    info = ra.align_fasta(m, os.path.join(TMP, 'dup.fa'), os.path.join(TMP, 'dup.msa'), workdir=TMP)
    assert info['status'] == 'error'


def test_plugin_spec_forms():
    class Spec(object):
        align = staticmethod(fake_align)
        params = {'pad': 1}
    for spec in ({'align': fake_align, 'params': {'pad': 1}}, Spec(), fake_align):
        m = ra._as_method('x', spec)
        assert callable(m.align)


def random_haplotypes(rng, n, unit='ACGGT', copies=(3, 9)):
    left = ''.join(rng.choice('ACGT') for _ in range(40))
    right = ''.join(rng.choice('ACGT') for _ in range(40))
    out = []
    for i in range(n):
        arr = ''.join(unit if rng.random() > 0.2 else unit[:2] + 'A' + unit[3:]
                      for _ in range(rng.randint(*copies)))
        out.append(('h%d' % i, left + arr + right))
    return out


def test_mafft_end_to_end():
    if not HAVE_MAFFT:
        print('  (mafft not installed; skipped)')
        return
    rng = random.Random(11)
    recs = random_haplotypes(rng, 8)
    recs[3] = (recs[3][0], recs[3][1][:50] + 'N' * 30 + recs[3][1][50:])
    inp = fa(os.path.join(TMP, 'h.fa'), recs)
    for method in ('mafft_fftns2', 'mafft_linsi', 'mafft_einsi', 'mafft_ginsi'):
        out = os.path.join(TMP, method + '.msa.fa')
        info = ra.align_fasta(method, inp, out, threads=1, workdir=TMP, timeout=120)
        assert info['status'] == 'ok', info
        rows = mg.read_msa(out)
        for (n, s), (m, r) in zip(recs, rows):
            assert n == m and r.replace('-', '') == s
        st = mg.msa_to_gfa(out, inp, os.path.join(TMP, method + '.gfa'))
        assert st['N_nodes'] == 1 and st['N_bp_in_nodes'] == 30
    # the CLI's MSA-only mode
    out = os.path.join(TMP, 'cli.msa.fa')
    js = os.path.join(TMP, 'cli.json')
    p = subprocess.run([sys.executable, os.path.join(ra.TOOLS, 'realign.py'), 'mafft_fftns2', '--input', inp,
                        '--msa-out', out, '--threads', '1', '--json', js, '--workdir', TMP],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert json.load(open(js))['status'] == 'ok'
    assert len(mg.read_msa(out)) == len(recs)
    # a timeout is reported as such and leaves no MSA
    big = [('b%d' % i, s * 40) for i, (_, s) in enumerate(random_haplotypes(rng, 30, copies=(40, 80)))]
    inb = fa(os.path.join(TMP, 'big.fa'), big)
    outb = os.path.join(TMP, 'big.msa.fa')
    info = ra.align_fasta('mafft_ginsi', inb, outb, threads=1, workdir=TMP, timeout=1)
    assert info['status'] == 'timeout', info['status']
    assert not os.path.exists(outb)


def test_region_mode_on_a_synthetic_package():
    if not HAVE_MAFFT:
        return
    rng = random.Random(2)
    regions = os.path.join(TMP, 'regions')
    rd = os.path.join(regions, 'TEST1')
    os.makedirs(rd)
    recs = [('CHM13#0#chr1', random_haplotypes(rng, 1)[0][1])] + random_haplotypes(rng, 5)
    fa(os.path.join(rd, 'hap32.fa'), recs)
    fa(os.path.join(rd, 'hap32.fragments.fa'), [('h9:enters_L', recs[2][1][:60])])
    json.dump({'region_id': 'TEST1', 'stratum': 'test', 'span_start': 1, 'span_end': len(recs[0][1])},
              open(os.path.join(rd, 'region.json'), 'w'))
    cand = os.path.join(TMP, 'cand')
    old = ra.RUNTIME_TSV
    ra.RUNTIME_TSV = os.path.join(TMP, 'runtime.tsv')
    try:
        row = ra.realign_region('mafft_linsi', rd, threads=1, cand_root=cand, workroot=TMP, fragments=True)
        assert row['status'] == 'ok', row
        o = ra.outputs('mafft_linsi', 'TEST1', cand)
        for k in ('msa', 'gfa', 'json', 'frag'):
            assert os.path.exists(o[k]), k
        j = json.load(open(o['json']))
        assert j['graph']['paths_checked'] == 6 and j['fragments']['n_written'] == 1
        again = ra.realign_region('mafft_linsi', rd, threads=1, cand_root=cand, workroot=TMP)
        assert again['status'] == 'skipped'
        lines = open(ra.RUNTIME_TSV).read().splitlines()
        assert len(lines) == 2 and lines[0].startswith('method\tregion_id')
    finally:
        ra.RUNTIME_TSV = old


if __name__ == '__main__':
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith('test_') and callable(v)]
    try:
        for k, f in tests:
            f()
            print('ok  ', k)
    finally:
        shutil.rmtree(TMP, ignore_errors=True)
    print('all %d tests passed' % len(tests))
