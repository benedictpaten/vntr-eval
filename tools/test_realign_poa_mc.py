#!/usr/bin/env python3
"""Tests for tools/realign_poa_mc.py (poa_abpoa_mc): offline, a few seconds.

    python3 tools/test_realign_poa_mc.py

1. trim(): a constructed case with a known cut, then random MSAs against a brute-force search over
   every cut point (same objective and tie-breaks as poaBarAligner.c's trim()).
2. bar_align()'s two-end path with a stand-in aligner and a small banding limit: every row spells its
   sequence, the shared flanks are columns, a string shorter than 2 x limit is cut once, the middle
   of a longer one is left unaligned (and counted).
3. The matrix file and the command line are Cactus's.
4. With the vendored abPOA 1.5.4 (skipped when absent): on toy pairs, abPOA with Cactus's command
   returns an alignment whose score under Cactus's scheme (5x5 matrix, convex gaps 400+30g /
   1200+g) equals the optimum of a Python dynamic program, and abPOA's CLI defaults do not (the
   same dynamic program under abPOA's default scheme reproduces the default runs, which checks the
   scorer); then an end-to-end run with an N run and an empty interior.
"""
import os
import random
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import realign_poa_mc as mc  # noqa: E402

GAP = ord('-')


def rows_to_msa(rows):
    return [bytearray(r, 'ascii') for r in rows]


def brute_cut(msa1, cs1, row1, msa2, cs2, row2, overlap):
    """Best cut by exhaustive search: keep (len1 - overlap + c) bases of row1 and (len2 - c) of row2."""
    l1 = sum(1 for c in msa1[row1] if c != GAP)
    l2 = sum(1 for c in msa2[row2] if c != GAP)

    def kept(msa, cs, row, k):
        s, n = 0, 0
        for i, c in enumerate(msa[row]):
            if c != GAP:
                if n < k:
                    s += cs[i]
                n += 1
        return s
    scores = [kept(msa1, cs1, row1, l1 - overlap + c) + kept(msa2, cs2, row2, l2 - c) for c in range(overlap + 1)]
    best = 0
    for c in range(1, overlap):                 # strictly better than the cut-at-0 start
        if scores[c] > scores[best]:
            best = c
    if scores[overlap] > scores[best]:
        best = overlap
    return best


def test_trim_constructed():
    # One sequence ACGTACGT (8 bp) is in both MSAs in full (overlap 8). In msa1 (forward) its first
    # 6 bases sit in columns shared with 2 other rows, its last 2 alone; in msa2 (reverse complement,
    # ACGTACGT is its own reverse complement) its first 2 bases are shared with 1 other row, the rest
    # alone. Keeping 6 bases in msa1 and 2 in msa2 keeps every shared column: cut point 6.
    msa1 = rows_to_msa(['ACGTACGT', 'ACGTAC--', 'ACGTAC--'])
    msa2 = rows_to_msa(['ACGTACGT', 'AC------'])
    cs1, cs2 = mc.make_column_scores(msa1), mc.make_column_scores(msa2)
    assert cs1 == [2, 2, 2, 2, 2, 2, 0, 0] and cs2 == [1, 1, 0, 0, 0, 0, 0, 0]
    cut = mc.trim(0, msa1, cs1, 0, msa2, cs2, 8)
    assert cut == 6, cut
    assert msa1[0].decode() == 'ACGTAC--' and msa2[0].decode() == 'AC------'
    assert cs1 == [2, 2, 2, 2, 2, 2, 0, 0] and cs2 == [1, 1, 0, 0, 0, 0, 0, 0]
    # a tie (all columns private): the C code keeps the cut at 0, i.e. all of the overlap in msa2
    m1, m2 = rows_to_msa(['ACGT']), rows_to_msa(['ACGT'])
    c1, c2 = mc.make_column_scores(m1), mc.make_column_scores(m2)
    assert mc.trim(0, m1, c1, 0, m2, c2, 4) == 0
    assert m1[0].decode() == '----' and m2[0].decode() == 'ACGT'
    # partial overlap: a 14 bp sequence, 10 bp prefixes at both ends -> overlap 6
    s = 'AAAACCCCGGGGTT'
    msa1 = rows_to_msa([s[:10], s[:10]])
    msa2 = rows_to_msa([mc.revcomp(s)[:10], '-' * 10])
    cs1, cs2 = mc.make_column_scores(msa1), mc.make_column_scores(msa2)
    cut = mc.trim(0, msa1, cs1, 0, msa2, cs2, 6)
    assert cut == 6                              # msa1's shared columns win: all of the overlap kept there
    assert msa1[0].decode() == s[:10] and msa2[0].decode() == mc.revcomp(s)[:4] + '-' * 6
    print('ok trim, constructed cases')


def random_msa(rng, nrow, ncol, p_gap):
    rows = []
    for _ in range(nrow):
        rows.append(''.join('-' if rng.random() < p_gap else rng.choice('ACGT') for _ in range(ncol)))
    return rows


def test_trim_random():
    rng = random.Random(7)
    n = 0
    for trial in range(400):
        nrow = rng.randint(1, 5)
        rows1 = random_msa(rng, nrow, rng.randint(3, 25), rng.choice([0.1, 0.3, 0.6]))
        rows2 = random_msa(rng, nrow, rng.randint(3, 25), rng.choice([0.1, 0.3, 0.6]))
        r = rng.randrange(nrow)
        l1 = len(rows1[r].replace('-', ''))
        l2 = len(rows2[r].replace('-', ''))
        if min(l1, l2) == 0:
            continue
        ov = rng.randint(1, min(l1, l2))
        m1, m2 = rows_to_msa(rows1), rows_to_msa(rows2)
        c1, c2 = mc.make_column_scores(m1), mc.make_column_scores(m2)
        want = brute_cut(m1, c1, r, m2, c2, r, ov)
        cut = mc.trim(r, m1, c1, r, m2, c2, ov)
        assert cut == want, (trial, cut, want)
        assert len(m1[r].replace(b'-', b'')) == l1 - ov + cut
        assert len(m2[r].replace(b'-', b'')) == l2 - cut
        assert c1 == mc.make_column_scores(m1) and c2 == mc.make_column_scores(m2)
        n += 1
    print('ok trim, %d random cases match brute force' % n)


def fake_poa(strings, tag=None, calls=None, **_):
    """A stand-in aligner: left-justified rows (a valid MSA; alignment quality is irrelevant here)."""
    if calls is not None:
        calls.append({'tag': tag, 'n': len(strings), 'progressive': True, 'seconds': 0, 'peak_rss_mb': 0,
                      'sampled_rss_mb': 0, 'predicted_mb': 0, 'command': ['fake']})
    m = max(len(s) for s in strings)
    return [s + '-' * (m - len(s)) for s in strings]


def test_bar_ends():
    old = mc.BANDING_LIMIT
    mc.BANDING_LIMIT = 10
    try:
        rng = random.Random(3)
        L, R = 'GA', 'T'
        inner = [''.join(rng.choice('ACGT') for _ in range(k)) for k in (4, 9, 10, 14, 19, 20, 21, 33, 0)]
        seqs = [L + x + R for x in inner]
        rows, info = mc.bar_align(seqs, (2, 1), poa=fake_poa)
        assert info['bar_rule'] == 'ends' and info['abpoa_calls'] == 2
        assert all(r.replace('-', '') == s for r, s in zip(rows, seqs))
        assert len(set(len(r) for r in rows)) == 1
        assert all(r[:2] == L and r[-1] == R for r in rows)
        assert info['rows_ge_10kb'] == 6 and info['rows_gt_20kb'] == 2
        assert info['unaligned_bp'] == (21 - 20) + (33 - 20), info
        # all lengths below the limit: one call on whole strings
        rows, info = mc.bar_align([L + x + R for x in inner[:2]], (2, 1), poa=fake_poa)
        assert info['bar_rule'] == 'single' and info['abpoa_calls'] == 1 and info['unaligned_bp'] == 0
    finally:
        mc.BANDING_LIMIT = old
    print('ok BAR two-end assembly (stand-in aligner)')


def test_matrix_and_command():
    d = tempfile.mkdtemp()
    p = os.path.join(d, 'm.mtx')
    mc.write_matrix(p)
    lines = open(p).read().split('\n')
    assert lines[0] == '\tA\tC\tG\tT\tN'
    assert lines[1] == 'A\t91\t-114\t-61\t-123\t-100'
    assert lines[5] == 'N\t-100\t-100\t-100\t-100\t100'
    cmd = mc.abpoa_command('in.fa', 'm.mtx', True)
    assert cmd[1:] == ['in.fa', '-O', '400,1200', '-E', '30,1', '-b', '1000', '-f', '0.1', '-t', 'm.mtx', '-r', '1',
                       '-m', '0', '-k', '15', '-w', '5', '-n', '500', '-p'], cmd
    assert mc.use_progressive([10, 5, 1]) and not mc.use_progressive([1] * 5001)
    print('ok matrix file and command line')


# ---------------------------------------------------------------- toy pairs against a dynamic program

IDX = {c: i for i, c in enumerate('ACGTN')}
MC_SCHEME = (mc.SUB_MATRIX, (400, 30), (1200, 1))
DEFAULT_SCHEME = ([2 if (i == j and i < 4) else (0 if (i == 4 or j == 4) else -4) for i in range(5) for j in range(5)],
                  (4, 2), (24, 1))


def gap_cost(g, scheme):
    (o1, e1), (o2, e2) = scheme[1], scheme[2]
    return min(o1 + g * e1, o2 + g * e2)


def aln_score(r1, r2, scheme):
    mat = scheme[0]
    s, i, n = 0, 0, len(r1)
    while i < n:
        a, b = r1[i], r2[i]
        if a != '-' and b != '-':
            s += mat[IDX[a] * 5 + IDX[b]]
            i += 1
        elif a == '-' and b == '-':
            i += 1
        else:
            which = 0 if a == '-' else 1
            j = i
            while j < n and (r1[j] == '-') == (which == 0) and (r2[j] == '-') == (which == 1):
                j += 1
            s -= gap_cost(j - i, scheme)
            i = j
    return s


def dp_opt(a, b, scheme):
    """Global alignment optimum under a convex (two-piece affine) gap scheme (Gotoh, two gap models)."""
    mat = scheme[0]
    (o1, e1), (o2, e2) = scheme[1], scheme[2]
    NEG = -10 ** 12
    n, m = len(a), len(b)
    H = [[NEG] * (m + 1) for _ in range(n + 1)]
    E = [[[NEG] * (m + 1) for _ in range(n + 1)] for _ in range(2)]
    F = [[[NEG] * (m + 1) for _ in range(n + 1)] for _ in range(2)]
    ops = ((o1, e1), (o2, e2))
    H[0][0] = 0
    for j in range(1, m + 1):
        for k, (o, e) in enumerate(ops):
            E[k][0][j] = -(o + j * e)
        H[0][j] = max(E[0][0][j], E[1][0][j])
    for i in range(1, n + 1):
        for k, (o, e) in enumerate(ops):
            F[k][i][0] = -(o + i * e)
        H[i][0] = max(F[0][i][0], F[1][i][0])
        for j in range(1, m + 1):
            best = H[i - 1][j - 1] + mat[IDX[a[i - 1]] * 5 + IDX[b[j - 1]]]
            for k, (o, e) in enumerate(ops):
                E[k][i][j] = max(H[i][j - 1] - (o + e), E[k][i][j - 1] - e)
                F[k][i][j] = max(H[i - 1][j] - (o + e), F[k][i - 1][j] - e)
                best = max(best, E[k][i][j], F[k][i][j])
            H[i][j] = best
    return H[n][m]


def toy_pairs(rng, k=14):
    pairs = []
    for t in range(k):
        n = rng.randint(40, 90)
        if t % 2:
            unit = ''.join(rng.choice('ACGT') for _ in range(rng.randint(3, 7)))
            a = ''.join(rng.choice('ACGT') for _ in range(10)) + unit * rng.randint(5, 9) + \
                ''.join(rng.choice('ACGT') for _ in range(10))
        else:
            a = ''.join(rng.choice('ACGT') for _ in range(n))
        b = list(a)
        for _ in range(rng.randint(2, 6)):         # substitutions: transitions and transversions
            p = rng.randrange(len(b))
            b[p] = rng.choice('ACGT')
        for _ in range(rng.randint(1, 3)):         # indels
            p = rng.randrange(len(b))
            if rng.random() < 0.5:
                del b[p:p + rng.randint(1, 8)]
            else:
                b[p:p] = list(''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 8))))
        pairs.append((a, ''.join(b)))
    return pairs


def run_abpoa(a, b, cmd_tail, workdir):
    import subprocess
    fa = os.path.join(workdir, 'pair.fa')
    with open(fa, 'w') as f:
        f.write('>a\n%s\n>b\n%s\n' % (a, b))
    p = subprocess.run([mc.ABPOA_MC, fa] + cmd_tail, capture_output=True, text=True, check=True)
    rows = {}
    name = None
    for ln in p.stdout.splitlines():
        if ln.startswith('>'):
            name = ln[1:].split()[0]
            rows[name] = ''
        elif name:
            rows[name] += ln.strip()
    return rows['a'].upper(), rows['b'].upper()


def _write_mat(path, mat):
    with open(path, 'w') as f:
        f.write('\tA\tC\tG\tT\tN\n')
        for i in range(5):
            f.write('ACGTN'[i] + ''.join('\t%d' % mat[i * 5 + j] for j in range(5)) + '\n')


def test_toy_pairs_real_abpoa():
    """Cactus's command is optimal under Cactus's scheme on every toy pair, and the check can fire:
    each of three mis-specified commands (abPOA defaults; Cactus gaps with abPOA's default scores;
    Cactus's command with HOXD70's transition score -31 in place of -61) is caught on some pair."""
    if not os.path.exists(mc.ABPOA_MC):
        print('skip toy pairs: no %s' % mc.ABPOA_MC)
        return
    d = tempfile.mkdtemp()
    mat = os.path.join(d, 'm.mtx')
    mc.write_matrix(mat)
    alt = list(mc.SUB_MATRIX)
    for i, j in ((0, 2), (2, 0), (1, 3), (3, 1)):
        alt[i * 5 + j] = -31
    mat31 = os.path.join(d, 'ts31.mtx')
    _write_mat(mat31, alt)
    mc_tail = mc.abpoa_command('X', mat, True)[2:]
    wrong = {'abPOA defaults': ['-r', '1', '-m', '0'],
             'Cactus gaps, default scores': ['-O', '400,1200', '-E', '30,1', '-b', '1000', '-f', '0.1', '-r', '1', '-m', '0'],
             'transitions -31': [mat31 if x == mat else x for x in mc_tail]}
    caught = dict((k, 0) for k in wrong)
    pairs = toy_pairs(random.Random(12), 30)
    for a, b in pairs:
        r1, r2 = run_abpoa(a, b, mc_tail, d)
        assert r1.replace('-', '') == a and r2.replace('-', '') == b
        opt = dp_opt(a, b, MC_SCHEME)
        got = aln_score(r1, r2, MC_SCHEME)
        assert got == opt, ('Cactus scheme: abPOA %d vs optimum %d' % (got, opt), r1, r2)
        for k, tail in wrong.items():
            w1, w2 = run_abpoa(a, b, tail, d)
            if k == 'abPOA defaults':
                assert aln_score(w1, w2, DEFAULT_SCHEME) == dp_opt(a, b, DEFAULT_SCHEME), 'scorer disagrees with abPOA'
            caught[k] += aln_score(w1, w2, MC_SCHEME) < opt
    assert all(caught.values()), caught
    print('ok toy pairs: Cactus command optimal under Cactus scheme on %d/%d; mis-specified commands caught on %s'
          % (len(pairs), len(pairs), ', '.join('%s %d' % kv for kv in caught.items())))


def test_end_to_end_real_abpoa():
    if not os.path.exists(mc.ABPOA_MC):
        print('skip end to end: no %s' % mc.ABPOA_MC)
        return
    import msa_graph
    rng = random.Random(5)
    core = ''.join(rng.choice('ACGT') for _ in range(300))
    seqs = [('CHM13#0', 'G' + core + 'T'),
            ('h1', 'G' + core[:100] + core[150:] + 'T'),
            ('h2', 'G' + core[:120] + 'NNNNNNNNNN' + core[120:] + 'T'),     # an N run
            ('h3', 'GT'),                                                    # empty interior
            ('h4', 'G' + core + 'T')]                                        # a duplicate (kept as a row)
    d = tempfile.mkdtemp()
    fa = os.path.join(d, 'in.fa')
    with open(fa, 'w') as f:
        for n, s in seqs:
            f.write('>%s\n%s\n' % (n, s))
    out = os.path.join(d, 'out.msa.fa')
    info = mc.align_fasta(fa, out, flanks=(1, 1), workdir=d)
    assert info['status'] == 'ok', info
    rows = dict(msa_graph.read_msa(out))
    assert all(rows[n].replace('-', '') == s for n, s in seqs)
    assert info['n_distinct'] == 5 and info['aligner']['calls'][0]['n'] == 5     # no dedup
    assert info['aligner']['calls'][0]['empty'] == 1
    assert info['masked_runs'] == 1
    print('ok end to end (N run, empty interior, duplicate row)')


if __name__ == '__main__':
    test_trim_constructed()
    test_trim_random()
    test_bar_ends()
    test_matrix_and_command()
    test_toy_pairs_real_abpoa()
    test_end_to_end_real_abpoa()
    print('all tests passed')
