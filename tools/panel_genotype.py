#!/usr/bin/env python3
"""panel_genotype.py -- prototype: genotype a repeat region as one site whose alleles are the panel
haplotypes, scored with the region's reads.

For each region: the distinct hap32 panel sequences (anchor to anchor) are the candidate alleles. The
reads the genome-wide mapping placed in the region's window are kept if they carry a variable 21-mer
(one not shared by every distinct haplotype, on either strand); each kept read gets its semi-global edit
distance d to every haplotype (tools/c/semig.c; reverse complement included) and is dropped if it is
more than MAX_ERR of its length from every haplotype. A diploid genotype (h1, h2) scores
    sum over reads of log( 0.5 exp(-lam d_r1) + 0.5 exp(-lam d_r2) )
and the best pair is compared with HG002's truth haplotypes (edit distance, best pairing).

    python3 tools/panel_genotype.py [--regions work/iterate/testset.tsv] [--lam 3] [--out work/panelgt]

Needs VNTR_REGIONS (region packages) and VNTR_WORK (the call_local.py work dir holding the cached window
read queries, e.g. work/iterate/w) in the environment.
"""
import argparse
import csv
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import call_local  # noqa: E402
import score_haplotypes as sh  # noqa: E402

K = 21
MAX_ERR = 0.08
DEPTH_W = float(os.environ.get('PANELGT_DEPTH_W', '1'))
COMP = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(COMP)[::-1]


def read_fastq(path):
    out = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().strip().upper()
            f.readline(); f.readline()
            out.append(s)
    return out


def kmers(s):
    return {s[i:i + K] for i in range(len(s) - K + 1)}


def region_reads(rid):
    wd = call_local.fetch_window_reads(rid, 200000)
    wd = wd if isinstance(wd, str) and os.path.isdir(wd) else os.path.join(
        call_local.STAGE3, 'hybrid200k', 'reads', rid)
    reads = []
    for f in ('reads_1.fq', 'reads_2.fq', 'reads_se.fq'):
        reads += read_fastq(os.path.join(wd, f))
    return reads


def genotype(rid, lam, out_dir, semig):
    rd = os.path.join(os.environ['VNTR_REGIONS'], rid)
    reg = sh.load_region(rd)
    fa = sh.read_fasta(os.path.join(rd, 'hap32.fa'))
    seqs = list(fa.values()) if isinstance(fa, dict) else [s for _, s in fa]
    haps = sorted(set(s.upper() for s in seqs))
    t1, t2 = reg['truth']
    ed = sh._load_edit_distance()
    # variable k-mers: in some but not all distinct haplotypes
    ks = [kmers(h) for h in haps]
    allk = set.union(*ks)
    common = set.intersection(*ks) if len(ks) > 1 else allk
    var = allk - common
    varall = var | {rc(k) for k in var}
    FL = 2000
    left = call_local.fetch_ref(reg['contig'], max(0, reg['a1'] - 1 - FL), reg['a1'] - 1)
    right = call_local.fetch_ref(reg['contig'], reg['b1'], reg['b1'] + FL)
    regk = allk | {rc(k) for k in allk}
    flk = kmers(left) | kmers(right)
    flk = flk | {rc(k) for k in flk}
    allr = region_reads(rid)
    reads = [r for r in allr if any(r[i:i + K] in regk for i in range(0, len(r) - K + 1, 3))]
    nflank = sum(1 for r in allr if all(r[i:i + K] in flk for i in (0, len(r) - K)) and
                 not any(r[i:i + K] in regk for i in range(0, len(r) - K + 1, 3)))
    wd = os.path.join(out_dir, rid)
    os.makedirs(wd, exist_ok=True)
    with open(os.path.join(wd, 'haps.fa'), 'w') as f:
        for i, h in enumerate(haps):
            f.write('>h%d\n%s\n' % (i, h))
    with open(os.path.join(wd, 'reads.fa'), 'w') as f:
        for i, r in enumerate(reads):
            f.write('>r%d\n%s\n' % (i, r))
    p = subprocess.run([semig, os.path.join(wd, 'haps.fa'), os.path.join(wd, 'reads.fa')],
                       capture_output=True, text=True, check=True)
    D = []
    for line, r in zip(p.stdout.splitlines(), reads):
        d = [int(x) for x in line.split('\t')[1:]]
        if min(d) <= MAX_ERR * len(r):
            D.append(d)
    n = len(haps)
    R = 150.0
    kappa = nflank / (4.0 * (FL - R + 1))        # read starts per bp per haplotype, from the flanks
    N = len(D)
    L = [len(h) for h in haps]
    best, bp = -math.inf, (0, 0)
    for i in range(n):
        for j in range(i, n):
            ll = 0.0
            for d in D:
                a, b = -lam * d[i], -lam * d[j]
                m = max(a, b)
                ll += m + math.log(0.5 * math.exp(a - m) + 0.5 * math.exp(b - m))
            if DEPTH_W > 0 and kappa > 0:
                mu = kappa * (L[i] + L[j] + 2 * (R - 1))
                ll += DEPTH_W * (N * math.log(mu) - mu - math.lgamma(N + 1))
            if ll > best:
                best, bp = ll, (i, j)
    h1, h2 = haps[bp[0]], haps[bp[1]]
    called = min(ed(h1, t1) + ed(h2, t2), ed(h1, t2) + ed(h2, t1))
    d1 = [ed(h, t1) for h in haps]
    d2 = [ed(h, t2) for h in haps]
    ceiling = min(min(d1[i] + d2[j], d1[j] + d2[i]) for i in range(n) for j in range(i, n))
    return {'region_id': rid, 'haps': n, 'reads_kept': len(D), 'kappa': round(kappa, 4), 'n_flank': nflank,
            'called_len': '%d/%d' % (len(h1), len(h2)), 'truth_len': '%d/%d' % (len(t1), len(t2)),
            'panel_call_ed': called, 'ceiling_ed': ceiling}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--regions', default='work/iterate/testset.tsv')
    ap.add_argument('--lam', type=float, default=3.0)
    ap.add_argument('--out', default='work/panelgt')
    a = ap.parse_args()
    semig = os.path.join(a.out, 'semig')
    rows = [r for r in csv.DictReader((l for l in open(a.regions) if not l.startswith('#')), delimiter='\t')]
    res = []
    for r in rows:
        try:
            x = genotype(r['region_id'], a.lam, a.out, semig)
            x['group'] = r['group']
        except Exception as e:  # noqa: BLE001
            x = {'region_id': r['region_id'], 'group': r['group'], 'error': repr(e)[:200]}
        res.append(x)
        print('\t'.join('%s=%s' % kv for kv in x.items()), flush=True)
    keys = ['region_id', 'group', 'haps', 'reads_kept', 'kappa', 'n_flank', 'called_len', 'truth_len', 'panel_call_ed', 'ceiling_ed', 'error']
    with open(os.path.join(a.out, 'panelgt_lam%g_dw%g.tsv' % (a.lam, DEPTH_W)), 'w') as f:
        f.write('\t'.join(keys) + '\n')
        for x in res:
            f.write('\t'.join(str(x.get(k, '')) for k in keys) + '\n')


if __name__ == '__main__':
    main()
