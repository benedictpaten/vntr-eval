#!/usr/bin/env python3
"""hardset.py -- the hard tandem-repeat benchmark: the regions whose realigned alignment agrees worst
with the pairwise optima, as a committed region list plus each region's alleles as FASTA.

    python3 tools/hardset.py                       # chr20, pair F1 < 0.8 -> hardset/chr20.tsv
    python3 tools/hardset.py --max-f1 0.9 -o hardset/chr20_lt0.9.tsv

Inputs (all under the contig's pgrealign run directory, work/full/<contig>/ by default):
  qc/latest.tsv              pgrealign qc of the latest realigned graph, with --graph (so each row
                             has the original graph's counts too); a region is hard when its pair F1,
                             2 pair_tp / (pair_opt + pair_msa), is below --max-f1
  <contig>.regions.tsv.gz    the region table (coordinates, targets, classes)
  extract/                   extract's packages (alleles) and occurrences (which allele is CHM13's)
  msas.banded.jsonl.gz       the realignment qc/latest.tsv scored, for its per-region run time
  catalogs/chm13v2.0_tr.bed.gz   period and motif of each target repeat
  --switches                 a switch table from vg-call-eval's offpanel/switches/ (columns start, end,
                             scale, pair); the local mis-phased stretches overlapping a region are counted

Outputs:
  hardset/<contig>.tsv                  the region list (columns in hardset/README.md)
  work/hardset/<contig>/fasta/<id>.fa   every distinct allele of the region as pgrealign extracted
                                        it, in its forward orientation: '>a<k> weight=<w>', where the
                                        weight is the number of panel haplotypes spelling it
"""
import argparse
import csv
import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402

PGREALIGN_DIR = config._env('VNTR_PGREALIGN', os.path.join(config.HOME, 'PycharmProjects/pgrealign'))
try:
    import pgrealign  # noqa: F401
except ImportError:
    sys.path.insert(0, PGREALIGN_DIR)
from pgrealign.qc import reference_alleles  # noqa: E402

SWITCHES = os.path.join(config.EVAL_DIR, 'work/offpanel/switches/switches_chr20_ont_1c8f89b7a.tsv')
CATALOG = os.path.join(config.REPO, 'catalogs/chm13v2.0_tr.bed.gz')
FIELDS = ['id', 'contig', 'ref_start', 'ref_end', 'span', 'L', 'R', 'targets', 'classes', 'period', 'motif',
          'alleles', 'longest', 'ref_allele', 'pair_f1_graph', 'pair_f1_realigned', 'ref_f1_graph',
          'ref_f1_realigned', 'realign_seconds', 'flip_stretches']


def f1(tp, opt, msa):
    tp, opt, msa = int(tp or 0), int(opt or 0), int(msa or 0)
    return round(2 * tp / (opt + msa), 4) if opt + msa else None


def hard_ids(qc_path, max_f1):
    """{id: qc row} of the regions whose realigned pair F1 is below max_f1."""
    out = {}
    with open(qc_path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            v = f1(r['pair_tp'], r['pair_opt'], r['pair_msa'])
            if v is not None and v < max_f1:
                out[r['id']] = r
    return out


def local_stretches(path):
    """[(start, end)] of the local mis-phased stretches: each pair of local switches, start to end."""
    if not path or not os.path.exists(path):
        return None
    span = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['scale'] == 'local':
                s, e = int(r['start']), int(r['end'])
                a, b = span.get(r['pair'], (s, e))
                span[r['pair']] = (min(a, s), max(b, e))
    return sorted(span.values())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--contig', default='chr20')
    ap.add_argument('--run-dir', help='pgrealign run directory (default work/full/<contig>)')
    ap.add_argument('--max-f1', type=float, default=0.8)
    ap.add_argument('--switches', default=SWITCHES)
    ap.add_argument('-o', '--out', help='default hardset/<contig>.tsv')
    ap.add_argument('--fasta-dir', help='default work/hardset/<contig>/fasta')
    a = ap.parse_args(argv)
    run = a.run_dir or os.path.join(config.WORK_DIR, 'full', a.contig)
    out = a.out or os.path.join(config.REPO, 'hardset', a.contig + '.tsv')
    fasta_dir = a.fasta_dir or os.path.join(config.WORK_DIR, 'hardset', a.contig, 'fasta')

    qc = hard_ids(os.path.join(run, 'qc/latest.tsv'), a.max_f1)
    print('%d regions with pair F1 < %g' % (len(qc), a.max_f1), file=sys.stderr)
    regions = {}
    with gzip.open(os.path.join(run, a.contig + '.regions.tsv.gz'), 'rt') as f:
        head = None
        for line in f:
            x = line.rstrip('\n').split('\t')
            if head is None:
                head = [h.lstrip('#') for h in x]
            elif x[0] in qc:
                regions[x[0]] = dict(zip(head, x))
    catalog = {}
    targets = {t for r in regions.values() for t in r['targets'].split(',')}
    with gzip.open(CATALOG, 'rt') as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            if len(x) > 6 and x[3] in targets:
                catalog[x[3]] = (x[5], x[6], int(x[2]) - int(x[1]))
    seconds = {}
    with gzip.open(os.path.join(run, 'msas.banded.jsonl.gz'), 'rt') as f:
        for line in f:
            i = line.find('"id":"')
            if i >= 0 and line[i + 6:line.find('"', i + 6)] in qc:
                m = json.loads(line)
                seconds[m['id']] = m.get('seconds')
    ref = reference_alleles(os.path.join(run, 'extract/occurrences.tsv.gz'), 'CHM13')
    stretches = local_stretches(a.switches)
    os.makedirs(fasta_dir, exist_ok=True)
    rows = []
    with gzip.open(os.path.join(run, 'extract/packages.jsonl.gz'), 'rt') as f:
        for line in f:
            p = json.loads(line)
            if p['id'] not in qc:
                continue
            with open(os.path.join(fasta_dir, p['id'] + '.fa.tmp'), 'w') as o:
                for al in p['alleles']:
                    o.write('>%s weight=%d\n%s\n' % (al['id'], al['weight'], al['seq'].upper()))
            os.replace(os.path.join(fasta_dir, p['id'] + '.fa.tmp'), os.path.join(fasta_dir, p['id'] + '.fa'))
            g, q = regions[p['id']], qc[p['id']]
            t0 = max(g['targets'].split(','), key=lambda t: catalog.get(t, ('', '', 0))[2])
            s, e = int(g['ref_start']), int(g['ref_end'])
            rows.append({
                'id': p['id'], 'contig': g['contig'], 'ref_start': s, 'ref_end': e, 'span': g['span'],
                'L': g['L'], 'R': g['R'], 'targets': g['targets'], 'classes': g['classes'],
                'period': catalog.get(t0, ('', ''))[0], 'motif': catalog.get(t0, ('', ''))[1],  # the longest target's
                'alleles': len(p['alleles']), 'longest': max(len(al['seq']) for al in p['alleles']),
                'ref_allele': ref.get(p['id'], ''),
                'pair_f1_graph': f1(q['graph_pair_tp'], q['pair_opt'], q['graph_pair_msa']),
                'pair_f1_realigned': f1(q['pair_tp'], q['pair_opt'], q['pair_msa']),
                'ref_f1_graph': f1(q['graph_ref_tp'], q['ref_opt'], q['graph_ref_msa']),
                'ref_f1_realigned': f1(q['ref_tp'], q['ref_opt'], q['ref_msa']),
                'realign_seconds': seconds.get(p['id'], ''),
                'flip_stretches': '' if stretches is None else sum(1 for a0, b0 in stretches if a0 <= e and b0 >= s),
            })
    rows.sort(key=lambda r: r['pair_f1_realigned'])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as o:
        o.write('\t'.join(FIELDS) + '\n')
        for r in rows:
            o.write('\t'.join('' if r[k] is None else str(r[k]) for k in FIELDS) + '\n')
    print('wrote %s (%d regions) and %d FASTAs in %s' % (out, len(rows), len(rows), fasta_dir), file=sys.stderr)


if __name__ == '__main__':
    main()
