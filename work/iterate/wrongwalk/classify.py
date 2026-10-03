#!/usr/bin/env python3
"""Covered chr20 repeat regions (every overlapping top-level snarl inside the region, so --repeat-sites
offers only whole panel walks): what kind of wrong walk does vg pick?

For each region: the called pair (hapscore/rep_site/<id>.fa), the truth pair (truth.fa), the best panel
pair (the ceiling), and the whole-allele prototype's chosen lengths (atomic TSV). Per called haplotype,
the panel allele it spells and its panel count; per truth haplotype, its closest panel allele.

    python3 work/iterate/wrongwalk/classify.py   -> work/iterate/wrongwalk/classify.tsv + summary
"""
import collections, csv, itertools, json, os, sys
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
sys.path.insert(0, R + '/tools')
import msa_graph                                   # noqa: E402
from score_haplotypes import edit_distance as ED   # noqa: E402

cls = json.load(open(R + '/work/iterate/atomic/repeat_cover_class_chr20.json'))
rows = json.load(open(R + '/work/iterate/atomic/compare_chr20.json'))
atomic = {r['region_id']: r for r in csv.DictReader(open(R + '/work/iterate/atomic/atomic_chr20_dw1.tsv'), delimiter='\t')}


def fa(p):
    return [(n, s.upper()) for n, s in msa_graph.read_fasta(p)]


out = []
for rid in sorted(k for k, c in cls.items() if c == 'covered' and k in rows):
    d = R + '/work/stage4/regions/' + rid
    truth = [s for _, s in fa(d + '/truth.fa')]
    panel = fa(d + '/hap32.fa')
    alle = collections.OrderedDict()
    for n, s in panel:
        alle.setdefault(s, []).append(n)
    seqs = list(alle)
    called = [s for _, s in fa(R + '/work/iterate/hapscore/rep_site/%s.fa' % rid)]
    if len(truth) != 2 or len(called) != 2:
        continue
    dt = [[ED(a, t) for t in truth] for a in seqs]          # panel allele x truth hap
    best = min(((dt[i][0] + dt[j][1], i, j) for i in range(len(seqs)) for j in range(len(seqs))))
    ci = [seqs.index(c) if c in alle else -1 for c in called]
    ed_called = min(ED(called[0], truth[0]) + ED(called[1], truth[1]), ED(called[0], truth[1]) + ED(called[1], truth[0]))
    # pair called slots to truth haps the better way, and to the ceiling's alleles for those haps
    a = ED(called[0], truth[0]) + ED(called[1], truth[1])
    b = ED(called[0], truth[1]) + ED(called[1], truth[0])
    order = (0, 1) if a <= b else (1, 0)                   # called slot k is paired with truth order[k]
    ceil_for = {0: best[1], 1: best[2]}
    slots = []
    for k in range(2):
        t = order[k]
        cs, bs = called[k], seqs[ceil_for[t]]
        slots.append({'right': cs == bs, 'dlen': len(cs) - len(bs), 'called_n': len(alle.get(cs, [])),
                      'best_n': len(alle[bs]), 'ed': ED(cs, truth[t]), 'ceil_ed': dt[ceil_for[t]][t]})
    at = atomic.get(rid, {})
    out.append({'region': rid, 'n_alleles': len(seqs), 'truth_len': '/'.join(str(len(t)) for t in truth),
                'called_len': '/'.join(str(len(c)) for c in called), 'atomic_len': at.get('called_len', ''),
                'ceil_len': '%d/%d' % (len(seqs[best[1]]), len(seqs[best[2]])),
                'ed': ed_called, 'ceil': best[0], 'atomic': rows[rid]['atomic'], 'excess': ed_called - best[0],
                'n_wrong': sum(not s['right'] for s in slots),
                'max_abs_dlen': max(abs(s['dlen']) for s in slots if not s['right']) if any(not s['right'] for s in slots) else 0,
                'truth_het': int(best[1] != best[2]), 'called_het': int(called[0] != called[1]),
                'wrong_called_n': ','.join(str(s['called_n']) for s in slots if not s['right']),
                'wrong_best_n': ','.join(str(s['best_n']) for s in slots if not s['right']),
                'dlen': ','.join(str(s['dlen']) for s in slots if not s['right'])})
keys = list(out[0])
with open(R + '/work/iterate/wrongwalk/classify.tsv', 'w') as f:
    f.write('\t'.join(keys) + '\n')
    for o in out:
        f.write('\t'.join(str(o[k]) for k in keys) + '\n')

tot = sum(o['ed'] for o in out); exc = sum(o['excess'] for o in out)
print('covered regions scored: %d; edit distance %d = ceiling %d + excess %d' % (len(out), tot, tot - exc, exc))
w = [o for o in out if o['n_wrong']]
print('regions with a wrong walk: %d (excess %d); right walks but not exact (ceiling error only): %d' % (
    len(w), sum(o['excess'] for o in w), sum(1 for o in out if not o['n_wrong'] and o['ed'])))
w.sort(key=lambda o: -o['excess'])
cum = 0
for i, o in enumerate(w):
    cum += o['excess']
    if i + 1 in (5, 10, 20, 40) or i + 1 == len(w):
        print('  top %d regions hold %d (%.0f%%) of the wrong-walk excess' % (i + 1, cum, 100 * cum / max(1, sum(x['excess'] for x in w))))
def bucket(o):
    m = o['max_abs_dlen']
    return 'same length' if m == 0 else ('|dlen| 1-9' if m < 10 else ('|dlen| 10-99' if m < 100 else '|dlen| >=100'))
print('\nwrong-walk regions by the largest length error of a wrong slot:')
for b in ['same length', '|dlen| 1-9', '|dlen| 10-99', '|dlen| >=100']:
    xs = [o for o in w if bucket(o) == b]
    print('  %-13s %3d regions, excess %6d' % (b, len(xs), sum(o['excess'] for o in xs)))
print('\nby zygosity (truth by ceiling pair / called):')
z = collections.Counter(); ze = collections.Counter()
for o in w:
    k = ('het' if o['truth_het'] else 'hom') + '->' + ('het' if o['called_het'] else 'hom')
    z[k] += 1; ze[k] += o['excess']
for k in sorted(z):
    print('  %-9s %3d regions, excess %6d' % (k, z[k], ze[k]))
print('  one slot wrong: %d regions (excess %d); both: %d (excess %d)' % (
    sum(o['n_wrong'] == 1 for o in w), sum(o['excess'] for o in w if o['n_wrong'] == 1),
    sum(o['n_wrong'] == 2 for o in w), sum(o['excess'] for o in w if o['n_wrong'] == 2)))
lens = [int(x) for o in w for x in o['dlen'].split(',') if x]
print('  wrong slots: %d, called longer than the right allele in %d, shorter in %d' % (len(lens), sum(x > 0 for x in lens), sum(x < 0 for x in lens)))
print('\nprototype on the wrong-walk regions: %d vs vg %d (ceiling %d); prototype better in %d, worse in %d' % (
    sum(o['atomic'] for o in w), sum(o['ed'] for o in w), sum(o['ceil'] for o in w),
    sum(o['atomic'] < o['ed'] for o in w), sum(o['atomic'] > o['ed'] for o in w)))
print('\ntop 15:')
for o in w[:15]:
    print('  %(region)s alleles %(n_alleles)2d truth %(truth_len)s ceil %(ceil_len)s called %(called_len)s proto %(atomic_len)s | ed %(ed)d ceil %(ceil)d proto %(atomic)d | dlen %(dlen)s wrong-called n %(wrong_called_n)s right n %(wrong_best_n)s' % o)
