#!/usr/bin/env python3
"""Compare vg call's called haplotypes (hapscore_chr20.py) with whole-allele genotyping of each repeat
(atomic_chr20.py) on chr20: summed edit distance to HG002's truth haplotypes (best pairing), the share of
vg call's error that comes from haplotypes no hap32 sequence spells (stitched across sites), and a
paired bootstrap over regions of the per-region difference."""
import csv, json, os, random, sys
R = '/Users/benedictpaten/PycharmProjects/vntr-eval'
sys.path.insert(0, R + '/tools')
import msa_graph  # noqa: E402

ARMS = [a for a in ['mc', 'abpoa_h32', 'st_chm13', 'st_medoid', 'rep_site', 'rep_desc', 'rep_link']
        if os.path.isdir(R + '/work/iterate/hapscore/' + a)]
atomic = {r['region_id']: r for r in csv.DictReader(open(R + '/work/iterate/atomic/atomic_chr20_dw%s.tsv'
                                                          % os.environ.get('PANELGT_DEPTH_W', '1')), delimiter='\t')}


def panel_seqs(rid):
    return {s.upper() for _, s in msa_graph.read_fasta(R + '/work/stage4/regions/%s/hap32.fa' % rid)}


rows = {}
for rid, a in atomic.items():
    if a.get('error') or a.get('panel_call_ed') in ('', None):
        continue
    rec = {'atomic': int(a['panel_call_ed']), 'ceiling': int(a['ceiling_ed'])}
    H = None
    ok = True
    for arm in ARMS:
        js = R + '/work/iterate/hapscore/%s/%s.json' % (arm, rid)
        fa = R + '/work/iterate/hapscore/%s/%s.fa' % (arm, rid)
        if not os.path.exists(js) or not os.path.exists(fa):
            ok = False
            break
        d = json.load(open(js))
        if d.get('status') != 'ok':
            ok = False
            break
        if H is None:
            H = panel_seqs(rid)
        haps = [s.upper() for _, s in msa_graph.read_fasta(fa)]
        # per called haplotype: its error under the scorer's best pairing, and whether a hap32 sequence spells it
        pairing = d.get('pairing', 'slot1-h1,slot2-h2')
        e = [d.get('ed_h1', 0), d.get('ed_h2', 0)]
        slot_err = {}
        for part, err in zip(pairing.split(','), e):
            slot_err[part.split('-')[0]] = err
        off = [haps[i] not in H for i in range(len(haps))]
        rec[arm] = d['ed']
        rec[arm + '_off'] = sum(off)
        rec[arm + '_off_ed'] = sum(slot_err.get('slot%d' % (i + 1), 0) for i in range(len(haps)) if off[i])
    if ok:
        rows[rid] = rec

print('regions scored by every method: %d' % len(rows))
tot = {k: sum(r[k] for r in rows.values()) for k in ['ceiling', 'atomic'] + ARMS}
print('\nsummed edit distance to the truth haplotypes (lower is better), regions exact')
print('  %-34s %9d  %d exact' % ('ceiling: best pair of hap32 sequences', tot['ceiling'], sum(1 for r in rows.values() if r['ceiling'] == 0)))
print('  %-34s %9d  %d exact' % ('atomic: one site of whole alleles', tot['atomic'], sum(1 for r in rows.values() if r['atomic'] == 0)))
for arm in ARMS:
    off = sum(r[arm + '_off'] for r in rows.values())
    offed = sum(r[arm + '_off_ed'] for r in rows.values())
    print('  %-34s %9d  %d exact;  %d of %d called haplotypes off-panel, carrying %d (%.0f%%) of the error' % (
        'vg call, ' + arm, tot[arm], sum(1 for r in rows.values() if r[arm] == 0), off, 2 * len(rows), offed,
        100.0 * offed / tot[arm] if tot[arm] else 0))

random.seed(1)
def boot(a, b, excl=()):
    d = [r[a] - r[b] for k, r in rows.items() if k not in excl]
    n = len(d)
    bs = sorted(sum(random.choice(d) for _ in range(n)) for _ in range(2000))
    return sum(d), bs[50], bs[1949], sum(x < 0 for x in d), sum(x > 0 for x in d)
print('\nper region, paired bootstrap of the edit-distance difference (negative: first is better)')
for a, b in [('atomic', 'st_medoid'), ('atomic', 'st_chm13'), ('atomic', 'abpoa_h32'), ('atomic', 'mc'),
             ('st_medoid', 'abpoa_h32'), ('st_medoid', 'mc'), ('rep_site', 'st_medoid'), ('rep_desc', 'st_medoid'), ('rep_link', 'st_medoid'), ('atomic', 'rep_site'), ('rep_desc', 'rep_site')]:
    if a not in ARMS + ['atomic'] or b not in ARMS + ['atomic']:
        continue
    s, lo, hi, bt, wr = boot(a, b)
    print('  %-10s - %-10s %+8d  [%+d, %+d]  better/worse %d/%d' % (a, b, s, lo, hi, bt, wr))
json.dump(rows, open(R + '/work/iterate/atomic/compare_chr20.json', 'w'))
