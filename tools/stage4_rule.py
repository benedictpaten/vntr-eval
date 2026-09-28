#!/usr/bin/env python3
"""stage4_rule.py -- Stage 4 detection rule: fit on the packaged non-chr20 regions, apply to chr20.

    python3 tools/stage4_rule.py fit   [--exclude-contig chr20] [--out work/stage4/rule.json]
    python3 tools/stage4_rule.py apply --rule work/stage4/rule.json --scan work/stage4/scan.tsv \
        --results work/stage4/results/mc --out work/stage4/hotspots.tsv

fit: labels from regions/regions.tsv (hotspot = hotspot_vntr + hotspot_other, control = the three
control strata), features from results/mc/<id>.json (Stage 0 of the MC subgraph). Candidate rules
(fixed before any chr20 region is looked at):
  - one feature >= T, for each of all_cost_over_opt, kmer_frac_extra, kmer_extra_per_kb;
  - all_cost_over_opt >= T1 AND/OR kmer_frac_extra >= T2.
T is the midpoint between adjacent observed values; the rule with the highest balanced accuracy
wins, except that a two-feature rule must beat the best one-feature rule by >= 0.02 to be chosen.
AUC (Mann-Whitney, ties 1/2) is reported per feature. A region with a missing feature is 'not
detected'. apply: the same rule on each chr20 region's Stage 0 JSON -> hotspots.tsv (every
packaged region, with a 'detected' column). Standard library only.
"""
import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FEATS = {'all_cost_over_opt': ('alignment', 'all_cost_over_opt'),
         'kmer_frac_extra': ('redundancy', 'kmer_frac_extra'),
         'kmer_extra_per_kb': ('redundancy', 'kmer_extra_per_kb')}
HOT = {'hotspot_vntr', 'hotspot_other'}
CTRL = {'control_vntr_matched', 'control_vntr_correct', 'control_nontr_sv'}


def feats(path):
    if not os.path.exists(path):
        return None
    d = json.load(open(path))
    out = {}
    for k, (sec, key) in FEATS.items():
        v = d.get(sec, {}).get(key)
        out[k] = float(v) if isinstance(v, (int, float)) else None
    return out


def auc(pos, neg):
    s = 0.0
    for p in pos:
        for n in neg:
            s += 1.0 if p > n else 0.5 if p == n else 0.0
    return s / (len(pos) * len(neg))


def thresholds(vals):
    u = sorted(set(vals))
    return [u[0] - 1e-9] + [(a + b) / 2 for a, b in zip(u, u[1:])] + [u[-1] + 1e-9]


def ba(pred, lab):
    tp = sum(1 for p, y in zip(pred, lab) if p and y)
    tn = sum(1 for p, y in zip(pred, lab) if not p and not y)
    npos = sum(lab)
    nneg = len(lab) - npos
    return 0.5 * (tp / npos + tn / nneg), tp, npos - tp, nneg - tn, tn


def rule_pred(rule, f):
    def ge(k, t):
        return f.get(k) is not None and f[k] >= t
    if rule['kind'] == 'single':
        return ge(rule['feature'], rule['T'])
    a, b = ge('all_cost_over_opt', rule['T1']), ge('kmer_frac_extra', rule['T2'])
    return (a and b) if rule['kind'] == 'and' else (a or b)


def cmd_fit(a):
    rows = []
    with open(os.path.join(ROOT, 'regions', 'regions.tsv')) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['contig'] == a.exclude_contig or r['stratum'] not in HOT | CTRL:
                continue
            f = feats(os.path.join(ROOT, 'results', 'mc', r['region_id'] + '.json'))
            if f is None:
                continue
            rows.append((r['region_id'], r['stratum'], r['stratum'] in HOT, f))
    lab = [y for _, _, y, _ in rows]
    rep = {'n': len(rows), 'n_hotspot': sum(lab), 'n_control': len(lab) - sum(lab), 'features': {}}
    best = None
    for k in FEATS:
        pos = [f[k] for _, _, y, f in rows if y and f[k] is not None]
        neg = [f[k] for _, _, y, f in rows if not y and f[k] is not None]
        fb = None
        for t in thresholds(pos + neg):
            rule = {'kind': 'single', 'feature': k, 'T': round(t, 6)}
            s = ba([rule_pred(rule, f) for *_, f in rows], lab)
            if fb is None or s[0] > fb[0][0]:
                fb = (s, rule)
        rep['features'][k] = {'auc': round(auc(pos, neg), 4), 'n_missing': len(rows) - len(pos) - len(neg),
                              'best_rule': fb[1], 'best_ba': round(fb[0][0], 4),
                              'tp_fn_fp_tn': fb[0][1:]}
        if best is None or fb[0][0] > best[0][0]:
            best = fb
    c1 = thresholds([f['all_cost_over_opt'] for *_, f in rows if f['all_cost_over_opt'] is not None])
    c2 = thresholds([f['kmer_frac_extra'] for *_, f in rows if f['kmer_frac_extra'] is not None])
    combo = None
    for kind in ('and', 'or'):
        for t1 in c1:
            for t2 in c2:
                rule = {'kind': kind, 'T1': round(t1, 6), 'T2': round(t2, 6)}
                s = ba([rule_pred(rule, f) for *_, f in rows], lab)
                if combo is None or s[0] > combo[0][0]:
                    combo = (s, rule)
    rep['best_single'] = {'rule': best[1], 'ba': round(best[0][0], 4), 'tp_fn_fp_tn': best[0][1:]}
    rep['best_combo'] = {'rule': combo[1], 'ba': round(combo[0][0], 4), 'tp_fn_fp_tn': combo[0][1:]}
    chosen = combo if combo[0][0] >= best[0][0] + 0.02 else best
    rep['rule'] = chosen[1]
    rep['rule_ba'] = round(chosen[0][0], 4)
    rep['rule_tp_fn_fp_tn'] = chosen[0][1:]
    rep['per_stratum_detected'] = {}
    for st in sorted(HOT | CTRL):
        sel = [f for _, s, _, f in rows if s == st]
        rep['per_stratum_detected'][st] = '%d/%d' % (sum(rule_pred(chosen[1], f) for f in sel), len(sel))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(rep, open(a.out, 'w'), indent=1)
    print(json.dumps(rep, indent=1))


def cmd_apply(a):
    rule = json.load(open(a.rule))['rule']
    scan = {}
    with open(a.scan) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['prefilter'] == '1':
                scan[r['region_id']] = r
    n = nd = 0
    with open(a.out, 'w') as fo:
        cols = ['region_id', 'contig', 'start', 'end', 'tr_class', 'period', 'chm13_len', 'max_diff',
                'n_walks', 'status'] + list(FEATS) + ['detected']
        fo.write('\t'.join(cols) + '\n')
        for rid, r in scan.items():
            f = feats(os.path.join(a.results, rid + '.json'))
            st = 'no_stage0' if f is None else 'ok'
            det = int(f is not None and rule_pred(rule, f))
            n += 1
            nd += det
            fo.write('\t'.join([rid, r['contig'], r['start'], r['end'], r['tr_class'], r['period'],
                                r['chm13_len'], r['max_diff'], r['n_walks'], st] +
                               ['' if f is None or f[k] is None else '%g' % f[k] for k in FEATS] +
                               [str(det)]) + '\n')
    print('%d prefiltered regions, %d detected (rule %s)' % (n, nd, json.dumps(rule)), file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    f = sub.add_parser('fit')
    f.add_argument('--exclude-contig', default='chr20')
    f.add_argument('--out', default=os.path.join(ROOT, 'work', 'stage4', 'rule.json'))
    p = sub.add_parser('apply')
    p.add_argument('--rule', required=True)
    p.add_argument('--scan', required=True)
    p.add_argument('--results', required=True)
    p.add_argument('--out', required=True)
    a = ap.parse_args(argv)
    {'fit': cmd_fit, 'apply': cmd_apply}[a.cmd](a)


if __name__ == '__main__':
    main()
