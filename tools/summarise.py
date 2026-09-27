#!/usr/bin/env python3
"""summarise.py -- compare candidate graphs across regions: results/<METHOD>/<region_id>.json
(written by tools/evaluate.py) into per-stratum tables.

    python3 tools/summarise.py [--results results] [--baseline mc] [--methods mc,mafft_linsi,...]
                               [--regions regions/regions.tsv] [--metrics headline|all]
                               [--md results/summary.md] [--tsv results/summary.tsv]
                               [--per-region results/per_region.tsv]

For every stratum (from each JSON's `stratum`, else regions.tsv) and every metric:
  - per method: n regions, median [q1-q3], mean;
  - paired against the baseline method (default 'mc') on the regions both have: the median and
    IQR of the per-region difference (method - baseline), how many regions got better / worse /
    stayed equal (the metric's direction is fixed in METRICS below), and a two-sided Wilcoxon
    signed-rank p-value (exact permutation distribution of the observed ranks, ties averaged,
    zero differences dropped; normal approximation above 60 pairs);
  - truvari: pooled F1 over the stratum (TP/FP/FN summed over regions), raw and refined, besides
    the per-region medians.
Only regions whose graph was valid (status 'ok') are used; invalid graphs are counted per method.
Markdown goes to stdout unless --md is given. Every metric is defined in tools/METRICS.md.
"""
import argparse
import collections
import csv
import glob
import json
import math
import os
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402

# (key path in the evaluate.py JSON, short label, direction: -1 lower is better, +1 higher is
#  better, 0 no preferred direction), headline flag
METRICS = [
    ('size.nodes_per_kb', 'nodes / CHM13 kb', -1, True),
    ('size.node_len_mean', 'mean node bp', +1, False),
    ('size.node_len_median', 'median node bp', +1, False),
    ('size.node_frac_1bp', 'frac 1-bp nodes', -1, False),
    ('size.edges_per_kb', 'edges / CHM13 kb', -1, False),
    ('size.bubbles', 'top-level bubbles', 0, False),
    ('size.largest_bubble_alleles', 'alleles, largest bubble', 0, False),
    ('size.distinct_walks', 'distinct walks', 0, False),
    ('alignment.all_cost_over_opt', 'cost/opt, all pairs', -1, True),
    ('alignment.all_pair_ratio_median', 'cost/opt, median pair', -1, False),
    ('alignment.all_excess_per_kb', 'excess edits / kb', -1, True),
    ('alignment.all_frac_pairs_optimal', 'frac pairs optimal', +1, False),
    ('alignment.ref_cost_over_opt', 'cost/opt vs CHM13', -1, False),
    ('alignment.affine_all_cost_over_opt', 'affine cost/opt, all', -1, True),
    ('alignment.affine_ref_cost_over_opt', 'affine cost/opt vs CHM13', -1, False),
    ('alignment.affine_all_excess_per_kb', 'affine excess / kb', -1, False),
    ('alignment.all_U_per_kb', 'unaligned homology bp / kb', -1, True),
    ('inflation.sv_pieces_median', 'SV pieces / path (median)', -1, True),
    ('inflation.sv_pieces_mean', 'SV pieces / path (mean)', -1, False),
    ('inflation.anchored_over_free', 'anchored/free vs CHM13', -1, False),
    ('inflation.indel_bp_ratio', 'indel bp / net length change', -1, True),
    ('inflation.indel_bp_excess_per_kb', 'excess indel bp / kb', -1, False),
    ('redundancy.kmer_frac_extra', 'k-mer extra positions (frac)', -1, True),
    ('redundancy.kmer_extra_per_kb', 'k-mer extra positions / kb', -1, False),
    ('truth.d_graph_sum', 'truth edits to best graph path (h1+h2)', -1, True),
    ('truth.h1_d_graph', 'truth h1 edits to graph', -1, False),
    ('truth.h2_d_graph', 'truth h2 edits to graph', -1, False),
    ('truth.d_panel_sum', 'truth edits to closest panel path (h1+h2)', 0, False),
    ('truth.truvari.raw_f1', 'truth-by-graph truvari F1, raw', +1, True),
    ('truth.truvari.refined_f1', 'truth-by-graph truvari F1, refined', +1, True),
    ('truth.truvari.phab_f1', 'truth-by-graph truvari F1, whole-span phab', +1, False),
    ('truth.truvari.records_sv', 'truth-by-graph SV records', 0, False),
    ('reads.all_unique', 'reads: unique placement', +1, False),
    ('reads.all_redundant_only', 'reads: redundant (cross-walk) only', -1, True),
    ('reads.all_tandem', 'reads: tandem', 0, False),
    ('reads.all_ratio_mean', 'reads: placements / max per path', -1, True),
    ('reads.mqlt5_redundant_only', 'reads MAPQ<5: redundant only', -1, False),
    ('reads.mqlt5_ratio_mean', 'reads MAPQ<5: placements / max per path', -1, False),
]
STRATA_ORDER = ['hotspot_vntr', 'control_vntr_matched', 'control_vntr_correct', 'hotspot_other',
                'control_nontr_sv']


def get(d, path):
    for k in path.split('.'):
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    if isinstance(d, bool):
        return int(d)
    return d if isinstance(d, (int, float)) else None


def quantile(s, p):
    k = (len(s) - 1) * p
    lo = int(math.floor(k))
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def qstats(xs):
    s = sorted(xs)
    if not s:
        return None
    return {'n': len(s), 'median': quantile(s, 0.5), 'q1': quantile(s, 0.25), 'q3': quantile(s, 0.75),
            'mean': sum(s) / len(s)}


def wilcoxon(diffs):
    """two-sided Wilcoxon signed-rank p-value; zeros dropped, tied |d| get average ranks.
    Exact (permutation distribution of the observed ranks) for n <= 60, else normal approx."""
    d = [x for x in diffs if x != 0]
    n = len(d)
    if n == 0:
        return None
    order = sorted(range(n), key=lambda i: abs(d[i]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(d[order[j + 1]]) == abs(d[order[i]]):
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2.0 + 1
        i = j + 1
    wplus = sum(r for r, x in zip(ranks, d) if x > 0)
    if n <= 60:
        r2 = [int(round(2 * r)) for r in ranks]           # doubled ranks are integers
        tot = sum(r2)
        dist = [0] * (tot + 1)
        dist[0] = 1
        for r in r2:
            for s in range(tot, r - 1, -1):
                dist[s] += dist[s - r]
        allc = float(2 ** n)
        w2 = int(round(2 * wplus))
        lo = sum(dist[:w2 + 1]) / allc
        hi = sum(dist[w2:]) / allc
        return min(1.0, 2 * min(lo, hi))
    mu = n * (n + 1) / 4.0
    ties = collections.Counter(ranks)
    var = n * (n + 1) * (2 * n + 1) / 24.0 - sum(t ** 3 - t for t in ties.values()) / 48.0
    if var <= 0:
        return None
    z = (abs(wplus - mu) - 0.5) / math.sqrt(var)
    return min(1.0, math.erfc(z / math.sqrt(2)))


def load(results, methods, regions_tsv):
    strat = {}
    if regions_tsv and os.path.exists(regions_tsv):
        with open(regions_tsv) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                strat[r['region_id']] = r.get('stratum')
    data = collections.defaultdict(dict)       # method -> region -> json
    invalid = collections.Counter()
    for m in methods:
        for p in sorted(glob.glob(os.path.join(results, m, '*.json'))):
            try:
                d = json.load(open(p))
            except (ValueError, OSError):
                invalid[m] += 1
                continue
            rid = d.get('region_id') or os.path.basename(p)[:-5]
            if d.get('status') != 'ok':
                invalid[m] += 1
                continue
            if not d.get('stratum'):
                d['stratum'] = strat.get(rid, 'unknown')
            data[m][rid] = d
    return data, invalid


def pooled_f1(ds, kind):
    t = collections.Counter()
    n = 0
    for d in ds:
        tv = (d.get('truth') or {}).get('truvari') or {}
        vals = [tv.get('%s_%s' % (kind, k)) for k in ('TP-base', 'FN', 'TP-comp', 'FP')]
        if None in vals:
            continue
        n += 1
        for k, v in zip(('tpb', 'fn', 'tpc', 'fp'), vals):
            t[k] += v
    rec = t['tpb'] / (t['tpb'] + t['fn']) if t['tpb'] + t['fn'] else None
    prec = t['tpc'] / (t['tpc'] + t['fp']) if t['tpc'] + t['fp'] else None
    if rec is None or prec is None:
        return n, None, t
    return n, (2 * rec * prec / (rec + prec) if rec + prec else 0.0), t


def fmt(x, nd=3):
    """integers as integers; |x| >= 100 to 0 decimals, >= 10 to 1, else up to `nd` (trailing zeros cut)"""
    if x is None:
        return '-'
    if isinstance(x, int) or float(x).is_integer():
        return str(int(x))
    a = abs(x)
    if a >= 100:
        return '%.0f' % x
    if a >= 10:
        return '%.1f' % x
    t = ('%.' + str(nd) + 'f') % x
    return t.rstrip('0').rstrip('.') if '.' in t else t


def fmtp(p):
    if p is None:
        return '-'
    return '%.2g' % p


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', default=config.RESULTS_DIR)
    ap.add_argument('--baseline', default='mc')
    ap.add_argument('--methods', default=None, help='comma list (default: every results/<METHOD>/ with JSON)')
    ap.add_argument('--regions', default=os.path.join(config.REGIONS_DIR, 'regions.tsv'))
    ap.add_argument('--metrics', default='headline', choices=('headline', 'all'))
    ap.add_argument('--strata', default=None, help='comma list of strata to report (default all)')
    ap.add_argument('--md', default=None, help='write the markdown here instead of stdout')
    ap.add_argument('--tsv', default=None, help='long-format table of every (stratum, metric, method)')
    ap.add_argument('--per-region', default=None, help='wide per-region table (region x method)')
    a = ap.parse_args(argv)

    if a.methods:
        methods = [m for m in a.methods.split(',') if m]
    else:
        methods = sorted(x for x in os.listdir(a.results)
                         if glob.glob(os.path.join(a.results, x, '*.json'))) if os.path.isdir(a.results) else []
    if a.baseline in methods:
        methods = [a.baseline] + [m for m in methods if m != a.baseline]
    if not methods:
        raise SystemExit('no results under %s' % a.results)
    data, invalid = load(a.results, methods, a.regions)
    base = data.get(a.baseline, {})
    mets = [m for m in METRICS if a.metrics == 'all' or m[3]]
    strata = sorted(set(d['stratum'] for m in methods for d in data[m].values()),
                    key=lambda s: (STRATA_ORDER.index(s) if s in STRATA_ORDER else 99, s))
    if a.strata:
        keep = set(a.strata.split(','))
        strata = [s for s in strata if s in keep]

    out = []
    rows = []
    out.append('# Candidate graphs by stratum\n')
    out.append('Results: `%s`. Baseline: `%s`. Methods: %s.\n' % (a.results, a.baseline,
                                                                ', '.join('`%s`' % m for m in methods)))
    out.append('| method | valid regions | invalid |')
    out.append('|---|---|---|')
    for m in methods:
        out.append('| %s | %d | %d |' % (m, len(data[m]), invalid[m]))
    out.append('')
    out.append('Cells: median [q1-q3] over regions. Paired columns compare each method with `%s` on the '
               'regions both have: median difference (method - baseline), regions better/worse/equal '
               '(by the metric\'s direction; \'.\' = no preferred direction), Wilcoxon signed-rank '
               'two-sided p.\n' % a.baseline)
    for st in strata:
        by = {m: {r: d for r, d in data[m].items() if d['stratum'] == st} for m in methods}
        nreg = len(set(r for m in methods for r in by[m]))
        out.append('## %s (%d regions)\n' % (st, nreg))
        # headline: medians only, methods as rows
        hl = [x for x in mets if x[3]]
        out.append('| method | n | ' + ' | '.join(x[1] for x in hl) +
                   ' | regions with cost/opt <= 1.1 | pooled truvari F1 raw / refined |')
        out.append('|---|---|' + '---|' * len(hl) + '---|---|')
        for m in methods:
            cells = []
            for key, lab, dr, _ in hl:
                q = qstats([v for v in (get(d, key) for d in by[m].values()) if v is not None])
                cells.append(fmt(q['median']) if q else '-')
            n1, f_raw, _ = pooled_f1(by[m].values(), 'raw')
            n2, f_ref, _ = pooled_f1(by[m].values(), 'refined')
            co = [get(d, 'alignment.all_cost_over_opt') for d in by[m].values()]
            co = [x for x in co if x is not None]
            tgt = '%d/%d' % (sum(1 for x in co if x <= 1.1), len(co)) if co else '-'
            out.append('| %s | %d | %s | %s | %s / %s |' % (m, len(by[m]), ' | '.join(cells), tgt, fmt(f_raw),
                                                          fmt(f_ref)))
        out.append('')
        # detail: one row per (metric, method)
        out.append('| metric | method | n | median [q1-q3] | mean | paired n | median diff [q1-q3] | '
                   'better / worse / equal | p |')
        out.append('|---|---|---|---|---|---|---|---|---|')
        for key, lab, dr, _ in mets:
            for m in methods:
                vals = {r: get(d, key) for r, d in by[m].items()}
                vals = {r: v for r, v in vals.items() if v is not None}
                q = qstats(list(vals.values()))
                row = {'stratum': st, 'metric': key, 'label': lab, 'direction': dr, 'method': m,
                       'n': q['n'] if q else 0}
                if q:
                    row.update({k: q[k] for k in ('median', 'q1', 'q3', 'mean')})
                pair = ''
                if m != a.baseline:
                    bv = {r: get(d, key) for r, d in base.items() if d['stratum'] == st}
                    common = [r for r in vals if bv.get(r) is not None]
                    diffs = [vals[r] - bv[r] for r in common]
                    if diffs:
                        qd = qstats(diffs)
                        if dr:
                            better = sum(1 for x in diffs if x * dr > 0)
                            worse = sum(1 for x in diffs if x * dr < 0)
                        else:
                            better = worse = None
                        equal = sum(1 for x in diffs if x == 0)
                        p = wilcoxon(diffs)
                        row.update({'n_paired': len(diffs), 'diff_median': qd['median'], 'diff_q1': qd['q1'],
                                    'diff_q3': qd['q3'], 'better': better, 'worse': worse, 'equal': equal,
                                    'wilcoxon_p': p})
                        bwe = ('%d / %d / %d' % (better, worse, equal)) if dr else ('. / . / %d' % equal)
                        pair = '%d | %s [%s-%s] | %s | %s' % (len(diffs), fmt(qd['median']), fmt(qd['q1']),
                                                              fmt(qd['q3']), bwe, fmtp(p))
                if not pair:
                    pair = '- | - | - | -'
                rows.append(row)
                if q:
                    out.append('| %s | %s | %d | %s [%s-%s] | %s | %s |' % (
                        lab, m, q['n'], fmt(q['median']), fmt(q['q1']), fmt(q['q3']), fmt(q['mean']), pair))
                else:
                    out.append('| %s | %s | 0 | - | - | %s |' % (lab, m, pair))
        out.append('')
        # pooled truvari
        out.append('| method | regions scored (raw / refined / phab) | raw TP-base / FN / TP-comp / FP | raw F1 | '
                   'refined TP-base / FN / TP-comp / FP | refined F1 | phab TP-base / FN / TP-comp / FP | '
                   'phab F1 |')
        out.append('|---|---|---|---|---|---|---|---|')
        for m in methods:
            n1, f1r, t1 = pooled_f1(by[m].values(), 'raw')
            n2, f1f, t2 = pooled_f1(by[m].values(), 'refined')
            n3, f1p, t3 = pooled_f1(by[m].values(), 'phab')
            out.append('| %s | %d / %d / %d | %d / %d / %d / %d | %s | %d / %d / %d / %d | %s | %d / %d / %d / %d | %s |' % (
                m, n1, n2, n3, t1['tpb'], t1['fn'], t1['tpc'], t1['fp'], fmt(f1r),
                t2['tpb'], t2['fn'], t2['tpc'], t2['fp'], fmt(f1f), t3['tpb'], t3['fn'], t3['tpc'], t3['fp'],
                fmt(f1p)))
            for kind, nn, ff, tt in (('raw', n1, f1r, t1), ('refined', n2, f1f, t2), ('phab', n3, f1p, t3)):
                rows.append({'stratum': st, 'metric': 'truth.truvari.pooled_%s_f1' % kind,
                             'label': 'pooled truvari F1 ' + kind, 'direction': 1, 'method': m, 'n': nn,
                             'median': ff, 'mean': ff, 'tp_base': tt['tpb'], 'fn': tt['fn'],
                             'tp_comp': tt['tpc'], 'fp': tt['fp']})
        out.append('')
    text = '\n'.join(out) + '\n'
    if a.md:
        os.makedirs(os.path.dirname(os.path.abspath(a.md)), exist_ok=True)
        open(a.md, 'w').write(text)
        print('wrote %s' % a.md, file=sys.stderr)
    else:
        sys.stdout.write(text)
    if a.tsv:
        cols = ['stratum', 'metric', 'label', 'direction', 'method', 'n', 'median', 'q1', 'q3', 'mean',
                'n_paired', 'diff_median', 'diff_q1', 'diff_q3', 'better', 'worse', 'equal', 'wilcoxon_p',
                'tp_base', 'fn', 'tp_comp', 'fp']
        with open(a.tsv, 'w') as f:
            f.write('\t'.join(cols) + '\n')
            for r in rows:
                f.write('\t'.join('' if r.get(c) is None else str(r.get(c)) for c in cols) + '\n')
        print('wrote %s' % a.tsv, file=sys.stderr)
    if a.per_region:
        keys = [m[0] for m in METRICS]
        with open(a.per_region, 'w') as f:
            f.write('\t'.join(['region_id', 'stratum', 'method'] + keys) + '\n')
            for m in methods:
                for r, d in sorted(data[m].items()):
                    f.write('\t'.join([r, d['stratum'], m] +
                                      ['' if get(d, k) is None else str(get(d, k)) for k in keys]) + '\n')
        print('wrote %s' % a.per_region, file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
