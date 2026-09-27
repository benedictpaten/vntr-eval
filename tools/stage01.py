#!/usr/bin/env python3
"""stage01.py -- the Stage 0-1 comparison of every realigner against Minigraph-Cactus (mc).

    python3 tools/stage01.py [--md results/stage01_tables.md] [--tsv results/stage01_pairs.tsv]
                             [--per-region FILE]

Reads the evaluate.py results of both arms and of the full-panel graphs:
  results/<method>/<id>.json        m(hap32) and, as <method>__all, m(all) projected onto hap32
  results/full/<method>/<id>.json   full-panel graphs judged on the full panel (evaluate.py --panel)
and writes tables that answer, per stratum and per method:
  A  do realigned hotspot VNTRs look like the matched controls, and how much of the gap closes;
  B  are the controls left unchanged (regressions against mc);
  C  which method is best, on regions every method of the arm has;
  D  the panel question: m(hap32) vs m__all vs mc (and the full-panel graphs vs full MC);
  E  truth in the graph: best-path edit distance, truth-by-graph truvari F1 (raw / refined / phab);
  F  Stage 1: placement classes of the existing MAPQ < 5 reads.
Every comparison is paired over the regions where all the graphs it names exist, and states n.
A region missing for a method is listed with the reason from the runtime tables. Metric
definitions: tools/METRICS.md.
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
from summarise import get, qstats, wilcoxon, quantile  # noqa: E402

HAP32 = ['mafft_fftns2', 'mafft_linsi', 'mafft_einsi', 'mafft_ginsi', 'poa_abpoa', 'poa_spoa',
         'poa_abpoa_mc', 'unit_aware']
ALL = [m + '__all' for m in HAP32[:4]] + ['mafft_fftnsi__all'] + [m + '__all' for m in HAP32[4:]]
FULL = [m for m in HAP32[:4]] + ['mafft_fftnsi'] + HAP32[4:]
FAMILY = {'mafft': 'progressive MSA (mafft)', 'poa': 'partial order (POA)', 'unit': 'repeat-unit-aware'}
STRATA = ['hotspot_vntr', 'control_vntr_matched', 'control_vntr_correct', 'hotspot_other', 'control_nontr_sv']
SHORT = {'hotspot_vntr': 'hotspot VNTR', 'control_vntr_matched': 'matched control',
         'control_vntr_correct': 'correct control', 'hotspot_other': 'hotspot other',
         'control_nontr_sv': 'non-TR SV control'}
CONTROLS = ['control_vntr_matched', 'control_vntr_correct', 'control_nontr_sv']

# key, label, direction (-1 lower is better), tolerance for better/worse counts
M_COST = ('alignment.all_cost_over_opt', 'cost/opt', -1, 0.005)
METRICS_A = [
    M_COST,
    ('alignment.all_excess_per_kb', 'excess edits/kb', -1, 0.0),
    ('alignment.affine_all_cost_over_opt', 'affine cost/opt', -1, 0.005),
    ('size.nodes_per_kb', 'nodes/kb', -1, 0.0),
    ('size.node_len_mean', 'mean node bp', +1, 0.0),
    ('inflation.sv_pieces_mean', 'SV pieces/path', -1, 0.0),
    ('redundancy.kmer_frac_extra', 'k-mer extra', -1, 0.0),
    ('alignment.all_U_per_kb', 'parallel homologous bp/kb', -1, 0.0),
]
METRICS_EXTRA = [
    ('alignment.ref_cost_over_opt', 'cost/opt vs CHM13', -1, 0.005),
    ('alignment.affine_ref_cost_over_opt', 'affine cost/opt vs CHM13', -1, 0.005),
    ('inflation.indel_bp_ratio', 'indel bp / net', -1, 0.0),
    ('inflation.anchored_over_free', 'anchored/free', -1, 0.005),
    ('truth.d_graph_sum', 'truth edits to graph (h1+h2)', -1, 0.0),
    ('reads.mqlt5_redundant_only', 'MAPQ<5 reads redundant only', -1, 0.0),
    ('reads.mqlt5_redundant', 'MAPQ<5 reads with excess placements', -1, 0.0),
    ('reads.mqlt5_ratio_mean', 'MAPQ<5 placements / max per path', -1, 0.0),
]
TARGET = 1.1


def link(rid):
    return '[%s](../regions/%s/)' % (rid, rid)


def f(x, nd=3):
    if x is None:
        return '-'
    if isinstance(x, float) and math.isnan(x):
        return '-'
    if isinstance(x, int) or float(x).is_integer() and abs(x) >= 10:
        return '%d' % x
    a = abs(x)
    if a >= 100:
        return '%.0f' % x
    if a >= 10:
        return '%.1f' % x
    return ('%.' + str(nd) + 'f') % x


def pct(x):
    return '-' if x is None else '%.0f%%' % (100 * x)


def med(vals):
    v = [x for x in vals if x is not None]
    return quantile(sorted(v), 0.5) if v else None


# ---------------------------------------------------------------- loading
def load_dir(d):
    out, bad = {}, {}
    for p in sorted(glob.glob(os.path.join(d, '*.json'))):
        try:
            j = json.load(open(p))
        except (ValueError, OSError) as e:
            bad[os.path.basename(p)[:-5]] = 'unreadable: %s' % e
            continue
        rid = j.get('region_id') or os.path.basename(p)[:-5]
        if j.get('status') != 'ok':
            bad[rid] = j.get('status')
            continue
        out[rid] = j
    return out, bad


def load_regions():
    regs = {}
    with open(os.path.join(config.REGIONS_DIR, 'regions.tsv')) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            r['span_bp'] = int(r['span_end']) - int(r['span_start']) + 1
            regs[r['region_id']] = r
    return regs


def load_status():
    """(method, region) -> 'status: note' for every region a method did not finish"""
    st = {}
    R = config.RESULTS_DIR

    def rd(name):
        p = os.path.join(R, name)
        return list(csv.DictReader(open(p), delimiter='\t')) if os.path.exists(p) else []
    for r in rd('realign_runtime.tsv'):
        st[(r['method'], r['region_id'])] = (r['status'], r.get('note', ''))
    for name in ('realign_runtime.all.mafft.tsv', 'realign_runtime.all.poa.tsv', 'realign_runtime.all.units.tsv'):
        for r in rd(name):
            st[(r['method'] + '__all', r['region_id'])] = (r['status'], r.get('note', ''))
    for r in rd('mcpoa_runtime.tsv'):
        st[('poa_abpoa_mc', r['region_id'])] = (r['status'], r.get('note', ''))
    for r in rd('mcpoa_runtime.all.tsv'):
        m = r['method'] if r['method'].endswith('__all') else r['method'] + '__all'
        st[(m, r['region_id'])] = (r['status'], r.get('note', ''))
    return st


class Data:
    def __init__(self):
        self.regs = load_regions()
        self.status = load_status()
        self.h = {}      # method -> rid -> json  (hap32-scored graphs: mc, m, m__all)
        self.full = {}   # method -> rid -> json  (full-panel graphs)
        self.bad = {}
        R = config.RESULTS_DIR
        for m in ['mc'] + HAP32 + ALL:
            self.h[m], b = load_dir(os.path.join(R, m))
            self.bad.update({(m, k): v for k, v in b.items()})
        for m in ['mc'] + FULL:
            self.full[m], b = load_dir(os.path.join(R, 'full', m))
            self.bad.update({('full/' + m, k): v for k, v in b.items()})

    def stratum(self, rid):
        return self.regs[rid]['stratum']

    def rids(self, stratum, *sets):
        """regions of a stratum present in every one of the given dicts"""
        out = [r for r in self.regs if self.regs[r]['stratum'] == stratum]
        for s in sets:
            out = [r for r in out if r in s]
        return sorted(out)

    def why_missing(self, method, rid):
        base = method[5:] if method.startswith('full/') else method
        if base.startswith('full/'):
            base = base[5:]
        s = self.status.get((base, rid))
        if s is None and not method.endswith('__all') and method.startswith('full/'):
            s = self.status.get((base + '__all', rid))
        if s is None:
            if (method, rid) in self.bad:
                return 'evaluate: %s' % self.bad[(method, rid)]
            return 'no graph'
        return s[0] + ((': ' + s[1][:90]) if s[1] else '')


def vals(D, rids, key):
    return [get(D[r], key) for r in rids]


def pair_stats(A, B, rids, key, direction, tol):
    """B relative to A over rids: n, med A, med B, med diff, better, worse, p"""
    a = [get(A[r], key) for r in rids]
    b = [get(B[r], key) for r in rids]
    pr = [(x, y) for x, y in zip(a, b) if x is not None and y is not None]
    if not pr:
        return None
    d = [y - x for x, y in pr]
    better = sum(1 for x in d if x * direction > tol)
    worse = sum(1 for x in d if -x * direction > tol)
    return {'n': len(pr), 'a': med([x for x, _ in pr]), 'b': med([y for _, y in pr]), 'diff': med(d),
            'better': better, 'worse': worse, 'equal': len(pr) - better - worse,
            'p': wilcoxon(d) if len(d) > 1 else None}


def n_target(D, rids, key='alignment.all_cost_over_opt'):
    v = [get(D[r], key) for r in rids]
    v = [x for x in v if x is not None]
    return sum(1 for x in v if x <= TARGET), len(v)


# ---------------------------------------------------------------- tables
class Out:
    def __init__(self):
        self.lines = []
        self.tsv = []

    def __call__(self, s=''):
        self.lines.append(s)

    def table(self, header, rows):
        self('| ' + ' | '.join(header) + ' |')
        self('|' + '---|' * len(header))
        for r in rows:
            self('| ' + ' | '.join(str(x) for x in r) + ' |')
        self()

    def rec(self, **kw):
        self.tsv.append(kw)


def table_A(D, o, arm_methods, arm_label, src='h'):
    """hotspot VNTRs vs matched controls; gap closure"""
    H = D.h if src == 'h' else D.full
    mc = H['mc']
    o('### %s\n' % arm_label)
    ctrl_all = D.rids('control_vntr_matched', mc)
    hot_all = D.rids('hotspot_vntr', mc)
    hdr = ['method', 'n hot / ctrl'] + [lab for _, lab, _, _ in METRICS_A] + ['hot <= 1.1', 'ctrl <= 1.1']
    rows = []
    cells = ['mc', '%d / %d' % (len(hot_all), len(ctrl_all))]
    for key, lab, dr, tol in METRICS_A:
        cells.append('%s / %s' % (f(med(vals(mc, hot_all, key))), f(med(vals(mc, ctrl_all, key)))))
    k1, n1 = n_target(mc, hot_all)
    k2, n2 = n_target(mc, ctrl_all)
    cells += ['%d/%d' % (k1, n1), '%d/%d' % (k2, n2)]
    rows.append(cells)
    clos_rows = []
    for m in arm_methods:
        if m not in H or not H[m]:
            continue
        hot = D.rids('hotspot_vntr', mc, H[m])
        ctrl = D.rids('control_vntr_matched', mc, H[m])
        cells = ['`%s`' % m, '%d / %d' % (len(hot), len(ctrl))]
        crow = ['`%s`' % m, str(len(hot))]
        for key, lab, dr, tol in METRICS_A:
            mh, mcn = med(vals(H[m], hot, key)), med(vals(H[m], ctrl, key))
            cells.append('%s / %s' % (f(mh), f(mcn)))
            base_h = med(vals(mc, hot, key))
            base_c = med(vals(mc, ctrl_all, key))
            clo = None
            if None not in (mh, base_h, base_c) and base_h != base_c:
                clo = (base_h - mh) / (base_h - base_c)
            crow.append('%s -> %s (%s)' % (f(base_h), f(mh), pct(clo)))
            o.rec(table='A_gap', arm=arm_label, method=m, stratum='hotspot_vntr', metric=key, n=len(hot),
                  mc_hot=base_h, m_hot=mh, mc_ctrl_all=base_c, m_ctrl=mcn, gap_closed=clo)
        k1, n1 = n_target(H[m], hot)
        k2, n2 = n_target(H[m], ctrl)
        kb, _ = n_target(mc, hot)
        cells += ['%d/%d' % (k1, n1), '%d/%d' % (k2, n2)]
        crow.append('%d -> %d of %d' % (kb, k1, n1))
        rows.append(cells)
        clos_rows.append(crow)
    o('Medians, hotspot VNTRs / matched VNTR controls, each method on the regions it has (n given; '
      'mc on all of them).\n')
    o.table(hdr, rows)
    o('Gap closure at hotspot VNTRs: mc median -> method median on the same hotspots, and in brackets the '
      'share of the distance from mc\'s hotspot median to mc\'s matched-control median (all %d) that it '
      'covers (100%% = hotspots reach the matched controls as they are in mc; >100%% = beyond). Last '
      'column: hotspots at cost/opt <= 1.1, mc -> method.\n' % len(ctrl_all))
    o.table(['method', 'n hot'] + [lab for _, lab, _, _ in METRICS_A] + ['<= 1.1'], clos_rows)


def table_B(D, o, methods, arm_label):
    mc = D.h['mc']
    o('### %s\n' % arm_label)
    rows = []
    worst = []
    for m in methods:
        if not D.h.get(m):
            continue
        for st in CONTROLS + ['hotspot_vntr', 'hotspot_other']:
            rids = D.rids(st, mc, D.h[m])
            ps = pair_stats(mc, D.h[m], rids, *M_COST[:1], M_COST[2], M_COST[3])
            if not ps:
                continue
            d = {r: get(D.h[m][r], M_COST[0]) - get(mc[r], M_COST[0]) for r in rids}
            reg01 = sorted([r for r in rids if d[r] > 0.01], key=lambda r: -d[r])
            reg05 = [r for r in reg01 if d[r] > 0.05]
            aff = pair_stats(mc, D.h[m], rids, 'alignment.affine_all_cost_over_opt', -1, 0.005)
            rows.append(['`%s`' % m, SHORT[st], ps['n'], f(ps['a']), f(ps['b']), f(ps['diff'], 4),
                         '%d / %d / %d' % (ps['better'], ps['worse'], ps['equal']),
                         '%.2g' % ps['p'] if ps['p'] is not None else '-', len(reg01), len(reg05),
                         '%d / %d' % (aff['better'], aff['worse']) if aff else '-',
                         ', '.join('%s %+.3f' % (link(r), d[r]) for r in reg01[:3]) or '-'])
            o.rec(table='B_controls', arm=arm_label, method=m, stratum=st, metric=M_COST[0], **ps,
                  regress_gt_0p01=len(reg01), regress_gt_0p05=len(reg05))
            if st in CONTROLS:
                for r in reg05:
                    worst.append((m, st, r, get(mc[r], M_COST[0]), get(D.h[m][r], M_COST[0])))
    o('Paired against mc on the regions both have. b/w/e = better / worse / equal by more than 0.005 '
      'cost/opt; p = Wilcoxon signed-rank; regressions = regions worse than mc by more than 0.01 / 0.05; '
      'affine b/w by more than 0.005.\n')
    o.table(['method', 'stratum', 'n', 'mc', 'method', 'median diff', 'b / w / e', 'p', 'regress > 0.01',
             '> 0.05', 'affine b / w', 'worst (diff)'], rows)
    return worst


def table_C(D, o, methods, arm_label, extra_keys=True):
    mc = D.h['mc']
    sets = [mc] + [D.h[m] for m in methods if D.h.get(m)]
    methods = [m for m in methods if D.h.get(m)]
    o('### %s\n' % arm_label)
    common_all = [r for r in D.regs if all(r in s for s in sets)]
    o('Regions where mc and all %d methods have a graph: %d of %d.\n' % (len(methods), len(common_all), len(D.regs)))
    keys = [M_COST, METRICS_A[2], METRICS_A[3], METRICS_A[5], METRICS_A[6], METRICS_A[7],
            METRICS_EXTRA[4], METRICS_EXTRA[6]]
    for key, lab, dr, tol in keys:
        rows = []
        for st in STRATA:
            rids = D.rids(st, *sets)
            if not rids:
                continue
            cells = [SHORT[st], len(rids)]
            for m in ['mc'] + methods:
                v = vals(mc if m == 'mc' else D.h[m], rids, key)
                if key == 'truth.d_graph_sum':
                    vv = [x for x in v if x is not None]
                    cells.append('%d' % sum(vv) if vv else '-')
                elif key == M_COST[0]:
                    k, n = n_target(mc if m == 'mc' else D.h[m], rids)
                    cells.append('%s (%d)' % (f(med(v)), k))
                else:
                    cells.append(f(med(v)))
            rows.append(cells)
        what = 'sum over regions' if key == 'truth.d_graph_sum' else 'median'
        extra = ' (k) = regions <= 1.1.' if key == M_COST[0] else ''
        o('**%s** (%s).%s\n' % (lab, what, extra))
        o.table(['stratum', 'n', 'mc'] + ['`%s`' % m for m in methods], rows)
    # best method per region
    rows = []
    for st in STRATA:
        rids = D.rids(st, *sets)
        if not rids:
            continue
        wins = collections.Counter()
        for r in rids:
            v = {m: get(D.h[m][r], M_COST[0]) for m in methods}
            v = {m: x for m, x in v.items() if x is not None}
            if not v:
                continue
            best = min(v.values())
            for m, x in v.items():
                if x <= best + 0.005:
                    wins[m] += 1
        rows.append([SHORT[st], len(rids)] + [wins[m] for m in methods])
    o('**Regions where each method is best or within 0.005 of the best cost/opt** (ties count for every '
      'method within 0.005).\n')
    o.table(['stratum', 'n'] + ['`%s`' % m for m in methods], rows)
    return common_all


def table_D(D, o):
    mc = D.h['mc']
    o('### m(hap32) vs m__all (full panel aligned, projected onto hap32) vs mc, all scored on hap32\n')
    o('Paired over the regions where mc, m(hap32) and m__all all exist. Coverage: regions with an m__all '
      'graph out of 149, and hotspot VNTRs among them. Retained = (mc - m__all) / (mc - m(hap32)) on the '
      'medians, shown where m(hap32) gained at least 0.02 on mc. all-vs-hap32 = regions where m__all is better / worse than m(hap32) by more than 0.005; '
      'all-vs-mc likewise against mc; p = Wilcoxon, m__all vs m(hap32).\n')
    methods = [m for m in HAP32] + ['mafft_fftnsi']
    rows = []
    for m in methods:
        ma = m + '__all'
        if not D.h.get(ma):
            continue
        cov = len(D.h[ma])
        covh = len([r for r in D.h[ma] if D.stratum(r) == 'hotspot_vntr'])
        for st in STRATA:
            if m == 'mafft_fftnsi':
                continue
            rids = D.rids(st, mc, D.h.get(m, {}), D.h[ma])
            if not rids:
                continue
            key, lab, dr, tol = M_COST
            a = med(vals(mc, rids, key))
            b = med(vals(D.h[m], rids, key))
            c = med(vals(D.h[ma], rids, key))
            # only meaningful where m(hap32) gained something to retain (>= 0.02 cost/opt)
            ret = (a - c) / (a - b) if None not in (a, b, c) and a - b >= 0.02 else None
            p1 = pair_stats(D.h[m], D.h[ma], rids, key, dr, tol)
            p2 = pair_stats(mc, D.h[ma], rids, key, dr, tol)
            aff = [med(vals(X, rids, 'alignment.affine_all_cost_over_opt')) for X in (mc, D.h[m], D.h[ma])]
            npk = [med(vals(X, rids, 'size.nodes_per_kb')) for X in (mc, D.h[m], D.h[ma])]
            k = [n_target(X, rids)[0] for X in (mc, D.h[m], D.h[ma])]
            rows.append(['`%s`' % m, '%d (%d hot)' % (cov, covh), SHORT[st], len(rids),
                         '%s / %s / %s' % (f(a), f(b), f(c)), pct(ret),
                         '%d / %d' % (p1['better'], p1['worse']), '%.2g' % p1['p'] if p1['p'] is not None else '-',
                         '%d / %d' % (p2['better'], p2['worse']),
                         '%d / %d / %d' % tuple(k),
                         '%s / %s / %s' % tuple(f(x) for x in aff), '%s / %s / %s' % tuple(f(x) for x in npk)])
            o.rec(table='D_panel', method=m, stratum=st, n=len(rids), mc=a, hap32=b, all=c, retained=ret,
                  all_vs_hap32_better=p1['better'], all_vs_hap32_worse=p1['worse'], p=p1['p'],
                  all_vs_mc_better=p2['better'], all_vs_mc_worse=p2['worse'])
    o.table(['method', 'm__all coverage', 'stratum', 'n', 'cost/opt mc / m / m__all', 'gain retained',
             'all-vs-hap32 b / w', 'p', 'all-vs-mc b / w', '<= 1.1 mc / m / m__all',
             'affine mc / m / m__all', 'nodes/kb mc / m / m__all'], rows)
    # FFT-NS-i: no hap32 arm; against mc only
    rows = []
    ma = 'mafft_fftnsi__all'
    for st in STRATA:
        rids = D.rids(st, mc, D.h.get(ma, {}))
        if not rids:
            continue
        p2 = pair_stats(mc, D.h[ma], rids, *M_COST[:1], M_COST[2], M_COST[3])
        rows.append([SHORT[st], len(rids), f(p2['a']), f(p2['b']), '%d / %d' % (p2['better'], p2['worse'])])
    if rows:
        o('`mafft_fftnsi__all` has no hap32 arm; against mc only:\n')
        o.table(['stratum', 'n', 'mc', 'mafft_fftnsi__all', 'b / w vs mc'], rows)


def table_D_full(D, o):
    mc_h, mc_f = D.h['mc'], D.full['mc']
    o('### Full MC on the full panel vs mc.gfa on hap32 (same regions, paired)\n')
    o('Full-panel scores use 300 sampled haplotype pairs (seeded, the same pairs for every graph of a '
      'region) and CHM13 x 300 haplotypes; hap32 scores use all pairs.\n')
    keys = [M_COST, ('alignment.ref_cost_over_opt', 'cost/opt vs CHM13', -1, 0.005)] + METRICS_A[1:]
    rows = []
    for st in STRATA:
        rids = D.rids(st, mc_h, mc_f)
        cells = [SHORT[st], len(rids)]
        for key, lab, dr, tol in keys:
            cells.append('%s / %s' % (f(med(vals(mc_h, rids, key))), f(med(vals(mc_f, rids, key)))))
        k1, _ = n_target(mc_h, rids)
        k2, _ = n_target(mc_f, rids)
        p = pair_stats(mc_h, mc_f, rids, *M_COST[:1], M_COST[2], M_COST[3])
        cells += ['%d / %d' % (k1, k2), '%d / %d' % (p['better'], p['worse'])]
        rows.append(cells)
    o.table(['stratum', 'n'] + ['%s hap32 / full' % lab for _, lab, _, _ in keys] +
            ['<= 1.1 hap32 / full', 'full better / worse'], rows)
    o('### Realigned full-panel graphs vs the full MC graph (both judged on the full panel)\n')
    o('Paired over regions where both full-panel graphs exist; b / w by more than 0.005 cost/opt.\n')
    rows = []
    for m in FULL:
        F = D.full.get(m)
        if not F:
            continue
        for st in STRATA:
            rids = D.rids(st, mc_f, F)
            if not rids:
                continue
            p = pair_stats(mc_f, F, rids, *M_COST[:1], M_COST[2], M_COST[3])
            aff = pair_stats(mc_f, F, rids, 'alignment.affine_all_cost_over_opt', -1, 0.005)
            ref = pair_stats(mc_f, F, rids, 'alignment.ref_cost_over_opt', -1, 0.005)
            npk = [med(vals(X, rids, 'size.nodes_per_kb')) for X in (mc_f, F)]
            u = [med(vals(X, rids, 'alignment.all_U_per_kb')) for X in (mc_f, F)]
            km = [med(vals(X, rids, 'redundancy.kmer_frac_extra')) for X in (mc_f, F)]
            sv = [med(vals(X, rids, 'inflation.sv_pieces_mean')) for X in (mc_f, F)]
            k = [n_target(X, rids)[0] for X in (mc_f, F)]
            rows.append(['`%s`' % m, SHORT[st], p['n'], '%s -> %s' % (f(p['a']), f(p['b'])),
                         '%d / %d' % (p['better'], p['worse']), '%.2g' % p['p'] if p['p'] is not None else '-',
                         '%d -> %d' % tuple(k), '%s -> %s' % (f(ref['a']), f(ref['b'])),
                         '%s -> %s' % (f(aff['a']), f(aff['b'])),
                         '%s -> %s' % tuple(f(x) for x in npk), '%s -> %s' % tuple(f(x) for x in u),
                         '%s -> %s' % tuple(f(x) for x in km), '%s -> %s' % tuple(f(x) for x in sv)])
            o.rec(table='D_full', method=m, stratum=st, **p)
    o.table(['method', 'stratum', 'n', 'cost/opt full MC -> m', 'b / w', 'p', '<= 1.1', 'vs CHM13',
             'affine', 'nodes/kb', 'parallel hom. bp/kb', 'k-mer extra', 'SV pieces/path'], rows)
    # gap closure on the full panel
    table_A(D, o, FULL, 'Hotspot VNTRs vs matched controls on the full panel (results/full)', src='full')


def truth_block(D, o, methods, arm_label):
    mc = D.h['mc']
    o('### %s\n' % arm_label)
    rows = []
    for m in methods:
        if not D.h.get(m):
            continue
        for st in STRATA:
            rids = D.rids(st, mc, D.h[m])
            if not rids:
                continue
            spell0 = kept0 = 0
            dsum_mc = dsum_m = 0
            worse = better = 0
            new_err = []
            for r in rids:
                tm, tn = mc[r].get('truth') or {}, D.h[m][r].get('truth') or {}
                for h in ('h1', 'h2'):
                    a, b = tm.get(h + '_d_graph'), tn.get(h + '_d_graph')
                    if a is None or b is None:
                        continue
                    dsum_mc += a
                    dsum_m += b
                    if a == 0:
                        spell0 += 1
                        if b == 0:
                            kept0 += 1
                        else:
                            new_err.append((r, h, b))
                    worse += b > a
                    better += b < a
            # pooled truvari over regions where both have counts
            pooled = {}
            for kind in ('raw', 'refined', 'phab'):
                t = {'mc': collections.Counter(), 'm': collections.Counter()}
                nreg = 0
                for r in rids:
                    ok = True
                    got = {}
                    for lab, X in (('mc', mc), ('m', D.h[m])):
                        tv = ((X[r].get('truth') or {}).get('truvari') or {})
                        v = [tv.get('%s_%s' % (kind, k)) for k in ('TP-base', 'FN', 'TP-comp', 'FP')]
                        if None in v:
                            ok = False
                        got[lab] = v
                    if not ok:
                        continue
                    nreg += 1
                    for lab in ('mc', 'm'):
                        for k, v in zip(('tpb', 'fn', 'tpc', 'fp'), got[lab]):
                            t[lab][k] += v
                pooled[kind] = (nreg, pf1(t['mc']), pf1(t['m']))
            rows.append(['`%s`' % m, SHORT[st], len(rids), '%d / %d' % (kept0, spell0),
                         '%d -> %d' % (dsum_mc, dsum_m), '%d / %d' % (better, worse),
                         '%s -> %s' % (f(pooled['raw'][1]), f(pooled['raw'][2])),
                         '%s -> %s' % (f(pooled['refined'][1]), f(pooled['refined'][2])),
                         '%s -> %s' % (f(pooled['phab'][1]), f(pooled['phab'][2])),
                         '%d/%d/%d' % (pooled['raw'][0], pooled['refined'][0], pooled['phab'][0]),
                         ', '.join('%s %s +%d' % (link(r), h, b) for r, h, b in sorted(new_err, key=lambda x: -x[2])[:3]) or '-'])
            o.rec(table='E_truth', arm=arm_label, method=m, stratum=st, n=len(rids), spellable_mc=spell0,
                  still_spellable=kept0, dsum_mc=dsum_mc, dsum_m=dsum_m, hap_better=better, hap_worse=worse,
                  raw_f1_mc=pooled['raw'][1], raw_f1_m=pooled['raw'][2], refined_f1_mc=pooled['refined'][1],
                  refined_f1_m=pooled['refined'][2], phab_f1_mc=pooled['phab'][1], phab_f1_m=pooled['phab'][2])
    o('Per method and stratum, over regions both graphs have. Spellable kept = HG002 haplotypes that mc '
      'spells exactly (d_graph = 0) and the method still spells exactly / all haplotypes mc spells. '
      'Edits = summed best-path edit distance over both haplotypes, mc -> method; haplotypes b / w = '
      'haplotypes whose distance fell / rose. Truvari F1 pooled (TP/FN/FP summed) over regions where both '
      'graphs were scored (count raw/refined/phab), mc -> method. New errors: the largest distances where '
      'mc spelled the haplotype exactly.\n')
    o.table(['method', 'stratum', 'n', 'spellable kept', 'edits mc -> m', 'haplotypes b / w',
             'raw F1', 'refined F1', 'phab F1', 'truvari regions', 'new errors (largest)'], rows)


def pf1(t):
    rec = t['tpb'] / (t['tpb'] + t['fn']) if t['tpb'] + t['fn'] else None
    prec = t['tpc'] / (t['tpc'] + t['fp']) if t['tpc'] + t['fp'] else None
    if rec is None or prec is None:
        return None
    return 2 * rec * prec / (rec + prec) if rec + prec else 0.0


def reads_block(D, o, methods, arm_label):
    mc = D.h['mc']
    o('### %s\n' % arm_label)
    rows = []
    classes = ['unique', 'tandem', 'redundant_only', 'redundant', 'absent']
    for m in methods:
        if not D.h.get(m):
            continue
        for st in STRATA:
            rids = [r for r in D.rids(st, mc, D.h[m])
                    if (mc[r].get('reads') or {}).get('mqlt5_n') and (D.h[m][r].get('reads') or {}).get('mqlt5_n')]
            if not rids:
                continue
            tot = {'mc': collections.Counter(), 'm': collections.Counter()}
            n = 0
            for r in rids:
                for lab, X in (('mc', mc), ('m', D.h[m])):
                    rd = X[r]['reads']
                    k = rd['mqlt5_n']
                    for c in classes:
                        tot[lab][c] += rd.get('mqlt5_' + c, 0) * k
                    tot[lab]['ratio'] += rd.get('mqlt5_ratio_mean', 0) * k
                n += mc[r]['reads']['mqlt5_n']
            ps = pair_stats(mc, D.h[m], rids, 'reads.mqlt5_redundant_only', -1, 0.0)
            ps2 = pair_stats(mc, D.h[m], rids, 'reads.mqlt5_redundant', -1, 0.0)
            rows.append(['`%s`' % m, SHORT[st], len(rids), n] +
                        ['%s -> %s' % (pct(tot['mc'][c] / n), pct(tot['m'][c] / n)) for c in classes] +
                        ['%s -> %s' % (f(tot['mc']['ratio'] / n, 2), f(tot['m']['ratio'] / n, 2)),
                         '%d / %d' % (ps['better'], ps['worse']), '%d / %d' % (ps2['better'], ps2['worse'])])
            o.rec(table='F_reads', arm=arm_label, method=m, stratum=st, n=len(rids), reads=n,
                  **{'mc_' + c: tot['mc'][c] / n for c in classes}, **{'m_' + c: tot['m'][c] / n for c in classes})
    o('MAPQ < 5 core reads of each region (MAPQ from the hap32 mapping), pooled over reads, mc -> method, '
      'over regions both graphs have. Classes: unique (1 placement), tandem (some haplotype holds the '
      'sequence twice), redundant only (several placements, no haplotype holds it twice: removable by '
      'alignment), redundant (placements > max(1, copies in one haplotype), which includes tandem reads with '
      'excess placements), absent (no panel path spells it). Ratio = placements / most copies in one '
      'path. b / w = regions whose fraction fell / rose.\n')
    o.table(['method', 'stratum', 'n', 'reads'] + classes + ['placement ratio', 'redundant only b / w',
                                                             'redundant b / w'], rows)


def missing_block(D, o, methods):
    o('### Regions without a graph, per method\n')
    rows = []
    for m in methods:
        miss = [r for r in D.regs if r not in D.h.get(m, {})]
        if not miss:
            continue
        by = collections.Counter(D.stratum(r) for r in miss)
        why = collections.Counter(D.why_missing(m, r).split(':')[0] for r in miss)
        rows.append(['`%s`' % m, len(D.h.get(m, {})), len(miss),
                     ', '.join('%s %d' % (SHORT[s], by[s]) for s in STRATA if by[s]),
                     ', '.join('%s %d' % kv for kv in why.most_common()),
                     ', '.join(sorted(miss, key=lambda r: -D.regs[r]['span_bp'])[:6]) + (' ...' if len(miss) > 6 else '')])
    o.table(['method', 'graphs', 'missing', 'missing by stratum', 'reason', 'largest missing'], rows)


def per_region(D, path):
    keys = [k for k, _, _, _ in METRICS_A + METRICS_EXTRA] + ['truth.truvari.raw_f1', 'truth.truvari.refined_f1',
                                                               'truth.truvari.phab_f1']
    with open(path, 'w') as fh:
        fh.write('\t'.join(['region_id', 'stratum', 'span_bp', 'panel', 'method'] + keys) + '\n')
        for src, H in (('hap32', D.h), ('full', D.full)):
            for m, X in H.items():
                for r in sorted(X):
                    fh.write('\t'.join([r, D.stratum(r), str(D.regs[r]['span_bp']), src, m] +
                                       ['' if get(X[r], k) is None else str(get(X[r], k)) for k in keys]) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--md', default=os.path.join(config.RESULTS_DIR, 'stage01_tables.md'))
    ap.add_argument('--tsv', default=os.path.join(config.RESULTS_DIR, 'stage01_pairs.tsv'))
    ap.add_argument('--per-region', default=None, help='also write a long per-region table (both arms, full panel)')
    a = ap.parse_args()
    D = Data()
    o = Out()
    o('# Stage 0-1 tables\n')
    o('Generated by `python3 tools/stage01.py` from `results/<method>/` (hap32-scored graphs) and '
      '`results/full/<method>/` (full-panel graphs). The narrative is in [stage01.md](stage01.md); metric '
      'definitions in [../tools/METRICS.md](../tools/METRICS.md).\n')
    o('Evaluated graphs: ' + ', '.join('`%s` %d' % (m, len(D.h[m])) for m in ['mc'] + HAP32 + ALL if D.h.get(m)) + '.\n')
    o('Full-panel graphs: ' + ', '.join('`%s` %d' % (m, len(D.full[m])) for m in ['mc'] + FULL if D.full.get(m)) + '.\n')
    if D.bad:
        o('Not ok: ' + ', '.join('%s %s (%s)' % (m, r, s) for (m, r), s in sorted(D.bad.items())) + '.\n')
    o('## Coverage\n')
    missing_block(D, o, HAP32 + ALL)
    o('## A. Hotspot VNTRs against matched controls\n')
    table_A(D, o, HAP32, 'hap32 arm: m(hap32)')
    table_A(D, o, ALL, 'full-panel arm projected onto hap32: m__all')
    o('## B. Controls (and hotspots) against mc, cost/opt\n')
    worst_h = table_B(D, o, HAP32, 'hap32 arm')
    worst_a = table_B(D, o, ALL, 'full-panel arm, projected')
    o('## C. Methods compared on common regions\n')
    table_C(D, o, HAP32, 'hap32 arm, all eight methods')
    table_C(D, o, ['mafft_fftns2__all', 'poa_abpoa__all', 'poa_spoa__all', 'poa_abpoa_mc__all', 'unit_aware__all'],
            'full-panel arm projected, the five methods that ran on most regions')
    table_C(D, o, ['mafft_linsi__all', 'mafft_einsi__all', 'mafft_ginsi__all', 'unit_aware__all', 'poa_abpoa__all'],
            'full-panel arm projected, with the all-pairs mafft modes (small and medium regions only)')
    o('## D. The panel question\n')
    table_D(D, o)
    table_D_full(D, o)
    o('## E. Truth in the graph\n')
    truth_block(D, o, HAP32, 'hap32 arm')
    truth_block(D, o, ALL, 'full-panel arm, projected')
    o('## F. Stage 1: existing reads\n')
    reads_block(D, o, HAP32, 'hap32 arm')
    reads_block(D, o, ALL, 'full-panel arm, projected')
    open(a.md, 'w').write('\n'.join(o.lines) + '\n')
    cols = sorted(set(k for r in o.tsv for k in r), key=lambda k: (k not in ('table', 'arm', 'method', 'stratum', 'metric', 'n'), k))
    with open(a.tsv, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in o.tsv:
            fh.write('\t'.join('' if r.get(c) is None else str(r.get(c)) for c in cols) + '\n')
    if a.per_region:
        per_region(D, a.per_region)
    print('wrote %s, %s%s' % (a.md, a.tsv, (', ' + a.per_region) if a.per_region else ''), file=sys.stderr)


if __name__ == '__main__':
    main()
