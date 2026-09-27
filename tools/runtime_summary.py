#!/usr/bin/env python3
"""Summarise aligner runtime and peak memory per method and panel from the runtime tables.

    python3 tools/runtime_summary.py            # writes results/runtime_summary.{tsv,md}

Sources (one row per region x method; written by the realign drivers):
    results/realign_runtime.tsv              hap32 arm: mafft, abPOA, spoa, unit_aware (tools/realign.py)
    results/mcpoa_runtime.tsv                hap32 arm: abPOA with Cactus BAR settings (tools/realign_poa_mc.py)
    results/realign_runtime.all.mafft.tsv    full panel: mafft (tools/mafft_panel.py)
    results/realign_runtime.all.poa.tsv      full panel: abPOA, spoa (tools/poa_panel.py)
    results/realign_runtime.all.units.tsv    full panel: unit_aware (tools/units_panel.py)
    results/mcpoa_runtime.all.tsv            full panel: abPOA with Cactus BAR settings (tools/poa_panel_mc.py)

Time is the aligner's wall-clock seconds (align_s), not CPU time, measured on a shared, often heavily
loaded 10-core machine. Peak memory is exact (/usr/bin/time -l) for the Cactus-settings abPOA drivers
and sampled from `ps` about once a second over the process tree for every other driver, so those are
lower bounds (short spikes are missed, and macOS compresses memory under pressure).

Outcome classes: ok; refused = not run because a predicted time or memory was over the cap;
timeout / memout = killed while running. Refused regions are the largest ones, so every
per-method statistic over 'ok' rows is censored from above; the 'common' columns compare methods on
the regions every compared method finished.
"""
import csv
import os
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), 'results')

SOURCES = [
    ('hap32', 'realign_runtime.tsv', 'sampled'),
    ('hap32', 'mcpoa_runtime.tsv', 'exact'),
    ('full', 'realign_runtime.all.mafft.tsv', 'sampled'),
    ('full', 'realign_runtime.all.poa.tsv', 'sampled'),
    ('full', 'realign_runtime.all.units.tsv', 'sampled'),
    ('full', 'mcpoa_runtime.all.tsv', 'exact'),
]
ORDER = ['mafft_fftns2', 'mafft_fftnsi', 'mafft_linsi', 'mafft_einsi', 'mafft_ginsi',
         'poa_abpoa', 'poa_spoa', 'poa_abpoa_mc', 'unit_aware']
SINGLE_THREADED = {'poa_abpoa', 'poa_spoa', 'poa_abpoa_mc'}
MIN_SAMPLED_S = 2.0


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def outcome(r):
    st = r['status']
    note = (r.get('note') or '').lower()
    if st == 'ok':
        return 'ok'
    if st == 'skipped' or note.startswith('not run') or note.startswith('not attempted'):
        return 'refused'
    return st  # timeout, memout (killed while running), error


def pct(xs, q):
    xs = sorted(xs)
    if not xs:
        return None
    k = (len(xs) - 1) * q
    lo = int(k)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def load():
    rows = []
    for panel, fn, rss_kind in SOURCES:
        for r in csv.DictReader(open(os.path.join(RES, fn)), delimiter='\t'):
            m = r['method'].replace('__all', '')
            secs = fnum(r.get('align_s'))
            if secs is None:
                secs = fnum(r.get('seconds'))
            rows.append(dict(panel=panel, method=m, region=r['region_id'], stratum=r.get('stratum', ''),
                             outcome=outcome(r), secs=secs, rss=fnum(r.get('peak_rss_mb')),
                             pred_mb=fnum(r.get('predicted_mb')), pred_s=fnum(r.get('predicted_s')),
                             threads=('1' if m in SINGLE_THREADED else (r.get('threads') or '')),
                             rss_kind=rss_kind, note=r.get('note') or ''))
    return rows


def summarise(rows):
    out = []
    for panel in ('hap32', 'full'):
        prs = [r for r in rows if r['panel'] == panel]
        methods = [m for m in ORDER if any(r['method'] == m for r in prs)]
        ok = {m: {r['region']: r for r in prs if r['method'] == m and r['outcome'] == 'ok'} for m in methods}
        # common set: regions every method with >= 100 finished regions finished (the all-pairs mafft
        # modes cover too few full-panel regions to join a common set)
        broad = [m for m in methods if len(ok[m]) >= 100]
        common = set.intersection(*[set(ok[m]) for m in broad]) if broad else set()
        for m in methods:
            mr = [r for r in prs if r['method'] == m]
            good = list(ok[m].values())
            secs = [r['secs'] for r in good if r['secs'] is not None]
            # a ~1 s ps poll cannot see a sub-second run: sampled peaks count only runs of >= MIN_SAMPLED_S
            rss = [r['rss'] for r in good if r['rss'] and (r['rss_kind'] == 'exact' or (r['secs'] or 0) >= MIN_SAMPLED_S)]
            hs = [r['secs'] for r in good if r['stratum'] == 'hotspot_vntr' and r['secs'] is not None]
            csecs = [ok[m][g]['secs'] for g in common if g in ok[m] and ok[m][g]['secs'] is not None]
            crss = [ok[m][g]['rss'] for g in common if g in ok[m] and ok[m][g]['rss']
                    and (ok[m][g]['rss_kind'] == 'exact' or (ok[m][g]['secs'] or 0) >= MIN_SAMPLED_S)]
            refused = [r for r in mr if r['outcome'] == 'refused']
            thr = sorted({r['threads'] for r in mr if r['threads']})
            out.append(dict(
                panel=panel, method=m, regions=len(mr), ok=len(good), refused=len(refused),
                timeout=sum(r['outcome'] == 'timeout' for r in mr),
                memout=sum(r['outcome'] == 'memout' for r in mr),
                threads='/'.join(thr), rss_measure=mr[0]['rss_kind'],
                sec_median=statistics.median(secs) if secs else None, sec_p90=pct(secs, 0.9),
                sec_max=max(secs) if secs else None, sec_total_h=sum(secs) / 3600 if secs else None,
                sec_median_hotspot=statistics.median(hs) if hs else None,
                rss_n=len(rss), rss_median_mb=statistics.median(rss) if rss else None, rss_p90_mb=pct(rss, 0.9),
                rss_max_mb=max(rss) if rss else None,
                common_n=len(common) if m in broad else 0,
                common_sec_median=statistics.median(csecs) if (csecs and m in broad) else None,
                common_sec_total_h=sum(csecs) / 3600 if (csecs and m in broad) else None,
                common_rss_median_mb=statistics.median(crss) if (crss and m in broad) else None,
                common_rss_max_mb=max(crss) if (crss and m in broad) else None,
                refused_msa_too_large=sum('msa' in r['note'].lower() and 'cells' in r['note'].lower() for r in refused),
                refused_pred_mb_max=max([r['pred_mb'] for r in refused if r['pred_mb'] and 'cells' not in r['note']], default=None),
                refused_pred_s_max=max([r['pred_s'] for r in refused if r['pred_s']], default=None)))
    return out


def fmt(x, d=0):
    if x is None:
        return '-'
    return f'{x:,.{d}f}'


def gb(x):
    return '-' if x is None else f'{x / 1024:.1f}'


def main():
    rows = load()
    summ = summarise(rows)
    keys = list(summ[0].keys())
    with open(os.path.join(RES, 'runtime_summary.tsv'), 'w') as fh:
        fh.write('\t'.join(keys) + '\n')
        for s in summ:
            fh.write('\t'.join('' if s[k] is None else (f'{s[k]:.3f}' if isinstance(s[k], float) else str(s[k]))
                               for k in keys) + '\n')
    lines = ['# Aligner runtime and memory', '',
             'Generated by `tools/runtime_summary.py` from the runtime tables in `results/`; see its docstring',
             'for sources and definitions. Time is aligner wall-clock (not CPU) on a shared, often heavily',
             'loaded 10-core Mac. Peak RSS is exact for Cactus-settings abPOA (`/usr/bin/time -l`) and',
             'sampled about once a second for every other method (a lower bound). "refused" = not run because',
             'a predicted time or memory exceeded the cap; timeout/memout = killed while running.',
             f'Sampled peaks count only runs of >= {MIN_SAMPLED_S:.0f} s (a ~1 s poll misses shorter ones); n is given.',
             'unit_aware on hap32 ran without memory sampling, so it has no peak. Threads differ between methods',
             '(mafft 1-3, unit_aware 2-3, POA 1), so wall time is not CPU time. Caps: hap32 900 s per mafft call and',
             '12 GB; full panel 1800 s per region and 10-12 GB. The largest regions were refused, so the per-method',
             'statistics are censored from above; the "common" tables compare methods on shared regions.', '']
    for panel, title in (('hap32', 'Aligned on hap32 (34 sequences)'),
                         ('full', 'Aligned on the full panel (distinct sequences; projection and graph building not included)')):
        ps = [s for s in summ if s['panel'] == panel]
        lines += [f'## {title}', '',
                  '| method | threads | ok / refused / timeout / memout | time median (hotspot VNTR) | p90 | max | total h | peak RSS median / p90 / max GB | RSS | largest refused, predicted |',
                  '|---|---|---|---|---|---|---|---|---|---|']
        for s in ps:
            ref = []
            if s['refused_pred_mb_max']:
                ref.append(f"{s['refused_pred_mb_max'] / 1024:.0f} GB")
            if s['refused_pred_s_max']:
                ref.append(f"{s['refused_pred_s_max'] / 3600:.1f} h")
            if s['refused_msa_too_large']:
                ref.append(f"{s['refused_msa_too_large']} MSA too large to build")
            lines.append(
                f"| {s['method']} | {s['threads'] or '-'} | {s['ok']} / {s['refused']} / {s['timeout']} / {s['memout']} "
                f"| {fmt(s['sec_median'], 1)} s ({fmt(s['sec_median_hotspot'], 1)} s) | {fmt(s['sec_p90'])} s "
                f"| {fmt(s['sec_max'])} s | {fmt(s['sec_total_h'], 2)} "
                f"| {gb(s['rss_median_mb'])} / {gb(s['rss_p90_mb'])} / {gb(s['rss_max_mb'])} (n={s['rss_n']}) | {s['rss_measure']} "
                f"| {', '.join(ref) or '-'} |")
        cm = [s for s in ps if s['common_n']]
        if cm:
            lines += ['', f"On the {cm[0]['common_n']} regions every method with >= 100 finished regions finished:", '',
                      '| method | time median | total h | peak RSS median / max GB |', '|---|---|---|---|']
            for s in cm:
                lines.append(f"| {s['method']} | {fmt(s['common_sec_median'], 1)} s | {fmt(s['common_sec_total_h'], 2)} "
                             f"| {gb(s['common_rss_median_mb'])} / {gb(s['common_rss_max_mb'])} |")
        lines.append('')
    open(os.path.join(RES, 'runtime_summary.md'), 'w').write('\n'.join(lines))
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
