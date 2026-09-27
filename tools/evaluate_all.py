#!/usr/bin/env python3
"""evaluate_all.py -- run tools/evaluate.py over every region of one or more methods.

    python3 tools/evaluate_all.py mc mafft_linsi [--reads] [--jobs 2] [--threads 2] [--force]
                                  [--regions L012184,L009656] [--stratum hotspot_vntr]
                                  [-- <extra evaluate.py options>]
    python3 tools/evaluate_all.py mc <method> --panel [...]   # full-panel graphs (tools/panel.py)

For the method 'mc' the graph is regions/<id>/mc.gfa; for any other method it is
candidates/<method>/<id>.gfa (regions without one are listed as missing). Results go to
results/<method>/<id>.json; an existing result is kept unless --force. Regions come from
regions/regions.tsv (or every regions/*/region.json), smallest first. Keep jobs x threads within
the machine's share (default 2 x 2).

--panel: the full-panel arm. The graph is work/panel/<method>/<id>.gfa (for 'mc' the full
Minigraph-Cactus graph written by `tools/panel.py mc-graph`), judged with
`evaluate.py --panel regions/<id>/hprc.fa.gz`; results go to results/<panel-subdir>/<method>/<id>.json
(--panel-subdir, default 'panel'; the Stage 0-1 evaluation used 'full').
"""
import argparse
import concurrent.futures
import csv
import glob
import json
import os
import subprocess
import sys
import time

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import config  # noqa: E402


def regions():
    tsv = os.path.join(config.REGIONS_DIR, 'regions.tsv')
    out = []
    if os.path.exists(tsv):
        with open(tsv) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                out.append(r)
    else:
        for p in sorted(glob.glob(os.path.join(config.REGIONS_DIR, '*', 'region.json'))):
            out.append(json.load(open(p)))
    for r in out:
        try:
            r['_bp'] = int(r['span_end']) - int(r['span_start']) + 1
        except (KeyError, ValueError, TypeError):
            r['_bp'] = 0
    return sorted(out, key=lambda r: r['_bp'])


def run(job):
    method, rid, gfa, out, a, extra = job
    rd = os.path.join(config.REGIONS_DIR, rid)
    cmd = [sys.executable, os.path.join(TOOLS, 'evaluate.py'), rd, gfa, '--name', method, '--out', out,
           '--threads', str(a.threads)] + (['--reads'] if a.reads else []) + extra
    t = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    logd = os.path.join(config.WORK_DIR, 'evaluate', *([a.panel_subdir] if '--panel' in extra else []), method, rid)
    os.makedirs(logd, exist_ok=True)
    with open(os.path.join(logd, 'evaluate.log'), 'w') as f:
        f.write(p.stderr)
    return method, rid, p.returncode, time.time() - t, (p.stderr.strip().splitlines() or [''])[-1][-200:]


def main():
    argv = sys.argv[1:]
    extra = []
    if '--' in argv:
        k = argv.index('--')
        argv, extra = argv[:k], argv[k + 1:]
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('methods', nargs='+')
    ap.add_argument('--reads', action='store_true')
    ap.add_argument('--jobs', type=int, default=2)
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--regions', default=None, help='comma list of region ids')
    ap.add_argument('--stratum', default=None, help='comma list of strata')
    ap.add_argument('--panel', action='store_true',
                    help='full-panel graphs: work/panel/<method>/<id>.gfa vs regions/<id>/hprc.fa.gz')
    ap.add_argument('--panel-subdir', default='panel',
                    help='with --panel: results go to results/<this>/<method>/<id>.json (default panel)')
    a = ap.parse_args(argv)
    regs = regions()
    if a.regions:
        keep = set(a.regions.split(','))
        regs = [r for r in regs if r['region_id'] in keep]
    if a.stratum:
        keep = set(a.stratum.split(','))
        regs = [r for r in regs if r.get('stratum') in keep]
    jobs, missing = [], []
    for m in a.methods:
        for r in regs:
            rid = r['region_id']
            if a.panel:
                gfa = os.path.join(config.WORK_DIR, 'panel', m, rid + '.gfa')
                fa = os.path.join(config.REGIONS_DIR, rid, 'hprc.fa.gz')
                if not os.path.exists(fa):
                    missing.append((m, rid))
                    continue
                ex = ['--panel', fa]
                out = os.path.join(config.RESULTS_DIR, a.panel_subdir, m, rid + '.json')
            else:
                gfa = (os.path.join(config.REGIONS_DIR, rid, 'mc.gfa') if m == 'mc'
                       else os.path.join(config.CANDIDATES_DIR, m, rid + '.gfa'))
                ex = []
                out = os.path.join(config.RESULTS_DIR, m, rid + '.json')
            if not os.path.exists(gfa):
                missing.append((m, rid))
                continue
            if os.path.exists(out) and not a.force:
                continue
            jobs.append((m, rid, gfa, out, a, ex + extra))
    print('%d evaluations to run, %d graphs missing' % (len(jobs), len(missing)), file=sys.stderr)
    for m, rid in missing[:20]:
        print('  missing: %s %s' % (m, rid), file=sys.stderr)
    bad = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for m, rid, rc, sec, last in ex.map(run, jobs):
            bad += rc != 0
            print('%s\t%s\trc=%d\t%.0fs\t%s' % (m, rid, rc, sec, last), flush=True)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
