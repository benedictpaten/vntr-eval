#!/usr/bin/env python3
"""replicate_noise.py -- each graph's own Stage 3 noise: how far its call moves under neutral changes.

For every graph and region, the Stage 3 call (call_local.py FINAL_ARM, hybrid200k) is compared with
the same graph's REPLICATE calls, made with call_local.py's replicate arms:
  flank50k  hybrid50k      the hybrid graph with a 50 kb instead of a 200 kb genome-wide flank, the
                           window's reads fetched and re-mapped to it (other competitors, other
                           minimizer statistics, other depth-rate windows at the window's edges)
  ids       hybrid200kids  the Stage 3 hybrid graph with the span's node IDs re-laid (half a spread
                           step on; mc's native IDs spread too), re-mapped: the same graph and the
                           same alignments (giraffe is ID-invariant), other depth-rate windows
  mcids     hybrid200kmcids  the same with the span's nodes on mc's own native span IDs, i.e. at mc's
                           ID density (not for mc, whose layout this is; not for mc_relabel, for
                           which it IS the mc call; impossible for a graph with more span nodes
                           than mc)
For mc, 'ids' moves its native layout to a spread one; for the others 'ids' shifts a spread layout
by half a step and 'mcids' moves it to mc's native density -- so every graph gets a numbering
change of the kind that moves mc (native density <-> spread).
Neither changes the graph's alignment or the reads' source, so any change of the called haplotypes
is noise of the pipeline, not a property of the graph.

The distance between two calls is the pair edit distance of their span haplotypes, built by the
Stage 3 scoring rule (score_haplotypes.Prepared: same records, same overlap and nesting rules), best
pairing; a replicate's phase blocks are oriented to the Stage 3 call's pair, never to the truth. So
the noise uses NO TRUTH: it can be measured for candidates whose scores nobody has seen.

  noise(graph, region) = max over the graph's replicates of d(replicate, Stage 3 call)

By the triangle inequality every measured replicate's ED to the truth lies within ED +- noise, so a
paired difference larger than noise(label) + noise(baseline) cannot be closed by any measured
replicate of either graph (score_haplotypes.py summarise --noise). The noise is a lower bound on the
real spread: two perturbations, one draw each.

With --truth-labels, the replicates of those graphs are also scored against the truth (ED per
replicate); by default only graphs whose Stage 3 scores are already known (mc and the null graphs,
the truth graph), never a candidate.

    python3 tools/replicate_noise.py --graphs mc,unit_aware__all,mafft_linsi,unit_aware,truth,mc_relabel,mc_unchop \
        --regions pilot --out results/stage3_noise

-> <out>.tsv (one row per graph and region), <out>.md (per graph and stratum), and
   work/stage3/noise/noise.json (every distance, the replicate VCFs used, flags)
"""
import argparse
import collections
import json
import os
import statistics
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import call_local as cl  # noqa: E402
import config  # noqa: E402
import score_haplotypes as sh  # noqa: E402

REPLICATES = collections.OrderedDict([('flank50k', 'hybrid50k'), ('ids', 'hybrid200kids'),
                                      ('mcids', 'hybrid200kmcids')])
# replicates an arm cannot give: mc_relabel's 'ids' arm would be mc spread with the half-step offset,
# i.e. exactly mc's own 'ids' replicate; 'mcids' is mc's own layout (for mc_relabel: the mc call)
NOT_DEFINED = {('mc_relabel', 'ids'), ('mc', 'mcids'), ('mc_relabel', 'mcids')}
KNOWN_SCORES = ('mc', 'mc_relabel', 'mc_unchop', 'truth')
GRAPHS = ['mc', 'unit_aware__all', 'mafft_linsi', 'unit_aware', 'truth', 'mc_relabel', 'mc_unchop']


def haplotypes(region, vcf, target=None):
    """The span haplotype pair a VCF writes under the scoring rule; phase blocks in the VCF's own
    orientation, or oriented to minimise the pair ED to `target` (another call's pair)."""
    recs, _, _ = sh.read_vcf_records(vcf, region['contig'], region['a1'], region['b1'])
    ploidy = sh.truth_ploidy(region['contig'], region['a1'], region['b1'])
    prep = sh.Prepared(recs, region['a1'], region['b1'], region['ref'], ploidy)
    if target is None or prep.n_blocks <= 1:
        seqs = prep.build([False] * prep.n_blocks)[0]
        return seqs, prep.n_blocks
    ph = sh.best_phase(prep, sh.EdCache(target))
    return ph['seqs'], prep.n_blocks


def pair_distance(a, b):
    return min(sh.edit_distance(a[0], b[0]) + sh.edit_distance(a[1], b[1]),
               sh.edit_distance(a[0], b[1]) + sh.edit_distance(a[1], b[0]))


def truth_ed(region, vcf):
    recs, _, _ = sh.read_vcf_records(vcf, region['contig'], region['a1'], region['b1'])
    return sh.score_haplotypes(region, recs)['ed']


def one(graph, rid, truth_labels):
    region = sh.load_region(os.path.join(config.REGIONS_DIR, rid))
    row = {'graph': graph, 'region_id': rid, 'stratum': region['stratum'], 'flags': []}
    main = cl.calls_path(graph, rid)
    if not os.path.exists(main):
        row['status'] = 'no Stage 3 call'
        return row
    ref_pair, nb = haplotypes(region, main)
    row['vcf'] = main
    row['blocks'] = nb
    if nb > 1:
        row['flags'].append('stage3_call_has_%d_phase_blocks' % nb)
    score = graph in truth_labels
    if score:
        row['ed'] = truth_ed(region, main)
    dists = {}
    for name, arm in REPLICATES.items():
        if (graph, name) in NOT_DEFINED:
            row['d_' + name] = None
            row['why_' + name] = 'not defined (see NOT_DEFINED)'
            continue
        v = cl.calls_path(graph, rid, arm)
        if not os.path.exists(v):
            row['d_' + name] = None
            row['why_' + name] = 'no replicate call'
            continue
        pair, nbr = haplotypes(region, v, target=ref_pair)
        d = 0 if tuple(pair) == tuple(ref_pair) else pair_distance(pair, ref_pair)
        dists[name] = d
        row['d_' + name] = d
        row['vcf_' + name] = v
        if score:
            row['ed_' + name] = truth_ed(region, v)
    row['n_replicates'] = len(dists)
    row['noise'] = max(dists.values()) if dists else None
    row['status'] = 'ok'
    return row


def markdown(rows, graphs, out_tsv):
    L = ['# Stage 3 replicate noise: how far each graph\'s own call moves under neutral changes (pilot)', '',
         'Written by `tools/replicate_noise.py` (the method is in its docstring) into `%s`. For each graph and '
         'region, d is the pair edit distance between the span haplotypes of the Stage 3 call (hybrid 200 kb, '
         'call_local.py) and of a replicate call of the same graph: `flank50k` (50 kb genome-wide flank, reads '
         're-fetched and re-mapped), `ids` (the Stage 3 hybrid graph with its span node IDs re-laid by half a '
         'spread step; for mc, its native IDs spread) and `mcids` (the span nodes on mc\'s native span IDs; not '
         'for mc or mc_relabel), each re-mapped. Noise = the largest. No truth is used, so candidate graphs are '
         'measured blind.'
         % os.path.relpath(out_tsv, config.REPO), '']
    by = collections.defaultdict(list)
    for r in rows:
        if r.get('status') == 'ok':
            by[(r['graph'], r['stratum'])].append(r)
    L += ['| graph | stratum | n | summed noise | median | max (region) | regions noise > 0 / > 10 | '
          'summed d flank50k / ids / mcids |', '|---|---|---|---|---|---|---|---|']
    for g in graphs:
        for st in sh.STRATA_ORDER:
            rs = by.get((g, st))
            if not rs:
                continue
            nz = [r['noise'] for r in rs if r['noise'] is not None]
            if not nz:
                continue
            mx = max(rs, key=lambda r: r['noise'] or -1)
            def sd(k):
                v = [r[k] for r in rs if r.get(k) is not None]
                return '%d (%d)' % (sum(v), len(v)) if v else '-'
            L.append('| %s | %s | %d | %d | %s | %d (%s) | %d / %d | %s / %s / %s |' % (
                g, st, len(nz), sum(nz), sh._fmt(statistics.median(nz) * 1.0), mx['noise'], mx['region_id'],
                sum(1 for x in nz if x > 0), sum(1 for x in nz if x > 10), sd('d_flank50k'), sd('d_ids'),
                sd('d_mcids')))
    L += ['', '## Per region', '', 'd flank50k / d ids / d mcids per graph; "-" where the replicate does not exist '
          '(summed columns above: sum (regions)).', '',
          '| region | stratum | ' + ' | '.join(graphs) + ' |', '|' + '---|' * (len(graphs) + 2)]
    cell = {(r['graph'], r['region_id']): r for r in rows}
    for rid in sh.PILOT if set(r['region_id'] for r in rows) <= set(sh.PILOT) else sorted({r['region_id'] for r in rows}):
        rr = [cell.get((g, rid)) for g in graphs]
        st = next((r['stratum'] for r in rr if r), '')

        def f(r):
            if not r or r.get('status') != 'ok':
                return '-'
            return '%s / %s / %s' % (sh._fmt(r.get('d_flank50k')), sh._fmt(r.get('d_ids')), sh._fmt(r.get('d_mcids')))
        L.append('| %s | %s | %s |' % (rid, st, ' | '.join(f(r) for r in rr)))
    known = [r for r in rows if r.get('status') == 'ok' and 'ed' in r]
    if known:
        L += ['', '## Replicate ED to the truth (graphs whose Stage 3 scores are already known only)', '',
              '| region | graph | ED Stage 3 call | ED flank50k | ED ids | ED mcids |', '|---|---|---|---|---|---|']
        for r in known:
            if not any(r.get(k) for k in ('ed', 'ed_flank50k', 'ed_ids', 'ed_mcids')):
                continue
            L.append('| %s | %s | %s | %s | %s | %s |' % (r['region_id'], r['graph'], r['ed'],
                                                          sh._fmt(r.get('ed_flank50k')), sh._fmt(r.get('ed_ids')),
                                                          sh._fmt(r.get('ed_mcids'))))
        L += ['', 'Rows where every ED is 0 are left out.']
    return '\n'.join(L) + '\n'


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--graphs', default=','.join(GRAPHS))
    ap.add_argument('--regions', default='pilot')
    ap.add_argument('--truth-labels', default=','.join(KNOWN_SCORES),
                    help='graphs whose replicates are also scored against the truth (default: never a candidate)')
    ap.add_argument('--out', default=os.path.join(config.RESULTS_DIR, 'stage3_noise'))
    ap.add_argument('--jobs', type=int, default=4)
    a = ap.parse_args(argv)
    graphs = [g for g in a.graphs.split(',') if g]
    truth_labels = set(x for x in a.truth_labels.split(',') if x)
    candidates = set(graphs) - set(KNOWN_SCORES)
    if truth_labels & candidates:
        print('warning: scoring candidate replicates against the truth: %s' % sorted(truth_labels & candidates),
              file=sys.stderr)
    ids = sh.select_regions(a.regions)
    import concurrent.futures
    tasks = [(g, rid) for g in graphs for rid in ids]
    with concurrent.futures.ProcessPoolExecutor(a.jobs) as ex:
        rows = list(ex.map(one, [t[0] for t in tasks], [t[1] for t in tasks], [truth_labels] * len(tasks)))
    cols = ['graph', 'region_id', 'stratum', 'status', 'n_replicates', 'd_flank50k', 'd_ids', 'd_mcids', 'noise',
            'blocks', 'ed', 'ed_flank50k', 'ed_ids', 'ed_mcids', 'why_mcids']
    with open(a.out + '.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join('' if r.get(c) is None else str(r.get(c)) for c in cols) + '\n')
    with open(a.out + '.md', 'w') as f:
        f.write(markdown(rows, graphs, a.out + '.tsv'))
    od = os.path.join(cl.STAGE3, 'noise')
    os.makedirs(od, exist_ok=True)
    with open(os.path.join(od, 'noise.json'), 'w') as f:
        json.dump(rows, f, indent=1)
    for r in rows:
        print(r['graph'], r['region_id'], r.get('status'), 'd50k', r.get('d_flank50k'), 'dids', r.get('d_ids'),
              'dmcids', r.get('d_mcids'), 'noise', r.get('noise'),
              ('ED %s/%s/%s/%s' % (r.get('ed'), r.get('ed_flank50k'), r.get('ed_ids'), r.get('ed_mcids')))
              if 'ed' in r else '')


if __name__ == '__main__':
    main()
