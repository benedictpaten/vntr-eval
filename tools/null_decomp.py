#!/usr/bin/env python3
"""null_decomp.py -- where does a Stage 3 null graph's call difference come from: the mapping or the caller?

The null graphs (call_local.py NULL_GRAPHS) carry MC's own alignment, written the way a candidate is:
  mc_relabel  mc.gfa with the span's node IDs spread as a candidate's are (order kept)
  mc_unchop   mc.gfa with every non-branching run of nodes merged
Each is re-mapped and called exactly like mc (work/stage3/diag/<null>@hybrid200k/). This tool takes one
null and, per region, against the Stage 3 mc call (work/stage3/calls/mc/):

  relabel ID...   1. checks the relabelled hybrid graph is the mc hybrid graph up to node IDs (nodes,
                     sequences, edges; the ID map comes from zipping the two hybrids' W-line walks);
                  2. compares the two GAFs read by read after translating mc's into the null's IDs;
                  3. calls each graph with the OTHER graph's alignments (translated);
                  4. calls both with --depth-term 0;
                  5. compares DR at records matched by their translated snarl ID, inside and outside
                     the span.
  unchop ID...    1. translates mc's GAF base-exactly onto the unchopped hybrid's nodes (the walks of
                     the two hybrids spell the same sequences, so every mc node is a set of pieces of
                     unchopped chunks) and compares it read by read with giraffe's own alignments;
                  2. calls the unchopped graph with mc's alignments;
                  3. calls both with --depth-term 0.

Haplotypes and ED come from the Stage 3 scoring rule (score_haplotypes.score_haplotypes), so every
ED here is comparable with results/stage3/<label>/<id>.json. (Before scorer version 2 this tool used
the gate's raw-rule builder, which shared the scorer's POS-order flaw: see score_haplotypes.Prepared.)
-> work/stage3/null/decomp/<null>/{<id>/..., decomp.json}

    python3 tools/null_decomp.py relabel L011138 L001909 L002013 L012272 L014297
    python3 tools/null_decomp.py unchop L009656 L014297 L005990 L012184 L002013 L011138 L012272 L001909
    python3 tools/null_decomp.py rescore-relabel      # recompute every ED of decomp.json from its kept VCFs
    python3 tools/null_decomp.py rescore-unchop       # (no vg run; alignment and VCF-identity checks kept)
"""
import argparse
import gzip
import json
import os
import re
import statistics
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import call_local as cl  # noqa: E402
import score_haplotypes as sh  # noqa: E402

STEP = re.compile(r'([<>])([^<>]+)')
FLIP = {'>': '<', '<': '>'}
H = os.path.join(cl.STAGE3, 'hybrid200k')
OUT = os.path.join(cl.STAGE3, 'null', 'decomp')


def load_gfa(path):
    """({node: seq length}, [((sample, hap, contig, start), [(orient, node)])], edges)"""
    lens, walks, edges = {}, [], set()
    with open(path) as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            if x[0] == 'S':
                lens[x[1]] = len(x[2])
            elif x[0] == 'L':
                edges.add((x[1], x[2], x[3], x[4]))
            elif x[0] == 'W':
                walks.append(((x[1], x[2], x[3], x[4]), STEP.findall(x[6])))
    return lens, walks, edges


def canon_edge(e):
    a, ao, b, bo = e
    r = (b, '-' if bo == '+' else '+', a, '-' if ao == '+' else '+')
    return min(e, r)


def id_map(gfa_a, gfa_b):
    """Node map a -> b from walks of identical names and step counts (the relabel null)."""
    la, wa, ea = load_gfa(gfa_a)
    lb, wb, eb = load_gfa(gfa_b)
    if [k for k, _ in wa] != [k for k, _ in wb]:
        raise SystemExit('walk names differ')
    m = {}
    for (_, sa), (_, sb) in zip(wa, wb):
        if len(sa) != len(sb):
            raise SystemExit('walk lengths differ')
        for (oa, na), (ob, nb) in zip(sa, sb):
            if oa != ob or m.setdefault(na, nb) != nb:
                raise SystemExit('inconsistent node map at %s' % na)
    iso = {'nodes_a': len(la), 'nodes_b': len(lb), 'unmapped_nodes': sum(1 for n in la if n not in m),
           'length_mismatch': sum(1 for n in la if n in m and lb.get(m[n]) != la[n]),
           'edges_only_a': len({canon_edge((m.get(a, a), ao, m.get(b, b), bo)) for a, ao, b, bo in ea}
                               - {canon_edge(e) for e in eb}),
           'edges_only_b': len({canon_edge(e) for e in eb}
                               - {canon_edge((m.get(a, a), ao, m.get(b, b), bo)) for a, ao, b, bo in ea})}
    return m, iso


def piece_map(gfa_a, gfa_b):
    """Each node of a as pieces of b's nodes, in a's forward orientation: [(b node, orient, start, end)]
    with start/end in b's forward coordinates. From walking both hybrids' W lines base by base."""
    la, wa, _ = load_gfa(gfa_a)
    lb, wb, _ = load_gfa(gfa_b)
    if [k for k, _ in wa] != [k for k, _ in wb]:
        raise SystemExit('walk names differ')
    m = {}
    for (_, sa), (_, sb) in zip(wa, wb):
        # b's steps as intervals on the walk
        ivs, q = [], 0
        for o, n in sb:
            ivs.append((q, q + lb[n], o, n))
            q += lb[n]
        j, p = 0, 0
        for om, n in sa:
            x0, x1 = p, p + la[n]
            while ivs[j][1] <= x0:
                j += 1
            pcs, k = [], j
            while k < len(ivs) and ivs[k][0] < x1:
                q0, q1, oc, c = ivs[k]
                lo, hi = max(x0, q0), min(x1, q1)
                M = q1 - q0
                cs, ce = (lo - q0, hi - q0) if oc == '>' else (M - (hi - q0), M - (lo - q0))
                pcs.append((c, oc, cs, ce))
                k += 1
            if om == '<':
                pcs = [(c, FLIP[o], cs, ce) for c, o, cs, ce in reversed(pcs)]
            if m.setdefault(n, pcs) != pcs:
                raise SystemExit('inconsistent pieces for node %s' % n)
            p = x1
        if p != q:
            raise SystemExit('walks spell different lengths')
    return m, lb


def translate_gaf(src, dst, step_fn):
    """Rewrite every GAF line's path, path length, start and end with step_fn(steps, start, end)."""
    bad = 0
    with open(src) as f, open(dst, 'w') as fo:
        for line in f:
            if line.startswith('@'):
                fo.write(line)
                continue
            x = line.split('\t')
            if x[5] != '*':
                t = step_fn(STEP.findall(x[5]), int(x[7]), int(x[8]))
                if t is None:
                    bad += 1
                    continue
                x[5], x[6], x[7], x[8] = t
            fo.write('\t'.join(x))
    return bad


def relabel_steps(m, lens):
    def fn(steps, a0, a1):
        st = [(o, m.get(n, n)) for o, n in steps]
        return ''.join(o + n for o, n in st), str(sum(lens[n] for _, n in st)), str(a0), str(a1)
    return fn


def unchop_steps(m, lb):
    def plen(st):
        return sum(lb[c] for _, c in st)

    def fn(steps, a0, a1):
        pcs = []
        for o, n in steps:
            x = m[n] if o == '>' else [(c, FLIP[oo], cs, ce) for c, oo, cs, ce in reversed(m[n])]
            for pc in x:
                if pcs and pcs[-1][:2] == pc[:2] and \
                        ((pc[1] == '>' and pcs[-1][3] == pc[2]) or (pc[1] == '<' and pcs[-1][2] == pc[3])):
                    c, oo, cs, ce = pcs[-1]
                    pcs[-1] = (c, oo, min(cs, pc[2]), max(ce, pc[3]))
                else:
                    pcs.append(pc)
        if any((cs, ce) != (0, lb[c]) for c, _, cs, ce in pcs[1:-1]):
            return None
        c, oo, cs, ce = pcs[0]
        lead = cs if oo == '>' else lb[c] - ce
        st = [(oo, c) for c, oo, _, _ in pcs]
        a0, a1 = a0 + lead, a1 + lead
        while len(st) > 1 and lb[st[0][1]] <= a0:      # a GAF path holds only the nodes it touches
            a0 -= lb[st[0][1]]
            a1 -= lb[st[0][1]]
            st = st[1:]
        while len(st) > 1 and plen(st[:-1]) >= a1:
            st = st[:-1]
        return ''.join(o + c for o, c in st), str(plen(st)), str(a0), str(a1)
    return fn


def gaf_records(path):
    d = {}
    with open(path) as f:
        for line in f:
            if line.startswith('@'):
                continue
            x = line.split('\t')
            mate = 'fn' if '\tfn:Z:' in line else 'fp' if '\tfp:Z:' in line else ''
            d[(x[0], mate, x[2])] = (x[5], x[7], x[8], x[11])
    return d


def compare_gafs(a, b, span_nodes):
    da, db = gaf_records(a), gaf_records(b)
    differ = [k for k in da if db.get(k) != da[k]]
    return {'reads': len(da), 'identical': len(da) - len(differ),
            'same_placement': sum(1 for k in da if k in db and db[k][:3] == da[k][:3]),
            'differ_touching_span': sum(1 for k in differ
                                        if any(n in span_nodes for _, n in STEP.findall(da[k][0])))}


def call_with(reg, graph_dir, gaf, out_prefix, extra=()):
    cmd = cl.vg_call_command(reg, os.path.join(graph_dir, 'call.gbz'), gaf, out_prefix + '.mosaic.tsv', 2, extra)
    cl.run(cmd, stdout=out_prefix + '.raw.vcf', stderr=out_prefix + '.log')
    with open(os.path.join(graph_dir, 'build.json')) as f:
        pj = json.load(f)
    cl.finish_vcf(out_prefix + '.raw.vcf', out_prefix + '.vcf.gz', reg, pj)
    return out_prefix + '.vcf.gz'


def vcf_body_equal(a, b):
    def body(p):
        with gzip.open(p, 'rt') as f:
            return [x for x in f if not x.startswith('#')]
    ba, bb = body(a), body(b)
    if not ba or not bb:
        raise SystemExit('empty VCF body: %s %s' % (a, b))
    return ba == bb


def haps(vcf, reg):
    """(ED to the truth, the sorted haplotype pair) under the Stage 3 scoring rule."""
    region = sh.load_region(reg['dir'])
    recs, _, _ = sh.read_vcf_records(vcf, region['contig'], region['a1'], region['b1'])
    res = sh.score_haplotypes(region, recs)
    return res['ed'], sorted(res['_seqs'])


def rescore_region(null, row):
    """Recompute every ED and haplotype identity of one decomp.json row from the VCFs run_region kept
    (the alignment comparisons and VCF-identity checks do not depend on the builder and are kept)."""
    rid = row['region_id']
    reg = cl.load_region(rid)
    od = os.path.join(OUT, null, rid)
    new = dict(row)
    ed_mc, h_mc = haps(cl.calls_path('mc', rid), reg)
    ed_null, h_null = haps(cl.calls_path(null, rid), reg)
    new.update({'ed_mc': ed_mc, 'ed_null': ed_null, 'haps_equal': h_mc == h_null})
    if null == 'mc_unchop':
        ed, h = haps(os.path.join(od, 'nullgraph_mcreads.vcf.gz'), reg)
        new.update({'ed_nullgraph_mcreads': ed, 'nullgraph_mcreads_haps_equal_mc': h == h_mc,
                    'nullgraph_mcreads_haps_equal_null': h == h_null})
    nd = {g: haps(os.path.join(od, g + '.nodepth.vcf.gz'), reg) for g in ('mc', null)}
    new['nodepth_ed_mc'], new['nodepth_ed_null'] = nd['mc'][0], nd[null][0]
    new['nodepth_haps_equal'] = nd['mc'][1] == nd[null][1]
    new['scorer_version'] = sh.SCORER_VERSION
    new['previous_eds'] = {k: row.get(k) for k in ('ed_mc', 'ed_null', 'ed_nullgraph_mcreads', 'nodepth_ed_mc',
                                                   'nodepth_ed_null') if k in row}
    return new


def dr_ratios(reg, vcf_a, vcf_b, m):
    """DR(b) / DR(a) at records matched by snarl ID (a's IDs translated), inside / outside the span."""
    def recs(p, mm):
        d = {}
        with gzip.open(p, 'rt') as f:
            for line in f:
                if line[0] == '#':
                    continue
                x = line.rstrip('\n').split('\t')
                key = ''.join(o + mm.get(n, n) for o, n in STEP.findall(x[2]))
                fm = dict(zip(x[8].split(':'), x[9].split(':')))
                d.setdefault(key, []).append((int(x[1]), fm.get('GT'), fm.get('DR')))
        return d
    a, b = recs(vcf_a, m), recs(vcf_b, {})
    out = {'span': [], 'flank': []}
    gt_diff = 0
    for k, va in a.items():
        vb = b.get(k)
        if not vb or len(va) != 1 or len(vb) != 1 or va[0][2] in (None, '.') or vb[0][2] in (None, '.'):
            continue
        if float(va[0][2]) == 0:
            continue
        inside = reg['span_start'] <= va[0][0] <= reg['span_end']
        out['span' if inside else 'flank'].append(float(vb[0][2]) / float(va[0][2]))
        gt_diff += inside and va[0][1] != vb[0][1]

    def q(v):
        return [round(min(v), 3), round(statistics.median(v), 3), round(max(v), 3), len(v)] if v else None
    return {'span_min_median_max_n': q(out['span']), 'flank_min_median_max_n': q(out['flank']),
            'span_matched_gt_differs': gt_diff}


def span_node_set(reg, hybrid_gfa_mc, hybrid_gfa_null):
    mc_span = set(load_gfa(os.path.join(reg['dir'], 'mc.gfa'))[0])
    flank = set(load_gfa(hybrid_gfa_mc)[0]) - mc_span
    return set(load_gfa(hybrid_gfa_null)[0]) - flank


def run_region(null, rid):
    reg = cl.load_region(rid)
    a, b = os.path.join(H, 'mc', rid), os.path.join(H, null, rid)
    ga, gb = os.path.join(a, 'hybrid.gfa'), os.path.join(b, 'hybrid.gfa')
    od = os.path.join(OUT, null, rid)
    os.makedirs(od, exist_ok=True)
    mc_vcf = cl.calls_path('mc', rid)
    null_vcf = cl.calls_path(null, rid)
    ed_mc, h_mc = haps(mc_vcf, reg)
    ed_null, h_null = haps(null_vcf, reg)
    row = {'region_id': rid, 'null': null, 'ed_mc': ed_mc, 'ed_null': ed_null, 'haps_equal': h_mc == h_null,
           'scorer_version': sh.SCORER_VERSION}
    lb = load_gfa(gb)[0]
    if null == 'mc_relabel':
        m, iso = id_map(ga, gb)
        row['isomorphism'] = iso
        inv = {v: k for k, v in m.items()}
        la = load_gfa(ga)[0]
        t = os.path.join(od, 'mcreads_translated.gaf')
        row['untranslatable_mc_reads'] = translate_gaf(os.path.join(a, 'reads.gaf'), t, relabel_steps(m, lb))
        row['gaf'] = compare_gafs(t, os.path.join(b, 'reads.gaf'), span_node_set(reg, ga, gb))
        v = call_with(reg, b, t, os.path.join(od, 'nullgraph_mcreads'))
        row['nullgraph_mcreads_equals_null_vcf'] = vcf_body_equal(v, null_vcf)
        os.remove(t)
        t = os.path.join(od, 'nullreads_translated.gaf')
        translate_gaf(os.path.join(b, 'reads.gaf'), t, relabel_steps(inv, la))
        v = call_with(reg, a, t, os.path.join(od, 'mcgraph_nullreads'))
        row['mcgraph_nullreads_equals_mc_vcf'] = vcf_body_equal(v, mc_vcf)
        os.remove(t)
        row['dr_null_over_mc'] = dr_ratios(reg, mc_vcf, null_vcf, m)
    else:
        m, _ = piece_map(ga, gb)
        t = os.path.join(od, 'mcreads_translated.gaf')
        row['untranslatable_mc_reads'] = translate_gaf(os.path.join(a, 'reads.gaf'), t, unchop_steps(m, lb))
        row['gaf'] = compare_gafs(t, os.path.join(b, 'reads.gaf'), span_node_set(reg, ga, gb))
        v = call_with(reg, b, t, os.path.join(od, 'nullgraph_mcreads'))
        ed, h = haps(v, reg)
        row['ed_nullgraph_mcreads'] = ed
        row['nullgraph_mcreads_haps_equal_mc'] = h == h_mc
        row['nullgraph_mcreads_haps_equal_null'] = h == h_null
        os.remove(t)
    nd = {}
    for g, gd in (('mc', a), (null, b)):
        v = call_with(reg, gd, os.path.join(gd, 'reads.gaf'), os.path.join(od, g + '.nodepth'), ['--depth-term', '0'])
        nd[g] = haps(v, reg)
    row['nodepth_ed_mc'], row['nodepth_ed_null'] = nd['mc'][0], nd[null][0]
    row['nodepth_haps_equal'] = nd['mc'][1] == nd[null][1]
    return row


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('kind', choices=('relabel', 'unchop', 'rescore-relabel', 'rescore-unchop'))
    ap.add_argument('ids', nargs='*')
    a = ap.parse_args(argv)
    if a.kind.startswith('rescore-'):
        null = 'mc_' + a.kind.split('-')[1]
        path = os.path.join(OUT, null, 'decomp.json')
        with open(path) as f:
            rows = [rescore_region(null, r) for r in json.load(f)]
        for r in rows:
            print(json.dumps({k: r.get(k) for k in ('region_id', 'ed_mc', 'ed_null', 'ed_nullgraph_mcreads',
                                                     'nodepth_ed_mc', 'nodepth_ed_null', 'previous_eds')}))
        with open(path, 'w') as f:
            json.dump(rows, f, indent=1)
        return
    null = 'mc_' + a.kind
    rows = []
    for rid in a.ids:
        r = run_region(null, rid)
        rows.append(r)
        print(json.dumps(r))
    with open(os.path.join(OUT, null, 'decomp.json'), 'w') as f:
        json.dump(rows, f, indent=1)


if __name__ == '__main__':
    main()
