"""4r: split the graph-implied cost/opt of each variant's hap32 pairs into pairs whose two alleles fall in one
bbt cluster (K = 32 or 64 partition, as bbt_plan) and pairs across clusters. Run from the repository root:
    python3 work/iterate/bbt_decomp.py"""
import json, os, subprocess, sys, concurrent.futures, collections
REPO = os.path.expanduser('~/PycharmProjects/vntr-eval')
sys.path.insert(0, os.path.join(REPO, 'tools'))
os.environ.setdefault('VNTR_REGIONS', os.path.join(REPO, 'work/stage4/regions'))
import iterate as it, panel, poa_panel, realign, msa_graph
OUT = os.path.join(REPO, 'work', 'iterate', 'bbt_decomp')     # per-pair costs (evaluate.py, VNTR_EVAL_PAIRS_OUT)
VS = ['poa_abpoa', 'poa_abpoa__all', 'mst', 'bbt32', 'bbt64', 'bbt32m', 'bbt64m']

def pairs(v, rid):
    o = os.path.join(OUT, 'pairs', v, rid + '.json')
    if not os.path.exists(o):
        os.makedirs(os.path.dirname(o), exist_ok=True)
        subprocess.run([sys.executable, os.path.join(REPO, 'tools/evaluate.py'), os.path.join(it.REGIONS, rid),
                        it.cand_gfa(v, rid), '--name', 'dc_' + v, '--out', o + '.eval.json', '--skip', 'truth',
                        '--threads', '2'], env=dict(os.environ, VNTR_EVAL_PAIRS_OUT=o), stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, cwd=REPO)
    return json.load(open(o))

def clusters(rid, K):
    fa, mp = poa_panel.ensure_union(os.path.join(it.REGIONS, rid), rid)
    rows = panel.read_map(mp)
    seqs = [s for _, s in msa_graph.read_fasta(fa)]
    recs = [('s%d' % k, realign.mask_runs(s.upper(), 1)[0]) for k, s in enumerate(seqs)]
    pan = [i for i, r in enumerate(rows) if int(r['weight']) > 0]
    g = [i for i in pan if any(x.startswith('GRCh38#') for x in rows[i]['members'])]
    reps, steps, _ = it.bbt_plan(recs, pan, g[0] if g else None, K)
    root = {r: r for r in reps}
    for x, p in steps:
        root[x] = root[p]
    name2c = {}
    for i, r in enumerate(rows):
        for m in r['members']:
            name2c[m] = root.get(i)
    return name2c

def one(rid):
    out = {}
    for K in (32, 64):
        cl = clusters(rid, K)
        for v in VS:
            acc = collections.Counter()
            for a, b, c, o in pairs(v, rid):
                t = 'in' if cl.get(a) is not None and cl.get(a) == cl.get(b) else 'x'
                acc[t + '_c'] += c; acc[t + '_o'] += o; acc[t + '_n'] += 1
            out['%s@%d' % (v, K)] = dict(acc)
    return rid, out

def main():
    ids = it.region_ids()
    res = {}
    with concurrent.futures.ProcessPoolExecutor(5) as ex:
        for rid, o in ex.map(one, ids):
            res[rid] = o
    json.dump(res, open(OUT + '.json', 'w'), indent=1)
    for K in (32, 64):
        print('partition K=%d' % K)
        for v in VS:
            t = collections.Counter()
            for rid in ids:
                t.update(res[rid]['%s@%d' % (v, K)])
            print('  %-15s within: %5d pairs cost/opt %.3f (opt %d)   across: %5d pairs cost/opt %.3f (opt %d)' % (
                v, t['in_n'], t['in_c'] / max(t['in_o'], 1), t['in_o'], t['x_n'], t['x_c'] / max(t['x_o'], 1), t['x_o']))


if __name__ == '__main__':
    main()
