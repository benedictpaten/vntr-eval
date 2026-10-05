#!/usr/bin/env python3
"""Interior sequence of every panel walk of star2/sampled.gbz over every region of test/patch.json.
Adapted from vntr-eval work/full/chr20/investigate/residual_sv/sampled.py (anchors by CHM13 coordinate:
the interior of a 1-based span s..e is the walk's sequence after ref base s-1 and before ref base e+1).
Output panel_star2.tsv.gz: region, walk name (sample#hap), seqstart of the W line, orient, seq."""
import collections, gzip, json, re, subprocess, sys

C = '/Users/benedictpaten/PycharmProjects/vntr-eval/work/full/chr20'
OUT = '/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-2026-07-vg-call-refactor/e061d8cf-7996-49a4-986b-3686374ef250/scratchpad/linkage'
VG = '/Users/benedictpaten/PycharmProjects/vntr-eval/work/bin/vg-91d38c802'
GBZ = C + '/test/banded/sampled.gbz'
STEP = re.compile(r'([<>])(\d+)')
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
A = {rid: tuple(v['span']) for rid, v in json.load(open(C + '/test/patch.json'))['regions'].items()}


def stream():
    p = subprocess.Popen([VG, 'convert', '-f', '-t', '1', GBZ], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for l in p.stdout:
        yield l
    p.wait()


seq = {}
ref_steps = None
for l in stream():
    if l[0] == 'S':
        _, n, s = l.rstrip('\n').split('\t')[:3]
        seq[int(n)] = s
    elif l[0] == 'W':
        x = l.split('\t', 7)
        if x[1] == 'CHM13' and x[3] == 'chr20':
            ref_steps = [(o, int(n)) for o, n in STEP.findall(x[6])]
            ref_start = int(x[4])
print('nodes', len(seq), 'ref steps', len(ref_steps), 'ref start', ref_start, file=sys.stderr, flush=True)

want = collections.defaultdict(list)
for r, (s, e) in A.items():
    want[s - 1].append((r, 'L'))
    want[e + 1].append((r, 'R'))
targets = sorted(want)
anch = collections.defaultdict(dict)
pos, ti = ref_start + 1, 0
for o, n in ref_steps:
    ln = len(seq[n])
    while ti < len(targets) and targets[ti] < pos + ln:
        t = targets[ti]
        if t >= pos:
            off = t - pos if o == '>' else ln - 1 - (t - pos)
            for r, side in want[t]:
                anch[r][side] = (n, off, o)
        ti += 1
    pos += ln
del ref_steps
bad = [r for r in A if 'L' not in anch[r] or 'R' not in anch[r] or anch[r]['L'][0] == anch[r]['R'][0]]
print('regions', len(A), 'without two distinct anchors', len(bad), bad[:5], file=sys.stderr, flush=True)

startnode = collections.defaultdict(list)
for r, d in anch.items():
    if r in bad:
        continue
    startnode[d['L'][0]].append(r)
    startnode[d['R'][0]].append(r)

out = gzip.open('/Users/benedictpaten/PycharmProjects/vntr-eval/work/full/chr20/repeat_sites/panel_banded.tsv.gz', 'wt', compresslevel=1)
nw = 0
for l in stream():
    if l[0] != 'W':
        continue
    x = l.split('\t', 7)
    name = '%s#%s' % (x[1], x[2])
    steps = STEP.findall(x[6])
    open_ = {}
    for i, (o, n) in enumerate(steps):
        n = int(n)
        for r in startnode.get(n, ()):
            Ln, Lo, _ = anch[r]['L']
            Rn, Ro, _ = anch[r]['R']
            if r not in open_:
                open_[r] = (i, n)
                continue
            j, m = open_.pop(r)
            if m == n:
                open_[r] = (i, n)
                continue
            sub = steps[j:i + 1]
            if m == Ln:
                fw, orient = sub, '+'
            else:
                fw, orient = [('<' if a == '>' else '>', b) for a, b in reversed(sub)], '-'
            parts = [seq[int(b)] if a == '>' else seq[int(b)].translate(COMP)[::-1] for a, b in fw]
            a0, a1 = fw[0][0], fw[-1][0]
            lcut = (Lo + 1) if a0 == '>' else (len(seq[Ln]) - Lo)
            rcut = Ro if a1 == '>' else (len(seq[Rn]) - 1 - Ro)
            s = ''.join(parts)
            interior = s[lcut:len(s) - (len(parts[-1]) - rcut)]
            out.write('%s\t%s\t%s\t%s\t%s\n' % (r, name, x[4], orient, interior))
    nw += 1
out.close()
print('walks', nw, file=sys.stderr)
