import csv, collections, random, sys
rows = list(csv.DictReader(open('region_metrics.tsv'), delimiter='\t'))
by = collections.defaultdict(dict)
for r in rows: by[r['region']][r['arm']] = int(r['fp']) + int(r['fn'])
random.seed(1)
pairs = [p.split(':') for p in sys.argv[1].split(',')]
for excl in ((), ('TR756034',)):
    print('excluding', excl or 'nothing')
    for a, ref in pairs:
        d = [by[r][a] - by[r][ref] for r in by if r not in excl]
        n = len(d); bs = sorted(sum(random.choice(d) for _ in range(n)) for _ in range(2000))
        b = sum(x < 0 for x in d); w = sum(x > 0 for x in d)
        print('  %-13s - %-13s %+5d  [%+d, %+d]  better/worse %d/%d  totals %d vs %d' % (a, ref, sum(d), bs[50], bs[1949], b, w,
              sum(by[r][a] for r in by if r not in excl), sum(by[r][ref] for r in by if r not in excl)))
