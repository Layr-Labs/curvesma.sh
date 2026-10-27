"""Emit the test vectors of the Python norm-bound reference (embedding_bounds.py) for `nbound12 check`, one 16-aligned block
per canonical pool row, from the regenerated pool (a, b from the per-index seeds of pool12-1.b85 + pool12-2.b85).  Run inside the verifier's sandbox or any
isolated Python: argv: package-dir out.txt [count]"""
import sys, json
from fractions import Fraction as F
sys.path.insert(0, sys.argv[1])
from embedding_bounds import root_disks, NormBounds, LOG_LOWER, LOG_UPPER, Interval
from common12 import integer_interval, integer_metric
from model12 import load_pool
from evidence12 import regenerate
cfg = json.load(open(sys.argv[1] + '/config.json'))
hdr, pool = load_pool(open(sys.argv[1] + '/pool12-1.b85').read() + open(sys.argv[1] + '/pool12-2.b85').read())
count = int(sys.argv[3]) if len(sys.argv) > 3 else len(pool)
rows = []
for r in pool[:count]:
    a, b, aq, bq = regenerate(cfg, r['seed'], r['index'])
    rows.append(dict(a=a, b=b, qa=F(aq) / integer_metric(cfg['matrix_a'])[1], seed=r['seed'], index=r['index']))
disks = root_disks()
out = []
places = [(r, i, rad) for (r, i), rad in disks if not i < -rad]
out.append('places %d' % len(places))
for r, i, rad in places:
    re = Interval.rational(r - rad, r + rad); im = Interval.rational(i - rad, i + rad)
    out.append('place %d %s %s %s %s' % (int(i == 0), re.lo.hex(), re.hi.hex(), im.lo.hex(), im.hi.hex()))
out.append('logtab ' + ' '.join(x.hex() for x in LOG_LOWER) + ' ' + ' '.join(x.hex() for x in LOG_UPPER))
qb, db, ib = integer_metric(cfg['matrix_b'])
step = 1
tot_w = tot_e = 0.0
for j, row in enumerate(rows[::step][:count]):
    a, b = row['a'], row['b']; radius = (cfg['production_A_squared'] - F(row['qa'])) * db
    suffix = b[1:]; t = sum(qb[0][i] * b[i] for i in range(1, 6))
    c = sum(qb[i][k] * b[i] * b[k] for i in range(1, 6) for k in range(1, 6)) - radius
    iv = integer_interval(qb[0][0], t, c)
    lo = max(iv.start, (b[0] // 16) * 16); hi = min(iv.stop - 1, (b[0] // 16) * 16 + 15)
    est = NormBounds(a, disks); leaves = est.adaptive_row(a, suffix, lo, hi, singleton_policy='reject')
    w = float((1 - F(row['qa']) / cfg['production_A_squared']) ** 3)
    tot_w += w; tot_e += w * est.counts['block_evaluations'] / (hi - lo + 1)
    out.append('row %d %d %d %.17g %d %d' % (iv.start, iv.stop - 1, b[0], w, row['seed'], row['index']))
    out.append('case %s %s %d %d %d %d' % (' '.join(map(str, a)), ' '.join(map(str, suffix)), lo, hi, est.counts['block_evaluations'], len(leaves)))
    for l, h, logs in leaves:
        if isinstance(logs, dict):
            out.append('leaf %d %d R' % (l, h))
        else:
            out.append('leaf %d %d %d %d %d %d' % (l, h, logs[0][0], logs[0][1], logs[1][0], logs[1][1]))
out.append('end weighted_evaluations_per_point %.12f' % (tot_e / tot_w))
open(sys.argv[2], 'w').write('\n'.join(out) + '\n')
print('cases', min(count, len(rows[::step])), 'weighted evals/point', tot_e / tot_w)
