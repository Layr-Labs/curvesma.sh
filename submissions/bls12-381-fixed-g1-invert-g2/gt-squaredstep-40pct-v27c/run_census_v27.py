#!/usr/bin/env python3
"""Driver: the DP-protocol census at four report depths, one JSON, one merge."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v27_measure as M

PLAN = {1: 60000, 2: 60000, 3: 30000, 4: 30000}


def main():
    out = {}
    tot_n = tot_rec = 0
    for d, n in sorted(PLAN.items()):
        r = M.dp_census(trials=n, dp_bits=d, ratio=1.0, workers=24)
        out[str(d)] = r
        tot_n += r['trials']
        tot_rec += r['recovered']
        print('dp=%d n=%d S=%d reports/trial %.2f stored/trial %.2f collide %.4f '
              'recover %.4f CPlo %.4f model %.4f'
              % (d, r['trials'], r['budget'], r['reports_per_trial'],
                 r['stored_per_trial'], r['p_collide_per_trial'], r['p_recover'],
                 r['cp95_lower'], r['model_success']), flush=True)
    out['pooled'] = dict(trials=tot_n, recovered=tot_rec, p_recover=tot_rec / tot_n,
                         cp95_lower=float(M.cp95_lower(tot_rec, tot_n)),
                         note='pooled over four report depths; NOT a production rate')
    out['production_like'] = dict(
        depth=1, trials=out['1']['trials'], recovered=out['1']['recovered'],
        p_recover=out['1']['p_recover'], cp95_lower=out['1']['cp95_lower'],
        S_over_theta=out['1']['budget'] / out['1']['theta'],
        note='largest realizable S/theta in the fixture; production S/theta = 2^45.8')
    out['plan'] = PLAN
    with open('evidence/v27-measure-census.json', 'w') as f:
        json.dump(out, f, sort_keys=True, default=str, separators=(',', ':'))
    print('WROTE evidence/v27-measure-census.json pooled n=%d p=%.4f CPlo=%.4f'
          % (tot_n, tot_rec / tot_n, out['pooled']['cp95_lower']))


if __name__ == '__main__':
    main()
