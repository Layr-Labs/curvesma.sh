#!/usr/bin/env python3
"""v31 EXTENDED RETENTION-TAIL ARM.  Route (b): the merge simulator is run over
the WHOLE production age range, not the old 16-theta fixture support.

The v30 tail row evaluated an asserted exponential envelope at 2^-3.589688: a
0.083-class constant charged into a claim whose error budget is 2^-100-class, and
the envelope was extrapolated far past the fitted support.  v31 replaces it with
data:

  * the fixture trail cap is raised from 16*theta to 128*theta == the PRODUCTION
    cap, and the toy budget from 32 to 8*32 = 256 steps, so raw merge ages are
    sampled on [1, 128] theta with no extrapolation inside the production range;
  * every raw age is recorded in a 1-step histogram (no bins), and the survival
    counts n(A > a) are the evidence the fit is made on;
  * the residual is only the part of the production range where no age was
    observed: (a_hi, 128 theta], closed-form integrated against an exponential
    whose rate is FITTED to the measured tail decades (per-decade slopes and
    least-squares residuals are both shipped) and whose anchor is a rule-of-three
    upper bound, not a measured point.

No solve/verify gate here: the twelve-shift gate literals come from v30's arm
(v30-age-law.json, 41426 merges, 0.999976 verified) and are unchanged in v31.
Skipping the 12-modexp verify keeps this arm cheap enough to run long.
"""
import json
import math
import multiprocessing as mp
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K
import classkey as CK
import routing
import toy_census as TC

BITS = 5
THETA = 2                          # dp_bits = 1 -> theta = 2 steps
CAP_MULT = 128                     # PRODUCTION cap, not the 16-theta fixture
BUDGET = 8 * K.step_budget(n=785)  # 256 steps = 8x the old support
WORKERS = 12
PER_WORKER = 2000


def batch(args):
    bits, seed, trials, dp_bits, cap_mult, budget = args
    toy = TC.Toy(bits=bits, seed=seed)
    rng = random.Random(seed ^ 0x51A7)
    theta = 1 << dp_bits
    cap = cap_mult * theta
    rec = dict(trials=trials, dp_bits=dp_bits, theta=theta, cap=cap,
               budget=budget, n_trails=0, merges=0, censored=0,
               age_hist={}, trail_len_hist={}, sum_l=0, sum_l2=0, max_age=0,
               steps=0, max_trail_len=0)
    for _ in range(trials):
        x = rng.randrange(1, toy.r)
        z = toy.T.pow(toy.g, x)
        store = {}
        used = 0
        end = 'budget'
        L = 0
        while used < budget:
            a = rng.randrange(1, toy.r)
            b = rng.randrange(1, toy.r)
            y = toy.T.mul(toy.T.pow(toy.g, a), toy.T.pow(z, b))
            L = 0
            end = 'budget'
            for _s in range(cap):
                if used >= budget:
                    end = 'budget'
                    break
                y, a, b, kind, _i = CK.step_v26(toy.T, y, a, b, toy.TBL, toy.U,
                                                toy.r, toy.tbits, toy.rows)
                used += 1
                L += 1
                k = toy.key(y)[0]
                is_dp = routing.dp_hit(k, dp_bits)
                if k in store:
                    # v30 semantics: an internal merge does NOT end the trail; the
                    # age is recorded and the trail keeps walking to its report.
                    rec['merges'] += 1
                    rec['age_hist'][str(L)] = rec['age_hist'].get(str(L), 0) + 1
                    rec['max_age'] = max(rec['max_age'], L)
                    continue
                store[k] = 1
                if is_dp:
                    end = 'report'
                    break
            rec['n_trails'] += 1
            rec['sum_l'] += L
            rec['sum_l2'] += L * L
            rec['steps'] += L
            rec['max_trail_len'] = max(rec['max_trail_len'], L)
            rec['trail_len_hist'][str(L)] = rec['trail_len_hist'].get(str(L), 0) + 1
            if end == 'budget':
                rec['censored'] += 1
    return rec


def merge(parts):
    out = dict(parts[0])
    for p in parts[1:]:
        for k, v in p.items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    out[k][kk] = out[k].get(kk, 0) + vv
            elif k in ('max_age', 'max_trail_len'):
                out[k] = max(out[k], v)
            elif isinstance(v, int) and k not in ('dp_bits', 'theta', 'cap', 'budget'):
                out[k] += v
    return out


def main():
    os.makedirs('evidence', exist_ok=True)
    jobs = [(BITS, 4400 + 13 * i, PER_WORKER, 1, CAP_MULT, BUDGET)
            for i in range(WORKERS)]
    with mp.Pool(WORKERS) as pool:
        recs = pool.map(batch, jobs)
    m = merge(recs)
    surv = {str(k): sum(v for a, v in m['age_hist'].items() if int(a) > k * m['theta'])
            for k in range(1, CAP_MULT + 1)}
    out = dict(rule='raw merge ages over the FULL production range [1,128] theta '
                    '(cap_mult=128 == production cap; budget 8x the old support); '
                    '1-step histogram, no bins, no solve/verify in this arm',
               theta=m['theta'], cap=m['cap'], cap_mult=CAP_MULT,
               budget=m['budget'], trials=m['trials'], n_trails=m['n_trails'],
               merges=m['merges'], censored=m['censored'],
               max_age=m['max_age'], max_age_over_theta=m['max_age'] / m['theta'],
               mean_trail_len=m['sum_l'] / m['n_trails'],
               sum_l=str(m['sum_l']), sum_l2=str(m['sum_l2']), steps=m['steps'],
               age_hist=m['age_hist'], trail_len_hist=m['trail_len_hist'],
               survival_gt_k_theta=surv,
               workers=WORKERS, per_worker=PER_WORKER)
    json.dump(out, open('evidence/v31-tail-raw.json', 'w'), indent=1,
              sort_keys=True)
    print('trials %d trails %d merges %d censored %d max_age %d (%.1f theta) '
          'mean_trail %.3f' % (m['trials'], m['n_trails'], m['merges'],
                               m['censored'], m['max_age'],
                               m['max_age'] / m['theta'],
                               m['sum_l'] / m['n_trails']))


if __name__ == '__main__':
    main()
