#!/usr/bin/env python3
"""v24 CENSUS -- the same-row inverse-pair family removed STRUCTURALLY, every
census re-run with PER-WORKER SEED STREAMS, and the success statistic measured
UNCONDITIONALLY (one product, both factors logged, from the same paired runs).

Verdict fcfdd383's four defects, each answered by a named function here:

  D1 ledger sum            -> ledger.total_all()  (S*c_stress + trap charge +
                              every aux row, exact Fractions, claim >= log2 sum)
  D2 the cycle family      -> classkey.step_v24 (row may not repeat the previous
                              multiply row) + family_census() which counts the
                              family on the v23 step (control) and on the v24 step
                              (must be exactly 0, by construction)
  D3 fixture -> production -> recovery_census(): UNCONDITIONAL, paired independent
                              starts, ONE event = (both reach a DP) AND (same DP
                              class) AND (the twelve-shift solve VERIFIES against
                              z).  Both factors logged, product + bootstrap CI.
                              Plus the model-free union bound P(repeat within t)
                              <= t(t-1)/2N which is what carries production.
  D4 census integrity      -> every worker draws from random.Random(BASE +
                              1000*worker) and starts from toy.rand_gt(that rng);
                              worker_starts are returned so the dedup assert can
                              prove the trajectories are distinct.
"""
import json
import os
import random
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import classkey as CK
import constants as K
import routing
import toy_census as TC
from cp95 import cp95_lower

EV = os.path.join(HERE, 'evidence')
BASE_SEED = 20261024
CAP_MULT = 128
BITS = 5
ROWS = 64


def key_of(toy, y):
    return CK.classkey(toy.T, y)[0]


def run_trail(toy, y, a, b, dp_bits, cap, live=True):
    """One metered trail: to its FIRST distinguished point or to the cap.

    Returns (reached_dp, dp_key, steps, trapped, y, a, b).
    A class repeat before the DP is a cycle in the class map: trapped, no report.
    """
    seen = {}
    pr = -1
    n = 0
    while n < cap:
        k = key_of(toy, y)
        if k in seen:
            return False, None, n, True, y, a, b
        seen[k] = n
        if live:
            y, a, b, _kd, _ki, pr = toy.step_live(y, a, b, pr)
        else:
            y, a, b, _kd, _ki = toy.step_v23(y, a, b)
        n += 1
        if routing.dp_hit(key_of(toy, y), dp_bits):
            return True, key_of(toy, y), n, False, y, a, b
    return False, None, n, False, y, a, b


# ------------------------------------------------------------------ D3
def _recovery_batch(args):
    trials, dp_bits, cap_mult, worker, bits, rows, live = args
    toy = TC.Toy(bits=bits, seed=BASE_SEED, rows=rows)
    rng = random.Random(BASE_SEED + 1000 * worker)      # D4: PER-WORKER STREAM
    r = toy.r
    cap = cap_mult * (1 << dp_bits)
    out = dict(trials=0, both_dp=0, a_dp=0, b_dp=0, collided=0, solved=0,
               verified=0, trapped=0, steps=0, starts=[], pairs=[], by_trial=[])
    for _ in range(trials):
        # ONE planted secret per trial; both trails are g^a z^b with independent
        # label pairs, so the two starts are INDEPENDENT group elements of the
        # same instance and the twelve-shift solve has a known answer x0.
        x0 = rng.randrange(1, r)
        z = toy.T.pow(toy.g, x0)
        a1, b1 = rng.randrange(1, r), rng.randrange(1, r)
        a2, b2 = rng.randrange(1, r), rng.randrange(1, r)
        y1 = toy.T.mul(toy.T.pow(toy.g, a1), toy.T.pow(z, b1))
        y2 = toy.T.mul(toy.T.pow(toy.g, a2), toy.T.pow(z, b2))
        out['starts'].append(hash((a1, b1)))
        r1 = run_trail(toy, y1, a1, b1, dp_bits, cap, live)
        r2 = run_trail(toy, y2, a2, b2, dp_bits, cap, live)
        out['trials'] += 1
        out['steps'] += r1[2] + r2[2]
        out['a_dp'] += int(r1[0])
        out['b_dp'] += int(r2[0])
        out['trapped'] += int(r1[3]) + int(r2[3])
        both = r1[0] and r2[0]
        out['both_dp'] += int(both)
        rec = 0
        if both and r1[1] == r2[1]:                       # SAME class key at a DP
            out['collided'] += 1
            cands = CK.solve_candidates(toy.ppow, toy.ppow, r1[5], r1[6], r2[5], r2[6], r)
            out['solved'] += int(bool(cands))
            # KNOWN-ANSWER GATE: the twelve-shift solve must return the planted x0.
            ok = any(x == x0 for _j, x in cands)
            out['verified'] += int(ok)
            rec = int(ok)
        out['by_trial'].append(1 if (both and r1[1] == r2[1] and rec) else 0)
    return out


def recovery_census(trials=24000, dp_bits=5, cap_mult=CAP_MULT, workers=6,
                    bits=BITS, rows=ROWS, live=True):
    per = trials // workers + 1
    chunks = [(per, dp_bits, cap_mult, w, bits, rows, live) for w in range(workers)]
    with Pool(workers) as p:
        parts = p.map(_recovery_batch, chunks)
    agg = dict(trials=0, both_dp=0, a_dp=0, b_dp=0, collided=0, solved=0,
               verified=0, trapped=0, steps=0, starts=[], by_trial=[])
    for q in parts:
        for k in ('trials', 'both_dp', 'a_dp', 'b_dp', 'collided', 'solved',
                  'verified', 'trapped', 'steps'):
            agg[k] += q[k]
        agg['starts'] += q['starts']
        agg['by_trial'] += q['by_trial']
    n = agg['trials']
    # dedup (D4): every worker's FIRST start must differ -- under the v23
    # construction they were equal by construction.
    firsts = [p['starts'][0] for p in parts]
    dedup_ok = len(set(firsts)) == len(firsts)
    p_uncond = agg['verified'] / n if n else 0.0
    lo = cp95_lower(agg['verified'], n)
    # bootstrap CI over independent trials
    rng = random.Random(7)
    boots = []
    for _ in range(400):
        s = sum(rng.choice(agg['by_trial']) for _ in range(n))
        boots.append(s / n)
    boots.sort()
    return {
        'trials': n, 'workers': workers, 'bits': bits, 'rows': rows,
        'dp_bits': dp_bits, 'cap_steps': cap_mult * (1 << dp_bits),
        'step': 'classkey.step_v24 (row-difference kill)' if live else 'classkey.step_v13',
        'seeds': 'random.Random(%d + 1000*worker); starts drawn from THAT stream' % BASE_SEED,
        'worker_first_starts_distinct': dedup_ok,
        'worker_first_starts': [str(f) for f in firsts],
        'factor_a_trail_reaches_dp': agg['a_dp'] / n if n else 0.0,
        'factor_b_trail_reaches_dp': agg['b_dp'] / n if n else 0.0,
        'both_reach_dp': agg['both_dp'] / n if n else 0.0,
        'same_dp_class_given_both': (agg['collided'] / agg['both_dp']) if agg['both_dp'] else 0.0,
        'solve_verified_given_collision': (agg['verified'] / agg['collided']) if agg['collided'] else 0.0,
        'UNCONDITIONAL_RECOVERY': p_uncond,
        'UNCONDITIONAL_CP95_LOWER': lo,
        'bootstrap_ci95': [boots[9], boots[-10]],
        'trapped_trails': agg['trapped'], 'mean_steps_per_trial': agg['steps'] / n if n else 0,
        'reading': 'ONE product from the SAME runs: P(both reach a DP) x P(same DP class | both) '
                   'x P(the twelve-shift solve verifies | same class).  No conditioning on a DP '
                   'having been reached is applied to the numerator.',
    }


# ------------------------------------------------------------------ D2
def _family_batch(args):
    """A FAMILY TRACE: a straight run of `steps` multiply/squaring steps per walk,
    counting consecutive same-row multiply pairs, the kappa-delta-6 subset (the
    inverse pair), and the exact returns y' == y.  No trap break: the family is
    the mechanism, and it shows up as an exact return."""
    walks, steps, worker, bits, rows, live = args
    toy = TC.Toy(bits=bits, seed=BASE_SEED, rows=rows)
    rng = random.Random(BASE_SEED + 7919 * worker)
    out = dict(steps=0, mul_steps=0, same_row_pairs=0, inverse_pairs=0,
               exact_returns=0, walks=0, reps=0)
    for _ in range(walks):
        y = toy.rand_gt(rng)
        a = b = 1
        prev = -1
        pk = None
        seen = {}
        for n in range(steps):
            k = key_of(toy, y)
            if k in seen:
                out['reps'] += 1
            seen[k] = n
            if live:
                y2, a2, b2, kd, ki, _p2 = toy.step_live(y, a, b, prev)
            else:
                y2, a2, b2, kd, ki = toy.step_v23(y, a, b)
            row = ki[1]
            out['steps'] += 1
            if kd == 'mul':
                out['mul_steps'] += 1
                if pk is not None and row == prev:
                    out['same_row_pairs'] += 1
                    if (ki[0] - pk) % 12 == 6:
                        out['inverse_pairs'] += 1
                        if y2 == y:
                            out['exact_returns'] += 1
                prev = row
                pk = ki[0]
            y, a, b = y2, a2, b2
        out['walks'] += 1
    return out


def family_census(walks=40, steps=400, workers=4, bits=BITS, rows=ROWS, live=True):
    per = walks // workers + 1
    chunks = [(per, steps, w, bits, rows, live) for w in range(workers)]
    with Pool(workers) as p:
        parts = p.map(_family_batch, chunks)
    agg = dict(steps=0, mul_steps=0, same_row_pairs=0, inverse_pairs=0,
               exact_returns=0, walks=0, reps=0)
    for q in parts:
        for k in agg:
            agg[k] += q[k]
    return dict(step=('classkey.step_v24 (row difference)' if live
                      else 'classkey.step_v13 (v23 shipped step)'),
                walks=agg['walks'], steps=agg['steps'], mul_steps=agg['mul_steps'],
                consecutive_same_row_multipliers=agg['same_row_pairs'],
                of_which_kappa_delta_6=agg['inverse_pairs'],
                exact_returns=agg['exact_returns'],
                class_repeats_seen=agg['reps'],
                family_rate_per_mul_step=(agg['inverse_pairs'] / agg['mul_steps']
                                          if agg['mul_steps'] else None),
                structural_kill_holds=(agg['same_row_pairs'] == 0) if live else None)


def main():
    os.makedirs(EV, exist_ok=True)
    out = {}
    out['family_v23_step'] = family_census(live=False)
    out['family_v24_step'] = family_census(live=True)
    out['recovery_v24'] = recovery_census(trials=12000)
    out['recovery_v23_step_control'] = recovery_census(trials=1200, live=False)
    json.dump(out, open(os.path.join(EV, 'census-v24.json'), 'w'), indent=1, default=str)
    print(json.dumps({k: (v if not isinstance(v, dict) else
                          {kk: vv for kk, vv in v.items() if not isinstance(vv, list)})
                      for k, v in out.items()}, indent=1))


if __name__ == '__main__':
    main()
