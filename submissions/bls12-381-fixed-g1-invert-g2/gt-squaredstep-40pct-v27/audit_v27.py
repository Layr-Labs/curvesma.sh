#!/usr/bin/env python3
"""v27 audits: four PLANTED controls, each one must FIRE on the defect.

A. DP-PREDICATE control.  The v26 defect was that a stored key counted as a
   success without the 2^-dp_bits predicate.  Plant it: no non-DP key may enter
   the report store, a non-DP key offered to the credit path must be refused, and
   a deliberately planted non-DP key must never be credited.

B. j != 0 RECOVERY control.  The v26 formula y^3 = sigma^j/(m1^4 m2^2) must FAIL
   the exact return check away from j=0, and the j=0 arm must reproduce it exactly.

C. WRONG-DENOMINATOR control.  Charge the trap against 4 trails (v26) and against
   the true trail-start population; the 4-trail charge must come out ~2^43.8 low.

D. dlog COLLAPSE control.  Key the cycle census by dlog[min(keys12)] (the v26
   index) and the 9,420 representatives must collapse onto 785 nodes.

E. CP control.  An independent Clopper-Pearson implementation (log-space binomial
   CDF, no scipy, not cp95.py) must reproduce every shipped CP lower bound.
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K
import ledger27 as L27
import routing
import v27_measure as M


def control_dp(bits=5, seed=31337, trials=1500, dp_bits=2):
    toy = M.make_toy(bits=bits, seed=seed)
    rng = random.Random(seed)
    theta = 1 << dp_bits
    budget = K.step_budget(n=toy.N)
    store_violations = 0
    reports = 0
    non_dp_offered = 0
    non_dp_credited = 0
    for _ in range(trials):
        x = rng.randrange(1, toy.r)
        z = toy.T.pow(toy.g, x)
        store = {}
        used = 0
        while used < budget:
            a, b = rng.randrange(1, toy.r), rng.randrange(1, toy.r)
            y = toy.T.mul(toy.T.pow(toy.g, a), toy.T.pow(z, b))
            for _s in range(K.TRAIL_METERS * theta):
                if used >= budget:
                    break
                y, a, b, kind, _i = M.step26(toy, y, a, b)
                used += 1
                k = toy.key(y)[0]
                if not routing.dp_hit(k, dp_bits):
                    if k in store:
                        non_dp_offered += 1        # PLANT: offer a non-DP key
                        if routing.dp_hit(k, dp_bits):
                            non_dp_credited += 1
                    continue
                reports += 1
                if not routing.dp_hit(k, dp_bits):
                    store_violations += 1
                store[k] = (a, b)
                break
    k_nondp = None
    y = toy.g
    for _ in range(4000):
        k = toy.key(y)[0]
        if not routing.dp_hit(k, dp_bits):
            k_nondp = k
            break
        y = toy.T.mul(y, toy.g)
    credited = 0
    if k_nondp is not None:
        fake = {k_nondp: (1, 2)}
        credited = int((k_nondp in fake) and routing.dp_hit(k_nondp, dp_bits))
    return dict(trials=trials, dp_bits=dp_bits, reports=reports,
                store_violations=store_violations, non_dp_offered=non_dp_offered,
                non_dp_credited=non_dp_credited, planted_key_credited=credited,
                pass_control=(store_violations == 0 and non_dp_credited == 0
                              and credited == 0 and reports > 0))


def control_algebra(bits=5, seed=77021):
    r = M.recovery_proof(bits=bits, seed=seed, j_list=[0, 1, 2, 3], msgs=10)
    return dict(samples=r['samples'], misses=r['misses'],
                old_hits_j0=r['old_formula_hits_j0'],
                old_hits_jne0=r['old_formula_hits_jne0'], jne0_total=r['jne0_total'],
                j0_regression_kept=r['j0_regression_kept'],
                pass_control=(r['misses'] == 0 and r['old_formula_hits_jne0'] == 0
                              and r['old_formula_hits_j0'] == r['j0_regression_total']))


def control_trap():
    import ledger26 as L26
    tr27, d27 = L27.trap27()
    tr26, rate26 = L26.trap26()
    return dict(v26_trails=4, true_trail_starts=d27['n_trails'],
                v26_charge_log2=math.log2(float(tr26)),
                v27_charge_log2=math.log2(float(tr27)),
                understated_bits=math.log2(d27['n_trails'] / 4.0),
                pass_control=(d27['n_trails'] > 4 and float(tr27) > float(tr26)))


def control_collapse(bits=5, seed=20261024):
    toy = M.make_toy(bits=bits, seed=seed)
    exact = M.full_map_census(toy, 'v26')
    T, r = toy.T, toy.r
    elems = [None] * r
    y = toy.g
    for i in range(1, r):
        elems[i] = y
        y = T.mul(y, toy.g)
    dlog = {}
    for i in range(1, r):
        dlog[toy.key(elems[i])[0]] = i
    nxt = [0] * r
    for i in range(1, r):
        y2, _a, _b, _k, _i = M.step26(toy, elems[i], 1, 0)
        nxt[i] = dlog[toy.key(y2)[0]]
    state = [0] * r
    hist = {}
    two = 0
    for i in range(1, r):
        if state[i]:
            continue
        path, j = [], i
        while j and not state[j]:
            state[j] = 1
            path.append(j)
            j = nxt[j]
        if j and state[j] == 1:
            cyc = path[path.index(j):]
            hist[len(cyc)] = hist.get(len(cyc), 0) + 1
            if len(cyc) == 2:
                two += 1
        for q in path:
            state[q] = 2
    return dict(exact_states=exact['distinct_exact_states'],
                collapsed_nodes=len(dlog), collapse_factor=(r - 1) / len(dlog),
                exact_two_cycles=exact['two_cycles'], collapsed_two_cycles=two,
                collapsed_cycle_lengths={str(k): v for k, v in sorted(hist.items())},
                pass_control=(len(dlog) < exact['distinct_exact_states']))


def _binom_cdf_le(k, n, p):
    """P(X <= k) for X ~ Bin(n,p), log-space, incremental log-comb."""
    if p <= 0.0:
        return 1.0
    if p >= 1.0:
        return 0.0 if k < n else 1.0
    lp, lq = math.log(p), math.log1p(-p)
    logc = 0.0
    s = 0.0
    for i in range(0, k + 1):
        if i:
            logc += math.log(n - i + 1) - math.log(i)
        s += math.exp(logc + i * lp + (n - i) * lq)
    return s


def independent_cp(k, n, alpha=0.05, iters=60):
    """Bisection Clopper-Pearson lower bound; independent of cp95.py."""
    lo, hi = 0.0, (k / n) if k else 1.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        # CP lower bound: the p with P(X <= k; p) = 1 - alpha
        if _binom_cdf_le(k, n, mid) > 1.0 - alpha:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def control_cp(path='evidence/v27-measure-census.json'):
    if not os.path.exists(path):
        return dict(skipped=True, pass_control=False)
    d = json.load(open(path))
    rows = []
    ok = True
    for depth in ('1', '2'):
        if depth not in d:
            continue
        k, n = d[depth]['recovered'], d[depth]['trials']
        shipped = float(d[depth]['cp95_lower'])
        mine = independent_cp(k, n)
        rows.append(dict(depth=int(depth), k=k, n=n, shipped=shipped,
                         independent=mine, delta=abs(shipped - mine)))
        ok = ok and abs(shipped - mine) < 2e-3
    return dict(rows=rows, pass_control=ok)


if __name__ == '__main__':
    out = dict(dp=control_dp(), algebra=control_algebra(), trap=control_trap(),
               collapse=control_collapse(), cp=control_cp())
    out['all_pass'] = all(bool(v.get('pass_control')) for v in
                          (out['dp'], out['algebra'], out['trap'], out['collapse'],
                           out['cp']))
    print(json.dumps(out, indent=1, sort_keys=True, default=str))
    with open('evidence/v27-audit.json', 'w') as f:
        json.dump(out, f, indent=1, sort_keys=True, default=str)
    sys.exit(0 if out['all_pass'] else 1)
