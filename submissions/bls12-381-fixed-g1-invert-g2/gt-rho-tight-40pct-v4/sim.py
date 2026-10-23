#!/usr/bin/env python3
"""Exponent-domain simulation of the same walk at the production constants.

The walk is the production rule reduced to its exponent structure: a state in
Z/M, an injective step (x -> 2x, the exponent image of a squaring, injective
because the class fold makes squaring injective and 2 is invertible mod the odd
modulus) and a table step (x -> x + t_j).  The step type and the table index come
from routing.fields(), i.e. the same bit layout the production step uses.

The budget S is computed by constants.step_budget() -- the single source of the
walk constants -- with the modulus M in place of N, so the simulation is the
production budget at scale, not a restatement of it.
"""
import json, random, sys
from math import log2, sqrt, exp
import constants as K
import routing


def odd_modulus(q, seed):
    """An odd modulus near 2^q (2 must be invertible for the injective step)."""
    rng = random.Random(seed)
    m = (1 << q) | 1
    while True:
        m += 2
        if all(m % p for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)):
            if pow(2, m-1, m) == 1:      # Fermat witness, cheap
                return m
    return m


def make_map(M, R, seed):
    rng = random.Random(seed)
    tbl = [rng.randrange(1, M) for _ in range(R)]
    def step(x):
        t, tw, dp, st = routing.fields(x)
        if st < 3:
            return (2*x) % M
        return (x + tbl[t % R]) % M
    return step


def first_repeat(M, R, limit, seed):
    """Steps until a state repeats; also returns whether a DP was hit."""
    step = make_map(M, R, seed)
    seen = set()
    x = seed % M or 1
    for i in range(1, limit+1):
        seen.add(x)
        x = step(x)
        if x in seen:
            return i, True
    return limit, False


def run(q, R, groups, seed, mu=K.MU_WALK, f=K.FRAC_SQR):
    import importlib
    # budget from the same function, at scale M and with the same constants
    st = K.step_budget(mu=mu, target=K.TARGET, f=f, n=(1 << q))
    M = odd_modulus(q, seed)
    limit = 3*int(st)
    tot = 0; within = 0; used = 0
    for g in range(groups):
        tc, tracked = first_repeat(M, R, limit, seed*1000003 + g)
        tot += tc; used += 1
        if tc <= st:
            within += 1
    return dict(q=q, M=M, R=R, groups=groups, S=int(st), mean_tc=tot/used,
                mean_over_sqrtM=(tot/used)/sqrt(M), p_within=within/used,
                mean_over_S=(tot/used)/st)


def pooled(rows):
    w = sum(r['groups'] for r in rows)
    p = sum(r['p_within']*r['groups'] for r in rows)/w
    var = sum(r['groups']*(r['p_within']-p)**2 for r in rows)/w
    return p, sqrt(var/w)


def ideal_at_budget(mu, f, target):
    return 1 - exp(-(1-f*f)*log2(1/(1-float(target)))*float(mu)/float(mu))


def main():
    rows = []
    for (q, R, groups) in ((12, 256, 1500), (13, 512, 1500), (14, 1024, 1500)):
        rows.append(run(q, R, groups, seed=q*7919))
    # control: no squaring steps at the same budget, to exhibit the injectivity effect
    ctrl = []
    for (q, R, groups) in ((12, 256, 1000), (13, 512, 1000)):
        ctrl.append(run(q, R, groups, seed=q*104729, f=K.Fraction(0, 1)))
    p, se = pooled(rows)
    pc, sec = pooled(ctrl)
    model = 1 - exp(-(1-float(K.FRAC_SQR)**2)*log2(1/(1-float(K.TARGET)))/float(K.MU_WALK))
    ideal = 1 - exp(-(1-float(K.FRAC_SQR)**2)*log2(1/(1-float(K.TARGET))))
    ratio = 0
    for r in rows:
        pass
    out = dict(rows=rows, control=ctrl, pooled_p=p, pooled_se=se,
               control_p=pc, control_se=sec, model_at_mu=model, ideal_at_mu_one=ideal,
               requirement=float(K.SUCCESS), f=float(K.FRAC_SQR), mu=float(K.MU_WALK),
               target=float(K.TARGET),
               mean_tc_over_sqrtM=sum(r['mean_over_sqrtM']*r['groups'] for r in rows)/sum(r['groups'] for r in rows),
               control_mean_tc_over_sqrtM=sum(r['mean_over_sqrtM']*r['groups'] for r in ctrl)/sum(r['groups'] for r in ctrl))
    out['injectivity_ratio'] = out['mean_tc_over_sqrtM']/out['control_mean_tc_over_sqrtM']
    tot = sum(r['groups'] for r in rows)
    out['pooled_se_binomial'] = sqrt(p*(1-p)/tot)
    out['model_rows'] = [1 - exp(-r['S']*r['S']*(1-float(K.FRAC_SQR)**2)/(2*float(K.MU_WALK)*r['M'])) for r in rows]
    out['model_pooled'] = sum(m*r['groups'] for m, r in zip(out['model_rows'], rows))/tot
    out['predicted_ratio'] = 1/sqrt(1-float(K.FRAC_SQR)**2)
    out['pass'] = bool(p >= float(K.SUCCESS))
    json.dump(out, open("sim-out.json", "w"), indent=1)
    for r in rows:
        print("q=%2d R=%4d M=%7d S=%5d mean Tc/sqrt(M)=%.4f P(Tc<=S)=%.4f"
              % (r['q'], r['R'], r['M'], r['S'], r['mean_over_sqrtM'], r['p_within']))
    for r in ctrl:
        print("control f=0 q=%2d R=%4d mean Tc/sqrt(M)=%.4f P=%.4f"
              % (r['q'], r['R'], r['mean_over_sqrtM'], r['p_within']))
    print("pooled P = %.4f +- %.4f (binomial)   model at integer S = %.4f   requirement %.2f"
          % (p, out['pooled_se_binomial'], out['model_pooled'], float(K.SUCCESS)))
    print("injectivity ratio measured %.4f vs predicted %.4f"
          % (out['injectivity_ratio'], out['predicted_ratio']))
    print("SIM", "PASS" if out['pass'] else "FAIL")
    return 0 if out['pass'] else 1


if __name__ == "__main__":
    sys.exit(main())
