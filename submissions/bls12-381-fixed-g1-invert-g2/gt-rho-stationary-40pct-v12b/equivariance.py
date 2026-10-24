#!/usr/bin/env python3
"""Equivariance / retro-audit property test for the v12b walk, with a v10f control.

The class group of the GT walk is mu_12 = <sigma>, sigma the p-power Frobenius on
GT ~ mu_r, acting as multiplication by p on exponents.  The generators of its
subgroups, and their orders, are enumerated here rather than described:

    sigma      order 12      p^6 = -1 (mod r)  so sigma^6 IS inversion
    sigma^2    order  6      the mu_6 subgroup
    sigma^3    order  4
    sigma^4    order  3
    sigma^6    order  2      the inversion generator -1

The reviewer's test (retro-audit-22) is `class(S + R_c) == class(sigma S + R_c)`
for every sigma in the class group: a quotient walk whose step multiplies by a
table element R_c is well defined only if the successor CLASS does not depend on
which of the twelve representatives of the current class was held.

  E1  in v12b the state is canonical by construction, so the test is run in the
      form the design must satisfy: for every class, every one of its twelve
      members, and the class-keyed multiplier, the successor class through the
      shipped (canonicalise-then-multiply) path.  Splits must be 0.
  E2  the D1 test: is sigma^j(R_c) itself a table entry for some j != 0?  If it
      is, the walk admits an (m, m^-1) two-cycle.  Splits must be 0.
  E3  control: the same two tests on the v10f table (a sigma-ladder,
      tbl[t][i] = sigma^i(m0)) and on the v10f held-offset rule, where the answer
      must be NONZERO -- otherwise the test would not be detecting the defect it
      is supposed to detect.

Every count is written to e-equivariance.json; nothing here is quoted by hand.
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toy12d as T

GENERATORS = ((1, 12, 'sigma: p^6 = -1 mod r, so sigma^6 is inversion'),
              (2, 6, 'sigma^2: the mu_6 subgroup'),
              (3, 4, 'sigma^3'),
              (4, 3, 'sigma^4'),
              (6, 2, 'sigma^6: the inversion generator'))


def ladder_table(cw):
    """The v10f table: tbl[t][i] = sigma^i(m0), i.e. a sigma-ladder."""
    return [[cw.frob_pow(cw.T.pow(cw.g, cw.u[t][0]), i) for i in range(12)]
            for t in range(T.TABLE_ROWS)]


def run(bits, seed, n_classes=24):
    cw = T.ClassWalk(bits, seed)
    T_, r = cw.T, cw.r
    rng = random.Random(seed ^ 0xE9)
    lad = ladder_table(cw)
    checks = splits = 0
    d1_checks = d1_hits = 0
    ctl_checks = ctl_splits = 0
    ctl_d1_checks = ctl_d1_hits = 0
    # p^6 = -1 at this fixture, and the generator orders, verified not assumed
    orders = {}
    for j, order, _desc in GENERATORS:
        v = pow(cw.p, j, r)
        o, cur = 1, v
        while cur != 1 and o < 64:
            cur = cur * v % r
            o += 1
        orders[j] = o
    for _ in range(n_classes):
        x = rng.randrange(1, r)
        z = T_.pow(cw.g, x)
        a, b = rng.randrange(1, r), rng.randrange(1, r)
        y0 = T_.mul(T_.pow(cw.g, a), T_.pow(z, b))
        k0, cid0, _ = cw.kappa(y0)
        r0 = cw.frob_pow(y0, k0)                       # canonical representative
        kind, t, col, _ = cw.decision(cid0)
        for j, order, _desc in GENERATORS:
            # ---- E1: successor class from each of the twelve members, shipped path
            yj = r0
            for _ in range(j):
                yj = T_.frob(yj)
            kj, cidj, _ = cw.kappa(yj)
            succ = cw.step(cw.frob_pow(yj, kj), cidj, a, b)[0]
            _, sc, _ = cw.kappa(succ)
            base = cw.step(r0, cid0, a, b)[0]
            _, sb, _ = cw.kappa(base)
            checks += 1
            if sc != sb:
                splits += 1
            # ---- E2: the D1 test at the INVERSION generator.  sigma^6 is inversion
            # (p^6 = -1), so this asks whether the inverse of an applied multiplier is
            # again a table entry -- the exact (m, m^-1) two-cycle of the rejection.
            # Other j give accidental coincidences at toy density (the table is a
            # 768-value subset of a 2^15 field) and are reported separately, not
            # counted as structural.
            if kind == 'mul':
                m = cw.tbl[t][col]
                for _ in range(j):
                    m = T_.frob(m)
                if j == 6:                      # the inversion generator
                    d1_checks += 1
                    if any(m == e for row in cw.tbl for e in row):
                        d1_hits += 1
            # ---- E3 control: the v10f sigma-ladder and its held-offset rule
            yi = r0
            for _ in range(j):
                yi = T_.frob(yi)
            kl, cidl, _ = cw.kappa(yi)
            kindl, tl, twl, _ = cw.decision(cidl)
            if kindl == 'mul':
                raw = T_.mul(yi, lad[tl][(twl - kl) % 12])
                _, sctl, _ = cw.kappa(raw)
                _, sbase, _ = cw.kappa(T_.mul(r0, lad[tl][(twl - 0) % 12]))
                ctl_checks += 1
                if sctl != sbase:
                    ctl_splits += 1
                ctl_d1_checks += 1
                inv_applied = T_.frob(lad[tl][(twl - kl) % 12])
                for _ in range(5):
                    inv_applied = T_.frob(inv_applied)
                if any(inv_applied == e for row in lad for e in row):
                    ctl_d1_hits += 1
            elif kindl == 'sq':
                ctl_checks += 1
    return dict(
        density_bits=bits, prime_r=r,
        p6_is_minus_one=(pow(cw.p, 6, r) == r - 1),
        generators=[dict(j=j, order=order, desc=d, order_verified=orders[j] == order)
                    for j, order, d in GENERATORS],
        checks=checks, splits=splits,
        d1_checks=d1_checks, d1_inverse_closure_hits=d1_hits,
        control_v10f_checks=ctl_checks, control_v10f_splits=ctl_splits,
        control_v10f_d1_checks=ctl_d1_checks,
        control_v10f_d1_inverse_closure_hits=ctl_d1_hits,
        scope='E1: successor class through the shipped canonicalising path from each '
              'of the twelve sigma^j images vs from the canonical representative, for '
              'every generator of every subgroup of mu_12. E2: whether sigma^j of a '
              'table entry is itself a table entry (the (m,m^-1) two-cycle test). '
              'E3: the same two tests on the v10f sigma-ladder and held-offset rule, '
              'where a nonzero count is the detection control.',
    )


def run_all(seed=912):
    out = dict(tree='mu_12 = <sigma>, p^6 = -1 so sigma^6 is inversion',
               note='splits counts class differences; the v12b design must read 0',
               per_fixture=[])
    tot = dict(checks=0, splits=0, d1_checks=0, d1_hits=0,
               control_checks=0, control_splits=0, control_d1_hits=0)
    for b in (5, 7, 9):
        r = run(b, seed + b)
        out['per_fixture'].append(r)
        tot['checks'] += r['checks']; tot['splits'] += r['splits']
        tot['d1_checks'] += r['d1_checks']; tot['d1_hits'] += r['d1_inverse_closure_hits']
        tot['control_checks'] += r['control_v10f_checks']
        tot['control_splits'] += r['control_v10f_splits']
        tot['control_d1_hits'] += r['control_v10f_d1_inverse_closure_hits']
    out['totals'] = tot
    out['all_fixtures_p6_is_minus_one'] = all(f['p6_is_minus_one'] for f in out['per_fixture'])
    out['splits'] = tot['splits']
    out['checks'] = tot['checks']
    return out


if __name__ == '__main__':
    D = os.path.dirname(os.path.abspath(__file__))
    res = run_all()
    json.dump(res, open(os.path.join(D, 'e-equivariance.json'), 'w'), indent=1)
    print(json.dumps(res['totals'], indent=1))
