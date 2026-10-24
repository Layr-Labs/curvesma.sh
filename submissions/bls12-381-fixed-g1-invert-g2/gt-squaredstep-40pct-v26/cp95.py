#!/usr/bin/env python3
"""Clopper-Pearson one-sided lower bound, Decimal, bisection.

The FUNCTION BODY IS COPIED VERBATIM from
    curvesmash/submissions/secp256k1-ecdlp/sigma-table-40pct-v15k/scripts/gen_cert_v15k.py
(the calibrated helper that certifies 0.398285, NOT the v12c float helper that
returned 3.1e-61 for 9/10).  The ONLY change is the bisection bracket: the v15k
copy hard-codes [0.30, 0.45] around its own 0.408 census, and a bracket that does
not contain the answer silently returns a bracket endpoint.  Here the bracket is
data-driven and the function asserts it brackets the root.

Calibration (asserted in toy_census.py): 0/10 -> 0.0, 10/10 -> 0.7411,
9/10 -> 0.6056, 8/10 -> 0.4901, 4/10 -> 0.1684.
"""
from decimal import Decimal


def cp95_lower(k, n, alpha=Decimal('0.05')):
    """one-sided Clopper-Pearson lower bound: smallest p with P[X>=k | p,n] >= alpha.
    P[X>=k] is strictly INCREASING in p, so bisection converges on the lower end."""
    k = int(k); n = int(n)
    if k == 0:
        return Decimal(0)
    lnC = sum(Decimal(i).ln() for i in range(n - k + 1, n + 1)) - sum(Decimal(i).ln() for i in range(1, k + 1))

    def tail(p):
        p = Decimal(p); q = Decimal(1) - p
        t = (lnC + k * p.ln() + (n - k) * q.ln()).exp()   # term_k
        s = t
        for j in range(k, n):
            t *= (n - j) * p / ((j + 1) * q)
            s += t
            if t < Decimal('1e-40') * s:
                break
        return s

    lo, hi = Decimal(0), Decimal('0.999999999')   # p = 1 would divide by q = 0
    if tail(hi) <= alpha:                      # interval still cannot bracket
        return Decimal(1)
    for _ in range(90):
        mid = (lo + hi) / 2
        if tail(mid) > alpha:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal('1e-13'):
            break
    assert tail(lo) <= alpha, "bracket lost the root"
    return lo   # conservative: guaranteed below the true bound


if __name__ == '__main__':
    for k, n in ((0, 10), (4, 10), (8, 10), (9, 10), (10, 10), (100, 300)):
        print(k, n, cp95_lower(k, n))
