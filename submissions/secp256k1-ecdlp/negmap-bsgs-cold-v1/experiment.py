"""Finite-scope checks of the signed-window (negation-map) BSGS in proof.md.

Evidence only at the stated small orders. It does not execute production
secp256k1 arithmetic, allocate the exponential table, time anything, or
establish the production cost bound, which is analytic (proof.md).
"""
import json
import sys


def plan(n, m):
    # Stride s=2m-1, offset c=m-1; window i covers exactly [i*s, (i+1)*s-1].
    s = 2 * m - 1
    c = m - 1
    h = n // s  # largest h with h*s <= n, so every covered scalar is distinct
    return s, c, h


def additive(n, m):
    """Z/nZ: class key of v is min(v, n-v), the analogue of an x-coordinate."""
    s, c, h = plan(n, m)
    key = lambda v: min(v % n, (-v) % n)
    table = {}
    for j in range(m):
        k = key(j)
        assert k not in table, "duplicate key would break the class argument"
        table[k] = j
    covered = h * s
    recovered = 0
    for x in range(n):
        out = None
        for i in range(h):
            t = (x - i * s - c) % n
            k = key(t)
            if k in table:
                j = table[k]
                for a in ((i * s + c + j) % n, (i * s + c - j) % n):
                    if a == x:  # independent verification [a]P=Q
                        out = a
                break  # first key hit stops the scan
        if x < covered:
            assert out == x, ("covered scalar not recovered", n, m, x)
            recovered += 1
        else:
            assert out is None, ("uncovered scalar returned a value", n, m, x)
    return {"order": n, "m": m, "stride": s, "h": h,
            "covered": covered, "recovered": recovered,
            "exact_coverage": recovered == covered}


def curve(p, m):
    """Real j=0 curve y^2=x^3+7 over F_p of prime order; key = x-coordinate."""
    def add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        if P[0] == Q[0] and (P[1] + Q[1]) % p == 0:
            return None
        if P == Q:
            l = 3 * P[0] * P[0] * pow(2 * P[1], p - 2, p) % p
        else:
            l = (Q[1] - P[1]) * pow(Q[0] - P[0], p - 2, p) % p
        x3 = (l * l - P[0] - Q[0]) % p
        return (x3, (l * (P[0] - x3) - P[1]) % p)

    def mul(k, P):
        R = None
        while k:
            if k & 1:
                R = add(R, P)
            P = add(P, P)
            k >>= 1
        return R

    pts = [(x, y) for x in range(p) for y in range(p) if (y * y - x ** 3 - 7) % p == 0]
    n = len(pts) + 1
    if any(n % d == 0 for d in range(2, int(n ** 0.5) + 1)):
        return None  # skip non-prime order
    G = pts[0]
    key = lambda R: ("inf",) if R is None else ("x", R[0])
    s, c, h = plan(n, m)
    table = {}
    R = None
    for j in range(m):
        k = key(R)
        assert k not in table, "two distinct classes shared an x-coordinate"
        table[k] = j
        R = add(R, G)
    covered = h * s
    recovered = 0
    for x in range(n):
        Q = mul(x, G)
        out = None
        for i in range(h):
            gamma = add(Q, mul((-(i * s + c)) % n, G))
            k = key(gamma)
            if k in table:
                j = table[k]
                for a in ((i * s + c + j) % n, (i * s + c - j) % n):
                    if mul(a, G) == Q:
                        out = a
                break
        if x < covered:
            assert out == x, ("curve: covered scalar not recovered", p, x)
            recovered += 1
        else:
            assert out is None or mul(out, G) == Q
    return {"field": p, "group_order": n, "m": m, "stride": s, "h": h,
            "covered": covered, "recovered": recovered,
            "exact_coverage": recovered == covered}


def main():
    cfg = json.load(sys.stdin)
    mode = cfg.get("mode", "all")
    res = {"scope": "toy groups only; does not validate production cost bounds"}
    if mode in ("all", "additive"):
        runs = []
        for n in (101, 257, 521, 1009):
            for m in (2, 4, 8, 16):
                if 2 * m - 1 <= n:
                    runs.append(additive(n, m))
        res["additive"] = runs
        res["additive_all_exact"] = all(x["exact_coverage"] for x in runs)
    if mode in ("all", "curve"):
        runs = []
        for p in (43, 67, 79, 97, 103, 109, 127, 139, 151):
            if p % 3 == 1:
                out = curve(p, 4)
                if out:
                    runs.append(out)
        res["curve"] = runs
        res["curve_all_exact"] = all(x["exact_coverage"] for x in runs)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
