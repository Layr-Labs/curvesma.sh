#!/usr/bin/env python3
"""Supporting evidence for the interleaved-chain BSGS candidate (v3).

The construction in proof.md is: W independent affine chains advanced in lockstep,
one Montgomery batch inversion per round over the W denominators of that round, and
a deterministic balanced-search-tree lookup on the baby-step x-coordinate.

This file implements exactly that and checks, at toy scale only:

  1. the interleaved chains recover the scalar on every generated instance;
  2. the batched inversion used by the round returns exactly the per-element inverses;
  3. chain advancement handles the point at infinity, doubling and inverse points,
     including the chain that starts at [0]G = infinity.

It is not evidence about production cost, memory, timing or behaviour at 2^128 steps.
"""
import json
import random
import sys
from math import isqrt

# Chain states used by one lockstep round.
NORMAL = 0     # distinct x: ordinary incomplete addition
FROM_INF = 1   # current point is infinity: O + inc = inc
DOUBLING = 2   # current point equals inc: doubling formula
TO_INF = 3     # current point is minus inc: result is infinity
SKIP = 4       # index past the chain bound


def inv(a, p):
    return pow(a, -1, p)


def batch_inv(values, p):
    """Montgomery: one inversion, 3(n-1) multiplications for n denominators."""
    n = len(values)
    prefix = [1] * n
    acc = 1
    for i in range(n):
        prefix[i] = acc
        acc = acc * values[i] % p
    inv_acc = inv(acc, p)
    out = [0] * n
    for i in range(n - 1, -1, -1):
        out[i] = inv_acc * prefix[i] % p
        inv_acc = inv_acc * values[i] % p
    return out


def add(P, Q, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        s = (3 * x1 * x1) * inv(2 * y1, p) % p
    else:
        s = (y2 - y1) * inv(x2 - x1, p) % p
    x3 = (s * s - x1 - x2) % p
    return (x3, (s * (x1 - x3) - y1) % p)


def mul(k, P, p):
    R = None
    while k:
        if k & 1:
            R = add(R, P, p)
        P = add(P, P, p)
        k >>= 1
    return R


def curve_order(p):
    n = 1
    for x in range(p):
        rhs = (x * x * x + 7) % p
        if rhs == 0:
            n += 1
        elif pow(rhs, (p - 1) // 2, p) == 1:
            n += 2
    return n


def factor(n):
    fs, d = [], 2
    while d * d <= n:
        if n % d == 0:
            fs.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        fs.append(n)
    return fs


def order_of(P, n, p):
    d = n
    for f in factor(n):
        while d % f == 0 and mul(d // f, P, p) is None:
            d //= f
    return d


def find_generator(p, n):
    best, best_ord, tried = None, 0, 0
    for x in range(p):
        rhs = (x * x * x + 7) % p
        if rhs == 0:
            P = (x, 0)
        elif pow(rhs, (p - 1) // 2, p) == 1:
            if p % 4 == 3:
                y = pow(rhs, (p + 1) // 4, p)
            else:
                y = next((c for c in range(p) if c * c % p == rhs), None)
                if y is None:
                    continue
            P = (x, y)
        else:
            continue
        tried += 1
        d = order_of(P, n, p)
        if d > best_ord:
            best, best_ord = P, d
        if best_ord == n or tried >= 256:
            break
    return best


def interleaved_bsgs(Q, G, r, p, W, stats=None):
    """W interleaved affine chains, one batched inversion per round."""
    if Q is None:
        return 0
    m = isqrt(r) + 1
    h = -(-r // m)
    Kb = -(-m // W)
    Kg = -(-h // W)

    Wm = mul(W, G, p)          # [W]G   : fixed baby increment
    D = mul(m, G, p)           # [m]G
    WD = mul(W, D, p)          # [W][m]G : fixed giant increment
    negWD = None if WD is None else (WD[0], (-WD[1]) % p)

    table = {}

    def advance(points, inc, count, start, record):
        """One lockstep round over all W chains with a single batched inversion."""
        kinds, den = [], []
        for c in range(W):
            idx = start + c
            P = points[c]
            if idx >= count:
                kinds.append(SKIP); den.append(1); continue
            if P is None:
                # O + inc = inc: the chain starting at [0]G enters the walk here.
                kinds.append(FROM_INF); den.append(1); continue
            if inc is None:
                kinds.append(TO_INF); den.append(1); continue
            dx = (inc[0] - P[0]) % p
            if dx != 0:
                kinds.append(NORMAL); den.append(dx)
            elif P == inc:
                dy = (2 * P[1]) % p
                if dy == 0:
                    kinds.append(TO_INF); den.append(1)
                else:
                    kinds.append(DOUBLING); den.append(dy)
            else:
                kinds.append(TO_INF); den.append(1)
        invs = batch_inv(den, p)
        if stats is not None:
            stats["from_inf_steps"] = stats.get("from_inf_steps", 0) + kinds.count(FROM_INF)
            stats["doubling_steps"] = stats.get("doubling_steps", 0) + kinds.count(DOUBLING)
            stats["to_inf_steps"] = stats.get("to_inf_steps", 0) + kinds.count(TO_INF)
        out = []
        for c in range(W):
            idx = start + c
            P = points[c]
            kind = kinds[c]
            if kind == SKIP:
                out.append(P)
                continue
            if kind == FROM_INF:
                out.append(inc)
                continue
            if kind == TO_INF:
                out.append(None)
                continue
            if kind == NORMAL:
                lam = (inc[1] - P[1]) * invs[c] % p
            else:  # DOUBLING
                lam = (3 * P[0] * P[0]) * invs[c] % p
            x3 = (lam * lam - P[0] - inc[0]) % p
            y3 = (lam * (P[0] - x3) - P[1]) % p
            if record is not None:
                record(P[0], idx)
            out.append((x3, y3))
        return out

    # baby chains: chain c holds [c + k*W]G, increment +[W]G
    points = [mul(c, G, p) for c in range(W)]
    for k in range(Kb):
        points = advance(points, Wm, m, k * W, lambda x, j: table.__setitem__(x, j))
        if k == 0 and stats is not None:
            # chain 0 starts at [0]G = O; after one round it must hold [W]G.
            stats["chain0_after_round0_is_WG"] = (points[0] == Wm)

    # giant chains: chain c holds Q - [c + k*W][m]G, increment -[W][m]G
    gpoints = []
    for c in range(W):
        cD = mul(c, D, p)
        gpoints.append(Q if cD is None else add(Q, (cD[0], (-cD[1]) % p), p))
    for k in range(Kg):
        start = k * W
        for c in range(W):
            idx = start + c
            if idx >= h:
                continue
            P = gpoints[c]
            if P is None:
                # giant chain reached infinity: Q = [idx*m]P, hence b = 0
                cand = (idx * m) % r
                if mul(cand, G, p) == Q:
                    return cand
                continue
            if P[0] in table:
                j = table[P[0]]
                for cand in ((idx * m + j) % r, (idx * m - j) % r):
                    if mul(cand, G, p) == Q:
                        return cand
        gpoints = advance(gpoints, negWD, h, start, None)

    return None


def main():
    req = json.load(sys.stdin)
    p = int(req.get("prime", 10007))
    trials = int(req.get("trials", 16))
    W = int(req.get("width", 8))
    if not (7 <= p <= 100003):
        raise SystemExit("prime out of range")
    if not (2 <= W <= 64):
        raise SystemExit("width out of range")

    n = curve_order(p)
    G = find_generator(p, n)
    r = order_of(G, n, p)

    rng = random.Random(20261105)
    successes, checks = 0, []
    stats = {}
    for _ in range(trials):
        x = rng.randrange(r)
        Q = mul(x, G, p)
        got = interleaved_bsgs(Q, G, r, p, W, stats)
        ok = got is not None and mul(got, G, p) == Q
        successes += 1 if ok else 0
        checks.append({"x": x, "recovered": got, "ok": ok})

    vals = list(range(2, 66))
    batched = batch_inv(vals, p)
    batch_ok = all(batched[i] * vals[i] % p == 1 for i in range(len(vals)))

    json.dump({
        "prime": p,
        "curve": "y^2 = x^3 + 7",
        "group_order": n,
        "generator_order_r": r,
        "width_W": W,
        "trials": trials,
        "successes": successes,
        "success_rate": successes / trials,
        "batch_inversion_matches_pow": batch_ok,
        "chain_state_counts": {
            "infinity_start_steps": stats.get("from_inf_steps", 0),
            "doubling_steps": stats.get("doubling_steps", 0),
            "inverse_point_steps": stats.get("to_inf_steps", 0),
            "chain0_after_round0_is_WG": stats.get("chain0_after_round0_is_WG", False),
        },
        "construction": "interleaved affine chains with explicit infinity, doubling and inverse-point handling",
        "scope": "toy scale only; does not establish production cost or memory",
        "sample": checks[:4],
    }, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
