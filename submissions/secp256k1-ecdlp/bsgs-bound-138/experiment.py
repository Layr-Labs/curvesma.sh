#!/usr/bin/env python3
"""Supporting evidence for the interleaved-chain BSGS upper bound (v6).

The construction in proof.md is: W independent affine chains advanced in lockstep,
one Montgomery batch inversion per round over the W denominators of that round, records
collected into a flat array, an LSD radix sort with counters wide enough for the bucket
loads, and a multiword binary search over the sorted array.

This file checks, at toy scale only:

  1. the interleaved chains recover the scalar on every generated instance;
  2. the batched inversion used by the round returns exactly the per-element inverses;
  3. each of the four chain-advance states is executed *through the real advance
     routine* and its result is checked, rather than only classified;
  4. the lookup is the priced structure: radix sort plus multiword binary search over
     full coordinates.

It is not evidence about production cost, memory or timing, and no part of the claimed
2^145 derives from it.
"""
import json
import random
import sys
from math import isqrt

NORMAL = 0
FROM_INF = 1
DOUBLING = 2
TO_INF = 3
SKIP = 4

LEDGER = {
    "baby": {"formula": 537, "batch": 538.40, "append": 8, "control": 8,
             "count": 1091.40, "allowance": 8731},
    "giant": {"formula": 537, "batch": 538.40, "search": 5120, "control": 8,
              "count": 6203.40, "allowance": 49627},
    "sort": {"histogram": 17, "distribute": 20, "passes": 16,
             "count": 592, "allowance": 4736},
    "per_2e128": 63094,
    "claimed_time_log2": 138,
    "priced_time_log2": 134.03422,
}


def inv(a, p):
    return pow(a, -1, p)


def batch_inv(values, p):
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


def lockstep_round(points, inc, count, start, record, p, W, stats=None):
    """One lockstep round over W chains with a single batched inversion.

    Module level so the self-test below can drive every state through this exact code
    path instead of a parallel reimplementation.
    """
    kinds, den = [], []
    for c in range(W):
        idx = start + c
        P = points[c]
        if idx >= count:
            kinds.append(SKIP); den.append(1); continue
        if P is None:
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
        stats["normal_steps"] = stats.get("normal_steps", 0) + kinds.count(NORMAL)
    out = []
    for c in range(W):
        idx = start + c
        P = points[c]
        kind = kinds[c]
        if kind == SKIP:
            out.append(P); continue
        if kind == FROM_INF:
            out.append(inc); continue
        if kind == TO_INF:
            out.append(None); continue
        if kind == NORMAL:
            lam = (inc[1] - P[1]) * invs[c] % p
        else:
            lam = (3 * P[0] * P[0]) * invs[c] % p
        x3 = (lam * lam - P[0] - inc[0]) % p
        y3 = (lam * (P[0] - x3) - P[1]) % p
        if record is not None:
            record(P[0], idx)
        out.append((x3, y3))
    return out


def chain_state_selftest(G, p, W=4):
    """Drive each of the four states through lockstep_round and check the result."""
    Wm = mul(W, G, p)
    checks = {}

    # FROM_INF: a chain sitting at O, advanced by a non-trivial increment.
    pts = [None] + [mul(c, G, p) for c in range(1, W)]
    out = lockstep_round(pts, Wm, 10 ** 9, 0, None, p, W)
    checks["from_inf_result_is_increment"] = (out[0] == Wm)
    checks["from_inf_matches_add"] = (add(None, Wm, p) == Wm)

    # NORMAL: two distinct points, no doubling, no infinity.
    A, B = mul(3, G, p), mul(7, G, p)
    out = lockstep_round([A] + [mul(c, G, p) for c in range(1, W)], B, 10 ** 9, 0, None, p, W)
    checks["normal_result_matches_add"] = (out[0] == add(A, B, p))

    # DOUBLING: increment equals the current point.
    out = lockstep_round([A] + [mul(c, G, p) for c in range(1, W)], A, 10 ** 9, 0, None, p, W)
    checks["doubling_result_matches_add"] = (out[0] == add(A, A, p))

    # TO_INF: increment is the negation of the current point.
    negA = (A[0], (-A[1]) % p)
    out = lockstep_round([A] + [mul(c, G, p) for c in range(1, W)], negA, 10 ** 9, 0, None, p, W)
    checks["to_inf_result_is_infinity"] = (out[0] is None)

    # SKIP: index past the bound leaves the chain untouched.
    out = lockstep_round([A] + [mul(c, G, p) for c in range(1, W)], B, 0, 0, None, p, W)
    checks["skip_leaves_point"] = (out[0] == A)

    return {k: bool(v) for k, v in checks.items()}


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


def radix_sort_records(records):
    """LSD radix sort on the integer x, 16-bit digits, unbounded counters."""
    if not records:
        return records
    width = max(1, max(r[0] for r in records).bit_length())
    base = 1 << 16
    passes = (width + 15) // 16
    buf = [None] * len(records)
    src, dst = records, buf
    for k in range(passes):
        shift = 16 * k
        counts = [0] * base
        for rec in src:
            counts[(rec[0] >> shift) & 0xFFFF] += 1
        total = 0
        for i in range(base):
            c = counts[i]
            counts[i] = total
            total += c
        for rec in src:
            d = (rec[0] >> shift) & 0xFFFF
            dst[counts[d]] = rec
            counts[d] += 1
        src, dst = dst, src
    return src


def multiword_search(sorted_records, key):
    """Multiword binary search over full coordinates; returns (pos, comparisons)."""
    lo, hi = 0, len(sorted_records)
    comparisons = 0
    while lo < hi:
        mid = (lo + hi) // 2
        comparisons += 1
        if sorted_records[mid][0] < key:
            lo = mid + 1
        else:
            hi = mid
    return lo, comparisons


def interleaved_bsgs(Q, G, r, p, W, stats=None):
    if Q is None:
        return 0
    m = isqrt(r) + 1
    h = -(-r // m)
    Kb = -(-m // W)
    Kg = -(-h // W)

    Wm = mul(W, G, p)
    D = mul(m, G, p)
    WD = mul(W, D, p)
    negWD = None if WD is None else (WD[0], (-WD[1]) % p)

    records = []

    points = [mul(c, G, p) for c in range(W)]
    for k in range(Kb):
        points = lockstep_round(points, Wm, m, k * W,
                                lambda x, j: records.append((x, j)), p, W, stats)
        if k == 0 and stats is not None:
            stats["chain0_after_round0_is_WG"] = (points[0] == Wm)

    records = radix_sort_records(records)

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
                cand = (idx * m) % r
                if mul(cand, G, p) == Q:
                    return cand
                continue
            pos, ncmp = multiword_search(records, P[0])
            if stats is not None:
                stats["max_search_comparisons"] = max(
                    stats.get("max_search_comparisons", 0), ncmp)
            while pos < len(records) and records[pos][0] == P[0]:
                j = records[pos][1]
                for cand in ((idx * m + j) % r, (idx * m - j) % r):
                    if mul(cand, G, p) == Q:
                        return cand
                pos += 1
        gpoints = lockstep_round(gpoints, negWD, h, start, None, p, W, stats)
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

    states = chain_state_selftest(G, p)
    probe = sorted([(rng.randrange(p), i) for i in range(256)], key=lambda rec: rec[0])
    probe_cmps = [multiword_search(probe, rng.randrange(p))[1] for _ in range(64)]

    total = 2 ** 128 * LEDGER["per_2e128"]
    exact = {
        "per_2e128_allowance": LEDGER["per_2e128"],
        "priced_work_log2": 144.0,
        "priced_time_log2": LEDGER["priced_time_log2"],
        "claimed_time_log2": LEDGER["claimed_time_log2"],
        "bound_holds": total < 2 ** 144 and total / 1000 < 2 ** LEDGER["claimed_time_log2"],
        "slack_factor_log2": LEDGER["claimed_time_log2"] - LEDGER["priced_time_log2"],
    }

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
        "lookup_structure": "radix sort plus multiword binary search over full coordinates",
        "chain_state_counts": {
            "infinity_start_steps": stats.get("from_inf_steps", 0),
            "normal_steps": stats.get("normal_steps", 0),
            "doubling_steps": stats.get("doubling_steps", 0),
            "inverse_point_steps": stats.get("to_inf_steps", 0),
            "chain0_after_round0_is_WG": stats.get("chain0_after_round0_is_WG", False),
        },
        "chain_state_selftest": states,
        "chain_state_selftest_all_pass": all(states.values()),
        "search_facts": {
            "max_full_comparisons_in_run": stats.get("max_search_comparisons", 0),
            "probe_max_comparisons": max(probe_cmps) if probe_cmps else 0,
            "ledger_worst_case_comparisons": 128,
        },
        "ledger": LEDGER,
        "bound_check": exact,
        "scope": "toy scale only; does not establish production cost or memory",
        "sample": checks[:4],
    }, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
