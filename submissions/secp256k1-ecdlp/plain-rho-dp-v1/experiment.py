"""Toy-scale check of the parallel r-adding Pollard rho with distinguished points.

Runs the same algorithm as proof.md (shared r-adding walk, distinguished-point
segments with fresh uniform restarts, segment cap, fixed global step budget,
sort-and-match of distinguished points, final verification) on small
prime-order curves y^2 = x^3 + 7 over F_p.  It measures the success rate at a
budget of c * sqrt(n / (1 - 2/R)) phase-1 steps.  It reports (a) the rate of a
repeated group element among phase-1 visited points against the heuristic-H
prediction 1 - exp(-k1^2 / (2 n / (1 - 2/R))), k1 = measured phase-1 points,
(b) the full algorithm's success rate, and (c) how often a phase-1 repeat was
detected and solved.  Every returned logarithm is verified.

Scope: toy curves only (n <= 2^32).  This is evidence for the random-walk
heuristic at small sizes, not a proof of production complexity, and it does
not measure instruction counts.  Standard library only.
"""
import json
import math
import random
import sys

# (p, n, G): y^2 = x^3 + 7 over F_p has prime order n; G is a finite point.
CURVES = {
    20: (1048783, 1050337, (59448, 508188)),
    22: (4194523, 4197601, (2899945, 1174969)),
    24: (16777723, 16773667, (14729108, 5846834)),
    26: (67109191, 67092871, (22863882, 31249551)),
    28: (268435987, 268414693, (37824420, 106902900)),
    30: (1073743417, 1073808937, (457478005, 379384870)),
    32: (4294968199, 4295016211, (1350594305, 260962323)),
}


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def add(P, Q, p):
    """Complete affine addition (None is the identity); used outside the walk."""
    if P is None:
        return Q
    if Q is None:
        return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % p == 0:
            return None
        lam = 3 * P[0] * P[0] * pow(2 * P[1], -1, p) % p
    else:
        lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
    x = (lam * lam - P[0] - Q[0]) % p
    return (x, (lam * (P[0] - x) - P[1]) % p)


def mul(k, P, p):
    R = None
    while k:
        if k & 1:
            R = add(R, P, p)
        P = add(P, P, p)
        k >>= 1
    return R


def check_curve(p, n, G):
    assert is_prime(p) and is_prime(n)
    assert (G[1] * G[1] - G[0] ** 3 - 7) % p == 0
    assert abs(p + 1 - n) <= 2 * math.isqrt(p) + 2
    assert n > 4 * math.isqrt(p) + 4          # so ord(G) = n forces #E = n
    assert mul(n, G, p) is None


def one_attack(p, n, G, Q, R, theta, walks, cap, budget, rng, stats):
    """One attack as in proof.md. Returns (log or None, phase1_points, phase1_repeat).

    Phase 1: ceil(budget / walks) steps per walk, round-robin; phase 2: cap more
    steps per walk under the same rules.  For measurement only, the toy also
    records every visited group element to tell whether a repeat (two visits of
    the same element) already happened among phase-1 visits: this is exactly
    the event that heuristic H bounds in proof.md."""
    table = []
    for _ in range(R):
        while True:                           # redraw M_j = O (as in proof.md)
            c, d = rng.randrange(n), rng.randrange(n)
            M = add(mul(c, G, p), mul(d, Q, p), p)
            if M is not None:
                table.append((M, c, d))
                break
    mask_part = R - 1
    part_bits = R.bit_length() - 1
    dp_mask = (1 << theta) - 1
    dps = []                                  # (x, start a, start b, length)
    seen = set()
    probe = {"phase1": True, "points1": 0, "repeat1": False}

    def visit(X):
        if X in seen:
            if probe["phase1"]:
                probe["repeat1"] = True
        else:
            seen.add(X)
        if probe["phase1"]:
            probe["points1"] += 1

    def fresh():
        while True:
            a, b = rng.randrange(n), rng.randrange(n)
            X = add(mul(a, G, p), mul(b, Q, p), p)
            if X is None:
                continue
            stats["starts"] += 1
            visit(X)
            if (X[0] >> part_bits) & dp_mask == 0:   # DP test also on start points
                dps.append((X[0], a, b, 0))
                continue
            return [X, a, b, 0]               # point, start a, start b, length

    state = [fresh() for _ in range(walks)]
    per_walk = -(-budget // walks)            # ceil, as in proof.md
    for t in range(per_walk + cap):
        if t == per_walk:
            probe["phase1"] = False
        for w in range(walks):
            X, a0, b0, ln = state[w]
            Rj = table[X[0] & mask_part][0]
            if X[0] == Rj[0]:
                stats["exceptional"] += 1     # incomplete addition: restart
                state[w] = fresh()
                continue
            lam = (Rj[1] - X[1]) * pow(Rj[0] - X[0], -1, p) % p
            x3 = (lam * lam - X[0] - Rj[0]) % p
            X = (x3, (lam * (X[0] - x3) - X[1]) % p)
            ln += 1
            stats["steps"] += 1
            visit(X)
            if (x3 >> part_bits) & dp_mask == 0:
                dps.append((x3, a0, b0, ln))
                state[w] = fresh()
            elif ln >= cap:
                stats["abandoned"] += 1
                state[w] = fresh()
            else:
                state[w] = [X, a0, b0, ln]
    dps.sort()
    pairs = 0
    for i in range(len(dps) - 1):
        if dps[i][0] != dps[i + 1][0]:
            continue
        pairs += 1
        if pairs > 16:                        # same 16-pair limit as proof.md
            break
        coeffs = []
        for (_, a0, b0, ln) in (dps[i], dps[i + 1]):
            X, a, b = add(mul(a0, G, p), mul(b0, Q, p), p), a0, b0
            for _ in range(ln):
                Rj, c, d = table[X[0] & mask_part]
                X = add(X, Rj, p)
                a, b = (a + c) % n, (b + d) % n
            coeffs.append((a, b, X))
        (a1, b1, X1), (a2, b2, X2) = coeffs
        if X1 != X2:                          # same x, opposite y: a1+b1x = -(a2+b2x)
            stats["negated_matches"] += 1
            a2, b2 = (-a2) % n, (-b2) % n
        if (b1 - b2) % n == 0:
            stats["degenerate"] += 1
            continue
        x = (a2 - a1) * pow(b1 - b2, -1, n) % n
        if mul(x, G, p) == Q:                 # charged final verification
            return x, probe["points1"], probe["repeat1"]
        stats["bad_verify"] += 1
    return None, probe["points1"], probe["repeat1"]


def main():
    cfg = json.load(sys.stdin)
    bits, trials, c_milli = cfg["bits"], cfg["trials"], cfg["c_milli"]
    R, theta, walks, seed = cfg["partitions"], cfg["theta"], cfg["walks"], cfg["seed"]
    cap_mult = cfg["cap_mult"]
    p, n, G = CURVES[bits]
    check_curve(p, n, G)
    c = c_milli / 1000.0
    n_eff = n / (1 - 2 / R)
    budget = int(c * math.sqrt(n_eff))
    cap = cap_mult << theta
    rng = random.Random(seed)
    stats = {"steps": 0, "starts": 0, "exceptional": 0, "abandoned": 0,
             "degenerate": 0, "bad_verify": 0, "negated_matches": 0}
    wins = repeats = solved_given_repeat = 0
    points1 = 0
    for _ in range(trials):
        x = rng.randrange(1, n)
        Q = mul(x, G, p)
        # input validation as in proof.md: canonical, on curve, [n]Q = O
        assert 0 <= Q[0] < p and 0 <= Q[1] < p
        assert (Q[1] * Q[1] - Q[0] ** 3 - 7) % p == 0 and mul(n, Q, p) is None
        got, k1, rep = one_attack(p, n, G, Q, R, theta, walks, cap, budget, rng, stats)
        points1 += k1
        repeats += rep
        if got is not None:
            assert got == x                   # every returned log is checked
            wins += 1
            solved_given_repeat += rep
    k1_mean = points1 / trials
    pred1 = 1 - math.exp(-k1_mean * k1_mean / (2 * n_eff))
    se1 = math.sqrt(pred1 * (1 - pred1) / trials)
    k_all = (stats["steps"] + stats["starts"]) / trials
    print(json.dumps({
        "curve": {"p": p, "n": n, "log2_n": round(math.log2(n), 3)},
        "partitions": R, "theta": theta, "walks": walks, "segment_cap": cap,
        "seed": seed, "trials": trials, "c": c, "phase1_budget_steps": budget,
        "mean_phase1_points": round(k1_mean, 1),
        "mean_all_points": round(k_all, 1),
        "restart_point_fraction": round(stats["starts"] / (stats["steps"] + stats["starts"]), 4),
        "H_test": {"phase1_repeat_rate": round(repeats / trials, 4),
                   "prediction_1_minus_exp": round(pred1, 4),
                   "z": round((repeats / trials - pred1) / se1, 2)},
        "algorithm_success_rate": round(wins / trials, 4),
        "solved_given_phase1_repeat": f"{solved_given_repeat}/{repeats}",
        "counters": stats,
        "scope": "toy curves only; tests heuristic H and DP detection at small n, "
                 "not production cost or instruction counts",
    }))

if __name__ == "__main__":
    main()
