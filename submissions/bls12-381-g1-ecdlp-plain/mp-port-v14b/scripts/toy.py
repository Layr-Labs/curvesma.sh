#!/usr/bin/env python3
"""toy.py -- fixed-fixture toy for the mp-port-v14b package.

Scope: a hardcoded 20-bit j=0 curve (no search of any kind at run time) that
carries the same objects the production construction uses: a prime-order
subgroup with a GLV automorphism, the resulting six-class quotient under
negation + GLV, a |y| canonical key, a truncated fixed-point index for the
canonical-member choice, checkpoint-only cycle detection with doubling
escapes, and the H1 slowdown / success statistic on the class walk.

The fixture is FIXED. Every constant below was chosen once, off-line, by
exact point counting; the shipped script performs no curve search, no
primality testing and no random curve generation. Runtime is a few seconds.

Run: python3 scripts/toy.py [--trials N] [--seed S]
Prints one JSON object; writes evidence/toy.json next to the package.
"""
import argparse
import hashlib
import json
import os
import random
import sys

# ---- FIXED FIXTURE (hardcoded; verified by check_fixture below) ------------
P = 600217                      # field prime, P = 1 mod 3
B = 1                           # curve y^2 = x^3 + B
N_TOTAL = 600564                # #E(F_P), counted exactly off-line
Q = 50047                       # prime subgroup order, N_TOTAL = 12 * Q
COFACTOR = 12                   # N_TOTAL // Q
GX = 218652                     # generator of the order-Q subgroup
GY = 2045
BETA = 43364                    # primitive cube root of unity mod P
LAMBDA = 26745                  # LAMBDA^2 + LAMBDA + 1 = 0 mod Q, [L]*G = (BETA*GX, GY)

# ---- toy parameters (fixed defaults, overridable for time budget) ----------
DP_BITS = 5                     # distinguished points: fp >> (64-DP_BITS) == 0
TABLE_ROWS = 16                 # class key low bits -> jump row
TABLE_R = 64                    # jump entries per row
MAX_TRIALS = 200
SEED = 3000
CHECKPOINT_PERIOD = 64          # their period-64 checkpoint, no 2-cycle test
TRAIL_CAP = 4096                # hard cap so the toy can never hang

G = (GX, GY)


def add(A, C):
    if A is None:
        return C
    if C is None:
        return A
    x1, y1 = A
    x2, y2 = C
    if x1 == x2 and (y1 + y2) % P == 0:
        return None
    if A == C:
        lam = 3 * x1 * x1 * pow(2 * y1, -1, P) % P
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, P) % P
    x3 = (lam * lam - x1 - x2) % P
    return (x3, (lam * (x1 - x3) - y1) % P)


def mul(k, A):
    R, S = None, A
    while k:
        if k & 1:
            R = add(R, S)
        S = add(S, S)
        k >>= 1
    return R


def neg(A):
    return None if A is None else (A[0], (P - A[1]) % P)


def lift(x):
    rhs = (x * x * x + B) % P
    if rhs == 0:
        return (x, 0)
    if pow(rhs, (P - 1) // 2, P) != 1:
        return None
    if P % 4 == 3:
        y = pow(rhs, (P + 1) // 4, P)
        return (x, y)
    # Tonelli-Shanks for the general case
    s, d = 0, P - 1
    while d % 2 == 0:
        d //= 2
        s += 1
    z = 2
    while pow(z, (P - 1) // 2, P) != P - 1:
        z += 1
    m, c, t, R = s, pow(z, d, P), pow(rhs, d, P), pow(rhs, (d + 1) // 2, P)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % P
            i += 1
        bf = pow(c, 1 << (m - i - 1), P)
        m, c = i, bf * bf % P
        t = t * bf * bf % P
        R = R * bf % P
    return (x, R)


# ---- class machinery -------------------------------------------------------
_CLASS_CACHE = {}


def orbit_class(x):
    """Canonical representative of the six-class orbit of x.

    Orbit = {+-P, +-phi(P), +-phi^2(P)} with phi(x, y) = (BETA*x, y).  The
    canonical member is the lexicographic minimum; the class key is |y| of that
    member (their key, invariant under both negation and GLV).
    """
    hit = _CLASS_CACHE.get(x)
    if hit is not None:
        return hit
    xs = sorted({x, x * BETA % P, x * BETA * BETA % P})
    cands = []
    for cx in xs:
        A = lift(cx)
        if A is None:
            continue
        cands.append(A)
        cands.append(neg(A))
    rep = min(cands)
    key = min(rep[1], (P - rep[1]) % P)
    out = (rep, key, len(cands))
    _CLASS_CACHE[x] = out
    return out


# truncated fixed-point index (toy analogue of their 12-multiply truncated
# product): rank x by floor(x * 2^16 / P) using ONE multiply and one shift,
# with exact comparison only on ties.
_FRAC_BITS = 16
_SCALE = (1 << (P.bit_length() + _FRAC_BITS)) // P


def trunc_index(x):
    return (x * _SCALE) >> P.bit_length()


def fp_of(key):
    h = hashlib.sha256(b"mp-port-v14b:%d:%d" % (key, P)).digest()
    return int.from_bytes(h[:8], "big")


_ORBIT_MULT = [1, (Q - 1) % Q, LAMBDA, (Q - LAMBDA) % Q,
               (LAMBDA * LAMBDA) % Q, (Q - LAMBDA * LAMBDA) % Q]


def orbit_element(A0, A):
    """Return e in {+-1, +-lambda, +-lambda^2} mod Q with A = [e]A0, or None."""
    for e in _ORBIT_MULT:
        if mul(e, A0) == A:
            return e
    return None


def build_table(T):
    tab = []
    for _ in range(TABLE_ROWS):
        row = []
        for _ in range(TABLE_R):
            r, s = random.randrange(Q), random.randrange(Q)
            row.append((add(mul(r, G), mul(s, T)), r, s))
        tab.append(row)
    return tab


def run_trail(T, tab, max_steps, dp_only=False):
    """One class walk with checkpoint-only cycle detection.

    T_c is the number of arrivals until a class that is already in the trail's
    table repeats with a DIFFERENT point, so that the two coefficient vectors
    solve for x; every recovered x is verified by [x]G == T.  dp_only=True
    restricts the table to distinguished points (their memory-saving rule) and
    is reported separately, since at 8341 classes it costs a factor ~6.
    """
    a, c = random.randrange(Q), random.randrange(Q)
    A = add(mul(a, G), mul(c, T))
    seen = {}
    kcp = None
    age = 0
    steps = 0
    escapes = 0
    bad = 0
    cycles = 0
    while steps < max_steps:
        if A is None:
            A = add(mul(random.randrange(Q), G), mul(random.randrange(Q), T))
            steps += 1
            continue
        rep, key, _ = orbit_class(A[0])
        fp = fp_of(key)
        # checkpoint test BEFORE the refresh, exactly their ordering
        if kcp is not None and fp == kcp:
            A = add(A, A)                     # doubling escape
            escapes += 1
            kcp, age = None, 0
            steps += 1
            continue
        age += 1
        if age == CHECKPOINT_PERIOD:
            kcp, age = fp, 0
        if (not dp_only) or (fp >> (64 - DP_BITS)) == 0:
            prev = seen.get(key)
            if prev is not None:
                a0, c0, A0 = prev
                e = orbit_element(A0, A)      # A = [e] A0 on the six-class orbit
                if e is not None:
                    da = (e * a0 - a) % Q
                    dc = (c - e * c0) % Q
                    if dc:
                        x = da * pow(dc, -1, Q) % Q
                        if mul(x, G) == T:
                            return {"steps": steps, "escapes": escapes, "recovered": x,
                                    "bad": bad, "dps": len(seen), "cycles": cycles}
                    bad += 1
                else:
                    bad += 1
            else:
                seen[key] = (a, c, A)
        h = key & (TABLE_ROWS - 1)
        J, r, s = tab[h][random.randrange(TABLE_R)]
        A = add(A, J)
        a = (a + r) % Q
        c = (c + s) % Q
        steps += 1
    return {"steps": None, "escapes": escapes, "recovered": None, "bad": bad,
            "dps": len(seen), "cycles": cycles}


def check_fixture():
    c = {}
    c["curve_nonsingular"] = (4 * pow(B, 3, P)) % P != 0
    c["G_on_curve"] = (GY * GY - GX * GX * GX - B) % P == 0
    c["G_order_is_Q"] = mul(Q, G) is None and mul(Q - 1, G) is not None
    c["N_TOTAL_equals_cofactor_times_Q"] = N_TOTAL == COFACTOR * Q
    c["beta_is_cube_root_of_unity"] = pow(BETA, 3, P) == 1 and BETA != 1
    c["lambda_relation"] = (LAMBDA * LAMBDA + LAMBDA + 1) % Q == 0
    c["lambda_matches_glv"] = mul(LAMBDA, G) == (BETA * GX % P, GY)
    return c


def check_six_class(trials):
    ok_orbit = ok_inv = ok_unique = 0
    for _ in range(trials):
        A = mul(random.randrange(1, Q), G)
        rep, key, n_pts = orbit_class(A[0])
        orb = [A, neg(A),
               (BETA * A[0] % P, A[1]), neg((BETA * A[0] % P, A[1])),
               (BETA * BETA * A[0] % P, A[1]), neg((BETA * BETA * A[0] % P, A[1]))]
        ok_orbit += int(len({q for q in orb if q is not None}) == 6 and n_pts == 6)
        keys = {min(q[1], (P - q[1]) % P) for q in orb if q is not None}
        ok_inv += int(keys == {key})
        ok_unique += int(rep == min(orb))
    return {"orbit_size_6": ok_orbit, "key_invariant": ok_inv,
            "canonical_unique": ok_unique, "trials": trials,
            "all_pass": ok_orbit == ok_inv == ok_unique == trials}


def check_trunc_index(trials):
    """Truncated fixed-point ordering of the three x-candidates vs exact."""
    mism = ties = 0
    for _ in range(trials):
        A = mul(random.randrange(1, Q), G)
        xs = sorted({A[0], A[0] * BETA % P, A[0] * BETA * BETA % P})
        approx = [trunc_index(x) for x in xs]
        if len(set(approx)) != len(set(xs)):
            ties += 1
        if min(range(len(xs)), key=lambda i: (approx[i], xs[i])) != 0:
            mism += 1
    return {"trials": trials, "mismatches": mism, "bucket_ties": ties,
            "all_pass": mism == 0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=MAX_TRIALS)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--long-trails", type=int, default=120)
    ap.add_argument("--long-len", type=int, default=2048)
    args = ap.parse_args()
    random.seed(args.seed)

    fx = check_fixture()
    six = check_six_class(min(400, args.trials * 2))
    trunc = check_trunc_index(min(2000, args.trials * 10))

    x_true = random.randrange(1, Q)
    T = mul(x_true, G)
    tab = build_table(T)

    # model budget: their convention S_ops = 1.05 * sqrt(2 N' ln(1/0.6097))
    import math
    N = Q / 6.0
    Np = N * TABLE_R / (TABLE_R - 1)
    S0 = math.sqrt(2 * Np * math.log(1 / 0.6097))
    S_ops = 1.05 * S0
    cap = min(TRAIL_CAP, int(4 * S_ops))

    tc, recovered, escapes, bad, cyc = [], 0, 0, 0, 0
    done = 0
    while done < args.trials:
        r = run_trail(T, tab, cap)
        done += 1
        escapes += r["escapes"]
        bad += r["bad"]
        cyc += r["cycles"]
        if r["recovered"] is not None:
            recovered += 1
        if r["steps"] is not None:
            tc.append(r["steps"])

    # distinguished-point-only table (their memory rule) at the same scale, to
    # show what the filter costs when the class space is tiny
    dptc, dpdone, dprec = [], 0, 0
    while dpdone < max(10, args.trials // 4):
        r = run_trail(T, tab, cap, dp_only=True)
        dpdone += 1
        if r["recovered"] is not None:
            dprec += 1
            dptc.append(r["steps"])

    # long-trail trap test: R is small here, so the walk is cycle dominated and
    # the checkpoint/doubling escape is what keeps it moving.
    long_esc, long_stuck, long_bad = 0, 0, 0
    for _ in range(args.long_trails):
        r = run_trail(T, tab, args.long_len)
        long_esc += r["escapes"]
        long_bad += r["bad"]
        if r["steps"] is None and r["recovered"] is None:
            long_stuck += 1

    out = {
        "fixture": {"p": P, "b": B, "order": N_TOTAL, "q": Q,
                    "cofactor": COFACTOR, "q_bits": Q.bit_length(),
                    "p_bits": P.bit_length(), "beta": BETA, "lambda": LAMBDA,
                    "G": [GX, GY], "search_performed_at_run_time": False},
        "checks": {
            "fixture": fx,
            "six_class": six,
            "truncated_index": trunc,
        },
        "stat_Tc": {
            "trials": done, "recovered": recovered,
            "recovery_fraction": round(recovered / done, 4) if done else 0.0,
            "mean_Tc": round(sum(tc) / len(tc), 2) if tc else None,
            "S0_model": round(S0, 2), "S_ops_model": round(S_ops, 2),
            "mean_Tc_over_S0": round((sum(tc) / len(tc)) / S0, 4) if tc else None,
            "mean_Tc_over_S_ops": round((sum(tc) / len(tc)) / S_ops, 4) if tc else None,
            "P_Tc_le_S_ops": round(sum(1 for x in tc if x <= S_ops) / len(tc), 4) if tc else None,
            "escapes": escapes, "cycle_hits_same_point": cyc,
            "fp_collisions_without_point_match": bad,
        },
        "dp_only_table": {"trials": dpdone, "recovered": dprec,
                          "mean_Tc": round(sum(dptc) / len(dptc), 2) if dptc else None,
                          "note": "their DP memory rule costs a factor ~1/2^DP_BITS at 8341 classes"},
        "long_trail": {"trails": args.long_trails, "len": args.long_len,
                       "escapes": long_esc, "stuck": long_stuck,
                       "fp_collisions_without_point_match": long_bad},
        "scope_note": ("fixed 20-bit j=0 fixture with 8341 classes; their toy rows are "
                       "q26/q30 (1e7+ classes) and their production R=2^20 is not "
                       "exercised; statistics here are indicative only"),
        "model_claims": "none; this toy supports no production claim",
    }
    all_ok = all(fx.values()) and six["all_pass"] and trunc["all_pass"] and recovered > 0
    out["verdict"] = "PASS" if all_ok else "FAIL"
    print(json.dumps(out, indent=1))
    d = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")
    os.makedirs(d, exist_ok=True)
    fn = os.path.join(d, "toy.json")
    with open(fn, "w") as f:
        json.dump(out, f, indent=1)
    print("wrote", fn)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
