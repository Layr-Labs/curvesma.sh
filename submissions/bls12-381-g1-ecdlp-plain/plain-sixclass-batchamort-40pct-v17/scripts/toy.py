#!/usr/bin/env python3
"""toy.py -- v17 fixed-fixture execution of the PRODUCTION spec, line by line.

The v17 spec (argument.md):
  state = (S, (a, c)) with S an affine G1 point and S = [a]G + [c]T;
  step  = S <- S + J_m, m = cell(|y(S)|), a <- a + r_m, c <- c + s_m.  There is
          NO per-step canonicalisation and NO per-step coefficient twist: the
          state is the point, the class key |y| only chooses the column and the
          distinguished point, and the orbit element e is resolved AT REPORT
          TIME from the two stored points (two Fp multiplies per hit).
  collision = two records with the same class key |y| and different points;
          with S1 = [e]S2 the report x = (e*a2 - a1)/(c1 - e*c2) is exact.

What this toy measures, and why each item is here:
  * invariant sweep -- S == [a]G + [c]T asserted after every transition, escape
    and restart, in a dedicated pass (count reported as `invariant_checks`);
  * `collision_report` -- for all six orbit elements e and for random pairs, the
    report path must recover e from the points and then x exactly.  This is the
    v17 replacement for v16's merge test; it is the invariant that the v16 crash
    violated (a coefficient pair that was twisted on one side only), and it is
    exact rather than statistical: 0 failures is required;
  * trail-length sweep at R = 2048 -- the H1 slowdown measurement: mean T_c over
    the random-mapping law against mean trail / sqrt(N), for dp_bits 0, 1, 2;
  * recovery rate -- rho at several kappa, on the pooled sweep.

Run: nice -n 10 python3 scripts/toy.py [--quick]
Writes evidence/toy.json next to the package root.
"""
import argparse
import hashlib
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import selector  # noqa: E402

P = 600217
B = 1
N_TOTAL = 600564
Q = 50047
GX, GY = 218652, 2045
BETA = 43364
LAMBDA = 26745
G = (GX, GY)
NCLS = (Q - 1) // 6

CHECKPOINT_PERIOD = 64
DP_BITS = 1
ORBIT_E = [1, Q - 1, LAMBDA, (Q - LAMBDA) % Q,
           LAMBDA * LAMBDA % Q, (Q - LAMBDA * LAMBDA) % Q]
_BETA_POW = [1, BETA, BETA * BETA % P]


def add(A, C):
    if A is None:
        return C
    if C is None:
        return A
    x1, y1 = A
    x2, y2 = C
    if x1 == x2 and (y1 + y2) % P == 0:
        return None
    if x1 == x2 and y1 == y2:
        if y1 == 0:
            return None
        lam = 3 * x1 * x1 % P * pow(2 * y1, -1, P) % P
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, P) % P
    x3 = (lam * lam - x1 - x2) % P
    return (x3, (lam * (x1 - x3) - y1) % P)


def mul(k, A):
    R, S = None, A
    k %= Q
    while k:
        if k & 1:
            R = add(R, S)
        S = add(S, S)
        k >>= 1
    return R


def class_key(S):
    """|y|: invariant under all six orbit elements -- the class key."""
    y = S[1]
    return y if y <= P - y else P - y


def orbit_apply(e, S):
    """[e]S for e in {+-1, +-lam, +-lam^2}, as the six explicit maps.

    The (k, sign) decomposition is read off ORBIT_E itself, the same table
    resolve_e returns indices into, so the two functions are exact inverses by
    construction.  Deriving the sign from `e >= Q // 2` is WRONG: lam = 26745
    exceeds Q // 2 = 25023 while lam is a positive representative, so that test
    applied the lam^2 map to lam and to -lam and broke e-resolution on exactly
    those two of the six elements.
    """
    if S is None:
        return None
    if e not in ORBIT_E:
        raise ValueError("orbit element not in the six-class orbit")
    k, neg = divmod(ORBIT_E.index(e), 2)
    x, y = S
    if k:
        x = x * _BETA_POW[k] % P
    if neg:
        y = (P - y) % P
    return (x, y)


def resolve_e(S1, S2):
    """The report-time orbit element: the e in S1 = [e]S2, or None if not mates.

    Two Fp multiplies (beta*x2, beta^2*x2) plus three compares -- charged once
    per record hit, not per step.
    """
    if S1 is None or S2 is None:
        return None
    x1, y1 = S1
    x2, y2 = S2
    bx = BETA * x2 % P
    if x1 == x2:
        k = 0
    elif x1 == bx:
        k = 1
    elif x1 == (BETA * BETA % P) * x2 % P:
        k = 2
    else:
        return None
    same = (y1 & 1) == (y2 & 1)
    return ORBIT_E[2 * k + (0 if same else 1)]


def recover_x(rec, S, a, c):
    """rec = (key, a0, c0, S0); current state (S, a, c) is its class mate."""
    _k, a0, c0, S0 = rec
    e = resolve_e(S0, S)
    if e is None:
        return None, None
    den = (c0 - e * c) % Q
    if den == 0:
        return e, None
    return e, (e * a - a0) * pow(den, -1, Q) % Q


def is_dp(K):
    return (K * 2654435761 % (1 << 20)) < (1 << (20 - DP_BITS))


def run(pool, tab, T, budget, A0_pair, stream, R, check_invariant=False,
        solve=True):
    a, c = A0_pair
    S = add(mul(a, G), mul(c, T))
    if S is None:
        a, c = stream.next_pair()
        S = add(mul(a, G), mul(c, T))
    st = {"steps": 0, "escapes": 0, "restarts": 0, "dps": 0, "hits": 0,
          "merges": 0, "false_keys": 0, "inv": 0, "checks": 1}
    kcp, cp = None, None
    while st["steps"] < budget:
        K = class_key(S)
        m = selector.cell(K, R, P)
        st["checks"] += 1
        J, rj, sj = tab[m]
        if is_dp(K):
            st["dps"] += 1
            prev = pool.get(K)
            if prev is not None:
                if prev[3] == S:
                    # MERGE: the two walks coincide in the point, so the pairs
                    # differ by a multiple of (x, -1).  Solving this is the same
                    # report as a class collision and must NOT be dropped -- the
                    # v16 crash was exactly a merge path that lost the pair.
                    st["merges"] += 1
                    den = (c - prev[2]) % Q
                    if den:
                        xm = (prev[1] - a) * pow(den, -1, Q) % Q
                        if solve and mul(xm, G) == T:
                            return st["steps"], st, xm
                else:
                    e, x = recover_x(prev, S, a, c)
                    st["hits"] += 1
                    if e is None:
                        st["false_keys"] += 1      # same |y| but a different class
                    elif solve and x is not None and mul(x, G) == T:
                        return st["steps"], st, x
            else:
                pool[K] = (K, a, c, S)
        if kcp is not None and K == kcp and cp == S:
            cp2 = add(cp, cp)                       # doubling escape: point cycle
            if cp2 is None:
                a, c = stream.next_pair()
                S = add(mul(a, G), mul(c, T))
                st["restarts"] += 1
            else:
                S, a, c = cp2, 2 * a % Q, 2 * c % Q
                st["escapes"] += 1
            kcp, cp = None, None
            st["steps"] += 1
            if check_invariant:
                assert add(mul(a, G), mul(c, T)) == S, "invariant broken at escape"
                st["inv"] += 1
            continue
        if st["steps"] % CHECKPOINT_PERIOD == CHECKPOINT_PERIOD - 1:
            kcp, cp = K, S
        S = add(S, J)
        if S is None:
            a, c = stream.next_pair()               # FRESH seeded draw
            S = add(mul(a, G), mul(c, T))
            st["restarts"] += 1
        else:
            a, c = (a + rj) % Q, (c + sj) % Q
        st["steps"] += 1
        if check_invariant:
            assert add(mul(a, G), mul(c, T)) == S, "invariant broken"
            st["inv"] += 1
    return None, st, None


def build_table(T, R, seed, n=64):
    h = hashlib.sha256(b"v17pool:%d:%d" % (seed, Q)).digest()
    gs, ts = [], []
    for _ in range(n):
        h = hashlib.sha256(h).digest()
        gs.append((mul(int.from_bytes(h[:4], "big") % Q, G),
                   int.from_bytes(h[:4], "big") % Q))
        ts.append((mul(int.from_bytes(h[4:8], "big") % Q, T),
                   int.from_bytes(h[4:8], "big") % Q))
    tab = []
    h2 = hashlib.sha256(b"v17:%d:%d" % (seed, Q)).digest()
    for _m in range(R):
        h2 = hashlib.sha256(h2).digest()
        gi = int.from_bytes(h2[:2], "big") % n
        ti = int.from_bytes(h2[2:4], "big") % n
        tab.append((add(gs[gi][0], ts[ti][0]), gs[gi][1], ts[ti][1]))
    return tab


def collision_report(tab, T, x, trials, seed):
    """The report invariant, exact: for every e and random pairs, recover e and x.

    Build two states that are class mates by construction: S2 = [a2]G + [c2]T,
    S1 = [e]S2 with a1 = e*a2, c1 = e*c2.  Then the resolver must return e and
    x = (e*a2 - a1)/(c1 - e*c2).  A coefficient pair twisted on one side only --
    the v16 toy crash -- fails this immediately.
    """
    out = {"trials": 0, "e_resolved": 0, "x_recovered": 0, "failures": [],
           "failed": 0, "per_e": {}}
    total_failed = 0
    for ei, e in enumerate(ORBIT_E):
        ok_e = ok_x = seen = 0
        for t in range(trials):
            # two INDEPENDENT pairs whose points are class mates: pick a1, a2, c2
            # at random and solve c1 from a1 + c1*x = e*(a2 + c2*x).  (Setting
            # a1,c1 = e*a2,e*c2 instead would make the denominator vanish -- that
            # degenerate construction is the v16 crash pattern.)
            a2 = (t * 92821 + 17 + ei) % Q
            c2 = (t * 68917 + 23 + 7 * ei) % Q
            a1 = (t * 15127 + 3 * ei + 91) % Q
            inv_x = pow(x, -1, Q)
            c1 = ((e * a2 - a1) * inv_x + e * c2) % Q
            S2 = add(mul(a2, G), mul(c2, T))
            if S2 is None:
                continue
            S1 = orbit_apply(e, S2)
            seen += 1
            if add(mul(a1, G), mul(c1, T)) != S1:
                total_failed += 1
                out["failures"].append({"e": e, "trial": t, "stage": "state"})
                continue
            # stored record is S1's, the live state is S2's (and vice versa)
            e_fwd = resolve_e(S1, S2)
            e_bwd = resolve_e(S2, S1)
            x_fwd = recover_x((0, a1, c1, S1), S2, a2, c2)[1]
            x_bwd = recover_x((0, a2, c2, S2), S1, a1, c1)[1]
            out["trials"] += 1
            if e_fwd == e and e_bwd == pow(e, -1, Q):
                ok_e += 1
            else:
                total_failed += 1
                out["failures"].append({"e": e, "trial": t, "stage": "resolve",
                                        "got": [e_fwd, e_bwd]})
            if x_fwd == x and x_bwd == x:
                ok_x += 1
            else:
                total_failed += 1
                out["failures"].append({"e": e, "trial": t, "stage": "x",
                                        "got": [x_fwd, x_bwd]})
        out["e_resolved"] += ok_e
        out["x_recovered"] += ok_x
        out["per_e"][str(e)] = {"trials": seen, "e_ok": ok_e, "x_ok": ok_x}
    out["failed"] = total_failed           # counted BEFORE the sample is truncated
    out["failures"] = out["failures"][:4]
    return out


def s_ref(R):
    Np = NCLS * R / (R - 1.0)
    return (2.0 * Np * (-math.log(0.6097))) ** 0.5


def sweep(R, dplist, trials, seed0, check_invariant=False):
    """Trail-length sweep at fixed R: mean T_c / S_ref and rho at several kappa."""
    global DP_BITS
    Sref = s_ref(R)
    kappas = [1.0, 1.005, 1.02, 1.025, 1.03]
    rows = []
    for dpb in dplist:
        DP_BITS = dpb
        hits, checks = [], 0
        for t in range(trials):
            x = ((seed0 + t) * 7919 + 13) % (Q - 1) + 1
            T = mul(x, G)
            tab = build_table(T, R, seed0 + t)
            pool = {}
            stream = selector.RestartStream(b"toy17-%d-%d" % (seed0 + t, R), Q)
            budget = int(round(3.0 * Sref)) + (2 << dpb)
            step, st, _ = run(pool, tab, T, budget,
                              ((t * 104729 + 7 + dpb) % Q, (t * 15485863 + 11) % Q),
                              stream, R, check_invariant)
            checks += st["checks"]
            hits.append(step if step is not None else budget)
        n = len(hits)
        mean_tc = sum(hits) / n / Sref
        var = sum((h / Sref - mean_tc) ** 2 for h in hits) / max(1, n - 1)
        row = {"R": R, "dp_bits": dpb, "trials": n,
               "mean_trail_over_sqrtN": round((2 ** dpb) / math.sqrt(NCLS), 6),
               "mean_tc_over_Sref": round(mean_tc, 5),
               "se": round(math.sqrt(var / n), 5),
               "slowdown_pct": round((mean_tc - 1) * 100, 3),
               "recovery_rate": {("kappa_%.3f" % k): round(
                   sum(1 for h in hits if h <= k * Sref) / n, 5) for k in kappas},
               "selector_checks": checks}
        rows.append(row)
    return {"S_ref": round(Sref, 5), "N_classes": NCLS, "rows": rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    global DP_BITS, CHECKPOINT_PERIOD
    x = 12345 % (Q - 1) + 1
    T = mul(x, G)
    R = 2048
    out = {"fixture": {"p": P, "q": Q, "classes": NCLS, "beta": BETA,
                       "lambda": LAMBDA, "order_total": N_TOTAL},
           "selector": "scripts/selector.py",
           "reachable_columns": {"R": R, "reachable": selector.reachable_columns(R, P)},
           "checkpoint_period": CHECKPOINT_PERIOD}

    # 1. invariant pass: every transition asserted
    tab = build_table(T, R, 11)
    pool = {}
    stream = selector.RestartStream(b"toy17-inv", Q)
    DP_BITS = 0
    # solve=False: the invariant pass must execute the WHOLE budget, so the
    # assertion count equals the transition count instead of stopping at the
    # first collision (the earlier form could never exceed ~100 checks).
    step, st, _ = run(pool, tab, T, 4000, (1234567 % Q, 7654321 % Q), stream, R,
                      check_invariant=True, solve=False)
    out["invariant_checks"] = st["inv"]
    out["invariant_pass_steps"] = st["steps"]
    out["selector_checks"] = st["checks"]

    # 2. report invariant (all six e), exact
    DP_BITS = 1
    cr = collision_report(tab, T, x, 25 if args.quick else 60, 5)
    out["collision_report"] = cr

    # 3. trail-length sweep
    dplist = [0] if args.quick else [0, 1, 2]
    trials = 120 if args.quick else 1400
    out["sweep"] = sweep(R, dplist, trials, 6000)
    out["A_measured"] = None
    best = None
    # smallest printed kappa whose pooled dp0 recovery is >= 0.39
    dp0 = [r for r in out["sweep"]["rows"] if r["dp_bits"] == 0]
    if dp0:
        rec = dp0[0]["recovery_rate"]
        for k in sorted(rec, key=lambda s: float(s.split("_")[1])):
            if rec[k] >= 0.39:
                best = float(k.split("_")[1])
                break
        out["A_measured"] = best
    out["verdict"] = ("PASS" if cr["failed"] == 0 and cr["x_recovered"] == cr["trials"]
                      and st["inv"] == st["steps"] and st["inv"] > 1000 else "FAIL")
    out["pass_rule"] = ("PASS iff (i) the report invariant recovers e and x on every "
                        "trial for all six e (failed == 0, x_recovered == trials), "
                        "(ii) invariant_checks == invariant_pass_steps > 1000, i.e. the "
                        "invariant was asserted after every single transition, escape "
                        "and restart of the invariant pass, and (iii) the sweep prints "
                        "the slowdown law and the recovery rate at the charged budget "
                        "kappa = A_H1 = 1.005")
    pkg = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs(os.path.join(pkg, "evidence"), exist_ok=True)
    with open(os.path.join(pkg, "evidence", "toy.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(json.dumps({"verdict": out["verdict"], "collision_failed": cr["failed"],
                      "collision_x": [cr["x_recovered"], cr["trials"]],
                      "invariant_checks": out["invariant_checks"],
                      "sweep": [(r["dp_bits"], r["slowdown_pct"], r["se"],
                                 r["recovery_rate"]["kappa_1.020"]) for r in out["sweep"]["rows"]],
                      "A_measured": out["A_measured"]}))
    return 0 if out["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
