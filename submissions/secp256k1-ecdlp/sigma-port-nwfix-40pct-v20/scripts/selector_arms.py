#!/usr/bin/env python3
"""selector_arms.py -- v20: selector dimension, MEASURED.

TWO ROWS, both requested by review:

(A) R/2 REACH DEFECT (shared with p256 v11).  qwalk.col feeds selector.cell_vec
    K uniform on [0, 2^39) (qwalk.py:186-187) and shifts by bitlen(2^40)-1 = 40
    (qwalk.py:193), so the largest reachable bucket is (2^39-1)*R >> 40 = R/2 - 1:
    only HALF the R columns are reachable from a uniform 39-bit key.  The package
    already books this as R_eff = reachable_columns(R, 2^40) (qwalk.py:145) and
    charges the trap rate 1/(6*R_eff), i.e. it priced the DEFECT rather than
    fixing it.  This script (i) measures the coverage of both variants, (ii)
    re-runs the floor-bearing cell at W=16 with a FULL-REACH selector (40-bit key,
    same map, mask fixed) and reports the diff.

(B) EC-COORDINATE SELECTOR vs HASH SURROGATE (verdict-v18 point 1: "hash-derived
    surrogate column selector" vs "elliptic-curve coordinates").  The surrogate
    selects on sm64(class).  The production selector selects on the y-coordinate
    of a real curve point, K = min(y, p-y).  This script builds a REAL curve
    y^2 = x^3 + 7 over a 24-bit prime p = 1 mod 3 (j=0, mu_6 automorphism, the
    secp256k1 shape), walks REAL points P -> P + u[m]G with the SAME table
    construction, the SAME cell(K, R, p) map with full reach, the SAME two-step
    return repair, the SAME fingerprint DP gate and the SAME same-challenge two-run
    pairing, and reports p_set against the surrogate at matched x, R, W, d, sets.

Output: evidence/selector-swap.json
"""
import hashlib
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk                  # noqa: E402
import transfer_arm as T      # noqa: E402
import census_v20n as C       # noqa: E402
import selector               # noqa: E402

OUT = os.path.join(D, "evidence", "selector-swap.json")
M40 = (1 << 40) - 1


def col_full(self, canon):
    """qwalk.Universe.col with the mask FIXED: a 40-bit key reaches all R columns."""
    K = qwalk.sm64(canon ^ np.uint64(0x243F6A8885A308D3)) >> np.uint64(24)
    K = K.astype(np.int64)                       # uniform on [0, 2^40)
    assert int(K.max(initial=0)) < (1 << 40)
    assert self.R.bit_length() + 40 < 63, (self.R, "K*R would overflow int64")
    return selector.cell_vec(K, self.R, 1 << 40).astype(np.int64)


def coverage(draws=400000, R=64, seed=3):
    rng = np.random.default_rng(seed)
    out = {}
    for name, bits in (("v18_39bit", 39), ("v20_40bit", 40)):
        K = rng.integers(0, 1 << bits, size=draws, dtype=np.int64)
        m = selector.cell_vec(K, R, 1 << 40)
        u = np.unique(m)
        out[name] = {"distinct_columns": int(u.size), "R": R,
                     "coverage_fraction": float(u.size / R),
                     "R_eff_measured": int(u.size),
                     "min": int(u.min()), "max": int(u.max())}
    out["diff"] = {"columns_gained": out["v20_40bit"]["distinct_columns"]
                   - out["v18_39bit"]["distinct_columns"],
                   "trap_rate_ratio": (out["v18_39bit"]["R_eff_measured"]
                                       / out["v20_40bit"]["R_eff_measured"]),
                   "bits_of_trap_row": math.log2(out["v18_39bit"]["R_eff_measured"]
                                                / out["v20_40bit"]["R_eff_measured"])}
    return out


# ---- (B) real EC walk -------------------------------------------------------
def prime_1mod3(minbits):
    n = (1 << minbits) + 1
    while True:
        if n % 3 == 1 and qwalk.is_prime(n):
            return n
        n += 1


class EC:
    """y^2 = x^3 + 7 over F_p, p = 1 mod 3, j = 0 (mu_6 automorphism)."""

    def __init__(self, p, R, seed=20250):
        self.p = p
        self.R = R
        self.G = self._find_point()
        u = np.empty(R, dtype=np.int64)
        for m in range(R):
            dd = hashlib.sha256(b"v20b-col:%d:%d:%d" % (seed, p, m)).digest()
            u[m] = (int.from_bytes(dd[:8], "big") % (p - 1)) + 1
        self.u = u.tolist()
        self.T = [self.mul(self.G, k) for k in self.u]        # increment table
        # class count: K = min(y, p-y) in [0, p/2]
        self.classes = (p - 1) // 2

    def _find_point(self):
        p = self.p
        for x in range(1, 10000):
            rhs = (x * x * x + 7) % p
            if pow(rhs, (p - 1) // 2, p) != 1:
                continue
            if p % 4 == 3:
                y = pow(rhs, (p + 1) // 4, p)
            else:
                y = None
                for _ in range(64):
                    pass
                continue                      # only p = 3 mod 4 curves are used
            if (y * y - rhs) % p == 0:
                return (x, y)
        raise RuntimeError("no point")

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        p = self.p
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1) * pow(2 * y1, p - 2, p) % p
        else:
            lam = (y2 - y1) * pow((x2 - x1) % p, p - 2, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)

    def mul(self, P, k):
        R, Q = None, P
        while k:
            if k & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            k >>= 1
        return R

    def key(self, P):
        y = P[1]
        return y if y <= self.p - y else self.p - y

    def dp(self, K, d):
        if d <= 0:
            return True
        return int(qwalk.sm64(np.array([K], dtype=np.uint64))[0]) < ((1 << 64) >> d)

    def cell(self, K):
        return int(selector.cell_vec(np.array([K], dtype=np.int64), self.R, self.p)[0])


def ec_run(E, sets, W, per, d, seed):
    """One run: `sets` challenges x W walkers x `per` steps, real EC points."""
    rng = np.random.default_rng(seed)
    solved = np.zeros(sets, dtype=bool)
    fires = restarts = 0
    for s in range(sets):
        rec = {}
        common = E.mul(E.G, int(rng.integers(1, E.p)))         # shared start (challenge)
        for w in range(W):
            P = E.add(common, E.mul(E.G, int(rng.integers(1, E.p))))
            prevK = None
            run = 0
            for _t in range(per):
                K = E.key(P)
                m = E.cell(K)
                if prevK is not None and K == prevK:            # two-step return repair
                    m = (m + 1) % E.R
                    fires += 1
                    run += 1
                    if run > qwalk.K_ESC_MAX:                   # priced restart
                        P = E.add(common, E.mul(E.G, int(rng.integers(1, E.p))))
                        run = 0
                        restarts += 1
                        prevK = None
                        continue
                else:
                    run = 0
                if E.dp(K, d):
                    if K in rec and rec[K] != w:
                        solved[s] = True
                    rec.setdefault(K, w)
                prevK = K
                P = E.add(P, E.T[m])
    return solved, fires, restarts


def main():
    t0 = time.time()
    cov = coverage()
    print("coverage:", json.dumps(cov["diff"]), flush=True)

    # (A) full-reach re-run of the floor-bearing cell
    saved = qwalk.Universe.col
    qwalk.Universe.col = col_full
    try:
        full = T.run_arm("W16F", 1 << 28, 2400, R=64, d=5, W=16, seed=20250)
        fullR = T.run_arm("R1024F", 1 << 28, 1600, R=1024, d=5, W=16, seed=20250)
    finally:
        qwalk.Universe.col = saved
    for o in (full, fullR):
        o["selector"] = "full_reach_40bit"

    # (B) EC-coordinate selector arm
    p = prime_1mod3(24)
    while p % 4 != 3:
        p = prime_1mod3(p.bit_length())
        p += 6
        while not (qwalk.is_prime(p) and p % 4 == 3 and p % 3 == 1):
            p += 6
    E = EC(p, 64)
    W, d, sets = 16, 5, 600
    x_shape = 1.0094576719251949          # the shipped shape (L28/L32 value)
    per = int(round(x_shape * math.sqrt(E.classes) / W))
    assert per > 0
    tA = time.time()
    solA, fires, restarts = ec_run(E, sets, W, per, d, seed=20250)
    kA = int(solA.sum())
    solB, f2, r2 = ec_run(E, sets, W, per, d, seed=99991)
    kB = int(solB.sum())
    kS = int((solA | solB).sum())
    ec = {"q": p, "q_log2": p.bit_length() - 1, "curve": "y^2=x^3+7",
          "classes": int(E.classes), "R": E.R, "dp_bits": d, "W": W,
          "sets": sets, "per_walker_steps": per,
          "visits_over_sqrtN": W * per / math.sqrt(E.classes),
          "selector": "EC_coordinate K=min(y,p-y), full reach",
          "p_set_measured": kS / sets, "runs": {"A": kA, "B": kB, "set_solved": kS,
                                                "sets": sets},
          "cp95_set_lower_f": float(C.cp95_lower(kS, sets)),
          "lag3_fires": fires + f2, "restarts": restarts + r2,
          "walks_executed": sets * W * per,
          "elapsed_s": round(time.time() - tA, 1),
          "floor_met_by_cp95": bool(float(C.cp95_lower(kS, sets)) >= 0.39)}

    # matched surrogate arm at the same W/d/SETS and the same shape
    sur = T.run_arm("SURmatch", 1 << 24, sets, R=64, d=d, W=W, seed=20250)
    pA, pB, n = ec["p_set_measured"], sur["p_set_measured"], sets
    se = math.sqrt(pA * (1 - pA) / n + pB * (1 - pB) / n)
    delta = pA - pB

    out = {"test": "selector_dimension", "coverage": cov,
           "full_reach_arms": [full, fullR],
           "ec_selector_arm": ec, "surrogate_matched_arm": sur,
           "swap": {"p_ec": pA, "p_surrogate": pB, "delta": delta,
                    "se_delta": se, "delta_ci95_normal": [delta - 1.96 * se,
                                                          delta + 1.96 * se],
                    "matched": ("same R, W, d, sets, same-challenge w=2 pairing, "
                                "equal dimensionless shape x; EC arm uses 24-bit "
                                "prime j=0 curve, surrogate 24-bit q")},
           "elapsed_s": round(time.time() - t0, 1)}
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("W16F p_set %.5f cp95 %.5f (v18 same cell %.5f)" %
          (full["p_set_measured"], full["cp95_set_lower_f"], 0.54000), flush=True)
    print("EC p_set %.5f cp95 %.5f | surrogate %.5f | delta %+.4f +- %.4f" %
          (pA, ec["cp95_set_lower_f"], pB, delta, 1.96 * se), flush=True)
    print("EC fires %d restarts %d %s" % (ec["lag3_fires"], ec["restarts"],
                                          "floor MET" if ec["floor_met_by_cp95"] else "BELOW"))


if __name__ == "__main__":
    main()
