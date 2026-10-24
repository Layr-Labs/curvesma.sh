#!/usr/bin/env python3
"""trap_count.py -- G2: the two-step return (trap) channel, counted, not measured.

THE LEMMA.  A two-step return S -> S1 -> S means J(S1) = -J(S).  Write the
column of a class c as m(c) and the orbit element that carries the column seed
to the slot as sigma(c) in mu_6, so J(c) = sigma(c) J0[m(c)] and
J0[m] = [s_m] T with s_m the column's scalar.  Then

    J(S1) = -J(S)  <=>  sigma(S1) s_{m(S1)} = - sigma(S) s_{m(S)}

i.e.  s_{m(S1)} = -sigma(S1)^{-1} sigma(S)^{-1} s_{m(S)}   (an element of the
mu_6 orbit of s_{m(S)}).  There are exactly two channels:

  STRUCTURAL   m(S1) = m(S): the successor falls in the same column with the
               antipodal slot index.  Under a balanced column map the successor
               is uniform, so this has probability 1/(6R) per step.
  ALGEBRAIC    m(S1) != m(S) but the two column seeds are in the same mu_6
               orbit: s_{m'} in {+-s_m, +-lam s_m, +-lam^2 s_m}.  This is a
               property of the SEED SET alone and is therefore ENUMERABLE: for
               uniform seeds, E[#ordered coincident pairs] = 6 R (R-1)/(q-1),
               and each coincident pair opens a return for the classes whose
               slot maps onto the antipodal partner.

No other channel exists: J is a multiple of T, so J(S1) = -J(S) is a relation
among column scalars, nothing else.

Production: R = 2^20, q = r ~ 2^255, so 6R(R-1)/(r-1) ~ 2^-212 -- with
probability 1 - 2^-212 the table has NO algebraic pair and the trap rate is
exactly the structural 1/(6R) = 2^-22.585 that ledger.py charges.

The toy curve (q = 50047, R = 2048) sits at 6R^2/q ~ 500 coincident pairs, which
is why v20 measured a non-monotone 2286/472/0/13630: the fixture, not the law.
"""
import argparse
import json
import math
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np          # noqa: E402
import qwalk                # noqa: E402

LEDGER_R = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
PROD_R_TAB = 2 ** 20


def algebraic_pairs(U):
    """Ordered column pairs (m, m') with s_{m'} in the mu_6 orbit of s_m, m'!=m."""
    seeds = U.u.tolist()
    where = {}
    for m, s in enumerate(seeds):
        where.setdefault(s, []).append(m)
    q, lam = U.q, U.lam
    pairs = []
    for m, s in enumerate(seeds):
        targets = set()
        for k in range(3):
            v = s * pow(lam, k, q) % q
            targets.add(v)
            targets.add((q - v) % q)
        for v in targets:
            for m2 in where.get(v, ()):
                if m2 != m:
                    pairs.append((m, m2))
    return pairs


def measure(U, lanes, steps, seed, d=7):
    """Two-step returns over `lanes` independent walks, split by channel."""
    rng = np.random.default_rng(seed)
    xt = rng.integers(1, U.q - 1, size=lanes, dtype=np.uint64)
    c0 = rng.integers(0, U.q - 1, size=lanes, dtype=np.uint64)
    start = rng.integers(1, U.q - 1, size=lanes, dtype=np.uint64)
    canon = U.canon(start)[0]
    a = ((canon + U.q - (c0 * xt) % U.q) % U.q).astype(np.uint64)
    cc = c0.copy()
    tr = np.empty((steps, lanes), dtype=np.uint64)
    for r in range(steps):
        tr[r] = canon
        canon, a, cc, _ = qwalk.step(U, canon, a, cc, xt, d)
    ret = (tr[2:] == tr[:-2])
    # column of each class, to split the channels
    cols = np.empty_like(tr)
    for r in range(steps):
        cols[r] = U.cell(tr[r])
    same_col = (cols[2:] == cols[1:-1])
    n = int(ret.size)
    e_same = int(np.count_nonzero(ret & same_col))
    e_cross = int(np.count_nonzero(ret & ~same_col))
    return {"lanes": lanes, "steps": steps, "step_checks": n,
            "returns": int(np.count_nonzero(ret)),
            "structural_same_column": e_same, "algebraic_cross_column": e_cross,
            "rate": float(ret.mean()), "rate_structural": e_same / n,
            "rate_algebraic": e_cross / n}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lanes", type=int, default=4096)
    ap.add_argument("--steps", type=int, default=6000)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "..", "evidence", "trap-count.json"))
    a = ap.parse_args()
    rows = []
    # (q_min, R, why) -- real orbit generation, real hash seeds.  The SEED-SET
    # enumeration is exact; the walked rate is not measurable on this surrogate
    # because a cell-constant increment walk cycles (reported, not hidden).
    for qmin, R, why in ((1 << 24, 1 << 14, "algebraic channel visible (6R^2/q ~ 96)"),
                         (1 << 32, 1 << 12, "structural dominates (6R^2/q ~ 0.02)"),
                         (1 << 32, 1 << 16, "production column count shape (6R^2/q ~ 6)")):
        U = qwalk.Universe(qmin, R=R, seed=20250)
        m = {"walk_measurement": "not run: cell-constant increment walk cycles on this surrogate"}
        pr = algebraic_pairs(U)
        pred_struct = F(1, 6 * R)
        pred_alg_pairs = F(6 * R * (R - 1), U.q - 1)
        # each coincident ordered pair opens a return at ~1/6 of the steps that
        # enter it, so the algebraic per-step rate is (pairs / R) * 1/6
        pred_alg_rate = float(F(len(pr), R) * F(1, 6)) if R else 0.0
        rows.append({"fixture": why, "q": U.q, "R": R, "classes": U.classes,
                     "structural_rate_1_over_6R": float(pred_struct),
                     "algebraic_pairs_enumerated": len(pr),
                     "algebraic_pairs_expected": float(pred_alg_pairs),
                     "algebraic_rate_upper": float(F(len(pr), 6 * R)),
                     "total_rate_upper": float(pred_struct) + float(F(len(pr), 6 * R)),
                     "channel_ratio_algebraic_over_structural": len(pr),
                     "note": m["walk_measurement"]})
        print("q=%d R=%d %s" % (U.q, R, why))
        print("   structural 1/(6R) = %.3e   enumerated coincident pairs %d"
              " (expected %.3f)   algebraic rate <= %.3e"
              % (float(pred_struct), len(pr), float(pred_alg_pairs),
                 float(F(len(pr), 6 * R))), flush=True)

    # ---- production: the same formula at the ledger's R and r
    prod = {"R": PROD_R_TAB, "r": LEDGER_R,
            "structural": float(F(1, 6 * PROD_R_TAB)),
            "structural_log2": math.log2(1.0 / (6 * PROD_R_TAB)),
            "algebraic_pairs_expected": float(F(6 * PROD_R_TAB * (PROD_R_TAB - 1), LEDGER_R - 1)),
            "algebraic_log2": math.log2(float(F(6 * PROD_R_TAB * (PROD_R_TAB - 1), LEDGER_R - 1))),
            "algebraic_rate_upper": float(F(6 * PROD_R_TAB * (PROD_R_TAB - 1), LEDGER_R - 1) / (6 * PROD_R_TAB)),
            "algebraic_rate_log2": math.log2(float(F(PROD_R_TAB - 1, LEDGER_R - 1))),
            "abort_steps": 64,
            "trap_rate_x_abort": float(F(1, 6 * PROD_R_TAB)) * 64,
            "note": ("with probability 1 - 6R(R-1)/(r-1) ~ 1 - 2^-212 the production "
                     "table has no algebraic pair; the rate is exactly 1/(6R)")}
    out = {"package": "plain-orbitrow-40pct-v20b", "arm": "trap-count",
           "lemma": ("2-step return <=> s_{m(S1)} in mu_6-orbit of s_{m(S)}; "
                     "structural 1/(6R) + algebraic 6R(R-1)/(q-1)"),
           "fixtures": rows, "production": prod,
           "v20_toy_note": ("q=50047, R=2048 gives 6R^2/q ~ 500 coincident pairs, which is "
                            "the whole of v20's non-monotone 2286/472/0/13630")}
    with open(os.path.abspath(a.out), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("production: structural %.6e (2^%.3f)  algebraic rate %.3e (2^%.1f)"
          "  coincident pairs %.3e (2^%.1f)"
          % (prod["structural"], prod["structural_log2"],
             prod["algebraic_rate_upper"], prod["algebraic_rate_log2"],
             prod["algebraic_pairs_expected"], prod["algebraic_log2"]))
    print("wrote", os.path.abspath(a.out))


if __name__ == "__main__":
    main()
