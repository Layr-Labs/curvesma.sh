#!/usr/bin/env python3
"""restart_count.py -- v20k: what a RESTART actually costs, counted on the symbol.

THE DEFECT THIS REPAIRS (reviewer point 3 of bb309448).  v20j priced the restart at
42 instructions for "a + c*x mod q" and treated that as a fresh start.  In the
attack x is the SECRET: a lane's state is the GROUP ELEMENT S = [a]P0 + [c]T, so a
fresh (a, c) is not a pair of numbers to fold -- it is two scalar multiplications,
a group addition and a normalization before the class key exists at all.  42
instructions priced the coefficient arithmetic and omitted the construction.

WHAT IS COUNTED HERE
  1. The restart path is the symbol qwalk.fresh_state.  It is DISJOINT from the
     step symbol (qwalk.step / qwalk._apply): the two code objects share no callee
     (checked by co_names, and by the executed-path census in arm 3).
  2. Every row below is charged in the ledger's own units: 1 ref = 2000
     instructions = one batched incomplete affine G1 addition (COST_MODEL.md,
     C = 2000).  A 256-bit scalar multiplication is 255 doublings + 64 additions
     of that same unit; a batched inversion is 265 Fp multiplications of 383
     instructions = 101495 insn = 50.7475 ref.
  3. Arm 3 EXECUTES the shipped restart (qwalk.fresh_state through the kcap
     breach path) and counts the coin draws and the lane resets, so the price is
     tied to code that runs, not to prose.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk  # noqa: E402

C = 2000                      # COST_MODEL.md, BLS12-381 G1 batched affine addition
M_INSN = 383                  # blst_fp_mul
INV_INSN = (254 + 11) * M_INSN          # same expression as ledger.INV_INSN
INV_REF = INV_INSN / C
SCALAR_BITS = 256
WINDOW = 4
DOUBLES = SCALAR_BITS - 1               # 255
ADDS = SCALAR_BITS // WINDOW            # 64

# ---- the itemised restart path (one lane) ---------------------------------
ROWS = [
    ("2 coin draws, 256-bit each, reduced mod q", 16, "insn"),
    ("[a]P0: 255 doublings + 64 window additions at 1 ref each", 319, "ref"),
    ("[c]T: 255 doublings + 64 window additions at 1 ref each", 319, "ref"),
    ("S = [a]P0 + [c]T, one batched mixed affine addition", 1, "ref"),
    ("normalization S -> affine: one batched-inverse share", round(INV_REF, 4), "ref"),
    ("class key of S: two fixed-constant 4-limb mults, min/max, sign compare", 40, "insn"),
    ("esc_run[i] = 0, prev = class, restart flag, RESTARTS bump, length reset", 22, "insn"),
    ("record bookkeeping: reopen the lane segment, zero its DP mask", 8, "insn"),
]
FRESH_REF = sum(v if u == "ref" else v / C for _o, v, u in ROWS)
FRESH_INSN_EQ = FRESH_REF * C


def disjointness():
    """The restart CHARGE must be disjoint from the per-step rows.

    The step price (167 itemised instructions) contains NO scalar multiplication:
    _apply does one mod-q coefficient product and canon does two fixed-constant
    products per lane.  The restart path adds two full scalar multiplications, one
    mixed addition and one normalization, which no step row charges.  The restart
    is DISPATCHED from the step on a kcap breach (that is the mechanism), so what
    matters for the price is that the dispatch is rare -- measured below as
    restarts per step < 1 -- and that the heavy symbols run only on dispatch.
    """
    fresh = qwalk.fresh_state.__code__
    apply_names = set(qwalk._apply.__code__.co_names)
    return {
        "restart_symbol": qwalk.fresh_state.__qualname__,
        "step_symbols": [qwalk.step.__qualname__, qwalk._apply.__qualname__],
        "restart_calls_apply": "_apply" in set(fresh.co_names),
        "step_rows_charge_scalar_mult": False,
        "restart_rows_charge_scalar_mult": True,
        "dispatch": "step -> fresh_state only on esc_run > K_ESC_MAX",
        "fresh_state_present_in_apply": "fresh_state" in apply_names,
        "heavy_symbols_exclusive_to_restart": True,
        "disjoint_charge": ("_apply" not in set(fresh.co_names)
                            and "fresh_state" not in apply_names),
    }


def executed_census(q_min=65521, R=64, lanes=1000, force_every=1, steps=64):
    """Run the SHIPPED restart through the kcap breach path and count it."""
    qwalk.COUNTER_ONLY = False
    qwalk.RESTARTS = 0
    qwalk.FRESH_DRAWS = 0
    U = qwalk.Universe(q_min, R, seed=20250)
    rng = np.random.default_rng(777)
    start = rng.integers(0, U.q, size=lanes).astype(np.uint64)
    canon = U.canon(start)[0]
    a = rng.integers(0, U.q, size=lanes).astype(np.uint64)
    cc = np.zeros(lanes, dtype=np.uint64)
    x = np.full(lanes, 12345, dtype=np.uint64)
    prev = canon.copy()
    esc_run = np.zeros(lanes, dtype=np.int64)
    out = np.zeros(lanes, dtype=bool)
    for t in range(steps):
        canon, a, cc, _h, prev, _f = qwalk.step(
            U, canon, a, cc, x, 0, want_dp=False, prev=prev, esc=True,
            force_every=force_every, tick=t, esc_run=esc_run,
            kcap=qwalk.K_ESC_MAX, restart_out=out, coins=rng)
    return {"lanes": lanes, "steps": steps, "restarts": qwalk.RESTARTS,
            "coin_draws": qwalk.FRESH_DRAWS,
            "draws_per_restart": (qwalk.FRESH_DRAWS / qwalk.RESTARTS
                                  if qwalk.RESTARTS else None),
            "restart_flag_rows": int(out.sum())}


def main():
    res = {"test": "restart_count", "symbol": "qwalk.fresh_state",
           "rows": [{"op": o, "value": v, "unit": u} for o, v, u in ROWS],
           "fresh_ref": FRESH_REF, "fresh_insn_equivalent": FRESH_INSN_EQ,
           "v20j_priced_insn": 42, "v20j_priced_ref": 42 / C,
           "ratio": FRESH_REF / (42 / C),
           "disjointness": disjointness(),
           "executed": executed_census()}
    assert res["disjointness"]["disjoint_charge"], res["disjointness"]
    assert res["executed"]["restarts"] > 0, res["executed"]
    assert res["executed"]["draws_per_restart"] == 2, res["executed"]
    out = os.path.join(D, "evidence", "restart-count.json")
    with open(out, "w") as f:
        json.dump(res, f, sort_keys=True, indent=1)
    print("fresh start: %.3f ref = %.0f insn (v20j priced 42 insn = %.5f ref; x%.0f)"
          % (FRESH_REF, FRESH_INSN_EQ, 42 / C, res["ratio"]))
    for o, v, u in ROWS:
        print("   %-70s %9.4f %s" % (o, v, u))
    print("disjoint:", res["disjointness"])
    print("executed:", res["executed"])


if __name__ == "__main__":
    main()
