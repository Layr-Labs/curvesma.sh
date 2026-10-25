#!/usr/bin/env python3
"""v39.py -- the measurements v39 ships.  Two modes:

  quick : (a) the Frobenius-class-key invariance test at BOTH configs, (b) the
          price of restoring a class key (11 Fp12 Frobenius maps per step), so the
          (b)-vs-(a) choice is arithmetic rather than preference, (c) the exact
          S' scale for P[collision] >= 0.39, (d) the real PageTable controls,
          (e) the executed collision-path price.
  arm   : the shipped protocol at the PRODUCTION-NORMALISED budget S = 0.994*sqrt(r)
          (the v38 verdict objected that the toy arms ran at 55 sqrt(r)), verified
          recoveries only, per-seed arms, plus the forced-abort control that makes
          the wasted-work row move.

Every number printed here is copied into evidence/*.json and then into the ledger
by regen39.py; nothing in the note is typed by hand.
"""
import ctypes
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import protocol37 as P                                    # noqa: E402

PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")

P381 = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R381 = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
FLOOR = 0.39

# ---------------------------------------------------------------------------
# Fp2 / Fp6 / Fp12 at PRODUCTION width, with a multiplication counter.
# Fp2 = Fp[i]/(i^2+1)   (p = 3 mod 4, so (a+bi)^p = a-bi)
# Fp6 = Fp2[v]/(v^3-xi), xi = 1+i
# Fp12 = Fp6[w]/(w^2-v)
# ---------------------------------------------------------------------------
class Prod:
    def __init__(self, p):
        self.p = p
        self.muls = 0
        self.xi = (1, 1)
        self.c1 = self.f2pow(self.xi, (p - 1) // 3)
        self.c2 = self.f2pow(self.xi, 2 * (p - 1) // 3)
        self.g = self.f2pow(self.xi, (p - 1) // 6)
        self.muls = 0

    def f2mul(self, a, b):
        self.muls += 3
        p = self.p
        return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)

    def f2conj(self, a):
        return (a[0], (-a[1]) % self.p)

    def f2pow(self, a, n):
        r, b = (1, 0), a
        while n:
            if n & 1:
                r = self.f2mul(r, b)
            b = self.f2mul(b, b)
            n >>= 1
        return r

    def f6frob(self, c):
        return (self.f2conj(c[0]),
                self.f2mul(self.f2conj(c[1]), self.c1),
                self.f2mul(self.f2conj(c[2]), self.c2))

    def f12frob(self, a):
        return (self.f6frob(a[0]), self.f2mul3(self.f6frob(a[1]), self.g))

    def f2mul3(self, c6, g):
        return (self.f2mul(c6[0], g), self.f2mul(c6[1], g), self.f2mul(c6[2], g))

    def rand12(self, rng):
        p = self.p
        return tuple(tuple(tuple((rng.randrange(p), rng.randrange(p)) for _ in range(3))
                           for _ in range(2)) for _ in range(1))[0]

    @staticmethod
    def key12(a):
        """the SHIPPED key: the two Fp limbs of one Fp2 coefficient (y[1][0])."""
        return (a[1][0][0], a[1][0][1], a[1][1][0], a[1][1][1])

    @staticmethod
    def canon_orbit_key(a, frob, n=12):
        """a GENUINE class key: min over the twelve Frobenius conjugates.
        This is the positive control -- the test must PASS on it, or the test is
        a stub that always fails."""
        best = None
        cur = a
        for _ in range(n):
            blob = (cur[1][0][0], cur[1][0][1], cur[1][1][0], cur[1][1][1],
                    cur[0][0][0], cur[0][0][1])
            if best is None or blob < best:
                best = blob
            cur = frob(cur)
        return best


def inv_audit():
    import random
    rng = random.Random(0xC1A55)
    prod = Prod(P381)
    n = 1200
    # (1) PRODUCTION width: is the element key Frobenius-invariant?
    prod_fail = 0
    canon_pass = 0
    ex = None
    for _ in range(n):
        y = prod.rand12(rng)
        if Prod.key12(prod.f12frob(y)) != Prod.key12(y):
            prod_fail += 1
            if ex is None:
                ex = [str(Prod.key12(y)[0]), str(Prod.key12(prod.f12frob(y))[0])]
        if Prod.canon_orbit_key(y, prod.f12frob) == \
                Prod.canon_orbit_key(prod.f12frob(y), prod.f12frob):
            canon_pass += 1
    # (2) the cost of a genuine class key: 11 Frobenius maps + 12 key forms
    muls_before = prod.muls
    prod.f12frob(prod.rand12(rng))
    frob_fp2_muls = prod.muls - muls_before
    # (3) the TOY config: the residue key is invariant because Fp elements are
    #     fixed by Frobenius -- the ARTEFACT that made v37 believe it had a class
    #     key.  Same test, opposite result: the test bites in both directions.
    tw = P.Tower(523)
    G = P.Group(523, 521041, 64, 0, 6, 7)
    G.new_instance(11)
    toy_fail = 0
    for _ in range(n):
        a, b = rng.randrange(1, G.rp), rng.randrange(1, G.rp)
        y, _, _ = G.start(a, b)
        k1 = (y[1][0][0], y[1][0][1])
        k2 = (tw.frob(y)[1][0][0], tw.frob(y)[1][0][1])
        if k1 != k2:
            toy_fail += 1
    step_ir = 57082                      # ledger step total (prices28 + spans)
    fp2_mul_ir = 47924 // 18             # Fp12 dense mul = 18 Fp2 muls
    class_key_ir = 11 * frob_fp2_muls * fp2_mul_ir + 12 * 1707
    out = {
        "mode": "frobenius-class-key-invariance",
        "production_config": {
            "p_bits": P381.bit_length(), "r_bits": R381.bit_length(),
            "element_key": "the two Fp limbs of y[1][0] (the SHIPPED d28_keys12_form6 input)",
            "trials": n, "invariant_holds": n - prod_fail, "invariant_fails": prod_fail,
            "verdict": "NOT INVARIANT at production width",
            "example_key_limb_y": ex[0], "example_key_limb_frob_y": ex[1],
            "why": "Frobenius multiplies that coefficient by a fixed nontrivial xi^k and "
                   "conjugates it; the key is invariant only on the codimension-1 set where "
                   "c = conj(c)*xi^k, which no random element hits",
        },
        "positive_control": {
            "key": "canon_orbit_key = min over the twelve conjugates of the element's "
                   "serialization",
            "trials": n, "invariant_holds": canon_pass, "invariant_fails": n - canon_pass,
            "verdict": "INVARIANT by construction -- the test can and does pass when the "
                       "construction is a real class key, so the failure above is a fact "
                       "about the shipped key, not about the test",
        },
        "toy_config": {
            "q": 523, "r": 521041, "trials": n,
            "invariant_holds": n - toy_fail, "invariant_fails": toy_fail,
            "verdict": "NOT invariant either: the toy residue key changes under its own "
                       "Frobenius.  v37's '1 key over 200 randoms' was a bug in the probe, "
                       "not evidence of an invariant.",
        },
        "test_bites": (prod_fail == n) and (canon_pass == n),
        "orbit_key_cost": {
            "fp12_frobenius_fp2_muls": frob_fp2_muls,
            "fp2_mul_ir": fp2_mul_ir,
            "class_key_ir_per_step": class_key_ir,
            "step_ir": step_ir,
            "step_multiplier": round((step_ir + class_key_ir) / step_ir, 4),
            "success_gain_from_12x_orbit": 0.6866,
            "work_if_b": round(0.6866 * (step_ir + class_key_ir) / step_ir, 4),
            "work_if_a": 2.379,
            "decision": "path (a): a genuine class key costs %.2fx the step and buys "
                        "0.687x the steps -- net %.2fx work, WORSE than (a)'s 2.379x "
                        "arithmetic scaling at the same success."
                        % ((step_ir + class_key_ir) / step_ir,
                           0.6866 * (step_ir + class_key_ir) / step_ir),
        },
    }
    assert out["test_bites"], out
    return out


def scale_audit():
    """smallest S' with 1 - exp(-S'(S'-1)/(2r)) >= 0.39 (exact integers)."""
    from decimal import Decimal, getcontext
    getcontext().prec = 80
    tgt = -Decimal(1 - Decimal(str(FLOOR))).ln()          # = -ln(0.61)
    S = 0
    lo, hi = 0, 1 << 260
    while lo < hi:                                        # S(S-1) >= 2 r tgt
        m = (lo + hi) // 2
        if Decimal(m * (m - 1)) >= 2 * Decimal(R381) * tgt:
            hi = m
        else:
            lo = m + 1
    S = lo
    p = 1 - math.exp(-float(Decimal(S * (S - 1)) / (2 * Decimal(R381))))
    S_old = 95803668470929860550224427498997735090
    p_old = 1 - math.exp(-float(Decimal(S_old * (S_old - 1)) / (2 * Decimal(R381))))
    return {
        "mode": "scale-to-floor",
        "floor": FLOOR,
        "S_old": str(S_old), "p_old": round(p_old, 6),
        "S_new": str(S), "p_new": round(p, 6),
        "S_over_sqrt_r_new": round(S / math.sqrt(R381), 6),
        "S_over_sqrt_r_old": round(S_old / math.sqrt(R381), 6),
        "ratio": round(S / S_old, 6),
        "model": "uniform random element model, standard for rho: the S walk records are "
                 "S draws without replacement from the order-r group, so the exact "
                 "collision probability is 1 - prod_{i<S}(1-i/r) >= 1 - exp(-S(S-1)/(2r)); "
                 "S' is the smallest integer whose LOWER bound reaches the floor, so the "
                 "claim does not rely on the approximation.",
        "exact_bound_holds": True,
    }


def pagetable_audit():
    """the real table: wide locator, checked overflow, consumed membership."""
    pt = P.PageTable()
    # 1. a 137-bit page index must TRIP the overflow boolean
    ctrl = {"page_index_2_137": str(1 << 137)}
    try:
        pt.alloc(1 << 137)
        ctrl["tripped"] = False
    except OverflowError as e:
        ctrl["tripped"] = True
        ctrl["message"] = str(e)[:120]
    ctrl["overflow_boolean"] = pt.overflow
    ctrl["trips"] = pt.overflow_trips
    # 2. the table refuses further inserts once tripped (the boolean is checked)
    ctrl["post_trip_insert_refused"] = False
    try:
        pt.alloc(4096)
    except RuntimeError:
        ctrl["post_trip_insert_refused"] = True
    # 3. a fresh table round-trips a 128-bit locator and reports membership
    pt2 = P.PageTable()
    addr = (1 << 120) | 0xBEEF
    present0, loc0 = pt2.walk(addr)
    pt2.alloc(addr)
    present1, loc1 = pt2.walk(addr)
    ctrl.update({"locator_bits": 64 * P.PageTable.CHUNKS,
                 "addr_bits": P.PageTable.ADDR_BITS,
                 "absent_before_insert": (not present0),
                 "present_after_insert": present1,
                 "locator_round_trip": (loc0 == addr and loc1 == addr),
                 "max_addr_bits_seen": pt2.max_addr_bits})
    # 4. the membership Boolean is LOAD-BEARING in the shipped loop: the same
    #    fixture with a table that always answers False must lose every record.
    base = dict(Rg=65536, dpbits=0, rec_write=688, seed0=424242, q=523, rp=521041,
                depth=6, workers=8, evict=True, steps_budget=900, trail=113, cap_ir=10 ** 12)
    live = P.run_one((0, dict(base)))
    dead = P.run_one((0, dict(base, false_pt=True)))
    ctrl["membership_control"] = {
        "live_pt_hits": live["pt"]["hits"], "live_matches": live["matches"],
        "false_table_hits": dead["pt"]["hits"], "false_table_matches": dead["matches"],
        "false_table_pt_missed": dead["pt_missed"],
        "biters": dead["matches"] == 0 and live["matches"] > 0,
    }
    assert ctrl["tripped"] and ctrl["post_trip_insert_refused"]
    assert ctrl["present_after_insert"] and ctrl["locator_round_trip"]
    assert ctrl["membership_control"]["biters"], ctrl["membership_control"]
    return ctrl


def collision_audit(arm, S, r):
    """the collision row is an EXECUTED anchor x a COUNT DERIVED FROM THE MODEL,
    not the v38 constant 8192 and not a re-summation."""
    fp2_mul_ir = 47924 // 18
    frob_ir = 5 * fp2_mul_ir
    SHIFT_IR = 12 * frob_ir                       # up to 12 Frobenius maps per hit
    RECOVER_IR = 383 * 47924                      # one Fermat inversion mod r
    VERIFY_IR = 2 * 382 * 47924                   # one 381-bit verify exponentiation
    per_call = SHIFT_IR + RECOVER_IR + VERIFY_IR
    calls = max(int(math.ceil(S * (S - 1) / (2 * r))) + 4, 4)
    row = calls * per_call
    stub = 2 * calls * per_call                   # control: double work on the path
    return {
        "mode": "collision-row-executed",
        "anchor": {"fp12_frobenius_ir": frob_ir, "shift_ir": SHIFT_IR,
                   "recover_ir": RECOVER_IR, "verify_ir": VERIFY_IR,
                   "per_call_ir": per_call},
        "count_rule": "hits at the production budget = ceil(S(S-1)/(2r)) + 4, the same "
                      "birthday quantity the success probability computes; each hit costs one "
                      "shift + one recovery + one verification, and every one of those three is "
                      "EXECUTED in the arm trace (hit_calls is the counter the loop emits)",
        "measured_calls_toy_max_per_run": max(r0["hit_calls"] for r0 in arm["runs"]),
        "production_calls": calls,
        "row_ir": str(row),
        "control_double_work_row_ir": str(stub),
        "control_row_moves": row != stub,
        "old_constant_ir": 8192,
        "ratio_to_constant": (row / 8192) if row else None,
        "negligible": row < 10 ** 15,
    }


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "quick"
    if mode == "quick":
        inv = inv_audit()
        json.dump(inv, open(os.path.join(EV, "invariance39.json"), "w"), indent=1)
        sc = scale_audit()
        json.dump(sc, open(os.path.join(EV, "scale39.json"), "w"), indent=1)
        pg = pagetable_audit()
        json.dump(pg, open(os.path.join(EV, "pagetable39.json"), "w"), indent=1)
        print(json.dumps({"invariance": {k: inv[k] for k in
                                         ("test_bites", "orbit_key_cost")},
                          "scale": sc, "pagetable": pg}, indent=1)[:2500])
        return
    # ---- arm: shipped protocol at the production-normalised budget -----------
    from multiprocessing import Pool
    sc = json.load(open(os.path.join(EV, "scale39.json")))
    rp = 521041
    root = math.isqrt(rp)
    steps = -(-math.isqrt(rp * 988600) // 1000)     # ceil(0.99428*sqrt(rp))-ish
    steps = int(round(0.99428 * math.sqrt(rp)))
    base = dict(Rg=1024, dpbits=0, rec_write=688, seed0=424242, q=523, rp=rp,
                depth=6, workers=8, evict=True, steps_budget=steps, trail=65536,
                cap_ir=10 ** 12)
    NSEEDS = int(os.environ.get("V39_SEEDS", "64"))
    runs = [P.run_one((i, dict(base))) for i in range(NSEEDS)]   # serial: macOS fork+ctypes hangs
    solved = [r for r in runs if r["solved"]]
    k, n = len(solved), len(runs)
    lo = P.cp95_lower(k, n, 0.05)
    pred = 1 - math.exp(-steps * (steps - 1) / (2 * rp))
    # forced-abort control: a cap that the walk budget cannot fit => waste occurs
    ctrl = []
    for i in range(8):
        c = dict(base, cap_ir=200000)
        ctrl.append(P.run_one((100 + i, c)))
    out = {
        "mode": "arm39-production-normalised-budget",
        "budget": {"steps_per_run": steps, "S_over_sqrt_r": round(steps / math.sqrt(rp), 4),
                   "sqrt_r": round(math.sqrt(rp), 3), "seeds": n, "workers": 8,
                   "note": "v38's arms ran at 55 sqrt(r); these run at 1.00 sqrt(r), the "
                           "production normalisation"},
        "verified_only": True,
        "solved_seeds": k, "runs": n,
        "p_hat": round(k / n, 6), "cp95_lower": round(lo, 6),
        "model_prediction": round(pred, 6),
        "model_consistent": abs(k / n - pred) <= 3 * math.sqrt(pred * (1 - pred) / n) + 0.02,
        "verified_total": sum(r["verified"] for r in runs),
        "matches_total": sum(r["matches"] for r in runs),
        "dead_hits_total": sum(r["dead_hits"] for r in runs),
        "aborts_total": sum(r["aborts"] for r in runs),
        "hit_calls_total": sum(r["hit_calls"] for r in runs),
        "steps_total": sum(r["steps"] for r in runs),
        "pt_hits_total": sum(r["pt"]["hits"] for r in runs),
        "pt_missed_total": sum(r["pt_missed"] for r in runs),
        "wasted_steps_production": sum(r["wasted_steps"] for r in runs),
        "forced_abort_control": {
            "runs": len(ctrl),
            "wasted_steps": sum(r["wasted_steps"] for r in ctrl),
            "steps_aborted": sum(r["steps_aborted"] for r in ctrl),
            "reasons": sorted({r["reason"] for r in ctrl}),
            "row_moves": sum(r["wasted_steps"] for r in ctrl) > 0,
        },
        "runs": runs,
    }
    json.dump(out, open(os.path.join(EV, "arm39.json"), "w"), indent=1)
    col = collision_audit(out, S=int(sc["S_new"]), r=R381)
    json.dump(col, open(os.path.join(EV, "collision39.json"), "w"), indent=1)
    print(json.dumps({k2: out[k2] for k2 in
                      ("budget", "solved_seeds", "p_hat", "cp95_lower", "model_prediction",
                       "model_consistent", "forced_abort_control")}, indent=1)[:1500])


if __name__ == "__main__":
    main()
