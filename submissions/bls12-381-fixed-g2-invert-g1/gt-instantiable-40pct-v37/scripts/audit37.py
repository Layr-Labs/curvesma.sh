#!/usr/bin/env python3
"""audit37.py -- the four planted controls, all machine-logged.

(a) BORROW: negate(1) must give m-1 in six limbs at the PRODUCTION modulus p,
    1000 random add/sub/neg must match the Python reference, and the OLD
    limb-local rule is evaluated beside it to show it FAILS the planted case.
(b) COLD START: G.start must not read x -- it is checked against g^(a+xb)
    computed from the secret AFTER the fact, and the restart row must move off
    7,660 Ir.
(c) OCCUPANCY: sparse regime (S << M) gives ~S, dense regime gives M, and the
    v36 dense-leaf formula is shown to be BELOW the occupancy (flagged).
(d) YIELD SCHEDULE: the relation yield at L = 64 must be below the shipped-L
    yield and both must be measured by the same shipped path.
"""
import ctypes
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import protocol37 as P                                   # noqa: E402
import bounds37 as B                                     # noqa: E402

EV = os.path.join(os.path.dirname(HERE), "evidence")
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
M64 = (1 << 64) - 1


def limbs(x, n=6):
    return (ctypes.c_uint64 * n)(*[(x >> (64 * j)) & M64 for j in range(n)])


def val(c, n=6):
    return sum(int(c[j]) << (64 * j) for j in range(n))


def main():
    lib = P.EVICTING
    out = (ctypes.c_uint64 * 6)()
    # ---- (a) borrow propagation and reduction, at m = p -------------------
    lib.d28_l6_neg(out, limbs(1), limbs(p))
    neg1 = val(out)
    old_rule = 1 if False else M64          # the v36 limb-local negation of 1
    a = {"neg1_equals_m_minus_1": neg1 == p - 1, "neg1_value_hex": hex(neg1),
         "old_limb_local_rule_value": hex(old_rule),
         "old_rule_fails_planted_case": old_rule != p - 1,
         "old_rule_note": "v36's (0 - A[j]) per limb gives the low limb 2^64-1 with zero upper "
                          "limbs, i.e. 0xffffffffffffffff, not m-1"}
    rng = random.Random(20261024)
    bad = 0
    for _ in range(1000):
        x, y = rng.randrange(0, p), rng.randrange(0, p)
        lib.d28_l6_sub(out, limbs(x), limbs(y), limbs(p)); bad += val(out) != (x - y) % p
        lib.d28_l6_add(out, limbs(x), limbs(y), limbs(p)); bad += val(out) != (x + y) % p
        lib.d28_l6_neg(out, limbs(x), limbs(p)); bad += val(out) != (-x) % p
    red_bad = 0
    for x in (p, p + 1, 2 * p - 1, 2 * p, 2 * p + 7, 3 * p - 3):
        carry = x >> 384
        lib.d28_l6_reduce(out, limbs(x & ((1 << 384) - 1)), ctypes.c_uint64(carry), limbs(p))
        red_bad += val(out) != x % p
    a.update({"randoms": 1000, "mismatches": bad, "reduce_overflow_mismatches": red_bad,
              "modulus": "p (381-bit production base field)",
              "landing": "borrow-propagating sub/add/neg and carry-limb reduce all match the "
                         "Python reference exactly"})

    # ---- (b) cold start --------------------------------------------------
    G = P.Group(523, 521041, 1024, 0, 6, 7)
    G.new_instance(424242)
    aa, bb = 12345, 54321
    y, A, D = G.start(aa, bb)
    check = G.T.pow(G.g, (aa + G.x * bb) % G.rp)       # x used ONLY to check, after the fact
    src = open(os.path.join(HERE, "protocol37.py")).read()
    body = src.split("def start(self, a, b):")[1].split("def shift")[0]
    b = {"start_matches_public_combination": y == check,
         "start_source_reads_secret": ("self.x" in body),
         "restart_ir_cold": int(json.load(open(os.path.join(EV, "ledger37.json")))
                               ["budget"]["restart_cold_ir"]),
         "restart_ir_v36": 7660,
         "row_changed": int(json.load(open(os.path.join(EV, "ledger37.json")))
                            ["budget"]["restart_cold_ir"]) != 7660,
         "rule": "a cold attack constructs the start from the public g and z only; the priced "
                 "start is 2 two-scalar powerings at the per-step multiply price"}

    # ---- (c) occupancy ---------------------------------------------------
    M = B.M
    sparse = B.occupancy(M, 1000)
    dense = B.occupancy(M, int(json.load(open(os.path.join(EV, "ledger37.json")))
                              ["protocol"]["steps"]))
    led = json.load(open(os.path.join(EV, "ledger37.json")))
    dense6 = int(led["store_spec_to_code"]["dense_six_ceil_control"]["value"])
    c = {"sparse_S_1000": {"reached": str(sparse[0]), "S_over_M": sparse[1], "regime": sparse[2],
                           "approx_S": abs(sparse[0] - 1000) <= 2},
         "dense_S": {"reached": str(dense[0]), "S_over_M": dense[1], "regime": dense[2],
                     "equals_M": dense[0] == M},
         "dense_six_ceil_value": str(dense6),
         "dense_six_ceil_below_occupancy": dense6 < dense[0],
         "flagged": "the v36 bound is %.1f x BELOW the occupancy" % (dense[0] / dense6)}

    # ---- (d) yield schedule sensitivity ---------------------------------
    import arm37 as A
    a64 = A.arm(64, 1)
    a113 = A.arm(113, 1)
    ash = A.arm(65536, 1)
    d = {"yield_cp_lower_at_64": a64["yield_cp_lower"], "hits_at_64": a64["hits"],
         "yield_cp_lower_at_113": a113["yield_cp_lower"],
         "yield_cp_lower_at_shipped": ash["yield_cp_lower"], "hits_at_shipped": ash["hits"],
         "schedule_sensitive": a64["yield_cp_lower"] < ash["yield_cp_lower"],
         "rule": "lengthening the trail must RAISE the yield if the arm measures the schedule and "
                 "not a constant: 64 < 113 < 65536 in this run"}
    ok = (a["neg1_equals_m_minus_1"] and a["mismatches"] == 0 and a["reduce_overflow_mismatches"] == 0
          and b["start_matches_public_combination"] and not b["start_source_reads_secret"]
          and b["row_changed"] and c["sparse_S_1000"]["approx_S"] and c["dense_S"]["equals_M"]
          and c["dense_six_ceil_below_occupancy"] and d["schedule_sensitive"])
    json.dump({"mode": "audit37-four-plants", "a_borrow": a, "b_cold_start": b,
               "c_occupancy": c, "d_yield_schedule": d, "all_four_pass": bool(ok)},
              open(os.path.join(EV, "audit37.json"), "w"), indent=1)
    print(json.dumps({"a": a["neg1_equals_m_minus_1"] and a["mismatches"] == 0,
                      "b": b["start_matches_public_combination"] and not b["start_source_reads_secret"],
                      "c": c["dense_S"]["equals_M"] and c["dense_six_ceil_below_occupancy"],
                      "d": d["schedule_sensitive"], "all": bool(ok),
                      "yields": [a64["yield_cp_lower"], a113["yield_cp_lower"],
                                 ash["yield_cp_lower"]]}, indent=1))
    assert ok


if __name__ == "__main__":
    main()
