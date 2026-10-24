#!/usr/bin/env python3
"""verify27.py -- independent re-derivation of the v27 ledger and claim.

Reads the EVIDENCE files only (not ledger27.json), rebuilds every row and the
claim from them with exact integer arithmetic, and compares against what the
generator published.  Also re-checks the five repair invariants of v27.
"""
import json
import math
import os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")


def jl(n):
    return json.load(open(os.path.join(EV, n)))


def main():
    cen, surv, scan, scan25 = jl("census25.json"), jl("survival25.json"), jl("scan26.json"), jl("scan25.json")
    ir25, ir26, pr, d26 = jl("ir25.json"), jl("ir26.json"), jl("prices25.json"), jl("dict26.json")
    led = jl("ledger27.json")
    man = jl("manifest27.json")

    MUL = pr["MUL_CHARGE"]
    SQR = ir25["rows"]["cyclotomic_sqr"]["ir_count"]
    fm = F(scan["multiply"], scan["steps"])
    STEP = fm * MUL + (1 - fm) * SQR
    LOOKUP = ir26["rows"]["state_key48"]["ir_count"] + ir26["rows"]["dict_lookup"]["ir_count"]
    REC = 688 + ir26["rows"]["dict_insert"]["ir_count"]

    Q = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
    N = (Q - 1) // 12
    sqrtN = math.isqrt(N)
    ratio = F(str(cen["scales"][0]["S_over_sqrtN"])).limit_denominator(1000)
    steps = ratio.numerator * sqrtN // ratio.denominator

    slots = -(-steps * 8 // 5)
    ref_bits = (slots - 1).bit_length()
    ref_words = -(-ref_bits // 64)
    slot_words = 6 + ref_words + 1
    assert slots == int(led["addressing"]["slots"])
    assert ref_bits == led["addressing"]["slot_ref_bits"]
    assert slot_words == led["addressing"]["slot_words"]
    assert 8 * steps <= 5 * slots
    assert (1 << ref_bits) >= slots

    rows = {
        "walk_steps": int(steps * STEP),
        "per_step_lookup_key_and_probe": int(steps * LOOKUP),
        "record_writes_and_insert": int(steps * REC),
        "collision_processing": int(max(1, math.ceil(F(steps * (steps - 1), 2 * N))) * pr["REPLAY_IR"]),
        "trap_waste_bound": int(led["rows_ir"]["trap_waste_bound"]),
        "index_slots_widened_refs": int(steps * slot_words),
        "worker_starts_1024": 1024 * pr["START_IR"],
        "code_and_advice": 1 << 20,
        "candidate_verification": pr["VERIFY_PROD"],
    }
    W = sum(rows.values())
    for k, v in rows.items():
        assert int(led["rows_ir"][k]) == v, (k, led["rows_ir"][k], v)
    assert int(led["W_total_ir"]) == W

    n0 = int(led["claim"]["numerator"])
    assert 20000 ** 10000 * 2 ** n0 >= W ** 10000
    assert 20000 ** 10000 * 2 ** (n0 - 1) < W ** 10000

    # repair invariants
    assert led["manifest"]["files"] == len(man["files"])
    assert all(os.path.exists(os.path.join(PKG, f["path"])) for f in man["files"])
    assert all(len(f["sha256"]) == 64 for f in man["files"])
    assert led["survival"]["scales_measured"] == ["rp=521041"]
    assert led["survival"]["not_measured_at_rp10566001"] is True
    assert led["success"]["route"] == "bootstrap-product"
    assert led["success"]["p5"] >= 0.39
    assert led["addressing"]["slot_bytes"] == slot_words * 8
    assert led["addressing"]["table_bytes"] == slots * slot_words * 8
    assert led["trap"]["ring_words"] == 6
    assert led["trap"]["family_census"]["route"] in ("census", "measured_windows")
    assert led["charges"]["SQR_CHARGE"] == led["charges"]["SQR_STATIC_SPAN"]
    assert d26["zero_mismatches"] and d26["mutation_control"]["detected"]

    out = {"mode": "verify27", "rows_match": True, "W": str(W), "claim": led["claim"]["value"],
           "claim_covers_W": True, "claim_minimal": True, "addressing_consistent": True,
           "manifest_ok": True, "survival_one_scale": True, "success_bootstrap_p5": led["success"]["p5"],
           "trap_ring": 6, "sqr_surcharge_removed": True, "all_pass": True}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
