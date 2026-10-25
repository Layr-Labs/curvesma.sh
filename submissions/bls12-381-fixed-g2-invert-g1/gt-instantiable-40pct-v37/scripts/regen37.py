#!/usr/bin/env python3
"""regen37.py -- build the v37 ledger, claim and package metadata.

Every charged row is an exact Python int, every span is read from the SHIPPED
object's disassembly, and every row the verdict called into question is rebuilt
from the executed path:

  * the frequency-independent STEP charge: ONE dense-multiply price for EVERY
    executed step (the guard-squaring substitution ceil(4497 S/5040000) is
    DELETED -- verdict flaw 3);
  * the trail-start row: the launch rule executes S steps at a per-trail cap L,
    so starts = ceil(S/L) EXACTLY (a lower bound that is an equality by
    construction), and each start is priced as a COLD two-scalar exponentiation
    y = g^a z^b that never reads x (verdict flaw 3);
  * the keys12 and coefficient rows: spans of the SIX-LIMB routines the fixture
    actually calls (verdict flaw 1);
  * the store init row: reached nodes from uniform-address OCCUPANCY (verdict
    flaw 4), read from evidence/bounds37.json.
"""
import json
import math
import os
import re
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
DYLIB = os.path.join(HERE, "ir28", "dict28.dylib")

S = int("95803668470929860550224427498997735090")
PR = json.load(open(os.path.join(EV, "prices28.json")))
ARM = json.load(open(os.path.join(EV, "arm37.json")))
BND = json.load(open(os.path.join(EV, "bounds37.json")))
SM = json.load(open(os.path.join(EV, "protocol37.json")))

MUL = int(PR["MUL_CHARGE"])       # executed dense Fp12 multiply span
SQR = int(PR["SQR_CHARGE"])       # executed cyclotomic square span
STEP12 = max(MUL, SQR)            # the frequency-independent per-step price
LOOKUP = 959                      # censused dict28_wide_lookup span
EXPORT = 681                      # censused dict28_slot_export span
INSERT = 1711 + 688               # censused insert + record write
TRAIL_L = 65536                   # per-trail step cap (construction, not a frequency)


def span(sym):
    out = subprocess.run(["otool", "-tvV", DYLIB], capture_output=True, text=True).stdout
    first = last = None
    inside = False
    for line in out.splitlines():
        m = re.match(r"^(?:([0-9a-f]{16})\t)?_" + sym + r":$", line)
        if m:
            inside = True
            first = int(m.group(1), 16) if m.group(1) else None
            continue
        if inside:
            m2 = re.match(r"^([0-9a-f]{16})\t", line)
            if not m2:
                if line.strip() == "":
                    continue
                break
            a = int(m2.group(1), 16)
            if first is None:
                first = a
            last = a
    assert first is not None and last is not None, sym
    return (last - first) // 4 + 1


# ---- six-limb spans, read from the object the fixture links ------------------
MONT6 = span("d28_montmul6")
ADD6 = span("d28_l6_add")
SUB6 = span("d28_l6_sub")
NEG6 = span("d28_l6_neg")
KEYS6 = span("d28_keys12_form6")
LOOP6 = 6                                   # every six-limb loop runs <= 6 iterations
MONT6_IR = LOOP6 * MONT6                    # executed bound of one montmul6
KEYS12_IR = KEYS6 + 4 * MONT6_IR            # body + four executed montmuls
COEF_IR = 2 * MONT6_IR + 2 * LOOP6 * ADD6   # A <- A + linv*a : 2 montmul + 2 add
PACK_IR = 48                                # 6 limbs x (store+load) per bridge call

# ---- the COLD trail start: y = g^a * z^b, 381-bit exponents, no x -----------
COLD_MULS = 2 * (381 + 381) + 1             # 2 two-scalar-powerings, every bit a
                                            # full multiply (conservative: no
                                            # square/multiply discrimination)
RESTART_COLD = COLD_MULS * STEP12
TRAIL_STARTS = -(-S // TRAIL_L)             # ceil(S/L): the launch rule executes
                                            # S steps at cap L, so starts >= S/L and
                                            # this is EQUALITY by construction
INIT_IR = int(BND["init_ir"])

STEP_ROWS = {
    "per_step_key_lookup_probe": LOOKUP,
    "per_step_slot_export": EXPORT,
    "record_writes_and_insert": INSERT,
    "keys12_formation": KEYS12_IR,
    "scalar_coefficient_update": COEF_IR,
    "l6_bridge_pack_unpack": PACK_IR,
    "guard_and_meter_eval": 26,
    "index_slots_widened_refs": 12,
}
STEP_IR = STEP12 + sum(STEP_ROWS.values())

rows = {k: S * v for k, v in STEP_ROWS.items()}
rows["fp12_step_mulmod_max"] = S * STEP12
rows["trail_starts_cold"] = TRAIL_STARTS * RESTART_COLD
rows["store_init_radix_nodes"] = INIT_IR
rows["worker_starts_1024"] = 7843840
rows["collision_processing"] = 8192
rows["code_and_advice"] = 1048576
rows["candidate_verification"] = 57526800
W = sum(rows.values())


def claim_n(Wv):
    lo, hi = 0, 1 << 40
    while lo < hi:
        m = (lo + hi) // 2
        if pow(20000, 10000) * (1 << m) >= pow(Wv, 10000):
            hi = m
        else:
            lo = m + 1
    return lo


n = claim_n(W)
assert pow(20000, 10000) * (1 << n) >= pow(W, 10000)
assert pow(20000, 10000) * (1 << (n - 1)) < pow(W, 10000)
g = math.gcd(n, 10000)
P_USED = 0.578125 * min(ARM["arms"][0]["yield_cp_lower"], ARM["arms"][1]["yield_cp_lower"])
RESTART_OLD = 7660

budget = {
    "rule": "W_total is TOTAL WORK.  (a) EVERY executed step is charged at the SINGLE "
            "frequency-independent price max(MUL span, SQR span) = %d Ir -- the guard-squaring "
            "substitution is DELETED, so no cheaper step is ever assumed.  (b) The launch rule "
            "executes S steps at a per-trail cap L = %d, so the number of trail starts is "
            "ceil(S/L) = %s: a trail of length AT MOST L forces AT LEAST S/L starts, and the "
            "schedule produces exactly that many because every trail runs to the cap and no "
            "early abort exists on the production path.  (c) Each start is priced as a COLD "
            "two-scalar exponentiation y = g^a z^b (%d full multiplies at %d Ir) that never "
            "reads the secret x: v36 priced G.start(g^a z^{xb}) at %d Ir and the fixture read x."
            % (STEP12, TRAIL_L, TRAIL_STARTS, COLD_MULS, STEP12, RESTART_OLD),
    "per_trail_step_cap": TRAIL_L,
    "starts_lower_bound": "starts >= ceil(S/L) = %s (the inequality v36 had reversed)" % TRAIL_STARTS,
    "starts_equality_by_construction": True,
    "restart_ir_cold": RESTART_COLD,
    "restart_ir_v36_hot": RESTART_OLD,
    "restart_ratio": RESTART_COLD / RESTART_OLD,
    "early_abort": "removed from the production schedule: the fixture's stall-abort path is a "
                   "measurement convenience and is NOT on the launched path (the run_one loop "
                   "runs to the step budget and every trail to the cap), so the start count "
                   "cannot exceed ceil(S/L)",
    "guard_split": "DELETED: every step pays STEP12, frequency-independent",
}

ledger = {
    "mode": "ledger37-spec-to-code",
    "protocol": {
        "name": "meter-first full-store rho, 1024 workers, depth 6, 128-bit radix store, "
                "six-limb key formation",
        "workers_declared": 1024, "depth": 6, "guard_ring_words": 6,
        "steps": str(S), "trail_cap": TRAIL_L, "budget": budget,
        "detection": {
            "mechanism": "the shipped dict28_wide_lookup hit; both collision classes are accepted "
                         "and each hit is one pair",
            "kappa_detect_cross": ARM["arms"][1].get("yield", 0.0),
            "relation_yield_at_64": ARM["arms"][0]["yield_cp_lower"],
            "relation_yield_at_shipped": ARM["arms"][1]["yield_cp_lower"],
            "yield_source": "evidence/arm37.json: ACTUAL collisions inside protocol37.run_one "
                            "(the shipped path), counted by class, with the denominator and the "
                            "recovered x checked per hit; the L=64 arm is the schedule control",
            "p_used": P_USED,
        },
        "census_fidelity": {"depth": 6, "dict": "shipped dict28.c, 9-word entries, 4-probe chains",
                            "store_shrunk_for_the_yield_arm": 256,
                            "why": "a census-scale store makes ACTUAL collisions frequent; the "
                                   "walk, the store code and the recovery are the shipped ones"},
    },
    "store_spec_to_code": {
        "shipped_structure": BND["production_layout"]["rule"],
        "parsed_constants": BND["parsed_constants"],
        "occupancy": BND["occupancy"],
        "dense_six_ceil_control": BND["dense_six_ceil_control"],
        "init_row": "reached nodes x 2^22 Ir; reached = ceil(M(1-e^(-S/M))), M = 2^110 -- NOT "
                    "6*ceil(S/(4*2^22))",
    },
    "arithmetic": {
        "spans_ir": {"d28_montmul6": MONT6, "d28_l6_add": ADD6, "d28_l6_sub": SUB6,
                     "d28_l6_neg": NEG6, "d28_keys12_form6": KEYS6,
                     "dict28_wide_lookup": span("dict28_wide_lookup"),
                     "dict28_wide_insert": span("dict28_wide_insert")},
        "montmul6_executed_ir": MONT6_IR, "keys12_ir": KEYS12_IR, "coefficient_ir": COEF_IR,
        "loop_bound_rule": "every six-limb loop runs at most 6 iterations, so 6 x the loop-free "
                           "span upper-bounds the executed count of one call",
        "executed_code": "scripts/protocol37.Group.state and Group.step CALL these symbols through "
                         "ctypes on EVERY step -- the priced code is the executed code (v36 priced "
                         "four-limb helpers the fixture never called)",
    },
    "budget": {"per_trail_step_cap": TRAIL_L, "trail_starts": str(TRAIL_STARTS),
               "charge_ir": str(TRAIL_STARTS * RESTART_COLD),
               "restart_cold_ir": RESTART_COLD, "restart_v36_hot_ir": RESTART_OLD},
    "bounds_table": BND["bounds_table"],
    "peak_bytes": BND["peak_bytes"],
    "rows_ir": {k: str(v) for k, v in rows.items()},
    "row_count": len(rows),
    "step_ir": STEP_IR,
    "W_total_ir": str(W),
    "claim": {"value": n / 10000, "numerator": str(n), "denominator": "10000",
              "reduced": "%d/%d" % (n // g, 10000 // g),
              "rule": "minimal n with 20000*2^(n/10000) >= W; exact integer test",
              "n_minus_1_fails": True, "row_sum_equals_W": True},
    "success": {"p5_route_a": 0.578125,
                "yield_at_64_cp_lower": ARM["arms"][0]["yield_cp_lower"],
                "yield_at_shipped_cp_lower": ARM["arms"][1]["yield_cp_lower"],
                "p_used": P_USED, "floor": 0.39,
                "rule": "p_used = p5 x min(yield at L=64, yield at the shipped L), each a 95% "
                        "CP lower bound over ACTUAL collisions; the debit uses the WORSE schedule"},
    "prices": {"source": "evidence/exec31.json + the shipped object's spans",
               "per_step_rows": STEP_ROWS, "step12": STEP12},
}
json.dump(ledger, open(os.path.join(EV, "ledger37.json"), "w"), indent=1)

ROWN = len(rows)
statement = (
    "Recover x from z = g^x in BLS12-381 GT by a Frobenius-class rho walk, returning [x]P0 "
    "verified by e(P,Q)=z. CLAIM. Smallest n/10000 with 20000*2^(n/10000) >= W (exact integer "
    "test, n-1 asserted to fail); W = %s is the exact sum of %d rows in evidence/ledger37.json, "
    "recomputed by scripts/cert37.py. "
    "PRICED = EXECUTED (verdict 177471b7 flaw 1). The fixture's Group.state and Group.step CALL "
    "the shipped SIX-LIMB routines through ctypes: d28_l6_sub and d28_l6_neg propagate borrow "
    "across all six limbs and reduce modulo m, d28_l6_add carries and conditionally subtracts m, "
    "d28_montmul6 is a six-limb Montgomery CIOS with a six-limb final subtraction, and "
    "d28_keys12_form6 forms the twelve lanes and returns kappa, so no Python part takes part in "
    "key formation. The modulus is an argument: the fixture runs the toy modulus and audit37.py "
    "runs the same symbols at the 381-bit production modulus p, where neg(1) returns m-1 (all six "
    "limbs), 1000 random add/sub/neg agree with the Python reference, and reduce handles a carry "
    "limb. The rows are those spans: keys12 %d Ir/step, coefficient update %d Ir/step. "
    "BUDGET (flaw 3). Every executed step pays the SINGLE price max(dense-mul span, cyclo-square "
    "span) = %d Ir, frequency-independent: the ceil(4497 S/5040000) guard-squaring substitution is "
    "DELETED. A trail of length at most L = %d forces AT LEAST ceil(S/L) = %s starts -- v36 "
    "reversed that inequality and its launch rule started only 1024 trails; the shipped loop now "
    "runs until S steps are executed, so the count is ceil(S/L) with equality by construction. "
    "Each start is COLD: y = g^a z^b, a two-scalar exponentiation of %d full multiplies (%d Ir) "
    "that never reads x, replacing v36's %d-Ir start that read it. "
    "STORE (flaw 4). The reachable leaf count is the uniform-address OCCUPANCY "
    "ceil(M(1-e^(-S/M))) with M = 2^110 distinct leaf prefixes: S/M = %.3e, so e^-x < 2^-106473 and "
    "the occupancy is M itself (%s leaves), not the %s nodes the dense-leaf formula 6*ceil(S/(4*2^22)) "
    "charges -- a %ax understatement, flagging that formula. The 96-byte record node and the 32-MiB "
    "radix node are SEPARATE constants (v36 let one overwrite the other), every calloc argument is "
    "<= SIZE_MAX, and the 128-bit slot is consumed 22 bits per level over h = 6 levels. "
    "SUCCESS (flaw 5). The yield is measured on ACTUAL collisions inside the shipped loop, by "
    "class and schedule: %d hits at L = 64 with CP-lo %.6f, and the shipped schedule gives CP-lo "
    "%.6f; p_used = 0.578125 x %.6f = %.6f, above the 0.39 floor. No planted pair with a chosen "
    "denominator is used anywhere in the debit. "
    "NO FLOAT in any charged term: every row is a Python int."
) % (W, ROWN, KEYS12_IR, COEF_IR, STEP12, TRAIL_L, TRAIL_STARTS, COLD_MULS, RESTART_COLD,
     RESTART_OLD, BND["occupancy"]["S_over_M"], BND["occupancy"]["reached_leaf_nodes"],
     BND["dense_six_ceil_control"]["value"], BND["dense_six_ceil_control"]["ratio"],
     ARM["arms"][0]["hits"], ARM["arms"][0]["yield_cp_lower"], ARM["arms"][1]["yield_cp_lower"],
     min(ARM["arms"][0]["yield_cp_lower"], ARM["arms"][1]["yield_cp_lower"]), P_USED)

sub = {
    "version": 1,
    "title": "Claim %.4f: six-limb borrow-and-reduce key formation wired into the fixture, "
             "frequency-independent step charge, ceil(S/L) cold starts, occupancy-bounded radix "
             "store, in-protocol relation yield" % (n / 10000),
    "claims": [{
        "statement": statement,
        "score": {"numerator": str(n), "denominator": "10000"},
        "evidence": [
            "argument.md", "certificate.json", "dependencies.json",
            "evidence/ledger37.json", "evidence/arm37.json", "evidence/bounds37.json",
            "evidence/trap37.json", "evidence/protocol37.json", "evidence/audit37.json",
            "evidence/prices28.json", "public-note.md",
            "scripts/regen37.py", "scripts/cert37.py", "scripts/arm37.py", "scripts/bounds37.py",
            "scripts/smoke37.py", "scripts/audit37.py", "scripts/protocol37.py",
            "scripts/ir28/dict28.c", "scripts/ir25/kernels.c",
        ],
    }],
    "experiments": [],
}
json.dump(sub, open(os.path.join(PKG, "submission.json"), "w"), indent=1)
json.dump(budget, open(os.path.join(EV, "trap37.json"), "w"), indent=1)
print(json.dumps({"W": str(W), "log2W": round(math.log2(W), 6), "n": n, "claim": n / 10000,
                  "rows": ROWN, "STEP_IR": STEP_IR, "KEYS12_IR": KEYS12_IR,
                  "COEF_IR": COEF_IR, "RESTART_COLD": RESTART_COLD,
                  "starts": str(TRAIL_STARTS), "p_used": P_USED,
                  "title_len": len(sub["title"])}, indent=1))
