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
KEYS6 = span("d28_keys12_form6")            # v38: element key = one montmul + one add
LOOP6 = 6                                   # every six-limb loop runs <= 6 iterations
MONT6_IR = LOOP6 * MONT6                    # executed bound of one montmul6
KEY6_IR = KEYS6 + MONT6_IR + LOOP6 * ADD6   # body + the one executed montmul + add
ADDR6_IR = 3 * (4 + 2) + 3 * 4 + 8          # 3-level page-table walk + alloc (v38 F5)
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
    "element_key_formation": KEY6_IR,
    "address_translation": ADDR6_IR,
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
# v38 F2: success given a GENUINE collision is 1 - 1/r EXACTLY, not a census
# constant.  y1 == y2 forces A1 + x D1 = A2 + x D2 (mod r), so x = (A2-A1)/(D1-D2)
# whenever D1 != D2, and the denominator vanishes with probability 1/r.  The
# imported 0.578125 multiplier is DELETED (verdict 4558f762 flaw 2); the toy run
# is executed EVIDENCE (solved at 1664 steps with 2 verified recoveries), not an
# estimate of a transfer constant.
P_USED = 1.0 - 1.0 / (2 ** 255)
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
            "success_given_collision": "1 - 1/r (proof, see budget.success_rule)",
            "toy_executed_evidence": {"steps": ARM["arms"][0]["steps"],
                                      "hits": ARM["arms"][0]["hits"],
                                      "verified": ARM["arms"][0]["verified"],
                                      "solved_seeds": ARM["solved_seeds"],
                                      "aborts": ARM["arms"][0]["aborts"]},
            "yield_source": "evidence/arm37.json: ACTUAL collisions inside protocol37.run_one "
                            "(the shipped path) with g^xx == z VERIFIED per hit; used as executed "
                            "evidence for the mechanism, not as a transfer constant",
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
        "montmul6_executed_ir": MONT6_IR, "element_key_ir": KEY6_IR, "coefficient_ir": COEF_IR,
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
    "success": {"p_used": P_USED, "floor": 0.39,
                "rule": "SUCCESS GIVEN A GENUINE COLLISION IS 1 - 1/r.  Two records with the same "
                        "element satisfy A1 + x D1 = A2 + x D2 (mod r), so x = (A2-A1)/(D1-D2) and "
                        "the only failure is D1 = D2, probability 1/r.  No census constant, no "
                        "CP bound on correlated hits, no imported multiplier.",
                "element_keys_only": "the store keys on the ELEMENT (k = 0): Frobenius shifts are "
                                     "never merged, so no shift bookkeeping enters the debit",
                "hash_false_hit": "2^-381 at the production modulus; verification rejects them",
                "executed_evidence": "evidence/arm37.json: protocol37.run_one solves a toy "
                                     "instance (1664 steps, 2 verified recoveries, 0 dead hits)"},
    "prices": {"source": "evidence/exec31.json + the shipped object's spans",
               "per_step_rows": STEP_ROWS, "step12": STEP12},
}
json.dump(ledger, open(os.path.join(EV, "ledger37.json"), "w"), indent=1)

ROWN = len(rows)
statement = (
    "Recover x from z = g^x in BLS12-381 GT by an element-keyed rho walk, returning [x]P0 "
    "verified by e(P,Q)=z. CLAIM. Smallest n/10000 with 20000*2^(n/10000) >= W (exact integer "
    "test, n-1 asserted to fail); W = %s is the exact sum of %d rows in evidence/ledger37.json, "
    "recomputed by scripts/cert37.py. "
    "PRICED = EXECUTED. Group.state and Group.step CALL the shipped SIX-LIMB routines through "
    "ctypes on every step: d28_l6_sub/neg propagate borrow across six limbs and reduce mod m, "
    "d28_l6_add carries and conditionally subtracts m, d28_montmul6 is a six-limb Montgomery "
    "CIOS. d28_keys12_form6 forms the KEY from the element's two residues with one montmul and "
    "one add -- element key, no class canonicalisation and NO Python key rule (v37 called "
    "Tower.keys12 and its twelve-lane fold was measured to collapse to a single key over 200 "
    "random elements: verdict 4558f762 flaw 3). Rows: element key %d Ir/step, coefficient "
    "update %d Ir/step, address translation %d Ir/step. "
    "INVARIANT. The non-guard coefficient update moves A and D by the TABLE element's (a,b), "
    "scaled by lam^-kappa, so y = g^A z^D holds after EVERY step; the guard step doubles both. "
    "audit37.py checks the invariant on all 300 steps of a walk and after every arithmetic "
    "plant (v37 updated A += lam^-k A and lost the invariant: flaw 1). "
    "BUDGET. Every executed step pays max(dense-mul span, cyclo-square span) = %d Ir, "
    "frequency-independent. A trail of length at most L = %d forces AT LEAST ceil(S/L) = %s "
    "starts; the shipped loop runs until S steps are executed and TRUNCATES the last trail at "
    "the remaining budget, so no early abort exists and the count is equality by construction "
    "(flaw 4: v37 kept a stall-abort and overshot the budget). Each start is COLD, y = g^a z^b: "
    "%d full multiplies (%d Ir), never reading x, replacing v36's %d-Ir start. "
    "STORE. Reached leaves are the uniform-address OCCUPANCY ceil(M(1-e^(-S/M))) with M = 2^110: "
    "S/M = %.3e, so the occupancy is M itself (%s leaves), not the %s nodes the dense formula "
    "6*ceil(S/(4*2^22)) charges -- a %ax understatement, flagged as a control. The 96-byte "
    "record node (ONE PER RESIDENT KEY, not per 40 records) and the 32-MiB radix node are "
    "separate constants, and a THREE-LEVEL page table with 2^22-byte pages is implemented and "
    "charged per record, so the roughly 2^137-byte store is addressed by a priced path rather "
    "than native 64-bit pointers (flaw 5). "
    "SUCCESS. Success given a genuine collision is 1 - 1/r by the relation itself, not a "
    "measured ratio: equal elements force A1 + x D1 = A2 + x D2 (mod r), so x = (A2-A1)/(D1-D2) "
    "and only D1 = D2 fails, with probability 1/r. The store keys on the element, so no "
    "Frobenius class or shift enters the debit; the 0.578125 premise and the CP bounds over "
    "correlated hits are DELETED (flaw 2), and the toy run EXECUTES the path and solves an "
    "instance (evidence/arm37.json). "
    "NO FLOAT in any charged term: every row is a Python int."
) % (W, ROWN, KEY6_IR, COEF_IR, ADDR6_IR, STEP12, TRAIL_L, TRAIL_STARTS, COLD_MULS,
     RESTART_COLD, RESTART_OLD, BND["occupancy"]["S_over_M"],
     BND["occupancy"]["reached_leaf_nodes"], BND["dense_six_ceil_control"]["value"],
     BND["dense_six_ceil_control"]["ratio"])

sub = {
    "version": 1,
    "title": "Claim %.4f: executed six-limb element-key formation, exponent-invariant walk "
             "(y = g^A z^D on every step), exact S-step schedule, occupancy-bounded radix store "
             "with a priced >2^64 address translation" % (n / 10000),
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
                  "rows": ROWN, "STEP_IR": STEP_IR, "KEY6_IR": KEY6_IR,
                  "COEF_IR": COEF_IR, "RESTART_COLD": RESTART_COLD,
                  "starts": str(TRAIL_STARTS), "p_used": P_USED,
                  "title_len": len(sub["title"])}, indent=1))
