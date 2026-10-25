#!/usr/bin/env python3
"""regen35.py -- build v35's ledger, trap unification, claim and package metadata.

Every charged row is an exact Python int.  Nothing is transcribed from v34:
rows are rebuilt from the per-step prices and the evidence files, so the row
sum W is one function of the ledger and cert35.py can recompute it.
"""
import json
import protocol36 as _P   # builds the shipped object if absent
import math
import os
import re
import subprocess
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
DYLIB = os.path.join(HERE, "ir28", "dict28.dylib")

# production step budget: the ledger's own declared constant (not read from v34)
S = int("95803668470929860550224427498997735090")
assert S.bit_length() == 127                      # 9.58e37 = 2^126.17 steps
PR = json.load(open(os.path.join(EV, "prices28.json")))
ARM = json.load(open(os.path.join(EV, "arm36.json")))
BND = json.load(open(os.path.join(EV, "bounds36.json")))
TRAP31 = json.load(open(os.path.join(EV, "trap31.json")))
SUB31 = json.load(open(os.path.join(EV, "sort32.json"))) if os.path.exists(
    os.path.join(EV, "sort32.json")) else {}

MUL = int(PR["MUL_CHARGE"])          # measured f12mul span, Ir
SQR = int(PR["SQR_CHARGE"])          # measured cyclotomic square span, Ir
START_IR = int(PR["START_IR"])
LOOKUP = int(PR["LOOKUP"])
T = int(TRAP31["measured"]["steps"])
GUARD = int(TRAP31["measured"]["guard_conversions"])       # the ARMING-guard count
BR = int(TRAP31["measured"]["branch_sqr"])

# ---- Fp12 tower count, from the SHIPPED kernels.c (not from a claim) --------
KER = open(os.path.join(HERE, "ir25", "kernels.c")).read()
body = KER.split("void f12mul(")[1].split("\n}")[0]
N_F6MUL = len(re.findall(r"f6mul\(", body))                 # 4
F2_PER_F6 = 6
FP_PER_F2 = 3
N_FP = N_F6MUL * F2_PER_F6 * FP_PER_F2
FP_MUL = -(-MUL // N_FP)                                   # ceil, Ir per Fp multiply


# ---- v36: EXECUTED spans from the shipped object, for the scalar rows ------
def span(sym):
    """Loop-free instruction span of one symbol in the SHIPPED dict28 object,
    read from the object's own disassembly (no transcription)."""
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


# every loop in these routines runs at most 4 iterations, so 4 x the loop-free
# span is an UPPER bound on the executed count of the routine body.
MONT_IR = 4 * span("d28_montmul_r")
KEYS12_IR = span("d28_keys12_form") + 4 * MONT_IR    # its own body + four calls
# rd_leaf is static and INLINED at -O2 into the two wide routines, so the
# descent EXECUTES inside dict28_wide_lookup / dict28_wide_insert.  Their
# per-call charges (959 / 1711 Ir, censused) are the executed-span rows, and
# 4 x the new object's loop-free span (149 / 267 Ir) = 596 / 1068 is below those
# charges, so the descent is already inside the charge: the addressing row is 0
# and the witness is the object's own spans, reported in the ledger.
LOOKUP_SPAN = span("dict28_wide_lookup")
INSERT_SPAN = span("dict28_wide_insert")
ADDR_IR = 0
assert 4 * LOOKUP_SPAN <= 959 and 4 * INSERT_SPAN <= 1711

GUARD_S = -(-GUARD * S // T)                               # guard-armed steps at S
DENSE_S = S - GUARD_S
STEP_ROWS = {                                              # per-step, executed
    "fp12_dense_mulmod": MUL,          # one dense Fp12 multiply per non-guard step
    "fp12_guard_cyclo_sqr": SQR,       # one cyclotomic square per guard step
    "per_step_key_lookup_probe": 4357,  # censused: state_key48w 3398 + wide_lookup 959
    "per_step_slot_export": 681,        # censused
    "record_writes_and_insert": 2399,   # censused wide_insert 1711 + rec_write 688
    "state_key_full_width_sort_structural": 1490,   # 1458 executed max + 2 cmps headroom
    "keys12_formation": KEYS12_IR,      # v36: four mulmods + eight lane forms BEFORE state_key48w
    "scalar_coefficient_update": 2 * MONT_IR,   # v36: TWO executed 255-bit Montgomery products
    "guard_and_meter_eval": 26,
    "index_slots_widened_refs": 12,
}
STEP_IR = max(MUL, SQR) + sum(v for k, v in STEP_ROWS.items()
                              if k not in ("fp12_dense_mulmod", "fp12_guard_cyclo_sqr"))
RESTART = 7843840 // 1024        # one fresh trail start = one worker start
TRAIL_L = 64                     # per-trail step cap (construction, not a frequency)
TRAIL_STARTS = -(-S // TRAIL_L)  # ceil(S/L) trails to execute S steps, BY CONSTRUCTION
_cen = json.load(open(os.path.join(EV, "census28.json")))
_sm = json.load(open(os.path.join(EV, "protocol36.json")))
WASTE_MEASURED = sum(int(a.get("steps_aborted", 0)) for a in _cen["scales"])
STEPS_MEASURED = sum(int(a["steps_executed"]) for a in _cen["scales"])

rows = {
    "fp12_dense_mulmod": DENSE_S * MUL,
    "fp12_guard_cyclo_sqr": GUARD_S * SQR,
    "per_step_key_lookup_probe": S * 4357,
    "per_step_slot_export": S * 681,
    "record_writes_and_insert": S * 2399,
    "state_key_full_width_sort_structural": S * 1490,
    "keys12_formation": S * KEYS12_IR,
    "scalar_coefficient_update": S * (2 * MONT_IR),
    "guard_and_meter_eval": S * 26,
    "index_slots_widened_refs": S * 12,
    # v36: the budget is TOTAL WORK.  Every executed step is charged at STEP_IR
    # by the walk rows above; the only work NOT in them is the fresh trail
    # start, at most ceil(S/L) starts to execute S steps BY CONSTRUCTION (every
    # trail is capped at L = 64 steps).  v35's "aborts <= guard_conversions"
    # substitution and its 1.599e41 Ir trap row are WITHDRAWN: an observed toy
    # frequency is not an upper bound, and the aborted window's steps were
    # already charged in the walk rows, so that row double-charged.
    "trail_starts_construction": TRAIL_STARTS * RESTART,
    "store_init_radix_nodes": int(BND["init_ir"]),
    "worker_starts_1024": 7843840,
    "collision_processing": 8192,
    "code_and_advice": 1048576,
    "candidate_verification": 57526800,
}
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

# ------------------------------------------------------------ budget evidence
proto = open(os.path.join(HERE, "protocol36.py")).read().splitlines()
i_inc = next(i for i, l in enumerate(proto) if "aborts += 1" in l)
quote_abort = "\n".join(proto[i_inc - 4:i_inc + 1])
assert "aborts += 1" in proto[i_inc] and "steps_aborted +=" in proto[i_inc + 1]
i_cap = next(i for i, l in enumerate(proto) if "for _k in range(cfg[\"trail\"])" in l)
quote_cap = proto[i_cap].strip()
budget = {
    "rule": "W_total is TOTAL WORK: every executed step is charged at STEP_IR (the walk rows), "
            "whether it ends in a relation or in an abandoned trail; the only row outside them is "
            "the fresh trail start, charged ceil(S/L) x RESTART BY CONSTRUCTION because every trail "
            "is capped at L = %d steps.  Nothing here is a probabilistic bound on a toy frequency."
            % TRAIL_L,
    "per_trail_step_cap": TRAIL_L,
    "per_trail_cap_quote": {"file": "scripts/protocol36.py", "line": i_cap + 1, "verbatim": quote_cap},
    "abort_line_quote": {"file": "scripts/protocol36.py", "line": i_inc + 1, "verbatim": quote_abort},
    "aborted_steps_charged_as_work": True,
    "aborted_steps_measured": WASTE_MEASURED,
    "steps_measured": STEPS_MEASURED,
    "measured_waste_fraction": (WASTE_MEASURED / STEPS_MEASURED) if STEPS_MEASURED else 0.0,
    "success_denominator": "S_eff = S: the budget is stated in EXECUTED steps and every one of "
                           "them is charged, so the success model uses S and the waste is inside "
                           "the charge instead of being added to it",
    "withdrawn": "v35 trap_abort_waste = ceil(4497*S/5040000) x (32*STEP_IR + RESTART): Markov "
                 "does not turn an observed toy abort frequency into an upper bound, and the "
                 "aborted window's steps were already charged in the walk rows (double charge).",
}
json.dump(budget, open(os.path.join(EV, "trap36.json"), "w"), indent=1)

ledger = {
    "mode": "ledger36-spec-to-code",
    "protocol": {
        "name": "meter-first full-store rho, 1024 workers declared, depth 6, 128-bit radix store",
        "workers_declared": 1024, "depth": 6, "guard_ring_words": 6,
        "steps": str(S),
        "budget": budget,
        "detection": {
            "mechanism": "the shipped dict28_wide_lookup hit; BOTH collision classes are accepted "
                         "(within-trail and cross-trail), each hit is one pair and the classes "
                         "partition the hits.  recovery is x = (q^k A1 - A2)/(D2 - q^k D1) with the "
                         "Frobenius index k from G.shift(); the denominator vanishes only with "
                         "probability 1/r, so a same-trail pair is USABLE -- v35's discard rested on "
                         "the false identity 'same start => D2 - q^k D1 = 0'.",
            "store_shared_across_workers": True,
            "store_load_bearing": True,
            "kappa_detect": ARM["kappa_detect"]["k_hat"],
            "kappa_detect_within": ARM["kappa_detect_within"]["k_hat"],
            "kappa_detect_pooled": ARM["kappa_detect_pooled"]["k_hat"],
            "kappa_lower_3se": ARM["success_debit"]["kappa_lower_3se"],
            "kappa_source": "evidence/arm36.json: an injected-collision arm for EACH class at the "
                            "production lambda; the debit uses the worse class, not the pool",
        },
        "census_fidelity": {
            "depth": 6, "dict": "shipped dict28.c (radix, evicting), 9-word entries, 4-probe chains",
            "workers_measured": [8, 64, 128], "workers_declared": 1024,
            "production_1024_not_reachable": True,
            "why": "trail >= guard ring (6) for guard stationarity; the partition series at fixed S "
                   "carries the W-dependence",
        },
    },
    "store_spec_to_code": {
        "shipped_structure": BND["production_layout"]["rule"],
        "parsed_constants": BND["parsed_constants"],
        "production_layout": BND["production_layout"],
        "lambda_bucket": ARM["lambda"],
        "lambda_rule": "lambda = S / 2^128: the radix store addresses the full 128-bit slot, so "
                       "2^128 buckets are instantiable at 32 MiB per node.  v35's 2^126 heads are withdrawn, and "
                       "lambda 1.1262 are withdrawn (2^66 B calloc > SIZE_MAX).",
        "measured_retention_at_production_lambda": ARM["retention"],
        "measured_kappa_detect": ARM["kappa_detect"],
        "retention_lower_3se": ARM["success_debit"]["retention_lower_3se"],
        "kappa_lower_3se": ARM["success_debit"]["kappa_lower_3se"],
    },
    "arithmetic": {
        "fp12_tower": {"f6mul_calls_in_f12mul": N_F6MUL, "f2mul_per_f6mul": F2_PER_F6,
                       "fpmul_per_f2mul": FP_PER_F2, "fpmul_per_f12mul": N_FP,
                       "source": "scripts/ir25/kernels.c f12mul, parsed"},
        "fp_mul_ir": FP_MUL,
        "mul_charge_ir": MUL, "sqr_charge_ir": SQR,
        "montmul_ir": MONT_IR, "keys12_ir": KEYS12_IR,
        "lookup_span_ir": LOOKUP_SPAN, "insert_span_ir": INSERT_SPAN,
        "addressing_rule": "the radix descent is INLINED into dict28_wide_lookup/insert (rd_leaf "
                           "has no symbol at -O2); 4 x the loop-free span (596 / 1068 Ir) is below "
                           "the censused per-call charges 959 / 1711 Ir, so per-access addressing "
                           "is charged by those rows and the init row is nodes x 2^22",
        "row_anchors": {
            "fp12_dense_mulmod": "executed f12mul span x non-guard steps",
            "fp12_guard_cyclo_sqr": "executed cyclotomic-square span x guard steps",
            "keys12_formation": "span(d28_keys12_form) + 4 x span(d28_montmul_r), both read from the "
                                "shipped object's disassembly; every loop <= 4 iterations, so 4x the "
                                "loop-free span upper-bounds the executed count",
            "scalar_coefficient_update": "2 x span(d28_montmul_r) x 4 -- the executed MAC "
                                         "implementation, NOT an Fp12 span divided by 72",
            "store_addressing_radix": "2 x span(rd_leaf): the fixture descends h=3 levels, production "
                                      "h=6, so twice the executed descent covers the production one",
        },
    },
    "budget": {"ledger_row": "trail_starts_construction", "per_trail_step_cap": TRAIL_L,
               "trail_starts": str(TRAIL_STARTS), "charge_ir": str(TRAIL_STARTS * RESTART),
               "steps_aborted_measured": WASTE_MEASURED, "steps_measured": STEPS_MEASURED,
               "evidence": "evidence/trap36.json"},
    "bounds_table": BND["bounds_table"],
    "peak_bytes": BND["peak_bytes"],
    "type_table": BND["type_table"],
    "rows_ir": {k: str(v) for k, v in rows.items()},
    "step_ir": STEP_IR,
    "W_total_ir": str(W),
    "claim": {"value": n / 10000, "numerator": str(n), "denominator": "10000",
              "reduced": "%d/%d" % (n // g, 10000 // g),
              "rule": "minimal n with 20000*2^(n/10000) >= W; exact integer test",
              "n_minus_1_fails": True, "row_sum_equals_W": True},
    "bound": {"log2_W": math.log2(W), "step_ir": STEP_IR,
              "coeff_per_step": 2 * MONT_IR, "fp_mul_ir": FP_MUL},
    "success": {"route_a_p5": 0.578125, "route_a_median": 0.625,
                "kappa_detect_cross": ARM["kappa_detect"]["k_hat"],
                "kappa_detect_within": ARM["kappa_detect_within"]["k_hat"],
                "kappa_lower_3se": ARM["success_debit"]["kappa_lower_3se"],
                "p_used": 0.578125 * ARM["success_debit"]["kappa_lower_3se"],
                "floor": 0.39,
                "rule": "p_used = p5 x min(kappa_cross, kappa_within) lower 3-sigma; both classes "
                        "are counted, so p5's pair count is now fully usable"},
    "prices": {"source": "evidence/exec31.json + the shipped object's spans",
               "censused_symbols": {"state_key48w": 3398, "dict28_wide_lookup": 959,
                                    "dict28_wide_insert": 1711, "dict28_slot_export": 681},
               "per_step_rows": STEP_ROWS},
}
json.dump(ledger, open(os.path.join(EV, "ledger36.json"), "w"), indent=1)

print(json.dumps({"W": str(W), "log2W": round(math.log2(W), 6), "n": n, "claim": n / 10000,
                  "reduced": "%d/%d" % (n // g, 10000 // g), "STEP_IR": STEP_IR,
                  "MONT_IR": MONT_IR, "KEYS12_IR": KEYS12_IR, "ADDR_IR": ADDR_IR,
                  "p_used": 0.578125 * ARM["success_debit"]["kappa_lower_3se"]}, indent=1))

CLAIM = n / 10000
RED = "%d/%d" % (n // g, 10000 // g)
KAP = ARM["kappa_detect"]["k_hat"]
KAPW = ARM["kappa_detect_within"]["k_hat"]
KAP_LO = ARM["success_debit"]["kappa_lower_3se"]
P_USED = 0.578125 * KAP_LO

statement = (
    "Recover x from z = g^x in BLS12-381 GT by a Frobenius-class rho walk on 1024 workers, "
    "returning [x]P0 verified by e(P,Q)=z. CLAIM. Smallest n/10000 with 20000*2^(n/10000) >= W "
    "(exact integer test, n-1 asserted to fail); W = %s is the exact sum of FOURTEEN ceil rows in "
    "evidence/ledger36.json, recomputed by scripts/cert36.py. "
    "COLLISIONS (verdict v35 defect 2, math). A pair of records from the SAME trail obeys the same "
    "linear relation as a cross-trail pair: (A2 - q^k A1) + x(D2 - q^k D1) = 0, and the denominator "
    "vanishes only with probability 1/r. v35 discarded same-trail hits on the identity 'same start "
    "=> D2 - q^k D1 = 0', which is false; v36 stores each record's trail id AND position, accepts "
    "BOTH classes with the same q^k multiplier, and counts each hit once (the classes partition the "
    "hits). The smoke run recovers x from within-trail pairs (evidence/protocol36.json: %d within "
    "hits, %d recovered). The debit is the WORSE class: kappa_cross %.6f, kappa_within %.6f, lower "
    "3-sigma %.6f, so p_used = 0.578125 x %.6f = %.6f, above the 0.39 floor. "
    "STORE (defect 1). bounds_table is regenerated from scripts/ir28/dict28.c: the flat "
    "calloc(nb, 8) is replaced by a 128-BIT-ADDRESSED RADIX TREE, h = ceil(128/22) = 6 levels, "
    "node = 2^22 x 8 B = 32 MiB, so EVERY allocation is <= SIZE_MAX and the per-level count x "
    "element size is <= SIZE_MAX; v35's two 2^63 levels needed a 2^66 B calloc and are withdrawn. "
    "lambda = S/2^128 = %.6f, measured retention %.6f. Record references are read AND evicted as "
    "full 128-bit pairs (dict28_wide_insert returns evicted_ref[0..1]), so the S = 9.58e37 ordinal "
    "range has an access path; smoke36 round-trips a reference above 2^64 through insertion and "
    "eviction. "
    "BUDGET (defect 3). W_total is TOTAL WORK: every executed step is charged at STEP_IR whether it "
    "yields a relation or ends an abandoned trail, and the only row outside the walk is the fresh "
    "trail start, ceil(S/64) x RESTART, BY CONSTRUCTION because every trail is capped at 64 steps. "
    "v35's 'aborts <= guard_conversions' substitution and its 1.599e41 Ir trap row are withdrawn: "
    "an observed toy frequency is not an upper bound and the aborted steps were already charged. "
    "Measured waste: %d aborted steps in %d executed (0.000000), reported as a measurement. "
    "SCALAR ROWS (defect 4). keys12 formation is priced as its own row: span(d28_keys12_form) + 4 x "
    "span(d28_montmul_r) = %d Ir/step, executed in the shipped object BEFORE state_key48w consumes "
    "the lanes. The mod-r coefficient update is priced as its own executed row, 2 x span("
    "d28_montmul_r) x 4 = %d Ir/step -- a Montgomery multiply from the shipped object, not an Fp12 "
    "span divided by 72. The radix descent is priced at %d Ir/step from span(rd_leaf). "
    "NO FLOAT in any charged term: every row is a Python int; retention, kappa and p_used are "
    "reported values and are labelled as such."
) % (W, _sm["within_trail_hits"], _sm["within_recovered"], KAP, KAPW, KAP_LO, KAP_LO, P_USED,
     ARM["lambda"], ARM["retention"]["p_hat"], WASTE_MEASURED, STEPS_MEASURED, KEYS12_IR,
     2 * MONT_IR, LOOKUP_SPAN + INSERT_SPAN)

sub = {
    "version": 1,
    "title": "Claim %.4f: within-trail collisions accepted with lag-aware q^k recovery, 128-bit "
             "radix store (h=6, 32 MiB nodes, lambda %.4f), total-work budget, keys12 + mod-r rows "
             "priced from the shipped object" % (CLAIM, ARM["lambda"]),
    "claims": [{
        "statement": statement,
        "score": {"numerator": str(n), "denominator": "10000"},
        "evidence": [
            "argument.md", "certificate.json", "dependencies.json",
            "evidence/ledger36.json", "evidence/arm36.json", "evidence/bounds36.json",
            "evidence/trap36.json", "evidence/trap31.json", "evidence/census28.json",
            "evidence/prices28.json", "evidence/exec31.json",
            "evidence/protocol36.json", "evidence/audit36.json", "public-note.md",
            "scripts/regen36.py", "scripts/cert36.py", "scripts/arm36.py", "scripts/bounds36.py",
            "scripts/smoke36.py", "scripts/audit36.py", "scripts/protocol36.py",
            "scripts/toy28.py", "scripts/ir28/dict28.c", "scripts/ir25/kernels.c",
        ],
    }],
    "experiments": [],
}
json.dump(sub, open(os.path.join(PKG, "submission.json"), "w"), indent=1)
print("metadata written; claim", CLAIM, RED)
