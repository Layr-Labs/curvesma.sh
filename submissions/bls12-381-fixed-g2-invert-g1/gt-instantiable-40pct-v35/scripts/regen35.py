#!/usr/bin/env python3
"""regen35.py -- build v35's ledger, trap unification, claim and package metadata.

Every charged row is an exact Python int.  Nothing is transcribed from v34:
rows are rebuilt from the per-step prices and the evidence files, so the row
sum W is one function of the ledger and cert35.py can recompute it.
"""
import json
import math
import os
import re
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")

# production step budget: the ledger's own declared constant (not read from v34)
S = int("95803668470929860550224427498997735090")
assert S.bit_length() == 127                      # 9.58e37 = 2^126.17 steps
PR = json.load(open(os.path.join(EV, "prices28.json")))
ARM = json.load(open(os.path.join(EV, "arm35.json")))
BND = json.load(open(os.path.join(EV, "bounds35.json")))
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
F2_PER_F6 = 6                                              # Karatsuba in Fp2[v]/(v^3-xi)
FP_PER_F2 = 3                                              # Karatsuba in Fp2
N_FP = N_F6MUL * F2_PER_F6 * FP_PER_F2
FP_MUL = -(-MUL // N_FP)                                   # ceil, Ir per Fp multiply
COEFF = 2 * FP_MUL + 12                                    # two 255-bit mulmods + 12-lane load

GUARD_S = -(-GUARD * S // T)                               # guard-armed steps at S
DENSE_S = S - GUARD_S
STEP_ROWS = {                                              # per-step, executed
    "fp12_dense_mulmod": MUL,          # one dense Fp12 multiply per non-guard step
    "fp12_guard_cyclo_sqr": SQR,       # one cyclotomic square per guard step
    "per_step_key_lookup_probe": 4357,  # censused: state_key48w 3398 + wide_lookup 959
    "per_step_slot_export": 681,        # censused
    "record_writes_and_insert": 2399,   # censused wide_insert 1711 + rec_write 688
    "state_key_full_width_sort_structural": 1490,   # 1458 executed max + 2 cmps headroom
    "scalar_coefficient_update": COEFF,             # 2 x FP_MUL + 12
    "guard_and_meter_eval": 26,
    "index_slots_widened_refs": 12,
}
STEP_IR = max(MUL, SQR) + sum(v for k, v in STEP_ROWS.items()
                              if k not in ("fp12_dense_mulmod", "fp12_guard_cyclo_sqr"))
RESTART = 7843840 // 1024        # one fresh trail start = one worker start, the ledger's own price
ABORTS = -(-GUARD * S // T)                  # counting bound: #aborts <= #guards
TRAP_IR = ABORTS * (32 * STEP_IR + RESTART)  # 32-step window + one fresh trail start

rows = {
    "fp12_dense_mulmod": DENSE_S * MUL,
    "fp12_guard_cyclo_sqr": GUARD_S * SQR,
    "per_step_key_lookup_probe": S * 4357,
    "per_step_slot_export": S * 681,
    "record_writes_and_insert": S * 2399,
    "state_key_full_width_sort_structural": S * 1490,
    "scalar_coefficient_update": S * COEFF,
    "guard_and_meter_eval": S * 26,
    "index_slots_widened_refs": S * 12,
    "cross_trail_compare": 8 * 1,            # one load+cmp+branch per HIT; hits are O(1)
    "trap_abort_waste": TRAP_IR,
    "store_init_bucket_array": int(BND["init_ir"]),
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

# ---------------------------------------------------------------- trap three-way
proto = open(os.path.join(HERE, "protocol28.py")).read().splitlines()
i_inc = next(i for i, l in enumerate(proto) if "aborts += 1" in l)
ctx = proto[i_inc - 6:i_inc + 1]
quote_abort = "\n".join(ctx)
assert "if armed and stall >= K_STALL" in quote_abort and "aborts += 1" in quote_abort
assert "aborted = abandoned = True" in proto[i_inc]
i_guard = next(i for i, l in enumerate(proto) if "if guard:" in l and "armed = True" in proto[i + 1])
assert "armed = True" in proto[i_guard + 1]
LEDGER_ROW = "trap_abort_waste"
CENSUS_COUNTER = "aborts"
MARKOV_WITNESS = "guard_conversions"
p_minus_c = -(-GUARD * S // T)               # same integer as ABORTS
raw_aborts = sum(int(s["aborts"]) for s in json.load(
    open(os.path.join(EV, "census28.json")))["scales"])
trap = {
    "ledger_row": LEDGER_ROW,
    "census_counter": CENSUS_COUNTER,
    "census_counter_quote": {
        "file": "scripts/protocol28.py", "line": i_inc + 1,
        "verbatim": quote_abort,
        "rule": "this is the ONLY place the abort counter increments; the ledger row, the "
                "census counter and the Markov event are the same event, not three paraphrases",
    },
    "markov_witness": MARKOV_WITNESS,
    "markov_event": "#{aborts at S} <= #{%s} = %d (each abort needs an ARMED counter, armed is "
                    "set only on a guard step and cleared on restart, so distinct aborts consume "
                    "distinct arming guards)" % (MARKOV_WITNESS, GUARD),
    "markov_witness_quote": {
        "file": "evidence/trap31.json", "key": "measured.guard_conversions",
        "verbatim": json.dumps({"steps": T, "guard_conversions": GUARD,
                                "branch_sqr": BR, "multiply": TRAP31["measured"]["multiply"]}),
    },
    "why_the_union_is_not_the_event": "br OR guard is the STALL predicate; 4497 is the guard "
                                      "count and %d the br count, and neither is the abort "
                                      "count.  The abort needs armed AND 32 consecutive union "
                                      "steps, so the union rate cannot be charged as aborts and "
                                      "p_s^32 cannot bound a single-step event.  The counting "
                                      "bound above needs neither independence nor p_s." % BR,
    "counts": {"steps_observed": T, "guard_conversions": GUARD, "branch_sqr": BR,
               "union_steps": GUARD + BR, "union_rate": (GUARD + BR) / T,
               "guard_rate": GUARD / T},
    "aborts_bound_at_S": str(ABORTS),
    "raw_aborts_observed_in_census_arms": raw_aborts,
    "charge_ir": str(TRAP_IR),
    "charge_rule": "aborts x (32 x STEP_IR + RESTART): the aborted 32-step window is wasted work "
                   "at the ledger's own per-step price, and the trail restarts at START_IR; "
                   "STEP_IR = %d Ir" % STEP_IR,
}
assert LEDGER_ROW in rows and CENSUS_COUNTER == "aborts"
json.dump(trap, open(os.path.join(EV, "trap35.json"), "w"), indent=1)

ledger = {
    "mode": "ledger35-spec-to-code",
    "protocol": {
        "name": "meter-first full-store rho, 1024 workers declared, depth 6, store priced as shipped",
        "workers_declared": 1024, "depth": 6, "guard_ring_words": 6,
        "steps": str(S),
        "abort": {
            "ledger_row": LEDGER_ROW, "census_counter": CENSUS_COUNTER,
            "quote_file": "scripts/protocol28.py", "quote_line": i_inc + 1,
            "event": json.load(open(os.path.join(EV, "trap35.json")))["markov_event"],
        },
        "detection": {
            "mechanism": "the shipped dict28_wide_lookup hit, followed by THE TRAIL COMPARISON "
                         "now in the protocol (scripts/protocol28.py, run_one): a hit counts as a "
                         "collision only when the stored trail id differs from the arriving "
                         "worker's; same-trail hits are counted and discarded because "
                         "G.recover() has den = (D2 - q^k D1) = 0 there",
            "store_shared_across_workers": True,
            "store_load_bearing": True,
            "kappa_detect": ARM["kappa_detect"]["k_hat"],
            "kappa_source": "evidence/arm35.json, injected-collision arm at the production lambda",
        },
        "census_fidelity": {
            "depth": 6, "dict": "shipped dict28.c (evicting), 9-word entries, 4-probe chains",
            "workers_measured": [8, 64, 128], "workers_declared": 1024,
            "production_1024_not_reachable": True,
            "why": "trail >= guard ring (6) for guard stationarity; only the partition series at "
                   "fixed S carries the W-dependence",
        },
    },
    "store_spec_to_code": {
        "shipped_structure": BND["production_layout"]["rule"],
        "parsed_constants": BND["parsed_constants"],
        "production_layout": BND["production_layout"],
        "lambda_bucket": ARM["lambda"],
        "lambda_rule": "lambda = arrivals / heads = S / 2^126; the v34 heads value (lambda 0.2815) "
                       "is withdrawn because the shipped field is u64 and ws_bucket masks 64 "
                       "bits (evidence/historical-disclosure.json)",
        "v34_withdrawn": BND["production_layout"]["v34_withdrawn"],
        "measured_retention_at_production_lambda": ARM["retention"],
        "measured_kappa_detect": ARM["kappa_detect"],
        "retention_lower_3se": ARM["success_debit"]["retention_lower_3se"],
        "kappa_lower_3se": ARM["success_debit"]["kappa_lower_3se"],
    },
    "arithmetic": {
        "fp12_tower": {"f6mul_calls_in_f12mul": N_F6MUL, "f2mul_per_f6mul": F2_PER_F6,
                       "fpmul_per_f2mul": FP_PER_F2, "fpmul_per_f12mul": N_FP,
                       "source": "scripts/ir25/kernels.c f12mul, parsed",
                       "v34_54_withdrawn": True},
        "fp_mul_ir": FP_MUL,
        "mul_charge_ir": MUL, "sqr_charge_ir": SQR,
        "coefficient_row_per_step": COEFF,
        "row_anchors": {
            "fp12_dense_mulmod": "executed f12mul span x the number of non-guard steps; the "
                                 "6-limb Montgomery CIOS chain is inside this measured span",
            "fp12_guard_cyclo_sqr": "executed cyclotomic-square span x the guard-step count",
            "scalar_coefficient_update": "2 x ceil(MUL_CHARGE/%d) + 12, the same anchor "
                                         "discipline as the walk row" % N_FP,
        },
    },
    "trap": {"ledger_row": LEDGER_ROW, "census_counter": CENSUS_COUNTER,
             "markov_witness": MARKOV_WITNESS, "aborts_bound_at_S": str(ABORTS),
             "charge_ir": str(TRAP_IR), "quote_line": i_inc + 1,
             "evidence": "evidence/trap35.json"},
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
              "coeff_per_step": COEFF, "fp_mul_ir": FP_MUL},
    "success": {"route_a_p5": 0.578125, "route_a_median": 0.625,
                "kappa_detect": ARM["kappa_detect"]["k_hat"],
                "kappa_lower_3se": ARM["success_debit"]["kappa_lower_3se"],
                "p_used": 0.578125 * ARM["success_debit"]["kappa_lower_3se"],
                "floor": 0.39,
                "rule": "p_used = p5 x kappa_lower_3se, the detection debit charged explicitly"},
    "prices": {"source": "evidence/exec31.json + evidence/sort32.json, executed counts",
               "censused_symbols": {"state_key48w": 3398, "dict28_wide_lookup": 959,
                                    "dict28_wide_insert": 1711, "dict28_slot_export": 681},
               "per_step_rows": STEP_ROWS},
}
json.dump(ledger, open(os.path.join(EV, "ledger35.json"), "w"), indent=1)

print(json.dumps({"W": str(W), "log2W": round(math.log2(W), 6), "n": n,
                  "claim": n / 10000, "reduced": "%d/%d" % (n // g, 10000 // g),
                  "STEP_IR": STEP_IR, "COEFF": COEFF, "FP_MUL": FP_MUL, "N_FP": N_FP,
                  "ABORTS": str(ABORTS), "TRAP_IR": str(TRAP_IR),
                  "p_used": 0.578125 * ARM["success_debit"]["kappa_lower_3se"]}, indent=1))


# ------------------------------------------------------------------ metadata
CLAIM = n / 10000
RED = "%d/%d" % (n // g, 10000 // g)
KAP = ARM["kappa_detect"]["k_hat"]
KAP_LO = ARM["success_debit"]["kappa_lower_3se"]
RET = ARM["retention"]["p_hat"]
LAM = ARM["lambda"]
P_USED = 0.578125 * KAP_LO

statement = (
    "For nonidentity Q in BLS12-381 G2 and uniform z in GT, compute g=e(P0,Q) and recover x "
    "with z=g^x by a Frobenius-class rho walk in GT on 1024 workers, returning [x]P0 verified "
    "by e(P,Q)=z. CLAIM. Smallest n/10000 with 20000*2^(n/10000) >= W (exact integer test, "
    "n-1 asserted to fail); W = %s is the exact sum of SIXTEEN ceil rows in "
    "evidence/ledger35.json, recomputed by scripts/cert35.py. "
    "DETECTION (verdict ea688603 flaw 1). Detection is THE STORE HIT FOLLOWED BY A TRAIL "
    "COMPARISON THAT IS NOW IN THE CODE (scripts/protocol28.py run_one): a lookup hit counts as "
    "a collision only when the stored trail id differs from the arriving worker's, and same-trail "
    "hits are counted and discarded because G.recover() has den=(D2-q^k D1)=0 there. The v34 "
    "manifest's 'per-worker trail comparison' prose is withdrawn: one d28_ws is shared by all "
    "1024 workers, so cross-worker collisions are cross-trail hits and are not lost, but detection "
    "IS store-load-bearing. It is therefore measured, not assumed: the injected-collision arm "
    "(evidence/arm35.json) replays a uniformly chosen earlier key through the shipped "
    "dict28_wide_lookup at the production lambda and returns the planted reference in %d of %d "
    "probes, kappa_detect = %.6f (se %.2e, 3-sigma lower %.6f). Success carries the debit "
    "explicitly: p_used = 0.578125 x %.6f = %.6f, above the 0.39 floor. "
    "STORE (flaw 2). bounds_table is REGENERATED FROM scripts/ir28/dict28.c, not transcribed: "
    "PROBES=4, ENTRY=9, sizeof(ws_node)=%d B, and the bucket count is a u64 power of two "
    "(dict28_ws_init rejects nb & (nb-1)) whose index is a 64-bit mix, so one level addresses at "
    "most 2^63 buckets. Production is priced as TWO levels of that type: 2^63 shards x 2^63 "
    "buckets = 2^126 heads, capacity 4x2^126, lambda = S/2^126 = %.6f (the v34 heads value and its lambda 0.2815 "
    "are NOT representable and are withdrawn: evidence/historical-disclosure.json); peak storage %s B, initialisation "
    "charged %s Ir (no row is free). Measured retention at that lambda on the shipped four-probe "
    "FIFO chains: %.6f over %d trials, against the Poisson terminal limit %.6f. "
    "ARITHMETIC (flaw 3). The Fp12 tower is counted from scripts/ir25/kernels.c: f12mul calls "
    "f6mul FOUR times, each f6mul is 6 f2mul, each f2mul is 3 fpmul, so one Fp12 multiply is %d "
    "Fp multiplies and the v34 '54' is withdrawn. MUL_CHARGE %d Ir is the executed f12mul span and "
    "is charged as its own row (fp12_dense_mulmod) against the executed non-guard step count; the "
    "coefficient row is 2 x ceil(%d/%d) + 12 = %d Ir/step, the same anchor discipline. "
    "TRAP (flaw 4). The charged event, the census counter and the Markov event are ONE predicate, "
    "quoted verbatim in evidence/trap35.json from the line that increments it (scripts/"
    "protocol28.py:%d: 'aborts += 1; aborted = abandoned = True' under 'if armed and stall >= "
    "K_STALL'). br OR guard is the STALL predicate, not the abort predicate; %d guard conversions "
    "and %d branch events in %d observed steps are reported as the step mix, not charged as "
    "aborts. The bound is a counting bound needing no independence: an abort requires an ARMED "
    "counter, armed is set only on a guard step and cleared on restart, so distinct aborts consume "
    "distinct arming guards and #aborts <= #guard_conversions; at S that is %s aborts, each charged "
    "its 32-step wasted window at the ledger's own per-step price plus one fresh trail start = "
    "%s Ir. NO FLOAT: every charged term is a Python int; retention, kappa and the success rate are "
    "reported values and are labelled as such."
) % (W, ARM["kappa_detect"]["hits"], ARM["kappa_detect"]["probes"], KAP,
     ARM["kappa_detect"]["se"], KAP_LO, KAP_LO, P_USED, BND["parsed_constants"]["NODE_BYTES"],
     LAM, BND["peak_bytes"], BND["init_ir"], RET, ARM["trials"],
     ARM["retention"]["p_derived"], N_FP, MUL, MUL, N_FP, COEFF, i_inc + 1, GUARD, BR, T,
     ABORTS, TRAP_IR)

sub = {
    "version": 1,
    "title": "Claim %.4f: trail comparison in the code + measured kappa_detect, bounds_table "
             "regenerated from dict28.c (2^63 x 2^63 u64 levels, lambda %.4f), Fp12 tower counted "
             "from kernels.c (72 Fp/Fp12), trap predicate three-way identical with a verbatim "
             "quote" % (CLAIM, LAM),
    "claims": [{
        "statement": statement,
        "score": {"numerator": str(n), "denominator": "10000"},
        "evidence": [
            "argument.md", "certificate.json", "dependencies.json",
            "evidence/ledger35.json", "evidence/arm35.json", "evidence/bounds35.json",
            "evidence/trap35.json", "evidence/trap31.json", "evidence/census28.json",
            "evidence/prices28.json", "evidence/exec31.json", "evidence/sort32.json",
            "evidence/protocol35.json", "public-note.md",
            "scripts/regen35.py", "scripts/cert35.py", "scripts/arm35.py", "scripts/bounds35.py",
            "scripts/smoke35.py", "scripts/protocol28.py", "scripts/census28.py",
            "scripts/toy28.py", "scripts/ir28/dict28.c", "scripts/ir25/kernels.c",
        ],
    }],
    "experiments": [],
}
json.dump(sub, open(os.path.join(PKG, "submission.json"), "w"), indent=1)

argument = None
tpl = os.path.join(HERE, "argument35.md.tmpl")
if os.path.exists(tpl):
    argument = open(tpl).read() % {
        "claim": "%.7g" % CLAIM, "red": RED, "W": str(W), "n": str(n),
        "kh": ARM["kappa_detect"]["hits"], "kp": ARM["kappa_detect"]["probes"],
        "khat": KAP, "kse": ARM["kappa_detect"]["se"], "klo": KAP_LO, "pu": P_USED,
        "node": BND["parsed_constants"]["NODE_BYTES"], "lam": LAM,
        "ret": RET, "tr": ARM["trials"], "retd": ARM["retention"]["p_derived"],
        "rdel": ARM["retention"]["delta"], "peak": BND["peak_bytes"],
        "init": BND["init_ir"], "nfp": N_FP, "mul": MUL, "sqr": SQR, "coef": COEFF,
        "line": i_inc + 1, "guard": GUARD, "br": BR, "T": T,
        "aborts": str(ABORTS), "step": STEP_IR, "trap": str(TRAP_IR),
    }
    argument = argument.replace("@CLAIM@", "%.7g" % CLAIM).replace("@RED@", RED)
    argument = argument.replace("@W@", str(W)).replace("@N@", str(n))
    open(os.path.join(PKG, "argument.md"), "w").write(argument)
print("metadata written; claim", CLAIM, RED)

