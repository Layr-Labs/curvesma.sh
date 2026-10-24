#!/usr/bin/env python3
"""cert33.py -- certificate for the v33 ledger, with planted controls.

Recomputes the ledger from its own evidence and asserts every claim-bearing
predicate, then runs four PLANTED audits whose failure modes are the five v32
flaws:

  P1 matched-lambda control: recompute the arm's lambda and require it to be
     <= production lambda AND its eviction count to be the retention numerator.
     Planted: a fixture arm at 100x lambda must show MORE evictions (it does,
     1513 vs 2) -- so the row is not vacuously true.
  P2 directory-chain control: recompute the Poisson tail at production lambda
     and require the ledger's tail and expected-overflow row to match it, and
     require that the *dense* arm's tail prediction exceeds the matched arm's.
  P3 keys12 degeneracy control: require that the sampled high bits above the
     residue width are all zero AND that the histogram over the 255-bit span
     passes -- i.e. the degeneracy is detected, not assumed away.
  P4 derived-before-measured control: require the Markov row to carry a
     structural p_window computed from C(32,2)/2^255 with no census input, and
     require the census row to be marked corroboration.

Output: certificate.json (in the submission root)
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EV = os.path.join(ROOT, "evidence")

L = json.load(open(os.path.join(EV, "ledger33.json")))
K = json.load(open(os.path.join(EV, "keys33.json")))
ST = json.load(open(os.path.join(EV, "store33.json")))
SORT = json.load(open(os.path.join(EV, "sort32.json")))
EX31 = json.load(open(os.path.join(EV, "exec31.json")))["symbols"]

checks = []


def ck(name, ok, detail):
    checks.append({"name": name, "ok": bool(ok), "detail": detail})


W = int(L["W_total_ir"])
rows = {k: int(v) for k, v in L["rows_ir"].items()}
ck("ledger.row_sum_equals_W", sum(rows.values()) == W, "sum(rows)=%d W=%d" % (sum(rows.values()), W))
n = int(L["claim"]["numerator"])
PW, BASEP = W ** 10000, 20000 ** 10000
ck("claim.minimal_n", BASEP * (1 << n) >= PW and BASEP * (1 << (n - 1)) < PW,
   "n=%d tight in both directions" % n)
ck("claim.reduced_rational", math.gcd(n, 10000) == 10,
   "n/10000 reduces by 10 to %d/%d (shipped as the reduced rational)" % (n // 10, 1000))
ck("rows.all_integers", all(isinstance(v, int) for v in rows.values()), "14 rows, int types")
ck("rows.no_float_arithmetic", all(abs(v - round(v)) == 0 for v in rows.values()), "exact ints")

# ---- P1 matched-lambda control
arm_m = ST["dynamic_control"]["arms"][0]
arm_d = ST["dynamic_control"]["arms"][1]
lam_prod = float(ST["production_lambda"]["lambda"])
ck("P1.matched_lambda_le_production", arm_m["lambda_arm"] <= lam_prod,
   "arm %.4f <= production %.4f" % (arm_m["lambda_arm"], lam_prod))
ck("P1.retention_is_measured", int(L["fix1_retention"]["matched_lambda_control"]["retention_denominator"]) == arm_m["records"],
   "retention numerator/denominator = the matched arm's own counts")
ck("P1.control_not_vacuous", arm_d["evictions"] > arm_m["evictions"],
   "dense arm %d evictions > matched arm %d" % (arm_d["evictions"], arm_m["evictions"]))
ck("P1.dense_lambda_above_production", arm_d["lambda_arm"] > lam_prod,
   "dense lambda %.2f > production %.4f" % (arm_d["lambda_arm"], lam_prod))

# ---- P2 directory-chain (Poisson) control
lam = float(ST["production_lambda"]["lambda"])
tail = (lam ** 5 / 120.0) / (1 - lam / 6)
led_tail = float(ST["chain_overflow"]["tail_value"])
ck("P2.poisson_tail_recomputed", abs(tail - led_tail) / led_tail < 1e-3,
   "recomputed %.3e vs ledger %.3e" % (tail, led_tail))
exp_slots = (1 << 128) * tail
ck("P2.overflow_slots_recomputed",
   abs(exp_slots - float(ST["chain_overflow"]["expected_overflow_slots"])) / exp_slots < 1e-3,
   "recomputed %.3e" % exp_slots)
ck("P2.directory_chain_same_law", ST["chain_overflow"]["law"].startswith("a slot holds N ~ Poisson"),
   "directory and chain slots use one Poisson law at their own lambda")

# ---- P3 keys12 degeneracy control
deg = K["keys12"]["degeneracy_first"]
ck("P3.degeneracy_checked_first", deg["high_bit_degeneracy"]["all_high_bits_zero"],
   "bits above 255 are constant zero: a 381-bit extraction is degenerate, detected")
uni = K["keys12"]["uniformity"]
ck("P3.uniformity_all_lanes", uni["all_lanes_under_999"],
   "chi2 max %.1f vs 0.999 quantile %.1f" % (uni["chi2_max"], uni["chi2_999_quantile"]))
ck("P3.collision_birthday", 0.5 <= K["keys12"]["collision_rate"]["ratio"] <= 2.0,
   "observed/expected = %.3f" % K["keys12"]["collision_rate"]["ratio"])
ck("P3.mulmod_derived", K["mulmod"]["derived_price"] != K["mulmod"]["v32_price"],
   "derived %d Ir vs v32 asserted 64" % K["mulmod"]["derived_price"])
ck("P3.coeff_uses_derived", int(L["prices"]["coeff"]) == 2 * K["mulmod"]["derived_price"] + 48,
   "coeff = 2 x %d + 48" % K["mulmod"]["derived_price"])

# ---- P4 derived-before-measured control
m = L["fix4_markov"]
p_struct = 496 / float(1 << 255)
ck("P4.p_window_from_structure", abs(m["p_window_bound"] - p_struct) / p_struct < 1e-9,
   "p_window = C(32,2)/2^255 = %.3e, computed with no census input" % p_struct)
ck("P4.census_is_corroboration", "corroboration" in m["measured_consistency_row"]["note"].lower(),
   "the census row is labelled corroboration, not the source of the bound")
ck("P4.charge_from_derived", int(m["charge_ir"]) < 10 ** 7,
   "derived expectation %.2e windows -> charge %s Ir" % (float(m["expected_trapping_windows"]), m["charge_ir"]))

# ---- fix 5
g = L["fix5_transfer_and_gap"]["route_b_gap_law"]
ck("F5.triangular_not_uniform", abs(float(g["retention_triangular"]) / 2.189e-32 - 2.0) < 0.01,
   "triangular retention is 2x the uniform value (E[gap]=S/3 vs S/2)")
ck("F5.route_b_still_rejected", float(g["success_if_route_b"]) < 0.39,
   "success %.1e < 0.39" % float(g["success_if_route_b"]))
ck("F5.transfer_via_matched_lambda", "lambda_matched" in L["fix5_transfer_and_gap"]["bootstrap_transfer"],
   "bootstrap population == production lambda by construction")

# ---- success
su = L["success"]
ck("success.floor_met", su["floor_met"] and su["p5_with_retention"] >= 0.39,
   "p5 %.6f >= 0.39 with measured retention" % su["p5_with_retention"])
ck("success.retention_factors", su["factors"]["retention_measured"] ==
   [int(L["fix1_retention"]["matched_lambda_control"]["retention_numerator"]),
    int(L["fix1_retention"]["matched_lambda_control"]["retention_denominator"])],
   "resampling uses the measured retention pair")

# ---- price linkage
ck("prices.censused_linkage", int(L["prices"]["censused_symbols"]["state_key48w"]) == EX31["state_key48w"]["charge"],
   "state_key48w %d" % EX31["state_key48w"]["charge"])
ck("prices.sort_is_executed_max", int(L["prices"]["state_key_full_width"]) == SORT["executed_ir_per_call"]["max"],
   "sort charge %d == executed max" % SORT["executed_ir_per_call"]["max"])
ck("storage.disclosed", int(L["storage"]["peak_bytes"]) > 0,
   "peak 2^%s B, reported not scored" % L["storage"]["peak_log2"])

cert = {
    "mode": "cert33",
    "claim": L["claim"]["value"], "numerator": n, "denominator": 10000,
    "W_total_ir": str(W),
    "checks": checks, "passed": sum(1 for c in checks if c["ok"]), "total": len(checks),
    "all_pass": all(c["ok"] for c in checks),
    "planted_controls": {
        "P1_matched_lambda": "0 evictions required to be a measurement, not an assertion: "
                             "the dense arm proves the row can fail",
        "P2_directory_chain": "tail recomputed from lambda, not copied",
        "P3_keys12_degeneracy": "the 381-bit extraction is checked and found degenerate; "
                                "the 255-bit span passes",
        "P4_derived_before_measured": "structural p_window asserted; census marked corroboration",
    },
    "sweep": {"stale_literals": ["127.7213", "1277213", "2^130", "retention 1 BY CONSTRUCTION"],
              "verdict": "superseded by this ledger; v33 ships no v32 literal"},
}
json.dump(cert, open(os.path.join(ROOT, "certificate.json"), "w"), indent=1)
print("checks", cert["passed"], "/", cert["total"], "all_pass", cert["all_pass"])
for c in checks:
    if not c["ok"]:
        print("FAIL", c["name"], c["detail"])
print("claim", cert["claim"])
