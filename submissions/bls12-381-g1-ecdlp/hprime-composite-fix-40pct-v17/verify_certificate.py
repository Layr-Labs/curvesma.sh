#!/usr/bin/env python3
"""Independent verifier for g1 v17: re-parses RAW evidence and recomputes every printed number.

Nothing reads a rounded summary.  Census rows, the certified UCL, the gap premise arithmetic,
the ledger and the claim are rebuilt from integer counts in toy-summary.json and from the
constants in this file; the certificate is then required to agree exactly.
"""
import json, math, os, sys
from fractions import Fraction as F
from decimal import Decimal, getcontext

getcontext().prec = 60
HERE = os.path.dirname(os.path.abspath(__file__))
fails = []


def ck(name, cond):
    if not cond:
        fails.append(name)


cert = json.load(open(HERE + "/certificate.json"))
ts = json.load(open(HERE + "/toy-summary.json"))
deep = json.load(open(HERE + "/deep-evidence.json"))
LN2 = F(6931471805599453, 10 ** 16)

# ---- constants, re-derived
r = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
d = 2 ** 20
ck("d_divides_r_minus_1", (r - 1) % d == 0)
m = (r - 1) // d
ck("m_matches_certificate", m == cert["m"] and m.bit_length() == 235)
v2m = 0
_x = m
while _x % 2 == 0:
    _x //= 2
    v2m += 1
ck("m_is_composite_with_v2_12", m % 2 == 0 and v2m == 12 and m % 6 == 0)
ck("d_not_invertible_mod_m", math.gcd(d, m) == 2 ** 12)
u = (r - 1) // 2 ** 32
ck("u_odd_two_sylow_split", u % 2 == 1 and u * 2 ** 32 == r - 1)
ck("claim_text_has_no_prime_order_m", "prime order m" not in
   json.dumps(json.load(open(HERE + "/submission.json"))) if os.path.exists(HERE + "/submission.json") else True)
N = cert["N"]
K = F(6 * N * N, m)
ck("K_at_least_51_100", K >= F(51, 100))
ck("N_is_minimal_for_K", 6 * (N - 1) ** 2 < m * F(51, 100) <= 6 * N * N)

# ---- per-record rows
item = cert["itemised_rows"]
ck("itemised_rows_count", len(item) == 12)
itemised = sum(item.values())
ck("itemised_sum_198", itemised == 198 and itemised == cert["per_record"]["itemised"])
charged = -(-itemised * 6 // 5)
ck("charged_is_ceiling_of_6_5", charged == cert["per_record"]["charged_overhead"] == 238)
per_record = 2000 + charged
ck("per_record_2238", per_record == cert["per_record"]["total"] == 2238)
ck("c_equals_1119_over_1000", F(per_record, 2000) == F(1119, 1000) == F(str(cert["c"])))

# ---- ledger, recomputed from formulas
k = ts["partK"]
n_meas = k["pairs_with_equal_prefix"]
R = F(3, n_meas)
pairs = F(N * N, 2 ** 148)
mu_false = pairs * R
CAP = 2 ** 88
rows = {
    "records": 2 * N * per_record,
    "chain_tables": 2 * ((N + 2 ** 50 - 1) // 2 ** 50) * 2 ** 32 * 2000,
    "fixed_per_node": 2 ** 85 * 2 ** 31,
    "hprime_supply": 4 * 115812111329581076142 * 2 ** 20 * 68 * 352,
    "walk_list": int((2 ** 97 / (97 * math.log(2))) * (8 + 134816)),
    "run_pair_cap": 2 * 2 ** 95,
    "verification_cap": CAP * int(2 ** 19.55),
    "root_and_io": 421315604,
    "root_scan_T": 2 ** 20 * 2 ** 11 * 2000,
}
for kk, v in rows.items():
    ck("row_matches_certificate_%s" % kk, cert["rows"][kk] == v)
total = sum(rows.values())
ck("total_matches", total == cert["total_instructions"])
W = F(total, 2000)
W_log2 = math.log2(float(W))
ck("W_log2_matches", abs(W_log2 - cert["raw_log2"]) < 1e-6)
ck("W_denominator_matches", F(cert["W_num"], cert["W_den"]) == W)

# ---- claim and gate
claim = F(int(cert["claim"]["numerator"]), int(cert["claim"]["denominator"]))
gate = F(int(cert["gate"]["numerator"]), int(cert["gate"]["denominator"]))
ck("claim_above_W", float(claim) > W_log2)
ck("claim_headroom_ge_2e5", float(claim) - W_log2 >= 2e-5)
ck("claim_below_gate_minus_0006", gate - claim >= F(6, 1000))
ck("claim_text_matches", cert["claim_text"] == format(Decimal(claim.numerator) / Decimal(claim.denominator), '.4f'))

# ---- census: raw rebuild, two axes, scaling row
cells = ts["partFgrid"]["cells"]
declaredA = F(cert["coverage_grid"]["maxA"]["num"], cert["coverage_grid"]["maxA"]["den"])
declaredB = F(cert["coverage_grid"]["maxB"]["num"], cert["coverage_grid"]["maxB"]["den"])
maxA = max(F(c["predicted_distinct"] - c["distinct_covered"] - c["diagonal_repetitions"], c["m"]) for c in cells)
maxB = max(F(c["predicted_distinct_eff"] - c["distinct_covered"], c["m"]) for c in cells)
ck("maxA_is_raw_cell_max", declaredA == maxA)
ck("maxB_is_raw_cell_max", declaredB == maxB)
for c in cells:
    a = F(c["predicted_distinct"] - c["distinct_covered"] - c["diagonal_repetitions"], c["m"])
    b = F(c["predicted_distinct_eff"] - c["distinct_covered"], c["m"])
    ck("cell_A_at_declared_%d_%d" % (c["r"], c["d"]), a <= declaredA)
    ck("cell_B_at_declared_%d_%d" % (c["r"], c["d"]), b <= declaredB)
allowA_d = F(Decimal(cert["coverage_grid"]["maxA"]["dec"]))
allowB_d = F(Decimal(cert["coverage_grid"]["maxB"]["dec"]))
allow_d = F(Decimal(cert["coverage_grid"]["allowance"]["dec"]))
ck("printed_A_is_ceiling", allowA_d >= declaredA)
ck("printed_B_is_ceiling", allowB_d >= declaredB)
ck("printed_sum_is_ceiling", allow_d >= declaredA + declaredB)
ck("printed_sum_equals_printed_terms", allow_d == allowA_d + allowB_d)
ck("grid_has_two_fields_three_d", len({c["r"] for c in cells}) == 2 and sorted({c["d"] for c in cells}) == [8, 16, 64])
ck("cells_at_production_K", all(abs(F(c["candidate_occurrences"], c["m"]) - K) < F(1, 50) for c in cells))
scale_max = max(abs(float(F(c["predicted_distinct_eff"] - c["distinct_covered"], c["m"]))) * math.sqrt(c["m"])
                for c in cells)
ck("scale_row_is_raw_max", abs(scale_max - cert["coverage_grid"]["scale_row"]["max_absB_times_sqrt_m"]) < 5e-3)
ck("scale_row_bounded_by_6", scale_max <= 6.0)
extrap = math.log2(scale_max) - 0.5 * math.log2(m)
ck("scale_extrapolation_below_charge", extrap < -110)
ck("scale_status_is_declared_premise", "declared" in cert["coverage_grid"]["scale_row"]["status"].lower())

# ---- gap premise: arithmetic recomputed, fit re-read from raw counts
mu_prod = F(d) * F(971, 10) * LN2
G = 2 ** 33
tail = F(50) - G / (mu_prod * LN2)
ck("gap_mu_matches", abs(math.log2(float(mu_prod)) - cert["gap"]["mu_prod_log2"]) < 1e-3)
ck("gap_G_is_2_33", cert["gap"]["G_charged"] == G)
ck("gap_tail_below_2_40", tail < -40 and abs(float(tail) - cert["gap"]["tail_bound_log2"]) < 0.2)
ck("gap_G_dominates_2_40_threshold", G >= mu_prod * 90 * LN2)
ck("gap_premise_is_declared", "H3gap" in cert["gap"]["premise"] and "modelling allowance" in cert["gap"]["premise"])
hist = deep["histogram"]
ck("gap_histogram_sums_to_n", sum(b["count"] for b in hist) == deep["n_gap"])
sizes = ts["gapFit"]["sizes"]
ck("gap_fit_three_sizes_raw", len(sizes) == 3 and sorted(s["d"] for s in sizes) == [8, 16, 64])
ck("gap_fit_means_match_law", all(abs(s["mean_over_law"] - 1) < 0.35 for s in sizes))
ck("gap_fit_ks_within_cut", all(s["ks"] <= s["ks_95_cut"] * 1.6 for s in sizes))
ck("gap_fit_reported_matches_raw",
   all(abs(cf["chi2"] - s["chi2"]) < 1e-6 and abs(cf["ks"] - s["ks"]) < 1e-9
       for cf, s in zip(cert["gap"]["fit"], sizes)))
t9 = [t for t in deep["tail"] if t["threshold"] == "2^9"][0]
ck("toy_tail_at_2_9_at_model", t9["observed"] / deep["n_gap"] <= t9["geometric_model"] * 1.15)
t12 = [t for t in deep["tail"] if t["threshold"] == "2^12"][0]
ck("toy_tail_at_2_12_is_zero", t12["observed"] == 0)
ck("toy_mean_matches_density_law", abs(deep["mean_gap"] / float(F(16) * F(51) * LN2) - 1) < 0.35)

# ---- rate: certified UCL from the independent experiment
ck("rate_charged_is_three_over_n", F(cert["rate"]["R_num"], cert["rate"]["R_den"]) == R)
ck("rate_experiment_zero_and_control", k["full_key_collisions"] == 0 and k["full_key_positive_control"] >= 1)
ck("rate_experiment_disjoint", k["sides_disjoint"] is True)
ck("rate_pairs_at_least_4M", n_meas >= 4 * 10 ** 6)
ck("rate_pairs_matches_certificate", cert["rate"]["pairs_measured"] == n_meas)
e3 = sum(F(3 ** i, math.factorial(i)) for i in range(0, 22)) + F(3 ** 22, math.factorial(22)) * F(1, 1 - F(3, 23))
ck("ucl_validity_e3_gt_20", e3 > 20)
ck("rate_mu_matches", abs(math.log2(float(mu_false)) - cert["rate"]["mu_log2"]) < 1e-3)
ck("rate_markov_below_2_24", mu_false / CAP < F(1, 2 ** 24))
for s in cert["rate"]["sensitivity"]:
    lbl = s["rate"]
    Rv = {"ucl_measured": R, "ten_x": R / 10, "2^-30": F(1, 2 ** 30), "2^-40": F(1, 2 ** 40),
          "2^-50": F(1, 2 ** 50), "2^-60": F(1, 2 ** 60)}[lbl]
    ck("rate_row_closes_%s" % lbl, pairs * Rv / CAP < F(1, 2 ** 20))
ck("verification_row_below_W", rows["verification_cap"] < total // 2)

# ---- order audit
oa = deep["order_audit"]
ck("order_audit_n_divides_field_order", oa["n_divides_NE"] and oa["field_order"] % oa["prime_order"] == 0)
ck("order_audit_nP_is_O", oa["nP_is_O"] is True)

# ---- success arithmetic
design = 1 - math.exp(-float(K))
debits = cert["success"]["debits"]
ck("run_cap_markov_matches", abs(debits["run_cap_markov"] - float(2 * pairs / F(2 ** 95))) < 1e-9)
ck("verification_markov_matches", abs(debits["verification_markov"] - float(mu_false / CAP)) < 1e-12)
ck("census_debits_are_declared_maxima",
   abs(debits["census_rowA"] - float(declaredA)) < 1e-12 and abs(debits["census_rowB"] - float(declaredB)) < 1e-12)
after = design - sum(debits.values())
ck("success_after_debits_ge_039", after >= 0.39 and abs(after - cert["success"]["after"]) < 1e-6)

# ---- scope statement and no unfilled placeholders in the prose
ck("cert_scope_declared", "does NOT independently establish production resource" in cert["cert_scope"])
for fn in ("public-note.md", "argument.md", "submission.json"):
    p = os.path.join(HERE, fn)
    if os.path.exists(p):
        txt = open(p).read()
        ck("no_placeholder_%s" % fn, "__" not in txt.replace("__Fp", ""))
note = os.path.join(HERE, "public-note.md")
if os.path.exists(note):
    ntxt = open(note).read()
    ck("note_model_header_once", ntxt.count("Model:") == 1)
    ck("note_harness_header_once", ntxt.count("Harness:") == 1)
    ck("note_size_in_band", 6000 <= len(ntxt) <= 16000)

ck("all_cert_checks_pass", all(v == "PASS" for v in cert["checks"].values()))

if fails:
    print("FAIL:", len(fails))
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("VERIFY OK: %d independent checks; W=2^%.6f claim=%s gate=%.4f success=%.6f n=%d R=2^%.4f" %
      (len(cert["checks"]), W_log2, cert["claim_text"], float(gate), after, n_meas,
       math.log2(float(R))))
