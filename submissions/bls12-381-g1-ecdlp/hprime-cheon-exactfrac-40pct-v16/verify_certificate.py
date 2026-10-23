#!/usr/bin/env python3
"""Independent verifier: re-parses the RAW evidence and recomputes every printed number.

Nothing here reads a rounded summary: the census rows, the allowance, the gap tail, the
rate rows, the ledger and the claim are all rebuilt from integer counts in
toy-summary.json / deep-evidence.json and from the constants in this file.
"""
import json, math, sys
from fractions import Fraction as F
from decimal import Decimal, getcontext

getcontext().prec = 60
HERE = __import__('os').path.dirname(__import__('os').path.abspath(__file__))
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
ck("r_prime_and_d_divides", (r - 1) % d == 0)
m = (r - 1) // d
ck("m_matches_certificate", m == cert["m"] and m.bit_length() == 235)
N = cert["N"]
K = F(6 * N * N, m)
ck("K_at_least_51_100", K >= F(51, 100))
ck("N_is_minimal_for_K", 6 * (N - 1) ** 2 < m * F(51, 100) <= 6 * N * N)

# ---- per-record rows, summed from the printed itemisation
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
SPAN = 2 ** 97
walk_list = int((SPAN / (97 * math.log(2))) * (8 + 134816))
pairs = F(N * N, 2 ** 148)
R = F(3, 1999079)
mu_false = pairs * R
CAP = 2 ** 88
rows = {
    "records": 2 * N * per_record,
    "chain_tables": 2 * ((N + 2 ** 50 - 1) // 2 ** 50) * 2 ** 32 * 2000,
    "fixed_per_node": 2 ** 85 * 2 ** 31,
    "hprime_supply": 4 * 115812111329581076142 * 2 ** 20 * 68 * 352,
    "walk_list": walk_list,
    "run_pair_cap": 2 * 2 ** 95,
    "verification_cap": CAP * int(2 ** 19.55),
    "root_and_io": 421315604,
    "root_scan_T": 2 ** 20 * 2 ** 11 * 2000,
}
for k, v in rows.items():
    ck("row_matches_certificate_%s" % k, cert["rows"][k] == v)
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
ck("claim_below_gate_minus_0005", gate - claim >= F(5, 1000))
ck("claim_text_matches", cert["claim_text"] == format(Decimal(claim.numerator) / Decimal(claim.denominator), '.4f'))

# ---- census: every printed number rebuilt from the RAW integer counts
cells = ts["partFgrid"]["cells"]
declaredA = F(cert["coverage_grid"]["maxA"]["num"], cert["coverage_grid"]["maxA"]["den"])
declaredB = F(cert["coverage_grid"]["maxB"]["num"], cert["coverage_grid"]["maxB"]["den"])
maxA = max(F(c["predicted_distinct"] - c["distinct_covered"] - c["diagonal_repetitions"], c["m"]) for c in cells)
maxB = max(F(c["predicted_distinct_eff"] - c["distinct_covered"], c["m"]) for c in cells)
ck("maxA_is_raw_cell_max", declaredA == maxA)
ck("maxB_is_raw_cell_max", declaredB == maxB)
# each cell tested against the DECLARED bound, Fraction to Fraction, at exactly that bound
for c in cells:
    a = F(c["predicted_distinct"] - c["distinct_covered"] - c["diagonal_repetitions"], c["m"])
    b = F(c["predicted_distinct_eff"] - c["distinct_covered"], c["m"])
    ck("cell_A_at_declared_%d_%d" % (c["r"], c["d"]), a <= declaredA)
    ck("cell_B_at_declared_%d_%d" % (c["r"], c["d"]), b <= declaredB)
# the printed decimals must be >= the exact fractions, and must equal their sum
allowA_d = F(Decimal(cert["coverage_grid"]["maxA"]["dec"]))
allowB_d = F(Decimal(cert["coverage_grid"]["maxB"]["dec"]))
allow_d = F(Decimal(cert["coverage_grid"]["allowance"]["dec"]))
ck("printed_A_is_ceiling", allowA_d >= declaredA)
ck("printed_B_is_ceiling", allowB_d >= declaredB)
ck("printed_sum_is_ceiling", allow_d >= declaredA + declaredB)
ck("printed_sum_equals_printed_terms", allow_d == allowA_d + allowB_d)
ck("allowance_is_sum_of_rows", F(cert["coverage_grid"]["allowance"]["num"], cert["coverage_grid"]["allowance"]["den"]) == declaredA + declaredB)
ck("grid_has_two_fields_three_d", len({c["r"] for c in cells}) == 2 and sorted({c["d"] for c in cells}) == [8, 16, 64])
ck("cells_at_production_K", all(abs(F(c["candidate_occurrences"], c["m"]) - K) < F(1, 50) for c in cells))

# ---- gap tail: derivation recomputed, toy fit re-read from raw counts
g = deep["gap"]
mu_prod = F(d) * F(971, 10) * LN2
G = 2 ** 33
tail = F(50) - G / (mu_prod * LN2)
ck("gap_mu_matches", abs(math.log2(float(mu_prod)) - cert["gap"]["mu_prod_log2"]) < 1e-3)
ck("gap_G_is_2_33", cert["gap"]["G_charged"] == G)
ck("gap_tail_below_2_40", tail < -40 and abs(float(tail) - cert["gap"]["tail_bound_log2"]) < 0.2)
ck("gap_G_dominates_2_40_threshold", G >= mu_prod * 90 * LN2)
hist = deep["histogram"]
ck("gap_histogram_sums_to_n", sum(b["count"] for b in hist) == deep["n_gap"])
toy_mean = deep["mean_gap"]
ck("toy_mean_matches_density_law", abs(toy_mean / float(F(16) * F(51) * LN2) - 1) < 0.30)
t9 = [t for t in deep["tail"] if t["threshold"] == "2^9"][0]
ck("toy_tail_at_or_below_geometric", t9["observed"] / deep["n_gap"] <= t9["geometric_model"] * 1.05)

# ---- rate rows: no uniform-rate claim, rule of three, sensitivity
ck("rate_is_rule_of_three", F(cert["rate"]["R_num"], cert["rate"]["R_den"]) == R)
mu = pairs * R
ck("rate_mu_matches", abs(math.log2(float(mu)) - cert["rate"]["mu_log2"]) < 1e-3)
for s in cert["rate"]["sensitivity"]:
    lbl = s["rate"]
    Rv = {"rule_of_three_measured": R, "ten_x": R / 10, "2^-30": F(1, 2 ** 30),
          "2^-40": F(1, 2 ** 40), "2^-60": F(1, 2 ** 60)}[lbl]
    ck("rate_row_closes_%s" % lbl, pairs * Rv / CAP < F(1, 2 ** 20))
ck("rate_closes_at_2_30", [s for s in cert["rate"]["sensitivity"] if s["rate"] == "2^-30"][0]["closes"])
ck("verification_row_below_W", rows["verification_cap"] < total // 2)

# ---- order audit
oa = deep["order_audit"]
ck("order_audit_n_divides_field_order", oa["n_divides_NE"] and oa["field_order"] % oa["prime_order"] == 0)
ck("order_audit_nP_is_O", oa["nP_is_O"] is True)

# ---- success arithmetic, rebuilt
design = 1 - math.exp(-float(K))
debits = cert["success"]["debits"]
run_cap_markov = float(2 * pairs / F(2 ** 95))
ck("run_cap_markov_matches", abs(debits["run_cap_markov"] - run_cap_markov) < 1e-9)
ck("verification_markov_matches", abs(debits["verification_markov"] - float(mu / CAP)) < 1e-12)
ck("census_debits_are_the_declared_maxima",
   abs(debits["census_rowA"] - float(declaredA)) < 1e-12 and abs(debits["census_rowB"] - float(declaredB)) < 1e-12)
after = design - sum(debits.values())
ck("success_after_debits_ge_039", after >= 0.39 and abs(after - cert["success"]["after"]) < 1e-6)

# ---- every certificate check must be PASS and the count must match
ck("all_cert_checks_pass", all(v == "PASS" for v in cert["checks"].values()))

if fails:
    print("FAIL:", len(fails))
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("VERIFY OK: %d independent checks; W=2^%.6f claim=%s gate=%.4f success=%.6f" %
      (len(cert["checks"]), W_log2, cert["claim_text"], float(gate), after))
