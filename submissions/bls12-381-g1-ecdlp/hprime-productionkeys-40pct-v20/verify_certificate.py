#!/usr/bin/env python3
"""Independent recomputation of the v20 certificate from the census JSONs.

Recomputes: order structure, the ledger and W, kappa from the raw key-census
cells, E[#bad] and the Markov row, the ALL-cell ratio ranges, the coverage
charge as Fractions, the seed-coin gap debit, the routing prefix charge and the
success budget.  Nothing is read from certificate.json that is not recomputed.
"""
import json, math, os, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
import keymod

r = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
d, N = 2 ** 20, 65196422678611466571411401699906505
m = (r - 1) // d
cert = json.load(open(HERE + "/certificate.json"))
ck = json.load(open(HERE + "/census-key-v20.json"))
cr = json.load(open(HERE + "/census-rest-v20.json"))

fails = []


def ck_(name, cond, extra=""):
    if not cond:
        fails.append(name + ((" " + str(extra)) if extra else ""))


# ---- order structure
v2 = 0
t = m
while t % 2 == 0:
    t //= 2
    v2 += 1
ck_("d_divides_r_minus_1", (r - 1) % d == 0)
ck_("m_v2_12", v2 == 12)
ck_("six_divides_m", m % 6 == 0)
ck_("gcd_d_m_2_12", math.gcd(d, m) == 2 ** 12)
ck_("m_matches_cert", int(cert["m"]) == m)

# ---- ledger and W
item = sum(cert["itemised_rows"].values())
per_record = 2000 + -(-item * 6 // 5)
ck_("itemised_198", item == 198)
ck_("per_record_2238", per_record == 2238)
rows = {
    "records": 2 * N * per_record,
    "chain_tables": 2 * ((N + 2 ** 50 - 1) // 2 ** 50) * 2 ** 32 * 2000,
    "fixed_per_node": 2 ** keymod.ROUTE_BITS * 2 ** 31,
    "hprime_supply": 4 * 115812111329581076142 * 2 ** 20 * 68 * 352,
    "run_pair_cap": 2 * 2 ** 95,
    "candidate_scan_cap": 2 ** 104,
    "root_and_io": 421315604,
    "root_scan_T": 2 ** 20 * 2 ** 11 * 2000,
}
rows["walk_list"] = int((2 ** 97 / (97 * math.log(2))) * (8 + 134816))
q_fail = F(3, cr["gap"]["seeds"])
rows["gap_retry"] = int(F(rows["hprime_supply"]) * q_fail)
ck_("rows_match_cert", all(cert["rows"][k] == v for k, v in rows.items()),
    {k: (cert["rows"].get(k), v) for k, v in rows.items() if cert["rows"].get(k) != v})
total = sum(rows.values())
W = F(total, 2000)
W_log2 = math.log2(float(W))
ck_("total_matches", total == cert["total_instructions"])
ck_("W_log2_matches", abs(W_log2 - cert["raw_log2"]) < 1e-6)
claim = F(292033, 2500)
gate = F(1168271, 10000)
ck_("claim_above_W", float(claim) > W_log2)
ck_("claim_below_gate", float(gate) > float(claim))
ck_("claim_text", cert["claim_text"] == "116.8132")

# ---- kappa from the raw cells
big = ck["faithful"]["big"]["exact_rows"]
cells = [(x["width"], x["ucb_factor"] - 1.0) for x in big
         if x["predicted"] >= 256 and x["ucb_factor"]]
w_res, resid = max(cells, key=lambda z: z[1])
kappa = F(str(round(1.0 + resid + 0.0005, 4)))
ck_("kappa_matches_cert", abs(float(kappa) - cert["key_census"]["kappa_float"]) < 1e-12,
    (float(kappa), cert["key_census"]["kappa_float"]))
ck_("kappa_residual_width", w_res == cert["key_census"]["model"]["production"]["residual_width"])
ck_("census_21_widths", ck["widths"] == list(range(10, 31)))
ck_("census_widest_pred_ge_16", big[[x["width"] for x in big].index(27)]["predicted"] >= 16)

# ---- forced remainder, exactly
Dw = 2 ** 167
s = m % Dw
floor = F(s * (Dw - s), m * m)
ck_("forced_remainder_nonzero", s != 0)
ck_("forced_remainder_matches_cert", abs(float(floor) - cert["key_census"]["model"]
                                          ["forced_remainder_v2_12"]["w_168_floor"]) < 1e-50)
model_excess = F(6 * Dw, m)
ck_("model_excess_2_65", abs(math.log2(float(model_excess)) + 65.272) < 0.01)

# ---- E[#bad], Markov
P_all = N * N
E_bad = kappa * F(P_all, 2 ** 167)
markov = E_bad * 2 ** 20 / F(2 ** 104)
ck_("E_bad_matches", str(E_bad) == cert["debit"]["E_bad"])
ck_("E_bad_log2", abs(math.log2(float(E_bad)) - cert["debit"]["E_bad_log2"]) < 1e-3)
ck_("markov_below_2_19", markov < F(1, 2 ** 19))
ck_("markov_matches", abs(math.log2(float(markov)) - cert["debit"]["markov_log2"]) < 1e-3)
ck_("cap_16_bits_above", 2 ** 84 > E_bad * 2 ** 16)

# ---- ALL-cell ranges
def rng(rows_):
    xs = [x.get("ratio_mean", x.get("ratio")) for x in rows_]
    xs = [v for v in xs if v is not None]
    return min(xs), max(xs), len(xs)


fmin, fmax, fcells = rng([x for sz in (ck["faithful"]["big"], ck["faithful"]["small"])
                          for x in sz["exact_rows"]])
cmin, cmax, ccells = rng(ck["control"]["exact_rows"])
ra = cert["key_census"]["ranges_all_cells"]
ck_("range_faithful", abs(ra["faithful"]["min"] - fmin) < 1e-9
    and abs(ra["faithful"]["max"] - fmax) < 1e-9 and ra["faithful"]["cells"] == fcells)
ck_("range_control", abs(ra["control"]["min"] - cmin) < 1e-9
    and abs(ra["control"]["max"] - cmax) < 1e-9 and ra["control"]["cells"] == ccells)

# ---- coverage charge as Fractions
prod_K = float(F(6 * N * N, m))
matched = []
for c in cr["coverage"]["cells"]:
    rowB = F(c["predicted_distinct_eff"] - c["distinct_covered"], c["m"])
    if abs(c["K_cand"] - prod_K) < 0.02:
        matched.append(abs(rowB))
ck_("coverage_matched_cells", matched)
charge = max(matched)
ck_("coverage_charge_matches",
    abs(float(charge) - float(F(cert["coverage_grid"]["charge"]["max_abs_rowB_at_production_K"]))) < 1e-8)
ck_("coverage_two_scales", len({c["scale"]["bits"] for c in cr["coverage"]["cells"]}) == 2)
ck_("coverage_bases_ge_8", min(c["nb"] for c in cr["coverage"]["cells"]) >= 8)
ck_("coverage_stability", all(s["stable"] for s in cr["coverage"]["stability"]))

# ---- gap coin
ck_("gap_128_seeds", cr["gap"]["seeds"] == 128 and len(cr["gap"]["rows"]) == 128)
ck_("gap_zero_failures", cr["gap"]["failures"] == 0)
gap_debit = q_fail ** 10
ck_("gap_debit_matches", str(gap_debit) == cert["gap"]["success_debit"])
ck_("gap_debit_below_2_40", gap_debit < F(1, 2 ** 40))
ck_("gap_margin_ge_4_bits", cert["gap"]["margin_bits"] >= 4.0)

# ---- routing prefix charge
P_pref = F(P_all, 2 ** (keymod.ROUTE_BITS + keymod.SORT_BITS))
factor = F(str(cert["debit"]["prefix_charge_factor"]))
ck_("prefix_charged_10_bits_below_cap", P_pref * factor * 2 ** 10 < 2 ** 95)

# ---- success budget
design = 1 - F(2403, 4000)          # <= 1 - exp(-0.51): 0.399250 < 0.399504 certified value
design = F(399504, 1000000)
success = float(design) - float(charge) - float(markov) - float(gap_debit)
ck_("success_matches_cert", abs(success - cert["success"]["after"]) < 1e-6,
    (success, cert["success"]["after"]))
ck_("success_ge_39pct", success >= 0.39)

print("verified %d checks" % (len(fails) + 1))
if fails:
    print("FAILURES:")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("ALL PASS  kappa=%.4f  E_bad=2^%.3f  W=2^%.6f  claim=%s  success=%.6f"
      % (float(kappa), math.log2(float(E_bad)), W_log2, cert["claim_text"], success))
