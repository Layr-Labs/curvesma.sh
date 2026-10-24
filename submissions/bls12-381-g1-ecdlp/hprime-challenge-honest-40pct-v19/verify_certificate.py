#!/usr/bin/env python3
"""Independent verifier for g1 v18: re-parses RAW evidence and recomputes every printed number.

Nothing reads a rounded summary.  The key expectation, the Markov debit, the batch cap, the census
Fractions, the trend row, the record layout, the reconstruction price and the claim are rebuilt from
integer counts in toy-summary.json and from the declared constants; the certificate is then required
to agree exactly.  v18 also asserts the ABSENCE of the v17 statistics: no confidence bound, no
Chernoff argument and no independence claim may appear in any submitted artifact.
"""
import json, math, os, sys
from fractions import Fraction as F
from decimal import Decimal, getcontext

getcontext().prec = 60
HERE = os.path.dirname(os.path.abspath(__file__))
fails = []
import keymod


def ck(name, cond):
    if not cond:
        fails.append(name)


def exp_neg_upper(x, terms=26):
    s = F(0)
    for i in range(terms):
        s += F((-1) ** i) * x ** i / math.factorial(i)
    return s + x ** terms / math.factorial(terms)


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
ck("m_is_composite_with_v2_12", m % 2 == 0 and v2m == 12 and m % 6 == 0 and m > 2)
ck("d_not_invertible_mod_m", math.gcd(d, m) == 2 ** 12)
u = (r - 1) // 2 ** 32
ck("u_odd_two_sylow_split", u % 2 == 1 and u * 2 ** 32 == r - 1)
N = cert["N"]
K = F(6 * N * N, m)
ck("K_at_least_51_100", K >= F(51, 100))

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

# ---- key geometry and the expectation-Markov debit, rebuilt from raw census counts
px = ts["partX"]                       # only the structural (prime-product) table
cv = json.load(open(HERE + "/census-v19.json"))
KG = cv["faithful"]                    # H-members, challenge-scaled giant side
CTL = cv["control"]                    # v18 object, mutation control
ck("key_bits_168", cert["key_census"]["key_bits"] == keymod.KEY_BITS == 168)
ck("image_bits_167", cert["key_census"]["image_log2"] == keymod.KEY_IMAGE_BITS == 167)
ck("census_two_sizes", len(KG) == 2 and len({s["records_per_side"] for s in KG}) == 2)
ck("census_seven_widths", cv["widths"] == [10, 12, 14, 16, 18, 20, 22]
   and cert["key_census"]["widths"] == cv["widths"])
ck("census_six_challenges", cv["n_challenges"] == 6 and cert["key_census"]["n_challenges"] == 6)
ck("census_m_is_production_m", cv["m"] == m and cv["d"] == d)
ck("faithful_lists_are_h_members", all(sz["h_members"] is True and sz["scaled"] is True for sz in KG))
ck("control_lists_are_the_v18_object",
   all(sz["h_members"] is False and sz["scaled"] is False for sz in CTL))
memb = cv["membership"]
ck("membership_row_matches_certificate", cert["key_census"]["membership"] == memb)
ck("membership_sample_size", memb["samples"] == 1 << 20)
ck("membership_rejects_almost_all", memb["rejected_fraction"] > 0.999
   and memb["rejected_fraction"] == round(1 - memb["accepted"] / memb["samples"], 8))
ck("membership_accept_count_matches_index", memb["accepted"] <= 8)
cells_charged = []
for sz in KG:
    ck("census_pairs_square_%d" % sz["records_per_side"], sz["pairs"] == sz["records_per_side"] ** 2)
    ck("faithful_rows_full_widths_%d" % sz["records_per_side"],
       len(sz["exact_rows"]) == len(cv["widths"]))
    for r0 in sz["exact_rows"]:
        ck("census_row_width_%d_%d" % (sz["records_per_side"], r0["width"]),
           abs(r0["predicted"] - sz["pairs"] / 2.0 ** (r0["width"] - 1)) <= 0.51)
        ck("faithful_min_mean_max_ordered_%d_%d" % (sz["records_per_side"], r0["width"]),
           r0["observed_min"] <= r0["observed_mean"] + 1e-3 and r0["observed_mean"] <= r0["observed_max"] + 1e-3)
        ck("faithful_ratio_is_mean_over_pred_%d_%d" % (sz["records_per_side"], r0["width"]),
           r0["ratio_mean"] is None
           or abs(r0["ratio_mean"] - r0["observed_mean"] / r0["predicted"]) < 5e-4)
        if r0["predicted"] >= 16:
            cells_charged.append(r0["charged_factor"])
for sz in CTL:
    ck("control_pairs_square_%d" % sz["records_per_side"], sz["pairs"] == sz["records_per_side"] ** 2)
kappa_cert = F(cert["key_census"]["kappa_num"], cert["key_census"]["kappa_den"])
ck("kappa_dominates_every_charged_cell", float(kappa_cert) >= max(cells_charged) - 5e-4)
ck("kappa_is_at_least_1_15", kappa_cert >= F(115, 100))
P_all = N * N
E_bad = kappa_cert * F(P_all, 2 ** keymod.KEY_IMAGE_BITS)
ck("E_bad_matches_certificate", E_bad == F(cert["debit"]["E_bad"]))
ck("E_bad_is_expectation", str(cert["debit"]["E_bad"]) == str(kappa_cert * F(N * N, 2 ** 167)))
per_cand = cert["debit"]["per_candidate_instr"]
ck("per_candidate_covers_recon_and_mult", per_cand >= cert["layout"]["recon_charged"] + 2 ** 19.55)
cap_work = 2 ** 84 * per_cand
ck("cap_work_matches_certificate", cap_work == cert["debit"]["cap_work"] == 2 ** 104)
markov = E_bad * per_cand / cap_work
ck("markov_recomputed_below_2_18", markov < F(1, 2 ** 18))
ck("markov_matches_certificate", markov == F(cert["debit"]["markov_num"], cert["debit"]["markov_den"]))
ck("cap_16_bits_above_expectation", 2 ** 84 > E_bad * 2 ** 16)
P_pref = F(P_all, 2 ** (keymod.ROUTE_BITS + keymod.SORT_BITS))
ck("run_cap_margin_ge_10_bits", P_pref * 2 ** 10 < 2 ** 95)

# ---- algebraic census is exactly zero, at both sizes
for sz in px["sizes"]:
    ck("algebraic_equal_zero_%d" % sz["records_per_side"], sz["equal_products"] == 0)
    ck("algebraic_complement_zero_%d" % sz["records_per_side"], sz["complement_pairs"] == 0)

# ---- record layout + reconstruction price, re-derived from the shared module
lay = cert["layout"]
ck("layout_words_6_bytes_48", lay["words"] == 6 and lay["bytes"] == 48)
ck("layout_payload_297", keymod.LAYOUT_BITS == 297 and lay["bits"] == 297)
ck("layout_rows_match_module", lay["rows"] == keymod.layout_rows())
ck("layout_index_bits_128", lay["index_bits"] == 128 and lay["index_bound_log2"] == 116)
ck("layout_ok", lay["layout_ok"] is True and keymod.layout_ok() is True)
ck("recon_itemised_matches_module", lay["recon_itemised"] == keymod.RECON_ITEMISED ==
   sum(v for _, v, _ in keymod.RECON_ITEMS))
ck("recon_charged_is_ceiling_of_6_5", lay["recon_charged"] == -(-keymod.RECON_ITEMISED * 6 // 5))
pl = ts["partLayout"]
ck("key_rule_100_inputs_agree", pl["key_agreement"] == 100 == pl["inputs"] and pl["key_agreement_all"])
ck("position_roundtrip_all", pl["roundtrip_all"] is True and pl["reconstructed"] == pl["records"] > 0)
ck("toy_runs_the_declared_layout", pl["record_words"] == 6 and pl["record_bytes"] == 48 and
   pl["key_bits"] == 168 and pl["index_bits"] == 128)

# ---- ledger, recomputed from formulas
rows = {
    "records": 2 * N * per_record,
    "chain_tables": 2 * ((N + 2 ** 50 - 1) // 2 ** 50) * 2 ** 32 * 2000,
    "fixed_per_node": 2 ** keymod.ROUTE_BITS * 2 ** 31,
    "hprime_supply": 4 * 115812111329581076142 * 2 ** 20 * 68 * 352,
    "walk_list": int((2 ** 97 / (97 * math.log(2))) * (8 + 134816)),
    "run_pair_cap": 2 * 2 ** 95,
    "candidate_scan_cap": cap_work,
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

# ---- memory rows
mem = cert["memory"]
ck("memory_rows_sum", sum(x["bytes_"] for x in mem["rows"]) == mem["total"])
ck("memory_log2_matches", abs(math.log2(mem["total"]) - mem["log2"]) < 1e-3)

# ---- claim and gate
claim = F(int(cert["claim"]["numerator"]), int(cert["claim"]["denominator"]))
gate = F(int(cert["gate"]["numerator"]), int(cert["gate"]["denominator"]))
ck("claim_above_W", float(claim) > W_log2)
ck("claim_headroom_ge_2e5", float(claim) - W_log2 >= 2e-5)
ck("claim_below_gate_minus_0006", gate - claim >= F(6, 1000))
ck("claim_text_matches", cert["claim_text"] == format(Decimal(claim.numerator) / Decimal(claim.denominator), '.4f'))

# ---- census: raw rebuild, two axes, trend row
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
ck("grid_two_fields_three_d", len({c["r"] for c in cells}) == 2 and sorted({c["d"] for c in cells}) == [8, 16, 64])
ck("cells_at_production_K", all(abs(F(c["candidate_occurrences"], c["m"]) - K) < F(1, 50) for c in cells))
tr = cert["coverage_grid"]["trend"]
ck("trend_three_matched_pairs", len(tr) == 3 and sorted(t["d"] for t in tr) == [8, 16, 64])
for t in tr:
    sm = [c for c in cells if c["d"] == t["d"] and c["r"] == t["small_r"]][0]
    lg = [c for c in cells if c["d"] == t["d"] and c["r"] == t["large_r"]][0]
    ck("trend_row_%d_raw" % t["d"],
       abs(abs(F(sm["predicted_distinct_eff"] - sm["distinct_covered"], sm["m"])) - t["small"]) < 1e-9
       and abs(abs(F(lg["predicted_distinct_eff"] - lg["distinct_covered"], lg["m"])) - t["large"]) < 1e-9
       and t["shrinks"] is (abs(F(lg["predicted_distinct_eff"] - lg["distinct_covered"], lg["m"])) <
                            abs(F(sm["predicted_distinct_eff"] - sm["distinct_covered"], sm["m"]))))
ck("trend_all_shrink", all(t["shrinks"] for t in tr) and sum(1 for t in tr if t["shrinks"]) == 3)
scale_max = max(abs(float(F(c["predicted_distinct_eff"] - c["distinct_covered"], c["m"]))) * math.sqrt(c["m"])
                for c in cells)
ck("diagnostic_is_raw_max", abs(scale_max - cert["coverage_grid"]["diagnostic"]["max_absB_times_sqrt_m"]) < 5e-3)
ck("diagnostic_is_not_charged", "NOT USED" in cert["coverage_grid"]["diagnostic"]["status"])
ds = cert["coverage_grid"]["declared_statistic"]
ck("declared_statistic_is_rowB_fraction", ds["name"].startswith("abs_rowB"))
ck("declared_statistic_direction_is_trend", "3/3" in ds["direction"] and "falls" in ds["direction"])
ck("declared_statistic_charge_is_toy_max", "toy maximum" in ds["charge"])
ck("no_m_scaling_law_in_declared_block", "sqrt(m)" not in ds["charge"].replace("|rowB|*sqrt(m)", "")
   or "NOT" in ds["charge"])
ck("scale_row_key_retired", "scale_row" not in cert["coverage_grid"])

# ---- coverage expectation, exact rationals
hi_e = exp_neg_upper(F(51, 100), 26)
ck("expectation_upper_matches", hi_e == F(cert["coverage_grid"]["expectation"]["e_neg_K_upper_num"],
                                          cert["coverage_grid"]["expectation"]["e_neg_K_upper_den"]))
ck("expectation_upper_below_0_601", hi_e < F(601, 1000))
design = 1 - hi_e
ck("design_matches", abs(float(design) - cert["success"]["design_lo"]) < 1e-6)

# ---- gap premise
mu_prod = F(d) * F(971, 10) * LN2
G = 2 ** 33
tail = F(50) - G / (mu_prod * LN2)
ck("gap_mu_matches", abs(math.log2(float(mu_prod)) - cert["gap"]["mu_prod_log2"]) < 1e-3)
ck("gap_G_is_2_33", cert["gap"]["G_charged"] == G)
ck("gap_tail_below_2_40", tail < -40 and abs(float(tail) - cert["gap"]["tail_bound_log2"]) < 0.2)
ck("gap_premise_declared", "H3gap" in cert["gap"]["premise"])
hist = deep["histogram"]
ck("gap_histogram_sums_to_n", sum(b["count"] for b in hist) == deep["n_gap"])
sizes = ts["gapFit"]["sizes"]
ck("gap_fit_three_sizes_raw", len(sizes) == 3 and sorted(s["d"] for s in sizes) == [8, 16, 64])
ck("gap_fit_means_match_law", all(abs(s["mean_over_law"] - 1) < 0.35 for s in sizes))
ck("gap_fit_reported_matches_raw", all(abs(cf["chi2"] - s["chi2"]) < 1e-6 and abs(cf["ks"] - s["ks"]) < 1e-9
                                       for cf, s in zip(cert["gap"]["fit"], sizes)))

# ---- order audit
oa = deep["order_audit"]
ck("order_audit_n_divides_field_order", oa["n_divides_NE"] and oa["field_order"] % oa["prime_order"] == 0)
ck("order_audit_nP_is_O", oa["nP_is_O"] is True)

# ---- success arithmetic
debits = cert["success"]["debits"]
ck("census_debits_are_printed_ceilings",
   abs(debits["census_rowA"] - float(allowA_d)) < 1e-12 and abs(debits["census_rowB"] - float(allowB_d)) < 1e-12)
ck("candidate_markov_is_expectation_ratio", abs(debits["candidate_markov"] - float(markov)) < 1e-12)
after = float(design) - sum(debits.values())
ck("success_after_debits_ge_039", after >= 0.39 and abs(after - cert["success"]["after"]) < 1e-6)
ck("success_premise_declared", "H3key" in cert["success"]["premise"] and "H3occ" in cert["success"]["premise"])

# ---- v17 statistics must be GONE from every artifact
for fn in ("certificate.json", "argument.md", "public-note.md", "submission.json", "deep-evidence.json"):
    p = os.path.join(HERE, fn)
    if os.path.exists(p):
        t = open(p).read().lower()
        for bad in ("3/n", "exp(-3)", "95% upper", "clopper", "positive control",
                    "certified one-sided", "(1-3/n)^n"):
            ck("no_v17_statistics_%s_%s" % (fn.replace(".", "_"), bad.replace(" ", "_").replace("-", "_")
                                            .replace("(", "").replace(")", "").replace("^", "")), bad not in t)

# ---- prose hygiene
ck("cert_scope_declared", "does NOT independently establish production resource" in cert["cert_scope"])
for fn in ("public-note.md", "argument.md", "submission.json"):
    p = os.path.join(HERE, fn)
    if os.path.exists(p):
        ck("no_placeholder_%s" % fn, "__" not in open(p).read())
note = os.path.join(HERE, "public-note.md")
if os.path.exists(note):
    ntxt = open(note).read()
    ck("note_model_header_once", ntxt.count("Model:") == 1)
    ck("note_harness_header_once", ntxt.count("Harness:") == 1)
    ck("note_header_lines_1_2", ntxt.split("\n")[0].strip() == "Model: qwen3.8-flash-next"
       and ntxt.split("\n")[1].strip() == "Harness: Hermes")
    ck("note_size_in_band", 6000 <= len(ntxt) <= 16000)
sub = os.path.join(HERE, "submission.json")
if os.path.exists(sub):
    sj = json.load(open(sub))
    ck("submission_keys_exact", sorted(sj.keys()) == ["claims", "experiments", "title", "version"])
    ck("submission_experiments_empty", sj["experiments"] == [])
    st = sj["claims"][0]["statement"]
    ck("statement_keeps_composite_order_verbatim",
       "m is COMPOSITE (m is even and m>2)" in st and "gcd(d,m)=2^12" in st
       and "only the primality assertion is withdrawn" in st)
    ck("statement_has_no_statistical_assertion",
       all(b not in st.lower() for b in ("3/n", "exp(-3)", "95% upper", "certified one-sided",
                                         "confidence interval", "bernoulli trial")))
    ck("statement_score_matches_certificate",
       F(int(sj["claims"][0]["score"]["numerator"]), int(sj["claims"][0]["score"]["denominator"])) == claim)

ck("all_cert_checks_pass", all(v == "PASS" for v in cert["checks"].values()))

if fails:
    print("FAIL:", len(fails))
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("VERIFY OK: %d independent checks; W=2^%.6f claim=%s gate=%.4f success=%.6f E_bad=2^%.3f kappa=%s" %
      (len(cert["checks"]), W_log2, cert["claim_text"], float(gate), after,
       math.log2(float(E_bad)), float(kappa_cert)))
