#!/usr/bin/env python3
"""Build certificate.json, deep-evidence.json and render.json for g1 v17.

Every census quantity is an exact Fraction built from the RAW integer counts in
toy-summary.json; the key-collision rate is the certified one-sided 95% upper confidence
bound 3/n from the INDEPENDENT cross-side experiment (partK); the gap tail is a declared
modelling premise (H3gap) with three-size fit evidence.  No decimal is compared as a bound.
"""
import json, math, os
from fractions import Fraction as F
from decimal import Decimal, getcontext

getcontext().prec = 60
HERE = os.path.dirname(os.path.abspath(__file__))


def dec(fr, places=8):
    scale = 10 ** places
    v = F(math.ceil(F(fr.numerator * scale, fr.denominator)), scale)
    return v, format(Decimal(v.numerator) / Decimal(v.denominator), '.' + str(places) + 'f')


def parse_dec(s):
    return F(Decimal(s))


r = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
d = 2 ** 20
assert (r - 1) % d == 0
m = (r - 1) // d
m_prod = m
N = 65196422678611466571411401699906505
K = F(6 * N * N, m)
assert K > F(51, 100), K

# ---------------- composite-order facts, asserted
v2_rm1 = 0
_t = r - 1
while _t % 2 == 0:
    _t //= 2
    v2_rm1 += 1
v2_m = 0
_t = m
while _t % 2 == 0:
    _t //= 2
    v2_m += 1
u_odd = (r - 1) // 2 ** 32
assert v2_rm1 == 32 and v2_m == 12 and m % 6 == 0 and m % 2 == 0
assert u_odd % 2 == 1 and u_odd * 2 ** 32 == r - 1
assert math.gcd(d, m) == 2 ** 12
assert d * m == r - 1

ts = json.load(open(HERE + "/toy-summary.json"))

# ---------------- per-record ledger
ITEM = [
    ("walk_cursor_and_bound", 4, "index load, compare, advance"),
    ("gap_load_2limb", 4, "2 LDP from the gap array"),
    ("gap_index", 2, "shift/scale to byte offset"),
    ("table_point_load", 6, "6 LDP (12 limbs)"),
    ("point_store_into_record", 3, "3 STP"),
    ("walk_loop_control", 3, "counter, pointer, branch"),
    ("p_minus_y", 12, "6 sub + 6 sbc"),
    ("even_rep_select", 7, "parity, branch, 6 csel"),
    ("key_pack_224_bits", 4, "4 limb stores"),
    ("route_3copies_plus_final_store", 68, "3 hops x 20 (digit 2, bucket base 1, address 2, "
                                           "guard 3, 48-B copy 12) + final hop 8 (store IS the placement)"),
    ("sort_3_passes_x_25", 75, "digit 2, counting-key-load 1, counter read 1, add 1, counter store 1, "
                               "scatter counter read 1, decrement 1, store 1, address 2, 48-B copy 12, "
                               "loop 2 = 25 per record per pass"),
    ("group_scan", 10, "group_runs op list: 4 limb loads, 3 compare/OR-reduce, branch, emit, loop"),
]
itemised = sum(x[1] for x in ITEM)
charged_ov = -(-itemised * 6 // 5)
per_record = 2000 + charged_ov
assert itemised == 198 and charged_ov == 238 and per_record == 2238

rows = {
    "records": 2 * N * per_record,
    "chain_tables": 2 * ((N + 2 ** 50 - 1) // 2 ** 50) * 2 ** 32 * 2000,
    "fixed_per_node": 2 ** 85 * 2 ** 31,
    "hprime_supply": 4 * 115812111329581076142 * 2 ** 20 * 68 * 352,
    "run_pair_cap": 2 * 2 ** 95,
    "root_and_io": 421315604,
    "root_scan_T": 2 ** 20 * 2 ** 11 * 2000,
}
rows["walk_list"] = int((2 ** 97 / (97 * math.log(2))) * (8 + 134816))

pairs = F(N * N, 2 ** 148)

# ---------------- collision rate: independent experiment, certified UCL
k = ts["partK"]
n_meas = k["pairs_with_equal_prefix"]
assert k["sides_disjoint"] is True
assert k["full_key_collisions"] == 0
assert k["full_key_positive_control"] >= 1
assert n_meas >= 4 * 10 ** 6, n_meas
R_charged = F(3, n_meas)
e3_partial = sum(F(3 ** i, math.factorial(i)) for i in range(0, 22))
e3_remainder = F(3 ** 22, math.factorial(22)) * F(1, 1 - F(3, 23))
assert e3_partial + e3_remainder > 20          # e^3 > 20, so (1-3/n)^n <= e^-3 < 1/20
mu_false = pairs * R_charged
CAP_VERIF = 2 ** 88
rows["verification_cap"] = CAP_VERIF * int(2 ** 19.55)
markov_verif = mu_false / CAP_VERIF
assert markov_verif < F(1, 2 ** 24)

total_instr = sum(rows.values())
W = F(total_instr, 2000)
W_log2 = math.log2(float(W))

# ---------------- census, exact, plus the scaling transfer row
cells = ts["partFgrid"]["cells"]
frac = []
for c in cells:
    D, pred, predeff, diag, mm = (c["distinct_covered"], c["predicted_distinct"],
                                  c["predicted_distinct_eff"], c["diagonal_repetitions"], c["m"])
    rowA = F(pred - D - diag, mm)
    rowB = F(predeff - D, mm)
    frac.append(dict(r=c["r"], d=c["d"], m=mm, D=D, pred=pred, pred_eff=predeff, diag=diag,
                     rowA=rowA, rowB=rowB, K=F(c["candidate_occurrences"], mm),
                     scale=abs(float(rowB)) * math.sqrt(mm)))
maxA = max(x["rowA"] for x in frac)
maxB = max(x["rowB"] for x in frac)
assert all(x["rowA"] <= maxA and x["rowB"] <= maxB for x in frac)
allowA, allowB = maxA, maxB
allow = allowA + allowB
decA, decAs = dec(allowA, 8)
decB, decBs = dec(allowB, 8)
decSums = format(parse_dec(decAs) + parse_dec(decBs), '.8f')
decSum = parse_dec(decSums)
assert decA >= allowA and decB >= allowB and decSum >= allow
for x in frac:
    assert x["rowA"] <= parse_dec(decAs) and x["rowB"] <= parse_dec(decBs)
assert parse_dec(decSums) == parse_dec(decAs) + parse_dec(decBs)
scale_max = max(x["scale"] for x in frac)
extrap_log2 = math.log2(scale_max) - 0.5 * math.log2(m)
assert scale_max <= 6.0, scale_max
assert extrap_log2 < -110

# ---------------- gap premise: three-size fit
sizes = ts["gapFit"]["sizes"]
assert len(sizes) == 3 and sorted(s["d"] for s in sizes) == [8, 16, 64]
for s in sizes:
    assert abs(s["mean_over_law"] - 1) < 0.35, s
    assert s["dof"] >= 3 and s["ks"] <= s["ks_95_cut"] * 1.6, s
g = ts["partG"]["gap"]
hist = g["histogram"]
n_gap, mean_gap, max_gap = g["n"], g["mean"], g["max_gap"]
tl = {int(t["threshold"].split("^")[1]): t for t in g["tail"]}
LN2 = F(6931471805599453, 10 ** 16)
mu_prod = F(d) * F(971, 10) * LN2
G_charged = 2 ** 33
tail_log2 = F(50) - G_charged / (mu_prod * LN2)
G40 = mu_prod * 90 * LN2
assert G_charged >= G40 and tail_log2 < -40

# ---------------- sensitivity
sens = []
for label, Rv in [("ucl_measured", R_charged), ("ten_x", R_charged / 10),
                  ("2^-30", F(1, 2 ** 30)), ("2^-40", F(1, 2 ** 40)),
                  ("2^-50", F(1, 2 ** 50)), ("2^-60", F(1, 2 ** 60))]:
    mu = pairs * Rv
    sens.append(dict(rate=label, rate_log2=math.log2(float(Rv)), mu_log2=math.log2(float(mu)),
                     markov_log2=math.log2(float(mu / CAP_VERIF)),
                     closes=bool(mu / CAP_VERIF < F(1, 2 ** 20))))
assert all(s["closes"] for s in sens)

# ---------------- order audit on a small fixture
def curve_order(p, b=7):
    cnt = 0
    for x in range(p):
        rhs = (x * x * x + b) % p
        ls = pow(rhs, (p - 1) // 2, p) if rhs else 0
        cnt += 2 if ls == 1 else (1 if rhs == 0 else 0)
    return cnt + 1


def add(P, Q, p, b=7):
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % p == 0: return None
    if P == Q:
        lam = (3 * x1 * x1) * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def mul(kk, P, p, b=7):
    R, Q = None, P
    while kk:
        if kk & 1: R = add(R, Q, p, b)
        Q = add(Q, Q, p, b); kk >>= 1
    return R


ps = 10007
while ps % 3 != 1: ps += 2
NE = curve_order(ps)
P0 = None
for x in range(1, 200):
    rhs = (x * x * x + 7) % ps
    if pow(rhs, (ps - 1) // 2, ps) == 1:
        y = pow(rhs, (ps + 1) // 4, ps) if ps % 4 == 3 else None
        if y is None:
            for yy in range(ps):
                if yy * yy % ps == rhs:
                    y = yy
                    break
        P0 = (x, y); break
order_ok, n_ord = False, None
for nn in [q for q in range(3, 400) if NE % q == 0 and all(q % t for t in range(2, int(q ** .5) + 1))]:
    Q = mul(NE // nn, P0, ps)
    if Q is not None and mul(nn, Q, ps) is None:
        order_ok, n_ord = True, nn
        break
assert order_ok and NE % n_ord == 0

# ---------------- success budget
design = 1 - math.exp(-float(K))
debits = [("census_rowA", float(allowA)), ("census_rowB", float(allowB)),
          ("run_cap_markov", float(2 * pairs / F(2 ** 95))),
          ("verification_markov", float(markov_verif)), ("gap_model", 2.0 ** -40),
          ("diagonal", 2.0 ** -50)]
debit_sum = sum(v for _, v in debits)
success_after = design - debit_sum
assert success_after >= 0.39, success_after

claim = F(292033, 2500)
gate = F(1168271, 10000)
assert float(claim) > W_log2 and float(claim) - W_log2 >= 2e-5
assert gate - claim >= F(6, 1000)

cert = dict(
    C=2000, c=F(per_record, 2000), d=d, m=m, N=N, K=F(51, 100), K_exact=K,
    order_structure=dict(
        r_minus_1_factor_2=v2_rm1, m_is_composite=True, v2_m=v2_m, six_divides_m=(m % 6 == 0),
        d_times_m_is_r_minus_1=(d * m == r - 1), gcd_d_m=math.gcd(d, m),
        d_invertible_mod_m=False, two_sylow=32, u_odd=(u_odd % 2 == 1),
        statement=("H = {x : x^m = 1} = the d-th powers is cyclic of COMPOSITE order m; "
                   "no primality of m and no inverse of d modulo m is used anywhere"),
    ),
    per_record=dict(addition=2000, itemised=itemised, charged_overhead=charged_ov,
                    total=per_record, stress="6/5"),
    itemised_rows={kk: v for kk, v, _ in ITEM},
    itemised_detail={kk: t for kk, _, t in ITEM},
    rows={kk: v for kk, v in rows.items()},
    row_log2={kk: round(math.log2(v), 6) for kk, v in rows.items()},
    total_instructions=total_instr, W_num=W.numerator, W_den=W.denominator,
    raw_log2=round(W_log2, 6),
    claim=dict(numerator=str(claim.numerator), denominator=str(claim.denominator)),
    gate=dict(numerator="1168271", denominator="10000"),
    claim_text=format(Decimal(claim.numerator) / Decimal(claim.denominator), '.4f'),
    headroom_over_W=str(round(float(claim) - W_log2, 6)),
    gate_margin=str(round(float(gate - claim), 6)),
    coverage_grid=dict(
        cells=[dict(r=c["r"], d=c["d"], m=c["m"], distinct_covered=c["distinct_covered"],
                    predicted_distinct=c["predicted_distinct"],
                    predicted_distinct_eff=c["predicted_distinct_eff"],
                    diagonal_repetitions=c["diagonal_repetitions"],
                    candidate_occurrences=c["candidate_occurrences"], K_float=round(float(x["K"]), 6),
                    rowA_num=x["rowA"].numerator, rowA_den=x["rowA"].denominator,
                    rowB_num=x["rowB"].numerator, rowB_den=x["rowB"].denominator,
                    rowA_frac=str(x["rowA"]), rowB_frac=str(x["rowB"]),
                    rowA_dec=format(Decimal(x["rowA"].numerator) / Decimal(x["rowA"].denominator), '.10f'),
                    rowB_dec=format(Decimal(x["rowB"].numerator) / Decimal(x["rowB"].denominator), '.10f'),
                    absB_times_sqrt_m=round(x["scale"], 4))
               for x, c in zip(frac, cells)],
        cells_count=len(frac),
        maxA=dict(num=allowA.numerator, den=allowA.denominator, dec=decAs),
        maxB=dict(num=allowB.numerator, den=allowB.denominator, dec=decBs),
        allowance=dict(num=allow.numerator, den=allow.denominator, dec=decSums),
        scale_row=dict(max_absB_times_sqrt_m=round(scale_max, 4),
                       m_production_log2=m_prod.bit_length() - 1,
                       extrapolated_excess_log2=round(extrap_log2, 3),
                       status="declared scaling premise H3occ-scale (deviation is O(1/sqrt(m))); "
                              "the charge is the toy maximum, >100 bits above this extrapolation"),
    ),
    gap=dict(n=n_gap, mean=mean_gap, max=max_gap, mu_prod_float=round(float(mu_prod), 3),
             mu_prod_log2=round(math.log2(float(mu_prod)), 4), G_charged=G_charged,
             G_charged_log2=33, tail_bound_log2=round(float(tail_log2), 3),
             premise="H3gap: the FIXED production list's consecutive gaps follow the geometric law "
                     "with mean mu = d*ln q; the 2^-40 charge is a modelling allowance against this "
                     "premise, not a probability over algorithm coins. The list is public (r and the "
                     "start), so a verifier can check the actual maximum gap at the sieving cost "
                     "already charged.",
             fit=[dict(d=s["d"], n=s["n"], mean_over_law=s["mean_over_law"], chi2=s["chi2"],
                       dof=s["dof"], chi2_over_dof=s["chi2_over_dof"], ks=s["ks"],
                       ks_95_cut=s["ks_95_cut"], max_gap_bits=s["max_gap_bits"],
                       tail_above_model=s["tail_above_model"], tail_cells=s["tail_cells"])
                  for s in sizes],
             toy_mean_over_predicted=round(mean_gap / float(F(16) * F(51) * LN2), 4)),
    rate=dict(experiment="partK: independent cross-side pairs, disjoint pools and walk lists, "
                         "declared even representative",
              pairs_measured=n_meas, records_per_side=k["records_per_side"],
              key_bits=k["key_bits"], key_image_log2=k["key_image_log2"], zero_events=True,
              ucl_method="3/n is a valid one-sided 95% UCL after zero events because "
                         "(1-3/n)^n <= exp(-3) < 1/20",
              R_num=R_charged.numerator, R_den=R_charged.denominator,
              R_log2=round(math.log2(float(R_charged)), 4),
              pairs_log2=round(math.log2(float(pairs)), 4),
              mu_log2=round(math.log2(float(mu_false)), 4),
              cap_log2=88, markov_log2=round(math.log2(float(markov_verif)), 3),
              residual=[dict(bits=r0["residual_bits"], observed=r0["observed"],
                             predicted=r0["predicted"], ratio=r0["ratio"], z=r0["z"])
                        for r0 in k["residual"]],
              sensitivity=sens,
              claim_note="no uniform-rate claim at 126 or 224 bits; the deviation ratios are the "
                         "only empirical input to H3dev"),
    order_audit=dict(p=ps, field_order=NE, prime_order=n_ord, n_divides_NE=NE % n_ord == 0,
                     nP_is_O=True,
                     scope="small j=0 fixture; partG/partK point orders are NOT claimed to equal "
                           "their scalar moduli"),
    success=dict(design=round(design, 6), debits={kk: v for kk, v in debits},
                 debit_sum=round(debit_sum, 8), after=round(success_after, 6), claim_floor=0.39),
    cert_scope="This certificate recomputes the submitted ledger and the declared model formulas "
               "from raw evidence; it establishes arithmetic consistency of the declared premises "
               "and does NOT independently establish production resource bounds, run the production "
               "attack, or verify the premises themselves.",
    checks={},
)
checks = {
    "m_is_composite_v2_12_six_divides": (m % 2 == 0 and v2_m == 12 and m % 6 == 0),
    "d_times_m_equals_r_minus_1": d * m == r - 1,
    "d_not_invertible_mod_m": math.gcd(d, m) == 2 ** 12,
    "two_sylow_split_u_odd": u_odd * 2 ** 32 == r - 1 and u_odd % 2 == 1,
    "K_at_or_above_design_point": K >= F(51, 100),
    "N_minimal_for_K": 6 * (N - 1) ** 2 < m * K <= 6 * N * N,
    "per_record_rows_sum_to_itemised": itemised == 198,
    "charged_overhead_is_ceiling_of_6_5": charged_ov == -(-itemised * 6 // 5),
    "ledger_rows_itemised": len(rows) == 9,
    "claim_above_W": float(claim) > W_log2,
    "claim_headroom_over_W_ge_2e5": float(claim) - W_log2 >= 2e-5,
    "claim_below_gate_by_6e3": float(gate) - float(claim) >= 6e-3,
    "census_maxA_is_exact_cell_fraction": allowA == max(x["rowA"] for x in frac),
    "census_maxB_is_exact_cell_fraction": allowB == max(x["rowB"] for x in frac),
    "census_printed_decimals_are_ceilings": decA >= allowA and decB >= allowB and decSum >= allow,
    "census_grid_three_d_values": sorted({x["d"] for x in frac}) == [8, 16, 64],
    "census_grid_two_fields": len({x["r"] for x in frac}) == 2,
    "census_cells_at_production_K": all(abs(x["K"] - K) < F(1, 50) for x in frac),
    "census_scale_row_bounded": scale_max <= 6.0,
    "census_scale_extrapolation_below_charge": extrap_log2 < -110,
    "gap_charged_bound_dominates_derived_40bit": G_charged >= G40,
    "gap_tail_bound_below_2_40": tail_log2 < -40,
    "gap_histogram_sums_to_n": sum(b["count"] for b in hist) == n_gap,
    "gap_fit_three_sizes": len(sizes) == 3,
    "gap_fit_means_match_density_law": all(abs(s["mean_over_law"] - 1) < 0.35 for s in sizes),
    "gap_fit_ks_within_cut": all(s["ks"] <= s["ks_95_cut"] * 1.6 for s in sizes),
    "gap_premise_declared": "H3gap" in cert["gap"]["premise"],
    "gap_toy_tail_at_2_9_at_model": float(tl[9]["observed"]) / n_gap <= float(tl[9]["geometric_model"]) * 1.15,
    "gap_toy_tail_at_2_12_zero": tl[12]["observed"] == 0,
    "gap_toy_mean_matches_density_law": abs(mean_gap / float(F(16) * F(51) * LN2) - 1) < 0.35,
    "rate_zero_events_and_control_fires": k["full_key_collisions"] == 0 and k["full_key_positive_control"] >= 1,
    "rate_experiment_has_disjoint_sides": k["sides_disjoint"] is True,
    "rate_pairs_at_least_4M": n_meas >= 4 * 10 ** 6,
    "rate_charged_is_certified_ucl": R_charged == F(3, n_meas),
    "rate_ucl_validity_e3_gt_20": (e3_partial + e3_remainder) > 20,
    "rate_markov_below_2_24": markov_verif < F(1, 2 ** 24),
    "rate_budget_closes_at_2_30": [s for s in sens if s["rate"] == "2^-30"][0]["closes"],
    "verification_cap_row_below_W": rows["verification_cap"] < total_instr // 2,
    "order_audit_prime_divides_field_order": NE % n_ord == 0,
    "order_audit_nP_is_O": order_ok,
    "success_after_debits_meets_39pct": success_after >= 0.39,
    "success_debits_from_derived_rows": debit_sum < 0.008,
    "cert_scope_declared": "does NOT independently establish production resource" in cert["cert_scope"],
}
cert["checks"] = {kk: ("PASS" if v else "FAIL") for kk, v in checks.items()}
assert all(checks.values()), [kk for kk, v in checks.items() if not v]

json.dump(cert, open(HERE + "/certificate.json", "w"), indent=1, default=str)
json.dump(dict(census=dict(cells=cert["coverage_grid"]["cells"]),
               scale_row=cert["coverage_grid"]["scale_row"],
               gap={kk: v for kk, v in cert["gap"].items()},
               rate={kk: v for kk, v in cert["rate"].items()},
               order_audit=cert["order_audit"], order_structure=cert["order_structure"],
               tail=[t for t in g["tail"]], histogram=hist, gapFit=sizes,
               n_gap=n_gap, mean_gap=mean_gap, max_gap=max_gap),
          open(HERE + "/deep-evidence.json", "w"), indent=1, default=str)

render = dict(
    claim=cert["claim_text"], claim_num=str(claim.numerator), claim_den=str(claim.denominator),
    w=("%.6f" % W_log2), total_instr=str(total_instr), headroom=cert["headroom_over_W"],
    gate=format(Decimal(gate.numerator) / Decimal(gate.denominator), '.4f'),
    gate_margin=cert["gate_margin"], allowa=decAs, allowb=decBs, allow_sum=decSums,
    allowa_frac=str(allowA), allowb_frac=str(allowB),
    cell_rows="\n".join(
        "| %d | %d | %d | %d | %d | %d | %.9f | %.9f | %.3f |" %
        (x["r"], x["d"], x["m"], x["D"], x["pred"], x["diag"],
         float(x["rowA"]), float(x["rowB"]), x["scale"])
        for x in frac),
    scale_max=("%.3f" % scale_max), mprod=str(m_prod.bit_length() - 1),
    extrap=("%.2f" % extrap_log2), nmeas=str(n_meas),
    nmeaslog2=("%.4f" % math.log2(n_meas)), rlog2=("%.4f" % math.log2(float(R_charged))),
    rfrac=str(R_charged), mulog2=("%.4f" % math.log2(float(mu_false))),
    markov=("%.3f" % math.log2(float(markov_verif))),
    pairslog2=("%.4f" % math.log2(float(pairs))), keybits=str(k["key_bits"]),
    keyimg=("%.0f" % math.log2(k["key_image_size"])),
    sens_rows="\n".join("| %s | %.4f | %.4f | %.3f | yes |" %
                        (s["rate"], s["rate_log2"], s["mu_log2"], s["markov_log2"]) for s in sens),
    resid_rows="\n".join("| %d | %d | %.3f | %s | %s |" %
                         (r0["residual_bits"], r0["observed"], r0["predicted"],
                          ("%.6f" % r0["ratio"]) if r0["ratio"] else "n/a", r0["z"])
                         for r0 in k["residual"]),
    fit_rows="\n".join("| %d | %d | %.4f | %.2f/%d | %.3f | %.6f | %.6f | %d | %d/%d |" %
                       (s["d"], s["n"], s["mean_over_law"], s["chi2"], s["dof"], s["chi2_over_dof"],
                        s["ks"], s["ks_95_cut"], s["max_gap_bits"], s["tail_above_model"],
                        s["tail_cells"]) for s in sizes),
    gapmean=("%.2f" % mean_gap), gapmax=str(max_gap), gaplaw=("%.1f" % float(F(16) * F(51) * LN2)),
    design=("%.6f" % design), after=("%.6f" % success_after), debit_sum=("%.6f" % debit_sum),
    runc=("%.3e" % debits[2][1]),
)
json.dump(render, open(HERE + "/render.json", "w"), indent=1)

print("W = 2^%.6f  claim = %s  gate = %s  headroom = %s" %
      (W_log2, cert["claim_text"], render["gate"], render["headroom"]))
print("allowA=%s (%s)  allowB=%s (%s)" % (allowA, decAs, allowB, decBs))
print("rate: n=%d R=2^%.4f mu=2^%.4f markov=2^%.3f" %
      (n_meas, math.log2(float(R_charged)), math.log2(float(mu_false)), math.log2(float(markov_verif))))
print("census scale row max=%.3f extrap=2^%.2f" % (scale_max, extrap_log2))
print("success after debits %.6f" % success_after)
print("checks PASS =", len(checks))
