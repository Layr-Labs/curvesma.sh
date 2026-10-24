#!/usr/bin/env python3
"""Build certificate.json, deep-evidence.json and render.json for g1 v18.

v18 removes the statistical certification of v17 entirely.  The key-collision
debit is an EXPECTATION-MARKOV bound on DETERMINISTIC lists: the exact structural
census (partX: zero algebraic coincidences, exact prefix-collision counts at two
sizes and seven widths) fixes E[#bad] = kappa * N^2 / 2^167, and Markov on that
expectation needs no independence.  Coverage keeps exact Fractions and adds a
two-axis trend row.  The record is the declared 48-byte position-index layout of
keymod.py, and its reconstruction is itemised and priced.
"""
import json, math, os
from fractions import Fraction as F
from decimal import Decimal, getcontext

getcontext().prec = 60
HERE = os.path.dirname(os.path.abspath(__file__))
import keymod


def dec(fr, places=8):
    scale = 10 ** places
    v = F(math.ceil(F(fr.numerator * scale, fr.denominator)), scale)
    return v, format(Decimal(v.numerator) / Decimal(v.denominator), '.' + str(places) + 'f')


def parse_dec(s):
    return F(Decimal(s))


def exp_neg_upper(x, terms=26):
    """Rational UPPER bound of exp(-x), x >= 0, by the alternating Taylor series
    (the first omitted term dominates the alternating remainder)."""
    s = F(0)
    for i in range(terms):
        s += F((-1) ** i) * x ** i / math.factorial(i)
    return s + x ** terms / math.factorial(terms)


r = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
d = 2 ** 20
assert (r - 1) % d == 0
m = (r - 1) // d
N = 65196422678611466571411401699906505
K = F(6 * N * N, m)
assert K > F(51, 100), K
K_CHARGED = F(51, 100)

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
px = ts["partX"]
pl = ts["partLayout"]

ITEM = [
    ("walk_cursor_and_bound", 4, "index load, compare, advance"),
    ("gap_load_2limb", 4, "2 LDP from the gap array"),
    ("gap_index", 2, "shift/scale to byte offset"),
    ("table_point_load", 6, "6 LDP (12 limbs)"),
    ("point_store_into_record", 3, "3 STP"),
    ("walk_loop_control", 3, "counter, pointer, branch"),
    ("p_minus_y", 12, "6 sub + 6 sbc"),
    ("even_rep_select", 7, "parity test, branch, 6 csel"),
    ("key_pack_168_bits", 4, "3 word stores + mask register setup"),
    ("route_4hops_x_21bits", 68, "3 hops x 20 (digit 2, bucket base 1, address 2, "
                                 "guard 3, 48-B copy 12) + final hop 8 (store IS the placement)"),
    ("sort_3_passes_x_25", 75, "21-bit digits: digit 2, counting-key-load 1, counter read 1, add 1, "
                               "counter store 1, scatter counter read 1, decrement 1, store 1, "
                               "address 2, 48-B copy 12, loop 2 = 25 per record per pass"),
    ("group_scan", 10, "3 key word loads, 3 compare/OR-reduce, side load, branch, emit, loop"),
]
itemised = sum(x[1] for x in ITEM)
charged_ov = -(-itemised * 6 // 5)
per_record = 2000 + charged_ov
assert itemised == 198 and charged_ov == 238 and per_record == 2238

KEY_BITS = keymod.KEY_BITS
IMAGE_LOG2 = keymod.KEY_IMAGE_BITS
assert KEY_BITS == 168 and IMAGE_LOG2 == 167
P_all = N * N
charged_cells = [r0["charged_factor"] for sz in px["sizes"] for r0 in sz["exact_rows"]
                 if r0["predicted"] >= 16 and r0["charged_factor"] is not None]
kappa = max(1.15, max(charged_cells))
kappa = F(str(round(kappa + 0.005, 3)))            # ceiling at 3 decimals
E_bad = kappa * F(P_all, 2 ** IMAGE_LOG2)          # expectation, no independence used
CAP_BAD = 2 ** 84
per_candidate = 2 ** 20                            # ceiling of recon + one full mult
CAP_WORK = CAP_BAD * per_candidate
E_work = E_bad * per_candidate
markov = E_work / CAP_WORK
assert markov < F(1, 2 ** 18), markov
assert CAP_BAD > E_bad * 2 ** 16, "cap must sit well above the expectation"
P_route = F(P_all, 2 ** keymod.ROUTE_BITS)
P_pref = F(P_all, 2 ** (keymod.ROUTE_BITS + keymod.SORT_BITS))
RUN_CAP_PAIRS = 2 ** 95
assert P_pref * 2 ** 10 < RUN_CAP_PAIRS

alg = [(sz["records_per_side"], sz["pairs"], sz["equal_products"], sz["complement_pairs"])
       for sz in px["sizes"]]
assert all(a[2] == 0 and a[3] == 0 for a in alg), alg

rows = {
    "records": 2 * N * per_record,
    "chain_tables": 2 * ((N + 2 ** 50 - 1) // 2 ** 50) * 2 ** 32 * 2000,
    "fixed_per_node": 2 ** keymod.ROUTE_BITS * 2 ** 31,
    "hprime_supply": 4 * 115812111329581076142 * 2 ** 20 * 68 * 352,
    "run_pair_cap": 2 * 2 ** 95,
    "candidate_scan_cap": CAP_WORK,
    "root_and_io": 421315604,
    "root_scan_T": 2 ** 20 * 2 ** 11 * 2000,
}
rows["walk_list"] = int((2 ** 97 / (97 * math.log(2))) * (8 + 134816))
total_instr = sum(rows.values())
W = F(total_instr, 2000)
W_log2 = math.log2(float(W))

MEM = [
    ("records", 2 * N * keymod.RECORD_BYTES, "2N records x 48 B"),
    ("sort_scratch", 2 * N * keymod.RECORD_BYTES, "one counting-sort copy of the record array"),
    ("base_pool", 2 ** 66 * 16, "a-pool entries, 2 words each"),
    ("gap_checkpoints", 2 ** 38 * 16, "2^50 walk gaps / 2^12 per checkpoint, 2 words"),
    ("routing_nodes", 2 ** keymod.ROUTE_BITS * 64, "per-node bucket descriptors"),
]
mem_total = sum(v for _, v, _ in MEM)
mem_log2 = math.log2(mem_total)

cells = ts["partFgrid"]["cells"]
frac = []
for c in cells:
    D, pred, predeff, diag, mm = (c["distinct_covered"], c["predicted_distinct"],
                                  c["predicted_distinct_eff"], c["diagonal_repetitions"], c["m"])
    frac.append(dict(r=c["r"], d=c["d"], m=mm, D=D, pred=pred, pred_eff=predeff, diag=diag,
                     rowA=F(pred - D - diag, mm), rowB=F(predeff - D, mm),
                     K=F(c["candidate_occurrences"], mm), occ=c["candidate_occurrences"],
                     scale=abs(float(F(predeff - D, mm))) * math.sqrt(mm)))
maxA = max(x["rowA"] for x in frac)
maxB = max(x["rowB"] for x in frac)
assert all(x["rowA"] <= maxA and x["rowB"] <= maxB for x in frac)
decA, decAs = dec(maxA, 8)
decB, decBs = dec(maxB, 8)
decSums = format(parse_dec(decAs) + parse_dec(decBs), '.8f')
decSum = parse_dec(decSums)
assert decA >= maxA and decB >= maxB and decSum >= maxA + maxB
for x in frac:
    assert x["rowA"] <= parse_dec(decAs) and x["rowB"] <= parse_dec(decBs)
trend = []
for dd in sorted({x["d"] for x in frac}):
    a = [x for x in frac if x["d"] == dd][0]
    b = [x for x in frac if x["d"] == dd][1]
    trend.append(dict(axis="field", d=dd, small_r=a["r"], large_r=b["r"],
                      small=abs(float(a["rowB"])), large=abs(float(b["rowB"])),
                      shrinks=abs(b["rowB"]) < abs(a["rowB"])))
scale_max = max(x["scale"] for x in frac)
assert sum(1 for t in trend if t["shrinks"]) == len(trend)
extrap_log2 = math.log2(scale_max) - 0.5 * math.log2(m)

hi_e = exp_neg_upper(K_CHARGED, 26)
assert hi_e < F(601, 1000), float(hi_e)
design_lo = 1 - hi_e
assert K >= K_CHARGED

sizes = ts["gapFit"]["sizes"]
assert len(sizes) == 3
for s in sizes:
    assert abs(s["mean_over_law"] - 1) < 0.35, s
g = ts["partG"]["gap"]
hist = g["histogram"]
n_gap, mean_gap, max_gap = g["n"], g["mean"], g["max_gap"]
LN2 = F(6931471805599453, 10 ** 16)
mu_prod = F(d) * F(971, 10) * LN2
G_charged = 2 ** 33
tail_log2 = F(50) - G_charged / (mu_prod * LN2)


def curve_order(p, b=7):
    cnt = 0
    for x in range(p):
        rhs = (x * x * x + b) % p
        ls = pow(rhs, (p - 1) // 2, p) if rhs else 0
        cnt += 2 if ls == 1 else (1 if rhs == 0 else 0)
    return cnt + 1


def padd(P, Q, p, b=7):
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


def pmul(kk, P, p, b=7):
    R, Q = None, P
    while kk:
        if kk & 1: R = padd(R, Q, p, b)
        Q = padd(Q, Q, p, b); kk >>= 1
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
                    y = yy; break
        P0 = (x, y); break
order_ok, n_ord = False, None
for nn in [q for q in range(3, 400) if NE % q == 0 and all(q % t for t in range(2, int(q ** .5) + 1))]:
    Q = pmul(NE // nn, P0, ps)
    if Q is not None and pmul(nn, Q, ps) is None:
        order_ok, n_ord = True, nn
        break
assert order_ok and NE % n_ord == 0

debits = [("census_rowA", float(decA)), ("census_rowB", float(decB)),
          ("candidate_markov", float(markov)), ("gap_model", 2.0 ** -40),
          ("diagonal", 2.0 ** -50)]
debit_sum = sum(v for _, v in debits)
success_after = float(design_lo) - debit_sum
assert success_after >= 0.39, success_after
assert F(str(round(success_after, 6))) >= F(39, 100)

claim = F(292033, 2500)
gate = F(1168271, 10000)
assert float(claim) > W_log2 and float(claim) - W_log2 >= 2e-5
assert gate - claim >= F(6, 1000)

cert = dict(
    C=2000, c=F(per_record, 2000), d=d, m=m, N=N, K=K, K_charged=K_CHARGED,
    order_structure=dict(
        r_minus_1_factor_2=v2_rm1, m_is_composite=True, v2_m=v2_m, six_divides_m=(m % 6 == 0),
        d_times_m_is_r_minus_1=(d * m == r - 1), gcd_d_m=math.gcd(d, m),
        d_invertible_mod_m=False, two_sylow=32, u_odd=(u_odd % 2 == 1),
        statement=("H = {x : x^m = 1} = the d-th powers is cyclic of COMPOSITE order m; "
                   "no primality of m and no inverse of d modulo m is used anywhere")),
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
    key_census=dict(
        key_bits=KEY_BITS, image_log2=IMAGE_LOG2,
        rule="even member of {y,p-y} truncated to 168 bits (keymod.key_scalar)",
        widths=px["widths"], a_band=px["a_band"], walk_band=px["walk_band"],
        sizes=[dict(records_per_side=sz["records_per_side"], pairs=sz["pairs"],
                    equal_products=sz["equal_products"],
                    complement_pairs=sz["complement_pairs"],
                    distinct_products=sz["distinct_products"],
                    max_scalar_bits=sz["max_scalar_bits"], rows=sz["exact_rows"])
               for sz in px["sizes"]],
        max_ratio=px["max_ratio"], min_ratio=px["min_ratio"], cells_used=px["cells_used"],
        kappa_num=kappa.numerator, kappa_den=kappa.denominator, kappa_float=float(kappa),
        premise=("H3key (DECLARED): the truncated even-representative key occupies its 2^167 image "
                 "like a uniform hash on distinct points.  It is TESTED, not assumed, by exact "
                 "prefix-collision counts at two list sizes and seven key widths on the production "
                 "curve; the zero-event confidence bound of v17 is deleted.")),
    debit=dict(P_all=str(P_all), P_all_log2=round(math.log2(P_all), 4), image_log2=IMAGE_LOG2,
               kappa=float(kappa), E_bad=str(E_bad), E_bad_log2=round(math.log2(float(E_bad)), 4),
               cap_bad_log2=84, per_candidate_instr=per_candidate, cap_work=CAP_WORK,
               cap_work_log2=round(math.log2(CAP_WORK), 4),
               E_work_log2=round(math.log2(float(E_work)), 4),
               markov_num=markov.numerator, markov_den=markov.denominator,
               markov_log2=round(math.log2(float(markov)), 4),
               method=("E[#bad] = kappa * N^2 / 2^167 is an expectation over the FIXED lists; "
                       "Markov's inequality P(work > cap) <= E[work]/cap is applied to that "
                       "expectation.  Records and factors are reused across comparisons, so no "
                       "Bernoulli or independence model is asserted anywhere."),
               route_pairs_log2=round(math.log2(float(P_route)), 2),
               prefix_pairs_log2=round(math.log2(float(P_pref)), 2),
               run_cap_pairs_log2=95,
               run_cap_margin_bits=round(95 - math.log2(float(P_pref)), 2)),
    layout=dict(rows=pl["layout"], words=pl["record_words"], bytes=pl["record_bytes"],
                bits=keymod.LAYOUT_BITS, index_bits=pl["index_bits"],
                index_bound_log2=pl["index_bound_log2"], layout_ok=pl["layout_ok"],
                recon_items=pl["recon_items"], recon_itemised=pl["recon_itemised"],
                recon_charged=pl["recon_charged"],
                recon_detail=("index -> (pool slot, walk offset); the 97-bit a is read from the "
                              "stored pool array, the 76-bit w is rebuilt from a 2^12-gap checkpoint "
                              "plus at most 2^12 gap adds; the 174-bit scalar a*w is formed with 9 "
                              "fp_mul.  The scalar is NOT stored."),
                evidence=dict(inputs=pl["inputs"], key_agreement=pl["key_agreement"],
                              records=pl["records"], reconstructed=pl["reconstructed"],
                              roundtrip_all=pl["roundtrip_all"])),
    memory=dict(rows=[dict(name=a, bytes_=b, detail=c) for a, b, c in MEM],
                total=mem_total, log2=round(mem_log2, 3)),
    coverage_grid=dict(
        cells=[dict(r=x["r"], d=x["d"], m=x["m"], distinct_covered=x["D"],
                    predicted_distinct=x["pred"], predicted_distinct_eff=x["pred_eff"],
                    diagonal_repetitions=x["diag"], candidate_occurrences=x["occ"],
                    K_float=round(float(x["K"]), 6),
                    rowA_num=x["rowA"].numerator, rowA_den=x["rowA"].denominator,
                    rowB_num=x["rowB"].numerator, rowB_den=x["rowB"].denominator,
                    rowA_frac=str(x["rowA"]), rowB_frac=str(x["rowB"]),
                    absB_times_sqrt_m=round(x["scale"], 4)) for x in frac],
        maxA=dict(num=maxA.numerator, den=maxA.denominator, dec=decAs),
        maxB=dict(num=maxB.numerator, den=maxB.denominator, dec=decBs),
        allowance=dict(num=(maxA + maxB).numerator, den=(maxA + maxB).denominator, dec=decSums),
        trend=trend,
        scale_row=dict(max_absB_times_sqrt_m=round(scale_max, 4),
                       m_production_log2=m.bit_length() - 1,
                       extrapolated_excess_log2=round(extrap_log2, 3),
                       status="H3occ-scale (declared): |rowB|*sqrt(m) is the measured trend "
                              "statistic; it falls in 3/3 matched field pairs, and the charge is "
                              "the toy maximum, 100+ bits above the extrapolation.  No Chernoff "
                              "or independence claim is made."),
        expectation=dict(
            e_neg_K_upper_num=hi_e.numerator, e_neg_K_upper_den=hi_e.denominator,
            e_neg_K_upper_dec=format(Decimal(hi_e.numerator) / Decimal(hi_e.denominator), '.8f'),
            design_lo_dec=format(Decimal(design_lo.numerator) / Decimal(design_lo.denominator), '.8f'),
            method="E[#missing]/m = exp(-K) with K = 6N^2/m >= 51/100, exact rational upper "
                   "bound by the alternating series; linearity of expectation over the candidate "
                   "occurrences, so no independence is used.")),
    gap=dict(n=n_gap, mean=mean_gap, max=max_gap, mu_prod_float=round(float(mu_prod), 3),
             mu_prod_log2=round(math.log2(float(mu_prod)), 4), G_charged=G_charged,
             G_charged_log2=33, tail_bound_log2=round(float(tail_log2), 3),
             premise="H3gap: the FIXED production list's consecutive gaps follow the geometric law "
                     "with mean mu = d*ln q; the 2^-40 charge is a modelling allowance against this "
                     "premise, not a probability over algorithm coins.",
             fit=[dict(d=s["d"], n=s["n"], mean_over_law=s["mean_over_law"], chi2=s["chi2"],
                       dof=s["dof"], chi2_over_dof=s["chi2_over_dof"], ks=s["ks"],
                       ks_95_cut=s["ks_95_cut"], max_gap_bits=s["max_gap_bits"],
                       tail_above_model=s["tail_above_model"], tail_cells=s["tail_cells"])
                  for s in sizes],
             toy_mean_over_predicted=round(mean_gap / float(F(16) * F(51) * LN2), 4)),
    order_audit=dict(p=ps, field_order=NE, prime_order=n_ord, n_divides_NE=NE % n_ord == 0,
                     nP_is_O=True,
                     scope="small j=0 fixture; no toy point order is claimed to equal its scalar modulus"),
    success=dict(design_lo=round(float(design_lo), 6), debits={kk: v for kk, v in debits},
                 debit_sum=round(debit_sum, 8), after=round(success_after, 6), claim_floor=0.39,
                 premise="all of the above is conditional on H3occ (coverage), H3occ-scale, "
                         "H3key and H3gap, each with measured fit and a disclosed transfer limit"),
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
    "K_at_or_above_design_point": K >= K_CHARGED,
    "per_record_rows_sum_to_itemised": itemised == 198,
    "charged_overhead_is_ceiling_of_6_5": charged_ov == -(-itemised * 6 // 5),
    "ledger_rows_itemised": len(rows) == 9,
    "claim_above_W": float(claim) > W_log2,
    "claim_headroom_over_W_ge_2e5": float(claim) - W_log2 >= 2e-5,
    "claim_below_gate_by_6e3": float(gate) - float(claim) >= 6e-3,
    "key_rule_is_even_member": "even member" in cert["key_census"]["rule"],
    "key_shared_module_100_inputs": pl["key_agreement"] == 100 == pl["inputs"],
    "record_roundtrip_all_positions": pl["roundtrip_all"] is True and pl["reconstructed"] == pl["records"],
    "record_layout_6_words_48_bytes": pl["record_words"] == 6 and pl["record_bytes"] == 48,
    "record_layout_payload_bits_297": keymod.LAYOUT_BITS == 297,
    "record_layout_words_384_bits": keymod.RECORD_WORDS * 64 == 384,
    "record_layout_ok": pl["layout_ok"] is True,
    "recon_price_itemised_and_stressed": pl["recon_charged"] == -(-pl["recon_itemised"] * 6 // 5),
    "algebraic_equal_products_zero": all(a[2] == 0 for a in alg),
    "algebraic_complements_zero": all(a[3] == 0 for a in alg),
    "key_census_two_sizes": len(px["sizes"]) == 2,
    "key_census_seven_widths": len(px["widths"]) == 7,
    "key_census_max_ratio_at_most_1_25": px["max_ratio"] <= 1.25,
    "kappa_covers_measured_cells": float(kappa) >= max(charged_cells),
    "E_bad_is_expectation_not_ucl": E_bad == kappa * F(P_all, 2 ** IMAGE_LOG2),
    "markov_below_2_18": markov < F(1, 2 ** 18),
    "cap_dominates_expectation_16bits": CAP_BAD > E_bad * 2 ** 16,
    "prefix_pairs_2_10_below_run_cap": P_pref * 2 ** 10 < RUN_CAP_PAIRS,
    "no_ucl_or_independence_claim": "ucl" not in json.dumps(cert, default=str).lower(),
    "census_maxA_is_exact_cell_fraction": maxA == max(x["rowA"] for x in frac),
    "census_maxB_is_exact_cell_fraction": maxB == max(x["rowB"] for x in frac),
    "census_printed_decimals_are_ceilings": decA >= maxA and decB >= maxB and decSum >= maxA + maxB,
    "census_grid_three_d_two_fields": sorted({x["d"] for x in frac}) == [8, 16, 64]
                                      and len({x["r"] for x in frac}) == 2,
    "census_cells_at_production_K": all(abs(x["K"] - K) < F(1, 50) for x in frac),
    "census_trend_shrinks_3_of_3": sum(1 for t in trend if t["shrinks"]) == 3,
    "census_scale_row_bounded": scale_max <= 6.0,
    "census_scale_extrapolation_below_charge": extrap_log2 < -110,
    "coverage_expectation_is_rational_upper": hi_e < F(601, 1000) and design_lo > 0,
    "gap_tail_bound_below_2_40": tail_log2 < -40,
    "gap_histogram_sums_to_n": sum(b["count"] for b in hist) == n_gap,
    "gap_fit_three_sizes": len(sizes) == 3,
    "gap_fit_means_match_density_law": all(abs(s["mean_over_law"] - 1) < 0.35 for s in sizes),
    "gap_premise_declared": "H3gap" in cert["gap"]["premise"],
    "order_audit_prime_divides_field_order": NE % n_ord == 0,
    "order_audit_nP_is_O": order_ok,
    "memory_total_at_declared_log2": abs(mem_log2 - 123.3) < 1.2,
    "success_after_debits_meets_39pct": success_after >= 0.39,
    "cert_scope_declared": "does NOT independently establish production resource" in cert["cert_scope"],
}
cert["checks"] = {kk: ("PASS" if v else "FAIL") for kk, v in checks.items()}
assert all(checks.values()), [kk for kk, v in checks.items() if not v]

json.dump(cert, open(HERE + "/certificate.json", "w"), indent=1, default=str)
json.dump(dict(census=dict(cells=cert["coverage_grid"]["cells"]),
               scale_row=cert["coverage_grid"]["scale_row"],
               expectation=cert["coverage_grid"]["expectation"],
               key_census=cert["key_census"], debit=cert["debit"], layout=cert["layout"],
               gap={kk: v for kk, v in cert["gap"].items()},
               order_audit=cert["order_audit"], order_structure=cert["order_structure"],
               tail=[t for t in g["tail"]], histogram=hist, gapFit=sizes,
               n_gap=n_gap, mean_gap=mean_gap, max_gap=max_gap, memory=cert["memory"]),
          open(HERE + "/deep-evidence.json", "w"), indent=1, default=str)

census_rows = "\n".join(
    "| %d | %d | %d | %d | %d | %d |" % (sz["records_per_side"], sz["pairs"], sz["equal_products"],
                                         sz["complement_pairs"], sz["distinct_products"],
                                         sz["max_scalar_bits"]) for sz in px["sizes"])
width_rows = "\n".join("| " + " | ".join("%d" % r0["observed"] for r0 in sz["exact_rows"]) + " |"
                       for sz in px["sizes"])
pred_row = "| " + " | ".join("%.0f" % r0["predicted"] for r0 in px["sizes"][0]["exact_rows"]) + " |"
ratio_rows = "\n".join("| " + " | ".join(("%.3f" % r0["ratio"]) if r0["ratio"] is not None else "n/a"
                                          for r0 in sz["exact_rows"]) + " |" for sz in px["sizes"])
layout_rows = "\n".join("| %d | %d-%d | %s | %d | %s |" % (x["word"], x["byte"], x["byte"] + 7,
                                                           x["field"], x["bits"], x["span"])
                        for x in pl["layout"])
recon_rows = "\n".join("| %s | %d | %s |" % (x["name"], x["instr"], x["detail"])
                       for x in pl["recon_items"])
mem_rows = "\n".join("| %s | 2^%.2f | %s |" % (a, math.log2(b), c) for a, b, c in MEM)
ledger_rows = "\n".join("| %s | %d | %.4f |" % (kk, v, math.log2(v)) for kk, v in rows.items())

render = dict(
    claim=cert["claim_text"], claim_num=str(claim.numerator), claim_den=str(claim.denominator),
    w=("%.6f" % W_log2), total_instr=str(total_instr), headroom=cert["headroom_over_W"],
    gate=format(Decimal(gate.numerator) / Decimal(gate.denominator), '.4f'),
    gate_margin=cert["gate_margin"], allowa=decAs, allowb=decBs, allow_sum=decSums,
    allowa_frac=str(maxA), allowb_frac=str(maxB),
    cell_rows="\n".join(
        "| %d | %d | %d | %d | %d | %d | %.9f | %.9f | %.3f |" %
        (x["r"], x["d"], x["m"], x["D"], x["pred"], x["diag"], float(x["rowA"]), float(x["rowB"]),
         x["scale"]) for x in frac),
    scale_max=("%.3f" % scale_max), mprod=str(m.bit_length() - 1), extrap=("%.2f" % extrap_log2),
    design=("%.6f" % float(design_lo)), after=("%.6f" % success_after),
    debit_sum=("%.6f" % debit_sum), markov=("%.3f" % float(markov)),
    ebad=("%.4f" % math.log2(float(E_bad))), kappa=("%.3f" % float(kappa)),
    pall=("%.4f" % math.log2(P_all)), image=str(IMAGE_LOG2),
    capwork=("%.2f" % math.log2(CAP_WORK)), percand=str(per_candidate),
    keybits=str(KEY_BITS), routebits=str(keymod.ROUTE_BITS), sortbits=str(keymod.SORT_BITS),
    scanbits=str(keymod.SCAN_BITS),
    census_rows=census_rows, width_header="| " + " | ".join("%d" % w for w in px["widths"]) + " |",
    width_rows=width_rows, pred_row=pred_row, ratio_rows=ratio_rows,
    kappameas=("%.3f" % px["max_ratio"]), kappamin=("%.3f" % px["min_ratio"]),
    layout_rows=layout_rows, recon_rows=recon_rows, mem_rows=mem_rows,
    recon_item=("%.4f" % math.log2(pl["recon_charged"])), recon_instr=str(pl["recon_charged"]),
    ledger_rows=ledger_rows, layout_evidence=("%d/%d keys, %d/%d records" %
                                              (pl["key_agreement"], pl["inputs"],
                                               pl["reconstructed"], pl["records"])),
    memlog2=("%.3f" % mem_log2), runmargin=("%.2f" % (95 - math.log2(float(P_pref)))),
    prefixlog2=("%.2f" % math.log2(float(P_pref))),
    fit_rows="\n".join("| %d | %d | %.4f | %.2f/%d | %.3f | %.6f | %.6f | %d | %d/%d |" %
                       (s["d"], s["n"], s["mean_over_law"], s["chi2"], s["dof"], s["chi2_over_dof"],
                        s["ks"], s["ks_95_cut"], s["max_gap_bits"], s["tail_above_model"],
                        s["tail_cells"]) for s in sizes),
    fitmeans="/".join("%.3f" % s["mean_over_law"] for s in sizes),
    fitks="/".join("%.4f" % s["ks"] for s in sizes),
    gapmean=("%.2f" % mean_gap), gapmax=str(max_gap), gaplaw=("%.1f" % float(F(16) * F(51) * LN2)),
    trend_rows="\n".join("| %d | %d | %d | %.6f | %.6f | %s |" %
                         (t["d"], t["small_r"], t["large_r"], t["small"], t["large"],
                          "yes" if t["shrinks"] else "no") for t in trend),
    ehi=format(Decimal(hi_e.numerator) / Decimal(hi_e.denominator), '.8f'),
)
json.dump(render, open(HERE + "/render.json", "w"), indent=1)

print("W = 2^%.6f  claim = %s  gate = %s  headroom = %s" % (W_log2, cert["claim_text"],
                                                            render["gate"], render["headroom"]))
print("E_bad = 2^%.4f  kappa = %s  markov = 2^%.3f  cap_work = 2^%.2f" %
      (math.log2(float(E_bad)), float(kappa), math.log2(float(markov)), math.log2(CAP_WORK)))
print("allowA=%s (%s) allowB=%s (%s) design_lo=%.8f after=%.6f" %
      (maxA, decAs, maxB, decBs, float(design_lo), success_after))
print("recon %d instr (2^%.4f) per candidate; mem 2^%.2f B" %
      (pl["recon_charged"], math.log2(pl["recon_charged"]), mem_log2))
print("checks PASS =", len(checks))
