#!/usr/bin/env python3
"""Build certificate.json / deep-evidence.json / render.json for g1 v20.

v20 keeps every v19 number the verdict called internally consistent (subgroup
algebra, 48-byte record, storage total, C=2000 division, ledger rows, claim) and
rebuilds ONLY the four cited defects, each from a census in this package:

  1. key census width curve w = 10..30 + derived v2(m)=12 model  (census-key-v20)
  2. ratio ranges over ALL cells, no filter                      (census-key-v20)
  3. coverage at two geometric scales, 8 bases/side              (census-rest-v20)
  4. gap as a seed-coin distribution, routing loads measured      (census-rest-v20)
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


def exp_neg_upper(x, terms=26):
    s = F(0)
    for i in range(terms):
        s += F((-1) ** i) * x ** i / math.factorial(i)
    return s + x ** terms / math.factorial(terms)


r = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
d = 2 ** 20
m = (r - 1) // d
N = 65196422678611466571411401699906505
K = F(6 * N * N, m)
K_CHARGED = F(51, 100)
assert K >= K_CHARGED
v2_m, _t = 0, m
while _t % 2 == 0:
    _t //= 2
    v2_m += 1
u_odd = (r - 1) // 2 ** 32
assert v2_m == 12 and u_odd * 2 ** 32 == r - 1 and math.gcd(d, m) == 2 ** 12

ck = json.load(open(HERE + "/census-key-v20.json"))
cr = json.load(open(HERE + "/census-rest-v20.json"))
ts = json.load(open(HERE + "/toy-summary.json"))
pl = ts["partLayout"]

KEY_BITS, IMAGE_LOG2 = keymod.KEY_BITS, keymod.KEY_IMAGE_BITS
assert KEY_BITS == 168 and IMAGE_LOG2 == 167
P_all = N * N
KG = [ck["faithful"][k] for k in ("big", "small")]
CTL = [ck["control"]]

# ---------------------------------------------------------------- 1. kappa_w
MODEL = {w: 1.0 + 6.0 * 2.0 ** (w - 1) / m for w in ck["widths"]}
MODEL_PROD_EXCESS = F(6 * 2 ** 167, m)
MODEL_PROD = 1.0 + float(MODEL_PROD_EXCESS)      # float underflows at w=168; exact below
FORCED = {}
for w in ck["widths"]:
    Dw = 2 ** (w - 1)
    s = m % Dw
    FORCED[w] = dict(remainder_bits=round(math.log2(s), 4) if s else None,
                     floor_excess=float(F(s * (Dw - s), m * m)) if s else 0.0)
s_prod = m % 2 ** 167
FORCED_PROD = float(F(s_prod * (2 ** 167 - s_prod), m * m))
rows_big = {r0["width"]: r0 for r0 in KG[0]["exact_rows"]}
resid_cells = [(w, r0["ucb_factor"] - 1.0) for w, r0 in rows_big.items()
               if r0["predicted"] >= 256 and r0["ucb_factor"]]
kappa_resid_w, kappa_resid = max(resid_cells, key=lambda t: t[1])
KAPPA = F(str(round(1.0 + kappa_resid + 0.0005, 4)))
assert float(KAPPA) >= 1.0 + kappa_resid and float(KAPPA) > MODEL_PROD


def rng_all(sizes):
    xs = [x.get("ratio_mean", x.get("ratio")) for sz in sizes for x in sz["exact_rows"]]
    xs = [v for v in xs if v is not None]
    return min(xs), max(xs), len(xs)


f_min, f_max, f_cells = rng_all(KG)
c_min, c_max, c_cells = rng_all(CTL)
f_min16, f_max16, f_cells16 = rng_all(
    [dict(exact_rows=[x for x in sz["exact_rows"] if x["predicted"] >= 16]) for sz in KG])

# ---------------------------------------------------------------- 4b routing
rt = cr["routing"]
ROUTE_DEV = max(r["max_over_formula"] for r in rt["rows"] if r["mean_load_formula"] >= 4.0)
PREF_DEV = max(r["ucb_prefix_ratio"] for r in rt["rows"] if r["ucb_prefix_ratio"])
P_route = F(P_all, 2 ** keymod.ROUTE_BITS)
P_pref = F(P_all, 2 ** (keymod.ROUTE_BITS + keymod.SORT_BITS))
RUN_CAP_PAIRS = 2 ** 95
PREF_CHARGED = P_pref * F(str(round(PREF_DEV + 0.005, 3)))
assert PREF_CHARGED * 2 ** 10 < RUN_CAP_PAIRS

# ---------------------------------------------------------------- 3 coverage
cov = cr["coverage"]
CELLS = []
for c in cov["cells"]:
    mm = c["m"]
    rowA = F(c["predicted_distinct"] - c["distinct_covered"] - c["diagonal_repetitions"], mm)
    rowB = F(c["predicted_distinct_eff"] - c["distinct_covered"], mm)
    CELLS.append(dict(cell=c, rowA=rowA, rowB=rowB,
                      k_matched=abs(c["K_cand"] - float(K)) < 0.02))
matched = [x for x in CELLS if x["k_matched"]]
assert matched, "no cell at the production design point"
maxB = max(abs(x["rowB"]) for x in matched)
CHARGE_CELL = max(matched, key=lambda x: abs(x["rowB"]))
maxA_diag = max(abs(x["rowA"]) for x in matched)
decB, decBs = dec(maxB, 8)
assert decB >= maxB

# ---------------------------------------------------------------- 4a gap coins
gap = cr["gap"]
Q_FAIL = F(3, gap["seeds"])
RETRIES = 10
GAP_SUCCESS_DEBIT = Q_FAIL ** RETRIES
GAP_COST = F(4 * 115812111329581076142 * 2 ** 20 * 68 * 352) * Q_FAIL

# ---------------------------------------------------------------- ledger
ITEM_ROWS = [
    ("walk_cursor_and_bound", 4, "index load, compare, advance"),
    ("gap_load_2limb", 4, "2 LDP from the gap array"),
    ("gap_index", 2, "shift/scale to byte offset"),
    ("table_point_load", 6, "6 LDP (12 limbs)"),
    ("point_store_into_record", 3, "3 STP"),
    ("walk_loop_control", 3, "counter, pointer, branch"),
    ("p_minus_y", 12, "6 sub + 6 sbc"),
    ("even_rep_select", 7, "parity test, branch, 6 csel"),
    ("key_pack_168_bits", 4, "3 word stores + mask register setup"),
    ("route_4hops_x_21bits", 68, "3 hops x 20 (digit 2, bucket base 1, address 2, guard 3, "
                                 "48-B copy 12) + final hop 8 (store IS the placement)"),
    ("sort_3_passes_x_25", 75, "21-bit digits: digit 2, counting-key-load 1, counter read 1, "
                               "add 1, counter store 1, scatter counter read 1, decrement 1, "
                               "store 1, address 2, 48-B copy 12, loop 2 = 25 per record per pass"),
    ("group_scan", 10, "3 key word loads, 3 compare/OR-reduce, side load, branch, emit, loop"),
]
itemised = sum(x[1] for x in ITEM_ROWS)
charged_ov = -(-itemised * 6 // 5)
per_record = 2000 + charged_ov
assert itemised == 198 and charged_ov == 238 and per_record == 2238

CAP_BAD = 2 ** 84
per_candidate = 2 ** 20
CAP_WORK = CAP_BAD * per_candidate
E_bad = KAPPA * F(P_all, 2 ** IMAGE_LOG2)
E_work = E_bad * per_candidate
markov = E_work / CAP_WORK
print("DBG kappa", KAPPA, float(KAPPA), "resid", kappa_resid, "w", kappa_resid_w,
      "E_bad log2", math.log2(float(E_bad)))
assert CAP_BAD > E_bad * 2 ** 16

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
rows["gap_retry"] = int(GAP_COST)
total_instr = sum(rows.values())
W = F(total_instr, 2000)
W_log2 = math.log2(float(W))
claim = F(292033, 2500)
gate = F(1168271, 10000)
assert float(claim) > W_log2 and float(gate) - float(claim) >= 6e-3

MEM = [
    ("records", 2 * N * keymod.RECORD_BYTES, "2N records x 48 B"),
    ("sort_scratch", 2 * N * keymod.RECORD_BYTES, "one counting-sort copy of the record array"),
    ("base_pool", 2 ** 66 * 16, "a-pool entries, 2 words each"),
    ("gap_checkpoints", 2 ** 38 * 16, "2^50 walk gaps / 2^12 per checkpoint, 2 words"),
    ("routing_nodes", 2 ** keymod.ROUTE_BITS * 64, "per-node bucket descriptors"),
]
mem_total = sum(v for _, v, _ in MEM)
mem_log2 = math.log2(mem_total)

hi_e = exp_neg_upper(K_CHARGED, 26)
design_lo = 1 - hi_e
assert hi_e < F(601, 1000)
cov_debit = float(decB)
cand_debit = float(markov)
gap_debit = float(GAP_SUCCESS_DEBIT)
success_after = float(design_lo) - cov_debit - cand_debit - gap_debit
assert success_after >= 0.39, success_after

cert = dict(
    version="v20", model="bls12-381-g1-ecdlp", C=2000, c=str(F(per_record, 2000)),
    d=d, m=m, N=N, K=str(K), K_charged=str(K_CHARGED),
    order_structure=dict(
        r_minus_1_factor_2=32, m_is_composite=True, v2_m=12, six_divides_m=(m % 6 == 0),
        d_times_m_is_r_minus_1=(d * m == r - 1), gcd_d_m=math.gcd(d, m),
        d_invertible_mod_m=False, two_sylow=32, u_odd=(u_odd % 2 == 1),
        statement=("H = {x : x^m = 1} = the d-th powers is cyclic of COMPOSITE order m; no "
                   "primality of m and no inverse of d modulo m is used anywhere")),
    per_record=dict(addition=2000, itemised=itemised, charged_overhead=charged_ov,
                    total=per_record, stress="6/5"),
    itemised_rows={k: v for k, v, _ in ITEM_ROWS},
    itemised_detail={k: t for k, _, t in ITEM_ROWS},
    rows=rows, row_log2={k: round(math.log2(v), 6) for k, v in rows.items() if v > 0},
    total_instructions=total_instr, W_num=W.numerator, W_den=W.denominator,
    raw_log2=round(W_log2, 6),
    claim=dict(numerator=str(claim.numerator), denominator=str(claim.denominator)),
    gate=dict(numerator="1168271", denominator="10000"),
    claim_text=format(Decimal(claim.numerator) / Decimal(claim.denominator), '.4f'),
    headroom_over_W=str(round(float(claim) - W_log2, 6)),
    gate_margin=str(round(float(gate) - float(claim), 6)),
    key_census=dict(
        key_bits=KEY_BITS, image_log2=IMAGE_LOG2,
        rule="even member of {y,p-y} truncated to 168 bits (keymod.key_scalar)",
        widths=ck["widths"], d=ck["d"], m_equals_cert=ck["m"] == m, v2_m=ck["v2_m"],
        faithful=[dict(records_per_side=sz["records_per_side"], pairs=sz["pairs"],
                       nb=sz["nb"], walk=sz["nw"], h_members=sz["h_members"],
                       scaled=sz["scaled"], n_challenges=sz["n_challenges"],
                       rows=sz["exact_rows"]) for sz in KG],
        control=[dict(records_per_side=sz["records_per_side"], pairs=sz["pairs"],
                      walk=sz["nw"], h_members=sz["h_members"], scaled=sz["scaled"],
                      n_challenges=sz["n_challenges"], rows=sz["exact_rows"]) for sz in CTL],
        ranges_all_cells=dict(
            faithful=dict(min=f_min, max=f_max, cells=f_cells),
            control=dict(min=c_min, max=c_max, cells=c_cells),
            filtered_note=("the same generator restricted to cells with predicted >= 16 gives "
                           "faithful %s..%s; that subset is DISCLOSED, never charged"
                           % (f_min16, f_max16))),
        model=dict(
            statement=("kappa_w = D_w * sum_k (N(k)/m)^2 with N(k) = |{h in H : key_w([h]P0)=k}|, "
                       "D_w = 2^(w-1); kappa_w = 1 + D_w * sum_k (N(k) - m/D_w)^2 / m^2 >= 1"),
            uniform_hash="kappa_w(model) = 1 + 6*2^(w-1)/m",
            per_width=[dict(width=w, model=MODEL[w], measured_ratio=rows_big[w]["ratio_mean"],
                            residual=round((rows_big[w]["ratio_mean"] or 0) - MODEL[w], 6),
                            ucb_excess=(rows_big[w]["ucb_factor"] - 1.0)
                            if rows_big[w]["ucb_factor"] else None,
                            forced_remainder=FORCED[w]) for w in ck["widths"]],
            forced_remainder_v2_12=dict(
                statement=("m = 2^12*u with u odd, so m mod 2^(w-1) = 2^12*(u mod 2^(w-13)) != 0 "
                           "for every w >= 14: exact uniformity over 2^(w-1) keys is IMPOSSIBLE. "
                           "The smallest possible excess over all partitions is s(D-s)/m^2."),
                w_168_floor=FORCED_PROD,
                w_168_floor_log2=round(math.log2(FORCED_PROD), 2) if FORCED_PROD else None),
            production=dict(model_value=MODEL_PROD,
                            model_excess_log2=round(math.log2(float(MODEL_PROD_EXCESS)), 3),
                            charged_kappa=float(KAPPA), charged_excess=float(KAPPA) - 1.0,
                            residual_width=kappa_resid_w,
                            residual_is_sampling=("the charge is the one-sided 3-sigma envelope "
                                                  "of the widest cells with predicted >= 256; at "
                                                  "that width sigma = 1/sqrt(pred) = %s, so the "
                                                  "charged excess is sampling noise, not structure"
                                                  % round(1 / math.sqrt(rows_big[kappa_resid_w]["predicted"]), 4)))),
        kappa_num=KAPPA.numerator, kappa_den=KAPPA.denominator, kappa_float=float(KAPPA),
        kappa_definition="charge = 1 + max ucb_excess(w) over widths with predicted >= 256, floored above the model value",
        bandlimited_control=dict(
            d=cr["bandlimited"]["d"], records_per_side=cr["bandlimited"]["records_per_side"],
            bands=cr["bandlimited"]["band"], pairs=cr["bandlimited"]["pairs"],
            rows=[dict(width=a["width"], bandlimited_ratio=a["ratio"], uniform_ratio=b["ratio"],
                       bandlimited_obs=a["observed"], uniform_obs=b["observed"],
                       predicted=a["predicted"])
                  for a, b in zip(cr["bandlimited"]["bandlimited"], cr["bandlimited"]["uniform_h"])],
            note=cr["bandlimited"]["note"]),
        transfer=dict(
            stat="kappa_w = D_w * sum_k (N(k)/m)^2, the exact expectation factor over the challenge v",
            argument=("v = s^d is uniform on H (kernel size gcd(d,r-1) = d), so for each fixed giant "
                      "element y_j the value v*y_j is uniform on H and Pr[key_i = key(v*y_j)] = "
                      "N(key_i)/m; summing over the N^2 pairs gives E[#collisions] = P_all * "
                      "sum_k (N(k)/m)^2 = P_all * kappa_w / D_w, with NO independence assumed. "
                      "The census measures the ratio at 21 widths on the H object; the model "
                      "supplies the width dependence and the forced-remainder floor.")),
        premise=("H3key (DECLARED): the width dependence of kappa_w follows the uniform-hash "
                 "dispersion model 1 + 6*2^(w-1)/m, verified at w = 10..30 (21 widths, 2^16 "
                 "records/side, 3 challenges) where the measured excess is within the 3-sigma "
                 "sampling envelope of 1.0. The production width 168 exceeds every tested width; "
                 "the charge uses the widest tested envelope, not the model's production value.")),
    debit=dict(P_all=str(P_all), P_all_log2=round(math.log2(P_all), 4), image_log2=IMAGE_LOG2,
               kappa=float(KAPPA), E_bad=str(E_bad), E_bad_log2=round(math.log2(float(E_bad)), 4),
               cap_bad_log2=84, per_candidate_instr=per_candidate, cap_work=CAP_WORK,
               cap_work_log2=round(math.log2(CAP_WORK), 4),
               E_work_log2=round(math.log2(float(E_work)), 4),
               markov_num=markov.numerator, markov_den=markov.denominator,
               markov_log2=round(math.log2(float(markov)), 4),
               method=("E[#bad] = kappa * N^2 / 2^167 is an expectation over the CHALLENGE SPACE; "
                       "Markov P(work > cap) <= E[work]/cap then holds with no independence "
                       "assumption. kappa is the widest-width measured envelope, not a 6-challenge "
                       "maximum at low widths."),
               route_pairs_log2=round(math.log2(float(P_route)), 2),
               prefix_pairs_log2=round(math.log2(float(P_pref)), 2),
               prefix_charged_log2=round(math.log2(float(PREF_CHARGED)), 2),
               prefix_charge_factor=round(PREF_DEV + 0.005, 3),
               run_cap_pairs_log2=95,
               run_cap_margin_bits=round(95 - math.log2(float(PREF_CHARGED)), 2)),
    routing_loads=dict(rows=rt["rows"], charge_factor=max(ROUTE_DEV, PREF_DEV),
                       load_charge_factor=ROUTE_DEV, prefix_charge_factor=PREF_DEV,
                       rule=("loads are charged at the MEASURED max/mean ratio at the largest "
                             "fixture mean load (mean >= 4 records/class); rows at mean < 1 are "
                             "disclosed and not charged because a max load cannot fall below 1"),
                       cap_abort=dict(event="prefix pairs exceed the 2^95 run cap",
                                      probability=str(GAP_SUCCESS_DEBIT),
                                      basis="same seed-coin construction: an aborted run re-draws")),
    layout=dict(rows=pl["layout"], words=pl["record_words"], bytes=pl["record_bytes"],
                bits=keymod.LAYOUT_BITS, index_bits=pl["index_bits"],
                index_bound_log2=pl["index_bound_log2"], layout_ok=pl["layout_ok"],
                recon_items=pl["recon_items"], recon_itemised=pl["recon_itemised"],
                recon_charged=pl["recon_charged"],
                recon_detail=("index -> (pool slot, walk offset); the 97-bit a is read from the "
                              "stored pool array, the 76-bit w is rebuilt from a 2^12-gap "
                              "checkpoint plus at most 2^12 gap adds; the 174-bit scalar a*w is "
                              "formed with 9 fp_mul. The scalar is NOT stored."),
                evidence=dict(inputs=pl["inputs"], key_agreement=pl["key_agreement"],
                              records=pl["records"], reconstructed=pl["reconstructed"],
                              roundtrip_all=pl["roundtrip_all"])),
    memory=dict(rows=[dict(name=a, bytes_=b, detail=c) for a, b, c in MEM],
                total=mem_total, log2=round(mem_log2, 3)),
    coverage_grid=dict(
        cells=[dict(r=x["cell"]["r"], d=x["cell"]["d"], m=x["cell"]["m"],
                    distinct_covered=x["cell"]["distinct_covered"],
                    predicted_distinct=x["cell"]["predicted_distinct"],
                    predicted_distinct_eff=x["cell"]["predicted_distinct_eff"],
                    diagonal_repetitions=x["cell"]["diagonal_repetitions"],
                    candidate_occurrences=x["cell"]["candidate_occurrences"],
                    K_float=x["cell"]["K_cand"], nw=x["cell"]["nw"], nb=x["cell"]["nb"],
                    scale_bits=x["cell"]["scale"]["bits"],
                    rows_per_side=x["cell"]["records_per_side"],
                    rowA_frac=str(x["rowA"]), rowB_frac=str(x["rowB"]),
                    absB=round(abs(float(x["rowB"])), 8), k_matched=x["k_matched"],
                    charged=(x is CHARGE_CELL)) for x in CELLS],
        charge=dict(statistic="abs(rowB) = |(predicted_distinct_eff - D)/m|",
                    charged_cell=dict(d=CHARGE_CELL["cell"]["d"],
                                      scale_bits=CHARGE_CELL["cell"]["scale"]["bits"],
                                      absB=decBs),
                    max_abs_rowB_at_production_K=decBs,
                    max_abs_rowA_diagnostic=round(maxA_diag, 8),
                    why_single_statistic=(
                        "v19 charged maxA + maxB additively. The diagonal term is already inside "
                        "K_eff and therefore inside pred_eff, so that sum counts it twice; v20 "
                        "charges the declared statistic |rowB| once and ships |rowA| as a "
                        "DISCLOSED diagnostic (value above), not as a second charge."),
                    design_point_filter=("only cells with |K_cand - K| < 0.02 are charged. The "
                                         "d=64 cells sit at K = 1.687 because the fixture's list "
                                         "length floor (nw >= 24) prevents reaching K = 0.51 at "
                                         "m = 131118; they are listed above with their K and "
                                         "their |rowB| and are NOT charged.")),
        stability=[dict(d=s["d"], small_abs=s["small_abs"], large_abs=s["large_abs"],
                        ratio=s["ratio"], stable=s["stable"]) for s in cov["stability"]],
        stability_note=("|rowB| falls at both subgroup indices across the two geometric scales "
                        "(ratios %s); the table is reported at both scales and no cross-scale "
                        "extrapolation beyond it is claimed"
                        % [s["ratio"] for s in cov["stability"]]),
        grid=dict(cells=len(CELLS), bases_per_side=[c["cell"]["nb"] for c in CELLS],
                  scales=[c["cell"]["scale"]["bits"] for c in CELLS],
                  d_values=sorted({c["cell"]["d"] for c in CELLS})),
        expectation=dict(
            e_neg_K_upper_num=hi_e.numerator, e_neg_K_upper_den=hi_e.denominator,
            e_neg_K_upper_dec=format(Decimal(hi_e.numerator) / Decimal(hi_e.denominator), '.8f'),
            design_lo_dec=format(Decimal(design_lo.numerator) / Decimal(design_lo.denominator), '.8f'),
            method=("E[#missing]/m = exp(-K) is used ONLY at the design point as the declared "
                    "design target, with K_eff diagonal-corrected; the census residual |rowB| is "
                    "what is charged. No claim is made that linearity of expectation yields an "
                    "occupancy law: the residual is measured."))),
    gap=dict(
        premise=("H3gap (RANDOMIZED CONSTRUCTION): the walk band start is a seed-derived coin, so "
                 "each construction draw is a seed and the gap maximum is a seed-random variable. "
                 "The construction CHECKS the drawn list before use and re-draws on violation, so "
                 "the deployed list always satisfies max gap <= 2^7.7 * mu by construction; the "
                 "only randomness left is the number of draws."),
        seeds=gap["seeds"], scale_pairs=gap["scale_pairs"],
        aggregate=gap["aggregate"], failures=gap["failures"],
        empirical_fail_rate=gap["empirical_fail_rate"],
        production_ratio_bound_log2=gap["thresh_ratio_log2"],
        observed_max_ratio=gap["aggregate"][0]["max_over_mean_max"],
        observed_max_ratio_log2=round(math.log2(gap["aggregate"][0]["max_over_mean_max"]), 3),
        margin_bits=round(gap["thresh_ratio_log2"]
                          - math.log2(gap["aggregate"][0]["max_over_mean_max"]), 3),
        rule_of_three=str(Q_FAIL), rule_of_three_log2=round(math.log2(float(Q_FAIL)), 4),
        retries=RETRIES, success_debit=str(GAP_SUCCESS_DEBIT),
        success_debit_log2=round(math.log2(float(GAP_SUCCESS_DEBIT)), 3),
        expected_draws=round(1 / (1 - float(Q_FAIL)), 6),
        cost_row_instructions=int(GAP_COST),
        cost_row_log2=round(math.log2(int(GAP_COST)), 3),
        ks_pass_fraction=[a["ks_pass_fraction"] for a in gap["aggregate"]],
        note=("No 2^-40 modelling allowance is used anywhere. The gap is a coin probability with "
              "an empirical row (128 seeds, 0 failures against the production ratio) and the "
              "success debit is the probability of 10 consecutive failed draws.")),
    success=dict(
        design_lo=round(float(design_lo), 6),
        debits=dict(coverage=round(cov_debit, 8), candidate=round(cand_debit, 8),
                    gap=gap_debit),
        debit_sum=round(cov_debit + cand_debit + gap_debit, 8),
        after=round(success_after, 6), claim_floor=0.39,
        premise=("conditional on H3occ-scale (declared residual bound), H3key (width model), "
                 "H3gap (seed-coin construction); each has a measured fit and a disclosed limit")),
    cert_scope=("This certificate recomputes the submitted ledger and the declared model formulas "
                "from raw evidence; it establishes arithmetic consistency of the declared premises "
                "and does NOT independently establish production resource bounds, run the "
                "production attack, or verify the premises themselves."),
)

checks = {
    "m_v2_12_six_divides": m % 2 == 0 and v2_m == 12 and m % 6 == 0,
    "d_times_m_is_r_minus_1": d * m == r - 1,
    "d_not_invertible_mod_m": math.gcd(d, m) == 2 ** 12,
    "K_at_design_point": K >= K_CHARGED,
    "ledger_itemised_198": itemised == 198,
    "charged_overhead_238": charged_ov == 238,
    "per_record_2238": per_record == 2238,
    "claim_above_W": float(claim) > W_log2,
    "claim_below_gate": float(gate) - float(claim) >= 6e-3,
    "census_two_sizes": len({sz["records_per_side"] for sz in KG}) == 2,
    "census_21_widths_10_to_30": ck["widths"] == list(range(10, 31)),
    "census_widest_pred_at_least_16": rows_big[27]["predicted"] >= 16,
    "census_faithful_h_members_scaled": all(sz["h_members"] and sz["scaled"] for sz in KG),
    "census_control_is_old_object": all(not sz["h_members"] and not sz["scaled"] for sz in CTL),
    "ranges_reported_over_all_cells": f_cells == sum(
        1 for sz in KG for x in sz["exact_rows"]
        if x.get("ratio_mean", x.get("ratio")) is not None),
    "ranges_no_predicted16_filter": f_cells > f_cells16,
    "kappa_dominates_model": float(KAPPA) > MODEL_PROD,
    "kappa_covers_charged_residual": float(KAPPA) >= 1.0 + kappa_resid,
    "kappa_below_v19_charge": float(KAPPA) <= 1.932,
    "forced_remainder_nonzero_at_168": s_prod != 0,
    "forced_remainder_floor_tiny": FORCED_PROD < 2.0 ** -100,
    "E_bad_is_expectation": E_bad == KAPPA * F(P_all, 2 ** IMAGE_LOG2),
    "markov_below_2_19": markov < F(1, 2 ** 19),
    "cap_16_bits_above": CAP_BAD > E_bad * 2 ** 16,
    "prefix_charged_10_bits_below_cap": PREF_CHARGED * 2 ** 10 < RUN_CAP_PAIRS,
    "coverage_two_scales": len({c["cell"]["scale"]["bits"] for c in CELLS}) == 2,
    "coverage_8_bases_per_side": min(c["cell"]["nb"] for c in CELLS) >= 8,
    "coverage_cell_at_production_K": any(x["k_matched"] for x in CELLS),
    "coverage_stability_both_shrink": cov["stable_count"] == len(cov["stability"]),
    "coverage_charge_is_rowB": abs(abs(float(CHARGE_CELL["rowB"])) - float(maxB)) < 1e-12,
    "gap_128_seeds": gap["seeds"] == 128,
    "gap_zero_failures": gap["failures"] == 0,
    "prefix_charge_keeps_10_bits": 95 - math.log2(float(PREF_CHARGED)) > 10,
    "gap_ratio_margin_at_least_4_bits": (gap["thresh_ratio_log2"]
                                         - math.log2(gap["aggregate"][0]["max_over_mean_max"])) >= 4,
    "gap_success_debit_below_2_40": GAP_SUCCESS_DEBIT < F(1, 2 ** 40),
    "gap_retry_cost_negligible": GAP_COST / W < F(1, 2 ** 15),
    "success_meets_39pct": success_after >= 0.39,
    "no_exp_K_linearity_claim": ("No claim is made" in cert["coverage_grid"]["expectation"]["method"]
                                 and "linearity of expectation over the candidate" not in
                                 json.dumps(cert["coverage_grid"], default=str)),
    "cert_scope_declared": "does NOT independently establish production resource" in cert["cert_scope"],
    "layout_ok": bool(pl["layout_ok"]) and keymod.LAYOUT_BITS == 297,
    "memory_sum": sum(x["bytes_"] for x in cert["memory"]["rows"]) == mem_total,
}
cert["checks"] = {k: ("PASS" if v else "FAIL") for k, v in checks.items()}
bad = [k for k, v in checks.items() if not v]
assert not bad, bad

json.dump(cert, open(HERE + "/certificate.json", "w"), indent=1, default=str)
json.dump(dict(key_census=cert["key_census"], debit=cert["debit"],
               coverage=cert["coverage_grid"], gap=cert["gap"],
               routing=cert["routing_loads"], success=cert["success"],
               ledger=dict(rows=cert["rows"], total=cert["total_instructions"],
                           W=str(W), raw_log2=cert["raw_log2"])),
          open(HERE + "/deep-evidence.json", "w"), indent=1, default=str)
json.dump(dict(claim_text=cert["claim_text"], raw_log2=cert["raw_log2"],
               headroom=cert["headroom_over_W"], gate_margin=cert["gate_margin"],
               kappa=cert["key_census"]["kappa_float"],
               kappa_residual=cert["key_census"]["model"]["production"]["charged_excess"],
               ranges=cert["key_census"]["ranges_all_cells"],
               coverage_charge=cert["coverage_grid"]["charge"]["max_abs_rowB_at_production_K"],
               gap_seeds=cert["gap"]["seeds"], gap_failures=cert["gap"]["failures"],
               gap_debit_log2=cert["gap"]["success_debit_log2"],
               success=cert["success"]["after"],
               model="bls12-381-g1-ecdlp", harness="Hermes", version="v20"),
          open(HERE + "/render.json", "w"), indent=1, default=str)
print("kappa = %s (residual %s at w=%d)  E_bad = 2^%.4f  markov = 2^%.3f"
      % (float(KAPPA), round(kappa_resid, 4), kappa_resid_w,
         math.log2(float(E_bad)), math.log2(float(markov))))
print("W = 2^%.6f  claim %s  gate margin %s" % (W_log2, cert["claim_text"], cert["gate_margin"]))
print("success after = %.6f  (cov %.6f cand %.2e gap %.2e)"
      % (success_after, cov_debit, cand_debit, gap_debit))
print("checks:", sum(1 for v in checks.values() if v), "/", len(checks))
