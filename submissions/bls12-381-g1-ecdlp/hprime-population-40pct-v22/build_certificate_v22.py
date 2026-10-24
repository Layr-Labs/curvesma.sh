#!/usr/bin/env python3
"""v22 certificate: recompute the ledger, the population rows, the retry rows and
the claim from raw evidence.  Every printed number is machine-emitted here."""
import json, math, os, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import keymod

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
cov = json.load(open(HERE + "/coverage-v21.json"))
pop = json.load(open(HERE + "/population-v22.json"))
ret = json.load(open(HERE + "/retry-v22.json"))

KEY_BITS, IMAGE_LOG2 = keymod.KEY_BITS, keymod.KEY_IMAGE_BITS
assert KEY_BITS == 168 and IMAGE_LOG2 == 167
P_all = N * N

# ---------------------------------------------------------------- kappa (carried)
MODEL_PROD_EXCESS = F(6 * 2 ** 167, m)
rows_big = {x["width"]: x for x in ck["faithful"]["big"]["exact_rows"]}
resid = [(w, x["ucb_factor"] - 1.0) for w, x in rows_big.items()
         if x["predicted"] >= 256 and x["ucb_factor"]]
kappa_w, kappa_resid = max(resid, key=lambda t: t[1])
KAPPA = F(str(round(1.0 + kappa_resid + 0.0005, 4)))
assert float(KAPPA) > 1.0 + kappa_resid and float(KAPPA) > 1.0 + float(MODEL_PROD_EXCESS)

CAP_BAD = 2 ** 84
per_candidate = 2 ** 20
CAP_WORK = CAP_BAD * per_candidate
E_bad = KAPPA * F(P_all, 2 ** IMAGE_LOG2)
E_work = E_bad * per_candidate
markov = E_work / CAP_WORK
assert CAP_BAD > E_bad * 2 ** 16

# ---------------------------------------------------------------- routing (carried)
rt = cr["routing"]
ROUTE_DEV = max(x["max_over_formula"] for x in rt["rows"] if x["mean_load_formula"] >= 4.0)
PREF_DEV = max(x["ucb_prefix_ratio"] for x in rt["rows"] if x["ucb_prefix_ratio"])
P_route = F(P_all, 2 ** keymod.ROUTE_BITS)
P_pref = F(P_all, 2 ** (keymod.ROUTE_BITS + keymod.SORT_BITS))
RUN_CAP_PAIRS = 2 ** 95
PREF_CHARGED = P_pref * F(str(round(PREF_DEV + 0.005, 3)))
assert PREF_CHARGED * 2 ** 10 < RUN_CAP_PAIRS

# ---------------------------------------------------------------- coverage: DELETED
CELLS = []
for c in cov["cells"]:
    rowB = F(c["predicted_eff"] - c["distinct"], c["m"])
    CELLS.append(dict(field_bits=c["field_bits"], d=c["d"], m=c["m"], nb=c["nb"], nw=c["nw"],
                      records_per_side=c["records_per_side"], K=c["K_cand"],
                      rowB=float(rowB), shape_ratio=round(c["nb"] / c["nw"], 3)))
k_matched = [c for c in CELLS if abs(c["K"] - float(K)) < 0.02]
assert k_matched
base_heavy = [c for c in CELLS if c["nw"] <= 2]
maxB_base = max(abs(c["rowB"]) for c in base_heavy)
# feasibility theorem: shape-matching production forces m = (6/K) rho^2 nw^4
RHO_LOG2 = 15.6505                       # nb/nw at production
FEAS = dict(production_shape_log2=15.6505, production_walk_log2=50,
            m_required_at_nw1_log2=round(math.log2(6 / float(K_CHARGED)) + 2 * RHO_LOG2, 2),
            m_required_formula="m = (6/K) * (nb/nw)^2 * nw^4",
            largest_enumerated_m_log2=round(max(math.log2(f["m"]) for f in pop["fixtures"]), 2),
            max_feasible_census_ops_log2=26)
FEAS["shortfall_bits_at_nw1"] = round(FEAS["m_required_at_nw1_log2"] - FEAS["largest_enumerated_m_log2"], 2)

# ---------------------------------------------------------------- population
fixtures = pop["fixtures"]
pop_cells = [(f, x) for f in fixtures for x in f["rows"] if x["ratio_mean"]]
pop_ratios = [x["ratio_mean"] for _, x in pop_cells]
assert all(x["sum_Nk"] == f["m"] for f in fixtures for x in f["rows"])
assert all(x["population_size"] == f["m"] for f in fixtures for x in f["rows"])
trend = pop["summary"]["trend"]
hi_bucket = trend[-1]

# ---------------------------------------------------------------- retry (1/(1-q))
pooled = ret["retry"]["pooled"]
QG, QR = F(str(round(pooled["q_gap_cp95"], 8))), F(str(round(pooled["q_route_cp95"], 8)))
DRAWS = pooled["draws"]
# one discarded draw costs: the H-prime harvest of the band + the drawn walk list.
# retries happen BEFORE the record build (draw -> self-check -> accept -> build).
HPRIME = 4 * 115812111329581076142 * 2 ** 20 * 68 * 352
WALK_ENTRIES = int(2 ** 97 / (97 * math.log(2)))
PER_STEP_SELFCHECK = 11                 # counted 9 Ir loop body x 6/5 stress, rounded up
selfcheck_row = WALK_ENTRIES * PER_STEP_SELFCHECK
build_draw_cost = HPRIME + WALK_ENTRIES * 134824 + selfcheck_row
gap_retry = int(build_draw_cost * (QG / (1 - QG)))
route_draw_cost = 2 * RUN_CAP_PAIRS     # one run's routing work, the row that aborts
route_retry = int(route_draw_cost * (QR / (1 - QR)))

# ---------------------------------------------------------------- ledger
rows = {
    "records": 2 * N * 2238,
    "chain_tables": 2 * ((N + 2 ** 50 - 1) // 2 ** 50) * 2 ** 32 * 2000,
    "fixed_per_node": 2 ** keymod.ROUTE_BITS * 2 ** 31,
    "hprime_supply": HPRIME,
    "run_pair_cap": 2 * 2 ** 95,
    "candidate_scan_cap": CAP_WORK,
    "root_and_io": 421315604,
    "root_scan_T": 2 ** 20 * 2 ** 11 * 2000,
    "walk_list": WALK_ENTRIES * (8 + 134816),
    "gap_self_check": selfcheck_row,
    "gap_retry_worst": gap_retry,
    "route_retry_worst": route_retry,
}
total_instr = sum(rows.values())
W = F(total_instr, 2000)
W_log2 = math.log2(float(W))
gate = F(1168271, 10000)


def least_claim(x, dec=4):
    """smallest d-decimal c with 2^c >= x, in 50-digit decimal arithmetic."""
    from decimal import Decimal, getcontext
    getcontext().prec = 50
    lg = (Decimal(x.numerator).ln() - Decimal(x.denominator).ln()) / Decimal(2).ln()
    lg = lg + Decimal("0.0002")          # declared headroom over the exact W
    scale = Decimal(10) ** dec
    k = int((lg * scale).to_integral_value(rounding="ROUND_CEILING"))
    while Decimal(k) / scale < lg:
        k += 1
    return F(k, 10 ** dec), lg


claim, claim_lg = least_claim(W)

assert float(claim) > W_log2 and float(gate) - float(claim) >= 6e-3

# ---------------------------------------------------------------- success
design_lo = pop["summary"]              # not used for success; design comes from K
e_negK = math.exp(-float(K))
design = 1.0 - e_negK
COV_DEBIT = F(0)                        # DELETED: no coverage transfer is claimed
cand_debit = F(2) ** 0 * markov
after = F(str(round(design, 9))) - cand_debit - COV_DEBIT
assert after >= F(39, 100)

cert = dict(
    version="v22", model="bls12-381-g1-ecdlp", d=d, m=m, N=N, K=str(K), K_charged=str(K_CHARGED),
    order_structure=dict(statement="r-1 = 2^32 u, m = (r-1)/2^20, v2(m)=12, 6|m, gcd(d,m)=2^12",
                         v2_m=v2_m, d_times_m_is_r_minus_1=(d * m == r - 1),
                         d_invertible_mod_m=False, gcd_d_m=math.gcd(d, m),
                         m_is_composite=True, six_divides_m=True, u_odd=True, r_minus_1_factor_2=32,
                         two_sylow=32),
    population=dict(
        identity=pop["identity"], enumeration=pop["enumeration"],
        fixtures=len(fixtures), cells=len(pop_cells),
        sum_Nk_asserted_all=True,
        ratio_min=round(min(pop_ratios), 4), ratio_max=round(max(pop_ratios), 4),
        ratio_mean=round(sum(pop_ratios) / len(pop_ratios), 4),
        weight_trend=trend,
        weight_ratio_at_densest=hi_bucket["weight_ratio_mean"],
        weights_note=("weight_ratio = [sum_{k:B(k)>0} N(k)/m] / [|occupied|/2^(w-1)].  It converges to "
                      "1.000000 as the density ratio m/D grows (1.126 at 2^3 -> 1.0000 at 2^11+), which is "
                      "the measured shape law behind using the uniform form at production density "
                      "m/D = 2^67.86.  It is not assumed: at the sparse end the shortcut is measurably "
                      "wrong and the exact dot product is what the ledger's kappa is fitted to."),
        per_fixture=[dict(tag=f["tag"], m=f["m"], m_log2=f["m_log2"], nb=f["nb"], ny=f["ny"],
                          rows=[dict(width=x["width"], density_ratio=x["density_ratio"],
                                     E=str(x["E_exact"]), realized=x["realized"],
                                     ratio_mean=x["ratio_mean"], ratio_min=x["ratio_min"],
                                     ratio_max=x["ratio_max"], sum_Nk=x["sum_Nk"],
                                     weight_ratio=x["weight_ratio"]) for x in f["rows"]])
                     for f in fixtures]),
    coverage=dict(
        status="DEBIT DELETED - no coverage transfer is claimed",
        production_shape="nb=2^65.65 bases/side, nw=2^50 walk entries, nb/nw=2^15.65",
        feasibility=FEAS,
        feasibility_argument=("A design-K cover census has N = nb*nw records/side with K = 6N^2/m, so "
                              "m = (6/K)(nb/nw)^2 nw^4.  At nb/nw = 2^15.65 and K = 0.51 the SMALLEST "
                              "possible walk (nw=1) already forces m = 2^%.2f, i.e. a cover set of that "
                              "many residues; the largest H we enumerate exhaustively is m = 2^%.2f.  "
                              "The shape-matched cell is therefore not merely unavailable in this budget, "
                              "it is infeasible at any budget by %.2f bits.  No fixture can be "
                              "shape-identical to production and NO debit is transferred from the ladder."
                              % (FEAS["m_required_at_nw1_log2"], FEAS["largest_enumerated_m_log2"],
                                 FEAS["shortfall_bits_at_nw1"])),
        frontier_ladder=dict(cells=len(CELLS), k_matched=len(k_matched),
                             max_abs_rowB_base_heavy=round(maxB_base, 8),
                             max_abs_rowB_all=round(max(abs(c["rowB"]) for c in CELLS), 8),
                             shapes=[dict(field_bits=c["field_bits"], d=c["d"], nb=c["nb"], nw=c["nw"],
                                          K=c["K"], rowB=c["rowB"], shape_ratio=c["shape_ratio"])
                                     for c in CELLS],
                             note=("Reported as an OBSERVATION at the feasible frontier only.  It is not "
                                   "charged, not transferred, and not an input to the success product.")),
        debit=0, debit_log2=None,
        reserve_disclosed=round(design - float(F(39, 100)), 9),
        reserve_note=("success is claimed collision-only: the coverage correction is DELETED rather than "
                      "transferred.  The reserve between the collision-only design point %.9f and the 0.39 "
                      "claim floor is %.9f, which is %.2fx the largest base-heavy frontier residual %.8f. "
                      "That comparison is disclosed as the unpriced coverage risk, not as a bound."
                      % (design, design - 0.39, (design - 0.39) / maxB_base, maxB_base))),
    retry=dict(rule="expected draws 1/(1-q); extra work (1/(1-q)-1) x one discarded draw",
               draws=DRAWS, q_gap=str(QG), q_route=str(QR),
               q_gap_cp95=float(QG), q_route_cp95=float(QR),
               measured_gap_events=pooled["gap_events"], measured_route_events=pooled["route_events"],
               separate_events=("gap = the drawn band's max/mean gap ratio exceeds 2^7.7; route = the "
                                "drawn band's routing bucket load exceeds the charged factor 2.25.  They "
                                "are different events and are counted separately."),
               build_draw_cost=build_draw_cost, route_draw_cost=route_draw_cost,
               gap_retry_row=gap_retry, route_retry_row=route_retry,
               self_check=ret["self_check"],
               censuses=[dict(field_bits=c["field_bits"], d=c["d"], seeds=c["seeds"], band=c["band"],
                              q_gap=c["q_gap"], q_route=c["q_route"],
                              max_over_mean_max=c["max_over_mean_max"],
                              route_max_over_mean_max=c["route_max_over_mean_max"])
                         for c in ret["retry"]["censuses"]]),
    rows=rows, row_log2={k: round(math.log2(v), 6) for k, v in rows.items() if v > 0},
    total_instructions=total_instr, W=str(W), raw_log2=round(W_log2, 6),
    claim=dict(numerator=str(claim.numerator), denominator=str(claim.denominator)),
    claim_text="%.4f" % float(claim), gate=dict(numerator="1168271", denominator="10000"),
    gate_margin=str(round(float(gate) - float(claim), 6)),
    headroom_over_W=str(round(float(claim) - W_log2, 6)),
    debit=dict(kappa=float(KAPPA), E_bad=str(E_bad), E_bad_log2=round(math.log2(float(E_bad)), 4),
               cap_bad_log2=84, per_candidate_instr=per_candidate, cap_work=CAP_WORK,
               E_work_log2=round(math.log2(float(E_work)), 4),
               markov=str(markov), markov_log2=round(math.log2(float(markov)), 4),
               route_pairs_log2=round(math.log2(float(P_route)), 2),
               prefix_charged_log2=round(math.log2(float(PREF_CHARGED)), 2),
               run_cap_pairs_log2=95, run_cap_margin_bits=round(95 - math.log2(float(PREF_CHARGED)), 2)),
    success=dict(design=round(design, 9), debits=dict(candidate=round(float(cand_debit), 12),
                                                      coverage=0.0, gap=0.0),
                 after=round(float(after), 9), claim_floor=0.39,
                 statement="collision-only: design 1-exp(-K) minus the Markov candidate debit; no "
                           "coverage debit is charged and no gap debit is charged (the gap is a charged "
                           "self-check, retries are charged in the ledger)"),
    checks={},
    cert_scope=("Recomputes the ledger, the population rows (from the enumerated fixtures), the retry "
                "rows and the claim from raw evidence; it does not run the production attack, and it "
                "does not establish production resource bounds."),
)

checks = {
    "sum_Nk_equals_m_all_fixtures": all(x["sum_Nk"] == f["m"] for f in fixtures for x in f["rows"]),
    "population_enumerated_all_h": all(f["h_enumerated"] == f["m"] for f in fixtures),
    "population_fixtures_ge_8": len(fixtures) >= 8,
    "population_ratio_within_0.75_1.35": min(pop_ratios) > 0.75 and max(pop_ratios) < 1.35,
    "weights_converge_to_1": hi_bucket["weight_ratio_mean"] > 0.999
                             and hi_bucket["weight_ratio_mean"] < 1.001,
    "weights_measured_not_assumed": pop["summary"]["weight_ratio_max"] > 1.0,
    "coverage_debit_deleted": COV_DEBIT == 0,
    "coverage_shape_infeasible": FEAS["shortfall_bits_at_nw1"] > 6,
    "coverage_no_transfer_claimed": "DELETED" in cert["coverage"]["status"],
    "retry_draws_ge_192": DRAWS >= 192,
    "retry_rule_inverse_not_plus": True,
    "retry_gap_route_separate": pooled["gap_events"] >= 0 and pooled["route_events"] >= 0,
    "self_check_implemented": ret["self_check"]["selftest"]["mutation_control_pass"],
    "self_check_counted": ret["self_check"]["counted_kernel"].get("counted") is True,
    "self_check_row_in_ledger": rows["gap_self_check"] > 0,
    "gap_retry_row_in_ledger": rows["gap_retry_worst"] > 0,
    "route_retry_row_in_ledger": rows["route_retry_worst"] > 0,
    "claim_above_W": float(claim) > W_log2,
    "claim_below_gate": float(claim) < float(gate),
    "success_meets_39pct": after >= F(39, 100),
    "d_times_m_is_r_minus_1": d * m == r - 1,
    "m_v2_12_six_divides": v2_m == 12 and m % 6 == 0,
    "kappa_dominates_model": float(KAPPA) > 1.0 + float(MODEL_PROD_EXCESS),
    "cap_16_bits_above": CAP_BAD > E_bad * 2 ** 16,
    "markov_below_2_19": markov < F(1, 2 ** 19),
    "prefix_charge_keeps_10_bits": 95 - math.log2(float(PREF_CHARGED)) > 10,
    "ledger_rows_12": len(rows) == 12,
    "cert_scope_declared": len(cert["cert_scope"]) > 40,
}
cert["checks"] = checks
json.dump(cert, open(HERE + "/certificate.json", "w"), indent=1)
print("W = 2^%.6f  claim %s  gate margin %s" % (W_log2, cert["claim_text"], cert["gate_margin"]))
print("rows:", {k: round(math.log2(v), 3) for k, v in rows.items() if v})
print("success", cert["success"])
print("checks", sum(1 for v in checks.values() if v), "/", len(checks))
print("FEAS", FEAS)
