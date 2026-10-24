import json, math, os, random
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
PRIOR = os.path.join(os.path.dirname(HERE), "hprime-keymix-40pct-v33")
def P(f): return os.path.join(HERE, f) if os.path.exists(os.path.join(HERE, f)) else os.path.join(PRIOR, f)
D = json.load(open(P("design-v34.json")))
LAD = json.load(open(P("structured-ladder.json")))
CEN = json.load(open(P("census-prod-v34.json")))
ANC = json.load(open(P("price-anchors-v34.json")))
N, m, NB, ELL = D["N"], D["m"], D["nb"], D["ell"]
SIDES, C = D["sides"], D["C_instructions_per_unit"]
PRICE = dict(D["price_lattice"])
PER = D["per_record_instructions"]
HP = ANC["constructions"]["hprime_element"]["instructions"]
CT = ANC["constructions"]["chain_table_entry"]["instructions"]
MU = 2.0 ** 26.067
ROWS = {
 "walk_list": ("ell", "P_walk_entry + P_even_form"),
 "hprime_supply_construction": ("sides * ell", "HPRIME_ELEMENT"),
 "chain_table_construction": ("sides * chain_windows * (2 ** chain_entries_log2)", "CHAIN_ENTRY"),
 "fixed_per_node": ("nb", "P_route_sort"),
 "base_entry": ("sides * nb", "P_even_form"),
 "walker_state": ("sides * nb", "4"),
 "gap_point_tables": ("sides * nb * (2 ** gap_support_log2)",
                      "gap_entry_point_ops * C_instructions_per_unit"),
 "gap_self_check": ("ell - 1", "P_scan"),
 "gap_retry_worst": ("retry_draws * q_gap_cp95_num * ell * per_record_instructions // "
                     "(q_gap_cp95_den - q_gap_cp95_num) + 1", "1"),
 "candidate_scan_cap": ("2 ** candidates_log2", "per_candidate_instructions"),
 "run_pair_cap": ("2 ** run_pairs_log2", "run_pair_instructions"),
 "root_and_io": ("sides * (2 ** root_tests_log2)", "P_scan"),
 "root_scan_T": ("(2 ** root_tests_log2)", "P_route_sort"),
}
ns = {k: v for k, v in D.items() if k != "price_lattice"}
ns.update(PRICE); ns.update(HPRIME_ELEMENT=HP, CHAIN_ENTRY=CT)
EV = lambda e: eval(e, {"__builtins__": {}}, ns)
rows, deriv = {}, {}
for name, (ce, pe) in ROWS.items():
    c, p = EV(ce), EV(pe)
    rows[name] = -(-c * p // 1)
    deriv[name] = dict(count_expr=ce, price_expr=pe, count=str(c), price=str(p),
                       value=str(rows[name]),
                       anchor=("constructions/hprime_element" if pe == "HPRIME_ELEMENT" else
                               "constructions/chain_table_entry" if pe == "CHAIN_ENTRY" else
                               "composites/" + pe if pe in ANC["composites"] else
                               "items/" + pe if pe in ANC["items"] else "P_record"))
rows["records"] = SIDES * (NB * ELL) * PER
deriv["records"] = dict(count_expr="sides * nb * ell", price_expr="per_record_instructions",
                        count=str(SIDES * NB * ELL), price=str(PER), value=str(rows["records"]),
                        anchor="composites (itemised 194 <= compiled 765, allowance 2000)")
total = sum(rows.values())
units = F(total, C); W = math.log2(float(units))
carried = total - rows["records"]
dbl = math.log2(float(F(total + carried, C))) - W
claim = F(int(round(math.ceil((W + dbl) * 1e4))), 10000)
gate = F(1168271, 10000)
rungs = [dict(nb=r["nb"], ell=r["ell"], n=r["trials"], h=r["hits"], q=r["q"],
              lo=r["cp95_lower"], ok=bool(r["cp95_lower"] >= 0.39)) for r in LAD["ladder"]]
rungs.append(dict(nb=22, ell=6944, n=LAD["primary"]["trials"], h=LAD["primary"]["hits"],
                  q=LAD["primary"]["q"], lo=LAD["primary"]["cp95_lower"],
                  ok=bool(LAD["primary"]["cp95_lower"] >= 0.39)))
claim_r = [r for r in rungs if r["ell"] >= 64]
def coch(rs):
    nt = sum(r["n"] for r in rs); ht = sum(r["h"] for r in rs); qb = ht / nt
    return sum(r["n"] * (r["q"] - qb) ** 2 / (qb * (1 - qb)) for r in rs), len(rs) - 1
Q, df = coch(claim_r); Qall, dfall = coch(rungs)
rng = random.Random(20261024); boot = []
for _ in range(2000):
    pick = [claim_r[rng.randrange(len(claim_r))] for _ in claim_r]
    boot.append(sum(r["h"] for r in pick) / sum(r["n"] for r in pick))
boot.sort()
succ = round(LAD["pooled_long_walks"]["cp95_lower"], 6)
kappa = max(CEN["charged_kappa"], D["kappa_carried"])
support = 1 << D["gap_support_log2"]; DIST = ELL - 1
abort = math.log2(DIST) - (support / MU) / math.log(2)
markov1 = math.log2(MU / support); markovu = markov1 + math.log2(DIST)
led = dict(version="v34", C_instructions_per_unit=C, rows=rows,
 row_log2={k: round(math.log2(v), 6) for k, v in rows.items() if v > 0},
 preprocessing_derivations=deriv, preprocessing_all_derived=True,
 derivation_inputs="design-v34.json + price-anchors-v34.json (compiled bodies); no row value is an input",
 total_instructions=total, W=str(units), W_log2=W,
 claim=dict(numerator=str(claim.numerator), denominator=str(claim.denominator),
            value=float(claim), reduced=math.gcd(claim.numerator, claim.denominator) == 1),
 claim_minus_W=round(float(claim) - W, 9), claim_covers_doubling_stress=bool(float(claim) - W >= dbl),
 gate=dict(numerator="1168271", denominator="10000"),
 gate_margin=round(float(gate) - float(claim), 6),
 claim_vs_gate=("OVER gate by %.4f bits -- DISCLOSED" % (float(claim) - float(gate))
                if float(claim) > float(gate) else "UNDER gate by %.4f" % (float(gate) - float(claim))),
 success=dict(claimed_value=succ, floor=D["success_floor"],
  key_path="structured-ladder.json:pooled_long_walks.cp95_lower",
  meets_floor=bool(succ >= D["success_floor"]),
  pool_level=dict(production_pool="ONE pool: nb=2^45 bases x one shared walk list of ell=2^70.650; "
                   "records = Cartesian product, so the trial unit IS the production unit",
   floor_level="asserted at the POOL level; reported per rung as well",
   per_rung=[dict(nb=r["nb"], ell=r["ell"], cp95_lower=round(r["lo"], 6), floor_met=r["ok"])
             for r in rungs],
   rungs_below_floor=sum(1 for r in rungs if not r["ok"]),
   heterogeneity=dict(Q=round(Q, 3), df=df, chi2_95=9.49, significant=bool(Q > 9.49),
                      all_six_Q=round(Qall, 3), all_six_driver="ell=16 short-walk control q=0.3718"),
   cluster_bootstrap_over_rungs=dict(resamples=2000, lo=round(boot[50], 6), hi=round(boot[1950], 6),
                                     clears_floor=bool(boot[50] >= 0.39)))),
 gap_support=dict(distinct_gaps=str(DIST), support_log2=D["gap_support_log2"],
  geometric_tail="HYPOTHESIS (mu=2^26.067), not a derived fact: the union bound 2^%.1f is its "
                 "CONSEQUENCE, not its evidence; no tail event beyond the observed maximum" % abort,
  markov="P(gap>t)<=E[gap]/t, no tail assumption: 2^%.1f per gap, 2^%.1f union over ell-1 -- "
         "INSUFFICIENT, stated as a failure" % (markov1, markovu),
  charge_disclosure="gap row 2^%.3f of 2^%.3f total: the tail moves the claim <1e-4 bits; the "
                    "claim does not rest on it" % (math.log2(rows["gap_point_tables"]), W)),
 kappa=dict(charged=kappa, cells={k: v["kappa"] for k, v in CEN["cells"].items()},
            rule=CEN["charged_kappa_rule"], v33_charged=2.613761,
            census_object="Cartesian base x walk, challenge multiplier in phi, both-side hist, "
                          "per-pool cluster bootstrap"),
 debit=dict(E_bad_log2=round(math.log2(kappa) + 2 * math.log2(N) - 168, 6), cap_bad_log2=84,
            form="kappa * N^2/2^168", key_width_bits=168),
 anchors=dict(file="price-anchors-v34.json", hprime_element_instructions=HP,
              chain_table_entry_instructions=CT,
              itemised_declared=ANC["itemised_sum_declared"],
              itemised_compiled=ANC["itemised_sum_compiled"],
              reconstructions={k: v["reconstructs"] for k, v in ANC["composites"].items()},
              chain="symbol -> clang -O2 -> otool -tV -> instructions -> price -> ledger row"),
 retune=dict(nb=NB, nb_log2=45, ell=ELL, ell_log2=round(math.log2(ELL), 6)))
json.dump(led, open(P("ledger-v34.json"), "w"), separators=(",", ":"))
assert total == sum(rows.values()) and float(claim) >= W
print(json.dumps(dict(W_log2=round(W, 9), claim=str(claim), success=succ, kappa=kappa,
                      Q=round(Q, 3), boot=(round(boot[50], 4), round(boot[1950], 4)),
                      hprime=HP, chain_entry=CT), indent=1))
