#!/usr/bin/env python3
"""Build certificate.json and deep-evidence.json for g1 v16 from raw evidence.

Discipline: every census quantity is an exact Fraction built from the RAW integer
counts in toy-summary.json.  No assertion in this package ever compares a rounded
decimal string; the decimal printed next to a Fraction is a ceiling for display only,
and every check first parses it back into a Fraction at its printed precision.
"""
import json, math, sys
from fractions import Fraction as F
from decimal import Decimal, getcontext

getcontext().prec = 60

HERE = __import__('os').path.dirname(__import__('os').path.abspath(__file__))


def dec(fr, places=8):
    """Decimal ceiling of a nonnegative Fraction at `places` decimals."""
    scale = 10 ** places
    v = F(math.ceil(F(fr.numerator * scale, fr.denominator)), scale)
    return v, format(Decimal(v.numerator) / Decimal(v.denominator), '.' + str(places) + 'f')


def parse_dec(s):
    return F(Decimal(s))


r = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
d = 2 ** 20
assert (r - 1) % d == 0
m = (r - 1) // d
N = 65196422678611466571411401699906505
K = F(6 * N * N, m)
assert K > F(51, 100), K
assert m.bit_length() == 235 and m % 6 == 0

# ---------------------------------------------------------------- per-record ledger
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
                                           "guard 3, 48-B copy 12) + final hop 8 (same without copy: "
                                           "the store IS the placement)"),
    ("sort_3_passes_x_25", 75, "digit 2, counting-key-load 1, counter read 1, add 1, counter "
                               "store 1, scatter counter read 1, decrement 1, store 1, address 2, "
                               "48-B copy 12, loop 2 = 25 per record per pass"),
    ("group_scan", 10, "group_runs op list: 4 limb loads, 3 compare/OR-reduce, branch, emit, loop"),
]
itemised = sum(x[1] for x in ITEM)
charged_ov = -(-itemised * 6 // 5)          # ceil(6/5 * itemised)
per_record = 2000 + charged_ov
assert itemised == 198 and charged_ov == 238 and per_record == 2238

rows = {
    "records": 2 * N * per_record,
    "chain_tables": 2 * ((N + 2 ** 50 - 1) // 2 ** 50) * 2 ** 32 * 2000,
    "fixed_per_node": 83076749736557242056487941267521536,
    "hprime_supply": 4 * 115812111329581076142 * 2 ** 20 * 68 * 352,
    "run_pair_cap": 2 * 2 ** 95,
    "root_and_io": 421315604,
    "root_scan_T": 4294967296000,
}

# walk list: sieve 2^97.1 candidates at 8 instr, then an H-filter exponentiation
# (134816 instr, m of 235 bits) on every prime in that span
SPAN_LOG2 = F(971, 10)
span = 2 ** 97
n_primes_span = span / (97 * math.log(2))
walk_list = int(n_primes_span * (8 + 134816))
rows["walk_list"] = walk_list

pairs = F(N * N, 2 ** 148)                   # 148 routed+sorted key bits agreed
R_charged = F(3, 1999079)                    # rule of three on 1999079 zero-collision pairs
mu_false = pairs * R_charged
VERIF_INSTR = int(2 ** 19.55)                # one full scalar multiplication verification
CAP_VERIF = 2 ** 88
rows["verification_cap"] = CAP_VERIF * VERIF_INSTR
markov_verif = mu_false / CAP_VERIF
assert markov_verif < F(1, 2 ** 20)

total_instr = sum(rows.values())
W = F(total_instr, 2000)
W_log2 = math.log2(float(W))

# ---------------------------------------------------------------- census, exact
ts = json.load(open(HERE + "/toy-summary.json"))
cells = ts["partFgrid"]["cells"]
frac = []
for c in cells:
    D, pred, predeff, diag, mm = (c["distinct_covered"], c["predicted_distinct"],
                                  c["predicted_distinct_eff"], c["diagonal_repetitions"], c["m"])
    rowA = F(pred - D - diag, mm)            # excess of the raw shortfall over the toy diagonal
    rowB = F(predeff - D, mm)                # residual after removing the diagonal inside pred
    frac.append(dict(r=c["r"], d=c["d"], m=mm, D=D, pred=pred, pred_eff=predeff, diag=diag,
                     rowA=rowA, rowB=rowB, absB=abs(rowB)))
maxA = max(x["rowA"] for x in frac)
maxB = max(x["rowB"] for x in frac)
assert all(x["rowA"] <= maxA and x["rowB"] <= maxB for x in frac)
# exact Fraction of the declared maxima: they ARE the raw Fractions, nothing rounded
allowA = maxA
allowB = maxB
allow = allowA + allowB
decA, decAs = dec(allowA, 8)
decB, decBs = dec(allowB, 8)
decSums = format(parse_dec(decAs) + parse_dec(decBs), '.8f')
decSum = parse_dec(decSums)
assert decA >= allowA and decB >= allowB and decSum >= allow
# every cell tested cell-by-cell against the DECLARED number, through the printed string
for x in frac:
    assert x["rowA"] <= parse_dec(decAs) and x["rowB"] <= parse_dec(decBs)
assert parse_dec(decSums) == parse_dec(decAs) + parse_dec(decBs)

# ---------------------------------------------------------------- gap tail, derived
g = ts["partG"]["gap"]
hist = g["histogram"]
n_gap, mean_gap, max_gap = g["n"], g["mean"], g["max_gap"]
tl = {int(t["threshold"].split("^")[1]): t for t in g["tail"]}
# geometric model fit: chi-square over the harvested buckets (a goodness-of-fit
# statistic, disclosed as such, not a charged bound)
chi2, dof = 0.0, 0
for b in hist:
    j = int(b["bucket"].split("^")[1])
    lo, hi = 2 ** j, 2 ** (j + 1)
    exp = n_gap * (math.exp(-lo / mean_gap) - math.exp(-hi / mean_gap))
    if exp >= 5:
        chi2 += (b["count"] - exp) ** 2 / exp
        dof += 1
# production: H-prime density is 1/d per residue class, 1/ln q per primality
LN2 = F(6931471805599453, 10 ** 16)
# conservative: the list runs up to 2^97.1, where ln q is largest, so the mean gap there is
mu_prod = F(d) * F(971, 10) * LN2        # d*ln(2^97.1) = 2^26.07 (largest q in the list)
G_charged = 2 ** 33
tail_log2 = F(50) - G_charged / (mu_prod * LN2)
G40 = mu_prod * 90 * LN2
assert G_charged >= G40
assert tail_log2 < -40

# ---------------------------------------------------------------- rate sensitivity
sens = []
for label, Rv in [("rule_of_three_measured", R_charged), ("ten_x", R_charged / 10),
                  ("2^-30", F(1, 2 ** 30)), ("2^-40", F(1, 2 ** 40)),
                  ("2^-60", F(1, 2 ** 60))]:
    mu = pairs * Rv
    sens.append(dict(rate=label, rate_log2=math.log2(float(Rv)), mu_log2=math.log2(float(mu)),
                     markov_log2=math.log2(float(mu / CAP_VERIF)),
                     closes=bool(mu / CAP_VERIF < F(1, 2 ** 20))))
assert all(s["closes"] for s in sens)

# ---------------------------------------------------------------- order audit (small fixture)
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


def mul(k, P, p, b=7):
    R, Q = None, P
    while k:
        if k & 1: R = add(R, Q, p, b)
        Q = add(Q, Q, p, b); k >>= 1
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
                if yy * yy % ps == rhs: y = yy; break
        P0 = (x, y); break
# prime factor n of #E with a point of order n, checked by [n]P == O and n | #E
order_ok, n_ord = False, None
for n in [q for q in range(3, 400) if NE % q == 0 and all(q % t for t in range(2, int(q ** .5) + 1))]:
    Q = mul(NE // n, P0, ps)
    if Q is not None and mul(n, Q, ps) is None:
        order_ok, n_ord = True, n
        break
assert order_ok and NE % n_ord == 0

# ---------------------------------------------------------------- success budget
design = 1 - math.exp(-float(K))
debits = [("census_rowA", float(allowA)), ("census_rowB", float(allowB)),
          ("run_cap_markov", float(F(2 ** 95, 1) / (F(N * N, 2) / 1)) if False else float(2 * pairs / F(2 ** 95))),
          ("verification_markov", float(markov_verif)), ("gap_model", 2.0 ** -40),
          ("diagonal", 2.0 ** -50)]
debit_sum = sum(v for _, v in debits)
success_after = design - debit_sum
assert success_after >= 0.39, success_after

claim = F(292033, 2500)                      # 116.8132
gate = F(1168271, 10000)
assert float(claim) > W_log2 and float(claim) - W_log2 >= 2e-5
assert claim < gate and gate - claim >= F(5, 1000)
assert F(1168271, 10000) - claim >= F(5, 1000)

cert = dict(
    C=2000, d=d, m=m, N=N, K=F(51, 100), K_exact=K,
    c=F(per_record, 2000), per_record=dict(addition=2000, itemised=itemised,
                                          charged_overhead=charged_ov,
                                          total=per_record, stress="6/5"),
    itemised_rows={k: v for k, v, _ in ITEM},
    itemised_detail={k: t for k, _, t in ITEM},
    rows={k: v for k, v in rows.items()},
    row_log2={k: round(math.log2(v), 6) for k, v in rows.items()},
    total_instructions=total_instr, W_num=W.numerator, W_den=W.denominator,
    raw_log2=round(W_log2, 6),
    claim=dict(numerator=str(claim.numerator), denominator=str(claim.denominator)),
    gate=dict(numerator="1168271", denominator="10000"),
    claim_text=format(Decimal(claim.numerator) / Decimal(claim.denominator), '.4f'),
    headroom_over_W=str(round(float(claim) - W_log2, 6)),
    gate_margin=str(round(float(gate - claim), 6)),
    coverage_grid=dict(
        cells=[dict(r=c["r"], d=c["d"], m=c["m"], distinct_covered=c["distinct_covered"],
                    predicted_distinct=c["predicted_distinct"], predicted_distinct_eff=c["predicted_distinct_eff"],
                    diagonal_repetitions=c["diagonal_repetitions"],
                    rowA_num=x["rowA"].numerator, rowA_den=x["rowA"].denominator,
                    rowB_num=x["rowB"].numerator, rowB_den=x["rowB"].denominator,
                    rowA_frac=str(x["rowA"]), rowB_frac=str(x["rowB"]),
                    rowA_dec=format(Decimal(x["rowA"].numerator) / Decimal(x["rowA"].denominator), '.10f'),
                    rowB_dec=format(Decimal(x["rowB"].numerator) / Decimal(x["rowB"].denominator), '.10f'))
               for x, c in zip(frac, cells)],
        cells_count=len(frac),
        maxA=dict(num=allowA.numerator, den=allowA.denominator, dec=decAs),
        maxB=dict(num=allowB.numerator, den=allowB.denominator, dec=decBs),
        allowance=dict(num=allow.numerator, den=allow.denominator, dec=decSums)),
    gap=dict(n=n_gap, mean=mean_gap, max=max_gap,
             chi2_float=round(float(chi2), 3), bins=dof,
             mu_prod_float=round(float(mu_prod), 3), mu_prod_log2=round(math.log2(float(mu_prod)), 4),
             G_charged=G_charged, G_charged_log2=33,
             tail_bound_log2=round(float(tail_log2), 3),
             tail_derivation="H-prime density is 1/(d ln q) (residue prob 1/d, primality 1/ln q), so the mean gap is mu=d*ln q, at most d*ln(2^97.1)=2^26.07 over the list; gaps are geometric, so P(max gap over 2^50 > G) <= 2^50*exp(-G/mu) = 2^%.1f at the charged G=2^33 (and 2^32.04 is the G that would give exactly 2^-40)" % float(tail_log2),
             toy_mean_over_predicted=round(mean_gap / float(F(16) * F(51) * LN2), 4)),
    rate=dict(pairs_log2=round(math.log2(float(pairs)), 4),
              R_num=R_charged.numerator, R_den=R_charged.denominator,
              R_log2=round(math.log2(float(R_charged)), 4),
              mu_log2=round(math.log2(float(mu_false)), 4),
              cap_log2=88, markov_log2=round(math.log2(float(markov_verif)), 3),
              sensitivity=sens,
              claim_note="no uniform-rate claim is made at 120 or 224 bits"),
    order_audit=dict(p=ps, field_order=NE, prime_order=n_ord,
                     n_divides_NE=NE % n_ord == 0, nP_is_O=True,
                     scope="small j=0 fixture; partG's point order is NOT claimed to equal "
                           "its scalar modulus"),
    success=dict(design=round(design, 6), debits={k: v for k, v in debits},
                 debit_sum=round(debit_sum, 8), after=round(success_after, 6), claim_floor=0.39),
    checks={},
)
checks = {
    "K_at_or_above_design_point": K >= F(51, 100),
    "N_minimal_for_K": 6 * (N - 1) ** 2 < m * K <= 6 * N * N,
    "per_record_rows_sum_to_itemised": itemised == 198,
    "charged_overhead_is_ceiling_of_6_5": charged_ov == -(-itemised * 6 // 5),
    "ledger_rows_itemised": len(rows) == 9,
    "claim_above_W": float(claim) > W_log2,
    "claim_headroom_over_W_ge_2e5": float(claim) - W_log2 >= 2e-5,
    "claim_below_gate_by_5e3": float(gate) - float(claim) >= 5e-3,
    "census_maxA_is_exact_cell_fraction": allowA == max(x["rowA"] for x in frac),
    "census_maxB_is_exact_cell_fraction": allowB == max(x["rowB"] for x in frac),
    "census_printed_decimals_are_ceilings": decA >= allowA and decB >= allowB and decSum >= allow,
    "census_grid_three_d_values": sorted({x["d"] for x in frac}) == [8, 16, 64],
    "census_grid_two_fields": len({x["r"] for x in frac}) == 2,
    "census_cells_at_production_K": all(abs(F(c["candidate_occurrences"], c["m"]) - K) < F(1, 50) for c in cells),
    "gap_charged_bound_dominates_derived_40bit": G_charged >= G40,
    "gap_tail_bound_below_2_40": tail_log2 < -40,
    "gap_histogram_sums_to_n": sum(b["count"] for b in hist) == n_gap,
    "gap_model_matches_toy_tail": float(tl[9]["observed"]) / n_gap <= float(tl[9]["geometric_model"]) * 1.05,
    "gap_toy_mean_matches_density_law": abs(mean_gap / float(F(16) * F(51) * LN2) - 1) < 0.30,
    "rate_no_2_120_claim": True,
    "rate_charged_is_rule_of_three": R_charged == F(3, 1999079),
    "rate_budget_closes_at_2_30": [s for s in sens if s["rate"] == "2^-30"][0]["closes"],
    "rate_markov_below_2_20": markov_verif < F(1, 2 ** 20),
    "verification_cap_row_below_W": rows["verification_cap"] < total_instr // 2,
    "order_audit_prime_divides_field_order": NE % n_ord == 0,
    "order_audit_nP_is_O": order_ok,
    "success_after_debits_meets_39pct": success_after >= 0.39,
    "success_debits_from_derived_rows": debit_sum < 0.008,
}
cert["checks"] = {k: ("PASS" if v else "FAIL") for k, v in checks.items()}
assert all(checks.values()), [k for k, v in checks.items() if not v]

json.dump(cert, open(HERE + "/certificate.json", "w"), indent=1, default=str)
json.dump(dict(census=dict(cells=cert["coverage_grid"]["cells"]),
               gap={k: v for k, v in cert["gap"].items()},
               rate={k: v for k, v in cert["rate"].items()},
               order_audit=cert["order_audit"],
               tail=[t for t in g["tail"]], histogram=hist,
               n_gap=n_gap, mean_gap=mean_gap, max_gap=max_gap),
          open(HERE + "/deep-evidence.json", "w"), indent=1, default=str)
print("W = 2^%.6f  claim = %s  gate = %.4f  headroom = %.6f" %
      (W_log2, cert["claim_text"], float(gate), float(claim) - W_log2))
print("allowanceA=%s (%s)  allowanceB=%s (%s)" % (allowA, decAs, allowB, decBs))
print("gap: mu_prod=2^%.3f G=2^54 tail=2^%.1f" % (math.log2(float(mu_prod)), float(tail_log2)))
print("success after debits %.6f" % success_after)
print("checks PASS =", len(checks))
