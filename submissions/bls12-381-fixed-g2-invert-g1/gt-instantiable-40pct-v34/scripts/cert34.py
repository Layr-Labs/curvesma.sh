#!/usr/bin/env python3
"""cert34.py -- independent re-derivation of every charged number in v34.

Prints one JSON line. Each check recomputes from evidence/*.json; none reads a
scalar out of submission.json. Exit 1 on any failure.
"""
import json, math, os, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
L = json.load(open(os.path.join(EV, "ledger34.json")))
ARM = json.load(open(os.path.join(EV, "arm34.json")))
TRAP = json.load(open(os.path.join(EV, "trap34.json")))
SUB = json.load(open(os.path.join(PKG, "submission.json")))
ok = {}
rows = {k: int(v) for k, v in L["rows_ir"].items()}

# 1. row sum == W
W = sum(rows.values())
ok["row_sum_equals_W"] = (W == int(L["W_total_ir"]))

# 2. claim is the minimal n, exact integer test, and n-1 fails
n = int(L["claim"]["numerator"])
ok["claim_exact"] = (pow(20000, 10000) * (1 << n) >= pow(W, 10000))
ok["claim_minimal"] = (pow(20000, 10000) * (1 << (n - 1)) < pow(W, 10000))

# 3. reduced rational, gcd 1, and submission.json carries the same rational
g = math.gcd(n, 10000)
ok["reduced"] = (L["claim"]["reduced"] == "%d/%d" % (n // g, 10000 // g))
ok["sub_score_matches"] = (SUB["claims"][0]["score"]["numerator"] == str(n)
                           and SUB["claims"][0]["score"]["denominator"] == "10000")

# 4. the fantasy rows are gone
ok["no_route_a_init_row"] = ("store_init_route_a" not in rows)
blob = SUB["claims"][0]["statement"]
ok["no_2_135_footprint"] = all("peak_bytes" not in json.dumps(v) for v in [L]) and "bucket_array_bytes" not in json.dumps(L)
ok["no_1490e39"] = not any("init" in k for k in rows) and not any("store_init" in k for k in rows)
ok["no_2_120_token"] = ("2^120" not in blob)

# 5. coefficient row is a ceil of the measured span, and 179 is not cited
MUL = int(json.load(open(os.path.join(EV, "prices28.json")))["MUL_CHARGE"])
fp = -(-MUL // 54)
S = int(L["protocol"]["steps"])
ok["coeff_row_ceiled"] = (rows["scalar_coefficient_update"] == (2 * fp + 12) * S)
ok["no_179_row"] = all(str(v) != "179" and not str(v).startswith("179") for v in rows.values())

# 6. trap event is the census predicate and the charge is the exact integer
trap_at_S = -(-TRAP["measured_n"] * S // TRAP["measured_d"])
ok["trap_count"] = (str(trap_at_S) == TRAP["trap_at_S"])
ok["trap_charge"] = (int(TRAP["charge_ir"]) == trap_at_S * 32 * 48)
ok["trap_event_is_union"] = ("br OR guard" in TRAP["event"])

# 7. retention pair: derived recomputed from lambda, measured reported with se
lam = Fraction(int(ARM["lambda_num"]), int(ARM["lambda_den"]))
der = sum(math.exp(-float(lam)) * float(lam) ** k / math.factorial(k) for k in range(5))  # = 1 - P(X>=5)
ok["retention_derived_matches"] = (abs(der - ARM["p_derived"]) < 1e-9)
ok["measured_arm_trials"] = (ARM["trials"] == 1485)
ok["lambda_int_error_small"] = (abs(Fraction(ARM["arrivals"], ARM["nb"]) - lam) / lam < 1e-3)
ok["wrong_lambda_control_does_not_match"] = (not ARM["control_matches"]
                                             and abs(ARM["control_2417"] - ARM["p_hat"]) > 3 * ARM["se"])

# 8. every evidence file named by the claim exists
missing = [e for e in SUB["claims"][0]["evidence"] if not os.path.exists(os.path.join(PKG, e))]
ok["evidence_present"] = (missing == [])

print(json.dumps({"ok": all(ok.values()), "checks": ok, "claim": str(n), "W_log2": math.log2(W),
                  "coeff_per_step": 2 * fp + 12, "missing": missing}, default=str))
sys.exit(0 if all(ok.values()) else 1)
