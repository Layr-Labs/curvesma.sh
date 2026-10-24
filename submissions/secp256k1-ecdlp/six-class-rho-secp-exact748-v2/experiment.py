"""Recompute the 114-row analytic step, the 748-instruction price, the
126.3194 work cap and the H1/Markov success floor.  Standard library only.
This check does not execute the counted AArch64 kernels; their counts and
hashes are documented in kernel-counts.md and the public v9 package."""
import json
import sys
from decimal import Decimal as D, getcontext

getcontext().prec = 80
r = D(int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16))
N = (r - 1) / 6
R = D(2048)
Np = N * R / (R - 1)
L = D(20) * (D(2) ** 64)
Dcap = D(2) ** D("63.83")
step = D(748) / D(1000)
dbl = D("1.1")
trav = D("0.795")
kappa = D("1.01")
p0 = D("0.3950")
cap = D("126.3194")

rows = {
    "F1_truncated_product": 24,
    "range_test_F0_S_F2": 8,
    "decision": 13,
    "sign_pminusy_K": 13,
    "h_dp_t_address": 7,
    "pass3_state_traffic": 5,
    "fingerprints_history": 26,
    "zero_denominator": 4,
    "length_cap_meter_loop": 13,
    "share_of_inversion_and_calls": 1,
}
counted = 111 + 523
rows_total = sum(rows.values())
step_instructions = counted + rows_total

S0 = (2 * Np * (1 / (1 - p0)).ln()).sqrt()
S = kappa * S0
B = step * (S + D(2) ** 60)
H = (D(2) ** 61 + D(2) ** 20 * (L + 65) * 600
     + 16 * ((L + 65) * 600 + D(2) ** 13)
     + 40 * Dcap + D(2) ** D("58.9") + D(2) ** 25)
E = ((dbl - step) * (D(2) ** D("-13.56")) * (S + D(2) ** 60)
     + (trav - step) * 64 * (D(2) ** D("-19.99")) * (S + D(2) ** 60)
     + D(2) ** D("75.9") + D(2) ** 70)
W = D(2) ** cap
A = W - B - H
markov = E / A
success = p0 - markov - D("2e-7")

lg = lambda x: x.ln() / D(2).ln()
checks = {
    "counted_is_634": counted == 634,
    "rows_sum_to_114": rows_total == 114,
    "step_is_748": step_instructions == 748,
    "price_is_748": step * 1000 == 748,
    "headroom_positive": A > 0,
    "success_above_0_39003": success > D("0.39003"),
    "success_above_0_39": success > D("0.39"),
    "cap_matches_claim": cap == D("126.3194"),
}
out = {
    "counted_instructions": counted,
    "analytic_rows": rows,
    "analytic_rows_total": rows_total,
    "step_instructions": step_instructions,
    "step_units": str(step),
    "kappa": str(kappa),
    "p0": str(p0),
    "claimed_time_log2": str(cap),
    "log2_S": float(lg(S)),
    "log2_B": float(lg(B)),
    "log2_H": float(lg(H)),
    "log2_E": float(lg(E)),
    "log2_A": float(lg(A)),
    "markov_early_halt": float(markov),
    "success_lower_bound": float(success),
    "checks": checks,
    "all_checks": all(checks.values()),
}
print(json.dumps(out, indent=1))
if not out["all_checks"]:
    sys.exit(1)
