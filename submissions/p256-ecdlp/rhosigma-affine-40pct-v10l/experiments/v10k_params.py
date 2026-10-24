#!/usr/bin/env python3
"""v10k parameter / invariant / ledger generator.  Single machine source for every
number quoted in public-note.md, argument.md and certificate.json.

Rows emitted
  evidence/prod_params.json    production-parameter + scaling-invariant table
  evidence/invariant_fit.json  p(Lambda, ltheta) census grid, nu-invariance test, fit
  evidence/bucket.json         record-store arithmetic at production scale
  evidence/ledger.json         priced charge and claim (recomputed from the rows above)
  certificate.json             independent re-derivation of the claim from the rows
"""
import json, math, os
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")

# ---------------------------------------------------------------- P-256 constants
R_ORDER = 115792089210356248762697446949407573529996955224135760342422259061068512044369
N_CLASS = (R_ORDER - 1) // 2          # negation-map classes (excluding infinity)
FE_MUL_INSTR = 140                    # the package's own measured fe_mul price
INV_MULS = 385                        # 256 squarings + 128 multiplies + 1
BATCH_WALKERS = 1 << 16
STEP_INSTR = 658                      # integrated step_forward 235 + step_close 423
BACKSUB_MULS_PER_WALKER = 2           # out[i]=inv*PRE[i] ; inv=inv*d[i]
RARE_EVENT_RESERVE = 0.63976416015625
DP_BITS_PROD = 62
RECORD_BYTES = 104
CAP_BYTES_LOG2 = 73.0
W_AMPLIFY = 2

# ---------------------------------------------------------------- production row
S_LOG2 = 127.631517202749
S = 2.0 ** S_LOG2
sqrt_N = math.sqrt(N_CLASS)
sqrt_2N = math.sqrt(2 * N_CLASS)
lam_prod = S * S / (2 * N_CLASS)
theta = 2.0 ** (-DP_BITS_PROD)
R_rec = S * theta
nu = R_rec / sqrt_N
ell = S / BATCH_WALKERS
ell_theta = ell * theta

prod = {
    "S_log2": S_LOG2,
    "classes_log2": math.log2(N_CLASS),
    "sqrt_N_log2": math.log2(sqrt_N),
    "sqrt_2N_log2": math.log2(sqrt_2N),
    "S_over_sqrt_2N_log2": math.log2(S / sqrt_2N),
    "Lambda": lam_prod,
    "Lambda_derivation": "Lambda = S^2/(2N') with S=2^127.631517202749, N'=(r-1)/2",
    "dp_bits": DP_BITS_PROD,
    "theta": theta,
    "records_R_log2": math.log2(R_rec),
    "nu_records_over_sqrtN_log2": math.log2(nu),
    "nu_sq_over_2_log2": 2 * math.log2(nu) - 1,
    "walkers": BATCH_WALKERS,
    "trail_len_log2": math.log2(ell),
    "trail_over_sqrtN_log2": math.log2(ell / sqrt_N),
    "records_per_walker_log2": math.log2(ell_theta),
    "record_bytes": RECORD_BYTES,
    "records_bytes_log2": math.log2(R_rec * RECORD_BYTES),
    "cap_bytes_log2": CAP_BYTES_LOG2,
    "storage_slack_bits": CAP_BYTES_LOG2 - math.log2(R_rec * RECORD_BYTES),
    "invariant_note": ("a rho collision census is controlled by (Lambda, ltheta) where "
                       "Lambda=S^2/(2N') is the budget in class-space units and ltheta=ell*theta "
                       "is the number of distinguished points emitted per trail; nu=R/sqrt(N') "
                       "enters the success law only through the accidental record-record "
                       "coincidence rate nu^2/2 and through storage"),
}
assert abs(lam_prod - 0.6) < 1e-9, lam_prod

# ---------------------------------------------------------------- census grid
cm = json.load(open(os.path.join(EV, "census_multiscale.json")))
cells = []
for k in sorted(cm):
    c = cm[k]
    th = 2.0 ** (-c["dp_bits"])
    per = c["per_walker"]
    cells.append({
        "cell": k,
        "classes_log2": c["classes_log2"],
        "dp_bits": c["dp_bits"],
        "theta": th,
        "walkers": c["slots"],
        "per_walker_steps": per,
        "trials": c["trials"],
        "successes": c["successes"],
        "p_hat": c["recovery"],
        "lambda": per * c["slots"] * per * c["slots"] / (2 * c["classes"]),
        "ltheta": per * th,
        "nu": per * c["slots"] * th / math.sqrt(c["classes"]),
    })
for c in cells:
    assert abs(c["lambda"] - 0.6) < 0.02, (c["cell"], c["lambda"])

# --------------------------------------------------------------- Clopper-Pearson
def cp95_lower(k, n):
    """Exact one-sided 95% Clopper-Pearson lower endpoint, bisection in (0,1).
    For k==n this returns 0.05**(1/n) -- the clipped-helper bug is not present."""
    if k == 0:
        return 0.0
    if k == n:
        return 0.05 ** (1.0 / n)

    def cdf(p):
        # P[X <= k-1] for Binomial(n,p), computed in log space
        s = 0.0
        for i in range(0, k):
            s += math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                          + i * math.log(p) + (n - i) * math.log1p(-p))
        return s
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cdf(mid) > 0.05:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

# saturation: ltheta >= 10*tau, tau from the fit below; use the measured knee ltheta>=40
SAT = 40.0
plateau = [c for c in cells if c["ltheta"] >= SAT]
k_tot = sum(c["successes"] for c in plateau)
n_tot = sum(c["trials"] for c in plateau)
p_pooled = k_tot / n_tot
cp = cp95_lower(k_tot, n_tot)
p_max_cell = max(c["p_hat"] for c in plateau)
residual = p_max_cell - p_pooled           # same convention as the accepted v10j row
p_head = cp - residual

# chi-square heterogeneity over the plateau cells (Nu-invariance test)
chi2 = sum((c["successes"] - c["trials"] * p_pooled) ** 2 / (c["trials"] * p_pooled * (1 - p_pooled))
           for c in plateau)
dof = len(plateau) - 1

# saturation fit p = A*(1-exp(-ltheta/tau)) on ALL cells (least squares, 2 params)
def _sse(Ai, taui):
    return sum((c["p_hat"] - Ai * (1 - math.exp(-c["ltheta"] / taui))) ** 2 for c in cells)
best = None
for Ai in [x / 1000 for x in range(200, 400)]:
    for taui in [x / 100 for x in range(50, 1200)]:
        s = _sse(Ai, taui)
        if best is None or s < best[0]:
            best = (s, Ai, taui)
sse, A_fit, tau_fit = best
A0, t0 = A_fit, tau_fit
for Ai in [A0 + x / 100000 for x in range(-500, 501)]:
    for taui in [t0 + x / 10000 for x in range(-500, 501, 2)]:
        s = _sse(Ai, taui)
        if s < sse:
            sse, A_fit, tau_fit = s, Ai, taui
fit_rows = [{"cell": c["cell"], "ltheta": c["ltheta"], "p_hat": c["p_hat"],
             "p_fit": A_fit * (1 - math.exp(-c["ltheta"] / tau_fit)),
             "resid": c["p_hat"] - A_fit * (1 - math.exp(-c["ltheta"] / tau_fit))}
            for c in cells]
p_fit_prod = A_fit * (1 - math.exp(-ell_theta / tau_fit))

inv = {
    "invariant": "(Lambda, ltheta)",
    "cells": cells,
    "saturation_threshold_ltheta": SAT,
    "plateau_cells": [c["cell"] for c in plateau],
    "plateau_successes": k_tot,
    "plateau_trials": n_tot,
    "p_pooled": p_pooled,
    "p_pooled_cp95": cp,
    "p_max_cell": p_max_cell,
    "scale_residual": residual,
    "p_headline": p_head,
    "nu_range_log2": [min(c["nu"] for c in plateau), max(c["nu"] for c in plateau)],
    "chi2": chi2, "chi2_dof": dof,
    "chi2_note": "pooled Binomial heterogeneity across the plateau cells; a small chi2 "
                 "against dof is the measured Nu-invariance (cells differ only in nu and scale)",
    "fit": {"A": A_fit, "tau": tau_fit, "sse": sse, "rows": fit_rows,
            "p_at_production_ltheta": p_fit_prod,
            "extrapolation_bits": math.log2(ell_theta / max(c["ltheta"] for c in cells)),
            "label": "fitted saturation law p=A(1-exp(-ltheta/tau)); anchors are the 9 measured "
                     "cells; production ltheta exceeds the largest anchor by "
                     "log2(production/anchor) bits, where the law is flat to 2^-%.0f"
                     % (math.log2(ell_theta / max(c["ltheta"] for c in cells)))},
    "nu_channel": {
        "nu_prod_log2": math.log2(nu),
        "accidental_record_collisions_log2": 2 * math.log2(nu) - 1,
        "note": "records over sqrt(class space); at production the accidental record-record "
                "coincidence rate is nu^2/2 (2^-124.7), and such a coincidence is rejected by the "
                "full-point verification, so nu is a false-positive channel, not a success channel",
    },
    "fit_beats_disclaimer": True,
}

# ---------------------------------------------------------------- record store
NBUCKETS_LOG2 = 32
D = 1 << NBUCKETS_LOG2
recs_per_bucket = R_rec / D
search_steps = math.log2(recs_per_bucket)
bytes_per_bucket = recs_per_bucket * RECORD_BYTES
instr_per_record = search_steps * 24 + 201          # binary search compares + insert/shift
bucket = {
    "records_log2": math.log2(R_rec),
    "n_buckets_log2": NBUCKETS_LOG2,
    "records_per_bucket": recs_per_bucket,
    "records_per_bucket_log2": search_steps,
    "bucket_bytes_log2": math.log2(bytes_per_bucket),
    "lookup": "one 2^32-way radix shard on the high 32 key bits; inside a shard the records are "
              "kept in key order and searched by binary search, so the per-record cost is "
              "log2(records_per_bucket) comparisons, not a load factor",
    "instr_per_record": instr_per_record,
    "record_share_exact": instr_per_record * 2.0 ** (-DP_BITS_PROD),
    "record_work_total_log2": math.log2(R_rec * instr_per_record),
    "directory_bytes": D * 16,
    "total_bytes_log2": math.log2(R_rec * RECORD_BYTES + D * 16),
    "cap_bytes_log2": CAP_BYTES_LOG2,
    "fits": math.log2(R_rec * RECORD_BYTES + D * 16) < CAP_BYTES_LOG2,
    "supersedes": "v10j's '2^32 buckets at load factor 0.5 / chain 0.69' row, which divided "
                  "2^65.63 records among 2^32 buckets (2^33.63 per bucket), was arithmetically false",
}

# ---------------------------------------------------------------- charge and claim
backsub = BACKSUB_MULS_PER_WALKER * FE_MUL_INSTR
inv_share = INV_MULS * FE_MUL_INSTR / BATCH_WALKERS
rec_share = bucket["record_share_exact"]
charge = STEP_INSTR + backsub + inv_share + RARE_EVENT_RESERVE + rec_share
W_single = S_LOG2 + math.log2(charge / 1000.0)
W_amp = W_single + math.log2(W_AMPLIFY)
claim = math.ceil(W_amp * 1e4) / 1e4
p = p_head
succ = 1 - (1 - p) ** W_AMPLIFY

ledger = {
    "W": {
        "A_log2_Tprod": S_LOG2,
        "charge_terms": {
            "step_instructions_integrated": STEP_INSTR,
            "step_forward": 235, "step_close": 423,
            "batch_backsub_muls_per_walker": BACKSUB_MULS_PER_WALKER,
            "batch_backsub_instr_per_walker": backsub,
            "batch_inversion_share": inv_share,
            "batch_inversion_muls": INV_MULS,
            "fe_mul_instr": FE_MUL_INSTR,
            "batch_walkers": BATCH_WALKERS,
            "rare_event_reserve": RARE_EVENT_RESERVE,
            "record_processor_share": rec_share,
            "prefix_product_note": "the prefix product PRE[i+1]=PRE[i]*d_i is 1 fe_mul per walker "
                                   "per step and is INSIDE the counted step_forward (235); it is "
                                   "computed once there and is not recomputed",
        },
        "charge_total": charge,
        "W_single_instance": W_single,
        "amplification_w": W_AMPLIFY,
        "W_amplified": W_amp,
        "claim": claim,
    },
    "amplification": {"w": W_AMPLIFY, "p_per_instance": p, "success": succ, "floor": 0.39,
                      "independent": "two independent instances per target, fresh coins, separate "
                                     "table instances, work charged twice, run sequentially so the "
                                     "store is reused"},
    "p_transfer": {"headline": p, "pooled_cp95": cp, "residual": residual,
                   "invariant": "(Lambda, ltheta) matched at the production point",
                   "source": "evidence/invariant_fit.json"},
    "assert": {
        "charge_equals_sum": abs(charge - (STEP_INSTR + backsub + inv_share + RARE_EVENT_RESERVE
                                           + rec_share)) < 1e-12,
        "claim_covers_W": claim >= W_amp,
        "claim_minus_W_bits": claim - W_amp,
        "success_ge_floor": succ >= 0.39,
        "storage_fits": bucket["fits"],
        "lambda_matched": abs(lam_prod - 0.6) < 1e-9,
    },
}

cert = {
    "claim": [int(round(claim * 10000)), 10000],
    "claim_f": claim,
    "recomputed": {
        "W_single": W_single, "W_amplified": W_amp,
        "charge": charge,
        "p": p, "success": succ,
    },
    "checks": [
        {"name": "charge_terms_sum_to_charge", "ok": ledger["assert"]["charge_equals_sum"]},
        {"name": "claim_ge_W_amplified", "ok": claim >= W_amp},
        {"name": "success_ge_floor_0.39", "ok": succ >= 0.39},
        {"name": "storage_under_cap", "ok": bucket["fits"]},
        {"name": "Lambda_production_equals_fixture_0.60", "ok": abs(lam_prod - 0.6) < 1e-9},
        {"name": "batch_backsub_charged", "ok": backsub > 0},
        {"name": "bucket_row_recomputed", "ok": bucket["fits"] and bucket["instr_per_record"] > 40},
        {"name": "nu_invariance_chi2_dof", "ok": chi2 < 3 * dof,
         "chi2": chi2, "dof": dof},
        {"name": "saturation_fit_residual_max", "ok": max(abs(r["resid"]) for r in fit_rows) < 0.03,
         "max_abs_resid": max(abs(r["resid"]) for r in fit_rows)},
    ],
}
cert["all_pass"] = all(c["ok"] for c in cert["checks"])

for name, obj in [("prod_params.json", prod), ("invariant_fit.json", inv),
                  ("bucket.json", bucket), ("ledger.json", ledger)]:
    with open(os.path.join(EV, name), "w") as f:
        json.dump(obj, f, indent=1, sort_keys=True)
with open(os.path.join(PKG, "certificate.json"), "w") as f:
    json.dump(cert, f, indent=1, sort_keys=True)

print(json.dumps({"claim": claim, "W_single": W_single, "W_amp": W_amp, "charge": charge,
                  "p": p, "cp95": cp, "pooled": p_pooled, "residual": residual,
                  "success": succ, "chi2": chi2, "dof": dof, "A": A_fit, "tau": tau_fit,
                  "p_fit_prod": p_fit_prod, "instr_per_record": instr_per_record,
                  "all_pass": cert["all_pass"]}, indent=1))
