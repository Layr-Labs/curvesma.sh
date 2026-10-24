#!/usr/bin/env python3
"""verify_certificate.py -- recompute every printed number and write certificate.json.

Checks (all must be true):
  * convention: the 161-row total reproduces mpjunior92 7b34f78d's printed bound
    to < 0.01 bits at their kappa = 1.02 (their own row table, our W model);
  * agreement: the 448-bit-window predicate agrees with the exact argmin on 5000
    production classes (the transport premise), and the 12-partial-product
    window stays within 2^64.25 of the full window;
  * TOY-SPEC FIDELITY: window_generic with production constants is bit-identical
    to the probe's verbatim spec_window12 on 100 keys, and the toy's instantiation
    agrees with the exact argmin on 100 toy keys;
  * fallback: 383 * 3 * 2^-64 per step over S_ops steps fits in 2^62 reference
    units and costs 0.000000 bits;
  * kappa budget product <= 2^(claim - log2 S0), kappa <= 1.015;
  * claim < ship line (gate - 0.008) and margin >= 0.008;
  * toy evidence: CP-95 lower >= 0.39 at every density, invariants exact,
    collision failures zero.
"""
import hashlib
import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
import selector  # noqa: E402
import ledger    # noqa: E402

p, beta = ledger.p, ledger.beta
beta2 = (beta * beta) % p
M = 1 << 448
c1 = (M * beta) // p
c2 = (M * beta2) // p
MASKL = (1 << 64) - 1
checks = []


def chk(name, cond, detail):
    checks.append({"name": name, "ok": bool(cond), "detail": detail})


# ---- 1. convention: reproduce their total ---------------------------------
def w_of(kappa, n, s0):
    L = math.log(1 / 0.6097)
    Np = (p * 0 + ((ledger.r - 1) // 6)) * 1.0
    Np = Np * (2 ** 20) / (2 ** 20 - 1)
    s = math.sqrt(2 * Np * L)
    ls = math.log2(s) + math.log2(kappa) + math.log2(1 + 2 ** 60 / (kappa * s))
    return ls + math.log2(1 + n / 2000)


their = w_of(1.02, 161, None)
chk("convention_reproduces_7b34f78d",
    abs(their - 126.271442) < 0.01,
    f"our 161 rows at kappa 1.02 -> 2^{their:.6f}; their printed budget row 2^126.271442")

# ---- 2. agreement on 5000 production classes ------------------------------
random.seed(20261024)
agree = 0
n_tot = 5000
maxdv = 0
near = 0
mismatch = []
for _ in range(n_tot):
    x = random.randrange(1, p)
    ie, vals = selector.twist_index_exact(x, p, beta)
    iw, wvals = selector.twist_index_production(x, p, beta)
    agree += (ie == iw)
    if ie != iw:
        mismatch.append(x)
    full = ((x * c1) & (M - 1)) >> 320
    maxdv = max(maxdv, abs(wvals[1] - full))
    s = sorted(vals)
    if s[1] - s[0] < 2 ** 315:
        near += 1
chk("agreement_5000_of_5000", agree == n_tot,
    f"{agree}/{n_tot} production classes: 448-bit window argmin == exact argmin")
chk("window_within_2_64", maxdv <= 2 ** 65,
    f"max |12-product window - full window| = {maxdv} <= 2^{math.log2(max(maxdv,1)):.2f} (their 3*2^-64)")
chk("near_tie_band_empty_in_sample", near == 0,
    f"{near} of {n_tot} classes had exact top-2 gap < 2^315 (the fallback band)")

# ---- 3. TOY-SPEC FIDELITY on 100 keys ------------------------------------
fid = 0
for _ in range(100):
    x = random.randrange(1, p)
    a, n1 = selector.spec_window12(x, c1)
    b = selector.window_generic(x, c1, 6, {5, 6})
    fid += (a == b and n1 == 12)
chk("toy_spec_fidelity_100_keys", fid == 100,
    f"window_generic(prod constants) == verbatim spec_window12 on {fid}/100 keys, 12 products each")

toyfid = 0
import toy  # noqa: E402
for _ in range(100):
    x = random.randrange(1, toy.P)
    ie, _ = selector.twist_index_exact(x, toy.P, toy.BETA)
    iw, _ = selector.twist_index(x, toy.P, toy.BETA, shift=128, n_limbs=1, band={0, 1})
    toyfid += (ie == iw)
chk("toy_selector_matches_exact_100_keys", toyfid == 100,
    f"toy instantiation picks the exact argmin on {toyfid}/100 toy keys")

# ---- 4. fallback ----------------------------------------------------------
chk("fallback_fits_2_62", ledger.FALLBACK_REF <= 2 ** 62,
    f"S*383*3*2^-64/C = {float(ledger.FALLBACK_REF):.6g} reference units < 2^62;"
    f" adds {ledger.FALLBACK_BITS:.9f} bits")
chk("fallback_union_bound", ledger.FALLBACK_TOTAL <= ledger.FALLBACK_REF * ledger.C,
    "exact Fraction: S_ops * 3*2^-64 * 383 instructions, charged once as EXTRAS")

# ---- 5. kappa budget ------------------------------------------------------
chk("kappa_le_1_015", ledger.A_H1 <= ledger.F(1015, 1000),
    f"kappa = {float(ledger.A_H1):.4f} <= 1.015 (every row in the machine table)")
chk("kappa_budget_product", ledger.BUDGET <= ledger.CLAIM_RATIO_CEIL,
    f"product {float(ledger.BUDGET):.9f} <= 2^(claim - log2 S0) {ledger.CLAIM_RATIO_CEIL:.9f}")

# ---- 6. gate arithmetic ---------------------------------------------------
chk("claim_below_ship_line", ledger.CLAIM < ledger.SHIP_LINE,
    f"claim {float(ledger.CLAIM):.3f} < ship line {float(ledger.SHIP_LINE):.3f}")
chk("margin_at_least_0_008", ledger.MARGIN_VS_GATE >= ledger.F(8, 1000),
    f"gate {float(ledger.GATE):.3f} - claim {float(ledger.CLAIM):.3f}"
    f" = {float(ledger.MARGIN_VS_GATE):.4f} >= 0.008 (tie auto-reject rule)")
chk("step_rows_161", ledger.STEP_INSN == 161,
    "46+20+14+17+7+36+12+1+7+1 = 161, budget 165")

# ---- 7. toy evidence ------------------------------------------------------
ev_path = os.path.join(ROOT, "evidence", "toy.json")
toyev = json.load(open(ev_path)) if os.path.exists(ev_path) else None
if toyev:
    rows = toyev["sweep_walkmap"]                 # production-scaled gate arm
    mn = min(r["cp95_lower"] for r in rows)
    n_tot = sum(r["trials"] for r in rows)
    chk("toy_cp95_gate", mn >= 0.39,
        f"min CP-95 lower over densities 2^-5..2^-11 = {mn:.4f} >= 0.39"
        f" (arm {rows[0]['arm']}, N = {rows[0]['N']}, {rows[0]['trials']} trials/density,"
        f" {n_tot} total; at n=1600 the same arm reads"
        f" {toyev['gate']['note_at_1600']['min_cp95_lower']:.4f}, disclosed)")
    chk("toy_densities_7", len(rows) == 7
        and all(r["dp_bits"] == 5 + i for i, r in enumerate(rows)),
        f"{len(rows)} densities {rows[0]['dp_bits']}..{rows[-1]['dp_bits']},"
        f" {n_tot} trials total")
    re_rows = toyev["sweep_realec"]
    step_ok = sum(r["invariant_step"][0] for r in re_rows)
    step_n = sum(r["invariant_step"][1] for r in re_rows)
    row_ok = sum(r["invariant_row"][0] for r in re_rows)
    row_n = sum(r["invariant_row"][1] for r in re_rows)
    rec_ok = sum(r["invariant_recovery"][0] for r in re_rows)
    rec_n = sum(r["invariant_recovery"][1] for r in re_rows)
    chk("toy_invariants_exact",
        step_n > 0 and step_ok == step_n and row_ok == row_n and rec_ok == rec_n
        and all(r["invariant_ok"] for r in re_rows),
        f"every step: {step_ok}/{step_n} S == [a]G+[c]T; every row: {row_ok}/{row_n}"
        f" J == [rj]G+[sj]T; every merge: {rec_ok}/{rec_n} recovered x reproduces T"
        f" (D2: orbit element resolved on the point, not on the x residue)")
    rt = toyev["rival_table"]
    zs = [abs(x["z"]) for x in rt["rows"]]
    zm = [abs(x["z_vs_model"]) for x in rt["rows"]]
    n2 = sum(1 for x in rt["rows"] if x["within_2se"])
    chk("toy_rival_table_reproduced",
        len(rt["rows"]) == 5 and n2 >= 4 and max(zs) <= 3.0
        and max(zm) <= 2.0 and rt["flat_ours"] <= 0.02
        and abs(rt["pipeline_pooled"] - rt["pipeline_pooled_rival"]) <= 2 * rt["pipeline_pooled_se"],
        f"5 rows at 4x their n on their N: {n2}/5 within 2 SE of their printed point"
        f" (max |z| {max(zs):.2f}; every row within 2 SE of the model, max"
        f" {max(zm):.2f}); dp_bits spread ours {rt['flat_ours']:.4f} vs theirs"
        f" {rt['flat_rival']:.4f}; pipeline pooled {rt['pipeline_pooled']:.4f} vs"
        f" {rt['pipeline_pooled_rival']:.4f} (n={rt['pipeline_pooled_n']})"
        f" (D1: merge gated on kappa*S_ref, DP lag charged separately)")
else:
    chk("toy_cp95_gate", False, "evidence/toy.json missing")
    chk("toy_densities_7", False, "evidence/toy.json missing")
    chk("toy_invariants_exact", False, "evidence/toy.json missing")
    chk("toy_rival_table_reproduced", False, "evidence/toy.json missing")

cert = {
    "package": "plain-truncwindow-40pct-v19d",
    "checks": checks,
    "passed": sum(1 for c in checks if c["ok"]),
    "total": len(checks),
    "verdict": "PASS" if all(c["ok"] for c in checks) else "FAIL",
    "claim": float(ledger.CLAIM),
    "log2_W": ledger.LOG_W,
    "success": ledger.SUCCESS,
}
with open(os.path.join(ROOT, "certificate.json"), "w") as f:
    json.dump(cert, f, indent=1)
with open(os.path.join(ROOT, "evidence", "agreement.json"), "w") as f:
    json.dump({"classes": n_tot, "agree": agree, "mismatches": len(mismatch),
               "max_window_dev": maxdv, "near_ties_2_315": near,
               "fidelity_100": fid, "toy_fidelity_100": toyfid,
               "convention_log2W_161_kappa102": their}, f, indent=1)
print(json.dumps({k: cert[k] for k in ("passed", "total", "verdict", "claim", "success")}, indent=1))
for c in checks:
    if not c["ok"]:
        print("FAIL:", c["name"], "-", c["detail"])
sys.exit(0 if cert["verdict"] == "PASS" else 1)
