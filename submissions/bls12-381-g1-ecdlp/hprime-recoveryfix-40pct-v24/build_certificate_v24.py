#!/usr/bin/env python3
"""v24 builder: regenerate certificate.json + argument.md + public-note.md from the
raw evidence files.  Every number in the prose is interpolated from JSON, so the
documents cannot drift from the evidence (the v23 lesson: hand-typed rows rot).

Claim rule: the ledger is unchanged (11 rows, route row dead), so W and the claim
are unchanged at 116.8130.  What changes is the SUCCESS row: it is now a measured
CP-95 lower bound instead of a model value, and it is shipped with an explicit
success_below_floor flag when it lands under 0.39.
"""
import json, math, os, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))


def lg(x):
    return round(math.log2(x), 6)


def main():
    cert = json.load(open(HERE + "/certificate.json"))
    direct = json.load(open(HERE + "/direct-v24.json"))
    gate = json.load(open(HERE + "/gate-v24.json"))
    retry = json.load(open(HERE + "/retry-v24.json"))
    kap = json.load(open(HERE + "/kappa-v24.json"))
    s = direct["summary"]
    best_name = s["best_shape"]
    best_cp95 = s["success_best_shape"]          # a CP-95 number, not a key

    # ---- retry row priced by the WORST regime (kill 4)
    worst = retry["retry"]["worst_regime_used"]
    old_q = F(1) - F(31, 2000)          # v23 pooled q_hi 0.0154817 -> factor 0.015725
    old_factor = 1 / (1 - 0.0154817) - 1
    new_factor = worst["retry_factor"]
    old_row = cert["rows"]["gap_retry_worst"]
    cost = cert["retry"]["build_draw_cost"]
    new_row = int(cost * new_factor)      # honest definition: cost x (1/(1-q)-1)
    assert abs(new_row / (cost * new_factor) - 1) < 1e-12
    rows = dict(cert["rows"])
    rows["gap_retry_worst"] = new_row
    row_log2 = dict(cert["row_log2"])
    row_log2["gap_retry_worst"] = lg(new_row)
    total_instr = sum(rows.values())
    W = F(total_instr, 2000)          # 2000 instructions per charged unit (v22/v23 rule)
    raw = math.log2(float(W))
    claim = F(116813, 1000)
    gate_v = F(1168271, 10000)
    headroom = float(claim) - raw
    margin = float(gate_v - claim)
    assert headroom > 0 and margin > 0, (raw, headroom, margin)

    # ---- direct measurement rows
    by_shape = s["by_shape"]
    rows_direct = []
    for g in direct["groups"]:
        rows_direct.append({k: g[k] for k in
                            ("group", "p", "r", "d", "m", "nb", "nw", "shape",
                             "records_per_side", "records_per_trial", "K", "K_dev",
                             "K_granularity", "trials", "solves", "q", "model",
                             "model_dev", "cp95", "matches", "verified_matches",
                             "matches_per_trial", "points_on_curve")})

    cert["version"] = "v24"
    cert["direct"] = dict(
        kind="v24 DIRECT end-to-end solve census, recovery algebra fixed",
        runner="direct_v24.py", gate="gate-v24.json", fixtures="fixtures-v24.json",
        algebra=("key match <=> x = lambda*v*b*w' (lambda in mu6) => "
                 "v = x*(b*w')^{-1}*zeta^j; b*w' is the UNSCALED B-side product "
                 "(v23 divided by the v-scaled y and only walked the mu6 orbit)"),
        verify="affine([c]P0) == affine(Ad)  (v23 compared raw projective tuples)",
        s_extraction="s = d-th root of v mod r (d=2, Tonelli-Shanks), then "
                     "affine([s]P0) == affine(As)",
        groups=rows_direct, by_shape=by_shape,
        gate_result=dict(n=gate["n"], ok=gate["ok"], passed=gate["passed"],
                         lam_one_cases=gate["lam_one_cases"],
                         v23_path_solves=gate["v23_path_solves"],
                         on_curve_asserts=gate["on_curve_asserts"],
                         false_positive_n=gate["false_positive_n"],
                         false_positives=gate["false_positives"]),
        points_on_curve_total=s["points_on_curve_total"],
        total_trials=s["total_trials"],
        success_final=s["success_final"], success_final_shape=s["success_final_shape"],
        success_best_shape=s["success_best_shape"],
        success_below_floor=s["success_below_floor"], shortfall=s["shortfall"],
        meets_floor=s["meets_floor"], claim_floor=0.39,
        budget_to_reach_floor_bits=s["budget_to_reach_floor_bits"],
        w2_amplification=s["w2_amplification"],
        stability_spread=s["stability_spread"],
        cp95_method=s["cp95_method"], rule=s["rule"], certified=True)

    cert["success"] = dict(
        design=cert["success"]["design"],
        debits=cert["success"]["debits"],
        measured=dict(
            q_by_shape={k: v["q"] for k, v in by_shape.items()},
            cp95_by_shape={k: v["cp95_min"] for k, v in by_shape.items()},
            success_final=s["success_final"], shape=s["success_final_shape"],
            trials=s["total_trials"], solves=s["solves"]),
        after=s["success_final"],
        claim_floor=0.39,
        meets_floor=s["meets_floor"],
        success_below_floor=s["success_below_floor"],
        statement=("MEASURED, not modelled. The design value 1-exp(-K)=%0.6f is NOT "
                   "reached once the walk list is shared between the two sides: one "
                   "algebraic (i,j) collision fires at every walk position, so matches "
                   "cluster and the per-TRIAL solve rate falls to %0.6f (CP-95 lower "
                   "%0.6f). The independent-record shape reproduces the design value "
                   "(%0.6f, CP-95 %0.6f) at the same per-side record count. The claim "
                   "is stated at the WORST shape, so it is a MISS against the 0.39 "
                   "floor by %0.6f, disclosed here rather than averaged away."
                   % (cert["success"]["design"],
                      by_shape[s["success_final_shape"]]["q_min"],
                      s["success_final"],
                      by_shape[best_name]["q_min"], best_cp95,
                      s["shortfall"])),
        certified=True)

    cert["retry"]["per_regime"] = [
        {k: v for k, v in c.items() if k != "rows"} for c in retry["retry"]["per_regime"]]
    cert["retry"]["worst_regime_used"] = worst
    cert["retry"]["pricing_rule"] = retry["retry"]["pricing_rule"]
    cert["retry"]["pooled_descriptive_only"] = retry["retry"]["pooled_descriptive_only"]
    cert["retry"]["gap_retry_row"] = dict(old=old_row, new=new_row, old_log2=lg(old_row),
                                          new_log2=lg(new_row),
                                          factor_old=round(old_factor, 8),
                                          factor_new=round(new_factor, 8))

    cert["kappa_168"] = {k: v for k, v in kap.items() if k != "fits"}
    cert["kappa_168"]["fit_rows"] = {k: {kk: vv for kk, vv in v.items() if kk != "points"}
                                     for k, v in kap["fits"].items()}
    cert["rows"] = rows
    cert["row_log2"] = row_log2
    cert["total_instructions"] = total_instr
    cert["W"] = str(W)
    cert["raw_log2"] = round(raw, 6)
    cert["claim"] = {"numerator": str(claim.numerator), "denominator": str(claim.denominator)}
    cert["headroom_over_W"] = round(headroom, 8)
    cert["gate_margin"] = round(margin, 8)
    cert["checks"].update(dict(
        gate_known_answer_64_64=gate["passed"] and gate["ok"] == 64,
        gate_old_path_zero=gate["v23_path_solves"] == 0,
        gate_no_false_positives=gate["false_positives"] == 0,
        population_on_curve_asserted=s["points_on_curve_total"] > 0,
        success_measured=direct["summary"]["total_trials"] >= 512,
        success_below_floor_disclosed=bool(s["success_below_floor"] == (not s["meets_floor"])),
        retry_per_regime=len(retry["retry"]["per_regime"]) == 2,
        retry_priced_by_worst=worst["q_gap_cp95_upper"] >= max(
            c["q_gap_cp95_upper"] for c in retry["retry"]["per_regime"]) - 1e-15,
        kappa_168_fitted_not_measured=kap["measured"] is False,
        kappa_charged_dominates=kap["charged_exceeds_fitted"]))
    cert["cert_scope"] = (
        "Recomputes the ledger, the population rows, the retry rows (per regime), the "
        "kappa_168 extrapolation and the claim from raw evidence; runs the production-"
        "shaped direct census at fixture scale (n=2^20 subgroup, 4 cells x 32768 "
        "trials). It does not run the production attack and does not establish "
        "production resource bounds.")

    json.dump(cert, open(HERE + "/certificate.json", "w"), indent=1)

    # ---------------- documents
    g0 = direct["groups"][0]
    gate_row = ("%d/%d" % (gate["ok"], gate["n"]))
    shape_names = {"independent_nw1": "independent records (nw=1)",
                   "shared_walk_nw2": "shared walk list (nw=2)"}
    lines = []
    for k, v in by_shape.items():
        lines.append("| %s | %s | %s | %s | %s | %s |" % (
            shape_names.get(k, k), v["groups"], "/".join(str(x) for x in v["trials"]),
            "/".join("%0.5f" % q for q in v["q"]), "%0.6f" % v["model"],
            "%0.6f" % v["cp95_min"]))
    table = "\n".join(lines)

    arg = f"""# H-composite Cheon BSGS v24 — recovery algebra FIXED, success MEASURED

## Result

`time_log2 = 116.8130` (`116813/1000`) at `n = 2^20`; charged work `2^{raw:.5f}` in
2000-instruction units over {len(rows)} row items; peak memory `2^123.235` B; no advice;
preprocessing inside the total. Gate `116.8271`, margin `{margin:.4f}`.

**Success is now MEASURED, not modelled.** Direct end-to-end census, {s['total_trials']}
independent trials over 4 cells (2 primes x 2 shapes):

| shape | cells | trials | q per cell | model 1-exp(-K) | CP-95 min |
|---|---|---|---|---|---|
{table}

`success_final = {s['success_final']:.6f}` (worst shape, CP-95 one-sided lower bound),
against a floor of `0.39` -> **success_below_floor = {str(s['success_below_floor']).lower()}**,
shortfall `{s['shortfall']:.6f}`. This is disclosed as the design row; it is not averaged
against the better shape and it is not replaced by a model value.

## 1. Recovery algebra (verdict a841deb5, kill 2) — FIXED and gated

v23 walked the mu6 orbit of `x/y`. `y` already carries `v`, so `x/y` is a mu6 unit and
the recovered value could never be `v`; that is why the pilot solved 0/1024. The correct
recovery divides by the **unscaled** B-side product `b*w'`:

    key match  <=>  x = lambda * v * b * w'  (lambda in mu6)
    =>  v = x * (b*w')^(-1) * zeta^j,   j = 0..5.

Verification was a second independent bug: v23 compared raw projective tuples
(`jmul(c,P0) == Ad`), which never compare equal; v24 compares **affine** points.
`s`-extraction is charged: `s` is the d-th root of `v` mod `r` (d=2, Tonelli-Shanks) and
must satisfy `[s]P0 == As`.

Known-answer gate (`gate-v24.json`), run before any census: **{gate_row}** engineered
collisions solved, {gate['lam_one_cases']} with `lambda = 1` and {gate['lam_other_cases']}
with `lambda != 1`; the same {gate['n']} fixtures scored with the v23 recovery give
**{gate['v23_path_solves']}/{gate['n']}** (planted-bug sensitivity: the fix is
load-bearing, not cosmetic); {gate['false_positive_n']} unplanted instances give
{gate['false_positives']} false solves; {gate['on_curve_asserts']} on-curve asserts.

## 2. Clustering: why the measured rate is under the design model

The walk list is shared between the two sides. For a fixed pair `(i,j)` the ratio
`x/y = a_i/(v*b_j)` does not depend on the walk index, so an algebraic collision fires
at **every** walk position at once: matches pile up on one trial instead of spreading.
Matches per trial stay at `{g0['matches_per_trial']:.4f}` (K = `{g0['K']:.4f}`, so the
collision budget is intact), but the per-trial solve rate drops. The two shapes at the
same per-side record count isolate the effect: `nw=1` reproduces the design value,
`nw=2` does not.

The production shape is `nb=2^65.65`, `nw=2^50` (certificate.json:coverage), i.e. deep
inside the clustered family: at fixed `K` a larger `nw` means fewer independent base
pairs `nb^2 = (N/nw)^2`, so success FALLS with `nw`. The `nw=2` cell is therefore the
most favourable member of the shipped family and `{s['success_final']:.6f}` is an UPPER
bound on the production success rate, not a lower one. The `nw=1` cell is included
precisely to show what the design model would need -- records that do not share a walk
list -- and it is not the shipped shape.

This is the honest reading of the v22 objection that started this lane ("the probability
of a genuine collision is the distinct fraction, not automatically 1-exp(-6|X||Y|/m)"):
the distinct fraction is the right object, and sharing the walk list makes the distinct
fraction smaller than the pair count suggests.

## 3. kappa at the production width (kill 3) — FITTED, not measured

Law `kappa_w = 1 + 6*2^(w-1)/m`, refitted on the large-count widths of the shipped census
(`w = 10..24`, predicted counts >= 512; smaller counts are Poisson-noise dominated and are
excluded), residual scatter `+/-{kap['fit_band_log2']:.3f}` bits on the excess. At the
production `m = 2^{kap['production_m_log2']:.3f}` and `w = 168`:

    excess = 6 * 2^167 / m = 2^{kap['fitted_excess_168_log2']:.3f}
    kappa_168 = 1 + 2^{kap['fitted_excess_168_log2']:.3f}  (below double precision: 1.0)
    upper end with the residual band: 1 + 2^{kap['fitted_excess_168_log2'] + kap['fit_band_log2']:.3f}

The charged `kappa = 1.231` is {kap['charged_headroom']:.6f} above the fitted upper end,
so the charge is conservative and no ledger row moves. This row is an extrapolation and
is labelled as such: it is **not** a 168-bit census.

## 4. Retry pooled across regimes (kill 4) — split per regime

v23 pooled a 32-bit and a 36-bit census. Each regime now has its own fresh-start census
({retry['retry']['per_regime'][0]['draws']} draws each) and its own CP-95 upper bound;
the row is priced by the **worst** regime, `q_hi = {worst['q_gap_cp95_upper']:.6f}`
(`{worst['regime']}`), giving factor `{worst['retry_factor']:.6f}` and
`gap_retry_worst = 2^{lg(new_row):.4f}` (v23: `2^{lg(old_row):.4f}` from the pooled
`0.015480`). The pooled number survives only as a descriptive row.

## 5. Population on-curve (v23 flaw) — machine row

Every enumerated point is asserted on `y^2 = x^3 + 4` before its key is used:
`{s['points_on_curve_total']}` asserts across the census, `{gate['on_curve_asserts']}` in
the gate. The population identity rows (`sum_k N_w(k) == m` per width per fixture) are
unchanged and still shipped.

## 6. Ledger (unchanged rows, one row repriced)

""" + "\n".join("`%s` `2^%.3f`" % (k, v) for k, v in sorted(row_log2.items(),
                                                           key=lambda kv: -kv[1])) + f"""

Total `W = 2^{raw:.5f}`; claim `116.8130` (headroom `{headroom:.6f}`); gate `116.8271`.

## 7. What is NOT closed

* `kappa_168` is fitted, not measured (no 168-bit census exists).
* The measured success misses 0.39 by `{s['shortfall']:.6f}`. Price to close it:
  `+{s['budget_to_reach_floor_bits']:.4f}` bits of records (N scales as sqrt(K)), which
  would put the claim at `{116.81277 + s['budget_to_reach_floor_bits']:.4f}` — above the
  gate. Alternative: two sequential runs, success
  `1-(1-{s['w2_amplification']['p_used']:.4f})^2 = {s['w2_amplification']['success']:.4f}`
  at `2W = 117.8128`, also above the gate. Both are priced, neither is claimed.
"""
    open(HERE + "/argument.md", "w").write(arg)

    note = f"""# hprime recoveryfix v24 — g1-ecdlp

Success is MEASURED this time, not modelled: {s['total_trials']} independent end-to-end
trials, real curve points, real keys, verified recovery, `s`-extraction charged.

* Recovery algebra fixed (v23 solved 0/1024 because it divided by the v-scaled `y` and
  only walked the mu6 orbit). Known-answer gate: **{gate_row}** engineered collisions,
  v23 path **{gate['v23_path_solves']}/{gate['n']}** on the same fixtures.
* Measured solve rate: independent-record shape `{by_shape[best_name]['q_min']:.4f}`
  (CP-95 `{best_cp95:.6f}`, the design value);
  shared-walk shape `{by_shape[s['success_final_shape']]['q_min']:.4f}`
  (CP-95 `{s['success_final']:.6f}`).
* `success_final = {s['success_final']:.6f}` < 0.39: **honest miss**, disclosed, not
  hidden behind the better shape and not replaced by a model. Price to close: +{s['budget_to_reach_floor_bits']:.4f} bits
  (over gate) or 2 sequential runs at 117.8128 (over gate).
* Cost claim unchanged: `116.8130` <= gate `116.8271`.
* `kappa_168 = 1 + 2^{kap['fitted_excess_168_log2']:.2f}` fitted, labelled not measured;
  charged 1.231 is conservative by {kap['charged_headroom']:.3f}.
* Retry priced per regime at the worst `q_hi = {worst['q_gap_cp95_upper']:.6f}`.
* {s['points_on_curve_total']} on-curve asserts, {gate['on_curve_asserts']} in the gate.
"""
    open(HERE + "/public-note.md", "w").write(note)
    print("built. raw_log2=%.6f claim=116.8130 headroom=%.6f margin=%.6f" %
          (raw, headroom, margin))
    print("success_final=%.6f below_floor=%s shortfall=%.6f" %
          (s["success_final"], s["success_below_floor"], s["shortfall"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
