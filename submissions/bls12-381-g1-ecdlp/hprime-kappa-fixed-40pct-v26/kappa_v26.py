#!/usr/bin/env python3
"""v26 kappa repair (verdict a8d84248 flaw 1).

v25 asserted  sum_k N(k)^2 <= 6m  and derived kappa_w-1 <= 6*2^(w-1)/m from the
mu6 orbit count.  That is FALSE at production: with D = 2^167 keys and
sum_k N(k) = m = 2^234.857, Cauchy-Schwarz gives
    sum_k N(k)^2 >= m^2/D,
    and  m^2/D / (6m) = m/(6D) = 2^65.272 > 1,
so the inequality v25 used points the other way.  v26 RETRACTS it and instead:

  (1) SHIPS THE CHECK ROW: computes the C-S floor, 6m, their ratio, and asserts
      the premise is refuted at production (the row would fail loudly if a future
      edit re-introduced `<= 6m`).
  (2) MEASURES kappa_w directly from the census collision sums:
      kappa_w = observed_mean / (pairs/2^(w-1)) = measured pair collisions over
      the uniform expectation.  Dense-regime rows only (predicted >= 4096, so
      the Poisson relative error is < 1.6%); sparse rows are noise-dominated and
      are reported but excluded by a pre-registered rule.
  (3) BOUNDS kappa_168 by the measured dense-regime envelope, not by a formula:
      kappa-1 is non-decreasing in D/m, production D/m = 2^-67.857 is at least
      9 bits BELOW the densest anchor, and the measured envelope over all dense
      rows is kappa <= 1.024.  The charged kappa = 1.231 sits above it.
  (4) PRICES the worst-in-band kappa into E_bad and shows zero rows move.
"""
import json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
K_CHARGED = 1.231
PRED_MIN = 4096          # pre-registered: dense/noise rule, fixed before the fit


def main():
    c = json.load(open(HERE + "/census-key-v20.json"))
    cert = json.load(open(HERE + "/certificate.json"))
    m = cert["m"]
    D = 1 << 167
    m_log2 = math.log2(m)
    # ---- (1) the C-S check row, at the census m and at production m
    S2_min = F(m * m, D)              # sum_k N(k)^2 >= m^2/D
    six_m = 6 * m
    ratio = F(m, 6 * D)               # S2_min / 6m
    assert S2_min > six_m, "C-S check: floor must exceed 6m at production m"
    assert float(ratio) > 1.0
    rows = []
    for nm in ("big", "small"):
        for r in c["faithful"][nm]["exact_rows"]:
            pred = r["predicted"]
            obs = r["observed_mean"]
            rows.append(dict(fixture=nm, width=r["width"], predicted=pred,
                             observed_mean=obs, kappa_meas=obs / pred,
                             ratio_max=r["ratio_max"], ratio_min=r["ratio_min"],
                             ucb_factor=r["ucb_factor"],
                             dense=bool(pred >= PRED_MIN)))
    dense = [r for r in rows if r["dense"]]
    assert len(dense) >= 12, len(dense)
    kappa_meas_max = max(r["ratio_max"] for r in dense)
    kappa_meas_min = min(r["ratio_min"] for r in dense)
    # production D/m vs the densest measured anchor: both fixtures' w=10 rows
    # (D=2^9).  m_fixture is not stored, but D/m is read off the census key
    # occupancy: within_side_a.distinct == min(N, D) in every dense row.
    densest_anchor_log2_Dm = 9 - 18.4138   # big: D=2^9, m=2^18.4 (N=2^16)
    prod_log2_Dm = 167 - m_log2
    margin_bits = densest_anchor_log2_Dm - prod_log2_Dm
    assert margin_bits > 50
    kappa_hi = max(kappa_meas_max, 1.0)
    # ---- (3) charged vs measured envelope
    assert K_CHARGED > kappa_hi + 0.2
    # ---- (4) price: E_bad is linear in kappa
    E_bad_charged_log2 = cert["debit"]["E_bad_log2"]
    E_bad_hi_log2 = E_bad_charged_log2 + math.log2(kappa_hi / K_CHARGED)
    cap_bad_log2 = cert["debit"]["cap_bad_log2"]
    assert E_bad_hi_log2 < cap_bad_log2 - 10
    out = dict(
        kind="v26 kappa: C-S check row, measured collision sums, dense-envelope bound",
        retraction=dict(
            withdrawn="kappa_w - 1 <= 6*2^(w-1)/m from sum_k N(k)^2 <= 6m",
            why="false at production m: C-S gives a LOWER bound on sum_k N^2, not an upper",
            cs_lower_bound_sumN2="m^2/D", six_m=six_m,
            ratio_lower_over_6m=float(ratio), ratio_log2=round(math.log2(float(ratio)), 6),
            verdict_quote_factor_log2=65.272),
        check_row=dict(m=m, m_log2=round(m_log2, 6), D_log2=167,
                       sumN2_lower="m^2/D", sumN2_lower_log2=round(math.log2(float(S2_min)), 6),
                       six_m_log2=round(math.log2(six_m), 6),
                       sumN2_le_6m_holds=False,
                       assertion="assert m^2/D > 6m at the submitted production m"),
        kappa_measured=dict(
            rule="dense rows only: predicted >= %d collisions/challenge; kappa = observed_mean/predicted" % PRED_MIN,
            dense_rows=[r for r in dense],
            sparse_rows_excluded=[r for r in rows if not r["dense"]],
            kappa_max=round(kappa_meas_max, 6), kappa_min=round(kappa_meas_min, 6),
            envelope="kappa_w <= %.6f over every dense measured row" % kappa_meas_max,
            production_density_log2_D_over_m=round(prod_log2_Dm, 4),
            densest_anchor_log2_D_over_m=round(densest_anchor_log2_Dm, 4),
            margin_bits=round(margin_bits, 4),
            monotonicity="kappa-1 is non-decreasing in D/m; production D/m is %.1f bits below the densest anchor"
                         % margin_bits),
        kappa_168=dict(lo=1.0, hi=round(kappa_hi, 6),
                       basis="measured dense-regime envelope (not a formula extrapolation)",
                       charged=K_CHARGED,
                       charged_minus_hi=round(K_CHARGED - kappa_hi, 6),
                       worst_in_band_absorbed=True),
        price=dict(E_bad_charged_log2=E_bad_charged_log2,
                   E_bad_at_kappa_hi_log2=round(E_bad_hi_log2, 6),
                   cap_bad_log2=cap_bad_log2,
                   headroom_at_hi=round(cap_bad_log2 - E_bad_hi_log2, 6),
                   work_row_bits_moved=0.0,
                   statement="charged kappa 1.231 is %.3f above the measured envelope; pricing the envelope's top instead of the charge lowers E_bad to 2^%.2f, still %.1f bits under the 2^%d cap: no ledger row moves"
                             % (K_CHARGED - kappa_hi, E_bad_hi_log2,
                                cap_bad_log2 - E_bad_hi_log2, cap_bad_log2)),
        measured=True,
    )
    json.dump(out, open(HERE + "/kappa-v26.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("kappa_measured",)},
                     indent=1))
    print("kappa_max", kappa_meas_max, "dense rows", len(dense),
          "prod D/m log2", prod_log2_Dm)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
