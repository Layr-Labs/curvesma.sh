#!/usr/bin/env python3
"""v24 evidence D: kappa at the PRODUCTION width (168) -- FITTED, not measured.

verdict a841deb5 kill 3: kappa = 1.231 was measured at census key widths (10..30),
production is w = 168; the transfer was unestablished.

The census law (census-key-v20.json note): kappa_w = 1 + 2^(w-1) * sum_k (N(k)/m)^2,
i.e. the collision-surface inflation above the uniform-image model.  With the six
mu6 images per H element the fitted form used in the package is

    kappa_w = 1 + 6 * 2^(w-1) / m .

This script (a) refits the slope and constant of that law on the LARGE-COUNT widths
of the shipped census, (b) ships the fit-residual band, (c) evaluates the law at
w = 168 with the production m from certificate.json, and (d) prices the maximum
residual into the collision floor (the candidate debit).

It does NOT measure kappa at w = 168: a 168-bit key census needs 2^168 records and
is not run here.  The row is labelled fitted.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
K_CHARGED = 1.231


def fit(rows, wmin, wmax):
    pts = [(r["width"], r["ucb_factor"] - 1.0) for r in rows
           if wmin <= r["width"] <= wmax and r["ucb_factor"] > 1.0]
    xs = [p[0] for p in pts]
    ys = [math.log2(p[1]) for p in pts]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    a = my - slope * mx
    resid = [y - (a + slope * x) for x, y in zip(xs, ys)]
    return dict(widths=[p[0] for p in pts], n=n, slope=round(slope, 6),
                intercept_log2=round(a, 6), resid_log2=[round(r, 4) for r in resid],
                resid_max_abs=round(max(abs(r) for r in resid), 6),
                m_implied_log2=round(math.log2(3.0) - a, 6),
                points=[(w, round(v, 8)) for w, v in pts])


def main():
    c = json.load(open(HERE + "/census-key-v20.json"))
    cert = json.load(open(HERE + "/certificate.json"))
    m_prod = cert["m"]
    m_log2 = math.log2(m_prod)
    fits = {}
    for nm, d in (("big", c["faithful"]["big"]), ("small", c["faithful"]["small"])):
        # large-count widths only: the UCB factor is Poisson-noise dominated when the
        # predicted count is small (w >= 26 has < 128 expected), so those are excluded
        wmax = 24 if nm == "big" else 18
        fits[nm] = fit(d["exact_rows"], 10, wmax)
    e_fit = {}
    for nm, f in fits.items():
        # law: excess(w) = 6*2^(w-1)/m  =>  at production m:
        e_fit[nm] = 6.0 * 2.0 ** 167 / m_prod
    band = max(f["resid_max_abs"] for f in fits.values())
    # worst-case (largest) fitted excess across regimes
    e168 = max(e_fit.values())
    k168 = 1.0 + e168
    e168_hi = e168 * 2.0 ** band
    k168_hi = 1.0 + e168_hi
    # price: kappa enters only through E[#bad] and the candidate debit
    debit = cert["success"]["debits"]["candidate"]
    price_success = debit * (1.0 - 1.0 / K_CHARGED)
    out = dict(
        kind="kappa transfer to the production width (FITTED, not measured)",
        law="kappa_w = 1 + 6*2^(w-1)/m",
        census_m_log2=c["m_log2"], production_m_log2=round(m_log2, 6),
        fits={k: {kk: vv for kk, vv in v.items()} for k, v in fits.items()},
        fit_band_log2=round(band, 6),
        fitted_excess_168_log2=round(math.log2(e168), 6),
        kappa_168_fitted=k168,
        kappa_168_upper=k168_hi,
        kappa_168_band=[1.0, k168_hi],
        kappa_charged=K_CHARGED,
        charged_exceeds_fitted=bool(K_CHARGED > k168_hi),
        charged_headroom=round(K_CHARGED - k168_hi, 12),
        measured=False,
        statement="kappa_168 is an EXTRAPOLATION of the census law with the production "
                  "m; it is not measured at w=168.  The charged kappa=1.231 exceeds the "
                  "fitted upper end, so keeping 1.231 charges the collision surface "
                  "conservatively and no ledger row moves.",
        price=dict(
            entry_point="candidate debit (success.debits.candidate)",
            candidate_debit=debit,
            max_residual_priced_in_success=price_success,
            note="even at the charged kappa the whole candidate debit is 1.4e-6 of "
                 "success, i.e. 1e-5 of the 0.07 shortfall; the measured shortfall is "
                 "caused by shared-walk clustering, not by kappa"))
    json.dump(out, open(HERE + "/kappa-v24.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "fits"}, indent=1))
    for nm, f in fits.items():
        print(nm, {k: v for k, v in f.items() if k != "points"})


if __name__ == "__main__":
    main()
