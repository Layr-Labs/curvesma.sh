"""Heuristic H1+ (proof.md section 3.4): the drift the ledger's 1.25 factor absorbs, the drift table, and the drift
at which the 0.39 floor or the 0.01 allowance is exhausted. Loaded by experiment.py part D."""

Dec = __import__("decimal").Decimal


def h1plus(Pm, S, eps, Npd, sumq2, loss_sum, p2, eps2, eps_long, margin, D, lg, floor):
    model_waste = D(eps2 + eps_long)
    allow = Dec("0.01")

    def bisect(f, lo, hi):
        for _ in range(80):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if f(mid) else (lo, mid)
        return lo

    def success_at_drift(points):
        s_eff = S * (1 + eps) / (1 + model_waste + Dec(points) / 100)
        return 1 - (-(s_eff * s_eff / (2 * Npd))).exp() - loss_sum

    tol_points = (eps - model_waste) * 100
    assert abs(success_at_drift(tol_points) - (Pm - loss_sum)) < Dec("1e-30")
    drift_lo = bisect(lambda k: success_at_drift(k) >= Dec("0.39"), Dec(0), Dec(100))
    assert drift_lo > tol_points
    table = {k: float(floor(success_at_drift(k))) for k in ("0", "0.25", "1", "1.5", "2", "3", "5")}
    sd_2cycle = (D(p2) * (1 - D(p2)) / (S * (1 + eps))).sqrt()
    sd_log2 = lg(D(eps2) * (D(margin) - 1) / sd_2cycle)
    assert sd_log2 > 59
    assert 66 < -lg(sd_2cycle) < 68 and 0.01 < D(eps2) * (D(margin) - 1) * 100 < 1

    def succ_idx(delta):
        n_eff = Npd * (1 - D(sumq2)) / (1 - D(sumq2) - delta)
        return 1 - (-((S * (1 + eps)) ** 2) / (2 * n_eff)).exp() - loss_sum

    succ0 = succ_idx(Dec(0))
    idx_lo = bisect(lambda d: succ0 - succ_idx(d) <= allow, Dec(0), Dec(1))
    assert idx_lo > Dec("0.02") and succ0 - succ_idx(Dec(1)) > allow
    return {"model_waste": model_waste, "allow": allow, "tol_points": tol_points, "drift_lo": drift_lo,
            "table": table, "sd_log2": sd_log2, "sd_2cycle_log2": -lg(sd_2cycle),
            "excess_2cycle_points": D(eps2) * (D(margin) - 1) * 100, "idx_lo": idx_lo}
