#!/usr/bin/env python3
"""v31 retention-tail arithmetic: measured support integral + SLOPE-FITTED residual.

Routes, stated before the numbers:

  ROUTE (b) MEASURED   I_meas = mean over raw merges of e^{-(128 - A/theta)} on the
                       extended arm (cap_mult 128 == the production cap, budget 8x
                       the old support).  No bins, no envelope, no extrapolation.

  RESIDUAL             the only part of the production range with no observed age
                       is (a_hi, 128 theta].  There the charge is
                         P(A > a) <= p_hi * exp(-kappa*(a - a_hi))
                       with kappa FITTED by least squares on the raw survival
                       counts in [a_lo, a_hi] (window: survival >= MIN_COUNT),
                       p_hi a rule-of-three + 2x upper bound at a_hi, and the
                       integral taken in closed form up to the production cap.
                       Per-decade slopes and the fit residuals are shipped, so
                       kappa is visibly data, not a constant: truncating the
                       decades MOVES it (control (a) in audit_v31).

  PRODUCTION CHARGE    the fixture age law is toy-scope (785 classes); it is NOT
                       used as a production envelope.  The production retention
                       miss is derived from pair-uniformity on the 128-theta cap:
                       with the meeting point uniform over the trail pair, the
                       younger trail's age has density 2(L-u)/L^2 on [0,L],
                       L = 128 theta, and the distance to its report is
                       D ~ Exp(theta) with E[D] = theta, so Markov in D gives

                         delta <= E[ theta / (L - age) ] = 2/128 = 2^-6

                       exactly (computed below), with the exact-exponential row
                       2/128^2 = 2^-13 shipped as the tighter law.  delta = 2^-6
                       is what v31 charges; it is 155 bits LARGER than the
                       fixture support integral, i.e. strictly more conservative
                       than v30's measured support row, and it needs no
                       extrapolation beyond any support.
"""
import json
import math

PROD_CAP = 128.0          # production cap in theta units
MIN_COUNT = 20


def _surv(hist):
    """survival counts: {age_in_steps: #merges with age > age_in_steps}."""
    ages = sorted(hist)
    total = sum(hist.values())
    out = {}
    run = 0
    for a in reversed(ages):
        out[a] = run
        run += hist[a]
    out[0] = total
    return out, total


def fit(raw, truncate=None, min_count=MIN_COUNT):
    """Least-squares slope of log2 survival on the raw tail decades."""
    hist = {int(k): int(v) for k, v in raw['age_hist'].items()}
    if truncate is not None:
        hist = {a: v for a, v in hist.items() if a <= int(truncate)}
    theta = raw['theta']
    s, merges = _surv(hist)
    ages = sorted(a for a in s if a > 0)
    counted = [a for a in ages if s[a] >= min_count]
    if len(counted) < 8:
        raise SystemExit('fit window too short: %r' % counted)
    a_hi = max(counted)
    # TAIL decades only: the last 4 doublings of the observed support.  Fitting
    # the whole curve would let the steep early part set the rate.
    a_lo = max(min(counted), a_hi // 16)
    win = [a for a in counted if a >= a_lo]
    if len(win) < 8:
        raise SystemExit('tail window too short: %r' % win)
    xs = [float(a) for a in ages if a_lo <= a <= a_hi]
    ys = [math.log2(s[a] / float(merges)) for a in xs]
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    icpt = my - slope * mx
    resid = [y - (icpt + slope * x) for x, y in zip(xs, ys)]
    ss = sum(r * r for r in resid)
    se_slope = math.sqrt(ss / max(1, n - 2) / sxx)
    kappa_step = -slope
    kappa = kappa_step * theta
    kappa_se = se_slope * theta
    # per-decade (doubling) slopes: the slope is measured, not assumed
    decades = []
    for k in range(1, 6):
        hi = a_hi // (2 ** (k - 1))
        lo = max(a_lo, hi // 2)
        if hi <= lo:
            continue
        if s.get(hi, 0) <= 0 or s.get(lo, 0) <= 0:
            continue
        decades.append(dict(window=[lo, hi],
                            slope_per_step=(math.log2(s[hi] / float(merges))
                                            - math.log2(s[lo] / float(merges)))
                            / float(hi - lo)))
    k_hi = s[a_hi]
    p_hi = min(1.0, 2.0 * (k_hi + 2.0) / (merges + 2.0))
    u_hi = a_hi / float(theta)
    k_lo = max(1e-9, kappa - 2.0 * kappa_se)
    # closed-form residual integral on (u_hi, 128]
    if abs(k_lo - 1.0) < 1e-12:
        inner = math.exp(-u_hi) * (PROD_CAP - u_hi)
    else:
        inner = (math.exp((1.0 - k_lo) * PROD_CAP)
                 - math.exp((1.0 - k_lo) * u_hi)) / (1.0 - k_lo)
    resid_charge = p_hi * math.exp(-PROD_CAP) * math.exp(k_lo * u_hi) * inner
    # the fit's own worst deviation is charged as an inflation of the residual
    # (2^resid_max_bits), so a bad exponential fit cannot hide a tail
    infl = 2.0 ** max(abs(r) for r in resid) if resid else 1.0
    resid_charge = resid_charge * infl
    i_meas = 0.0
    for a, v in hist.items():
        i_meas += v * math.exp(-(PROD_CAP - a / float(theta)))
    i_meas /= merges
    return dict(theta=theta, merges=merges, truncate=truncate,
                a_lo_steps=a_lo, a_hi_steps=a_hi, u_lo=a_lo / float(theta),
                u_hi=u_hi, window_points=n, kappa_per_theta=kappa,
                kappa_se=kappa_se, kappa_lo=k_lo, slope_per_step=slope,
                intercept_log2=icpt, resid_max_bits=max(abs(r) for r in resid),
                resid_rms_bits=math.sqrt(ss / n), decade_slopes=decades,
                p_hi=p_hi, p_hi_counts=k_hi,
                integral_measured=i_meas,
                integral_measured_log2=(math.log2(i_meas) if i_meas > 0 else None),
                residual_charge=resid_charge,
                residual_log2=(math.log2(resid_charge) if resid_charge > 0 else None),
                delta_fixture=i_meas + resid_charge,
                delta_fixture_log2=math.log2(i_meas + resid_charge)
                if (i_meas + resid_charge) > 0 else None)


def production_delta(cap=PROD_CAP):
    """Pair-uniform production retention miss, computed, no fitted envelope.

    age density f(u) = 2 (L-u)/L^2 ; Markov in D (E[D] = theta) gives
    delta = int f(u) * theta/(L-u) du = 2/L = 2^-6.  Exact-exponential law:
    int f(u) e^{-(L-u)} du = 2 (1 - e^-L (1+L)) / L^2 = 2^-13.
    """
    L = float(cap)
    exact = 2.0 * (1.0 - math.exp(-L) * (1.0 + L)) / (L * L)
    markov = 2.0 / L
    return dict(cap_theta=L, delta_markov=markov,
                delta_markov_log2=math.log2(markov),
                delta_exact=exact, delta_exact_log2=math.log2(exact),
                charged=max(markov, exact), rule='pair-uniform age law; Markov in '
                'the report distance D (E[D]=theta); no tail envelope')


if __name__ == '__main__':
    raw = json.load(open('evidence/v31-tail-raw.json'))
    f = fit(raw)
    p = production_delta()
    print('window [%d,%d] steps = [%.1f,%.1f] theta  points %d'
          % (f['a_lo_steps'], f['a_hi_steps'], f['u_lo'], f['u_hi'],
             f['window_points']))
    print('kappa %.6f +/- %.6f per theta (lo %.6f)  resid max %.3f bits rms %.3f'
          % (f['kappa_per_theta'], f['kappa_se'], f['kappa_lo'],
             f['resid_max_bits'], f['resid_rms_bits']))
    print('decades %s' % [(d['window'], round(d['slope_per_step'], 4))
                          for d in f['decade_slopes']])
    print('I_meas 2^%.2f  residual 2^%s  fixture delta 2^%.2f'
          % (f['integral_measured_log2'], f['residual_log2'],
             f['delta_fixture_log2']))
    print('production delta markov 2^%.4f exact 2^%.4f charged 2^%.4f'
          % (p['delta_markov_log2'], p['delta_exact_log2'],
             math.log2(p['charged'])))
