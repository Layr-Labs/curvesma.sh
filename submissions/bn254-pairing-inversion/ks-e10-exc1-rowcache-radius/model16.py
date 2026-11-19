"""Evaluator source is model16core.py inside lib.b126.

Local edit (this submission): the b-tail row-cache debit [5c]/[5d] is charged
per radius bin instead of at the uniform worst case.  row_cache() of
model16core computes the strong-row ratio N6/Nstrong once, as the maximum over
256 certified radius bins, and charges the per-cell saving 138.4/600 at that
worst case for every cell.  The population is a sum over slices whose b-side
radius is R^2 = (lam - Q_a) A^2, and the same lattice_count_bounds.strong_rows
bound is certified on every one of those bins, so each slice can be charged at
its own ratio.  The weighted mean below uses the model's own slice_cells
weights (Geometry.slice_cells ~ (R^2/A2)^3), which are the weights that define
the cell population, and it is rounded down.  Nothing else changes: the
per-row credit (155.4 -> 64 -> 17), the short-row penalty, the build charges
and f are exactly as in model16core.
"""
import common12
from model16core import *
from model16core import _ctx
import model16core as _m16c
from fractions import Fraction as _F


def _strong_bins(geo, bins=256):
    """strong_rows' certified (N6/Nstrong, (N5-Nstrong)/N6) on each radius bin, with the model shell weights."""
    import lattice_count_bounds as _lc
    q = [[_F(v, geo.db) for v in r] for r in geo.qb]
    al = q[0][0]
    S = [[q[i][j] - q[i][0] * q[0][j] / al for j in range(1, 6)] for i in range(1, 6)]
    S4 = [r[:4] for r in S[:4]]
    mats = (q, S, S4, [r[:5] for r in q[:5]])
    det = [_lc.determinant_rational(m) for m in mats]
    assert det[0] == al * det[1]
    p = _lc.pi_interval(); e = 1 << 128
    p = _F(int(p[0] * e), e), _F(-(-p[1] * e // 1), e)
    s6, s5, s4, sz = map(_lc.sqrt_interval, det)
    t6, t5, t4, tz = (_lc.sqrt_interval(_lc.cell_radius(m)[1])[1] for m in mats)
    V5 = lambda r, k: 8 * p[k] ** 2 / 15 * r ** 5

    def ratio(a, b):
        rr = _lc.sqrt_interval(a - 64 * al); A = _lc.sqrt_interval(a); B = _lc.sqrt_interval(b)[1]
        assert rr[0] > t5 and A[0] > t6
        low = (V5(rr[0] - t5, 0) / s5[1] - p[1] ** 2 / 2 * (rr[1] + t4) ** 4 / s4[0]) / 2
        n6 = (p[0] ** 3 / 6 * (A[0] - t6) ** 6 / s6[1] - V5(A[1] + tz, 1) / sz[0]) / 2
        assert low > 0 and n6 > 0
        return (float(p[1] ** 3 / 12 * (B + t6) ** 6 / s6[0] / low),
                float((V5(B + t5, 1) / s5[0] / 2 - low) / n6))

    A2 = _F(_ctx_a2); lam = _F(_ctx_lam); th = _F(_ctx_th)
    lo = (lam - th) * A2; hi = lam * A2
    out = []
    for j in range(bins):
        a = lo + (hi - lo) * j / bins; b = lo + (hi - lo) * (j + 1) / bins
        u, v = ratio(a, b)
        out.append((float((b / A2) ** 3 - (a / A2) ** 3), u, v))   # Geometry.slice_cells / vb / upper
    return out


_ctx_a2 = _ctx_lam = _ctx_th = 0


def row_cache(cfg, geo):
    global _ctx_a2, _ctx_lam, _ctx_th
    _ctx_a2 = geo.A2
    _ctx_lam = cfg['production_region_lam']
    _ctx_th = cfg['production_theta_region']
    rc = _m16c.row_cache(cfg, geo)
    bins = _strong_bins(geo)
    den = sum(w for w, _, _ in bins)
    inv = sum(w / u for w, u, _ in bins) / den      # cells-weighted mean of 1/(N6/Nstrong)
    xs = sum(w * v for w, _, v in bins) / den       # cells-weighted mean short-row share
    d = (C_NB_ROW - C_ROW_LOOKUP) * inv - (C_ROW_LOOKUP - C_NB_ROW / 16) * xs + (C_ROW_LOOKUP - C_DESC) * inv
    rc['inv_strong'] = inv
    rc['xshort_w'] = xs
    rc['d'] = int(d * 1e6) / 1e6                     # rounded down: under-credit the saving
    return rc
