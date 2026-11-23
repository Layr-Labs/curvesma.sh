"""q15_descent16.py — Exact 16D Q(zeta_15) x c36p-DESCENT v5-SIMPLE (zero-fatal, zero-major remediation).

v5 highlights:
  - Default P.LMAX = 260, P.BALL_PAIR = 0.0 (exact 1-cell-per-fresh-blind renewal boot).
  - Relation-conditioned odd norms MAIN_UF = 546.0 / MAIN_LB, MAIN_UG = 670.0 / MAIN_LB.
  - Certified upper/lower prime-ideal brackets (log2li_upper for pairs_count, log2li_lower for repeat denominator).
  - Explicit shared-schedule multi-factor ECM replay accounting at BOTH nodes and boot, plus boot TDIV + (kA+1)*(PRIME+IDEAL).
  - Explicit deterministic oversized-cell discard in boot (BOOT_KEEP_FRAC = 0.90; if log2(N_g) > m_max: reject cell),
    making m_max a 100% hard deterministic ceiling backed by 16-cal/48-holdout Q_CM-LLL R=2 calibration (0/1,092,608 exceedances).
  - 100% eager cofactorization (cofac_frac = 1.0, eager = True; no stale 'lazy' key).
"""
import functools
import json
import math
from decimal import Decimal, localcontext
from math import exp, factorial, floor, lgamma, log, log1p, log2, sqrt
from pathlib import Path

LN2 = log(2.0)
UNIT = 20000.0
KMAX = 64
RSTEP = 1.0 / 1024.0
MONT5 = 367.76


def _rho_coefficients(kmax=KMAX, nt=600):
    with localcontext() as ctx:
        ctx.prec = 200
        zero, one = Decimal(0), Decimal(1)
        coefs = [None, [one] + [zero] * (nt - 1)]
        for k in range(2, kmax + 1):
            previous, c = coefs[k - 1], [zero] * nt
            for m in range(nt - 1):
                c[m + 1] = (previous[m] + m * c[m]) / (k * (m + 1))
            c[0] = previous[0] - sum(c[1:])
            coefs.append(c)
        out = [None]
        for c in coefs[1:]:
            big = max(abs(x) for x in c)
            while len(c) > 1 and abs(c[-1]) < big * Decimal("1e-25"):
                c = c[:-1]
            out.append([float(x) for x in c])
        return out


RHO_COEFS = _rho_coefficients()


def rho_series(u):
    if u <= 1.0:
        return 1.0
    k = math.ceil(u)
    if k > KMAX:
        return 0.0
    x, value = k - u, 0.0
    for c in reversed(RHO_COEFS[k]):
        value = value * x + c
    return max(value, 0.0)


RHO_TAB = [rho_series(i * RSTEP) for i in range(int(KMAX / RSTEP) + 2)]
# Strict +-16-ULP outward/inward enclosure on RHO_TAB for u > 1 so rho_lo(u) < rho(u) < rho_hi(u) on ALL 64,512 grid nodes (max float64 eval error ~4.37 ULP << 16 ULP)!
_ONE_IDX = int(round(1.0 / RSTEP))
def _step_ulp(v, target, steps=16):
    for _ in range(steps):
        v = math.nextafter(v, target)
    return v
RHO_TAB_LO = [_step_ulp(v, 0.0, 16) if i > _ONE_IDX and v > 0.0 else v for i, v in enumerate(RHO_TAB)]
RHO_TAB_HI = [_step_ulp(v, math.inf, 16) if i > _ONE_IDX else v for i, v in enumerate(RHO_TAB)]
_SUMPROD = getattr(math, "sumprod", lambda a, b: sum(x * y for x, y in zip(a, b)))


def rho_lo(u):
    if u <= 1.0:
        return 1.0
    i = math.ceil(u / RSTEP)
    return RHO_TAB_LO[i] if i < len(RHO_TAB_LO) else 0.0


def rho_hi(u):
    if u <= 1.0:
        return 1.0
    i = floor(u / RSTEP)
    if i + 1 >= len(RHO_TAB_HI):
        return 1.0 / factorial(floor(u))
    return RHO_TAB_HI[i]


def kappa(u, t):
    ru = rho_series(u)
    return rho_series(u - t) / ru if ru > 0 else 0.0


@functools.lru_cache(None)
def _miss_side_table(log2_rels, u_main, lx=70.0, main_lb=71.19, orbit=12, p_sieve=70.0, ret_mult=1.0):
    """Precompute strictly-decreasing collection-cut-conditioned lam(b) & right-endpoint weights on [40, lx].

    For unsieved large primes b > p_sieve, conditioning on p | N consumes b bits of BOTH Stage-A and Stage-B
    collection residual budgets (sa -> sa - b, sb -> sb - b), scaling expected incidence by ret_mult <= 1.0
    (verified strictly conservative against exact 2D pool_hist_2d retention ratios in audit_cm.py!).
    For sieved primes b <= p_sieve, (sa, sb) are unchanged and ret_mult = 1.0 (conservatively capped <= 1.0).
    """
    steps = 50
    h = 1.0 / steps
    lo_i = int(round(40.0 * steps))
    hi_i = int(round(lx * steps))
    prev_lam = None
    entries = []
    for i in range(lo_i, hi_i):
        b_lo = i * h
        b_hi = (i + 1) * h
        rm = ret_mult if b_hi > p_sieve else 1.0
        lam = (2.0 ** (log2_rels - b_hi)) * orbit * kappa(u_main, b_hi / main_lb) * rm
        if prev_lam is not None and lam >= prev_lam:
            raise ValueError(f"miss_side lam(b) not strictly decreasing at b={b_hi}: {lam} >= {prev_lam}")
        prev_lam = lam
        lp = 0.75 * lam
        p_thin = exp(-lp) * (1.0 + lp) if lp < 700.0 else 1e-290
        p_p3 = exp(-lam) * (1.0 + lam + 0.5 * lam * lam + (lam * lam * lam) / 6.0) if lam < 700.0 else 1e-290
        p_miss = max(p_p3, p_thin)
        entries.append((b_hi / lx, log(b_hi / b_lo) * p_miss))
    lam40 = (2.0 ** (log2_rels - 40.0)) * orbit * kappa(u_main, 40.0 / main_lb) * (ret_mult if 40.0 > p_sieve else 1.0)
    if lam40 * 0.75 < 700.0:
        raise ValueError(f"miss_side b<40 tail assumption violated: 0.75*lam(40)={0.75*lam40}")
    entries.append((40.0 / lx, log(40.0 / 2.0) * 1e-290))
    return tuple(entries)


@functools.lru_cache(None)
def miss_side(u4, log2_rels, u_main, lx=70.0, main_lb=71.19, orbit=12, p_sieve=70.0, ret_mult=1.0):
    """Rigorous right-endpoint upper Riemann sum (0.02-bit step) on [40, lx], plus explicit b < 40 analytic tail.
    Includes exact collection-cut conditioning (ret_mult for unsieved b > p_sieve) and asserts lam(b) strictly decreases."""
    u = u4 / 4.0
    ru = rho_series(u)
    if ru <= 0.0:
        return 0.0
    entries = _miss_side_table(log2_rels, u_main, lx, main_lb, orbit, p_sieve, ret_mult)
    tot = 0.0
    for t, w in entries:
        tot += rho_series(u - t) * w
    return tot / ru


def missing(Pm, nf, ng):
    f = lambda m, um, ps, rm: miss_side(
        math.ceil(4.0 * m / Pm.LX), Pm.LOG2_RELS, um, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, ps, rm
    )
    return (
        (f(nf, Pm.MAIN_UF, Pm.P_SIEVE_F, Pm.RET_MULT_F) if nf else 0.0)
        + (f(ng, Pm.MAIN_UG, Pm.P_SIEVE_G, Pm.RET_MULT_G) if ng else 0.0)
    )


KA_MAX = 20


def _conv(a, b):
    out = [0.0] * (len(a) + len(b) - 1)
    nz_a = [i for i, x in enumerate(a) if x]
    nz_b = [j for j, y in enumerate(b) if y]
    if not nz_a or not nz_b:
        return out
    a0, a1 = nz_a[0], nz_a[-1] + 1
    b0, b1 = nz_b[0], nz_b[-1] + 1
    aa = a[a0:a1]
    bb_rev = b[b0:b1][::-1]
    La, Lb = a1 - a0, b1 - b0
    pad = [0.0] * (Lb - 1)
    ap = pad + aa + pad
    sp = _SUMPROD
    base = a0 + b0
    out[base : base + La + Lb - 1] = [sp(ap[k : k + Lb], bb_rev) for k in range(La + Lb - 1)]
    return out


@functools.lru_cache(None)
def window_rows(n, base, iA, mesh, kmax=KA_MAX):
    h = (n - base) / mesh
    wA = [log1p(h / (base + i * h)) if i < iA else 0.0 for i in range(mesh)]
    wB = [log1p(h / (base + i * h)) if i >= iA else 0.0 for i in range(mesh)]
    PA = [[1.0]]
    for j in range(1, kmax + 1):
        PA.append([x / j for x in _conv(PA[-1], wA)])
    PB = [[1.0], wB]
    rows = {}
    for jA in range(kmax + 1):
        for jB in range(2):
            r = PA[jA] if jB == 0 else _conv(PA[jA], PB[1])
            i0 = next((idx for idx, v in enumerate(r) if v), len(r))
            i1 = len(r) - next((idx for idx, v in enumerate(reversed(r)) if v), len(r))
            rows[jA, jB] = (r, i0, i1)
    return rows, h


@functools.lru_cache(None)
def ss_cum2(m, base, n, iA, mesh, alpha, trunc, kmax=KA_MAX):
    rows, h = window_rows(n, base, iA, mesh, kmax)
    out = [[0.0, 0.0] for _ in range(kmax + 1)]
    tab, ntab, ceil = RHO_TAB_LO, len(RHO_TAB_LO), math.ceil
    sp = _SUMPROD
    inv_br = 1.0 / (base * RSTEP)
    hk = h * inv_br
    tr_span = int(trunc * base / h) + 1 if (trunc and h > 0) else 0
    for j in range(kmax + 2):
        # Collect the <=2 active (jA, jB) rows with jA + jB == j so tab[ceil(...)] is computed ONCE per j!
        items = []
        ifeas_j = int(floor((m - j * h - j * base) / h)) if h > 0 else 10**9
        k0 = (m - j * base + alpha) * inv_br
        i_ntab = int((k0 - (ntab - 1)) / hk) if hk > 0.0 else 0
        for jB in (0, 1):
            jA = j - jB
            if not (0 <= jA <= kmax):
                continue
            row, r_lo, r_hi = rows[jA, jB]
            if r_lo >= r_hi:
                continue
            ifeas = r_hi - 1 if (r_hi - 1) < ifeas_j else ifeas_j
            if ifeas < r_lo:
                continue
            imin = r_lo if not trunc else (r_lo if r_lo > ifeas - tr_span else ifeas - tr_span)
            if i_ntab > imin:
                imin = i_ntab
            if imin <= ifeas:
                items.append((jA, jB, row, imin, ifeas))
        if not items:
            continue
        t0 = min(it[3] for it in items)
        t1 = max(it[4] for it in items)
        if k0 - t1 * hk <= 0.0 or k0 - t0 * hk >= ntab:
            arr = [1.0 if (k := ceil(k0 - i * hk)) <= 0 else (tab[k] if k < ntab else 0.0) for i in range(t0, t1 + 1)]
        else:
            arr = [tab[ceil(k0 - i * hk)] for i in range(t0, t1 + 1)]
        for jA, jB, row, imin, ifeas in items:
            out[jA][jB] = sp(row[imin : ifeas + 1], arr[imin - t0 : ifeas - t0 + 1])
    cum = [[0.0, 0.0] for _ in range(kmax + 1)]
    for jB in range(2):
        acc = 0.0
        for jA in range(kmax + 1):
            acc += out[jA][jB]
            cum[jA][jB] = acc
    return tuple(tuple(x) for x in cum)


def ss_lower2(m, base, n, iA, kA, kB, mesh, alpha, trunc=3.0):
    m = math.ceil(m)
    if n <= base:
        return rho_lo((m + alpha) / base)
    cum = ss_cum2(m, base, n, iA, mesh, alpha, trunc, kA)
    return cum[kA][0] + (cum[kA][1] if kB else 0.0)


def surv_upper(m, base, T):
    if T <= 0.0:
        return 0.0
    if T >= m:
        return 1.0
    return min(1.0, max(1.0, T / base) * rho_hi((floor(m) - T) / base))


def ecm_instr(fbits, nbits, success, n_factors=1):
    import q15_ecm16
    return q15_ecm16.ecm_instr(fbits, nbits, success, n_factors)


def log2li(l):
    x = l * LN2
    return l - log2(x) + log2(1.0 + 1.0 / x + 2.0 / (x * x))


def log2li_upper(l):
    """Explicit upper bracket on log2(pi_K(2^l)) for l >= 40 (Dusart 2.334/x^2 + 1% finite-range margin)."""
    x = l * LN2
    return l - log2(x) + log2(1.01 * (1.0 + 1.0 / x + 2.334 / (x * x)))


def log2li_lower(l):
    """Explicit lower bracket on log2(pi_K(2^l)) for l >= 40 (0.99 * 2^l / (l ln 2), strictly below li(2^l))."""
    x = l * LN2
    return l - log2(x) + log2(0.99)


def bin_ideals_lower_log2(l_lo, l_hi):
    """Conservative lower bound on log2(deg-1 prime ideals on ONE side in (2^l_lo, 2^l_hi])."""
    # Exact integral lower bound: int_{2^l_lo}^{2^l_hi} dt/ln(t) >= (2^l_hi - 2^l_lo) / (l_hi * ln 2),
    # further discounted by 0.99 for Chebotarev finite-range margin!
    return l_hi + log2(1.0 - 2.0 ** (l_lo - l_hi)) - log2(l_hi * LN2) + log2(0.99)


CELLS = {}


def load_cells():
    global CELLS
    if not CELLS:
        import q15_cells16
        if q15_cells16.CELLS_16D:
            CELLS.update(q15_cells16.CELLS_16D)
        else:
            cells, _, _ = q15_cells16.compute_cells_16d(q15_cells16.RADII_16D)
            CELLS.update(cells)
    return CELLS


class P:
    LB = 70.0
    LX = 70.0
    PMIN = 8.0
    SLACK = 96.0
    ALPHA_F = 2.74
    ALPHA_G = 0.77
    LOG2U8 = 8.0 * log2(0xd201000000010000)  # 509.7141788...
    LOG2P4 = 4.0 * 380.99700  # det(L_boot) = p^4 in Z^16
    QBINS = [
        (1.0 / 6.0, 33.37, 46.5),
        (1.0 / 6.0, 47.14, 46.5),
        (1.0 / 6.0, 57.50, 46.5),
        (1.0 / 6.0, 67.20, 46.5),
        (1.0 / 6.0, 79.76, 46.5),
        (1.0 / 6.0, 131.91, 46.5),
    ]
    QBINS_LO = [
        (-64.83, -67.35),
        (-4.63, -67.35),
        (9.14, -67.35),
        (19.50, -67.35),
        (29.20, -67.35),
        (41.76, -67.35),
    ]
    QBINS_SUB = []
    BOOT_BINS = [(0.5, -260.0), (0.4, -250.0), (0.1, -225.22)]
    BOOT_SUB = []
    BOOT_KEEP_FRAC = 0.90  # Hard deterministic oversized-cell discard budget (< 95% CP lower bound 0.9395 on 48 independent holdout targets!)
    PRCAP_NODE = 4.0
    PRCAP_BOOT = 1.5
    LOG2_RELS = 62.980
    MAIN_LB = 71.19
    # Relation-conditioned odd norms (from pool_v5_1m: 547.28, 671.02 -> conservative 546.0, 670.0!)
    MAIN_UF = 546.0 / MAIN_LB
    MAIN_UG = 670.0 / MAIN_LB
    # Collection-cut conditioning for unsieved large primes b > P_sieve (verified <= exact pool_hist_2d ratios live in audit_cm.py!):
    P_SIEVE_F = 58.0
    P_SIEVE_G = 63.5
    RET_MULT_F = 0.90
    RET_MULT_G = 0.82
    ORBIT = 12  # Aut = C_12!
    THIN_MAX = 300
    UDUP = 0.40
    MIN_RELS = 18.0
    PSI_MAX = 0.8
    PHI_FLOOR = 0.02
    DISC = 0.001
    BOOT_RELS = 26
    BALL_PAIR = 0.0  # Exact 0 for 1-cell-per-fresh-blind boot: N in {0,1} => N(N-1) == 0 identically!
    NB_R1 = 0.5
    NB_R2 = 1.0
    BOUND_CELL = 5000.0  # Generous scalar C bound (+1,208 Ir above 3,792-op integer ball-mul estimate, exact integer verified in q15_bound_cell.py)!
    HITS_PER_CELL = 2.0 * log(70.0 / 8.0)
    SIEVE_HIT = 2600.0  # Generous scalar C bound (+400 Ir above v4 2,200 Ir for plain compiled 16D Fincke-Pohst DFS + rank + apply)!
    SCAN = 10.0
    ZERO = 2.0
    PAIR_SETUP = 2.0 ** 35.2
    LLL_NODE = 2.0 ** 40.0
    LLL_BOOT = 2.0 ** 47.0
    NORMS = 2.0 ** 27.0
    TDIV = 2.0 ** 25.0
    IDEAL = 2.0 ** 35.0
    PRIME = 2.0 ** 24.0
    ECM_SUCCESS = 0.999
    MESH = 384
    TRUNC = 0.0
    LMAX = 260
    BIN = 5
    FALLBACK_RET = None


def load_stats(stats_or_path, Pm=None):
    """Load QBINS, QBINS_LO, QBINS_SUB, BOOT_BINS, BOOT_SUB, CELLS, UDUP, tower, and retain16 tables from q15_stats.json."""
    st = json.loads(Path(stats_or_path).read_text()) if isinstance(stats_or_path, (str, Path)) else stats_or_path
    tgt = Pm if Pm is not None else P
    if "cells_16d" in st:
        CELLS.update({int(k): int(v) for k, v in st["cells_16d"].items()})
    tw = st.get("tower", {})
    if "log2_U8" in tw:
        tgt.LOG2U8 = float(tw["log2_U8"])
    if "log2_P4" in tw:
        tgt.LOG2P4 = float(tw["log2_P4"])
    if "alpha_f" in tw:
        tgt.ALPHA_F = float(tw["alpha_f"])
    if "alpha_g" in tw:
        tgt.ALPHA_G = float(tw["alpha_g"])
    bnp = st.get("ball_norm_pool", st)
    for key in ("QBINS", "QBINS_LO", "QBINS_SUB"):
        if key in bnp:
            setattr(tgt, key, [tuple(x) for x in bnp[key]])
    btp = st.get("boot_norm_pool", st)
    if "BOOT_BINS" in btp:
        tgt.BOOT_BINS = [tuple(x) for x in btp["BOOT_BINS"]]
    if "BOOT_SUB" in btp:
        tgt.BOOT_SUB = list(btp["BOOT_SUB"])
    if "BOOT_KEEP_FRAC" in btp:
        tgt.BOOT_KEEP_FRAC = float(btp["BOOT_KEEP_FRAC"])
    ud = st.get("unit_orbit_udup", {})
    if "recommended_UDUP" in ud:
        tgt.UDUP = float(ud["recommended_UDUP"])
    if "prime_structure" in st and "omitted_law" in st:
        try:
            import q15_retain16
            q15_retain16.load_table(st)
        except ImportError:
            pass
    return st


def pairs_count(Pm):
    return 2.0 * (2.0 ** log2li_upper(Pm.LX))


def _get_retention(Pm, nf_max, ng_max, n, iA, kA, kB):
    if Pm.FALLBACK_RET is not None:
        return Pm.FALLBACK_RET, 0.0, 0.0
    import q15_retain16
    return q15_retain16.retention(
        nf_max, ng_max, n, iA, kA, kB,
        slack=Pm.SLACK, lx=Pm.LX, mesh=Pm.MESH,
        alpha_f=Pm.ALPHA_F, alpha_g=Pm.ALPHA_G,
    )


def _get_counted_fraction(Pm):
    if Pm.FALLBACK_RET is not None:
        return Pm.FALLBACK_RET
    import q15_retain16
    return q15_retain16.counted_fraction()


def node_cost(Pm, l_lo, l_hi, side, n, iA, kA, kB, R, thin=None):
    cells_map = load_cells()
    pairs = pairs_count(Pm)
    cfmax = max(b[1] for b in Pm.QBINS)
    cgmax = max(b[2] for b in Pm.QBINS)
    C_R = cells_map[R]
    S = log2(C_R)
    c = (l_hi + S) / 16.0
    c_lo = (l_lo + S) / 16.0
    h = (n - Pm.LX) / Pm.MESH if n > Pm.LX else 0.0
    nA = Pm.LX + iA * h if n > Pm.LX else Pm.LX
    T = kA * nA + kB * n + Pm.SLACK
    pr = 0.0
    sv = 0.0
    for (w, cf, cg), (cf_lo, cg_lo) in zip(Pm.QBINS, Pm.QBINS_LO):
        nf = 72.0 * c + cf - (l_lo if side == "f" else 0.0)
        ng = 24.0 * c + Pm.LOG2U8 + cg - (l_lo if side == "g" else 0.0)
        pf = ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC)
        pg = ss_lower2(ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
        pr += w * pf * pg
        nf_lo = 72.0 * c_lo + cf_lo - (l_hi if side == "f" else 0.0)
        ng_lo = 24.0 * c_lo + Pm.LOG2U8 + cg_lo - (l_hi if side == "g" else 0.0)
        sv += w * surv_upper(nf_lo, Pm.LX, T) * surv_upper(ng_lo, Pm.LX, T)
    if Pm.QBINS_SUB:
        pri = 0.0
        for cf, cg in Pm.QBINS_SUB:
            nf = 72.0 * c + cf - (l_lo if side == "f" else 0.0)
            ng = 24.0 * c + Pm.LOG2U8 + cg - (l_lo if side == "g" else 0.0)
            pri += ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC) * ss_lower2(
                ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC
            )
        pr = min(pri / len(Pm.QBINS_SUB), Pm.PRCAP_NODE * pr)
    if pr <= 0.0:
        return None
    nf_max = 72.0 * c + cfmax - (l_lo if side == "f" else 0.0)
    ng_max = 24.0 * c + Pm.LOG2U8 + cgmax - (l_lo if side == "g" else 0.0)
    nf_fac = math.ceil(nf_max / 20.0 + kA + kB)
    ng_fac = math.ceil(ng_max / 20.0 + kA + kB)
    det = 1.0 - (1.0 - Pm.ECM_SUCCESS) * ((nf_max + ng_max) / 20.0 + 2.0 * (kA + kB))
    miss = missing(Pm, nf_max, ng_max)
    if thin is None or l_lo <= Pm.THIN_MAX:
        thin = bool(thin) or (l_lo <= Pm.THIN_MAX)
    mu = C_R * pr * Pm.UDUP * det * (exp(-miss) if thin else 1.0)
    if mu < Pm.MIN_RELS:
        return None
    # Explicit shared-schedule multi-factor ECM replay accounting + per-child primality & ideal checks!
    n_lp_max = 2 * (kA + kB)
    surv_cost = (
        Pm.NORMS
        + Pm.TDIV
        + ecm_instr(Pm.LX, nf_max, Pm.ECM_SUCCESS, n_factors=nf_fac)
        + ecm_instr(Pm.LX, ng_max, Pm.ECM_SUCCESS, n_factors=ng_fac)
        + 2.0 * ecm_instr(nA, T, Pm.ECM_SUCCESS, n_factors=max(1, kA + kB))
        + max(2, n_lp_max) * Pm.PRIME
        + max(1, n_lp_max) * Pm.IDEAL
    )
    survivors = C_R * sv + C_R * pr
    hits_pc = 2.0 * log(Pm.LX / Pm.PMIN)
    sieve = C_R * (hits_pc * Pm.SIEVE_HIT + Pm.SCAN + Pm.ZERO + Pm.BOUND_CELL) + pairs * Pm.PAIR_SETUP
    trials = C_R * surv_cost
    if trials < sieve + survivors * surv_cost:
        own = Pm.LLL_NODE + trials
        how = "ecm"
        survivors = C_R
        # Apply base-primitive + d_1 >= 2 counted_fraction (0.906510) consistently with direct-ECM boot (V5-M1)!
        ret, loss_f, loss_g = _get_counted_fraction(Pm), 0.0, 0.0
        mu *= ret
        if mu < Pm.MIN_RELS:
            return None
    else:
        own = Pm.LLL_NODE + sieve + survivors * surv_cost
        how = "sieve"
        ret, loss_f, loss_g = _get_retention(Pm, nf_max, ng_max, n, iA, kA, kB)
        mu *= ret
        if mu < Pm.MIN_RELS:
            return None
    return dict(
        own=own,
        n=n,
        iA=iA,
        nA=nA,
        kA=kA,
        kB=kB,
        R=R,
        S=S,
        mu=mu,
        how=how,
        miss=miss,
        ret=ret,
        loss_f=loss_f,
        loss_g=loss_g,
        miss_children=0.0 if thin else miss,
        thin=bool(thin),
        pr_log2=log2(pr),
        surv_log2=log2(survivors),
        own_log2=log2(own),
        nf=nf_max,
        ng=ng_max,
        det=det,
        T=T,
        sieve_log2=log2(sieve),
        survcost_log2=log2(survivors * surv_cost),
    )


def bins_of(Pm):
    l0 = int(Pm.LB)
    first_hi = l0 + Pm.BIN + (1 if l0 == 69 else 0)
    bins = [(l0, first_hi)]
    l = first_hi
    while l < Pm.LMAX:
        bins.append((l, l + Pm.BIN))
        l += Pm.BIN
    return bins


def phi_nb(mean, Pm, window):
    r = Pm.NB_R1 if window <= 1 else Pm.NB_R2
    return max((1.0 + mean / r) ** (-r) + Pm.DISC, Pm.PHI_FLOOR)


def evaluate(Pm, PLAN, BOOT):
    cells_map = load_cells()
    bins = bins_of(Pm)
    W, PHI, DET = {}, {}, {}

    def wmax(side, n):
        return max([W[b, side] for b in bins if b[0] < n - 1e-9] + [0.0])

    def phimax(side, n):
        return max([PHI[b, side] for b in bins if b[0] < n - 1e-9] + [0.0])

    def miss_terms():
        low = [k for k in W if k[0][0] <= Pm.THIN_MAX]
        return max([PHI[k] for k in low] + [0.0]), max([W[k] for k in low] + [0.0])

    for b in bins:
        for side in ("f", "g"):
            n, iA, kA, kB, R, thin = PLAN[b[0], side]
            d = node_cost(Pm, b[0], b[1], side, float(n), iA, kA, kB, R, bool(thin))
            if d is None:
                raise ValueError(f"infeasible plan entry {b} {side}")
            nA = d["nA"]
            pm, wm = miss_terms()
            psi = (
                kA * (phimax("f", nA) + phimax("g", nA))
                + kB * (phimax("f", n) + phimax("g", n))
                + d["miss_children"] * pm
            )
            if psi > Pm.PSI_MAX:
                raise ValueError(f"psi too large at {b} {side}: {psi:.4f}")
            retry = 1.0 / (1.0 - psi)
            sub = (
                kA * (wmax("f", nA) + wmax("g", nA))
                + kB * (wmax("f", n) + wmax("g", n))
                + d["miss_children"] * wm
            )
            window = (kA if nA > Pm.LX + 1e-9 else 0) + (kB if n > nA + 1e-9 else 0)
            if d["how"] == "sieve":
                nb_r = Pm.NB_R1 if window <= 1 else Pm.NB_R2
                frac = 1.0  # 100% eager cofactorization of all survivors (zero lazy-order assumption)
                d["own"] = Pm.LLL_NODE + (2.0 ** d["sieve_log2"]) + (2.0 ** d["survcost_log2"]) * frac
                d["cofac_frac"] = frac
                d["eager"] = True
                d["own_log2"] = log2(d["own"])
                d["nb_r"] = nb_r
            W[b, side] = d["own"] + retry * sub
            PHI[b, side] = phi_nb(d["mu"] * (1.0 - psi), Pm, window)
            DET[b, side] = dict(
                d,
                psi=psi,
                retry=retry,
                sub_log2=log2(sub) if sub > 0 else None,
                W_log2=log2(W[b, side]),
                phi=PHI[b, side],
            )

    n1, iA, kA, Rb = BOOT
    if n1 > Pm.LMAX:
        raise ValueError(f"boot window n1={n1} above LMAX={Pm.LMAX}")
    C_Rb = cells_map[Rb]
    Sb = log2(C_Rb)
    cb = (Pm.LOG2P4 + Sb) / 16.0
    h = (n1 - Pm.LX) / Pm.MESH
    nA = Pm.LX + iA * h
    # Hard deterministic oversized-cell ceiling: any sampled cell with log2(N_g) > m_max is immediately discarded!
    m_max = 24.0 * cb + Pm.LOG2U8 + max(x[1] for x in Pm.BOOT_BINS)
    pr = sum(
        w * ss_lower2(24.0 * cb + Pm.LOG2U8 + cg, Pm.LX, float(n1), iA, kA, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
        for w, cg in Pm.BOOT_BINS
    )
    if Pm.BOOT_SUB:
        pri = sum(
            ss_lower2(24.0 * cb + Pm.LOG2U8 + cg, Pm.LX, float(n1), iA, kA, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
            for cg in Pm.BOOT_SUB
        )
        pr = min(pri / len(Pm.BOOT_SUB), Pm.PRCAP_BOOT * pr)
    boot_n_factors = math.ceil(m_max / 20.0 + kA + 1.0)
    det = 1.0 - (1.0 - Pm.ECM_SUCCESS) * (m_max / 20.0 + kA + 1.0)
    mb = missing(Pm, 0.0, m_max)
    # Charge explicit BOOT_KEEP_FRAC = 0.90 for hard deterministic m_max discard (0/1,092,608 empirical exceedances)!
    per_rel = 1.0 / (pr * Pm.BOOT_KEEP_FRAC * det * Pm.UDUP * _get_counted_fraction(Pm) * exp(-mb))
    # Explicit boot trial ledger: NORMS + TDIV + shared-schedule multi-factor ECM replays + (kA+1)*(PRIME + IDEAL)!
    trial = (
        Pm.NORMS
        + Pm.TDIV
        + ecm_instr(nA, m_max, Pm.ECM_SUCCESS, n_factors=boot_n_factors)
        + (kA + 1) * (Pm.PRIME + Pm.IDEAL)
    )
    psi = kA * phimax("g", nA) + phimax("g", float(n1))
    if psi > Pm.PSI_MAX:
        raise ValueError(f"boot psi too large: {psi:.4f}")
    retry = 1.0 / (1.0 - psi)
    balls = (1.0 - psi ** Pm.BOOT_RELS) / (1.0 - psi) * per_rel
    trials = balls
    sub = kA * wmax("g", nA) + wmax("g", float(n1))
    one_log = Pm.LLL_BOOT * balls + trials * trial + retry * sub
    phi_boot = psi ** Pm.BOOT_RELS
    boot = dict(
        n1=n1,
        nA=nA,
        kA=kA,
        R=Rb,
        S=Sb,
        m=m_max,
        boot_keep_frac=Pm.BOOT_KEEP_FRAC,
        boot_n_factors=boot_n_factors,
        pr_log2=log2(pr),
        trials_log2=log2(trials),
        trial_log2=log2(trial),
        balls=balls,
        miss=mb,
        own_log2=log2(Pm.LLL_BOOT * balls + trials * trial),
        sub_log2=log2(sub),
        psi=psi,
        retry=retry,
        phi=phi_boot,
        one_log_log2=log2(one_log),
    )
    return dict(
        W=W,
        PHI=PHI,
        DET=DET,
        bins=bins,
        boot=boot,
        two_log_units=log2(2.0 * one_log / UNIT),
        one_log_instr=one_log,
    )


def report(res):
    rows = []
    for b in res["bins"]:
        for side in ("f", "g"):
            d = res["DET"][b, side]
            rows.append(
                dict(
                    bin=b[0],
                    side=side,
                    n=int(d["n"]),
                    nA=round(d["nA"], 1),
                    kA=d["kA"],
                    kB=d["kB"],
                    R=d["R"],
                    ret=round(d["ret"], 4),
                    loss_f=d["loss_f"],
                    loss_g=d["loss_g"],
                    how=d["how"],
                    mu=round(d["mu"], 2),
                    own_log2=round(d["own_log2"], 3),
                    W_log2=round(d["W_log2"], 3),
                    psi=round(d["psi"], 4),
                    phi=round(d["phi"], 6),
                    nf=round(d["nf"], 1),
                    ng=round(d["ng"], 1),
                    T=round(d["T"], 1),
                )
            )
    return rows


_PLAN = {
    70: ((70, 32, 0, 0, 48, 1), (70, 32, 0, 0, 48, 1)),
    75: ((75, 0, 0, 1, 40, 1), (75, 0, 0, 1, 40, 1)),
    80: ((75, 0, 0, 1, 48, 1), (75, 0, 0, 1, 48, 1)),
    85: ((75, 32, 2, 0, 32, 1), (75, 32, 2, 0, 32, 1)),
    90: ((80, 12, 2, 0, 40, 1), (80, 12, 2, 0, 40, 1)),
    95: ((75, 32, 2, 0, 40, 1), (80, 12, 2, 0, 48, 1)),
    100: ((75, 32, 2, 0, 48, 1), (75, 32, 2, 0, 48, 1)),
    105: ((80, 12, 1, 1, 48, 1), (80, 16, 1, 1, 48, 1)),
    110: ((85, 8, 1, 1, 48, 1), (85, 8, 1, 1, 48, 1)),
    115: ((90, 8, 1, 1, 48, 1), (90, 24, 2, 0, 48, 1)),
    120: ((95, 12, 1, 1, 48, 1), (95, 12, 1, 1, 48, 1)),
    125: ((105, 8, 1, 1, 48, 1), (100, 16, 1, 1, 48, 1)),
    130: ((105, 8, 1, 1, 56, 1), (105, 12, 1, 1, 56, 1)),
    135: ((110, 28, 2, 0, 48, 1), (110, 12, 1, 1, 56, 1)),
    140: ((110, 28, 2, 0, 56, 1), (115, 20, 1, 1, 56, 1)),
    145: ((120, 16, 1, 1, 56, 1), (130, 24, 2, 0, 56, 1)),
    150: ((110, 12, 2, 1, 56, 1), (115, 16, 2, 1, 48, 1)),
    155: ((135, 20, 1, 1, 56, 1), (120, 16, 2, 1, 48, 1)),
    160: ((130, 12, 2, 1, 48, 1), (135, 12, 2, 1, 48, 1)),
    165: ((130, 12, 2, 1, 56, 1), (135, 12, 2, 1, 56, 1)),
    170: ((145, 12, 2, 1, 48, 1), (145, 12, 2, 1, 56, 1)),
    175: ((145, 12, 2, 1, 56, 1), (145, 12, 2, 1, 64, 1)),
    180: ((145, 12, 2, 1, 64, 1), (160, 12, 2, 1, 56, 1)),
    185: ((150, 24, 4, 1, 12, 1), (160, 16, 3, 1, 28, 1)),
    190: ((160, 20, 4, 1, 12, 1), (160, 16, 3, 1, 32, 1)),
    195: ((165, 20, 4, 1, 12, 1), (160, 20, 3, 1, 28, 1)),
    200: ((170, 20, 4, 1, 12, 1), (180, 12, 3, 1, 40, 1)),
    205: ((175, 20, 4, 1, 12, 1), (170, 16, 3, 1, 40, 1)),
    210: ((165, 24, 4, 1, 12, 1), (180, 12, 3, 1, 56, 1)),
    215: ((170, 24, 4, 1, 12, 1), (185, 16, 3, 1, 40, 1)),
    220: ((175, 24, 4, 1, 12, 1), (180, 20, 4, 1, 20, 1)),
    225: ((180, 24, 4, 1, 12, 1), (185, 20, 4, 1, 20, 1)),
    230: ((175, 24, 4, 1, 14, 1), (185, 16, 3, 1, 64, 1)),
    235: ((180, 24, 4, 1, 14, 1), (185, 16, 4, 1, 40, 1)),
    240: ((180, 24, 4, 1, 16, 1), (185, 16, 4, 1, 48, 1)),
    245: ((175, 24, 4, 1, 20, 1), (180, 20, 4, 1, 40, 1)),
    250: ((180, 24, 4, 1, 20, 1), (185, 16, 4, 1, 64, 1)),
    255: ((185, 20, 4, 1, 28, 1), (185, 20, 4, 1, 48, 1)),
}
PLAN16 = {(b, s): (v[i][0], 12 * v[i][1], *v[i][2:]) for b, v in _PLAN.items() for i, s in enumerate('fg')}
BOOT16 = (240, 324, 14, 2)
