"""q15_descent16.py - 16D Q(zeta_15) descent model: one certified trial per slot, Boole over slots,
i.i.d. fresh-e Wald boot renewals (DESCENT.md).

  - LB = LX = 68.5, LMAX = 260, 38 bins ([68.5, 75), [75, 80), ..., [255, 260)) -> 76 nodes.
  - Each slot draws p_rand from the 4,740 known-log f_9 ideals in [2^14, 2^16), reduces L_{Q, p_rand} and keeps it
    only if the yield certificate (q15_cert16.py) gives mu(L) >= MU_CERT at error DELTA_CERT; draws are charged
    (LLL_NODE + CERT_DRAW) / CERT_PASS_LO. One sieve trial per slot (K_TRIALS = 1), first usable relation, children once.
  - Slot failure <= P_NODE = exp(-MU_CERT) + DELTA_CERT; Boole recurrence W = own + sum_children W,
    PHI = P_NODE + sum_children PHI (no independence across slots).
  - Linear missing-leaf retention max(0.0, 1.0 - miss); BOUND_CELL charged on sieve and direct-ECM branches;
    branch-tree-capped L_min = 2, 12-bit geometric ladder ECM.
  - Boot BOOT16 = (230, 312, 15, 2), BOOT_KEEP_FRAC = 0.90, BALL_PAIR = 0.0; E[runs/log] = 1/(1 - Phi_tree).
"""
import functools
import json
import math
from decimal import Decimal, localcontext
from math import exp, factorial, floor, log, log1p, log2
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
def _miss_side_table(log2_rels, u_main, lx=68.5, main_lb=71.19, orbit=12, p_sieve=68.5, ret_mult=1.0):
    """Precompute strictly-decreasing collection-cut-conditioned lam(b) & right-endpoint weights on [40, lx]."""
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
def miss_side(u4, log2_rels, u_main, lx=68.5, main_lb=71.19, orbit=12, p_sieve=68.5, ret_mult=1.0):
    """Rigorous right-endpoint upper Riemann sum (0.02-bit step) on [40, lx], plus explicit b < 40 analytic tail."""
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
    LB = 68.5
    LX = 68.5
    BP_LO = 14.0
    BP_HI = 16.0
    K_TRIALS = 1
    MU_MIN = 15.0
    MU_CERT = 15.0
    DELTA_CERT = 1.0e-9
    CERT_DRAW = 2.0 ** 34
    CERT_PASS_LO = 2.0 ** -20
    P_EMPTY = 1.0e-4
    P_NODE = P_EMPTY + DELTA_CERT
    PMIN = 8.0
    SLACK = 96.0
    ALPHA_F = 2.74
    ALPHA_G = 0.77
    LOG2U8 = 8.0 * log2(0xd201000000010000)
    LOG2P4 = 4.0 * 380.99700
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
    BOOT_KEEP_FRAC = 0.90
    PRCAP_NODE = 4.0
    PRCAP_BOOT = 1.5
    LOG2_RELS = 62.980
    MAIN_LB = 71.19
    MAIN_UF = 546.0 / MAIN_LB
    MAIN_UG = 670.0 / MAIN_LB
    P_SIEVE_F = 58.0
    P_SIEVE_G = 63.5
    RET_MULT_F = 0.90
    RET_MULT_G = 0.82
    ORBIT = 12
    UDUP = 0.40
    BALL_PAIR = 0.0
    BOUND_CELL = 5000.0
    SIEVE_HIT = 2600.0
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


def bins_of(Pm):
    """Bins from LX up to LMAX: first bin [LX, 75), then 5-bit bins [75, 80), ..., [255, 260)."""
    bins = [(float(Pm.LX), 75.0)]
    l = 75
    while l < Pm.LMAX:
        bins.append((float(l), float(l + Pm.BIN)))
        l += Pm.BIN
    return bins


def node_cost(Pm, l_lo, l_hi, side, n, iA, kA, kB, R):
    """Evaluate ONE single-ball composite L_{Q, p_rand} trial with LINEAR missing-leaf retention
    max(0.0, 1.0 - miss), 100% eager cofactorization, BOUND_CELL on BOTH sieve and direct-ECM branches,
    and branch-tree-capped 12-bit geometric ladder ECM."""
    cells_map = load_cells()
    pairs = pairs_count(Pm)
    cfmax = max(b[1] for b in Pm.QBINS)
    cgmax = max(b[2] for b in Pm.QBINS)
    C_R = cells_map[R]
    S = log2(C_R)
    c_hi = (l_hi + Pm.BP_HI + S) / 16.0
    c_lo = (l_lo + Pm.BP_LO + S) / 16.0
    s_f_hi = (l_lo if side == "f" else 0.0) + Pm.BP_LO
    s_g_hi = l_lo if side == "g" else 0.0
    s_f_lo = (l_hi if side == "f" else 0.0) + Pm.BP_HI
    s_g_lo = l_hi if side == "g" else 0.0

    h = (n - Pm.LX) / Pm.MESH if n > Pm.LX else 0.0
    nA = Pm.LX + iA * h if n > Pm.LX else Pm.LX
    T = kA * nA + kB * n + Pm.SLACK

    pr = 0.0
    sv = 0.0
    for (w, cf, cg), (cf_lo, cg_lo) in zip(Pm.QBINS, Pm.QBINS_LO):
        nf = 72.0 * c_hi + cf - s_f_hi
        ng = 24.0 * c_hi + Pm.LOG2U8 + cg - s_g_hi
        pf = ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC)
        pg = ss_lower2(ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
        pr += w * pf * pg
        nf_lo = 72.0 * c_lo + cf_lo - s_f_lo
        ng_lo = 24.0 * c_lo + Pm.LOG2U8 + cg_lo - s_g_lo
        sv += w * surv_upper(nf_lo, Pm.LX, T) * surv_upper(ng_lo, Pm.LX, T)
    if Pm.QBINS_SUB:
        pri = 0.0
        for cf, cg in Pm.QBINS_SUB:
            nf = 72.0 * c_hi + cf - s_f_hi
            ng = 24.0 * c_hi + Pm.LOG2U8 + cg - s_g_hi
            pri += (
                ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC)
                * ss_lower2(ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
            )
        pr = min(pri / len(Pm.QBINS_SUB), Pm.PRCAP_NODE * pr)
    if pr <= 0.0:
        return None

    nf_max = 72.0 * c_hi + cfmax - s_f_hi
    ng_max = 24.0 * c_hi + Pm.LOG2U8 + cgmax - s_g_hi
    nf_fac = math.ceil(nf_max / 16.0 + kA + kB)
    ng_fac = math.ceil(ng_max / 16.0 + kA + kB)
    det = 1.0 - (1.0 - Pm.ECM_SUCCESS) * ((nf_max + ng_max) / 16.0 + 2.0 * (kA + kB))
    miss = missing(Pm, nf_max, ng_max)
    lin_keep = max(0.0, 1.0 - miss)
    if lin_keep <= 0.0:
        return None
    mu = C_R * pr * Pm.UDUP * det * lin_keep

    n_lp_max = 2 * (kA + kB)
    surv_cost_core = (
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
    trials_direct = C_R * (surv_cost_core + Pm.BOUND_CELL)
    if trials_direct < sieve + survivors * surv_cost_core:
        own_1 = (Pm.LLL_NODE + Pm.CERT_DRAW) / Pm.CERT_PASS_LO + trials_direct
        how = "ecm"
        survivors = C_R
        ret, loss_f, loss_g = _get_counted_fraction(Pm), 0.0, 0.0
        mu *= ret
    else:
        own_1 = (Pm.LLL_NODE + Pm.CERT_DRAW) / Pm.CERT_PASS_LO + sieve + survivors * surv_cost_core
        how = "sieve"
        ret, loss_f, loss_g = _get_retention(Pm, nf_max, ng_max, n, iA, kA, kB)
        mu *= ret
    if mu < Pm.MU_MIN:
        return None
    own = float(Pm.K_TRIALS) * own_1
    return dict(
        own_1=own_1,
        own=own,
        own_log2=log2(own),
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
        lin_keep=lin_keep,
        ret=ret,
        loss_f=loss_f,
        loss_g=loss_g,
        pr_log2=log2(pr),
        surv_log2=log2(survivors),
        sieve_log2=log2(sieve),
        survcost_log2=log2(survivors * surv_cost_core),
        nf=nf_max,
        ng=ng_max,
        det=det,
        T=T,
        cofac_frac=1.0,
        eager=True,
    )


def evaluate(Pm, PLAN, BOOT):
    """Exact Pure k=1 Zero-Retry Boole + Unbounded i.i.d. Fresh-e Wald Boot evaluation."""
    cells_map = load_cells()
    bins = bins_of(Pm)
    W, PHI, NSL, DET = {}, {}, {}, {}

    def wmax(side, lim):
        return max([W[b, side] for b in bins if (b, side) in W and b[0] < lim - 1e-9] + [0.0])

    def phimax(side, lim):
        return max([PHI[b, side] for b in bins if (b, side) in PHI and b[0] < lim - 1e-9] + [0.0])

    def nsmax(side, lim):
        return max([NSL[b, side] for b in bins if (b, side) in NSL and b[0] < lim - 1e-9] + [0])

    for b in bins:
        for side in ("f", "g"):
            n, iA, kA, kB, R = PLAN[b[0], side]
            h_tmp = (n - Pm.LX) / Pm.MESH if n > Pm.LX else 0.0
            nA_tmp = Pm.LX + iA * h_tmp if n > Pm.LX else Pm.LX
            if not ((kA == 0 or nA_tmp < b[0] - 1e-9) and (kB == 0 or n < b[0] - 1e-9)):
                raise ValueError(f"non-decreasing plan at {b} {side}: n={n}, nA={nA_tmp}")
            d = node_cost(Pm, b[0], b[1], side, float(n), iA, kA, kB, R)
            if d is None:
                raise ValueError(f"infeasible plan entry {b} {side}")
            nA = d["nA"]
            sub_w = kA * (wmax("f", nA) + wmax("g", nA)) + kB * (wmax("f", float(n)) + wmax("g", float(n)))
            sub_p = kA * (phimax("f", nA) + phimax("g", nA)) + kB * (phimax("f", float(n)) + phimax("g", float(n)))
            sub_n = kA * (nsmax("f", nA) + nsmax("g", nA)) + kB * (nsmax("f", float(n)) + nsmax("g", float(n)))
            W[b, side] = d["own"] + sub_w
            PHI[b, side] = Pm.P_NODE + sub_p
            NSL[b, side] = 1 + sub_n
            DET[b, side] = dict(
                d,
                p0=Pm.P_NODE,
                sub_log2=log2(sub_w) if sub_w > 0 else None,
                W_log2=log2(W[b, side]),
                phi_dag=PHI[b, side],
                n_slots=NSL[b, side],
            )

    n1, iA_b, kA_b, Rb = BOOT
    if n1 > Pm.LMAX:
        raise ValueError(f"boot window n1={n1} above LMAX={Pm.LMAX}")
    C_Rb = cells_map[Rb]
    Sb = log2(C_Rb)
    cb = (Pm.LOG2P4 + Sb) / 16.0
    hb = (n1 - Pm.LX) / Pm.MESH
    nA_b = Pm.LX + iA_b * hb
    m_max = 24.0 * cb + Pm.LOG2U8 + max(x[1] for x in Pm.BOOT_BINS)
    pr_b = sum(
        w * ss_lower2(24.0 * cb + Pm.LOG2U8 + cg, Pm.LX, float(n1), iA_b, kA_b, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
        for w, cg in Pm.BOOT_BINS
    )
    if Pm.BOOT_SUB:
        pri_b = sum(
            ss_lower2(24.0 * cb + Pm.LOG2U8 + cg, Pm.LX, float(n1), iA_b, kA_b, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
            for cg in Pm.BOOT_SUB
        )
        pr_b = min(pri_b / len(Pm.BOOT_SUB), Pm.PRCAP_BOOT * pr_b)
    if pr_b <= 0.0:
        raise ValueError("infeasible boot pr_b <= 0")
    boot_n_factors = math.ceil(m_max / 16.0 + kA_b + 1.0)
    det_b = 1.0 - (1.0 - Pm.ECM_SUCCESS) * (m_max / 16.0 + kA_b + 1.0)
    mb = missing(Pm, 0.0, m_max)
    lin_keep_b = max(0.0, 1.0 - mb)
    if lin_keep_b <= 0.0:
        raise ValueError("infeasible boot lin_keep_b <= 0")
    prob_b = pr_b * Pm.BOOT_KEEP_FRAC * det_b * Pm.UDUP * _get_counted_fraction(Pm) * lin_keep_b
    trial_b = (
        Pm.LLL_BOOT
        + Pm.NORMS
        + Pm.TDIV
        + ecm_instr(nA_b, m_max, Pm.ECM_SUCCESS, n_factors=boot_n_factors)
        + (kA_b + 1) * (Pm.PRIME + Pm.IDEAL)
    )
    boot_own_1run = trial_b / prob_b
    sub_w_1run = kA_b * wmax("g", nA_b) + wmax("g", float(n1))
    phi_tree = kA_b * phimax("g", nA_b) + phimax("g", float(n1))
    n_slots_tree = kA_b * nsmax("g", nA_b) + nsmax("g", float(n1))
    if phi_tree >= 0.99:
        raise ValueError(f"phi_tree >= 0.99: {phi_tree}")
    e_runs = 1.0 / (1.0 - phi_tree)
    one_log = e_runs * (boot_own_1run + sub_w_1run)
    two_log_u = log2(2.0 * one_log / UNIT)
    boot_u = log2(2.0 * e_runs * boot_own_1run / UNIT)
    sub_u = log2(2.0 * e_runs * sub_w_1run / UNIT) if sub_w_1run > 0 else 0.0
    boot = dict(
        n1=n1,
        nA=nA_b,
        kA=kA_b,
        R=Rb,
        S=Sb,
        m=m_max,
        boot_keep_frac=Pm.BOOT_KEEP_FRAC,
        boot_n_factors=boot_n_factors,
        pr_log2=log2(pr_b),
        trials_1run_log2=log2(1.0 / prob_b),
        trial_log2=log2(trial_b),
        balls=e_runs / prob_b,
        miss=mb,
        lin_keep=lin_keep_b,
        own_log2=log2(e_runs * boot_own_1run),
        sub_log2=log2(e_runs * sub_w_1run),
        boot_u=boot_u,
        sub_u=sub_u,
        phi_tree=phi_tree,
        e_runs=e_runs,
        n_slots_tree=n_slots_tree,
        phi=0.0,
        one_log_log2=log2(one_log),
    )
    return dict(
        W=W,
        PHI=PHI,
        NSL=NSL,
        DET=DET,
        bins=bins,
        boot=boot,
        two_log_units=two_log_u,
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
                    n=d["n"],
                    nA=round(d["nA"], 2),
                    kA=d["kA"],
                    kB=d["kB"],
                    R=d["R"],
                    ret=round(d["ret"], 4),
                    loss_f=d["loss_f"],
                    loss_g=d["loss_g"],
                    how=d["how"],
                    mu=round(d["mu"], 2),
                    lin_keep=round(d["lin_keep"], 4),
                    own_log2=round(d["own_log2"], 3),
                    W_log2=round(d["W_log2"], 3),
                    phi_dag=round(d["phi_dag"], 6),
                    n_slots=d["n_slots"],
                    nf=round(d["nf"], 1),
                    ng=round(d["ng"], 1),
                    T=round(d["T"], 1),
                )
            )
    return rows


# Shipped plan PLAN16 (MESH=384, TRUNC=0.0):
_PLAN = {
    68.5: ((68.5, 384, 0, 0, 56), (68.5, 384, 0, 0, 56)),
    75.0: ((68.5, 384, 0, 0, 64), (68.5, 384, 0, 0, 64)),
    80.0: ((80.0, 192, 1, 0, 48), (75.0, 0, 0, 1, 48)),
    85.0: ((75.0, 0, 0, 1, 56), (75.0, 0, 0, 1, 56)),
    90.0: ((90.0, 96, 2, 0, 56), (90.0, 96, 2, 0, 56)),
    95.0: ((80.0, 192, 2, 0, 64), (80.0, 192, 2, 0, 64)),
    100.0: ((90.0, 96, 1, 1, 48), (90.0, 96, 1, 1, 48)),
    105.0: ((90.0, 96, 1, 1, 56), (85.0, 144, 1, 1, 64)),
    110.0: ((95.0, 144, 1, 1, 56), (90.0, 96, 2, 1, 56)),
    115.0: ((95.0, 144, 1, 1, 64), (110.0, 240, 2, 0, 56)),
    120.0: ((120.0, 192, 2, 0, 64), (95.0, 144, 2, 1, 56)),
    125.0: ((110.0, 288, 2, 0, 64), (95.0, 240, 2, 1, 56)),
    130.0: ((100.0, 192, 2, 1, 56), (100.0, 240, 2, 1, 56)),
    135.0: ((100.0, 192, 2, 1, 64), (100.0, 240, 2, 1, 64)),
    140.0: ((105.0, 192, 2, 1, 64), (110.0, 192, 2, 1, 64)),
    145.0: ((115.0, 192, 2, 1, 56), (110.0, 240, 2, 1, 64)),
    150.0: ((120.0, 192, 2, 1, 64), (130.0, 144, 2, 1, 64)),
    155.0: ((130.0, 144, 2, 1, 64), (135.0, 144, 2, 1, 64)),
    160.0: ((130.0, 192, 2, 1, 64), (135.0, 192, 2, 1, 64)),
    165.0: ((135.0, 192, 2, 1, 64), (145.0, 192, 2, 1, 64)),
    170.0: ((145.0, 192, 2, 1, 64), (150.0, 192, 2, 1, 64)),
    175.0: ((160.0, 144, 2, 1, 64), (150.0, 192, 3, 1, 40)),
    180.0: ((150.0, 288, 4, 1, 12), (150.0, 192, 4, 1, 28)),
    185.0: ((155.0, 288, 4, 1, 12), (160.0, 192, 3, 1, 40)),
    190.0: ((160.0, 288, 4, 1, 12), (160.0, 192, 3, 1, 48)),
    195.0: ((155.0, 336, 4, 1, 12), (160.0, 192, 3, 1, 56)),
    200.0: ((165.0, 288, 4, 1, 12), (165.0, 192, 4, 1, 32)),
    205.0: ((165.0, 288, 4, 1, 14), (165.0, 240, 3, 1, 48)),
    210.0: ((165.0, 288, 4, 1, 16), (160.0, 288, 4, 1, 24)),
    215.0: ((160.0, 288, 4, 1, 20), (165.0, 288, 4, 1, 24)),
    220.0: ((165.0, 288, 4, 1, 20), (165.0, 288, 4, 1, 28)),
    225.0: ((180.0, 144, 3, 1, 112), (170.0, 240, 4, 1, 40)),
    230.0: ((165.0, 288, 4, 1, 24), (175.0, 192, 4, 1, 64)),
    235.0: ((170.0, 192, 3, 1, 128), (170.0, 240, 4, 1, 56)),
    240.0: ((175.0, 192, 3, 1, 128), (170.0, 240, 4, 1, 64)),
    245.0: ((175.0, 240, 3, 1, 96), (170.0, 240, 4, 1, 80)),
    250.0: ((175.0, 240, 3, 1, 128), (170.0, 240, 4, 1, 96)),
    255.0: ((165.0, 288, 4, 1, 48), (170.0, 240, 4, 1, 112)),
}
PLAN16 = {(float(b), s): v[i] for b, v in _PLAN.items() for i, s in enumerate("fg")}
BOOT16 = (230, 312, 15, 2)
