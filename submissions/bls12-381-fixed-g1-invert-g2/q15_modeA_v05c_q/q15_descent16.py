"""q15_descent16.py — Pure Mode-A k=1 Zero-Retry Boole + Unbounded i.i.d. Fresh-e Wald Boot
16D Q(zeta_15) Descent Evaluator.

Pure Mode-A architecture (zero internal retries, single special-q lattice L_Q, unconditional Boole):
  - Single special-q lattice L_Q (index q in Z^16) per slot; LB = LX, LMAX = 260, 5-bit bins.
  - Pure k = 1 zero-retry node failure bound:
      p_0(b, s) = max(P0_FLOOR, (1.0 + NB_ALPHA * mu(b, s) / NB_R) ** (-NB_R))
    with NB_R = 1.0, NB_ALPHA = 0.5, P0_FLOOR = 0.001, generous MU_MIN >= 500.
  - Bottom-up unconditional Boole tree recursion over c_det deterministic slots:
      c_det(b, s) = 1 + kA*(cdet_max(f, nA) + cdet_max(g, nA)) + kB*(cdet_max(f, n) + cdet_max(g, n))
      phi(b, s)   = p_0(b, s) + kA*(phi_max(f, nA) + phi_max(g, nA)) + kB*(phi_max(f, n) + phi_max(g, n))
      W(b, s)     = own_1(b, s) + kA*(w_max(f, nA) + w_max(g, nA)) + kB*(w_max(f, n) + w_max(g, n))
  - Unbounded i.i.d. fresh-e Wald boot (BOOT16 = (230, 288, 12, 2), ZERO B_max/T_boot_max):
      c_det_one   = kA_b * cdet_max(g, nA_b) + cdet_max(g, n1)
      Phi_tree    = kA_b * phi_max(g, nA_b) + phi_max(g, n1) <= MAX_PHI_TREE = 0.60
      e_boot_runs = 1.0 / (1.0 - Phi_tree),  boot_phi = 0.0,  p_descent_fail = 0.0
      one_log     = e_boot_runs * (boot_own_1run + ch_boot)
  - Explicit BOUND_CELL charged on BOTH sieve and direct-ECM code paths, linear union-bound missing-leaf
    retention max(0.0, 1.0 - miss), and branch-tree-capped L_min = 2, 12-bit geometric ladder ECM.
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
_RHO_COEFS_REV = tuple(tuple(reversed(c)) if c is not None else None for c in RHO_COEFS)


def rho_series(u):
    if u <= 1.0:
        return 1.0
    k = math.ceil(u)
    if k > KMAX:
        return 0.0
    x, value = k - u, 0.0
    for c in _RHO_COEFS_REV[k]:
        value = value * x + c
    return value if value > 0.0 else 0.0


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
def ss_cum2(m, base, n, iA, mesh, alpha, trunc, kmax=KA_MAX, need_b=1):
    rows, h = window_rows(n, base, iA, mesh, kmax)
    out = [[0.0, 0.0] for _ in range(kmax + 1)]
    tab, ntab, ceil = RHO_TAB_LO, len(RHO_TAB_LO), math.ceil
    sp = _SUMPROD
    inv_br = 1.0 / (base * RSTEP)
    hk = h * inv_br
    tr_span = int(trunc * base / h) + 1 if (trunc and h > 0) else 0
    jBs = (0, 1) if need_b else (0,)
    j_lim = kmax + (2 if need_b else 1)
    # base * inv_br == 1 / RSTEP == 1024 exactly; precompute j=0 ceiling indices once across all j.
    step_j = int(round(1.0 / RSTEP))
    k0_0 = (m + alpha) * inv_br
    plan = []
    g_lo, g_hi = 10**9, -1
    for j in range(j_lim):
        items = []
        ifeas_j = int(floor((m - j * h - j * base) / h)) if h > 0 else 10**9
        k0 = k0_0 - j * step_j
        i_ntab = int((k0 - (ntab - 1)) / hk) if hk > 0.0 else 0
        for jB in jBs:
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
        if t0 < g_lo:
            g_lo = t0
        if t1 > g_hi:
            g_hi = t1
        plan.append((j, k0, t0, t1, items))
    if plan:
        idx0 = [ceil(k0_0 - i * hk) for i in range(g_lo, g_hi + 1)]
        for j, k0, t0, t1, items in plan:
            sh = j * step_j
            sl = idx0[t0 - g_lo : t1 - g_lo + 1]
            if k0 - t1 * hk <= 0.0 or k0 - t0 * hk >= ntab:
                arr = [1.0 if (k := x - sh) <= 0 else (tab[k] if k < ntab else 0.0) for x in sl]
            else:
                arr = [tab[x - sh] for x in sl]
            for jA, jB, row, imin, ifeas in items:
                out[jA][jB] = sp(row[imin : ifeas + 1], arr[imin - t0 : ifeas - t0 + 1])
    cum = [[0.0, 0.0] for _ in range(kmax + 1)]
    for jB in jBs:
        acc = 0.0
        for jA in range(kmax + 1):
            acc += out[jA][jB]
            cum[jA][jB] = acc
    return tuple(tuple(x) for x in cum)


def ss_lower2(m, base, n, iA, kA, kB, mesh, alpha, trunc=0.0):
    m = math.ceil(m)
    if n <= base:
        return rho_lo((m + alpha) / base)
    nb = 1 if kB else 0
    cum = ss_cum2(m, base, n, iA, mesh, alpha, trunc, kA, nb)
    return cum[kA][0] + (cum[kA][1] if nb else 0.0)


@functools.lru_cache(None)
def _ss_lower2_tau_c(m, base, n, iA, kA, kB, mesh, alpha, tau):
    rows, h = window_rows(n, base, iA, mesh, kA)
    tab, ntab, ceil = RHO_TAB_LO, len(RHO_TAB_LO), math.ceil
    inv_br = 1.0 / (base * RSTEP)
    hk = h * inv_br
    sp = _SUMPROD
    step_j = int(round(1.0 / RSTEP))
    k0_0 = (m + alpha) * inv_br
    cap = m if m < tau else tau
    jBs = (0, 1) if kB else (0,)
    j_lim = kA + (2 if kB else 1)
    plan = []
    g_lo, g_hi = 10**9, -1
    for j in range(j_lim):
        items = []
        ifeas_j = int(floor((cap - j * h - j * base) / h)) if h > 0 else 10**9
        k0 = k0_0 - j * step_j
        i_ntab = int((k0 - (ntab - 1)) / hk) if hk > 0.0 else 0
        for jB in jBs:
            jA = j - jB
            if not (0 <= jA <= kA):
                continue
            row, r_lo, r_hi = rows[jA, jB]
            if r_lo >= r_hi:
                continue
            ifeas = r_hi - 1 if (r_hi - 1) < ifeas_j else ifeas_j
            if ifeas < r_lo:
                continue
            imin = r_lo if r_lo > i_ntab else i_ntab
            if imin <= ifeas:
                items.append((row, imin, ifeas))
        if not items:
            continue
        t0 = min(it[1] for it in items)
        t1 = max(it[2] for it in items)
        if t0 < g_lo:
            g_lo = t0
        if t1 > g_hi:
            g_hi = t1
        plan.append((j, k0, t0, t1, items))
    tot = 0.0
    if plan:
        idx0 = [ceil(k0_0 - i * hk) for i in range(g_lo, g_hi + 1)]
        for j, k0, t0, t1, items in plan:
            sh = j * step_j
            sl = idx0[t0 - g_lo : t1 - g_lo + 1]
            if k0 - t1 * hk <= 0.0 or k0 - t0 * hk >= ntab:
                arr = [1.0 if (k := x - sh) <= 0 else (tab[k] if k < ntab else 0.0) for x in sl]
            else:
                arr = [tab[x - sh] for x in sl]
            for row, imin, ifeas in items:
                tot += sp(row[imin : ifeas + 1], arr[imin - t0 : ifeas - t0 + 1])
    return tot


def ss_lower2_tau(m, base, n, iA, kA, kB, mesh, alpha, tau):
    """Exact ss_lower2 restricted to window-prime sum j*base + (i+j)*h <= tau (2-stage boot cofactor cap)."""
    m = math.ceil(m)
    if n <= base:
        return rho_lo((m + alpha) / base)
    return _ss_lower2_tau_c(m, base, n, iA, kA, kB, mesh, alpha, tau)


def surv_upper(m, base, T):
    """Strict upper bound on P(residual cofactor <= 2^T after stripping <=base factors) for all m >= 0, T >= 0.
    At T == 0 (pure base-smooth check R == 1), returns rho_hi(floor(m)/base), NEVER 0.0."""
    if T < 0.0:
        return 0.0
    if T >= m:
        return 1.0
    if T == 0.0:
        return min(1.0, rho_hi(floor(m) / base))
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


CELLS = {}


def load_cells():
    global CELLS
    import q15_cells16
    if q15_cells16.CELLS_16D:
        for k, v in q15_cells16.CELLS_16D.items():
            if k in CELLS:
                assert CELLS[k] == v, f"CELLS mismatch at R={k}"
            else:
                CELLS[k] = v
    elif not CELLS:
        cells, _, _ = q15_cells16.compute_cells_16d(q15_cells16.RADII_16D)
        CELLS.update(cells)
    return CELLS


class P:
    # Mode-A pure k=1 zero-retry Boole + unbounded i.i.d. fresh-e Wald boot parameters
    LB = 68.5
    LX = 68.5
    NB_R = 1.0
    NB_ALPHA = 0.5
    P0_FLOOR = 0.001
    DISC = 0.0001
    MU_MIN = 800.0
    MAX_PHI_TREE = 0.60
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
    TAU_BOOT = 2000
    BOOT_LO_OFF = -270.0
    PRCAP_NODE = 4.0
    PRCAP_BOOT = 1.5
    LOG2_RELS = 62.952
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
    BOUND_CELL = 4000.0
    SIEVE_HIT = 1100.0
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


def p0_modeA(mu, Pm=P):
    """Pure Mode-A k = 1 zero-retry slot failure probability: max(P0_FLOOR, (1 + NB_ALPHA * mu / NB_R)^(-NB_R)) + DISC."""
    return max(Pm.P0_FLOOR, (1.0 + Pm.NB_ALPHA * mu / Pm.NB_R) ** (-Pm.NB_R)) + Pm.DISC


def load_stats(stats_or_path, Pm=None):
    """Load QBINS, QBINS_LO, QBINS_SUB, BOOT_BINS, BOOT_SUB, CELLS, UDUP, tower, and retain16 tables from q15_stats."""
    st = json.loads(Path(stats_or_path).read_text()) if isinstance(stats_or_path, (str, Path)) else stats_or_path
    tgt = Pm if Pm is not None else P
    if "cells_16d" in st:
        CELLS.update({int(k): int(v) for k, v in st["cells_16d"].items()})
    load_cells()
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
    first_hi = 75.0 if Pm.LX < 75.0 else float(int(Pm.LX // Pm.BIN) * Pm.BIN + Pm.BIN)
    bins = [(float(Pm.LX), first_hi)]
    l = first_hi
    while l < Pm.LMAX - 1e-9:
        bins.append((float(l), float(min(Pm.LMAX, l + Pm.BIN))))
        l += Pm.BIN
    return bins


def node_cost(Pm, l_lo, l_hi, side, n, iA, kA, kB, R):
    """Evaluate ONE pure k=1 zero-retry single-special-q (L_Q) descent slot with linear missing-leaf retention
    max(0.0, 1.0 - miss), 100% eager cofactorization, BOUND_CELL on BOTH sieve and direct-ECM branches,
    and branch-tree-capped L_min=2, 12-bit geometric ladder ECM."""
    cells_map = load_cells()
    pairs = pairs_count(Pm)
    cfmax = max(b[1] for b in Pm.QBINS)
    cgmax = max(b[2] for b in Pm.QBINS)
    C_R = cells_map[R]
    S = log2(C_R)
    c_hi = (l_hi + S) / 16.0
    c_lo = (l_lo + S) / 16.0
    s_f_hi = l_lo if side == "f" else 0.0
    s_g_hi = l_lo if side == "g" else 0.0
    s_f_lo = l_hi if side == "f" else 0.0
    s_g_lo = l_hi if side == "g" else 0.0

    h = (n - Pm.LX) / Pm.MESH if n > Pm.LX else 0.0
    nA = Pm.LX + iA * h if n > Pm.LX else Pm.LX
    T = kA * nA + kB * n + Pm.SLACK
    T0 = kA * nA + kB * n  # Exact post-ecm(LX) zero-slack cofactor bound (ZERO sieve slack!)

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
    # Single global pre-ECM lower-norm floor (from min(QBINS_LO): 15b below observed_min cf, 9b below observed_min cg):
    # bound_cell_check_and_cap enforces the upper caps (nf_max, ng_max); exact integer norms computed in Pm.NORMS
    # (already charged as the first step of e_s1, BEFORE Pm.TDIV and Stage-1 ecm(LX)) deterministically enforce
    # the lower-norm floors log2(|N_f|/q_f) >= nf_floor and log2(|N_g|/q_g) >= ng_floor with zero approximation gap.
    # Cells below floor are rejected deterministically before Stage-1 ecm(LX) (0.0 acceptance below floor);
    # on admitted cells above floor, under declared work heuristics H1 / H-EA, sv0 = surv_upper(nf_floor, LX, T0) * surv_upper(ng_floor, LX, T0)
    # monotonically bounds Stage-1 survival to exact residual cofactor (R_f <= 2^T0, R_g <= 2^T0).
    nf_floor = 72.0 * c_lo + min(x[0] for x in Pm.QBINS_LO) - s_f_lo
    ng_floor = 24.0 * c_lo + Pm.LOG2U8 + min(x[1] for x in Pm.QBINS_LO) - s_g_lo
    sv0 = surv_upper(nf_floor, Pm.LX, T0) * surv_upper(ng_floor, Pm.LX, T0)
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
    e_s1 = (
        Pm.NORMS
        + Pm.TDIV
        + ecm_instr(Pm.LX, nf_max, Pm.ECM_SUCCESS, n_factors=nf_fac)
        + ecm_instr(Pm.LX, ng_max, Pm.ECM_SUCCESS, n_factors=ng_fac)
    )
    # 2-stage post-ecm(LX) zero-slack cofactor gate:
    #   - If kA+kB == 0 (T0 == 0): only R_f == 1 and R_g == 1 survives -> ZERO Stage-2 ECM (e_s2 = 0.0).
    #   - Else: only survivors with log2(R_f) <= T0 and log2(R_g) <= T0 (rate <= sv0 + pr under H1 / H-EA) run
    #           deterministic 2.0 * ecm_instr(nA, T0, ECM_SUCCESS) (lo_bits = 16.0) + lp*(PRIME + IDEAL).
    if kA + kB == 0:
        e_s2 = 0.0
    else:
        e_s2 = (
            2.0 * ecm_instr(nA, T0, Pm.ECM_SUCCESS, n_factors=kA + kB)
            + max(2, n_lp_max) * Pm.PRIME
            + max(1, n_lp_max) * Pm.IDEAL
        )
    survivors = C_R * sv + C_R * pr
    survivors_s2 = C_R * sv0 + C_R * pr
    surv_cost_tot = survivors * e_s1 + survivors_s2 * e_s2
    hits_pc = 2.0 * log(Pm.LX / Pm.PMIN)
    sieve = C_R * (hits_pc * Pm.SIEVE_HIT + Pm.SCAN + Pm.ZERO + Pm.BOUND_CELL) + pairs * Pm.PAIR_SETUP
    trials_direct = C_R * (Pm.BOUND_CELL + e_s1 + (sv0 + pr) * e_s2)
    if trials_direct < sieve + surv_cost_tot:
        own_1 = Pm.LLL_NODE + trials_direct
        how = "ecm"
        survivors = C_R
        ret, loss_f, loss_g = _get_counted_fraction(Pm), 0.0, 0.0
        mu *= ret
    else:
        own_1 = Pm.LLL_NODE + sieve + surv_cost_tot
        how = "sieve"
        ret, loss_f, loss_g = _get_retention(Pm, nf_max, ng_max, n, iA, kA, kB)
        mu *= ret
    if mu < Pm.MU_MIN:
        return None
    p0 = p0_modeA(mu, Pm)
    return dict(
        own_1=own_1,
        own_1_log2=log2(own_1),
        n=n,
        iA=iA,
        nA=nA,
        kA=kA,
        kB=kB,
        R=R,
        S=S,
        mu=mu,
        p0=p0,
        how=how,
        miss=miss,
        lin_keep=lin_keep,
        ret=ret,
        loss_f=loss_f,
        loss_g=loss_g,
        pr_log2=log2(pr),
        surv_log2=log2(survivors),
        sieve_log2=log2(sieve),
        survcost_log2=log2(max(1.0, surv_cost_tot)),
        nf=nf_max,
        ng=ng_max,
        nf_floor=nf_floor,
        ng_floor=ng_floor,
        det=det,
        T=T,
        T0=T0,
        sv0=sv0,
        cofac_frac=1.0,
        eager=True,
    )


def evaluate(Pm, PLAN, BOOT):
    """Pure Mode-A k=1 zero-retry unconditional Boole DAG + Unbounded i.i.d. Fresh-e Wald Boot."""
    cells_map = load_cells()
    bins = bins_of(Pm)
    W_NODE, PHI_NODE, C_DET, DET = {}, {}, {}, {}

    def w_max(side, lim):
        return max([W_NODE[b, side] for b in bins if (b, side) in W_NODE and b[0] < lim - 1e-9] + [0.0])

    def phi_max(side, lim):
        return max([PHI_NODE[b, side] for b in bins if (b, side) in PHI_NODE and b[0] < lim - 1e-9] + [0.0])

    def cdet_max(side, lim):
        return max([C_DET[b, side] for b in bins if (b, side) in C_DET and b[0] < lim - 1e-9] + [0])

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
            c_det = (
                1
                + (kA * (cdet_max("f", nA) + cdet_max("g", nA)) if kA else 0)
                + (kB * (cdet_max("f", float(n)) + cdet_max("g", float(n))) if kB else 0)
            )
            phi_n = (
                d["p0"]
                + (kA * (phi_max("f", nA) + phi_max("g", nA)) if kA else 0.0)
                + (kB * (phi_max("f", float(n)) + phi_max("g", float(n))) if kB else 0.0)
            )
            w_n = (
                d["own_1"]
                + (kA * (w_max("f", nA) + w_max("g", nA)) if kA else 0.0)
                + (kB * (w_max("f", float(n)) + w_max("g", float(n))) if kB else 0.0)
            )
            C_DET[b, side] = c_det
            PHI_NODE[b, side] = phi_n
            W_NODE[b, side] = w_n
            DET[b, side] = dict(
                d,
                c_det=c_det,
                phi_node=phi_n,
                w_node_log2=log2(w_n / UNIT),
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
    tau_b = getattr(Pm, "TAU_BOOT", None)
    boot_lo_off = getattr(Pm, "BOOT_LO_OFF", -270.0)
    m_min = 24.0 * cb + Pm.LOG2U8 + boot_lo_off
    if tau_b is None:
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
    else:
        pr_b = sum(
            w * ss_lower2_tau(24.0 * cb + Pm.LOG2U8 + cg, Pm.LX, float(n1), iA_b, kA_b, 1, Pm.MESH, Pm.ALPHA_G, tau_b)
            for w, cg in Pm.BOOT_BINS
        )
        if Pm.BOOT_SUB:
            pri_b = sum(
                ss_lower2_tau(24.0 * cb + Pm.LOG2U8 + cg, Pm.LX, float(n1), iA_b, kA_b, 1, Pm.MESH, Pm.ALPHA_G, tau_b)
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
    if tau_b is None:
        sv_b = 1.0
        trial_b = (
            Pm.LLL_BOOT
            + Pm.NORMS
            + Pm.TDIV
            + ecm_instr(nA_b, m_max, Pm.ECM_SUCCESS, n_factors=boot_n_factors)
            + (kA_b + 1) * (Pm.PRIME + Pm.IDEAL)
        )
    else:
        # 2-stage tau-capped boot with explicit two-sided pre-ECM norm window [m_min, m_max]:
        # Before Stage-1 ECM, reject any cell with log2(|N_g|) < m_min or log2(|N_g|) > m_max
        # (acceptance probability = 0.0 identically for all m < m_min, including 0 < m < 515.01).
        # On admitted cells m in [m_min, m_max], Stage-1 ecm(LX, m_max) strips <=68.5b factors,
        # and Stage-2 runs only if residual cofactor log2(R_g) <= tau_b = 2000; since m - tau_b >= m_min - tau_b = 515.01,
        # Stage-1 survival on admitted norms is monotonically bounded by surv_upper(m_min, LX, tau_b) under H1 / H-EA
        # (and is deterministically 0 for m < m_min by the exact-NORMS pre-screen).
        sv_b = surv_upper(m_min, Pm.LX, tau_b)
        ecm_s1 = ecm_instr(Pm.LX, m_max, Pm.ECM_SUCCESS, n_factors=boot_n_factors)
        ecm_s2 = ecm_instr(nA_b, tau_b, Pm.ECM_SUCCESS, n_factors=kA_b + 1)
        trial_b = (
            Pm.LLL_BOOT
            + Pm.NORMS
            + Pm.TDIV
            + ecm_s1
            + (sv_b + pr_b) * (ecm_s2 + (kA_b + 1) * (Pm.PRIME + Pm.IDEAL))
        )
    boot_own_1run = trial_b / prob_b
    ch_boot = kA_b * w_max("g", nA_b) + w_max("g", float(n1))
    phi_tree = kA_b * phi_max("g", nA_b) + phi_max("g", float(n1))
    c_det_one = kA_b * cdet_max("g", nA_b) + cdet_max("g", float(n1))
    if phi_tree > Pm.MAX_PHI_TREE + 1e-12:
        raise ValueError(f"Phi_tree={phi_tree:.6f} exceeds MAX_PHI_TREE={Pm.MAX_PHI_TREE}")
    e_boot_runs = 1.0 / (1.0 - phi_tree)
    one_log = e_boot_runs * (boot_own_1run + ch_boot)
    two_log_u = log2(2.0 * one_log / UNIT)
    boot_u = log2(2.0 * e_boot_runs * boot_own_1run / UNIT)
    sub_u = log2(2.0 * e_boot_runs * ch_boot / UNIT) if ch_boot > 0 else 0.0
    boot = dict(
        n1=n1,
        nA=nA_b,
        kA=kA_b,
        R=Rb,
        S=Sb,
        m=m_max,
        m_min=m_min,
        tau_boot=tau_b,
        boot_lo_off=boot_lo_off if tau_b is not None else None,
        surv_b=sv_b,
        boot_keep_frac=Pm.BOOT_KEEP_FRAC,
        boot_n_factors=boot_n_factors,
        pr_log2=log2(pr_b),
        trials_1run_log2=log2(1.0 / prob_b),
        trial_log2=log2(trial_b),
        e_trials_log=e_boot_runs / prob_b,
        miss=mb,
        lin_keep=lin_keep_b,
        own_log2=log2(e_boot_runs * boot_own_1run),
        sub_log2=log2(e_boot_runs * ch_boot),
        boot_u=boot_u,
        sub_u=sub_u,
        phi_tree=phi_tree,
        c_det_one=c_det_one,
        e_boot_runs=e_boot_runs,
        p0_max=max(d["p0"] for d in DET.values()),
        mu_min=min(d["mu"] for d in DET.values()),
        phi=0.0,
        one_log_log2=log2(one_log),
    )
    return dict(
        W_NODE=W_NODE,
        PHI_NODE=PHI_NODE,
        C_DET=C_DET,
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
                    p0=round(d["p0"], 6),
                    c_det=d["c_det"],
                    phi_node=round(d["phi_node"], 6),
                    lin_keep=round(d["lin_keep"], 4),
                    own_1_log2=round(d["own_1_log2"], 3),
                    w_node_log2=round(d["w_node_log2"], 3),
                    nf=round(d["nf"], 1),
                    ng=round(d["ng"], 1),
                    T=round(d["T"], 1),
                )
            )
    return rows


# Shipped Mode-A pure k=1 zero-retry Boole + 2-stage post-ecm(LX) zero-slack T0 + 2-stage tau=2000 Wald boot plan:
PLAN_TAG = "MODEA_V05C_LO270_SH1100_DISC_MU800"
_PLAN = {
    68.5: ((68.5, 384, 0, 0, 56), (68.5, 384, 0, 0, 56)),
    75.0: ((68.5, 384, 0, 0, 64), (68.5, 384, 0, 0, 64)),
    80.0: ((68.5, 384, 0, 0, 72), (68.5, 384, 0, 0, 72)),
    85.0: ((68.5, 384, 0, 0, 80), (80.0, 192, 1, 0, 60)),
    90.0: ((90.0, 192, 1, 0, 60), (80.0, 0, 0, 1, 60)),
    95.0: ((95.0, 144, 1, 0, 72), (90.0, 192, 1, 0, 72)),
    100.0: ((80.0, 0, 0, 1, 80), (90.0, 288, 1, 0, 76)),
    105.0: ((90.0, 192, 2, 0, 68), (80.0, 384, 2, 0, 68)),
    110.0: ((80.0, 384, 2, 0, 76), (80.0, 384, 2, 0, 76)),
    115.0: ((85.0, 384, 2, 0, 72), (85.0, 384, 2, 0, 76)),
    120.0: ((85.0, 384, 2, 0, 80), (85.0, 240, 2, 1, 72)),
    125.0: ((85.0, 240, 2, 1, 76), (85.0, 240, 2, 1, 80)),
    130.0: ((85.0, 384, 3, 0, 76), (85.0, 384, 3, 0, 80)),
    135.0: ((100.0, 192, 2, 1, 68), (100.0, 192, 2, 1, 72)),
    140.0: ((100.0, 192, 2, 1, 76), (100.0, 192, 2, 1, 80)),
    145.0: ((110.0, 144, 2, 1, 72), (110.0, 144, 2, 1, 80)),
    150.0: ((110.0, 144, 2, 1, 84), (100.0, 384, 3, 0, 80)),
    155.0: ((100.0, 384, 3, 0, 84), (105.0, 384, 3, 0, 80)),
    160.0: ((105.0, 384, 3, 0, 80), (115.0, 288, 2, 1, 88)),
    165.0: ((115.0, 288, 2, 1, 84), (125.0, 240, 2, 1, 80)),
    170.0: ((125.0, 240, 2, 1, 80), (115.0, 240, 3, 1, 72)),
    175.0: ((135.0, 192, 2, 1, 84), (115.0, 288, 3, 1, 72)),
    180.0: ((130.0, 288, 2, 1, 84), (120.0, 240, 3, 1, 80)),
    185.0: ((135.0, 288, 2, 1, 84), (125.0, 240, 3, 1, 80)),
    190.0: ((145.0, 240, 2, 1, 88), (140.0, 192, 3, 1, 80)),
    195.0: ((150.0, 240, 2, 1, 88), (130.0, 288, 3, 1, 80)),
    200.0: ((150.0, 288, 2, 1, 88), (140.0, 240, 3, 1, 84)),
    205.0: ((145.0, 240, 3, 1, 64), (145.0, 240, 3, 1, 84)),
    210.0: ((150.0, 240, 3, 1, 64), (150.0, 240, 3, 1, 84)),
    215.0: ((165.0, 288, 2, 1, 92), (150.0, 288, 3, 1, 76)),
    220.0: ((170.0, 192, 3, 1, 68), (155.0, 288, 3, 1, 72)),
    225.0: ((190.0, 240, 2, 1, 84), (165.0, 240, 3, 1, 84)),
    230.0: ((195.0, 240, 2, 1, 88), (170.0, 240, 3, 1, 92)),
    235.0: ((200.0, 240, 2, 1, 92), (170.0, 288, 3, 1, 76)),
    240.0: ((205.0, 240, 2, 1, 96), (170.0, 288, 3, 1, 92)),
    245.0: ((215.0, 240, 2, 1, 88), (190.0, 240, 3, 1, 84)),
    250.0: ((215.0, 240, 2, 1, 104), (195.0, 240, 3, 1, 88)),
    255.0: ((240.0, 192, 2, 1, 104), (200.0, 240, 3, 1, 92)),
}
PLAN16 = {}
for _b0, (_pf, _pg) in _PLAN.items():
    PLAN16[float(_b0), "f"] = _pf
    PLAN16[float(_b0), "g"] = _pg
BOOT16 = (230, 288, 12, 2)
