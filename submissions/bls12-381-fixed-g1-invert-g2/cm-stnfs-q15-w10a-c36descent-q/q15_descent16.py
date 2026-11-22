"""q15_descent16.py — Exact 16D Q(zeta_15) port of c36p/descent16.py.

Uncapped 16D special-q lattice sieving (zero per-prime caps, zero survivor caps, zero Q*q' fallback lattices, 100% eager cofactorization).
EXACT c36p architecture:
  - Uncapped dsieve16 + BOUND_CELL live yield certificate + q15_retain16.retention
  - NegBin r = 0.5 (<= 1 effective window prime) / r = 1.0 (>= 2), DISC = 0.001, PHI_FLOOR = 0.02
  - Lazy cofactorization (lazy_fraction exact NegBin sum)
  - Whole-ball renewal boot with BOOT_RELS = 18, BALL_PAIR = 2^24, phi_boot = psi^18
  - Exact 16D norms: c = (l + S)/16, nf = 72*c + cf, ng = 24*c + 8*log2(|u|) + cg
  - Boot 16D lattice (det p^4): cb = (4*log2(p) + S_b)/16, ng_boot = 24*cb + 8*log2(|u|) + cg
  - Missing FB ideals under C_12 orbit (ORBIT = 12), MAIN_LB = 71.5, LOG2_RELS = 63.36
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


def rho_lo(u):
    if u <= 1.0:
        return 1.0
    i = math.ceil(u / RSTEP)
    return RHO_TAB[i] if i < len(RHO_TAB) else 0.0


def rho_hi(u):
    if u <= 1.0:
        return 1.0
    i = floor(u / RSTEP)
    if i + 1 >= len(RHO_TAB):
        return 1.0 / factorial(floor(u))
    return RHO_TAB[i]


def kappa(u, t):
    ru = rho_series(u)
    return rho_series(u - t) / ru if ru > 0 else 0.0


@functools.lru_cache(None)
def miss_side(u4, log2_rels, u_main, lx=71.0, main_lb=71.5, orbit=12):
    u = u4 / 4.0
    tot = 0.0
    hi10 = int(round(lx * 10.0))
    for b10 in range(400, hi10):
        b = b10 / 10.0
        lam = (2.0 ** (log2_rels - b)) * orbit * kappa(u_main, b / main_lb)
        lp = 0.75 * lam
        p_thin = exp(-lp) * (1.0 + lp) if lp < 700.0 else 0.0
        p_p3 = exp(-lam) * (1.0 + lam + 0.5 * lam * lam + (lam * lam * lam) / 6.0) if lam < 700.0 else 0.0
        p_miss = max(p_p3, p_thin)
        tot += kappa(u, b / lx) * log((b + 0.1) / b) * p_miss
    return tot


def missing(Pm, nf, ng):
    f = lambda m, um: miss_side(
        math.ceil(4.0 * m / Pm.LX), Pm.LOG2_RELS, um, Pm.LX, Pm.MAIN_LB, Pm.ORBIT
    )
    return (f(nf, Pm.MAIN_UF) if nf else 0.0) + (f(ng, Pm.MAIN_UG) if ng else 0.0)


KA_MAX = 20


def _conv(a, b):
    out = [0.0] * (len(a) + len(b) - 1)
    nz = [j for j, y in enumerate(b) if y]
    if nz:
        j0, j1 = nz[0], nz[-1] + 1
        bb = b[j0:j1]
        for i, x in enumerate(a):
            if x:
                out[i + j0 : i + j1] = [o + x * y for o, y in zip(out[i + j0 : i + j1], bb)]
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
    tab, ntab, ceil = RHO_TAB, len(RHO_TAB), math.ceil
    sp = math.sumprod
    inv_br = 1.0 / (base * RSTEP)
    hk = h * inv_br
    for (jA, jB), (row, r_lo, r_hi) in rows.items():
        if r_lo >= r_hi:
            continue
        j = jA + jB
        ifeas = min(r_hi - 1, int(floor((m - j * h - j * base) / h))) if h > 0 else (r_hi - 1)
        if ifeas < r_lo:
            continue
        imin = r_lo if not trunc else max(r_lo, ifeas - int(trunc * base / h) - 1)
        k0 = (m - j * base + alpha) * inv_br
        # Skip indices where ceil(k0 - i * hk) >= ntab (tab[k] == 0.0)
        if hk > 0.0:
            i_ntab = int((k0 - (ntab - 1)) / hk)
            if i_ntab > imin:
                imin = i_ntab
        if imin > ifeas:
            continue
        # Fast C-level sumprod over non-zero slice [imin..ifeas]
        if k0 - ifeas * hk <= 0.0 or k0 - imin * hk >= ntab:
            s = 0.0
            for i in range(imin, ifeas + 1):
                k = ceil(k0 - i * hk)
                s += row[i] * (1.0 if k <= 0 else tab[k] if k < ntab else 0.0)
        else:
            s = sp(row[imin : ifeas + 1], [tab[ceil(k0 - i * hk)] for i in range(imin, ifeas + 1)])
        out[jA][jB] = s
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


def ecm_instr(fbits, nbits, success):
    import q15_ecm16
    return q15_ecm16.ecm_instr(fbits, nbits, success)


def log2li(l):
    x = l * LN2
    return l - log2(x) + log2(1.0 + 1.0 / x + 2.0 / (x * x))


# Exact 16D ball half-counts loaded/computed via q15_cells16
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
    LB = 71.0
    LX = 71.0
    PMIN = 8.0
    SLACK = 96.0
    ALPHA_F = 2.74
    ALPHA_G = 0.77
    LOG2U8 = 8.0 * log2(0xd201000000010000)  # 509.7141788...
    LOG2P4 = 4.0 * 380.99700  # det(L_boot) = p^4 in Z^16
    # Default placeholders (overwritten verbatim from q15_stats.json by load_stats()!)
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
    BOOT_BINS = [(0.5, -260.0), (0.4, -250.0), (0.1, -240.0)]
    BOOT_SUB = []
    PRCAP_NODE = 4.0
    PRCAP_BOOT = 1.5
    LOG2_RELS = 63.658
    MAIN_LB = 71.5
    MAIN_UF = 549.0 / MAIN_LB  # Q(zeta_15): N_f ~ 549..556 bits
    MAIN_UG = 678.0 / MAIN_LB  # Q(zeta_15): N_g ~ 678..681 bits
    ORBIT = 12  # Aut = C_12!
    THIN_MAX = 300
    UDUP = 0.65
    MIN_RELS = 18.0
    PSI_MAX = 0.8
    PHI_FLOOR = 0.02
    DISC = 0.001
    BOOT_RELS = 18
    BALL_PAIR = 2.0 ** 24
    NB_R1 = 0.5
    NB_R2 = 1.0
    BOUND_CELL = 3600.0
    HITS_PER_CELL = 2.0 * log(71.0 / 8.0)
    SIEVE_HIT = 2200.0
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
    LMAX = 256
    BIN = 5
    # Optional fallback retention when q15_retain16 is not yet present
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
    return 2.0 * (2.0 ** log2li(Pm.LX))


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
    det = 1.0 - (1.0 - Pm.ECM_SUCCESS) * ((nf_max + ng_max) / 20.0 + 2.0 * (kA + kB))
    miss = missing(Pm, nf_max, ng_max)
    if thin is None or l_lo <= Pm.THIN_MAX:
        thin = bool(thin) or (l_lo <= Pm.THIN_MAX)
    mu = C_R * pr * Pm.UDUP * det * (exp(-miss) if thin else 1.0)
    if mu < Pm.MIN_RELS:
        return None
    surv_cost = (
        Pm.NORMS
        + Pm.TDIV
        + ecm_instr(Pm.LX, nf_max, Pm.ECM_SUCCESS)
        + ecm_instr(Pm.LX, ng_max, Pm.ECM_SUCCESS)
        + 2.0 * (ecm_instr(nA, T, Pm.ECM_SUCCESS) + Pm.PRIME)
        + Pm.IDEAL
    )
    survivors = C_R * sv + C_R * pr
    hits_pc = 2.0 * log(Pm.LX / Pm.PMIN)
    sieve = C_R * (hits_pc * Pm.SIEVE_HIT + Pm.SCAN + Pm.ZERO + Pm.BOUND_CELL) + pairs * Pm.PAIR_SETUP
    trials = C_R * surv_cost
    if trials < sieve + survivors * surv_cost:
        own = Pm.LLL_NODE + trials
        how = "ecm"
        survivors = C_R
        ret, loss_f, loss_g = 1.0, 0.0, 0.0
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
    # Match c36p convention: first bin (LB, LB + BIN) or (69, 75) -> for LB=71, 5-bit bins 71..256!
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


def lazy_fraction(mu, r, psi, tol=1e-13, m_cap=20000):
    retry = (1.0 + 1e-12) / (1.0 - psi)
    lq = log(mu / (r + mu))
    lp0 = r * log(r / (r + mu))
    tot = mass = 0.0
    m = 0
    pm = 1.0
    while (m <= mu or 1.0 - mass > tol) and m < m_cap:
        pm = exp(lgamma(m + r) - lgamma(r) - lgamma(m + 1) + lp0 + m * lq)
        tot += pm * min(1.0, retry / (m + 1) + psi ** m)
        mass += pm
        m += 1
        if m > mu and pm < 1e-20:
            break
    rem = max(0.0, 1.0 - mass)
    if m >= m_cap and rem > 0.0:
        # Since min(1, retry/(m+1) + psi^m) is strictly decreasing in m, its value at m_cap upper-bounds every term m >= m_cap!
        tail_val = min(1.0, retry / (m_cap + 1) + (psi ** m_cap if m_cap < 500 else 0.0))
        return min(1.0, tot + rem * tail_val)
    return min(1.0, tot + rem + pm * (r + mu) / r)


def evaluate(Pm, PLAN, BOOT):
    cells_map = load_cells()
    if Pm.FALLBACK_RET is None:
        import q15_retain16
        pass  # q15_retain16.loss is a pure function of its explicit arguments
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
                frac = 1.0  # 100% eager cofactorization of all survivors (no E[s]*E[frac(M)] assumption!)
                d["own"] = Pm.LLL_NODE + (2.0 ** d["sieve_log2"]) + (2.0 ** d["survcost_log2"]) * frac
                d["lazy"] = frac
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
    det = 1.0 - (1.0 - Pm.ECM_SUCCESS) * (m_max / 20.0 + kA + 1.0)
    mb = missing(Pm, 0.0, m_max)
    per_rel = 1.0 / (pr * det * Pm.UDUP * _get_counted_fraction(Pm) * exp(-mb))
    trial = Pm.NORMS + ecm_instr(nA, m_max, Pm.ECM_SUCCESS) + Pm.PRIME + Pm.IDEAL
    psi = kA * phimax("g", nA) + phimax("g", float(n1))
    if psi > Pm.PSI_MAX:
        raise ValueError(f"boot psi too large: {psi:.4f}")
    retry = 1.0 / (1.0 - psi)
    balls = (1.0 - psi ** Pm.BOOT_RELS) / (1.0 - psi) * (per_rel / C_Rb + Pm.BALL_PAIR)
    trials = balls * C_Rb
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
    71: ((71, 32, 0, 0, 40, 1), (71, 32, 0, 0, 40, 1)),
    76: ((76, 0, 0, 1, 32, 1), (76, 0, 0, 1, 32, 1)),
    81: ((76, 0, 0, 1, 40, 1), (76, 0, 0, 1, 40, 1)),
    86: ((76, 32, 2, 0, 28, 1), (76, 32, 2, 0, 28, 1)),
    91: ((76, 32, 2, 0, 32, 1), (76, 32, 2, 0, 32, 1)),
    96: ((76, 32, 2, 0, 40, 1), (76, 32, 2, 0, 40, 1)),
    101: ((76, 32, 2, 0, 48, 1), (76, 32, 2, 0, 48, 1)),
    106: ((81, 12, 1, 1, 48, 1), (81, 16, 1, 1, 48, 1)),
    111: ((86, 12, 1, 1, 40, 1), (86, 12, 1, 1, 48, 1)),
    116: ((91, 24, 2, 0, 40, 1), (91, 24, 2, 0, 48, 1)),
    121: ((96, 12, 1, 1, 48, 1), (96, 16, 1, 1, 48, 1)),
    126: ((101, 16, 1, 1, 40, 1), (111, 20, 2, 0, 48, 1)),
    131: ((111, 12, 1, 1, 40, 1), (111, 12, 1, 1, 48, 1)),
    136: ((111, 16, 1, 1, 48, 1), (111, 20, 1, 1, 48, 1)),
    141: ((121, 24, 2, 0, 48, 1), (101, 16, 2, 1, 48, 1)),
    146: ((111, 12, 2, 1, 48, 1), (111, 12, 2, 1, 48, 1)),
    151: ((111, 16, 2, 1, 48, 1), (111, 12, 2, 1, 56, 1)),
    156: ((121, 12, 2, 1, 48, 1), (126, 12, 2, 1, 48, 1)),
    161: ((126, 12, 2, 1, 48, 1), (131, 12, 2, 1, 48, 1)),
    166: ((126, 12, 2, 1, 56, 1), (131, 16, 2, 1, 48, 1)),
    171: ((131, 12, 2, 1, 56, 1), (131, 16, 2, 1, 56, 1)),
    176: ((141, 12, 2, 1, 56, 1), (146, 12, 2, 1, 56, 1)),
    181: ((146, 12, 2, 1, 56, 1), (151, 16, 2, 1, 48, 1)),
    186: ((156, 12, 2, 1, 56, 1), (151, 16, 3, 1, 32, 1)),
    191: ((161, 12, 2, 1, 56, 1), (166, 16, 2, 1, 48, 1)),
    196: ((161, 12, 2, 1, 64, 1), (166, 16, 2, 1, 56, 1)),
    201: ((171, 12, 2, 1, 64, 1), (156, 20, 3, 1, 32, 1)),
    206: ((186, 12, 2, 1, 56, 1), (171, 20, 3, 1, 28, 1)),
    211: ((186, 12, 2, 1, 64, 1), (176, 20, 3, 1, 28, 1)),
    216: ((181, 16, 2, 1, 64, 1), (176, 20, 3, 1, 32, 1)),
    221: ((216, 24, 3, 1, 16, 1), (191, 24, 3, 1, 20, 1)),
    226: ((191, 16, 2, 1, 64, 1), (196, 20, 3, 1, 24, 1)),
    231: ((221, 24, 3, 1, 16, 1), (201, 20, 3, 1, 28, 1)),
    236: ((181, 24, 4, 1, 16, 1), (201, 20, 3, 1, 32, 1)),
    241: ((186, 24, 4, 1, 16, 1), (206, 24, 3, 1, 24, 1)),
    246: ((191, 24, 4, 1, 16, 1), (211, 24, 3, 1, 24, 1)),
    251: ((181, 28, 4, 1, 16, 1), (216, 24, 3, 1, 24, 1)),
    256: ((186, 28, 4, 1, 16, 1), (216, 24, 3, 1, 28, 1)),
}
PLAN16 = {(b, s): (v[i][0], 12 * v[i][1], *v[i][2:]) for b, v in _PLAN.items() for i, s in enumerate('fg')}
BOOT16 = (261, 336, 14, 2)
