"""stnfs descent model (DESCENT.md): ball-region special-q nodes with two-tier
windows, semismooth lower sums, survivor upper bounds at each bin's lowest norms,
missing factor-base logs (discard or descend), negative-binomial node failure and
the finite-tree work recurrence; the boot by ECM trials on re-randomised certified
p^6 lattices.  Stdlib only; the shipped PLAN (offline grid search over n, nA, kA,
kB, R and the missing-log option per bin, node failure <= 0.05) is evaluated and
checked."""
import math, functools
from decimal import Decimal, localcontext

LN2 = math.log(2); UNIT = 20000.0
KMAX = 64
RSTEP = 1.0 / 128


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
        # series truncated at nt terms (exact to 1e-7 relative through u = 64 at
        # nt = 600; nt = 170 was not: it plateaued beyond u ~ 32); trailing terms
        # below 1e-25 of the largest are dropped
        out = [None]
        for c in coefs[1:]:
            big = max(abs(x) for x in c)
            while len(c) > 1 and abs(c[-1]) < big * Decimal('1e-25'):
                c = c[:-1]
            out.append([float(x) for x in c])
        return out


RHO_COEFS = _rho_coefficients()


def rho_series(u):
    if u <= 1:
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
    """lower bracket of rho(u): the table value at the next grid point."""
    if u <= 1:
        return 1.0
    i = math.ceil(u / RSTEP)
    return RHO_TAB[i] if i < len(RHO_TAB) else 0.0


def rho_hi(u):
    """upper bracket: the previous grid point; beyond the table 1/floor(u)!."""
    if u <= 1:
        return 1.0
    i = math.floor(u / RSTEP)
    if i + 1 >= len(RHO_TAB):
        return 1.0 / math.factorial(math.floor(u))
    return RHO_TAB[i]


def kappa(u, t):
    return rho_series(u - t) / rho_series(u) if rho_series(u) > 0 else 0.0


@functools.lru_cache(None)
def miss_side(u4, log2_rels, u_main):
    """missing-log ideals on one m-bit side, u = m/69 = u4/4 (DESCENT.md): sum over
    1/10-bit bands b < 69 of kappa(u, b/69) ln((b+.1)/b) P(Poisson(lambda) <= 2),
    lambda = 2^log2_rels 2 kappa(u_main, b/69)/2^b."""
    u = u4 / 4.0; tot = 0.0
    for b10 in range(400, 690):
        b = b10 / 10.0
        lam = 2.0 ** (log2_rels - b) * 2 * kappa(u_main, b / 69.0)
        tot += kappa(u, b / 69.0) * math.log((b + 0.1) / b) * math.exp(-lam) * (1 + lam + lam * lam / 2)
    return tot


def missing(Pm, nf, ng):
    """missing-log ideals per relation: f side at nf bits, g side at ng bits."""
    f = lambda m, um: miss_side(math.ceil(4 * m / 69.0), Pm.LOG2_RELS, um)
    return (f(nf, Pm.MAIN_UF) if nf else 0.0) + (f(ng, Pm.MAIN_UG) if ng else 0.0)


# ---------------- two-window semismooth lower sums ----------------
KA_MAX = 16


def _conv(a, b):
    out = [0.0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


@functools.lru_cache(None)
def window_rows(n, base, iA, mesh):
    """rows[(jA, jB)]: measure/(jA! jB!) of jA A-type and jB B-type interval
    products, indexed by the sum of interval indices; the grid is h = (n-base)/mesh,
    A-type intervals are those with index < iA (so nA = base + iA h), jA <= KA_MAX,
    jB <= 1.  Interval weight = log(right/left) = int dt/t."""
    h = (n - base) / mesh
    wA = [math.log1p(h / (base + i * h)) if i < iA else 0.0 for i in range(mesh)]
    wB = [math.log1p(h / (base + i * h)) if i >= iA else 0.0 for i in range(mesh)]
    PA = [[1.0]]
    for j in range(1, KA_MAX + 1):
        PA.append([x / j for x in _conv(PA[-1], wA)])
    PB = [[1.0], wB]
    rows = {}
    for jA in range(KA_MAX + 1):
        for jB in range(2):
            rows[(jA, jB)] = PA[jA] if jB == 0 else _conv(PA[jA], PB[1])
    return rows, h


@functools.lru_cache(None)
def ss_cum2(m, base, n, iA, mesh, alpha, trunc):
    """lower sums L[jA][jB] (cumulative in jA, exact in jB) at integer m."""
    rows, h = window_rows(n, base, iA, mesh)
    out = [[0.0, 0.0] for _ in range(KA_MAX + 1)]
    for (jA, jB), row in rows.items():
        j = jA + jB
        imax = len(row) - 1
        # feasible cells have left + j h <= m; keep those within trunc*base of the
        # largest feasible left endpoint sum (dropping cells keeps a lower sum)
        ifeas = min(imax, int(math.floor((m - j * h - j * base) / h))) if h > 0 else imax
        if ifeas < 0:
            continue
        imin = 0 if not trunc else max(0, ifeas - int(trunc * base / h) - 1)
        s = 0.0
        tab, ntab, ceil = RHO_TAB, len(RHO_TAB), math.ceil
        for i in range(imin, ifeas + 1):
            weight = row[i]
            if not weight:
                continue
            left = j * base + i * h
            if left + j * h <= m:
                u = (m - left + alpha) / base   # rho_lo(u), inlined
                if u <= 1:
                    s += weight
                else:
                    k = ceil(u / RSTEP)
                    if k < ntab:
                        s += weight * tab[k]
        out[jA][jB] = s
    cum = [[0.0, 0.0] for _ in range(KA_MAX + 1)]
    for jB in range(2):
        acc = 0.0
        for jA in range(KA_MAX + 1):
            acc += out[jA][jB]
            cum[jA][jB] = acc
    return tuple(tuple(x) for x in cum)


def ss_lower2(m, base, n, iA, kA, kB, mesh, alpha, trunc=3.0):
    """P(all primes <= 2^n, <= kA in (2^base, 2^nA], <= kB in (2^nA, 2^n]) lower sum."""
    m = math.ceil(m)
    if n <= base:
        return rho_lo((m + alpha) / base)
    cum = ss_cum2(m, base, n, iA, mesh, alpha, trunc)
    return cum[kA][0] + (cum[kA][1] if kB else 0.0)


def surv_upper(m, base, T):
    """P(product of the prime factors above 2^base of an m-bit integer <= 2^T)
    <= (T/base) rho((m-T)/base) for base <= T < m."""
    if T <= 0:
        return 0.0
    if T >= m:
        return 1.0
    return min(1.0, max(1.0, T / base) * rho_hi((math.floor(m) - T) / base))


# ---------------- ECM price (GMP-ECM table; measured ecm12 Montgomery multiplication) ----
ECM_TAB = {20: (11e3, 86), 25: (5e4, 214), 30: (25e4, 430), 35: (1e6, 904), 40: (3e6, 2350), 45: (11e6, 4480),
           50: (43e6, 7553), 55: (11e7, 17769), 60: (26e7, 42017), 65: (85e7, 69408), 70: (29e8, 102212),
           75: (76e8, 188056), 80: (25e9, 265557)}
MONT5 = 367.76


def ecm_instr(fbits, nbits, success):
    """instructions of an ECM schedule finding every prime of <= fbits bits in an
    nbits-bit integer with probability >= success per prime: GMP-ECM (B1,
    expected curves) at the next table digit size (geometric extrapolation
    beyond 80 digits), curves x ln(1/(1-success)), stage 1 = 1.44 B1 ladder bits
    at 11 modular multiplications, stage 2 charged as a second stage 1, one
    L-limb Montgomery multiplication = 367.76 (L/5)^2 instructions."""
    d = fbits / 3.3219
    ks = sorted(ECM_TAB)
    if d <= ks[-1]:
        k = min(x for x in ks if x >= d); b1, cur = ECM_TAB[k]
    else:
        b1, cur = ECM_TAB[ks[-1]]; e = (d - ks[-1]) / 5; b1 *= 3.3 ** e; cur *= 1.5 ** e
    curves = cur * math.log(1 / (1 - success))
    L = max(5, math.ceil(nbits / 64)); mulmod = MONT5 * (L / 5.0) ** 2
    return curves * 2 * 1.44 * b1 * 11 * mulmod


def log2li(l):
    x = l * LN2
    return l - math.log2(x) + math.log2(1 + 1 / x + 2 / x / x)


# exact half counts (N_12(R^2)-1)/2 of the ball regions (dsieve16 count)
CELLS = {4: 13315520, 6: 1573359464, 8: 47982017060, 10: 687193400292, 12: 6073556677752, 14: 38414493101628,
         16: 190058478539564, 20: 2754504316057176, 24: 24505531366261200, 28: 155610653675916932,
         32: 771909602591773000, 40: 11221331971811657968, 48: 99995402971660463340,
         56: 635616373364584622852, 64: 3155035383949694198996, 80: 45900097950907999286464,
         96: 409193506988359583551012, 112: 2601677688965738307448548, 128: 12916156661164168796080328, 160: 187942845067274337246511324, 192: 1675659512501707896389090868}


class P:
    """the parameter set (proof.md / DESCENT.md)."""
    LB = 69.0; LX = 69.0; PMIN = 8.0; SLACK = 96.0
    ALPHA_F = 5.0; ALPHA_G = 1.0; LOG2U6 = 6 * 63.714
    # cf sextile maxima of the pooled ball norms, each paired with the overall cg
    # maximum; both shifted by the worst sampled per-(ideal, radius) mean minus the
    # pooled mean (cf 20.68 + 1.58 -> 23, cg -2.53 + 7.81 -> 6), stats16.json
    QBINS = [(1 / 6, -26.4 + 23, 14.5 + 6), (1 / 6, -11.2 + 23, 14.5 + 6), (1 / 6, -0.1 + 23, 14.5 + 6),
             (1 / 6, 10.5 + 23, 14.5 + 6), (1 / 6, 23.2 + 23, 14.5 + 6), (1 / 6, 86.6 + 23, 14.5 + 6)]
    BOOT_BINS = [(0.5, 1.4 + 4), (0.4, 6.8 + 4), (0.1, 12.7 + 4)]
    # lower norm edges per sextile for the survivor upper bound: previous sextile's
    # raw cf maximum (observed minimum -140.05 for the first), the observed cg
    # minimum -54.56 throughout (stats16.json), both shifted down by the margin
    QBINS_LO = [(-140.05 - 23, -54.56 - 6), (-26.4 - 23, -54.56 - 6), (-11.2 - 23, -54.56 - 6),
                (-0.1 - 23, -54.56 - 6), (10.5 - 23, -54.56 - 6), (23.2 - 23, -54.56 - 6)]
    # missing factor-base logarithms: main-phase relations 2^LOG2_RELS with mean
    # sizes u_f, u_g; an ideal with <= 2 occurrences is treated as having no log
    LOG2_RELS = 64.162; MAIN_UF = 750.0 / 69; MAIN_UG = 498.5 / 69   # relation-weighted means, rounded down
    THIN_MAX = 75               # nodes with q <= 2^80 discard relations touching a missing log
    UDUP = 0.25                 # distinct relations per raw relation (unit duplicates; see DESCENT.md)
    MIN_RELS = 16.0; PSI_MAX = 0.8; PHI_FLOOR = 0.02
    NB_R1 = 0.5; NB_R2 = 1.0    # negative-binomial shapes: windows of <= 1 / >= 2 primes (toy fits 0.60 / 1.86-12.2)
    BOUND_CELL = 400.0          # certified 64-cell block norm bounds, both sides, per cell
    HITS_PER_CELL = 2 * math.log(69.0 / 8.0)
    SIEVE_HIT = 1050.0; SCAN = 10.0; ZERO = 2.0   # separate f and g score arrays
    PAIR_SETUP = 2.0 ** 25.1; LLL_NODE = 2.0 ** 36; NORMS = 2.0 ** 26; TDIV = 2.0 ** 24.5; IDEAL = 2.0 ** 34; PRIME = 2.0 ** 24
    ECM_SUCCESS = 0.999
    KS = (1, 2, 3, 4, 6, 8, 12, 16)
    RADII = (4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192)
    MESH = 48; TRUNC = 0.0
    LOG2P6 = 6 * 380.99700
    LMAX = 310; BIN = 5; TSTEP = 10


def pairs_count(Pm):
    return 2.0 * 2 ** log2li(Pm.LX)


def node_cost(Pm, l_lo, l_hi, side, n, iA, kA, kB, R, thin=None):
    """own cost (instructions) and statistics of one node plan; None if infeasible."""
    pairs = pairs_count(Pm)
    cfmax = max(b[1] for b in Pm.QBINS); cgmax = max(b[2] for b in Pm.QBINS)
    S = math.log2(CELLS[R]); c = (l_hi + S) / 12
    h = (n - Pm.LX) / Pm.MESH if n > Pm.LX else 0.0
    nA = Pm.LX + iA * h if n > Pm.LX else Pm.LX
    T = kA * nA + kB * n + Pm.SLACK
    c_lo = (l_lo + S) / 12
    pr = 0.0; sv = 0.0
    for (w, cf, cg), (cf_lo, cg_lo) in zip(Pm.QBINS, Pm.QBINS_LO):
        nf = 72 * c + cf - (l_lo if side == 'f' else 0)
        ng = 12 * c + Pm.LOG2U6 + cg - (l_lo if side == 'g' else 0)
        pf = ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC)
        pg = ss_lower2(ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
        pr += w * pf * pg
        # survivors: upper bound at the bin's lowest norms (surv_upper decreases in m)
        nf_lo = 72 * c_lo + cf_lo - (l_hi if side == 'f' else 0)
        ng_lo = 12 * c_lo + Pm.LOG2U6 + cg_lo - (l_hi if side == 'g' else 0)
        sv += w * surv_upper(nf_lo, Pm.LX, T) * surv_upper(ng_lo, Pm.LX, T)
    if pr <= 0:
        return None
    nf_max = 72 * c + cfmax - (l_lo if side == 'f' else 0)
    ng_max = 12 * c + Pm.LOG2U6 + cgmax - (l_lo if side == 'g' else 0)
    det = 1.0 - (1 - Pm.ECM_SUCCESS) * ((nf_max + ng_max) / 20.0 + 2 * (kA + kB))
    miss = missing(Pm, nf_max, ng_max)
    if thin is None or l_lo <= Pm.THIN_MAX:
        thin = thin or l_lo <= Pm.THIN_MAX
    mu = CELLS[R] * pr * Pm.UDUP * det * (math.exp(-miss) if thin else 1.0)
    if mu < Pm.MIN_RELS:
        return None
    surv_cost = (Pm.NORMS + Pm.TDIV + ecm_instr(Pm.LX, nf_max, Pm.ECM_SUCCESS) + ecm_instr(Pm.LX, ng_max, Pm.ECM_SUCCESS)
                 + 2 * (ecm_instr(nA, T, Pm.ECM_SUCCESS) + Pm.PRIME) + Pm.IDEAL)
    survivors = CELLS[R] * sv + CELLS[R] * pr
    sieve = CELLS[R] * (Pm.HITS_PER_CELL * Pm.SIEVE_HIT + Pm.SCAN + Pm.ZERO + Pm.BOUND_CELL) + pairs * Pm.PAIR_SETUP
    trials = CELLS[R] * surv_cost
    if trials < sieve + survivors * surv_cost:
        own = Pm.LLL_NODE + trials; how = 'ecm'; survivors = CELLS[R]
    else:
        own = Pm.LLL_NODE + sieve + survivors * surv_cost; how = 'sieve'
    return dict(own=own, n=n, iA=iA, nA=nA, kA=kA, kB=kB, R=R, S=S, mu=mu, how=how, miss=miss,
                miss_children=0.0 if thin else miss, thin=bool(thin), pr_log2=math.log2(pr),
                surv_log2=math.log2(survivors), own_log2=math.log2(own), nf=nf_max, ng=ng_max, det=det, T=T,
                sieve_log2=math.log2(sieve), survcost_log2=math.log2(survivors * surv_cost))


def bins_of(Pm):
    bins = []
    l = int(Pm.LB) + 1
    while l < Pm.LMAX:
        bins.append((l, l + Pm.BIN)); l += Pm.BIN
    return bins


def phi_nb(mean, Pm, window):
    """P(no distinct relation) for a negative-binomial count of the given mean."""
    r = Pm.NB_R1 if window <= 1 else Pm.NB_R2
    return max((1.0 + mean / r) ** (-r), Pm.PHI_FLOOR)


def evaluate(Pm, PLAN, BOOT):
    """evaluate the shipped plan: per-bin W (instructions), phi, details, and the
    two-log expected work; raises on an infeasible plan entry."""
    bins = bins_of(Pm)
    W = {}; PHI = {}; DET = {}

    def wmax(side, n):
        return max([W[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])

    def phimax(side, n):
        return max([PHI[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])

    def miss_terms():
        """failure probability and work of a missing-log child: the worst node with
        q <= 2^(THIN_MAX + 5) (those nodes discard relations touching missing logs)."""
        low = [k for k in W if k[0][0] <= Pm.THIN_MAX]
        return (max([PHI[k] for k in low] + [0.0]), max([W[k] for k in low] + [0.0]))
    for b in bins:
        for side in ('f', 'g'):
            n, iA, kA, kB, R, thin = PLAN[(b[0], side)]
            d = node_cost(Pm, b[0], b[1], side, float(n), iA, kA, kB, R, bool(thin))
            if d is None:
                raise ValueError('infeasible plan entry %s %s' % (b, side))
            nA = d['nA']
            pm, wm = miss_terms()
            psi = kA * (phimax('f', nA) + phimax('g', nA)) + kB * (phimax('f', n) + phimax('g', n)) + d['miss_children'] * pm
            if psi > Pm.PSI_MAX:
                raise ValueError('psi too large at %s %s' % (b, side))
            retry = 1.0 / (1.0 - psi)
            sub = kA * (wmax('f', nA) + wmax('g', nA)) + kB * (wmax('f', n) + wmax('g', n)) + d['miss_children'] * wm
            W[(b, side)] = d['own'] + retry * sub
            PHI[(b, side)] = phi_nb(d['mu'] * (1 - psi), Pm, kA + kB)
            DET[(b, side)] = dict(d, psi=psi, retry=retry, sub_log2=math.log2(sub) if sub > 0 else None,
                                  W_log2=math.log2(W[(b, side)]), phi=PHI[(b, side)])
    n1, iA, kA, Rb = BOOT
    Sb = math.log2(CELLS[Rb]); cb = (Pm.LOG2P6 + Sb) / 12
    h = (n1 - Pm.LX) / Pm.MESH; nA = Pm.LX + iA * h
    m_max = 12 * cb + Pm.LOG2U6 + max(x[1] for x in Pm.BOOT_BINS)
    pr = sum(w * ss_lower2(12 * cb + Pm.LOG2U6 + cg, Pm.LX, float(n1), iA, kA, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
             for w, cg in Pm.BOOT_BINS)
    det = 1.0 - (1 - Pm.ECM_SUCCESS) * (m_max / 20.0 + kA + 1)
    # distinct relations (unit mates collapse, UDUP); when a ball is exhausted the
    # boot draws a fresh exponent e and a new certified lattice (LLL_NODE each)
    trials = Pm.MIN_RELS / (pr * det * Pm.UDUP)
    balls = 1.0 + trials / CELLS[Rb]
    trial = Pm.NORMS + ecm_instr(nA, m_max, Pm.ECM_SUCCESS) + Pm.PRIME + Pm.IDEAL
    pm, wm = miss_terms(); mb = missing(Pm, 0, m_max)
    psi = kA * phimax('g', nA) + phimax('g', float(n1)) + mb * pm
    if psi > Pm.PSI_MAX:
        raise ValueError('boot psi too large')
    retry = 1 / (1 - psi)
    sub = kA * wmax('g', nA) + wmax('g', float(n1)) + mb * wm
    one_log = Pm.LLL_NODE * balls + trials * trial + retry * sub
    phi_boot = psi ** Pm.MIN_RELS   # all 16 kept relations fail (each fails with probability <= psi)
    boot = dict(n1=n1, nA=nA, kA=kA, R=Rb, S=Sb, m=m_max, pr_log2=math.log2(pr), trials_log2=math.log2(trials),
                trial_log2=math.log2(trial), balls=balls, miss=mb, own_log2=math.log2(trials * trial), sub_log2=math.log2(sub),
                psi=psi, retry=retry, phi=phi_boot, one_log_log2=math.log2(one_log))
    return dict(W=W, PHI=PHI, DET=DET, bins=bins, boot=boot, two_log_units=math.log2(2 * one_log / UNIT),
                one_log_instr=one_log)

# the shipped plan (one entry per (bin lower edge, side):
# (n, iA, kA, kB, R), nA = 69 + iA (n - 69)/48)
PLAN16 = {
    (70, 'f'): (69, 48, 2, 0, 192, 1),
    (70, 'g'): (69, 48, 2, 0, 192, 1),
    (75, 'f'): (75, 48, 2, 0, 96, 1),
    (75, 'g'): (75, 48, 2, 0, 112, 1),
    (80, 'f'): (75, 48, 2, 0, 128, 1),
    (80, 'g'): (75, 48, 2, 0, 160, 1),
    (85, 'f'): (75, 48, 2, 0, 192, 1),
    (85, 'g'): (75, 48, 2, 0, 192, 0),
    (90, 'f'): (85, 18, 1, 1, 160, 1),
    (90, 'g'): (85, 18, 1, 1, 192, 1),
    (95, 'f'): (85, 18, 1, 1, 192, 0),
    (95, 'g'): (85, 18, 2, 1, 192, 1),
    (100, 'f'): (85, 18, 2, 1, 192, 0),
    (100, 'g'): (85, 48, 3, 0, 192, 1),
    (105, 'f'): (85, 48, 3, 0, 192, 1),
    (105, 'g'): (85, 48, 4, 0, 192, 1),
    (110, 'f'): (95, 30, 2, 1, 160, 1),
    (110, 'g'): (95, 30, 3, 1, 160, 1),
    (115, 'f'): (105, 21, 2, 1, 192, 1),
    (115, 'g'): (95, 30, 3, 1, 192, 0),
    (120, 'f'): (95, 30, 3, 1, 192, 1),
    (120, 'g'): (95, 48, 4, 0, 192, 0),
    (125, 'f'): (105, 21, 3, 1, 192, 0),
    (125, 'g'): (105, 35, 3, 1, 192, 0),
    (130, 'f'): (105, 35, 3, 1, 160, 1),
    (130, 'g'): (105, 35, 4, 1, 192, 1),
    (135, 'f'): (105, 35, 3, 1, 192, 0),
    (135, 'g'): (115, 27, 4, 1, 192, 1),
    (140, 'f'): (115, 27, 3, 1, 192, 0),
    (140, 'g'): (125, 22, 4, 1, 192, 0),
    (145, 'f'): (115, 27, 4, 1, 160, 0),
    (145, 'g'): (115, 38, 4, 1, 192, 0),
    (150, 'f'): (125, 22, 4, 1, 192, 1),
    (150, 'g'): (135, 26, 4, 1, 192, 1),
    (155, 'f'): (135, 26, 3, 1, 192, 0),
    (155, 'g'): (135, 33, 4, 1, 192, 1),
    (160, 'f'): (125, 31, 4, 1, 160, 0),
    (160, 'g'): (135, 41, 4, 1, 160, 1),
    (165, 'f'): (135, 26, 4, 1, 192, 0),
    (165, 'g'): (145, 35, 4, 1, 192, 1),
    (170, 'f'): (125, 39, 4, 1, 192, 0),
    (170, 'g'): (155, 31, 4, 1, 192, 0),
    (175, 'f'): (135, 33, 4, 1, 192, 0),
    (175, 'g'): (135, 33, 6, 1, 192, 0),
    (180, 'f'): (145, 35, 4, 1, 128, 0),
    (180, 'g'): (155, 42, 4, 1, 192, 1),
    (185, 'f'): (145, 35, 4, 1, 192, 1),
    (185, 'g'): (165, 38, 4, 1, 192, 0),
    (190, 'f'): (155, 31, 4, 1, 192, 0),
    (190, 'g'): (165, 43, 4, 1, 192, 0),
    (195, 'f'): (155, 37, 4, 1, 160, 0),
    (195, 'g'): (175, 43, 4, 1, 192, 1),
    (200, 'f'): (155, 42, 4, 1, 160, 0),
    (200, 'g'): (165, 38, 6, 1, 96, 0),
    (205, 'f'): (165, 38, 4, 1, 160, 0),
    (205, 'g'): (185, 31, 6, 1, 112, 0),
    (210, 'f'): (165, 38, 4, 1, 192, 0),
    (210, 'g'): (195, 40, 4, 1, 192, 0),
    (215, 'f'): (185, 31, 4, 1, 192, 0),
    (215, 'g'): (185, 31, 6, 1, 192, 0),
    (220, 'f'): (175, 39, 4, 1, 192, 0),
    (220, 'g'): (205, 34, 6, 1, 64, 0),
    (225, 'f'): (185, 40, 4, 1, 160, 1),
    (225, 'g'): (185, 48, 8, 0, 32, 1),
    (230, 'f'): (185, 40, 4, 1, 192, 0),
    (230, 'g'): (205, 34, 6, 1, 112, 0),
    (235, 'f'): (185, 40, 6, 1, 40, 0),
    (235, 'g'): (195, 48, 8, 0, 28, 1),
    (240, 'f'): (195, 40, 6, 1, 32, 0),
    (240, 'g'): (195, 48, 8, 0, 40, 1),
    (245, 'f'): (195, 40, 6, 1, 40, 0),
    (245, 'g'): (205, 48, 8, 0, 28, 1),
    (250, 'f'): (195, 44, 6, 1, 32, 0),
    (250, 'g'): (205, 48, 8, 0, 32, 1),
    (255, 'f'): (185, 48, 8, 0, 32, 1),
    (255, 'g'): (215, 48, 8, 0, 24, 1),
    (260, 'f'): (205, 44, 6, 1, 32, 1),
    (260, 'g'): (215, 48, 8, 0, 32, 1),
    (265, 'f'): (195, 48, 8, 0, 28, 1),
    (265, 'g'): (215, 48, 8, 0, 40, 1),
    (270, 'f'): (195, 48, 8, 0, 32, 1),
    (270, 'g'): (225, 48, 8, 0, 28, 1),
    (275, 'f'): (195, 48, 8, 0, 40, 1),
    (275, 'g'): (225, 48, 8, 0, 32, 1),
    (280, 'f'): (205, 48, 8, 0, 28, 1),
    (280, 'g'): (235, 48, 8, 0, 28, 1),
    (285, 'f'): (205, 48, 8, 0, 32, 1),
    (285, 'g'): (235, 48, 8, 0, 32, 1),
    (290, 'f'): (215, 48, 8, 0, 24, 1),
    (290, 'g'): (245, 48, 8, 0, 24, 1),
    (295, 'f'): (215, 48, 8, 0, 28, 1),
    (295, 'g'): (245, 48, 8, 0, 28, 1),
    (300, 'f'): (215, 48, 8, 0, 32, 1),
    (300, 'g'): (255, 48, 8, 0, 24, 1),
    (305, 'f'): (215, 48, 8, 0, 40, 1),
    (305, 'g'): (255, 48, 8, 0, 28, 1),
}
# boot plan: (n1, iA, kA, R)
BOOT16 = (235, 42, 12, 4)


def report(res):
    """compact per-bin table for the experiment report."""
    rows = []
    for b in res['bins']:
        for side in ('f', 'g'):
            d = res['DET'][(b, side)]
            rows.append(dict(bin=b[0], side=side, n=int(d['n']), nA=round(d['nA'], 1), kA=d['kA'], kB=d['kB'], R=d['R'],
                             S=round(d['S'], 2), how=d['how'], mu=round(d['mu'], 2), pr_log2=round(d['pr_log2'], 2),
                             surv_log2=round(d['surv_log2'], 2), own_log2=round(d['own_log2'], 3), W_log2=round(d['W_log2'], 3),
                             psi=round(d['psi'], 4), phi=round(d['phi'], 6), det=round(d['det'], 4),
                             nf=round(d['nf'], 1), ng=round(d['ng'], 1), T=round(d['T'], 1)))
    return rows


if __name__ == '__main__':
    import time, sys, json
    Pm = P()
    PLAN, BOOT = PLAN16, BOOT16
    t = time.time()
    res = evaluate(Pm, PLAN, BOOT)
    print(f"two-log 2^{res['two_log_units']:.3f} units; boot {res['boot']}; evaluate {time.time()-t:.2f}s")
    for r in report(res):
        print(r)
