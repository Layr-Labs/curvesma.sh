"""stnfs v16 descent model: ball-region special-q nodes with two-tier windows and
an explicit finite-tree success recurrence (DESCENT.md).

Node (special-q ideal of l bits on side f or g):
  region   = integer points of the c-space ball of radius R (one of each +-pair,
             exact counts CELLS[R]) in the LLL-reduced basis of L_q (det q);
  sieve    = every factor-base prime p in [2^8, 2^69] and every relative-degree-one
             root pair, by exact hit-lattice enumeration (dsieve16.c, measured);
  survivor = both residual masses (norm bound minus score) <= T_side = kA nA + kB n + SLACK;
  relation = after ECM at the factor-base level on both exact norms and ECM at
             level nA on the residuals, every residual prime is <= 2^n with at
             most kA primes in (2^69, 2^nA] and at most kB (0 or 1) prime in
             (2^nA, 2^n] per side (a kB=1 prime is identified by a primality test);
  children = those window primes (kA per side at bins <= nA, kB at bins <= n).
Boot: ECM trials at level nA on random points of the p^6 lattice's ball
(g-side norm only), the same two-tier predicate with n1 and kB = 1.

Yield (per cell): the standard semismooth heuristic with the two windows,
  P = sum_{jA<=kA, jB<=kB} (1/jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t
  (A-type t in [69, nA], B-type t in (nA, n]), as numerical LOWER sums on a mesh
  (rho decreasing: left endpoints, cells crossing sum t = m dropped, rows
  truncated), averaged over the measured (cf, cg) sextile bins of the ball norm
  offsets (QBINS; per-bin maxima, shifted by the worst sampled per-ideal mean),
  times UDUP (distinct relations per raw relation) and the ECM detection factor.
Survivors (per cell): the monotone UPPER bound (T/69) rho((m - T)/69) per side.
Success: negative-binomial(mu, r) distinct relations, independent relation subtrees,
  phi(node) = max((1 + mu (1 - psi)/r)^-r, PHI_FLOOR) (negative binomial, r = NB_R1 for windows of one prime, NB_R2 otherwise), psi = sum over children of
  phi(child bin maximum), expected relations tried <= 1/(1 - psi), children
  charged at the most expensive bin at or below their window.
Stdlib only; the shipped PLAN (chosen by an offline grid search over n, nA, kA, kB
and R per bin, each node's failure probability at most 0.05) is evaluated and checked.
"""
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
    QBINS = [(1 / 6, -26.4 + 13, -3.3 + 5), (1 / 6, -11.2 + 13, -0.5 + 5), (1 / 6, -0.1 + 13, 1.7 + 5),
             (1 / 6, 10.5 + 13, 3.4 + 5), (1 / 6, 23.2 + 13, 4.9 + 5), (1 / 6, 86.6 + 13, 14.5 + 5)]
    BOOT_BINS = [(0.5, 1.4 + 4), (0.4, 6.8 + 4), (0.1, 12.7 + 4)]
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


def node_cost(Pm, l_lo, l_hi, side, n, iA, kA, kB, R):
    """own cost (instructions) and statistics of one node plan; None if infeasible."""
    pairs = pairs_count(Pm)
    cfmax = max(b[1] for b in Pm.QBINS); cgmax = max(b[2] for b in Pm.QBINS)
    S = math.log2(CELLS[R]); c = (l_hi + S) / 12
    h = (n - Pm.LX) / Pm.MESH if n > Pm.LX else 0.0
    nA = Pm.LX + iA * h if n > Pm.LX else Pm.LX
    T = kA * nA + kB * n + Pm.SLACK
    pr = 0.0; sv = 0.0
    for w, cf, cg in Pm.QBINS:
        nf = 72 * c + cf - (l_lo if side == 'f' else 0)
        ng = 12 * c + Pm.LOG2U6 + cg - (l_lo if side == 'g' else 0)
        pf = ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC)
        pg = ss_lower2(ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
        pr += w * pf * pg
        sv += w * surv_upper(nf, Pm.LX, T) * surv_upper(ng, Pm.LX, T)
    if pr <= 0:
        return None
    nf_max = 72 * c + cfmax - (l_lo if side == 'f' else 0)
    ng_max = 12 * c + Pm.LOG2U6 + cgmax - (l_lo if side == 'g' else 0)
    det = 1.0 - (1 - Pm.ECM_SUCCESS) * ((nf_max + ng_max) / 20.0 + 2 * (kA + kB))
    mu = CELLS[R] * pr * Pm.UDUP * det
    if mu < Pm.MIN_RELS:
        return None
    surv_cost = (Pm.NORMS + Pm.TDIV + ecm_instr(Pm.LX, nf_max, Pm.ECM_SUCCESS) + ecm_instr(Pm.LX, ng_max, Pm.ECM_SUCCESS)
                 + 2 * (ecm_instr(nA, T, Pm.ECM_SUCCESS) + Pm.PRIME) + Pm.IDEAL)
    survivors = CELLS[R] * sv + mu / (Pm.UDUP * det)
    sieve = CELLS[R] * (Pm.HITS_PER_CELL * Pm.SIEVE_HIT + Pm.SCAN + Pm.ZERO + Pm.BOUND_CELL) + pairs * Pm.PAIR_SETUP
    trials = CELLS[R] * surv_cost
    if trials < sieve + survivors * surv_cost:
        own = Pm.LLL_NODE + trials; how = 'ecm'; survivors = CELLS[R]
    else:
        own = Pm.LLL_NODE + sieve + survivors * surv_cost; how = 'sieve'
    return dict(own=own, n=n, iA=iA, nA=nA, kA=kA, kB=kB, R=R, S=S, mu=mu, how=how, pr_log2=math.log2(pr),
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
    for b in bins:
        for side in ('f', 'g'):
            n, iA, kA, kB, R = PLAN[(b[0], side)]
            d = node_cost(Pm, b[0], b[1], side, float(n), iA, kA, kB, R)
            if d is None:
                raise ValueError('infeasible plan entry %s %s' % (b, side))
            nA = d['nA']
            psi = kA * (phimax('f', nA) + phimax('g', nA)) + kB * (phimax('f', n) + phimax('g', n))
            if psi > Pm.PSI_MAX:
                raise ValueError('psi too large at %s %s' % (b, side))
            retry = 1.0 / (1.0 - psi)
            sub = kA * (wmax('f', nA) + wmax('g', nA)) + kB * (wmax('f', n) + wmax('g', n))
            W[(b, side)] = d['own'] + retry * sub
            PHI[(b, side)] = phi_nb(d['mu'] * (1 - psi), Pm, kA + kB)
            DET[(b, side)] = dict(d, psi=psi, retry=retry, sub_log2=math.log2(sub) if sub > 0 else None,
                                  W_log2=math.log2(W[(b, side)]), phi=PHI[(b, side)])
    n1, iA, kA, Rb = BOOT
    Sb = math.log2(CELLS[Rb]); cb = (Pm.LOG2P6 + Sb) / 12
    h = (n1 - Pm.LX) / Pm.MESH; nA = Pm.LX + iA * h
    m_max = 12 * cb + Pm.LOG2U6 + max(x[1] for x in Pm.BOOT_BINS)
    pr = sum(w * ss_lower2(12 * cb + Pm.LOG2U6 + cg, Pm.LX, float(n1), iA, kA, 1, Pm.MESH, 0.0, Pm.TRUNC)
             for w, cg in Pm.BOOT_BINS)
    det = 1.0 - (1 - Pm.ECM_SUCCESS) * (m_max / 20.0 + kA + 1)
    # distinct relations (unit mates collapse, UDUP); when a ball is exhausted the
    # boot draws a fresh exponent e and a new certified lattice (LLL_NODE each)
    trials = Pm.MIN_RELS / (pr * det * Pm.UDUP)
    balls = 1.0 + trials / CELLS[Rb]
    trial = Pm.NORMS + ecm_instr(nA, m_max, Pm.ECM_SUCCESS) + Pm.PRIME + Pm.IDEAL
    psi = kA * phimax('g', nA) + phimax('g', float(n1))
    if psi > Pm.PSI_MAX:
        raise ValueError('boot psi too large')
    retry = 1 / (1 - psi)
    sub = kA * wmax('g', nA) + wmax('g', float(n1))
    one_log = Pm.LLL_NODE * balls + trials * trial + retry * sub
    phi_boot = psi ** Pm.MIN_RELS   # all 16 kept relations fail (each fails with probability <= psi)
    boot = dict(n1=n1, nA=nA, kA=kA, R=Rb, S=Sb, m=m_max, pr_log2=math.log2(pr), trials_log2=math.log2(trials),
                trial_log2=math.log2(trial), balls=balls, own_log2=math.log2(trials * trial), sub_log2=math.log2(sub),
                psi=psi, retry=retry, phi=phi_boot, one_log_log2=math.log2(one_log))
    return dict(W=W, PHI=PHI, DET=DET, bins=bins, boot=boot, two_log_units=math.log2(2 * one_log / UNIT),
                one_log_instr=one_log)

# the shipped plan (one entry per (bin lower edge, side):
# (n, iA, kA, kB, R), nA = 69 + iA (n - 69)/48)
PLAN16 = {
    (70, 'f'): (69, 48, 2, 0, 112),
    (70, 'g'): (69, 48, 2, 0, 128),
    (75, 'f'): (69, 48, 2, 0, 160),
    (75, 'g'): (75, 48, 2, 0, 80),
    (80, 'f'): (75, 48, 2, 0, 96),
    (80, 'g'): (75, 48, 2, 0, 112),
    (85, 'f'): (85, 18, 1, 1, 80),
    (85, 'g'): (75, 48, 3, 0, 128),
    (90, 'f'): (85, 18, 2, 1, 80),
    (90, 'g'): (85, 18, 1, 1, 128),
    (95, 'f'): (85, 18, 3, 1, 96),
    (95, 'g'): (85, 48, 2, 0, 160),
    (100, 'f'): (85, 48, 2, 0, 160),
    (100, 'g'): (85, 48, 3, 0, 128),
    (105, 'f'): (85, 48, 3, 0, 128),
    (105, 'g'): (95, 11, 4, 1, 160),
    (110, 'f'): (95, 30, 2, 1, 112),
    (110, 'g'): (95, 48, 3, 0, 128),
    (115, 'f'): (95, 30, 2, 1, 160),
    (115, 'g'): (95, 30, 3, 1, 160),
    (120, 'f'): (95, 48, 3, 0, 160),
    (120, 'g'): (95, 30, 4, 1, 160),
    (125, 'f'): (105, 35, 2, 1, 160),
    (125, 'g'): (115, 27, 2, 1, 192),
    (130, 'f'): (105, 21, 3, 1, 192),
    (130, 'g'): (115, 27, 3, 1, 160),
    (135, 'f'): (105, 35, 3, 1, 160),
    (135, 'g'): (125, 22, 4, 1, 112),
    (140, 'f'): (115, 27, 3, 1, 160),
    (140, 'g'): (125, 22, 4, 1, 160),
    (145, 'f'): (105, 48, 4, 0, 192),
    (145, 'g'): (115, 38, 4, 1, 160),
    (150, 'f'): (115, 38, 3, 1, 192),
    (150, 'g'): (145, 42, 4, 1, 40),
    (155, 'f'): (125, 31, 3, 1, 192),
    (155, 'g'): (135, 41, 3, 1, 192),
    (160, 'f'): (125, 39, 3, 1, 192),
    (160, 'g'): (145, 35, 4, 1, 96),
    (165, 'f'): (135, 33, 3, 1, 192),
    (165, 'g'): (145, 48, 6, 0, 40),
    (170, 'f'): (135, 26, 4, 1, 192),
    (170, 'g'): (155, 31, 4, 1, 160),
    (175, 'f'): (135, 33, 4, 1, 160),
    (175, 'g'): (175, 25, 4, 1, 160),
    (180, 'f'): (145, 42, 3, 1, 192),
    (180, 'g'): (135, 48, 6, 0, 160),
    (185, 'f'): (155, 37, 3, 1, 192),
    (185, 'g'): (165, 38, 4, 1, 160),
    (190, 'f'): (175, 34, 4, 1, 48),
    (190, 'g'): (145, 48, 6, 0, 160),
    (195, 'f'): (175, 34, 3, 1, 160),
    (195, 'g'): (175, 48, 6, 0, 40),
    (200, 'f'): (175, 34, 4, 1, 80),
    (200, 'g'): (155, 48, 6, 0, 160),
    (205, 'f'): (155, 42, 4, 1, 160),
    (205, 'g'): (185, 40, 4, 1, 192),
    (210, 'f'): (165, 38, 4, 1, 160),
    (210, 'g'): (165, 48, 6, 0, 160),
    (215, 'f'): (165, 43, 4, 1, 160),
    (215, 'g'): (215, 35, 4, 1, 160),
    (220, 'f'): (175, 39, 4, 1, 160),
    (220, 'g'): (205, 44, 6, 1, 20),
    (225, 'f'): (175, 43, 4, 1, 160),
    (225, 'g'): (215, 41, 4, 1, 160),
    (230, 'f'): (195, 48, 4, 0, 192),
    (230, 'g'): (205, 41, 6, 1, 40),
    (235, 'f'): (185, 40, 4, 1, 192),
    (235, 'g'): (225, 42, 6, 1, 20),
    (240, 'f'): (195, 40, 4, 1, 160),
    (240, 'g'): (235, 48, 8, 0, 8),
    (245, 'f'): (205, 41, 6, 1, 20),
    (245, 'g'): (235, 42, 6, 1, 20),
    (250, 'f'): (195, 44, 4, 1, 192),
    (250, 'g'): (245, 43, 4, 1, 160),
    (255, 'f'): (225, 48, 4, 0, 192),
    (255, 'g'): (225, 48, 6, 0, 80),
    (260, 'f'): (215, 41, 4, 1, 160),
    (260, 'g'): (245, 48, 6, 0, 40),
    (265, 'f'): (235, 48, 8, 0, 8),
    (265, 'g'): (235, 48, 6, 0, 80),
    (270, 'f'): (235, 48, 8, 0, 8),
    (270, 'g'): (255, 48, 6, 0, 40),
    (275, 'f'): (225, 42, 4, 1, 192),
    (275, 'g'): (265, 48, 8, 0, 8),
    (280, 'f'): (235, 42, 4, 1, 160),
    (280, 'g'): (275, 46, 6, 1, 12),
    (285, 'f'): (215, 48, 6, 0, 80),
    (285, 'g'): (265, 46, 6, 1, 20),
    (290, 'f'): (235, 45, 4, 1, 192),
    (290, 'g'): (255, 45, 6, 1, 40),
    (295, 'f'): (255, 48, 6, 0, 20),
    (295, 'g'): (295, 40, 6, 1, 20),
    (300, 'f'): (225, 45, 6, 1, 40),
    (300, 'g'): (285, 44, 6, 1, 20),
    (305, 'f'): (255, 45, 4, 1, 160),
    (305, 'g'): (265, 48, 6, 0, 160),
}
# boot plan: (n1, iA, kA, R)
BOOT16 = (245, 43, 12, 4)


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
