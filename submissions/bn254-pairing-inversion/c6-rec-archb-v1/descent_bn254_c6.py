"""
descent_bn254_c6.py

Pure-stdlib live 62-row bottom-up 12D descent recurrence & boot evaluator for
C6-m184_w216-rec at LB = 59, ORBIT = 6, with exact C6 base_primitive = 0.683542,
exact stored-prefix ECM extraction pricing, exact scaled-remainder-tree strip pricing,
exact PGF-envelope failure certificates (phi <= 0.04, psi < 0.8), exact E[N] recurrence,
and strict RT single-side tau-orbit repeated-ideal bound (li(2^59) / 6).
"""

import math, functools
from math import exp, factorial, floor, ceil, log, log1p, log2, sqrt
from fractions import Fraction as F
from operator import mul as _mul

LN2 = log(2.0)
UNIT = 20000.0
RSTEP = 1.0 / 1024.0
KA_MAX = 16
MONT5 = 367.76
TORSION = 12
WHEELS = (30030,)
SEPARATION = 0.5
RET_BOOT = 0.9
SURV_X = 2.0
SURV_DM = 10.0
CHILD_X = 1.25
LB_WAVE = 44

# BN254 12D norm-envelope bins (conservatively retaining bn15r's 12D boot & special-q bins
# with ZERO credit for PARI/GP-measured -1.71b boot and -8.33b special-q advantages of C6-m184_w216-rec!)
LOG2P6 = 6.0 * log2(21888242871839275222246405745257275088696311157297823662689037894645226208583)
LOG2U6 = 6.0 * log2(4965661367192848881)
QBINS = [(1/6, 13.6, 14.6), (1/6, 23.7, 14.6), (1/6, 31.5, 14.6), (1/6, 38.7, 14.6), (1/6, 47.6, 14.6), (1/6, 98.3, 14.6)]
QBINS_LO = [(-92.7, -53.4), (-6.5, -53.4), (3.6, -53.4), (11.4, -53.4), (18.6, -53.4), (27.5, -53.4)]
QBINS_SUB = [[-1.9, 5.4], [5.4, 4.7], [10.0, 6.0], [13.6, 6.4], [16.6, 6.7], [19.2, 7.6], [21.6, 8.5], [23.7, 8.3],
             [25.8, 9.4], [27.8, 10.0], [29.6, 9.4], [31.5, 10.2], [33.2, 10.4], [35.0, 10.2], [36.8, 11.3], [38.7, 11.3],
             [40.7, 11.5], [42.8, 12.6], [45.0, 11.2], [47.6, 12.2], [50.4, 13.6], [53.9, 12.7], [59.2, 11.6], [98.3, 14.6]]
BOOT_BINS = [(0.5, -178.78), (0.4, -175.78), (0.1, -171.93)]
BOOT_SUB = [-189.12, -186.97, -185.57, -184.44, -183.37, -182.43, -181.41, -180.27, -178.78, -171.93]

DF = 48
ALPHA_F = 0.0
ALPHA_G = 3.4
SLACK = 96.0
THIN_MAX = 72
UDUP = 0.6869
MESH = 96
PRCAP_NODE = 4.0
PRCAP_BOOT = 1.5
ECM_SUCCESS = 0.999
BOUND_CELL = 2400.0
HITS_PER_CELL = 2.0 * log(58.0 / 8.0)
SIEVE_HIT = 1100.0
SCAN = 10.0
ZERO = 2.0
PAIR_SETUP = 2.0 ** 25.1
LLL_NODE = 2.0 ** 36
NORMS = 2.0 ** 26
TDIV = 2.0 ** 32
IDEAL = 2.0 ** 34
PRIME = 2.0 ** 24
GRID = 4
SL_RANGE = (328, 380)

def _log2li(l):
    x = l * LN2
    return l - log2(x) + log2(1.0 + 1.0 / x + 2.0 / (x * x))

def _mulmod_instr(nbits):
    L = max(5, ceil(nbits / 64.0))
    return MONT5 * (L / 5.0) ** 2, L

def _gcd_instr(nbits):
    L = max(5, ceil(nbits / 64.0))
    return 2.0 * nbits * (12.0 * L + 16.0)

def _pi_upper(x):
    return 1.26 * x / log(x) if x >= 17 else 7.0

def _stage_ops(b1, b2, wheel):
    s1 = ceil(1.5 * b1) * 11
    phi = wheel
    for q in (2, 3, 5, 7, 11, 13):
        if wheel % q == 0:
            phi = phi // q * (q - 1)
    J = phi // 2
    G = ceil(b2 / wheel) + 1
    s2 = 6 * (wheel // 2 + G) + 3 * J * G if b2 > b1 else 0
    return 1.1 * (s1 + s2), s1, s2, J, G

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

def _convd(A, B):
    o0 = _conv(A[0], B[0])
    out = [o0]
    for a, b in zip(A[1:], B[1:]):
        x = _conv(a, B[0])
        y = _conv(A[0], b)
        L = max(len(x), len(y))
        out.append([(x[i] if i < len(x) else 0.0) + (y[i] if i < len(y) else 0.0) for i in range(L)])
    return out

def failure_bound(mu, shape, psi, kmax=64):
    if mu == 0:
        return dict(phi=1.0, raw=1.0, k=1, t=0.0)
    beta = 0.5 * mu / shape
    best = dict(raw=1.0, k=None, t=None)
    for k in range(1, kmax + 1):
        if k == 1:
            t = 0.0
            logtail = -shape * log1p(beta)
        else:
            t = min(0.8, ((k - 1) / (k - 1 + shape)) * (1.0 + 1.0 / beta))
            logtail = -shape * log1p(beta * (1.0 - t)) - (k - 1) * log(t)
        tail = exp(min(0.0, logtail))
        raw = min(1.0, math.nextafter(tail + psi ** k + 1e-10, math.inf))
        if raw < best["raw"]:
            best = dict(raw=raw, k=k, t=t)
    return dict(best, phi=max(0.01, best["raw"]))

def evaluate_descent(desc_data, rho_fn, fft_F_fn, gbatch_cost_fn, pbits_fn,
                     log2_rels, peel_h, nf_main=405.0, ng_main=469.0):
    """Execute the COMPLETE 62-row LB=59 bottom-up descent recurrence & boot LIVE!"""
    LB = 59.0
    ORBIT = 6
    CELLS = {int(k): v for k, v in desc_data["cells"].items()}
    OMIT = desc_data["omitted_law"]
    BASE_PRIM = desc_data["base_prim_c6"]
    PLAN = desc_data["plan59"]
    BOOT = tuple(desc_data["boot"])
    EA = tuple(desc_data["ea"])

    # Precompute high-res rho table on RSTEP = 1/1024 up to u = 32 (fast, < 0.03s!)
    RHO_TAB = [rho_fn(i * RSTEP) for i in range(32 * 1024 + 2)]
    def rho_lo(u):
        if u <= 1.0: return 1.0
        i = ceil(u / RSTEP)
        return RHO_TAB[i] if i < len(RHO_TAB) else 0.0
    def rho_hi(u):
        if u <= 1.0: return 1.0
        i = floor(u / RSTEP)
        if i + 1 >= len(RHO_TAB): return 1.0 / factorial(floor(u))
        return RHO_TAB[i]
    def _rho(rem, base):
        k = ceil(rem / (base * RSTEP))
        return 1.0 if k <= 0 else (RHO_TAB[k] if k < len(RHO_TAB) else 0.0)
    def kappa(u, t):
        ru = rho_fn(u)
        return rho_fn(u - t) / ru if ru > 0.0 else 0.0

    @functools.lru_cache(None)
    def semismooth(fbits, b1, b2, steps=256):
        m = fbits - log2(TORSION)
        base, top = log2(b1), log2(b2)
        acc = rho_lo(m / base)
        if top > base:
            for i in range(steps):
                lo = base + (top - base) * i / steps
                hi = base + (top - base) * (i + 1) / steps
                acc += rho_lo((m - lo) / base) * log(hi / lo)
        return min(acc, 1.0 - 1e-15)

    @functools.lru_cache(None)
    def power_corr(fbits, b1, b2):
        m = fbits - log2(TORSION)
        u = m / log2(b1)
        tilt = rho_hi(max(u - 2.0, 1.0)) / max(rho_lo(u), 1e-300)
        y = sqrt(b1)
        return max(0.0, 1.0 - tilt * 2.0 / (y * log(y)))

    @functools.lru_cache(None)
    def extraction_price(fbits, nbits, success=0.999):
        nb = ceil(nbits)
        mul, L = _mulmod_instr(nb)
        g = _gcd_instr(nb)
        F_ch = max(1, ceil(nb / 20.0))
        best = None
        for k in range(10, 37, 2):
            b1 = 1 << k
            Q = ceil(_pi_upper(b1))
            queries = min(Q, 2 * F_ch * ceil(log2(Q)) + 2)
            for ratio in (1, 10, 30, 100, 300):
                b2 = b1 * ratio
                p = semismooth(fbits, b1, b2) * power_corr(fbits, b1, b2)
                if p <= 1e-12: continue
                curves = ceil(log1p(-success) / log1p(-SEPARATION * p))
                for wheel in (WHEELS if ratio > 1 else (0,)):
                    if wheel and wheel > b2: continue
                    if wheel: ops, _, _, _, G = _stage_ops(b1, b2, wheel)
                    else: s1 = ceil(1.5 * b1) * 11; G = 0; ops = 1.1 * s1
                    cost = curves * (ops * mul + 8 * L * (Q + G) + 2 * g + 5 * (queries + G) * g)
                    if best is None or cost < best:
                        best = cost
        return best

    @functools.lru_cache(None)
    def batch_per(fbits, nbits, logS):
        leaf = max(1, ceil(nbits / 64.0))
        pb = pbits_fn(0, fbits)
        maxL = max(0, min(ceil(log2((ceil(pb / 64.0) + 1) / leaf)), int(logS)))
        best_c, best_L = None, 0
        for L_i in range(0, maxL + 1):
            c = gbatch_cost_fn(pb, leaf=leaf, L=L_i)
            if best_c is None or c < best_c:
                best_c, best_L = c, L_i
        return best_c + _gcd_instr(nbits) + 12.0 * _mulmod_instr(nbits)[0], best_L

    def strip_price(fbits, nbits, Ns):
        e = extraction_price(fbits, nbits, 0.999)
        if Ns < 1.0: return e
        b, L = batch_per(fbits, nbits, log2(Ns))
        return min(e, b * ceil(Ns / (2.0**L)) * (2.0**L) / Ns)

    def ptree_cost(x):
        return sum((2.0**(x - 6.0 + 0.53 - j - 1)) * (fft_F_fn(2.0**(j + 1)) + 3.0 * (2.0**(j + 1))) for j in range(ceil(x - 6.0 + 0.53)))

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
                rows[jA, jB] = PA[jA] if jB == 0 else _conv(PA[jA], PB[1])
        return rows, h

    @functools.lru_cache(None)
    def ss_cum2(m, base, n, iA, mesh, alpha, kmax=KA_MAX):
        rows, h = window_rows(n, base, iA, mesh, kmax)
        out = [[0.0, 0.0] for _ in range(kmax + 1)]
        ntab = len(RHO_TAB)
        for (jA, jB), row in rows.items():
            j = jA + jB
            ifeas = min(len(row) - 1, int(floor((m - j * h - j * base) / h))) if h > 0 else len(row) - 1
            if ifeas < 0: continue
            k0 = (m - j * base + alpha) / (base * RSTEP)
            hk = h / (base * RSTEP)
            s = 0.0
            for i in range(ifeas + 1):
                if row[i]:
                    k = ceil(k0 - i * hk)
                    s += row[i] * (1.0 if k <= 0 else (RHO_TAB[k] if k < ntab else 0.0))
            out[jA][jB] = s
        cum = [[0.0, 0.0] for _ in range(kmax + 1)]
        for jB in range(2):
            acc = 0.0
            for jA in range(kmax + 1):
                acc += out[jA][jB]
                cum[jA][jB] = acc
        return tuple(tuple(x) for x in cum)

    def ss_lower2(m, base, n, iA, kA, kB, mesh, alpha):
        m = ceil(m)
        if n <= base: return rho_lo((m + alpha) / base)
        cum = ss_cum2(m, base, n, iA, mesh, alpha, kA)
        return cum[kA][0] + (cum[kA][1] if kB else 0.0)

    def ss_upper2(m, base, n, iA, kA, kB, mesh, alpha):
        m = floor(m)
        if n <= base: return rho_hi((m + alpha) / base)
        rows, h = window_rows(n, base, iA, mesh, kA)
        tot = 0.0
        for (jA, jB), row in rows.items():
            if jA > kA or (jB and not kB): continue
            j = jA + jB
            for i, w in enumerate(row):
                if w and j * base + i * h <= m:
                    u = (m - (j * base + i * h) - j * h + alpha) / base
                    tot += w * (1.0 if u <= 1.0 else rho_hi(u))
        return min(tot, 1.0)

    def surv_upper(m, base, T):
        if T <= 0: return 0.0
        if T >= m: return 1.0
        return min(1.0, max(1.0, T / base) * rho_hi((floor(m) - T) / base))

    # Exact C6 base_primitive = 0.683542 in retain16 counted_fraction & retention!
    ord_both = max(0.0, OMIT["f"][3] + OMIT["g"][3] - 1.0)
    COUNTED_FRAC = ord_both * BASE_PRIM

    def loss_side(side, m, n, iA, kA, kB, alpha, slack):
        s = int(floor((slack - 1.0 - m / 128.0) * GRID))
        b0, fine, blk, _ = OMIT[side]
        hi = lambda xg: ss_upper2(max(m - xg, 1.0), LB, n, iA, kA, kB, MESH, alpha) if m - xg > 1.0 else 1.0
        num = sum(fine[max(0, s - SL_RANGE[0]) :]) * hi(128.0) + sum(v * hi(64.0 * b) for b, v in enumerate(blk, 2))
        return num / (ss_lower2(m, LB, n, iA, kA, kB, MESH, alpha) * b0)

    def retention_node(nf, ng, n, iA, kA, kB, slack):
        lf = loss_side("f", ceil(nf), n, iA, kA, kB, ALPHA_F, slack)
        lg = loss_side("g", ceil(ng), n, iA, kA, kB, ALPHA_G, slack)
        return max(0.0, 1.0 - lf - lg) * COUNTED_FRAC

    @functools.lru_cache(None)
    def miss_side(u4, u_main):
        u = u4 / 4.0
        tot = 0.0
        for b10 in range(300, int(10 * LB)):
            b = b10 / 10.0
            c = b + 0.1
            lam = (2.0 ** (log2_rels - c)) * ORBIT * kappa(u_main, c / LB)
            tot += kappa(u, c / LB) * log(c / b) * exp(-lam * (1.0 - (peel_h if b10 >= 544 else 0.4)))
        return tot

    u_main_f, u_main_g = nf_main / LB, ng_main / LB
    def missing(nf, ng):
        return (miss_side(ceil(4.0 * nf / LB), u_main_f) if nf else 0.0) + (miss_side(ceil(4.0 * ng / LB), u_main_g) if ng else 0.0)

    pairs = 2.0 * (2.0 ** _log2li(LB))
    cfmax = max(b[1] for b in QBINS)
    cgmax = max(b[2] for b in QBINS)

    # Execute ALL 62 rows bottom-up!
    W = {}; PHI = {}; DET = {}; CC = {}
    def wm(s, n): return max([v for (b, ss), v in W.items() if ss == s and b[0] < n - 1e-9] + [0.0])
    def pm(s, n): return max([v for (b, ss), v in PHI.items() if ss == s and b[0] < n - 1e-9] + [0.0])
    def cm(s, n): return max([v for (b, ss), v in CC.items() if ss == s and b[0] < n - 1e-9] + [0.0])

    for b_lo, b_hi, side, n_i, iA, kA, kB, R, thin in PLAN:
        n = float(n_i)
        S = log2(CELLS[R])
        c = (b_hi + S) / 12.0
        c_lo = (b_lo + S) / 12.0
        h_w = (n - LB) / MESH if n > LB else 0.0
        nA = LB + iA * h_w if n > LB else LB
        T = kA * nA + kB * n + SLACK
        T0 = ceil(kA * (F(58) + (F(iA * (n_i - 58), 96) if n > LB else 0)) + kB * n_i)

        pr = sv = sv0 = 0.0
        for (w, cf, cg), (cf_lo, cg_lo) in zip(QBINS, QBINS_LO):
            nf = DF * c + cf - (b_lo if side == "f" else 0.0)
            ng = 12.0 * c + LOG2U6 + cg - (b_lo if side == "g" else 0.0)
            pr += w * ss_lower2(nf, LB, n, iA, kA, kB, MESH, ALPHA_F) * ss_lower2(ng, LB, n, iA, kA, kB, MESH, ALPHA_G)
            nf_lo = DF * c_lo + cf_lo - (b_hi if side == "f" else 0.0)
            ng_lo = 12.0 * c_lo + LOG2U6 + cg_lo - (b_hi if side == "g" else 0.0)
            sv += w * surv_upper(nf_lo, LB, T) * surv_upper(ng_lo, LB, T)
            sv0 += w * surv_upper(nf_lo, LB, T0) * surv_upper(ng_lo, LB, T0)
        pri = 0.0
        for cf, cg in QBINS_SUB:
            nf = DF * c + cf - (b_lo if side == "f" else 0.0)
            ng = 12.0 * c + LOG2U6 + cg - (b_lo if side == "g" else 0.0)
            pri += ss_lower2(nf, LB, n, iA, kA, kB, MESH, ALPHA_F) * ss_lower2(ng, LB, n, iA, kA, kB, MESH, ALPHA_G)
        pr = min(pri / 24.0, PRCAP_NODE * pr)

        nf_max = DF * c + cfmax - (b_lo if side == "f" else 0.0)
        ng_max = 12.0 * c + LOG2U6 + cgmax - (b_lo if side == "g" else 0.0)
        det = 1.0 - (1.0 - ECM_SUCCESS) * ((nf_max + ng_max) / 20.0 + 2.0 * (kA + kB))
        miss = missing(nf_max, ng_max)
        is_thin = bool(thin) or (b_lo <= THIN_MAX)
        mu = CELLS[R] * pr * UDUP * det * (exp(-miss) if is_thin else 1.0)
        ef = extraction_price(LB, nf_max, ECM_SUCCESS)
        eg = extraction_price(LB, ng_max, ECM_SUCCESS)
        sc_ecm = NORMS + TDIV + ef + eg + 2.0 * (extraction_price(nA, T, ECM_SUCCESS) + PRIME) + IDEAL
        sieve = CELLS[R] * (HITS_PER_CELL * SIEVE_HIT + SCAN + ZERO + BOUND_CELL) + pairs * PAIR_SETUP
        if CELLS[R] * sc_ecm < sieve + CELLS[R] * (sv + pr) * sc_ecm:
            how = "ecm"
            mu *= COUNTED_FRAC
            Ns = CELLS[R]
        else:
            how = "sieve"
            mu *= retention_node(nf_max, ng_max, n, iA, kA, kB, SLACK)
            Ns = CELLS[R] * (sv + pr)
        strip = NORMS + TDIV + strip_price(LB, nf_max, Ns) + strip_price(LB, ng_max, Ns) + IDEAL
        ext = 2.0 * (extraction_price(nA, T0, ECM_SUCCESS) + PRIME) + ef + eg
        gated = min(CELLS[R], CELLS[R] * (sv0 + pr))
        own = (LLL_NODE + CELLS[R] * strip + gated * ext if how == "ecm" else LLL_NODE + sieve + Ns * strip + gated * ext) + 2**40

        miss_ch = 0.0 if is_thin else miss
        lowp = max([v for (bb, ss), v in PHI.items() if bb[0] <= THIN_MAX] + [0.0])
        loww = max([v for (bb, ss), v in W.items() if bb[0] <= THIN_MAX] + [0.0])
        psi = kA * (pm("f", nA) + pm("g", nA)) + kB * (pm("f", n) + pm("g", n)) + miss_ch * lowp
        assert psi < 0.8, f"psi={psi} >= 0.8 at ({b_lo},{side})"
        win_cnt = (kA if nA > LB + 1e-9 else 0) + (kB if n > nA + 1e-9 else 0)
        r_sh = 0.5 if win_cnt <= 1 else 1.0
        cert = failure_bound(mu, r_sh, psi)
        phi = cert["phi"]
        assert phi <= 0.04 + 1e-12, f"phi={phi} > 0.04 at ({b_lo},{side})"
        sub = kA * (wm("f", nA) + wm("g", nA)) + kB * (wm("f", n) + wm("g", n)) + miss_ch * loww
        work = own + sub / (1.0 - psi)
        bkey = (b_lo, b_hi)
        W[bkey, side] = work
        PHI[bkey, side] = phi
        retry = 1.0 / (1.0 - psi)
        low_c = max([v for (bb, ss), v in CC.items() if bb[0] <= THIN_MAX] + [0.0])
        CC[bkey, side] = 1.0 + retry * (kA * (cm("f", nA) + cm("g", nA)) + kB * (cm("f", n) + cm("g", n)) + miss_ch * low_c)
        DET[bkey, side] = dict(how=how, R=R, mu=mu, psi=psi, phi=phi, W_log2=log2(work / UNIT))

    # Live Boot Evaluation (Gate A + Gate E)!
    n_b, iA_b, kA_b, Rb = BOOT
    kx, tau = EA
    hb = (n_b - LB) / MESH
    nA_b = LB + iA_b * hb
    x_ea = LB + kx * hb
    c_b = (LOG2P6 + log2(CELLS[Rb])) / 12.0
    m_max = 12.0 * c_b + LOG2U6 + max(v[1] for v in BOOT_BINS)
    m_min = 12.0 * c_b + LOG2U6 + min(BOOT_SUB)
    tops = [LB + (k + 1) * hb for k in range(MESH)]
    sc_w = max(W.values())
    fs = (tuple(wm("g", t) / sc_w for t in tops), tuple(pm("g", t) for t in tops))

    @functools.lru_cache(None)
    def ea_rows(base, n, iA, kA, mesh, kx_i, tau_i, fs_t):
        h = (n - base) / mesh
        w = [log1p(h / (base + k * h)) for k in range(mesh)]
        C = len(fs_t)
        def part(sel):
            return [[w[k] if sel(k) else 0.0 for k in range(mesh)]] + [[w[k] * f[k] if sel(k) else 0.0 for k in range(mesh)] for f in fs_t]
        wl = part(lambda k: k < kx_i); wh = part(lambda k: kx_i <= k < iA); wb = part(lambda k: k >= iA)
        one = [[1.0]] + [[0.0]] * C
        PL = [one]; PH = [one]
        for j in range(1, kA + 1):
            PL.append([[x / j for x in v] for v in _convd(PL[-1], wl)])
            PH.append([[x / j for x in v] for v in _convd(PH[-1], wh)])
        def add(X, Y):
            return [[(x[i] if i < len(x) else 0.0) + (y[i] if i < len(y) else 0.0) for i in range(max(len(x), len(y)))] for x, y in zip(X, Y)]
        F_tab = {}
        for jj in range(kA + 2):
            lim = floor((tau_i + 1e-9 - jj * base) / h - jj)
            if lim < 0: continue
            cut = lambda X: [v[: max(0, min(len(v), lim + 1))] for v in X]
            H0 = cut(PH[jj]) if jj <= kA else None
            H1 = cut(_convd(PH[jj - 1], wb)) if jj >= 1 else None
            for jl in range(kA - jj + 2):
                ps = ([H0] if H0 is not None and jl <= kA - jj else []) + ([H1] if H1 is not None and jl <= kA - jj + 1 else [])
                if not ps: continue
                H = ps[0] if len(ps) == 1 else add(*ps)
                if not any(H[0]): continue
                r = _convd(PL[jl], H)
                f = F_tab.setdefault(jl + jj, [[] for _ in r])
                for fc, rc in zip(f, r):
                    if len(fc) < len(rc): fc.extend([0.0] * (len(rc) - len(fc)))
                    for i, x in enumerate(rc): fc[i] += x
        return F_tab, h

    def ea_mass(m, fs_t):
        F_tab, h = ea_rows(LB, float(n_b), iA_b, kA_b, MESH, kx, tau, fs_t)
        tot = [0.0] * (len(fs_t) + 1)
        for j, rows in F_tab.items():
            L = min(max(len(r) for r in rows) - 1, floor((m - j * h - j * LB) / h)) + 1
            if L <= 0: continue
            rr = [_rho(m - j * LB - i * h + ALPHA_G, LB) for i in range(L)]
            for c_i, row in enumerate(rows):
                tot[c_i] += sum(map(_mul, row[:L], rr))
        return tot

    def psurv(m, x, tau_v, mesh=192):
        h = (tau_v - x) / mesh
        w = [log1p(h / (x + k * h)) for k in range(mesh)]
        tot = rho_hi(m / x)
        P = [1.0]
        for j in range(1, int(tau_v // x) + 1):
            P = [y / j for y in _conv(P, w)]
            tot += sum(p * rho_hi(max(m - j * x - (i + j) * h, 0.0) / x) for i, p in enumerate(P) if p and j * x + i * h <= tau_v + 1e-9)
        return min(1.0, tot)

    pr_b = sum(w * ea_mass(ceil(12.0 * c_b + LOG2U6 + cg), fs)[0] for w, cg in BOOT_BINS)
    pri_b = sum(ea_mass(ceil(12.0 * c_b + LOG2U6 + cg), fs)[0] for cg in BOOT_SUB) / len(BOOT_SUB)
    pr_b = min(pri_b, PRCAP_BOOT * pr_b)
    det_b = 1.0 - 0.001 * (m_max / 20.0 + 2.0 * (kA_b + 1))
    prob_b = RET_BOOT * pr_b * det_b * UDUP * COUNTED_FRAC * 0.5
    lowp_b = max(v for (bb, s), v in PHI.items() if bb[0] <= THIN_MAX)
    loww_b = max(v for (bb, s), v in W.items() if bb[0] <= THIN_MAX)
    miss_b = missing(0.0, m_max)
    ew = ep = 0.0
    for mc in sorted({ceil(12.0 * c_b + LOG2U6 + cg) for cg in [v[1] for v in BOOT_BINS] + list(BOOT_SUB)}):
        t_m = ea_mass(mc, fs)
        ew = max(ew, (t_m[1] / t_m[0]) * sc_w)
        ep = max(ep, t_m[2] / t_m[0])
    psi_b = CHILD_X * ep + miss_b * lowp_b
    sub_b = CHILD_X * ew + miss_b * loww_b
    assert psi_b < 0.8
    surv_b = min(1.0, SURV_X * psurv(m_min - SURV_DM, x_ea, tau))
    sb_b = strip_price(x_ea, m_max, 2.0 ** LB_WAVE)
    trial_b = 2**40 + NORMS + sb_b + surv_b * (extraction_price(x_ea, m_max, 0.999) + extraction_price(nA_b, tau, 0.999)) + PRIME + IDEAL
    over_b = (2.0 ** LB_WAVE) * trial_b / (1.0 - psi_b)
    tree_b = 2.0 * ptree_cost(x_ea)
    one_b = (trial_b / prob_b + sub_b) / (1.0 - psi_b) + over_b + tree_b
    two_log_u = log2(2.0 * one_b / UNIT)

    # Live E[N] and strict RT repeated-ideal collision bound (single-side li(2^59)/6)!
    low_c_b = max(v for (bb, s), v in CC.items() if bb[0] <= THIN_MAX)
    en = 2.0 * (1.0 / (1.0 - psi_b)) * (kA_b * cm("g", nA_b) + cm("g", float(n_b)) + miss_b * low_c_b)
    logen = log2(en)
    small_orbits = _log2li(LB) - log2(6.0)
    budget = min((x / 1000.0 for x in range(18000, 26001)),
                 key=lambda x: 2.0**(logen - x) + 2.0**(2.0 * x - 1.0 - small_orbits))

    return dict(
        ok=True,
        nodes=len(DET),
        lb=LB,
        base_prim=BASE_PRIM,
        counted_frac=COUNTED_FRAC,
        two_log=two_log_u,
        boot_psi=psi_b,
        boot_cfg=list(BOOT),
        phi_max=max(PHI.values()),
        psi_max=max(d["psi"] for d in DET.values()),
        mu_min=min(d["mu"] for d in DET.values()),
        log2_EN=logen,
        budget=budget,
        small_orbits_rt=small_orbits,
        abort_log2=logen - budget,
        repeat_log2=2.0 * budget - 1.0 - small_orbits,
    )
