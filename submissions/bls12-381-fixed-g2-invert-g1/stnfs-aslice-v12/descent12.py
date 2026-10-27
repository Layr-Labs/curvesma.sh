"""Descent cost by dynamic programming over 1-bit size levels (stnfs v5; v12: ECM priced with the measured ecm12 Montgomery multiplication through model12.ecm_instr, and no silently uncosted child bins)."""
import math, sys
import model12
import filtsim
rho = filtsim.rho   # tabulated Dickman rho (fast; the DP evaluates it millions of times)
LN2 = math.log(2)
ECM_TAB = model12.ECM_TAB
ecm_instr = model12.ecm_instr
def proba(log2n, alpha, log2b):
    u = (log2n + alpha) / log2b
    return rho(u) + (1 - 0.5772) * rho(u - 1) / (log2n * LN2)
def kappa(u, t):
    r = rho(u); return rho(u - t) / r if r > 0 else 0.0
def log2li(l):
    x = l * LN2; return l - math.log2(x) + math.log2(1 + 1 / x + 2 / x / x)
def pois_cdf(lam, k):
    s = 0.0; term = math.exp(-lam)
    for i in range(k + 1): s += term; term *= lam / (i + 1)
    return s

class P:
    lB = 68.5; alpha_f = 4.87; alpha_g = 0.84; CF = 54.0; CG = 2.0; log2U6 = 6 * 63.714   # v12: CF = 54 from box-uniform coefficients (72 c + 50..54 measured); v10 used 6
    BOOT_NORM = 2700.0; LLL12 = 2.0 ** 30; LLL6 = 2.0 ** 13; NORM = 2.0 ** 21
    SIEVE_PP = 250.0; LX = 65.0; NSIEVE_MAX = 80.0; PAD_CHILD = 1.5; PAD_NODE = 2.0
    LMAX = 360

def children_bins(nf, ng, n, lB, pf_alpha, pg_alpha):
    """expected children per 1-bit bin [b, b+1) for b in [ceil(lB), n), per side"""
    out = []
    uf = (nf + pf_alpha) / n; ug = (ng + pg_alpha) / n
    b = math.floor(lB)
    while b < n:
        lo, hi = max(b, lB), min(b + 1, n)
        if hi > lo:
            t = 0.5 * (lo + hi) / n; w = (hi - lo) / n / t
            out.append((b, kappa(uf, t) * w, kappa(ug, t) * w))
        b += 1
    return out

# ---------------- the DP over 2-bit bins (the sandbox and the claim use solve_fast) ----------------
def solve_fast(Pm, verbose=False, s_step=1):
    lB = Pm.lB
    bins = list(range(int(math.floor(lB)), Pm.LMAX + 1, 2))
    cost = {}; detail = {}
    S_grid = [float(x) for x in range(28, 100, s_step)]
    for b in bins:
        lo_b = max(float(b), lB)
        l = 0.5 * (lo_b + b + 2.0)
        for side in ('f', 'g'):
            if b + 2.0 <= lB + 1e-9:
                cost[(b, side)] = 0.0; continue
            best = None

            targets = [lB] + [float(x) for x in range(bins[0] + 2, b + 1, 2) if x > lB]
            for nn in targets:
                if nn >= l: continue
                opt = node_opt(l, side, nn, Pm, S_grid)
                if opt is None: continue
                c, S, nf, ng, keep, how = opt
                kids = children_bins(nf, ng, nn, lB, Pm.alpha_f, Pm.alpha_g)
                sub = 0.0
                for (bb, kf, kg) in kids:
                    key = bb - ((bb - bins[0]) % 2)
                    if (key, 'f') not in cost:
                        raise AssertionError('uncosted child bin %s for node bin %s' % (key, b))
                    sub += Pm.PAD_CHILD * (kf * cost[(key, 'f')] + kg * cost[(key, 'g')])
                tot = Pm.PAD_NODE * c + sub
                if best is None or tot < best[0]:
                    best = (tot, nn, c, S, nf, ng, keep, how, sum(Pm.PAD_CHILD * (kf + kg) for (_, kf, kg) in kids))
            cost[(b, side)] = best[0]; detail[(b, side)] = best
    bestb = None
    for lB1 in range(150, Pm.LMAX, 10):
        pr = proba(Pm.BOOT_NORM, 0.0, lB1)
        if pr <= 0: continue
        cand = 1 / pr
        own = cand * (Pm.NORM + ecm_instr(lB1, Pm.BOOT_NORM)) + Pm.LLL12 * max(1, cand / 2 ** 30)
        kids = children_bins(Pm.BOOT_NORM, 1.0, lB1, lB, 0.0, 0.0)
        sub = 0.0
        for (bb, kf, kg) in kids:
            key = bb - ((bb - bins[0]) % 2); sub += Pm.PAD_CHILD * kf * cost[(key, 'g')]
        tot = Pm.PAD_NODE * own + sub
        if bestb is None or tot < bestb[0]: bestb = (tot, lB1, own, sub, cand, sum(Pm.PAD_CHILD * kf for (_, kf, _) in kids))
    return math.log2(bestb[0]), bestb, detail

def node_opt(l, side, n, Pm, S_grid):
    best = None
    c0 = l / 12
    for S in S_grid:
        c = c0 + S / 12
        nf = 72 * c + Pm.CF; ng = 12 * c + Pm.log2U6 + Pm.CG
        if side == 'f': nf -= l
        else: ng -= l
        pr = proba(nf, Pm.alpha_f, n) * proba(ng, Pm.alpha_g, n)
        if pr <= 0: continue
        if 2 ** S * pr >= 4:
            per = Pm.NORM + ecm_instr(n, nf) + ecm_instr(n, ng)
            c1 = 2 ** S * per + Pm.LLL12
            if best is None or c1 < best[0]: best = (c1, S, nf, ng, 1.0, 'ecm')
            if n > Pm.NSIEVE_MAX: break
        if n <= Pm.NSIEVE_MAX:
            for k in (3, 4):
                lam_f = lam_top(nf, Pm.alpha_f, n, Pm.LX); lam_g = lam_top(ng, Pm.alpha_g, n, Pm.LX)
                keep = pois_cdf(lam_f, k) * pois_cdf(lam_g, k)
                if 2 ** S * pr * keep < 4: continue
                surv = rho(max(nf - k * n, 0) / Pm.LX) * rho(max(ng - k * n, 0) / Pm.LX)
                pairs = 2 ** log2li(Pm.LX) / 2
                c2 = pairs * Pm.LLL6 + 2 ** S * (Pm.SIEVE_PP + surv * (Pm.NORM + ecm_instr(n, nf) + ecm_instr(n, ng))) + Pm.LLL12
                if best is None or c2 < best[0]: best = (c2, S, nf, ng, keep, 'sieve%d' % k)
    return best

_LAM = {}
def lam_top(nbits, alpha, n, LX):
    key = (round(nbits, 1), n)
    if key in _LAM: return _LAM[key]
    u = (nbits + alpha) / n; t0 = LX / n
    if t0 >= 1: v = 0.0
    else:
        m = 40; s = 0.0
        for i in range(m):
            t = t0 + (1 - t0) * (i + 0.5) / m; s += kappa(u, t) * ((1 - t0) / m) / t
        v = s
    _LAM[key] = v
    return v
