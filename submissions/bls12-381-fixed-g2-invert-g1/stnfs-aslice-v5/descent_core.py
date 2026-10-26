"""Descent cost by dynamic programming over 1-bit size levels (stnfs v5).

A node = prime ideal of norm 2^l on side s.  Its special-q lattice (12-dim, det 2^l) is reduced once
(2^LLL12 instructions).  Candidates = 2^S lattice points in a ball; coordinates 2^(l/12 + S/12);
   N_f = 72 c + CF,  N_g = 12 c + 6 log2 U + CG,  the q-side norm divided by 2^l.
The node picks a target bound 2^n (n < l): both norms must be 2^n-smooth; S is the least with 2^S P >= 4.
Candidate test:
   (ECM)   exact norms + full ECM to 2^n on both cofactors (no sieving: for n >> 65 the factor-base sieve
           carries no information), cost per candidate = NORM + ecm(n, N_f) + ecm(n, N_g);
   (SIEVE) only when n <= NSIEVE_MAX: sieve ideals <= 2^65 (per-(node, ideal) 6-dim lattice setup charged:
           2^LLL6 x pairs), per-point SIEVE_PP, then cofactor survivors (slack S_s = k log2 n per side) with
           ECM(n); relations with more than k unsieved primes per side are lost (keep factor).
Children: every prime factor of norm in (2^lB, 2^n] of either norm is a child; expected count per 1-bit bin
from the Dickman-conditioned density kappa(t)/t, times PAD_CHILD; each child costs cost(bin, side).
Boot: lattice L_T (det p^6), N_g ~ 2^BOOT_NORM; B1 optimised; candidates 1/P(B1-smooth); ECM(B1) each;
children g-side in (lB, B1].
"""
import math, sys
import filtsim
rho = filtsim.rho
LN2 = math.log(2)
ECM_TAB = {20: (11e3, 86), 25: (5e4, 214), 30: (25e4, 430), 35: (1e6, 904), 40: (3e6, 2350), 45: (11e6, 4480),
           50: (43e6, 7553), 55: (11e7, 17769), 60: (26e7, 42017), 65: (85e7, 69408), 70: (29e8, 102212), 75: (76e8, 188056), 80: (25e9, 265557)}
def ecm_instr(fbits, nbits, success=0.99):
    d = fbits / 3.3219
    ks = sorted(ECM_TAB)
    if d <= ks[-1]:
        k = min(x for x in ks if x >= d); b1, cur = ECM_TAB[k]
    else:   # extrapolate: B1 x3.3 and curves x1.5 per 5 digits (matches the table's trend)
        b1, cur = ECM_TAB[ks[-1]]; e = (d - ks[-1]) / 5; b1 *= 3.3 ** e; cur *= 1.5 ** e
    curves = cur * math.log(1 / (1 - success))
    L = max(1, math.ceil(nbits / 64)); mulmod = 2.2 * L * L + 40
    return curves * 2 * 1.44 * b1 * 11 * mulmod
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
    lB = 68.5; alpha_f = 4.87; alpha_g = 0.84; CF = 6.0; CG = 2.0; log2U6 = 6 * 63.714
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

def node_options(l, side, n, Pm):
    """best (cost, S, nf, ng, keep) over S for a node at level l targeting 2^n"""
    best = None
    for S in [x * 0.5 for x in range(30, 170)]:
        c = l / 12 + S / 12
        nf = 72 * c + Pm.CF; ng = 12 * c + Pm.log2U6 + Pm.CG
        if side == 'f': nf -= l
        else: ng -= l
        pf = proba(nf, Pm.alpha_f, n); pg = proba(ng, Pm.alpha_g, n); pr = pf * pg
        if pr <= 0: continue
        # ECM option
        if 2 ** S * pr >= 4:
            per = Pm.NORM + ecm_instr(n, nf) + ecm_instr(n, ng)
            c1 = 2 ** S * per + Pm.LLL12
            if best is None or c1 < best[0]: best = (c1, S, nf, ng, 1.0, 'ecm')
        # SIEVE option
        if n <= Pm.NSIEVE_MAX:
            for k in (2, 3, 4):
                # unsieved primes per side in (LX, n]: Poisson(lambda); keep = P(count <= k) both sides
                def lam(nbits, alpha):
                    u = (nbits + alpha) / n; t0 = Pm.LX / n
                    if t0 >= 1: return 0.0
                    m = 100; s = 0.0
                    for i in range(m):
                        t = t0 + (1 - t0) * (i + 0.5) / m; s += kappa(u, t) * ((1 - t0) / m) / t
                    return s
                keep = pois_cdf(lam(nf, Pm.alpha_f), k) * pois_cdf(lam(ng, Pm.alpha_g), k)
                if 2 ** S * pr * keep < 4: continue
                Sf = k * n; Sg = k * n
                surv = rho(max(nf - Sf, 0) / Pm.LX) * rho(max(ng - Sg, 0) / Pm.LX)
                pairs = 2 ** log2li(Pm.LX) / 2
                c2 = pairs * Pm.LLL6 + 2 ** S * (Pm.SIEVE_PP + surv * (Pm.NORM + ecm_instr(n, nf) + ecm_instr(n, ng))) + Pm.LLL12
                if best is None or c2 < best[0]: best = (c2, S, nf, ng, keep, 'sieve%d' % k)
    return best

def solve(Pm, verbose=True):
    lB = Pm.lB
    levels = list(range(math.floor(lB), Pm.LMAX + 1))      # bin b covers [b, b+1)
    cost = {}   # (b, side) -> expected total instructions to finish a node whose size is in bin b
    detail = {}
    for b in levels:
        l = b + 0.5 if b >= lB else lB            # representative size
        if b < lB: l = lB
        for side in ('f', 'g'):
            if b <= lB - 1e-9 or (b == math.floor(lB) and l <= lB):
                cost[(b, side)] = 0.0; continue
            best = None
            for n in range(math.floor(lB), b + 1):      # target bound 2^n with n < l
                nn = n + 0.999 if n >= lB else lB
                nn = min(nn, l - 0.5)
                if nn <= lB and n > lB: continue
                opt = node_options(l, side, max(nn, lB), Pm)
                if opt is None: continue
                c, S, nf, ng, keep, how = opt
                kids = children_bins(nf, ng, max(nn, lB), lB, Pm.alpha_f, Pm.alpha_g)
                sub = sum(Pm.PAD_CHILD * (kf * cost[(bb, 'f')] + kg * cost[(bb, 'g')]) for (bb, kf, kg) in kids)
                tot = Pm.PAD_NODE * c + sub
                if best is None or tot < best[0]:
                    best = (tot, n, c, S, nf, ng, keep, how, sum(Pm.PAD_CHILD * (kf + kg) for (_, kf, kg) in kids))
            cost[(b, side)] = best[0]; detail[(b, side)] = best
    if verbose:
        for b in (69, 70, 72, 75, 80, 90, 100, 120, 150, 180, 220, 260, 300, 340):
            for side in ('f', 'g'):
                if (b, side) in detail:
                    tot, n, c, S, nf, ng, keep, how, nk = detail[(b, side)]
                    print("  node 2^%d %s: target 2^%d via %s S=%.1f Nf=%.0f Ng=%.0f keep=%.2f own 2^%.1f children %.1f subtree 2^%.1f" % (b, side, n, how, S, nf, ng, keep, math.log2(c), nk, math.log2(tot)))
    # boot
    bestb = None
    for lB1 in range(150, Pm.LMAX, 5):
        pr = proba(Pm.BOOT_NORM, 0.0, lB1)
        if pr <= 0: continue
        cand = 1 / pr
        own = cand * (Pm.NORM + ecm_instr(lB1, Pm.BOOT_NORM)) + Pm.LLL12 * max(1, cand / 2 ** 30)
        kids = children_bins(Pm.BOOT_NORM, 1.0, lB1, lB, 0.0, 0.0)
        sub = sum(Pm.PAD_CHILD * kf * cost[(bb, 'g')] for (bb, kf, kg) in kids)
        tot = Pm.PAD_NODE * own + sub
        if bestb is None or tot < bestb[0]: bestb = (tot, lB1, own, sub, cand, sum(Pm.PAD_CHILD * kf for (_, kf, _) in kids))
    tot, lB1, own, sub, cand, nk = bestb
    if verbose:
        print("  boot: B1=2^%d candidates 2^%.1f own 2^%.1f (x%g pad) tree 2^%.1f children %.1f -> per log 2^%.2f instr = 2^%.2f units" % (
            lB1, math.log2(cand), math.log2(own), Pm.PAD_NODE, math.log2(sub), nk, math.log2(tot), math.log2(tot) - math.log2(20000)))
    return math.log2(tot), bestb, detail

if __name__ == "__main__":
    Pm = P()
    if len(sys.argv) > 1: Pm.lB = float(sys.argv[1])
    r = solve(Pm)
    print("descent per log: 2^%.2f units; two logs: 2^%.2f units" % (r[0] - math.log2(20000), r[0] + 1 - math.log2(20000)))
