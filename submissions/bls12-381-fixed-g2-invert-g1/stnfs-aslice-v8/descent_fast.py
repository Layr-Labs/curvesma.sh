"""Fast descent DP (sandbox version of descent2): 2-bit size bins, S step 1, target bound stepped by 2 bits,
rho via the shared table.  Same model, coarser grid -> cost is an upper bound of the fine DP's optimum."""
import math, sys
import descent_core as d2
rho = d2.rho; LN2 = d2.LN2

def solve_fast(Pm, verbose=False):
    lB = Pm.lB
    bins = list(range(int(math.floor(lB)), Pm.LMAX + 1, 2))        # bin b covers [b, b+2)
    cost = {}; detail = {}
    S_grid = [float(x) for x in range(28, 100)]
    for b in bins:
        lo_b = max(float(b), lB)                       # the bin's part above the factor-base bound
        l = 0.5 * (lo_b + b + 2.0)                     # representative size of the nodes in the bin
        for side in ('f', 'g'):
            if b + 2.0 <= lB + 1e-9:                   # entirely inside the factor base: no descent
                cost[(b, side)] = 0.0; continue
            best = None
            targets = sorted({lB} | {float(n) for n in range(int(math.floor(lB)), b + 1, 2) if lB <= n < l})
            for n in targets:
                nn = max(min(n + 1.0, l - 0.5), lB)    # target bound 2^nn with lB <= nn < l
                if nn >= l: continue
                opt = node_opt(l, side, nn, Pm, S_grid)
                if opt is None: continue
                c, S, nf, ng, keep, how = opt
                kids = d2.children_bins(nf, ng, nn, lB, Pm.alpha_f, Pm.alpha_g)
                sub = 0.0
                for (bb, kf, kg) in kids:
                    key = bb - ((bb - bins[0]) % 2)
                    if key not in cost: continue       # bins above the current one cannot occur (children are smaller)
                    sub += Pm.PAD_CHILD * (kf * cost[(key, 'f')] + kg * cost[(key, 'g')])
                tot = Pm.PAD_NODE * c + sub
                if best is None or tot < best[0]:
                    best = (tot, n, c, S, nf, ng, keep, how, sum(Pm.PAD_CHILD * (kf + kg) for (_, kf, kg) in kids))
            cost[(b, side)] = best[0]; detail[(b, side)] = best
    bestb = None
    for lB1 in range(150, Pm.LMAX, 10):
        pr = d2.proba(Pm.BOOT_NORM, 0.0, lB1)
        if pr <= 0: continue
        cand = 1 / pr
        own = cand * (Pm.NORM + d2.ecm_instr(lB1, Pm.BOOT_NORM)) + Pm.LLL12 * max(1, cand / 2 ** 30)
        kids = d2.children_bins(Pm.BOOT_NORM, 1.0, lB1, lB, 0.0, 0.0)
        sub = 0.0
        for (bb, kf, kg) in kids:
            key = bb - ((bb - bins[0]) % 2); sub += Pm.PAD_CHILD * kf * cost.get((key, 'g'), 0.0)
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
        pr = d2.proba(nf, Pm.alpha_f, n) * d2.proba(ng, Pm.alpha_g, n)
        if pr <= 0: continue
        if 2 ** S * pr >= 4:
            per = Pm.NORM + d2.ecm_instr(n, nf) + d2.ecm_instr(n, ng)
            c1 = 2 ** S * per + Pm.LLL12
            if best is None or c1 < best[0]: best = (c1, S, nf, ng, 1.0, 'ecm')
            if n > Pm.NSIEVE_MAX: break          # larger S only costs more once feasible
        if n <= Pm.NSIEVE_MAX:
            for k in (3, 4):
                lam_f = lam_top(nf, Pm.alpha_f, n, Pm.LX); lam_g = lam_top(ng, Pm.alpha_g, n, Pm.LX)
                keep = d2.pois_cdf(lam_f, k) * d2.pois_cdf(lam_g, k)
                if 2 ** S * pr * keep < 4: continue
                surv = rho(max(nf - k * n, 0) / Pm.LX) * rho(max(ng - k * n, 0) / Pm.LX)
                pairs = 2 ** d2.log2li(Pm.LX) / 2
                c2 = pairs * Pm.LLL6 + 2 ** S * (Pm.SIEVE_PP + surv * (Pm.NORM + d2.ecm_instr(n, nf) + d2.ecm_instr(n, ng))) + Pm.LLL12
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
            t = t0 + (1 - t0) * (i + 0.5) / m; s += d2.kappa(u, t) * ((1 - t0) / m) / t
        v = s
    _LAM[key] = v
    return v

if __name__ == "__main__":
    import time
    class P2(d2.P): LMAX = 265
    Pm = P2(); Pm.lB = float(sys.argv[1]) if len(sys.argv) > 1 else 69.5
    t = time.time()
    r = solve_fast(Pm)
    tot, lB1, own, sub, cand, nk = r[1]
    print("fast DP: boot B1=2^%d candidates 2^%.1f own 2^%.1f tree 2^%.1f -> per log 2^%.2f units; both 2^%.2f  (%.1fs)" % (
        lB1, math.log2(cand), math.log2(own), math.log2(sub), r[0] - math.log2(20000), r[0] + 1 - math.log2(20000), time.time() - t))
