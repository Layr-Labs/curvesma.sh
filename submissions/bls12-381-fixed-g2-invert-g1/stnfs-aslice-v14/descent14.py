"""stnfs-aslice-v14 descent: every node, including the boot, is a window-sieved
special-q lattice node (DESCENT.md). This replaces the v13 recurrence: the exact
norm screen at 2^26 instructions per trial is replaced by (i) certified c-space
block bounds enforcing the envelopes, (ii) a bucket sieve over the sieving primes
with per-(prime, root pair) hit-lattice walks, (iii) a threshold on the residual
mass with a bounded child window (k primes above the sieve bound), and (iv) ECM
on survivors' residuals. The no-sieve ECM-per-trial option is retained as the
honest control at every bin.

Envelope constants (CF, CG) are enforced by the certified bounds; the yield model
applies proba() to the envelope and multiplies by the measured retained fraction
(ENV_FACTOR). SIEVE_CELL and PAIR_SETUP are kernel-measured allowances
(KERNELS.md/measurements14.txt).
"""
import math
import filtsim
rho = filtsim.rho
LN2 = math.log(2)
UNIT = 20000.0

# --- GMP-ECM table model (measured ecm12 Montgomery multiplication) ---
ECM_TAB = {20: (11e3, 86), 25: (5e4, 214), 30: (25e4, 430), 35: (1e6, 904), 40: (3e6, 2350), 45: (11e6, 4480),
           50: (43e6, 7553), 55: (11e7, 17769), 60: (26e7, 42017), 65: (85e7, 69408), 70: (29e8, 102212),
           75: (76e8, 188056), 80: (25e9, 265557)}
MONT5 = 367.76
ECM_FIND = 0.99


def ecm_instr(fbits, nbits, success=ECM_FIND):
    d = fbits / 3.3219
    ks = sorted(ECM_TAB)
    if d <= ks[-1]:
        k = min(x for x in ks if x >= d); b1, cur = ECM_TAB[k]
    else:
        b1, cur = ECM_TAB[ks[-1]]; e = (d - ks[-1]) / 5; b1 *= 3.3 ** e; cur *= 1.5 ** e
    curves = cur * math.log(1 / (1 - success))
    L = max(5, math.ceil(nbits / 64)); mulmod = MONT5 * (L / 5.0) ** 2
    return curves * 2 * 1.44 * b1 * 11 * mulmod


def proba(log2n, alpha, log2b):
    u = (log2n + alpha) / log2b
    return rho(u) + (1 - 0.5772) * rho(u - 1) / (log2n * LN2)


def kappa(u, t):
    r = rho(u)
    return rho(u - t) / r if r > 0 else 0.0


def log2li(l):
    x = l * LN2
    return l - math.log2(x) + math.log2(1 + 1 / x + 2 / x / x)


def pois_cdf(lam, k):
    s = 0.0; term = math.exp(-lam)
    for i in range(k + 1):
        s += term; term *= lam / (i + 1)
    return s


_LAM = {}
def lam_top(nbits, alpha, n, LX, m=40):
    key = (round(nbits, 1), n)
    if key in _LAM:
        return _LAM[key]
    u = (nbits + alpha) / n; t0 = LX / n
    if t0 >= 1:
        v = 0.0
    else:
        s = 0.0
        for i in range(m):
            t = t0 + (1 - t0) * (i + 0.5) / m
            s += kappa(u, t) * ((1 - t0) / m) / t
        v = s
    _LAM[key] = v
    return v


def children_bins(nf, ng, n, lB, pf_alpha, pg_alpha, scale=(1.0, 1.0)):
    """expected children per 1-bit bin [b, b+1) for b in [floor(lB), n), per side,
    scaled by the conditional means given the <=k window (sieve nodes)."""
    out = []
    uf = (nf + pf_alpha) / n; ug = (ng + pg_alpha) / n
    b = math.floor(lB)
    while b < n:
        lo, hi = max(b, lB), min(b + 1, n)
        if hi > lo:
            t = 0.5 * (lo + hi) / n; w = (hi - lo) / n / t
            out.append((b, kappa(uf, t) * w * scale[0], kappa(ug, t) * w * scale[1]))
        b += 1
    return out


class P:
    lB = 69.0; alpha_f = 5.0; alpha_g = 1.0
    CF = 96.0; CG = 17.0; log2U6 = 6 * 63.714
    BOOT_NORM = 2800.0; LLL12 = 2.0 ** 36; NORM = 2.0 ** 26
    LX = 69.0; PAD_CHILD = 2.0; PAD_NODE = 4.0
    LMAX = 305; MIN_RELS = 8.0
    ENV_FACTOR = 0.999          # measured retained fraction at (CF, CG)
    SIEVE_CELL = 3268000.0    # = 2^21.64: measured 2946 Ir/leaf x the smallest
                               # measured leaves/hit (254, at p=2^20 where the hit
                               # mass lives) + record ~20, times 4.32 hits/cell,
                               # + 10.9 zero/apply/scan + 1200 c-space certified
                               # bound screen (2x the measured a-space eval at a
                               # doubled eval rate; see KERNELS.md)
    PAIR_SETUP = 2.0 ** 24.8   # measured sound-build 12-dim LLL (11.6M/21.0M/30.1M
                               # Ir at p=2^20/2^40/2^62) + functional per (p, root)
    KS = (2, 3, 4, 6)
    SLACK = 96.0                # threshold residual slack (sub-2^8 mass, powers, rounding)
    STOP_FACTOR = 4.0


def node(l, side, n, Pm, S_grid):
    """Both options' costs are increasing in S (the per-trial/per-cell constants
    are non-decreasing in the norm sizes, which grow with S), and feasibility
    (2^S x pr x keep) is increasing in S up to the small keep drift, so only the
    minimal feasible S per option/window needs evaluation; the S scan ascends and
    records each (option, k) the first time it becomes feasible."""
    c0 = l / 12
    LXe = min(n, Pm.LX)
    pairs = 2.0 * 2 ** log2li(LXe)
    best = None
    done_ecm = False
    done_k = {k: False for k in Pm.KS}
    for S in S_grid:
        c = c0 + S / 12
        nf = 72 * c + Pm.CF; ng = 12 * c + Pm.log2U6 + Pm.CG
        if side == 'f':
            nf -= l
        else:
            ng -= l
        pr = proba(nf, Pm.alpha_f, n) * proba(ng, Pm.alpha_g, n) * Pm.ENV_FACTOR
        if pr <= 0:
            continue
        if 2 ** S * pr < Pm.MIN_RELS:
            continue
        # option A: no sieve, full ECM schedule per trial (control; also the only
        # option when the sieve model is disabled)
        if not done_ecm:
            done_ecm = True
            per = Pm.NORM + ecm_instr(n, nf) + ecm_instr(n, ng)
            c1 = 2 ** S * per + Pm.LLL12
            if best is None or c1 < best[0]:
                best = (c1, S, nf, ng, 1.0, 'ecm', n, None)
        if all(done_k.values()):
            break
        lam_f = lam_top(nf, Pm.alpha_f, n, LXe)
        lam_g = lam_top(ng, Pm.alpha_g, n, LXe)
        for k in Pm.KS:
            if done_k[k]:
                continue
            kf = pois_cdf(lam_f, k); kg = pois_cdf(lam_g, k)
            keep = kf * kg
            if 2 ** S * pr * keep < Pm.MIN_RELS:
                continue
            done_k[k] = True
            surv = (rho(max(nf - k * n - Pm.SLACK, 0) / LXe) *
                    rho(max(ng - k * n - Pm.SLACK, 0) / LXe))
            surv_instr = Pm.NORM + 2 * ecm_instr(n, k * n + Pm.SLACK)
            surv_rate = surv + pr * keep
            need_surv = Pm.STOP_FACTOR * Pm.MIN_RELS * (1.0 + surv / max(pr * keep, 1e-300))
            ecm_cells = min(2 ** S * surv_rate, need_surv)
            c2 = (pairs * Pm.PAIR_SETUP + 2 ** S * Pm.SIEVE_CELL +
                  ecm_cells * surv_instr + Pm.LLL12)
            sf = pois_cdf(lam_f, k - 1) / kf if kf > 0 else 0.0
            sg = pois_cdf(lam_g, k - 1) / kg if kg > 0 else 0.0
            if best is None or c2 < best[0]:
                best = (c2, S, nf, ng, keep, 'sieve%d' % k, n, (sf, sg))
    return best


def solve(Pm, s_step=2):
    bins = list(range(int(math.floor(Pm.lB)), Pm.LMAX + 1, s_step))
    cost = {}; detail = {}
    S_grid = [12 * math.log2(2 * k + 1) for k in (2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128)]
    for b in bins:
        lo_b = max(float(b), Pm.lB)
        l = 0.5 * (lo_b + b + 2.0)
        for side in ('f', 'g'):
            if b + 2.0 <= Pm.lB + 1e-9:
                cost[(b, side)] = 0.0; continue
            best = None
            # target grid: step 2 below bin 150, step 4 above (the coarser upper
            # grid overestimates the optimum by +0.05 bits, measured against the
            # full step-2 grid; the shipped number is the conservative one)
            tstep = 2 if b < 150 else 4
            targets = [Pm.lB] + [float(x) for x in range(bins[0] + 2, b + 1, tstep) if x > Pm.lB]
            for nn in targets:
                if nn >= l:
                    continue
                opt = node(l, side, nn, Pm, S_grid)
                if opt is None:
                    continue
                c, S, nf, ng, keep, how, nt, cscale = opt
                kids = children_bins(nf, ng, nt, Pm.lB, Pm.alpha_f, Pm.alpha_g,
                                      scale=cscale or (1.0, 1.0))
                sub = 0.0; nkids = 0.0
                for (bb, kf, kg) in kids:
                    key = bb - ((bb - bins[0]) % 2)
                    if (key, 'f') not in cost:
                        raise AssertionError('uncosted child bin %s for node bin %s' % (key, b))
                    sub += Pm.PAD_CHILD * (kf * cost[(key, 'f')] + kg * cost[(key, 'g')])
                    nkids += Pm.PAD_CHILD * (kf + kg)
                tot = Pm.PAD_NODE * c + sub
                if best is None or tot < best[0]:
                    best = (tot, nt, c, S, nf, ng, keep, how, nkids)
            cost[(b, side)] = best[0] if best else float('inf')
            if best:
                detail[(b, side)] = best
    # boot: the p^6-lattice as one window-sieved g-side node
    bestb = None
    for lB1 in range(150, Pm.LMAX, 10):
        pr = proba(Pm.BOOT_NORM, 0.0, lB1) * Pm.ENV_FACTOR
        if pr <= 0:
            continue
        LXe = min(lB1, Pm.LX)
        for k in Pm.KS:
            lam = lam_top(Pm.BOOT_NORM, 0.0, lB1, LXe)
            keep = pois_cdf(lam, k)
            for Sb in S_grid:
                if 2 ** Sb * pr * keep < Pm.MIN_RELS:
                    continue
                surv = rho(max(Pm.BOOT_NORM - k * lB1 - Pm.SLACK, 0) / LXe)
                pairs = 2.0 * 2 ** log2li(LXe)
                surv_instr = Pm.NORM + ecm_instr(lB1, k * lB1 + Pm.SLACK)
                surv_rate = surv + pr * keep
                need_surv = Pm.STOP_FACTOR * Pm.MIN_RELS * (1.0 + surv / max(pr * keep, 1e-300))
                ecm_cells = min(2 ** Sb * surv_rate, need_surv)
                own = (pairs * Pm.PAIR_SETUP + 2 ** Sb * Pm.SIEVE_CELL +
                       ecm_cells * surv_instr + Pm.LLL12)
                kids = children_bins(Pm.BOOT_NORM, 1.0, lB1, Pm.lB, 0.0, 0.0)
                sf = pois_cdf(lam, k - 1) / keep if keep > 0 else 0.0
                sub = 0.0
                for (bb, kf, kg) in kids:
                    key = bb - ((bb - bins[0]) % 2)
                    sub += Pm.PAD_CHILD * kf * sf * cost[(key, 'g')]
                tot = Pm.PAD_NODE * own + sub
                if bestb is None or tot < bestb[0]:
                    bestb = (tot, lB1, own, sub, (Sb, k, keep),
                             sum(Pm.PAD_CHILD * kf * sf for (_, kf, _) in kids))
    if bestb is None:
        raise ValueError('no feasible boot')
    one_log = math.log2(bestb[0])
    return one_log, bestb, detail


def solve_units(Pm):
    one_log, bestb, detail = solve(Pm)
    two_log_units = one_log - math.log2(UNIT) + 1.0
    return two_log_units, bestb, detail
