"""two_bound_prod.py - live production-scale check of the corrected two-bound product tree.

The shipped recursion is V20.enum_two_bound in v20_lat_b2.py: it intersects the full 12D
Q38 interval with the two rank-6 weighted forms F_a = 1764*Q(a) <= RA2W and
F_b = 1000*Q(b) <= RB2W (CORRECT b coordinates Bp[i][6+u]), relaxes a form/level whose
alpha < 2^20 (monotone: dropping the non-negative term can only add leaves), and keeps
the package's chord split.  This module measures that product tree at production
parameters with the package's 32 x 4 = 128-draw exact rejection-sampled Haldane-UMVUE
Horvitz-Thompson estimator (6 level-1 subtrees per draw = 768 sampled subtrees), and
separately exercises the native Q38/Q32 recursion over the same 128 draws with the
strict u64/i64 checks.

Exact per-draw C (the conservative bound):
    C = n0*75 + n1*115 + (chords/leaves)*sp_per_chord + (fast*13.35+bnd*53.05)/leaves + 2.85
with n0 = level-0 nodes / product leaf, n1 = reached level-1 nodes / product leaf,
n0/n1 per sampled level-1 subtree.  The HT ratio over the 32 equal-log strata uses the
exact multiplicity m_p and exact 1/p weights, N_hat_s = M_s*(K-1)/(S_T-1), and the
linearized stratified SE including Var(N_hat_s)=N_hat_s^2/(K-2).
"""
import math, random, time
from fractions import Fraction as Fr
from math import sqrt, log2

import model_bn254_v6 as M
import descent_bn254_v6 as DESC
import v20_lat_b2 as V20

WA, WB = 1764, 1000
Q32 = 1 << 32
Q38 = 1 << 38
N_STRATA = 32
M_CH = 6
C0, C1, APPLY = 75.0, 115.0, 2.85
FAST_IR, BND_IR = 13.35, 53.05
THETA = 1 << 20
SEED_ORB = 20261122
SEED_BOX = 0xB0B0B0
# regression pins: the corrected 128-draw x6 exact estimator (2026-11-22)
PIN_MEAN = 36.39157215192703
PIN_SE = 1.294421324205492
PIN_P3SE = 40.274836124543505
PIN_103 = 41.48308120827981


def _geom():
    D = M.load_data(); T4 = D["t4_geom"]; PQ = D["pool_quad"]
    global _S, _S2
    G_sym, D_mats, UMATS, S, S2, G_SC = M.build_exact_int_geom(T4)
    _S, _S2 = S, S2
    A = int(round(WA * (PQ["ra_op"] ** 2) * G_SC))
    B = int(round(WB * (PQ["rb_ga_op"] ** 2) * G_SC))
    T = int(PQ["t12_box"])
    assert A + B == T
    return dict(G_sym=G_sym, S=S, S2=S2, G_SC=G_SC, A=A, B=B, T=T,
                p_lat=M.OP["p_lat"], lx2=M.OP["lx2"])


_S = _S2 = None


def _in_swar_box(v):
    a0, b0 = v[:6], v[6:]
    A3M = M.A_TAB_SLOTS_A3; B0M = M.B0_ARENA_MAX
    for ta, tb in ((a0, b0), (V20._mv6(_S, a0), V20._mv6(_S, b0)), (V20._mv6(_S2, a0), V20._mv6(_S2, b0))):
        if not (-510 <= ta[0] <= 510 and -510 <= ta[1] <= 510 and -510 <= ta[2] <= 510 and -A3M < ta[3] < A3M and -510 <= ta[4] <= 510 and -510 <= ta[5] <= 510):
            return False
        if not (-B0M < tb[0] < B0M and -1022 <= tb[1] <= 1022 and -510 <= tb[2] <= 510 and -510 <= tb[3] <= 510 and -510 <= tb[4] <= 510 and -510 <= tb[5] <= 510):
            return False
    return True


def _box_interval(base, B0, l0, h0):
    lo = Fr(l0); hi = Fr(h0)
    A3M = M.A_TAB_SLOTS_A3; B0M = M.B0_ARENA_MAX
    def apply(c, slope, bl, br, slo=False, shi=False):
        nonlocal lo, hi
        if slope == 0:
            return (bl < c < br) if (slo and shi) else (bl <= c <= br)
        x1 = Fr(bl - c, slope); x2 = Fr(br - c, slope)
        if x1 > x2:
            x1, x2 = x2, x1
        lx = (math.floor(x1) + 1) if slo else -((-x1.numerator) // x1.denominator)
        hx = (math.ceil(x2) - 1) if shi else (x2.numerator // x2.denominator)
        if lx > lo:
            lo = Fr(lx)
        if hx < hi:
            hi = Fr(hx)
        return True
    for vbase, vslope in ((base[:6], B0[:6]), (V20._mv6(_S, base[:6]), V20._mv6(_S, B0[:6])), (V20._mv6(_S2, base[:6]), V20._mv6(_S2, B0[:6]))):
        for c in range(6):
            if c == 3:
                if not apply(vbase[c], vslope[c], -A3M, A3M, True, True):
                    return None
            else:
                if not apply(vbase[c], vslope[c], -510, 510):
                    return None
    for vbase, vslope in ((base[6:], B0[6:]), (V20._mv6(_S, base[6:]), V20._mv6(_S, B0[6:])), (V20._mv6(_S2, base[6:]), V20._mv6(_S2, B0[6:]))):
        for c in range(6):
            if c == 0:
                if not apply(vbase[c], vslope[c], -B0M, B0M, True, True):
                    return None
            elif c == 1:
                if not apply(vbase[c], vslope[c], -1022, 1022):
                    return None
            else:
                if not apply(vbase[c], vslope[c], -510, 510):
                    return None
    return max(l0, math.floor(lo)), min(h0, math.floor(hi))


def _fstats(base, B0, l0, h0, clamp=False):
    clips0 = 0
    if clamp:
        bi = _box_interval(base, B0, l0, h0)
        if bi is None:
            return dict(lvs=0, fast=0, bnd=0, segs=0, bspl=0, fspl=0, cross=0, clips=0)
        bl, bh = bi
        if bl > l0:
            clips0 += 1; l0 = bl
        if bh < h0:
            clips0 += 1; h0 = bh
    clips = clips0
    if not clamp:
        while l0 <= h0 and not _in_swar_box([base[c] + l0 * B0[c] for c in range(12)]):
            clips += 1; l0 += 1
        while l0 <= h0 and not _in_swar_box([base[c] + h0 * B0[c] for c in range(12)]):
            clips += 1; h0 -= 1
    if l0 > h0:
        return dict(lvs=0, fast=0, bnd=0, segs=0, bspl=0, fspl=0, cross=0, clips=clips)
    length = h0 - l0 + 1
    vl = [base[c] + l0 * B0[c] for c in range(12)]
    vh = [base[c] + h0 * B0[c] for c in range(12)]
    clo, chi = V20._cone_of(vl), V20._cone_of(vh)
    if clo >= 0 and clo == chi:
        return dict(lvs=length, fast=length, bnd=0, segs=1, bspl=0, fspl=0, cross=0, clips=clips)
    db0_k = (B0[6], B0[6] + B0[8] + B0[10], B0[6] + B0[7] + B0[10])
    vlb0_k = (vl[6], vl[6] + vl[8] + vl[10], vl[6] + vl[7] + vl[10])
    if clo >= 0 and chi >= 0 and (clo >> 1) == (chi >> 1):
        k = clo >> 2
        tf = V20._floor_div(-vlb0_k[k], db0_k[k]); tc = V20._ceil_div(-vlb0_k[k], db0_k[k])
        if tf == tc:
            return dict(lvs=length, fast=length - 1, bnd=1,
                        segs=(1 if tf > 0 else 0) + (1 if tf + 1 < length else 0),
                        bspl=1, fspl=0, cross=1, clips=clips)
        return dict(lvs=length, fast=length, bnd=0, segs=2, bspl=1, fspl=0, cross=1, clips=clips)
    umax = length - 1
    a3, a4, a5 = vl[3], vl[4], vl[5]; d3, d4, d5 = B0[3], B0[4], B0[5]
    bps = {0, umax}; cross = 0
    for c, m in ((a3 - a4, d3 - d4), (a3 + a4, d3 + d4), (a3 - a5, d3 - d5), (a3 + a5, d3 + d5), (a4 - a5, d4 - d5), (a4 + a5, d4 + d5)):
        if m != 0:
            v0, v1 = c, c + umax * m
            if v0 == 0 or v1 == 0 or ((v0 ^ v1) < 0):
                cross += 1
                for rx in (V20._floor_div(-c, m), V20._ceil_div(-c, m)):
                    if 0 <= rx <= umax:
                        bps.add(rx)
    sbps = sorted(bps); k_mask = 0
    for u in sbps:
        x3, x4, x5 = abs(a3 + u * d3), abs(a4 + u * d4), abs(a5 + u * d5)
        if x3 >= x4 and x3 >= x5:
            k_mask |= 1
        if x4 >= x3 and x4 >= x5:
            k_mask |= 2
        if x5 >= x3 and x5 >= x4:
            k_mask |= 4
    for ki in range(3):
        if (k_mask & (1 << ki)) and db0_k[ki] != 0:
            c, m = vlb0_k[ki], db0_k[ki]
            v0, v1 = c, c + umax * m
            if v0 == 0 or v1 == 0 or ((v0 ^ v1) < 0):
                cross += 1
                for rx in (V20._floor_div(-c, m), V20._ceil_div(-c, m)):
                    if 0 <= rx <= umax:
                        bps.add(rx)
    sbps = sorted(bps); cur_c = -2; cur_R = -1; fast = bnd = segs = 0
    for idx, u in enumerate(sbps):
        cu = V20._cone_of([vl[c] + u * B0[c] for c in range(12)])
        if cu >= 0:
            if cur_c == cu and cur_R + 1 == u:
                cur_R = u
            else:
                if cur_c >= 0:
                    segs += 1
                cur_c, cur_R = cu, u
            fast += 1
        else:
            if cur_c >= 0:
                segs += 1; cur_c = -2
            bnd += 1
        if idx + 1 < len(sbps) and sbps[idx] + 1 <= sbps[idx + 1] - 1:
            uA, uB = sbps[idx] + 1, sbps[idx + 1] - 1
            cA = V20._cone_of([vl[c] + uA * B0[c] for c in range(12)])
            if cA >= 0:
                if cur_c == cA and cur_R + 1 == uA:
                    cur_R = uB
                else:
                    if cur_c >= 0:
                        segs += 1
                    cur_c, cur_R = cA, uB
                fast += (uB - uA + 1)
            else:
                if cur_c >= 0:
                    segs += 1; cur_c = -2
                bnd += (uB - uA + 1)
    if cur_c >= 0:
        segs += 1
    return dict(lvs=length, fast=fast, bnd=bnd, segs=segs, bspl=0, fspl=1, cross=cross, clips=clips)


def _splus(st):
    return st["segs"] * 12.0 + st["bspl"] * 12.0 + st["fspl"] * 22.0 + st["cross"] * 10.0 + st["clips"] * 60.0 + 16.0


def _mat_gram_full(B, G_sym):
    G = [[0] * 12 for _ in range(12)]
    for i in range(12):
        for j in range(12):
            G[i][j] = WA * sum(B[i][u] * G_sym[u][v] * B[j][v] for u in range(6) for v in range(6)) \
                + WB * sum(B[i][6 + u] * G_sym[u][v] * B[j][6 + v] for u in range(6) for v in range(6))
    return G


def _ldl12(G):
    n = 12; Dv = [1] + [0] * n; lam = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            u = G[i][j]
            for k in range(j):
                u = (Dv[k + 1] * u - lam[i][k] * lam[j][k]) // Dv[k]
            if j < i:
                lam[i][j] = u
            else:
                Dv[i + 1] = u
    alpha = [Fr(Dv[k + 1], Dv[k]) for k in range(n)]
    mu = [[Fr(lam[j][k], Dv[k + 1]) if j > k else Fr(0) for k in range(n)] for j in range(n)]
    return alpha, mu


def _load_orbits(g, K_DRAW, seed=SEED_ORB):
    rng_orb = random.Random(seed); orbits = []
    G_sym = g["G_sym"]
    for s_idx in range(N_STRATA):
        lo_b = g["p_lat"] + (g["lx2"] - g["p_lat"]) * s_idx / N_STRATA
        hi_b = g["p_lat"] + (g["lx2"] - g["p_lat"]) * (s_idx + 1) / N_STRATA
        L_i = int(2.0 ** lo_b) | 1; H_i = int(2.0 ** hi_b)
        if (H_i & 1) == 0:
            H_i -= 1
        M_s = (H_i - L_i) // 2 + 1; S_T = 0; draws = []
        for _ in range(K_DRAW):
            tries = 0
            while True:
                tries += 1
                pp4 = L_i + 2 * rng_orb.randrange(M_s)
                if (pp4 & 7) not in (1, 3) or not M._isp(pp4):
                    continue
                xp4 = M._ppow_h([0, 1, 0, 0, 0, 0], pp4, pp4); xp4[1] = (xp4[1] - 1) % pp4
                gp4 = M.H_COEFFS if not any(xp4) else M._pgcd(xp4, M.H_COEFFS, pp4)
                if len(gp4) <= 1:
                    continue
                rs4 = M._roots_mod(gp4, pp4)
                if not rs4:
                    continue
                vp4 = []; u_mod4 = M.U_BN254 % pp4
                for sr4 in rs4:
                    ws4 = V20._eval_poly(M.W_COEFFS, sr4, pp4)
                    if ws4 != 0:
                        rh_pos = V20._sqrt_mod(pow((u_mod4 * ws4) % pp4, -1, pp4), pp4)
                        if rh_pos > 0:
                            vp4.append((sr4, rh_pos))
                if vp4:
                    m_p4 = len(vp4) // 3
                    s04, rh_pos = rng_orb.choice(vp4)
                    rh4 = rh_pos if rng_orb.randrange(2) == 0 else pp4 - rh_pos
                    break
            S_T += tries
            i34 = pow(3, -1, pp4); sp4 = [1]
            for _ in range(5):
                sp4.append((sp4[-1] * s04) % pp4)
            es4 = [(sum(M.ZK_NUM[k][j] * sp4[j] for j in range(6)) * i34) % pp4 for k in range(6)]
            lf4 = es4 + [((-rh4) * e) % pp4 for e in es4]; i04 = pow(lf4[0], -1, pp4)
            B4 = [[pp4 if c == 0 else 0 for c in range(12)]]
            for r1 in range(1, 12):
                c0 = (-lf4[r1] * i04) % pp4
                if c0 > pp4 // 2:
                    c0 -= pp4
                B4.append([c0 if c == 0 else (1 if c == r1 else 0) for c in range(12)])
            _, D4, lm4 = DESC.bareiss_lll_12d(B4, G_sym, WA, WB, False, 9, 10, (99, 100), True)
            draws.append(dict(pp4=pp4, m_p=m_p4, s0=s04, rho=rh4, B4=B4, D4=D4, lm4=lm4, tries=tries))
        orbits.append(dict(M_s=M_s, S_T=S_T, draws=draws))
    return orbits


def _prep(d, g):
    B4 = d["B4"]; D4 = d["D4"]; lm4 = d["lm4"]; T = g["T"]
    a12, m12 = _ldl12(_mat_gram_full(B4, g["G_sym"]))
    GA = V20._tb_gram(B4, g["G_sym"], 0, WA); GB = V20._tb_gram(B4, g["G_sym"], 6, WB)
    aa, ma = V20._tb_ldl6(GA); ab, mb = V20._tb_ldl6(GB)
    Xm = [0] * 12
    for k in range(11, -1, -1):
        rk = math.isqrt(((D4[k] * T) + D4[k + 1] - 1) // D4[k + 1]) + 1
        Xm[k] = rk + sum((Xm[j] + 1) // 2 for j in range(k + 1, 12))
    ec = [sum(Xm[j] for j in range(k + 1, 12)) + 1 for k in range(12)]
    b2 = [min((1 << 64) - 1, (D4[k + 1] << 64) // (D4[k] * T)) for k in range(12)]
    inv = [((D4[k] * T << 38) + D4[k + 1] - 1) // D4[k + 1] for k in range(12)]
    mu32 = [[((lm4[j][k] << 32) // D4[k + 1]) if j > k else 0 for k in range(12)] for j in range(12)]
    Bn = [D4[k + 1] / D4[k] for k in range(12)]
    ctx = V20._tb_setup(B4, g["G_sym"], WA, WB, g["A"], g["B"], T, mu32, b2, inv, ec, prod=True)
    d.update(a12=a12, m12=m12, aa=aa, ma=ma, ab=ab, mb=mb, ec=ec, b2=b2, inv=inv,
             mu32=mu32, Bn=Bn, ctx=ctx, GA=GA, GB=GB)
    return d


def _sample_l1(d, g, rng_box):
    D4 = d["D4"]; mu32 = d["mu32"]; ec = d["ec"]; b2 = d["b2"]; inv = d["inv"]
    Bn = d["Bn"]; R_f = sqrt(float(g["T"])); inv_sB = [1.0 / sqrt(Bn[k]) for k in range(12)]
    mu_c = [[mu32[j][k] for j in range(k + 1, 12)] for k in range(12)]
    gv = [rng_box.gauss(0.0, 1.0) for _ in range(11)]
    rg = R_f * (rng_box.random() ** (1.0 / 10.0)) / sqrt(math.sumprod(gv, gv))
    xp = [0] * 12; rm = Q38
    for k in range(11, 1, -1):
        ck = -math.sumprod(mu_c[k], xp[k + 1:]) if k < 11 else 0
        if abs(ck) > 10 ** 12:
            return None
        xk = int(round(ck * (1.0 / Q32) + gv[k - 1] * rg * inv_sB[k])); xp[k] = xk
        if abs(xk) > 10 ** 9:
            return None
        ad = abs((xk << 32) - ck)
        if ad > ec[k]:
            d19 = (ad - ec[k]) >> 13; sub = ((d19 * d19) * b2[k]) >> 64
            if sub > rm:
                return None
            rm -= sub
    return xp, rm


def _prod_interval(base, B0, g, l0e, h0e):
    G_sym = g["G_sym"]; A = g["A"]; B = g["B"]
    a0 = base[:6]; b0 = base[6:]; u = B0[:6]; w = B0[6:]
    Aa = V20._qf6(G_sym, u); Ab = V20._qf6(G_sym, w)
    Ba = 2 * sum(a0[i] * G_sym[i][j] * u[j] for i in range(6) for j in range(6))
    Bb = 2 * sum(b0[i] * G_sym[i][j] * w[j] for i in range(6) for j in range(6))
    Ca = V20._qf6(G_sym, a0); Cb = V20._qf6(G_sym, b0)
    def iv(Ac, Bv, Cv, T):
        if Ac <= 0:
            return None
        disc = Bv * Bv - 4 * Ac * (Cv - T)
        if disc < 0:
            return None
        s = math.isqrt(disc)
        lo = math.ceil(Fr(-Bv - s, 2 * Ac)); hi = math.floor(Fr(-Bv + s, 2 * Ac))
        while lo <= hi and Ac * lo * lo + Bv * lo + Cv > T:
            lo += 1
        while Ac * (lo - 1) * (lo - 1) + Bv * (lo - 1) + Cv <= T:
            lo -= 1
        while lo <= hi and Ac * hi * hi + Bv * hi + Cv > T:
            hi -= 1
        while Ac * (hi + 1) * (hi + 1) + Bv * (hi + 1) + Cv <= T:
            hi += 1
        return lo, hi
    ia = iv(WA * Aa, WA * Ba, WA * Ca, A); ib = iv(WB * Ab, WB * Bb, WB * Cb, B)
    if ia is None or ib is None:
        return None
    lo = max(ia[0], ib[0], l0e); hi = min(ia[1], ib[1], h0e)
    return (lo, hi) if lo <= hi else None


def _ell_x0_interval(rem0, c0, g):
    """Exact full-form level-0 interval, as in the package prototype."""
    from math import floor
    alpha = None  # filled by caller via closure; kept for API clarity
    return None


def _exact_subtree(d, xp, g):
    a12 = d["a12"]; m12 = d["m12"]; aa = d["aa"]; ma = d["ma"]; ab = d["ab"]; mb = d["mb"]
    T = g["T"]; A = g["A"]; B = g["B"]
    x = [0] * 12
    for k in range(2, 12):
        x[k] = xp[k]
    S = Fr(0)
    for k in range(2, 12):
        c = sum((m12[j][k] * x[j] for j in range(k + 1, 12)), Fr(0))
        S += a12[k] * (x[k] + c) ** 2
    rem = T - S
    if rem < 0:
        return None
    Sa = Fr(0); Sb = Fr(0)
    for k in range(2, 6):
        ca = sum((ma[j][k] * x[j] for j in range(k + 1, 12)), Fr(0)); Sa += aa[k] * (x[k] + ca) ** 2
        cb = sum((mb[j][k] * x[j] for j in range(k + 1, 12)), Fr(0)); Sb += ab[k] * (x[k] + cb) ** 2
    rema = A - Sa; remb = B - Sb
    if rema < 0 or remb < 0:
        return None
    def root_iv(rem, alpha, c):
        if rem < 0:
            return 0, -1
        m = math.floor(-c)
        if alpha * (m + c) ** 2 <= rem:
            seed = m
        elif alpha * (m + 1 + c) ** 2 <= rem:
            seed = m + 1
        else:
            return 0, -1
        lo = seed
        while alpha * (lo - 1 + c) ** 2 <= rem:
            lo -= 1
        hi = seed
        while alpha * (hi + 1 + c) ** 2 <= rem:
            hi += 1
        return lo, hi
    c1 = sum((m12[j][1] * x[j] for j in range(2, 12)), Fr(0))
    lo1, hi1 = root_iv(rem, a12[1], c1)
    ca1 = sum((ma[j][1] * x[j] for j in range(2, 12)), Fr(0))
    la1, ha1 = root_iv(rema, aa[1], ca1)
    cb1 = sum((mb[j][1] * x[j] for j in range(2, 12)), Fr(0))
    lb1, hb1 = root_iv(remb, ab[1], cb1)
    lo = max(lo1, la1, lb1); hi = min(hi1, ha1, hb1)
    if lo > hi:
        return None
    enodes = eleaves = 0
    for x1 in range(lo1, hi1 + 1):
        t1 = x1 + c1; rem0e = rem - a12[1] * t1 * t1
        if rem0e < 0:
            continue
        c0e = sum((m12[j][0] * x[j] for j in range(1, 12)), Fr(0))
        l0e, h0e = root_iv(rem0e, a12[0], c0e)
        if l0e <= h0e:
            enodes += 1; eleaves += (h0e - l0e + 1)
    nodes0 = leaves = chords = 0; fast = bnd = 0; sp = 0.0
    for x1 in range(lo, hi + 1):
        x[1] = x1
        t1 = x1 + c1; rem0 = rem - a12[1] * t1 * t1
        if rem0 < 0:
            continue
        ra0 = rema - aa[1] * (x1 + ca1) ** 2
        rb0 = remb - ab[1] * (x1 + cb1) ** 2
        if ra0 < 0 or rb0 < 0:
            continue
        nodes0 += 1
        base = [sum(x[j] * d["B4"][j][c] for j in range(1, 12)) for c in range(12)]
        pv = _prod_interval(base, d["B4"][0], g, -10 ** 12, 10 ** 12)
        if pv is None:
            continue
        c0 = sum((m12[j][0] * x[j] for j in range(1, 12)), Fr(0))
        l0f, h0f = root_iv(rem0, a12[0], c0)
        l = max(pv[0], l0f); h = min(pv[1], h0f)
        if l > h:
            continue
        chords += 1; leaves += (h - l + 1)
        st = _fstats(base, d["B4"][0], l, h, clamp=True)
        fast += st["fast"]; bnd += st["bnd"]; sp += _splus(st)
    x[1] = 0
    return nodes0, leaves, chords, enodes, eleaves, fast, bnd, sp


def _fp_subtree(d, xp, rm):
    ctx = d["ctx"]; B4 = d["B4"]; a12 = d["a12"]; m12 = d["m12"]
    T = d["_T"]
    x = [0] * 12
    for k in range(2, 12):
        x[k] = xp[k]
    S = Fr(0)
    for k in range(2, 12):
        c = sum((m12[j][k] * x[j] for j in range(k + 1, 12)), Fr(0))
        S += a12[k] * (x[k] + c) ** 2
    rem = T - S
    if rem < 0:
        return None
    rema = remb = Q38
    for k in (5, 4, 3, 2):
        if ctx["oka"][k]:
            ck = V20._tb_ck(ctx["mua"], x, k)
            rema = V20._tb_sub(rema, ck, ctx["ecab"][k], ctx["b2a"], k, x[k])
            if rema is None:
                return None
        if ctx["okb"][k]:
            ck = V20._tb_ck(ctx["mub"], x, k)
            remb = V20._tb_sub(remb, ck, ctx["ecab"][k], ctx["b2b"], k, x[k])
            if remb is None:
                return None
    def root_iv(rem, alpha, c):
        if rem < 0:
            return 0, -1
        m = math.floor(-c)
        if alpha * (m + c) ** 2 <= rem:
            seed = m
        elif alpha * (m + 1 + c) ** 2 <= rem:
            seed = m + 1
        else:
            return 0, -1
        lo = seed
        while alpha * (lo - 1 + c) ** 2 <= rem:
            lo -= 1
        hi = seed
        while alpha * (hi + 1 + c) ** 2 <= rem:
            hi += 1
        return lo, hi
    c1 = sum((m12[j][1] * x[j] for j in range(2, 12)), Fr(0))
    lo, hi = root_iv(rem, a12[1], c1)
    if ctx["oka"][1]:
        la, ha = V20._tb_iv(rema, ctx["inva"], V20._tb_ck(ctx["mua"], x, 1), ctx["ecab"][1], 1)
        lo = max(lo, la); hi = min(hi, ha)
    if ctx["okb"][1]:
        lb, hb = V20._tb_iv(remb, ctx["invb"], V20._tb_ck(ctx["mub"], x, 1), ctx["ecab"][1], 1)
        lo = max(lo, lb); hi = min(hi, hb)
    if lo > hi:
        return None
    nodes0 = leaves = chords = 0; fast = bnd = 0; sp = 0.0
    for x1 in range(lo, hi + 1):
        x[1] = x1
        t1 = x1 + c1; rem0 = rem - a12[1] * t1 * t1
        if rem0 < 0:
            continue
        ra0 = rema; rb0 = remb
        if ctx["oka"][1]:
            ra0 = V20._tb_sub(rema, V20._tb_ck(ctx["mua"], x, 1), ctx["ecab"][1], ctx["b2a"], 1, x1)
            if ra0 is None:
                continue
        if ctx["okb"][1]:
            rb0 = V20._tb_sub(remb, V20._tb_ck(ctx["mub"], x, 1), ctx["ecab"][1], ctx["b2b"], 1, x1)
            if rb0 is None:
                continue
        c0 = sum((m12[j][0] * x[j] for j in range(1, 12)), Fr(0))
        l, h = root_iv(rem0, a12[0], c0)
        if l > hi + 10 ** 12:
            continue
        base = [sum(x[j] * B4[j][c] for j in range(1, 12)) for c in range(12)]
        pv = _prod_interval(base, B4[0], d["_g"], -10 ** 12, 10 ** 12)
        if pv is None:
            continue
        l = max(l, pv[0]); h = min(h, pv[1])
        if ctx["oka"][0]:
            la, ha = V20._tb_iv(ra0, ctx["inva"], V20._tb_ck(ctx["mua"], x, 0), ctx["ecab"][0], 0)
            l = max(l, la); h = min(h, ha)
        if ctx["okb"][0]:
            lb, hb = V20._tb_iv(rb0, ctx["invb"], V20._tb_ck(ctx["mub"], x, 0), ctx["ecab"][0], 0)
            l = max(l, lb); h = min(h, hb)
        if l > h:
            continue
        nodes0 += 1
        st = _fstats(base, B4[0], l, h, clamp=True)
        if st["lvs"] == 0:
            continue
        chords += 1; leaves += st["lvs"]
        fast += st["fast"]; bnd += st["bnd"]; sp += _splus(st)
    x[1] = 0
    return nodes0, leaves, chords, 0, 0, fast, bnd, sp


def _measure(d, g, rng_box, mode, m_ch=M_CH):
    reach = nodes0 = leaves = chords = enodes = eleaves = fast = bnd = 0; sp = 0.0
    relaxed = []
    for _ in range(m_ch):
        s = _sample_l1(d, g, rng_box)
        if s is None:
            continue
        xp, rm = s
        r = _exact_subtree(d, xp, g) if mode == "exact" else _fp_subtree(d, xp, rm)
        if r is None:
            continue
        reach += 1; nodes0 += r[0]; leaves += r[1]; chords += r[2]
        enodes += r[3]; eleaves += r[4]; fast += r[5]; bnd += r[6]; sp += r[7]
    relaxed = list(d["ctx"]["relaxed"]) if mode != "exact" else []
    d["ctx"]["relaxed"] = []
    if leaves == 0:
        return None
    leaf = (fast * FAST_IR + bnd * BND_IR) / leaves
    spc = sp / chords if chords else 0.0
    chl = chords / float(leaves)
    n0 = nodes0 / float(leaves); n1 = reach / float(leaves)
    return dict(reach=reach, nodes0=nodes0, leaves=leaves, chords=chords,
                enodes=enodes, eleaves=eleaves, fast=fast, bnd=bnd, sp=sp,
                leaf=leaf, spc=spc, chl=chl, n0=n0, n1=n1,
                ne0=(enodes / float(eleaves)) if eleaves else None,
                ir=n0 * C0 + n1 * C1 + chl * spc + leaf + APPLY,
                relaxed=relaxed)


def _ht(orbits, key, K_DRAW):
    Y = X = 0.0
    for o in orbits:
        N_hat = o["M_s"] * (K_DRAW - 1.0) / (o["S_T"] - 1.0); sy = sx = 0.0
        for d in o["draws"]:
            v = key(d)
            if v is None:
                continue
            sy += d["m_p"] * v / d["pp4"]; sx += d["m_p"] / d["pp4"]
        Y += N_hat * sy / K_DRAW; X += N_hat * sx / K_DRAW
    return Y / X if X else None


def _htse(orbits, key, mean, K_DRAW):
    X = 0.0; var = 0.0
    for o in orbits:
        N_hat = o["M_s"] * (K_DRAW - 1.0) / (o["S_T"] - 1.0); sx = 0.0; us = []
        for d in o["draws"]:
            v = key(d)
            if v is None:
                continue
            w = d["m_p"] / d["pp4"]; sx += w; us.append(w * (v - mean))
        X += N_hat * sx / K_DRAW
        ub = sum(us) / K_DRAW; s2 = sum((u - ub) ** 2 for u in us) / (K_DRAW - 1)
        var_N = (N_hat * N_hat) / float(K_DRAW - 2)
        var += (N_hat * N_hat + var_N) * (s2 / K_DRAW) + var_N * (ub * ub)
    return sqrt(var) / X



WITNESS_ORBITS = (
    (94113356107938059, 92020840183298046, 22911412904947162, 1, "a", 5),
    (641242830317571683, 446788972188432864, 377812225690998158, 1, "a", 5),
    (1738579582004560859, 522942446165567004, 1060315281056501781, 1, "b", 5),
)


def _orbit_draw(pp4, s0, rho, m_p, g):
    i3 = pow(3, -1, pp4); sp = [1]
    for _ in range(5):
        sp.append((sp[-1] * s0) % pp4)
    es = [(sum(M.ZK_NUM[k][j] * sp[j] for j in range(6)) * i3) % pp4 for k in range(6)]
    lf = es + [((-rho) * e) % pp4 for e in es]; i0 = pow(lf[0], -1, pp4)
    B4 = [[pp4 if c == 0 else 0 for c in range(12)]]
    for r1 in range(1, 12):
        c0 = (-lf[r1] * i0) % pp4
        if c0 > pp4 // 2:
            c0 -= pp4
        B4.append([c0 if c == 0 else (1 if c == r1 else 0) for c in range(12)])
    _, D4, lm4 = DESC.bareiss_lll_12d(B4, g["G_sym"], WA, WB, False, 9, 10, (99, 100), True)
    return dict(pp4=pp4, m_p=m_p, s0=s0, rho=rho, B4=B4, D4=D4, lm4=lm4, tries=1)


EXTRA_CAL_DRAWS = 19
SEED_CAL = 0xC0FFEE


def alpha_scan_calibration(orbits, g):
    """Live alpha calibration: the estimator's 128 draws, 19 extra production-prime draws
    (fixed calibration seed), and the 3 relaxed draws found by the full 1,024-draw scan
    (all at k=5).  Each relaxed draw is rebuilt from its disclosed (p, s0, rho), reduced
    live, run guarded, then forced onto the fast path to show the u64 check fails."""
    extra = [d for o in _load_orbits(g, 1, SEED_CAL) for d in o["draws"]][:EXTRA_CAL_DRAWS]
    for d in extra:
        d["_g"] = g; d["_T"] = g["T"]
    scan = _alpha_scan(orbits + [dict(M_s=0, S_T=0, draws=extra)])
    wits = []
    for pp4, s0, rho, m_p, form, k in WITNESS_ORBITS:
        d = _orbit_draw(pp4, s0, rho, m_p, g); _prep(d, g); d["_T"] = g["T"]; d["_g"] = g
        aa, ab = d["aa"], d["ab"]
        rel = [(f, kk) for f, al in (("a", aa), ("b", ab)) for kk in range(6) if al[kk] < THETA]
        rec = dict(pp4=pp4, relaxed=rel, guarded_ok=False)
        rng = random.Random(SEED_BOX)
        for _ in range(8):
            smp = _sample_l1(d, g, rng)
            if smp is None:
                continue
            d["ctx"]["relaxed"] = []
            try:
                r = _fp_subtree(d, smp[0], smp[1])
            except OverflowError as ex:
                rec["guarded_error"] = str(ex); r = None
            if r is not None:
                rec["guarded_ok"] = True; break
        ctx = d["ctx"]; saved = {}
        for f, kk in rel:
            key = "oka" if f == "a" else "okb"
            saved[key] = list(ctx[key]); ctx[key][kk] = True
        for f, kk in rel:
            inv = ctx["inva"] if f == "a" else ctx["invb"]
            try:
                V20._tb_iv(Q38, inv, 0, 0, kk)
            except OverflowError as ex:
                rec["unguarded_interval_probe"] = str(ex)
        for key, val in saved.items():
            ctx[key] = val
        wits.append(rec)
    est_draws = sum(len(o["draws"]) for o in orbits)
    total = est_draws + len(extra) + len(wits)
    n_rel = sum(len(w["relaxed"]) for w in wits)
    return dict(estimator_draws=est_draws, extra_cal_draws=len(extra), calibration_draws=total,
                relaxed=n_rel, witnesses=wits,
                fastpath_frac_k5=1.0 - n_rel / float(2 * total),
                fastpath_frac_all=1.0 - n_rel / float(6 * total),
                scan=scan,
                note="full 1,024-draw scan found 3 relaxed draws, all k=5; this live subset re-runs those 3 plus the estimator draws")


def alpha_scan_full():
    """Static alpha scan over the full 32 x 32 = 1,024-draw production sample.
    Returns the per-level fast-path coverage and live witnesses for every relaxed draw:
    with the alpha guard the Q38/Q32 subtree runs; with the guard forced on (the old
    unconditional plan) the native u64 operand check raises, so the audit fails closed."""
    g = _geom(); orbits = _load_orbits(g, 32)
    for o in orbits:
        for d in o["draws"]:
            d["_g"] = g; d["_T"] = g["T"]
    scan = _alpha_scan(orbits)
    witnesses = []
    for o in orbits:
        for d in o["draws"]:
            rel = d.get("rel") or []
            if not rel:
                continue
            _prep(d, g); d["_T"] = g["T"]; d["_g"] = g
            rec = dict(pp4=d["pp4"], relaxed=[list(x) for x in rel])
            rng = random.Random(SEED_BOX)
            ok_guarded = False
            for _ in range(8):
                smp = _sample_l1(d, g, rng)
                if smp is None:
                    continue
                d["ctx"]["relaxed"] = []
                try:
                    r = _fp_subtree(d, smp[0], smp[1])
                except OverflowError as ex:
                    rec["guarded_error"] = str(ex); r = None
                if r is not None:
                    ok_guarded = True
                    break
            rec["guarded_ok"] = ok_guarded
            # negative control: force the relaxed level onto the fast path
            ctx = d["ctx"]; saved = {}
            for form, k in rel:
                key = "oka" if form == "a" else "okb"
                saved[key] = list(ctx[key]); ctx[key][k] = True
            forced = None
            try:
                rng2 = random.Random(SEED_BOX)
                for _ in range(8):
                    smp = _sample_l1(d, g, rng2)
                    if smp is None:
                        continue
                    _fp_subtree(d, smp[0], smp[1])
                forced = "no_overflow"
            except OverflowError as ex:
                forced = str(ex)
            finally:
                for key, val in saved.items():
                    ctx[key] = val
            rec["unguarded"] = forced
            probe = None
            for form, k in rel:
                inv = ctx["inva"] if form == "a" else ctx["invb"]
                try:
                    V20._tb_iv(Q38, inv, 0, 0, k)
                except OverflowError as ex:
                    probe = str(ex)
            rec["unguarded_interval_probe"] = probe
            witnesses.append(rec)
    return dict(scan=scan, witnesses=witnesses,
                fastpath_frac_k5=scan["fastpath_frac_k5"],
                fastpath_frac_all=scan["fastpath_frac_all"],
                relaxed=scan["relaxed_draws"], elapsed_s=None)

def _alpha_scan(orbits):
    per_k = {k: dict(total=0, relaxed=0) for k in range(6)}
    min_a = min_b = None; relaxed_draws = 0; relaxed_k5 = 0; recs = []
    for o in orbits:
        for d in o["draws"]:
            if "aa" not in d:
                Gs = d["_g"]["G_sym"]
                d["aa"], _ = V20._tb_ldl6(V20._tb_gram(d["B4"], Gs, 0, WA))
                d["ab"], _ = V20._tb_ldl6(V20._tb_gram(d["B4"], Gs, 6, WB))
            aa = d["aa"]; ab = d["ab"]; rel = []
            for k in range(6):
                for form, al in (("a", aa), ("b", ab)):
                    per_k[k]["total"] += 1
                    if al[k] < THETA:
                        per_k[k]["relaxed"] += 1; rel.append((form, k))
                if aa[k] < THETA or ab[k] < THETA:
                    pass
            m = min(min(aa), min(ab))
            if min_a is None or m < min_a:
                min_a = m
            d["rel"] = rel
            if rel:
                relaxed_draws += 1
                if ("a", 5) in rel or ("b", 5) in rel:
                    relaxed_k5 += 1
                recs.append(dict(pp4=d["pp4"], relaxed=rel, min_bits=round(log2(float(m)), 3)))
    return dict(per_k=per_k, min_alpha_bits=round(log2(float(min_a)), 3),
                relaxed_draws=relaxed_draws, relaxed_k5=relaxed_k5,
                fastpath_frac_k5=1.0 - per_k[5]["relaxed"] / float(per_k[5]["total"]),
                fastpath_frac_all=1.0 - sum(v["relaxed"] for v in per_k.values()) /
                float(sum(v["total"] for v in per_k.values())), records=recs)


def _bform_check(d, g):
    B4 = d["B4"]; G_sym = g["G_sym"]; GB = d["GB"]; bad = 0
    for i in (0, 1, 5, 6, 11):
        for j in (0, 1, 5, 6, 11):
            q = sum(B4[i][6 + u] * G_sym[u][v] * B4[j][6 + v] for u in range(6) for v in range(6))
            if GB[i][j] != WB * q:
                bad += 1
    wrong = V20._tb_gram(B4, G_sym, 0, WB)
    ndiff = sum(1 for i in range(12) for j in range(12) if GB[i][j] != wrong[i][j])
    return dict(checked=25, mismatches=bad, entries_differing_from_a_form=ndiff,
                bform_ok=(bad == 0 and ndiff > 0))


def _wrong_bform_control(orbits, g):
    """Negative control: rebuild the b-form from the a-coordinates (the 43.03 bug) on a
    32-draw control subset and show the level-0 node density changes.  The check fails
    if the correction is reverted."""
    rng = random.Random(SEED_BOX)
    cn = cl = wn = wl = 0
    sample = [d for o in orbits for d in o["draws"]][:32]
    for d in sample:
        smp = _sample_l1(d, g, rng)
        if smp is None:
            continue
        r1 = _exact_subtree(d, smp[0], g)
        saved = (d["ab"], d["mb"])
        d["ab"], d["mb"] = V20._tb_ldl6(V20._tb_gram(d["B4"], g["G_sym"], 0, WB))
        r2 = _exact_subtree(d, smp[0], g)
        d["ab"], d["mb"] = saved
        if r1 and r2 and r1[1] and r2[1]:
            cn += r1[0]; cl += r1[1]; wn += r2[0]; wl += r2[1]
    ratio = (wn / float(wl)) / (cn / float(cl)) if cl and wl else None
    return dict(draws=len(sample), n0_ratio_mean=(round(ratio, 4) if ratio else None),
                ok=(ratio is not None and ratio > 1.05 and cl > 0))


def measure(K_DRAW=4, seed_orb=SEED_ORB, m_ch=M_CH):
    t0 = time.time(); g = _geom(); orbits = _load_orbits(g, K_DRAW, seed_orb)
    V20.TB_MAX.update(dict(inv=0, b2=0, r2=0, d19=0, d2=0, ck=0, rem=0))
    for o in orbits:
        for d in o["draws"]:
            _prep(d, g); d["_T"] = g["T"]; d["_g"] = g
    rng_e = random.Random(SEED_BOX)
    for o in orbits:
        for d in o["draws"]:
            d["exact"] = _measure(d, g, rng_e, "exact", m_ch=m_ch)
    rng_f = random.Random(SEED_BOX)
    overflow = None
    for o in orbits:
        for d in o["draws"]:
            try:
                d["fp"] = _measure(d, g, rng_f, "fp", m_ch=max(1, m_ch // 4))
            except OverflowError as ex:
                overflow = str(ex); d["fp"] = None
    ex_ir = _ht(orbits, lambda d: d["exact"]["ir"] if d.get("exact") else None, K_DRAW)
    ex_se = _htse(orbits, lambda d: d["exact"]["ir"] if d.get("exact") else None, ex_ir, K_DRAW)
    fp_ir = _ht(orbits, lambda d: d["fp"]["ir"] if d.get("fp") else None, K_DRAW)
    fp_se = _htse(orbits, lambda d: d["fp"]["ir"] if d.get("fp") else None, fp_ir, K_DRAW)
    n0 = _ht(orbits, lambda d: d["exact"]["n0"] if d.get("exact") else None, K_DRAW)
    n1 = _ht(orbits, lambda d: d["exact"]["n1"] if d.get("exact") else None, K_DRAW)
    chl = _ht(orbits, lambda d: d["exact"]["chl"] if d.get("exact") else None, K_DRAW)
    spc = _ht(orbits, lambda d: d["exact"]["spc"] if d.get("exact") else None, K_DRAW)
    leaf = _ht(orbits, lambda d: d["exact"]["leaf"] if d.get("exact") else None, K_DRAW)
    tot = dict(draws=sum(len(o["draws"]) for o in orbits),
               completed=sum(1 for o in orbits for d in o["draws"] if d.get("exact")),
               fp_completed=sum(1 for o in orbits for d in o["draws"] if d.get("fp")),
               nodes0=sum(d["exact"]["nodes0"] for o in orbits for d in o["draws"] if d.get("exact")),
               leaves=sum(d["exact"]["leaves"] for o in orbits for d in o["draws"] if d.get("exact")),
               chords=sum(d["exact"]["chords"] for o in orbits for d in o["draws"] if d.get("exact")),
               ell_nodes=sum(d["exact"]["enodes"] for o in orbits for d in o["draws"] if d.get("exact")),
               ell_leaves=sum(d["exact"]["eleaves"] for o in orbits for d in o["draws"] if d.get("exact")))
    bf = _bform_check(orbits[0]["draws"][0], g)
    wrong = _wrong_bform_control(orbits, g)
    cal = alpha_scan_calibration(orbits, g)
    inv_max = 0; ec_max = 0
    for o in orbits:
        for d in o["draws"]:
            ctx = d["ctx"]
            for k in range(6):
                if ctx["oka"][k]:
                    inv_max = max(inv_max, ctx["inva"][k]); ec_max = max(ec_max, ctx["ecab"][k])
                if ctx["okb"][k]:
                    inv_max = max(inv_max, ctx["invb"][k]); ec_max = max(ec_max, ctx["ecab"][k])
    r2_cond = ((Q38 * inv_max) >> 64) + 1
    rad_cond = V20._isqrt_u40_lut9_q32(r2_cond)
    d19_cond = (rad_cond + ec_max) >> 13
    cond = dict(inv_bits=inv_max.bit_length(), r2_bits=r2_cond.bit_length(),
                d19_bits=d19_cond.bit_length(), d2_bits=(d19_cond * d19_cond).bit_length(),
                r2_lt_u40=(r2_cond < (1 << 40)))
    # fixed-path operand maxima (bit lengths)
    max_inv = max(max(d["ctx"]["inva"][k].bit_length() for k in range(6) if d["ctx"]["oka"][k]) for o in orbits for d in o["draws"])
    max_b2 = max(max(d["ctx"]["b2a"][k].bit_length() for k in range(6) if d["ctx"]["oka"][k]) for o in orbits for d in o["draws"])
    upper = ex_ir + 3.0 * ex_se
    return dict(k_draw=K_DRAW, draws=tot["draws"], completed=tot["completed"],
                fp_completed=tot["fp_completed"], nodes0=tot["nodes0"], leaves=tot["leaves"],
                chords=tot["chords"], ell_nodes=tot["ell_nodes"], ell_leaves=tot["ell_leaves"],
                wmean_ir=ex_ir, se_ir=ex_se, wmean_plus_3se_ir=upper, upper_103=1.03 * upper,
                n0=n0, n1=n1, chords_per_leaf=chl, sp_per_chord=spc, leaf_ir=leaf,
                fp_wmean_ir=fp_ir, fp_se_ir=fp_se, fp_wmean_plus_3se_ir=fp_ir + 3.0 * fp_se,
                alpha=cal, alpha_scan_estimator=cal["scan"],
                bform=bf, wrong_bform=wrong, cond_bounds=cond,
                correction_ok=(bf["bform_ok"] and wrong["ok"]),
                overflow=overflow, overflow_free=(overflow is None),
                max_inv_bits=V20.TB_MAX["inv"].bit_length(), max_b2_bits=V20.TB_MAX["b2"].bit_length(),
                max_r2_bits=V20.TB_MAX["r2"].bit_length(), max_d19_bits=V20.TB_MAX["d19"].bit_length(),
                max_d2_bits=V20.TB_MAX["d2"].bit_length(), max_ck_bits=V20.TB_MAX["ck"].bit_length(),
                pins=dict(mean=PIN_MEAN, se=PIN_SE, p3se=PIN_P3SE, upper_103=PIN_103),
                pin_ok=(abs(ex_ir - PIN_MEAN) <= 2e-6 and abs(ex_se - PIN_SE) <= 2e-6),
                elapsed_s=round(time.time() - t0, 3))


if __name__ == "__main__":
    r = measure()
    print({k: r[k] for k in ("k_draw", "draws", "completed", "fp_completed", "nodes0",
                             "leaves", "wmean_ir", "se_ir", "wmean_plus_3se_ir", "upper_103",
                             "fp_wmean_ir", "fp_se_ir", "overflow", "max_inv_bits", "max_b2_bits",
                             "pin_ok", "elapsed_s")})
    print("alpha", r["alpha"])
    print("bform", r["bform"])
