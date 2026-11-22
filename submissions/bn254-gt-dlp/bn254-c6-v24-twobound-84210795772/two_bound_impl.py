"""two_bound_impl.py — production-scale two-bound product-region Fincke-Pohst
measurement for the BN254 product region.

The tree keeps two rank-6 weighted-form budgets (F_a = W_a Q(a), F_b = W_b Q(b)) in
the package's 12D LLL basis and prunes on either budget; the level-0 chord is the
intersection of the two exact level-0 intervals.  This module measures the tree at
production parameters by sampling level-1 nodes (x_2..x_11) from the ellipsoid
volume measure, checking the level-5..2 a/b ancestors, and enumerating x_1 and the
x_0 chord exactly.  It uses exact Fractions for the two form LDLs and the package's
integer quadratic for the level-0 product chord; the chord split uses the package's
fast chord statistics.

Pure stdlib.  Measurement only; the shipped kernel would use Q38/Q32 fixed point.
"""
import math, random
from fractions import Fraction as Fr
from math import sqrt, log2, ceil, floor
import model_bn254_v6 as M
import descent_bn254_v6 as DESC
import v20_lat_b2 as V20

WA, WB = 1764, 1000
Q32 = 1 << 32
N_STRATA = 32
DEF_K_DRAW = 32
DEF_M_CH = 8


def _mat_gram(B, G_sym, part):
    n = 12
    G = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            q = sum(B[i][u] * G_sym[u][v] * B[j][v] for u in range(6) for v in range(6))
            G[i][j] = (WA if part == 'a' else WB) * q
    return G


def _ldl(G, npos):
    n = len(G)
    Dv = [1] + [0] * n
    lam = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(min(i, npos - 1) + 1):
            u = G[i][j]
            for k in range(j):
                u = (Dv[k + 1] * u - lam[i][k] * lam[j][k]) // Dv[k]
            if j < i:
                lam[i][j] = u
            else:
                Dv[i + 1] = u
    alpha = [Fr(Dv[k + 1], Dv[k]) for k in range(npos)]
    mu = [[Fr(0) for _ in range(n)] for _ in range(n)]
    for k in range(npos):
        for j in range(k + 1, n):
            mu[j][k] = Fr(lam[j][k], Dv[k + 1])
    return alpha, mu


def _root_interval(rem, alpha, c):
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


def _in_swar_box(v, S, S2):
    a0, b0 = v[:6], v[6:]
    A3M = M.A_TAB_SLOTS_A3; B0M = M.B0_ARENA_MAX
    for ta, tb in ((a0, b0), (V20._mv6(S, a0), V20._mv6(S, b0)), (V20._mv6(S2, a0), V20._mv6(S2, b0))):
        if not (-510 <= ta[0] <= 510 and -510 <= ta[1] <= 510 and -510 <= ta[2] <= 510 and -A3M < ta[3] < A3M and -510 <= ta[4] <= 510 and -510 <= ta[5] <= 510):
            return False
        if not (-B0M < tb[0] < B0M and -1022 <= tb[1] <= 1022 and -510 <= tb[2] <= 510 and -510 <= tb[3] <= 510 and -510 <= tb[4] <= 510 and -510 <= tb[5] <= 510):
            return False
    return True


def _fstats(base, B0, l0, h0, S, S2):
    """Package fast chord stats with the SWAR box-convexity endpoint clip."""
    clips = 0
    while l0 <= h0 and not _in_swar_box([base[c] + l0 * B0[c] for c in range(12)], S, S2):
        clips += 1; l0 += 1
    while l0 <= h0 and not _in_swar_box([base[c] + h0 * B0[c] for c in range(12)], S, S2):
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
    for c, m in ((a3 - a4, d3 - d4), (a3 + a4, d3 + d4), (a3 - a5, d3 - d5),
                 (a3 + a5, d3 + d5), (a4 - a5, d4 - d5), (a4 + a5, d4 + d5)):
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
            c, m = vlb0_k[ki], db0_k[ki]; v0, v1 = c, c + umax * m
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


def _fstats_product(base, B0, l0, h0, S, S2):
    return _fstats(base, B0, l0, h0, S, S2)


def _splus(st):
    return st["segs"] * 12.0 + st["bspl"] * 12.0 + st["fspl"] * 22.0 + \
        st["cross"] * 10.0 + st["clips"] * 60.0 + 16.0


def measure_two_bound_tree(seed=20261122, K_DRAW=DEF_K_DRAW, M_CH=DEF_M_CH):
    """Return production two-bound tree statistics and the HT per-leaf cost."""
    dat = M.load_data(); T4 = dat["t4_geom"]; PQ = dat["pool_quad"]
    G_sym, D_mats, UMATS, S, S2, G_SC = M.build_exact_int_geom(T4)
    RA = PQ["ra_op"]; RB = PQ["rb_ga_op"]
    RA2W = int(round(WA * (RA ** 2) * G_SC)); RB2W = int(round(WB * (RB ** 2) * G_SC))
    QA_B = RA2W // WA; QB_B = RB2W // WB; T_BOX = int(PQ["t12_box"])
    p_lat_f = M.OP["p_lat"]; lx2_f = M.OP["lx2"]

    def _load_orbits():
        rng_orb = random.Random(seed); orbits = []
        for s_idx in range(N_STRATA):
            lo_b = p_lat_f + (lx2_f - p_lat_f) * s_idx / N_STRATA
            hi_b = p_lat_f + (lx2_f - p_lat_f) * (s_idx + 1) / N_STRATA
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
                draws.append(dict(pp4=pp4, m_p=m_p4, B4=B4, D4=D4, lm4=lm4))
            orbits.append(dict(M_s=M_s, S_T=S_T, draws=draws))
        return orbits

    def _prep(d):
        B4 = d["B4"]; D4 = d["D4"]; lm4 = d["lm4"]
        GA = _mat_gram(B4, G_sym, 'a'); GB = _mat_gram(B4, G_sym, 'b')
        aa, ma = _ldl(GA, 6); ab, mb = _ldl(GB, 6)
        Xm = [0] * 12
        for k in range(11, -1, -1):
            rk = math.isqrt(((D4[k] * T_BOX) + D4[k + 1] - 1) // D4[k + 1]) + 1
            Xm[k] = rk + sum((Xm[j] + 1) // 2 for j in range(k + 1, 12))
        d["aa"] = aa; d["ma"] = ma; d["ab"] = ab; d["mb"] = mb
        d["ec"] = [sum(Xm[j] for j in range(k + 1, 12)) + 1 for k in range(12)]
        d["b2"] = [min((1 << 64) - 1, (D4[k + 1] << 64) // (D4[k] * T_BOX)) for k in range(12)]
        d["inv"] = [((D4[k] * T_BOX << 38) + D4[k + 1] - 1) // D4[k + 1] for k in range(12)]
        d["mu32"] = [[((lm4[j][k] << 32) // D4[k + 1]) if j > k else 0 for k in range(12)] for j in range(12)]
        d["Bn"] = [D4[k + 1] / D4[k] for k in range(12)]

    rng_box = random.Random(0xB0B0B0)

    def _sample_l1(d):
        D4 = d["D4"]; mu32 = d["mu32"]; ec = d["ec"]; b2 = d["b2"]; inv = d["inv"]
        Bn = d["Bn"]; R_f = sqrt(float(T_BOX)); inv_sB = [1.0 / sqrt(Bn[k]) for k in range(12)]
        mu_c = [[mu32[j][k] for j in range(k + 1, 12)] for k in range(12)]
        gv = [rng_box.gauss(0.0, 1.0) for _ in range(11)]
        rg = R_f * (rng_box.random() ** (1.0 / 10.0)) / sqrt(math.sumprod(gv, gv))
        xp = [0] * 12; rm = 1 << 38
        for k in range(11, 1, -1):
            ck = -math.sumprod(mu_c[k], xp[k + 1:]) if k < 11 else 0
            if abs(ck) > 10 ** 12:
                return None
            xk = int(round(ck * (1.0 / (1 << 32)) + gv[k - 1] * rg * inv_sB[k])); xp[k] = xk
            if abs(xk) > 10 ** 9:
                return None
            ad = abs((xk << 32) - ck)
            if ad > ec[k]:
                d19 = (ad - ec[k]) >> 13; sub = ((d19 * d19) * b2[k]) >> 64
                if sub > rm:
                    return None
                rm -= sub
        return xp, rm

    def _prod_interval(base, B0, l0e, h0e):
        a0 = base[:6]; b0 = base[6:]; u = B0[:6]; w = B0[6:]
        Aa = V20._qf6(G_sym, u); Ab = V20._qf6(G_sym, w)
        Ba = 2 * sum(a0[i] * G_sym[i][j] * u[j] for i in range(6) for j in range(6))
        Bb = 2 * sum(b0[i] * G_sym[i][j] * w[j] for i in range(6) for j in range(6))
        Ca = V20._qf6(G_sym, a0); Cb = V20._qf6(G_sym, b0)
        def iv(A, Bv, Cv, T):
            if A <= 0:
                return None
            disc = Bv * Bv - 4 * A * (Cv - T)
            if disc < 0:
                return None
            s = math.isqrt(disc)
            lo = math.ceil(Fr(-Bv - s, 2 * A)); hi = math.floor(Fr(-Bv + s, 2 * A))
            while lo <= hi and A * lo * lo + Bv * lo + Cv > T:
                lo += 1
            while A * (lo - 1) * (lo - 1) + Bv * (lo - 1) + Cv <= T:
                lo -= 1
            while lo <= hi and A * hi * hi + Bv * hi + Cv > T:
                hi -= 1
            while A * (hi + 1) * (hi + 1) + Bv * (hi + 1) + Cv <= T:
                hi += 1
            return lo, hi
        ia = iv(Aa, Ba, Ca, QA_B); ib = iv(Ab, Bb, Cb, QB_B)
        if ia is None or ib is None:
            return None
        lo = max(ia[0], ib[0], l0e); hi = min(ia[1], ib[1], h0e)
        return (lo, hi) if lo <= hi else None

    def _l1_subtree(d, xp, rm):
        aa = d["aa"]; ma = d["ma"]; ab = d["ab"]; mb = d["mb"]
        mu32 = d["mu32"]; ec = d["ec"]; inv = d["inv"]; b2 = d["b2"]
        x = [0] * 12
        for k in range(2, 12):
            x[k] = xp[k]
        Sa = Fr(0); Sb = Fr(0)
        for k in range(2, 6):
            ca = sum((ma[j][k] * x[j] for j in range(k + 1, 12)), Fr(0))
            Sa += aa[k] * (x[k] + ca) ** 2
            cb = sum((mb[j][k] * x[j] for j in range(k + 1, 12)), Fr(0))
            Sb += ab[k] * (x[k] + cb) ** 2
        rema = RA2W - Sa; remb = RB2W - Sb
        if rema < 0 or remb < 0:
            return None
        # full form level 1: Q38 fixed point, c1 in Q32
        c1 = -sum((mu32[j][1] * x[j] for j in range(2, 12)), 0)
        r2 = ((rm * inv[1]) >> 64) + 1
        if r2 < 1 or r2 >= (1 << 40):
            return None
        rad = V20._isqrt_u40_lut9_q32(r2)
        lo1 = (c1 - ec[1] - rad + Q32 - 1) >> 32; hi1 = (c1 + ec[1] + rad) >> 32
        ca1 = sum((ma[j][1] * x[j] for j in range(2, 12)), Fr(0))
        la1, ha1 = _root_interval(rema, aa[1], ca1)
        cb1 = sum((mb[j][1] * x[j] for j in range(2, 12)), Fr(0))
        lb1, hb1 = _root_interval(remb, ab[1], cb1)
        lo = max(lo1, la1, lb1); hi = min(hi1, ha1, hb1)
        if lo > hi:
            return None
        nodes0 = leaves = chords = 0; fast = bnd = 0; sp = 0.0
        for x1 in range(lo, hi + 1):
            x[1] = x1
            t1 = (x1 << 32) - c1
            ad = -t1 if t1 < 0 else t1
            rem0 = rm
            if ad > ec[1]:
                d19 = (ad - ec[1]) >> 13; sub = ((d19 * d19) * b2[1]) >> 64
                if sub > rm:
                    continue
                rem0 = rm - sub
            ra0 = rema - aa[1] * (x1 + ca1) ** 2
            rb0 = remb - ab[1] * (x1 + cb1) ** 2
            if ra0 < 0 or rb0 < 0:
                continue
            nodes0 += 1
            base = [sum(x[j] * d["B4"][j][c] for j in range(1, 12)) for c in range(12)]
            pv = _prod_interval(base, d["B4"][0], -10 ** 12, 10 ** 12)
            if pv is None:
                continue
            c0 = -sum((mu32[j][0] * x[j] for j in range(1, 12)), 0)
            r20 = ((rem0 * inv[0]) >> 64) + 1
            if r20 < 1 or r20 >= (1 << 40):
                continue
            rad0 = V20._isqrt_u40_lut9_q32(r20)
            l0f = (c0 - ec[0] - rad0 + Q32 - 1) >> 32; h0f = (c0 + ec[0] + rad0) >> 32
            l = max(pv[0], l0f); h = min(pv[1], h0f)
            if l > h:
                continue
            chords += 1; leaves += (h - l + 1)
            st = _fstats_product(base, d["B4"][0], l, h, S, S2)
            fast += st["fast"]; bnd += st["bnd"]; sp += _splus(st)
        x[1] = 0
        return nodes0, leaves, chords, fast, bnd, sp

    def _measure(d):
        reach = 0; nodes0 = leaves = chords = fast = bnd = 0; sp = 0.0
        for _ in range(M_CH):
            sres = _sample_l1(d)
            if sres is None:
                continue
            xp, rm = sres
            r = _l1_subtree(d, xp, rm)
            if r is None:
                continue
            reach += 1; nodes0 += r[0]; leaves += r[1]; chords += r[2]
            fast += r[3]; bnd += r[4]; sp += r[5]
        if leaves == 0:
            return None
        leaf = (fast * 13.35 + bnd * 53.05) / leaves
        spc = sp / chords if chords else 0.0
        chl = chords / float(leaves)
        n0 = nodes0 / float(leaves); n1 = reach / float(leaves)
        C0 = 74.0; C1 = 115.0; APPLY = 2.85
        return dict(reach=reach, nodes0=nodes0, leaves=leaves, chords=chords,
                    n0=n0, n1=n1, chl=chl, spc=spc, leaf=leaf,
                    ir=n0 * C0 + n1 * C1 + chl * spc + leaf + APPLY)

    t0 = __import__("time").time()
    orbits = _load_orbits()
    for o in orbits:
        for d in o["draws"]:
            _prep(d)
    for o in orbits:
        for d in o["draws"]:
            d["l1"] = _measure(d)

    def _ht(key):
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

    def _htse(key, mean):
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

    ir = _ht(lambda d: d["l1"]["ir"] if d.get("l1") else None)
    se = _htse(lambda d: d["l1"]["ir"] if d.get("l1") else None, ir)
    n0 = _ht(lambda d: d["l1"]["n0"] if d.get("l1") else None)
    n1 = _ht(lambda d: d["l1"]["n1"] if d.get("l1") else None)
    chl = _ht(lambda d: d["l1"]["chl"] if d.get("l1") else None)
    spc = _ht(lambda d: d["l1"]["spc"] if d.get("l1") else None)
    leaf = _ht(lambda d: d["l1"]["leaf"] if d.get("l1") else None)
    tot_draws = sum(1 for o in orbits for d in o["draws"] if d.get("l1"))
    tot_nodes = sum(d["l1"]["nodes0"] for o in orbits for d in o["draws"] if d.get("l1"))
    tot_leaves = sum(d["l1"]["leaves"] for o in orbits for d in o["draws"] if d.get("l1"))
    return dict(k_draw=K_DRAW, m_ch=M_CH, strata=N_STRATA,
                draws=sum(len(o["draws"]) for o in orbits), completed=tot_draws,
                nodes0=tot_nodes, leaves=tot_leaves,
                wmean_ir=ir, se_ir=se, wmean_plus_3se_ir=ir + 3.0 * se,
                n0=n0, n1=n1, chords_per_leaf=chl, sp_per_chord=spc, leaf_ir=leaf,
                C0=74.0, C1=115.0, apply=2.85,
                elapsed_s=round(__import__("time").time() - t0, 3))
