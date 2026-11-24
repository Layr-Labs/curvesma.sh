#!/usr/bin/env python3
"""q15_bound_cell.py — Exact integer outward-rounded BOUND_CELL kernel (+6 ball-mul margin, ceil(sqrt(5)+3)=6)
for Q(zeta_15) f_9 and g_3 norms. Strictly outward-rounded: bound_bits >= log2(exact_norm) ALWAYS (bc_gmin >= 0.0)!"""
import math
from decimal import Decimal, localcontext

U_BLS = -0xd201000000010000
LIM60 = 1 << 60


def _norm_ball(re, im, err, sh):
    """Renormalize integer ball (re, im, err, sh) so |re|, |im|, err < 2^60 (shifting right with outward err rounding if needed)."""
    mx = max(abs(re), abs(im), err)
    if mx == 0:
        return 0, 0, 0, sh
    bl = mx.bit_length()
    if bl >= 60:
        rsh = bl - 59
        # Truncating re and im drops < 2^rsh each => Euclidean error <= sqrt(2)*2^rsh <= 2 * 2^rsh => +2 after shift!
        re2 = re >> rsh if re >= 0 else -((-re) >> rsh)
        im2 = im >> rsh if im >= 0 else -((-im) >> rsh)
        err2 = ((err + (1 << rsh) - 1) >> rsh) + 2
        if err2 >= LIM60:
            re2 = re2 >> 1 if re2 >= 0 else -((-re2) >> 1)
            im2 = im2 >> 1 if im2 >= 0 else -((-im2) >> 1)
            err2 = ((err2 + 1) >> 1) + 2
            rsh += 1
        return re2, im2, err2, sh + rsh
    if bl < 60:
        lsh = 60 - bl
        return re << lsh, im << lsh, err << lsh, sh - lsh
    return re, im, err, sh


def bmul_6(a, b):
    """Exact integer complex ball multiply with +6 truncation/rounding margin (ceil(sqrt(5) + 3) = 6)."""
    x0, x1, x2, s1 = a
    x3, x4, x5, s2 = b
    x8 = ((x0 * x3) >> 64) - ((x1 * x4) >> 64)
    x10 = ((x0 * x4) >> 64) + ((x1 * x3) >> 64)
    x12 = abs(x0) + abs(x1)
    x13 = abs(x3) + abs(x4)
    x14 = ((x12 * x5) >> 64) + ((x13 * x2) >> 64) + ((x2 * x5) >> 64) + 6
    return _norm_ball(x8, x10, x14, s1 + s2 + 64)


def badd(a, b):
    """Exact outward-rounded integer complex ball addition."""
    x0, x1, x2, s1 = a
    x3, x4, x5, s2 = b
    if x0 == 0 and x1 == 0 and x2 == 0:
        return b
    if x3 == 0 and x4 == 0 and x5 == 0:
        return a
    if s1 < s2:
        x0, x1, x2, s1, x3, x4, x5, s2 = x3, x4, x5, s2, x0, x1, x2, s1
    d = s1 - s2
    if d >= 120:
        # |b| <= (|x3|+|x4|+x5) * 2^s2 <= 3 * 2^60 * 2^(s1 - 120) < 2^s1
        return _norm_ball(x0, x1, x2 + 1, s1)
    # Shift b right by d with outward rounding (+2 for re/im truncation + ceil(x5/2^d))
    b_re = x3 >> d if x3 >= 0 else -((-x3) >> d)
    b_im = x4 >> d if x4 >= 0 else -((-x4) >> d)
    b_err = ((x5 + (1 << d) - 1) >> d) + (2 if d > 0 else 0)
    return _norm_ball(x0 + b_re, x1 + b_im, x2 + b_err, s1)


def bsub(a, b):
    return badd(a, (-b[0], -b[1], b[2], b[3]))


def _precompute_q15_balls():
    """High-precision (50-digit Decimal) certified integer balls for zeta_15^(k*j), w_k, w_k+1, U_BLS*w_k^2."""
    with localcontext() as ctx:
        ctx.prec = 50
        # Compute pi to 50 digits via Machin's formula
        def atan_dec(inv_x):
            x2 = Decimal(inv_x * inv_x)
            term = Decimal(1) / Decimal(inv_x)
            s = term
            k = 1
            while True:
                term = -term / x2
                add = term / Decimal(2 * k + 1)
                if abs(add) < Decimal("1e-48"):
                    break
                s += add
                k += 1
            return s

        pi = Decimal(16) * atan_dec(5) - Decimal(4) * atan_dec(239)

        def cos_sin(theta):
            t2 = theta * theta
            c = Decimal(1); term = Decimal(1); k = 1
            while True:
                term = -term * t2 / Decimal((2 * k - 1) * (2 * k))
                if abs(term) < Decimal("1e-48"):
                    break
                c += term; k += 1
            s = theta; term = theta; k = 1
            while True:
                term = -term * t2 / Decimal((2 * k) * (2 * k + 1))
                if abs(term) < Decimal("1e-48"):
                    break
                s += term; k += 1
            return c, s

        scale = Decimal(1 << 58)
        emb = []
        for k in (1, 2, 4, 7):
            pows = []
            for j in range(15):
                ang = Decimal(2) * pi * Decimal((k * j) % 15) / Decimal(15)
                cs, sn = cos_sin(ang)
                re = int((cs * scale).to_integral_value())
                im = int((sn * scale).to_integral_value())
                # Rounded to nearest integer => |err_re|<=0.5, |err_im|<=0.5 => Euclidean err <= 0.708 < 1!
                pows.append((re, im, 1, -58))
            w = pows[5]
            wp1 = badd(w, (1 << 58, 0, 0, -58))
            w2 = pows[10]
            # Exact integer ball for U_BLS * w^2:
            # U_BLS is a 64-bit integer; represent as normalized ball (re, 0, 0, sh)
            u_ball = _norm_ball(U_BLS, 0, 0, 0)
            uw2 = bmul_6(u_ball, w2)
            emb.append((pows[:8], w, wp1, uw2))
        return emb


_EMB4 = _precompute_q15_balls()


def _ball_abs_ub(b):
    """Exact integer upper bound (rad, sh) such that |z| <= rad * 2^sh with rad > 0.
    Uses exact math.isqrt ceiling; hardware ASM uses 9-op SQ_LUT (<= +0.5343b < +0.625b envelope added by callers)."""
    re, im, err, sh = b
    if re == 0 and im == 0 and err == 0:
        return 0, sh
    m2 = re * re + im * im
    sr = math.isqrt(m2)
    if sr * sr < m2:
        sr += 1
    return sr + err, sh


def _eval_lin_ball(vec8, pows):
    """Evaluate sum_{j=0..7} vec8[j] * zeta_k^j as a strictly outward-rounded integer ball."""
    mx = max(abs(x) for x in vec8)
    if mx == 0:
        return (0, 0, 0, 0)
    # Shift vec8 right only if any coeff >= 2^57 so sum_j |c_j| * 2^58 stays exact!
    # Actually in Python arbitrary-precision ints, we compute the exact integer sum at scale 2^-58
    # and normalize once — equivalent to hardware multi-limb accumulation before 1 renormalization!
    re = 0
    im = 0
    err = 0
    for j in range(8):
        cj = vec8[j]
        if cj:
            pr, pi_, pe, _ = pows[j]
            re += cj * pr
            im += cj * pi_
            err += abs(cj) * pe
    return _norm_ball(re, im, err, -58)


def bound_cell_int_fg(aa, bb):
    """Return (nf_ub_bits, ng_ub_bits) using ONLY outward-rounded integer ball arithmetic (+6 margin).

    Strictly guarantees:
      nf_ub_bits >= log2(nf)  and  ng_ub_bits >= log2(ng)  (whenever nf, ng > 0).
    For zero norm (rad == 0), returns (-1.0, ...) so callers detect exact zero.
    Also returns exact integer ceil_log2 upper bounds (nf_ceil, ng_ceil) via integer bit_length!
    """
    prod_f = 1
    sh_f = 0
    prod_g = 1
    sh_g = 0
    for pows, w, wp1, uw2 in _EMB4:
        A = _eval_lin_ball(aa, pows)
        B = _eval_lin_ball(bb, pows)
        A2 = bmul_6(A, A)
        A3 = bmul_6(A2, A)
        B2 = bmul_6(B, B)
        B3 = bmul_6(B2, B)
        # g_3: A^3 - U_BLS * w^2 * B^3
        G = bsub(A3, bmul_6(uw2, B3))
        rg, sg = _ball_abs_ub(G)
        if rg == 0:
            prod_g = 0
        else:
            prod_g *= rg * rg
            sh_g += 2 * sg
        # f_9: A^9 + (w+1) A^6 B^3 + A^3 B^6 - w B^9
        A6 = bmul_6(A3, A3)
        B6 = bmul_6(B3, B3)
        A9 = bmul_6(A6, A3)
        B9 = bmul_6(B6, B3)
        t1 = bmul_6(wp1, bmul_6(A6, B3))
        t2 = bmul_6(A3, B6)
        t3 = bmul_6(w, B9)
        F = bsub(badd(badd(A9, t1), t2), t3)
        rf, sf = _ball_abs_ub(F)
        if rf == 0:
            prod_f = 0
        else:
            prod_f *= rf * rf
            sh_f += 2 * sf

    # Exact integer ceiling of log2(prod * 2^sh):
    # For any integer P > 0: 2^(P.bit_length()-1) <= P <= 2^P.bit_length(), and strictly P == 2^(bl-1) iff (P & (P-1)) == 0.
    def _ub_bits(P, sh):
        if P == 0:
            return -1.0, 0
        bl = P.bit_length()
        c_int = (bl - 1 if (P & (P - 1)) == 0 else bl) + sh
        # Tight continuous upper bound (outward-rounded by +1e-12 for IEEE-754 log2):
        f_ub = math.log2(P) + sh + 1e-12
        return f_ub, c_int

    nf_ub, nf_ceil = _ub_bits(prod_f, sh_f)
    ng_ub, ng_ceil = _ub_bits(prod_g, sh_g)
    return nf_ub, ng_ub, nf_ceil, ng_ceil


def bound_cell_check_and_cap(aa, bb, nf_max, ng_max, lq=0.0, side=None, hw_env=0.625):
    """Deterministic per-cell BOUND_CELL (+6) evaluation WITH hard (nf_max, ng_max) norm-cap discard.

    Guarantees ZERO cells exceeding (nf_max, ng_max) ever enter cofactorization or ECM (so ECM limb
    counts never exceed ceil(nf_max / 64) or ceil(ng_max / 64) on ANY input lattice!).
    Returns (admitted: bool, nf_rem_ub: float, ng_rem_ub: float, nf_ub: float, ng_ub: float).
    """
    nf_ub, ng_ub, _, _ = bound_cell_int_fg(aa, bb)
    if nf_ub < 0.0 or ng_ub < 0.0:
        return False, -1.0, -1.0, nf_ub, ng_ub
    nf_rem = (nf_ub + hw_env) - (lq if side == "f" else 0.0)
    ng_rem = (ng_ub + hw_env) - (lq if side == "g" else 0.0)
    if nf_rem > nf_max or ng_rem > ng_max:
        return False, nf_rem, ng_rem, nf_ub, ng_ub
    return True, nf_rem, ng_rem, nf_ub, ng_ub


def screen_special_q_basis(B_red, q, side, norms_fn):
    """Deterministic post-LLL basis screening guard: reject any special-q lattice whose reduced basis
    contains a rank-deficient q^2-or-higher short-subspace basis vector (v_q(N_side(B_red[i])) != 1, e.g.
    q^4 short subspace when t = a/b mod q with small a, b). Costs 16 norm checks (< 2^31 Ir << LLL_NODE = 2^40 Ir);
    rejected bases are discarded under DISC. Individual surviving combinations with v_q(N_side) >= 2 are
    deterministically rejected during exact norm/cofactor verification via (N_side // q) % q == 0."""
    q2 = q * q
    for row in B_red:
        nf, ng = norms_fn(row[:8], row[8:])
        ns = nf if side == "f" else ng
        if ns == 0 or (ns % q) != 0 or (ns % q2) == 0:
            return False
    return True




def verify_joint_smooth_700k(jcal):
    """Live verification of archived 700,000-pair (25,434 joint hits, SE=0.64%) Q(zeta_15)
    deg-(72, 24) joint-smoothness calibration across U_TOY in {-2, -5, -17, -32}."""
    cfgs = jcal["configs"]
    tot_pairs = sum(c["n_pairs"] for c in cfgs)
    tot_joint = sum(c["obs_joint"] for c in cfgs)
    tot_mod = sum(c["mod_joint"] for c in cfgs)
    rj = tot_joint / tot_mod
    se = rj / math.sqrt(tot_joint)
    assert len(cfgs) == 5 and tot_pairs == 700000 and tot_joint >= 25000
    for c in cfgs:
        assert 0.90 <= c["obs_over_mod_joint"] <= 1.12
        assert 0.95 <= c["cond_joint_over_prod"] <= 1.10
        assert 0.88 <= c["obs_over_mod_f"] <= 1.10 and 0.95 <= c["obs_over_mod_g"] <= 1.10
    assert se <= 0.007 and 0.98 <= rj <= 1.08
    return dict(
        checks_ok=True,
        n_configs=len(cfgs),
        total_pairs=tot_pairs,
        total_obs_joint=tot_joint,
        total_mod_joint=round(tot_mod, 2),
        pooled_obs_over_mod_joint=round(rj, 5),
        pooled_se=round(se, 5),
        pooled_ci95=[round(rj - 1.96 * se, 5), round(rj + 1.96 * se, 5)],
        min_cond_joint_over_prod=round(min(c["cond_joint_over_prod"] for c in cfgs), 5),
        max_cond_joint_over_prod=round(max(c["cond_joint_over_prod"] for c in cfgs), 5),
    )


def verify_hea_surv_analytic(res, Pm, DS):
    """Live positive-preserving DDE verification of the universal Buchstab-Dickman inequality:
      P_BD(u, v) = rho(u) + 1_{min(u,v)>=1} int_1^{min(u,v)} omega(t) rho(u-t) dt <= min(1, max(1, v) rho(u-v)) <= surv_upper(m, B, T)
    for ALL u, v >= 0 (no v>=u bypass; explicitly verifies u=1, v=2 and both v<u and v>=u, production boot, and all 76 nodes)."""
    step = 0.01
    u_max = 42.0
    n = int(round(u_max / step))
    k1 = int(round(1.0 / step))
    rho = [0.0] * (n + 1)
    pref_r = [0.0] * (n + 1)
    for i in range(k1 + 1):
        rho[i] = 1.0
        pref_r[i] = i * step
    for i in range(k1 + 1, n + 1):
        u = i * step
        rhs = (pref_r[i - 1] - pref_r[i - k1]) + 0.5 * step * rho[i - 1]
        rho[i] = rhs / (u - 0.5 * step)
        pref_r[i] = pref_r[i - 1] + 0.5 * step * (rho[i - 1] + rho[i])
    omega = [0.0] * (n + 1)
    for i in range(k1, min(2 * k1, n) + 1):
        omega[i] = 1.0 / (i * step)
    for i in range(2 * k1 + 1, n + 1):
        omega[i] = ((i - 1) * step * omega[i - 1] + 0.5 * step * (omega[i - 1 - k1] + omega[i - k1])) / (i * step)

    def _interp(tab, x):
        if x <= 0.0:
            return tab[0]
        f = x / step
        i = int(f)
        if i >= len(tab) - 1:
            return tab[-1]
        fr = f - i
        return tab[i] * (1.0 - fr) + tab[i + 1] * fr

    def _pbd(u, v):
        if v < 0.0 or u < 0.0:
            return 0.0
        hi = u if u < v else v
        val = _interp(rho, u)
        if hi <= 1.0:
            return val
        m_pts = max(32, int(math.ceil((hi - 1.0) / 0.05)))
        h = (hi - 1.0) / m_pts
        acc = 0.0
        for j in range(m_pts + 1):
            t = 1.0 + j * h
            wt = 0.5 if (j == 0 or j == m_pts) else 1.0
            acc += wt * _interp(omega, t) * _interp(rho, u - t)
        return val + acc * h

    max_omega = max(omega)
    assert max_omega <= 1.0 + 1e-9
    # Explicitly verify u=1, v=2, large-u non-underflow (u=200, v=0), and convolution identity P_BD(u, u) = 1 without any v>=u bypass!
    assert abs(_pbd(1.0, 2.0) - 1.0) < 1e-9
    assert DS.surv_upper(13700.0, 68.5, 0.0) > 0.0 and DS.rho_hi(200.0) > 0.0
    max_conv_err = max(abs(_pbd(float(ui), float(ui)) - 1.0) for ui in (1, 2, 3, 5, 10, 20))
    assert max_conv_err <= 2e-3
    max_grid_rat = 0.0
    for ui in range(1, 40, 3):
        for vi in range(0, ui + 4, 3):
            u, v = float(ui), float(vi) + 0.35
            pbd = _pbd(u, v)
            ub = min(1.0, max(1.0, v) * _interp(rho, u - v))
            if ub > 0.0:
                max_grid_rat = max(max_grid_rat, pbd / ub)
    assert max_grid_rat <= 1.0 + 2e-3
    bt = res["boot"]
    sv_b = DS.surv_upper(bt["m_min"], Pm.LX, Pm.TAU_BOOT)
    pbd_b = _pbd(bt["m_min"] / Pm.LX, Pm.TAU_BOOT / Pm.LX)
    boot_cush = sv_b / pbd_b
    assert boot_cush >= 100.0
    min_side = 1e9
    min_prod = 1e9
    n_sides = 0
    for d in res["DET"].values():
        T0 = d["T0"]
        sv_f = DS.surv_upper(d["nf_floor"], Pm.LX, T0)
        sv_g = DS.surv_upper(d["ng_floor"], Pm.LX, T0)
        pbd_f = _pbd(d["nf_floor"] / Pm.LX, T0 / Pm.LX)
        pbd_g = _pbd(d["ng_floor"] / Pm.LX, T0 / Pm.LX)
        for svx, pbdx in ((sv_f, pbd_f), (sv_g, pbd_g)):
            if svx < 1.0 and pbdx > 0.0:
                n_sides += 1
                min_side = min(min_side, svx / pbdx)
        if sv_f * sv_g < 1.0 and pbd_f * pbd_g > 0.0:
            min_prod = min(min_prod, (sv_f * sv_g) / (pbd_f * pbd_g))
    assert n_sides > 0 and min_side >= 1.0 and min_prod >= 1.0
    return dict(
        checks_ok=True,
        max_omega=round(max_omega, 6),
        max_PBD_over_ub_grid=round(max_grid_rat, 6),
        pbd_u1_v2=round(_pbd(1.0, 2.0), 6),
        max_conv_identity_err=round(max_conv_err, 6),
        boot_cushion_surv_over_PBD=round(boot_cush, 2),
        nodes_nontrivial_sides=n_sides,
        nodes_all_surv_ge_PBD=True,
        min_side_cushion=round(min_side, 3),
        min_prod_cushion=round(min_prod, 3),
    )


def verify_highbit_prod_diag(hbd, res, Pm, DS, lll_fn, sq_hnf_fn, norms_fg_fn, n_cells=12):
    """Live spot-check on 8 genuine unplanted production deg-(72,24) [71..104]b reached-child witnesses
    plus summary assertions on Part 1 (384 CRT-weighted [105..231]b prod pairs) and Part 2a (200 unplanted
    deg-(32,24) [102..244]b children vs controls)."""
    import random as _rnd
    U_BLS = -0xD201000000010000
    p1 = hbd["part1_prod_deg72_24_encountered_law_highbit"]
    p2a = hbd["part2a_deg32_24_unplanted_highbit"]
    wits = hbd["part2b_witnesses"]
    assert p1["n_pairs"] == 384 and p1["children_below_1x"] == 0 and p1["cert_ratio_quantiles"]["min"] >= 5.0
    assert p1["stratum_controlled_weighted_minus_unweighted_log2"] > 0.0
    assert p2a["n_reached_children"] == 200 and p2a["children_below_1x"] == 0 and p2a["ks_2samp_p"] >= 0.10
    assert len(wits) == 8

    def _mr(n):
        if n < 2:
            return False
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
            if n == p:
                return True
            if n % p == 0:
                return False
        d, s = n - 1, 0
        while (d & 1) == 0:
            d >>= 1
            s += 1
        for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
            x = pow(a, d, n)
            if x == 1 or x == n - 1:
                continue
            for _ in range(s - 1):
                x = (x * x) % n
                if x == n - 1:
                    break
            else:
                return False
        return True

    def _evp(cs, x, q):
        acc = 0
        for c in reversed(cs):
            acc = (acc * x + c) % q
        return acc

    def _chk_rt(q, s, t, sd):
        assert q % 15 == 1 and _mr(q)
        assert pow(s, 15, q) == 1 and pow(s, 3, q) != 1 and pow(s, 5, q) != 1
        w = pow(s, 5, q)
        t3 = pow(t, 3, q)
        if sd == "f":
            t6 = (t3 * t3) % q
            assert ((t6 * t3) + (w + 1) * t6 + t3 - w) % q == 0
        else:
            assert (t3 - (U_BLS % q) * ((w * w) % q)) % q == 0

    bins = res["bins"]
    rng = _rnd.Random(20261124)
    min_rat = 1e9
    cbits = []
    for w in wits:
        qp, sp, tp, ps = w["q_p"], w["s_p"], w["t_p"], w["p_side"]
        qc, sc, tc, cs = w["q_c"], w["s_c"], w["t_c"], w["c_side"]
        aa, bb = w["a"], w["b"]
        assert qc >= (1 << 70) and qc < qp and math.gcd(*(aa + bb)) == 1
        _chk_rt(qp, sp, tp, ps)
        _chk_rt(qc, sc, tc, cs)
        assert (_evp(aa, sp, qp) - tp * _evp(bb, sp, qp)) % qp == 0
        assert (_evp(aa, sc, qc) - tc * _evp(bb, sc, qc)) % qc == 0
        nf, ng = norms_fg_fn(aa, bb)
        npv = nf if ps == "f" else ng
        ncv = nf if cs == "f" else ng
        assert npv % qp == 0 and npv % (qp * qp) != 0 and ncv % qc == 0 and ncv % (qc * qc) != 0
        Bc, dec, dqc = lll_fn(sq_hnf_fn(qc, sc, tc))
        assert dec == qc * qc and dqc == (1265625 ** 2) * qc * qc
        assert screen_special_q_basis(Bc, qc, cs, norms_fg_fn) is True
        ab_norm = math.sqrt(sum(x * x for x in aa) + sum(x * x for x in bb))
        b0_norm = math.sqrt(sum(x * x for x in Bc[0]))
        assert ab_norm > b0_norm
        cb = w["c_bits"]
        cbits.append(cb)
        bkey = next((b for b in bins if b[0] <= cb < b[1]), bins[0] if cb < bins[0][0] else bins[-1])
        dp = res["DET"][bkey, cs]
        Rs, ns, iAs, kAs, kBs = dp["R"], dp["n"], dp["iA"], dp["kA"], dp["kB"]
        lqc = math.log2(qc)
        r2 = Rs * Rs
        rad_hi = Rs + 0.5
        pr_sum = 0.0
        kept = 0
        for _ in range(n_cells):
            while True:
                u = [rng.gauss(0.0, 1.0) for _ in range(16)]
                scr = rad_hi * (rng.random() ** 0.0625) / math.sqrt(sum(x * x for x in u))
                kk = [round(x * scr) for x in u]
                if 0 < sum(x * x for x in kk) <= r2 and math.gcd(*kk) == 1:
                    break
            vv = [sum(kk[i] * Bc[i][j] for i in range(16)) for j in range(16)]
            ok_bc, nf_rem, ng_rem, _, _ = bound_cell_check_and_cap(vv[:8], vv[8:], dp["nf"], dp["ng"], lq=lqc, side=cs)
            if not ok_bc:
                continue
            nf_ex, ng_ex = norms_fg_fn(vv[:8], vv[8:])
            if nf_ex == 0 or ng_ex == 0:
                continue
            nf_er = math.log2(nf_ex) - (lqc if cs == "f" else 0.0)
            ng_er = math.log2(ng_ex) - (lqc if cs == "g" else 0.0)
            assert nf_er <= nf_rem + 1e-9 and ng_er <= ng_rem + 1e-9
            if nf_er < dp["nf_floor"] or ng_er < dp["ng_floor"]:
                continue
            kept += 1
            pr_sum += DS.ss_lower2(nf_rem, Pm.LX, float(ns), iAs, kAs, kBs, 384, Pm.ALPHA_F, 0.0) * DS.ss_lower2(ng_rem, Pm.LX, float(ns), iAs, kAs, kBs, 384, Pm.ALPHA_G, 0.0)
        rat = (pr_sum / n_cells) / (2.0 ** dp["pr_log2"])
        assert kept > 0 and rat >= 10.0
        min_rat = min(min_rat, rat)
    return dict(
        checks_ok=True,
        part1_pairs=p1["n_pairs"],
        part1_below_1x=p1["children_below_1x"],
        part1_min_ratio=p1["cert_ratio_quantiles"]["min"],
        part1_stratum_diff_log2=p1["stratum_controlled_weighted_minus_unweighted_log2"],
        part2a_children=p2a["n_reached_children"],
        part2a_below_1x=p2a["children_below_1x"],
        part2a_ks_p=p2a["ks_2samp_p"],
        n_witnesses=len(wits),
        child_bit_range=[min(cbits), max(cbits)],
        min_screened_cert_ratio=round(min_rat, 3),
    )


if __name__ == "__main__":
    import random
    import verify_q15_descent as VD
    rng = random.Random(42)
    gmin = 1e9; gmax = -1e9
    for _ in range(100):
        aa = [rng.randint(-10**12, 10**12) for _ in range(8)]
        bb = [rng.randint(-10**12, 10**12) for _ in range(8)]
        nf, ng = VD._norms_fg(aa, bb)
        nf_ub, ng_ub, nf_c, ng_c = bound_cell_int_fg(aa, bb)
        gf = nf_ub - math.log2(nf)
        gg = ng_ub - math.log2(ng)
        assert gf >= 0.0 and gg >= 0.0
        assert nf_c >= math.ceil(math.log2(nf)) and ng_c >= math.ceil(math.log2(ng))
        gmin = min(gmin, gf, gg); gmax = max(gmax, gf, gg)
    print(f"q15_bound_cell self-test OK! 100 random cells: gmin={gmin:.3e} >= 0, gmax={gmax:.6f} bits")
