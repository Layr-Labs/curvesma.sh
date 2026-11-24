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
