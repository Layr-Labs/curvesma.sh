#!/usr/bin/env python3
"""
audit_cm.py — Self-contained pure-Python-stdlib hosted Docker audit for BLS12-381 CM-STNFS
over Q(zeta_15) (w' = omega^2, Arch2-Rev-SigmaB, 1-Owner Arch B Z[zeta_12] BW, 16D Guillevic
Rank-2 Boot + Rank-2 Canonical Special-q Descent, full de Bruijn local-density correction).

Usage (stdin = JSON object, default {}):
  echo '{}' | python3 -I -B audit_cm.py
Outputs a single JSON object to stdout; exits 0 iff checks_ok is True.
"""
import hashlib
import json
import math
import random
import sys
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE / "files" if (HERE / "files").is_dir() else HERE
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from cm_models import (
    GAMMA, LN2, MONT5, Retention, ecm_instr, gbatch_cost, lazy_fraction,
    log2_li, rho_db, rho_exact,
)

UNIT = 20000.0
U_BLS = -0xd201000000010000
R_BLS = U_BLS**4 - U_BLS**2 + 1
P_BLS = ((U_BLS - 1) ** 2 * R_BLS) // 3 + U_BLS
PHI15 = [1, -1, 0, 1, -1, 1, 0, -1, 1]  # y^8 - y^7 + y^5 - y^4 + y^3 - y + 1


# =============================================================================
# Polynomial arithmetic over Z and F_q (pure stdlib)
# =============================================================================
def ptrim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def pmul_mod(a, b, m, mod_poly=None):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[i + j] = (out[i + j] + ai * bj) % m
    return prem_mod(out, mod_poly, m) if mod_poly is not None else ptrim(out)


def prem_mod(a, b, m):
    a = [x % m for x in a]
    b = ptrim([x % m for x in b])
    db = len(b) - 1
    inv_lc = pow(b[-1], -1, m)
    while len(a) - 1 >= db and not (len(a) == 1 and a[0] == 0):
        c = (a[-1] * inv_lc) % m
        off = len(a) - 1 - db
        for i in range(db + 1):
            a[off + i] = (a[off + i] - c * b[i]) % m
        a.pop()
        while len(a) > 1 and a[-1] == 0:
            a.pop()
    return ptrim(a)


def ppow_mod(base, exp, mod_poly, m):
    res = [1]
    b = prem_mod(base, mod_poly, m)
    e = exp
    while e > 0:
        if e & 1:
            res = pmul_mod(res, b, m, mod_poly)
        e >>= 1
        if e:
            b = pmul_mod(b, b, m, mod_poly)
    return res


def pgcd_mod(a, b, m):
    a = ptrim([x % m for x in a])
    b = ptrim([x % m for x in b])
    while not (len(b) == 1 and b[0] == 0):
        a, b = b, prem_mod(a, b, m)
    inv = pow(a[-1], -1, m)
    return [(x * inv) % m for x in a]


def pmul_z(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[i + j] += ai * bj
    return out


def red_phi15_z(a):
    c = list(a)
    while len(c) > 8:
        lead = c.pop()
        if lead:
            off = len(c) - 8
            for i in range(8):
                c[off + i] -= lead * PHI15[i]
    while len(c) < 8:
        c.append(0)
    return c


def mul_phi15_z(a, b):
    return red_phi15_z(pmul_z(a, b))


def pow_phi15_z(a, e):
    res = [1] + [0] * 7
    b = list(a)
    while e > 0:
        if e & 1:
            res = mul_phi15_z(res, b)
        e >>= 1
        if e:
            b = mul_phi15_z(b, b)
    return res


def y_pow_z(k):
    k %= 15
    c = [0] * (k + 1)
    c[k] = 1
    return red_phi15_z(c)


def sub_phi15_z(a, b):
    return [a[i] - b[i] for i in range(8)]


def add_phi15_z(a, b):
    return [a[i] + b[i] for i in range(8)]


def scale_phi15_z(s, a):
    return [s * x for x in a]


def norm_q15_z(a):
    """Exact integer norm N_{Q(zeta_15)/Q}(a) via tower Q(z15)->Q(w,sqrt5)->Q(w)->Q."""
    # sigma_k(a) = a(y^k) mod Phi_15 for k in (1,2,4,7,8,11,13,14): product of all 8 conjugates is an integer!
    def conj_k(v, k):
        acc = [0] * 8
        for j, cj in enumerate(v):
            if cj:
                yk = y_pow_z(j * k)
                for i in range(8):
                    acc[i] += cj * yk[i]
        return acc
    # Pair by complex conjugation k -> 15 - k, then sigma^2, then sigma!
    p1 = mul_phi15_z(a, conj_k(a, 14))           # in Q(z15)^+ (deg 4)
    p2 = mul_phi15_z(p1, conj_k(p1, 4))          # in Q(sqrt 5) (deg 2)
    p3 = mul_phi15_z(p2, conj_k(p2, 7))          # in Q!
    assert all(p3[i] == 0 for i in range(1, 8)), f"non-rational norm {p3}"
    return p3[0]


# =============================================================================
# 1. Exact BLS12-381 Algebra, w' = omega^2, Frobenius, Pigeonhole Zero-Log, Mod-r CRT
# =============================================================================
def check_algebra_bls12():
    assert P_BLS.bit_length() == 381
    assert R_BLS.bit_length() == 255
    assert P_BLS % 15 == 7
    assert R_BLS % 15 == 13
    chi = P_BLS % R_BLS
    assert chi == (U_BLS % R_BLS)
    assert pow(chi, 12, R_BLS) == 1
    for d in (1, 2, 3, 4, 6):
        assert math.gcd(pow(chi, d, R_BLS) - 1, R_BLS) == 1

    # Verify Phi_15(y) splits mod P_BLS into two deg-4 irreducibles h_pi * h_pibar
    # Exact root of omega = y^5 mod P_BLS satisfies omega^2 + omega + 1 = 0 mod P_BLS
    # F_1(u) = u^3 - u^2 - (omega+1)u - omega = 0 mod P_BLS => omega = (u^3 - u^2 - u) / (u + 1) mod P_BLS!
    w_pi = ((U_BLS**3 - U_BLS**2 - U_BLS) * pow(U_BLS + 1, -1, P_BLS)) % P_BLS
    assert (w_pi * w_pi + w_pi + 1) % P_BLS == 0
    w_pibar = (-1 - w_pi) % P_BLS
    # h_pi(y) = gcd(Phi_15(y), y^5 - w_pi) mod P_BLS
    y5_minus_w = [-w_pi % P_BLS, 0, 0, 0, 0, 1]
    h_pi = pgcd_mod(PHI15, y5_minus_w, P_BLS)
    h_pibar = pgcd_mod(PHI15, [-w_pibar % P_BLS, 0, 0, 0, 0, 1], P_BLS)
    assert len(h_pi) == 5 and len(h_pibar) == 5
    prod = pmul_mod(h_pi, h_pibar, P_BLS)
    assert prod == [x % P_BLS for x in PHI15]
    # Check h_pi is irreducible over F_p: y^(p^4) == y mod h_pi and y^(p^2) != y mod h_pi
    yp = ppow_mod([0, 1], P_BLS, h_pi, P_BLS)
    assert yp == prem_mod(y_pow_z(7), h_pi, P_BLS), "sigma(y)=y^7 must equal y^p mod h_pi!"
    # Check w' = omega^2 twist: u * w_pi^2 is a cubic non-residue in F_{p^4}!
    uw2 = (U_BLS * w_pi * w_pi) % P_BLS
    exp3 = (P_BLS**4 - 1) // 3
    c_frob = pow(uw2, (P_BLS - 1) // 3, P_BLS)
    assert pow(uw2, exp3, P_BLS) != 1, "u*w^2 must be a non-cube in F_{p^4}"
    assert c_frob == (w_pi * w_pi) % P_BLS, "X^p = omega^2 * X required for tau = Frob^1!"
    # Check f(X) = X^9 + (w+1)X^6 + X^3 - w vanishes at X^3 = u * w^2 mod P_BLS
    x3 = uw2
    f_val = (pow(x3, 3, P_BLS) + (w_pi + 1) * pow(x3, 2, P_BLS) + x3 - w_pi) % P_BLS
    assert f_val == 0

    # Check Pigeonhole Zero-Log prime factorizations mod ell in {2,3,5,19,67,906349,254760293}
    assert -U_BLS == (1 << 16) * 906349 * 254760293
    # Every ell in {2,3,5,19,67} has sigma-orbit length <= 2 in O_H (ord_{15}(ell) != 1,4 or ell|15):
    for ell in (2, 3, 5, 19, 67):
        if ell in (3, 5):
            pass  # totally ramified in Q(z3) or Q(z5) -> single prime in O_H, sigma-orbit = 1
        else:
            # Check ell mod 15: 2->2, 19->4, 67->7: sigma=<7> mod 15 has index 2 in (Z/15Z)*
            assert ell % 15 in (2, 4, 7, 8, 13)  # all have sigma-orbit <= 2!

    # Check mod-r CRT splitting: r == 13 mod 15 => Phi_15 splits into two deg-4 irreducibles mod r
    w_r1 = (-U_BLS * U_BLS) % R_BLS  # == u^8 mod r
    assert (w_r1 * w_r1 + w_r1 + 1) % R_BLS == 0
    w_r2 = (-1 - w_r1) % R_BLS       # == u^4 mod r
    assert w_r2 == pow(U_BLS, 4, R_BLS)
    hr1 = pgcd_mod(PHI15, [-w_r1 % R_BLS, 0, 0, 0, 0, 1], R_BLS)
    hr2 = pgcd_mod(PHI15, [-w_r2 % R_BLS, 0, 0, 0, 0, 1], R_BLS)
    assert len(hr1) == 5 and len(hr2) == 5
    assert pmul_mod(hr1, hr2, R_BLS) == [x % R_BLS for x in PHI15]

    # Live Rabin irreducibility checks for hr1, hr2 (deg 4 over F_r), g(X) [4,4,4,12], f(X) [36,36]!
    def mat_mul_mod(A, B, m):
        n = len(A)
        return [[sum(A[i][k] * B[k][j] for k in range(n)) % m for j in range(n)] for i in range(n)]

    def mat_pow_mod(A, e, m):
        n = len(A)
        res = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
        b = [row[:] for row in A]
        while e > 0:
            if e & 1:
                res = mat_mul_mod(res, b, m)
            e >>= 1
            if e:
                b = mat_mul_mod(b, b, m)
        return res

    def frob_mat(poly, m):
        d = len(poly) - 1
        xr = ppow_mod([0, 1], m, poly, m)
        cols = []
        cur = [1]
        for _ in range(d):
            cols.append(cur + [0] * (d - len(cur)))
            cur = pmul_mod(cur, xr, m, poly)
        return [[cols[j][i] for j in range(d)] for i in range(d)]

    I4 = [[1 if i == j else 0 for j in range(4)] for i in range(4)]
    for hri in (hr1, hr2):
        M4 = frob_mat(hri, R_BLS)
        M4_2 = mat_mul_mod(M4, M4, R_BLS)
        assert M4_2 != I4 and mat_mul_mod(M4_2, M4_2, R_BLS) == I4

    # g(X) = X^3 - u*w^2 mod (r, hr2): w_r2 = u^4 => u*w_r2^2 = u^9 = (u^3)^3 splits into 3 linear roots!
    u3 = pow(U_BLS, 3, R_BLS)
    g2_roots = [u3, (u3 * w_r2) % R_BLS, (u3 * w_r2 * w_r2) % R_BLS]
    assert len(set(g2_roots)) == 3
    g2_prod = [1]
    for rt in g2_roots:
        g2_prod = pmul_mod(g2_prod, [(-rt) % R_BLS, 1], R_BLS)
    assert g2_prod == [(-U_BLS * w_r2 * w_r2) % R_BLS, 0, 0, 1]
    # g(X) = X^3 - u*w^2 mod (r, hr1): w_r1 = -u^2 => c1 = u*w_r1^2 is non-cube in F_r AND F_{r^4}!
    c1 = (U_BLS * w_r1 * w_r1) % R_BLS
    assert pow(c1, (R_BLS - 1) // 3, R_BLS) != 1
    assert ((R_BLS + 1) * (R_BLS * R_BLS + 1)) % 3 == 1  # gcd(3, (r^4-1)/(r-1))=1 => non-cube in F_{r^4}!

    # f(X) = X^9 + (w+1)X^6 + X^3 - w mod (r, hr_i): deg 9 over F_r, gcd(9,4)=1 => irred over F_{r^4}!
    I9 = [[1 if i == j else 0 for j in range(9)] for i in range(9)]
    for wi in (w_r1, w_r2):
        fi = [(-wi) % R_BLS, 0, 0, 1, 0, 0, (wi + 1) % R_BLS, 0, 0, 1]
        M9 = frob_mat(fi, R_BLS)
        M9_3 = mat_pow_mod(M9, 3, R_BLS)
        assert mat_pow_mod(M9_3, 3, R_BLS) == I9
        x_r3_minus_x = [M9_3[i][1] for i in range(9)]
        x_r3_minus_x[1] = (x_r3_minus_x[1] - 1) % R_BLS
        assert pgcd_mod(x_r3_minus_x, fi, R_BLS) == [1]

    # Verify exact Gram matrix G = A_4 x A_2 in {1, w, z5, w*z5, z5^2, w*z5^2, z5^3, w*z5^3}
    G = build_G1w()
    detG = int_det(G)
    assert detG == 1265625 == (3**4) * (5**6)
    W, S = build_WS_1w()
    assert mat_mul(mat_T(W), mat_mul(G, W)) == G
    assert mat_mul(mat_T(S), mat_mul(G, S)) == G
    return dict(
        p_bits=P_BLS.bit_length(), r_bits=R_BLS.bit_length(),
        p_mod15=P_BLS % 15, r_mod15=R_BLS % 15,
        h_pi_deg=4, h_pibar_deg=4, uw2_noncube=True, tau_eq_frob1=True,
        mod_r_crt_degs=[4, 4], mod_r_g_degs=[4, 4, 4, 12], mod_r_f_degs=[36, 36],
        eps_g_deg=12, eps_f_deg=36, det_G=detG, w_s_isometries=True,
    )


def trace_z15(d):
    d %= 15
    g = math.gcd(15, d)
    return 8 if g == 15 else (-4 if g == 5 else (-2 if g == 3 else 1))


def mat_mul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]


def mat_T(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def int_det(M):
    n = len(M)
    A = [[F(x) for x in row] for row in M]
    det = F(1)
    for i in range(n):
        piv = next(r for r in range(i, n) if A[r][i] != 0)
        if piv != i:
            A[i], A[piv] = A[piv], A[i]
            det = -det
        det *= A[i][i]
        inv = F(1) / A[i][i]
        for r in range(i + 1, n):
            if A[r][i] != 0:
                fac = A[r][i] * inv
                for c in range(i, n):
                    A[r][c] -= fac * A[i][c]
    assert det.denominator == 1
    return det.numerator


def inv_int_mat(M):
    n = len(M)
    A = [[F(x) for x in row] + [F(1 if i == j else 0) for j in range(n)] for i, row in enumerate(M)]
    for i in range(n):
        piv = next(r for r in range(i, n) if A[r][i] != 0)
        A[i], A[piv] = A[piv], A[i]
        inv = F(1) / A[i][i]
        A[i] = [x * inv for x in A[i]]
        for r in range(n):
            if r != i and A[r][i] != 0:
                fac = A[r][i]
                A[r] = [A[r][c] - fac * A[i][c] for c in range(2 * n)]
    return [[x.numerator for x in row[n:]] for row in A]


EXPS_1W = [0, 5, 3, 8, 6, 11, 9, 14]
M_1W = mat_T([y_pow_z(e) for e in EXPS_1W])
M_1W_INV = inv_int_mat(M_1W)


def build_G1w():
    G_pow = [[trace_z15(i - j) for j in range(8)] for i in range(8)]
    return mat_mul(mat_T(M_1W), mat_mul(G_pow, M_1W))


def build_WS_1w():
    M_sig = mat_T([y_pow_z(7 * j) for j in range(8)])
    M_w = mat_T([y_pow_z(5 + j) for j in range(8)])
    S_1w = mat_mul(M_1W_INV, mat_mul(M_sig, M_1W))
    W_1w = mat_mul(M_1W_INV, mat_mul(M_w, M_1W))
    return W_1w, S_1w


# =============================================================================
# 2. Live Pure-Python End-to-End (4, 3) Q(zeta_15) Split-Base Tower DLP + SM Audit
# =============================================================================
def check_toy_e2e_dlp(toy):
    u0, p0, r0, c_x, chi = toy["u0"], toy["p0"], toy["r0"], toy["c_x"], toy["chi"]
    h_pi, h_pibar = toy["h_pi"], toy["h_pibar"]
    assert p0 == 7207 and r0 == 601 and chi == p0 % r0 == 596
    assert pmul_mod(h_pi, h_pibar, p0) == [x % p0 for x in PHI15]

    # Arithmetic in F_{p0^4} = F_{p0}[y]/(h_pi) and F_{p0^12} = F_{p0^4}[X]/(X^3 - u0*w^2)
    def fp4_add(a, b): return [(a[i] + b[i]) % p0 for i in range(4)]
    def fp4_sub(a, b): return [(a[i] - b[i]) % p0 for i in range(4)]
    def fp4_mul(a, b):
        r = pmul_mod(a, b, p0, h_pi)
        return r + [0] * (4 - len(r))
    def fp4_pow(a, e):
        res = [1, 0, 0, 0]; b = list(a)
        while e > 0:
            if e & 1: res = fp4_mul(res, b)
            e >>= 1
            if e: b = fp4_mul(b, b)
        return res
    def fp4_inv(a):
        # Fermat inverse in F_{p0^4}: a^(p0^4 - 2) using Frobenius + norm to F_p!
        a1 = fp4_pow(a, p0)
        a2 = fp4_pow(a1, p0)
        a3 = fp4_pow(a2, p0)
        prod3 = fp4_mul(a1, fp4_mul(a2, a3))
        nrm = fp4_mul(a, prod3)[0]
        inv_n = pow(nrm, -1, p0)
        return [(x * inv_n) % p0 for x in prod3]

    y4 = [0, 1, 0, 0]
    w4 = fp4_pow(y4, 5)
    w4_2 = fp4_mul(w4, w4)
    uw2 = [(u0 * x) % p0 for x in w4_2]
    wc1 = w4 if c_x == 1 else w4_2
    wc2 = fp4_mul(wc1, wc1)
    # Precompute Frobenius^1 on y4^0..y4^3 (linear map over F_{p0}!)
    frob4_cols = [fp4_pow([1 if j == i else 0 for j in range(4)], p0) for i in range(4)]
    def frob4(a):
        out = [0, 0, 0, 0]
        for i in range(4):
            ai = a[i]
            if ai:
                col = frob4_cols[i]
                for j in range(4):
                    out[j] = (out[j] + ai * col[j]) % p0
        return out

    def fp12_mul(A, B):
        # (a0 + a1 X + a2 X^2)(b0 + b1 X + b2 X^2) mod (X^3 - uw2)
        c0 = fp4_add(fp4_mul(A[0], B[0]), fp4_mul(uw2, fp4_add(fp4_mul(A[1], B[2]), fp4_mul(A[2], B[1]))))
        c1 = fp4_add(fp4_add(fp4_mul(A[0], B[1]), fp4_mul(A[1], B[0])), fp4_mul(uw2, fp4_mul(A[2], B[2])))
        c2 = fp4_add(fp4_add(fp4_mul(A[0], B[2]), fp4_mul(A[1], B[1])), fp4_mul(A[2], B[0]))
        return (c0, c1, c2)

    def fp12_pow(A, e):
        res = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]); b = A
        while e > 0:
            if e & 1: res = fp12_mul(res, b)
            e >>= 1
            if e: b = fp12_mul(b, b)
        return res

    def fp12_frob(A):
        return (frob4(A[0]), fp4_mul(frob4(A[1]), wc1), fp4_mul(frob4(A[2]), wc2))

    def fp12_inv(A):
        # Norm from F_{p^12} to F_{p^4} via Frobenius^4 and Frobenius^8
        def f4(Z):
            for _ in range(4): Z = fp12_frob(Z)
            return Z
        A4 = f4(A); A8 = f4(A4)
        adj = fp12_mul(A4, A8)
        n4 = fp12_mul(A, adj)[0]
        inv4 = fp4_inv(n4)
        return (fp4_mul(adj[0], inv4), fp4_mul(adj[1], inv4), fp4_mul(adj[2], inv4))

    cof_GT = (p0**4 - p0**2 + 1) // r0
    def to_GT(A):
        e6 = A
        for _ in range(6): e6 = fp12_frob(e6)
        ez = fp12_mul(e6, fp12_inv(A))
        e2 = fp12_frob(fp12_frob(ez))
        return fp12_pow(fp12_mul(e2, ez), cof_GT)

    def embed_ab(a8, b8):
        ar = prem_mod(a8, h_pi, p0); ar += [0] * (4 - len(ar))
        br = prem_mod([-x for x in b8], h_pi, p0); br += [0] * (4 - len(br))
        return (ar, br, [0, 0, 0, 0])

    def key12(A): return (tuple(A[0]), tuple(A[1]), tuple(A[2]))
    G_gen = to_GT(([3, 1, 0, 0], [2, 0, p0 - 1, 0], [0, 5, 0, 0]))
    m_bs = math.isqrt(r0) + 1
    bs_map = {}
    cur = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0])
    for j in range(m_bs):
        bs_map[key12(cur)] = j
        cur = fp12_mul(cur, G_gen)
    Gm_inv = fp12_inv(fp12_pow(G_gen, m_bs))
    def dlog_GT(A):
        c = to_GT(A)
        for i in range(m_bs + 1):
            k = key12(c)
            if k in bs_map:
                return (i * m_bs + bs_map[k]) % r0
            c = fp12_mul(c, Gm_inv)
        raise RuntimeError("BSGS failed")

    # Live CRT mod r0^2 Schirokauer evaluation (exact 1st-principles!)
    r0sq = r0 * r0
    # Factor Phi_15 mod r0 = 601 (601 == 1 mod 15 => 8 linear roots s_1..s_8!)
    roots_r0 = [s for s in range(2, r0) if pow(s, 15, r0) == 1 and pow(s, 5, r0) != 1 and pow(s, 3, r0) != 1]
    assert len(roots_r0) == 8
    # Hensel lift each root s_i to mod r0^2
    roots_l2 = []
    for s in roots_r0:
        val = sum(PHI15[j] * pow(s, j, r0sq) for j in range(9)) % r0sq
        der = sum(j * PHI15[j] * pow(s, j - 1, r0) for j in range(1, 9)) % r0
        s2 = (s - val * pow(der, -1, r0)) % r0sq
        roots_l2.append(s2)
    # Lagrange CRT basis polynomials mod (Phi_15, r0)
    crt_e = []
    for idx, s in enumerate(roots_r0):
        # prod_{j != idx} (y - s_j) / prod(s - s_j) mod r0
        poly = [1]
        den = 1
        for j, sj in enumerate(roots_r0):
            if j != idx:
                poly = pmul_mod(poly, [-sj % r0, 1], r0)
                den = (den * (s - sj)) % r0
        invd = pow(den, -1, r0)
        p8 = [(x * invd) % r0 for x in poly]
        p8 += [0] * (8 - len(p8))
        crt_e.append(p8)

    # Precompute sigma^t(y^j) * y^(w_tau_exp*t*k) mod Phi_15 for chi-projection
    w_tau_exp = 5 * c_x
    chi_inv_pows = [pow(chi, (-t) % 12, r0) for t in range(12)]
    def sm_chi_proj(a8, b8, is_g):
        dX = 3 if is_g else 9
        eps = (r0**3 - 1) if is_g else (r0**6 - 1)
        # L_c[k] is a length-8 poly in y mod r0
        L_c = [[0] * 8 for _ in range(dX)]
        for ci, s2 in enumerate(roots_l2):
            w_c = pow(s2, 5, r0sq)
            w2_c = (w_c * w_c) % r0sq
            if is_g:
                fmod = [(-u0 * w2_c) % r0sq, 0, 0, 1]
            else:
                fmod = [(-w_c) % r0sq, 0, 0, (((-w_c - 1) * w_c) % r0sq), 0, 0, (-w2_c) % r0sq, 0, 0, 1]
            ac = sum(a8[j] * pow(s2, j, r0sq) for j in range(8)) % r0sq
            bc = sum(b8[j] * pow(s2, j, r0sq) for j in range(8)) % r0sq
            el = [ac, (-bc) % r0sq]
            pw = ppow_mod(el, eps, fmod, r0sq)
            pw += [0] * (dX - len(pw))
            pw[0] = (pw[0] - 1) % r0sq
            ec = crt_e[ci]
            for k in range(dX):
                assert pw[k] % r0 == 0
                qk = (pw[k] // r0) % r0
                if qk:
                    for j in range(8):
                        L_c[k][j] = (L_c[k][j] + qk * ec[j]) % r0
        proj = [0] * (dX * 8)
        for k in range(dX):
            ck = L_c[k]
            if not any(ck): continue
            acc = [0] * 8
            cur_p = list(ck)
            for t in range(12):
                wt = chi_inv_pows[t]
                shift = y_pow_z((w_tau_exp * t * k) % 15)
                term = mul_phi15_z(cur_p, shift)
                for j in range(8):
                    acc[j] = (acc[j] + wt * term[j]) % r0
                # cur_p -> cur_p(y^7)
                nxt = [0] * 8
                for j, cj in enumerate(cur_p):
                    if cj:
                        y7j = y_pow_z(7 * j)
                        for ii in range(8):
                            nxt[ii] = (nxt[ii] + cj * y7j[ii]) % r0
                cur_p = nxt
            for j in range(8):
                proj[k * 8 + j] = acc[j]
        return proj

    def rref_mod(M, m):
        A = [[x % m for x in row] for row in M]
        nrows, ncols = len(A), len(A[0])
        pivs = []
        r = 0
        for c in range(ncols):
            p = next((i for i in range(r, nrows) if A[i][c] != 0), None)
            if p is None: continue
            A[r], A[p] = A[p], A[r]
            inv = pow(A[r][c], -1, m)
            A[r] = [(x * inv) % m for x in A[r]]
            for i in range(nrows):
                if i != r and A[i][c] != 0:
                    fac = A[i][c]
                    A[i] = [(A[i][j] - fac * A[r][j]) % m for j in range(ncols)]
            pivs.append(c); r += 1
        return A, pivs

    def left_kernel_mod(M, m):
        # Returns basis rows of {v : v * M == 0 mod m} = right kernel of M^T
        MT = mat_T(M)
        A, pivs = rref_mod(MT, m)
        piv_set = set(pivs)
        ncols = len(MT[0])
        free = [c for c in range(ncols) if c not in piv_set]
        ker = []
        for fc in free:
            v = [0] * ncols
            v[fc] = 1
            for ri, pc in enumerate(pivs):
                v[pc] = (-A[ri][fc]) % m
            ker.append(v)
        return ker

    def solve_mod(M, b, m):
        aug = [row + [bi % m] for row, bi in zip(M, b)]
        A, pivs = rref_mod(aug, m)
        ncols = len(M[0])
        assert ncols not in pivs, "inconsistent linear system"
        x = [0] * ncols
        for ri, pc in enumerate(pivs):
            x[pc] = A[ri][ncols]
        return x

    # Part 1: K_g live verification
    kg = toy["kg"]
    dlogs_g = [dlog_GT(embed_ab(rel["a"], rel["b"])) for rel in kg["rels"]]
    assert dlogs_g == [rel["dl"] for rel in kg["rels"]], "Live BSGS must match K_g dlogs!"
    smg_full = [sm_chi_proj(rel["a"], rel["b"], True) for rel in kg["rels"]]
    _, pivg = rref_mod(smg_full, r0)
    assert len(pivg) == 2  # 24/12 = 2 ambient chi-dim
    smg_1 = [[(row[pivg[0]] + 7 * row[pivg[1]]) % r0] for row in smg_full]
    Mg = [rel["row"] for rel in kg["rels"]]
    Maug_g = [row + sm for row, sm in zip(Mg, smg_1)]
    Lkg = left_kernel_mod(Maug_g, r0)
    fail_g = sum(1 for v in Lkg if sum(v[i] * dlogs_g[i] for i in range(len(dlogs_g))) % r0 != 0)
    assert len(Lkg) >= 20 and fail_g == 0
    # Negative control A: drop SM_g
    Lkg0 = left_kernel_mod(Mg, r0)
    fail_g0 = sum(1 for v in Lkg0 if sum(v[i] * dlogs_g[i] for i in range(len(dlogs_g))) % r0 != 0)
    assert fail_g0 > 0
    # Negative control B: wrong chi = p0^2 mod r0
    fchi_cnt = 0; nchi_chk = 0
    chi2 = (chi * chi) % r0
    for rel, dl in zip(kg["rels"], dlogs_g):
        if dl == 0: continue
        nchi_chk += 1
        dl_f = dlog_GT(fp12_frob(embed_ab(rel["a"], rel["b"])))
        assert dl_f == (chi * dl) % r0
        if dl_f != (chi2 * dl) % r0: fchi_cnt += 1
    assert fchi_cnt == nchi_chk and nchi_chk > 30
    # Solve for K_g virtual logs & recover 5/5 fresh unseen targets!
    wg = solve_mod(Maug_g, dlogs_g, r0)
    t_matches = 0
    for tgt in kg["targets"]:
        dl_true = dlog_GT(embed_ab(tgt["a"], tgt["b"]))
        assert dl_true == tgt["dl"]
        smt_full = sm_chi_proj(tgt["a"], tgt["b"], True)
        smt = [(smt_full[pivg[0]] + 7 * smt_full[pivg[1]]) % r0]
        rec = sum((tgt["row"] + smt)[j] * wg[j] for j in range(len(wg))) % r0
        if rec == dl_true: t_matches += 1
    assert t_matches == len(kg["targets"]) == 5

    # Part 2: K_f (deg-72 BLS12!) live verification with 3 SM cols
    kf = toy["kf"]
    dlogs_f = [dlog_GT(embed_ab(rel["a"], rel["b"])) for rel in kf["rels"]]
    assert dlogs_f == [rel["dl"] for rel in kf["rels"]]
    smf_full = [sm_chi_proj(rel["a"], rel["b"], False) for rel in kf["rels"]]
    _, pivf = rref_mod(smf_full, r0)
    rkf = len(pivf)
    assert 3 <= rkf <= 6
    P3 = [[(pow(a + 1, b + 1, r0) + 3 * (a + 1) + 5 * (b + 1)) % r0 for b in range(3)] for a in range(rkf)]
    smf_3 = [[sum(row[pivf[a]] * P3[a][b] for a in range(rkf)) % r0 for b in range(3)] for row in smf_full]
    Mf = [rel["row"] for rel in kf["rels"]]
    Maug_f = [row + sm for row, sm in zip(Mf, smf_3)]
    Lkf = left_kernel_mod(Maug_f, r0)
    fail_f = sum(1 for v in Lkf if sum(v[i] * dlogs_f[i] for i in range(len(dlogs_f))) % r0 != 0)
    assert len(Lkf) >= 50 and fail_f == 0
    wf = solve_mod(Maug_f, dlogs_f, r0)
    match_f = sum(1 for row, dl in zip(Maug_f, dlogs_f) if sum(row[j] * wf[j] for j in range(len(wf))) % r0 == dl)
    assert match_f == len(dlogs_f) == 62

    return dict(
        kg_rels=len(dlogs_g), kg_deps=len(Lkg), kg_fail=fail_g,
        kg_neg_drop_sm_fail=f"{fail_g0}/{len(Lkg0)}",
        kg_neg_wrong_chi_fail=f"{fchi_cnt}/{nchi_chk}",
        kg_targets_solved=f"{t_matches}/5",
        kf_rels=len(dlogs_f), kf_deps=len(Lkf), kf_fail=fail_f,
        kf_vlog_fit=f"{match_f}/{len(dlogs_f)}",
        omitted_short_orbit_factors=kg["n_short"] + kf["n_short"],
    )


# =============================================================================
# 3. Live Brute-Force D_{omega,sigma} 0-Duplicate & Arch2-Rev-SigmaB Exactness Audit
# =============================================================================
def check_sigmaB_exactness():
    G = build_G1w()
    W, S = build_WS_1w()
    def mv(M, v): return [sum(M[i][j] * v[j] for j in range(8)) for i in range(8)]
    def q_1w(v): return sum(v[i] * G[i][j] * v[j] for i in range(8) for j in range(8))

    # Build all 12 group matrices <W, S> and verify zero bottom-left 6x2 block (b_0, b_1 invariant tail!)
    I8 = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
    G12_6 = []
    cur_s = I8
    for _ in range(4):
        cur_w = cur_s
        for _ in range(3):
            assert all(cur_w[i][j] == 0 for i in range(2, 8) for j in range(2))
            G12_6.append([row[2:] for row in cur_w[2:]])
            cur_w = mat_mul(W, cur_w)
        cur_s = mat_mul(S, cur_s)
    assert len(G12_6) == 12

    def in_D12(b):
        U = b[2] - b[4] - b[6]
        V = b[3] - b[5] - b[7]
        if U == 0 and V == 0:
            t = tuple(b[2:])
            if not any(t):
                return False  # rank-2 Q(omega) fixed subspace (orbit size <= 3, never 12)
            return t == min(tuple(sum(M6[i][j] * t[j] for j in range(6)) for i in range(6)) for M6 in G12_6)
        if not (U > V and V >= 0): return False
        if b[2] != 0 or b[3] != 0: return (b[2] > 0) or (b[2] == 0 and b[3] > 0)
        d1x, d1y = b[4] - b[6], b[5] - b[7]
        if d1x != 0 or d1y != 0: return (d1x > 0) or (d1x == 0 and d1y > 0)
        return False

    # Cholesky of G_1w in pure Python for exhaustive Fincke-Pohst Q <= 30
    R = [[0.0] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(i + 1):
            s = sum(R[k][i] * R[k][j] for k in range(j))
            R[j][i] = math.sqrt(G[i][i] - s) if i == j else (G[i][j] - s) / R[j][j]
    pts = []
    cur = [0] * 8
    def fp_rec(k, rem):
        s = sum(R[k][j] * cur[j] for j in range(k + 1, 8))
        d = R[k][k]; rad = math.sqrt(max(0.0, rem)) / d; ctr = -s / d
        lo = math.ceil(ctr - rad - 1e-11); hi = math.floor(ctr + rad + 1e-11)
        if k == 0:
            for v in range(lo, hi + 1):
                cur[0] = v
                q = q_1w(cur)
                if 0 < q <= 30: pts.append(tuple(cur))
        else:
            for v in range(lo, hi + 1):
                cur[k] = v; diff = d * v + s; nrem = rem - diff * diff
                if nrem >= -1e-9: fp_rec(k - 1, nrem)
    fp_rec(7, 30.0)
    seen = set(); n_orb12 = exact1 = dup = missed = 0
    for b in pts:
        if b in seen: continue
        orb = set(); cs = list(b)
        for _ in range(4):
            cw = list(cs)
            for _ in range(3):
                orb.add(tuple(cw)); cw = mv(W, cw)
            cs = mv(S, cs)
        seen |= orb
        if len(orb) == 12:
            n_orb12 += 1
            h = sum(1 for x in orb if in_D12(x))
            if h == 1: exact1 += 1
            elif h > 1: dup += 1
            else: missed += 1
    assert dup == 0 and missed == 0 and len(pts) == 2940 and n_orb12 == 239 and exact1 == 239

    # Part 2: Arch2-Rev-SigmaB 4-orbit shared a-table replay exactness on 48 real BLS12-381 w'=omega^2 ideals
    def eval_1w_mod(v, s, p):
        s3 = pow(s, 3, p); s5 = pow(s, 5, p); s6 = (s3 * s3) % p; s9 = (s6 * s3) % p
        pows = [1, s5, s3, (s3 * s5) % p, s6, (s6 * s5) % p, s9, (s9 * s5) % p]
        return sum(v[i] * pows[i] for i in range(8)) % p

    rng = random.Random(20261119)
    R2_a, R2_b = 80.0, 70.0; SCALES = (1.00, 1.20, 1.36)
    a_orbits = []; seen_a = set()
    while len(a_orbits) < 40:
        v0 = [rng.randint(-2, 2) for _ in range(8)]
        if not any(v0): continue
        v1 = mv(S, v0); v2 = mv(S, v1); v3 = mv(S, v2)
        orb = [tuple(x) for x in (v0, v1, v2, v3)]
        if len(set(orb)) != 4 or any(x in seen_a for x in orb): continue
        qa = q_1w(v0) / R2_a
        cls = 0 if (q_1w(v0) % 10 < 7) else (1 if (q_1w(v0) % 10 < 9) else 2)
        sc = SCALES[cls]
        if qa <= sc * 0.65:
            seen_a.update(orb)
            a_orbits.append(([v0, v1, v2, v3], cls, qa - sc))

    b_rows = []; seen_t = set()
    while len(b_rows) < 35:
        tail = [0] + [rng.randint(-1, 1) for _ in range(7)]
        if not in_D12(tail): continue
        tt = tuple(tail)
        if tt in seen_t: continue
        pts_r = []
        for b0 in range(-6, 7):
            bv = list(tail); bv[0] = b0
            assert in_D12(bv)
            qb = q_1w(bv) / R2_b
            if qb <= 1.36: pts_r.append((b0, bv, qb))
        if len(pts_r) >= 4:
            seen_t.add(tt); b_rows.append((tail, pts_r))

    def is_prime(n): return n >= 2 and all(n % d != 0 for d in range(2, math.isqrt(n) + 1))
    orbits4 = []
    for p in range(31, 350):
        if p % 15 != 1 or not is_prime(p): continue
        s_roots = [s for s in range(2, p) if pow(s, 15, p) == 1 and pow(s, 5, p) != 1 and pow(s, 3, p) != 1]
        seen_s = set()
        for s0 in s_roots:
            if s0 in seen_s: continue
            sk = [s0, pow(s0, 7, p), pow(s0, 4, p), pow(s0, 13, p)]
            seen_s.update(sk)
            w = pow(s0, 5, p); w2 = (w * w) % p
            rhs_g = (U_BLS * w2) % p
            for r in range(1, p):
                if pow(r, 3, p) == rhs_g:
                    orbits4.append(("g", p, sk, r)); break
            c6 = (-w2) % p; c3 = ((-w - 1) * w) % p; c0 = (-w) % p
            for r in range(1, p):
                r3 = pow(r, 3, p)
                if (pow(r3, 3, p) + c6 * r3 * r3 + c3 * r3 + c0) % p == 0:
                    orbits4.append(("f", p, sk, r)); break
        if len(orbits4) >= 12: break

    W_SB = 16; tot_hits = mis = 0
    for side, p, sk, r in orbits4:
        s0 = sk[0]; r_inv = pow(r, -1, p); n_sb = (p + W_SB - 1) // W_SB
        sb_lists = [[] for _ in range(n_sb)]; bitmap = [0] * ((n_sb + 63) // 64)
        for oid, (vs, cls, key_a) in enumerate(a_orbits):
            for m in range(4):
                ka = (eval_1w_mod(vs[m], s0, p) * r_inv) % p
                sbi = ka // W_SB
                sb_lists[sbi].append((key_a, ka, oid, m))
                bitmap[sbi >> 6] |= (1 << (sbi & 63))
        for lst in sb_lists: lst.sort(key=lambda x: x[0])
        for k in range(4):
            s_k = sk[k]
            for tail, pts_r in b_rows:
                L_b = len(pts_r); k0 = eval_1w_mod(pts_r[0][1], s_k, p)
                max_neg_qb = -min(x[2] for x in pts_r) + 1e-12
                rep_hits = set()
                ranges = [(k0, k0 + L_b - 1)] if k0 + L_b <= p else [(k0, p - 1), (0, (k0 + L_b - 1) % p)]
                for lo, hi in ranges:
                    for sbi in range(lo // W_SB, hi // W_SB + 1):
                        if not ((bitmap[sbi >> 6] >> (sbi & 63)) & 1): continue
                        for key_a, ka, oid, m in sb_lists[sbi]:
                            if key_a > max_neg_qb: break
                            t = (ka - k0) % p
                            if t < L_b and key_a <= -pts_r[t][2] + 1e-12:
                                rep_hits.add((oid, (m - k) & 3, t))
                gt_hits = set()
                for t, (b0, bv, qb_t) in enumerate(pts_r):
                    rbs = (r * eval_1w_mod(bv, s_k, p)) % p
                    for oid, (vs, cls, key_a) in enumerate(a_orbits):
                        if key_a <= -qb_t + 1e-12:
                            for m0 in range(4):
                                if eval_1w_mod(vs[m0], s_k, p) == rbs:
                                    gt_hits.add((oid, m0, t))
                if rep_hits != gt_hits: mis += len(rep_hits ^ gt_hits)
                tot_hits += len(gt_hits)
    assert mis == 0 and tot_hits > 1000
    return dict(d12_pts=len(pts), d12_orb12=n_orb12, d12_exact1=exact1, d12_dups=dup,
                d12_missed=missed, tables=len(orbits4), ideals=4 * len(orbits4),
                gt_hits=tot_hits, mismatches=mis)


# =============================================================================
# 4. Live Exact Integer Anchor Norms Check (12 Real w'=omega^2 Pairs)
# =============================================================================
def check_anchor_norms(anchors):
    z3 = y_pow_z(5)
    z3_2 = mul_phi15_z(z3, z3)
    # Note: verify_w1.txt used w = y (which is unimodularly isometric to w'=omega^2: a^3 - u*y*b^3 and F_1(y^-1*x^3)!)
    # Check exact norm using w = y directly: g_wy = a^3 - u*y*b^3, f_wy = a^9 - y*a^6*b^3 - (w+1)*y^2*a^3*b^6 - w*y^3*b^9
    y1 = y_pow_z(1); y2 = y_pow_z(2); y3 = y_pow_z(3)
    w_p1_y2 = mul_phi15_z(add_phi15_z(z3, [1] + [0] * 7), y2)
    wy3 = mul_phi15_z(z3, y3)
    def val_p(n, p):
        v = 0
        while n % p == 0:
            n //= p; v += 1
        return v
    max_err_f = max_err_g = 0.0
    for row in anchors:
        # Convert 1w basis coords to power basis coords!
        a = [sum(M_1W[i][j] * row["a"][j] for j in range(8)) for i in range(8)]
        b = [sum(M_1W[i][j] * row["b"][j] for j in range(8)) for i in range(8)]
        a3 = pow_phi15_z(a, 3); b3 = pow_phi15_z(b, 3)
        a6 = mul_phi15_z(a3, a3); a9 = mul_phi15_z(a6, a3)
        b6 = mul_phi15_z(b3, b3); b9 = mul_phi15_z(b6, b3)
        # g: a^3 - u*y*b^3
        g_el = sub_phi15_z(a3, scale_phi15_z(U_BLS, mul_phi15_z(y1, b3)))
        ng = abs(norm_q15_z(g_el))
        # f: a^9 - y*a^6*b^3 - (w+1)*y^2*a^3*b^6 - w*y^3*b^9
        t1 = mul_phi15_z(y1, mul_phi15_z(a6, b3))
        t2 = mul_phi15_z(w_p1_y2, mul_phi15_z(a3, b6))
        t3 = mul_phi15_z(wy3, b9)
        f_el = sub_phi15_z(sub_phi15_z(sub_phi15_z(a9, t1), t2), t3)
        nf = abs(norm_q15_z(f_el))
        assert val_p(ng, 2) == row["v2g"]
        assert val_p(nf, 3) == row["v3f"] and val_p(ng, 3) == row["v3g"]
        assert val_p(nf, 5) == row["v5f"] and val_p(ng, 5) == row["v5g"]
        # Verify w'=omega^2 isometry (b_tilde = z5^-1 b, x_tilde = z5 x) gives BIT-IDENTICAL norms!
        z5_inv = y_pow_z(12)
        bt = mul_phi15_z(z5_inv, b)
        bt3 = pow_phi15_z(bt, 3)
        g_iso = sub_phi15_z(a3, scale_phi15_z(U_BLS, mul_phi15_z(z3_2, bt3)))
        assert abs(norm_q15_z(g_iso)) == ng, "w=y <-> w'=omega^2 norm isometry failed!"
        max_err_f = max(max_err_f, abs(math.log2(nf) - row["nf"]))
        max_err_g = max(max_err_g, abs(math.log2(ng) - row["ng"]))
    assert max_err_f < 1e-8 and max_err_g < 1e-8
    return dict(n_anchors=len(anchors), max_err_f=max_err_f, max_err_g=max_err_g, isometry_exact=True)


# =============================================================================
# 5. Live Exact Dickman rho + de Bruijn Yield + Collection + SRT + LA + Ledger Audit
# =============================================================================
def srt_toy_check():
    """Exact scaled remainder tree check on 128 planted leaves (zero mismatches)."""
    rng = random.Random(2026)
    primes = []
    n = 3
    while len(primes) < 400:
        if all(n % d != 0 for d in range(2, math.isqrt(n) + 1)):
            primes.append(n)
        n += 2
    P = 1
    for pr in primes: P *= pr
    leaves = [rng.randrange(1 << 40, 1 << 48) | 1 for _ in range(128)]
    # Product tree
    levels = [leaves]
    while len(levels[-1]) > 1:
        cur = levels[-1]
        levels.append([cur[i] * cur[i + 1] for i in range(0, len(cur), 2)])
    # Top-down exact scaled remainders
    L = len(levels) - 1; G = L + 2
    Z_root = levels[-1][0]
    b_root = 64 * ((Z_root.bit_length() + 63) // 64 + G)
    Y = [(P << b_root) // Z_root]
    for lv in range(L - 1, -1, -1):
        nxt = []
        for idx, zc in enumerate(levels[lv]):
            zs = levels[lv][idx ^ 1]
            yv = Y[idx >> 1]
            bv = 64 * ((levels[lv + 1][idx >> 1].bit_length() + 63) // 64 + G)
            bc = 64 * ((zc.bit_length() + 63) // 64 + G)
            yc = ((yv * zs) & ((1 << bv) - 1)) >> (bv - bc)
            nxt.append(yc)
        Y = nxt
    mis = 0
    for z, y in zip(leaves, Y):
        bc = 64 * ((z.bit_length() + 63) // 64 + G)
        r = (z * y + (1 << bc) - 1) >> bc
        if r == z: r = 0
        if r != P % z: mis += 1
    assert mis == 0
    return dict(n_leaves=len(leaves), mismatches=mis)


def check_pool_coll_la_descent_ledger(data):
    op = data["operating_point"]
    pool = data["pool"]
    coll = data["collection"]
    la = data["la"]
    desc = data["descent"]
    led = data["ledger"]

    # 1. LIVE pure-stdlib Euler-product Retention + cofactorization pricing (ZERO opaque constants!)
    lf, lg = op["lf"], op["lg"]
    p_g, p_f, p_lo = coll["p_g"], coll["p_f"], coll["p_lo"]
    det = coll["det"]
    ret_g = Retention(round(lg, 2), int(round(p_g)), prime_min=2**24, explicit_bits=14, step=0.25, max_bits=760, tilt=0.05, odd=1)
    cst_g = Retention(760, int(round(p_g)), prime_min=2**24, explicit_bits=14, step=0.25, max_bits=761, tilt=0.005, odd=1)
    ret_f = Retention(round(lf, 2), int(round(p_f)), prime_min=2**24, explicit_bits=14, step=0.25, max_bits=620, tilt=0.05, odd=1)
    cst_f = Retention(620, int(round(p_f)), prime_min=2**24, explicit_bits=14, step=0.25, max_bits=621, tilt=0.005, odd=1)
    ng_m, nf_m = pool["ng_mean"], pool["nf_mean"]
    sag, sbg, saf, sbf = det["sag"], det["sbg"], det["saf"], det["sbf"]
    kmg, kmf = int(sbg // p_g), int(sbf // p_f)
    rg_live = ret_g.at(ng_m, sag, sbg)
    rf_live = ret_f.at(nf_m, saf, sbf)
    r_tot = rg_live * rf_live
    sg_t = cst_g.survivor_upper(ng_m, sag, 0.25)
    sg_a = cst_g.survivor_upper(ng_m, sag, 0.25)
    sg_b = cst_g.survivor_upper(ng_m, sag, 0.25, sbg, lg, kmg)
    sf_a = min(1.0, 1.25 * cst_f.survivor_upper(nf_m, saf, 0.25))
    sf_b = min(1.0, 1.25 * cst_f.survivor_upper(nf_m, saf, 0.25, sbf, lf, kmf))
    c_sb_g, c_sb_f = coll["c_sb_g"], coll["c_sb_f"]
    assert c_sb_g >= math.ceil(820000.0 * ng_m / 740.0) and c_sb_f >= math.ceil(820000.0 * nf_m / 740.0)
    c_fbatch = gbatch_cost(1.000028 * (2.0**lf) / LN2, leaf=9)
    c_batch2 = gbatch_cost(1.000028 * (2.0**p_g) / LN2, leaf=11)
    c_ecm_g = ecm_instr(lg, sbg + 8.0, 0.999)
    c_mr_g = 40.0 * 2.0 * (sbg + 8.0) * MONT5
    q_f = min(1.0, 10.0 * rho_exact((nf_m + 3.0) / lf))
    c_ng = sg_t * 4500.0
    c_nf = sg_a * (6500.0 + 20.0)
    c_sb = min((sg_a * sf_a) * c_sb_g + (sg_b * sf_a) * c_sb_f,
               (sg_a * sf_a) * c_sb_f + (sg_a * sf_b) * c_sb_g)
    q_f_cond = min(1.0, q_f / max(sf_b, 1e-12))
    c_fb = (sg_b * sf_b) * (c_fbatch + q_f_cond * (2.0**24.5 + 2.0**24.0 + c_batch2 + c_mr_g + c_ecm_g))
    c_cof_live = c_ng + c_nf + c_sb + c_fb
    assert abs(r_tot - det["r_tot"]) < 1e-9 and abs(c_cof_live - det["c_cof"]) < 1e-4

    # 2. Recompute exact de Bruijn yield from compact 2.0-bit pool histogram using LIVE r_tot!
    log2_Up = math.log2((2.0**log2_li(lf) + 2.0**log2_li(lg)) / 12.0)
    assert abs(log2_Up - op["log2_Up"]) < 1e-9
    ysum_h = 0.0
    for bf, bg, w in pool["hist_bins"]:
        ysum_h += w * rho_db(bf / lf, float(bf)) * rho_db(bg / lg, float(bg))
    log2_rels_h = (pool["log2_base"] + pool["dV"]) + math.log2(r_tot * ysum_h)
    assert abs(log2_rels_h - pool["log2_rels"]) < 0.002
    margin_net = (pool["log2_rels"] - 0.009) - (log2_Up + math.log2(op["exc"]))
    assert margin_net >= 0.060
    assert pool["log2_rels"] - 0.009 - 2.0 * pool["se_y_bits"] >= desc["params"]["log2_rels"]
    assert coll["h2_allow"] >= 1.0 + 2.0 * pool["se_c_rel"]

    # 3. LIVE 16D ellipsoid tiled integral + Band-1 bsearch integral + P(occ)=1.0 SB overhead + owner RAM!
    V_lin = 2.0 ** pool["V_kept"]
    p2_mult = pool["p2_mult"]
    N_a_tbl = 4.0 * (2.0 ** pool["log2_Va"])
    N_7 = 2.0 ** pool["log2_N7_d12"]
    Lb = pool["Lb_avg"]
    w_tile, w_sb = coll["w_tile"], coll["w_sb"]
    s_vals = tuple(pool["s_vals"])
    s_max = max(s_vals)
    w_cls = (225.0 / 256.0, 225.0 / 2048.0, 23.0 / 2048.0)
    Rb_unit2 = pool["A2_prod"] / pool["r_ab"]
    L_max_live = 2.0 * math.sqrt(s_max * Rb_unit2 / 8.0)
    assert abs(L_max_live - pool["L_max"]) < 1e-4

    def phi_a_s(q):
        return sum(wc * (min(0.65 * sc, max(0.0, sc - q)) ** 4) for wc, sc in zip(w_cls, s_vals))

    nx, nm = 400, 48
    h_tot = c_tot = nt_row = sb_row = l_avg_live = 0.0
    for ix in range(nx):
        x = (ix + 0.5) / nx
        w7 = 3.5 * (x ** 2.5) / nx
        L = L_max_live * math.sqrt(1.0 - x)
        l_avg_live += w7 * L
        du = L / nm
        for im in range(nm):
            u = ((im + 0.5) / nm - 0.5) * L
            h_tot += w7 * du * phi_a_s(s_max * (x + (2.0 * u / L_max_live) ** 2))
        K_t = max(1, math.ceil(L / w_tile))
        nt_row += w7 * K_t
        u_lo = -0.5 * L
        for _ in range(K_t):
            u_hi = min(0.5 * L, u_lo + w_tile)
            Lj_int = max(1.0, float(math.ceil(u_hi - u_lo)))
            u_min = 0.0 if (u_lo <= 0.0 <= u_hi) else (u_lo if abs(u_lo) < abs(u_hi) else u_hi)
            q_min = s_max * (x + (2.0 * u_min / L_max_live) ** 2)
            sp_sb = 1.0 + (Lj_int - 1.0) / w_sb
            sb_row += w7 * sp_sb
            c_tot += w7 * (sp_sb * w_sb) * phi_a_s(q_min)
            u_lo = u_hi
    cand_ratio_live = c_tot / h_tot
    assert abs(l_avg_live - Lb) / Lb < 1e-3
    assert cand_ratio_live * 1.01 <= coll["cand_ratio_charged"]
    assert abs(nt_row - coll["nt_per_row"]) < 1e-2 and abs(sb_row - coll["sb_per_row"]) < 1e-2

    # Band-1 integral over [24, p_lo)
    b1_bs = b1_n1 = 0.0
    dl = (p_lo - 24.0) / 200.0
    for i in range(200):
        l0 = 24.0 + i * dl
        l1 = l0 + dl
        dn = 1.01 * (2.0**log2_li(l1) - 2.0**log2_li(l0))
        mu = N_a_tbl / (2.0 ** (0.5 * (l0 + l1)))
        b1_bs += dn * (coll["c_bsearch_hdr"] + coll["c_bsearch_step"] * (math.log2(mu + 1.0) + 1.0))
        b1_n1 += dn

    def price_b(p_hi):
        bld = 12.0 * N_a_tbl * math.log2(N_a_tbl)
        sp1 = 1.01 * (2.0**(2.0 * p_lo)) / (2.0 * p_lo * LN2 - 1.0) / 4.0
        bld += (b1_n1 / 4.0) * coll["c_list_a"] * N_a_tbl + coll["c_lbkt"] * sp1
        mem = (b1_n1 / 4.0) * 12.0 * N_a_tbl + 8.0 * sp1
        ipc1 = math.log(p_lo / 24.0) * coll["c_hit_lo"] * p2_mult + (b1_bs * N_7 * Lb + b1_n1 * N_7 * 36.0) / V_lin
        n2 = 1.01 * (2.0**log2_li(p_hi) - 2.0**log2_li(p_lo))
        sp2 = 1.01 * (2.0**(2.0 * p_hi) - 2.0**(2.0 * p_lo)) / (2.0 * p_hi * LN2 - 1.0) / 4.0
        bld += (n2 / 4.0) * coll["c_list_a"] * N_a_tbl + coll["c_lbkt"] * (sp2 / w_sb)
        mem += (n2 / 4.0) * 16.0 * N_a_tbl + 8.125 * (sp2 / w_sb)
        ipc_pr = (n2 * N_7 * (nt_row * coll["c_tile_hdr"] + sb_row * coll["c_sb_hdr"])) / V_lin
        ipc_hi = math.log(p_hi / p_lo) * (coll["c_cand"] * coll["cand_ratio_charged"] + coll["c_hit_hi"]) * p2_mult
        return bld, mem, ipc1 + ipc_pr + ipc_hi

    bg, mg, ipc_g = price_b(p_g)
    bf, mf, ipc_f = price_b(p_f)
    ipc_tot = ipc_g + ipc_f + coll["c_scan_zero"] + c_cof_live
    coll_u = pool["V_kept"] + math.log2(ipc_tot * coll["h2_allow"] / UNIT)
    build_u = math.log2((bg + bf) / UNIT)
    rows_d = dict(led["rows"])
    assert coll_u <= rows_d["relation_collection"]
    assert build_u <= rows_d["replay_a_tables"]
    assert 4.0 * N_a_tbl < 2.0**64, "Column owner S_g+S_f must fit < 2^64 B!"
    assert 16.0 * N_a_tbl + ((2.0**p_g) / w_sb) * 8.125 + w_tile * (1 << 20) < 2.0**64, "Max Band-2 table must fit < 2^64 B!"

    # 4. LIVE band_split & first-principles ecm_split + Chebyshev tail < 2^-70!
    np_band = 1.25506 * (2.0**48) / (48.0 * LN2)
    c_rel_band = np_band * (6 * (10 + 12) + 44 + 4) + (13 + 16) * 1024 + (1 << 26)
    n_rels_upper = 2.0 ** (pool["log2_rels"] + 0.55)
    band_u_live = math.log2(n_rels_upper * c_rel_band / UNIT)
    assert band_u_live <= rows_d["band_split"]
    ef_1 = ecm_instr(lf, 640, 0.999)
    eg_1 = ecm_instr(lg, 768, 0.999)
    c_ecm_rel = (13.0 * 1.002 * ef_1 + 14.0 * 40.0 * 2.0 * 640.0 * MONT5 + 13.0 * (1 << 26)
               + 16.0 * 1.002 * eg_1 + 17.0 * 40.0 * 2.0 * 768.0 * MONT5 + 16.0 * (1 << 26))
    ew_ecm_live = (c_ecm_rel / UNIT) * n_rels_upper
    ecm_sched_u = max(ef_1, eg_1) / UNIT
    ecm_var_live = 0.00201 * ecm_sched_u * ew_ecm_live
    ecm_gap = (2.0 ** rows_d["ecm_split_expected"]) - ew_ecm_live
    assert ecm_gap > 0
    ecm_cheb_p = ecm_var_live / (ecm_gap * ecm_gap)
    assert ecm_cheb_p < 2.0**-70

    # 5. Recompute 2^16..2^23 filtsim_z12 monotonicity, exact k5 from exp_unit, 1-owner LA RAM (x+y+tab)!
    lad = la["ladder"]
    assert [x["log2U"] for x in lad] == list(range(16, 24))
    for i in range(1, len(lad)):
        assert lad[i]["nl"] < lad[i - 1]["nl"] and lad[i]["sw"] < lad[i - 1]["sw"]
    r23 = lad[-1]
    cap_add = (0.045 / 2.0) * (2.0**(-log2_Up / 4.0)) / math.sqrt(0.0001)
    assert (r23["nl"] + 2.0 * r23["se_l"]) * 1.0159 + cap_add <= la["nl_cap"]
    assert (r23["sw"] + 2.0 * r23["se_sw"]) * 1.0159 + cap_add <= la["sw_cap"]
    assert r23["unit_frac"] >= 0.974 and r23["exp_unit"] >= 0.930
    eu_lo = la["eu_23"] - 2.0 * la["se_eu_23"]
    k5_live = eu_lo * 7.80 + (1.0 - eu_lo) * 33.0
    assert k5_live <= la["k5_charged"]
    tab_e = {2: 28.0, 3: 18.25, 4: 14.0}
    n_tb = {g: int(2.0 * (2.0**log2_li(x)) / 12.0) for g, x in tab_e.items()}
    te = sum((12**g) // 2 * math.comb(n_tb[g] + g - 1, g) for g in tab_e) + 12 * n_tb[2]
    N = 2.0**(log2_Up + math.log2(la["nl_cap"]))
    N_sw = 2.0**(log2_Up + math.log2(la["sw_cap"]))
    m, n, n_sm = la["m"], la["n"], la["n_sm"]
    la_mem_full = 2.0 * N_sw * 160.0 + (n / 5.0) * 160.0 * te
    assert la_mem_full < 2.0**64, "1-owner LA RAM (input x + output y + dense table) must fit < 2^64 B!"
    c_prod = (la["sw_cap"] / la["nl_cap"]) * (201.0 * la["k5_charged"] * la["tab_ratio_charged"] + la["c_row_chi"]) + 18.6 * te / N
    bw_sig = 2.0**-7
    kap = ((1.0 + bw_sig) * (N + n_sm) + 6.0 / (1.0 / m + 1.0 / n)) / N
    per_n2 = (c_prod + 60.0) * (n / m + 1.0) * kap + (1.0 / m + 1.0 / n) * (c_prod + n * 97.0 + 270.0)
    la_u = math.log2((N * N * per_n2 + N * (m + n) * 270.0) / UNIT)
    L_bm = (1.0 + bw_sig) * (N + n_sm) * (1.0 / m + 1.0 / n) + 6.0
    fft_L = 8.5 * (18.0 * L_bm) * (math.log2(18.0 * L_bm)**2) + 1200.0 + 504.0 * (18.0 * L_bm)
    gen_u = math.log2((4.0 * ((m + n)**3) * math.log2(L_bm) * (fft_L + 640.0 * L_bm)) / UNIT)
    assert la_u <= rows_d["linear_algebra"]
    assert gen_u <= rows_d["matrix_generator"]

    # 6. LIVE 78-node bottom-up descent recurrence + boot derivation + [11]/[7e] checks!
    dp = desc["params"]
    plan = desc["plan"]
    assert len(plan) == 2 * desc["n_bins"] == 78
    W_d, PHI_d, CN_d = {}, {}, {}
    bins_seen = []
    for row in plan:
        lo, hi, side, how, n_c, nA, kA, kB, R_n, mu, mc, sv_l2, sc_l2, own_l2_ref, W_l2_ref, phi_ref, en_l2_ref = row
        bn = (lo, hi)
        if not bins_seen or bins_seen[-1] != bn:
            bins_seen.append(bn)
        wmax = lambda s, lim: max([W_d[b, s] for b in bins_seen if b[0] < lim - 1e-9 and (b, s) in W_d] + [0.0])
        pmax = lambda s, lim: max([PHI_d[b, s] for b in bins_seen if b[0] < lim - 1e-9 and (b, s) in PHI_d] + [0.0])
        cmax = lambda s, lim: max([CN_d[b, s] for b in bins_seen if b[0] < lim - 1e-9 and (b, s) in CN_d] + [0.0])
        low_k = [k for k in W_d if k[0][0] <= dp["thin_max"]]
        pm = max([PHI_d[k] for k in low_k] + [0.0])
        wm = max([W_d[k] for k in low_k] + [0.0])
        cm = max([CN_d[k] for k in low_k] + [0.0])

        psi = kA * (pmax("f", nA) + pmax("g", nA)) + kB * (pmax("f", n_c) + pmax("g", n_c)) + mc * pm
        sub = kA * (wmax("f", nA) + wmax("g", nA)) + kB * (wmax("f", n_c) + wmax("g", n_c)) + mc * wm
        sub_c = kA * (cmax("f", nA) + cmax("g", nA)) + kB * (cmax("f", n_c) + cmax("g", n_c)) + mc * cm
        win = (kA if nA > dp["lx"] + 1e-9 else 0) + (kB if n_c > nA + 1e-9 else 0)
        r_nb = dp["nb_r1"] if win <= 1 else dp["nb_r2"]
        if how == "sieve":
            frac = lazy_fraction(mu, r_nb, psi)
            own = dp["lll_node"] + (2.0**sv_l2) + (2.0**sc_l2) * frac
        else:
            own = 2.0**own_l2_ref
        W_d[bn, side] = own + sub / (1.0 - psi)
        PHI_d[bn, side] = max((1.0 + (mu * (1.0 - psi)) / r_nb) ** (-r_nb) + dp["disc"], dp["phi_floor"])
        CN_d[bn, side] = 1.0 + sub_c / (1.0 - psi)
        assert abs(math.log2(W_d[bn, side]) - W_l2_ref) < 2e-3
        assert abs(PHI_d[bn, side] - phi_ref) < 1e-4
        assert abs(math.log2(CN_d[bn, side]) - en_l2_ref) < 2e-3

    bt = desc["boot"]
    assert bt["m_max"] >= desc["audit_11"]["ball_max"], "Boot m_max must strictly bound live [11] ball_max!"
    low_k = [k for k in W_d if k[0][0] <= dp["thin_max"]]
    pm = max(PHI_d[k] for k in low_k)
    wm = max(W_d[k] for k in low_k)
    cm = max(CN_d[k] for k in low_k)
    pmax_g = lambda lim: max(PHI_d[b, "g"] for b in bins_seen if b[0] < lim - 1e-9)
    wmax_g = lambda lim: max(W_d[b, "g"] for b in bins_seen if b[0] < lim - 1e-9)
    cmax_g = lambda lim: max(CN_d[b, "g"] for b in bins_seen if b[0] < lim - 1e-9)
    psi_b = bt["kA"] * pmax_g(bt["nA"]) + pmax_g(float(bt["n1"])) + bt["miss"] * pm
    phi_boot_live = psi_b ** dp["boot_rels"]
    det_b = 1.0 - (1.0 - dp["ecm_success"]) * (bt["m_max"] / 20.0 + bt["kA"] + 1.0)
    per_rel_b = 1.0 / ((2.0**bt["pr_log2"]) * det_b * dp["udup_boot"] * dp["ret_boot"])
    trial_b = dp["norms_boot"] + ecm_instr(bt["nA"], bt["m_max"], dp["ecm_success"]) + dp["prime_boot"] + dp["ideal_boot"]
    balls_b = (1.0 - phi_boot_live) / (1.0 - psi_b) * (per_rel_b / bt["cells_b"] + dp["ball_pair"])
    sub_b = bt["kA"] * wmax_g(bt["nA"]) + wmax_g(float(bt["n1"])) + bt["miss"] * wm
    one_log_live = dp["lll_node"] * balls_b + (balls_b * bt["cells_b"]) * trial_b + sub_b / (1.0 - psi_b)
    two_log_live = math.log2(2.0 * one_log_live / UNIT)
    EN_two_live = 2.0 * (bt["kA"] * cmax_g(bt["nA"]) + cmax_g(float(bt["n1"])) + bt["miss"] * cm) / (1.0 - psi_b)
    log2_EN_live = math.log2(EN_two_live)
    assert abs(two_log_live - desc["two_log_meas"]) < 2e-3
    assert abs(phi_boot_live - desc["boot_phi_meas"]) / desc["boot_phi_meas"] < 1e-3
    assert abs(log2_EN_live - desc["log2_EN_meas"]) < 2e-3
    assert two_log_live < desc["two_log_upper"]
    assert phi_boot_live <= desc["boot_phi_upper"]
    assert log2_EN_live - desc["node_cap_log2"] <= desc["counter_abort_log2_upper"]
    assert desc["audit_11"]["kappa_hat"] < 10.0 and desc["audit_11"]["ball_pair_charged"] == 2**24
    assert desc["audit_7e"]["delta_cf"] <= 10.0 and desc["audit_7e"]["delta_cg"] <= 4.0

    # 7. Exact 80-digit Decimal + Fraction ledger, STRICT ADDITIVE UNION BOUND, gate, and submission.json!
    with localcontext() as ctx:
        ctx.prec = 80
        two = Decimal(2); ln2 = two.ln()
        p2 = lambda v: (Decimal(str(v)) * ln2).exp()
        lg2 = lambda v: v.ln() / ln2

        core_cap = sum(p2(v) for _, v in led["rows"])
        core_cap_l2 = lg2(core_cap)
        assert core_cap_l2 <= Decimal(led["core_cap_claim"])
        exp_core = core_cap - p2(desc["descent_cap"]) + p2(desc["two_log_upper"])
        exp_core_l2 = lg2(exp_core)
        assert exp_core_l2 <= Decimal(led["exp_run_claim"])

        p_cap = p2(Decimal(str(desc["two_log_upper"])) - Decimal(str(desc["descent_cap"])))
        p_cnt = p2(Decimal(str(desc["counter_abort_log2_upper"])))
        p_rep = p2(Decimal(str(desc["repeat_log2_upper"])))
        order_d = Decimal(str(R_BLS))
        fail_core = (
            two**-66 + Decimal(str(led["cap_abort"])) + two**-70
            + Decimal(2) * Decimal(str(desc["boot_phi_upper"]))
            + p_cap + p_cnt + p_rep + Decimal(2) / order_d
        )
        fail_hr = (
            two**-66 + Decimal(str(led["cap_abort"])) + two**-70
            + Decimal(2) * Decimal(str(desc["boot_phi_upper"]))
            + Decimal("1.1") * p_cap + p_cnt + p_rep + Decimal(2) / order_d + two**-48
        )
        succ_core = Decimal(1) - fail_core
        succ_hr = Decimal(1) - fail_hr
        lo = F(led["p_floor"])
        assert F(str(succ_hr)) >= lo > F(39, 100)
        K = int(led["K_gate"])
        assert 0 < K < (1 << 64)
        assert K * lo >= F(39, 100) * (1 << 64) > (K - 1) * lo

        ew_p = lg2(Decimal(1) + Decimal(K) * (p2(led["exp_run_claim"]) + p2(led["wrapper_log2"])) / (two**64))
        ew_g1 = lg2(Decimal(1) + Decimal(K) * Decimal(10) * (p2(led["exp_run_claim"]) + p2(led["wrapper_log2"])) / (two**64))
        assert ew_p < Decimal(led["claim_p"])
        assert ew_g1 < Decimal(led["claim_g1"])

        # Check submission.json consistency (supports both C=20,000 and C=2,000 tracks!)
        manifest_ok = True
        track_claim = None
        for mp in (ROOT / "submission.json", HERE / "submission.json"):
            if mp.is_file():
                mf = json.loads(mp.read_text())
                sc = mf["claims"][0]["score"]
                declared = F(int(sc["numerator"]), int(sc["denominator"]))
                if declared == F(led["claim_p"]):
                    track_claim = led["claim_p"]
                    assert ("2^" + led["claim_p"]) in mf["title"]
                    assert ("2^" + led["exp_run_claim"]) in mf["claims"][0]["statement"]
                elif declared == F(led["claim_g1"]):
                    track_claim = led["claim_g1"]
                    assert ("2^" + led["claim_g1"]) in mf["title"]
                    assert ("2^" + led["exp_run_claim"]) in mf["claims"][0]["statement"]
                else:
                    manifest_ok = False
                break
        assert manifest_ok

    return dict(
        r_tot_live=round(r_tot, 6), c_cof_live=round(c_cof_live, 4),
        log2_rels_hist=round(log2_rels_h, 6), margin_net_bits=round(margin_net, 5),
        L_max_live=round(L_max_live, 2), cand_ratio_live=round(cand_ratio_live, 4),
        ipc_tot=round(ipc_tot, 4), coll_u=round(coll_u, 6), build_u=round(build_u, 4),
        ecm_cheb_log2=round(math.log2(ecm_cheb_p), 2),
        k5_live=round(k5_live, 4), la_u=round(la_u, 6), gen_u=round(gen_u, 4),
        la_mem_full_log2=round(math.log2(la_mem_full), 3),
        desc_nodes_live=len(plan), desc_two_log_live=round(two_log_live, 6),
        desc_phi_boot_live=float(f"{phi_boot_live:.4e}"), desc_EN_log2_live=round(log2_EN_live, 4),
        core_cap_log2=str(core_cap_l2)[:18], exp_run_log2=str(exp_core_l2)[:18],
        succ_headroom=str(succ_hr)[:10], p_floor=led["p_floor"], K_gate=K,
        ew_p_log2=str(ew_p)[:18], ew_g1_log2=str(ew_g1)[:18], track_claim=track_claim,
    )


# =============================================================================
# 6. Inline Fast H-X Smoke Check (Z[zeta_12] mod 1000033, 4 SM cols, 1-block pad neg control)
# =============================================================================
def check_hx_smoke():
    Q = 1000033  # == 1 mod 12
    chi = 385966
    assert (pow(chi, 4, Q) - pow(chi, 2, Q) + 1) % Q == 0
    units12 = [1, chi, pow(chi, 2, Q), pow(chi, 3, Q), (pow(chi, 2, Q) - 1) % Q, (pow(chi, 3, Q) - chi) % Q]
    units12 += [(-u) % Q for u in units12]
    rng = random.Random(12345)
    D, s_sm, m, n = 48, 4, 12, 8
    N = D - s_sm
    # Build sparse Z[zeta_12] live matrix M (D x N) + 4 dense SM cols S (D x 4)
    M_rows = []
    for i in range(D):
        cols = rng.sample(range(N), 14)
        M_rows.append([(c, rng.choice(units12)) for c in cols])
    S_cols = [[rng.randrange(1, Q) for _ in range(D)] for _ in range(s_sm)]

    def matvec_B(v):
        return [sum(val * v[c] for c, val in row) % Q for row in M_rows]

    # Build Krylov controllability C_d = [B z, B Y, ..., B^d Y] and observability O_e A
    Y = [list(col) for col in S_cols] + [matvec_B([rng.randrange(1, Q) for _ in range(D)]) for _ in range(n - s_sm)]
    d = math.ceil(D / n) + 1
    e = math.ceil(D / m) + 1

    def rank_mod(vecs):
        # vecs = list of column vectors (length D) -> rank over F_Q
        basis = {}
        for v in vecs:
            cur = [x % Q for x in v]
            for i in range(D):
                if cur[i]:
                    if i in basis:
                        b = basis[i]
                        f = cur[i]
                        for j in range(i, D):
                            cur[j] = (cur[j] - f * b[j]) % Q
                    else:
                        inv = pow(cur[i], -1, Q)
                        basis[i] = [(x * inv) % Q for x in cur]
                        break
        return len(basis)

    # Controllability C_d (columns in Im(B) of dim N)
    cols_C = []
    cur_Y = [matvec_B(y) for y in Y]
    for _ in range(d):
        cols_C.extend(cur_Y)
        cur_Y = [matvec_B(y) for y in cur_Y]
    rk_Cd = rank_mod(cols_C)

    # Observability O_e A for (a) block-diagonal X with 4 pad coords in 4 distinct blocks, (b) 1-block pad neg control!
    def obs_rank(one_block_pad):
        # Pad coordinates are N..N+3 (the 4 zero columns of B). X is m x D with m disjoint blocks!
        blocks = [[] for _ in range(m)]
        perm = list(range(N)); rng.shuffle(perm)
        for idx, c in enumerate(perm):
            blocks[idx % m].append((c, rng.randrange(1, Q)))
        for k in range(s_sm):
            blk = 0 if one_block_pad else k
            blocks[blk].append((N + k, rng.randrange(1, Q)))
        # Compute O_e * A where A = [M | S] (D x D)
        # Rows of O_e are X_l, X_l B, ..., X_l B^(e-1): multiply on right by A!
        # Easier: apply X^T B^t to every column of A (D columns) -> (m*e) x D matrix!
        A_cols = []
        for c in range(N):
            col = [0] * D
            for r in range(D):
                for cc, vv in M_rows[r]:
                    if cc == c: col[r] = vv
            A_cols.append(col)
        A_cols.extend(S_cols)
        obs_vecs = []
        for col in A_cols:
            out = []
            cur = col
            for _ in range(e):
                for blk in blocks:
                    out.append(sum(w * cur[idx] for idx, w in blk) % Q)
                cur = matvec_B(cur)
            obs_vecs.append(out)
        # Rank of obs_vecs (D vectors of length m*e)
        L = m * e; basis = {}
        for v in obs_vecs:
            cur = list(v)
            for i in range(L):
                if cur[i]:
                    if i in basis:
                        b = basis[i]; f = cur[i]
                        for j in range(i, L): cur[j] = (cur[j] - f * b[j]) % Q
                    else:
                        inv = pow(cur[i], -1, Q)
                        basis[i] = [(x * inv) % Q for x in cur]
                        break
        return len(basis)

    rk_obs = obs_rank(False)
    rk_neg = obs_rank(True)
    assert rk_Cd == N and rk_obs == D and rk_neg == D - (s_sm - 1)
    return dict(D=D, N=N, s_sm=s_sm, rank_Cd=rk_Cd, rank_OeA=rk_obs, neg_1blk_rank=rk_neg)


def main():
    raw = sys.stdin.read().strip()
    req = json.loads(raw) if raw else {}
    if req != {}:
        raise ValueError("audit_cm.py takes empty JSON object {}")
    data = json.loads((ROOT / "cm_data.json").read_text())
    toy = json.loads((ROOT / "toy_e2e.json").read_text())
    out = {}
    ck = {}

    out["algebra_bls12"] = check_algebra_bls12(); ck["algebra_bls12"] = True
    out["toy_e2e_dlp"] = check_toy_e2e_dlp(toy); ck["toy_e2e_dlp"] = True
    out["sigmaB_exactness"] = check_sigmaB_exactness(); ck["sigmaB_exactness"] = True
    out["anchor_norms"] = check_anchor_norms(data["pool"]["anchors"]); ck["anchor_norms"] = True
    out["srt_toy"] = srt_toy_check(); ck["srt_toy"] = True
    out["ledger_gate_claim"] = check_pool_coll_la_descent_ledger(data); ck["ledger_gate_claim"] = True
    out["hx_smoke"] = check_hx_smoke(); ck["hx_smoke"] = True

    res = dict(checks_ok=all(ck.values()), checks=ck, summary=out)
    print(json.dumps(res, separators=(",", ":")))
    if not res["checks_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
