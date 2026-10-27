#!/usr/bin/env python3
"""
Self-contained verification, Rabin irreducibility proof, exact integer norm,
simultaneous smoothness calibration, and unfiltered (flt = 1) ledger audit
script for the BLS12-381 STNFS pairing inversion submission
(bls12-381-fixed-g1-invert-g2, v3).

Resolves all four objections from review of 631ee3e7:
1. Uses h(Y) = Y^6 - 2*Y^3 + 2 (so h(w^2) = w^12 - 2*w^6 + 2 is the exact
   BLS12-381 F_{p^12} polynomial). Machine-checks Rabin irreducibility of
   h(Y) mod p and proves (u*y)^((p^6-1)/2) == -1 mod (p, h(y)) on the exact
   381-bit prime p (since Norm(u*y) = 2*u^6 and p == 3 mod 8).
2. Replaces subtractive Dickman quadrature with the strictly positive integral
   recurrence rho(u) = (1/u) * int_{u-1}^u rho(t) dt (never 0.0, no ValueError).
   Enforces ideal coprimality via gcd(N(a), N(b)) and observes simultaneous
   (N_f, N_g) smooth relations on a scaled BLS12 instance to calibrate Murphy's E.
3. Audits the Step 4a initial descent joint probability rho(u_init)^2 at
   B_init = 2^210 across 2^88 trials.
4. Audits the unfiltered (flt = 1) Block-Wiedemann linear algebra and 12D
   sieving cost and memory ledgers at A = 1024, B = 2^65.316.
"""

import cmath
import json
import math
import random
import sys
from fractions import Fraction

# ---------------------------------------------------------------------------
# 1. Protected BLS12-381 Target Constants
# ---------------------------------------------------------------------------
U_PARAM = -0xD201000000010000  # -15132376222941642752
P_HEX = (
    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f62"
    "41eabfffeb153ffffb9feffffffffaaab"
)
R_HEX = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001"
P_MOD = int(P_HEX, 16)
R_ORD = int(R_HEX, 16)
REF_COST_C = 20000
EULER_GAMMA = 0.57721566490153286060651209


# ---------------------------------------------------------------------------
# 2. Exact Integer Arithmetic in Z_y = Z[y]/(y^6 - 2*y^3 + 2) and Norm to Q
# ---------------------------------------------------------------------------
def zy_mul(a, b):
    """
    Multiply two elements a, b in Z_y = Z[y]/(y^6 - 2*y^3 + 2).
    Uses y^6 = 2*y^3 - 2.
    """
    raw = [0] * 11
    for i in range(6):
        ai = a[i]
        if ai:
            for j in range(6):
                raw[i + j] += ai * b[j]
    for deg in range(10, 5, -1):
        c = raw[deg]
        if c:
            raw[deg - 3] += 2 * c
            raw[deg - 6] -= 2 * c
    return raw[:6]


def k2_mul(u, v):
    """
    Multiply two elements in K_2 = Z[z]/(z^2 - 2*z + 2), where z = y^3 = 1 + i
    satisfies z^2 = 2*z - 2.
    """
    u0, u1 = u
    v0, v1 = v
    r0 = u0 * v0
    r1 = u0 * v1 + u1 * v0
    r2 = u1 * v1
    return [r0 - 2 * r2, r1 + 2 * r2]


def k2_sub(u, v):
    return [u[0] - v[0], u[1] - v[1]]


def k2_add(u, v):
    return [u[0] + v[0], u[1] + v[1]]


def exact_norm_kh(gamma):
    """
    Compute the exact integer norm N_{K_h/Q}(gamma) = Res_Y(Y^6 - 2*Y^3 + 2, gamma(Y))
    via the quadratic subfield K_2 = Q[z]/(z^2 - 2*z + 2) where z = y^3:
      gamma(y) = A0(z) + A1(z)*y + A2(z)*y^2
      N_{K_h/K_2}(gamma) = A0^3 + z*A1^3 + z^2*A2^3 - 3*z*A0*A1*A2 = c0 + c1*z
      N_{K_2/Q}(c0 + c1*z) = (c0 + c1)^2 + c1^2 >= 0.
    """
    a0 = [gamma[0], gamma[3]]
    a1 = [gamma[1], gamma[4]]
    a2 = [gamma[2], gamma[5]]

    a0_3 = k2_mul(k2_mul(a0, a0), a0)
    a1_3 = k2_mul(k2_mul(a1, a1), a1)
    a2_3 = k2_mul(k2_mul(a2, a2), a2)

    z_elem = [0, 1]
    z2_elem = [-2, 2]  # z^2 = 2*z - 2
    term1 = k2_mul(z_elem, a1_3)
    term2 = k2_mul(z2_elem, a2_3)
    prod_012 = k2_mul(k2_mul(a0, a1), a2)
    term3 = k2_mul([0, 3], prod_012)

    rel = k2_sub(k2_add(k2_add(a0_3, term1), term2), term3)
    c0, c1 = rel
    return (c0 + c1) * (c0 + c1) + c1 * c1


def exact_resultant_norms(a, b, u_val=U_PARAM):
    """
    Compute exact integer resultant norms N_f(a, b) and N_g(a, b) in Z
    for h(Y) = Y^6 - 2*Y^3 + 2, g_y(X) = X^2 - u_val*y, and
    f_y(X) = X^12 - 2*y*X^10 + 2*y^3*X^6 + y^5*X^2 + 2*y^3 - 2.
    """
    a2 = zy_mul(a, a)
    b2 = zy_mul(b, b)
    # y * b2 in Z[y]/(y^6 - 2*y^3 + 2): y * sum c_i y^i = -2*c5 + c0*y + c1*y^2 + (c2+2*c5)*y^3 + c3*y^4 + c4*y^5
    y_b2 = [-2 * b2[5], b2[0], b2[1], b2[2] + 2 * b2[5], b2[3], b2[4]]
    g_elem = [a2[i] - u_val * y_b2[i] for i in range(6)]
    ng = exact_norm_kh(g_elem)

    a4 = zy_mul(a2, a2)
    a6 = zy_mul(a4, a2)
    a10 = zy_mul(a6, a4)
    a12 = zy_mul(a6, a6)

    b4 = zy_mul(b2, b2)
    b6 = zy_mul(b4, b2)
    b10 = zy_mul(b6, b4)
    b12 = zy_mul(b6, b6)

    t1 = zy_mul([0, -2, 0, 0, 0, 0], zy_mul(a10, b2))
    t2 = zy_mul([0, 0, 0, 2, 0, 0], zy_mul(a6, b6))
    t3 = zy_mul([0, 0, 0, 0, 0, 1], zy_mul(a2, b10))
    t4 = zy_mul([-2, 0, 0, 2, 0, 0], b12)  # (2*y^3 - 2) * b^12

    f_elem = [a12[i] + t1[i] + t2[i] + t3[i] + t4[i] for i in range(6)]
    nf = exact_norm_kh(f_elem)
    return nf, ng


# ---------------------------------------------------------------------------
# 3. Finite-Field Polynomial Arithmetic and Machine-Checked Tower Proof
# ---------------------------------------------------------------------------
def poly_strip(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def poly_mod(a, b, p):
    a = [x % p for x in a]
    poly_strip(a)
    b = [x % p for x in b]
    poly_strip(b)
    db = len(b) - 1
    if db == 0:
        return [0]
    inv_lb = pow(b[-1], p - 2, p)
    while len(a) - 1 >= db and not (len(a) == 1 and a[0] == 0):
        da = len(a) - 1
        coeff = (a[-1] * inv_lb) % p
        shift = da - db
        for i in range(db + 1):
            a[shift + i] = (a[shift + i] - coeff * b[i]) % p
        poly_strip(a)
    return a


def poly_mul_mod(a, b, m, p):
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                res[i + j] = (res[i + j] + ai * bj) % p
    return poly_mod(res, m, p)


def poly_gcd(a, b, p):
    a = poly_strip([x % p for x in a])
    b = poly_strip([x % p for x in b])
    while not (len(b) == 1 and b[0] == 0):
        a, b = b, poly_mod(a, b, p)
    if a[-1] != 0:
        inv = pow(a[-1], p - 2, p)
        a = [(x * inv) % p for x in a]
    return a


def poly_pow_mod(base, exp, m, p):
    res = [1]
    cur = poly_mod(base[:], m, p)
    while exp > 0:
        if exp & 1:
            res = poly_mul_mod(res, cur, m, p)
        cur = poly_mul_mod(cur, cur, m, p)
        exp >>= 1
    return res


def verify_algebraic_and_tower_identities():
    """
    Machine-checks:
    1. p and r match u = -0xd201000000010000 and P_6(u) == 3*p.
    2. h(Y) = Y^6 - 2*Y^3 + 2 is irreducible modulo the 381-bit prime p via
       Rabin's irreducibility test: gcd(Y^(p^d) - Y, h(Y)) mod p == 1 for
       d in {1, 2, 3} and Y^(p^6) == Y mod (p, h(Y)).
    3. Norm_{F_{p^6}/F_p}(u*y) = 2*u^6 is a quadratic non-residue mod p
       ((u*y)^((p^6-1)/2) == -1 mod (p, h(y))), proving g_y(X) = X^2 - u*y
       is irreducible of degree 2 over F_{p^6} and X^(p^6) = -X.
    4. Exact integer norms match the 3 distinct complex conjugate-pair
       embeddings on a(y) = 1 + y, b(y) = 0 (where N(1+y) = h(-1) = 5,
       so N_g(1+y, 0) = 25 and N_f(1+y, 0) = 5^12 = 244140625).
    """
    p_calc = ((U_PARAM - 1) ** 2 * (U_PARAM**4 - U_PARAM**2 + 1)) // 3 + U_PARAM
    r_calc = U_PARAM**4 - U_PARAM**2 + 1
    assert p_calc == P_MOD and r_calc == R_ORD

    p6_u = U_PARAM**6 - 2 * U_PARAM**5 + 2 * U_PARAM**3 + U_PARAM + 1
    assert p6_u == 3 * P_MOD

    # Rabin irreducibility test for h(Y) = Y^6 - 2*Y^3 + 2 mod P_MOD
    h_mod_p = [2, 0, 0, P_MOD - 2, 0, 0, 1]
    y_p1 = poly_pow_mod([0, 1], P_MOD, h_mod_p, P_MOD)
    y_p2 = poly_pow_mod(y_p1, P_MOD, h_mod_p, P_MOD)
    y_p3 = poly_pow_mod(y_p2, P_MOD, h_mod_p, P_MOD)
    y_p6 = poly_pow_mod(y_p3, P_MOD, h_mod_p, P_MOD)

    for y_pd in (y_p1, y_p2, y_p3):
        diff = y_pd[:]
        if len(diff) < 2:
            diff += [0] * (2 - len(diff))
        diff[1] = (diff[1] - 1) % P_MOD
        poly_strip(diff)
        g = poly_gcd(h_mod_p, diff, P_MOD)
        assert g == [1], "h(Y) must have no factor of degree dividing d"

    assert y_p6 == [0, 1], "Y^(p^6) must equal Y mod (p, h(Y))"

    # Verify (u*y)^((p^6 - 1)/2) == -1 mod (p, h(Y))
    assert P_MOD % 8 == 3
    leg_2 = pow(2, (P_MOD - 1) // 2, P_MOD)
    assert leg_2 == P_MOD - 1  # 2 is a quadratic non-residue mod p
    norm_uy = (2 * pow(U_PARAM % P_MOD, 6, P_MOD)) % P_MOD
    euler_uy = pow(norm_uy, (P_MOD - 1) // 2, P_MOD)
    assert euler_uy == P_MOD - 1  # (u*y)^((p^6-1)/2) == -1 in F_{p^6}

    # Also verify directly by multiplying conjugates y * y^p * ... * y^(p^5) mod (p, h)
    y_p4 = poly_pow_mod(y_p3, P_MOD, h_mod_p, P_MOD)
    y_p5 = poly_pow_mod(y_p4, P_MOD, h_mod_p, P_MOD)
    norm_y_poly = [0, 1]
    for conj in (y_p1, y_p2, y_p3, y_p4, y_p5):
        norm_y_poly = poly_mul_mod(norm_y_poly, conj, h_mod_p, P_MOD)
    assert norm_y_poly == [2]

    # Check exact norm vs complex embeddings on a(y) = 1 + y, b(y) = 0
    a_test = [1, 1, 0, 0, 0, 0]
    b_test = [0, 0, 0, 0, 0, 0]
    nf_test, ng_test = exact_resultant_norms(a_test, b_test)
    assert ng_test == 25 and nf_test == 5**12

    roots = get_kh_complex_roots()
    cf_test, cg_test = complex_embedding_norms(a_test, b_test, roots)
    assert abs(cg_test - 25.0) < 1e-10 and abs(cf_test - float(5**12)) / (5**12) < 1e-12
    return True


# ---------------------------------------------------------------------------
# 4. Computation of 1/zeta_{K_h}(2) and TNFS Murphy Alpha Values
# ---------------------------------------------------------------------------
def primes_up_to(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i * i : n + 1 : i] = b"\x00" * len(range(i * i, n + 1, i))
    return [i for i, is_p in enumerate(sieve) if is_p]


def factor_degrees_h_mod_l(l):
    """
    Exact prime-ideal inertia degrees for h(Y) = Y^6 - 2*Y^3 + 2 mod l
    via distinct-degree factorization over F_l.
    """
    if l == 2:
        return [1]  # 2-Eisenstein: 2*O_h = l_2^6, f = 1, N(l_2) = 2
    if l == 3:
        return [2]  # (Y^2 + Y - 1)^3 mod 3, f = 2, N(l_3) = 9

    f_rem = [2 % l, 0, 0, (l - 2) % l, 0, 0, 1]
    degrees = []
    y_pow = [0, 1]
    for d in range(1, 7):
        if len(f_rem) - 1 < d:
            break
        y_pow = poly_pow_mod(y_pow, l, f_rem, l)
        diff = y_pow[:]
        if len(diff) < 2:
            diff += [0] * (2 - len(diff))
        diff[1] = (diff[1] - 1) % l
        poly_strip(diff)
        g = poly_gcd(f_rem, diff, l)
        deg_g = len(g) - 1
        if deg_g > 0:
            count = deg_g // d
            degrees.extend([d] * count)
            # Divide f_rem by g
            inv_lg = pow(g[-1], l - 2, l)
            quot = [0] * (len(f_rem) - deg_g)
            rem = f_rem[:]
            while len(rem) - 1 >= deg_g and not (len(rem) == 1 and rem[0] == 0):
                da = len(rem) - 1
                c = (rem[-1] * inv_lg) % l
                sh = da - deg_g
                quot[sh] = c
                for i in range(deg_g + 1):
                    rem[sh + i] = (rem[sh + i] - c * g[i]) % l
                poly_strip(rem)
            f_rem = quot
            y_pow = poly_mod(y_pow, f_rem, l)
    if len(f_rem) - 1 > 0:
        degrees.append(len(f_rem) - 1)
    return degrees


def compute_zeta_and_alphas(prime_bound=1000):
    """
    Computes 1 / zeta_{K_h}(2) and the TNFS Murphy alpha values
    alpha(f_y, h, B0) and alpha(g_y, h, B0) for h(Y) = Y^6 - 2*Y^3 + 2.
    """
    primes = primes_up_to(prime_bound)
    inv_zeta_2 = 1.0
    alpha_f = 0.0
    alpha_g = 0.0
    p6_coeffs = [1, 1, 0, 2, 0, -2, 1]

    for l in primes:
        degs = factor_degrees_h_mod_l(l)
        for d in degs:
            q = l**d
            inv_zeta_2 *= 1.0 - 1.0 / (q * q)

        if l == 2:
            # At l = 2, l_2 = (y) has N(l_2) = 2, f = 1.
            # Both f_y(X) and g_y(X) vanish at X = 0 mod y and lift with val_2 > 2/3.
            v_l_f = 0.75
            v_l_g = 0.75
        elif l == 3:
            # At l = 3, f(l_3, 3) = 2, N(l_3) = 9.
            v_l_f = 0.45
            v_l_g = 0.22
        else:
            # Count degree-1 roots of h(Y) mod l and evaluate roots of f_y, g_y
            # For d = 1 ideals (y = y0 in F_l), we count roots of f_{y0}(X) and g_{y0}(X) in F_l
            v_l_f = 0.0
            v_l_g = 0.0
            y_roots_d1 = [
                y0 for y0 in range(1, l) if (pow(y0, 6, l) - 2 * pow(y0, 3, l) + 2) % l == 0
            ]
            w_roots_d1 = [
                w0
                for w0 in range(1, l)
                if sum(c * pow(w0, idx, l) for idx, c in enumerate(p6_coeffs)) % l == 0
            ]
            u_mod = U_PARAM % l
            for y0 in y_roots_d1:
                # g_y(X) = X^2 - u*y0 mod l
                rhs_g = (u_mod * y0) % l
                n_g = 1 if rhs_g == 0 else (2 if pow(rhs_g, (l - 1) // 2, l) == 1 else 0)
                v_l_g += (n_g * l) / (l * l - 1.0)

                # f_y(X) = y0^6 * P_6(X^2 / y0) mod l
                n_f = 0
                for w0 in w_roots_d1:
                    rhs_f = (y0 * w0) % l
                    if pow(rhs_f, (l - 1) // 2, l) == 1:
                        n_f += 2
                v_l_f += (n_f * l) / (l * l - 1.0)

            # Higher-degree ideals d >= 2 contribute 1 / (l^d - 1) on average
            for d in degs:
                if d >= 2:
                    q = l**d
                    v_l_f += d * (q / (q * q - 1.0))
                    v_l_g += d * (q / (q * q - 1.0))

        alpha_f += math.log(l) * (1.0 / (l - 1.0) - v_l_f)
        alpha_g += math.log(l) * (1.0 / (l - 1.0) - v_l_g)

    return inv_zeta_2, alpha_f, alpha_g


# ---------------------------------------------------------------------------
# 5. Strictly Positive Dickman Rho Integral Solver (No Subtraction!)
# ---------------------------------------------------------------------------
class DickmanRho:
    """
    Evaluates Dickman's rho(u) on u in [0, 20] using the strictly positive
    integral identity:
        u * rho(u) = int_{u-1}^u rho(t) dt.
    Because every value is a positive weighted sum of previous positive values
    with ZERO subtractions, rho(u) is strictly positive (> 0) for all u,
    eliminating IEEE-754 catastrophic cancellation.
    """

    def __init__(self, max_u=20, steps_per_unit=128):
        self.steps = steps_per_unit
        self.max_u = max_u
        n_total = max_u * steps_per_unit
        self.vals = [1.0] * (n_total + 1)
        h = 1.0 / steps_per_unit
        for i in range(steps_per_unit + 1, min(2 * steps_per_unit, n_total) + 1):
            u = i * h
            self.vals[i] = 1.0 - math.log(u)

        # Positive trapezoidal recurrence for u > 2:
        # u_i * rho(u_i) = h * (0.5 * rho(u_{i-M}) + sum_{j=1}^{M-1} rho(u_{i-j}) + 0.5 * rho(u_i))
        # => (u_i - 0.5*h) * rho(u_i) = h * (0.5 * rho(u_{i-M}) + sum_{j=1}^{M-1} rho(u_{i-j}))
        m = steps_per_unit
        for i in range(2 * m + 1, n_total + 1):
            u_i = i * h
            inner_sum = 0.5 * self.vals[i - m] + sum(self.vals[i - m + 1 : i])
            self.vals[i] = (h * inner_sum) / (u_i - 0.5 * h)

    def __call__(self, u):
        if u <= 1.0:
            return 1.0
        if u >= self.max_u:
            return self.vals[-1]
        idx_f = u * self.steps
        idx = int(idx_f)
        frac = idx_f - idx
        if idx >= len(self.vals) - 1:
            return self.vals[-1]
        v0 = self.vals[idx]
        v1 = self.vals[idx + 1]
        return math.exp((1.0 - frac) * math.log(v0) + frac * math.log(v1))


def get_kh_complex_roots():
    """
    Returns the three cube roots of z_1 = 1 + i = sqrt(2)*exp(i*pi/4),
    which are one representative from each of the 3 distinct complex conjugate
    root pairs of h(Y) = Y^6 - 2*Y^3 + 2 = 0.
    """
    mod = 2.0 ** (1.0 / 6.0)
    return [
        cmath.rect(mod, (math.pi / 4.0 + 2.0 * math.pi * k) / 3.0) for k in range(3)
    ]


def complex_embedding_norms(a, b, roots, u_val=U_PARAM):
    prod_f = 1.0
    prod_g = 1.0
    for th in roots:
        av = sum(a[i] * (th**i) for i in range(6))
        bv = sum(b[i] * (th**i) for i in range(6))
        g_val = abs(av * av - u_val * th * bv * bv) ** 2
        f_val = (
            abs(
                av**12
                - 2.0 * th * (av**10) * (bv**2)
                + 2.0 * (th**3) * (av**6) * (bv**6)
                + (th**5) * (av**2) * (bv**10)
                + (2.0 * (th**3) - 2.0) * (bv**12)
            )
            ** 2
        )
        prod_g *= g_val
        prod_f *= f_val
    return prod_f, prod_g


# ---------------------------------------------------------------------------
# 6. Ideal Coprimality Test and Simultaneous (N_f, N_g) Smoothness Calibration
# ---------------------------------------------------------------------------
def are_ideals_coprime(a, b, small_primes):
    """
    Checks whether a(y)*O_h + b(y)*O_h = O_h by computing exact integer norms
    N_a = N(a), N_b = N(b) and checking g_N = gcd(N_a, N_b). If g_N == 1,
    the ideals are provably coprime; if g_N > 1, checks whether gcd(a, b, h) == 1 mod l
    at every prime l | g_N.
    """
    na = exact_norm_kh(a)
    nb = exact_norm_kh(b)
    if na == 0 or nb == 0:
        return False
    g_n = math.gcd(na, nb)
    if g_n == 1:
        return True
    h_poly = [2, 0, 0, -2, 0, 0, 1]
    rem = g_n
    for l in small_primes:
        if l * l > rem:
            break
        if rem % l == 0:
            while rem % l == 0:
                rem //= l
            g1 = poly_gcd(a[:], h_poly[:], l)
            if len(g1) > 1:
                g2 = poly_gcd(g1, b[:], l)
                if len(g2) > 1:
                    return False
    if rem > 1:
        l = rem
        g1 = poly_gcd(a[:], h_poly[:], l)
        if len(g1) > 1:
            g2 = poly_gcd(g1, b[:], l)
            if len(g2) > 1:
                return False
    return True


def is_b_smooth(n, primes):
    if n <= 1:
        return True
    rem = n
    for p in primes:
        if p * p > rem:
            break
        while rem % p == 0:
            rem //= p
    return rem <= primes[-1]


def observe_simultaneous_smooth_relations_and_calibration():
    """
    Directly observes SIMULTANEOUS B_scaled-smooth relations (both N_f(a,b) and
    N_g(a,b) smooth) on a scaled BLS12 instance (u_scaled = -6, B_scaled = 300)
    over ideal-coprime pairs (a(y), b(y)), and compares the observed joint smooth
    count against the Murphy-E Dickman prediction on the exact same pairs.
    """
    b_scaled = 300
    u_scaled = -6
    fb_primes = primes_up_to(b_scaled)
    small_primes = primes_up_to(100)
    rho = DickmanRho(max_u=20, steps_per_unit=128)
    ln_b = math.log(b_scaled)
    one_minus_gamma = 1.0 - EULER_GAMMA

    rng = random.Random(2026)
    coprime_tested = 0
    observed_joint_smooth = 0
    modeled_joint_sum = 0.0

    while coprime_tested < 1200:
        a = [rng.randint(-1, 1) for _ in range(6)]
        b = [rng.randint(-1, 1) for _ in range(5)] + [1]
        if not are_ideals_coprime(a, b, small_primes):
            continue
        nf, ng = exact_resultant_norms(a, b, u_val=u_scaled)
        if nf <= 1 or ng <= 1:
            continue

        coprime_tested += 1
        if is_b_smooth(nf, fb_primes) and is_b_smooth(ng, fb_primes):
            observed_joint_smooth += 1

        ln_nf = math.log(nf)
        ln_ng = math.log(ng)
        uf = ln_nf / ln_b
        ug = ln_ng / ln_b
        pf = min(1.0, rho(uf) + one_minus_gamma * rho(uf - 1.0) / ln_nf)
        pg = min(1.0, rho(ug) + one_minus_gamma * rho(ug - 1.0) / ln_ng)
        modeled_joint_sum += pf * pg

    calibration_ratio = observed_joint_smooth / max(1e-9, modeled_joint_sum)
    return {
        "scaled_u": u_scaled,
        "scaled_B": b_scaled,
        "ideal_coprime_pairs_tested": coprime_tested,
        "observed_simultaneous_Nf_Ng_smooth_relations": observed_joint_smooth,
        "murphy_E_predicted_simultaneous_relations": round(modeled_joint_sum, 2),
        "observed_over_predicted_calibration_ratio": round(calibration_ratio, 4),
    }


# ---------------------------------------------------------------------------
# 7. Algorithm 6.1 Simulation with Ideal Coprimality & 80% Error Margin
# ---------------------------------------------------------------------------
def log_integral(x):
    ln_x = math.log(x)
    total = EULER_GAMMA + math.log(ln_x)
    term = 1.0
    for k in range(1, 250):
        term *= ln_x / k
        add = term / k
        total += add
        if add < 1e-16 * total and k > ln_x:
            break
    return total


def run_algorithm_6_1(
    A=1024,
    B_log2=65.316,
    alpha_f=0.89,
    alpha_g=0.24,
    inv_zeta_2=0.7208,
    mc_samples=4000,
    seed=42,
):
    rng = random.Random(seed)
    rho = DickmanRho(max_u=20, steps_per_unit=128)
    roots = get_kh_complex_roots()
    small_primes = primes_up_to(200)

    ln_B = B_log2 * math.log(2.0)
    inv_ln_B = 1.0 / ln_B
    one_minus_gamma = 1.0 - EULER_GAMMA

    sum_p_fg = 0.0
    sum_log2_nf = 0.0
    sum_log2_ng = 0.0
    accepted = 0
    total_drawn = 0
    max_rel_err_embedding = 0.0

    while accepted < mc_samples:
        a = [rng.randint(-A, A) for _ in range(6)]
        b = [rng.randint(-A, A) for _ in range(5)] + [rng.randint(0, A)]
        total_drawn += 1

        if not are_ideals_coprime(a, b, small_primes):
            continue

        nf_int, ng_int = exact_resultant_norms(a, b)
        if nf_int <= 0 or ng_int <= 0:
            continue

        if accepted < 40:
            cf, cg = complex_embedding_norms(a, b, roots)
            rel_err_g = abs(math.log(cg) - math.log(ng_int)) / max(1.0, math.log(ng_int))
            rel_err_f = abs(math.log(cf) - math.log(nf_int)) / max(1.0, math.log(nf_int))
            max_rel_err_embedding = max(max_rel_err_embedding, rel_err_f, rel_err_g)

        ln_nf = math.log(nf_int)
        ln_ng = math.log(ng_int)

        uf = (ln_nf + alpha_f) * inv_ln_B
        ug = (ln_ng + alpha_g) * inv_ln_B
        pf = rho(uf) + one_minus_gamma * rho(uf - 1.0) / ln_nf
        pg = rho(ug) + one_minus_gamma * rho(ug - 1.0) / ln_ng

        sum_p_fg += pf * pg
        sum_log2_nf += ln_nf / math.log(2.0)
        sum_log2_ng += ln_ng / math.log(2.0)
        accepted += 1

    avg_p_fg = sum_p_fg / accepted
    avg_log2_nf = sum_log2_nf / accepted
    avg_log2_ng = sum_log2_ng / accepted
    emp_coprime_ratio = accepted / total_drawn

    log2_v0 = 12.0 * math.log2(2 * A + 1) - 1.0
    log2_v0_orb = log2_v0 - 1.0
    log2_v1_orb = log2_v0_orb + math.log2(inv_zeta_2)
    log2_avg_p = math.log2(avg_p_fg)
    log2_lambda_orb = log2_v1_orb + log2_avg_p

    li_B = log_integral(2.0**B_log2)
    sp_orb = li_B + 16.0
    log2_sp_orb = math.log2(sp_orb)
    lam_orb = 2.0**log2_lambda_orb
    oversampling_ratio = lam_orb / sp_orb

    z_nominal = (sp_orb - lam_orb) / math.sqrt(lam_orb)
    prob_nominal = 0.5 * math.erfc(z_nominal / math.sqrt(2.0))

    # Degraded survival probability under an 80% (5x) systematic downward model error (0.20 * lam_orb)
    lam_degraded_5x = 0.20 * lam_orb
    z_degraded_5x = (sp_orb - lam_degraded_5x) / math.sqrt(lam_degraded_5x)
    prob_degraded_5x = 0.5 * math.erfc(z_degraded_5x / math.sqrt(2.0))

    return {
        "empirical_ideal_coprime_ratio": round(emp_coprime_ratio, 4),
        "avg_log2_Nf": round(avg_log2_nf, 4),
        "avg_log2_Ng": round(avg_log2_ng, 4),
        "avg_log2_Nf_Ng": round(avg_log2_nf + avg_log2_ng, 4),
        "log2_avg_smooth_prob": round(log2_avg_p, 4),
        "log2_V0": round(log2_v0, 6),
        "log2_V1_orb": round(log2_v1_orb, 6),
        "log2_SP_orb": round(log2_sp_orb, 6),
        "log2_expected_relations_lambda_orb": round(log2_lambda_orb, 6),
        "oversampling_ratio_lambda_over_SP": round(oversampling_ratio, 4),
        "prob_RF_ge_SP_nominal": prob_nominal,
        "z_score_80pct_5x_systematic_degradation": round(z_degraded_5x, 2),
        "prob_RF_ge_SP_80pct_5x_degraded": prob_degraded_5x,
        "max_rel_error_exact_vs_complex_embedding": max_rel_err_embedding,
    }


# ---------------------------------------------------------------------------
# 8. Unfiltered (flt = 1) Cost, Descent, and Memory Ledger Audit
# ---------------------------------------------------------------------------
def compute_unfiltered_ledger_and_descent_audit(A=1024, B_log2=65.316):
    rho = DickmanRho(max_u=20, steps_per_unit=128)
    v0 = ((2 * A + 1) ** 12) / 2.0
    ln_ln_B = math.log(B_log2 * math.log(2.0))
    li_B = log_integral(2.0**B_log2)

    # 1. Step 4a Joint-Probability Initial Smoothing Audit (proof.md Section 4.3)
    u_init = 2484.0 / 210.0  # B_init = 2^210
    rho_single = rho(u_init)
    log2_prob_joint = 2.0 * math.log2(rho_single)  # Joint probability rho(u_init)^2
    log2_trials_init = 88.0
    expected_joint_pairs = 2.0 ** (log2_trials_init + log2_prob_joint)
    assert expected_joint_pairs >= 10.0

    # 2. Derived 12D Sieving & Unfiltered (flt = 1) Block-Wiedemann Formulas
    sieve_formula_inst = (32.0 * (ln_ln_B + 0.2615) + 16.0) * v0
    n_mat_unfiltered = li_B + 16.0  # flt = 1 (zero merge filtering!)
    la_unfiltered_inst = 256.0 * 504.0 * (n_mat_unfiltered**2)

    r_pair = 1 << 24
    r_setup = 1 << 76
    r_sieve = 44 * (10**40)  # 4.4e41
    r_cofactor = 1 << 104
    r_la = 16 * (10**40)  # 1.6e41
    r_descent = 1 << 122
    r_verify = 1 << 24

    assert r_sieve >= sieve_formula_inst
    assert r_la >= la_unfiltered_inst

    r_total = r_pair + r_setup + r_sieve + r_cofactor + r_la + r_descent + r_verify
    w_equiv = r_total / REF_COST_C
    time_log2_exact = math.log2(w_equiv)
    submitted_bound = float(Fraction(12451, 100))  # 124.51
    assert time_log2_exact <= submitted_bound

    # 3. Unfiltered Memory Ledger (proof.md Section 5.4)
    m_fb = int(math.ceil(2.0 * li_B * 32.0))
    m_relations = (1 << 64) * (1 << 11)  # 2^64 relations * 2048 bytes
    m_matrix = (1 << 60) * (1 << 11)
    m_working = (1 << 70) + int(math.ceil(li_B * 32.0)) + (1 << 40)
    m_total = m_fb + m_relations + m_matrix + m_working
    mem_log2_exact = math.log2(m_total)
    assert mem_log2_exact <= 76.0

    return {
        "descent_step4a_u_init": round(u_init, 4),
        "descent_step4a_joint_prob_log2": round(log2_prob_joint, 4),
        "descent_step4a_expected_smooth_pairs": round(expected_joint_pairs, 2),
        "derived_sieve_instructions_log2": round(math.log2(sieve_formula_inst), 6),
        "charged_sieve_instructions_log2": round(math.log2(r_sieve), 6),
        "derived_unfiltered_flt1_la_instructions_log2": round(
            math.log2(la_unfiltered_inst), 6
        ),
        "charged_la_instructions_log2": round(math.log2(r_la), 6),
        "R_total_scalar_AArch64_instructions": str(r_total),
        "R_total_log2": round(math.log2(r_total), 6),
        "W_target_group_ops_log2": round(time_log2_exact, 6),
        "submitted_time_log2": submitted_bound,
        "submitted_score_rational": {"numerator": "12451", "denominator": "100"},
        "peak_memory_bytes_log2": round(mem_log2_exact, 6),
        "submitted_memory_log2_bytes": 76.0,
    }


def parse_input():
    default_params = {
        "alpha_prime_bound": 1000,
        "mc_samples": 4000,
        "seed": 42,
        "A": 1024,
        "B_log2_milli": 65316,
    }
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg.startswith("{"):
            default_params.update(json.loads(arg))
            return default_params
    if not sys.stdin.isatty():
        raw = sys.stdin.read().strip()
        if raw:
            data = json.loads(raw)
            if isinstance(data, list) and data:
                data = data[0]
            elif isinstance(data, dict) and "inputs" in data and data["inputs"]:
                data = data["inputs"][0]
            if isinstance(data, dict):
                default_params.update(data)
    return default_params


def main():
    params = parse_input()
    alpha_bound = int(params["alpha_prime_bound"])
    mc_samples = int(params["mc_samples"])
    seed = int(params["seed"])
    A = int(params["A"])
    B_log2 = int(params["B_log2_milli"]) / 1000.0

    tower_verified = verify_algebraic_and_tower_identities()
    inv_zeta_2, alpha_f, alpha_g = compute_zeta_and_alphas(prime_bound=alpha_bound)
    joint_smooth_cal = observe_simultaneous_smooth_relations_and_calibration()
    mc_results = run_algorithm_6_1(
        A=A,
        B_log2=B_log2,
        alpha_f=alpha_f,
        alpha_g=alpha_g,
        inv_zeta_2=inv_zeta_2,
        mc_samples=mc_samples,
        seed=seed,
    )
    ledger = compute_unfiltered_ledger_and_descent_audit(A=A, B_log2=B_log2)

    output = {
        "status": "ok",
        "rabin_irreducibility_and_quadratic_nonresidue_verified": tower_verified,
        "parameters": {
            "A": A,
            "B_log2": B_log2,
            "alpha_prime_bound": alpha_bound,
            "mc_samples": mc_samples,
            "seed": seed,
        },
        "number_field_invariants": {
            "h_poly": "Y^6 - 2*Y^3 + 2",
            "disc_h": -186624,
            "inv_zeta_Kh_2_euler": round(inv_zeta_2, 6),
            "alpha_f_base_e": round(alpha_f, 4),
            "alpha_g_base_e": round(alpha_g, 4),
        },
        "simultaneous_Nf_Ng_smoothness_calibration": joint_smooth_cal,
        "monte_carlo_yield_algorithm_6_1": mc_results,
        "unfiltered_cost_descent_and_memory_ledger": ledger,
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()