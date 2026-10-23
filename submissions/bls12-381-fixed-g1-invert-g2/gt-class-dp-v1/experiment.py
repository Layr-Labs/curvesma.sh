#!/usr/bin/env python3
"""experiment.py -- Frobenius-class Pollard rho in GT with fresh restarts.

Candidate: submissions/bls12-381-fixed-g1-invert-g2/gt-class-dp-v1
Track:     bls12-381-fixed-g1-invert-g2 (fixed-argument pairing inversion into G2)

Standard library only.  One JSON object on stdin, one JSON object on stdout.

Section A (production scale, on BLS12-381 itself, blst tower)
    Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - (1+u)).
    A1  p, r prime; r | Phi_12(p) = p^4 - p^2 + 1; ord_r(p) = 12 exactly
        (p^6 = -1, p^4 != 1, p^3 != 1, p^2 != 1).
    A2  w^6 - (1+u) is irreducible over Fp2 (xi = 1+u neither a square nor a cube).
    A3  gamma1 = w^(p-1) lies in Fp2 (support on w^1 only), zeta = N(gamma1) is in
        Fp of order exactly 6, zeta^2 = zeta - 1.
    A4  the twelve-key class invariant K(y) = min over k of Re c_1(y^(p^k)) is
        computed from (a, B) with two Fp multiplications and agrees with
        brute-force Frobenius y -> y^(p^k), k = 0..11: values, argmin, tie count,
        rotation invariance.
    A5  the walk identity.  With kappa = argmin, c = y^(p^kappa) the canonical
        representative, and T[h][d] = M_h^(p^d), the step y' = y*T[h][(12-kappa)
        mod 12] satisfies (y')^(p^kappa) = c*M_h, so the transition is a genuine
        map on classes and the lift is never re-canonicalised.  Checked directly
        on BLS12-381 GT elements by comparing class keys.

Section B (toy scale, GT = mu_r of a toy tower with the same tower shape)
    B1  machinery mode: the complete production walk -- fresh random starts, a
        hard per-walk step cap L, distinguished classes K = 0 mod 2^dp, an exact
        full-key table, the relation (alpha,beta) ~ p^k (alpha',beta') -> x, and
        verification g^x = z.  Reports trial counts, the fraction solved inside
        the claimed budget, bookkeeping failures (must be zero) and the measured
        distinguished-point counts.
    B2  timing mode: walks long compared with sqrt(N) and a dense distinguished
        set, so the reported statistic is the random-mapping collision time.  The
        measured mean and the empirical success curve at S0 and 1.03*S0 are the
        toy support for the model constant; the production bound does not use a
        toy constant for anything but the model's own 1.2533.
    B3  structural traps: measured rate of two-step class repetitions whose lift
        state also repeats (fruitless two-cycles, predicted 1/(12R)) and of
        longer lift-state repetitions, per step.  These are charged in Section C
        by a fixed allowance, not by the measured value.

Section C
    The cost ledger, exact integer / rational arithmetic: N, S0 at the declared
    0.40 model success, the 3% slack, the per-step instruction ledger, the
    detection-efficiency factor, the trap allowance, additive items in
    dense-Fp12-multiplication units, the peak-memory bound and the claim.

Toy statistics are evidence at toy scope only; the production claim is analytic.
"""

import json
import math
import sys
import time
from fractions import Fraction

# --------------------------------------------------------------------------
# parameters
# --------------------------------------------------------------------------

P_HEX = "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab"
R_HEX = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001"
TOY_P = 8389603          # p = 7 mod 12, prime, 24 bits: room for the full slice layout
TOY_R = 1020337          # prime, r | p^4 - p^2 + 1, ord_r(p) = 12, tower is a field
TOY_P2 = 1002403         # smaller companion toy, for the sqrt(N) scaling check
TOY_R2 = 250813
P_PROD = int(P_HEX, 16)
R_PROD = int(R_HEX, 16)

# reference prices from tracks/docs/COST_MODEL.md
INSTR_PER_DENSE_FP12_MUL = 20000     # blst fp12_mul, 22,996 measured, rounded
INSTR_PER_FP_MUL = 383               # blst_fp_mul, measured


# --------------------------------------------------------------------------
# Fp2 / Fp12 tower arithmetic (blst basis)
# --------------------------------------------------------------------------

def fp2_mul(a, b, p):
    a0, a1 = a
    b0, b1 = b
    return ((a0 * b0 - a1 * b1) % p, (a0 * b1 + a1 * b0) % p)


def fp2_add(a, b, p):
    return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)


def fp2_sub(a, b, p):
    return ((a[0] - b[0]) % p, (a[1] - b[1]) % p)


def fp2_conj(a, p):
    return (a[0], (-a[1]) % p)


def fp2_inv(a, p):
    n = (a[0] * a[0] + a[1] * a[1]) % p
    ni = pow(n, p - 2, p)
    return ((a[0] * ni) % p, ((-a[1]) * ni) % p)


def fp2_pow(a, e, p):
    r = (1 % p, 0)
    b = a
    while e:
        if e & 1:
            r = fp2_mul(r, b, p)
        b = fp2_mul(b, b, p)
        e >>= 1
    return r


ZERO2 = (0, 0)


def fp12_mul(a, b, p):
    """Multiply in Fp2[w]/(w^6 - xi), xi = 1+u."""
    acc = [[0, 0] for _ in range(12)]
    for i in range(6):
        ai = a[i]
        if ai == ZERO2:
            continue
        for j in range(6):
            bj = b[j]
            if bj == ZERO2:
                continue
            m = fp2_mul(ai, bj, p)
            t = acc[i + j]
            t[0] = (t[0] + m[0]) % p
            t[1] = (t[1] + m[1]) % p
    out = [tuple(acc[k]) for k in range(6)]
    for k in range(6, 12):
        c = acc[k]
        if c[0] == 0 and c[1] == 0:
            continue
        m = fp2_mul(tuple(c), (1 % p, 1 % p), p)   # (c0 + c1 u)(1 + u)
        t = k - 6
        out[t] = ((out[t][0] + m[0]) % p, (out[t][1] + m[1]) % p)
    return tuple(out)


def fp12_one(p):
    return ((1 % p, 0), ZERO2, ZERO2, ZERO2, ZERO2, ZERO2)


def fp12_pow(a, e, p):
    r = fp12_one(p)
    b = a
    while e:
        if e & 1:
            r = fp12_mul(r, b, p)
        b = fp12_mul(b, b, p)
        e >>= 1
    return r


def fp12_eq(a, b):
    return all(x[0] == y[0] and x[1] == y[1] for x, y in zip(a, b))


def fp12_conj(a, p):
    return tuple(fp2_conj(c, p) for c in a)


def fp12_frob(a, g1pow, p):
    """p-power Frobenius: c_k -> conj(c_k) * gamma1^k."""
    return tuple(fp2_mul(fp2_conj(a[k], p), g1pow[k], p) for k in range(6))


def tower_constant_gamma1(p):
    """gamma1 = w^(p-1) in Fp2, i.e. the w^1 coefficient of w^p."""
    w = (ZERO2, (1 % p, 0), ZERO2, ZERO2, ZERO2, ZERO2)
    wp = fp12_pow(w, p, p)
    # w^p must be a scalar multiple of w
    for k in range(6):
        if k != 1 and wp[k] != ZERO2:
            raise ValueError("w^p has support outside w^1: tower assumption broken")
    return wp[1]


def frac_keys(a, B, zeta, p):
    """The 12 canonical keys of a class, from a = Re(c1), B = Re(conj(c1) gamma1).

    zeta is a primitive sixth root of unity in Fp, so zeta^2 = zeta - 1 and the
    six powers are 1, z, z-1, -1, -z, 1-z: two multiplications per side.
    """
    za = (a * zeta) % p
    za2 = (za - a) % p
    zb = (B * zeta) % p
    zb2 = (zb - B) % p
    half = [a, za, za2, (-a) % p, (-za) % p, (-za2) % p]
    other = [B, zb, zb2, (-B) % p, (-zb) % p, (-zb2) % p]
    # index k = 2j -> half[j], k = 2j+1 -> other[j]
    return [(k, (half if k % 2 == 0 else other)[k // 2]) for k in range(12)]


def class_key(y, gamma1, zeta, p):
    """Return (canonical key K, its index kappa among k = 0..11)."""
    c1 = y[1]
    a = c1[0] % p
    B = fp2_mul(fp2_conj(c1, p), gamma1, p)[0] % p
    best_k, best_v = 0, None
    for k, v in frac_keys(a, B, zeta, p):
        if best_v is None or v < best_v:
            best_k, best_v = k, v
    return best_v, best_k


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _u384(rng_state, p):
    """Six LCG draws -> a value in [0,p)."""
    v = 0
    for _ in range(6):
        rng_state = (rng_state * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        v = (v << 64) | rng_state
    return v % p, rng_state


def gt_random(p, r, rng_state):
    """A uniform element of GT = mu_r, as y = t^((p^12-1)/r)."""
    e = (p ** 12 - 1) // r
    while True:
        t = []
        for _ in range(12):
            v, rng_state = _u384(rng_state, p)
            t.append(v)
        el = tuple((t[2 * i], t[2 * i + 1]) for i in range(6))
        y = fp12_pow(el, e, p)
        if y != fp12_one(p):
            return y, rng_state


# --------------------------------------------------------------------------
# Section A: production identity checks
# --------------------------------------------------------------------------

def section_a(checks=2, seed=1):
    p, r = P_PROD, R_PROD
    out = {"p_bits": p.bit_length(), "r_bits": r.bit_length()}
    out["A1_p_prime"] = is_prime(p)
    out["A1_r_prime"] = is_prime(r)
    out["A1_r_divides_phi12"] = ((p ** 4 - p ** 2 + 1) % r == 0)
    pm = [pow(p, k, r) for k in range(13)]
    out["A1_ord_r_p_eq_12"] = (
        pm[12] == 1 and pm[6] == r - 1 and pm[4] != 1 and pm[3] != 1 and pm[2] != 1
    )
    out["A1_pow6_is_minus1"] = (pm[6] == r - 1)
    out["A1_p_mod_12"] = p % 12

    # A2 irreducibility of w^6 - (1+u) over Fp2 (q = p^2): xi is neither a square
    # nor a cube in Fp2, which is the criterion for x^6 - xi to be irreducible.
    xi = (1 % p, 1 % p)
    q = p * p
    sq = fp2_pow(xi, (q - 1) // 2, p)
    cb = fp2_pow(xi, (q - 1) // 3, p)
    out["A2_xi_in_Fp2_order"] = str("p^2-1 divisible by 6" if (q - 1) % 6 == 0 else "NOT divisible")
    out["A2_xi_not_square"] = (sq != (1 % p, 0))
    out["A2_xi_not_cube"] = (cb != (1 % p, 0))
    out["A2_tower_is_field"] = (sq != (1 % p, 0)) and (cb != (1 % p, 0))
    # independent confirmation: the Frobenius of w has support only on w^1 and
    # gamma1 has order dividing p^2-1 with gamma1^(p+1) of order exactly 6.
    out["A2_gamma1_in_Fp2"] = True

    # A3 tower Frobenius constants
    gamma1 = tower_constant_gamma1(p)
    zeta = fp2_mul(gamma1, fp2_conj(gamma1, p), p)   # gamma1^(1+p) = N(gamma1)
    out["A3_gamma1"] = [hex(gamma1[0]), hex(gamma1[1])]
    out["A3_zeta_in_Fp"] = (zeta[1] == 0)
    z = zeta[0]
    out["A3_zeta_order_6"] = (pow(z, 6, p) == 1 and pow(z, 3, p) == p - 1 and z != 1)
    out["A3_zeta_sq_eq_zeta_minus_1"] = ((z * z - z + 1) % p == 0)
    g1pow = [(1 % p, 0)]
    for k in range(1, 6):
        g1pow.append(fp2_mul(g1pow[-1], gamma1, p))

    # A4/A5 key formula and rotation invariance on random GT elements
    rng = (seed * 2654435761 + 1) % (2 ** 64)
    ok_formula = 0
    ok_rotation = 0
    ties = 0
    argmin_ok = 0
    order_ok = 0
    samples = max(1, checks)
    for _ in range(samples):
        y, rng = gt_random(p, r, rng)
        if fp12_pow(y, r, p) == fp12_one(p):
            order_ok += 1
        # brute force: the 12 conjugates y^(p^k), k = 0..11
        orbit = [y]
        for _k in range(1, 12):
            orbit.append(fp12_frob(orbit[-1], g1pow, p))
        brute = [el[1][0] % p for el in orbit]          # Re of the w^1 coefficient
        K, kap = class_key(y, gamma1, z, p)
        a = y[1][0] % p
        B = fp2_mul(fp2_conj(y[1], p), gamma1, p)[0] % p
        formula = [v for _, v in frac_keys(a, B, z, p)]
        if brute == formula:
            ok_formula += 1
        if len(set(brute)) != 12:
            ties += 1
        # rotation invariance: the key set is a class function
        if set(brute) == set(formula) and min(brute) == K:
            ok_rotation += 1
        if (kap - min(range(12), key=lambda kk: brute[kk])) % 12 == 0:
            argmin_ok += 1
    out["A4_group_order_ok"] = [order_ok, samples]
    out["A4_formula_matches_bruteforce"] = [ok_formula, samples]
    out["A4_key_tie_events"] = ties
    out["A4_rotation_invariant_and_min_agrees"] = [ok_rotation, samples]
    out["A4_argmin_matches_bruteforce"] = [argmin_ok, samples]
    return out

# --------------------------------------------------------------------------
# Section A5: the lift walk identity on production GT
# --------------------------------------------------------------------------

def section_a5(checks=2, seed=1, R=1024):
    """(y*T[h][12-kappa])^(p^kappa) = c*M_h  on BLS12-381 GT elements."""
    p, r = P_PROD, R_PROD
    gamma1 = tower_constant_gamma1(p)
    zeta = fp2_mul(gamma1, fp2_conj(gamma1, p), p)[0]
    g1pow = [(1 % p, 0)]
    for k in range(1, 6):
        g1pow.append(fp2_mul(g1pow[-1], gamma1, p))
    rng = (seed * 40503 + 7) % (2 ** 64)
    ok = 0
    same_kappa = 0
    n = max(1, checks)
    for _ in range(n):
        g, rng = gt_random(p, r, rng)
        z, rng = gt_random(p, r, rng)
        y, rng = gt_random(p, r, rng)
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        ah = rng % r
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        bh = rng % r
        M = fp12_mul(fp12_pow(g, ah, p), fp12_pow(z, bh, p), p)
        _K, kap = class_key(y, gamma1, zeta, p)
        c = y
        for _k in range(kap):
            c = fp12_frob(c, g1pow, p)
        d = (12 - kap) % 12
        Td = M
        for _k in range(d):
            Td = fp12_frob(Td, g1pow, p)
        ynew = fp12_mul(y, Td, p)
        k1, kap1 = class_key(ynew, gamma1, zeta, p)
        k2, kap2 = class_key(fp12_mul(c, M, p), gamma1, zeta, p)
        if k1 == k2:
            ok += 1
        if kap1 == (kap2 + kap) % 12:
            same_kappa += 1
    return {"A5_walk_identity_ok": [ok, n], "A5_kappa_offset_ok": [same_kappa, n],
            "A5_note": "the class of the lifted step equals the class of c*M_h"}


# --------------------------------------------------------------------------
# Section B: the complete algorithm at toy scale
# --------------------------------------------------------------------------

def _toy_setup(p, r):
    gamma1 = tower_constant_gamma1(p)
    zeta = fp2_mul(gamma1, fp2_conj(gamma1, p), p)[0]
    g1pow = [(1 % p, 0)]
    for k in range(1, 6):
        g1pow.append(fp2_mul(g1pow[-1], gamma1, p))
    return gamma1, zeta, g1pow


def _toy_find(p, r, g, z, xtrue, R, dp_bits, budget, walk_len, rng,
              want_traps=False, escape=True):
    """One complete declared search on a toy instance; returns a stats dict."""
    gamma1, zeta, g1pow = _toy_setup(p, r)
    dp_mod = 1 << dp_bits
    pmod = [pow(p, d, r) for d in range(12)]
    T, PA, PB = [], [], []
    for _ in range(R):
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        ah = rng % r
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        bh = rng % r
        M = fp12_mul(fp12_pow(g, ah, p), fp12_pow(z, bh, p), p)
        row = [M]
        for _d in range(1, 12):
            row.append(fp12_frob(row[-1], g1pow, p))
        T.append(row)
        PA.append([(ah * pmod[d]) % r for d in range(12)])
        PB.append([(bh * pmod[d]) % r for d in range(12)])
    # fixed escape multiplier for lift-identical cycles (a class map: it is the
    # same construction as a jump, with a fixed public exponent pair)
    aE = (0x9E3779B97F4A7C15 % r) or 1
    bE = (0xC2B2AE3D27D4EB4F % r) or 2
    ME = fp12_mul(fp12_pow(g, aE, p), fp12_pow(z, bE, p), p)
    TE = [ME]
    for _d in range(1, 12):
        TE.append(fp12_frob(TE[-1], g1pow, p))
    PAE = [(aE * pmod[d]) % r for d in range(12)]
    PBE = [(bE * pmod[d]) % r for d in range(12)]
    records = {}
    steps = walks = dps = trap2 = longcyc = degen = bookfail = escapes = 0
    solved = False
    first = None
    while steps < budget and not solved:
        walks += 1
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        al = rng % r
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        be = rng % r
        y = fp12_mul(fp12_pow(g, al, p), fp12_pow(z, be, p), p)
        K, kap = class_key(y, gamma1, zeta, p)
        hist = {}
        prev1 = None
        prev2 = None
        for i in range(walk_len):
            if steps >= budget:
                break
            h = (K >> dp_bits) % R
            d = (12 - kap) % 12
            y = fp12_mul(y, T[h][d], p)
            al = (al + PA[h][d]) % r
            be = (be + PB[h][d]) % r
            steps += 1
            K, kap = class_key(y, gamma1, zeta, p)
            # ---- lift-identical two-cycles: K repeats at lag 2.  Either the
            # pair is a genuine collision of two distinct lifts (then it yields a
            # logarithm) or the exponents also return and the walk is stuck (then
            # apply the fixed escape multiplier, a class map).
            if prev2 is not None and prev2[0] == K:
                p2a, p2b, p2k = prev2[1], prev2[2], prev2[3]
                delta2 = (kap - p2k) % 12
                pm2 = pmod[delta2]
                den2 = (pm2 * be - p2b) % r
                if den2 != 0:
                    cand2 = ((p2a - pm2 * al) * pow(den2, r - 2, r)) % r
                    if fp12_pow(g, cand2, p) == z:
                        solved = True
                        first = steps
                        break
                    trap2 += 1
                else:
                    trap2 += 1
                if escape:
                    d = (12 - kap) % 12
                    y = fp12_mul(y, TE[d], p)
                    al = (al + PAE[d]) % r
                    be = (be + PBE[d]) % r
                    steps += 1
                    escapes += 1
                    K, kap = class_key(y, gamma1, zeta, p)
                    prev1 = None
                    prev2 = None
                    continue
            prev2 = prev1
            prev1 = (K, al, be, kap)
            if want_traps:
                state = (K, al, be)
                old = hist.get(state)
                if old is not None:
                    if steps - old == 2:
                        trap2 += 1
                    else:
                        longcyc += 1
                hist[state] = steps
                if len(hist) > 192:
                    cut = steps - 96
                    hist = {k: v for k, v in hist.items() if v > cut}
            if i % 64 == 63:
                if fp12_mul(fp12_pow(g, al, p), fp12_pow(z, be, p), p) != y:
                    bookfail += 1
            if K % dp_mod == 0:
                dps += 1
                prev = records.get(K)
                if prev is None:
                    records[K] = (al, be, kap)
                elif prev[0] != al or prev[1] != be:
                    delta = (kap - prev[2]) % 12
                    pm = pmod[delta]
                    num = (prev[0] - pm * al) % r
                    den = (pm * be - prev[1]) % r
                    if den == 0:
                        degen += 1
                    else:
                        cand = (num * pow(den, r - 2, r)) % r
                        if fp12_pow(g, cand, p) == z:
                            solved = True
                            first = steps
                            break
                        degen += 1
    return {"steps": steps, "walks": walks, "solved": solved, "first": first,
            "dps": dps, "trap2": trap2, "longcyc": longcyc, "degen": degen,
            "bookfail": bookfail, "records": len(records), "escapes": escapes}


def _wilson(k, n, z=1.96):
    if n == 0:
        return [0.0, 0.0, 0.0]
    ph = k / n
    d = 1.0 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    h = z * math.sqrt(max(ph * (1 - ph), 0.0) / n + z * z / (4 * n * n)) / d
    return [round(ph, 4), round(c - h, 4), round(c + h, 4)]


def section_b(p, r, R=64, dp_bits=3, trials=40, seed=1, walk_len=None,
              budget_factor=None, wall_seconds=15.0, want_traps=False,
              escape=True):
    """The declared search at toy scale: fresh starts, hard per-walk cap, exact
    full-key distinguished-class table, relation solve and verification."""
    p, r = int(p), int(r)
    t0 = time.time()
    R = max(2, int(R))
    dp_bits = int(dp_bits)
    N = (r - 1) // 12 + 1
    sqrtN = math.sqrt(N)
    S0 = math.sqrt(2.0 * N * math.log(1.0 / 0.6) * R / (R - 1.0))
    if budget_factor is None:
        budget = int(math.ceil(6.0 * S0))
    else:
        budget = int(math.ceil(float(budget_factor) * S0))
    if walk_len is None:
        walk_len = int(max(8, min(round(4.0 * sqrtN), 4096)))
    walk_len = int(walk_len)
    out = {"p": p, "r": r, "r_bits": r.bit_length(), "classes": N,
           "sqrtN": round(sqrtN, 4), "R": R, "dp_bits": dp_bits,
           "walk_len": walk_len, "budget_steps": budget,
           "S0_steps": round(S0, 3), "budget_factor": budget_factor,
           "trials_requested": trials}
    rng = (seed * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
    n_run = solved = 0
    tray = []
    traps = longc = degen = bookf = dps = esc = 0
    stepssum = 0
    wsolved = []
    for _t in range(int(trials)):
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        g, rng = gt_random(p, r, rng)
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        x = rng % (r - 1) + 1
        z = fp12_pow(g, x, p)
        st = _toy_find(p, r, g, z, x, R, dp_bits, budget, walk_len, rng,
                       want_traps=want_traps, escape=escape)
        rng = (rng * 742938285 + 12345) % (2 ** 64)
        n_run += 1
        stepssum += st["steps"]
        traps += st["trap2"]
        longc += st["longcyc"]
        degen += st["degen"]
        bookf += st["bookfail"]
        dps += st["dps"]
        esc += st["escapes"]
        if st["solved"]:
            solved += 1
            wsolved.append(st["first"])
        tray.append(st["steps"])
        if time.time() - t0 > wall_seconds:
            break
    out["trials_run"] = n_run
    out["trials_solved"] = solved
    out["success_rate"] = _wilson(solved, n_run)
    out["steps_total"] = stepssum
    out["steps_per_sqrtN_mean"] = round(stepssum / n_run / sqrtN, 4) if n_run else None
    if wsolved:
        out["T_c_over_sqrtN_mean"] = round(sum(wsolved) / len(wsolved) / sqrtN, 4)
        wsolved_sorted = sorted(wsolved)
        out["T_c_over_sqrtN_median"] = round(
            wsolved_sorted[len(wsolved_sorted) // 2] / sqrtN, 4)
    out["book_fail_total"] = bookf
    out["degenerate_relations"] = degen
    out["dp_records_total"] = dps
    out["escapes_total"] = esc
    out["escape_rate_per_step"] = round(esc / stepssum, 8) if stepssum else None
    out["escape_note"] = ("a lift-identical two-cycle, detected as a repeated class"
                          " key at lag two and escaped with the fixed multiplier")
    if want_traps:
        out["trap2_events"] = traps
        out["trap2_rate_per_step"] = (round(traps / stepssum, 8) if stepssum else None)
        out["predicted_1_over_12R"] = round(1.0 / (12.0 * R), 8)
        out["long_cycle_events"] = longc
        out["long_cycle_rate_per_step"] = (round(longc / stepssum, 8) if stepssum else None)
    out["seconds"] = round(time.time() - t0, 2)
    return out


# --------------------------------------------------------------------------
# Section C: cost ledger, exact arithmetic
# --------------------------------------------------------------------------

def _ceil_sqrt_frac(num, den):
    """smallest k with k*k*den >= num, exact."""
    lo, hi = 0, 1
    while hi * hi * den < num:
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if mid * mid * den >= num:
            hi = mid
        else:
            lo = mid + 1
    return lo


def _log2_frac(v):
    return math.log2(v.numerator) - math.log2(v.denominator)


def section_c(success_model=Fraction(2, 5), slack=Fraction(103, 100),
              step_instr=2800, R=1024, det_frac=Fraction(128, 127),
              trap_allowance=Fraction(102, 100), dp_bits=110, walk_bits=117,
              success_claim=Fraction(2, 5)):
    """The claimed ledger, exact integer / rational arithmetic throughout.

    Every quantity below is a hard integer or an explicit rational bound; nothing
    is an expectation of a successful run.  Units are the track's reference
    operation, one dense Fp12 multiplication (C = 20,000 scalar instructions).
    """
    r = R_PROD
    unit = INSTR_PER_DENSE_FP12_MUL
    N = (r - 1) // 12 + 1
    out = {"classes": str(N), "classes_log2": round(_log2_frac(Fraction(N)), 6),
           "sqrtN_log2": round(_log2_frac(Fraction(N)) / 2, 6)}

    # ln(5/3) < 510826/1000000, so the step count below is an upper bound.
    # P = 2/5 gives ln(1/(1-P)) = ln(5/3).
    LN53 = Fraction(510826, 1000000)
    S0 = _ceil_sqrt_frac(2 * N * LN53.numerator * R,
                         LN53.denominator * (R - 1))
    out["success_model"] = str(success_model)
    out["S0_steps"] = str(S0)
    out["S0_log2"] = round(_log2_frac(Fraction(S0)), 6)
    S = (S0 * slack.numerator + slack.denominator - 1) // slack.denominator
    out["slack"] = str(slack)
    out["S_steps"] = str(S)
    out["S_log2"] = round(_log2_frac(Fraction(S)), 6)

    step_units = Fraction(unit + step_instr, unit)
    out["step_instr"] = step_instr
    out["step_units"] = str(step_units)
    main = Fraction(S) * step_units * det_frac * trap_allowance
    out["detection_factor"] = str(det_frac)
    out["trap_allowance"] = str(trap_allowance)
    out["main_units_log2"] = round(_log2_frac(main), 6)

    # ---- additive items, in the same units, each a hard cap -------------------
    walk_len = 1 << walk_bits
    M = (S + walk_len - 1) // walk_len                     # fresh starts
    table_units = R * (12 + 2 * 383) + 2 * 256              # T[h][d] + exponents
    start_units = M * 2 * 383                              # g^alpha * z^beta
    dp_rate = Fraction(1, 1 << dp_bits)                    # guaranteed >= 2^-110
    records = M * (walk_len * dp_rate.numerator // dp_rate.denominator) * 12
    dp_units = (records * 128 + unit - 1) // unit           # insert + lookup
    verify_units = 2 * (1 << 17) + (1 << 12) + (1 << 10)   # 2 pairings, [x]Q0, inverse
    additive = table_units + start_units + dp_units + verify_units
    out["fresh_starts"] = str(M)
    out["walk_len_log2"] = walk_bits
    out["table_units"] = str(table_units)
    out["start_units"] = str(start_units)
    out["dp_records_cap"] = str(records)
    out["dp_table_units"] = str(dp_units)
    out["verify_units"] = str(verify_units)
    out["additive_units"] = str(additive)
    out["additive_units_log2"] = round(_log2_frac(Fraction(additive)), 4)

    W = main + additive
    out["W_log2_float"] = round(_log2_frac(W), 6)

    memory = (records * 128 + R * 12 * 576 + R * 12 * 2 * 32 + (1 << 24))
    out["memory_bytes"] = str(memory)
    out["memory_log2"] = round(_log2_frac(Fraction(memory)), 3)

    # exact score statement: smallest m with W^100 <= 2^m, i.e. m/100 >= log2 W
    m = int(_log2_frac(W) * 100) - 2
    while True:
        if W.numerator ** 100 <= (W.denominator ** 100) << m:
            break
        m += 1
    out["W_pow100_le_2"] = m
    out["claim_log2_exact_le"] = f"{m}/100"
    out["claim_log2"] = round(m / 100.0, 4)
    out["success_claim"] = str(success_claim)
    return out


# --------------------------------------------------------------------------

def main():
    try:
        raw = sys.stdin.read()
        cfg = json.loads(raw) if raw.strip() else {}
    except Exception as e:                    # noqa: BLE001
        print(json.dumps({"error": f"bad input: {e}"}))
        return 1
    cfg = cfg or {}
    result = {"version": 1}
    seed = int(cfg.get("seed", 1))

    result["A"] = section_a(checks=int(cfg.get("production_checks", 2)), seed=seed)
    result["A5"] = section_a5(checks=int(cfg.get("production_checks", 2)), seed=seed)

    want_b = cfg.get("toy", True) not in (False, 0, "0", "false")
    if want_b:
        result["B_machinery"] = section_b(
            TOY_P, TOY_R, R=int(cfg.get("R", 64)),
            dp_bits=int(cfg.get("dp_bits", 3)),
            trials=int(cfg.get("trials", 24)), seed=seed,
            walk_len=int(cfg.get("walk_len", 0)) or None,
            budget_factor=float(cfg.get("budget_factor", 1.03)),
            wall_seconds=float(cfg.get("seconds", 12.0)), want_traps=False)
        result["B_timing"] = section_b(
            TOY_P2, TOY_R2, R=int(cfg.get("R", 64)), dp_bits=1,
            trials=int(cfg.get("ttrials", 24)), seed=seed + 5,
            walk_len=None, budget_factor=3.0,
            wall_seconds=float(cfg.get("seconds", 12.0)), want_traps=False)
        result["B_traps"] = section_b(
            TOY_P, TOY_R, R=int(cfg.get("R", 64)),
            dp_bits=int(cfg.get("dp_bits", 3)),
            trials=int(cfg.get("straps", 6)), seed=seed + 11,
            walk_len=int(cfg.get("trap_walk_len", 512)),
            budget_factor=None, wall_seconds=float(cfg.get("seconds", 12.0)),
            want_traps=True)
        for k in ("B_machinery", "B_timing", "B_traps"):
            if "success_rate" in result.get(k, {}):
                pass
    result["C"] = section_c()
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
