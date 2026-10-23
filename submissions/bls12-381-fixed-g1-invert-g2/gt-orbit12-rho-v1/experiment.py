#!/usr/bin/env python3
"""experiment.py -- orbit-12 Frobenius-class distinguished-point rho in GT.

Candidate: submissions/bls12-381-fixed-g1-invert-g2/gt-orbit12-rho-v1
Track:     bls12-381-fixed-g1-invert-g2 (fixed-argument pairing inversion)

Standard library only.  One JSON object on stdin, one JSON object on stdout.

Section A (production scale, on BLS12-381 itself, blst tower)
    Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - (1+u)).
    A1  p, r prime; r | Phi_12(p) = p^4 - p^2 + 1; ord_r(p) = 12 exactly
        (p^6 = -1, p^4 != 1, p^3 != 1, p^2 != 1).
    A2  w^6 - (1+u) is irreducible over Fp2, so the tower is the field F_{p^12}.
    A3  Frobenius constant gamma1 = w^(p-1) in Fp2, zeta = gamma1^(p+1) = N(gamma1)
        in Fp is a primitive sixth root of unity (zeta^6 = 1, zeta^2 = zeta - 1).
    A4  The 12-key canonicalisation formula
            a = Re(c_1(y)),  B = Re(conj(c_1(y)) * gamma1),
            key_{2j}(y) = a*zeta^j,  key_{2j+1}(y) = B*zeta^j   (j = 0..5)
        is checked against brute-force Frobenius y -> y^(p^k) for k = 0..11 on
        random GT elements (y = t^((p^12-1)/r), y != 1, y^r = 1).
    A5  Rotation invariance: the 12-element key *set* of y and of y^p coincide,
        so min over the set is a class function; argmin(y^p) = argmin(y) - 1.

Section B (toy scale, same code path as production)
    A toy prime p with p = 7 (mod 12) and a prime r | p^4 - p^2 + 1, ord_r(p) = 12,
    built in the same tower.  GT = mu_r.  The complete algorithm runs there:
    lift/kappa class walk, T[h][d] = M_h^(p^d) jump table, exponent bookkeeping
    g^alpha z^beta == y checked at every distinguished point, the declared
    dictionary (bit-slice-disjoint DP test / bucket / walk index / probe step,
    odd multiplicative probe, four-probe cap, drop on fourth miss, hard cap) and
    the 12-candidate collision solve.  Measured: steps/sqrt(N) at first
    collision, success rate, bookkeeping mismatches, dictionary statistics and
    the realised bucket distribution.

Section C
    The cost ledger in exact integer / rational arithmetic.

Toy statistics are evidence at toy scope only; the production claim is analytic.
"""

import json
import sys
import time
from fractions import Fraction

# --------------------------------------------------------------------------
# parameters
# --------------------------------------------------------------------------

P_HEX = "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab"
R_HEX = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001"
TOY_P = 1002403          # p = 7 mod 12, prime
TOY_R = 250813           # prime, r | p^4 - p^2 + 1, ord_r(p) = 12, tower is a field
TOY_P2 = 1009483         # smaller companion toy, for the sqrt(N) scaling check
TOY_R2 = 51109
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
# Section B: the complete algorithm at toy scale
# --------------------------------------------------------------------------

def _tower_is_field(p):
    """w^6 - (1+u) irreducible over Fp2: xi is neither a square nor a cube."""
    xi = (1 % p, 1 % p)
    q = p * p
    return (fp2_pow(xi, (q - 1) // 2, p) != (1 % p, 0)) and \
           (fp2_pow(xi, (q - 1) // 3, p) != (1 % p, 0))


def find_toy(lo_hi=(101, 4000), r_lo=40000, r_hi=4000000):
    """A toy prime p = 7 mod 12 and a prime r | p^4-p^2+1 with ord_r(p) = 12."""
    primes = []
    for n in range(lo_hi[0], lo_hi[1], 12):
        p = n + (7 - n % 12) % 12
        if p < lo_hi[0]:
            continue
        if is_prime(p) and p % 12 == 7:
            primes.append(p)
    for p in primes:
        M = p ** 4 - p ** 2 + 1
        # small trial division
        m = M
        d = 2
        while d * d <= m and d < 200000:
            if m % d == 0:
                while m % d == 0:
                    m //= d
                if r_lo <= d <= r_hi and is_prime(d) and pow(p, 6, d) == d - 1 and pow(p, 4, d) != 1:
                    if _tower_is_field(p):
                        return p, d
            d += 1 if d == 2 else 2
        if r_lo <= m <= r_hi and is_prime(m) and pow(p, 6, m) == m - 1 and pow(p, 4, m) != 1:
            if _tower_is_field(p):
                return p, m
    return None, None


def toy_instances(p, r, n):
    """n - 1 toy GT generators g and targets z = g^x."""
    e = (p ** 12 - 1) // r
    res = []
    h = 12345
    attempts = 0
    while len(res) < n and attempts < 200 * n + 200:
        attempts += 1
        h = (h * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        el = []
        for i in range(12):
            h = (h * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
            el.append(h % p)
        t = tuple((el[2 * i], el[2 * i + 1]) for i in range(6))
        g = fp12_pow(t, e, p)
        if g == fp12_one(p) or fp12_pow(g, r, p) != fp12_one(p):
            continue
        h = (h * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        x = h % r
        if x == 0:
            continue
        res.append((g, x, fp12_pow(g, x, p)))
    return res


def section_b(p, r, R=8, wander_max=None, trials=12, walks=4, dp_bits=2,
              bucket_bits=3, slot_bits=6, h_bits=3, seed=1, cap_factor=6,
              budget_factor=None):
    """Faithful run of the declared algorithm on GT = mu_r of a toy tower."""
    p = int(p)
    r = int(r)
    gamma1 = tower_constant_gamma1(p)
    zeta2 = fp2_mul(gamma1, fp2_conj(gamma1, p), p)
    zeta = zeta2[0]
    total_bits = p.bit_length()
    dp_bits = int(dp_bits)
    # bit-slice layout, mirroring production: [0,dp) DP test, then the 128-bit
    # class identity, then bucket, walk index and probe step on disjoint slices.
    buck_lo = dp_bits
    buck_hi = min(total_bits, buck_lo + bucket_bits)
    h_lo = buck_hi
    h_hi = min(total_bits, h_lo + h_bits)
    step_lo = h_hi
    step_hi = min(total_bits, step_lo + slot_bits - 1)
    key_lo = step_hi
    key_hi = total_bits
    if key_lo >= key_hi:
        raise ValueError("toy p too small for the declared disjoint slice layout")
    C_BUCK = (0x9E3779B97F4A7C15C2B2AE3D27D4EB4F % (2 ** 64)) | 1

    def bucket_of(K):
        v = (K >> buck_lo) & ((1 << (buck_hi - buck_lo)) - 1)
        return (v * C_BUCK) % (1 << slot_bits)

    def hindex_of(K):
        return (K >> h_lo) & ((1 << (h_hi - h_lo)) - 1)

    def step_of(K):
        m = (1 << (slot_bits - 1)) - 1
        return (1 + 2 * ((K >> step_lo) & m)) % (1 << slot_bits) | 1

    def ident_of(K):
        return (K >> key_lo) & ((1 << (key_hi - key_lo)) - 1)

    stats = {
        "p": p, "r": r, "r_bits": r.bit_length(),
        "classes": (r - 1) // 12 + 1, "dp_bits": dp_bits,
        "R": R, "walks": walks, "trials": trials,
        "slot_bits": slot_bits, "bucket_bits": bucket_bits,
        "layout": {"dp": [0, dp_bits], "ident": [key_lo, key_hi],
                   "bucket": [buck_lo, buck_hi], "h": [h_lo, h_hi],
                   "step": [step_lo, step_lo + slot_bits - 1]},
    }
    N = (r - 1) // 12 + 1
    import math
    sqrtN = math.sqrt(N)

    # jump table: M_h = g^ah z^bh, T[h][d] = M_h^(p^d)
    def build_tables(g, z, rng):
        T = []
        PA = []
        PB = []
        pmod = [pow(p, d, r) for d in range(12)]
        for _ in range(R):
            rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
            ah = rng % r
            rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
            bh = rng % r
            M = fp12_mul(fp12_pow(g, ah, p), fp12_pow(z, bh, p), p)
            g1pow = [(1 % p, 0)]
            for k in range(1, 6):
                g1pow.append(fp2_mul(g1pow[-1], gamma1, p))
            row = [M]
            for d in range(1, 12):
                row.append(fp12_frob(row[-1], g1pow, p))
            T.append(row)
            PA.append([(ah * pmod[d]) % r for d in range(12)])
            PB.append([(bh * pmod[d]) % r for d in range(12)])
        return T, PA, PB, rng

    CHECK = 24

    def run_trial(g, xtrue, z, rng, cap_steps):
        T, PA, PB, rng = build_tables(g, z, rng)
        nslots = 1 << slot_bits
        slots = [None] * nslots
        records = []
        probes = inserts = drops = lookups = dp_count = steps = 0
        escapes = 0
        book_fail = 0
        bucket_hist = {}
        solved = None
        starts = []
        for _ in range(walks):
            rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
            al = rng % r
            rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
            be = rng % r
            y = fp12_mul(fp12_pow(g, al, p), fp12_pow(z, be, p), p)
            K, kap = class_key(y, gamma1, zeta, p)
            starts.append([y, al, be, K, kap, []])
        active = [True] * walks
        while steps < cap_steps and solved is None:
            moved = False
            for wi in range(walks):
                if solved is not None:
                    break
                y, al, be, K, kap, hist = starts[wi]
                h = hindex_of(K) % R
                d = (12 - kap) % 12
                y = fp12_mul(y, T[h][d], p)
                al = (al + PA[h][d]) % r
                be = (be + PB[h][d]) % r
                steps += 1
                moved = True
                K, kap = class_key(y, gamma1, zeta, p)
                # ---- cycle escape: a repeat inside the recent window means the
                # walk is trapped; squaring is a class map and breaks the cycle.
                if K in hist:
                    y = fp12_mul(y, y, p)
                    al = (2 * al) % r
                    be = (2 * be) % r
                    steps += 1
                    K, kap = class_key(y, gamma1, zeta, p)
                    hist = []
                    escapes += 1
                else:
                    hist.append(K)
                    if len(hist) > CHECK:
                        hist.pop(0)
                if K % (1 << dp_bits) != 0:
                    starts[wi] = [y, al, be, K, kap, hist]
                    continue
                dp_count += 1
                if fp12_mul(fp12_pow(g, al, p), fp12_pow(z, be, p), p) != y:
                    book_fail += 1
                ident = ident_of(K)
                b = bucket_of(K)
                bucket_hist[b] = bucket_hist.get(b, 0) + 1
                stepv = step_of(K)
                pos = b
                found = free = None
                for _pr in range(4):
                    probes += 1
                    cur = slots[pos]
                    if cur is None:
                        free = pos
                        break
                    if records[cur][0] == ident:
                        found = cur
                        break
                    pos = (pos + stepv) % nslots
                if found is not None:
                    lookups += 1
                    a1, b1 = records[found][1], records[found][2]
                    for j in range(12):
                        pj = pow(p, j, r)
                        num = (al - a1 * pj) % r
                        den = (b1 * pj - be) % r
                        if den == 0:
                            continue
                        cand = num * pow(den, r - 2, r) % r
                        if fp12_pow(g, cand, p) == z:
                            solved = cand
                            break
                elif free is not None:
                    records.append([ident, al, be])
                    slots[free] = len(records) - 1
                    inserts += 1
                else:
                    drops += 1
                starts[wi] = [y, al, be, K, kap, hist]
            if not moved:
                break
        return {
            "steps": steps, "solved": solved == xtrue, "solved_value": solved,
            "dp": dp_count, "probes": probes, "inserts": inserts,
            "drops": drops, "book_fail": book_fail, "escapes": escapes,
            "lookups": lookups,
            "buckets_used": len(bucket_hist),
            "bucket_max": max(bucket_hist.values()) if bucket_hist else 0,
        }

    rng = (seed * 2654435761 + 12345) % (2 ** 64)
    results = []
    inst = toy_instances(p, r, trials)
    if budget_factor is not None:
        cap = max(2, int(budget_factor * sqrtN))
        stats["budget_factor"] = budget_factor
        stats["budget_steps"] = cap
    else:
        cap = int(cap_factor * sqrtN * max(1, walks))
    t0 = time.time()
    for g, xtrue, z in inst:
        res = run_trial(g, xtrue, z, rng, cap)
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        results.append(res)
        if time.time() - t0 > 40:
            break
    ok = [x for x in results if x["solved"]]
    steps = sorted(x["steps"] for x in ok)
    stats["trials_run"] = len(results)
    stats["trials_solved"] = len(ok)
    stats["success_rate"] = (Fraction(len(ok), len(results)) if results else 0)
    stats["steps_over_sqrtN_median"] = (
        (Fraction(steps[len(steps) // 2], 1) / Fraction(int(sqrtN * 1000), 1000)) if steps else None)
    stats["steps_median"] = steps[len(steps) // 2] if steps else None
    stats["sqrtN"] = round(sqrtN, 3)
    stats["steps_list"] = [x["steps"] for x in results]
    stats["book_fail_total"] = sum(x["book_fail"] for x in results)
    stats["dp_total"] = sum(x["dp"] for x in results)
    stats["probes_total"] = sum(x["probes"] for x in results)
    stats["inserts_total"] = sum(x["inserts"] for x in results)
    stats["drops_total"] = sum(x["drops"] for x in results)
    stats["escapes_total"] = sum(x["escapes"] for x in results)
    stats["lookups_total"] = sum(x["lookups"] for x in results)
    stats["buckets_used_max"] = max((x["buckets_used"] for x in results), default=0)
    stats["bucket_max_occupancy"] = max((x["bucket_max"] for x in results), default=0)
    stats["slots"] = 1 << slot_bits
    stats["ident_bits"] = key_hi - key_lo
    stats["theory_sqrt_pi2_sqrtN"] = round(1.2533 * sqrtN, 1)
    if budget_factor is not None:
        stats["budget_success_rate"] = round(len(ok) / max(1, len(results)), 4)
    return stats


# --------------------------------------------------------------------------
# Section C: cost ledger, exact arithmetic
# --------------------------------------------------------------------------

def _ceil_sqrt(num, den):
    """smallest k with k^2 * den >= num, exactly."""
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


def _log2_lower_bound(v):
    """math.log2 of a positive Fraction, in float (reporting only)."""
    import math
    return math.log2(v.numerator) - math.log2(v.denominator)


def section_c(success_claim=Fraction(40, 100), success_model=Fraction(42, 100),
              efficiency=Fraction(113, 100), cycle_factor=Fraction(101, 100),
              step_instr=2600, R=1024,
              dp_bits=100, slot_bits=30, record_bytes=96, handle_bytes=16):
    """The claimed ledger, exact integer / rational arithmetic throughout."""
    p, r = P_PROD, R_PROD
    N = (r - 1) // 12 + 1                        # classes of size 12, plus {1}
    out = {"classes": str(N), "classes_log2": round(_log2_lower_bound(Fraction(N)), 4)}

    # ln(1/(1-s)) <= s + s^2/(2(1-s))  (Taylor remainder bound), all exact
    s = success_model

    def ln_inv_upper(t):
        return t + t * t / (2 * (1 - t))

    L = ln_inv_upper(s)
    out["ln_inv_bound_model"] = str(L)
    out["ln_inv_bound_claim"] = str(ln_inv_upper(success_claim))
    # k0 = ceil(sqrt(2 N L))
    k0 = _ceil_sqrt(2 * N * L.numerator, L.denominator)
    out["k0_steps"] = str(k0)
    out["k0_log2"] = round(_log2_lower_bound(Fraction(k0)), 6)
    # the model success probability at k0 is at least the claim
    out["model_success_at_k0_ge_claim"] = (Fraction(k0 * k0, 1) / (2 * N) >= ln_inv_upper(success_claim))
    out["model_success_flag_note"] = "1-exp(-k^2/(2N)) >= claim because k^2/(2N) >= ln(1/(1-claim))"

    out["efficiency_charge"] = str(efficiency)
    out["k0_charged_steps"] = str(-(-k0 * efficiency.numerator // efficiency.denominator))
    out["k0_charged_log2"] = round(_log2_lower_bound(Fraction(
        -(-k0 * efficiency.numerator // efficiency.denominator))), 6)
    k0 = -(-k0 * efficiency.numerator // efficiency.denominator)
    total_steps = -(-k0 * cycle_factor.numerator // cycle_factor.denominator)
    out["cycle_factor"] = str(cycle_factor)
    out["total_steps"] = str(total_steps)
    step_units = Fraction(INSTR_PER_DENSE_FP12_MUL + step_instr, INSTR_PER_DENSE_FP12_MUL)
    out["step_instr_extra"] = step_instr
    out["step_units"] = str(step_units)
    main = total_steps * step_units

    dps = -(-total_steps // (1 << dp_bits))
    slots = 1 << slot_bits
    dict_instr = 2 * dps * 4 * 120            # lookup + insert, <=4 probes, 120 instr/probe
    init_instr = slots * 2                    # clear one handle slot
    prep_instr = R * 12 * INSTR_PER_FP_MUL + R * 2 * 400 * INSTR_PER_FP_MUL
    verify_instr = 12 * 400 * INSTR_PER_FP_MUL + 200 * INSTR_PER_DENSE_FP12_MUL
    overhead_units = Fraction(dict_instr + init_instr + prep_instr + verify_instr,
                              INSTR_PER_DENSE_FP12_MUL)
    W = main + overhead_units
    out["dictionary_instr"] = str(dict_instr)
    out["cold_start_instr"] = str(init_instr)
    out["preprocessing_instr"] = str(prep_instr)
    out["verification_instr"] = str(verify_instr)
    out["overhead_units"] = str(overhead_units)
    out["main_units_log2"] = round(_log2_lower_bound(main), 4)
    out["W_log2_float"] = round(_log2_lower_bound(W), 6)
    out["dp_stored"] = str(dps)
    mem = slots * handle_bytes + dps * record_bytes
    out["memory_bytes"] = str(mem)
    out["memory_log2"] = round(_log2_lower_bound(Fraction(mem)), 2)

    # exact score statement: smallest m with W^100 <= 2^m
    m = int(_log2_lower_bound(W) * 100) - 2
    while True:
        if int(W.numerator) ** 100 <= (int(W.denominator) ** 100) << m:
            break
        m += 1
    out["W_pow100_le_2"] = m
    out["claim_log2_exact_le"] = f"{m}/100"
    out["claim_log2"] = round(m / 100.0, 4)
    return out


# --------------------------------------------------------------------------

def main():
    try:
        raw = sys.stdin.read()
        cfg = json.loads(raw) if raw.strip() else {}
    except Exception as e:            # noqa: BLE001
        print(json.dumps({"error": f"bad input: {e}"}))
        return 1
    cfg = cfg or {}
    t0 = time.time()
    result = {"version": 1}

    prod_checks = int(cfg.get("production_checks", 2))
    result["A"] = section_a(checks=prod_checks, seed=int(cfg.get("seed", 1)))

    want_b = cfg.get("toy", True) not in (False, 0, "0", "false")
    if want_b:
        p_toy = cfg.get("p")
        r_toy = cfg.get("r")
        if p_toy is None or r_toy is None:
            p_toy, r_toy = TOY_P, TOY_R
        if p_toy is None:
            result["B"] = {"error": "no toy pair found"}
        else:
            result["B"] = section_b(
                p_toy, r_toy,
                R=int(cfg.get("R", 8)),
                trials=int(cfg.get("trials", 12)),
                walks=int(cfg.get("walks", 4)),
                dp_bits=int(cfg.get("dp_bits", 2)),
                bucket_bits=int(cfg.get("bucket_bits", 3)),
                slot_bits=int(cfg.get("slot_bits", 6)),
                h_bits=int(cfg.get("h_bits", 3)),
                seed=int(cfg.get("seed", 1)),
                cap_factor=float(cfg.get("cap_factor", 6)),
            )
    if want_b and cfg.get("scaling", True) not in (False, 0, "0", "false"):
        try:
            result["B_small"] = section_b(
                TOY_P2, TOY_R2,
                R=int(cfg.get("R", 8)), trials=int(cfg.get("trials", 12)),
                walks=int(cfg.get("walks", 4)), dp_bits=int(cfg.get("dp_bits", 2)),
                bucket_bits=int(cfg.get("bucket_bits", 3)),
                slot_bits=int(cfg.get("slot_bits", 6)), h_bits=int(cfg.get("h_bits", 3)),
                seed=int(cfg.get("seed", 1)) + 7, cap_factor=float(cfg.get("cap_factor", 6)))
        except Exception as e:                       # noqa: BLE001
            result["B_small"] = {"error": str(e)}
    if want_b:
        # success-rate evidence at the claimed budget: the model's k0 for the
        # 0.42 and 0.45 success targets, in units of sqrt(N_classes)
        for bf, tag in ((1.0415, "budget_0.42"), (1.0934, "budget_0.45"), (1.10, "budget_0.455")):
            try:
                result[tag] = section_b(
                    TOY_P, TOY_R, R=int(cfg.get("R", 8)), trials=int(cfg.get("btrials", 120)),
                    walks=int(cfg.get("walks", 4)), dp_bits=int(cfg.get("dp_bits", 2)),
                    bucket_bits=int(cfg.get("bucket_bits", 3)),
                    slot_bits=int(cfg.get("slot_bits", 6)), h_bits=int(cfg.get("h_bits", 3)),
                    seed=int(cfg.get("seed", 1)) + 31, cap_factor=6, budget_factor=bf)
                result[tag] = {k: result[tag][k] for k in
                               ("classes", "sqrtN", "budget_factor", "budget_steps",
                                "trials_run", "trials_solved", "budget_success_rate",
                                "book_fail_total")}
            except Exception as e:                    # noqa: BLE001
                result[tag] = {"error": str(e)}
    result["C"] = section_c()
    result["elapsed_seconds"] = round(time.time() - t0, 2)
    print(json.dumps(result, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
