#!/usr/bin/env python3
"""experiment.py -- O(1)-memory staggered-checkpoint Pollard rho in GT for fixed-G1 pairing inversion.

Candidate: submissions/bls12-381-fixed-g1-invert-g2/gt-stagger-rho-v1
Track:     bls12-381-fixed-g1-invert-g2 (fixed-argument pairing inversion into G2)

Standard library only.  One JSON object on stdin, one JSON object on stdout.

Section A (production scale, on BLS12-381 itself)
    A1  p and r prime; r | Phi_12(p); ord_r(p) = 12 exactly with p^6 = -1.
    A2  the blst tower Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6-(1+u)) is a field
        (xi = 1+u neither a square nor a cube) and gamma1 = w^(p-1) lies in Fp2.
    A3  a random element of GT = mu_r computed as t^((p^12-1)/r) is nonidentity and
        has order r, so GT is cyclic of prime order r.
    A4  the collision identity used by the attack, on real 381-bit GT elements:
        with z = g^x and X_i = g^(a_i) z^(b_i), two equal states with different b
        return x = (a1-a2)/(b2-b1) mod r.  Checked by constructing two equal states
        from chosen exponents and recovering x.
    A5  sampled step-table entries are nonidentity and of order r.

Section B (schedule assertions and toy scale)
    B1  the live-checkpoint count is ENUMERATED, not asserted: every hare position
        in [1, 2^20] is scanned, plus the probe positions 1024, 4096 and 2^20, and
        the maximum is compared with the bound k*log2(W) = 12.  (The rejected
        predecessor priced k*(1-1/W) = 3.5 states here; that was wrong, and this
        section is why the new claim cannot repeat it.)
    B2  the closed-form detection law versus a direct simulation of the comparison
        decisions, exhaustively over every (mu, lam) with mu+lam <= 60.
    B3  the complete algorithm on a toy cyclic group of prime order q: success
        fraction inside the declared hard cap ceil(1.16*sqrt(q)), degenerate
        relations, bookkeeping failures (must be zero), live-bound violations
        (must be zero), and the observed working-set size.

Section C
    The cost ledger in exact integer / rational arithmetic: the hard cap
    B = ceil((116/100)*sqrt(r)) computed by integer square root, the per-step
    instruction ledger in units of the track's reference dense Fp12
    multiplication (C = 20,000), the additive items, the peak-memory bound and
    the exact power-of-ten bound on the claim.

Toy statistics are evidence at toy scope only; the production claim is analytic,
and its success probability is a computed quantile of the random-map model, not a
work bound.
"""

import json
import math
import sys
import time
from fractions import Fraction
from bisect import bisect_left
from math import isqrt

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



# --------------------------------------------------------- Section A: production

def section_a(checks=2, seed=1, table_checks=2):
    t0 = time.time()
    p, r = P_PROD, R_PROD
    out = {"p_bits": p.bit_length(), "r_bits": r.bit_length()}
    out["A1_p_prime"] = is_prime(p)
    out["A1_r_prime"] = is_prime(r)
    out["A1_r_divides_phi12"] = ((p ** 4 - p ** 2 + 1) % r == 0)
    pm = [pow(p, k, r) for k in range(13)]
    out["A1_ord_r_p_eq_12"] = (pm[12] == 1 and pm[6] == r - 1 and pm[4] != 1
                              and pm[3] != 1 and pm[2] != 1)
    out["A1_pow6_is_minus1"] = (pm[6] == r - 1)

    xi = (1 % p, 1 % p)
    q2 = p * p
    sq = fp2_pow(xi, (q2 - 1) // 2, p)
    cb = fp2_pow(xi, (q2 - 1) // 3, p)
    out["A2_tower_is_field"] = (sq != (1 % p, 0)) and (cb != (1 % p, 0))
    gamma1 = tower_constant_gamma1(p)
    out["A2_gamma1_in_Fp2"] = True

    rng = (seed * 2654435761 + 7) % (2 ** 64)
    # A3: a random GT element is nonidentity of order r
    y, rng = gt_random(p, r, rng)
    out["A3_gt_random_order_r"] = (y != fp12_one(p)) and (fp12_pow(y, r, p) == fp12_one(p))

    # A4: the collision identity on real GT elements
    ok = 0
    n = max(1, checks)
    for _ in range(n):
        g, rng = gt_random(p, r, rng)
        z, rng = gt_random(p, r, rng)
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        x = rng % (r - 1) + 1
        z = fp12_pow(g, x, p)                     # z = g^x, so log_g z = x
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        a1 = rng % r
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        b1 = rng % r
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        delta = rng % (r - 1) + 1
        a2 = (a1 - x * delta) % r
        b2 = (b1 + delta) % r
        X1 = fp12_mul(fp12_pow(g, a1, p), fp12_pow(z, b1, p), p)
        X2 = fp12_mul(fp12_pow(g, a2, p), fp12_pow(z, b2, p), p)
        den = (b2 - b1) % r
        rec = None
        if fp12_eq(X1, X2) and den != 0:
            rec = (a1 - a2) * pow(den, r - 2, r) % r
        if fp12_eq(X1, X2) and rec == x:
            ok += 1
    out["A4_collision_identity_ok"] = [ok, n]

    # A5: table entries g^(c_h)
    ok_t = 0
    nt = max(1, table_checks)
    for _ in range(nt):
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        c = rng % (r - 1) + 1
        Th = fp12_pow(y, c, p)
        if (not fp12_eq(Th, fp12_one(p))) and fp12_pow(Th, r, p) == fp12_one(p):
            ok_t += 1
    out["A5_table_in_GT_ok"] = [ok_t, nt]
    out["seconds"] = round(time.time() - t0, 2)
    return out


def _sophie_germain(bits=21):
    """prime q with p = 2q+1 prime, q about 2^(bits-1)."""
    q = (1 << (bits - 1)) + 1
    if q % 2 == 0:
        q += 1
    while True:
        if is_prime(q) and is_prime(2 * q + 1):
            return q, 2 * q + 1
        q += 2


MIX = 0x9E3779B97F4A7C15

# ---------------------------------------------------------------- the schedule
#
# Grid P_k = { ceil(2^j (k+i)/k) : j >= 0, 0 <= i < k }: exactly k points in every
# dyadic level [2^j, 2^{j+1}).  A checkpoint created at p is compared with the hare
# at every position in the time window [p+1, W*p], then dropped.
#
K_SCHED = 4
W_SCHED = 8


def sched_grid(k, hi):
    """Sorted distinct grid points of P_k that are <= hi."""
    P = set()
    j = 0
    while (1 << j) <= hi:
        p2 = 1 << j
        for i in range(k):
            v = -(-(p2 * (k + i)) // k)
            if v <= hi:
                P.add(v)
        j += 1
    return sorted(P)


def sched_live(h, P, W):
    """Checkpoints alive at hare position h: those p in P with p < h <= W*p."""
    lo = -(-h // W)
    return bisect_left(P, h) - bisect_left(P, lo)


def check_actives(k=K_SCHED, W=W_SCHED, probe=(1024, 4096, 1 << 20),
                  hi=1 << 20, n_max=None):
    """Enumerate the live-checkpoint count: the bound is asserted, not assumed.

    The comparison window has multiplicative span W = 2^b, and the grid carries k
    points per dyadic level, so at most b*k checkpoints are simultaneously live;
    this walks every hare position up to `hi` plus the probe points.
    """
    P = sched_grid(k, hi)
    probes = {}
    mx = 0
    arg = 0
    for h in range(1, hi + 1):
        n = sched_live(h, P, W)
        if n > mx:
            mx, arg = n, h
        if h in probe:
            probes[h] = n
    for h in probe:
        if h > hi:
            probes[h] = sched_live(h, sched_grid(k, h), W)
    if n_max is None:
        n_max = k * (W.bit_length() - 1)
    return {"k": k, "W": W, "grid_points_le_2^20": len(P),
            "probe_counts": {str(h): probes[h] for h in sorted(probes)},
            "max_over_1..2^20": mx, "argmax": arg,
            "n_max_formula": n_max,
            "assert_ok": mx <= n_max and max(probes.values()) <= n_max,
            "old_claim_k_times_1_minus_1_over_W": round(k * (1 - 1.0 / W), 2)}


def sched_law(mu, lam, P, W):
    """Detection index, verified exhaustively in section B2.

    m = p* + lam with p* = min{ p in P : p >= max(mu, ceil(lam/(W-1))) }.
    """
    lam = max(1, lam)
    target = max(mu, -(-lam // (W - 1)))
    return P[bisect_left(P, target)] + lam


def sched_law_simulate(mu, lam, k, W, horizon):
    """Direct simulation: first hare position that matches a live checkpoint."""
    P = [p for p in sched_grid(k, horizon) if p < horizon]
    for h in range(1, horizon + 1):
        lo = -(-h // W)
        for p in P:
            if p >= h:
                break
            if p < lo:
                continue
            if p >= mu and (h - p) % lam == 0:
                return h
    return None


def check_law(k=K_SCHED, W=W_SCHED, max_rho=60):
    """B2: the closed-form law against a direct simulation, exhaustively."""
    horizon = k * (1 << 8)
    P = sched_grid(k, horizon * 2)
    bad = []
    n = 0
    for mu in range(0, max_rho + 1):
        for lam in range(1, max_rho + 1 - mu):
            if mu + lam > max_rho:
                continue
            n += 1
            got = sched_law(mu, lam, P, W)
            want = sched_law_simulate(mu, lam, k, W, horizon)
            if want is None or got != want:
                bad.append([mu, lam, got, want])
    return {"pairs": n, "mismatches": len(bad), "first_bad": bad[:4],
            "max_ratio_m_over_rho": round(max(
                [sched_law(mu, lam, P, W) / (mu + lam)
                 for mu in range(0, max_rho + 1)
                 for lam in range(1, max_rho + 1 - mu) if mu + lam > 0]), 4)}


def _hash_index(X, R, bits=64):
    v = (X * MIX) % (1 << 64)
    v ^= v >> 29
    v = (v * 0xBF58476D1CE4E5B9) % (1 << 64)
    return (v >> (64 - (R - 1).bit_length())) & (R - 1)


def stagger_walk(q, p, g, T, R, cap, rng, k=K_SCHED, W=W_SCHED, word_mask=(1 << 16) - 1):
    """One complete staggered-checkpoint rho search on the cyclic group <g>.

    States are X = g^a * T^b with exponent pair (a, b); the walk is
    X <- X * T_h, T_h = g^(A_h) T^(B_h), h a hash of X.  A checkpoint stores
    (X & word_mask, a, b) at every grid position and is compared with the hare on
    the window [p+1, W*p] only; the stored count is asserted against k*log2(W).
    """
    A, B, Th = [], [], []
    for _ in range(R):
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        ah = rng % q
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        bh = rng % q
        A.append(ah)
        B.append(bh)
        Th.append(pow(g, ah, p) * pow(T, bh, p) % p)

    rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
    a = rng % q
    rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
    b = rng % q
    X = pow(g, a, p) * pow(T, b, p) % p

    grid = sched_grid(k, cap + 1)
    gi = 0
    live = []                      # list of (p, word, a, b, X) -- O(1) length
    max_live = 0
    compares = 0
    n_max = k * (W.bit_length() - 1)
    res = {"evals": 0, "capped": False, "solved": False, "x": None,
           "degenerate": 0, "bookfail": 0, "live_violation": 0,
           "solve_correct": False,
           "max_live": 0, "compares": compares}
    h = 0
    while h < cap:
        h += 1
        hh = _hash_index(X, R)
        X = X * Th[hh] % p
        a = (a + A[hh]) % q
        b = (b + B[hh]) % q
        # drop checkpoints whose window has closed (expiry is a single test on
        # the oldest entry: p is increasing in the list)
        while live and W * live[0][0] < h:
            live.pop(0)
        w = X & word_mask
        hit = None
        for rec in live:
            compares += 1
            if rec[1] == w and rec[4] == X:
                hit = rec
                break
        if len(live) > max_live:
            max_live = len(live)
        if max_live > n_max:
            res["live_violation"] += 1
        res["evals"] = h
        if hit is not None:
            res["solved"] = True
            ta, tb = hit[2], hit[3]
            if pow(g, ta, p) * pow(T, tb, p) % p != hit[4]:
                res["bookfail"] += 1
            den = (b - tb) % q
            if den == 0:
                res["degenerate"] += 1
                res["solved"] = False
            else:
                res["x"] = (ta - a) * pow(den, q - 2, q) % q
                res["solve_correct"] = (pow(g, res["x"], p) == T)
            res["max_live"] = max_live
            res["compares"] = compares
            return res
        # create the checkpoint at this grid position, after the comparison
        if gi < len(grid) and grid[gi] == h:
            live.append((h, w, a, b, X))
            gi += 1
    res["capped"] = True
    res["max_live"] = max_live
    res["compares"] = compares
    return res


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    ph = k / n
    d = 1 + z * z / n
    c = ph + z * z / (2 * n)
    r = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    return [round((c - r) / d, 4), round((c + r) / d, 4)]


def section_b(bits=21, trials=200, R=1024, seed=1, cap_factor=1.16,
              wall_seconds=10.0, k=K_SCHED, W=W_SCHED):
    """Toy-scale run of the complete algorithm plus the schedule assertions."""
    t0 = time.time()
    q, p = _sophie_germain(bits)
    # a generator of the order-q subgroup of F_p^* : a non-1 square, q prime
    grng = (seed * 2654435761 + 7) % (2 ** 64)
    while True:
        grng = (grng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        hh = grng % p
        gg = hh * hh % p
        if gg != 1 and pow(gg, q, p) == 1:
            break
    sqrtq = math.sqrt(q)
    cap = int(math.ceil(cap_factor * sqrtq))
    out = {"B1_actives": check_actives(k, W),
           "B2_law": check_law(k, W),
           "B3_toy": {"k": k, "W": W, "toy_q_bits": q.bit_length(), "toy_q": q,
                      "R": R, "sqrt_q": round(sqrtq, 3),
                      "cap_evals": cap, "cap_factor": cap_factor,
                      "trials_requested": trials,
                      "n_max_asserted": k * (W.bit_length() - 1)}}
    rng = (seed * 2654435761 + 7) % (2 ** 64)
    solved = n_run = degen = book = capped = livebad = solvok = 0
    evals_all, evals_ok, maxlive = [], [], 0
    for _ in range(int(trials)):
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        x = rng % (q - 1) + 1
        T = pow(gg, x, p)
        st = stagger_walk(q, p, gg, T, R, cap, rng, k, W)
        rng = (rng * 7439281 + 1) % (2 ** 64)
        n_run += 1
        evals_all.append(st["evals"])
        degen += st["degenerate"]
        book += st["bookfail"]
        livebad += st["live_violation"]
        maxlive = max(maxlive, st["max_live"])
        capped += 1 if st["capped"] else 0
        if st["solved"]:
            solved += 1
            solvok += 1 if st.get("solve_correct") else 0
            evals_ok.append(st["evals"])
        if time.time() - t0 > wall_seconds:
            break
    b3 = out["B3_toy"]
    b3["trials_run"] = n_run
    b3["trials_solved"] = solved
    b3["success_rate"] = round(solved / n_run, 4) if n_run else None
    b3["success_wilson95"] = wilson(solved, n_run)
    b3["capped_trials"] = capped
    b3["solved_plus_capped_eq_run"] = (solved + capped == n_run)
    b3["solve_correct_total"] = solvok
    b3["degenerate_relations"] = degen
    b3["book_fail_total"] = book
    b3["live_bound_violations"] = livebad
    b3["observed_max_live"] = maxlive
    if evals_all:
        srt = sorted(evals_all)
        b3["evals_over_sqrtq_median"] = round(srt[len(srt) // 2] / sqrtq, 4)
        b3["evals_over_sqrtq_p90"] = round(srt[int(0.9 * (len(srt) - 1))] / sqrtq, 4)
        b3["evals_over_sqrtq_max"] = round(srt[-1] / sqrtq, 4)
    if evals_ok:
        b3["evals_solved_over_sqrtq_mean"] = round(
            sum(evals_ok) / len(evals_ok) / sqrtq, 4)
    out["seconds"] = round(time.time() - t0, 2)
    return out


# ------------------------------------------------------------- Section C: ledger

def _ceil_sqrt_frac(num, den):
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


def section_c(cap_num=116, cap_den=100, step_instr=328,
              table_size=1024, pairing_units=1 << 17, verify_units=1 << 12,
              k=K_SCHED, W=W_SCHED):
    """Exact ledger.  Units: one dense Fp12 multiplication at the fixed C = 20,000."""
    unit = INSTR_PER_DENSE_FP12_MUL
    r = R_PROD
    out = {"r": str(r), "r_log2": round(_log2_frac(Fraction(r)), 6),
           "sqrt_r_log2": round(_log2_frac(Fraction(r)) / 2, 6),
           "reference_price": unit, "table_size": table_size,
           "cap_factor": f"{cap_num}/{cap_den}",
           "success_claim_lower": "39/100"}
    # hard cap B = ceil(c*sqrt(r)): least integer with (cap_den*B)^2 >= cap_num^2*r
    num, den = cap_num * cap_num * r, cap_den * cap_den
    B = isqrt(num // den)
    while den * B * B < num:
        B += 1
    out["cap_evals"] = str(B)
    out["cap_log2"] = round(_log2_frac(Fraction(B)), 6)
    out["cap_tight"] = (den * (B - 1) * (B - 1) < num)   # B-1 would fail
    step_units = Fraction(unit + step_instr, unit)
    out["step_instr"] = step_instr
    out["step_units"] = str(step_units)
    out["step_units_log2"] = round(_log2_frac(step_units), 6)
    main = Fraction(B) * step_units
    out["main_units_log2"] = round(_log2_frac(main), 6)
    table = table_size * 2 * 383            # g^(A_h) and z^(B_h) per entry
    reduction = pairing_units               # one pairing for g = e(P,Q0)
    verify = pairing_units + verify_units   # e(P,Q) = z and [a]Q0
    # expected cost of 64-bit prefix false positives, in reference units:
    # 2^127.7 steps * 12 live states * 2^-64 = 2^67.3 false hits, each resolved by
    # a 120-instruction full compare = 2^74.2 instructions = 2^59.9 units
    prefix_fp = 1 << 60
    add = table + reduction + verify + prefix_fp
    out["table_units"] = str(table)
    out["reduction_units"] = str(reduction)
    out["verify_units"] = str(verify)
    out["prefix_false_positive_units_log2"] = round(
        _log2_frac(Fraction(prefix_fp)), 2)
    out["additive_units"] = str(add)
    out["additive_units_log2"] = round(_log2_frac(Fraction(add)), 4)
    tot = main + add
    out["total_log2_float"] = round(_log2_frac(tot), 6)
    n_states = k * (W.bit_length() - 1)
    memory = (table_size * 576 + n_states * 112 + (1 << 16))
    out["states"] = n_states
    out["state_bytes"] = 112
    out["memory_bytes"] = str(memory)
    out["memory_log2"] = round(_log2_frac(Fraction(memory)), 3)
    m = int(_log2_frac(tot) * 100) - 3
    while True:
        if tot.numerator ** 100 <= (tot.denominator ** 100) << m:
            break
        m += 1
    out["claim_log2_exact_le"] = f"{m}/100"
    out["claim_log2"] = round(m / 100.0, 4)
    return out


def main():
    try:
        raw = sys.stdin.read()
        cfg = json.loads(raw) if raw.strip() else {}
    except Exception as e:                        # noqa: BLE001
        print(json.dumps({"error": f"bad input: {e}"}))
        return 1
    cfg = cfg or {}
    result = {"version": 2, "candidate": "gt-stagger-rho-v1",
              "schedule": {"k": K_SCHED, "W": W_SCHED,
                           "states": K_SCHED * (W_SCHED.bit_length() - 1)}}
    seed = int(cfg.get("seed", 1))
    result["A"] = section_a(checks=int(cfg.get("production_checks", 2)), seed=seed,
                            table_checks=int(cfg.get("table_checks", 2)))
    result["B"] = section_b(bits=int(cfg.get("toy_bits", 21)),
                            trials=int(cfg.get("trials", 200)),
                            R=int(cfg.get("R", 1024)), seed=seed,
                            cap_factor=float(cfg.get("cap_factor", 1.16)),
                            wall_seconds=float(cfg.get("seconds", 20.0)))
    result["C"] = section_c()
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
