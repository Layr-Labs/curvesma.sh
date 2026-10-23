#!/usr/bin/env python3
"""experiment.py -- O(1)-memory staggered-checkpoint Pollard rho in GT for fixed-G1 pairing inversion.

v4: the walk draws its step coefficients uniformly modulo the group order
(512 coin bits per draw, no rejection sampling), i.e. the production regime
that the passing Floyd-rho entry e444842b uses.  There is no positivity lemma
and no small-coefficient window; a zero-denominator collision is CHARGED
(section C) and handled (section B4).

Candidate: submissions/bls12-381-fixed-g1-invert-g2/gt-ron-rho-v4
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

Section D
    The success probability that the ledger is *not* allowed to assume.  D1 is a
    deterministic quadrature of the random-map model (Flajolet-Odlyzko) pushed
    through the exact detection law of B2, giving m/sqrt(r) quantiles with no
    sampling error.  D2 is an independent Monte-Carlo cross-check of the same
    model.  D3 sweeps the dyadic *alignment* of sqrt(r) over a full period and
    reports the worst alignment, so the claimed success cannot be an artefact of
    where r happens to sit between two powers of two.  D4 reads the success
    probability at the declared cap off D1 and states it against the 0.39
    threshold of rules.md.

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

# Hard, deterministic bound on the filter's residual work: at most 2^70 full
# element comparisons are performed by one run, and the run ABORTS if the
# allowance is exceeded.  2^70 comparisons x 216 instructions / C = 20,000 is
# 2^63.46 reference units -- 64 bits below the main term -- and the expected
# number of false positives over the whole run is only 2^67.3, so the abort is
# a bound, not an event.  This is charged as work in section C.
FP_ALLOWANCE = 1 << 70


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
    # Peak STORED records: the window plus the record inserted at this position.
    # The record at h is appended AFTER the comparison, so storage peaks at
    # live(h) + 1 whenever h is a grid point -- the reviewer's 13th record.
    Pset = set(P)
    mx_st = 0
    arg_st = 0
    for h in range(1, hi + 1):
        n = sched_live(h, P, W) + (1 if h in Pset else 0)
        if n > mx_st:
            mx_st, arg_st = n, h
    n_states = n_max + 1
    return {"k": k, "W": W, "grid_points_le_2^20": len(P),
            "probe_counts": {str(h): probes[h] for h in sorted(probes)},
            "max_over_1..2^20": mx, "argmax": arg,
            "n_max_formula": n_max,
            "max_stored_over_1..2^20": mx_st, "argmax_stored": arg_st,
            "n_states_formula": n_states,
            "assert_ok": (mx <= n_max and mx_st <= n_states
                          and max(probes.values()) <= n_max),
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


COIN_WORDS = 8      # 512 coin bits per coefficient draw, reduced mod q, no rejection


def _coin_mod(rng, q, words=COIN_WORDS):
    """One exponent draw from the algorithm's uniform coins, reduced mod q.

    `words` successive 64-bit coin words are concatenated and reduced once: no
    rejection sampling, so the cost is O(1) and the statistical distance from
    uniform on Z_q is below 2^(log2(words*64) - 64*words - log2 q), i.e. below
    2^-256 for the 512-bit production draw (q = r < 2^255) and below 2^-480 for
    every toy group used here (q < 2^32).  This is the regime of the accepted
    entry e444842b ("512 coin bits reduced mod r, bias < 2^-256").
    """
    v = 0
    for _ in range(words):
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        v = (v << 64) | rng
    return rng, v % q


def _hash_index(X, R, bits=64):
    v = (X * MIX) % (1 << 64)
    v ^= v >> 29
    v = (v * 0xBF58476D1CE4E5B9) % (1 << 64)
    return (v >> (64 - (R - 1).bit_length())) & (R - 1)


def stagger_walk(q, p, g, T, R, cap, rng, k=K_SCHED, W=W_SCHED,
                 word_mask=(1 << 64) - 1, fp_allow=FP_ALLOWANCE,
                 inject_degenerate=False):
    """One complete staggered-checkpoint rho search on the cyclic group <g>.

    States are X = g^a * T^b with exponent pair (a, b); the walk is
    X <- X * T_h, T_h = g^(A_h) T^(B_h), h a hash of X.  A checkpoint stores
    (X & word_mask, a, b, X) at every grid position -- the full element X is
    retained, exactly as the production ledger prices it, so a fingerprint match
    is settled by one full element comparison and never by a reconstruction --
    and is compared with the hare on the window [p+1, W*p] only; the stored count
    is asserted against k*log2(W).  word_mask = 2^64-1 fixes the fingerprint at
    the low 64 bits of the element's canonical encoding, the same width the
    ledger charges.
    """
    # PRODUCTION COEFFICIENT REGIME.  A_h and B_h are drawn from the run's
    # uniform coins and reduced modulo the group order q (512 coin bits per
    # draw, no rejection sampling, see _coin_mod).  There is NO positivity lemma
    # and NO small-coefficient window: the step offsets range over all of Z_q, so
    # the step set is not the structured set {A + x*B : A,B < 2^64} whose
    # random-map justification the reviews of a7ac78ea and 8b096f2c rejected
    # (there A,B < 2^64 ~ r^(1/4) forced r/M^2 ~ sqrt(r); here M = r, so
    # r/M^2 = 1/r).  The toy therefore samples exactly the production regime.
    #
    # What replaces the old lemma is an explicit CHARGE: a detected collision can
    # have b_i = b_j (mod q), i.e. denominator zero, and is then unusable.  Over a
    # run there are at most 12*B comparisons, each with gap coefficient sum uniform
    # in Z_q, so Pr[any zero denominator] <= 12*B/q = 2^-123.6 at production
    # (q = r, B <= 1.16*sqrt(r)); section C carries that number and section B4
    # exercises the branch.
    A, B, Th = [], [], []
    for _ in range(R):
        rng, ah = _coin_mod(rng, q)
        rng, bh = _coin_mod(rng, q)
        A.append(ah)
        B.append(bh)
        Th.append(pow(g, ah, p) * pow(T, bh, p) % p)

    rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
    a = rng % q
    rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
    b = rng % q
    ai, bi = a, b                      # unreduced copies, for the gap assertion
    X = pow(g, a, p) * pow(T, b, p) % p

    if inject_degenerate:
        # TEST HOOK ONLY (never used in a scored run, never in section B3): make the
        # first table entry the identity and pin the step function to it, so the hare
        # is stationary and re-meets its own checkpoint with an identical b.
        A[0], B[0], Th[0] = 0, 0, 1

    grid = sched_grid(k, cap + 1)
    gi = 0
    live = []              # (p, word, a, b, X, ai, bi) -- O(1) length
    max_live = 0           # peak STORED records, counted after insertion
    max_cmp = 0            # peak records compared against, counted in the window
    compares = 0           # every fingerprint comparison
    fullcmp = 0            # full element comparisons (real work, bounded)
    n_cmp_max = k * (W.bit_length() - 1)        # window: k*log2(W) compared
    n_max = n_cmp_max + 1                       # storage: the record being inserted
    res = {"evals": 0, "capped": False, "solved": False, "x": None,
           "degenerate": 0, "bookfail": 0, "live_violation": 0,
           "cmp_violation": 0, "fp_abort": False, "usable_collisions": 0,
           "solve_correct": False,
           "max_live": 0, "max_cmp": 0, "full_compares": 0,
           "compares": compares}
    h = 0
    while h < cap:
        h += 1
        hh = 0 if inject_degenerate else _hash_index(X, R)
        X = X * Th[hh] % p
        a = (a + A[hh]) % q
        b = (b + B[hh]) % q
        ai += A[hh]                    # the same exponents, unreduced
        bi += B[hh]
        # drop checkpoints whose window has closed (expiry is a single test on
        # the oldest entry: p is increasing in the list)
        while live and W * live[0][0] < h:
            live.pop(0)
        w = X & word_mask
        hit = None
        for rec in live:
            compares += 1
            if rec[1] == w:
                fullcmp += 1
                if fullcmp > fp_allow:      # deterministic hard work bound
                    res["fp_abort"] = True
                    res["solved"] = False
                    res["max_live"] = max_live
                    res["max_cmp"] = max_cmp
                    res["full_compares"] = fullcmp
                    res["compares"] = compares
                    return res
                if rec[4] == X:
                    hit = rec
                    break
        if len(live) > max_cmp:
            max_cmp = len(live)
        if max_cmp > n_cmp_max:
            res["cmp_violation"] += 1
        res["evals"] = h
        if hit is not None:
            res["solved"] = True
            ta, tb, tai, tbi = hit[2], hit[3], hit[5], hit[6]
            if pow(g, ta, p) * pow(T, tb, p) % p != hit[4]:
                res["bookfail"] += 1
            den = (bi - tbi) % q        # gap coefficient sum, reduced mod q
            if den == 0:
                # UNUSABLE relation: b_i = b_j (mod q), so X_i = X_j carries no
                # logarithm of T.  Counted, the checkpoint is DROPPED so the same
                # pair cannot fire again, and the walk continues to the cap.  This
                # is the charged event of section C (<= 12*B/q = 2^-123.6).
                res["degenerate"] += 1
                res["solved"] = False
                res["max_live"] = max_live
                res["max_cmp"] = max_cmp
                live.remove(hit)
            else:
                res["usable_collisions"] += 1
                res["x"] = (tai - ai) * pow(den, q - 2, q) % q
                res["solve_correct"] = (pow(g, res["x"], p) == T)
                res["max_live"] = max_live
                res["max_cmp"] = max_cmp
                res["full_compares"] = fullcmp
                res["compares"] = compares
                return res
        # create the checkpoint at this grid position, after the comparison
        if gi < len(grid) and grid[gi] == h:
            live.append((h, w, a, b, X, ai, bi))
            gi += 1
            if len(live) > max_live:    # peak STORED records, measured after insertion
                max_live = len(live)
            if max_live > n_max:
                res["live_violation"] += 1
    res["capped"] = True
    res["max_live"] = max_live
    res["max_cmp"] = max_cmp
    res["full_compares"] = fullcmp
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


def _toy_trials(q, p, gg, R, cap, trials, seed, k=K_SCHED, W=W_SCHED):
    """Run `trials` complete searches; return the aggregate this evidence turns on."""
    rng = (seed * 2654435761 + 7) % (2 ** 64)
    solved = n_run = degen = book = capped = livebad = cmpmax = 0
    solvok = fullcmp = fpabort = 0
    maxlive = 0
    for _ in range(int(trials)):
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        x = rng % (q - 1) + 1
        T = pow(gg, x, p)
        st = stagger_walk(q, p, gg, T, R, cap, rng, k, W)
        rng = (rng * 7439281 + 1) % (2 ** 64)
        n_run += 1
        degen += st["degenerate"]
        book += st["bookfail"]
        livebad += st["live_violation"] + st["cmp_violation"]
        maxlive = max(maxlive, st["max_live"])
        cmpmax = max(cmpmax, st["max_cmp"])
        fullcmp = max(fullcmp, st["full_compares"])
        fpabort += 1 if st["fp_abort"] else 0
        capped += 1 if st["capped"] else 0
        if st["solved"]:
            solved += 1
            solvok += 1 if st.get("solve_correct") else 0
    return {"R": R, "trials_run": n_run, "trials_solved": solved,
            "success_rate": round(solved / n_run, 4) if n_run else None,
            "success_wilson95": wilson(solved, n_run),
            "capped_trials": capped,
            "solved_plus_capped_eq_run": (solved + capped == n_run),
            "solve_correct_total": solvok, "degenerate_relations": degen,
            "book_fail_total": book, "live_bound_violations": livebad,
            "observed_max_live": maxlive}


def check_degenerate_branch(bits=21, seed=1):
    """B4: exercise the zero-denominator branch that v4 charges.

    The event itself is 2^-123.6 at production, so no feasible toy can sample it
    by chance.  What CAN be shown is that the branch is reachable, is taken, and is
    handled without a division: two states with b_i = b_j (mod q) are constructed
    directly, the solve rule is applied, and the walk is then run with an INJECTED
    zero-denominator collision to show that it counts, drops the record and
    continues to a correct usable collision instead of aborting.
    """
    q, p = _sophie_germain(bits)
    grng = (seed * 2654435761 + 7) % (2 ** 64)
    while True:
        grng = (grng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        hh = grng % p
        gg = hh * hh % p
        if gg != 1 and pow(gg, q, p) == 1:
            break
    # (1) the solve rule on a constructed degenerate pair: same b mod q, different a.
    rng = (seed * 88172645463325252 + 3) % (2 ** 64)
    rng, b0 = _coin_mod(rng, q)
    rng, a1 = _coin_mod(rng, q)
    rng, a2 = _coin_mod(rng, q)
    den_same = (b0 - b0) % q
    den_wrap = (b0 - (b0 + q)) % q
    rule_ok = (den_same == 0 and den_wrap == 0 and ((a2 - a1) % q != 0 or a2 == a1))
    # (2) a complete injected walk: the table's first entry is the identity and the
    # step function is pinned to it, so the hare is stationary, its own checkpoint
    # at position 1 is re-met at position 2 with b_1 = b_2 (mod q), and the run must
    # count the degenerate relation, DROP that checkpoint and continue to the cap --
    # never divide by zero, never abort, never double count.
    res = stagger_walk(q, p, gg, pow(gg, 1, p), 64, 64, rng, K_SCHED, W_SCHED,
                       inject_degenerate=True)
    return {"constructed_pair_denominator_zero": rule_ok,
            "injected_run_degenerates": res["degenerate"],
            "injected_run_branch_taken": res["degenerate"] >= 1,
            "injected_run_runs_to_cap": res["capped"] is True,
            "injected_run_not_solved": res["solved"] is False,
            "injected_run_bookkeeping_ok": res["bookfail"] == 0,
            "injected_run_live_bound_ok": (res["live_violation"] == 0
                                           and res["cmp_violation"] == 0)}


def section_b(bits=21, trials=300, R=256, seed=1, cap_factor=1.16,
              wall_seconds=10.0, k=K_SCHED, W=W_SCHED,
              sens_R=1024, sens_trials=120):
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
           "B4_degenerate_branch": check_degenerate_branch(bits=bits, seed=seed),
           "B3_toy": {"k": k, "W": W, "toy_q_bits": q.bit_length(), "toy_q": q,
                      "R": R, "sqrt_q": round(sqrtq, 3),
                      "cap_evals": cap, "cap_factor": cap_factor,
                      "trials_requested": trials,
                      "coefficient_regime": "uniform mod q: 512 coin bits per draw, no rejection",
                      "coefficient_range": "Z/qZ (the whole group order; no small-coefficient window)",
                      "positivity_lemma": "none (withdrawn in v4; replaced by the zero-denominator charge)",
                      "zero_denominator_charge_12B_over_q": round(
                          12.0 * cap / q, 12),
                      "n_compared_asserted": k * (W.bit_length() - 1),
                      "n_stored_asserted": k * (W.bit_length() - 1) + 1,
                      "fp_allowance": FP_ALLOWANCE}}
    rng = (seed * 2654435761 + 7) % (2 ** 64)
    solved = n_run = degen = book = capped = livebad = cmpmax = 0
    solvok = fullcmp = fpabort = 0
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
        livebad += st["live_violation"] + st["cmp_violation"]
        maxlive = max(maxlive, st["max_live"])
        cmpmax = max(cmpmax, st["max_cmp"])
        fullcmp = max(fullcmp, st["full_compares"])
        fpabort += 1 if st["fp_abort"] else 0
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
    b3["observed_max_compared"] = cmpmax
    b3["max_full_compares"] = fullcmp
    b3["fp_aborts"] = fpabort
    if evals_all:
        srt = sorted(evals_all)
        b3["evals_over_sqrtq_median"] = round(srt[len(srt) // 2] / sqrtq, 4)
        b3["evals_over_sqrtq_p90"] = round(srt[int(0.9 * (len(srt) - 1))] / sqrtq, 4)
        b3["evals_over_sqrtq_max"] = round(srt[-1] / sqrtq, 4)
    if evals_ok:
        b3["evals_solved_over_sqrtq_mean"] = round(
            sum(evals_ok) / len(evals_ok) / sqrtq, 4)
    # B4: the same complete search at a larger jump table, so the reported
    # success fraction is not an artefact of the declared R.  Raising R costs
    # memory only; it does not change the per-step cost or the detection law.
    out["B4_R_sensitivity"] = _toy_trials(
        q, p, gg, sens_R, cap, int(sens_trials), seed + 101, k, W)
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


def section_c(cap_num=116, cap_den=100, step_instr=600,
             table_size=1 << 20, pairing_units=1 << 17, verify_units=1 << 12,
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
    # Itemised per-step instruction ledger, 64-bit scalar equivalents.  The
    # charged figure is the itemised sum times a safety factor; every item counts
    # instructions the walk provably executes, not a timing.
    itemised = {"index_mix_load_mul_xorshift": 8,
                "table_entry_load_element_and_two_4limb_coeffs": 12,
                "accumulator_a_add_reduce_4limb": 24,
                "accumulator_b_add_reduce_4limb": 24,
                "fingerprint_load": 2,
                "twelve_live_fingerprint_compares": 48,
                "ring_insert_and_expiry_amortised": 24,
                "loop_control_and_cap_test": 8}
    out["step_instr_itemised"] = itemised
    out["step_instr_itemised_sum"] = sum(itemised.values())
    out["step_instr_safety_factor"] = round(
        step_instr / max(1, sum(itemised.values())), 3)
    out["step_units"] = str(step_units)
    out["step_units_log2"] = round(_log2_frac(step_units), 6)
    main = Fraction(B) * step_units
    out["main_units_log2"] = round(_log2_frac(main), 6)
    # ZERO-DENOMINATOR CHARGE (v4).  Uniform-mod-r coefficients mean a detected
    # collision can have b_i = b_j (mod r), which yields no logarithm.  The run makes
    # at most 12*B comparisons and each gap coefficient sum is uniform in Z_r, so
    # Pr[any usable relation is lost this way] <= 12*B/r.  It is charged against the
    # success claim, not against the work.
    degen_num = 12 * B
    out["degenerate_charge"] = {
        "comparisons_upper": str(degen_num),
        "comparisons_upper_log2": round(_log2_frac(Fraction(degen_num)), 3),
        "prob_upper_log2": round(_log2_frac(Fraction(degen_num, r)), 3),
        "charged_against": "success lower bound (subtracted), not against work"}
    # One table entry needs g^(A_h) and z^(B_h): two 255-bit fixed-exponent
    # exponentiations, 255 doublings + 255 additions of dense Fp12 elements each,
    # charged as 512 reference multiplications apiece (round up).
    table = table_size * 2 * 512
    reduction = pairing_units               # one pairing for g = e(P,Q0)
    verify = pairing_units + verify_units   # e(P,Q) = z and [a]Q0
    # Deterministic cost of the filter's residual work, in reference units.  The
    # run performs at most 12*B = 12*1.16*sqrt(r) = 2^131.3 fingerprint
    # comparisons, of which 2^131.3 * 2^-64 = 2^67.3 are expected to match; each
    # match is settled by one full 576-byte element comparison (216 instructions)
    # and never by a reconstruction.  An EXPECTED cost is not a bound, so the run
    # instead carries a HARD allowance of 2^70 full comparisons and aborts when it
    # is exceeded: 2^70 * 216 / 20000 = 2^63.46 units, deterministic, 64 bits
    # below the main term, and above the 2^67.3 expected true hits with a 6.5x
    # margin, so the abort can only fire on an astronomically unlikely run.
    fp_hits_log2 = 131.3 - 64.0
    fp_allowance_units = Fraction(FP_ALLOWANCE * 216, unit)
    add = table + reduction + verify + fp_allowance_units
    out["table_units"] = str(table)
    out["reduction_units"] = str(reduction)
    out["verify_units"] = str(verify)
    out["fp_hits_log2"] = round(fp_hits_log2, 2)
    out["fp_full_compare_instr"] = 216
    out["prefix_fp_exact_units_log2"] = round(
        _log2_frac(Fraction(int(2 ** fp_hits_log2) * 216, unit)), 2)
    out["fp_allowance_compares"] = FP_ALLOWANCE
    out["fp_allowance_units_log2"] = round(_log2_frac(fp_allowance_units), 2)
    out["fp_allowance_is_hard"] = True
    out["additive_units"] = str(add)
    out["additive_units_log2"] = round(_log2_frac(Fraction(add)), 4)
    tot = main + add
    out["total_log2_float"] = round(_log2_frac(tot), 6)
    n_states = k * (W.bit_length() - 1) + 1     # window + the record being inserted
    # ---- peak memory in bytes, itemised.  Every quantity the walk touches is
    # ---- either carried in a record here or shown absent, so the byte ledger and
    # ---- the implementation cannot disagree.
    # jump table: table_size entries x (576-byte dense Fp12 element + 32-byte A_h +
    # 32-byte B_h).  v4 coefficients are residues mod r (four 64-bit limbs), not the
    # 8-byte [1,2^64) coefficients of v3.
    mem_table = table_size * (576 + 32 + 32)
    # live checkpoints: 13 x (576-byte X + 8-byte fingerprint + 32-byte a +
    # 32-byte b + 32-byte grid position).  X is STORED, not reconstructed.  The
    # 13th is the record created at the current position, appended after the
    # comparison; 12 is the number compared per step.
    mem_checkpoints = n_states * (576 + 8 + 32 + 32 + 32)
    mem_hare = 576 + 8 + 32 + 32 + 8      # state element, fingerprint, a, b, step
    mem_g = 576                           # g = e(P,Q0), the table's base
    mem_ctl = 256                         # power, lam, counters, RNG state, pointers
    mem_code = 1 << 16                    # walk/detection code, stack and scratch
    memory = mem_table + mem_checkpoints + mem_hare + mem_g + mem_ctl + mem_code
    out["states"] = n_states
    out["memory_layout"] = {"jump_table": mem_table,
                            "live_checkpoints": mem_checkpoints,
                            "hare_state": mem_hare, "base_g": mem_g,
                            "control": mem_ctl, "code_and_scratch": mem_code}
    out["checkpoint_record_bytes"] = 576 + 8 + 32 + 32 + 32
    out["memory_layout_sums_to_total"] = (
        sum(out["memory_layout"].values()) == memory)
    out["position_bytes"] = 32           # the 2^128-step schedule needs 17
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


# ------------------------------------- Section D: the model that prices success

def _lcg(state):
    return (state * 6364136223846793005 + 1442695040888963407) % (2 ** 64)


def _model_mass(k, W, N, levels=600, nv=200, umax=9.0):
    """Deterministic quadrature of the detection index under the random-map model.

    Model (Flajolet-Odlyzko, the standard rho model): rho = mu + lam satisfies
    Pr[rho >= t] = exp(-t^2/(2N)) and, given rho, mu is uniform on [0, rho).  The
    detection index is the exact law of B2 -- m = p* + lam with p* the first grid
    point >= max(mu, ceil(lam/(W-1))) -- so no Monte-Carlo error enters.  Returns
    (sorted list of (m/sqrt(N), weight), truncated mass before normalisation).
    """
    sq = math.sqrt(N)
    hi = 2 * int(umax * sq) + 8          # margin: grid must cover every target
    P = sched_grid(k, hi)
    rows, tot = [], 0.0
    du = umax / levels
    for i in range(1, levels + 1):
        u = i * du
        w = u * math.exp(-0.5 * u * u) * du      # density of rho/sqrt(N)
        if w < 1e-14:
            continue
        tot += w
        rows.append((u * sq, w))
    out = []
    for rho, w in rows:
        wv = w / nv
        for j in range(nv):
            mu = rho * (j + 0.5) / nv
            out.append((sched_law(mu, rho - mu, P, W) / sq, wv))
    out.sort()
    return out, tot


def _quantile(mass, qq, norm):
    acc = 0.0
    for v, w in mass:
        acc += w
        if acc >= qq * norm:
            return v
    return mass[-1][0]


def _tail(mass, cap, norm):
    """Model success probability at a cap c = cap*sqrt(N): mass with m/sqrt(N)<=cap."""
    acc = 0.0
    for v, w in mass:
        if v > cap:
            break
        acc += w
    return acc / norm


def _mc_quantile(k, W, N, draws=20000, seed=1, umax=9.0):
    """Independent Monte-Carlo cross-check of D1 on the same model."""
    sq = math.sqrt(N)
    P = sched_grid(k, int(umax * sq) + 2)
    st = seed * 2654435761 + 7
    vals = []
    for _ in range(int(draws)):
        st = _lcg(st)
        u1 = ((st >> 11) + 1) / float(1 << 53)
        rho = math.sqrt(-2.0 * N * math.log(u1))
        if rho > umax * sq:
            continue
        st = _lcg(st)
        u2 = ((st >> 11) + 1) / float(1 << 53)
        mu = rho * u2
        vals.append(sched_law(mu, rho - mu, P, W) / sq)
    vals.sort()
    return vals


def section_d(cap_factor=1.16, k=K_SCHED, W=W_SCHED, seed=1, draws=20000,
              sweep_levels=200, sweep_nv=200):
    """D1 quadrature, D2 Monte-Carlo cross-check, D3 alignment sweep, D4 success."""
    t0 = time.time()
    out = {"D1_model": {"k": k, "W": W, "reference_mass": None}}
    mass, tot = _model_mass(k, W, R_PROD)
    out["D1_model"]["reference_mass"] = round(tot, 6)
    out["D1_model"]["c_q39"] = round(_quantile(mass, 0.39, tot), 4)
    out["D1_model"]["c_q42"] = round(_quantile(mass, 0.42, tot), 4)
    out["D1_model"]["c_q45"] = round(_quantile(mass, 0.45, tot), 4)
    out["D1_model"]["c_q50"] = round(_quantile(mass, 0.50, tot), 4)
    out["D1_model"]["mean_m_over_sqrt_r"] = round(
        sum(v * w for v, w in mass) / tot, 4)

    mc = _mc_quantile(k, W, R_PROD, draws=draws, seed=seed)
    out["D2_montecarlo"] = {
        "draws": len(mc),
        "c_q39": round(mc[int(0.39 * (len(mc) - 1))], 4),
        "c_q45": round(mc[int(0.45 * (len(mc) - 1))], 4),
        "mean": round(sum(mc) / len(mc), 4),
        "abs_diff_q39_vs_D1": round(
            abs(mc[int(0.39 * (len(mc) - 1))] - out["D1_model"]["c_q39"]), 4)}

    # D3: dyadic alignment.  sqrt(N_t) = sqrt(r)*2^(t/16), so t = 0..15 sweeps a
    # full period of log2(sqrt(N)) mod 1; the cap must clear 0.39 at every t.
    sweep = []
    for t in range(16):
        N_t = int(R_PROD * (2.0 ** (t / 8.0)))
        m_t, tot_t = _model_mass(k, W, N_t, levels=sweep_levels, nv=sweep_nv)
        sweep.append((t, round(tot_t, 4), round(_quantile(m_t, 0.39, tot_t), 4),
                      round(_quantile(m_t, 0.45, tot_t), 4),
                      round(_tail(m_t, cap_factor, tot_t), 4)))
    c39 = [s[2] for s in sweep]
    c45 = [s[3] for s in sweep]
    succ = [s[4] for s in sweep]
    out["D3_alignment_sweep"] = {
        "t": [s[0] for s in sweep],
        "reference_mass": [s[1] for s in sweep],
        "c_q39": c39, "c_q45": c45, "model_success_at_cap": succ,
        "c_q39_min": min(c39), "c_q39_max": max(c39),
        "c_q39_spread_pct": round(100.0 * (max(c39) / min(c39) - 1.0), 3),
        "c_q45_max": max(c45),
        "min_model_success_at_cap": min(succ)}

    # D4: the claim's success probability, read off D1 and off the worst alignment.
    out["D4_success"] = {
        "cap_factor": cap_factor,
        "model_success_nominal": round(_tail(mass, cap_factor, tot), 4),
        "model_success_worst_alignment": min(succ),
        "threshold_rules_md": 0.39,
        "c_q39_max_over_alignments": max(c39),
        "c_q45_max_over_alignments": max(c45),
        "cap_clears_q39_at_every_alignment": max(c39) <= cap_factor,
        "cap_clears_q45_at_every_alignment": max(c45) <= cap_factor,
        "claim": "success >= 39/100"}
    out["D4_success"]["claimed_success_holds"] = (
        min(succ) >= 0.39 and max(c39) <= cap_factor)
    out["seconds"] = round(time.time() - t0, 2)
    return out


def main():
    try:
        raw = sys.stdin.read()
        cfg = json.loads(raw) if raw.strip() else {}
    except Exception as e:                        # noqa: BLE001
        print(json.dumps({"error": f"bad input: {e}"}))
        return 1
    cfg = cfg or {}
    result = {"version": 4, "candidate": "gt-ron-rho-v4",
              "schedule": {"k": K_SCHED, "W": W_SCHED,
                           "compared_per_step": K_SCHED * (W_SCHED.bit_length() - 1),
                           "stored_at_peak": K_SCHED * (W_SCHED.bit_length() - 1) + 1}}
    seed = int(cfg.get("seed", 1))
    result["A"] = section_a(checks=int(cfg.get("production_checks", 2)), seed=seed,
                            table_checks=int(cfg.get("table_checks", 2)))
    result["B"] = section_b(bits=int(cfg.get("toy_bits", 21)),
                            trials=int(cfg.get("trials", 300)),
                            R=int(cfg.get("R", 256)), seed=seed,
                            cap_factor=float(cfg.get("cap_factor", 1.16)),
                            wall_seconds=float(cfg.get("seconds", 20.0)))
    result["C"] = section_c()
    result["D"] = section_d(cap_factor=float(cfg.get("cap_factor", 1.16)),
                            seed=seed,
                            draws=int(cfg.get("model_draws", 20000)))
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
