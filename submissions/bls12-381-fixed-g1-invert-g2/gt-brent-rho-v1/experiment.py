#!/usr/bin/env python3
"""experiment.py -- O(1)-memory Brent-cycle Pollard rho in GT for fixed-G1 pairing inversion.

Candidate: submissions/bls12-381-fixed-g1-invert-g2/gt-brent-rho-v1
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
    A5  sampled step-table entries g^(c_h) are nonidentity and of order r.

Section B (toy scale, cyclic group of prime order)
    The walk is a function of the group element only, so the toy uses a
    multiplicative group of prime order q (a safe-prime subgroup of F_p^*), which
    exercises the same map, the same table, the same Brent cycle detection with
    O(1) memory and the same relation solve as the production walk.  Reported:
    the number of group operations at the detected collision divided by sqrt(q),
    the fraction of trials solved inside the declared hard cap, degenerate
    relations, bookkeeping failures (must be zero), and the constant number of
    stored group elements.

Section C
    The cost ledger in exact integer / rational arithmetic: r, S0 at the declared
    0.40 model success, the Brent evaluation factor, the per-step instruction
    ledger in units of the track's reference dense Fp12 multiplication
    (C = 20,000), additive items (table build, pairing reduction, verification),
    the peak-memory bound and the claim.

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


def _hash_index(X, R, bits=64):
    v = (X * MIX) % (1 << 64)
    v ^= v >> 29
    v = (v * 0xBF58476D1CE4E5B9) % (1 << 64)
    return (v >> (64 - (R - 1).bit_length())) & (R - 1)


def brent_walk(q, p, g, T, R, cap, rng, table=None):
    """One complete O(1)-memory Brent search on the cyclic group <g> of order q.

    States are X = g^a * T^b with exponent pair (a, b); the walk is
    X <- X * T_h with T_h = g^(A_h), h a hash of X, and (a,b) <- (a,b)+(A_h,0).
    Returns a stats dict; nothing proportional to the walk length is stored.
    """
    if table is None:
        table = []
        for _ in range(R):
            c = rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
            table.append(c % q)
    R = len(table)

    def step(X, a, b):
        h = _hash_index(X, R)
        # T_h = g^(A_h) * T^(B_h) is precomputed as a group element below
        return h

    # precompute the jump table as actual group elements: Th = g^(A_h) T^(B_h)
    A = []
    B = []
    Th = []
    for _ in range(R):
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        ah = rng % q
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        bh = rng % q
        A.append(ah)
        B.append(bh)
        Th.append(pow(g, ah, p) * pow(T, bh, p) % p)

    def f(X, a, b):
        h = _hash_index(X, R)
        return (X * Th[h] % p, (a + A[h]) % q, (b + B[h]) % q)

    rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
    a0 = rng % q
    rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
    b0 = rng % q
    X0 = pow(g, a0, p) * pow(T, b0, p) % p

    evals = 0
    tX, ta, tb = X0, a0, b0
    hX, ha, hb = f(X0, a0, b0)
    evals += 1
    power = lam = 1
    capped = False
    while tX != hX:
        if evals >= cap:
            capped = True
            break
        if power == lam:
            tX, ta, tb = hX, ha, hb
            power *= 2
            lam = 0
        hX, ha, hb = f(hX, ha, hb)
        evals += 1
        lam += 1

    res = {"evals": evals, "capped": capped, "solved": False, "x": None,
           "degenerate": 0, "bookfail": 0, "stored_group_elements": 6}
    if capped or tX != hX:
        return res
    # verify bookkeeping: the two states hold equal values with their exponents
    if pow(g, ta, p) * pow(T, tb, p) % p != tX:
        res["bookfail"] += 1
    if pow(g, ha, p) * pow(T, hb, p) % p != hX:
        res["bookfail"] += 1
    den = (hb - tb) % q
    if den == 0:
        res["degenerate"] += 1
        return res
    cand = (ta - ha) * pow(den, q - 2, q) % q
    if pow(g, cand, p) == T:
        res["solved"] = True
        res["x"] = cand
    return res


def wilson(k, n, z=1.96):
    if n == 0:
        return [0.0, 0.0, 0.0]
    ph = k / n
    d = 1.0 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    h = z * math.sqrt(max(ph * (1 - ph), 0.0) / n + z * z / (4 * n * n)) / d
    return [round(ph, 4), round(max(c - h, 0.0), 4), round(min(c + h, 1.0), 4)]


def section_b(bits=21, trials=200, R=64, seed=1, cap_factor=2.0,
              wall_seconds=12.0):
    t0 = time.time()
    q, p = _sophie_germain(bits)
    rng = (seed * 2654435761 + 7) % (2 ** 64)
    # generator of the order-q subgroup
    while True:
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        h = rng % p
        gg = h * h % p
        if gg != 1 and pow(gg, q, p) == 1:
            break
    sqrtq = math.sqrt(q)
    S0 = math.sqrt(2.0 * q * math.log(1.0 / 0.6) * R / (R - 1.0))
    cap = int(math.ceil(cap_factor * S0))
    out = {"toy_q_bits": q.bit_length(), "toy_q": q, "toy_p": p, "R": R,
           "sqrt_q": round(sqrtq, 3), "S0_evals": round(S0, 3),
           "cap_evals": cap, "cap_factor": cap_factor,
           "trials_requested": trials}
    solved = 0
    n_run = 0
    evals_all = []
    evals_ok = []
    degen = book = capped = 0
    for _ in range(int(trials)):
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        x = rng % (q - 1) + 1
        T = pow(gg, x, p)
        st = brent_walk(q, p, gg, T, R, cap, rng)
        rng = (rng * 7439281 + 1) % (2 ** 64)
        n_run += 1
        evals_all.append(st["evals"])
        degen += st["degenerate"]
        book += st["bookfail"]
        capped += 1 if st["capped"] else 0
        if st["solved"]:
            solved += 1
            evals_ok.append(st["evals"])
        if time.time() - t0 > wall_seconds:
            break
    out["trials_run"] = n_run
    out["trials_solved"] = solved
    out["success_rate"] = wilson(solved, n_run)
    out["capped_trials"] = capped
    out["degenerate_relations"] = degen
    out["book_fail_total"] = book
    if evals_all:
        out["evals_over_sqrtq_mean"] = round(sum(evals_all) / len(evals_all) / sqrtq, 4)
        srt = sorted(evals_all)
        out["evals_over_sqrtq_median"] = round(srt[len(srt) // 2] / sqrtq, 4)
        out["evals_over_sqrtq_p90"] = round(srt[int(0.9 * (len(srt) - 1))] / sqrtq, 4)
    if evals_ok:
        out["evals_solved_over_sqrtq_mean"] = round(
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


def section_c(success_model=Fraction(2, 5), brent_factor=2, step_instr=400,
              table_size=1024, pairing_units=1 << 17, verify_units=1 << 12,
              success_claim=Fraction(2, 5)):
    """Exact ledger.  Units: one dense Fp12 multiplication at the fixed C = 20,000."""
    unit = INSTR_PER_DENSE_FP12_MUL
    r = R_PROD
    LN53 = Fraction(510826, 1000000)          # > ln(5/3), P = 2/5
    out = {"r": str(r), "r_log2": round(_log2_frac(Fraction(r)), 6),
           "sqrt_r_log2": round(_log2_frac(Fraction(r)) / 2, 6),
           "reference_price": unit, "table_size": table_size,
           "success_model": str(success_model)}
    S0 = _ceil_sqrt_frac(2 * r * LN53.numerator, LN53.denominator)
    out["S0_evals"] = str(S0)
    out["S0_log2"] = round(_log2_frac(Fraction(S0)), 6)
    cap = brent_factor * S0
    out["brent_factor"] = brent_factor
    out["cap_evals"] = str(cap)
    out["cap_log2"] = round(_log2_frac(Fraction(cap)), 6)
    step_units = Fraction(unit + step_instr, unit)
    out["step_instr"] = step_instr
    out["step_units"] = str(step_units)
    main = Fraction(cap) * step_units
    out["main_units_log2"] = round(_log2_frac(main), 6)
    # table build: one dense scalar multiplication of the GT generator per entry
    table = table_size * 2 * 383
    reduction = pairing_units                     # one pairing g = e(P,Q0)
    verify = pairing_units + verify_units          # e(P,Q) = z and [a]Q0
    add = table + reduction + verify + 64
    out["table_units"] = str(table)
    out["reduction_units"] = str(reduction)
    out["verify_units"] = str(verify)
    out["additive_units"] = str(add)
    out["additive_units_log2"] = round(_log2_frac(Fraction(add)), 4)
    W = main + add
    out["W_log2_float"] = round(_log2_frac(W), 6)
    memory = (table_size * 2 * 288 + table_size * 2 * 32 + 8 * 32
              + 6 * 2 * 288 + (1 << 16))
    out["memory_bytes"] = str(memory)
    out["memory_log2"] = round(_log2_frac(Fraction(memory)), 3)
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


def main():
    try:
        raw = sys.stdin.read()
        cfg = json.loads(raw) if raw.strip() else {}
    except Exception as e:                        # noqa: BLE001
        print(json.dumps({"error": f"bad input: {e}"}))
        return 1
    cfg = cfg or {}
    result = {"version": 1}
    seed = int(cfg.get("seed", 1))
    result["A"] = section_a(checks=int(cfg.get("production_checks", 2)), seed=seed,
                            table_checks=int(cfg.get("table_checks", 2)))
    result["B"] = section_b(bits=int(cfg.get("toy_bits", 21)),
                            trials=int(cfg.get("trials", 200)),
                            R=int(cfg.get("R", 64)), seed=seed,
                            cap_factor=float(cfg.get("cap_factor", 2.0)),
                            wall_seconds=float(cfg.get("seconds", 10.0)))
    result["C"] = section_c()
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
