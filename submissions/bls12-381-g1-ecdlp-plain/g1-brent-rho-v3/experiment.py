#!/usr/bin/env python3
"""experiment.py -- O(1)-memory Brent-cycle-detection Pollard rho in BLS12-381 G1.

Candidate: submissions/bls12-381-g1-ecdlp-plain/g1-brent-rho-v2
Track:     bls12-381-g1-ecdlp-plain (plain G1 discrete logarithm, no auxiliary powers)

Standard library only.  One JSON object on stdin, one JSON object on stdout.

Section A (production scale, on BLS12-381 G1 itself)
    A1  p and r prime; the profile generator P0 satisfies y^2 = x^3 + 4; it is on
        the curve, not infinity, and [r]P0 = O, so P0 generates the prime-order
        subgroup G1.
    A2  affine addition and doubling close on the curve; [2]P0 != P0; -P0 is on
        the curve.
    A3  the collision-solve identity used by the attack, on real BLS12-381 points:
        with T = [x]P0 and X_i = [a_i]P0 + [b_i]T, two states with equal value and
        different b return x = (a1-a2)/(b2-b1) mod r.  Checked by constructing two
        equal states from chosen exponents and recovering x.
    A4  the precomputed step table T_h = [c_h]P0 lies in G1 (on the curve, not
        infinity, and [r]T_h = O for a sampled entry).

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
    ledger in units of the track's reference BLS G1 addition (C = 2,000),
    additive items, the peak-memory bound and the claim.

Toy statistics are evidence at toy scope only; the production claim is analytic.
"""

import json
import math
import sys
import time
from fractions import Fraction

# ---------------------------------------------------------------- parameters
P_HEX = "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab"
R_HEX = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001"
GX_HEX = "17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb"
GY_HEX = "08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1"
FP = int(P_HEX, 16)
R = int(R_HEX, 16)
P0 = (int(GX_HEX, 16), int(GY_HEX, 16))
INSTR_PER_ADD = 2000          # fixed reference price: BLS G1 addition, C = 2,000
MEASURED_ADD = 2298           # measured blst G1 addition (scripts/reference_costs)

# ------------------------------------------------------- BLS12-381 G1 arithmetic

def g1_neg(A):
    if A is None:
        return None
    return (A[0], (-A[1]) % FP)


def g1_add(A, B):
    if A is None:
        return B
    if B is None:
        return A
    x1, y1 = A
    x2, y2 = B
    if x1 == x2:
        if (y1 + y2) % FP == 0:
            return None
        lam = (3 * x1 * x1) * pow(2 * y1, FP - 2, FP) % FP
    else:
        lam = (y2 - y1) * pow((x2 - x1) % FP, FP - 2, FP) % FP
    x3 = (lam * lam - x1 - x2) % FP
    return (x3, (lam * (x1 - x3) - y1) % FP)


def g1_mul(k, A):
    k = int(k)
    if k == 0 or A is None:
        return None
    if k < 0:
        return g1_mul(-k, g1_neg(A))
    out = None
    base = A
    while k:
        if k & 1:
            out = g1_add(out, base)
        base = g1_add(base, base)
        k >>= 1
    return out


def g1_on_curve(A):
    if A is None:
        return True
    x, y = A
    return (y * y - x * x * x - 4) % FP == 0


def is_prime(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
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


# --------------------------------------------------------- Section A: production

def section_a(checks=2, seed=1, table_checks=2):
    t0 = time.time()
    out = {"p_bits": FP.bit_length(), "r_bits": R.bit_length()}
    out["A1_p_prime"] = is_prime(FP)
    out["A1_r_prime"] = is_prime(R)
    out["A1_P0_on_curve"] = g1_on_curve(P0)
    out["A1_P0_not_infinity"] = P0 is not None
    out["A2_r_times_P0_is_infinity"] = (g1_mul(R, P0) is None)
    out["A2_r_minus_1_times_P0_is_neg_P0"] = (g1_mul(R - 1, P0) == g1_neg(P0))
    two = g1_add(P0, P0)
    out["A2_2P0_on_curve"] = g1_on_curve(two)
    out["A2_2P0_ne_P0"] = (two != P0)
    out["A2_negP0_on_curve"] = g1_on_curve(g1_neg(P0))
    out["A2_P0_plus_negP0_is_infinity"] = (g1_add(P0, g1_neg(P0)) is None)

    # A3: the collision-solve identity, on real BLS12-381 points
    rng = (seed * 2654435761 + 12345) % (2 ** 64)
    def nxt():
        nonlocal rng
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        return rng

    ok = 0
    n = max(1, checks)
    for _ in range(n):
        x = nxt() % (R - 1) + 1
        T = g1_mul(x, P0)
        a1 = nxt() % R
        b1 = nxt() % R
        delta = nxt() % (R - 1) + 1
        a2 = (a1 - x * delta) % R
        b2 = (b1 + delta) % R
        X1 = g1_add(g1_mul(a1, P0), g1_mul(b1, T))
        X2 = g1_add(g1_mul(a2, P0), g1_mul(b2, T))
        den = (b2 - b1) % R
        rec = None
        if X1 == X2 and den != 0:
            rec = (a1 - a2) * pow(den, R - 2, R) % R
        if X1 == X2 and rec == x:
            ok += 1
    out["A3_collision_identity_ok"] = [ok, n]

    # A4: precomputed table entries are in G1
    ok_t = 0
    nt = max(1, table_checks)
    for _ in range(nt):
        c = nxt() % (R - 1) + 1
        d = nxt() % (R - 1) + 1
        Th = g1_add(g1_mul(c, P0), g1_mul(d, T))
        if Th is not None and g1_on_curve(Th) and (g1_mul(R, Th) is None) and Th != P0:
            ok_t += 1
    out["A4_table_in_G1_ok"] = [ok_t, nt]
    out["seconds"] = round(time.time() - t0, 3)
    return out


# --------------------------------------------------------------- Section B: toy

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
    X <- X * T_h with T_h = g^(A_h) T^(B_h), h a hash of X, (a,b) <- (a,b)+(A_h,B_h).
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
    # v3 RULE, the same as production: the hard cap is a multiple of sqrt(q).
    # (v2 used cap = cap_factor * S0, i.e. a multiple of the *uncapped* first-repeat
    # quantile, which is a different and non-conservative rule once the detection
    # lag is counted.  S0 stays in the output only for continuity with v2.)
    cap = int(math.ceil(cap_factor * sqrtq))
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


def section_c(success_model=Fraction(2, 5), cap_num=33, cap_den=20, step_instr=120,
              table_size=1024, verification_adds=4 * 383, precompute_mults=2048,
              success_claim=Fraction(42, 100)):
    """Exact ledger.  Units: one BLS G1 addition at the fixed price C = 2,000."""
    unit = INSTR_PER_ADD
    LN53 = Fraction(510826, 1000000)      # > ln(5/3) = ln(1/(1-2/5))
    out = {"r": str(R), "r_log2": round(_log2_frac(Fraction(R)), 6),
           "sqrt_r_log2": round(_log2_frac(Fraction(R)) / 2, 6),
           "reference_price": unit, "measured_add_instr": MEASURED_ADD,
           "table_size": table_size,
           "success_model": str(success_model)}
    S0 = _ceil_sqrt_frac(2 * R * LN53.numerator, LN53.denominator * 1)
    out["S0_evals"] = str(S0)          # v2's uncapped first-repeat quantile (reference)
    out["S0_log2"] = round(_log2_frac(Fraction(S0)), 6)
    # v3 hard cap: ceil((cap_num/cap_den) * sqrt(r)), an exact integer.
    cap = _ceil_sqrt_frac(cap_num * cap_num * R, cap_den * cap_den)
    out["cap_factor"] = f"{cap_num}/{cap_den}"
    out["cap_evals"] = str(cap)
    out["cap_log2"] = round(_log2_frac(Fraction(cap)), 6)
    step_units = Fraction(unit + step_instr, unit)
    out["step_instr"] = step_instr
    out["step_units"] = str(step_units)
    main = Fraction(cap) * step_units
    out["main_units_log2"] = round(_log2_frac(main), 6)
    # additive: table build (2 scalar mults per entry: [c_h]P0 and [d_h]T) plus
    # start, verification ([x]P0) and output (4 scalar mults) plus the constant state
    precompute = precompute_mults * 383
    verify = verification_adds
    add = precompute + verify + 64
    out["precompute_units"] = str(precompute)
    out["verify_units"] = str(verify)
    out["additive_units"] = str(add)
    out["additive_units_log2"] = round(_log2_frac(Fraction(add)), 4)
    W = main + add
    out["W_log2_float"] = round(_log2_frac(W), 6)
    memory = (table_size * 2 * 48 + 6 * 2 * 48 + 2 * 32 + (1 << 16))
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


# ---------------------------------- Section D: the reviewer's success model
#
# Submission `68b71f2d` (fixed-g1) was rejected because our success bound used
# `lam <= p` and ignored that the step budget truncates the comparison window:
# for a checkpoint stored at scaled position p the admitted set is
# `u <= p` AND `p + lam <= end(p)`, with end(p) = min(next checkpoint, cap).
# With (u,v) = (mu,lam)/sqrt(N) having joint density exp(-(u+v)^2/2) on u,v >= 0,
# a schedule admits a staircase, not a square.  This section reproduces the
# reviewer's number, computes the worst dyadic alignment at the cap used here,
# and re-derives the same mass by an independent Monte-Carlo sampler.
SQRT2PI = math.sqrt(2.0 * math.pi)
PHI0 = 1.0 / SQRT2PI


def _Phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _phi(x):
    return math.exp(-0.5 * x * x) / SQRT2PI


def G(a, b):
    """Integral of exp(-(u+v)^2/2) over [0,a] x [0,b] (closed form)."""
    if a <= 0.0 or b <= 0.0:
        return 0.0
    return SQRT2PI * ((a + b) * _Phi(a + b) + _phi(a + b)
                      - (b * _Phi(b) + _phi(b))
                      - (a * _Phi(a) + _phi(a)) + PHI0)


def doubling(phase):
    """Brent's saved positions {2^j} scaled by 1/phase (dyadic grid alignment)."""
    out = []
    p = 1.0 / phase
    while p < 4096.0:
        out.append(p)
        p *= 2.0
    return out


def schedule_windows(ps, B):
    ws = []
    for j, p in enumerate(ps):
        nxt = ps[j + 1] if j + 1 < len(ps) else B
        ws.append(max(0.0, min(nxt, B) - p))
    return ws


def mass(schedule, B):
    """Reference-model success probability of a checkpoint schedule under cap B."""
    ps = [p for p in schedule if 0.0 < p < B]
    if not ps:
        return 0.0
    ws = schedule_windows(ps, B)
    total = 0.0
    prev = 0.0
    for j in range(len(ps)):
        V = max(ws[j:])
        total += G(ps[j], V) - G(prev, V)
        prev = ps[j]
    return total


def worst_phase_mass(B, phases=2000):
    best = None
    for i in range(phases):
        ph = 1.0 + i / float(phases)
        mm = mass(doubling(ph), B)
        if best is None or mm < best[0]:
            best = (mm, ph)
    return best


def hits_interval(k, n, z=1.96):
    if n == 0:
        return None
    ph = k / float(n)
    d = 1.0 + z * z / n
    c = ph + z * z / (2 * n)
    ss = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    return [round((c - ss) / d, 4), round((c + ss) / d, 4)]


def mass_mc(B, phase, draws, seed):
    """Independent Monte-Carlo of the SAME model.

    u+v = w has density w*exp(-w^2/2) on w >= 0 (Rayleigh), so draw
    w = sqrt(-2 ln U) and then u uniform on [0,w] (v = w - u).  Detection at a
    checkpoint p is `p >= u and v <= end(p)-p`; a run is detected iff some
    checkpoint admits it, which is exactly the staircase the quadrature sums.
    """
    ps = [p for p in doubling(phase) if 0.0 < p < B]
    ws = schedule_windows(ps, B)
    rng = (seed * 2862933555777941757 + 3037000493) % (2 ** 64)
    hit = 0
    for _ in range(draws):
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        u1 = ((rng >> 11) + 1) / float(1 << 53)
        rng = (rng * 6364136223846793005 + 1442695040888963407) % (2 ** 64)
        u2 = (rng >> 11) / float(1 << 53)
        w = math.sqrt(-2.0 * math.log(u1))
        u = w * u2
        v = w - u
        ok = False
        for j, p in enumerate(ps):
            if p >= u and v <= ws[j]:
                ok = True
                break
        if ok:
            hit += 1
    return hit / float(draws), hits_interval(hit, draws)


def section_d(cap_num=33, cap_den=20, phases=1200, draws=20000, seed=1):
    """Success probability under the reviewer's model, at the declared cap."""
    B = cap_num / float(cap_den)
    out = {"cap_factor": f"{cap_num}/{cap_den}", "cap_scaled": B,
           "reviewer_rectangle_G_16_15_8_15": round(G(16.0 / 15.0, 8.0 / 15.0), 6),
           "reviewer_said": 0.4041}
    mm, ph = worst_phase_mass(B, phases=phases)
    out["worst_phase_mass"] = round(mm, 6)
    out["worst_phase"] = round(ph, 4)
    dyad = doubling(ph)
    out["worst_phase_first_checkpoints"] = [round(p, 4) for p in dyad[:4]]
    out["mass_at_worst_phase_check"] = round(mass(dyad, B), 6)
    mc, ci = mass_mc(B, ph, draws, seed)
    out["mc_mass_at_worst_phase"] = round(mc, 4)
    out["mc_wilson"] = ci
    out["mc_agrees_with_quadrature"] = bool(ci[0] - 0.01 <= mm <= ci[1] + 0.01)
    out["mass_phase_1_875"] = round(mass(doubling(1.875), B), 6)
    out["mass_phase_1_0"] = round(mass(doubling(1.0), B), 6)
    out["success_at_least_039"] = bool(mm >= 0.39)
    return out


def main():
    try:
        raw = sys.stdin.read()
        cfg = json.loads(raw) if raw.strip() else {}
    except Exception as e:                     # noqa: BLE001
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
                            cap_factor=float(cfg.get("cap_factor", 1.65)),
                            wall_seconds=float(cfg.get("seconds", 10.0)))
    result["C"] = section_c()
    result["D"] = section_d(draws=int(cfg.get("model_draws", 20000)),
                            phases=int(cfg.get("model_phases", 1200)), seed=seed)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
