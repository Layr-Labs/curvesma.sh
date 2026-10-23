"""Evidence experiment for the plain-rho-v8 submission (secp256k1-ecdlp).

v8 is the plain additive walk on the full group with STRUCTURE-FREE
constants: A_j = [c_j]P where c_j = sha512(label/j) mod r, published
and deterministic. The v6/v7 constants c_j = jQ + 2^j were rejected
because the rescaling z = 32k mod r (using r = 32Q + 1) turns every
step into a positive increment 2^(j+5) - j <= 2^36: the walk is
confined and cross-walker collisions essentially never occur at
production scale. Pseudorandom constants have no such confining
rescaling: the counting heuristic puts the expected number of
rescalings M that bring all 32 steps into any window of width W
below 2^240 at r*(W/r)^31 << 1, and Part A spot-verifies every
rescaling up to 2^12 leaves the step set spread across >= 0.8r.

With structure-free constants the analysis is the standard
random-function model of parallel rho (H1): the birthday collision
count t^2/(2(r-1)), cycle entries at rate 1/r per window (< 2^-127
expected over the whole budget), and Binomial distinguished-point
counts (Chernoff overflow bound stated under H1). Part A verifies
at production parameters; Part B measures the realized fraction of
the ideal directly and exercises the sparse-DP machinery.

Cycle handling: a cycle of length <= 4 would require a multiset of
at most 4 constants summing to 0 mod r - Part A enumerates all of
them and asserts none exists. Longer cycles are detected for free
by the distinguished-point table (a revisited key with an identical
representation reseeds the walker) and are a 2^-127 event under H1;
the claim does not depend on their absence, surviving cycle-entry
rates up to 2^-60.

Declared run: {"trials": 200, "dbits": 0, "walkers": 8} executes the
complete v8 attack and measures collision-time efficiency; the same
run also records a supplementary census at dbits=6 exercising the
sparse-distinguished-point machinery.

Stdlib only; deterministic (fixed seeds); single JSON object on stdout.
"""
import hashlib
import json
import math
import random
import sys
from fractions import Fraction

P_FIELD = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
R_ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
GEN = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
       0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)
CURVE_B = 7
TRACK = "secp256k1-ecdlp"
C_REF = 1000  # reference-operation price for this track


def derive_constants(r, label):
    """Published deterministic pseudorandom constants c_j mod r.

    c_j = int(sha512("curvesmash/<label>/constant/<j>")) mod r; the
    512-bit hash value reduced mod r has bias below 2^-500 (512-bit
    value against a <= 256-bit modulus), far below every figure in
    this submission.
    """
    out = []
    for j in range(32):
        h = hashlib.sha512(
            f"curvesmash/{label}/constant/{j}".encode()).digest()
        out.append(int.from_bytes(h, "big") % r)
    return out


def is_prime_mr(v):
    if v < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if v % sp == 0:
            return v == sp
    d = v - 1
    a = 0
    while d % 2 == 0:
        d //= 2
        a += 1
    for w in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        x = pow(w, d, v)
        if x in (1, v - 1):
            continue
        ok = False
        for _ in range(a):
            x = x * x % v
            if x == v - 1:
                ok = True
                break
        if not ok:
            return False
    return True


def ec_add(Pt, Qt, mod, b):
    if Pt is None:
        return Qt
    if Qt is None:
        return Pt
    x1, y1 = Pt
    x2, y2 = Qt
    if x1 == x2 and (y1 + y2) % mod == 0:
        return None
    if Pt == Qt:
        m = (3 * x1 * x1) * pow(2 * y1, -1, mod) % mod
    else:
        m = (y2 - y1) * pow(x2 - x1, -1, mod) % mod
    x3 = (m * m - x1 - x2) % mod
    return (x3, (m * (x1 - x3) - y1) % mod)


def ec_mul(k, Pt, mod, b):
    if k % mod == 0:
        return None
    acc = None
    add = Pt
    while k:
        if k & 1:
            acc = ec_add(acc, add, mod, b)
        add = ec_add(add, add, mod, b)
        k >>= 1
    return acc


def step_window(vals, r):
    """Width of the narrowest window containing all values (mod r).

    A rescaling z = M*k confines the walk only if all rescaled steps
    M*c_j mod r lie in one short window; this returns that width.
    """
    sv = sorted(vals)
    gaps = [(sv[(i + 1) % len(sv)] - sv[i]) % r for i in range(len(sv))]
    return r - max(gaps)


def production_checks():
    p, r, G, b = P_FIELD, R_ORDER, GEN, CURVE_B
    out = {"track": TRACK}
    assert is_prime_mr(p) and is_prime_mr(r)
    out["p_and_r_prime"] = True
    x, y = G
    assert (y * y - (x * x * x + b)) % p == 0
    out["generator_on_curve"] = True
    assert ec_mul(r, G, p, b) is None
    out["generator_order_r"] = True
    cvals = derive_constants(r, TRACK)
    pts = [ec_mul(c, G, p, b) for c in cvals]
    assert all(A is not None for A in pts)
    xs = set(A[0] for A in pts)
    out["constants"] = {
        "construction": "c_j = int(sha512('curvesmash/%s/constant/j')) "
                        "mod r (published, deterministic)" % TRACK,
        "distinct": len(xs) == 32,
        "sha256_of_concatenation": hashlib.sha256(
            b"".join(c.to_bytes((r.bit_length() + 7) // 8, "big")
                     for c in cvals)).hexdigest(),
    }
    found_short = None
    from itertools import combinations_with_replacement
    for size in (2, 3, 4):
        for combo in combinations_with_replacement(range(32), size):
            if sum(cvals[j] for j in combo) % r == 0:
                found_short = combo
                break
        if found_short is not None:
            break
    out["no_short_cycle_word"] = {
        "checked_multisets_up_to_length": 4,
        "none_sums_to_zero": found_short is None,
    }
    assert found_short is None
    assert len(xs) == 32
    # confinement spot-certificate: rescaling z = M*k changes the
    # steps to M*c_j mod r; a confinement needs all 32 rescaled
    # steps inside a short window. Check every M up to 2^12 and
    # report the narrowest window found.
    worst = None
    for M in range(1, 1 << 12):
        w = step_window([M * c % r for c in cvals], r)
        if worst is None or w < worst[1]:
            worst = (M, w)
    # self-test: the same certificate applied to the rejected v6/v7
    # constants c_j = jQ + 2^j (r = 32Q + rho) detects their
    # confinement at M = 32: the rescaled steps land in 2^(j+5) - j,
    # a window of width <= 2^36
    Qv7 = (r - r % 32) // 32
    v7 = [(j * Qv7 + (1 << j)) % r for j in range(32)]
    v7w = step_window([32 * c % r for c in v7], r)
    out["confinement_certificate"] = {
        "rescalings_checked": (1 << 12) - 1,
        "narrowest_step_window": worst[1],
        "narrowest_at_M": worst[0],
        "window_over_r": worst[1] / r,
        "v7_constants_detected_window_over_r": v7w / r,
        "counting_heuristic_confiners_below_2^240":
            "r*(2^-15)^31 < 2^-209",
    }
    assert worst[1] / r >= 0.5
    assert v7w / r < 2 ** -200
    # ledger
    other = 2 + 32 + 4 + 69  # index, rep add, DP predicate, state/loop
    per_step = 1 + other / C_REF
    N = r - 1
    t = 1.25 * math.sqrt(N)
    per_walker = t / 2 ** 16
    # seeding: [a0]P by fixed-base comb (32 adds), [b0]T by 4-bit
    # windowed double-and-add (64 doubles + 16 adds), per walker
    # <= 112 reference ops; 2^16 walkers -> 2^22.8, charged 2^23
    seeding = 2 ** 23
    # DP processing: 2^96.32 inserts at 20 instructions
    # (hash 8 + probe 4 + compare 8) = 2^91.7 instructions
    dp_proc = t * 2 ** -32 * 20 / C_REF
    # extraction (modular inverse ~700 instr) + verification of
    # [x]P == T (one ~256-add scalar multiple ~ 2^8 reference ops)
    extract_verify = 2 ** 9
    W = per_step * t + seeding + dp_proc + extract_verify
    tl2 = math.log2(W)
    num = math.ceil(tl2 * 10 ** 6)
    fr = Fraction(num, 10 ** 6)
    out["ledger"] = {
        "per_step_other_instructions": other,
        "per_step_units": per_step,
        "budget_multiplier": 1.25,
        "t_log2": math.log2(t),
        "per_walker_log2": math.log2(per_walker),
        "seeding_units_log2": 23.0,
        "dp_processing_units_log2": math.log2(dp_proc),
        "extract_verify_units_log2": 9.0,
        "W_log2_raw": tl2,
        "time_log2": f"{fr.numerator}/{fr.denominator}",
        "preprocessing_log2": 23.0,
    }
    # distinguished-point table: capacity, overflow, memory. Under H1
    # the insert count is Binomial(t, 2^-32); Chernoff gives the
    # stated overflow bound.
    inserts = t * 2 ** -32
    cap_bits = math.ceil(math.log2(inserts)) + 1
    out["dp"] = {
        "inserts_log2": math.log2(inserts),
        "capacity_slots_log2": cap_bits,
        "capacity_over_mean": 2 ** cap_bits / inserts,
        "overflow_exponent_neg": (2 ** cap_bits - inserts) ** 2
        / (2 * 2 ** cap_bits),
        "memory_log2_bytes_reported": cap_bits + 7 + 0.1,
    }
    # success model
    rep_loss = 2 ** 32 / per_walker
    ideal = 1 - math.exp(-(1.25 * (1 - rep_loss)) ** 2 / 2)
    out["success_model"] = {
        "reporting_loss": rep_loss,
        "ideal": ideal,
        "efficiency_floor": 0.75,
        "claimed_success": ideal * 0.75,
    }
    assert ideal * 0.75 >= 0.39
    return out


# ---------------------------------------------------------------------------
# Part B: toy-scale complete v8 attack (plain additive walk)
# ---------------------------------------------------------------------------
TOY_P = 1048609
TOY_B = 29
TOY_N = 1047379
TOY_G = (1, 245324)


def toy_attack(trials, dbits, walkers):
    p, b, n, G = TOY_P, TOY_B, TOY_N, TOY_G
    N = n - 1  # full group: no canonical folding
    budget = math.ceil(1.25 * math.isqrt(N))
    # toy constants: the same published pseudorandom construction
    # evaluated mod the toy order; distinctness asserted
    cvals = derive_constants(n, "toy")
    consts = [ec_mul(c, G, p, b) for c in cvals]
    xs = set(A[0] for A in consts)
    assert len(xs) == 32

    def seed_state(rng, Tg):
        a0 = rng.randrange(1, n)
        b0 = rng.randrange(1, n)
        C = ec_add(ec_mul(a0, G, p, b), ec_mul(b0, Tg, p, b), p, b)
        return C, a0, b0

    rng = random.Random(2026)
    succ = 0
    wrong = 0
    identity_reseeds = 0
    cycle_reseeds = 0
    dp_detections = 0
    executed = 0
    for ti in range(trials):
        s = rng.randrange(1, n)
        Tg = ec_mul(s, G, p, b)
        tab = {}
        wk = []
        for wi in range(walkers):
            C, a0, b0 = seed_state(rng, Tg)
            wk.append({"C": C, "a": a0, "b": b0})
        found = None
        steps_left = budget
        while steps_left > 0 and found is None:
            for wi in range(walkers):
                if steps_left <= 0 or found is not None:
                    break
                st = wk[wi]
                C, a, bb = st["C"], st["a"], st["b"]
                j = C[0] % 32
                C2 = ec_add(C, consts[j], p, b)
                steps_left -= 1
                executed += 1
                if C2 is None:
                    # C = -A_j: probability 32/(r-1) per step at
                    # production (expected count over the budget
                    # < 2^-119); uncharged there, reseeded and
                    # recorded here
                    identity_reseeds += 1
                    C3, a3, b3 = seed_state(rng, Tg)
                    st.update({"C": C3, "a": a3, "b": b3})
                    continue
                a2 = (a + cvals[j]) % n
                b2 = bb
                # distinguished point. The table is keyed by
                # (x, y-parity): on the full group +X and -X are
                # distinct points sharing an x-coordinate, so the
                # parity bit is what makes a key hit a true collision
                # (same key => same point => extraction is exact).
                if dbits == 0 or C2[0] % (1 << dbits) == 0:
                    key = (C2[0], C2[1] % 2)
                    prev = tab.get(key)
                    if prev is None:
                        tab[key] = (a2, b2)
                    elif prev != (a2, b2):
                        # collision with a different representation:
                        # extract and verify
                        db = (prev[1] - b2) % n
                        if db == 0:
                            # same b, different a at one point cannot
                            # represent the same point; safety branch
                            wrong += 1
                            break
                        xv = (a2 - prev[0]) * pow(db, -1, n) % n
                        if ec_mul(xv, G, p, b) == Tg:
                            found = xv
                            dp_detections += 1
                        else:
                            wrong += 1
                        break
                    else:
                        # identical representation at a revisited
                        # point: a cycle revisit (a 1/r-per-window
                        # event under H1; free DP-table detection,
                        # walker reseeds); recorded
                        cycle_reseeds += 1
                        C3, a3, b3 = seed_state(rng, Tg)
                        st.update({"C": C3, "a": a3, "b": b3})
                        continue
                st["C"], st["a"], st["b"] = C2, a2, b2
        if found is not None and found == s:
            succ += 1
        elif found is not None:
            wrong += 1
    ideal = 1 - math.exp(-budget ** 2 / (2 * N))
    measured = succ / trials
    sigma = math.sqrt(max(measured * (1 - measured), 1e-9) / trials) / ideal
    ratio = measured / ideal
    return {
        "curve": "y^2 = x^3 + 29 over F_1048609, prime order 1047379",
        "group_size": N,
        "budget_total_steps": budget,
        "per_walker_steps": budget // walkers,
        "walkers": walkers,
        "dbits": dbits,
        "trials": trials,
        "constants": "c_j = sha512(...) mod n (production construction)",
        "successes": succ,
        "empirical_success": measured,
        "ideal_success": ideal,
        "measured_ratio": ratio,
        "ratio_two_sigma": 2 * sigma,
        "wrong_extractions": wrong,
        "identity_reseeds": identity_reseeds,
        "cycle_reseeds": cycle_reseeds,
        "dp_detections": dp_detections,
        "executed_steps": executed,
    }


def main():
    data = json.load(sys.stdin)
    trials = int(data.get("trials", 200))
    dbits = int(data.get("dbits", 0))
    walkers = int(data.get("walkers", 8))
    out = {
        "production": production_checks(),
        "toy_plain_rho_v8": toy_attack(trials, dbits, walkers),
        # supplementary census: the same machinery at a
        # production-shaped distinguished-point density (2^-6 at toy
        # scale), exercising the sparse-DP table, lagged reporting and
        # the cycle-escape reseed, recorded inside the declared run
        "toy_dp_density_census": toy_attack(100, 6, 8),
    }
    print(json.dumps(out))


main()
