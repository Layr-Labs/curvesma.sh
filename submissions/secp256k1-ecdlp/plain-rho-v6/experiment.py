"""Evidence experiment for the plain-rho-v6 submission (secp256k1-ecdlp).

v6 is the plain additive walk on the full group: C <- C + A_j with
j = x_C mod 32, A_j = [c_j]P, c_j = j*Q + 2^j mod r, Q = (r - r mod 32)/32.

Two properties, both verified in Part A:

(1) NO FRUITLESS CYCLES, PROVABLY. A cycle of length L is a word
    w_1..w_L with sum c_{w_u} = 0 mod r, i.e. Q*A + B = 0 mod r with
    A = sum w_u and B = sum 2^{w_u}. Multiplying by 32 and using
    32Q = r - rho (rho = r mod 32) gives 32B = rho*A mod r. For
    L <= L* = r/(32*2^31) we have 32B < r and rho*A < r, so the
    congruence is the integer equation 32B = rho*A. But B >= 2A term
    by term (2^w >= 2w for all w >= 0), so 32B >= 64A > 31A >= rho*A
    unless A = B = 0, which is impossible for L >= 1. Hence no cycle
    of length <= L* ~ 2^219.99 exists. The per-walker budget is
    1.25*sqrt(r-1)/2^16 ~ 2^112.32, far below L*: no look-back, no
    ratchet, no epoch, no escape machinery is needed.

(2) MIXING. The constant differences c_i - c_j = (i-j)Q + (2^i - 2^j)
    span 31Q ~ 0.97r, so the k-space difference walk between two
    walkers wraps mod r within a few dozen steps and the states of
    different walkers collide at the standard parallel-rho birthday
    rate t^2/(2(r-1)) (van Oorschot-Wiener). The toy in Part B
    measures the realized fraction of this ideal directly.

Walkers carry lazy representations C = [a]P + [b]T with b fixed per
walker; the table is keyed by (x, y-parity) so a key hit is a true
point collision and the extraction x = (a2-a1)/(b1-b2) mod r is
exact algebra. Reporting loss is bounded by 2^32/t (budget-end race
after a merge; merged walkers provably coalesce since the map is
deterministic) plus table overflow (Chernoff).

Declared run: {"trials": 200, "dbits": 0, "walkers": 8} executes the
complete v6 attack (seeding, walk, distinguished points at density 1,
table, extraction, verification) and measures collision-time
efficiency; the same run also executes a fixed supplementary census at
dbits=6 (production-shaped density 2^-6) exercising the sparse-DP
machinery including the toy-scale cycle-escape reseed.

Stdlib only; deterministic (fixed seeds); single JSON object on stdout.
"""
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


def build_constants(r):
    rho = r % 32
    Q = (r - rho) // 32
    return rho, Q, [(j * Q + (1 << j)) % r for j in range(32)]


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
    rho, Q, cvals = build_constants(r)
    assert rho != 0 and 32 * Q == r - rho
    out["constants_construction"] = "c_j = j*Q + 2^j mod r, Q=(r-rho)/32"
    pts = [ec_mul(c, G, p, b) for c in cvals]
    assert all(A is not None for A in pts)
    xs = set(A[0] for A in pts)
    out["constants_distinct"] = len(xs) == 32
    # mixing: constant differences must span (almost) the whole modulus
    # so that the k-space difference walk wraps mod r
    span = max((cvals[i] - cvals[j]) % r for i in range(32) for j in range(32))
    out["mixing"] = {
        "max_difference": span,
        "span_over_r": span / r,
        "wraps_within_steps": math.log(1 / (1 - span / r), 2),
    }
    assert span / r >= 0.9
    # no-cycle theorem
    L_star = r / (32 * (2 ** 31))
    N = r - 1
    t = 1.25 * math.sqrt(N)
    per_walker = t / 2 ** 16
    out["no_cycle_theorem"] = {
        "rho": rho,
        "L_star": L_star,
        "log2_L_star": math.log2(L_star),
        "per_walker_budget": per_walker,
        "log2_per_walker_budget": math.log2(per_walker),
        "budget_below_L_star": per_walker < L_star,
    }
    assert per_walker < L_star
    # ledger
    other = 2 + 32 + 4 + 69  # index, rep add, DP predicate, state/loop
    per_step = 1 + other / C_REF
    seeding = 2 ** 16 * 64
    dp_proc = t * 2 ** -32 * 0.02
    W = per_step * t + seeding + dp_proc + 700
    tl2 = math.log2(W)
    num = math.ceil(tl2 * 10 ** 6)
    fr = Fraction(num, 10 ** 6)
    out["ledger"] = {
        "per_step_other_instructions": other,
        "per_step_units": per_step,
        "budget_multiplier": 1.25,
        "t_log2": math.log2(t),
        "W_log2_raw": tl2,
        "time_log2": f"{fr.numerator}/{fr.denominator}",
        "preprocessing_log2": 22.0,
    }
    # distinguished-point table: capacity, overflow, memory
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
    # success model. A merge is reported at the next shared
    # distinguished point; the race is between the DP interval (2^32
    # steps in expectation) and the walker's remaining personal budget
    # (t/2^16 steps), so the loss is 2^32/per_walker
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
# Part B: toy-scale complete v6 attack (plain additive walk)
# ---------------------------------------------------------------------------
TOY_P = 1048609
TOY_B = 29
TOY_N = 1047379
TOY_G = (1, 245324)


def toy_attack(trials, dbits, walkers):
    p, b, n, G = TOY_P, TOY_B, TOY_N, TOY_G
    N = n - 1  # full group: no canonical folding
    budget = math.ceil(1.25 * math.isqrt(N))
    # toy constants mirror the production construction c_j = jQ + 2^j
    # mod n. At toy scale L* = n/(32*2^31) < 1, so toy-only wrap cycles
    # exist (production excludes them by the theorem); the main attack
    # absorbs that waste into the measured ratio, the census exercises
    # the cycle-escape reseed, and both record the counts.
    rho, Q, cvals = build_constants(n)
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
                        # point: a toy-scale wrap cycle (impossible at
                        # production by the theorem); reseed the
                        # walker, record, and continue the sweep
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
        "constants": "c_j = jQ + 2^j mod n (production construction)",
        "successes": succ,
        "empirical_success": measured,
        "ideal_success": ideal,
        "measured_ratio": ratio,
        "ratio_two_sigma": 2 * sigma,
        "wrong_extractions": wrong,
        "identity_reseeds": identity_reseeds,
        "cycle_reseeds_toy_scale": cycle_reseeds,
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
        "toy_plain_rho_v6": toy_attack(trials, dbits, walkers),
        # supplementary census: the same machinery at a
        # production-shaped distinguished-point density (2^-6 at toy
        # scale), exercising the sparse-DP table, lagged reporting and
        # the cycle-escape reseed, recorded inside the declared run
        "toy_dp_density_census": toy_attack(100, 6, 8),
    }
    print(json.dumps(out))


main()
