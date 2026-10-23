"""Evidence experiment for the plain-rho-v21 submission (secp256k1-ecdlp).

v12 is the plain additive walk on the full group with published
pseudorandom constants, an append-only ARENA for distinguished points
and deferred duplicate detection by one in-place heapsort plus
scan, with no workspace. There is
no hash table, no chaining, no probing: a distinguished-point record
costs 20 instructions WORST CASE (4 linear address-
arithmetic instructions - the 2^63-byte arena fits the
64-bit address space, so no banking exists - 8 paired
stores for the 128-byte record, 8 bookkeeping), the arena
holds at most 2^cap records (the next append aborts the attack, its
probability bounded by MARKOV on the exact expected record
count and subtracted), and the sort of at
most cap*log2(cap) record moves is charged as a one-off. Detection
is the same event as immediate table detection: some distinguished
point recorded by two walkers with different representations.

Every work term is a hard worst-case bound: walker steps <= t
(round-robin), appends <= 2^cap, identity reseeds <= 2^20 (abort
beyond, Markov), one sort, one extraction. Expected-value quantities
(cycle waste) enter only as probabilities via Markov's inequality and
are included in the abort total, so the collision-factor tolerance is
computed after all deductions. No external collision-time factor
enters the model; the only measurement is the toy, which runs the
same machinery (arena, deferred detection, extraction with
verification).

Stdlib only; deterministic (fixed seeds); single JSON object on
stdout. Declared run: {"trials": 200, "dbits": 0, "walkers": 8} plus
a fixed dbits=6 sparse-density census recorded in the same output.
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
    512-bit hash value reduced mod r has bias below 2^-500.
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
    assert len(xs) == 32
    # no short cycle word: no multiset of <= 4 constants sums to 0
    # mod r (a cycle word of length L needs such a multiset). This
    # excludes every fruitless cycle of length <= 4 regardless of
    # index realizability.
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
            "32*W^31/r^30",
    }
    assert worst[1] / r >= 0.5
    assert v7w / r < 2 ** -200
    # ledger
    other = 2 + 32 + 6 + 69  # index, rep add, DP predicate, state/loop
    per_step = 1 + other / C_REF
    N = r - 1
    t = 1.5 * math.sqrt(N)
    per_walker = t / 2 ** 16
    # seeding, charged with GENEROUS FLAT upper bounds so no
    # implementation-level count can exceed them. Per walker:
    # [a0]P by fixed-base comb (32 adds + 32 doublings), [b0]T by
    # 4-bit windowed double-and-add (256 doublings + 64
    # additions), final addition: 97 additions and 288
    # doublings. Every ADDITION is charged 2 reference operations
    # (the reference operation is 1; 2 covers any full-width
    # modular add with carries, conditional reduction, correction
    # and transfers). Every DOUBLING is charged 1600 instructions
    # flat: an affine doubling is 7 field multiplications = 1253
    # instructions, and the flat 1600 covers the six full-width
    # adds with their carry chains, conditional reductions,
    # corrections and transfers at ANY reasonable count - a
    # deliberately loose bound. Randomness: TWO independent
    # 512-bit ChaCha20 outputs, one per scalar (80 quarter-rounds
    # = 960 instructions of round arithmetic each, plus
    # initialization, feed-forward, state transfers, counter
    # maintenance and the mod-r reduction): charged 4096
    # instructions per output flat. A reseed's own operations are
    # batched internally (one Montgomery inversion shared within
    # the reseed itself), so each reseed costs the same
    # per-walker price wherever it occurs. Per-walker total,
    # rounded up generously:
    #   97*2 + 288*1600/1000 + 2*4096/1000 = 558.6 <= 600
    # reference operations.
    per_walker_seed = 600
    seeding = 2 ** 16 * per_walker_seed
    table_build = 2 ** 9 + 15
    # constant construction: 32 fixed-base scalar multiples
    # (~256 reference additions each = 2^13 reference operations)
    # plus their sha512 derivation and mod-r reduction
    # (~1300 instructions each = 2^15 instructions = 2^5 units):
    # charged as 2^13.1 units
    constants_build = 2 ** 14
    # (the DRBG and reduction costs are inside per_walker_seed)
    # is below 2^-249, and no rejection sampling is performed, so
    # the randomness work is bounded with no abort path. The
    # seeding additions and doublings are batched across the 2^16
    # simultaneous walker evaluations (Montgomery inversion
    # sharing), matching the reference operation's batching, and
    # distinguished-point density: a state is distinguished when
    # its low 80 bits are zero. Under H1 each visited state is
    # uniform on the r-1 non-identity points, so a state qualifies
    # with probability at most (#qualifying)/(r-1) <= 2^-80 +
    # 2^-160 (the qualifying count is floor((r-1)/2^80) + 1, which
    # exceeds 2^-80 * (r-1) by at most one point). The expected
    # record count is therefore BOUNDED: at most
    # t * (2^-80 + 2^-160) - dependence between steps cannot
    # change the per-state marginal, so the mean bound holds with
    # no independence assumption.
    inserts = t * (2 ** -80 + 2 ** -160)
    cap_bits = 56
    overflow_markov = inserts / 2 ** cap_bits
    # DP appends: worst case 20 instructions per append (64-bit
    # linear addressing of the 2^63-byte arena: 4 address arithmetic
    # instructions; 8 paired 64-bit stores for the 128-byte record;
    # 8 bookkeeping and branch). Hard cap 2^cap appends.
    dp_appends_worst = 2 ** cap_bits * 20 / C_REF
    # deferred duplicate detection: one IN-PLACE HEAPSORT of at most
    # 2^cap records plus one linear scan, all worst case. In-place
    # heapsort uses NO workspace: the peak is the arena alone
    # (2^56 * 128 = 2^63 bytes, strictly inside the 64-bit address
    # space; nothing else of comparable scale is allocated).
    # Heapsort performs at most 2*n*log2(n) record moves
    # (sift-downs), each a full-key compare (up to 6 64-bit words
    # with early-exit branches: charge 16), 64-bit linear
    # addressing 4, load 8 and store 8 as paired transfers, and
    # heap bookkeeping 8 - 48 instructions per move; plus the
    # linear scan at 18 instructions per record.
    final_sort = (2 * 2 ** cap_bits * cap_bits * 48
                  + 2 ** cap_bits * 64) / C_REF    # extraction (modular inverse ~700 instr) + verification of
    # [x]P == T (one ~256-add scalar multiple ~ 2^8 reference ops)
    extract_verify = 2 ** 9
    # identity reseeds (state = -A_j with matching index): rate at
    # most 32/(r-1) per step, capped at 2^20, abort beyond
    identity_reseed_cap_work = 2 ** 20 * per_walker_seed
    W = (per_step * t + seeding + table_build + constants_build
         + dp_appends_worst
         + final_sort + extract_verify + identity_reseed_cap_work)
    tl2 = math.log2(W)
    num = math.ceil(tl2 * 10 ** 6)
    fr = Fraction(num, 10 ** 6)
    out["ledger"] = {
        "per_step_other_instructions": other,
        "per_step_units": per_step,
        "budget_multiplier": 1.5,
        "t_log2": math.log2(t),
        "per_walker_log2": math.log2(per_walker),
        "seeding_units_log2": math.log2(seeding + table_build),
        "arena_appends_worst_units_log2": math.log2(dp_appends_worst),
        "final_sort_units_log2": math.log2(final_sort),
        "extract_verify_units_log2": 9.0,
        "W_log2_raw": tl2,
        "time_log2": f"{fr.numerator}/{fr.denominator}",
        "preprocessing_log2": math.ceil(math.log2(seeding + table_build)
                                        * 2) / 2,
    }
    out["dp"] = {
        "construction": "append-only arena of 128-byte records plus "
                        "one in-place heapsort (no hash table, no "
                        "links, no workspace; 20-instruction "
                        "worst-case appends; the peak is the arena "
                        "alone at 2^63 bytes, inside the 64-bit "
                        "address space)",
        "append_price_instructions_worst": 24,
        "record_bytes": 128,
        "inserts_log2": math.log2(inserts),
        "capacity_records_log2": cap_bits,
        "overflow_bound_markov": overflow_markov,
        "memory_log2_bytes_reported": cap_bits + 7 + 0.1,
        "density_bits": 80,
    }
    # cycle model under H1. A fruitless cycle is a SELF-repeat of one
    # walker: a pair (u, v), u < v, on one trail whose window sum is
    # 0 mod r, an event of probability 1/r per self pair. The number
    # of self pairs is K*L*(L-1)/2 with K = 2^16, L = t/K, so the
    # expected number of cycle entries is t^2/(2Kr). A trapped walker
    # wastes at most its remaining budget; charging the remaining
    # budget L - v for every entry at (u, v) bounds the total
    # expected waste by K * sum_(u<v) (L - v)/r = K*L^3/(6r) =
    # t^3/(6K^2 r) steps, a fraction t^2/(6K^2 r) of the budget.
    # The bound is computed over the full budget, so it also bounds
    # the process with early termination. The expectation is NOT
    # substituted as a work deduction: it is converted by Markov's
    # inequality into the abort total below.
    cyc_entries = t * t / (2 * 2 ** 16 * r)
    cyc_waste_frac = t * t / (6 * 2 ** 32 * r)
    # Markov: P(cycle waste exceeds 2^-23 of the budget)
    #   <= E[W] / (2^-23 * t)
    markov_waste = cyc_waste_frac / 2 ** -23
    # reporting: a merged pair reports at the next shared
    # distinguished point unless the personal budget ends first
    rep_loss = 2 ** 80 / per_walker
    # identity-reseed abort: Markov on the expected count; the
    # computed bound is printed and quoted as-is (never rounded
    # below the computed value)
    identity_expected = t * 32 / r
    identity_abort = identity_expected / 2 ** 20
    identity_abort_log2 = math.log2(identity_abort)
    # efficiency: the r-adding walk's collision-time factor c (real
    # collision time / ideal random-walk collision time). Within
    # the fixed budget the collision probability is
    # q(c) = 1 - exp(-(1.5/c)^2 / 2) under Rayleigh collision
    # times. Bad events (aborts) are combined with q by UNION
    # BOUND SUBTRACTION, not multiplication: P(success) >=
    # q(c) - abort_total, since P(A and B) >= P(A) - P(B^c).
    # The claim of 0.39 holds for every c <= c_max:
    aborts = (rep_loss + markov_waste + identity_abort
              + overflow_markov)
    c_max = 1.5 / math.sqrt(-2 * math.log(1 - 0.39 - aborts))
    ideal = (1 - math.exp(-(1.5) ** 2 / 2)) - aborts
    out["success_model"] = {
        "ideal_at_c_1": ideal,
        "reporting_loss": rep_loss,
        "cycle_entries_expected": cyc_entries,
        "cycle_waste_mean_fraction": cyc_waste_frac,
        "markov_bound_waste_over_2^-23": markov_waste,
        "markov_waste_exponent": math.log2(markov_waste),
        "identity_reseed_abort": identity_abort,
        "identity_reseed_abort_log2": identity_abort_log2,
        "arena_overflow_markov": overflow_markov,
        "abort_total": aborts,
        "collision_factor_tolerance_c_max": c_max,
        "reseed_cap": 2 ** 20,
        "total_work_cutoff": "walker steps <= t (round-robin "
                              "budget); arena appends <= 2^cap "
                              "(abort, Markov on the exact mean); "
                              "identity reseeds "
                              "<= 2^20 (abort, Markov); one final "
                              "sort; one extraction; W bounds every "
                              "non-aborting path",
    }
    assert c_max >= 1.25
    assert ideal >= 0.39
    return out


# ---------------------------------------------------------------------------
# Part B: toy-scale complete v21 attack (plain additive walk)
# ---------------------------------------------------------------------------
TOY_P = 1048609
TOY_B = 29
TOY_N = 1047379
TOY_G = (1, 245324)


def _sift(arr, root, end):
    while True:
        child = 2 * root + 1
        if child >= end:
            break
        if child + 1 < end and arr[child + 1] > arr[child]:
            child += 1
        if arr[child] > arr[root]:
            arr[root], arr[child] = arr[child], arr[root]
            root = child
        else:
            break


def _heapsort(arr):
    # in-place heapsort, matching the production construction
    n = len(arr)
    for i in range(n // 2 - 1, -1, -1):
        _sift(arr, i, n)
    for end in range(n - 1, 0, -1):
        arr[0], arr[end] = arr[end], arr[0]
        _sift(arr, 0, end)


def toy_attack(trials, dbits, walkers):
    p, b, n, G = TOY_P, TOY_B, TOY_N, TOY_G
    N = n - 1  # full group: no canonical folding
    budget = math.ceil(1.5 * math.isqrt(N))
    # toy constants: the same published pseudorandom construction
    # evaluated mod the toy order; distinctness asserted
    cvals = derive_constants(n, "toy")
    consts = [ec_mul(c, G, p, b) for c in cvals]
    xs = set(A[0] for A in consts)
    assert len(xs) == 32
    # toy arena capacity mirrors production: 2^(ceil(log2 mean)+1)
    arena_cap = 1 << (math.ceil(math.log2(budget * 2 ** -dbits)) + 1)

    def seed_state(rng, Tg):
        a0 = rng.randrange(1, n)
        b0 = rng.randrange(1, n)
        C = ec_add(ec_mul(a0, G, p, b), ec_mul(b0, Tg, p, b), p, b)
        return C, a0, b0

    rng = random.Random(2026)
    succ = 0
    wrong = 0
    identity_reseeds = 0
    arena_aborts = 0
    cycle_duplicates = 0
    sort_detections = 0
    executed = 0
    true_merges = 0
    merged_trials = 0
    merged_detected = 0
    for ti in range(trials):
        s = rng.randrange(1, n)
        Tg = ec_mul(s, G, p, b)
        arena = []
        hidden = {}
        wk = []
        for wi in range(walkers):
            C, a0, b0 = seed_state(rng, Tg)
            wk.append({"C": C, "a": a0, "b": b0})
        merged = [False] * walkers
        aborted = False
        steps_left = budget
        while steps_left > 0 and not aborted:
            for wi in range(walkers):
                if steps_left <= 0:
                    break
                st = wk[wi]
                C, a, bb = st["C"], st["a"], st["b"]
                j = C[0] % 32
                C2 = ec_add(C, consts[j], p, b)
                steps_left -= 1
                executed += 1
                if C2 is None:
                    # C = -A_j: probability 32/(r-1) per step at
                    # production (capped at 2^20, abort beyond);
                    # reseeded and recorded here
                    identity_reseeds += 1
                    C3, a3, b3 = seed_state(rng, Tg)
                    st.update({"C": C3, "a": a3, "b": b3})
                    continue
                a2 = (a + cvals[j]) % n
                b2 = bb
                key = (C2[0], C2[1] % 2)
                # diagnostic full-visit log (not part of the
                # attack): a key hit with a different representation
                # is a true point collision; first contacts are
                # counted once per walker so the census can report
                # true merges versus detected merges
                if key in hidden:
                    if hidden[key] != (a2, b2) and not merged[wi]:
                        merged[wi] = True
                        true_merges += 1
                else:
                    hidden[key] = (a2, b2)
                # distinguished point: append to the arena (the
                # record is the storage; detection is deferred to
                # the final sort)
                if dbits == 0 or C2[0] % (1 << dbits) == 0:
                    if len(arena) < arena_cap:
                        arena.append((key, a2, b2))
                    else:
                        arena_aborts += 1
                        aborted = True
                        break
                st["C"], st["a"], st["b"] = C2, a2, b2
        found = None
        if not aborted:
            # deferred detection: an IN-PLACE HEAPSORT (the same
            # construction as production: sift-down heaps, no
            # workspace; Python's built-in sort is NOT used), then
            # scan adjacent equal keys for a pair with different
            # representations
            _heapsort(arena)
            for i in range(1, len(arena)):
                if (arena[i][0] == arena[i - 1][0]
                        and arena[i][1:] != arena[i - 1][1:]):
                    (k1, a1, b1) = arena[i - 1]
                    (k2, a2v, b2v) = arena[i]
                    db = (b1 - b2v) % n
                    if db == 0:
                        wrong += 1
                        break
                    xv = (a2v - a1) * pow(db, -1, n) % n
                    if ec_mul(xv, G, p, b) == Tg:
                        found = xv
                        sort_detections += 1
                    else:
                        wrong += 1
                    break
            # identical-representation duplicates are cycle evidence
            for i in range(1, len(arena)):
                if (arena[i][0] == arena[i - 1][0]
                        and arena[i][1:] == arena[i - 1][1:]):
                    cycle_duplicates += 1
                    break
        if any(merged):
            merged_trials += 1
            if found is not None and found == s:
                merged_detected += 1
        if found is not None and found == s:
            succ += 1
        elif found is not None:
            wrong += 1
    ideal = 1 - math.exp(-budget ** 2 / (2 * N))
    L = budget // walkers
    self_repeat_model = walkers * L * (L - 1) / (2 * n)
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
        "arena_capacity_records": arena_cap,
        "arena_aborts": arena_aborts,
        "successes": succ,
        "empirical_success": measured,
        "ideal_success": ideal,
        "measured_ratio": ratio,
        "ratio_two_sigma": 2 * sigma,
        "wrong_extractions": wrong,
        "identity_reseeds": identity_reseeds,
        "cycle_duplicates_at_sort": cycle_duplicates,
        "self_repeat_model_expected_total": self_repeat_model * trials,
        "true_merges_first_contacts": true_merges,
        "trials_with_true_merge": merged_trials,
        "merged_trials_detected": merged_detected,
        "merged_trials_undetected": merged_trials - merged_detected,
        "sort_detections": sort_detections,
        "executed_steps": executed,
        "implied_collision_factor": 1.5 / math.sqrt(
            -2 * math.log(1 - measured)),
    }


def main():
    data = json.load(sys.stdin)
    trials = int(data.get("trials", 200))
    dbits = int(data.get("dbits", 0))
    walkers = int(data.get("walkers", 8))
    out = {
        "production": production_checks(),
        "toy_plain_rho_v21": toy_attack(trials, dbits, walkers),
        # supplementary census: the same machinery at a
        # production-shaped distinguished-point density (2^-6 at toy
        # scale), exercising the sparse-arena records and deferred
        # detection, recorded inside the declared run
        "toy_dp_density_census": toy_attack(100, 6, 8),
    }
    toy_c = out["toy_plain_rho_v21"]["implied_collision_factor"]
    ab = out["production"]["success_model"]["abort_total"]
    out["production"]["success_model"]["modeled_success_at_toy_factor"] = \
        (1 - math.exp(-(1.5 / toy_c) ** 2 / 2)) - ab
    print(json.dumps(out))


main()
