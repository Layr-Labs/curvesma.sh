"""Evidence experiment for the plain-rho-v38 submission (secp256k1-ecdlp).

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


def fermat_chain_instr(p, mult_instr):
    """Instruction count of the FIXED binary left-to-right chain
    for z^(p-2) plus affine normalization. The chain is
    constructive and VERIFIED by verify_inversion_chain below:
    (bitlen(p-2)-1) squarings + (popcount(p-2)-1) multiplies,
    then 4 normalization multiplies (u = z^-1; x = X*u^2,
    y = Y*u^2*u), every multiplication charged mult_instr
    (the pinned per-multiplication price), 20 control/addressing
    instructions per chain step, 8 transfers."""
    e = p - 2
    steps = (e.bit_length() - 1) + (e.bit_count() - 1)
    return (steps + 4) * mult_instr + steps * 20 + 8


def verify_inversion_chain():
    """Execute the priced chain and verify it computes a field
    inverse. This is the exact construction the ledger prices:
    left-to-right binary exponentiation of z^(p-2) starting at
    z (the MSB is consumed by the start)."""
    p = P_FIELD
    e = p - 2
    z = 0x123456789ABCDEF1
    assert 0 < z < p
    r_ = z
    sq = 0
    mu = 0
    for i in range(e.bit_length() - 2, -1, -1):
        r_ = (r_ * r_) % p
        sq += 1
        if (e >> i) & 1:
            r_ = (r_ * z) % p
            mu += 1
    assert (r_ * z) % p == 1, "chain must compute the inverse"
    assert r_ == pow(z, p - 2, p)
    return {"squarings": sq, "multiplies": mu, "verified": True}


def verify_canonicalization():
    """Execute and verify the charged canonicalization construction:
    ANY value inside the multiplication calibration's output
    container (5x52 = 260 bits) is reduced to the canonical
    residue in [0,p) by the FOLD construction: repack V as
    L + H*2^256 (H < 2^4), fold with the fixed constant
    K = 2^256 mod p, carry-fold (at most 2 iterations; after
    them W < 2^256 + 2^33 < 2p by the explicit inequality
    2p = 2^257 - 2^33 - 1954), then conditional subtracts.
    The construction assumes NOTHING about the calibration's
    output beyond the container size."""
    p = P_FIELD
    base = 256 if p.bit_length() <= 257 else 384
    folds = 2 if base == 256 else 12
    subtracts = 4 if base == 256 else 10
    container = base + (4 if base == 256 else 32)
    K = (1 << base) % p
    rng = random.Random(0xC0FFEE)
    for _ in range(100):
        V = rng.randrange(0, 1 << container)
        W = V
        for _ in range(folds):
            hi = W >> base
            W = (W & ((1 << base) - 1)) + hi * K
        assert (W - V % p) % p == 0
        for _ in range(subtracts):
            if W >= p:
                W -= p
        assert W < p and W == V % p
    return {"verified": True, "base": base, "folds": folds,
            "container_bits": container, "constant_K": K,
            "trials": 100}


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
    # every main-loop step is a genuine incomplete affine
    # addition. WHAT THE 6M REFERENCE CHARGE COVERS (identified
    # exactly): 3 multiplication instructions for the addition
    # arithmetic (slope lambda = (y2-y1)*(x2-x1)^-1 after the
    # batch inversion; x3 = lambda^2 - x1 - x2; y3 =
    # lambda*(x1-x3) - y1) plus 3 multiplication instructions
    # for the walker's share of the Montgomery batch passes
    # (two forward sweep multiplications and one backward sweep
    # multiplication per step). NOTHING ELSE is inside the
    # reference charge. The exceptional states S = +A_j and
    # S = -A_j both have a ZERO slope denominator, are detected
    # by the full-width zero test, and RESEED - a doubling never
    # executes in the main loop. The stored state is made
    # CANONICAL by an explicit charged construction each step:
    # the multiplication calibration (179/383 instructions)
    # computes in 5x52-bit limbs, so each output coordinate is
    # repacked to the stored 4/6x64-bit limb form and once
    # conditionally reduced against p - see the
    # canonicalization_repack schedule entry.
    FIELD_MULT_INSTR = 179  # pinned instructions per field
    # multiplication on this track (BLS clone pins 383)
    # canonicalization per step, by the EXECUTED AND VERIFIED
    # fold construction (see verify_canonicalization; the
    # construction assumes NOTHING about the multiplication
    # calibration's output beyond its container size):
    # any value inside the 260-bit container is repacked as
    # L + H*2^256 with H < 2^4 (carry-propagating repack over
    # all five limbs, 32 instructions), folded with the fixed
    # constant K = 2^256 mod p = 2^32 + 977 = 4294968273
    # (one word multiply
    # 16 + add-with-carry 16), carry-folded at most twice
    # (2 x 24: test the bit above 2^256 + conditional add of
    # K; after two folds W < 2^256 + 2^33 < 2p by the explicit
    # inequality), then up to four conditional subtracts
    # (4 x 12: 4-limb compare + conditional 4-limb subtract)
    # - 160 instructions per value. FOUR values per step are
    # normalized conservatively (both stored coordinates plus
    # the two intermediate products that feed them):
    canon_check = verify_canonicalization()
    assert canon_check["verified"]
    CANON_INSTR = 4 * 160
    LOAD_STORE_INSTR = 64  # 16 64-bit state words x 4
    # instructions each (load or store, generous); the BLS
    # clone charges 96 (24 words)
    ZERO_TEST_INSTR = 4  # full-width zero test: 3 limb ORs +
    # branch (4 limbs); BLS charges 8 for 6 limbs
    # Fermat inversion + affine normalization, COMPUTED from the
    # fixed binary chain (see fermat_chain_instr / the verified
    # inversion_chain output): 255 squarings + 248 multiplies +
    # 4 normalization multiplies on this track, 20 control
    # instructions per chain step, 8 transfers
    NORM_INV_INSTR = fermat_chain_instr(P_FIELD, FIELD_MULT_INSTR)
    inv_chain = verify_inversion_chain()
    assert inv_chain["verified"]
    # the itemized worst-case per-step schedule (every
    # instruction outside the 6M reference charge):
    psched = {
        "state_loads_stores": LOAD_STORE_INSTR,
        # x,y in/out: 16 words x 4 instructions each
        "index_extract": 2,      # j = x low word & 31
        "rep_add": 32,           # (a,c_j) mod r: limb adds +
        # conditional subtract (generous at any limb width)
        "dp_predicate": 6,       # low-80-bit zero test + branch
        "denominator_test": ZERO_TEST_INSTR,
        "montgomery_pass": 24,   # per-step share of the batch
        # sweeps: 3 passes x (2 loads + 1 multiply-linked
        # store + pointer arithmetic) - generous
        "round_robin": 8,       # walker index, loop control,
        # batch trigger
        "reseed_branch": 6,     # not-taken branch bookkeeping
        "canonicalization_repack": CANON_INSTR,
        # repack both coordinates from the multiplication's
        # 5x52-bit internal limbs to canonical 4x64-bit stored
        # limbs + one conditional reduction against p each
    }
    other = sum(psched.values())
    per_step = 1 + other / C_REF
    N = r - 1
    t = 1.5 * math.sqrt(N)
    per_walker = t / 2 ** 16
    # flat component charges (instructions)
    DBL_FLAT_INSTR = 1600
    DRBG_FLAT_INSTR = 4096
    JADD_INSTR = 3401  # Jacobian-Jacobian add, ~19 mults
    JDBL_INSTR = 1432  # Jacobian doubling, ~8 mults
    # per-op unit ceilings used by the table build and the
    # reseed price (all unbatched Jacobian / Fermat prices)
    jadd_unit = -(-JADD_INSTR // C_REF)
    jdbl_unit = -(-JDBL_INSTR // C_REF)
    inv_unit = -(-NORM_INV_INSTR // C_REF)
    # seeding, charged with GENEROUS FLAT upper bounds so no
    # implementation-level count can exceed them. Per walker:
    # [a0]P by fixed-base comb (32 adds + 32 doublings), [b0]T by
    # 4-bit windowed double-and-add (256 doublings + 64
    # additions), final addition: 97 additions and 288
    # doublings. Every ADDITION is charged 2 reference operations
    # (the reference operation is 1; 2 covers any full-width
    # modular add with carries, conditional reduction, correction
    # and transfers). Every DOUBLING is charged 1600/3000
    # instructions flat (secp/BLS; an affine doubling is 7 field
    # multiplications = 1253/2681 instructions, so the flat
    # allowance covers the adds with their carry chains,
    # conditional reductions, corrections and transfers at any
    # reasonable count). Randomness: TWO independent 512-bit
    # ChaCha20 outputs per walker, charged 4096 instructions per
    # output flat (covers 80 quarter-rounds = 960 instructions
    # plus initialization, feed-forward, state transfers,
    # counters and the mod-r reduction). Exact component sum:
    #   97*2 + 288*1600/C_REF + 2*4096/C_REF
    #     = 662.992 (secp) / 630.096 (BLS) reference operations.
    # The INITIAL seeding batch-shares the Montgomery inversion
    # across the 2^16 INDEPENDENT walkers (a batch of
    # independent scalar multiplications), which is what the
    # batched affine prices assume. The per-walker price is
    # DERIVED and ASSERTED here, not merely assigned:
    add_unit = 2
    dbl_unit = -(-DBL_FLAT_INSTR // C_REF)      # ceil to whole ops
    drbg_unit = -(-2 * DRBG_FLAT_INSTR // C_REF)
    per_walker_computed = 97 * add_unit + 288 * dbl_unit + drbg_unit
    per_walker_seed = -(-per_walker_computed // 100) * 100
    assert per_walker_seed >= per_walker_computed
    # seeding inversions: the scalar multiplications have 385
    # SEQUENTIAL layers per walker (97 additions + 288 doublings;
    # each layer's affine output feeds the next), and each layer
    # is one Montgomery batch across the 2^16 INDEPENDENT walkers:
    # 385 layer-batches, one Fermat inversion each (the
    # batch-multiplication share is already inside the
    # per-operation prices)
    seeding = (2 ** 16 * per_walker_seed
               + 385 * NORM_INV_INSTR / C_REF)
    # FIXED-BASE COMB TABLE for [a_0]P, complete construction,
    # charged with NO batching anywhere: 224 setup doublings
    # (8 teeth at 32-bit spacing: [2^32]P ... [2^224]P via
    # sequential doublings), 247 subset additions (2^8 - 1
    # nonzero subsets minus the 8 single-tooth bases), and 255
    # affine normalizations (8 bases + 247 subsets) - every
    # operation an unbatched Jacobian operation and every
    # affine result normalized by its OWN charged Fermat
    # inversion; plus the 15-entry window table for [b_0]T:
    # 14 additions + 15 normalizations. COMPUTED:
    table_build = (224 * jdbl_unit + 247 * jadd_unit
                   + 255 * inv_unit
                   + 14 * jadd_unit + 15 * inv_unit)
    # constant construction: 32 fixed-base scalar multiples
    # (~256 reference additions each) plus their sha512
    # derivation and mod-r reduction, charged generously as
    # 2^14 reference operations flat
    constants_build = 2 ** 14
    # distinguished-point density: the predicate tests the low 80
    # bits of the field x-COORDINATE (not a scalar mod r). Under
    # H1 each visited state is a uniform subgroup point. At most
    # TWO curve points share one admissible x (the +/- y pair),
    # and at most floor(P_FIELD/2^80)+1 <= P_FIELD/2^80 + 1 field
    # values pass the predicate, so at most
    # 2*(P_FIELD/2^80 + 1) subgroup points qualify out of r-1;
    # on secp256k1 (P_FIELD ~ r) this conditional-on-H1 bound is
    # 2^-79 + 2^-253, on BLS the full-field count is vacuous and
    # the 2^-80 figure is the standard equidistribution heuristic
    # H2 (subgroup x-coordinates equidistributed; checked
    # empirically by the toy census below). A CONSERVATIVE
    # uniform bound, charged on BOTH tracks:
    #   density <= 2^-79 + 2^-160
    # (the secp rigorous figure and the BLS heuristic rounded up
    # by a full factor of 2). Dependence between steps cannot
    # change the per-state marginal, so Markov applies to this
    # mean bound with no independence assumption.
    inserts = t * (2 ** -79 + 2 ** -160)
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
    # linear scan at 64 instructions per record flat.
    final_sort = (2 * 2 ** cap_bits * cap_bits * 48
                  + 2 ** cap_bits * 64) / C_REF    # extraction (modular inverse ~700 instr) + verification of
    # [x]P == T (one ~256-add scalar multiple ~ 2^8 reference ops)
    extract_verify = 2 ** 9
    # special-state reseeds (S = +/-A_j, both detected by the
    # zero-denominator branch): rate at most 64/(r-1) per step
    # under H1, capped at 2^20, abort beyond. Reseed WORK is a
    # COMPLETE UNBATCHED construction with no batching
    # assumption: both scalar multiplications in Jacobian
    # coordinates (mixed additions ~11 mults, Jacobian-Jacobian
    # combine ~19 mults - ALL 97 additions charged at the
    # Jacobian-Jacobian price ~3401/7277 instructions; Jacobian
    # doublings ~8 mults ~1432/3064 instructions), ONE Fermat
    # normalization inversion making the affine state explicit
    # (~271 multiplications, 48509/103793 instructions) because
    # the walk needs x = X/Z^2 for j = x mod 32 and the
    # predicate, and TWO fresh ChaCha20 outputs; DERIVED and
    # ASSERTED:
    drbg_unit = -(-2 * DRBG_FLAT_INSTR // C_REF)
    reseed_computed = (97 * jadd_unit + 288 * jdbl_unit
                       + inv_unit + drbg_unit)
    reseed_price = -(-reseed_computed // 100) * 100
    assert reseed_price >= reseed_computed
    identity_reseed_cap_work = 2 ** 20 * reseed_price
    # recurring MAIN-LOOP batch inversions: the affine walk
    # completes one round-robin round of 2^16 walker steps per
    # Montgomery batch (the next step of a walker needs its
    # affine x immediately, but the 2^16 walkers are
    # independent), so the attack executes ceil(t / 2^16)
    # batches, each with ONE Fermat inversion:
    batch_inversions = math.ceil(t / 2 ** 16) * NORM_INV_INSTR / C_REF
    W = (per_step * t + seeding + table_build + constants_build
         + dp_appends_worst
         + final_sort + extract_verify + identity_reseed_cap_work
         + batch_inversions)
    tl2 = math.log2(W)
    num = math.ceil(tl2 * 10 ** 6)
    fr = Fraction(num, 10 ** 6)
    out["ledger"] = {
        "per_step_other_instructions": other,
        "per_step_schedule": psched,
        "per_step_units": per_step,
        "batch_inversions_units_log2": math.log2(batch_inversions),
        "budget_multiplier": 1.5,
        "t_log2": math.log2(t),
        "per_walker_log2": math.log2(per_walker),
        "seeding_units_log2": math.log2(seeding + table_build),
        "canonicalization_check": canon_check,
        "inversion_chain": {"squarings": inv_chain["squarings"], "multiplies": inv_chain["multiplies"], "instructions": NORM_INV_INSTR, "verified": True},
        "per_walker_seed": per_walker_seed,
        "per_walker_computed": per_walker_computed,
        "reseed_price": reseed_price,
        "reseed_computed": reseed_computed,
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
    # special-state reseed abort (S = +/-A_j, BOTH exceptional
    # states, rate 64/(r-1) per step under H1): Markov on the
    # expected count; the computed bound is printed and quoted
    # as-is (never rounded below the computed value)
    identity_expected = t * 64 / (r - 1)
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
# Part B: toy-scale complete v38 attack (plain additive walk)
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
                if C2 is None or C == consts[j]:
                    # exceptional states S = -A_j (identity: zero
                    # denominator, C2 is None) and S = +A_j (the
                    # doubling case, also a zero denominator) BOTH
                    # reseed in production - a doubling never
                    # executes in the main loop. Rate 64/(r-1) per
                    # step under H1 (capped at 2^20, abort beyond);
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
        "toy_plain_rho_v38": toy_attack(trials, dbits, walkers),
        # supplementary census: the same machinery at a
        # production-shaped distinguished-point density (2^-6 at toy
        # scale), exercising the sparse-arena records and deferred
        # detection, recorded inside the declared run
        "toy_dp_density_census": toy_attack(100, 6, 8),
    }
    toy_c = out["toy_plain_rho_v38"]["implied_collision_factor"]
    ab = out["production"]["success_model"]["abort_total"]
    out["production"]["success_model"]["modeled_success_at_toy_factor"] = \
        (1 - math.exp(-(1.5 / toy_c) ** 2 / 2)) - ab
    print(json.dumps(out))


main()
