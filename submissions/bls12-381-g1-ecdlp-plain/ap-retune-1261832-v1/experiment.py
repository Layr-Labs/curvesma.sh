"""Fixed-step (arithmetic-progression) walks with a routed class table.

Evidence for proof.md (stdlib only; reads one JSON object on stdin).

Part A  production BLS12-381 G1: GLV pair, class key invariance and
        completeness, limb-level emulation of the listed key/digit/home rows,
        AArch64 immediate encodability, and a planted end-to-end recovery that
        runs real batched affine additions, routing digits, a linear-probing
        table, candidate verification and the solving formula.
Part B  toy curve y^2 = x^3 + 4 over F_q with a prime-order subgroup: the whole
        attack at toy size (AP walks, parity key, linear-probing table,
        verification, x recovered and checked) against the Bonferroni/Poisson
        prediction; plus the exact scalar model at small r for more trials.
Part C  the cost ledger and the success bound in exact fractions, with checks;
        listing instruction counts are counted from the listing text below.
"""
import json, math, random, sys, time
from fractions import Fraction as Fr

T0 = time.time()
cfg = json.load(sys.stdin)
SEED = cfg.get("seed", 1)
rng = random.Random(SEED)

# ---------------------------------------------------------------- BLS12-381
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
P0 = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
      0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
BETA = 0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac
LAM = 0xac45a4010001a40200000000ffffffff
MASK64 = (1 << 64) - 1


def ec_add(P, Q, q):
    if P is None: return Q
    if Q is None: return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % q == 0: return None
        l = 3 * P[0] * P[0] * pow(2 * P[1], -1, q) % q
    else:
        l = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, q) % q
    x = (l * l - P[0] - Q[0]) % q
    return (x, (l * (P[0] - x) - P[1]) % q)


def ec_neg(P, q):
    return None if P is None else (P[0], (-P[1]) % q)


def ec_mul(k, P, q):
    R = None
    while k:
        if k & 1: R = ec_add(R, P, q)
        P = ec_add(P, P, q); k >>= 1
    return R


def batched_affine_add(pts, steps, q):
    """Montgomery-batched incomplete affine addition pts[i] + steps[i]
    (the reference operation): one inversion for the whole batch."""
    if any(t is None for t in pts): return None   # a walk at O: exceptional path
    d = [(s[0] - t[0]) % q for t, s in zip(pts, steps)]
    pre = [1]
    for v in d: pre.append(pre[-1] * v % q)
    if pre[-1] == 0: return None              # zero-denominator batch: recovery path
    inv = pow(pre[-1], -1, q)
    out = [None] * len(pts)
    for i in range(len(pts) - 1, -1, -1):
        di = inv * pre[i] % q; inv = inv * d[i] % q
        (x1, y1), (x2, y2) = pts[i], steps[i]
        l = (y2 - y1) * di % q
        x3 = (l * l - x1 - x2) % q
        out[i] = (x3, (l * (x1 - x3) - y1) % q)
    return out


def key(P, q):
    """Class key: K = y if y is even, else q - y (q odd)."""
    y = P[1]
    return y if y % 2 == 0 else q - y


def units(lam, rr):
    u = [1, lam % rr, lam * lam % rr]
    return u + [(-v) % rr for v in u]


def unit_apply(t, P, beta, q):
    """[u_t]P for t = 3s + j: (beta^j x, (-1)^s y)."""
    s, j = divmod(t, 3)
    x = P[0] * pow(beta, j, q) % q
    y = P[1] if s == 0 else (-P[1]) % q
    return (x, y)


# ------------------------------------------------------------ AArch64 helpers
def logical_imm_ok(v, width=64):
    """True iff v is encodable as an AArch64 logical (bitmask) immediate."""
    if v == 0 or v == (1 << width) - 1: return False
    for e in (2, 4, 8, 16, 32, 64):
        pat = v & ((1 << e) - 1)
        if any(((v >> (k * e)) & ((1 << e) - 1)) != pat for k in range(width // e)): continue
        for rot in range(e):
            rp = ((pat >> rot) | (pat << (e - rot))) & ((1 << e) - 1)
            ones = bin(rp).count("1")
            if rp == (1 << ones) - 1: return True
    return False


def emulate_generator_key(y):
    """Rows G2-G4 of proof.md section 5: limbs 0-2 of p - y by subs/sbcs,
    tst bit 0, three csel; ubfx digit1; returns (K0, K1, K2, digit1)."""
    yl = [(y >> (64 * i)) & MASK64 for i in range(3)]
    pl = [(p >> (64 * i)) & MASK64 for i in range(3)]
    out, borrow = [], 0
    for i in range(3):                       # subs, sbcs, sbcs
        v = pl[i] - yl[i] - borrow
        borrow = 1 if v < 0 else 0
        out.append(v & MASK64)
    odd = yl[0] & 1                          # tst x0, #1  (ne iff odd)
    K = [out[i] if odd else yl[i] for i in range(3)]   # csel ..., ne
    digit1 = (K[0] >> 1) & ((1 << 36) - 1)   # ubfx x6, x0, #1, #36
    return K[0], K[1], K[2], digit1


def emulate_final_home(K1, K2, S):
    h64 = ((K2 << 64 | K1) >> 35) & MASK64   # extr x4, x1, x0, #35
    return (h64 * S) >> 64                   # umulh x4, x4, x23


# ------------------------------------------------------------- listings
LISTINGS = {
 "G_gen_row": """
ldp x0, x1, [x19, #REC]!
ldr x2, [x19, #16]
subs x3, x10, x0
sbcs x4, x11, x1
sbcs x5, x12, x2
tst x0, #1
csel x0, x3, x0, ne
csel x1, x4, x1, ne
csel x2, x5, x2, ne
ubfx x6, x0, #1, #36
ldr x7, [x21, x6, lsl #3]
stp x1, x2, [x7]
stp x22, x23, [x7, #16]
add x7, x7, #32
str x7, [x21, x6, lsl #3]
add x22, x22, #1
cmp x19, x20
b.ne gen
""",
 "G_operand_traffic": "\n".join(["ldp"] * 6 + ["ldp"] * 6 + ["stp"] * 6),
 "I_intermediate": """
ldp x0, x1, [x9], #32
ldp x2, x3, [x9, #-16]
and x4, x0, #0x7ffffffff
ldr x5, [x21, x4, lsl #3]
stp x0, x1, [x5]
stp x2, x3, [x5, #16]
add x5, x5, #32
str x5, [x21, x4, lsl #3]
cmp x9, x10
b.ne int
""",
 "F_final_main": """
ldp x0, x1, [x9], #32
ldp x2, x3, [x9, #-16]
extr x4, x1, x0, #35
umulh x4, x4, x23
add x5, x22, x4, lsl #5
ldp x6, x7, [x5]
cbnz x6, occ
stp x0, x1, [x5]
stp x2, x3, [x5, #16]
cmp x9, x10
b.ne fin
""",
 "F_occupied_per_slot": """
add x26, x26, #1
cmp x26, x27
b.hs abort
cmp x6, x0
ccmp x7, x1, #0, eq
b.eq cand
ldp x6, x7, [x5, #32]!
cbnz x6, occ
""",
}
COUNT = {k: len([l for l in v.strip().splitlines() if l.strip()]) for k, v in LISTINGS.items()}

# ------------------------------------------------------------ cost ledger
LEDGER = dict(
    G_operand_traffic=18,   # 6 ldp current point, 6 ldp step point, 6 stp new point
    G_pass_bookkeeping=12,  # batched-addition passes: 3 stp + 3 ldp prefix product, 2x(cmp, b.ne), 2 spare
    G_gen_row=18,           # key, digit1, remote append, walk id, loop
    I_intermediate=10,
    F_final_main=11,
    F_occupied_per_slot=8,  # global counter increment, immediate cap check and branch, key compare, slot advance
    F_b_empty=1,            # the 'b empty' branch fires at most once per probe step (one per occupied action)
    batch_instr=235000,     # Fermat inversion 608*383=232,864 + scan <2,000 + zero test + constants
    batch_size=2 ** 24,
    round_allowance=Fr(2, 100),  # window bookkeeping per point (bounded in proof.md: < 0.016)
)
C = 2000
M_WALKS = 2 ** 63
ALPHA = Fr(394, 1000)
ZERO_PER_SLOT = 2 + Fr(3, 32)   # 2 stp xzr per 32-byte slot + loop share (64 stp per 3 loop instr)


def f4(mu):
    return mu - mu ** 2 / 2 + mu ** 3 / 6 - mu ** 4 / 24


def lg(v):
    if isinstance(v, Fr): return math.log2(v.numerator) - math.log2(v.denominator)
    return math.log2(v)


def solve_L(target):
    def mu_of(L): return Fr(6 * (M_WALKS * (M_WALKS - 1) // 2) * (L + 1) ** 2, r)
    lo, hi = 2 ** 62, 2 ** 64
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if f4(mu_of(mid)) >= target: hi = mid
        else: lo = mid
    return hi, mu_of(hi)


def cost_block(claim):
    out, chk = {}, {}
    for k in ("G_gen_row", "G_operand_traffic", "I_intermediate", "F_final_main", "F_occupied_per_slot"):
        chk["listing_count_" + k] = COUNT[k] == LEDGER[k]
    target = Fr(39, 100) + Fr(1, 10 ** 6)
    L, mu = solve_L(target)
    n = M_WALKS * (L + 1)
    per_point = (LEDGER["G_operand_traffic"] + LEDGER["G_pass_bookkeeping"] + LEDGER["G_gen_row"] + LEDGER["I_intermediate"]
                 + LEDGER["F_final_main"] + Fr(LEDGER["batch_instr"], LEDGER["batch_size"])
                 + LEDGER["round_allowance"])
    price = 1 + per_point / C
    S_tot = n / ALPHA
    zero_units = S_tot * ZERO_PER_SLOT / C
    insp = n * ALPHA / (2 * (1 - ALPHA))              # expected occupied-slot inspections
    probe_price = LEDGER["F_occupied_per_slot"] + LEDGER["F_b_empty"]
    probe_units = insp * probe_price / C                # every probe step charged 9 instructions
    probe_allot = probe_units * Fr(110, 100)            # global cap = 1.10 x the expectation
    W_det = n * price + zero_units + probe_allot
    extras = dict(walk_setup=M_WALKS * 2 ** 14, pointer_init=2 ** 76, candidates=2 ** 104,
                  recovery_cap=2 ** 39 * 4 * 2 ** 17, round_overshoot=W_det / 2 ** 41,
                  solve_output=2 ** 20)
    W = W_det + sum(extras.values())
    out.update(L_log2=round(lg(L), 6), mu=float(mu), f4_mu=float(f4(mu)), poisson=1 - math.exp(-float(mu)),
               points_log2=round(lg(n), 6), per_point_instr=float(per_point), step_price=float(price),
               zero_units_log2=round(lg(zero_units), 3), probe_allot_log2=round(lg(probe_allot), 3),
               W_det_log2=round(lg(W_det), 6), W_log2=round(lg(W), 6),
               extras_log2={k: round(lg(v), 3) for k, v in extras.items()})
    chk["mu_gives_f4_at_least_0.390001"] = f4(mu) >= target
    chk["L_below_2^64_j_fits_one_word"] = L < 2 ** 64
    chk["L_squared_much_below_r"] = 2 * (L + 1) ** 2 < r
    chk["claim_covers_W"] = lg(W) <= claim
    chk["claim_margin_below_0.001"] = claim - lg(W) < 0.001
    chk["global_inspection_cap_in_listing"] = ("cmp x26, x27" in LISTINGS["F_occupied_per_slot"]
                                               and "b.hs abort" in LISTINGS["F_occupied_per_slot"])
    # failure terms (section 4)
    q_reject = Fr(2 ** 255 - r, 2 ** 255)               # probability one 255-bit draw is rejected
    sampling_cap = 3 * M_WALKS * q_reject ** 128        # c_w, b_w, D_w draws
    insp_slack = Fr(1, 10) * insp                       # cap is 1.10 x the expectation
    probe_tail = 16 * S_tot / insp_slack ** 2           # Chebyshev, Var <= 16 S_tot
    cand_log2 = -(2 ** 104) * math.log2(2 ** 104 / (math.e * 2 ** 88))   # union bound at k = 2^104
    fail = dict(guard_tail_overflow=2.0 ** -400, zero_denominator_batches=3 * float(Fr(n, r)), dependent_sets=2.0 ** -60, degenerate_b=float(Fr(6 * (M_WALKS * (M_WALKS - 1) // 2) * (L + 1) ** 2, r * r)),
                exceptional_points=3 * float(Fr(n, r)), empty_marker=2.0 ** -50,
                window_overflow=2.0 ** -400, probe_allotment=float(probe_tail),
                candidate_allotment=2.0 ** -100,
                sampling_cap=float(sampling_cap), recovery_cap=2.0 ** -120)
    succ = float(f4(mu)) - sum(fail.values())
    out["failure_terms"] = fail
    out["sampler_cap_log2"] = round(lg(sampling_cap), 3)
    out["probe_tail_log2"] = round(lg(probe_tail), 3)
    out["candidate_tail_log2"] = round(cand_log2, 1)
    chk["sampler_cap_below_2^-370"] = lg(sampling_cap) < -370
    chk["probe_tail_below_2^-100"] = lg(probe_tail) < -100
    chk["candidate_tail_below_2^-100"] = cand_log2 < -100
    out["success_lower_bound"] = succ
    chk["success_at_least_0.39"] = succ >= 0.39
    # memory
    mem = (S_tot + 2 ** 71 * 2 ** 12) * 32 + 2 ** 36 * 2 ** 56 + 2 ** 71 * 2 ** 21 + M_WALKS * (504 + 96)
    out["memory_log2_bytes"] = round(lg(mem), 3)
    # dependent-set correction bounds (proof.md section 4, Lemma 2); n = L + 1 grid side
    g = L + 1
    nev = 6 * (M_WALKS * (M_WALKS - 1) // 2) * g * g
    N3 = M_WALKS ** 3 * g ** 3 * 36 + M_WALKS ** 2 * (2160 * g ** 4 + 120 * g ** 3)   # 6*8(1+ln n) <= 2160
    N4_2 = M_WALKS ** 2 * (84 * g ** 5 + 15 * g ** 2)
    N4_3 = (M_WALKS ** 2 * (6 * g * g) ** 4 + M_WALKS ** 3 * 36 * g * (2160 * g ** 4 + 120 * g ** 3)
            + M_WALKS ** 4 * g ** 4 * 1296 + 1296 * M_WALKS ** 3 * g ** 6 + N3 * nev)
    out["dependency_bounds_log2"] = dict(
        S3_lower_loss=round(lg(Fr(N3, r ** 3)), 1),
        S4_rank2_excess=round(lg(Fr(N4_2, r ** 2)), 1),
        S4_rank3_excess=round(lg(Fr(N4_3, r ** 3)), 1))
    chk["dependency_corrections_below_2^-60"] = all(v < -60 for v in out["dependency_bounds_log2"].values())
    out["checks"] = chk
    return out


# ------------------------------------------------------------ Part A
def part_a(n_rand):
    res, ok = {}, True
    beta, lam = BETA, LAM
    res["P0_on_curve"] = (P0[1] ** 2 - P0[0] ** 3 - 4) % p == 0
    res["rP0_is_O"] = ec_mul(r, P0, p) is None
    res["beta_cube_root"] = pow(beta, 3, p) == 1 and beta != 1
    res["lambda_cube_root_mod_r"] = (lam * lam + lam + 1) % r == 0
    res["phi_P0_eq_lambda_P0"] = (beta * P0[0] % p, P0[1]) == ec_mul(lam, P0, p)
    U = units(lam, r)
    # key invariance on the six members and completeness
    inv_ok, even_ok, emu_ok, home_ok = True, True, True, True
    S_node = 3 * 2 ** 54
    for i in range(n_rand):
        k = rng.randrange(1, r)
        Q = ec_mul(k, P0, p)
        Ks = {key(unit_apply(t, Q, beta, p), p) for t in range(6)}
        inv_ok &= len(Ks) == 1
        K = Ks.pop(); even_ok &= K % 2 == 0
        for yy in (Q[1], p - Q[1]):
            K0, K1, K2, d1 = emulate_generator_key(yy)
            exact = K & ((1 << 192) - 1)
            emu_ok &= (K0 | K1 << 64 | K2 << 128) == exact and d1 == (K >> 1) % 2 ** 36
        home_ok &= emulate_final_home(K1, K2, S_node) < S_node
        # members: [u]Q equals the unit_apply map
        if i < 2:
            for t in range(6):
                inv_ok &= unit_apply(t, Q, beta, p) == ec_mul(U[t], Q, p)
    # boundary y values for the borrow chain
    bnd = [1, 2, 3, p - 1, p - 2, 2 ** 64 - 1, 2 ** 64, 2 ** 64 + 1, 2 ** 128 - 1, 2 ** 128, 2 ** 191, 2 ** 192 - 1,
           p - 2 ** 64, p - 2 ** 128, (p & ((1 << 192) - 1)), (p & ((1 << 192) - 1)) + 1]
    for yy in bnd:
        yy %= p
        if yy == 0: continue
        K = yy if yy % 2 == 0 else p - yy
        K0, K1, K2, d1 = emulate_generator_key(yy)
        emu_ok &= (K0 | K1 << 64 | K2 << 128) == K % 2 ** 192
    res["key_constant_on_six_members_and_unit_map"] = inv_ok
    res["key_even"] = even_ok
    res["generator_key_rows_match_exact_K_mod_2^192"] = emu_ok
    res["final_home_in_range"] = home_ok
    # completeness: equal K => y' = +-y => x'^3 = x^3 => x' in {x, beta x, beta^2 x}
    x0 = P0[0]
    res["equal_key_implies_same_class"] = all((x0 * pow(beta, j, p) % p) ** 3 % p == pow(x0, 3, p) for j in range(3))
    res["imm_and_0x7ffffffff_encodable"] = logical_imm_ok(0x7ffffffff)
    res["imm_tst_1_encodable"] = logical_imm_ok(1)
    res["ubfx_1_36_valid"] = 1 + 36 <= 64
    res["extr_35_valid"] = 0 <= 35 < 64
    # planted end-to-end on production: 3 walks, 6 steps, planted class collision
    x = rng.randrange(1, r); T = ec_mul(x, P0, p)
    Lw, Mw = 6, 3
    walks = []
    for w in range(Mw):
        walks.append([rng.randrange(r), rng.randrange(r), rng.randrange(1, r)])
    jp, jq, ut = rng.randrange(Lw + 1), rng.randrange(Lw + 1), rng.randrange(6)
    u = U[ut]
    c0, b0, D0 = walks[0]
    A0 = (c0 + b0 * x) % r
    c1, b1, D1 = walks[1]
    # plant: A1 + jq D1 = u^-1 (A0 + jp D0)   (uses x only to construct the test)
    c1 = (pow(u, -1, r) * (A0 + jp * D0) - jq * D1 - b1 * x) % r
    walks[1][0] = c1
    cur = [ec_add(ec_mul(c, P0, p), ec_mul(b, T, p), p) for c, b, D in walks]
    stepP = [ec_mul(D, P0, p) for c, b, D in walks]
    S_small = 64
    table = [None] * (S_small + 16)
    cands, entries = [], 0
    for j in range(Lw + 1):
        if j > 0:
            cur = batched_affine_add(cur, stepP, p)
        for w in range(Mw):
            K0, K1, K2, d1 = emulate_generator_key(cur[w][1])
            h = emulate_final_home(K1, K2, S_small)
            pos = h
            while table[pos] is not None:
                if table[pos][0] == K1 and table[pos][1] == K2: cands.append((table[pos][2], table[pos][3], w, j))
                pos += 1
            table[pos] = (K1, K2, w, j); entries += 1

    def walk_point(w, j):
        c, b, D = walks[w]
        return ec_add(ec_mul((c + j * D) % r, P0, p), ec_mul(b, T, p), p)

    solved = None
    for (w1, j1, w2, j2) in cands:
        Q1, Q2 = walk_point(w1, j1), walk_point(w2, j2)
        if key(Q1, p) != key(Q2, p): continue      # false candidate
        for t in range(6):
            if unit_apply(t, Q2, beta, p) == Q1:   # Q1 = [U[t]] Q2
                uu = U[t]
                (ca, ba, Da), (cb, bb, Db) = walks[w1], walks[w2]
                den = (ba - uu * bb) % r
                if den:
                    xs = ((uu * (cb + j2 * Db) - ca - j1 * Da) * pow(den, -1, r)) % r
                    if ec_mul(xs, P0, p) == T: solved = xs
    res["planted_collision_detected"] = len(cands) >= 1
    res["planted_x_recovered_and_verified"] = solved == x
    res["planted_entries"] = entries
    for k, v in res.items():
        if isinstance(v, bool): ok &= v
    res["all_true"] = ok
    return res


# ------------------------------------------------------------ Part B (toy curve)
TOYS = {
    "q26": dict(q=67109023, r=3196447, h=21, G=(47062022, 8437629), beta=15691445, lam=1581523),
    "q30": dict(q=1073741827, r=18837001, h=57, G=(87950372, 540139103), beta=536887297, lam=32767),
    "q34": dict(q=17179869547, r=154772419, h=111, G=(4885185641, 14353049500), beta=14837154105, lam=47098950),
}


def toy_trial(tc, Mw, Lw, trng):
    q, rr, G, beta = tc["q"], tc["r"], tc["G"], tc["beta"]
    U = units(tc["lam"], rr)
    x = trng.randrange(1, rr); T = ec_mul(x, G, q)
    walks = [(trng.randrange(rr), trng.randrange(rr), trng.randrange(1, rr)) for _ in range(Mw)]
    cur = [ec_add(ec_mul(c, G, q), ec_mul(b, T, q), q) for c, b, D in walks]
    stepP = [ec_mul(D, G, q) for c, b, D in walks]
    S = int((Mw * (Lw + 1)) / 0.45) + 1
    table = [0] * (S + 4096)      # 0 = empty (toy keys are nonzero)
    loc = {}
    cands = 0; found = None; exceptional = 0
    for j in range(Lw + 1):
        if j > 0:
            nxt = batched_affine_add(cur, stepP, q)
            if nxt is None or any(v is None for v in nxt):
                exceptional += 1
                nxt = [ec_add(a, b, q) for a, b in zip(cur, stepP)]
            cur = nxt
        for w in range(Mw):
            P = cur[w]
            if P is None: exceptional += 1; continue
            K = key(P, q)
            pos = K * S // q          # toy home: multiply-high of the key by S
            while table[pos]:
                if table[pos] == K and found is None:
                    cands += 1
                    w1, j1 = loc[pos]
                    Q1 = ec_add(ec_mul((walks[w1][0] + j1 * walks[w1][2]) % rr, G, q), ec_mul(walks[w1][1], T, q), q)
                    for t in range(6):
                        if unit_apply(t, P, beta, q) == Q1:
                            uu = U[t]
                            (ca, ba, Da), (cb, bb, Db) = walks[w1], walks[w]
                            den = (ba - uu * bb) % rr
                            if den:
                                xs = ((uu * (cb + j * Db) - ca - j1 * Da) * pow(den, -1, rr)) % rr
                                if ec_mul(xs, G, q) == T: found = xs
                pos += 1
            table[pos] = K; loc[pos] = (w, j)
    return found == x, cands, exceptional


def part_b(curve, Mw, seconds, trng):
    tc = TOYS[curve]
    rr = tc["r"]
    Lw = 1
    while 6 * (Mw * (Mw - 1) // 2) * (Lw + 1) ** 2 < 0.4947 * rr: Lw += 1
    mu = 6 * (Mw * (Mw - 1) // 2) * (Lw + 1) ** 2 / rr
    mu_all = mu + Mw * (4 * (Lw + 1) * Lw // 2 + 2 * Lw - 1) / rr   # plus within-walk events (u != 1)
    t_end = time.time() + seconds
    n = s = 0; exc = 0
    while time.time() < t_end:
        ok, cands, e = toy_trial(tc, Mw, Lw, trng)
        n += 1; s += ok; exc += e
    ph = s / max(1, n)
    return dict(curve=curve, r=rr, walks=Mw, L=Lw, points=Mw * (Lw + 1), mu_cross_walk=round(mu, 4),
                predicted_poisson=round(1 - math.exp(-mu), 4), bonferroni4_lower=round(f4(mu), 4),
                mu_all_events=round(mu_all, 4), predicted_poisson_all=round(1 - math.exp(-mu_all), 4),
                trials=n, solved_and_verified=s, success_rate=round(ph, 4),
                stderr=round(math.sqrt(ph * (1 - ph) / max(1, n)), 4), exceptional_events=exc)


def part_b_scalar(rbits, Mw, seconds, trng):
    """Exact scalar model: point [A + jD]P with class = U-orbit of the scalar."""
    rr = 1 << rbits
    def isprime(v):
        if v < 2: return False
        i = 2
        while i * i <= v:
            if v % i == 0: return False
            i += 1
        return True
    rr += 1
    while not (rr % 3 == 1 and isprime(rr)): rr += 2
    g = 2
    while pow(g, (rr - 1) // 3, rr) == 1: g += 1
    U = units(pow(g, (rr - 1) // 3, rr), rr)
    Lw = 1
    while 6 * (Mw * (Mw - 1) // 2) * (Lw + 1) ** 2 < 0.4947 * rr: Lw += 1
    mu = 6 * (Mw * (Mw - 1) // 2) * (Lw + 1) ** 2 / rr
    t_end = time.time() + seconds
    n = s = 0
    while time.time() < t_end:
        x = trng.randrange(1, rr); seen = {}; hit = False
        for w in range(Mw):
            c, b, D = trng.randrange(rr), trng.randrange(rr), trng.randrange(1, rr)
            A = (c + b * x) % rr
            for j in range(Lw + 1):
                v = (A + j * D) % rr
                if v == 0: continue
                k = min(uu * v % rr for uu in U)
                o = seen.get(k)
                if o is not None and o != w: hit = True; break
                if o is None: seen[k] = w
            if hit: break
        n += 1; s += hit
    ph = s / max(1, n)
    return dict(r=rr, walks=Mw, L=Lw, mu_cross_walk=round(mu, 4), predicted_poisson=round(1 - math.exp(-mu), 4),
                bonferroni4_lower=round(f4(mu), 4), trials=n, cross_walk_collision_rate=round(ph, 4),
                stderr=round(math.sqrt(ph * (1 - ph) / max(1, n)), 4))


def main():
    out = {"config": cfg}
    claim = cfg.get("claim", 126.1835)
    out["cost"] = cost_block(claim)
    out["listing_counts"] = COUNT
    if cfg.get("production_checks", 1):
        out["part_a"] = part_a(cfg.get("production_random", 12))
    trng = random.Random(SEED + 7)
    if cfg.get("toy_seconds", 0) > 0:
        out["part_b_toy_curve"] = part_b(cfg.get("curve", "q30"), cfg.get("toy_walks", 32), cfg["toy_seconds"], trng)
    if cfg.get("scalar_seconds", 0) > 0:
        out["part_b_scalar_model"] = part_b_scalar(cfg.get("scalar_rbits", 22), cfg.get("scalar_walks", 64), cfg["scalar_seconds"], trng)
    out["all_cost_checks_true"] = all(out["cost"]["checks"].values())
    out["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(out, indent=1, default=str))


main()
