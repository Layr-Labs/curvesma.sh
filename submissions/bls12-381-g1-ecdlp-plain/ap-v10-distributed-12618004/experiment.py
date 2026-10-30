"""Fixed-step (arithmetic-progression) walks with a distributed class table.

Evidence for proof.md (stdlib only; reads one JSON object on stdin).

Part A  production BLS12-381 G1: GLV pair, class key invariance and
        completeness, limb-level emulation of the listed key/digit rows,
        AArch64 immediate encodability, and a planted end-to-end recovery that
        runs real batched affine additions, a split generator/message/network/
        receive/server pipeline (all addresses 64-bit local), a chained table,
        candidate verification and the solving formula.
Part A2 production-bit heuristic-H' consistency evidence: real G1 AP walks,
        the exact production extraction code, bit-balance z-scores over the 192
        key and 71 routing bits, pairwise bit chi-square, small-table Poisson
        occupancy, and truncated 164-bit fingerprint collisions on high/middle/
        low bit slices.  This is consistency evidence on a scaled sample, not a
        proof of H'.
Part B  toy curve y^2 = x^3 + 4 over F_q with a prime-order subgroup: the whole
        attack at toy size (AP walks, parity key, chained table,
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
    h64 = K2 & MASK64                        # h64 = limb 2 of the class key
    return (h64 * S) >> 64                   # umulh x8, x2, x20


# ------------------------------------------------------------- listings
LISTINGS = {
 "GEN_key": """
ldp x0, x1, [x19, #REC]!
ldp x2, x22, [x19, #16]
subs x3, x10, x0
sbcs x4, x11, x1
sbcs x5, x12, x2
tst x0, #1
csel x0, x3, x0, ne
csel x1, x4, x1, ne
csel x2, x5, x2, ne
ubfx x6, x0, #1, #36
ubfx x7, x1, #0, #35
""",
 "MSG": """
stp x1, x2, [x14]
stp x22, x23, [x14, #16]
stp x6, x7, [x14, #32]
add x14, x14, #56
""",
 "NET": """
ldr x9, [x13, #32]
ldr x11, [x12, x9, lsl #3]
br x11
ldr x9, [x13, #40]
ldr x11, [x12, x9, lsl #3]
br x11
""",
 "RCV": """
ldp x1, x2, [x9]
ldp x22, x23, [x9, #16]
add x9, x9, #56
cmp x9, x13
b.ne next
mov x9, x10
""",
 "SRV_insert": """
umulh x8, x2, x20
add x8, x25, x8, lsl #3
ldp x15, x13, [x25]
add x16, x15, #40
cmp x16, x13
b.hs abort
ldr x9, [x8, #32]
str x9, [x15]
str x15, [x8, #32]
stp x1, x2, [x15, #8]
stp x22, x23, [x15, #24]
str x16, [x25]
cbz x9, next
""",
 "F_occupied_per_slot": """
ldp x3, x4, [x9, #8]
adds x26, x26, #1
tst x26, #0xffff
b.ne scanc
ldadd x27, [x30], x31
mov x26, xzr
cmp x27, x29
b.lo scanc
ldadd x28, [x32], x33
cmp x28, x18
b.hs abort
cmp x3, x1
ccmp x4, x2, #0, eq
b.eq cand
ldr x9, [x9]
cbnz x9, scan
""",
 "LOOP_CONTROL": """
subs x17, x17, #1
b.ne group
""",
 "CAND_COUNTER": """
ldadd x29, [x30], x31
cmp x29, x32
b.hi abort
ldadd x27, [x33], x31
cmp x27, x18
b.hs abort
""",
}
COUNT = {k: len([l for l in v.strip().splitlines() if l.strip()]) for k, v in LISTINGS.items()}

# ------------------------------------------------------------ cost ledger
UNROLL = 128
LEDGER = dict(
    GEN_operand_traffic=18,  # 6 ldp current point, 6 ldp step point, 6 stp new point (reference interface)
    GEN_bookkeeping=12,      # batched-addition passes: 3 stp + 3 ldp prefix product, 2x(cmp,b.ne), 2 spare
    GEN_key=11,              # 2 record loads, subs/sbcs/sbcs, tst, 3 csel, ubfx digit1, ubfx digit2
    MSG=4,                   # 2 stp payload, 1 stp descriptor, 1 post-index enqueue
    NET=6,                   # 2 levels x (load destination, load next hop, branch)
    RCV=6,                   # 2 ldp payload, load next pointer, load/add/store queue head
    SRV_insert=13,           # umulh, head address, header ldp, bump add, cmp, b.hs, head ldr,
                             # 2 links, 2 stp payload, bump str, cbz (loop control separate)
    F_occupied_per_slot=9,   # common path: head load, local counter, block test, key compare, chain step
    F_block=4,               # per 2^16 comparisons: global low-word ldadd, local reset, carry test, branch
    F_carry=3,               # per 2^64 comparisons: high-word ldadd, cap compare, abort branch
    CAND_COUNTER=6,          # shared 128-bit candidate counter: two ldadd's plus two cap checks
    LOOP_CONTROL=2,          # subs/b.ne once per UNROLL points, on the generator and on the server
    batch_instr=235000,      # Fermat inversion 608*383=232,864 + scan <2,000 + zero test + constants
    batch_size=2 ** 26,      # 2^37 generator nodes x 2^26 walk records each
    round_allowance=Fr(2, 1000),  # per-point batch dispatch / queue tail publish allowance
)
BLK = 2 ** 16                 # inspections per shared-counter update (mask 0xffff is encodable)
C = 2000
M_WALKS = 2 ** 63
N_GEN = 2 ** 37           # generator nodes; each holds 2^26 walk records locally
S_TOT = 2 ** 128          # 2^71 table servers x 2^57 head slots
N_NODES = 2 ** 71
S_NODE = 2 ** 57
ZERO_PER_HEAD = Fr(500244, 10 ** 6)  # 8-byte heads: 1 stp per 2 heads + 3/64 loop share
ARENA_MARGIN = Fr(105, 100)           # per-server entry arena capacity = 1.05 x the mean


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


def written_regs(text):
    """Destination registers written by a bare listing (labels/comments absent)."""
    out = set()
    no_dest = {"stp", "str", "cmp", "tst", "ccmp", "b.hs", "b.eq", "b.ne",
               "cbz", "cbnz", "dsb", "br", "b.cc", "b.lo", "b.hi"}
    dst1 = {"ldr", "subs", "sbcs", "csel", "ubfx", "extr", "umulh", "add",
            "adc", "adds", "mov", "ldadd"}
    for line in text.strip().splitlines():
        line = line.split("//")[0].strip()
        if not line or line.endswith(":"): continue
        parts = line.replace(",", " ").split()
        op = parts[0]
        if op in no_dest:
            continue
        if op == "ldp":
            out.update(parts[1:3])
        elif op in dst1:
            out.add(parts[1])
        else:
            raise ValueError("unparsed listing opcode: " + op)
    return out


def emulate_listing_GF(points, S, cap_entries):
    """Execute the v10 distributed data flow on one table-server bank.

    Server header at x25: bump pointer at offset 0, arena end at offset 8,
    receive-queue head at offset 16, S head cells from offset 32.  Entries are
    40 bytes: next at 0, K1 at 8, K2 at 16, walk id at 24, j at 32.  x8 = bank
    base plus the home offset (never written by the scan); x9 is the old head,
    then the traversal pointer; x10-x12 are the generator's p limbs, set once
    and never written; x15/x16 are the bump pointer and the new bump; x13 is
    the arena end fetched per point by ldp x15,x13,[x25].  The generator and
    the server are separate nodes; only the 32-byte key/locator payload (plus
    the 16-byte routing descriptor) crosses the network.  Returns the head
    cells, arena, candidates and a data-flow consistency flag.
    """
    pl = [(p >> (64 * i)) & MASK64 for i in range(3)]
    heads = [0] * S
    arena = []
    B = 32 + 8 * S                 # arena starts after the server header and heads
    header = [B, B + cap_entries * 40]
    cands, ok = [], True
    for (yy, w, j) in points:
        K0, K1, K2, d1 = emulate_generator_key(yy)
        K1 &= MASK64
        K2 &= MASK64
        h = emulate_final_home(K1, K2, S)
        x8 = h * 8                              # bank-relative head-cell offset
        x15, x13 = header[0], header[1]         # ldp x15, x13, [x25]
        x16 = x15 + 40                          # add x16, x15, #40
        ok &= (x13 == B + cap_entries * 40)     # arena end reloaded per point
        if x16 > x13:
            return dict(heads=heads, arena=arena, cands=cands, ok=False, full=True, B=B)
        x9 = heads[h]                           # ldr x9, [x8, #32]
        arena.append(dict(off=x15, next=x9, K1=K1, K2=K2, w=w, j=j))
        heads[h] = x15                          # str x15, [x8, #32]
        header[0] = x16                         # str x16, [x25]
        while x9:                               # listing F
            en = arena[(x9 - B) // 40]
            ok &= (en["off"] == x9)
            if en["K1"] == K1 and en["K2"] == K2:
                cands.append((en["w"], en["j"], w, j))
            x9 = en["next"]                     # ldr x9, [x9]
        ok &= (pl == [(p >> (64 * i)) & MASK64 for i in range(3)])
    return dict(heads=heads, arena=arena, cands=cands, ok=ok, full=False, B=B)


def cost_block(claim):
    out, chk = {}, {}
    for k in ("GEN_key", "MSG", "NET", "RCV", "SRV_insert", "LOOP_CONTROL", "CAND_COUNTER"):
        chk["listing_count_" + k] = COUNT[k] == LEDGER[k]
    chk["listing_F_common_path_is_9"] = (COUNT["F_occupied_per_slot"] == 16
                                        and LEDGER["F_occupied_per_slot"] == 9)
    chk["F_block_and_carry_paths_priced"] = (LEDGER["F_block"] == 4 and LEDGER["F_carry"] == 3
                                            and COUNT["F_occupied_per_slot"]
                                            - LEDGER["F_occupied_per_slot"] == 7)
    # Distributed protocol listings (proof.md sections 2 and 5).  The generator
    # and the table server are separate processors; every address in the server
    # body is a 64-bit local pointer, there is no selector register, no barrier
    # and no >64-bit address construction anywhere.
    GEN = LISTINGS["GEN_key"]; F = LISTINGS["F_occupied_per_slot"]
    MSG = LISTINGS["MSG"]; NET = LISTINGS["NET"]; RCV = LISTINGS["RCV"]
    SRV = LISTINGS["SRV_insert"]
    wGEN, wF = written_regs(GEN), written_regs(F)
    chk["generator_preserves_p_limbs_x10_x11_x12"] = not ({"x10", "x11", "x12"} & wGEN)
    chk["generator_reloads_walk_id_x22_per_point"] = ("x22" in wGEN
                                                      and "ldp x2, x22, [x19, #16]" in GEN)
    chk["generator_forms_key_and_digits"] = ("ubfx x6, x0, #1, #36" in GEN
                                             and "ubfx x7, x1, #0, #35" in GEN)
    chk["message_stores_payload_and_descriptor"] = ("stp x1, x2, [x14]" in MSG
                                                    and "stp x22, x23, [x14, #16]" in MSG
                                                    and "stp x6, x7, [x14, #32]" in MSG)
    chk["message_advances_queue_pointer"] = "add x14, x14, #56" in MSG
    chk["network_is_two_levels_of_three"] = (COUNT["NET"] == 6 and NET.count("br x11") == 2
                                             and NET.count("ldr x11, [x12, x9, lsl #3]") == 2)
    chk["receiver_loads_payload_and_advances_queue"] = (
        COUNT["RCV"] == 6 and "ldp x1, x2, [x9]" in RCV
        and "ldp x22, x23, [x9, #16]" in RCV and "add x9, x9, #56" in RCV
        and "cmp x9, x13" in RCV and "mov x9, x10" in RCV)
    chk["server_insert_is_13_instructions"] = COUNT["SRV_insert"] == 13
    chk["server_links_before_scan"] = (SRV.index("str x9, [x15]")
                                       < SRV.index("str x15, [x8, #32]")
                                       < SRV.index("cbz x9, next"))
    chk["server_head_cell_offset_is_32"] = ("ldr x9, [x8, #32]" in SRV
                                            and "str x15, [x8, #32]" in SRV)
    chk["server_loads_arena_end_x13_per_point"] = "ldp x15, x13, [x25]" in SRV
    chk["F_traverses_x9_and_preserves_x8"] = ("ldp x3, x4, [x9, #8]" in F
                                               and "ldr x9, [x9]" in F
                                               and "x8" not in wF
                                               and not ({"x10", "x11", "x12"} & wF))
    alltext = "\n".join(LISTINGS[k] for k in LISTINGS)
    chk["no_mmio_selector_or_barrier_in_listings"] = ("dsb" not in alltext
                                                      and "selector" not in alltext.lower())
    chk["server_addresses_are_64_bit_locals"] = ("movk" not in alltext
                                                 and "add x8, x25, x8, lsl #3" in SRV
                                                 and "ldp x15, x13, [x25]" in SRV
                                                 and "ldr x9, [x8, #32]" in SRV)
    chk["two_generator_batches_cover_all_walks"] = N_GEN * LEDGER["batch_size"] == M_WALKS
    target = Fr(39, 100) + Fr(1, 10 ** 6)
    L, mu = solve_L(target)
    n = M_WALKS * (L + 1)
    # loop control is paid on the generator and on the table server, each
    # unrolled UNROLL times (the accepted batched-block convention).
    lc = Fr(2 * LEDGER["LOOP_CONTROL"], UNROLL)
    per_point = (LEDGER["GEN_operand_traffic"] + LEDGER["GEN_bookkeeping"]
                 + LEDGER["GEN_key"] + LEDGER["MSG"] + LEDGER["NET"]
                 + LEDGER["RCV"] + LEDGER["SRV_insert"] + lc
                 + Fr(LEDGER["batch_instr"], LEDGER["batch_size"])
                 + LEDGER["round_allowance"])
    price = 1 + per_point / C
    zero_units = Fr(S_TOT) * ZERO_PER_HEAD / C
    insp = Fr(n * (n - 1), 2 * S_TOT)                # expected chain comparisons E[X]
    probe_price = (Fr(LEDGER["F_occupied_per_slot"]) + Fr(LEDGER["F_block"], BLK)
                   + Fr(LEDGER["F_carry"], 2 ** 64))
    cap_frac = Fr(1005, 1000) * insp
    cap_hi = int(cap_frac // 2 ** 64)       # high-word cap: abort when the new high word reaches it
    cap_eff = Fr(cap_hi) * 2 ** 64          # comparisons allowed before the abort
    hidden = N_NODES * BLK                  # <= 2^87 unflushed per-server comparisons at the abort
    probe_allot = (cap_eff + hidden) * probe_price / C
    W_det = n * price + zero_units + probe_allot
    extras = dict(walk_setup=M_WALKS * 2 ** 14, network_setup=2 ** 90,
                  server_headers=N_NODES * 32,
                  candidates=2 ** 90 * 2 ** 14,
                  recovery_cap=2 ** 39 * 4 * 2 ** 17, solve_output=2 ** 20)
    W = W_det + sum(extras.values())
    out.update(L_log2=round(lg(L), 6), mu=float(mu), f4_mu=float(f4(mu)), poisson=1 - math.exp(-float(mu)),
               points_log2=round(lg(n), 6), per_point_instr=float(per_point), step_price=float(price),
               zero_units_log2=round(lg(zero_units), 3), probe_allot_log2=round(lg(probe_allot), 3),
               W_det_log2=round(lg(W_det), 6), W_log2=round(lg(W), 6),
               extras_log2={k: round(lg(v), 3) for k, v in extras.items()},
               inspections_log2=round(lg(insp), 6), taylor=round(float(insp) / n, 6))
    chk["mu_gives_f4_at_least_0.390001"] = f4(mu) >= target
    chk["L_below_2^64_j_fits_one_word"] = L < 2 ** 64
    chk["L_squared_much_below_r"] = 2 * (L + 1) ** 2 < r
    chk["claim_covers_W"] = lg(W) <= claim
    chk["claim_margin_below_0.001"] = claim - lg(W) < 0.001
    chk["claim_is_smallest_1e-5_value"] = math.ceil(lg(W) * 10 ** 5) / 10 ** 5 == claim
    chk["inspection_counter_is_shared_128_bit"] = (
        "ldadd x27, [x30], x31" in LISTINGS["F_occupied_per_slot"]
        and "ldadd x28, [x32], x33" in LISTINGS["F_occupied_per_slot"]
        and "tst x26, #0xffff" in LISTINGS["F_occupied_per_slot"])
    chk["candidate_counter_is_shared_128_bit"] = (
        COUNT["CAND_COUNTER"] == 6
        and "ldadd x29, [x30], x31" in LISTINGS["CAND_COUNTER"]
        and "ldadd x27, [x33], x31" in LISTINGS["CAND_COUNTER"])
    chk["unrolled_loop_control_two_sided"] = (UNROLL == 128 and LEDGER["LOOP_CONTROL"] == 2
                                              and lc == Fr(4, 128))
    chk["candidate_cap_costs_at_most_2^104_units"] = 2 ** 90 * 2 ** 14 == 2 ** 104
    chk["network_setup_at_most_2^90_units"] = extras["network_setup"] == 2 ** 90
    chk["server_headers_at_most_2^76_bytes"] = extras["server_headers"] == 2 ** 76
    chk["no_address_mapping_term"] = "address_mapping" not in extras
    # Multiplicity-aware inspection moments (proof.md section 4).  Under H' the
    # buckets are iid uniform per *class*; condition on the class partition
    # c_1..c_m (which depends only on the walk coins, not on the buckets).  Then
    # X = X_same + X_distinct, X_same = sum_a C(c_a,2) forced same-class
    # comparisons, E[X_same] = 6 C(n,2)/r = mu, and
    # X_distinct = sum_{a<b} c_a c_b 1[B_a = B_b] has
    # E = (n^2 - sum c_a^2)/(2S) <= insp and Var = p(1-p) sum (c_a c_b)^2
    # <= c_max^2 E, because P(B_a=B_b=B_c) = 1/S^2 = P(B_a=B_b)P(B_a=B_c) makes
    # the three-class covariance exactly 0 and disjoint pairs are independent.
    EX = insp
    varX = 64 * insp                        # c_max <= 8 except with prob < 2^-771
    cap_frac = Fr(1005, 1000) * insp        # 1.005 x E[X_distinct]
    cap_hi = int(cap_frac // 2 ** 64)       # high-word cap in x28
    cap_eff = Fr(cap_hi) * 2 ** 64          # effective cap enforced by the listing
    slack_eff = cap_eff - insp - 2 ** 30    # X_same < 2^30 with failure 2^-2^30*29
    probe_tail_eff = varX / slack_eff ** 2
    mu_same = Fr(6 * (n * (n - 1) // 2), r)
    cmax_fail = Fr(n ** 8, math.factorial(8)) * Fr(6 ** 7, r ** 7)
    chk["inspection_cap_fits_128_bit_counter"] = cap_hi < 2 ** 64
    chk["inspection_effective_slack_at_least_0.49pct"] = slack_eff >= Fr(49, 10000) * insp
    chk["inspection_overrun_below_2^-30"] = lg(probe_tail_eff + mu_same / 2 ** 30) < -30
    chk["class_multiplicity_at_most_8_at_2^-400"] = lg(cmax_fail) < -400
    # failure terms (section 4)
    q_reject = Fr(2 ** 255 - r, 2 ** 255)               # probability one 255-bit draw is rejected
    sampling_cap = 3 * M_WALKS * q_reject ** 128        # c_w, b_w, D_w draws
    # Candidate tail (proof.md section 4): multiplicity-weighted Chebyshev.
    # Z = sum_{a<b} c_a c_b 1[F_a=F_b] over 164-bit class fingerprints;
    # E[Z] <= C(n,2) 2^-164, Var[Z] <= c_max^2 E[Z] <= 256 E[Z]; cap 2^90.
    Ecan = Fr(n * (n - 1), 2 * 2 ** 164)
    var_can = 64 * Ecan
    cap_can = 2 ** 90 - 2 ** 30
    cand_tail = var_can / (cap_can - Ecan) ** 2
    same_tail = mu_same / 2 ** 30
    log2_same_tail = lg(same_tail)
    fail = dict(zero_denominator_batches=3 * float(Fr(n, r)),
                dependent_sets=2.0 ** -60,
                degenerate_b=float(Fr(6 * (M_WALKS * (M_WALKS - 1) // 2) * (L + 1) ** 2, r * r)),
                exceptional_points=3 * float(Fr(n, r)),
                arena_overflow=2.0 ** -400,
                probe_allotment=float(probe_tail_eff) + float(same_tail),
                candidate_allotment=float(cand_tail) + float(same_tail),
                sampling_cap=float(sampling_cap), recovery_cap=2.0 ** -120)
    succ = float(f4(mu)) - sum(fail.values())
    out["failure_terms"] = fail
    out["sampler_cap_log2"] = round(lg(sampling_cap), 3)
    out["probe_tail_log2"] = round(min(lg(probe_tail_eff), log2_same_tail), 3)
    out["inspection_overrun_log2"] = round(lg(probe_tail_eff + mu_same / 2 ** 30), 3)
    out["candidate_tail_log2"] = round(lg(cand_tail), 3)
    out["same_class_tail_log2"] = round(log2_same_tail, 3)
    out["inspection_cap_high_word_log2"] = round(math.log2(cap_hi), 4)
    chk["sampler_cap_below_2^-370"] = lg(sampling_cap) < -370
    chk["candidate_chebyshev_tail_below_2^-84"] = lg(cand_tail) < -84
    chk["same_class_markov_tail_below_2^-31"] = log2_same_tail < -31
    out["success_lower_bound"] = succ
    chk["success_at_least_0.39"] = succ >= 0.39
    # memory: 2^71 table servers, each with a local 64-bit-addressable bank
    n_node_mean = Fr(n, N_NODES)
    arena_cap = int(ARENA_MARGIN * n_node_mean) + 1
    bank = 32 + S_NODE * 8 + arena_cap * 40
    mem = (S_TOT * 8 + N_NODES * 32 + N_NODES * arena_cap * 40
           + N_GEN * 32 + M_WALKS * (504 + 96) + 2 ** 71 * 64)
    out["memory_log2_bytes"] = round(lg(mem), 3)
    out["arena_capacity_log2_entries"] = round(math.log2(arena_cap), 3)
    out["per_server_bank_log2_bytes"] = round(math.log2(bank), 3)
    chk["per_server_arena_below_2^63_bytes"] = arena_cap * 40 < 2 ** 63
    chk["per_server_heads_below_2^63_bytes"] = S_NODE * 8 < 2 ** 63
    chk["per_server_bank_below_2^61.4_bytes"] = bank < 2 ** 61.4
    chk["generator_record_store_below_2^63_bytes"] = LEDGER["batch_size"] * 600 < 2 ** 63
    chk["memory_below_2^132.3"] = lg(mem) < 132.3
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
    S_node = 2 ** 57
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
    res["imm_tst_ffff_encodable"] = logical_imm_ok(0xffff)
    res["ubfx_1_36_valid"] = 1 + 36 <= 64
    res["final_home_hashes_K2_directly"] = ("extr" not in LISTINGS["GEN_key"]
                                              and "umulh x8, x2, x20" in LISTINGS["SRV_insert"])
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
    points = []
    for j in range(Lw + 1):
        if j > 0:
            cur = batched_affine_add(cur, stepP, p)
        for w in range(Mw):
            points.append((cur[w][1], w, j))
    GF = emulate_listing_GF(points, S_small, len(points) + 4)
    heads, arena, cands = GF["heads"], GF["arena"], GF["cands"]
    entries = len(arena)
    # every entry is reachable exactly once from its bucket head; a prepended
    # entry's next is an earlier arena entry (index < i) or 0
    reach = 0
    for h in range(S_small):
        e, seen = heads[h], 0
        while e:
            reach += 1
            seen += 1
            if seen > entries: break
            e = arena[(e - GF["B"]) // 40]["next"]
    res["listing_GF_dataflow_simulation"] = GF["ok"] and not GF["full"]
    base = GF["B"]
    res["chained_insert_links_consistent"] = (reach == entries
        and all(arena[i]["next"] == 0
                or ((arena[i]["next"] - base) % 40 == 0
                    and arena[i]["next"] < arena[i]["off"])
                for i in range(entries)))
    res["listing_stores_locator_per_point"] = all(
        arena[i]["w"] == points[i][1] and arena[i]["j"] == points[i][2]
        for i in range(entries))

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


# ------------------------------------------------------------ Part A2 (H')
def part_a2(n_walks=64, n_steps=32):
    """Heuristic-H' consistency evidence at production bit positions.

    Real BLS12-381 G1 AP walks Q_{w,j} = [c_w + j D_w]P0 + [b_w]T (T = [x]P0),
    so the key bits examined here are the production bits of actual subgroup
    points.  Reports (a) bit balance over the 199 bits that H' treats as
    uniform (K1, K2, digit1, digit2), (b) pairwise bit-independence chi-square,
    (c) small-table Poisson occupancy and same-bucket pair counts, (d) counts of
    truncated 164-bit fingerprint collisions on high/middle/low bit slices, and
    (e) within-walk vs cross-walk collision rates.  This is consistency
    evidence on a scaled sample, NOT a proof of H'; H' is only used for
    cost-side expectations and concentration in proof.md section 4.
    """
    trng = random.Random(20261030)
    x = trng.randrange(1, r)
    T = ec_mul(x, P0, p)
    n_walks = max(1, min(128, int(n_walks)))
    n_steps = max(1, min(64, int(n_steps)))
    walks, cur, stepP = [], [], []
    for _ in range(n_walks):
        c, b, D = trng.randrange(r), trng.randrange(r), trng.randrange(1, r)
        walks.append((c, b, D))
        cur.append(ec_add(ec_mul(c, P0, p), ec_mul(b, T, p), p))
        stepP.append(ec_mul(D, P0, p))
    keys, degenerate = [], 0
    for j in range(n_steps + 1):
        if j:
            nxt = batched_affine_add(cur, stepP, p)
            if nxt is None:
                degenerate += 1
                nxt = [ec_add(a, s, p) for a, s in zip(cur, stepP)]
            cur = nxt
        for w in range(n_walks):
            if cur[w] is None:
                degenerate += 1
                continue
            K0, K1, K2, d1 = emulate_generator_key(cur[w][1])
            keys.append((K0, K1, K2, d1, w, j))
    m = len(keys)
    res = dict(walks=n_walks, steps=n_steps, points=m, degenerate_batches=degenerate)

    # (a) bit balance of the 199 bits H' treats as uniform: K1 (64), K2 (64),
    # digit1 (36), digit2 (35).  K0 bit 0 is the forced parity bit, not tested.
    rows = []
    for (K0, K1, K2, d1, w, j) in keys:
        d2 = K1 & ((1 << 35) - 1)
        rows.append([(K1 >> t) & 1 for t in range(64)]
                    + [(K2 >> t) & 1 for t in range(64)]
                    + [(d1 >> t) & 1 for t in range(36)]
                    + [(d2 >> t) & 1 for t in range(35)])
    nbits = len(rows[0])
    zs = [(2 * sum(r[b] for r in rows) - m) / math.sqrt(m) for b in range(nbits)]
    res["bit_balance_tests"] = nbits
    res["bit_balance_max_abs_z"] = round(max(abs(z) for z in zs), 3)

    # (b) pairwise bit-independence chi-square (2x2 tables, 1 dof) on K1||K2
    pairs = [((11 * i) % 128, (97 * i + 13) % 128) for i in range(64)]
    chi = []
    for a, b in pairs:
        if a == b:
            b = (b + 1) % 128
        n11 = sum(1 for r in rows if r[a] and r[b])
        n10 = sum(1 for r in rows if r[a] and not r[b])
        n01 = sum(1 for r in rows if not r[a] and r[b])
        n00 = m - n11 - n10 - n01
        den = (n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00)
        if den:
            chi.append(m * (n11 * n00 - n10 * n01) ** 2 / den)
    res["pairwise_tests"] = len(chi)
    res["pairwise_max_chi2"] = round(max(chi), 3)

    # (c) occupancy of a scaled S-bucket table; home = umulh(K2, S) exactly as
    # in the server listing, Poisson(lambda = m/S) prediction
    for S in (64, 1024):
        cnt = [0] * S
        for (K0, K1, K2, d1, w, j) in keys:
            cnt[(K2 * S) >> 64] += 1
        lam = m / S
        hist = {}
        for c in cnt:
            hist[c] = hist.get(c, 0) + 1
        cells, acc_o, acc_e, chi2 = [], 0, 0.0, 0.0
        for c in range(0, max(cnt) + 2):
            acc_o += hist.get(c, 0)
            acc_e += S * (lam ** c) * math.exp(-lam) / math.factorial(c)
            if acc_e >= 5:
                cells.append((acc_o, acc_e)); acc_o, acc_e = 0, 0.0
        if acc_e > 0:
            cells.append((acc_o, acc_e))
        for o, e in cells:
            if e > 0:
                chi2 += (o - e) ** 2 / e
        dof = max(1, len(cells) - 1)
        pairs_obs = sum(c * (c - 1) // 2 for c in cnt)
        res["occupancy_S%d" % S] = dict(
            lambda_=round(lam, 3), chi2=round(chi2, 2), dof=dof,
            pairs_observed=pairs_obs, pairs_expected=round(m * (m - 1) / (2 * S), 1))
        res["occupancy_S%d_ok" % S] = chi2 < dof + 6 * math.sqrt(2 * dof) + 10

    # (d) truncated 164-bit fingerprint F = K1 || K2 || digit1 collision counts
    F = [((K1 << 100) | (K2 << 36) | d1) for (K0, K1, K2, d1, w, j) in keys]

    def npairs(tags):
        d = {}
        for t in tags:
            d[t] = d.get(t, 0) + 1
        return sum(c * (c - 1) // 2 for c in d.values())

    fpc, maxz = {}, 0.0
    for b in (12, 14, 16, 18, 20):
        for name, tags in (("high", [f >> (164 - b) for f in F]),
                           ("mid", [(f >> 60) & ((1 << b) - 1) for f in F]),
                           ("low", [f & ((1 << b) - 1) for f in F])):
            obs = npairs(tags)
            exp = m * (m - 1) / 2 / 2 ** b
            z = (obs - exp) / math.sqrt(exp) if exp > 0 else 0.0
            fpc["b%d_%s" % (b, name)] = dict(observed=obs, expected=round(exp, 2),
                                             z=round(z, 2))
            maxz = max(maxz, abs(z))
    res["fingerprint_collisions"] = fpc
    res["fingerprint_max_abs_z"] = round(maxz, 3)

    # (e) within-walk vs cross-walk collision rates at b = 16 (serial
    # correlation check); both are compared with 2^-16.
    b = 16
    tags = [f >> (164 - b) for f in F]
    wpairs = wobs = 0
    for w in range(n_walks):
        tw = [tags[i] for i, k in enumerate(keys) if k[4] == w]
        wpairs += len(tw) * (len(tw) - 1) // 2
        wobs += npairs(tw)
    tpairs = m * (m - 1) // 2
    cobs, cpairs = npairs(tags) - wobs, tpairs - wpairs
    wrate = (wobs / wpairs) if wpairs else 0.0
    crate = (cobs / cpairs) if cpairs else 0.0
    wexp, cexp = wpairs * 2.0 ** -b, cpairs * 2.0 ** -b
    wz = (wobs - wexp) / math.sqrt(wexp) if wexp > 0 else 0.0
    cz = (cobs - cexp) / math.sqrt(cexp) if cexp > 0 else 0.0
    res["within_cross_b16"] = dict(within_pairs=wpairs, within_observed=wobs,
                                   within_expected=round(wexp, 3), within_z=round(wz, 3),
                                   cross_pairs=cpairs, cross_observed=cobs,
                                   cross_expected=round(cexp, 3), cross_z=round(cz, 3),
                                   within_rate=round(wrate, 8),
                                   cross_rate=round(crate, 8),
                                   expected_rate=round(2.0 ** -b, 8),
                                   ratio=round(wrate / crate, 3) if crate else None)

    res["points_at_least_2000"] = m >= 2000
    res["bit_balance_ok"] = res["bit_balance_max_abs_z"] < 5.0
    res["pairwise_ok"] = res["pairwise_max_chi2"] < 30.0
    res["fingerprint_ok"] = maxz < 5.0
    res["within_cross_ok"] = abs(wz) < 5.0 and abs(cz) < 5.0
    res["all_true"] = all(bool(res[k]) for k in
                          ("points_at_least_2000", "bit_balance_ok", "pairwise_ok",
                           "fingerprint_ok", "within_cross_ok"))
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
    S = int((Mw * (Lw + 1)) / 0.30) + 1
    heads, arena = {}, []           # direct buckets with chained heads
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
            h = K * S // q            # toy home: multiply-high of the key by S
            e = heads.get(h)
            while e is not None:
                if arena[e][0] == K and found is None:
                    cands += 1
                    w1, j1 = arena[e][1]
                    Q1 = ec_add(ec_mul((walks[w1][0] + j1 * walks[w1][2]) % rr, G, q), ec_mul(walks[w1][1], T, q), q)
                    for t in range(6):
                        if unit_apply(t, P, beta, q) == Q1:
                            uu = U[t]
                            (ca, ba, Da), (cb, bb, Db) = walks[w1], walks[w]
                            den = (ba - uu * bb) % rr
                            if den:
                                xs = ((uu * (cb + j * Db) - ca - j1 * Da) * pow(den, -1, rr)) % rr
                                if ec_mul(xs, G, q) == T: found = xs
                e = arena[e][2]
            arena.append((K, (w, j), heads.get(h)))
            heads[h] = len(arena) - 1
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
    claim = cfg.get("claim", 126.18004)
    out["cost"] = cost_block(claim)
    out["listing_counts"] = COUNT
    if cfg.get("production_checks", 1):
        out["part_a"] = part_a(cfg.get("production_random", 12))
        out["part_a2"] = part_a2(cfg.get("a2_walks", 64), cfg.get("a2_steps", 32))
    trng = random.Random(SEED + 7)
    if cfg.get("toy_seconds", 0) > 0:
        out["part_b_toy_curve"] = part_b(cfg.get("curve", "q30"), cfg.get("toy_walks", 32), cfg["toy_seconds"], trng)
    if cfg.get("scalar_seconds", 0) > 0:
        out["part_b_scalar_model"] = part_b_scalar(cfg.get("scalar_rbits", 22), cfg.get("scalar_walks", 64), cfg["scalar_seconds"], trng)
    out["all_cost_checks_true"] = all(out["cost"]["checks"].values())
    out["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(out, indent=1, default=str))


main()
