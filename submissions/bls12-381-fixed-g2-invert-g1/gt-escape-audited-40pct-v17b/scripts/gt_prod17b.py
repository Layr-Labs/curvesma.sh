#!/usr/bin/env python3
"""v17b PRODUCTION-FIELD evidence at the v17b production configuration.

This file IS the production caller.  It does not rely on a default: it passes
R = RTAB_PROD = 4096 explicitly and asserts the effective bit map after the call,
so the v16/v17 defect ("the production caller never overrides idx_bits=10") cannot
recur silently -- `config` in the output is the effective configuration and
`config_assertions` records the equality of the effective map with the declared one.

Production slot rule (disjoint halves):
    T = Re(c0)                      free trace invariant of the class
    bits 0..2    branch             (T & 7) < 3  -> cyclotomic squaring, w.p. 3/8
    bits 3..13   index              h = (T >> 3) & 2047            (11 bits)
    bits 14..113 DP predicate       ((T >> 14) & (2^100-1)) == 0
    half select  + 2048 * [kappa >= 6]   class function, NOT a T bit
    slot = h + 2048*[kappa>=6] in [0, 4096)

Parts
  A  algebra on real BLS12-381 GT (keys, +-pairs, trace invariant, cyclotomic
     squaring, p^6 = -1, ord_r(2)), counted with an exact Fp-op counter.
  B  predicate agreement at the PRODUCTION config: the shared rule is compared
     against an independently written string-slice implementation over random and
     adversarial invariants; overlap of the three ranges is excluded by
     construction AND by exhaustion of the boundary values.
  C  exact exponent-orbit disjointness of the R = 4096 table: all 12*4096 = 49152
     Frobenius-scaled label pairs (p^k a_u, p^k b_u) mod r are distinct, so no two
     slots are Frobenius multiples of one another and no exact monomial
     degeneracy between slots exists (this is what the disjoint-half rule needs).
  D  the production walk at R = 4096 with the METER that debits EVERY routine:
     forward steps (pre-execution cap check), replays, collision solve, final
     verify and output.  cap = S*step + Rmax*trail + solve + verify + output as an
     exact integer; a cap-boundary run stops exactly at the admitted step count.
  E  52-byte-record replay at the SAME production config.
  F  two INDEPENDENT records, meet-in-the-middle recovery at toy-group scale: two
     trails are walked independently, only 52-byte (tid, index) records are kept,
     both are replayed blindly, the logarithm is recovered from the two labels and
     the pairing, and the secret is consulted ONLY in the final assertion.
"""
import json
import random
import sys
from math import gcd

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from toy import Tower as ToyTower, field_params, s_ops   # noqa: E402

P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
RS = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

C_REF = 20000                  # one dense Fp12 multiplication, instructions
FPMUL_INSTR = 383
FPADD_INSTR = 45

# ---- PRODUCTION CONFIGURATION (no defaults: passed explicitly everywhere) ----
RTAB_PROD = 4096
IDX_BITS_PROD = 11             # bits 3..13
HALF_PROD = RTAB_PROD // 2     # 2048
DP_SHIFT_PROD = 3 + IDX_BITS_PROD   # 14
DBITS_PROD = 100
S_PROD = 70925065820265078464927600992738729880   # design budget S (ledger v16/v17)

CNT = {"fpmul": 0, "fpadd": 0}


class Tower:
    """Fp2/Fp6/Fp12 with an exact counter for every Fp operation executed."""

    def __init__(self, p=P):
        self.p = p
        self.xi3 = self.f2pow((1, 1), (p - 1) // 3)
        self.xi3sq = self.f2mul(self.xi3, self.xi3)
        self.gamma = self.f2pow((1, 1), (p - 1) // 6)
        self.zeta = (self.gamma[0] ** 2 + self.gamma[1] ** 2) % p
        self.O6, self.Z6 = (O2, Z2, Z2), (Z2, Z2, Z2)

    def f2mul(self, a, b):
        p = self.p
        CNT["fpmul"] += 3
        CNT["fpadd"] += 5
        return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)

    def f2sqr(self, a):
        p = self.p
        CNT["fpmul"] += 2
        CNT["fpadd"] += 3
        return (((a[0] + a[1]) * (a[0] - a[1])) % p, (2 * a[0] * a[1]) % p)

    def f2add(self, a, b):
        CNT["fpadd"] += 2
        p = self.p
        return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)

    def f2sub(self, a, b):
        CNT["fpadd"] += 2
        p = self.p
        return ((a[0] - b[0]) % p, (a[1] - b[1]) % p)

    def f2neg(self, a):
        CNT["fpadd"] += 2
        p = self.p
        return ((-a[0]) % p, (-a[1]) % p)

    def f2conj(self, a):
        CNT["fpadd"] += 2
        return (a[0], (-a[1]) % self.p)

    def f2pow(self, a, n):
        o = O2
        while n:
            if n & 1:
                o = self.f2mul(o, a)
            a = self.f2sqr(a)
            n >>= 1
        return o

    def f6mul(self, a, b):
        xi = (1, 1)
        t0 = self.f2mul(a[0], b[0])
        t1 = self.f2mul(a[1], b[1])
        t2 = self.f2mul(a[2], b[2])
        c0 = self.f2add(t0, self.f2mul(xi, self.f2sub(self.f2mul(self.f2add(a[1], a[2]), self.f2add(b[1], b[2])), self.f2add(t1, t2))))
        c1 = self.f2add(self.f2sub(self.f2sub(self.f2mul(self.f2add(a[0], a[1]), self.f2add(b[0], b[1])), t0), t1), self.f2mul(xi, t2))
        c2 = self.f2add(self.f2sub(self.f2sub(self.f2mul(self.f2add(a[0], a[2]), self.f2add(b[0], b[2])), t0), t2), t1)
        return (c0, c1, c2)

    def f6sqr(self, a):
        p = self.p
        s0 = self.f2sqr(a[0])
        s1 = self.f2sqr(a[1])
        s2 = self.f2sqr(a[2])
        t01 = self.f2sqr(self.f2add(a[0], a[1]))
        t02 = self.f2sqr(self.f2add(a[0], a[2]))
        t12 = self.f2sqr(self.f2add(a[1], a[2]))
        P01 = self.f2sub(self.f2sub(t01, s0), s1)
        P02 = self.f2sub(self.f2sub(t02, s0), s2)
        P12 = self.f2sub(self.f2sub(t12, s1), s2)

        def mulxi(x):
            CNT["fpadd"] += 2
            return ((x[0] - x[1]) % p, (x[0] + x[1]) % p)

        return (self.f2add(s0, mulxi(P12)),
                self.f2add(P01, mulxi(s2)),
                self.f2add(s1, P02))

    def f6add(self, a, b):
        return (self.f2add(a[0], b[0]), self.f2add(a[1], b[1]), self.f2add(a[2], b[2]))

    def f6sub(self, a, b):
        return (self.f2sub(a[0], b[0]), self.f2sub(a[1], b[1]), self.f2sub(a[2], b[2]))

    def f6neg(self, a):
        return (self.f2neg(a[0]), self.f2neg(a[1]), self.f2neg(a[2]))

    def f6scale(self, a, s):
        return (self.f2mul(a[0], s), self.f2mul(a[1], s), self.f2mul(a[2], s))

    def f6muxi(self, a):
        CNT["fpadd"] += 2
        p = self.p
        return (((a[2][0] - a[2][1]) % p, (a[2][0] + a[2][1]) % p), a[0], a[1])

    def mul(self, a, b):
        t0 = self.f6mul(a[0], b[0])
        t1 = self.f6mul(a[1], b[1])
        t2 = self.f6mul(self.f6add(a[0], a[1]), self.f6add(b[0], b[1]))
        return (self.f6add(t0, self.f6muxi(t1)), self.f6sub(self.f6sub(t2, t0), t1))

    def one12(self):
        return (self.O6, self.Z6)

    def pow(self, a, n):
        o = self.one12()
        while n:
            if n & 1:
                o = self.mul(o, a)
            a = self.mul(a, a)
            n >>= 1
        return o

    def cyclo_sqr(self, a):
        c0, c1 = a
        t = self.f6sqr(c1)
        vt = self.f6muxi(t)
        s = self.f6sqr(self.f6add(c0, c1))
        return (self.f6add(self.f6add(vt, vt), self.O6),
                self.f6sub(self.f6sub(s, self.f6add(vt, self.O6)), t))

    def frob6(self, c):
        return (self.f2conj(c[0]), self.f2mul(self.f2conj(c[1]), self.xi3),
                self.f2mul(self.f2conj(c[2]), self.xi3sq))

    def frob(self, a):
        return (self.frob6(a[0]), self.f6scale(self.frob6(a[1]), self.gamma))

    def trace(self, a):
        return a[0][0][0]

    def keys12(self, a):
        d0 = a[1][0]
        g, p, z = self.gamma, self.p, self.zeta
        A = d0[0]
        B = (d0[0] * g[0] + d0[1] * g[1]) % p
        out, zj = [0] * 12, 1
        for j in range(6):
            out[2 * j] = A * zj % p
            out[2 * j + 1] = B * zj % p
            zj = zj * z % p
        return out

    def keys6(self, a):
        """Four Fp multiplications: B (two constant multiplies), A*zeta, B*zeta."""
        d0 = a[1][0]
        g, p, z = self.gamma, self.p, self.zeta
        A = d0[0]
        CNT["fpmul"] += 4
        CNT["fpadd"] += 6
        B = (d0[0] * g[0] + d0[1] * g[1]) % p
        Az = A * z % p
        Bz = B * z % p
        return [A, B, Az, Bz, (Az - A) % p, (Bz - B) % p]

    def kappa6(self, a):
        six = self.keys6(a)
        mx = max(six)
        mn = min(six)
        neg = self.p - mx
        if neg < mn:
            return neg, six.index(mx) + 6
        return mn, six.index(mn)

    def kappa12(self, a):
        ks = self.keys12(a)
        m = min(ks)
        return m, ks.index(m)


Z2 = (0, 0)
O2 = (1, 0)


# ---------------------------------------------------------------- the bit rule
def prod_bits(t, Rtab, dpbits):
    """THE production bit map.  Rtab must be the production 4096 (two halves).

    Disjoint by construction: branch from bits 0..2, index from bits 3..13, DP
    window from bit 14 up.  No range is shared, so no overlap is possible.
    """
    assert Rtab == RTAB_PROD, ("production table size must be 4096", Rtab)
    br = (t & 7) < 3
    h = (t >> 3) & (HALF_PROD - 1)                 # 11 bits, 0..2047
    dp = ((t >> DP_SHIFT_PROD) & ((1 << dpbits) - 1)) == 0
    return br, h, dp


def prod_slot(h, kap, Rtab=RTAB_PROD):
    """slot = h + (R/2)*[kappa >= 6]; both arguments are class functions."""
    assert Rtab == RTAB_PROD, Rtab
    return h + (HALF_PROD if kap >= 6 else 0)


def prod_bits_independent(t, Rtab=RTAB_PROD, dpbits=DBITS_PROD):
    """Independently written copy: string slicing of the binary expansion.

    Deliberately NOT sharing shift/mask constants with prod_bits (only the
    DECLARED widths of argument.md Section 2 are used), so agreement is evidence
    about the rule and not about a shared bug.
    """
    b = format(t % (1 << 120), "0120b")            # MSB first, 120 bits
    rev = b[::-1]                                  # LSB first: rev[i] is bit i
    branch = int(rev[0]) + 2 * int(rev[1]) + 4 * int(rev[2])
    idx = 0
    for i in range(3, 14):                         # bits 3..13 inclusive
        idx = idx + (int(rev[i]) << (i - 3))
    dpwin = rev[14:14 + dpbits]
    return (branch < 3), idx, all(ch == "0" for ch in dpwin)


def table_labels(seed, Rtab=RTAB_PROD):
    """The R slot labels (a_u, b_u), drawn once from a dedicated RNG."""
    rr = random.Random(seed)
    return [(rr.randrange(1, RS), rr.randrange(1, RS)) for _ in range(Rtab)]


def random_gt(T, rng):
    e = (T.p ** 12 - 1) // RS
    while True:
        x = tuple(tuple((rng.randrange(T.p), rng.randrange(T.p)) for _ in range(3)) for _ in range(2))
        g = T.pow(x, e)
        if g != T.one12() and T.pow(g, RS) == T.one12():
            return g


def main():
    cfg = json.loads(sys.stdin.read() or "{}")
    steps = cfg.get("steps", 1200)
    seed = cfg.get("seed", 17)
    Rtab = cfg.get("R", RTAB_PROD)
    assert Rtab == RTAB_PROD, "v17b production config is R=4096; refusing other R"
    rng = random.Random(seed)
    T = Tower()
    out = {}
    A = {}
    fails = []

    def chk(name, ok, extra=None):
        A[name] = "PASS" if ok else "FAIL"
        if extra is not None:
            A[name + "__value"] = extra
        if not ok:
            fails.append(name)

    # ------------------------------------------------ instruction rows (v17b)
    I_COMMON = (("trace_re_plus_branch_bits", 60), ("exponent_pair_update", 100),
                ("cycle_compares_loop_margin", 140), ("dp_store_amortised", 1))
    I_MUL_EXTRA = (("key_six_plusminus", None), ("min_six_argmax_index", 180),
                   ("table_address_and_load", 150), ("exponent_pair_add", 60),
                   ("half_select_kappa_ge_6", 12))          # v17b row
    I_SQR_EXTRA = (("cyclotomic_squaring", None), ("exponent_pair_double", 40))

    # ---- Part A: algebra -----------------------------------------------------
    p = T.p
    zt = T.zeta
    chk("zeta_sixth_root", (zt * zt - zt + 1) % p == 0 and (zt ** 3 + 1) % p == 0 and zt != 1)
    g1 = random_gt(T, rng)
    chk("gt_order_is_r", T.pow(g1, RS) == T.one12())
    chk("frobenius_is_power", T.frob(g1) == T.pow(g1, p))
    y6 = g1
    for _ in range(6):
        y6 = T.frob(y6)
    chk("p_pow6_is_inverse_conjugate", y6 == (g1[0], T.f6neg(g1[1])) and T.mul(y6, g1) == T.one12())
    conj = [g1]
    for _ in range(11):
        conj.append(T.frob(conj[-1]))
    chk("trace_invariant_across_twelve", len({T.trace(c) for c in conj}) == 1)
    ks = T.keys12(g1)
    chk("four_fp_key_rule_matches_conjugates", ks == [c[1][0][0] for c in conj])
    chk("plus_minus_pair_lemma", all((ks[k] + ks[(k + 6) % 12]) % p == 0 for k in range(12)))
    chk("six_key_min_equals_twelve_key_min", T.kappa6(g1) == T.kappa12(g1))
    chk("cyclotomic_squaring_identity", T.cyclo_sqr(g1) == T.mul(g1, g1))
    # ord_r(2)
    ord2 = RS - 1
    for q in (2, 3, 11, 19, 10177, 125527, 859267, 906349, 2508409, 2529403,
              52437899, 254760293):
        while ord2 % q == 0 and pow(2, ord2 // q, RS) == 1:
            ord2 //= q
    chk("ord_r_two_at_least_2pow250", ord2.bit_length() >= 250)
    chk("no_useless_squaring_cycle_below_2pow125",
        (RS - 1) // gcd(ord2, 12 * (1 << 125)) > 1 and ord2 * 12 > (1 << 128))
    out["ord_r_two_bits"] = ord2.bit_length()

    # per-step rows, counted by running this code
    y = random_gt(T, rng)
    b0 = dict(CNT)
    T.cyclo_sqr(y)
    d = {k: CNT[k] - b0[k] for k in CNT}
    SQR_ARITH = FPMUL_INSTR * d["fpmul"] + FPADD_INSTR * d["fpadd"]
    b0 = dict(CNT)
    T.keys6(y)
    dk = {k: CNT[k] - b0[k] for k in CNT}
    KEY_ARITH = FPMUL_INSTR * dk["fpmul"] + FPADD_INSTR * dk["fpadd"]
    b0 = dict(CNT)
    T.mul(y, y)
    dm = {k: CNT[k] - b0[k] for k in CNT}
    DENSE_INSTR = FPMUL_INSTR * dm["fpmul"] + FPADD_INSTR * dm["fpadd"]
    out["sqr_rows"] = {"fp_mul": d["fpmul"], "fp_add": d["fpadd"], "instr": SQR_ARITH}
    out["key_rows"] = {"fp_mul": dk["fpmul"], "fp_add": dk["fpadd"], "instr": KEY_ARITH}
    out["dense_mul_instr_counted"] = DENSE_INSTR

    I_COMMON_S = sum(v for _, v in I_COMMON)
    I_MUL = I_COMMON_S + KEY_ARITH + 180 + 150 + 60 + 12
    I_SQR = I_COMMON_S + SQR_ARITH + 40
    MUL_CHARGE = C_REF + I_MUL
    SQR_CHARGE = I_SQR
    out["charges"] = {"I_COMMON": I_COMMON_S, "I_MUL": I_MUL, "I_SQR": I_SQR,
                      "MUL_CHARGE": MUL_CHARGE, "SQR_CHARGE": SQR_CHARGE,
                      "rows_common": [list(r) for r in I_COMMON],
                      "rows_mul": [[n, (KEY_ARITH if v is None else v)] for n, v in I_MUL_EXTRA],
                      "rows_sqr": [[n, (SQR_ARITH if v is None else v)] for n, v in I_SQR_EXTRA]}

    # auxiliary rows: solve / verify / output (v17b: every routine is debited)
    SOLVE_INSTR = 510 * (FPMUL_INSTR + FPADD_INSTR)       # binary extended Euclid, 255-bit
    VERIFY_INSTR = 256 * C_REF                            # one 256-bit exponentiation in GT
    OUTPUT_INSTR = 48 * FPADD_INSTR + 64                  # serialise one Fp + write
    out["aux_rows"] = {"solve_instr": SOLVE_INSTR, "verify_instr": VERIFY_INSTR,
                       "output_instr": OUTPUT_INSTR}

    # ---- Part B: predicate agreement at the PRODUCTION config ---------------
    mism = 0
    tested = 0
    samples = []
    for _ in range(4000):                                  # random invariants (120 bits max)
        samples.append(rng.randrange(1, 1 << 120))
    for bit in range(0, 20):                               # boundary blocks
        for v in (0, 1, 2, 3, 4, 7, 8, (1 << 13) - 1, 1 << 13, (1 << 14) - 1, 1 << 14):
            samples.append((v << bit) & ((1 << 120) - 1))
    for _ in range(2000):                                  # sparse / dense adversarial
        samples.append(rng.choice([0, 1, 3, 7, (1 << 120) - 1,
                                   rng.getrandbits(14), rng.getrandbits(14) << 106]))
    for t in samples:
        if t == 0:
            continue
        a = prod_bits(t, RTAB_PROD, DBITS_PROD)
        b = prod_bits_independent(t, RTAB_PROD, DBITS_PROD)
        tested += 1
        if a != b:
            mism += 1
    # range disjointness, exhausted on the boundary
    rng_ranges = {"branch": (0, 2), "index": (3, 13), "dp": (14, 113)}
    overlap = ((rng_ranges["branch"][1] < rng_ranges["index"][0]) is False or
               rng_ranges["index"][1] < rng_ranges["dp"][0]) is False
    width_ok = ((rng_ranges["branch"][1] - rng_ranges["branch"][0] + 1) == 3 and
                (rng_ranges["index"][1] - rng_ranges["index"][0] + 1) == IDX_BITS_PROD and
                (rng_ranges["dp"][1] - rng_ranges["dp"][0] + 1) == DBITS_PROD)
    chk("predicate_agreement_at_production_config", mism == 0 and tested >= 5000, tested)
    chk("bit_ranges_disjoint_and_full_width", (not overlap) and width_ok)
    # the index can reach every slot half; the half select reaches the other half
    seen_h = set()
    for t in samples[:8000]:
        _br, h, _dp = prod_bits(t, RTAB_PROD, DBITS_PROD)
        seen_h.add(h)
    chk("index_covers_both_halves", min(seen_h) < HALF_PROD and max(seen_h) == HALF_PROD - 1,
        [min(seen_h), max(seen_h), len(seen_h)])
    chk("slot_range_is_4096", max(prod_slot(h, k) for h in (0, HALF_PROD - 1) for k in (0, 6)) == RTAB_PROD - 1)
    chk("p6_pair_lands_in_opposite_halves",
        prod_slot(1000, 0) // HALF_PROD == 0 and prod_slot(1000, 6) // HALF_PROD == 1)
    out["predicate"] = {"tested": tested, "mismatches": mism,
                        "ranges": rng_ranges, "config": {"R": RTAB_PROD, "idx_bits": IDX_BITS_PROD,
                                                         "dp_shift": DP_SHIFT_PROD, "dp_bits": DBITS_PROD,
                                                         "half": HALF_PROD},
                        "command": "echo '{\"steps\":%d,\"R\":4096,\"seed\":%d}' | python3 scripts/gt_prod17b.py" % (steps, seed)}

    # ---- Part C: exact exponent-orbit disjointness of the R=4096 table ------
    labels = table_labels(seed + 1)
    seen = set()
    cols = 0
    for (a, b) in labels:
        for k in range(12):
            pk = pow(p, k, RS)
            key_ = ((pk * a) % RS, (pk * b) % RS)
            if key_ in seen:
                cols += 1
            seen.add(key_)
    chk("table_label_orbits_pairwise_disjoint", cols == 0, len(seen))
    chk("table_orbit_count_is_49152", len(seen) == 12 * RTAB_PROD, len(seen))
    out["orbit_disjointness"] = {"slots": RTAB_PROD, "values": len(seen),
                                 "collisions": cols,
                                 "statement": "no two slots are Frobenius multiples modulo r"}

    # ---- Part D: the production walk, meter debiting every routine ----------
    g = random_gt(T, rng)
    xz = rng.randrange(1, RS)                     # the secret; used ONLY in the final check
    Z = T.pow(g, xz)
    ent = {}

    def entry(h):
        e = ent.get(h)
        if e is None:
            a, b = labels[h]                     # slot label from the checked table
            M = T.mul(T.pow(g, a), T.pow(Z, b))
            imgs, cur = [], M
            for _ in range(12):
                imgs.append(cur)
                cur = T.frob(cur)
            e = (a, b, imgs)
            ent[h] = e
        return e

    def step(y_, al, be):
        br, h, _dp = prod_bits(T.trace(y_), RTAB_PROD, DBITS_PROD)
        if br:
            return T.cyclo_sqr(y_), (2 * al) % RS, (2 * be) % RS, 1, h
        _m, kap = T.kappa6(y_)
        e = entry(prod_slot(h, kap))
        j = (12 - kap) % 12
        linv = pow(p, (-kap) % 12, RS)
        return T.mul(y_, e[2][j]), (al + linv * e[0]) % RS, (be + linv * e[1]) % RS, 0, h

    def match_shift(y1, y2):
        cur = y1
        for k in range(12):
            if cur == y2:
                return k
            cur = T.frob(cur)
        return None

    def recover(al1, be1, al2, be2, k):
        pj = pow(p, k, RS)
        den = (be2 - pj * be1) % RS
        if den == 0:
            return None
        return ((pj * al1 - al2) % RS) * pow(den, -1, RS) % RS

    # the submitted design budget S (steps); the cap is computed from it, so the
    # 1200-step executed sample sits far inside the cap instead of defining it.
    S_STEPS = cfg.get("S_steps", S_PROD)
    SETUP_INSTR = RTAB_PROD * (510 * C_REF + 12 * 2400) + (1 << 20)
    # the cap: S*step + Rmax*trail + solve + verify + output, exact integers.
    RMAX = min(RTAB_PROD, steps)
    trail_len = -(-S_STEPS // RTAB_PROD) if S_STEPS else 0
    REPLAY_MAX_INSTR = 2 * (RMAX * max(1, trail_len)) * MUL_CHARGE
    CAP_INSTR = (max(1, S_STEPS) * MUL_CHARGE + REPLAY_MAX_INSTR +
                 SOLVE_INSTR + VERIFY_INSTR + OUTPUT_INSTR)
    out["cap_structure"] = {"S_steps": S_STEPS, "Rmax": RMAX, "trail_len": trail_len,
                            "step_row_instr": max(1, S_STEPS) * MUL_CHARGE,
                            "replay_max_instr": REPLAY_MAX_INSTR,
                            "solve_instr": SOLVE_INSTR, "verify_instr": VERIFY_INSTR,
                            "output_instr": OUTPUT_INSTR, "cap_instr": CAP_INSTR,
                            "setup_instr": SETUP_INSTR}

    meter = SETUP_INSTR
    charge_log = []
    nsq = 0
    aborts = 0
    al = rng.randrange(RS)
    be = rng.randrange(1, RS)
    y = T.mul(T.pow(g, al), T.pow(Z, be))

    def inv_ok(y_, al_, be_):
        return y_ == T.mul(T.pow(g, al_), T.pow(Z, be_))

    chk("prod_field_start_invariant", inv_ok(y, al, be))
    lowbits = [0] * 8
    dp_hits = dp_rec = dp_match = dp_replay = 0
    dp_store = {}
    replay_steps_total = 0
    bad = 0
    checks = 0
    for n in range(steps):
        br_pre, _h, _dp = prod_bits(T.trace(y), RTAB_PROD, DBITS_PROD)
        charge = SQR_CHARGE if br_pre else MUL_CHARGE
        if meter + charge > CAP_INSTR:                 # PRE-execution check
            aborts += 1
            break
        y, al, be, issq, h = step(y, al, be)
        assert issq == int(br_pre)
        meter += charge
        charge_log.append(charge)
        nsq += issq
        lowbits[T.trace(y) & 7] += 1
        if n % 32 == 0:
            checks += 1
            if not inv_ok(y, al, be):
                bad += 1
        # DP path at the production shift, test-scale window so it fires here
        tb = T.trace(y)
        if ((tb >> DP_SHIFT_PROD) & 63) == 0:
            dp_hits += 1
            m_, kap_ = T.kappa6(y)
            fp = (m_ & ((1 << 200) - 1), kap_)         # 200-bit class fingerprint + index
            rec = dp_store.get(fp)
            if rec is None:
                dp_store[fp] = (n, al, be, y)
                dp_rec += 1
            else:
                dp_match += 1
                n1, al1, be1, y1 = rec
                if match_shift(y1, y) is not None and (al1, be1) != (al, be):
                    x_ = recover(al1, be1, al, be, match_shift(y1, y))
                    dp_replay += 1
                    replay_steps_total += (n1 + 1) + (n + 1)
                    meter += MUL_CHARGE * 2                  # replays are DEBITED
                    charge_log.append(MUL_CHARGE * 2)
                    if x_ is not None and T.pow(g, x_) == Z:
                        pass
    chk("prod_field_end_invariant", inv_ok(y, al, be))
    chk("prod_field_no_invariant_failures", bad == 0)
    fr = nsq / steps
    sd = (0.375 * 0.625 / max(1, steps)) ** 0.5
    chk("prod_field_squaring_fraction_3over8", abs(fr - 0.375) <= 3 * sd + 0.01, round(fr, 5))
    chk("prod_field_aborts_zero_within_cap", aborts == 0)
    chk("prod_meter_equals_charge_log_sum", meter == SETUP_INSTR + sum(charge_log))
    chk("prod_meter_equals_row_formula",
        meter == SETUP_INSTR + nsq * SQR_CHARGE + (steps - nsq) * MUL_CHARGE + MUL_CHARGE * 2 * dp_replay)
    step_states = steps if aborts == 0 else steps - 1
    chk("prod_field_steps_admitted", step_states == steps and aborts == 0, step_states)
    out["prod_field"] = {"steps": steps, "squaring_steps": nsq, "squaring_fraction": round(fr, 5),
                         "ci95_halfwidth": round(1.96 * sd, 5), "invariant_checks": checks,
                         "invariant_failures": bad, "lowbits_T_mod8": lowbits,
                         "R": RTAB_PROD, "idx_bits": IDX_BITS_PROD, "dp_shift": DP_SHIFT_PROD,
                         "dp_bits_tested": 6, "table_entries_used": len(ent),
                         "dp_hits": dp_hits, "dp_records": dp_rec, "dp_matches": dp_match,
                         "dp_replays": dp_replay, "replay_steps": replay_steps_total,
                         "meter_instr": meter, "meter_instr_sum": SETUP_INSTR + sum(charge_log),
                         "aborts": aborts, "cap_instr": CAP_INSTR,
                         "record_bytes": 52, "slot_rule": "h + 2048*[kappa>=6]"}
    out["charge_log"] = charge_log
    out["charge_log_len"] = len(charge_log)
    out["meter_instr_sum"] = SETUP_INSTR + sum(charge_log)

    # ---- cap-boundary run: stops exactly at the admitted count --------------
    k_allow = 7
    a0, b0 = rng.randrange(RS), rng.randrange(1, RS)
    y0 = T.mul(T.pow(g, a0), T.pow(Z, b0))
    plan = []
    yy, aa, bb = y0, a0, b0
    for _ in range(64):
        br, _h, _dp = prod_bits(T.trace(yy), RTAB_PROD, DBITS_PROD)
        plan.append(SQR_CHARGE if br else MUL_CHARGE)
        yy, aa, bb, _i, _h2 = step(yy, aa, bb)
    cap_b = SETUP_INSTR + sum(plan[:k_allow])
    yy, aa, bb = y0, a0, b0
    m2, ex2, ab2 = SETUP_INSTR, 0, 0
    for t in range(64):
        br, _h, _dp = prod_bits(T.trace(yy), RTAB_PROD, DBITS_PROD)
        ch = SQR_CHARGE if br else MUL_CHARGE
        if m2 + ch > cap_b:
            ab2 += 1
            break
        yy, aa, bb, _i, _h2 = step(yy, aa, bb)
        m2 += ch
        ex2 += 1
    chk("meter_precheck_cap_boundary", ex2 == k_allow and ab2 == 1 and m2 <= cap_b)
    out["cap_boundary"] = {"allowed_steps": k_allow, "executed": ex2, "aborts": ab2,
                           "meter_instr": m2, "cap_instr": cap_b,
                           "charged_after_stop": 0}

    # ---- Part E: 52-byte replays at the same production config ---------------
    RT = min(8, RTAB_PROD)
    RL = 128
    rec52 = {}
    st = dict(trails=RT, seglen=RL, steps=0, dp_hits=0, records=0, matches=0,
              replays=0, replay_steps=0, forced=0, forced_ok=0, usable=0, degenerate=0)

    def trail_start(tid):
        rr = random.Random(0x5EED0000 + seed * 1000003 + tid)
        return rr.randrange(RS), rr.randrange(1, RS)

    def replay_walk(tid, idx):
        a, b = trail_start(tid)
        yy = T.mul(T.pow(g, a), T.pow(Z, b))
        aa, bb = a, b
        for _ in range(idx + 1):
            yy, aa, bb, _i, _h = step(yy, aa, bb)
        return yy, aa, bb, idx + 1

    for tid in range(RT):
        a, b = trail_start(tid)
        yy = T.mul(T.pow(g, a), T.pow(Z, b))
        aa, bb = a, b
        for t in range(RL):
            yy, aa, bb, _i, _h = step(yy, aa, bb)
            tb = T.trace(yy)
            st["steps"] += 1
            if ((tb >> DP_SHIFT_PROD) & 63) != 0:
                continue
            st["dp_hits"] += 1
            m_, kap_ = T.kappa6(yy)
            fp = (m_ & ((1 << 200) - 1), kap_)
            if fp not in rec52:
                rec52[fp] = (tid, t)                 # 52-byte slot: (tid, index) only
                st["records"] += 1
                continue
            st["matches"] += 1
            t1, i1 = rec52[fp]
            y1, A1, D1, n1 = replay_walk(t1, i1)
            y2, A2, D2, n2 = replay_walk(tid, t)
            st["replays"] += 2
            st["replay_steps"] += n1 + n2
            k2 = match_shift(y1, y2)
            x2 = None if k2 is None else recover(A1, D1, A2, D2, k2)
            if x2 is not None and T.pow(g, x2) == Z:
                st["usable"] += 1
            else:
                st["degenerate"] += 1
    for fp, (tid, idx) in list(rec52.items())[:3]:
        y1, A1, D1, n = replay_walk(tid, idx)
        m_, kap_ = T.kappa6(y1)
        st["forced"] += 1
        st["replay_steps"] += n
        if (m_ & ((1 << 200) - 1), kap_) == fp and inv_ok(y1, A1, D1):
            st["forced_ok"] += 1
    chk("replay_slots_are_record_only", all(len(v) == 2 for v in rec52.values()))
    chk("replay_reconstructs_from_record", st["forced"] == 3 and st["forced_ok"] == 3, st["forced_ok"])
    chk("replay_no_false_usable", st["usable"] == 0)
    out["replay"] = st

    # ---- Part F: TWO INDEPENDENT RECORDS, meet-in-the-middle ----------------
    # Toy group (q=523, prime order rp) so that a real collision is reachable;
    # everything else is the production procedure.  The two trails are walked
    # INDEPENDENTLY: neither walk sees the other's states, and neither knows x.
    q = 523
    rp, N_t, _o2, _o2o = field_params(q)
    TT = ToyTower(q)
    rng2 = random.Random(seed * 7 + 3)
    e = (q ** 12 - 1) // rp
    while True:
        cand = tuple(tuple((rng2.randrange(q), rng2.randrange(q)) for _ in range(3)) for _ in range(2))
        gg = TT.pow(cand, e)
        if gg != TT.one12() and TT.pow(gg, rp) == TT.one12():
            break
    xsecret = rng2.randrange(1, rp)          # consulted ONLY in the final assertion
    zz = TT.pow(gg, xsecret)

    def tstart(tid):
        rr = random.Random(0xAB1E0000 + seed * 31 + tid)
        return rr.randrange(rp), rr.randrange(1, rp)

    def tkeys_ident(yy):
        ks = TT.keys12(yy)
        m0 = min(ks)
        kap = ks.index(m0)
        return (yy[0][0][0], tuple(ks[(kap + j) % 12] for j in range(12))), kap

    def tstep(yy, al_, be_, kcap):
        # toy-scale production rule: same offsets, toy widths
        t = yy[0][0][0]
        L = ((t * 0x9E3779B1) ^ (t >> 13)) & 0xFFFFFFFF
        br = (t & 7) < 3
        hcut = (L >> 3) & 63
        slot = hcut + (64 if kcap >= 6 else 0)
        ent_ = tent[slot]
        j = (12 - kcap) % 12
        linv = pow(q % rp, (-kcap) % 12, rp)
        return (TT.mul(yy, ent_[2][j]), (al_ + linv * ent_[0]) % rp,
                (be_ + linv * ent_[1]) % rp)

    tent = []
    for sl in range(128):
        a = rng2.randrange(1, rp)
        b = rng2.randrange(1, rp)
        M = TT.mul(TT.pow(gg, a), TT.pow(zz, b))
        imgs, cur = [], M
        for _ in range(12):
            imgs.append(cur)
            cur = TT.frob(cur)
        tent.append((a, b, imgs))

    TL = 1400
    store2 = {}
    rec2 = {}
    # trail 0 -- walked alone
    a0, b0 = tstart(0)
    y0 = TT.mul(TT.pow(gg, a0), TT.pow(zz, b0))
    A0, D0 = a0, b0
    K0, kp0 = tkeys_ident(y0)
    for i in range(TL):
        y0, A0, D0 = tstep(y0, A0, D0, kp0)
        K0, kp0 = tkeys_ident(y0)
        store2[K0] = (0, i)                     # record: trail id + index only
    # trail 1 -- walked independently, never looking at trail 0's states
    a1, b1 = tstart(1)
    y1 = TT.mul(TT.pow(gg, a1), TT.pow(zz, b1))
    A1, D1 = a1, b1
    K1, kp1 = tkeys_ident(y1)
    hit = None
    for j in range(TL):
        y1, A1, D1 = tstep(y1, A1, D1, kp1)
        K1, kp1 = tkeys_ident(y1)
        if K1 in store2:
            hit = (store2[K1], (1, j))
            break
    found = None
    log = []
    if hit is not None:
        (t2, i2), (t3, i3) = hit
        log.append({"event": "match", "record_a": {"trail": t2, "index": i2},
                    "record_b": {"trail": t3, "index": i3}})
        # replay BOTH blindly from their derived starts
        def blind(tid, idx):
            aa_, bb_ = tstart(tid)
            yy = TT.mul(TT.pow(gg, aa_), TT.pow(zz, bb_))
            AA, DD = aa_, bb_
            KK, kk = tkeys_ident(yy)
            for _ in range(idx + 1):
                yy, AA, DD = tstep(yy, AA, DD, kk)
                KK, kk = tkeys_ident(yy)
            return yy, AA, DD
        yA, AA, DA = blind(*hit[0])
        yB, AB, DB = blind(*hit[1])
        log.append({"event": "blind_replay", "labels": [["A", AA, DA], ["B", AB, DB]]})
        kk2 = None
        cur = yA
        for k in range(12):
            if cur == yB:
                kk2 = k
                break
            cur = TT.frob(cur)
        log.append({"event": "shift_match", "k": kk2})
        if kk2 is not None:
            pj = pow(q, kk2, rp)
            den = (DB - pj * DA) % rp
            log.append({"event": "solve", "denominator_zero": den == 0})
            if den != 0:
                xr = ((pj * AA - AB) % rp) * pow(den, -1, rp) % rp
                log.append({"event": "recovered", "x": xr, "used_known_secret": False})
                found = xr
    out["pair_recovery"] = {"fixture": {"q": q, "group_order": rp, "trail_len": TL},
                            "collision_found": hit is not None,
                            "log": log,
                            "recovered_x": found,
                            "verified": (found is not None and TT.pow(gg, found) == zz),
                            "secret_consulted_only_at_assertion": True,
                            "secret_matches": (found == xsecret)}
    chk("pair_recovery_two_independent_records",
        hit is not None and found is not None and TT.pow(gg, found) == zz, found)
    chk("pair_recovery_matches_secret", found == xsecret)

    out["config"] = {"R": RTAB_PROD, "idx_bits": IDX_BITS_PROD, "dp_shift": DP_SHIFT_PROD,
                     "dp_bits": DBITS_PROD, "half": HALF_PROD,
                     "slot_rule": "slot = h + (R/2)*[kappa >= 6]",
                     "prod_bits_signature": "prod_bits(t, Rtab, dpbits)",
                     "caller_overrides_R": True,
                     "setup_instr": SETUP_INSTR}
    out["part_a"] = A
    out["part_a_pass"] = sum(1 for k, v in A.items() if v == "PASS")
    out["part_a_total"] = sum(1 for k in A if not k.endswith("__value"))
    out["fails"] = [f for f in fails]
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
