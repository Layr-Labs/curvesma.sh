"""Exact baby-step giant-step (negation keys) for BLS12-381 G1: production ledger and exhaustive toy runs.

Part 1 recomputes every production number quoted in proof.md (parameters, exact success count, itemised
instruction charges, W, the rounded-up score and peak memory) with exact integer and Fraction arithmetic.
Part 2 runs the same algorithm (same chain layout, same incomplete batched affine additions, same identity
and exceptional-case checks, a real LSD radix sort and merge) on small prime-order curves y^2 = x^3 + b and
recovers every scalar of each toy group exhaustively, checking the success count against the formula.
Toy runs are evidence at toy scale only. Standard library only.
"""
import json
import sys
import time
from fractions import Fraction

# ----------------------------------------------------------------------------------------------------------
# Part 1: production ledger (single source of every number in proof.md, claim.json and submission.json)
# ----------------------------------------------------------------------------------------------------------
P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
C = 2000            # fixed reference price (instructions per reference operation) for this track
MUL = 383           # pinned blst_fp_mul instruction count
K = 2 ** 20         # parallel chains per list (batch size of every batched addition round)
SUCCESS_TARGET = Fraction(2, 5)
DIGIT_BITS = 32     # LSD radix digit width
DIGITS = 12         # 12 x 32 = 384 >= 381 key bits
REC_BYTES = 64      # 48-byte key (6 limbs) + 16-byte index

# Itemised instruction charges (scalar AArch64 instructions). Big arrays (the two lists and the scatter
# buffer) are charged 2 instructions per 64-bit word accessed: the access plus one for the high part of the
# 128-bit address. Small arrays (chain state, counters) are charged 1 per word.
SUB = 6 + 1 + 6 + 6 + 6                      # modular 384-bit subtraction: subs/sbcs, borrow mask, and, adds/adcs, loads of p
ADD_EXTRA = [                                # per batched incomplete affine addition, beyond the reference unit
    ("6 modular subtractions outside the 6M (dx, dy, lam^2-x1, -x2, x1-x3, -y1)", 6 * SUB),
    ("load chain point x1,y1 (12 words)", 12),
    ("load constant addend x2,y2 (12 words)", 12),
    ("store dx and prefix product, reload both in the backward pass (24 words)", 24),
    ("store result x3,y3 (12 words)", 12),
    ("loop control", 3),
]
GIANT_CHECK = [("exceptional check x(A)==x(S): 6 eor, 5 orr, 1 branch", 12)]
INVERSION = 240000   # Fermat a^(p-2): 380 squarings + 228 multiplications = 608 M, plus 5 per exponent bit
REC = [                                      # per emitted record (each baby and each giant point)
    ("canonical reduction of x: 6 limb loads of p, 6 subs/sbcs, 6 csel", 18),
    ("histogram update for 12 digits: extract, address, ldp, adds/adc, stp (8 each)", DIGITS * 8),
    ("store 8-word record (2 per word), 128-bit pointer advance, 128-bit index update", 16 + 2 + 2),
]
SCATTER = [                                  # per record per LSD pass
    ("load 8-word record (2 per word)", 16),
    ("extract 32-bit digit", 1),
    ("counter address", 1),
    ("load 128-bit counter (2 words)", 2),
    ("destination address: 128-bit shift by 6 (3) and 128-bit add (2)", 5),
    ("increment counter (adds/adc) and store it (2 words)", 4),
    ("store 8-word record (2 per word)", 16),
    ("128-bit source pointer advance (2) and 128-bit end test with branch (3)", 5),
]
COUNTER_PASS = 3 + 9                         # per counter per digit: zeroing (stp + loop 2) and prefix sum (ldp, adds/adc, stp, loop 3)
MERGE = [                                    # per merge step (each step advances one pointer)
    ("load both 6-word keys (2 per word)", 24),
    ("compare 6 limbs from the top: cmp + branch each, final branch", 13),
    ("128-bit pointer advance", 2),
    ("two 128-bit bounds tests with branch", 6),
]
GLAW = 250000                                # one textbook affine group-law operation with all cases (O, P=-Q,
                                             # P=Q doubling, generic): 1 inversion, <= 4 M, <= 10 sub/compare, 50 control
SMUL = 128000000                             # double-and-add, 255 iterations of 2 group-law ops + 24 control
FINISH = 2 * SMUL + 1000                     # two candidate scalars verified by scalar multiplication; input
                                             # decoding, -G, -S, candidate arithmetic and output (1000)


def total(items):
    return sum(n for _, n in items)


def checks():
    assert SUB == 25
    assert INVERSION + 4 * MUL + 10 * SUB + 50 <= GLAW and 255 * (2 * GLAW + 24) <= SMUL
    assert (380 + 228) * MUL + 381 * 5 <= INVERSION
    assert bin(P - 2).count("1") == 229 and (P - 2).bit_length() == 381
    assert DIGITS * DIGIT_BITS >= P.bit_length()


def instr_per_element(giant):
    add = total(ADD_EXTRA) + (total(GIANT_CHECK) if giant else 0)
    return add + total(REC) + DIGITS * total(SCATTER) + total(MERGE)


def work(n_order, kk, lb, lg):
    """Exact W = H + Rins/C for chain count kk, chain lengths lb (babies) and lg (giants)."""
    m, g = kk * lb, kk * lg
    H = kk * (lb - 1) + kk * (lg - 1)                       # batched incomplete affine additions
    parts = [
        ("baby additions: extra instructions", kk * (lb - 1) * total(ADD_EXTRA)),
        ("giant additions: extra instructions incl. exceptional check",
         kk * (lg - 1) * (total(ADD_EXTRA) + total(GIANT_CHECK))),
        ("one inversion per batched round", ((lb - 1) + (lg - 1)) * INVERSION),
        ("record emission (m + g records)", (m + g) * total(REC)),
        ("LSD scatter passes (12 per record)", (m + g) * DIGITS * total(SCATTER)),
        ("counter zeroing and prefix sums", 2 * DIGITS * (2 ** DIGIT_BITS) * COUNTER_PASS),
        ("merge (m + g steps charged)", (m + g) * total(MERGE)),
        ("chain starts and S: 2K+1 scalar mults, K group-law ops", (2 * kk + 1) * SMUL + kk * GLAW),
        ("final verification and I/O", FINISH),
    ]
    Rins = sum(v for _, v in parts)
    work.parts = parts
    return m, g, H, Rins, H + Fraction(Rins, C)


def choose(n_order, kk):
    """Smallest W over chain lengths with g*(2m+1) >= 2/5 * n and g*(2m+1) + m < n - s."""
    a = 1 + Fraction(instr_per_element(False), C)
    b = 1 + Fraction(instr_per_element(True), C)
    need = SUCCESS_TARGET * n_order
    best = None
    # continuous optimum m* = sqrt(b*need/(2a)); scan chain lengths around it
    import math
    m0 = math.isqrt(int(b * need / (2 * a)) + 1)
    lb0 = max(2, m0 // kk)
    for lb in range(max(2, lb0 - 3), lb0 + 4):
        m = kk * lb
        lg = max(2, -(-int(-(-need // (2 * m + 1))) // kk))
        while kk * lg * (2 * m + 1) < need:
            lg += 1
        g = kk * lg
        s = 2 * m + 1
        assert g * s + m < n_order - s and 2 * m < n_order
        w = work(n_order, kk, lb, lg)[4]
        if best is None or w < best[0]:
            best = (w, lb, lg)
    return best[1], best[2]


def score_ceiling(w):
    """Smallest k with W <= 2^(k/100), exactly."""
    k = 100 * (w.numerator.bit_length() - w.denominator.bit_length() - 1)
    while w.numerator ** 100 > 2 ** k * w.denominator ** 100:
        k += 1
    while k > 0 and w.numerator ** 100 <= 2 ** (k - 1) * w.denominator ** 100:
        k -= 1
    return k


def log2f(fr):
    import math
    return math.log2(fr.numerator) - math.log2(fr.denominator)


def production():
    checks()
    lb, lg = choose(R, K)
    m, g, H, Rins, W = work(R, K, lb, lg)
    s = 2 * m + 1
    covered = g * s
    k = score_ceiling(W)
    rounds = (lb - 1) + (lg - 1)
    counters = DIGITS * (2 ** DIGIT_BITS) * 16
    state = K * 256
    code = 2 ** 20
    peak = REC_BYTES * m + max(REC_BYTES * m, 2 * REC_BYTES * g) + counters + state + code
    out = {
        "K": K, "Lb": lb, "Lg": lg, "m": m, "g": g, "s": s,
        "covered": covered, "success_count": covered + 1,
        "success_lower": float(Fraction(covered, R)), "success_exact": float(Fraction(covered + 1, R)),
        "success_ok": Fraction(covered, R) >= SUCCESS_TARGET,
        "disjoint_ok": g * s + m < R - s,
        "H_additions": H, "log2_H": round(log2f(Fraction(H)), 6),
        "R_instructions": Rins, "R_over_C": round(Rins / C, 1),
        "instr_per_baby_element": instr_per_element(False),
        "instr_per_giant_element": instr_per_element(True),
        "add_extra": total(ADD_EXTRA), "rec": total(REC), "scatter_per_pass": total(SCATTER),
        "merge": total(MERGE), "rounds": rounds,
        "W_float": float(W), "log2_W": round(log2f(W), 6),
        "score": "%d/100" % k, "score_decimal": k / 100,
        "peak_bytes": peak, "log2_peak_bytes": round(log2f(Fraction(peak)), 4),
        "breakdown_instructions": [[name, v, round(log2f(Fraction(v, C)), 4)] for name, v in work.parts],
        "baby_phase_units_log2": round(log2f(K * (lb - 1) + Fraction(
            K * (lb - 1) * total(ADD_EXTRA) + (lb - 1) * INVERSION + m * (total(REC) + DIGITS * total(SCATTER))
            + DIGITS * (2 ** DIGIT_BITS) * COUNTER_PASS + (K + 1) * SMUL, C)), 4),
        "setup_units_log2": round(log2f(Fraction((2 * K + 1) * SMUL + K * GLAW, C)), 4),
    }
    assert out["success_ok"] and out["disjoint_ok"]
    return out


# ----------------------------------------------------------------------------------------------------------
# Part 2: toy runs of the same algorithm
# ----------------------------------------------------------------------------------------------------------
TOYS = [(1009, 11, 967), (4003, 2, 4111), (12043, 5, 11827)]   # (p, b, prime group order n) of y^2=x^3+b
TOY_K = 4
TOY_DIGIT_BITS = 4


class Exceptional(Exception):
    pass


class Toy:
    def __init__(self, p, b, n):
        self.p, self.b, self.n = p, b, n
        cnt = 1
        for x in range(p):
            v = (x * x * x + b) % p
            cnt += 1 if v == 0 else (2 if pow(v, (p - 1) // 2, p) == 1 else 0)
        assert cnt == n and all(n % d for d in range(2, int(n ** 0.5) + 1)), "toy order is not prime"
        self.G = self.first_point()
        self.counts = {}

    def bump(self, key, by=1):
        self.counts[key] = self.counts.get(key, 0) + by

    def first_point(self):
        p = self.p
        for x in range(1, p):
            v = (x * x * x + self.b) % p
            for y in range(1, p):
                if y * y % p == v:
                    return (x, y)
        raise AssertionError

    def gadd(self, P1, P2):
        """Textbook affine group law with every case: None is the identity."""
        p = self.p
        if P1 is None:
            return P2
        if P2 is None:
            return P1
        x1, y1 = P1
        x2, y2 = P2
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None                               # P2 = -P1 (covers 2-torsion, absent here)
            lam = 3 * x1 * x1 * pow(2 * y1, p - 2, p) % p  # doubling, a = 0
        else:
            lam = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
        x3 = (lam * lam - x1 - x2) % p
        self.bump("group_law")
        return (x3, (lam * (x1 - x3) - y1) % p)

    def smul(self, k, A):
        """Double-and-add with the textbook group law; k reduced mod n."""
        k %= self.n
        Rp = None
        for bit in bin(k)[2:] if k else "":
            Rp = self.gadd(Rp, Rp)
            if bit == "1":
                Rp = self.gadd(Rp, A)
        self.bump("scalar_mult")
        return Rp

    def batch_add(self, pts, addend):
        """Batched incomplete affine addition pts[c] + addend with one inversion; raises on x1 == x2."""
        p = self.p
        x2, y2 = addend
        pref, acc = [], 1
        for x1, y1 in pts:
            dx = (x2 - x1) % p
            if dx == 0:
                raise Exceptional
            pref.append(acc)
            acc = acc * dx % p
        inv = pow(acc, p - 2, p)
        self.bump("inversion")
        out = [None] * len(pts)
        for c in range(len(pts) - 1, -1, -1):
            x1, y1 = pts[c]
            dx = (x2 - x1) % p
            dinv = inv * pref[c] % p
            inv = inv * dx % p
            lam = (y2 - y1) * dinv % p
            x3 = (lam * lam - x1 - x2) % p
            out[c] = (x3, (lam * (x1 - x3) - y1) % p)
        self.bump("batched_add", len(pts))
        return out

    def neg(self, A):
        return (A[0], (-A[1]) % self.p)

    def radix_sort(self, recs):
        """LSD radix sort of (key, index) records on the full key, digit by digit; stable counting sort."""
        digits = -(-self.p.bit_length() // TOY_DIGIT_BITS)
        mask = (1 << TOY_DIGIT_BITS) - 1
        for d in range(digits):
            shift = d * TOY_DIGIT_BITS
            cnt = [0] * (mask + 1)
            for key, _ in recs:
                cnt[(key >> shift) & mask] += 1
            pos, run = [0] * (mask + 1), 0
            for v in range(mask + 1):
                pos[v], run = run, run + cnt[v]
            out = [None] * len(recs)
            for rec in recs:
                v = (rec[0] >> shift) & mask
                out[pos[v]] = rec
                pos[v] += 1
            recs = out
            self.bump("scatter_moves", len(recs))
        return recs

    def setup(self):
        """Q-independent part: baby list (sorted), S = sG, giant start offsets P_c = [c Lg s]G."""
        n = self.n
        self.Lb, self.Lg = choose(n, TOY_K)
        self.m, self.g = TOY_K * self.Lb, TOY_K * self.Lg
        self.s = 2 * self.m + 1
        negG = self.neg(self.G)
        B = [self.smul((c + 1) * self.Lb, self.G) for c in range(TOY_K)]   # chain tops, j = (c+1)Lb
        recs = [(B[c][0], (c + 1) * self.Lb) for c in range(TOY_K)]
        for t in range(1, self.Lb):
            B = self.batch_add(B, negG)          # descending chain j -> j-1, j >= 2: never exceptional
            recs += [(B[c][0], (c + 1) * self.Lb - t) for c in range(TOY_K)]
        assert sorted(j for _, j in recs) == list(range(1, self.m + 1))
        self.baby = self.radix_sort(recs)
        self.S = self.smul(self.s, self.G)
        self.negS = self.neg(self.S)
        self.offsets = [self.smul(n - c * self.Lg * self.s, self.G) for c in range(TOY_K)]   # -P_c
        self.setup_counts = dict(self.counts)
        self.counts = {}

    def verify(self, cands, Q):
        for xc in cands:
            if self.smul(xc % self.n, self.G) == Q:
                return xc % self.n
        return None

    def solve(self, Q):
        n, s, Lg = self.n, self.s, self.Lg
        A = []
        for c in range(TOY_K):
            start = self.gadd(Q, self.offsets[c])        # Q - c*Lg*s*G, all cases handled
            if start is None:                             # the identity: x = c*Lg*s
                return self.verify([c * Lg * s], Q), "start-identity"
            A.append(start)
        recs = [(A[c][0], c * Lg) for c in range(TOY_K)]
        for t in range(1, Lg):
            for c in range(TOY_K):
                if A[c][0] == self.S[0]:                  # A_i = +-S: x = (i+1)s or (i-1)s
                    i = c * Lg + t - 1
                    return self.verify([(i + 1) * s, (i - 1) * s], Q), "exceptional"
            A = self.batch_add(A, self.negS)
            recs += [(A[c][0], c * Lg + t) for c in range(TOY_K)]
        recs = self.radix_sort(recs)
        a = b = 0
        steps = 0
        while a < len(self.baby) and b < len(recs):
            steps += 1
            ka, kb = self.baby[a][0], recs[b][0]
            if ka == kb:
                j, i = self.baby[a][1], recs[b][1]
                self.bump("merge_steps", steps)
                return self.verify([i * s + j, i * s - j], Q), "match"
            if ka < kb:
                a += 1
            else:
                b += 1
        self.bump("merge_steps", steps)
        self.bump("max_merge", 0)
        self.counts["max_merge"] = max(self.counts["max_merge"], steps)
        return None, "fail"


def run_toy(spec):
    p, b, n = spec
    toy = Toy(p, b, n)
    t0 = time.time()
    toy.setup()
    m, g, s = toy.m, toy.g, toy.s
    G = toy.G
    predicted = set()
    for i in range(g):
        for e in range(-m, m + 1):
            predicted.add((i * s + e) % n)
    extra = (n - s) % n
    assert extra not in predicted
    predicted.add(extra)
    # enumerate every scalar x: Q = xG by a running sum (identity at x = 0)
    wrong = 0
    found = set()
    kinds = {}
    Qs = [None]
    for x in range(1, n):
        Qs.append(toy.gadd(Qs[-1], G))
    toy.counts = {}
    for x in range(n):
        out, kind = toy.solve(Qs[x])
        kinds[kind] = kinds.get(kind, 0) + 1
        if out is not None:
            if out != x:
                wrong += 1
            found.add(x)
    res = {
        "p": p, "b": b, "n": n, "K": TOY_K, "Lb": toy.Lb, "Lg": toy.Lg, "m": m, "g": g, "s": s,
        "formula_success_count": g * s + 1, "observed_success_count": len(found),
        "success_set_equals_formula_set": found == predicted, "wrong_outputs": wrong,
        "success_fraction": "%d/%d" % (len(found), n), "formula_fraction": "%d/%d" % (g * s + 1, n),
        "outcomes": kinds, "max_merge_steps": toy.counts.get("max_merge", 0), "merge_bound": m + g - 1,
        "baby_batched_adds": toy.setup_counts.get("batched_add", 0), "baby_formula": TOY_K * (toy.Lb - 1),
        "baby_list_sorted_and_complete": [r[0] for r in toy.baby] == sorted(r[0] for r in toy.baby),
        "seconds": round(time.time() - t0, 2),
    }
    return res, toy, predicted


def check_giant_counts(toy, predicted):
    """Re-run the first uncovered target (a full-length run ending in a failed merge) and compare its counted
    work with the ledger formula: every giant addition, every scatter move and at most m+g-1 merge steps."""
    x = min(v for v in range(1, toy.n) if v not in predicted)
    Q = toy.smul(x, toy.G)
    toy.counts = {}
    out, kind = toy.solve(Q)
    c = toy.counts
    want_scatter = (-(-toy.p.bit_length() // TOY_DIGIT_BITS)) * toy.g
    return {
        "x": x, "outcome": kind, "full_length_run": kind == "fail" and out is None,
        "giant_batched_adds": c.get("batched_add", 0), "formula": TOY_K * (toy.Lg - 1),
        "scatter_moves": c.get("scatter_moves", 0), "formula_scatter": want_scatter,
        "merge_steps": c.get("merge_steps", 0), "merge_bound": toy.m + toy.g - 1,
    }


def main():
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    args = json.loads(raw) if raw.strip() else {}
    if len(sys.argv) > 1:
        args = json.loads(sys.argv[1])
    groups = int(args.get("groups", 3))
    if not 1 <= groups <= 3:
        raise SystemExit("groups must be 1..3")
    prod = production()
    print("PRODUCTION LEDGER")
    print(json.dumps(prod, indent=1, default=str))
    ok = True
    for spec in TOYS[:groups]:
        res, toy, predicted = run_toy(spec)
        gc = check_giant_counts(toy, predicted)
        good = (res["success_set_equals_formula_set"] and res["wrong_outputs"] == 0
                and res["observed_success_count"] == res["formula_success_count"]
                and res["max_merge_steps"] <= res["merge_bound"]
                and res["baby_batched_adds"] == res["baby_formula"] and res["baby_list_sorted_and_complete"]
                and gc["full_length_run"] and gc["giant_batched_adds"] == gc["formula"]
                and gc["scatter_moves"] == gc["formula_scatter"] and gc["merge_steps"] <= gc["merge_bound"])
        res["full_run_count_check"] = gc
        res["all_checks_pass"] = good
        ok = ok and good
        print("TOY GROUP n=%d" % spec[2])
        print(json.dumps(res, indent=1))
    print("RESULT: %s" % ("all toy checks pass" if ok else "TOY CHECK FAILED"))
    print(json.dumps({"production_score": prod["score"], "log2_W": prod["log2_W"], "toy_checks_pass": ok}))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
