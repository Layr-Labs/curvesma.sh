"""Six-class (negation + GLV) parallel Pollard rho on secp256k1.

Port of the construction of submission 1c600df4 (mpjunior92, bls12-381-g1-ecdlp-plain,
glv-six-class-rho-v11) to secp256k1: same class walk, escapes, checkpoint and
re-walk solving; secp256k1 constants, a 4-limb truncated canonical-index rule,
toy curves y^2 = x^3 + 7, and the resource rules of this package's proof.md.

Part A checks, on secp256k1 itself, every algebraic identity the attack uses
(including the truncated-product canonical-index rule of proof.md section 1.1,
against the exact rule on random and adversarial near-tie inputs):
the GLV pair (beta, lambda), the size-6 orbits, the |y| class key, the
canonical-member rule, the twisted-table step identity, coefficient tracking
along real production walk steps, the doubling escape as a class map, and the
collision-solving formula (including its sign) on two production trails that
are made to merge.

Part B runs the complete algorithm of proof.md (twisted R-adding class walk,
2-cycle escape, checkpoint escape, trail cap, distinguished classes, trail
restart after every distinguished class, re-walk solving, verification) on the
prime-order subgroup of a toy curve y^2 = x^3 + 7 over F_q, q = 1 mod 3, and
measures the number of class visits (trail starts + walk operations) until the first collision between
different trails, against the random-mapping model on N = (r'-1)/6 classes
(or (r'-1)/classes for the classes = 1, 2 control walks).
Toy evidence only: it does not validate production cost bounds.
"""
import json, math, random, sys, time

T0 = time.time()
cfg = json.load(sys.stdin)
curve = cfg.get("curve", "q25")
trials = cfg.get("trials", 100)
R = cfg.get("R", 256)
dp_bits = cfg.get("dp_bits", 4)
seed = cfg.get("seed", 1)
prod_checks = cfg.get("production_checks", 2)
budget = cfg.get("seconds", 14)
dump = cfg.get("dump", False)
pipeline = cfg.get("pipeline", False)
long_trails = cfg.get("long_trails", 0)
long_len = cfg.get("long_len", 4096)
long_escape_add = cfg.get("long_escape_add", False)
canon_random = cfg.get("canon_random", None)
canon_stress = cfg.get("canon_stress", None)
run_toy = cfg.get("toy", True)
seq_n = cfg.get("sequence_check", 200000)
CLS = cfg.get("classes", 6)
rng = random.Random(seed)

# ---------------------------------------------------------------- curves
PROD = dict(
    q=int("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16),
    r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16),
    h=1,
    gx=int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16),
    gy=int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16),
    beta=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16),
    lam=int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72", 16),
)
CB = 7  # curve coefficient b (y^2 = x^3 + 7) for production and toys
TOYS = {
    "q25": (33556867, 33563419, 1),
    "q29": (536873251, 536839203, 3),
}


class Curve:
    """Affine arithmetic on y^2 = x^3 + 7 over F_q (points are tuples, None = infinity)."""

    def __init__(self, q, r, h):
        self.q, self.r, self.h = q, r, h

    def on_curve(self, P):
        return P is None or (P[1] * P[1] - P[0] ** 3 - CB) % self.q == 0

    def add(self, P, Q):
        q = self.q
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % q == 0:
                return None
            l = 3 * x1 * x1 * pow(2 * y1, -1, q) % q
        else:
            l = (y2 - y1) * pow(x2 - x1, -1, q) % q
        x3 = (l * l - x1 - x2) % q
        return (x3, (l * (x1 - x3) - y1) % q)

    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.q)

    def mul(self, k, P):
        k %= self.r
        Rr = None
        while k:
            if k & 1:
                Rr = self.add(Rr, P)
            P = self.add(P, P)
            k >>= 1
        return Rr

    def fixed_base(self, P):
        """Precompute j * 2^(8i) * P for a fixed base P (toy curves only)."""
        self.fb = []
        Q = P
        for i in range((self.r.bit_length() + 7) // 8):
            row = [None]
            for j in range(1, 256):
                row.append(self.add(row[-1], Q))
            self.fb.append(row)
            Q = self.add(row[-1], Q)

    def fmul(self, k):
        """[k]G via the fixed-base table: at most 4 additions on the toy curves."""
        k %= self.r
        Rr, i = None, 0
        while k:
            Rr = self.add(Rr, self.fb[i][k & 255])
            k >>= 8
            i += 1
        return Rr

    def phi(self, P, j=1):
        if P is None:
            return None
        return (P[0] * pow(self.beta, j % 3, self.q) % self.q, P[1])


def sqrt_mod(a, q):
    a %= q
    if a == 0:
        return 0
    if pow(a, (q - 1) // 2, q) != 1:
        return None
    s, t = 0, q - 1
    while t % 2 == 0:
        t //= 2
        s += 1
    z = 2
    while pow(z, (q - 1) // 2, q) != q - 1:
        z += 1
    m, c, tt, x = s, pow(z, t, q), pow(a, t, q), pow(a, (t + 1) // 2, q)
    while tt != 1:
        i, t2 = 1, tt * tt % q
        while t2 != 1:
            t2 = t2 * t2 % q
            i += 1
        b = pow(c, 1 << (m - i - 1), q)
        m, c, tt, x = i, b * b % q, tt * b * b % q, x * b % q
    return x


def cube_roots_of_unity(m):
    for g in range(2, 1000):
        c = pow(g, (m - 1) // 3, m)
        if c != 1:
            return [c, c * c % m]


def make_toy(name):
    q, n, h = TOYS[name]
    r = n // h
    E = Curve(q, r, h)
    while True:
        x = rng.randrange(q)
        y = sqrt_mod(x ** 3 + CB, q)
        if y is None:
            continue
        G = None  # cofactor multiple [h](x, y), computed without reducing h mod r
        P, k = (x, y), h
        while k:
            if k & 1:
                G = E.add(G, P)
            P = E.add(P, P)
            k >>= 1
        if G is not None and E.mul(r, G) is None:
            break
    E.G = G
    E.fixed_base(G)
    for b in cube_roots_of_unity(q):
        for l in cube_roots_of_unity(r):
            if E.mul(l, G) == (G[0] * b % q, G[1]):
                E.beta, E.lam = b, l
    return E


# ---------------------------------------------------------------- the class walk
class Walker:
    """Exactly the walk of proof.md section 2 (production rules), on curve E."""

    def __init__(self, E, R, dp_bits, T, rng, start_point=None):
        self.E, self.R, self.dp_bits = E, R, dp_bits
        self.rbits = R.bit_length() - 1
        assert 1 << self.rbits == R
        self.dpmask = (1 << dp_bits) - 1
        self.Lmax = 20 << dp_bits
        r = E.r
        lam = E.lam
        self.units = [(-1) ** s * pow(lam, j, r) % r for s in (0, 1) for j in (0, 1, 2)]
        self.uinv = [pow(u, -1, r) for u in self.units]
        self.T = T
        self.table = []  # table[h][t] = (point u_t^{-1} M_h, coeff c, coeff d)
        while len(self.table) < R:
            c, d = rng.randrange(r), rng.randrange(r)
            M = start_point(c, d) if start_point else E.add(E.mul(c, E.G), E.mul(d, T))
            if M is None:  # M_h = O: redraw (probability 1/r)
                continue
            row = []
            for t in range(6):
                s, j = divmod(t, 3)
                Mt = E.phi(M, -j)
                if s:
                    Mt = E.neg(Mt)
                ui = self.uinv[t]
                row.append((Mt, ui * c % r, ui * d % r))
            self.table.append(row)
        self.escape_add = False
        self.esc_row = None
        self.stats = dict(steps=0, doublings=0, esc2=0, esccp=0, abandoned=0, exceptional=0)

    def escape(self, P, a, b, track):
        """Escape from the class of P: doubling (the attack). escape_add = True is a
        CONTROL only: an additive escape P + [u_t^{-1}]M_esc, which proof.md rejects."""
        r = self.E.r
        self.stats["doublings"] += 1
        if self.escape_add:
            t = self.canon(P)[1]
            TP, tc, td = self.esc_row[t]
            Q = self.add_checked(P, TP)
            if track:
                a, b = (a + tc) % r, (b + td) % r
            return Q, a, b
        Q = self.E.add(P, P)  # doubling of a G1 point never hits y = 0 (odd order)
        if track:
            a, b = 2 * a % r, 2 * b % r
        return Q, a, b

    def add_checked(self, P, Q):
        """Incomplete affine addition: None signals an exceptional operand pair
        (equal x-coordinates, i.e. P = +-Q), on which the trail is abandoned."""
        if P[0] == Q[0]:
            return None
        return self.E.add(P, Q)

    def canon(self, P):
        """Return (K, t, canonical point): canonical = units[t] * P."""
        q, beta = self.E.q, self.E.beta
        x, y = P
        if CLS == 1:  # control: point itself, key injective on points
            return x + q * y, 0, P
        if CLS == 2:  # control: negation classes, key x
            if q - y < y:
                return x, 3, (x, q - y)
            return x, 0, P
        bx = beta * x % q
        b2x = (-x - bx) % q  # beta^2 x = -(1+beta) x
        j = 0
        cx = x
        if bx < cx:
            cx, j = bx, 1
        if b2x < cx:
            cx, j = b2x, 2
        ny = q - y
        if ny < y:
            return ny, 3 + j, (cx, ny)
        return y, j, (cx, y)

    def step_target(self, P, K, t):
        h = K & (self.R - 1)
        return self.table[h][t]

    def is_dp(self, K):
        return ((K >> self.rbits) & self.dpmask) == 0

    def run_trail(self, P, coeffs=None, on_class=None):
        """Walk from P until a distinguished class or the trail cap.

        Returns (status, K, length, P_end, coeffs_end). If coeffs=(a,b) is given,
        coefficients with P = aG + bT are tracked (used only when re-walking).
        on_class(K) is called on every class arrival (instrumentation only)."""
        E, r = self.E, self.E.r
        k1 = k2 = None  # keys one and two arrivals back
        cp = None  # checkpoint key
        length = 0
        pending_double = False
        a = b = None
        if coeffs is not None:
            a, b = coeffs
        while True:
            K, t, c = self.canon(P)
            if on_class is not None:
                on_class(K)
            if self.is_dp(K):
                if coeffs is not None:
                    u = self.units[t]
                    return "dp", K, length, c, (u * a % r, u * b % r)
                return "dp", K, length, c, None
            if length >= self.Lmax:
                self.stats["abandoned"] += 1
                return "cap", K, length, c, None
            if pending_double:
                P, a, b = self.escape(P, a, b, coeffs is not None)
                if P is None:
                    self.stats["exceptional"] += 1
                    return "exc", K, length, None, None
                pending_double = False
                k1 = k2 = cp = None
                length += 1
                continue
            if k2 is not None and K == k2:
                self.stats["esc2"] += 1
                if K < k1:
                    P, a, b = self.escape(P, a, b, coeffs is not None)
                    if P is None:
                        self.stats["exceptional"] += 1
                        return "exc", K, length, None, None
                    k1 = k2 = cp = None
                    length += 1
                    continue
                pending_double = True  # take one ordinary step to the smaller class, then double
            elif cp is not None and K == cp:
                self.stats["esccp"] += 1
                best = (K, P, a, b)
                Q, qa, qb = P, a, b
                trav = 0
                while True:
                    trav += 1
                    if trav > 64:  # not a cycle of length <= 64: abandon the trail
                        self.stats["abandoned"] += 1
                        return "trav", K, length, None, None
                    TP, tc, td = self.step_target(Q, *self.canon(Q)[:2])
                    Q = self.add_checked(Q, TP)
                    if Q is None:
                        self.stats["exceptional"] += 1
                        return "exc", K, length, None, None
                    if coeffs is not None:
                        qa, qb = (qa + tc) % r, (qb + td) % r
                    length += 1
                    self.stats["steps"] += 1
                    KQ = self.canon(Q)[0]
                    if on_class is not None:
                        on_class(KQ)  # traversal arrivals are counted visits too
                    if KQ == K:
                        break
                    if KQ < best[0]:
                        best = (KQ, Q, qa, qb)
                _, P, a, b = best
                P, a, b = self.escape(P, a, b, coeffs is not None)
                if P is None:
                    self.stats["exceptional"] += 1
                    return "exc", K, length, None, None
                k1 = k2 = cp = None
                length += 1
                continue
            if length % 64 == 0:
                cp = K
            TP, tc, td = self.step_target(P, K, t)
            P = self.add_checked(P, TP)
            if P is None:
                self.stats["exceptional"] += 1
                return "exc", K, length, None, None
            if coeffs is not None:
                a, b = (a + tc) % r, (b + td) % r
            self.stats["steps"] += 1
            length += 1
            k2, k1 = k1, K


# ---------------------------------------------------------------- fast canonical index (production)
M64 = (1 << 64) - 1
GQ = (PROD["beta"] << 384) // PROD["q"]  # floor(beta/p * 2^384): 6 limbs
GL = [(GQ >> (64 * j)) & M64 for j in range(6)]
TOL = 40


def exact_j(x):
    q, b = PROD["q"], PROD["beta"]
    v = [x, b * x % q, b * b % q * x % q]
    return min(range(3), key=lambda j: v[j])


def fast_j(x):
    """Return j, or None for 'fall back to the exact rule'. Only 8 limb products are
    formed: the low halves (mul) of the four column-5 products x_i*G_{5-i} and the high
    halves (umulh) of the four column-4 products x_i*G_{4-i}; their sum mod 2^64 is F1.
    The true bits 320..383 of x*G, i.e. floor(2^64 frac(x beta/p)) up to < 1, lie in
    [F1, F1 + 14]. F0 = top limb of x approximates 2^64 x/p within [0, 2). With
    x + beta x + beta^2 x = p or 2p, F2 = 2^64 - F0 - F1 or 2^65 - F0 - F1."""
    xl = [(x >> (64 * i)) & M64 for i in range(4)]
    F1 = 0
    for i in range(4):
        F1 += (xl[i] * GL[5 - i]) & M64
        F1 += (xl[i] * GL[4 - i]) >> 64
    F1 &= M64
    if F1 < TOL or F1 > M64 - TOL:  # beta x mod p within ~2^-58 p of 0 or p
        return None
    F0 = xl[3]
    S = F0 + F1
    if abs(S - (1 << 64)) <= TOL:  # cannot decide whether x + beta x wraps mod p
        return None
    F2 = (1 << 64) - S if S < (1 << 64) else (1 << 65) - S
    vals = [F0, F1, F2]
    jm = min(range(3), key=lambda j: vals[j])
    for j in range(3):
        if j != jm and vals[j] - vals[jm] <= TOL:
            return None
    return jm


def check_fast_j(n_random, n_stress):
    q, b = PROD["q"], PROD["beta"]
    rng = random.Random(seed + 10 ** 7)  # separate stream: Part B is unaffected by this check
    res = dict(random=n_random, random_fallbacks=0, random_mismatches=0,
               stress=0, stress_fallbacks=0, stress_mismatches=0)
    for _ in range(n_random):
        if time.time() - T0 > budget:  # stay inside the runner's time limit
            break
        res["random_done"] = res.get("random_done", 0) + 1
        x = rng.randrange(1, q)
        f = fast_j(x)
        if f is None:
            res["random_fallbacks"] += 1
        elif f != exact_j(x):
            res["random_mismatches"] += 1
    mults = [b - 1, b * b - 1, b - b * b, 1 + b, b]
    for _ in range(n_stress):
        if time.time() - T0 > budget:
            break
        k = rng.randrange(150, 250)
        d = rng.randrange(-(1 << k), 1 << k)
        mlt = mults[rng.randrange(len(mults))]
        x = d * pow(mlt % q, -1, q) % q
        if mlt == 1 + b:
            x = (q - d) * pow(mlt % q, -1, q) % q
        if x == 0:
            continue
        res["stress"] += 1
        f = fast_j(x)
        if f is None:
            res["stress_fallbacks"] += 1
        elif f != exact_j(x):
            res["stress_mismatches"] += 1
    return res


# ---------------------------------------------------------------- per-step sequence simulator
# Section 5 prices each per-step row as the exact length of a written-out sequence.
# seq_rows() holds the rows as lists of ARM-like operations, SEQ_SEMANTICS gives every
# operation's exact semantics, and check_sequences() executes each row against fast_j
# on random production inputs and on boundary inputs at its own branch, reporting the
# length of the sequence it ran and its mismatch count.
MASK64 = M64
INDEX_MASK = (1 << 20) - 1          # h = K mod R, R = 2^20
TWO64_M81 = M64 - 80                # 2^64 - 81: the register the range test compares with
I_INSTR = 1 << 17                   # I < 2^17, section 5's inversion bound
M_INSTR = 227                       # M <= 227, section 5's multiplication bound
BATCH_MIN = 1 << 15                 # section 3: every walk this table prices is in a batch of L >= 2^15
BATCH_MAX = 1 << 16

# Condition codes, as the flags (N, Z, C) each one tests; `lo` is `cc` and `hs` is `cs`.
SEQ_CONDS = {
    "eq": lambda n, z, c: z,
    "ne": lambda n, z, c: not z,
    "hs": lambda n, z, c: c,
    "lo": lambda n, z, c: not c,
    "hi": lambda n, z, c: c and not z,
    "ls": lambda n, z, c: not (c and not z),
}
SEQ_CONDS["cs"] = SEQ_CONDS["hs"]
SEQ_CONDS["cc"] = SEQ_CONDS["lo"]

SEQ_SEMANTICS = {
    "add": "d = (a + b) mod 2^64",
    "adds": "d = (a + b) mod 2^64; flags from the true sum: C = carry out of bit 63, N and Z from d",
    "adc": "d = (a + b + C) mod 2^64",
    "add_lsl6": "d = (a + (b << 6)) mod 2^64",
    "sub": "d = (a - b) mod 2^64",
    "subs": "d = (a - b) mod 2^64; C = NOT borrow, so C = 1 exactly when a >= b, N and Z from d",
    "sbcs": "d = (a - b - (1 - C)) mod 2^64, with the flags of that true difference",
    "neg": "d = (-a) mod 2^64",
    "and": "d = a & b",
    "mul": "d = (a * b) mod 2^64, the low half of the product",
    "madd": "d = (a + b * c) mod 2^64",
    "umulh": "d = floor(a * b / 2^64), the high half of the 128-bit product",
    "cmp": "the flags of (a - b), no destination",
    "tst": "the flags of (a & b), no destination",
    "ccmp": "if cond holds: the flags of (a - b); otherwise the immediate's flags",
    "csel": "d = a if cond else b",
    "csinc": "d = a if cond else b + 1",
    "csneg": "d = (-b) mod 2^64 if cond else a",
    "cset": "d = 1 if cond else 0",
    "ldp": "d1 = M[base + off], d2 = M[base + off + 8]; a post-index adds off to base",
    "stp": "M[base + off] = a, M[base + off + 8] = b; a post-index adds off to base",
    "ldr": "d = M[base + off]; a post-index adds off to base",
    "str": "M[base + off] = a; a post-index adds off to base",
    "b.cond": "if cond holds: execution continues at the target named",
    "cbz": "if a == 0: execution continues at the target named",
    "cbnz": "if a != 0: execution continues at the target named",
}


def seq_expand(seq):
    """Expand a row's sequence into a straight-line Python function.

    Each operation becomes the one statement SEQ_SEMANTICS describes, with registers as
    local variables, so this is a transcription of the row rather than a second
    interpretation.  It exists because 95 operations times 200,000 inputs must fit the
    runner's 20 s; the interpreter above is the definition and check_sequences requires
    the two to agree register for register on a sample and on every boundary input."""
    reads, outs, defined, lines = [], [], [], []

    def rd(*names):
        for x in names:
            if x is not None and x not in outs and x not in reads:
                reads.append(x)

    def wr(*names):
        for x in names:
            if x not in outs:
                outs.append(x)
            if x not in defined:
                defined.append(x)

    def conds(c):
        return SEQ_EXPR[c]

    def exit_parts():
        writes = " ; ".join("regs[%r] = %s" % (x, x) for x in defined)
        return (writes + " ; ") if writes else ""

    def emit(statement, written=()):
        lines.append("    " + statement)
        for x in written:
            wr(x)

    for mn, r1, r2, r3, ex in seq:
        if mn == "csel":
            rd(r2, r3)
            emit("%s = %s if %s else %s" % (r1, r2, conds(ex), r3), [r1])
        elif mn == "csneg":
            rd(r2, r3)
            emit("%s = ((-%s) & MASK64) if %s else %s" % (r1, r3, conds(ex), r2), [r1])
        elif mn == "csinc":
            rd(r2, r3)
            emit("%s = %s if %s else %s + 1" % (r1, r2, conds(ex), r3), [r1])
        elif mn == "cset":
            emit("%s = 1 if %s else 0" % (r1, conds(r2)), [r1])
        elif mn == "cmp":
            rd(r2, r3)
            emit("_d = (%s - %s) & MASK64; n = _d >> 63; z = int(_d == 0); c = int(%s >= %s)"
                 % (r2, r3, r2, r3))
        elif mn == "tst":
            rd(r2, r3)
            emit("_d = %s & %s; n = _d >> 63; z = int(_d == 0); c = 0" % (r2, r3))
        elif mn == "ccmp":
            imm, cond = ex
            rd(r2, r3)
            emit("if %s:" % conds(cond))
            emit("    _d = (%s - %s) & MASK64; n = _d >> 63; z = int(_d == 0); c = int(%s >= %s)"
                 % (r2, r3, r2, r3))
            emit("else:")
            emit("    n = %d; z = %d; c = %d" % ((imm >> 3) & 1, (imm >> 2) & 1, (imm >> 1) & 1))
        elif mn == "b.cond":
            cond, label, leave = ex
            if leave:
                emit("if %s: %sreturn fired + [%r]" % (conds(cond), exit_parts(), label))
            else:
                emit("if %s: fired.append(%r)" % (conds(cond), label))
        elif mn == "cbz" or mn == "cbnz":
            label, leave = ex
            rd(r2)
            test = "%s == 0" % r2 if mn == "cbz" else "%s != 0" % r2
            if leave:
                emit("if %s: %sreturn fired + [%r]" % (test, exit_parts(), label))
            else:
                emit("if %s: fired.append(%r)" % (test, label))
        elif mn in ("ldp", "ldr", "stp", "str"):
            mode, off = ex
            base = r3
            rd(base)
            if mn in ("ldp", "stp"):
                rd(r1, r2) if mn == "stp" else None
            else:
                rd(r1) if mn == "str" else None
            if mode == "post":       # a post-index reaches memory first, then adds
                at, reach = "_b", 0
                lines.append("    _b = %s" % base)
            else:
                at, reach = base, off
            if mn == "ldp":
                emit("%s = mem[%s + %d]; %s = mem[%s + %d]"
                     % (r1, at, reach, r2, at, reach + 8), [r1, r2])
            elif mn == "ldr":
                emit("%s = mem[%s + %d]" % (r1, at, reach), [r1])
            elif mn == "stp":
                emit("mem[%s + %d] = %s; mem[%s + %d] = %s"
                     % (at, reach, r1, at, reach + 8, r2))
            else:
                emit("mem[%s + %d] = %s" % (at, reach, r1))
            if mode == "post":
                emit("%s = %s + %d" % (r3, "_b", off), [r3])
        elif mn == "add":
            rd(r2, r3)
            emit("%s = (%s + %s) & MASK64" % (r1, r2, r3), [r1])
        elif mn == "add_lsl6":
            rd(r2, r3)
            emit("%s = (%s + (%s << 6)) & MASK64" % (r1, r2, r3), [r1])
        elif mn == "adds":
            rd(r2, r3)
            emit("_s = %s + %s; c = int(_s > MASK64); %s = _s & MASK64; n = %s >> 63; z = int(%s == 0)"
                 % (r2, r3, r1, r1, r1), [r1])
        elif mn == "adc":
            rd(r2, r3)
            emit("%s = (%s + %s + c) & MASK64" % (r1, r2, r3), [r1])
        elif mn == "sub":
            rd(r2, r3)
            emit("%s = (%s - %s) & MASK64" % (r1, r2, r3), [r1])
        elif mn == "subs":
            rd(r2, r3)
            emit("_d = (%s - %s) & MASK64; c = int(%s >= %s); n = _d >> 63; z = int(_d == 0); %s = _d"
                 % (r2, r3, r2, r3, r1), [r1])
        elif mn == "sbcs":
            rd(r2, r3)
            emit("_d = (%s - %s - 1 + c) & MASK64; c = int(%s >= %s + 1 - c); n = _d >> 63; "
                 "z = int(_d == 0); %s = _d" % (r2, r3, r2, r3, r1), [r1])
        elif mn == "neg":
            rd(r2)
            emit("%s = (-%s) & MASK64" % (r1, r2), [r1])
        elif mn == "mul":
            rd(r2, r3)
            emit("%s = (%s * %s) & MASK64" % (r1, r2, r3), [r1])
        elif mn == "madd":
            rd(ex, r2, r3)
            emit("%s = (%s + %s * %s) & MASK64" % (r1, ex, r2, r3), [r1])
        elif mn == "umulh":
            rd(r2, r3)
            emit("%s = (%s * %s) >> 64" % (r1, r2, r3), [r1])
        elif mn == "and":
            rd(r2, r3)
            emit("%s = %s & %s" % (r1, r2, r3), [r1])
        else:
            raise AssertionError("no expansion for " + mn)
    lines.append("    " + exit_parts() + "return fired")
    src = ["def _f(regs, mem):",
           "    " + ("; ".join("%s = regs[%r]" % (x, x) for x in reads) if reads else "pass"),
           "    n = z = c = 0",
           "    fired = []"]
    src += lines
    ns = {"MASK64": MASK64, "SEQ_EXPR": SEQ_EXPR}
    exec(compile("\n".join(src), "<seq>", "exec"), ns)
    fn = ns["_f"]
    fn.source = "\n".join(src)
    fn.reads, fn.outs = reads, outs
    return fn


SEQ_EXPR = {"eq": "z", "ne": "not z", "hs": "c", "cs": "c", "lo": "not c", "cc": "not c",
            "hi": "c and not z", "ls": "not (c and not z)"}


def seq_imm(x):
    if isinstance(x, int):
        return "#0x%X" % x if x > (1 << 32) else "#%d" % x
    return "#" + x[1:] if x[:1] == "k" and x[1:].isdigit() else x


def seq_mem(base, ex):
    mode, off = ex
    return "[%s], #%d" % (base, off) if mode == "post" else "[%s, #%d]" % (base, off)


def seq_text(seq):
    """Render a row's sequence the way section 5 writes it out."""
    out = []
    for mn, r1, r2, r3, ex in seq:
        if mn in ("csel", "csinc", "csneg"):
            out.append("%s %s, %s, %s, %s" % (mn, r1, r2, r3, ex))
        elif mn == "cset":
            out.append("cset %s, %s" % (r1, r2))
        elif mn == "neg":
            out.append("neg %s, %s" % (r1, r2))
        elif mn in ("cmp", "tst"):
            out.append("%s %s, %s" % (mn, r2, seq_imm(r3)))
        elif mn == "ccmp":
            out.append("ccmp %s, %s, #%d, %s" % (r2, r3, ex[0], ex[1]))
        elif mn == "b.cond":
            out.append("b.%s %s" % (ex[0], ex[1]))
        elif mn in ("cbz", "cbnz"):
            out.append("%s %s, %s" % (mn, r2, ex[0]))
        elif mn == "ldp":
            out.append("ldp %s, %s, %s" % (r1, r2, seq_mem(r3, ex)))
        elif mn == "ldr":
            out.append("ldr %s, %s" % (r1, seq_mem(r3, ex)))
        elif mn == "stp":
            out.append("stp %s, %s, %s" % (r1, r2, seq_mem(r3, ex)))
        elif mn == "str":
            out.append("str %s, %s" % (r1, seq_mem(r3, ex)))
        elif mn == "madd":
            out.append("madd %s, %s, %s, %s" % (r1, r2, seq_imm(r3), ex))
        elif mn == "add_lsl6":
            out.append("add %s, %s, %s, lsl #6" % (r1, r2, r3))
        elif ex is None:
            out.append("%s %s, %s, %s" % (mn, r1, r2, seq_imm(r3)))
        else:
            raise AssertionError("no rendering for " + mn)
    return " ; ".join(out)


def seq_rows():
    """Every priced per-step row of section 5, as an executable sequence.

    The values the proof keeps in registers across the batch loop (p0, the all-ones
    word, 2^64 - 81, 2, 3, 6, the meter price, the table base) are registers here and an
    inline literal is an immediate; both are loaded once, inside setup's 2^24 reserve."""
    f1 = [("ldp", "G1", "G2", "gmem", ("off", 0)),
          ("ldp", "G3", "G4", "gmem", ("off", 16)),
          ("ldr", "G5", None, "gmem", ("off", 32)),
          ("mul", "F1", "x0", "G5", None),
          ("madd", "F1", "x1", "G4", "F1"),
          ("madd", "F1", "x2", "G3", "F1"),
          ("madd", "F1", "x3", "G2", "F1"),
          ("umulh", "t0", "x0", "G4", None), ("add", "F1", "F1", "t0", None),
          ("umulh", "t0", "x1", "G3", None), ("add", "F1", "F1", "t0", None),
          ("umulh", "t0", "x2", "G2", None), ("add", "F1", "F1", "t0", None),
          ("umulh", "t0", "x3", "G1", None), ("add", "F1", "F1", "t0", None)]
    wrap = [("sub", "t0", "F1", 40, None),
            ("cmp", None, "t0", "m64m81", None),
            ("b.cond", None, None, None, ("hi", "range", True)),
            ("adds", "S", "F0", "F1", None),
            ("csneg", "T", "S", "S", "cc"),
            ("cmp", None, "T", 40, None),
            ("b.cond", None, None, None, ("ls", "wrap", True)),
            ("neg", "F2", "S", None, None)]
    decide = [("cmp", None, "F0", "F1", None),
              ("csel", "m1", "F0", "F1", "lo"),
              ("csel", "M1", "F1", "F0", "lo"),
              ("cset", "jm", "hs", None, None),
              ("cmp", None, "m1", "F2", None),
              ("csel", "mn", "m1", "F2", "lo"),
              ("csel", "u1", "F2", "m1", "lo"),
              ("csel", "j", "two", "jm", "hi"),
              ("cmp", None, "M1", "u1", None),
              ("csel", "mid", "M1", "u1", "lo"),
              ("sub", "gap", "mid", "mn", None),
              ("cmp", None, "gap", 40, None),
              ("b.cond", None, None, None, ("ls", "gap", True))]
    sign = [("subs", "n0", "p0", "y0", None),
            ("sbcs", "n1", "ones", "y1", None),
            ("sbcs", "n2", "ones", "y2", None),
            ("sbcs", "n3", "ones", "y3", None),
            ("subs", "t0", "n0", "y0", None),
            ("sbcs", "t1", "n1", "y1", None),
            ("sbcs", "t2", "n2", "y2", None),
            ("sbcs", "t3", "n3", "y3", None),
            ("cset", "s", "cc", None, None),
            ("csel", "K0", "n0", "y0", "cc"),
            ("csel", "K1", "n1", "y1", "cc"),
            ("csel", "K2", "n2", "y2", "cc"),
            ("csel", "K3", "n3", "y3", "cc")]
    route = [("and", "h", "K0", INDEX_MASK, None),
             ("cbz", None, "K1", None, ("dp", True)),
             ("madd", "t", "s", "three", "j"),
             ("madd", "w", "h", "six", "t"),
             ("add_lsl6", "xT", "base", "w", None)]
    traffic = [("ldp", "px0", "px1", "rec", ("off", 0)),
               ("ldp", "px2", "px3", "rec", ("off", 16)),
               ("ldp", "py0", "py1", "rec", ("off", 32)),
               ("ldp", "py2", "py3", "rec", ("off", 48)),
               ("ldp", "c0", "c1", "xT", ("off", 0)),
               ("ldp", "c2", "c3", "xT", ("off", 16)),
               ("ldp", "c4", "c5", "xT", ("off", 32)),
               ("ldp", "c6", "c7", "xT", ("off", 48)),
               ("stp", "nx0", "nx1", "rec", ("off", 0)),
               ("stp", "nx2", "nx3", "rec", ("off", 16)),
               ("stp", "ny0", "ny1", "rec", ("off", 32)),
               ("stp", "ny2", "ny3", "rec", ("off", 48))]
    hist = [("cbnz", None, "flag", None, ("esc", True)),
            ("ldp", "f0", "f1", "hp", ("post", 32)),
            ("ldr", "f2", None, "hp", ("off", -16)),
            ("cmp", None, "f0", "K0", None),
            ("ccmp", None, "f1", "K2", (0, "eq")),
            ("ccmp", None, "f2", "K3", (0, "eq")),
            ("b.cond", None, None, None, ("eq", "cycle2", True)),
            ("stp", "K0", "K2", "hp", ("off", -32)),
            ("str", "K3", None, "hp", ("off", -16)),
            ("ldp", "g0", "g1", "cp", ("off", 0)),
            ("ldr", "g2", None, "cp", ("off", 16)),
            ("cmp", None, "g0", "K0", None),
            ("ccmp", None, "g1", "K2", (0, "eq")),
            ("ccmp", None, "g2", "K3", (0, "eq")),
            ("b.cond", None, None, None, ("eq", "checkpoint", True)),
            ("tst", None, "len0", 63, None),
            ("b.cond", None, None, None, ("ne", "skip_cp", True)),
            ("stp", "K0", "K2", "cp", ("off", 0)),
            ("str", "K3", None, "cp", ("off", 16))]
    meter = [("adds", "len0", "len0", "one", None),
             ("adc", "len1", "len1", "zr", None),
             ("cmp", None, "len1", "twenty", None),
             ("b.cond", None, None, None, ("hs", "cap", True)),
             ("subs", "mtr", "mtr", "one", None),
             ("add", "ops", "ops", "price", None),
             ("b.cond", None, None, None, ("lo", "refill", False)),
             ("ldr", "qp", None, "qp", ("post", 8)),
             ("subs", "cnt", "cnt", "one", None),
             ("b.cond", None, None, None, ("ne", "loop", True))]
    batch_price = -(-(I_INSTR - 3 * M_INSTR + 4) // BATCH_MIN)   # ceil of section 5's bound
    return [
        dict(key="F1", row="1.1 F1: 8 limb products, 3 loads of G", ops=f1),
        dict(key="wrap", row="1.1 range test, S, wrap test, F2", ops=wrap),
        dict(key="decide", row="1.1 decision: min, max, argmin, median, gap", ops=decide),
        dict(key="sign", row="sign, p - y, select limbs of K", ops=sign),
        dict(key="route", row="h, distinguished test, t, table address", ops=route),
        dict(key="traffic", row="loads of table point and state, store of state", ops=traffic),
        dict(key="hist", row="fingerprint tests, history and checkpoint stores", ops=hist),
        dict(key="zero", row="zero-denominator test, folded into the batch", ops=[]),
        dict(key="meter", row="length, cap test, meter, queue and loop", ops=meter),
        dict(key="batch", row="batch inversion: I/L - 3M/L < I/2^15", ops=[],
             price=batch_price, text="(an upper bound, not a sequence: checked over every L)"),
    ]


def seq_immify(seq):
    """Move every inline literal into a constant register (k40 holds 40), so a row spends
    no instruction on one and the interpreter sees register operands only."""
    out = []
    for mn, r1, r2, r3, ex in seq:
        if isinstance(r3, int):
            r3 = "k%d" % r3
        out.append((mn, r1, r2, r3, ex))
    return out


SEQ_ROWS = seq_rows()
for _r in SEQ_ROWS:
    _r["ops"] = seq_immify(_r["ops"])
SEQ_OPS = dict((r["key"], r["ops"]) for r in SEQ_ROWS)
SEQ_CONSTS = dict((_c, int(_c[1:])) for _r in SEQ_ROWS for (_mn, _a, _b, _c, _e) in _r["ops"]
                  if isinstance(_c, str) and _c[:1] == "k" and _c[1:].isdigit())
SEQ_FAST = dict((_r["key"], seq_expand(_r["ops"])) for _r in SEQ_ROWS if _r["ops"])
# Every mnemonic a row uses must have a stated semantics: nothing runs undocumented.
SEQ_UNDOCUMENTED = sorted(set(_mn for _r in SEQ_ROWS for (_mn, _a, _b, _c, _e) in _r["ops"]
                               if _mn not in SEQ_SEMANTICS))


SEQ_ORDER = ("F1", "wrap", "decide", "sign", "route", "traffic", "hist", "zero", "meter", "batch")


def _limbs(v):
    return [(v >> (64 * i)) & MASK64 for i in range(4)]


def _fp(k):
    """The package's fingerprint: limbs 0, 2 and 3 of the class key."""
    return [(k >> 0) & MASK64, (k >> 128) & MASK64, (k >> 192) & MASK64]


def _ref_f1(xl):
    """F1 as fast_j's own lines form it: eight limb products summed mod 2^64."""
    F1 = 0
    for i in range(4):
        F1 += (xl[i] * GL[5 - i]) & MASK64
        F1 += (xl[i] * GL[4 - i]) >> 64
    return F1 & MASK64


def _ref_decide(F0, F1, F2):
    """fast_j's decision: the argmin and the second-smallest-minus-smallest gap test."""
    if F0 < F1:
        m1, M1, jm = F0, F1, 0
    else:
        m1, M1, jm = F1, F0, 1
    if m1 < F2:
        mn, u1 = m1, F2
    else:
        mn, u1 = F2, m1
    j = 2 if m1 > F2 else jm
    mid = sorted((F0, F1, F2))[1]
    return dict(m1=m1, M1=M1, jm=jm, mn=mn, u1=u1, j=j, mid=mid, gap=(mid - mn) & MASK64)


def _chk_f1(xl, regs, mem):
    """F1 from the eight limb products of x's limbs and G's five limbs."""
    regs["x0"], regs["x1"], regs["x2"], regs["x3"] = xl
    regs["F0"] = xl[3]
    SEQ_FAST["F1"](regs, mem)
    return int(regs["F1"] != _ref_f1(xl))


def _chk_wrap(regs, F0, F1):
    """The range test, S, the wrap distance and F2, on one (F0, F1) pair."""
    regs["F0"], regs["F1"] = F0, F1
    fired = SEQ_FAST["wrap"](regs, {})
    in_range = F1 < TOL or F1 > M64 - TOL
    bad = int(("range" in fired) != in_range)
    if in_range:
        return bad, True
    S = F0 + F1
    wrapped = abs(S - (1 << 64)) <= TOL
    bad += int(regs["S"] != (S & MASK64))
    bad += int(regs["T"] != abs(S - (1 << 64)))
    bad += int(("wrap" in fired) != wrapped)
    if not wrapped:
        bad += int(regs["F2"] != ((-S) & MASK64))
    return bad, wrapped


def _chk_decide(regs, F0, F1, F2):
    """min, max, argmin, median, gap and the gap test, on one triple."""
    regs["F0"], regs["F1"], regs["F2"] = F0, F1, F2
    fired = SEQ_FAST["decide"](regs, {})
    ref = _ref_decide(F0, F1, F2)
    got = (regs["m1"], regs["M1"], regs["jm"], regs["mn"], regs["u1"], regs["j"],
           regs["mid"], regs["gap"])
    want = (ref["m1"], ref["M1"], ref["jm"], ref["mn"], ref["u1"], ref["j"],
            ref["mid"], ref["gap"])
    bad = int(got != want) + int(("gap" in fired) != (ref["gap"] <= TOL))
    return bad, bool("gap" in fired)


def _chk_sign(yl, regs):
    """p - y, the cost-model sign and the limbs of the key K = min(y, p - y)."""
    p = PROD["q"]
    regs["y0"], regs["y1"], regs["y2"], regs["y3"] = yl
    regs["p0"], regs["ones"] = p & MASK64, MASK64
    SEQ_FAST["sign"](regs, {})
    y = yl[0] | (yl[1] << 64) | (yl[2] << 128) | (yl[3] << 192)
    ny = p - y
    neg = 1 if ny < y else 0
    got = (regs["n0"], regs["n1"], regs["n2"], regs["n3"], regs["s"],
           regs["K0"], regs["K1"], regs["K2"], regs["K3"])
    want = tuple(_limbs(ny) + [neg] + _limbs(ny if neg else y))
    return int(got != want)


def _chk_route(regs, kl, K, s, j):
    """h = K mod 2^20, the distinguished test, t = 3s + j and the table address."""
    regs["K0"], regs["K1"], regs["K2"], regs["K3"] = kl
    regs["s"], regs["j"] = s, j
    regs["two"], regs["three"], regs["six"], regs["base"] = 2, 3, 6, 0x100000
    fired = SEQ_FAST["route"](regs, {})
    h = K & INDEX_MASK
    dp = ((K >> 64) & MASK64) == 0      # section 2: limb 1 of K is the distinguished field
    bad = int(regs["h"] != h) + int(("dp" in fired) != dp)
    if not dp:
        w = 6 * h + 3 * s + j
        bad += int((regs["t"], regs["w"], regs["xT"]) != (3 * s + j, w, 0x100000 + (w << 6)))
    return bad, dp


_TBL = [(0xA5A50000 + 0x1111 * i) & MASK64 for i in range(8)]


def _chk_traffic(regs, xl, yl, x, y, xT):
    """The 12 moves: state and table loaded, the new point stored, no word dropped."""
    p = PROD["q"]
    nx, ny = (3 * x + 1) % p, (5 * y + 2) % p      # the priced reference addition's output
    nxl, nyl = _limbs(nx), _limbs(ny)
    tbl = _TBL
    mem = {}
    for i in range(4):
        mem[0x200000 + 8 * i] = xl[i]
        mem[0x200000 + 32 + 8 * i] = yl[i]
        mem[xT + 8 * i] = tbl[i]
        mem[xT + 32 + 8 * i] = tbl[4 + i]
    regs["rec"], regs["xT"] = 0x200000, xT
    regs["nx0"], regs["nx1"], regs["nx2"], regs["nx3"] = nxl
    regs["ny0"], regs["ny1"], regs["ny2"], regs["ny3"] = nyl
    SEQ_FAST["traffic"](regs, mem)
    got = (regs["px0"], regs["px1"], regs["px2"], regs["px3"],
           regs["py0"], regs["py1"], regs["py2"], regs["py3"],
           regs["c0"], regs["c1"], regs["c2"], regs["c3"],
           regs["c4"], regs["c5"], regs["c6"], regs["c7"])
    want = tuple(xl) + tuple(yl) + tuple(tbl)
    bad = int(got != want)
    got = (mem[0x200000], mem[0x200008], mem[0x200010], mem[0x200018],
           mem[0x200020], mem[0x200028], mem[0x200030], mem[0x200038])
    return bad + int(got != tuple(nxl) + tuple(nyl))


def _chk_hist(regs, K, kl, slot, cp, flag, length):
    """The pending-escape flag, the 2-cycle fingerprint test, the checkpoint test and
    the checkpoint write, on one history state."""
    mem = {}
    for i in range(3):
        mem[0x200040 + 8 * i] = slot[i]
        mem[0x200080 + 8 * i] = cp[i]
    regs["hp"], regs["cp"] = 0x200040, 0x200080
    regs["K0"], regs["K1"], regs["K2"], regs["K3"] = kl
    regs["flag"], regs["len0"] = flag, length & MASK64
    fired = SEQ_FAST["hist"](regs, mem)
    fp = _fp(K)
    if flag:
        return int(fired != ["esc"]), True
    cycle2 = list(slot) == fp
    bad = int(("cycle2" in fired) != cycle2)
    if cycle2:
        return bad, True
    bad += int((mem[0x200040], mem[0x200048], mem[0x200050]) != tuple(fp))
    hit = list(cp) == fp
    bad += int(("checkpoint" in fired) != hit)
    if hit:
        return bad, True
    write = (length % 64) == 0
    bad += int(("skip_cp" in fired) == write)
    want = fp if write else cp
    bad += int((mem[0x200080], mem[0x200088], mem[0x200090]) != tuple(want))
    return bad, False


def _chk_meter(regs, length, mtr, cnt):
    """The 2-limb length, the trail cap test, the meter block and the loop tail."""
    mem = {0x400000: 0x0123456789ABCDEF}
    regs["len1"], regs["len0"] = (length >> 64) & MASK64, length & MASK64
    regs["mtr"], regs["cnt"], regs["qp"], regs["ops"] = mtr, cnt, 0x400000, 0
    regs["zr"], regs["one"], regs["twenty"], regs["price"] = 0, 1, 20, 1099
    fired = SEQ_FAST["meter"](regs, mem)
    newlen = length + 1
    bad = int(((regs["len1"] << 64) | regs["len0"]) != newlen & ((1 << 128) - 1))
    cap = newlen >= 20 * (1 << 64)
    bad += int(("cap" in fired) != cap)
    if cap:
        return bad
    bad += int((regs["mtr"], regs["ops"], regs["qp"], regs["cnt"])
               != ((mtr - 1) & MASK64, 1099, 0x400008, (cnt - 1) & MASK64))
    bad += int(("refill" in fired) != (mtr == 0))
    bad += int(("loop" in fired) != (((cnt - 1) & MASK64) != 0))
    return bad


def _chk_batch(rng, L, zero_at):
    """The batch's zero-product test: the product of the L denominators is zero exactly
    when some member's x-coordinates are equal, and redoing the inversions one by one
    locates that member."""
    p = PROD["q"]
    xs = [rng.randrange(1, p) for _ in range(L)]
    ts = [rng.randrange(1, p) for _ in range(L)]
    if zero_at is not None:
        ts[zero_at] = xs[zero_at]
    prod = 1
    for i in range(L):
        prod = prod * ((xs[i] - ts[i]) % p) % p
    limbs = _limbs(prod)
    test = (limbs[0] | limbs[1]) | (limbs[2] | limbs[3])     # the three ORRs
    bad = int((test == 0) != (zero_at is not None))
    if zero_at is not None:                                   # the 600-unit locate action
        bad += int([i for i in range(L) if xs[i] == ts[i]] != [zero_at])
        bad += sum(1 for i in range(L) if i != zero_at and (xs[i] - ts[i]) % p == 0)
    return bad


def _chk_bound():
    """The batch-inversion row's price: (I - 3M)/L + 4/L <= 4 for every L in [2^15, 2^16]."""
    bad = 0
    for L in range(BATCH_MIN, BATCH_MAX + 1):
        bad += int((I_INSTR - 3 * M_INSTR) / L + 4.0 / L > 4)
    return bad


def check_sequences(n_random, deadline=None):
    """Execute every priced row of section 5; report its count and its mismatches.

    `deadline` (a time.time() value) stops the random loop early on a machine much
    slower than the reference one, so the declared experiment stays inside the runner's
    20 s; the count it reached is reported either way."""
    st = {}
    for r in SEQ_ROWS:
        st[r["key"]] = dict(row=r["row"],
                            instructions=r.get("price", len(r["ops"])),
                            ops=r.get("text", seq_text(r["ops"])),
                            runs=0, mismatches=0)
    runs, mis = [0] * 10, [0] * 10
    rng = random.Random(seed + 991)
    regs = dict(zr=0, one=1, two=2, three=3, six=6, twenty=20, price=1099,
                m64m81=TWO64_M81, p0=PROD["q"] & MASK64, ones=MASK64, gmem=0x300000)
    regs.update(SEQ_CONSTS)
    gmem = dict((0x300000 + 8 * i, GL[i + 1]) for i in range(5))
    composed = dict(runs=0, mismatches=0, fallbacks=0, dp=0)
    q = PROD["q"]
    done = 0
    for i in range(n_random):
        if (i & 2047) == 0 and deadline is not None and time.time() > deadline:
            break
        done = i + 1
        x = rng.randrange(1, q)
        y = rng.randrange(1, q)
        fj = fast_j(x)
        xl = _limbs(x)
        runs[0] += 1
        mis[0] += _chk_f1(xl, regs, gmem)
        F0 = xl[3]
        F1 = _ref_f1(xl)
        runs[1] += 1
        bad, wrapped = _chk_wrap(regs, F0, F1)
        mis[1] += bad
        F2 = (-(F0 + F1)) & MASK64
        fallback = wrapped or F1 < TOL or F1 > M64 - TOL
        if not fallback:
            runs[2] += 1
            bad, fallback = _chk_decide(regs, F0, F1, F2)
            mis[2] += bad
        composed["runs"] += 1
        composed["mismatches"] += int(fallback != (fj is None))
        if fallback:
            composed["fallbacks"] += 1
            continue
        composed["mismatches"] += int(regs["j"] != fj)
        ny = q - y
        K = y if y <= ny else ny
        yl = _limbs(y)
        runs[3] += 1
        mis[3] += _chk_sign(yl, regs)
        runs[4] += 1
        kl = _limbs(K)
        bad, dp = _chk_route(regs, kl, K, 1 if ny < y else 0, fj)
        mis[4] += bad
        if dp:
            composed["dp"] += 1
            continue
        runs[5] += 1
        mis[5] += _chk_traffic(regs, xl, yl, x, y, regs["xT"])
        fp = _fp(K)
        if i % 6 == 0:
            slot = list(fp)
        elif i % 6 == 1:
            slot = [fp[0], fp[1] ^ 1, fp[2]]
        elif i % 6 == 2:
            slot = _fp(K ^ (1 << 64))      # differs only in limb 1: the fingerprint agrees
        else:
            slot = [x & MASK64, (x >> 128) & MASK64, (y >> 64) & MASK64]
        if i % 5 == 0:
            cp = list(fp)
        elif i % 5 == 1:
            cp = _fp(K ^ (1 << 128))
        else:
            cp = [y & MASK64, (y >> 96) & MASK64, (y >> 192) & MASK64]
        length = (((x >> 128) % 22) << 64) | ((x * 3) & MASK64)
        runs[6] += 1
        mis[6] += _chk_hist(regs, K, kl, slot, cp, 1 if i % 8 == 0 else 0, length)[0]
        runs[8] += 1
        mis[8] += _chk_meter(regs, length, (y >> 192) % 5, (y >> 128) % 4)

    # ---- boundary inputs aimed at each row's own branch
    E = [0, 1, 2, 39, 40, 41, 42, M64 - 42, M64 - 41, M64 - 40, M64 - 1, M64, 1 << 63]
    for F0 in E:
        for F1 in E:
            runs[1] += 1
            mis[1] += _chk_wrap(regs, F0, F1)[0]
    for F0 in E:
        for F1 in E:
            for F2 in E:
                runs[2] += 1
                mis[2] += _chk_decide(regs, F0, F1, F2)[0]
    Y = [0, 1, 2, (1 << 32) - 1, 1 << 32, (1 << 32) + 1, M64, 1 << 64, (1 << 64) + 1,
         (1 << 128) - 1, 1 << 128, 1 << 255, q - 2, q - 1, (q - 1) // 2, (q + 1) // 2,
         q - 1 - (1 << 64), 3 * M64]
    for y in Y:
        runs[3] += 1
        mis[3] += _chk_sign(_limbs(y), regs)
    KS = [0, 1, INDEX_MASK, INDEX_MASK + 1, M64, 1 << 64, (1 << 64) + 1, (1 << 64) | INDEX_MASK,
          1 << 128, M64 << 64, (M64 << 64) | M64, (1 << 192) - 1,
          (M64 << 192) | (1 << 128) | (5 << 64) | 7]
    for K in KS:
        for s in (0, 1):
            for j in (0, 1, 2):
                runs[4] += 1
                mis[4] += _chk_route(regs, _limbs(K), K, s, j)[0]
    for (x, y) in [(0, 0), (1, 1), (q - 1, q - 1), (M64 << 192, M64), ((1 << 128) - 1, (1 << 192) - 1),
                   (1 << 255, 3 * M64)]:
        runs[5] += 1
        mis[5] += _chk_traffic(regs, _limbs(x), _limbs(y), x, y,
                                                    0x100000 + (0x1ABCD << 6))
    KH = (M64 << 192) | (0x1234 << 128) | (M64 << 64) | 0x5678
    fph = _fp(KH)
    for flag in (0, 1):
        for slot in (list(fph), [fph[0] ^ 1, fph[1], fph[2]], [fph[0], fph[1] ^ 1, fph[2]],
                     [fph[0], fph[1], fph[2] ^ 1], _fp(KH ^ (1 << 64)), [0, 0, 0], [M64] * 3):
            for cp in (list(fph), [fph[0] ^ 1, fph[1], fph[2]], _fp(KH ^ (1 << 128)),
                       [fph[0], fph[1], fph[2] ^ (1 << 63)]):
                for length in (0, 1, 63, 64, 65, 128):
                    runs[6] += 1
                    mis[6] += _chk_hist(regs, KH, _limbs(KH), slot, cp, flag, length)[0]
    for len1 in (0, 19, 20, 21):
        for len0 in (0, 1, M64):
            for mtr in (0, 1, 2):
                for cnt in (0, 1, 2):
                    runs[8] += 1
                    mis[8] += _chk_meter(regs, (len1 << 64) | len0, mtr, cnt)
    for b in range(120):
        runs[7] += 1
        mis[7] += _chk_batch(rng, 64, None if b % 2 else rng.randrange(64))
    runs[9] = BATCH_MAX - BATCH_MIN + 1
    mis[9] = _chk_bound()
    for i, key in enumerate(SEQ_ORDER):
        st[key]["runs"], st[key]["mismatches"] = runs[i], mis[i]
    total = sum(st[k]["instructions"] for k in st)
    return dict(rows=[st[k] for k in ("F1", "wrap", "decide", "sign", "route", "traffic",
                                      "hist", "zero", "meter", "batch")],
                itemised_total=total,
                random_inputs=done,
                random_requested=n_random,
                deadline_hit=done < n_random,
                composed=composed,
                semantics_check=dict(mnemonics=len(SEQ_SEMANTICS), undocumented=SEQ_UNDOCUMENTED),
                total_mismatches=sum(st[k]["mismatches"] for k in st) + composed["mismatches"])


# ---------------------------------------------------------------- Part A: production checks
def part_a(nchecks):
    out = {}
    E = Curve(PROD["q"], PROD["r"], PROD["h"])
    q, r = E.q, E.r
    E.beta, E.lam = PROD["beta"], PROD["lam"]
    G = (PROD["gx"], PROD["gy"])
    E.G = G
    out["P0_on_curve"] = E.on_curve(G)
    out["r_P0_is_infinity"] = E.mul(r - 1, G) == E.neg(G)
    out["q_mod_3"] = q % 3
    out["beta_cube_root_of_unity"] = (E.beta * E.beta + E.beta + 1) % q == 0
    out["lambda_cube_root_of_unity"] = (E.lam * E.lam + E.lam + 1) % r == 0
    out["phi_P0_equals_lambda_P0"] = E.mul(E.lam, G) == E.phi(G)
    out["cofactor_1_r_prime_order_group"] = PROD["h"] == 1
    out["log2_r"] = math.log2(r)
    out["log2_N_classes"] = math.log2((r - 1) // 6)
    rm1, fac = r - 1, {}
    for pr in [2, 3, 149, 631, 107361793816595537, 174723607534414371449, 341948486974166000522343609283189]:
        while rm1 % pr == 0:
            fac[pr] = fac.get(pr, 0) + 1
            rm1 //= pr
    order = r - 1
    for pr, ex in fac.items():
        for _ in range(ex):
            if pow(2, order // pr, r) == 1:
                order //= pr
            else:
                break
    out["r_minus_1_fully_factored"] = rm1 == 1 and all(pow(2, order // pr, r) != 1 for pr in fac if order % pr == 0)
    out["log2_order_of_2_mod_r"] = round(math.log2(order), 4)
    xs = rng.randrange(1, r)
    T = E.mul(xs, G)
    W = Walker(E, 4, 0, T, rng)  # tiny table, for identity checks only
    ok_orbit = ok_key = ok_canon = ok_twist = ok_double = True
    for _ in range(nchecks * 4):
        k = rng.randrange(1, r)
        P = E.mul(k, G)
        members = [E.mul(u * k % r, G) for u in W.units]
        for idx, u in enumerate(W.units):
            s, j = divmod(idx, 3)
            Pj = E.phi(P, j)
            if s:
                Pj = E.neg(Pj)
            if members[idx] != Pj:
                ok_orbit = False
        if len(set(members)) != 6:
            ok_orbit = False
        cans = [W.canon(Mm) for Mm in members]
        if len({c[0] for c in cans}) != 1 or cans[0][0] != min(P[1], q - P[1]):
            ok_key = False
        if len({c[2] for c in cans}) != 1:
            ok_canon = False
        for Mm, (K, t, c) in zip(members, cans):
            if E.mul(W.units[t], Mm) != c:
                ok_canon = False
            h = K & (W.R - 1)
            TP = W.table[h][t][0]
            Mh = W.table[h][0][0]
            n1 = W.canon(E.add(Mm, TP))
            n2 = W.canon(E.add(c, Mh))
            if n1[2] != n2[2]:
                ok_twist = False
            if W.canon(E.add(Mm, Mm))[2] != W.canon(E.add(c, c))[2]:
                ok_double = False
    out["orbit_size_6_and_unit_action"] = ok_orbit
    out["key_abs_y_is_class_invariant"] = ok_key
    out["canonical_member_unique"] = ok_canon
    out["twisted_step_is_class_map"] = ok_twist
    out["doubling_is_class_map"] = ok_double
    ok_tab = all(W.table[h][t][0] == E.mul(W.uinv[t], W.table[h][0][0]) for h in range(W.R) for t in range(6))
    out["twisted_table_entries_equal_u_inverse_M"] = ok_tab

    Wp = Walker(E, 16, 64, T, rng)  # small table (R does not affect these identities)
    a, b = rng.randrange(r), rng.randrange(r)
    P = E.add(E.mul(a, G), E.mul(b, T))
    steps = 6 * nchecks
    trailA = []
    for i in range(steps):
        K, t, c = Wp.canon(P)
        TP, tc, td = Wp.step_target(P, K, t)
        P = E.add(P, TP)
        a, b = (a + tc) % r, (b + td) % r
        trailA.append((P, a, b))
    out["production_steps_coefficients_ok"] = E.add(E.mul(a, G), E.mul(b, T)) == P
    i0 = steps // 3
    PA, aA, bA = trailA[i0]
    v = rng.randrange(6)
    u = Wp.units[v]
    bB = rng.randrange(r)
    aB = (u * (aA + bA * xs) - bB * xs) % r
    PB = E.mul(u, PA)
    assert E.add(E.mul(aB, G), E.mul(bB, T)) == PB
    for i in range(i0 + 1, steps):
        K, t, c = Wp.canon(PB)
        TP, tc, td = Wp.step_target(PB, K, t)
        PB = E.add(PB, TP)
        aB, bB = (aB + tc) % r, (bB + td) % r
    PA, aA, bA = trailA[-1]
    KA, tA, cA = Wp.canon(PA)
    KB, tB, cB = Wp.canon(PB)
    out["merged_trails_same_canonical"] = cA == cB
    uA, uB = Wp.units[tA], Wp.units[tB]
    a1, b1 = uA * aA % r, uA * bA % r  # canonical coefficients of trail A
    a2, b2 = uB * aB % r, uB * bB % r
    xsol = (a2 - a1) * pow(b1 - b2, -1, r) % r
    out["solve_x_equals_(aB-aA)/(bA-bB)"] = xsol == xs and E.mul(xsol, G) == T
    nr = 15000 * nchecks if canon_random is None else canon_random
    ns = 5000 * nchecks if canon_stress is None else canon_stress
    out["fast_canonical_index"] = check_fast_j(nr, ns)
    out["part_a_seconds"] = round(time.time() - T0, 2)
    return out


# ---------------------------------------------------------------- Part B: toy run
def toy_trial(E, rng, R, dp_bits):
    r = E.r
    xs = rng.randrange(1, r)
    T = E.fmul(xs)

    def point(a, b):
        return E.fmul(a + b * xs)

    W = Walker(E, R, dp_bits, T, rng, point)
    N = (r - 1) // CLS
    first_seen = {}  # class key -> trail id (instrumentation only)
    state = {"visits": 0, "tc": None, "trail": None}

    def on_class(K):
        state["visits"] += 1
        tid = state["trail"]
        prev = first_seen.get(K)
        if prev is None:
            first_seen[K] = tid
        elif prev != tid and state["tc"] is None:
            state["tc"] = state["visits"]

    dp_table = {}
    trail_id = 0
    solved = None
    Np = N * R / (R - 1)
    B = 1.05 * math.sqrt(2 * Np * math.log(1 / 0.605))
    D_cap = int(2 * (B + W.Lmax) / (1 << dp_bits)) + 1
    NODES = cfg.get("nodes", 16)
    node_shift = W.rbits + dp_bits
    node_cap = -(-2 * D_cap // NODES)
    node_load = [0] * NODES
    matches = 0
    while solved is None:
        if pipeline and state["visits"] >= B:
            solved = False  # phase-1 budget exhausted and the open trail is finished
            break
        a, b = rng.randrange(r), rng.randrange(r)
        P = point(a, b)
        if P is None:
            continue
        trail_id += 1
        state["trail"] = trail_id
        status, K, length, c, _ = W.run_trail(P, None, on_class)
        if status != "dp":
            continue
        if K in dp_table:
            matches += 1
            if pipeline and matches > 8:
                solved = False
                break
            (a2, b2, l2) = dp_table[K]
            s1 = W.stats.copy()
            _, K1, _, c1, (A1, B1) = W.run_trail(point(a, b), (a, b))
            _, K2, _, c2, (A2, B2) = W.run_trail(point(a2, b2), (a2, b2))
            W.stats.update(s1)  # re-walk cost is accounted separately
            if c1 == c2 and (B1 - B2) % r:
                x = (A2 - A1) * pow(B1 - B2, -1, r) % r
                if E.fmul(x) == T:
                    solved = (x == xs)
                    break
        else:
            if pipeline and len(dp_table) >= D_cap:
                solved = False
                break
            if pipeline:
                nd = (K >> node_shift) % NODES
                if node_load[nd] >= node_cap:
                    solved = False
                    W.stats["node_overflow"] = W.stats.get("node_overflow", 0) + 1
                    break
                node_load[nd] += 1
            dp_table[K] = (a, b, length)
        if time.time() - T0 > budget:
            return None
    nl = dict(max_node_load_over_mean=(max(node_load) * NODES / max(1, sum(node_load))) if pipeline and sum(node_load) else None,
              max_node_load_over_cap=(max(node_load) / node_cap) if pipeline else None)
    return dict(tc=state["tc"], solve_visits=state["visits"], ok=solved, trails=trail_id, **nl, **W.stats)


def long_trail_test(E, rng, R, L, ntrails, escape_add):
    """Many escapes per trail: walks with no distinguished points run L operations
    each. Separates repeats inside fruitless cycles of length <= 64 (same escape
    epoch, <= 64 arrivals since the previous visit) from all other returns of a trail to its own past ('non-
    fruitless returns': across an escape, or after a longer gap). H1 says the
    latter are no more frequent than ordinary self-collisions (a trail has one with
    probability about 1 - exp(-L^2/(2N))); a structural trap would make them frequent. escape_add = True runs the rejected additive escape as a control."""
    r = E.r
    xs = rng.randrange(1, r)
    T = E.fmul(xs)
    W = Walker(E, R, 0, T, rng, lambda c, d: E.fmul(c + d * xs))
    W.is_dp = lambda K: False
    W.Lmax = L
    if escape_add:
        W.escape_add = True
        W.esc_row = [W.table[0][t] for t in range(6)]
        W.table = [W.table[R - 1]] + W.table[1:]  # M_esc not used as a partition multiplier
    out = dict(trails=0, operations=0, escapes=0, nonfruitless_returns=0, fruitless_repeats=0,
               trails_with_nonfruitless_returns=0)
    N = (r - 1) // 6
    for _ in range(ntrails):
        if time.time() - T0 > budget:
            break
        P = E.fmul(rng.randrange(1, r))
        seen = {}
        cnt = [0, 0]  # [non-fruitless returns, fruitless-cycle repeats]
        arr = [0]

        def on_class(K):
            ep = W.stats["doublings"]
            arr[0] += 1
            if K in seen:
                pep, parr = seen[K]
                if pep < ep or arr[0] - parr > 64:
                    cnt[0] += 1
                else:
                    cnt[1] += 1
            seen[K] = (ep, arr[0])

        s0 = dict(W.stats)
        W.run_trail(P, None, on_class)
        out["trails"] += 1
        out["operations"] += W.stats["steps"] + W.stats["doublings"] - s0["steps"] - s0["doublings"]
        out["escapes"] += W.stats["doublings"] - s0["doublings"]
        out["nonfruitless_returns"] += cnt[0]
        out["fruitless_repeats"] += cnt[1]
        out["trails_with_nonfruitless_returns"] += cnt[0] > 0
    out["expected_self_collisions_per_trail_L2_over_2N"] = round(L * L / (2 * N), 4)
    return out


def main():
    out = {"config": cfg}
    if seq_n:
        # section 5's rows, executed; leave the toy phase and the output some of the budget
        out["sequence_check"] = check_sequences(int(seq_n), deadline=T0 + budget - 4)
    if prod_checks:
        out["part_a"] = part_a(prod_checks)
    if long_trails:
        out["long_trail_test"] = dict(curve=curve, R=R, long_len=long_len, escape="additive (control)" if long_escape_add else "doubling",
                                      **long_trail_test(make_toy(curve), rng, R, long_len, long_trails, long_escape_add))
    if not run_toy:
        out["seconds"] = round(time.time() - T0, 2)
        print(json.dumps(out, indent=1))
        return
    E = make_toy(curve)
    r = E.r
    N = (r - 1) // CLS
    Np = N * R / (R - 1)
    S0 = math.sqrt(2 * Np * math.log(1 / 0.605))  # pure model: probability 0.395
    S1 = 1.05 * S0  # claimed phase-1 budget of proof.md section 4 (H1: >= 0.395)
    out["toy"] = dict(q=E.q, r=r, h=E.h, beta=E.beta, lam=E.lam, N=N, R=R, dp_bits=dp_bits,
                      mean_trail_over_sqrtN=round((1 << dp_bits) / math.sqrt(N), 4),
                      S0=round(S0, 1), S1=round(S1, 1),
                      model_mean_tc=round(math.sqrt(math.pi * Np / 2), 1))
    res = []
    for _ in range(trials):
        tr = toy_trial(E, rng, R, dp_bits)
        if tr is None:
            break
        res.append(tr)
    n = len(res)
    if n:
        tcs = [x["tc"] for x in res if x["tc"] is not None]
        steps = sum(x["steps"] for x in res)
        cens = pipeline  # in pipeline mode T_c is censored at the budget
        summ = dict(
            trials_completed=n,
            mode="pipeline (scaled budget, phase 2, report cap, 8-match cap)" if pipeline else "run until solved",
            pipeline_success_rate=round(sum(1 for x in res if x["ok"]) / n, 4) if pipeline else None,
            all_solved_correctly=None if pipeline else all(x["ok"] for x in res),
            mean_tc_over_sqrtN=None if cens else round(sum(tcs) / len(tcs) / math.sqrt(N), 4),
            stderr_mean_tc_over_sqrtN=None if cens else round(math.sqrt(sum((t / math.sqrt(N) - sum(tcs) / len(tcs) / math.sqrt(N)) ** 2 for t in tcs) / max(1, len(tcs) - 1) / len(tcs)), 4),
            model_mean_over_sqrtN=round(math.sqrt(math.pi * Np / 2) / math.sqrt(N), 4),
            frac_tc_le_S0=round(sum(t <= S0 for t in tcs) / n, 4),
            frac_tc_le_S1=round(sum(t <= S1 for t in tcs) / n, 4),
            binomial_stderr_at_0p39=round(math.sqrt(0.39 * 0.61 / n), 4),
            model_P_S0=0.395,
            H1_floor_P_S1=0.395,
            model_P_S1=round(1 - math.exp(-S1 * S1 / (2 * Np)), 4),
            mean_solve_over_tc=None if cens else round(sum(x["solve_visits"] / x["tc"] for x in res if x["tc"]) / len(tcs), 3),
            visits_per_step=round(sum(x["solve_visits"] for x in res) / max(1, sum(x["steps"] + x["doublings"] for x in res)), 4),
            total_steps=steps,
            doublings_per_step=sum(x["doublings"] for x in res) / max(1, steps),
            esc2_events=sum(x["esc2"] for x in res),
            esc2_expected_events=round(steps / (6 * R), 1),
            esc2_per_step=sum(x["esc2"] for x in res) / max(1, steps),
            esc2_model_1_over_6R=1 / (6 * R),
            checkpoint_escapes_per_step=sum(x["esccp"] for x in res) / max(1, steps),
            abandoned_trails=sum(x["abandoned"] for x in res),
            exceptional_additions=sum(x["exceptional"] for x in res),
            node_overflows=sum(x.get("node_overflow", 0) for x in res),
            max_node_load_over_cap=max((x["max_node_load_over_cap"] or 0) for x in res) if pipeline else None,
            mean_max_node_load_over_mean=(sum(x["max_node_load_over_mean"] or 0 for x in res) / n) if pipeline else None,
        )
        out["toy_summary"] = summ
        if dump:
            out["tc_list"] = tcs
    out["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(out, indent=1))


main()
