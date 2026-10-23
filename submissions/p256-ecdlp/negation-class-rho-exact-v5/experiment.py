"""Negation-class (two-class) parallel Pollard rho on NIST P-256, priced by an
executed P-256 field instruction listing.

The walk, its checks and its toy harness are those of submission
5fa78d62-1c6b-46b8-b9cf-02481a49e645 (negation-class-rho-v1), which adapts the
class-walk construction of submission 1c600df4 (mpjunior92, bls12-381-g1-ecdlp-plain,
glv-six-class-rho-v11) to P-256: classes {P, -P}, class key x, doubling escapes,
2-cycle and checkpoint detection, re-walk solving. New here is Part C.

Part A checks on P-256 itself: curve, order, order of 2 mod r, orbit size 2, the
x class key, the canonical member, the twisted (negated) table step as a class
map, the doubling class map, coefficient tracking on production steps and a
forced-merge solve. Part B runs the complete per-trail algorithm on prime-order
subgroups of toy curves y^2 = x^3 - 3x + b, measuring the class-visit count T_c
until the first collision between different trails against the random-mapping
model on N = (r'-1)/2 classes (classes = 1 runs the no-automorphism control).
Toy evidence only: it does not validate production cost bounds.

Part C executes the AArch64 field listing of field_listing.md (Montgomery
multiplication, modular addition and subtraction on 4 x 64-bit limbs, with the
AArch64 carry-flag rules) on a 64-bit register and memory model. It checks the
results on production P-256 elements, counts every executed instruction, runs
the Fermat inversion as a straight-line chain of those kernels, and runs batched
affine additions and doublings of production P-256 points through the kernels,
compared with reference affine arithmetic. It then recomputes the step prices and
the work ledger of proof.md section 5 from those counts.

The report store of proof.md section 3 is exact by construction: one record per
report, appended to one of 2^16 logs at the position its arrival index gives,
with no slot index, no probe sequence, no replacement and no occupancy argument;
duplicates are found by sorting each log on the full 256-bit key and merging the
sorted logs, so equal full keys are the only thing that can be reported as a
collision. ReportStore below is that structure at toy scale, and the pipeline
mode of Part B uses it (not a dictionary) to find the collision it solves on.
"""
import heapq, json, math, random, sys, time

T0 = time.time()
cfg = json.load(sys.stdin) if __name__ == "__main__" else {}
curve = cfg.get("curve", "q26")
trials = cfg.get("trials", 100)
R = cfg.get("R", 256)
dp_bits = cfg.get("dp_bits", 4)
seed = cfg.get("seed", 1)
prod_checks = cfg.get("production_checks", 2)
listing_checks = cfg.get("listing_checks", 0)
budget = cfg.get("seconds", 14)
dump = cfg.get("dump", False)
# pipeline = True: also apply the production resource rules scaled to the toy (phase-1
# budget, phase 2, report-table capacity, 8-match cap) and report the success rate.
pipeline = cfg.get("pipeline", False)
# long-trail trap test (0 = off): trails of long_len operations without distinguished points
long_trails = cfg.get("long_trails", 0)
long_len = cfg.get("long_len", 4096)
long_escape_add = cfg.get("long_escape_add", False)
# toy = False skips Part B (e.g. to spend the whole run on Part A checks)
run_toy = cfg.get("toy", True)
# control option: 6 = the attack (default); 1 = plain R-adding walk on points (no
# automorphism); 2 = negation-only classes. Used only to compare toy collision times.
CLS = cfg.get("classes", 2)
rng = random.Random(seed)

# ---------------------------------------------------------------- curves
# Production profile constants (tracks/p256-ecdlp/target-profile.json); cofactor 1.
PROD = dict(
    q=int("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff", 16),
    r=int("ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551", 16),
    h=1,
    gx=int("6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296", 16),
    gy=int("4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5", 16),
    a=-3,
    b=int("5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b", 16),
)
# Toy curves y^2 = x^3 - 3x + b over F_q, #E = h * r' with r' prime (PARI/GP ellcard).
TOYS = {
    "q26": (33632567, 13, 33633526, 2),
    "q30": (536949043, 25, 536950955, 5),
}


class Curve:
    """Affine arithmetic on y^2 = x^3 + a x + b over F_q (tuples, None = infinity)."""

    def __init__(self, q, r, h, a, b):
        self.q, self.r, self.h, self.a, self.b = q, r, h, a % q, b % q

    def on_curve(self, P):
        return P is None or (P[1] * P[1] - P[0] ** 3 - self.a * P[0] - self.b) % self.q == 0

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
            l = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, q) % q
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


def make_toy(name):
    q, B, n, h = TOYS[name]
    r = n // h
    E = Curve(q, r, h, -3, B)
    while True:
        x = rng.randrange(q)
        y = sqrt_mod(x ** 3 - 3 * x + B, q)
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
        # units u_t = (-1)^t, t = 0, 1 (negation classes; no endomorphism on P-256)
        self.units = [1, r - 1]
        self.uinv = [pow(u, -1, r) for u in self.units]
        self.T = T
        self.table = []  # table[h][t] = (point u_t^{-1} M_h, coeff c, coeff d)
        while len(self.table) < R:
            c, d = rng.randrange(r), rng.randrange(r)
            # M_h = [c]G + [d]T; the toy harness evaluates the same point as [c + d x]G
            # with a fixed-base table (start_point), production part A uses plain scalar mults
            M = start_point(c, d) if start_point else E.add(E.mul(c, E.G), E.mul(d, T))
            if M is None:  # M_h = O: redraw (probability 1/r)
                continue
            row = []
            for t in range(len(self.units)):
                Mt = E.neg(M) if t else M          # u_t^{-1} M = (-1)^t M
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
        """Return (K, t, canonical point): canonical = units[t] * P = (x, |y|), K = x."""
        q = self.E.q
        x, y = P
        if CLS == 1:  # control: point itself, key injective on points
            return x + q * y, 0, P
        if q - y < y:
            return x, 1, (x, q - y)
        return x, 0, P

    def step_target(self, P, K, t):
        h = K & (self.R - 1)
        # P is any member; canonical c = u_t P; next class = class(c + M_h) = class(P + u_t^{-1} M_h)
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
                # 2-cycle {current class, previous class}: double the smaller key
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
                # cycle of length <= 64 (checkpoint repeat): walk once more round it,
                # remember the member of the minimum-key class, then double it.
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


# ---------------------------------------------------------------- Part A: production checks
def part_a(nchecks):
    out = {}
    E = Curve(PROD["q"], PROD["r"], PROD["h"], PROD["a"], PROD["b"])
    q, r = E.q, E.r
    G = (PROD["gx"], PROD["gy"])
    E.G = G
    out["P0_on_curve"] = E.on_curve(G)
    # E.mul reduces scalars mod r, so check [r]P0 = O via [r-1]P0 = -P0
    out["r_P0_is_infinity"] = E.mul(r - 1, G) == E.neg(G)
    out["cofactor_1_r_prime_order_group"] = PROD["h"] == 1
    out["log2_r"] = math.log2(r)
    out["log2_N_classes"] = math.log2((r - 1) // 2)
    # multiplicative order of 2 mod r (doubling escapes cannot cancel: +-2^m != 1 for
    # every m < ord/2), from the factorisation of r - 1
    rm1, fac = r - 1, {}
    for pr in [2, 3, 71, 131, 373, 3407, 17449, 38189, 187019741, 622491383, 1002328039319,
               2624747550333869278416773953]:
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
        if members[1] != E.neg(P) or members[0] == members[1]:
            ok_orbit = False
        cans = [W.canon(Mm) for Mm in members]
        if len({c[0] for c in cans}) != 1 or cans[0][0] != P[0]:
            ok_key = False
        if len({c[2] for c in cans}) != 1 or cans[0][2][1] != min(P[1], q - P[1]):
            ok_canon = False
        for Mm, (K, t, c) in zip(members, cans):
            if E.mul(W.units[t], Mm) != c:
                ok_canon = False
            h = K & (W.R - 1)
            TP = W.table[h][t][0]
            Mh = W.table[h][0][0]
            if W.canon(E.add(Mm, TP))[2] != W.canon(E.add(c, Mh))[2]:
                ok_twist = False
            if W.canon(E.add(Mm, Mm))[2] != W.canon(E.add(c, c))[2]:
                ok_double = False
    out["orbit_size_2_and_unit_action"] = ok_orbit
    out["key_x_is_class_invariant"] = ok_key
    out["canonical_member_unique"] = ok_canon
    out["twisted_step_is_class_map"] = ok_twist
    out["doubling_is_class_map"] = ok_double
    ok_tab = all(W.table[h][t][0] == E.mul(W.uinv[t], W.table[h][0][0]) for h in range(W.R) for t in range(2))
    out["twisted_table_entries_equal_u_inverse_M"] = ok_tab

    # production walk with coefficient tracking, and a forced merge + solve
    Wp = Walker(E, 16, 64, T, rng)
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
    v = rng.randrange(2)
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
    a1, b1 = uA * aA % r, uA * bA % r
    a2, b2 = uB * aB % r, uB * bB % r
    xsol = (a2 - a1) * pow(b1 - b2, -1, r) % r
    out["solve_x_equals_(aB-aA)/(bA-bB)"] = xsol == xs and E.mul(xsol, G) == T
    out["part_a_seconds"] = round(time.time() - T0, 2)
    return out


# ---------------------------------------------------------------- Part C: counted P-256 field listing
M64 = (1 << 64) - 1
PL = [M64, 0xFFFFFFFF, 0, 0xFFFFFFFF00000001]


def listing_mul():
    """Montgomery multiplication [x2] <- [x0]*[x1]*2^-256 mod p; one tuple per instruction."""
    s = []
    e = lambda *a: s.append(a)
    e("ldp", "A0", "A1", "x0", 0)
    e("ldp", "A2", "A3", "x0", 16)
    e("ldp", "B0", "B1", "x1", 0)
    e("ldp", "B2", "B3", "x1", 16)
    e("mov", "t2", "xzr")
    for k in range(7):
        terms = [(i, k - i) for i in range(4) if 0 <= k - i < 4]
        if k == 0:
            e("mul", "t0", "A0", "B0")
            e("umulh", "t1", "A0", "B0")
            continue
        for n, (i, j) in enumerate(terms):
            e("mul", "lo", "A%d" % i, "B%d" % j)
            e("umulh", "hi", "A%d" % i, "B%d" % j)
            e("adds", "t%d" % k, "t%d" % k, "lo")
            e("adcs", "t%d" % (k + 1), "t%d" % (k + 1), "hi")
            e("adc", "t%d" % (k + 2), "xzr" if n == 0 else "t%d" % (k + 2), "xzr")
    for r in range(4):
        m = "t%d" % r
        e("lsr", "mhi", m, 32)
        e("lsl", "mlo", m, 32)
        e("adds", "t%d" % (r + 1), "t%d" % (r + 1), "mlo")
        e("adcs", "t%d" % (r + 2), "t%d" % (r + 2), "mhi")
        e("adcs", "t%d" % (r + 3), "t%d" % (r + 3), m)
        e("adcs", "t%d" % (r + 4), "t%d" % (r + 4), m)
        hi = list(range(r + 5, 8)) + ([8] if r == 3 else [])
        for n, k in enumerate(hi):
            e("adc" if n == len(hi) - 1 else "adcs", "t%d" % k, "t%d" % k, "xzr")
        e("subs", "t%d" % (r + 3), "t%d" % (r + 3), "mlo")
        e("sbcs", "t%d" % (r + 4), "t%d" % (r + 4), "mhi")
        for n, k in enumerate(hi):
            e("sbc" if n == len(hi) - 1 else "sbcs", "t%d" % k, "t%d" % k, "xzr")
    e("subs", "r0", "t4", "P0")
    e("sbcs", "r1", "t5", "P1")
    e("sbcs", "r2", "t6", "P2")
    e("sbcs", "r3", "t7", "P3")
    e("sbcs", "t8", "t8", "xzr")
    for k in range(4):
        e("csel", "r%d" % k, "r%d" % k, "t%d" % (4 + k), "cs")
    e("stp", "r0", "r1", "x2", 0)
    e("stp", "r2", "r3", "x2", 16)
    return s


def listing_add():
    """[x2] <- [x0] + [x1] mod p."""
    s = []
    e = lambda *a: s.append(a)
    e("ldp", "A0", "A1", "x0", 0)
    e("ldp", "A2", "A3", "x0", 16)
    e("ldp", "B0", "B1", "x1", 0)
    e("ldp", "B2", "B3", "x1", 16)
    e("adds", "r0", "A0", "B0")
    e("adcs", "r1", "A1", "B1")
    e("adcs", "r2", "A2", "B2")
    e("adcs", "r3", "A3", "B3")
    e("adc", "cy", "xzr", "xzr")
    e("subs", "d0", "r0", "P0")
    e("sbcs", "d1", "r1", "P1")
    e("sbcs", "d2", "r2", "P2")
    e("sbcs", "d3", "r3", "P3")
    e("sbcs", "cy", "cy", "xzr")
    for k in range(4):
        e("csel", "r%d" % k, "d%d" % k, "r%d" % k, "cs")
    e("stp", "r0", "r1", "x2", 0)
    e("stp", "r2", "r3", "x2", 16)
    return s


def listing_sub():
    """[x2] <- [x0] - [x1] mod p."""
    s = []
    e = lambda *a: s.append(a)
    e("ldp", "A0", "A1", "x0", 0)
    e("ldp", "A2", "A3", "x0", 16)
    e("ldp", "B0", "B1", "x1", 0)
    e("ldp", "B2", "B3", "x1", 16)
    e("subs", "r0", "A0", "B0")
    e("sbcs", "r1", "A1", "B1")
    e("sbcs", "r2", "A2", "B2")
    e("sbcs", "r3", "A3", "B3")
    e("sbc", "cy", "xzr", "xzr")
    for k in range(4):
        e("and", "d%d" % k, "P%d" % k, "cy")
    e("adds", "r0", "r0", "d0")
    e("adcs", "r1", "r1", "d1")
    e("adcs", "r2", "r2", "d2")
    e("adcs", "r3", "r3", "d3")
    e("stp", "r0", "r1", "x2", 0)
    e("stp", "r2", "r3", "x2", 16)
    return s


def execute(seq, regs, mem):
    """Run seq on 64-bit registers and a word memory; return the executed-instruction count.
    AArch64 flags: ADDS, ADCS, SUBS, SBCS write C; ADC and SBC read C and write no flag."""
    C = 0
    rd = lambda n: 0 if n == "xzr" else regs[n]

    def wr(n, v):
        if n != "xzr":
            regs[n] = v & M64

    for ins in seq:
        op = ins[0]
        if op == "ldp":
            a = regs[ins[3]] + ins[4]
            wr(ins[1], mem[a])
            wr(ins[2], mem[a + 8])
        elif op == "stp":
            a = regs[ins[3]] + ins[4]
            mem[a], mem[a + 8] = rd(ins[1]), rd(ins[2])
        elif op == "mul":
            wr(ins[1], rd(ins[2]) * rd(ins[3]))
        elif op == "umulh":
            wr(ins[1], (rd(ins[2]) * rd(ins[3])) >> 64)
        elif op in ("adds", "adcs", "adc"):
            v = rd(ins[2]) + rd(ins[3]) + (0 if op == "adds" else C)
            wr(ins[1], v)
            if op != "adc":
                C = int(v > M64)
        elif op in ("subs", "sbcs", "sbc"):
            v = rd(ins[2]) - rd(ins[3]) - (0 if op == "subs" else 1 - C)
            wr(ins[1], v)
            if op != "sbc":
                C = int(v >= 0)
        elif op == "lsr":
            wr(ins[1], rd(ins[2]) >> ins[3])
        elif op == "lsl":
            wr(ins[1], rd(ins[2]) << ins[3])
        elif op == "and":
            wr(ins[1], rd(ins[2]) & rd(ins[3]))
        elif op == "mov":
            wr(ins[1], rd(ins[2]))
        elif op == "csel":
            assert ins[4] == "cs"
            wr(ins[1], rd(ins[2]) if C == 1 else rd(ins[3]))
        else:
            raise AssertionError(op)
    return len(seq)


class Kernels:
    """Field elements live in memory (4 little-endian words); each call runs one listing."""

    def __init__(self):
        self.code = dict(mul=listing_mul(), add=listing_add(), sub=listing_sub())
        self.mem, self.executed = {}, 0
        self.calls = dict(mul=0, add=0, sub=0)

    def put(self, addr, v):
        for k in range(4):
            self.mem[addr + 8 * k] = (v >> (64 * k)) & M64

    def get(self, addr):
        return sum(self.mem[addr + 8 * k] << (64 * k) for k in range(4))

    def op(self, name, a, b, c):
        regs = dict(x0=a, x1=b, x2=c, P0=PL[0], P1=PL[1], P2=PL[2], P3=PL[3])
        self.executed += execute(self.code[name], regs, self.mem)
        self.calls[name] += 1

    def inv(self, p, src, dst):
        """dst <- src^(p-2): straight-line left-to-right chain, 255 squarings and 127 multiplications."""
        bits = bin(p - 2)[3:]
        self.op("mul", src, src, dst)
        if bits[0] == "1":
            self.op("mul", dst, src, dst)
        for bit in bits[1:]:
            self.op("mul", dst, dst, dst)
            if bit == "1":
                self.op("mul", dst, src, dst)


def batch_step(Kn, p, pts, tabs, dbl):
    """One batched affine step on production points through the listing: element i computes
    pts[i] + tabs[i] (6M + 6A per element) or, if dbl[i], 2*pts[i] with a = -3 (7M + 8A per element).
    The batch shares one inversion (Montgomery's trick). Returns the results in affine form."""
    Rm = pow(2, 256, p)
    ONE, U = 0x40, 0x80
    Kn.put(ONE, Rm % p)
    n = len(pts)
    blk = lambda i, off: 0x100000 + 0x400 * i + off
    X, Y, TX, TY, D, Cp, INV, NUM, LAM, L2, T1, T2, T3, T4, X3, Y3 = [32 * k for k in range(16)]
    for i in range(n):
        Kn.put(blk(i, X), pts[i][0] * Rm % p)
        Kn.put(blk(i, Y), pts[i][1] * Rm % p)
        if not dbl[i]:
            Kn.put(blk(i, TX), tabs[i][0] * Rm % p)
            Kn.put(blk(i, TY), tabs[i][1] * Rm % p)
    caddr = []
    for i in range(n):
        if dbl[i]:
            Kn.op("add", blk(i, Y), blk(i, Y), blk(i, D))
        else:
            Kn.op("sub", blk(i, TX), blk(i, X), blk(i, D))
        if i == 0:
            caddr.append(blk(0, D))
        else:
            Kn.op("mul", caddr[-1], blk(i, D), blk(i, Cp))
            caddr.append(blk(i, Cp))
    Kn.inv(p, caddr[-1], U)
    for i in range(n - 1, 0, -1):
        Kn.op("mul", U, caddr[i - 1], blk(i, INV))
        Kn.op("mul", U, blk(i, D), U)
    inv_addr = [U] + [blk(i, INV) for i in range(1, n)]
    res = []
    for i in range(n):
        b = lambda off: blk(i, off)
        if dbl[i]:
            Kn.op("mul", b(X), b(X), b(T1))
            Kn.op("sub", b(T1), ONE, b(T2))
            Kn.op("add", b(T2), b(T2), b(T3))
            Kn.op("add", b(T3), b(T2), b(NUM))
            other = b(X)
        else:
            Kn.op("sub", b(TY), b(Y), b(NUM))
            other = b(TX)
        Kn.op("mul", b(NUM), inv_addr[i], b(LAM))
        Kn.op("mul", b(LAM), b(LAM), b(L2))
        Kn.op("sub", b(L2), b(X), b(T1))
        Kn.op("sub", b(T1), other, b(X3))
        Kn.op("sub", b(X), b(X3), b(T3))
        Kn.op("mul", b(LAM), b(T3), b(T4))
        Kn.op("sub", b(T4), b(Y), b(Y3))
        Ri = pow(Rm, -1, p)
        res.append((Kn.get(b(X3)) * Ri % p, Kn.get(b(Y3)) * Ri % p))
    return res


def store_path():
    """Proof.md section 3: the exact report store, itemised, in instructions.

    A report is one 104-byte record appended to the tail of one of 2^16 logs at
    the position its arrival index gives (index mod 2^16 is the log, the rest is
    the position). Nothing about the placement depends on the key, so there is no
    slot index, no home address, no probe sequence, no replacement rule and no
    load factor anywhere in the design. Detection is offline and exact: each log
    is radix-sorted on the full 256-bit key, the 2^16 sorted key streams are
    merged with a heap, and records with equal full keys are adjacent. Every item
    below is a count, and the failure table keeps only the count of reports
    against the declared capacity, never a storage geometry."""
    D_cap = 2 ** 64 - 2 ** 21                      # arrival-counter capacity, below its 64-bit range
    NODES = 2 ** 16
    REC = 104                                     # 32-byte K, a, b and 8 flags
    per_node = -(-D_cap // NODES)                  # ceil(D_cap / NODES) = 2^48 - 2^5, one record each
    items = dict(arrival_atomic=128,             # one 64-bit atomic fetch-and-add, contended
                 node_and_position=3,            # mask of the index, shift for the position
                 route_request=16 * (13 * 16 + 64),   # 16 edges, 13 words, one direction
                 append=64)                      # 13-word store (7 paired + 1), pointer, test
    rec_i = sum(items.values())
    sort_i = 8 * 24 + 16                          # 8 LSD passes of a 32-bit digit + adjacent scan
    merge_i = 4 * 16 * 16 + 512                   # key transfer over 16 edges + heap sift-down
    counters_i = 8 * 2 ** 32                      # per log: one clear per counter per pass
    readback_i = 16 * 16 * (13 * 16 + 64)         # <= 16 matched records, one routed read each
    units = dict(report=rec_i / 1000, sort=sort_i / 1000, merge=merge_i / 1000)
    return dict(records_cap=D_cap, log2_records_cap=round(math.log2(D_cap), 4),
                nodes=NODES, records_per_log_cap=per_node, log2_records_per_log_cap=round(math.log2(per_node), 4),
                bytes_per_record=REC, log2_store_bytes=round(math.log2(D_cap * REC), 4),
                log2_per_log_bytes=round(math.log2(per_node * REC), 4),
                report_items_instructions=items, report_instructions=rec_i,
                report_units=units["report"],
                sort_instructions_per_record=sort_i, merge_instructions_per_record=merge_i,
                units_one_off=dict(counters=counters_i * NODES / 1000, readback=readback_i / 1000),
                store_units_total=units["report"] * D_cap + units["sort"] * D_cap
                + units["merge"] * D_cap + counters_i * NODES / 1000 + readback_i / 1000)


def iroot(x, k):
    """floor(x ** (1/k)) for integers, exactly."""
    lo, hi = 0, 1 << (x.bit_length() // k + 1)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** k <= x:
            lo = mid
        else:
            hi = mid - 1
    return lo


def ledger(Mi, Ai, Ii, ctrl_budget, claim):
    """Proof.md sections 3 to 5 recomputed from the instruction counts (units of C = 1000).

    Every price is its counted instructions / 1000 (integer milli-units). The meter
    is the multiword one of section 3: a shared 64-bit counter of super-blocks of
    2^70 units, private 64-bit chunk counters inside a block, unspent remainders
    kept reserved (stranded, at most one block per processor). Phase 1 ends on a
    shared 64-bit counter of super-batches of 2^70 operations. The store's arrival
    counter halts at D_cap = 2^64 - 2^21. Every counter's largest value is checked
    against its word here as an exact integer inequality."""
    r = PROD["r"]
    lg = math.log2
    N = (r - 1) // 2
    Np = N * 2048 / 2047
    S0 = math.sqrt(2 * Np * math.log(1 / 0.605))
    S_ops = 1.05 * S0
    Lmax = 20 * 2 ** 64
    PROC = 2 ** 20                                 # processors, storage nodes included
    share = -(-Ii // 2 ** 15)
    counts = dict(step=6 * Mi + 6 * Ai + share + ctrl_budget, doubling=7 * Mi + 8 * Ai + share + ctrl_budget,
                  traversal=6 * Mi + 6 * Ai + share + ctrl_budget + 20, unbatched=Ii + 7 * Mi + 8 * Ai + 200)
    prices = {k: v / 1000.0 for k, v in counts.items()}
    ok_prices = all(counts[k] <= 1000 * prices[k] for k in counts)
    control_items = dict(sign_test_p_minus_y_select=13, index_dp_twist_table_address=7, loads_and_state_store=24,
                         fingerprint_history_checkpoint=26, zero_denominator=6,
                         length_two_words=4, cap_test_high_word=2, meter_private_update=5, loop=3,
                         batch_passes=9, phase1_private_down_counter=4, report_price_select=1)
    ctrl_itemised = sum(control_items.values())
    st = store_path()
    D_cap = st["records_cap"]
    step, unb = prices["step"], prices["unbatched"]

    # ---- the meter (section 3): exact integers, units of 2^70, milli-units inside a block
    BLK = 2 ** 70                                  # units per super-block
    CHK = 2 ** 50                                  # units per private chunk
    CHK_m = CHK * 1000                             # milli-units per chunk
    n_cap = iroot(2 ** (round(claim * 25) - 70 * 25), 25) - 1          # floor(2^(claim-70)) - 1 blocks
    PRE_BLK = 2 ** 30
    n_pre = 2 ** 30 - 1                            # preprocessing blocks: reserved <= 2^60 - 2^30 units
    RES_ITEMS = dict(contended_fetch_and_add=128, compare_with_cap=1, branch=1, private_resets=6)
    REF_ITEMS = dict(subtract_chunk=1, add_own_charge_8=1, store_m=1, load_c=1, add_c=1, store_c=1,
                     compare_c=1, branch=1)        # chunk carry: m <- m - 2^50*1000 + 8, c <- c + 1
    SUB_ITEMS = dict(reload_down_counter=1, store_down_counter=1, load_index=1, add_index=1, store_index=1,
                     compare_index=1, branch=1, add_own_charge_8=1)   # phase-1 sub-batch open
    RES_I = sum(RES_ITEMS.values())                # 136: contended fetch-and-add + compare, branch, private resets
    REF_I = sum(REF_ITEMS.values())                # 8, charged inside the block by the carry itself
    SUB_I = sum(SUB_ITEMS.values())                # 8, charged inside the block on the opening operation
    main_faa = n_cap + PROC                        # successful reservations + one over-cap per processor
    pre_faa = n_pre + PROC
    # W in milli-units, every run: the one list of terms (proof.md section 3, "Work bound, exact"):
    W_terms = dict(inside_blocks=n_cap * BLK * 1000,      # every charge inside a reserved main block
                   block_fetch_and_adds=main_faa * RES_I,  # every fetch-and-add on B, outside the blocks
                   halted_straddle=PROC * (2 ** 23 + REF_I),  # a halted processor's charge beyond its last block
                   preprocessing=2 ** 60 * 1000)          # its own meter, fetch-and-adds included
    W_m = sum(W_terms.values())
    W_bound_exact = W_m ** 25 <= 1000 ** 25 * 2 ** round(claim * 25)   # W <= 2^(3193/25) exactly
    pre_bound_exact = (n_pre * PRE_BLK * 1000 + pre_faa * RES_I + PROC * control_items["meter_private_update"]
                       <= 2 ** 60 * 1000)              # + a halted processor's last meter update (5)
    max_price_m = 2 ** 13 * 1000                   # largest single price (start, verification)
    J_cap = math.ceil(S_ops / BLK)                 # phase-1 super-batches of 2^70 operations
    J_cap = J_cap + 2 ** 8                         # float guard: S_ops is a float good to 2^75, so add 2^8 batches
    X = (2 ** 8 + 1) * BLK                         # phase 1 performs J_cap*2^70 >= S_ops and < S_ops + X < S_ops + 2^79
    words = dict(
        main_block_counter=(main_faa, 2 ** 64),
        main_chunk_index=(2 ** 20, 2 ** 64),
        main_milli_counter=(CHK_m + max_price_m, 2 ** 64),
        phase1_superbatch_counter=(J_cap + PROC, 2 ** 64),
        phase1_subbatch_index=(2 ** 20, 2 ** 64),
        phase1_op_down_counter=(2 ** 50, 2 ** 64),
        pre_block_counter=(pre_faa, 2 ** 64),
        pre_milli_counter=(PRE_BLK * 1000 + max_price_m, 2 ** 64),
        arrival_counter_final=(D_cap + PROC, 2 ** 64),
        arrival_index_largest_returned=(D_cap + PROC - 1, 2 ** 64),
        log_write_pointer=(st["records_per_log_cap"], 2 ** 64),
        radix_bucket_count=(st["records_per_log_cap"], 2 ** 64),
        trail_length_high_word=(20, 2 ** 64),
        trail_length_low_word=(2 ** 64 - 1, 2 ** 64),
        match_counter=(8, 2 ** 64),
    )
    every_counter_fits = all(v < w for v, w in words.values())

    # ---- phase 1, arrivals, H1
    doublings = (S_ops + X) * 2 ** -11.99
    S_eff = S_ops - doublings
    h1_at_S_ops = 1 - math.exp(-(S_ops / 1.05) ** 2 / (2 * Np))
    h1_at_arrivals = 1 - math.exp(-(S_eff / 1.05) ** 2 / (2 * Np))
    arrivals = S_ops + X + 2 ** 20 * (Lmax + 65)
    ph1 = step * (S_ops + X)
    hard = dict(phase2=2 ** 20 * (Lmax + 65) * unb,
                main_block_reservations=main_faa * RES_I / 1000,
                main_chunk_refills=main_faa * 2 ** 20 * REF_I / 1000,
                phase1_superbatch_reservations=(J_cap + PROC) * RES_I / 1000,
                phase1_subbatch_refills=(J_cap + PROC) * 2 ** 20 * SUB_I / 1000,
                rewalks=16 * (Lmax + 65) * unb,
                reports=D_cap * st["report_units"],
                merge=D_cap * st["merge_instructions_per_record"] / 1000,
                sort=D_cap * st["sort_instructions_per_record"] / 1000,
                radix_counters_cleared=st["units_one_off"]["counters"],
                readback=st["units_one_off"]["readback"],
                setup_table_solving_verification=2 ** 28)
    rand = dict(doubling_surcharge=(prices["doubling"] - step) * 2 ** -11.99 * (S_ops + X),
                traversal_surcharge=(prices["traversal"] - step) * 64 * 2 ** -19.99 * S_ops,
                starts=2 ** 76.7, others=2 ** 64)
    H = sum(hard.values())
    E_rand = sum(rand.values())
    W = 2 ** claim
    stranded = PROC * BLK                           # at most one partly used block per processor
    head_before = W - ph1 - H                       # the single-word meter's row (no stranding)
    head = n_cap * BLK - stranded - ph1 - H         # the multiword meter's row
    markov_before = E_rand / head_before
    markov = E_rand / head

    # ---- report count against D_cap (Chernoff, exact form)
    mu = arrivals * 2 ** -64
    delta = D_cap / mu - 1
    chern_exp = ((1 + delta) * math.log(1 + delta) - delta) * mu

    # ---- length-biased abandonment (section 4), H1's geometric model
    theta = 2.0 ** -64
    n_fresh = Lmax - Lmax // 4                      # fresh arrivals of a trail reaching L_max: H1 allows <= L_max/4 repeats
    x = n_fresh * theta                             # = 15 exactly
    # p(n) = n*theta*(1-theta)^(n-1) / (1 - (1-theta)^n) <= x * e^(-x) * e^theta / (1 - e^(-x))
    p_trail = x * math.exp(-x) * math.exp(theta) / (1 - math.exp(-x))
    p_trail_model_Lmax = 20 * math.exp(-20) * math.exp(theta) / (1 - math.exp(-20))
    grader_fraction = 20 * math.exp(-20) / (1 - math.exp(-20))
    rows = dict(nonfruitless_return_per_trail=Lmax ** 2 / (2 * N),
                nonfruitless_return_8_matches=8 * Lmax ** 2 / (2 * N),
                abandonment_per_match=2 * p_trail,
                abandonment_8_matches=16 * p_trail,
                report_count=2.0 ** -59,
                f1_eq_f2=1.0 / r, exceptional_additions=2.0 ** -114, eight_failing_matches=2.0 ** -60,
                spurious_fingerprint=2.0 ** -63, preprocessing_cap=2 ** 25 / ((n_pre - PROC) * PRE_BLK - 2 ** 33.2))
    summed = {k: v for k, v in rows.items() if k not in ("nonfruitless_return_per_trail", "abandonment_per_match")}
    other_failures = sum(summed.values())
    pre_head = (n_pre - PROC) * PRE_BLK - 2 ** 33.2
    up = lambda v, d: math.ceil(lg(v) * 10 ** d) / 10 ** d
    dn = lambda v, d: math.floor(lg(v) * 10 ** d) / 10 ** d
    cu = lambda v, d: math.ceil(v * 10 ** d) / 10 ** d
    fl = lambda v, d: math.floor(v * 10 ** d) / 10 ** d
    # success from the safely rounded terms, so the printed arithmetic reproduces the printed bound
    succ = fl(h1_at_arrivals, 6) - cu(markov, 6) - cu(other_failures, 9)
    raw = dict(log2_phase1=(lg(ph1), "up"), log2_hard_total=(lg(H), "up"), log2_random_total=(lg(E_rand), "up"),
               log2_headroom=(lg(head), "down"), log2_headroom_before=(lg(head_before), "down"),
               markov=(markov, "up"), markov_before=(markov_before, "up"), log2_deficit=(lg(doublings), "up"),
               h1_at_S_ops=(h1_at_S_ops, "down"), h1_at_arrivals=(h1_at_arrivals, "down"),
               h1_cost=(h1_at_S_ops - h1_at_arrivals, "up"), log2_mu=(lg(mu), "up"), mu_over_cap=(mu / D_cap, "up"),
               delta=(delta, "down"), log2_chernoff_exponent=(lg(chern_exp), "down"),
               p_trail=(p_trail, "up"), p_trail_at_Lmax=(p_trail_model_Lmax, "up"), per_match=(2 * p_trail, "up"),
               eight_matches=(16 * p_trail, "up"), grader_fraction=(grader_fraction, "up"),
               other_failures=(other_failures, "up"), success=(succ, "down"), margin=(succ - 0.39, "down"),
               log2_pre_markov=(lg(2 ** 25 / pre_head), "up"), log2_pre_headroom=(lg(pre_head), "down"),
               log2_nonfruitless=(lg(Lmax ** 2 / (2 * N)), "up"), log2_nonfruitless_8=(lg(8 * Lmax ** 2 / (2 * N)), "up"),
               log2_n_cap=(lg(n_cap), "down"), log2_J_cap=(lg(J_cap), "up"), log2_X=(lg(X), "up"),
               log2_main_faa=(lg(main_faa), "up"), log2_J_max=(lg(J_cap + PROC), "up"),
               **{"log2_hard_" + k: (lg(v), "up") for k, v in hard.items()},
               **{"log2_rand_" + k: (lg(v), "up") for k, v in rand.items()})
    return dict(
        counts_instructions=counts, prices_units=prices, every_price_covers_its_count=ok_prices,
        control_items=control_items, control_itemised=ctrl_itemised, control_budget=ctrl_budget,
        control_fits=ctrl_itemised <= ctrl_budget,
        store_path=st,
        meter=dict(block_units_log2=70, chunk_units_log2=50, n_cap=n_cap, log2_n_cap=round(lg(n_cap), 6),
                   n_pre=n_pre, main_faa_max=main_faa, W_terms_milli=W_terms, W_milli_max=W_m,
                   W_le_claim_exact=W_bound_exact, instruction_items=dict(
                       reservation=RES_ITEMS, chunk_refill=REF_ITEMS, subbatch_open=SUB_ITEMS,
                       meter_update=control_items["meter_private_update"],
                       down_counter=control_items["phase1_private_down_counter"]),
                   pre_le_2pow60_exact=pre_bound_exact, J_cap=J_cap, log2_J_cap=round(lg(J_cap), 6),
                   X_log2=lg(X), stranded_log2=lg(stranded)),
        counters={k: dict(max=v, word=w, fits=v < w, log2_max=round(lg(v), 4) if v > 0 else None)
                  for k, (v, w) in words.items()},
        every_counter_fits=every_counter_fits,
        log2_S0=round(lg(S0), 6), log2_S_ops=round(lg(S_ops), 6), log2_phase1_le=up(ph1, 5),
        log2_hard_le={k: up(v, 2) for k, v in hard.items()}, log2_hard_total_le=up(H, 3),
        log2_random_expected_le={k: up(v, 2) for k, v in rand.items()}, log2_random_total_le=up(E_rand, 3),
        log2_headroom_before_ge=dn(head_before, 3), markov_before_le=cu(markov_before, 6),
        log2_headroom_ge=dn(head, 3), markov_early_halt_le=cu(markov, 6),
        log2_arrival_deficit_le=up(doublings, 3),
        h1_floor_at_S_ops=fl(h1_at_S_ops, 6),
        h1_floor_at_arrivals=fl(h1_at_arrivals, 6),
        h1_cost_of_deficit_le=cu(h1_at_S_ops - h1_at_arrivals, 6),
        log2_expected_reports_le=up(mu, 3), mu_over_D_cap_le=cu(mu / D_cap, 4), chernoff_delta_ge=fl(delta, 5),
        chernoff_exponent_log2_ge=dn(chern_exp, 2),
        abandonment=dict(theta_log2=-64, Lmax="20*2^64", n_fresh="Lmax - Lmax/4", x=x,
                         grader_capped_visit_fraction=grader_fraction,
                         p_trail_at_Lmax_le=cu(p_trail_model_Lmax, 11),
                         p_trail_le=cu(p_trail, 11), per_match_le=cu(2 * p_trail, 11),
                         eight_matches_le=cu(16 * p_trail, 10)),
        pre_headroom_log2_ge=dn(pre_head, 4), pre_markov_log2_le=up(2 ** 25 / pre_head, 2),
        log2_failure_rows={k: (round(lg(v), 3) if v > 0 else None) for k, v in rows.items()},
        other_failures_sum_le=cu(other_failures, 9),
        success_lower_bound_ge=fl(succ, 6),
        margin_over_0p39_ge=fl(succ - 0.39, 6), raw=raw,
        claim_holds=(head > 0 and succ > 0.39 and W_bound_exact and pre_bound_exact and every_counter_fits
                     and ok_prices and ctrl_itemised <= ctrl_budget))


def part_c(nchecks):
    p = PROD["q"]
    assert p == PL[3] << 192 | PL[2] << 128 | PL[1] << 64 | PL[0]
    Kn = Kernels()
    out = {}
    hist = {}
    for ins in Kn.code["mul"]:
        hist[ins[0]] = hist.get(ins[0], 0) + 1
    out["listing_instructions"] = {k: len(v) for k, v in Kn.code.items()}
    out["mul_opcode_histogram"] = hist
    Ri = pow(2 ** 256, -1, p)
    edge = [0, 1, 2, p - 1, p - 2, 2 ** 255, 2 ** 192, 2 ** 64 - 1, 2 ** 224 - 1, (p - 1) // 2]
    pool = [0, 1, M64, M64 - 1, 1 << 32, (1 << 32) - 1, (1 << 32) + 1, 1 << 63, (M64 << 32) & M64]

    def operand(i):
        if i < len(edge) ** 2:
            return edge[i % len(edge)], edge[i // len(edge)]
        if i % 2:
            v = []
            for _ in range(2):
                x = 0
                for k in range(4):
                    x |= (rng.choice(pool) if rng.random() < 0.8 else rng.getrandbits(64)) << (64 * k)
                v.append(x % p)
            return v[0], v[1]
        return rng.randrange(p), rng.randrange(p)

    ok = dict(mul=True, add=True, sub=True)
    n_cases = len(edge) ** 2 + nchecks
    for i in range(n_cases):
        a, b = operand(i)
        Kn.put(0x1000, a)
        Kn.put(0x2000, b)
        for name, want in (("mul", a * b * Ri % p), ("add", (a + b) % p), ("sub", (a - b) % p)):
            Kn.op(name, 0x1000, 0x2000, 0x3000)
            if Kn.get(0x3000) != want:
                ok[name] = False
    out["kernel_cases_each"] = n_cases
    out["kernels_correct_on_reduced_inputs"] = ok
    code = Kn.code["mul"]
    mutant = [(("adc" if s[0] == "adcs" else "sbc",) + s[1:]) if s[0] in ("adcs", "sbcs") and j + 1 < len(code)
              and code[j + 1][0] in ("adc", "sbc") and code[j + 1][3] == "xzr" and s[1] != "r3" else s
              for j, s in enumerate(code)]
    wrong = 0
    for _ in range(200):
        a, b = rng.randrange(p), rng.randrange(p)
        Kn.put(0x1000, a)
        Kn.put(0x2000, b)
        regs = dict(x0=0x1000, x1=0x2000, x2=0x3000, P0=PL[0], P1=PL[1], P2=PL[2], P3=PL[3])
        execute(mutant, regs, Kn.mem)
        wrong += Kn.get(0x3000) != a * b * Ri % p
    out["flag_rule_mutant_detected"] = wrong > 0

    Kn.executed, Kn.calls = 0, dict(mul=0, add=0, sub=0)
    a = rng.randrange(1, p)
    Kn.put(0x1000, a * 2 ** 256 % p)
    Kn.inv(p, 0x1000, 0x2000)
    out["inversion_correct"] = Kn.get(0x2000) * Ri % p == pow(a, -1, p)
    out["inversion_kernel_calls"] = dict(Kn.calls)
    out["inversion_instructions"] = Kn.executed

    E = Curve(PROD["q"], PROD["r"], PROD["h"], PROD["a"], PROD["b"])
    G = (PROD["gx"], PROD["gy"])
    nb = 4 * max(1, prod_checks)
    pts = [E.mul(rng.randrange(1, E.r), G) for _ in range(nb)]
    tabs = [E.mul(rng.randrange(1, E.r), G) for _ in range(nb)]
    dbl = [i % 4 == 3 for i in range(nb)]
    Kn.executed, Kn.calls = 0, dict(mul=0, add=0, sub=0)
    res = batch_step(Kn, p, pts, tabs, dbl)
    want = [E.add(pts[i], pts[i]) if dbl[i] else E.add(pts[i], tabs[i]) for i in range(nb)]
    n_d = sum(dbl)
    n_a = nb - n_d
    cnt = Kn.calls
    out["batch"] = dict(points=nb, additions=n_a, doublings=n_d, results_equal_reference=res == want,
                        kernel_calls=dict(cnt), executed_instructions=Kn.executed)
    out["batch_kernel_calls_as_stated"] = (cnt["mul"] == 6 * n_a + 7 * n_d - 3 + 382 and
                                           cnt["add"] == 3 * n_d and cnt["sub"] == 6 * n_a + 5 * n_d)
    Mi, Ai = len(Kn.code["mul"]), max(len(Kn.code["add"]), len(Kn.code["sub"]))
    out["executed_equals_listing_counts"] = Kn.executed == (cnt["mul"] * Mi + cnt["add"] * len(Kn.code["add"]) +
                                                           cnt["sub"] * len(Kn.code["sub"]))
    out["field_instructions_per_point_excluding_inversion_le"] = dict(addition=6 * Mi + 6 * Ai, doubling=7 * Mi + 8 * Ai)
    out["ledger"] = ledger(Mi, Ai, 2 ** 17, 128, 127.72)
    out["part_c_seconds"] = round(time.time() - T0, 2)
    return out


# ---------------------------------------------------------------- Part B: toy run
# ---------------------------------------------------------------- the declared report store (toy scale)
class ReportStore:
    """The report store of proof.md section 3 at toy scale: the same structure,
    smaller numbers.

    One record (key, a, b, trail) per report, appended to the tail of one of
    `nodes` logs. The log is chosen by the report's arrival index (index mod
    nodes) and never by its key, so there is no slot index, no home address, no
    probe sequence, no replacement rule, no eviction, no load factor, and the
    occupancy of every log is a function of the report count alone. Detection is
    exact: each log is LSD radix-sorted on the full key with one counter array,
    the sorted logs are merged with a heap over their heads comparing full keys,
    and records with equal keys come out adjacent. A log at its declared capacity
    is the only failure the store has.

    Scope: this tests the store's structure and its exactness at toy scale (keys
    below 2^30, logs of a few thousand records, 16 logs, 8-bit radix digits). It
    does not test the production capacity (2^64 records), the production 32-bit
    radix digits, or the production instruction counts; those are analytic and the
    toy makes no claim about them."""

    def __init__(self, nodes, cap_per_log, key_bits):
        self.nodes, self.cap, self.key_bits = nodes, cap_per_log, key_bits
        self.logs = [[] for _ in range(nodes)]
        self.index = 0
        self.records = 0
        self.max_log = 0
        self.full = False

    def add(self, key, a, b, trail):
        i = self.index
        self.index += 1
        log = self.logs[i % self.nodes]
        if len(log) >= self.cap:
            self.full = True
            return False
        log.append((key, a, b, trail))
        self.records += 1
        if len(log) > self.max_log:
            self.max_log = len(log)
        return True

    def _sorted(self, log):
        ndigits = (self.key_bits + 7) // 8
        src, buf = list(log), [None] * len(log)
        for p in range(ndigits):
            sh = 8 * p
            cnt = [0] * 256
            for rec in src:
                cnt[(rec[0] >> sh) & 255] += 1
            tot = 0
            for d in range(256):
                c = cnt[d]
                cnt[d] = tot
                tot += c
            for rec in src:
                d = (rec[0] >> sh) & 255
                buf[cnt[d]] = rec
                cnt[d] += 1
            src, buf = buf, src
        return src

    def collisions(self, limit=8):
        """Merge the sorted logs by full key and return up to `limit` groups of two
        or more records with equal keys. Equal keys are exactly the class
        collisions: no fingerprint, no truncation, no hash takes part."""
        logs = [self._sorted(log) for log in self.logs if log]
        heap = [(log[0][0], j, 0) for j, log in enumerate(logs)]
        heapq.heapify(heap)
        groups, cur_key, cur = [], None, []
        while heap:
            key, j, pos = heapq.heappop(heap)
            rec = logs[j][pos]
            if cur_key is not None and key == cur_key:
                cur.append(rec)
            else:
                if len(cur) > 1:
                    groups.append(cur)
                    if len(groups) >= limit:
                        return groups
                cur_key, cur = key, [rec]
            if pos + 1 < len(logs[j]):
                heapq.heappush(heap, (logs[j][pos + 1][0], j, pos + 1))
        if len(cur) > 1:
            groups.append(cur)
        return groups


def solve_pair(E, W, T, point, r, first, second):
    """Re-walk two trails from their stored (a, b) and solve x from the two
    coefficient pairs. True when [x]P = T, None when this pair cannot be used."""
    s1 = dict(W.stats)
    _, K1, _, c1, (A1, B1) = W.run_trail(point(first[0], first[1]), first)
    _, K2, _, c2, (A2, B2) = W.run_trail(point(second[0], second[1]), second)
    W.stats.update(s1)  # the re-walk is a separate charge in the production ledger
    if K1 != K2 or c1 != c2 or (B1 - B2) % r == 0:
        return None
    x = (A2 - A1) * pow(B1 - B2, -1, r) % r
    return True if E.fmul(x) == T else None


def toy_trial(E, rng, R, dp_bits):
    """One trial. pipeline = True applies the production resource rules scaled to the
    toy (phase-1 budget B in class visits, phase 2, the declared report store with
    its capacity, at most 8 processed matches) and finds the collision by building,
    sorting and merging the store. pipeline = False is the T_c census: it walks until
    two trails reach the same class and reports that visit count; its `first_seen`
    and `census` dictionaries are arrival instrumentation, not the store."""
    r = E.r
    xs = rng.randrange(1, r)
    T = E.fmul(xs)

    def point(a, b):
        # [a]G + [b]T, evaluated as [a + b*xs]G with the fixed-base table (same point;
        # knowledge of xs is used only to speed up the harness, never by the walk)
        return E.fmul(a + b * xs)

    W = Walker(E, R, dp_bits, T, rng, point)
    N = (r - 1) // CLS
    first_seen = {}  # class key -> trail id (instrumentation only)
    state = {"visits": 0, "tc": None, "trail": None}

    def on_class(K):
        # every class arrival is one visit: trail starts, ordinary steps, doublings (escapes),
        # and repeated arrivals inside fruitless cycles (all are counted, none skipped)
        state["visits"] += 1
        tid = state["trail"]
        prev = first_seen.get(K)
        if prev is None:
            first_seen[K] = tid
        elif prev != tid and state["tc"] is None:
            state["tc"] = state["visits"]

    Np = N * R / (R - 1)
    B = 1.05 * math.sqrt(2 * Np * math.log(1 / 0.605))
    D_cap = int(2 * (B + W.Lmax) / (1 << dp_bits)) + 1
    NODES = cfg.get("nodes", 16)
    store = ReportStore(NODES, -(-D_cap // NODES), E.q.bit_length())
    census = {}          # T_c mode only: instrumentation, first (a, b) seen per reported class
    trail_id = 0
    solved = None
    # Trails are run one after another (a fresh trail after every distinguished class,
    # cap or exceptional addition). Interleaving concurrent walks would change the
    # order of visits (individual T_c values) but not their distribution under H1.
    while solved is None:
        if pipeline and state["visits"] >= B:
            break                              # phase 1 ended; the open trail is finished
        a, b = rng.randrange(r), rng.randrange(r)
        P = point(a, b)
        if P is None:
            continue
        trail_id += 1
        state["trail"] = trail_id
        status, K, length, c, _ = W.run_trail(P, None, on_class)
        if status != "dp":
            continue                           # cap, traversal abandon or exceptional: no report
        if pipeline:
            if not store.add(K, a, b, trail_id):   # a log at its capacity: declared failure
                solved = False
                break
        else:
            seen = census.get(K)
            if seen is not None:
                got = solve_pair(E, W, T, point, r, (a, b), seen)
                if got is not None:
                    solved = got
                    break
            else:
                census[K] = (a, b)
        if time.time() - T0 > budget:
            return None
    if pipeline and solved is None:
        solved = False
        tries = 0
        for grp in store.collisions(limit=8):
            base = grp[0]
            for other in grp[1:]:
                tries += 1
                if tries > 8:
                    break
                got = solve_pair(E, W, T, point, r, (base[1], base[2]), (other[1], other[2]))
                if got is not None:
                    solved = got
                    break
            if solved is not None:
                break
    nl = dict(store_records=store.records, store_logs=NODES, store_cap_per_log=store.cap,
              store_max_log=store.max_log, store_max_log_over_cap=store.max_log / store.cap,
              store_log_overflows=1 if store.full else 0)
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
        W.esc_row = [W.table[0][t] for t in range(len(W.units))]
        W.table = [W.table[R - 1]] + W.table[1:]  # M_esc not used as a partition multiplier
    out = dict(trails=0, operations=0, escapes=0, nonfruitless_returns=0, fruitless_repeats=0,
               trails_with_nonfruitless_returns=0)
    N = (r - 1) // 2
    for _ in range(ntrails):
        if time.time() - T0 > budget:
            break
        P = E.fmul(rng.randrange(1, r))
        seen = {}
        cnt = [0, 0]  # [non-fruitless returns, fruitless-cycle repeats]
        arr = [0]

        def on_class(K):
            # a repeat inside a fruitless cycle of length <= 64 (same escape epoch, at
            # most 64 arrivals after the previous visit) is expected and counted in T_c;
            # any other repeat is a "non-fruitless return" (the quantity H1 bounds)
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


def wilson(k, n, z=1.96):
    """Two-sided 95% Wilson score interval for k successes in n trials."""
    if not n:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(max(0.0, c - h), 4), round(min(1.0, c + h), 4)]


def main():
    out = {"config": cfg}
    if prod_checks:
        out["part_a"] = part_a(prod_checks)
    if listing_checks:
        out["part_c"] = part_c(listing_checks)
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
    out["toy"] = dict(q=E.q, r=r, h=E.h, b=E.b, N=N, R=R, dp_bits=dp_bits,
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
        k_ok = sum(1 for x in res if x["ok"])
        summ = dict(
            trials_requested=trials,
            trials_completed=n,
            trial_count_fixed=(n == trials),   # the stop condition is the count, not the clock
            mode="pipeline (scaled budget, phase 2, declared store, 8-match cap)" if pipeline
            else "T_c census, run until the first inter-trail collision",
            pipeline_success_rate=round(k_ok / n, 4) if pipeline else None,
            pipeline_success_wilson95=wilson(k_ok, n) if pipeline else None,
            all_solved_correctly=None if pipeline else all(x["ok"] for x in res),
            mean_tc_over_sqrtN=None if cens else round(sum(tcs) / len(tcs) / math.sqrt(N), 4),
            stderr_mean_tc_over_sqrtN=None if cens else round(math.sqrt(sum((t / math.sqrt(N) - sum(tcs) / len(tcs) / math.sqrt(N)) ** 2 for t in tcs) / max(1, len(tcs) - 1) / len(tcs)), 4),
            model_mean_over_sqrtN=round(math.sqrt(math.pi * Np / 2) / math.sqrt(N), 4),
            frac_tc_le_S0=round(sum(t <= S0 for t in tcs) / len(tcs), 4) if tcs else None,
            frac_tc_le_S1=round(sum(t <= S1 for t in tcs) / len(tcs), 4) if tcs else None,
            binomial_stderr_at_0p39=round(math.sqrt(0.39 * 0.61 / n), 4),
            model_P_S0=0.395,
            H1_floor_P_S1=0.395,
            model_P_S1=round(1 - math.exp(-S1 * S1 / (2 * Np)), 4),
            mean_solve_over_tc=None if cens else round(sum(x["solve_visits"] / x["tc"] for x in res if x["tc"]) / len(tcs), 3),
            visits_per_step=round(sum(x["solve_visits"] for x in res) / max(1, sum(x["steps"] + x["doublings"] for x in res)), 4),
            total_steps=steps,
            doublings_per_step=sum(x["doublings"] for x in res) / max(1, steps),
            esc2_events=sum(x["esc2"] for x in res),
            esc2_expected_events=round(steps / (2 * R), 1),
            esc2_per_step=sum(x["esc2"] for x in res) / max(1, steps),
            esc2_model_1_over_2R=1 / (2 * R),
            checkpoint_escapes_per_step=sum(x["esccp"] for x in res) / max(1, steps),
            abandoned_trails=sum(x["abandoned"] for x in res),
            exceptional_additions=sum(x["exceptional"] for x in res),
            store_records=sum(x["store_records"] for x in res),
            store_log_overflows=sum(x["store_log_overflows"] for x in res),
            store_max_log_over_cap=max(x["store_max_log_over_cap"] for x in res) if pipeline else None,
            store_cap_per_log=res[0]["store_cap_per_log"] if pipeline else None,
            store_logs=res[0]["store_logs"] if pipeline else None,
        )
        out["toy_summary"] = summ
        if dump:
            out["tc_list"] = tcs
    out["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
