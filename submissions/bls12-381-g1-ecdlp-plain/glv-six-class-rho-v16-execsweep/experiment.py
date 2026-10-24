"""Six-class (negation + GLV) parallel Pollard rho on BLS12-381 G1, parity key.

Part A checks, on BLS12-381 G1 itself, every algebraic identity the attack uses
(including the truncated-product canonical-index rule of proof.md section 1.1,
against the exact rule on random and adversarial near-tie inputs):
the GLV pair (beta, lambda), the size-6 orbits, the parity class key
K = y if y is even else p - y (invariant under y -> p - y, hence under phi, so a
complete class key; bit 0 of K is identically 0), the canonical-member rule, the
twisted-table step identity, coefficient tracking along real production walk
steps, the doubling escape as a class map, and the collision-solving formula
(including its sign) on two production trails that are made to merge. The older
exact |y| = min(y, p - y) key is kept as the reference rule and checked to
induce the same partition of those orbits.

Part B runs the complete algorithm of proof.md (twisted R-adding class walk,
2-cycle escape, checkpoint escape, trail cap, distinguished classes, trail
restart after every distinguished class, re-walk solving, verification) on the
prime-order subgroup of a toy curve y^2 = x^3 + 4 over F_q, q = 1 mod 3, and
measures the number of class visits (trail starts + walk operations) until the first collision between
different trails, against the random-mapping model on N = (r'-1)/6 classes
(or (r'-1)/classes for the classes = 1, 2 control walks).
"sweep": true reports, for the same walk, the empirical T_c distribution and the
fractions P[T_c <= kappa S0] on the kappa grid of proof.md section 4 (the
measurement that chooses the H1 walk margin), with the declared sweep table.
Toy evidence only: it does not validate production cost bounds.
"""
import json, math, random, sys, time

T0 = time.time()
cfg = json.load(sys.stdin)
curve = cfg.get("curve", "q26")
trials = cfg.get("trials", 100)
R = cfg.get("R", 256)
dp_bits = cfg.get("dp_bits", 4)
seed = cfg.get("seed", 1)
prod_checks = cfg.get("production_checks", 2)
budget = cfg.get("seconds", 14)
dump = cfg.get("dump", False)
# pipeline = True: also apply the production resource rules scaled to the toy (phase-1
# budget, phase 2, report-table capacity, 8-match cap) and report the success rate.
pipeline = cfg.get("pipeline", False)
# sweep = True: report the measured T_c distribution and the per-kappa fractions
run_sweep = cfg.get("sweep", False)
# class key of the toy walk: "parity" (the attack) or "abs" (the reference rule
# |y| = min(y, p - y) of the earlier versions, kept as a control)
class_key = cfg.get("class_key", "parity")
KAPPA_GRID = (1.05, 1.04, 1.03, 1.02, 1.015)
GATE = 0.40  # the measured P[T_c <= kappa S0] must reach this at the chosen kappa
KAPPA = 1.02  # H1 walk margin chosen by the sweep (proof.md section 4)
H1_FLOOR = 0.3903  # H1 floor at kappa: P(T_c <= S_ops) >= 1 - exp(-(S_ops/kappa)^2/(2N'))
# long-trail trap test (0 = off): trails of long_len operations without distinguished points
long_trails = cfg.get("long_trails", 0)
long_len = cfg.get("long_len", 4096)
long_escape_add = cfg.get("long_escape_add", False)
# section 1.1 check sizes (defaults: 15,000 random and 5,000 stress inputs per production check)
canon_random = cfg.get("canon_random", None)
canon_stress = cfg.get("canon_stress", None)
# toy = False skips Part B (e.g. to spend the whole run on Part A checks)
run_toy = cfg.get("toy", True)
# control option: 6 = the attack (default); 1 = plain R-adding walk on points (no
# automorphism); 2 = negation-only classes. Used only to compare toy collision times.
CLS = cfg.get("classes", 6)
rng = random.Random(seed)

# ---------------------------------------------------------------- curves
# Production profile constants (tracks/bls12-381-g1-ecdlp-plain/target-profile.json).
PROD = dict(
    q=int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab", 16),
    r=int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16),
    h=int("396c8c005555e1568c00aaab0000aaab", 16),
    gx=int("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb", 16),
    gy=int("08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1", 16),
    beta=int("1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac", 16),
    lam=int("ac45a4010001a40200000000ffffffff", 16),
)
# Toy curves y^2 = x^3 + 4 over F_q, q = 1 mod 3, #E = h * r' with r' prime.
# (q, #E, h) found offline by trial; generator, beta, lambda derived below.
TOYS = {
    "q26": (58493047, 58489293, 3),
    "q30": (861101287, 861140127, 3),
    "q34": (10146855109, 10146690129, 3),
}


class Curve:
    """Affine arithmetic on y^2 = x^3 + 4 over F_q (points are tuples, None = infinity)."""

    def __init__(self, q, r, h):
        self.q, self.r, self.h = q, r, h

    def on_curve(self, P):
        return P is None or (P[1] * P[1] - P[0] ** 3 - 4) % self.q == 0

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
        y = sqrt_mod(x ** 3 + 4, q)
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
        # the parity key has bit 0 identically 0, so the partition index and the
        # distinguished test read bits 1.. of K (control walks do not shift)
        self.shift = 1 if (CLS == 6 and class_key == "parity") else 0
        r = E.r
        lam = E.lam
        # units u = (-1)^s * lam^j, index t = 3*s + j
        self.units = [(-1) ** s * pow(lam, j, r) % r for s in (0, 1) for j in (0, 1, 2)]
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
            for t in range(6):
                s, j = divmod(t, 3)
                # u_t^{-1} M = (-1)^s phi^{-j}(M): one multiplication by beta^{-j}, no scalar mult
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
        """Return (K, t, canonical point): canonical = units[t] * P.

        The attack key is the PARITY key K = y if y is even else q - y. q is
        odd, so exactly one of y, q - y is even: K is invariant under y -> -y
        and unchanged by phi(x, y) = (beta x, y), so it is a complete class key
        (equal K gives y' = +-y and x'^3 = x^3, hence x' = beta^j x), and bit 0
        of K is identically 0. The canonical member is c = (beta^j x, K) with
        j = argmin beta^j x, so c = [u_t]P with t = 3s + j and s = [y odd].
        class_key = "abs" runs the earlier reference rule K = min(y, q - y)
        instead (a control that proof.md section 1 quotes)."""
        if class_key == "abs" and CLS == 6:
            return self.canon_ref(P)
        q, beta = self.E.q, self.E.beta
        x, y = P
        if CLS == 1:  # control: point itself, key injective on points
            return x + q * y, 0, P
        if CLS == 2:  # control: negation classes, key x
            if y & 1:
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
        if y & 1:
            ny = q - y
            return ny, 3 + j, (cx, ny)
        return y, j, (cx, y)

    def canon_ref(self, P):
        """The earlier reference key K = |y| = min(y, q - y) (control only)."""
        q, beta = self.E.q, self.E.beta
        x, y = P
        bx = beta * x % q
        b2x = (-x - bx) % q
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

    def idx(self, K):
        """partition index: bits 1..20 of K (bit 0 of K is always 0)"""
        return self.kshift(K) & (self.R - 1)

    def kshift(self, K):
        return (K >> 1) if self.shift else K

    def step_target(self, P, K, t):
        # P is any member; canonical c = u_t P; next class = class(c + M_h) = class(P + u_t^{-1} M_h)
        return self.table[self.idx(K)][t]

    def is_dp(self, K):
        """distinguished: the dp_bits above the partition index are zero"""
        return (self.kshift(K) >> self.rbits) & self.dpmask == 0

    def run_trail(self, P, coeffs=None, on_class=None):
        """Walk from P until a distinguished class or the trail cap.

        Returns (status, K, length, P_end, coeffs_end). If coeffs=(a,b) is given,
        coefficients with P = aG + bT are tracked (used only when re-walking).
        on_class(K) is called on every class arrival (instrumentation only)."""
        E, r = self.E, self.E.r
        cp = None  # checkpoint key
        length = 0
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
            if cp is not None and K == cp:
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
                cp = None
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


# ---------------------------------------------------------------- fast canonical index (production)
# proof.md section 1.1: decide j = argmin_j (beta^j x mod p) from a truncated fixed-point
# product, falling back to the exact rule on near-ties. Emulates the 64-bit limb algorithm.
M64 = (1 << 64) - 1
GQ = (PROD["beta"] << 448) // PROD["q"]  # floor(beta/p * 2^448): 7 limbs
GL = [(GQ >> (64 * j)) & M64 for j in range(7)]
PH0, PH1 = (PROD["q"] >> 256) & M64, PROD["q"] >> 320  # top 125 bits of p
PTOP = PROD["q"] >> 320
TOL = 16


def exact_j(x):
    q, b = PROD["q"], PROD["beta"]
    v = [x, b * x % q, b * b % q * x % q]
    return min(range(3), key=lambda j: v[j])


def fast_j(x):
    """Return j, or None for 'fall back to the exact rule'. Only 12 limb products are
    formed: the high halves (umulh) of the six column-5 products x_i*G_{5-i} and the low
    halves (mul) of the six column-6 products x_i*G_{6-i}; their sum mod 2^64 is L'."""
    xl = [(x >> (64 * i)) & M64 for i in range(6)]
    L = 0
    for i in range(6):
        L += (xl[i] * GL[5 - i]) >> 64
        L += (xl[i] * GL[6 - i]) & M64
    L &= M64  # true bits 384..447 of x*G mod 2^448 are in [L, L + 12]
    if L < 16 or L > M64 - 16:  # beta*x mod p within ~2^-59 p of 0 or p
        return None
    V1 = ((PH0 * L >> 64) + PH1 * L) >> 64  # 0 <= (beta x mod p)/2^320 - V1 < 3
    X = x >> 320  # top limb of x, error < 1
    S = X + V1
    if abs(S - PTOP) <= TOL:  # cannot decide whether x + beta x wraps mod p
        return None
    V2 = PTOP - S if S < PTOP else 2 * PTOP - S  # within 5 of (beta^2 x mod p)/2^320
    vals = [X, V1, V2]
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
    # adversarial inputs: two of {x, beta x, beta^2 x} (or x + beta x and p, or beta x and 0)
    # within 2^k of each other, k = 300..380, so that the approximate comparison is stressed
    mults = [b - 1, b * b - 1, b - b * b, 1 + b, b]
    for _ in range(n_stress):
        if time.time() - T0 > budget:
            break
        k = rng.randrange(300, 381)
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


# ---------------------------------------------------------------- Part A: production checks
def part_a(nchecks):
    out = {}
    E = Curve(PROD["q"], PROD["r"], PROD["h"])
    q, r = E.q, E.r
    E.beta, E.lam = PROD["beta"], PROD["lam"]
    G = (PROD["gx"], PROD["gy"])
    E.G = G
    out["P0_on_curve"] = E.on_curve(G)
    # E.mul reduces scalars mod r, so check [r]P0 = O via [r-1]P0 = -P0
    out["r_P0_is_infinity"] = E.mul(r - 1, G) == E.neg(G)
    out["q_mod_3"] = q % 3
    out["beta_cube_root_of_unity"] = (E.beta * E.beta + E.beta + 1) % q == 0
    out["lambda_cube_root_of_unity"] = (E.lam * E.lam + E.lam + 1) % r == 0
    out["phi_P0_equals_lambda_P0"] = E.mul(E.lam, G) == E.phi(G)
    out["points_x0_are_3_torsion_not_G1"] = E.add(E.add((0, 2), (0, 2)), (0, 2)) is None
    out["log2_r"] = math.log2(r)
    out["log2_N_classes"] = math.log2((r - 1) // 6)
    # multiplicative order of 2 mod r (doubling escapes cannot cancel: u*2^m != 1 for
    # every unit u and every m < ord/6), from the factorisation of r - 1
    rm1, fac = r - 1, {}
    for pr in [2, 3, 11, 19, 10177, 125527, 859267, 906349, 2508409, 2529403, 52437899, 254760293]:
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
    ok_orbit = ok_key = ok_canon = ok_twist = ok_double = ok_ref = True
    for _ in range(nchecks * 4):
        k = rng.randrange(1, r)
        P = E.mul(k, G)
        members = [E.mul(u * k % r, G) for u in W.units]
        # unit action matches (-1)^s phi^j
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
        keven = P[1] if P[1] % 2 == 0 else q - P[1]
        if len({c[0] for c in cans}) != 1 or cans[0][0] != keven:
            ok_key = False
        if (cans[0][0] & 1) or cans[0][0] != W.canon((P[0], q - P[1]))[0]:
            ok_key = False  # bit 0 identically 0; invariant under y -> -y
        # the earlier exact reference rule K = min(y, q - y) induces the same
        # partition: both keys are class-invariant and lie in the same +-y pair
        refs = [W.canon_ref(Mm) for Mm in members]
        if len({c[0] for c in refs}) != 1 or len({c[2] for c in refs}) != 1:
            ok_ref = False
        if refs[0][0] not in (cans[0][0], q - cans[0][0]):
            ok_ref = False
        if len({c[2] for c in cans}) != 1:
            ok_canon = False
        for Mm, (K, t, c) in zip(members, cans):
            if E.mul(W.units[t], Mm) != c:
                ok_canon = False
            # twisted step: class(Mm + u_t^-1 M_h) == class(c + M_h)
            h = W.idx(K)
            TP = W.table[h][t][0]
            Mh = W.table[h][0][0]
            n1 = W.canon(E.add(Mm, TP))
            n2 = W.canon(E.add(c, Mh))
            if n1[2] != n2[2]:
                ok_twist = False
            if W.canon(E.add(Mm, Mm))[2] != W.canon(E.add(c, c))[2]:
                ok_double = False
    out["orbit_size_6_and_unit_action"] = ok_orbit
    out["key_parity_even_representative_is_class_invariant"] = ok_key
    out["key_abs_y_reference_is_same_partition"] = ok_ref
    out["canonical_member_unique"] = ok_canon
    out["twisted_step_is_class_map"] = ok_twist
    out["doubling_is_class_map"] = ok_double
    ok_tab = all(W.table[h][t][0] == E.mul(W.uinv[t], W.table[h][0][0]) for h in range(W.R) for t in range(6))
    out["twisted_table_entries_equal_u_inverse_M"] = ok_tab

    # production walk with coefficient tracking, and a forced merge + solve
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
    # trail B starts at a random member of the class of A's point at step i0, with its
    # own random coefficients (constructed with knowledge of x: test harness only)
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

    dp_table = {}
    trail_id = 0
    solved = None
    # scaled production rules (pipeline mode): phase-1 budget in class visits
    # B = KAPPA*S0, S0 = sqrt(2 N' ln(1/0.6097)) class arrivals (trail starts are counted, which is conservative), report
    # capacity 2*theta*(B + L_max) + 1, at most 8 processed matches. With one walk,
    # phase 2 is the completion of the trail that is open when B is reached.
    Np = N * R / (R - 1)
    B = KAPPA * math.sqrt(2 * Np * math.log(1 / 0.6097))
    D_cap = int(2 * (B + W.Lmax) / (1 << dp_bits)) + 1
    matches = 0
    # Trails are run one after another (a fresh trail after every distinguished class,
    # cap or exceptional addition). Interleaving concurrent walks would change the
    # order of visits (individual T_c values) but not their distribution under H1.
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
            # re-walk both trails with coefficients
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
            dp_table[K] = (a, b, length)
        if time.time() - T0 > budget:
            return None
    return dict(tc=state["tc"], solve_visits=state["visits"], ok=solved, trails=trail_id, **W.stats)


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



SWEEP_TABLE = [
    # measured with this file's sweep mode (T_c mode, R = 2048, q26 toy of section 6):
    # dp, trials, seeds, trail/sqrt(N), mean T_c/sqrt(N) +- se, and
    # P[T_c <= kappa S0] +- binomial se on the kappa grid (1.05, 1.04, 1.03, 1.02, 1.015)
    dict(dp_bits=0, trials=22000, seeds=[7501, 7510], trail_over_sqrtN=0.000555,
         mean=1.2525, se=0.0044,
         fracs=[[0.4170, 0.0033], [0.4107, 0.0033], [0.4052, 0.0033], [0.3995, 0.0033], [0.3968, 0.0033]]),
    dict(dp_bits=1, trials=22000, seeds=[7601, 7610], trail_over_sqrtN=0.001110,
         mean=1.2533, se=0.0044,
         fracs=[[0.4160, 0.0033], [0.4101, 0.0033], [0.4047, 0.0033], [0.3989, 0.0033], [0.3961, 0.0033]]),
    dict(dp_bits=2, trials=44000, seeds=[7001, 7701], trail_over_sqrtN=0.002219,
         mean=1.2565, se=0.0031,
         fracs=[[0.4188, 0.0024], [0.4130, 0.0023], [0.4069, 0.0023], [0.4012, 0.0023], [0.3983, 0.0023]]),
    dict(dp_bits=3, trials=22000, seeds=[7101, 7110], trail_over_sqrtN=0.004438,
         mean=1.2541, se=0.0044,
         fracs=[[0.4209, 0.0033], [0.4146, 0.0033], [0.4079, 0.0033], [0.4013, 0.0033], [0.3982, 0.0033]]),
    dict(dp_bits="pooled", trials=110000, seeds=[7001, 7806], trail_over_sqrtN=None,
         mean=1.2546, se=0.0019,
         fracs=[[0.4183, 0.0015], [0.4123, 0.0015], [0.4063, 0.0015], [0.4004, 0.0015], [0.3975, 0.0015]]),
]

# ---------------------------------------------------------------- cost arithmetic (proof.md sections 3-5)
# The arithmetic consequences of the section-5 prices, with the rounding the proof uses
# (upper bounds rounded up, lower bounds rounded down), except the final ceiling, which is
# exact integer arithmetic: W_up^1000 <= 2^k and W_up^1000 > 2^(k-1) for the claim k/1000.
# The prices themselves rest on the section-5 instruction listing, not on this block.
CLAIM = (126241, 1000)  # time_log2 = 126.241 (the smallest thousandth covering W_up)
LEDGER = [("1.1 v1", 24), ("1.1 range test, V1, X, S, wrap test, V2", 15), ("1.1 decision", 14),
          ("parity bit, p - y limbs 0-3, select limbs 0-3 of K", 11),
          ("h (bits 1-20), distinguished test, twist index, table address", 8),
          ("loads of table point and state, store of state", 19), ("checkpoint test, refresh test, refresh", 12),
          ("length increment (69-bit), cap test, loop", 9), ("zero-denominator test", 1),
          ("per-batch meter charge", 1), ("batch inversion share", 1)]
SLACK = 2
C_INSTR, M_INSTR, ADDSUB_INSTR = 2000, 383, 32
OTHER_RANDOM_LOG2 = 74.4  # section 5: starts 2^62.2 * 2^12, fallbacks 2^72.2, words/attempts/probes 2^63
D_CAP = 1 << 63  # exact report capacity: 1.79 x the mean report count (Chernoff < 2^-60)
METER_MAX = 0.000199  # meter_3x ceiling that keeps success_3x >= 0.3901
LN_UP = (494788247, 10 ** 9)  # ceil(ln(1/0.6097) * 1e9) / 1e9, an upper bound
TRAVERSAL_EXTRA = 18  # copy the remembered member (12) + fingerprint comparison (6)
TRAVERSAL_SURCHARGE = 0.009  # 18 instructions / 2000


def up(v, d):
    return (math.floor(v * 10 ** d) + 1) / 10 ** d  # smallest d-decimal number above v


def dn(v, d):
    return (math.ceil(v * 10 ** d) - 1) / 10 ** d  # largest d-decimal number below v


def ceil_sqrt_frac(num, den):
    """Smallest s with s*s*den >= num (num, den > 0 integers): an exact ceil of sqrt(num/den)."""
    s = math.isqrt(num // den)
    while s * s * den < num:
        s += 1
    while s > 0 and (s - 1) * (s - 1) * den >= num:
        s -= 1
    return s


def chain(kappa, budget, claim_k):
    Rt = 1 << 20
    sub = sum(n for _, n in LEDGER)
    step_num = C_INSTR + budget                      # step price = step_num / 2000
    dbl_num = 7 * M_INSTR + 7 * ADDSUB_INSTR + budget  # doubling price = dbl_num / 2000
    # S_ops upper bound, exact: S_ops^2 = kappa^2 * 2 * N' * ln(1/0.6097)
    N = (PROD["r"] - 1) // 6
    kn, kd = int(round(kappa * 1000)), 1000
    num = kn * kn * 2 * N * Rt * LN_UP[0]
    den = kd * kd * (Rt - 1) * LN_UP[1]
    S_ops_up = ceil_sqrt_frac(num, den)
    S_ops = math.sqrt(num / den)
    S0 = S_ops / kappa
    ops_up = S_ops_up + (1 << 60)
    # expected random surcharges (units), from the section-4 rates
    E_dbl = (1 / (6 * Rt) + 2 ** -38) * ops_up
    E_trav = (2 / (6 * Rt) + 64 * 2 ** -37.99) * ops_up
    S_dbl = E_dbl * (dbl_num - step_num) / 2000
    S_trav = E_trav * TRAVERSAL_SURCHARGE
    S_other = 2 ** OTHER_RANDOM_LOG2
    random_sum = S_dbl + S_trav + S_other
    # hard extras (integer units, section 5)
    Lmax = 20 << 64
    phase2 = (1 << 20) * (Lmax + 65) * 119
    rewalk = 8 * 2 * ((Lmax + 65) * 119 + (1 << 12))
    reports = D_CAP // 2
    hash_init = math.ceil(2 ** 57.1)
    table = math.ceil(2 ** 32.2)
    blocks = math.ceil(2 ** 82.3)
    H_up = phase2 + rewalk + reports + hash_init + table + blocks + (1 << 61) + (1 << 13)
    A_up = math.ceil(3 * random_sum / METER_MAX)  # the meter's random allowance
    # exact integer work bound in 1/2000 units, and the exact ceiling k/1000
    W_2000 = step_num * ops_up + 2000 * (H_up + A_up)
    k = claim_k
    ceil_ok = W_2000 ** 1000 <= (2 ** k) * 2000 ** 1000
    ceil_tight = W_2000 ** 1000 > (2 ** (k - 1)) * 2000 ** 1000
    head_up = 2 ** (k / 1000) - (step_num * ops_up + 2000 * H_up) / 2000
    meter1 = up(random_sum / head_up, 6)
    meter3 = up(3 * random_sum / head_up, 6)
    success = dn(H1_FLOOR - meter1 - 0.00000011, 5)
    success3 = dn(H1_FLOOR - meter3 - 0.00000011, 5)
    return dict(sub=sub, budget=budget, kappa=kappa, step=step_num / 2000, dbl_num=dbl_num,
                dbl=dbl_num / 2000, S0=S0, S_ops=S_ops, S_ops_up=S_ops_up, ops_up=ops_up,
                log2_S0=math.log2(S0), log2_S_ops=math.log2(S_ops),
                log2_U1=math.log2(step_num * ops_up / 2000), log2_step=math.log2(step_num / 2000),
                E_dbl=E_dbl, S_dbl=S_dbl, E_trav=E_trav, S_trav=S_trav, S_other=S_other,
                random_sum=random_sum, log2_random_3x=math.log2(3 * random_sum),
                H_up=H_up, log2_H=math.log2(H_up), A_up=A_up, log2_A=math.log2(A_up),
                W_2000=W_2000, log2_W=math.log2(W_2000 / 2000), k=k, claim=k / 1000,
                headroom=head_up, log2_headroom=math.log2(head_up),
                meter1=meter1, meter3=meter3, success=success, success3=success3,
                ceil_ok=ceil_ok, ceil_tight=ceil_tight)


def cost():
    k = chain(KAPPA, sum(n for _, n in LEDGER) + SLACK, CLAIM[0])
    below = chain(KAPPA, sum(n for _, n in LEDGER) + SLACK, CLAIM[0] - 1)
    trav_max = k["budget"] + TRAVERSAL_EXTRA
    inv = (380 + bin(PROD["q"] - 2).count("1") - 1) * M_INSTR
    small = inv + 2000
    unbatched = inv + 2000 + 3 * M_INSTR + 6 * ADDSUB_INSTR + k["budget"] + 100
    return {
        "ledger_rows": [[name, n] for name, n in LEDGER],
        "ledger_subtotal": k["sub"], "slack": SLACK, "step_budget_instr": k["budget"],
        "price_step": "%.4f" % k["step"],
        "price_traversal": "%.4f" % (k["step"] + TRAVERSAL_SURCHARGE), "traversal_max_instr": trav_max,
        "price_doubling": "%.4f" % k["dbl"], "doubling_instr": k["dbl_num"],
        "doubling_surcharge": "%.4f" % (k["dbl"] - k["step"]),
        "fermat_inversion_instr": inv, "small_batch_inversion_instr_bound": small,
        "unbatched_operation_instr_bound": unbatched,
        "kappa": KAPPA, "kappa_grid": list(KAPPA_GRID), "kappa_gate": GATE,
        "log2_S0": "%.5f" % k["log2_S0"], "log2_S_ops": "%.5f" % k["log2_S_ops"],
        "S_ops_upper_integer": str(k["S_ops_up"]),
        "log2_U1_upper": "%.6f" % k["log2_U1"], "log2_step": "%.6f" % k["log2_step"],
        "log2_expected_doublings": "%.2f" % math.log2(k["E_dbl"]),
        "log2_doubling_surcharge_upper": "%.2f" % math.log2(k["S_dbl"]),
        "log2_traversal_steps_upper": "%.2f" % math.log2(k["E_trav"]),
        "log2_traversal_surcharge_upper": "%.2f" % math.log2(k["S_trav"]),
        "log2_other_random_upper": OTHER_RANDOM_LOG2,
        "log2_random_total_upper_3x": "%.2f" % k["log2_random_3x"],
        "D_cap_log2": 63, "D_cap": str(D_CAP), "expected_reports_log2_upper": "%.2f" % 62.17,
        "log2_hard_extras_upper": "%.2f" % k["log2_H"], "hard_extras_upper": str(k["H_up"]),
        "meter_allowance_upper": str(k["A_up"]), "log2_meter_allowance": "%.2f" % k["log2_A"],
        "W_upper_in_2000ths": str(k["W_2000"]), "log2_W_upper": "%.6f" % k["log2_W"],
        "claim": "%.3f" % k["claim"], "claim_fraction": list(CLAIM),
        "log2_headroom_lower": "%.2f" % k["log2_headroom"],
        "meter_upper": "%.6f" % k["meter1"], "meter_upper_3x_rate": "%.6f" % k["meter3"],
        "success_lower": "%.5f" % k["success"], "success_lower_3x_rate": "%.5f" % k["success3"],
        "checks": {
            "ledger_subtotal_115_budget_117": k["sub"] == 115 and k["budget"] == 117,
            "traversal_price_covers_listing": trav_max == 135 and abs((k["step"] + TRAVERSAL_SURCHARGE - 1) * 2000 - 135) < 1e-9,
            "doubling_price_is_listing_over_C": k["dbl_num"] == round(k["dbl"] * C_INSTR),
            "small_batch_inversion_le_118_units": small <= 118 * C_INSTR,
            "unbatched_operation_le_119_units": unbatched <= 119 * C_INSTR,
            "fermat_inversion_below_2_18": inv + 2000 < 1 << 18,
            "success_margin_at_least_1e-4": k["success"] >= 0.39 + 1e-4,
            "success_with_3x_cycle_rate_at_least_0.3901": k["success3"] >= 0.39 + 1e-4,
            "claim_is_exact_ceiling_thousandth": k["ceil_ok"] and k["ceil_tight"],
            "claim_is_minimal_thousandth": not (below["ceil_ok"]),
            "claim_above_phase1_plus_hard": k["claim"] > math.log2(2 ** k["log2_U1"] + 2 ** k["log2_H"]),
        },
    }


def main():
    out = {"config": cfg}
    out["cost"] = cost()
    pooled = [r for r in SWEEP_TABLE if r["dp_bits"] == "pooled"][0]
    out["sweep"] = dict(chosen_kappa=KAPPA, gate=GATE, grid=list(KAPPA_GRID), table=SWEEP_TABLE,
                        measured_pooled_at_kappa=pooled["fracs"][KAPPA_GRID.index(KAPPA)],
                        measured_pooled_at_1_015=pooled["fracs"][KAPPA_GRID.index(1.015)],
                        seats_above_floor_at_kappa=round((pooled["fracs"][KAPPA_GRID.index(KAPPA)][0] - H1_FLOOR)
                                                         / pooled["fracs"][KAPPA_GRID.index(KAPPA)][1], 1),
                        source="measured with this experiment's sweep mode (T_c mode: run until the "
                               "first inter-trail collision, every trial's x verified)")
    if prod_checks:
        out["part_a"] = part_a(prod_checks)
    if long_trails:
        out["long_trail_test"] = dict(curve=curve, R=R, long_len=long_len, escape="additive (control)" if long_escape_add else "doubling",
                                      **long_trail_test(make_toy(curve), rng, R, long_len, long_trails, long_escape_add))
    if run_sweep:
        E = make_toy(curve)
        N = (E.r - 1) // CLS
        Np = N * R / (R - 1)
        S0 = math.sqrt(2 * Np * math.log(1 / 0.6097))
        sq = math.sqrt(N)
        tcs = []
        verified = True
        for _ in range(trials):
            tr = toy_trial(E, rng, R, dp_bits)
            if tr is None:
                break
            verified = verified and bool(tr["ok"])
            if tr["tc"] is not None:
                tcs.append(tr["tc"])
        n = len(tcs)
        fr = {}
        for kap in KAPPA_GRID:
            c = sum(1 for t in tcs if t <= kap * S0)
            p = c / n
            fr["%.3f" % kap] = [round(p, 5), round(math.sqrt(p * (1 - p) / n), 5), c]
        mean = sum(tcs) / n / sq
        sd = math.sqrt(sum((t / sq - mean) ** 2 for t in tcs) / (n - 1))
        out["sweep_run"] = dict(curve=curve, seed=seed, R=R, dp_bits=dp_bits, trials=n,
                                class_key=class_key, all_verified=verified,
                                trail_over_sqrtN=round((1 << dp_bits) / sq, 6),
                                mean_tc_over_sqrtN=round(mean, 5), se_mean=round(sd / math.sqrt(n), 5),
                                model_mean_over_sqrtN=round(math.sqrt(math.pi * Np / 2) / sq, 5),
                                fracs=fr, gate_reached=GATE)
        out["seconds"] = round(time.time() - T0, 2)
        print(json.dumps(out, indent=1))
        return
    if not run_toy:
        out["seconds"] = round(time.time() - T0, 2)
        print(json.dumps(out, indent=1))
        return
    E = make_toy(curve)
    r = E.r
    N = (r - 1) // CLS
    Np = N * R / (R - 1)
    S0 = math.sqrt(2 * Np * math.log(1 / 0.6097))  # pure model: probability 0.3903
    S1 = KAPPA * S0  # claimed phase-1 budget of proof.md section 4 (H1: >= 0.3903)
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
            model_P_S0=0.3903,
            H1_floor_P_S1=0.3903,
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
        )
        out["toy_summary"] = summ
        if dump:
            out["tc_list"] = tcs
    out["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(out, indent=1))


main()
