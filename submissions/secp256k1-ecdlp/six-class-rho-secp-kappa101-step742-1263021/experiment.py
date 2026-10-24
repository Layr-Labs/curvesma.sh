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

Kernels mode ("kernels": true) checks the cost ledger of proof.md section 5: limb-for-limb Python
ports of the field kernels of kernels/secp_walk.c (fe_mul, fe_sqr, fe_sub, fe_canon, and the batched
affine addition built from them) against big integers and affine secp256k1 arithmetic, the instruction
counts parsed from the shipped AArch64 assembly kernels/secp_loops.s and kernels/secp_fe.s
(instructions per loop iteration of the single-block loops pass1 and pass2), the per-step budget, and
the work-meter ledger recomputed to 50 digits.
"""
import json, math, random, sys, time

KAPPA = 1.01   # slowdown allowance of H1 (proof.md section 4)
P0 = 0.3939086   # probability floor of H1 at S0 = S_ops / KAPPA

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
CLS = cfg.get("classes", 6)
kernels = cfg.get("kernels", False)
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
        return self.mul_raw(k, P)

    def mul_raw(self, k, P):
        """Scalar multiplication without reducing by a presumed point order."""
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
        if G is not None and E.mul_raw(r, G) is None:
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
    B = KAPPA * math.sqrt(2 * Np * math.log(1 / (1 - P0)))
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


# ---------------------------------------------------------------- kernels (proof.md section 5)
M64 = (1 << 64) - 1
C0 = (1 << 32) + 977                     # 2^256 = C0 (mod p)
KST = {}                                 # statistics of the limb ports (extreme values met)


def _kmax(key, v):
    if v > KST.get(key, -1):
        KST[key] = v


def _addc(a, b, c):
    s = a + b + c
    return s & M64, s >> 64


def _subc(a, b, c):
    d = a - b - c
    return d & M64, 1 if d < 0 else 0


def _inc(a, c):
    """a + c where the C code adds a carry bit without a carry out (asserted)."""
    assert a + c <= M64
    return a + c


def _limbs(v):
    return [(v >> (64 * i)) & M64 for i in range(4)]


def _val(l):
    return sum(x << (64 * i) for i, x in enumerate(l))


def _lohi(x, y):
    p = x * y
    return p & M64, p >> 64


def p_mul(a, b):
    """Limb-for-limb port of fe_mul in kernels/secp_walk.c (a, b: lists of 4 limbs)."""
    z = [0] * 8
    l, h = zip(*[_lohi(a[0], b[j]) for j in range(4)])
    z[0] = l[0]
    z[1], c = _addc(l[1], h[0], 0)
    z[2], c = _addc(l[2], h[1], c)
    z[3], c = _addc(l[3], h[2], c)
    z[4] = _inc(h[3], c)
    for i in range(1, 4):
        l, h = zip(*[_lohi(a[i], b[j]) for j in range(4)])
        p1, c = _addc(l[1], h[0], 0)
        p2, c = _addc(l[2], h[1], c)
        p3, c = _addc(l[3], h[2], c)
        p4 = _inc(h[3], c)
        z[i], c = _addc(z[i], l[0], 0)
        z[i + 1], c = _addc(z[i + 1], p1, c)
        z[i + 2], c = _addc(z[i + 2], p2, c)
        z[i + 3], c = _addc(z[i + 3], p3, c)
        z[i + 4] = _inc(p4, c)
    return _reduce(z)


def _reduce(z):
    g, k = zip(*[_lohi(z[4 + j], C0) for j in range(4)])
    q1, c = _addc(g[1], k[0], 0)
    q2, c = _addc(g[2], k[1], c)
    q3, c = _addc(g[3], k[2], c)
    q4 = _inc(k[3], c)
    s0, c = _addc(z[0], g[0], 0)
    s1, c = _addc(z[1], q1, c)
    s2, c = _addc(z[2], q2, c)
    s3, c = _addc(z[3], q3, c)
    s4 = _inc(q4, c)
    _kmax("s4", s4)
    assert s4 <= C0
    w0, w1 = _lohi(s4, C0)
    _kmax("w1", w1)
    assert w1 <= 1
    t0, c = _addc(s0, w0, 0)
    t1, c = _addc(s1, w1, c)
    t2, c = _addc(s2, 0, c)
    t3, c = _addc(s3, 0, c)
    m = C0 if c else 0
    KST["final_fold_taken"] = KST.get("final_fold_taken", 0) + (1 if c else 0)
    r0, c = _addc(t0, m, 0)
    r1, c = _addc(t1, 0, c)
    r2, c = _addc(t2, 0, c)
    r3 = _inc(t3, c)
    return [r0, r1, r2, r3]


def p_sqr(a):
    """Limb-for-limb port of fe_sqr."""
    a0, a1, a2, a3 = a
    l01, h01 = _lohi(a0, a1); l02, h02 = _lohi(a0, a2); l03, h03 = _lohi(a0, a3)
    l12, h12 = _lohi(a1, a2); l13, h13 = _lohi(a1, a3); l23, h23 = _lohi(a2, a3)
    o1 = l01
    o2, c = _addc(h01, l02, 0)
    o3, c = _addc(h02, l03, c)
    o4 = _inc(h03, c)
    p3 = l12
    p4, c = _addc(h12, l13, 0)
    p5 = _inc(h13, c)
    o3, c = _addc(o3, p3, 0)
    o4, c = _addc(o4, p4, c)
    o5, c = _addc(p5, l23, c)
    o6 = _inc(h23, c)
    d7 = o6 >> 63
    d6 = ((o6 << 1) & M64) | (o5 >> 63)
    d5 = ((o5 << 1) & M64) | (o4 >> 63)
    d4 = ((o4 << 1) & M64) | (o3 >> 63)
    d3 = ((o3 << 1) & M64) | (o2 >> 63)
    d2 = ((o2 << 1) & M64) | (o1 >> 63)
    d1 = (o1 << 1) & M64
    z0, g0 = _lohi(a0, a0); y2, g1 = _lohi(a1, a1); y4, g2 = _lohi(a2, a2); y6, g3 = _lohi(a3, a3)
    z1, c = _addc(d1, g0, 0)
    z2, c = _addc(d2, y2, c)
    z3, c = _addc(d3, g1, c)
    z4, c = _addc(d4, y4, c)
    z5, c = _addc(d5, g2, c)
    z6, c = _addc(d6, y6, c)
    z7, c = _addc(d7, g3, c)
    assert c == 0
    return _reduce([z0, z1, z2, z3, z4, z5, z6, z7])


def p_sub(a, b):
    """fe_sub: a any 256-bit value, b canonical."""
    assert _val(b) < PROD["q"]
    d0, bw = _subc(a[0], b[0], 0)
    d1, bw = _subc(a[1], b[1], bw)
    d2, bw = _subc(a[2], b[2], bw)
    d3, bw = _subc(a[3], b[3], bw)
    m = C0 if bw else 0
    r0, bw2 = _subc(d0, m, 0)
    r1, bw2 = _subc(d1, 0, bw2)
    r2, bw2 = _subc(d2, 0, bw2)
    r3, bw2 = _subc(d3, 0, bw2)
    assert bw2 == 0   # no second borrow: the wrapped value is above C0 whenever b is canonical
    return [r0, r1, r2, r3]


def p_canon(r):
    t0, c = _addc(r[0], C0, 0)
    t1, c = _addc(r[1], 0, c)
    t2, c = _addc(r[2], 0, c)
    t3, c = _addc(r[3], 0, c)
    return [t0, t1, t2, t3] if c else list(r)


def _kernel_port_check(nvec):
    q = PROD["q"]
    rr = random.Random(777)
    edge = [0, 1, 2, q - 1, q, q + 1, 2**256 - 1, 2**256 - 2, 2**255, 2**192, 2**64 - 1, 2**128 - 1,
            2**33, C0, q - 2**33, 2**256 - 2**33, 2**128, 2**64, (2**64 - 1) << 192]
    vals = edge + [rr.getrandbits(256) for _ in range(nvec)] + [q - rr.getrandbits(40) for _ in range(nvec // 8)] \
        + [2**256 - 1 - rr.getrandbits(40) for _ in range(nvec // 8)] + [rr.getrandbits(64) << (64 * rr.randrange(4)) for _ in range(nvec // 8)]
    bad = 0
    n = 0
    for a in edge:
        for b in edge:
            r = _val(p_mul(_limbs(a), _limbs(b)))
            n += 1
            bad += not (r < 2**256 and r % q == a * b % q)
    for _ in range(nvec):
        a, b = rr.choice(vals), rr.choice(vals)
        r = _val(p_mul(_limbs(a), _limbs(b)))
        n += 1
        bad += not (r < 2**256 and r % q == a * b % q)
    ns = 0
    for a in vals:
        r = _val(p_sqr(_limbs(a)))
        ns += 1
        bad += not (r < 2**256 and r % q == a * a % q)
    nsub = 0
    for _ in range(nvec):
        a, b = rr.choice(vals), rr.choice(vals) % q
        r = _val(p_sub(_limbs(a), _limbs(b)))
        nsub += 1
        bad += not (r < 2**256 and r % q == (a - b) % q)
        # when a is canonical the result is canonical
        if a < q:
            bad += not (r < q)
        c = _val(p_canon(_limbs(rr.choice(vals))))
        bad += not (c < q)
    for a in vals[:2000]:
        bad += _val(p_canon(_limbs(a))) != a % q
    return dict(mul_checked=n, sqr_checked=ns, sub_checked=nsub, mismatches=bad, extremes=dict(KST))


def _batch_port_check(nwalks):
    """Port of pass1/pass2 (shared inversion by pow) on production secp256k1 points."""
    q = PROD["q"]
    rr = random.Random(4242)
    Gp = (PROD["gx"], PROD["gy"])

    def aadd(A, B):
        if A[0] == B[0]:
            lam = 3 * A[0] * A[0] * pow(2 * A[1], -1, q) % q
        else:
            lam = (B[1] - A[1]) * pow(B[0] - A[0], -1, q) % q
        x = (lam * lam - A[0] - B[0]) % q
        return (x, (lam * (A[0] - x) - A[1]) % q)

    def amul(k, A):
        R = None
        while k:
            if k & 1:
                R = A if R is None else aadd(R, A)
            A = aadd(A, A)
            k >>= 1
        return R

    tab = [amul(rr.randrange(1, PROD["r"]), Gp) for _ in range(4)]
    walks = [amul(rr.randrange(1, PROD["r"]), Gp) for _ in range(nwalks)]
    tps = []
    for w in walks:
        t = tab[rr.randrange(4)]
        if t[0] == w[0]:
            t = tab[(tab.index(t) + 1) % 4]
        if rr.random() < 0.5:
            t = (t[0], q - t[1])
        tps.append(t)
    W = [(_limbs(x), _limbs(y)) for x, y in walks]
    T = [(_limbs(x), _limbs(y)) for x, y in tps]
    n = len(W)
    acc = [None] * n
    acc[0] = p_sub(T[0][0], W[0][0])
    for i in range(1, n):
        dx = p_sub(T[i][0], W[i][0])
        acc[i] = p_mul(acc[i - 1], dx)
    R = _limbs(pow(_val(acc[n - 1]), -1, q))

    def finish(w, t, inv):
        dy = p_sub(t[1], w[1])
        lam = p_mul(dy, inv)
        l2 = p_sqr(lam)
        x3 = p_sub(p_sub(l2, w[0]), t[0])
        x3 = p_canon(x3)
        tt = p_sub(w[0], x3)
        uu = p_mul(lam, tt)
        y3 = p_canon(p_sub(uu, w[1]))
        return (x3, y3)

    out = [None] * n
    for i in range(n - 1, 0, -1):
        inv = p_mul(R, acc[i - 1])
        dx = p_sub(T[i][0], W[i][0])
        R = p_mul(R, dx)
        out[i] = finish(W[i], T[i], inv)
    out[0] = finish(W[0], T[0], R)
    bad = 0
    for i in range(n):
        want = aadd(walks[i], tps[i])
        got = (_val(out[i][0]), _val(out[i][1]))
        bad += got != want
    return dict(walks=n, mismatches=bad)


def _asm_counts(text):
    """Per function: instructions, and for single-block loops the instructions per iteration.
    Directive lines (.text, .globl, ...) are skipped; labels are `name:` (functions), `.Lname:` or `N:`
    (local); a branch to `Nb` is a backward branch to the nearest preceding label N."""
    funcs = {}
    cur = None
    for line in text.split("\n"):
        if not line.strip():
            continue
        if not line[0].isspace():
            name = line.split(":")[0]
            if name.startswith(".L") or name.isdigit():
                funcs[cur]["items"].append(("label", name))
            else:
                cur = name
                funcs[cur] = {"items": []}
            continue
        if line.strip().startswith("."):
            continue
        funcs[cur]["items"].append(("ins", line.strip()))
    res = {}
    for name, f in funcs.items():
        items = f["items"]
        pos = {}
        idx = 0
        for kind, val in items:
            if kind == "label":
                pos.setdefault(val, []).append(idx)
            else:
                idx += 1
        total = idx
        loops = []
        idx = 0
        for kind, val in items:
            if kind == "ins":
                mn = val.split()[0]
                tgt = val.split()[-1]
                key = tgt[:-1] if (tgt[:-1].isdigit() and tgt.endswith("b")) else tgt
                if (mn in ("b", "cbz", "cbnz", "tbz", "tbnz") or mn.startswith("b.")) and key in pos and pos[key][0] <= idx:
                    start = max(p for p in pos[key] if p <= idx)
                    body = [v for k, v in items if k == "ins"][start:idx + 1]
                    inner = [b for b in body[:-1] if b.split()[0] in ("b", "bl", "cbz", "cbnz", "tbz", "tbnz", "ret") or b.split()[0].startswith("b.")]
                    loops.append(dict(label=tgt, iteration=len(body), inner_branches=len(inner)))
                idx += 1
        res[name] = dict(instructions=total, loops=loops,
                         loop_iteration_max=max([l["iteration"] for l in loops]) if loops else 0,
                         calls=sum(1 for k, v in items if k == "ins" and v.split()[0] == "bl"))
    return res


# ---- a small AArch64 interpreter: executes the shipped kernels/secp_loops.s (proof.md section 5.2)
def _a64_ops(s):
    ops, cur, depth = [], "", 0
    for ch in s:
        if ch == "[":
            depth += 1
        if ch == "]":
            depth -= 1
        if ch == "," and depth == 0:
            ops.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        ops.append(cur.strip())
    return ops


def a64_load(text):
    prog, labels = [], {}
    for raw in text.split("\n"):
        line = raw.split("//")[0].rstrip()
        if not line.strip():
            continue
        if not line[0].isspace():
            labels.setdefault(line.rstrip(":"), []).append(len(prog))
            continue
        ins = line.strip()
        if ins.startswith("."):
            continue
        mn, _, rest = ins.partition(" ")
        prog.append((mn, _a64_ops(rest)))
    return prog, labels


def a64_run(prog, labels, entry, regs, mem, limit=10**7):
    """Run from label `entry` until RET; returns the number of executed instructions. regs: dict x0..x30, sp."""
    M = M64
    R = {("x%d" % i): 0 for i in range(31)}
    R["sp"] = 0
    R.update(regs)
    flags = {"C": 0, "Z": 0}
    pc = labels[entry][0]
    steps = 0

    def rd(o):
        if o == "xzr":
            return 0
        return R[o]

    def wr(o, v):
        if o != "xzr":
            R[o] = v & M

    def imm(o):
        return int(o.lstrip("#"), 0)

    def ea(op, post):
        # returns (address, new_base_value or None, base_register)
        inner = op.strip()
        pre = inner.endswith("!")
        inner = inner.rstrip("!")[1:-1]
        parts = [p.strip() for p in inner.split(",")]
        base = parts[0]
        off = imm(parts[1]) if len(parts) > 1 else 0
        if post is not None:
            return rd(base), (rd(base) + imm(post)) & M, base
        if pre:
            return (rd(base) + off) & M, (rd(base) + off) & M, base
        return (rd(base) + off) & M, None, base

    while True:
        mn, o = prog[pc]
        steps += 1
        assert steps < limit
        pc += 1
        if mn in ("mul", "umulh"):
            a, b = rd(o[1]), rd(o[2])
            wr(o[0], (a * b) & M if mn == "mul" else (a * b) >> 64)
        elif mn in ("adds", "adc", "adcs", "add") and not (mn == "add" and (o[2].startswith("#") or len(o) == 4)):
            a, b = rd(o[1]), rd(o[2])
            cin = flags["C"] if mn in ("adc", "adcs") else 0
            t = a + b + cin
            wr(o[0], t)
            if mn in ("adds", "adcs"):
                flags["C"] = t >> 64
                flags["Z"] = 1 if (t & M) == 0 else 0
        elif mn == "add":     # add xd, xn, #imm | add xd, xn, xm, lsl #k
            if o[2].startswith("#"):
                wr(o[0], rd(o[1]) + imm(o[2]))
            else:
                wr(o[0], rd(o[1]) + (rd(o[2]) << imm(o[3].split()[1])))
        elif mn in ("sub", "subs", "sbc", "sbcs"):
            a = rd(o[1])
            b = imm(o[2]) if o[2].startswith("#") else rd(o[2])
            cin = flags["C"] if mn in ("sbc", "sbcs") else 1
            t = a + ((~b) & M) + cin
            wr(o[0], t)
            if mn in ("subs", "sbcs"):
                flags["C"] = t >> 64
                flags["Z"] = 1 if (t & M) == 0 else 0
        elif mn == "csel":
            cond = o[3]
            ok = (flags["C"] == 0) if cond in ("lo", "cc") else (flags["C"] == 1)
            wr(o[0], rd(o[1]) if ok else rd(o[2]))
        elif mn == "extr":
            wr(o[0], ((rd(o[1]) << 64 | rd(o[2])) >> imm(o[3])))
        elif mn == "lsr":
            wr(o[0], rd(o[1]) >> imm(o[2]))
        elif mn == "lsl":
            wr(o[0], rd(o[1]) << imm(o[2]))
        elif mn == "mov":
            wr(o[0], imm(o[1]) if o[1].startswith("#") else rd(o[1]))
        elif mn == "movk":
            sh = imm(o[2].split()[1]) if len(o) > 2 else 0
            wr(o[0], (rd(o[0]) & ~(0xFFFF << sh)) | ((imm(o[1]) & 0xFFFF) << sh))
        elif mn in ("ldp", "stp"):
            post = o[3] if len(o) > 3 else None
            addr, newb, base = ea(o[2], post)
            if mn == "ldp":
                wr(o[0], mem[addr])
                wr(o[1], mem[addr + 8])
            else:
                mem[addr] = rd(o[0])
                mem[addr + 8] = rd(o[1])
            if newb is not None:
                R[base] = newb
        elif mn == "ldr":
            post = o[2] if len(o) > 2 else None
            addr, newb, base = ea(o[1], post)
            wr(o[0], mem[addr])
            if newb is not None:
                R[base] = newb
        elif mn == "b.ne":
            if flags["Z"] == 0:
                tgt = o[0]
                assert tgt.endswith("b")      # numeric backward label
                cands = [p for p in labels[tgt[:-1]] if p <= pc - 1]
                pc = cands[-1]
        elif mn == "ret":
            return steps
        else:
            raise ValueError("unsupported instruction: " + mn)


def asm_exec_check(text, sizes=(2, 3, 5, 9, 17)):
    """Execute pass1/pass2 of the shipped assembly on production secp256k1 points and compare with affine
    arithmetic; the executed-instruction difference between batch sizes n and n+1 is the per-walk count."""
    q = PROD["q"]
    rr = random.Random(2468)
    Gp = (PROD["gx"], PROD["gy"])

    def aadd(A, B):
        if A[0] == B[0]:
            lam = 3 * A[0] * A[0] * pow(2 * A[1], -1, q) % q
        else:
            lam = (B[1] - A[1]) * pow(B[0] - A[0], -1, q) % q
        x = (lam * lam - A[0] - B[0]) % q
        return (x, (lam * (A[0] - x) - A[1]) % q)

    def amul(k, A):
        R = None
        while k:
            if k & 1:
                R = A if R is None else aadd(R, A)
            A = aadd(A, A)
            k >>= 1
        return R

    prog, labels = a64_load(text)
    tab = [amul(rr.randrange(1, PROD["r"]), Gp) for _ in range(3)]
    tab += [(x, q - y) for x, y in tab]
    res = dict(walks_checked=0, mismatches=0, executed={})
    for n in sizes:
        walks = [amul(rr.randrange(1, PROD["r"]), Gp) for _ in range(n)]
        tps = [tab[rr.randrange(len(tab))] for _ in range(n)]
        mem = {}
        WB, TB, PB, AB, RB, SB = 0x10000, 0x20000, 0x30000, 0x40000, 0x50000, 0x90000
        for i in range(n):
            for k, v in enumerate(_limbs(walks[i][0]) + _limbs(walks[i][1])):
                mem[WB + 64 * i + 8 * k] = v
            for k, v in enumerate(_limbs(tps[i][0]) + _limbs(tps[i][1])):
                mem[TB + 64 * i + 8 * k] = v
            mem[PB + 8 * i] = TB + 64 * i
        for a in range(0, 80 + 64, 8):
            mem[SB - 112 + a] = 0
        s1 = a64_run(prog, labels, "pass1", {"x0": n, "x1": WB, "x2": PB, "x3": AB, "sp": SB}, mem)
        last = sum(mem[AB + 64 * n + 32 + 8 * k] << (64 * k) for k in range(4))
        rinv = pow(last, -1, q)
        for k, v in enumerate(_limbs(rinv)):
            mem[RB + 8 * k] = v
        s2 = a64_run(prog, labels, "pass2", {"x0": n, "x1": WB, "x2": PB, "x3": AB, "x4": RB, "sp": SB}, mem)
        for i in range(n):
            got = (sum(mem[WB + 64 * i + 8 * k] << (64 * k) for k in range(4)),
                   sum(mem[WB + 64 * i + 32 + 8 * k] << (64 * k) for k in range(4)))
            res["walks_checked"] += 1
            res["mismatches"] += got != aadd(walks[i], tps[i])
        res["executed"][n] = dict(pass1=s1, pass2=s2)
    ns = sorted(res["executed"])
    res["per_extra_walk"] = dict(
        pass1=[res["executed"][b]["pass1"] - res["executed"][a]["pass1"] for a, b in zip(ns, ns[1:])],
        pass2=[res["executed"][b]["pass2"] - res["executed"][a]["pass2"] for a, b in zip(ns, ns[1:])])
    res["per_extra_walk_divided"] = dict(
        pass1=[(res["executed"][b]["pass1"] - res["executed"][a]["pass1"]) / (b - a) for a, b in zip(ns, ns[1:])],
        pass2=[(res["executed"][b]["pass2"] - res["executed"][a]["pass2"]) / (b - a) for a, b in zip(ns, ns[1:])])
    return res


# prices of proof.md section 5 in instructions (C = 1000 instructions per unit).  Each price is the
# derived bound itself: the ordinary step has no margin over its table, the doubling escape is its
# 950-instruction bound and the traversal step its 766-instruction bound (proof.md sections 0, 5.3).
PRICE_STEP, PRICE_TRAVERSAL, PRICE_DOUBLING = 742, 766, 950
# Analytic bookkeeping rows of proof.md section 5.3: the instructions the pseudocode of section 2
# performs per class arrival outside the two counted loops.  Each is the sum of an explicit sequence,
# with no allowance on top.  No row repeats an action of the counted 111 + 521: those loops load and
# store the walk point and consume the table pointer, and the zero-denominator test is charged here
# because the counted pass-1 iteration does not contain it.
ROWS = dict(F1_truncated_product=20, range_test_F0_S_F2=8, decision=13, sign_pminusy_K=17,
            h_dp_t_address=7, arrival_traffic=3, fingerprints_history=26, zero_denominator=4,
            length_cap_meter_loop=13, share_of_inversion_and_calls=1)


def kernels_check():
    import os
    from decimal import Decimal as D, getcontext
    getcontext().prec = 60
    here = os.path.dirname(os.path.abspath(__file__))
    res = {}
    res["python_port_vs_bigint"] = _kernel_port_check(2500)
    res["python_batch_port_vs_affine"] = _batch_port_check(24)
    asm = {}
    texts = {}
    for fname in ("secp_loops.s", "secp_fe.s"):
        for base in (os.path.join(here, "kernels"), "/snapshot/files/kernels", os.path.join(here, "files", "kernels")):
            path = os.path.join(base, fname)
            if os.path.exists(path):
                with open(path, encoding="utf-8") as fh:
                    texts[fname] = fh.read()
                asm.update(_asm_counts(texts[fname]))
                break
        else:
            res["asm"] = fname + " not found"
            return res
    res["asm_executed_vs_affine"] = asm_exec_check(texts["secp_loops.s"])
    res["asm_per_function"] = {k: dict(instructions=v["instructions"], loop_iteration_max=v["loop_iteration_max"],
                                       inner_branches=[l["inner_branches"] for l in v["loops"]], calls=v["calls"])
                               for k, v in sorted(asm.items())}
    loops_single_block = all(l["inner_branches"] == 0 for k in ("pass1", "pass2") for l in asm[k]["loops"]) \
        and len(asm["pass1"]["loops"]) == 1 and len(asm["pass2"]["loops"]) == 1
    arith = asm["pass1"]["loop_iteration_max"] + asm["pass2"]["loop_iteration_max"]
    rows = sum(ROWS.values())
    step = arith + rows
    mulf, sqrf = asm["fe_mul_f"]["instructions"], asm["fe_sqr_f"]["instructions"]
    glue = asm["fe_inv"]["instructions"] - asm["fe_inv"]["calls"]
    inv_bound = 255 * (sqrf + 1) + 15 * (mulf + 1) + glue + 255 * 4
    fixed_per_batch = (asm["pass1"]["instructions"] - asm["pass1"]["loop_iteration_max"]
                       + asm["pass2"]["instructions"] - asm["pass2"]["loop_iteration_max"] + 64)
    batch_overhead = inv_bound + fixed_per_batch
    res.update(
        loops_are_single_block=loops_single_block,
        arithmetic_instructions_per_step=arith, rows=ROWS, rows_total=rows, step_instructions_bound=step,
        fe_mul_f_instructions=mulf, fe_sqr_f_instructions=sqrf,
        inversion_instructions_bound=inv_bound, per_batch_fixed_instructions=batch_overhead,
        per_batch_overhead_below_2_pow_15=batch_overhead < 2**15,
        priced=dict(step=PRICE_STEP, traversal=PRICE_TRAVERSAL, doubling=PRICE_DOUBLING),
        traversal_instructions_bound=step + 24,
        asm_executed_matches_counts=bool(res["asm_executed_vs_affine"]["mismatches"] == 0 and set(res["asm_executed_vs_affine"]["per_extra_walk_divided"]["pass1"]) == {float(asm["pass1"]["loop_iteration_max"])} and set(res["asm_executed_vs_affine"]["per_extra_walk_divided"]["pass2"]) == {float(asm["pass2"]["loop_iteration_max"])}),
        price_covers_bound=bool(step <= PRICE_STEP and step + 24 <= PRICE_TRAVERSAL and batch_overhead < 2**15),
        margin_step_instructions=PRICE_STEP - step)

    # operation costs of the generated assembly (proof.md 5.1) and the identities that tie them to the loop counts
    OPC = dict(mul=12 + 3 * 17 + 28, sqr=12 + 9 + 7 + 15 + 28, sub=9, canon=8, add=9, add_weak_twice=14)
    ident1 = OPC["mul"] + OPC["sub"] + 5 + 4 + 2
    ident2 = 4 * OPC["mul"] + OPC["sqr"] + 5 * OPC["sub"] + 2 * OPC["canon"] + 15 + 6 + 2
    # doubling: 5 mul (acc, inv, R, lambda, lambda*(x-x3)), 2 sqr (x^2, lambda^2), dx = 2y (1 add), 3x^2 = v+v+v
    # (2 adds with a second carry fold), 4 sub, 2 canon, <= 72 loads/stores/pointers/loop control
    dbl_arith = 5 * OPC["mul"] + 2 * OPC["sqr"] + OPC["add"] + 2 * OPC["add_weak_twice"] + 4 * OPC["sub"] + 2 * OPC["canon"] + 72
    dbl_total = dbl_arith + rows + 50 + 30      # + bookkeeping rows, list moves, escape logic
    res["doubling"] = dict(op_costs=OPC, identity_pass1=ident1, identity_pass2=ident2,
                           identities_hold=bool(ident1 == asm["pass1"]["loop_iteration_max"] and ident2 == asm["pass2"]["loop_iteration_max"]),
                           arithmetic_bound=dbl_arith, total_bound=dbl_total, priced=PRICE_DOUBLING,
                           price_covers_bound=bool(dbl_total <= PRICE_DOUBLING))
    two = D(2)
    lg = lambda x: x.ln() / two.ln()
    r = PROD["r"]
    N = D((r - 1) // 6)
    Rr = D(2048)
    Np = N * Rr / (Rr - 1)
    kappa, p0 = D(str(KAPPA)), D(str(P0))
    S0 = (2 * Np * (D(1) / (D(1) - p0)).ln()).sqrt()
    Sops = kappa * S0
    step_u, dbl_u, trav_u = D(PRICE_STEP) / 1000, D(PRICE_DOUBLING) / 1000, D(PRICE_TRAVERSAL) / 1000
    ph1 = step_u * (Sops + two ** 60)
    # report store of proof.md section 3: integer D_cap, bucket count B, probe depth 8,
    # overflow list of B 8-byte addresses, report price 40 + 8*0.05
    from decimal import ROUND_CEILING
    Dcap = int((two ** D("63.83")).to_integral_value(rounding=ROUND_CEILING))
    B = -((-(2 * Dcap)) // 2 ** 16)
    report_price = D(40) + 8 * D("0.05")
    ovf_insert = 6 * B + 35
    ovf_lookup = 16 * B
    claim = D("126.3021")
    Wmax = two ** claim
    refills = int((Wmax / two ** 40).to_integral_value(rounding=ROUND_CEILING)) + 2 ** 20
    Lmax = 20 * two ** 64
    capped = dict(
        meter_remainders=two ** 61,
        meter_block_refills=D(refills) * 30 / 1000,
        phase2=two ** 20 * (Lmax + 65) * 600,
        rewalks=8 * 2 * ((Lmax + 65) * 600 + two ** 13),
        reports=D(Dcap) * report_price,
        overflow_worst_case=D(Dcap) * (ovf_insert + ovf_lookup) / 1000,
        store_init=D(Dcap) * 2 * 16 / 1000,
        table_verification=two ** 25,
        node_setup=two ** 26,
    )
    hard = sum(capped.values())
    rate = two ** D("-13.5681")        # upper bound on 1/(6R) + 2^-19.9954 = 2^-13.5681008
    random_items = dict(
        doubling=rate * (Sops + two ** 60) * (dbl_u - step_u),
        traversal=(trav_u - step_u) * 64 * two ** D("-19.99") * (Sops + two ** 60),
        starts=two ** 13 * (Sops / two ** 64 + two ** 20),
        other=two ** 70,
    )
    extras = sum(random_items.values())
    head = Wmax - ph1 - hard
    markov = extras / head
    allowance = (extras / D("0.0038883")).to_integral_value(rounding=ROUND_CEILING)
    w_ledger = ph1 + hard + allowance
    success = p0 - D("0.0038883") - D("2e-7")
    res["ledger"] = dict(
        log2_S0=float(lg(S0)), log2_Sops=float(lg(Sops)), log2_phase1=float(lg(ph1)),
        log2_hard_extras={k: round(float(lg(v)), 2) for k, v in capped.items()},
        log2_hard_total=float(lg(hard)),
        log2_random_items={k: round(float(lg(v)), 2) for k, v in random_items.items()},
        log2_expected_random_total=float(lg(extras)),
        log2_headroom=float(lg(head)), markov_early_halt_probability=float(markov),
        success_lower_bound=float(success),
        success_with_unrounded_ratio=float(p0 - markov - D("2e-7")),
        claimed_time_log2=float(claim),
        store=dict(buckets_per_node=B, bucket_bytes=104, overflow_address_bytes=8,
                   probe_instructions=39, probe_price=float(D("0.05")),
                   report_price=float(report_price), overflow_insert=ovf_insert,
                   overflow_lookup=ovf_lookup, blocks_per_refill=2 ** 40, refills=refills),
        log2_allowance=float(lg(allowance)), log2_w_ledger=float(lg(w_ledger)),
        w_ledger_below_claim=bool(w_ledger <= Wmax),
        claim_is_three_decimal_ceiling=bool(w_ledger <= Wmax and lg(w_ledger) > claim - D("0.001")),
        claim_above_phase1_plus_hard=bool(Wmax > ph1 + hard),
        success_at_least_0_39=bool(success > D("0.39")),
        success_at_least_0_39002=bool(success > D("0.39002")),
        model_probability_at_S_ops=float(1 - (1 - p0) ** (kappa * kappa)),
        hard_below_2_pow_107_156=bool(hard < two ** D("107.156")),
        random_below_2_pow_111_01=bool(extras < two ** D("111.01")),
        kappa=float(kappa), p0=float(p0),
        markov_below_0_0038883=bool(markov < D("0.0038883")),
    )
    return res



def main():
    out = {"config": cfg}
    if kernels:
        out["kernels"] = kernels_check()
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
    S0 = math.sqrt(2 * Np * math.log(1 / (1 - P0)))  # pure model: probability P0
    S1 = KAPPA * S0  # claimed phase-1 budget of proof.md section 4 (H1: >= P0)
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
            model_P_S0=P0,
            H1_floor_P_S1=P0,
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
