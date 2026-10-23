"""Negation-class (two-class) parallel Pollard rho on NIST P-256.

Adaptation of the class-walk construction of submission 1c600df4 (mpjunior92,
bls12-381-g1-ecdlp-plain, glv-six-class-rho-v11; ported to secp256k1 in our
a188e335) to P-256, which has no efficient endomorphism: classes {P, -P}, class
key x, doubling escapes, 2-cycle and checkpoint detection, re-walk solving.

Part A checks on P-256 itself: curve, order, order of 2 mod r, orbit size 2, the
x class key, the canonical member, the twisted (negated) table step as a class
map, the doubling class map, coefficient tracking on production steps and a
forced-merge solve. Part B runs the complete per-trail algorithm on prime-order
subgroups of toy curves y^2 = x^3 - 3x + b, measuring the class-visit count T_c
until the first collision between different trails against the random-mapping
model on N = (r'-1)/2 classes (classes = 1 runs the no-automorphism control).
Toy evidence only: it does not validate production cost bounds.

Kernels mode ("kernels": true) checks the cost ledger of proof.md section 5
(package negation-class-rho-counted): a pure-Python port of the 4x64-bit
Montgomery kernels of kernels/p256_walk.c against big integers, the instruction
counts parsed from the shipped AArch64 assembly kernels/p256_field.s and
kernels/p256_loops.s (instructions per loop iteration = longest path through the
loop), and the exact arithmetic of the work and success ledger.
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
# long-trail trap test (0 = off): trails of long_len operations without distinguished points
long_trails = cfg.get("long_trails", 0)
long_len = cfg.get("long_len", 4096)
long_escape_add = cfg.get("long_escape_add", False)
# toy = False skips Part B (e.g. to spend the whole run on Part A checks)
run_toy = cfg.get("toy", True)
# kernels = True: run the cost-ledger checks of proof.md section 5 (see kernels_check)
kernels = cfg.get("kernels", False)
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
    # B = 1.05*S0, S0 = sqrt(2 N' ln(1/0.605)) class arrivals (trail starts are counted, which is conservative), report
    # capacity 2*theta*(B + L_max) + 1, at most 8 processed matches. With one walk,
    # phase 2 is the completion of the trail that is open when B is reached.
    Np = N * R / (R - 1)
    B = 1.05 * math.sqrt(2 * Np * math.log(1 / 0.605))
    D_cap = int(2 * (B + W.Lmax) / (1 << dp_bits)) + 1
    # storage nodes (pipeline mode): reports are routed by the key bits just above the
    # partition and distinguished bits (production: bits 128..143 of K) to NODES nodes,
    # each with capacity ceil(2 * D_cap / NODES) as in production (2 D_cap slots in
    # total); a report to a full node is a failure, as in production.
    NODES = cfg.get("nodes", 16)
    node_shift = W.rbits + dp_bits
    node_cap = -(-2 * D_cap // NODES)
    node_load = [0] * NODES
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



# ------------------------------------------------------------ cost ledger
def _kernel_port_check(ncases):
    """Port of mulm/subm/addm of kernels/p256_walk.c with 64-bit limb semantics."""
    P = PROD["q"]
    M64 = (1 << 64) - 1
    P3 = 0xFFFFFFFF00000001
    R256 = 1 << 256
    Rinv = pow(R256, -1, P)
    pl = [(P >> (64 * i)) & M64 for i in range(4)]

    def limbs(x):
        return [(x >> (64 * i)) & M64 for i in range(4)]

    def val(l):
        return sum(l[i] << (64 * i) for i in range(len(l)))

    def mulm(a, b):
        t = [0] * 8
        c = 0
        for j in range(4):
            w = a[0] * b[j] + c
            t[j] = w & M64
            c = w >> 64
        t[4] = c
        for i in (1, 2, 3):
            c = 0
            for j in range(4):
                w = a[i] * b[j] + t[i + j] + c
                t[i + j] = w & M64
                c = w >> 64
            t[i + 4] = c
        top = 0
        for i in range(4):
            m = t[i]
            s1 = t[i + 1] + (m << 32)
            t[i + 1] = s1 & M64
            c = s1 >> 64
            s1 = t[i + 2] + c
            t[i + 2] = s1 & M64
            c = s1 >> 64
            s1 = t[i + 3] + c + m * P3
            t[i + 3] = s1 & M64
            c = s1 >> 64
            for j in range(i + 4, 8):
                s1 = t[j] + c
                t[j] = s1 & M64
                c = s1 >> 64
            top += c
        res = t[4:8]
        diff, bw = [], 0
        for i in range(4):
            d = res[i] - pl[i] - bw
            diff.append(d & M64)
            bw = 1 if d < 0 else 0
        return res if (bw and not top) else diff

    def subm(a, b):
        d, bw = [], 0
        for i in range(4):
            x = a[i] - b[i] - bw
            d.append(x & M64)
            bw = 1 if x < 0 else 0
        mk = M64 if bw else 0
        out, c = [], 0
        for i in range(4):
            x = d[i] + (pl[i] & mk) + c
            out.append(x & M64)
            c = x >> 64
        return out

    def addm(a, b):
        sm, c = [], 0
        for i in range(4):
            x = a[i] + b[i] + c
            sm.append(x & M64)
            c = x >> 64
        diff, bw = [], 0
        for i in range(4):
            d = sm[i] - pl[i] - bw
            diff.append(d & M64)
            bw = 1 if d < 0 else 0
        return sm if (bw and not c) else diff

    rr = random.Random(20260101)
    edge = [0, 1, 2, P - 1, P - 2, P - 3, (P + 1) // 2, 1 << 255, 1 << 224, 1 << 192, 1 << 96,
            1 << 64, (1 << 64) - 1, 1 << 32, R256 % P, Rinv, P - (1 << 224), P - (1 << 192)]
    bad = n = 0
    for x in edge:
        for y in edge:
            bad += val(mulm(limbs(x), limbs(y))) != x * y * Rinv % P
            bad += val(addm(limbs(x), limbs(y))) != (x + y) % P
            bad += val(subm(limbs(x), limbs(y))) != (x - y) % P
            n += 1
    for _ in range(ncases):
        x, y = rr.randrange(P), rr.randrange(P)
        bad += val(mulm(limbs(x), limbs(y))) != x * y * Rinv % P
        bad += val(addm(limbs(x), limbs(y))) != (x + y) % P
        bad += val(subm(limbs(x), limbs(y))) != (x - y) % P
        n += 1
    return dict(cases=n, mismatches=bad)


def _asm_counts(text):
    """Per function of a compiler-generated .s file: instruction count, whether
    the control-flow graph has a cycle, and for looping functions the longest
    path (in instructions) of one loop iteration, taken over all back edges.
    The longest path includes the cold slow-path block (its `bl` counts as one
    instruction; the callee is priced separately), so it bounds every iteration."""
    import re
    funcs, cur = {}, None
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("//"):
            continue
        m = re.match(r"^([A-Za-z_.$][\w.$]*):$", line)
        if m:
            name = m.group(1)
            if name.startswith(".L"):
                funcs[cur]["labels"][name] = len(funcs[cur]["ins"])
            else:
                cur = name
                funcs[cur] = dict(ins=[], labels={})
            continue
        funcs[cur]["ins"].append(line)
    out = {}
    for name, f in funcs.items():
        ins, labels = f["ins"], f["labels"]
        n = len(ins)
        succ = [[] for _ in range(n)]
        calls = rets = 0
        mix = {}
        for pos, line in enumerate(ins):
            op = line.split()[0]
            mix[op] = mix.get(op, 0) + 1
            calls += op in ("bl", "blr")
            rets += op == "ret"
            tgt = line.split()[-1].split(",")[-1]
            if op == "ret":
                continue
            if op == "b":
                succ[pos].append(labels[tgt])
                continue
            if op.startswith("b.") or op in ("cbz", "cbnz", "tbz", "tbnz"):
                succ[pos].append(labels[tgt])
            if pos + 1 < n:
                succ[pos].append(pos + 1)
        # depth-first search from the entry: back edges are edges into a node on the stack
        color, back = [0] * n, []
        stack = [(0, 0)]
        color[0] = 1
        while stack:
            v, k = stack.pop()
            if k < len(succ[v]):
                stack.append((v, k + 1))
                w = succ[v][k]
                if color[w] == 1:
                    back.append((v, w))
                elif color[w] == 0:
                    color[w] = 1
                    stack.append((w, 0))
            else:
                color[v] = 2
        loop_iter = None
        if back:
            dag = [[w for w in succ[v] if (v, w) not in set(back)] for v in range(n)]
            order, seen, stack = [], [False] * n, []
            for root in range(n):
                if seen[root]:
                    continue
                seen[root] = True
                stack.append((root, 0))
                while stack:
                    v, k = stack.pop()
                    if k < len(dag[v]):
                        stack.append((v, k + 1))
                        w = dag[v][k]
                        if not seen[w]:
                            seen[w] = True
                            stack.append((w, 0))
                    else:
                        order.append(v)
            loop_iter = 0
            for u, h in set(back):
                best = [0] * n   # longest path v -> u in instructions, 0 if u is unreachable
                best[u] = 1
                for v in order:   # post-order: successors first
                    if v == u:
                        continue
                    cand = [best[w] for w in dag[v] if best[w] > 0]
                    best[v] = 1 + max(cand) if cand else 0
                loop_iter = max(loop_iter, best[h])
        out[name] = dict(
            instructions=n, cfg_has_cycle=bool(back), loop_iteration_max=loop_iter,
            calls=calls, rets=rets,
            mul_umulh=mix.get("mul", 0) + mix.get("umulh", 0),
            loads_stores=sum(v for k, v in mix.items() if k.startswith(("ldr", "ldp", "ldu", "str", "stp", "stu"))),
        )
    return out


# prices of proof.md section 5.3 (instructions), checked against the counts
PRICE_STEP, PRICE_TRAVERSAL, PRICE_DOUBLING, PRICE_INV = 1350, 1450, 1800, 1 << 17


def kernels_check():
    import os
    from decimal import Decimal as D, getcontext
    getcontext().prec = 60
    here = os.path.dirname(os.path.abspath(__file__))
    res = {}
    res["python_port_vs_bigint"] = _kernel_port_check(2500)
    asm = {}
    for fname in ("p256_field.s", "p256_loops.s", "p256_slow.s"):
        for base in (os.path.join(here, "kernels"), "/snapshot/files/kernels", os.path.join(here, "files", "kernels")):
            path = os.path.join(base, fname)
            if os.path.exists(path):
                with open(path, encoding="utf-8") as fh:
                    asm.update(_asm_counts(fh.read()))
                break
        else:
            res["asm"] = fname + " not found"
            return res
    res["asm_per_function"] = {k: v for k, v in sorted(asm.items())}
    batch_share = 4 + 1   # I / 2^15 with I = 2^17; plus < 2^12 other instructions per batch of 2^15
    step = asm["pass1"]["loop_iteration_max"] + asm["pass2"]["loop_iteration_max"] + batch_share
    dbl = asm["pass1_dbl"]["loop_iteration_max"] + asm["pass2_dbl"]["loop_iteration_max"] + batch_share
    popc = bin(PROD["q"] - 2).count("1")
    # Fermat inversion: 255 loop iterations (a squaring and a conditional multiplication each), plus entry and exit
    inv_bound = 255 * asm["fe_inv"]["loop_iteration_max"] + asm["fe_inv"]["instructions"]
    res["popcount_p_minus_2"] = popc
    res["step_instructions_bound"] = step
    res["doubling_instructions_bound"] = dbl
    res["inversion_instructions_bound"] = inv_bound
    res["priced"] = dict(step=PRICE_STEP, traversal=PRICE_TRAVERSAL, doubling=PRICE_DOUBLING, inversion=PRICE_INV)
    handoff = 60   # flag test, call, return and the decision in the slow path (allowance)
    trav = step + asm["trav_update"]["instructions"] + handoff
    dbl_total = dbl + asm["move_to_dbl"]["instructions"] + asm["move_back"]["instructions"] + handoff
    res["traversal_instructions_bound"] = trav
    res["doubling_with_moves_instructions_bound"] = dbl_total
    res["price_covers_bound"] = (step <= PRICE_STEP and dbl_total <= PRICE_DOUBLING and inv_bound <= PRICE_INV
                                 and trav <= PRICE_TRAVERSAL)
    res["margin_step"] = round(PRICE_STEP / step, 4)
    res["margin_doubling"] = round(PRICE_DOUBLING / dbl_total, 4)
    res["margin_traversal"] = round(PRICE_TRAVERSAL / trav, 4)

    # ledger (proof.md sections 4 and 5), exact to 50 digits. Starts, "others" and the hard
    # extras (2^97.4) are inputs taken from the ledger of proof.md section 5.4.
    two = D(2)
    lg = lambda x: x.ln() / two.ln()
    r = PROD["r"]
    N = D((r - 1) // 2)
    Rr = D(2048)
    Np = N * Rr / (Rr - 1)
    S0 = (2 * Np * (D(1) / D("0.605")).ln()).sqrt()
    Sops = D("1.05") * S0
    step_u, dbl_u, trav_u = D(PRICE_STEP) / 1000, D(PRICE_DOUBLING) / 1000, D(PRICE_TRAVERSAL) / 1000
    ph1 = step_u * (Sops + two ** 64)
    hard = two ** D("97.4")
    rate = two ** D("-11.99")
    extras = (rate * Sops * (dbl_u - step_u) + (trav_u - step_u) * 64 * two ** D("-19.99") * Sops
              + two ** D("76.7") + two ** 70)
    claim = D("128.06")
    Wmax = two ** claim
    head = Wmax - ph1 - hard
    markov = extras / head
    res["ledger"] = dict(
        log2_S0=float(lg(S0)), log2_Sops=float(lg(Sops)), log2_phase1=float(lg(ph1)),
        log2_hard_extras=float(lg(hard)), log2_expected_random_extras=float(lg(extras)),
        log2_headroom=float(lg(head)), markov_early_halt_probability=float(markov),
        success_lower_bound=float(D("0.395") - markov - D("2e-7")),
        claimed_time_log2=float(claim),
        claim_above_phase1_plus_hard=bool(Wmax > ph1 + hard),
        success_at_least_0_39=bool(D("0.395") - markov - D("2e-7") > D("0.39")),
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
            esc2_expected_events=round(steps / (2 * R), 1),
            esc2_per_step=sum(x["esc2"] for x in res) / max(1, steps),
            esc2_model_1_over_2R=1 / (2 * R),
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
