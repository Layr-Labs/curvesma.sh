"""Six-class (negation + GLV) parallel Pollard rho on BLS12-381 G1.

Part A checks on BLS12-381 G1 itself every identity the attack uses: the GLV
pair (beta, lambda), the six-element orbits under negation and phi, the |y|
class key, the unique canonical member, the truncated-product twist index of
proof.md section 1.1 against the exact rule (random and adversarial near-tie
inputs), the twisted-table step as a class map, the doubling escape as a class
map, the order of 2 mod r, coefficient tracking along production steps, and
recovery of a planted x from two production trails forced to merge.

Part B runs the complete walk of proof.md (twisted R-adding class walk,
2-cycle escape, checkpoint escape, trail cap, distinguished classes, re-walk
solving, verification) on the prime-order subgroup of a toy curve
y^2 = x^3 + 4 over F_q and measures class visits to the first collision
between different trails, against the random-mapping model.

The "cost" block prints the whole cost model of proof.md sections 3-5 from
exact integers (r, N, R) and 80-digit decimals, including the exact-integer
check of the claimed score. Toy evidence only; cost bounds are analytic.
"""
import json, math, random, sys, time
from decimal import Decimal, getcontext

getcontext().prec = 80
T0 = time.time()
cfg = json.load(sys.stdin)
curve = cfg.get("curve", "q26")
trials = cfg.get("trials", 120)
R = cfg.get("R", 256)
dp_bits = cfg.get("dp_bits", 4)
seed = cfg.get("seed", 1)
prod_checks = cfg.get("production_checks", 1)
budget = cfg.get("seconds", 14)
dump = cfg.get("dump", False)
pipeline = cfg.get("pipeline", True)
run_toy = cfg.get("toy", True)
canon_random = cfg.get("canon_random", 2000)
canon_stress = cfg.get("canon_stress", 1000)
CLS = cfg.get("classes", 6)
rng = random.Random(seed)

PROD = dict(
    q=int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab", 16),
    r=int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16),
    h=int("396c8c005555e1568c00aaab0000aaab", 16),
    gx=int("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb", 16),
    gy=int("08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1", 16),
    beta=int("1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac", 16),
    lam=int("ac45a4010001a40200000000ffffffff", 16),
)
TOYS = {"q26": (58493047, 58489293, 3), "q30": (861101287, 861140127, 3)}


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
        self.fb = []
        Q = P
        for i in range((self.r.bit_length() + 7) // 8):
            row = [None]
            for j in range(1, 256):
                row.append(self.add(row[-1], Q))
            self.fb.append(row)
            Q = self.add(row[-1], Q)

    def fmul(self, k):
        """[k]G through the fixed-base table (toy harness only)."""
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
        G, P, k = None, (x, y), h
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


class Walker:
    """Exactly the walk of proof.md section 2, with the production rules."""

    def __init__(self, E, R, dp_bits, T, rng, start_point=None):
        self.E, self.R, self.dp_bits = E, R, dp_bits
        self.rbits = R.bit_length() - 1
        self.dpmask = (1 << dp_bits) - 1
        self.Lmax = 20 << dp_bits
        r, lam = E.r, E.lam
        self.units = [(-1) ** s * pow(lam, j, r) % r for s in (0, 1) for j in (0, 1, 2)]
        self.uinv = [pow(u, -1, r) for u in self.units]
        self.T = T
        self.table = []
        while len(self.table) < R:
            c, d = rng.randrange(r), rng.randrange(r)
            M = start_point(c, d) if start_point else E.add(E.mul(c, E.G), E.mul(d, T))
            if M is None:
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
        self.stats = dict(steps=0, doublings=0, esc2=0, esccp=0, abandoned=0, exceptional=0)

    def escape(self, P, a, b, track):
        r = self.E.r
        self.stats["doublings"] += 1
        Q = self.E.add(P, P)
        if track:
            a, b = 2 * a % r, 2 * b % r
        return Q, a, b

    def add_checked(self, P, Q):
        if P[0] == Q[0]:
            return None
        return self.E.add(P, Q)

    def canon(self, P):
        q, beta = self.E.q, self.E.beta
        x, y = P
        if CLS == 1:
            return x + q * y, 0, P
        if CLS == 2:
            if q - y < y:
                return x, 3, (x, q - y)
            return x, 0, P
        bx = beta * x % q
        b2x = (-x - bx) % q
        j, cx = 0, x
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
        E, r = self.E, self.E.r
        k1 = k2 = cp = None
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
                pending_double = True
            elif cp is not None and K == cp:
                self.stats["esccp"] += 1
                best = (K, P, a, b)
                Q, qa, qb = P, a, b
                trav = 0
                while True:
                    trav += 1
                    if trav > 64:
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
                        on_class(KQ)
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


# --------------------------------------------------- section 1.1 truncated product
M64 = (1 << 64) - 1
GQ = (PROD["beta"] << 448) // PROD["q"]
GL = [(GQ >> (64 * j)) & M64 for j in range(7)]
PH0, PH1 = (PROD["q"] >> 256) & M64, PROD["q"] >> 320
PTOP = PROD["q"] >> 320
TOL = 16


def exact_j(x):
    q, b = PROD["q"], PROD["beta"]
    v = [x, b * x % q, b * b % q * x % q]
    return min(range(3), key=lambda j: v[j])


def fast_j(x):
    xl = [(x >> (64 * i)) & M64 for i in range(6)]
    L = 0
    for i in range(6):
        L += (xl[i] * GL[5 - i]) >> 64
        L += (xl[i] * GL[6 - i]) & M64
    L &= M64
    if L < 16 or L > M64 - 16:
        return None
    V1 = ((PH0 * L >> 64) + PH1 * L) >> 64
    X = x >> 320
    S = X + V1
    if abs(S - PTOP) <= TOL:
        return None
    V2 = PTOP - S if S < PTOP else 2 * PTOP - S
    vals = [X, V1, V2]
    jm = min(range(3), key=lambda j: vals[j])
    for j in range(3):
        if j != jm and vals[j] - vals[jm] <= TOL:
            return None
    return jm


def check_fast_j(n_random, n_stress):
    q, b = PROD["q"], PROD["beta"]
    rng2 = random.Random(seed + 10 ** 7)
    res = dict(random=n_random, random_fallbacks=0, random_mismatches=0,
               stress=0, stress_fallbacks=0, stress_mismatches=0)
    for _ in range(n_random):
        if time.time() - T0 > budget:
            break
        x = rng2.randrange(1, q)
        f = fast_j(x)
        if f is None:
            res["random_fallbacks"] += 1
        elif f != exact_j(x):
            res["random_mismatches"] += 1
    mults = [b - 1, b * b - 1, b - b * b, 1 + b, b]
    for _ in range(n_stress):
        if time.time() - T0 > budget:
            break
        k = rng2.randrange(300, 381)
        d = rng2.randrange(-(1 << k), 1 << k)
        mlt = mults[rng2.randrange(len(mults))]
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


# ------------------------------------------------------------------ Part A
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
    out["log2_r"] = round(math.log2(r), 10)
    out["log2_N_classes"] = round(math.log2((r - 1) // 6), 10)
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
    out["r_minus_1_fully_factored"] = rm1 == 1
    out["log2_order_of_2_mod_r"] = round(math.log2(order), 4)
    xs = rng.randrange(1, r)
    T = E.mul(xs, G)
    W = Walker(E, 4, 0, T, rng)
    ok_orbit = ok_key = ok_canon = ok_twist = ok_double = True
    for _ in range(nchecks * 2):
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
            if W.canon(E.add(Mm, TP))[2] != W.canon(E.add(c, Mh))[2]:
                ok_twist = False
            if W.canon(E.add(Mm, Mm))[2] != W.canon(E.add(c, c))[2]:
                ok_double = False
    out["orbit_size_6_and_unit_action"] = ok_orbit
    out["key_abs_y_is_class_invariant"] = ok_key
    out["canonical_member_unique"] = ok_canon
    out["twisted_step_is_class_map"] = ok_twist
    out["doubling_is_class_map"] = ok_double
    out["twisted_table_entries_equal_u_inverse_M"] = all(
        W.table[h][t][0] == E.mul(W.uinv[t], W.table[h][0][0]) for h in range(W.R) for t in range(6))
    # production walk with coefficient tracking, forced merge, solve
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
    v = rng.randrange(6)
    u = Wp.units[v]
    bB = rng.randrange(r)
    aB = (u * (aA + bA * xs) - bB * xs) % r
    PB = E.mul(u, PA)
    out["second_trail_start_is_same_class"] = E.add(E.mul(aB, G), E.mul(bB, T)) == PB
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
    out["fast_canonical_index"] = check_fast_j(canon_random, canon_stress)
    return out


# ------------------------------------------------------------------ Part B
def toy_trial(E, rng, R, dp_bits):
    r = E.r
    xs = rng.randrange(1, r)
    T = E.fmul(xs)

    def point(a, b):
        return E.fmul(a + b * xs)

    W = Walker(E, R, dp_bits, T, rng, point)
    N = (r - 1) // CLS
    first_seen = {}
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
    B = 1.05 * math.sqrt(2 * Np * math.log(1 / 0.6095))
    D_cap = int(2 * (B + W.Lmax) / (1 << dp_bits)) + 1
    matches = 0
    while solved is None:
        if pipeline and state["visits"] >= B:
            solved = False
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
            W.stats.update(s1)
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


# ------------------------------------------------------------------ cost model
def cost():
    r = PROD["r"]
    N = (r - 1) // 6
    R = 1 << 16
    Np = Decimal(N) * Decimal(R) / Decimal(R - 1)
    l2 = lambda x: Decimal(x).ln() / Decimal(2).ln()
    L, C = 128, 2000
    p_step = Decimal(1) + Decimal(L) / C
    p_trav = Decimal(1) + Decimal(L + 24) / C
    p_dbl = Decimal("1.5")
    kappa, floor, ph = Decimal("1.05"), Decimal("0.3905"), Decimal("0.00048")
    S0 = (2 * Np * (Decimal(1) / (1 - floor)).ln()).sqrt()
    Sops = kappa * S0
    S60 = Sops + Decimal(2) ** 60
    U1 = p_step * S60
    rate2, short = Decimal(1) / (6 * R), Decimal(4) / (R * R)
    E_dbl = (p_dbl - p_step) * (rate2 + short) * S60
    E_trav = (p_trav - p_step) * 64 * short * 16 * S60
    E_rand = Decimal(2) ** Decimal("74.8")
    E_tot = E_dbl + E_trav + E_rand
    head = E_tot / ph
    hard = sum([Decimal(2) ** 61, Decimal(2) ** Decimal("95.4"), Decimal(2) ** Decimal("79.4"),
                Decimal(2) ** Decimal("62.3"), Decimal(2) ** Decimal("57.1"), Decimal(2) ** 28,
                Decimal(2) ** 13])
    W = U1 + hard + head
    Wup = int(W.to_integral_value(rounding="ROUND_CEILING")) + 2
    k = math.ceil(math.log2(Wup) * 1000)
    head_at = Decimal(2) ** (Decimal(k) / 1000) - U1 - hard
    ph_bound = E_tot / head_at
    Dcap = Decimal(2) ** Decimal("63.25")
    mean_rep = (S60 + Decimal(2) ** Decimal("88.32")) / Decimal(2) ** 64
    mem = Dcap * 104 * 2 + 6 * R * 160 + Decimal(2) ** 20 * 512
    pre = Decimal(2) ** Decimal("57.1") + Decimal(2) ** 28 + Decimal(2) ** 32 + Decimal(2) ** 21
    f = lambda x, n=10: round(float(l2(x)), n)
    return {
        "C": C, "R": R, "table_rows_x_units": 6 * R,
        "r": str(r), "N": str(N), "log2_N": f(N), "log2_N_prime": f(Np),
        "kappa": "1.05", "H1_floor": "0.3905",
        "S0_log2": f(S0), "S_ops_log2": f(Sops), "S_ops_int": str(int(Sops)),
        "step_instructions": L, "step_price": str(p_step),
        "traversal_price": str(p_trav), "doubling_price": str(p_dbl),
        "U1_log2": f(U1),
        "doubling_rate_log2": round(float(l2(rate2)), 6),
        "short_cycle_rate_log2": round(float(l2(short)), 6),
        "E_dbl_log2": round(float(l2(E_dbl)), 6),
        "E_trav_log2": round(float(l2(E_trav)), 6),
        "E_rand_log2": 74.8,
        "E_total_log2": round(float(l2(E_tot)), 6),
        "p_halt_target": "0.00048", "headroom_log2": round(float(l2(head)), 6),
        "headroom_at_W_max_log2": round(float(l2(head_at)), 6),
        "p_halt_log2": round(float(l2(ph_bound)), 4),
        "p_halt_bound": round(float(ph_bound), 9),
        "hard_extras_log2": round(float(l2(hard)), 6),
        "W_log2": f(W), "W_up": str(Wup), "claim": str(k) + "/1000",
        "check_2k_ge_W_up_1000": pow(2, k) >= pow(Wup, 1000),
        "check_2k_minus_1_lt_W_up_1000": pow(2, k - 1) < pow(Wup, 1000),
        "success_bound": str(round(0.3905 - round(float(ph_bound), 9) - 0.00000012, 7)),
        "report_capacity_log2": "63.25",
        "mean_reports_log2": round(float(l2(mean_rep)), 4),
        "report_capacity_over_mean_log2": round(float(l2(Dcap / mean_rep)), 4),
        "memory_log2_bytes": round(float(l2(mem)), 4),
        "table_bytes": 6 * R * 160,
        "preprocessing_expected_lt": "57.11",
        "preprocessing_cap_log2": "57.2",
        "phase2_log2": round(float(l2(Decimal(2) ** 20 * (20 * Decimal(2) ** 64 + 65) * 133)), 4),
        "rewalk_log2": round(float(l2(8 * 2 * ((20 * Decimal(2) ** 64 + 65) * 133 + Decimal(2) ** 12))), 4),
    }


def main():
    out = {"config": cfg}
    if prod_checks:
        out["part_a"] = part_a(prod_checks)
    out["cost"] = cost()
    if not run_toy:
        out["seconds"] = round(time.time() - T0, 2)
        print(json.dumps(out, indent=1))
        return
    E = make_toy(curve)
    r = E.r
    N = (r - 1) // CLS
    Np = N * R / (R - 1)
    S0 = math.sqrt(2 * Np * math.log(1 / 0.6095))
    S1 = 1.05 * S0
    out["toy"] = dict(q=E.q, r=r, h=E.h, beta=E.beta, lam=E.lam, N=N, R=R, dp_bits=dp_bits,
                      L_max=20 << dp_bits,
                      mean_trail_over_sqrtN=round((1 << dp_bits) / math.sqrt(N), 4),
                      S0=round(S0, 1), S_ops=round(S1, 1),
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
        cens = pipeline
        out["toy"]["summary"] = dict(
            trials_completed=n,
            mode="pipeline (scaled phase-1 budget, phase 2, report cap, 8-match cap)" if pipeline else "run until solved",
            pipeline_success_rate=round(sum(1 for x in res if x["ok"]) / n, 4) if pipeline else None,
            all_solved_correctly=None if pipeline else all(x["ok"] for x in res),
            mean_tc_over_sqrtN=None if cens else round(sum(tcs) / len(tcs) / math.sqrt(N), 4),
            model_mean_over_sqrtN=round(math.sqrt(math.pi * Np / 2) / math.sqrt(N), 4),
            frac_tc_le_S_ops=round(sum(t <= S1 for t in tcs) / n, 4),
            binomial_stderr_at_0p3905=round(math.sqrt(0.3905 * 0.6095 / n), 4),
            model_P_S_ops=0.3905,
            total_steps=steps,
            esc2_events=sum(x["esc2"] for x in res),
            esc2_expected_events=round(steps / (6 * R), 1),
            esc2_per_step=round(sum(x["esc2"] for x in res) / max(1, steps), 6),
            esc2_model_1_over_6R=round(1 / (6 * R), 6),
            checkpoint_escapes=sum(x["esccp"] for x in res),
            abandoned_trails=sum(x["abandoned"] for x in res),
            exceptional_additions=sum(x["exceptional"] for x in res),
        )
        if dump:
            out["tc_list"] = tcs
    out["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(out, indent=1))


main()
