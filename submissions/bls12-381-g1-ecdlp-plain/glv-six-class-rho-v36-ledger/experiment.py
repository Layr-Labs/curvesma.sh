import json, math, random, sys, time

T0 = time.time()
cfg = json.load(sys.stdin)
curve = cfg.get("curve", "q26")
trials = cfg.get("trials", 100)
R = cfg.get("R", 256)
dp_bits = cfg.get("dp_bits", 4)
seed = cfg.get("seed", 1)
pchk = cfg.get("production_checks", 2)
budget = cfg.get("seconds", 14)
dump = cfg.get("dump", False)
pipeline = cfg.get("pipeline", False)
KAPPA = 1.005
H1F = 0.3903
c_rnd = cfg.get("c_rnd", None)
c_str = cfg.get("c_str", None)
rtoy = cfg.get("toy", True)
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
TOYS = {
    "q26": (58493047, 58489293, 3),
    "q30": (861101287, 861140127, 3),
    "q34": (10146855109, 10146690129, 3),
}


class Curve:

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
        G = None
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


class Walker:

    def __init__(self, E, R, dp_bits, T, rng, start_point=None):
        self.E, self.R, self.dp_bits = E, R, dp_bits
        self.rbits = R.bit_length() - 1
        assert 1 << self.rbits == R
        self.shift = 1
        self.dpmask = (1 << dp_bits) - 1
        self.Lmax = 20 << dp_bits
        r = E.r
        lam = E.lam
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
        bx = beta * x % q
        b2x = (-x - bx) % q
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


    def step_target(self, P, K, t):
        h = (K >> self.shift) & (self.R - 1)
        return self.table[h][t]

    def is_dp(self, K):
        return ((K >> (self.shift + self.rbits)) & self.dpmask) == 0

    def run_trail(self, P, coeffs=None, on_class=None):
        E, r = self.E, self.E.r
        cp = None
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


M64 = (1 << 64) - 1
GQ = (PROD["beta"] << 448) // PROD["q"]
GL = [(GQ >> (64 * j)) & M64 for j in range(7)]
PH0, PH1 = (PROD["q"] >> 256) & M64, PROD["q"] >> 320
PTOP = PROD["q"] >> 320
TOL = 16


def xj(x):
    q, b = PROD["q"], PROD["beta"]
    v = [x, b * x % q, b * b % q * x % q]
    return min(range(3), key=lambda j: v[j])


def fj(x):
    xl = [(x >> (64 * i)) & M64 for i in range(6)]
    L = 0
    for i in range(6):
        L += (xl[i] * GL[5 - i]) >> 64
        L += (xl[i] * GL[6 - i]) & M64
    L &= M64
    if L < 16 or L > M64 - 15:
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


def chk_row45(n_random, n_boundary):
    q = PROD["q"]
    pl = [(q >> (64 * i)) & M64 for i in range(4)]
    rng = random.Random(seed + 2 * 10 ** 7)
    res = dict(random=0, boundary=0, dp_keys=0, mismatches=0, k_b0=0,
               neg_var=0, addr_n=0, addr_mm=0)
    cases = [rng.randrange(1, q) for _ in range(n_random)]
    half = q >> 1
    cases += [1, 2, 3, q - 1, q - 2, q - 3] + [half + d for d in range(-3, 4)]
    picks = [0, 1, M64, M64 - 1] + pl + [v + 1 for v in pl] + [v - 1 for v in pl]
    while len(cases) < n_random + 13 + n_boundary:
        y = rng.randrange(1, q) >> 256 << 256
        for i in range(4):
            y |= (picks[rng.randrange(len(picks))] if rng.random() < 0.7 else rng.getrandbits(64)) << (64 * i)
        if 0 < y < q:
            cases.append(y)
    n_dp = max(1, n_boundary // 4)
    while len(cases) < n_random + 13 + n_boundary + n_dp:
        K = (rng.randrange(2, q) >> 128 << 128) | (rng.getrandbits(64) & ~1)
        if 0 < K < q:
            cases.append(K if rng.random() < 0.5 else q - K)
    table = (rng.getrandbits(32) << 31) & M64
    for i, y in enumerate(cases):
        res["random" if i < n_random else "boundary" if i < n_random + 13 + n_boundary
            else "dp_keys"] += 1
        yl = [(y >> (64 * k)) & M64 for k in range(3)]
        d, borrow = [], 0
        for k in range(3):
            v = pl[k] - yl[k] - borrow
            borrow = 1 if v < 0 else 0
            d.append(v & M64)
        s_bit = yl[0] & 1
        K3 = [d[k] if s_bit else yl[k] for k in range(3)]
        h = (K3[0] >> 1) & 0x3FFFFF
        dist = K3[1] == 0
        Kx = y if y % 2 == 0 else q - y
        if Kx & 1:
            res["k_b0"] += 1
        yn = q - y
        if (yn if yn % 2 == 0 else q - yn) != Kx:
            res["neg_var"] += 1
        fp = K3[0] | K3[1] << 64 | K3[2] << 128
        if (s_bit != (y & 1) or K3 != [(Kx >> (64 * k)) & M64 for k in range(3)]
                or fp != Kx % (1 << 192) or h != (Kx >> 1) % (1 << 22)
                or dist != ((Kx >> 64) % (1 << 64) == 0)):
            res["mismatches"] += 1
        for j in range(3):
            x5 = (j + (s_bit << 2)) & M64
            x5 = (x5 + (h << 3)) & M64
            x5 = (table + (x5 << 8)) & M64
            t = 3 * (y % 2) + j
            want = table + 256 * (8 * ((Kx >> 1) % (1 << 22)) + 4 * (t // 3) + t % 3)
            res["addr_n"] += 1
            if x5 != want:
                res["addr_mm"] += 1
    idx = set()
    for hh in (0, 1, 2, 12345, (1 << 22) - 2, (1 << 22) - 1):
        for t in range(6):
            idx.add(8 * hh + 4 * (t // 3) + t % 3)
    res["addr_injective"] = len(idx) == 36 and max(idx) < 1 << 25
    return res


def _flags_sub(a, b):
    res = (a - b) & M64
    return (res >> 63, int(res == 0), int(a >= b), (((a ^ b) & (a ^ res)) >> 63) & 1)


def _flags_add(a, b):
    tot = a + b
    res = tot & M64
    return (res >> 63, int(res == 0), tot >> 64, ((~(a ^ b) & (a ^ res)) >> 63) & 1)


def _nzcv(k):
    return ((k >> 3) & 1, (k >> 2) & 1, (k >> 1) & 1, k & 1)


_COND = dict(eq=lambda f: f[1] == 1, hs=lambda f: f[2] == 1, lo=lambda f: f[2] == 0,
             hi=lambda f: f[2] == 1 and f[1] == 0, ls=lambda f: f[2] == 0 or f[1] == 1)


def listing_rows_1_3(x):
    xl = [(x >> (64 * i)) & M64 for i in range(6)]
    x9 = 0
    for i in range(6):
        x9 = (x9 + ((xl[i] * GL[5 - i]) >> 64) + ((xl[i] * GL[6 - i]) & M64)) & M64
    x5 = xl[5]
    x8 = (x9 - 16) & M64
    x0 = (PH0 * x9) >> 64
    x1 = (PH1 * x9) & M64
    x2 = (PH1 * x9) >> 64
    f = _flags_add(x1, x0)
    x1 = (x2 + f[2]) & M64
    x2 = (x5 + x1) & M64
    f = _flags_sub(x2, PTOP)
    x3 = (x2 - PTOP) & M64
    x4 = PTOP if _COND["hs"](f) else 0
    x13 = (x3 + 16) & M64
    x2 = (x4 - x3) & M64
    f = _flags_sub(x1, x5)
    lo = _COND["lo"](f)
    x3, x4, x6 = (x1 if lo else x5), (x5 if lo else x1), (1 if lo else 0)
    f = _flags_sub(x2, x4)
    x7 = x2 if _COND["lo"](f) else x4
    f = _flags_sub(x2, x3)
    x7 = x3 if _COND["lo"](f) else x7
    x3 = x2 if _COND["lo"](f) else x3
    x6 = 2 if _COND["lo"](f) else x6
    x7 = (x7 - x3) & M64
    f = _flags_sub(x13, 32)
    f = _flags_add(x8, 31) if _COND["hi"](f) else _nzcv(2)
    f = _flags_sub(x7, 16) if _COND["lo"](f) else _nzcv(0)
    fallback = _COND["ls"](f)
    L = 0
    for i in range(6):
        L += (xl[i] * GL[5 - i]) >> 64
        L += (xl[i] * GL[6 - i]) & M64
    L &= M64
    V1 = ((PH0 * L >> 64) + PH1 * L) >> 64
    X = x >> 320
    S = X + V1
    V2 = PTOP - S if S < PTOP else 2 * PTOP - S
    v = sorted([X, V1, V2])
    tests = (L < 16 or L > M64 - 15, abs(S - PTOP) <= TOL, v[1] - v[0] <= TOL)
    return fallback, x6, tests


def _listing_compare(x, f, res, kind):
    fb, j, tests = listing_rows_1_3(x)
    if fb != any(tests) or fb != (f is None):
        res["listing_fb_mm"] += 1
    if not fb and j != f:
        res["listing_j_mm"] += 1
    if kind == "stress":
        for name, hit in zip(("A_range", "B_wrap", "C_gap"), tests):
            res["stress_test_" + name] += hit


def check_fj(n_random, n_stress):
    q, b = PROD["q"], PROD["beta"]
    rng = random.Random(seed + 10 ** 7)
    res = dict(random=n_random, random_fallbacks=0, random_mismatches=0,
               stress=0, stress_fallbacks=0, stress_mismatches=0,
               listing_merged_fallback_mismatches=0, listing_j_mismatches=0,
               stress_test_A_range=0, stress_test_B_wrap=0, stress_test_C_gap=0)
    for _ in range(n_random):
        if time.time() - T0 > budget:
            break
        res["rand_done"] = res.get("rand_done", 0) + 1
        x = rng.randrange(1, q)
        f = fj(x)
        if f is None:
            res["random_fallbacks"] += 1
        elif f != xj(x):
            res["rand_mm"] += 1
        _listing_compare(x, f, res, "random")
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
        f = fj(x)
        if f is None:
            res["stress_fallbacks"] += 1
        elif f != xj(x):
            res["stress_mm"] += 1
        _listing_compare(x, f, res, "stress")
    return res


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
    out["points_x0_are_3_torsion_not_G1"] = E.add(E.add((0, 2), (0, 2)), (0, 2)) is None
    out["log2_r"] = math.log2(r)
    out["log2_N_classes"] = math.log2((r - 1) // 6)
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
    W = Walker(E, 4, 0, T, rng)
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
        Kp = P[1] if P[1] % 2 == 0 else q - P[1]
        absy = min(P[1], q - P[1])
        if (len({c[0] for c in cans}) != 1 or cans[0][0] != Kp or Kp & 1
                or Kp != (absy if absy % 2 == 0 else q - absy)
):
            ok_key = False
        if len({c[2] for c in cans}) != 1:
            ok_canon = False
        for Mm, (K, t, c) in zip(members, cans):
            if E.mul(W.units[t], Mm) != c:
                ok_canon = False
            h = (K >> W.shift) & (W.R - 1)
            TP = W.table[h][t][0]
            Mh = W.table[h][0][0]
            n1 = W.canon(E.add(Mm, TP))
            n2 = W.canon(E.add(c, Mh))
            if n1[2] != n2[2]:
                ok_twist = False
            if W.canon(E.add(Mm, Mm))[2] != W.canon(E.add(c, c))[2]:
                ok_double = False
    out["orbit_size_6_and_unit_action"] = ok_orbit
    out["key_parity_is_class_invariant_bit0_zero"] = ok_key
    out["canonical_member_unique"] = ok_canon
    out["twisted_step_is_class_map"] = ok_twist
    out["doubling_is_class_map"] = ok_double
    ok_tab = all(W.table[h][t][0] == E.mul(W.uinv[t], W.table[h][0][0]) for h in range(W.R) for t in range(6))
    out["twisted_table_entries_equal_u_inverse_M"] = ok_tab

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
    nr = 15000 * nchecks if c_rnd is None else c_rnd
    ns = 5000 * nchecks if c_str is None else c_str
    out["fastcan"] = check_fj(nr, ns)
    out["row45"] = chk_row45(20000 * nchecks, 2000 * nchecks)
    out["part_a_seconds"] = round(time.time() - T0, 2)
    return out


def toy_trial(E, rng, R, dp_bits):
    r = E.r
    xs = rng.randrange(1, r)
    T = E.fmul(xs)

    def point(a, b):
        return E.fmul(a + b * xs)

    W = Walker(E, R, dp_bits, T, rng, point)
    N = (r - 1) // 6
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
    B = KAPPA * math.sqrt(2 * Np * math.log(1 / 0.6097))
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





from fractions import Fraction as Fr

CLAIM = (1262053, 10000)
LEDGER = [("1.1 v1 (first load pre-indexed)", 24), ("1.1 L' - 16, V1, X, S, S - PTOP + 16, V2", 11),
          ("1.1 decision and the one merged fallback branch", 15),
          ("parity bit of y (ands), p - y limbs 0-2, select limbs 0-2 of K", 9),
          ("h, distinguished test, table address 8h + 4s + j", 6),
          ("loads of table point and state, store of state", 19),
          ("checkpoint test (limbs 0-2, validity via limb 1), refresh test", 8),
          ("length increment (69-bit), loop (cap test in refresh)", 6)]
SLACK = 0
CC, MM, AS = 2000, 383, 32
FINV = (380 + bin(PROD["q"] - 2).count("1") - 1) * MM
B_MN = 1 << 22
B_IT = [("zero-denominator test on the batch product", 9), ("batch meter charge", 10),
               ("Fermat inversion", FINV), ("exponent scan (< 2,000)", 1999),
               ("per-batch constants and pointers (x10-x12, x16, x17, x29, x19-x23)", 16)]
B_SH = Fr(561, 10000)
OVS = 44
RF_L, RF_I, RF_P = 11, 11, 64
ORND, HRD = 76.5, 97.3
TRV_X = 25
TRV_S = Fr(135, 10000)
RECORD = 504
R_F = [("x", 0, 48), ("y", 48, 48), ("next table-entry address", 96, 8),
                 ("K_cp limbs 0, 1", 104, 16), ("K_cp limb 2", 120, 8), ("length (2 limbs)", 128, 16),
                 ("start coefficient a", 144, 32), ("start coefficient b", 176, 32),
                 ("remembered member x, y (traversal)", 208, 96), ("its fingerprint", 304, 24),
                 ("traversal counter", 328, 8), ("batch denominator", 336, 48),
                 ("batch prefix product", 384, 48)]
L_OFF = [("ldp", 504), ("ldp", 16), ("ldp", 32), ("ldp", 48), ("ldr", 64), ("str", 96),
                   ("ldp", 104), ("ldr", 120), ("ldp", 128), ("stp", 128), ("stp", 104), ("str", 120)]
T_N, T_SL, T_BY = (1 << 22) * 8, 256, 160


def up(v, d):
    return (math.floor(v * 10 ** d) + 1) / 10 ** d


def dn(v, d):
    return (math.ceil(v * 10 ** d) - 1) / 10 ** d


def budgets():
    listed = sum(n for _, n in LEDGER)
    budget = Fr(listed) + B_SH + SLACK
    step = 1 + budget / CC
    refresh = Fr(RF_I, CC)
    per_op = step + refresh / RF_P
    return listed, budget, step, refresh, per_op


def cost_chain(c):
    listed, budget, step, refresh, per_op = budgets()
    dbl_instr = 7 * MM + 7 * AS + budget
    dbl = dbl_instr / CC
    Rt = 1 << 22
    Np = (PROD["r"] - 1) // 6 * Rt / (Rt - 1)
    S0 = math.sqrt(2 * Np * math.log(1 / (1 - H1F)))
    S_ops = KAPPA * S0
    lU1 = up(math.log2(float(per_op) * (S_ops + 2 ** OVS)), 6)
    l_rate = up(math.log2(1 / (6 * Rt) + 2 ** -38), 3)
    l_S = up(math.log2(S_ops), 4)
    l_dbl = up(math.log2(float(dbl - step)) + l_rate + l_S, 2)
    l_travn = up(math.log2((2 / (6 * Rt) + 64 * 2 ** -38.99) * (S_ops + 2 ** OVS)), 2)
    l_trav = up(math.log2(float(TRV_S)) + l_travn, 2)
    l_sum = math.log2(2 ** l_dbl + 2 ** l_trav + 2 ** ORND)
    l_head = dn(lU1 + math.log2(math.expm1((c - lU1) * math.log(2)) - 2 ** (HRD - lU1)), 2)
    l_meter = up(l_sum - l_head, 2)
    meter = up(2 ** l_meter, 6)
    meter3 = up((3 * 2 ** l_dbl + 3 * 2 ** l_trav + 2 ** ORND) / 2 ** l_head, 5)
    return dict(listed=listed, budget=budget, step=step, refresh=refresh, per_op=per_op,
                dbl_instr=dbl_instr, dbl=dbl, S0=S0, S_ops=S_ops, lU1=lU1,
                l_rate=l_rate, l_S=l_S, l_dbl=l_dbl, l_travn=l_travn, l_trav=l_trav,
                l_tot=up(l_sum, 2), l_head=l_head, l_meter=l_meter, meter=meter,
                success=dn(H1F - meter - 0.00000011, 5), meter3=meter3,
                success3=dn(H1F - meter3 - 0.00000011, 5))


def dec(f, d):
    assert (f * 10 ** d).denominator == 1
    return ("%." + str(d) + "f") % float(f)


def cost():
    c = CLAIM[0] / CLAIM[1]
    k = cost_chain(c)
    below = cost_chain((CLAIM[0] - 1) / CLAIM[1])
    trav_max = k["budget"] + TRV_X
    trav_price = k["step"] + TRV_S
    batch_work = sum(n for _, n in B_IT)
    small = FINV + 2000
    unbatched = FINV + 2000 + 3 * MM + 6 * AS + k["budget"] + RF_I + 100
    Lmax = 20 * 2 ** 64
    hard = [2 ** 63, 2 ** 82.3, 2 ** 22 * (Lmax + 65) * 119, 8 * 2 * ((Lmax + 65) * 119 + 2 ** 12),
            2 ** 63.25 * 0.5, 2 ** 57.1, 2 ** 34.3]
    rand = [2 ** 64.4 * (2 ** 12 + float(k["refresh"])), 2 ** 63, 2 ** 8, 2 ** 72.3, 2 ** -92]
    fields = sorted(R_F, key=lambda f: f[1])
    layout_ok = (all(o % 8 == 0 and n % 8 == 0 for _, o, n in fields)
                 and all(fields[i][1] + fields[i][2] <= fields[i + 1][1] for i in range(len(fields) - 1))
                 and fields[-1][1] + fields[-1][2] <= RECORD
                 and all(o % 8 == 0 and (o <= 504 if op in ("ldp", "stp") else o <= 32760)
                         for op, o in L_OFF) and RECORD == 504)
    ok = lambda kk: kk["success"] >= 0.39 + 1e-4 and kk["success3"] >= 0.39 + 1e-4
    return {
        "rows": [[name, n] for name, n in LEDGER],
        "listed": k["listed"], "batch_items": [[name, n] for name, n in B_IT],
        "batch_work": batch_work, "batch_share": "< " + dec(B_SH, 4),
        "step_budget": dec(k["budget"], 4),
        "price_step": dec(k["step"], 8),
        "price_rf": dec(k["refresh"], 4), "rf_listed": RF_L,
        "per_op": dec(k["per_op"], 10),
        "price_trv": dec(trav_price, 8), "trav_max": dec(trav_max, 4),
        "price_dbl": dec(k["dbl"], 8), "dbl_instr": dec(k["dbl_instr"], 4),
        "dbl_surch": dec(k["dbl"] - k["step"], 4),
        "finv": FINV, "small_inv": small,
        "unbatch": dec(unbatched, 4),
        "layout": [[name, o, n] for name, o, n in fields], "rec_used": fields[-1][1] + fields[-1][2],
        "tab_log2": math.log2(T_N * T_SL),
        "hard_sum": "%.4f" % math.log2(sum(hard)),
        "rand_sum": "%.4f" % math.log2(sum(rand)),
        "l2_S0": "%.5f" % math.log2(k["S0"]), "l2_Sops": "%.5f" % math.log2(k["S_ops"]),
        "l2_U1": "%.6f" % k["lU1"], "l_rate": "%.3f" % k["l_rate"],
        "l_dbl": "%.2f" % k["l_dbl"],
        "l_travn": "%.2f" % k["l_travn"], "l_trav": "%.2f" % k["l_trav"],
        "l_ornd": ORND, "l_tot": "%.2f" % k["l_tot"],
        "l_hrd": HRD,
        "claim": "%.4f" % c, "claim_fr": list(CLAIM), "l_head": "%.2f" % k["l_head"],
        "l_meter": "%.2f" % k["l_meter"], "meter": "%.6f" % k["meter"],
        "success": "%.5f" % k["success"],
        "meter_3x": "%.5f" % k["meter3"], "success_3x": "%.5f" % k["success3"],
        "checks": {
            "listed_98": k["listed"] == 98 and k["budget"] == Fr(980561, 10000),
            "batch_share": batch_work < B_SH * B_MN,
            "refresh_ok": RF_L <= RF_I,
            "trav_ok": trav_max + 2 <= (trav_price - 1) * CC,
            "dbl_ok": k["dbl_instr"] == k["dbl"] * CC,
            "small_118": small <= 118 * CC,
            "unbatch_119": unbatched <= 119 * CC,
            "finv_18": FINV + 2000 < 1 << 18,
            "layout_ok": layout_ok,
            "tab_ok": T_N * T_SL == 1 << 33 and T_BY <= T_SL,
            "hard_ok": math.log2(sum(hard)) < HRD,
            "rand_ok": math.log2(sum(rand)) < ORND,
            "succ_ok": k["success"] >= 0.39 + 1e-4,
            "succ3_ok": k["success3"] >= 0.39 + 1e-4,
            "claim_min": ok(k) and not ok(below),
            "claim_hard": c > math.log2(2 ** k["lU1"] + 2 ** HRD),
            "ovs_ok": OVS == 44,
        },
    }


def main():
    out = {"config": cfg}
    out["cost"] = cost()
    if pchk:
        out["part_a"] = part_a(pchk)
    if not rtoy:
        out["seconds"] = round(time.time() - T0, 2)
        print(json.dumps(out, indent=1))
        return
    E = make_toy(curve)
    r = E.r
    N = (r - 1) // 6
    Np = N * R / (R - 1)
    S0 = math.sqrt(2 * Np * math.log(1 / 0.6097))
    S1 = KAPPA * S0
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
        cens = pipeline
        m = sum(tcs) / len(tcs) / math.sqrt(N)
        summ = dict(
            trials_done=n, mode="pipeline" if pipeline else "tc",
            pipe_success=round(sum(1 for x in res if x["ok"]) / n, 4) if pipeline else None,
            all_solved=None if pipeline else all(x["ok"] for x in res),
            tc_over_sqrtN=None if cens else round(m, 4),
            frac_le_S1=round(sum(t <= S1 for t in tcs) / n, 4),
            binom_se_0p39=round(math.sqrt(0.39 * 0.61 / n), 4),
            H1_P_S1=0.3903, model_P_S1=round(1 - math.exp(-S1 * S1 / (2 * Np)), 4),
            steps_tot=steps, dbl_per_step=sum(x["doublings"] for x in res) / max(1, steps),
            esc2_exp=round(steps / (6 * R), 1), esc2_model=1 / (6 * R),
            esccp_rate=sum(x["esccp"] for x in res) / max(1, steps),
            abandoned=sum(x["abandoned"] for x in res), exceptional=sum(x["exceptional"] for x in res),
        )
        out["toy_summary"] = summ
        if dump:
            out["tc_list"] = tcs
    out["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(out, indent=1))


main()
