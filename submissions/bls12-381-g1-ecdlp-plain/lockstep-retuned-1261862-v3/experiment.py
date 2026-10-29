"""Lockstep six-class (negation + GLV) Pollard rho on BLS12-381 G1.

Part A checks, on BLS12-381 G1 itself, every algebraic identity and every listed
instruction row that the cost argument of proof.md relies on:
the GLV pair (beta, lambda), size-6 orbits and the unit action, the parity class
key, the fraction-domain canonical-index rule of proof.md section 1.1 against the
exact rule (random and adversarial near-tie inputs), rows 1-3 of the section-5
listing executed instruction by instruction on 64-bit registers with NZCV flags
(fallback branch against the OR of the three section-1.1 tests, j against the
rule), rows 4-5 and the refresh-batch row 4' limb by limb (key limbs, parity,
h = bits 1..40, distinguished test, table address, checkpoint compare and the
invalid marker), the twisted-table step identity, coefficient tracking along
production steps, the doubling escape as a class map and the collision-solving
formula on two production trails made to merge.

Part B runs the complete algorithm of proof.md on the prime-order subgroup of a
toy curve y^2 = x^3 + 4 over F_q: M walks in lockstep batches, global
checkpoint-test batches every k batches and refresh batches every 64k batches,
synchronous traversal and doubling escape, trail cap checked at refresh
batches, distinguished classes, restarts, phase-1 budget, phase 2, report
capacity, 8-match cap, replay re-walks from (a, b, start batch) and
verification. It reports T_c (class arrivals until the first inter-trail
collision) and the pipeline success rate against the random-mapping model.
Toy evidence only: it does not validate production cost bounds.
"""
import json, math, random, sys, time
from fractions import Fraction as Fr

T0 = time.time()
cfg = json.load(sys.stdin)
curve = cfg.get("curve", "q26")
trials = cfg.get("trials", 100)
R = cfg.get("R", 2048)
dp_bits = cfg.get("dp_bits", 3)
seed = cfg.get("seed", 1)
prod_checks = cfg.get("production_checks", 2)
budget = cfg.get("seconds", 14)
M_WALKS = cfg.get("walks", 16)
K_TEST = cfg.get("k_test", 1)
canon_random = cfg.get("canon_random", None)
canon_stress = cfg.get("canon_stress", None)
run_toy = cfg.get("toy", True)
dump = cfg.get("dump", False)
KAPPA = 1.005
H1_FLOOR = 0.39003
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


class ClassWalk:
    """The class walk of proof.md section 2 on curve E: twisted R-adding table,
    parity key, exact canonical twist index (section 1.1 reproduces it, Part A),
    h = bits 1..log2(R) of K, distinguished test on the next dp bits (the toy
    analogue of the production test on limb 1 of K)."""

    def __init__(self, E, R, dp_bits, T, rng, point=None):
        self.E, self.R, self.dp_bits = E, R, dp_bits
        self.rbits = R.bit_length() - 1
        assert 1 << self.rbits == R
        self.dpmask = (1 << dp_bits) - 1
        r = E.r
        lam = E.lam
        self.units = [(-1) ** s * pow(lam, j, r) % r for s in (0, 1) for j in (0, 1, 2)]
        self.uinv = [pow(u, -1, r) for u in self.units]
        self.T = T
        self.table = []
        while len(self.table) < R:
            c, d = rng.randrange(r), rng.randrange(r)
            M = point(c, d) if point else E.add(E.mul(c, E.G), E.mul(d, T))
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

    def canon(self, P):
        """(K, t, c): parity key K = y if y even else q - y; c = units[t] P is the
        canonical member (beta^j x minimal, y = K), t = 3s + j, s = [y odd]."""
        q, beta = self.E.q, self.E.beta
        x, y = P
        bx = beta * x % q
        b2x = (-x - bx) % q
        j, cx = 0, x
        if bx < cx:
            cx, j = bx, 1
        if b2x < cx:
            cx, j = b2x, 2
        if y & 1:
            return q - y, 3 + j, (cx, q - y)
        return y, j, (cx, y)

    def target(self, K, t):
        return self.table[(K >> 1) & (self.R - 1)][t]

    def is_dp(self, K):
        return ((K >> (1 + self.rbits)) & self.dpmask) == 0


class Trail:
    __slots__ = ("P", "a", "b", "ts", "kcp", "tid", "A", "B", "active")


def lockstep_body(W, w, t, k, Pref, Lmax, on_arrival, stats):
    """One walk's work in batch t, up to (not including) its batched addition.

    Returns ("step", target) or ("dp", K, canonical point, canonical coefficients)
    or ("cap",) / ("exc",) / ("trav",). Exactly the order of the section-5 loop:
    distinguished test, checkpoint test (test batches, t = 0 mod k), refresh
    (refresh batches, t = 0 mod Pref: cap test, then store K_cp). A detected cycle
    is traversed synchronously, its minimum-key member doubled, and the body is
    re-entered for the escape point in the same batch. Coefficients (w.A, w.B) are
    tracked when not None (replay re-walks)."""
    E, r = W.E, W.E.r
    track = w.A is not None
    while True:
        K, tt, c = W.canon(w.P)
        if W.is_dp(K):
            if track:
                u = W.units[tt]
                return ("dp", K, c, (u * w.A % r, u * w.B % r))
            return ("dp", K, c, None)
        if t % k == 0 and w.kcp is not None and w.kcp == K:
            stats["detections"] += 1
            best = (K, w.P, w.A, w.B)
            Q, qa, qb = w.P, w.A, w.B
            n = 0
            while True:
                n += 1
                if n > Pref:
                    return ("trav",)
                KQ, tq, _ = W.canon(Q)
                TP, tc, td = W.target(KQ, tq)
                if Q[0] == TP[0]:
                    return ("exc",)
                Q = E.add(Q, TP)
                if track:
                    qa, qb = (qa + tc) % r, (qb + td) % r
                stats["traversal_ops"] += 1
                KQ = W.canon(Q)[0]
                on_arrival(KQ)
                if KQ == K:
                    break
                if KQ < best[0]:
                    best = (KQ, Q, qa, qb)
            _, Pm, am, bm = best
            w.P = E.add(Pm, Pm)
            stats["doublings"] += 1
            if w.P is None:
                return ("exc",)
            if track:
                w.A, w.B = 2 * am % r, 2 * bm % r
            on_arrival(W.canon(w.P)[0])
            w.kcp = None
            continue
        if t % Pref == 0:
            if t - w.ts >= Lmax:
                return ("cap",)
            w.kcp = K
        return ("step", W.target(K, tt))


def replay(W, point, a, b, ts, k, Pref, Lmax, stats):
    """Re-walk a trail from its stored (a, b, start batch) with the same schedule:
    the same batch indices give exactly the same sequence of actions, so the replay
    ends at the same distinguished class with tracked coefficients."""
    w = Trail()
    w.P, w.A, w.B, w.ts, w.kcp = point(a, b), a, b, ts, None
    t = ts
    E, r = W.E, W.E.r
    while True:
        res = lockstep_body(W, w, t, k, Pref, Lmax, lambda K: None, stats)
        if res[0] != "step":
            return res
        TP, tc, td = res[1]
        if w.P[0] == TP[0]:
            return ("exc",)
        w.P = E.add(w.P, TP)
        w.A, w.B = (w.A + tc) % r, (w.B + td) % r
        t += 1


def toy_trial(E, rng, R, dp_bits, M, k):
    """One complete run of the scaled production pipeline on the toy curve."""
    r = E.r
    xs = rng.randrange(1, r)
    T = E.fmul(xs)

    def point(a, b):
        return E.fmul(a + b * xs)

    W = ClassWalk(E, R, dp_bits, T, rng, point)
    Pref = 64 * k
    Lmax = 20 << dp_bits
    N = (r - 1) // 6
    Np = N * R / (R - 1)
    S1 = KAPPA * math.sqrt(2 * Np * math.log(1 / (1 - H1_FLOOR)))
    B1 = math.ceil(S1 / M)
    D_cap = math.ceil(2.16 * (B1 * M + M * (1 << dp_bits)) / (1 << dp_bits))
    stats = dict(detections=0, traversal_ops=0, doublings=0, abandoned=0, exceptional=0,
                 batched_ops=0, trails=0, matches=0, probes=0)
    first_seen = {}
    st = {"arrivals": 0, "tc": None, "tid": 0}
    cur = [0]

    def on_arrival(K):
        st["arrivals"] += 1
        prev = first_seen.get(K)
        if prev is None:
            first_seen[K] = cur[0]
        elif prev != cur[0] and st["tc"] is None:
            st["tc"] = st["arrivals"]

    def start(w, t):
        while True:
            a, b = rng.randrange(r), rng.randrange(r)
            P = point(a, b)
            if P is not None:
                break
        st["tid"] += 1
        w.P, w.a, w.b, w.ts, w.kcp, w.tid, w.A, w.B, w.active = P, a, b, t, None, st["tid"], None, None, True
        cur[0] = w.tid
        stats["trails"] += 1
        on_arrival(W.canon(P)[0])

    # Report log + linear-probing index, as in production (proof.md section 2):
    # m slots (a power of two, load <= 0.149 at capacity), slot = (fingerprint, log index),
    # probe cap 256 (overflow = failure). Toy keys are short, so the toy index and
    # fingerprint come from a multiplicative hash of K instead of separate key bits.
    m = 1
    while m * 0.149 < D_cap:
        m *= 2
    slots = [None] * m
    rlog = []
    outcome = [None]

    def report(w, K, t):
        """Distinguished class: probe the index; on a match re-walk both trails and solve."""
        hk = (K * 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
        i, fp = hk & (m - 1), hk >> 40
        for _ in range(256):
            stats["probes"] += 1
            sl = slots[i]
            if sl is None:
                break
            if sl[0] == fp and rlog[sl[1]][0] == K:
                stats["matches"] += 1
                if stats["matches"] > 8:
                    outcome[0] = False
                    return
                _, a2, b2, t2 = rlog[sl[1]]
                r1 = replay(W, point, w.a, w.b, w.ts, k, Pref, Lmax, {"detections": 0, "traversal_ops": 0, "doublings": 0})
                r2 = replay(W, point, a2, b2, t2, k, Pref, Lmax, {"detections": 0, "traversal_ops": 0, "doublings": 0})
                if r1[0] == "dp" and r2[0] == "dp" and r1[2] == r2[2]:
                    (A1, B1_), (A2, B2_) = r1[3], r2[3]
                    if (B1_ - B2_) % r:
                        x = (A2 - A1) * pow(B1_ - B2_, -1, r) % r
                        if E.fmul(x) == T:
                            outcome[0] = (x == xs)
                return
            i = (i + 1) & (m - 1)
        else:
            outcome[0] = False
            return
        if len(rlog) >= D_cap:
            outcome[0] = False
            return
        slots[i] = (fp, len(rlog))
        rlog.append((K, w.a, w.b, w.ts))

    walks = []
    for _ in range(M):
        w = Trail()
        start(w, 0)
        walks.append(w)
    t = 0
    phase1 = True
    while outcome[0] is None:
        if phase1 and t >= B1:
            phase1 = False
        live = 0
        for w in walks:
            if not w.active:
                continue
            cur[0] = w.tid
            while True:
                res = lockstep_body(W, w, t, k, Pref, Lmax, on_arrival, stats)
                if res[0] == "step":
                    break
                if res[0] == "dp":
                    report(w, res[1], t)
                    if outcome[0] is not None:
                        break
                else:
                    stats["abandoned" if res[0] in ("cap", "trav") else "exceptional"] += 1
                if phase1:
                    start(w, t)
                else:
                    w.active = False
                    break
            if outcome[0] is not None:
                break
            if not w.active:
                continue
            TP = res[1][0]
            if w.P[0] == TP[0]:
                stats["exceptional"] += 1
                if phase1:
                    start(w, t + 1)
                else:
                    w.active = False
                continue
            w.P = W.E.add(w.P, TP)
            stats["batched_ops"] += 1
            on_arrival(W.canon(w.P)[0])
            live += 1
        if outcome[0] is not None:
            break
        if not phase1 and live == 0:
            outcome[0] = False
            break
        t += 1
        if time.time() - T0 > budget:
            return None
    return dict(ok=outcome[0], tc=st["tc"], arrivals=st["arrivals"], S1=S1, B1=B1, **stats)


M64 = (1 << 64) - 1
GQ = (PROD["beta"] << 448) // PROD["q"]
GL = [(GQ >> (64 * j)) & M64 for j in range(7)]
DQ = (1 << 448) // PROD["q"]
DINT, DFRAC = DQ >> 64, DQ & M64
THR_GAP = 41
THR_B = 31
STRIDE = 304
DP_MASK0 = 0xff80000000000000


def exact_j(x):
    q, b = PROD["q"], PROD["beta"]
    v = [x, b * x % q, b * b % q * x % q]
    return min(range(3), key=lambda j: v[j])


def fast_j(x):
    """Section 1.1 rule. Returns j, or None for 'fall back to the exact rule'.
    F1 = L' (12 limb products), F0 = 9 X + umulh(X, Dfrac), F2 = -(F0 + F1) mod 2^64."""
    xl = [(x >> (64 * i)) & M64 for i in range(6)]
    L = 0
    for i in range(6):
        L += (xl[i] * GL[5 - i]) >> 64
        L += (xl[i] * GL[6 - i]) & M64
    L &= M64
    X = xl[5]
    F0 = (DINT * X + ((X * DFRAC) >> 64)) & M64
    F2 = (-(F0 + L)) & M64
    if L + 14 > M64 or F2 < THR_B:
        return None
    vals = [F0, L, F2]
    jm = min(range(3), key=lambda j: vals[j])
    gap = min(vals[j] - vals[jm] for j in range(3) if j != jm)
    if gap < THR_GAP:
        return None
    return jm


def fractions_exact(x):
    """Exact 2^64 * v_j / p for j = 0, 1, 2 (as Fractions), for error statistics."""
    q, b = PROD["q"], PROD["beta"]
    v = [x, b * x % q, b * b % q * x % q]
    return [Fr(vj << 64, q) for vj in v]


def _flags_sub(a, b):
    """NZCV of cmp a, b (subs xzr, a, b) on 64-bit registers."""
    res = (a - b) & M64
    return (res >> 63, int(res == 0), int(a >= b), (((a ^ b) & (a ^ res)) >> 63) & 1)


def _flags_add(a, b):
    """NZCV of cmn a, b (adds xzr, a, b) on 64-bit registers."""
    tot = a + b
    res = tot & M64
    return (res >> 63, int(res == 0), tot >> 64, ((~(a ^ b) & (a ^ res)) >> 63) & 1)


def _nzcv(k):
    return ((k >> 3) & 1, (k >> 2) & 1, (k >> 1) & 1, k & 1)


_COND = dict(eq=lambda f: f[1] == 1, hs=lambda f: f[2] == 1, lo=lambda f: f[2] == 0,
             hi=lambda f: f[2] == 1 and f[1] == 0, ls=lambda f: f[2] == 0 or f[1] == 1)




def listing_rows_1_3(x):
    """Rows 1-3 of the section-5 listing executed instruction by instruction on 64-bit
    registers with NZCV flags. Resident registers: x22..x28 = G0..G6, x29 = Dfrac,
    x30 = Dint = 9, x14 = 41, x15 = 2. Returns (fallback taken, j in x6, (A, B, C))
    with A, B, C the three section-1.1 tests evaluated on exact integers."""
    reg = {"x22": GL[0], "x23": GL[1], "x24": GL[2], "x25": GL[3], "x26": GL[4], "x27": GL[5],
           "x28": GL[6], "x29": DFRAC, "x30": DINT, "x14": THR_GAP, "x15": 2}
    xl = [(x >> (64 * i)) & M64 for i in range(6)]
    for i in range(6):
        reg["x%d" % i] = xl[i]
    umulh = lambda a, b: (reg[a] * reg[b]) >> 64
    # row 1 (after the three ldp of x limbs into x0..x5)
    reg["x9"] = umulh("x5", "x22")
    reg["x8"] = umulh("x4", "x23"); reg["x9"] = (reg["x9"] + reg["x8"]) & M64
    reg["x9"] = (reg["x5"] * reg["x23"] + reg["x9"]) & M64
    for (xa, ga, xb, gb) in (("x3", "x24", "x4", "x24"), ("x2", "x25", "x3", "x25"),
                             ("x1", "x26", "x2", "x26"), ("x0", "x27", "x1", "x27")):
        reg["x8"] = umulh(xa, ga); reg["x9"] = (reg["x9"] + reg["x8"]) & M64
        reg["x9"] = (reg[xb] * reg[gb] + reg["x9"]) & M64
    reg["x9"] = (reg["x0"] * reg["x28"] + reg["x9"]) & M64
    # row 2
    reg["x0"] = umulh("x5", "x29")
    reg["x0"] = (reg["x5"] * reg["x30"] + reg["x0"]) & M64
    reg["x2"] = (reg["x0"] + reg["x9"]) & M64
    reg["x2"] = (-reg["x2"]) & M64
    # row 3
    f = _flags_sub(reg["x9"], reg["x0"])
    lo = _COND["lo"](f)
    reg["x3"] = reg["x9"] if lo else reg["x0"]
    reg["x4"] = reg["x0"] if lo else reg["x9"]
    reg["x6"] = 1 if lo else 0
    f = _flags_sub(reg["x2"], reg["x4"])
    reg["x7"] = reg["x2"] if _COND["lo"](f) else reg["x4"]
    f = _flags_sub(reg["x2"], reg["x3"])
    lo = _COND["lo"](f)
    reg["x7"] = reg["x3"] if lo else reg["x7"]
    reg["x3"] = reg["x2"] if lo else reg["x3"]
    reg["x6"] = reg["x15"] if lo else reg["x6"]
    reg["x7"] = (reg["x7"] - reg["x3"]) & M64
    f = _flags_add(reg["x9"], 14)                                   # cmn x9,#14
    f = _flags_sub(reg["x2"], 31) if _COND["lo"](f) else _nzcv(0)    # ccmp x2,#31,#0,cc
    f = _flags_sub(reg["x7"], reg["x14"]) if _COND["hs"](f) else _nzcv(0)  # ccmp x7,x14,#0,hs
    fallback = _COND["lo"](f)                                        # b.lo fallback
    # the three tests on exact integers
    L = 0
    for i in range(6):
        L += (xl[i] * GL[5 - i]) >> 64
        L += (xl[i] * GL[6 - i]) & M64
    L &= M64
    X = x >> 320
    F0 = DINT * X + ((X * DFRAC) >> 64)
    F2 = (-(F0 + L)) % (1 << 64)
    v = sorted([F0, L, F2])
    tests = (L + 14 >= 1 << 64, F2 < THR_B, v[1] - v[0] < THR_GAP)
    return fallback, reg["x6"], tests


def check_fast_j(n_random, n_stress):
    """fast_j and the listing against the exact rule; error statistics of F0, F1, F2
    against the section-1.1 bounds e0 in (-11.85, 0], e1 in (-12.13, 0], e2 = -e0-e1."""
    q, b = PROD["q"], PROD["beta"]
    rng2 = random.Random(seed + 10 ** 7)
    res = dict(random=0, random_fallbacks=0, random_mismatches=0,
               stress=0, stress_fallbacks=0, stress_mismatches=0,
               listing_fallback_mismatches=0, listing_j_mismatches=0,
               stress_test_A_range=0, stress_test_B_wrap=0, stress_test_C_gap=0,
               err_bound_violations=0)
    emin, emax = [0.0] * 3, [0.0] * 3

    def one(x, kind):
        f = fast_j(x)
        if f is None:
            res[kind + "_fallbacks"] += 1
        elif f != exact_j(x):
            res[kind + "_mismatches"] += 1
        fb, j, tests = listing_rows_1_3(x)
        if fb != any(tests) or fb != (f is None):
            res["listing_fallback_mismatches"] += 1
        if not fb and j != f:
            res["listing_j_mismatches"] += 1
        if kind == "stress":
            for name, hit in zip(("A_range", "B_wrap", "C_gap"), tests):
                res["stress_test_" + name] += hit
        if kind == "random":
            T = fractions_exact(x)
            xl5 = x >> 320
            L = 0
            xl = [(x >> (64 * i)) & M64 for i in range(6)]
            for i in range(6):
                L += (xl[i] * GL[5 - i]) >> 64
                L += (xl[i] * GL[6 - i]) & M64
            L &= M64
            F0 = DINT * xl5 + ((xl5 * DFRAC) >> 64)
            e0, e1 = F0 - T[0], L - T[1]
            if T[1] < 13:
                return
            F2 = (-(F0 + L)) % (1 << 64)
            e2 = F2 - T[2]
            if e2 < -(1 << 63):
                e2 += 1 << 64
            for i, e in enumerate((e0, e1, e2)):
                emin[i], emax[i] = min(emin[i], float(e)), max(emax[i], float(e))
            if not (Fr(-1185, 100) < e0 <= 0 and Fr(-1213, 100) < e1 <= 0 and e2 == -e0 - e1):
                res["err_bound_violations"] += 1

    for _ in range(n_random):
        if time.time() - T0 > budget:
            break
        res["random"] += 1
        one(rng2.randrange(1, q), "random")
    mults = [b - 1, b * b - 1, b - b * b, 1 + b, b]
    for _ in range(n_stress):
        if time.time() - T0 > budget:
            break
        kk = rng2.randrange(300, 381)
        d = rng2.randrange(-(1 << kk), 1 << kk)
        mlt = mults[rng2.randrange(len(mults))]
        x = d * pow(mlt % q, -1, q) % q
        if mlt == 1 + b:
            x = (q - d) * pow(mlt % q, -1, q) % q
        if x == 0:
            continue
        res["stress"] += 1
        one(x, "stress")
    # inputs aimed at each fallback window: F1 near 0 / 2^64, F2 near 0 / 2^64
    res["error_ranges_random"] = [[round(emin[i], 3), round(emax[i], 3)] for i in range(3)]
    return res


def check_rows_4_5(n_random, n_boundary):
    """Rows 4 and 5 (ordinary batches) and row 4' (refresh batches) emulated limb by
    limb, with the test-batch checkpoint compare and the invalid marker.
    Row 4:  ldp x0,x1,[x19,#48]; subs x2,x10,x0; sbcs x3,x11,x1; ands x8,x0,#1;
            csel x0,x2,x0,ne; csel x1,x3,x1,ne
    Row 4': ldp x0,x1,[x19,#48]; ldr x7,[x19,#64]; subs x2,x10,x0; sbcs x3,x11,x1;
            sbcs x9,x12,x7; ands x8,x0,#1; csel x0,x2,x0,ne; csel x1,x3,x1,ne; csel x7,x9,x7,ne
    Row 5:  ubfx x4,x0,#1,#40; cbz x1,report; add x5,x6,x8,lsl #2; add x5,x5,x4,lsl #3;
            add x5,x21,x5,lsl #8; str x5,[x19,#96]
    Compared with K = y if y even else p - y: limbs 0-1 (rows 4/5), limbs 0-2 (row 4'),
    s = [y odd], h = bits 1..40 of K, distinguished iff limb 1 of K = 0, address
    table + 256 (8h + 4s + j); the test-batch compare (ldr x2,[x19,#104]; cmp x0,x2;
    b.eq) never matches the invalid marker 1 because bit 0 of K is 0."""
    q = PROD["q"]
    pl = [(q >> (64 * i)) & M64 for i in range(4)]
    rng3 = random.Random(seed + 2 * 10 ** 7)
    res = dict(random=0, boundary=0, distinguished_keys=0, mismatches=0, row4p_mismatches=0,
               k_bit0_nonzero=0, invalid_marker_matches=0, address_checks=0, address_mismatches=0)
    cases = [rng3.randrange(1, q) for _ in range(n_random)]
    half = q >> 1
    cases += [1, 2, 3, q - 1, q - 2, q - 3] + [half + d for d in range(-3, 4)]
    picks = [0, 1, M64, M64 - 1] + pl + [v + 1 for v in pl] + [v - 1 for v in pl]
    while len(cases) < n_random + 13 + n_boundary:
        y = rng3.randrange(1, q) >> 256 << 256
        for i in range(4):
            y |= (picks[rng3.randrange(len(picks))] if rng3.random() < 0.7 else rng3.getrandbits(64)) << (64 * i)
        if 0 < y < q:
            cases.append(y)
    n_dp = max(1, n_boundary // 4)
    while len(cases) < n_random + 13 + n_boundary + n_dp:
        K = (rng3.randrange(2, q) >> 128 << 128) | (rng3.getrandbits(55) & ~1)
        if 0 < K < q:
            cases.append(K if rng3.random() < 0.5 else q - K)
    table = (rng3.getrandbits(12) << 52) & M64
    for i, y in enumerate(cases):
        res["random" if i < n_random else "boundary" if i < n_random + 13 + n_boundary
            else "distinguished_keys"] += 1
        yl = [(y >> (64 * k)) & M64 for k in range(3)]
        d, borrow = [], 0
        for k in range(3):
            v = pl[k] - yl[k] - borrow
            borrow = 1 if v < 0 else 0
            d.append(v & M64)
        s_bit = yl[0] & 1
        K2 = [d[k] if s_bit else yl[k] for k in range(2)]
        K3 = [d[k] if s_bit else yl[k] for k in range(3)]
        h = (K2[0] >> 1) & ((1 << 40) - 1)
        dist = K2[1] == 0 and (K2[0] & DP_MASK0) == 0
        Kx = y if y % 2 == 0 else q - y
        if Kx & 1 or K2[0] & 1:
            res["k_bit0_nonzero"] += 1
        if K2[0] == 1:
            res["invalid_marker_matches"] += 1
        if (s_bit != (y & 1) or K2 != [(Kx >> (64 * k)) & M64 for k in range(2)]
                or h != (Kx >> 1) % (1 << 40) or dist != ((Kx >> 55) % (1 << 73) == 0)):
            res["mismatches"] += 1
        if K3[0] | K3[1] << 64 | K3[2] << 128 != Kx % (1 << 192):
            res["row4p_mismatches"] += 1
        for j in range(3):
            x5 = (j + (s_bit << 2)) & M64
            x5 = (x5 + (h << 3)) & M64
            x5 = (table + (x5 << 8)) & M64
            want = table + 256 * (8 * ((Kx >> 1) % (1 << 40)) + 4 * (y % 2) + j)
            res["address_checks"] += 1
            if x5 != want:
                res["address_mismatches"] += 1
    idx = set()
    for hh in (0, 1, 2, 12345, (1 << 40) - 2, (1 << 40) - 1):
        for t in range(6):
            idx.add(8 * hh + 4 * (t // 3) + t % 3)
    res["address_map_injective_in_2_43_slots"] = len(idx) == 36 and max(idx) < 1 << 43
    res["table_end_offset_fits_64_bits"] = 256 * (1 << 43) < 1 << 64
    return res


LISTING = {
    "row1": ["ldp x0,x1,[x19,#304]!", "ldp x2,x3,[x19,#16]", "ldp x4,x5,[x19,#32]",
             "umulh x9,x5,x22", "umulh x8,x4,x23", "add x9,x9,x8", "madd x9,x5,x23,x9",
             "umulh x8,x3,x24", "add x9,x9,x8", "madd x9,x4,x24,x9",
             "umulh x8,x2,x25", "add x9,x9,x8", "madd x9,x3,x25,x9",
             "umulh x8,x1,x26", "add x9,x9,x8", "madd x9,x2,x26,x9",
             "umulh x8,x0,x27", "add x9,x9,x8", "madd x9,x1,x27,x9", "madd x9,x0,x28,x9"],
    "row2": ["umulh x0,x5,x29", "madd x0,x5,x30,x0", "add x2,x0,x9", "neg x2,x2"],
    "row3": ["cmp x9,x0", "csel x3,x9,x0,lo", "csel x4,x0,x9,lo", "cset x6,lo",
             "cmp x2,x4", "csel x7,x2,x4,lo", "cmp x2,x3", "csel x7,x3,x7,lo", "csel x3,x2,x3,lo",
             "csel x6,x15,x6,lo", "sub x7,x7,x3",
             "cmn x9,#14", "ccmp x2,#31,#0,cc", "ccmp x7,x14,#0,hs", "b.lo fallback"],
    "row4": ["ldp x0,x1,[x19,#48]", "subs x2,x10,x0", "sbcs x3,x11,x1", "ands x8,x0,#1",
             "csel x0,x2,x0,ne", "csel x1,x3,x1,ne"],
    "row5": ["ubfx x4,x0,#1,#40", "cbz x1,dpcheck", "add x5,x6,x8,lsl #2", "add x5,x5,x4,lsl #3",
             "add x5,x21,x5,lsl #8", "str x5,[x19,#96]"],
    "row7": ["cmp x19,x20", "b.ne loop"],
    "test_extra": ["ldr x2,[x19,#104]", "cmp x0,x2", "b.eq detect"],
    "row4_refresh": ["ldp x0,x1,[x19,#48]", "ldr x7,[x19,#64]", "subs x2,x10,x0", "sbcs x3,x11,x1",
                     "sbcs x9,x12,x7", "ands x8,x0,#1", "csel x0,x2,x0,ne", "csel x1,x3,x1,ne",
                     "csel x7,x9,x7,ne"],
    "refresh_extra": ["ldp x2,x3,[x19,#128]", "cmp x13,x2", "sbcs xzr,x16,x3", "b.hs abandon",
                      "stp x0,x1,[x19,#104]", "str x7,[x19,#120]"],
    "dpcheck": ["tst x0,#0xff80000000000000", "b.ne back_to_row5", "b report"],
}
REGS = {"x%d" % i for i in range(31)} | {"xzr"}
CONDS = {"eq", "ne", "hs", "cs", "lo", "cc", "mi", "pl", "vs", "vc", "hi", "ls", "ge", "lt", "gt", "le"}


def bitmask_imm(v):
    """True iff v is an AArch64 64-bit logical (bitmask) immediate."""
    if v in (0, M64):
        return False
    for e in (2, 4, 8, 16, 32, 64):
        mask = (1 << e) - 1
        chunk = v & mask
        if any(((v >> (e * i)) & mask) != chunk for i in range(64 // e)):
            continue
        for rot in range(e):
            c = ((chunk >> rot) | (chunk << (e - rot))) & mask
            if c and c != mask and (c & (c + 1)) == 0:
                return True
    return False


def check_encodable(ins):
    """Operand-range check of one listed instruction (the forms used in proof.md)."""
    op, _, rest = ins.partition(" ")
    if op.startswith("b."):
        return op[2:] in CONDS
    if op in ("b", "bl"):
        return True
    args = []
    depth, cur = 0, ""
    for ch in rest:
        if ch == "[":
            depth += 1
        if ch == "]":
            depth -= 1
        if ch == "," and depth == 0:
            args.append(cur.strip())
            cur = ""
        else:
            cur += ch
    args.append(cur.strip())
    imm = lambda a: int(a[1:], 0)
    if op in ("ldp", "stp", "ldr", "str"):
        regs, mem = args[:-1], args[-1]
        if not all(r in REGS and r != "xzr" or (op.startswith("st") and r in REGS) for r in regs):
            return False
        pre = mem.endswith("!")
        base, _, off = mem.strip("[]!").partition(",")
        off = imm(off.strip()) if off else 0
        if base.strip() not in REGS or off % 8:
            return False
        if op in ("ldp", "stp"):
            return -512 <= off <= 504 and len(regs) == 2
        return (0 <= off <= 32760) if not pre else (-256 <= off <= 255)
    if op in ("umulh", "mul"):
        return len(args) == 3 and all(a in REGS for a in args)
    if op == "madd":
        return len(args) == 4 and all(a in REGS for a in args)
    if op in ("add", "sub", "adds", "subs", "cmp", "cmn"):
        if args[-1].startswith("#"):
            v = imm(args[-1])
            return all(a in REGS for a in args[:-1]) and (0 <= v <= 4095)
        if "lsl" in args[-1]:
            sh = int(args[-1].split("#")[1])
            return all(a in REGS for a in args[:-1]) and 0 <= sh <= 63
        return all(a in REGS for a in args)
    if op in ("sbcs", "adcs", "adc", "sbc", "csel", "cset", "neg"):
        regs = [a for a in args if a in REGS]
        conds = [a for a in args if a in CONDS]
        return len(regs) + len(conds) == len(args)
    if op in ("ccmp", "ccmn"):
        a, b, nzcv, cond = args
        if a not in REGS or cond not in CONDS or not (0 <= imm(nzcv) <= 15):
            return False
        return b in REGS or (b.startswith("#") and 0 <= imm(b) <= 31)
    if op in ("ands", "tst", "and"):
        return all(a in REGS for a in args[:-1]) and args[-1].startswith("#") and bitmask_imm(imm(args[-1]))
    if op == "ubfx":
        lsb, w = imm(args[2]), imm(args[3])
        return args[0] in REGS and args[1] in REGS and 0 <= lsb <= 63 and 1 <= w <= 64 - lsb
    if op in ("cbz", "cbnz"):
        return args[0] in REGS
    return False


def check_listing():
    bad = [(row, ins) for row, seq in LISTING.items() for ins in seq if not check_encodable(ins)]
    counts = {row: len(seq) for row, seq in LISTING.items()}
    negative = [ins for ins in ("ccmp x2,#32,#0,cc", "ands x8,x0,#0", "ldp x0,x1,[x19,#508]!",
                                "ubfx x4,x0,#30,#40", "cmp x9,#4096")
                if check_encodable(ins)]
    return dict(rows=counts, ordinary_listed=sum(counts[r] for r in ("row1", "row2", "row3", "row4", "row5", "row7")) + 19,
                unencodable=bad, negative_controls_wrongly_accepted=negative)


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
    out["Dint_Dfrac"] = [DINT, hex(DFRAC)]
    out["log2_2^384_over_p"] = math.log2(2 ** 384 / q)
    rm1, fac = r - 1, {}
    for pr in [2, 3, 11, 19, 10177, 125527, 859267, 906349, 2508409, 2529403, 52437899, 254760293]:
        while rm1 % pr == 0:
            fac[pr] = fac.get(pr, 0) + 1
            rm1 //= pr
    order = r - 1
    for pr in fac:
        while order % pr == 0 and pow(2, order // pr, r) == 1:
            order //= pr
    out["r_minus_1_fully_factored"] = rm1 == 1
    out["log2_order_of_2_mod_r"] = round(math.log2(order), 4)
    xs = rng.randrange(1, r)
    T = E.mul(xs, G)
    W = ClassWalk(E, 4, 0, T, rng)
    ok_orbit = ok_key = ok_canon = ok_twist = ok_double = True
    for _ in range(nchecks * 4):
        k = rng.randrange(1, r)
        P = E.mul(k, G)
        members = [E.mul(u * k % r, G) for u in W.units]
        for idx in range(6):
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
        if len({c[0] for c in cans}) != 1 or cans[0][0] != Kp or Kp & 1:
            ok_key = False
        if len({c[2] for c in cans}) != 1:
            ok_canon = False
        for Mm, (K, t, c) in zip(members, cans):
            if E.mul(W.units[t], Mm) != c or fast_j(Mm[0]) not in (None, t % 3):
                ok_canon = False
            TP = W.target(K, t)[0]
            Mh = W.table[(K >> 1) & (W.R - 1)][0][0]
            if W.canon(E.add(Mm, TP))[2] != W.canon(E.add(c, Mh))[2]:
                ok_twist = False
            if W.canon(E.add(Mm, Mm))[2] != W.canon(E.add(c, c))[2]:
                ok_double = False
    out["orbit_size_6_and_unit_action"] = ok_orbit
    out["key_parity_is_class_invariant_bit0_zero"] = ok_key
    out["canonical_member_unique_and_fast_rule_agrees"] = ok_canon
    out["twisted_step_is_class_map"] = ok_twist
    out["doubling_is_class_map"] = ok_double
    out["twisted_table_entries_equal_u_inverse_M"] = all(
        W.table[h][t][0] == E.mul(W.uinv[t], W.table[h][0][0]) for h in range(W.R) for t in range(6))
    Wp = ClassWalk(E, 16, 64, T, rng)
    a, b = rng.randrange(r), rng.randrange(r)
    P = E.add(E.mul(a, G), E.mul(b, T))
    steps = 6 * nchecks
    trailA = []
    for _ in range(steps):
        K, t, c = Wp.canon(P)
        TP, tc, td = Wp.target(K, t)
        P = E.add(P, TP)
        a, b = (a + tc) % r, (b + td) % r
        trailA.append((P, a, b))
    out["production_steps_coefficients_ok"] = E.add(E.mul(a, G), E.mul(b, T)) == P
    i0 = steps // 3
    PA, aA, bA = trailA[i0]
    u = Wp.units[rng.randrange(6)]
    bB = rng.randrange(r)
    aB = (u * (aA + bA * xs) - bB * xs) % r
    PB = E.mul(u, PA)
    for _ in range(i0 + 1, steps):
        K, t, c = Wp.canon(PB)
        TP, tc, td = Wp.target(K, t)
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
    out["solve_x_equals_(e2-e1)/(f1-f2)"] = xsol == xs and E.mul(xsol, G) == T
    nr = 12000 * nchecks if canon_random is None else canon_random
    ns = 4000 * nchecks if canon_stress is None else canon_stress
    out["fast_canonical_index"] = check_fast_j(nr, ns)
    out["rows_4_5_and_4prime"] = check_rows_4_5(15000 * nchecks, 1500 * nchecks)
    out["part_a_seconds"] = round(time.time() - T0, 2)
    return out


CLAIM = (1261862, 10000)
C_INSTR, M_INSTR, ADDSUB_INSTR = 2000, 383, 32
R_PROD, M_PROD, K_PROD, P_PROD = 1 << 40, 1 << 22, 256, 1 << 14
L_MAX = 20 * 2 ** 73
THETA_LOG2 = -73
D_CAP_LOG2 = 54.35
SLOTS_LOG2 = 57
PROBE_CAP = 256
PROBE_PRICE = 0.015
LEDGER = [("row 1: L' (3 ldp of x, first pre-indexed; 6 umulh, 5 add, 6 madd; G resident)", 20),
          ("row 2: F0 = 9X + umulh(X, Dfrac) (umulh, madd), F2 = -(F0 + L') (add, neg)", 4),
          ("row 3: argmin and gap (11), merged fallback test cmn/ccmp/ccmp/b.lo (4)", 15),
          ("row 4: limbs 0-1 of y, p - y limbs 0-1, parity, select limbs 0-1 of K", 6),
          ("row 5: h (ubfx), distinguished test (cbz), table address, str", 6),
          ("row 6: entry address, table point and state loads, state store", 19),
          ("loop control (cmp, b.ne)", 2)]
TEST_EXTRA, REFRESH_EXTRA = 3, 12
FERMAT_INV = (380 + bin(PROD["q"] - 2).count("1") - 1) * M_INSTR
BATCH_ITEMS = [("zero-denominator test on the batch product", 9), ("batch meter charge", 10),
               ("Fermat inversion", FERMAT_INV), ("exponent scan (< 2,000)", 1999),
               ("resident constants and pointers (16 registers + 3 pointers)", 32),
               ("batch index increment and variant dispatch", 8)]
BATCH_SHARE = Fr(561, 10000)
UNBATCHED = 119
PREPROC_CAP_LOG2 = 56
PREPROC_CAP_FAIL = 2 ** -18
HASH_INIT = 2 ** SLOTS_LOG2 * 1.25 / 2000
HARD_EXTRAS_LOG2, RANDOM_LOG2 = 106.22, 92.2
RECORD = STRIDE
RECORD_FIELDS = [("x", 0, 48), ("y", 48, 48), ("next table-entry address", 96, 8),
                 ("K_cp limbs 0, 1", 104, 16), ("K_cp limb 2", 120, 8),
                 ("trail expiry batch index (2 limbs)", 128, 16),
                 ("start coefficient a", 144, 32), ("start coefficient b", 176, 32),
                 ("batch denominator", 208, 48), ("batch prefix product", 256, 48)]
LISTING_OFFSETS = [("ldp", 304), ("ldp", 16), ("ldp", 32), ("ldp", 48), ("ldr", 64), ("str", 96),
                   ("ldr", 104), ("stp", 104), ("str", 120), ("ldp", 128)]


def up(v, d):
    return (math.floor(v * 10 ** d) + 1) / 10 ** d


def dn(v, d):
    return (math.ceil(v * 10 ** d) - 1) / 10 ** d


def dec(f, d):
    assert (f * 10 ** d).denominator == 1
    return ("%." + str(d) + "f") % float(f)


_S2 = None


def stirling2_log2(n, k):
    """log2 S2(n, k) from one table of the exact recurrence, n <= 300, k <= 150."""
    global _S2
    if _S2 is None:
        _S2 = [[0] * 151 for _ in range(301)]
        _S2[0][0] = 1
        for i in range(1, 301):
            for j in range(1, min(i, 150) + 1):
                _S2[i][j] = j * _S2[i - 1][j] + _S2[i - 1][j - 1]
    v = _S2[n][k]
    return math.log2(v) if v else -math.inf


def cycle_rates():
    """Per-step rates: 2-cycles 1/(6R); other fruitless cycles of length L <= 64 need
    every partition index to occur >= 2 times: <= sum_L sum_b S2(L, b) R^-(L-b)."""
    lr = math.log2(R_PROD)
    other = 0.0
    for L in range(3, 65):
        for b_ in range(1, L // 2 + 1):
            other += 2 ** (stirling2_log2(L, b_) - lr * (L - b_))
    two = 1 / (6 * R_PROD)
    long_max = max(stirling2_log2(L, L // 2) - lr * (L - L // 2) for L in range(65, 301))
    return two, other, long_max


def cost_chain(c):
    listed = sum(n for _, n in LEDGER)
    step_instr = Fr(listed) + BATCH_SHARE
    step = 1 + step_instr / C_INSTR
    per_op = step + (Fr(TEST_EXTRA, K_PROD) + Fr(REFRESH_EXTRA, P_PROD)) / C_INSTR
    N = (PROD["r"] - 1) // 6
    Np = N * R_PROD / (R_PROD - 1)
    S0 = math.sqrt(2 * Np * math.log(1 / (1 - H1_FLOOR)))
    S_ops = KAPPA * S0
    lU1 = up(math.log2(float(per_op) * (S_ops + M_PROD)), 6)
    room = 2 ** c - 2 ** lU1 - 2 ** HARD_EXTRAS_LOG2
    l_head = dn(math.log2(room), 2) if room > 0 else -math.inf
    meter = up(min(1.0, 2 ** (RANDOM_LOG2 - l_head)), 7)
    meter3 = up(min(1.0, 3 * 2 ** (RANDOM_LOG2 - l_head)), 7)
    return dict(listed=listed, step_instr=step_instr, step=step, per_op=per_op, S0=S0, S_ops=S_ops,
                lU1=lU1, l_head=l_head, meter=meter, meter3=meter3,
                success=dn(H1_FLOOR - meter - PREPROC_CAP_FAIL - 0.00000011, 5),
                success3=dn(H1_FLOOR - meter3 - PREPROC_CAP_FAIL - 0.00000011, 5))


def cost():
    c = CLAIM[0] / CLAIM[1]
    k = cost_chain(c)
    below = cost_chain((CLAIM[0] - 1) / CLAIM[1])
    batch_work = sum(n for _, n in BATCH_ITEMS)
    unbatched = FERMAT_INV + 2000 + 5 * M_INSTR + 6 * ADDSUB_INSTR + 98 + 100 + 84
    two, other, long_max = cycle_rates()
    rate = two + other
    ops_total = 2 ** k["lU1"]
    detections = rate * ops_total + 2 ** -63 * ops_total / K_PROD
    trav_cost = detections * (0.05 + UNBATCHED * 3 + 1)
    mean_reports = (ops_total + M_PROD * (L_MAX + P_PROD)) * 2 ** THETA_LOG2
    starts = (mean_reports + M_PROD) * (2 ** 12 + 0.5)
    alpha = 2 ** (D_CAP_LOG2 - SLOTS_LOG2)
    xr = math.e ** (1 - alpha) * alpha
    probes_mean = 1 + sum(L * L * xr ** L for L in range(1, 4000))
    probe_tail_log2 = math.log2(PROBE_CAP * xr ** PROBE_CAP / (1 - xr))
    reports_random = (mean_reports + 16) * probes_mean * PROBE_PRICE
    dpcheck = 2 ** -64 * ops_total * 0.01
    fallbacks = 2 ** -55 * ops_total * 0.25
    rand = trav_cost + starts + fallbacks + reports_random + dpcheck + 2 ** 8
    hard = [M_PROD * (L_MAX + P_PROD) * UNBATCHED, 8 * 2 * ((L_MAX + P_PROD) * UNBATCHED + 2 ** 12),
            2 ** D_CAP_LOG2 * (0.5 + PROBE_CAP * PROBE_PRICE), HASH_INIT, 2 ** 40 * 2 ** 12 + 2 ** 41,
            M_PROD * 2 ** 12, M_PROD * (TEST_EXTRA + REFRESH_EXTRA) / C_INSTR, 2 ** 13 + 2 ** 9]
    fields = sorted(RECORD_FIELDS, key=lambda f: f[1])
    layout_ok = (all(o % 8 == 0 and n % 8 == 0 for _, o, n in fields)
                 and all(fields[i][1] + fields[i][2] <= fields[i + 1][1] for i in range(len(fields) - 1))
                 and fields[-1][1] + fields[-1][2] <= RECORD
                 and all(o % 8 == 0 and (o <= 504 if op in ("ldp", "stp") else o <= 32760)
                         for op, o in LISTING_OFFSETS))
    ok = lambda kk: kk["success"] >= 0.39 + 1e-5 and kk["success3"] >= 0.39 + 1e-5
    mem = 2 ** SLOTS_LOG2 * 16 + 2 ** D_CAP_LOG2 * 128 + (1 << 40) * 8 * 256 + M_PROD * RECORD + 2 ** 20
    return {
        "ledger_rows": [[name, n] for name, n in LEDGER], "listed_per_step": k["listed"],
        "per_batch_items": [[name, n] for name, n in BATCH_ITEMS], "per_batch_instr": batch_work,
        "per_batch_share_per_step": "< " + dec(BATCH_SHARE, 4),
        "step_instr": dec(k["step_instr"], 4), "price_step": dec(k["step"], 8),
        "test_batch_extra_per_walk": TEST_EXTRA, "refresh_batch_extra_per_walk": REFRESH_EXTRA,
        "phase1_rate_per_operation": str(float(k["per_op"])),
        "unbatched_operation_instr_bound": unbatched,
        "fermat_inversion_instr": FERMAT_INV,
        "log2_two_cycle_rate": round(math.log2(two), 3), "log2_other_short_cycle_rate": round(math.log2(other), 2),
        "log2_max_term_cycles_65_to_300": round(long_max, 1),
        "log2_expected_detections": round(math.log2(detections), 2),
        "log2_expected_traversal_units": round(math.log2(trav_cost), 2),
        "log2_expected_repeat_arrivals_from_detection_delay": round(math.log2(rate * ops_total * 2 * P_PROD), 2),
        "log2_random_total": round(math.log2(rand), 3),
        "log2_mean_reports": round(math.log2(mean_reports), 3), "log2_report_capacity": D_CAP_LOG2,
        "report_index_slots_log2": SLOTS_LOG2, "report_index_load_at_capacity": round(alpha, 4),
        "linear_probing_tail_base_e^(1-a)a": round(xr, 4), "probe_run_ge_256_log2_prob_per_op": round(probe_tail_log2, 1),
        "expected_probes_bound": round(probes_mean, 3), "log2_expected_restart_units": round(math.log2(starts), 2), "log2_hard_extras": round(math.log2(sum(hard)), 4),
        "record_layout": [[name, o, n] for name, o, n in fields],
        "log2_memory_bytes": round(math.log2(mem), 3),
        "log2_S0": "%.5f" % math.log2(k["S0"]), "log2_S_ops": "%.5f" % math.log2(k["S_ops"]),
        "log2_phase1_deterministic_upper": "%.6f" % k["lU1"], "claim": "%.4f" % c, "claim_fraction": list(CLAIM),
        "log2_headroom_lower": "%.2f" % k["l_head"], "meter_halt_upper": "%.7f" % k["meter"],
        "success_lower": "%.5f" % k["success"], "meter_halt_upper_3x_random": "%.7f" % k["meter3"],
        "success_lower_3x_random": "%.5f" % k["success3"],
        "checks": {
            "listed_per_step_72_budget_72.0561": k["listed"] == 72 and k["step_instr"] == Fr(720561, 10000),
            "per_batch_work_below_share_times_2_22": batch_work < BATCH_SHARE * M_PROD,
            "unbatched_operation_le_119_units": unbatched <= UNBATCHED * C_INSTR,
            "record_layout_fits_stride_and_offsets_encodable": layout_ok,
            "hard_extras_below_stated": math.log2(sum(hard)) < HARD_EXTRAS_LOG2,
            "random_below_2_92.2": math.log2(rand) < RANDOM_LOG2,
            "memory_below_2_62.3": mem < 2 ** 62.3,
            "preprocessing_cap_gap_vs_sampling": 2 ** 37 / (2 ** PREPROC_CAP_LOG2 - (HASH_INIT + 2 ** 52 + 2 ** 41 + 2 ** 34)) < PREPROC_CAP_FAIL,
            "report_capacity_2.16_times_mean": 2 ** D_CAP_LOG2 >= 2.16 * mean_reports,
            "probe_overflow_below_2^-300_total": probe_tail_log2 + math.log2(2 ** D_CAP_LOG2) < -300,
            "all_arrays_below_2^63_bytes": max(2 ** SLOTS_LOG2 * 16, 2 ** D_CAP_LOG2 * 128, (1 << 43) * 256) < 2 ** 63,
            "listing_encodable_and_72": not check_listing()["unencodable"] and check_listing()["ordinary_listed"] == 72,
            "long_cycle_terms_below_2_-900": long_max < -900,
            "refresh_period_multiple_of_test_period_and_64k": P_PROD % K_PROD == 0 and P_PROD == 64 * K_PROD,
            "lmax_multiple_of_refresh_period": L_MAX % P_PROD == 0,
            "success_margin_at_least_1e-5": ok(k),
            "claim_is_smallest_4_decimal_value": ok(k) and not ok(below),
        },
    }


def main():
    out = {"config": cfg}
    out["cost"] = cost()
    out["listing_encoding"] = check_listing()
    if prod_checks:
        out["part_a"] = part_a(prod_checks)
    if not run_toy:
        out["seconds"] = round(time.time() - T0, 2)
        print(json.dumps(out, indent=1))
        return
    E = make_toy(curve)
    r = E.r
    N = (r - 1) // 6
    Np = N * R / (R - 1)
    S0 = math.sqrt(2 * Np * math.log(1 / (1 - H1_FLOOR)))
    out["toy"] = dict(q=E.q, r=r, N=N, R=R, dp_bits=dp_bits, walks=M_WALKS, k_test=K_TEST,
                      refresh_period=64 * K_TEST, S0=round(S0, 1), S1=round(KAPPA * S0, 1),
                      mean_trail_over_sqrtN=round((1 << dp_bits) / math.sqrt(N), 4))
    res = []
    for _ in range(trials):
        tr = toy_trial(E, rng, R, dp_bits, M_WALKS, K_TEST)
        if tr is None:
            break
        res.append(tr)
    n = len(res)
    if n:
        S1 = KAPPA * S0
        tcs = [x["tc"] for x in res if x["tc"] is not None]
        ops = sum(x["batched_ops"] for x in res)
        out["toy_summary"] = dict(
            trials_completed=n,
            pipeline_success_rate=round(sum(1 for x in res if x["ok"]) / n, 4),
            binomial_stderr_at_0p39=round(math.sqrt(0.39 * 0.61 / n), 4),
            frac_tc_le_S1=round(sum(t <= S1 for t in tcs) / n, 4),
            frac_tc_le_S0=round(sum(t <= S0 for t in tcs) / n, 4),
            model_P_S1=round(1 - math.exp(-S1 * S1 / (2 * Np)), 4), H1_floor=H1_FLOOR,
            batched_ops=ops,
            detections_per_op=sum(x["detections"] for x in res) / max(1, ops),
            two_cycle_model_1_over_6R=1 / (6 * R),
            traversal_ops=sum(x["traversal_ops"] for x in res),
            abandoned=sum(x["abandoned"] for x in res),
            exceptional=sum(x["exceptional"] for x in res),
            matches=sum(x["matches"] for x in res),
            index_probes_per_report=round(sum(x["probes"] for x in res) / max(1, sum(x["trails"] for x in res)), 3),
        )
        if dump:
            out["tc_list"] = tcs
    out["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(out, indent=1))


main()
