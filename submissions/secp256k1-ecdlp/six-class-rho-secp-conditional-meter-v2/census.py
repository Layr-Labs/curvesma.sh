"""Independent toy-scale census of the six-class rho walk used by the
secp256k1-ecdlp submission.  This is a stand-alone reimplementation of the
toy part of the public v9 evidence: it does not execute any downloaded
candidate file.  Stdlib only."""
import argparse, json, math, random, time

CB = 7
TOYS = {
    "q25": (33556867, 33563419, 1),
    "q29": (536873251, 536839203, 3),
}

class Curve:
    def __init__(self, q, r, h):
        self.q, self.r, self.h = q, r, h
    def on_curve(self, P):
        return P is None or (P[1] * P[1] - P[0] ** 3 - CB) % self.q == 0
    def add(self, P, Q):
        q = self.q
        if P is None: return Q
        if Q is None: return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % q == 0: return None
            l = 3 * x1 * x1 * pow(2 * y1, -1, q) % q
        else:
            l = (y2 - y1) * pow(x2 - x1, -1, q) % q
        x3 = (l * l - x1 - x2) % q
        return (x3, (l * (x1 - x3) - y1) % q)
    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.q)
    def mul_raw(self, k, P):
        R = None
        while k:
            if k & 1: R = self.add(R, P)
            P = self.add(P, P); k >>= 1
        return R
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
        R, i = None, 0
        while k:
            R = self.add(R, self.fb[i][k & 255])
            k >>= 8; i += 1
        return R
    def phi(self, P, j=1):
        if P is None: return None
        return (P[0] * pow(self.beta, j % 3, self.q) % self.q, P[1])

def sqrt_mod(a, q):
    a %= q
    if a == 0: return 0
    if pow(a, (q - 1) // 2, q) != 1: return None
    s, t = 0, q - 1
    while t % 2 == 0: t //= 2; s += 1
    z = 2
    while pow(z, (q - 1) // 2, q) != q - 1: z += 1
    m, c, tt, x = s, pow(z, t, q), pow(a, t, q), pow(a, (t + 1) // 2, q)
    while tt != 1:
        i, t2 = 1, tt * tt % q
        while t2 != 1: t2 = t2 * t2 % q; i += 1
        b = pow(c, 1 << (m - i - 1), q)
        m, c, tt, x = i, b * b % q, tt * b * b % q, x * b % q
    return x

def cube_roots_of_unity(m):
    for g in range(2, 1000):
        c = pow(g, (m - 1) // 3, m)
        if c != 1: return [c, c * c % m]

def make_toy(name, rng):
    q, n, h = TOYS[name]
    r = n // h
    E = Curve(q, r, h)
    while True:
        x = rng.randrange(q)
        y = sqrt_mod(x ** 3 + CB, q)
        if y is None: continue
        G, P, k = None, (x, y), h
        while k:
            if k & 1: G = E.add(G, P)
            P = E.add(P, P); k >>= 1
        if G is not None and E.mul_raw(r, G) is None: break
    E.G = G; E.fixed_base(G)
    for b in cube_roots_of_unity(q):
        for l in cube_roots_of_unity(r):
            if E.mul_raw(l, G) == (G[0] * b % q, G[1]):
                E.beta, E.lam = b, l
    return E

class Walker:
    def __init__(self, E, R, dp_bits, T, rng, start_point):
        self.E, self.R, self.dp_bits = E, R, dp_bits
        self.rbits = R.bit_length() - 1
        assert 1 << self.rbits == R
        self.dpmask = (1 << dp_bits) - 1
        self.Lmax = 20 << dp_bits
        r = E.r; lam = E.lam
        self.units = [(-1) ** s * pow(lam, j, r) % r for s in (0, 1) for j in (0, 1, 2)]
        self.uinv = [pow(u, -1, r) for u in self.units]
        self.T = T; self.table = []
        while len(self.table) < R:
            c, d = rng.randrange(r), rng.randrange(r)
            M = start_point(c, d)
            if M is None: continue
            row = []
            for t in range(6):
                s, j = divmod(t, 3)
                Mt = E.phi(M, -j)
                if s: Mt = E.neg(Mt)
                ui = self.uinv[t]
                row.append((Mt, ui * c % r, ui * d % r))
            self.table.append(row)
        self.stats = dict(steps=0, doublings=0, esc2=0, esccp=0, abandoned=0, exceptional=0)
    def escape(self, P):
        self.stats["doublings"] += 1
        return self.E.add(P, P)
    def add_checked(self, P, Q):
        if P[0] == Q[0]: return None
        return self.E.add(P, Q)
    def canon(self, P):
        q, beta = self.E.q, self.E.beta
        x, y = P
        bx = beta * x % q
        b2x = (-x - bx) % q
        j, cx = 0, x
        if bx < cx: cx, j = bx, 1
        if b2x < cx: cx, j = b2x, 2
        ny = q - y
        if ny < y: return ny, 3 + j, (cx, ny)
        return y, j, (cx, y)
    def step_target(self, P, K, t):
        h = K & (self.R - 1)
        return self.table[h][t]
    def is_dp(self, K):
        return ((K >> self.rbits) & self.dpmask) == 0
    def run_trail(self, P, on_class=None):
        E = self.E
        k1 = k2 = None; cp = None; length = 0; pending_double = False
        while True:
            K, t, c = self.canon(P)
            if on_class is not None: on_class(K)
            if self.is_dp(K): return "dp", K, length
            if length >= self.Lmax:
                self.stats["abandoned"] += 1; return "cap", K, length
            if pending_double:
                P = self.escape(P)
                if P is None:
                    self.stats["exceptional"] += 1; return "exc", K, length
                pending_double = False; k1 = k2 = cp = None; length += 1; continue
            if k2 is not None and K == k2:
                self.stats["esc2"] += 1
                if K < k1:
                    P = self.escape(P)
                    if P is None:
                        self.stats["exceptional"] += 1; return "exc", K, length
                    k1 = k2 = cp = None; length += 1; continue
                pending_double = True
            elif cp is not None and K == cp:
                self.stats["esccp"] += 1
                best = (K, P); Q = P; trav = 0
                while True:
                    trav += 1
                    if trav > 64:
                        self.stats["abandoned"] += 1; return "trav", K, length
                    TP, tc, td = self.step_target(Q, *self.canon(Q)[:2])
                    Q = self.add_checked(Q, TP)
                    if Q is None:
                        self.stats["exceptional"] += 1; return "exc", K, length
                    length += 1; self.stats["steps"] += 1
                    KQ = self.canon(Q)[0]
                    if on_class is not None: on_class(KQ)
                    if KQ == K: break
                    if KQ < best[0]: best = (KQ, Q)
                _, P = best
                P = self.escape(P)
                if P is None:
                    self.stats["exceptional"] += 1; return "exc", K, length
                k1 = k2 = cp = None; length += 1; continue
            if length % 64 == 0: cp = K
            TP, tc, td = self.step_target(P, K, t)
            P = self.add_checked(P, TP)
            if P is None:
                self.stats["exceptional"] += 1; return "exc", K, length
            self.stats["steps"] += 1; length += 1; k2, k1 = k1, K

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curve", default="q29", choices=["q25", "q29"])
    ap.add_argument("--trials", type=int, default=1000)
    ap.add_argument("--R", type=int, default=2048)
    ap.add_argument("--dp", type=int, default=5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--kappa", type=float, default=1.01)
    ap.add_argument("--p0", type=float, default=0.3950)
    args = ap.parse_args()
    t0 = time.time()
    rng = random.Random(args.seed)
    E = make_toy(args.curve, rng)
    r = E.r
    xs = rng.randrange(1, r)
    T = E.fmul(xs)
    N = (r - 1) // 6
    Np = N * args.R / (args.R - 1)
    S0 = math.sqrt(2 * Np * math.log(1 / (1 - args.p0)))
    W = Walker(E, args.R, args.dp, T, rng, lambda c, d: E.fmul(c + d * xs))
    tcs = []
    visits_total = 0
    first_seen = {}
    for trial in range(args.trials):
        if time.time() - t0 > 900:
            break
        state = {"visits": 0, "tc": None, "trail": 0}
        def on_class(K):
            state["visits"] += 1
            prev = first_seen.get(K)
            if prev is None:
                first_seen[K] = state["trail"]
            elif prev != state["trail"] and state["tc"] is None:
                state["tc"] = state["visits"]
        # fresh collision table for each trial: T_c is a per-instance quantity
        first_seen.clear()
        while state["tc"] is None:
            a, b = rng.randrange(r), rng.randrange(r)
            P = E.fmul(a + b * xs)
            if P is None: continue
            state["trail"] += 1
            W.run_trail(P, on_class)
        tcs.append(state["tc"])
        visits_total += state["visits"]
    n = len(tcs)
    mean = sum(tcs) / n / math.sqrt(N)
    var = sum((t / math.sqrt(N) - mean) ** 2 for t in tcs) / max(1, n - 1)
    stderr = math.sqrt(var / n)
    model = 1 - math.exp(-S0 * S0 / (2 * Np))
    within = sum(t <= S0 for t in tcs) / n
    print(json.dumps({
        "curve": args.curve, "R": args.R, "dp_bits": args.dp, "seed": args.seed,
        "trials": n, "q": E.q, "r": r, "N": N,
        "mean_tc_over_sqrtN": mean, "stderr": stderr,
        "slowdown_vs_model": mean / math.sqrt(math.pi * Np / 2 / N) - 1,
        "model_mean_over_sqrtN": math.sqrt(math.pi * Np / 2) / math.sqrt(N),
        "S0_over_sqrtN": S0 / math.sqrt(N),
        "frac_tc_le_S0": within, "model_P_S0": model,
        "doublings_per_step": W.stats["doublings"] / max(1, W.stats["steps"]),
        "seconds": round(time.time() - t0, 2),
    }, indent=1))

if __name__ == "__main__":
    main()
