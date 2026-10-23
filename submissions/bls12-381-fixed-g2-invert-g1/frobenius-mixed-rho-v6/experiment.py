"""Toy v11 budget run (stdlib only): python3 -I experiment.py TAU SEED N. Protocol of RESULT-toy-v11c: ToyGT(2^18,
201), r_a 1024, S 458, DP spacing 4, PRF restarts, no cap, stop at verified output, max_steps 13000, branch key per trial
Random(SEED*100003+i). T = sqrt(2 N TAU / V), V = 1 - s^2 - (1-s)^2/r_a, the withdrawn §9.3 expression at s = 458/1024, kept as an assumed lower bound (not derived). Prints one JSON object."""
import hashlib, json, math, random, sys

M64, MP, MA, MB = (1 << 64) - 1, (1 << 127) - 1, 0x9E3779B97F4A7C15D1B54A32D192ED03, 0x632BE59BD9B4E019F1357BD8A2C6D4E3
KP = hashlib.sha256(b"gt-v11-toy").digest()
DEG, R_A, S, AB, DPM, MAXC, MAXS = 12, 1024, 458, 10, 3, 64, 13000


def is_prime(n):
    return n > 1 and all(n % p for p in range(2, int(n ** 0.5) + 1))


def factors(n):
    out, d = [], 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1
    return out + ([n] if n > 1 else [])


def toy_params(min_order):
    q = 2
    while True:
        q += 1
        if is_prime(q):
            big = [p for p in factors(q ** 4 - q * q + 1) if p >= min_order and p not in (2, 3)]
            if big:
                return q, min(big)


def trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def pmod(a, b, q):
    a = list(a)
    inv = pow(b[-1], q - 2, q)
    while len(a) >= len(b):
        c, sh = a[-1] * inv % q, len(a) - len(b)
        for i, bi in enumerate(b):
            a[sh + i] = (a[sh + i] - c * bi) % q
        a.pop()
        trim(a)
    return a


def pgcd(a, b, q):
    a, b = trim([c % q for c in a]), trim([c % q for c in b])
    while b:
        a, b = b, pmod(a, b, q)
    return a


class Fq12:
    def __init__(self, q, f):
        self.q, self.f = q, f
        self.one, self.x = [1] + [0] * 11, [0, 1] + [0] * 10

    def mul(self, a, b):
        q, f = self.q, self.f
        p = [0] * 23
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    p[i + j] += ai * bj
        for k in range(22, 11, -1):
            c = p[k] % q
            if c:
                for i in range(12):
                    p[k - 12 + i] -= c * f[i]
            p[k] = 0
        return [c % q for c in p[:12]]

    def pow(self, a, e):
        r, b = self.one, a
        while e:
            if e & 1:
                r = self.mul(r, b)
            e >>= 1
            if e:
                b = self.mul(b, b)
        return r

    def irreducible(self):
        q, fr = self.q, [self.x]
        for _ in range(12):
            fr.append(self.pow(fr[-1], q))
        if fr[12] != self.x:
            return False
        for p in (2, 3):
            t = list(fr[12 // p])
            t[1] = (t[1] - 1) % q
            if len(pgcd(t, self.f + [1], q)) != 1:
                return False
        return True


class FixedBase:
    def __init__(self, field, base, order):
        self.field, self.table, b = field, [], base
        for _ in range((order.bit_length() + 3) // 4):
            row = [field.one]
            for _ in range(15):
                row.append(field.mul(row[-1], b))
            self.table.append(row)
            b = field.mul(row[-1], b)

    def pow(self, e):
        r = None
        for i, row in enumerate(self.table):
            d = (e >> (4 * i)) & 15
            if d:
                r = row[d] if r is None else self.field.mul(r, row[d])
        return self.field.one if r is None else r


def mv(M, v, q):
    return [sum(x * y for x, y in zip(row, v)) % q for row in M]


def mm(A, B, q):
    cols = list(zip(*B))
    return [[sum(x * y for x, y in zip(row, col)) % q for col in cols] for row in A]


def mix64(x):
    z = x & M64
    z ^= z >> 30
    z = z * 0xBF58476D1CE4E5B9 & M64
    z ^= z >> 27
    z = z * 0x94D049BB133111EB & M64
    return z ^ (z >> 31)


def prf(msg, m):
    h0 = hashlib.sha256(KP + msg + b"\x00").digest()
    h1 = hashlib.sha256(KP + msg + b"\x01").digest()
    return int.from_bytes(h0 + h1, "big") % m


class Group:
    def __init__(self, min_order, seed=201):
        self.q, self.m = toy_params(min_order)
        q, m, rng = self.q, self.m, random.Random(seed)
        while True:
            c = [rng.randrange(q) for _ in range(12)]
            if c[0] and Fq12(q, c).irreducible():
                break
        self.F = F = Fq12(q, c)
        g = F.one
        while g == F.one:
            g = F.pow([rng.randrange(q) for _ in range(12)], (q ** 12 - 1) // m)
        self.gtab = FixedBase(F, g, m)
        self.qp = [q ** i for i in range(12)]
        xq = F.pow(F.x, q)
        cols = [F.pow(xq, i) if i else F.one for i in range(12)]
        self.sig = [[cols[i][j] for i in range(12)] for j in range(12)]
        self.sinv = self.sig
        for _ in range(10):
            self.sinv = mm(self.sinv, self.sig, q)
        self.lam = q % m
        self.lp = [pow(self.lam, k, m) for k in range(12)]
        self.li = [pow(self.lam, -k, m) for k in range(12)]
        self.started = 0

    def key(self, x):
        best, kb, v = None, 0, x
        for n in range(12):
            e = sum(a * b for a, b in zip(v, self.qp))
            if best is None or e < best:
                best, kb = e, n
            v = mv(self.sig, v, self.q)
        return (best * MA + MB) % MP, kb

    def canon(self, st):
        k = st[4]
        return st[3], self.lp[k] * st[1] % self.m, self.lp[k] * st[2] % self.m


class Walk:
    def __init__(self, G, zt, rng, kb):
        m, F = G.m, G.F
        self.G, self.m, self.kb = G, m, kb
        self.T, self.Ec, self.Ed = [], [], []
        for _ in range(R_A):
            c, d = rng.randrange(m), rng.randrange(m)
            v = F.mul(G.gtab.pow(c), zt.pow(d))
            rt, rc, rd = [], [], []
            for k in range(12):
                rt.append(v)
                rc.append(c * G.li[k] % m)
                rd.append(d * G.li[k] % m)
                v = mv(G.sinv, v, self.G.q)
            self.T.append(rt)
            self.Ec.append(rc)
            self.Ed.append(rd)

    def step(self, st):
        x, a, b, K, k = st
        G, m = self.G, self.m
        u = mix64((K & M64) ^ self.kb)
        if (u >> (64 - AB)) < S:
            x2 = G.F.mul(x, x)
            K2, k2 = G.key(x2)
            return (x2, 2 * a % m, 2 * b % m, K2, k2)
        h = u & (R_A - 1)
        x2 = G.F.mul(x, self.T[h][k])
        K2, k2 = G.key(x2)
        return (x2, (a + self.Ec[h][k]) % m, (b + self.Ed[h][k]) % m, K2, k2)

    def restart(self, st, e):
        x, a, b = st[0], st[1], st[2]
        x2 = self.G.F.mul(x, self.G.gtab.pow(e))
        K2, k2 = self.G.key(x2)
        return (x2, (a + e) % self.m, b, K2, k2)

    def escape(self, members):
        G, m = self.G, self.m
        s = min(members, key=lambda t: t[3])
        x = s[0]
        for _ in range(s[4]):
            x = mv(G.sig, x, G.q)
        Ks, As, Bs = G.canon(s)
        e = prf(b"\x02" + Ks.to_bytes(48, "big"), m)
        x2 = G.F.mul(x, G.gtab.pow(e))
        K2, k2 = G.key(x2)
        return (x2, (As + e) % m, Bs, K2, k2)


def walk_cycle(walk, N1, first=None):
    pts, s = [], (walk.step(N1) if first is None else first)
    while True:
        pts.append(s)
        if s[3] == N1[3]:
            return pts
        if len(pts) >= MAXC:
            return None
        s = walk.step(s)


def solve(m, a1, b1, a2, b2):
    return None if (b1 - b2) % m == 0 else (a2 - a1) * pow(b1 - b2, -1, m) % m


def trial(G, rng, kb, max_steps):
    m = G.m
    y = rng.randrange(1, m)
    z = G.gtab.pow(y)
    zt = FixedBase(G.F, z, m)
    walk = Walk(G, zt, rng, kb)
    w0 = G.started
    G.started += 2
    cur = []
    for _ in range(2):
        a, b = rng.randrange(m), rng.randrange(m)
        x = G.F.mul(G.gtab.pow(a), zt.pow(b))
        K, k = G.key(x)
        cur.append((x, a, b, K, k))
    isdp = lambda K: ((K >> AB) & DPM) == 0  # noqa: E731
    prev, anc, tc, rounds = [None, None], [None, None], [0, 0], [0, 0]
    visited, dps, fu, sa, steps = {}, {}, None, None, 0
    for st in cur:
        K, A, B = G.canon(st)
        visited.setdefault(K, (A, B))
    while steps < max_steps:
        for w in range(2):
            W = cur[w]
            N1 = walk.step(W)
            steps += 1
            pv, an = prev[w], anc[w]
            prev[w] = W[3]
            rounds[w] += 1
            due = rounds[w] % 64 == 0
            out = None
            while True:
                K1, A, B = G.canon(N1)
                seen = visited.get(K1)
                if seen is not None:
                    if fu is None and solve(m, A, B, *seen) is not None:
                        fu = steps
                else:
                    visited[K1] = (A, B)
                if isdp(K1):
                    old = dps.get(K1)
                    if old is None:
                        dps[K1] = (A, B)
                    else:
                        yc = solve(m, A, B, *old)
                        if yc is not None:
                            if fu is None:
                                fu = steps
                            if G.gtab.pow(yc) == z:
                                if sa is None:
                                    sa = steps
                                return fu, sa
                    steps += 1
                    tc[w] += 1
                    N1 = walk.restart(N1, prf(b"\x01" + (w0 + w).to_bytes(8, "big") + tc[w].to_bytes(8, "big"), m))
                    prev[w], anc[w], pv, an, due, out = None, None, None, None, False, "restart"
                    continue
                det = None
                if pv == N1[3]:
                    N2 = walk.step(N1)
                    if N2[3] == W[3]:
                        det = walk_cycle(walk, N1, first=N2)
                elif an == N1[3]:
                    det = walk_cycle(walk, N1)
                if det is None:
                    break
                steps += len(det)
                hit = next((P for P in det if isdp(P[3])), None)
                if hit is not None:
                    N1 = hit
                    continue
                E = walk.escape([N1] + det[:-1])
                steps += 1
                prev[w], anc[w], pv, an, out = None, None, None, None, "escape"
                N1 = E
            cur[w] = N1
            if due and out is None:
                anc[w] = N1[3]
    return fu, sa


def wilson(k, n, z=1.96):
    ph = k / n
    den = 1 + z * z / n
    mid = (ph + z * z / (2 * n)) / den
    half = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return [round(max(0.0, mid - half), 4), round(min(1.0, mid + half), 4)]


def main():
    tau, seed, n = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    G = Group(1 << 18, 201)
    N = (G.m - 1) // 12
    s = S / R_A
    V = 1 - s * s - (1 - s) ** 2 / R_A
    T = math.sqrt(2 * N * tau / V)
    rng = random.Random(seed)
    rows = []
    for i in range(n):
        kb = random.Random(seed * 100003 + i).getrandbits(64)
        f, v = trial(G, rng, kb, MAXS)
        rows.append([f, v])
    kw = sum(1 for f, v in rows if v is not None and v <= T)
    fw = sum(1 for f, v in rows if f is not None and f <= T)
    fn = sum(1 for f, v in rows if v is not None)
    print(json.dumps({"tau": tau, "T_steps": round(T, 3), "N": N, "V": round(V, 6), "n": n, "trial_seed": seed,
                      "verified_within_T": kw, "verified_rate": round(kw / n, 4), "verified_wilson95": wilson(kw, n),
                      "first_useful_within_T": fw, "first_useful_rate": round(fw / n, 4),
                      "first_useful_wilson95": wilson(fw, n), "verified_eventually": fn}))


if __name__ == "__main__":
    main()
