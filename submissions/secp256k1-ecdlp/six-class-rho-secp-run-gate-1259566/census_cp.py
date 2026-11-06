"""Toy census of the six-class rho walk with checkpoint-only cycle detection (proof.md section 5).

At each class arrival: a distinguished class ends the trail; at length >= L_max the trail is abandoned; if the
key equals the checkpoint K_cp, the cycle is traversed once more (<= 64 steps, else abandoned) and the walk
doubles the member of its minimum-key class and clears K_cp; at length mod 64 == 0 the key becomes K_cp;
otherwise P <- P + Tab[h][t]. T_c counts all class arrivals of all trails up to the first arrival at a class
already visited by another trail, and is compared with the random-mapping model sqrt(pi*N'/2),
N' = N*R/(R-1). Every chunk of trials draws a fresh table. Stdlib only.

Each detected cycle is classified by the partition indices h of its steps. Index-repeating cycles (every h
used at least twice) are the cycles of the ledger, at most 1/(6R) + sum_{L=3..64} sum_b S2(L,b) R^-(L-b) per
operation. Coincidental cycles (some h used once) close by chance: the next class is one of the trail's last
min(length + 1, 64) classes with probability min(length + 1, 64)/N, which the census sums as their model.

usage: python3 census_cp.py --curve q29 --R 2048 65536 1048576 --trials 64000 --chunk 4000 --jobs 8 --out x.json
"""
import argparse, json, math, random, time
from array import array
from multiprocessing import Pool

TOYS = {"q25": (33556867, 33563419, 1), "q29": (536873251, 536839203, 3)}
B = 7   # y^2 = x^3 + 7


def sqrt_mod(a, q):
    """Square root modulo q = 3 (mod 4), or None."""
    a %= q
    y = pow(a, (q + 1) // 4, q)
    return y if y * y % q == a else None


def repeating_model(R):
    """1/(6R) plus the sum over cycle lengths 3..64 of S2(L, b) R^-(L-b) (associated Stirling numbers)."""
    S2 = [[0] * 33 for _ in range(65)]
    S2[0][0] = 1
    for L in range(2, 65):
        for b in range(1, L // 2 + 1):
            S2[L][b] = b * S2[L - 1][b] + (L - 1) * S2[L - 2][b - 1]
    return 1 / (6 * R) + sum(S2[L][b] / R ** (L - b) for L in range(3, 65) for b in range(1, L // 2 + 1))


class Toy:
    def __init__(self, name, rng):
        q, n, h = TOYS[name]
        self.q, self.r = q, n // h
        while True:
            x = rng.randrange(q)
            y = sqrt_mod(x ** 3 + B, q)
            if y is None:
                continue
            G = self.mul((x, y), h)
            if G is not None and self.mul(G, self.r) is None:
                break
        self.G = G
        beta = next(c for c in (pow(g, (q - 1) // 3, q) for g in range(2, 50)) if c != 1)
        lam = next(c for c in (pow(g, (self.r - 1) // 3, self.r) for g in range(2, 50)) if c != 1)
        for b in (beta, beta * beta % q):
            for l in (lam, lam * lam % self.r):
                if self.mul(G, l) == (G[0] * b % q, G[1]):
                    self.beta, self.lam = b, l
        self.fb = []
        P = G
        for _ in range((self.r.bit_length() + 7) // 8):
            row = [None]
            for _ in range(255):
                row.append(self.add(row[-1], P))
            self.fb.append(row)
            P = self.add(row[-1], P)

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

    def mul(self, P, k):
        R = None
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R

    def fmul(self, k):
        k %= self.r
        R, i = None, 0
        while k:
            R = self.add(R, self.fb[i][k & 255])
            k >>= 8
            i += 1
        return R


def run_chunk(args):
    name, R, dp_bits, trials, seed, p0 = args
    rng = random.Random(seed)
    E = Toy(name, random.Random(12345))
    q, r, beta = E.q, E.r, E.beta
    xs = rng.randrange(1, r)
    rbits = R.bit_length() - 1
    assert 1 << rbits == R
    dpmask = (1 << dp_bits) - 1
    Lmax = 20 << dp_bits
    binv = (1, beta * beta % q, beta)           # beta^-j: phi^-j(x, y) = (beta^-j x, y)
    tx, ty = array("q"), array("q")             # M_h; Tab[h][3s + j] = (-1)^s phi^-j(M_h) on the fly
    while len(tx) < R:
        M = E.fmul(rng.randrange(r) + rng.randrange(r) * xs)
        if M is None:
            continue
        tx.append(M[0])
        ty.append(M[1])
    add = E.add

    def canon(P):
        x, y = P
        bx = beta * x % q
        b2x = (-x - bx) % q
        j, cx = 0, x
        if bx < cx:
            cx, j = bx, 1
        if b2x < cx:
            cx, j = b2x, 2
        ny = q - y
        if ny < y:
            return ny, 1, j
        return y, 0, j

    def target(K, s, j):
        h = K & (R - 1)
        return (tx[h] * binv[j] % q, (q - ty[h]) % q if s else ty[h])

    N = (r - 1) // 6
    Np = N * R / (R - 1)
    S0 = math.sqrt(2 * Np * math.log(1 / (1 - p0)))
    s1 = s2 = 0.0
    within = 0
    doublings = steps = trav_steps = rep = coin = 0
    near = 0                                    # sum over ordinary steps of min(length + 1, 64)
    for _ in range(trials):
        seen = {}
        visits = 0
        trail = 0
        tc = None
        while tc is None:
            P = E.fmul(rng.randrange(r) + rng.randrange(r) * xs)
            if P is None:
                continue
            trail += 1
            cp = None
            length = 0
            while True:
                K, s, j = canon(P)
                visits += 1
                prev = seen.get(K)
                if prev is None:
                    seen[K] = trail
                elif prev != trail:
                    tc = visits
                    break
                if ((K >> rbits) & dpmask) == 0 or length >= Lmax:
                    break
                if cp is not None and K == cp:
                    best = (K, P)
                    Q = P
                    trav = 0
                    hit = False
                    ok = True
                    Kq, sq, jq = K, s, j
                    hs = {}
                    while True:
                        trav += 1
                        if trav > 64:
                            ok = False
                            break
                        hs[Kq & (R - 1)] = hs.get(Kq & (R - 1), 0) + 1
                        TP = target(Kq, sq, jq)
                        if Q[0] == TP[0]:
                            ok = False
                            break
                        Q = add(Q, TP)
                        trav_steps += 1
                        length += 1
                        Kq, sq, jq = canon(Q)
                        visits += 1
                        prev = seen.get(Kq)
                        if prev is None:
                            seen[Kq] = trail
                        elif prev != trail:
                            tc = visits
                            hit = True
                            break
                        if Kq == K:
                            break
                        if Kq < best[0]:
                            best = (Kq, Q)
                    if hit or not ok:
                        break
                    if min(hs.values()) >= 2:
                        rep += 1
                    else:
                        coin += 1
                    P = add(best[1], best[1])
                    doublings += 1
                    if P is None:
                        break
                    cp = None
                    length += 1
                    continue
                if length % 64 == 0:
                    cp = K
                TP = target(K, s, j)
                if P[0] == TP[0]:
                    break
                P = add(P, TP)
                steps += 1
                near += min(length + 1, 64)
                length += 1
        z = tc / math.sqrt(N)
        s1 += z
        s2 += z * z
        within += tc <= S0
    return dict(R=R, seed=seed, trials=trials, s1=s1, s2=s2, within=within, doublings=doublings, steps=steps,
                trav_steps=trav_steps, rep=rep, coin=coin, near=near)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curve", default="q29", choices=sorted(TOYS))
    ap.add_argument("--R", type=int, nargs="+", default=[65536, 1048576])
    ap.add_argument("--dp", type=int, default=5)
    ap.add_argument("--trials", type=int, default=16000)
    ap.add_argument("--chunk", type=int, default=2000)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--seed", type=int, default=31001)
    ap.add_argument("--p0", type=float, default=0.3950)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    t0 = time.time()
    E = Toy(a.curve, random.Random(12345))
    r = E.r
    N = (r - 1) // 6
    tasks = [(a.curve, R, a.dp, a.chunk, a.seed + 1000003 * R.bit_length() + c, a.p0)
             for R in a.R for c in range(a.trials // a.chunk)]
    with Pool(a.jobs) as pool:
        out = pool.map(run_chunk, tasks, chunksize=1)
    res = dict(walk="checkpoint-only", curve=a.curve, q=E.q, r=r, N=N, dp_bits=a.dp, chunk=a.chunk, p0=a.p0,
               per_R={})
    for R in a.R:
        rows = [o for o in out if o["R"] == R]
        n = sum(o["trials"] for o in rows)
        mean = sum(o["s1"] for o in rows) / n
        var = sum(o["s2"] for o in rows) / n - mean * mean
        stderr = math.sqrt(var / n)

        def spread(xs):
            m = sum(xs) / len(xs)
            return math.sqrt(sum((x - m) ** 2 for x in xs) / max(1, len(xs) - 1) / len(xs))
        Np = N * R / (R - 1)
        model = math.sqrt(math.pi * Np / 2) / math.sqrt(N)
        opr = [o["steps"] + o["trav_steps"] + o["doublings"] for o in rows]
        ops = sum(opr)
        res["per_R"][str(R)] = dict(
            trials=n, tables=len(rows), mean_tc_over_sqrtN=mean, stderr=stderr, model_mean_over_sqrtN=model,
            slowdown_vs_model=mean / model - 1, slowdown_stderr=stderr / model,
            stderr_from_table_spread=spread([o["s1"] / o["trials"] for o in rows]) / model,
            frac_tc_le_S0=sum(o["within"] for o in rows) / n,
            operations=ops,
            doublings_per_op=sum(o["doublings"] for o in rows) / ops,
            traversal_steps_per_op=sum(o["trav_steps"] for o in rows) / ops,
            index_repeating_cycles_per_op=sum(o["rep"] for o in rows) / ops,
            index_repeating_model=repeating_model(R),
            coincidental_cycles_per_op=sum(o["coin"] for o in rows) / ops,
            coincidental_stderr_from_table_spread=spread([o["coin"] / x for o, x in zip(rows, opr)]),
            coincidental_model=sum(o["near"] for o in rows) / N / ops)
    res["seconds"] = round(time.time() - t0, 1)
    s = json.dumps(res, indent=1)
    print(s)
    if a.out:
        open(a.out, "w").write(s + "\n")


if __name__ == "__main__":
    main()
