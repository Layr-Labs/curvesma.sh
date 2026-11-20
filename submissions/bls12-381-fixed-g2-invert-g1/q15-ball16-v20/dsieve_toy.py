"""Descent sieve hit enumeration (toy): points of the index-p sublattice L_p = {v in Z^16 : c.v = 0 mod p} in the
ball |v| <= R, by Schnorr-Euchner enumeration on an LLL-reduced basis of L_p; counts tree nodes per hit and checks
the hit set against brute force on small radii."""
import math, random, sys


def lll(B, delta=0.99):
    B = [list(map(int, r)) for r in B]; n = len(B)
    def gso(B):
        Bs = []; mu = [[0.0] * n for _ in range(n)]; nb = []
        for i in range(n):
            v = [float(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = sum(a * b for a, b in zip(B[i], Bs[j])) / nb[j]
                v = [a - mu[i][j] * b for a, b in zip(v, Bs[j])]
            Bs.append(v); nb.append(sum(a * a for a in v))
        return Bs, mu, nb
    k = 1; Bs, mu, nb = gso(B)
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [a - q * b for a, b in zip(B[k], B[j])]; Bs, mu, nb = gso(B)
        if nb[k] >= (delta - mu[k][k - 1] ** 2) * nb[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu, nb = gso(B); k = max(k - 1, 1)
    return B


def sublattice_basis(c, p):
    """rows: p e_0 and e_i - (c_i/c_0) e_0 (c_0 invertible mod p)."""
    n = len(c); inv = pow(c[0], -1, p); B = [[p] + [0] * (n - 1)]
    for i in range(1, n):
        r = [0] * n; r[0] = (-c[i] * inv) % p; r[i] = 1; B.append(r)
    return B


def enumerate_ball(B, R2):
    """Schnorr-Euchner over all lattice points v = x B with |v|^2 <= R2; returns (points, nodes)."""
    n = len(B)
    Bs = []; mu = [[0.0] * n for _ in range(n)]; nb = []
    for i in range(n):
        v = [float(x) for x in B[i]]
        for j in range(i):
            mu[i][j] = sum(a * b for a, b in zip(B[i], Bs[j])) / nb[j]
            v = [a - mu[i][j] * b for a, b in zip(v, Bs[j])]
        Bs.append(v); nb.append(sum(a * a for a in v))
    x = [0] * n; pts = []; nodes = [0]
    def rec(k, rem):
        c = -sum(mu[j][k] * x[j] for j in range(k + 1, n))
        r = math.sqrt(max(rem, 0.0) / nb[k])
        for t in range(math.ceil(c - r - 1e-12), math.floor(c + r + 1e-12) + 1):
            nodes[0] += 1
            d = (t - c) ** 2 * nb[k]
            if d <= rem + 1e-9:
                x[k] = t
                if k == 0:
                    if any(x):
                        pts.append(tuple(sum(x[i] * B[i][j] for i in range(n)) for j in range(n)))
                else:
                    rec(k - 1, rem - d)
        x[k] = 0
    rec(n - 1, float(R2))
    return pts, nodes[0]


def brute(c, p, R):
    n = len(c); out = []; v = [0] * n; R2 = R * R
    def rec(i, rem):
        if i == n:
            if any(v) and sum(a * b for a, b in zip(c, v)) % p == 0: out.append(tuple(v))
            return
        r = math.isqrt(rem)
        for t in range(-r, r + 1):
            v[i] = t; rec(i + 1, rem - t * t)
        v[i] = 0
    rec(0, R2)
    return out


if __name__ == "__main__":
    rng = random.Random(7)
    # correctness against brute force (dimension 16, small radius)
    for p in (101, 1009):
        c = [rng.randrange(1, p) for _ in range(16)]
        pts, nodes = enumerate_ball(lll(sublattice_basis(c, p)), 4)
        bf = brute(c, p, 2)
        print("check p", p, "hits", len(pts), "brute", len(bf), "equal", sorted(pts) == sorted(bf), "nodes/hit", round(nodes / max(1, len(pts)), 2))
    # nodes per hit as the expected hit count varies
    for R, p in ((3, 1009), (4, 10007), (4, 100003), (5, 1000003), (6, 10000019)):
        c = [rng.randrange(1, p) for _ in range(16)]
        pts, nodes = enumerate_ball(lll(sublattice_basis(c, p)), R * R)
        print("R", R, "p", p, "hits", len(pts), "nodes", nodes, "nodes/hit", round(nodes / max(1, len(pts)), 2))
        sys.stdout.flush()
