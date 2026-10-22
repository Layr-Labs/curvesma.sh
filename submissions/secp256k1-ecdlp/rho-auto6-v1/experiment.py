"""Toy-scale check of truncated Pollard rho on equivalence classes of the
order-6 automorphism group <-1, phi> of a j=0 curve y^2 = x^3 + 7 (the same
curve family and endomorphism as secp256k1), using the r-adding walk,
canonical class representatives and fruitless-cycle escape described in
proof.md.  Python standard library only.  Evidence at toy scope only; it does
not establish production cost."""
import json, math, random, sys

def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def sqrt_mod(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1: tt = tt * tt % p; i += 1
        b = pow(c, 1 << (m - i - 1), p); m = i; c = b * b % p; t = t * c % p; r = r * b % p
    return r

def add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        l = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
    else:
        l = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)

def mul(k, P, p):
    R = None
    while k:
        if k & 1: R = add(R, P, p)
        P = add(P, P, p); k >>= 1
    return R

def find_curve(bits, rng):
    while True:
        p = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if p % 3 != 1 or not is_prime(p): continue
        r0 = sqrt_mod(p - 3, p)
        a, b, lim = p, r0, math.isqrt(p)
        while b > lim: a, b = b, a % b
        x = b; y2 = (p - x * x) // 3; y = math.isqrt(y2)
        if y * y != y2: continue
        while True:
            X = rng.randrange(p); Y = sqrt_mod(X ** 3 + 7, p)
            if Y: break
        G = (X, Y)
        for t in (2 * x, -2 * x, x + 3 * y, -(x + 3 * y), x - 3 * y, -(x - 3 * y)):
            n = p + 1 - t
            if is_prime(n) and mul(n, G, p) is None:
                return p, n, G

def run(bits, trials, K, seed):
    rng = random.Random(seed)
    p, n, P = find_curve(bits, rng)
    # beta: primitive cube root of unity mod p; lam: matching one mod n
    g = 2
    while pow(g, (p - 1) // 3, p) == 1: g += 1
    beta = pow(g, (p - 1) // 3, p)
    h = 2
    while pow(h, (n - 1) // 3, n) == 1: h += 1
    lam = pow(h, (n - 1) // 3, n)
    if mul(lam, P, p) != (beta * P[0] % p, P[1]):
        lam = lam * lam % n
    assert mul(lam, P, p) == (beta * P[0] % p, P[1])
    b2 = beta * beta % p; l2 = lam * lam % n
    def canon(R, a, b):
        # orbit {+-phi^k R}; choose min x, then y <= p-1-y
        x, y = R
        cands = ((x, 1), (beta * x % p, lam), (b2 * x % p, l2))
        xm, f = min(cands)
        if y > p - y: y = p - y; f = (-f) % n
        return (xm, y), a * f % n, b * f % n
    Bp = math.sqrt(2 * (n / 6) * math.log(1 / 0.61))
    # budget rule used in proof.md (before the separate cycle/restart factors)
    Bprod = 1.05 * math.sqrt(2 * (n / 6) * math.log(1 / 0.58))
    steps_list, fruitless_total, ok, okprod = [], 0, 0, 0
    for _ in range(trials):
        x = rng.randrange(1, n); Q = mul(x, P, p)
        table = []
        for _ in range(K):
            c, d = rng.randrange(n), rng.randrange(n)
            table.append((add(mul(c, P, p), mul(d, Q, p), p), c, d))
        a0, b0 = rng.randrange(n), rng.randrange(1, n)
        R, a, b = canon(add(mul(a0, P, p), mul(b0, Q, p), p), a0, b0)
        seen = {}; steps = 0; prev_j = -1
        while True:
            key = R
            if key in seen:
                a2, b2c = seen[key]
                if (b - b2c) % n:
                    xr = (a2 - a) * pow(b - b2c, -1, n) % n
                    assert xr == x
                    break
                # fruitless cycle: escape by doubling
                fruitless_total += 1
                R, a, b = canon(add(R, R, p), 2 * a, 2 * b); steps += 1
                continue
            seen[key] = (a, b)
            j = R[0] % K
            if j == prev_j: j = (j + 1) % K  # 2-cycle reduction
            M, c, d = table[j]; prev_j = j
            S = add(R, M, p)
            if S is None:
                a0, b0 = rng.randrange(n), rng.randrange(1, n)
                R, a, b = canon(add(mul(a0, P, p), mul(b0, Q, p), p), a0, b0); steps += 1; continue
            R, a, b = canon(S, a + c, b + d); steps += 1
        steps_list.append(steps)
        if steps <= Bp: ok += 1
        if steps <= Bprod: okprod += 1
    mean = sum(steps_list) / trials
    return {"p": p, "n": n, "K": K, "trials": trials, "seed": seed,
            "budget_B_for_0.39": round(Bp, 1),
            "empirical_success_within_B": round(ok / trials, 4),
            "predicted_success_within_B": 0.39,
            "budget_rule_proof_md": round(Bprod, 1),
            "empirical_success_within_budget_rule": round(okprod / trials, 4),
            "mean_steps": round(mean, 1),
            "predicted_mean_sqrt_pi_n_over_12": round(math.sqrt(math.pi * n / 12), 1),
            "mean_ratio": round(mean / math.sqrt(math.pi * n / 12), 4),
            "fruitless_escapes_total": fruitless_total,
            "scope": "toy j=0 curves only; supports the class-walk collision model, not production cost"}

if __name__ == "__main__":
    args = json.load(sys.stdin)
    print(json.dumps(run(args["bits"], args["trials"], args.get("partitions", 64), args.get("seed", 1))))
