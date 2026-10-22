"""Toy end-to-end run of the class-rho method: group Z_n (n prime, n = 1 mod 12) with the order-12 automorphism
u -> lam*u standing in for z -> z^p. Parallel walks, distinguished points, 2-cycle escape. Toy only; it validates the
logic (class canonicalisation, exponent recovery, success probability), not production costs."""
import json
import math
import random
import sys


def prime(n):
    return n > 1 and all(n % q for q in range(2, math.isqrt(n) + 1))


def run(n, trials):
    lam = next(pow(g, (n - 1) // 12, n) for g in range(2, n) if pow(pow(g, (n - 1) // 12, n), 6, n) != 1
               and pow(pow(g, (n - 1) // 12, n), 4, n) != 1)           # element of order exactly 12
    pw = [pow(lam, i, n) for i in range(12)]
    classes = (n - 1) // 12 + 1
    k = math.sqrt(-2 * math.log(0.61) * classes)
    R, theta, walks = 64, 4.0 / math.sqrt(classes), 6
    ok = 0
    steps_used = 0
    for t in range(trials):
        rnd = random.Random(1000 + t)
        x = rnd.randrange(1, n)
        table = [rnd.randrange(1, n) for _ in range(R)]          # step j multiplies by g^table[j]: u += table[j]
        dps = {}
        budget = 1.1 * k + walks * 1.0 / theta
        used = 0
        found = None
        W = [[rnd.randrange(n), rnd.randrange(n)] for _ in range(walks)]   # exponents (a,b): u = a + b*x
        hist = [[] for _ in range(walks)]
        length = [0] * walks
        while used < budget and found is None:
            for w in range(walks):
                a, b = W[w]
                u = (a + b * x) % n
                i = min(range(12), key=lambda s: pw[s] * u % n)       # canonical conjugate
                a, b = pw[i] * a % n, pw[i] * b % n
                u = (a + b * x) % n
                j = (u * 2654435761 >> 8) % R
                if u in hist[w][-3:]:
                    j = (j + 1) % R                                    # escape a fruitless short cycle
                hist[w] = (hist[w] + [u])[-4:]
                length[w] += 1
                used += 1
                if u % int(1 / theta) == 0 or length[w] > 20 / theta:  # distinguished point (class invariant)
                    if u in dps and dps[u][1] != b:
                        a2, b2 = dps[u]
                        found = (a - a2) * pow(b2 - b, -1, n) % n
                        break
                    dps[u] = (a, b)
                    W[w] = [rnd.randrange(n), rnd.randrange(n)]
                    hist[w], length[w] = [], 0
                else:
                    W[w] = [(a + table[j]) % n, b]
        ok += found == x
        steps_used += used
    return {"order": n, "classes": classes, "trials": trials, "recovered": ok, "success_rate": ok / trials,
            "mean_steps_over_sqrt_classes": steps_used / trials / math.sqrt(classes)}


if __name__ == "__main__":
    n = json.load(sys.stdin)["order"]
    n = next(m for m in range(n, 2 * n) if m % 12 == 1 and prime(m))
    res = run(n, 400)
    assert res["success_rate"] >= 0.39, res
    res["scope"] = "toy group only; does not validate production cost bounds"
    print(json.dumps(res))
