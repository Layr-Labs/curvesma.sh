"""Toy end-to-end run of the class-walk scheme: prime-order subgroup of Z_q^* with an order-12 exponent
automorphism y -> y^lam, state kept un-canonicalised, conjugated jump tables, distinguished points.
Only a small-group sanity check; it says nothing about production costs."""
import json
import random
import sys


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def setup(order):
    r = order
    while not (r % 12 == 1 and is_prime(r)):
        r -= 1
    m = 2
    while not is_prime(m * r + 1):
        m += 1
    q = m * r + 1
    g = next(pow(h, m, q) for h in range(2, q) if pow(h, m, q) != 1)
    lam = next(pow(h, (r - 1) // 12, r) for h in range(2, r) if len({pow(pow(h, (r - 1) // 12, r), i, r) for i in range(12)}) == 12)
    return r, q, g, lam


def solve(r, q, g, lam, z, rng, tab=64, dpbits=3):
    H = [pow(lam, i, r) for i in range(12)]
    Hinv = [pow(h, -1, r) for h in H]
    jumps = [(rng.randrange(r), rng.randrange(r)) for _ in range(tab)]
    M = [pow(g, a, q) * pow(z, b, q) % q for a, b in jumps]
    Mc = [[pow(Mj, hi, q) for hi in Hinv] for Mj in M]
    ec = [[(a * hi % r, b * hi % r) for hi in Hinv] for a, b in jumps]
    classes = (r - 1) // 12
    budget = int(3.0 * classes ** 0.5) + 200
    seen, steps = {}, 0
    while steps < budget:
        al, be = rng.randrange(r), rng.randrange(r)
        y = pow(g, al, q) * pow(z, be, q) % q
        for _ in range(64 * (1 << dpbits)):
            imgs = [pow(y, h, q) for h in H]
            idx = min(range(12), key=imgs.__getitem__)
            key = imgs[idx]
            j = (key * 2654435761 >> 7) % tab
            steps += 1
            if (key * 40503 >> 5) % (1 << dpbits) == 0:
                ca, cb = al * H[idx] % r, be * H[idx] % r
                if key in seen:
                    a2, b2 = seen[key]
                    if b2 != cb:
                        return (a2 - ca) * pow(cb - b2, -1, r) % r, steps
                else:
                    seen[key] = (ca, cb)
                break
            y = y * Mc[j][idx] % q
            al, be = (al + ec[j][idx][0]) % r, (be + ec[j][idx][1]) % r
    return None, steps


def main():
    spec = json.load(sys.stdin)
    r, q, g, lam = setup(spec["order"])
    rng = random.Random(spec.get("seed", 1))
    ok = bad = 0
    for _ in range(spec.get("trials", 20)):
        x = rng.randrange(1, r)
        got, _ = solve(r, q, g, lam, pow(g, x, q), rng)
        if got is None:
            continue
        ok += got == x
        bad += got != x
    print(json.dumps({"subgroup_order": r, "modulus": q, "automorphism_order": 12,
                      "recovered": ok, "wrong": bad, "scope": "toy group only"}))
    assert bad == 0


main()
