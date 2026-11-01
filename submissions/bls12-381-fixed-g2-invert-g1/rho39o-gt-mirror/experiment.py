"""Toy-scale validation of the GT parallel distinguished-points rho attack.

Same walk mechanics as the curve version but in multiplicative notation
over a prime-order subgroup of F_p^* (mirroring GT: partitions and
distinguished tests read the element representation, identity is never
distinguished). Runs P = 8 walks with a shared table, collision solve,
independent verification. Mechanism check at toy scale only; it does
not execute or prove production cost bounds.
"""
import json
import math
import random
import sys

R_PARTITIONS = 64
WALKS = 8
DISTINGUISHED_MOD = 4  # toy density 1/4 (production uses 2^-32)
TAIL_FACTOR = 1  # minimal tails: mechanism check only
OFFSET = 2**40  # toy per-walk coefficient offset (production: 2**151)
TRIALS = 2000
MIN_RATE = 0.34  # nominal prediction ~0.40; margin for noise


def find_group(n):
    """Find prime p = k*n+1 and generator g of the order-n subgroup."""
    for k in range(1, 200):
        p = k * n + 1
        if p < 5:
            continue
        if all(p % d for d in range(2, int(p**0.5) + 1)):
            for h in range(2, p):
                g = pow(h, k, p)
                if g != 1 and pow(g, n, p) == 1:
                    return p, g
    raise AssertionError(f"no toy group for order {n}")


def walk_length(n):
    # Same formula as production: smallest T with T(T-1) >= 2*n*0.5113.
    target = 2 * n * 0.5113
    t = int(math.sqrt(target)) + 1
    while t * (t - 1) < target:
        t += 1
    return t


def partition(Y):
    return Y % R_PARTITIONS


def distinguished(Y):
    # Window disjoint from the partition bits (mirrors production v15).
    return Y != 1 and (Y >> 6) % DISTINGUISHED_MOD == 0


DEGENS = [0]
DBS = []

def rho_attack(n, g, z, p, rng):
    t_walk = walk_length(n)
    per_walk = t_walk // WALKS + 1
    tail = TAIL_FACTOR * DISTINGUISHED_MOD
    table = {}
    u = [rng.randrange(1, 2**32) for _ in range(R_PARTITIONS)]
    v = [rng.randrange(1, 2**32) for _ in range(R_PARTITIONS)]
    M = [pow(g, u[j], p) * pow(z, v[j], p) % p
         for j in range(R_PARTITIONS)]
    for w in range(WALKS):
        ca, cb = w * OFFSET + rng.randrange(2**32), rng.randrange(2**32)
        Y = pow(g, ca, p) * pow(z, cb, p) % p
        for _ in range(per_walk + tail):
            f = partition(Y)
            ca = ca + u[f]
            cb = cb + v[f]
            Y = Y * M[f] % p
            if distinguished(Y):
                if Y in table:
                    a0, b0 = table[Y]
                    da = (ca - a0) % n
                    db = (cb - b0) % n
                    DBS.append(db)
                    if db == 0:
                        DEGENS[0] += 1
                    if db != 0:
                        rec = (-da * pow(db, -1, n)) % n
                        if pow(g, rec, p) == z:
                            return rec
                        raise AssertionError("wrong scalar from solve")
                else:
                    table[Y] = (ca, cb)
    return None


def main():
    order = json.load(sys.stdin)["order"]
    assert isinstance(order, int) and 2 <= order <= 10000
    p, g = find_group(order)
    assert pow(g, order, p) == 1 and g != 1
    rng = random.Random(0xC0FFEE + order)
    t_walk = walk_length(order)
    pred = 1.0 - math.exp(-t_walk * (t_walk - 1) / (2 * order))
    ok = 0
    for _ in range(TRIALS):
        secret = rng.randrange(order)
        z = pow(g, secret, p)
        if rho_attack(order, g, z, p, rng) is not None:
            ok += 1
    rate = ok / TRIALS
    assert rate >= MIN_RATE, f"success rate {rate:.3f} below {MIN_RATE}"
    BINS = 25
    bins = [0] * BINS
    for db in DBS:
        bins[(db * BINS) // order] += 1
    exp_b = len(DBS) / BINS if DBS else 1
    chi_b = sum((h - exp_b) ** 2 / exp_b for h in bins) if DBS else 0
    assert chi_b < 40, f"db chi-square {chi_b:.2f} above 40 (uniformity break)"
    print(json.dumps({
        "group": f"order-{order} subgroup of F_{p}^*",
        "order": order,
        "walks": WALKS,
        "trials": TRIALS,
        "successes": ok,
        "success_rate": round(rate, 4),
        "predicted_walk_collision": round(pred, 4),
        "t_walk": t_walk,
        "all_recovered_correct": True,
        "degenerate_encounters": DEGENS[0],
        "db_solves": len(DBS),
        "db_chi2_24df": round(chi_b, 2),
        "scope": "toy group only; does not validate production cost bounds",
    }))


main()
