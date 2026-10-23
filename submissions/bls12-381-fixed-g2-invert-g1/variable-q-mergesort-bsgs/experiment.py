"""Toy mergesort BSGS. Not a production cost or curve-break check."""
import json
import math
import sys


def prime_factors(n):
    fac = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            fac.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        fac.append(n)
    return fac


def is_prime(n):
    if n < 2:
        return False
    return prime_factors(n) == [n]


def is_primitive(g, q):
    order = q - 1
    return all(pow(g, order // f, q) != 1 for f in prime_factors(order))


def giant_inverse(h, order, mod):
    """h^order = 1 implies h^(order-1) is the inverse. order-2 is not."""
    if pow(h, order, mod) != 1:
        raise AssertionError("element not in the subgroup")
    inv = pow(h, order - 1, mod)
    if (inv * h) % mod != 1:
        raise AssertionError("r-1 did not invert")
    if order > 2 and (pow(h, order - 2, mod) * h) % mod == 1:
        raise AssertionError("r-2 inverted; group law is wrong")
    return inv


def recover(mod, order, g, secret):
    m = math.isqrt(order) + 1
    step = pow(g, m, mod)
    giant = giant_inverse(step, order, mod)
    babies = []
    acc = 1
    for j in range(m):
        babies.append((acc, j))
        acc = (acc * g) % mod
    babies.sort()
    keys = [row[0] for row in babies]
    target = pow(g, secret, mod)
    hops = (order + m - 1) // m
    for i in range(hops):
        pos = bisect_left(keys, target)
        if pos < len(keys) and keys[pos] == target:
            return (i * m + babies[pos][1]) % order
        target = (target * giant) % mod
    raise AssertionError("discrete log not recovered")


def bisect_left(keys, target):
    lo = 0
    hi = len(keys)
    while lo < hi:
        mid = (lo + hi) // 2
        if keys[mid] < target:
            lo = mid + 1
        else:
            hi = mid
    return lo


def witness_29():
    mod, order, g = 233, 29, 2
    if pow(g, order, mod) != 1:
        raise AssertionError("witness order")
    for secret in range(order):
        if recover(mod, order, g, secret) != secret:
            raise AssertionError("witness mismatch")
    return order


def main():
    spec = json.load(sys.stdin)
    q = spec["modulus"]
    g = spec["generator"]
    if not is_prime(q):
        raise SystemExit("modulus must be prime")
    if g <= 1 or g >= q or not is_primitive(g, q):
        raise SystemExit("generator must be a primitive root")
    order = q - 1
    checked = 0
    for secret in range(order):
        if recover(q, order, g, secret) != secret:
            raise SystemExit("recovery mismatch")
        checked += 1
    witness = witness_29()
    json.dump(
        {
            "group": "multiplicative group of a prime field",
            "modulus": q,
            "order": order,
            "recovered": checked,
            "inverse_exponent": "order-1",
            "witness_order": witness,
            "witness_rejects_exponent": "order-2",
            "scope": "toy groups only; no production cost",
        },
        sys.stdout,
    )


if __name__ == "__main__":
    main()
