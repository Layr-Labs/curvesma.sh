"""Toy six-class checkpoint walk, the linear unit identity, and phi(P0).

The toy group is Z/nZ, n prime, n ≡ 1 (mod 6), with multiplication by a
root of X^2+X+1 as the order-6 action. A checkpoint collision must return
the planted scalar. That run is not a production instruction count.

The production checks are integer identities: the six unit maps on (a, b)
agree with multiplication by ±λ^k, β^2 + β + 1 = 0 in Fp, and
phi(P0) = [λ]P0. proof.md extends the generator identity by the
homomorphism property.
"""

from __future__ import annotations

import json
import random
import sys

import reference

P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
BETA = 0x1A0111EA397FE699EC02408663D4DE85AA0D857D89759AD4897D29650FB85F9B409427EB4F49FFFD8BFD00000000AAAC
P0 = (
    0x17F1D3A73197D7942695638C4FA9AC0FC3688C4F9774B905A14E3A3F171BAC586C55E83FF97A1AEFFB3AF00ADB22C6BB,
    0x08B3F481E3AAA0F1A09E30ED741D8AE4FCF5E095D5D00AF600DB18CB2C04B3EDD03CC744A2888AE40CAA232946C5E7E1,
)


def _add(p, q):
    if p is None:
        return q
    if q is None:
        return p
    x1, y1 = p
    x2, y2 = q
    if x1 == x2 and (y1 + y2) % P == 0:
        return None
    if p != q:
        slope = ((y2 - y1) * pow(x2 - x1, -1, P)) % P
    else:
        slope = ((3 * x1 * x1) * pow(2 * y1, -1, P)) % P
    x3 = (slope * slope - x1 - x2) % P
    y3 = (slope * (x1 - x3) - y1) % P
    return x3, y3


def _mul(k, point):
    acc = None
    base = point
    while k:
        if k & 1:
            acc = _add(acc, base)
        base = _add(base, base)
        k >>= 1
    return acc


def generator_identity():
    lam = reference.LAMBDA
    if lam * lam + lam + 1 != reference.GROUP_ORDER:
        raise AssertionError("lambda polynomial")
    if (BETA * BETA + BETA + 1) % P != 0:
        raise AssertionError("beta polynomial")
    # β² x = -x - β x, checked on the generator x-coordinate.
    x = P0[0]
    beta_x = (BETA * x) % P
    beta2_x = (BETA * beta_x) % P
    if beta2_x != (-x - beta_x) % P:
        raise AssertionError("beta^2 x identity")
    image = (beta_x, P0[1])
    if _mul(lam, P0) != image:
        raise AssertionError("phi(P0) != [lambda]P0")
    if not reference.unit_images_match_multiplication():
        raise AssertionError("linear unit map")
    return True


def _miller(n: int, rounds: int = 8) -> bool:
    if n < 3 or n % 2 == 0:
        return False
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19)[:rounds]:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def toy_prime(bits: int) -> int:
    n = (1 << (bits - 1)) + 1
    n += (1 - n) % 6
    while n.bit_length() == bits:
        if _miller(n):
            return n
        n += 6
    raise RuntimeError("no prime")


def tonelli(n: int, p: int) -> int:
    n %= p
    if pow(n, (p - 1) // 2, p) != 1:
        raise ValueError("non-residue")
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)
    while t != 1:
        i = 1
        t2 = (t * t) % p
        while t2 != 1:
            i += 1
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r


def toy_lambda(n: int) -> int:
    root = tonelli(-3 % n, n)
    lam = ((-1 + root) * pow(2, -1, n)) % n
    if (lam * lam + lam + 1) % n != 0 or pow(lam, 3, n) != 1 or lam == 1:
        raise AssertionError("toy lambda")
    other = (-1 - lam) % n
    return lam if lam.bit_length() <= other.bit_length() else other


def orbit(point: int, n: int, lam: int) -> set[int]:
    lam2 = (-1 - lam) % n
    return {(m * point) % n for m in (1, n - 1, lam, (-lam) % n, lam2, (-lam2) % n)}


def canonical(point: int, n: int, lam: int) -> tuple[int, int]:
    lam2 = (-1 - lam) % n
    best = n
    best_m = 1
    for m in (1, n - 1, lam, (-lam) % n, lam2, (-lam2) % n):
        img = (m * point) % n
        if img < best:
            best = img
            best_m = m
    return best, best_m


def solve_toy(n: int, lam: int, secret: int, trials_seed: int) -> dict:
    rng = random.Random(trials_seed)
    n_jumps = 32
    table = [(rng.randrange(1, n), rng.randrange(n)) for _ in range(n_jumps)]
    theta = 8
    classes = (n - 1) // 6
    budget = reference.sigma_upper_int(classes) + reference.WALK_SLACK + 8 * theta
    store: dict[int, tuple[int, int]] = {}
    steps = 0
    restarts = 0
    while steps < budget and restarts < 64:
        a = rng.randrange(n)
        b = rng.randrange(n)
        if a == 0 and b == 0:
            b = 1
        restarts += 1
        escaped = False
        while steps < budget and not escaped:
            steps += 1
            point = (a + b * secret) % n
            if point == 0:
                if b % n == 0:
                    escaped = True
                    break
                return {"ok": True, "steps": steps, "restarts": restarts, "budget": budget}
            can, m = canonical(point, n, lam)
            a = (a * m) % n
            b = (b * m) % n
            if can % theta == 0:
                prev = store.get(can)
                if prev is None:
                    store[can] = (a, b)
                else:
                    a0, b0 = prev
                    db = (b - b0) % n
                    if db == 0:
                        escaped = True
                        break
                    got = ((a0 - a) * pow(db, -1, n)) % n
                    if got != secret % n:
                        raise AssertionError("wrong scalar")
                    return {"ok": True, "steps": steps, "restarts": restarts, "budget": budget}
            du, dv = table[can % n_jumps]
            a = (a + du) % n
            b = (b + dv) % n
    return {"ok": False, "steps": steps, "restarts": restarts, "budget": budget}


def run(bits: int, trials: int) -> dict:
    n = toy_prime(bits)
    lam = toy_lambda(n)
    if len(orbit(2, n, lam)) != 6:
        raise AssertionError("orbit size")
    hits = 0
    step_sum = 0
    for i in range(trials):
        secret = random.Random(1000 + bits * 100 + i).randrange(n)
        out = solve_toy(n, lam, secret, trials_seed=50_000 + i)
        hits += int(out["ok"])
        step_sum += out["steps"]
    return {
        "bits": bits,
        "prime": n,
        "classes": (n - 1) // 6,
        "trials": trials,
        "recoveries": hits,
        "mean_steps": step_sum / trials,
        "sigma_upper": reference.sigma_upper_int((n - 1) // 6),
    }


def main() -> None:
    data = json.load(sys.stdin)
    bits = data["bits"]
    trials = data["trials"]
    if (
        isinstance(bits, bool)
        or isinstance(trials, bool)
        or not isinstance(bits, int)
        or not isinstance(trials, int)
        or not (12 <= bits <= 16)
        or not (1 <= trials <= 8)
    ):
        raise SystemExit("input out of range")
    reference.self_check()
    identity = generator_identity()
    toy = run(bits, trials)
    if toy["recoveries"] < 1:
        raise SystemExit("toy walk recovered nothing")
    print(json.dumps({
        "generator_endomorphism_identity": identity,
        "linear_unit_action": True,
        "toy": toy,
        "scope": "toy checkpoint recoveries, linear unit maps, one generator identity; not a production cost measurement",
    }))


if __name__ == "__main__":
    main()
