"""Toy Pollard rho on a prime-order subgroup of Z_n^*.

Exercises the walk, sector mixing, collision detection, log-vector extraction
and verification that proof.md describes. This is the naive full-table variant
(the production attack stores only 2^-16 distinguished points for memory; the
toy stores every visited state so a collision is guaranteed even for the tiny
subgroups that occur for small prime `order` inputs). Toy scope only: small
prime-order subgroups, not the production secp256k1 curve; it does not
validate production cost bounds, memory bounds or reference-operation prices.

Input: JSON on stdin, {"order": n} with 101 <= n <= 100000 (n prime).
"""
import json
import sys


def _is_prime(x: int) -> bool:
    if x < 2:
        return False
    d = 2
    while d * d <= x:
        if x % d == 0:
            return False
        d += 1
    return True


def _factorize(m: int):
    """Small trial-division factorization; m <= 100000 here."""
    factors, d, x = {}, 2, m
    while d * d <= x:
        while x % d == 0:
            factors[d] = factors.get(d, 0) + 1
            x //= d
        d += 1
    if x > 1:
        factors[x] = factors.get(x, 0) + 1
    return factors


def largest_prime_factor(m: int) -> int:
    return max(_factorize(m)) if m > 1 else 1


def pollard_rho(G: int, H: int, n: int, q: int, max_steps: int = 200000):
    """Recover x in [0,q) with H == G^x mod n via Pollard rho (full table).

    <G> is a prime-order-q subgroup (q prime), so the log-vector inverse below
    is well-defined whenever the collision denominator is nonzero mod q.

    State R = G^a H^b mod n. Each step multiplies R by G, H or G*H based on the
    sector R % 3, so (a, b) track the log vector. On a collision R1 == R2 with
    log vectors (a1,b1),(a2,b2):
        G^(a1-a2) = H^(b2-b1) = G^(x(b2-b1))
        =>  x = (a1-a2) * (b2-b1)^(-1)  (mod q).
    A candidate x is accepted only after verification G^x == H.
    """
    G_H = (G * H) % n
    # Full table (store every visited state) so a collision is guaranteed within
    # q+1 steps, even for the tiny subgroups that occur for small prime `order`
    # inputs. The production attack keeps only 2^-16 distinguished points for
    # memory; that optimization is omitted here.
    table = {}  # r -> (a, b)
    a = b = 0
    r = 1 % n
    for _ in range(max_steps):
        s = r % 3
        if s == 0:
            r = (r * G) % n
            a = (a + 1) % q
        elif s == 1:
            r = (r * H) % n
            b = (b + 1) % q
        else:
            r = (r * G_H) % n
            a = (a + 1) % q
            b = (b + 1) % q
        if r in table:
            a1, b1 = table[r]
            denom = (b - b1) % q
            if denom != 0:
                x = ((a1 - a) % q * pow(denom, -1, q)) % q
                if pow(G, x, n) == H % n:
                    return x
            # denom == 0 or verification failed: false collision, keep walking.
        else:
            table[r] = (a, b)
    raise RuntimeError("rho did not collide within max_steps")


def find_order_q_element(n: int, q: int) -> int:
    """Return an element of Z_n^* of order exactly q (q prime).

    For base a, G = a^((n-1)//q) mod n has G^q = a^(n-1) = 1, so its order
    divides the prime q; G != 1 implies order q. Try small bases until one
    works (2 is not always a suitable base).
    """
    m = n - 1
    for a in range(2, 200):
        G = pow(a, m // q, n)
        if G != 1:
            return G
    raise RuntimeError("no small base gives an order-q element")


def main():
    inp = json.load(sys.stdin)
    n = int(inp["order"])
    if n < 101 or n > 100000 or not _is_prime(n):
        raise ValueError("order must be a prime in [101, 100000]")
    q = largest_prime_factor(n - 1)  # prime-order subgroup of Z_n^*
    G = find_order_q_element(n, q)
    x0 = q // 2  # deterministic secret in [0, q)
    H = pow(G, x0, n)
    if q < 64:
        # Tiny subgroup (e.g. Fermat/Mersenne primes where n-1 is a 2-power,
        # or primes whose largest prime factor is small): rho degenerates
        # (periodic log-vector, only false collisions). Recover exhaustively
        # instead; O(q) is trivial at this size and the goal is exercising the
        # recovery/verification, not rho at q < 64.
        x = next((cand for cand in range(q) if pow(G, cand, n) == H % n), None)
        if x is None:
            raise RuntimeError("exhaustive recovery failed")
        method = "exhaustive (tiny subgroup, q < 64)"
    else:
        x = pollard_rho(G, H, n, q)
        method = "pollard-rho (full table)"
    assert x % q == x0 % q, "recovery mismatch"
    print(json.dumps({
        "group": "prime-order subgroup of Z_n^*",
        "order": n,
        "subgroup_order": q,
        "generator": G,
        "secret": x0,
        "recovered": x,
        "method": method,
        "all_recovered": True,
        "scope": "toy only; does not validate production cost bounds",
    }))


if __name__ == "__main__":
    main()
