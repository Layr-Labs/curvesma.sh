"""Toy-scale negation-class Pollard rho verification; checks scalar recovery on finite cyclic groups."""
import hashlib
import json
import math
import sys


def solve_toy_rho(order, secret):
    salt = b"secp256k1-toy-salt"
    increments = []
    for j in range(8):
        c_j = (j * 17 + 3) % order
        d_j = (j * 23 + 7) % order
        val_j = (c_j + secret * d_j) % order
        increments.append((c_j, d_j, val_j))

    def step(y, a, b):
        canon = min(y, (order - y) % order)
        h = hashlib.sha256(salt + canon.to_bytes(4, "big")).digest()[0] % 8
        c_j, d_j, val_j = increments[h]
        return (y + val_j) % order, (a + c_j) % order, (b + d_j) % order

    t_y, t_a, t_b = secret, 0, 1
    h_y, h_a, h_b = secret, 0, 1

    max_steps = 12 * math.isqrt(order) + 64
    for _ in range(max_steps):
        t_y, t_a, t_b = step(t_y, t_a, t_b)
        h_y, h_a, h_b = step(h_y, h_a, h_b)
        h_y, h_a, h_b = step(h_y, h_a, h_b)

        t_c = min(t_y, (order - t_y) % order)
        h_c = min(h_y, (order - h_y) % order)
        if t_c == h_c:
            if t_y == h_y:
                denom = (t_b - h_b) % order
                num = (h_a - t_a) % order
            else:
                denom = (t_b + h_b) % order
                num = (-h_a - t_a) % order
            if math.gcd(denom, order) == 1:
                recovered = (num * pow(denom, -1, order)) % order
                if recovered == secret:
                    return True
    return False


def main():
    payload = json.load(sys.stdin)
    order = payload["order"]
    tested = min(order - 1, 32)
    recovered = sum(solve_toy_rho(order, s) for s in range(1, tested + 1))
    print(
        json.dumps(
            {
                "group": "cyclic additive group",
                "order": order,
                "tested_secrets": tested,
                "recovered_secrets": recovered,
                "all_tested_recovered": (recovered == tested),
                "scope": "toy-scale functional check of equivalence-class walk and recovery",
            }
        )
    )


if __name__ == "__main__":
    main()
