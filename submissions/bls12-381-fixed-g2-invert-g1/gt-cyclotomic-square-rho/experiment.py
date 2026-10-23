"""Toy Pollard three-way walk on a prime-order additive group.

Doubling stands in for the cyclotomic squaring branch. Adding the generator
stands in for multiplication by g. Adding the secret stands in for
multiplication by the target. This checks that a distinguished-point
collision recovers the discrete log. It does not speak to production
field costs.
"""
import json
import sys


def main():
    cfg = json.load(sys.stdin)
    n = int(cfg["order"])
    secret = int(cfg["secret"]) % n
    theta = int(cfg["theta"])
    budget = int(cfg["budget"])
    if n < 8 or not (0 < secret < n):
        raise SystemExit("bad input")
    if theta < 2 or budget < 1:
        raise SystemExit("bad budget")
    # X = a * secret + b  (mod n), starting from the target itself.
    x = secret
    a = 1
    b = 0
    seen = {}
    found = None
    steps = 0
    dps = 0
    squares = 0
    muls = 0
    for steps in range(1, budget + 1):
        branch = x % 3
        if branch == 0:
            x = (x + x) % n
            a = (a + a) % n
            b = (b + b) % n
            squares += 1
        elif branch == 1:
            x = (x + 1) % n
            b = (b + 1) % n
            muls += 1
        else:
            x = (x + secret) % n
            a = (a + 1) % n
            muls += 1
        if x % theta == 0:
            dps += 1
            if x in seen:
                pa, pb = seen[x]
                da = (a - pa) % n
                db = (pb - b) % n
                if da == 0:
                    seen[x] = (a, b)
                    continue
                try:
                    inv = pow(da, -1, n)
                except ValueError:
                    seen[x] = (a, b)
                    continue
                guess = (db * inv) % n
                if guess == secret:
                    found = guess
                    break
            seen[x] = (a, b)
    out = {
        "group": "integers modulo order under addition",
        "order": n,
        "theta": theta,
        "budget": budget,
        "steps": steps,
        "squarings": squares,
        "multiplications": muls,
        "distinguished_stored": dps,
        "recovered": found is not None,
        "match": found == secret,
        "scope": "toy group only; does not validate production cost bounds",
    }
    print(json.dumps(out))


if __name__ == "__main__":
    main()
