"""Toy check of the Cheon rectangle and its mu_6 lift.

The production attack searches the exponent k = log_eta(A_d) in the cyclic
group of N-th roots of unity, stores the mu_6 orbit of every baby point,
and lifts a hit k to the scalar a = k + t*N.  This script replays exactly
that logic in the multiplicative group Z_l^* for a small prime l with
d = 8, d | l-1 and 6 | N = (l-1)/d.  The group operation is multiplication
mod l; the mu_6 orbit is the six-unit set {+-1, +-lambda, +-lambda^2}.
Nothing here is production-size evidence: it checks the covering argument
and the lift rule only, at the scope of the small group.

Python standard library only.
"""

import json
import math
import sys


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True


def find_prime(bits):
    """Smallest prime l in [2^(bits-1), 2^bits) with l = 1 (mod 48)."""
    n = 1 << (bits - 1)
    while n % 48 != 1:
        n += 1
    limit = 1 << bits
    while n < limit:
        if is_prime(n):
            return n
        n += 48
    raise ValueError("no prime in range")


def prime_factors(n):
    out = set()
    q = 2
    while q * q <= n:
        while n % q == 0:
            out.add(q)
            n //= q
        q += 1
    if n > 1:
        out.add(n)
    return out


def primitive_root(l):
    order = l - 1
    factors = prime_factors(order)
    for g in range(2, l):
        if all(pow(g, order // q, l) != 1 for q in factors):
            return g
    raise ValueError("no primitive root")


def run(bits):
    d = 8
    l = find_prime(bits)
    N = (l - 1) // d
    if N % 6 != 0:
        raise ValueError("6 must divide N")
    g = primitive_root(l)
    eta = pow(g, d, l)                       # order N
    zeta = pow(g, N, l)                      # order d
    lam = pow(g, (l - 1) // 3, l)            # order 3

    units = [1, l - 1, lam, l - lam, lam * lam % l, l - lam * lam % l]
    if len(set(units)) != 6:
        raise AssertionError("mu_6 is not free")
    # The six units must be exactly eta^(t*N/6); that identity is what
    # turns an orbit hit into the shift t*N/6 in the exponent.
    shifts = sorted(pow(eta, t * N // 6, l) for t in range(6))
    if sorted(units) != shifts:
        raise AssertionError("orbit shift identity failed")

    # Rectangle side chosen so 6*G^2/N >= 0.395 and 6*G^2 < N.
    G = 1
    while 6 * G * G * 1000 < 395 * N:
        G += 1
    if 6 * G * G >= N:
        raise AssertionError("rectangles not disjoint")

    # Baby side: orbit of [eta^j] for j < G, stored as a map unit->j.
    orbit_index = {}
    for j in range(G):
        base = pow(eta, j, l)
        for u in units:
            orbit_index[u * base % l] = j

    hits = 0
    recovered = 0
    for a in range(l - 1):
        s = pow(g, a, l)
        A_d = pow(s, d, l)                   # = eta^a
        found = None
        for i in range(G):
            giant = pow(eta, (a - G * i) % N, l)
            if giant in orbit_index:
                found = (i, orbit_index[giant])
                break
        if found is None:
            continue
        hits += 1
        i, j = found
        ok = False
        for t in range(6):
            k = (G * i + j + t * N // 6) % N
            for lift in range(d):
                cand = k + lift * N
                if cand < l - 1 and pow(g, cand, l) == s:
                    ok = True
        if not ok:
            raise AssertionError("hit did not lift")
        recovered += 1

    expected = 6 * G * G / N
    observed = hits / (l - 1)
    if hits != recovered:
        raise AssertionError("recovery mismatch")
    return {
        "prime": l,
        "d": d,
        "N": N,
        "G": G,
        "scalars": l - 1,
        "hits": hits,
        "recovered": recovered,
        "expected_fraction": expected,
        "observed_fraction": observed,
        "covered_set_size": len(orbit_index) * G,
        "delta": abs(expected - observed),
    }


def main():
    params = json.load(sys.stdin)
    bits = int(params["bits"])
    result = run(bits)
    if result["recovered"] != result["hits"]:
        raise SystemExit("lift failures")
    if result["delta"] > 1e-6:
        raise SystemExit("coverage mismatch")
    print(json.dumps({
        "prime": result["prime"],
        "d": result["d"],
        "N": result["N"],
        "G": result["G"],
        "scalars": result["scalars"],
        "hits": result["hits"],
        "recovered": result["recovered"],
        "expected_fraction": result["expected_fraction"],
        "observed_fraction": result["observed_fraction"],
        "scope": "toy multiplicative group only; orbit shift and lift rule, "
                 "not production instruction counts",
    }))


if __name__ == "__main__":
    main()
