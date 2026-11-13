"""Toy 16-walk parallel distinguished-point rho experiment (multiplicative
prime-order subgroup of F_p* — the faithful abstract-group model for
ECDLP-style walks). Evidence for H1 only.

Runs the claimed construction exactly — W=16 interleaved walks with
power-of-two multiplier bookkeeping, low5 selector on the group element's
integer value, shared DP table, inter-walk match recovery with assertion —
on a prime-order subgroup of F_p*, and reports the total visit count
(summed over walks) at the first DP-table match against the corrected
union birthday model sqrt(pi*N/2)*sqrt(W/(W-1)). Small-scope evidence only.
"""
import json
import math
import random
import sys

W = 16


def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def find_group(bits: int, rng: random.Random):
    """Safe-ish prime p = 2q+1 with q prime near `bits` bits; returns
    (p, q, g) with g a generator of the order-q subgroup of F_p*."""
    while True:
        q = rng.getrandbits(bits - 1) | (1 << (bits - 2)) | 1
        if not is_probable_prime(q):
            continue
        p = 2 * q + 1
        if is_probable_prime(p):
            g = 2
            while True:
                if is_probable_prime(g) and pow(g, 2, p) != 1 and pow(g, q, p) != 1:
                    return p, q, pow(g, 2, p)
                g += 1


def run_trial(p: int, q: int, g: int, dp_bits: int, rng: random.Random, cap: int):
    """One 16-walk run in the order-q subgroup <g> of F_p*.
    Secret x in [0,q); public Q = g^x. Walk i holds (a_i, b_i) with point
    P_i = g^(a_i*x + b_i) mod p, updated multiplicatively by
    R_j = g^(2^aj*x + 2^bj). DP on the integer value of P_i.
    Returns (total_visits, recovered)."""
    x = rng.randrange(q)
    Q = pow(g, x, p)
    m = 32
    aj = rng.sample(range(0, 2 * q.bit_length()), m)
    bj = rng.sample(range(0, 2 * q.bit_length()), m)
    addpoints = [Q << 0 for _ in range(m)]
    addpoints = [pow(Q, pow(2, aj[j], q), p) * pow(g, pow(2, bj[j], q), p) % p for j in range(m)]
    dp_mask = (1 << dp_bits) - 1
    table = {}
    state = [(rng.randrange(q), rng.randrange(q)) for _ in range(W)]
    points = [pow(g, a * x % q + b, p) for (a, b) in state]
    visits = 0
    per_walk = [0] * W
    w = 0
    while True:
        i = w % W
        if per_walk[i] >= cap:
            w += 1
            if w >= W * 4:
                return visits, False
            continue
        a, b = state[i]
        P = points[i]
        visits += 1
        per_walk[i] += 1
        if P & dp_mask == 0:
            rec = table.get(P)
            if rec is not None:
                (sa, sb, si) = rec
                denom = (a - sa) % q
                assert denom != 0, "degenerate collision"
                xhat = ((sb - b) * pow(denom, -1, q)) % q
                assert xhat == x, "recovery failed"
                return visits, True
            table[P] = (a % q, b % q, i)
        j = P & 31
        P = P * addpoints[j] % p
        a += 1 << aj[j]
        b += 1 << bj[j]
        a %= 4 * q
        b %= 4 * q
        state[i] = (a, b)
        points[i] = P
        w += 1


def main() -> None:
    req = json.load(sys.stdin)
    bits = int(req["bits"])
    trials = int(req["trials"])
    dp_bits = 6
    rng = random.Random(0xC0FFEE ^ bits)
    p, q, g = find_group(bits, rng)
    total = 0
    maxt = 0
    completed = 0
    capped = 0
    model = math.sqrt(math.pi * q / 2) * math.sqrt(W / (W - 1))
    cap = int(30 * model) + (1 << 16)
    for _ in range(trials):
        t, ok = run_trial(p, q, g, dp_bits, rng, cap)
        if ok:
            total += t
            maxt = max(maxt, t)
            completed += 1
        else:
            capped += 1
    mean = total / completed if completed else 0.0
    out = {
        "group": "prime-order subgroup of F_p* (multiplicative)",
        "walks": W,
        "order_bits": q.bit_length(),
        "order": q,
        "dp_bits": dp_bits,
        "trials": trials,
        "completed": completed,
        "capped_failures": capped,
        "mean_total_visits": mean,
        "union_birthday_model": model,
        "ratio_mean_over_model": mean / model,
        "max_total_visits": maxt,
        "visit_cap_per_walk": cap,
        "scope": "toy group only; validates the value-sequence model H1 of the claimed 16-walk construction, not production costs",
    }
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
