
"""Recompute the declared GT work bound and run the same small adding walk."""
import json
import math
import sys
from decimal import Decimal, localcontext

TRACK = "bls12-381-fixed-g1-invert-g2"
ORDER = 52435875175126190479447740508185965837690552500527637822603658699938581184513
C = 20000
EXPECTED_W = 2395039763148592466895916532144209349628
EXPECTED_H = 2060899291373846202583682192643547253372
EXPECTED_R = 6682809435494925286244686790013241925107712
EXPECTED_NUM = 130815246
EXPECTED_DEN = 1000000
SCALAR_INS = 21299200
PAIRING_ONE = 9343


def ceil_sqrt(n):
    s = math.isqrt(n)
    return s if s * s == n else s + 1


def ceil_div(a, b):
    return (a + b - 1) // b


def recompute():
    iterations = 3 * ceil_sqrt(ORDER)
    walk = 9 * ceil_sqrt(ORDER)
    setup = 64 * 510
    verify = 510
    pairings = 2 * PAIRING_ONE
    h = walk + setup + verify + pairings
    eq = iterations * 2048
    coeff = iterations * 3 * 2 * 256
    hash_ins = iterations * 3 * 2048
    reserve = 1 << 40
    r_ins = eq + coeff + hash_ins + SCALAR_INS + reserve
    w = h + ceil_div(r_ins, C)
    return {"W": w, "H": h, "R": r_ins, "iterations": iterations}


def log_ok(value, num, den):
    with localcontext() as ctx:
        ctx.prec = 80
        log2 = Decimal(value).ln() / Decimal(2).ln()
        claimed = Decimal(num) / Decimal(den)
    if claimed < log2:
        raise SystemExit("declared score is below log2(W)")
    if claimed - log2 > Decimal("0.000002"):
        raise SystemExit("declared score is more than two ticks above log2(W)")


def floyd_once(order, secret, seed):
    def lcg(x):
        return (x * 6364136223846793005 + 1) & ((1 << 64) - 1)

    state = seed & ((1 << 64) - 1)
    jumps = []
    for _ in range(32):
        state = lcg(state)
        c = (state % (order - 1)) + 1
        state = lcg(state)
        d = state % order
        jumps.append((c, d))

    def step(point, a, b):
        h = ((point * 0x9E3779B97F4A7C15 + seed) >> 32) & 31
        c, d = jumps[h]
        return (point + c + d * secret) % order, (a + c) % order, (b + d) % order

    tp = ta = tb = 0
    hp = ha = hb = 0
    for _ in range(3 * ceil_sqrt(order)):
        tp, ta, tb = step(tp, ta, tb)
        hp, ha, hb = step(*step(hp, ha, hb))
        if tp == hp and (ta, tb) != (ha, hb):
            db = (tb - hb) % order
            if db == 0:
                return False
            da = (ha - ta) % order
            return da * pow(db, -1, order) % order == secret
    return False


def main():
    payload = json.load(sys.stdin)
    samples = int(payload["samples"])
    if not 1 <= samples <= 30:
        raise SystemExit("samples out of declared range")
    got = recompute()
    if got["W"] != EXPECTED_W or got["H"] != EXPECTED_H or got["R"] != EXPECTED_R:
        raise SystemExit(f"recomputed bound {got} does not match the claim")
    log_ok(got["W"], EXPECTED_NUM, EXPECTED_DEN)
    toys = []
    for order in (10007, 100003, 1000003):
        hits = 0
        for i in range(samples):
            secret = 1 + (i * 17) % (order - 1)
            seed = 0xA5A5A5A5 ^ (i * 0x9E3779B9)
            hits += bool(floyd_once(order, secret, seed))
        if hits * 10 < samples * 9:
            raise SystemExit(f"walk success {hits}/{samples} on order {order} is below 0.90")
        toys.append({"order": order, "hits": hits, "samples": samples})
    print(json.dumps({
        "track": TRACK,
        "W": str(got["W"]),
        "H": str(got["H"]),
        "score": [EXPECTED_NUM, EXPECTED_DEN],
        "toys": toys,
        "scope": "toy groups and arithmetic identity only",
    }))


if __name__ == "__main__":
    main()
