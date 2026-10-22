
"""Recompute the declared work bound and run a small adding walk.

The walk is evidence about three prime-order additive groups only.
"""
import json
import math
import sys
from decimal import Decimal, localcontext

TRACK = "p256-ecdlp"
ORDER = 115792089210356248762697446949407573529996955224135760342422259061068512044369
C = 1000
M_INS = 1024
A_INS = 128
GM, GA = 16, 28
EXPECTED_W = 72912983316395138025578212591182278225982
EXPECTED_NUM = 135743300
EXPECTED_DEN = 1000000
EXPECTED_ITER = 1020847100643973146639479902990327873536
EXPECTED_SQRT = 340282366881324382213159967663442624512


def ceil_sqrt(n):
    s = math.isqrt(n)
    return s if s * s == n else s + 1


def ceil_div(a, b):
    return (a + b - 1) // b


def recompute():
    iterations = 3 * ceil_sqrt(ORDER)
    walk = 9 * ceil_sqrt(ORDER)
    total_ops = walk + 64 * 520 + 520
    g_ins = GM * M_INS + GA * A_INS
    eq = iterations * (8 * M_INS + 8 * A_INS)
    coeff = iterations * 3 * 2 * 256
    hash_ins = iterations * 3 * 256
    reserve = 1 << 40
    total_ins = total_ops * g_ins + eq + coeff + hash_ins + reserve
    w = ceil_div(total_ins, C)
    return {
        "iterations": iterations,
        "sqrt": ceil_sqrt(ORDER),
        "W": w,
        "group_ins": g_ins,
    }


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


def toys(samples):
    report = []
    for order in (10007, 100003, 1000003):
        hits = 0
        for i in range(samples):
            secret = 1 + (i * 17) % (order - 1)
            seed = 0xA5A5A5A5 ^ (i * 0x9E3779B9)
            hits += bool(floyd_once(order, secret, seed))
        if hits * 10 < samples * 9:
            raise SystemExit(f"walk success {hits}/{samples} on order {order} is below 0.90")
        report.append({"order": order, "hits": hits, "samples": samples})
    return report


def main():
    payload = json.load(sys.stdin)
    samples = int(payload["samples"])
    if not 1 <= samples <= 30:
        raise SystemExit("samples out of declared range")
    got = recompute()
    if got["W"] != EXPECTED_W or got["iterations"] != EXPECTED_ITER or got["sqrt"] != EXPECTED_SQRT:
        raise SystemExit(f"recomputed bound {got} does not match the claim")
    log_ok(got["W"], EXPECTED_NUM, EXPECTED_DEN)
    report = toys(samples)
    print(json.dumps({
        "track": TRACK,
        "W": str(got["W"]),
        "iterations": str(got["iterations"]),
        "score": [EXPECTED_NUM, EXPECTED_DEN],
        "toys": report,
        "scope": "toy groups and arithmetic identity only",
    }))


if __name__ == "__main__":
    main()
