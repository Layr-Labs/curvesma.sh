"""Toy additive BSGS, plus the integer ledger in proof.md.

The production table is not built. The ledger check only re-evaluates the
published integer bounds.
"""
import json
import math
import sys
from decimal import Decimal, getcontext


def toy_bsgs(order):
    m = math.isqrt(order) + 1
    babies = {j % order: j for j in range(m)}
    for secret in range(order):
        value = secret
        found = False
        for i in range((order + m - 1) // m):
            if value in babies:
                recovered = (i * m + babies[value]) % order
                if recovered != secret:
                    raise AssertionError("recovered scalar mismatch")
                found = True
                break
            value = (value - m) % order
        if not found:
            raise AssertionError("discrete logarithm not recovered")


def ledger():
    r = int(
        "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
        16,
    )
    m = 143186686509166558328661515423460682036
    c = 20000
    side = c + 2304 + 576 + 128 + 3072
    if side != 26080:
        raise AssertionError("side")
    n = (391 * r + 999) // 1000
    if not ((m - 1) * (m - 1) < n <= m * m <= r):
        raise AssertionError("coverage window")
    if not m < 2**127:
        raise AssertionError("counter width")
    b = 2**84
    t = b * 2**193
    prod = m * m
    if 100 * prod * (t - r) < 39 * r * t:
        raise AssertionError("success inequality")
    fixed_one = 16 * (64 + 32) * 65536
    variable = (side + 64) * (2 * m)
    setup = (16384 + 3072 + 512) * c + 2**16
    verify = (16384 + 3072) * c
    cap = b * 2048
    total = variable + 2 * fixed_one + cap + setup + verify
    pre = setup + side * m + fixed_one
    mem = 1296 * m + 2**21
    published_total = 7486945464191340615970310450630709904145664
    published_pre = 3734308784159063841211492322243855087587712
    published_mem = 185569945715879859593945323988805046015808
    if total != published_total or pre != published_pre or mem != published_mem:
        raise AssertionError("ledger mismatch")
    getcontext().prec = 80
    ln2 = Decimal(2).ln()

    def is_upper(k, numer, denom):
        # numer/denom <= 2^(k/10**6) iff ln(numer)*1e6 <= ln(denom)*1e6 + ln(2)*k
        left = Decimal(numer).ln() * 1000000
        right = Decimal(denom).ln() * 1000000 + ln2 * k
        return left <= right

    if not is_upper(128137645, total, c):
        raise AssertionError("time ceiling")
    if is_upper(128137644, total, c):
        raise AssertionError("time ceiling not tight")
    if not is_upper(127134109, pre, c):
        raise AssertionError("preprocessing ceiling")
    if is_upper(127134108, pre, c):
        raise AssertionError("preprocessing ceiling not tight")
    if not is_upper(137091015, mem, 1):
        raise AssertionError("memory ceiling")
    if is_upper(137091014, mem, 1):
        raise AssertionError("memory ceiling not tight")
    # Miniature collision lemma: d*s ≡ C (mod 2^w) has 0 or 2^{v2(d)} solutions.
    w = 8
    modulus = 1 << w
    d = 12
    t = 0
    dd = d
    while dd % 2 == 0:
        dd //= 2
        t += 1
    for residue in range(modulus):
        sols = sum(1 for s in range(modulus) if (d * s - residue) % modulus == 0)
        if sols not in (0, 1 << t):
            raise AssertionError("small collision count")


def main():
    order = json.load(sys.stdin)["order"]
    toy_bsgs(order)
    ledger()
    print(
        json.dumps(
            {
                "group": "integers modulo order under addition",
                "order": order,
                "checked": order,
                "all_recovered": True,
                "ledger": True,
                "scope": "toy group and integer bounds only; not a production execution",
            }
        )
    )


if __name__ == "__main__":
    main()
