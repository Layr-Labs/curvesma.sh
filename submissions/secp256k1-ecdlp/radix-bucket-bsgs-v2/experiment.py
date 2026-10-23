"""Toy check of truncated radix-sorted BSGS on a small additive cyclic group.

Not evidence of production curve arithmetic, memory, or instruction totals.
The production bound is the integer ledger in proof.md.
"""
import json
import math
import sys

BITS = 8


def counting_sort(records, shift):
    mask = (1 << BITS) - 1
    count = [0] * (mask + 1)
    for key, _index, _sign in records:
        count[(key >> shift) & mask] += 1
    total = 0
    for digit, card in enumerate(count):
        count[digit] = total
        total += card
    out = [None] * len(records)
    for record in records:
        digit = (record[0] >> shift) & mask
        out[count[digit]] = record
        count[digit] += 1
    return out


def radix_sort(records):
    if not records:
        return []
    width = max(key for key, _index, _sign in records).bit_length()
    shift = 0
    out = list(records)
    while True:
        out = counting_sort(out, shift)
        shift += BITS
        if shift >= max(width, 1):
            break
    return out


def shape(order):
    covered = (39 * order + 99) // 100
    m = math.isqrt(covered)
    if m * m < covered:
        m += 1
    if m < 1:
        m = 1
    if m > order:
        m = order
    intervals = (covered + m - 1) // m
    if m * intervals > order:
        intervals = (order + m - 1) // m
    return covered, m, intervals


def lookup(babies, value):
    lo = 0
    hi = len(babies)
    while lo < hi:
        mid = (lo + hi) // 2
        if babies[mid][0] < value:
            lo = mid + 1
        else:
            hi = mid
    while lo < len(babies) and babies[lo][0] == value:
        yield babies[lo][1]
        lo += 1


def recover(secret, order, m, intervals, babies):
    value = secret
    for i in range(intervals):
        for index in lookup(babies, value):
            scalar = (i * m + index) % order
            if scalar == secret:
                return scalar
        value = (value - m) % order
    return None


def main():
    order = json.load(sys.stdin)["order"]
    covered, m, intervals = shape(order)
    # The toy group has one point per value, so the sign byte is always 0.
    babies = radix_sort([(j % order, j, 0) for j in range(m)])
    if babies != sorted((j % order, j, 0) for j in range(m)):
        raise AssertionError("radix order")
    limit = min(order, m * intervals)
    checked = 0
    for secret in range(order):
        got = recover(secret, order, m, intervals, babies)
        if secret < limit:
            if got != secret:
                raise AssertionError("missed covered scalar")
            checked += 1
        elif got is not None and got != secret:
            raise AssertionError("wrong scalar")
    if checked * 100 < 39 * order:
        raise AssertionError("success below 39/100")
    print(json.dumps({
        "group": "integers modulo order under addition",
        "order": order,
        "m": m,
        "intervals": intervals,
        "covered": checked,
        "success_on_this_order": checked / order,
        "all_covered_recovered": True,
        "scope": "toy group only; does not validate production cost bounds",
    }))


if __name__ == "__main__":
    main()
