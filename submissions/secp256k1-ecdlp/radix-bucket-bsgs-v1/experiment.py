"""Truncated BSGS coverage check on a toy cyclic group.

Not evidence of production curve arithmetic, memory, or instruction totals.
The production bound is the integer ledger in proof.md.
"""
import json
import math
import sys


def counting_sort(items, shift, bits=8):
    mask = (1 << bits) - 1
    count = [0] * (mask + 1)
    for key, _payload in items:
        count[(key >> shift) & mask] += 1
    total = 0
    for digit, card in enumerate(count):
        count[digit] = total
        total += card
    out = [None] * len(items)
    for item in items:
        digit = (item[0] >> shift) & mask
        out[count[digit]] = item
        count[digit] += 1
    return out


def radix_sort(items, bits=8):
    if not items:
        return []
    width = max(key for key, _payload in items).bit_length()
    shift = 0
    out = list(items)
    while True:
        out = counting_sort(out, shift, bits)
        shift += bits
        if shift >= max(width, 1):
            break
    return out


def shape(order):
    covered_goal = (39 * order + 99) // 100
    m = math.isqrt(covered_goal)
    if m * m < covered_goal:
        m += 1
    if m < 1:
        m = 1
    if m > order:
        m = order
    intervals = (covered_goal + m - 1) // m
    if m * intervals > order:
        intervals = (order + m - 1) // m
    return m, intervals


def recover(secret, order, m, intervals, babies):
    value = secret
    probes = max(len(babies).bit_length(), 1)
    for i in range(intervals):
        lo = 0
        hi = len(babies)
        found = None
        for _ in range(probes):
            if lo >= hi:
                break
            mid = (lo + hi) // 2
            key = babies[mid][0]
            if key == value:
                found = babies[mid][1]
                break
            if key < value:
                lo = mid + 1
            else:
                hi = mid
        if found is not None:
            scalar = i * m + found
            if scalar % order == secret:
                return scalar % order
        value = (value - m) % order
    return None


def main():
    order = json.load(sys.stdin)["order"]
    m, intervals = shape(order)
    babies = radix_sort([(j % order, j) for j in range(m)])
    if babies != sorted((j % order, j) for j in range(m)):
        raise AssertionError("radix order")
    limit = min(order, m * intervals)
    covered = 0
    for secret in range(order):
        got = recover(secret, order, m, intervals, babies)
        if secret < limit:
            if got != secret:
                raise AssertionError("missed covered scalar")
            covered += 1
        elif got is not None and got != secret:
            raise AssertionError("wrong scalar")
    if covered * 100 < 39 * order:
        raise AssertionError("success below 39/100")
    print(json.dumps({
        "group": "integers modulo order under addition",
        "order": order,
        "m": m,
        "intervals": intervals,
        "covered": covered,
        "success_on_this_order": covered / order,
        "all_covered_recovered": True,
        "scope": "toy group only; does not validate production cost bounds",
    }))


if __name__ == "__main__":
    main()
