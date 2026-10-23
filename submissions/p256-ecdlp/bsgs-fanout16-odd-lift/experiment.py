"""Odd-lift coverage on a toy additive group. Not a production cost measurement."""
import json
import sys


def recover(order, secret):
    # Largest power of two whose signed odds stay inside the toy order.
    bound = 2
    while 2 * ((bound * 2) - 1) < order:
        bound *= 2
    length = bound
    half = length // 2
    giants = (order + length - 1) // length
    lift = secret if secret % 2 == 1 else secret + order
    matches = 0
    recovered = None
    for index in range(giants):
        center = (2 * index + 1) * length
        delta = lift - center
        if not (-length < delta < length and delta % 2 != 0):
            continue
        baby = (abs(delta) - 1) // 2
        if not 0 <= baby < half:
            continue
        raw = center + delta
        reduced = raw
        guards = 0
        while reduced >= order:
            reduced -= order
            guards += 1
            if guards > 4:
                break
        matches += 1
        if reduced == secret:
            recovered = reduced
    return recovered, matches, length, half, giants


def main():
    order = json.load(sys.stdin)["order"]
    if order % 2 == 0 or order < 5:
        raise SystemExit("order must be an odd integer >= 5")
    checked = 0
    max_matches = 0
    for secret in range(order):
        recovered, matches, length, half, giants = recover(order, secret)
        if recovered != secret or matches != 1:
            raise SystemExit(f"miss {secret} -> {recovered} matches={matches}")
        if matches > max_matches:
            max_matches = matches
        checked += 1
    # Every leaf index of a complete 16-ary tree fits in ceil(log16(n)) digits.
    leaves = 4096
    width = 16
    digits = 0
    span = 1
    while span < leaves:
        span *= width
        digits += 1
    for leaf in range(leaves):
        value = leaf
        count = 0
        while value or count == 0:
            value //= width
            count += 1
            if value == 0:
                break
        if count > digits:
            raise SystemExit(f"path longer than {digits} for leaf {leaf}")
    print(json.dumps({
        "group": "integers modulo odd order under addition",
        "order": order,
        "checked": checked,
        "identity_recovered": True,
        "matches_per_scalar": max_matches,
        "all_recovered": True,
        "fanout_digits_checked": digits,
        "scope": "toy odd-lift identity only; does not validate production instruction costs",
    }))


if __name__ == "__main__":
    main()
