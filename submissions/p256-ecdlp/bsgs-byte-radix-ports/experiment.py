"""Odd-lift coverage on a toy additive group, plus an in-page byte radix.

Not a production cost measurement. The radix uses 256 ports of capacity 4
and must match a stable sort on the same integer keys.
"""
import json
import sys

PORT_CAP = 4
NPORTS = 256


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


def _byte(key, shift):
    return (key >> shift) & 0xFF


def radix_pass(keys, shift):
    """Stable counting pass. Ports spill at PORT_CAP in arrival order."""
    hist = [0] * NPORTS
    for key in keys:
        hist[_byte(key, shift)] += 1
    start = [0] * NPORTS
    running = 0
    for digit in range(NPORTS):
        start[digit] = running
        running += hist[digit]
    out = [0] * len(keys)
    cursor = start[:]
    ports = [[] for _ in range(NPORTS)]

    def flush(digit):
        for key in ports[digit]:
            out[cursor[digit]] = key
            cursor[digit] += 1
        ports[digit].clear()

    for key in keys:
        digit = _byte(key, shift)
        ports[digit].append(key)
        if len(ports[digit]) == PORT_CAP:
            flush(digit)
    for digit in range(NPORTS):
        if ports[digit]:
            flush(digit)
    return out


def radix_sort(keys, nbytes):
    out = list(keys)
    for pass_index in range(nbytes):
        out = radix_pass(out, 8 * pass_index)
    return out


def check_radix():
    keys = [
        0x10F3, 0x00F3, 0x10F4, 0x0001, 0xABCD, 0x00F3, 0x7FFF, 0x0100,
        0x0000, 0xFF00, 0x00FF, 0x8000, 0x7FFF, 0x0002,
    ]
    got = radix_sort(keys, 2)
    stable = sorted(keys)
    if got != stable:
        raise SystemExit(f"radix mismatch {got} != {stable}")
    # A full page of 32 records, 4-byte keys, the production page width.
    page = [(i * 0x9E3779B1 + 0x00F3) & 0xFFFFFFFF for i in range(32)]
    page[3] = page[7] = page[0]
    if radix_sort(page, 4) != sorted(page):
        raise SystemExit("page radix mismatch")


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
    check_radix()
    print(json.dumps({
        "group": "integers modulo odd order under addition",
        "order": order,
        "checked": checked,
        "identity_recovered": True,
        "matches_per_scalar": max_matches,
        "all_recovered": True,
        "radix_ports": NPORTS,
        "port_cap": PORT_CAP,
        "page_records_checked": 32,
        "scope": "toy odd-lift identity and in-page radix only; does not validate production instruction costs",
    }))


if __name__ == "__main__":
    main()
