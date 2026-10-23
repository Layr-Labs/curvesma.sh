"""Odd-lift BSGS with cross-page planted-base MSD byte distribution.

No heap, no port queue, no 16-ary merge. Each MSD pass is one histogram,
one exclusive-prefix planted base, and one stable scatter whose destination
index is global, so a record may change page. Unary buckets are not
rewritten. The toy curve is exhaustive; the synthetic table stresses pages.
Production instruction totals are not measured here.
"""

from __future__ import annotations

import json
import sys

PAGE = 32
RADIX = 256

# y^2 = x^3 + 2x + 2 over F_17. Affine points: 18. Group order: 19.
P = 17
A = 2
B = 2
ORDER = 19
G = (0, 6)
T_BITS = 4


def inv(x: int) -> int:
    x %= P
    if x == 0:
        raise ZeroDivisionError
    return pow(x, P - 2, P)


def neg(point):
    if point is None:
        return None
    return (point[0], (-point[1]) % P)


def add(p1, p2):
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    x1, y1 = p1
    x2, y2 = p2
    if x1 == x2 and (y1 + y2) % P == 0:
        return None
    if p1 == p2:
        lam = ((3 * x1 * x1 + A) * inv(2 * y1)) % P
    else:
        lam = ((y2 - y1) * inv(x2 - x1)) % P
    x3 = (lam * lam - x1 - x2) % P
    y3 = (lam * (x1 - x3) - y1) % P
    return (x3, y3)


def mul(k: int, point):
    if k < 0:
        return mul(-k, neg(point))
    acc = None
    base = point
    while k:
        if k & 1:
            acc = add(acc, base)
        base = add(base, base)
        k >>= 1
    return acc


def point_count() -> int:
    total = 1
    for x in range(P):
        rhs = (x * x * x + A * x + B) % P
        if pow(rhs, (P - 1) // 2, P) == 0:
            total += 1
        elif pow(rhs, (P - 1) // 2, P) == 1:
            total += 2
    return total


def n_babies(t_bits: int) -> int:
    return 1 << (t_bits - 1)


def n_giants(order: int, t_bits: int) -> int:
    modulus = 1 << t_bits
    need = (2 * order + modulus - 1) // modulus
    return (need + 2) // 2


def odd_representative(scalar: int, order: int) -> int:
    if scalar == 0:
        return 0
    if scalar % 2 == 1:
        return scalar
    return scalar + order


def representations(z: int, t_bits: int, giants: int):
    modulus = 1 << t_bits
    babies = n_babies(t_bits)
    hits = []
    for i in range(giants):
        base = (2 * i + 1) * modulus
        for j in range(babies):
            odd = 2 * j + 1
            if base + odd == z:
                hits.append((i, j, 1))
            if base - odd == z:
                hits.append((i, j, -1))
    return hits


def planted_msd(records, key_of, key_len: int, counter: dict):
    """Recursive MSD. `records` is a list of payloads. Returns a new list.

    Planted bases are the exclusive prefix of the 256-bin histogram.
    Scatter destinations are absolute indexes in the bucket buffer, which
    is the global order once the top call returns, so writes are cross-page
    whenever source and destination pages differ. No heap is built.
    """

    def rec(items, depth: int, origin: int):
        if len(items) <= 1 or depth == key_len:
            return items
        counts = [0] * RADIX
        for index, item in enumerate(items):
            counts[key_of(item)[depth]] += 1
            counter["histogram_increments"] += 1
            src = origin + index
            counter["reads"] += 1
            if src // PAGE != index // PAGE:
                # local index is not yet a global destination; real
                # cross-page events are counted on the scatter below.
                pass
        occupied = [b for b in range(RADIX) if counts[b]]
        if len(occupied) == 1:
            counter["unary_skips"] += 1
            return rec(items, depth + 1, origin)
        base = [0] * RADIX
        running = 0
        for b in range(RADIX):
            base[b] = running
            running += counts[b]
            counter["prefix_steps"] += 1
        counter["planted_directories"] += 1
        buf = [None] * len(items)
        cursor = base[:]
        for index, item in enumerate(items):
            digit = key_of(item)[depth]
            dest = cursor[digit]
            cursor[digit] += 1
            buf[dest] = item
            counter["scatter_writes"] += 1
            src_page = (origin + index) // PAGE
            dst_page = (origin + dest) // PAGE
            if src_page != dst_page:
                counter["cross_page_writes"] += 1
            span = counts[digit]
            leaf_pages = (span + PAGE - 1) // PAGE
            if leaf_pages > 1:
                counter["multi_page_buckets"] += 1
        out = []
        for b in occupied:
            start = base[b]
            stop = start + counts[b]
            out.extend(rec(buf[start:stop], depth + 1, origin + start))
        return out

    counter.setdefault("histogram_increments", 0)
    counter.setdefault("scatter_writes", 0)
    counter.setdefault("cross_page_writes", 0)
    counter.setdefault("planted_directories", 0)
    counter.setdefault("prefix_steps", 0)
    counter.setdefault("unary_skips", 0)
    counter.setdefault("multi_page_buckets", 0)
    counter.setdefault("reads", 0)
    return rec(records, 0, 0)


def binary_search(table, x: int):
    lo = 0
    hi = len(table)
    probes = 0
    while lo < hi:
        mid = (lo + hi) // 2
        probes += 1
        if table[mid][0] < x:
            lo = mid + 1
        else:
            hi = mid
    if lo < len(table) and table[lo][0] == x:
        return table[lo], probes
    return None, probes


def build_babies(t_bits: int, counter: dict):
    babies = n_babies(t_bits)
    step = add(G, G)
    counter["doublings"] += 1
    point = G
    rows = []
    for j in range(babies):
        odd = 2 * j + 1
        rows.append((point[0], point[1], odd))
        if j + 1 != babies:
            point = add(point, step)
            counter["baby_additions"] += 1
    xs = [row[0] for row in rows]
    if len(xs) != len(set(xs)):
        raise AssertionError("odd babies share an x-coordinate")
    stat = {}
    ordered = planted_msd(rows, lambda row: bytes([row[0]]), 1, stat)
    for key, value in stat.items():
        counter[key] = counter.get(key, 0) + value
    counter["sorts"] += 1
    return ordered


def recover(public, t_bits: int, counter: dict):
    if public is None:
        return 0, 0
    table = build_babies(t_bits, counter)
    modulus = 1 << t_bits
    giants = n_giants(ORDER, t_bits)
    # Sequential doubling yields 2^t G and the stride 2^{t+1} G.
    pow2 = G
    for _ in range(t_bits):
        pow2 = add(pow2, pow2)
        counter["doublings"] += 1
    offset = pow2  # 2^t G
    stride = add(pow2, pow2)
    counter["doublings"] += 1
    point = add(public, neg(offset))
    counter["giant_additions"] += 1
    probes = 0
    for i in range(giants):
        if point is not None:
            hit, used = binary_search(table, point[0])
            probes += used
            counter["lookups"] += 1
            counter["lookup_probes"] += used
            if hit is not None:
                _x, y, odd = hit
                if point[1] == y:
                    sign = 1
                elif point[1] == (-y) % P:
                    sign = -1
                else:
                    sign = 0
                if sign != 0:
                    base = (2 * i + 1) * modulus
                    return (base + sign * odd) % ORDER, probes
        if i + 1 != giants:
            point = add(point, neg(stride))
            counter["giant_additions"] += 1
    raise AssertionError("odd-lift search missed a covered scalar")


def check_curve():
    if point_count() != ORDER:
        raise AssertionError(point_count())
    if mul(ORDER, G) is not None:
        raise AssertionError("order")
    for k in range(1, ORDER):
        if mul(k, G) is None:
            raise AssertionError(k)
    if add(G, neg(G)) is not None:
        raise AssertionError("negation")
    if mul(2, G) != add(G, G):
        raise AssertionError("double")


def check_cover():
    giants = n_giants(ORDER, T_BITS)
    if n_babies(T_BITS) != 8 or giants != 2:
        raise AssertionError((n_babies(T_BITS), giants))
    for z in range(1, 2 * ORDER, 2):
        hits = representations(z, T_BITS, giants)
        if len(hits) != 1:
            raise AssertionError((z, hits))
    for z in range(0, 2 * ORDER, 2):
        if representations(z, T_BITS, giants):
            raise AssertionError(z)
    for scalar in range(ORDER):
        z = odd_representative(scalar, ORDER)
        if scalar == 0:
            if z != 0:
                raise AssertionError(z)
            continue
        if z % 2 != 1 or not (1 <= z < 2 * ORDER):
            raise AssertionError((scalar, z))
        if len(representations(z, T_BITS, giants)) != 1:
            raise AssertionError((scalar, z))


def check_recover():
    counter = {
        "doublings": 0,
        "baby_additions": 0,
        "giant_additions": 0,
        "lookups": 0,
        "lookup_probes": 0,
        "sorts": 0,
    }
    found, _ = recover(None, T_BITS, counter)
    if found != 0:
        raise AssertionError(found)
    for scalar in range(ORDER):
        got, _ = recover(mul(scalar, G), T_BITS, counter)
        if got != scalar:
            raise AssertionError((scalar, got))
    return counter


def check_distribution():
    records = []
    for i in range(4096):
        key = bytes([i % 64]) + i.to_bytes(3, "big")
        records.append((key, i))
    counter = {}
    ordered = planted_msd(records, lambda row: row[0], 4, counter)
    keys = [row[0] for row in ordered]
    if keys != sorted(keys):
        raise AssertionError("msd order")
    if sorted(row[1] for row in ordered) != list(range(4096)):
        raise AssertionError("permutation")
    if counter["cross_page_writes"] <= 0:
        raise AssertionError("expected cross-page scatters")
    if counter["multi_page_buckets"] <= 0:
        raise AssertionError("expected a bucket spanning pages")
    if counter["planted_directories"] <= 0:
        raise AssertionError("missing planted base")
    if counter["scatter_writes"] > 4096 * 4:
        raise AssertionError(counter["scatter_writes"])
    if counter["histogram_increments"] > 4096 * 4:
        raise AssertionError(counter["histogram_increments"])
    # Identical keys: one planted directory chain, unary until the last byte
    # is unnecessary because the single bucket never splits.
    same = [(b"\x00\x00\x00\x00", i) for i in range(128)]
    same_counter = {}
    same_ordered = planted_msd(same, lambda row: row[0], 4, same_counter)
    if [row[1] for row in same_ordered] != list(range(128)):
        raise AssertionError("stable unary")
    if same_counter["scatter_writes"] != 0:
        raise AssertionError(same_counter["scatter_writes"])
    return counter


def production_shape():
    """Closed forms for the P-256 instance. Not executed as a search."""
    babies = 1 << 127
    giants = 1 << 128
    key_len = 32
    record_bytes = 80
    return {
        "babies": babies,
        "giants": giants,
        "modulus": 1 << 128,
        "key_len": key_len,
        "page": PAGE,
        "radix": RADIX,
        "record_bytes": record_bytes,
        "baby_additions": babies - 1,
        "doublings": 129,
        "giant_additions": giants,
        "lookups": giants,
        "scatter_writes_bound": babies * key_len,
        "histogram_increments_bound": babies * key_len,
        "memory_bytes_bound": babies * record_bytes * 2,
    }


def run_experiment():
    check_curve()
    check_cover()
    dlp = check_recover()
    dist = check_distribution()
    shape = production_shape()
    if shape["memory_bytes_bound"].bit_length() - 1 >= 143:
        raise AssertionError("memory")
    return {
        "ok": True,
        "success": 1,
        "advice_bits": 0,
        "dlp": {k: dlp[k] for k in ("doublings", "baby_additions", "giant_additions", "lookups", "sorts")},
        "distribution": {
            "cross_page_writes": dist["cross_page_writes"],
            "planted_directories": dist["planted_directories"],
            "scatter_writes": dist["scatter_writes"],
            "multi_page_buckets": dist["multi_page_buckets"],
        },
        "production": {
            "babies": "2^127",
            "giants": "2^128",
            "scatter_writes_bound": "2^127 * 32",
            "memory_bytes_bound": "80 * 2^128",
        },
    }


def run_mode(mode: str):
    if mode == "cover":
        check_curve()
        check_cover()
        return {"ok": True, "mode": mode, "scope": "toy curve only"}
    if mode == "distribution":
        dist = check_distribution()
        return {
            "ok": True,
            "mode": mode,
            "cross_page_writes": dist["cross_page_writes"],
            "scope": "synthetic keys only",
        }
    if mode == "all":
        result = run_experiment()
        result["mode"] = mode
        result["scope"] = "toy curve and synthetic keys; not a production cost measurement"
        return result
    raise AssertionError(mode)


def main():
    spec = json.load(sys.stdin)
    result = run_mode(spec["mode"])
    if not result["ok"]:
        raise SystemExit(1)
    json.dump(result, sys.stdout)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
