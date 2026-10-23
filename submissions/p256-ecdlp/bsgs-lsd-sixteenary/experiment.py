"""Odd-lift BSGS with a global LSD byte distribution and an implicit
sixteen-ary lookup. No heap, no per-bucket planted base, no port queue.
The toy curve is exhaustive. Production instruction totals are not measured.
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
        leg = pow(rhs, (P - 1) // 2, P)
        if leg == 0:
            total += 1
        elif leg == 1:
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


def lsd_distribute(records, key_of, key_len: int, counter: dict):
    """Stable global LSD. One 256-bin histogram per key byte, then one
    scatter into a second buffer. Destination indexes are global, so a
    record may change page. A byte that is constant on the whole table
    is not rewritten. No recursion and no per-bucket directory.
    """
    counter.setdefault("histogram_increments", 0)
    counter.setdefault("scatter_writes", 0)
    counter.setdefault("cross_page_writes", 0)
    counter.setdefault("prefix_steps", 0)
    counter.setdefault("unary_skips", 0)
    counter.setdefault("passes", 0)
    items = list(records)
    for depth in range(key_len - 1, -1, -1):
        counter["passes"] += 1
        counts = [0] * RADIX
        for item in items:
            counts[key_of(item)[depth]] += 1
            counter["histogram_increments"] += 1
        occupied = [b for b in range(RADIX) if counts[b]]
        if len(occupied) == 1:
            counter["unary_skips"] += 1
            continue
        base = [0] * RADIX
        running = 0
        for b in range(RADIX):
            base[b] = running
            running += counts[b]
            counter["prefix_steps"] += 1
        cursor = base[:]
        buf = [None] * len(items)
        for index, item in enumerate(items):
            digit = key_of(item)[depth]
            dest = cursor[digit]
            cursor[digit] += 1
            buf[dest] = item
            counter["scatter_writes"] += 1
            if index // PAGE != dest // PAGE:
                counter["cross_page_writes"] += 1
        items = buf
    return items


def sixteen_ary_lower(table, x: int):
    """Lower bound of x in a sorted table. Each stage splits the live
    index interval into 16 contiguous parts and binary-searches the part.
    Pivot indexes are arithmetic; no node is allocated.
    """
    lo = 0
    hi = len(table)
    probes = 0
    stages = 0
    while hi - lo >= 16:
        stages += 1
        span = hi - lo
        base, extra = divmod(span, 16)

        def start(k: int, origin: int = lo, width: int = base, rest: int = extra) -> int:
            return origin + k * width + (k if k < rest else rest)

        low, high = 0, 15
        while low < high:
            mid = (low + high + 1) // 2
            probes += 1
            if table[start(mid)][0] <= x:
                low = mid
            else:
                high = mid - 1
        lo = start(low)
        hi = start(low + 1) if low < 15 else hi
    while lo < hi:
        mid = (lo + hi) // 2
        probes += 1
        if table[mid][0] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo, probes, stages


def sixteen_ary_search(table, x: int):
    lo, probes, stages = sixteen_ary_lower(table, x)
    if lo < len(table) and table[lo][0] == x:
        return table[lo], probes, stages
    return None, probes, stages


def worst_binary_probes(span: int) -> int:
    if span <= 0:
        return 0
    mid = span // 2
    return 1 + max(worst_binary_probes(mid), worst_binary_probes(span - (mid + 1)))


def production_probe_cap() -> int:
    """Exact control cost on a table of 2^127 sorted keys."""
    span = 1 << 127
    stages = 0
    probes = 0
    while span >= 16:
        stages += 1
        probes += 4
        span //= 16
    probes += worst_binary_probes(span)
    if stages != 31 or span != 8:
        raise AssertionError((stages, span, probes))
    return probes


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
    ordered = lsd_distribute(rows, lambda row: bytes([row[0]]), 1, stat)
    for key, value in stat.items():
        counter[key] = counter.get(key, 0) + value
    keys = [row[0] for row in ordered]
    if keys != sorted(keys):
        raise AssertionError("baby lsd order")
    counter["sorts"] += 1
    return ordered


def recover(public, t_bits: int, counter: dict):
    if public is None:
        return 0, 0
    table = build_babies(t_bits, counter)
    modulus = 1 << t_bits
    giants = n_giants(ORDER, t_bits)
    pow2 = G
    for _ in range(t_bits):
        pow2 = add(pow2, pow2)
        counter["doublings"] += 1
    offset = pow2
    stride = add(pow2, pow2)
    counter["doublings"] += 1
    point = add(public, neg(offset))
    counter["giant_additions"] += 1
    probes = 0
    for i in range(giants):
        if point is not None:
            hit, used, _stages = sixteen_ary_search(table, point[0])
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
    ordered = lsd_distribute(records, lambda row: row[0], 4, counter)
    keys = [row[0] for row in ordered]
    if keys != sorted(keys):
        raise AssertionError("lsd order")
    if sorted(row[1] for row in ordered) != list(range(4096)):
        raise AssertionError("permutation")
    if counter["cross_page_writes"] <= 0:
        raise AssertionError("expected cross-page scatters")
    if counter["passes"] != 4:
        raise AssertionError(counter["passes"])
    if counter["scatter_writes"] > 4096 * 4:
        raise AssertionError(counter["scatter_writes"])
    if counter["histogram_increments"] != 4096 * 4:
        raise AssertionError(counter["histogram_increments"])
    same = [(b"\x00\x00\x00\x00", i) for i in range(128)]
    same_counter = {}
    same_ordered = lsd_distribute(same, lambda row: row[0], 4, same_counter)
    if [row[1] for row in same_ordered] != list(range(128)):
        raise AssertionError("stable unary")
    if same_counter["scatter_writes"] != 0:
        raise AssertionError(same_counter["scatter_writes"])
    probes = []
    stages = []
    for i, row in enumerate(ordered):
        found, used, nstages = sixteen_ary_search(ordered, row[0])
        if found[1] != row[1]:
            raise AssertionError((i, found, row))
        probes.append(used)
        stages.append(nstages)
        missing, _used, _n = sixteen_ary_search(ordered, b"\xff" + row[0][1:])
        if missing is not None:
            raise AssertionError("false hit")
    if max(stages) < 1:
        raise AssertionError("sixteen-ary did not run")
    if max(probes) > 64:
        raise AssertionError(max(probes))
    counter["max_lookup_probes"] = max(probes)
    counter["max_stages"] = max(stages)
    return counter


def production_shape():
    babies = 1 << 127
    giants = 1 << 128
    probes = production_probe_cap()
    if probes > 128:
        raise AssertionError(probes)
    return {
        "babies": babies,
        "giants": giants,
        "modulus": 1 << 128,
        "key_len": 32,
        "page": PAGE,
        "radix": RADIX,
        "arity": 16,
        "record_bytes": 80,
        "lookup_probe_cap": probes,
        "scatter_passes": 32,
    }


def run_experiment():
    check_curve()
    check_cover()
    dlp = check_recover()
    dist = check_distribution()
    shape = production_shape()
    return {
        "ok": True,
        "success": 1,
        "advice_bits": 0,
        "dlp": {
            k: dlp[k]
            for k in ("doublings", "baby_additions", "giant_additions", "lookups", "sorts")
        },
        "distribution": {
            "cross_page_writes": dist["cross_page_writes"],
            "scatter_writes": dist["scatter_writes"],
            "unary_skips": dist["unary_skips"],
            "max_lookup_probes": dist["max_lookup_probes"],
            "max_stages": dist["max_stages"],
        },
        "production": {
            "babies": "2^127",
            "giants": "2^128",
            "lookup_probe_cap": shape["lookup_probe_cap"],
            "scatter_passes": shape["scatter_passes"],
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
            "max_lookup_probes": dist["max_lookup_probes"],
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
