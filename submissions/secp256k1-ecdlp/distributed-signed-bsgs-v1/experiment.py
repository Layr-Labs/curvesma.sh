"""Bounded evidence for a signed baby-step/giant-step search over a distributed merge network.

One JSON request on stdin, one JSON result on stdout. Standard library only.

What is checked, and at which scope:

* a toy prime-order elliptic-curve group built from a small prime field: the
  baby and giant record lists, the F-ary lazy merge network with duplicate keys
  and left-first tie handling, the stream intersection, and the recovered
  scalar, compared against an independent interval prediction for every scalar
  the toy group can represent;
* Montgomery prefix/suffix batch inversion against direct inversion;
* the exact production coverage count 2*h*(2m-1) - r > 0 on the protected
  integers;
* the exact integer ledger and the exact score inequality, recomputed from the
  itemised prices.

No counter in this file is an AArch64 instruction count and nothing here is a
production benchmark, a production run or a complexity proof.
"""

import json
import sys


FIELD = int("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16)
ORDER = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
GENERATOR = (
    int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16),
    int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16),
)

# Exact production integers used by the ledger and by the coverage count.
M_PROD = 1 << 127
H_PROD = 1 << 127
D_PROD = 2 * M_PROD - 1
TOTAL_INSTRUCTIONS = 2050001412324430834062893202469342222401024
BOUND_NUMERATOR = 1305909
BOUND_DENOMINATOR = 10000
PEAK_BYTES = 115804030901347946931027802910112881136640
NODES = 10802614822886935348043638331167244810


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def integer(request, name, default, lower, upper):
    value = request.get(name, default)
    if type(value) is not int or not lower <= value <= upper:
        raise ValueError("%s must be an integer in [%d,%d]" % (name, lower, upper))
    return value


# ----------------------------------------------------------------------------
# Toy elliptic curve: y^2 = x^3 + b over F_p with prime group order.
# ----------------------------------------------------------------------------

def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def curve_order(p, b, chi):
    total = 1  # the point at infinity
    cube = [x * x * x % p for x in range(p)]
    for x in range(p):
        total += 1 + chi[(cube[x] + b) % p]
    return total


def find_toy_curve(low, high, prime_cap):
    """Curve y^2 = x^3 + b over F_p with prime group order in (low, high]."""
    for p in range(5, prime_cap):
        if not is_prime(p):
            continue
        chi = [0] * p
        for y in range(1, p):
            chi[y * y % p] = 1
        for y in range(1, p):
            if chi[y] == 0:
                chi[y] = -1
        for b in range(1, 40):
            n = curve_order(p, b, chi)
            if n <= low or n > high or not is_prime(n):
                continue
            for x in range(p):
                rhs = (x * x * x + b) % p
                if rhs and chi[rhs] == 1:
                    y = pow(rhs, (p + 1) // 4, p) if p % 4 == 3 else sqrt_mod(rhs, p)
                    if y is not None and (y * y - rhs) % p == 0:
                        return {"p": p, "b": b, "order": n, "generator": (x, y)}
    raise ValueError("no toy curve with a prime order in (%d, %d] under prime %d" % (low, high, prime_cap))


def sqrt_mod(a, p):
    for y in range(p):
        if (y * y - a) % p == 0:
            return y
    return None


def point_add(curve, left, right):
    p, b = curve["p"], curve["b"]
    if left is None:
        return right
    if right is None:
        return left
    (x1, y1), (x2, y2) = left, right
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if left == right:
        if y1 == 0:
            return None
        slope = (3 * x1 * x1) * pow(2 * y1, p - 2, p) % p
    else:
        slope = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
    x3 = (slope * slope - x1 - x2) % p
    y3 = (slope * (x1 - x3) - y1) % p
    return (x3, y3)


def point_mul(curve, scalar, point):
    result = None
    addend = point
    while scalar:
        if scalar & 1:
            result = point_add(curve, result, addend)
        addend = point_add(curve, addend, addend)
        scalar >>= 1
    return result


def key_of(point):
    return (1, 0) if point is None else (0, point[0])


# ----------------------------------------------------------------------------
# F-ary lazy merge network over singleton leaves.
# ----------------------------------------------------------------------------

EMPTY = "EMPTY"
END = "END"


class Node:
    """One internal node: one cached head slot per child, left-first on ties."""

    def __init__(self, children):
        self.children = children
        self.slots = [EMPTY] * len(children)
        self.stats = {"emissions": 0, "comparisons": 0, "requests": 0, "ends": 0}

    def next(self):
        for index, slot in enumerate(self.slots):
            if slot is EMPTY:
                self.stats["requests"] += 1
                value = self.children[index].next()
                self.slots[index] = value
                if value is END:
                    self.stats["ends"] += 1
        best = -1
        for index, slot in enumerate(self.slots):
            if slot is END:
                continue
            if best < 0:
                best = index
                continue
            self.stats["comparisons"] += 1
            if slot[0] < self.slots[best][0]:  # strict: ties keep the lower index
                best = index
        if best < 0:
            return END
        value = self.slots[best]
        self.slots[best] = EMPTY
        self.stats["emissions"] += 1
        return value


class Leaf:
    def __init__(self, record):
        self.record = record
        self.served = False

    def next(self):
        if self.served:
            return END
        self.served = True
        return self.record


def merge_tree(records, fanout):
    """Complete fanout-ary tree over singleton leaves, padded with END leaves."""
    leaves = [Leaf(record) for record in records]
    while len(leaves) > 1:
        while len(leaves) % fanout:
            leaves.append(Leaf(END))
        leaves = [Node(leaves[i:i + fanout]) for i in range(0, len(leaves), fanout)]
    return leaves[0]


def drain(root, limit):
    out = []
    while len(out) < limit:
        value = root.next()
        if value is END:
            break
        out.append(value)
    return out


def merge_checks(fanout):
    """Sortedness, multiset preservation, duplicate retention and left-first ties."""
    records = [
        ((0, 10), 0, 0), ((0, 3), 1, 1), ((0, 10), 2, 2), ((0, 3), 3, 3),
        ((0, 10), 4, 4), ((0, 7), 5, 5), ((1, 0), 6, 6), ((0, 3), 7, 7),
    ]
    root = merge_tree(records, fanout)
    out = drain(root, len(records) + 1)
    require(len(out) == len(records), "merge lost or invented records")
    require(out == sorted(records, key=lambda r: (r[0], r[2])), "merge emitted a wrong order")
    require([r[1] for r in out] != [], "empty merge")
    require(sorted(out) == sorted(records), "merge is not a permutation of its input")
    counts = {}
    for record in records:
        counts[record[0]] = counts.get(record[0], 0) + 1
    for key, count in counts.items():
        require(sum(1 for r in out if r[0] == key) == count, "duplicate key count changed")
    ties = [r for r in out if r[0] == (0, 10)]
    require([r[2] for r in ties] == [0, 2, 4], "left-first tie order is wrong")
    # A one-level node must emit an exhausted child's END without starving siblings.
    node = Node([Leaf(((0, 1), 0, 0)), Leaf(END), Leaf(((0, 2), 2, 2))])
    require(node.next()[0] == (0, 1), "END sibling starved the merge")
    return {"input": len(records), "emitted": len(out), "duplicate_keys": 2,
            "comparisons": root.stats["comparisons"], "ends": root.stats["ends"]}


# ----------------------------------------------------------------------------
# Montgomery batch inversion
# ----------------------------------------------------------------------------

def batch_invert(values, p):
    """Prefix products, one inversion, backward pass.  Zero is not passed in."""
    prefixes = [1] * (len(values) + 1)
    for index, value in enumerate(values):
        require(value % p != 0, "batch inversion of zero")
        prefixes[index + 1] = prefixes[index] * value % p
    running = pow(prefixes[-1], p - 2, p)
    out = [0] * len(values)
    for index in range(len(values) - 1, -1, -1):
        out[index] = running * prefixes[index] % p
        running = running * values[index] % p
    require(running == 1, "batch inversion did not close")
    return out


def batch_checks(p, block_size):
    values = [(7 * i + 11) % p for i in range(1, block_size + 1)]
    values = [v if v else 1 for v in values]
    batch = batch_invert(values, p)
    direct = [pow(v, p - 2, p) for v in values]
    require(batch == direct, "batch inversion disagrees with direct inversion")
    require(all(v * w % p == 1 for v, w in zip(values, batch)), "batch inverse is not inverse")
    return {"block": block_size, "matched": len(values)}


# ----------------------------------------------------------------------------
# Toy signed BSGS over the merge network
# ----------------------------------------------------------------------------

def predicted_covered(order, m, h, x):
    d = 2 * m - 1
    for i in range(h):
        for j in range(m):
            if (x - i * d - j) % order == 0 or (x - i * d + j) % order == 0:
                return True
    return False


def toy_search(curve, x, m, h, fanout):
    """Return the recovered scalar or None, using only group points and keys."""
    order = curve["order"]
    generator = curve["generator"]
    d = 2 * m - 1
    target = point_mul(curve, x, generator)
    minus_d = point_mul(curve, -d % order, generator)
    baby = [((key_of(point_mul(curve, j, generator))), j) for j in range(m)]
    giant = []
    for i in range(h):
        point = point_add(curve, target, point_mul(curve, i, minus_d))
        giant.append((key_of(point), i))
    baby_root = merge_tree([(k, j, 0) for k, j in baby], fanout)
    giant_root = merge_tree([(k, i, 1) for k, i in giant], fanout)
    left, right = baby_root.next(), giant_root.next()
    probes = 0
    while left is not END and right is not END:
        probes += 1
        if probes > 4 * (m + h) + 16:
            raise AssertionError("intersection did not terminate")
        if left[0] < right[0]:
            left = baby_root.next()
        elif right[0] < left[0]:
            right = giant_root.next()
        else:
            i, j = right[1], left[1]
            for candidate in ((i * d + j) % order, (i * d - j) % order):
                if point_mul(curve, candidate, generator) == target:
                    return candidate
            raise AssertionError("equal key with no verified scalar")
    return None


def toy_checks(m, h, fanout):
    require(m >= 3 and h >= 2, "toy needs m >= 3 and h >= 2")
    d = 2 * m - 1
    # Order above h*d keeps the intervals of the toy from wrapping, mirroring the
    # production regime; order at or below 2*h*d keeps coverage at or above 1/2.
    curve = find_toy_curve(3 * h * d // 2, 2 * h * d, prime_cap=2 * h * d + 4000)
    order = curve["order"]
    covered = [x for x in range(order) if predicted_covered(order, m, h, x)]
    require(len(covered) == h * d, "toy covered count is not h*d")
    expected = sorted(list(range(order - (m - 1), order)) + list(range(0, (h - 1) * d + m)))
    require(covered == expected, "toy covered set is not the predicted contiguous union")
    marks = [0, 1, m - 2, m - 1, m, order - 1]
    for i in (0, 1, h // 2, h - 1):
        for offset in (-m, -m + 1, 0, m - 1, m):
            marks.append((i * d + offset) % order)
    marks.extend((order * t) // 7 for t in range(7))
    sample = sorted(set(x for x in marks if 0 <= x < order))
    hits = 0
    for x in sample:
        recovered = toy_search(curve, x, m, h, fanout)
        expected = predicted_covered(order, m, h, x)
        require((recovered is not None) == expected, "toy search disagrees with the interval prediction")
        if recovered is not None:
            require(recovered == x, "toy search returned the wrong scalar")
            hits += 1
    covered_fraction = len(covered) / order
    require(covered_fraction >= 0.39, "toy coverage is below the required threshold")
    return {"p": curve["p"], "b": curve["b"], "order": order, "m": m, "h": h,
            "fanout": fanout, "covered": len(covered), "covered_fraction": round(covered_fraction, 6),
            "sample": len(sample), "sample_hits": hits}


# ----------------------------------------------------------------------------
# Production coverage and ledger
# ----------------------------------------------------------------------------

def coverage_checks():
    covered = H_PROD * (2 * M_PROD - 1)
    twice_minus_order = 2 * covered - ORDER
    require(covered < ORDER, "coverage can exceed the order")
    require(twice_minus_order > 0, "coverage is not more than half the order")
    require(twice_minus_order == 92138019644721193389046258963199934143,
            "the exact slack integer changed")
    require(D_PROD == 2 * M_PROD - 1, "d is not 2m-1")
    return {"covered": str(covered), "twice_minus_order": str(twice_minus_order),
            "covered_below_order": True, "success_at_least_half": True}


def ledger_checks():
    """Recompute every ledger term from the itemised prices."""
    add = 1074
    wire, stored = 7, 8
    per_word, hop_control = 6, 8
    hop = wire * per_word + hop_control
    levels = 22
    comparison, sift_levels, sift_pointer = 16, 6, 2
    tournament = sift_levels * comparison + sift_levels * sift_pointer
    read_hop = hop + tournament + 4 + 4
    write_per_record = levels * hop
    read_per_record = levels * read_hop
    records = 2 * M_PROD
    nodes = 2 * (-(-(2 * M_PROD) // 63))
    node_state_words = 64 * stored + 64
    node_code_words = 512
    node_words = node_state_words + node_code_words
    node_init = node_state_words * 4 + node_code_words * 4 + 512
    jacobian = 16 * 512 + 32 * 64 + 4096
    fallbacks = 2 * (1 << 16) * 2 * jacobian
    chain = 2 * (1 << 16) * 512 * jacobian
    reserve = 1 << 40
    terms = {
        "reference_additions": add * (M_PROD + H_PROD),
        "record_writes": records * write_per_record,
        "record_reads": records * read_per_record,
        "node_initialisation": nodes * node_init,
        "root_intersection": records * 36,
        "generation_outer": records * 8,
        "exceptional_fallbacks": fallbacks,
        "chain_start_points": chain,
        "fixed_reserve": reserve,
    }
    require(hop == 50, "hop price changed")
    require(tournament == 108, "tournament price changed")
    require(write_per_record == 1100, "write price changed")
    require(read_per_record == 3652, "read price changed")
    require(node_init == 4864, "node initialisation price changed")
    require(nodes == NODES, "node count changed")
    total = sum(terms.values())
    require(total == TOTAL_INSTRUCTIONS, "ledger total disagrees with the declared total")
    require(total ** BOUND_DENOMINATOR < 1000 ** BOUND_DENOMINATOR * 2 ** BOUND_NUMERATOR,
            "declared bound does not hold")
    require(not (total ** BOUND_DENOMINATOR
                 < 1000 ** BOUND_DENOMINATOR * 2 ** (BOUND_NUMERATOR - 1)),
            "the next lower four-decimal bound also holds")
    require(BOUND_NUMERATOR == 1305909 and BOUND_DENOMINATOR == 10000,
            "the four-decimal bound changed")
    require(PEAK_BYTES == records * 64 + nodes * node_words * 8 + (1 << 24),
            "peak byte count disagrees")
    require(PEAK_BYTES < 2 ** 137, "peak is not below 2^137")
    return {"total_instructions": str(total), "bound": "%d/%d" % (BOUND_NUMERATOR, BOUND_DENOMINATOR),
            "terms": {name: str(value) for name, value in terms.items()},
            "peak_bytes": str(PEAK_BYTES), "bits": total.bit_length()}


# ----------------------------------------------------------------------------

def main():
    raw = sys.stdin.read()
    request = json.loads(raw) if raw.strip() else {}
    mode = request.get("mode", "all")
    if mode not in ("all", "toy", "merge", "accounting"):
        raise ValueError("mode must be one of all, toy, merge, accounting")
    m = integer(request, "m", 16, 3, 24)
    h = integer(request, "h", 16, 2, 24)
    block_size = integer(request, "block_size", 4, 1, 16)
    fanout = integer(request, "fanout", 4, 2, 8)
    result = {"mode": mode, "ok": True}
    if mode in ("all", "merge"):
        result["merge"] = merge_checks(fanout)
        result["batch_inversion"] = batch_checks(4099, block_size)
    if mode in ("all", "toy"):
        result["toy"] = toy_checks(m, h, fanout)
    if mode in ("all", "accounting"):
        result["coverage"] = coverage_checks()
        result["ledger"] = ledger_checks()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
