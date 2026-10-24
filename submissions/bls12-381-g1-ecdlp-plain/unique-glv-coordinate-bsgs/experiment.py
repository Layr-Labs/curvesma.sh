"""Toy rectangle coverage, plus one production generator identity.

Cardinality uses Z/rZ with r = lam^2 + lam + 1. Multiplication by lam stands
in for the GLV endomorphism. Every feasible rectangle 1 <= h <= m-1 with
2*m-1 <= lam is built as four explicit blocks and checked against
h*(4*m-h-2). Orbit-minimum recovery is checked on the h ~ (2/3)m slice.
Neither check is a production cost measurement.

The generator check is phi(P0) = [lambda]P0 on the profile generator.
proof.md extends it to G1 because both maps are group homomorphisms.
"""

import json
import sys

P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
LAM = 0xAC45A4010001A40200000000FFFFFFFF
BETA = 0x1A0111EA397FE699EC02408663D4DE85AA0D857D89759AD4897D29650FB85F9B409427EB4F49FFFD8BFD00000000AAAC
P0 = (
    0x17F1D3A73197D7942695638C4FA9AC0FC3688C4F9774B905A14E3A3F171BAC586C55E83FF97A1AEFFB3AF00ADB22C6BB,
    0x08B3F481E3AAA0F1A09E30ED741D8AE4FCF5E095D5D00AF600DB18CB2C04B3EDD03CC744A2888AE40CAA232946C5E7E1,
)


def _add(p, q):
    if p is None:
        return q
    if q is None:
        return p
    x1, y1 = p
    x2, y2 = q
    if x1 == x2 and (y1 + y2) % P == 0:
        return None
    if p != q:
        slope = ((y2 - y1) * pow(x2 - x1, -1, P)) % P
    else:
        slope = ((3 * x1 * x1) * pow(2 * y1, -1, P)) % P
    x3 = (slope * slope - x1 - x2) % P
    y3 = (slope * (x1 - x3) - y1) % P
    return x3, y3


def _mul(k, point):
    acc = None
    base = point
    while k:
        if k & 1:
            acc = _add(acc, base)
        base = _add(base, base)
        k >>= 1
    return acc


def _generator_identity():
    assert (BETA * BETA + BETA + 1) % P == 0
    assert LAM * LAM + LAM + 1 == 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
    image = ((BETA * P0[0]) % P, P0[1])
    scaled = _mul(LAM, P0)
    if image != scaled:
        raise AssertionError("phi(P0) != [lambda]P0")
    return True


def _orbit_min(a, mod, lam):
    if a % mod == 0:
        return 0
    best = mod
    for sigma in (1, -1, lam, -lam, 1 + lam, -(1 + lam)):
        best = min(best, (sigma * a) % mod)
    return best


def _blocks(m, h, lam, order):
    """Four blocks. A repeated residue is a failed disjointness check."""
    seen = {}

    def put(tag, residue, meta):
        residue %= order
        prev = seen.get(residue)
        if prev is not None:
            raise AssertionError(
                "overlap %s %s vs %s %s at %d" % (prev[0], prev[1], tag, meta, residue)
            )
        seen[residue] = (tag, meta)

    for b in range(h):
        for a in range(m):
            put("P+", a + b * lam, (a, b))
    for b in range(h):
        for a in range(1, m):
            put("P-", -a + b * lam, (a, b))
    for a in range(1, m):
        for b in range(h):
            if a + b < h:
                continue
            put("Q+", (1 + lam) * a + b * lam, (a, b))
    for a in range(1, m):
        limit = a if a < h else h
        for b in range(limit):
            put("Q-", -(1 + lam) * a + b * lam, (a, b))
    return seen


def _recover_all(seen, m, h, lam, order):
    babies = {}
    for a in range(m):
        babies.setdefault(_orbit_min(a, order, lam), []).append(a)
    if any(len(v) != 1 for v in babies.values()):
        raise AssertionError("baby orbits collided")
    units = (1, -1, lam, -lam, 1 + lam, -(1 + lam))
    hit = 0
    for target in seen:
        found = False
        for b in range(h):
            key = _orbit_min((target - b * lam) % order, order, lam)
            for a in babies.get(key, ()):
                for sigma in units:
                    if (b * lam + sigma * a) % order == target:
                        found = True
                        break
                if found:
                    break
            if found:
                break
        if not found:
            raise AssertionError("unrecovered residue")
        hit += 1
    return hit


def _coverage(lam):
    order = lam * lam + lam + 1
    m_max = (lam + 1) // 2
    checked = 0
    recovered = 0
    densest = None
    for m in range(2, m_max + 1):
        if 2 * m - 1 > lam:
            break
        for h in range(1, m):
            if (m - 1) + (h - 1) > lam - 1:
                break
            seen = _blocks(m, h, lam, order)
            formula = h * (4 * m - h - 2)
            if len(seen) != formula:
                raise AssertionError(
                    "|S|=%d formula=%d m=%d h=%d" % (len(seen), formula, m, h)
                )
            checked += 1
            # Production band: one recovery sweep per m, at h = min(m-1, floor(2m/3)).
            if h == min(m - 1, max(1, (2 * m) // 3)):
                recovered += _recover_all(seen, m, h, lam, order)
                densest = {"m": m, "h": h, "cardinality": len(seen)}
    if checked < 1 or densest is None:
        raise AssertionError("no rectangle")
    return {
        "lam": lam,
        "order": order,
        "rectangles": checked,
        "recovered_residues": recovered,
        "two_thirds": densest,
    }


def main():
    lam = json.load(sys.stdin)["lam"]
    if not isinstance(lam, int) or isinstance(lam, bool) or not (6 <= lam <= 80):
        raise SystemExit("lam out of range")
    identity = _generator_identity()
    toy = _coverage(lam)
    print(json.dumps({
        "generator_endomorphism_identity": identity,
        "toy": toy,
        "scope": "toy rectangle cardinality plus one generator identity; not a production cost measurement",
    }))


if __name__ == "__main__":
    main()
