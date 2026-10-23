"""Toy orbit coverage, plus one production generator identity.

The coverage check uses the additive group Z/rZ, r = lam^2 + lam + 1, with
multiplication by lam standing in for the GLV endomorphism. It checks the
cardinality identity and that canonical-orbit matching recovers every residue
of the explicit four-multiplier set. It does not speak to instruction costs.

The generator check is the single finite identity phi(P0) = [lambda]P0 on the
profile generator. proof.md extends it to all of G1 because both maps are
group homomorphisms.
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


def _boundary(lam):
    # Largest m with h = m-1 inside the proof constraints.
    m = min((lam + 1) // 2, (lam + 2) // 2)
    while m >= 2:
        h = m - 1
        if 2 * m - 1 <= lam and (m - 1) + (h - 1) <= lam - 1:
            return m, h
        m -= 1
    raise AssertionError("no boundary parameters")


def _orbit_min(a, mod, lam):
    if a % mod == 0:
        return 0
    best = mod
    for sigma in (1, -1, lam, -lam, 1 + lam, -(1 + lam)):
        best = min(best, (sigma * a) % mod)
    return best


def _coverage(lam):
    order = lam * lam + lam + 1
    m, h = _boundary(lam)
    babies = {}
    for a in range(m):
        babies.setdefault(_orbit_min(a, order, lam), []).append(a)
    if any(len(v) != 1 for v in babies.values()):
        raise AssertionError("baby orbits collided")
    claimed = set()
    for b in range(h):
        for a in range(m):
            for sigma in (1, -1, 1 + lam, -(1 + lam)):
                if a == 0 and sigma != 1:
                    continue
                claimed.add((b * lam + sigma * a) % order)
    expect = (m - 1) * (3 * m - 1)
    if len(claimed) != expect:
        raise AssertionError("cardinality %d != %d" % (len(claimed), expect))
    recovered = 0
    sigmas = (1, -1, lam, -lam, 1 + lam, -(1 + lam))
    for x in range(order):
        found = None
        for b in range(h):
            hit = babies.get(_orbit_min((x - b * lam) % order, order, lam))
            if not hit:
                continue
            for a in hit:
                for sigma in sigmas:
                    if (b * lam + sigma * a) % order == x:
                        found = x
        if x in claimed and found != x:
            raise AssertionError("missed claimed residue %d" % x)
        if found is not None:
            recovered += 1
    if recovered < len(claimed):
        raise AssertionError("recovered %d claimed %d" % (recovered, len(claimed)))
    return {
        "lam": lam,
        "order": order,
        "m": m,
        "h": h,
        "claimed": len(claimed),
        "recovered": recovered,
        "formula": expect,
    }


def main():
    lam = json.load(sys.stdin)["lam"]
    if not (6 <= lam <= 80):
        raise SystemExit("lam out of range")
    identity = _generator_identity()
    toy = _coverage(lam)
    print(json.dumps({
        "generator_endomorphism_identity": identity,
        "toy": toy,
        "scope": "toy cardinality plus one generator identity; not a production cost measurement",
    }))


if __name__ == "__main__":
    main()
