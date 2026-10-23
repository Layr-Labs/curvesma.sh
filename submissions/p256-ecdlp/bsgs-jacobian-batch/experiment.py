"""Toy-scale run of the production P-256 BSGS update.

The same Jacobian mixed-addition walk, the same block batch inversion
(prefix products, one inversion, back substitution), the same open-addressed
lookup with the same probe cap and the same coverage loop as the production
claim, with m = ceil(sqrt(n)) instead of 2^128.  It is toy-scope evidence for
the update and the coverage argument only; it is not evidence for the
production instruction counts.
"""
import json
import math
import sys

P_MOD = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
CURVE_B = 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B
GX = 0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296
GY = 0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5
ORDER = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
BLOCK = 8
PROBES = 64
LOAD_SHIFT = 4
HASH_MULT = 0x9E3779B97F4A7C15

INF = (0, 1, 0)


def jdbl(pt):
    x, y, z = pt
    if z == 0 or y == 0:
        return INF
    delta = z * z % P_MOD
    gamma = y * y % P_MOD
    beta = x * gamma % P_MOD
    alpha = 3 * (x - delta) * (x + delta) % P_MOD
    x3 = (alpha * alpha - 8 * beta) % P_MOD
    z3 = ((y + z) * (y + z) - gamma - delta) % P_MOD
    y3 = (alpha * (4 * beta - x3) - 8 * gamma * gamma) % P_MOD
    return (x3, y3, z3)


def jadd_mixed(pt, ax, ay):
    """add-2007-bl: Jacobian + affine, no inversion."""
    x1, y1, z1 = pt
    if z1 == 0:
        return (ax, ay, 1)
    z1z1 = z1 * z1 % P_MOD
    u2 = ax * z1z1 % P_MOD
    s2 = ay * z1 % P_MOD * z1z1 % P_MOD
    h = (u2 - x1) % P_MOD
    rr = 2 * (s2 - y1) % P_MOD
    if h == 0:
        if rr == 0:
            return jdbl(pt)
        return INF
    hh = h * h % P_MOD
    i4 = 4 * hh % P_MOD
    j = h * i4 % P_MOD
    v = x1 * i4 % P_MOD
    x3 = (rr * rr - j - 2 * v) % P_MOD
    y3 = (rr * (v - x3) - 2 * y1 * j) % P_MOD
    z3 = ((z1 + h) * (z1 + h) - z1z1 - hh) % P_MOD
    return (x3, y3, z3)


def batch_invert(vals):
    n = len(vals)
    pre = [1] * (n + 1)
    for i, v in enumerate(vals):
        pre[i + 1] = pre[i] * v % P_MOD
    inv = pow(pre[n], P_MOD - 2, P_MOD)
    out = [0] * n
    for i in range(n - 1, -1, -1):
        out[i] = inv * pre[i] % P_MOD
        inv = inv * vals[i] % P_MOD
    return out


def normalize(buf):
    """One batch inversion for the whole block; None marks infinity."""
    zs = [z if z else 1 for (_, _, z) in buf]
    invs = batch_invert(zs)
    out = []
    for (x, y, z), zi in zip(buf, invs):
        if z == 0:
            out.append(None)
            continue
        zi2 = zi * zi % P_MOD
        out.append((x * zi2 % P_MOD, y * zi2 * zi % P_MOD))
    return out


def blocks(start, ax, ay, count):
    cur = start
    buf = []
    for _ in range(count):
        buf.append(cur)
        cur = jadd_mixed(cur, ax, ay)
        if len(buf) == BLOCK:
            yield normalize(buf)
            buf = []
    if buf:
        yield normalize(buf)


def slot_of(x, shift):
    return ((x * HASH_MULT) & 0xFFFFFFFFFFFFFFFF) >> (64 - shift)


def insert(table, mask, x, y, idx):
    pos = slot_of(x, mask.bit_length())
    for step in range(PROBES):
        cell = table[pos]
        if cell is None:
            table[pos] = (x, y, idx)
            return True
        if cell[0] == x and cell[1] == y:
            return True
        pos = (pos + 1) & mask
    return False


def lookup(table, mask, x, y):
    pos = slot_of(x, mask.bit_length())
    for step in range(PROBES):
        cell = table[pos]
        if cell is None:
            return None
        if cell[0] == x and cell[1] == y:
            return cell[2]
        pos = (pos + 1) & mask
    return None


def scalar_mul(k, pt):
    acc = INF
    for bit in bin(k)[2:]:
        acc = jdbl(acc)
        if bit == "1":
            acc = jadd_mixed(acc, pt[0], pt[1])
    return acc


def to_affine(pt):
    x, y, z = pt
    if z == 0:
        return None
    zi = pow(z, P_MOD - 2, P_MOD)
    zi2 = zi * zi % P_MOD
    return (x * zi2 % P_MOD, y * zi2 * zi % P_MOD)


def recover(qx, qy, m, h):
    slots = 1
    while slots < (m << LOAD_SHIFT):
        slots <<= 1
    mask = slots - 1
    table = [None] * slots
    base = 0
    for blk in blocks(INF, GX, GY, m):
        for off, pt in enumerate(blk):
            idx = base + off
            if idx == 0:
                continue
            if not insert(table, mask, pt[0], pt[1], idx):
                raise AssertionError("probe cap exceeded")
        base += len(blk)
    gm = scalar_mul(m, (GX, GY))
    gma = to_affine(gm)
    dx, dy = gma[0], (-gma[1]) % P_MOD
    start = (qx, qy, 1) if qy is not None else INF
    i = 0
    for blk in blocks(start, dx, dy, h):
        for off, pt in enumerate(blk):
            step = i + off
            if pt is None:
                return (step * m) % ORDER
            j = lookup(table, mask, pt[0], pt[1])
            if j is not None:
                return (step * m + j) % ORDER
        i += len(blk)
    return None


def main():
    n = json.load(sys.stdin)["scalars"]
    if not isinstance(n, int) or n < 2 or n > 4096:
        raise SystemExit("scalars must be an integer in [2, 4096]")
    if pow(GY, 2, P_MOD) != (GX * GX % P_MOD * GX - 3 * GX + CURVE_B) % P_MOD:
        raise AssertionError("generator is not on the curve")
    if to_affine(scalar_mul(ORDER, (GX, GY))) is not None:
        raise AssertionError("[r]P is not infinity")
    m = math.isqrt(n) + 1
    h = (n + m - 1) // m
    for x in range(n):
        q = to_affine(scalar_mul(x, (GX, GY)))
        qx, qy = (0, None) if q is None else q
        got = recover(qx, qy, m, h)
        if got != x:
            raise AssertionError("scalar %d recovered as %r" % (x, got))
    print(json.dumps({
        "curve": "NIST P-256 (secp256r1), a=-3",
        "update": "Jacobian mixed addition, block batch inversion, 64-probe lookup",
        "checked": n,
        "m": m,
        "h": h,
        "all_recovered": True,
        "scope": "toy size on the production curve; not a production cost bound",
    }))


main()
