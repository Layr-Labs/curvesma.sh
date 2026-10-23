"""Build the binary input file for step.c (blst GT Frobenius-class walk).

Run with the py_ecc venv:  work-gt/.venv/bin/python gen_inputs.py --seed 1 --s 0.46 --out DIR

Writes DIR/inputs.bin (kernel input) and DIR/instance.json (seed, s, S, the instance
exponents for replica invariant checks). The table R_h = g^{c_h} z^{d_h} is built with
the instance exponent (a simulation shortcut); the kernel never reads that exponent.

Binary layout, all integers little-endian (u64 unless stated), Fp as 48-byte big-endian
canonical integers, Fp12 as 12 Fp in blst struct order (fp6[0].fp2[0], fp6[1].fp2[0],
...): for i in 0..5 the coefficients of w^i are A_i = (c_i + c_{i+6}) + c_{i+6}·u for the
py_ecc coefficients c_0..c_11 (u^2 = -1); struct order is A0, A2, A4, A1, A3, A5.

  header : magic, version, RA, ADDER_BITS, DP_BITS, S, W, r[4] (u256)
  consts : zeta (Fp), g0, g1 (Fp, gamma_1 = g0 + g1 u)
  fam    : 12 x (family, j, eps) for n = 0..11 (family 0 = a-value, 1 = B-value)
  table  : RA*12 entries, index h*12 + k: Fp12 T[h][k], ec[4], ed[4] (u256)
  walks  : W entries: Fp12 x, a[4], b[4] (u256)
"""
import argparse
import json
import os
import random
import struct
import sys

from py_ecc.optimized_bls12_381 import FQ12, G1, G2, curve_order, field_modulus, multiply, pairing

P = field_modulus
R = curve_order
MAGIC = int.from_bytes(b"GTRHO001", "little")
VERSION = 1
ADDER_BITS = 10
DP_BITS = 4
WALKS = 64


def fp2_mul(x, y):
    return ((x[0] * y[0] - x[1] * y[1]) % P, (x[0] * y[1] + x[1] * y[0]) % P)


def fp2_pow(x, e):
    out = (1, 0)
    while e:
        if e & 1:
            out = fp2_mul(out, x)
        x = fp2_mul(x, x)
        e >>= 1
    return out


def fp2_inv(x):
    n = pow((x[0] * x[0] + x[1] * x[1]) % P, -1, P)
    return (x[0] * n % P, -x[1] * n % P)


def u64(x):
    return struct.pack("<Q", x)


def u256(x):
    return b"".join(u64((x >> (64 * i)) & (2 ** 64 - 1)) for i in range(4))


def fp(x):
    return (x % P).to_bytes(48, "big")


STRUCT_ORDER = [0, 2, 4, 1, 3, 5]


def fp12_struct(coeffs):
    """py_ecc coefficients (12 ints) -> 12 Fp in blst struct order, as bytes."""
    out = b""
    for i in STRUCT_ORDER:
        re = (coeffs[i] + coeffs[i + 6]) % P
        im = coeffs[i + 6] % P
        out += fp(re) + fp(im)
    return out


def gamma_table():
    """gamma_n = xi^((p^n - 1)/6) in Fp2, xi = 1 + u, for n = 0..11."""
    return [fp2_pow((1, 1), (P ** n - 1) // 6) for n in range(12)]


def family_table(gam):
    """For each n: (family, j, eps) with V_n = value_family * (-1)^eps * zeta^j.

    Even n: V_n = a * gamma_n, family 0 (a = real part of the w^1 coefficient).
    Odd n:  V_n = B * s_n, family 1, B = Re(conj(A1) * gamma_1), s_n = Re(gamma_n / gamma_1).
    """
    for n in range(0, 12, 2):
        assert gam[n][1] == 0, "even gamma must lie in Fp"
    zeta = None
    for cand in [gam[2], gam[4], gam[8], gam[10]]:
        z = cand[0]
        if (z * z - z + 1) % P == 0:
            zeta = z
            break
    assert zeta is not None, "no primitive sixth root among gamma_even"
    fams = []
    for n in range(12):
        if n % 2 == 0:
            s = gam[n][0]
            fam = 0
        else:
            ratio = fp2_mul(gam[n], fp2_inv(gam[1]))
            assert ratio[1] == 0, "odd gamma ratio must lie in Fp"
            s = ratio[0]
            fam = 1
        found = None
        for j in range(3):
            for eps in (0, 1):
                v = pow(zeta, j, P) * (P - 1 if eps else 1) % P
                if v == s:
                    found = (fam, j, eps)
        assert found is not None, "gamma value is not a sixth root of unity times zeta^j"
        fams.append(found)
    for fam in (0, 1):
        pairs = {(j, e) for (f, j, e) in fams if f == fam}
        assert len(pairs) == 6, "each family must cover all six sixth roots"
    return zeta, fams


def frobenius_matrix():
    """M[row][col]: coefficients of (w^col)^p in py_ecc basis, so sigma(c) = M c."""
    cols = []
    for i in range(12):
        w_i = FQ12([0] * i + [1] + [0] * (11 - i))
        cols.append((w_i ** P).coeffs)
    return [[cols[c][r] % P for c in range(12)] for r in range(12)]


def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) % P for j in range(n)] for i in range(n)]


def mat_vec(A, v):
    return [sum(A[i][k] * v[k] for k in range(len(v))) % P for i in range(len(A))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--s", type=float, required=True, help="square probability (0 = add-only)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    RA = 1 << ADDER_BITS
    S = round(args.s * RA)
    rng = random.Random(args.seed)

    b = rng.randrange(1, R)
    y = rng.randrange(R)
    g = pairing(G2, multiply(G1, b))
    z = g ** y
    gcoef = g.coeffs
    assert (g ** R) == FQ12.one(), "g must lie in the order-r subgroup"

    gam = gamma_table()
    zeta, fams = family_table(gam)
    lam = P % R
    M = frobenius_matrix()
    powers = [[[int(i == j) for j in range(12)] for i in range(12)]]
    for _ in range(11):
        powers.append(mat_mul(powers[-1], M))

    c = [rng.randrange(R) for _ in range(RA)]
    d = [rng.randrange(R) for _ in range(RA)]
    table = []
    for h in range(RA):
        Rh = (g ** ((c[h] + y * d[h]) % R)).coeffs
        for k in range(12):
            Tk = mat_vec(powers[(12 - k) % 12], Rh)
            ec = pow(lam, -k, R) * c[h] % R
            ed = pow(lam, -k, R) * d[h] % R
            table.append((Tk, ec, ed))
    # spot-check the twist convention against exponentiation
    h0, k0 = 5 % RA, 3
    expect = (g ** (pow(lam, -k0, R) * (c[h0] + y * d[h0]) % R)).coeffs
    assert table[h0 * 12 + k0][0] == list(expect), "twisted table does not match exponentiation"

    walks = []
    for _ in range(WALKS):
        aw = rng.randrange(R)
        bw = rng.randrange(R)
        xw = (g ** ((aw + y * bw) % R)).coeffs
        walks.append((xw, aw, bw))

    head = b"".join([
        u64(MAGIC), u64(VERSION), u64(RA), u64(ADDER_BITS), u64(DP_BITS), u64(S), u64(WALKS), u256(R),
    ])
    consts = fp(zeta) + fp(gam[1][0]) + fp(gam[1][1])
    fam_bytes = b"".join(u64(f) + u64(j) + u64(e) for (f, j, e) in fams)
    os.makedirs(args.out, exist_ok=True)
    path = os.path.join(args.out, "inputs.bin")
    with open(path, "wb") as f:
        f.write(head + consts + fam_bytes)
        for Tk, ec, ed in table:
            f.write(fp12_struct(Tk) + u256(ec) + u256(ed))
        for xw, aw, bw in walks:
            f.write(fp12_struct(xw) + u256(aw) + u256(bw))
    meta = {"seed": args.seed, "s": args.s, "S": S, "RA": RA, "adder_bits": ADDER_BITS,
            "dp_bits": DP_BITS, "walks": WALKS, "y": y, "b": b, "zeta": zeta, "gamma1": gam[1],
            "families": fams}
    with open(os.path.join(args.out, "instance.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print("wrote", path, os.path.getsize(path), "bytes; S =", S, file=sys.stderr)


if __name__ == "__main__":
    main()
