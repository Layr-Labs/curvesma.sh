"""Python reference for measurement/step.c on BLS12-381 G1: inputs, and an exact replica of its walk.

  python3 -I walk_ref.py gen <inputs.txt> [seed]
  python3 -I walk_ref.py check <inputs.txt> <batch-steps> <dump.txt>

`gen` writes NA adders [c_i]P0+[d_i]T and B start points [a_w]P0+[b_w]T for a test
target T=[x]P0 (x comes from the seed; the walk never uses it). `check` replays every
walk with Python integers and compares the final points and event counts with the C
dump, so the measured binary is shown to compute the walk specified in proof.md.

Walk (proof.md section 1). M(v) = v*2^384 mod p is blst's Montgomery representative.
  key(W)   = M(y) if M(y) is even, else p - M(y)   (the class key; even)
  canon(W) = (beta^k x, (-1)^s y): k minimises M(beta^k x), s = M(y) mod 2
  f(W)     = W + sigma(R_i), i = (key >> 1) mod NA, sigma = (x -> beta^-k x, y -> (-1)^s y)
The walk stores actual points; the class of f(W) is the class of canon(W) + R_i.
"""
import sys

p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
BETA = 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe
BETA2 = BETA * BETA % p
G = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
     0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
RM = pow(2, 384, p)
ADDER_BITS, B, ANCHOR_PERIOD, MAX_CYCLE, DP_BITS = 16, 16384, 64, 64, 60
NA = 1 << ADDER_BITS
LIMB = (1 << 64) - 1


def M(v):
    return v * RM % p


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % p == 0:
            return None
        m = 3 * P[0] * P[0] * pow(2 * P[1], -1, p) % p
    else:
        m = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
    x = (m * m - P[0] - Q[0]) % p
    return (x, (m * (P[0] - x) - P[1]) % p)


def mul(k, P):
    R = None
    while k:
        if k & 1:
            R = add(R, P)
        P = add(P, P)
        k >>= 1
    return R


def comb(P):
    """rows[j][d] = [d * 256^j]P for 8-bit windows (fixed-base multiplication)."""
    rows = []
    for _ in range(32):
        row = [None, P]
        for _ in range(254):
            row.append(add(row[-1], P))
        rows.append(row)
        P = add(row[255], P)
    return rows


def mul_comb(rows, k):
    R, j = None, 0
    while k:
        R = add(R, rows[j][k & 255])
        k >>= 8
        j += 1
    return R


def canon_info(P):
    x, y = P
    bx = BETA * x % p
    ms = [M(x), M(bx), M((-x - bx) % p)]
    k = ms.index(min(ms))
    my = M(y)
    s = my & 1
    return k, s, (p - my if s else my)


def canon(P):
    k, s, _ = canon_info(P)
    return ((P[0], BETA * P[0] % p, BETA2 * P[0] % p)[k], p - P[1] if s else P[1])


def gen(path, seed):
    import random
    rng = random.Random(seed)
    gp = comb(G)
    gq = comb(mul_comb(gp, rng.randrange(1, r)))
    with open(path, "w") as out:
        for _ in range(NA + B):
            S = add(mul_comb(gp, rng.randrange(r)), mul_comb(gq, rng.randrange(r)))
            out.write("%096x %096x\n" % S)


def check(path, steps, dump):
    vals = open(path).read().split()
    pts = [(int(vals[2 * i], 16), int(vals[2 * i + 1], 16)) for i in range(len(vals) // 2)]
    T, W = pts[:NA], pts[NA:NA + B]

    def key(a):
        return canon_info(a)[2]

    def f(a):
        k, s, kk = canon_info(a)
        X, Y = T[(kk >> 1) & (NA - 1)]
        return add(a, ((X, BETA2 * X % p, BETA * X % p)[k], p - Y if s else Y))

    def escape(cyc):
        c = canon(min(cyc, key=key))
        return add(c, c)

    prev = [None] * B
    anchor = [None] * B
    two = longc = dps = 0
    for t in range(1, steps + 1):
        for w in range(B):
            nw = f(W[w])
            kn = key(nw)
            if kn & LIMB == prev[w]:
                if key(f(nw)) == key(W[w]):
                    nw = escape([W[w], nw])
                    two += 1
            elif kn & LIMB == anchor[w]:
                cyc, s = [nw], f(nw)
                while key(s) != kn and len(cyc) <= MAX_CYCLE:
                    cyc.append(s)
                    s = f(s)
                if key(s) == kn:
                    nw = escape(cyc)
                    longc += 1
            prev[w] = key(W[w]) & LIMB
            W[w] = nw
            kn = key(nw)
            if t % ANCHOR_PERIOD == 0:
                anchor[w] = kn & LIMB
            if (kn >> (ADDER_BITS + 1)) & ((1 << DP_BITS) - 1) == 0:
                dps += 1
    got = [tuple(int(v, 16) for v in line.split()) for line in open(dump) if line.strip()]
    for P in W:
        assert (P[1] * P[1] - P[0] ** 3 - 4) % p == 0
    mismatches = sum(1 for a, b in zip(got, W) if a != b) + abs(len(got) - B)
    print({"walks": B, "batch_steps": steps, "mismatches": mismatches,
           "two_cycles": two, "long_cycles": longc, "dps": dps})
    return mismatches == 0


if __name__ == "__main__":
    if sys.argv[1] == "gen":
        gen(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    else:
        sys.exit(0 if check(sys.argv[2], int(sys.argv[3]), sys.argv[4]) else 1)
