"""Python reference for measurement/step.c on secp256k1: inputs, and an exact replica of its walk.

  python3 -I walk_ref.py gen <inputs.txt> [seed]
  python3 -I walk_ref.py check <inputs.txt> <batch-steps> <dump.txt>

`gen` writes NA adders [c_i]G+[d_i]Q and B start points [a_w]G+[b_w]Q for a test
target Q=[x]G (x comes from the seed; the walk never uses it). `check` replays every
walk with Python integers and compares the final points and event counts with the C
dump, so the measured binary is shown to compute the walk specified in proof.md.

Walk (proof.md section 1). Values are the reduced integers in [0,p).
  key(W)   = y if y is even, else p - y            (the class key; even)
  canon(W) = (beta^k x, (-1)^s y): k minimises beta^k x, s = y mod 2
  f(W)     = W + sigma(R_i), i = (key >> 1) mod NA, sigma = (x -> beta^-k x, y -> (-1)^s y)
The walk stores actual points; the class of f(W) is the class of canon(W) + R_i.
"""
import os
import sys

p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
r = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
BETA = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
BETA2 = BETA * BETA % p
G = (0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
     0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8)
RM = 1
# Production values; a test build overrides them as step.c's -D flags do.
ADDER_BITS = int(os.environ.get("ADDER_BITS", 16))
B = int(os.environ.get("B", 65536))
DP_BITS = int(os.environ.get("DP_BITS", 60))
CAP_PERIODS = int(os.environ.get("CAP_PERIODS", 1 << DP_BITS))
ANCHOR_PERIOD, MAX_CYCLE, SPARES = 64, 64, 1024
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
        for _ in range(NA + B + SPARES):
            S = add(mul_comb(gp, rng.randrange(r)), mul_comb(gq, rng.randrange(r)))
            out.write("%064x %064x\n" % S)


def check(path, steps, dump):
    vals = open(path).read().split()
    pts = [(int(vals[2 * i], 16), int(vals[2 * i + 1], 16)) for i in range(len(vals) // 2)]
    T, W, spares = pts[:NA], pts[NA:NA + B], pts[NA + B:NA + B + SPARES]

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
    tper = [0] * B
    two = longc = dps = abandoned = 0
    nxt = 0

    def new_trail(w):
        nonlocal nxt
        W[w] = spares[nxt % SPARES]
        nxt += 1
        prev[w] = anchor[w] = None
        tper[w] = 0

    for t in range(1, steps + 1):
        for w in reversed(range(B)):  # the kernel's order, which fixes the order spares are used
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
            kn = key(nw)
            if (kn >> (ADDER_BITS + 1)) & ((1 << DP_BITS) - 1) == 0:
                dps += 1
                new_trail(w)
                continue
            if t % ANCHOR_PERIOD == 0:
                anchor[w] = kn & LIMB
                tper[w] += 1
                if tper[w] > CAP_PERIODS:
                    abandoned += 1
                    new_trail(w)
                    continue
            prev[w] = key(W[w]) & LIMB
            W[w] = nw
    got = [tuple(int(v, 16) for v in line.split()) for line in open(dump) if line.strip()]
    for P in W:
        assert (P[1] * P[1] - P[0] ** 3 - 7) % p == 0
    mismatches = sum(1 for a, b in zip(got, W) if a != b) + abs(len(got) - B)
    print({"walks": B, "batch_steps": steps, "mismatches": mismatches,
           "two_cycles": two, "long_cycles": longc, "dps": dps, "abandoned": abandoned})
    return mismatches == 0


if __name__ == "__main__":
    if sys.argv[1] == "gen":
        gen(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    else:
        sys.exit(0 if check(sys.argv[2], int(sys.argv[3]), sys.argv[4]) else 1)
