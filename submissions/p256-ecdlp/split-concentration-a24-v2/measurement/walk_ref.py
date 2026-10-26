"""Python reference for measurement/step.c on P-256: inputs, and an exact replica of its walk.

  python3 -I walk_ref.py gen <inputs.txt> [seed]
  python3 -I walk_ref.py check <inputs.txt> <batch-steps> <dump.txt> [kernel-stderr]

`gen` writes NA adders [c_i]P+[d_i]Q, B start points [a_w]P+[b_w]Q and SPARES further start
points, all of that form, for a test target Q=[x]P (x comes from the seed; the walk never uses it). `check` replays every
walk with Python integers and compares the final points and event counts with the C
dump, so the measured binary is shown to compute the walk specified in proof.md. With the
kernel's stderr file it also compares the event counts the kernel printed and dp_digest, a
digest of the class key of every distinguished point recorded (look-ahead ones included):
the four 64-bit limbs of each key, least significant first, folded by h = h*1099511628211 + limb.

Walk (proof.md section 1). M(v) = v*2^256 mod p is the kernel's Montgomery representative.
  key(W)   = M(x)                                  (the class key of {W, -W})
  canon(W) = (x, (-1)^s y), s = M(y) mod 2
  f(W)     = W + sigma(R_i), i = key mod NA, sigma = (y -> (-1)^s y)
The walk stores actual points; the class of f(W) is the class of canon(W) + R_i.
A fruitless cycle is left by 2*canon(W), W the member of least key. The last RING_K escape keys of a
trail are kept for a diagnostic count (loops) that changes nothing in the walk.
"""
import os
import sys

p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
r = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
A = p - 3
G = (0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,
     0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)
RM = pow(2, 256, p)
ADDER_BITS = int(os.environ.get("ADDER_BITS", 16))
B = int(os.environ.get("B", 65536))
DP_BITS = int(os.environ.get("DP_BITS", 60))
CAP_PERIODS = int(os.environ.get("CAP_PERIODS", 1 << DP_BITS))
FMASK = (1 << int(os.environ.get("FILTER_BITS", 64))) - 1
RING_K = 4
ANCHOR_PERIOD, MAX_CYCLE, SPARES = 81, 64, 1024
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
        m = (3 * P[0] * P[0] + A) * pow(2 * P[1], -1, p) % p
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
    """(k, s, key) as in the BLS replica; k is always 0 here."""
    return 0, M(P[1]) & 1, M(P[0])


def canon(P):
    _, s, _ = canon_info(P)
    return (P[0], p - P[1] if s else P[1])


def gen(path, seed):
    import random
    rng = random.Random(seed)
    gp = comb(G)
    gq = comb(mul_comb(gp, rng.randrange(1, r)))
    with open(path, "w") as out:
        for _ in range(NA + B + SPARES):
            S = add(mul_comb(gp, rng.randrange(r)), mul_comb(gq, rng.randrange(r)))
            out.write("%064x %064x\n" % S)


def check(path, steps, dump, kernel_err=None):
    vals = open(path).read().split()
    pts = [(int(vals[2 * i], 16), int(vals[2 * i + 1], 16)) for i in range(len(vals) // 2)]
    T, W, spares = pts[:NA], pts[NA:NA + B], pts[NA + B:NA + B + SPARES]

    def key(a):
        return canon_info(a)[2]

    def f(a):
        _, s, kk = canon_info(a)
        X, Y = T[kk & (NA - 1)]
        return add(a, (X, p - Y if s else Y))

    def escape(cyc):
        c = canon(min(cyc, key=key))
        return add(c, c)

    ring = [[] for _ in range(B)]

    def esc_seen(w, kk):
        lst = ring[w]
        if kk in lst:
            return True
        lst.append(kk)
        if len(lst) > RING_K:
            lst.pop(0)
        return False

    prev = [None] * B
    anchor = [None] * B
    tper = [0] * B
    two = longc = dps = abandoned = ahead = loops = 0
    nxt = 0
    dig = 0

    def new_trail(w):
        nonlocal nxt
        W[w] = spares[nxt % SPARES]
        nxt += 1
        prev[w] = anchor[w] = None
        tper[w] = 0
        ring[w] = []

    for t in range(1, steps + 1):
        for w in reversed(range(B)):
            nw = f(W[w])
            kn = key(nw)
            dp = ((kn >> ADDER_BITS) & ((1 << DP_BITS) - 1) == 0)
            esc = False
            if not dp and kn & FMASK == prev[w]:
                if key(f(nw)) == key(W[w]):
                    nw = escape([W[w], nw])
                    two += 1
                    esc = True
                    kn = key(nw)
                    dp = ((kn >> ADDER_BITS) & ((1 << DP_BITS) - 1) == 0)
                    if not dp and esc_seen(w, kn & LIMB):
                        loops += 1
            if not dp and not esc and kn & FMASK == anchor[w]:
                cyc, s = [nw], f(nw)
                while key(s) != kn and len(cyc) <= MAX_CYCLE:
                    if ((key(s) >> ADDER_BITS) & ((1 << DP_BITS) - 1) == 0):
                        nw = s
                        kn = key(nw)
                        dp = True
                        ahead += 1
                        break
                    cyc.append(s)
                    s = f(s)
                if not dp and key(s) == kn:
                    nw = escape(cyc)
                    longc += 1
                    kn = key(nw)
                    dp = ((kn >> ADDER_BITS) & ((1 << DP_BITS) - 1) == 0)
                    if not dp and esc_seen(w, kn & LIMB):
                        loops += 1
            if dp:
                dps += 1
                k = key(nw)
                for i in range(4):
                    dig = (dig * 1099511628211 + ((k >> (64 * i)) & LIMB)) & LIMB
                new_trail(w)
                continue
            if t % ANCHOR_PERIOD == 0:
                anchor[w] = kn & FMASK
                tper[w] += 1
                if tper[w] > CAP_PERIODS:
                    abandoned += 1
                    new_trail(w)
                    continue
            prev[w] = key(W[w]) & FMASK
            W[w] = nw
    got = [tuple(int(v, 16) for v in line.split()) for line in open(dump) if line.strip()]
    for P in W:
        assert (P[1] * P[1] - P[0] ** 3 - A * P[0] - 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b) % p == 0
    mismatches = sum(1 for a, b in zip(got, W) if a != b) + abs(len(got) - B)
    mine = {"two_cycles": two, "long_cycles": longc, "dps": dps, "ahead_dps": ahead, "abandoned": abandoned, "loops": loops}
    if kernel_err:
        lines = [l for l in open(kernel_err).read().splitlines() if l.strip()]
        kern = dict(tok.partition("=")[::2] for tok in lines[-1].split())
        for name, val in mine.items():
            if int(kern.get(name, -1)) != val:
                mismatches += 1
        if kern.get("dp_digest") != "%016x" % dig:
            mismatches += 1
    print({"walks": B, "batch_steps": steps, "mismatches": mismatches, **mine, "dp_digest": "%016x" % dig})
    return mismatches == 0


if __name__ == "__main__":
    if sys.argv[1] == "gen":
        gen(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    else:
        sys.exit(0 if check(sys.argv[2], int(sys.argv[3]), sys.argv[4],
                            sys.argv[5] if len(sys.argv) > 5 else None) else 1)
