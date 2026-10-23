"""Python reference for measurement/step.c on BLS12-381 G1: inputs, and an exact replica of its walk.

  python3 -I walk_ref.py gen <inputs.txt> [seed]
  python3 -I walk_ref.py check <inputs.txt> <batch-steps> <dump.txt> [<kernel-stderr>]

`gen` writes NA adders [c_i]P0+[d_i]T, B start points [a_w]P0+[b_w]T and SPARES further start
points of that form, for a test target T=[x]P0 (x comes from the seed; the walk never uses it).
`check` replays every walk with Python integers and compares the final points and event counts
with the C dump, so the measured binary is shown to compute the walk specified in proof.md. With
the kernel's stderr file it also compares the kernel's own event counts, dp_digest (a digest of the
key of every distinguished point a trail records, look-ahead points included) and loops (a
diagnostic count of escapes whose key repeats one of the trail's last four escape keys).

Walk (proof.md section 1). M(v) = v*2^384 mod p is blst's Montgomery representative.
  key(W)   = M(y) if M(y) is even, else p - M(y)   (the class key; even)
  canon(W) = (beta^k x, (-1)^s y): k minimises M(beta^k x), s = M(y) mod 2
  f(W)     = W + sigma(R_i), i = (key >> 1) mod NA, sigma = (x -> beta^-k x, y -> (-1)^s y)
The walk stores actual points; the class of f(W) is the class of canon(W) + R_i.
"""
import os
import sys

p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
BETA = 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe
BETA2 = BETA * BETA % p
G = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
     0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
RM = pow(2, 384, p)
# Production values; a test build overrides them as step.c's -D flags do.
ADDER_BITS = int(os.environ.get("ADDER_BITS", 16))
B = int(os.environ.get("B", 16384))
DP_BITS = int(os.environ.get("DP_BITS", 60))
CAP_PERIODS = int(os.environ.get("CAP_PERIODS", 1 << DP_BITS))
# Filters compare the low 64 key bits; a test build narrows them (-DFILTER_MASK) to force
# false filter matches.
FMASK = (1 << int(os.environ.get("FILTER_BITS", 64))) - 1
ANCHOR_PERIOD, MAX_CYCLE, SPARES = 64, 64, 1024
NA = 1 << ADDER_BITS
LIMB = (1 << 64) - 1
# Diagnostic only: each trail keeps its last RING_K escape keys; an escape whose key repeats one is counted
# in loops (the walk is not changed).
RING_K = 4


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
            out.write("%096x %096x\n" % S)


def check(path, steps, dump, err=None):
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

    ring = [[] for _ in range(B)]

    def esc_seen(w, kk):
        # Records an escape key of walk w; counts it in loops if it repeats one of the last RING_K keys.
        nonlocal loops
        lst = ring[w]
        if kk in lst:
            loops += 1
            return
        lst.append(kk)
        if len(lst) > RING_K:
            lst.pop(0)

    prev = [None] * B
    anchor = [None] * B
    tper = [0] * B
    two = longc = dps = abandoned = ahead = loops = 0
    nxt = 0
    dig = 0  # digest of the keys of the recorded distinguished points, as the kernel's record_dp

    def new_trail(w):
        nonlocal nxt
        W[w] = spares[nxt % SPARES]
        nxt += 1
        prev[w] = anchor[w] = None
        tper[w] = 0
        ring[w] = []

    for t in range(1, steps + 1):
        for w in reversed(range(B)):  # the kernel's order, which fixes the order spares are used
            nw = f(W[w])
            kn = key(nw)
            # The distinguished-point test comes first, on every point the walk produces,
            # including look-ahead points and escape points, so an escape never skips one.
            dp = ((kn >> (ADDER_BITS + 1)) & ((1 << DP_BITS) - 1) == 0)
            esc = False
            if not dp and kn & FMASK == prev[w]:
                if key(f(nw)) == key(W[w]):
                    nw = escape([W[w], nw])
                    two += 1
                    esc = True
                    kn = key(nw)
                    dp = ((kn >> (ADDER_BITS + 1)) & ((1 << DP_BITS) - 1) == 0)
                    if not dp:
                        esc_seen(w, kn & LIMB)
            # not elif: a false predecessor match must not suppress the anchor check
            if not dp and not esc and kn & FMASK == anchor[w]:
                cyc, s = [nw], f(nw)
                while key(s) != kn and len(cyc) <= MAX_CYCLE:
                    if ((key(s) >> (ADDER_BITS + 1)) & ((1 << DP_BITS) - 1) == 0):
                        # The walk reaches s next, so s is the trail's distinguished point and its record.
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
                    dp = ((kn >> (ADDER_BITS + 1)) & ((1 << DP_BITS) - 1) == 0)
                    if not dp:
                        esc_seen(w, kn & LIMB)
            if dp:
                dps += 1
                for i in range(6):  # the six key limbs, least significant first
                    dig = (dig * 1099511628211 + ((kn >> (64 * i)) & LIMB)) & LIMB
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
        assert (P[1] * P[1] - P[0] ** 3 - 4) % p == 0
    mismatches = sum(1 for a, b in zip(got, W) if a != b) + abs(len(got) - B)
    if err is not None:
        # The kernel's last stderr line: its event counts and the record digest must agree.
        lines = [l for l in open(err).read().splitlines() if l.strip()]
        kv = dict(tok.split("=", 1) for tok in lines[-1].split())
        for name, mine in (("dps", dps), ("two_cycles", two), ("long_cycles", longc),
                           ("ahead_dps", ahead), ("abandoned", abandoned), ("loops", loops)):
            mismatches += int(kv[name]) != mine
        mismatches += kv["dp_digest"] != "%016x" % dig
    print({"walks": B, "batch_steps": steps, "mismatches": mismatches,
           "two_cycles": two, "long_cycles": longc, "dps": dps, "ahead_dps": ahead, "abandoned": abandoned,
           "loops": loops, "dp_digest": "%016x" % dig})
    return mismatches == 0


if __name__ == "__main__":
    if sys.argv[1] == "gen":
        gen(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    else:
        err = sys.argv[5] if len(sys.argv) > 5 else None
        sys.exit(0 if check(sys.argv[2], int(sys.argv[3]), sys.argv[4], err) else 1)
