"""Check Lemma 1 of proof.md: the 64-bit fraction approximations used by step.c.

Input (JSON on stdin, or a seed as the only argument): {"seed": int, "samples": int}.

For a = M(x) and C_0 = floor(2^512/p), C_1 = floor(2^512*beta/p) (limbs C_{j,m}):
  g_j = (sum_{i+m=7} lo(a_i C_{j,m}) + sum_{i+m=6} hi(a_i C_{j,m})) mod 2^64.
Lemma 1: if g_j <= 2^64 - 14 then floor(2^64 * M(beta^j x) / p) - g_j is in [0, 12].
Samples are uniform a in [1,p) plus edge cases (small a, a near p, a with a short top).
It also checks the kernel's decision rule for k against the exact minimum.
"""
import json
import random
import sys

p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
BETA = 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe
MASK = (1 << 64) - 1
C = [(1 << 512) // p, (BETA << 512) // p]
ERR = 13


def limbs(v, n):
    return [(v >> (64 * i)) & MASK for i in range(n)]


CL = [limbs(c, 8) for c in C]


def g(a, j):
    A = limbs(a, 6)
    s = 0
    for i in range(6):
        if 0 <= 7 - i < 8:
            s += (A[i] * CL[j][7 - i]) & MASK
        if 0 <= 6 - i < 8:
            s += (A[i] * CL[j][6 - i]) >> 64
    return s & MASK


def decide(a):
    """The kernel's rule: k, or None for the exact fallback."""
    g0, g1 = g(a, 0), g(a, 1)
    g2 = (-g0 - g1 - 2 * ERR) & MASK
    gs = [g0, g1, g2]
    if max(gs) > MASK - 2 * ERR - 1:
        return None
    k = gs.index(min(gs))
    if all(gs[l] - gs[k] > 2 * ERR for l in range(3) if l != k):
        return k
    return None


def main():
    if len(sys.argv) > 1:
        cfg = {"seed": int(sys.argv[1])}
    else:
        cfg = json.load(sys.stdin)
    rng = random.Random(cfg.get("seed", 1))
    samples = min(int(cfg.get("samples", 200000)), 400000)
    worst = 0
    decided = fallback = 0
    for t in range(samples):
        if t % 4:
            a = rng.randrange(1, p)
        else:
            a = rng.choice([rng.randrange(1, 1 << rng.randrange(1, 381)), p - rng.randrange(1, 1 << 70)])
        for j in (0, 1):
            gv = g(a, j)
            if gv > MASK - 13:
                continue
            d = ((a * pow(BETA, j, p) % p) << 64) // p - gv
            assert 0 <= d <= 12, (hex(a), j, d)
            worst = max(worst, d)
        k = decide(a)
        ms = [a, BETA * a % p, BETA * BETA * a % p]
        if k is None:
            fallback += 1
        else:
            assert k == ms.index(min(ms)), hex(a)
            decided += 1
    print(json.dumps({"samples": samples, "lemma1_holds": True, "worst_floor_gap": worst,
                      "decisions_checked_exact": decided, "fallbacks": fallback}))


if __name__ == "__main__":
    main()
