"""Check Lemma 1 of proof.md for the 4 x 64-bit kernel (step.c, canon_k).

Input (JSON on stdin, or a seed as the only argument): {"seed": int, "samples": int}.

x is fully reduced, held as limbs a_0..a_3 (2^64 radix). With
  C1 = floor(2^384 * beta / p) and limbs C1_1..C1_5 (C1_0 is not used),
  g_1 = (sum_{i+j=5} lo(a_i C1_j) + sum_{i+j=4} hi(a_i C1_j)) mod 2^64,
  g_0 = a_3.
Lemma 1: if g_j <= 2^64 - 10, then floor(2^64 * (beta^j x mod p) / p) - g_j is in [0, 8]
(for g_0 in [0, 1]). The script also checks the kernel's decision rule for k (FRAC_ERR = 11)
against the exact minimum of x, beta x, beta^2 x mod p.
"""
import json
import random
import sys

p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
BETA = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
MASK = (1 << 64) - 1
ERR = 11
C1 = (BETA << 384) // p
K = [(C1 >> (64 * j)) & MASK for j in range(6)]
assert K[1:] == [0x7a4a36aebcbb3d53, 0xc1396c28ec7e6dee, 0x9cf0497512f58995,
                 0x6e64479eac3434e9, 0x7ae96a2b657c0710]   # FRAC_K in step.c


def limbs(x):
    return [(x >> (64 * i)) & MASK for i in range(4)]


def g1(x):
    a = limbs(x)
    s = sum((a[i] * K[5 - i]) & MASK for i in range(4) if 1 <= 5 - i <= 5)
    s += sum((a[i] * K[4 - i]) >> 64 for i in range(4) if 1 <= 4 - i <= 5)
    return s & MASK


def decide(x):
    g0, gb = x >> 192, g1(x)
    g2 = (-g0 - gb - 2 * ERR) & MASK
    gs = [g0, gb, g2]
    if max(gs) > MASK - 2 * ERR - 1:
        return None
    k = gs.index(min(gs))
    if all(gs[l] - gs[k] > 2 * ERR for l in range(3) if l != k):
        return k
    return None


def main():
    cfg = {"seed": int(sys.argv[1])} if len(sys.argv) > 1 else json.load(sys.stdin)
    rng = random.Random(cfg.get("seed", 1))
    samples = min(int(cfg.get("samples", 200000)), 400000)
    worst = [0, 0]
    decided = fallback = fallback_uniform = 0
    for t in range(samples):
        if t % 4:
            x = rng.randrange(1, p)
        else:
            x = rng.choice([rng.randrange(1, 1 << rng.randrange(1, 256)), p - rng.randrange(1, 1 << 70)])
        for j, gv in ((0, x >> 192), (1, g1(x))):
            if gv > MASK - 9:
                continue
            d = ((x * pow(BETA, j, p) % p) << 64) // p - gv
            assert 0 <= d <= (1 if j == 0 else 8), (hex(x), j, d)
            worst[j] = max(worst[j], d)
        k = decide(x)
        ms = [x, BETA * x % p, BETA * BETA * x % p]
        if k is None:
            fallback += 1
            fallback_uniform += t % 4 != 0
        else:
            assert k == ms.index(min(ms)), hex(x)
            decided += 1
    print(json.dumps({"samples": samples, "lemma1_holds": True, "worst_floor_gap_g0": worst[0],
                      "worst_floor_gap_g1": worst[1], "decisions_checked_exact": decided,
                      "fallbacks": fallback,
                      "fallbacks_uniform_x": fallback_uniform}))


if __name__ == "__main__":
    main()
