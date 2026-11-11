"""Lattice points per volume of the half-ellipsoids by an acceptance experiment (uniform integer points of the frame box of
common12.uniform_half_ellipsoid, exact hit test).  Trial i uses the six words splitmix64(seed, 6 i + j), y_j = word mod (2 cap_j + 1) - cap_j."""
from fractions import Fraction as F
from math import isqrt, pi, sqrt
import common12 as C
M64 = (1 << 64) - 1
def splitmix(seed, n):
    z = (seed + (n + 1) * 0x9E3779B97F4A7C15) & M64; z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64; z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64; return z ^ (z >> 31)
def frame(q, iq, radius2):
    V, fi = C._frame(q); cap = []
    for j in range(6): c = F(radius2) * fi[j][j]; cap.append(isqrt(c.numerator // c.denominator))
    return V, cap
def replay(q, iq, radius2, seed, i0, n, cuts):
    """Hits of trials i0..i0+n-1 below each cut (Q <= cut, x5 >= 1), pure Python."""
    V, cap = frame(q, iq, radius2); out = [0] * len(cuts)
    for i in range(i0, i0 + n):
        y = [splitmix(seed, 6 * i + j) % (2 * c + 1) - c for j, c in enumerate(cap)]; x = [sum(y[k] * V[k][j] for k in range(6)) for j in range(6)]
        if x[5] >= 1:
            v = C.quadratic(q, x)
            for k, cut in enumerate(cuts):
                if v <= cut: out[k] += 1
    return out
def ratio(q, den, radius2, cap, hits, n):
    """Points per volume of {x5 >= 1, Q <= radius2/den} from hits of n trials, and its relative SE."""
    import ideal_canonical as IC
    det = F(IC.determinant(q), den ** 6); box = 1
    for c in cap: box *= 2 * c + 1
    vol = pi ** 3 / 6 * (radius2 / den) ** 3 / (2 * sqrt(float(det))); return hits / n * box / vol, 1 / sqrt(hits)
