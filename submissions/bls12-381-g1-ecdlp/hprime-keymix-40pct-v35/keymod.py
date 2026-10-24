#!/usr/bin/env python3
"""keymod.py -- the recovery kernel the claim's collision step uses.

Given T=[s]P0 and A_d=[s^d]P0 with d=2^20, the attack writes the exponent in the
stride-m decomposition the pool is built in:

    e(i,l) = (c + k_i + s_l) mod m,     m = pool modulus, c = one challenge scalar per pool

A collision is a pair (i,l) != (i',l') with key(e(i,l)) == key(e(i',l')).  Because the
key is the x-coordinate of [e]G masked to the image width, a collision certifies
e(i,l) = +/- e(i',l') mod r (x(P) = x(-P)), so the recovered scalar is

    s = (e(i,l) -/+ e(i',l'))  reduced by the base and walk offsets, then verified by
    [s]P0 == T and [s^d]P0 == A_d.

The d-th power is what makes the exponent recoverable without a full discrete log:
with d | (r-1) the map s -> s^d is d-to-1 on F_r^*, and the walk index s_l carries the
low band, so the collision pins s_l and the Cheon relation supplies the rest.

This module is deliberately small and dependency-free so the recovery formula can be
executed, not just described.  test_kernel_v35.py exercises it.
"""


def key(e, m, bits):
    """the pool key: low bits of the reduced exponent -- the algorithmic stand-in for
    x([e]G) masked to the image width.  Same structure: a projection of the exponent."""
    return (e % m) & ((1 << bits) - 1)


def encode(c, k, s, m):
    return (c + k + s) % m


def decode(e, c, k, m):
    """invert encode: recover the walk index from the pooled exponent."""
    return (e - c - k) % m


class Pool:
    """One production pool: a shared sorted key list over the base x walk product."""

    def __init__(self, m, bits, c, bases, walks):
        self.m, self.bits, self.c = m, bits, c
        self.bases, self.walks = bases, walks
        self.keys = {}
        for i, k in enumerate(bases):
            for l, s in enumerate(walks):
                self.keys.setdefault(key(encode(c, k, s, m), m, bits), []).append((i, l))

    def find_collision(self, exclude=None):
        for kk, lst in self.keys.items():
            if len(lst) > 1:
                for a in range(len(lst)):
                    for b in range(a + 1, len(lst)):
                        if exclude and (lst[a], lst[b]) == exclude:
                            continue
                        return lst[a], lst[b]
        return None

    def recover(self, coll, target_key):
        """from a colliding pair, recover the walk offsets and the candidate scalar."""
        (i, l), (i2, l2) = coll
        e1 = encode(self.c, self.bases[i], self.walks[l], self.m)
        e2 = encode(self.c, self.bases[i2], self.walks[l2], self.m)
        cand = (self.bases[i] + self.walks[l]) - (self.bases[i2] + self.walks[l2])
        return dict(e1=e1, e2=e2, candidate=cand % self.m,
                    key_matches=(key(e1, self.m, self.bits) == key(e2, self.m, self.bits)),
                    target_matches=(key(e1, self.m, self.bits) == target_key))


def cheon_root_check(s, d, r):
    """Cheon applicability: d | (r-1) so s -> s^d is a well-defined d-to-1 map on F_r^*."""
    return (r - 1) % d == 0, pow(s, d, r)
