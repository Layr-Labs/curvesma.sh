#!/usr/bin/env python3
"""v32 toy fixture -- the SAME 785-class order-12 tower the census uses.

Extracted verbatim from v27 toy_census.py (class Toy and FIX); the walk itself
lives in protocol.py and is NOT duplicated here.  The fixture builds the twelve
pre-expanded columns TBL[t][k] = sigma^-k(M_t) and their discrete logs U[t][k],
so step_v26 can be executed exactly as shipped.
"""
import random
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tower as TW
import classkey as CK
import routing
FIX = {
    5: dict(r=9421,   p=1152921504607943563, N=785),
    7: dict(r=37633,  p=1152921504613142059, N=3136),
}
BITS = 5
TABLE_ROWS = 64
SMALL_ROWS = 8            # the control table on which the 1/(12R) family is visible
DP_DENSITIES = (5, 7, 9)
TRIALS = 6400
TRAIL_METERS = 16
TRAILS = 4


class Toy:
    def __init__(self, bits=BITS, seed=20261024, rows=TABLE_ROWS):
        f = FIX[bits]
        self.bits, self.r, self.p, self.N = bits, f['r'], f['p'], f['N']
        self.T = TW.Tower(self.p)
        self.rng = random.Random(seed)
        e = (self.p ** 12 - 1) // self.r
        while True:
            y = self.T.pow(self.rand_raw(), e)
            if y != self.T.one12() and self.T.pow(y, self.r) == self.T.one12():
                self.g = y
                break
        self.stride = self.rng.randrange(1, self.r)
        self.c0 = self.rng.randrange(1, self.r)
        self.us = [((self.c0 + t * self.stride) % self.r) or 1 for t in range(rows)]
        assert all(1 <= u < self.r for u in self.us)
        self.tbl = [self.T.pow(self.g, u) for u in self.us]
        self.tbits = rows.bit_length() - 1
        self.rows = rows
        self.ppow = [pow(self.p, i, self.r) for i in range(12)]
        self._build_frob_pow()
        # the TWELVE PRE-EXPANDED COLUMNS: TBL[t][k] = sigma^{-k}(M_t), U[t][k]
        self.TBL = [[self.frob_pow(self.tbl[t], (-k) % 12) for k in range(12)]
                    for t in range(rows)]
        self.U = [[self.us[t] * self.ppow[(-k) % 12] % self.r for k in range(12)]
                  for t in range(rows)]
        chk = 0
        for t in (0, 1, rows - 1):
            for k in (0, 1, 5, 6, 11):
                if self.T.pow(self.g, self.U[t][k]) != self.TBL[t][k]:
                    chk += 1
                if self.frob_pow(self.tbl[t], (-k) % 12) != self.TBL[t][k]:
                    chk += 1
        self.sigma6_is_inversion = all(
            self.frob_pow(y, 6) == self.T.pow(y, self.r - 1) for y in
            [self.tbl[i] for i in range(0, rows, 3)])
        self.column_identity_failures = chk
        # NO COLUMN IS AN INVOLUTION: TBL[t][k]^2 != 1 for every k, which is what
        # makes "same column twice" a non-return (2*u*p^-k != 0 mod r).
        self.column_involution_failures = sum(
            1 for t in range(rows) for k in range(12)
            if self.T.mul(self.TBL[t][k], self.TBL[t][k]) == self.T.one12())
        # and TBL[t][k+6] IS the inverse of TBL[t][k]: the v13 family exists here
        self.column_inverse_pairs = sum(
            1 for t in range(rows) for k in range(6)
            if self.T.mul(self.TBL[t][k], self.TBL[t][k + 6]) == self.T.one12())

    def rand_raw(self):
        p, rng = self.p, self.rng
        return tuple(tuple((rng.randrange(p), rng.randrange(p)) for _ in range(3))
                     for _ in range(2))

    def rand_gt(self, rng=None):
        """A fresh GT element.  rng lets a WORKER draw from its OWN stream: with no
        rng every worker replays toy.rng, so worker i's j-th start equals worker
        0's j-th start and the 'independent repeats' are the same trajectories."""
        src = self.rng if rng is None else rng
        return self.T.pow(self.g, src.randrange(1, self.r))

    # ---- THE v24 shipped step: guard-free + the same-row family removed
    def step_live(self, y, a, b, prev_row=-1):
        return CK.step_v24(self.T, y, a, b, self.TBL, self.U, self.r,
                           self.tbits, self.rows, prev_row)

    # ---- the v23 shipped step, kept as the control that SHOWS the family
    def step_v23(self, y, a, b):
        return CK.step_v13(self.T, y, a, b, self.TBL, self.U, self.r, self.tbits)


    def _build_frob_pow(self):
        T = self.T
        self.S = [[None] * 3 for _ in range(12)]
        for k in range(12):
            for i in range(3):
                v = [(0, 0), (0, 0), (0, 0)]
                v[i] = (1, 0)
                v = (v[0], v[1], v[2])
                for _ in range(k):
                    v = T.frob6(v)
                self.S[k][i] = v[i]
        self.G = [T.one6[0]]
        for _ in range(1, 12):
            self.G.append(T.f2mul(T.f2conj(self.G[-1]), T.gamma))

    def frob_pow(self, y, k):
        k %= 12
        if k == 0:
            return y
        T, S, G = self.T, self.S[k], self.G[k]
        out = []
        for half in y:
            if k % 2:
                half = (T.f2conj(half[0]), T.f2conj(half[1]), T.f2conj(half[2]))
            out.append(tuple(T.f2mul(half[i], S[i]) for i in range(3)))
        return (out[0], T.f6scale(out[1], G))

    # ---- THE shipped step (guarded)
    def step(self, y, a, b, guard=None):
        return CK.step(self.T, y, a, b, self.TBL, self.U, self.r, self.tbits, guard)

    # ---- the v13 rule (depth 0): the control that must show the family
    def step_v13(self, y, a, b):
        return CK.step_v13(self.T, y, a, b, self.TBL, self.U, self.r, self.tbits)

    def key(self, y):
        return CK.classkey(self.T, y)

    def rot(self, y):
        return CK.classkey(self.T, y)[1]

    def dp(self, y, dp_bits):
        return routing.dp_hit(self.key(y)[0], dp_bits)

    # ---- CONTROLS
    def step_v12b(self, y, a, b):
        return CK.step_v10f(self.T, y, a, b, self.tbl, self.us, self.r,
                            self.frob_pow, self.ppow, {}, self.tbits)

    def step_classkeyed(self, y, a, b):
        return CK.step_classkeyed(self.T, y, a, b, self.tbl, self.us, self.r,
                                  self.tbits)

