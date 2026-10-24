#!/usr/bin/env python3
"""v14 faithful order-12 tower toy.  Everything executes THE SHIPPED STEP
(`classkey.step` with the depth-D ring / `classkey.classkey` /
`classkey.solve_candidates` / `routing.dp_hit`); the control rules in `classkey`
are the only other code paths and every property test runs on them too, so each
test is seen to FAIL on a mechanism that has the defect.

T1  RETURN-LAG LADDER, exhaustive over all r-1 elements, from an empty ring:
    for guard depth D in {0 (the v13 rule), 1, 2, 12} count elements whose image
    returns to the start at lag L.  The prediction is min lag = 2D+2 exactly:
    D=0 -> 2-cycles (the family the v13 verdict exhibits), D=1 -> nothing below 4,
    D=2 -> nothing below 6, D=12 -> nothing below 26.  This is the census in the
    space the trap lives in: it counts RETURNS OF THE WHOLE STATE, and the ring is
    part of the state.
T1j EXACT JOINT-SPACE CYCLE CENSUS at toy depth 1: the joint state is
    (element, ring) with ring in {empty} U (rows x 12), 97 x (r-1) states; the map
    is built for ALL of them and the full cycle-length histogram is reported.
    Prediction: 0 two-cycles, 0 three-cycles, four-cycles present.  The v13 rule
    (depth 0) on the same table is the control and must show its two-cycles.
T2  equivariance of the FULL step on SUCCESSOR ELEMENTS with matched rings:
    step(sigma^j y) must equal sigma^j step(y) exactly, for every j in Z/12.
    Controls: v12c class-keyed and v12b canonicalising (both split).
T2b the +-pair / label-bit algebra, unchanged from v13 (the guard does not touch
    the involution branch of the solve).
T3  coalescence: two sigma-partners of one class walked independently with fresh
    rings, keys compared at every checkpoint.  Plus the measured memory-divergence
    row (independent rings) that the model restatement in argument.md section 5
    refers to.
T4  census: 6400 trials per density, 4 trails, trails end at their first
    distinguished point, one shared store per trial, UNCONDITIONAL statistic,
    every recovery verified against the planted secret.
T5  end-to-end solve from (a,b) records across all twelve sigma-shifts.

Scope is stated, not glossed: the toy step budget is supercritical
(S/sqrt(N) ~ 2.3 at the fixture against 0.28 at production) and its DP depth is
2^-5..2^-9 against production 2^-80, so the census measures the MERGE -> SOLVE ->
VERIFY mechanism and its density trend, not a production recovery rate.
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tower as TW
import classkey as CK
import routing
import constants as K
from cp95 import cp95_lower

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


def enumerate_states(toy):
    out = [None] * toy.r
    y = toy.g
    for i in range(1, toy.r):
        out[i] = y
        y = toy.T.mul(y, toy.g)
    return out


def cycles_of(nxt, n):
    state = [0] * n
    hist = {}
    cyc_states = 0
    for s in range(n):
        if state[s]:
            continue
        path, pos = [], {}
        v = s
        while v is not None and not state[v] and v not in pos:
            pos[v] = len(path)
            path.append(v)
            v = nxt[v]
        if v is not None and v in pos:
            L = len(path) - pos[v]
            hist[L] = hist.get(L, 0) + 1
            cyc_states += L
        for u in path:
            state[u] = 1
    return dict(zip(map(str, sorted(hist)), [hist[k] for k in sorted(hist)])), cyc_states


def t1_return_ladder(toy64, toy8, idx, depth, maxlag):
    """Exhaustive return-lag census from an empty ring at guard depth `depth`.

    A return of the whole state at lag L is a cycle of length L, because the ring
    is part of the state and the start ring is empty: if y_L == y_0 the walk has
    come back to the exact state it left.
    """
    T, r = toy64.T, toy64.r
    out = dict(depth=depth, maxlag=maxlag, elements=r-1, lag_counts={})
    counts = {L: 0 for L in range(2, maxlag + 1)}
    for i in range(1, r):
        g = CK.Guard(depth) if depth else None
        y, a, b = idx[i], 1, 0
        ys = [y]
        for _s in range(maxlag):
            if depth:
                y, a, b, _k, _i = toy64.step(y, a, b, g)
            else:
                y, a, b, _k, _i = toy64.step_v13(y, a, b)
            ys.append(y)
        for L in range(2, maxlag + 1):
            if ys[L] == ys[0]:
                counts[L] += 1
    out['lag_counts'] = {str(L): counts[L] for L in counts}
    nz = [L for L in counts if counts[L]]
    out['min_lag'] = min(nz) if nz else None
    out['total_returns_below_26'] = sum(v for L, v in counts.items() if L <= 25)
    return out


def t1j_joint_census(toy8, idx, depth=1):
    """EXACT joint-space cycle census: every (element, ring) state, all of them."""
    T, r, rows = toy8.T, toy8.r, len(toy8.us)
    ring_states = [None]
    for t in range(rows):
        for c in range(12):
            ring_states.append((t, c))
    nr = len(ring_states)
    n = (r - 1) * nr
    ridx = {rs: j for j, rs in enumerate(ring_states)}
    nxt = [0] * n
    pos = {idx[i]: i for i in range(1, r)}
    for i in range(1, r):
        y0 = idx[i]
        base = (i - 1) * nr
        for j, rs in enumerate(ring_states):
            g = CK.Guard(depth) if depth else None
            if depth and rs is not None:
                g.slots[0] = rs
                for q in range(1, depth):
                    g.slots[q] = (-1, 0)
            y, a, b, kind, info = (toy8.step(y0, 1, 0, g) if depth
                                   else toy8.step_v13(y0, 1, 0))
            newring = None
            if depth:
                newring = g.slots[0] if g.slots[0][0] >= 0 else None
            if y in pos:
                nxt[base + j] = (pos[y] - 1) * nr + ridx[newring]
            else:
                nxt[base + j] = n          # the identity (and only it) is a sink
    nxt.append(n)
    n = n + 1
    hist, mass = cycles_of(nxt, n)
    return dict(depth=depth, states=n, ring_states=nr, cycles=hist, on_cycle=mass,
                two_cycles=hist.get('2', 0), three_cycles=hist.get('3', 0),
                four_cycles=hist.get('4', 0),
                min_cycle=(min(int(k) for k in hist) if hist else None))


def guard_fire_rate(rows, steps=4000, seed=5):
    """Measured guard fire rate per multiply step at table size `rows`.

    Production has rows = 2^14, so the rate there is ~D/R = 12/16384; this row is
    a fixture-scale measurement and the production retention is the analytic
    1 - 2*D*D/R, not a transfer of this number.
    """
    toy = Toy(rows=rows, seed=seed)
    rng = random.Random(seed)
    fires = muls = 0
    for w in range(steps // 40):
        y = toy.rand_gt()
        g = CK.Guard(K.GUARD_DEPTH)
        for _s in range(40):
            y, a, b, kind, info = toy.step(y, 1, 0, g)
            if kind == 'mul':
                muls += 1
                fires += info[3]
    return dict(rows=rows, mul_steps=muls, fires=fires,
                rate_per_mul_step=(fires / muls if muls else None),
                d_over_r=K.GUARD_DEPTH / rows)


def t1_control_cycles(toy8, idx, rows_label):
    """The v13 rule (depth 0) as an exact map on elements, small table."""
    T, r = toy8.T, toy8.r
    pos = {idx[i]: i for i in range(1, r)}
    nxt = [0] * (r + 1)
    for i in range(1, r):
        y, a, b, kind, info = toy8.step_v13(idx[i], 1, 0)
        nxt[i] = pos.get(y, r)
    nxt[r] = r
    hist, mass = cycles_of(nxt, r + 1)
    kinds = {'sq': 0, 'mul': 0}
    for i in range(1, r):
        _y, _a, _b, kind, _i = toy8.step_v13(idx[i], 1, 0)
        kinds[kind] += 1
    return dict(table_rows=rows_label, states=r-1, cycles=hist, on_cycle=mass,
                two_cycles=hist.get('2', 0), three_cycles=hist.get('3', 0),
                sqr_fraction=kinds['sq'] / (r - 1),
                predicted_2cycles_expectation=(r - 1) / (12.0 * rows_label))


def t2_equivariance(toy, idx, samples=64, steps=1):
    """step(sigma^j y) == sigma^j step(y) on ELEMENTS with matched rings."""
    T, r = toy.T, toy.r
    out = dict(samples=samples, shipped_split=0, shipped_key_split=0, checks=0,
               control_split_classkeyed=0, control_split_v12b=0,
               shipped_steps_sq=0, shipped_steps_mul=0, guard_fired=0)
    for s in range(samples):
        y = idx[1 + (s * 97) % (r - 1)]
        for j in range(12):
            ga, gb = CK.Guard(K.GUARD_DEPTH), CK.Guard(K.GUARD_DEPTH)
            ya, yb = y, toy.frob_pow(y, j)
            for _st in range(steps):
                ya, _a1, _b1, k1, i1 = toy.step(ya, 1, 0, ga)
                yb, _a2, _b2, k2, i2 = toy.step(yb, 1, 0, gb)
                out['guard_fired'] += i1[3] + i2[3]
                if k1 == 'sq':
                    out['shipped_steps_sq'] += 1
                else:
                    out['shipped_steps_mul'] += 1
            out['checks'] += 1
            if yb != toy.frob_pow(ya, j):
                out['shipped_split'] += 1
            if toy.key(yb)[0] != toy.key(ya)[0]:
                out['shipped_key_split'] += 1
            c1, _x, _y, _z, _w = toy.step_classkeyed(toy.frob_pow(y, j), 1, 0)
            cb, _x2, _y2, _z2, _w2 = toy.step_classkeyed(y, 1, 0)
            if toy.key(c1)[0] != toy.key(cb)[0]:
                out['control_split_classkeyed'] += 1
            v1, _x3, _y3 = toy.step_v12b(toy.frob_pow(y, j), 1, 0)
            vb, _x4, _y4 = toy.step_v12b(y, 1, 0)
            if v1 != toy.frob_pow(vb, j):
                out['control_split_v12b'] += 1
    out['shipped_clean'] = (out['shipped_split'] == 0 and out['shipped_key_split'] == 0)
    out['control_fires'] = (out['control_split_classkeyed'] > 0
                            and out['control_split_v12b'] > 0)
    return out


def t2b_pair_algebra(toy, samples=32, steps=8):
    """The +-pair / label-bit case, on a planted target (unchanged from v13)."""
    T, r, rng = toy.T, toy.r, random.Random(11)
    out = dict(samples=0, steps=steps, inverse_checks=0, inverse_holds=0,
               label_bit_attempts=0, label_bit_verified=0,
               control_attempts=0, control_wrong_sign_verified=0,
               zero_denominator=0, class_splits=0)
    for _ in range(samples):
        x = rng.randrange(1, r)
        z = T.pow(toy.g, x)
        a1, b1, b2 = rng.randrange(1, r), rng.randrange(1, r), rng.randrange(1, r)
        y1 = T.mul(T.pow(toy.g, a1), T.pow(z, b1))
        a2 = (-a1 - (b1 + b2) * x) % r
        y2 = T.mul(T.pow(toy.g, a2), T.pow(z, b2))
        if T.mul(y1, y2) != T.one12():
            continue
        out['samples'] += 1
        g1, g2 = CK.Guard(K.GUARD_DEPTH), CK.Guard(K.GUARD_DEPTH)
        for _s in range(steps):
            if toy.key(y1)[0] != toy.key(y2)[0]:
                out['class_splits'] += 1
            y1, a1, b1, _k, _i = toy.step(y1, a1, b1, g1)
            y2, a2, b2, _k2, _i2 = toy.step(y2, a2, b2, g2)
            out['inverse_checks'] += 1
            if T.mul(y1, y2) == T.one12():
                out['inverse_holds'] += 1
        xs = CK.label_bit_solve(a1, b1, a2, b2, r)
        out['label_bit_attempts'] += 1
        if xs is None:
            out['zero_denominator'] += 1
        elif T.pow(toy.g, xs) == z:
            out['label_bit_verified'] += 1
        den = (b1 - b2) % r
        if den:
            xc = (a2 - a1) % r * pow(den, -1, r) % r
            out['control_attempts'] += 1
            if T.pow(toy.g, xc) == z:
                out['control_wrong_sign_verified'] += 1
    out['inverse_clean'] = (out['inverse_holds'] == out['inverse_checks'])
    out['class_clean'] = (out['class_splits'] == 0)
    out['label_bit_ok'] = (out['label_bit_verified'] == out['label_bit_attempts']
                           and out['label_bit_attempts'] > 0)
    out['control_fires'] = (out['control_wrong_sign_verified'] == 0
                            and out['control_attempts'] > 0)
    return out


def t3_coalescence(toy, idx, pairs=64, steps=6):
    """Two sigma-partners of one class, independent labels, fresh rings."""
    T, r = toy.T, toy.r
    byclass = {}
    for i in range(1, r):
        byclass.setdefault(toy.key(idx[i])[0], []).append(i)
    keys = [k for k, v in byclass.items() if len(v) > 1]
    out = dict(pairs_available=len(keys), pairs=0, steps=steps, shipped_splits=0,
               control_splits=0, shipped_compared=0, control_compared=0,
               checkpoints=[], guard_fires=0)
    rng = random.Random(7)
    for _ in range(min(pairs, len(keys))):
        k = rng.choice(keys)
        i, j = rng.sample(byclass[k], 2)
        ya, yb = idx[i], idx[j]
        aa, ab, ba, bb = 1, 1, 0, 0
        ga, gb = CK.Guard(K.GUARD_DEPTH), CK.Guard(K.GUARD_DEPTH)
        for s in range(steps):
            if toy.key(ya)[0] != toy.key(yb)[0]:
                out['shipped_splits'] += 1
            out['shipped_compared'] += 1
            ya, aa, ba, _k, i1 = toy.step(ya, aa, ba, ga)
            yb, ab, bb, _k2, i2 = toy.step(yb, ab, bb, gb)
            out['guard_fires'] += i1[3] + i2[3]
        out['checkpoints'].append(1 if toy.key(ya)[0] == toy.key(yb)[0] else 0)
        ya, yb = idx[i], idx[j]
        for _s in range(steps):
            ya, _x, _y, _z, _w = toy.step_classkeyed(ya, 1, 0)
            yb, _x2, _y2, _z2, _w2 = toy.step_classkeyed(yb, 1, 0)
            out['control_compared'] += 1
            if toy.key(ya)[0] != toy.key(yb)[0]:
                out['control_splits'] += 1
    out['pairs'] = min(pairs, len(keys))
    out['shipped_clean'] = (out['shipped_splits'] == 0)
    out['control_fires'] = (out['control_splits'] > 0)
    out['coalesced_at_end'] = sum(out['checkpoints'])
    # memory-divergence row: the SAME class with deliberately different rings.
    # The applied column differs whenever exactly one of the two rings holds the
    # current row, and then the two successors are different classes.  This is
    # the priced consequence of a per-trail ring; argument.md section 5 restates
    # the birthday model over (class, ring) states and shows the class marginal is
    # unaffected.
    md = dict(pairs=0, steps=0, mul_steps=0, column_disagreements=0,
              merges_surviving_window=0, window=K.GUARD_DEPTH)
    rngm = random.Random(99)
    rows_t = len(toy.us)
    for _ in range(min(pairs, len(keys))):
        k = rngm.choice(keys)
        i, j = rngm.sample(byclass[k], 2)
        ya, yb = idx[i], idx[j]
        ga = CK.Guard(K.GUARD_DEPTH)
        gb = CK.Guard(K.GUARD_DEPTH)
        for q in range(K.GUARD_DEPTH):
            ga.slots[q] = (rngm.randrange(rows_t), rngm.randrange(12))
            gb.slots[q] = (rngm.randrange(rows_t), rngm.randrange(12))
        aa = ab = 1
        ba = bb = 0
        same = 0
        for _s in range(K.GUARD_DEPTH):
            ya, aa, ba, k1, i1 = toy.step(ya, aa, ba, ga)
            yb, ab, bb, k2, i2 = toy.step(yb, ab, bb, gb)
            md['steps'] += 1
            if k1 == 'mul':
                md['mul_steps'] += 1
                if i1[3] != i2[3]:
                    md['column_disagreements'] += 1
            if toy.key(ya)[0] != toy.key(yb)[0]:
                same = 1
        md['pairs'] += 1
        if not same:
            md['merges_surviving_window'] += 1
    out['memory_divergence'] = md
    return out


def t5_solver_end_to_end(toy, samples=12):
    """End-to-end from (a,b) records: for EVERY j in Z/12."""
    T, r, rng = toy.T, toy.r, random.Random(13)
    out = dict(samples=0, j_hits=0, j6_hits=0, wrong_j_accepted=0,
               zero_denominators=0, candidate_solves=0, verified=0, per_j={},
               unique_accept=0)
    for _ in range(samples):
        x = rng.randrange(1, r)
        z = T.pow(toy.g, x)
        a1, b1 = rng.randrange(1, r), rng.randrange(1, r)
        y1 = T.mul(T.pow(toy.g, a1), T.pow(z, b1))
        out['samples'] += 1
        for j in range(12):
            pj = pow(toy.p, j, r)
            b2 = rng.randrange(1, r)
            a2 = (pj * a1 + (pj * b1 - b2) * x) % r
            y2 = T.mul(T.pow(toy.g, a2), T.pow(z, b2))
            assert y2 == toy.frob_pow(y1, j)
            assert toy.key(y1)[0] == toy.key(y2)[0]
            cands = CK.solve_candidates([pow(toy.p, q, r) for q in range(12)],
                                        None, a1, b1, a2, b2, r)
            out['candidate_solves'] += 1
            out['per_j'][str(j)] = out['per_j'].get(str(j), 0) + len(cands)
            if not cands:
                out['zero_denominators'] += 1
                continue
            qtrue = (-j) % 12
            acc = [q for q, xx in cands if T.pow(toy.g, xx) == z]
            out['verified'] += len(acc)
            out['unique_accept'] += (1 if len(acc) == 1 else 0)
            if qtrue in acc:
                out['j_hits'] += 1
                if j == 6:
                    out['j6_hits'] += 1
            out['wrong_j_accepted'] += len([q for q in acc if q != qtrue])
    out['all_j_solved'] = (out['j_hits'] == 12 * out['samples'])
    out['no_false_accept'] = (out['wrong_j_accepted'] == 0)
    return out


def _census_batch(args):
    dp_bits, trials, seed, bits = args
    toy = Toy(bits=bits, seed=seed)
    T, r, rng = toy.T, toy.r, random.Random(seed ^ 0x5EED)
    theta = 1 << dp_bits
    cap = TRAIL_METERS * theta
    rec = dict(dp_bits=dp_bits, trials=trials, theta=theta, cap=cap, recovered=0,
               collisions=0, solved=0, zero_denominator=0, invariant_fail=0,
               records=0, steps=0, candidate_solves=0, verified=0,
               steps_first_report=0, first_report_trials=0, invariant_checks=0,
               attribution_fail=0, planted_keys=0, collisions_at_j6=0,
               collisions_at_j0=0, solve_candidates_used=0, guard_fires=0,
               mul_steps=0)
    for t in range(trials):
        x = rng.randrange(1, r)
        z = T.pow(toy.g, x)
        store = {}
        hit = got = False
        for trail in range(TRAILS):
            a = rng.randrange(1, r)
            b = rng.randrange(1, r)
            y = T.mul(T.pow(toy.g, a), T.pow(z, b))
            gr = CK.Guard(K.GUARD_DEPTH)
            for s in range(cap):
                y, a, b, kind, info = toy.step(y, a, b, gr)
                rec['steps'] += 1
                if kind == 'mul':
                    rec['mul_steps'] += 1
                    rec['guard_fires'] += info[3]
                if toy.dp(y, dp_bits):
                    if T.mul(T.pow(toy.g, a), T.pow(z, b)) != y:
                        rec['invariant_fail'] += 1
                    rec['invariant_checks'] += 1
                    rec['steps_first_report'] += s + 1
                    rec['first_report_trials'] += 1
                    k = toy.key(y)[0]
                    rec['records'] += 1
                    rec['planted_keys'] += 1
                    if k != toy.key(y)[0]:
                        rec['attribution_fail'] += 1
                    if k in store:
                        rec['collisions'] += 1
                        a2, b2 = store[k]
                        cands = CK.solve_candidates(
                            [pow(toy.p, j, r) for j in range(12)], None,
                            a, b, a2, b2, r)
                        rec['candidate_solves'] += len(cands)
                        rec['solve_candidates_used'] += 1
                        if not cands:
                            rec['zero_denominator'] += 1
                        for _j, xx in cands:
                            if T.pow(toy.g, xx) == z:
                                rec['verified'] += 1
                                if _j == 6:
                                    rec['collisions_at_j6'] += 1
                                if _j == 0:
                                    rec['collisions_at_j0'] += 1
                                got = True
                                break
                        hit = True
                        break
                    store[k] = (a, b)
                    break
            if hit:
                break
        if hit:
            rec['solved'] += 1
        if got:
            rec['recovered'] += 1
    return rec


def _merge(parts):
    out = dict(parts[0])
    for p in parts[1:]:
        for k, v in p.items():
            if isinstance(v, int) and k not in ('dp_bits', 'theta', 'cap'):
                out[k] += v
    return out


def t4_census(dp_bits, trials=TRIALS, bits=BITS, workers=16):
    import multiprocessing as mp
    nw = max(1, min(workers, trials))
    per = [trials // nw] * nw
    for i in range(trials - sum(per)):
        per[i] += 1
    args = [(dp_bits, n, 1000 + 17 * i, bits) for i, n in enumerate(per) if n]
    with mp.Pool(len(args)) as pool:
        parts = pool.map(_census_batch, args)
    rec = _merge(parts)
    rec['trials'] = trials
    rec['success'] = rec['recovered'] / trials
    rec['cp95_lower'] = float(cp95_lower(rec['recovered'], trials))
    rec['solve_rate_given_collision'] = (rec['solved'] / rec['collisions']
                                        if rec['collisions'] else None)
    rec['record_co_fraction'] = (rec['collisions'] / rec['records']
                                 if rec['records'] else None)
    rec['mean_steps_per_trial'] = rec['steps'] / trials
    rec['mean_steps_to_first_report'] = (rec['steps_first_report']
                                         / max(1, rec['first_report_trials']))
    rec['guard_fire_rate_per_mul_step'] = (rec['guard_fires'] / rec['mul_steps']
                                           if rec['mul_steps'] else None)
    rec['workers'] = len(args)
    return rec


def quick():
    toy = Toy()
    toy8 = Toy(rows=SMALL_ROWS, seed=777)
    idx = enumerate_states(toy)
    idx8 = enumerate_states(toy8)
    out = dict(
        fixture=dict(bits=toy.bits, r=toy.r, p_bits=toy.p.bit_length(), N=toy.N,
                     table_rows=len(toy.us), table_cols=12, guard_depth=K.GUARD_DEPTH,
                     ord_r_p=[i for i in range(1, 13) if pow(toy.p, i, toy.r) == 1][0],
                     sigma6_is_inversion=toy.sigma6_is_inversion,
                     column_identity_failures=toy.column_identity_failures,
                     column_involution_failures=toy.column_involution_failures,
                     column_inverse_pairs_visible=toy.column_inverse_pairs,
                     small_table_rows=SMALL_ROWS),
        kappa_bits=K.KAPPA_BITS,
        routing=routing.selftest(n=20000),
        t1=dict(
            d12=t1_return_ladder(toy, toy8, idx, K.GUARD_DEPTH, 26),
            d2=t1_return_ladder(toy8, toy8, idx8, 2, 8),
            d1=t1_return_ladder(toy8, toy8, idx8, 1, 6),
            d0_v13_rule=t1_return_ladder(toy8, toy8, idx8, 0, 8),
            control_v13_rule_small_table=t1_control_cycles(toy8, idx8, SMALL_ROWS),
            control_v13_rule_v13_table=t1_control_cycles(toy, idx, TABLE_ROWS),
        ),
        t1j=dict(
            guarded_d1=t1j_joint_census(toy8, idx8, depth=1),
            v13_rule=t1j_joint_census(toy8, idx8, depth=0),
        ),
        guard_rate=dict(r64=guard_fire_rate(64), r512=guard_fire_rate(512)),
        t2=t2_equivariance(toy, idx),
        t2b=t2b_pair_algebra(toy),
        t3=t3_coalescence(toy, idx),
        t5=t5_solver_end_to_end(toy),
        t4=dict(),
        cp_selftest=dict(case_0_of_10=str(cp95_lower(0, 10)),
                         case_9_of_10=str(cp95_lower(9, 10)),
                         case_10_of_10=str(cp95_lower(10, 10)),
                         case_8_of_10=str(cp95_lower(8, 10))),
    )
    os.makedirs('evidence', exist_ok=True)
    with open('evidence/toy-census.json', 'w') as f:
        json.dump(out, f, sort_keys=True, default=str, separators=(',', ':'))
    t1 = out['t1']
    print("fixture r=%d N=%d classes=%d sigma6=inv %s colidfail %d involutionfail %d invpairs %d"
          % (toy.r, toy.N, t1['d12']['elements'], toy.sigma6_is_inversion,
             toy.column_identity_failures, toy.column_involution_failures,
             toy.column_inverse_pairs))
    for tag in ('d0_v13_rule', 'd1', 'd2', 'd12'):
        t = t1[tag]
        print("T1 %-12s min_lag %s  returns<=25 %d  lag_counts %s"
              % (tag, t['min_lag'], t['total_returns_below_26'],
                 {k: v for k, v in t['lag_counts'].items() if v}))
    c = t1['control_v13_rule_small_table']
    print("T1 control v13 rule R=8: cycles %s (2c %d, expected %.1f) sqr_frac %.4f"
          % (c['cycles'], c['two_cycles'], c['predicted_2cycles_expectation'],
             c['sqr_fraction']))
    c64 = t1['control_v13_rule_v13_table']
    print("T1 control v13 rule R=64: cycles %s (2c %d) <- the v13 fixture"
          % (c64['cycles'], c64['two_cycles']))
    for tag in ('guarded_d1', 'v13_rule'):
        j = out['t1j'][tag]
        print("T1j %-12s joint states %d  cycles %s  min_cycle %s  2c %d 3c %d 4c %d"
              % (tag, j['states'], j['cycles'], j['min_cycle'], j['two_cycles'],
                 j['three_cycles'], j['four_cycles']))
    print("T2", out['t2'])
    print("T2b", out['t2b'])
    print("T3", {k: v for k, v in out['t3'].items() if k not in ('checkpoints',)})
    print("T5", out['t5'])
    return out


def census_one(d):
    out = json.load(open('evidence/toy-census.json'))
    r4 = t4_census(d)
    out['t4'][str(d)] = r4
    with open('evidence/toy-census.json', 'w') as f:
        json.dump(out, f, sort_keys=True, default=str, separators=(',', ':'))
    with open('evidence/t4-dp%d.json' % d, 'w') as f:
        json.dump(r4, f, sort_keys=True, default=str, separators=(',', ':'))
    print("T4 dp2^-%d: success %.4f CP95 %s collisions %d verified %d inv_fail %d "
          "zeroden %d attr_fail %d rec_co %.4f mean_steps %.1f mean_to_report %.1f "
          "guard_rate %.5f"
          % (d, r4['success'], r4['cp95_lower'], r4['collisions'], r4['verified'],
             r4['invariant_fail'], r4['zero_denominator'], r4['attribution_fail'],
             r4['record_co_fraction'], r4['mean_steps_per_trial'],
             r4['mean_steps_to_first_report'], r4['guard_fire_rate_per_mul_step']))
    return r4


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'census':
        for dd in [int(v) for v in sys.argv[2:]]:
            census_one(dd)
    else:
        quick()
