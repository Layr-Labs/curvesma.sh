#!/usr/bin/env python3
"""v13 faithful order-12 tower toy.  Everything executes THE SHIPPED STEP
(`classkey.step` / `classkey.classkey` / `classkey.solve_candidates` /
`routing.dp_hit`); the control rules in `classkey` are the only other code paths
and every property test runs on them too, so each test is seen to FAIL on a
mechanism that has the defect.

T1  exhaustive cycle census over ALL r-1 states and ALL (r-1)/12 classes for the
    shipped rule, plus the v12b canonicalising control (its 2-cycles must appear)
    and the v12c class-keyed control; step-kind census = measured f_sqr.
T2  equivariance of the FULL step on SUCCESSOR ELEMENTS: step(sigma^j y) must
    equal sigma^j step(y) exactly, for every j in Z/12, with a counter incremented
    per mismatch.  Controls: v12c class-keyed and v12b canonicalising (both split).
T2b the +-pair / label-bit algebra: two walks started at y and at y^{-1} with
    INDEPENDENT label pairs stay exact inverses; the j=6 solve recovers the
    planted x and the j=0 control (sign bit flipped) FAILS verification.
T3  coalescence: two sigma-partners of one class walked independently, keys
    compared at every checkpoint.  Control: v12c class-keyed.
T4  census: 6400 trials per density, 4 trails, trails end at their first
    distinguished point, one shared store per trial, UNCONDITIONAL success
    statistic, every recovery verified against the planted secret by
    exponentiation, per-report invariant g^a z^b == y asserted, record
    co-fraction measured.

Scope is stated, not glossed: the toy step budget is supercritical
(S/sqrt(N) ~ 2.3 at the fixture against 0.28 at production) and its DP depth is
2^-5..2^-9 against production 2^-80, so the census measures the MERGE -> SOLVE ->
VERIFY mechanism and its density trend, not a production recovery rate.  The
production recovery model is the analytic birthday expression in
ledger.model_success.
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
        self.us = [(self.c0 + t * self.stride) % self.r for t in range(rows)]
        assert all(1 <= u < self.r for u in self.us)
        self.tbl = [self.T.pow(self.g, u) for u in self.us]
        self.tbits = rows.bit_length() - 1
        self.ppow = [pow(self.p, i, self.r) for i in range(12)]
        self._build_frob_pow()
        # the TWELVE PRE-EXPANDED COLUMNS: TBL[t][k] = sigma^{-k}(M_t), U[t][k]
        self.TBL = [[self.frob_pow(self.tbl[t], (-k) % 12) for k in range(12)]
                    for t in range(rows)]
        self.U = [[self.us[t] * self.ppow[(-k) % 12] % self.r for k in range(12)]
                  for t in range(rows)]
        # column identity check: TBL[t][k] must be g^{U[t][k]} and the sigma-algebra
        chk = 0
        for t in (0, 1, 7, rows - 1):
            for k in (0, 1, 5, 6, 11):
                g_u = self.T.pow(self.g, self.U[t][k])
                if g_u != self.TBL[t][k]:
                    chk += 1
                if self.frob_pow(self.tbl[t], (-k) % 12) != self.TBL[t][k]:
                    chk += 1
        # and sigma^6 IS inversion on this fixture
        self.sigma6_is_inversion = all(
            self.frob_pow(y, 6) == self.T.pow(y, self.r - 1) for y in
            [self.tbl[i] for i in range(0, rows, 7)])
        self.column_identity_failures = chk

    def rand_raw(self):
        p, rng = self.p, self.rng
        return tuple(tuple((rng.randrange(p), rng.randrange(p)) for _ in range(3))
                     for _ in range(2))

    def rand_gt(self):
        return self.T.pow(self.g, self.rng.randrange(1, self.r))

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

    # ---- THE shipped step
    def step(self, y, a, b):
        return CK.step(self.T, y, a, b, self.TBL, self.U, self.r, self.tbits)

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


def idx_inv(toy, idx, y):
    d = toy.__dict__.setdefault('_inv', {})
    if not d:
        for i in range(1, toy.r):
            d[idx[i]] = i
    return d.get(y)


def t1_cycles(toy, idx):
    T, r = toy.T, toy.r
    n = r - 1
    nxt_state = [None] * r
    kinds = {'sq': 0, 'mul': 0}
    for i in range(1, r):
        y2, a2, b2, kind, info = toy.step(idx[i], 1, 0)
        kinds[kind] += 1
        nxt_state[i] = idx_inv(toy, idx, y2)
    st_state, cyc_state = cycles_of(nxt_state, r)
    # class level
    cid, rep = {}, {}
    for i in range(1, r):
        k = toy.key(idx[i])[0]
        c = cid.setdefault(k, len(cid))
        rep.setdefault(c, i)
    ncls = len(cid)
    nxt_cls = [None] * ncls
    for c, i in rep.items():
        y2, _a, _b, _k, _i = toy.step(idx[i], 1, 0)
        nxt_cls[c] = cid[toy.key(y2)[0]]
    st_cls, cyc_cls = cycles_of(nxt_cls, ncls)
    # controls
    nxt_v12b = [None] * r
    for i in range(1, r):
        y2, _a, _b = toy.step_v12b(idx[i], 1, 0)
        nxt_v12b[i] = idx_inv(toy, idx, y2)
    st_v12b, cyc_v12b = cycles_of(nxt_v12b, r)
    nxt_v12c = [None] * r
    for i in range(1, r):
        y2, _a, _b, _k, _i = toy.step_classkeyed(idx[i], 1, 0)
        nxt_v12c[i] = idx_inv(toy, idx, y2)
    st_v12c, cyc_v12c = cycles_of(nxt_v12c, r)
    return dict(states=n, classes=ncls, classes_exact=(ncls == toy.N),
                shipped_state_cycles=st_state, shipped_state_on_cycles=cyc_state,
                shipped_state_2cycles=st_state.get('2', 0),
                shipped_class_cycles=st_cls, shipped_class_on_cycles=cyc_cls,
                shipped_class_2cycles=st_cls.get('2', 0),
                step_kinds=kinds, sqr_fraction=kinds['sq'] / n,
                control_v12b_state_cycles=st_v12b,
                control_v12b_state_2cycles=st_v12b.get('2', 0),
                control_v12b_on_cycles=cyc_v12b,
                control_v12c_state_cycles=st_v12c,
                control_v12c_state_2cycles=st_v12c.get('2', 0),
                control_v12c_on_cycles=cyc_v12c)


def t2_equivariance(toy, idx, samples=64):
    """step(sigma^j y) == sigma^j step(y) on ELEMENTS, every j, counter per check."""
    T, r = toy.T, toy.r
    out = dict(samples=samples, shipped_split=0, shipped_key_split=0, checks=0,
               control_split_classkeyed=0, control_split_v12b=0,
               shipped_steps_sq=0, shipped_steps_mul=0)
    for s in range(samples):
        y = idx[1 + (s * 97) % (r - 1)]
        base, _a, _b, _k, _i = toy.step(y, 1, 0)
        kb = toy.key(base)[0]
        for j in range(12):
            yj = toy.frob_pow(y, j)
            sj, _a2, _b2, kind2, _i2 = toy.step(yj, 1, 0)
            out['checks'] += 1
            if sj != toy.frob_pow(base, j):            # ELEMENT equality
                out['shipped_split'] += 1
            if toy.key(sj)[0] != kb:                   # key equality
                out['shipped_key_split'] += 1
            if kind2 == 'sq':
                out['shipped_steps_sq'] += 1
            else:
                out['shipped_steps_mul'] += 1
            c1, _x, _y, _z, _w = toy.step_classkeyed(yj, 1, 0)
            cb, _x2, _y2, _z2, _w2 = toy.step_classkeyed(y, 1, 0)
            if toy.key(c1)[0] != toy.key(cb)[0]:
                out['control_split_classkeyed'] += 1
            v1, _x3, _y3 = toy.step_v12b(yj, 1, 0)
            vb, _x4, _y4 = toy.step_v12b(y, 1, 0)
            if v1 != toy.frob_pow(vb, j):
                out['control_split_v12b'] += 1
    out['shipped_clean'] = (out['shipped_split'] == 0 and out['shipped_key_split'] == 0)
    out['control_fires'] = (out['control_split_classkeyed'] > 0
                            and out['control_split_v12b'] > 0)
    return out


def t2b_pair_algebra(toy, samples=32, steps=8):
    """The +-pair / label-bit case, on a planted target.

    Two walks start at y1 = g^a1 z^b1 and y2 = y1^{-1} = g^a2 z^b2 with labels that
    are NOT simple negations (b2 is drawn independently and a2 is placed so that
    the pair is exact).  They stay exact inverses (the +-pair is a fixed point of
    the equivariant map), and the j=6 label-bit solve recovers the planted x.
    Control: flip the sign bit (assume j=0) and verification FAILS.
    """
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
        for _s in range(steps):
            if toy.key(y1)[0] != toy.key(y2)[0]:
                out['class_splits'] += 1
            y1, a1, b1, _k, _i = toy.step(y1, a1, b1)
            y2, a2, b2, _k2, _i2 = toy.step(y2, a2, b2)
            out['inverse_checks'] += 1
            if T.mul(y1, y2) == T.one12():
                out['inverse_holds'] += 1
        # the label-bit solve (j = 6)
        xs = CK.label_bit_solve(a1, b1, a2, b2, r)
        out['label_bit_attempts'] += 1
        if xs is None:
            out['zero_denominator'] += 1
        elif T.pow(toy.g, xs) == z:
            out['label_bit_verified'] += 1
        # control: assume y1 == y2 (sign bit flipped, j = 0)
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
    T, r = toy.T, toy.r
    byclass = {}
    for i in range(1, r):
        byclass.setdefault(toy.key(idx[i])[0], []).append(i)
    keys = [k for k, v in byclass.items() if len(v) > 1]
    out = dict(pairs_available=len(keys), pairs=0, steps=steps, shipped_splits=0,
               control_splits=0, shipped_compared=0, control_compared=0,
               checkpoints=[])
    rng = random.Random(7)
    for _ in range(min(pairs, len(keys))):
        k = rng.choice(keys)
        i, j = rng.sample(byclass[k], 2)
        ya, yb = idx[i], idx[j]
        aa, ab, ba, bb = 1, 1, 0, 0
        for s in range(steps):
            if toy.key(ya)[0] != toy.key(yb)[0]:
                out['shipped_splits'] += 1
            out['shipped_compared'] += 1
            ya, aa, ba, _k, _i = toy.step(ya, aa, ba)
            yb, ab, bb, _k2, _i2 = toy.step(yb, ab, bb)
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
    return out


def t5_solver_end_to_end(toy, samples=12):
    """End-to-end from (a,b) records: for EVERY j in Z/12, two records that are
    sigma^j-partners of one class must solve to the planted x by the shipped
    twelve-way formula, and every wrong j must be rejected by verification.
    j = 6 is the +-pair case and is counted separately."""
    T, r, rng = toy.T, toy.r, random.Random(13)
    out = dict(samples=0, j_hits=0, j6_hits=0, wrong_j_accepted=0,
               zero_denominators=0, candidate_solves=0, verified=0, per_j={})
    for _ in range(samples):
        x = rng.randrange(1, r)
        z = T.pow(toy.g, x)
        a1, b1 = rng.randrange(1, r), rng.randrange(1, r)
        y1 = T.mul(T.pow(toy.g, a1), T.pow(z, b1))
        out['samples'] += 1
        for j in range(12):
            pj = pow(toy.p, j, r)
            # INDEPENDENT label pair placed exactly on the sigma^j-partner:
            # a2 = p^j a1 + (p^j b1 - b2) x, so y2 = sigma^j(y1) and b2 is free.
            b2 = rng.randrange(1, r)
            a2 = (pj * a1 + (pj * b1 - b2) * x) % r
            y2 = T.mul(T.pow(toy.g, a2), T.pow(z, b2))
            assert y2 == toy.frob_pow(y1, j)          # the sigma-relation is exact
            assert toy.key(y1)[0] == toy.key(y2)[0]   # same class, both records
            cands = CK.solve_candidates([pow(toy.p, q, r) for q in range(12)],
                                        None, a1, b1, a2, b2, r)
            out['candidate_solves'] += 1
            out['per_j'][str(j)] = out['per_j'].get(str(j), 0) + len(cands)
            if not cands:
                out['zero_denominators'] += 1
                continue
            # y2 = sigma^j(y1) means y1 = sigma^{-j}(y2): the accepted shift is
            # q = -j mod 12, which for j = 6 is j itself (the +-pair case).
            qtrue = (-j) % 12
            acc = [q for q, xx in cands if T.pow(toy.g, xx) == z]
            out['verified'] += len(acc)
            out['unique_accept'] = out.get('unique_accept', 0) + (1 if len(acc) == 1 else 0)
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
               collisions_at_j0=0, solve_candidates_used=0)
    for t in range(trials):
        x = rng.randrange(1, r)
        z = T.pow(toy.g, x)
        store = {}
        hit = got = False
        for trail in range(TRAILS):
            a = rng.randrange(1, r)
            b = rng.randrange(1, r)
            y = T.mul(T.pow(toy.g, a), T.pow(z, b))
            for s in range(cap):
                y, a, b, kind, _info = toy.step(y, a, b)
                rec['steps'] += 1
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
    rec['workers'] = len(args)
    return rec


def quick():
    toy = Toy()
    idx = enumerate_states(toy)
    out = dict(
        fixture=dict(bits=toy.bits, r=toy.r, p_bits=toy.p.bit_length(), N=toy.N,
                     table_rows=len(toy.us), table_cols=12,
                     ord_r_p=[i for i in range(1, 13) if pow(toy.p, i, toy.r) == 1][0],
                     sigma6_is_inversion=toy.sigma6_is_inversion,
                     column_identity_failures=toy.column_identity_failures),
        kappa_bits=K.KAPPA_BITS,
        routing=routing.selftest(n=20000),
        t1=t1_cycles(toy, idx),
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
        json.dump(out, f, indent=1, sort_keys=True, default=str)
    t1 = out['t1']
    print("fixture r=%d N=%d classes=%d exact=%s sigma6=inv %s colidfail %d"
          % (toy.r, toy.N, t1['classes'], t1['classes_exact'],
             toy.sigma6_is_inversion, toy.column_identity_failures))
    print("T1 shipped: state cycles %s (2c %d) class cycles %s (2c %d) kinds %s sqr_frac %.5f"
          % (t1['shipped_state_cycles'], t1['shipped_state_2cycles'],
             t1['shipped_class_cycles'], t1['shipped_class_2cycles'],
             t1['step_kinds'], t1['sqr_fraction']))
    print("T1 control v12b cycles %s (2c %d)  v12c cycles %s (2c %d)"
          % (t1['control_v12b_state_cycles'], t1['control_v12b_state_2cycles'],
             t1['control_v12c_state_cycles'], t1['control_v12c_state_2cycles']))
    print("T2", out['t2'])
    print("T2b", out['t2b'])
    print("T3", {k: v for k, v in out['t3'].items() if k != 'checkpoints'})
    print("T5", out['t5'])
    return out


def census_one(d):
    out = json.load(open('evidence/toy-census.json'))
    r4 = t4_census(d)
    out['t4'][str(d)] = r4
    with open('evidence/toy-census.json', 'w') as f:
        json.dump(out, f, indent=1, sort_keys=True, default=str)
    with open('evidence/t4-dp%d.json' % d, 'w') as f:
        json.dump(r4, f, indent=1, sort_keys=True, default=str)
    print("T4 dp2^-%d: success %.4f CP95 %s collisions %d verified %d inv_fail %d "
          "zeroden %d attr_fail %d rec_co %.4f mean_steps %.1f mean_to_report %.1f"
          % (d, r4['success'], r4['cp95_lower'], r4['collisions'], r4['verified'],
             r4['invariant_fail'], r4['zero_denominator'], r4['attribution_fail'],
             r4['record_co_fraction'], r4['mean_steps_per_trial'],
             r4['mean_steps_to_first_report']))
    return r4


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'census':
        for dd in [int(v) for v in sys.argv[2:]]:
            census_one(dd)
    else:
        quick()
