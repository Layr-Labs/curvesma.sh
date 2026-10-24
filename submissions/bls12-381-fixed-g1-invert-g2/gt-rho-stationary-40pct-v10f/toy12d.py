#!/usr/bin/env python3
"""v10d faithful order-12 tower toy: the PRODUCTION class walk, and what v10c got wrong.

v10c's audit compared `keys12(y)` with `keys12(sigma^j y)` and read 528/528 FAIL.  That
comparison is against the *un-rotated* tuple, which is not the shipped class id.  The
shipped id is `kappa` = the twelve keys rotated so the minimum is first, and that id is
Frobenius invariant -- audit A below reads 0 failures.  So the v10c root cause ("the
class id is not Frobenius invariant") is WRONG, and the spec (`classwalk.twist_index`,
i = (tw - k) mod 12) is correct as written.

The defect that actually produced recovery 0.000 is the accelerated replay: v10c cached
the successor *element* per class and replayed it.  On a cache hit the held point is a
different member of the class than the one that materialised the record, so the stored
y2 is sigma^c(y * g^delta) while the labels are updated by delta alone.  Thereafter the
labels describe sigma^-c(y), every reported triple is off by a Frobenius power, and every
merge equation is unsatisfiable -- 0/289 solves, exactly as measured, with the class map
itself intact.

v10d keeps the class *decision* cached (legitimate: the decision is a class function,
audit B) and recomputes the field operation every step, so the invariant g^a z^b == y
holds at every step (audit C) and a merge solve verifies.
"""
import hashlib
import json
import random
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tower as TW

FIX = {
    5:  dict(r=9421,     lam=8943,   p=1152921504607943563, N=785,   spacing=8),
    7:  dict(r=37633,    lam=26672,  p=1152921504613142059, N=3136,  spacing=16),
    9:  dict(r=150589,   lam=100177, p=1152921504615061627, N=12549, spacing=32),
    11: dict(r=602137,   lam=117029, p=1152921504660984043, N=50178, spacing=64),
}
FRAC_SQR = 3
MU_WALK = 1.1
TARGET = 0.40
TABLE_ROWS = 64


def step_budget(N, mu=MU_WALK, target=TARGET, f=FRAC_SQR / 8.0):
    from math import log, isqrt
    s2 = 2 * mu * N * log(1 / (1 - target)) / (1 - f * f)
    s = isqrt(int(s2))
    while s * s < s2:
        s += 1
    while s > 1 and (s - 1) * (s - 1) >= s2:
        s -= 1
    return s


def fnv(vals):
    h = 0x811C9DC5
    for v in vals:
        x = v & 0xFFFFFFFFFFFFFFFF
        while True:
            h = ((h ^ (x & 0xFF)) * 0x01000193) & 0xFFFFFFFF
            x >>= 8
            if not x:
                break
    return h


class ClassWalk:
    """The production step rule.  Only the *decision* is cached, never an element."""

    def __init__(self, bits, seed):
        f = FIX[bits]
        self.bits, self.r, self.p, self.N, self.spacing = bits, f['r'], f['p'], f['N'], f['spacing']
        self.T = TW.Tower(self.p)
        self.rng = random.Random(seed)
        e = (self.p ** 12 - 1) // self.r
        while True:
            y = self.T.pow(self.rand_gt(), e)
            if y != self.T.one12() and self.T.pow(y, self.r) == self.T.one12():
                self.g = y
                break
        self.u = [self.rng.randrange(1, self.r) for _ in range(TABLE_ROWS)]
        self.ppow = [pow(self.p, i, self.r) for i in range(12)]
        self.tbl = []
        for t in range(TABLE_ROWS):
            y = self.T.pow(self.g, self.u[t])
            row = [y]
            for i in range(1, 12):
                row.append(self.T.frob(row[-1]))
            self.tbl.append(row)
        self.cache = {}          # class id -> decision (a class function of the id alone)
        self.materialised = 0
        self.deg_kappa = 0       # orbits where the minimum is not attained once

    def rand_gt(self):
        p, rng = self.p, self.rng
        return tuple(tuple((rng.randrange(p), rng.randrange(p)) for _ in range(3))
                     for _ in range(2))

    # ---- class key: the twelve keys rotated so the minimum is first (spec `kappa`)
    def kappa(self, y):
        ks = self.T.keys12(y)
        m = min(ks)
        k = ks.index(m)
        if ks.count(m) != 1:
            self.deg_kappa += 1
        return k, tuple(ks[k:] + ks[:k]), m

    def decision(self, cid):
        rec = self.cache.get(cid)
        if rec is not None:
            return rec
        h = fnv(list(cid))
        w = h % 8
        if w < FRAC_SQR:
            rec = ('sq', None, None, None)
        else:
            rec = ('mul', (h >> 8) % TABLE_ROWS, (h >> 20) % 12, None)
        self.cache[cid] = rec
        self.materialised += 1
        return rec

    def step(self, y, cid, k):
        """One real field operation: returns (y', delta or None, kind).  The twist is
        (tw - k) mod 12 with k the held element's own offset -- classwalk.twist_index."""
        kind, t, tw, _ = self.decision(cid)
        if kind == 'sq':
            return self.T.cyclo_sqr_fast(y), None, 'sq'
        i = (tw - k) % 12
        return self.T.mul(y, self.tbl[t][i]), (self.u[t] * self.ppow[i]) % self.r, 'mul'


# ---------------------------------------------------------------- audits
def audit(bits, seed):
    """A: the shipped class id is Frobenius invariant (48 x 11 = 528 checks).
       B: the step is a class function -- successor class identical on the orbit (528).
       Both counted, neither filtered."""
    f = FIX[bits]
    cw = ClassWalk(bits, seed)
    T = cw.T
    A_ok = A_bad = B_ok = B_bad = 0
    for el in range(48):
        x = cw.rng.randrange(1, cw.r)
        z = T.pow(cw.g, x)
        a, b = cw.rng.randrange(cw.r), cw.rng.randrange(cw.r)
        y0 = T.mul(T.pow(cw.g, a), T.pow(z, b))
        k0, cid0, _ = cw.kappa(y0)
        y2, d0, kind0 = cw.step(y0, cid0, k0)
        _, succ0, _ = cw.kappa(y2)
        yj = y0
        for j in range(1, 12):
            yj = T.frob(yj)
            kj, cidj, _ = cw.kappa(yj)
            if cidj == cid0:
                A_ok += 1
            else:
                A_bad += 1
            yj2, dj, kindj = cw.step(yj, cidj, kj)
            _, succj, _ = cw.kappa(yj2)
            if succj == succ0 and kindj == kind0:
                B_ok += 1
            else:
                B_bad += 1
    return dict(density_bits=bits, prime_r=cw.r,
                audit_A_key_invariance_checks=A_ok + A_bad, audit_A_pass=A_ok, audit_A_fail=A_bad,
                audit_B_step_class_function_checks=B_ok + B_bad, audit_B_pass=B_ok, audit_B_fail=B_bad,
                kappa_degenerate=cw.deg_kappa,
                scope='A: class id of a materialised element vs its 11 Frobenius conjugates; '
                      'B: successor class and step type from each conjugate vs the base')


def sweep(bits, trials, seed, inv_trials=4):
    """Protocol of record (v10b sweep4x.run_group): four trails of
    w = max(S//4, 8*spacing) steps sharing one DP store; a trail walks its whole
    length, every distinguished class is stored, and a repeat solves.  The trail is
    NOT restarted at its first report (v10c did that, and burned 3/4 of its budget
    cycling through non-DP classes)."""
    f = FIX[bits]
    r, spacing, N = f['r'], f['spacing'], f['N']
    S = step_budget(N)
    w_trail = max(S // 4, 8 * spacing)
    n_trails = 4
    budget = n_trails * w_trail
    cw = ClassWalk(bits, seed)
    T = cw.T
    rng = random.Random(seed ^ 0x5EED)
    st = dict(trials=0, success=0, steps=0, reports=0, trail_starts=0, merges=0,
              solve_success=0, degenerate_den=0, stranded_steps=0, false_confirm=0,
              inv_checks=0, inv_failures=0, repeated_class_revisits=0)
    for ti in range(trials):
        x = rng.randrange(1, r)
        z = T.pow(cw.g, x)
        store = {}
        ok = False
        for seg in range(n_trails):
            if ok:
                break
            st['trail_starts'] += 1
            a, b = rng.randrange(1, r), rng.randrange(1, r)
            y = T.mul(T.pow(cw.g, a), T.pow(z, b))
            k, cid, m = cw.kappa(y)
            for s_i in range(w_trail):
                y2, delta, kind = cw.step(y, cid, k)
                if kind == 'mul':
                    a = (a + delta) % r
                else:
                    a, b = (2 * a) % r, (2 * b) % r
                y = y2
                k, cid, m = cw.kappa(y)
                st['steps'] += 1
                if ti < inv_trials and seg == 0:          # exact label invariant, counted
                    st['inv_checks'] += 1
                    if T.mul(T.pow(cw.g, a), T.pow(z, b)) != y:
                        st['inv_failures'] += 1
                if (m ^ (m >> 32)) % spacing == 0:
                    st['reports'] += 1
                    prev = store.get(cid)
                    if prev is None:
                        store[cid] = (a, b, k)
                    elif (prev[0], prev[1]) != (a, b):
                        st['repeated_class_revisits'] += 1
                        a1, b1, k1 = prev
                        a2, b2, k2 = a, b, k
                        for j in ((k2 - k1) % 12, (k1 - k2) % 12):
                            pj = cw.ppow[j]
                            den = (b2 * pj - b1) % r
                            if den == 0:
                                st['degenerate_den'] += 1
                                continue
                            xg = (a1 - a2 * pj) % r * pow(den, -1, r) % r
                            st['merges'] += 1
                            if T.pow(cw.g, xg) == z:
                                st['solve_success'] += 1
                                st['success'] += 1
                                ok = True
                                break
                if ok:
                    break
        if not ok:
            st['stranded_steps'] += budget
        st['trials'] += 1
    out = dict(density_bits=bits, prime_r=r, field_p=f['p'], lam=f['lam'], N=N,
               spacing=spacing, S_ref=S, w_trail=w_trail, n_trails=n_trails,
               budget=budget, table_rows=TABLE_ROWS, frac_sqr='3/8', mu=MU_WALK,
               target=TARGET, classes_materialised=cw.materialised,
               kappa_degenerate=cw.deg_kappa, **st)
    out['recovery'] = '%d/%d' % (st['success'], st['trials'])
    out['recovery_frac'] = st['success'] / st['trials']
    return out


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'audit':
        print(json.dumps(audit(int(sys.argv[2]), 900 + int(sys.argv[2])), indent=1))
    else:
        bits = int(sys.argv[2])
        trials = int(sys.argv[3]) if len(sys.argv) > 3 else 1600
        print(json.dumps(sweep(bits, trials, 900 + bits), indent=1))
