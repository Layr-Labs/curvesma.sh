#!/usr/bin/env python3
"""v12b faithful order-12 tower toy: the CANONICALISING walk, with a cycle census.

Three things this toy does that v10f's did not:

  1. it executes the shipped v12b step verbatim -- canonical representative held,
     label renormalisation (a,b) <- ((a+delta) p^k, b p^k) on every step, table
     built from independent lower-half exponents -- and imports the same rule
     functions from `classwalk`;
  2. it censuses ORBITS, not outcomes: every trail records exact state repeats
     (class id, a, b) and flags any orbit whose length is below the declared
     bound, so a stationary cycle cannot hide behind a recovery fraction;
  3. it solves with the v12b solve (x = (a1-a2)/(b2-b1) mod r, no k-search) and
     verifies every recovered x against the planted secret by exponentiation.

The two-to-three order-of-magnitude gap between toy and production cycle rates is
stated, not hidden: at r ~ 2^15 a lap that uses one row's complete index set
returns element and labels exactly, because p has order 12 and p != 1 gives
sum_{i<12} p^i = 0 (mod r).  That family is negligible at production (2^14 * 12
multipliers against r ~ 2^255) and dominant here, so toy cycle LENGTHS are
fixture-limited evidence; the production statement about the (m, m^-1) family is
the algebraic identity in argument.md, not a toy rate.
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tower as TW
import classwalk as CW
import constants as K

FIX = {
    5:  dict(r=9421,     lam=8943,   p=1152921504607943563, N=785,   spacing=8),
    7:  dict(r=37633,    lam=26672,  p=1152921504613142059, N=3136,  spacing=16),
    9:  dict(r=150589,   lam=100177, p=1152921504615061627, N=12549, spacing=32),
    11: dict(r=602137,   lam=117029, p=1152921504660984043, N=50178, spacing=64),
}
FRAC_SQR = 3
MU_WALK = 1.1
TARGET = 0.40
TABLE_ROWS = 64          # toy: 64 rows x 12 columns; production 2^14 x 12
ORBIT_BOUND = 8          # declared: no orbit may close below this many steps


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
    """The production step rule: canonical representative, inverse-free table.

    Only the *decision* is cached (it is a class function); the field element is
    recomputed every step, so g^a z^b == y holds at every step including the
    renormalisation.
    """

    def __init__(self, bits, seed):
        f = FIX[bits]
        self.bits, self.r, self.p = bits, f['r'], f['p']
        self.N, self.spacing = f['N'], f['spacing']
        self.T = TW.Tower(self.p)
        self.rng = random.Random(seed)
        e = (self.p ** 12 - 1) // self.r
        while True:
            y = self.T.pow(self.rand_gt(), e)
            if y != self.T.one12() and self.T.pow(y, self.r) == self.T.one12():
                self.g = y
                break
        self.ppow = [pow(self.p, i, self.r) for i in range(12)]
        self._build_frob_pow()
        self._build_table()
        self.cache = {}
        self.materialised = 0
        self.deg_kappa = 0

    def rand_gt(self):
        p, rng = self.p, self.rng
        return tuple(tuple((rng.randrange(p), rng.randrange(p)) for _ in range(3))
                     for _ in range(2))

    # ---- fast sigma^k: frob^k on the tower is (conjugate or not) times a constant
    def _build_frob_pow(self):
        """Precompute frob^k as one conjugation decision plus one Fp2 scaling.

        frob6(c)_i = conj(c_i) * xi3^i, so frob6^k(c)_i = E_k(c_i) * S[k][i] with
        E_k = conj when k is odd and identity when k is even.  S and the Fp12
        gamma-factor G_k = gamma^(1+p+...+p^(k-1)) are built by running the
        literal tower.frob and then asserted against it on random elements.
        """
        T, p = self.T, self.p
        self.S = [[None] * 3 for _ in range(12)]
        for k in range(12):
            for i in range(3):
                v = [(0, 0), (0, 0), (0, 0)]
                v[i] = (1, 0)
                v = (v[0], v[1], v[2])
                # literal frob6^k on limb i
                for _ in range(k):
                    v = T.frob6(v)
                self.S[k][i] = v[i]
        self.G = [T.one6[0]]
        for k in range(1, 12):
            self.G.append(T.f2mul(T.f2conj(self.G[-1]), T.gamma))

    def frob_pow(self, y, k):
        k %= 12
        if k == 0:
            return y
        T = self.T
        par = k & 1
        S = self.S[k]
        def f6(c):
            out = []
            for i in range(3):
                ci = c[i]
                if par:
                    ci = T.f2mul(T.f2conj(ci), S[i])
                else:
                    ci = T.f2mul(ci, S[i])
                out.append(ci)
            return (out[0], out[1], out[2])
        return (f6(y[0]), T.f6scale(f6(y[1]), self.G[k]))

    def check_frob_pow(self, n=24):
        """Differential: fast frob_pow vs the literal tower.frob applied k times."""
        ok = bad = 0
        for _ in range(n):
            y = self.T.pow(self.g, self.rng.randrange(1, self.r))
            for k in range(12):
                lit = y
                for _ in range(k):
                    lit = self.T.frob(lit)
                if self.frob_pow(y, k) == lit:
                    ok += 1
                else:
                    bad += 1
        return dict(checks=ok + bad, pass_=ok, fail=bad)

    # ---- table: independent lower-half exponents, no inverse is an entry
    def _build_table(self):
        r, rng = self.r, self.rng
        half = (r - 1) // 2
        self.u = [[rng.randrange(1, half + 1) for _ in range(12)]
                  for _ in range(TABLE_ROWS)]
        flat = [v for row in self.u for v in row]
        self.table_inverse_free = CW.table_is_inverse_free(flat, r)
        self.table_all_lower_half = all(1 <= v <= half for v in flat)
        self.tbl = [[self.T.pow(self.g, v) for v in row] for row in self.u]

    # ---- class key: the twelve keys rotated so the minimum is first
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
        if h % 8 < FRAC_SQR:
            rec = ('sq', None, None, None)
        else:
            rec = ('mul', (h >> 8) % TABLE_ROWS, (h >> 20) % 12, None)
        self.cache[cid] = rec
        self.materialised += 1
        return rec

    def step(self, y, cid, a, b):
        """One real field operation on the CANONICAL representative, then the
        renormalisation.  Returns (y', a', b', kind, k', kappa)."""
        kind, t, col, _ = self.decision(cid)
        T, r, p = self.T, self.r, self.p
        if kind == 'sq':
            raw = T.cyclo_sqr_fast(y)
            a2, b2 = (2 * a) % r, (2 * b) % r
            delta = None
        else:
            raw = T.mul(y, self.tbl[t][col])
            a2, b2 = (a + self.u[t][col]) % r, b % r
            delta = self.u[t][col]
        k, _, _ = self.kappa(raw)
        y2 = self.frob_pow(raw, k)
        pk = self.ppow[k]
        return y2, a2 * pk % r, b2 * pk % r, kind, k, delta


# ---------------------------------------------------------------- audits
def audit(bits, seed, reps=48):
    """A: class id Frobenius invariant.  B: the step is a class function.

    B is the equivariance audit in the form this design needs, and it is not the
    v10f form.  Because the state is held canonical, the step is a class function
    iff (B1) the renormalisation lands on ONE representative -- canon(sigma^j(r))
    == r for all twelve members of the class -- and (B2) the successor class
    computed from each member, through the shipped path, is the same class.  Both
    are counted over all twelve members.

    Control C is the v10f index rule (i = (tw-k) mod 12) evaluated on the NEW
    table: it splits, because the new table is not a sigma-ladder.  The shipped
    step does not use it; reporting the split is the evidence that the old rule
    is gone rather than merely unused.
    """
    cw = ClassWalk(bits, seed)
    T = cw.T
    A_ok = A_bad = B_ok = B_bad = E_ok = E_bad = 0
    C_ok = C_bad = 0
    for _ in range(reps):
        x = cw.rng.randrange(1, cw.r)
        z = T.pow(cw.g, x)
        a, b = cw.rng.randrange(cw.r), cw.rng.randrange(cw.r)
        y0 = T.mul(T.pow(cw.g, a), T.pow(z, b))
        k0, cid0, _ = cw.kappa(y0)
        yc = cw.frob_pow(y0, k0)                  # canonical representative
        succ_c, _, _, kind_c, _, _ = cw.step(yc, cid0, a, b)
        _, succ_cid, _ = cw.kappa(succ_c)
        yj = y0
        for j in range(1, 12):
            yj = T.frob(yj)
            kj, cidj, _ = cw.kappa(yj)
            if cidj == cid0:
                A_ok += 1
            else:
                A_bad += 1
            # B1: renormalsation lands on the same representative
            # B2: successor class through the shipped (canonicalising) path
            yjc = cw.frob_pow(yj, kj)
            sj, _, _, kindj, _, _ = cw.step(yjc, cidj, a, b)
            _, succj, _ = cw.kappa(sj)
            if yjc == yc and succj == succ_cid and kindj == kind_c:
                B_ok += 1
            else:
                B_bad += 1
            # C (control): the v10f (tw-k) index rule on the new table
            kind, t, tw, _ = cw.decision(cidj)
            if kind == 'mul':
                _, succ_old, _ = cw.kappa(T.mul(yj, cw.tbl[t][(tw - kj) % 12]))
                if succ_old == succ_cid:
                    C_ok += 1
                else:
                    C_bad += 1
            # E: the D1 test -- is the inverse of an applied multiplier also a
            # table entry?  In v12b the applied multiplier is a plain entry, so
            # this asks whether sigma^6(entry) is in the table (it never is).
            if kind == 'sq':
                E_ok += 1
            else:
                inv_m = cw.frob_pow(cw.tbl[t][tw], 6)
                in_table = any(inv_m == e for row in cw.tbl for e in row)
                if in_table:
                    E_bad += 1
                else:
                    E_ok += 1
    return dict(density_bits=bits, prime_r=cw.r,
                audit_A_key_invariance_checks=A_ok + A_bad, audit_A_pass=A_ok,
                audit_A_fail=A_bad,
                audit_B_equivariance_checks=B_ok + B_bad, audit_B_pass=B_ok,
                audit_B_fail=B_bad,
                audit_E_inverse_free_checks=E_ok + E_bad, audit_E_pass=E_ok,
                audit_E_fail=E_bad,
                control_C_old_rule_checks=C_ok + C_bad, control_C_old_rule_splits=C_bad,
                kappa_degenerate=cw.deg_kappa,
                table_entries=TABLE_ROWS * 12,
                table_inverse_free=cw.table_inverse_free,
                table_all_lower_half=cw.table_all_lower_half,
                frob_pow_differential=cw.check_frob_pow(),
                scope='A: class id vs the eleven conjugates. B1: canon(sigma^j r) == r. '
                      'B2: successor class and step type from each member through the '
                      'shipped canonicalising path vs from the canonical rep. E: '
                      'sigma^6(entry), the inverse under p^6 = -1, is never a table entry. '
                      'C: the v10f (tw-k) rule on the new table, as a control.')


def sweep(bits, trials, seed, inv_trials=4):
    """Four trails of w = max(S//4, 8*spacing) steps sharing one DP store.

    Every reported class is stored; a second, different label pair at the same
    class solves by x = (a1-a2)/(b2-b1) mod r with no k-search, and the recovered
    x is VERIFIED by exponentiation against the planted secret -- there is no
    seen[]-style proxy anywhere in this file.  A trail walks its whole length (it
    is not restarted at a report) and records exact state repeats.
    """
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
              inv_checks=0, inv_failures=0, renorm_checks=0, renorm_failures=0,
              repeated_class_revisits=0, exact_state_repeats=0, orbit_censor=0,
              min_orbit=0, self_repeat_diff=0)
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
            y = cw.frob_pow(y, k)                  # trail starts canonical ...
            a, b = a * cw.ppow[k] % r, b * cw.ppow[k] % r   # ... with its labels
            seen = {}
            for s_i in range(w_trail):
                y2, a2, b2, kind, kk, delta = cw.step(y, cid, a, b)
                inv = T.mul(T.pow(cw.g, a2), T.pow(z, b2)) == y2
                y, a, b = y2, a2, b2
                k, cid, m = cw.kappa(y)
                st['steps'] += 1
                if ti < inv_trials and seg == 0:
                    st['inv_checks'] += 1
                    st['renorm_checks'] += 1
                    if not inv:
                        st['inv_failures'] += 1
                        st['renorm_failures'] += 1
                # orbit census: exact state repeat (class id AND labels)
                key = (cid, a, b)
                if key in seen:
                    st['exact_state_repeats'] += 1
                    L = s_i - seen[key]
                    if L < ORBIT_BOUND:
                        st['orbit_censor'] += 1
                    st['min_orbit'] = L if st['min_orbit'] == 0 else min(st['min_orbit'], L)
                seen[key] = s_i
                if (m ^ (m >> 32)) % spacing == 0:
                    st['reports'] += 1
                    prev = store.get(cid)
                    if prev is None:
                        store[cid] = (a, b)
                    elif (prev[0], prev[1]) != (a, b):
                        st['repeated_class_revisits'] += 1
                        if seg == 0:
                            st['self_repeat_diff'] += 1
                        xg = CW.solve_x(prev[0], prev[1], a, b, r)
                        if xg is None:
                            st['degenerate_den'] += 1
                            continue
                        st['merges'] += 1
                        if T.pow(cw.g, xg) == z:
                            st['solve_success'] += 1
                            st['success'] += 1
                            ok = True
                        else:
                            st['false_confirm'] += 1
                if ok:
                    break
        if not ok:
            st['stranded_steps'] += budget
        st['trials'] += 1
    out = dict(density_bits=bits, prime_r=r, field_p=f['p'], lam=f['lam'], N=N,
               spacing=spacing, S_ref=S, w_trail=w_trail, n_trails=n_trails,
               budget=budget, table_rows=TABLE_ROWS, table_columns=12,
               frac_sqr='3/8', mu=MU_WALK, target=TARGET,
               orbit_bound=ORBIT_BOUND,
               classes_materialised=cw.materialised, kappa_degenerate=cw.deg_kappa,
               table_inverse_free=cw.table_inverse_free, **st)
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
