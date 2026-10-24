#!/usr/bin/env python3
"""transfer_arm.py -- v18: the fixture->production TRANSFER TEST (verdict cd18d785).

The last verdict rejects the package on one premise: the success law is measured on
a scalar quotient surrogate at q ~ 4.29e9 (N_classes 7.15e8) and the claim transfers
it to secp256k1 (N_classes ~ 2^254) with no support.  This script measures the law at
a LADDER OF ORDERS instead of asserting the transfer.

Dimensionless shape held fixed across arms:  x = visits / sqrt(N_classes), the
birthday scaling variable, at the shipped value 1.009444.  If success is a function
of x alone, it is free of the order and the production value is the measured value,
with the measured order-trend reported as a disclosed band.

ARMS  (R=64, d=5, W=16, same-challenge w=2 as pair_census.py)
    L24  q ~ 2^24   N ~ 2.8e6
    L28  q ~ 2^28   N ~ 4.5e7
    L32  q ~ 2^32   N ~ 7.2e8   the shipped surrogate (reproduction arm)
    L36  q ~ 2^36   N ~ 1.1e10  16x the shipped class count

The shipped qwalk.py keeps every word in np.uint64, exact only while q < 2^32 (the
product u*x must fit 2^64).  For L36 this file supplies a 12-bit-limb constant
multiply that stays inside uint64 for q < 2^44:

    k*c mod q = sum_i ((c >> 12i) & 0xFFF) * ((k << 12i) mod q)  mod q

every term < 2^12 * 2^36 = 2^48, three terms < 2^50, so no word overflows.  The
general product u[m]*x is a TABLE built once per run (R x n_sets Python-int
products, looked up by np.searchsorted over the sorted challenge block), and the
coefficient update multiplies by one of the six mu_6 units -- a CONSTANT -- so it
uses the same limb path.

VALIDATION, mandatory before any extended arm: `validate_equiv` drives the shipped
qwalk step and this path side by side at q ~ 2^32 and asserts identical canonical
classes, DP hits and escape fires on every step of a 400-step trace.  A failure
stops the ladder at 2^32 (evidence/transfer-equiv.json records it).

Each arm writes evidence/transfer-arm-<name>.json when it finishes.
"""
import hashlib
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk                     # noqa: E402
import census_v20n as C          # noqa: E402

M12 = (1 << 12) - 1
OUTDIR = os.path.join(D, "evidence")


class Universe64(qwalk.Universe):
    """qwalk.Universe with a uint64-exact modular path for q < 2^44."""

    def __init__(self, q_min, R, seed=20250, row_mode="orbit"):
        super().__init__(q_min, R, seed, row_mode)
        assert self.q < (1 << 44), self.q
        self.nl = max(1, (self.q.bit_length() + 11) // 12)
        self._c = {}

    def cmul(self, k, c):
        """k*c mod q, k a CONSTANT, c a uint64 vector.  Exact for q < 2^44."""
        terms = self._c.get(k)
        if terms is None:
            terms = [((k << (12 * i)) % self.q, 12 * i) for i in range(self.nl)]
            self._c[k] = terms
        acc = None
        for w, sh in terms:
            t = ((c >> np.uint64(sh)) & np.uint64(M12)) * np.uint64(w)
            acc = t if acc is None else acc + t
        return (acc % np.uint64(self.q)).astype(np.uint64)

    def canon(self, c):
        """Same signature and same value as qwalk.Universe.canon."""
        q = np.uint64(self.q)
        c = np.asarray(c, dtype=np.uint64)
        u1 = self.cmul(self.lam, c)
        u2 = self.cmul(self.lam2, c)
        w = np.minimum(np.minimum(c, u1), u2)
        W = np.maximum(np.maximum(c, u1), u2)
        pos = w <= q - W
        cn = np.where(pos, w, q - W)
        hp = np.where(c <= u1, np.where(c <= u2, 0, 2), np.where(u1 <= u2, 1, 2))
        hq = np.where(c >= u1, np.where(c >= u2, 0, 2), np.where(u1 >= u2, 1, 2))
        e = np.where(pos, self.LAMH[hp], q - self.LAMH[hq])
        return cn.astype(np.uint64), e.astype(np.uint64), hp.astype(np.int64)

    def scale_e(self, e, v):
        """e*v mod q for e in the six mu_6 units (a vector), v arbitrary."""
        q = np.uint64(self.q)
        v = np.asarray(v, dtype=np.uint64)
        lam = np.uint64(self.lam)
        lam2 = np.uint64(self.lam2)
        p1 = self.cmul(self.lam, v)
        p2 = self.cmul(self.lam2, v)
        one = e == np.uint64(1)
        two = e == lam
        three = e == lam2
        mone = e == q - np.uint64(1)
        ml = e == q - lam
        ml2 = e == q - lam2
        r = np.where(one, v, np.where(two, p1, p2))
        r = np.where(mone, (q - v) % q, r)
        r = np.where(ml, (q - p1) % q, r)
        r = np.where(ml2, (q - p2) % q, r)
        assert bool(np.all(one | two | three | mone | ml | ml2)), "e not a mu_6 unit"
        return r.astype(np.uint64)


class Ext:
    """Shim so census_v20n.walk_arm can drive the extended path unchanged."""

    Universe = Universe64
    K_ESC = qwalk.K_ESC
    K_ESC_MAX = qwalk.K_ESC_MAX
    COUNTER_ONLY = False
    sm64 = staticmethod(qwalk.sm64)

    @staticmethod
    def _inc(U, m, x):
        rows = np.searchsorted(U.xt_sorted, np.asarray(x, dtype=np.uint64))
        assert bool(np.all(U.xt_sorted[rows] == np.asarray(x, dtype=np.uint64)))
        return U.inc_sorted[m, rows]

    @staticmethod
    def _apply(U, canon, a, cc, x, m):
        q = np.uint64(U.q)
        u = U.u[m]
        t = canon + Ext._inc(U, m, x)
        t = np.where(t >= q, t - q, t)
        cn, e, _hp = U.canon(t)
        a = U.scale_e(e, a)
        cc2 = cc + u
        cc2 = np.where(cc2 >= q, cc2 - q, cc2)
        cc = U.scale_e(e, cc2)
        return cn, a, cc, e

    @staticmethod
    def step(U, canon, a, cc, x, d, want_dp=False, prev=None, esc=False,
             force_every=0, tick=0, esc_run=None, kcap=None, restart_out=None,
             coins=None):
        m = U.col(canon)
        cn, a2, cc2, _e = Ext._apply(U, canon, a, cc, x, m)
        fired = np.zeros(canon.size, dtype=bool)
        if esc and prev is not None:
            fire = cn == prev
            if force_every:
                fire = fire | (((np.arange(canon.size) + tick) % force_every) == 0)
            if fire.any():
                idx = np.flatnonzero(fire)
                ma = (m[idx] + 1) % U.R
                cn2, a3, cc3, _e2 = Ext._apply(U, canon[idx], a[idx], cc[idx], x[idx], ma)
                cn = cn.copy(); a2 = a2.copy(); cc2 = cc2.copy()
                cn[idx] = cn2; a2[idx] = a3; cc2[idx] = cc3
                fired = fire
        if esc_run is not None:
            esc_run[:] = np.where(fired, esc_run + 1, 0)
            if kcap is not None:
                if restart_out is not None:
                    restart_out[:] = False
                over = esc_run > kcap
                if int(np.count_nonzero(over)):
                    idx = np.flatnonzero(over)
                    if coins is None:
                        raise AssertionError("kcap breach must resample coins")
                    cn3, a4, cc4, _a0, _c0 = Ext.fresh_state(U, x[idx], coins, idx.size)
                    cn = cn.copy(); a2 = a2.copy(); cc2 = cc2.copy()
                    cn[idx] = cn3; a2[idx] = a4; cc2[idx] = cc4
                    esc_run[idx] = 0
        hit = U.dp(cn, d) if want_dp else None
        prev_out = canon.copy()
        if restart_out is not None and restart_out.any():
            prev_out = prev_out.copy()
            prev_out[restart_out] = cn[restart_out]
        return cn, a2.astype(np.uint64), cc2.astype(np.uint64), hit, prev_out, fired

    @staticmethod
    def fresh_state(U, x, coins, n):
        a0 = coins.integers(0, U.q, size=n).astype(np.uint64)
        c0 = coins.integers(0, U.q, size=n).astype(np.uint64)
        t = np.array([(int(ai) + int(ci) * int(xi)) % U.q
                      for ai, ci, xi in zip(a0.tolist(), c0.tolist(),
                                            np.asarray(x).tolist())], dtype=np.uint64)
        cn, e, _hp = U.canon(t)
        return cn, U.scale_e(e, a0), U.scale_e(e, c0), a0, c0


def validate_equiv(steps=400, seed=20250):
    U0 = qwalk.Universe(4_290_000_000, 64, seed=seed)
    U1 = Universe64(4_290_000_000, 64, seed=seed)
    assert U0.q == U1.q
    rng = np.random.default_rng(seed + 1)
    lanes = 48
    x = np.repeat(rng.integers(1, U0.q, size=lanes // 16).astype(np.uint64), 16)
    a = rng.integers(0, U0.q, size=lanes).astype(np.uint64)
    cc = rng.integers(0, U0.q, size=lanes).astype(np.uint64)
    c0 = U0.canon(a)[0]
    c1 = U1.canon(a)[0]
    assert np.array_equal(c0, c1), "canon mismatch at start"
    U1.xt_sorted = np.sort(np.unique(x))
    U1.inc_sorted = np.array([[int(v) % U1.q for v in
                               (U1.u.astype(object)[:, None] *
                                U1.xt_sorted.astype(object)[None, :]).tolist()[i]]
                              for i in range(U1.R)], dtype=np.uint64)
    prev0 = c0.copy()
    prev1 = c1.copy()
    nf = 0
    for t in range(steps):
        c0, a, cc, hit0, prev0, f0 = qwalk.step(
            U0, c0, a, cc, x, 5, True, prev0, True, 0, t, None, None, None, None)
        c1, a, cc, hit1, prev1, f1 = Ext.step(
            U1, c1, a, cc, x, 5, True, prev1, True, 0, t, None, None, None, None)
        assert np.array_equal(c0, c1), ("class divergence at step %d" % t)
        assert np.array_equal(hit0, hit1), ("dp divergence at step %d" % t)
        assert np.array_equal(f0, f1), ("fire divergence at step %d" % t)
        nf += int(f0.sum())
    return {"test": "equiv_extended_vs_shipped", "q": int(U0.q), "steps": steps,
            "lanes": lanes, "class_divergences": 0, "dp_divergences": 0,
            "fire_divergences": 0, "trace_equal": True, "fires_in_trace": nf,
            "note": ("shipped qwalk.step and the 12-bit-limb extended step agree "
                     "on canonical class, DP gate and escape fire at every step")}


def run_arm(name, q_min, sets, R=64, d=5, W=16, seed=20250, use_ext=False, shape=None):
    shape = C.SHAPE if shape is None else shape
    t0 = time.time()
    U = Universe64(q_min, R, seed=seed) if use_ext else qwalk.Universe(q_min, R, seed=seed)
    shim = Ext if use_ext else qwalk
    per = int(round(float(shape) * math.sqrt(U.classes) / W))
    xt, _xv0, _a0, _c0, _cn0 = C.make_start(U, sets, W, seed)
    xt_sha = hashlib.sha256(xt.tobytes()).hexdigest()
    if use_ext:
        U.xt_sorted = np.sort(xt)
        U.inc_sorted = np.array(
            [[int(v) % U.q for v in row] for row in
             (U.u.astype(object)[:, None] * U.xt_sorted.astype(object)[None, :]).tolist()],
            dtype=np.uint64)
    saved = C.qwalk
    C.qwalk = shim
    try:
        sols = []
        fires3 = inv_bad = restarts = 0
        for r in (0, 1):
            rng = np.random.default_rng((seed * 1000 + r + 1) & 0x7FFFFFFF)
            xv = np.repeat(xt, W).astype(np.uint64)
            lanes = xt.size * W
            aa = rng.integers(0, U.q, size=lanes).astype(np.uint64)
            cq = rng.integers(0, U.q, size=lanes).astype(np.uint64)
            t0v = np.array([(int(ai) + int(ci) * int(xi)) % U.q
                            for ai, ci, xi in zip(aa.tolist(), cq.tolist(), xv.tolist())],
                           dtype=np.uint64)
            cn, e, _hp = U.canon(t0v) if use_ext else (U.canon(t0v)[0], U.canon(t0v)[1], 0)
            if use_ext:
                a_all = U.scale_e(e, aa)
                c_all = U.scale_e(e, cq)
            else:
                a_all = (e * aa) % U.q
                c_all = (e * cq) % U.q
            res = C.walk_arm(U, cn, a_all, c_all, xv, d, per, True, 3, W, sets,
                             coins_seed=0xC0FFEE + r)
            sols.append(res["succ_step"] >= 0)
            fires3 += res["fires_lag3"]
            inv_bad += res["invariant_violations"]
            restarts += res["restarts"]
    finally:
        C.qwalk = saved
    A, B = sols
    sset = A | B
    n = sets
    kA, kB, kS = int(A.sum()), int(B.sum()), int(sset.sum())
    p = (kA + kB) / (2 * n)
    out = {"arm": name, "q_min": q_min, "q": int(U.q), "q_log2": U.q.bit_length() - 1,
           "classes": int(U.classes), "classes_log2": math.log2(U.classes),
           "R": R, "dp_bits": d, "walks_per_set": W, "sets": n, "per_walker_steps": per,
           "visits_per_set": W * per,
           "visits_over_sqrtN": W * per / math.sqrt(U.classes),
           "extended_arithmetic": bool(use_ext),
           "challenge_block_sha256": xt_sha, "seed": seed,
           "runs": {"A_solved": kA, "B_solved": kB, "both": int((A & B).sum()),
                    "set_solved": kS, "sets": n},
           "p_run_pooled": p, "p_set_measured": kS / n,
           "p_set_formula_2p_minus_p2": 2 * p - p ** 2,
           "cp95_set_lower": str(C.cp95_lower(kS, n)),
           "cp95_set_lower_f": float(C.cp95_lower(kS, n)),
           "cp95_run_A_lower": str(C.cp95_lower(kA, n)),
           "lag3_fires": fires3, "restarts": restarts,
           "invariant_violations_sampled": inv_bad,
           "elapsed_s": round(time.time() - t0, 1), "floor": 0.39,
           "floor_met_by_cp95": bool(float(C.cp95_lower(kS, n)) >= 0.39)}
    with open(os.path.join(OUTDIR, "transfer-arm-%s.json" % name), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("%-4s q=2^%-3d N=2^%6.2f per=%-6d sets=%-5d p_run %.5f p_set %.5f cp95 "
          "%.5f  %.0fs" % (name, out["q_log2"], out["classes_log2"], per, n, p,
                           out["p_set_measured"], out["cp95_set_lower_f"],
                           out["elapsed_s"]), flush=True)
    return out


SPEC = {
    "L24": dict(q_min=1 << 24, sets=1600),
    "L28": dict(q_min=1 << 28, sets=1200),
    "L32": dict(q_min=4_290_000_000, sets=800),
    "L36": dict(q_min=1 << 36, sets=400, use_ext=True),
    "L40": dict(q_min=1 << 40, sets=128, use_ext=True),
}


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    which = os.environ.get("ARMS", "L24,L28,L32,L36").split(",")
    if any(SPEC[a].get("use_ext") for a in which):
        eq = validate_equiv()
        with open(os.path.join(OUTDIR, "transfer-equiv.json"), "w") as f:
            json.dump(eq, f, indent=1, sort_keys=True)
        print("equiv ok:", eq["q"], eq["steps"], "fires", eq["fires_in_trace"], flush=True)
        assert eq["trace_equal"]
    for a in which:
        run_arm(a, **SPEC[a])


if __name__ == "__main__":
    main()
