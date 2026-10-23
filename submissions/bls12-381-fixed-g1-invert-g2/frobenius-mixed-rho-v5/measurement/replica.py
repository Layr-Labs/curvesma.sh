"""Python replica of step.c (walk-spec section 6, v10): reruns the blst GT walk and compares the kernel's dump exactly.

  python3 replica.py <input-dir> <rounds> [--direct N] [--cap C]
  python3 replica.py <input-dir> 0 --replay [--pairs N] [--t2-steps K] [--workers W] [--pool-rounds R]

<input-dir> holds inputs.bin and instance.json from gen_inputs.py, plus the kernel's
dump-<rounds>.txt and native-<rounds>.stderr when present (measure.sh output). Writes
replica-<rounds>.txt next to them and prints a report.

--cap C sets the trail cap of a test build (step.c compiled with -DCAP_STEPS_OVERRIDE=C); the
production cap is 2^66 (spec section 6).

Per round, each walk (v10): produce the new point; the DP test first (record, restart, chained while the restarted
state is a DP); then the detector (2-cycle, anchor): a detected cycle is escaped to c_s * Y_j (the least-key state s,
its canonical representative c_s = sigma^kn(s), canonical exponents, j = mix(low64(K_s)) mod TR); the trail cap restarts
a trail that runs CAP steps without a DP. The anchor refreshes on the global round unless the walk restarted.

Checks:
  1. the walk reproduces the kernel's dump line for line, and the event counts match the kernel's final stderr line;
  2. the 4-multiplication key (SPEC section 2) equals the direct minimum over the 12 conjugates (--direct N states);
  3. exponent invariant: every final state equals g^(a + y b), with y the instance exponent;
  4. full-visit check: every restart output key is compared with all visited keys (main-step states of every walk and
     the start states) and with all earlier restart outputs; both counts must be 0;
  5. --replay: replay check (T1 first-DP equality, T2 detection and escape synchronisation) of two walks at the same
     element with different detector states and exponents (DP_BITS must be 2).
"""
import argparse
import json
import multiprocessing
import os
import random
import struct
import sys

from py_ecc.optimized_bls12_381 import FQ12, G1, G2, curve_order, field_modulus, multiply, pairing

P = field_modulus
R = curve_order
MONT = (1 << 384) % P
ANCHOR_PERIOD = 64
MAX_CYCLE = 64
MAGIC = 0x3130304f48525447
STRUCT_ORDER = [0, 2, 4, 1, 3, 5]
GAMMA = 0x9E3779B97F4A7C15  # restart index multiplier (spec section 6)
TR = 1024  # restart table size (spec section 6, v10)
M64 = (1 << 64) - 1
DP_RING = 65536  # distinguished-point ring size of step.c


def mix64(z):
    """mix of spec section 6, all arithmetic mod 2^64."""
    z &= M64
    z ^= z >> 30
    z = (z * 0xBF58476D1CE4E5B9) & M64
    z ^= z >> 27
    z = (z * 0x94D049BB133111EB) & M64
    z ^= z >> 31
    return z


class Inputs:
    def __init__(self, path):
        data = open(path, "rb").read()
        self.off = 0
        self.data = data
        assert self.u64() == MAGIC and self.u64() == 5, "bad magic or version"
        self.RA = self.u64()
        self.ADDER_BITS = self.u64()
        self.DP_BITS = self.u64()
        self.S = self.u64()
        self.NW = self.u64()
        assert self.u256() == R, "group order mismatch"
        self.LP = [self.u256() for _ in range(12)]  # lambda^k
        self.DP_MASK = (1 << self.DP_BITS) - 1
        self.ZETA = self.fp()
        self.G0 = self.fp()
        self.G1 = self.fp()
        self.FAM, self.JJ, self.EPS = [], [], []
        for _ in range(12):
            self.FAM.append(self.u64())
            self.JJ.append(self.u64())
            self.EPS.append(self.u64())
        self.TAB, self.EC, self.ED = [], [], []
        for _ in range(self.RA * 12):
            self.TAB.append(self.fp12())
            self.EC.append(self.u256())
            self.ED.append(self.u256())
        self.walks = []
        for _ in range(self.NW):
            x = self.fp12()
            a = self.u256()
            b = self.u256()
            self.walks.append((x, a, b))
        # restart table block (walk-spec section 6, v10): TR, then Y_j = g^(f_j), then f_j
        assert self.u64() == TR, "restart table size"
        self.RY = [self.fp12() for _ in range(TR)]
        self.RF = [self.u256() for _ in range(TR)]

    def u64(self):
        v = struct.unpack_from("<Q", self.data, self.off)[0]
        self.off += 8
        return v

    def u256(self):
        return sum(self.u64() << (64 * i) for i in range(4))

    def fp(self):
        v = int.from_bytes(self.data[self.off:self.off + 48], "big")
        self.off += 48
        return v

    def fp12(self):
        """12 Fp in blst struct order -> py_ecc coefficients c_0..c_11."""
        A = [None] * 6
        for i in STRUCT_ORDER:
            re = self.fp()
            im = self.fp()
            A[i] = (re, im)
        c = [0] * 12
        for i in range(6):
            re, im = A[i]
            c[i + 6] = im
            c[i] = (re - im) % P
        return FQ12(c)


_FROB = None


def frob_matrix():
    """Matrix of the p-power map on the py_ecc coefficients: sigma(c) = M c (as in gen_inputs.py)."""
    global _FROB
    if _FROB is None:
        cols = []
        for i in range(12):
            w_i = FQ12([0] * i + [1] + [0] * (11 - i))
            cols.append((w_i ** P).coeffs)
        _FROB = [[cols[c][r] % P for c in range(12)] for r in range(12)]
    return _FROB


def sigma_pow(x, n):
    """sigma^n(x) = x^(p^n), n < 12 (spec section 2)."""
    if n == 0:
        return x
    M = frob_matrix()
    c = list(x.coeffs)
    for _ in range(n):
        c = [sum(M[r][k] * c[k] for k in range(12)) % P for r in range(12)]
    return FQ12(c)


def struct_hex(x):
    """py_ecc FQ12 -> the 12 canonical Fp of blst struct order, as 96-char hex strings."""
    c = x.coeffs
    out = []
    A = {i: ((c[i] + c[i + 6]) % P, c[i + 6] % P) for i in range(6)}
    for i in STRUCT_ORDER:
        re, im = A[i]
        out.append("%096x" % re)
        out.append("%096x" % im)
    return out


def fresh(s, phase, rnd, restarted):
    """The state after a restart or escape: detector cleared; the anchor refreshes on the round boundary unless the
    walk restarted in this round (spec section 6)."""
    st = dict(s, prevK=None, prevValid=False, anchorK=None, anchorValid=False, age=0)
    if not restarted and (rnd - phase + 1) % ANCHOR_PERIOD == 0:
        st["anchorK"] = s["K"]
        st["anchorValid"] = True
    return st


class Walk:
    def __init__(self, inp, cap=1 << 66, track=True):
        self.inp = inp
        self.RA = inp.RA
        self.S = inp.S
        self.dp_mask = inp.DP_MASK
        self.cap = cap  # trail cap in steps (production 2^66)
        self.track = track  # full-visit bookkeeping (main runs)
        # steps: main walk steps; adds, sqrs: every group operation (cycle checks, escapes' walk-throughs included)
        # restarts: at a DP (chained ones included); cap_restarts: at the trail cap; chained: restarts to a DP state
        # escapes: one per detected cycle; escape_dps: escape targets that are DPs
        self.counts = dict(steps=0, adds=0, sqrs=0, dps=0, cyc2=0, cycl=0, restarts=0, cap_restarts=0, chained=0,
                           escapes=0, escape_dps=0)
        self.dp_lines = []  # distinguished-point records in order: fingerprint, canonical exponents
        self.rc = [0] * inp.NW  # per-walk restart counter t (spec section 6)
        self.seen_walk = [set() for _ in range(inp.NW)] if track else None
        self.seen_all = set() if track else None
        self.rep = dict(walk=0, all=0)
        self.visited = set()
        self.restart_keys = set()
        self.fv = dict(restart_in_visited=0, restart_repeat=0)
        self.last_event = None
        self.last_dp = None
        self.last_least = None
        self.last_out = None

    def key(self, x):
        """4-multiplication formula (SPEC section 2). Returns (mont, canonical value, index)."""
        inp = self.inp
        c = x.coeffs
        a0 = (c[1] + c[7]) % P
        a1 = c[7] % P
        B = (a0 * inp.G0 + a1 * inp.G1) % P
        az = a0 * inp.ZETA % P
        bz = B * inp.ZETA % P
        fam = [[a0, az, (az - a0) % P], [B, bz, (bz - B) % P]]
        best = None
        for n in range(12):
            v = fam[inp.FAM[n]][inp.JJ[n]]
            if inp.EPS[n]:
                v = (-v) % P
            m = v * MONT % P
            if best is None or m < best[0]:
                best = (m, v, n)
        return best

    def step(self, s):
        """s = dict(x, a, b, K, Kc, kn). Mirrors step() in step.c (one group operation)."""
        inp = self.inp
        h = s["K"] & (self.RA - 1)
        if h < self.S:
            x2 = s["x"] * s["x"]
            a2 = 2 * s["a"] % R
            b2 = 2 * s["b"] % R
            self.counts["sqrs"] += 1
        else:
            idx = h * 12 + s["kn"]
            x2 = s["x"] * inp.TAB[idx]
            a2 = (s["a"] + inp.EC[idx]) % R
            b2 = (s["b"] + inp.ED[idx]) % R
            self.counts["adds"] += 1
        m, v, n = self.key(x2)
        return dict(x=x2, a=a2, b=b2, K=m, Kc=v, kn=n)

    def escape(self, s, w):
        """Escape from the least-key state s (spec section 6, v10): x <- sigma^kn(s) Y_j, a <- lambda^kn a_s + f_j,
        b <- lambda^kn b_s, j = mix(low64(K_s)) mod TR. No restart counter is used."""
        inp = self.inp
        self.last_least = s
        cs = sigma_pow(s["x"], s["kn"])
        j = mix64(s["K"] & M64) & (TR - 1)
        x = cs * inp.RY[j]
        lam = inp.LP[s["kn"]]
        A = s["a"] * lam % R
        B = s["b"] * lam % R
        m, v, n = self.key(x)
        self.counts["escapes"] += 1
        return dict(x=x, a=(A + inp.RF[j]) % R, b=B, K=m, Kc=v, kn=n)

    def dp_hit(self, K):
        return (K >> self.inp.ADDER_BITS) & self.dp_mask == 0

    def restart(self, s, w):
        """Restart of walk w (spec section 6, v10): t = the walk's counter after increment, j = mix((w + 1) GAMMA + t)
        mod TR; the restart element Y_j of the restart table: x <- x Y_j, a <- a + f_j, b unchanged."""
        inp = self.inp
        self.rc[w] += 1
        t = self.rc[w]
        j = mix64(((w + 1) * GAMMA + t) & M64) & (TR - 1)
        x = s["x"] * inp.RY[j]
        m, v, n = self.key(x)
        out = dict(x=x, a=(s["a"] + inp.RF[j]) % R, b=s["b"], K=m, Kc=v, kn=n)
        if self.track:
            sig = (tuple(x.coeffs), out["a"], out["b"])
            if sig in self.seen_walk[w]:
                self.rep["walk"] += 1
            if sig in self.seen_all:
                self.rep["all"] += 1
            self.seen_walk[w].add(sig)
            self.seen_all.add(sig)
            if m in self.visited:
                self.fv["restart_in_visited"] += 1
            if m in self.restart_keys:
                self.fv["restart_repeat"] += 1
            self.restart_keys.add(m)
            self.visited.add(m)
        return out

    def dp_chain(self, s, w):
        """A distinguished point s: store its record, restart; chained while the restarted state is a DP."""
        cur = s
        while True:
            self.counts["dps"] += 1
            m = cur["K"]
            lam = self.inp.LP[cur["kn"]]
            # fingerprint: Montgomery key bits 70..197 (bits 10..69 are zero for every distinguished point)
            self.dp_lines.append("dp %016x %016x %064x %064x" % (
                (m >> 70) & (2**64 - 1), (m >> 134) & (2**64 - 1), cur["a"] * lam % R, cur["b"] * lam % R))
            cur = self.restart(cur, w)
            self.counts["restarts"] += 1
            if not self.dp_hit(cur["K"]):
                break
            self.counts["chained"] += 1
        return cur

    def advance(self, c, w, rnd, phase=0, allow_dp=True):
        """One round of walk w (spec sections 5-6, v10). c: state dict with the detector fields. Returns the new state.
        self.last_event: None, 'dp', 'cyc2' or 'cycl' (the event of this round)."""
        self.last_event = None
        nxt = self.step(c)
        self.counts["steps"] += 1
        if self.track:
            self.visited.add(nxt["K"])
        if allow_dp and self.dp_hit(nxt["K"]):
            self.last_event = "dp"
            self.last_dp = nxt
            return fresh(self.dp_chain(nxt, w), phase, rnd, True)
        two = False
        if c["prevValid"] and c["prevK"] == nxt["K"]:
            nn = self.step(nxt)
            if nn["K"] == c["K"]:
                two = True
        ev = None
        if two:
            self.counts["cyc2"] += 1
            self.last_event = "cyc2"
            ev = self.escape(nxt if nxt["K"] < c["K"] else c, w)
        elif c["anchorValid"] and c["anchorK"] == nxt["K"]:
            mn = nxt
            t = self.step(nxt)
            length = 1
            while t["K"] != nxt["K"] and length <= MAX_CYCLE:
                if t["K"] < mn["K"]:
                    mn = t
                t = self.step(t)
                length += 1
            if t["K"] == nxt["K"]:
                self.counts["cycl"] += 1
                self.last_event = "cycl"
                ev = self.escape(mn, w)
        if ev is not None:
            self.last_out = ev
            if self.track:
                self.visited.add(ev["K"])
            if allow_dp and self.dp_hit(ev["K"]):
                self.counts["escape_dps"] += 1
                return fresh(self.dp_chain(ev, w), phase, rnd, True)
            return fresh(ev, phase, rnd, False)
        age = c["age"] + 1
        if age >= self.cap:
            out = self.restart(nxt, w)
            self.counts["cap_restarts"] += 1
            if allow_dp and self.dp_hit(out["K"]):
                out = self.dp_chain(out, w)
            return fresh(out, phase, rnd, True)
        st = dict(nxt, prevK=c["K"], prevValid=True, anchorK=c["anchorK"], anchorValid=c["anchorValid"], age=age)
        if (rnd - phase + 1) % ANCHOR_PERIOD == 0:
            st["anchorK"] = nxt["K"]
            st["anchorValid"] = True
        return st

    def initial_states(self):
        states = []
        for w, (x, a, b) in enumerate(self.inp.walks):
            m, v, n = self.key(x)
            states.append(dict(x=x, a=a, b=b, K=m, Kc=v, kn=n, prevK=None, prevValid=False, anchorK=None,
                               anchorValid=False, age=0))
            if self.track:
                self.visited.add(m)
        return states

    def run(self, rounds, pool=None, pool_every=4):
        inp = self.inp
        states = self.initial_states()
        for rnd in range(rounds):
            for w in range(inp.NW):
                states[w] = self.advance(states[w], w, rnd)
            if pool is not None and rnd % pool_every == 0:
                for s in states:
                    pool.append(dict(x=s["x"], a=s["a"], b=s["b"], K=s["K"], kn=s["kn"]))
        return states


def dump_lines(states):
    lines = []
    for w, s in enumerate(states):
        a = "%064x" % s["a"]
        b = "%064x" % s["b"]
        parts = ["%d" % w, a, b, "%d" % s["kn"], "%096x" % s["Kc"]] + struct_hex(s["x"])
        lines.append(" ".join(parts))
    return lines


def parse_counts(path):
    last = [l for l in open(path).read().splitlines() if "=" in l][-1]
    out = {}
    for tok in last.split():
        k, _, v = tok.partition("=")
        out[k] = int(v)
    return out


def direct_key(x):
    """Direct min over V_n(sigma^n x), sigma^n x = x^(p^n), by exponentiation."""
    best = None
    for n in range(12):
        y = x if n == 0 else x ** (P ** n)
        v = (y.coeffs[1] + y.coeffs[7]) % P
        m = v * MONT % P
        if best is None or m < best[0]:
            best = (m, v, n)
    return best


# ---------------------------------------------------------------------------------------------------------------------
# Replay check (v10): two walks at the same element with different detector states and exponents.
REPLAY = {}


def _elem_eq(x1, x2):
    return tuple(x1.coeffs) == tuple(x2.coeffs)


def _pair_states(walk, pool, rng, y):
    """Two walk states at the same element: pool element, different exponents (a + b y fixed), different detectors."""
    pa = pool[rng.randrange(len(pool))]
    bB = rng.randrange(1, R)
    if bB == pa["b"]:
        bB = (bB + 1) % R
    aB = (pa["a"] + (pa["b"] - bB) * y) % R
    A = dict(x=pa["x"], a=pa["a"], b=pa["b"], K=pa["K"], kn=pa["kn"], Kc=0)
    B = dict(x=pa["x"], a=aB, b=bB, K=pa["K"], kn=pa["kn"], Kc=0)
    phiA = rng.randrange(ANCHOR_PERIOD)
    phiB = (phiA + rng.randrange(1, ANCHOR_PERIOD)) % ANCHOR_PERIOD
    for st, phi in ((A, phiA), (B, phiB)):
        st.update(prevK=pool[rng.randrange(len(pool))]["K"], prevValid=True,
                  anchorK=pool[rng.randrange(len(pool))]["K"], anchorValid=True, age=0)
        st["phi"] = phi
    return A, B


def replay_pair(i):
    """T1 (first-DP equality, DP test on) and T2 (detection and escape synchronisation, DP test off) for pair i."""
    ctx = REPLAY
    walk, pool, y = ctx["walk"], ctx["pool"], ctx["y"]
    rng = random.Random(ctx["seed"] * 1000003 + i)
    A0, B0 = _pair_states(walk, pool, rng, y)
    res = dict(i=i)
    # ---- T1: lockstep with the DP test on, until the first DP or detection in either walk
    A, B = dict(A0), dict(B0)
    phiA, phiB = A0["phi"], B0["phi"]
    t1 = "none"
    t1_mismatch = 0
    t1_round = None
    for rnd in range(ctx["t1_limit"]):
        A = walk.advance(A, 0, rnd, phase=phiA, allow_dp=True)
        evA = walk.last_event
        dpA = walk.last_dp
        B = walk.advance(B, 1, rnd, phase=phiB, allow_dp=True)
        evB = walk.last_event
        dpB = walk.last_dp
        if evA is None and evB is None:
            if not _elem_eq(A["x"], B["x"]):
                t1_mismatch += 1
            continue
        t1_round = rnd
        if evA == "dp" and evB == "dp":
            t1 = "both_dp"
            if not _elem_eq(dpA["x"], dpB["x"]):
                t1_mismatch += 1
        elif evA == "dp" or evB == "dp":
            t1 = "one_dp"
            t1_mismatch += 1
        else:
            t1 = "detect"
        break
    res.update(t1=t1, t1_mismatch=t1_mismatch, t1_round=t1_round)
    # ---- T2: lockstep with the DP test off, for t2_steps rounds; pair the detections of the two walks
    A, B = dict(A0), dict(B0)
    phi = {"A": A0["phi"], "B": B0["phi"]}
    ev_lists = {"A": [], "B": []}
    t2_mismatch = 0
    for rnd in range(ctx["t2_steps"]):
        for name, st, w in (("A", A, 0), ("B", B, 1)):
            nst = walk.advance(st, w, rnd, phase=phi[name], allow_dp=False)
            if walk.last_event in ("cyc2", "cycl"):
                least, out = walk.last_least, walk.last_out
                ev_lists[name].append(dict(
                    round=rnd, kind=walk.last_event, Ks=least["K"], kn=least["kn"],
                    la=least["a"], lb=least["b"], x=tuple(out["x"].coeffs), a=out["a"], b=out["b"]))
            if name == "A":
                A = nst
            else:
                B = nst
    # pair the detections of the same cycle (same least key)
    byK = {}
    for e in ev_lists["B"]:
        byK.setdefault(e["Ks"], []).append(e)
    matched = 0
    zero_off = 0
    elem_mismatch = 0
    elem_offset_mismatch = 0
    canon_mismatch = 0
    unmatched_A = 0
    for e in ev_lists["A"]:
        lst = byK.get(e["Ks"])
        if not lst:
            unmatched_A += 1
            continue
        f = lst.pop(0)
        matched += 1
        if e["x"] != f["x"]:
            elem_mismatch += 1
        da = (e["a"] - f["a"]) % R
        db = (e["b"] - f["b"]) % R
        if da == 0 and db == 0:
            zero_off += 1
        if (da + db * y) % R != 0:
            elem_offset_mismatch += 1
        # the escape adds the same f_j to both: the offset equals the offset of the canonical least states
        canon_a = (walk.inp.LP[e["kn"]] * e["la"] - walk.inp.LP[f["kn"]] * f["la"]) % R
        canon_b = (walk.inp.LP[e["kn"]] * e["lb"] - walk.inp.LP[f["kn"]] * f["lb"]) % R
        if (canon_a, canon_b) != (da, db):
            canon_mismatch += 1
    res.update(t2=dict(
        steps=ctx["t2_steps"], detections_A=len(ev_lists["A"]), detections_B=len(ev_lists["B"]),
        matched=matched, unmatched_A=unmatched_A, unmatched_B=sum(len(v) for v in byK.values()),
        escape_element_mismatches=elem_mismatch, zero_offsets=zero_off,
        offset_element_mismatches=elem_offset_mismatch, offset_vs_canonical_mismatches=canon_mismatch))
    res["t2_mismatch"] = elem_mismatch + zero_off + elem_offset_mismatch + canon_mismatch
    return res


def replay_main(inp, meta, args):
    global REPLAY
    assert inp.DP_BITS == 2, "replay check needs DP_BITS 2 inputs"
    walk_main = Walk(inp, track=True)
    pool = []
    walk_main.run(args.pool_rounds, pool=pool, pool_every=4)
    walk = Walk(inp, track=False)
    REPLAY = dict(walk=walk, pool=pool, y=meta["y"], seed=args.seed, t1_limit=args.t1_limit, t2_steps=args.t2_steps)
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(args.workers) as pool_proc:
        results = pool_proc.map(replay_pair, range(args.pairs), chunksize=8)
    t1 = {}
    for r in results:
        t1[r["t1"]] = t1.get(r["t1"], 0) + 1
    t1_mismatch = sum(r["t1_mismatch"] for r in results)
    t2 = {k: sum(r["t2"][k] for r in results) for k in results[0]["t2"]}
    t2_mismatch = sum(r["t2_mismatch"] for r in results)
    return {"pairs": args.pairs, "pool_states": len(pool), "T1": dict(categories=t1, mismatches=t1_mismatch,
                                                                       limit=args.t1_limit),
            "T2": dict(t2, mismatches=t2_mismatch), "seed": args.seed}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("rounds", type=int)
    ap.add_argument("--direct", type=int, default=2)
    ap.add_argument("--cap", type=int, default=1 << 66, help="trail cap (test builds only)")
    ap.add_argument("--replay", action="store_true")
    ap.add_argument("--pairs", type=int, default=2000)
    ap.add_argument("--t2-steps", type=int, default=4096)
    ap.add_argument("--t1-limit", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--pool-rounds", type=int, default=256)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    inp = Inputs(os.path.join(args.dir, "inputs.bin"))
    meta = json.load(open(os.path.join(args.dir, "instance.json")))
    y = meta["y"]
    if args.replay:
        rep = replay_main(inp, meta, args)
        out = dict(replay=rep, dir=args.dir)
        with open(os.path.join(args.dir, "replay.json"), "w") as f:
            json.dump(out, f, indent=1)
        print(out)
        return 0
    walk = Walk(inp, args.cap)
    states = walk.run(args.rounds)
    lines = dump_lines(states)
    out_path = os.path.join(args.dir, "replica-%d.txt" % args.rounds)
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    report = {"rounds": args.rounds, "walks": inp.NW, "S": inp.S, "RA": inp.RA,
              "DP_BITS": inp.DP_BITS, "cap": args.cap}

    dump_path = os.path.join(args.dir, "dump-%d.txt" % args.rounds)
    if os.path.exists(dump_path):
        kernel = open(dump_path).read().splitlines()
        k_walks = [l for l in kernel if not l.startswith("dp ")]
        k_dps = [l for l in kernel if l.startswith("dp ")]
        mismatch = [i for i, (a, b) in enumerate(zip(lines, k_walks)) if a != b]
        report["dump_lines"] = len(k_walks)
        report["dump_exact_match"] = (lines == k_walks)
        report["dump_mismatch_walks"] = mismatch[:10]
        report["dp_records"] = len(k_dps)
        # step.c dumps the first min(dp_head, DP_RING) slots of its ring; record j is stored in slot j % DP_RING
        ring = [None] * DP_RING
        for j, line in enumerate(walk.dp_lines):
            ring[j % DP_RING] = line
        report["dp_records_total"] = len(walk.dp_lines)
        report["dp_exact_match"] = (ring[:min(len(walk.dp_lines), DP_RING)] == k_dps)
    else:
        report["dump_exact_match"] = None

    native = os.path.join(args.dir, "native-%d.stderr" % args.rounds)
    if os.path.exists(native):
        kc = parse_counts(native)
        report["kernel_counts"] = kc
        report["replica_counts"] = dict(walk.counts)
        report["counts_match"] = all(kc.get(k) == walk.counts[k] for k in walk.counts)

    # Key check: formula vs direct exponentiation on a few states (initial and final).
    checks = []
    sample = [0, len(states) - 1][: args.direct]
    for w in sample:
        x = states[w]["x"] if w == len(states) - 1 else inp.walks[w][0]
        f = walk.key(x)
        d = direct_key(x)
        checks.append({"walk": w, "formula": [f[0] == d[0], f[2], d[2]], "value_equal": f[1] == d[1]})
    report["key_checks"] = checks
    report["key_formula_matches_direct"] = all(c["formula"][0] and c["value_equal"] for c in checks)

    # Exponent invariant: x = g^(a + y b). g = e(P, Q0) with P = [b]P0 from instance.json.
    g = pairing(G2, multiply(G1, meta["b"]))
    inv = []
    for s in states:
        inv.append(g ** ((s["a"] + y * s["b"]) % R) == s["x"])
    report["exponent_invariant_holds"] = all(inv)
    # Canonical exponents: sigma^kn(x) = g^(lambda^kn (a + y b)) for the sampled final states.
    canon = []
    for w in sample:
        st = states[w]
        xc = sigma_pow(st["x"], st["kn"])
        canon.append(xc == g ** (inp.LP[st["kn"]] * (st["a"] + y * st["b"]) % R))
    report["canonical_invariant_holds"] = all(canon)
    report["walks_checked"] = len(inv)
    report["restart_repeats"] = dict(walk.rep, restart_states=walk.counts["restarts"] + walk.counts["cap_restarts"])
    report["full_visit"] = dict(walk.fv, visited_keys=len(walk.visited), restart_outputs=len(walk.restart_keys))
    print(report)
    return 0 if report.get("dump_exact_match") in (True, None) else 1


if __name__ == "__main__":
    sys.exit(main())
