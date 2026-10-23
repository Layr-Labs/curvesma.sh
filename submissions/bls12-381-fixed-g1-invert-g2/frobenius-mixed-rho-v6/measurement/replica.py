"""Python replica of step.c (walk-spec v11). Reruns the walk and compares the kernel dump exactly.

  python3 replica.py <input-dir> <rounds> [--cap C] [--escape-zero]
  python3 replica.py <input-dir> 0 --replay [--pairs N] [--t2-steps K] [--workers W] [--pool-rounds R]
"""
import argparse
import hashlib
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
VERSION = 6
KB = 0x5851F42D4C957F2D  # branch key kappa_b (spec section 5, measurement value; step.c BRANCH_KEY)
PRF_KEY_TEXT = b"gt-v11-measurement"  # step.c PRF_KEY_TEXT
KP = hashlib.sha256(PRF_KEY_TEXT).digest()  # PRF key (spec section 6)
COMB_BYTES = 32
COMB_VALUES = 255
M64 = (1 << 64) - 1
DP_RING = 65536


def mix64(z):
    """mix of spec section 5 (splitmix64 finaliser), all arithmetic mod 2^64."""
    z &= M64
    z ^= z >> 30
    z = (z * 0xBF58476D1CE4E5B9) & M64
    z ^= z >> 27
    z = (z * 0x94D049BB133111EB) & M64
    z ^= z >> 31
    return z


def prf(m):
    """PRF(m) = (SHA256(KP || m || 0x00) || SHA256(KP || m || 0x01)) as a."""
    h = hashlib.sha256(KP + m + b"\x00").digest() + hashlib.sha256(KP + m + b"\x01").digest()
    return int.from_bytes(h, "big") % R


def restart_msg(w, t):
    return b"\x01" + w.to_bytes(8, "big") + t.to_bytes(8, "big")


def escape_msg(K):
    """0x02 || K_s: the 48-byte big-endian Montgomery integer of the key."""
    return b"\x02" + K.to_bytes(48, "big")


class Inputs:
    def __init__(self, path):
        data = open(path, "rb").read()
        self.off = 0
        self.data = data
        assert self.u64() == MAGIC and self.u64() == VERSION, "bad magic or version"
        self.RA = self.u64()
        self.ADDER_BITS = self.u64()
        self.DP_BITS = self.u64()
        self.S = self.u64()
        self.NW = self.u64()
        assert self.u256() == R, "group order mismatch"
        self.LP = [self.u256() for _ in range(12)]
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
        self.G = self.fp12()
        assert self.off == len(data), "trailing bytes in the input"
        self.comb = build_comb(self.G)

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


def build_comb(g):
    """COMB[i][v-1] = g^(v 2^(8i)), i < 32, v = 1..255 (spec section 6): the."""
    table = []
    base = g
    for i in range(COMB_BYTES):
        row = [base]
        for _ in range(2, COMB_VALUES + 1):
            row.append(row[-1] * base)
        table.append(row)
        if i + 1 < COMB_BYTES:
            for _ in range(8):
                base = base * base
    return table


def comb_mul(comb, x, e):
    """x * g^e by the comb: one multiplication per nonzero byte of e."""
    n = 0
    for i in range(COMB_BYTES):
        b = (e >> (8 * i)) & 0xFF
        if b:
            x = x * comb[i][b - 1]
            n += 1
    return x, n


_FROB = None


def frob_matrix():
    """Matrix of the p-power map on the py_ecc coefficients: sigma(c) = M c."""
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
    """py_ecc FQ12 -> the 12 canonical Fp of blst struct order, as 96-char."""
    c = x.coeffs
    out = []
    A = {i: ((c[i] + c[i + 6]) % P, c[i + 6] % P) for i in range(6)}
    for i in STRUCT_ORDER:
        re, im = A[i]
        out.append("%096x" % re)
        out.append("%096x" % im)
    return out


def fresh(s, phase, rnd, restarted, age=0):
    """The state after a restart or an escape: detector cleared, trail age."""
    st = dict(s, prevK=None, prevValid=False, anchorK=None, anchorValid=False, age=age)
    if not restarted and (rnd - phase + 1) % ANCHOR_PERIOD == 0:
        st["anchorK"] = s["K"]
        st["anchorValid"] = True
    return st


class Walk:
    def __init__(self, inp, cap=1 << 66, track=True, escape_zero=False):
        self.inp = inp
        self.RA = inp.RA
        self.S = inp.S
        self.dp_mask = inp.DP_MASK
        self.cap = cap
        self.track = track
        self.escape_zero = escape_zero
        self.counts = dict(steps=0, adds=0, sqrs=0, dps=0, cyc2=0, cycl=0, restarts=0, cap_restarts=0, chained=0,
                           escapes=0, escape_dps=0, comb_mults=0)
        self.diag = dict(cap_after_escape=0, escapes_in_capped_trails=0)
        self.trail_escapes = [0] * inp.NW
        self.dp_lines = []
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
        """4-multiplication formula (SPEC section 2)."""
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
        """s = dict(x, a, b, K, Kc, kn)."""
        inp = self.inp
        u = mix64((s["K"] & M64) ^ KB)
        h = u & (self.RA - 1)
        if (u >> (64 - inp.ADDER_BITS)) < self.S:
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
        """Escape from the least-key state s (spec section 6, v11): x <-."""
        inp = self.inp
        self.last_least = s
        cs = sigma_pow(s["x"], s["kn"])
        e = 0 if self.escape_zero else prf(escape_msg(s["K"]))
        x, n = comb_mul(inp.comb, cs, e)
        self.counts["comb_mults"] += n
        lam = inp.LP[s["kn"]]
        A = s["a"] * lam % R
        B = s["b"] * lam % R
        m, v, n = self.key(x)
        self.counts["escapes"] += 1
        return dict(x=x, a=(A + e) % R, b=B, K=m, Kc=v, kn=n)

    def dp_hit(self, K):
        return (K >> self.inp.ADDER_BITS) & self.dp_mask == 0

    def restart(self, s, w):
        """Restart of walk w (spec section 6, v11): t = the walk's counter after."""
        inp = self.inp
        self.rc[w] += 1
        t = self.rc[w]
        e = prf(restart_msg(w, t))
        x, n = comb_mul(inp.comb, s["x"], e)
        self.counts["comb_mults"] += n
        m, v, k = self.key(x)
        out = dict(x=x, a=(s["a"] + e) % R, b=s["b"], K=m, Kc=v, kn=k)
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
        """A distinguished point s: store its record, restart; chained while the."""
        cur = s
        while True:
            self.counts["dps"] += 1
            m = cur["K"]
            lam = self.inp.LP[cur["kn"]]
            self.dp_lines.append("dp %016x %016x %064x %064x" % (
                (m >> 70) & (2**64 - 1), (m >> 134) & (2**64 - 1), cur["a"] * lam % R, cur["b"] * lam % R))
            cur = self.restart(cur, w)
            self.counts["restarts"] += 1
            if not self.dp_hit(cur["K"]):
                break
            self.counts["chained"] += 1
        return cur

    def advance(self, c, w, rnd, phase=0, allow_dp=True):
        """One round of walk w (spec sections 5-6, v11)."""
        self.last_event = None
        nxt = self.step(c)
        self.counts["steps"] += 1
        if self.track:
            self.visited.add(nxt["K"])
        if allow_dp and self.dp_hit(nxt["K"]):
            self.last_event = "dp"
            self.last_dp = nxt
            self.trail_escapes[w] = 0
            return fresh(self.dp_chain(nxt, w), phase, rnd, True)
        age = c["age"] + 1
        if age >= self.cap:
            self.counts["cap_restarts"] += 1
            if self.trail_escapes[w]:
                self.diag["cap_after_escape"] += 1
                self.diag["escapes_in_capped_trails"] += self.trail_escapes[w]
            self.trail_escapes[w] = 0
            out = self.restart(nxt, w)
            if allow_dp and self.dp_hit(out["K"]):
                out = self.dp_chain(out, w)
            return fresh(out, phase, rnd, True)
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
            self.trail_escapes[w] += 1
            if allow_dp and self.dp_hit(ev["K"]):
                self.counts["escape_dps"] += 1
                self.trail_escapes[w] = 0
                return fresh(self.dp_chain(ev, w), phase, rnd, True)
            return fresh(ev, phase, rnd, False, age)
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
    """Direct min over V_n(sigma^n x), sigma^n x = x^(p^n), by."""
    best = None
    for n in range(12):
        y = x if n == 0 else x ** (P ** n)
        v = (y.coeffs[1] + y.coeffs[7]) % P
        m = v * MONT % P
        if best is None or m < best[0]:
            best = (m, v, n)
    return best


REPLAY = {}


def _elem_eq(x1, x2):
    return tuple(x1.coeffs) == tuple(x2.coeffs)


def _pair_states(walk, pool, rng, y):
    """Two walk states at the same element: pool element, different."""
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
    """T1 (first-DP equality, DP test on) and T2 (detection and escape."""
    ctx = REPLAY
    walk, pool, y = ctx["walk"], ctx["pool"], ctx["y"]
    rng = random.Random(ctx["seed"] * 1000003 + i)
    A0, B0 = _pair_states(walk, pool, rng, y)
    res = dict(i=i)
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
    A, B = dict(A0), dict(B0)
    phi = {"A": A0["phi"], "B": B0["phi"]}
    ev_lists = {"A": [], "B": []}
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
    ap.add_argument("--escape-zero", action="store_true", help="test builds only: escape exponent 0")
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
    walk = Walk(inp, args.cap, escape_zero=args.escape_zero)
    states = walk.run(args.rounds)
    lines = dump_lines(states)
    out_path = os.path.join(args.dir, "replica-%d.txt" % args.rounds)
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    report = {"rounds": args.rounds, "walks": inp.NW, "S": inp.S, "RA": inp.RA,
              "DP_BITS": inp.DP_BITS, "cap": args.cap, "escape_zero": args.escape_zero}

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

    checks = []
    sample = [0, len(states) - 1][: args.direct]
    for w in sample:
        x = states[w]["x"] if w == len(states) - 1 else inp.walks[w][0]
        f = walk.key(x)
        d = direct_key(x)
        checks.append({"walk": w, "formula": [f[0] == d[0], f[2], d[2]], "value_equal": f[1] == d[1]})
    report["key_checks"] = checks
    report["key_formula_matches_direct"] = all(c["formula"][0] and c["value_equal"] for c in checks)

    g = pairing(G2, multiply(G1, meta["b"]))
    report["generator_is_input"] = (g == inp.G)
    inv = []
    for s in states:
        inv.append(g ** ((s["a"] + y * s["b"]) % R) == s["x"])
    report["exponent_invariant_holds"] = all(inv)
    canon = []
    for w in sample:
        st = states[w]
        xc = sigma_pow(st["x"], st["kn"])
        canon.append(xc == g ** (inp.LP[st["kn"]] * (st["a"] + y * st["b"]) % R))
    report["canonical_invariant_holds"] = all(canon)
    report["walks_checked"] = len(inv)
    report["restart_repeats"] = dict(walk.rep, restart_states=walk.counts["restarts"] + walk.counts["cap_restarts"])
    report["full_visit"] = dict(walk.fv, visited_keys=len(walk.visited), restart_outputs=len(walk.restart_keys))
    report["diag"] = walk.diag
    print(report)
    return 0 if report.get("dump_exact_match") in (True, None) else 1


if __name__ == "__main__":
    sys.exit(main())
