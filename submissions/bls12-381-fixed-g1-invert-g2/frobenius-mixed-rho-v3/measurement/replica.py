"""Python replica of step.c: reruns the blst GT walk and compares the kernel's dump exactly.

  python3 replica.py <input-dir> <rounds> [--direct N]

<input-dir> holds inputs.bin and instance.json from gen_inputs.py, plus the kernel's
dump-<rounds>.txt and native-<rounds>.stderr when present (measure.sh output). Writes
replica-<rounds>.txt next to them and prints a report.

Checks:
  1. the walk (same rules as step.c, including per-walk anchors and escapes) reproduces the
     kernel's dump line for line, and the event counts match the kernel's final stderr line;
  2. the 4-multiplication key (SPEC section 2) equals the direct minimum over the 12
     conjugates V_n(sigma^n(x)), where sigma^n(x) = x^(p^n) is computed by exponentiation
     (for --direct N states, default 2);
  3. exponent invariant: every final state equals g^(a + y b), with y the instance exponent.
"""
import argparse
import json
import os
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


class Inputs:
    def __init__(self, path):
        data = open(path, "rb").read()
        self.off = 0
        self.data = data
        assert self.u64() == MAGIC and self.u64() == 2, "bad magic or version"
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


class Walk:
    def __init__(self, inp):
        self.inp = inp
        self.RA = inp.RA
        self.S = inp.S
        self.dp_mask = inp.DP_MASK
        # steps: main walk steps. adds, sqrs: every group operation, including the extra
        # steps of cycle checks and the escapes (squarings), as counted by step.c.
        self.counts = dict(steps=0, adds=0, sqrs=0, dps=0, cyc2=0, cycl=0)
        self.dp_lines = []  # distinguished-point records in order: fingerprint, canonical exponents

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

    def escape(self, m):
        x2 = m["x"] * m["x"]
        mm, v, n = self.key(x2)
        self.counts["sqrs"] += 1
        return dict(x=x2, a=2 * m["a"] % R, b=2 * m["b"] % R, K=mm, Kc=v, kn=n)

    def dp_hit(self, K):
        return (K >> self.inp.ADDER_BITS) & self.dp_mask == 0

    def run(self, rounds):
        inp = self.inp
        NW = inp.NW
        states, prevK, prevValid = [], [None] * NW, [False] * NW
        anchorK, anchorValid = [None] * NW, [False] * NW
        for w, (x, a, b) in enumerate(inp.walks):
            m, v, n = self.key(x)
            states.append(dict(x=x, a=a, b=b, K=m, Kc=v, kn=n))
        for rnd in range(rounds):
            for w in range(NW):
                c = states[w]
                nxt = self.step(c)
                self.counts["steps"] += 1
                two = False
                if prevValid[w] and prevK[w] == nxt["K"]:
                    nn = self.step(nxt)
                    if nn["K"] == c["K"]:
                        two = True
                if two:
                    self.counts["cyc2"] += 1
                    nxt = self.escape(nxt if nxt["K"] < c["K"] else c)
                elif anchorValid[w] and anchorK[w] == nxt["K"]:
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
                        nxt = self.escape(mn)
                prevK[w] = c["K"]
                prevValid[w] = True
                if (rnd + 1) % ANCHOR_PERIOD == 0:
                    anchorK[w] = nxt["K"]
                    anchorValid[w] = True
                states[w] = nxt
                if self.dp_hit(nxt["K"]):
                    self.counts["dps"] += 1
                    m = nxt["K"]
                    lam = self.inp.LP[nxt["kn"]]
                    # fingerprint: Montgomery key bits 70..197 (bits 10..69 are zero for every distinguished point)
                    self.dp_lines.append("dp %016x %016x %064x %064x" % (
                        (m >> 70) & (2**64 - 1), (m >> 134) & (2**64 - 1), nxt["a"] * lam % R, nxt["b"] * lam % R))
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("rounds", type=int)
    ap.add_argument("--direct", type=int, default=2)
    args = ap.parse_args()
    inp = Inputs(os.path.join(args.dir, "inputs.bin"))
    meta = json.load(open(os.path.join(args.dir, "instance.json")))
    y = meta["y"]
    walk = Walk(inp)
    states = walk.run(args.rounds)
    lines = dump_lines(states)
    out_path = os.path.join(args.dir, "replica-%d.txt" % args.rounds)
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    report = {"rounds": args.rounds, "walks": inp.NW, "S": inp.S, "RA": inp.RA}

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
        report["dp_exact_match"] = (walk.dp_lines == k_dps)
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
        xc = st["x"] if st["kn"] == 0 else st["x"] ** (P ** st["kn"])
        canon.append(xc == g ** (inp.LP[st["kn"]] * (st["a"] + y * st["b"]) % R))
    report["canonical_invariant_holds"] = all(canon)
    report["walks_checked"] = len(inv)
    print(report)
    return 0 if report.get("dump_exact_match") in (True, None) else 1


if __name__ == "__main__":
    sys.exit(main())
