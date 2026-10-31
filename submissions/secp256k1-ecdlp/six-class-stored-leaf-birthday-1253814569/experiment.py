"""Stored-leaf birthday search on secp256k1 (proof.md). Part A runs kernels/leaf.s and secp_fe.s on production points
against exact arithmetic and the exact class rule; part B runs the emits, both hops (every message delivered by
route.s dlv) and the grouping on a 4 x 4 machine through kernels/leaf.s and route.s; part C runs the toy pipeline at
q25; part D runs certificate.py. Input JSON on stdin, output JSON."""
import contextlib
import io
import json
import math
import os
import random
import sys
import time

T0 = time.time()
cfg = json.load(sys.stdin)
SEED = cfg.get("seed", 1)
TRIALS = cfg.get("trials", 6000)
CURVE = cfg.get("curve_trials", 150)
BUDGET = cfg.get("seconds", 15)
DO_KERNELS = cfg.get("kernels", True)
RESERVE = 0.5                                   # seconds kept for the output
T_END = T0 + BUDGET - RESERVE
MIN_RANDOM, MIN_NEAR = 2000, 600                # class-rule checks that always run
MIN_MULTI, MIN_CURVE = 400, 10                 # toy trials that always run (each at most its requested count)


def package_root():
    """The package directory: this file's, files/ beside it, or /snapshot/files (the hosted runner's layout).

    returns: absolute path of the directory that contains kernels/a64.py.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    for base in (here, os.path.join(here, "files"), "/snapshot/files"):
        if os.path.isfile(os.path.join(base, "kernels", "a64.py")):
            return base
    raise FileNotFoundError("kernels/a64.py is not beside " + here + " nor under files/ or /snapshot/files")


HERE = package_root()
sys.path.insert(0, os.path.join(HERE, "kernels"))
import a64  # noqa: E402

M64 = 2 ** 64 - 1
p = 2 ** 256 - 2 ** 32 - 977
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
GX = int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16)
GY = int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16)
BETA = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)
LAM = int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72", 16)
C0 = 2 ** 32 + 977
GQ = (BETA << 384) // p
GL = [(GQ >> (64 * k)) & M64 for k in range(6)]
OUT = {}


class Curve:
    """Short Weierstrass y^2 = x^3 + 7 over F_q with prime group order r."""

    def __init__(self, q, order):
        """q: field prime; order: group order. returns: None."""
        self.q, self.r = q, order

    def add(self, P, Q):
        """Affine addition. P, Q: points or None (infinity). returns: P + Q."""
        q = self.q
        if P is None:
            return Q
        if Q is None:
            return P
        if P[0] == Q[0]:
            if (P[1] + Q[1]) % q == 0:
                return None
            l = 3 * P[0] * P[0] * pow(2 * P[1], -1, q) % q
        else:
            l = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, q) % q
        x = (l * l - P[0] - Q[0]) % q
        return (x, (l * (P[0] - x) - P[1]) % q)

    def neg(self, P):
        """Negation. P: point. returns: -P."""
        return None if P is None else (P[0], (-P[1]) % self.q)

    def mul(self, k, P):
        """Double-and-add. k: integer; P: point. returns: [k]P."""
        k %= self.r
        acc = None
        while k:
            if k & 1:
                acc = self.add(acc, P)
            P = self.add(P, P)
            k >>= 1
        return acc


PROD = Curve(p, r)
GP = (GX, GY)


def limbs(v, n=4):
    """Little-endian 64-bit words. v: integer; n: count. returns: list."""
    return [(v >> (64 * k)) & M64 for k in range(n)]


def put(mem, a, v, n=4):
    """Store v as n words at address a. mem: dict; a: address; v: integer; n: words. returns: None."""
    for k, w in enumerate(limbs(v, n)):
        mem[a + 8 * k] = w


def get(mem, a, n=4):
    """Read n words at address a. mem: dict; a: address; n: words. returns: integer."""
    return sum(mem[a + 8 * k] << (64 * k) for k in range(n))


def exact_key(x):
    """Canonical class value and index. x: field element. returns: (min(x, beta x, beta^2 x), j)."""
    v = [x, BETA * x % p, BETA * BETA % p * x % p]
    j = min(range(3), key=lambda t: v[t])
    return v[j], j


def rule_j(x):
    """Bit-exact model of the KEY macro (proof.md 3.3). x: integer < 2^256. returns: j, or None when a window fires."""
    xl = limbs(x)
    f1 = sum(((xl[i] * GL[5 - i]) & M64) + ((xl[i] * GL[4 - i]) >> 64) for i in range(4)) & M64
    f0 = xl[3]
    s, d = (f0 + f1) & M64, (f0 - f1) & M64
    m = f1 if f0 >= f1 else f0
    e = m + s
    if f1 + 16 > M64 or s + 17 > M64 or (e & M64) + 33 > M64 or d < 16:
        return None
    return 2 if e > M64 else 1 if f0 >= f1 else 0


def deadline(frac):
    """Deadline of a budget-limited loop: a share of the budget that remains when the loop starts.

    frac: share of the remaining budget. returns: absolute time.
    """
    now = time.time()
    return now + frac * max(0.0, T_END - now)


def near_tie_x(rg):
    """A near-tie input for the class rule: a small value divided by one of the six linear forms.

    rg: Random. returns: field element.
    """
    mlt = rg.choice([BETA - 1, BETA * BETA - 1, BETA - BETA * BETA, 1 + BETA, 2 + BETA, 1 + 2 * BETA])
    dd = rg.randrange(-(1 << rg.randrange(150, 250)), 1 << rg.randrange(150, 250))
    return dd * pow(mlt % p, -1, p) % p or 1


def part_a_identities(rng):
    """Production identities and the class rule against the exact rule. rng: Random. returns: dict."""
    res = {}
    res["beta_cube_root_of_1_mod_p"] = pow(BETA, 3, p) == 1 and BETA != 1
    res["lambda_cube_root_of_1_mod_r"] = pow(LAM, 3, r) == 1 and LAM != 1
    res["lambda_G_equals_beta_x"] = PROD.mul(LAM, GP) == (BETA * GX % p, GY)
    res["seven_nonresidue_no_x0_point"] = pow(7, (p - 1) // 2, p) == p - 1
    res["window_forms_nonzero"] = all(c % p for c in (1 - BETA, BETA, 1 + BETA, 2 + BETA, 1 + 2 * BETA))
    bad = drops = done = 0
    kinds = {"random": 0, "near_tie": 0}
    dl = deadline(0.10)
    while done < 26000 and (done < MIN_RANDOM or time.time() < dl):
        if done < 20000:
            x, kind = rng.randrange(1, p), "random"
        else:
            x, kind = near_tie_x(rng), "near_tie"
        kinds[kind] += 1
        done += 1
        j = rule_j(x)
        if j is None:
            drops += 1
        elif j != exact_key(x)[1]:
            bad += 1
    extra = max(0, MIN_NEAR - kinds["near_tie"])
    rng2 = random.Random(SEED + 7919)  # the guaranteed near ties come from their own stream, so rng is unchanged
    for _ in range(extra):
        x = near_tie_x(rng2)
        kinds["near_tie"] += 1
        j = rule_j(x)
        if j is None:
            drops += 1
        elif j != exact_key(x)[1]:
            bad += 1
    res["class_rule"] = dict(checked=kinds, windows=drops, mismatches=bad, near_ties_from_second_stream=extra)
    return res


SENTINEL = 0x5E5E5E5E5E5E5E5E
M36 = 2 ** 36 - 1


class Port:
    """Capture of a 48-byte port (proof.md 2.1): destination +0, address +8, four data words +16..+47; a message is
    complete when its last word, at +40, is written."""

    def __init__(self):
        """Empty capture. returns: None."""
        self.words, self.msgs = [0] * 6, []

    def store(self, k, v):
        """Device hook. k: byte offset in the block; v: word stored. returns: None."""
        self.words[k // 8] = v
        if k == 40:
            self.msgs.append(list(self.words))


LINK = 1 << 40   # where dlv finds the message in transit: an address of the link, not of a node


class LinkView:
    """Address space of one delivery: the message's 48 bytes at LINK, the destination's memory everywhere else."""

    def __init__(self, mem, msg):
        """mem: the destination's memory; msg: [destination, address, w0, w1, w2, w3]. returns: None."""
        self.mem, self.msg = mem, msg

    def __getitem__(self, a):
        """A load. a: byte address. returns: the word."""
        return self.msg[(a - LINK) // 8] if LINK <= a < LINK + 48 else self.mem[a]

    def __setitem__(self, a, v):
        """A store, into the destination's memory. a: byte address; v: word. returns: None."""
        assert not LINK <= a < LINK + 48
        self.mem[a] = v


def link(prog, labels, mem, msg, calls, cert):
    """Deliver one message (proof.md 2.1) by running route.s dlv on it, so its four words reach the destination through
    counted stores; record the count against the certificate's DLV.

    prog, labels: route.s; mem: the destination's memory; msg: [destination, address, w0, w1, w2, w3]; calls: list
    receiving (entry, executed count, formula). returns: None."""
    calls.append(("dlv", a64.a64_run(prog, labels, "dlv", dict(x0=LINK, x30=SENTINEL), LinkView(mem, msg)), cert.DLV))


class LazyTable(dict):
    """Memory in which the unwritten words of one cursor table read as its initial entry (cursor, slot end), standing in
    for tinit's 2^36 entries in the production-size runs of part A."""

    def __init__(self, lo, hi, base, end):
        """lo, hi: the table's address range; base, end: the initial entry. returns: None."""
        super().__init__()
        self.lo, self.hi, self.base, self.end = lo, hi, base, end

    def __missing__(self, a):
        """An unwritten word. a: address. returns: base or end inside the table (KeyError elsewhere)."""
        if self.lo <= a < self.hi:
            return self.base if (a - self.lo) % 16 == 0 else self.end
        raise KeyError(a)


def header_values(X):
    """The eight header values of X (proof.md 3.1). X: point. returns: list of field elements, in header order."""
    x, y = X
    b2, i2, i4 = BETA * BETA % p, pow(2, -1, p), pow(4, -1, p)
    return [x, y, BETA * x % p, b2 * x % p, x * i2 % p, BETA * x * i2 % p, 3 * x * x * i4 % p, 3 * b2 * x * x * i4 % p]


def y_entry(Y):
    """The Y table entry (proof.md 3.1): y^2, y, x, 3x^2 and x + C0. Y: point. returns: list of five integers."""
    x, y = Y
    return [y * y % p, y, x, 3 * x * x % p, x + C0]


def six_leaves(X, Y):
    """Exact leaves x(phi^t(X) + Y) and x(phi^t(X) - Y), in entry order c = 2t + s. X, Y: points.
    returns: list of six field elements."""
    out = []
    for t in range(3):
        PX = (pow(BETA, t, p) * X[0] % p, X[1])
        out += [PROD.add(PX, Y)[0], PROD.add(PX, PROD.neg(Y))[0]]
    return out


def header(mem, Hb, X=None):
    """Header of the batch kernels: G1..G5, the rows of beta and beta^2, and X's eight values when X is given.

    mem: memory; Hb: header address; X: point or None. returns: None.
    """
    if X is not None:
        for k, v in enumerate(header_values(X)):
            put(mem, Hb + 32 * k, v)
    for k in range(1, 6):
        mem[Hb + 408 + 8 * k] = GL[k]
    for jj, base in ((1, Hb + 456), (2, Hb + 584)):
        b = pow(BETA, jj, p)
        for k in range(4):
            put(mem, base + 32 * k, b * pow(2, 64 * k, p) % p)


def batch_tables(mem, Hb, Yb, X, Ys):
    """Write X's header and the Y table. mem: memory; Hb, Yb: addresses; X: point; Ys: points. returns: None."""
    header(mem, Hb, X)
    for i, Y in enumerate(Ys):
        for k, v in enumerate(y_entry(Y)):
            put(mem, Yb + 160 * i + 32 * k, v)


def emit_reference(leaves, W2, start, end, cert):
    """Exact model of keys' emits (proof.md 3.2) in keys' order: entries from the last down, leaves c = 0..5 (type c + 1).

    leaves: six leaves per incidence; W2: index word; start, end: initial cursor and slot end of d1; cert: certificate.
    returns: (messages, final cursors {d1: cursor}, products, instructions saved by drops, window drops, slot drops).
    """
    msgs, cur, prods, saved, wins, full = [], {}, 0, 0, 0, 0
    for l in range(len(leaves) - 1, -1, -1):
        for c, x in enumerate(leaves[l]):
            j = rule_j(x)
            if j is None:
                wins += 1
                saved += 1 + cert.EMIT
                continue
            K, jj = exact_key(x)
            assert jj == j
            prods += j != 0
            kl = limbs(K)
            d1 = kl[0] & M36
            cc = cur.get(d1, start(d1))
            if cc >= end(d1):
                full += 1
                saved += cert.EMIT - cert.EMIT_CHECK + (j != 0)
                continue
            msgs.append([d1, cc, kl[2], kl[3], W2, (kl[1] & M36) + ((c + 1) << 36) + (l << 39)])
            cur[d1] = cc + 32
    return msgs, cur, prods, saved, wins, full


def leaf_batch(prog, labels, X, Ys, cert, W2=0x123456789, entries_override=None):
    """Run fwd, bwd and keys of kernels/leaf.s on one batch; each call must return through the sentinel x30.

    prog, labels: loaded kernel; X: point; Ys: points; cert: certificate module (count constants); W2: index word;
    entries_override: optional six leaves per incidence written into the entries before the key pass (window tests).
    returns: dict with counts, checks and the messages.
    """
    q2 = len(Ys)
    Hb, Yb, Sb, Pb, SP = 0x100000, 0x200000, 1 << 34, 0x40000000, 0x7000000
    BASE, END = 1 << 62, (1 << 62) + (1 << 40)
    mem = LazyTable(Pb + 64, Pb + 64 + 16 * 2 ** 36, BASE, END)
    batch_tables(mem, Hb, Yb, X, Ys)
    regs = dict(x0=q2, x1=Yb, x3=Sb, x25=1, x26=0, x27=0, x28=0, sp=SP, x30=SENTINEL)
    regs.update(as_regs(X[1] * X[1] % p, ("x19", "x20", "x21", "x22")))
    out = {}
    nf = a64.a64_run(prog, labels, "fwd", regs, mem, out_regs=out)
    acc = reg_int(out, ("x25", "x26", "x27", "x28"))
    want = 1
    for Y in Ys:
        want = want * (Y[1] * Y[1] - X[1] * X[1]) % p
    R = pow(acc, p - 2, p)
    regs = dict(x1=Yb + 160 * q2 + 32, x2=Hb, x3=Sb + 256 * (q2 - 1), sp=SP, x30=SENTINEL)
    regs.update(as_regs(R, ("x25", "x26", "x27", "x28")))
    out = {}
    nb = a64.a64_run(prog, labels, "bwd", regs, mem, out_regs=out)
    final_R = reg_int(out, ("x25", "x26", "x27", "x28"))   # the last R update, R_0 D_0 = 1, left unread (3.4)
    slot_bad = 0
    leaves = []
    for i, Y in enumerate(Ys):
        six = six_leaves(X, Y)
        slot_bad += sum(get(mem, Sb + 256 * i + 32 * c) != six[c] for c in range(6))
        leaves.append(six)
    if entries_override:
        leaves = entries_override
        for i, six in enumerate(leaves):
            for c, x in enumerate(six):
                put(mem, Sb + 256 * i + 32 * c, x)
    port = Port()
    regs = dict(x0=Pb, x2=Hb, x3=Sb + 256 * (q2 - 1), x17=W2, x18=0, sp=SP, x30=SENTINEL)
    out = {}
    nk = a64.a64_run(prog, labels, "keys", regs, mem, out_regs=out, dev=(Pb, 48, port.store))
    exp, cur, prods, saved, wins, _ = emit_reference(leaves, W2, lambda d: BASE, lambda d: END, cert)
    cursors_ok = all(mem[Pb + 64 + 16 * d] == c for d, c in cur.items())
    return dict(fwd=nf, bwd=nb, keys=nk, acc_ok=acc % p == want, leaf_mismatches=slot_bad,
                messages_ok=port.msgs == exp and cursors_ok, messages=len(port.msgs), products=out["x18"],
                products_ref=prods, windows=wins, port_fixed=out["x0"] == Pb,
                keys_formula=nk == cert.KEY_FIX + cert.KEY_INC * q2 + cert.KMUL_PATH * out["x18"] - saved,
                final_R_is_one=final_R % p == 1)


def window_x(rng, form):
    """A field element whose class rule hits a window. rng: Random; form: multiplier of the small linear form.
    returns: x with rule_j(x) None."""
    while True:
        x = rng.randrange(1, 2 ** 150) * pow(form % p, -1, p) % p
        if rule_j(x) is None:
            return x


SITE_KERNELS = """
t_shx:
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tSHX x5,x25,x28,x26,x27
\tret
t_fdx:
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tFDX x10,x11,x7,x24,x10,x11,x5,x6,x12,x13
\tret
t_squ:
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tstr x30,[sp,#-16]!
\tSQ5 x14,x15,x16,x17,x5,x6,x7,x8,x9,x0,x18,x29,x30,x19,x20,x21,x22,x23,x24
\tF5 x8,x5,x6,x9,x0,x19,x20,x21,x22,x18,x29
\tldr x30,[sp],#16
\tret
t_sqv:
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tstr x30,[sp,#-16]!
\tSQ5 x10,x11,x12,x13,x5,x6,x7,x8,x9,x14,x15,x16,x17,x23,x24,x0,x18,x29,x30
\tF5 x8,x5,x6,x9,x14,x10,x11,x12,x13,x15,x16
\tldr x30,[sp],#16
\tret
t_sqc:
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tstr x30,[sp,#-16]!
\tSQ5 x10,x11,x12,x13,x5,x6,x7,x8,x9,x14,x15,x16,x17,x23,x24,x0,x18,x29,x30
\tldp x15,x16,[x2,#192]
\tldp x17,x18,[x2,#208]
\tAC4 x8,x5,x6,x9,x14,x15,x16,x17,x18
\tF5 x8,x5,x6,x9,x14,x10,x11,x12,x13,x15,x16
\tldr x30,[sp],#16
\tret
t_ltf:
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tSB6 x14,x15,x16,x17,x5,x29,x6,x7,x10,x11
\tFDX x14,x15,x5,x29,x14,x15,x18,x19,x16,x17
\tret
t_kmul:
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tstr x30,[sp,#-16]!
\tKMUL
\tldr x30,[sp],#16
\tret
t_fs:
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tFS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
\tret
t_faddw:
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tFADDW x23,x24,x25,x26,x5,x6,x7,x8,x14,x15,x16,x17,x9
\tret
"""


def run_site(prog, labels, entry, regs, mem=None):
    """Run one site-test entry with a sentinel x30.

    prog, labels: loaded program; entry: label; regs: input registers; mem: memory dict or None.
    returns: the final registers.
    """
    out = {}
    a64.a64_run(prog, labels, entry, dict(regs, sp=0x7000000, x30=SENTINEL), {} if mem is None else mem, out_regs=out)
    return out


def reg_int(regs, names):
    """Integer from registers, low limb first. regs: register values; names: register names. returns: integer."""
    return sum(regs[n] << (64 * k) for k, n in enumerate(names))


def as_regs(v, names):
    """Registers holding the limbs of v, low limb first. v: integer < 2^256; names: register names. returns: dict."""
    return {n: w for n, w in zip(names, limbs(v, len(names)))}


def macro_counts(src, entry, names):
    """Invocations of the named macros, through nested macros, in one kernel of leaf.s (its label to its first ret).

    src: leaf.s text; entry: kernel label; names: macro names. returns: {name: count}.
    """
    macros, cur, lines = {}, None, [ln.split("//")[0].strip() for ln in src.split("\n")]
    for s in lines:
        if s.startswith(".macro"):
            cur = s.split()[1].rstrip(",")
            macros[cur] = []
        elif s.startswith(".endm"):
            cur = None
        elif cur is not None:
            macros[cur].append(s)
    out = dict.fromkeys(names, 0)

    def walk(seq):
        """Count invocations in a list of lines. seq: lines. returns: None."""
        for s in seq:
            w = s.split()[0] if s else ""
            if w in out:
                out[w] += 1
            if w in macros:
                walk(macros[w])
    start = lines.index(entry + ":")
    walk(lines[start + 1:lines.index("ret", start)])
    return out


def sqrt_mod_p(v):
    """A square root mod p, or None. v: residue. returns: integer or None."""
    s_ = pow(v, (p + 1) // 4, p)
    return s_ if s_ * s_ % p == v % p else None


def site_checks(rng, fprog, flabels):
    """Each site type of proof.md 3.4 and F's field routines against exact arithmetic: inputs built to meet a site's
    predicate must deviate, boundary and random ones must not, and each must match the bit-level model; FS, FADDW,
    fe_mul_f and fe_sqr_f must be exact. rng: Random; fprog, flabels: loaded secp_fe.s. returns: dict.
    """
    prog, labels = a64.a64_load(open(os.path.join(HERE, "kernels", "leaf.s")).read() + SITE_KERNELS)
    res = {}

    def tally(name, cases):
        """Record (predicate, deviates, matches model) triples. name: site; cases: list. returns: None."""
        res[name] = dict(cases=len(cases), predicate_holds=sum(c[0] for c in cases),
                         deviates_iff_predicate=all(c[0] == c[1] for c in cases), model=all(c[2] for c in cases))

    def fold(v):
        """One fold of a 5- or 6-word value: (v mod 2^256) + floor(v / 2^256) C0, unreduced. returns: integer."""
        return (v % 2 ** 256) + (v >> 256) * C0

    cases = []
    for kind in (0, 1, 2):          # built to carry, just below the carry, random
        for _ in range(48):
            b3 = rng.randrange(2 ** 40, 2 ** 64)
            hi = (b3 * C0) >> 64
            low = (2 ** 192 - hi + rng.randrange(hi), 2 ** 192 - hi - 1 - rng.randrange(2 ** 20),
                   rng.getrandbits(192))[kind]
            B = low + (b3 << 192)
            o = reg_int(run_site(prog, labels, "t_shx", as_regs(B, ("x25", "x26", "x27", "x28"))),
                        ("x5", "x25", "x26", "x27"))
            want = B * 2 ** 64 - b3 * p
            cases.append((low + hi >= 2 ** 192, o != want, o == want % 2 ** 256))
    tally("row_shift", cases)
    cases = []
    for kind in (0, 1, 2):
        for _ in range(48):
            u = rng.randrange(1, 2 ** 66)
            m = u * C0
            low = (2 ** 256 - m + rng.randrange(m), 2 ** 256 - m - 1 - rng.randrange(2 ** 20), rng.getrandbits(256))[kind]
            Z = low + (u << 256)
            regs = as_regs(low, ("x10", "x11", "x12", "x13"))
            regs.update(x7=u & M64, x24=u >> 64)
            o = reg_int(run_site(prog, labels, "t_fdx", regs), ("x10", "x11", "x12", "x13"))
            t = fold(Z)
            cases.append((t >= 2 ** 256, o != t, o == t % 2 ** 256))
    tally("fold", cases)
    cases = []
    for kind in (0, 1, 2):          # the fused fold: built to carry, just below the carry, random
        for _ in range(48):
            u = rng.randrange(1, 2 ** 66 + 2)
            m = u * C0
            low = (2 ** 256 - m + rng.randrange(m), 2 ** 256 - m - 1 - rng.randrange(2 ** 20), rng.getrandbits(256))[kind]
            sub = rng.getrandbits(256)
            Z = low + (u << 256) + sub - 2 ** 256      # the accumulator before SB6, which subtracts sub and adds 2^256
            regs = as_regs(Z, ("x14", "x15", "x16", "x17", "x5", "x29"))
            regs.update(as_regs(sub, ("x6", "x7", "x10", "x11")))
            o = reg_int(run_site(prog, labels, "t_ltf", regs), ("x14", "x15", "x16", "x17"))
            t = fold(low + (u << 256))
            cases.append((t >= 2 ** 256, o != t, o == t % 2 ** 256))
    tally("fused_fold", cases)
    Hb = 0x100000
    for name, entry, a_regs, o_regs in (("square_fold_u", "t_squ", ("x14", "x15", "x16", "x17"), ("x19", "x20", "x21", "x22")),
                                        ("square_fold_v", "t_sqv", ("x10", "x11", "x12", "x13"), ("x10", "x11", "x12", "x13")),
                                        ("square_plus_constant", "t_sqc", ("x10", "x11", "x12", "x13"),
                                         ("x10", "x11", "x12", "x13"))):
        cases = []
        for kind in (0, 1, 2):      # true value below 2^64, in (2^65, 2^80), random
            for _ in range(32):
                while True:
                    T = (rng.randrange(2 ** 64), rng.randrange(2 ** 65, 2 ** 80), rng.randrange(p))[kind]
                    if entry != "t_sqc":
                        a = sqrt_mod_p(T)
                        if a is None:
                            continue
                        a += rng.choice((0, p)) if a < 2 ** 256 - p else 0
                        c = 0
                    else:
                        a = rng.getrandbits(256)
                        c = (T - a * a) % p
                    break
                mem = {}
                put(mem, Hb + 192, c)
                regs = as_regs(a, a_regs)
                regs.update(x2=Hb)
                o = reg_int(run_site(prog, labels, entry, regs, mem), o_regs)
                t = fold(fold(a * a) + c)
                cases.append((t >= p, o != T, o == t % 2 ** 256))
        tally(name, cases)
    cases = []
    for kind in (0, 1, 2):          # FS(a, b) with b a fold output: b - a > p, b >= p without it, b random below p
        for _ in range(32):
            b = (rng.randrange(p, 2 ** 256), rng.randrange(p, 2 ** 256), rng.randrange(p))[kind]
            a = (rng.randrange(b - p) if b > p else 0, rng.randrange(b - p, p), rng.randrange(p))[kind]
            regs = as_regs(a, ("x5", "x6", "x7", "x8"))
            regs.update(as_regs(b, ("x19", "x20", "x21", "x22")))
            o = reg_int(run_site(prog, labels, "t_fs", regs), ("x9", "x10", "x11", "x12"))
            d = (a - b) % 2 ** 256
            mdl = (d - C0) % 2 ** 256 if a < b else d
            cases.append((b - a > p, o % p != (a - b) % p, o == mdl))
    tally("fs_on_fold_output", cases)
    cases = []
    for jj, base in ((1, Hb + 456), (2, Hb + 584)):
        b = pow(BETA, jj, p)
        rows = [b * pow(2, 64 * k, p) % p for k in range(4)]
        binv = pow(b, -1, p)
        for kind in (0, 1, 2):      # beta^j x below 2^66 C0, above 2^99, random
            for _ in range(16):
                T = (rng.randrange(2 ** 66 * C0), rng.randrange(2 ** 99, 2 ** 120), rng.randrange(p))[kind]
                x = T * binv % p
                mem = {}
                for k in range(4):
                    put(mem, base + 32 * k, rows[k])
                regs = as_regs(x, ("x20", "x21", "x22", "x23"))
                regs.update(x29=base)
                o = reg_int(run_site(prog, labels, "t_kmul", regs, mem), ("x24", "x25", "x26", "x27"))
                Z = sum(w * rows[k] for k, w in enumerate(limbs(x)))
                t = fold(Z)
                cases.append((t >= p, o != T, o == t % 2 ** 256))
    tally("kmul", cases)
    edge = [0, 1, p - 1, C0, 2 ** 255]
    pairs = [(a, b) for a in edge for b in edge] + [(rng.randrange(p), rng.randrange(p)) for _ in range(48)]
    weak_b = [p, p + 1, 2 ** 256 - 1, 2 ** 256 - C0] + [rng.randrange(p - C0, p) + C0 for _ in range(12)]
    fs = fa = fw = True
    for a, b in pairs:
        regs = as_regs(a, ("x5", "x6", "x7", "x8"))
        regs.update(as_regs(b, ("x19", "x20", "x21", "x22")))
        fs &= reg_int(run_site(prog, labels, "t_fs", regs), ("x9", "x10", "x11", "x12")) == (a - b) % p
        regs = as_regs(a, ("x5", "x6", "x7", "x8"))
        regs.update(as_regs(b, ("x14", "x15", "x16", "x17")))
        fa &= reg_int(run_site(prog, labels, "t_faddw", regs), ("x23", "x24", "x25", "x26")) % p == (a + b) % p
    for a in edge + [rng.randrange(p) for _ in range(8)]:
        for b in weak_b:
            regs = as_regs(a, ("x5", "x6", "x7", "x8"))
            regs.update(as_regs(b, ("x14", "x15", "x16", "x17")))
            fw &= reg_int(run_site(prog, labels, "t_faddw", regs), ("x23", "x24", "x25", "x26")) % p == (a + b) % p
    res["fs_canonical_exact"] = dict(cases=len(pairs), exact=fs)
    res["faddw_exact"] = dict(cases=len(pairs), exact=fa)
    res["faddw_canonical_plus_weak_exact"] = dict(cases=len(weak_b) * (len(edge) + 8), exact=fw)
    weak = [0, 1, C0, p - 1, p, p + 1, 2 ** 256 - 1] + [rng.getrandbits(256) for _ in range(12)]
    mul_in = [(a, b) for a in weak[:7] for b in weak[:7]] + [(a, rng.getrandbits(256)) for a in weak[7:]]
    for _ in range(12):             # products with a small residue, which take the final carry of the reduction
        a = rng.randrange(1, p)
        mul_in.append((a, rng.randrange(1, 2 ** 64) * pow(a, -1, p) % p))
    sqr_in = weak + [p - rng.randrange(1, 2 ** 100) for _ in range(12)] + [2 ** 256 - rng.randrange(1, 2 ** 40)
                                                                         for _ in range(6)]
    while len(sqr_in) < len(weak) + 42:     # square roots of small residues reach the final +C0 of the reduction
        s_ = rng.randrange(1, 2 ** 64)
        if pow(s_, (p - 1) // 2, p) == 1:
            sqr_in.append(pow(s_, (p + 1) // 4, p))
    ok_m = ok_s = True
    for a, b in mul_in:
        mem = {}
        put(mem, 0x1000, a)
        put(mem, 0x1100, b)
        a64.a64_run(fprog, flabels, "fe_mul_f", dict(x0=0x2000, x1=0x1000, x2=0x1100, sp=0x100000, x30=SENTINEL), mem)
        ok_m &= get(mem, 0x2000) % p == a * b % p
    for a in sqr_in:
        mem = {}
        put(mem, 0x1000, a)
        a64.a64_run(fprog, flabels, "fe_sqr_f", dict(x0=0x2000, x1=0x1000, sp=0x100000, x30=SENTINEL), mem)
        ok_s &= get(mem, 0x2000) % p == a * a % p
    def final_carry(v):
        """Whether reducing the 512-bit v takes the routines' final +C0 (the second fold reaches 2^256). returns: bool."""
        t = (v % 2 ** 256) + (v >> 256) * C0
        return (t % 2 ** 256) + (t >> 256) * C0 >= 2 ** 256
    res["fe_mul_f_exact"] = dict(cases=len(mul_in), exact=ok_m, final_carry_cases=sum(final_carry(a * b) for a, b in mul_in))
    res["fe_sqr_f_exact"] = dict(cases=len(sqr_in), exact=ok_s, final_carry_cases=sum(final_carry(a * a) for a in sqr_in))
    res["all_ok"] = all(v.get("deviates_iff_predicate", True) and v.get("model", True) and v.get("exact", True)
                        and v.get("predicate_holds", 1) > 0 for v in res.values() if isinstance(v, dict))
    return res


def degenerate_batch(prog, labels, X, Ys):
    """fwd and bwd on a batch in which some Y_l = +-phi^t(X) (D_l = 0), with R = acc^(p-2) as fe_inv returns it.

    prog, labels: loaded leaf.s; X: point; Ys: points. returns: (acc, number of exact leaves among the other Y_l).
    """
    mem = {}
    Hb, Yb, Sb, SP = 0x100000, 0x200000, 1 << 34, 0x7000000
    batch_tables(mem, Hb, Yb, X, Ys)
    regs = dict(x0=len(Ys), x1=Yb, x3=Sb, x25=1, x26=0, x27=0, x28=0, sp=SP, x30=SENTINEL)
    regs.update(as_regs(X[1] * X[1] % p, ("x19", "x20", "x21", "x22")))
    out = {}
    a64.a64_run(prog, labels, "fwd", regs, mem, out_regs=out)
    acc = reg_int(out, ("x25", "x26", "x27", "x28"))
    regs = dict(x1=Yb + 160 * len(Ys) + 32, x2=Hb, x3=Sb + 256 * (len(Ys) - 1), sp=SP, x30=SENTINEL)
    regs.update(as_regs(pow(acc, p - 2, p), ("x25", "x26", "x27", "x28")))
    a64.a64_run(prog, labels, "bwd", regs, mem)
    good = 0
    for i, Y in enumerate(Ys):
        if (Y[1] * Y[1] - X[1] * X[1]) % p:
            six = six_leaves(X, Y)
            good += sum(get(mem, Sb + 256 * i + 32 * c) == six[c] for c in range(6))
    return acc, good


def part_a_kernels(rng, cert):
    """Execute leaf.s and fe_inv on production points; derive the per-incidence and fixed counts.

    rng: Random; cert: certificate module. returns: dict.
    """
    res = {}
    prog, labels = a64.a64_load(open(os.path.join(HERE, "kernels", "leaf.s")).read())
    X = PROD.mul(rng.randrange(1, r), GP)
    runs = {}
    for q2 in (1, 2, 5, 9):
        Ys = [PROD.mul(rng.randrange(1, r), GP) for _ in range(q2)]
        runs[q2] = leaf_batch(prog, labels, X, Ys, cert)
    res["runs"] = runs
    res["fwd_per_incidence"] = sorted({(runs[b]["fwd"] - runs[a]["fwd"]) // (b - a) for a, b in ((1, 5), (5, 9), (2, 9))})
    res["bwd_per_incidence"] = sorted({(runs[b]["bwd"] - runs[a]["bwd"]) // (b - a) for a, b in ((1, 5), (5, 9), (2, 9))})
    res["fwd_fixed"] = sorted({d["fwd"] - res["fwd_per_incidence"][0] * q for q, d in runs.items()})
    res["bwd_fixed"] = sorted({d["bwd"] - res["bwd_per_incidence"][0] * q for q, d in runs.items()})
    res["keys_formula"] = ("keys = KEY_FIX + KEY_INC*q2 + KMUL_PATH*products - (1 + EMIT per window-dropped leaf) - "
                           "(EMIT - EMIT_CHECK + [j != 0] per leaf dropped at a full slot)")
    res["keys_formula_holds"] = all(d["keys_formula"] for d in runs.values())
    res["all_correct"] = all(d["acc_ok"] and d["leaf_mismatches"] == 0 and d["messages_ok"] and d["products"] ==
                             d["products_ref"] and d["port_fixed"] for d in runs.values())
    forms = [BETA, 1 + BETA, 1 - BETA, 2 + BETA, 1 + 2 * BETA]
    Ys = [PROD.mul(rng.randrange(1, r), GP) for _ in range(3)]
    over = [[window_x(rng, forms[(3 * i + c) % 5]) if (i + c) % 3 == 0 else rng.randrange(1, p) for c in range(6)]
            for i in range(3)]
    w = leaf_batch(prog, labels, X, Ys, cert, entries_override=over)
    res["window_batch"] = dict(leaves=18, windows=w["windows"], messages=w["messages"], messages_ok=w["messages_ok"],
                               keys=w["keys"], keys_formula=w["keys_formula"])
    fprog, flabels = a64.a64_load(open(os.path.join(HERE, "kernels", "secp_fe.s")).read())
    cnts, ok = set(), True
    for _ in range(3):
        a = rng.randrange(1, p)
        mem = {}
        put(mem, 0x1000, a)
        cnts.add(a64.a64_run(fprog, flabels, "fe_inv", dict(x0=0x2000, x1=0x1000, sp=0x100000, x30=SENTINEL), mem,
                             limit=10 ** 6))
        ok &= get(mem, 0x2000) % p == pow(a, p - 2, p)
    res["fe_inv"] = dict(counts=sorted(cnts), correct=ok)
    src = open(os.path.join(HERE, "kernels", "leaf.s")).read()
    names = ("SHX", "FDX", "F5")
    counts = {k: macro_counts(src, k, names) for k in ("fwd", "bwd", "keys")}
    tot = {n: sum(c[n] for c in counts.values()) for n in names}
    res["site_counts"] = dict(counts, per_incidence=tot)
    res["sites_match_certificate"] = res["site_counts"]["per_incidence"] == cert.SITES_PER_INCIDENCE
    srng = random.Random(SEED + 15485863)   # own stream: the census is unchanged
    res["site_checks"] = site_checks(srng, fprog, flabels)
    Yd = [PROD.mul(srng.randrange(1, r), GP) for _ in range(3)]
    phX = (BETA * X[0] % p, X[1])
    acc0, good = degenerate_batch(prog, labels, X, [Yd[0], phX, Yd[1], PROD.neg(X), Yd[2]])
    res["d_zero_batch"] = dict(acc_zero=acc0 % p == 0, exact_leaves=good, all_other_leaves_wrong=good == 0)
    res["final_R_one_and_unused"] = all(d["final_R_is_one"] and d["messages_ok"] for d in runs.values())
    prog_r, lab_r = a64.a64_load(open(os.path.join(HERE, "kernels", "route.s")).read())
    res["program_instructions"] = dict(leaf=len(prog), route=len(prog_r), secp_fe=len(fprog),
                                       total=len(prog) + len(prog_r) + len(fprog))
    return res


def scan_model(keys, cap):
    """Path model of route.s scan (paired comparisons, one end test per pair). keys: (K2, K3) in pointer order;
    cap: candidate-list entries. returns: (instruction count, candidate index pairs)."""
    n = len(keys)
    cands = []
    cnt = 5 + 2                      # ldr, ldp, sub, add, tbz; cmp, b.hs
    if n == 1:
        return cnt + 1, cands        # ret
    i = 0
    if (n - 1) % 2 == 1:
        cnt += 5                     # odd first comparison: ldr, ldp, cmp, ccmp, b.eq
        if keys[0] == keys[1]:
            if len(cands) >= cap:
                return cnt + 3, cands
            cands.append((0, 1))
            cnt += 3                 # cmp, b.hs, stp
        cnt += 3                     # mov, mov, mov
        if keys[0] != keys[1]:
            cnt += 1                 # b 2f
        i = 1
    while i + 2 <= n - 1:
        cnt += 5                     # ldr, ldp, cmp, ccmp, b.eq
        if keys[i] == keys[i + 1]:
            if len(cands) >= cap:
                return cnt + 3, cands
            cands.append((i, i + 1))
            cnt += 4                 # cmp, b.hs, stp, b
        cnt += 5
        if keys[i + 1] == keys[i + 2]:
            if len(cands) >= cap:
                return cnt + 3, cands
            cands.append((i + 1, i + 2))
            cnt += 4
        cnt += 2                     # cmp, b.lo
        i += 2
    return cnt + 1, cands


class Node:
    """One node of the 4 x 4 machine: its memory, at the same addresses on every node (proof.md 5)."""
    P1, P2, HB, SB, SP = 0x10000, 0x20000, 0x100000, 1 << 34, 0x7000000
    R1, E1, R2, E2 = 0x1000000, 0x2000000, 0x3000000, 0x4000000
    CT, C4, BK, PA, CL = 0x5000000, 0x5100000, 0x6000000, 0x6100000, 0x6200000

    def __init__(self):
        """Empty memory. returns: None."""
        self.mem = {}


def planted_x(rg, d1, d2, K=None):
    """A field element whose key has digits d1 = K0 mod 2^36 and d2 = K1 mod 2^36 (or is the given key K), with a class
    index that the KEY rule gets right.

    rg: Random; d1, d2: digits; K: a key, or None for a fresh one.
    returns: (x, K).
    """
    while True:
        Kc = K
        if Kc is None:
            Kc = d1 | rg.getrandbits(28) << 36 | (d2 | rg.getrandbits(28) << 36) << 64 | rg.getrandbits(128) << 128
            if Kc >= p or exact_key(Kc)[0] != Kc:
                continue
        j = rg.randrange(3)
        x = Kc * pow(BETA, (3 - j) % 3, p) % p
        if rule_j(x) == j and exact_key(x) == (Kc, j):
            return x, Kc


def run_k(prog, labels, entry, regs, nd, calls, formula, port=None):
    """Run one kernel on a node (outermost RET through the sentinel x30) and record (entry, count, formula).

    prog, labels: program; entry: label; regs: registers; nd: Node; calls: list; formula: the certificate's count;
    port: (address, Port) or None. returns: final registers.
    """
    o = {}
    regs = dict(regs, sp=nd.SP, x30=SENTINEL)
    dev = (port[0], 48, port[1].store) if port else None
    calls.append((entry, a64.a64_run(prog, labels, entry, regs, nd.mem, out_regs=o, dev=dev), formula))
    return o


def part_b(cert):
    """4 x 4 machine (proof.md 2.3, 4): keys, fin, h2s, fin, bkt, prefix4, pass4 and scan against an exact model of the
    same order, each message delivered by running dlv; every call's count against its formula. cert: certificate.
    returns: dict."""
    rg = random.Random(SEED + 32452843)   # own stream: the census is unchanged
    prog, labels = a64.a64_load(open(os.path.join(HERE, "kernels", "route.s")).read())
    lprog, llabels = a64.a64_load(open(os.path.join(HERE, "kernels", "leaf.s")).read())
    side, q2, cap1, cap2, b3, b4, cap3, ncl = 4, 4, 7, 7, 2, 4, 8, 64
    n3, n4, st1, st2, bkb = 2 ** b3, 2 ** b4, 32 * cap1, 32 * cap2, 16 * cap3
    ids = [(i, j) for i in range(side) for j in range(side)]
    nodes = {n: Node() for n in ids}
    calls = []
    for (i, j), nd in nodes.items():
        tf = cert.TINIT_ENTRY * side + 1
        run_k(prog, labels, "tinit", dict(x0=side, x1=nd.P1 + 64, x2=nd.R1 + i * st1, x3=nd.R1 + (i + 1) * st1), nd,
              calls, tf)
        run_k(prog, labels, "tinit", dict(x0=side, x1=nd.P2 + 64, x2=nd.R2 + j * st2, x3=nd.R2 + (j + 1) * st2), nd,
              calls, tf)
        run_k(prog, labels, "init3", dict(x0=n3, x1=nd.CT, x2=nd.BK, x3=bkb), nd, calls, cert.INIT3_BKT * n3 + 1)
        run_k(prog, labels, "zero4", dict(x0=n4, x1=nd.C4), nd, calls, n4 // 4 * cert.ZERO4_PER4 + 1)
        header(nd.mem, nd.HB)
    pool = [planted_x(rg, rg.randrange(side), rg.randrange(side))[1] for _ in range(10)]
    batches, shared = {}, 0
    for n in ids:
        lst = []
        for _ in range(6 * q2):
            if rg.random() < 0.15:
                lst.append(planted_x(rg, 0, 0, rg.choice(pool))[0])
                shared += 1
            else:
                lst.append(planted_x(rg, rg.randrange(side), rg.randrange(side))[0])
        batches[n] = [lst[6 * l:6 * l + 6] for l in range(q2)]
    # hop 1: keys' emits write each record into the pair slot of its sender row at (d1, j); then the final cursors
    msgs_ok, ref1, drops = True, {n: {i: [] for i in range(side)} for n in ids}, dict(hop1=0, hop2=0, bucket=0)
    fin_f = 3 + cert.FIN_DEST * side
    for (i, j), nd in nodes.items():
        for l, six in enumerate(batches[(i, j)]):
            for c_, x in enumerate(six):
                put(nd.mem, nd.SB + 256 * l + 32 * c_, x)
        W2 = i | (side * i + j) << 36
        exp, cur, prods, saved, wins, full = emit_reference(batches[(i, j)], W2, lambda d: nd.R1 + i * st1,
                                                            lambda d: nd.R1 + (i + 1) * st1, cert)
        port = Port()
        o = run_k(lprog, llabels, "keys", dict(x0=nd.P1, x2=nd.HB, x3=nd.SB + 256 * (q2 - 1), x17=W2, x18=0), nd, calls,
                  cert.KEY_FIX + cert.KEY_INC * q2 + cert.KMUL_PATH * prods - saved, (nd.P1, port))
        msgs_ok &= port.msgs == exp and o["x18"] == prods and wins == 0
        drops["hop1"] += full
        for m in port.msgs:
            link(prog, labels, nodes[(m[0], j)].mem, m, calls, cert)
        for m in exp:
            ref1[(m[0], j)][i].append(tuple(m[2:]))
        port = Port()
        run_k(prog, labels, "fin", dict(x0=nd.P1, x1=side, x2=nd.E1 + 32 * i), nd, calls, fin_f, (nd.P1, port))
        msgs_ok &= port.msgs == [[d, nd.E1 + 32 * i, cur.get(d, nd.R1 + i * st1), 0, 0, 0] for d in range(side)]
        for m in port.msgs:
            link(prog, labels, nodes[(m[0], j)].mem, m, calls, cert)
    # barrier; hop 2: h2s streams the pair slots through the row port to (own row, d2)
    ref2 = {n: {j: [] for j in range(side)} for n in ids}
    slots_ok = True
    for (a, j), nd in nodes.items():
        for i in range(side):
            got = [tuple(nd.mem[x + 8 * k] for k in range(4)) for x in range(nd.R1 + i * st1, nd.mem[nd.E1 + 32 * i], 32)]
            slots_ok &= got == ref1[(a, j)][i]
        used, stored, nrec = {}, 0, 0
        for i in range(side):
            for rec in ref1[(a, j)][i]:
                nrec += 1
                d2 = rec[3] & M36
                if used.get(d2, 0) < cap2:
                    used[d2] = used.get(d2, 0) + 1
                    ref2[(a, d2)][j].append(rec)
                    stored += 1
                else:
                    drops["hop2"] += 1
        port = Port()
        loops = sum(len(ref1[(a, j)][i]) // 8 for i in range(side))
        run_k(prog, labels, "h2s", dict(x0=nd.P2, x3=side, x13=nd.E1, x14=nd.R1, x15=st1), nd, calls,
              cert.H2S_FIX + cert.H2S_SLOT * side + cert.H2S_REC * stored + cert.H2S_DROP * (nrec - stored)
              + cert.LOOP8 * loops,
              (nd.P2, port))
        msgs_ok &= len(port.msgs) == stored
        for m in port.msgs:
            link(prog, labels, nodes[(a, m[0])].mem, m, calls, cert)
        port = Port()
        run_k(prog, labels, "fin", dict(x0=nd.P2, x1=side, x2=nd.E2 + 32 * j), nd, calls, fin_f, (nd.P2, port))
        for m in port.msgs:
            link(prog, labels, nodes[(a, m[0])].mem, m, calls, cert)
    # barrier; grouping at each home
    found, exp_c, scan_ok, tags_ok = set(), set(), True, True
    for n, nd in nodes.items():
        for j in range(side):
            got = [tuple(nd.mem[x + 8 * k] for k in range(4)) for x in range(nd.R2 + j * st2, nd.mem[nd.E2 + 32 * j], 32)]
            slots_ok &= got == ref2[n][j]
            tags_ok &= all(rec[2] >> 36 == side * (rec[2] & M36) + j and rec[3] & M36 == n[1] for rec in got)
        bk, nrec = [[] for _ in range(n3)], 0
        for j in range(side):
            for k, rec in enumerate(ref2[n][j]):
                nrec += 1
                t = rec[0] & (n3 - 1)
                if len(bk[t]) < cap3:
                    bk[t].append((nd.R2 + j * st2 + 32 * k, rec))
                else:
                    drops["bucket"] += 1
        seq = [e for t in range(n3) for e in bk[t]]
        loops = sum(len(ref2[n][j]) // 8 for j in range(side))
        run_k(prog, labels, "bkt", dict(x0=side, x13=nd.E2, x14=nd.R2, x15=st2, x3=n3 - 1, x4=nd.CT, x5=n4 - 1,
                                        x6=nd.C4), nd, calls,
              cert.BKT_FIX + cert.BKT_SLOT * side + cert.BKT_REC * len(seq) + cert.BKT_DROP * (nrec - len(seq))
              + cert.LOOP8 * loops)
        scan_ok &= all([nd.mem[nd.BK + t * bkb + 16 * k], nd.mem[nd.BK + t * bkb + 16 * k + 8]] ==
                       [e[0], e[1][1] & (n4 - 1)] for t in range(n3) for k, e in enumerate(bk[t]))
        run_k(prog, labels, "prefix4", dict(x0=n4, x1=nd.C4, x2=nd.PA), nd, calls, n4 // 4 * cert.PREFIX4_PER4 + 1)
        run_k(prog, labels, "pass4", dict(x0=n3, x1=nd.CT, x2=nd.BK, x3=bkb, x4=nd.C4), nd, calls,
              cert.PASS4_BKT * n3 + cert.PASS4_REC * len(seq) + cert.LOOP8 * sum(len(b) // 8 for b in bk) + 1)
        order = sorted(seq, key=lambda e: e[1][1] & (n4 - 1))
        scan_ok &= [nd.mem[nd.PA + 8 * k] for k in range(len(seq))] == [e[0] for e in order]
        if len(seq) >= 2:
            mc, cands = scan_model([e[1][:2] for e in order], ncl)
            o = run_k(prog, labels, "scan", dict(x0=len(seq), x1=nd.PA, x2=nd.CL, x3=nd.CL + 16 * ncl), nd, calls, mc)
            for c in range(nd.CL, o["x2"], 16):
                found.add((n, nd.mem[c], nd.mem[c + 8]))
            exp_c |= {(n, order[a][0], order[b][0]) for a, b in cands}
    groups = {}
    for n in ids:
        for six in batches[n]:
            for x in six:
                groups.setdefault(exact_key(x)[0], []).append(n)
    cnt, form = {}, {}
    for e, c, f in calls:
        cnt[e], form[e] = cnt.get(e, 0) + c, form.get(e, 0) + f
    res = dict(nodes=len(ids), leaves=6 * q2 * len(ids), shared_key_leaves=shared,
               equal_key_groups=sum(len(v) >= 2 for v in groups.values()), candidates=len(found),
               candidates_equal_exact_order_reference=found == exp_c and len(found) > 0,
               messages_equal_emit_reference=msgs_ok, slots_equal_reference=slots_ok, source_tags_ok=tags_ok,
               buckets_pointer_order_and_scan_ok=scan_ok, drops=drops, kernel_calls=len(calls), counts=cnt,
               deliveries=sum(e == "dlv" for e, _, _ in calls), formulas_ok=all(c == f for _, c, f in calls))
    res["full_paths"] = full_paths(prog, labels, lprog, llabels, rg, cert)
    res["unrolled_paths"] = unrolled_paths(prog, labels, rg, cert)
    return res


def full_paths(prog, labels, lprog, llabels, rg, cert):
    """Paths forced on small inputs: keys and h2s with every slot full, bkt with a full bucket, scan with a full list.

    prog, labels: route.s; lprog, llabels: leaf.s; rg: Random; cert: certificate module. returns: dict of checks.
    """
    out, calls = {}, []
    nd = Node()
    header(nd.mem, nd.HB)
    leaves = [[planted_x(rg, rg.randrange(4), 0)[0] for _ in range(6)] for _ in range(3)]
    for l, six in enumerate(leaves):
        for c_, x in enumerate(six):
            put(nd.mem, nd.SB + 256 * l + 32 * c_, x)
    run_k(prog, labels, "tinit", dict(x0=4, x1=nd.P1 + 64, x2=nd.R1, x3=nd.R1), nd, calls, cert.TINIT_ENTRY * 4 + 1)
    exp, _, prods, saved, _, full = emit_reference(leaves, 5, lambda d: nd.R1, lambda d: nd.R1, cert)
    port = Port()
    o = run_k(lprog, llabels, "keys", dict(x0=nd.P1, x2=nd.HB, x3=nd.SB + 256 * 2, x17=5, x18=0), nd, calls,
              cert.KEY_FIX + cert.KEY_INC * 3 + cert.KMUL_PATH * prods - saved, (nd.P1, port))
    out["keys_every_slot_full"] = dict(count=calls[-1][1], formula=calls[-1][2], messages=len(port.msgs), dropped=full,
                                       products_ok=o["x18"] == prods)
    nd = Node()
    recs = [[rg.getrandbits(64), rg.getrandbits(64), 7, rg.randrange(4)] for _ in range(5)]
    for k, rec in enumerate(recs):
        for w in range(4):
            nd.mem[nd.R1 + 32 * k + 8 * w] = rec[w]
    nd.mem[nd.E1] = nd.R1 + 32 * len(recs)
    run_k(prog, labels, "tinit", dict(x0=4, x1=nd.P2 + 64, x2=nd.R2, x3=nd.R2), nd, calls, cert.TINIT_ENTRY * 4 + 1)
    port = Port()
    run_k(prog, labels, "h2s", dict(x0=nd.P2, x3=1, x13=nd.E1, x14=nd.R1, x15=32 * 8), nd, calls,
          cert.H2S_FIX + cert.H2S_SLOT + cert.H2S_DROP * 5, (nd.P2, port))
    out["h2s_every_slot_full"] = dict(count=calls[-1][1], formula=calls[-1][2], messages=len(port.msgs))
    nd = Node()
    run_k(prog, labels, "init3", dict(x0=4, x1=nd.CT, x2=nd.BK, x3=16), nd, calls, cert.INIT3_BKT * 4 + 1)
    run_k(prog, labels, "zero4", dict(x0=4, x1=nd.C4), nd, calls, cert.ZERO4_PER4 + 1)
    for k in range(3):
        rec = [4 * rg.getrandbits(40) + 1, 4 * rg.getrandbits(40) + 2, rg.getrandbits(64), rg.getrandbits(64)]
        for w in range(4):
            nd.mem[nd.R2 + 32 * k + 8 * w] = rec[w]
    nd.mem[nd.E2] = nd.R2 + 96
    run_k(prog, labels, "bkt", dict(x0=1, x13=nd.E2, x14=nd.R2, x15=32 * 8, x3=3, x4=nd.CT, x5=3, x6=nd.C4), nd, calls,
          cert.BKT_FIX + cert.BKT_SLOT + cert.BKT_REC + 2 * cert.BKT_DROP)
    out["bkt_full_bucket"] = dict(count=calls[-1][1], formula=calls[-1][2],
                                  d4_counted_once=sum(nd.mem[nd.C4 + 8 * k] for k in range(4)) == 1,
                                  entry_ok=[nd.mem[nd.BK + 16], nd.mem[nd.BK + 24]] == [nd.R2, 2])
    nd = Node()
    keys = [(1, 1), (1, 1), (2, 2), (3, 3), (3, 3), (4, 4)]
    for k, kk in enumerate(keys):
        nd.mem[nd.R2 + 32 * k], nd.mem[nd.R2 + 32 * k + 8] = kk
        nd.mem[nd.PA + 8 * k] = nd.R2 + 32 * k
    mc, mcand = scan_model(keys, 1)
    o = run_k(prog, labels, "scan", dict(x0=6, x1=nd.PA, x2=nd.CL, x3=nd.CL + 16), nd, calls, mc)
    out["scan_full_list"] = dict(count=calls[-1][1], model=mc, listed=(o["x2"] - nd.CL) // 16, model_listed=len(mcand))
    k1, h2, b3 = out["keys_every_slot_full"], out["h2s_every_slot_full"], out["bkt_full_bucket"]
    out["ok"] = (all(c == f for _, c, f in calls) and k1["messages"] == 0 and k1["dropped"] == 18 and k1["products_ok"]
                 and h2["messages"] == 0 and b3["d4_counted_once"] and b3["entry_ok"]
                 and out["scan_full_list"]["listed"] == 1)
    return out


def unrolled_paths(prog, labels, rg, cert):
    """h2s, bkt and pass4 on slots and buckets of 0 to 19 records (every 1-, 2- and 4-record block combination, up to two
    8-record passes), some full; outputs against exact models and counts against the certificate's formulas.

    prog, labels: route.s; rg: Random; cert: certificate module. returns: dict of checks.
    """
    calls, out, ns = [], {}, list(range(20))
    rg.shuffle(ns)
    nd = Node()                     # h2s: 20 slots of 32*19 bytes; destination d takes at most 24 + 12d records
    st, ends = 32 * 19, {d: nd.R2 + 32 * (24 + 12 * d) for d in range(4)}
    recs = []
    for s_, n in enumerate(ns):
        cur = []
        for k in range(n):
            rec = [rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(28) << 36 | rg.randrange(4)]
            for w in range(4):
                nd.mem[nd.R1 + st * s_ + 32 * k + 8 * w] = rec[w]
            cur.append(rec)
        nd.mem[nd.E1 + 32 * s_] = nd.R1 + st * s_ + 32 * n
        recs.append(cur)
    for d in range(4):
        nd.mem[nd.P2 + 64 + 16 * d], nd.mem[nd.P2 + 72 + 16 * d] = nd.R2, ends[d]
    curs, exp = {d: nd.R2 for d in range(4)}, []
    for cur in recs:
        for rec in cur:
            d = rec[3] & M36
            if curs[d] < ends[d]:
                exp.append([d, curs[d]] + rec)
                curs[d] += 32
    port = Port()
    run_k(prog, labels, "h2s", dict(x0=nd.P2, x3=len(ns), x13=nd.E1, x14=nd.R1, x15=st), nd, calls,
          cert.H2S_FIX + cert.H2S_SLOT * len(ns) + cert.H2S_REC * len(exp) + cert.H2S_DROP * (sum(ns) - len(exp))
          + cert.LOOP8 * sum(n // 8 for n in ns), (nd.P2, port))
    out["h2s"] = dict(records=sum(ns), sent=len(exp), messages_ok=port.msgs == exp,
                      cursors_ok=all(nd.mem[nd.P2 + 64 + 16 * d] == curs[d] for d in range(4)))
    nd = Node()                     # bkt: 20 slots; d3 bucket t holds at most 24 + 12t entries
    n3, n4, bkb = 4, 4, 16 * 64
    run_k(prog, labels, "init3", dict(x0=n3, x1=nd.CT, x2=nd.BK, x3=bkb), nd, calls, cert.INIT3_BKT * n3 + 1)
    for t in range(n3):
        nd.mem[nd.CT + 16 * t + 8] = nd.BK + bkb * t + 16 * (24 + 12 * t)
    run_k(prog, labels, "zero4", dict(x0=n4, x1=nd.C4), nd, calls, cert.ZERO4_PER4 + 1)
    bk, cnt4, stored = [[] for _ in range(n3)], [0] * n4, 0
    for s_, n in enumerate(ns):
        for k in range(n):
            rec = [rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(64)]
            a = nd.R2 + st * s_ + 32 * k
            for w in range(4):
                nd.mem[a + 8 * w] = rec[w]
            t = rec[0] & (n3 - 1)
            if len(bk[t]) < 24 + 12 * t:
                bk[t].append((a, rec[1] & (n4 - 1)))
                cnt4[rec[1] & (n4 - 1)] += 1
                stored += 1
        nd.mem[nd.E2 + 32 * s_] = nd.R2 + st * s_ + 32 * n
    run_k(prog, labels, "bkt", dict(x0=len(ns), x13=nd.E2, x14=nd.R2, x15=st, x3=n3 - 1, x4=nd.CT, x5=n4 - 1, x6=nd.C4),
          nd, calls, cert.BKT_FIX + cert.BKT_SLOT * len(ns) + cert.BKT_REC * stored + cert.BKT_DROP * (sum(ns) - stored)
          + cert.LOOP8 * sum(n // 8 for n in ns))
    ok_b = all([nd.mem[nd.BK + bkb * t + 16 * k], nd.mem[nd.BK + bkb * t + 16 * k + 8]] == list(e)
               for t in range(n3) for k, e in enumerate(bk[t]))
    ok_b &= all(nd.mem[nd.CT + 16 * t] == nd.BK + bkb * t + 16 * len(bk[t]) for t in range(n3))
    ok_b &= [nd.mem[nd.C4 + 8 * v] for v in range(n4)] == cnt4
    out["bkt"] = dict(records=sum(ns), bucketed=stored, entries_cursors_and_counts_ok=ok_b)
    nd = Node()                     # pass4: 20 buckets of 0 to 19 entries
    ents, bkb = [], 16 * 19
    for t, n in enumerate(ns):
        cur = [(rg.getrandbits(64), rg.randrange(n4)) for _ in range(n)]
        for k, (a, v) in enumerate(cur):
            nd.mem[nd.BK + bkb * t + 16 * k], nd.mem[nd.BK + bkb * t + 16 * k + 8] = a, v
        nd.mem[nd.CT + 16 * t] = nd.BK + bkb * t + 16 * n
        ents.append(cur)
    tot = [sum(v == w for c in ents for _, v in c) for w in range(n4)]
    for w in range(n4):
        nd.mem[nd.C4 + 8 * w] = nd.PA + 8 * sum(tot[:w])
    run_k(prog, labels, "pass4", dict(x0=len(ns), x1=nd.CT, x2=nd.BK, x3=bkb, x4=nd.C4), nd, calls,
          cert.PASS4_BKT * len(ns) + cert.PASS4_REC * sum(ns) + cert.LOOP8 * sum(n // 8 for n in ns) + 1)
    order = [a for w in range(n4) for c in ents for a, v in c if v == w]
    out["pass4"] = dict(entries=sum(ns), pointer_order_ok=[nd.mem[nd.PA + 8 * k] for k in range(sum(ns))] == order)
    out["counts_equal_formulas"] = all(c == f for _, c, f in calls)
    out["ok"] = (out["counts_equal_formulas"] and out["h2s"]["messages_ok"] and out["h2s"]["cursors_ok"]
                 and out["bkt"]["entries_cursors_and_counts_ok"] and out["pass4"]["pointer_order_ok"]
                 and out["h2s"]["sent"] < sum(ns) and out["bkt"]["bucketed"] < sum(ns))
    return out


def toy_setup(rng):
    """The q25 toy curve of package F (q = 33556867, r = 33563419, j = 0): generator, beta, lambda.
    rng: Random. returns: (Curve, G, beta, lambda)."""
    q, order = 33556867, 33563419
    E = Curve(q, order)
    while True:
        x = rng.randrange(q)
        y = pow(x ** 3 + 7, (q + 1) // 4, q)
        if (y * y - x ** 3 - 7) % q == 0:
            G = (x, y)
            break
    assert E.mul(order, G) is None
    cr = lambda m: [c for g in range(2, 50) for c in [pow(g, (m - 1) // 3, m)] if c != 1][0]
    for b in (cr(q), cr(q) ** 2 % q):
        for l in (cr(order), cr(order) ** 2 % order):
            if E.mul(l, G) == (G[0] * b % q, G[1]):
                return E, G, b, l
    raise AssertionError("no beta/lambda pair")


def toy_trial_logs(rng, rt, lam, nodes, q1, q2):
    """One toy trial in discrete-log form (proof.md 8): class of [t]G is t^6 mod r; leaf phi^t(X) + sY has logarithm
    lam^t x + s y and Q coefficient lam^t b.

    rng: Random; rt: group order; lam: cube root of unity; nodes, q1, q2: machine shape.
    returns: (scalars, success over all pairs, success over cross-node pairs, colliding pairs).
    """
    xi = rng.randrange(rt)
    sc = []
    for g in range(nodes):
        bg = rng.randrange(rt)
        xs = [(rng.randrange(rt), bg) for _ in range(q1)]
        ys = [rng.randrange(1, rt) for _ in range(q2)]
        sc.append((xs, ys))
    units = [pow(lam, t, rt) for t in range(3)]
    keys, live = [], []
    for g, (xs, ys) in enumerate(sc):
        for jx, (a, b) in enumerate(xs):
            tx = (a + b * xi) % rt
            if tx == 0 or any((c - e * w * tx) % rt == 0 for c in ys for w in units for e in (1, -1)):
                continue  # X = O or D = 0 (Y = +-phi^t(X)): the batch is a defect row (proof.md 6.4)
            live.append((g, b, tx, ys))
            keys += [pow((w * tx + e * c) % rt, 6, rt) for c in ys for w in units for e in (1, -1)]
    groups = {}
    if len(set(keys)) < len(keys):
        for g, b, tx, ys in live:
            for c in ys:
                for w in units:
                    for e in (1, -1):
                        t = (w * tx + e * c) % rt
                        groups.setdefault(pow(t, 6, rt), []).append((g, w * b % rt, t))
    succ = succ_x = False
    pairs = 0
    for v in groups.values():
        if len(v) < 2:
            continue
        for u_ in range(len(v)):
            for w_ in range(u_ + 1, len(v)):
                (ga, ba, ta), (gb, bb, tb) = v[u_], v[w_]
                pairs += 1
                u = tb * pow(ta, -1, rt) % rt
                if (bb - u * ba) % rt:
                    succ = True
                    succ_x |= ga != gb
    return (xi, sc), succ, succ_x, pairs


TOY_TAB = {}


def toy_table(E, G):
    """Fixed-base table of G: rows of [d * 256^w]G. E: curve; G: generator. returns: list of rows (cached)."""
    if G not in TOY_TAB:
        rows, P = [], G
        for _ in range(4):
            row = [None]
            for _ in range(255):
                row.append(E.add(row[-1], P))
            rows.append(row)
            P = E.add(row[-1], P)
        TOY_TAB[G] = rows
    return TOY_TAB[G]


def toy_fmul(E, tab, k):
    """[k]G by the fixed-base table. E: curve; tab: toy_table rows; k: scalar. returns: point."""
    k %= E.r
    acc, w = None, 0
    while k:
        acc = E.add(acc, tab[w][k & 255])
        k >>= 8
        w += 1
    return acc


def toy_trial_curve(E, G, beta, lam, scalars, nodes, q1, q2):
    """The same trial on the toy curve: six leaves per incidence from one batched inverse of D = y_Y^2 - y_X^2 (proof.md
    3.1), keys min(x, bx, b^2x), each equal-key pair resolved with the points, solved and checked.

    E, G, beta, lam: toy curve data; scalars: from toy_trial_logs; nodes, q1, q2: shape.
    returns: (success, equal-key pairs).
    """
    xi, sc = scalars
    q, rt = E.q, E.r
    tab = toy_table(E, G)
    fm = lambda k: toy_fmul(E, tab, k)
    Q = fm(xi)
    b2 = beta * beta % q
    groups = {}
    for g, (xs, ys) in enumerate(sc):
        Ypts = [fm(c) for c in ys]
        for jx, (a, b) in enumerate(xs):
            X = fm(a + b * xi)
            if X is None or any((Y[1] * Y[1] - X[1] * X[1]) % q == 0 for Y in Ypts):
                continue
            Ds = [(Y[1] * Y[1] - X[1] * X[1]) % q for Y in Ypts]
            pre = [1]
            for d in Ds:
                pre.append(pre[-1] * d % q)
            R = pow(pre[-1], -1, q)
            for iy in range(q2 - 1, -1, -1):
                inv = R * pre[iy] % q
                R = R * Ds[iy] % q
                Y = Ypts[iy]
                u, v = (Y[1] - X[1]) * inv % q, (Y[1] + X[1]) * inv % q
                for t in range(3):
                    bt = pow(beta, t, q)
                    Qt = (Y[0] * Y[0] + bt * X[0] * Y[0] + bt * bt * X[0] * X[0]) % q
                    for s, w in ((1, u), (-1, v)):
                        x = (Qt * Qt * w * w - bt * X[0] - Y[0]) % q
                        key = min(x, x * beta % q, x * b2 % q)
                        groups.setdefault(key, []).append((g, t, s, a, b, ys[iy]))
    succ, pairs = False, 0
    for v in groups.values():
        for u_ in range(len(v)):
            for w_ in range(u_ + 1, len(v)):
                pairs += 1
                pts = []
                for (g, t, s, a, b, c) in (v[u_], v[w_]):
                    lt = pow(lam, t, rt)
                    P = E.add(E.add(fm(lt * a), E.mul(lt * b % rt, Q)), fm(s * c))
                    pts.append((P, (lt * a + s * c) % rt, lt * b % rt))
                (PA, aA, bA), (PB, aB, bB) = pts
                for k in range(3):
                    for e in (1, -1):
                        img = (PA[0] * pow(beta, k, q) % q, PA[1] if e == 1 else (-PA[1]) % q)
                        if img == PB:
                            u = e * pow(lam, k, rt) % rt
                            den = (bB - u * bA) % rt
                            if den:
                                cand = (u * aA - aB) * pow(den, -1, rt) % rt
                                assert E.mul(cand, G) == Q
                                succ = True
    return succ, pairs


def mu_of(nodes, n0, rt):
    """Cross-node collision mean |Gamma| P1. nodes, n0: shape; rt: order. returns: float."""
    s = nodes * n0
    return (s * s - nodes * n0 * n0) / 2 * (6 * (rt - 1) + 1) / rt ** 2


def part_c(rng):
    """Toy success near the production mu, and curve/log equivalence. rng: Random. returns: dict."""
    E, G, beta, lam = toy_setup(rng)
    rt = E.r
    res = dict(q=E.q, r=rt)
    nodes, q1, q2 = 25, 4, 4
    ok = okx = done = pairs = 0
    dl = deadline(0.70)
    while done < TRIALS and (done < MIN_MULTI or time.time() < dl):
        _, s1, s2, pr = toy_trial_logs(rng, rt, lam, nodes, q1, q2)
        ok += s1
        okx += s2
        pairs += pr
        done += 1
    n0 = 6 * q1 * q2
    S = nodes * n0
    mu_all = S * (S - 1) / 2 * (6 * (rt - 1) + 1) / rt ** 2
    mu = mu_of(nodes, n0, rt)
    ph = ok / max(done, 1)
    p1t = (6 * (rt - 1) + 1) / rt ** 2
    gam = (S * S - nodes * n0 * n0) / 2
    cs = gam * 12 * (q1 + q2) * S * (p1t * p1t + 36 / (rt - 1) ** 2)
    res["multi_node_grid"] = dict(nodes=nodes, grid=[q1, q2], leaves_per_incidence=6, leaves=S, trials=done,
                                  success=round(ph, 5), se=round(math.sqrt(ph * (1 - ph) / max(done, 1)), 5),
                                  cross_node_success=round(okx / max(done, 1), 5),
                                  mean_colliding_pairs=round(pairs / max(done, 1), 5), mu_cross=round(mu, 6),
                                  mu_all_pairs=round(mu_all, 6), poisson=round(1 - math.exp(-mu), 5),
                                  chen_stein_b1_b2_at_this_scale=round(cs, 5))
    agree = n = curve_succ = 0
    dl = deadline(1.0)
    while n < CURVE and (n < MIN_CURVE or time.time() < dl):
        scal, s1, _, pr = toy_trial_logs(rng, rt, lam, nodes, q1, q2)
        s2, pr2 = toy_trial_curve(E, G, beta, lam, scal, nodes, q1, q2)
        agree += (s1 == s2) and (pr == pr2)
        curve_succ += s2
        n += 1
    res["curve_pipeline_vs_logs"] = dict(trials=n, identical_outcomes=agree, curve_successes=curve_succ)
    return res


def load_certificate():
    """Import certificate.py with its printout captured. returns: the module."""
    buf = io.StringIO()
    sys.path.insert(0, HERE)
    with contextlib.redirect_stdout(buf):
        import certificate as cert
    return cert


def part_d(cert):
    """The certificate's verdict, claim and the count constants checked against parts A and B.

    cert: certificate module. returns: dict.
    """
    names = ("LEAVES", "FWD_INC", "FWD_FIX", "BWD_INC", "BWD_FIX", "KEY_INC", "KEY_FIX", "EMIT", "EMIT_INC",
             "EMIT_CHECK", "KMUL_PATH", "INV", "TINIT_ENTRY", "FIN_DEST", "H2S_REC", "H2S_DROP", "H2S_SLOT", "H2S_FIX",
             "BKT_REC", "BKT_DROP", "BKT_SLOT", "BKT_FIX", "PASS4_REC", "PASS4_BKT", "SCAN_REC", "SCAN_FIX", "LOOP8",
             "ZERO4_PER4", "PREFIX4_PER4", "INIT3_BKT", "DLV", "DONE")
    return dict(passed=not cert.FAILED, failed=cert.FAILED, claim=format(float(cert.CLAIM), ".7f"), q1=cert.Q1,
                mu=float(cert.T["mu"]), success_lower_bound=float(cert.T["lo"]), per_leaf=float(cert.PER_LEAF),
                d_batch=float(cert.T["d_batch"]), g_def=float(cert.T["g_def"]),
                counts={k: getattr(cert, k) for k in names})


def main():
    """Run the parts within the time budget and print one JSON object. returns: None."""
    rng = random.Random(SEED)
    cert = load_certificate()
    OUT["input"] = cfg
    OUT["budget"] = dict(seconds=BUDGET, reserve=RESERVE, minimums=dict(class_rule_random=MIN_RANDOM,
                         class_rule_near_tie=MIN_NEAR, multi_node_trials=MIN_MULTI,
                         curve_trials=MIN_CURVE), note="each count runs at least its minimum (or its requested count if "
                         "smaller), then continues while its share of the remaining budget lasts")
    OUT["A_identities"] = part_a_identities(rng)
    if DO_KERNELS:
        OUT["A_kernels"] = a = part_a_kernels(rng, cert)
        OUT["B_routing"] = b = part_b(cert)
        OUT["counts_match_certificate"] = (
            a["fwd_per_incidence"] == [cert.FWD_INC] and a["bwd_per_incidence"] == [cert.BWD_INC]
            and a["fwd_fixed"] == [cert.FWD_FIX] and a["bwd_fixed"] == [cert.BWD_FIX] and a["keys_formula_holds"]
            and a["window_batch"]["keys_formula"] and a["fe_inv"]["counts"] == [cert.INV] and b["formulas_ok"]
            and b["full_paths"]["ok"] and b["unrolled_paths"]["ok"] and cert.SCAN_REC == 6
            and a["sites_match_certificate"])
    OUT["C_toy"] = part_c(rng)
    OUT["D_certificate"] = part_d(cert)
    OUT["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(OUT, indent=1))


if __name__ == "__main__":
    main()
