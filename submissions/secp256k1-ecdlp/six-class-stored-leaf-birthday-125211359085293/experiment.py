"""Stored-leaf birthday search on secp256k1 (proof.md). Part A runs kernels/leaf.s and secp_fe.s on production points
against exact arithmetic, the exact class rule and the fraction class key, with the rare-carry fixups and window drops
forced by planted points and leaves and each site type run alone (sites.py); part B runs the emits, both hops (every
message delivered by route.s dlv) and the grouping on a 4 x 4 machine through kernels/leaf.s and route.s; part C
(toy.py) runs the toy pipeline at q25; part D runs certificate.py. Input JSON on stdin, output JSON."""
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
MIN_CLASSES = 3000                              # class-key checks that always run (a third at window edges)
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
sys.path.insert(0, HERE)
import toy  # noqa: E402
import sites  # noqa: E402
Curve = toy.Curve

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
GW = (BETA << 448) // p                               # G' of the key paths (proof.md 3.3): limbs g0..g6, g1..g6 = GL
GG = [(GW >> (64 * k)) & M64 for k in range(7)]
M48 = 2 ** 48 - 1
OUT = {}


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


def kparts(x):
    """Bit-exact model of the KEY macro (proof.md 3.3). x: integer < 2^256. returns: (j, whether x is in the window
    set W1-W4, where j may differ from the exact index)."""
    xl = limbs(x)
    f1 = sum(((xl[i] * GL[5 - i]) & M64) + ((xl[i] * GL[4 - i]) >> 64) for i in range(4)) & M64
    f0 = xl[3]
    s, d = (f0 + f1) & M64, (f0 - f1) & M64
    m = f1 if f0 >= f1 else f0
    e = m + s
    win = f1 + 16 > M64 or s + 17 > M64 or (e & M64) + 33 > M64 or d < 16
    return (2 if e > M64 else 1 if f0 >= f1 else 0), win


def rule_j(x):
    """KEY's class index. x: integer < 2^256. returns: j in 0, 1, 2."""
    return kparts(x)[0]


def ext_a(x):
    """Bit-exact model of the key extension (proof.md 3.3): columns 4 to 6 of x*G' without column 3's carry, column 6
    being KEY's F1 plus the carry out of column 5. x: integer < 2^256. returns: A1 < 2^192."""
    a = limbs(x)
    f1 = sum(((a[i] * GL[5 - i]) & M64) + ((a[i] * GL[4 - i]) >> 64) for i in range(4)) & M64
    c4 = sum(((a[i] * GG[4 - i]) & M64) + ((a[i] * GG[3 - i]) >> 64) for i in range(4))
    c5 = sum(((a[i] * GG[5 - i]) & M64) + ((a[i] * GG[4 - i]) >> 64) for i in range(4))
    return (c4 + (c5 << 64) + (f1 << 128)) % 2 ** 192


def key_path(x, j):
    """Bit-exact model of keys' path j (proof.md 3.3). x: integer < 2^256; j: KEY's index.
    returns: A, the 192-bit value whose top 176 bits are the key, or None when the path's window drops the leaf."""
    if j == 0:
        return x >> 64
    a1 = ext_a(x)
    if j == 1:
        return None if (a1 + 7) % 2 ** 16 < 8 else a1
    a2 = -((x >> 64) + a1) % 2 ** 192
    return None if a2 % 2 ** 16 < 16 else a2


def key_model(x):
    """The value keys emits for leaf x: its path's A, or None if a window drops the leaf. x: integer < 2^256.
    returns: (A or None, j)."""
    j = rule_j(x)
    return key_path(x, j), j


def class_key(x):
    """The class key: the top 176 bits of the canonical value v* (proof.md 3.3). x: field element. returns: integer."""
    return exact_key(x)[0] >> 80


def frac192(v):
    """floor(2^192 v / p). v: field element. returns: integer."""
    return (v << 192) // p


def deadline(frac):
    """Deadline of a budget-limited loop. frac: share of the remaining budget. returns: absolute time."""
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
    res["key_constants"] = GW >> 64 == GQ and GG[1:] == GL and GG[0] == 0x03fde1631c4b80af
    bad = drops = wrong = done = 0
    kinds = {"random": 0, "near_tie": 0}
    dl = deadline(0.10)
    while done < 26000 and (done < MIN_RANDOM or time.time() < dl):
        if done < 20000:
            x, kind = rng.randrange(1, p), "random"
        else:
            x, kind = near_tie_x(rng), "near_tie"
        kinds[kind] += 1
        done += 1
        j, win = kparts(x)
        drops += win
        if j != exact_key(x)[1]:
            wrong += win
            bad += not win
    extra = max(0, MIN_NEAR - kinds["near_tie"])
    rng2 = random.Random(SEED + 7919)  # the guaranteed near ties come from their own stream, so rng is unchanged
    for _ in range(extra):
        x = near_tie_x(rng2)
        kinds["near_tie"] += 1
        j, win = kparts(x)
        drops += win
        if j != exact_key(x)[1]:
            wrong += win
            bad += not win
    res["class_rule"] = dict(checked=kinds, in_window_set=drops, wrong_index_in_window_set=wrong,
                             mismatches_outside_window_set=bad, near_ties_from_second_stream=extra)
    res["class_keys"] = class_keys(random.Random(SEED + 104729))   # own stream: the census is unchanged
    return res


def edge_class(rg):
    """A canonical value v* whose bits 64..79 sit at a window edge (0..16 or 2^16 - 16..2^16 - 1), some with bits 0..63
    above 2^64 - C0 - 5 (so that floor(2^192 v*/p) carries into bit 64). rg: Random. returns: field element."""
    while True:
        low = rg.choice(list(range(17)) + list(range(2 ** 16 - 16, 2 ** 16)))
        tail = rg.choice([rg.getrandbits(64), 2 ** 64 - 1 - rg.randrange(C0 + 5)])
        v = (rg.getrandbits(176) << 80 | low << 64 | tail) % p
        if v and exact_key(v)[0] == v:
            return v


def class_keys(rg):
    """Every member of a class whose KEY index is exact emits the top 176 bits of v* unless its window drops it, on random
    and window-edge classes; the path errors stay in their bounds (proof.md 3.3). rg: Random. returns: dict."""
    out = dict(classes=0, edge_classes=0, members=[0, 0, 0], window_drops=[0, 0, 0], wrong_index=0, mismatches=0,
               error_j1=[99, -99], error_j2=[99, -99])
    dl = deadline(0.04)
    while out["classes"] < MIN_CLASSES or (out["classes"] < 4 * MIN_CLASSES and time.time() < dl):
        edge = out["classes"] % 3 == 2
        v = edge_class(rg) if edge else exact_key(rg.randrange(1, p))[0]
        out["classes"] += 1
        out["edge_classes"] += edge
        for m in range(3):              # the member x with beta^m x = v*: exact index m
            x = v * pow(BETA, (3 - m) % 3, p) % p
            j = rule_j(x)
            if j != m:
                out["wrong_index"] += 1
                continue
            out["members"][j] += 1
            if j:
                a1 = ext_a(x)
                e1 = frac192(BETA * x % p) - a1
                if j == 1:
                    out["error_j1"] = [min(out["error_j1"][0], e1), max(out["error_j1"][1], e1)]
                else:
                    e2 = -((x >> 64) + a1) % 2 ** 192 - frac192(BETA * BETA % p * x % p)
                    out["error_j2"] = [min(out["error_j2"][0], e2), max(out["error_j2"][1], e2)]
            a = key_path(x, j)
            if a is None:
                out["window_drops"][j] += 1
            elif a >> 16 != v >> 80:
                out["mismatches"] += 1
    out["ok"] = (out["mismatches"] == 0 and out["wrong_index"] == 0 and out["window_drops"][0] == 0
                 and 0 <= out["error_j1"][0] and out["error_j1"][1] <= 7 and 0 <= out["error_j2"][0]
                 and out["error_j2"][1] <= 10)
    return out


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
    """Deliver one message (proof.md 2.1) by running route.s dlv on it; record its count against cert.DLV.

    prog, labels: route.s; mem: destination memory; msg: [dest, address, w0..w3]; calls: list. returns: None."""
    calls.append(("dlv", a64.a64_run(prog, labels, "dlv", dict(x0=LINK, x30=SENTINEL), LinkView(mem, msg)), cert.DLV))


class LazyTable(dict):
    """Memory in which the unwritten words of one cursor table read as its initial cursor, standing in for tinit's 2^36
    entries in the production-size runs of part A."""

    def __init__(self, lo, hi, base):
        """lo, hi: the table's address range; base: the initial cursor. returns: None."""
        super().__init__()
        self.lo, self.hi, self.base = lo, hi, base

    def __missing__(self, a):
        """An unwritten word. a: address. returns: the initial cursor inside the table (KeyError elsewhere)."""
        if self.lo <= a < self.hi:
            return self.base
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
    """Header of the batch kernels: G1..G5, g0 and g1 of G', and X's eight values when X is given.

    mem: memory; Hb: header address; X: point or None. returns: None.
    """
    if X is not None:
        for k, v in enumerate(header_values(X)):
            put(mem, Hb + 32 * k, v)
    for k in range(1, 6):
        mem[Hb + 408 + 8 * k] = GL[k]
    mem[Hb + 456], mem[Hb + 464] = GG[0], GG[1]


def batch_tables(mem, Hb, Yb, X, Ys):
    """Write X's header and the Y table. mem: memory; Hb, Yb: addresses; X: point; Ys: points. returns: None."""
    header(mem, Hb, X)
    for i, Y in enumerate(Ys):
        for k, v in enumerate(y_entry(Y)):
            put(mem, Yb + 160 * i + 32 * k, v)


def emit_reference(leaves, W2, start, end, cert):
    """Exact model of keys' emits (proof.md 3.2) in keys' order: entries from the last down, leaves c = 0..5 (type c + 1).

    leaves: six leaves per incidence; W2: index word; start, end: initial cursor and slot end of d1; cert: certificate.
    returns: (messages, final cursors {d1: cursor}, the count keys adds to x18, instructions saved by drops, slot drops).
    """
    msgs, cur, meter, saved, full = [], {}, 0, 0, 0
    for l in range(len(leaves) - 1, -1, -1):
        for c, x in enumerate(leaves[l]):
            A, j = key_model(x)
            meter += (0, cert.J1_PATH, cert.J2_PATH)[j]
            if A is None:               # a window dropped the leaf: no emit and no branch back
                saved += cert.EMIT + 1
                continue
            lo, mid, hi = limbs(A, 3)
            d1 = hi >> 28
            cc = cur.get(d1, start(d1))
            cur[d1] = cc + 32               # LDADD advances the cursor whether or not the slot is full
            if cc >= end(d1):
                full += 1
                saved += cert.EMIT - cert.EMIT_CHECK
                continue
            k2 = hi & (2 ** 28 - 1) | (mid & M36) << 28
            msgs.append([d1, cc, k2, lo, W2, (mid >> 28) + ((c + 1) << 36) + (l << 39)])
    return msgs, cur, meter, saved, full


def d4_of(k3):
    """bkt's d4 digit of a record: bits 16..63 of its K3 word (proof.md 4). k3: word. returns: integer."""
    return k3 >> 16 & M48


def fixup_extra(prog):
    """The fixups of a loaded leaf.s (proof.md 3.4). prog: from a64_load.
    returns: {index of a fixup's STADD: instructions the fixup adds, its chain and STADD after the branch}."""
    return {i: i - max(k for k in range(i) if prog[k][0] in ("b.lo", "b.hs", "cbnz")) for i, (mn, _) in enumerate(prog)
            if mn == "stadd"}


def flags_dead_after_fixups(prog, labels):
    """Whether every path from each fixup's join point in leaf.s sets the flags before reading them, or returns (kernels
    may leave the flags changed, proof.md 3.2): there a taken branch and a fixup differ only in the flags.
    prog, labels: loaded leaf.s. returns: bool."""
    reads = ("adc", "adcs", "sbc", "sbcs", "csel", "cinc", "cset", "ccmp", "ccmn", "b.eq", "b.ne", "b.lo", "b.hs")

    def target(t, pc):
        """Branch target index, as a64_run resolves it. t: label; pc: index after the branch. returns: index."""
        if t[:-1].isdigit():
            c_ = labels[t[:-1]]
            return [q for q in c_ if q < pc][-1] if t[-1] == "b" else [q for q in c_ if q >= pc][0]
        return labels[t][0]
    for i in fixup_extra(prog):
        todo, seen = [i + 1], set()
        while todo:
            j = todo.pop()
            if j in seen:
                continue
            seen.add(j)
            mn, o = prog[j]
            if mn in reads:
                return False
            if mn in ("adds", "subs", "cmp", "cmn", "tst", "ret"):
                continue
            if mn in ("b", "tbz", "tbnz", "cbz", "cbnz"):
                todo.append(target(o[-1], j + 1))
                if mn == "b":
                    continue
            todo.append(j + 1)
    return True


def leaf_batch(prog, labels, X, Ys, cert, W2=0x123456789, entries_override=None):
    """Run fwd, bwd and keys of leaf.s on one batch, each returning through the sentinel x30, counting each fixup's
    executions.

    prog, labels: kernel; X: point; Ys: points; cert: certificate; W2: index word; entries_override: six leaves per
    incidence written before the key pass (window and drop tests). returns: dict of counts and checks.
    """
    fx = fixup_extra(prog)
    vis = {k: dict.fromkeys(fx, 0) for k in ("fwd", "bwd")}
    q2 = len(Ys)
    Hb, Yb, Sb, Pb, SP = 0x100000, 0x200000, 1 << 34, 0x40000000, 0x7000000
    BASE, END = 1 << 62, (1 << 62) + (1 << 40)
    mem = LazyTable(Pb + 64, Pb + 64 + 8 * 2 ** 36, BASE)
    batch_tables(mem, Hb, Yb, X, Ys)
    regs = dict(x0=q2, x1=Yb, x3=Sb, x25=1, x26=0, x27=0, x28=0, sp=SP, x30=SENTINEL)
    regs.update(as_regs(X[1] * X[1] % p, ("x19", "x20", "x21", "x22")))
    out = {}
    nf = a64.a64_run(prog, labels, "fwd", regs, mem, out_regs=out, visits=vis["fwd"])
    words = [out["x0"]]
    acc = reg_int(out, ("x25", "x26", "x27", "x28"))
    want = 1
    for Y in Ys:
        want = want * (Y[1] * Y[1] - X[1] * X[1]) % p
    R = pow(acc, p - 2, p)
    regs = dict(x1=Yb + 160 * q2 + 32, x2=Hb, x3=Sb + 256 * (q2 - 1), sp=SP, x30=SENTINEL)
    regs.update(as_regs(R, ("x25", "x26", "x27", "x28")))
    out = {}
    nb = a64.a64_run(prog, labels, "bwd", regs, mem, out_regs=out, visits=vis["bwd"])
    words.append(out["x0"])
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
    regs = dict(x0=Pb, x2=Hb, x3=Sb + 256 * (q2 - 1), x12=END, x17=W2, x18=0, sp=SP, x30=SENTINEL)
    out = {}
    nk = a64.a64_run(prog, labels, "keys", regs, mem, out_regs=out, dev=(Pb, 48, port.store))
    exp, cur, meter, saved, _ = emit_reference(leaves, W2, lambda d: BASE, lambda d: END, cert)
    cursors_ok = all(mem[Pb + 64 + 8 * d] == c for d, c in cur.items())
    flat = [x for six in leaves for x in six]
    win = [kparts(x)[1] for x in flat]
    wrong = sum(w and rule_j(x) != exact_key(x)[1] for w, x in zip(win, flat))
    canon = sum(key_model(x)[0] is not None and key_model(x)[0] >> 16 != class_key(x) for w, x in zip(win, flat) if not w)
    nfix = {k: [sum(v for i, v in vis[k].items() if fx[i] == e) for e in (3, 4)] for k in vis}
    ext = {k: 3 * nfix[k][0] + 4 * nfix[k][1] for k in vis}
    return dict(fwd=nf, bwd=nb, keys=nk, acc_ok=acc % p == want, leaf_mismatches=slot_bad,
                fixups_3_4=nfix, fixup_extra=ext, fixup_words_ok=words == [C0 * sum(nfix[k]) for k in ("fwd", "bwd")],
                fwd_bwd_formulas=nf == cert.FWD_FIX + cert.FWD_INC * q2 + ext["fwd"]
                and nb == cert.BWD_FIX + cert.BWD_INC * q2 + ext["bwd"],
                window_drops=sum(key_model(x)[0] is None for x in flat),
                messages_ok=port.msgs == exp and cursors_ok, messages=len(port.msgs), meter=out["x18"],
                meter_ref=meter, windows=sum(win), wrong_index=wrong, keys_not_canonical=canon,
                port_fixed=out["x0"] == Pb + 64,
                keys_formula=nk == cert.KEY_FIX + cert.KEY_INC * q2 + out["x18"] - saved,
                final_R_is_one=final_R % p == 1)


def window_x(rng, form):
    """A field element in KEY's window set. rng: Random; form: multiplier of the small linear form.
    returns: x with kparts(x)[1] true."""
    while True:
        x = rng.randrange(1, 2 ** 150) * pow(form % p, -1, p) % p
        if kparts(x)[1]:
            return x


def reg_int(regs, names):
    """Integer from registers, low limb first. regs: register values; names: register names. returns: integer."""
    return sum(regs[n] << (64 * k) for k, n in enumerate(names))


def as_regs(v, names):
    """Registers holding the limbs of v, low limb first. v: integer < 2^256; names: register names. returns: dict."""
    return {n: w for n, w in zip(names, limbs(v, len(names)))}


def macro_counts(src, entry, names):
    """Macro invocations, nested included, in one kernel of leaf.s (label to first ret).

    src: leaf.s text; entry: label; names: macro names. returns: {name: count}.
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


def cube_root(v):
    """A cube root mod p, or None (p = 7 mod 9, so a cube's root is its power (p + 2)/9). v: residue.
    returns: integer or None."""
    c_ = pow(v, (p + 2) // 9, p)
    return c_ if c_ * c_ * c_ % p == v % p else None


def fixup_points(rg, X):
    """Points Y whose batch takes chosen fixups (proof.md 3.4): fwd's FS forming D_l (and the next incidence's first row
    shift of acc = D), bwd's FS of y_Y - y_X, its FADDW of y_Y + y_X, of x_Y + x_X/2 and LT's of x_X + (x_Y + C0).

    rg: Random; X: point. returns: list of five points, in that order."""
    xX, yX = X
    b, h, out = yX * yX % p, xX * pow(2, -1, p) % p, []
    T = 2 ** 64 - 1
    while len(out) < 5:
        k, e, n = len(out), rg.randrange(C0), rg.randrange(1, 2 ** 120) << 64
        if k == 0:
            D = p - rg.randrange(1, b)
            D += T - D % 2 ** 64 - rg.randrange(2 ** 20)
            y = sqrt_mod_p((D + b) % p) if p - b <= D < p and D >> 192 >= 2 ** 52 else None
            ok = y is not None and (y * y % p - b) % 2 ** 64 < C0 and y * y % p < b
        elif k == 1:
            y = yX - (rg.randrange(1, yX >> 64) << 64) + e
            ok = 0 <= y < yX and (y - yX) % 2 ** 64 < C0
        elif k == 2:
            y = 2 ** 256 - yX + n + T - e
            ok = y < p and y + yX >= 2 ** 256 and (y + yX) % 2 ** 64 >= 2 ** 64 - C0
        else:
            x = 2 ** 256 - (h if k == 3 else xX + C0) + n + T - e
            s_ = x + (h if k == 3 else xX + C0)
            ok = x < p and s_ >= 2 ** 256 and s_ % 2 ** 64 >= 2 ** 64 - C0
            y = sqrt_mod_p((x ** 3 + 7) % p) if ok else None
            ok = ok and y is not None
        if ok and k < 3:
            x = cube_root((y * y - 7) % p)
            ok = x is not None
        if ok:
            out.append((x, y))
    return out


def drop_leaf(rg, j):
    """A leaf outside KEY's window set whose index is j (1 or 2) and whose path-j window drops it (proof.md 3.3): a member
    of a class whose v* has bits 64..79 near 0. rg: Random; j: path. returns: field element."""
    while True:
        v = (rg.getrandbits(176) << 80 | rg.choice([0, 1, 2, 2 ** 16 - 1]) << 64 | rg.getrandbits(64)) % p
        if not v or exact_key(v)[0] != v:
            continue
        x = v * pow(BETA, 3 - j, p) % p
        if kparts(x) == (j, False) and key_model(x) == (None, j):
            return x


def degenerate_batch(prog, labels, X, Ys):
    """fwd and bwd on a batch with some Y_l = +-phi^t(X) (D_l = 0), R = acc^(p-2) as fe_inv returns it.

    prog, labels: leaf.s; X: point; Ys: points. returns: (acc, exact leaves among the other Y_l).
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
    net = {q: {k: d[k] - d["fixup_extra"][k] for k in ("fwd", "bwd")} for q, d in runs.items()}
    for k in ("fwd", "bwd"):
        res[k + "_per_incidence"] = sorted({(net[b][k] - net[a][k]) // (b - a) for a, b in ((1, 5), (5, 9), (2, 9))})
        res[k + "_fixed"] = sorted({net[q][k] - res[k + "_per_incidence"][0] * q for q in runs})
    res["fwd_bwd_formula"] = "fwd (bwd) = FWD_FIX + FWD_INC*q2 (BWD_...) + 3 or 4 per fixup (FIX_EXTRA)"
    res["fwd_bwd_formulas_hold"] = all(d["fwd_bwd_formulas"] and d["fixup_words_ok"] for d in runs.values())
    res["fixups"] = dict(sites=len(fixup_extra(prog)), flags_dead_at_every_join=flags_dead_after_fixups(prog, labels))
    res["keys_formula"] = ("keys = KEY_FIX + KEY_INC*q2 + (the count keys adds to x18: J1_PATH per j = 1 leaf, J2_PATH "
                           "per j = 2 leaf) - (EMIT - EMIT_CHECK per leaf dropped at a full slot) - (EMIT + 1 per leaf "
                           "dropped by its window)")
    res["keys_formula_holds"] = all(d["keys_formula"] for d in runs.values())
    res["all_correct"] = all(d["acc_ok"] and d["leaf_mismatches"] == 0 and d["messages_ok"] and d["meter"] ==
                             d["meter_ref"] and d["port_fixed"] and d["keys_not_canonical"] == 0 for d in runs.values())
    forms = [BETA, 1 + BETA, 1 - BETA, 2 + BETA, 1 + 2 * BETA]
    Ys = [PROD.mul(rng.randrange(1, r), GP) for _ in range(3)]
    over = [[window_x(rng, forms[(3 * i + c) % 5]) if (i + c) % 3 == 0 else rng.randrange(1, p) for c in range(6)]
            for i in range(3)]
    w = leaf_batch(prog, labels, X, Ys, cert, entries_override=over)
    res["window_batch"] = dict(leaves=18, windows=w["windows"], wrong_index=w["wrong_index"], messages=w["messages"],
                               messages_ok=w["messages_ok"], keys=w["keys"], keys_formula=w["keys_formula"])
    frng = random.Random(SEED + 49979687)   # own stream: the census is unchanged
    Yf = fixup_points(frng, X)
    f = leaf_batch(prog, labels, X, [Yf[0], PROD.mul(frng.randrange(1, r), GP)] + Yf[1:], cert)
    res["fixup_batch"] = dict(incidences=6, fixups_3_4=f["fixups_3_4"], leaf_mismatches=f["leaf_mismatches"],
                              fixup_words_ok=f["fixup_words_ok"], formulas=f["fwd_bwd_formulas"],
                              messages_ok=f["messages_ok"] and f["keys_formula"],
                              ok=f["leaf_mismatches"] == 0 and f["fixup_words_ok"] and f["fwd_bwd_formulas"]
                              and f["messages_ok"] and f["keys_formula"] and f["fixups_3_4"]["fwd"][0] >= 1
                              and f["fixups_3_4"]["fwd"][1] >= 1 and f["fixups_3_4"]["bwd"][1] >= 4)
    over = [[drop_leaf(frng, 1 + (i + c) % 2) if c in (0, 3) else frng.randrange(1, p) for c in range(6)]
            for i in range(3)]
    d = leaf_batch(prog, labels, X, Ys, cert, entries_override=over)
    split = [sum(key_model(x) == (None, j) for six in over for x in six) for j in range(3)]
    res["window_drop_batch"] = dict(leaves=18, planted=6, drops_by_path=split, messages=d["messages"],
                                    messages_ok=d["messages_ok"], meter=d["meter"], meter_ok=d["meter"] == d["meter_ref"],
                                    keys=d["keys"], formula=cert.KEY_FIX + cert.KEY_INC * 3 + d["meter"]
                                    - (cert.EMIT + 1) * sum(split), keys_formula=d["keys_formula"],
                                    ok=split[0] == 0 and split[1] >= 3 and split[2] >= 3
                                    and d["messages"] == 18 - sum(split) and d["messages_ok"]
                                    and d["meter"] == d["meter_ref"] and d["keys_formula"])
    src = open(os.path.join(HERE, "kernels", "leaf.s")).read()
    gp, gl = a64.a64_load(src.replace("\tb.lo 9f\n", "\tcbnz xzr,9f\n").replace("\tb.hs 9f\n", "\tcbnz xzr,9f\n"))
    g = leaf_batch(gp, gl, X, Ys[:2], cert)     # every fixup branch never taken: each fixup runs in place, as an inline chain
    tot = sum(cert.FIX_SITES.values())
    ext = sum(cert.FIX_EXTRA[m] * cert.FIX_SITES[m] for m in cert.FIX_SITES)
    res["forced_fixups"] = dict(incidences=2, fixups_3_4=g["fixups_3_4"], fixup_extra=g["fixup_extra"],
                                leaf_mismatches=g["leaf_mismatches"], formulas=g["fwd_bwd_formulas"],
                                ok=g["leaf_mismatches"] == 0 and g["fwd_bwd_formulas"] and g["fixup_words_ok"]
                                and g["messages_ok"] and sum(map(sum, g["fixups_3_4"].values())) == 2 * tot
                                and sum(g["fixup_extra"].values()) == 2 * ext)
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
    names = ("SHX", "FDX", "F5", "FS", "FADDW")
    counts = {k: macro_counts(src, k, names) for k in ("fwd", "bwd", "keys")}
    tot = {n: sum(c[n] for c in counts.values()) for n in names}
    res["site_counts"] = dict(counts, per_incidence=tot)
    res["sites_match_certificate"] = ({n: tot[n] for n in cert.SITES_PER_INCIDENCE} == cert.SITES_PER_INCIDENCE
                                      and {n: counts["fwd"][n] + counts["bwd"][n] for n in cert.FIX_SITES}
                                      == cert.FIX_SITES and sum(counts["keys"].values()) == 0)
    srng = random.Random(SEED + 15485863)   # own stream: the census is unchanged
    res["site_checks"] = sites.site_checks(sys.modules[__name__], srng, fprog, flabels)
    Yd = [PROD.mul(srng.randrange(1, r), GP) for _ in range(3)]
    phX = (BETA * X[0] % p, X[1])
    acc0, good = degenerate_batch(prog, labels, X, [Yd[0], phX, Yd[1], PROD.neg(X), Yd[2]])
    res["d_zero_batch"] = dict(acc_zero=acc0 % p == 0, exact_leaves=good, all_other_leaves_wrong=good == 0)
    res["final_R_one_and_unused"] = all(d["final_R_is_one"] and d["messages_ok"] for d in runs.values())
    prog_r, lab_r = a64.a64_load(open(os.path.join(HERE, "kernels", "route.s")).read())
    res["program_instructions"] = dict(leaf=len(prog), route=len(prog_r), secp_fe=len(fprog),
                                       total=len(prog) + len(prog_r) + len(fprog))
    return res


def range_scan_model(groups, cap):
    """Path model of route.s scan_ranges (the paired scanner of 892a3df8 and b00400fa): per d4 range the first record,
    one single comparison when n - 1 is odd, the TBZ tree into the 32-pair block for ((n - 1)/2) mod 32 pairs, then the
    block once per 64 records; a pair costs 6 (LDP, two LDRs, CMP, CCMP, B.EQ), and a hit re-tests both comparisons.

    groups: the K2 of each range's records, in pointer order; cap: candidate-list entries.
    returns: (instruction count, candidate index pairs in the concatenated order).
    """
    count, pairs, offset = 1, [], 0                  # cbz
    if not groups:
        return 2, pairs

    class Full(Exception):
        """The candidate list was full at an append."""

    def append(i, j):
        """One append after its list-end test (CMP, B.HS, STP). i, j: indices in the range. returns: None."""
        nonlocal count
        count += 2
        if len(pairs) >= cap:
            raise Full
        pairs.append((offset + i, offset + j))
        count += 1
    try:
        for keys in groups:
            count += 3                               # ldr end, cmp, b.hs
            if not keys:
                count += 2                           # subs, b.ne
                continue
            m, i = len(keys) - 1, 1
            count += 5                               # first pointer and K2, sub, lsr, tbz
            if m & 1:
                count += 4                           # ldr, ldr, cmp, b.eq
                if keys[0] == keys[1]:
                    append(0, 1)
                    count += 1                       # b back
                count += 2                           # mov, mov
                i += 1
            r = (m // 2) & 31
            count += 12 + (r & 1)                    # lsr, and, five movs, the TBZ tree (and a B when r is odd)

            def pair_step():
                """The pair of comparisons (i - 1, i) and (i, i + 1). returns: None."""
                nonlocal count, i
                count += 6
                first, second = keys[i - 1] == keys[i], keys[i] == keys[i + 1]
                if first or second:
                    count += 2                       # cmp, b.ne
                    if first:
                        append(i - 1, i)
                    count += 2                       # cmp, b.ne
                    if second:
                        append(i, i + 1)
                        count += 1                   # b back
                i += 2
            for _ in range(r):
                pair_step()
            count += 3                               # mov, cmp, b.lo
            while i < len(keys):
                for _ in range(32):
                    pair_step()
                count += 3
            count += 2                               # subs, b.ne
            offset += len(keys)
    except Full:
        return count + 1, pairs                      # ret
    return count + 1, pairs


class Node:
    """One node of the 4 x 4 machine: its memory, at the same addresses on every node (proof.md 5)."""
    P1, P2, HB, SB, SP = 0x10000, 0x20000, 0x100000, 1 << 34, 0x7000000
    R1, E1, R2, E2 = 0x1000000, 0x2000000, 0x3000000, 0x4000000
    CT, C4, BK, PA, CL = 0x5000000, 0x5100000, 1 << 38, 0x6100000, 0x6200000
    STRIDE = 1 << 38                 # d3 bucket t's entries end at BK + (t + 1) * STRIDE

    def __init__(self):
        """Empty memory. returns: None."""
        self.mem = {}


def planted_x(rg, d1, d2, v=None, d4=None):
    """A member of a class whose key has digits d1, d2 and d4 < 16 (or of class v*), outside KEY's window set and not
    dropped by its path (proof.md 3.3, 4). rg: Random; d1, d2: digits; v: canonical value or None; d4: digit or None.
    returns: (x, v*)."""
    while True:
        vc = v
        if vc is None:
            dd = rg.randrange(16) if d4 is None else d4
            vc = (d1 << 220 | rg.getrandbits(28) << 192 | d2 << 156 | rg.getrandbits(28) << 128 | dd << 80
                  | rg.randrange(16, 2 ** 16 - 16) << 64 | rg.getrandbits(64))
            if exact_key(vc)[0] != vc:
                continue
        j = rg.randrange(3)
        x = vc * pow(BETA, (3 - j) % 3, p) % p
        if kparts(x) == (j, False) and key_model(x)[0] is not None and key_model(x)[0] >> 16 == vc >> 80:
            return x, vc


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
    """4 x 4 machine (proof.md 2.3, 4): keys, fin, h2s, fin, bkt, prefix4, pass4 and scan_ranges against an exact model of the
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
        run_k(prog, labels, "tinit", dict(x0=side, x1=nd.P1 + 64, x2=nd.R1 + i * st1), nd, calls, tf)
        run_k(prog, labels, "tinit", dict(x0=side, x1=nd.P2 + 64, x2=nd.R2 + j * st2), nd, calls, tf)
        run_k(prog, labels, "init3", dict(x0=n3, x1=nd.CT, x2=nd.BK + nd.STRIDE - bkb, x3=nd.STRIDE), nd, calls,
              cert.INIT3_BKT * n3 + 1)
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
        exp, cur, meter, saved, full = emit_reference(batches[(i, j)], W2, lambda d: nd.R1 + i * st1,
                                                      lambda d: nd.R1 + (i + 1) * st1, cert)
        port = Port()
        o = run_k(lprog, llabels, "keys", dict(x0=nd.P1, x2=nd.HB, x3=nd.SB + 256 * (q2 - 1), x12=nd.R1 + (i + 1) * st1,
                                               x17=W2, x18=0), nd, calls,
                  cert.KEY_FIX + cert.KEY_INC * q2 + meter - saved, (nd.P1, port))
        msgs_ok &= port.msgs == exp and o["x18"] == meter
        drops["hop1"] += full
        for m in port.msgs:
            link(prog, labels, nodes[(m[0], j)].mem, m, calls, cert)
        for m in exp:
            ref1[(m[0], j)][i].append(tuple(m[2:]))
        port = Port()
        end1 = nd.R1 + (i + 1) * st1
        run_k(prog, labels, "fin", dict(x0=nd.P1, x1=side, x2=nd.E1 + 32 * i, x6=end1), nd, calls, fin_f, (nd.P1, port))
        msgs_ok &= port.msgs == [[d, nd.E1 + 32 * i, min(cur.get(d, nd.R1 + i * st1), end1), 0, 0, 0]
                                 for d in range(side)]
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
        run_k(prog, labels, "h2s", dict(x0=nd.P2, x3=side, x13=nd.E1, x14=nd.R1, x15=st1, x12=nd.R2 + (j + 1) * st2),
              nd, calls,
              cert.H2S_FIX + cert.H2S_SLOT * side + cert.H2S_REC * stored + cert.H2S_DROP * (nrec - stored)
              + cert.LOOP8 * loops,
              (nd.P2, port))
        msgs_ok &= len(port.msgs) == stored
        for m in port.msgs:
            link(prog, labels, nodes[(a, m[0])].mem, m, calls, cert)
        port = Port()
        run_k(prog, labels, "fin", dict(x0=nd.P2, x1=side, x2=nd.E2 + 32 * j, x6=nd.R2 + (j + 1) * st2), nd, calls,
              fin_f, (nd.P2, port))
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
        b0 = nd.BK + nd.STRIDE - bkb
        scan_ok &= all([nd.mem[b0 + t * nd.STRIDE + 16 * k], nd.mem[b0 + t * nd.STRIDE + 16 * k + 8]] ==
                       [e[0], d4_of(e[1][1])] for t in range(n3) for k, e in enumerate(bk[t]))
        run_k(prog, labels, "prefix4", dict(x0=n4, x1=nd.C4, x2=nd.PA), nd, calls, n4 // 4 * cert.PREFIX4_PER4 + 1)
        run_k(prog, labels, "pass4", dict(x0=n3, x1=nd.CT, x2=b0, x3=nd.STRIDE, x4=nd.C4), nd, calls,
              cert.PASS4_BKT * n3 + cert.PASS4_REC * len(seq) + cert.LOOP8 * sum(len(b) // 8 for b in bk) + 1)
        order = sorted(seq, key=lambda e: d4_of(e[1][1]))
        scan_ok &= [nd.mem[nd.PA + 8 * k] for k in range(len(seq))] == [e[0] for e in order]
        scan_ok &= all(d4_of(e[1][1]) < n4 for e in seq)
        groups = [[e[1][0] for e in order if d4_of(e[1][1]) == t] for t in range(n4)]
        mc, cands = range_scan_model(groups, ncl)
        o = run_k(prog, labels, "scan_ranges", dict(x0=n4, x1=nd.PA, x2=nd.C4, x3=nd.CL, x4=nd.CL + 16 * ncl), nd, calls,
                  mc)
        scan_ok &= mc <= cert.SCAN_REC * len(seq) + cert.SCAN_RANGE * n4 + 8 * (ncl + 1) + 2
        for c in range(nd.CL, o["x3"], 16):
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
    """Forced paths: keys and h2s with every slot full, bkt with a full bucket, scan_ranges with a full list.

    prog, labels: route.s; lprog, llabels: leaf.s; rg: Random; cert: certificate. returns: dict of checks.
    """
    out, calls = {}, []
    nd = Node()
    header(nd.mem, nd.HB)
    leaves = [[planted_x(rg, rg.randrange(4), 0)[0] for _ in range(6)] for _ in range(3)]
    for l, six in enumerate(leaves):
        for c_, x in enumerate(six):
            put(nd.mem, nd.SB + 256 * l + 32 * c_, x)
    run_k(prog, labels, "tinit", dict(x0=4, x1=nd.P1 + 64, x2=nd.R1), nd, calls, cert.TINIT_ENTRY * 4 + 1)
    exp, _, meter, saved, full = emit_reference(leaves, 5, lambda d: nd.R1, lambda d: nd.R1, cert)
    port = Port()
    o = run_k(lprog, llabels, "keys", dict(x0=nd.P1, x2=nd.HB, x3=nd.SB + 256 * 2, x12=nd.R1, x17=5, x18=0), nd, calls,
              cert.KEY_FIX + cert.KEY_INC * 3 + meter - saved, (nd.P1, port))
    out["keys_every_slot_full"] = dict(count=calls[-1][1], formula=calls[-1][2], messages=len(port.msgs), dropped=full,
                                       meter_ok=o["x18"] == meter)
    nd = Node()
    recs = [[rg.getrandbits(64), rg.getrandbits(64), 7, rg.randrange(4)] for _ in range(5)]
    for k, rec in enumerate(recs):
        for w in range(4):
            nd.mem[nd.R1 + 32 * k + 8 * w] = rec[w]
    nd.mem[nd.E1] = nd.R1 + 32 * len(recs)
    run_k(prog, labels, "tinit", dict(x0=4, x1=nd.P2 + 64, x2=nd.R2), nd, calls, cert.TINIT_ENTRY * 4 + 1)
    port = Port()
    run_k(prog, labels, "h2s", dict(x0=nd.P2, x3=1, x13=nd.E1, x14=nd.R1, x15=32 * 8, x12=nd.R2), nd, calls,
          cert.H2S_FIX + cert.H2S_SLOT + cert.H2S_DROP * 5, (nd.P2, port))
    out["h2s_every_slot_full"] = dict(count=calls[-1][1], formula=calls[-1][2], messages=len(port.msgs))
    nd = Node()
    run_k(prog, labels, "init3", dict(x0=4, x1=nd.CT, x2=nd.BK + nd.STRIDE - 16, x3=nd.STRIDE), nd, calls,
          cert.INIT3_BKT * 4 + 1)
    run_k(prog, labels, "zero4", dict(x0=4, x1=nd.C4), nd, calls, cert.ZERO4_PER4 + 1)
    for k in range(3):
        rec = [4 * rg.getrandbits(40) + 1, 2 << 16 | rg.getrandbits(16), rg.getrandbits(64), rg.getrandbits(64)]
        for w in range(4):
            nd.mem[nd.R2 + 32 * k + 8 * w] = rec[w]
    nd.mem[nd.E2] = nd.R2 + 96
    run_k(prog, labels, "bkt", dict(x0=1, x13=nd.E2, x14=nd.R2, x15=32 * 8, x3=3, x4=nd.CT, x5=3, x6=nd.C4), nd, calls,
          cert.BKT_FIX + cert.BKT_SLOT + cert.BKT_REC + 2 * cert.BKT_DROP)
    out["bkt_full_bucket"] = dict(count=calls[-1][1], formula=calls[-1][2],
                                  d4_counted_once=sum(nd.mem[nd.C4 + 8 * k] for k in range(4)) == 1,
                                  entry_ok=[nd.mem[nd.BK + 2 * nd.STRIDE - 16], nd.mem[nd.BK + 2 * nd.STRIDE - 8]]
                                  == [nd.R2, 2], cursor_at_end=nd.mem[nd.CT + 8] == nd.BK + 2 * nd.STRIDE)
    nd = Node()
    keys = [(1, 1), (1, 1), (2, 2), (3, 3), (3, 3), (4, 4)]
    for k, kk in enumerate(keys):
        nd.mem[nd.R2 + 32 * k], nd.mem[nd.R2 + 32 * k + 8] = kk
        nd.mem[nd.PA + 8 * k] = nd.R2 + 32 * k
    groups = [[kk[0] for kk in keys if kk[1] == t] for t in range(8)]
    cursor = nd.PA
    for t, g in enumerate(groups):
        cursor += 8 * len(g)
        nd.mem[nd.C4 + 8 * t] = cursor
    mc, mcand = range_scan_model(groups, 1)
    o = run_k(prog, labels, "scan_ranges", dict(x0=8, x1=nd.PA, x2=nd.C4, x3=nd.CL, x4=nd.CL + 16), nd, calls, mc)
    out["scan_full_list"] = dict(count=calls[-1][1], model=mc, listed=(o["x3"] - nd.CL) // 16, model_listed=len(mcand))
    k1, h2, b3 = out["keys_every_slot_full"], out["h2s_every_slot_full"], out["bkt_full_bucket"]
    out["ok"] = (all(c == f for _, c, f in calls) and k1["messages"] == 0 and k1["dropped"] == 18 and k1["meter_ok"]
                 and h2["messages"] == 0 and b3["d4_counted_once"] and b3["entry_ok"] and b3["cursor_at_end"]
                 and out["scan_full_list"]["listed"] == 1)
    return out


def unrolled_paths(prog, labels, rg, cert):
    """h2s, bkt and pass4 on slots and buckets of 0 to 19 records (every 1-, 2- and 4-record block combination, up to two
    8-record passes), some full; outputs against exact models and counts against the certificate's formulas.

    prog, labels: route.s; rg: Random; cert: certificate. returns: dict of checks.
    """
    calls, out, ns = [], {}, list(range(20))
    rg.shuffle(ns)
    nd = Node()                     # h2s: 20 slots of 32*19 bytes; destination d takes at most 24 + 12d records
    st, end = 32 * 19, nd.R2 + 2 ** 20      # one slot end for every destination; d's cursor starts 24 + 12d below it
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
        nd.mem[nd.P2 + 64 + 8 * d] = end - 32 * (24 + 12 * d)
    curs, exp = {d: end - 32 * (24 + 12 * d) for d in range(4)}, []
    for cur in recs:
        for rec in cur:
            d = rec[3] & M36
            if curs[d] < end:
                exp.append([d, curs[d]] + rec)
            curs[d] = min(curs[d] + 32, end)    # a drop stores the end back
    port = Port()
    run_k(prog, labels, "h2s", dict(x0=nd.P2, x3=len(ns), x13=nd.E1, x14=nd.R1, x15=st, x12=end), nd, calls,
          cert.H2S_FIX + cert.H2S_SLOT * len(ns) + cert.H2S_REC * len(exp) + cert.H2S_DROP * (sum(ns) - len(exp))
          + cert.LOOP8 * sum(n // 8 for n in ns), (nd.P2, port))
    out["h2s"] = dict(records=sum(ns), sent=len(exp), messages_ok=port.msgs == exp,
                      cursors_ok=all(nd.mem[nd.P2 + 64 + 8 * d] == curs[d] for d in range(4)))
    nd = Node()                     # bkt: 20 slots; d3 bucket t holds at most 24 + 12t entries
    n3, n4 = 4, 4
    run_k(prog, labels, "init3", dict(x0=n3, x1=nd.CT, x2=nd.BK + nd.STRIDE - 16 * 64, x3=nd.STRIDE), nd, calls,
          cert.INIT3_BKT * n3 + 1)
    start = [nd.BK + (t + 1) * nd.STRIDE - 16 * (24 + 12 * t) for t in range(n3)]
    for t in range(n3):
        nd.mem[nd.CT + 8 * t] = start[t]
    run_k(prog, labels, "zero4", dict(x0=n4, x1=nd.C4), nd, calls, cert.ZERO4_PER4 + 1)
    bk, cnt4, stored = [[] for _ in range(n3)], [0] * n4, 0
    for s_, n in enumerate(ns):
        for k in range(n):
            rec = [rg.getrandbits(64), rg.randrange(n4) << 16 | rg.getrandbits(16), rg.getrandbits(64), rg.getrandbits(64)]
            a = nd.R2 + st * s_ + 32 * k
            for w in range(4):
                nd.mem[a + 8 * w] = rec[w]
            t = rec[0] & (n3 - 1)
            if len(bk[t]) < 24 + 12 * t:
                bk[t].append((a, d4_of(rec[1])))
                cnt4[d4_of(rec[1])] += 1
                stored += 1
        nd.mem[nd.E2 + 32 * s_] = nd.R2 + st * s_ + 32 * n
    run_k(prog, labels, "bkt", dict(x0=len(ns), x13=nd.E2, x14=nd.R2, x15=st, x3=n3 - 1, x4=nd.CT, x5=n4 - 1, x6=nd.C4),
          nd, calls, cert.BKT_FIX + cert.BKT_SLOT * len(ns) + cert.BKT_REC * stored + cert.BKT_DROP * (sum(ns) - stored)
          + cert.LOOP8 * sum(n // 8 for n in ns))
    ok_b = all([nd.mem[start[t] + 16 * k], nd.mem[start[t] + 16 * k + 8]] == list(e)
               for t in range(n3) for k, e in enumerate(bk[t]))
    ok_b &= all(nd.mem[nd.CT + 8 * t] == start[t] + 16 * len(bk[t]) for t in range(n3))
    ok_b &= [nd.mem[nd.C4 + 8 * v] for v in range(n4)] == cnt4
    out["bkt"] = dict(records=sum(ns), bucketed=stored, entries_cursors_and_counts_ok=ok_b)
    nd = Node()                     # pass4: 20 buckets of 0 to 19 entries
    ents, b0 = [], nd.BK + nd.STRIDE - 16 * 19
    for t, n in enumerate(ns):
        cur = [(rg.getrandbits(64), rg.randrange(n4)) for _ in range(n)]
        for k, (a, v) in enumerate(cur):
            nd.mem[b0 + nd.STRIDE * t + 16 * k], nd.mem[b0 + nd.STRIDE * t + 16 * k + 8] = a, v
        nd.mem[nd.CT + 8 * t] = b0 + nd.STRIDE * t + 16 * n
        ents.append(cur)
    tot = [sum(v == w for c in ents for _, v in c) for w in range(n4)]
    for w in range(n4):
        nd.mem[nd.C4 + 8 * w] = nd.PA + 8 * sum(tot[:w])
    run_k(prog, labels, "pass4", dict(x0=len(ns), x1=nd.CT, x2=b0, x3=nd.STRIDE, x4=nd.C4), nd, calls,
          cert.PASS4_BKT * len(ns) + cert.PASS4_REC * sum(ns) + cert.LOOP8 * sum(n // 8 for n in ns) + 1)
    order = [a for w in range(n4) for c in ents for a, v in c if v == w]
    out["pass4"] = dict(entries=sum(ns), pointer_order_ok=[nd.mem[nd.PA + 8 * k] for k in range(sum(ns))] == order)
    out["counts_equal_formulas"] = all(c == f for _, c, f in calls)
    out["scan_ranges"] = scan_paths(prog, labels, rg, cert)
    out["ok"] = (out["counts_equal_formulas"] and out["h2s"]["messages_ok"] and out["h2s"]["cursors_ok"]
                 and out["bkt"]["entries_cursors_and_counts_ok"] and out["pass4"]["pointer_order_ok"]
                 and out["h2s"]["sent"] < sum(ns) and out["bkt"]["bucketed"] < sum(ns) and out["scan_ranges"]["ok"])
    return out


def scan_paths(prog, labels, rg, cert):
    """scan_ranges alone (proof.md 4): ranges of 0 to 129 records (every remainder ((n - 1)/2) mod 32 with either parity,
    up to two 64-record blocks) with K2 from three values or random; every range all equal (both comparisons of each pair
    hit); a full list at a pair's first and at its second append; no range. Candidates and counts against
    range_scan_model, counts within the certificate's bound, and every instruction of scan_ranges must run.

    prog, labels: route.s; rg: Random; cert: certificate. returns: dict.
    """
    first = labels["scan_ranges"][0]
    visits = {i: 0 for i in range(first, len(prog))}
    runs = []

    def run(groups, cap):
        """One call on the given ranges (K2 lists) with a list of cap entries. returns: (ok, listed, count)."""
        nd = Node()
        a, cursor = nd.R2, nd.PA
        for t, g in enumerate(groups):
            for k2 in g:
                nd.mem[a], nd.mem[a + 8] = k2, rg.getrandbits(64)     # K2, and a K3 the scan never reads
                nd.mem[cursor] = a
                a, cursor = a + 32, cursor + 8
            nd.mem[nd.C4 + 8 * t] = cursor
        o = {}
        n = a64.a64_run(prog, labels, "scan_ranges", dict(x0=len(groups), x1=nd.PA, x2=nd.C4, x3=nd.CL,
                        x4=nd.CL + 16 * cap, sp=nd.SP, x30=SENTINEL), nd.mem, out_regs=o, visits=visits)
        mc, pairs = range_scan_model(groups, cap)
        got = [(nd.mem[c], nd.mem[c + 8]) for c in range(nd.CL, o["x3"], 16)]
        want = [(nd.R2 + 32 * i, nd.R2 + 32 * j) for i, j in pairs]
        bound = cert.SCAN_REC * sum(map(len, groups)) + cert.SCAN_RANGE * len(groups) + 8 * (len(pairs) + 1) + 2
        runs.append((n == mc and got == want and n <= bound, len(pairs), n, sum(map(len, groups))))
    ns = list(range(130))
    rg.shuffle(ns)
    run([[rg.randrange(3) for _ in range(n)] for n in ns], 4096)
    run([[rg.getrandbits(64) for _ in range(n)] for n in ns], 4096)
    run([[5, 5, 7, 7, 9]], 1)                          # the second equal pair, in the next pair, finds the list full
    run([[1, 2, 2, 4, 4] + list(range(10, 70))], 1)    # 65 records: the second equal pair, in the block, ends it
    run([], 1)
    run([[0] * n for n in range(0, 130, 3)], 2 ** 14)  # every range all equal: both comparisons of every pair hit
    hit = sum(v > 0 for v in visits.values())
    return dict(calls=len(runs), records=sum(r[3] for r in runs), listed=[r[1] for r in runs], counts=[r[2] for r in runs],
                instructions=len(visits), instructions_executed=hit,
                ok=all(r[0] for r in runs) and hit == len(visits) and runs[2][1] == 1 and runs[3][1] == 1)


def load_certificate():
    """Import certificate.py with its printout captured. returns: the module."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        import certificate as cert
    return cert


def part_d(cert):
    """The certificate's verdict, claim and the count constants checked against parts A and B.

    cert: certificate module. returns: dict.
    """
    names = ("LEAVES", "FWD_INC", "FWD_FIX", "BWD_INC", "BWD_FIX", "KEY_INC", "KEY_FIX", "EMIT", "EMIT_INC",
             "EMIT_CHECK", "KROW0", "KROW", "J1_PATH", "J2_PATH", "FIX_SITES", "FIX_EXTRA", "DELTA_F", "INV", "TINIT_ENTRY", "FIN_DEST", "H2S_REC", "H2S_DROP", "H2S_SLOT", "H2S_FIX",
             "BKT_REC", "BKT_DROP", "BKT_SLOT", "BKT_FIX", "PASS4_REC", "PASS4_BKT", "SCAN_REC", "SCAN_RANGE", "SCAN_FIX", "LOOP8",
             "ZERO4_PER4", "PREFIX4_PER4", "INIT3_BKT", "DLV", "DONE")
    return dict(passed=not cert.FAILED, failed=cert.FAILED, claim="%d.%012d" % divmod(int(cert.CLAIM * 10 ** 12), 10 ** 12), q1=cert.Q1,
                mu=float(cert.T["mu"]), success_lower_bound=float(cert.T["lo"]), per_leaf=float(cert.PER_LEAF),
                d_batch=float(cert.T["d_batch"]), g_def=float(cert.T["g_def"]), d_fixup_meter=float(cert.T["d_f"]),
                d_window_drop=float(cert.T["d_win"]), d_key_meter=float(cert.T["d_w"]),
                counts={k: str(getattr(cert, k)) if not isinstance(getattr(cert, k), (int, tuple)) else getattr(cert, k)
                        for k in names})


def main():
    """Run the parts within the time budget and print one JSON object. returns: None."""
    rng = random.Random(SEED)
    cert = load_certificate()
    OUT["input"] = cfg
    OUT["budget"] = dict(seconds=BUDGET, reserve=RESERVE, minimums=dict(class_rule_random=MIN_RANDOM,
                         class_rule_near_tie=MIN_NEAR, class_keys=MIN_CLASSES, multi_node_trials=MIN_MULTI,
                         curve_trials=MIN_CURVE), note="each count runs at least its minimum (or its requested count if "
                         "smaller), then continues while its share of the remaining budget lasts")
    OUT["A_identities"] = part_a_identities(rng)
    if DO_KERNELS:
        OUT["A_kernels"] = a = part_a_kernels(rng, cert)
        OUT["B_routing"] = b = part_b(cert)
        OUT["counts_match_certificate"] = (
            a["fwd_per_incidence"] == [cert.FWD_INC] and a["bwd_per_incidence"] == [cert.BWD_INC]
            and a["fwd_fixed"] == [cert.FWD_FIX] and a["bwd_fixed"] == [cert.BWD_FIX] and a["keys_formula_holds"]
            and a["fwd_bwd_formulas_hold"] and a["fixup_batch"]["ok"] and a["window_drop_batch"]["ok"]
            and a["forced_fixups"]["ok"] and a["window_drop_batch"]["keys"] == a["window_drop_batch"]["formula"]
            and OUT["A_identities"]["class_keys"]["ok"]
            and a["site_checks"]["fixups_ok"] and a["fixups"]["flags_dead_at_every_join"]
            and a["window_batch"]["keys_formula"] and a["fe_inv"]["counts"] == [cert.INV] and b["formulas_ok"]
            and b["full_paths"]["ok"] and b["unrolled_paths"]["ok"]
            and a["sites_match_certificate"])
    OUT["C_toy"] = toy.part_c(rng, deadline, TRIALS, CURVE, MIN_MULTI, MIN_CURVE)
    OUT["D_certificate"] = part_d(cert)
    OUT["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(OUT, indent=1))


if __name__ == "__main__":
    main()
