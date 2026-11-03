"""Stored-leaf birthday search on secp256k1 (proof.md): exact kernel counts, grouped paths, success checks."""
import contextlib
import io
import json
import os
import random
import sys
import time
from fractions import Fraction

T0 = time.time()
cfg = json.load(sys.stdin)
SEED = cfg.get("seed", 1)
BUDGET = cfg.get("seconds", 15)
DO_KERNELS = cfg.get("kernels", True)
RESERVE = 0.5
T_END = T0 + BUDGET - RESERVE
MIN_RANDOM, MIN_NEAR = 2000, 600
MIN_CLASSES = 3000


def package_root():
    here = os.path.dirname(os.path.abspath(__file__))
    for base in (here, os.path.join(here, "files"), "/snapshot/files"):
        if os.path.isfile(os.path.join(base, "kernels", "a64.py")):
            return base
    raise FileNotFoundError("kernels/a64.py is not beside " + here + " nor under files/ or /snapshot/files")


HERE = package_root()
sys.path.insert(0, os.path.join(HERE, "kernels"))
import a64
sys.path.insert(0, HERE)
import sites
import machine

M64 = 2 ** 64 - 1
M36 = 2 ** 36 - 1
M128 = 2 ** 128
p = 2 ** 256 - 2 ** 32 - 977
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
GX = int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16)
GY = int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16)
BETA = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)
LAM = int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72", 16)
C0 = 2 ** 32 + 977
GQ = (BETA << 384) // p
GL = [(GQ >> (64 * k)) & M64 for k in range(6)]
FB = int("e4437ed6010e88286f547fa90abfe4c3", 16)
FA = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16)
COMP_OFF = (160, 128, 96, 64, 32, 0)
ERR = ((-73, 73), (-73, 206), (-212, 146))
D_EDGE = ((73, 73), (206, 73), (146, 212))
UNROLL = 2
OUT, CERT_LINES = {}, []


class Curve:

    def __init__(self, q, order):
        self.q, self.r = q, order

    def add(self, P, Q):
        q = self.q
        if P is None:
            return Q
        if Q is None:
            return P
        if P[0] == Q[0]:
            if (P[1] + Q[1]) % q == 0:
                return None
            lam = 3 * P[0] * P[0] * pow(2 * P[1], -1, q) % q
        else:
            lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, q) % q
        x = (lam * lam - P[0] - Q[0]) % q
        return (x, (lam * (P[0] - x) - P[1]) % q)

    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.q)

    def mul(self, k, P):
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
    return [(v >> (64 * k)) & M64 for k in range(n)]


def put(mem, a, v, n=4):
    for k, w in enumerate(limbs(v, n)):
        mem[a + 8 * k] = w


def get(mem, a, n=4):
    return sum(mem[a + 8 * k] << (64 * k) for k in range(n))


def exact_key(x):
    v = [x, BETA * x % p, BETA * BETA % p * x % p]
    j = min(range(3), key=lambda t: v[t])
    return v[j], j


def t128(v):
    return (v << 128) // p


def comp(v):
    return (v << 100) // p


def label(v):
    c1, c2 = comp(v), comp(BETA * v % p)
    return c1 & M36, c2 & M36, c1 >> 36, c2 >> 36


def class_label(x):
    return label(exact_key(x)[0])


def neg_t(v):
    return -t128(v % p) % M128


def corr(v):
    return neg_t(v), neg_t(BETA * v)


def u_words(c, cw=(0, 0)):
    a = limbs(c)
    f1 = sum(((a[i] * GL[5 - i]) & M64) + ((a[i] * GL[4 - i]) >> 64) for i in range(4)) & M64
    words = [(a[0] * GL[4]) & M64, (a[1] * GL[3]) & M64, (a[1] * GL[2]) >> 64, (a[2] * GL[2]) & M64,
             (a[2] * GL[1]) >> 64, (a[3] * GL[1]) & M64, (a[3] * GL[0]) >> 64]
    sm = ((a[0] * (GL[3] >> 3)) >> 64) + sum(w >> 3 for w in words)
    u0 = ((c >> 128) + cw[0]) % M128
    u1 = ((((f1 + (sm >> 61)) & M64) << 64 | (sm << 3) & M64) + cw[1]) % M128
    return [u0, u1, -(u0 + u1) % M128]


def leaf_w(c3, t, ca, cb):
    ua, ub, ug = u_words(c3[0], ca), u_words(c3[1], cb), u_words(c3[2])
    return [(ua[m] + ub[(m + t) % 3] + ug[(m + 2 * t) % 3]) % M128 for m in (0, 1)]


def record_fields(w):
    f0, f1 = w[0] >> 64, w[1] >> 64
    s = (f0 + f1) & M64
    j = (2 if f1 + s > M64 else 1) if f0 >= f1 else (2 if f0 + s > M64 else 0)
    n = ~(w[0] + w[1]) % M128
    a, b = ((w[0], w[1]), (w[1], n), (n, w[0]))[j]
    return j, a >> 28 & M36, b >> 28 & M36, a >> 64, b >> 64


def in_W(x):
    v = [x << 64, (BETA * x % p) << 64, (BETA * BETA % p * x % p) << 64]
    top = p << 64
    return (any(w < p or w > top - p for w in v[:2]) or v[2] < p or v[2] > top - 4 * p or abs(v[0] - v[1]) < 2 * p
            or -p < v[0] - v[2] < 4 * p or -p < v[1] - v[2] < 4 * p)


def in_D(x):
    v = (x, BETA * x % p, BETA * BETA % p * x % p)
    return any(t128(w) % 2 ** 28 < lo or t128(w) % 2 ** 28 >= 2 ** 28 - hi for w, (lo, hi) in zip(v, D_EDGE))


def dft_leaf(c3, t, a, b):
    return (c3[0] - b + pow(BETA, t, p) * (c3[1] - a) + pow(BETA, 2 * t, p) * c3[2]) % p


def inverse_dft(L):
    i3, b2 = pow(3, -1, p), BETA * BETA % p
    return ((L[0] + L[1] + L[2]) * i3 % p, (L[0] + b2 * L[1] + BETA * L[2]) * i3 % p,
            (L[0] + BETA * L[1] + b2 * L[2]) * i3 % p)


def represent(x, t, a, b, rg):
    c = [rg.randrange(p) if rg.random() < 0.85 else rg.randrange(C0) + p for _ in range(3)]
    bt, b2t, k = pow(BETA, t, p), pow(BETA, 2 * t, p), rg.randrange(3)
    if k == 0:
        c[0] = (x + b - bt * (c[1] - a) - b2t * c[2]) % p
    elif k == 1:
        c[1] = (a + pow(bt, -1, p) * (x + b - c[0] - b2t * c[2])) % p
    else:
        c[2] = pow(b2t, -1, p) * (x + b - c[0] - bt * (c[1] - a)) % p
    return c


def edge_leaf(rg, kind):
    m = rg.randrange(3)
    if kind == "D":
        r28 = rg.choice((rg.randrange(260), 2 ** 28 - 1 - rg.randrange(260), rg.randrange(2 ** 28)))
        v = -(-(rg.getrandbits(100) << 28 | r28) * p // M128)
        return v * pow(BETA, -m, p) % p or 1
    off = rg.randrange(6 << 192)
    w = off if rg.random() < 0.5 else p - 1 - off
    if kind == "Xm":
        return w * pow(BETA, -m, p) % p or 1
    return w * pow((1 - BETA, 1 - BETA * BETA, BETA - BETA * BETA)[m] % p, -1, p) % p or 1


def deadline(frac):
    now = time.time()
    return now + frac * max(0.0, T_END - now)


def near_tie_x(rg):
    mlt = rg.choice([1, BETA, BETA * BETA, 1 - BETA, 1 - BETA * BETA, BETA - BETA * BETA])
    dd = rg.randrange(-(1 << rg.randrange(150, 250)), 1 << rg.randrange(150, 250))
    return dd * pow(mlt % p, -1, p) % p or 1


def index_check(x, rg):
    t, a, b = rg.randrange(3), rg.randrange(1, p), rg.randrange(1, p)
    j = record_fields(leaf_w(represent(x, t, a, b, rg), t, corr(b), corr(a)))[0]
    return in_W(x), j != exact_key(x)[1]


def part_a_identities(rng):
    res = {}
    res["beta_cube_root_of_1_mod_p"] = pow(BETA, 3, p) == 1 and BETA != 1
    res["lambda_cube_root_of_1_mod_r"] = pow(LAM, 3, r) == 1 and LAM != 1
    res["lambda_G_equals_beta_x"] = PROD.mul(LAM, GP) == (BETA * GX % p, GY)
    res["seven_nonresidue_no_x0_point"] = pow(7, (p - 1) // 2, p) == p - 1
    res["window_forms_nonzero"] = all(c % p for c in (1 - BETA, 1 - BETA * BETA, BETA - BETA * BETA, 3))
    res["label_constants"] = (FA * FA - FA * FB + FB * FB == p and FB * BETA % p == FA - FB
                              and (FA - FB) * BETA % p == p - FA and GL[3] % 8 == 5 and (GL[3] >> 3) < 2 ** 61)
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
        win, err = index_check(x, rng)
        drops += win
        wrong += win and err
        bad += err and not win
    extra = max(0, MIN_NEAR - kinds["near_tie"])
    rng2 = random.Random(SEED + 7919)
    for _ in range(extra):
        win, err = index_check(near_tie_x(rng2), rng2)
        kinds["near_tie"] += 1
        drops += win
        wrong += win and err
        bad += err and not win
    res["class_rule"] = dict(checked=kinds, in_window_set=drops, wrong_index_in_window_set=wrong,
                             mismatches_outside_window_set=bad, near_ties_from_second_stream=extra)
    res["class_labels"] = class_labels(random.Random(SEED + 104729))
    return res


def class_labels(rg):
    out = dict(classes=0, edge_classes=0, members=[0, 0, 0], in_W=0, in_D=0, wrong_index=0, mismatches=0,
               errors=[[99, -99] for _ in range(3)])
    dl = deadline(0.05)
    while out["classes"] < MIN_CLASSES or (out["classes"] < 4 * MIN_CLASSES and time.time() < dl):
        edge = out["classes"] % 3 == 2
        x0 = edge_leaf(rg, ("D", "Xm", "diff")[out["classes"] // 3 % 3]) if edge else rg.randrange(1, p)
        v = exact_key(x0)[0]
        out["classes"] += 1
        out["edge_classes"] += edge
        lab = label(v)
        for m in range(3):
            x = v * pow(BETA, (3 - m) % 3, p) % p
            t, a, b = rg.randrange(3), rg.randrange(1, p), rg.randrange(1, p)
            w = leaf_w(represent(x, t, a, b, rg), t, corr(b), corr(a))
            for k, c in enumerate((w[0], w[1], ~(w[0] + w[1]) % M128)):
                e = (t128(pow(BETA, k, p) * x % p) - c + 2 ** 127) % M128 - 2 ** 127
                out["errors"][k] = [min(out["errors"][k][0], e), max(out["errors"][k][1], e)]
            f = record_fields(w)
            if in_W(x):
                out["in_W"] += 1
                continue
            if f[0] != m:
                out["wrong_index"] += 1
                continue
            out["members"][m] += 1
            if in_D(x):
                out["in_D"] += 1
            elif f[1:] != lab:
                out["mismatches"] += 1
    out["ok"] = (out["mismatches"] == 0 and out["wrong_index"] == 0 and out["in_D"] > 0 and out["in_W"] > 0
                 and all(lo <= e[0] and e[1] <= hi for e, (lo, hi) in zip(out["errors"], ERR)))
    return out


SENTINEL = 0x5E5E5E5E5E5E5E5E
A_HB, A_KB, A_YB, A_SB, A_SP = (1 << 41) + 0x100000, (1 << 41) + 0x180000, (1 << 41) + 0x200000, 1 << 42, \
    (1 << 41) + 0x7000000


class LazyTable(dict):

    def __init__(self, tables):
        super().__init__()
        self.tables = tables

    def __missing__(self, a):
        for lo, hi, first, stride in self.tables:
            if lo <= a < hi:
                return first + (a - lo) // 8 * stride
        raise KeyError(a)


def rows(c):
    return [c * pow(2, 64 * k, p) % p for k in range(4)]


def kblock(mem, Kb, a=None):
    for k, w in enumerate((GL[1], GL[2], GL[3], GL[4], GL[5], GL[0], GL[3] >> 3, GL[0])):
        mem[Kb + 8 * k] = w
    for k, w in enumerate(corr(a) if a is not None else (0, 0)):
        put(mem, Kb + 64 + 16 * k, w, 2)


def batch_tables(mem, Hb, Kb, Yb, Sb, X, Ys):
    a, y = X
    for i, c in enumerate((2 * a ** 3, 2 * a, a * a)):
        for k, v in enumerate(rows(c % p)):
            put(mem, Hb + 128 * i + 32 * k, v)
    put(mem, Hb + 384, y)
    put(mem, Hb + 416, pow(a, 4, p))
    kblock(mem, Kb, a)
    for l, (b, yy) in enumerate(Ys):
        for k, v in enumerate((yy * yy, yy, b, b ** 3, 3 * b * b, b ** 4)):
            put(mem, Yb + 192 * l + 32 * k, v % p)
        for k, w in enumerate(corr(b)):
            put(mem, Sb + 256 * l + 192 + 16 * k, w, 2)


def components(X, Y):
    a, b = X[0], Y[0]
    E, F, G = (b ** 4 + 2 * a ** 3 * b) % p, (a ** 4 + 2 * a * b ** 3) % p, 3 * a * a * b * b % p
    out = []
    for w in (pow(Y[1] + X[1], -2, p), pow(Y[1] - X[1], -2, p)):
        out += [E * w % p, F * w % p, G * w % p]
    return out


def six_leaves(X, Y):
    out = []
    for t in range(3):
        PX = (pow(BETA, t, p) * X[0] % p, X[1])
        out += [PROD.add(PX, Y)[0], PROD.add(PX, PROD.neg(Y))[0]]
    return out


def fixup_extra(prog):
    return {i: i - max(k for k in range(i) if prog[k][0] in ("b.lo", "b.hs", "cbnz")) for i, (mn, _) in enumerate(prog)
            if mn == "stadd"}


def flags_dead_after_fixups(prog, labels):
    reads = ("adc", "adcs", "sbc", "sbcs", "csel", "cinc", "cset", "ccmp", "ccmn", "b.eq", "b.ne", "b.lo", "b.hs")

    def target(t, pc):
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


def leaf_batch(prog, labels, X, Ys, cert, W2=0x123456789, override=None):
    fx = fixup_extra(prog)
    vis = {k: dict.fromkeys(fx, 0) for k in ("fwd", "bwd")}
    q2 = len(Ys)
    Hb, Kb, Yb, Sb, SP = A_HB, A_KB, A_YB, A_SB, A_SP
    BASE, END = 1 << 62, (1 << 62) + (1 << 40)
    mem = LazyTable([(0, 8 << 36, machine.PB1 + 32, machine.BLK), (machine.RC, machine.RC + (8 << 36), BASE, 0)])
    batch_tables(mem, Hb, Kb, Yb, Sb, X, Ys)
    regs = dict(x0=q2, x1=Yb, x3=Sb, x25=1, x26=0, x27=0, x28=0, sp=SP, x30=SENTINEL)
    regs.update(as_regs(X[1] * X[1] % p, ("x19", "x20", "x21", "x22")))
    out = {}
    nf = a64.a64_run(prog, labels, "fwd", regs, mem, out_regs=out, visits=vis["fwd"])
    words = [out["x0"]]
    acc = reg_int(out, ("x25", "x26", "x27", "x28"))
    want = 1
    for Y in Ys:
        want = want * (Y[1] * Y[1] - X[1] * X[1]) % p
    regs = dict(x1=Yb + 192 * q2 + 32, x3=Sb + 256 * (q2 - 1), sp=Hb, x30=SENTINEL)
    regs.update(as_regs(pow(acc, p - 2, p), ("x25", "x26", "x27", "x28")))
    out = {}
    nb = a64.a64_run(prog, labels, "bwd", regs, mem, out_regs=out, visits=vis["bwd"])
    words.append(out["x0"])
    final_R = reg_int(out, ("x25", "x26", "x27", "x28"))
    sp_kept = out["sp"] == Hb
    bad_comp = bad_leaf = 0
    leaves = []
    for l, Y in enumerate(Ys):
        got = [get(mem, Sb + 256 * l + o) for o in COMP_OFF]
        bad_comp += sum(g % p != w for g, w in zip(got, components(X, Y)))
        six = six_leaves(X, Y)
        bad_leaf += sum(dft_leaf(got[3 * s:3 * s + 3], t, X[0], Y[0]) != six[2 * t + s] for s in (0, 1) for t in range(3))
        leaves.append(six)
    for (l, s), (x, t, rg) in (override or {}).items():
        c3 = represent(x, t, X[0], Ys[l][0], rg)
        for o, c in zip(COMP_OFF[3 * s:3 * s + 3], c3):
            put(mem, Sb + 256 * l + o, c)
        for t_ in range(3):
            leaves[l][2 * t_ + s] = dft_leaf(c3, t_, X[0], Ys[l][0])
    port = machine.Port(mem, machine.MB1, machine.PB1)
    regs = dict(x3=Sb + 256 * (q2 - 1) + 160, x12=END, x17=W2, x18=0, sp=Kb, x30=SENTINEL)
    out = {}
    nk = a64.a64_run(prog, labels, "keys", regs, mem, out_regs=out, dev=port.dev)
    sp_kept &= out["sp"] == Kb
    snd = machine.Sender(lambda d: BASE, END, False, cert)
    meter, fields = machine.keys_reference(sys.modules[__name__], mem, Sb, Kb, q2, W2, snd, cert)
    flat = [leaves[l][2 * t + s] for l in range(q2 - 1, -1, -1) for s in (0, 1) for t in range(3)]
    win = [in_W(x) for x in flat]
    dset = [in_D(x) for x in flat]
    wrong = sum(not w and f[0] != exact_key(x)[1] for w, f, x in zip(win, fields, flat))
    canon = sum(f[1:] != class_label(x) for w, d, f, x in zip(win, dset, fields, flat) if not w and not d)
    nfix = {k: [sum(v for i, v in vis[k].items() if fx[i] == e) for e in (3, 4)] for k in vis}
    ext = {k: 3 * nfix[k][0] + 4 * nfix[k][1] for k in vis}
    return dict(fwd=nf, bwd=nb, keys=nk, acc_ok=acc % p == want, component_mismatches=bad_comp,
                identity_mismatches=bad_leaf, fixups_3_4=nfix, fixup_extra=ext,
                fixup_words_ok=words == [C0 * sum(nfix[k]) for k in ("fwd", "bwd")],
                fwd_bwd_formulas=nf == cert.FWD_FIX + cert.FWD_INC * q2 + ext["fwd"]
                and nb == cert.BWD_FIX + cert.BWD_INC * q2 + ext["bwd"],
                records_ok=port.blocks == snd.msgs and not port.msgs and machine.state_ok(mem, snd, 0, machine.PB1),
                records=len(snd.seq), meter=out["x18"], meter_ref=meter, in_W=sum(win),
                in_D=sum(d and not w for w, d in zip(win, dset)), wrong_index=wrong, keys_not_canonical=canon,
                keys_formula=nk == cert.KEY_FIX + cert.KEY_INC * q2 + out["x18"] + snd.extra,
                final_R_is_one=final_R % p == 1, sp_kept=sp_kept)


def reg_int(regs, names):
    return sum(regs[n] << (64 * k) for k, n in enumerate(names))


def as_regs(v, names):
    return {n: w for n, w in zip(names, limbs(v, len(names)))}


def macro_counts(src, entry, names):
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
    s_ = pow(v, (p + 1) // 4, p)
    return s_ if s_ * s_ % p == v % p else None


def cube_root(v):
    c_ = pow(v, (p + 2) // 9, p)
    return c_ if c_ * c_ * c_ % p == v % p else None


def fixup_points(rg, X):
    yX = X[1]
    b, out = yX * yX % p, []
    T = 2 ** 64 - 1
    while len(out) < 3:
        k, e, n = len(out), rg.randrange(C0), rg.randrange(1, 2 ** 120) << 64
        if k == 0:
            D = p - rg.randrange(1, b)
            D += T - D % 2 ** 64 - rg.randrange(2 ** 20)
            y = sqrt_mod_p((D + b) % p) if p - b <= D < p and D >> 192 >= 2 ** 52 else None
            ok = y is not None and (y * y % p - b) % 2 ** 64 < C0 and y * y % p < b
        elif k == 1:
            y = yX - (rg.randrange(1, yX >> 64) << 64) + e
            ok = 0 <= y < yX and (y - yX) % 2 ** 64 < C0
        else:
            y = 2 ** 256 - yX + n + T - e
            ok = y < p and y + yX >= 2 ** 256 and (y + yX) % 2 ** 64 >= 2 ** 64 - C0
        x = cube_root((y * y - 7) % p) if ok else None
        if x is not None:
            out.append((x, y))
    return out


def degenerate_batch(prog, labels, X, Ys):
    mem = {}
    Hb, Kb, Yb, Sb, SP = 0x100000, 0x180000, 0x200000, 1 << 34, 0x7000000
    batch_tables(mem, Hb, Kb, Yb, Sb, X, Ys)
    regs = dict(x0=len(Ys), x1=Yb, x3=Sb, x25=1, x26=0, x27=0, x28=0, sp=SP, x30=SENTINEL)
    regs.update(as_regs(X[1] * X[1] % p, ("x19", "x20", "x21", "x22")))
    out = {}
    a64.a64_run(prog, labels, "fwd", regs, mem, out_regs=out)
    acc = reg_int(out, ("x25", "x26", "x27", "x28"))
    regs = dict(x1=Yb + 192 * len(Ys) + 32, x3=Sb + 256 * (len(Ys) - 1), sp=Hb, x30=SENTINEL)
    regs.update(as_regs(pow(acc, p - 2, p), ("x25", "x26", "x27", "x28")))
    a64.a64_run(prog, labels, "bwd", regs, mem)
    good = 0
    for i, Y in enumerate(Ys):
        if (Y[1] * Y[1] - X[1] * X[1]) % p:
            good += sum(get(mem, Sb + 256 * i + o) % p == c for o, c in zip(COMP_OFF, components(X, Y)))
    return acc, good


def edge_batch(prog, labels, X, Ys, cert, rg, kinds, need):
    over = {}
    for l in range(len(Ys)):
        for s in (0, 1):
            while True:
                x = edge_leaf(rg, kinds[(2 * l + s) % len(kinds)])
                if len(over) >= 4 or (in_W(x) if need == "W" else in_D(x) and not in_W(x)):
                    break
            over[(l, s)] = (x, rg.randrange(3), rg)
    return leaf_batch(prog, labels, X, Ys, cert, override=over)


def part_a_kernels(rng, cert):
    res = {}
    prog, labels = a64.a64_load(open(os.path.join(HERE, "kernels", "leaf.s")).read())
    X = PROD.mul(rng.randrange(1, r), GP)
    runs = {}
    for q2 in (2, 4, 6, 10):
        Ys = [PROD.mul(rng.randrange(1, r), GP) for _ in range(q2)]
        runs[q2] = leaf_batch(prog, labels, X, Ys, cert)
    res["runs"] = runs
    net = {q: {k: d[k] - d["fixup_extra"][k] for k in ("fwd", "bwd")} for q, d in runs.items()}
    slope = {}
    for k in ("fwd", "bwd"):
        slope[k] = sorted({Fraction(net[b][k] - net[a][k], b - a) for a, b in ((2, 6), (6, 10), (4, 10))})
        res[k + "_per_incidence"] = [float(v) for v in slope[k]]
        res[k + "_fixed"] = sorted({float(net[q][k] - slope[k][0] * q) for q in runs})
    res["fwd_bwd_formula"] = "fwd (bwd) = FWD_FIX + FWD_INC*q2 (BWD_...) + 3 or 4 per fixup (FIX_EXTRA), q2 even"
    res["fwd_bwd_formulas_hold"] = all(d["fwd_bwd_formulas"] and d["fixup_words_ok"] for d in runs.values())
    res["counts_ok"] = (slope["fwd"] == [cert.FWD_INC] and slope["bwd"] == [cert.BWD_INC]
                        and res["fwd_fixed"] == [cert.FWD_FIX] and res["bwd_fixed"] == [cert.BWD_FIX])
    res["fixups"] = dict(sites=len(fixup_extra(prog)), flags_dead_at_every_join=flags_dead_after_fixups(prog, labels))
    tbzs = [o for mn, o in prog if mn in ("tbz", "tbnz")]
    res["loop_ends_on_bit_33"] = [o for o in tbzs if o[0] == "x3"] == [["x3", "#33", "1b"]] * 2 and len(tbzs) == 38 \
        and all(o[1] == "#13" for o in tbzs if o[0] == "x16")
    res["keys_formula"] = ("keys = KEY_FIX + KEY_INC*q2 + (the count keys adds to x18: JPATH per j = 1 or 2 leaf) + "
                           "(FLUSH per block sent, KEY_BLOCK_DROP per block dropped)")
    res["keys_formula_holds"] = all(d["keys_formula"] for d in runs.values())
    res["all_correct"] = all(d["acc_ok"] and d["component_mismatches"] == 0 and d["identity_mismatches"] == 0
                             and d["records_ok"] and d["meter"] == d["meter_ref"] and d["keys_not_canonical"] == 0
                             and d["wrong_index"] == 0 and d["sp_kept"] for d in runs.values())
    erng = random.Random(SEED + 49979687)
    Ys = [PROD.mul(rng.randrange(1, r), GP) for _ in range(4)]
    for name, kinds, need in (("window_batch", ("Xm", "diff"), "W"), ("boundary_batch", ("D",), "D")):
        w = edge_batch(prog, labels, X, Ys, cert, erng, kinds, need)
        res[name] = dict(leaves=24, in_W=w["in_W"], in_D_outside_W=w["in_D"], wrong_index=w["wrong_index"],
                         keys_not_canonical=w["keys_not_canonical"], records=w["records"], records_ok=w["records_ok"],
                         meter_ok=w["meter"] == w["meter_ref"], keys=w["keys"], keys_formula=w["keys_formula"],
                         ok=(w[{"W": "in_W", "D": "in_D"}[need]] >= 4 and w["wrong_index"] == 0
                             and w["keys_not_canonical"] == 0 and w["records"] == 24 and w["records_ok"]
                             and w["meter"] == w["meter_ref"] and w["keys_formula"]))
    Yf = fixup_points(erng, X)
    f = leaf_batch(prog, labels, X, [Yf[0], PROD.mul(erng.randrange(1, r), GP)] + Yf[1:], cert)
    res["fixup_batch"] = dict(incidences=4, fixups_3_4=f["fixups_3_4"], component_mismatches=f["component_mismatches"],
                              fixup_words_ok=f["fixup_words_ok"], formulas=f["fwd_bwd_formulas"],
                              records_ok=f["records_ok"] and f["keys_formula"],
                              ok=f["component_mismatches"] == 0 and f["identity_mismatches"] == 0
                              and f["fixup_words_ok"] and f["fwd_bwd_formulas"] and f["records_ok"]
                              and f["keys_formula"] and f["fixups_3_4"]["fwd"][0] >= 1
                              and f["fixups_3_4"]["fwd"][1] >= 1 and f["fixups_3_4"]["bwd"][1] >= 2)
    src = open(os.path.join(HERE, "kernels", "leaf.s")).read()
    gp, gl = a64.a64_load(src.replace("\tb.lo 9f\n", "\tcbnz xzr,9f\n").replace("\tb.hs 9f\n", "\tcbnz xzr,9f\n"))
    g = leaf_batch(gp, gl, X, Ys[:2], cert)
    tot = sum(cert.FIX_SITES.values())
    ext = sum(cert.FIX_EXTRA[m] * cert.FIX_SITES[m] for m in cert.FIX_SITES)
    res["forced_fixups"] = dict(incidences=2, fixups_3_4=g["fixups_3_4"], fixup_extra=g["fixup_extra"],
                                component_mismatches=g["component_mismatches"], formulas=g["fwd_bwd_formulas"],
                                ok=g["component_mismatches"] == 0 and g["identity_mismatches"] == 0
                                and g["fwd_bwd_formulas"] and g["fixup_words_ok"] and g["records_ok"]
                                and sum(map(sum, g["fixups_3_4"].values())) == 2 * tot
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
    names = ("SHX", "FDX", "F5", "FS", "FADDW")
    counts = {k: macro_counts(src, k, names) for k in ("fwd", "bwd", "keys")}
    per = {n: sum(c[n] for c in counts.values()) // UNROLL for n in names}
    res["site_counts"] = dict(counts, per_incidence=per, incidences_per_loop_pass=UNROLL)
    res["sites_match_certificate"] = (all(sum(c[n] for c in counts.values()) % UNROLL == 0 for n in names)
                                      and {n: per[n] for n in cert.SITES_PER_INCIDENCE} == cert.SITES_PER_INCIDENCE
                                      and {n: per[n] for n in cert.FIX_SITES} == cert.FIX_SITES
                                      and sum(counts["keys"].values()) == 0)
    srng = random.Random(SEED + 15485863)
    res["site_checks"] = sites.site_checks(sys.modules[__name__], srng, fprog, flabels)
    Yd = [PROD.mul(srng.randrange(1, r), GP) for _ in range(4)]
    phX = (BETA * X[0] % p, X[1])
    acc0, good = degenerate_batch(prog, labels, X, [Yd[0], phX, Yd[1], PROD.neg(X), Yd[2], Yd[3]])
    res["d_zero_batch"] = dict(acc_zero=acc0 % p == 0, exact_components=good, all_other_components_wrong=good == 0)
    res["final_R_one_and_unused"] = all(d["final_R_is_one"] and d["records_ok"] for d in runs.values())
    prog_r, lab_r = a64.a64_load(open(os.path.join(HERE, "kernels", "route.s")).read())
    res["program_instructions"] = dict(leaf=len(prog), route=len(prog_r), secp_fe=len(fprog),
                                       total=len(prog) + len(prog_r) + len(fprog))
    return res


def lattice_point(tv, tu):
    c1 = (2 * (FA * tv + (FA - FB) * tu) + p) // (2 * p)
    c2 = (2 * ((FA - FB) * tv - FB * tu) + p) // (2 * p)
    return c1 * FB + c2 * (FA - FB), c1 * (FA - FB) - c2 * FA


def planted_x(rg, d1, d2, v=None, j=None):
    while True:
        vc = v
        if vc is None:
            c1 = (rg.getrandbits(40) << 24 | rg.randrange(4)) << 36 | d1
            c2 = (rg.getrandbits(16) << 48 | rg.randrange(16)) << 36 | d2
            vv, uu = lattice_point((4 * c1 + 1) * p >> 102, (4 * c2 + 1) * p >> 102)
            if not (0 < vv < p and uu == BETA * vv % p and comp(vv) == c1 and comp(uu) == c2 and exact_key(vv)[0] == vv):
                continue
            vc = vv
        m = rg.randrange(3) if j is None else j
        x = vc * pow(BETA, (3 - m) % 3, p) % p
        if exact_key(x)[1] == m and not in_W(x) and not in_D(x):
            return x, vc


def load_certificate():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        import certificate as cert
    lines = [" ".join(ln.split()) for ln in buf.getvalue().splitlines()]
    return cert, [ln for ln in lines if " ok" in ln or "FAIL" in ln or ln.startswith(("q1 ", "R_total", "declared:"))]


def part_d(cert):
    names = ("LEAVES", "FWD_INC", "FWD_FIX", "BWD_INC", "BWD_FIX", "KEY_INC", "KEY_FIX", "EMIT", "EMIT_INC",
             "F1_COL", "COL4", "JPATH", "FIX_SITES", "FIX_EXTRA", "DELTA_F", "INV", "TINIT_ENTRY", "K", "FLUSH",
             "KEY_BLOCK_DROP", "H2S_BLOCK_DROP", "FIN_DEST", "FIN_SENT", "FIN_HELD", "FIN_FIX", "H2S_REC", "H2S_SLOT",
             "H2S_FIX", "BKT_REC", "BKT_DROP", "BKT_SLOT", "BKT_FIX", "PASS4_REC", "PASS4_BKT", "PASS4_FIX", "SCAN_REC",
             "SCAN_RANGE", "SCAN_FIX", "LOOP8", "ZERO4_PER4", "PREFIX4_PER4", "INIT3_BKT", "DLV", "DLVK", "DONE",
             "G_GROUPS", "HB", "HCHG_NEW", "HCHG_SLOT", "HCHG_FIXG", "HCHG_PROBE", "HCHG_VER", "HCHG_ADV",
             "HCHG_END", "HCHG_LS0", "HCHG_LS1", "ZERO8_PER16", "ZW", "ZERO_PASS", "STOP_CHECK", "NBAR_ALLOW", "HOME_MAX")
    return dict(passed=not cert.FAILED, failed=cert.FAILED, checks=CERT_LINES, claim="%d.%012d" % divmod(int(cert.CLAIM_E * 10 ** 12), 10 ** 12), q1=cert.Q1GD,
                base_claim="%d.%012d" % divmod(int(cert.CLAIM * 10 ** 12), 10 ** 12),
                mu=float(cert.TFULL["mu"]), success_lower_bound=float(cert.TFULL["lo"]), per_leaf=float(cert.PER_LEAF_E),
                groups=cert.G_GROUPS, phi_16=float(cert.PHI[16]),
                d_batch=float(cert.T["d_batch"]), g_def=float(cert.T["g_def"]), d_fixup_meter=float(cert.T["d_f"]),
                d_label_error=float(cert.T["d_lab"]), d_key_meter=float(cert.T["d_j"]),
                counts={k: str(getattr(cert, k)) if not isinstance(getattr(cert, k), (int, tuple)) else getattr(cert, k)
                        for k in names})


def main():
    global CERT_LINES
    rng = random.Random(SEED)
    cert, CERT_LINES = load_certificate()
    OUT["input"] = cfg
    OUT["budget"] = dict(seconds=BUDGET, reserve=RESERVE, minimums=dict(class_rule_random=MIN_RANDOM,
                         class_rule_near_tie=MIN_NEAR, class_labels=MIN_CLASSES), note="each count runs at least its "
                         "minimum, then continues while its share of the remaining budget lasts, up to its cap (26,000 "
                         "class-rule inputs, 12,000 classes)")
    OUT["A_identities"] = part_a_identities(rng)
    if DO_KERNELS:
        OUT["A_kernels"] = a = part_a_kernels(rng, cert)
        OUT["B_routing"] = b = machine.part_b(sys.modules[__name__], cert)
        OUT["counts_match_certificate"] = (
            a["counts_ok"] and a["keys_formula_holds"] and a["fwd_bwd_formulas_hold"] and a["all_correct"]
            and a["fixup_batch"]["ok"] and a["window_batch"]["ok"] and a["boundary_batch"]["ok"]
            and a["forced_fixups"]["ok"] and OUT["A_identities"]["class_labels"]["ok"]
            and OUT["A_identities"]["class_rule"]["mismatches_outside_window_set"] == 0
            and a["site_checks"]["fixups_ok"] and a["site_checks"]["all_ok"] and a["fixups"]["flags_dead_at_every_join"]
            and a["fe_inv"]["counts"] == [cert.INV] and a["fe_inv"]["correct"] and b["formulas_ok"]
            and b["block_paths"]["ok"] and b["unrolled_paths"]["ok"] and b["grouped_paths"]["ok"]
            and a["sites_match_certificate"]
            and a["d_zero_batch"]["acc_zero"] and a["d_zero_batch"]["all_other_components_wrong"]
            and a["final_R_one_and_unused"] and a["loop_ends_on_bit_33"])
    # Same bounded sandbox: run cap-edge component and retain its structured result.
    import runpy
    saved_stdin = sys.stdin
    capture = io.StringIO()
    try:
        sys.stdin = io.StringIO('{}')
        with contextlib.redirect_stdout(capture):
            runpy.run_path(os.path.join(HERE, 'guard_probe.py'))
    finally:
        sys.stdin = saved_stdin
    OUT['guard_edge'] = json.loads(capture.getvalue())
    OUT['counts_match_certificate'] &= OUT['guard_edge']['all_checks']
    OUT["D_certificate"] = part_d(cert)
    OUT["D_certificate"]["breakdown"] = {"leaf":cert.lg(sum(cert.PHI[g]*cert.E_LEAF[g] for g in range(cert.G_GROUPS))), "route":cert.lg(sum(cert.PHI[g]*cert.E_ROUTE[g] for g in range(cert.G_GROUPS))), "group":cert.lg(sum(cert.PHI[g]*cert.E_GROUP[g] for g in range(cert.G_GROUPS))), "sameX":cert.lg(cert.SXV), "K2":cert.lg(cert.EK2*cert.HCHG_VER), "list":cert.lg(cert.ELIST*cert.HCHG_LS1), "defect":cert.lg(cert.DEFV+cert.DEFK2+cert.DEFL), "coins":cert.lg(cert.NODES*(cert.Q1GD+cert.Q2+1)*cert.COIN_SCALAR), "S":cert.lg(cert.SG_FULL), "worst":cert.lg(cert.R_WORST), "total":cert.lg(cert.E_TOTAL), "exchange":cert.lg(cert.NODES*cert.G_GROUPS*4*cert.EXCH), "control":cert.lg(cert.NODES*cert.PHASE_CTRL_G)}
    OUT["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(OUT, indent=1))


if __name__ == "__main__":
    main()
