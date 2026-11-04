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
FB = int("e4437ed6010e88286f547fa90abfe4c3", 16)
FA = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16)
FRAC_OFF = ((320, 336, 288, 304, 256, 272), (64, 80, 96, 112, 128, 144))
ENTRY, ACC, Y_OFF, CORR, KEYS0 = 512, 256, (192, 224, 352, 384, 416, 448), 160, 320
BT = (1, BETA, BETA * BETA % p)
TERMS = dict(e=1, c_low=7, d=56, total=64)
ERR = ((-132, 194), (-132, 194), (-389, 264))
D_EDGE = ((194, 132), (194, 132), (264, 389))
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
def fold(z):
    return (z % 2 ** 256 + (z >> 256) * C0) % 2 ** 256
def shx_rows(c):
    R = [c]
    for _ in range(3):
        R.append((R[-1] % 2 ** 192 * 2 ** 64 + (R[-1] >> 192) * C0) % 2 ** 256)
    return R
def tfp(R, M, tm=None, sform=False):
    f = s = lo = d = 0
    for k in range(4):
        m, w = M >> 64 * k & M64, limbs(R[k])
        f += (m * w[3] & M64) + (m * w[2] >> 64)
        s += ((m * w[2] & M64) >> 3) + ((m * w[1] >> 64) >> 3)
        lo += m * w[0] + ((m * w[1] & M64) << 64)
        d += (m * w[2] & M64) % 8 + (m * w[1] >> 64) % 8
    A = ((f << 64) + 8 * s) % M128
    if tm is not None:
        Z = sum((M >> 64 * k & M64) * R[k] for k in range(4))
        e, tot = (t128(Z % p) - (Z >> 128)) % M128, (t128(Z % p) - A) % M128
        for k, v in (("e", e), ("c_low", lo >> 128), ("d", d), ("total", tot)):
            tm[k] = max(tm.get(k, 0), v)
        tm["sum_ok"] = tm.get("sum_ok", True) and tot == e + (lo >> 128) + d and s < 2 ** 64
    return A >> 64 << 64 | s if sform else A
def fracs(mem, e, s):
    v = [get(mem, e + o, 2) for o in FRAC_OFF[s]]
    v[5] = v[5] >> 64 << 64 | (v[5] << 3 & M64)
    return v
def frac_errs(mem, e, X, Y):
    cs, out = components(X, Y), []
    for s in (0, 1):
        ex = ((cs[3 * s] - Y[0]) % p, (cs[3 * s + 1] - X[0]) % p, cs[3 * s + 2])
        out += [(t128(BT[i % 2] * ex[i // 2] % p) - v + 2 ** 127) % M128 - 2 ** 127 for i, v in enumerate(fracs(mem, e, s))]
    return out
def leaf_w(f6, t):
    u = [[f6[2 * i], f6[2 * i + 1], -(f6[2 * i] + f6[2 * i + 1]) % M128] for i in range(3)]
    return [(u[0][m] + u[1][(m + t) % 3] + u[2][(m + 2 * t) % 3]) % M128 for m in (0, 1)]
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
def weak(rg):
    return rg.randrange(p) if rg.random() < 0.85 else rg.randrange(C0) + p
def represent(x, t, a, b, rg, tm=None):
    w, c, k = weak(rg) or 1, [rg.randrange(p) for _ in range(3)], rg.randrange(3)
    c[k] = 0
    v = (x + b + BT[t] * a) * pow(w, -1, p) - sum(BT[t * i % 3] * c[i] for i in range(3))
    c[k] = v * pow(BT[t * k % 3], -1, p) % p
    bw = BETA * w % p
    bw += p if bw < C0 and rg.random() < 0.5 else 0
    out = []
    for ci, y in zip(c, (b, a, None)):
        R, cw = shx_rows(ci), corr(y) if y else (0, 0)
        out += [(tfp(R, w, tm) + cw[0]) % M128, (tfp(R, bw, tm) + cw[1]) % M128]
    return out, ((c[0] * w - b) % p, (c[1] * w - a) % p, c[2] * w % p)
def edge_leaf(rg, kind):
    m = rg.randrange(3)
    if kind == "D":
        r28 = rg.choice((rg.randrange(420), 2 ** 28 - 1 - rg.randrange(420), rg.randrange(2 ** 28)))
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
    j = record_fields(leaf_w(represent(x, t, a, b, rg)[0], t))[0]
    return in_W(x), j != exact_key(x)[1]
def part_a_identities(rng):
    res = {}
    res["beta_cube_root_of_1_mod_p"] = pow(BETA, 3, p) == 1 and BETA != 1
    res["lambda_cube_root_of_1_mod_r"] = pow(LAM, 3, r) == 1 and LAM != 1
    res["lambda_G_equals_beta_x"] = PROD.mul(LAM, GP) == (BETA * GX % p, GY)
    res["seven_nonresidue_no_x0_point"] = pow(7, (p - 1) // 2, p) == p - 1
    res["window_forms_nonzero"] = all(c % p for c in (1 - BETA, 1 - BETA * BETA, BETA - BETA * BETA, 3))
    res["label_constants"] = (FA * FA - FA * FB + FB * FB == p and FB * BETA % p == FA - FB
                              and (FA - FB) * BETA % p == p - FA)
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
               errors=[[99, -99] for _ in range(3)], tfp_terms={})
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
            w = leaf_w(represent(x, t, a, b, rg, out["tfp_terms"])[0], t)
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
    tm = out["tfp_terms"]
    out["ok"] = (out["mismatches"] == 0 and out["wrong_index"] == 0 and out["in_D"] > 0 and out["in_W"] > 0
                 and all(lo <= e[0] and e[1] <= hi for e, (lo, hi) in zip(out["errors"], ERR))
                 and tm["sum_ok"] and all(tm[k] <= v for k, v in TERMS.items()))
    return out
SENTINEL = 0x5E5E5E5E5E5E5E5E
A_HB, A_SB, A_SP = (1 << 41) + 0x100000, 1 << 42, (1 << 41) + 0x7000000
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
def batch_tables(mem, Hb, Sb, X, Ys):
    a, y = X
    for i, c in enumerate((2 * a ** 3, 2 * a, a * a, BETA)):
        for k, v in enumerate(rows(c % p)):
            put(mem, Hb + 128 * i + 32 * k - 528 * (i == 3), v)
    put(mem, Hb + 384, y)
    put(mem, Hb + 416, pow(a, 4, p))
    for k, w in enumerate(corr(a)):
        put(mem, Hb + 448 + 16 * k, w, 2)
    for l, (b, yy) in enumerate(Ys):
        for o, v in zip(Y_OFF, (yy * yy, yy, b, b ** 3, 3 * b * b, b ** 4)):
            put(mem, Sb + ENTRY * l + o, v % p)
        for k, w in enumerate(corr(b)):
            put(mem, Sb + ENTRY * l + CORR + 16 * k, w, 2)
def components(X, Y):
    a, b = X[0], Y[0]
    E, F, G = (b ** 4 + 2 * a ** 3 * b) % p, (a ** 4 + 2 * a * b ** 3) % p, 3 * a * a * b * b % p
    out = []
    for w in (pow(Y[1] + X[1], -2, p), pow(Y[1] - X[1], -2, p)):
        out += [E * w % p, F * w % p, G * w % p]
    return out
def bwd_model(X, Ys, acc, D, R, tm=None):
    a, yx = X
    out, prod = {}, lambda m, rr: sum((m >> 64 * k & M64) * rr[k] for k in range(4))
    for l in range(len(Ys) - 1, -1, -1):
        b, yy = Ys[l]
        RR = shx_rows(R)
        inv, R = fold(prod(acc[l], RR)), fold(prod(D[l], RR))
        II = shx_rows(inv)
        w = []
        for y2 in ((yy - yx) % p, yy + yx - (p if yy + yx >> 256 else 0)):
            z = fold(prod(y2, II))
            w.append(fold(z * z % 2 ** 256 + (z * z >> 256) * C0))
        bw = [fold(prod(v, rows(BETA))) for v in w]
        cs = (fold(prod(b, rows(2 * a ** 3 % p)) + pow(b, 4, p)), fold(prod(pow(b, 3, p), rows(2 * a % p)) + pow(a, 4, p)),
              fold(prod(3 * b * b % p, rows(a * a % p))))
        fr = {}
        for i, (c, y) in enumerate(zip(cs, (b, a, None))):
            Rc, cw = shx_rows(c), corr(y) if y else (0, 0)
            for s in (0, 1):
                fr[s, 2 * i] = (tfp(Rc, w[s], tm) + cw[0]) % M128
                fr[s, 2 * i + 1] = tfp(Rc, bw[s], tm, i == 2) if i == 2 else (tfp(Rc, bw[s], tm) + cw[1]) % M128
        out[l] = [fr[s, i] for s in (0, 1) for i in range(6)]
    return out, R
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
def leaf_batch(prog, labels, X, Ys, cert, W2=0x123456789, override=None, entry="bwd"):
    fx = fixup_extra(prog)
    vis = {k: dict.fromkeys(fx, 0) for k in ("fwd", "bwd")}
    q2, Sb, SP = len(Ys), A_SB, A_SP
    BASE, END = 1 << 62, (1 << 62) + (1 << 40)
    mem = LazyTable([(0, 8 << 36, machine.PB1 + 32, machine.BLK), (machine.RC, machine.RC + (8 << 36), BASE, 0)])
    batch_tables(mem, A_HB, Sb, X, Ys)
    regs = dict(x0=q2, x3=Sb, x25=1, x26=0, x27=0, x28=0, sp=SP, x30=SENTINEL)
    regs.update(as_regs(X[1] * X[1] % p, ("x19", "x20", "x21", "x22")))
    out = {}
    nf = a64.a64_run(prog, labels, "fwd", regs, mem, out_regs=out, visits=vis["fwd"])
    words, acc = [out["x0"]], reg_int(out, ("x25", "x26", "x27", "x28"))
    accs, Ds, R0 = [get(mem, Sb + ENTRY * l + ACC) for l in range(q2)], [get(mem, Sb + ENTRY * l) for l in range(q2)], \
        pow(acc, p - 2, p)
    regs = dict(x3=Sb + ENTRY * (q2 - 1), sp=A_HB, x30=SENTINEL)
    regs.update(as_regs(R0, ("x25", "x26", "x27", "x28")))
    out = {}
    nb = a64.a64_run(prog, labels, entry, regs, mem, out_regs=out, visits=vis["bwd"])
    words.append(out["x0"])
    want = 1
    for Y in Ys:
        want = want * (Y[1] * Y[1] - X[1] * X[1]) % p
    final_R = reg_int(out, ("x25", "x26", "x27", "x28"))
    sp_kept = out["sp"] == A_HB
    tm = {}
    model, mR = bwd_model(X, Ys, accs, Ds, R0, tm)
    bad_model, bad_bound, bad_leaf, err, leaves = int(mR != final_R), 0, 0, [[99, -99], [99, -99]], []
    for l, Y in enumerate(Ys):
        e, cs = Sb + ENTRY * l, components(X, Y)
        bad_model += sum(get(mem, e + o, 2) != w for o, w in zip(FRAC_OFF[0] + FRAC_OFF[1], model[l]))
        for i, d in enumerate(frac_errs(mem, e, X, Y)):
            g = i % 6 // 4
            err[g] = [min(err[g][0], d), max(err[g][1], d)]
            bad_bound += not g - 1 <= d <= TERMS["total"]
        six = six_leaves(X, Y)
        bad_leaf += sum(dft_leaf(cs[3 * s:3 * s + 3], t, X[0], Y[0]) != six[2 * t + s] for s in (0, 1) for t in range(3))
        leaves.append(six)
    for (l, s), (x, t, rg) in (override or {}).items():
        f6, c3 = represent(x, t, X[0], Ys[l][0], rg)
        f6[5] = f6[5] >> 64 << 64 | (f6[5] & M64) >> 3
        for o, v in zip(FRAC_OFF[s], f6):
            put(mem, Sb + ENTRY * l + o, v, 2)
        for t_ in range(3):
            leaves[l][2 * t_ + s] = (c3[0] + BT[t_] * c3[1] + BT[2 * t_ % 3] * c3[2]) % p
    port = machine.Port(mem, machine.MB1, machine.PB1)
    regs = dict(x3=Sb + ENTRY * (q2 - 1) + KEYS0, x12=END, x17=W2, x18=0, sp=SP, x30=SENTINEL)
    out = {}
    nk = a64.a64_run(prog, labels, "keys", regs, mem, out_regs=out, dev=port.dev)
    sp_kept &= out["sp"] == SP
    snd = machine.Sender(lambda d: BASE, END, False, cert)
    meter, fields = machine.keys_reference(sys.modules[__name__], mem, Sb, q2, W2, snd, cert)
    flat = [leaves[l][2 * t + s] for l in range(q2 - 1, -1, -1) for s in (0, 1) for t in range(3)]
    win = [in_W(x) for x in flat]
    dset = [in_D(x) for x in flat]
    wrong = sum(not w and f[0] != exact_key(x)[1] for w, f, x in zip(win, fields, flat))
    canon = sum(f[1:] != class_label(x) for w, d, f, x in zip(win, dset, fields, flat) if not w and not d)
    nfix = {k: [sum(v for i, v in vis[k].items() if fx[i] == e) for e in (3, 4)] for k in vis}
    ext = {k: 3 * nfix[k][0] + 4 * nfix[k][1] for k in vis}
    return dict(fwd=nf, bwd=nb, keys=nk, acc_ok=acc % p == want, fraction_mismatches=bad_model,
                fractions_outside_bounds=bad_bound, fraction_errors=err, tfp_terms=tm, identity_mismatches=bad_leaf,
                fixups_3_4=nfix, fixup_extra=ext,
                fixup_words_ok=words == [C0 * sum(nfix[k]) for k in ("fwd", "bwd")],
                fwd_bwd_formulas=nf == cert.fwd_count(q2) + ext["fwd"]
                and nb == (cert.BWD_FIX if entry == "bwd4" else cert.BWD2_FIX) + (cert.BWD_INC if entry == "bwd4" else cert.BWD2_INC) * q2 + ext["bwd"],
                records_ok=port.blocks == snd.msgs and not port.msgs and machine.state_ok(mem, snd, 0, machine.PB1),
                records=len(snd.seq), key_paths=meter, in_W=sum(win),
                in_D=sum(d and not w for w, d in zip(win, dset)), wrong_index=wrong, keys_not_canonical=canon,
                keys_formula=nk == cert.KEY_FIX + cert.KEY_INC * q2 + meter + snd.extra,
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
    Hb, Sb, SP = 0x100000, 1 << 35, 0x7000000
    batch_tables(mem, Hb, Sb, X, Ys)
    regs = dict(x0=len(Ys), x3=Sb, x25=1, x26=0, x27=0, x28=0, sp=SP, x30=SENTINEL)
    regs.update(as_regs(X[1] * X[1] % p, ("x19", "x20", "x21", "x22")))
    out = {}
    a64.a64_run(prog, labels, "fwd", regs, mem, out_regs=out)
    acc = reg_int(out, ("x25", "x26", "x27", "x28"))
    regs = dict(x3=Sb + ENTRY * (len(Ys) - 1), sp=Hb, x30=SENTINEL)
    regs.update(as_regs(pow(acc, p - 2, p), ("x25", "x26", "x27", "x28")))
    a64.a64_run(prog, labels, "bwd", regs, mem)
    good = 0
    for i, Y in enumerate(Ys):
        if (Y[1] * Y[1] - X[1] * X[1]) % p:
            good += sum(-1 <= d <= TERMS["total"] for d in frac_errs(mem, Sb + ENTRY * i, X, Y))
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
    src = open(os.path.join(HERE, "kernels", "leaf.s")).read()
    prog, labels = a64.a64_load(src)
    X = PROD.mul(rng.randrange(1, r), GP)
    runs = {}
    for q2 in (4, 8):
        Ys = [PROD.mul(rng.randrange(1, r), GP) for _ in range(q2)]
        runs[q2] = leaf_batch(prog, labels, X, Ys, cert, entry="bwd4")
    res["runs"] = runs
    net = {q: {k: d[k] - d["fixup_extra"][k] for k in ("fwd", "bwd")} for q, d in runs.items()}
    slope = {}
    for k in ("fwd", "bwd"):
        slope[k] = sorted({Fraction(net[b][k] - net[a][k], b - a) for a, b in ((4, 8),)})
        res[k + "_per_incidence"] = [float(v) for v in slope[k]]
    res["bwd_fixed"] = sorted({float(net[q]["bwd"] - slope["bwd"][0] * q) for q in runs})
    res["fwd_counts"] = {q: net[q]["fwd"] for q in runs}
    res["fwd_bwd_formula"] = "fwd = 9+97 q2+2 fl(q2/8)+(q2 mod 8); bwd4 = BWD_FIX+BWD_INC*q2; q2 multiple of 4"
    res["fwd_bwd_formulas_hold"] = all(d["fwd_bwd_formulas"] and d["fixup_words_ok"] for d in runs.values())
    res["counts_ok"] = (all(net[q]["fwd"] == cert.fwd_count(q) for q in runs) and slope["bwd"] == [cert.BWD_INC]
                        and res["bwd_fixed"] == [cert.BWD_FIX])
    res["fixups"] = dict(sites=len(fixup_extra(prog)), flags_dead_at_every_join=flags_dead_after_fixups(prog, labels))
    tbzs = [o for mn, o in prog if mn in ("tbz", "tbnz")]
    res["loop_ends_on_bit_34"] = [o for o in tbzs if o[0] == "x3"] == [["x3", "#34", "1b"]] * 3 and len(tbzs) == 39 \
        and all(o[1] == "#13" for o in tbzs if o[0] == "x16")
    res["branch_encoding_ok"] = (4 * 1366 * 4 <= 32764 and "\ttbz x3,#34,1b\n\tldp x0,x30,[sp],#16\n\tret\n" in src)
    res["keys_formula"] = ("keys = KEY_FIX + KEY_INC*q2 + (JPATH per j = 1 or 2 leaf of the bit-exact model) + "
                           "(FLUSH per block sent, KEY_BLOCK_DROP per block dropped)")
    res["keys_formula_holds"] = all(d["keys_formula"] for d in runs.values())
    res["tfp_terms"] = dict(largest={k: max(d["tfp_terms"][k] for d in runs.values()) for k in TERMS}, bounds=TERMS,
                            sums_exact=all(d["tfp_terms"]["sum_ok"] for d in runs.values()))
    res["fraction_errors"] = [[min(d["fraction_errors"][g][0] for d in runs.values()), max(d["fraction_errors"][g][1]
                              for d in runs.values())] for g in (0, 1)]
    res["all_correct"] = all(d["acc_ok"] and d["fraction_mismatches"] == 0 and d["fractions_outside_bounds"] == 0
                             and d["identity_mismatches"] == 0
                             and d["records_ok"] and d["keys_not_canonical"] == 0
                             and d["wrong_index"] == 0 and d["sp_kept"] for d in runs.values())
    erng = random.Random(SEED + 49979687)
    Ys = [PROD.mul(rng.randrange(1, r), GP) for _ in range(4)]
    for name, kinds, need in (("window_batch", ("Xm", "diff"), "W"), ("boundary_batch", ("D",), "D")):
        w = edge_batch(prog, labels, X, Ys, cert, erng, kinds, need)
        res[name] = dict(leaves=24, in_W=w["in_W"], in_D_outside_W=w["in_D"], wrong_index=w["wrong_index"],
                         keys_not_canonical=w["keys_not_canonical"], records=w["records"], records_ok=w["records_ok"],
                         key_paths=w["key_paths"], keys=w["keys"], keys_formula=w["keys_formula"],
                         ok=(w[{"W": "in_W", "D": "in_D"}[need]] >= 4 and w["wrong_index"] == 0
                             and w["keys_not_canonical"] == 0 and w["records"] == 24 and w["records_ok"]
                             and w["keys_formula"] and w["fwd_bwd_formulas"]
                             and w["fixup_words_ok"]))
    Yf = fixup_points(erng, X)
    f = leaf_batch(prog, labels, X, [Yf[0], PROD.mul(erng.randrange(1, r), GP)] + Yf[1:], cert)
    res["fixup_batch"] = dict(incidences=4, fixups_3_4=f["fixups_3_4"], fraction_mismatches=f["fraction_mismatches"],
                              fixup_words_ok=f["fixup_words_ok"], formulas=f["fwd_bwd_formulas"],
                              records_ok=f["records_ok"] and f["keys_formula"],
                              ok=f["fraction_mismatches"] == 0 and f["fractions_outside_bounds"] == 0
                              and f["identity_mismatches"] == 0
                              and f["fixup_words_ok"] and f["fwd_bwd_formulas"] and f["records_ok"]
                              and f["keys_formula"] and f["fixups_3_4"]["fwd"][0] >= 1
                              and f["fixups_3_4"]["fwd"][1] >= 1 and f["fixups_3_4"]["bwd"][1] >= 2)
    gp, gl = a64.a64_load(src.replace("\tb.lo 9f\n", "\tcbnz xzr,9f\n").replace("\tb.hs 9f\n", "\tcbnz xzr,9f\n"))
    g = leaf_batch(gp, gl, X, Ys[:2], cert)
    tot = sum(cert.FIX_SITES.values())
    ext = sum(cert.FIX_EXTRA[m] * cert.FIX_SITES[m] for m in cert.FIX_SITES)
    res["forced_fixups"] = dict(incidences=2, fixups_3_4=g["fixups_3_4"], fixup_extra=g["fixup_extra"],
                                fraction_mismatches=g["fraction_mismatches"], formulas=g["fwd_bwd_formulas"],
                                ok=g["fraction_mismatches"] == 0 and g["fractions_outside_bounds"] == 0
                                and g["identity_mismatches"] == 0
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
    names = ("SHX", "SHXT", "FDX", "F5", "FS", "FADDW")
    counts = {k: macro_counts(src, k, names + ("FWDI", "BWDI")) for k in ("fwd", "bwd", "bwd4", "keys")}
    copies = {"fwd": counts["fwd"]["FWDI"], "bwd": counts["bwd"]["BWDI"], "bwd4": counts["bwd4"]["BWDI"]}
    per = {n: counts["fwd"][n] // copies["fwd"] + counts["bwd4"][n] // copies["bwd4"] for n in names}
    per["SHX"] += per.pop("SHXT")
    res["site_counts"] = dict(counts, per_incidence=per, incidence_copies_in_listing=copies)
    res["sites_match_certificate"] = (all(counts[k][n] % copies[k] == 0 for k in copies for n in names)
                                      and {n: per[n] for n in cert.SITES_PER_INCIDENCE} == cert.SITES_PER_INCIDENCE
                                      and {n: per[n] for n in cert.FIX_SITES} == cert.FIX_SITES
                                      and sum(counts["keys"].values()) == 0)
    srng = random.Random(SEED + 15485863)
    res["site_checks"] = sites.site_checks(sys.modules[__name__], srng, fprog, flabels)
    Yd = [PROD.mul(srng.randrange(1, r), GP) for _ in range(4)]
    phX = (BETA * X[0] % p, X[1])
    acc0, good = degenerate_batch(prog, labels, X, [Yd[0], phX, Yd[1], PROD.neg(X), Yd[2], Yd[3]])
    res["d_zero_batch"] = dict(acc_zero=acc0 % p == 0, fractions_within_bounds=good,
                               all_other_fractions_wrong=good == 0)
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
             "JPATH", "FIX_SITES", "FIX_EXTRA", "DELTA_F", "INV", "TINIT_ENTRY", "K", "FLUSH",
             "KEY_BLOCK_DROP", "H2S_BLOCK_DROP", "FIN_DEST", "FIN_SENT", "FIN_HELD", "FIN_FIX", "H2S_REC", "H2S_SLOT",
             "H2S_FIX", "H2S_ARR", "H2S_ARR_FIX", "BKT_REC", "BKT_DROP", "BKT_SLOT", "BKT_FIX", "PASS4_REC", "PASS4_BKT",
             "PASS4_FIX", "SCAN_REC",
             "SCAN_RANGE", "SCAN_FIX", "LOOP8", "ZERO4_PER4", "PREFIX4_PER4", "INIT3_BKT", "DLV", "DLVK", "DONE",
             "G_GROUPS", "HB", "HCHG_NEW", "HCHG_SLOT", "HCHG_FIXG", "HCHG_VISIT", "HCHG_VER", "VISIT",
             "HCHG_END", "HCHG_LS0", "HCHG_LS1", "ZW", "ZERO_PASS", "STOP_CHECK", "NBAR_ALLOW", "HOME_MAX", "HEADER")
    return dict(passed=not cert.FAILED, failed=cert.FAILED, checks=CERT_LINES,
                claim="%d.%09d" % divmod(int(cert.DECLARED["time_log2"] * 10 ** 9), 10 ** 9),
                least_15_grid="%d.%015d" % divmod(int(cert.CLAIM_E * 10 ** 15), 10 ** 15), q1=cert.Q1GD,
                per_leaf_declared=format(float(cert.DECLARED["per_leaf"]), ".4f"),
                base_claim="%d.%012d" % divmod(int(cert.CLAIM * 10 ** 12), 10 ** 12),
                mu=float(cert.TFULL["mu"]), success_lower_bound=float(cert.TFULL["lo"]), per_leaf=float(cert.PER_LEAF_E),
                groups=cert.G_GROUPS, phi_16=float(cert.PHI[16]),
                d_batch=float(cert.T["d_batch"]), g_def=float(cert.T["g_def"]), d_fixup_meter=float(cert.T["d_f"]),
                d_label_error=float(cert.T["d_lab"]),
                counts={k: str(getattr(cert, k)) if not isinstance(getattr(cert, k), (int, tuple)) else getattr(cert, k)
                        for k in names})
def main():
    global CERT_LINES
    rng = random.Random(SEED)
    cert, CERT_LINES = load_certificate()
    OUT["input"] = cfg
    OUT["budget"] = dict(seconds=BUDGET, reserve=RESERVE, minimums=dict(class_rule_random=MIN_RANDOM,
                         class_rule_near_tie=MIN_NEAR, class_labels=MIN_CLASSES), note="minimums first, then shares of the remaining budget up to the caps")
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
            and a["tfp_terms"]["sums_exact"]
            and all(a["tfp_terms"]["largest"][k] <= v for k, v in TERMS.items())
            and a["d_zero_batch"]["acc_zero"] and a["d_zero_batch"]["all_other_fractions_wrong"]
            and a["final_R_one_and_unused"] and a["loop_ends_on_bit_34"] and a["branch_encoding_ok"])
    OUT["D_certificate"] = part_d(cert)
    OUT["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(OUT, indent=1))
if __name__ == "__main__":
    main()
