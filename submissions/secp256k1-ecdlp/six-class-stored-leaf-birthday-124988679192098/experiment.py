"""Stored-leaf birthday search on secp256k1 (proof.md). Part A runs kernels/leaf.s and secp_fe.s on production points:
components and 3.1's identity against exact arithmetic, records against a bit-exact model of keys, each leaf's index
and label against the exact ones outside W and D' (3.3), fixups forced by planted points and a never-branching copy,
and each site type alone (sites.py); part B (machine.py) runs the emits, both hops (every message delivered by route.s
dlv or dlvk) and the grouping on a 4 x 4 machine, and the record blocks on every path; part D runs certificate.py.
Input JSON on stdin, output JSON."""
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
RESERVE = 0.5                                   # seconds kept for the output
T_END = T0 + BUDGET - RESERVE
MIN_RANDOM, MIN_NEAR = 2000, 600                # class-rule checks that always run
MIN_CLASSES = 3000                              # class-label checks that always run (a third at edges of W and D')


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
import sites  # noqa: E402
import machine  # noqa: E402

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
FB = int("e4437ed6010e88286f547fa90abfe4c3", 16)      # Frobenius pi = FA + FB*omega, N(pi) = p (proof.md 3.3)
FA = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16)
COMP_OFF = (160, 128, 96, 64, 32, 0)    # alpha', beta'', gamma of sign +, then of sign -, in a batch entry (3.1)
ERR = ((-73, 73), (-73, 206), (-212, 146))   # T(v_m) - W_m, m = 0, 1, and T(v_2) - N', over every t (3.3)
D_EDGE = ((73, 73), (206, 73), (146, 212))   # D': T(v_m) mod 2^28 below the first or within the second of 2^28
UNROLL = 2                              # incidences per pass of the fwd, bwd and keys loops
OUT, CERT_LINES = {}, []


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
            lam = 3 * P[0] * P[0] * pow(2 * P[1], -1, q) % q
        else:
            lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, q) % q
        x = (lam * lam - P[0] - Q[0]) % q
        return (x, (lam * (P[0] - x) - P[1]) % q)

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


def t128(v):
    """T(v) = floor(2^128 v/p). v: field element. returns: integer."""
    return (v << 128) // p


def comp(v):
    """A label component: floor(2^100 v/p). v: field element. returns: integer below 2^100."""
    return (v << 100) // p


def label(v):
    """The class label of canonical value v as record fields (proof.md 3.3). returns: (d1, d2, K2, K3)."""
    c1, c2 = comp(v), comp(BETA * v % p)
    return c1 & M36, c2 & M36, c1 >> 36, c2 >> 36


def class_label(x):
    """The label of x's class. x: field element. returns: (d1, d2, K2, K3)."""
    return label(exact_key(x)[0])


def neg_t(v):
    """A correction word, -T(v) mod 2^128 (proof.md 3.3). v: field element. returns: integer."""
    return -t128(v % p) % M128


def corr(v):
    """The correction words of v: -T(v) and -T(beta v) mod 2^128. v: field element. returns: tuple."""
    return neg_t(v), neg_t(BETA * v)


def u_words(c, cw=(0, 0)):
    """keys' three fractions of a stored component c < 2^256 (proof.md 3.3): U0 = c >> 128 and U1 = c6:(8s mod 2^64)
    from F1C, each plus its correction word mod 2^128, and U2 = -U0 - U1. c: integer; cw: correction words.
    returns: [U0, U1, U2]."""
    a = limbs(c)
    f1 = sum(((a[i] * GL[5 - i]) & M64) + ((a[i] * GL[4 - i]) >> 64) for i in range(4)) & M64
    words = [(a[0] * GL[4]) & M64, (a[1] * GL[3]) & M64, (a[1] * GL[2]) >> 64, (a[2] * GL[2]) & M64,
             (a[2] * GL[1]) >> 64, (a[3] * GL[1]) & M64, (a[3] * GL[0]) >> 64]
    sm = ((a[0] * (GL[3] >> 3)) >> 64) + sum(w >> 3 for w in words)
    u0 = ((c >> 128) + cw[0]) % M128
    u1 = ((((f1 + (sm >> 61)) & M64) << 64 | (sm << 3) & M64) + cw[1]) % M128
    return [u0, u1, -(u0 + u1) % M128]


def leaf_w(c3, t, ca, cb):
    """W0 and W1 of leaf t of a sign (proof.md 3.3): W_m = U_m(alpha') + U_(m+t)(beta'') + U_(m+2t)(gamma) mod 2^128.
    c3: the stored (alpha', beta'', gamma); ca, cb: the correction words of b (added to alpha') and of a (to beta'').
    returns: [W0, W1]."""
    ua, ub, ug = u_words(c3[0], ca), u_words(c3[1], cb), u_words(c3[2])
    return [(ua[m] + ub[(m + t) % 3] + ug[(m + 2 * t) % 3]) % M128 for m in (0, 1)]


def record_fields(w):
    """KEY's index on the top words F0, F1 of W0, W1, and the emitted pair (W0, W1), (W1, N') or (N', W0) with
    N' = ~(W0 + W1) (proof.md 3.2, 3.3). w: [W0, W1]. returns: (j, d1, d2, K2, K3)."""
    f0, f1 = w[0] >> 64, w[1] >> 64
    s = (f0 + f1) & M64
    j = (2 if f1 + s > M64 else 1) if f0 >= f1 else (2 if f0 + s > M64 else 0)
    n = ~(w[0] + w[1]) % M128
    a, b = ((w[0], w[1]), (w[1], n), (n, w[0]))[j]
    return j, a >> 28 & M36, b >> 28 & M36, a >> 64, b >> 64


def in_W(x):
    """The window set W of proof.md 3.3, outside which KEY's index is argmin: with X_m = 2^64 v_m/p, (W1) X_0 or X_1
    below 1 or above 2^64 - 1, (W2) X_2 below 1 or above 2^64 - 4, (W3) |X_0 - X_1| < 2, (W4) X_0 - X_2 or X_1 - X_2 in
    (-1, 4). x: leaf value in [1, p). returns: bool."""
    v = [x << 64, (BETA * x % p) << 64, (BETA * BETA % p * x % p) << 64]     # p X_m
    top = p << 64
    return (any(w < p or w > top - p for w in v[:2]) or v[2] < p or v[2] > top - 4 * p or abs(v[0] - v[1]) < 2 * p
            or -p < v[0] - v[2] < 4 * p or -p < v[1] - v[2] < 4 * p)


def in_D(x):
    """The label-error set D' of proof.md 3.3: some T(v_m) mod 2^28 within its fraction's error range of a block
    boundary. x: leaf value. returns: bool."""
    v = (x, BETA * x % p, BETA * BETA % p * x % p)
    return any(t128(w) % 2 ** 28 < lo or t128(w) % 2 ** 28 >= 2 ** 28 - hi for w, (lo, hi) in zip(v, D_EDGE))


def dft_leaf(c3, t, a, b):
    """Leaf t of a sign from its components (proof.md 3.1): (alpha' - b) + beta^t (beta'' - a) + beta^(2t) gamma.
    c3: (alpha', beta'', gamma); a, b: x_X and x_Y. returns: field element."""
    return (c3[0] - b + pow(BETA, t, p) * (c3[1] - a) + pow(BETA, 2 * t, p) * c3[2]) % p


def inverse_dft(L):
    """Components (alpha, beta', gamma) whose leaves alpha + beta^t beta' + beta^(2t) gamma are L_0, L_1, L_2
    (proof.md 3.1). L: three field elements. returns: tuple."""
    i3, b2 = pow(3, -1, p), BETA * BETA % p
    return ((L[0] + L[1] + L[2]) * i3 % p, (L[0] + b2 * L[1] + BETA * L[2]) * i3 % p,
            (L[0] + BETA * L[1] + b2 * L[2]) * i3 % p)


def represent(x, t, a, b, rg):
    """Stored components (alpha', beta'', gamma) of a sign whose leaf t is x for x_X = a and x_Y = b: two of them
    random below p or a small value plus p (as a fold output can be), the third solved mod p.
    x: leaf value; t: 0, 1, 2; a, b: field elements; rg: Random. returns: list of three integers below 2^256."""
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
    """A leaf value at an edge of D' ("D": some T(v_m) within 260 of a multiple of 2^28) or of W ("Xm": some v_m within
    6*2^192 of 0 or p; "diff": some (beta^i - beta^j) x within 6*2^192 of 0 or p). rg: Random; kind: text.
    returns: field element."""
    m = rg.randrange(3)
    if kind == "D":
        r28 = rg.choice((rg.randrange(260), 2 ** 28 - 1 - rg.randrange(260), rg.randrange(2 ** 28)))
        v = -(-(rg.getrandbits(100) << 28 | r28) * p // M128)      # the least v with that T(v)
        return v * pow(BETA, -m, p) % p or 1
    off = rg.randrange(6 << 192)
    w = off if rg.random() < 0.5 else p - 1 - off
    if kind == "Xm":
        return w * pow(BETA, -m, p) % p or 1
    return w * pow((1 - BETA, 1 - BETA * BETA, BETA - BETA * BETA)[m] % p, -1, p) % p or 1


def deadline(frac):
    """Deadline of a budget-limited loop. frac: share of the remaining budget. returns: absolute time."""
    now = time.time()
    return now + frac * max(0.0, T_END - now)


def near_tie_x(rg):
    """A leaf near an edge of W: a small value divided by one of the six forms whose small values make W (3.3).
    rg: Random. returns: field element."""
    mlt = rg.choice([1, BETA, BETA * BETA, 1 - BETA, 1 - BETA * BETA, BETA - BETA * BETA])
    dd = rg.randrange(-(1 << rg.randrange(150, 250)), 1 << rg.randrange(150, 250))
    return dd * pow(mlt % p, -1, p) % p or 1


def index_check(x, rg):
    """KEY's index from random stored components of x against the exact index. x: leaf; rg: Random.
    returns: (whether x is in W, whether the index is wrong)."""
    t, a, b = rg.randrange(3), rg.randrange(1, p), rg.randrange(1, p)
    j = record_fields(leaf_w(represent(x, t, a, b, rg), t, corr(b), corr(a)))[0]
    return in_W(x), j != exact_key(x)[1]


def part_a_identities(rng):
    """Production identities and KEY's index against the exact class rule. rng: Random. returns: dict."""
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
    rng2 = random.Random(SEED + 7919)  # the guaranteed near ties come from their own stream, so rng is unchanged
    for _ in range(extra):
        win, err = index_check(near_tie_x(rng2), rng2)
        kinds["near_tie"] += 1
        drops += win
        wrong += win and err
        bad += err and not win
    res["class_rule"] = dict(checked=kinds, in_window_set=drops, wrong_index_in_window_set=wrong,
                             mismatches_outside_window_set=bad, near_ties_from_second_stream=extra)
    res["class_labels"] = class_labels(random.Random(SEED + 104729))   # own stream
    return res


def class_labels(rg):
    """Every member of a class, given random stored components and correction words, emits the class label unless it
    lies in W or D', on random classes and classes with a member at an edge of W or D'; the computed fractions' errors
    stay in their bounds (proof.md 3.3). rg: Random. returns: dict."""
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
        for m in range(3):              # the member x with beta^m x = v*: exact index m
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
    """Memory in which the unwritten entries of a table read as their initial values, standing in for init3's and
    tinit's 2^36 entries in the production-size runs of part A (proof.md 2.3)."""

    def __init__(self, tables):
        """tables: (lo, hi, first, stride) per table, entry t reading first + t*stride. returns: None."""
        super().__init__()
        self.tables = tables

    def __missing__(self, a):
        """An unwritten word. a: address. returns: its table entry's initial value (KeyError outside the tables)."""
        for lo, hi, first, stride in self.tables:
            if lo <= a < hi:
                return first + (a - lo) // 8 * stride
        raise KeyError(a)


def rows(c):
    """The four tabled rows of c (proof.md 3.1): c * 2^(64k) mod p, k = 0..3. c: field element. returns: list."""
    return [c * pow(2, 64 * k, p) % p for k in range(4)]


def kblock(mem, Kb, a=None):
    """keys' block K (proof.md 3.2): G1 G2, G3 G4, G5 G0, G3 >> 3 G0, then -T(a) and -T(beta a) (zero when a is None).
    mem: memory; Kb: address; a: x_X or None. returns: None."""
    for k, w in enumerate((GL[1], GL[2], GL[3], GL[4], GL[5], GL[0], GL[3] >> 3, GL[0])):
        mem[Kb + 8 * k] = w
    for k, w in enumerate(corr(a) if a is not None else (0, 0)):
        put(mem, Kb + 64 + 16 * k, w, 2)


def batch_tables(mem, Hb, Kb, Yb, Sb, X, Ys):
    """X's header H (rows of 2a^3, 2a and a^2, y_X, a^4) and keys block, the Y table (y^2, y, b, b^3, 3b^2, b^4) and
    the batch entries' correction words of b (proof.md 3.1). mem: memory; Hb, Kb, Yb, Sb: addresses; X: point; Ys:
    points. returns: None."""
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
    """The six components of an incidence in entry order (proof.md 3.1): alpha' = E w, beta'' = F w and gamma = G w,
    E = b^4 + 2a^3 b, F = a^4 + 2ab^3, G = 3a^2 b^2, for w = u^2 then v^2, u = 1/(y_Y + y_X), v = 1/(y_Y - y_X).
    X, Y: points with D != 0. returns: list of six field elements."""
    a, b = X[0], Y[0]
    E, F, G = (b ** 4 + 2 * a ** 3 * b) % p, (a ** 4 + 2 * a * b ** 3) % p, 3 * a * a * b * b % p
    out = []
    for w in (pow(Y[1] + X[1], -2, p), pow(Y[1] - X[1], -2, p)):
        out += [E * w % p, F * w % p, G * w % p]
    return out


def six_leaves(X, Y):
    """Exact leaves x(phi^t(X) + Y) and x(phi^t(X) - Y), in the order 2t + s. X, Y: points.
    returns: list of six field elements."""
    out = []
    for t in range(3):
        PX = (pow(BETA, t, p) * X[0] % p, X[1])
        out += [PROD.add(PX, Y)[0], PROD.add(PX, PROD.neg(Y))[0]]
    return out


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


def leaf_batch(prog, labels, X, Ys, cert, W2=0x123456789, override=None):
    """Run fwd, bwd and keys of leaf.s on one batch (RET through the sentinel x30), counting fixup executions; check the
    components against exact arithmetic, 3.1's identity against exact affine addition, the records against the model
    of keys, and each leaf's index and label against the exact ones outside W and D'. prog, labels: kernel; X: point;
    Ys: an even number of points; cert: certificate; W2: index word; override: {(l, s): (leaf, t, rg)}, components of
    sign s of incidence l whose leaf t is the given one, stored before keys. returns: dict of counts and checks."""
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
    final_R = reg_int(out, ("x25", "x26", "x27", "x28"))   # the last R update, R_0 D_0 = 1, left unread (3.4)
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
    shift of acc = D), bwd's FS of y_Y - y_X and its FADDW of y_Y + y_X.

    rg: Random; X: point. returns: list of three points, in that order."""
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
    """fwd and bwd on a batch with some Y_l = +-phi^t(X) (D_l = 0), R = acc^(p-2) as fe_inv returns it.

    prog, labels: leaf.s; X: point; Ys: an even number of points. returns: (acc, exact components among the other
    Y_l)."""
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
    """leaf_batch with each sign's components replaced by a random representation of an edge leaf of the given kinds;
    the first four such leaves lie in W (need "W") or D' (need "D"). returns: leaf_batch's dict."""
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
    """Execute leaf.s and fe_inv on production points; derive the per-incidence and fixed counts.

    rng: Random; cert: certificate module. returns: dict.
    """
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
        and all(o[1] == "#13" for o in tbzs if o[0] == "x16")   # bwd, keys; the 36 emits test slot k by TBZ #13
    res["keys_formula"] = ("keys = KEY_FIX + KEY_INC*q2 + (the count keys adds to x18: JPATH per j = 1 or 2 leaf) + "
                           "(FLUSH per block sent, KEY_BLOCK_DROP per block dropped)")
    res["keys_formula_holds"] = all(d["keys_formula"] for d in runs.values())
    res["all_correct"] = all(d["acc_ok"] and d["component_mismatches"] == 0 and d["identity_mismatches"] == 0
                             and d["records_ok"] and d["meter"] == d["meter_ref"] and d["keys_not_canonical"] == 0
                             and d["wrong_index"] == 0 and d["sp_kept"] for d in runs.values())
    erng = random.Random(SEED + 49979687)   # own stream
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
    g = leaf_batch(gp, gl, X, Ys[:2], cert)     # every fixup branch never taken: each fixup runs in place, as an inline chain
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
    srng = random.Random(SEED + 15485863)   # own stream
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
    """The point (v, u) of {(v, u): u = beta v mod p} found by rounding (tv, tu) in the basis (B, A - B), (A - B, -A)
    of proof.md 3.3 (determinant -p), within about 2^128 of it in each coordinate. returns: (v, u)."""
    c1 = (2 * (FA * tv + (FA - FB) * tu) + p) // (2 * p)
    c2 = (2 * ((FA - FB) * tv - FB * tu) + p) // (2 * p)
    return c1 * FB + c2 * (FA - FB), c1 * (FA - FB) - c2 * FA


def planted_x(rg, d1, d2, v=None, j=None):
    """A member of a class whose label has digits d1 and d2, d3 < 4 and d4 < 16 (or of class v*), outside W and D': the
    class's v* is the lattice point nearest the quarter points of the two blocks of 2^100 v*/p and 2^100 (beta v* mod
    p)/p it needs, so all three of its fractions lie far from block boundaries.
    rg: Random; d1, d2: digits; v: canonical value or None; j: the member's index, or None for a random one.
    returns: (x, v*)."""
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
    """Import certificate.py with its printout captured. returns: (the module, its printed check lines)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        import certificate as cert
    lines = [" ".join(ln.split()) for ln in buf.getvalue().splitlines()]
    return cert, [ln for ln in lines if " ok" in ln or "FAIL" in ln or ln.startswith(("q1 ", "R_total", "declared:"))]


def part_d(cert):
    """The certificate's verdict, claim and the count constants checked against parts A and B.

    cert: certificate module. returns: dict.
    """
    names = ("LEAVES", "FWD_INC", "FWD_FIX", "BWD_INC", "BWD_FIX", "KEY_INC", "KEY_FIX", "EMIT", "EMIT_INC",
             "F1_COL", "COL4", "JPATH", "FIX_SITES", "FIX_EXTRA", "DELTA_F", "INV", "TINIT_ENTRY", "K", "FLUSH",
             "KEY_BLOCK_DROP", "H2S_BLOCK_DROP", "FIN_DEST", "FIN_SENT", "FIN_HELD", "FIN_FIX", "H2S_REC", "H2S_SLOT",
             "H2S_FIX", "BKT_REC", "BKT_DROP", "BKT_SLOT", "BKT_FIX", "PASS4_REC", "PASS4_BKT", "PASS4_FIX", "SCAN_REC",
             "SCAN_RANGE", "SCAN_FIX", "LOOP8", "ZERO4_PER4", "PREFIX4_PER4", "INIT3_BKT", "DLV", "DLVK", "DONE")
    return dict(passed=not cert.FAILED, failed=cert.FAILED, checks=CERT_LINES, claim="%d.%012d" % divmod(int(cert.CLAIM * 10 ** 12), 10 ** 12), q1=cert.Q1,
                mu=float(cert.T["mu"]), success_lower_bound=float(cert.T["lo"]), per_leaf=float(cert.PER_LEAF),
                d_batch=float(cert.T["d_batch"]), g_def=float(cert.T["g_def"]), d_fixup_meter=float(cert.T["d_f"]),
                d_label_error=float(cert.T["d_lab"]), d_key_meter=float(cert.T["d_j"]),
                counts={k: str(getattr(cert, k)) if not isinstance(getattr(cert, k), (int, tuple)) else getattr(cert, k)
                        for k in names})


def main():
    """Run the parts within the time budget and print one JSON object. returns: None."""
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
            and b["block_paths"]["ok"] and b["unrolled_paths"]["ok"] and a["sites_match_certificate"]
            and a["d_zero_batch"]["acc_zero"] and a["d_zero_batch"]["all_other_components_wrong"]
            and a["final_R_one_and_unused"] and a["loop_ends_on_bit_33"])
    OUT["D_certificate"] = part_d(cert)
    OUT["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(OUT, indent=1))


if __name__ == "__main__":
    main()
