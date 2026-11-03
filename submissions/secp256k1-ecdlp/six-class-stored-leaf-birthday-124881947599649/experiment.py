"""Stored-leaf birthday search on secp256k1 (proof.md 8): part A runs leaf.s and secp_fe.s on production points
against exact models (fractions, records, labels, fixups, sites.py); part B (machine.py) the 4 x 4 machine; part D
certificate.py."""
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
FB = int("e4437ed6010e88286f547fa90abfe4c3", 16)      # Frobenius pi = FA + FB*omega, N(pi) = p (proof.md 3.3)
FA = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16)
FRAC_OFF = ((160, 176, 128, 144, 96, 112), (0, 16, 32, 48, 64, 80))   # V0, V1 of alpha', beta'', gamma in a batch
                                        # entry, sign + then - (3.1); gamma's V1 is stored as (s, hi), low word 8s
BT = (1, BETA, BETA * BETA % p)
TERMS = dict(e=1, c_low=7, d=56, total=64)      # a truncated product's error terms and their sum, at most (3.3)
ERR = ((-132, 194), (-132, 194), (-389, 264))   # T(v_m) - W_m, m = 0, 1, and T(v_2) - N', over every t (3.3)
D_EDGE = ((194, 132), (194, 132), (264, 389))   # D': T(v_m) mod 2^28 below the first or within the second of 2^28
UNROLL = 4                              # incidences per pass of the fwd, bwd and keys loops
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


def fold(z):
    """FDX as executed (proof.md 3.4): (z mod 2^256) + floor(z/2^256) C0, kept mod 2^256. returns: integer."""
    return (z % 2 ** 256 + (z >> 256) * C0) % 2 ** 256


def sqf(S):
    """SQF on S = c0..c7 as executed (proof.md 3.4): c7 C0 into words 3..6 (the carry out of word 6 dropped), then
    (c4, c5, c6) C0 into words 0..3, C0 added on a carry out of word 3. returns: (output, step 1 dropped, fixup)."""
    mid = (S >> 192) % 2 ** 256 + (S >> 448) * C0
    S1 = S % 2 ** 192 + mid % 2 ** 256 * 2 ** 192
    t = S1 % 2 ** 256 + (S1 >> 256) * C0
    return t % 2 ** 256 + C0 * (t >> 256), mid >> 256 > 0, t >> 256 > 0


def shx_rows(c):
    """Rows c, SHX(c), SHX^2(c), SHX^3(c) as executed (proof.md 3.4). returns: four integers."""
    R = [c]
    for _ in range(3):
        R.append((R[-1] % 2 ** 192 * 2 ** 64 + (R[-1] >> 192) * C0) % 2 ** 256)
    return R


def tfp(R, M, tm=None, sform=False):
    """TR0/TRK's truncated product on rows R and multiplier M (3.3): A = f 2^64 + 8s mod 2^128. tm: None or a dict of
    the largest terms and their sum checks; sform: return (top word, s). returns: A."""
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
    """Sign s's six stored fractions at entry e as values (gamma's V1 from its s-form). returns: list."""
    v = [get(mem, e + o, 2) for o in FRAC_OFF[s]]
    v[5] = v[5] >> 64 << 64 | (v[5] << 3 & M64)
    return v


def frac_errs(mem, e, X, Y):
    """T - stored value, signed, of the twelve fractions at entry e (FRAC_OFF order). returns: list."""
    cs, out = components(X, Y), []
    for s in (0, 1):
        ex = ((cs[3 * s] - Y[0]) % p, (cs[3 * s + 1] - X[0]) % p, cs[3 * s + 2])
        out += [(t128(BT[i % 2] * ex[i // 2] % p) - v + 2 ** 127) % M128 - 2 ** 127 for i, v in enumerate(fracs(mem, e, s))]
    return out


def leaf_w(f6, t):
    """W0 and W1 of leaf t of a sign (proof.md 3.3): W_m = U_m(alpha) + U_(m+t)(beta') + U_(m+2t)(gamma) mod 2^128,
    U_2 = -U_0 - U_1. f6: V0, V1 of alpha, beta' and gamma. returns: [W0, W1]."""
    u = [[f6[2 * i], f6[2 * i + 1], -(f6[2 * i] + f6[2 * i + 1]) % M128] for i in range(3)]
    return [(u[0][m] + u[1][(m + t) % 3] + u[2][(m + 2 * t) % 3]) % M128 for m in (0, 1)]


def record_fields(w):
    """KEY's index on the top words of W0, W1 and the emitted pair (W0, W1), (W1, N') or (N', W0), N' = ~(W0 + W1)
    (proof.md 3.2, 3.3). w: [W0, W1]. returns: (j, d1, d2, K2, K3)."""
    f0, f1 = w[0] >> 64, w[1] >> 64
    s = (f0 + f1) & M64
    j = (2 if f1 + s > M64 else 1) if f0 >= f1 else (2 if f0 + s > M64 else 0)
    n = ~(w[0] + w[1]) % M128
    a, b = ((w[0], w[1]), (w[1], n), (n, w[0]))[j]
    return j, a >> 28 & M36, b >> 28 & M36, a >> 64, b >> 64


def in_W(x):
    """The window set W of proof.md 3.3, outside which KEY's index is argmin: with X_m = 2^64 v_m/p, (W1) X_0 or X_1
    below 1 or above 2^64 - 1, (W2) X_2 below 1 or above 2^64 - 4, (W3) |X_0 - X_1| < 2, (W4) X_0 - X_2 or X_1 - X_2
    in (-1, 4). x: leaf value in [1, p). returns: bool."""
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
    """Leaf t of a sign (proof.md 3.1): (alpha' - b) + beta^t (beta'' - a) + beta^(2t) gamma. c3: (alpha', beta'',
    gamma); a, b: x_X, x_Y. returns: field element."""
    return (c3[0] - b + pow(BETA, t, p) * (c3[1] - a) + pow(BETA, 2 * t, p) * c3[2]) % p


def inverse_dft(L):
    """Components (alpha, beta', gamma) whose leaves alpha + beta^t beta' + beta^(2t) gamma are L_0, L_1, L_2
    (proof.md 3.1). L: three field elements. returns: tuple."""
    i3, b2 = pow(3, -1, p), BETA * BETA % p
    return ((L[0] + L[1] + L[2]) * i3 % p, (L[0] + b2 * L[1] + BETA * L[2]) * i3 % p,
            (L[0] + BETA * L[1] + b2 * L[2]) * i3 % p)


def weak(rg):
    """A stored value as a fold can leave it: random below p or (15%) a value below C0 plus p. returns: integer."""
    return rg.randrange(p) if rg.random() < 0.85 else rg.randrange(C0) + p


def represent(x, t, a, b, rg, tm=None):
    """A sign's six fractions for leaf t = x through the truncated products on a random representation (3.3):
    w weak, two of E, F, G random below p, the third solved, beta w, rows by SHX, corrections. returns: (fractions,
    (alpha, beta', gamma))."""
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
    """A leaf at an edge of D' ("D": a T(v_m) within 420 of a multiple of 2^28) or W ("Xm": a v_m within 6*2^192 of 0
    or p; "diff": a (beta^i - beta^j) x so). rg: Random; kind: text. returns: field element."""
    m = rg.randrange(3)
    if kind == "D":
        r28 = rg.choice((rg.randrange(420), 2 ** 28 - 1 - rg.randrange(420), rg.randrange(2 ** 28)))
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
    """A leaf near an edge of W: a small value divided by one of the six forms whose small values make W (3.3). rg:
    Random. returns: field element."""
    mlt = rg.choice([1, BETA, BETA * BETA, 1 - BETA, 1 - BETA * BETA, BETA - BETA * BETA])
    dd = rg.randrange(-(1 << rg.randrange(150, 250)), 1 << rg.randrange(150, 250))
    return dd * pow(mlt % p, -1, p) % p or 1


def index_check(x, rg):
    """KEY's index from a random representation of x (represent) against the exact index. x: leaf; rg: Random.
    returns: (whether x is in W, whether the index is wrong)."""
    t, a, b = rg.randrange(3), rg.randrange(1, p), rg.randrange(1, p)
    j = record_fields(leaf_w(represent(x, t, a, b, rg)[0], t))[0]
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
    """Every member of a class, randomly represented, emits the class label unless in W or D'; errors and terms in
    their bounds (3.3). rg: Random. returns: dict."""
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
        for m in range(3):              # the member x with beta^m x = v*: exact index m
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
A_HB, A_YB, A_SB, A_SP = (1 << 41) + 0x100000, (1 << 41) + 0x200000, 1 << 42, (1 << 41) + 0x7000000


class LazyTable(dict):
    """Memory whose unwritten table entries read as their initial values (init3's and tinit's 2^36 entries in part A's
    production-size runs, proof.md 2.3)."""

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


def batch_tables(mem, Hb, Yb, Sb, X, Ys):
    """X's header H with beta's rows at H - 144, the Y table and the entries' correction words (3.1).
    returns: None."""
    a, y = X
    for i, c in enumerate((2 * a ** 3, 2 * a, a * a, BETA)):
        for k, v in enumerate(rows(c % p)):
            put(mem, Hb + 128 * i + 32 * k - 528 * (i == 3), v)
    put(mem, Hb + 384, y)
    put(mem, Hb + 416, pow(a, 4, p))
    for k, w in enumerate(corr(a)):
        put(mem, Hb + 448 + 16 * k, w, 2)
    for l, (b, yy) in enumerate(Ys):
        for k, v in enumerate((yy * yy, yy, b, b ** 3, 3 * b * b, b ** 4)):
            put(mem, Yb + 192 * l + 32 * k, v % p)
        for k, w in enumerate(corr(b)):
            put(mem, Sb + 256 * l + 192 + 16 * k, w, 2)


def components(X, Y):
    """The exact components alpha' = E w, beta'' = F w, gamma = G w (proof.md 3.1), w = u^2 then v^2. X, Y: points
    with D != 0. returns: six field elements."""
    a, b = X[0], Y[0]
    E, F, G = (b ** 4 + 2 * a ** 3 * b) % p, (a ** 4 + 2 * a * b ** 3) % p, 3 * a * a * b * b % p
    out = []
    for w in (pow(Y[1] + X[1], -2, p), pow(Y[1] - X[1], -2, p)):
        out += [E * w % p, F * w % p, G * w % p]
    return out


def bwd_model(X, Ys, acc, D, R, tm=None):
    """bwd's twelve stored words per incidence, bit for bit (proof.md 3.1, 3.3, 3.4), from fwd's acc and D and R. tm:
    as for tfp. returns: ({l: the words in FRAC_OFF order}, the last R)."""
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
            w.append(sqf(z * z)[0])
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
    """Exact leaves x(phi^t(X) + Y) and x(phi^t(X) - Y), in the order 2t + s. X, Y: points. returns: list of six field
    elements."""
    out = []
    for t in range(3):
        PX = (pow(BETA, t, p) * X[0] % p, X[1])
        out += [PROD.add(PX, Y)[0], PROD.add(PX, PROD.neg(Y))[0]]
    return out


def fixup_extra(prog):
    """The fixups of a loaded leaf.s (proof.md 3.4). prog: from a64_load. returns: {index of a fixup's STADD:
    instructions the fixup adds, its chain and STADD after the branch}."""
    return {i: i - max(k for k in range(i) if prog[k][0] in ("b.lo", "b.hs", "cbnz")) for i, (mn, _) in enumerate(prog)
            if mn == "stadd"}


def flags_dead_after_fixups(prog, labels):
    """Whether every path from each fixup's join sets the flags before reading them, or returns (proof.md 3.4), where
    a taken branch and a fixup differ only in the flags. prog, labels: loaded leaf.s. returns: bool."""
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
    """fwd, bwd and keys on one batch, counting fixups; fractions against bwd's model and their bounds, 3.1's
    identity, records against keys' model, indices and labels outside W and D'. Ys: a multiple of 4 of points;
    override: {(l, s): (leaf, t, rg)} replaces sign s of incidence l by represent's fractions. returns: dict."""
    fx = fixup_extra(prog)
    vis = {k: dict.fromkeys(fx, 0) for k in ("fwd", "bwd")}
    q2, Sb, SP = len(Ys), A_SB, A_SP
    BASE, END = 1 << 62, (1 << 62) + (1 << 40)
    mem = LazyTable([(0, 8 << 36, machine.PB1 + 32, machine.BLK), (machine.RC, machine.RC + (8 << 36), BASE, 0)])
    batch_tables(mem, A_HB, A_YB, Sb, X, Ys)
    regs = dict(x0=q2, x1=A_YB, x3=Sb, x25=1, x26=0, x27=0, x28=0, sp=SP, x30=SENTINEL)
    regs.update(as_regs(X[1] * X[1] % p, ("x19", "x20", "x21", "x22")))
    out = {}
    nf = a64.a64_run(prog, labels, "fwd", regs, mem, out_regs=out, visits=vis["fwd"])
    words, acc = [out["x0"]], reg_int(out, ("x25", "x26", "x27", "x28"))
    accs, Ds, R0 = [get(mem, Sb + 256 * l + 32) for l in range(q2)], [get(mem, Sb + 256 * l) for l in range(q2)], \
        pow(acc, p - 2, p)
    regs = dict(x1=A_YB + 192 * q2 + 32, x3=Sb + 256 * (q2 - 1), sp=A_HB, x30=SENTINEL)
    regs.update(as_regs(R0, ("x25", "x26", "x27", "x28")))
    out = {}
    nb = a64.a64_run(prog, labels, "bwd", regs, mem, out_regs=out, visits=vis["bwd"])
    words.append(out["x0"])
    want = 1
    for Y in Ys:
        want = want * (Y[1] * Y[1] - X[1] * X[1]) % p
    final_R = reg_int(out, ("x25", "x26", "x27", "x28"))   # the last R update, R_0 D_0 = 1, left unread (3.4)
    sp_kept = out["sp"] == A_HB
    tm = {}
    model, mR = bwd_model(X, Ys, accs, Ds, R0, tm)
    bad_model, bad_bound, bad_leaf, err, leaves = int(mR != final_R), 0, 0, [[99, -99], [99, -99]], []
    for l, Y in enumerate(Ys):
        e, cs = Sb + 256 * l, components(X, Y)
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
        f6[5] = f6[5] >> 64 << 64 | (f6[5] & M64) >> 3          # gamma's V1 in its s-form (0 mod 8 as computed)
        for o, v in zip(FRAC_OFF[s], f6):
            put(mem, Sb + 256 * l + o, v, 2)
        for t_ in range(3):
            leaves[l][2 * t_ + s] = (c3[0] + BT[t_] * c3[1] + BT[2 * t_ % 3] * c3[2]) % p
    port = machine.Port(mem, machine.MB1, machine.PB1)
    regs = dict(x3=Sb + 256 * (q2 - 1) + 160, x12=END, x17=W2, x18=0, sp=SP, x30=SENTINEL)
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
    nfix = {k: [sum(v for i, v in vis[k].items() if fx[i] == e) for e in (3, 4, 6)] for k in vis}
    ext = {k: 3 * nfix[k][0] + 4 * nfix[k][1] + 6 * nfix[k][2] for k in vis}
    return dict(fwd=nf, bwd=nb, keys=nk, acc_ok=acc % p == want, fraction_mismatches=bad_model,
                fractions_outside_bounds=bad_bound, fraction_errors=err, tfp_terms=tm, identity_mismatches=bad_leaf, fixups_3_4_6=nfix, fixup_extra=ext,
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
    """Macro invocations, nested included, in one kernel of leaf.s (label to first ret). src: leaf.s text; entry:
    label; names: macro names. returns: {name: count}."""
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
    """A cube root mod p, or None (p = 7 mod 9, so a cube's root is its power (p + 2)/9). v: residue. returns: integer
    or None."""
    c_ = pow(v, (p + 2) // 9, p)
    return c_ if c_ * c_ * c_ % p == v % p else None


def fixup_points(rg, X):
    """Points Y whose batch takes chosen fixups (proof.md 3.4): fwd's FS forming D_l (and the next incidence's row
    shift of acc = D), bwd's FS of y_Y - y_X and its FADDW of y_Y + y_X. rg: Random; X: point. returns: three points."""
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


def edge_batch(prog, labels, X, Ys, cert, rg, kinds, need):
    """leaf_batch with each sign's fractions replaced by a representation of an edge leaf of the given kinds, the
    first four in W (need "W") or D' (need "D"). returns: leaf_batch's dict."""
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
    """Execute leaf.s and fe_inv on production points; derive the per-incidence and fixed counts. rng: Random; cert:
    certificate module. returns: dict."""
    res = {}
    prog, labels = a64.a64_load(open(os.path.join(HERE, "kernels", "leaf.s")).read())
    X = PROD.mul(rng.randrange(1, r), GP)
    runs = {}
    for q2 in (4, 8, 12):
        Ys = [PROD.mul(rng.randrange(1, r), GP) for _ in range(q2)]
        runs[q2] = leaf_batch(prog, labels, X, Ys, cert)
    res["runs"] = runs
    net = {q: {k: d[k] - d["fixup_extra"][k] for k in ("fwd", "bwd")} for q, d in runs.items()}
    slope = {}
    for k in ("fwd", "bwd"):
        slope[k] = sorted({Fraction(net[b][k] - net[a][k], b - a) for a, b in ((4, 8), (8, 12), (4, 12))})
        res[k + "_per_incidence"] = [float(v) for v in slope[k]]
        res[k + "_fixed"] = sorted({float(net[q][k] - slope[k][0] * q) for q in runs})
    res["fwd_bwd_formula"] = "fwd (bwd) = FWD_FIX + FWD_INC*q2 (BWD_...) + 3, 4 or 6 per fixup (FIX_EXTRA), q2 = 0 mod 4"
    res["fwd_bwd_formulas_hold"] = all(d["fwd_bwd_formulas"] and d["fixup_words_ok"] for d in runs.values())
    res["counts_ok"] = (slope["fwd"] == [cert.FWD_INC] and slope["bwd"] == [cert.BWD_INC]
                        and res["fwd_fixed"] == [cert.FWD_FIX] and res["bwd_fixed"] == [cert.BWD_FIX])
    res["fixups"] = dict(sites=len(fixup_extra(prog)), flags_dead_at_every_join=flags_dead_after_fixups(prog, labels))
    res["loop_ends_on_bit_33"] = [o for mn, o in prog if mn in ("tbz", "tbnz")] == [["x3", "#33", "1b"]] * 2   # bwd, keys
    res["keys_formula"] = "keys = KEY_FIX + KEY_INC*q2 + x18's count (JPATH per j = 1, 2 leaf) + FLUSH per block sent" \
                          " + KEY_BLOCK_DROP per block dropped"
    res["keys_formula_holds"] = all(d["keys_formula"] for d in runs.values())
    res["tfp_terms"] = dict(largest={k: max(d["tfp_terms"][k] for d in runs.values()) for k in TERMS}, bounds=TERMS,
                            sums_exact=all(d["tfp_terms"]["sum_ok"] for d in runs.values()))
    res["fraction_errors"] = [[min(d["fraction_errors"][g][0] for d in runs.values()), max(d["fraction_errors"][g][1]
                              for d in runs.values())] for g in (0, 1)]     # alpha and beta' (bound [-1, 64]); gamma
    res["all_correct"] = all(d["acc_ok"] and d["fraction_mismatches"] == 0 and d["fractions_outside_bounds"] == 0
                             and d["identity_mismatches"] == 0
                             and d["records_ok"] and d["meter"] == d["meter_ref"] and d["keys_not_canonical"] == 0
                             and d["wrong_index"] == 0 and d["sp_kept"] for d in runs.values())
    erng = random.Random(SEED + 49979687)   # own stream
    Ys = [PROD.mul(rng.randrange(1, r), GP) for _ in range(4)]
    for name, kinds, need in (("window_batch", ("Xm", "diff"), "W"), ("boundary_batch", ("D",), "D")):
        w = edge_batch(prog, labels, X, Ys, cert, erng, kinds, need)
        res[name] = dict(leaves=24, in_W=w["in_W"], in_D_outside_W=w["in_D"], wrong_index=w["wrong_index"],
                         keys_not_canonical=w["keys_not_canonical"], records_ok=w["records_ok"], ok=(w[{"W": "in_W", "D": "in_D"}[need]] >= 4 and w["wrong_index"] == 0
                             and w["keys_not_canonical"] == 0 and w["records"] == 24 and w["records_ok"]
                             and w["meter"] == w["meter_ref"] and w["keys_formula"] and w["fwd_bwd_formulas"]
                             and w["fixup_words_ok"]))
    Yf = fixup_points(erng, X)
    f = leaf_batch(prog, labels, X, [Yf[0], PROD.mul(erng.randrange(1, r), GP)] + Yf[1:], cert)
    res["fixup_batch"] = dict(incidences=4, fixups_3_4_6=f["fixups_3_4_6"], fraction_mismatches=f["fraction_mismatches"],
                              formulas=f["fwd_bwd_formulas"], ok=f["fraction_mismatches"] == 0 and f["fractions_outside_bounds"] == 0
                              and f["fixup_words_ok"] and f["fwd_bwd_formulas"] and f["records_ok"]
                              and f["keys_formula"] and f["fixups_3_4_6"]["fwd"][0] >= 1
                              and f["fixups_3_4_6"]["fwd"][1] >= 1 and f["fixups_3_4_6"]["bwd"][1] >= 2)
    src = open(os.path.join(HERE, "kernels", "leaf.s")).read()
    gp, gl = a64.a64_load(src.replace("\tb.lo 9f\n", "\tcbnz xzr,9f\n").replace("\tb.hs 9f\n", "\tcbnz xzr,9f\n"))
    g = leaf_batch(gp, gl, X, Ys, cert)     # every fixup branch never taken: each fixup runs in place, as an inline chain
    tot = sum(cert.FIX_SITES.values())
    ext = sum(cert.FIX_EXTRA[m] * cert.FIX_SITES[m] for m in cert.FIX_SITES)
    res["forced_fixups"] = dict(incidences=4, fixups_3_4_6=g["fixups_3_4_6"], fixup_extra=g["fixup_extra"],
                                fraction_mismatches=g["fraction_mismatches"], formulas=g["fwd_bwd_formulas"],
                                ok=g["fraction_mismatches"] == 0 and g["fractions_outside_bounds"] == 0
                                and g["fwd_bwd_formulas"] and g["fixup_words_ok"] and g["records_ok"]
                                and sum(map(sum, g["fixups_3_4_6"].values())) == 4 * tot
                                and sum(g["fixup_extra"].values()) == 4 * ext)
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
    names = ("SHX", "SHXT", "FDX", "SQF", "FS", "FADDW")
    counts = {k: macro_counts(src, k, names) for k in ("fwd", "bwd", "keys")}
    per = {n: sum(c[n] for c in counts.values()) // UNROLL for n in names}
    per["SHX"] += per.pop("SHXT")      # SHXT is a row shift without its low word (3.4)
    res["site_counts"] = dict(counts, per_incidence=per, incidences_per_loop_pass=UNROLL)
    res["sites_match_certificate"] = (all(sum(c[n] for c in counts.values()) % UNROLL == 0 for n in names)
                                      and {n: per[n] for n in cert.SITES_PER_INCIDENCE} == cert.SITES_PER_INCIDENCE
                                      and {n: per[n] for n in cert.FIX_SITES} == cert.FIX_SITES
                                      and sum(counts["keys"].values()) == 0)
    srng = random.Random(SEED + 15485863)   # own stream
    res["site_checks"] = sites.site_checks(sys.modules[__name__], srng, fprog, flabels)
    res["final_R_one_and_unused"] = all(d["final_R_is_one"] and d["records_ok"] for d in runs.values())
    prog_r, lab_r = a64.a64_load(open(os.path.join(HERE, "kernels", "route.s")).read())
    res["program_instructions"] = dict(leaf=len(prog), route=len(prog_r), secp_fe=len(fprog),
                                       total=len(prog) + len(prog_r) + len(fprog))
    return res


def lattice_point(tv, tu):
    """The point (v, u), u = beta v mod p, from rounding (tv, tu) in 3.3's basis (B, A - B), (A - B, -A), within about
    2^128 of it in each coordinate. returns: (v, u)."""
    c1 = (2 * (FA * tv + (FA - FB) * tu) + p) // (2 * p)
    c2 = (2 * ((FA - FB) * tv - FB * tu) + p) // (2 * p)
    return c1 * FB + c2 * (FA - FB), c1 * (FA - FB) - c2 * FA


def planted_x(rg, d1, d2, v=None, j=None):
    """A member, outside W and D', of a class with digits d1, d2 and a cell below 16 (or of class v), v* the
    lattice point nearest its blocks' quarter points. returns: (x, v*)."""
    while True:
        vc = v
        if vc is None:
            c1 = rg.getrandbits(64) << 36 | d1
            c2 = (rg.getrandbits(7) << 57 | rg.randrange(16)) << 36 | d2
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
    """The certificate's verdict, claim and the count constants checked against parts A and B. returns: dict."""
    names = ("LEAVES", "FWD_INC", "FWD_FIX", "BWD_INC", "BWD_FIX", "KEY_INC", "KEY_FIX", "EMIT", "EMIT_INC", "JPATH", "FIX_SITES", "FIX_EXTRA", "DELTA_F", "INV", "TINIT_ENTRY", "K", "FLUSH",
             "KEY_BLOCK_DROP", "H2S_BLOCK_DROP", "FIN_DEST", "FIN_SENT", "FIN_HELD", "FIN_FIX", "H2S_REC", "H2S_SLOT",
             "H2S_FIX", "H2G_REC", "H2G_SLOT", "H2G_FIX", "H2G_STEP", "H2G_HIT", "M_CELL", "ZERO_PASS", "ZERO_BYTES",
             "LOOP8", "INIT3_ENTRY", "DLV", "DLVK", "DONE")
    return dict(passed=not cert.FAILED, failed=cert.FAILED, checks=CERT_LINES, claim="%d.%012d" % divmod(int(cert.CLAIM * 10 ** 12), 10 ** 12), q1=cert.Q1,
                mu=float(cert.T["mu"]), success_lower_bound=float(cert.T["lo"]), per_leaf=float(cert.PER_LEAF),
                d_batch=float(cert.T["d_batch"]), g_def=float(cert.T["g_def"]), d_fixup_meter=float(cert.T["d_f"]),
                d_label_error=float(cert.T["d_lab"]), d_key_meter=float(cert.T["d_j"]), walk_budget=cert.T["cap_W"],
                d_walk_budget=float(cert.T["d_W"]), lemma_u_eps=[float(cert.EPS_GEN), float(cert.EPS_EDGE),
                float(cert.EPS_HOME)],
                counts={k: str(getattr(cert, k)) if not isinstance(getattr(cert, k), (int, tuple)) else getattr(cert, k)
                        for k in names})


def lemma_u_toy():
    """Lemma U at toy scale (3.3): on p' near 2^18 with secp256k1's A'/B', each digit set's largest cell count
    over p' 2^-f against the bound at this scale. returns: dict."""
    bt = 492
    while True:
        bt += 1
        at = round(bt * 0x114ca50f7a8e2f3f657c1108d9d44cfd9 / 0xe4437ed6010e88286f547fa90abfe4c3)
        pt = at * at - at * bt + bt * bt
        if all(pt % q for q in range(2, int(pt ** 0.5) + 1)):
            break
    be, k = (at - bt) * pow(bt, -1, pt) % pt, 6
    assert pow(be, 3, pt) == 1 != be
    labs = []
    for x in range(1, pt):
        u = be * x % pt
        v = min(x, u, be * u % pt)
        labs.append(((v << k) // pt, ((be * v % pt) << k) // pt))
    h = (pt >> k) + 1
    diam = int(max(at * at + bt * bt, (2 * bt - at) ** 2 + (2 * at - bt) ** 2) ** 0.5) + 2
    ew = (h + at) * (h + 2 * at - bt) * 4 ** k / pt ** 2 - 1
    rows = []
    for a, b in ((0, 5), (1, 4), (2, 3), (3, 2), (2, 4), (1, 5)):
        cnt = {}
        for c1, c2 in labs:
            key = (c1 % 2 ** a, c2 % 2 ** b)
            cnt[key] = cnt.get(key, 0) + 1
        bound = (1 + ew) * (1 + 1.75 * 2.0 ** (b - k) + 6 * 2.0 ** (a - k)) + 36 * diam * 2.0 ** (a + k) * (h + 2 * at) / pt ** 2
        rows.append(dict(a=a, b=b, max_over_share=round(max(cnt.values()) * 2 ** (a + b) / pt, 4), bound=round(bound, 4)))
    return dict(p=pt, K=k, digit_sets=rows, ok=all(r["max_over_share"] <= r["bound"] for r in rows))


def main():
    """Run the parts within the time budget and print one JSON object. returns: None."""
    global CERT_LINES
    rng = random.Random(SEED)
    cert, CERT_LINES = load_certificate()
    OUT["input"] = cfg
    OUT["budget"] = dict(seconds=BUDGET, reserve=RESERVE, minimums=dict(class_rule_random=MIN_RANDOM,
                         class_rule_near_tie=MIN_NEAR, class_labels=MIN_CLASSES), note="each runs its minimum, then its "
                         "share of the remaining budget, up to 26,000 class-rule inputs and 12,000 classes")
    OUT["A_identities"] = part_a_identities(rng)
    OUT["A_lemma_u_toy"] = lemma_u_toy()
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
            and a["tfp_terms"]["sums_exact"]
            and all(a["tfp_terms"]["largest"][k] <= v for k, v in TERMS.items())
            and a["final_R_one_and_unused"] and a["loop_ends_on_bit_33"])
    OUT["D_certificate"] = part_d(cert)
    OUT["seconds"] = round(time.time() - T0, 2)
    print(json.dumps(OUT, indent=1))


if __name__ == "__main__":
    main()
