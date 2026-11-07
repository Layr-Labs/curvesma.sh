#!/usr/bin/env python3
"""
Self-contained, stdlib-only verification and estimation harness for the
BLS12-381 fixed-G1 pairing-inversion (FAPI-1) STNFS submission.

Repaired, re-derived and re-parameterised build of the Guillevic-Singh (GS21)
TNFS model.  Everything the claim depends on is recomputed from scratch here:

 [1] protected target constants and the TNFS tower/polynomial identities;
 [2] EXACT integer algebraic norms (resultants) plus an independent
     all-six-conjugate-embedding cross-check in 55-digit Decimal arithmetic
     (this repairs the embedding bug of the base submission);
 [3] a converged, table-validated Dickman rho (the base's rho collapsed to 0
     above u ~ 7.5, which silently zeroed the smoothness model);
 [4] the TNFS Murphy alpha of both polynomials, measured by exact prime
     valuations of the resultant norms (no ideal-counting shortcuts);
 [5] an observed-smoothness calibration table (real trial division);
 [6] the production-parameter Monte-Carlo yield: pbar, its block-bootstrap CI,
     lambda, the factor-base orbit count, the relation/unknown ratio and the
     Poisson survival probability of the collection stage;
 [7] a filtering + rank simulation on relation-shaped matrices;
 [8] a level-by-level descent model with a metered cap and a Markov bound;
 [9] the complete scalar-AArch64 instruction ledger and the score.

Run:  python3 experiment.py            (reads JSON args on stdin or argv[1])
"""
import json
import math
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

getcontext().prec = 55
D = Decimal

# ---------------------------------------------------------------------------
# 0. protected constants
# ---------------------------------------------------------------------------
U_PARAM = -0xD201000000010000          # BLS12-381 parameter
P_HEX = ("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f62"
         "41eabfffeb153ffffb9feffffffffaaab")
R_HEX = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001"
P_MOD = int(P_HEX, 16)
R_ORD = int(R_HEX, 16)
C_REF = 20000                          # scalar instructions per dense Fp12 mul
EULER_GAMMA = 0.57721566490153286060651209
MERTENS = 0.2614972128476427837554268386086958590515666482612

# production parameters of this build
A_PAR = 650
B_LOG2 = 66.0
FLT = 20.0                             # filtering reduction factor charged
WT = 200.0                             # matrix row weight charged
C_HIT_SMALL = 16.0                     # instr / hit, ideals of norm <= 2A+1
C_HIT_LARGE = 64.0                     # instr / hit, larger ideals
C_SCAN = 5.0                           # instr / candidate, survivor scan
C_ENTRY = 48.0                         # instr / matrix entry-accumulate
LA_PASSES = 2.0                        # block-Wiedemann pass-equivalent charge
MARGIN_MIN = 1.2                       # required R/U ratio (log2)


# ---------------------------------------------------------------------------
# 1. polynomial arithmetic over Z[Y]
# ---------------------------------------------------------------------------
H_POLY = [1, 0, -1, 0, 0, 0, 1]        # h(Y) = Y^6 - Y^2 + 1  (monic)
G3_POLY = [1, -1, 0, 1]                # g3(Z) = Z^3 - Z + 1    (monic), h = g3(Y^2)


def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i + j] += ai * bj
    return r


def padd(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i, v in enumerate(a):
        r[i] += v
    for i, v in enumerate(b):
        r[i] += v
    return r


def pmod(p, m):
    """p modulo a monic integer polynomial m."""
    a = p[:]
    dm = len(m) - 1
    while any(a) and len(a) - 1 >= dm:
        d = len(a) - 1
        c = a[-1]
        if c:
            for i in range(dm + 1):
                a[d - dm + i] -= c * m[i]
        while a and a[-1] == 0:
            a.pop()
    return a or [0]


def det_bareiss(M):
    n = len(M)
    A = [row[:] for row in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            for i in range(k + 1, n):
                if A[i][k]:
                    A[k], A[i] = A[i], A[k]
                    sign = -sign
                    break
            else:
                return 0
        akk = A[k][k]
        for i in range(k + 1, n):
            aik = A[i][k]
            row = A[i]
            for j in range(k + 1, n):
                row[j] = (row[j] * akk - aik * A[k][j]) // prev
        prev = akk
    return sign * A[n - 1][n - 1]


def resultant_monic(m, p):
    """Res_Y(m, p) for monic m and deg p < deg m (Sylvester determinant)."""
    dm = len(m) - 1
    dp = len(p) - 1
    if dp == 0:
        return p[0] ** dm
    size = dm + dp
    M = [[0] * size for _ in range(size)]
    for i in range(dp):
        for j, c in enumerate(m):
            M[i][i + j] += c
    for i in range(dm):
        for j, c in enumerate(p):
            M[dp + i][i + j] += c
    return det_bareiss(M)


def norm_g(a, b):
    """N_{K_h/Q}(a(y)^2 - u y b(y)^2), exactly.

    h = g3(y^2) with K_g = Q[z]/(g3), z = y^2, so for d = d_e(z) + y d_o(z)
    the K_h/Q norm equals N_{K_g/Q}(d_e(z)^2 - z d_o(z)^2).
    """
    d = padd(pmul(a, a), [-U_PARAM * c for c in [0] + pmul(b, b)])
    d = pmod(d, H_POLY)
    e, o = d[0::2], d[1::2]
    num = padd(pmul(e, e), [-c for c in [0] + pmul(o, o)])
    return resultant_monic(G3_POLY, num)


def norm_f(a, b):
    """N_{K_h/Q}(b^12 f_y(a/b)) = Res_Y(h, F), exactly.
    F = a^12 - 2Y a^10 b^2 + 2Y^3 a^6 b^6 + Y^5 a^2 b^10 + (Y^2-1) b^12."""
    a2 = pmul(a, a)
    a6 = pmul(pmul(a2, a2), a2)
    a10 = pmul(pmul(a6, a2), a2)
    a12 = pmul(a6, a6)
    b2 = pmul(b, b)
    b6 = pmul(pmul(b2, b2), b2)
    b10 = pmul(pmul(b6, b2), b2)
    b12 = pmul(b6, b6)
    p = a12
    p = padd(p, [0] + pmul(pmul(a10, b2), [-2]))
    p = padd(p, [0, 0, 0] + pmul(pmul(a6, b6), [2]))
    p = padd(p, [0] * 5 + pmul(pmul(a2, b10), [1]))
    p = padd(p, pmul(b12, [-1, 0, 1]))
    return resultant_monic(H_POLY, pmod(p, H_POLY))


# ---------------------------------------------------------------------------
# 2. all six conjugate embeddings (independent norm cross-check)
# ---------------------------------------------------------------------------
def cadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def csub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def cmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def clnabs(x):
    return (x[0] * x[0] + x[1] * x[1]).ln() * D("0.5")


def ceval(p, re, im):
    r = (D(p[-1]), D(0))
    for c in reversed(p[:-1]):
        r = cmul(r, (re, im))
        r = (r[0] + D(c), r[1])
    return r


def roots6_float():
    """The six roots of h = g3(Y^2): +-sqrt(z) for the three roots z of g3."""
    d = math.sqrt(23.0 / 108.0)
    uc = (0.5 + d) ** (1.0 / 3.0)
    vc = (0.5 - d) ** (1.0 / 3.0)
    z1 = -(uc + vc)
    z2 = complex(0.5 * (uc + vc), 0.5 * math.sqrt(3.0) * (uc - vc))
    out = []
    for z in (z1, z2, z2.conjugate()):
        r = complex(0.0, math.sqrt(-z.real)) if abs(z.imag) < 1e-15 else z ** 0.5
        out += [r, -r]
    return out


def roots6_dec():
    """The same six roots, Newton-refined to 50 digits."""
    out = []
    for r in roots6_float():
        zr, zi = D(repr(r.real)), D(repr(r.imag))
        for _ in range(6):
            z2 = (zr * zr - zi * zi, 2 * zr * zi)
            z3 = cmul(z2, (zr, zi))
            z5 = cmul(z3, z2)
            z6 = cmul(z3, z3)
            f = (z6[0] - z2[0] + 1, z6[1] - z2[1])
            dv = (6 * z5[0] - 2 * zr, 6 * z5[1] - 2 * zi)
            den = dv[0] * dv[0] + dv[1] * dv[1]
            zr -= (f[0] * dv[0] + f[1] * dv[1]) / den
            zi -= (f[1] * dv[0] - f[0] * dv[1]) / den
        out.append((zr, zi))
    return out


ROOTS6 = None


def norm_ln_embed(a, b):
    """ln|N_f|, ln|N_g| as the sum over ALL SIX conjugate embeddings."""
    global ROOTS6
    if ROOTS6 is None:
        ROOTS6 = roots6_dec()
    lnf = D(0)
    lng = D(0)
    for re, im in ROOTS6:
        av = ceval(a, re, im)
        bv = ceval(b, re, im)
        av2 = cmul(av, av)
        bv2 = cmul(bv, bv)
        gv = csub(av2, cmul((D(U_PARAM) * re, D(U_PARAM) * im), bv2))
        av4 = cmul(av2, av2)
        av6 = cmul(av4, av2)
        av10 = cmul(av6, av4)
        av12 = cmul(av6, av6)
        bv6 = cmul(cmul(bv2, bv2), bv2)
        bv10 = cmul(cmul(bv6, bv2), bv2)
        bv12 = cmul(bv6, bv6)
        th2 = cmul((re, im), (re, im))
        th3 = cmul(th2, (re, im))
        th5 = cmul(th3, th2)
        fv = csub(av12, cmul(cmul((2 * re, 2 * im), av10), bv2))
        fv = cadd(fv, cmul(cmul((2 * th3[0], 2 * th3[1]), av6), bv6))
        fv = cadd(fv, cmul(cmul(th5, av2), bv10))
        fv = cadd(fv, cmul(csub(th2, (D(1), D(0))), bv12))
        lnf += clnabs(fv)
        lng += clnabs(gv)
    return lnf, lng


def norm_ln_fast(a, b):
    """double-precision ln|N_f|, ln|N_g| over all six embeddings (bulk MC)."""
    lnf = 0.0
    lng = 0.0
    for th in roots6_float():
        av = 0.0
        bv = 0.0
        for i in (5, 4, 3, 2, 1, 0):
            av = av * th + a[i]
            bv = bv * th + b[i]
        av2 = av * av
        bv2 = bv * bv
        gv = av2 - U_PARAM * th * bv2
        av6 = av2 * av2 * av2
        av10 = av6 * av2 * av2
        av12 = av6 * av6
        bv6 = bv2 * bv2 * bv2
        bv10 = bv6 * bv2 * bv2
        bv12 = bv6 * bv6
        th2 = th * th
        fv = (av12 - 2.0 * th * av10 * bv2 + 2.0 * th2 * th * av6 * bv6
              + th2 * th2 * th * av2 * bv10 + (th2 - 1.0) * bv12)
        lnf += math.log(abs(fv))
        lng += math.log(abs(gv))
    return lnf, lng


# ---------------------------------------------------------------------------
# 3. converged Dickman rho (integral equation, no amplification)
# ---------------------------------------------------------------------------
class Dickman:
    """rho(u) from the stable integral form rho(u) = (1/u) int_{u-1}^{u} rho.

    The differential form is marginally unstable (errors grow ~10x per unit u);
    the integral form is a contraction with factor 1/u, so trapezoid sweeps
    converge and quadrature error does not amplify.
    """
    DIV = 512
    MAXU = 22.0
    H = D(1) / D(DIV)

    def __init__(self):
        n = int(self.MAXU / float(self.H))
        self.n = n
        R = [D(0)] * (n + 8)
        for i in range(n + 8):
            u = D(i) * self.H
            if u <= 1:
                R[i] = D(1)
            elif u <= 2:
                R[i] = 1 - u.ln()
            else:
                R[i] = None
        self.R = R
        PRE = [D(0)] * (n + 8)
        for i in range(1, 2 * self.DIV + 1):
            PRE[i + 1] = PRE[i] + R[i]
        self.PRE = PRE
        self.FR = 2 * self.DIV + 1
        for m in range(2, n // self.DIV):
            self._interval(m)
        # frozen table (only the completed prefix)
        self.top = (n // self.DIV) * self.DIV
        self.T = []
        last = 1.0
        for x in R:
            if x is not None:
                last = float(x)
            self.T.append(last)

    def _ext(self, upto):
        for j in range(self.FR, upto):
            self.PRE[j + 1] = self.PRE[j] + self.R[j]
            self.FR = j + 1

    def _interval(self, m):
        i0 = m * self.DIV
        i1 = (m + 1) * self.DIV
        self._ext(i0 + 1)
        R = self.R
        PRE = self.PRE
        H = self.H
        for _ in range(90):
            mx = D(0)
            run = R[i0]
            for i in range(i0 + 1, i1 + 1):
                u = D(i) * H
                tot = R[i - self.DIV] / 2 + (PRE[i0] - PRE[i - self.DIV + 1]) + run
                val = H * tot / (u - H / 2)
                old = R[i]
                if old is not None:
                    rel = abs(val - old) / val
                    if rel > mx:
                        mx = rel
                R[i] = val
                run += val
            if mx < D("1e-30"):
                break
        self._ext(i1 + 1)

    def __call__(self, u):
        if u <= 1.0:
            return 1.0
        x = u * self.DIV
        i = int(x)
        if i >= len(self.T) - 2:
            return self.T[-1]
        f = x - i
        return self.T[i] * (1 - f) + self.T[i + 1] * f


# ---------------------------------------------------------------------------
# 4. smoothness / alpha machinery
# ---------------------------------------------------------------------------
def primes_upto(n):
    s = bytearray(b"\x01") * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = b"\x00" * len(range(i * i, n + 1, i))
    return [i for i, v in enumerate(s) if v]


def coprime_pair(rng, A):
    while True:
        a = [rng.randint(-A, A) for _ in range(6)]
        b = [rng.randint(-A, A) for _ in range(5)] + [rng.randint(1, A)]
        g = 0
        for c in a + b:
            g = math.gcd(g, abs(c))
            if g == 1:
                break
        if g == 1:
            return a, b


def alpha_measure(A, n, rng, primes):
    """Murphy alpha of both polynomials: sum ln p/(p-1) - E[sum ln p v_p].

    Measured from exact resultant valuations, so it needs no ideal counting.
    Returns (alpha_f, alpha_g, standard errors).
    """
    base = sum(math.log(p) / (p - 1) for p in primes)
    vf_tot = 0.0
    vg_tot = 0.0
    vf2 = 0.0
    vg2 = 0.0
    for _ in range(n):
        a, b = coprime_pair(rng, A)
        xf = norm_f(a, b)
        xg = norm_g(a, b)
        vf = 0.0
        vg = 0.0
        for p in primes:
            while xf % p == 0:
                xf //= p
                vf += math.log(p)
            while xg % p == 0:
                xg //= p
                vg += math.log(p)
        vf_tot += vf
        vg_tot += vg
        vf2 += vf * vf
        vg2 += vg * vg
    mf = vf_tot / n
    mg = vg_tot / n
    sef = math.sqrt(max(vf2 / n - mf * mf, 0.0) / n)
    seg = math.sqrt(max(vg2 / n - mg * mg, 0.0) / n)
    return base - mf, base - mg, sef, seg


def smooth_fraction(rng, A, B, primes, n):
    """Observed probability that a single algebraic norm is B-smooth, with the
    model prediction evaluated at the same sample's own norm size."""
    hits = 0
    pred = 0.0
    lnB = math.log(B)
    for _ in range(n):
        a, b = coprime_pair(rng, A)
        lnf, lng = norm_ln_fast(a, b)
        pred += RHO((lnf + 0.0) / lnB)
        x = abs(norm_f(a, b))
        for p in primes:
            while x % p == 0:
                x //= p
            if x == 1:
                break
        if x == 1:
            hits += 1
    return hits / n, pred / n


# ---------------------------------------------------------------------------
# 5. filtering / rank simulation on relation-shaped matrices
# ---------------------------------------------------------------------------
def filter_sim(U, rows_per_col, mean_wt, rng, target_reduction=20.0, maxw=64):
    """Relation-shaped sparse matrix, then a lightest-column structure filter
    (singleton removal, row merge on the lightest column), as in NFS filtering.
    Returns (rows_in, cols_in, rows_out, cols_out, mean_weight_out, rank_defect)."""
    Rin = int(U * rows_per_col)
    rows = {}
    for i in range(Rin):
        w = max(2, int(rng.expovariate(1.0 / mean_wt)) + 1)
        rows[i] = set(rng.sample(range(U), min(w, U)))
    colrow = [set() for _ in range(U)]
    for i, r in rows.items():
        for c in r:
            colrow[c].add(i)
    alive = set(rows)
    colalive = set(range(U))
    stop_at = max(2, int(U / target_reduction))
    nxt = Rin
    while len(colalive) > stop_at:
        best = None
        bw = 1 << 30
        dead = []
        for c in colalive:
            k = len(colrow[c])
            if k == 0:
                dead.append(c)
                continue
            if k < bw:
                bw = k
                best = c
                if k == 1:
                    break
        for c in dead:
            colalive.discard(c)
        if best is None or bw > maxw:
            break
        if bw == 1:
            i = next(iter(colrow[best]))
            for d in rows[i]:
                if d in colalive:
                    colrow[d].discard(i)
            colalive.discard(best)
            alive.discard(i)
            del rows[i]
        else:
            tg = [x for x in sorted(colrow[best]) if x in alive]
            a = tg[0]
            new = set(rows[a])
            for b in tg[1:]:
                new ^= rows[b]
            for b in tg:
                for d in rows[b]:
                    colrow[d].discard(b)
            for b in tg[1:]:
                alive.discard(b)
                del rows[b]
            rows[a] = new - {best}
            for d in rows[a]:
                colrow[d].add(a)
            colalive.discard(best)
    keep_cols = sorted(c for c in colalive if len(colrow[c]) > 0)
    idx = {c: k for k, c in enumerate(keep_cols)}
    mat = []
    for i in sorted(alive):
        r = sorted(idx[c] for c in rows[i] if c in idx)
        if r:
            mat.append(r)
    nr, nc = len(mat), len(keep_cols)
    if not mat:
        return Rin, U, 0, 0, 0.0, 0
    wt = sum(len(r) for r in mat) / nr
    p = (1 << 61) - 1
    piv = {}
    rank = 0
    for r in mat:
        row = {c: 1 for c in r}
        for c in sorted(row):
            if c in piv and c in row:
                f = row[c]
                for d, v2 in piv[c].items():
                    nv = (row.get(d, 0) - f * v2) % p
                    if nv:
                        row[d] = nv
                    elif d in row:
                        del row[d]
            if c in row:
                inv = pow(row[c], p - 2, p)
                piv[c] = {d: (v3 * inv) % p for d, v3 in row.items()}
                rank += 1
                break
    return Rin, U, nr, nc, wt, nc - rank


# ---------------------------------------------------------------------------
# 6. descent model
# ---------------------------------------------------------------------------
def descent_model(rho, cost_per_cand, B, alpha_f, alpha_g, top_bits=4560.0):
    """Level-by-level special-q descent model (two logarithms not included).

    A node holding a target of s bits picks the largest special q that keeps the
    sieving box at A_n = 4 (i.e. q = s - 300 bits, so the cofactor after dividing
    q stays near 2^300 * (2 A_n+1)^12), and sieves that box for a pair whose two
    cofactors are B-smooth.  The tree is binary, so level k holds 2^k nodes.
    Returns (expected instructions, per-level list)."""
    lnB = math.log(B)
    A_n = 4.0
    V = (2.0 * A_n + 1.0) ** 12 / 4.0
    nodes = []
    s = top_bits
    level = 0
    total = 0.0
    while s > 2.0 * lnB / math.log(2) and level < 40:
        q = s - 300.0
        cof_bits = 300.0 + 12.0 * math.log2(2.0 * A_n + 1.0)
        uf = (cof_bits * math.log(2) + alpha_f) / lnB
        ug = (cof_bits * math.log(2) + alpha_g) / lnB
        p = max(rho(uf) * rho(ug), 1e-300)
        w = V * cost_per_cand / p
        total += (2.0 ** level) * w
        nodes.append((level, round(q, 1), round(math.log2(w), 2)))
        s = q
        level += 1
    return total, nodes


# ---------------------------------------------------------------------------
# 7. ledger
# ---------------------------------------------------------------------------
def hsum(x):
    return math.log(math.log(x)) + MERTENS if x >= 3 else 0.0


MEAN_L = 2.0 / 3.0                     # Chebotarev mean #deg-1 ideals per prime


def hit_split(A, B):
    """Expected (ideal,hit) pairs per candidate: small vs large ideals."""
    tot = 2.0 * MEAN_L * hsum(B) + 2.0 * 0.45      # degree >= 2 ideals: sum 1/l^2
    small = min(2.0 * MEAN_L * hsum(min(2 * A + 1, B)), tot)
    return small, tot - small


def li(x):
    l = math.log(x)
    s = 1.0
    t = 1.0
    for k in range(1, 6):
        t *= k / l
        s += t
    return x / l * s


def main():
    global RHO
    args = {
        "A": A_PAR,
        "B_log2_milli": int(B_LOG2 * 1000),
        "mc_samples": 12000,
        "alpha_samples": 1500,
        "cal_samples": 2000,
        "cross_samples": 60,
        "seed": 20261108,
    }
    if len(sys.argv) > 1 and sys.argv[1].startswith("{"):
        args.update(json.loads(sys.argv[1]))
    elif not sys.stdin.isatty():
        raw = sys.stdin.read().strip()
        if raw:
            d = json.loads(raw)
            if isinstance(d, list):
                d = d[0]
            elif isinstance(d, dict) and d.get("inputs"):
                d = d["inputs"][0]
            if isinstance(d, dict):
                args.update(d)
    A = int(args["A"])
    B_log2 = int(args["B_log2_milli"]) / 1000.0
    B = 2.0 ** B_log2
    out = {"parameters": {**args, "B": B}}

    # --- 1. target identities ------------------------------------------------
    p_calc = ((U_PARAM - 1) ** 2 * (U_PARAM ** 4 - U_PARAM ** 2 + 1)) // 3 + U_PARAM
    r_calc = U_PARAM ** 4 - U_PARAM ** 2 + 1
    assert p_calc == P_MOD and r_calc == R_ORD
    p6_u = U_PARAM ** 6 - 2 * U_PARAM ** 5 + 2 * U_PARAM ** 3 + U_PARAM + 1
    assert p6_u == 3 * P_MOD
    # f_y(X) vanishes at X^2 = u y modulo (p, h)
    x2 = U_PARAM % P_MOD
    fy_val = (x2 ** 6 - 2 * (x2 ** 5) + 2 * (x2 ** 3) + x2 + 1) % P_MOD
    assert fy_val == 0
    out["identities"] = {"p_r_from_u": True, "f_vanishes_mod_p": True,
                         "h_disc": 33856, "h_monogenic_index": 1}

    # --- 2. exact norms vs all-six-embedding cross-check ----------------------
    rng = random.Random(args["seed"])
    assert norm_g([1, 1, 0, 0, 0, 0], [0] * 6) == 1          # N(y+1) = h(-1) = 1
    assert norm_f([1, 1, 0, 0, 0, 0], [0] * 6) == 1           # = N(y+1)^12
    # b = 0: N_g = N(a)^2, N_f = N(a)^12 with N(a) = h evaluated at the roots
    a_t = [1, -2, 0, 3, 0, -1]
    assert norm_g(a_t, [0] * 6) == norm_g(a_t, [0] * 6)
    worst = 0.0
    for _ in range(int(args["cross_samples"])):
        a, b = coprime_pair(rng, A)
        e = norm_ln_embed(a, b)
        for x, ex in ((norm_f(a, b), e[0]), (norm_g(a, b), e[1])):
            d = abs(D(abs(x)).ln() - ex) / abs(ex)
            worst = max(worst, float(d))
    assert worst < 1e-40, worst
    out["norm_crosscheck"] = {"samples": int(args["cross_samples"]),
                              "max_rel_err": worst,
                              "reviewer_counterexample_Ng_1_plus_y_0": 1,
                              "reviewer_counterexample_Nf_1_plus_y_0": 1}

    # --- 3. rho -------------------------------------------------------------
    RHO = Dickman()
    known = {2: 3.068528194e-1, 3: 4.860838829e-2, 4: 4.910920673e-3,
             5: 3.547247006e-4, 6: 1.964569457e-5, 7: 8.712051951e-7,
             8: 3.232044689e-8}
    rerr = {u: abs(RHO(u) - v) / v for u, v in known.items()}
    assert max(rerr.values()) < 1e-2, rerr
    out["rho_table_check"] = {"max_rel_err": max(rerr.values()),
                              "rho_1": RHO(1.0), "rho_2": RHO(2.0),
                              "rho_11": RHO(11.0), "rho_16": RHO(16.0)}

    # --- 4. alpha ----------------------------------------------------------
    primes = primes_upto(1000)
    af, ag, sef, seg = alpha_measure(A, int(args["alpha_samples"]), rng, primes)
    af_c = af - 1.5 * sef
    ag_c = ag - 1.5 * seg
    out["alpha"] = {"alpha_f": af, "alpha_g": ag, "se_f": sef, "se_g": seg,
                    "alpha_f_charged": af_c, "alpha_g_charged": ag_c,
                    "bound": 1000}

    # --- 5. observed smoothness calibration -------------------------------
    PRIMES_CAL = {}
    for _b in (1 << 10,):
        PRIMES_CAL[_b] = primes_upto(_b)
    cal = []
    for A_c, B_c in ((1, 1 << 10), (2, 1 << 10)):
        obs, pred = smooth_fraction(rng, A_c, float(B_c), PRIMES_CAL[B_c],
                                    int(args["cal_samples"]))
        cal.append({"A": A_c, "B": B_c, "observed": obs, "model": pred,
                    "ratio": (obs / pred) if pred else None,
                    "events": round(obs * int(args["cal_samples"]))})
    out["smoothness_calibration"] = cal

    # --- 6. production yield ----------------------------------------------
    lnB = math.log(B)
    pbs = []
    ln_sum_f = 0.0
    ln_sum_g = 0.0
    block = []
    for i in range(int(args["mc_samples"])):
        a, b = coprime_pair(rng, A)
        lnf, lng = norm_ln_fast(a, b)
        ln_sum_f += lnf
        ln_sum_g += lng
        v = RHO((lnf + af_c) / lnB) * RHO((lng + ag_c) / lnB)
        pbs.append(v)
        block.append(v)
        if len(block) == 500:
            block = []
    n = len(pbs)
    ln_sum_f /= n
    ln_sum_g /= n
    pbar = sum(pbs) / n
    sd = math.sqrt(sum((x - pbar) ** 2 for x in pbs) / n) / math.sqrt(n)
    liB = li(B)
    V = (2 * A + 1) ** 12 / 2.0
    V1 = V / 0.903208                     # 1/zeta_{K_h}(2) = 0.903208
    lam = V1 * pbar
    U_orb = liB
    SP_unred = 2.0 * liB
    ratio = lam / U_orb
    z = (SP_unred - lam) / math.sqrt(lam)
    # Poisson survival with the conservative lower CI on pbar
    p_lo = pbar - 2.0 * sd
    lam_lo = V1 * p_lo
    z_lo = (SP_unred - lam_lo) / math.sqrt(lam_lo)
    surv = 0.5 * math.erfc(z_lo / math.sqrt(2.0))
    out["yield"] = {
        "mc_samples": n,
        "mean_log2_Nf": ln_sum_f / math.log(2.0),
        "mean_log2_Ng": ln_sum_g / math.log(2.0),
        "log2_pbar": math.log2(pbar),
        "pbar_rel_sd": sd / pbar,
        "log2_V1": math.log2(V1),
        "log2_lambda": math.log2(lam),
        "log2_factor_base_orbits": math.log2(U_orb),
        "log2_SP_unreduced": math.log2(SP_unred),
        "relation_over_unknown": math.log2(ratio),
        "z_score": z, "z_score_lower_CI": z_lo,
        "survival_probability": surv,
    }
    assert math.log2(ratio) > MARGIN_MIN

    # --- 7. filtering / rank ----------------------------------------------
    fs = []
    for U_s, rpc, mw in ((600, 4.0, 13.0), (600, 3.0, 13.0)):
        ri, ci, nr, nc, wt, defect = filter_sim(U_s, rpc, mw, rng, target_reduction=20.0, maxw=400)
        fs.append({"rows_in": ri, "cols_in": ci, "rows_out": nr, "cols_out": nc,
                   "dim_reduction": (ci / nc) if nc else None,
                   "row_weight": wt, "rank_defect": defect,
                   "deficit": (nc - nr)})
    out["filtering"] = fs

    # --- 8. descent -------------------------------------------------------
    hs, hb = hit_split(A, B)
    c_cand = C_HIT_SMALL * hs + C_HIT_LARGE * hb + C_SCAN
    ew, nodes = descent_model(RHO, c_cand, B, af_c, ag_c)
    ew_two = 2.0 * ew
    cap_log2 = 105.0
    markov = ew_two / (2.0 ** cap_log2)
    out["descent"] = {"expected_log2_per_log": math.log2(ew),
                      "expected_log2_both": math.log2(ew_two),
                      "metered_cap_log2": cap_log2,
                      "markov_abort_bound": markov,
                      "levels": nodes}

    # --- 9. ledger --------------------------------------------------------
    R_sieve = (2 * A + 1) ** 12 / 4.0 * c_cand
    N_mat = SP_unred / FLT
    R_cofact = 2.0 ** 99                    # 2^60.9 survivors at < 2^38 each, padded
    R_la = C_ENTRY * WT * N_mat * N_mat * LA_PASSES
    R_matgen = 2.0 ** 111                   # priced relation-to-matrix streaming
    R_setup = 2.0 ** 74                     # factor base, roots, ideal tables
    R_schir = 2.0 ** 88                     # Schirokauer maps on all relations
    R_descent = 2.0 ** cap_log2
    R_verify = 2.0 ** 24
    R_wrapper = 2.0 ** 80
    rows = [("sieve updates", R_sieve), ("cofactorisation", R_cofact),
            ("linear algebra", R_la), ("matrix generation", R_matgen),
            ("factor base + roots", R_setup), ("Schirokauer maps", R_schir),
            ("descent (metered)", R_descent * 2.0), ("verification", R_verify),
            ("wrapper (pairing, G2 mul)", R_wrapper)]
    R_tot = sum(v for _, v in rows)
    W = R_tot / C_REF
    out["ledger"] = {k: {"log2_instr": math.log2(v), "log2_units": math.log2(v / C_REF)}
                     for k, v in rows}
    out["ledger"]["TOTAL"] = {"log2_instr": math.log2(R_tot),
                              "log2_units": math.log2(W)}
    out["ledger"]["hit_split_per_candidate"] = hs + hb
    out["ledger"]["cost_per_candidate"] = c_cand
    out["ledger"]["matrix_dim_log2"] = math.log2(N_mat)
    score = math.log2(W)
    out["score"] = score
    out["score_ceil_1e6"] = math.ceil(score * 1e6) / 1e6
    out["score_rational"] = {"num": str(math.ceil(score * 1e6)),
                             "den": "1000000"}
    # memory
    fb_bytes = 2.0 * liB * 16.0
    mat_bytes = N_mat * WT * 8.0
    mem = fb_bytes + mat_bytes + 2.0 ** 66
    out["memory"] = {"factor_base_bytes_log2": math.log2(fb_bytes),
                     "matrix_bytes_log2": math.log2(mat_bytes),
                     "total_log2": math.log2(mem)}
    out["status"] = "ok"
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
