"""Evidence for stnfs-aslice-v6 (bls12-381-fixed-g2-invert-g1).  Stdlib only; one run.

The attack solves the discrete logarithm of z = e(P0,Q)^a in the order-r subgroup of F_{p^12}^*
(the route rules.md admits) with a Special Tower Number Field Sieve, then returns P = [a]P0.

This revision is a re-derivation from the ground up after v1-v4 were rejected.  Its central
correction is the UNIT-ORBIT accounting of relation collection: the unit group O_h^* of the tower
field K_h acts freely on candidate pairs (a,b) -> (eps a, eps b) and all pairs of one orbit give
the SAME relation (same ideal factorisation; log eps = 0 mod r since eps lies in F_{p^6}).  v1-v4
credited every point of the sieve region although a region point has on average ~8 unit
multiples inside the region.  v5 sieves a fundamental domain F of the unit action (canonical a)
with an architecture built around it (per-a 6-dimensional slices), measures all kernel constants
as retired instructions, simulates filtering with the production ideal-occurrence profile, and
prices the descent by an explicit dynamic programme.

Parts (every number printed is recomputed here from the embedded data):
  [0] exact check Res_y(h, Res_x(f, g)) = 0 mod p (the two sides share a root in F_{p^12});
  [1] Dickman rho from Marsaglia's exact series, self-tested;
  [U] the unit group of K_h: all products of the two fundamental units (PARI) up to height 4,
      each verified exactly (norm +-1), count compared with the embedded brute-force count;
  [2] held-out pool: for the first `samples` seeds the point is regenerated, its exact norms and
      its canonical-a flag are recomputed and compared with the embedded inside set;
  [3] region geometry: volume fraction, canonical fraction f_c, slice-size distribution;
  [4] the ARCH-E cost model at the claimed operating point (and its neighbourhood): relation
      collection with the measured slice-sieve constants, filtering credit from [5], LA;
  [5] filtering simulation at reduced scale (production occurrence profile) + embedded
      larger-scale convergence table and the DGP21 record calibration;
  [6] individual-logarithm descent: dynamic programme over 1-bit size levels (boot + tree);
  [7] LA inner-loop instruction count re-derived from the shipped assembly;
  [8] boot lattice: exact membership and |det| = p^6 for the embedded example;
  ledger: total log2 W.
"""
import base64
import json
import math
import os
import random
import sys
import zlib
from decimal import Decimal, getcontext

_HERE = os.path.dirname(os.path.abspath(__file__))
# the verifier runs this file as /snapshot/entrypoint.py with the submission files under /snapshot/files/
_DIR = os.path.join(_HERE, "files") if os.path.exists(os.path.join(_HERE, "files", "filtsim.py")) else _HERE
sys.path.insert(0, _DIR)
import filtsim          # noqa: E402  (filtering simulator; also provides rho)
import descent_fast     # noqa: E402  (descent dynamic programme)
import descent_core     # noqa: E402

LN2 = math.log(2)
UNIT = 20000.0                     # 1 unit = one dense Fp12 multiplication = 20,000 AArch64 instructions

# ------------------------------------------------------------------ target and polynomials ----
P = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
U = 15132376222941642752            # |BLS12-381 seed u|; g = x^2 + |u| y
H = [-1, -1, 1, 1, 0, 0, 1]         # h(y) = y^6 + y^3 + y^2 - y - 1 (v1-v4 tower; h(0) = -1 is what makes the tower generate F_{p^12}, part [0b])
# f = y^6 * 3p(x^2/y) reduced mod h, p(X) the BLS12 characteristic polynomial: (x-degree, y-polynomial)
FTERMS = [(12, [1]), (10, [0, -2]), (6, [0, 0, 0, 2]), (2, [0, 0, 0, 0, 0, 1]), (0, [1, 1, -1, -1])]
ALPHA_F = 4.868                     # log2; GS alpha_tnfs_2d (2000 samples) on (f, h); v1 used 4.721
ALPHA_G = 0.844                     # log2; same code on (g, h)
INV_ZETA = 0.9168940567189809       # density of O_h-coprime pairs, 1/zeta_Kh(2) (v1-v4)
FUND_UNITS = ([0, -1], [1, 0, -1, -1, 1, -1], [-2, 0, 1, 1, -1, 1])   # PARI bnfinit: fundamental units of O_h (rank 3, torsion +-1)
UNIT_EXP = (30, 12, 6)              # exponent ranges that reach every unit of height <= 4 (max needed: 20, 9, 3)
UNITS_BRUTE = 595                   # brute-force count of units of height <= 4 (mod sign, excluding 1), canon.c
REGULATOR = 0.8578
NBAR = 1.947                        # mean number of degree-1 prime ideals per prime, both sides (v1-v4, exact to 2^16)

# ------------------------------------------------------------------ region (ellipsoid) --------
A0 = 1 << 10                        # pool coordinates are x/A0; the production region is the same shape scaled to A
Q_FRAC = 1.0 / 32                   # volume fraction of the ellipsoid in the box [-1,1]^12
# @@ELL_M@@
ELL_M = [
    [0.5528759885, -0.1260799858, 0.1447295864, -0.002614363166, -0.02018451304, 0.0956632886, 0, 0, 0, 0, 0, 0],
    [-0.1260799858, 0.4730025995, -0.06826093202, 0.08498958762, -0.001395387797, -0.06294781792, 0, 0, 0, 0, 0, 0],
    [0.1447295864, -0.06826093202, 0.4552066288, -0.05662522189, 0.04120452268, -0.02031101183, 0, 0, 0, 0, 0, 0],
    [-0.002614363166, 0.08498958762, -0.05662522189, 0.4367157553, 0.0001512689972, -0.02187465135, 0, 0, 0, 0, 0, 0],
    [-0.02018451304, -0.001395387797, 0.04120452268, 0.0001512689972, 0.4322681594, 0.02026484702, 0, 0, 0, 0, 0, 0],
    [0.0956632886, -0.06294781792, -0.02031101183, -0.02187465135, 0.02026484702, 0.4746279595, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0.5907186542, -0.1192505666, 0.1631112482, -0.003050921392, -0.04000276543, 0.1352719521],
    [0, 0, 0, 0, 0, 0, -0.1192505666, 0.4974396776, -0.08398310872, 0.09084451399, -0.01109118934, -0.08908672312],
    [0, 0, 0, 0, 0, 0, 0.1631112482, -0.08398310872, 0.4933090946, -0.0460977154, 0.03924332861, -0.03828634712],
    [0, 0, 0, 0, 0, 0, -0.003050921392, 0.09084451399, -0.0460977154, 0.4744115241, 0.02143547529, -0.04387611416],
    [0, 0, 0, 0, 0, 0, -0.04000276543, -0.01109118934, 0.03924332861, 0.02143547529, 0.4924220451, 0.05071856806],
    [0, 0, 0, 0, 0, 0, 0.1352719521, -0.08908672312, -0.03828634712, -0.04387611416, 0.05071856806, 0.5676976811],
]
# @@END_ELL_M@@

# ------------------------------------------------------------------ operating point ----------
# @@OPPOINT@@
OP = dict(la=10.50, lb=68.50, excess=1.20, theta=0.40, lx2=64.0, k_f=4, k_g=3)
# @@END_OPPOINT@@
CLAIM = 120.497520   # @@CLAIM@@
MARGIN = 0.35   # required log2 yield margin over excess*U (>= 2.5 bootstrap standard deviations of the embedded pool estimate)

# ------------------------------------------------------------------ measured kernel constants -
# slice7.c (per-a slice sieve, production-width: lines, segments, recentred cosets), retired instructions on AArch64
# (/usr/bin/time -l; KERNELS.md).  LINE path (p < 2^P_MID): per bucket update, by expected points per line coset
# h = (line points)/p; RECORD path (p >= 2^P_MID): per hit of the 6-dim walk incl. the 12-byte record and the apply pass,
# by expected points per slice coset h = (slice points)/p.
UPD_TAB = [(15.5, 15.6), (11.5, 23.7), (7.5, 56.6), (3.5, 191.0), (0.6, 1450.0)]   # line path, h = points per line coset (slice7 measure7.sh)
REC_TAB = [(21.0, 37.0), (13.0, 58.0), (6.0, 166.7), (3.0, 576.0), (0.0, 2100.0)]   # record path incl. record + apply, h = points per slice coset
C_COSET_SPARSE = 2130.0  # per (slice, pair) when h < 1 on the record path: the whole exact enumeration of an (almost always empty) coset incl. recentring
P_MID = 40.0             # primes below 2^P_MID are sieved line by line with direct byte updates; above, by records
C_LSET = 866.0           # per (pair, line, coset) setup of the line path: c_k, exact recentring, walk centre (measured)
C_SET2 = 866.0 + 75.0 + 72.0   # per (slice, pair) of the record path: recentring + walk centre (866, the measured line-path setup, 6-dim here) + a(s) mod p and c by Horner/Barrett (<= 75) + 36 streamed words
C_SCAN = 8.0             # per bucket byte at threshold time (load, compare, branch; reports are rare)
C_ENTRY = 35.0           # LA: la_mac.s 33/entry (part [7]) + row finish amortised
BW = 3.0                 # block Wiedemann: ~3 N single-vector products
W_ROW = 100.0            # charged row weight of the filtered matrix (the merge target; N^2 W is minimal there in the simulation)
PAD_N = 2.0              # pad on the simulated matrix dimension
C_NORM = 2.0 ** 21       # exact resultant norms of one survivor
C_TDIV = 2.0 ** 22       # per survivor: trial division by the primes < 2^20 (their hits are not recorded)
LX0 = 8.0                # primes below 2^8 are not sieved (their logs are inside the threshold slack)
DESCENT_CAP = 104.0      # metered cap on both descents together (units)
LOG2_SM = 88.0           # Schirokauer maps: <= 2^37 instr per relation x <= 2^64.5 relations
LOG2_FBC = 70.0          # factor-base construction: 2^62 ideals x (roots 2^20 + 6-dim reduction 2^13 + tables)
LOG2_CANON = 60.0        # canonical a-list: |E_a| x 235 units x ~60 instructions
LOG2_MISC = 60.0         # relation storage/dedupe/filtering passes, ellipsoid fit, challenge pairing
LOG2_VERIFY = 11.0

# ------------------------------------------------------------------ filtering tables ----------
# filtsim at scale 2^17 (kmax 64, W<=150): N/U by (u_f, u_g) and excess.  @@RATIO_GRID@@
RATIO_GRID = {
    (9.5, 6.5): {1.05: 0.0021, 1.1: 0.0020, 1.2: 0.0020, 1.3: 0.0020},
    (9.5, 7.3): {1.05: 0.0029, 1.1: 0.0027, 1.2: 0.0023, 1.3: 0.0022},
    (9.5, 8.0): {1.05: 0.0064, 1.1: 0.0049, 1.2: 0.0037, 1.3: 0.0031},
    (10.3, 6.5): {1.05: 0.0031, 1.1: 0.0029, 1.2: 0.0025, 1.3: 0.0023},
    (10.3, 7.3): {1.05: 0.0087, 1.1: 0.0069, 1.2: 0.0050, 1.3: 0.0039},
    (10.3, 8.0): {1.05: 0.0248, 1.1: 0.0197, 1.2: 0.0127, 1.3: 0.0087},
    (11.0, 6.5): {1.05: 0.0068, 1.1: 0.0057, 1.2: 0.0042, 1.3: 0.0034},
    (11.0, 7.3): {1.05: 0.0199, 1.1: 0.0166, 1.2: 0.0116, 1.3: 0.0086},
    (11.0, 8.0): {1.05: 0.0401, 1.1: 0.0340, 1.2: 0.0244, 1.3: 0.0178},
    (11.7, 6.5): {1.05: 0.0162, 1.1: 0.0132, 1.2: 0.0085, 1.3: 0.0064},
    (11.7, 7.3): {1.05: 0.0368, 1.1: 0.0310, 1.2: 0.0229, 1.3: 0.0174},
    (11.7, 8.0): {1.05: 0.0587, 1.1: 0.0504, 1.2: 0.0388, 1.3: 0.0303},
}
# @@CONV_TABLE@@
CONV_TABLE = [[1.2, 15, 1, 0.0534, 100, 25], [1.2, 15, 2, 0.0532, 100, 24], [1.2, 16, 1, 0.0491, 100, 27], [1.2, 16, 2, 0.0508, 100, 20], [1.2, 17, 1, 0.0518, 100, 22], [1.2, 17, 2, 0.0499, 100, 25], [1.2, 18, 1, 0.0479, 100, 24], [1.3, 15, 1, 0.0461, 100, 22], [1.3, 15, 2, 0.0471, 100, 20], [1.3, 16, 1, 0.0459, 100, 26], [1.3, 16, 2, 0.0452, 100, 21], [1.3, 17, 1, 0.0466, 100, 28], [1.3, 17, 2, 0.0446, 100, 25], [1.3, 18, 1, 0.0419, 100, 26], [1.4, 15, 1, 0.0468, 100, 30], [1.4, 15, 2, 0.0457, 100, 32], [1.4, 16, 1, 0.0441, 100, 21], [1.4, 16, 2, 0.0444, 100, 23], [1.4, 17, 1, 0.0417, 100, 26], [1.4, 17, 2, 0.0412, 100, 30], [1.4, 18, 1, 0.0392, 100, 30]]
# @@END_CONV_TABLE@@

_REPORT = []


def say(*args):
    """progress lines go to stderr; stdout carries only the final JSON object (the verifier's protocol)."""
    line = " ".join(str(a) for a in args)
    _REPORT.append(line)
    print(line, file=sys.stderr, flush=True)


# ------------------------------------------------------------------ exact polynomial arithmetic
def ptrim(p):
    while p and p[-1] == 0:
        p.pop()
    return p


def pmul(p, q):
    if not p or not q:
        return []
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                r[i + j] += a * b
    return r


def padd(p, q):
    n = max(len(p), len(q))
    return ptrim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def pmod_h(p):
    d = len(H) - 1
    p = ptrim(p[:])
    while len(p) - 1 >= d and p:
        shift = len(p) - 1 - d
        f = p[-1]
        if f:
            for i in range(d + 1):
                p[shift + i] -= f * H[i]
        ptrim(p)
    return p


def rmul(p, q):
    return pmod_h(pmul(p, q))


MAT = [[0] * 6 for _ in range(6)]        # multiplication-by-y matrix of O_h = Z[y]/h
for _i in range(5):
    MAT[_i][_i + 1] = 1
for _j in range(6):
    MAT[5][_j] = -H[_j]


def matmul(X, Y):
    return [[sum(X[i][k] * Y[k][j] for k in range(6)) for j in range(6)] for i in range(6)]


def matpoly_eval(Q):
    X = [[0] * 6 for _ in range(6)]
    for d in range(5, -1, -1):
        if any(any(r) for r in X):
            X = matmul(X, MAT)
        qd = Q[d] if d < len(Q) else 0
        if qd:
            for i in range(6):
                X[i][i] += qd
    return X


def det6(X):
    n = 6
    M = [row[:] for row in X]
    sign, prev = 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    M[k], M[i] = M[i], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def res_y_h(Q):
    """Res_y(h, Q) = N_{K_h/Q}(Q(y)) = det of multiplication by Q(y) on Z[y]/h."""
    return det6(matpoly_eval(Q))


def norm_f(a, b):
    apow, bpow = [[1]], [[1]]
    for _ in range(12):
        apow.append(rmul(apow[-1], a))
        bpow.append(rmul(bpow[-1], b))
    acc = []
    for (xd, fy) in FTERMS:
        acc = padd(acc, rmul(pmod_h(fy[:]), rmul(apow[xd], bpow[12 - xd])))
    return abs(res_y_h(pmod_h(acc)))


def norm_g(a, b):
    a2 = rmul(a, a)
    b2 = rmul(b, b)
    return abs(res_y_h(padd(a2, rmul([0, U], b2))))


def gcd_ints(xs):
    g = 0
    for x in xs:
        g = math.gcd(g, x)
    return g


def log2_exact(n):
    k = n.bit_length()
    return k - 60 + math.log2(n >> (k - 60)) if k > 60 else math.log2(n)


# ------------------------------------------------------------------ [0] resultant check ---------
def check_resultant():
    """Res_y(h, Res_x(f, g)) == 0 mod p.  Res_x(f, g) with g = x^2 + U y monic in x: substitute x^2 -> -U y
    into f (f is even in x), giving a polynomial F0(y); Res_x(f, g) = F0(y)^2 up to sign (two roots of g), so
    Res_y(h, Res_x(f,g)) = Res_y(h, F0)^2 and it suffices that Res_y(h, F0) = 0 mod p."""
    F0 = []
    for (xd, fy) in FTERMS:
        term = pmod_h(fy[:])
        for _ in range(xd // 2):
            term = rmul(term, [0, -U])
        F0 = padd(F0, term)
    r = res_y_h(pmod_h(F0))
    ok = r != 0 and r % P == 0
    say(f"[0] Res_y(h, f(x^2 -> -U y)) has {abs(r).bit_length()} bits and is {'0' if ok else 'NOT 0'} mod p "
        f"-> f and g share a root in F_{{p^12}}: {'ok' if ok else 'FAIL'}")
    return ok


def f6_mulp(A, B):
    return [c % P for c in pmod_h(pmul([c % P for c in A], [c % P for c in B]))]


def f6_pow(A, n):
    r = [1]
    base = [c % P for c in A]
    while n:
        if n & 1:
            r = f6_mulp(r, base)
        base = f6_mulp(base, base)
        n >>= 1
    return (r + [0] * 6)[:6]


def check_tower():
    """(i) h is irreducible mod p: y^(p^6) = y and y^(p^k) != y for k = 1, 2, 3 in F_p[y]/h (so the minimal polynomial of y
    has degree 6); (ii) u*y is a non-square in F_{p^6}: (u y)^((p^6-1)/2) = -1 (explicit) -- equivalently N(u y) = u^6 h(0)
    = -u^6 is a non-square in F_p because p = 3 mod 4; so g = x^2 + |u| y = x^2 - u y is irreducible over F_{p^6} and
    F_p[y,x]/(h,g) = F_{p^12}."""
    yp = f6_pow([0, 1], P)                      # y^p = Frobenius(y)
    frob = yp
    ok_irr = True
    cur = [0, 1, 0, 0, 0, 0]
    for k in range(1, 7):
        # cur = y^(p^k): apply Frobenius by composition cur(y) -> cur(y^p)
        acc = []
        for coef in reversed(cur):
            acc = f6_mulp(acc, frob) if acc else []
            acc = [c % P for c in padd(acc, [coef % P])]
        cur = (acc + [0] * 6)[:6]
        is_y = cur == [0, 1, 0, 0, 0, 0]
        if k in (1, 2, 3) and is_y:
            ok_irr = False
        if k == 6 and not is_y:
            ok_irr = False
    uy = [(-U) % P, 0, 0, 0, 0, 0]
    uy = f6_mulp(uy, [0, 1])                    # u*y with u = -U
    chi = f6_pow(uy, (P ** 6 - 1) // 2)
    ok_sq = chi == [(P - 1), 0, 0, 0, 0, 0]
    say(f"[0b] tower generates F_{{p^12}}: h irreducible mod p (y^(p^6)=y, y^(p^k)!=y for k<=3): {'ok' if ok_irr else 'FAIL'}; "
        f"p mod 4 = {P % 4}, N(y) = h(0) = {H[0]}; (u y)^((p^6-1)/2) = -1 in F_p[y]/h: {'ok' if ok_sq else 'FAIL'} -> g irreducible over F_{{p^6}}")
    return ok_irr and ok_sq


# ------------------------------------------------------------------ [1] Dickman rho -----------
rho = filtsim.rho


def check_rho():
    ref = {2.0: 0.30685281944, 3.0: 0.04860838829, 5.0: 3.5472470045e-4, 10.0: 2.77017183e-11, 20.0: 2.46e-29}
    worst = 0.0
    for u, v in ref.items():
        worst = max(worst, abs(rho(u) / v - 1))
    say(f"[1] Dickman rho table: max relative deviation from reference values {worst:.2e} {'ok' if worst < 3e-3 else 'FAIL'}")
    return worst < 3e-3


def proba(log2n, alpha, log2b):
    """Murphy/Dickman B-smoothness probability of an integer of log2n bits with alpha correction."""
    u = (log2n + alpha) / log2b
    return rho(u) + (1 - 0.5772156649015329) * rho(u - 1.0) / (log2n * LN2)


# ------------------------------------------------------------------ [U] units -------------------
def unit_inverse(u):
    """inverse of a unit in O_h by exact linear algebra on the multiplication matrix."""
    from fractions import Fraction
    cols = [(rmul(u, [1 if j == i else 0 for j in range(6)]) + [0] * 6)[:6] for i in range(6)]
    A = [[Fraction(cols[j][i]) for j in range(6)] + [Fraction(1 if i == 0 else 0)] for i in range(6)]
    for k in range(6):
        piv = next(i for i in range(k, 6) if A[i][k] != 0)
        A[k], A[piv] = A[piv], A[k]
        d = A[k][k]
        A[k] = [v / d for v in A[k]]
        for i in range(6):
            if i != k and A[i][k]:
                f = A[i][k]
                A[i] = [vi - f * vk for vi, vk in zip(A[i], A[k])]
    inv = [A[i][6] for i in range(6)]
    assert all(c.denominator == 1 for c in inv)
    return ptrim([int(c) for c in inv])


def build_units(height=4):
    gens = [list(x) for x in FUND_UNITS]
    invs = [unit_inverse(g) for g in gens]
    for g, gi in zip(gens, invs):
        assert ptrim(rmul(g, gi)) == [1]

    def powers(u, ui, n):
        out = {0: [1]}
        for k in range(1, n + 1):
            out[k] = rmul(out[k - 1], u)
            out[-k] = rmul(out[-(k - 1)], ui)
        return out
    P = [powers(g, gi, n) for g, gi, n in zip(gens, invs, UNIT_EXP)]
    units = {}
    for b in range(-UNIT_EXP[1], UNIT_EXP[1] + 1):
        for c in range(-UNIT_EXP[2], UNIT_EXP[2] + 1):
            pbc = rmul(P[1][b], P[2][c])
            for a in range(-UNIT_EXP[0], UNIT_EXP[0] + 1):
                u = (rmul(P[0][a], pbc) + [0] * 6)[:6]
                if max(abs(x) for x in u) <= height:
                    key = tuple(u)
                    nkey = tuple(-x for x in u)
                    if key in units or nkey in units:
                        continue
                    units[key] = (a, b, c)
    bad = sum(1 for k in units if abs(res_y_h(list(k))) != 1)
    one = tuple([1, 0, 0, 0, 0, 0])
    lst = [list(k) for k in units if k != one]
    say(f"[U] units of O_h = products of the three fundamental units with |exponents| <= {UNIT_EXP}: {len(lst)} of height <= {height} "
        f"(mod sign, excluding 1), all with exact norm +-1: {'ok' if bad == 0 else 'FAIL'}; brute-force count (canon.c): {UNITS_BRUTE} "
        f"{'ok' if len(lst) == UNITS_BRUTE else 'MISMATCH'}; rank 3, torsion +-1, regulator {REGULATOR}")
    return lst


# ------------------------------------------------------------------ region / pool ---------------
def quad_blocks(xa, xb):
    qa = 0.0
    qb = 0.0
    for i in range(6):
        ra = ELL_M[i]
        rb = ELL_M[6 + i]
        ta = 0.0
        tb = 0.0
        for j in range(6):
            ta += ra[j] * xa[j]
            tb += rb[6 + j] * xb[j]
        qa += ta * xa[i]
        qb += tb * xb[i]
    return qa, qb


POOL_TAG = "stnfs-v5-ours-tpool:%d"
THETA_MAX = 0.6     # @@THETA_MAX@@  the pool is uniform over D = {Q_a(a) <= THETA_MAX, Q_a + Q_b <= 1, b5 >= 1}
POOL_N = 85000     # @@POOL_N@@
POOL_VOL = 6.982594352970235e+37      # @@POOL_VOL@@  integer-point volume of D at A0 (computed in region_volume(), checked against the embedded value)


def inv_matrix(Mx):
    n = len(Mx)
    A = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(Mx)]
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(A[i][k]))
        A[k], A[p] = A[p], A[k]
        d = A[k][k]
        A[k] = [v / d for v in A[k]]
        for i in range(n):
            if i != k and A[i][k]:
                f = A[i][k]
                A[i] = [vi - f * vk for vi, vk in zip(A[i], A[k])]
    return [row[n:] for row in A]


def _pool_boxes():
    MA = [r[:6] for r in ELL_M[:6]]
    MB = [r[6:] for r in ELL_M[6:]]
    IA, IB = inv_matrix(MA), inv_matrix(MB)
    ra = [math.sqrt(THETA_MAX * IA[i][i]) for i in range(6)]
    rb = [math.sqrt(IB[i][i]) for i in range(6)]
    return ra, rb


def sample_coords(i, boxes):
    """the pool sampler: uniform integer points of D (rejection from bounding boxes), coprime content."""
    ra, rb = boxes
    rng = random.Random(POOL_TAG % i)
    while True:
        a = [rng.randint(-int(ra[k] * A0), int(ra[k] * A0)) for k in range(6)]
        qa, _ = quad_blocks([c / A0 for c in a], [0.0] * 6)
        if qa > THETA_MAX:
            continue
        s = math.sqrt(1 - qa)
        b = [rng.randint(-int(rb[k] * s * A0), int(rb[k] * s * A0)) for k in range(5)] + [rng.randint(1, max(1, int(rb[5] * s * A0)))]
        _, qb = quad_blocks([0.0] * 6, [c / A0 for c in b])
        if qb > 1 - qa:
            continue
        if math.gcd(gcd_ints(a), gcd_ints(b)) != 1:
            continue
        return a, b, qa


def slice_weight_norm():
    """E[(1-q)^3] for q = Q_a(a), a uniform in {Q_a <= THETA_MAX}: the pool draws a uniformly there and b uniformly
    in the a-slice, so a sample represents (1-q)^3 / E[(1-q)^3] points of D (uniform over D after weighting)."""
    m = 4000
    return sum(3 * ((i + .5) / m * THETA_MAX) ** 2 * (1 - (i + .5) / m * THETA_MAX) ** 3 for i in range(m)) * THETA_MAX / m / THETA_MAX ** 3


def region_volume():
    """integer points of D at A0: sum over a with Q_a <= THETA_MAX of |{b: Q_b <= 1-Q_a, b5>=1}| = vol_a * 3 int_0^T q^2 (1-q)^3 dq * vol_b / 2."""
    MA = [r[:6] for r in ELL_M[:6]]
    MB = [r[6:] for r in ELL_M[6:]]
    vol_a = vball(6) / math.sqrt(det(MA))
    vol_b = vball(6) / math.sqrt(det(MB))
    m = 4000
    integ = sum(3 * ((i + .5) / m * THETA_MAX) ** 2 * (1 - (i + .5) / m * THETA_MAX) ** 3 for i in range(m)) * THETA_MAX / m
    return vol_a * integ * vol_b / 2 * A0 ** 12


_UNIT_FORMS = None


def unit_forms(units):
    """Q_a(eps a) = a^T G_eps a with G_eps = U_eps^T M_a U_eps, U_eps the (integer) matrix of multiplication by eps."""
    global _UNIT_FORMS
    if _UNIT_FORMS is None:
        forms = []
        for u in units:
            cols = [(rmul(u, [1 if j == i else 0 for j in range(6)]) + [0] * 6)[:6] for i in range(6)]   # eps * y^i
            Uq = [[cols[j][i] for j in range(6)] for i in range(6)]                                      # Uq[i][j] = coeff i of eps*y^j
            MU = [[sum(ELL_M[i][k] * Uq[k][j] for k in range(6)) for j in range(6)] for i in range(6)]
            G = [[sum(Uq[k][i] * MU[k][j] for k in range(6)) for j in range(6)] for i in range(6)]
            forms.append([[G[i][j] + (G[j][i] if j != i else 0.0) for j in range(i + 1)] for i in range(6)])   # lower triangle, off-diagonals doubled
        _UNIT_FORMS = forms
    return _UNIT_FORMS


def canonical_a(a, units, qa):
    """a is canonical iff Q_a(a) <= Q_a(eps a) for every unit eps (mod sign); exact in the integers up to float rounding
    of the form (the pool flags were computed with the same real forms in canon.c)."""
    x = [c / A0 for c in a]
    thr = qa - 1e-12
    for G in unit_forms(units):
        q = 0.0
        for i in range(6):
            Gi = G[i]
            xi = x[i]
            t = 0.0
            for j in range(i + 1):
                t += Gi[j] * x[j]
            q += t * xi
        if q < thr:
            return False
    return True


def load_pool():
    """embedded canonical rows of the targeted pool: (idx, log2 nf, log2 ng, Q_a).  Non-canonical samples are not
    sieved and are not embedded; their number is POOL_N - count."""
    b85 = ""
    N = count = None
    for p in (1, 2, 3):
        fn = os.path.join(_DIR, "pool_canon_%d.txt" % p)
        if not os.path.exists(fn):
            continue
        lines = open(fn).read().splitlines()
        meta = dict(tok.split("=") for tok in lines[0].split() if "=" in tok)
        if p == 1:
            N, count = int(meta["N"]), int(meta["count"])
        b85 += "".join(lines[2:])
    raw = zlib.decompress(base64.b85decode(b85))
    recs = []
    i, idx = 0, -1
    while i < len(raw):
        delta, shift = 0, 0
        while True:
            byte = raw[i]
            i += 1
            delta |= (byte & 0x7F) << shift
            if not byte & 0x80:
                break
            shift += 7
        idx += delta
        nf = int.from_bytes(raw[i:i + 2], "big") / 256 + 600
        ng = int.from_bytes(raw[i + 2:i + 4], "big") / 256 + 400
        qa = int.from_bytes(raw[i + 4:i + 6], "big") / 65535.0 * THETA_MAX
        i += 6
        recs.append((idx, nf, ng, qa))
    assert len(recs) == count and N == POOL_N
    return recs, N


def verify_pool(samples, units):
    recs, N = load_pool()
    emb = {r[0]: r for r in recs}
    boxes = _pool_boxes()
    hi = min(samples, N)
    ncan = mism = flagmism = 0
    for i in range(hi):
        a, b, qa = sample_coords(i, boxes)
        flag = canonical_a(a, units, qa)
        e = emb.get(i)
        if flag != (e is not None):
            flagmism += 1
            continue
        if flag:
            ncan += 1
            nf = log2_exact(norm_f(a, b))
            ng = log2_exact(norm_g(a, b))
            if abs(e[1] - nf) > 1.0 / 256 + 1e-9 or abs(e[2] - ng) > 1.0 / 256 + 1e-9 or abs(e[3] - qa) > THETA_MAX / 65535 + 1e-9:
                mism += 1
    vol = region_volume()
    say(f"[2] targeted pool over D(theta_max={THETA_MAX}): N={N} samples, {len(recs)} canonical (f_c={len(recs)/N:.4f}); "
        f"seeds [0,{hi}) regenerated: {ncan} canonical, {flagmism} flag mismatches, {mism} norm/Q_a mismatches "
        f"{'ok' if mism == 0 and flagmism == 0 else 'FAIL'}; |D| = 2^{math.log2(vol):.4f} points (embedded 2^{math.log2(POOL_VOL):.4f})")
    return recs, N, vol, mism == 0 and flagmism == 0 and abs(math.log2(vol) - math.log2(POOL_VOL)) < 1e-6


# ------------------------------------------------------------------ [3] geometry ----------------
def det(A):
    n = len(A)
    M = [r[:] for r in A]
    d = 1.0
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k]))
        if abs(M[p][k]) < 1e-300:
            return 0.0
        if p != k:
            M[k], M[p] = M[p], M[k]
            d = -d
        d *= M[k][k]
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n):
                M[i][j] -= f * M[k][j]
    return d


def vball(d):
    return math.pi ** (d / 2) / math.gamma(d / 2 + 1)


def log2_li(log2x):
    x = log2x * LN2
    return log2x - math.log2(x) + math.log2(1 + 1 / x + 2 / x / x)


def mertens_sum(log2x):           # sum_{p <= x} 1/p
    return math.log(log2x * LN2) + 0.2615


def kappa(u, t):
    r = rho(u)
    return rho(u - t) / r if r > 0 else 0.0


def poisson_cdf(lam, k):
    s, term = 0.0, math.exp(-lam)
    for i in range(k + 1):
        s += term
        term *= lam / (i + 1)
    return s


def ratio_lookup(uf, ug, excess):
    """bilinear in (u_f, u_g) on log(N/U), linear in excess; clamped to the grid (edges are conservative here)."""
    ufs = sorted(set(k[0] for k in RATIO_GRID))
    ugs = sorted(set(k[1] for k in RATIO_GRID))

    def interp1(xs, x):
        x = min(max(x, xs[0]), xs[-1])
        for a, b in zip(xs, xs[1:]):
            if x <= b:
                return a, b, (x - a) / (b - a)
        return xs[-2], xs[-1], 1.0
    f0, f1, tf = interp1(ufs, uf)
    g0, g1, tg = interp1(ugs, ug)

    def at(f, g):
        d = RATIO_GRID[(f, g)]
        es = sorted(d)
        e0, e1, te = interp1(es, excess)
        return math.log(d[e0]) + (math.log(d[e1]) - math.log(d[e0])) * te
    v = (1 - tf) * ((1 - tg) * at(f0, g0) + tg * at(f0, g1)) + tf * ((1 - tg) * at(f1, g0) + tg * at(f1, g1))
    return math.exp(v)


# ------------------------------------------------------------------ [4] ARCH-E cost model -------
def lines_per_slice(la):
    """number of b5 lines of the half-slice: floor(sqrt(R^2 (M_b^-1)_55)) at scale A (R = 1 in box-normalised coordinates)"""
    MB = [r[6:] for r in ELL_M[6:]]
    IB = inv_matrix(MB)
    return max(1.0, math.floor(math.sqrt(IB[5][5]) * 2 ** la))


def c_tab(tab, log2h):
    if log2h >= tab[0][0]:
        return tab[0][1]
    for (h1, c1), (h2, c2) in zip(tab, tab[1:]):
        if log2h >= h2:
            return c2 + (c1 - c2) * (log2h - h2) / (h1 - h2)
    return tab[-1][1]


def evaluate(rows, N, vol, la, lb, k_f, k_g, lx2, theta, excess, scale_consts=1.0, pad_n=None):
    """rows: canonical pool rows (nf, ng, qa) at A0; N pool samples over the region of |D| = vol points.  Returns dict or None."""
    pad_n = PAD_N if pad_n is None else pad_n
    s72 = 72 * (la - math.log2(A0))
    s12 = 12 * (la - math.log2(A0))
    MB = [r[6:] for r in ELL_M[6:]]
    vol_b = vball(6) / math.sqrt(det(MB))
    nF = 0
    Y = 0.0
    setup_w = 0.0
    nf_mean = ng_mean = 0.0
    EW = slice_weight_norm()
    for (nf, ng, qa) in rows:
        if qa > theta:
            continue
        wt = (1 - qa) ** 3 / EW                                          # sample weight (see slice_weight_norm)
        nF += wt
        nfA, ngA = nf + s72, ng + s12
        Y += wt * proba(nfA, ALPHA_F, lb) * proba(ngA, ALPHA_G, lb)
        nf_mean += wt * nfA
        ng_mean += wt * ngA
        setup_w += wt / (vol_b * (1 - qa) ** 3 / 2 * (2 ** la) ** 6)   # wt / (points of this a-slice, b5>0)
    if nF <= 0:
        return None
    nf_mean /= nF
    ng_mean /= nF
    uf, ug = nf_mean / lb, ng_mean / lb
    scale = vol * (2 ** la / A0) ** 12 / N / 2                     # pool sample -> production points; /2 for a5 > 0
    log2_F = math.log2(nF * scale)
    # cofactorisation tier: primes in (2^lx2, 2^lb] are not sieved; a relation is kept iff each side has <= k such primes
    t3 = lx2 / lb

    def lam_top(u):
        m, s = 200, 0.0
        for i in range(m):
            t = t3 + (1 - t3) * (i + 0.5) / m
            s += kappa(u, t) * ((1 - t3) / m) / t
        return s
    lam_f, lam_g = lam_top(uf), lam_top(ug)
    keep = poisson_cdf(lam_f, k_f) * poisson_cdf(lam_g, k_g)
    log2_rels = math.log2(Y * scale * INV_ZETA * keep)
    log2_U = log2_li(lb)
    log2_need = log2_U + math.log2(excess)
    feasible = log2_rels >= log2_need + MARGIN
    # ---- cost per F-point (instructions) ----
    mean_slice_log2 = math.log2(nF / setup_w)
    lines_log2 = math.log2(lines_per_slice(la))
    mean_line_log2 = mean_slice_log2 - lines_log2
    upd = 0.0
    sparse = 0.0
    lp = LX0
    while lp < lx2:
        lp2 = min(lp + 0.25, lx2)
        dsum = mertens_sum(lp2) - mertens_sum(lp)
        if lp2 <= P_MID:                                   # line path: h = line points / p
            upd += NBAR * dsum * c_tab(UPD_TAB, mean_line_log2 - 0.5 * (lp + lp2))
        else:                                              # record path: h = slice points / p
            log2h = mean_slice_log2 - 0.5 * (lp + lp2)
            if log2h >= 0:
                upd += NBAR * dsum * c_tab(REC_TAB, log2h)
            else:
                dpairs = (NBAR / 2) * (2 ** log2_li(lp2) - 2 ** log2_li(lp))
                sparse += dpairs * C_COSET_SPARSE * (setup_w / nF)
        lp = lp2
    # line-path setups: 2 cosets per (pair below P_MID, line)
    lset = C_LSET * 2 * lines_per_slice(la) * (NBAR / 2) * 2 ** log2_li(min(lx2, P_MID)) * (setup_w / nF)
    pairs2 = (NBAR / 2) * 2 ** log2_li(lx2)
    setup = C_SET2 * pairs2 * (setup_w / nF)
    Sf, Sg = k_f * lb, k_g * lb
    surv = rho(max(nf_mean - Sf, 0) / lx2) * rho(max(ng_mean - Sg, 0) / lx2)
    # survivors: exact norms, trial division by the small primes, the recorded large-prime hits divided out, then full
    # ECM (primes <= B) on each cofactor, which is at most S_f / S_g bits by the threshold
    c_cof = C_NORM + C_TDIV + descent_core.ecm_instr(lb, min(nf_mean, Sf)) + descent_core.ecm_instr(lb, min(ng_mean, Sg))
    cof = surv * c_cof
    per_point = (upd + sparse + setup + lset + C_SCAN) * scale_consts + cof
    log2_relcol = log2_F + math.log2(per_point) - math.log2(UNIT)
    direct = [d for (ex, lu, sd, d, w, mc) in CONV_TABLE if abs(ex - excess) < 1e-9 and lu >= 17]
    ratio = (sum(direct) / len(direct)) if direct else ratio_lookup(uf, ug, excess)
    log2_N = log2_U + math.log2(ratio * pad_n)
    log2_la = math.log2(BW) + 2 * log2_N + math.log2(W_ROW * C_ENTRY / UNIT)
    return dict(la=la, lb=lb, k_f=k_f, k_g=k_g, lx2=lx2, theta=theta, excess=excess, feasible=feasible,
                log2_F=log2_F, log2_rels=log2_rels, log2_need=log2_need, uf=uf, ug=ug, nf=nf_mean, ng=ng_mean,
                per_point=per_point, upd=upd, sparse=sparse, lset=lset, setup=setup, cof=cof, surv=surv, c_cof=c_cof, keep=keep,
                lam_f=lam_f, lam_g=lam_g, log2_relcol=log2_relcol, ratio=ratio, log2_N=log2_N, log2_la=log2_la,
                mean_slice_log2=mean_slice_log2, nF=nF, pairs2=pairs2)


def total_log2(e, log2_desc):
    return math.log2(2 ** e["log2_relcol"] + 2 ** e["log2_la"] + 2 ** log2_desc + 2 ** LOG2_SM + 2 ** LOG2_FBC
                     + 2 ** LOG2_CANON + 2 ** LOG2_MISC + 2 ** LOG2_VERIFY)


# ------------------------------------------------------------------ [7] LA kernel ---------------
def check_la_kernel():
    lines = open(os.path.join(_DIR, "la_mac.s")).read().splitlines()
    on, body = False, []
    for ln in lines:
        s = ln.strip()
        if s.startswith("LBB0_2:"):
            on = True
            continue
        if on:
            if s.startswith("LBB") or s.startswith(";") and "bb." in s:
                break
            if s and not s.startswith(";") and not s.startswith("."):
                body.append(s.split()[0])
    branches = [m for m in body if m.startswith("b.") or m == "b"]
    return len(body), branches


# ------------------------------------------------------------------ [8] boot lattice ------------
def f6_mul(A, B):
    return [c % P for c in pmod_h(pmul([c % P for c in A], [c % P for c in B]))]


def f6_inv(A):
    def ptrimmod(p):
        return ptrim([c % P for c in p])

    def pdivmod(a, b):
        a = ptrimmod(a[:])
        b = ptrimmod(b[:])
        db = len(b) - 1
        inv = pow(b[-1], P - 2, P)
        q = [0] * max(1, len(a) - db)
        while len(a) - 1 >= db and any(a):
            d = len(a) - 1 - db
            c = a[-1] * inv % P
            q[d] = c
            for i in range(len(b)):
                a[d + i] = (a[d + i] - c * b[i]) % P
            a = ptrimmod(a)
        return q, a
    r0, r1 = ptrimmod(H[:]), ptrimmod(A[:])
    t0, t1 = [], [1]
    while any(r1):
        q, r = pdivmod(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, ptrimmod(padd(t0, [-c for c in pmul(q, t1)]))
    c = r0[0] if r0 else 1
    return [x * pow(c, P - 2, P) % P for x in t0]


def f12_inv(zu, zv):
    """(zu + zv theta)^-1 = (zu - zv theta)/N with N = zu^2 + U y zv^2 in F_{p^6} (theta^2 = -U y)."""
    nrm = padd(f6_mul(zu, zu), f6_mul([0, U], f6_mul(zv, zv)))
    ninv = f6_inv(nrm)
    return f6_mul(zu, ninv), [(-c) % P for c in f6_mul(zv, ninv)]


def boot_membership(x12, u, v):
    a, b = x12[:6], x12[6:]
    t = padd(f6_mul(a, v), [(-c) % P for c in f6_mul(b, u)])
    return all(c % P == 0 for c in t)


def det12(Bm):
    n = 12
    M = [row[:] for row in Bm]
    sign, prev = 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    M[k], M[i] = M[i], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def check_boot():
    bd = json.loads(zlib.decompress(base64.b85decode(open(os.path.join(_DIR, "boot_data.b85")).read())))
    assert bd.get("h", H) == H
    u, v = f12_inv(bd["z_u"], bd["z_v"])
    bad = sum(0 if boot_membership(bv, u, v) else 1 for bv in bd["basis"])
    d = det12(bd["basis"])
    ok = bad == 0 and abs(d) == P ** 6
    lens = sorted(math.log2(math.sqrt(sum(c * c for c in bv))) for bv in bd["basis"])
    say(f"[8] boot lattice L_z (det p^6): {12 - bad}/12 basis vectors satisfy the exact subfield condition; |det| == p^6: "
        f"{abs(d) == P ** 6} ({d.bit_length()} bits); LLL vector lengths 2^{lens[0]:.1f}..2^{lens[-1]:.1f} (p^(1/2) = 2^{math.log2(P)/2:.1f}); "
        f"measured mean log2 N_g by ball radius: " + ", ".join(f"x{m}: {q[3]:.0f}" for m, q in sorted(bd['quants'].items(), key=lambda kv: float(kv[0]))))
    return ok, bd


# ------------------------------------------------------------------ main ------------------------
def main():
    cfg = json.load(sys.stdin) if not sys.stdin.isatty() else {}
    samples = int(cfg.get("samples", 2000))
    seed = int(cfg.get("seed", 1))
    say(f"evidence: stnfs-aslice-v6 (bls12-381-fixed-g2-invert-g1); samples={samples} seed={seed}; 1 unit = {UNIT:.0f} instructions")
    ok0 = check_resultant()
    ok0b = check_tower()
    ok1 = check_rho()
    units = build_units()
    recs, N, vol, ok2 = verify_pool(samples, units)
    rows = [(nf, ng, qa) for (idx, nf, ng, qa) in recs]
    EW = slice_weight_norm()
    f_c = sum((1 - r[2]) ** 3 / EW for r in rows) / N                 # weighted canonical fraction over D
    qv = vball(12) / math.sqrt(det(ELL_M)) / 2 ** 12
    MA = [r[:6] for r in ELL_M[:6]]
    MB = [r[6:] for r in ELL_M[6:]]
    cross = max(abs(ELL_M[i][j]) for i in range(6) for j in range(6, 12))
    MBm = [r[6:] for r in ELL_M[6:]]
    IBm = inv_matrix(MBm)
    r5 = [math.sqrt(IBm[i][i]) * 2 ** OP["la"] for i in range(5)]
    L5 = 1.0
    for x in r5:
        L5 *= 2 * x + 1
    nlines = lines_per_slice(OP["la"])
    slice_pts = vball(6) / math.sqrt(det(MBm)) / 2 * 2 ** (6 * OP['la'])          # largest half-slice (Q_a = 0)
    nrec = slice_pts * NBAR * (mertens_sum(OP['lx2']) - mertens_sum(P_MID))
    nseg = math.ceil(nlines / max(1, math.floor(2 ** 62 / L5)))
    say(f"[3m] slice layout at A=2^{OP['la']}: largest half-slice 2^{math.log2(slice_pts):.2f} points, {nlines:.0f} lines (2^{math.log2(nlines):.2f}), "
        f"5-dim box per line 2^{math.log2(L5):.2f} bytes, slice 2^{math.log2(L5 * nlines):.2f} bytes (box factor {L5 * nlines / slice_pts:.1f}), "
        f"{nseg} segments of <= 2^62 bytes; record path (2^{P_MID:.0f} <= p <= 2^{OP['lx2']:.0f}): 2^{math.log2(nrec):.2f} records x 12 bytes = 2^{math.log2(nrec * 12):.2f} bytes")
    say(f"[3] ellipsoid: volume fraction {qv:.5f} (claimed {Q_FRAC:.5f}); a/b cross terms {cross:g} (sigma-invariant by construction); "
        f"canonical-a fraction over D (slice-weighted): f_c = {f_c:.4f}; Q_a <= theta={OP['theta']}: {sum(1 for r in rows if r[2] <= OP['theta'])} canonical pool points")

    # [4] cost model at the claimed operating point and its neighbourhood
    op = OP
    e0 = evaluate(rows, N, vol, op["la"], op["lb"], op["k_f"], op["k_g"], op["lx2"], op["theta"], op["excess"])
    say(f"[4] operating point A=2^{op['la']} B=2^{op['lb']} excess={op['excess']} theta={op['theta']} X2=2^{op['lx2']} k=({op['k_f']},{op['k_g']}):")
    say(f"    F-points 2^{e0['log2_F']:.3f} (mean slice 2^{e0['mean_slice_log2']:.1f} points); mean log2 N_f {e0['nf']:.1f}, N_g {e0['ng']:.1f} "
        f"(u_f {e0['uf']:.2f}, u_g {e0['ug']:.2f})")
    say(f"    relation orbits 2^{e0['log2_rels']:.3f} from {sum(1 for r in rows if r[2] <= op['theta'])} pool points (keep factor {e0['keep']:.3f}: unsieved primes per side Poisson({e0['lam_f']:.2f}), "
        f"Poisson({e0['lam_g']:.2f})) vs needed 2^{e0['log2_need']:.3f} + margin {MARGIN}: {'feasible' if e0['feasible'] else 'INFEASIBLE'}")
    say(f"    instructions per F-point {e0['per_point']:.1f} = updates (line path below 2^{P_MID:.0f}, record path above) {e0['upd']:.1f} + empty cosets {e0['sparse']:.1f} + record-path setups {e0['setup']:.1f} + line setups {e0['lset']:.2f} "
        f"(sigma-pairs <= X2: 2^{math.log2(e0['pairs2']):.1f}) + scan {C_SCAN} + cofactorisation {e0['cof']:.1f} "
        f"(survivors 2^{math.log2(e0['surv']):.1f} x 2^{math.log2(e0['c_cof']):.1f})")
    say(f"    relation collection 2^{e0['log2_relcol']:.3f} units; filtering N/U = {e0['ratio']:.4f} (mean of the embedded 2^17-2^18 runs at this profile and excess, W = {W_ROW:.0f}) x pad {PAD_N} -> N = 2^{e0['log2_N']:.2f}; "
        f"LA 2^{e0['log2_la']:.3f} units")
    # [5] filtering simulation at reduced scale, at the operating point's profile
    L = op["lb"]
    Us, res = filtsim.run(L, e0["uf"], e0["ug"], op["excess"], seed, 14, kmax=64, targets=(100, 150), verbose=False)
    r100 = res[100]
    say(f"[5] filtering simulation (production occurrence profile with the dense columns propagated through merge, scale 2^14 sampled columns, "
        f"seed {seed}): N/U = {1 / r100[0]:.4f} at full row weight {r100[3]:.0f} (max |coefficient| {r100[2]} <= 255); model uses {e0['ratio']:.4f} x pad {PAD_N}")
    if CONV_TABLE:
        say("    embedded runs at this profile, W = 100 (excess, log2 sampled columns, seed -> N/U, full row weight, max coef): " +
            "; ".join(f"({a},{b},{c}: {d:.4f},{w:.0f},{m})" for (a, b, c, d, w, m) in CONV_TABLE))
    cal = filtsim.calibrate_dgp21(13, seed, verbose=False)
    say(f"    DGP21 record calibration (F_{{p^6}} 521-bit, B=2^27, 0.896 relations per ideal; scale 2^13 here, 2^17 offline): ideals appearing "
        f"{cal['appearing']:.3f} (reported 0.87), rows kept by singleton removal {cal['purged_frac']:.3f} (reported 0.38), merge to density 100/150: "
        f"{cal['merge100']:.3f}/{cal['merge150']:.3f} of purged rows (reported 0.33/0.29); offline 2^17: 0.82, 0.61-0.65, 0.27-0.35 (proof.md H3)")

    # [6] descent
    Pm = descent_core.P()
    Pm.lB = op["lb"]
    Pm.alpha_f, Pm.alpha_g = ALPHA_F, ALPHA_G
    Pm.LMAX = 265
    Pm.BOOT_NORM = BOOT_NORM
    ld, bestb, detail = descent_fast.solve_fast(Pm)
    tot, lB1, own, sub, cand, nk = bestb
    log2_desc_units = ld - math.log2(UNIT)
    say(f"[6] descent DP (1 log): boot B1=2^{lB1}, 2^{math.log2(cand):.1f} candidates, boot 2^{math.log2(own):.1f} instr (x{Pm.PAD_NODE} pad), "
        f"tree 2^{math.log2(sub):.1f} instr; per log 2^{log2_desc_units:.2f} units; two logs 2^{log2_desc_units + 1:.2f}; "
        f"cap 2^{DESCENT_CAP:.0f} -> Markov P(abort) <= {2 ** (log2_desc_units + 1 - DESCENT_CAP):.4f}")
    for b in (70, 80, 100, 150, 220, 260):
        for side in ("f", "g"):
            key = (b - ((b - int(math.floor(op['lb']))) % 2), side)
            if key in detail:
                t, n, c, S, nf, ng, keep, how, nkids = detail[key]
                say(f"    node 2^{key[0]}-2^{key[0]+2} {side}: target 2^{n} via {how}, 2^{S:.0f} candidates, N_f {nf:.0f} N_g {ng:.0f} bits, "
                    f"keep {keep:.2f}, own 2^{math.log2(c):.1f}, children {nkids:.1f}, subtree 2^{math.log2(t):.1f} instr")

    # [7] LA kernel
    la_n, la_br = check_la_kernel()
    say(f"[7] la_mac.s inner loop: {la_n} instructions per entry, branches {la_br} ({'ok' if la_n == 33 and la_br == ['b.ne'] else 'MISMATCH'})")

    # [8] boot lattice
    ok8, bd = check_boot()

    # ledger
    tot = total_log2(e0, DESCENT_CAP)
    say(f"ledger (log2 units): relcol {e0['log2_relcol']:.3f}, LA {e0['log2_la']:.3f}, descent cap {DESCENT_CAP:.0f}, SM {LOG2_SM:.0f}, "
        f"FB {LOG2_FBC:.0f}, canonical list {LOG2_CANON:.0f}, misc {LOG2_MISC:.0f}, verify {LOG2_VERIFY:.0f}")
    say(f"total_log2W_full_ledger = {tot:.6f}  (claimed {CLAIM:.6f}: {'ok' if tot <= CLAIM else 'EXCEEDS CLAIM'})")
    # sensitivities
    for label, kw in (("sieve constants x1.5", dict(scale_consts=1.5)), ("matrix pad 4 instead of 2", dict(pad_n=4.0))):
        e = evaluate(rows, N, vol, op["la"], op["lb"], op["k_f"], op["k_g"], op["lx2"], op["theta"], op["excess"], **kw)
        say(f"    sensitivity [{label}]: relcol {e['log2_relcol']:.3f}, LA {e['log2_la']:.3f}, total {total_log2(e, DESCENT_CAP):.3f}")
    checks = dict(resultant=ok0, tower=ok0b, rho=ok1, pool=ok2, boot=ok8, feasible=e0["feasible"], within_claim=tot <= CLAIM)
    print(json.dumps({"checks_ok": all(checks.values()), "checks": checks, "total_log2W_full_ledger": round(tot, 6),
                      "claimed": CLAIM, "operating_point": OP, "f_c": round(f_c, 5),
                      "log2_relcol": round(e0["log2_relcol"], 4), "log2_la": round(e0["log2_la"], 4),
                      "descent_two_logs_log2_units": round(log2_desc_units + 1, 3), "report": _REPORT}, indent=1), flush=True)


BOOT_NORM = 2725.0   # @@BOOT_NORM@@  mean log2 N_g of boot candidates at the useful radius (measured, part [8])

if __name__ == "__main__":
    main()
