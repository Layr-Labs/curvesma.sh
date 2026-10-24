"""Evidence for the STNFS pairing-inversion claim (bls12-381-fixed-g2-invert-g1). Stdlib only.

The attack solves the discrete logarithm of the challenge target z = e(P0, Q)^a in the order-r
subgroup GT of F_{p^12}^* by a Special Tower Number Field Sieve, then returns P = [a]P0.
This script reproduces, in exact integer arithmetic where it matters, the pipeline behind the
claimed total work (1 unit = one dense Fp12 multiplication = 20,000 scalar AArch64
instructions, the track's fixed reference price):

  (0) Polynomial-selection validity: Res_y(h, Res_x(f, g)) is computed exactly and checked to
      be 0 mod p (the two NFS sides share a root mod p, i.e. generate the target field F_{p^12}).
  (1) Dickman rho by Marsaglia's exact piecewise power series; self-tested.
  (2) Held-out verification: eval_inside_1.txt, eval_inside_2.txt carry the inside-ellipsoid
      subset of a held-out pool of N per-sample-seeded candidates (the ellipsoid was fitted on
      a disjoint sample set). The run recomputes the exact norms for the first `samples`
      candidate indices from their seeds, re-checks ellipsoid membership and confirms the
      recomputed inside-set matches the embedded one.
  (3) Slab geometry of the embedded ellipsoid, computed from its matrix: volume fraction and
      the mean volume of the 3-dim slices (this fixes the number of slabs, which the sieve
      model charges per (small ideal, slab) setup).
  (4) Joint (A, B) optimizer under the ARCH-D sieve model. Every point of the sieved
      ellipsoid is charged (yield is credited only from its in-box part), the slab count is
      the true one from (3), and the small/large-ideal threshold is optimised. Constants are
      retired-instruction counts of correctness-checked kernels shipped in this package
      (C_HIT, C_ENTRY) or 1.25x pads of static counts (C_UPD, C_SCAN, C_SET, C_REP); see
      KERNELS.md. Prints box baseline, ellipsoid optimum, relcol/LA split and sensitivities.
  (5) Toy-scale exact-smoothness check of the rho model (real random integers, exact factoring).
  (6) Exact degree-1 prime-ideal root statistics of the real (h, f, g) at small primes.
  (7) Individual-log descent estimate (booting + special-q tree) with this package's rho and the
      sieve's own per-point cost, and the metered-cap success bound.
  (8) Instruction count of the shipped LA inner loop (la_mac.s) re-derived from the assembly.
  Full ledger: relcol + LA + descent cap + Schirokauer maps + factor-base construction + misc.

alpha_f, alpha_g and 1/zeta_Kh(2) are documented inputs (Guillevic-Singh alpha code, offline);
everything else is recomputed here. Toy parts are evidence at toy scale only.
Usage: python evidence.py  (reads {samples, seed} as JSON on stdin; prints JSON on stdout)
"""
import base64
import json
import math
import os
import random
import sys
import zlib
from decimal import Decimal, getcontext

# ---- fixed parameters (provenance: Guillevic-Singh TNFS polyselect, sage, offline) ----
P = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
U = 15132376222941642752          # |BLS12-381 seed u|; g = x^2 + u*y
H = [-1, -1, 1, 1, 0, 0, 1]       # h(y) = y^6+y^3+y^2-y-1, monic, index = degree
# f(x,y) = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (-y^3-y^2+y+1): (x-degree, y-poly)
FTERMS = [(12, [1]), (10, [0, -2]), (6, [0, 0, 0, 2]), (2, [0, 0, 0, 0, 0, 1]), (0, [1, 1, -1, -1])]
ALPHA_F = 4.7208                  # log2, GS alpha_tnfs_2d, 1000 samples
ALPHA_G = 0.8445
INV_ZETA = 0.9168940567189809     # 1/zeta_Kh(2), Kh = Q[y]/h
AUT = 2                           # x -> -x automorphism (Frob^6 conjugation), kappa=2
# ---- H4 sieve-cost model (ARCH-D, distributed slab buckets; KERNELS.md) --------------
# The region (fitted ellipsoid, volume fraction q) is cut into slabs by fixing 9 of the 12
# coordinates; each slab holds a bucket array over the bounding box of its 3-dim slice. Small
# ideals (p <= X) are sieved per slab (one setup per (ideal, slab), then stride updates);
# large ideals walk their lattice once over the whole ellipsoid and each hit is canonicalised
# and routed to the owner of its slab. Constants are instructions per event:
C_UPD = 17.5                      # small-ideal stride update + in-slice check: 14 static x 1.25
C_REP = 7.5                       # per sieved point: threshold/report: 6 static x 1.25
C_SCAN = 12.5                     # per scanned slice-bbox point: 10 static x 1.25
C_SET = 176.0                     # per (small ideal, slab) setup: 141 static x 1.25
#   1.25 pads the static-to-retired gap measured on the walk kernel (122.4 retired vs 103 static)
NBAR = 1.947                      # mean ideals per prime (measured, roots to 2^16)
S_CONST = -2.60                   # Mertens-type constant of S(B) (measured: S(2^16) = 2.5481)
C_ENTRY = 35.0                    # instructions per nonzero matrix entry per mat-vec: the shipped
                                  # sign-split kernel (la_mac.c/.s) runs 33 per entry (branch-free
                                  # loop, re-counted in part [8]) + ~100 per row (2 loop
                                  # prologues/epilogues + la_finish.c reduction) / WT=200 -> 33.5;
                                  # tested exact vs bigints (la_test.c)
BW = 3.0                          # block Wiedemann: ~3N single-vector products
WT = 200.0                        # row weight after filtering
FLT = 45.0                        # filtering factor (matrix/FB, record-calibrated)
EXCESS = 1.1                      # relations over factor base
MARGIN = 0.25                     # required log2 yield margin over EXCESS*FB (~2.8 sigma of the
                                  # held-out yield estimate's Monte-Carlo spread)
UNIT = 20000.0    # this track: 1 unit = one dense Fp12 mul = 20,000 AArch64 instructions
DESCENT_CAP = 106.0   # metered cap on both individual-log descents together (units); part [7]
LOG2_SM = 87.0        # Schirokauer maps: <= 2^37 instr per relation (exponent <= 72*255 bits in
                      # O_f/r^2, ~2^22 instr per product) x <= 2^63.4 relations = 2^86.1 units
LOG2_FBC = 72.0       # factor-base construction (2^61.9 primes x ~2^20 root finding) + per-ideal
                      # 12-dim lattice reduction (2^63.0 ideals x <= 2^22 instr) = 2^70.9 units
LOG2_MISC = 60.0      # filtering (sort/dedupe 2^63.4 relations), relation/matrix storage writes,
                      # ellipsoid fit, polynomial/alpha search, challenge pairing; each <= 2^56
LOG2_VERIFY = 11.0    # final pairing, subgroup checks, [a]P0, serialization
ETA = 6
LN2 = math.log(2)
GAMMA1 = 1.0 - 0.5772156649015329
RHO_FLOOR = 125.9618             # log2 of sqrt(pi*r/(2*12)): orbit-12 rho STEP-count floor;
                                  # per-step cost is 0.77-1.15 units (cyclotomic squaring vs dense
                                  # mul); the public frontier claims a total of 125.73

A0 = 1 << 10

# The harness parses stdout as a single JSON object, so all human-readable lines go to stderr
# and are also collected into the JSON "report" field. Only the final JSON object hits stdout.
_REPORT = []


def say(*args):
    line = " ".join(str(a) for a in args)
    _REPORT.append(line)
    print(line, file=sys.stderr)


Q_FRAC = 1.0 / 32                 # ellipsoid volume fraction
A_GRID = [9.3 + 0.05 * i for i in range(13)]     # log2 A (focused window; optimum verified
                                             # interior under the full 29x186 research grid)
B_GRID = [64.0 + 0.25 * i for i in range(22)]    # log2 B (focused window; optimum interior
                                             # under the extended-to-76.25 research grid)
C12 = math.pi ** 6 / 720          # volume of the unit 12-ball
# ---- Z[y] arithmetic, index = degree ------------------------------------------------


def ptrim(p):
    while p and p[-1] == 0:
        p.pop()
    return p


def pmul(p, q):
    if not p or not q:
        return []
    out = [0] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        if pi:
            for j, qj in enumerate(q):
                out[i + j] += pi * qj
    return ptrim(out)


def padd(p, q):
    n = max(len(p), len(q))
    out = [0] * n
    for i, pi in enumerate(p):
        out[i] += pi
    for i, qi in enumerate(q):
        out[i] += qi
    return ptrim(out)


def pscale(p, c):
    return ptrim([c * x for x in p])


def pmod_h(p):
    """Reduce modulo the monic h (deg 6)."""
    d = len(H) - 1
    while len(p) - 1 >= d and p:
        shift = len(p) - 1 - d
        f = p[-1]
        if f:
            for i in range(d + 1):
                p[shift + i] -= f * H[i]
        ptrim(p)
    return p


def rmul(p, q):
    """Multiplication in Z[y]/h with deg<6 representatives."""
    return pmod_h(pmul(p, q))


# ---- companion matrix and resultants -----------------------------------------------

# companion matrix of h: M[i][i+1] = 1 (i<5), M[5][j] = -H[j]
MAT = [[0] * 6 for _ in range(6)]
for _i in range(5):
    MAT[_i][_i + 1] = 1
for _j in range(6):
    MAT[5][_j] = -H[_j]


def matmul(X, Y):
    return [[sum(X[i][k] * Y[k][j] for k in range(6)) for j in range(6)] for i in range(6)]


def matpoly_eval(Q):
    """Q(M) for Q a deg<=5 integer poly (list, index=degree), via Horner."""
    X = [[0] * 6 for _ in range(6)]
    for d in range(5, -1, -1):
        # X = X*M + q_d*I
        if any(any(r) for r in X):
            X = matmul(X, MAT)
        qd = Q[d] if d < len(Q) else 0
        if qd:
            for i in range(6):
                X[i][i] += qd
    return X


def det6(X):
    """Exact determinant of a 6x6 integer matrix (fraction-free Bareiss)."""
    n = 6
    M = [row[:] for row in X]
    sign = 1
    prev = 1
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
    """Res_y(h, Q) for deg Q <= 5: det of Q(companion matrix) (h monic)."""
    return det6(matpoly_eval(Q))


def norm_f(a, b):
    """N_f = |Res_y(h, Res_x(f, a-bx))|; Res_x(f, a-bx) = sum_i f_i(y) a^i b^(12-i)
    computed exactly in Z[y]/h."""
    apow = [[1]]
    bpow = [[1]]
    for _ in range(12):
        apow.append(rmul(apow[-1], a))
        bpow.append(rmul(bpow[-1], b))
    acc = []
    for (xd, fy) in FTERMS:
        term = rmul(pmod_h(fy[:]), rmul(apow[xd], bpow[12 - xd]))
        acc = padd(acc, term)
    return abs(res_y_h(pmod_h(acc)))


def norm_g(a, b):
    """N_g = |Res_y(h, Res_x(g, a-bx))|; Res_x(g, a-bx) = a^2 + u*y*b^2."""
    a2 = rmul(a, a)
    b2 = rmul(b, b)
    term = rmul([0, U], b2)
    return abs(res_y_h(padd(a2, term)))


def gcd_ints(xs):
    g = 0
    for x in xs:
        g = math.gcd(g, abs(x))
    return g


def log2_exact(n):
    """log2(n) as float, accurate to ~1e-15 relative."""
    bl = n.bit_length()
    if bl <= 53:
        return math.log2(n)
    lead = n >> (bl - 53)
    return (bl - 53) + math.log2(lead)


def check_resultant():
    # divide f by g = x^2 + U*y (monic in x) over Z[y]: remainder R0 + R1 x.
    import copy
    fdict = {xd: fy[:] for (xd, fy) in FTERMS}
    gdict = {2: [1], 0: [0, U]}
    f = copy.deepcopy(fdict)
    while True:
        d = max(f)
        if d < 2:
            break
        lc = f[d]
        if lc:
            for j, gj in gdict.items():
                t = pmul(lc, gj)
                key = d - 2 + j
                f[key] = pmod_h(padd(f.get(key, []), pscale(t, -1)))
        del f[d]
    R1 = f.get(1, [])
    R0 = f.get(0, [])
    # Res_x(f, g) = Res_x(R0 + R1 x, x^2 + U y) = R0^2 + U y R1^2   (roots sum 0, prod U y)
    resx = padd(rmul(pmod_h(R0[:]), pmod_h(R0[:])), rmul([0, U], rmul(pmod_h(R1[:]), pmod_h(R1[:]))))
    val = res_y_h(pmod_h(resx))
    say(f"[0] Res_y(h, Res_x(f,g)) mod p == 0 ? {val % P == 0}  (resultant has {val.bit_length()} bits)")
    return val % P == 0


# Yield-weighted covariance ellipsoid (vol. fraction 1/32 of the box); inside iff x^T ELL_M x <= 1,
# x = (a0..a5,b0..b5)/A.  Fitted on 500k exact-norm samples (seeds 1,2,3,25,50), disjoint from the eval pool.
ELL_M = [
    [0.5560579791, -0.1423773008, 0.1521544956, -0.008625772368, -0.02973055294, 0.1003647259, 0.004511129827, -0.00227398145, 0.003383957204, 0.008243079155, -0.004030328247, 0.005972241557],
    [-0.1423773008, 0.4621393439, -0.09370194495, 0.08404246835, -0.002881408129, -0.07470611654, 0.007156692057, 0.008003986377, 0.002261093061, 0.006069441532, 0.001915111497, -0.004958311802],
    [0.1521544956, -0.09370194495, 0.4566787804, -0.0618933189, 0.04300443037, -0.02438387976, 0.0001938634248, 0.006766577467, -0.003077441838, 0.002771500882, -0.001683319262, -0.0001330399408],
    [-0.008625772368, 0.08404246835, -0.0618933189, 0.43327467, -0.002973151503, -0.02621867505, 0.001268798138, 0.0002935839084, 0.004432853397, 0.001128585067, 0.0004112728192, -0.003456690784],
    [-0.02973055294, -0.002881408129, 0.04300443037, -0.002973151503, 0.4293307361, 0.0225005934, 0.003127443711, -0.0003115580231, 3.486702869e-05, 0.002091222687, -0.001536013733, -0.004107703974],
    [0.1003647259, -0.07470611654, -0.02438387976, -0.02621867505, 0.0225005934, 0.4770210492, -0.005000558504, -5.286617558e-05, -0.002085356467, 0.0009341353118, 0.002839858152, 0.00412506036],
    [0.004511129827, 0.007156692057, 0.0001938634248, 0.001268798138, 0.003127443711, -0.005000558504, 0.6126805356, -0.1185127897, 0.1698813919, 0.002268118173, -0.04962560224, 0.1447363527],
    [-0.00227398145, 0.008003986377, 0.006766577467, 0.0002935839084, -0.0003115580231, -5.286617558e-05, -0.1185127897, 0.5025909011, -0.07234540809, 0.08816798656, -0.009979898828, -0.08509656347],
    [0.003383957204, 0.002261093061, -0.003077441838, 0.004432853397, 3.486702869e-05, -0.002085356467, 0.1698813919, -0.07234540809, 0.5012616259, -0.04865565647, 0.05348257839, -0.04265795094],
    [0.008243079155, 0.006069441532, 0.002771500882, 0.001128585067, 0.002091222687, 0.0009341353118, 0.002268118173, 0.08816798656, -0.04865565647, 0.4835948843, 0.01282656811, -0.04081343341],
    [-0.004030328247, 0.001915111497, -0.001683319262, 0.0004112728192, -0.001536013733, 0.002839858152, -0.04962560224, -0.009979898828, 0.05348257839, 0.01282656811, 0.4952216288, 0.04648801675],
    [0.005972241557, -0.004958311802, -0.0001330399408, -0.003456690784, -0.004107703974, 0.00412506036, 0.1447363527, -0.08509656347, -0.04265795094, -0.04081343341, 0.04648801675, 0.5731120139],
]

# ---- Dickman rho: exact piecewise power series (Marsaglia) -------------------------
# On (k-1, k], rho(k - xi) = sum_m c^(k)_m xi^m, with (k - xi) d/dxi rho = rho_{k-1}(xi):
#   c^(k)_{m+1} = (c^(k-1)_m + m c^(k)_m) / (k (m+1)),  c^(k)_0 = c^(k-1)_0 - sum_{m>=1} c^(k)_m.
# c_0 cancels catastrophically for large k, so the coefficients are built once in high-precision
# decimal (then stored as floats). This is exact up to float rounding for all u we need, with no
# large-u clamping (the old table-integration clamped rho to 0 above u~11, biasing the f-side tail).
_RHO_KMAX = 32       # rho(32) ~ 1e-49; covers the small-B regime where slab volume >= B
_RHO_NT = 220
_RHO_STEP = 2.0 ** -10


def _build_rho_coefs():
    getcontext().prec = 150
    zero, one = Decimal(0), Decimal(1)
    coefs = [None, [one] + [zero] * (_RHO_NT - 1)]          # k = 1: rho == 1 on [0, 1]
    for k in range(2, _RHO_KMAX + 1):
        prev = coefs[k - 1]
        c = [zero] * _RHO_NT
        for m in range(_RHO_NT - 1):
            c[m + 1] = (prev[m] + m * c[m]) / (k * (m + 1))
        c[0] = prev[0] - sum(c[1:])
        coefs.append(c)
    return [None] + [[float(x) for x in c] for c in coefs[1:]]


_RHO_COEFS = _build_rho_coefs()


def rho_series(u):
    """Exact Marsaglia series value of rho(u) (slow; used to build the table and self-test)."""
    if u <= 1.0:
        return 1.0
    k = math.ceil(u)
    if k > _RHO_KMAX or k < 1:
        return 0.0
    xi = k - u
    s = 0.0
    for ci in reversed(_RHO_COEFS[k]):
        s = s * xi + ci
    return s if s > 0.0 else 0.0


# dense table for the optimizer's hot path: O(1) linear interpolation, built from the exact
# series once at import (~0.1 s). Max relative error vs the series is ~2e-6 on (1, 18).
_RHO_NMAX = int((_RHO_KMAX - 0.001) / _RHO_STEP)
_RHO_TAB = [1.0 if i * _RHO_STEP <= 1.0 else rho_series(i * _RHO_STEP) for i in range(_RHO_NMAX + 2)]


def rho(u):
    if u <= 1.0:
        return 1.0
    x = u / _RHO_STEP
    if x >= _RHO_NMAX:
        return 0.0
    i = int(x)
    f = x - i
    return _RHO_TAB[i] * (1.0 - f) + _RHO_TAB[i + 1] * f


def proba(log2n, alpha, log2b):
    """Smoothness probability of an integer of size 2^log2n at bound 2^log2b, alpha-corrected,
    with the second-order term rho(u) + (1-gamma) rho(u-1)/ln N (Murphy/Dickman)."""
    u = (log2n + alpha) / log2b
    ln_n = log2n * LN2
    return rho(u) + GAMMA1 * rho(u - 1.0) / ln_n if ln_n > 0 else rho(u)


def check_rho():
    ok = True
    for u, want in ((1.5, 1 - math.log(1.5)), (2.0, 1 - math.log(2)),
                    (3.0, 0.048608388291131), (5.0, 3.547247e-04), (8.0, 3.232069e-08),
                    (10.0, 2.770172e-11), (12.0, 1.419713e-14), (15.0, 7.589908e-20)):
        got = rho(u)
        good = abs(got - want) < 1e-3 * want + 1e-300
        ok = ok and good
        say(f"[1] rho({u:g}) = {got:.6e}  expected {want:.6e}  {'ok' if good else 'MISMATCH'}")
    return ok


# ---- yield-weighted ellipsoid (fitted offline, embedded) ---------------------------
# ELL_M injected below from _ellm.py


def quad(x):
    """x^T ELL_M x, fixed summation order (identical in generation and audit)."""
    s = 0.0
    for i in range(12):
        row = ELL_M[i]
        t = 0.0
        for j in range(12):
            t += row[j] * x[j]
        s += t * x[i]
    return s


def sample_one(i):
    """Held-out candidate i: its own RNG seeded by a string (CPython seeds str/bytes via SHA-512,
    stable across versions/platforms). Returns (12 coeffs / A0, log2 N_f, log2 N_g)."""
    rng = random.Random("bls12-381-stnfs-eval:%d" % i)
    while True:
        a = [rng.randint(-A0, A0) for _ in range(6)]
        b = [rng.randint(-A0, A0) for _ in range(5)] + [rng.randint(1, A0)]
        if math.gcd(gcd_ints(a), gcd_ints(b)) != 1:
            continue
        nf = norm_f(a, b)
        ng = norm_g(a, b)
        if nf == 0 or ng == 0:
            continue
        return [x / A0 for x in a] + [x / A0 for x in b], log2_exact(nf), log2_exact(ng)


# ---- embedded held-out inside-set ---------------------------------------------------


def _data_dir():
    for d in ("/snapshot/files", os.path.dirname(os.path.abspath(__file__)), "."):
        if os.path.exists(os.path.join(d, "eval_inside_1.txt")):
            return d
    raise SystemExit("eval_inside_*.txt not found")


def load_inside():
    """Decode the embedded inside-set: list of (index, log2 N_f, log2 N_g), plus (N, count)."""
    d = _data_dir()
    b85, N, count = "", None, None
    for p in (1, 2):
        lines = open(os.path.join(d, "eval_inside_%d.txt" % p)).read().splitlines()
        meta = dict(tok.split("=") for tok in lines[0].split() if "=" in tok)
        if p == 1:
            N, count = int(meta["N"]), int(meta["count"])
        b85 += "".join(lines[3:])           # skip the 3-line header
    raw = zlib.decompress(base64.b85decode(b85))
    recs = []
    i, idx = 0, -1
    while i < len(raw):
        delta, shift = 0, 0
        while True:
            byte = raw[i]; i += 1
            delta |= (byte & 0x7F) << shift
            if not byte & 0x80:
                break
            shift += 7
        idx += delta
        nf = int.from_bytes(raw[i:i + 2], "big") / 256 + 600; i += 2
        ng = int.from_bytes(raw[i:i + 2], "big") / 256 + 400; i += 2
        recs.append((idx, nf, ng))
    return recs, N, count


# ---- parts 3-4: the optimizer -------------------------------------------------------


def log2_li(log2x):
    lnx = log2x * LN2
    return log2x - math.log2(lnx) + math.log2(1 + 1 / lnx + 2 / lnx ** 2 + 6 / lnx ** 3)


def grid(NF, NG, mass_in, qfrac, excess=EXCESS):
    """Feasible (A, B) grid rows for a region; the yield E[P] sum is done once per (A, B)."""
    rows = []
    n = len(NF)
    for la in A_GRID:
        sf = [v + ETA * 12 * (la - 10.0) for v in NF]
        sg = [v + ETA * 2 * (la - 10.0) for v in NG]
        log2V0 = 12 * math.log2(2 * 2.0 ** la + 1) - 1
        for lb in B_GRID:
            ep = 0.0
            for fv, gv in zip(sf, sg):
                ep += proba(fv, ALPHA_F, lb) * proba(gv, ALPHA_G, lb)
            ep /= n
            if ep <= 0:
                continue
            log2_rels = log2V0 + math.log2(mass_in) + math.log2(INV_ZETA) + math.log2(ep)
            log2_FB = 1 + log2_li(lb)
            if log2_rels < math.log2(excess) + log2_FB + MARGIN:
                continue
            rows.append((la, lb, log2V0, log2_FB, mass_in, qfrac))
    return rows


def s_model(log2b):
    """Expected FB-ideal updates per candidate at bound 2^log2b: NBAR*(ln ln B + M)
    + S_CONST, with the small-prime ideal-deficit constant measured exactly to
    2^16 (roots.py; validated against direct congruence tests at toy scale)."""
    return max(0.0, NBAR * (math.log(log2b * LN2) + 0.2615) + S_CONST)


C_HIT = 150.0    # instr per large-ideal hit, RETIRED (kernel archd_hit.c, KERNELS.md): 12-dim
                # Schnorr-Euchner walk carrying 7 ambient state words, (a,b) -> canonical
                # orbit representative under +-1 and the x -> -x automorphism, slab-descriptor
                # lookup, and send to the slab owner + apply there: 140.1 measured; + one
                # intermediate routing hop (~7, the measured apply cost) -> 147; claimed 150.
                # (Single shared memory would measure 131.1, but 2^122 B is not 64-bit
                # addressable, so the distributed variant is the one charged.)
_BALLBOX3 = 1.910                 # 2^3 / V_3(unit ball): bbox waste of 3-dim ellipsoid slices


def relcol_instr(la, lb, log2V0, qfrac, vbar):
    """ARCH-D relation-collection instructions (before the automorphism factor).
    E = qfrac*V0 points are sieved (the whole ellipsoid, although yield is credited only from
    its in-box part); #slabs = E / (vbar*A^3) with vbar the true mean slice volume (part [3]);
    small ideals p <= X = 2^lx pay one setup per (ideal, slab) plus stride updates, larger
    ideals pay C_HIT per hit; lx is optimised. Returns (instructions, lx)."""
    m = 2.0 ** la
    E = qfrac * 2 ** log2V0
    nslabs = E / (vbar * m ** 3)
    S_B = s_model(lb)
    best = None
    for t in range(40, int(4 * lb) + 1):
        lx = t / 4.0
        S_X = s_model(lx)
        c = (nslabs * NBAR * 2 ** lx / (lx * LN2) * C_SET + E * S_X * C_UPD
             + E * max(0.0, S_B - S_X) * C_HIT)
        if best is None or c < best[0]:
            best = (c, lx)
    return best[0] + E * C_REP + _BALLBOX3 * E * C_SCAN, best[1]


def best_cost(rows, vbar, flt=FLT):
    """Minimum log2 W over feasible rows: (log2W, logA, logB, relcol, LA, lx) in units."""
    best = None
    for la, lb, log2V0, log2_FB, mass_in, qfrac in rows:
        instr, lx = relcol_instr(la, lb, log2V0, qfrac, vbar)
        log2_Trc = math.log2(instr) - math.log2(AUT)
        log2_N = log2_FB - math.log2(flt * AUT)
        log2_Tla = math.log2(BW * WT * C_ENTRY) + 2 * log2_N
        log2_W = math.log2(2 ** (log2_Trc - 40) + 2 ** (log2_Tla - 40)) + 40 - math.log2(UNIT)
        cand = (log2_W, la, lb, log2_Trc - math.log2(UNIT), log2_Tla - math.log2(UNIT), lx)
        if best is None or cand[0] < best[0]:
            best = cand
    return best


# ---- part 3: slab geometry of the embedded ellipsoid ----------------------------------


def _det(A):
    A = [r[:] for r in A]
    n = len(A)
    d = 1.0
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(A[i][k]))
        if p != k:
            A[k], A[p] = A[p], A[k]
            d = -d
        d *= A[k][k]
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            for j in range(k, n):
                A[i][j] -= f * A[k][j]
    return d


def _vball(d):
    return math.pi ** (d / 2) / math.gamma(d / 2 + 1)


def slab_geometry():
    """Volume fraction of {x^T M x <= 1} in [-1,1]^12, and the best choice of 3 slab coordinates
    T with its mean slice volume vol(E)/vol(proj E) = (V12/V9)/sqrt(det M_TT), in units of A^3."""
    qv = _vball(12) / math.sqrt(_det(ELL_M)) / 2 ** 12
    best = None
    import itertools
    for T in itertools.combinations(range(12), 3):
        vbar = (_vball(12) / _vball(9)) / math.sqrt(_det([[ELL_M[i][j] for j in T] for i in T]))
        if best is None or vbar > best[0]:
            best = (vbar, T)
    return qv, best[0], best[1]


# ---- part 5: toy exact-smoothness validation ----------------------------------------


def toy_smoothness(trials, seed, bits=40, log2b=13.0):
    rng = random.Random(seed)
    lim = 1 << int(log2b)
    sieve = bytearray([1]) * lim
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(lim ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    primes = [i for i in range(lim) if sieve[i]]
    nums = [rng.getrandbits(bits) | (1 << (bits - 1)) for _ in range(trials)]
    logs = [log2_exact(v) for v in nums]
    smooth = []
    for v in nums:
        m = v
        for pr in primes:
            if pr * pr > m:
                break
            while m % pr == 0:
                m //= pr
        smooth.append(m == 1 or m < lim)
    order = sorted(range(trials), key=lambda i: logs[i])
    ter = trials // 3
    say(f"[5] toy exact smoothness: {bits}-bit ints, B=2^{log2b:.0f}, {trials} samples")
    say(f"    overall smooth rate: {sum(smooth)/trials:.4f}  rho-model: {proba(bits, 0.0, log2b):.4f}")
    for t in range(3):
        idx = order[t * ter:(t + 1) * ter]
        rate = sum(smooth[i] for i in idx) / len(idx)
        mlog = sum(logs[i] for i in idx) / len(idx)
        say(f"    tercile {t+1} (mean log2N={mlog:.1f}): rate={rate:.4f}  rho-model={proba(mlog, 0.0, log2b):.4f}")
    small = sum(smooth[i] for i in order[:ter]) / ter
    large = sum(smooth[i] for i in order[2 * ter:]) / ter
    if large > 0:
        say(f"    small/large tercile yield ratio: {small/large:.1f}x (concentration mechanism)")
    return True


def roots_smallprimes(pmax):
    """[6] Exact degree-1 prime-ideal statistics of the real (h, f, g) for p <= pmax
    (small primes only: cost is O(pi(pmax)*pmax)). Returns (#primes, mean ideals/prime,
    S(pmax)). Provenance of the full-scale constants: research/nfs/h4test/roots.py to
    2^16 gives NBAR = 1.947, S(2^16) = 2.5481; direct congruence tests on real toy
    candidates reproduce the S model (1.95-1.96 measured vs 1.98 model at B = 2^12)."""
    lim = 1 << 12
    sieve = bytearray([1]) * lim
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(lim ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    tot_n = 0
    tot_s = 0.0
    n = 0
    for p in range(2, lim):
        if not sieve[p] or p > pmax:
            continue
        # roots of h mod p by trial
        iotas = [x for x in range(p) if poly_eval_mod(H, x, p) == 0]
        if len(iotas) == p:
            continue
        npr = 0
        for iota in iotas:
            # f(x, iota) mod p as dense poly in x (degrees 12,10,6,2,0)
            fx = [0] * 13
            for (xd, fy) in FTERMS:
                fx[xd] = poly_eval_mod(fy, iota, p)
            npr += sum(1 for x in range(p) if poly_eval_mod(fx, x, p) == 0)
            # g = x^2 + u*iota: 2 roots iff -u*iota is QR
            if pow((-U * iota) % p, (p - 1) // 2, p) == 1:
                npr += 2
        tot_n += npr
        tot_s += npr / p
        n += 1
    return n, tot_n / n if n else 0.0, tot_s


def poly_eval_mod(poly, x, p):
    r = 0
    for c in reversed(poly):
        r = (r * x + (c % p)) % p
    return r


# ---- part 7: individual-log descent estimate (booting + special-q tree) --------------------
# Port of the research-record model with the package's own Dickman rho (rho := 0 beyond the
# table, which only lowers success probabilities, i.e. over-estimates cost). All costs in
# instructions; every node is charged a full factor-base setup and C_SV per sieved candidate.
D_B0 = 66.0                      # descent leaves: factor-base bound (FB is 2^66.75; 66 is stricter)
D_SETUP = 2.0 ** 70              # per node: walk setup for ~2^62.2 FB ideals x 2^7.8 instr
D_ECM = {20: (11e3, 90), 25: (5e4, 300), 30: (25e4, 790), 35: (1e6, 1700), 40: (3e6, 5000),
         45: (11e6, 10600), 50: (43e6, 19300), 55: (11e7, 49000), 60: (26e7, 124000),
         65: (85e7, 210000), 70: (29e8, 340000), 75: (76e8, 700000), 80: (25e9, 1500000)}
# (B1, expected curves) per factor size in decimal digits: conservative older GMP-ECM tables


def d_ecm_cost(b_bits, size_bits):
    """Instructions to strip all prime factors <= 2^b_bits from a size_bits cofactor: 3x the
    expected curve count for the digit level, each curve ~30*B1 group ops of ~5*(limbs^2)
    word-multiplies."""
    if b_bits <= D_B0:
        return 0.0
    d = b_bits / 3.3219
    ks = sorted(D_ECM)
    k = min([x for x in ks if x >= d] or [ks[-1]])
    b1, c = D_ECM[k]
    return 3 * c * 30 * b1 * (2 * 2.5 * max(size_bits / 64, 1) ** 2)


_D_FG = {}


def _d_fg(bc, tmax, dk):
    """F(k), G(k) densities on a dk-bit grid: F counts factorisations of a k-bit rough part into
    primes in (2^B0, 2^bc], G into primes > 2^B0 (k F(k) = 1[k<=bc]/ln2 + int F)."""
    key = (bc, dk)
    nk = int(tmax / dk) + 1
    if key in _D_FG and len(_D_FG[key][0]) >= nk:
        return _D_FG[key]
    F = [0.0] * nk
    G = [0.0] * nk
    j0 = int(round(D_B0 / dk))
    for i in range(nk):
        k = i * dk
        if k < D_B0:
            continue
        sF = (1 / LN2) if k <= bc else 0.0
        sG = 1 / LN2
        # t runs over the smallest-prime size in [B0, k-B0]; j = (k - t)/dk
        tj = j0
        while tj * dk <= k - D_B0 + 1e-9:
            j = i - tj
            if tj * dk <= bc:
                sF += F[j] * dk
            sG += G[j] * dk
            tj += 1
        F[i] = sF / k
        G[i] = sG / k
    _D_FG[key] = (F, G)
    return F, G


def d_side(n, bc, T, dk=4.0):
    """(P_success, P_survive) for an n-bit norm sieved to 2^B0 with rough-part threshold 2^T."""
    F, G = _d_fg(bc, max(T, 1.0), dk)
    ps = pv = rho(n / D_B0)
    for i in range(int(T / dk) + 1):
        k = i * dk
        if k < D_B0 or k > n:
            continue
        w = LN2 * rho((n - k) / D_B0) * dk
        ps += F[i] * w
        pv += G[i] * w
    return ps, pv


def d_children(n, bc):
    """Expected number of primes in (2^B0, 2^bc] in a smooth n-bit norm (Dickman saddle)."""
    if bc <= D_B0:
        return 0.0
    u = n / bc
    xi = 1.0
    for _ in range(60):
        xi = math.log(1 + u * xi)
    t0 = D_B0 / bc
    m = 200
    return sum(math.exp(xi * (t0 + (1 - t0) * (j + .5) / m)) / (t0 + (1 - t0) * (j + .5) / m)
               for j in range(m)) * (1 - t0) / m


def d_node(nf, ng, bc, c_sv):
    best = None
    tf_grid = [d_side(nf, bc, T) + (T,) for T in range(int(D_B0), int(min(nf, 2 * bc + 400)) + 1, 20)]
    tg_grid = [d_side(ng, bc, T) + (T,) for T in range(int(D_B0), int(min(ng, 2 * bc + 400)) + 1, 20)]
    for psf, pvf, Tf in tf_grid:
        for psg, pvg, Tg in tg_grid:
            ps = psf * psg
            if ps <= 0:
                continue
            c = (c_sv + pvf * pvg * (d_ecm_cost(bc, Tf) + d_ecm_cost(bc, Tg))) / ps
            if best is None or c < best[0]:
                best = (c, ps)
    return best


def d_special_q(bq, bc, c_sv):
    """Special-q node of size bq descending to bc: the q-lattice of index 2^bq in Z^12 is sieved
    over V candidates (coefficients 2^c, c = (bq + log2 V)/12); V is the fixed point giving one
    success. Norms from the held-out pool means at A = 2^10 (alpha included); q is removed from
    the g side (worst case: leaves the f norm intact)."""
    lv = 30.0
    for _ in range(25):
        c = (bq + lv) / 12
        nf = 772.6 + 4.72 + 72 * (c - 10)
        ng = max(503.9 + 0.84 + 12 * (c - 10) - bq, 1.0)
        b = d_node(nf, ng, bc, c_sv)
        new = -math.log2(b[1])
        if abs(new - lv) < 0.05:
            break
        lv = 0.5 * lv + 0.5 * new
    return b[0] + D_SETUP, d_children(nf, bc) + d_children(ng, bc)


def d_booting(b1, c_sv, nbits=2472):
    """Booting: randomised lifts t -> z*g^t of g-side norm ~2472 bits (LLL on the weighted
    splitting lattice) until the lift is 2^b1-smooth."""
    best = None
    for T in range(int(D_B0), nbits + 1, 40):
        ps, pv = d_side(nbits, b1, T)
        if ps <= 0:
            continue
        c = (c_sv + pv * d_ecm_cost(b1, T)) / ps + D_SETUP
        if best is None or c < best:
            best = c
    return best, d_children(nbits, b1)


def descent_estimate(c_sv):
    out = []
    for chain in ([233, 166, 120, 90, 66], [200, 150, 110, 85, 66], [266, 200, 150, 110, 85, 66],
                  [166, 125, 95, 75, 66], [233, 180, 140, 105, 80, 66], [200, 140, 100, 75, 66]):
        cb, kids = d_booting(chain[0], c_sv)
        tot, nodes = cb, kids
        for bq, bc in zip(chain, chain[1:]):
            c, k = d_special_q(bq, bc, c_sv)
            tot += nodes * c
            nodes *= k
        out.append((math.log2(tot) - math.log2(UNIT), chain))
    return out


# ---- part 8: the shipped LA inner loop, re-counted from its assembly ----------------------


def check_la_kernel():
    """Count the instructions of the inner loop of _mac_run in la_mac.s (clang -O2 arm64) and
    check that it is branch-free apart from its back-edge (so static = retired per entry)."""
    lines = open(os.path.join(_data_dir(), "la_mac.s")).read().splitlines()
    on, body = False, []
    for ln in lines:
        t = ln.split(";")[0].strip()
        if ln.startswith("LBB0_2:"):
            on = True
            continue
        if on and t:
            if t.endswith(":"):
                break
            body.append(t.split()[0])
            if t.split()[0] == "b.ne":
                break
    branches = [op for op in body if op.startswith("b") and op not in ("bic", "bics", "bfi", "bfxil")]
    return len(body), branches


def main():
    global C_HIT
    cfg = json.load(sys.stdin)
    samples = int(cfg["samples"])
    seed = int(cfg.get("seed", 1))
    say(f"evidence: STNFS pairing inversion (bls12-381-fixed-g2-invert-g1); samples={samples} seed={seed}")
    say(f"units: 1 unit = {UNIT:.0f} AArch64 instructions (one dense Fp12 mul); ARCH-D sieve model "
        f"(upd={C_UPD} rep={C_REP} scan={C_SCAN} set={C_SET:.0f} hit={C_HIT:.0f}); "
        f"C_ENTRY={C_ENTRY} FLT={FLT} EXCESS={EXCESS} margin={MARGIN} aut={AUT}")

    ok0 = check_resultant()
    ok1 = check_rho()

    recs, N, count = load_inside()
    inside_map = {i: (nf, ng) for (i, nf, ng) in recs}
    nin = mism = 0
    box_nf, box_ng = [], []          # ALL recomputed indices form a fresh box sample
    hi = min(samples, max(i for i, _, _ in recs) + 1)
    for i in range(hi):
        x, nf, ng = sample_one(i)
        box_nf.append(nf)
        box_ng.append(ng)
        if quad(x) <= 1.0:
            nin += 1
            emb = inside_map.get(i)
            if emb is None or abs(emb[0] - nf) > 1.0 / 256 + 1e-9 or abs(emb[1] - ng) > 1.0 / 256 + 1e-9:
                mism += 1
        elif i in inside_map:
            mism += 1
    emb_prefix = sum(1 for i in inside_map if i < hi)
    say(f"[2] held-out pool: N={N} candidates, {count} inside the ellipsoid (embedded mass_in={count/N:.5f})")
    say(f"[2] recomputed indices [0,{hi}) exactly: {nin} inside (mass_in={nin/hi:.5f}); "
        f"embedded inside in range={emb_prefix}; mismatches={mism} {'ok' if mism == 0 else 'FAIL'}")

    qv, vbar, tco = slab_geometry()
    say(f"[3] ellipsoid volume fraction {qv:.5f} (claimed {Q_FRAC:.5f}); slab coordinates {tco}: mean "
        f"3-dim slice volume {vbar:.3f} A^3 (a full (2A)^3 cube would be 8; #slabs = E/(vbar A^3))")
    ok3 = abs(qv - Q_FRAC) < 1e-4

    NF = [r[1] for r in recs]
    NG = [r[2] for r in recs]
    mass_in = count / N
    box_rows = grid(box_nf, box_ng, 1.0, 1.0)
    box = best_cost(box_rows, 8.0)               # box region: every slab is a full cube
    ell_rows = grid(NF, NG, mass_in, Q_FRAC)
    ell = best_cost(ell_rows, vbar)
    log2V0 = 12 * math.log2(2 * 2.0 ** ell[1] + 1) - 1
    instr, lx = relcol_instr(ell[1], ell[2], log2V0, Q_FRAC, vbar)
    eff_pt = instr / (Q_FRAC * 2 ** log2V0)
    eff_credited = instr / (mass_in * 2 ** log2V0)
    say(f"[4] box baseline   log2W={box[0]:.3f}  A=2^{box[1]:.2f} B=2^{box[2]:.2f} relcol={box[3]:.3f} LA={box[4]:.3f}")
    say(f"[4] ellipsoid 1/32 log2W={ell[0]:.3f}  A=2^{ell[1]:.2f} B=2^{ell[2]:.2f} relcol={ell[3]:.3f} "
        f"LA={ell[4]:.3f} small/large threshold X=2^{lx:.2f}  (region gains {box[0]-ell[0]:.2f} bits)")
    say(f"[4] relcol per sieved point {eff_pt:.0f} instr; per credited (in-box) point {eff_credited:.0f} instr")
    say("[4] sensitivity to C_HIT (instr per large-ideal hit):")
    base_chit = C_HIT
    sweep = {}
    for ch in (0, 103, 122, 131, 150, 180, 250):
        C_HIT = float(ch)
        sweep[ch] = best_cost(ell_rows, vbar)[0]
        say(f"      C_HIT={ch:4d}  log2W={sweep[ch]:.3f}")
    C_HIT = base_chit
    flt20 = best_cost(ell_rows, vbar, flt=20.0)[0]
    say(f"[4] FLT=20 instead of 45: log2W={flt20:.3f}")
    say(f"[4] comparison: orbit-12 rho step floor 2^{RHO_FLOOR}; model total is {RHO_FLOOR - ell[0]:.1f} bits below it")

    npr, nbar_s, s_s = roots_smallprimes(1 << 11)
    say(f"[6] root stats p<=2^11: {npr} primes, mean ideals/prime={nbar_s:.3f} (embedded NBAR={NBAR}), "
        f"S(2^11)={s_s:.3f}, S(2^66) model={s_model(66.0):.2f} updates/point")

    toy_smoothness(4000, seed)

    chains = descent_estimate(eff_pt)
    d_best = min(c for c, _ in chains)
    d_two = d_best + 1.0
    p_abort = 2 ** (d_two - DESCENT_CAP)
    say(f"[7] descent per log (units), chains: " + ", ".join(f"{c[0]}->: 2^{v:.1f}" for v, c in chains))
    say(f"[7] best chain 2^{d_best:.2f} per log; two logs 2^{d_two:.2f}; meter cap 2^{DESCENT_CAP:.0f} -> "
        f"Markov P(abort) <= {p_abort:.4f}, success >= {1 - p_abort:.4f}")

    la_n, la_br = check_la_kernel()
    say(f"[8] la_mac.s inner loop: {la_n} instructions per matrix entry, branches {la_br} "
        f"({'ok' if la_n == 33 and la_br == ['b.ne'] else 'MISMATCH'})")

    ledger = [ell[0], DESCENT_CAP, LOG2_SM, LOG2_FBC, LOG2_MISC, LOG2_VERIFY]
    total = math.log2(sum(2 ** v for v in ledger))
    say(f"[ledger] relcol+LA 2^{ell[0]:.6f} + descent cap 2^{DESCENT_CAP} + Schirokauer 2^{LOG2_SM} + "
        f"FB construction 2^{LOG2_FBC} + misc 2^{LOG2_MISC} + verification 2^{LOG2_VERIFY} = 2^{total:.6f}")
    ok = ok0 and ok1 and ok3 and mism == 0 and la_n == 33 and p_abort < 0.61
    say("[done] checks:", "ok" if ok else "FAILED")
    print(json.dumps({
        "box_log2W": round(box[0], 3), "ellipsoid_log2W": round(ell[0], 6),
        "region_bits": round(box[0] - ell[0], 3),
        "ellipsoid_A_log2": round(ell[1], 3), "ellipsoid_B_log2": round(ell[2], 3),
        "small_large_threshold_log2": lx,
        "relcol_log2": round(ell[3], 3), "LA_log2": round(ell[4], 3),
        "relcol_instr_per_sieved_point": round(eff_pt, 1),
        "relcol_instr_per_credited_point": round(eff_credited, 1),
        "slab_mean_slice_volume_A3": round(vbar, 4), "ellipsoid_volume_fraction": round(qv, 6),
        "mass_in_embedded": round(mass_in, 5), "mass_in_recomputed": round(nin / hi, 5),
        "recompute_indices": hi, "recompute_mismatches": mism,
        "resultant_ok": ok0, "rho_ok": ok1,
        "c_hit_sweep": {str(k): round(v, 3) for k, v in sweep.items()}, "flt20_log2W": round(flt20, 3),
        "descent_best_log2_per_log": round(d_best, 2), "descent_cap_log2": DESCENT_CAP,
        "descent_markov_abort_bound": round(p_abort, 5),
        "la_loop_instructions": la_n,
        "total_log2W_full_ledger": round(total, 6),
        "rho_orbit12_step_floor_log2": RHO_FLOOR,
        "rootstats_primes": npr, "rootstats_nbar": round(nbar_s, 3),
        "checks_ok": ok,
        "report": _REPORT,
    }, indent=1))


if __name__ == "__main__":
    main()
