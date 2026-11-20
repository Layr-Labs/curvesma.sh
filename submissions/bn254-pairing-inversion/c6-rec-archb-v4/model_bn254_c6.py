"""
model_bn254_c6.py

Pure-Python-stdlib verification & ledger engine for the hardened stacked BN254 STNFS attack:
  - Tower: C6-m184_w216-rec (content-reduced reciprocal rank-2 CM-sextic, order-6 tau = Frob^10)
  - 100% Pure Baseline Scalar AArch64 (scalar integer ALU only, no vector extensions, no FEAT_LSE):
      C_HIT_LO=5.0, C_PROBE=18.0, C_CAND=13.0, C_HIT_HI=9.0,
      C_BKT=25.0 (12 write + 8 xfer + 5 scalar apply), C_SCAN_ZERO=2.15 (64-bit scalar SWAR),
      C_LIST_A=35.0 (3-pass pure-C build).
  - Primary model: Full de Bruijn local-density correction `rho_dB` (live norm-histogram-weighted Retention
    and Stage B) + 4-parameter deterministic hard-cap `Z[omega]` LA `(N_sw_cap, N_l_cap, k5_cap, tab_r_cap)`
    with `lazy = 0` + self-contained Band-3 7D Fincke-Pohst coset walk (ZERO imported walk tables!) +
    un-discounted Schirokauer map + live 62-row LB=59 12D descent recurrence (`T0` at `LB=59`) + live
    12D lattice norm audit + live pure-Python `filtsim_z6` + strict RT/FIX-U!
"""

import json, math, random
from math import log, log2, exp, sqrt, ceil, floor, pi, log1p, gamma
from fractions import Fraction as F
from decimal import Decimal, localcontext
from collections import Counter
from array import array
from pathlib import Path

import descent_bn254_c6 as DESC

LN2 = math.log(2.0)
GAMMA = 0.5772156649015329
UNIT = 20000.0

U_BN254 = 4965661367192848881
P_BN254 = 36 * U_BN254**4 + 36 * U_BN254**3 + 24 * U_BN254**2 + 6 * U_BN254 + 1
R_BN254 = 36 * U_BN254**4 + 36 * U_BN254**3 + 18 * U_BN254**2 + 6 * U_BN254 + 1

H_COEFFS = [3, -4, 0, 2, 1, -2, 1]
W_COEFFS = [0, -2, -1, 0, 1, -1]
SIGMA_NUM = [9, -4, -5, 1, 3, -2]
ZK_NUM = [
    [3,  0,  0,  0,  0,  0],
    [0,  3,  0,  0,  0,  0],
    [9, -4, -5,  1,  3, -2],
    [0,  2,  1, -2,  0,  1],
    [6,  1, -1, -1,  3, -1],
    [-6, 6,  6,  0, -3,  3],
]
F_REC_COEFFS = {
    8: [3, -4, -2, 0, 2, -2],
    6: [0, 6, 3, 0, -3, 3],
    4: [-12, 0, 0, 0, 0, 0],
    2: [12, -6, -3, 0, 3, -3],
    0: [-5, 4, 2, 0, -2, 2],
}
CONTENT_BITS = 12.0 + 6.0 * log2(3.0)
EXACT_PRIM_T4 = 0.683542
ALPHA_F_RAW = -20.010 + 0.08
ALPHA_G_RAW = -0.119 + 0.08

_DATA = None
def load_data():
    global _DATA
    if _DATA is None:
        _DATA = json.loads((Path(__file__).resolve().parent / "c6_data.json").read_text())
    return _DATA

# Self-contained pure baseline scalar AArch64 integer kernel prices (ZERO imported walk tables!)
C_HIT_LO = 6.0        # Band 1 [2^14, 2^43) exact 6-instr pure scalar pointwise b-lookup hit loop
C_LOOK_LO = 24.0      # Band 1 per-(b, signed ideal) lookup in 256-shell a-table
C_PROBE = 18.0        # Band 2 [2^43, 2^53.5) 64-bit bitmap word probe (pure scalar)
C_OCC_SB = 13.0       # Band 2 occupied super-bucket entry/directory/terminator overhead
C_CAND = 13.0         # Band 2 candidate check (qa <= qa_max, t < L_b, a_slot < M[t])
C_HIT_HI = 9.0        # Band 2 pure scalar score update into single-owner (t, a_slot) arena
C_FP_NODE = 20.0      # Band 3 [2^53.5, 2^58.0) 7D Fincke-Pohst internal node (shift+fsqrt+bounds)
C_FP_LEAF = 20.0      # Band 3 7D Fincke-Pohst leaf step (7 adds + t-range + rem check)
FP_WT_T = 0.40        # Optimal 7D weighted-cylinder aspect ratio (analytic min 1/sqrt(6) = 0.408)
FP_LEAVES_PER_HIT = 2.15  # Conservative bound (+12% above exact analytic Vol(B7)/Vol(Cyl7) = 1.9213!)
FP_RUN_OVERHANG = 1.05    # Run-radius Q_a variation factor across an a-run
C_RUN_INIT = 500.0    # Band 3 per-(run, ideal) 7D Babai/GS coordinate projection

def fp7_vr(wt_t=FP_WT_T):
    """Exact analytical 7D enclosing ellipsoid to cylinder volume ratio Vol(B7(sqrt(1+wt^2))) / Vol(Cyl7)."""
    v7 = (pi ** 3.5) / gamma(4.5)
    v6 = (pi ** 3.0) / gamma(4.0)
    return (v7 * ((1.0 + wt_t * wt_t) ** 3.5)) / (v6 * 2.0 * wt_t)

def fp7_gauss_nph(h, wt_t=FP_WT_T):
    """Exact parameter-free 7D Gaussian-volume-ratio internal nodes/hit at dimensionless hit density h = V_b*L_a/p."""
    v7 = (pi ** 3.5) / gamma(4.5)
    vr = fp7_vr(wt_t)
    tot = 1.0  # root node (depth 0)
    for d in range(1, 7):
        vd = (pi ** (0.5 * d)) / gamma(0.5 * d + 1.0)
        tot += (vd / (v7 ** (d / 7.0))) * ((vr * h) ** (d / 7.0))
    return tot / h
C_LLL7 = 2.0 ** 20    # Band 3 per-ideal 7D LLL basis reduction budget
C_BKT_WRITE = 12.0    # Band 3 pre-walk bucket emit + append
C_BKT_XFER = 8.0      # Band 3 explicit owner-to-owner bucket staging transfer
C_BKT_APPLY = 5.0     # Band 3 pure scalar bucket apply
C_SCAN_ZERO = 2.50    # Fused 64-bit pure-scalar SWAR varying-T(a) + b-in-P2 scan+zero (9 ops/4 active cells + setup)
C_LIST_A = 35.0       # 3-pass pure-C a-table build per entry
C_UNIT_5 = 7.8        # Unit-entry Z[omega] cost per 5-lane BW block
C_GEN_1 = 33.0        # General non-unit Z[omega] entry cost per 5-lane BW block

C_NORMS = 27000.0     # Analytical scalar norm evaluation (12-coord degree-8/2 f_rec, g_rec norms)
C_SB = 400000.0       # Analytical full-entry-norm 8-limb (<=512b) Stage-B cofactorization (P14 + MR + rho/ECM, NO 119b entry assumption!)
C_TDIV = 2.0 ** 23.45
C_BATCH = 2.0 ** 22.12
C_BATCH2 = 2.0 ** 23.13
C_GBATCH = 2.0 ** 23.0
C_LBUCKET = 2.25
LIST_NS = 256
LI_PAD = 1.001
Q_F_SAFETY = 100.0
KMAX = 4

FFT_A, FFT_B, FFT_XFER, FFT_LOCAL, FFT_COPY, FFT_ROOT, FFT_LEAF = 8.5, 1200.0, 504.0, 67108864, 3.0, 6.0, 32768.0
NTT3 = (69.75, 119.0, 90.0, 151.5, 371.0, 139.25, 271.5)
NTT3_XFER = 3.0 * FFT_XFER + 128.0
SPLIT_Q, SPLIT_REC, SPLIT_C, SPLIT_IDEAL, NORM_BITS = 1099511627776.0, 56, 33024.0, 67108864.0, 576
SPLIT_CAP, SPLIT_RES = 112.0, 84.0

# 4-parameter deterministic post-filtering hard caps (abort if ANY exceeded; LA & gen & RAM priced AT CAP with lazy=0!)
W_ROW = 201.0         # Enforced by W_target = 200 merge stop
NU_SW_CAP = 0.00310   # +91% above U'=2^22 mean (0.001625), +65% above U'=2^21 mean (0.001875)
NU_L_CAP = 0.00175    # +87% above U'=2^22 mean (0.000936), +64% above U'=2^21 mean (0.001066)
K5_CAP = 11.50        # +40% above measured k5_ev (8.19) and +13.6% above explicit-only k5 (10.12)
TAB_R_CAP = 0.85      # +14.2% above measured tab_r (0.744)
P_CAP_H3 = 0.001000   # Explicit exploratory heuristic premise H3 (NOT derived from sample variances!)

C_PROJ = 28.0
C_FULL = 97.0
C_Z64 = 40.0
C_RED9 = 270.0
BW_SIGMA = 0.0078125
PREP_OWNER = 2.0 ** 80
PREP_GUARD = 2.0 ** 80
PREP_JOINT = 2.0 ** 80
LA_RESERVE = 2.0 ** 78
G2_OUTPUT_BOUND = 69456896.0

OP = dict(
    lf=60.0, lg=59.0, exc=0.95, sk=0.65,
    p_cut=43.0, p_mid=53.5, lx2=58.0, slack=114.0, slack_b=119.0,
    w_sb=128, c_tail_z6=420.0,
    rmass=5.52, fge16=0.989, margin=0.0666, work_allow=1.0033,
)

# ==============================================================================
# 1. EXACT DICKMAN RHO & DE BRUIJN RHO_DB (PURE STDLIB)
# ==============================================================================
_RHO_COEFS = None

def _build_rho_coefs(kmax=33, nt=140):
    with localcontext() as ctx:
        ctx.prec = 70
        zero, one = Decimal(0), Decimal(1)
        coefs = [None, [one] + [zero] * (nt - 1)]
        for k in range(2, kmax + 1):
            prev, c = coefs[k - 1], [zero] * nt
            for m in range(nt - 1):
                c[m + 1] = (prev[m] + m * c[m]) / (k * (m + 1))
            c[0] = prev[0] - sum(c[1:])
            coefs.append(c)
        out = [None]
        for c in coefs[1:]:
            big = max(abs(x) for x in c)
            while len(c) > 1 and abs(c[-1]) < big * Decimal("1e-18"):
                c = c[:-1]
            out.append([float(x) for x in c])
        return out

def rho(u):
    global _RHO_COEFS
    if u <= 1.0:
        return 1.0
    if _RHO_COEFS is None:
        _RHO_COEFS = _build_rho_coefs()
    k = ceil(u)
    if k >= len(_RHO_COEFS):
        return 0.0
    x = k - u
    v = 0.0
    for c in reversed(_RHO_COEFS[k]):
        v = v * x + c
    return max(0.0, v)

def rho_db(u, log2_N):
    if u <= 1.0:
        return 1.0
    return max(1e-300, rho(u) - GAMMA * rho(u - 1.0) / max(1.0, log2_N * LN2))

def kappa_rho(u, t):
    ru = rho(u)
    return rho(u - t) / ru if ru > 0.0 else 0.0

# ==============================================================================
# 2. LIVE EULER-PRODUCT RETENTION & SURVIVOR MODEL
# ==============================================================================
def _primes_upto(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(sqrt(n)) + 1):
        if sieve[p]:
            sieve[p * p :: p] = b"\x00" * ((n - p * p) // p + 1)
    return [p for p in range(2, n + 1) if sieve[p]]

def _series_coeffs(weights, normalizer, length):
    terms = [(j, j * v) for j, v in enumerate(weights) if j and abs(v) >= 1e-24]
    out = [exp(-normalizer)] + [0.0] * length
    for n in range(1, length + 1):
        val = 0.0
        for j, w in terms:
            if j > n: break
            val += w * out[n - j]
        out[n] = val / n
    return out

class Retention:
    def __init__(self, B_bits, X_bits, prime_min=2**14, step=0.25, max_bits=600, tilt=0.05, explicit_bits=14, power_terms=8, odd=1):
        self.B, self.X, self.minimum, self.step = B_bits, X_bits, prime_min, step
        self.tilt, self.length = tilt, ceil(max_bits / step)
        self.names = ("small", "sieve_geo", "large_geo", "sieve_sf", "large_sf")
        weights = {k: [0.0] * (self.length + 1) for k in self.names}
        totals = Counter()
        cutoff = min(float(explicit_bits), B_bits)
        for p in _primes_upto(int(2 ** cutoff) - 1)[odd:]:
            b = log2(p)
            bi = round(b / step)
            q = exp(tilt * bi * step - LN2 * b)
            groups = [("small", False)] if p < prime_min else [(("sieve" if b < X_bits else "large") + "_geo", False), (("sieve" if b < X_bits else "large") + "_sf", True)]
            for name, sf in groups:
                totals[name] += log1p(q) if sf else -log1p(-q)
                qk = q
                for k in range(1, self.length + 1):
                    idx = k * bi
                    if idx > self.length: break
                    weights[name][idx] += (-1 if sf and k % 2 == 0 else 1) * qk / k
                    qk *= q
                    if qk < 1e-25: break
        def continuum(left, right, pref):
            if left >= right: return
            w = step / 4.0
            cells = ceil((right - left) / w)
            w = (right - left) / cells
            for i in range(cells):
                b = left + (i + 0.5) * w
                bi = round(b / step)
                for k in range(1, power_terms + 1):
                    amt = exp(k * tilt * bi * step - (k - 1) * LN2 * b) * w / (k * b)
                    idx = k * bi
                    sign = -1 if k % 2 == 0 else 1
                    grps = ((pref, amt),) if pref == "small" else ((pref + "_geo", amt), (pref + "_sf", sign * amt))
                    for g, x in grps:
                        totals[g] += x
                        if idx <= self.length: weights[g][idx] += x
        lm = log2(prime_min)
        continuum(cutoff, min(lm, B_bits), "small")
        continuum(max(cutoff, lm), min(X_bits, B_bits), "sieve")
        continuum(max(cutoff, X_bits), B_bits, "large")
        self.weights, self.normalizers = weights, {k: totals[k] for k in self.names}
        all_w = [sum(weights[k][j] for k in ("small", "sieve_geo", "large_geo")) for j in range(self.length + 1)]
        self.denominator = _series_coeffs(all_w, sum(totals[k] for k in ("small", "sieve_geo", "large_geo")), self.length)
        self.selection_normalizer = exp(totals["sieve_sf"] + totals["large_sf"] - totals["sieve_geo"] - totals["large_geo"])
        self._c = {}

    def co(self, k):
        if k not in self._c:
            self._c[k] = _series_coeffs(self.weights[k], self.normalizers[k], self.length)
        return self._c[k]

    def two(self, n, a, b=9**9, win=None, sf=False):
        a = min(a, n)
        if a < 0 or b < 0: return 0.0
        key = (n, b, win, sf)
        if key not in self._c:
            if key[1:] not in self._c:
                st = self.step; sm = self.co("small"); D = min(self.length, 2048); m = [0.0] * (D + 1)
                for l, v in enumerate(self.co("large_sf" if sf else "large_geo")[: min(b, D) + 1]):
                    if v and (win is None or l * st <= win[2] or any(self.X * k - win[2] <= l * st <= win[1] * k + win[2] for k in range(1, win[0] + 1))):
                        for t in range(D + 1 - l): m[l + t] += v * sm[t]
                self._c[key[1:]] = m
            m = self._c[key[1:]]; sv = self.co("sieve_sf" if sf else "sieve_geo"); c = array("d"); t = 0.0
            for x in range(min(n, len(m) - 1) + 1):
                t += m[x] * sv[n - x]; c.append(t)
            self._c[key] = c
        c = self._c[key]
        return c[min(a, len(c) - 1)]

    def cap(self, nb): return floor(ceil(nb) / log2(self.minimum))
    def at(self, nb, sl, sb=None):
        n = round(nb / self.step); rb = (self.cap(nb) + 1) / 16.0
        r = self.selection_normalizer * self.two(n, floor((sl - rb) / self.step), 9**9 if sb is None else floor((sb - rb - 0.125) / self.step), None, True) / self.denominator[n]
        return max(0.0, min(1.0, r))
    def survivor_upper(self, nb, sl, ex=0.0, sb=None, lb=None, kmax=4):
        n = round(nb / self.step); a = floor((sl + ex) / self.step)
        return (self.two(n, a) if sb is None else self.two(n, a, floor((sb + 0.25) / self.step), (kmax, lb, (self.cap(nb) + 2) / 16.0 + 0.25))) / self.denominator[n]

# ==============================================================================
# 3. LIVE ALGEBRA + LIVE MOD-R CRT FACTORIZATION OF h, f_rec, AND g_rec ([0])
# ==============================================================================
def _ptrim(a):
    while len(a) > 1 and a[-1] == 0: a.pop()
    return a
def _pmul(a, b, mod):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y: c[i + j] = (c[i + j] + x * y) % mod
    return _ptrim(c)
def _pmod(a, b, mod):
    a = [x % mod for x in a]; b = _ptrim([x % mod for x in b])
    inv = pow(b[-1], -1, mod)
    while len(a) >= len(b) and not (len(a) == 1 and a[0] == 0):
        d = len(a) - len(b); fac = (a[-1] * inv) % mod
        for i in range(len(b)): a[d + i] = (a[d + i] - fac * b[i]) % mod
        _ptrim(a)
    return a
def _pdiv(a, b, mod):
    a = _ptrim([x % mod for x in a]); b = _ptrim([x % mod for x in b])
    q = [0] * max(1, len(a) - len(b) + 1); inv = pow(b[-1], -1, mod)
    while len(a) >= len(b) and not (len(a) == 1 and a[0] == 0):
        d = len(a) - len(b); fac = (a[-1] * inv) % mod; q[d] = fac
        for i in range(len(b)): a[d + i] = (a[d + i] - fac * b[i]) % mod
        _ptrim(a)
    return _ptrim(q), a
def _pgcd(a, b, mod):
    a = _ptrim([x % mod for x in a]); b = _ptrim([x % mod for x in b])
    while not (len(b) == 1 and b[0] == 0):
        a, b = b, _pmod(a, b, mod)
    inv = pow(a[-1], -1, mod)
    return [(x * inv) % mod for x in a]
def _ppowmod(base, e, m, mod):
    res = [1]; b = _pmod(base, m, mod)
    while e > 0:
        if e & 1: res = _pmod(_pmul(res, b, mod), m, mod)
        b = _pmod(_pmul(b, b, mod), m, mod); e >>= 1
    return res
def _ddf_degrees(f, mod):
    f = _ptrim([x % mod for x in f]); inv = pow(f[-1], -1, mod); f = [(x * inv) % mod for x in f]
    df = [(i * f[i]) % mod for i in range(1, len(f))]
    if len(_pgcd(f, df, mod)) != 1: return None
    degs = []; h = [0, 1]; rem = f[:]; d = 1
    while len(rem) - 1 >= 2 * d:
        h = _ppowmod(h, mod, rem, mod)
        hx = h[:] + [0] * max(0, 2 - len(h)); hx[1] = (hx[1] - 1) % mod
        g = _pgcd(hx, rem, mod); gd = len(g) - 1
        if gd > 0:
            degs.extend([d] * (gd // d))
            rem, _ = _pdiv(rem, g, mod)
            h = _pmod(h, rem, mod)
        d += 1
    if len(rem) - 1 > 0: degs.append(len(rem) - 1)
    return sorted(degs)
def _roots_mod(f, mod):
    rng = random.Random(12345)
    f = _ptrim([x % mod for x in f]); inv = pow(f[-1], -1, mod); f = [(x * inv) % mod for x in f]
    out, stack = [], [f]
    while stack:
        cur = stack.pop(); deg = len(cur) - 1
        if deg == 0: continue
        if deg == 1: out.append((-cur[0]) % mod); continue
        while True:
            a = rng.randrange(mod)
            pw = _ppowmod([a, 1], (mod - 1) // 2, cur, mod)
            pw[0] = (pw[0] - 1) % mod
            g = _pgcd(pw, cur, mod)
            if 0 < len(g) - 1 < deg:
                q, _ = _pdiv(cur, g, mod)
                stack.append(g); stack.append(q); break
    return sorted(out)

def _pmul_h(a, b, mod):
    r = _pmod(_pmul(a, b, mod), H_COEFFS, mod)
    return (r + [0] * 6)[:6]
def _ppow_h(a, e, mod):
    r = _ppowmod(a, e, H_COEFFS, mod)
    return (r + [0] * 6)[:6]
def _pinv_h(a, mod):
    r0, r1 = [x % mod for x in H_COEFFS], _ptrim([x % mod for x in a])
    t0, t1 = [0], [1]
    while not (len(r1) == 1 and r1[0] == 0):
        q, rem = _pdiv(r0, r1, mod)
        qt = _pmul(q, t1, mod)
        L = max(len(t0), len(qt))
        t2 = _ptrim([((t0[i] if i < len(t0) else 0) - (qt[i] if i < len(qt) else 0)) % mod for i in range(L)])
        r0, r1, t0, t1 = r1, rem, t1, t2
    inv = pow(r0[0], -1, mod)
    return ([(x * inv) % mod for x in t0] + [0] * 6)[:6]
def _norm_Kh_exact(coeffs_int):
    M = [[F(0)] * 6 for _ in range(6)]; cur = [F(c) for c in coeffs_int]
    for col in range(6):
        for r in range(6): M[r][col] = cur[r]
        top = cur[5]; nxt = [F(0)] + cur[:5]
        for j in range(6): nxt[j] -= top * H_COEFFS[j]
        cur = nxt
    det = F(1)
    for col in range(6):
        piv = next(r for r in range(col, 6) if M[r][col] != 0)
        if piv != col: M[col], M[piv] = M[piv], M[col]; det = -det
        pv = M[col][col]; det *= pv
        for r in range(col + 1, 6):
            if M[r][col] != 0:
                fac = M[r][col] / pv
                for j in range(col, 6): M[r][j] -= fac * M[col][j]
    return int(det)

def verify_tower_algebra_and_crt():
    """Live big-int tower algebra AND live mod-r CRT factorization of h, f_rec, AND g_rec!"""
    ord_ok = (pow(P_BN254, 12, R_BN254) == 1) and (pow(P_BN254, 6, R_BN254) != 1) and (pow(P_BN254, 4, R_BN254) != 1)
    y = [0, 1, 0, 0, 0, 0]
    inv3 = pow(3, -1, P_BN254)
    sig_p = [(c * inv3) % P_BN254 for c in SIGMA_NUM]
    y_p1 = _ppow_h(y, P_BN254, P_BN254)
    y_p2 = _ppow_h(y_p1, P_BN254, P_BN254)
    y_p3 = _ppow_h(y_p2, P_BN254, P_BN254)
    y_p4 = _ppow_h(y_p3, P_BN254, P_BN254)
    y_p6 = _ppow_h(_ppow_h(y_p4, P_BN254, P_BN254), P_BN254, P_BN254)
    h_irred = (y_p6 == y) and (y_p2 != y) and (y_p3 != y)
    sig_frob_ok = (y_p4 == sig_p)

    uw = [(U_BN254 * c) % P_BN254 for c in W_COEFFS]
    uw_half = _ppow_h(uw, (P_BN254 - 1) // 2, P_BN254)
    cur, acc = uw_half, [1, 0, 0, 0, 0, 0]
    val_p6 = val_p10 = None
    for j in range(10):
        acc = _pmul_h(acc, cur, P_BN254)
        if j == 5: val_p6 = acc[:]
        if j == 9: val_p10 = acc[:]
        cur = _ppow_h(cur, P_BN254, P_BN254)
    minus_one = [P_BN254 - 1, 0, 0, 0, 0, 0]
    g_irred = (val_p6 == minus_one)
    tau_frob10 = sig_frob_ok and (val_p10 == minus_one)

    chi = pow(P_BN254, 10, R_BN254)
    chi_expr = 36 * U_BN254**3 + 18 * U_BN254**2 + 6 * U_BN254 + 2
    chi_ok = (chi == chi_expr) and ((chi * chi - chi + 1) % R_BN254 == 0) and ((pow(chi, 3, R_BN254) + 1) % R_BN254 == 0)

    inv_uw = _pinv_h(uw, P_BN254)
    f_eval = [0] * 6; x2_pow = [1, 0, 0, 0, 0, 0]
    for d in (0, 2, 4, 6, 8):
        term = _pmul_h([c % P_BN254 for c in F_REC_COEFFS[d]], x2_pow, P_BN254)
        f_eval = [(a + b) % P_BN254 for a, b in zip(f_eval, term)]
        x2_pow = _pmul_h(x2_pow, inv_uw, P_BN254)
    f_vanishes = all(c == 0 for c in f_eval)
    nw, n_lead, n_const = abs(_norm_Kh_exact(W_COEFFS)), abs(_norm_Kh_exact(F_REC_COEFFS[8])), abs(_norm_Kh_exact(F_REC_COEFFS[0]))

    # LIVE MOD-r CRT FACTORIZATION OF h, f_rec, AND g_rec ACROSS ALL 6 ROOTS OF h mod r!
    h_mod_r_degs = _ddf_degrees(H_COEFFS, R_BN254)
    rts_r = _roots_mod(H_COEFFS, R_BN254)
    g_mod_r_degs = []
    f_mod_r_degs = []
    for s in rts_r:
        ws = sum(W_COEFFS[j] * pow(s, j, R_BN254) for j in range(6)) % R_BN254
        g_mod_r_degs.append(_ddf_degrees([(-1) % R_BN254, 0, (U_BN254 * ws) % R_BN254], R_BN254))
        fp = [0] * 9
        for d, cs in F_REC_COEFFS.items():
            fp[d] = sum(cs[j] * pow(s, j, R_BN254) for j in range(6)) % R_BN254
        f_mod_r_degs.append(_ddf_degrees(fp, R_BN254))
    crt_ok = (h_mod_r_degs == [1, 1, 1, 1, 1, 1] and len(rts_r) == 6
              and all(gd == [2] for gd in g_mod_r_degs)
              and all(fd == [2, 2, 2, 2] for fd in f_mod_r_degs))

    T4 = load_data()["t4_geom"]
    S = T4["SMAT"]
    def iso_err(G):
        return max(abs(sum(S[a][i] * G[a][b] * S[b][j] for a in range(6) for b in range(6)) - G[i][j]) for i in range(6) for j in range(6))
    iso_a, iso_b = iso_err(T4["GA"]), iso_err(T4["GB"])
    ok = ord_ok and h_irred and g_irred and tau_frob10 and chi_ok and f_vanishes and crt_ok and (nw == 216) and (n_lead == 729) and (n_const == 729) and (iso_a < 1e-10) and (iso_b < 1e-10)
    return dict(
        ok=ok, ord_r_p_12=ord_ok, h_irred_mod_p=h_irred, g_rec_irred_mod_p6=g_irred,
        tau_is_frob10=tau_frob10, chi_192bit_ok=chi_ok, f_rec_vanishes_at_inv_uw=f_vanishes,
        norm_w=nw, norm_f_lead=n_lead, norm_f_const=n_const,
        mod_r_h_degs=h_mod_r_degs, mod_r_g_degs=g_mod_r_degs, mod_r_f_degs=f_mod_r_degs,
        mod_r_crt_all_deg2_ok=crt_ok, sigma_iso_err_Ga=iso_a, sigma_iso_err_Gb=iso_b,
    )

# ==============================================================================
# 4. LIVE EXHAUSTIVE SMALL-RADIUS FUNDAMENTAL DOMAIN AUDIT (R^2 = 50) ([1])
# ==============================================================================
def verify_exhaustive_fd(R2=50.0):
    """Exhaustive Fincke-Pohst enumeration of ALL nonzero v with T2(v) <= R2: asserts 0 missed & 0 dups!"""
    T4 = load_data()["t4_geom"]
    EMB_ZK = [[complex(r, im) for r, im in row] for row in T4["EMB_ZK"]]
    FAC6, UMATS, S = T4["FAC6"], T4["UMATS"], T4["SMAT"]
    GT2 = [[sum(2.0 * (EMB_ZK[j][a] * EMB_ZK[j][b].conjugate()).real for j in range(3)) for b in range(6)] for a in range(6)]

    def mv(M, v): return tuple(sum(M[r][c] * v[c] for c in range(6)) for r in range(6))
    def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(6)) for j in range(6)] for i in range(6)]
    S2 = mm(S, S)
    def qf(v):
        s = 0.0
        for r in range(6):
            vr = v[r]
            if vr:
                s += GT2[r][r] * vr * vr
                for c in range(r + 1, 6):
                    if v[c]: s += 2.0 * GT2[r][c] * vr * v[c]
        return s

    L = [[0.0] * 6 for _ in range(6)]
    for r in range(6):
        for c in range(r + 1):
            s = GT2[r][c] - sum(L[r][k] * L[c][k] for k in range(c))
            L[r][c] = sqrt(s) if r == c else s / L[c][c]
    pts, v = [], [0] * 6
    def enum(idx, rem):
        shift = sum(L[j][idx] * v[j] for j in range(idx + 1, 6))
        diag = L[idx][idx]; rad = sqrt(max(0.0, rem))
        lo = ceil((-rad - shift) / diag - 1e-12); hi = floor((rad - shift) / diag + 1e-12)
        for x in range(lo, hi + 1):
            d = diag * x + shift; nrem = rem - d * d
            if nrem >= -1e-9:
                v[idx] = x
                if idx == 0:
                    if any(v): pts.append(tuple(v))
                else:
                    enum(idx - 1, max(0.0, nrem))
    enum(5, R2 + 1e-9)

    def xv_of(p):
        return [2.0 * ((z := sum(EMB_ZK[j][k] * p[k] for k in range(6))).real**2 + z.imag**2) for j in range(3)]
    def is_can(p):
        if tuple(-x for x in p) < p: return False
        xv = xv_of(p)
        for fi in range(6):
            dot = sum(FAC6[fi][j] * xv[j] for j in range(3))
            if dot < -1e-9: return False
            if abs(dot) <= 1e-9:
                uv = mv(UMATS[fi], p); nuv = tuple(-x for x in uv)
                if (nuv if nuv < uv else uv) < p: return False
        for Sm in (S, S2):
            vk = mv(Sm, p); nvk = tuple(-x for x in vk)
            if (nvk if nvk < vk else vk) < p: return False
        return True
    def orbit_rep(p):
        cur = p
        for _ in range(50):
            xv = xv_of(cur); bf = -1; bd = -1e-9
            for fi in range(6):
                dot = sum(FAC6[fi][j] * xv[j] for j in range(3))
                if dot < bd: bd = dot; bf = fi
            if bf < 0: break
            cur = mv(UMATS[bf], cur)
        mq = qf(cur); comp = {cur}; q = [cur]
        while q:
            u = q.pop()
            for w in [tuple(-x for x in u), mv(S, u), mv(S2, u)] + [mv(U, u) for U in UMATS]:
                if w not in comp and abs(qf(w) - mq) <= 1e-7:
                    comp.add(w); q.append(w)
        return min(comp)

    orbits = {}
    for p in pts:
        orbits.setdefault(orbit_rep(p), []).append(p)
    exact1 = missed = dups = 0
    for members in orbits.values():
        c = sum(1 for m in members if is_can(m))
        if c == 1: exact1 += 1
        elif c == 0: missed += 1
        else: dups += c - 1
    return dict(
        ok=(exact1 == len(orbits) and missed == 0 and dups == 0 and len(orbits) == 107),
        R2=R2, nonzero_pts=len(pts), total_orbits=len(orbits),
        exact1=exact1, missed_orbits=missed, duplicate_orbits=dups,
    )

# ==============================================================================
# 5. LIVE ARCH B-HYBRID PROTOTYPE + 64-BIT SCALAR SWAR SCAN CHECK ([2])
# ==============================================================================
def verify_archB_and_scalar_swar():
    T4 = load_data()["t4_geom"]
    GA, GB, S, FAC6, UMATS = T4["GA"], T4["GB"], T4["SMAT"], T4["FAC6"], T4["UMATS"]
    EMB_ZK = [[complex(r, im) for r, im in row] for row in T4["EMB_ZK"]]
    EMB_W = [complex(r, im) for r, im in T4["EMB_W"]]

    # 1. Verify 64-bit pure-scalar SWAR varying-T(a) + b-in-P2 scan+zero identity across 10,000 random 4-tuples
    rng = random.Random(20261120)
    swar_errs = 0
    MASK = 0x8000800080008000
    P2_BIAS = 0x0030003000300030  # +48 = +3*16 quantized bits when b(t) in P_2
    for _ in range(10000):
        t_a = [rng.randint(500, 12000) for _ in range(4)]
        b_in_p2 = (rng.random() < (1.0 / 9.0))
        sc = [rng.randint(0, 12500) for _ in range(4)]
        v64 = sc[0] | (sc[1] << 16) | (sc[2] << 32) | (sc[3] << 48)
        b4 = [(0x8000 - t) for t in t_a]
        bias64 = b4[0] | (b4[1] << 16) | (b4[2] << 32) | (b4[3] << 48)
        if b_in_p2:
            v64 = v64 + P2_BIAS
        swar_hit = (((v64 + bias64) & MASK) != 0)
        shift = 48 if b_in_p2 else 0
        exact_hit = any((sc[i] + shift) >= t_a[i] for i in range(4))
        if swar_hit != exact_hit:
            swar_errs += 1

    # 2. Check N_g,rec(a,b) == N_g,rec(a,0) invariance on 500 production-scale pairs
    max_inv = 0.0
    for _ in range(500):
        a = [rng.randint(-450, 450) for _ in range(6)]
        b = [rng.randint(-750, 750) for _ in range(6)]
        if not any(a) or not any(b): continue
        lab = la0 = 0.0
        for j in range(3):
            za = sum(EMB_ZK[j][k] * a[k] for k in range(6))
            zb = sum(EMB_ZK[j][k] * b[k] for k in range(6))
            uw_a2 = U_BN254 * EMB_W[j] * (za * za)
            vab = uw_a2 - zb * zb
            lab += log2(vab.real**2 + vab.imag**2)
            la0 += log2(uw_a2.real**2 + uw_a2.imag**2)
        max_inv = max(max_inv, abs(lab - la0))

    # 3. Reduced-scale Arch B-Hybrid exactness prototype on 24 real g_rec ideal pairs
    def mv(M, v): return [sum(M[r][c] * v[c] for c in range(6)) for r in range(6)]
    def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(6)) for j in range(6)] for i in range(6)]
    S2 = mm(S, S)
    def qf(G, v):
        s = 0.0
        for i in range(6):
            vi = v[i]
            if vi:
                s += G[i][i] * vi * vi
                for j in range(i + 1, 6):
                    if v[j]: s += 2.0 * G[i][j] * vi * v[j]
        return s
    def is_can(a):
        if tuple(-x for x in a) < tuple(a): return False
        xv = [2.0 * ((z := sum(EMB_ZK[j][k] * a[k] for k in range(6))).real**2 + z.imag**2) for j in range(3)]
        for fi in range(6):
            dot = sum(FAC6[fi][j] * xv[j] for j in range(3))
            if dot < -1e-9: return False
            if abs(dot) <= 1e-9:
                uv = mv(UMATS[fi], a); nuv = [-x for x in uv]
                if tuple(nuv if tuple(nuv) < tuple(uv) else uv) < tuple(a): return False
        for Sm in (S, S2):
            vk = mv(Sm, a); nvk = [-x for x in vk]
            if tuple(nvk if tuple(nvk) < tuple(vk) else vk) < tuple(a): return False
        return True
    def eval_poly_mod(cs, s, p):
        v = 0
        for c in reversed(cs): v = (v * s + c) % p
        return v
    def eval_zk_mod(vec, s, p):
        inv3 = pow(3, -1, p); sp = [1]
        for _ in range(5): sp.append((sp[-1] * s) % p)
        return (sum(vec[k] * (sum(ZK_NUM[k][j] * sp[j] for j in range(6)) % p) for k in range(6) if vec[k]) * inv3) % p
    def is_prime(n):
        return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))

    ideals = []
    for p in range(101, 1500, 2):
        if not is_prime(p): continue
        for s in range(p):
            if eval_poly_mod(H_COEFFS, s, p) == 0:
                uw = ((U_BN254 % p) * eval_poly_mod(W_COEFFS, s, p)) % p
                if uw == 0: continue
                inv_uw = pow(uw, -1, p)
                if pow(inv_uw, (p - 1) // 2, p) == 1:
                    for r in range(1, (p + 1) // 2):
                        if (r * r) % p == inv_uw: ideals.append((p, s, r)); break
        if len(ideals) >= 24: break

    RA2, RB2, KAP = 36.0, 110.0, 0.322
    a_pts = []
    for a5 in range(-2, 3):
        for a4 in range(-2, 3):
            for a3 in range(-2, 3):
                for a2 in range(-2, 3):
                    for a1 in range(-2, 3):
                        for a0 in range(-2, 3):
                            av = [a0, a1, a2, a3, a4, a5]
                            if any(av) and (qv := qf(GA, av)) <= RA2 and is_can(av):
                                a_pts.append((qv / RA2, av))
    a_pts.sort(key=lambda x: x[0]); a_pts = a_pts[:80]
    b_rows, b_pts = [], []
    for b5 in range(1, 3):
        for b4 in range(-2, 3):
            for b3 in range(-2, 3):
                for b2 in range(-2, 3):
                    for b1 in range(-2, 3):
                        rp = [(b0, bv, qv / RB2) for b0 in range(-6, 7) if (qv := qf(GB, (bv := [b0, b1, b2, b3, b4, b5]))) <= RB2]
                        if rp:
                            b_rows.append((len(rp), [x[1] for x in rp], [x[2] for x in rp], min(x[2] for x in rp)))
                            for _, bv, qn in rp: b_pts.append((bv, qn, len(b_pts)))
                        if len(b_rows) >= 50: break
                    if len(b_rows) >= 50: break
                if len(b_rows) >= 50: break
            if len(b_rows) >= 50: break
        if len(b_rows) >= 50: break

    gt = {(ai, bi): 0 for ai, (qan, _) in enumerate(a_pts) for bi, (_, qbn, _) in enumerate(b_pts) if qbn <= 1.0 - KAP * qan + 1e-12}
    tot_hits = 0
    for p, s, rho_r in ideals:
        ar = [eval_zk_mod(av, s, p) for _, av in a_pts]
        br = [eval_zk_mod(bv, s, p) for bv, _, _ in b_pts]
        for ai, bi in gt:
            rb = (rho_r * br[bi]) % p
            c = int(ar[ai] == rb) + int(ar[ai] == (-rb) % p)
            if c: gt[(ai, bi)] += c; tot_hits += c

    W_SB = 4; sc = {k: 0 for k in gt}; b_map = {tuple(bv): bi for bv, _, bi in b_pts}
    for p, s, rho_r in ideals:
        inv_r = pow(rho_r, -1, p); n_sb = (p + W_SB - 1) // W_SB
        tab = {}; bm = [0] * ((n_sb + 63) // 64 + 1)
        for ai, (qan, av) in enumerate(a_pts):
            kta = (eval_zk_mod(av, s, p) * inv_r) % p; q = kta // W_SB
            tab.setdefault(q, []).append((qan, kta, ai)); bm[q >> 6] |= (1 << (q & 63))
        for L_b, b_vecs, qb_norms, qb_min in b_rows:
            qa_max = (1.0 - qb_min) / KAP + 1e-12; k0 = eval_zk_mod(b_vecs[0], s, p)
            for sign in (1, -1):
                lo = k0 if sign == 1 else (-(k0 + L_b - 1)) % p
                hi = (k0 + L_b - 1) % p if sign == 1 else (-k0) % p
                rngs = [(lo // W_SB, hi // W_SB)] if lo <= hi else [(lo // W_SB, n_sb - 1), (0, hi // W_SB)]
                for q_lo, q_hi in rngs:
                    for q in range(q_lo, q_hi + 1):
                        if not ((bm[q >> 6] >> (q & 63)) & 1): continue
                        for qan, kta, ai in tab[q]:
                            if qan > qa_max: break
                            t = (kta - k0) % p if sign == 1 else (-kta - k0) % p
                            if t < L_b and qb_norms[t] <= 1.0 - KAP * qan + 1e-12:
                                sc[(ai, b_map[tuple(b_vecs[t])])] += 1
    mismatches = sum(1 for k in gt if sc[k] != gt[k])
    return dict(
        ok=(mismatches == 0 and swar_errs == 0 and max_inv <= 1e-12),
        scalar_swar_10k_errors=swar_errs,
        max_ng_invariance_bits=max_inv,
        n_ideals=len(ideals), n_a=len(a_pts), n_b=len(b_pts),
        active_pairs=len(gt), total_hits=tot_hits, mismatches=mismatches,
    )

# ==============================================================================
# 6. LIVE FULL 16-ROW LEDGER (100% SELF-CONTAINED, PURE SCALAR, FULL DE BRUIJN)
# ==============================================================================
def log2_li(log2x):
    x = log2x * LN2
    return log2x - log2(x) + log2(1.0 + 1.0 / x + 2.0 / (x * x))

def log2_U_z6(lf, lg):
    uf = (2.0 ** log2_li(lf)) / 6.0
    ug = (2.0 ** log2_li(lg)) / 6.0
    deg2 = (2.0 ** (0.5 * lf) / (0.5 * lf * LN2) + 2.0 ** (0.5 * lg) / (0.5 * lg * LN2)) / 6.0
    return log2(uf + ug + deg2)

def fft_F(N):
    return FFT_A * N * (log2(max(N, 2.0)) ** 2) + FFT_B + (FFT_XFER * N if N > FFT_LOCAL else 0.0)

def ntt_F(N):
    k = max(2, ceil(log2(N))); n = 2.0 ** k; o = k % 2
    return NTT3[0]*n*k + NTT3[5 if o else 1]*n + NTT3[2]*N + NTT3[3]*k + NTT3[6 if o else 4] + (NTT3_XFER*n if n > FFT_LOCAL else 0.0)

def mul_F(N):
    return fft_F(N) if max(2, ceil(log2(N))) > 55 else min(fft_F(N), ntt_F(N))

def pbits(lo, hi):
    return 1.000028 * (2.0 ** hi - (2.0 ** lo if lo else 0)) / LN2 if lo == 0 else 1.000028 * (2.0 ** hi) / LN2

def gbatch_cost(pb, leaf=9, L=None, c_leaf=FFT_LEAF):
    SP = ceil(pb / 64.0)
    L = L if L is not None else ceil(log2((SP + 1) / leaf))
    G = L + 2; n = lambda j: leaf * (2 ** j)
    PT = sum((2 ** (L - j)) * (mul_F(2 * n(j - 1)) + FFT_COPY * 2 * n(j - 1)) for j in range(1, L + 1))
    RT = sum((2 ** (L - j + 1)) * (mul_F(n(j) + G + n(j - 1)) + FFT_COPY * (n(j - 1) + G)) for j in range(1, L + 1))
    return (PT + RT + FFT_ROOT * mul_F(2 * (max(SP, n(L)) + G))) / (2 ** L) + c_leaf

def split_parts(R, lb):
    np_ = 1.25506 * (2.0 ** lb) / (lb * LN2); h = ceil(log2(np_)); D = h + 1
    SP = 2.0 ** ceil(log2(np_ * (lb / 64.0 + 1.0)))
    L = ceil(log2(SPLIT_REC * SPLIT_Q)); G = L + 2; C = SPLIT_C
    W = R / SPLIT_Q + 1.0; B = SPLIT_REC * R * D; occ = 2 * NORM_BITS
    return dict(
        trees=B * (2 * (C * (36 * L + 2 * G) + 3 * FFT_B) + 2**15),
        roots=6 * C * (2 * SP * W * D + 136 * B) + 6 * FFT_B * B,
        gcd=B * 2.0**24, small=R * 2.0**26, output=occ * R * 2.0**20,
        ideal=occ * R * SPLIT_IDEAL, control=B * 2.0**16,
        ptree=D * SP * C + 2 * FFT_B * (2.0**h), debit=128 * (W + 2),
    )

def _det_mat(M):
    n = len(M); A = [row[:] for row in M]; det = 1.0
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(A[r][col]))
        if piv != col: A[col], A[piv] = A[piv], A[col]; det = -det
        pv = A[col][col]; det *= pv
        for r in range(col + 1, n):
            fac = A[r][col] / pv
            for j in range(col, n): A[r][j] -= fac * A[col][j]
    return det

_RET_CACHE = {}
_PEEL_CACHE = {}
_DESC_CACHE = {}

def _peel_fixed_point(lf, lg, exc, nf=405.0, ng=469.0):
    lU = log2_U_z6(lf, lg)
    terms = []
    for m, lb in ((nf, lf), (ng, lg)):
        u = m / lb; steps = int(round((lb - 30.0) * 10))
        for i in range(steps):
            b = 30.0 + i * 0.1; c = b + 0.1
            k = kappa_rho(u, c / lb)
            terms.append((2.0 ** (lU - c) * 6.0 * k, k * log(c / b)))
    def run(e_val):
        def step(h):
            tot = 0.0; omh = 1.0 - h
            for lam0, wt in terms:
                lam = lam0 * e_val
                if lam < 100.0:
                    tot += wt * exp(-lam * omh)
            return 1.0 - exp(-tot)
        h = 1.0
        for _ in range(100):
            hn = step(h)
            if abs(hn - h) < 1e-10: return hn
            h = hn
        return h
    h_val = run(exc)
    lo, hi = 0.50, 1.00
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if run(mid) > 0.75: lo = mid
        else: hi = mid
    return dict(h=h_val, exc_crit=hi, stalls=(h_val > 0.75))

def compute_ledger(track="bn254-gt-dlp", use_db=True, c_tail_z6=None, cap_scale=1.0):
    D = load_data()
    PQ = D["pool_quad"]
    T4 = D["t4_geom"]
    FILT = D["filtsim_raw"]
    lf, lg, exc = OP["lf"], OP["lg"], OP["exc"]
    p_cut, pm, lx2, sl, sb, w_sb = OP["p_cut"], OP["p_mid"], OP["lx2"], OP["slack"], OP["slack_b"], OP["w_sb"]
    if c_tail_z6 is None: c_tail_z6 = OP["c_tail_z6"]

    # 1. Live Euler-product Retention & Survivor model (histogram-weighted across ALL 200k-pool norm bins!)
    rkey = (lg, lx2)
    if rkey not in _RET_CACHE:
        _RET_CACHE[rkey] = (
            Retention(lg, lx2, prime_min=2**14, explicit_bits=14, step=0.25, max_bits=600, tilt=0.05, odd=1),
            Retention(600.0, lx2, prime_min=2**14, explicit_bits=14, step=0.25, max_bits=601, tilt=0.005, odd=1),
        )
    ret_obj, cost_obj = _RET_CACHE[rkey]

    # 2. Live 200k-pool self-consistent yield (PRIMARY = FULL DE BRUIJN `rho_db` + bin-weighted Retention!)
    prim_tail_corr = EXACT_PRIM_T4 / PQ["p_prim"]
    lU = log2_U_z6(lf, lg)
    target_lr = lU + log2(exc) + OP["margin"]
    log2_V_sub0 = PQ["log2_V0"] + log2(PQ["frac_reg"])
    n_in = PQ["n_in"]
    bins = PQ["bins"]

    def eval_d(d):
        ys = ys2 = nf_w = ng_w = rg_w = 0.0
        for lnf0, lng0, cnt in bins:
            nf_raw = lnf0 + 48.0 * d; ng_raw = lng0 + 12.0 * d
            nf = nf_raw + ALPHA_F_RAW; ng = ng_raw + ALPHA_G_RAW
            ng_q = round(ng_raw * 2.0) / 2.0
            r_g_bin = (8.0 / 9.0) * ret_obj.at(ng_q, sl, sb) + (1.0 / 9.0) * ret_obj.at(ng_q - 3.0, sl - 3.0, sb)
            y = (rho_db(nf / lf, nf) * rho_db(ng / lg, ng) if use_db else rho(nf / lf) * rho(ng / lg)) * r_g_bin * (1.0 - 0.0014) * prim_tail_corr
            wy = cnt * y; ys += wy; ys2 += cnt * (y * y)
            nf_w += wy * nf_raw; ng_w += wy * ng_raw; rg_w += wy * r_g_bin
        my = ys / n_in
        return log2_V_sub0 + 12.0 * d + log2(max(my, 1e-300)), my, ys2 / n_in, nf_w / ys - CONTENT_BITS, ng_w / ys, rg_w / ys

    lo, hi = -2.5, +2.5
    for _ in range(32):
        mid = 0.5 * (lo + hi)
        if eval_d(mid)[0] < target_lr: lo = mid
        else: hi = mid
    d_opt = 0.5 * (lo + hi)
    lr, my, my2, nf_rel, ng_rel, r_g = eval_d(d_opt)
    se_bits = sqrt(max(my2 - my * my, 0.0) / n_in) / max(my, 1e-300) / LN2

    # 3. Exact C6 geometry + 100% SELF-CONTAINED PURE SCALAR Arch B-Hybrid per-cell & lists pricing
    log2_V = log2_V_sub0 + 12.0 * d_opt
    log2_Vb = PQ["log2_Vb"] + 6.0 * d_opt
    log2_Na = PQ["log2_Va"] + 6.0 * d_opt
    cells_use = 2.0 ** log2_V; V_b_max = 2.0 ** log2_Vb; na_use = 2.0 ** log2_Na
    GB = T4["GB"]; g00 = GB[0][0]; det_GB = _det_mat(GB)
    S5 = [[GB[r][c] - GB[r][0] * GB[0][c] / g00 for c in range(1, 6)] for r in range(1, 6)]
    det_S5 = _det_mat(S5)
    Rb_max = (V_b_max * sqrt(det_GB) / ((pi**3) / 12.0)) ** (1.0 / 6.0)
    L_max = 2.0 * Rb_max / sqrt(g00)
    n5_use = ((4.0 * pi**2) / 15.0) * (Rb_max ** 5) / sqrt(det_S5)
    L_b_avg = V_b_max / n5_use
    score_arena_log2 = log2((L_max + 2.0) * na_use * 2.0 + na_use * 2.0)

    kap = PQ["kap"]
    n_r, n_s = 400, 32
    tot_len = tot_act = tot_c2 = 0.0
    for i in range(n_r):
        r2 = (i + 0.5) / n_r; w5 = 2.5 * (r2 ** 1.5) / n_r; Lb = L_max * sqrt(1.0 - r2)
        tot_len += w5 * Lb
        act_s = sum((min(1.0, (1.0 - (r2 + (1.0 - r2) * ((j + 0.5) / n_s) ** 2)) / kap) ** 3) for j in range(n_s)) / n_s
        tot_act += w5 * Lb * act_s
        tot_c2 += w5 * (Lb + w_sb - 1.0) * (min(1.0, (1.0 - r2) / kap) ** 3)
    frac_reg_analytic = tot_act / tot_len
    cr2 = tot_c2 / tot_act
    avg_sb = 1.0 + max(0.0, L_b_avg - 1.0) / w_sb
    bmw = 1.0 + max(0.0, L_b_avg - 1.0) / (64.0 * w_sb)

    # Cofactorization: LIVE norm-histogram-weighted Retention survivors + analytical Stage-B kernels!
    c_gb = max(C_GBATCH, gbatch_cost(pbits(0, lf), leaf=9))
    c_b2 = max(C_BATCH2, gbatch_cost(pbits(0, lg), c_leaf=2**19))
    cc = (c_gb, C_TDIV + C_BATCH + c_b2)
    cof_sum = sf_sum = s0_sum = 0.0
    for lnf0, lng0, cnt in bins:
        nf_raw = lnf0 + 48.0 * d_opt - CONTENT_BITS; ng_raw = lng0 + 12.0 * d_opt
        ng_q = round(ng_raw * 2.0) / 2.0
        sf_b = (8.0 / 9.0) * cost_obj.survivor_upper(ng_q, sl, 0.25) + (1.0 / 9.0) * cost_obj.survivor_upper(ng_q - 3.0, sl - 3.0, 0.25)
        s0_b = (8.0 / 9.0) * cost_obj.survivor_upper(ng_q, sl, 0.25, sb, lg, KMAX) + (1.0 / 9.0) * cost_obj.survivor_upper(ng_q - 3.0, sl - 3.0, 0.25, sb, lg, KMAX)
        q_f_b = min(1.0, Q_F_SAFETY * rho((nf_raw + 1.58) / lf))
        cof_sum += cnt * (sf_b * (C_NORMS + C_SB) + s0_b * (cc[0] + q_f_b * cc[1]))
        sf_sum += cnt * sf_b; s0_sum += cnt * s0_b
    cof_b = cof_sum / n_in; sf = sf_sum / n_in; s0 = s0_sum / n_in

    # Band 1 [2^14, 2^43): Pointwise b-lookup (exact Q_b(b) -> exact 256-shell cutoff M(b))
    rep_lo = log(p_cut / 14.0) * (C_HIT_LO * 1.004 + 0.06) + LI_PAD * (2.0**log2_li(p_cut) - 2.0**log2_li(14.0)) * V_b_max * C_LOOK_LO / cells_use
    # Band 2 [2^43, 2^53.5): Exact 64-bit bitmap words + Poisson occupied SBs + exact 6D cr2 candidates
    probe = LI_PAD * (2.0**log2_li(pm) - 2.0**log2_li(p_cut)) * n5_use * bmw * C_PROBE / cells_use
    occ_tot = 0.0; lp_o = p_cut
    while lp_o < pm - 1e-9:
        lp_o2 = min(lp_o + 0.1, pm)
        dp_o = LI_PAD * (2.0**log2_li(lp_o2) - 2.0**log2_li(lp_o))
        occ_tot += dp_o * n5_use * avg_sb * (1.0 - exp(-na_use * w_sb / (2.0**lp_o)))
        lp_o = lp_o2
    occ_sb = (occ_tot / cells_use) * C_OCC_SB
    hi_ap = log(pm / p_cut) * (C_CAND * cr2 + C_HIT_HI)

    # Band 3 [2^53.5, 2^58.0): 100% SELF-CONTAINED 7D Fincke-Pohst coset walk with EXACT Gaussian-volume integral!
    S_avg = V_b_max * PQ["frac_reg"]
    hits_top = log(lx2 / pm)
    dp_top = LI_PAD * 0.5 * (2.0**log2_li(lx2) - 2.0**log2_li(pm))
    La_long = 2.0 ** OP["rmass"]; La_short = 2.0 ** 2.50
    f_hl = OP["fge16"] * La_long; f_hs = (1.0 - OP["fge16"]) * La_short
    w_fp_sum = nph_fp_sum = 0.0
    for ii in range(400):
        lp_b3 = pm + (lx2 - pm) * (ii + 0.5) / 400.0
        wt_b3 = 1.0 / lp_b3  # exact d(ln ln p) hit-density weight!
        nph_p = (f_hl * fp7_gauss_nph(V_b_max * La_long / (2.0 ** lp_b3))
                 + f_hs * fp7_gauss_nph(V_b_max * La_short / (2.0 ** lp_b3))) / (f_hl + f_hs)
        w_fp_sum += wt_b3; nph_fp_sum += wt_b3 * nph_p
    fp_nodes_int = nph_fp_sum / w_fp_sum
    fp_vr_exact = fp7_vr()
    c_fp_walk = hits_top * (FP_RUN_OVERHANG * (C_FP_NODE * fp_nodes_int + C_FP_LEAF * FP_LEAVES_PER_HIT))
    c_babai_init = (2.0 * dp_top) * C_RUN_INIT * (OP["fge16"] / (S_avg * La_long) + (1.0 - OP["fge16"]) / (S_avg * La_short))
    bkt_top = hits_top * (C_BKT_WRITE + C_BKT_XFER + C_BKT_APPLY)
    per_cell = rep_lo + probe + occ_sb + hi_ap + c_fp_walk + c_babai_init + bkt_top + C_SCAN_ZERO + cof_b
    coll_u = log2(cells_use * per_cell * OP["work_allow"] / UNIT)

    # Lists build at C_LIST_A = 35.0
    n_lo = LI_PAD * 0.5 * (2.0**log2_li(p_cut) - 2.0**log2_li(14.0))
    sp_lo = LI_PAD * 0.5 * (2.0**(2.0 * p_cut)) / (2.0 * p_cut * LN2 - 1.0)
    ins_lo = (n_lo + 1.0) * (C_LIST_A * na_use) + C_LBUCKET * LIST_NS * sp_lo
    byt_lo = 16.0 * na_use * n_lo + 8.0 * LIST_NS * sp_lo
    n_hi = LI_PAD * 0.5 * (2.0**log2_li(pm) - 2.0**log2_li(p_cut))
    sp_hi = LI_PAD * 0.5 * (2.0**(2.0 * pm)) / (2.0 * pm * LN2 - 1.0)
    ins_hi = n_hi * (C_LIST_A * na_use) + C_LBUCKET * (sp_hi / w_sb)
    byt_hi = 16.0 * na_use * n_hi + 8.125 * (sp_hi / w_sb)
    ins_pre = na_use * (1000.0 + 12.0 * log2_Na)
    lists_u = log2((ins_lo + ins_hi + ins_pre) / UNIT)
    ram_log2 = log2(byt_lo + byt_hi + 8.0 * cells_use * hits_top + 2.0 * na_use)

    # 4. 3-Scale W=200 Real-C Ladder (2^20 -> 2^21 -> 2^22) & 4-Parameter Deterministic Hard-Cap LA!
    def stats_of(key):
        rs = FILT[key]; n = len(rs)
        def ms(idx):
            vs = [x[idx] for x in rs]; m = sum(vs) / n
            s2 = sum((v - m) ** 2 for v in vs) / (n - 1) if n > 1 else 0.0
            return m, sqrt(s2 / n), max(vs)
        sw, sw_se, sw_mx = ms(1); nl, nl_se, nl_mx = ms(2); k5, _, k5_mx = ms(3); tr, _, tr_mx = ms(4); uf, _, _ = ms(5); eu, _, eu_mn = ms(6)
        k5_exp_mean = eu * C_UNIT_5 + (1.0 - eu) * C_GEN_1
        return dict(n=n, N_sw_U=sw, N_sw_se=sw_se, N_sw_max=sw_mx, N_l_U=nl, N_l_se=nl_se, N_l_max=nl_mx,
                    k5_ev=k5, k5_exp=k5_exp_mean, tab_r=tr, tab_r_max=tr_mx, unit_frac=uf, exp_unit=eu)

    cur20 = stats_of("m184_z6_60.0_59.0_e0.95_20")
    cur21 = stats_of("m184_z6_60.0_59.0_e0.95_21")
    cur22 = stats_of("m184_z6_60.0_59.0_e0.95_22")
    assert cur20["n"] == 4 and cur21["n"] == 8 and cur22["n"] == 4

    # Enforce ALL FOUR post-filtering hard caps deterministically (abort if ANY exceeded!)
    nu_sw_cap = NU_SW_CAP * cap_scale
    nu_l_cap = NU_L_CAP * cap_scale
    k5_cap = K5_CAP
    tab_r_cap = TAB_R_CAP
    p_cap_h3 = P_CAP_H3  # Explicit exploratory heuristic premise H3 (no sample-variance Chebyshev!)

    pf_rl_cap = min(1.0, nu_l_cap / nu_sw_cap)
    log2_Ns_cap = lU + log2(nu_sw_cap); log2_N_cap = lU + log2(nu_l_cap)
    n_sm = 5
    # N_l * c_prod = N_sw * (W_ROW * k5_cap * tab_r_cap + c_tail_z6), strictly monotonic in (N_sw, N_l, k5, tab_r, W_ROW)!
    c_prod = (W_ROW * k5_cap * tab_r_cap + c_tail_z6) / pf_rl_cap

    def bw_steps(lN, m, n): return (1.0 + BW_SIGMA) * (2.0**lN + n_sm) * (1.0 / m + 1.0 / n) + 6.0
    def gen_at_cap(lN_cap, m, n):
        L_cap = bw_steps(lN_cap, m, n); k = m + n
        return log2(4.0 * (k**3) * log2(L_cap) * (fft_F(18.0 * L_cap) + 640.0 * L_cap) / UNIT)

    best_la = None
    for m in (256, 512, 768, 1024):
        for n in (20, 25, 30, 35, 40):
            N = 2.0 ** log2_N_cap; L_st = bw_steps(log2_N_cap, m, n); s = N * (1.0 / m + 1.0 / n)
            kap_bw = L_st / s; hor = max(s, (1.0 + BW_SIGMA) * (N + n_sm) / n + 4.0) / s
            per_n2 = (c_prod + C_PROJ) * (n / m + 1.0) * kap_bw + (1.0 / m + 1.0 / n) * (c_prod + (n_sm * C_FULL + (n - n_sm) * C_Z64) + C_RED9) * hor
            per_n1 = (m + n) * C_RED9
            la_u = log2((N * N * per_n2 + N * per_n1 + PREP_OWNER + PREP_GUARD + PREP_JOINT + LA_RESERVE + 2.0**20) / UNIT)
            gen_u = gen_at_cap(log2_N_cap, m, n)
            both = log2(2.0**la_u + 2.0**gen_u)
            if best_la is None or both < best_la["both"]:
                ncap = 2.0 ** log2_Ns_cap; nslots = ceil(n / 5.0)
                mem_1own = log2(ncap * (160.0 * nslots + 160.0 * nslots + W_ROW * 8.0) + (2.0**55.04) * nslots + 65.0 * (2.0**53))
                best_la = dict(la=la_u, gen=gen_u, both=both, m=m, n=n, mem_1own=mem_1own, log2_N_cap=log2_N_cap, log2_Ns_cap=log2_Ns_cap, c_prod=c_prod)

    pkey = (lf, lg, exc)
    if pkey not in _PEEL_CACHE:
        _PEEL_CACHE[pkey] = _peel_fixed_point(lf, lg, exc)
    peel = _PEEL_CACHE[pkey]

    # 5. LIVE 62-ROW LB=59 DESCENT RECURRENCE (T0 at LB=59, ZERO stored two_log or log2_EN!)
    rels_low = lr - 2.25 * se_bits
    dkey = (round(rels_low, 2), round(peel["h"], 4))
    if dkey not in _DESC_CACHE:
        _DESC_CACHE[dkey] = DESC.evaluate_descent(
            D["descent"], rho, fft_F, gbatch_cost, pbits,
            log2_rels=rels_low, peel_h=peel["h"], nf_main=405.0, ng_main=469.0,
        )
    desc = _DESC_CACHE[dkey]

    # 6. Remaining ledger rows (FULL un-discounted Schirokauer sm_u & self-contained 7D LLL bases_u!)
    rels_up = lr + 2.25 * se_bits; u_log = log2(UNIT); lb_max = max(lf, lg)
    canon_u = log2_Na + 42.8 - u_log
    r_bits = R_BN254.bit_length()
    bit_f = 48*48*512 + 7*4*36*512 + 2**17 + 3*48*1024 + 2*48*6*8*40 + 48*1024
    bit_g = 12*12*512 + 1*6*6*512 + 2**15 + 3*12*1024 + 2*12*6*8*40 + 12*1024
    sm_u = rels_up + log2(48 * r_bits * bit_f + 12 * r_bits * bit_g) - u_log  # ZERO 1/10 discount!
    misc_u = log2(2.0**(rels_up + 30.4) + 2.0**(rels_up + 14.0)) - u_log
    pi_b = log2(1.25506 * (2.0**lb_max) / (lb_max * LN2))
    fb_u = pi_b + 23.4 - u_log
    bases_u = log2(dp_top * C_LLL7 / UNIT)
    split_u = log2(sum(split_parts(2.0**rels_up, lb=lb_max).values())) - u_log + 1.0
    L_gt = ceil(lb_max - 6.0 + 0.53)
    pg = sum((2.0**(lb_max - 6.0 + 0.53 - j - 1)) * (fft_F(2.0**(j + 1)) + 3.0 * (2.0**(j + 1))) for j in range(L_gt))
    gtree_u = log2(2.0 * pg) - u_log

    rows = dict(
        collection=coll_u, chi_setup=65.712288, la=best_la["la"], gen=best_la["gen"],
        lists=lists_u, two_logs=desc["two_log"], canon=canon_u, sm=sm_u, misc=misc_u,
        fb=fb_u, bases=bases_u, split=split_u, gtree=gtree_u,
        cert=50.000001, verify=11.000001, wrapper=53.712288,
    )
    if track == "bn254-pairing-inversion":
        rows["output"] = log2(G2_OUTPUT_BOUND / UNIT)

    exp_run = log2(sum(2.0**v for v in rows.values()))

    # 7. Strict Additive Union Bound (FIX-U) + Strict RT + Clean H3 Hard-Cap Bound!
    ncap = ceil(2.0 ** best_la["log2_Ns_cap"])
    fail_terms = dict(
        main_matrix=4.0 * ncap / (2.0**64),
        filt_cap_h3=p_cap_h3,
        split_abort=(2.0**split_u) / (2.0**SPLIT_CAP - (2.0**SPLIT_RES) / UNIT),
        log_phase_cap=2.0 * (2.0**(desc["two_log"] - 1.0 - 112.0)),
        node_abort=2.0 ** desc["abort_log2"],
        repeat_ideal=2.0 ** desc["repeat_log2"],
        group_order=2.0 / R_BN254,
        numerical=2.0 ** (-48.0),
    )
    succ_exact = 1.0 - sum(fail_terms.values())
    succ_decl = f"{floor(succ_exact * 1e6) / 1e6:.6f}"
    s_frac = F(succ_decl)
    K_gate = int(-(-(F(39, 100) * (2**64)) // s_frac))
    gate_passes = (F(K_gate, 2**64) * s_frac >= F(39, 100))
    gate_prev_fails = (F(K_gate - 1, 2**64) * s_frac < F(39, 100))
    exp_run_decl = ceil(exp_run * 1e9) / 1e9
    claim_exact = log2(1.0 + float(F(K_gate, 2**64)) * (2.0**exp_run_decl))
    claim_num = ceil(claim_exact * 1e9)
    claim_decl = f"{claim_num / 1e9:.9f}"

    return dict(
        track=track, use_db=use_db, c_tail_z6=c_tail_z6, cap_scale=cap_scale,
        op=OP,
        pool=dict(n_total=PQ["n_total"], n_in=n_in, n_prim_in=PQ["n_prim_in"], frac_reg=PQ["frac_reg"], frac_reg_analytic=frac_reg_analytic),
        geometry=dict(log2_U=lU, log2_V=log2_V, se_bits=se_bits, log2_Vb=log2_Vb, log2_Na=log2_Na,
                      log2_N5=log2(n5_use), L_b_avg=L_b_avg, L_max=L_max, avg_sb=avg_sb, bmw=bmw, cr2=cr2,
                      nf_rel=nf_rel, ng_rel=ng_rel, r_g=r_g, sf=sf, s0=s0,
                      score_arena_log2=score_arena_log2, ram_log2=ram_log2),
        per_cell=dict(total=per_cell, rep_lo=rep_lo, probe=probe, occ_sb=occ_sb, hi_ap=hi_ap,
                      fp_walk=c_fp_walk, fp_nodes_int=fp_nodes_int, fp_vr_exact=fp_vr_exact,
                      babai_init=c_babai_init, bkt_top=bkt_top,
                      scan_zero=C_SCAN_ZERO, cof=cof_b),
        la=dict(bw=[best_la["m"], best_la["n"]], log2_N_cap=best_la["log2_N_cap"], log2_Ns_cap=best_la["log2_Ns_cap"],
                nu_sw_cap=nu_sw_cap, nu_l_cap=nu_l_cap, k5_cap=k5_cap, tab_r_cap=tab_r_cap,
                pf_rl_cap=pf_rl_cap, c_prod=c_prod,
                mem_1own_log2=best_la["mem_1own"], fits_1own=(best_la["mem_1own"] < 64.0),
                p_cap_h3=p_cap_h3, p_cap_h3_decl=f"{p_cap_h3:.6f}",
                ladder_monotonic=(cur20["N_sw_U"] > cur21["N_sw_U"] > cur22["N_sw_U"] and cur20["N_l_U"] > cur21["N_l_U"] > cur22["N_l_U"]),
                all_runs_below_4caps=(max(cur20["N_sw_max"], cur21["N_sw_max"], cur22["N_sw_max"]) < nu_sw_cap
                                      and max(cur20["N_l_max"], cur21["N_l_max"], cur22["N_l_max"]) < nu_l_cap
                                      and max(cur20["k5_exp"], cur21["k5_exp"], cur22["k5_exp"]) < k5_cap
                                      and max(cur20["tab_r_max"], cur21["tab_r_max"], cur22["tab_r_max"]) < tab_r_cap),
                filtsim_cur20=cur20, filtsim_cur21=cur21, filtsim_cur22=cur22, peel=peel),
        descent=desc,
        rows=rows,
        exp_run=exp_run, exp_run_decl=f"{exp_run_decl:.9f}",
        fail_terms=fail_terms,
        succ_exact=succ_exact, succ_decl=succ_decl,
        K_gate=K_gate, gate_passes=gate_passes, gate_prev_fails=gate_prev_fails,
        claim_exact=claim_exact, claim_decl=claim_decl, score_numerator=str(claim_num),
    )
