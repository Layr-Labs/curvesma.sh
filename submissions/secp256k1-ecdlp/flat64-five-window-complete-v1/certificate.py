"""Unified exact-rational certificate for the Card65 five-window search and Cards66/67 repairs.
Computes every probability, Eisenstein 5-window bound, Lemma R/V bound, Phase-0 mining bound,
cap tail, 128-bit fixup meter, explicit scalar generic/callback service surcharge, 27-row E_RUN,
peak memory, preprocessing, and exact rational 15-decimal-place claim without floating-point
dependence.
"""
import functools, math, time
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, getcontext
from fractions import Fraction as Fr

t0 = time.time()
rn, en, dc = range, enumerate, dict
FAILED = []

def check(tag, ok, info=""):
    if not ok:
        FAILED.append((tag, info))
    return bool(ok)

def lg(x):
    f = Fr(x)
    return math.log2(f.numerator) - math.log2(f.denominator)

def cd(a, b):
    return -(-a // b)

def ceilf(x):
    x = Fr(x)
    return -(-x.numerator // x.denominator)

def rat(x):
    x = Fr(x)
    return f"{x.numerator}/{x.denominator}"

def dec_up(x, places):
    x = Fr(x); s = 10**places; n = ceilf(x * s)
    return ("-" if n < 0 else "") + f"{abs(n) // s}.{abs(n) % s:0{places}d}"

def dec_dn(x, places):
    x = Fr(x); s = 10**places; n = (x.numerator * s) // x.denominator
    return ("-" if n < 0 else "") + f"{abs(n) // s}.{abs(n) % s:0{places}d}"

def expi(x, terms=60, bits=400):
    assert abs(x) <= 1
    n, d, S = x.numerator, x.denominator, 1 << bits
    lo = hi = 0; I = J = S
    for k in rn(terms + 1):
        if k:
            a, b = sorted((I * n, J * n))
            I, J = a // (d * k), cd(b, d * k)
        lo += I; hi += J
    e = abs(n) * J // (d * (terms + 1)) + 2
    return Fr(lo - e, S), Fr(hi + e, S)

def _ln2():
    a = b = Fr(0)
    for k in rn(1, 200):
        t = Fr(1, (2 * k - 1) * (3 ** (2 * k - 1)))
        a += 2 * t; b += 2 * t
    rem = Fr(2, 399 * (3 ** 397))
    return a, b + rem

LN2 = _ln2()

@functools.lru_cache(None)
def P2(x):
    x = Fr(x); k = math.floor(x); f = x - k
    s = Fr(1 << k) if k >= 0 else Fr(1, 1 << -k)
    if f == 0:
        return s, s
    l0, l1 = expi(f * LN2[0], 55)[0], expi(f * LN2[1], 55)[1]
    return l0 * s, l1 * s

def Z2(n, d):
    return P2(Fr(n, d))[0]

def least_grid_up(val, scale):
    val = Fr(val); getcontext().prec = 80
    st = (Decimal(val.numerator).ln() - Decimal(val.denominator).ln()) / Decimal(2).ln() * scale
    c_ = Fr(int(st.to_integral_value(rounding=ROUND_CEILING)), scale)
    for _ in rn(400):
        if P2(c_)[0] >= val and P2(c_ - Fr(1, scale))[1] < val:
            return c_
        c_ += Fr(1, scale) if P2(c_)[0] < val else -Fr(1, scale)
    raise AssertionError("least_grid_up failed")

def enh(o):
    o = Fr(o)
    if o <= 0: return Fr(1)
    J, mv = 0, o
    while mv > 1:
        mv /= 2; J += 1
    mlo = Fr((mv.numerator << 72) // mv.denominator, 2**72)
    I = expi(-mlo, 18, 84)[1]
    for _ in rn(J):
        I = I * I
        I = Fr(cd(I.numerator << 96, I.denominator), 1 << 96)
    return I

# 1. Curve & topology constants
p = 2**256 - 2**32 - 977
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
C0 = 2**32 + 977
BETA_P = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)
LAM_R = int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72", 16)
FA = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16)
FB = int("e4437ed6010e88286f547fa90abfe4c3", 16)
C = 1000

NODES = 2**58
N_BANK = NODES  # 2^48 owners, each with 1024 bank endpoints
M_G = 1024
N_CLUST = NODES // M_G  # 2^48 = 2^33 (s) * 2^15 (c' mod 2^15)
K_S = 33
D1, D2, D3, D4 = 40, 42, 64, 22
F1, F2 = 1 << D1, 1 << D2
N_REL, N_HOME = F1, F1 * F2
ALLN = NODES + N_REL + N_HOME
MX, YPC, YROUNDS = 6, 33, 63520
L_CHAIN = YROUNDS * YPC
Q2 = MX * L_CHAIN
KEYBITS = 64 + D4
KBF = 700
K1, K2B = 128, 256
K = K1
HB = 50
THRESH = Fr(39, 100)
LEAVES = 6

# 2. Exact AArch64 kernel instruction counts (verified by experiment.py execution!)
# FWD: base 1437/24 minus 1/Y (283cdacf) = 1433/24
FWD_INC, FWD_FIX, FWD_REM = Fr(1433, 24), 9, (8, 366, 723, 1083)
def fwd_count(q2):
    assert q2 % 6 == 0
    m, rem = divmod(q2 // 6, 4)
    return 13 if q2 == 0 else FWD_FIX + FWD_REM[rem] + 1433 * m

# BWD: exact measured 18355 = 24 + 18331 for 24 incidences (4 Ys) => BWD_INC = 18331/24!
BWD_INC, BWD_FIX, BWD_REM = Fr(18331, 24), 24, (0, 4585, 9167, 13748)
def bwd_count(q2):
    assert q2 % 6 == 0
    n = q2 // 6
    return BWD_FIX + 18331 * (n // 4) + BWD_REM[n % 4]

CH, B_CAP, CHK = 2**22, 2**26 - 1, 2
FB_FIX = max(fwd_count(6 * c) - 6 * c * FWD_INC for c in rn(1, 5)) + BWD_FIX + max(BWD_REM[c] - 6 * c * BWD_INC for c in rn(4))
FIX_SITES_QUAD = dc(SHX=18, FDX=6, FS=2, BLK=2, RW=8)
FIX_EXTRA = dc(SHX=5, FDX=5, FS=6, BLK=8, RW=6)
EMIT = 7  # explicit tbnz #12 (no S1 permission fault)
EMIT_INC = LEAVES * EMIT
KEY_W3 = Fr(1, 1)
KFC = 10
KEC = Fr(50, 6)
KEY_INC, KEY_FIX = Fr(838, 6), 8
EXC_ENTRY = 0  # flat 64-bit analytical VA: zero exception entry!
FLUSH = 18     # FLUSH7 (17) + b (1)
FLUSH2 = 18    # hf cold stub (18)
XEMIT_HOT = 14 # +1 for tbnz #12
JPATH, JPATH2, SJOIN = 3, 2, 0
INV = 22549
INIT3_BKT, TINIT_ENTRY = 4, 3
FIN_DEST, FIN_SENT, FIN_HELD, FIN_FIX = 15, 4, 2, 2
FINH_DEST, FINH_SENT = 16, 5
H2S_ARR, H2S_ONE = 6, 1  # +1 for tbnz #13 (no S1 fault)
H2S_FULL, H2S_EMPTY, H2S_ARR_FIX = 1, 6, 18
H2S_BLK128 = 733  # exact measured full 128-record h2s count with dedicated SP source table + tbnz #13
DONE = 2; RCV, WAKE, WY, EXPECT, SWAP, CTX = 1, 5, 21, 39, 64, 0
DONE_C = DONE + RCV + WAKE  # 8
RCV1, FIN_CNT, H2S_DONE = 4, 6, 4
EXPECTH = 78; RDONE = 1 + H2S_DONE; DONE_C1 = DONE + RCV1 + WAKE  # 11
DLV = 6; MSG_LINK = DONE_C + DLV + CTX  # 14
FINAL = FIN_DEST + FIN_SENT + DONE_C1 + H2S_ARR_FIX + FIN_CNT  # 54
DLVHB = 272    # direct-R2 K2B=256 full-block dlvh non-HCHGI upper bound (measured: 264 insns)
BULK_H = 48    # explicit scalar bulk_h_receive per-block wrapper upper bound (measured: 24 insns)
def dsg(n, m):
    if m == 0: return 1
    I = n % 2; o = m - I
    return 5 + I * 5 + (o // 2) * 6 + (o % 2) * 4 - (o == 0)
def dlvh_copy(o, I, J=K2B):
    if o == 0 and I == J: return DLVHB
    return 23 - (o > 0) + dsg(0, o) + dsg(o, I) + dsg(o + I, J - o - I)
DLVHS = max(dlvh_copy(a, b) for a in rn(K2B + 1) for b in rn(K2B + 1 - a) if (a, b) != (0, K2B))
BLOCKH = FLUSH2 + EXC_ENTRY + DONE_C + DLVHB + BULK_H  # 18 + 0 + 8 + 272 + 48 = 346
FINALH = FINH_DEST + FINH_SENT + DONE_C + DLVHS + MSG_LINK + (DLVHS - DLVHB) + BULK_H

check("c1_h2s", H2S_ARR == 6 and FLUSH == 18 and FLUSH2 == 18 and EXC_ENTRY == 0)
check("c2_key", KEY_INC == 8 + 48 + KEY_W3 + 6 * (4 + EMIT) + KFC + Fr(1, 3) - 2 + KEC and KEY_FIX == 8)
check("c3_fb", all(fwd_count(6 * n) + bwd_count(6 * n) <= FB_FIX + 6 * n * (FWD_INC + BWD_INC) for n in rn(1, 9)))

# 3. 5-window Eisenstein K' lattice & Lemma F/K/R bounds
NCAND = dc(A=8, B=6, EB=7)
EISR3 = {k: 19 + 31 * v for k, v in NCAND.items()}
KCL = dc(E=("A",) * 5, F=("A",) * 3 + ("B",) * 2, G=("A",) * 3 + ("B",) * 2)
def tbldk_count(n, cl): return 4268 + n * (487 + 6 * EISR3[cl])
def tbldf_count(n, v): return 4 + n * (6 + 46 * v)

ROWS_52 = 2**52 + 4 * (2**51)
ROWS_F = 3 * (2**45) + 2 * (2**44)
XTBLD_P1_ONE = 5 * (4 + 2 * 11) + 20 * 2 * ROWS_52
XTBLD_FWD = 2 * (5 * 26 + 200 * ROWS_F)
XTBLD_E = tbldk_count(2**52, "A") + 4 * tbldk_count(2**51, "A")
XTBLD_G = tbldk_count(2**52, "A") + 2 * tbldk_count(2**51, "A") + 2 * tbldk_count(2**51, "B")
XTBLD_F = 3 * tbldk_count(2**45, "A") + 2 * tbldk_count(2**44, "B")
XTBLD_ALL = 6 * XTBLD_P1_ONE + XTBLD_FWD + XTBLD_E + XTBLD_F + XTBLD_G
TBLD = tbldf_count(2**52, 1) + tbldf_count(2**51, 4)

def _cvg():
    ds, u, v = [], BETA_P, p
    while u: ds += [v // u]; u, v = v % u, u
    o, qm, qk, pm, pk = [], 0, 1, 1, 0
    for d in ds:
        qn, pn = d * qk + qm, d * pk + pm
        o += [(qn, qn * BETA_P - pn * p, d)]
        qm, qk, pm, pk = qk, qn, pk, pn
    return o
CVG = _cvg()

def eis_cls(X2):
    for k in rn(1, len(CVG) - 1):
        (qk, dk, _), (qm, dm, _), (qn, _, a) = CVG[k], CVG[k - 1], CVG[k + 1]
        if qk + qm <= X2 - 1 < qn + qk:
            g = max(abs(dk), abs(dm) - (min((X2 - 1 - qm) // qk, a) - 1) * abs(dk))
    bb = None
    for k in rn(1, len(CVG)):
        (qk, dk, _), (qm, dm, _) = CVG[k], CVG[k - 1]
        if qk + qm > X2: break
        if (dk > 0) == (dm > 0): continue
        for j in rn((X2 - qk - qm) // qk + 1):
            dw, qw = dm + j * dk, qm + j * qk
            if (dw > 0) == (dk > 0) and dw: break
            if bb is None or abs(dk) + abs(dw) < bb[0]:
                bb = (abs(dk) + abs(dw), (-dk, qk), (-dw, qw))
    _, U, W = bb
    if U[0] * W[1] < U[1] * W[0]: U, W = W, U
    lo, hi = min(0, U[0]) + min(0, W[0]), max(0, U[0]) + max(0, W[0])
    o = (-(lo + hi) // 2, (X2 - U[1] - W[1]) // 2)
    C_ = (g + 1) // 2
    cs = sorted(((i * U[0] + j * W[0], i * U[1] + j * W[1]) for i in rn(-4, 5) for j in rn(-4, 5)
                 if o[0] + lo + i * U[0] + j * W[0] <= C_ and o[0] + hi + i * U[0] + j * W[0] >= -C_
                 and o[1] + i * U[1] + j * W[1] < X2 and o[1] + U[1] + W[1] + i * U[1] + j * W[1] > 0),
                key=lambda t: (abs(t[0]) + 64 * abs(t[1]), t))
    return dc(X2=X2, C=C_, U=U, W=W, o=o, c=cs, kS=(W[0] * o[1] - W[1] * o[0]) % p, kT=(U[1] * o[0] - U[0] * o[1]) % p)

def eis_red(c, v):
    U, W, o = c["U"], c["W"], c["o"]
    S = (W[1] * v + c["kS"]) % p
    T = (c["kT"] - U[1] * v) % p
    f1, f2 = o[0] + (S * U[0] + T * W[0]) // p, o[1] + (S * U[1] + T * W[1]) // p
    for l1, l2 in c["c"]:
        if 0 <= f2 + l2 < c["X2"] and abs(f1 + l1) <= c["C"]:
            return f1 + l1, f2 + l2
    return f1, f2

def eis_cover(c):
    U, W, (o1, o2), X2, C_ = c["U"], c["W"], c["o"], c["X2"], c["C"]
    V = [(o1, o2), (o1 + U[0], o2 + U[1]), (o1 + U[0] + W[0], o2 + U[1] + W[1]), (o1 + W[0], o2 + W[1])]
    E_ = list(zip(V, V[1:] + V[:1]))
    ys = sorted({y for y in [v[1] for v in V] + [e for _, l in c["c"] for e in (-l, X2 - l)] if o2 <= y <= o2 + U[1] + W[1]})
    def xs(y):
        z = [x0 + Fr(y - y0) * (x1 - x0) / (y1 - y0) for (x0, y0), (x1, y1) in E_ if y0 != y1 and min(y0, y1) <= y <= max(y0, y1)]
        return min(z), max(z)
    for ya, yb in zip(ys, ys[1:]):
        cm = []
        for a, b in sorted((-C_ - l1, C_ - l1) for l1, l2 in c["c"] if -l2 <= ya and yb <= X2 - l2):
            if cm and a <= cm[-1][1] + 1: cm[-1] = (cm[-1][0], max(cm[-1][1], b))
            else: cm += [(a, b)]
        (La, Ra), (Lb, Rb) = xs(ya), xs(yb)
        if not any(a <= min(La, Lb) and max(Ra, Rb) <= b for a, b in cm):
            return False
    return True

def isum(*rs): return sum(x[0] for x in rs), sum(x[1] for x in rs)

X2A = 45402545270420317080910594221470922860
X2B = 64502973549206556628585045361533709079
EIS = dc(A=eis_cls(X2A), B=eis_cls(X2B), EB=eis_cls(2**128 - 6 * (X2A - 1)))
check("c7_eis", all(eis_cover(EIS[k]) and len(EIS[k]["c"]) == NCAND[k] for k in EIS)
      and 5 * (X2A - 1) + (EIS["EB"]["X2"] - 1) < 2**128
      and 3 * (X2A - 1) + 2 * (X2B - 1) < 2**128)

C_E = sum(EIS[c]["C"] for c in KCL["E"]) + EIS["EB"]["C"]
C_FG = sum(EIS[c]["C"] for c in KCL["F"])
E0R, E1R = (0, 1 + Fr(1, 2**30)), (Fr(-5, 2), Fr(5, 2) + Fr(1, 2**180))

def _half(lo, hi, ti, e):
    t0, t1 = (lo >> 128, hi >> 128) if ti else (0, 0)
    s = [m * x / 2**128 for m in (lo, hi) for x in e]
    return (min(t0, 0) - 1, max(t1 - 1, 0) + 5 + ti), (math.floor(min(s)), 2 + math.ceil(max(s)))

def _kt(C_):
    es = (E0R[0] + E1R[0], E0R[1] + E1R[1])
    p1 = _half(-C_, C_, 1, E0R)
    p2 = _half(0, 2**128 - 1, 0, E1R)
    p3 = _half(-C_ - 2**128, C_, 1, es)
    return (isum(p1[0], p2[0]), isum(p1[1], p2[1])), ((p3[0][0] - p1[0][1], p3[0][1] - p1[0][0]), (p3[1][0] - p1[1][1], p3[1][1] - p1[1][0]))

KT = dc(E=_kt(C_E), F=_kt(C_FG), G=_kt(C_FG))
def _table3(rs):
    rot = lambda a, b: (a, b, (-2 - a[1] - b[1], -a[0] - b[0]))
    (e0, e1), (f0, f1), (g0, g1) = rs
    al = rot(isum(e0, (-1, 0)), isum(e1, (-1, 0)))
    bp = rot(*[(-f[1] - 2, -f[0] - 1) for f in (f0, f1)])
    ga = rot(g0, g1)
    W = lambda m, t: isum(al[m % 3], bp[(m + t) % 3], ga[(m + 2 * t + 1) % 3], (0, 2))
    wr = [[W(o + i, t) for t, o in en((2, 2, 0))] for i in (0, 1)]
    nr = [(-wr[0][t][1] - wr[1][t][1] - 1, -wr[0][t][0] - wr[1][t][0] + 1) for t in rn(3)]
    return [(min(x[0] for x in r_), max(x[1] for x in r_)) for r_ in (wr[0], wr[1], nr)]

ERRH = _table3([(KT[c][0][0], KT[c][1][0]) for c in "EFG"])
ERRL = _table3([(KT[c][0][1], KT[c][1][1]) for c in "EFG"])
WDH, WDL = sum(hi - lo for lo, hi in ERRH), sum(hi - lo for lo, hi in ERRL)
WH, WL = max(hi - lo for lo, hi in ERRH) + 1, max(hi - lo for lo, hi in ERRL) + 1
ELO, EHI = min(ERRH[0][0], ERRH[1][0]), max(ERRH[0][1], ERRH[1][1])
NWIN = 2 * (EHI - ELO) + (1 - 2 * ELO) + (3 + 2 * EHI) + 2 * (EHI + 1 - ELO) + 2 * (3 * EHI + 3 - 3 * ELO)
TWB = Fr(2, r - 1)
MH_A, MH_B, ML = 20, 21, 23
D_WJ = TWB * NWIN * 2**192
D_LAB = TWB * (WDH * 2**(256 - MH_A) + WDL * 2**(256 - ML))

def lemma_r(k, c_, w, vertical):
    N_, l, P_ = 2**(k - c_) + 1, w * cd(p, 2**k), p >> (k - c_)
    A_ = w * 2**(256 - c_)
    b = ((l + FA) * (Fr(p * p, 3 * P_) + 3 * (p + P_) + N_ * (8 * l + 2 * (2 * FA - FB))) / p
         if vertical else
         (Fr(p * p, 3 * P_) + Fr(3, 2) * (p + P_) + N_ * (4 * l + 2 * FA)) * (l + 2 * FA - FB) / p)
    return b, A_, 3 * b / A_ - 1

LEMMA_R = [lemma_r(64, MH_A, WH - 1, 1), lemma_r(64, MH_B, WH - 1, 0),
           lemma_r(128, ML + 1, WL - 1, 1), lemma_r(128, ML, WL - 1, 0)]
D_LAB_SYM = 3 * TWB * sum(b for b, _, _ in LEMMA_R)
B_HA, B_HB, B_L, B_LA = 2**MH_A + WH - 1, 2**MH_B + WH - 1, 2**ML + WL - 1, 2**(ML + 1) + WL - 1
check("c9_kb", B_HA * B_HB * (B_LA * 2**128 + 2 * FA) * (B_L * 2**128 + 4 * FA - 2 * FB) < p * 2**(MH_A + MH_B + 2 * ML + 1) * (1 + Fr(1, KBF)))

KB = functools.lru_cache(None)(lambda fixed: Fr(6 * 2**(256 - fixed) * (KBF + 1), (r - 1) * KBF))
KB_BUCKET = KB(D1 + D2 + HB)
PAD_C = max(3 * EHI + 3, EHI + 1 - ELO, 1 - 3 * ELO)
WRAPT_C = max(EHI, 1 - ELO, 1 - 2 * ELO, 2 + 2 * EHI)
def _npad(i, m, d=PAD_C):
    if not 0 <= i < m: return 0
    lo = max(0, i - d)
    A_, B_ = (lo, min(m - 1, m + d - 2 * i)), (max(lo, m - d - i), min(m - 1, 2 * m + d - 2 * i))
    return sum(max(0, y - x + 1) for x, y in (A_, B_)) - (max(0, min(A_[1], B_[1]) - max(A_[0], B_[0]) + 1) if A_[0] <= A_[1] and B_[0] <= B_[1] else 0)
_LC = [(0, 0), (1, -PAD_C), (0, 2**64 - 1), (-2, 2**64 + PAD_C), (-1, 2**64 - PAD_C), (-2, 2 * 2**64 + PAD_C)]
co = sorted({0, 2**64} | {t + e for a, b in _LC for cg, g in _LC if a != cg for t in [(g - b) // (a - cg)] for e in rn(-3, 4) if 0 <= t + e <= 2**64})
NC_C = sum((_npad(s, 2**64) + _npad(e - 1, 2**64)) * (e - s) // 2 for s, e in zip(co, co[1:]))
ROWMAX132 = cd(NC_C, 1 << 42) + (len(co) - 1) * 2**64
CNT132 = cd(B_HA * B_HB * (ROWMAX132 + (1 << 22) * 3 * B_HA), 1 << 50) + (1 << 22) * 4 * B_HA * B_HB + 2 * WRAPT_C * (1 << 14) * B_HB + (1 << 22) * B_HA * (6 * WRAPT_C + 8)
DPAR = max(FA, 2 * FA - FB)
KB_C132 = Fr(6 * CNT132 * ((B_LA * (1 << 128) + 2 * DPAR) * (B_L * (1 << 128) + 2 * DPAR)), p * (r - 1))
check("c19c_c132", KB_C132 < Fr(34, 100) * KB_BUCKET)

# 4. 5-window P1_5WIN WLT, batch defects, 128-bit fixup meter
WLT = [(12, C0 + 1, 2), (12, 5, 1), (12, 6, 1), (12, 2**66, 2), (6, 2, 1), (6, 2, 1), (6, 2**66, 2), (6, 2**63 + 2, 2)]
NW = [2**127 * (2**(32 * i + 1) * (1 + Fr(1, 2**20)) + (i < 3) * 2**65 * C0) + C0 * 2**192 for i in (1, 2, 3)]
D_WLT = sum(Fr(f * w * C0 * cd(p, 2**(64 * l)), r - 1) for f, w, l in WLT) + 6 * sum(NW) / (r - 1)
FIX_BATCH = sum(FIX_SITES_QUAD.values()) * CH * 6 // 4
QP = Q2 // 2
SHX_BAD, FOLD_BAD = 2**64 * 2**33, 2**66 * C0
SITES = [(3 * QP, SHX_BAD, 36), (QP, FOLD_BAD, 36), (3 * QP, SHX_BAD, 36), (QP, FOLD_BAD, 36), (QP, FOLD_BAD, 36), (3 * QP, SHX_BAD, 36)]
D_SITES = sum(Fr(n * bad * mult, r - 1) for n, bad, mult in SITES)
D_YGEN = Fr(YROUNDS * (YPC + 34), r - 1)
D_BATCH = Fr(CH, L_CHAIN) * (D_SITES + Q2 * Fr(6 * (C0 + 1), r - 1)) + CH * D_YGEN
# Pair-product canonicality is distinct from the existing individual-inverse exclusion.
D_PAIR_INVERSE_CHAIN = Fr(12 * (C0 + 1) * QP, r - 1)
D_PAIR_INVERSE = Fr(CH, L_CHAIN) * D_PAIR_INVERSE_CHAIN
D_BATCH += D_PAIR_INVERSE
check("pair_chunk_population", D_PAIR_INVERSE == Fr(36 * (C0+1) * CH, r-1)
      and CHK*CH >= L_CHAIN)
FIX_ROW, FIX_FOLD, FIX_LK7, FIX_SUB, FIX_RW = C0 * 2**192, (2**66 + 1) * 2**128 * C0, 21 * C0 * 2**192, C0 * 2**192, C0 * 2**192
FIXUPS = [("BLK", QP, FIX_LK7, 36), ("SHX", 3 * QP, FIX_ROW, 36), ("FDX", QP, FIX_FOLD, 36), ("SHX", 3 * QP, FIX_ROW, 36),
          ("FDX", QP, FIX_FOLD, 36), ("FDX", QP, FIX_FOLD, 36), ("FS", QP, FIX_SUB, 36), ("SHX", 3 * QP, FIX_ROW, 36), ("RW", 2 * Q2, FIX_RW, 6)]
E_FIX = sum(Fr(n * bad * mult, r - 1) for _, n, bad, mult in FIXUPS) + CHK * D_BATCH * FIX_BATCH + 2 * Q2 * D_WLT
DELTA_F = 1 << 71  # 128-bit cumulative fixup meter (integrated in ytransport.s:ychunk)
EF_Q = E_FIX * C0 / (2**32 * (DELTA_F - 33 * FIX_BATCH))

# 5. Phase 0 mining + Vandermonde GF(1024) + Phase 1-5 success
RHO = Fr(r, r - 1)**2
SLF = 6 * RHO * Fr(2**256, r)
SL6 = 6 * L_CHAIN * Fr(2**256, r)
D_DROP = 2 * SLF * L_CHAIN / (B_CAP - 2 * L_CHAIN - SL6)**2
CLUST_BATCH = M_G * LEAVES * Q2
FIXED_LEAVES = 2 * CLUST_BATCH
PARTNER_PAIRS = CLUST_BATCH * (CLUST_BATCH - 1)
P1 = Fr(6 * (r - 1) + 1, r * r)
P1_GOOD = Fr(6 * (r - 1), r * r) * (1 - 3 * D_WJ - D_LAB_SYM) - 12 * D_WLT / r
K_U_CLUST = LEAVES * M_G * MX
W_ELL = 1 + 47 * (3 * D_WJ + D_LAB_SYM)
SAMEX_FAC = Fr(CLUST_BATCH - K_U_CLUST, CLUST_BATCH - 1) * Fr(r, r - 2) * W_ELL + Fr(96 * (K_U_CLUST - 1), CLUST_BATCH - 1)
KB_M, KB_D1, KB64 = KB(KEYBITS + D1 + D2), KB(D1), KB(64)
LA = Fr(6, r) / P1_GOOD
D_M = Fr(184 * 13 + 8 * 10 + 2 * 192 * 13, 2)
DSY = (3 * 192)**2 * 624 * K_U_CLUST
DWS = 48 * RHO * CLUST_BATCH * D_WLT * (2 + D_M) + DSY
SIP = 240 * RHO * KB_M
R2_TIB = 38220
CAP_H = ((R2_TIB * (2**40)) // (32 * K2B)) * K2B
CAND_CAP = 2082000000
REQ_CAP = 256 * 125 * 2**39  # 1000 PiB of request records per G on flat64 VA
RELAY_CAP = 2**49 - 2**43 + 2**41

LOG2_POOLS = 14; N_POOLS = 1 << LOG2_POOLS; R_POOL = N_POOLS // M_G; K_GEN = NODES // N_POOLS
SAME_M_NODE_PAIRS = M_G * (N_CLUST * N_CLUST * (R_POOL - 1) // R_POOL) // 2
SHARED_BOTH_CLUST_PAIRS = (N_CLUST * N_CLUST // (R_POOL * R_POOL) - N_CLUST) // 2
ACP = N_CLUST * (N_CLUST - 1) // 2
GOOD_NODE_PAIRS = SAME_M_NODE_PAIRS + SHARED_BOTH_CLUST_PAIRS * (M_G * (M_G - 1) // 2) + (ACP - SHARED_BOTH_CLUST_PAIRS) * (M_G * (M_G - 1))
YB = 256; YR_MINE = 992; PTS_MINE_CHAIN = YR_MINE * YPC
NBATCH_MINE = 256 * 23750000 * YB
PTS_PER_GEN_MINE = NBATCH_MINE * PTS_MINE_CHAIN
M_POOL = K_GEN * PTS_PER_GEN_MINE
S_MINE_TOT = NODES * PTS_PER_GEN_MINE
N_CHAINS_POOL = K_GEN * NBATCH_MINE
PPMC = PTS_MINE_CHAIN * (PTS_MINE_CHAIN - 1) // 2
P_BUCKET_MAX = Fr(6 * (1 << 63), r - 1)
P_SAME_CHAIN_POOL = N_CHAINS_POOL * PPMC * P_BUCKET_MAX
P_CHAIN_COUPLE_POOL = N_CHAINS_POOL * PPMC * (M_POOL * M_POOL + (N_CHAINS_POOL - 1) * PPMC) * (P_BUCKET_MAX**2)
GAMMA_MINE_CROSS = M_POOL * (M_POOL - 1) // 2 - N_CHAINS_POOL * PPMC
P_BDAY_HIT_MIN = Fr((r - 1) * (r - 1), (1 << 193) * r * r) - Fr(6 * (r - 1), r * r)
MU_POOL_LO = GAMMA_MINE_CROSS * P_BDAY_HIT_MIN
CS_POOL = GAMMA_MINE_CROSS * (2 * PTS_MINE_CHAIN * M_POOL) * (2 * P_BUCKET_MAX * P_BUCKET_MAX)
_P_A102 = Fr(6 * (1 << 103), r - 1)
_P_A174 = Fr(6 * (1 << (256 - D1 - D2)), r - 1)
_P_A215 = Fr(6 * (1 << (256 - D1)), r - 1)
P_CELL154_MASK_POOL = (Fr(M_POOL * (M_POOL - 1), 2) * M_POOL + 2 * N_CHAINS_POOL * PPMC * M_POOL) * P_BUCKET_MAX * _P_A102
POIS_LE3 = enh(MU_POOL_LO) * sum((MU_POOL_LO**k_) / math.factorial(k_) for k_ in rn(4))
P_FAIL_POOL = POIS_LE3 + CS_POOL + P_CELL154_MASK_POOL + P_SAME_CHAIN_POOL + P_CHAIN_COUPLE_POOL
P_FAIL_ALL_POOLS = N_POOLS * P_FAIL_POOL
LAM_ORD = N_POOLS * (M_POOL * M_POOL) * P_BUCKET_MAX
POINTWISE_X_MAX = M_POOL * M_POOL * P_BUCKET_MAX
SAME_OFF_VAND = (2 * LEAVES * MX)**2 + (M_G // R_POOL) * (LEAVES * MX)**2 + Fr(K_U_CLUST * K_U_CLUST, N_CLUST)
P_BAD_X = (65536 * sum((LAM_ORD**(d_ - 1)) * (2 * LAM_ORD / r) for d_ in (2, 3, 4))
           + 65536 * sum((LAM_ORD**(d_ - 1)) * N_POOLS * (M_POOL * M_POOL) * Fr(30 * (1 << 63), r * r) for d_ in (2, 3)))
P_X_ORBIT_ALL_POOLS = 4 * LAM_ORD * (6 * LAM_ORD / r)
P_MINE_ZERO_ALL = Fr(S_MINE_TOT, r)
P_SAME_A_ALL = N_POOLS * Fr(M_POOL * (M_POOL - 1) // 2 * 6, r - 1)
D_XMINE_ALL = NODES * NBATCH_MINE * Fr(YR_MINE * (YPC + 34), r - 1)
ECMX = N_POOLS * Fr(M_POOL * (M_POOL - 1), 2) * _P_A102
E_CAND_MINE = cd(ECMX.numerator, ECMX.denominator)
E_VIS_MINE = cd((S_MINE_TOT * (S_MINE_TOT - 1) // 2) * 6 * (1 << (256 - 131)), r - 1)
E_K2_MINE = cd(N_POOLS * (M_POOL * (M_POOL - 1) // 2) * 6 * (1 << (256 - 131)), r - 1)
_OTH_H_OUT = N_POOLS * Fr(M_POOL * (M_POOL - 1), 2) * _P_A174 * _P_A102 + 2 * PTS_MINE_CHAIN * M_POOL * _P_A102
_OTH_G_OUT = 2 * PTS_PER_GEN_MINE * M_POOL * _P_A102
_M_R_OUT = 18 * N_POOLS * Fr(M_POOL * (M_POOL - 1), 2) * _P_A215 * _P_A102
_cc = 1 << 20
_RL_CR_MINE = 36 * PTS_MINE_CHAIN * M_POOL * _P_A102
_VAR_R_OUT = 324 * N_POOLS * M_POOL * _P_A215 * (1 + RHO * RHO * PTS_MINE_CHAIN * _P_A215) * (M_POOL * _P_A102)**2 + 18 * PTS_MINE_CHAIN * _M_R_OUT
P_MINE_CAPS = Fr(1, 1 << 200) + LAM_ORD * (_OTH_H_OUT / Fr(CAND_CAP - 2 * PTS_MINE_CHAIN) + 2 * _OTH_G_OUT / Fr(REQ_CAP - 4 * PTS_MINE_CHAIN) + _RL_CR_MINE / _cc + _VAR_R_OUT / Fr(RELAY_CAP - 36 * PTS_MINE_CHAIN - _cc - _M_R_OUT)**2)
_MINE_K_CH = CAP_H // PTS_MINE_CHAIN
_MINE_Y_MEAN = S_MINE_TOT * _P_A174
_MINE_HOME_EXPONENT = 82 + Fr(3, 2) * _MINE_Y_MEAN - _MINE_K_CH
check("c19b_home_cap", _MINE_HOME_EXPONENT < -200)
check("c19b_mine_denominators", CAND_CAP > 2*PTS_MINE_CHAIN
      and REQ_CAP > 4*PTS_MINE_CHAIN
      and RELAY_CAP > 36*PTS_MINE_CHAIN+_cc+_M_R_OUT)
check("c19b_unchanged_exclusions", P_SAME_A_ALL < Fr(1, 1<<39)
      and D_XMINE_ALL < Fr(1, 1<<140) and P_MINE_ZERO_ALL < Fr(1, 1<<140))
check("c19b_mine", MU_POOL_LO > 31 and P_FAIL_ALL_POOLS < Z2(-19, 1) and P_BAD_X < Z2(-140, 1) and P_MINE_CAPS < Z2(-19, 1))

# 6. Original 22000-group power-law schedule, with 2^27 slots per bank.
G_BASE, G_EXP = 22000, 17
BANK_SLOT_BITS = 27
BANK_SLOTS = 1 << BANK_SLOT_BITS
ARCHIVE_CAP = 1024 * BANK_SLOTS

def group_sizes(q1):
    ms = [0]; o = G_BASE**G_EXP; J = q1**32
    for g in rn(1, G_BASE):
        I = J * (g**G_EXP); m = int(q1 * ((g / G_BASE)**(G_EXP / 32.0))) + 2
        while m**32 * o > I: m -= 1
        while ((m + 1)**32) * o <= I: m += 1
        ms.append(m)
    ms.append(q1)
    qs = [b-a for a,b in zip(ms, ms[1:])]
    return qs, ms

Q1GD = 256 * 106000000000  # 27,136,000,000,000 < 2^45
check("c19d_pack", Q1GD < (1 << 45) and N_CLUST == (1 << 48) and NODES == (1 << 58) and NBATCH_MINE * PTS_MINE_CHAIN < (1 << 56))
Q1G, MQ1 = group_sizes(Q1GD)
G_GROUPS = len(Q1G)
check("c19e_archive_cap", max(Q1G) <= ARCHIVE_CAP and min(Q1G) >= 1 and sum(Q1G) == Q1GD and MQ1[-1] == Q1GD)
ck = P_FAIL_ALL_POOLS + P_SAME_A_ALL + P_BAD_X + D_XMINE_ALL + P_MINE_ZERO_ALL + P_X_ORBIT_ALL_POOLS + P_MINE_CAPS
d5 = P1 * P1 + Fr(36, (r - 1)**2)
cf = KB(D1 + D2)
cz = Fr(K_U_CLUST * L_CHAIN, r - 1) * cf**2
cp = SAMEX_FAC * PARTNER_PAIRS + CLUST_BATCH * CLUST_BATCH + 96 * K_U_CLUST * K_U_CLUST
_FDIFF0 = SAME_OFF_VAND + 320 * K_U_CLUST * K_U_CLUST * ((POINTWISE_X_MAX / MX)**4) * KB_M + 2 * CLUST_BATCH * CLUST_BATCH * RHO * KB_M
cr = 1 + Fr(1, r)
RL_SAME_PL = W_ELL * Fr(r, r - 2) * (RHO * Fr(CLUST_BATCH - 1, 2) * KB_D1 * KB_M + 96 * K_U_CLUST * KB_M + SIP)
cy = [NODES * LEAVES * q * Q2 for q in Q1G]
SG_MAX = max(cy)
SPROD_MAX = max(cy[g] * (NODES * LEAVES * MQ1[g + 1] * Q2) for g in rn(G_GROUPS))
M_VG_MAX = 18 * (RHO * Fr(SPROD_MAX) * KB_D1 * KB_M + Fr(SG_MAX) * RL_SAME_PL)
TH_RELAY = RELAY_CAP - 18 * FIXED_LEAVES - RELAY_CAP * 9 // 10
check("c_relay_positive_gap", TH_RELAY > M_VG_MAX)
cn = (Fr(TH_RELAY) - M_VG_MAX)**2

def spx(q1p):
    n0p = LEAVES * q1p * Q2; Sp = NODES * n0p; gamma = GOOD_NODE_PAIRS * (n0p * n0p)
    I = gamma * P1; mu = gamma * P1_GOOD; cs = gamma * (2 * CLUST_BATCH * Sp) * d5; e_hi = enh(mu)
    d_f = SLF * q1p * (EF_Q + CHK * D_BATCH)
    m_p2 = Sp * cf + 1; v_p2 = Sp * (K_U_CLUST * cf + cz) + 1
    d_p2 = v_p2 / Fr(CAP_H - FIXED_LEAVES - m_p2 - CAP_H // 4)**2 + LA * Sp * D_WLT * cf / (CAP_H // 4)
    lam_int = RHO * (Sp + FIXED_LEAVES) * KB_M + 2 * RHO * Sp * D_WLT * KB_M * (1 + D_M)
    dw = RHO * Sp * D_WLT * (2 + D_M) + DWS
    alpha_cross = (2 * K_U_CLUST + SAMEX_FAC * (CLUST_BATCH - 1)) * Sp + cp
    other_h = RHO * (Fr(Sp) * Sp * cf * KB_M + alpha_cross * KB_M + 2 * L_CHAIN * cf) + SIP * Sp
    d_cl = (other_h + KB_M * dw * (LA * cf * Sp + 2 * FIXED_LEAVES)) / (CAND_CAP - 76)
    other_g = RHO * KB_M * (Fr(n0p) * Sp * 2 + 2 * LEAVES * MX * Sp + 96 * K_U_CLUST * n0p + 96 * LEAVES * MX * K_U_CLUST)
    d_rq = (other_g + 2 * KB_M * n0p * dw) / (REQ_CAP - 2 * L_CHAIN - 600)
    rl_cap_cr = RELAY_CAP * 9 // 10
    rl_cr = 18 * KB_M * (FIXED_LEAVES * (Sp + 2 * dw) + LA * KB_D1 * (SPROD_MAX * RHO * D_WLT * (2 + D_M) + SG_MAX * DWS))
    o_ = RHO * SPROD_MAX * KB_D1 * KB_M
    rl_fdiff = _FDIFF0 + 2 * CLUST_BATCH * Fr(Sp) * KB_M
    v_rl = 648 * (o_ * rl_fdiff + CLUST_BATCH * Fr(Sp) * RL_SAME_PL)
    d_rl = rl_cr / rl_cap_cr + v_rl / cn
    per_leaf = D_BATCH + d_f + d_p2 + d_cl + d_rq + d_rl + D_DROP
    rows_ = I * cr * (Fr(1, r) + lam_int + 2 * per_leaf)
    g_def = NODES * q1p * CHK * D_BATCH
    glob = ck + NODES * (3 * q1p + 2 * NBATCH_MINE + 1) * Fr(1, 2**254) + g_def
    lo = 1 - e_hi - cs - rows_ - glob
    return dc(mu=mu, lo=lo, cs=cs, rows=rows_, glob=glob, Sp=Sp, d_p2=d_p2, d_f=d_f, d_cl=d_cl, d_rq=d_rq, d_rl=d_rl, g_def=g_def)

TFULL = spx(Q1GD)
check("c_home_positive_gap", CAP_H-FIXED_LEAVES-(TFULL["Sp"]*cf+1)-CAP_H//4 > 0)
check("c_success_unchanged_bounds", TFULL["cs"] < Fr(1, 1<<82) and D_DROP < Fr(1, 1<<26))
LO_RUN = TFULL["lo"]
GB, GATE_INSNS = 128, 256
K_GATE = cd(THRESH.numerator * LO_RUN.denominator * (1 << GB), THRESH.denominator * LO_RUN.numerator)
P_GATE = Fr(K_GATE, 1 << GB)
GATED_LO = P_GATE * LO_RUN
check("c20_c21_gate", LO_RUN >= Fr(9995, 10000) and GATED_LO >= THRESH and Fr(K_GATE - 1, 1 << GB) * LO_RUN < THRESH)
check("c28_c30_caps", TFULL["d_f"] < Z2(-188, 10) and max(TFULL["d_cl"], TFULL["d_rq"], TFULL["d_rl"]) < Z2(-152, 10) and TFULL["glob"] < Z2(-18, 1))

def pq(q):
    n = LEAVES * q * Q2; g = GOOD_NODE_PAIRS * n * n
    return g * P1_GOOD, 2 * cr * g * P1 * (NODES * n * (K_U_CLUST * cf + cz) + 1) / Fr(CAP_H - FIXED_LEAVES - NODES * n * cf - 1 - CAP_H // 4)**2

V = [1 - spx(k)["lo"] - enh(pq(k)[0]) - pq(k)[1] for k in rn(6)]
DF = []
while V: DF += [V[0]]; V = [b - a for a, b in zip(V, V[1:])]
PDR = math.lcm(*[x.denominator for x in DF]); PRN = [int(x * PDR) for x in DF]
def phig(q):
    m, t = pq(q); e = enh(m) + t
    return e.numerator * PDR + sum(c_ * math.comb(q, k) for k, c_ in en(PRN)) * e.denominator, e.denominator * PDR
TPN = [phig(m) for m in MQ1[:-1]]
PHI = [Fr(1)] + [Fr(1) if N >= D else Fr(cd(N << 256, D), 1 << 256) for N, D in TPN[1:]]
check("c32_mono", all(a >= b for a, b in zip(PHI, PHI[1:])))

# 7. Complete 27-row cost accounting including explicit scalar generic & callback surcharges!
MSG = 256
GEN_MUL = 1 << 20  # explicit scalar callback envelope (same as 283cdacf scalar_costs.py!)
YR_FIX, YR_BLK, YR_BLK_MINE = 22558, 62249, 21636
CHAIN_ADD = 2**15
CHAIN_SET = 2 * GEN_MUL + 17 * CHAIN_ADD
Y_TRANSPORT = 129; OFFER = 576
METER_CHUNK = 24  # conservative charge for inline 128-bit meter check + adds/adc in ytransport.s:ychunk

def y_group(q):
    o = cd(q, M_G)
    J = YROUNDS * (cd(o, YB) * YR_FIX + o * YR_BLK) + o * CHAIN_SET
    I = cd(q * L_CHAIN * Y_TRANSPORT, M_G)
    card54_y = 512 * cd(q, M_G) * L_CHAIN + 16 * cd(q, M_G)
    return J + I + card54_y

XSETUP = 2**17; ENTRY_X = 18; COIN_SCALAR = 64; BATCH_SETUP = 9544; INV_CALL = 16; HEADER = 80
PROGRAM_BYTES = 1 << 20
PROGRAM_SOURCE_BYTES = 330336  # seven expanded templates plus 15 extra PE generic bindings
PROGRAM_ANALYTIC_BYTES = 262144  # explicitly described shared-loop callback construction
check("c_program_construction", PROGRAM_SOURCE_BYTES + PROGRAM_ANALYTIC_BYTES <= PROGRAM_BYTES)
COLD_EXTRA = 2**22 + 2**21  # cold flat-64 node init + SP source tables + endpoint state init
NODE_INIT = (PROGRAM_BYTES // 32) * 2 * MSG + COLD_EXTRA
# Preserve and expose the leader's protocol/fabric cold reserve. The bound
# below separately covers endpoint contexts/maps, callback/collector setup,
# source-pointer helpers and a conservative regenerated controller schedule.
# Advice/table row construction is separately charged, never hidden here.
COLD_SCHEDULE = (G_BASE + 1) * 46 * 16 * 65536
check("c_cold_schedule", COLD_SCHEDULE < (1 << 41))
COLD_EXPLICIT_BOUND = ((ALLN + N_BANK + 1) * (1 << 24)
                       + NODES * (1 << 27) + N_REL * (1 << 20)
                       + N_HOME * (1 << 24) + (1 << 41))
COLD_PROTOCOL_ALLOWANCE = 1 << 112
check("c_cold_explicit", COLD_EXPLICIT_BOUND <= COLD_PROTOCOL_ALLOWANCE)
COLD_TOTAL = (ALLN + N_BANK + 1) * NODE_INIT + COLD_PROTOCOL_ALLOWANCE

def bank_control(q):
    # Begin/seal/release: one bounded command per active bank in each pass;
    # metadata initialization is 4/active slot plus fixed call/setup. Double
    # that per-slot bound and use 4096 per command/sweep for all caller glue.
    return N_CLUST * (3 * 4096 * cd(q, BANK_SLOTS) + 8 * q + 4096) + 256
ZW = 16384; ZERO_PASS = 130; DLVH0 = 15; STOP_CHECK = 5
NBAR_ALLOW = 418009; PHASE_CTRL_G = 18304000; OUTPUT = 2**20
MULX, SQRX, ADDX = 109, 83, 27
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX
CAND_LEAVES = 6 * CHAIN_ADD + 12 * MULX
RESOLVE_RESERVE = 1 << 21
REQ_RESERVE = 4096
CAND_ITEM = 4 * 8 * 2 * MSG + 2 * CAND_LEAVES + 36 * 600 + 2 * (384 * 200 + GEN_MUL)
RESOLVE = 2 * CAND_ITEM + 1 + RESOLVE_RESERVE
REQ_SERVE = 2**10 + GEN_MUL + REQ_RESERVE
BARRIER_MSGS = 2 * (N_HOME + NODES + 2 * (N_REL - 1))
BARRIER = max(2**12 * ALLN, 7 * (40 * ALLN + BARRIER_MSGS * MSG))
EXCH_ROW = ALLN * 4 * 2**10
HR_INIT = N_HOME * (INIT3_BKT + 1 + DLVH0)
HCHG_NEW, MAIN_VISIT, UNWEIGHTED_VISIT, HCHG_VER, HCHG_LS1 = 3, 10, 12, 5, 6
SG = [NODES * LEAVES * q * Q2 for q in Q1G]
SG_FULL = NODES * LEAVES * Q1GD * Q2
PAIRM = KB_C132 + Fr(1, r)
PAIRK2 = KB(D1 + D2 + 64) + Fr(1, r)
PJ_MEAN = Fr(62, 54) + JPATH * (2 * D_WJ + Fr(1, r - 1))

WGL, WGR, WGG, WGJ, WGK2, SAMEX = [], [], [], [], [], []
for g, q in en(Q1G):
    s = SG[g]
    leaf_ = (NODES * q * (CHK * (FB_FIX + COIN_SCALAR + 41 + INV + INV_CALL + KEY_FIX + METER_CHUNK) + BATCH_SETUP)
             + NODES * y_group(q) + bank_control(q)
             + cd(NODES * q * (FWD_INC + BWD_INC + KEY_INC - EMIT_INC) * Q2, 1)
             + NODES * q * (OFFER * 1023 // 1024 + 8 * 4096))
    card54_r = (45 * (s // K1)
                + 41 * NODES * F1 + 98 * N_HOME + 8 * N_REL * N_CLUST + 64 * N_CLUST + 2 * NODES)
    route_ = (s * (EMIT * K1 + FLUSH + EXC_ENTRY + DONE_C1) // K1
              + (s // K1) * H2S_BLK128 + (s % K1) * H2S_ARR
              + NODES * (FINAL * F1 + FIN_FIX) + s * BLOCKH // K2B + N_REL * (FINALH * F2 + FIN_FIX)
              + N_HOME * EXPECTH + N_REL * (EXPECT + RDONE) + NODES * (INIT3_BKT * F1 + 1) + card54_r)
    gp = NODES * LEAVES * MQ1[g] * Q2
    samex_ = N_CLUST * q * CLUST_BATCH * (CLUST_BATCH - 1) // 2
    main_ = s * gp + s * (s - 1) // 2 - samex_
    WGL.append(leaf_); WGR.append(route_)
    WGG.append(s * HCHG_NEW + ALLN * STOP_CHECK + MAIN_VISIT * main_ * PAIRM)
    WGK2.append(HCHG_VER * main_ * PAIRK2)
    WGJ.append(s * PJ_MEAN)
    SAMEX.append(samex_)

# Phase 0 mining (charged at full MSG=256 and GEN_MUL=2^20!)
GEN_MINE_Y = YR_MINE * ((NBATCH_MINE // YB) * YR_FIX + NBATCH_MINE * YR_BLK_MINE) + NBATCH_MINE * (CHAIN_SET + 2 * COIN_SCALAR + 64)
GEN_MINE_NODE = 64 + GEN_MINE_Y + PTS_PER_GEN_MINE * XEMIT_HOT
_P_MINE_DEF = D_XMINE_ALL + P_MINE_ZERO_ALL + NODES * (2 * NBATCH_MINE + 1) * Fr(1, 2**254)
W_MINE_DET_MAX = S_MINE_TOT * CAP_H * (UNWEIGHTED_VISIT + HCHG_VER) + N_HOME * CAND_CAP * (HCHG_LS1 + RESOLVE + 2 * REQ_SERVE + 8 * MSG)
DEF_MINE = ceilf(_P_MINE_DEF * W_MINE_DET_MAX)
PSEL_POOL = 16 + 1024 * 256
ROUTE_MINE = (S_MINE_TOT * (FLUSH + EXC_ENTRY + DONE_C1) // K1 + S_MINE_TOT * H2S_ARR + cd(S_MINE_TOT, K1) * H2S_FULL
              + NODES * (FINAL * F1 + FIN_FIX) + cd(S_MINE_TOT * BLOCKH, K2B) + N_REL * (FINALH * F2 + FIN_FIX)
              + N_HOME * EXPECTH + N_REL * (EXPECT + RDONE) + S_MINE_TOT * (HCHG_NEW + RCV)
              + UNWEIGHTED_VISIT * E_VIS_MINE + HCHG_VER * E_K2_MINE + HCHG_LS1 * E_CAND_MINE
              + E_CAND_MINE * (RESOLVE + 2 * REQ_SERVE + 8 * MSG) + DEF_MINE + EXCH_ROW
              + N_POOLS * (PSEL_POOL + 1024 * (5 * MSG + 32) + MSG + 4 + 8 * (REQ_SERVE + 5 * MSG))
              + ALLN * (MSG + STOP_CHECK) + NODES * MX * (4 * MSG) + 3 * N_HOME + 6 * S_MINE_TOT
              + HR_INIT + NODES * (INIT3_BKT * F1 + 2) + N_REL * ((INIT3_BKT + TINIT_ENTRY) * F2 + 64) + 5 * BARRIER)
E_MINE = NODES * GEN_MINE_NODE + ROUTE_MINE

ZERO_INIT = N_HOME * (2**HB // ZW * ZERO_PASS + 1)
E_PRE = (COLD_TOTAL + E_MINE + NODES * (3 * cd(Q1GD, M_G) + 1) * COIN_SCALAR + NODES * MX * XSETUP
         + NODES * 6 * CH * ENTRY_X + NODES * HEADER + N_REL * (INIT3_BKT * F2 + 1 + TINIT_ENTRY * F2 + 1)
         + NODES * (INIT3_BKT * F1 + 2) + ZERO_INIT + HR_INIT + NODES * XTBLD_ALL + NODES * TBLD)
_SIPT = 240 * RHO
PCFA = NODES * (3 * Q1GD + 2 * NBATCH_MINE + 1) * Fr(1, 2**254)
P_EXC = P_X_ORBIT_ALL_POOLS + P_BAD_X + D_XMINE_ALL + P_MINE_ZERO_ALL + PCFA
BATCH_EXTRA = (MSG - 128) + max(0, 2 * 7 * (40 + 2 * MSG) - 8192)
EXTRA8_BAR = 8 * max(0, cd(BARRIER, ALLN) - 4096) * (NODES + N_REL)
W_DET_MAX = (sum(WGL[g] + WGR[g] + SG[g] * HCHG_NEW + ALLN * STOP_CHECK + EXCH_ROW for g in rn(G_GROUPS))
             + NODES * Q1GD * (B_CAP // CH + 1) * FIX_BATCH * max(FIX_EXTRA.values()) + SG_FULL * JPATH
             + SG_FULL * CAP_H * (UNWEIGHTED_VISIT + HCHG_VER + HCHG_LS1)
             + N_HOME * CAND_CAP * (RESOLVE + 2 * REQ_SERVE)
             + NODES * Q1GD * BATCH_EXTRA + EXTRA8_BAR * Q1GD)
EXC_WORK = ceilf(P_EXC * W_DET_MAX)
PAIR_INVERSE_ALL = NODES * Q1GD * CHK * D_PAIR_INVERSE
PAIR_INVERSE_EXC_WORK = ceilf(PAIR_INVERSE_ALL * W_DET_MAX)

c4 = 2 * K_U_CLUST * SG_FULL + SAMEX_FAC * SG_FULL * (CLUST_BATCH - 1) + SAMEX_FAC * PARTNER_PAIRS + CLUST_BATCH * CLUST_BATCH + 96 * K_U_CLUST * K_U_CLUST
ELIST = RHO * (Fr(SG_FULL) * SG_FULL * KB_M + c4 * KB_M) + _SIPT * KB_M * SG_FULL + 72
EK2 = RHO * (c4 * KB64) + _SIPT * KB64 * SG_FULL + 72
PAIR_SAMEX = Fr(CLUST_BATCH - K_U_CLUST, CLUST_BATCH - 1) * Fr(r, r - 2) * W_ELL * PAIRM + Fr(96 * (K_U_CLUST - 1), CLUST_BATCH - 1) * (KB_BUCKET + Fr(1, r))
SXK2 = PAIR_SAMEX * sum(SAMEX)
SXV = PAIR_SAMEX * UNWEIGHTED_VISIT * sum(SAMEX) + EXC_WORK + PAIR_INVERSE_EXC_WORK
BRACKET = RHO * (Fr(SG_FULL) * SG_FULL * KB(D1 + D2) + c4) + _SIPT * SG_FULL + 72
DEFV = UNWEIGHTED_VISIT * D_BATCH * (CAP_H * SG_FULL + sum(SAMEX))
DEFL_PAIRS = 2 * D_BATCH * BRACKET + D_BATCH * CAP_H * SG_FULL
DEFK2, DEFL, DEFJ = HCHG_VER * DEFL_PAIRS, HCHG_LS1 * DEFL_PAIRS, JPATH * D_BATCH * SG_FULL
DEFWLT = SG_FULL * (RHO * SG_FULL * D_WLT * (2 + D_M) * (UNWEIGHTED_VISIT * KB_C132 + HCHG_VER * KB(D1 + D2 + 64) + (HCHG_LS1 + RESOLVE + 2 * REQ_SERVE) * KB_M)
                    + DWS * KB(132) * (UNWEIGHTED_VISIT + HCHG_VER + HCHG_LS1 + RESOLVE + 2 * REQ_SERVE) + JPATH * D_WLT)
E_FIX_WORK = NODES * Q1GD * E_FIX * max(FIX_EXTRA.values())
E_CAND_WORK = (ELIST + DEFL_PAIRS) * RESOLVE
E_REQ_WORK = (ELIST + DEFL_PAIRS) * (2 * REQ_SERVE)

def _ps(a): return sum(PHI[g] * a[g] for g in rn(G_GROUPS))
WQ = _ps(Q1G)
ROWS = [
    COLD_TOTAL,
    NODES * GEN_MUL,
    (2 + 2 * G_GROUPS) * ALLN * MSG,
    E_MINE,
    NODES * (3 * cd(Q1GD, M_G) + 1) * COIN_SCALAR,
    NODES * MX * XSETUP,
    NODES * 6 * CH * ENTRY_X + NODES * HEADER,
    NODES * XTBLD_ALL,
    NODES * TBLD,
    N_REL * ((INIT3_BKT + TINIT_ENTRY) * F2 + 2) + NODES * (INIT3_BKT * F1 + 2),
    ZERO_INIT,
    HR_INIT,
    _ps(WGL) + NODES * WQ * BATCH_EXTRA + EXTRA8_BAR * WQ,
    E_FIX_WORK,
    _ps(WGJ),
    _ps(WGR),
    _ps(WGG),
    SXV,
    _ps(WGK2) + SXK2 * HCHG_VER + EK2 * HCHG_VER,
    ELIST * HCHG_LS1,
    DEFV + DEFK2 + DEFL + DEFJ + DEFWLT,
    E_CAND_WORK,
    E_REQ_WORK,
    sum(PHI) * EXCH_ROW,
    (NBAR_ALLOW + 5) * BARRIER,
    ALLN * PHASE_CTRL_G,
    OUTPUT,
]
E_RUN = sum(ROWS)
E_TOTAL = GATE_INSNS + P_GATE * E_RUN
W_E = E_TOTAL / C
CLAIM = least_grid_up(W_E, 10**15)
assert P2(CLAIM)[0] >= W_E and P2(CLAIM - Fr(1, 10**15))[1] < W_E

# 8. Exact memory census (disjoint flat 64-bit address spaces, PEAK < 2^138 B)
TiB, PiB, EiB = 1 << 40, 1 << 50, 1 << 60
G_BATCH_META_BYTES = (2 * MX + 2 * cd(Q1GD, M_G) + Q1GD + 1) * 32 + (1 << 30)
HCC_BYTES = (1 << 20) + 8 * 2048 + 4096 + CAND_CAP * 528
G_MEM_REGIONS = [
    ("LC1", 0, 8 * TiB),
    ("YR", 32 * TiB, 33 * TiB),
    ("G34_CTRL_MMIO", 34 * TiB, 38 * TiB),
    ("YQ_YVIEW", 40 * TiB - 64, 40 * TiB + (1 << 36) + (1 << 31)),
    ("EB", 42 * TiB, 42 * TiB + (1 << 37)),
    ("CN", 84 * TiB, 84 * TiB + 64 * Q1GD),
    ("CNM", 2048 * TiB, 2048 * TiB + 64 * NBATCH_MINE),
    ("H1_CONTIG", 64 * PiB, 64 * PiB + F1 * 8192),
    ("F_TBL_5WIN", 512 * PiB, 512 * PiB + 32 * PiB),
    ("FWD_TBL_5WIN", 544 * PiB, 544 * PiB + 16 * PiB),
    ("P1_TBL_5WIN", 1 * EiB, 1 * EiB + 6 * 3 * (1 << 58)),
    ("E_TBL_5WIN", 6 * EiB, 9 * EiB),
    ("G_TBL_5WIN", 10 * EiB, 13 * EiB),
    ("G_BATCH_META", 13 * EiB, 13 * EiB + G_BATCH_META_BYTES),
    ("REQ_CURSOR_MONITOR", 13 * EiB + 512 * PiB - 2048, 13 * EiB + 512 * PiB),
    ("REQ_64B", 13 * EiB + 512 * PiB, 13 * EiB + 512 * PiB + REQ_CAP * 64 + 64),
    ("BT0_52", 0xf800_0000_0000_0000, 0xfc00_0000_0000_0000),
    ("BT1_51x4", 0xfc00_0000_0000_0000, 1 << 64),
]
R_MEM_REGIONS = [
    ("RCS_CTRL", 0x8000, 0x10000 + (1 << 20) + 16384),
    ("R_SLOTS_1024x8K", 0xc000_0000_0000, 0xc000_0000_0000 + 1024 * 8192),
    ("R_SLOT_TABLES", 192 * TiB + (16 << 20), 192 * TiB + (17 << 20)),
    ("R_GENERIC", 192 * TiB + (20 << 20), 192 * TiB + (22 << 20)),
    ("LC2", 32 * TiB, 64 * TiB),
    ("RC2", 96 * TiB, 128 * TiB),
    ("E2P", 128 * TiB, 192 * TiB),
    ("R_SP_TBL", (1 << 48) - (1 << 20), 1 << 48),
    ("DPW_RCB", 256 * TiB, 256 * TiB + 65536 + RELAY_CAP * 64),
    ("H2U_CONTIG", 1 * EiB, 1 * EiB + F2 * 16384),
]
H_MEM_REGIONS = [
    ("H_CTRL_CODE", 0x10000, 0x10000 + (1 << 24)),
    ("H_HCC_CAND", 0x7ffc_0000_0000_00, 0x7ffc_0000_0000_00 + HCC_BYTES),
    ("H_WIN_BUF", 0x7fff_0000_0000_00, 0x7fff_0000_0000_00 + (1 << 24)),
    ("H_SP_HDR", (1 << 55) - (1 << 20), 1 << 55),
    ("H_BULK_CTRL", (1 << 55) - (64 << 20), (1 << 55) - (64 << 20) + 65536),
    ("H_GENERIC", (1 << 55) - (32 << 20), (1 << 55) - (30 << 20)),
    ("R2_38220T", 1 << 55, (1 << 55) + CAP_H * 32),
    ("HH_2p50", 0xff80_0000_0000_0000, 0xff80_0000_0000_0000 + (1 << HB) * 8),
    ("NX_19110T", 0xffa0_0000_0000_0000, 0xffa0_0000_0000_0000 + CAP_H * 16),
]

def _check_map(regs):
    rs = sorted(regs, key=lambda x: x[1])
    for i, (nm, lo, hi) in en(rs):
        if not (0 <= lo < hi <= (1 << 64)): return False, 0
        if i > 0 and rs[i - 1][2] > lo: return False, 0
    return True, sum(hi - lo for _, lo, hi in rs) + PROGRAM_BYTES + (2 << 20)

B_MEM_REGIONS = [
    ("B_DATA", 0, 1 << 62),
    ("B_META", 1 << 62, (1 << 62) + (1 << 31)),
    ("B_INBOX", (1 << 62) + (1 << 32), (1 << 62) + (1 << 32) + (1 << 24)),
    ("B_CTRL", (1 << 62) + (1 << 33), (1 << 62) + (1 << 33) + (1 << 24)),
    ("B_CODE_GENERIC", (1 << 62) + (1 << 34), (1 << 62) + (1 << 34) + (1 << 24)),
]
# Explicit service bindings: programs, the five generic pages, dedicated PE
# stacks, callback/collector state and fixed source tables. Subintervals can
# nest inside reserved regions; distinct live callback/stack intervals are
# checked separately below. The extra _check_map allowance is unused reserve,
# not a substitute for these concrete addresses.
ROLE_MAPS = dict(G=G_MEM_REGIONS, R=R_MEM_REGIONS, H=H_MEM_REGIONS, B=B_MEM_REGIONS)
GENERIC_BASES = dict(G=34*TiB+0x7800000, R=192*TiB+(20<<20),
                     H=(1<<55)-(32<<20), B=(1<<62)+(1<<34))
PROGRAM_BASES = dict(G=34*TiB+(1<<20), R=0x10000, H=0x10000,
                    B=GENERIC_BASES["B"]+(1<<20))
SERVICE_SPANS = []
for role, base in GENERIC_BASES.items():
    for j, name in en(("generic_control", "generic_remote", "generic_ack", "generic_inbox", "generic_handlers")):
        SERVICE_SPANS.append((role, name, base+j*65536, base+(j+1)*65536))
    SERVICE_SPANS.append((role, "16_PE_stacks", base+3*65536+8192, base+3*65536+40960))
    SERVICE_SPANS.append((role, "program", PROGRAM_BASES[role], PROGRAM_BASES[role]+PROGRAM_BYTES))
G34 = 34*TiB
SERVICE_SPANS += [
    ("G", "chain_blocks", G34+0x2000000, G34+0x2000000+2304*256),
    ("G", "YSP", G34+0x4000000, G34+0x4000000+8192),
    ("G", "Y_leaf_state", G34+0x7100000, G34+0x7100000+256),
    ("G", "archive_source", G34+0x7200000, G34+0x7200000+8192),
    ("G", "member_offer", G34+0x7300000, G34+0x7300000+8192),
    ("G", "owner_receiver", G34+0x7500000, G34+0x7500000+8192),
    ("G", "owner_inbox", G34+0x7510000, G34+0x7510000+512),
    ("G", "owner_metadata", G34+0x7520000, G34+0x7520000+256),
    ("G", "command_source", G34+0x7600000, G34+0x7600000+8192),
    ("G", "command_inbox", G34+0x7610000, G34+0x7610000+128),
    ("G", "command_descriptor", G34+0x7620000, G34+0x7620000+128),
    ("G", "driver_state", G34+0x7640000, G34+0x7640000+256),
    ("G", "mining_callback", G34+0x7900000, G34+0x7900000+65536),
    ("G", "collector", G34+0x7a00000, G34+0x7a00000+263680),
    ("R", "member_source_tables", 192*TiB+(16<<20), 192*TiB+(17<<20)),
    ("R", "legacy_source_header", (1<<48)-4096, 1<<48),
    ("H", "bulk_control", (1<<55)-(64<<20), (1<<55)-(64<<20)+65536),
    ("H", "source_header", (1<<55)-4096, 1<<55),
    ("H", "resolve_scratch", GENERIC_BASES["H"]+4*65536+24576, GENERIC_BASES["H"]+4*65536+28672),
    ("H", "resolve_stack", GENERIC_BASES["H"]+5*65536-8192, GENERIC_BASES["H"]+5*65536),
    ("B", "payload_inbox", (1<<62)+(1<<32), (1<<62)+(1<<32)+512),
    ("B", "command_inbox", (1<<62)+(1<<32)+65536, (1<<62)+(1<<32)+65536+128),
    ("B", "archive_control", (1<<62)+(1<<33), (1<<62)+(1<<33)+8192),
]
def _backed(role, start, end):
    return any(lo<=start<end<=hi for _,lo,hi in ROLE_MAPS[role])
check("c_service_backing", all(_backed(role,start,end) for role,_,start,end in SERVICE_SPANS))
# Each program is disjoint from every other explicitly live service span.
check("c_program_placement", all(role!=r or name=="program" or end<=a or b<=start
      for role,name,start,end in SERVICE_SPANS
      for r,a,b in ((r,PROGRAM_BASES[r],PROGRAM_BASES[r]+PROGRAM_BYTES) for r in ROLE_MAPS)))
# These broad selectors are MMIO/readonly routing apertures, not local RAM.
# Per-endpoint latch/register state is counted in FABRIC_BYTES.
SERVICE_ALIASES = [
    ("G", "READY_R", 8*TiB, 16*TiB),
    ("G", "BULK_ACK", 24*TiB, 24*TiB+4096),
    ("G", "YVIEW", 38*TiB, 38*TiB+(1<<35)),
    ("R", "ACK_G", TiB, TiB+8192),
    ("R", "READY_H", 64*TiB, 96*TiB),
]
check("c_service_decode", all(0<=start<end<=1<<64 and
      all(end<=lo or hi<=start for _,lo,hi in ROLE_MAPS[role])
      for role,_,start,end in SERVICE_ALIASES))

okB, B_BYTES = _check_map(B_MEM_REGIONS)
okG, G_BYTES = _check_map(G_MEM_REGIONS)
okR, R_BYTES = _check_map(R_MEM_REGIONS)
okH, H_BYTES = _check_map(H_MEM_REGIONS)
# Reserve 4 * 64 KiB switches + 1 MiB bridge state per physical endpoint (same fabric census as 283cdacf!)
FABRIC_BYTES = (ALLN + N_BANK + 1) * (4 * 65536 + (1 << 20))
PEAK_BYTES = (NODES + 1) * G_BYTES + N_REL * R_BYTES + N_HOME * H_BYTES + N_BANK * B_BYTES + FABRIC_BYTES
check("c40_maps", okG and okR and okH and okB and PEAK_BYTES < (1 << 138))

LIVE_BEST = Fr(123186563709484764, 10**15)
TARGET = Fr(123156563709000000, 10**15)
check("c48_target", CLAIM <= TARGET and len(ROWS) == 27 and sum(ROWS) == E_RUN and all(v >= 0 for v in ROWS))

BUILD_SECONDS = round(time.time() - t0, 3)

def summary():
    return dict(
        all_checks_pass=not FAILED,
        failed=FAILED,
        exact_claim_rational=rat(CLAIM),
        exact_claim_decimal=dec_up(CLAIM, 15),
        live_best_decimal=dec_up(LIVE_BEST, 15),
        gain_vs_live_bits_lower=dec_dn(LIVE_BEST - CLAIM, 15),
        target_margin_bits_lower=dec_dn(TARGET - CLAIM, 15),
        run_log2_upper_6dp=dec_up(least_grid_up(E_RUN, 10**6), 6),
        run_units_log2_upper_6dp=dec_up(least_grid_up(E_RUN / C, 10**6), 6),
        total_log2_upper_6dp=dec_up(least_grid_up(E_TOTAL, 10**6), 6),
        preprocessing_units_log2_upper_4dp=dec_up(least_grid_up(E_PRE / C, 10**4), 4),
        mine_log2_upper_6dp=dec_up(least_grid_up(E_MINE, 10**6), 6),
        per_leaf_upper_4dp=dec_up(E_RUN / SG_FULL, 4),
        success_lower_12dp=dec_dn(LO_RUN, 12),
        gate_upper_15dp=dec_up(P_GATE, 15),
        exact_gate_numerator=str(K_GATE),
        gate_denominator=str(1 << GB),
        peak_memory_bytes=str(PEAK_BYTES),
        peak_memory_log2_upper_6dp=dec_up(least_grid_up(PEAK_BYTES, 10**6), 6),
        request_arena=dict(start=13*EiB+512*PiB, cursor=13*EiB+512*PiB-8,
                           monitor_granule=[13*EiB+512*PiB-2048,13*EiB+512*PiB],
                           end=13*EiB+512*PiB+64*REQ_CAP, capacity=REQ_CAP),
        memory_maps={role:[dict(region=name, begin=str(lo), end=str(hi), bytes=str(hi-lo))
                                for name,lo,hi in regs] for role,regs in ROLE_MAPS.items()},
        per_role_extra_reserved_bytes=str(PROGRAM_BYTES+(2<<20)),
        service_bindings=[dict(role=role, purpose=name, begin=str(lo), end=str(hi))
                          for role,name,lo,hi in SERVICE_SPANS],
        service_aliases=[dict(role=role, purpose=name, begin=str(lo), end=str(hi))
                          for role,name,lo,hi in SERVICE_ALIASES],
        role_bytes=dict(G=str(G_BYTES), R=str(R_BYTES), H=str(H_BYTES), B=str(B_BYTES)),
        archive=dict(groups=G_GROUPS, maximum_group=max(Q1G), minimum_group=min(Q1G),
                     banks_per_owner=1024, bank_slots=BANK_SLOTS, capacity=ARCHIVE_CAP,
                     maximum_active_banks=cd(max(Q1G), BANK_SLOTS), B_population=N_BANK),
        program=dict(reserved=PROGRAM_BYTES, assembled=PROGRAM_SOURCE_BYTES,
                     analytic=PROGRAM_ANALYTIC_BYTES, constructive_total=PROGRAM_SOURCE_BYTES+PROGRAM_ANALYTIC_BYTES),
        pair_inverse=dict(chain_budget=rat(D_PAIR_INVERSE_CHAIN), nominal_chunk_budget=rat(D_PAIR_INVERSE),
                          global_budget=rat(PAIR_INVERSE_ALL), exception_work=str(PAIR_INVERSE_EXC_WORK),
                          batch_budget=rat(D_BATCH), chain_covered_by_nominal_chunks=CHK*CH>=L_CHAIN),
        mining=dict(chains_per_generator=NBATCH_MINE, last_batch_start=NBATCH_MINE-YB,
                    points_per_chain=PTS_MINE_CHAIN, max_index=PTS_PER_GEN_MINE-1,
                    home_chain_threshold=_MINE_K_CH, home_union_exponent_upper=dec_up(_MINE_HOME_EXPONENT,6),
                    ancillary_exclusion_checks=dict(same_A_lt_2m39=P_SAME_A_ALL<Fr(1,1<<39),
                        mining_defect_lt_2m140=D_XMINE_ALL<Fr(1,1<<140),
                        mining_zero_lt_2m140=P_MINE_ZERO_ALL<Fr(1,1<<140),
                        cs_lt_2m82=TFULL["cs"]<Fr(1,1<<82),drop_lt_2m26=D_DROP<Fr(1,1<<26)),
                    pool_mean_lower_6dp=dec_dn(MU_POOL_LO, 6),
                    caps_log2_upper_4dp=dec_up(least_grid_up(P_MINE_CAPS, 10**4), 4)),
        cold=dict(program_bytes=PROGRAM_BYTES, node_init=NODE_INIT,
                  protocol_allowance=str(COLD_PROTOCOL_ALLOWANCE),
                  explicit_bound=str(COLD_EXPLICIT_BOUND), schedule_bound=COLD_SCHEDULE,
                  bank_control_full=str(sum(bank_control(q) for q in Q1G)),
                  bank_control_expected=str(_ps([bank_control(q) for q in Q1G]))),
        role_max_va_hex=dict(G=hex(G_MEM_REGIONS[-1][2]), R=hex(R_MEM_REGIONS[-1][2]), H=hex(H_MEM_REGIONS[-1][2]), B=hex(B_MEM_REGIONS[-1][2])),
        five_window_bounds=dict(ERRH=ERRH, ERRL=ERRL, WH=WH, WL=WL, NWIN=NWIN,
                                D_WJ_log2=dec_up(least_grid_up(D_WJ, 10**4), 4), D_LAB_SYM_log2=dec_up(least_grid_up(D_LAB_SYM, 10**4), 4),
                                D_LABEL_COMBINED_log2=dec_up(least_grid_up(3*D_WJ+D_LAB_SYM, 10**4), 4),
                                D_WLT_log2=dec_up(least_grid_up(D_WLT, 10**4), 4), KB_C132_ratio=round(float(KB_C132 / KB_BUCKET), 4)),
        cap_tails_log2=dict(d_f=dec_up(least_grid_up(TFULL["d_f"], 10**4), 4), d_p2=dec_up(least_grid_up(TFULL["d_p2"], 10**4), 4),
                            d_cl=dec_up(least_grid_up(TFULL["d_cl"], 10**4), 4), d_rq=dec_up(least_grid_up(TFULL["d_rq"], 10**4), 4),
                            d_rl=dec_up(least_grid_up(TFULL["d_rl"], 10**4), 4), glob=dec_up(least_grid_up(TFULL["glob"], 10**4), 4)),
        rows=[dict(row=i + 1, log2_upper_4dp=dec_up(least_grid_up(Fr(v), 10**4), 4),
                   share_pct=round(100.0 * float(Fr(v) / E_RUN), 4)) for i, v in en(ROWS)],
        cert_seconds=BUILD_SECONDS,
    )

if __name__ == "__main__":
    import json
    print(json.dumps(summary(), indent=2))
