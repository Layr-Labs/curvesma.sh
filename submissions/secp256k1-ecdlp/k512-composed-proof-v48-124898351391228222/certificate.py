"""Exact certificate for proof.md: parameters, success bound, work ledger, claim, memory (integers and Fractions only).
Run python3 certificate.py; exit status 0 and a final PASS mean every check holds."""
import math
import os
import sys
from fractions import Fraction

FAILED = []


def check(name, ok, detail=""):
    """Record one named check.

    name: description; ok: truth value; detail: optional text.
    returns: None (prints the line, remembers failures).
    """
    print("%-66s %s%s" % (name, "ok" if ok else "FAIL", ("  " + detail) if detail else ""))
    if not ok:
        FAILED.append(name)


def lg(x):
    """Floating log2 of a positive rational, for display only.

    x: int or Fraction.
    returns: float.
    """
    f = Fraction(x)
    return math.log2(f.numerator) - math.log2(f.denominator)


def cdiv(a, b):
    """Ceiling of a/b for integers, b > 0.

    a, b: integers.
    returns: integer.
    """
    return -((-a) // b)


def exp_interval(x, terms=60):
    """Rational enclosure of exp(x) for |x| <= 1.

    x: Fraction; terms: series length.
    returns: (lo, hi) Fractions with lo <= exp(x) <= hi.
    """
    assert abs(x) <= 1
    s, t = Fraction(1), Fraction(1)
    for k in range(1, terms + 1):
        t = t * x / k
        s += t
    rem = abs(t * x) / (terms + 1) * 3
    return s - rem, s + rem


def ln2_interval(bits=300):
    """Rational enclosure of ln 2 = 2 atanh(1/3).

    bits: target width 2^-bits.
    returns: (lo, hi).
    """
    z, z2 = Fraction(1, 3), Fraction(1, 9)
    tot, pw, k = Fraction(0), z, 0
    while True:
        tot += pw / (2 * k + 1)
        rest = pw * z2 / ((2 * k + 3) * (1 - z2))
        if rest < Fraction(1, 2 ** bits):
            return 2 * tot, 2 * (tot + rest)
        pw *= z2
        k += 1


LN2 = ln2_interval()


def pow2_interval(c):
    """Rational enclosure of 2^c for a rational c.

    c: Fraction.
    returns: (lo, hi).
    """
    n = c.numerator // c.denominator
    f = c - n
    lo, hi = exp_interval(f * LN2[0])[0], exp_interval(f * LN2[1])[1]
    sc = Fraction(2) ** n
    return lo * sc, hi * sc


def up2(x):
    """Least multiple of 1/100 at or above log2(x), exact (x^100 against powers of two), for the tables' display.

    x: positive int or Fraction.
    returns: float.
    """
    f = Fraction(x)
    a, b = f.numerator ** 100, f.denominator ** 100

    def le(k):
        """Whether x^100 <= 2^k. k: integer. returns: bool."""
        return a <= b << k if k >= 0 else a << -k <= b
    k = math.ceil(lg(f) * 100)
    while not le(k):
        k += 1
    while le(k - 1):
        k -= 1
    return k / 100



p = 2 ** 256 - 2 ** 32 - 977
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
C0 = 2 ** 32 + 977
BETA_P = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)   
C = 1000                      
NODES = 2 ** 72               
SIDE = 2 ** 36
Q2 = 2 ** 25                  
D1, D2, D3, D4 = 36, 36, 24, 46
KEYBITS = 110                 
THRESH = Fraction(390000001, 10 ** 9)   


LEAVES = 6                        
FWD_INC, FWD_FIX = 97, 5          
                                  
BWD_INC, BWD_FIX = 1285, 5   
                                  
FWD_CODE_PAIR, BWD_CODE_PAIR = 225, 2741   
                                  
KEYS_CODE_PAIR = 915          
FIX_SITES = dict(SHX=18, FDX=10, F5=2, FS=2, FADDW=1)   
                                                        
FIX_EXTRA = dict(SHX=3, FDX=3, F5=3, FS=4, FADDW=4)     
EMIT = 7                          
EMIT_INC = LEAVES * EMIT          
KEY_INC, KEY_FIX = 140, 5   
                                  
K = 512                           
PB1C, PB2C, D4BASE, BLKC = 1 << 51, 1 << 52, 1 << 49, 1 << 15   
MB1C, MB2C = (1 << 49) - 64, (1 << 49) - 128   
FLUSH = 8                         
                                  
FLUSH_KEYS = FLUSH + 1              
KEY_BLOCK_DROP, H2S_BLOCK_DROP = 8, 9   
                                  
JPATH = 5                         
                                  
INV = 23646                       
TINIT_ENTRY = 3                   
FIN_DEST, FIN_SENT, FIN_HELD, FIN_FIX = 14, 3, 2, 2   
                                  
                                  
H2S_REC, H2S_SLOT = 7, 11                 
SCAN_RANGE = 26                           
H2S_FIX, BKT_FIX, SCAN_FIX = 4, 4, SCAN_RANGE * 2 ** D4 + 2   
                                          
BKT_REC, BKT_DROP, BKT_SLOT = 8, 7, 12    
                                          
PASS4_REC, PASS4_BKT, PASS4_FIX = 3, 11, 2   
SCAN_REC = Fraction(195, 64)              
                                          
LOOP8 = 2                                 
ZERO4_PER4, PREFIX4_PER4, INIT3_BKT = 4, 10, 4   
                                                 
check("key pass: 140 = 2 signs * (6 LDP of the fractions + 28 sums + 3 leaves * (index ADD, CMP, B.HS, CMN, B.HS + "
      "emit 7)), straight-line (no loop TBZ)", KEY_INC == 2 * (6 + 28 + 3 * (5 + EMIT)))
check("j = 1 and j = 2 paths (B.HS taken; CMN, B.HS; ADDS, ADC, MVN of N'; emit; meter ADD; B back) exceed the j = 0 path "
      "(ADD, CMP, B.HS, CMN, B.HS; emit) by 5", (3 + 2 + 3 + EMIT + 2) - (5 + EMIT) == JPATH)
check("emit: AND or BIC (8*d1), LDADD of its LC entry (the slot), d2 and w3 (LSR or EOR, and ADD of the entry pointer "
      "shifted left by 31), two STPs, TBZ of bit 13 for slot k", EMIT == 1 + 1 + 2 + 2 + 1)
DONE = 2   
DLV = 6    
DLVK = 4 * K + 2   
MSG_LINK = DONE + DLV   
BLOCK = FLUSH + DONE + DLVK   
BLOCK_KEYS = FLUSH_KEYS + DONE + DLVK   
FINAL = FIN_DEST + FIN_SENT + DONE + DLVK + MSG_LINK   
L8 = Fraction(LOOP8, 8)   
HOP1, HOP2 = EMIT + Fraction(BLOCK_KEYS, K), H2S_REC + L8 + Fraction(BLOCK, K)
GROUP = BKT_REC + L8 + PASS4_REC + L8 + SCAN_REC
check("a full block: 2060 = flush 8 + in-flight count 2 + dlvk 2050 when sent (2061 at hop 1 with keys' skip B), "
      "more than 8 or 9 when dropped; fin per destination at most 2077 (14 + 3, the final block's 2052 and the final "
      "cursor's 8)",
      BLOCK == 2060 and BLOCK_KEYS == 2061 and BLOCK_KEYS > max(KEY_BLOCK_DROP, H2S_BLOCK_DROP) and FINAL == 2077
      and FINAL > FIN_DEST + FIN_HELD + MSG_LINK)
check("per record: hop 1 7 + 2061/512, hop 2 7.25 + 2060/512 (a block of k = 512 records), grouping 14.546875; "
      "36.845703125 in all", HOP1 + HOP2 + GROUP == Fraction(18865, 512))


GEN_MUL = 677570           
                           
                           
COIN_SCALAR = 64           
BATCH_SETUP = 160          
                           
INV_CALL = 16              
NODE_INIT = 2 ** 20        
                           
EXPAND_NODE = 2 ** 36     
HEADER = 72                
MSG = 40                   
                           
CAND_CAP = 2 ** 28         
RESOLVE = 1557745          
                           
                           
REQ_CAP = 2 ** 29          
REQ_SERVE = 2 ** 10
RELAY_CAP = 2 ** 34        
EXCH = 2 ** 41             
BARRIERS, BARRIER = 32, 2 ** 82   
OUTPUT = 2 ** 20
DELTA_J = 2 ** 42          
DELTA_F = 2 ** 45          
FIX_BATCH = sum(FIX_SITES.values()) * Q2   
CAP_P = 2 ** 22 - 2 ** 12  


P1 = Fraction(6 * (r - 1) + 1, r * r)    
TWO_BOUND = Fraction(2, r - 1)           



SHX_BAD = 2 ** 64 * 2 ** 33     
FOLD_BAD = 2 ** 66 * C0         
CP_BAD = (2 ** 66 + 1) * C0     
SQ_BAD = C0 * C0 + 2 * C0       
SITES = [  
    ("fwd: row shifts of acc_(l-1), l >= 1", 3 * (Q2 - 1), SHX_BAD, 6),
    ("fwd: fold of acc_l, l >= 1", Q2 - 1, FOLD_BAD, 6),
    ("bwd: row shifts of R_l", 3 * Q2, SHX_BAD, 6),
    ("bwd: fold of inv_l = R_l acc_(l-1), l >= 1", Q2 - 1, FOLD_BAD, 6),
    ("bwd: fold of R_(l-1) = R_l D_l, l >= 1", Q2 - 1, FOLD_BAD, 6),
    ("bwd: row shifts of inv_l", 3 * Q2, SHX_BAD, 6),
    ("bwd: folds of u = inv_l (y_Y - y_X) and v", 2 * Q2, FOLD_BAD, 3),
    ("bwd: square folds of u^2 and v^2", 2 * Q2, SQ_BAD, 6),
    ("bwd: folds of E, F, G (tabled rows; b^4, a^4 added to E, F)", 3 * Q2, CP_BAD, 8),
    ("bwd: row shifts of E, F, G", 9 * Q2, SHX_BAD, 8),
    ("bwd: folds of beta u^2 and beta v^2 (tabled rows of beta)", 2 * Q2, FOLD_BAD, 6),
]
SITES_PER_INCIDENCE = dict(SHX=18, FDX=10, F5=2)   
EXACT_AT_L0 = dict(SHX=3, FDX=3)   
D_SITES = sum(Fraction(n * bad * mult, r - 1) for _, n, bad, mult in SITES)


FIX_ROW = C0 * 2 ** 192         
FIX_FOLD = (2 ** 66 + 1) * 2 ** 128 * C0   
                                           
FIX_SQ = 2 ** 193               
FIX_SUB = C0 * 2 ** 192         
FIXUPS = [  
    ("fwd: FS forming D_l", "FS", Q2, FIX_SUB, 6),
    ("fwd: row shifts of acc_(l-1), l >= 1", "SHX", 3 * (Q2 - 1), FIX_ROW, 6),
    ("fwd: fold of acc_l, l >= 1", "FDX", Q2 - 1, FIX_FOLD, 6),
    ("bwd: row shifts of R_l", "SHX", 3 * Q2, FIX_ROW, 6),
    ("bwd: folds of inv_l and R_(l-1), l >= 1", "FDX", 2 * (Q2 - 1), FIX_FOLD, 6),
    ("bwd: FS of y_Y - y_X", "FS", Q2, FIX_SUB, 3),
    ("bwd: FADDW of y_Y + y_X", "FADDW", Q2, FIX_SUB, 3),
    ("bwd: row shifts of inv_l", "SHX", 3 * Q2, FIX_ROW, 6),
    ("bwd: folds of u and v", "FDX", 2 * Q2, FIX_FOLD, 3),
    ("bwd: square folds of u^2 and v^2", "F5", 2 * Q2, FIX_SQ, 6),
    ("bwd: folds of E, F, G", "FDX", 3 * Q2, FIX_FOLD, 8),
    ("bwd: row shifts of E, F, G", "SHX", 9 * Q2, FIX_ROW, 8),
    ("bwd: folds of beta u^2 and beta v^2", "FDX", 2 * Q2, FIX_FOLD, 6),
]
E_FIX = sum(Fraction(n * bad * mult, r - 1) for _, _, n, bad, mult in FIXUPS)
KB = lambda fixed: Fraction(6 * 2 ** (256 - fixed) * (2 ** 17 + 1), (r - 1) * 2 ** 17)  
                              
RHO = Fraction(r, r - 1) ** 4  
FIXED_LEAVES = 2 * LEAVES     
PARTNER_PAIRS = 2 * 10        
                              


def success_terms(q1):
    """Success lower bound pieces for q1 X generators per node.

    q1: integer.
    returns: dict of exact Fractions.
    """
    n0 = LEAVES * q1 * Q2
    S = NODES * n0
    gamma = Fraction(S * S - NODES * n0 * n0, 2)
    mu = gamma * P1
    deg = 2 * LEAVES * (q1 + Q2) * S       
    cs = gamma * deg * (P1 * P1 + Fraction(36, (r - 1) ** 2))
    mu_lo = Fraction((mu.numerator << 200) // mu.denominator, 2 ** 200)
    e_hi = exp_interval(-mu_lo, 40)[1]
    
    d_batch = D_SITES + Q2 * Fraction(6, r - 1) + Fraction(1, r)   
    d_wj = TWO_BOUND * 23 * 2 ** 192   
                                       
    var_j = (Fraction(2, 9) + 2 * d_wj) * n0 + 1
    d_j = var_j / (DELTA_J - n0 // q1 - 1 - 2 * n0 * d_wj) ** 2 + q1 * d_batch
    d_f = q1 * E_FIX / (DELTA_F - FIX_BATCH) + q1 * d_batch   
    d_lab = TWO_BOUND * 1305 * 2 ** 228   
                                          
    
    
    m_p = n0 * KB(D1) + 1        
    v_p = m_p + 1 + 60 * (q1 + Q2) * KB(D1)
    d_p = v_p / Fraction(CAP_P - FIXED_LEAVES - m_p) ** 2
    m_3 = S * KB(D1 + D2 + D3) + 1
    cap_3 = cdiv(m_3.numerator, m_3.denominator) + 2 ** 32
    v_3 = m_3 + 1 + 60 * (q1 + Q2) * KB(D1 + D2 + D3)
    d_3 = v_3 / Fraction(cap_3 - FIXED_LEAVES - m_3) ** 2
    
    
    kb_m = KB(KEYBITS + D1 + D2)   
    lam_int = RHO * (S + FIXED_LEAVES) * KB(D1 + D2 + D3 + D4)
    sip = RHO * Fraction(5, 2) * (1 + 11 * (d_wj + d_lab)) * kb_m   
                                   
    other_h = RHO * (Fraction(S) * S * KB(D1 + D2) * kb_m + (2 * S + PARTNER_PAIRS) * kb_m) + sip * S
    d_cl = other_h / (CAND_CAP - 4)
    other_g = RHO * Fraction(n0) * S * kb_m * 2 + sip * n0
    d_rq = other_g / (REQ_CAP - 8)
    other_v = RHO * (2 + 2 * 8) * (Fraction(S, SIDE) * S * KB(D2) * kb_m + (2 * S + PARTNER_PAIRS) * kb_m) \
        + 18 * sip * S  
    d_rl = other_v / (RELAY_CAP - 16)
    per_leaf = d_batch + d_wj + d_j + d_f + d_lab + 2 * d_p + d_3 + d_cl + d_rq + d_rl
    rows = (mu * (1 + Fraction(1, r)) * (Fraction(1, r) + lam_int + 2 * per_leaf))
    g_def = NODES * q1 * d_batch
    g_rndr = Fraction(0)  
    glob = NODES * (q1 + Q2 + 1) * Fraction(1, 2 ** 254) + g_def + g_rndr
    lo = 1 - e_hi - cs - rows - glob
    return dict(n0=n0, S=S, gamma=gamma, mu=mu, cs=cs, e_hi=e_hi, d_batch=d_batch, d_wj=d_wj, d_j=d_j, d_f=d_f, d_lab=d_lab,
                d_p=d_p, d_3=d_3, lam_int=lam_int, d_cl=d_cl, d_rq=d_rq, d_rl=d_rl, g_def=g_def, rows=rows, glob=glob,
                lo=lo, cap_3=cap_3, m_p=m_p, m_3=m_3)


lo_q, hi_q = 2 ** 27, 2 ** 28
while hi_q - lo_q > 1:
    mid = (lo_q + hi_q) // 2
    if success_terms(mid)["lo"] >= THRESH:
        hi_q = mid
    else:
        lo_q = mid
Q1 = hi_q
T = success_terms(Q1)
N0, S = T["n0"], T["S"]
FA, FB = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16), int("e4437ed6010e88286f547fa90abfe4c3", 16)
BOX = (p >> 100) + 1                      
WIDE = 654 * (p >> 128) + 4               
check("lattice lemma: basis (B, A - B), (A - B, -A) of {(v, u): u = beta v mod p}, |det| = p, parallelogram 2A - B wide "
      "in u and A in v; box factor below 1 + 2^-17, one side's below 1 + 2^-18",
      FA * FA - FA * FB + FB * FB == p and (FB * BETA_P - (FA - FB)) % p == 0 and ((FA - FB) * BETA_P + FA) % p == 0
      and (BOX + WIDE + FA) * (BOX + WIDE + 2 * FA - FB) * 2 ** 200 < p * p * (1 + Fraction(1, 2 ** 17))
      and (BOX + WIDE) * 2 ** 100 < p * (1 + Fraction(1, 2 ** 18)))
check("q1 = %d is the least with success bound >= 0.390000001" % Q1, T["lo"] >= THRESH and success_terms(Q1 - 1)["lo"] < THRESH)
check("leaves S = 2^72*6*q1*q2 = 2^%.6f" % lg(S), S == NODES * 6 * Q1 * Q2)
check("mu = |Gamma|*P1 = %.9f (cross-node pairs)" % float(T["mu"]), T["mu"] > Fraction(494, 1000))
check("Chen-Stein b1+b2 < 2^-95", T["cs"] < Fraction(1, 2 ** 95), "2^%.2f" % lg(T["cs"]))
check("intruder rate per colliding pair < 2^-12.7", T["lam_int"] < pow2_interval(Fraction(-127, 10))[0], "2^%.3f" % lg(T["lam_int"]))


def isum(*rs):
    """Interval sum. rs: (lo, hi) pairs. returns: (lo, hi)."""
    return sum(x[0] for x in rs), sum(x[1] for x in rs)


TFP_TERMS = ((0, 1), (0, 7), (0, 56))   
RAW = isum(*TFP_TERMS)                  
COR = isum(RAW, (-1, 0))                
U3, UC3 = (RAW, RAW, (-2 - 2 * RAW[1], -2 * RAW[0])), (COR, COR, (-2 - 2 * COR[1], -2 * COR[0]))   
check("3.3's fraction ranges from the terms e 0..1, c_low 0..7, d 0..56: T - estimate in [0, 64], corrected [-1, 64]; "
      "U_2 in [-130, 0], corrected [-130, 2]", RAW == (0, 64) and COR == (-1, 64) and U3[2] == (-130, 0)
      and UC3[2] == (-130, 2))
WR = [[isum(UC3[m], UC3[(m + t) % 3], U3[(m + 2 * t) % 3], (0, 2)) for t in range(3)] for m in (0, 1)]   
NR = [(-WR[0][t][1] - WR[1][t][1] - 1, -WR[0][t][0] - WR[1][t][0] + 1) for t in range(3)]   
LABEL_ERR = [(min(x[0] for x in rs), max(x[1] for x in rs)) for rs in (WR[0], WR[1], NR)]
check("3.3's error table from the fraction ranges: W0 per t [-2, 194], [-132, 130], [-131, 132]; W1 [-2, 194], "
      "[-131, 132], [-132, 130]; N' [-389, 5], [-263, 264], [-263, 264]; over t [-132, 194], [-132, 194], [-389, 264]: "
      "D' 1305 residues (d_lab's), widening 654 (WIDE's)",
      WR == [[(-2, 194), (-132, 130), (-131, 132)], [(-2, 194), (-131, 132), (-132, 130)]]
      and NR == [(-389, 5), (-263, 264), (-263, 264)] and LABEL_ERR == [(-132, 194), (-132, 194), (-389, 264)]
      and sum(hi - lo for lo, hi in LABEL_ERR) == 1305 and max(hi - lo for lo, hi in LABEL_ERR) + 1 == 654
      and WIDE == 654 * (p >> 128) + 4)
check("wrong class index rate per leaf < 2^-58.47", T["d_wj"] < pow2_interval(Fraction(-5847, 100))[0], "2^%.4f" % lg(T["d_wj"]))
for name, n, bad, mult in SITES:
    k = -(-n // Q2)
    print("  site %-52s %-6s bad <= 2^%.2f, x%d/(r-1), each <= 2^%.2f" % (name, ("%dq2" % k if k > 1 else "q2") + (
        "%+d" % (n - k * Q2) if n != k * Q2 else ""), up2(bad), mult, up2(Fraction(bad * mult, r - 1))))
check("site table: 18 row shifts (3 by SHXT), 10 folds, 2 square folds per incidence, less 3 + 3 exact at l = 0",
      sum(n for _, n, b, _ in SITES if b == SHX_BAD) == SITES_PER_INCIDENCE["SHX"] * Q2 - EXACT_AT_L0["SHX"]
      and sum(n for _, n, b, _ in SITES if b in (FOLD_BAD, CP_BAD)) == SITES_PER_INCIDENCE["FDX"] * Q2
      - EXACT_AT_L0["FDX"]
      and sum(n for _, n, b, _ in SITES if b == SQ_BAD) == SITES_PER_INCIDENCE["F5"] * Q2 and len(SITES) == 11)
check("bad sets: rows 2^97, folds below 2^98.02, square folds below 2^64.01",
      SHX_BAD == 2 ** 97 and all(b < pow2_interval(Fraction(e, 100))[0] for b, e in ((FOLD_BAD, 9802), (CP_BAD, 9802),
                                                                                     (SQ_BAD, 6401))))
for name, _, n, bad, mult in FIXUPS:
    k = -(-n // Q2)
    print("  fixup %-51s %-6s words <= 2^%.2f, x%d/(r-1), each <= 2^%.2f" % (name, ("%dq2" % k if k > 1 else "q2") + (
        "%+d" % (n - k * Q2) if n != k * Q2 else ""), up2(bad), mult, up2(Fraction(bad * mult, r - 1))))
FIX_RUNS = {m: sum(n for _, ms, n, _, _ in FIXUPS if ms == m) for m in FIX_SITES}
check("fixup table: per incidence 18 row shifts, 10 folds, 2 square folds, 2 FS and 1 FADDW, less 3 + 3 exact at l = 0",
      FIX_RUNS == {m: FIX_SITES[m] * Q2 - EXACT_AT_L0.get(m, 0) for m in FIX_SITES}
      and all(bad == dict(SHX=FIX_ROW, FDX=FIX_FOLD, F5=FIX_SQ, FS=FIX_SUB, FADDW=FIX_SUB)[m] for _, m, _, bad, _ in FIXUPS)
      and all(FIX_SITES[m] == SITES_PER_INCIDENCE[m] for m in ("SHX", "FDX", "F5")))
check("a call's fixup word (C0 per fixup, at most 33 q2 fixups) stays below 2^63", FIX_BATCH * C0 < 2 ** 63,
      "2^%.3f" % lg(FIX_BATCH * C0))
check("mean fixups per node (clean batches) q1*E_FIX = 2^%.3f < 2^28.68; cap 2^45 > 2^16 times it plus 33 q2" % lg(Q1 * E_FIX),
      Q1 * E_FIX < pow2_interval(Fraction(2868, 100))[0] and DELTA_F > 2 ** 16 * Q1 * E_FIX + FIX_BATCH)
check("fixup meter rate (Markov) < 2^-16.32", T["d_f"] < pow2_interval(Fraction(-1632, 100))[0], "2^%.3f" % lg(T["d_f"]))
check("label-error rate per leaf < 2^-16.65", T["d_lab"] < pow2_interval(Fraction(-1665, 100))[0],
      "2^%.4f" % lg(T["d_lab"]))
check("batch defect rate (site table, D = 0, X = O) < 2^-126.05", T["d_batch"] < pow2_interval(Fraction(-12605, 100))[0],
      "2^%.3f" % lg(T["d_batch"]))
check("j = 1, 2 meter rate < 2^-30", T["d_j"] < Fraction(1, 2 ** 30), "2^%.2f" % lg(T["d_j"]))
check("pair-slot overflow rate, each hop, < 2^-19.9", T["d_p"] < pow2_interval(Fraction(-199, 10))[0],
      "2^%.3f" % lg(T["d_p"]))
check("d3 bucket overflow rate < 2^-30", T["d_3"] < Fraction(1, 2 ** 30), "2^%.2f" % lg(T["d_3"]))
check("candidate cap rate < 2^-23.4, request cap rate < 2^-26.0",
      T["d_cl"] < pow2_interval(Fraction(-234, 10))[0] and T["d_rq"] < pow2_interval(Fraction(-260, 10))[0],
      "2^%.3f 2^%.3f" % (lg(T["d_cl"]), lg(T["d_rq"])))
check("relay cap rate (a reply counts 8 records) < 2^-25.2", T["d_rl"] < pow2_interval(Fraction(-252, 10))[0],
      "2^%.3f" % lg(T["d_rl"]))
check("whole-run row: a computational defect in any batch < 2^-26.5", T["g_def"] < pow2_interval(Fraction(-265, 10))[0],
      "2^%.3f" % lg(T["g_def"]))
check("whole-run rows (defective batch anywhere, coin rejection twice) < 2^-26.5",
      T["glob"] < pow2_interval(Fraction(-265, 10))[0], "2^%.3f" % lg(T["glob"]))
check("success >= 1 - e^-mu - (b1+b2) - rows - global >= 0.390000001", T["lo"] >= THRESH,
      "%.10f" % (math.floor(T["lo"] * 10 ** 10) / 10 ** 10))


CAP_J = cdiv(2 * N0, 3) + DELTA_J
GENS = Q1 + Q2 + 1
rows = []


def row(name, value, pre=False):
    """Add one ledger row.

    name: label; value: exact integer instructions; pre: True if the row is preprocessing (independent of Q).
    returns: None.
    """
    assert isinstance(value, int) and value >= 0
    rows.append((name, value, pre))


row("program delivery, node and port state: 2^72 * 2^20", NODES * NODE_INIT, True)
row("coins: (q1+q2+1) scalars * 64 per node", NODES * GENS * COIN_SCALAR, True)
row("Y generators [c]P with their table entries: q2 * 677570 per node", NODES * Q2 * GEN_MUL, True)
row("R_g = [b]Q and X = R_g + [a]P with their header values: (q1+1) * 677570 per node", NODES * (Q1 + 1) * GEN_MUL)
row("rows of beta: 72 per node", NODES * HEADER, True)
row("LC tables (init3) and RC tables (tinit): 2 * (4 * 2^36 + 1) + 2 * (3 * 2^36 + 1) per node",
    NODES * 2 * (INIT3_BKT * SIDE + 1 + TINIT_ENTRY * SIDE + 1), True)
row("Q broadcast: 2 messages to every node over two hops at 40", 2 * (SIDE + NODES) * MSG)
row("d4 counters zeroed: 2^46/4 * 4 + 1 per node", NODES * (2 ** D4 // 4 * ZERO4_PER4 + 1), True)
row("d3 cursor table: 4 * 2^24 + 1 per node", NODES * (INIT3_BKT * 2 ** D3 + 1), True)
row("forward: q1 * (97 q2 + 5)", NODES * Q1 * (FWD_INC * Q2 + FWD_FIX))
row("inversion: q1 * (23,646 + 16)", NODES * Q1 * (INV + INV_CALL))
row("backward: q1 * (1285 q2 + 5)", int(NODES * Q1 * (BWD_INC * Q2 + BWD_FIX)))
row("straight-line code expansion: 2^72 * 2^36", NODES * EXPAND_NODE, True)
row("arithmetic fixups: 2^45 * 4 per node", NODES * DELTA_F * max(FIX_EXTRA.values()))
row("key pass without emits: q1 * (98 q2 + 5)", int(NODES * Q1 * ((KEY_INC - EMIT_INC) * Q2 + KEY_FIX)))
row("j = 1 and j = 2 key paths: cap_j * 5 per node", NODES * CAP_J * JPATH)
row("batch set-up: q1 * 160", NODES * Q1 * BATCH_SETUP)
row("hop 1, the emits and full blocks (flush, in-flight count, dlvk): S * (7 + 2061/512)", S * (EMIT * K + BLOCK_KEYS) // K)
row("hop 2, h2s and full blocks: S * (7.25 + 2060/512) + 11 * 2^36 + 4 per node",
    S * (H2S_REC * K + LOOP8 * K // 8 + BLOCK) // K + NODES * (H2S_SLOT * SIDE + H2S_FIX))
row("final blocks and cursors, in-flight counts and deliveries, both hops: 2 * (2077 * 2^36 + 2) per node",
    NODES * 2 * (FINAL * SIDE + FIN_FIX))
row("bucket pass bkt: S * 8.25 + 12 * 2^36 + 4 per node", S * BKT_REC + S * LOOP8 // 8 + NODES * (BKT_SLOT * SIDE + BKT_FIX))
row("prefix sums of d4: 2^46/4 * 10 + 1 per node", NODES * (2 ** D4 // 4 * PREFIX4_PER4 + 1))
row("local pass pass4: S * 3.25 + 11 * 2^24 + 2 per node",
    S * PASS4_REC + S * LOOP8 // 8 + NODES * (PASS4_BKT * 2 ** D3 + PASS4_FIX))
row("scan_ranges: S * 195/64 + 26 * 2^46 + 2 + 8 * (2^28 + 1) per node", S * SCAN_REC.numerator // SCAN_REC.denominator
    + NODES * (SCAN_FIX + 8 * (CAND_CAP + 1)))
row("candidates: 2^28 * 1557745 per node", NODES * CAND_CAP * RESOLVE)
row("requests served: 2^29 * 2^10 per node", NODES * REQ_CAP * REQ_SERVE)
row("phase-4 count exchanges: 4 * 2^41 per node", NODES * 4 * EXCH)
row("barriers: 32 * 2^82", BARRIERS * BARRIER)
row("phase control: kernel calls between barriers, 2^10 per node", NODES * 2 ** 10)
row("output", OUTPUT)
R_TOTAL = sum(v for _, v, _ in rows)
R_PRE = sum(v for _, v, pre in rows if pre)
W = Fraction(R_TOTAL, C)
print()
for name, v, pre in rows:
    print("  %-58s %s 2^%.4f" % (name, "P" if pre else " ", lg(v)))
print("  %-58s   2^%.6f instructions" % ("total", lg(R_TOTAL)))
print()


def claim_of(w):
    """Smallest multiple of 10^-12 whose power of two is at least w, verified by interval bounds.

    w: Fraction > 1.
    returns: Fraction claim.
    """
    c = Fraction(math.ceil(lg(w) * 10 ** 15), 10 ** 15)
    for _ in range(400):
        if pow2_interval(c)[0] >= w and pow2_interval(c - Fraction(1, 10 ** 15))[1] < w:
            return c
        c += Fraction(1, 10 ** 15) if pow2_interval(c)[0] < w else -Fraction(1, 10 ** 15)
    raise AssertionError("claim not certified")


CLAIM = claim_of(W)
check("claim = least 15-decimal c with 2^c >= R/1000", pow2_interval(CLAIM)[0] >= W
      and pow2_interval(CLAIM - Fraction(1, 10 ** 15))[1] < W, "%d/10^15" % (CLAIM * 10 ** 15))
PER_LEAF = Fraction(R_TOTAL, S)
check("work per leaf = %.7f instructions" % float(PER_LEAF), 286 < PER_LEAF < 290)


CAP_3 = T["cap_3"]
SLOTS = SIDE * CAP_P * 32     
FIXED = [("code incl. straight-line fwd/bwd/keys (.rept Q2/2), constants, stack", 2 ** 38),
         ("barrier words: 5 granules of at most 2 KiB", 5 * 2048),
         ("header H and the rows of beta below it", 1024), ("X table: H's rows and values, -T(a), -T(beta a), y^2 per X", Q1 * 512),
         ("Y table: y^2, y, b, b^3, 3b^2, b^4", Q2 * 192), ("scalars a, b, c", GENS * 32), ("batch entries", Q2 * 256),
         ("the ports' message blocks", 2 * 64), ("LC tables at 0 and 2^39", 2 * SIDE * 8),
         ("RC tables at 2^43 and 2^43 + 2^39", 2 * SIDE * 8), ("record blocks at 2^51 and 2^52",
          2 * SIDE * 2 ** 15),
         ("final-cursor tables E1 and E2", 2 * SIDE * 32), ("d3 cursor table at 2^40", 2 ** D3 * 8),
         ("d4 counters at 2^49", 2 ** D4 * 8), ("candidate list", CAND_CAP * 16),
         ("request and reply buffers", REQ_CAP * 64 + 8 * 2 * CAND_CAP * 64), ("relay buffer", RELAY_CAP * 64),
         ("phase-4 count, offset and cursor tables", 3 * SIDE * 32)]
FIXED_BYTES = sum(b for _, b in FIXED)
BUCKET_STRIDE = 2 ** 38      
PHASE3 = [("d3 buckets: 2^24 strides of 2^38 bytes, cap_3 entries (address, d4 counter) at the end of each",
           2 ** D3 * BUCKET_STRIDE),
          ("pointer array", 2 ** D3 * CAP_3 * 8)]
PHASE3_BYTES = sum(b for _, b in PHASE3)
for name, b in FIXED + [("R1 at 2^53: hop-1 pair slots (phases 1-2)", SLOTS),
                        ("R2 at 2^53 + |R1|: hop-2 pair slots (phases 2-4)", SLOTS)] + PHASE3:
    print("  %-58s 2^%.4f bytes" % (name, lg(b)))
check("fixed regions lie below 2^53", FIXED_BYTES <= 2 ** 53, "2^%.4f" % lg(FIXED_BYTES))
ALIGNED = 4 * SIDE * 8 + 2 ** D3 * 8 + 2 ** D4 * 8 + 2 * SIDE * 2 ** 15
check("aligned tables: LC [0, 2^39), [2^39, 2^40), d3 [2^40, 2^40 + 2^27), RC [2^43, 2^44), d4 [2^49, 2^50), "
      "column blocks [2^51, 2^52), row blocks [2^52, 2^52 + 2^51); message blocks at 2^49 - 128 and 2^49 - 64 "
      "(48 bytes each, disjoint from d4 and the record blocks); the rest fits between 2^44 and 2^45",
      SIDE * 8 == 2 ** 39 and 2 ** D3 * 8 == 2 ** 27 and 2 ** D4 * 8 == 2 ** 49 and SIDE * 2 ** 15 == 2 ** 51
      and D4BASE == 1 << 49 and D4BASE + 2 ** D4 * 8 <= PB1C and PB1C == 1 << 51 and PB2C == 1 << 52
      and MB1C + 48 <= D4BASE and MB2C + 48 <= MB1C and (1 << 48) <= MB2C
      and PB2C + (1 << 51) <= 1 << 53 and FIXED_BYTES - ALIGNED <= (1 << 45) - (1 << 44),
      "2^%.4f" % lg(FIXED_BYTES - ALIGNED))
check("record block: header and k slots within 2^15; bit 14 set exactly at slot k (TBZ #14); bases/stride keep "
      "low 14 bits zero; cap_p a multiple of k (keep lemma)",
      32 + 32 * K <= 2 ** 15 and all((((32 * n) >> 14) & 1 == 1) == (n == K) for n in range(1, K + 1))
      and PB1C & 0x3FFF == 0 and PB2C & 0x3FFF == 0 and BLKC & 0x3FFF == 0 and CAP_P % K == 0
      and S % K == 0)
check("batch entries fill the lower 2^33 bytes of their 2^34 block (all three kernels straight-line .rept Q2/2)",
      Q2 * 256 == 2 ** 33)
check("record fields: k < 2^28 in w2 bits 36-63; the sign's code 2 or 5 in w3 bits 36-38 and l < 2^25 (l - 1 mod 2^25 "
      "for s = -1) in bits 39-63", Q1 <= 2 ** 28 and Q2 <= 2 ** 25)
check("R1 and R2 end below 2^64", 2 ** 53 + 2 * SLOTS < 2 ** 64, "2^%.6f" % lg(2 ** 53 + 2 * SLOTS))
check("phase 3: buckets and pointer array fit in R1's range (R1 is dead after hop 2)", PHASE3_BYTES <= SLOTS,
      "2^%.4f" % lg(PHASE3_BYTES))
check("d3 buckets: cap_3 entries of 16 bytes fit below a stride of 2^38, and R1 is a multiple of the stride",
      16 * CAP_3 < BUCKET_STRIDE and 2 ** 53 % BUCKET_STRIDE == 0, "16*cap_3 = 2^%.4f" % lg(16 * CAP_3))
check("hop-1 cursors cannot wrap: the last slot's base plus 32 per leaf of a node stays below 2^64 (hop-2 and bucket "
      "cursors stop at their end)", 2 ** 53 + 32 * CAP_P * (SIDE - 1) + 32 * N0 < 2 ** 64,
      "2^%.4f" % lg(2 ** 53 + 32 * CAP_P * (SIDE - 1) + 32 * N0))
NODE_BYTES = FIXED_BYTES + 2 * SLOTS   
PEAK = NODES * NODE_BYTES
check("peak memory 2^72 * (fixed + R1 + R2) = 2^%.4f bytes" % lg(PEAK), PEAK < 2 ** 136)
check("caps exceed crude means: cap_p > m_p + 2^20, cap_3 > m_3 + 2^31",
      CAP_P > T["m_p"] + 2 ** 20 and CAP_3 > T["m_3"] + 2 ** 31)
check("preprocessing (rows marked P) = 2^%.7f units, included in the total" % lg(Fraction(R_PRE, C)), R_PRE < R_TOTAL)


DECLARED = dict(time_log2=Fraction(124898351391228222, 10 ** 15), memory_log2_bytes=Fraction(1359990, 10 ** 4),
                preprocessing_log2=Fraction(1084542, 10 ** 4), success=Fraction(390000001, 10 ** 9),
                per_leaf=Fraction(2869510, 10 ** 4))


def least_at_or_above(d, x, scale):
    """Whether 2^d >= x and 2^(d - 1/scale) < x, both certified by the interval bounds of pow2_interval.

    d: Fraction; x: positive int or Fraction; scale: denominator of the decimal grid.
    returns: bool.
    """
    return pow2_interval(d)[0] >= x and pow2_interval(d - Fraction(1, scale))[1] < x


check("declared time_log2 %d/10^15: the claim, least 15-decimal value with 2^c >= R/1000" % (DECLARED["time_log2"] * 10 ** 15),
      DECLARED["time_log2"] == CLAIM and least_at_or_above(DECLARED["time_log2"], W, 10 ** 15))
check("declared memory 2^135.9989 bytes: least 4-decimal bound >= peak (2^%.7f)" % lg(PEAK),
      least_at_or_above(DECLARED["memory_log2_bytes"], PEAK, 10 ** 4))
check("declared preprocessing 2^%s units: least 4-decimal bound >= rows marked P / 1000" % format(float(DECLARED["preprocessing_log2"]), ".4f"),
      least_at_or_above(DECLARED["preprocessing_log2"], Fraction(R_PRE, C), 10 ** 4))
check("declared success 0.390000001 <= the bound", DECLARED["success"] <= T["lo"])
check("declared work per leaf 286.9510: least 4-decimal bound >= R/S (%.7f)" % float(PER_LEAF),
      DECLARED["per_leaf"] >= PER_LEAF > DECLARED["per_leaf"] - Fraction(1, 10 ** 4))
ANALYTIC = sum(v for n, v, _ in rows if n.split(":")[0] in ("program delivery, node and port state",
               "straight-line code expansion", "coins",
               "Y generators [c]P with their table entries", "R_g = [b]Q and X = R_g + [a]P with their header values",
               "rows of beta", "batch set-up",
               "Q broadcast", "candidates", "requests served", "phase-4 count exchanges", "barriers", "phase control",
               "output",
               "inversion")) - NODES * Q1 * INV



MULX, SQRX, ADDX = 115 + 4, 82 + 3, 9 + 2 * 9
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX

GEN_ITEM = 256 * (DBL + MADD_) + 2 * (INV + 3 * MULX + 2 * SQRX + 10 * ADDX) + 64
GEN_TEST = 256 * 2 * (9 + 4) + 64   
                                    
GEN_CTRL = 1600   
NEGT = 30   
            
            
RED = 9     
Y_ENTRY = 3 * SQRX + 2 * MULX + 2 * ADDX + 2 * NEGT + 7 * RED + 14   
                                            
X_ENTRY = MADD_ + 3 * SQRX + 11 * MULX + 2 * ADDX + 2 * NEGT + 17 * RED + 32   
                                            
                                            
                                            
check("generator (exact arithmetic): doubling 933, mixed addition 1,477; with the exceptional-case tests, control, "
      "the Y entry and the X header values at most %d, the allowance" % (GEN_ITEM + GEN_TEST + GEN_CTRL + max(Y_ENTRY,
      X_ENTRY)), DBL == 933 and MADD_ == 1477 and GEN_ITEM + GEN_TEST + GEN_CTRL + max(Y_ENTRY, X_ENTRY) <= GEN_MUL,
      "%.4fx" % (math.floor(10 ** 4 * GEN_MUL / (GEN_ITEM + GEN_TEST + GEN_CTRL + max(Y_ENTRY, X_ENTRY))) / 10 ** 4))



BARRIER_ITEM = 5 * 8 + 5 * (6 + MSG_LINK) + 1 + 20









BATCH_ITEM = (2 + 4) + 60 + (1 + 4 + 4 + 4 + 2 + 1) + (2 * 4 + 1) + (4 + 1 + 2 + 1 + 1 + 1) + 2 * (4 + 1) \
    + 2 * (2 + 1 + 3) + 2 * (1 + 4 + 1 + 1) + (2 + 2 + 1 + 1) + 2
check("batch set-up: at most %d instructions, within its allowance of 160" % BATCH_ITEM, BATCH_ITEM <= BATCH_SETUP)
check("barrier: at most 131 instructions per node, in-flight counts and deliveries included (2^10 per node is > 7x)",
      BARRIER_ITEM == 131 and BARRIER > 7 * NODES * BARRIER_ITEM)




EXCH_ITEM = 12 * SIDE + 4 * (14 + MSG_LINK) * (8 * REQ_CAP)
check("phase-4 count exchange: < 2^40.5 per node, within 2^41 (at most 8 * 2^29 messages)",
      EXCH_ITEM < pow2_interval(Fraction(405, 10))[0] and EXCH_ITEM < EXCH, "2^%.3f" % lg(EXCH_ITEM))


check("a broadcast message costs at most 16, in-flight count and delivery included, within its price of 40",
      2 + 3 + 3 + MSG_LINK <= MSG)




CAND_LEAVES = INV + 5 * SQRX + 18 * MULX + 20 * ADDX
CAND_ITEM = 4 * 8 * 2 * MSG + 2 * CAND_LEAVES + 9 * 600 + 384 * 200 + 2 * 256 * (DBL + MADD_) // 2 + INV
check("candidate: messages, both records' three candidate leaves, nine unit tests (600 each), solve mod r, [x]P = Q = "
      "%d; charged above twice that" % CAND_ITEM, RESOLVE == 2 * CAND_ITEM + 1,
      "%.3fx" % (math.floor(10 ** 3 * RESOLVE / CAND_ITEM) / 10 ** 3))
PROGRAM_BYTES = 2 ** 16
check("program delivery: 2^16 bytes = 2^11 messages over 2 hops at 40 = 2^17.3 < node set-up 2^20",
      PROGRAM_BYTES // 32 * 2 * MSG < NODE_INIT)
HERE = os.path.dirname(os.path.abspath(__file__))
SRC_BYTES = sum(os.path.getsize(os.path.join(HERE, "kernels", f)) for f in ("leaf.s", "route.s", "secp_fe.s"))
check("compact template: measured .s sources + 16 KiB (expansion, generator, resolver, phase control) <= 2^16",
      SRC_BYTES + 16384 <= 2 ** 16, "%d bytes" % (SRC_BYTES + 16384))
CODE_BYTES = 4 * ((FWD_CODE_PAIR + BWD_CODE_PAIR + KEYS_CODE_PAIR) * (Q2 // 2) + 10)
check("straight-line stream: 4*((225 + 2741 + 915)*2^24 + 10) bytes fits the 2^38 code reservation",
      CODE_BYTES <= 2 ** 38, "2^%.4f" % lg(CODE_BYTES))
EXPAND_ITEM = 3 * ((CODE_BYTES + 15) // 16) + 2 * ((CODE_BYTES + 63) // 64) + 8
check("expansion itemized per node (3 per 16-byte unit + DC CVAU and IC IVAU per line + 8) within 2^36",
      EXPAND_ITEM <= EXPAND_NODE, "%d" % EXPAND_ITEM)
check("straight-line .rept is exact: Q2/2 pairs, no residual trip", Q2 % 2 == 0)
check("analytic allowances are below 2^-13.7 of the total", ANALYTIC * pow2_interval(Fraction(137, 10))[1] < R_TOTAL,
      "2^%.3f" % lg(Fraction(ANALYTIC, R_TOTAL)))

print()
print("q1 %d  q2 2^25  n0 %d  S 2^%.6f  mu %.9f  success >= %.10f" % (Q1, N0, lg(S), float(T["mu"]),
                                                                    math.floor(T["lo"] * 10 ** 10) / 10 ** 10))
print("R_total %d instructions  W 2^%.10f  claim %d/10^15  per leaf %.4f" % (R_TOTAL, lg(W), CLAIM * 10 ** 15,
                                                                      float(PER_LEAF)))
print("declared: memory_log2_bytes %s (exact %.7f)  preprocessing_log2 %s (exact %.7f)" % (
    format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK), format(float(DECLARED["preprocessing_log2"]), ".4f"),
    lg(Fraction(R_PRE, C))))
print("FAIL: " + ", ".join(FAILED) if FAILED else "PASS")
if __name__ == "__main__":
    sys.exit(1 if FAILED else 0)
