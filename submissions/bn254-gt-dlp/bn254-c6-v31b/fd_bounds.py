#!/usr/bin/env python3
"""fd_bounds.py -- input-uniform bounds for the 12-D Q38/Q32 Fincke-Pohst tree of v26 Lattice Band 2 (stdlib, exact
Fraction arithmetic with outward-rounded square roots, ~0.2 s). Evaluates the constants of KERNELS.md section 1.3
Theorem A (completeness: every point with exact Q12 <= T is emitted; pad_1 slack lemma) and Theorem B (every emitted
leaf has Q12 <= (1 + eps_or) T), using only: min ||b_k*||^2 >= (13/20)^11 * 753,729,000 = 2^22.6531, T12 <= 2^46.803,
LLL(99/100) (|mu_jk| <= 1/2, B_{k+1} >= 0.74 B_k), and beta_k = B_k/T12 < 1 at every level that runs the fixed-point
body (levels with beta_k >= 1 run the exact-integer body; Lemma E: only k >= 3; KERNELS.md section 1.3). The radius bound uses the 7-op LUT's ratio Lambda (bucket-exhaustive certificate below); the
5-op LUT4 actually used is checked in v20_lat_b2.py to be outward and never above the 7-op LUT in any bucket, so
Lambda bounds it too."""
import math, random, sys, time
from fractions import Fraction as Fr
Q32 = 1 << 32
Q38 = 1 << 38
M64 = (1 << 64) - 1
BMIN = Fr(753729000) * Fr(13, 20) ** 11          # min_k ||b_k*||^2  (= 2^22.6531)
TMAX = int(2.0 ** 46.803) + 2                    # T12 upper bound (2^46.803), rounded up
KAPPA = Fr(99, 100) - Fr(1, 4)                   # Lovasz: B_{k+1} >= KAPPA * B_k  (|mu| <= 1/2)
BETA_LO = BMIN / TMAX                            # min beta_k = B_k / T   (= 2^-24.150)
LUT = [0] * 512 + [math.ceil(math.sqrt(i + 1) * (2.0 ** 52.5)) + (1 << 31) for i in range(512, 2048)]
MD0 = 2 * 1000000000039
def sqrt_up(q, p=60):
    """rational upper bound of sqrt(q), error <= 2^-p."""
    q = Fr(q)
    n = q.numerator << (2 * p)
    v = -((-n) // q.denominator)
    return Fr(math.isqrt(v) + 1, 1 << p)
def sqrt_dn(q, p=60):
    q = Fr(q)
    return Fr(math.isqrt((q.numerator << (2 * p)) // q.denominator), 1 << p)
def isqrt_lut(r2):
    if not (1 <= r2 < (1 << 40)):
        raise OverflowError("r2 outside [1, 2^40): %d" % r2)
    s = (64 - r2.bit_length()) & ~1
    return LUT[(r2 << s) >> 53] >> (s >> 1)
def lut_certificate():
    """For every (s, idx) bucket of r2 in [1, 2^40): (i) rad^2 > r2_max * 2^52 (strict outward on the whole
    bucket); (ii) Lambda = max rad^2 / (r2_min * 2^52) (exact squared overshoot factor)."""
    errs = 0
    lam = Fr(0)
    nb = 0
    for s in range(0, 63, 2):
        lo_s, hi_s = max(1, 1 << max(0, 62 - s)), min(1 << 40, 1 << (64 - s))
        if s == 62:
            lo_s = 1
        if lo_s >= hi_s:
            continue
        for idx in range(512, 2048):
            a = -((-(idx << 53)) >> s) if s <= 53 else -((-(idx << 53)) // (1 << s))
            b = -((-((idx + 1) << 53)) // (1 << s)) - 1
            r_min, r_max = max(a, lo_s), min(b, hi_s - 1)
            if r_min > r_max:
                continue
            for r in (r_min, r_max):
                ss = (64 - r.bit_length()) & ~1
                if ss != s or ((r << ss) >> 53) != idx:
                    errs += 1
            nb += 1
            rad = LUT[idx] >> (s >> 1)
            if rad * rad <= r_max << 52:
                errs += 1
            v = Fr(rad * rad, r_min << 52)
            if v > lam:
                lam = v
    return errs, lam, nb
def uniform_bounds(LAM, kappa=KAPPA, beta_lo=BETA_LO, grid=400):
    """Returns dict with z_k, e_k, F(b) evaluator, eps_or bound, completeness margins."""
    sL = sqrt_up(LAM)                                   # 1 + delta_L
    z = [None] * 12
    for k in range(0, 11):
        s = {j: sqrt_up(Fr(1) / kappa ** (j - k)) for j in range(k + 1, 12)}   # sqrt(beta_k / beta_j) bound
        m = {}
        for j in range(11, k, -1):
            m[j] = s[j] + Fr(3, 2) + sum((m[i] + 1) / 2 for i in range(j + 1, 12))
        h = {}
        for j in range(11, k, -1):
            ec_j = sum(m[i] for i in range(j + 1, 12)) + 1
            h[j] = sum(h[i] for i in range(j + 1, 12)) / 2 + ec_j / Q32 + sL * (s[j] + Fr(1, 64) + Fr(1, 1 << 20))
        ecs = sum(m.values()) + 1
        Es = sum(h.values())
        z[k] = (ecs + Es) / Q32                         # >= sqrt(beta_k) (ec_k + E_k) / 2^32
    z[11] = Fr(1, Q32)
    e = {}
    for k in range(2, 12):
        G = 1 + (1 << 19) * z[k]                        # >= sqrt(beta_k) * g_k
        e[k] = ((1 << 20) * G * (1 + Fr(1, 1 << 37)) + 2 * G * G + 2) / Q38
    D1 = sum(e.values())                                # Delta_1: levels 11..2
    W = sL * sqrt_up(1 + Fr(1, 1 << 12) + Fr(1, 1 << 38)) + z[1]
    def Fb(b):
        isb = sqrt_up(Fr(1) / b)
        ib = Fr(1) / b
        S = (2 * W * isb + 1) ** 2 + (1 + Fr(1, 1 << 26)) * W * W * ib + 2 * W * isb + 4 + (1 << 20) * W
        return S / Q38 + z[1] * (6 * W + 2 + z[1])
    def eps_at(Fv, b0):
        inner = (1 + Fv) * (1 + Fr(1, Q38)) + b0 / 4096
        return D1 + (sL * sqrt_up(inner) + z[0]) ** 2 - 1
    pts = [beta_lo * (Fr(1) / beta_lo) ** Fr(i, grid) if False else None for i in range(0)]  # placeholder
    pts = []
    r = sqrt_dn(sqrt_dn(sqrt_dn(Fr(1) / beta_lo)))      # (1/beta_lo)^(1/8)
    q = Fr(1) + Fr(1, grid // 8 if grid >= 16 else 2)
    b = beta_lo
    while b < 1:
        pts.append(b)
        b = b * q
    pts.append(Fr(1))
    worst = None
    for i in range(len(pts) - 1):
        val = eps_at(Fb(pts[i]), min(Fr(1), pts[i + 1] / kappa))
        if worst is None or val > worst[0]:
            worst = (val, pts[i], pts[i + 1])
    WR = {}
    WR["R1K1_max_log2"] = math.log2((1 << 39) * (2 * W + 1) * Q32 * z[1] + (1 << 33))   # R1*K1 + 2^33
    WR["r2_0_max_log2"] = math.log2((1 + Fb(beta_lo)) * (Fr(4096) / beta_lo + Fr(1, 1 << 26)) + 1)
    WR["dr1_max_log2"] = math.log2(Q38 * (2 * W + 2) + Q38)
    return dict(sL=sL, z=z, e=e, D1=D1, W=W, Fb=Fb, eps_at=eps_at, eps=worst[0], eps_where=worst[1:],
                F_lo=Fb(beta_lo), npts=len(pts), WR=WR)
G6 = [[753729, 251243, 251243, 251244, 251242, 251244], [251243, 1163069, -455913, -455913, 517661, 517662],
      [251243, -455913, 1163069, 189495, 189494, -455913], [251244, -455913, 189495, 1542060, 127747, -394165],
      [251242, 517661, 189494, 127747, 1542058, 455914], [251244, 517662, -455913, -394165, 455914, 1542060]]
KAP_L = Fr(74, 100)                              # LLL(99/100), |mu| <= 1/2: B_{j+1} >= 0.74 B_j
P_E_CHARGED = 3 << 12                           # Ir per exact node (charged; >= 2x the operation count below)
def _det(M):
    M = [[Fr(x) for x in r] for r in M]; d = Fr(1)
    for i in range(len(M)):
        d *= M[i][i]
        for r in range(i + 1, len(M)):
            f = M[r][i] / M[i][i]
            M[r] = [a - f * b for a, b in zip(M[r], M[i])]
    return d
def lemma_e():
    """Lemma E: for p < 2^61 and T >= T_lo, beta_k = B_k/T >= 1 is impossible at k <= 2.  If B_k >= T then
    B_j >= max(Bmin, 0.74^(j-k) T) for j >= k (Lovasz) and B_j >= Bmin for j < k, so prod_j B_j >= L_k; but
    prod_j B_j = det(Gram L_p) = p^2 1764^6 1000^6 det(G)^2 < U.  Returns exact log2(L_k / U) per k (> 0: impossible)."""
    U = Fr(1 << 122) * Fr(1764) ** 6 * Fr(1000) ** 6 * _det(G6) ** 2
    T_lo = Fr(1764) * Fr(58908, 100) ** 2 * 196608 - 1   # production R_a = 589.0871 >= 589.08
    out = []
    for k in range(12):
        L = BMIN ** k
        for j in range(k, 12):
            L *= max(BMIN, KAP_L ** (j - k) * T_lo)
        out.append(math.log2(L / U))
    return out
def exact_nodes_bound(UB, eta_max, lem):
    """Exact-body nodes per orbit.  At a visited prefix the stored remainder is <= R T, R = 1 + Delta_1 + 12 * 2^-38
    (Theorem B's per-level bounds; an exact level adds < 2^-38).  A fast level j places at most
    floor(2 sqrt(Lambda (R/beta_j + 2^-12 + R 2^-38)) + 2 eta_max + 2^-31) + 1 integers in [lo, hi]; an exact level keeps
    at most floor(2 sqrt(R/beta_j)) + 1.  If level k is exact (beta_k >= 1, so k >= 3 by Lemma E) then
    beta_j >= 0.74^(j-k) for j > k, so level k is entered at most prod_{j>k} c(0.74^(j-k)) times per orbit."""
    R = 1 + UB["D1"] + Fr(12, 1 << 38)
    def c(b):
        w = 2 * sqrt_up(UB["sL"] ** 2 * (R / b + Fr(1, 4096) + R / (1 << 38))) + 2 * eta_max + Fr(1, 1 << 31)
        return w.numerator // w.denominator + 1
    tot = 0
    for k in range(12):
        if lem[k] > 0:
            continue
        n = 1
        for j in range(k + 1, 12):
            n *= c(KAP_L ** (j - k))
        tot += n
    return tot
def fp_nodes_above_bound(UB, eta_max, lem):
    """Fixed-point nodes above an exact level, per orbit.  With level k exact (beta_k >= 1, k >= 3 by Lemma E), a level
    j > k is entered at most prod_{k < i ... } bounded as in exact_nodes_bound: level j is visited at most
    prod_{i > j} c(0.74^(i-k)) times.  Summed over j = k+1..11 (the root counts once) and over every k that can be exact."""
    R = 1 + UB["D1"] + Fr(12, 1 << 38)
    def c(b):
        w = 2 * sqrt_up(UB["sL"] ** 2 * (R / b + Fr(1, 4096) + R / (1 << 38))) + 2 * eta_max + Fr(1, 1 << 31)
        return w.numerator // w.denominator + 1
    tot = 0
    for k in range(12):
        if lem[k] > 0:
            continue
        for j in range(k + 1, 12):
            n = 1
            for i in range(j + 1, 12):
                n *= c(KAP_L ** (i - k))
            tot += n
    return tot
FP_NODE_MAX_IR = 105                            # largest FP node price counted from kernels_v31.s (FP_CL)
C_LAT_EXC_IR = int(2.0 ** 31.2)                 # floor(2^31.2) = 2,466,810,933, the model's per-orbit exact-level allowance
def swar21_max_coord_bound(eps_or):
    """Theorem C (3x21-bit biased SWAR cur zero-overflow bound across ALL FP levels k=0..11, input-uniform):
    At any visited FP node k (beta_j < 1 for all j >= k), cur_k = sum_{j=k}^{11} x_j b_j = u_orth + u_proj, where:
      - u_orth = sum_{m=k}^{11} y_m b_m^* satisfies Q12(u_orth) <= (1 + eps_or) * TMAX;
      - u_proj = sum_{m=0}^{k-1} (sum_{j=k}^{11} mu_{jm} x_j) b_m^* has orthogonal components bounded via
        |mu_{jm}| <= 1/2, |x_j| <= C_j * sqrt(TMAX / B_j) (C_j solving backward from Lovasz B_j >= KAPPA^{j-i} B_i,
        plus 1/sqrt(1/beta_j) <= 1 since beta_j < 1), and sqrt(B_m / B_j) <= KAPPA^{-(j-m)/2}.
    By Cauchy-Schwarz in the weighted G6 inner product, every coordinate of cur_k at every k=0..11 satisfies
    |cur_{k,c}| <= sqrt(Q12(cur_k) / 1000 * max_i (G6^{-1})_{ii}) < 2^20 = 1,048,576 (zero overflow in 3x21-bit SWAR)."""
    n = 6
    A = [[Fr(G6[i][j]) for j in range(n)] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]; A[r] = [A[r][j] - f * A[c][j] for j in range(2 * n)]
    max_ginv = max(A[i][n + i] for i in range(n))
    leaf_ub = sqrt_up((1 + eps_or) * Fr(TMAX, 1000) * max_ginv)
    # Input-uniform backward recurrence for C_j >= Xm[j] * sqrt(B_j / TMAX) on all FP levels (beta_j < 1):
    # rk_j <= sqrt(TMAX / B_j) + 2 -> rk_j * sqrt(B_j / TMAX) <= 1 + 2*sqrt(beta_j) < 3 (actually at fast levels beta_j < 1,
    # tighter: |y_j| <= sqrt((1+eps_or)*TMAX/B_j), so |x_j| <= sqrt((1+eps_or)*TMAX/B_j) + 1/2 + sum_{i>j} |x_i|/2).
    # For any visited node at level j: |y_j| <= sqrt((1+eps_or)/beta_j), and |x_j - (-sum_{i>j} mu_{ij} x_i)| <= |y_j| + 1/2.
    # Let D_j >= max |x_j| * sqrt(beta_j) (using beta_j < 1 => (1/2)*sqrt(beta_j) < 1/2):
    r_eff = sqrt_up(1 + eps_or)
    inv_sk = [sqrt_up(Fr(1) / (KAPPA ** d)) for d in range(12)]
    D = [Fr(0)] * 12
    for j in range(11, -1, -1):
        D[j] = r_eff + Fr(1, 2) + sum(Fr(1, 2) * D[i] * inv_sk[i - j] for i in range(j + 1, 12))
    # For each level k in 0..11, sum_{m=0}^{k-1} (1/2 * sum_{j=k}^{11} D[j] * inv_sk[j-m])^2 bounds Q12(u_proj) / TMAX!
    max_q_ratio = Fr(1) + eps_or
    for k in range(1, 12):
        proj_q = Fr(0)
        for m in range(k):
            comp = sum(Fr(1, 2) * D[j] * inv_sk[j - m] for j in range(k, 12))
            proj_q += comp * comp
        tot_q = (1 + eps_or) + proj_q
        if tot_q > max_q_ratio:
            max_q_ratio = tot_q
    # Also include k0_body step cur + x6*b0 (which is a leaf in Q12 <= (1+eps_or)*TMAX) and cur += b1 (next k=1 point):
    # cur_1 + b1 has norm <= sqrt(max_q_ratio*TMAX) + sqrt(TMAX) (since beta_1 < 1 and |mu_{10}| <= 1/2 => ||b1||^2 <= 5/4 TMAX).
    inter_ub = (sqrt_up(max_q_ratio) + sqrt_up(Fr(5, 4))) * sqrt_up(Fr(TMAX, 1000) * max_ginv)
    return leaf_ub, inter_ub
def exact_node_ops():
    """Operation count of one exact node from the operand widths (64-bit limbs): D_k <= det^2 < 2^485 (8 limbs),
    |lambda_jk| <= D_{k+1}/2 (8), S = sum_{j>k} lambda_jk x_j (<= 8 terms, 9 limbs), N = S + D_{k+1} x (9),
    Q_k = T D_k D_{k+1} < 2^1017 (16), rem Q_k and 2^38 N^2 (<= 19).  a x b limb product 4ab + 2(a+b)
    (mul, umulh, adds, adcs, loads/stores); n-limb add/sub/compare 3n; 64-bit-quotient division by an n-limb
    divisor (quotient estimated from the top limbs, multiply back, subtract, corrections) 4n + 2(n + 1) + 9n + 10."""
    mw = lambda a, b: 4 * a * b + 2 * (a + b)
    dv = lambda n: mw(1, n) + 9 * n + 10
    centre = 8 * (mw(8, 1) + 3 * 9) + dv(9)                      # S and floor(S / D_{k+1})
    per_cand = mw(1, 8) + 3 * 9 + mw(9, 9) + 3 * 19 + 3 * 19 + 10   # N, N^2, shift by 38, compare, loop
    child = 3 * 19 + dv(16) + 100                                  # remainder, ceil division, fast child state
    return centre + mw(1, 16) + 4 * per_cand + 3 * child + 20
def completeness_margin():
    """Lemma FD (completeness) reduces to   slack(R1, alpha, eta) = alpha*eta*(R1-1) - 2*R1*eta + 1
    + (3R1+1)(R1-1)/4 > 0  with alpha = 2^38 beta_1 >= 4 and eta = ec_1/2^32 in (0, 1/2).  It is affine in alpha
    and eta, so its minimum over the box is at a corner; evaluate exactly at all corners for R1 = 1..64 (for R1 >= 2 the slack is >= eta (2 R1 - 4) + 1 > 0 for every R1, KERNELS.md 1.3)."""
    worst = None
    rkm = sqrt_up(Fr(1) / BETA_LO) + Fr(3, 2)
    Xm = {}
    for j in range(11, 1, -1):
        Xm[j] = rkm + sum((Xm[i] + 1) / 2 for i in range(j + 1, 12))
    eta_max = (sum(Xm.values()) + 1) / Q32
    for R in range(1, 65):
        for alpha in (Fr(4), Fr(1 << 38)):
            for eta in (Fr(0), eta_max):
                sl = alpha * eta * (R - 1) - 2 * R * eta + 1 + Fr((3 * R + 1) * (R - 1), 4)
                if worst is None or sl < worst[0]:
                    worst = (sl, alpha, eta, R)
    return worst + (eta_max,)
def main():
    t0 = time.time()
    errs, LAMv, nb = lut_certificate()
    UB = uniform_bounds(LAMv)
    cm = completeness_margin()
    out = dict(lut7_buckets=nb, lut7_outward_errs=errs, lut7_lambda=float(LAMv), beta_lo_log2=math.log2(BETA_LO),
               z_max_log2=math.log2(max(UB["z"][:11])), delta1_levels_11_to_2=float(UB["D1"]), F_beta_lo=float(UB["F_lo"]),
               eps_or_bound=float(UB["eps"]), word_widths=UB["WR"], pad1_min_slack=float(cm[0]))
    lem = lemma_e()
    n_e = exact_nodes_bound(UB, cm[4], lem)
    n_fp = fp_nodes_above_bound(UB, cm[4], lem)
    ops_e = exact_node_ops()
    swar21_leaf_ub, swar21_inter_ub = swar21_max_coord_bound(UB["eps"])
    out.update(lemma_e_log2_margin_k0_k3=[round(v, 3) for v in lem[:4]], exact_levels_possible=[k for k in range(12) if lem[k] <= 0],
               exact_nodes_per_orbit_max=n_e, exact_nodes_log2=round(math.log2(n_e), 3), exact_node_ops=ops_e,
               exact_node_price_charged=P_E_CHARGED, exact_per_orbit_ir_log2=round(math.log2(n_e * P_E_CHARGED), 3),
               fp_nodes_above_exact_max=n_fp, exact_level_work_ir=n_e * P_E_CHARGED + n_fp * FP_NODE_MAX_IR,
               exact_level_allowance_ir=C_LAT_EXC_IR,
               swar21_leaf_coord_ub=round(float(swar21_leaf_ub), 2), swar21_max_coord_ub=round(float(swar21_inter_ub), 2),
               swar21_zero_overflow_ok=(swar21_inter_ub < (1 << 20)))
    out["ok"] = (errs == 0 and UB["eps"] <= Fr(23389, 10000000) and cm[0] > 0 and min(lem[:3]) > 0
                 and 2 * ops_e <= P_E_CHARGED and n_e * P_E_CHARGED + n_fp * FP_NODE_MAX_IR <= C_LAT_EXC_IR
                 and swar21_inter_ub < (1 << 20))
    out["seconds"] = round(time.time() - t0, 3)
    return out
if __name__ == "__main__":
    import json
    o = main(); print(json.dumps(o, default=str)); sys.exit(0 if o["ok"] else 1)
