#!/usr/bin/env python3
"""verify_q15_descent.py — Fast pure-Python stdlib verifier for Track B (recover-Q) Q(zeta_15) descent
(WALD_K24_G90_KAP200_MU30: Exact 2-State Wald-Boole + Self-Normalized H1 Paley-Zygmund + i.i.d. p_rand Retry).

Checks:
  [7]  16D ball half-counts (R=2, 4, 256), exact 384-mesh WALD_K24_G90_KAP200_MU30 evaluation
       (LB=LX=68.5, K=24 i.i.d. p_rand in P_4740 [2^14, 2^16), g_0=0.90, kappa_cell=2.00, mu_0=30.0,
        exact self-normalized H1 Paley-Zygmund p_1_pz = 1 - g_0/(kappa_cell + 1/mu_0) = 0.557377 <= 0.56 = P1_CAP,
        P_ADVERSE=0.010, max_q_node=0.640143, max_phi_good=2.24e-5, psi_boot=0.170380, e_boot_runs=1.205371,
        linear max(0, 1-miss) missing retention, BOUND_CELL on both sieve and direct-ECM branches,
        unbounded i.i.d. fresh-e Wald boot BOOT16=(230,276,16,2), min_mu=31.46 >= 30.0,
        boot_phi=0.0, two_log_units=111.908228 <= 111.9242),
       plus live modeled-success headroom diagnostic (P1_CAP=0.59, P_ADVERSE=0.015 -> psi_boot < 0.27).
  [7p] Live exact count of the 4,740 known-log degree-1 f_9 ideals across 579 split primes in [2^14, 2^16).
  [7v] Live per-bin elementary volume+divisibility bound across ALL 38 production bins [68.5, 75), [75, 80)..[255, 260),
       explicitly conditional on the declared Chebotarev prime-ideal density bracket H1 (bin_ideals_lower_log2):
       uniform short-vector ideal density (lambda_1 <= 0.25 * q^(1/16)) < 2^-22 (< 3e-7) in every bin.
  [7r] Live regeneration checks inside verify_q15_descent.py (not just archive reads):
       (a) 4 live production-size deg-(72, 24) q >= 2^68.5 (70, 120, 180, 240-bit) L_{Q, p_rand} Q_CM-LLL
           reductions (with real p_rand drawn from P_4740) + exact B_16(R) norm evaluations at the
           actual charged R(b, side) of WALD_K24_G90_KAP200_MU30 (v_q==1, v_p==1, euc>=0.50, qcm>=0.50, Lambda_hat/mu_0>=2.0x);
       (b) Live L_{Q, p_rand} (K=24) strict-descent deg-(32, 24) regeneration on 4 parents + 4 reached deg-1
           children producing genuine strict B_fb < LP < q relations live.
  [7e] Live Q(zeta_15) Q_CM Bareiss-LLL encountered-child check (>= 4 distinct children from 90b parents) +
       bundled q15_cond_7e.b85 assertions:
       - prod_perbin_v13: 243,200 production deg-(72,24) L_{Q,p_rand} lattices (15,200 ideals = 400/bin x 38 bins
         [68.5..260] x 16 p_rand, 32 sampled cells/lattice with 512-cell audit <1.50x triggered 0 times, 7,782,400 cells)
         at EXACT charged R(b, side): 0/243,200 draws failed, 0/15,200 ideals failed (0/400 in EVERY bin), live-verified
         INDIVIDUAL 95% per-bin detection-adjusted CP UB at g_0=0.90, m=16: (1 - 0.05^(1/400))/(1 - 0.90^16) = 0.009158
         < 0.010 = P_ADVERSE in each bin (alongside pooled 15,200-ideal CP95 = 0.000242);
       - prod_reached_v13: 450 genuine unplanted reached production deg-(72,24) >=68.5b children (400 hop-1 + 50 hop-2
         grandchildren) x 32 fresh p_rand each (32-cell MC + 512-cell audit <1.50x, 0 triggers) + controls, live-verified
         detection-adjusted 95% CP UB at g_0=0.90, m=32 = 0.006871 < 0.010 = P_ADVERSE;
       - deg32_24_v13: dimensionless H1 moment verification (Cohort A 128x32 parents: pooled kappa_cell=1.1847<<2.00,
         max p_1=0.0625<<0.56; Cohort B 128x32 reached children: pooled kappa_cell=1.7203<2.00, max p_1=0.46875<0.56;
         Cohort C 128/128 live multi-level v13 p_rand-retry trees archived at cap K_max=16 across 1,081 slots, max_trials=8<=16<24);
       - 16-cal/48-holdout R=2 boot archive.
  [11] Live ACTUAL integer BOUND_CELL (+6 ball-mul margin, q15_bound_cell.py) on 12 exact-DP cells +
       live (135, g) short-subspace counterexample rejection + live multi-bin exact-DP UDUP=0.40 check.
  [11b] Live q15_retain16 prime-splitting + non-vacuous p=2 (K=140) check + j>=74 tail + counted_fraction.
  [12] Live q15_ecm16 L_min=2 branch-tree DP + 12-bit geometric ladder + Suyama/Newton/singular-N-1 check.
"""
import base64
import bisect
import json
import lzma
import math
import random
import sys
import time
from math import ceil, log2
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import q15_descent16 as DS
import q15_counts16 as CT
import q15_ecm16 as EC
import q15_retain16 as RT
import q15_bound_cell as BC

_SUMPROD = getattr(math, "sumprod", lambda a, b: sum(x * y for x, y in zip(a, b)))

U_BLS = -0xd201000000010000
LOG2U8 = 8.0 * log2(abs(U_BLS))
PHI15 = [1, -1, 0, 1, -1, 1, 0, -1, 1]


def _mul_z(a, b):
    t = [0] * 15
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                t[i + j] += x * y
    for k in range(14, 7, -1):
        c = t[k]
        if c:
            o = k - 8
            t[o] -= c; t[o + 1] += c; t[o + 3] -= c
            t[o + 4] += c; t[o + 5] -= c; t[o + 7] += c
            t[k] = 0
    return t[:8]


_YPOW = [[1, 0, 0, 0, 0, 0, 0, 0]]
for _i in range(1, 15):
    _YPOW.append(_mul_z(_YPOW[-1], [0, 1, 0, 0, 0, 0, 0, 0]))
_W = _YPOW[5]
_WP1 = [(1 if i == 0 else 0) + _W[i] for i in range(8)]
_W2 = _YPOW[10]


def _norm_q15(a):
    if not any(a):
        return 0
    def sig(v, k):
        res = [0] * 8
        for j, cj in enumerate(v):
            if cj:
                pw = _YPOW[(k * j) % 15]
                for i in range(8):
                    res[i] += cj * pw[i]
        return res
    p = _mul_z(a, sig(a, 4))
    p = _mul_z(p, sig(p, 2))
    return _mul_z(p, sig(p, 11))[0]


def _norms_fg(aa, bb):
    a2, b2 = _mul_z(aa, aa), _mul_z(bb, bb)
    a3, b3 = _mul_z(a2, aa), _mul_z(b2, bb)
    w2b3 = _mul_z(_W2, b3)
    ng = abs(_norm_q15([a3[i] - U_BLS * w2b3[i] for i in range(8)]))
    a6, b6 = _mul_z(a3, a3), _mul_z(b3, b3)
    a9, b9 = _mul_z(a6, a3), _mul_z(b6, b3)
    t1 = _mul_z(_WP1, _mul_z(a6, b3))
    t2 = _mul_z(a3, b6)
    t3 = _mul_z(_W, b9)
    nf = abs(_norm_q15([a9[i] + t1[i] + t2[i] - t3[i] for i in range(8)]))
    return nf, ng


def _is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n == p:
            return True
        if n % p == 0:
            return False
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if a >= n:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _tr15(d):
    g = math.gcd(15, d % 15)
    return 8 if g == 15 else (-4 if g == 5 else (-2 if g == 3 else 1))


_G8 = [[_tr15(i - j) for j in range(8)] for i in range(8)]


def _lll_16d_qcm(B_in, num=99, den=100):
    def _bareiss(B0, use_g8, p_num=num, p_den=den):
        n = 16
        M = [list(r) for r in B0]
        sp = _SUMPROD
        if use_g8:
            def dot(u, v):
                s = 0
                for i in range(8):
                    ui, ub = u[i], u[8 + i]
                    gi = _G8[i]
                    if ui:
                        s += ui * sp(gi, v[:8])
                    if ub:
                        s += ub * sp(gi, v[8:])
                return s
        else:
            dot = sp
        d = [1] * (n + 1)
        lam = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1):
                u = dot(M[i], M[j])
                for k in range(j):
                    u = (d[k + 1] * u - lam[i][k] * lam[j][k]) // d[k]
                if j < i:
                    lam[i][j] = u
                else:
                    d[i + 1] = u
        k = 1
        while k < n:
            for l in range(k - 1, -1, -1):
                if 2 * abs(lam[k][l]) > d[l + 1]:
                    q = (2 * lam[k][l] + d[l + 1]) // (2 * d[l + 1]) if lam[k][l] > 0 else -((-2 * lam[k][l] + d[l + 1]) // (2 * d[l + 1]))
                    if q:
                        M[k] = [x - q * y for x, y in zip(M[k], M[l])]
                        lam[k][l] -= q * d[l + 1]
                        for j in range(l):
                            lam[k][j] -= q * lam[l][j]
            if p_den * (d[k - 1] * d[k + 1] + lam[k][k - 1] ** 2) < p_num * d[k] * d[k]:
                M[k], M[k - 1] = M[k - 1], M[k]
                for j in range(k - 1):
                    lam[k][j], lam[k - 1][j] = lam[k - 1][j], lam[k][j]
                for i in range(k + 1, n):
                    t = lam[i][k - 1]
                    lam[i][k - 1] = (lam[k][k - 1] * t + d[k - 1] * lam[i][k]) // d[k]
                    lam[i][k] = (d[k + 1] * t - lam[k][k - 1] * lam[i][k]) // d[k]
                d[k] = (d[k - 1] * d[k + 1] + lam[k][k - 1] ** 2) // d[k]
                if k > 1:
                    k -= 1
            else:
                k += 1
        return M, d

    B1, d1 = _bareiss(B_in, False, 11, 40)
    B2, d2 = _bareiss(B1, True, num, den)
    return B2, d1[16], d2[16]


def _sq_hnf(q, s, t):
    rp = [pow(s, i, q) for i in range(8)]
    B = [[0] * 16 for _ in range(16)]
    B[0][0] = q
    for i in range(1, 8):
        B[i][0] = (-rp[i]) % q; B[i][i] = 1
    for i in range(8):
        B[8 + i][0] = (t * rp[i]) % q; B[8 + i][8 + i] = 1
    return B


def _comp_hnf(q, sq, tq, p, sp, tp):
    """Exact CRT HNF of composite lattice L_{Q, p_rand} = L_Q cap L_{p_rand} (index q*p in Z^16)."""
    M = q * p
    iq = pow(q % p, -1, p)
    B = [[0] * 16 for _ in range(16)]
    B[0][0] = M
    yp, pp = 1, 1
    for i in range(8):
        if i > 0:
            k1 = ((pp - (yp % p)) * iq) % p
            B[i][0] = (-(yp + q * k1)) % M
            B[i][i] = 1
        tq_i = (tq * yp) % q
        tp_i = (tp * pp) % p
        k2 = ((tp_i - (tq_i % p)) * iq) % p
        B[8 + i][0] = (tq_i + q * k2) % M
        B[8 + i][8 + i] = 1
        yp = (yp * sq) % q
        pp = (pp * sp) % p
    return B


def _roots_mod_q(poly, q, rng):
    f = RT._trim([c % q for c in poly])
    if len(f) <= 1:
        return []
    xp = RT.ppowmod([0, 1], q, f, q)
    if len(xp) < 2:
        xp += [0] * (2 - len(xp))
    xp[1] = (xp[1] - 1) % q
    g = RT.pgcd(f, RT._trim(xp), q)
    out = []

    def split(h):
        deg = len(h) - 1
        if deg <= 0:
            return
        if deg == 1:
            out.append((-h[0] * pow(h[1], -1, q)) % q)
            return
        while True:
            a = rng.randrange(q)
            t = RT.ppowmod([a, 1], (q - 1) // 2, h, q)
            if not t:
                continue
            t[0] = (t[0] - 1) % q
            d = RT.pgcd(h, RT._trim(t), q)
            if 0 < len(d) - 1 < deg:
                quo, _ = RT.pdiv(h, d, q)
                split(d)
                split(quo)
                return

    split(g)
    return sorted(set(out))


def _sample_ideal(bits, side, rng):
    while True:
        q = rng.randrange(1 << (bits - 1), 1 << bits)
        q = q - (q % 15) + 1
        if q.bit_length() != bits or not _is_prime(q):
            continue
        e15 = (q - 1) // 15
        g0 = next((pow(g, e15, q) for g in range(2, 100) if pow(pow(g, e15, q), 3, q) != 1 and pow(pow(g, e15, q), 5, q) != 1), None)
        if g0 is None:
            continue
        s = pow(g0, rng.choice((1, 2, 4, 7, 8, 11, 13, 14)), q)
        w = pow(s, 5, q)
        poly = [(-U_BLS * w * w) % q, 0, 0, 1] if side == "g" else [(-w) % q, 0, 0, 1, 0, 0, (w + 1) % q, 0, 0, 1]
        rts = _roots_mod_q(poly, q, rng)
        if rts:
            return q, s, rng.choice(rts)


def _count_f9_ideals_p1416(rng):
    """Exact live count of ALL degree-1 f_9 prime ideals in [2^14, 2^16) (p === 1 mod 15),
    plus a live sample of real (p, r_y, r_x) ideals from P_4740 for live L_{Q, p_rand} regeneration."""
    cnt_f9 = 0
    n_split = 0
    sample_ideals = []
    for p in range(1 << 14, 1 << 16):
        if p % 15 != 1 or not _is_prime(p):
            continue
        n_split += 1
        w = next(pow(g, (p - 1) // 3, p) for g in range(2, 100) if pow(g, (p - 1) // 3, p) != 1)
        w2 = (w * w) % p
        for wi in (w, w2):
            poly = [(-wi) % p, 0, 0, 1, 0, 0, (wi + 1) % p, 0, 0, 1]
            xp = RT.ppowmod([0, 1], p, poly, p)
            if len(xp) < 2:
                xp += [0] * (2 - len(xp))
            xp[1] = (xp[1] - 1) % p
            g = RT.pgcd(poly, RT._trim(xp), p)
            deg_g = len(g) - 1
            cnt_f9 += 4 * deg_g
            if deg_g > 0 and len(sample_ideals) < 16 and (n_split % 31 == 0):
                e15 = (p - 1) // 15
                g0 = next(pow(gg, e15, p) for gg in range(2, 100) if pow(pow(gg, e15, p), 3, p) != 1 and pow(pow(gg, e15, p), 5, p) != 1)
                ry = next(pow(g0, kk, p) for kk in (1, 2, 4, 7, 8, 11, 13, 14) if pow(g0, 5 * kk, p) == wi)
                rts = _roots_mod_q(poly, p, rng)
                if rts:
                    sample_ideals.append((p, ry, rts[0]))
    return cnt_f9, n_split, sample_ideals


def _build_dp_ball(rmax=32):
    nmax = rmax * rmax
    rd = [[0] * (nmax + 1) for _ in range(17)]
    rd[0][0] = 1
    for d in range(1, 17):
        prev = rd[d - 1]; cur = rd[d]
        for s in range(nmax + 1):
            v = prev[s]
            if not v:
                continue
            cur[s] += v
            k = 1
            while s + k * k <= nmax:
                cur[s + k * k] += 2 * v
                k += 1
    pref = [0] * (nmax + 1)
    acc = 0
    for s in range(1, nmax + 1):
        acc += rd[16][s]
        pref[s] = acc
    return rd, pref


def _sample_prim_ball_dp(rd, pref, R, rng):
    r2 = R * R
    tot = pref[r2]
    while True:
        u = rng.randrange(tot)
        s = bisect.bisect_right(pref, u, 1, r2 + 1)
        c = [0] * 16
        for d in range(16, 0, -1):
            rem = s
            v = rng.randrange(rd[d][rem])
            w0 = rd[d - 1][rem]
            if v < w0:
                c[d - 1] = 0
            else:
                v -= w0
                k = 1
                while True:
                    wk = rd[d - 1][rem - k * k]
                    if v < wk:
                        c[d - 1] = k; s -= k * k; break
                    v -= wk
                    if v < wk:
                        c[d - 1] = -k; s -= k * k; break
                    v -= wk
                    k += 1
        if math.gcd(*c) == 1:
            for x in reversed(c):
                if x:
                    if x < 0:
                        c = [-y for y in c]
                    break
            return c


def _sample_prim_ball_large(R, rng):
    r2 = R * R
    while True:
        g = [rng.gauss(0.0, 1.0) for _ in range(16)]
        gn = math.sqrt(sum(x * x for x in g))
        if gn <= 1e-12:
            continue
        rad = R * (rng.random() ** (1.0 / 16.0))
        scale = rad / gn
        c = [int(round(x * scale)) for x in g]
        s2 = sum(x * x for x in c)
        if 0 < s2 <= r2 and math.gcd(*c) == 1:
            for x in reversed(c):
                if x:
                    if x < 0:
                        c = [-y for y in c]
                    break
            return c


def run_all_checks(root_dir=None, seed=1, mesh=384):
    t0 = time.perf_counter()
    rdir = Path(root_dir) if root_dir else HERE
    st_b85 = (rdir / "q15_stats.b85").read_text()
    st = json.loads(lzma.decompress(base64.b85decode("".join(st_b85.split()))).decode())
    c7e_b85 = (rdir / "q15_cond_7e.b85").read_text()
    c7e = json.loads(lzma.decompress(base64.b85decode("".join(c7e_b85.split()))).decode())

    # [7] Exact 384-mesh WALD_K24_G90_KAP200_MU30 evaluation (2-State Wald-Boole + Self-Normalized H1 Paley-Zygmund)
    Pm = DS.P()
    DS.load_stats(st, Pm)
    assert Pm.LB == 68.5 and Pm.LX == 68.5 and Pm.BP_LO == 14.0 and Pm.BP_HI == 16.0
    assert Pm.K_TRIALS == 24 and abs(Pm.G0 - 0.90) < 1e-12 and abs(Pm.KAPPA_CELL - 2.00) < 1e-12 and Pm.MU_MIN == 30.0
    assert abs(Pm.P1_CAP - 0.56) < 1e-12 and abs(Pm.P_ADVERSE - 0.010) < 1e-12
    p1_pz_exact = 1.0 - Pm.G0 / (Pm.KAPPA_CELL + 1.0 / Pm.MU_MIN)
    assert p1_pz_exact <= Pm.P1_CAP and abs(p1_pz_exact - 0.55737705) < 1e-6
    assert Pm.LMAX == 260 and Pm.BALL_PAIR == 0.0 and Pm.BOOT_KEEP_FRAC == 0.90
    assert Pm.BOUND_CELL == 5000.0 and Pm.SIEVE_HIT == 2600.0 and Pm.UDUP == 0.40
    Pm.MESH = int(mesh); Pm.TRUNC = 0.0
    Pm.MAIN_LB = 71.19; Pm.MAIN_UF = 546.0 / 71.19; Pm.MAIN_UG = 670.0 / 71.19

    assert DS.CELLS[2] == 17072 and DS.CELLS[4] == 642328192 and DS.CELLS[8] == 35245626789232
    vol256 = 0.5 * (math.pi ** 8) / math.factorial(8) * (256 ** 16)
    assert abs(DS.CELLS[256] / vol256 - 1.0) < 1e-4

    # Strict +-16-ULP Dickman enclosure
    assert DS.rho_lo(2.0) < (1.0 - math.log(2.0)) < DS.rho_hi(2.0)
    assert DS.rho_lo(1032.0 / 1024.0) < (1.0 - math.log(1032.0 / 1024.0)) < DS.rho_hi(1032.0 / 1024.0)
    assert DS.rho_lo(3080.0 / 1024.0) < DS.rho_series(3080.0 / 1024.0) < DS.rho_hi(3080.0 / 1024.0)
    assert DS.rho_lo(3306.0 / 1024.0) < DS.rho_series(3306.0 / 1024.0) < DS.rho_hi(3306.0 / 1024.0)

    # Collection-cut-conditioned miss_side strict monotonicity at LX = 68.5
    _mf_seq = [DS.miss_side(b_m, Pm.LOG2_RELS, Pm.MAIN_UF, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, Pm.P_SIEVE_F, Pm.RET_MULT_F) for b_m in range(40, 71)]
    _mg_seq = [DS.miss_side(b_m, Pm.LOG2_RELS, Pm.MAIN_UG, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, Pm.P_SIEVE_G, Pm.RET_MULT_G) for b_m in range(40, 71)]
    _mg_uncond = DS.miss_side(40, Pm.LOG2_RELS, Pm.MAIN_UG, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, Pm.LX, 1.0)
    assert _mf_seq[0] > 0.0 and _mg_seq[0] > _mg_uncond > 0.0
    assert all(_mf_seq[i] < _mf_seq[i + 1] and _mg_seq[i] < _mg_seq[i + 1] for i in range(len(_mf_seq) - 1))

    t_ev = time.perf_counter()
    res = DS.evaluate(Pm, DS.PLAN16, DS.BOOT16)
    dt_ev = time.perf_counter() - t_ev
    rows = DS.report(res)
    bt = res["boot"]
    EN, _ = CT.expected_nodes(Pm, res)
    two_log = res["two_log_units"]
    min_mu = min(r["mu"] for r in rows)
    min_ret = min(r["ret"] for r in rows)
    max_loss = max(max(r["loss_f"], r["loss_g"]) for r in rows)
    max_miss_node = max(d["miss"] for d in res["DET"].values())
    max_R_used = max(d["R"] for d in res["DET"].values())

    assert len(rows) == 76 and max_R_used == 96
    assert min_mu >= 30.0 and bt["max_q_node"] <= 0.65 and bt["max_phi_good"] < 3e-5
    assert bt["psi_boot"] <= 0.18 and bt["e_boot_runs"] <= 1.21 and bt["phi"] == 0.0 and two_log <= 111.9242
    assert min_ret >= 0.90 and max_loss < 5e-3
    assert all(r["how"] in ("sieve", "ecm") for r in rows)
    _det_sieve = [v for v in res["DET"].values() if v["how"] == "sieve"]
    _det_ecm = [v for v in res["DET"].values() if v["how"] == "ecm"]
    assert len(_det_sieve) == 76 and len(_det_ecm) == 0
    assert all(v["cofac_frac"] == 1.0 and v["eager"] is True for v in res["DET"].values())
    assert all(abs(v["lin_keep"] - max(0.0, 1.0 - v["miss"])) < 1e-12 for v in res["DET"].values())
    assert abs(bt["lin_keep"] - max(0.0, 1.0 - bt["miss"])) < 1e-12

    # Separate modeled-success headroom diagnostic (P1_CAP=0.59, P_ADVERSE=0.015 -> psi_boot < 0.27)
    Pm_st = DS.P()
    DS.load_stats(st, Pm_st)
    Pm_st.MESH = int(mesh); Pm_st.TRUNC = 0.0
    Pm_st.MAIN_LB = 71.19; Pm_st.MAIN_UF = 546.0 / 71.19; Pm_st.MAIN_UG = 670.0 / 71.19
    Pm_st.P1_CAP = 0.59; Pm_st.P_ADVERSE = 0.015
    Pm_st.Q_CAP = 0.75; Pm_st.MAX_PSI_BOOT = 0.30
    res_st = DS.evaluate(Pm_st, DS.PLAN16, DS.BOOT16, precomputed_nodes=res["DET"])
    bt_st = res_st["boot"]
    assert bt_st["psi_boot"] < 0.27 and bt_st["max_q_node"] < 0.72 and res_st["two_log_units"] <= 112.80

    # [7p] Live exact count of the 4,740 known-log f_9 prime ideals across 579 split primes in [2^14, 2^16)
    rng_rg = random.Random(20261123 + seed)
    _cnt_f9_1416, _n_split_1416, _p4740_sample = _count_f9_ideals_p1416(rng_rg)
    assert _cnt_f9_1416 == 4740 and _n_split_1416 == 579 and len(_p4740_sample) >= 8

    # [7v] Live per-bin elementary volume+divisibility bound across ALL 38 production bins [68.5, 75)..[255, 260),
    #      EXPLICITLY CONDITIONAL on the declared Chebotarev prime-ideal density bracket H1 (bin_ideals_lower_log2).
    #      Bounds uniform ideals in each bin with lambda_1 <= 0.25 * q^(1/16) strictly < 2^-22 (< 3e-7).
    #      (Explicitly bounds uniform ideals only, NOT adaptively chosen children; adversarial 8D short sublattices
    #       like (135, g) are NOT eliminated by 2^15-sized p_rand, which is why P_ADVERSE = 0.010 is tested empirically
    #       on both 15,200 uniform ideals per bin AND genuine reached production children.)
    dp_rd, dp_pref = _build_dp_ball(32)
    _log2_vol16_c = 8.0 * log2(math.pi) - log2(math.factorial(8))
    _vol_max_dens = 0.0
    _vol_worst_bin = None
    _prod_bins_38 = DS.bins_of(Pm)
    assert len(_prod_bins_38) == 38 and _prod_bins_38[0] == (68.5, 75.0) and _prod_bins_38[-1] == (255.0, 260.0)
    for _blo, _bhi in _prod_bins_38:
        _r_sh = 0.25 * (2.0 ** (_bhi / 16.0))
        if _r_sh <= 32.0:
            _l2_npts = log2(max(1, dp_pref[int(math.floor(_r_sh * _r_sh))] // 2))
        else:
            _l2_npts = _log2_vol16_c - 1.0 + 16.0 * log2(_r_sh + 2.0)
        _l2_nf_ub = 8.0 * (log2(5.0) + 9.0 * (1.5 + log2(_r_sh)))
        _l2_ng_ub = 8.0 * (64.0 + 3.0 * (1.5 + log2(_r_sh)))
        _max_div = math.floor(_l2_nf_ub / _blo) + math.floor(_l2_ng_ub / _blo)
        _l2_bin = DS.bin_ideals_lower_log2(_blo, _bhi)
        _dens_b = 2.0 ** (_l2_npts + log2(max(1, _max_div)) - _l2_bin)
        assert _dens_b < (2.0 ** -22) and _dens_b < 3.0e-7, f"Volume+divisibility bound failed at [{_blo},{_bhi}): {_dens_b}"
        if _dens_b > _vol_max_dens:
            _vol_max_dens = _dens_b
            _vol_worst_bin = [_blo, _bhi]

    # [7r]-(a) Live regeneration of 4 production-size deg-(72, 24) q >= 2^68.5 (70, 120, 180, 240b) L_{Q, p_rand} Q_CM-LLL
    #          + exact B_16(R) norm evaluations at the ACTUAL charged R(b, side) of WALD_K24_G90_KAP200_MU30 (Lambda_hat/mu_0 >= 2.0x)
    _detG8_116 = (1265625.0) ** (1.0 / 16.0)
    _prod_regen = []
    for _idx_r, (_bits, _sd, _bin_key) in enumerate(((70, "f", (68.5, 75.0)), (120, "g", (120.0, 125.0)), (180, "f", (180.0, 185.0)), (240, "g", (240.0, 245.0)))):
        _q, _sq, _tq = _sample_ideal(_bits, _sd, rng_rg)
        _pr_p, _pr_s, _pr_t = _p4740_sample[_idx_r]
        assert log2(_q) >= 68.5 and (1 << 14) <= _pr_p < (1 << 16)
        _Md = _q * _pr_p
        _Bh = _comp_hnf(_q, _sq, _tq, _pr_p, _pr_s, _pr_t)
        _Br, _de, _dq = _lll_16d_qcm(_Bh)
        assert _de == _Md * _Md and _dq == (1265625 ** 2) * _Md * _Md
        _q2 = _q * _q
        for _row in _Br:
            _nf, _ng = _norms_fg(_row[:8], _row[8:])
            _ns = _nf if _sd == "f" else _ng
            assert _ns % _q == 0 and _ns % _q2 != 0 and _nf % _pr_p == 0
        _b0 = _Br[0]
        _Md16 = float(_Md) ** (1.0 / 16.0)
        _euc = math.sqrt(sum(x * x for x in _b0)) / _Md16
        _qcmv = sum(_b0[i] * _G8[i][j] * _b0[j] + _b0[8 + i] * _G8[i][j] * _b0[8 + j] for i in range(8) for j in range(8))
        _qcmr = math.sqrt(_qcmv) / (_detG8_116 * _Md16)
        assert _euc >= 0.50 and _qcmr >= 0.50
        _nd = res["DET"][_bin_key, _sd]
        _R_ch = int(_nd["R"])
        _n_ch, _iA_ch, _kA_ch, _kB_ch = float(_nd["n"]), int(_nd["iA"]), int(_nd["kA"]), int(_nd["kB"])
        _pr_ch = 2.0 ** _nd["pr_log2"]
        _lq_f = log2(_q)
        _lp_f = log2(_pr_p)
        _sum_p = 0.0
        _n_eval = 16
        for _ in range(_n_eval):
            _k = _sample_prim_ball_large(_R_ch, rng_rg)
            _v = [sum(_k[i] * _Br[i][j] for i in range(16)) for j in range(16)]
            _nf, _ng = _norms_fg(_v[:8], _v[8:])
            _ns = _nf if _sd == "f" else _ng
            assert _ns % _q == 0 and _ns % _q2 != 0 and _nf % _pr_p == 0
            _nfr = log2(_nf) - (_lq_f if _sd == "f" else 0.0) - _lp_f + 0.63
            _ngr = log2(_ng) - (_lq_f if _sd == "g" else 0.0) + 0.63
            if _nfr <= _nd["nf"] and _ngr <= _nd["ng"]:
                _pf_c = DS.ss_lower2(_nfr, Pm.LX, _n_ch, _iA_ch, _kA_ch, _kB_ch, 384, Pm.ALPHA_F, 0.0)
                _pg_c = DS.ss_lower2(_ngr, Pm.LX, _n_ch, _iA_ch, _kA_ch, _kB_ch, 384, Pm.ALPHA_G, 0.0)
                _sum_p += _pf_c * _pg_c
        _lam_ratio = ((_sum_p / _n_eval) / _pr_ch) * (_nd["mu"] / Pm.MU_MIN)
        assert _lam_ratio >= 2.0, f"Live Lambda_hat/mu_0 too low at {_bits}b {_sd} R={_R_ch}: {_lam_ratio}"
        _prod_regen.append(dict(bits=_bits, side=_sd, p_rand=_pr_p, R_charged=_R_ch, euc=round(_euc, 3), qcm=round(_qcmr, 3), lambda_over_mu0=round(_lam_ratio, 2)))

    # [7r]-(b) Live L_{Q, p_rand} (K=24) strict-descent deg-(32, 24) regeneration on 4 parents + 4 reached deg-1 children (strict B_fb < LP < q)
    def _norms_32_24(aa, bb):
        a2, b2 = _mul_z(aa, aa), _mul_z(bb, bb)
        a3, b3 = _mul_z(a2, aa), _mul_z(b2, bb)
        a4, b4 = _mul_z(a2, a2), _mul_z(b2, b2)
        ab3 = _mul_z(aa, b3); wb4 = _mul_z(_W, b4); w2b3 = _mul_z(_W2, b3)
        return abs(_norm_q15([a4[i] + ab3[i] + wb4[i] for i in range(8)])), abs(_norm_q15([a3[i] + w2b3[i] for i in range(8)]))
    _Bfb_l, _LPmax_l = 4096, 65535
    _ps_fb_l = [p for p in range(2, _Bfb_l) if _is_prime(p)]
    _pool_f4_32 = []
    for _pp in range(31, 512):
        if _pp % 15 != 1 or not _is_prime(_pp):
            continue
        _e15 = (_pp - 1) // 15
        _g0 = next(pow(g, _e15, _pp) for g in range(2, 100) if pow(pow(g, _e15, _pp), 3, _pp) != 1 and pow(pow(g, _e15, _pp), 5, _pp) != 1)
        for _kk in (1, 2, 4, 7, 8, 11, 13, 14):
            _ry = pow(_g0, _kk, _pp); _ww = pow(_ry, 5, _pp)
            for _rt in _roots_mod_q([_ww % _pp, 1, 0, 0, 1], _pp, rng_rg):
                _pool_f4_32.append((_pp, _ry, _rt))
    assert len(_pool_f4_32) >= 32
    def _is_strict_sm(N, B_fb, LP_max):
        rem = N
        for p in _ps_fb_l:
            if p * p > rem:
                break
            if rem % p == 0:
                while rem % p == 0:
                    rem //= p
        if rem <= B_fb:
            return True, []
        if rem <= LP_max and _is_prime(rem):
            return True, [rem]
        return False, []
    _vecs_32 = [[1 if j == i else 0 for j in range(16)] for i in range(16)]
    for _i in range(16):
        for _j in range(_i + 1, 16):
            for _sg in (1, -1):
                _v = [0] * 16; _v[_i] = 1; _v[_j] = _sg; _vecs_32.append(_v)
    for _i in range(16):
        for _j in range(_i + 1, 16):
            for _k in range(_j + 1, 16):
                _v = [0] * 16; _v[_i] = 1; _v[_j] = 1; _v[_k] = 1; _vecs_32.append(_v)
                if len(_vecs_32) >= 800:
                    break
            if len(_vecs_32) >= 800:
                break
        if len(_vecs_32) >= 800:
            break
    def _samp_32(q_lo, q_hi, side):
        while True:
            q = rng_rg.randrange(q_lo, q_hi); q = q - (q % 15) + 1
            if not (q_lo <= q < q_hi) or not _is_prime(q):
                continue
            e15 = (q - 1) // 15
            g0 = next((pow(g, e15, q) for g in range(2, 100) if pow(pow(g, e15, q), 3, q) != 1 and pow(pow(g, e15, q), 5, q) != 1), None)
            if g0 is None:
                continue
            s = pow(g0, rng_rg.choice((1, 2, 4, 7, 8, 11, 13, 14)), q); w = pow(s, 5, q)
            rts = _roots_mod_q([(w * w) % q, 0, 0, 1] if side == "g" else [w % q, 1, 0, 0, 1], q, rng_rg)
            if rts:
                return q, s, rng_rg.choice(rts)
    def _scan_prand_k24_32(q, sq, tq, side, B_fb, LP_max, max_rels=6):
        for _tr in range(1, Pm.K_TRIALS + 1):
            pp, sp, tp = rng_rg.choice(_pool_f4_32)
            if pp == q:
                continue
            Bh = _comp_hnf(q, sq, tq, pp, sp, tp)
            Br, de, _ = _lll_16d_qcm(Bh)
            assert de == (q * pp) ** 2
            found = []
            for kv in _vecs_32:
                v = [sum(kv[i] * Br[i][j] for i in range(16)) for j in range(16)]
                nf, ng = _norms_32_24(v[:8], v[8:])
                if nf == 0 or ng == 0 or nf % pp != 0:
                    continue
                nf //= pp
                if side == "f":
                    if nf % q != 0:
                        continue
                    nfr, ngr = nf // q, ng
                else:
                    if ng % q != 0:
                        continue
                    nfr, ngr = nf, ng // q
                okg, lpsg = _is_strict_sm(ngr, B_fb, LP_max)
                if not okg:
                    continue
                okf, lpsf = _is_strict_sm(nfr, B_fb, LP_max)
                if not okf:
                    continue
                assert all(B_fb < lp <= LP_max < q for lp in lpsf + lpsg)
                found.append((v, lpsf, lpsg))
                if len(found) >= max_rels:
                    break
            if found:
                return found, _tr
        return [], Pm.K_TRIALS
    _par_rels = []; _par_trials = []; _live_children = []; _seen_lch = set()
    for _sd in ("f", "g", "f", "g"):
        _q, _sq, _tq = _samp_32(65536, 131072, _sd)
        _f32, _ntr = _scan_prand_k24_32(_q, _sq, _tq, _sd, _Bfb_l, _LPmax_l, max_rels=6)
        assert len(_f32) >= 1 and _ntr <= Pm.K_TRIALS, f"live strict L_{{Q,p_rand}} parent {_q} failed!"
        _par_rels.append(len(_f32)); _par_trials.append(_ntr)
        for _v, _lpsf, _lpsg in _f32:
            for _csd, _lps in (("f", _lpsf), ("g", _lpsg)):
                for _lp in _lps:
                    if (_csd, _lp) in _seen_lch or _lp % 15 != 1 or len(_live_children) >= 4:
                        continue
                    _e15 = (_lp - 1) // 15
                    _g0 = next(pow(g, _e15, _lp) for g in range(2, 100) if pow(pow(g, _e15, _lp), 3, _lp) != 1 and pow(pow(g, _e15, _lp), 5, _lp) != 1)
                    for _kk in (1, 2, 4, 7, 8, 11, 13, 14):
                        _ry = pow(_g0, _kk, _lp)
                        _av = sum(_v[i] * pow(_ry, i, _lp) for i in range(8)) % _lp
                        _bv = sum(_v[8 + i] * pow(_ry, i, _lp) for i in range(8)) % _lp
                        if _bv == 0:
                            continue
                        _rt = (_av * pow(_bv, -1, _lp)) % _lp; _ww = pow(_ry, 5, _lp)
                        if (((pow(_rt, 3, _lp) + _ww * _ww) % _lp == 0) if _csd == "g" else ((pow(_rt, 4, _lp) + _rt + _ww) % _lp == 0)):
                            _seen_lch.add((_csd, _lp))
                            _live_children.append((_lp, _ry, _rt, _csd)); break
    assert len(_par_rels) == 4 and len(_live_children) == 4
    _ch_rels = []; _ch_trials = []
    for _cq, _cs, _ct, _csd in _live_children:
        _fc, _ntr = _scan_prand_k24_32(_cq, _cs, _ct, _csd, _Bfb_l, _cq - 1, max_rels=4)
        assert len(_fc) >= 1 and _ntr <= Pm.K_TRIALS, f"live reached L_{{Q,p_rand}} child {_cq} failed!"
        _ch_rels.append(len(_fc)); _ch_trials.append(_ntr)

    # [11] Live ACTUAL integer BOUND_CELL (+6 ball-mul margin) check on 12 exact-DP sampled primitive cells
    rng = random.Random(20261122 + seed)
    q_bc, s_bc, t_bc = _sample_ideal(90, "g", rng)
    B_bc, det_e, det_qcm = _lll_16d_qcm(_sq_hnf(q_bc, s_bc, t_bc))
    assert det_e == q_bc * q_bc and det_qcm == (1265625 ** 2) * q_bc * q_bc
    bc_gmin = 1e9
    bc_gmax = -1e9
    for _ in range(12):
        k = _sample_prim_ball_dp(dp_rd, dp_pref, 32, rng)
        v = [sum(k[i] * B_bc[i][j] for i in range(16)) for j in range(16)]
        nf, ng = _norms_fg(v[:8], v[8:])
        assert ng % q_bc == 0 and nf > 0 and ng > 0
        nf_ub, ng_ub, nf_c, ng_c = BC.bound_cell_int_fg(v[:8], v[8:])
        gf = nf_ub - log2(nf)
        gg = ng_ub - log2(ng)
        assert gf >= 0.0 and gg >= 0.0, f"BOUND_CELL undershoot! gf={gf}, gg={gg}"
        assert nf_c >= ceil(log2(nf)) and ng_c >= ceil(log2(ng))
        bc_gmin = min(bc_gmin, gf, gg)
        bc_gmax = max(bc_gmax, gf, gg)
    assert bc_gmin >= 0.0 and bc_gmax <= 0.5343

    # [7e] Live encountered-child offset check (>= 4 distinct children across fresh 90-bit parents)
    cf_lim = max(x[1] for x in Pm.QBINS)
    cg_lim = max(x[2] for x in Pm.QBINS)
    live_enc_cf, live_enc_cg = [], []
    seen_live_ch = set()
    S16 = log2(DS.CELLS[16])
    probe_ps = (4111, 4201, 4231, 4261, 4391, 4441, 4561, 4591, 4621, 4651, 4751, 4801, 4831, 4861, 4951, 5011, 5101, 5281, 5431, 5521, 5581, 5641, 5701, 5791, 5821, 5851, 5881, 6091, 6121, 6151, 6211, 6271, 6301, 6361, 6421, 6451, 6481, 6571, 6661, 6691, 6781, 6841, 6871, 6961, 6991, 7057, 7151, 7211, 7321, 7351, 7411, 7451, 7541, 7561, 7591, 7621, 7681, 7741, 7841, 7951, 8011, 8101, 8161, 8191)
    for side_p in ("f", "g"):
        qp, sp, tp = _sample_ideal(90, side_p, rng)
        Bp, _, _ = _lll_16d_qcm(_sq_hnf(qp, sp, tp))
        got_this_p = 0
        for _ in range(120):
            if got_this_p >= 2:
                break
            k = _sample_prim_ball_dp(dp_rd, dp_pref, 32, rng)
            v = [sum(k[i] * Bp[i][j] for i in range(16)) for j in range(16)]
            nf, ng = _norms_fg(v[:8], v[8:])
            for c_sd, Nv in (("f", nf // qp if side_p == "f" else nf), ("g", ng // qp if side_p == "g" else ng)):
                if got_this_p >= 2:
                    break
                for pp in probe_ps:
                    if (c_sd, pp) in seen_live_ch or Nv % pp != 0 or pp % 15 != 1 or not _is_prime(pp):
                        continue
                    e15 = (pp - 1) // 15
                    g0 = next((pow(g, e15, pp) for g in range(2, 80) if pow(pow(g, e15, pp), 3, pp) != 1 and pow(pow(g, e15, pp), 5, pp) != 1), None)
                    if g0 is None:
                        continue
                    aa, bb = v[:8], v[8:]
                    found_rt = False
                    for kk in (1, 2, 4, 7, 8, 11, 13, 14):
                        ry = pow(g0, kk, pp)
                        av = sum(aa[i] * pow(ry, i, pp) for i in range(8)) % pp
                        bv = sum(bb[i] * pow(ry, i, pp) for i in range(8)) % pp
                        if bv == 0:
                            continue
                        rt = (av * pow(bv, -1, pp)) % pp
                        ww = pow(ry, 5, pp)
                        ok_r = ((pow(rt, 3, pp) - U_BLS * ww * ww) % pp == 0) if c_sd == "g" else ((pow(rt, 9, pp) + (ww + 1) * pow(rt, 6, pp) + pow(rt, 3, pp) - ww) % pp == 0)
                        if ok_r:
                            Bc, de_c, _ = _lll_16d_qcm(_sq_hnf(pp, ry, rt))
                            assert de_c == pp * pp
                            cc = (log2(pp) + S16) / 16.0
                            for _ in range(4):
                                kc = _sample_prim_ball_dp(dp_rd, dp_pref, 16, rng)
                                vc = [sum(kc[i] * Bc[i][j] for i in range(16)) for j in range(16)]
                                nfc, ngc = _norms_fg(vc[:8], vc[8:])
                                live_enc_cf.append(log2(nfc) - 72.0 * cc)
                                live_enc_cg.append(log2(ngc) - 24.0 * cc - LOG2U8)
                            seen_live_ch.add((c_sd, pp))
                            got_this_p += 1
                            found_rt = True
                            break
                    if found_rt and got_this_p >= 2:
                        break
    assert len(seen_live_ch) >= 4 and len(live_enc_cf) >= 16 and max(live_enc_cf) <= cf_lim and max(live_enc_cg) <= cg_lim

    # Assert bundled v13 Exact 2-State Wald-Boole + Self-Normalized H1 Paley-Zygmund archive:
    assert set(c7e.keys()) == {"description", "prod_perbin_v13", "prod_reached_v13", "deg32_24_v13", "boot_r2_32x48"}

    # Campaign 1: 243,200 production deg-(72, 24) lattices (400 ideals/bin x 38 bins x 16 p_rand)
    pb13 = c7e["prod_perbin_v13"]
    pb_p = pb13["WALD_K24_G90_KAP200_MU30"]
    det_pow_16 = 1.0 - (Pm.G0 ** 16)
    unpooled_cp95_g90 = (1.0 - (0.05 ** (1.0 / 400.0))) / det_pow_16
    pooled_cp95_g90 = (1.0 - (0.05 ** (1.0 / 15200.0))) / det_pow_16
    assert unpooled_cp95_g90 < Pm.P_ADVERSE and pooled_cp95_g90 < Pm.P_ADVERSE / 40.0
    assert pb13["pool_f9_14_16_exact_size"] == 4740 and pb13["n_bins"] == 38
    assert pb13["geometry_global"]["min_euc_b0_over_M16"] > 0.50 and pb13["geometry_global"]["min_qcm_b0_over_detG8_M16"] > 0.50
    assert pb_p["n_ideals_total"] == 15200 and pb_p["n_ideals_per_bin"] == 400 and pb_p["m_draws"] == 16
    assert pb_p["total_lattices"] == 243200 and pb_p["total_cells"] == 7782400
    assert pb_p["ideals_with_any_fail_mu0"] == 0 and pb_p["draws_below_mu0"] == 0 and pb_p["max_bin_ideals_fail_mu0"] == 0
    assert pb_p["cell_vq_failures"] == 0 and pb_p["cell_zeros"] == 0 and pb_p["cell_capped"] == 0 and pb_p["vq_basis_failures"] == 0
    assert pb_p["all_lattices_Lambda_over_mu0"]["min"] >= 6.0 and pb_p["per_ideal_min16_over_muplan"]["min"] >= 3.0
    assert len(pb_p["per_bin"]) == 38
    for _bi, _pb in enumerate(pb_p["per_bin"]):
        _blo = 68.5 if _bi == 0 else float(_pb["bin"][0])
        assert _pb["n"] == 400 and _pb["lattices"] == 6400
        assert _pb["Rf"] == DS.PLAN16[_blo, "f"][4] and _pb["Rg"] == DS.PLAN16[_blo, "g"][4]
        assert _pb["fail_mu0"] == 0 and _pb["fail_muplan"] == 0 and _pb["min_over_mu0"] >= 6.0

    # Campaign 2: Genuine unplanted reached production >=68.5b children (hop-1 + hop-2) x 32 p_rand + controls
    pr13 = c7e["prod_reached_v13"]
    pr_p = pr13["WALD_K24_G90_KAP200_MU30"]
    _n_rch = int(pr13["n_reached_total"])
    _m_rch = int(pr_p["m_draws"])
    det_pow_32 = 1.0 - (Pm.G0 ** _m_rch)
    rch_cp95_g90 = (1.0 - (0.05 ** (1.0 / _n_rch))) / det_pow_32
    assert _n_rch >= 450 and _m_rch == 32 and pr_p["total_child_lattices"] == _n_rch * _m_rch
    assert pr_p["children_vq_ok"] == _n_rch and pr_p["children_any_fail_mu0"] == 0 and pr_p["draws_below_mu0"] == 0
    assert pr_p["control_any_fail_mu0"] == 0 and pr_p["control_draws_below_mu0"] == 0
    assert rch_cp95_g90 < Pm.P_ADVERSE
    assert pr_p["reached_min32_over_mu0"]["min"] >= 4.0 and pr_p["reached_min32_over_muplan"]["min"] >= 2.0

    # Campaign 3: Dimensionless H1 moment verification (kappa_cell) + Paley-Zygmund + live v13 trees in deg-(32, 24)
    d32 = c7e["deg32_24_v13"]
    _cA = d32["cohortA_parents_H1_moments"]
    _cB = d32["cohortB_reached_children_H1_moments"]
    _cC = d32["cohortC_live_v13_prand_retry_trees"]
    assert _cA["n_ideals"] == 128 and _cA["draws_per_ideal"] == 32 and _cA["total_lattices"] == 4096
    assert _cA["pooled_dimensionless_H1_moments"]["kappa_cell_self_norm_E_M_M1_over_c1Lam2"] < Pm.KAPPA_CELL
    assert _cA["paley_zygmund_check"]["max_Q_p1_emp"] < Pm.P1_CAP and _cA["paley_zygmund_check"]["k16_complete_failures"] == 0
    assert _cB["n_ideals"] == 128 and _cB["draws_per_ideal"] == 32 and _cB["total_lattices"] == 4096
    assert _cB["pooled_dimensionless_H1_moments"]["kappa_cell_self_norm_E_M_M1_over_c1Lam2"] < Pm.KAPPA_CELL
    assert _cB["paley_zygmund_check"]["max_Q_p1_emp"] < Pm.P1_CAP and _cB["paley_zygmund_check"]["k16_complete_failures"] == 0
    assert _cC["roots_tested"] == 128 and _cC["roots_completed"] == 128 and _cC["roots_failed"] == 0
    assert _cC["K_max"] == 16 and _cC["total_slots_visited"] == 1081 and _cC["max_trials_slot"] <= _cC["K_max"] < Pm.K_TRIALS

    # 100%-fresh production-HNF 16-cal / 48-holdout boot archive with full 64-entry per-target array
    br2 = c7e["boot_r2_32x48"]
    bnp = st["boot_norm_pool"]
    assert br2["fresh_from_production_hnf_all_64"] is True and bnp["fresh_from_production_hnf_all_64"] is True
    cg_by_t = br2["cg_max_by_target"]
    assert len(cg_by_t) == 64 and cg_by_t == bnp["cg_max_by_target"]
    assert abs(br2["cal16_cg_max"] - max(cg_by_t[:16])) < 1e-9 and abs(br2["holdout48_cg_max"] - max(cg_by_t[16:])) < 1e-9
    assert abs(br2["cg_max"] - max(br2["cal16_cg_max"], br2["holdout48_cg_max"])) < 1e-9
    assert br2["n_exceed_m_max"] == 0 and br2["exhaustive_cells_total"] == 64 * 17072 and br2["holdout48_cells"] == 48 * 17072
    assert br2["holdout48_exceed_targets"] == 0 and br2["holdout48_exceed_cal_cut"] == 0 and br2["holdout48_exceed_boot_bins_max"] == 0
    cp95_h48_live = 0.05 ** (1.0 / 48.0)
    assert cp95_h48_live >= Pm.BOOT_KEEP_FRAC and abs(br2["cp95_lower_bound_holdout48"] - round(cp95_h48_live, 5)) < 1e-5
    assert len(br2["raw_deciles_10"]) == 10 and br2["cg_max"] <= max(x[1] for x in Pm.BOOT_BINS)
    assert br2["E_S2_over_ES2"] <= 2.0 and Pm.UDUP <= br2["udup_48cells_mean_inv_m"] <= 1.0

    # Live multi-bin UDUP=0.40 exact-DP recomputation under WALD_K24_G90_KAP200_MU30
    _ud_live = {}
    for _uk, _side_u in (
        ("sq_104f_R48_48cells", "f"), ("sq_104g_R48_48cells", "g"),
        ("sq_134f_R56_48cells", "f"), ("sq_134g_R56_48cells", "g"),
        ("sq_204f_R12_48cells", "f"), ("sq_204g_R40_48cells", "g"),
    ):
        _ue = st["unit_orbit_udup"][_uk]
        assert len(_ue["m_true"]) == 48 and len(_ue["cfs"]) == 48 and len(_ue["cgs"]) == 48
        _inv_ms = [1.0 / m for m in _ue["m_true"]]
        _m_inv = sum(_inv_ms) / 48.0
        _se_m = math.sqrt(sum((x - _m_inv) ** 2 for x in _inv_ms) / (47.0 * 48.0))
        _b_lo = float(_ue["bin_lo"]); _b_hi = _b_lo + 5.0; _R_u = int(_ue["R"])
        assert f"_R{_R_u}_" in _uk and abs(_m_inv - _ue["mean_inv_m"]) < 1e-3 and _m_inv - 2.0 * _se_m >= Pm.UDUP
        _SR_u = log2(DS.CELLS[_R_u]); _lq_u = log2(int(_ue["q"]))
        _c_u = (_lq_u + Pm.BP_HI + _SR_u) / 16.0
        _c_dp = (_b_hi + Pm.BP_HI + _SR_u) / 16.0
        _sf_u = (_lq_u if _side_u == "f" else 0.0) + Pm.BP_LO; _sg_u = (_lq_u if _side_u == "g" else 0.0)
        _sf_dp = (_b_lo if _side_u == "f" else 0.0) + Pm.BP_LO; _sg_dp = (_b_lo if _side_u == "g" else 0.0)
        _n_s, _iA_s, _kA_s, _kB_s, _R_plan = DS.PLAN16[(_b_lo, _side_u)]
        _ws = [
            DS.ss_lower2(72.0 * _c_u + cf_ - _sf_u, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_F, 0.0)
            * DS.ss_lower2(24.0 * _c_u + Pm.LOG2U8 + cg_ - _sg_u, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_G, 0.0)
            for cf_, cg_ in zip(_ue["cfs"], _ue["cgs"])
        ]
        _pr_c = sum(
            w * DS.ss_lower2(72.0 * _c_dp + cf_ - _sf_dp, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_F, 0.0)
              * DS.ss_lower2(24.0 * _c_dp + Pm.LOG2U8 + cg_ - _sg_dp, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_G, 0.0)
            for w, cf_, cg_ in Pm.QBINS
        )
        _pr_sub = sum(
            DS.ss_lower2(72.0 * _c_dp + cf_ - _sf_dp, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_F, 0.0)
            * DS.ss_lower2(24.0 * _c_dp + Pm.LOG2U8 + cg_ - _sg_dp, Pm.LX, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_G, 0.0)
            for cf_, cg_ in Pm.QBINS_SUB
        ) / len(Pm.QBINS_SUB)
        _pr_dp = min(_pr_sub, Pm.PRCAP_NODE * _pr_c)
        if _R_u == _R_plan:
            assert abs(log2(_pr_dp) - res["DET"][(_b_lo, _b_hi), _side_u]["pr_log2"]) < 1e-12
        _tw = sum(w / m for w, m in zip(_ws, _ue["m_true"])) / sum(_ws)
        _dy = (sum(w / m for w, m in zip(_ws, _ue["m_true"])) / 48.0) / _pr_dp
        assert Pm.UDUP <= _m_inv <= 1.0 and 0.0 < _tw <= 1.0 and _dy >= Pm.UDUP
        _ud_live[_uk] = dict(bin_lo=int(_b_lo), R=_R_u, mean_inv_m=round(_m_inv, 4), true_yield_weighted_inv_m=round(_tw, 4), dedup_yield_over_charged_pr=round(_dy, 4))
    _ub = st["unit_orbit_udup"]["boot_R2_48cells"]
    assert len(_ub["m_true"]) == 48 and len(_ub["cgs"]) == 48
    _cb2 = (Pm.LOG2P4 + log2(DS.CELLS[2])) / 16.0
    _n1b, _iAb, _kAb, _ = DS.BOOT16
    _wsb = [DS.ss_lower2(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LX, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, 0.0) for cg_ in _ub["cgs"]]
    _prb = sum(w * DS.ss_lower2(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LX, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, 0.0) for w, cg_ in Pm.BOOT_BINS)
    _prb = min(sum(DS.ss_lower2(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LX, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, 0.0) for cg_ in Pm.BOOT_SUB) / len(Pm.BOOT_SUB), Pm.PRCAP_BOOT * _prb)
    _inv_mb = [1.0 / m for m in _ub["m_true"]]
    _mb_inv = sum(_inv_mb) / 48.0
    _se_mb = math.sqrt(sum((x - _mb_inv) ** 2 for x in _inv_mb) / (47.0 * 48.0))
    _twb = sum(w / m for w, m in zip(_wsb, _ub["m_true"])) / sum(_wsb)
    _dyb = (sum(w / m for w, m in zip(_wsb, _ub["m_true"])) / 48.0) / _prb
    assert abs(_mb_inv - _ub["mean_inv_m"]) < 1e-3
    assert Pm.UDUP <= _mb_inv - 2.0 * _se_mb and _mb_inv <= 1.0 and Pm.UDUP <= _twb <= 1.0 and _dyb >= Pm.UDUP
    _ud_live["boot_R2_48cells"] = dict(mean_inv_m=round(_mb_inv, 4), true_yield_weighted_inv_m=round(_twb, 4), dedup_yield_over_charged_pr=round(_dyb, 4))

    # Live deterministic rejection test on production-bin (135, g) short-subspace counterexample
    _q_ce = 81009504249198611324243595142290546015061
    _s_ce = 14886707496769063264106121185902690387739
    _t_ce = (5 * pow(8, -1, _q_ce)) % _q_ce
    _B_ce, _, _ = _lll_16d_qcm(_sq_hnf(_q_ce, _s_ce, _t_ce))
    assert BC.screen_special_q_basis(_B_ce, _q_ce, "g", _norms_fg) is False, "Must reject q^4 short-subspace ideal!"
    _k_ce = [1 if i == 0 else (47 if i == 8 else 0) for i in range(16)]
    _v_ce = [sum(_k_ce[i] * _B_ce[i][j] for i in range(16)) for j in range(16)]
    _d135 = res["DET"][(135.0, 140.0), "g"]
    _ok_ce, _nf_ce_rem, _ng_ce_rem, _, _ = BC.bound_cell_check_and_cap(_v_ce[:8], _v_ce[8:], _d135["nf"], _d135["ng"], lq=log2(_q_ce), side="g")
    assert _ok_ce is False and _nf_ce_rem > _d135["nf"] and _ng_ce_rem > _d135["ng"]
    assert BC.screen_special_q_basis(B_bc, q_bc, "g", _norms_fg) is True

    # [11b] Live q15_retain16 prime-splitting + non-vacuous p=2 (K=140) check + counted_fraction
    tab = st["prime_structure"]["table"]
    for pp in (2, 3, 5, 19, 31, 61, 101, 257):
        assert RT.structure(pp) == tab[str(pp)]
    assert abs(RT.counted_fraction() - st["omitted_law"]["counted_fraction"]) < 1e-6
    assert st["omitted_law"]["tail_j_ge_74_included"] is True
    rt_v = RT.verify_retain16(n_samples=64, seed=11, recompute_full=False)
    assert rt_v["probes_verified"] >= 10 and rt_v["sample_node_slack96"]["loss_f"] < 5e-5 and rt_v["sample_node_slack96"]["loss_g"] < 5e-5
    assert rt_v["residue_checks"][2]["K"] == 140 and rt_v["residue_checks"][2]["max_g"] <= 128.0

    # [12] Live q15_ecm16 L_min=2 branch-tree DP + geometric ladder + Suyama/Newton/singular-N-1 check
    assert EC.L_MIN_DEFAULT == 2 and EC.B_MIN_DEFAULT == 32 and EC.SEPARATION_DEFAULT == 0.5 and EC.RETRY_PR_DESC == 0.0043 and EC.NEWTON_PRECHECK_INSTR == 2.0 ** 26
    l2_boot, _ = EC.tree_branch_caps(int(ceil(bt["m"])), EC.B_MIN_DEFAULT, EC.L_MIN_DEFAULT)
    assert l2_boot >= ceil(bt["m"] / 64.0) ** 2
    sch_node = EC.schedule_v9(Pm.LX, int(ceil(res["DET"][(68.5, 75.0), "f"]["nf"])), Pm.ECM_SUCCESS, n_factors=50)
    sch_boot = EC.schedule_v9(bt["nA"], int(ceil(bt["m"])), Pm.ECM_SUCCESS, n_factors=bt["boot_n_factors"])
    lad_boot = EC.ladder_v9_instr(16.0, bt["nA"], int(ceil(bt["m"])), Pm.ECM_SUCCESS, bt["boot_n_factors"])
    assert sch_node["curves"] > 0 and sch_boot["curves"] > 0 and lad_boot >= sch_boot["instructions"]
    suyama_retry = EC.verify_suyama_retry_counterexample()
    _suyama_ok = bool(suyama_retry["perfect_power_precheck_ok"] and suyama_retry["bad_gcd_eq_N"] and suyama_retry["singular_nm1_rejected"] and suyama_retry["split_factor"] > 1)
    assert _suyama_ok, "Suyama retry & perfect-power/singular precheck failed!"

    wall_s = time.perf_counter() - t0
    return dict(
        all_ok=True,
        plan_tag="WALD_K24_G90_KAP200_MU30",
        wall_s=round(wall_s, 3),
        eval_384_s=round(dt_ev, 3),
        mesh=mesh,
        n_nodes=len(rows),
        n_sieve=len(_det_sieve),
        n_ecm=len(_det_ecm),
        max_R_used=max_R_used,
        lx=Pm.LX,
        bp=[Pm.BP_LO, Pm.BP_HI],
        p4740_exact_count=_cnt_f9_1416,
        p4740_split_primes=_n_split_1416,
        k_trials=Pm.K_TRIALS,
        g0=Pm.G0,
        kappa_cell=Pm.KAPPA_CELL,
        mu_0=Pm.MU_MIN,
        p1_pz_exact=round(p1_pz_exact, 6),
        p1_cap=Pm.P1_CAP,
        p_adverse=Pm.P_ADVERSE,
        two_log_units=round(two_log, 6),
        boot_u=round(bt["boot_u"], 4),
        sub_u=round(bt["sub_u"], 4),
        psi_boot=round(bt["psi_boot"], 6),
        e_boot_runs=round(bt["e_boot_runs"], 6),
        max_q_node=round(bt["max_q_node"], 6),
        max_psi_node=round(bt["max_psi_node"], 6),
        max_phi_good=float(f"{bt['max_phi_good']:.4e}"),
        expected_slot_trials_2log=round(EN, 2),
        boot_phi=0.0,
        boot_m=round(bt["m"], 2),
        boot_miss=round(bt["miss"], 5),
        boot_lin_keep=round(bt["lin_keep"], 5),
        max_miss_node=round(max_miss_node, 5),
        boot_keep_frac=bt["boot_keep_frac"],
        boot_n_factors=bt["boot_n_factors"],
        boot_n1=int(bt["n1"]),
        boot_nA=round(bt["nA"], 2),
        boot_kA=int(bt["kA"]),
        boot_R=int(bt["R"]),
        boot_one_log_log2=round(bt["one_log_log2"], 4),
        wald_stress=dict(
            p1_cap=Pm_st.P1_CAP,
            p_adverse=Pm_st.P_ADVERSE,
            psi_boot=round(bt_st["psi_boot"], 6),
            e_boot_runs=round(bt_st["e_boot_runs"], 6),
            max_q_node=round(bt_st["max_q_node"], 6),
            max_phi_good=float(f"{bt_st['max_phi_good']:.4e}"),
            two_log_units=round(res_st["two_log_units"], 6),
        ),
        uniform_vol_div_bound_38bins_cond_H1=dict(
            conditional_on_H1_chebotarev_bracket=True,
            bin_range=[68.5, 260.0],
            n_bins=38,
            worst_bin=_vol_worst_bin,
            max_uniform_density=float(f"{_vol_max_dens:.3e}"),
            log2_max_uniform_density=round(log2(_vol_max_dens), 2),
            strictly_below_2_m22=True,
        ),
        v13_empirical_summary=dict(
            perbin_total_lattices=pb_p["total_lattices"],
            perbin_total_cells=pb_p["total_cells"],
            perbin_ideals_fail_mu0=pb_p["ideals_with_any_fail_mu0"],
            perbin_draws_fail_mu0=pb_p["draws_below_mu0"],
            perbin_det_power_g90_m16=round(det_pow_16, 6),
            perbin_unpooled_cp95_g90=round(unpooled_cp95_g90, 6),
            perbin_pooled_cp95_g90=round(pooled_cp95_g90, 7),
            perbin_min_over_mu0=pb_p["all_lattices_Lambda_over_mu0"]["min"],
            perbin_median_over_mu0=pb_p["all_lattices_Lambda_over_mu0"]["median"],
            reached_n_total=_n_rch,
            reached_n_hop1=int(pr13["n_hop1_children"]),
            reached_n_hop2=int(pr13["n_hop2_grandchildren"]),
            reached_m_draws=_m_rch,
            reached_fails_mu0=pr_p["children_any_fail_mu0"],
            reached_det_power_g90_m32=round(det_pow_32, 6),
            reached_det_adj_cp95_g90=round(rch_cp95_g90, 6),
            reached_min32_over_mu0=pr_p["reached_min32_over_mu0"]["min"],
            reached_med32_over_mu0=pr_p["reached_min32_over_mu0"]["median"],
            deg32_parents_kappa_cell=_cA["pooled_dimensionless_H1_moments"]["kappa_cell_self_norm_E_M_M1_over_c1Lam2"],
            deg32_parents_max_p1=_cA["paley_zygmund_check"]["max_Q_p1_emp"],
            deg32_children_kappa_cell=_cB["pooled_dimensionless_H1_moments"]["kappa_cell_self_norm_E_M_M1_over_c1Lam2"],
            deg32_children_max_p1=_cB["paley_zygmund_check"]["max_Q_p1_emp"],
            deg32_trees_completed=_cC["roots_completed"],
            deg32_trees_k_max_archived=_cC["K_max"],
            deg32_trees_slots=_cC["total_slots_visited"],
            deg32_trees_mean_trials=_cC["mean_trials_per_slot"],
            deg32_trees_max_trials=_cC["max_trials_slot"],
            deg32_trees_max_depth=_cC["max_depth"],
        ),
        live_regen_prod_q685_240=_prod_regen,
        live_regen_strict_deg32_24=dict(
            B_fb=_Bfb_l,
            LP_max=_LPmax_l,
            pool_f4_size=len(_pool_f4_32),
            parents_ok=len(_par_rels),
            parent_rels_found=_par_rels,
            parent_trials_used=_par_trials,
            reached_children_ok=len(_ch_rels),
            child_rels_found=_ch_rels,
            child_trials_used=_ch_trials,
        ),
        pm_log2_rels=Pm.LOG2_RELS,
        pm_main_lb=Pm.MAIN_LB,
        pm_main_nf=round(Pm.MAIN_UF * Pm.MAIN_LB, 4),
        pm_main_ng=round(Pm.MAIN_UG * Pm.MAIN_LB, 4),
        pm_alpha_f=Pm.ALPHA_F,
        pm_alpha_g=Pm.ALPHA_G,
        pm_pair_setup=Pm.PAIR_SETUP,
        pm_lll_node=Pm.LLL_NODE,
        pm_lll_boot=Pm.LLL_BOOT,
        pm_bound_cell=Pm.BOUND_CELL,
        pm_sieve_hit=Pm.SIEVE_HIT,
        pm_ret_mult_f=Pm.RET_MULT_F,
        pm_ret_mult_g=Pm.RET_MULT_G,
        f4_counterexample_rejected=True,
        suyama_retry_ok=_suyama_ok,
        udup_live_recomputed=_ud_live,
        eager_all_100pct=True,
        linear_miss_retention=True,
        min_mu=round(min_mu, 2),
        min_ret=round(min_ret, 5),
        max_loss=float(f"{max_loss:.3e}"),
        bound_cell_live_min_gap_bits=bc_gmin,
        bound_cell_live_max_gap_bits=bc_gmax,
        live_enc_distinct_children=len(seen_live_ch),
        live_enc_max_cf_cg=[round(max(live_enc_cf), 2), round(max(live_enc_cg), 2)],
        retain16_live_ok=True,
        retain16_probes=int(rt_v["probes_verified"]),
        retain16_p2_K=int(rt_v["residue_checks"][2]["K"]),
        retain16_p2_max_g=rt_v["residue_checks"][2]["max_g"],
        boot_r2=br2,
        ecm_schedules=dict(
            l_min=EC.L_MIN_DEFAULT,
            b_min=EC.B_MIN_DEFAULT,
            sep=EC.SEPARATION_DEFAULT,
            p_r=EC.RETRY_PR_DESC,
            node_685f=dict(curves=sch_node["curves"], n_factors=sch_node["n_factors"], log2_ir=round(sch_node["log2_instructions"], 2)),
            boot_nA=dict(curves=sch_boot["curves"], n_factors=sch_boot["n_factors"], log2_ir=round(sch_boot["log2_instructions"], 2), ladder_log2_ir=round(log2(lad_boot), 2)),
        ),
    )


if __name__ == "__main__":
    out = run_all_checks()
    print(json.dumps(out, indent=2))
