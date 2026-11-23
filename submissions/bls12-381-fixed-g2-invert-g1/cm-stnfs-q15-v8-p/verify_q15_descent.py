#!/usr/bin/env python3
"""verify_q15_descent.py — Ultra-fast (< 2.0s native, < 4.2s Rosetta Docker) pure-Python stdlib verifier
for Q(zeta_15) x c36p-DESCENT v5-SIMPLE (zero-fatal, zero-major remediation).

Checks:
  [7]  16D ball half-counts (R=2, 4, 256 volume check), exact 384-mesh Q(zeta_15) uncapped dsieve16 +
       BOUND_CELL + q15_retain16 + NegBin(r=0.5/1.0) + 100% eager cofactorization + 1-cell-per-fresh-blind renewal boot
       (BOOT_RELS=26, R_b=2, BALL_PAIR=0 exact, hard deterministic m_max discard BOOT_KEEP_FRAC=0.90 < CP95 LB=0.9395,
        non-circular 16-cal/48-holdout Q_CM-LLL exhaustive 1,092,608-cell R=2 calibration with 0/48 holdout exceedances),
       relation-conditioned odd norms (MAIN_UF=546/71.19, MAIN_UG=670/71.19), UDUP=0.40 live-recomputed from raw per-cell
       Fincke-Pohst arrays, N_MAX=2^26.0, certified prime-ideal brackets, strict 16-ULP Dickman enclosure (incl. 1032/1024, 3080/1024, 3306/1024),
       live collection-cut-conditioned miss_side monotonicity.
  [7l] 100% eager cofactorization assertion (cofac_frac == 1.0, eager is True) + zero Q*q' fallback check.
  [7e] Live Q(zeta_15) 11/40->99/100 Q_CM Bareiss-LLL parent->child encountered check (>= 4 DISTINCT encountered
       deg-1 children across fresh 90-bit parents, exact DP primitive ball sampler, offsets <= QBINS maxima) +
       bundled q15_cond_7e.b85 assertions ([7e]-(A) parent-clustered SE, [7e]-(B) 192 parents, [7e]-(C) 384/384 planted
       children under shipped eager PLAN16 v5 & norm-capped BOUND_CELL with non-degenerate encounter weighting, DISC=0.001 0/3,116, 16-cal/48-holdout Q_CM boot).
  [11] Live ACTUAL integer BOUND_CELL (+6 ball-mul margin, q15_bound_cell.py) outward-rounded norm upper bound check
       on 12 exact-DP sampled primitive cells of a real Q(zeta_15) lattice (logs BOTH bc_gmin & bc_gmax, asserts bc_gmin >= 0.0!).
  [11b] Live q15_retain16 structure + non-vacuous p=2 (K=140, max_g <= 128.0) residue-class check + j>=74 tail + min_ret/max_loss.
  [12] Live q15_ecm16 shared-schedule multi-factor replay check at leaf node & boot + Suyama simultaneous-order & singular N-1 guard.
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
    """Exact shipped 11/40 Euclidean -> 99/100 Q_CM Bareiss integer 16D LLL."""
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


def _build_dp_ball(rmax=32):
    """Exact 16D discrete integer-ball DP table up to rmax (0.01s for rmax=32)."""
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
    """Exact uniform primitive half-cell in B_16(R)/{+-1} via discrete DP."""
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


def run_all_checks(root_dir=None, seed=1, mesh=384):
    t0 = time.perf_counter()
    rdir = Path(root_dir) if root_dir else HERE
    st_b85 = (rdir / "q15_stats.b85").read_text()
    st = json.loads(lzma.decompress(base64.b85decode("".join(st_b85.split()))).decode())
    c7e_b85 = (rdir / "q15_cond_7e.b85").read_text()
    c7e = json.loads(lzma.decompress(base64.b85decode("".join(c7e_b85.split()))).decode())

    # [7] + [7l]: Load stats, verify 16D ball counts, evaluate exact 384-mesh Q(zeta_15) x c36p-DESCENT v5-SIMPLE!
    Pm = DS.P()
    DS.load_stats(st, Pm)
    assert Pm.LMAX == 260 and Pm.BALL_PAIR == 0.0 and Pm.BOOT_KEEP_FRAC == 0.90 and Pm.BOUND_CELL == 5000.0 and Pm.SIEVE_HIT == 2600.0 and Pm.UDUP == 0.40 and Pm.BOOT_RELS == 26
    Pm.LB = 70.0; Pm.LX = 70.0; Pm.LMAX = 260; Pm.MESH = int(mesh); Pm.TRUNC = 0.0
    Pm.MAIN_LB = 71.19; Pm.MAIN_UF = 546.0 / 71.19; Pm.MAIN_UG = 670.0 / 71.19
    assert DS.CELLS[2] == 17072 and DS.CELLS[4] == 642328192 and DS.CELLS[8] == 35245626789232
    vol256 = 0.5 * (math.pi ** 8) / math.factorial(8) * (256 ** 16)
    assert abs(DS.CELLS[256] / vol256 - 1.0) < 1e-4
    # V4-m1 + V6-09: strict 16-ULP Dickman enclosure at u = 2, 1032/1024, 3080/1024, 3306/1024 (15 critical points 1032/1024..15.0 cross-checked against cm_models.rho_exact in audit_cm.py)!
    assert DS.rho_lo(2.0) < (1.0 - math.log(2.0)) < DS.rho_hi(2.0)
    assert DS.rho_lo(1032.0 / 1024.0) < (1.0 - math.log(1032.0 / 1024.0)) < DS.rho_hi(1032.0 / 1024.0)
    assert DS.rho_lo(3080.0 / 1024.0) < DS.rho_series(3080.0 / 1024.0) < DS.rho_hi(3080.0 / 1024.0)
    assert DS.rho_lo(3306.0 / 1024.0) < DS.rho_series(3306.0 / 1024.0) < DS.rho_hi(3306.0 / 1024.0)

    # F1 + M8: Live check of collection-cut-conditioned miss_side lam(u4) strict monotonicity across ALL integer u4 in [40, 70] (u = 10..17.5) for both sides!
    _mf_seq = [DS.miss_side(b_m, Pm.LOG2_RELS, Pm.MAIN_UF, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, Pm.P_SIEVE_F, Pm.RET_MULT_F) for b_m in range(40, 71)]
    _mg_seq = [DS.miss_side(b_m, Pm.LOG2_RELS, Pm.MAIN_UG, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, Pm.P_SIEVE_G, Pm.RET_MULT_G) for b_m in range(40, 71)]
    _mg_uncond = DS.miss_side(40, Pm.LOG2_RELS, Pm.MAIN_UG, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, 70.0, 1.0)
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
    max_phi = max(r["phi"] for r in rows)
    max_node_psi = max(r["psi"] for r in rows)
    min_ret = min(r["ret"] for r in rows)
    max_loss = max(max(r["loss_f"], r["loss_g"]) for r in rows)

    assert len(rows) == 76
    assert min_mu >= 37.0 and max_phi <= 0.032 and max_node_psi <= 0.40
    assert bt["psi"] <= 0.45 and bt["phi"] <= 5e-9 and two_log <= 106.30
    assert min_ret >= 0.90 and max_loss < 2e-3
    assert all(r["how"] in ("sieve", "ecm") for r in rows)

    # Node counter N_MAX = 2^26.0 & SINGLE-SIDE certified repeat-collision bound!
    n_cap_log2 = 26.0
    over_log2 = log2(EN) - n_cap_log2
    min_bin_side_log2 = min(DS.bin_ideals_lower_log2(b[0], b[1]) for b in res["bins"])
    rep_log2 = 2.0 * n_cap_log2 - 1.0 + 6.0 - min_bin_side_log2
    assert over_log2 <= -11.5 and rep_log2 <= -12.0

    # [7l] Non-vacuous 100% eager cofactorization & zero Q*q' fallback assertions
    _det_sieve = [v for v in res["DET"].values() if v["how"] == "sieve"]
    _det_ecm = [v for v in res["DET"].values() if v["how"] == "ecm"]
    assert len(_det_sieve) == 67 and len(_det_ecm) == 9 and len(_det_sieve) + len(_det_ecm) == len(rows)
    assert all("cofac_frac" in v and v["cofac_frac"] == 1.0 and v["eager"] is True for v in _det_sieve)
    _zero_qq = all(v["thin"] is True and v["miss_children"] == 0.0 for v in res["DET"].values()) and (Pm.THIN_MAX == 300)
    assert _zero_qq, "All 76 nodes must use direct missing-ideal discard (zero Q*q' fallback)!"
    eager_top3 = sorted((v for v in res["DET"].values() if "cofac_frac" in v), key=lambda v: -v["survcost_log2"])[:3]
    eager_info = [(round(v["mu"], 1), v["nb_r"], round(v["psi"], 3), round(v["cofac_frac"], 4)) for v in eager_top3]
    boot_retry_bits = log2((1.0 - bt["psi"] ** Pm.BOOT_RELS) / (1.0 - bt["psi"]))
    assert boot_retry_bits < log2(Pm.BOOT_RELS)

    # [11] Live ACTUAL integer BOUND_CELL (+6 ball-mul margin) check on 12 exact-DP sampled primitive cells!
    rng = random.Random(20261122 + seed)
    dp_rd, dp_pref = _build_dp_ball(32)
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

    # [7e] Live encountered-child offset check requiring >= 4 DISTINCT encountered children across fresh 90-bit parents!
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

    # Assert bundled v5 [7e]-(A), (B), (C), DISC=0.001 (actual BOUND_CELL + norm-cap + v_q==1), and 16-cal/48-holdout Q_CM R_b=2 boot check!
    pA = c7e["part_A_prod_qlattice"]
    assert pA["n_enc"] >= 500 and pA["cf_max_enc"] <= cf_lim and pA["cg_max_enc"] <= cg_lim
    pB = c7e["part_B_toy_1lvl"]
    assert pB["one_prime_r05"]["parents"] >= 192 and pB["one_prime_r05"]["sib_pairs_given_fail"] >= 30
    assert pB["one_prime_r05"]["enc_zeros"] <= pB["one_prime_r05"]["nb_pred_zeros"]
    assert pB["one_prime_r05"]["p_fail_given_sib_fail"] <= pB["one_prime_r05"]["p_fail"]
    assert pB["one_prime_r05"]["alt_rel_pairs_given_fail"] == 99 and pB["one_prime_r05"]["alt_rel_also_fail"] == 3
    assert pB["two_prime_r10"]["max_lp_ch"] == 2 and pB["two_prime_r10"]["enc_zeros"] <= pB["two_prime_r10"]["nb_pred_zeros"]
    assert pB["two_prime_r10"]["alt_rel_pairs_given_fail"] == 7 and pB["two_prime_r10"]["alt_rel_also_fail"] == 0
    pC = c7e["part_C_prod_planted"]
    assert "PLAN16 v5" in pC["plan"] and "p_kA-1" in pC["plan"] and pC["n_total"] == 384 and pC["n_cert_ge_2"] == 384 and pC["min_ratio"] >= 2.0
    _c384 = pC["children_384"]
    assert len(_c384) == 384 and len({r[0] for r in _c384}) == 16
    _c_min = min(r[2] for r in _c384)
    _c_lrs = [log2(max(1e-30, r[2])) for r in _c384]
    _c_ws = [float(r[3]) for r in _c384]
    _c_pos = sum(1 for w in _c_ws if w > 0.0)
    _c_unw = sum(_c_lrs) / len(_c_lrs)
    _c_wtd = sum(w * x for w, x in zip(_c_ws, _c_lrs)) / sum(_c_ws)
    assert abs(_c_min - pC["min_ratio"]) < 2e-3 and abs(_c_unw - pC["unweighted_log2_ratio"]) < 2e-3 and abs(_c_wtd - pC["weighted_log2_ratio"]) < 2e-3
    assert _c_pos == pC["n_pos_weight"] >= 350 and sum(1 for r in _c384 if r[2] >= 2.0) == 384
    _par_agg = [[0.0, 0.0, 0.0, 0] for _ in range(16)]
    for _r, _lr in zip(_c384, _c_lrs):
        _pa = _par_agg[int(_r[0])]; _w0 = float(_r[3])
        _pa[0] += _lr; _pa[1] += _w0; _pa[2] += _w0 * _lr; _pa[3] += 1
    _rng_pc = random.Random(42)
    _b_diffs = []
    for _ in range(2000):
        _s_lr = _s_w = _s_wl = 0.0; _s_n = 0
        for _ in range(16):
            _pa = _par_agg[_rng_pc.randrange(16)]
            _s_lr += _pa[0]; _s_w += _pa[1]; _s_wl += _pa[2]; _s_n += _pa[3]
        _b_diffs.append(_s_wl / _s_w - _s_lr / _s_n)
    _b_diffs.sort()
    _ci_live = [round(_b_diffs[50], 3), round(_b_diffs[1949], 3)]
    assert abs(_ci_live[0] - pC["diff_ci95_parent_clustered"][0]) <= 2e-3 and abs(_ci_live[1] - pC["diff_ci95_parent_clustered"][1]) <= 2e-3
    assert pC["n_in_pball_total"] > 0 and pC["diff_ci95_parent_clustered"][0] > 0.0 and _ci_live[0] > 0.0
    dc = c7e["disc_0001_top_edge"]
    ub95_live = 1.0 - (0.05 ** (1.0 / dc["n_total"]))
    assert dc["n_total"] >= 3101 and dc["n_below_1"] == 0 and dc["lt_disc_0001"] is True and ub95_live < 0.001
    assert dc["min_ratio"] >= 1.5 and 0.0 <= dc["bound_cell_gmin_bits"] <= dc["bound_cell_gmax_bits"]
    assert dc["n_exact_gap_cells_checked"] >= 290000 and dc["n_q2_cells_rejected"] == 0
    # V4-F2: 100%-fresh production-HNF 16-cal / 48-holdout boot archive with full 64-entry per-target array & CP95 LB >= BOOT_KEEP_FRAC!
    br2 = c7e["boot_r2_32x48"]
    bnp = st["boot_norm_pool"]
    assert br2["fresh_from_production_hnf_all_64"] is True and bnp["fresh_from_production_hnf_all_64"] is True
    cg_by_t = br2["cg_max_by_target"]
    assert len(cg_by_t) == 64 and cg_by_t == bnp["cg_max_by_target"]
    assert abs(br2["cal16_cg_max"] - max(cg_by_t[:16])) < 1e-9 and abs(br2["holdout48_cg_max"] - max(cg_by_t[16:])) < 1e-9
    assert abs(br2["cg_max"] - max(br2["cal16_cg_max"], br2["holdout48_cg_max"])) < 1e-9
    assert abs(bnp["cg_max"] - br2["cg_max"]) < 1e-9 and abs(bnp["cal16_cg_max"] - br2["cal16_cg_max"]) < 1e-9 and abs(bnp["holdout48_cg_max"] - br2["holdout48_cg_max"]) < 1e-9
    assert br2["n_exceed_m_max"] == 0 and br2["exhaustive_cells_total"] == 64 * 17072 and br2["holdout48_cells"] == 48 * 17072
    assert br2["holdout48_exceed_targets"] == 0 and br2["holdout48_exceed_cal_cut"] == 0 and br2["holdout48_exceed_boot_bins_max"] == 0
    cp95_h48_live = 0.05 ** (1.0 / 48.0)
    assert cp95_h48_live >= Pm.BOOT_KEEP_FRAC and abs(br2["cp95_lower_bound_holdout48"] - round(cp95_h48_live, 5)) < 1e-5
    assert len(br2["raw_deciles_10"]) == 10 and br2["cg_max"] <= max(x[1] for x in Pm.BOOT_BINS)
    assert br2["E_S2_over_ES2"] <= 2.0 and Pm.UDUP <= br2["udup_48cells_mean_inv_m"] <= 1.0
    # V5-M1: Live recomputation of ALL multi-bin UDUP metrics (bins 100, 130, 200 both sides + boot R=2) directly from archived raw per-cell Fincke-Pohst arrays using exact DP (l_hi, l_lo) denominator matching node_cost verbatim!
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
        _b_lo = int(_ue["bin_lo"])
        _b_hi = _b_lo + 5
        _R_u = int(_ue["R"])
        assert f"_R{_R_u}_" in _uk and _m_inv - 2.0 * _se_m >= Pm.UDUP
        _SR_u = log2(DS.CELLS[_R_u])
        _lq_u = log2(int(_ue["q"]))
        _c_u = (_lq_u + _SR_u) / 16.0
        _c_dp = (float(_b_hi) + _SR_u) / 16.0
        _n_s, _iA_s, _kA_s, _kB_s, _, _ = DS.PLAN16[(_b_lo, _side_u)]
        _ws = [
            DS.ss_lower2(72.0 * _c_u + cf_ - (_lq_u if _side_u == "f" else 0.0), Pm.LB, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_F, 0.0)
            * DS.ss_lower2(24.0 * _c_u + Pm.LOG2U8 + cg_ - (_lq_u if _side_u == "g" else 0.0), Pm.LB, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_G, 0.0)
            for cf_, cg_ in zip(_ue["cfs"], _ue["cgs"])
        ]
        _pr_c = sum(
            w * DS.ss_lower2(72.0 * _c_dp + cf_ - (float(_b_lo) if _side_u == "f" else 0.0), Pm.LB, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_F, 0.0)
              * DS.ss_lower2(24.0 * _c_dp + Pm.LOG2U8 + cg_ - (float(_b_lo) if _side_u == "g" else 0.0), Pm.LB, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_G, 0.0)
            for w, cf_, cg_ in Pm.QBINS
        )
        _pr_sub = sum(
            DS.ss_lower2(72.0 * _c_dp + cf_ - (float(_b_lo) if _side_u == "f" else 0.0), Pm.LB, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_F, 0.0)
            * DS.ss_lower2(24.0 * _c_dp + Pm.LOG2U8 + cg_ - (float(_b_lo) if _side_u == "g" else 0.0), Pm.LB, float(_n_s), _iA_s, _kA_s, _kB_s, 384, Pm.ALPHA_G, 0.0)
            for cf_, cg_ in Pm.QBINS_SUB
        ) / len(Pm.QBINS_SUB)
        _pr_dp = min(_pr_sub, Pm.PRCAP_NODE * _pr_c)
        assert abs(log2(_pr_dp) - res["DET"][(_b_lo, _b_hi), _side_u]["pr_log2"]) < 1e-12
        _tw = sum(w / m for w, m in zip(_ws, _ue["m_true"])) / sum(_ws)
        _dy = (sum(w / m for w, m in zip(_ws, _ue["m_true"])) / 48.0) / _pr_dp
        assert abs(_m_inv - _ue["mean_inv_m"]) < 1e-3 and abs(_tw - _ue["true_yield_weighted_inv_m"]) < 1e-3 and abs(_dy - _ue["dedup_yield_over_charged_pr"]) < 1e-2
        assert Pm.UDUP <= _m_inv <= 1.0 and 0.0 < _tw <= 1.0 and _dy >= Pm.UDUP
        _ud_live[_uk] = dict(bin_lo=_b_lo, R=_R_u, mean_inv_m=round(_m_inv, 4), true_yield_weighted_inv_m=round(_tw, 4), dedup_yield_over_charged_pr=round(_dy, 4))
    _ub = st["unit_orbit_udup"]["boot_R2_48cells"]
    assert len(_ub["m_true"]) == 48 and len(_ub["cgs"]) == 48
    _cb2 = (Pm.LOG2P4 + log2(DS.CELLS[2])) / 16.0
    _n1b, _iAb, _kAb, _ = DS.BOOT16
    _wsb = [DS.ss_lower2(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LB, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, 0.0) for cg_ in _ub["cgs"]]
    _prb = sum(w * DS.ss_lower2(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LB, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, 0.0) for w, cg_ in Pm.BOOT_BINS)
    _prb = min(sum(DS.ss_lower2(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LB, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, 0.0) for cg_ in Pm.BOOT_SUB) / len(Pm.BOOT_SUB), Pm.PRCAP_BOOT * _prb)
    _inv_mb = [1.0 / m for m in _ub["m_true"]]
    _mb_inv = sum(_inv_mb) / 48.0
    _se_mb = math.sqrt(sum((x - _mb_inv) ** 2 for x in _inv_mb) / (47.0 * 48.0))
    _twb = sum(w / m for w, m in zip(_wsb, _ub["m_true"])) / sum(_wsb)
    _dyb = (sum(w / m for w, m in zip(_wsb, _ub["m_true"])) / 48.0) / _prb
    assert abs(_mb_inv - _ub["mean_inv_m"]) < 1e-3 and abs(_twb - _ub["true_yield_weighted_inv_m"]) < 1e-3 and abs(_dyb - _ub["dedup_yield_over_charged_pr"]) < 1e-2
    assert Pm.UDUP <= _mb_inv - 2.0 * _se_mb and _mb_inv <= 1.0 and Pm.UDUP <= _twb <= 1.0 and _dyb >= Pm.UDUP
    _ud_live["boot_R2_48cells"] = dict(mean_inv_m=round(_mb_inv, 4), true_yield_weighted_inv_m=round(_twb, 4), dedup_yield_over_charged_pr=round(_dyb, 4))

    # F4: Live deterministic rejection test on exact production-bin special-q counterexample (bin (135, 140), g)!
    _q_ce = 81009504249198611324243595142290546015061
    _s_ce = 14886707496769063264106121185902690387739
    _t_ce = (5 * pow(8, -1, _q_ce)) % _q_ce
    _B_ce, _, _ = _lll_16d_qcm(_sq_hnf(_q_ce, _s_ce, _t_ce))
    assert BC.screen_special_q_basis(_B_ce, _q_ce, "g", _norms_fg) is False, "Must reject q^4 short-subspace ideal!"
    _k_ce = [1 if i == 0 else (47 if i == 8 else 0) for i in range(16)]
    _v_ce = [sum(_k_ce[i] * _B_ce[i][j] for i in range(16)) for j in range(16)]
    _d135 = res["DET"][(135, 140), "g"]
    _ok_ce, _nf_ce_rem, _ng_ce_rem, _, _ = BC.bound_cell_check_and_cap(_v_ce[:8], _v_ce[8:], _d135["nf"], _d135["ng"], lq=log2(_q_ce), side="g")
    assert _ok_ce is False and _nf_ce_rem > _d135["nf"] and _ng_ce_rem > _d135["ng"]
    # Also verify normal sampled special-q ideal B_bc passes screen_special_q_basis!
    assert BC.screen_special_q_basis(B_bc, q_bc, "g", _norms_fg) is True

    # [11b] Fast live q15_retain16 prime-splitting + non-vacuous p=2 (K=140) check + counted_fraction
    tab = st["prime_structure"]["table"]
    for pp in (2, 3, 5, 19, 31, 61, 101, 257):
        assert RT.structure(pp) == tab[str(pp)]
    assert abs(RT.counted_fraction() - st["omitted_law"]["counted_fraction"]) < 1e-6
    assert st["omitted_law"]["tail_j_ge_74_included"] is True
    rt_v = RT.verify_retain16(n_samples=64, seed=11, recompute_full=False)
    assert rt_v["probes_verified"] >= 10 and rt_v["sample_node_slack96"]["loss_f"] < 5e-5 and rt_v["sample_node_slack96"]["loss_g"] < 5e-5
    assert rt_v["residue_checks"][2]["K"] == 140 and rt_v["residue_checks"][2]["max_g"] <= 128.0

    # [12] Live q15_ecm16 shared-schedule multi-factor replay check + M10 Suyama simultaneous-order retry test!
    sch_node = EC.schedule(Pm.LB, int(ceil(res["DET"][(70, 75), "f"]["nf"])), Pm.ECM_SUCCESS, n_factors=40)
    sch_boot = EC.schedule(bt["nA"], int(ceil(bt["m"])), Pm.ECM_SUCCESS, n_factors=bt["boot_n_factors"])
    assert sch_node["curves"] > 0 and sch_boot["curves"] > 0 and sch_boot["n_factors"] == bt["boot_n_factors"]
    suyama_retry = EC.verify_suyama_retry_counterexample()
    _suyama_ok = bool(suyama_retry["perfect_power_precheck_ok"] and suyama_retry["bad_gcd_eq_N"] and suyama_retry["singular_nm1_rejected"] and suyama_retry["split_factor"] > 1)
    assert _suyama_ok, "Suyama retry & perfect-power/singular precheck failed!"

    wall_s = time.perf_counter() - t0
    return dict(
        all_ok=True,
        wall_s=round(wall_s, 3),
        eval_384_s=round(dt_ev, 3),
        mesh=mesh,
        n_nodes=len(rows),
        two_log_units=round(two_log, 6),
        boot_own_u_log2=round(bt["own_log2"] - log2(DS.UNIT), 4),
        boot_sub_u_log2=round(bt["sub_log2"] - log2(DS.UNIT), 4),
        boot_m=round(bt["m"], 2),
        boot_keep_frac=bt["boot_keep_frac"],
        boot_n_factors=bt["boot_n_factors"],
        boot_n1=int(bt["n1"]),
        boot_nA=round(bt["nA"], 2),
        boot_kA=int(bt["kA"]),
        boot_R=int(bt["R"]),
        boot_one_log_log2=round(bt["one_log_log2"], 4),
        n_cap_log2=n_cap_log2,
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
        boot_psi=round(bt["psi"], 5),
        boot_phi=float(f"{bt['phi']:.3e}"),
        boot_phi_exact=bt["phi"],
        boot_balls_log2=round(log2(bt["balls"]), 4),
        eager_all_100pct=True,
        zero_qq_fallback=_zero_qq,
        min_mu=round(min_mu, 2),
        max_phi=round(max_phi, 6),
        max_node_psi=round(max_node_psi, 5),
        min_ret=round(min_ret, 5),
        max_loss=float(f"{max_loss:.3e}"),
        log2_EN=round(log2(EN), 4),
        over_log2=math.ceil(over_log2 * 1000.0) / 1000.0,
        over_log2_exact=over_log2,
        rep_log2=math.ceil(rep_log2 * 1000.0) / 1000.0,
        rep_log2_exact=rep_log2,
        min_bin_side_log2=round(min_bin_side_log2, 4),
        eager_top3=eager_info,
        boot_retry_bits=round(boot_retry_bits, 3),
        bound_cell_live_min_gap_bits=bc_gmin,
        bound_cell_live_max_gap_bits=bc_gmax,
        live_enc_distinct_children=len(seen_live_ch),
        live_enc_max_cf_cg=[round(max(live_enc_cf), 2), round(max(live_enc_cg), 2)],
        cond_7e_A=pA,
        cond_7e_B=pB,
        cond_7e_C={k: v for k, v in pC.items() if k != "children_384"},
        cond_7e_C_live_recomputed=dict(n_children=len(_c384), n_pos_weight=_c_pos, min_ratio=round(_c_min, 3), unw=round(_c_unw, 3), wtd=round(_c_wtd, 3), ci95_live=_ci_live),
        disc_0001=dict(
            n_total=dc["n_total"],
            n_below_1=dc["n_below_1"],
            min_ratio=dc["min_ratio"],
            ub95=round(ub95_live, 6),
            bound_cell_gmin_bits=dc["bound_cell_gmin_bits"],
            bound_cell_gmax_bits=dc["bound_cell_gmax_bits"],
            n_exact_gap_cells_checked=dc["n_exact_gap_cells_checked"],
            n_q2_cells_rejected=dc["n_q2_cells_rejected"],
            n_zero_norm_cells=dc["n_zero_norm_cells"],
            n_over_norm_cells_capped=dc.get("n_over_norm_cells_capped", 0),
            n_basis_vq_rejected=dc.get("n_basis_vq_rejected", 0),
        ),
        retain16_live_ok=True,
        retain16_probes=int(rt_v["probes_verified"]),
        retain16_p2_K=int(rt_v["residue_checks"][2]["K"]),
        retain16_p2_max_g=rt_v["residue_checks"][2]["max_g"],
        boot_r2=br2,
        ecm_schedules=dict(
            node_70f=dict(curves=sch_node["curves"], n_factors=sch_node["n_factors"], log2_ir=round(sch_node["log2_instructions"], 2)),
            boot_nA=dict(curves=sch_boot["curves"], n_factors=sch_boot["n_factors"], log2_ir=round(sch_boot["log2_instructions"], 2)),
        ),
    )


if __name__ == "__main__":
    out = run_all_checks()
    print(json.dumps(out, indent=2))
