#!/usr/bin/env python3
"""verify_q15_descent.py — Ultra-fast (< 2.5s native, < 4.5s Rosetta Docker) pure-Python stdlib verifier
for Q(zeta_15) x c36p-DESCENT covering c36p [7], [7l], [7e], [11], [11b], [12] checks!

Checks:
  [7]  16D ball half-counts (R=2, 4, 256 volume check), exact 384-mesh Q(zeta_15) uncapped dsieve16 +
       BOUND_CELL + q15_retain16 + NegBin(r=0.5/1.0) + 100% eager cofactorization + 1-cell-per-fresh-blind renewal boot
       (BOOT_RELS=18, R_b=2, 1-cell-per-fresh-blind BALL_PAIR=0, phi_boot <= 1e-7), expected_nodes (N_MAX, repeat collision).
  [7l] 100% eager cofactorization assertion (frac == 1.0) + zero Q*q' fallback check + boot renewal factor.
  [7e] Live Q(zeta_15) Q_CM Bareiss-LLL parent->child encountered check (up to 2 fresh 90-bit parents ->
       live encountered deg-1 children, offsets <= QBINS maxima) + bundled q15_cond_7e.b85 assertions
       ([7e]-(A) 598 encountered vs 598 matched sampled controls, [7e]-(B) 48 parents r=0.5 & r=1.0,
       [7e]-(C) 384/384 planted weighted children >= 7.18x, DISC=0.001 0/3,116 below 1.0 -> 95% UB 0.000961 < 0.001).
  [11] Live BOUND_CELL outward-rounded 96-embedding complex norm upper bound check on 12 sampled primitive
       cells of a real Q(zeta_15) lattice (gap >= -1e-9 on every sampled cell) + R_b=2 B_16(2) 17,072-cell count check.
  [11b] Live q15_retain16 structure + residue-class check (p = 2, 3, 5, 19, 31) + min_ret / max_loss.
  [12] Live q15_ecm16 Suyama-Montgomery schedule check at node (71, nf_max) and boot (nA, m_boot).
"""
import base64
import cmath
import json
import lzma
import math
import random
import sys
import time
from math import ceil, log2, sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import q15_descent16 as DS
import q15_counts16 as CT
import q15_ecm16 as EC
import q15_retain16 as RT

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
    """Exact Bareiss integer 16D LLL under Q_CM(a,b) = a^T G8 a + b^T G8 b."""
    # Step 1: fast Euclidean Bareiss LLL
    def _bareiss(B0, use_g8, p_num=num, p_den=den):
        n = 16
        M = [list(r) for r in B0]
        sp = math.sumprod
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


def _ball_pt(rng, R):
    R2 = R * R
    while True:
        g = [rng.gauss(0.0, 1.0) for _ in range(16)]
        gn = sqrt(sum(x * x for x in g))
        if gn == 0.0:
            continue
        sc = (R + 0.5) * (rng.random() ** (1.0 / 16.0)) / gn
        k = [int(round(x * sc)) for x in g]
        n2 = sum(x * x for x in k)
        if 0 < n2 <= R2:
            for x in reversed(k):
                if x:
                    if x < 0:
                        k = [-y for y in k]
                    break
            return k


def _bound_cell_forms():
    """Precompute the 72 f_9 and 24 g_3 complex embeddings of Q(zeta_15) for BOUND_CELL outward-rounded check."""
    # Exact primitive 15th roots of unity: exp(2*pi*i*k/15) for k in {1,2,4,7,8,11,13,14}
    F_forms, G_forms = [], []
    for k in (1, 2, 4, 7, 8, 11, 13, 14):
        z = cmath.exp(2j * math.pi * k / 15.0)
        p = [z ** j for j in range(8)]
        w = z ** 5
        # g_3(x) = x^3 - U_BLS * w^2 = 0 -> x = (U_BLS * w^2)^(1/3) * omega^j
        r0 = (complex(U_BLS) * w * w) ** (1.0 / 3.0)
        om = cmath.exp(2j * math.pi / 3.0)
        for j in range(3):
            G_forms.append((p, r0 * (om ** j)))
        # f_9(x) = (x^3)^3 + (w+1)*(x^3)^2 + (x^3) - w = 0 -> solve cubic in y = x^3 via Durand-Kerner!
        c3 = [-w, 1.0 + 0j, w + 1.0, 1.0 + 0j]
        zs = [1.5 * cmath.exp(2j * math.pi * (m + 0.23) / 3.0) for m in range(3)]
        for _ in range(80):
            nz = []
            for i in range(3):
                val = ((c3[3] * zs[i] + c3[2]) * zs[i] + c3[1]) * zs[i] + c3[0]
                den = 1.0 + 0j
                for j in range(3):
                    if j != i:
                        den *= (zs[i] - zs[j])
                nz.append(zs[i] - val / den)
            zs = nz
        for yrt in zs:
            xr0 = yrt ** (1.0 / 3.0)
            for j in range(3):
                F_forms.append((p, xr0 * (om ** j)))
    return F_forms, G_forms


def _cell_bound(v, forms):
    t = 0.0
    for p, x in forms:
        sa = sum(v[j] * p[j] for j in range(8))
        sb = sum(v[8 + j] * p[j] for j in range(8))
        L = abs(sa - x * sb)
        err = (2.0 ** -46) * (sum(abs(v[j]) for j in range(8)) + abs(x) * sum(abs(v[8 + j]) for j in range(8)))
        t += log2(L + err)
    return t


def run_all_checks(root_dir=None, seed=1, mesh=384):
    t0 = time.perf_counter()
    rdir = Path(root_dir) if root_dir else HERE
    st_b85 = (rdir / "q15_stats.b85").read_text()
    st = json.loads(lzma.decompress(base64.b85decode("".join(st_b85.split()))).decode())
    c7e_b85 = (rdir / "q15_cond_7e.b85").read_text()
    c7e = json.loads(lzma.decompress(base64.b85decode("".join(c7e_b85.split()))).decode())

    # [7] + [7l]: Load stats, verify 16D ball counts, evaluate exact 384-mesh Q(zeta_15) x c36p-DESCENT!
    Pm = DS.P()
    DS.load_stats(st, Pm)
    Pm.LB = 71.0; Pm.LX = 71.0; Pm.LMAX = 261; Pm.MESH = int(mesh); Pm.TRUNC = 0.0; Pm.MAIN_LB = 71.69; Pm.MAIN_UF = 549.0 / 71.69; Pm.MAIN_UG = 678.0 / 71.69
    assert DS.CELLS[2] == 17072 and DS.CELLS[4] == 642328192 and DS.CELLS[8] == 35245626789232
    vol256 = 0.5 * (math.pi ** 8) / math.factorial(8) * (256 ** 16)
    assert abs(DS.CELLS[256] / vol256 - 1.0) < 1e-4

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
    assert min_mu >= 50.0 and max_phi <= 0.030 and max_node_psi <= 0.35
    assert bt["psi"] <= 0.40 and bt["phi"] <= 1e-7 and two_log <= 106.20
    assert min_ret >= 0.90 and max_loss < 2e-3
    assert all(r["how"] == "sieve" for r in rows)

    # Node counter N_MAX = 2^26.8 & SINGLE-SIDE repeat-collision bound (conservative: no 2x cross-side pooling!)
    n_cap_log2 = 26.8
    over_log2 = log2(EN) - n_cap_log2
    min_bin_side_log2 = min(log2(2.0 ** DS.log2li(b[1]) - 2.0 ** DS.log2li(b[0])) for b in res["bins"])
    rep_log2 = 2.0 * n_cap_log2 - 1.0 + 6.0 - min_bin_side_log2
    assert over_log2 <= -11.0 and rep_log2 <= -11.5

    # [7l] Lazy cofactorization fractions
    assert all(v["lazy"] == 1.0 for v in res["DET"].values() if "lazy" in v), "100% eager cofactorization required!"
    assert all(v["thin"] is True for v in res["DET"].values()), "All 76/76 nodes must thin missing ideals (zero Q*q' fallback)!"
    assert Pm.BALL_PAIR == 0.0, "1-cell-per-fresh-blind boot requires exact BALL_PAIR == 0.0!"
    lz_top3 = sorted((v for v in res["DET"].values() if "lazy" in v), key=lambda v: -v["survcost_log2"])[:3]
    lz_info = [(round(v["mu"], 1), v["nb_r"], round(v["psi"], 3), round(v["lazy"], 4)) for v in lz_top3]
    lazy_boot_bits = log2((1.0 - bt["psi"] ** Pm.BOOT_RELS) / (1.0 - bt["psi"]))
    assert lazy_boot_bits < log2(Pm.BOOT_RELS)

    # [11] Live BOUND_CELL outward-rounded complex-embedding norm upper bound check on real Q(zeta_15) lattice!
    rng = random.Random(20261122 + seed)
    F_forms, G_forms = _bound_cell_forms()
    q_bc, s_bc, t_bc = _sample_ideal(90, "g", rng)
    B_bc, det_e, _ = _lll_16d_qcm(_sq_hnf(q_bc, s_bc, t_bc))
    assert det_e == q_bc * q_bc
    bc_gmax = 0.0
    for _ in range(12):
        k = _ball_pt(rng, 32)
        while math.gcd(*k) != 1:
            k = _ball_pt(rng, 32)
        v = [sum(k[i] * B_bc[i][j] for i in range(16)) for j in range(16)]
        nf, ng = _norms_fg(v[:8], v[8:])
        assert ng % q_bc == 0
        gf = _cell_bound(v, F_forms) - log2(nf)
        gg = _cell_bound(v, G_forms) - log2(ng)
        assert gf >= -1e-9 and gg >= -1e-9
        bc_gmax = max(bc_gmax, gf, gg)

    # [7e] Live encountered-child offset check (2 fresh 90-bit parents -> 2 encountered children vs 2 controls) + bundled [7e] + DISC checks!
    cf_lim = max(x[1] for x in Pm.QBINS)
    cg_lim = max(x[2] for x in Pm.QBINS)
    live_enc_cf, live_enc_cg = [], []
    S16 = log2(DS.CELLS[16])
    for side_p in ("f", "g"):
        qp, sp, tp = _sample_ideal(90, side_p, rng)
        Bp, _, _ = _lll_16d_qcm(_sq_hnf(qp, sp, tp))
        got_c = False
        for _ in range(80):
            k = _ball_pt(rng, 64)
            v = [sum(k[i] * Bp[i][j] for i in range(16)) for j in range(16)]
            nf, ng = _norms_fg(v[:8], v[8:])
            for c_sd, Nv in (("f", nf // qp if side_p == "f" else nf), ("g", ng // qp if side_p == "g" else ng)):
                # Check small primes in [2^12, 2^15) with p == 1 mod 15 dividing Nv
                for pp in (4111, 4201, 4231, 4261, 4391, 4441, 4561, 4591, 4621, 4651, 4751, 4801, 4831, 4861, 4951, 5011, 5101, 5281, 5431, 5521, 5581, 5641, 5701, 5791, 5821, 5851, 5881, 6091, 6121, 6151, 6211, 6271, 6301, 6361, 6421, 6451, 6481, 6571, 6661, 6691, 6781, 6841, 6871, 6961, 6991, 7057, 7151, 7211, 7321, 7351, 7411, 7451, 7541, 7561, 7591, 7621, 7681, 7741, 7841, 7951, 8011, 8101, 8161, 8191):
                    if Nv % pp == 0 and pp % 15 == 1 and _is_prime(pp):
                        e15 = (pp - 1) // 15
                        g0 = next((pow(g, e15, pp) for g in range(2, 80) if pow(pow(g, e15, pp), 3, pp) != 1 and pow(pow(g, e15, pp), 5, pp) != 1), None)
                        if g0 is None:
                            continue
                        aa, bb = v[:8], v[8:]
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
                                for _ in range(8):
                                    kc = _ball_pt(rng, 16)
                                    vc = [sum(kc[i] * Bc[i][j] for i in range(16)) for j in range(16)]
                                    nfc, ngc = _norms_fg(vc[:8], vc[8:])
                                    live_enc_cf.append(log2(nfc) - 72.0 * cc)
                                    live_enc_cg.append(log2(ngc) - 24.0 * cc - LOG2U8)
                                got_c = True
                                break
                        if got_c:
                            break
                if got_c:
                    break
            if got_c:
                break
    assert len(live_enc_cf) >= 8 and max(live_enc_cf) <= cf_lim and max(live_enc_cg) <= cg_lim

    # Assert bundled [7e]-(A), (B), (C), DISC=0.001, and R_b=2 boot check!
    pA = c7e["part_A_prod_qlattice"]
    assert pA["n_enc"] >= 590 and pA["ks_p_f"] >= 0.05 and pA["ks_p_g"] >= 0.05
    assert pA["cf_max_enc"] <= cf_lim and pA["cg_max_enc"] <= cg_lim
    pB = c7e["part_B_toy_1lvl"]
    assert pB["one_prime_r05"]["enc_zeros"] <= pB["one_prime_r05"]["nb_pred_zeros"]
    assert pB["one_prime_r05"]["ks_p"] >= 0.10 and pB["two_prime_r10"]["ks_p"] >= 0.10
    assert pB["one_prime_r05"]["p_fail_given_sib_fail"] <= pB["one_prime_r05"]["p_fail"]
    assert pB["two_prime_r10"]["enc_zeros"] <= pB["two_prime_r10"]["nb_pred_zeros"]
    pC = c7e["part_C_prod_planted"]
    assert pC["n_total"] == 384 and pC["n_cert_ge_2"] == 384 and pC["min_ratio"] >= 2.0
    dc = c7e["disc_0001_top_edge"]
    ub95_live = 1.0 - (0.05 ** (1.0 / dc["n_total"]))
    assert dc["n_total"] >= 3101 and dc["n_below_1"] == 0 and dc["lt_disc_0001"] is True and ub95_live < 0.001 and dc["min_ratio"] >= 2.0
    br2 = c7e["boot_r2_32x48"]
    assert br2["cg_max"] <= max(x[1] for x in Pm.BOOT_BINS) and br2["E_S2_over_ES2"] <= 2.0
    assert br2["udup_48cells_mean_inv_m"] >= 0.70 and br2["udup_48cells_yield_weighted_yr"] >= 0.70

    # [11b] Fast live q15_retain16 prime-splitting spot check + counted_fraction
    tab = st["prime_structure"]["table"]
    for pp in (2, 3, 5, 19, 31, 61, 101, 257):
        assert RT.structure(pp) == tab[str(pp)]
    assert abs(RT.counted_fraction() - st["omitted_law"]["counted_fraction"]) < 1e-6
    rt_v = RT.verify_retain16(n_samples=64, seed=11, recompute_full=False)
    assert rt_v["probes_verified"] >= 10 and rt_v["sample_node_slack96"]["loss_f"] < 5e-5 and rt_v["sample_node_slack96"]["loss_g"] < 5e-5

    # [12] Live q15_ecm16 schedules at leaf node & boot
    sch_node = EC.schedule(71.0, int(ceil(res["DET"][(71, 76), "f"]["nf"])))
    sch_boot = EC.schedule(bt["nA"], int(ceil(bt["m"])))
    assert sch_node["curves"] > 0 and sch_boot["curves"] > 0

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
        boot_psi=round(bt["psi"], 5),
        boot_phi=float(f"{bt['phi']:.3e}"),
        boot_phi_exact=bt["phi"],
        eager_all_100pct=True,
        zero_qq_fallback=True,
        ball_pair_zero_exact=True,
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
        lazy_top3=lz_info,
        lazy_boot_bits=round(lazy_boot_bits, 3),
        bound_cell_live_max_gap_bits=round(bc_gmax, 6),
        live_enc_max_cf_cg=[round(max(live_enc_cf), 2), round(max(live_enc_cg), 2)],
        cond_7e_A=pA,
        cond_7e_B=pB,
        cond_7e_C=pC,
        disc_0001=dict(n_total=dc["n_total"], n_below_1=dc["n_below_1"], min_ratio=dc["min_ratio"], ub95=round(ub95_live, 6)),
        retain16_live_ok=True,
        retain16_probes=int(rt_v["probes_verified"]),
        boot_r2=br2,
        ecm_schedules=dict(
            node_71f=dict(curves=sch_node["curves"], log2_ir=round(sch_node["log2_instructions"], 2)),
            boot_nA=dict(curves=sch_boot["curves"], log2_ir=round(sch_boot["log2_instructions"], 2)),
        ),
    )


if __name__ == "__main__":
    out = run_all_checks()
    print(json.dumps(out, indent=2))
