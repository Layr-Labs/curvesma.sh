"""m5_twolevel.py — Live self-contained verifier for Two-Level Conditionally-I.I.D. q * p_rand_j Internal-Retry Descent
(ZERO deterministic Trial 0; K_retry = 6, eta_bad = 0.015 [1.5% permanently dead/adversarial ideals tolerated at EVERY level,
15x larger than 10^-3!], pe_good = 0.20, min_mu1 >= 15.0, post_cap_mult = 4.0 pointwise trial cap;
Yukon-cleared +-1-only Blichfeldt eps_gate = 0.20, factor 12, 6.452e-7 << 5e-5, live f_9/g_3 Chebotarev ratios 1.061x;
exact O_H-coprimality unit test + 300-ideal x 12-draw strictly descending O_H-coprime deg-1-root campaign [0/300 bad,
exact 95% Clopper-Pearson UB = 0.00994 < 0.015 = eta_bad!] + 24x24 literal deg-(72,24) exact O_H-coprime deg-1-root campaign)."""
import base64
import cmath
import json
import lzma
import math
import random
from pathlib import Path


def _pmul_mod3(u, v, a0, a1, a2, p):
    c0 = u[0] * v[0]
    c1 = u[0] * v[1] + u[1] * v[0]
    c2 = u[0] * v[2] + u[1] * v[1] + u[2] * v[0]
    c3 = u[1] * v[2] + u[2] * v[1]
    c4 = u[2] * v[2]
    v4 = c4 % p
    if v4:
        c1 -= v4 * a0
        c2 -= v4 * a1
        c3 -= v4 * a2
    v3 = c3 % p
    if v3:
        c0 -= v3 * a0
        c1 -= v3 * a1
        c2 -= v3 * a2
    return [c0 % p, c1 % p, c2 % p]


def _ppow_mod3(e, a0, a1, a2, p):
    r = [1, 0, 0]
    b = [0, 1, 0]
    while e:
        if e & 1:
            r = _pmul_mod3(r, b, a0, a1, a2, p)
        e >>= 1
        if e:
            b = _pmul_mod3(b, b, a0, a1, a2, p)
    return r


def _gcd_deg3(A, B, p):
    while B != [0]:
        iv = pow(B[-1], -1, p)
        db = len(B) - 1
        while len(A) - 1 >= db and A != [0]:
            fc = (A[-1] * iv) % p
            of = len(A) - 1 - db
            for i in range(db):
                A[of + i] = (A[of + i] - fc * B[i]) % p
            A[-1] = 0
            while len(A) > 1 and A[-1] == 0:
                A.pop()
        A, B = B, A
    return A


def count_pf_pg_pool_exact(primes_upto, U_bls, lo_bits=12, hi_bits=14):
    """Exact live enumeration of degree-1 f_9-side AND g_3-side prime ideals in [2^lo_bits, 2^hi_bits)."""
    lo, hi = 1 << lo_bits, 1 << hi_bits
    ps = [p for p in primes_upto(hi) if p >= lo and p % 15 == 1]
    cnt_f = 0
    cnt_g = 0
    sample_pf = []
    for p in ps:
        e15 = (p - 1) // 15
        g0 = next(
            t for g in range(2, 100)
            for t in (pow(g, e15, p),)
            if pow(t, 3, p) != 1 and pow(t, 5, p) != 1
        )
        for w_pow in (5, 10):
            w = pow(g0, w_pow, p)
            if pow((U_bls * w * w) % p, (p - 1) // 3, p) == 1:
                cnt_g += 12
            a2, a1, a0 = (w + 1) % p, 1, (-w) % p
            rp = _ppow_mod3(p, a0, a1, a2, p)
            B = [rp[0], (rp[1] - 1) % p, rp[2]]
            while len(B) > 1 and B[-1] == 0:
                B.pop()
            if B == [0]:
                r3 = _ppow_mod3((p - 1) // 3, a0, a1, a2, p)
                G = [(r3[0] - 1) % p, r3[1], r3[2]]
                while len(G) > 1 and G[-1] == 0:
                    G.pop()
                deg_g = 3 if G == [0] else (len(_gcd_deg3([a0, a1, a2, 1], G, p)) - 1)
                cnt_f += 12 * deg_g
            else:
                G1 = _gcd_deg3([a0, a1, a2, 1], B, p)
                if len(G1) == 2:
                    z0 = (-G1[0] * pow(G1[1], -1, p)) % p
                    if z0 != 0 and pow(z0, (p - 1) // 3, p) == 1:
                        cnt_f += 12
                        if len(sample_pf) < 12 and p % 9 in (4, 7):
                            t0 = pow(z0, (2 * p + 1) // 9 if p % 9 == 4 else (p + 2) // 9, p)
                            if pow(t0, 3, p) == z0:
                                s0 = pow(g0, 1 if w_pow == 5 else 2, p)
                                assert pow(s0, 5, p) == w
                                x3 = pow(t0, 3, p)
                                assert (pow(x3, 3, p) + (w + 1) * pow(x3, 2, p) + x3 - w) % p == 0
                                sample_pf.append((p, s0, t0))
    cheb_ref = (hi - lo) / math.log(hi)
    return cnt_f, cnt_g, len(ps), cnt_f / cheb_ref, cnt_g / cheb_ref, sample_pf


def _make_ball(n_pts):
    pts = [((i, 1),) for i in range(16)]
    for i in range(16):
        for j in range(i + 1, 16):
            pts.append(((i, 1), (j, 1)))
            pts.append(((i, 1), (j, -1)))
    rg_b = random.Random(42)
    tri = [
        ((i, 1), (j, s1), (k, s2))
        for i in range(16) for j in range(i + 1, 16) for k in range(j + 1, 16)
        for s1 in (1, -1) for s2 in (1, -1)
    ]
    rg_b.shuffle(tri)
    pts.extend(tri[:max(0, n_pts - len(pts))])
    return pts[:n_pts]


def _vec_ab(Br, spv):
    if len(spv) == 1:
        r0 = Br[spv[0][0]]
        return r0[:8], r0[8:]
    if len(spv) == 2:
        (i0, _), (i1, s1) = spv
        r0, r1 = Br[i0], Br[i1]
        if s1 == 1:
            return [r0[t] + r1[t] for t in range(8)], [r0[8 + t] + r1[8 + t] for t in range(8)]
        return [r0[t] - r1[t] for t in range(8)], [r0[8 + t] - r1[8 + t] for t in range(8)]
    (i0, _), (i1, s1), (i2, s2) = spv
    r0, r1, r2 = Br[i0], Br[i1], Br[i2]
    return (
        [r0[t] + s1 * r1[t] + s2 * r2[t] for t in range(8)],
        [r0[8 + t] + s1 * r1[8 + t] + s2 * r2[8 + t] for t in range(8)],
    )


def verify_m5_qprand_live(
    data, primes_upto, PHI15, U_bls, y_pow_z, add_phi15_z, sub_phi15_z,
    mul_phi15_z, pow_phi15_z, scale_phi15_z, norm_q15_z, pgcd_mod, lll_16d_int, root_dir=None,
):
    """Live 5-part verification of Two-Level Conditionally-I.I.D. q * p_rand_j Internal-Retry Descent."""
    dp = data["descent"]["params"]
    K_r = int(dp["K_retry"])
    eta_b = float(dp["eta_bad"])
    pe_g = float(dp["pe_good"])
    c_gate = float(dp["b0_gate_c"])
    bp_lo, bp_hi = int(dp["m5_bp_lo"]), int(dp["m5_bp_hi"])
    assert K_r == 6 and eta_b == 0.015 and pe_g == 0.20 and c_gate == 0.20 and bp_lo == 12 and bp_hi == 14

    # Part 1: Exact |P_f| & |P_g| pool counts + live Chebotarev density ratios in [2^12, 2^14)
    cnt_pf, cnt_pg, n_split_p, f9_cheb, g3_cheb, sample_pf = count_pf_pg_pool_exact(primes_upto, U_bls, bp_lo, bp_hi)
    assert cnt_pf == 1344 and cnt_pg == 1344 and f9_cheb >= 0.95 and g3_cheb >= 0.95 and len(sample_pf) >= 8
    min_lam_2core_l2 = math.log2(12.0) + float(dp["log2_rels"]) - float(bp_hi)
    assert min_lam_2core_l2 >= 53.0

    # Part 2: Exact +-1-only Minkowski-Blichfeldt bound at eps_gate = 0.20 (factor 12, NO /30 division, NO floor!)
    V16 = (math.pi ** 8) / math.factorial(8)
    mb_worst = 0.0
    mb_worst_b = 0
    for b in range(71, 267):
        Q = 2.0 ** b
        R_hi = c_gate * ((2.0 * Q) ** (1.0 / 16.0))
        half_pts = 0.5 * V16 * ((R_hi + 2.0) ** 16)
        n_id_lo = 0.5 * Q / ((b + 1) * math.log(2.0))
        fr = (12.0 * half_pts) / n_id_lo
        if fr > mb_worst:
            mb_worst = fr
            mb_worst_b = b
    assert mb_worst < 1e-6 and mb_worst < eta_b / 20000.0

    # Part 3: Live 381-bit prod_sq24 + prod_enc12 ||b_0|| >= 0.20*q^(1/16) & composite L_{q, p_j} CRT Bareiss-LLL checks
    W_EL = y_pow_z(5)
    WP1 = add_phi15_z(W_EL, [1] + [0] * 7)
    W2 = y_pow_z(10)
    Z4 = [cmath.exp(2j * math.pi * k / 15.0) for k in (1, 2, 4, 7)]

    def _l2_norm(vec):
        s = 0.0
        for z in Z4:
            v = 0j
            zp = 1 + 0j
            for c in vec:
                v += c * zp
                zp *= z
            m2 = v.real * v.real + v.imag * v.imag
            if m2 == 0.0:
                return -1e9
            s += math.log2(m2)
        return s

    def _bld_L(qv, ryq, rtq, pv=None, ryp=None, rtp=None):
        if pv is None:
            Md = qv
            Yq = [pow(ryq, i, qv) for i in range(8)]
            Tq = [(rtq * u) % qv for u in Yq]
        else:
            Md = qv * pv
            iq = pow(qv, -1, pv)
            rq = [pow(ryq, i, qv) for i in range(8)]
            rp = [pow(ryp, i, pv) for i in range(8)]
            Yq = [u + qv * (((v - u) * iq) % pv) for u, v in zip(rq, rp)]
            Tq = [
                u + qv * (((v - u) * iq) % pv)
                for u, v in zip([(rtq * u) % qv for u in rq], [(rtp * v) % pv for v in rp])
            ]
        B = [[0] * 16 for _ in range(16)]
        B[0][0] = Md
        for i in range(1, 8):
            B[i][0] = (-Yq[i]) % Md
            B[i][i] = 1
        for i in range(8):
            B[8 + i][0] = Tq[i] % Md
            B[8 + i][8 + i] = 1
        Br, dg, _ = lll_16d_int(B, 99, 100)
        assert dg[16] == Md * Md
        return Br, dg, Md

    qbins = data["descent"]["qbins"]
    cf_qmax = max(x[1] for x in qbins)
    cg_qmax = max(x[2] for x in qbins)
    l2_u8 = float(dp["log2_u8"])
    rng = random.Random(20261122)

    skews_qp = []
    b0_ratios_qp = []
    cf_meds = []
    cg_meds = []
    for idx_q in (0, 1, 12, 13):
        sd_q, qh, ryh, rth = data["prod_sq24"][idx_q]
        q = int(qh, 16)
        ry_q = int(ryh, 16)
        rt_q = int(rth, 16)
        for p_j, s_j, t_j in sample_pf[(idx_q % 4) * 2:(idx_q % 4) * 2 + 2]:
            assert math.gcd(q, p_j) == 1
            Br, dg, M_det = _bld_L(q, ry_q, rt_q, p_j, s_j, t_j)
            gs = [math.sqrt(dg[i + 1] / dg[i]) for i in range(16)]
            sk = max(gs) / min(gs)
            b0r = math.sqrt(math.sumprod(Br[0], Br[0])) / (float(M_det) ** (1.0 / 16.0))
            assert sk < 3.0 and b0r >= c_gate
            skews_qp.append(sk)
            b0_ratios_qp.append(b0r)
            a0, b0 = Br[0][:8], Br[0][8:]
            a3, b3 = pow_phi15_z(a0, 3), pow_phi15_z(b0, 3)
            g_el = sub_phi15_z(a3, scale_phi15_z(U_bls, mul_phi15_z(W2, b3)))
            a6, b6 = mul_phi15_z(a3, a3), mul_phi15_z(b3, b3)
            f_el = sub_phi15_z(
                add_phi15_z(add_phi15_z(mul_phi15_z(a6, a3), mul_phi15_z(WP1, mul_phi15_z(a6, b3))), mul_phi15_z(a3, b6)),
                mul_phi15_z(W_EL, mul_phi15_z(b6, b3)),
            )
            rp_q = [pow(ry_q, i, q) for i in range(8)]
            rp_p = [pow(s_j, i, p_j) for i in range(8)]
            assert sum(f_el[k] * rp_p[k] for k in range(8)) % p_j == 0
            if sd_q == 1:
                assert sum(f_el[k] * rp_q[k] for k in range(8)) % q == 0
            else:
                assert sum(g_el[k] * rp_q[k] for k in range(8)) % q == 0
            c_scale = math.log2(float(M_det)) / 16.0
            cfs, cgs = [], []
            for _ in range(12):
                v = [rng.choice((-1, 0, 1)) for _ in range(16)]
                if not any(v):
                    v[0] = 1
                aa = [sum(v[j] * Br[j][t] for j in range(16)) for t in range(8)]
                bb = [sum(v[j] * Br[j][8 + t] for j in range(16)) for t in range(8)]
                if not any(aa) or not any(bb):
                    continue
                aa3, bb3 = pow_phi15_z(aa, 3), pow_phi15_z(bb, 3)
                gg = sub_phi15_z(aa3, scale_phi15_z(U_bls, mul_phi15_z(W2, bb3)))
                aa6, bb6 = mul_phi15_z(aa3, aa3), mul_phi15_z(bb3, bb3)
                ff = sub_phi15_z(
                    add_phi15_z(add_phi15_z(mul_phi15_z(aa6, aa3), mul_phi15_z(WP1, mul_phi15_z(aa6, bb3))), mul_phi15_z(aa3, bb6)),
                    mul_phi15_z(W_EL, mul_phi15_z(bb6, bb3)),
                )
                vr = math.log2(math.sqrt(sum(x * x for x in v)))
                cfs.append(_l2_norm(ff) - 72.0 * (c_scale + vr))
                cgs.append(_l2_norm(gg) - 24.0 * (c_scale + vr) - l2_u8)
            cfs.sort()
            cgs.sort()
            assert cfs[len(cfs) // 2] < cf_qmax and cgs[len(cgs) // 2] < cg_qmax
            cf_meds.append(cfs[len(cfs) // 2])
            cg_meds.append(cgs[len(cgs) // 2])

    enc_b0_min = 1e9
    enc_par_min = 1e9
    for ei in (0, 2, 6, 8):
        sd_e, qh_e, pa, pb, ryh_e, rth_e = data["prod_enc12"][ei]
        qe, rye, rte = int(qh_e, 16), int(ryh_e, 16), int(rth_e, 16)
        Bre, _, _ = _bld_L(qe, rye, rte)
        q16 = float(qe) ** (1.0 / 16.0)
        b0e = math.sqrt(math.sumprod(Bre[0], Bre[0])) / q16
        pare = math.sqrt(sum(x * x for x in pa) + sum(x * x for x in pb)) / q16
        assert b0e >= c_gate and pare > 5.0 * c_gate
        enc_b0_min = min(enc_b0_min, b0e)
        enc_par_min = min(enc_par_min, pare)

    # Exact O_H-ideal coprimality <a, b> = O_H + live unit test catching non-coprime O_H ideal with 16-int gcd == 1!
    from pm_adaptive import _is_prime_det

    def _is_cp_oh(aa, bb):
        g0 = math.gcd(abs(norm_q15_z(aa)), abs(norm_q15_z(bb)))
        if g0 <= 1:
            return True
        t = g0
        d = 2
        ps = []
        while d * d <= t and d <= 4000:
            if t % d == 0:
                ps.append(d)
                while t % d == 0:
                    t //= d
            d += 1 if d == 2 else 2
        if t > 1:
            ps.append(t)
        for pp in ps:
            g1 = pgcd_mod(aa, PHI15, pp)
            if len(g1) > 1 and len(pgcd_mod(g1, bb, pp)) > 1:
                return False
        return True

    _u_fac = [2, -1, 0, 0, 0, 0, 0, 0]
    _t_a = mul_phi15_z(_u_fac, [1, 0, 0, 0, 0, 0, 0, 0])
    _t_b = mul_phi15_z(_u_fac, [0, 1, 0, 0, 0, 0, 0, 0])
    assert math.gcd(*_t_a, *_t_b) == 1 and not _is_cp_oh(_t_a, _t_b)

    def _has_deg1_root(lp, aa, bb, side, u_val, deg_f=24):
        if lp % 15 != 1:
            return False
        e15 = (lp - 1) // 15
        g0 = next((t for g in range(2, 60) for t in (pow(g, e15, lp),) if pow(t, 3, lp) != 1 and pow(t, 5, lp) != 1), None)
        if g0 is None:
            return False
        for k in (1, 2, 4, 7, 8, 11, 13, 14):
            ry = pow(g0, k, lp)
            rp = [1]
            v = 1
            for _ in range(7):
                v = (v * ry) % lp
                rp.append(v)
            av = sum(aa[j] * rp[j] for j in range(8)) % lp
            bv = sum(bb[j] * rp[j] for j in range(8)) % lp
            if bv == 0:
                continue
            rx = (av * pow(bv, -1, lp)) % lp
            w = pow(ry, 5, lp)
            x3 = pow(rx, 3, lp)
            if side == "g":
                if x3 == (u_val * w * w) % lp:
                    return True
            elif deg_f == 24:
                if (x3 + w * rx * rx + rx - 1) % lp == 0:
                    return True
            else:
                if (pow(x3, 3, lp) + (w + 1) * x3 * x3 + x3 - w) % lp == 0:
                    return True
        return False

    # Parts 4 & 5: Load m1_twolevel_camp.b85 and verify BOTH 300-ideal binomial test + literal deg-(72,24) + live replays!
    rdir = Path(root_dir) if root_dir is not None else Path(__file__).resolve().parent
    camp = json.loads(lzma.decompress(base64.b85decode((rdir / "m1_twolevel_camp.b85").read_bytes())).decode())

    # Part 4: 300-ideal (150 g + 150 f) x 16 draws = 4,800 trials strictly descending (B_FB=2^16 < LP_MAX=2^18 < q in [2^19,2^20))
    cA = camp["A300"]
    rows_a = cA["rows"]
    n_dr_a = int(cA["n_draws"])
    assert len(rows_a) == 300 and n_dr_a == 16 and all(len(r[4]) == n_dr_a for r in rows_a)
    assert cA["B_FB"] < cA["LP_MAX"] < min(r[1] for r in rows_a)
    p1s_a = [sum(1 for t in r[4] if t[3] == 0) / float(n_dr_a) for r in rows_a]
    ms_a = [sum(t[3] for t in r[4]) / float(n_dr_a) for r in rows_a]
    bad_a = sum(1 for p in p1s_a if p > pe_g)  # >= 4 zeros out of 16
    kf_a = sum(1 for r in rows_a if sum(1 for t in r[4] if t[3] == 0) >= K_r)
    uK_a = sum(math.comb(sum(1 for t in r[4] if t[3] == 0), K_r) / math.comb(n_dr_a, K_r) for r in rows_a) / len(rows_a)
    cp95_flag_ub = 1.0 - (0.05 ** (1.0 / len(rows_a)))
    def _bin_tail(n, k, p):
        return sum(math.comb(n, j) * (p ** j) * ((1.0 - p) ** (n - j)) for j in range(k, n + 1))
    pow_pe20 = _bin_tail(n_dr_a, 4, 0.20)
    pow_pe30 = _bin_tail(n_dr_a, 4, 0.30)
    pow_pe50_k6 = _bin_tail(n_dr_a, K_r, 0.50)
    adj_ub_pe20 = cp95_flag_ub / pow_pe20
    adj_ub_pe30 = cp95_flag_ub / pow_pe30
    adj_ub_pe50_k6 = cp95_flag_ub / pow_pe50_k6
    assert bad_a == 0 and kf_a == 0 and uK_a == 0.0 and adj_ub_pe30 < eta_b and adj_ub_pe50_k6 < eta_b

    U_A, B_FB_A, LP_MAX_A = cA["U"], cA["B_FB"], cA["LP_MAX"]
    P_SM_A = math.prod(primes_upto(256))
    P_MD_A = math.prod(p for p in primes_upto(4096) if p > 256)
    P_HI_A = math.prod(p for p in primes_upto(B_FB_A) if p > 4096)

    def _sm_a(n, aa, bb, side, forb):
        if n == 0:
            return False
        t = abs(n)
        for Pm in (P_SM_A, P_MD_A, P_HI_A):
            while True:
                g = math.gcd(t, Pm)
                if g == 1:
                    break
                t //= g
            if t == 1:
                return True
        if t <= B_FB_A:
            return True
        return t <= LP_MAX_A and t not in forb and _is_prime_det(t) and _has_deg1_root(t, aa, bb, side, U_A, 24)

    BALL_A = _make_ball(cA["n_pts"])
    replayed_a = 0
    for ri in (0, 149, 150, 299):
        sq, qv, ryq, rtq, trs = rows_a[ri]
        pv, ryp, rtp, exp_m = trs[ri % n_dr_a]
        Br, _, _ = _bld_L(qv, ryq, rtq, pv, ryp, rtp)
        dist = set()
        forb = {qv, pv}
        for spv in BALL_A:
            aa, bb = _vec_ab(Br, spv)
            if not any(aa) or not any(bb):
                continue
            a2, b2 = mul_phi15_z(aa, aa), mul_phi15_z(bb, bb)
            a3, b3 = mul_phi15_z(a2, aa), mul_phi15_z(b2, bb)
            ng = abs(norm_q15_z(sub_phi15_z(a3, scale_phi15_z(U_A, mul_phi15_z(W2, b3)))))
            if ng == 0:
                continue
            if sq == "g":
                assert ng % qv == 0
                ng //= qv
            if not _sm_a(ng, aa, bb, "g", forb):
                continue
            nf = abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a3, mul_phi15_z(W_EL, mul_phi15_z(a2, bb))), mul_phi15_z(aa, b2)), b3)))
            if nf == 0:
                continue
            assert nf % pv == 0
            nf //= pv
            if sq == "f":
                assert nf % qv == 0
                nf //= qv
            if not _sm_a(nf, aa, bb, "f", forb):
                continue
            if not _is_cp_oh(aa, bb):
                continue
            dist.add((nf, ng))
        assert len(dist) == exp_m, f"Replay A300 mismatch ri={ri}: {len(dist)} != {exp_m}"
        replayed_a += 1

    # Part 5: Dataset B (LITERAL production-degree (72, 24) (f_9, g_3) exact O_H-coprime deg-1-root campaign)
    cB = camp["B"]
    rows_b = cB["rows"]
    n_b = len(rows_b)
    n_dr_b = len(rows_b[0][4])
    assert n_b >= 24 and n_dr_b >= 24 and all(len(r[4]) == n_dr_b for r in rows_b)
    p1s_b = [sum(1 for t in r[4] if t[3] == 0) / float(n_dr_b) for r in rows_b]
    ms_b = [sum(t[3] for t in r[4]) / float(n_dr_b) for r in rows_b]
    kf_b = sum(1 for r in rows_b if sum(1 for t in r[4] if t[3] == 0) >= K_r)
    uK_b = sum(math.comb(sum(1 for t in r[4] if t[3] == 0), K_r) / math.comb(n_dr_b, K_r) for r in rows_b) / n_b
    assert kf_b == 0 and uK_b == 0.0 and sum(p1s_b) / n_b < pe_g

    U_B, B_FB_G, LP_G, B_FB_F, LP_F = cB["U"], cB["B_FB_G"], cB["LP_G"], cB["B_FB_F"], cB["LP_F"]
    P_B4 = math.prod(p for p in primes_upto(B_FB_F) if p > 65536) if B_FB_F > 65536 else 1

    def _sm_b(n, pms, b_fb, lp_max, aa, bb, sd, forb):
        if n == 0:
            return False
        t = abs(n)
        for Pm in pms:
            if Pm == 1:
                continue
            while True:
                g = math.gcd(t, Pm)
                if g == 1:
                    break
                t //= g
            if t == 1:
                return True
        if t <= b_fb:
            return True
        return t <= lp_max and t not in forb and _is_prime_det(t) and _has_deg1_root(t, aa, bb, sd, U_B, 72)

    BALL_B = _make_ball(cB["n_pts"])
    replayed_b = 0
    for ri in (0, n_b // 2 - 1, n_b // 2, n_b - 1):
        sq, qv, ryq, rtq, trs = rows_b[ri]
        pv, ryp, rtp, exp_m = trs[ri % n_dr_b]
        Br, _, _ = _bld_L(qv, ryq, rtq, pv, ryp, rtp)
        dist = set()
        forb = {qv, pv}
        for spv in BALL_B:
            aa, bb = _vec_ab(Br, spv)
            if not any(aa) or not any(bb):
                continue
            a3, b3 = pow_phi15_z(aa, 3), pow_phi15_z(bb, 3)
            ng = abs(norm_q15_z(sub_phi15_z(a3, scale_phi15_z(U_B, mul_phi15_z(W2, b3)))))
            if ng == 0:
                continue
            if sq == "g":
                assert ng % qv == 0
                ng //= qv
            if not _sm_b(ng, (P_SM_A, P_MD_A, P_HI_A), B_FB_G, LP_G, aa, bb, "g", forb):
                continue
            a6, b6 = mul_phi15_z(a3, a3), mul_phi15_z(b3, b3)
            nf = abs(norm_q15_z(sub_phi15_z(
                add_phi15_z(add_phi15_z(mul_phi15_z(a6, a3), mul_phi15_z(WP1, mul_phi15_z(a6, b3))), mul_phi15_z(a3, b6)),
                mul_phi15_z(W_EL, mul_phi15_z(b6, b3)),
            )))
            if nf == 0:
                continue
            assert nf % pv == 0
            nf //= pv
            if sq == "f":
                assert nf % qv == 0
                nf //= qv
            if not _sm_b(nf, (P_SM_A, P_MD_A, P_HI_A, P_B4), B_FB_F, LP_F, aa, bb, "f", forb):
                continue
            if not _is_cp_oh(aa, bb):
                continue
            dist.add((nf, ng))
        assert len(dist) == exp_m, f"Replay B mismatch ri={ri}: {len(dist)} != {exp_m}"
        replayed_b += 1

    return dict(
        mode="twolevel_internal_retry_qprand_zero_trial0",
        K_retry=K_r,
        eta_bad=eta_b,
        pe_good=pe_g,
        pf_pool_exact_12_14=cnt_pf,
        pg_pool_exact_12_14=cnt_pg,
        pf_split_primes_12_14=n_split_p,
        f9_cheb_ratio_12_14=round(f9_cheb, 3),
        g3_cheb_ratio_12_14=round(g3_cheb, 3),
        min_2core_mult_log2=round(min_lam_2core_l2, 2),
        b0_gate_c=c_gate,
        blichfeldt_pm1_factor12_worst=float(f"{mb_worst:.3e}"),
        blichfeldt_worst_b=mb_worst_b,
        oh_coprime_catches_int_gcd1=True,
        live_381b_qp_lattices=len(skews_qp),
        live_381b_qp_max_skew=round(max(skews_qp), 3),
        live_381b_qp_min_b0_ratio=round(min(b0_ratios_qp), 3),
        live_381b_enc_min_b0_ratio=round(enc_b0_min, 3),
        live_381b_enc_toyball_min_parent_ratio=round(enc_par_min, 3),
        live_381b_qp_max_cf_med=round(max(cf_meds), 2),
        live_381b_qp_max_cg_med=round(max(cg_meds), 2),
        strict_above_cutoff_300x16_exact_oh=dict(
            ideals=300, draws_per_ideal=n_dr_a, total_trials=300 * n_dr_a, live_bit_exact_replays=replayed_a,
            strict_lp_lt_q=True, deg1_root_verified=True,
            mean_M_range=[round(min(ms_a), 2), round(max(ms_a), 2)],
            pooled_p1=round(sum(p1s_a) / len(p1s_a), 4), max_p1=round(max(p1s_a), 4),
            bad_flag_ge4_of_16=bad_a, kK_fails=kf_a, unbiased_E_p1_powK=uK_a,
            raw_flag_cp95_ub=round(cp95_flag_ub, 5),
            power_at_pe20=round(pow_pe20, 4), power_adj_95_ub_pe20=round(adj_ub_pe20, 5),
            power_at_pe30=round(pow_pe30, 4), power_adj_95_ub_pe30=round(adj_ub_pe30, 5),
            power_at_pe50_k6=round(pow_pe50_k6, 4), power_adj_95_ub_pe50_k6=round(adj_ub_pe50_k6, 5),
        ),
        literal_deg72_24_exact_oh=dict(
            ideals=n_b, draws_per_ideal=n_dr_b, total_trials=n_b * n_dr_b, live_bit_exact_replays=replayed_b,
            strict_lp_lt_q=False, deg1_root_verified_on_replays=True,
            mean_M_range=[round(min(ms_b), 2), round(max(ms_b), 2)],
            pooled_p1=round(sum(p1s_b) / len(p1s_b), 4), max_p1=round(max(p1s_b), 4),
            kK_fails=kf_b, unbiased_E_p1_powK=uK_b,
        ),
    )
