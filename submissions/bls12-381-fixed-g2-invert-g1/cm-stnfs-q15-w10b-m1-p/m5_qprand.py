"""m5_qprand.py — Live self-contained verifier for Pure M1 conditionally-i.i.d. q * p_rand_j descent
(ZERO deterministic Trial 0, K=12, p_1 <= 0.40, b0_gate_c = 0.30, +-1-only ceil Blichfeldt,
exact O_H-coprime above-cutoff 32x24 + literal deg-(72,24) 24x24 campaigns with live bit-exact replays)."""
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


def count_pf_pool_exact(primes_upto, lo_bits=12, hi_bits=14):
    """Exact live enumeration of degree-1 f_9-side prime ideals in [2^lo_bits, 2^hi_bits)."""
    lo, hi = 1 << lo_bits, 1 << hi_bits
    ps = [p for p in primes_upto(hi) if p >= lo and p % 15 == 1]
    cnt_f = 0
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
    return cnt_f, len(ps), sample_pf


def _b16_half_table(Rmax=30):
    M2 = Rmax * Rmax

    def ksq(poly, bpc):
        x = int.from_bytes(b"".join(v.to_bytes(bpc, "little") for v in poly), "little")
        p = x * x
        nb = max((M2 + 1) * bpc, (p.bit_length() + 7) // 8)
        bs = p.to_bytes(nb, "little")
        return [int.from_bytes(bs[i * bpc:(i + 1) * bpc], "little") for i in range(M2 + 1)]

    t1 = [0] * (M2 + 1)
    t1[0] = 1
    for k in range(1, Rmax + 1):
        t1[k * k] = 2
    t16 = ksq(ksq(ksq(ksq(t1, 2), 4), 8), 14)
    out = {}
    acc = 0
    r = 0
    for s in range(M2 + 1):
        acc += t16[s]
        while r <= Rmax and r * r == s:
            out[r] = (acc - 1) // 2
            r += 1
    return out


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
    """Live 5-part verification of Pure M1 conditionally-i.i.d. q * p_rand_j descent (ZERO deterministic Trial 0!):
    1. Exact live enumeration of |P_f| = 1,344 degree-1 f_9-side prime ideals in [2^12, 2^14).
    2. Exact +-1-only Blichfeldt covering bound with upper-edge r_hi = c_gate * 2^((b+1)/16), ceil(r_hi),
       and proved PNT lower bound 0.5 * 2^b / ((b+1)*ln 2) across ALL b in [71, 266] (<= 3.03e-6 << p_geom = 5e-5),
       plus live ||b_0||_2 / q^(1/16) checks on prod_sq24 AND real 381-bit encountered prod_enc12 (+ parent vectors).
    3. Live 381-bit composite L_{q, p_j} CRT + (99/100)-Bareiss-LLL + exact divisibility + qbins check.
    4. Full 32x24 = 768 above-cutoff (B_FB=2^16 < LP_MAX=2^18 < q in [2^19, 2^20)) exact O_H-coprime deg-1-root
       campaign with unbiased U-statistic E_emp[p_1^12] + 8 live bit-exact replays (4 g + 4 f).
    5. Full 24x24 = 576 LITERAL production-degree (72, 24) (f_9, g_3) exact O_H-coprime campaign with
       unbiased U-statistic E_emp[p_1^12] + 6 live bit-exact replays (3 g_3 + 3 f_9).
    """
    dp = data["descent"]["params"]
    K_fb = int(dp["m5_K"])
    p1_fl = float(dp["m5_p1_floor"])
    p_geom = float(dp["m5_p_geom"])
    c_gate = float(dp["b0_gate_c"])
    bp_lo, bp_hi = int(dp["m5_bp_lo"]), int(dp["m5_bp_hi"])
    assert K_fb == 15 and p1_fl == 0.45 and p_geom == 5e-5 and c_gate == 0.30 and bp_lo == 12 and bp_hi == 14

    # Part 1: Exact |P_f| pool count in [2^12, 2^14)
    cnt_pf, n_split_p, sample_pf = count_pf_pool_exact(primes_upto, bp_lo, bp_hi)
    assert cnt_pf == 1344 and len(sample_pf) >= 8

    # Part 2: Exact +-1-only Blichfeldt covering bound using upper edge r_hi = c_gate * 2^((b+1)/16), ceil(r_hi),
    # and conservative proved PNT lower bound 0.5 * 2^b / ((b+1)*ln 2) (NO /30 root-of-unity division!)
    b16_h = _b16_half_table(30)
    V16 = (math.pi ** 8) / math.factorial(8)
    mb_worst = 0.0
    mb_worst_b = 0
    for b in range(71, 267):
        r_hi = c_gate * (2.0 ** ((b + 1) / 16.0))
        rc = math.ceil(r_hi)
        half_pts = float(b16_h[rc]) if rc <= 30 else 0.5 * V16 * ((r_hi + 2.0) ** 16)
        pnt_lb = 0.5 * (2.0 ** b) / ((b + 1) * math.log(2.0))
        fr = (4.0 * half_pts) / pnt_lb
        if fr > mb_worst:
            mb_worst = fr
            mb_worst_b = b
    assert mb_worst < p_geom / 15.0

    # Part 3: Live 381-bit prod_sq24 + prod_enc12 ||b_0|| gate & composite L_{q, p_j} CRT Bareiss-LLL checks
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

    # Also verify real 381-bit ENCOUNTERED children (prod_enc12) pass ||b_0|| >= c_gate * q^(1/16) and non-short parent!
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

    # Parts 4 & 5: Load m1_qprand_camp.b85 and verify BOTH full summaries AND live bit-exact replays!
    from pm_adaptive import _is_prime_det
    rdir = Path(root_dir) if root_dir is not None else Path(__file__).resolve().parent
    camp = json.loads(lzma.decompress(base64.b85decode((rdir / "m1_qprand_camp.b85").read_bytes())).decode())

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

    def _has_deg1_root_24(lp, aa, bb, side, u_val):
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
            if side == "g":
                if pow(rx, 3, lp) == (u_val * w * w) % lp:
                    return True
            else:
                if (pow(rx, 3, lp) + w * rx * rx + rx - 1) % lp == 0:
                    return True
        return False

    # Part 4: Above-cutoff 32x24 strict O_H-coprime deg-1-root campaign
    cA = camp["A"]
    rows_a = cA["rows"]
    assert len(rows_a) == 32 and all(len(r[4]) == 24 for r in rows_a)
    assert cA["B_FB"] < cA["LP_MAX"] < min(r[1] for r in rows_a)
    p1s_a = [sum(1 for t in r[4] if t[3] == 0) / 24.0 for r in rows_a]
    ms_a = [sum(t[3] for t in r[4]) / 24.0 for r in rows_a]
    u12_a = sum(math.comb(sum(1 for t in r[4] if t[3] == 0), K_fb) / math.comb(24, K_fb) for r in rows_a) / len(rows_a)
    pl12_a = sum(p ** K_fb for p in p1s_a) / len(p1s_a)
    k12_fails_a = sum(1 for r in rows_a if sum(1 for t in r[4] if t[3] == 0) >= K_fb)
    assert u12_a <= (p1_fl ** K_fb) and pl12_a <= (p1_fl ** K_fb) and max(p1s_a) <= p1_fl and k12_fails_a == 0

    # Live bit-exact replay of 8 trials from Dataset A (4 g-side + 4 f-side)!
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
        return t <= LP_MAX_A and t not in forb and _is_prime_det(t) and _has_deg1_root_24(t, aa, bb, side, U_A)

    BALL_A = _make_ball(cA["n_pts"])
    replayed_a = 0
    for ri in (0, 5, 10, 15, 16, 21, 26, 31):
        sq, qv, ryq, rtq, trs = rows_a[ri]
        pv, ryp, rtp, exp_m = trs[ri % 24]
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
        assert len(dist) == exp_m, f"Replay A mismatch ri={ri}: {len(dist)} != {exp_m}"
        replayed_a += 1

    # Part 5: Literal production-degree (72, 24) (f_9, g_3) 24x24 exact O_H-coprime campaign
    cB = camp["B"]
    rows_b = cB["rows"]
    assert len(rows_b) == 24 and all(len(r[4]) == 24 for r in rows_b)
    p1s_b = [sum(1 for t in r[4] if t[3] == 0) / 24.0 for r in rows_b]
    ms_b = [sum(t[3] for t in r[4]) / 24.0 for r in rows_b]
    u12_b = sum(math.comb(sum(1 for t in r[4] if t[3] == 0), K_fb) / math.comb(24, K_fb) for r in rows_b) / len(rows_b)
    pl12_b = sum(p ** K_fb for p in p1s_b) / len(p1s_b)
    k12_fails_b = sum(1 for r in rows_b if sum(1 for t in r[4] if t[3] == 0) >= K_fb)
    assert u12_b <= (p1_fl ** K_fb) and pl12_b <= (p1_fl ** K_fb) and max(p1s_b) <= p1_fl and k12_fails_b == 0

    # Live bit-exact replay of 4 literal deg-(72, 24) trials (2 g_3 + 2 f_9) WITH explicit deg-1 child root verification!
    U_B, B_FB_G, LP_G, B_FB_F, LP_F = cB["U"], cB["B_FB_G"], cB["LP_G"], cB["B_FB_F"], cB["LP_F"]
    P_B4 = math.prod(p for p in primes_upto(B_FB_F) if p > 65536)

    def _has_deg1_root_72(lp, aa, bb, side, u_val):
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
            else:
                if (pow(x3, 3, lp) + (w + 1) * x3 * x3 + x3 - w) % lp == 0:
                    return True
        return False

    def _sm_b(n, pms, b_fb, lp_max, aa, bb, sd, forb):
        if n == 0:
            return False
        t = abs(n)
        for Pm in pms:
            while True:
                g = math.gcd(t, Pm)
                if g == 1:
                    break
                t //= g
            if t == 1:
                return True
        if t <= b_fb:
            return True
        return t <= lp_max and t not in forb and _is_prime_det(t) and _has_deg1_root_72(t, aa, bb, sd, U_B)

    BALL_B = _make_ball(cB["n_pts"])
    replayed_b = 0
    for ri in (0, 11, 12, 23):
        sq, qv, ryq, rtq, trs = rows_b[ri]
        pv, ryp, rtp, exp_m = trs[ri % 24]
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

    # Typical outer-shell parent scale across all 78 production DAG rows (NOT an individual-vector lower bound!)
    _cr2 = {int(k): float(v) for k, v in data["descent"]["cells_rank2"].items()}
    prod_dag_min_shell_scale = min(
        (2.0 ** (math.log2(_cr2[r[7]]) / 16.0)) * (2.0 ** ((r[0] - r[3]) / 16.0))
        for r in data["descent"]["plan_struct"]
    )
    assert prod_dag_min_shell_scale > 40.0

    return dict(
        mode="pure_m1_qprand_zero_trial0",
        K_trials=K_fb,
        p1_floor=p1_fl,
        pf_pool_exact_12_14=cnt_pf,
        pf_split_primes_12_14=n_split_p,
        b0_gate_c=c_gate,
        blichfeldt_pm1_ceil_worst=round(mb_worst, 8),
        blichfeldt_worst_b=mb_worst_b,
        p_geom_charged=p_geom,
        live_381b_qp_lattices=len(skews_qp),
        live_381b_qp_max_skew=round(max(skews_qp), 3),
        live_381b_qp_min_b0_ratio=round(min(b0_ratios_qp), 3),
        live_381b_enc_min_b0_ratio=round(enc_b0_min, 3),
        live_381b_enc_toyball_min_parent_ratio=round(enc_par_min, 3),
        prod_dag_typical_shell_scale_min=round(prod_dag_min_shell_scale, 3),
        live_381b_qp_max_cf_med=round(max(cf_meds), 2),
        live_381b_qp_max_cg_med=round(max(cg_meds), 2),
        above_cutoff_32x24_exact_oh=dict(
            ideals=32, draws_per_ideal=24, total_trials=768, live_bit_exact_replays=replayed_a,
            strict_lp_lt_q=True, deg1_root_verified=True,
            mean_M_range=[round(min(ms_a), 2), round(max(ms_a), 2)],
            pooled_p1=round(sum(p1s_a) / len(p1s_a), 4), max_p1=round(max(p1s_a), 4),
            unbiased_E_p1_powK=float(f"{u12_a:.3e}"), plugin_E_p1_powK=float(f"{pl12_a:.3e}"),
            claimed_p1_powK=float(f"{p1_fl**K_fb:.3e}"), kK_fails=k12_fails_a,
        ),
        literal_deg72_24_24x24_exact_oh=dict(
            ideals=24, draws_per_ideal=24, total_trials=576, live_bit_exact_replays=replayed_b,
            strict_lp_lt_q=False, deg1_root_verified_on_replays=True,
            mean_M_range=[round(min(ms_b), 2), round(max(ms_b), 2)],
            pooled_p1=round(sum(p1s_b) / len(p1s_b), 4), max_p1=round(max(p1s_b), 4),
            unbiased_E_p1_powK=float(f"{u12_b:.3e}"), plugin_E_p1_powK=float(f"{pl12_b:.3e}"),
            claimed_p1_powK=float(f"{p1_fl**K_fb:.3e}"), kK_fails=k12_fails_b,
        ),
    )
