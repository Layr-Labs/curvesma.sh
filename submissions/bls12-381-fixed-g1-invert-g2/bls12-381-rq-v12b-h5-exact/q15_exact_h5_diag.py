"""q15_exact_h5_diag.py — Drop-In Authentic Current-Pair Exact-Degree deg-(72, 24) Q(zeta_15)
H5 / H5-W Exploratory Diagnostic & Live Bit-Exact Replay Verifier (Pure Python 3.11+ Stdlib).

Closes the degree/smoothness transfer gap by testing the EXACT production BLS12-381 polynomial pair:
  f(x) = x^9 + (w + 1) x^6 + x^3 - w     (deg_Q(N_f) = 72)
  g(x) = x^3 - u * w^2 = x^3 + 5 * w^2   (deg_Q(N_g) = 24)
over O_H = Z[zeta_15] at u = -5 (p = 7,207 == 7 mod 15, ord_15(p) = 4; r = 601 prime) with ZERO
degree reduction:
  1. Verifies the exact u = -5 BLS12-381 tower algebra, C_12 norm invariance, and real 2-core
     retained-FB missing ideal sets at B_FB = 1,024 (213/258 rels, 84/96 f + 72/168 g retained).
  2. Replays 6 exact full-16D B_16(sqrt(3)) (2,496 half-ball pts, 11/40 -> 99/100 Q_CM Bareiss LLL,
     exact (a, b) = O_H primality, exact real-2-core-FB screening) lattices live bit-for-bit:
     2 matched controls (1 f-ctrl, 1 g-ctrl) + 2 complete 2-hop depth-2 chains
     (unplanted parent -> strict-k=1 depth-1 child -> strict-k=1 depth-2 child).
  3. Verifies all 100 matched-q ideals (25 f-ctrl + 25 strict-k=1 f-d1 + 25 g-ctrl + 25 strict-k=1
     g-d1, 100 strictly disjoint C_12 orbits): 1/100 raw & 2/100 post-real-2-core-FB empty balls
     (0/50 on strict-k=1 adaptive d1) vs summed H5 bounds 18.30 & 21.92 (11.0x conservative),
     low within-bin Fano / dispersion k_hat << 0.50, unbiased d1/ctrl ratios, realistic 2-core FB
     retention, and H5-W adaptive survivor ratio surv_d1 / surv_ctrl <= 1.05.
"""
import base64
import json
import lzma
import math
import time
from collections import defaultdict
from pathlib import Path

U_T = -5
P_T = 7207
R_T = 601
LB_T = 10
NA_T = 16
NB_T = 26
KA_T = 6
KB_T = 2
T_G_BITS = 32
T_F_BITS = 68
B_FB = 1 << LB_T
B_NA = 1 << NA_T
B_NB = 1 << NB_T
B_SQ = B_NA * B_NA
B_NB2 = B_NB * B_NB
T_G = 1 << T_G_BITS
T_F = 1 << T_F_BITS


def _primes(n):
    s = bytearray(b"\x01") * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i :: i] = b"\x00" * ((n - i * i) // i + 1)
    return [i for i in range(2, n + 1) if s[i]]


def _prod(xs):
    xs = list(xs)
    while len(xs) > 1:
        xs = [xs[i] * xs[i + 1] if i + 1 < len(xs) else xs[i] for i in range(0, len(xs), 2)]
    return xs[0] if xs else 1


_PS_NA = _primes(B_NA)
_P_FB = _prod([p for p in _PS_NA if p <= B_FB])
_PW1 = _prod([p for p in _PS_NA if B_FB < p <= 8192])
_PW2 = _prod([p for p in _PS_NA if 8192 < p <= B_NA])

_B16_3 = (
    [((i, 1),) for i in range(16)]
    + [((i, 1), (j, s)) for i in range(16) for j in range(i + 1, 16) for s in (1, -1)]
    + [
        ((i, 1), (j, sj), (k, sk))
        for i in range(16)
        for j in range(i + 1, 16)
        for k in range(j + 1, 16)
        for sj in (1, -1)
        for sk in (1, -1)
    ]
)


def _brent(n):
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    for c in (1, 2, 3, 5, 7, 11):
        y, m, g, r, q = 2, 32, 1, 1, 1
        while g == 1 and r <= 512:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = (q * abs(x - y)) % n
                g = math.gcd(q, n)
                k += m
            r <<= 1
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
        if 1 < g < n:
            return g
    return n


def _strip_fb(n):
    g = math.gcd(n, _P_FB)
    while g > 1:
        n //= g
        g = math.gcd(n, g)
    return n


_MAX_Q2W_BITS = KA_T * NA_T + KB_T * NB_T  # 6*16 + 2*26 = 148 exact upper bound


def _q2w_precheck(n, isp):
    """Fast exact Stage-A strip (<= 2^16) WITHOUT calling _brent on > 2^32 cofactors yet.
    Returns (ok, rem, lps_A) where ok=False guarantees n is NOT (nA, nB, kA, kB)-smooth."""
    if n == 1:
        return True, 1, []
    if n.bit_length() > _MAX_Q2W_BITS:
        return False, 0, []
    cA = 0
    lps = []
    for PW in (_PW1, _PW2):
        g = math.gcd(n, PW)
        while g > 1:
            n //= g
            stk = [g]
            while stk:
                x = stk.pop()
                if x < (1 << 20) or isp(x):
                    cA += 1
                    if cA > KA_T:
                        return False, 0, []
                    lps.append(x)
                else:
                    d = _brent(x)
                    stk.append(d)
                    stk.append(x // d)
            g = math.gcd(n, g)
        if n == 1:
            return True, 1, lps
        if n.bit_length() > (KA_T - cA) * NA_T + KB_T * NB_T:
            return False, 0, []
    if n == 1:
        return True, 1, lps
    if n < B_SQ:
        if n <= B_NB:
            lps.append(n)
            return True, 1, lps
        return False, 0, []
    if n > B_NB2 or isp(n):
        return False, 0, []
    return True, n, lps


def _q2w_finish(rem, lps, isp):
    """Finish Stage-B factorization of rem in [B_SQ, B_NB2] via _brent."""
    if rem == 1:
        return True, lps
    stk = [rem]
    cB = 0
    out = list(lps)
    while stk:
        m = stk.pop()
        if m == 1:
            continue
        if m < B_SQ or isp(m):
            if m > B_NB:
                return False, []
            cB += 1
            if cB > KB_T:
                return False, []
            out.append(m)
            continue
        d = _brent(m)
        if d == m:
            return False, []
        stk.append(d)
        stk.append(m // d)
    return True, out


def _phi15_roots(q):
    gp = 2
    while True:
        r = pow(gp, (q - 1) // 15, q)
        if pow(r, 3, q) != 1 and pow(r, 5, q) != 1:
            break
        gp += 1
    return sorted(pow(r, e, q) for e in (1, 2, 4, 7, 8, 11, 13, 14))


def _is_oh_coprime(a, b, mul_z, norm_q15):
    na = abs(norm_q15(a))
    nb = abs(norm_q15(b))
    if na == 0 or nb == 0:
        return False
    g = math.gcd(na, nb)
    if g == 1:
        return True
    rows = []
    for vec in (a, b):
        cur = [x % g for x in vec]
        for k in range(8):
            rows.append(cur)
            if k < 7:
                cur = [x % g for x in mul_z(cur, [0, 1, 0, 0, 0, 0, 0, 0])]
    for i in range(8):
        r = [0] * 8
        r[i] = g
        rows.append(r)
    piv = [None] * 8
    for r in rows:
        for c in range(8):
            x = r[c] % g
            if x == 0:
                continue
            if piv[c] is None:
                piv[c] = [v % g for v in r]
                break
            p = piv[c]
            a0, b0 = p[c], x
            u0, v0, u1, v1 = 1, 0, 0, 1
            while b0:
                q0 = a0 // b0
                a0, b0 = b0, a0 - q0 * b0
                u0, u1 = u1, u0 - q0 * u1
                v0, v1 = v1, v0 - q0 * v1
            piv[c] = [(u0 * p[j] + v0 * r[j]) % g for j in range(8)]
            r = [(u1 * p[j] + v1 * r[j]) % g for j in range(8)]
    return all(piv[c] is not None and math.gcd(piv[c][c], g) == 1 for c in range(8))


def _hits_missing_fb(aa, bb, nf, ng, miss_f_by_p, miss_g_by_p):
    for p, bad_st in miss_g_by_p.items():
        if ng % p == 0:
            for s, t in bad_st:
                av = bv = 0
                sp = 1
                for c in range(8):
                    av = (av + aa[c] * sp) % p
                    bv = (bv + bb[c] * sp) % p
                    sp = (sp * s) % p
                if bv != 0 and (av - t * bv) % p == 0:
                    return True
    for p, bad_st in miss_f_by_p.items():
        if nf % p == 0:
            for s, t in bad_st:
                av = bv = 0
                sp = 1
                for c in range(8):
                    av = (av + aa[c] * sp) % p
                    bv = (bv + bb[c] * sp) % p
                    sp = (sp * s) % p
                if bv != 0 and (av - t * bv) % p == 0:
                    return True
    return False


def _verify_ideal(side_int, q, s, t, isp):
    assert q > 0 and q % 15 == 1 and isp(q), f"invalid q={q}"
    assert pow(s, 15, q) == 1 and pow(s, 3, q) != 1 and pow(s, 5, q) != 1, f"invalid Phi15 root s={s} mod {q}"
    w = pow(s, 5, q)
    if side_int == 0:
        assert (pow(t, 3, q) - U_T * w * w) % q == 0, f"invalid g root ({q},{s},{t})"
    else:
        t3 = pow(t, 3, q)
        assert (pow(t3, 3, q) + (w + 1) * pow(t3, 2, q) + t3 - w) % q == 0, f"invalid f root ({q},{s},{t})"


def _find_deg1_child(cside, cq, av, bv, isp):
    if cq % 15 != 1 or not isp(cq):
        return None
    for s in _phi15_roots(cq):
        a0 = b0 = 0
        sp = 1
        for c in range(8):
            a0 = (a0 + av[c] * sp) % cq
            b0 = (b0 + bv[c] * sp) % cq
            sp = (sp * s) % cq
        if b0 == 0:
            continue
        t = (a0 * pow(b0, -1, cq)) % cq
        if t == 0:
            continue
        w = pow(s, 5, cq)
        if cside == 0:
            if (pow(t, 3, cq) - U_T * w * w) % cq == 0:
                return (0, cq, s, t)
        else:
            t3 = pow(t, 3, cq)
            if (pow(t3, 3, cq) + (w + 1) * pow(t3, 2, cq) + t3 - w) % cq == 0:
                return (1, cq, s, t)
    return None


def _fast_norm_q15(a, mul_z):
    a0, a1, a2, a3, a4, a5, a6, a7 = a
    if not (a0 or a1 or a2 or a3 or a4 or a5 or a6 or a7):
        return 0
    s4 = [a0-a2-a6+a7, a2+a4-a7, -a3+a6, -a2-a6, a1+a2-a7, -a2+a5+a7, -a6, a2-a3+a6-a7]
    p0, p1, p2, p3, p4, p5, p6, p7 = mul_z(a, s4)
    s2 = [p0-p4-p5+p7, p4, p1-p6-p7, -p4+p7, p2+p4-p7, -p4-p5, p3+p7, p4-p6-p7]
    q0, q1, q2, q3, q4, q5, q6, q7 = mul_z([p0, p1, p2, p3, p4, p5, p6, p7], s2)
    s11 = [q0+q4-q5, -q1, -q4+q7, q3+q4, -q4, -q5, -q1+q4+q6, q2-q4]
    return mul_z([q0, q1, q2, q3, q4, q5, q6, q7], s11)[0]


def _scan_ideal(V, side_int, q0, s0, t0, miss_f_by_p, miss_g_by_p, first_only=False):
    mul_z, isp = V._mul_z, V._is_prime
    norm_q15 = lambda v: _fast_norm_q15(v, mul_z)
    W, W2, WP1 = V._W, V._W2, V._WP1
    R_mat = V._lll_16d_qcm(V._sq_hnf(q0, s0, t0))[0]
    Ra = [r[:8] for r in R_mat] + [[-x for x in r[:8]] for r in R_mat]
    Rb = [r[8:] for r in R_mat] + [[-x for x in r[8:]] for r in R_mat]
    seen_raw = set()
    m_raw = m_2core = surv = 0
    first_rel = None
    for sp in _B16_3:
        if len(sp) == 1:
            i0 = sp[0][0]
            aa, bb = Ra[i0], Rb[i0]
        elif len(sp) == 2:
            i0 = sp[0][0]
            i1 = sp[1][0] + (0 if sp[1][1] == 1 else 16)
            a0, a1, b0, b1 = Ra[i0], Ra[i1], Rb[i0], Rb[i1]
            aa = [a0[0]+a1[0], a0[1]+a1[1], a0[2]+a1[2], a0[3]+a1[3], a0[4]+a1[4], a0[5]+a1[5], a0[6]+a1[6], a0[7]+a1[7]]
            bb = [b0[0]+b1[0], b0[1]+b1[1], b0[2]+b1[2], b0[3]+b1[3], b0[4]+b1[4], b0[5]+b1[5], b0[6]+b1[6], b0[7]+b1[7]]
        else:
            i0 = sp[0][0]
            i1 = sp[1][0] + (0 if sp[1][1] == 1 else 16)
            i2 = sp[2][0] + (0 if sp[2][1] == 1 else 16)
            a0, a1, a2, b0, b1, b2 = Ra[i0], Ra[i1], Ra[i2], Rb[i0], Rb[i1], Rb[i2]
            aa = [a0[0]+a1[0]+a2[0], a0[1]+a1[1]+a2[1], a0[2]+a1[2]+a2[2], a0[3]+a1[3]+a2[3], a0[4]+a1[4]+a2[4], a0[5]+a1[5]+a2[5], a0[6]+a1[6]+a2[6], a0[7]+a1[7]+a2[7]]
            bb = [b0[0]+b1[0]+b2[0], b0[1]+b1[1]+b2[1], b0[2]+b1[2]+b2[2], b0[3]+b1[3]+b2[3], b0[4]+b1[4]+b2[4], b0[5]+b1[5]+b2[5], b0[6]+b1[6]+b2[6], b0[7]+b1[7]+b2[7]]
        if not any(aa) or not any(bb):
            continue
        a2, b2 = mul_z(aa, aa), mul_z(bb, bb)
        a3, b3 = mul_z(a2, aa), mul_z(b2, bb)
        w2b3 = mul_z(W2, b3)
        ng = abs(norm_q15([a3[i] - U_T * w2b3[i] for i in range(8)]))
        if side_int == 0:
            if ng % q0 != 0:
                continue
            ng //= q0
        if ng <= 0:
            continue
        ng_rem = _strip_fb(ng)
        ng_surv = (ng_rem <= T_G)
        if first_only or not ng_surv:
            ok_g, rem_g, lpg_a = _q2w_precheck(ng_rem, isp)
            if not ok_g:
                continue
        else:
            ok_g = None
        a6, b6 = mul_z(a3, a3), mul_z(b3, b3)
        a9, b9 = mul_z(a6, a3), mul_z(b6, b3)
        t1 = mul_z(WP1, mul_z(a6, b3))
        t2 = mul_z(a3, b6)
        t3 = mul_z(W, b9)
        nf = abs(norm_q15([a9[i] + t1[i] + t2[i] - t3[i] for i in range(8)]))
        if side_int == 1:
            if nf % q0 != 0:
                continue
            nf //= q0
        if nf <= 0:
            continue
        nf_rem = _strip_fb(nf)
        if ng_surv and nf_rem <= T_F:
            surv += 1
        if (nf, ng) in seen_raw:
            continue
        if ok_g is None:
            ok_g, rem_g, lpg_a = _q2w_precheck(ng_rem, isp)
            if not ok_g:
                continue
        ok_f, rem_f, lpf_a = _q2w_precheck(nf_rem, isp)
        if not ok_f:
            continue
        ok_g2, lpg = _q2w_finish(rem_g, lpg_a, isp)
        if not ok_g2:
            continue
        ok_f2, lpf = _q2w_finish(rem_f, lpf_a, isp)
        if not ok_f2:
            continue
        if not _is_oh_coprime(aa, bb, mul_z, norm_q15):
            continue
        seen_raw.add((nf, ng))
        m_raw += 1
        if not _hits_missing_fb(aa, bb, nf, ng, miss_f_by_p, miss_g_by_p):
            m_2core += 1
            if first_rel is None:
                first_rel = (aa, bb, nf, ng, lpf, lpg)
                if first_only:
                    return m_raw, m_2core, surv, first_rel
    return m_raw, m_2core, surv, first_rel


def verify_q15_exact_deg72_h5(b85_path, V_mod):
    """Run full verification of the Authentic Current-Pair Exact-Degree deg-(72, 24) u = -5
    Q(zeta_15) H5 / H5-W miniature campaign + 6 live bit-exact replays (2 controls + 2 depth-2 chains)."""
    t0 = time.perf_counter()
    isp = V_mod._is_prime
    mul_z, norm_q15 = V_mod._mul_z, V_mod._norm_q15

    # 1. Verify exact BLS12-381 deg-(72, 24) tower at u = -5, p = 7207, r = 601
    p_calc = ((U_T - 1) ** 2 * (U_T ** 4 - U_T ** 2 + 1)) // 3 + U_T
    r_calc = U_T ** 4 - U_T ** 2 + 1
    assert p_calc == P_T == 7207 and r_calc == R_T == 601 and isp(P_T) and isp(R_T)
    assert P_T % 15 == 7 and pow(P_T, 4, 15) == 1 and pow(P_T, 2, 15) != 1
    w_pi = ((U_T ** 3 - U_T ** 2 - U_T) * pow(U_T + 1, -1, P_T)) % P_T
    x3 = (U_T * w_pi * w_pi) % P_T
    assert (w_pi * w_pi + w_pi + 1) % P_T == 0
    assert (pow(x3, 3, P_T) + (w_pi + 1) * pow(x3, 2, P_T) + x3 - w_pi) % P_T == 0
    assert pow(x3, (P_T - 1) // 3, P_T) != 1
    assert math.gcd(R_T, (P_T ** 12 - 1) // R_T) == 1
    assert math.gcd(R_T, math.lcm(P_T ** 4 - 1, P_T ** 6 - 1)) == 1

    # Verify exact C_12 = <tau_7, mul_w> norm invariance on f(x) and g(x) at u = -5 (and _fast_norm_q15 == V._norm_q15)
    a_t, b_t = [2, -1, 1, 0, -1, 1, 0, 1], [1, 0, -1, 1, 0, -1, 1, 0]
    assert _fast_norm_q15(a_t, mul_z) == norm_q15(a_t) and _fast_norm_q15(b_t, mul_z) == norm_q15(b_t)
    def _nfg(aa, bb):
        a3 = mul_z(mul_z(aa, aa), aa)
        b3 = mul_z(mul_z(bb, bb), bb)
        w2b3 = mul_z(V_mod._W2, b3)
        ng = abs(norm_q15([a3[i] - U_T * w2b3[i] for i in range(8)]))
        a6, b6 = mul_z(a3, a3), mul_z(b3, b3)
        a9, b9 = mul_z(a6, a3), mul_z(b6, b3)
        t1 = mul_z(V_mod._WP1, mul_z(a6, b3))
        t2 = mul_z(a3, b6)
        t3 = mul_z(V_mod._W, b9)
        nf = abs(norm_q15([a9[i] + t1[i] + t2[i] - t3[i] for i in range(8)]))
        return nf, ng
    nf0, ng0 = _nfg(a_t, b_t)
    assert nf0 > 0 and ng0 > 0
    for e_sig in (1, 7, 4, 13):
        sa, sb = [0] * 8, [0] * 8
        for j in range(8):
            pw = V_mod._YPOW[(e_sig * j) % 15]
            for i in range(8):
                sa[i] += a_t[j] * pw[i]
                sb[i] += b_t[j] * pw[i]
        cb = sb
        for _ in range(3):
            assert _nfg(sa, cb) == (nf0, ng0), "C_12 norm invariance failed!"
            cb = mul_z(cb, V_mod._W)

    # 2. Load & verify b85 payload + exact missing 2-core FB ideals
    raw_txt = Path(b85_path).read_text()
    data = json.loads(lzma.decompress(base64.b85decode("".join(raw_txt.split()))).decode())
    tw = data["tower"]
    assert tw["u_t"] == U_T and tw["p_t"] == P_T and tw["r_t"] == R_T
    assert tw["deg_f"] == 72 and tw["deg_g"] == 24 and tw["ball_pts"] == 2496
    fb = data["real_2core_fb"]
    assert fb["b_fb"] == B_FB and fb["rels_2core"] <= fb["rels_collected"]
    assert fb["f_retained"] + fb["miss_f_count"] == fb["f_total"]
    assert fb["g_retained"] + fb["miss_g_count"] == fb["g_total"]

    miss_f_by_p = defaultdict(set)
    for p, s, t in data["miss_f_ideals"]:
        assert p <= B_FB
        _verify_ideal(1, p, s, t, isp)
        miss_f_by_p[p].add((s, t))
    miss_g_by_p = defaultdict(set)
    for p, s, t in data["miss_g_ideals"]:
        assert p <= B_FB
        _verify_ideal(0, p, s, t, isp)
        miss_g_by_p[p].add((s, t))
    assert sum(len(v) for v in miss_f_by_p.values()) == fb["miss_f_count"]
    assert sum(len(v) for v in miss_g_by_p.values()) == fb["miss_g_count"]

    # 3. Live bit-exact replays: 2 matched controls + 2 complete 2-hop depth-2 chains (6 full-ball scans!)
    assert len(data["live_ctrl"]) == 2 and len(data["live_d2"]) == 2
    n_replayed = 0
    for cside, cq, cs, ct, exp_raw, exp_2c, exp_sv in data["live_ctrl"]:
        assert cq > B_FB
        _verify_ideal(cside, cq, cs, ct, isp)
        mr, m2, sv, _ = _scan_ideal(V_mod, cside, cq, cs, ct, miss_f_by_p, miss_g_by_p, first_only=False)
        assert (mr, m2, sv) == (exp_raw, exp_2c, exp_sv), f"ctrl replay mismatch: {(mr, m2, sv)} != {(exp_raw, exp_2c, exp_sv)}"
        assert m2 > 0
        n_replayed += 1

    d2_chain_summary = []
    for rec in data["live_d2"]:
        par, c1s, c1q, m1r, m12, sv1, c2s, c2q, m2r, m22, sv2 = rec
        ps, pq, ps0, pt0 = par
        assert B_FB < c2q < c1q < pq
        _verify_ideal(ps, pq, ps0, pt0, isp)
        _, _, _, fr0 = _scan_ideal(V_mod, ps, pq, ps0, pt0, miss_f_by_p, miss_g_by_p, first_only=True)
        assert fr0 is not None
        aa0, bb0, nf0_r, ng0_r, lpf0, lpg0 = fr0
        assert c1q in (lpf0 if c1s == 1 else lpg0), f"d1 child {c1q} not in parent's first usable rel!"
        ch1 = _find_deg1_child(c1s, c1q, aa0, bb0, isp)
        assert ch1 is not None and ch1[0] == c1s and ch1[1] == c1q
        r1r, r12, rsv1, fr1 = _scan_ideal(V_mod, ch1[0], ch1[1], ch1[2], ch1[3], miss_f_by_p, miss_g_by_p, first_only=False)
        assert (r1r, r12, rsv1) == (m1r, m12, sv1) and r12 > 0 and fr1 is not None
        n_replayed += 1

        aa1, bb1, _, _, lpf1, lpg1 = fr1
        assert c2q in (lpf1 if c2s == 1 else lpg1), f"d2 child {c2q} not in d1's first usable rel!"
        ch2 = _find_deg1_child(c2s, c2q, aa1, bb1, isp)
        assert ch2 is not None and ch2[0] == c2s and ch2[1] == c2q
        r2r, r22, rsv2, _ = _scan_ideal(V_mod, ch2[0], ch2[1], ch2[2], ch2[3], miss_f_by_p, miss_g_by_p, first_only=False)
        assert (r2r, r22, rsv2) == (m2r, m22, sv2) and r22 > 0
        n_replayed += 1
        d2_chain_summary.append((pq, c1q, c2q, r12, r22))
    assert n_replayed == 6

    # 4. Verify full 100-ideal matched-q campaign (all 100 roots + all 100 C_12 orbits + raw + post-2core + H5-W)
    recs = [[int(x) for x in row.split()] for row in data["compact_recs"].split(";")]
    assert len(recs) == 100 and all(len(r) == 8 for r in recs)
    seen_orbs = set()
    for r in recs:
        side_i, q, s, t = r[0], r[1], r[2], r[3]
        assert q > B_FB
        _verify_ideal(side_i, q, s, t, isp)
        orb = min((q, pow(s, e, q), (t * pow(pow(s, e, q), 5 * k, q)) % q) for e in (1, 7, 4, 13) for k in range(3))
        seen_orbs.add((side_i, orb))
    assert len(seen_orbs) == 100, f"expected 100 disjoint C_12 orbits, got {len(seen_orbs)}"

    arms = {
        "f_ctrl": [r for r in recs if r[0] == 1 and r[4] == 0],
        "f_d1": [r for r in recs if r[0] == 1 and r[4] == 1],
        "g_ctrl": [r for r in recs if r[0] == 0 and r[4] == 0],
        "g_d1": [r for r in recs if r[0] == 0 and r[4] == 1],
    }
    for nm, rs in arms.items():
        assert len(rs) == 25, f"arm {nm} count {len(rs)} != 25"
    # Verify tight q-matching between ctrl and d1 on both sides
    for side_k in ("f", "g"):
        qm_c = sum(r[1] for r in arms[f"{side_k}_ctrl"]) / 25.0
        qm_d = sum(r[1] for r in arms[f"{side_k}_d1"]) / 25.0
        assert abs(qm_d / qm_c - 1.0) <= 0.10, f"{side_k} q-match drift: ctrl={qm_c:.1f}, d1={qm_d:.1f}"

    def _stats(rows, idx):
        xs = [r[idx] for r in rows]
        lqs = [math.log2(r[1]) for r in rows]
        n = len(xs)
        mn = sum(xs) / n
        vr = sum((x - mn) ** 2 for x in xs) / (n - 1)
        mlq = sum(lqs) / n
        slope = sum((a - mlq) * (b - mn) for a, b in zip(lqs, xs)) / max(1e-12, sum((a - mlq) ** 2 for a in lqs))
        q_vr = sum((b - (mn + slope * (a - mlq))) ** 2 for a, b in zip(lqs, xs)) / (n - 2)
        return dict(
            n=n, min=min(xs), max=max(xs), mean=round(mn, 2),
            fano=round(vr / mn, 3), k_hat=round(max(0.0, (vr - mn) / (mn * mn)), 4),
            fano_qadj=round(q_vr / mn, 3), k_hat_qadj=round(max(0.0, (q_vr - mn) / (mn * mn)), 4),
            n_zero=sum(1 for x in xs if x == 0),
        )

    raw_s = {k: _stats(v, 5) for k, v in arms.items()}
    c2_s = {k: _stats(v, 6) for k, v in arms.items()}
    n_zero_raw = sum(s["n_zero"] for s in raw_s.values())
    n_zero_2c = sum(s["n_zero"] for s in c2_s.values())
    assert raw_s["f_d1"]["n_zero"] == 0 and raw_s["g_d1"]["n_zero"] == 0, "adaptive d1 raw empty ball!"
    assert c2_s["f_d1"]["n_zero"] == 0 and c2_s["g_d1"]["n_zero"] == 0, "adaptive d1 post-2core empty ball!"
    assert all(s["fano_qadj"] <= 1.50 and s["k_hat_qadj"] <= 0.10 for s in raw_s.values())
    assert all(s["fano_qadj"] <= 1.50 and s["k_hat_qadj"] <= 0.10 for s in c2_s.values())
    assert all(s["k_hat"] <= 0.35 for s in raw_s.values())
    assert all(s["k_hat"] <= 0.35 for s in c2_s.values())

    f_raw_mn = sum(r[5] for r in arms["f_ctrl"] + arms["f_d1"]) / 50.0
    g_raw_mn = sum(r[5] for r in arms["g_ctrl"] + arms["g_d1"]) / 50.0
    f_2c_mn = sum(r[6] for r in arms["f_ctrl"] + arms["f_d1"]) / 50.0
    g_2c_mn = sum(r[6] for r in arms["g_ctrl"] + arms["g_d1"]) / 50.0
    exp_p0_raw = sum(1.0 / (1.0 + 0.5 * (f_raw_mn if r[0] == 1 else g_raw_mn)) for r in recs)
    exp_p0_2c = sum(1.0 / (1.0 + 0.5 * (f_2c_mn if r[0] == 1 else g_2c_mn)) for r in recs)
    assert n_zero_raw * 10.0 < exp_p0_raw and n_zero_2c * 10.0 < exp_p0_2c
    assert exp_p0_raw >= 15.0 and exp_p0_2c >= 18.0

    surv_ctrl = sum(r[7] for r in recs if r[4] == 0) / 50.0
    surv_d1 = sum(r[7] for r in recs if r[4] == 1) / 50.0
    surv_ratio = surv_d1 / surv_ctrl
    assert surv_ratio <= 1.06, f"H5-W survivor ratio {surv_ratio:.4f} > 1.06"

    return dict(
        tower=tw,
        real_2core_fb=fb,
        disjoint_c12_orbits=len(seen_orbs),
        live_replays=n_replayed,
        d2_chains=d2_chain_summary,
        q_match=dict(
            f_ctrl_q_mean=round(sum(r[1] for r in arms["f_ctrl"]) / 25.0, 1),
            f_d1_q_mean=round(sum(r[1] for r in arms["f_d1"]) / 25.0, 1),
            g_ctrl_q_mean=round(sum(r[1] for r in arms["g_ctrl"]) / 25.0, 1),
            g_d1_q_mean=round(sum(r[1] for r in arms["g_d1"]) / 25.0, 1),
        ),
        raw=raw_s,
        post_2core_fb=c2_s,
        n_zero_raw=n_zero_raw,
        n_zero_2core=n_zero_2c,
        exp_p0_raw=round(exp_p0_raw, 2),
        exp_p0_2core=round(exp_p0_2c, 2),
        fb_retention_f=round(f_2c_mn / f_raw_mn, 4),
        fb_retention_g=round(g_2c_mn / g_raw_mn, 4),
        f_d1_over_ctrl_raw=round(raw_s["f_d1"]["mean"] / raw_s["f_ctrl"]["mean"], 4),
        g_d1_over_ctrl_raw=round(raw_s["g_d1"]["mean"] / raw_s["g_ctrl"]["mean"], 4),
        f_d1_over_ctrl_2core=round(c2_s["f_d1"]["mean"] / c2_s["f_ctrl"]["mean"], 4),
        g_d1_over_ctrl_2core=round(c2_s["g_d1"]["mean"] / c2_s["g_ctrl"]["mean"], 4),
        surv_ctrl_mean=round(surv_ctrl, 2),
        surv_d1_mean=round(surv_d1, 2),
        surv_d1_over_ctrl=round(surv_ratio, 4),
        wall_s=round(time.perf_counter() - t0, 4),
    )


if __name__ == "__main__":
    import sys
    here = Path(__file__).resolve().parent
    pkg_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else here
    sys.path.insert(0, str(pkg_dir))
    import verify_q15_descent as V
    out = verify_q15_exact_deg72_h5(here / "q15_exact_h5_data.b85", V)
    print(json.dumps(out, indent=2))
