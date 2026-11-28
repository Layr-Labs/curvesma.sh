"""q15_exact_h5_diag.py — Drop-In Authentic Current-Pair Exact-Degree deg-(72, 24) Q(zeta_15)
Fixed-2-Window EA-Gated H5 / H5-W Exploratory Diagnostic & Independent Production-Estimator
mu_pred Calibration Verifier (Pure Python 3.11+ Stdlib, v3).

Evaluates the EXACT production BLS12-381 polynomial pair with ZERO degree reduction:
  f(x) = x^9 + (w + 1) x^6 + x^3 - w     (deg_Q(N_f) = 72)
  g(x) = x^3 - u * w^2 = x^3 + 5 * w^2   (deg_Q(N_g) = 24)
over O_H = Z[zeta_15] at u = -5 (p = 7,207 == 7 mod 15, ord_15(p) = 4; r = 601 prime).

Key properties (v3 — addressing both Codex Tier-1/2 and live Yukon review):
  1. Independent production-estimator mu_pred calibration: For EVERY ONE of the 100 matched-q
     ideals (100 strictly disjoint C_12 orbits), computes mu_pred_raw_i and mu_pred_2core_i
     INDEPENDENTLY using the literal unmodified production estimator q15_descent16.ss_lower2_tau
     (with production ALPHA_F = 2.74, ALPHA_G = 0.77, kA = 6, kB = 1) + production DS.kappa
     linear missing-FB union bound max(0, 1 - miss_f - miss_g) over unique O_H-primitive
     (nf, ng) pairs in B_16(sqrt(3)) — ZERO peeking at observed relation counts. Directly
     calibrates obs_raw / mu_pred_raw and obs_2core / mu_pred_2core and compares actual empty
     balls against the INDEPENDENT summed H5 bound sum_i max(0.001, (1 + 0.5 * mu_pred_i)^{-1}).
  2. Strict deterministic EA gate (usable => surv): Fixed 2-window smoothness
     (B_FB = 2^10 = 1,024, B_A = 2^16, B_B = 2^26, kA = 6, kB = 1) is strictly gated inside
     (T_G = 2^40, T_F = 2^88) after B_FB strip, so 0 <= m_2core <= m_raw <= surv holds
     deterministically on every ideal and every first_rel passes EA and real-2-core-FB screening.
  3. Explicit fixed-2-window scoping: Because deg_Q(f)/16 = 4.5 on g-side makes toy-scale
     full-relation max_lp < q infeasible (requiring q >= 2^40, R >= 3), companion primes may
     exceed q up to B_B = 2^26 while every tracked d1/d2 branch strictly reduces
     B_FB < c2q < c1q < q_parent via a verified v_{cq}(N) == 1 split deg-1 child ideal.
     Full-relation strict max_lp < q multi-level tree descent is tested separately in the
     complementary Q(zeta_15) deg-(32, 24) campaigns (ADAPTIVE.md §2(B)-(C)).
  4. Fast CI default + full diagnostics mode:
       * Default (full_replays=False, <= 0.6s pinned linux/amd64 Docker): verifies exact u = -5
         tower + C_12 invariance, all 108 missing-FB ideal roots, all 100 campaign ideal roots
         & 100 disjoint C_12 orbits, live_ctrl/live_d2 cross-match against compact_recs,
         independent mu_pred calibration, and 1 complete 2-hop depth-2 chain (2 full-ball live
         bit-exact replays + live DS.ss_lower2_tau mu_pred recomputation).
       * Full mode (full_replays=True, ~1.5s Docker, for {"part": "diagnostics"}): replays ALL
         6 full-ball lattices (2 matched controls + 2 complete 2-hop depth-2 chains) + live
         DS.ss_lower2_tau mu_pred recomputation on all 6.
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
KB_T = 1
T_G_BITS = 40
T_F_BITS = 88
MU_MESH = 32
MU_IA = 12
ALPHA_F_PROD = 2.74
ALPHA_G_PROD = 0.77
B_FB = 1 << LB_T
B_NA = 1 << NA_T
B_NB = 1 << NB_T
B_SQ = B_NA * B_NA
B_NB2 = B_NB * B_NB
T_G = 1 << T_G_BITS
T_F = 1 << T_F_BITS
_MAX_Q2W_BITS = KA_T * NA_T + KB_T * NB_T


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


def _q2w_precheck(n, isp):
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
    if n < B_SQ:
        if n <= B_NB:
            lps.append(n)
            return True, 1, lps
        return False, 0, []
    if n > B_NB2 or isp(n):
        return False, 0, []
    return True, n, lps


def _q2w_finish(rem, lps, isp):
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


def _miss_kappa_tables(DS, miss_f_by_p, miss_g_by_p):
    mf_items = [(len(sts) / p, math.log2(p) / float(NB_T)) for p, sts in sorted(miss_f_by_p.items())]
    mg_items = [(len(sts) / p, math.log2(p) / float(NB_T)) for p, sts in sorted(miss_g_by_p.items())]
    mf_cache, mg_cache = {}, {}
    def get_mf(bf):
        v = mf_cache.get(bf)
        if v is None:
            u = bf / float(NB_T)
            v = sum(w * DS.kappa(u, t) for w, t in mf_items)
            mf_cache[bf] = v
        return v
    def get_mg(bg):
        v = mg_cache.get(bg)
        if v is None:
            u = bg / float(NB_T)
            v = sum(w * DS.kappa(u, t) for w, t in mg_items)
            mg_cache[bg] = v
        return v
    return get_mf, get_mg


def _scan_ideal(V, DS, side_int, q0, s0, t0, miss_f_by_p, miss_g_by_p, first_only=False, miss_fns=None):
    mul_z, isp = V._mul_z, V._is_prime
    norm_q15 = lambda v: _fast_norm_q15(v, mul_z)
    W, W2, WP1 = V._W, V._W2, V._WP1
    R_mat = V._lll_16d_qcm(V._sq_hnf(q0, s0, t0))[0]
    Ra = [r[:8] for r in R_mat] + [[-x for x in r[:8]] for r in R_mat]
    Rb = [r[8:] for r in R_mat] + [[-x for x in r[8:]] for r in R_mat]
    seen_all = set()
    seen_raw = set()
    m_raw = m_2core = surv = 0
    mu_raw = mu_2c = 0.0
    first_rel = None
    get_mf, get_mg = miss_fns if (miss_fns is not None and not first_only) else (None, None)
    ss_tau = DS.ss_lower2_tau if get_mf is not None else None
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
        ng_ok_ea = (ng_rem <= T_G)
        if first_only and not ng_ok_ea:
            continue
        if get_mf is None and not ng_ok_ea:
            continue
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
        is_coprime = None
        if get_mf is not None and (nf, ng) not in seen_all:
            seen_all.add((nf, ng))
            is_coprime = _is_oh_coprime(aa, bb, mul_z, norm_q15)
            if is_coprime:
                bf = math.ceil(2.0 * math.log2(nf)) * 0.5
                bg = math.ceil(2.0 * math.log2(ng)) * 0.5
                p_pt = (
                    ss_tau(bf, float(LB_T), float(NB_T), MU_IA, KA_T, KB_T, MU_MESH, ALPHA_F_PROD, float(T_F_BITS))
                    * ss_tau(bg, float(LB_T), float(NB_T), MU_IA, KA_T, KB_T, MU_MESH, ALPHA_G_PROD, float(T_G_BITS))
                )
                mu_raw += p_pt
                mu_2c += p_pt * max(0.0, 1.0 - get_mf(bf) - get_mg(bg))
        if not ng_ok_ea:
            continue
        nf_rem = _strip_fb(nf)
        if nf_rem > T_F:
            continue
        surv += 1
        if (nf, ng) in seen_raw:
            continue
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
        if is_coprime is None:
            is_coprime = _is_oh_coprime(aa, bb, mul_z, norm_q15)
        if not is_coprime:
            continue
        seen_raw.add((nf, ng))
        m_raw += 1
        if not _hits_missing_fb(aa, bb, nf, ng, miss_f_by_p, miss_g_by_p):
            m_2core += 1
            if first_rel is None:
                first_rel = (aa, bb, nf, ng, lpf, lpg)
                if first_only:
                    return m_raw, m_2core, surv, mu_raw, mu_2c, first_rel
    assert 0 <= m_2core <= m_raw <= surv
    return m_raw, m_2core, surv, mu_raw, mu_2c, first_rel


def verify_q15_exact_deg72_h5(b85_path, V_mod, DS_mod=None, full_replays=False):
    """Verify the Authentic Current-Pair Exact-Degree deg-(72, 24) u = -5 fixed-2-window
    EA-gated Q(zeta_15) H5 / H5-W campaign + independent production DS.ss_lower2_tau mu_pred
    calibration + live bit-exact full-ball replays (2 replays default <= 0.6s Docker; 6 with
    full_replays=True)."""
    t0 = time.perf_counter()
    if DS_mod is None:
        DS_mod = getattr(V_mod, "DS", None)
    if DS_mod is None:
        import q15_descent16 as DS_mod
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

    # Verify exact C_12 norm invariance and _fast_norm_q15 == V._norm_q15
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

    # 2. Load & verify b85 payload + all 108 archived missing 2-core FB ideal roots
    raw_txt = Path(b85_path).read_text() if b85_path is not None else _EMBEDDED_H5_DATA_B85
    data = json.loads(lzma.decompress(base64.b85decode("".join(raw_txt.split()))).decode())
    tw = data["tower"]
    assert tw["u_t"] == U_T and tw["p_t"] == P_T and tw["r_t"] == R_T
    assert tw["deg_f"] == 72 and tw["deg_g"] == 24 and tw["ball_pts"] == 2496
    assert tw["kA"] == KA_T and tw["kB"] == KB_T and tw["T_g_bits"] == T_G_BITS and tw["T_f_bits"] == T_F_BITS
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
    miss_fns = _miss_kappa_tables(DS_mod, miss_f_by_p, miss_g_by_p)

    # 3. Verify all 100 archived campaign ideals (roots, 100 disjoint C_12 orbits, usable => surv)
    recs = []
    for row in data["compact_recs"].split(";"):
        ps = row.split()
        assert len(ps) == 10
        recs.append([int(ps[0]), int(ps[1]), int(ps[2]), int(ps[3]), int(ps[4]),
                     int(ps[5]), int(ps[6]), int(ps[7]), float(ps[8]), float(ps[9])])
    assert len(recs) == 100
    rec_by_ideal = {}
    seen_orbs = set()
    for r in recs:
        side_i, q, s, t, is_d1, m_r, m_2, sv, mu_r, mu_2 = r
        assert q > B_FB and is_d1 in (0, 1) and 0 <= m_2 <= m_r <= sv and mu_r >= mu_2 > 0.0
        _verify_ideal(side_i, q, s, t, isp)
        orb = min((q, pow(s, e, q), (t * pow(pow(s, e, q), 5 * k, q)) % q) for e in (1, 7, 4, 13) for k in range(3))
        seen_orbs.add((side_i, orb))
        rec_by_ideal[(side_i, q, s, t)] = (is_d1, m_r, m_2, sv, mu_r, mu_2)
    assert len(seen_orbs) == 100, f"expected 100 disjoint C_12 orbits, got {len(seen_orbs)}"

    # Cross-match ALL archived live_ctrl & live_d2 c1 rows against compact_recs unconditionally
    assert len(data["live_ctrl"]) == 2 and len(data["live_d2"]) == 2
    for cside, cq, cs, ct, er, e2, esv, emr, em2 in data["live_ctrl"]:
        assert rec_by_ideal.get((cside, cq, cs, ct)) == (0, er, e2, esv, emr, em2)
    for rec in data["live_d2"]:
        c1s, c1q, c1s0, c1t0, m1r, m12, sv1, mu1r, mu12 = rec[1:10]
        assert rec_by_ideal.get((c1s, c1q, c1s0, c1t0)) == (1, m1r, m12, sv1, mu1r, mu12)

    # 4. Live bit-exact replays + live DS.ss_lower2_tau mu_pred verification
    #    Default (full_replays=False): 1 complete 2-hop depth-2 chain (2 full-ball replays, <= 0.6s Docker).
    #    Full (full_replays=True): 2 controls + 2 two-hop depth-2 chains (6 full-ball replays, ~1.5s Docker).
    n_replayed = 0
    if full_replays:
        for cside, cq, cs, ct, er, e2, esv, emr, em2 in data["live_ctrl"]:
            mr, m2, sv, mur, mu2, fr = _scan_ideal(V_mod, DS_mod, cside, cq, cs, ct, miss_f_by_p, miss_g_by_p, first_only=False, miss_fns=miss_fns)
            assert (mr, m2, sv) == (er, e2, esv) and abs(mur - emr) <= 2e-4 and abs(mu2 - em2) <= 2e-4
            assert 0 < m2 <= mr <= sv and fr is not None and _strip_fb(fr[3]) <= T_G and _strip_fb(fr[2]) <= T_F
            n_replayed += 1

    d2_chain_summary = []
    chains_to_run = data["live_d2"] if full_replays else data["live_d2"][:1]
    for rec in chains_to_run:
        par, c1s, c1q, c1s0, c1t0, m1r, m12, sv1, mu1r, mu12, c2s, c2q, c2s0, c2t0, m2r, m22, sv2, mu2r, mu22 = rec
        ps, pq, ps0, pt0 = par
        assert B_FB < c2q < c1q < pq
        _verify_ideal(ps, pq, ps0, pt0, isp)
        _, _, _, _, _, fr0 = _scan_ideal(V_mod, DS_mod, ps, pq, ps0, pt0, miss_f_by_p, miss_g_by_p, first_only=True)
        assert fr0 is not None
        aa0, bb0, nf0_r, ng0_r, lpf0, lpg0 = fr0
        assert _strip_fb(ng0_r) <= T_G and _strip_fb(nf0_r) <= T_F and not _hits_missing_fb(aa0, bb0, nf0_r, ng0_r, miss_f_by_p, miss_g_by_p)
        n1_val = nf0_r if c1s == 1 else ng0_r
        assert c1q in (lpf0 if c1s == 1 else lpg0) and n1_val % c1q == 0 and (n1_val // c1q) % c1q != 0
        ch1 = _find_deg1_child(c1s, c1q, aa0, bb0, isp)
        assert ch1 == (c1s, c1q, c1s0, c1t0)
        r1r, r12, rsv1, rmu1r, rmu12, fr1 = _scan_ideal(V_mod, DS_mod, ch1[0], ch1[1], ch1[2], ch1[3], miss_f_by_p, miss_g_by_p, first_only=False, miss_fns=miss_fns)
        assert (r1r, r12, rsv1) == (m1r, m12, sv1) and abs(rmu1r - mu1r) <= 2e-4 and abs(rmu12 - mu12) <= 2e-4
        assert 0 < r12 <= r1r <= rsv1 and fr1 is not None
        n_replayed += 1

        aa1, bb1, nf1_r, ng1_r, lpf1, lpg1 = fr1
        assert _strip_fb(ng1_r) <= T_G and _strip_fb(nf1_r) <= T_F and not _hits_missing_fb(aa1, bb1, nf1_r, ng1_r, miss_f_by_p, miss_g_by_p)
        n2_val = nf1_r if c2s == 1 else ng1_r
        assert c2q in (lpf1 if c2s == 1 else lpg1) and n2_val % c2q == 0 and (n2_val // c2q) % c2q != 0
        ch2 = _find_deg1_child(c2s, c2q, aa1, bb1, isp)
        assert ch2 == (c2s, c2q, c2s0, c2t0)
        r2r, r22, rsv2, rmu2r, rmu22, fr2 = _scan_ideal(V_mod, DS_mod, ch2[0], ch2[1], ch2[2], ch2[3], miss_f_by_p, miss_g_by_p, first_only=False, miss_fns=miss_fns)
        assert (r2r, r22, rsv2) == (m2r, m22, sv2) and abs(rmu2r - mu2r) <= 2e-4 and abs(rmu22 - mu22) <= 2e-4
        assert 0 < r22 <= r2r <= rsv2 and fr2 is not None and _strip_fb(fr2[3]) <= T_G and _strip_fb(fr2[2]) <= T_F
        n_replayed += 1
        d2_chain_summary.append((pq, c1q, c2q, r12, r22))
    assert n_replayed == (6 if full_replays else 2)

    # 4b. Live bit-exact STRICT-DECREASE (ALL primes on BOTH f and g sides strictly < q_parent,
    #     ZERO companion primes >= q_parent!) exact-degree deg-(72, 24) u = -5 witnesses
    _STRICT_U5_WITNESSES = (
        (1, 423931, 72377, 186995, 7, 251431),
        (1, 348811, 87983, 83727, 1013, 271651),
        (1, 3245611, 197891, 740795, 39, 1969831),
        (1, 3455821, 553313, 1115103, 84, 3113221),
    )
    def _factor_full_lt(n, bound):
        rem = n
        fs = []
        d = 2
        while d * d <= rem and d <= 4096:
            while rem % d == 0:
                fs.append(d)
                rem //= d
            d += 1 if d == 2 else 2
        stk = [rem] if rem > 1 else []
        while stk:
            m = stk.pop()
            if m == 1:
                continue
            if isp(m):
                assert m < bound, f"strict decrease violated: prime {m} >= q={bound}"
                fs.append(m)
            else:
                g = _brent(m)
                assert 1 < g < m, f"brent failed on {m}"
                stk.append(g)
                stk.append(m // g)
        prod_chk = 1
        for p in fs:
            prod_chk *= p
        assert prod_chk == n
        return sorted(fs)

    strict_rows = []
    for sw_side, sw_q, sw_s, sw_t, sw_idx, sw_exp_max in _STRICT_U5_WITNESSES:
        _verify_ideal(sw_side, sw_q, sw_s, sw_t, isp)
        B_sw, _, _ = V_mod._lll_16d_qcm(V_mod._sq_hnf(sw_q, sw_s, sw_t))
        comb = _B16_3[sw_idx]
        v = [0] * 16
        for bi, sgn in comb:
            row = B_sw[bi]
            if sgn == 1:
                for j in range(16):
                    v[j] += row[j]
            else:
                for j in range(16):
                    v[j] -= row[j]
        aa, bb = v[:8], v[8:]
        assert _is_oh_coprime(aa, bb, mul_z, norm_q15)
        nf_sw, ng_sw = _nfg(aa, bb)
        assert nf_sw > 0 and ng_sw > 0 and nf_sw % sw_q == 0 and (nf_sw // sw_q) % sw_q != 0
        ff_sw = _factor_full_lt(nf_sw // sw_q, sw_q)
        fg_sw = _factor_full_lt(ng_sw, sw_q)
        max_lp_sw = max(ff_sw[-1] if ff_sw else 1, fg_sw[-1] if fg_sw else 1)
        assert max_lp_sw == sw_exp_max and max_lp_sw < sw_q
        strict_rows.append(dict(
            q=sw_q, idx=sw_idx, max_lp=max_lp_sw,
            ratio_max_lp_over_q=round(max_lp_sw / sw_q, 4),
            n_primes_f_gt_fb=sum(1 for p in ff_sw if p > B_FB),
            n_primes_g_gt_fb=sum(1 for p in fg_sw if p > B_FB),
        ))
    strict_dec = dict(
        ok=True,
        strict_max_lp_lt_q=True,
        zero_companions_ge_q=True,
        n_strict_witnesses=len(strict_rows),
        witnesses=strict_rows,
    )

    # 5. Summary stats + INDEPENDENT production-estimator calibration across all 100 ideals
    arms = {
        "f_ctrl": [r for r in recs if r[0] == 1 and r[4] == 0],
        "f_d1": [r for r in recs if r[0] == 1 and r[4] == 1],
        "g_ctrl": [r for r in recs if r[0] == 0 and r[4] == 0],
        "g_d1": [r for r in recs if r[0] == 0 and r[4] == 1],
    }
    for nm, rs in arms.items():
        assert len(rs) == 25, f"arm {nm} count {len(rs)} != 25"

    def _stats(rows, idx, mu_idx):
        xs = [r[idx] for r in rows]
        mus = [r[mu_idx] for r in rows]
        lqs = [math.log2(r[1]) for r in rows]
        n = len(xs)
        mn = sum(xs) / n
        mu_mn = sum(mus) / n
        vr = sum((x - mn) ** 2 for x in xs) / (n - 1)
        # Residual variance controlling for independent per-ideal mu_pred_i (eliminates q-size drift!)
        res_vr = sum((x - (mn / mu_mn) * m) ** 2 for x, m in zip(xs, mus)) / (n - 1)
        mlq = sum(lqs) / n
        slope = sum((a - mlq) * (b - mn) for a, b in zip(lqs, xs)) / max(1e-12, sum((a - mlq) ** 2 for a in lqs))
        q_vr = sum((b - (mn + slope * (a - mlq))) ** 2 for a, b in zip(lqs, xs)) / (n - 2)
        ind_p0 = sum(max(0.001, 1.0 / (1.0 + 0.5 * m)) for m in mus)
        return dict(
            n=n, min=min(xs), max=max(xs), mean=round(mn, 2),
            mu_pred_mean=round(mu_mn, 2),
            obs_over_mu_pred=round(mn / mu_mn, 4),
            ind_p0_sum=round(ind_p0, 2),
            fano=round(vr / mn, 3), k_hat=round(max(0.0, (vr - mn) / (mn * mn)), 4),
            fano_mu_adj=round(res_vr / mn, 3),
            fano_qadj=round(q_vr / mn, 3), k_hat_qadj=round(max(0.0, (q_vr - mn) / (mn * mn)), 4),
            n_zero=sum(1 for x in xs if x == 0),
        )

    raw_s = {k: _stats(v, 5, 8) for k, v in arms.items()}
    c2_s = {k: _stats(v, 6, 9) for k, v in arms.items()}
    n_zero_raw = sum(s["n_zero"] for s in raw_s.values())
    n_zero_2c = sum(s["n_zero"] for s in c2_s.values())
    ind_p0_raw = sum(max(0.001, 1.0 / (1.0 + 0.5 * r[8])) for r in recs)
    ind_p0_2c = sum(max(0.001, 1.0 / (1.0 + 0.5 * r[9])) for r in recs)
    pooled_cal_raw = sum(r[5] for r in recs) / sum(r[8] for r in recs)
    pooled_cal_2c = sum(r[6] for r in recs) / sum(r[9] for r in recs)
    # Assert every arm has obs_over_mu_pred >= 0.90, n_zero < ind_p0_sum, and mu-adjusted Fano <= 1.60!
    assert all(s["obs_over_mu_pred"] >= 0.90 and s["n_zero"] < s["ind_p0_sum"] and s["fano_mu_adj"] <= 1.60 for s in raw_s.values())
    assert all(s["obs_over_mu_pred"] >= 0.90 and s["n_zero"] < s["ind_p0_sum"] and s["fano_mu_adj"] <= 1.60 for s in c2_s.values())
    assert pooled_cal_raw >= 1.00 and pooled_cal_2c >= 1.00
    assert n_zero_raw * 3.0 < ind_p0_raw and n_zero_2c * 3.0 < ind_p0_2c

    f_raw_tot = sum(r[5] for r in arms["f_ctrl"] + arms["f_d1"])
    g_raw_tot = sum(r[5] for r in arms["g_ctrl"] + arms["g_d1"])
    f_2c_tot = sum(r[6] for r in arms["f_ctrl"] + arms["f_d1"])
    g_2c_tot = sum(r[6] for r in arms["g_ctrl"] + arms["g_d1"])
    surv_ctrl = sum(r[7] for r in recs if r[4] == 0) / 50.0
    surv_d1 = sum(r[7] for r in recs if r[4] == 1) / 50.0
    surv_ratio = surv_d1 / surv_ctrl
    assert surv_ratio <= 1.08, f"H5-W survivor ratio {surv_ratio:.4f} > 1.08"

    return dict(
        protocol="fixed_2window_ea_gated_kB1_deg72_24",
        usable_implies_surv=True,
        full_replays=bool(full_replays),
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
        independent_h5_p0_raw=round(ind_p0_raw, 2),
        independent_h5_p0_2core=round(ind_p0_2c, 2),
        pooled_obs_over_mu_pred_raw=round(pooled_cal_raw, 4),
        pooled_obs_over_mu_pred_2core=round(pooled_cal_2c, 4),
        fb_retention_f=round(f_2c_tot / f_raw_tot, 4),
        fb_retention_g=round(g_2c_tot / g_raw_tot, 4),
        f_d1_over_ctrl_raw=round(raw_s["f_d1"]["mean"] / raw_s["f_ctrl"]["mean"], 4),
        g_d1_over_ctrl_raw=round(raw_s["g_d1"]["mean"] / raw_s["g_ctrl"]["mean"], 4),
        f_d1_over_ctrl_2core=round(c2_s["f_d1"]["mean"] / c2_s["f_ctrl"]["mean"], 4),
        g_d1_over_ctrl_2core=round(c2_s["g_d1"]["mean"] / c2_s["g_ctrl"]["mean"], 4),
        surv_ctrl_mean=round(surv_ctrl, 2),
        surv_d1_mean=round(surv_d1, 2),
        surv_d1_over_ctrl=round(surv_ratio, 4),
        strict_decrease_deg72_24=strict_dec,
        wall_s=round(time.perf_counter() - t0, 4),
    )


if __name__ == "__main__":
    import sys
    here = Path(__file__).resolve().parent
    args = [a for a in sys.argv[1:] if a != "--full"]
    full = "--full" in sys.argv[1:]
    pkg_dir = Path(args[0]) if args else here
    sys.path.insert(0, str(pkg_dir))
    import q15_descent16 as DS
    import verify_q15_descent as V
    out = verify_q15_exact_deg72_h5(None, V, DS, full_replays=full)
    print(json.dumps(out, indent=2))


_EMBEDDED_H5_DATA_B85 = """
{Wp48S^xk9=GL@E0stWa8~^|S5YJf5;25L{H(dZdh>C|aR|CXjM(rJQe>eg`*kK;RKW%*(W%qB0Xo8;8WNkJFpzI&K>@_KILs)IB
eIObNpxy|5jjQHCEa23+Qv25+Ae$tcA%{b@{1-I8r?ru~$y#ia<qC9oK}euYA<*uEe_w+F2Qh9Y1bgy@eU8*ep}J*ijOUV-m9njn
biGrSA@;WbhXSM+&2Y1O+);K4M;2jm6TP&zSFY>CR$;80$yC0w_nB2t?@jj;B)VC6b78`FPsi+6aPZ0`;6s6K5Y<oh>rHe+bOM0w
7^rk$<P1^|i!KL8<n@Uss#87zIlrt++M47Y@}nFWz%N`*sio<r-$o^;4XB`57lJe-zxWgADUP_Tlbeohb!~618S~f~3Fu<zLVOur
MY2DA1YI$T7CkgCv1J=28K&sxzn~G9PT_(2L;qaOgNBNK8z^?G+!c?N(1TRl%JGD3(v}hxC8%UNRJ^B;;Tkd?z#!(&dSyo(l$6(p
lvChws{yJ`pOM9%aw|SK!VeAqzTQU>dbjw#M014Dy7GWZ#*Q^iE%W5l#rtqGv^q{oV;KKE!0~e;s~t4#z%kr{*d@CiS&d;w!wGlU
zSqmIE7zNC6rPrce<w)mcMHKhdmU|Tq==`MN5jC^PTeGy8nqB96Q}7d>fcDmLxP=GHXrc9N|LDpVziQ=%$KDp=fc!Y>QaRs&||PB
Qm>kkw9AdyRz6|r$I6kX1RVH(e?MQYu(RE6p?d+9k^sahK^UdCh+RXe+aH#&^m@-!rHz@}!Ticb*;3Gkcd5`OheFfJahFE|D;2)s
8;&c|0q8}Yo6JDnvbg6a+vu3o=}K9hKzP4K;Z~`<!SCB;f}4C`E!<&l++ik!ki7WDT3)*G2${9kWw@rlKDV7%=lmvF5W_)w0bHq(
KC}D+X3sR}dFyKd1u7P{Fo7rf1NP^-5qYYMuHp`!ri0UO9RqdbIM$i(SisL0uI0(0TY|;vX~P16Pd6wa3J^`h)%Zsu4GgHQxe&(a
U)~Hsk?-bVDyLgdH@&=MNV2IkIU^R)wIKS%wgyk3^sf!>5+pDOY<_Q(EY)O>mGjE+ATOeEl3yN(iamE8-)-|cJZi_O8lT;F+S)z>
InAize&ON+gn+7%&RY&Dl)p=KktZXL{Uof2glfNZLm*TvNm}`7iw{JK-S<@xZ&!+5@xcVVbsWWBGFD;EgG#mvVG1cF-o$w1)59Gi
==rpvS9Y$nxeWmCqK_L|U&#tRv6mG`I=AdhBFRbv1tVk*jGfJvmI5?AqsI%E54R8j=-#S>LrKY@L^a6PCax+JKSNYF?~N4szObF9
?=C^Ie-LBBj0Z_wnXZv>nv%5ZGstLAY=pnMa}Ny_lU;2V%W;@WNl#GXQ#Tqjm<B?Fu5~cW&e5TJmx-E;Y)?sZ*#L?by!PM~*(Sx<
2|k*e=5)}QJ?{0zmy@j@gWlL)i%&WU&8Hw2a=IswEQGZPLz!4~qS=-><Df@aLr^5cccf29{Fcdu{>wKtXH6#=W-B?)>rer?P<h|F
h#IKrQJ1X*x;J4T9Q)`|Eo@rAF;DmkZ!}rVZ|u4APQGT@^if?l+u`d!$qb2zUczs*%PkQ_Zav+#UJUiZR!MedvTlYAGKq$qOd3$*
oolf+M7?ibg1@)=xZz8nKAq#2kq7Bcn^7Y}m!}84LB?t^p6$|#-3`_w+KXtlL8l6^-|@4c7FrDFTmN<!Y6CS-XBr-SZVjz3Hiw46
4OAKEi|_cC{cauwV%;n`G-Bx)V5c~kQYVb*UxQ;Uet{ks=KG7ea{BgV(a6vQk>0%v{_0NhW6g3&3NwBiZ!}CxgA6Uh3*kRmIf*dp
{7KwAgEk^Csf*Z2k?!JW15-@4X{N(g7N(6f3M#EBf$m_96#B0nK*rW9nW$SPRtGN!(u!koLO8wr!iU=uK4{}|Mw{LW51u=dh7-m4
&3m17#3@q9o1W~T!yzj(jggWNVRAkUf@npFrzXLX8|)2Vi2rP8Dz!$1wHd1&3gcXEzi{s{Pg1O(7HJ<j7E+QlkO3w!_MR9L*9g%q
ft|K5B1&?^+-788oGQbBYQ|{+7M9%R0sgEl3x=%$=6pM7$x!l1nU<=R8KW7Mj%sYE+P2e-7C-ueA6rGrkc=^QZ)-en);L9MiE2@X
U*nf7elqOO^O@VJo3P%RN=qP%1)W_(h9qjhh{Oz6)s-i@l%f(8r`z}Hi@MB*00lT=><#@JVu0o!<h1c*Xq`mJ&k7l+OSeR$15#Qm
QMQ|1coK9l6G`LTjqx;^YN<|48P0nN(BBlR*;ljiS+TuISjgM<1JO%7ALe}j)e)d9$FtTF(N($(cDHs*^7~E68rz=l%o_<%-hMj%
`R{vfj%}oDJG26R*sN~0iOz`3=n?c4uKoY8i=_4sE`+K#)c`rP46p~ZUXNEjD;bMV0$(1D@t>@5cS|+15z&5YA+|kAqd(mK3;u<0
t$FH11H7Ez!XX{>b2vs4L%*W&Rth6h<fv`d&~-~Fh{J9@TP4Q79-|9u$9dkJ&R44Xrq4i|6*O1!X41M?+Bb*bg4c41oP8CJgJiR5
XKsU^nso_X-?>(&Y{%k%C}9`(@?LE{%E)eWtqdH)h{R6(<d1ndY>F0|ss|tNFEOOmH>%J=tD6bFAH1UBX|3>T43~xz#@gB}WPM<Y
>a4B~2ID{TNF%(kgYTfm#O@U|qYQqy8x~EZBXH$`I76j_5=LRk23>P%J2#{z=nVwsj9Wl-=LvST9H7I|n(1ZBw-L&wWF>%tqj|*4
IC3|}<{VtEMTgKr@xa#iKsoKRUj}k)TIEL%1uJ(2X0L**U%E`lmRk;;za8ITAC1hNcTo$0{@@0rcqE-eBkXzct0f6;vcuoA-bbo-
xIgJ&&;<TT!WxK3WtW2_!s%sj1<xZy^sBY|kkyb_8r2(@b?3&CmE7wXa5LyP5ucLbNWzbNrxt-Lr77e(zX@tZco9{j35RE5%vlR?
Y{ohPC3+e2Rf&;3pJL9rur}rD!CQDPH)apw-4Z<m=k@{Ni+_5pkcmfY+!CWV{y1?ckXW(1e!-~H+;3i*f4NJEP8a|H1h<zfY?lzn
00Gkyr7-{iE*QfSvBYQl0ssI200dcD
"""
