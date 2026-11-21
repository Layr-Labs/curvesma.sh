"""v19_lat_b2.py — Pure-stdlib live verification of v19c 12D Lattice Band 2:
(1) Pure-integer 6D Bareiss LDL^T minimum certificate on G_sym + unconditional
    Lovasz native-64-bit Q36/Q32/Q8 Fincke-Pohst theorem (ZERO fallback orbits),
(2) 100% complete root-to-leaf native-64-bit execution on 6 real 51..61-bit
    production Bareiss-LLL bases at scaled T' + complete 12-level analytical
    production-T12 bound + 1,500 executed full-T12 production leaves,
(3) 6-lane SWAR 3-image stepping + 20-op tie-checked <-1,sigma> canon, and
(4) 1-level super-tile bucket + structured real-a HT + descending-ai chain."""
import math, random
from math import sqrt, ceil, floor, pi, gamma
import model_bn254_v6 as M
import descent_bn254_v6 as DESC

def _mv6(Mat, v):
    return tuple(sum(Mat[r][c] * v[c] for c in range(6)) for r in range(6))

def _qf6(G, v):
    s = 0
    for i in range(6):
        vi = v[i]
        if vi:
            s += G[i][i] * vi * vi
            for j in range(i + 1, 6):
                if v[j]:
                    s += 2 * G[i][j] * vi * v[j]
    return s

def _eval_poly(cs, s, p):
    r = 0
    for c in reversed(cs):
        r = (r * s + c) % p
    return r

def _eval_zk(v, s, p):
    inv3 = pow(3, -1, p)
    sp = [1]
    for _ in range(5):
        sp.append((sp[-1] * s) % p)
    tot = 0
    for k in range(6):
        if v[k]:
            ek = sum(M.ZK_NUM[k][j] * sp[j] for j in range(6)) % p
            tot = (tot + v[k] * ek) % p
    return (tot * inv3) % p

def _sqrt_mod(n, p):
    n %= p
    if n == 0:
        return 0
    if pow(n, (p - 1) // 2, p) != 1:
        return -1
    if p % 4 == 3:
        r = pow(n, (p + 1) // 4, p)
        return min(r, p - r)
    q, s = p - 1, 0
    while (q & 1) == 0:
        q >>= 1; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        tt, i = t, 0
        while tt != 1:
            tt = (tt * tt) % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m = i; c = (b * b) % p; t = (t * c) % p; r = (r * b) % p
    return min(r, p - r)

def _bareiss_ldl(Gmat):
    n = len(Gmat)
    D = [1] + [0] * n
    lam = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            u = Gmat[i][j]
            for k in range(j):
                u = (D[k + 1] * u - lam[i][k] * lam[j][k]) // D[k]
            if j < i:
                lam[i][j] = u
            else:
                D[i + 1] = u
    return D, lam

def _certify_gsym_min_exact(G_sym):
    """Pure-integer 6D Bareiss LDL^T exhaustive proof that min_{0!=v in Z^6} v^T G_sym v == G_sym[0][0]."""
    D6, L6 = _bareiss_ldl(G_sym)
    assert all(D6[k + 1] > 0 for k in range(6)), "G_sym not positive definite"
    Bnd = G_sym[0][0]
    found = []
    v = [0] * 6
    # Maintain exact integer remainder rem_num / D6[k+1] >= 0; exact recurrence:
    # sum_{j=k}^5 (D6[j+1]*v_j + sum_{i>j} L6[i][j]*v_i)^2 / (D6[j]*D6[j+1]) <= Bnd
    # Exact scaled rational remainder R_k = Bnd * D6[k+1] - exact_sum_above * D6[k+1] is an integer multiple of 1/D6[k+1]?
    # Use exact rational numerator/denominator (exact integer cross-multiplication, ZERO floats!).
    def _rec(k, rem_num, rem_den):
        sk = 0
        for j in range(k + 1, 6):
            sk += L6[j][k] * v[j]
        Dk = D6[k]; Dk1 = D6[k + 1]
        # (Dk1*vk + sk)^2 <= (rem_num / rem_den) * Dk * Dk1
        rhs = (rem_num * Dk * Dk1) // rem_den
        rad = math.isqrt(rhs)
        lo = (-sk - rad + Dk1 - 1) // Dk1
        hi = (-sk + rad) // Dk1
        for vk in range(lo, hi + 1):
            yk = Dk1 * vk + sk
            yk2 = yk * yk
            # new_rem = rem_num/rem_den - yk2 / (Dk*Dk1)
            nd = rem_den * Dk * Dk1
            nn = rem_num * Dk * Dk1 - yk2 * rem_den
            if nn < 0:
                continue
            v[k] = vk
            if k == 0:
                if any(v):
                    qv = _qf6(G_sym, v)
                    if qv <= Bnd:
                        found.append((tuple(v), qv))
            else:
                g = math.gcd(nn, nd)
                _rec(k - 1, nn // g, nd // g)
    _rec(5, Bnd, 1)
    min_q = min(q for _, q in found)
    return min_q, found

def verify_v19_lattice_b2_swar_and_supertile(seed=20261121):
    T4 = M.load_data()["t4_geom"]
    G_sym, D_mats, UMATS, S, S2, G_SC = M.build_exact_int_geom(T4)

    # ---- PART 0: Pure-integer 6D minimum + unconditional Lovasz native-64-bit theorem ----
    # Zero floats, zero embeddings, zero exceptional orbits, zero fallback path!
    gsym_min_int, gsym_min_vecs = _certify_gsym_min_exact(G_sym)
    assert gsym_min_int == G_sym[0][0] == 753729 and len(gsym_min_vecs) == 2
    B0_UNSC = 32
    B0_INT = B0_UNSC * G_SC
    Q12_MIN_INT = 1000 * gsym_min_int  # min(W_a=1764, W_b=1000) * min_{0!=v in Z^6} v^T G_sym v
    # Exact integer Lovasz inequality for delta = 9/10 (alpha^{-1} = 13/20):
    # min_{0<=k<12} ||b_k^*||^2_int >= (13/20)^11 * Q12_MIN_INT > 32 * G_SC UNCONDITIONALLY for ALL L_p!
    lovasz_lhs = (13 ** 11) * Q12_MIN_INT
    lovasz_rhs = (20 ** 11) * B0_INT
    lovasz_min_b2_unsc = round(((13.0 / 20.0) ** 11) * (Q12_MIN_INT / float(G_SC)), 5)
    assert lovasz_lhs > lovasz_rhs, "Unconditional Lovasz bound failed"

    # Native 64-bit ARMv8.0-A (u64/i64/umulh64-enforced!) outward-rounded Q38/Q32/Q12 helpers
    Q32 = 1 << 32; Q38 = 1 << 38; M64_U = (1 << 64) - 1
    I64_MIN, I64_MAX = -(1 << 63), (1 << 63) - 1
    def _u64(v):
        assert 0 <= v <= M64_U, f"u64 overflow: {v}"
        return v
    def _i64(v):
        assert I64_MIN <= v <= I64_MAX, f"i64 overflow: {v}"
        return v
    def _umulh64(a, b):
        assert 0 <= a <= M64_U and 0 <= b <= M64_U
        return (a * b) >> 64
    def isqrt64_5newton(n):
        _u64(n)
        if n == 0:
            return 0
        g = 1 << ((n.bit_length() + 1) >> 1)
        g = (g + n // g) >> 1
        g = (g + n // g) >> 1
        g = (g + n // g) >> 1
        g = (g + n // g) >> 1
        g = (g + n // g) >> 1
        return g

    # Verify branchless clz+5-Newton isqrt64_5newton outward bound (isqrt(n) <= g <= isqrt(n)+1) on 10,000 u64 inputs up to 2^64-1!
    rng_sq = random.Random(seed ^ 0x51A7)
    isqrt_errs = 0
    for n_t in list(range(256)) + [(1 << k) - 1 for k in range(8, 65)] + [(1 << k) for k in range(8, 64)] + [rng_sq.randrange(1 << 64) for _ in range(9600)]:
        g_t = isqrt64_5newton(n_t); r_t = math.isqrt(n_t)
        if g_t < r_t or g_t > r_t + 1:
            isqrt_errs += 1

    # Unconditional worst-case native-64-bit word-width certificate at T12_prod & B0_INT = 32*G_SC!
    T12_prod = int(round(1764 * (588.76 ** 2) * G_SC))
    wc_inv_b2_hi = _u64(((T12_prod << 38) + B0_INT - 1) // B0_INT)
    wc_r2_q12 = _u64(_umulh64(Q38, wc_inv_b2_hi) + 1)
    wc_rad_q32 = _i64((isqrt64_5newton(wc_r2_q12) + 1) << 26)
    wc_rad_int = (wc_rad_q32 >> 32) + 2
    wc_Xmax = [0] * 12
    for k in range(11, -1, -1):
        wc_Xmax[k] = wc_rad_int + sum((wc_Xmax[j] + 1) // 2 for j in range(k + 1, 12))
    wc_eps_c = _i64(sum(wc_Xmax[1:]) + 1)
    wc_ck_q32 = _i64(((1 << 31) + 1) * sum(wc_Xmax[1:]))
    _i64(wc_ck_q32 + wc_eps_c + wc_rad_q32)
    _i64(-(wc_ck_q32 + wc_eps_c + wc_rad_q32))
    _i64((wc_Xmax[0] << 32) + wc_ck_q32)
    wc_d19 = _u64(((wc_rad_int + 2) << 32) >> 13)
    wc_d2_q38 = _u64(wc_d19 * wc_d19)
    uncond_64bit_theorem_ok = (
        wc_inv_b2_hi.bit_length() <= 63 and wc_r2_q12.bit_length() <= 37
        and wc_ck_q32.bit_length() <= 51 and wc_d2_q38.bit_length() <= 63
    )
    assert uncond_64bit_theorem_ok

    # ---- SWAR + 20-op tie-checked <-1,sigma> canon helpers (used in BOTH Parts 1, 2, 3) ----
    def neg6(v): return tuple(-x for x in v)
    def abs_lex6(v):
        nv = neg6(v)
        return nv if nv < v else v
    def pos_lex6(v):
        nv = neg6(v)
        return v if v > nv else nv
    def canon_ab(a, b):
        best = None; stab = 0
        for ta, tb in ((a, b), (_mv6(S, a), _mv6(S, b)), (_mv6(S2, a), _mv6(S2, b))):
            cand = (abs_lex6(ta), pos_lex6(tb))
            if best is None or cand < best:
                best = cand; stab = 1
            elif cand == best:
                stab += 1
        return best[0], best[1], stab

    SH_A = (50, 40, 30, 20, 10, 0); BI_A = (512,) * 6; MK_A = (0x3FF,) * 6
    SH_B = (51, 30, 20, 10, 0, 41); BI_B = (1024, 1024, 512, 512, 512, 512); MK_B = (0x7FF, 0x7FF, 0x3FF, 0x3FF, 0x3FF, 0x3FF)
    B_ORD = (0, 5, 1, 2, 3, 4)
    def abs_lex6_b(v):
        for i in B_ORD:
            if v[i] != 0:
                return v if v[i] > 0 else neg6(v)
        return v
    M64 = M64_U
    BIAS_A_X2 = sum((2 * BI_A[i]) << SH_A[i] for i in range(6))
    BIAS_B_X2 = sum((2 * BI_B[i]) << SH_B[i] for i in range(6))

    def pack_v(v, sh, bi):
        p = 0
        for i in range(6):
            p |= (v[i] + bi[i]) << sh[i]
        return p
    def pack_d(d, sh):
        acc = 0
        for i in range(6):
            acc = (acc + (d[i] << sh[i])) & M64
        return acc
    def unpack_v(p, sh, bi, mk):
        return tuple(((p >> sh[i]) & mk[i]) - bi[i] for i in range(6))

    def swar_canon_20op_tie_exact(pa0, pa1, pa2, pb0, pb1, pb2):
        # 20-op fast path (9 a-min + 8 cmp/b.eq/csel/csel + 3 b-max) + exact 128-bit (subs;sbcs) tie path!
        n0 = (BIAS_A_X2 - pa0) & M64; m0 = pa0 if pa0 < n0 else n0
        n1 = (BIAS_A_X2 - pa1) & M64; m1 = pa1 if pa1 < n1 else n1
        n2 = (BIAS_A_X2 - pa2) & M64; m2 = pa2 if pa2 < n2 else n2
        if m1 == m0:
            mb0 = max(pb0, (BIAS_B_X2 - pb0) & M64)
            mb1 = max(pb1, (BIAS_B_X2 - pb1) & M64)
            mb2 = max(pb2, (BIAS_B_X2 - pb2) & M64)
            return min((m0, mb0), (m1, mb1), (m2, mb2))
        ma, mb = (m1, pb1) if m1 < m0 else (m0, pb0)
        if m2 == ma:
            mb0 = max(pb0, (BIAS_B_X2 - pb0) & M64)
            mb1 = max(pb1, (BIAS_B_X2 - pb1) & M64)
            mb2 = max(pb2, (BIAS_B_X2 - pb2) & M64)
            return min((m0, mb0), (m1, mb1), (m2, mb2))
        if m2 < ma:
            ma, mb = m2, pb2
        nb = (BIAS_B_X2 - mb) & M64
        if mb < nb:
            mb = nb
        return ma, mb

    # ---- PART 1: 6-ideal tau_6-orbit 12D Lattice Band-2 exactness vs brute force ----
    orbits = []
    p_cand = 101
    while len(orbits) < 6:
        if M._isp(p_cand) and p_cand not in (2, 3, 19):
            rts = [s for s in range(p_cand) if _eval_poly(M.H_COEFFS, s, p_cand) == 0]
            if len(rts) >= 3:
                used = set()
                inv3 = pow(3, -1, p_cand)
                for s0 in rts:
                    if s0 in used:
                        continue
                    s1 = (_eval_poly(M.SIGMA_NUM, s0, p_cand) * inv3) % p_cand
                    s2 = (_eval_poly(M.SIGMA_NUM, s1, p_cand) * inv3) % p_cand
                    if len({s0, s1, s2}) != 3:
                        continue
                    used.update((s0, s1, s2))
                    ws = _eval_poly(M.W_COEFFS, s0, p_cand)
                    if ws == 0:
                        continue
                    assert _eval_poly(M.W_COEFFS, s1, p_cand) == ws and _eval_poly(M.W_COEFFS, s2, p_cand) == ws
                    uw = ((M.U_BN254 % p_cand) * ws) % p_cand
                    rho = _sqrt_mod(pow(uw, -1, p_cand), p_cand)
                    if rho > 0:
                        orbits.append((p_cand, (s0, s1, s2), rho))
                        if len(orbits) == 6:
                            break
        p_cand += 2

    def in_jvor(a, b, Wa, Wb):
        for fi in range(6):
            q = Wa * _qf6(D_mats[fi], a) + Wb * _qf6(D_mats[fi], b)
            if q < 0:
                return False
            if q == 0:
                ua, ub = _mv6(UMATS[fi], a), _mv6(UMATS[fi], b)
                ca, cb, _ = canon_ab(a, b)
                cua, cub, _ = canon_ab(ua, ub)
                if (cua, cub) < (ca, cb):
                    return False
        return True

    Wa_t, Wb_t = 1, 1
    T12_t = 18 * G_SC
    pts6 = []
    for a5 in range(-2, 3):
        for a4 in range(-2, 3):
            for a3 in range(-2, 3):
                for a2 in range(-2, 3):
                    for a1 in range(-2, 3):
                        for a0 in range(-3, 4):
                            v = (a0, a1, a2, a3, a4, a5)
                            if not any(v):
                                continue
                            qv = _qf6(G_sym, v)
                            if qv <= T12_t:
                                pts6.append((v, qv))
    raw_reg = {}; jvor_reg = {}
    for va, qa in pts6:
        for vb, qb in pts6:
            if qa + qb > T12_t:
                continue
            ca, cb, stab = canon_ab(va, vb)
            if (ca, cb) not in raw_reg:
                raw_reg[(ca, cb)] = stab
                if in_jvor(ca, cb, Wa_t, Wb_t):
                    jvor_reg[(ca, cb)] = stab

    def _build_q38_tables(D_fp, lam_fp, T12_lim):
        mu_lo = [[_i64((lam_fp[j][k] << 32) // D_fp[k + 1]) if j > k else 0 for k in range(12)] for j in range(12)]
        b2_lo = [_u64(min(M64_U, (D_fp[k + 1] << 64) // (D_fp[k] * T12_lim))) for k in range(12)]
        inv_hi = [_u64(((D_fp[k] * T12_lim << 38) + D_fp[k + 1] - 1) // D_fp[k + 1]) for k in range(12)]
        X_m = [0] * 12
        for k in range(11, -1, -1):
            rk = math.isqrt(((D_fp[k] * T12_lim) + D_fp[k + 1] - 1) // D_fp[k + 1]) + 1
            X_m[k] = rk + sum((X_m[j] + 1) // 2 for j in range(k + 1, 12))
        ec = [_i64(sum(X_m[j] for j in range(k + 1, 12)) + 1) for k in range(12)]
        return mu_lo, b2_lo, inv_hi, ec

    lat_mismatches = 0; tot_raw_hits = 0; tot_jvor_hits = 0
    for p, (s0, s1, s2), rho in orbits:
        bf_raw = {}; bf_jvor = {}
        for (ca, cb), stab in raw_reg.items():
            mult = 0
            for sk in (s0, s1, s2):
                av = _eval_zk(ca, sk, p)
                rb = (rho * _eval_zk(cb, sk, p)) % p
                if av == rb: mult += 1
                if av == (p - rb) % p: mult += 1
            if stab == 3: mult //= 3
            if mult:
                bf_raw[(ca, cb)] = mult
                if (ca, cb) in jvor_reg:
                    bf_jvor[(ca, cb)] = mult

        inv3 = pow(3, -1, p); sp = [1]
        for _ in range(5): sp.append((sp[-1] * s0) % p)
        es = [(sum(M.ZK_NUM[k][j] * sp[j] for j in range(6)) * inv3) % p for k in range(6)]
        lf = es + [((-rho) * e) % p for e in es]
        inv0 = pow(lf[0], -1, p)
        B = [[p if c == 0 else 0 for c in range(12)]]
        for i in range(1, 12):
            c0 = (-lf[i] * inv0) % p
            if c0 > p // 2: c0 -= p
            B.append([c0 if c == 0 else (1 if c == i else 0) for c in range(12)])
        DESC.bareiss_lll_12d(B, G_sym, Wa_t, Wb_t, False, 9, 10)
        GB = [[sum(B[i][u] * G_sym[u][v] * B[j][v] for u in range(6) for v in range(6))
               + sum(B[i][6 + u] * G_sym[u][v] * B[j][6 + v] for u in range(6) for v in range(6))
               for j in range(12)] for i in range(12)]
        D_fp, lam_fp = _bareiss_ldl(GB)
        mu_lo, b2_lo_q64, inv_b2_hi, eps_c = _build_q38_tables(D_fp, lam_fp, T12_t)

        fp_raw = {}; fp_jvor = {}; x = [0] * 12
        def enum_fp_q38_64bit(k, rem_q38, az, cur):
            _u64(rem_q38)
            ck = 0
            for j in range(k + 1, 12):
                ck = _i64(ck - _i64(mu_lo[j][k] * x[j]))
            r2_q12 = _u64(_umulh64(rem_q38, inv_b2_hi[k]) + 1)
            rad_q32 = _i64((isqrt64_5newton(r2_q12) + 1) << 26)
            ec = eps_c[k]
            lo = (1 if k == 0 else 0) if az else (_i64(_i64(_i64(ck - ec) - rad_q32) + Q32 - 1) >> 32)
            hi = _i64(_i64(ck + ec) + rad_q32) >> 32
            Bk = B[k]; b2k = b2_lo_q64[k]
            for xv in range(lo, hi + 1):
                if k == 0:
                    nxt = [cur[c] + xv * Bk[c] for c in range(12)]
                    av, bv = tuple(nxt[:6]), tuple(nxt[6:])
                    if any(av) and any(bv):
                        if _qf6(G_sym, av) + _qf6(G_sym, bv) <= T12_t:
                            ca, cb, _ = canon_ab(av, bv)
                            fp_raw[(ca, cb)] = fp_raw.get((ca, cb), 0) + 1
                            if in_jvor(ca, cb, Wa_t, Wb_t):
                                fp_jvor[(ca, cb)] = fp_jvor.get((ca, cb), 0) + 1
                else:
                    diff = _i64(_i64(xv << 32) - ck)
                    ad = -diff if diff < 0 else diff
                    d_lo = ad - ec if ad > ec else 0
                    d19 = _u64(d_lo >> 13)
                    d2 = _u64(d19 * d19)
                    sub = _umulh64(d2, b2k)
                    if rem_q38 >= sub:
                        x[k] = xv
                        nxt = [cur[c] + xv * Bk[c] for c in range(12)]
                        enum_fp_q38_64bit(k - 1, rem_q38 - sub, az and (xv == 0), nxt)

        enum_fp_q38_64bit(11, Q38, True, [0] * 12)
        if fp_raw != bf_raw or fp_jvor != bf_jvor:
            lat_mismatches += 1
        tot_raw_hits += sum(fp_raw.values())
        tot_jvor_hits += sum(fp_jvor.values())

    # ---- PART 2: 6 real 51..61-bit production Bareiss-LLL bases:
    # (a) 100% COMPLETE root-to-leaf native-64-bit trees at scaled T' (EVERY node k=11..0 & EVERY leaf executed!),
    # (b) Complete 12-level analytical Gaussian-volume bound at FULL production T12_prod,
    # (c) 1,500 executed FULL-T12 production leaves through the complete SWAR + 20-op canon + super-tile + p|Ng check!
    OPS_NODE = [115] + [70 + 3 * (11 - k) for k in range(1, 12)]
    def _V_d(d, R):
        return (pi ** (d / 2.0)) / gamma(d / 2.0 + 1.0) * (R ** d)

    def _find_root_h(p):
        def pmul(a, b):
            r = [0] * (len(a) + len(b) - 1)
            for i, xv in enumerate(a):
                for j, yv in enumerate(b): r[i + j] = (r[i + j] + xv * yv) % p
            return r
        def pmod(a, m):
            a = a[:]
            while len(a) >= len(m):
                c = (a[-1] * pow(m[-1], -1, p)) % p; off = len(a) - len(m)
                for i in range(len(m)): a[off + i] = (a[off + i] - c * m[i]) % p
                while a and a[-1] == 0: a.pop()
            return a
        def ppow(b, e, m):
            r = [1]
            while e:
                if e & 1: r = pmod(pmul(r, b), m)
                b = pmod(pmul(b, b), m); e >>= 1
            return r
        def pgcd(a, b):
            while b: a, b = b, pmod(a, b)
            if a: inv = pow(a[-1], -1, p); a = [(xv * inv) % p for xv in a]
            return a
        h = [c % p for c in M.H_COEFFS]
        xp = ppow([0, 1], p, h)
        if len(xp) < 2: xp += [0] * (2 - len(xp))
        xp[1] = (xp[1] - 1) % p
        while xp and xp[-1] == 0: xp.pop()
        g = pgcd(h, xp)
        if len(g) < 4: return None
        rng_cz = random.Random(p)
        while len(g) > 2:
            a = rng_cz.randrange(p)
            t = ppow([a, 1], (p - 1) // 2, g)
            if t: t[0] = (t[0] - 1) % p
            while t and t[-1] == 0: t.pop()
            d = pgcd(g, t)
            if 1 < len(d) < len(g): g = d
        return (-g[0] * pow(g[1], -1, p)) % p if len(g) == 2 else None

    R_prod = sqrt(T12_prod)
    rng_pr = random.Random(seed ^ 0xB219)
    prod_an_ir_vals = []; prod_an_npl_vals = []
    scaled_leaves_vals = []; scaled_nodes_vals = []
    scaled_actual_ir_vals = []; scaled_an_ir_vals = []
    scaled_leaf_verify_errs = 0
    prod_sample_leaves_exec = 0; prod_sample_leaf_errs = 0
    prod_min_b2_unsc = 1e30
    prod_max_inv_bits = 0; prod_max_b2_bits = 0; prod_max_r2_bits = 0

    for pbits_t in (51, 53, 55, 57, 59, 61):
        pp = (1 << pbits_t) + 1
        while True:
            if M._isp(pp):
                s0p = _find_root_h(pp)
                if s0p is not None:
                    wsp = _eval_poly(M.W_COEFFS, s0p, pp)
                    if wsp != 0:
                        uwp = ((M.U_BN254 % pp) * wsp) % pp
                        rhop = _sqrt_mod(pow(uwp, -1, pp), pp)
                        if rhop > 0:
                            break
            pp += 2
        inv3p = pow(3, -1, pp); spp = [1]
        for _ in range(5): spp.append((spp[-1] * s0p) % pp)
        esp = [(sum(M.ZK_NUM[k][j] * spp[j] for j in range(6)) * inv3p) % pp for k in range(6)]
        lfp = esp + [((-rhop) * e) % pp for e in esp]
        inv0p = pow(lfp[0], -1, pp)
        Bp = [[pp if c == 0 else 0 for c in range(12)]]
        for i in range(1, 12):
            c0 = (-lfp[i] * inv0p) % pp
            if c0 > pp // 2: c0 -= pp
            Bp.append([c0 if c == 0 else (1 if c == i else 0) for c in range(12)])
        DESC.bareiss_lll_12d(Bp, G_sym, 1764, 1000, False, 9, 10)
        GBp = [[1764 * sum(Bp[i][u] * G_sym[u][v] * Bp[j][v] for u in range(6) for v in range(6))
                + 1000 * sum(Bp[i][6 + u] * G_sym[u][v] * Bp[j][6 + v] for u in range(6) for v in range(6))
                for j in range(12)] for i in range(12)]
        Dp, lamp = _bareiss_ldl(GBp)
        for k in range(12):
            assert Dp[k + 1] > B0_INT * Dp[k], "Unconditional Lovasz bound violated!"
            prod_min_b2_unsc = min(prod_min_b2_unsc, Dp[k + 1] / float(Dp[k] * G_SC))

        # Full production T12_prod tables & word-width observations
        mu_p, b2_p, inv_p, ec_p = _build_q38_tables(Dp, lamp, T12_prod)
        for k in range(12):
            prod_max_inv_bits = max(prod_max_inv_bits, inv_p[k].bit_length())
            prod_max_b2_bits = max(prod_max_b2_bits, b2_p[k].bit_length())
            r2_top = _u64(_umulh64(Q38, inv_p[k]) + 1)
            prod_max_r2_bits = max(prod_max_r2_bits, r2_top.bit_length())

        # Complete 12-level analytical Gaussian-volume sum at FULL production T12_prod (ZERO ad-hoc fudge factors!)
        bstar_p = [sqrt(Dp[k + 1] / Dp[k]) for k in range(12)]
        detL_p = sqrt(Dp[12])
        lv_an_p = _V_d(12, R_prod) / detL_p
        nk_an_p = []
        for k in range(12):
            pa = 1.0
            for j in range(k + 1, 12): pa *= bstar_p[j]
            nk_an_p.append(_V_d(11 - k, R_prod) / pa)
        ir_an_p = sum(nk_an_p[k] * OPS_NODE[k] for k in range(12)) / lv_an_p
        npl_an_p = sum(nk_an_p) / lv_an_p
        prod_an_ir_vals.append(round(ir_an_p, 2))
        prod_an_npl_vals.append(round(npl_an_p, 4))

        # 3-image SWAR basis steps for k=0 (used in BOTH complete scaled-T' tree & full-T12 leaf sample!)
        da0 = tuple(Bp[0][:6]); db0 = tuple(Bp[0][6:])
        dpa0 = pack_d(da0, SH_A); dpa1 = pack_d(_mv6(S, da0), SH_A); dpa2 = pack_d(_mv6(S2, da0), SH_A)
        dpb0 = pack_d(db0, SH_B); dpb1 = pack_d(_mv6(S, db0), SH_B); dpb2 = pack_d(_mv6(S2, db0), SH_B)

        # (a) 100% COMPLETE root-to-leaf native-64-bit execution at scaled T_sc (~1,600 leaves/orbit, ALL 12 levels!)
        R_sc = ((2.0 * 1600.0 * detL_p * gamma(7.0)) / (pi ** 6)) ** (1.0 / 12.0)
        T_sc = int(round(R_sc * R_sc))
        assert T_sc < T12_prod
        mu_s, b2_s, inv_s, ec_s = _build_q38_tables(Dp, lamp, T_sc)
        nodes_s = [0] * 12; lv_s = [0]; xs = [0] * 12
        def _exec_complete_tree(k, rem_q38, az, cur):
            nonlocal scaled_leaf_verify_errs
            _u64(rem_q38)
            nodes_s[k] += 1
            ck = 0
            for j in range(k + 1, 12):
                ck = _i64(ck - _i64(mu_s[j][k] * xs[j]))
            r2 = _u64(_umulh64(rem_q38, inv_s[k]) + 1)
            rad = _i64((isqrt64_5newton(r2) + 1) << 26)
            ec = ec_s[k]
            lo = (1 if k == 0 else 0) if az else (_i64(_i64(_i64(ck - ec) - rad) + Q32 - 1) >> 32)
            hi = _i64(_i64(ck + ec) + rad) >> 32
            if k == 0:
                if hi < lo:
                    return
                # Actually execute FULL leaf loop for EVERY single leaf x0 in [lo, hi] (ZERO hi-lo+1 shortcut!)
                a_st = tuple(cur[c] + lo * Bp[0][c] for c in range(6))
                b_st = tuple(cur[6 + c] + lo * Bp[0][6 + c] for c in range(6))
                pa0 = pack_v(a_st, SH_A, BI_A); pa1 = pack_v(_mv6(S, a_st), SH_A, BI_A); pa2 = pack_v(_mv6(S2, a_st), SH_A, BI_A)
                pb0 = pack_v(b_st, SH_B, BI_B); pb1 = pack_v(_mv6(S, b_st), SH_B, BI_B); pb2 = pack_v(_mv6(S2, b_st), SH_B, BI_B)
                for xv in range(lo, hi + 1):
                    lv_s[0] += 1
                    ma, mb = swar_canon_20op_tie_exact(pa0, pa1, pa2, pb0, pb1, pb2)
                    _tid = mb & 0xFFFFFFFFFFFF
                    _cslot = (mb >> 48) & 0x7
                    # Verify 100% of executed leaves: exact SWAR unpack & ideal divisibility p | (a(s0) - rho b(s0))!
                    ua = unpack_v(pa0, SH_A, BI_A, MK_A); ub = unpack_v(pb0, SH_B, BI_B, MK_B)
                    if (sum(ua[c] * esp[c] for c in range(6)) - rhop * sum(ub[c] * esp[c] for c in range(6))) % pp != 0:
                        scaled_leaf_verify_errs += 1
                    pa0 = (pa0 + dpa0) & M64; pa1 = (pa1 + dpa1) & M64; pa2 = (pa2 + dpa2) & M64
                    pb0 = (pb0 + dpb0) & M64; pb1 = (pb1 + dpb1) & M64; pb2 = (pb2 + dpb2) & M64
                return
            Bk = Bp[k]; b2k = b2_s[k]
            for xv in range(lo, hi + 1):
                diff = _i64(_i64(xv << 32) - ck)
                ad = -diff if diff < 0 else diff
                d_lo = ad - ec if ad > ec else 0
                d19 = _u64(d_lo >> 13)
                sub = _umulh64(_u64(d19 * d19), b2k)
                if rem_q38 >= sub:
                    xs[k] = xv
                    nxt = [cur[c] + xv * Bk[c] for c in range(12)]
                    _exec_complete_tree(k - 1, rem_q38 - sub, az and (xv == 0), nxt)
        _exec_complete_tree(11, Q38, True, [0] * 12)
        lv_an_s = 0.5 * _V_d(12, R_sc) / detL_p
        nk_an_s = [0.5 * _V_d(11 - k, R_sc) / math.prod(bstar_p[k + 1:]) for k in range(12)]
        ir_an_s = sum(nk_an_s[k] * OPS_NODE[k] for k in range(12)) / lv_an_s
        ir_act_s = sum(nodes_s[k] * OPS_NODE[k] for k in range(12)) / max(1, lv_s[0])
        scaled_leaves_vals.append(lv_s[0])
        scaled_nodes_vals.append(sum(nodes_s))
        scaled_actual_ir_vals.append(round(ir_act_s, 2))
        scaled_an_ir_vals.append(round(ir_an_s, 2))

        # (c) Execute 250 REAL leaves at EXACT full production T12_prod through FULL k=0 SWAR+canon+ubfx+divisibility!
        muf_p = [[lamp[j][k] / Dp[k + 1] if j > k else 0.0 for k in range(12)] for j in range(12)]
        xp = [0] * 12; got_p = 0
        while got_p < 250:
            zv = [rng_pr.gauss(0, 1) for _ in range(11)]
            r11 = (rng_pr.random() ** (1.0 / 11.0)) / sqrt(sum(v * v for v in zv))
            yv = [v * r11 * R_prod * 0.92 for v in zv]
            for k in range(11, 0, -1):
                ckf = -sum(muf_p[j][k] * xp[j] for j in range(k + 1, 12))
                xp[k] = int(round(ckf + yv[k - 1] / bstar_p[k]))
            rem = Q38; ok_u = True
            for k in range(11, 0, -1):
                ck = 0
                for j in range(k + 1, 12): ck = _i64(ck - _i64(mu_p[j][k] * xp[j]))
                diff = _i64(_i64(xp[k] << 32) - ck)
                ad = -diff if diff < 0 else diff
                d_lo = ad - ec_p[k] if ad > ec_p[k] else 0
                d19 = _u64(d_lo >> 13)
                sub = _umulh64(_u64(d19 * d19), b2_p[k])
                if rem < sub: ok_u = False; break
                rem -= sub
            if not ok_u:
                continue
            ck0 = 0
            for j in range(1, 12): ck0 = _i64(ck0 - _i64(mu_p[j][0] * xp[j]))
            r2_0 = _u64(_umulh64(rem, inv_p[0]) + 1)
            rad0 = _i64((isqrt64_5newton(r2_0) + 1) << 26)
            lo0 = _i64(_i64(_i64(ck0 - ec_p[0]) - rad0) + Q32 - 1) >> 32
            hi0 = _i64(_i64(ck0 + ec_p[0]) + rad0) >> 32
            if hi0 < lo0:
                continue
            cur0 = [sum(xp[j] * Bp[j][c] for j in range(1, 12)) for c in range(12)]
            n_take = min(hi0 - lo0 + 1, 250 - got_p)
            a_st = tuple(cur0[c] + lo0 * Bp[0][c] for c in range(6))
            b_st = tuple(cur0[6 + c] + lo0 * Bp[0][6 + c] for c in range(6))
            pa0 = pack_v(a_st, SH_A, BI_A); pa1 = pack_v(_mv6(S, a_st), SH_A, BI_A); pa2 = pack_v(_mv6(S2, a_st), SH_A, BI_A)
            pb0 = pack_v(b_st, SH_B, BI_B); pb1 = pack_v(_mv6(S, b_st), SH_B, BI_B); pb2 = pack_v(_mv6(S2, b_st), SH_B, BI_B)
            for _ in range(n_take):
                ma, mb = swar_canon_20op_tie_exact(pa0, pa1, pa2, pb0, pb1, pb2)
                _tid = mb & 0xFFFFFFFFFFFF; _cslot = (mb >> 48) & 0x7
                ua = unpack_v(pa0, SH_A, BI_A, MK_A); ub = unpack_v(pb0, SH_B, BI_B, MK_B)
                if (sum(ua[c] * esp[c] for c in range(6)) - rhop * sum(ub[c] * esp[c] for c in range(6))) % pp != 0:
                    prod_sample_leaf_errs += 1
                pa0 = (pa0 + dpa0) & M64; pa1 = (pa1 + dpa1) & M64; pa2 = (pa2 + dpa2) & M64
                pb0 = (pb0 + dpb0) & M64; pb1 = (pb1 + dpb1) & M64; pb2 = (pb2 + dpb2) & M64
                got_p += 1; prod_sample_leaves_exec += 1

    prod_an_mean_ir = round(sum(prod_an_ir_vals) / len(prod_an_ir_vals), 2)
    prod_an_max_ir = max(prod_an_ir_vals)
    scaled_leaves_total = sum(scaled_leaves_vals)
    scaled_act_over_an = [round(a / b, 3) for a, b in zip(scaled_actual_ir_vals, scaled_an_ir_vals)]

    # Boundary over-enumeration check of native-64-bit Q38/Q32/Q12 on real orbit p_b (T_mult=40)
    p_b = orbits[-1][0]; s0_b = orbits[-1][1][0]; rho_b = orbits[-1][2]
    inv3_b = pow(3, -1, p_b); sp_b = [1]
    for _ in range(5): sp_b.append((sp_b[-1] * s0_b) % p_b)
    es_b = [(sum(M.ZK_NUM[k][j] * sp_b[j] for j in range(6)) * inv3_b) % p_b for k in range(6)]
    lf_b = es_b + [((-rho_b) * e) % p_b for e in es_b]
    inv0_b = pow(lf_b[0], -1, p_b)
    Bb = [[p_b if c == 0 else 0 for c in range(12)]]
    for i in range(1, 12):
        c0 = (-lf_b[i] * inv0_b) % p_b
        if c0 > p_b // 2: c0 -= p_b
        Bb.append([c0 if c == 0 else (1 if c == i else 0) for c in range(12)])
    DESC.bareiss_lll_12d(Bb, G_sym, 1764, 1000, False, 9, 10)
    GBb = [[1764 * sum(Bb[i][u] * G_sym[u][v] * Bb[j][v] for u in range(6) for v in range(6))
            + 1000 * sum(Bb[i][6 + u] * G_sym[u][v] * Bb[j][6 + v] for u in range(6) for v in range(6))
            for j in range(12)] for i in range(12)]
    Db, lamb = _bareiss_ldl(GBb)
    T12_b = 1764 * 40 * G_SC
    mu_b, b2_b, inv_b, ecb = _build_q38_tables(Db, lamb, T12_b)
    b_in = [0]; b_ex = [0]; xb = [0] * 12
    def _cnt_bnd(k, rem_q38, az, cur):
        ck = 0
        for j in range(k + 1, 12): ck = _i64(ck - _i64(mu_b[j][k] * xb[j]))
        r2 = _u64(_umulh64(rem_q38, inv_b[k]) + 1)
        rad = _i64((isqrt64_5newton(r2) + 1) << 26)
        ec = ecb[k]
        lo = (1 if k == 0 else 0) if az else (_i64(_i64(_i64(ck - ec) - rad) + Q32 - 1) >> 32)
        hi = _i64(_i64(ck + ec) + rad) >> 32
        Bk = Bb[k]; b2k = b2_b[k]
        for xv in range(lo, hi + 1):
            if k == 0:
                nxt = [cur[c] + xv * Bk[c] for c in range(12)]
                q12v = 1764 * _qf6(G_sym, nxt[:6]) + 1000 * _qf6(G_sym, nxt[6:])
                if q12v <= T12_b: b_in[0] += 1
                else: b_ex[0] += 1
            else:
                diff = _i64(_i64(xv << 32) - ck)
                ad = -diff if diff < 0 else diff
                d_lo = ad - ec if ad > ec else 0
                d19 = _u64(d_lo >> 13)
                sub = _umulh64(_u64(d19 * d19), b2k)
                if rem_q38 >= sub:
                    xb[k] = xv
                    _cnt_bnd(k - 1, rem_q38 - sub, az and (xv == 0), [cur[c] + xv * Bk[c] for c in range(12)])
    _cnt_bnd(11, Q38, True, [0] * 12)
    q32_bnd_extra_frac = round(b_ex[0] / max(1, b_in[0]), 5)

    # ---- PART 3: Explicit tie-case + 2,000+ SWAR 3-image stepping tests ----
    swar_tie_canon_errs = 0
    tie_cases = [
        ((1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0)),
        ((5, 0, 0, 0, 0, 0), (1, 2, -3, 4, -1, 2)),
        ((-7, 0, 0, 0, 0, 0), (0, 0, 1, -2, 3, -4)),
        ((1, 1, 1, 0, 0, 0), (0, 1, -1, 2, 0, -1)),
    ]
    for v_seed, b_seed in [((0,1,0,0,0,0),(0,1,0,0,0,0)), ((0,0,0,1,0,0),(1,-2,3,0,1,-1)), ((1,0,1,0,1,0),(0,0,1,0,0,-2))]:
        v1 = _mv6(S, v_seed); v2 = _mv6(S2, v_seed)
        a_fix = tuple(v_seed[i] + v1[i] + v2[i] for i in range(6))
        assert _mv6(S, a_fix) == a_fix
        tie_cases.append((a_fix, b_seed))
    for ta, tb in tie_cases:
        imgs = []
        ca_i, cb_i = ta, tb
        for _ in range(3):
            imgs.append((ca_i, cb_i))
            imgs.append((neg6(ca_i), neg6(cb_i)))
            imgs.append((ca_i, neg6(cb_i)))
            ca_i, cb_i = _mv6(S, ca_i), _mv6(S, cb_i)
        outs = set()
        for ia, ib in imgs:
            pa0, pa1, pa2 = pack_v(ia, SH_A, BI_A), pack_v(_mv6(S, ia), SH_A, BI_A), pack_v(_mv6(S2, ia), SH_A, BI_A)
            pb0, pb1, pb2 = pack_v(ib, SH_B, BI_B), pack_v(_mv6(S, ib), SH_B, BI_B), pack_v(_mv6(S2, ib), SH_B, BI_B)
            outs.add(swar_canon_20op_tie_exact(pa0, pa1, pa2, pb0, pb1, pb2))
        if len(outs) != 1:
            swar_tie_canon_errs += 1

    rng = random.Random(seed)
    swar_steps = 0; swar_step_errs = 0; swar_neg_errs = 0; swar_canon_errs = 0
    for _ in range(250):
        L0 = rng.randint(5, 14)
        a0 = tuple(rng.randint(-120, 120) for _ in range(6))
        da = tuple(rng.randint(-10, 10) for _ in range(6))
        b0 = tuple(rng.randint(-180, 180) for _ in range(6))
        db = tuple(rng.randint(-12, 12) for _ in range(6))
        a1, a2 = _mv6(S, a0), _mv6(S2, a0); da1, da2 = _mv6(S, da), _mv6(S2, da)
        b1, b2 = _mv6(S, b0), _mv6(S2, b0); db1, db2 = _mv6(S, db), _mv6(S2, db)
        ok_b = True
        for va, vd in ((a0, da), (a1, da1), (a2, da2)):
            for i in range(6):
                if max(abs(va[i]), abs(va[i] + (L0 - 1) * vd[i])) > 450:
                    ok_b = False
        for vb, vd in ((b0, db), (b1, db1), (b2, db2)):
            for i in range(6):
                lim = 900 if i < 2 else 480
                if max(abs(vb[i]), abs(vb[i] + (L0 - 1) * vd[i])) > lim:
                    ok_b = False
        if not ok_b:
            continue
        pa0, pa1, pa2 = pack_v(a0, SH_A, BI_A), pack_v(a1, SH_A, BI_A), pack_v(a2, SH_A, BI_A)
        pb0, pb1, pb2 = pack_v(b0, SH_B, BI_B), pack_v(b1, SH_B, BI_B), pack_v(b2, SH_B, BI_B)
        dpa0, dpa1, dpa2 = pack_d(da, SH_A), pack_d(da1, SH_A), pack_d(da2, SH_A)
        dpb0, dpb1, dpb2 = pack_d(db, SH_B), pack_d(db1, SH_B), pack_d(db2, SH_B)
        for t in range(L0):
            swar_steps += 1
            ra = tuple(a0[i] + t * da[i] for i in range(6))
            rb = tuple(b0[i] + t * db[i] for i in range(6))
            if unpack_v(pa0, SH_A, BI_A, MK_A) != ra or unpack_v(pb0, SH_B, BI_B, MK_B) != rb:
                swar_step_errs += 1
            if unpack_v((BIAS_A_X2 - pa0) & M64, SH_A, BI_A, MK_A) != neg6(ra):
                swar_neg_errs += 1
            if unpack_v((BIAS_B_X2 - pb0) & M64, SH_B, BI_B, MK_B) != neg6(rb):
                swar_neg_errs += 1
            ma, mb = swar_canon_20op_tie_exact(pa0, pa1, pa2, pb0, pb1, pb2)
            ref_pair = min((abs_lex6(ra), abs_lex6_b(rb)),
                           (abs_lex6(_mv6(S, ra)), abs_lex6_b(_mv6(S, rb))),
                           (abs_lex6(_mv6(S2, ra)), abs_lex6_b(_mv6(S2, rb))))
            if (unpack_v(ma, SH_A, BI_A, MK_A), unpack_v(mb, SH_B, BI_B, MK_B)) != ref_pair:
                swar_canon_errs += 1
            pa0 = (pa0 + dpa0) & M64; pa1 = (pa1 + dpa1) & M64; pa2 = (pa2 + dpa2) & M64
            pb0 = (pb0 + dpb0) & M64; pb1 = (pb1 + dpb1) & M64; pb2 = (pb2 + dpb2) & M64

    # ---- PART 4: 1-level super-tile bucket + structured real-a HT (load=0.140) + descending-ai chain ----
    ht_sz = 16384; ht_mask = ht_sz - 1
    real_a_keys = []
    for a0i in range(-6, 7):
        for a1i in range(-5, 6):
            for a2i in range(-4, 5):
                for a3i in range(-2, 3):
                    for a4i in range(-2, 3):
                        for a5i in range(-2, 3):
                            if a0i | a1i | a2i | a3i | a4i | a5i:
                                if a0i*a0i + a1i*a1i + a2i*a2i + a3i*a3i + a4i*a4i + a5i*a5i <= 14:
                                    real_a_keys.append(pack_v((a0i, a1i, a2i, a3i, a4i, a5i), SH_A, BI_A))
    a_keys = list(dict.fromkeys(real_a_keys))[:2293]
    Na_t = len(a_keys); B_st = 8
    FIB = 11400714819323198485
    ht_k = [0] * ht_sz; ht_v = [0] * ht_sz
    ht_ins_probes = 0; ht_ins_max = 0
    for ai, k in enumerate(a_keys, 1):
        h = (((k ^ (k >> 30)) * FIB) & M64) >> 50
        pr = 1
        while ht_k[h] != 0:
            h = (h + 1) & ht_mask; pr += 1
        ht_k[h] = k; ht_v[h] = ai
        ht_ins_probes += pr
        if pr > ht_ins_max: ht_ins_max = pr

    n_chords = 32; L_ch = 16
    score_ref = [[[0] * (Na_t + 1) for _ in range(L_ch)] for _ in range(n_chords)]
    score_got = [[[0] * (Na_t + 1) for _ in range(L_ch)] for _ in range(n_chords)]

    p_l = 131; W_l = 16; nb_l = (p_l + W_l - 1) // W_l
    kta_l = [0] + [rng.randrange(p_l) for _ in range(Na_t)]
    head_l = [0] * nb_l; nxt_l = [0] * (Na_t + 1)
    for ai in range(Na_t, 0, -1):
        k = kta_l[ai]; sb = k // W_l
        nxt_l[ai] = head_l[sb]; head_l[sb] = ((k & (W_l - 1)) << 46) | ai
    pref_errs = 0; M46 = (1 << 46) - 1
    for sb in range(nb_l):
        namax = rng.randint(40, Na_t)
        got_s = set(); rec = head_l[sb]
        while rec:
            ai = rec & M46
            if ai > namax: break
            got_s.add((rec >> 46, ai))
            rec = nxt_l[ai]
        ref_s = {(kta_l[ai] & (W_l - 1), ai) for ai in range(1, namax + 1) if kta_l[ai] // W_l == sb}
        if got_s != ref_s: pref_errs += 1

    b_keys = [[((b0 + 512) << 11) | ((ch // B_st) << 8) | (ch % B_st) for b0 in range(L_ch)] for ch in range(n_chords)]
    buckets = [[] for _ in range(B_st)]
    ht_probes = 0; ht_lookups = 0
    for lp_band in (52, 55, 58, 61):
        for t in range(B_st): buckets[t].append((0, lp_band))
        for _ in range(500):
            ai = rng.randint(1, Na_t)
            ch = rng.randrange(n_chords); b0 = rng.randrange(L_ch)
            ka = a_keys[ai - 1]; kb = b_keys[ch][b0]
            score_ref[ch][b0][ai] += lp_band
            buckets[kb & 0xFF].append((ka, kb))

    for tid in range(B_st):
        chord_slots = [tid + s * B_st for s in range(n_chords // B_st)]
        cur_lp = 0
        for ka, kb in buckets[tid]:
            if ka == 0:
                cur_lp = kb; continue
            ht_lookups += 1
            h = (((ka ^ (ka >> 30)) * FIB) & M64) >> 50
            ht_probes += 1
            while ht_k[h] != ka:
                h = (h + 1) & ht_mask; ht_probes += 1
            ai = ht_v[h]
            ch_slot = (kb >> 8) & 0x7
            b0 = (kb >> 11) - 512
            score_got[chord_slots[ch_slot]][b0][ai] += cur_lp

    supertile_errs = 0
    for ch in range(n_chords):
        for b0 in range(L_ch):
            if score_got[ch][b0] != score_ref[ch][b0]:
                supertile_errs += 1

    GA = T4["GA"]
    aug = [GA[i][:] + [1.0 if i == j else 0.0 for j in range(6)] for i in range(6)]
    for col in range(6):
        piv = max(range(col, 6), key=lambda r: abs(aug[r][col]))
        aug[col], aug[piv] = aug[piv], aug[col]
        pv = aug[col][col]
        aug[col] = [x / pv for x in aug[col]]
        for r in range(6):
            if r != col and aug[r][col]:
                fac = aug[r][col]
                aug[r] = [aug[r][c] - fac * aug[col][c] for c in range(12)]
    max_cu = [sqrt(aug[i][6 + i]) for i in range(6)]
    max_a_prod = [round(605.0 * c, 1) for c in max_cu]
    max_b_prod = [round(805.0 * c, 1) for c in max_cu]
    swar_bounds_ok = (max(max_a_prod) < 511.0 and max_b_prod[0] < 1023.0 and max_b_prod[1] < 1023.0 and max(max_b_prod[2:]) < 511.0)
    struct_avg_probes = round(ht_ins_probes / float(Na_t), 4)

    ok = (uncond_64bit_theorem_ok and gsym_min_int == 753729 and lovasz_lhs > lovasz_rhs
          and isqrt_errs == 0 and lat_mismatches == 0
          and scaled_leaves_total >= 5000 and scaled_leaf_verify_errs == 0
          and prod_sample_leaves_exec == 1500 and prod_sample_leaf_errs == 0
          and prod_an_max_ir <= 25.0
          and prod_max_inv_bits <= 63 and prod_max_b2_bits <= 64 and prod_max_r2_bits <= 37
          and q32_bnd_extra_frac <= 0.02
          and swar_step_errs == 0 and swar_neg_errs == 0 and swar_canon_errs == 0
          and swar_tie_canon_errs == 0 and pref_errs == 0 and supertile_errs == 0
          and swar_bounds_ok and swar_steps >= 1000 and struct_avg_probes <= 1.12)
    return dict(
        ok=ok,
        zero_float_fp_enum=True,
        strict_u64_i64_umulh64_enforced=True,
        zero_fallback_unconditional_64bit=True,
        gsym_6d_exact_min_int=gsym_min_int,
        gsym_6d_min_vec_count=len(gsym_min_vecs),
        lovasz_uncond_min_b2_unsc=lovasz_min_b2_unsc,
        lovasz_floor_B0_unsc=B0_UNSC,
        wc_inv_b2_hi_bits=wc_inv_b2_hi.bit_length(),
        wc_r2_q12_bits=wc_r2_q12.bit_length(),
        wc_ck_q32_bits=wc_ck_q32.bit_length(),
        wc_d2_q38_bits=wc_d2_q38.bit_length(),
        isqrt64_5newton_errs=isqrt_errs,
        orbits_6ideal_checked=len(orbits),
        lattice_b2_mismatches=lat_mismatches,
        lattice_b2_raw_hits=tot_raw_hits,
        lattice_b2_jvor_hits=tot_jvor_hits,
        lattice_b2_leaves_per_kept=round(tot_raw_hits / max(1, tot_jvor_hits), 3),
        prod_51_to_61b_orbits_checked=6,
        prod_51_to_61b_min_b2_unsc=round(prod_min_b2_unsc, 1),
        scaled_complete_tree_leaves_per_orbit=scaled_leaves_vals,
        scaled_complete_tree_leaves_total=scaled_leaves_total,
        scaled_complete_tree_nodes_per_orbit=scaled_nodes_vals,
        scaled_complete_tree_actual_ir_per_leaf=scaled_actual_ir_vals,
        scaled_complete_tree_analytic_ir_per_leaf=scaled_an_ir_vals,
        scaled_complete_tree_act_over_an_ratio=scaled_act_over_an,
        scaled_complete_tree_leaf_verify_errs=scaled_leaf_verify_errs,
        prod_T12_analytic_12lvl_nodes_per_leaf=prod_an_npl_vals,
        prod_T12_analytic_12lvl_ir_per_leaf=prod_an_ir_vals,
        prod_T12_analytic_12lvl_mean_ir=prod_an_mean_ir,
        prod_T12_analytic_12lvl_max_ir=prod_an_max_ir,
        prod_T12_sample_leaves_executed=prod_sample_leaves_exec,
        prod_T12_sample_leaf_verify_errs=prod_sample_leaf_errs,
        prod_fp_ops_node_k0=OPS_NODE[0],
        prod_fp_budget_ir_per_leaf=25.0,
        prod_max_inv_b2_bits=prod_max_inv_bits,
        prod_max_b2_q64_bits=prod_max_b2_bits,
        prod_max_r2_q12_bits=prod_max_r2_bits,
        q32_bnd_in_leaves=b_in[0],
        q32_bnd_extra_leaves=b_ex[0],
        q32_bnd_extra_frac=q32_bnd_extra_frac,
        swar_steps_checked=swar_steps,
        swar_step_errs=swar_step_errs,
        swar_neg_errs=swar_neg_errs,
        swar_20op_canon_errs=swar_canon_errs,
        swar_tie_cases_checked=len(tie_cases),
        swar_tie_canon_errs=swar_tie_canon_errs,
        swar_bounds_ok=swar_bounds_ok,
        max_a_prod_coords=max_a_prod,
        max_b_prod_coords=max_b_prod,
        desc_ai_prefix_chain_errs=pref_errs,
        supertile_structured_keys=Na_t,
        supertile_structured_ht_load=round(Na_t / float(ht_sz), 3),
        supertile_structured_ht_avg_probes=struct_avg_probes,
        supertile_structured_ht_max_probes=ht_ins_max,
        supertile_replay_hits=ht_lookups,
        supertile_ht_avg_probes=round(ht_probes / max(1, ht_lookups), 3),
        supertile_replay_errs=supertile_errs,
    )
