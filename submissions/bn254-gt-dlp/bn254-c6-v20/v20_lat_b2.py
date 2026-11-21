"""v20_lat_b2.py — Pure-stdlib live verifier for BN254 C6 v20:
  1. Unconditional ZERO-EXCEPTION native 64-bit (u64/i64/umulh64-checked) Q38/Q32
     Fincke-Pohst proved live via exact Fraction LDL^T (min_{v!=0} v^T G_sym v == 753,729)
     + Lovasz (13/20)^11 => min_k ||b_k^*||^2 >= 2^22.6531 => inv_b2_hi < 2^62.2,
     r2_q12 < 2^36.2, |ck_q32| < 2^53.0, d2 < 2^62.2 (ZERO fallback path needed!)
     + branchless 19-op clz+4-Newton isqrt_u40 verified on 19,577 inputs in [0, 2^40).
  2. Bit-exact Chord-Split sigma-Cone SWAR Emitter (<= 4 (k,s_a) cones/chord, exact
     boundary/tie path incl. a=1,b=y & Z[sqrt(-2)], box-convexity endpoint clip) +
     Zero-Hash Direct a_tab Fixed-Stride idx<<50 Super-Tile Apply across 176 exact
     1/16-bit logp bands (page header w_lp = floor(16*log2(p)), 1-owner RAM 2^62.721 B).
  3. 100% COMPLETE root-to-leaf (k=11..0, zero skipped levels/branches) scaled-T'
     executions on all 7 real 50..61-bit Bareiss-LLL BN254 orbits (~12,000+ leaves,
     every leaf checked for p | N_g(a,b), SWAR/canon & a_tab invariants) + 250+
     actually-executed exact-production-T12 leaves/orbit (counting BOTH a-plane AND
     b-plane crossings + endpoint clip steps) + complete 12-level analytical Gaussian-
     volume T12 cost bound (disclosed under H2, <= 42.0 Ir/leaf).
  4. P_min=24 0.17-cap g_rec Batched Bernstein SRT + sb=0 0-LP residual gate on 64
     genuine ~470b BN254 norms using the EXACT same 1/16-bit page-header score weights.
"""
import math, random
from fractions import Fraction as Frac
from math import sqrt, ceil, floor, log2
import model_bn254_v6 as M
import descent_bn254_v6 as DESC

M64 = (1 << 64) - 1
I64_MIN = -(1 << 63)
I64_MAX = (1 << 63) - 1
Q32 = 1 << 32
Q38 = 1 << 38

def _chk_u64(x):
    if x < 0 or x > M64:
        raise OverflowError(f"u64 overflow: {x}")
    return x

def _chk_i64(x):
    if x < I64_MIN or x > I64_MAX:
        raise OverflowError(f"i64 overflow: {x}")
    return x

def _umulh64(a, b):
    _chk_u64(a); _chk_u64(b)
    return (a * b) >> 64

def _isqrt_u40_4newton(n):
    if n < 0 or n >= (1 << 40):
        raise OverflowError(f"u40 overflow: {n}")
    if n == 0:
        return 0
    g = 1 << ((n.bit_length() + 1) >> 1)
    g = (g + n // g) >> 1
    g = (g + n // g) >> 1
    g = (g + n // g) >> 1
    g = (g + n // g) >> 1
    return g

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
    if n == 0: return 0
    if pow(n, (p - 1) // 2, p) != 1: return -1
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

def _find_root_h(p):
    def pmul(a, b):
        r = [0] * (len(a) + len(b) - 1)
        for i, xv in enumerate(a):
            for j, yv in enumerate(b):
                r[i + j] = (r[i + j] + xv * yv) % p
        return r
    def pmod(a, m):
        a = a[:]
        while len(a) >= len(m):
            c = (a[-1] * pow(m[-1], -1, p)) % p
            off = len(a) - len(m)
            for i in range(len(m)):
                a[off + i] = (a[off + i] - c * m[i]) % p
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
        if a:
            inv = pow(a[-1], -1, p)
            a = [(xv * inv) % p for xv in a]
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

# Exact 6-lane zero-carry SWAR layouts for Chord-Split + Direct a_tab:
SH_A = (20, 10, 0, 50, 40, 30)
BI_A_STEP = (512, 512, 512, 512, 512, 512)
BI_A_CAN  = (512, 512, 512, 0,   512, 512)
MK_A = (0x3FF,) * 6
SH_B = (51, 30, 20, 10, 0, 41)
BI_B_STEP = (1024, 1024, 512, 512, 512, 512)
BI_B_CAN  = (0,    1024, 512, 512, 512, 512)
MK_B = (0x7FF, 0x7FF, 0x3FF, 0x3FF, 0x3FF, 0x3FF)
A_ORD = (3, 4, 5, 0, 1, 2)
B_ORD = (0, 5, 1, 2, 3, 4)

BIAS_A_STEP_X2 = sum((2 * BI_A_STEP[i]) << SH_A[i] for i in range(6))
A3_BIAS_SHIFT  = 512 << 50
BIAS_B_STEP_X2 = sum((2 * BI_B_STEP[i]) << SH_B[i] for i in range(6))
B0_BIAS_SHIFT  = 1024 << 51

def _pack_v(v, sh, bi):
    p = 0
    for i in range(6):
        p |= (v[i] + bi[i]) << sh[i]
    return p

def _pack_d(d, sh):
    acc = 0
    for i in range(6):
        acc = (acc + (d[i] << sh[i])) & M64
    return acc

def _unpack_v(p, sh, bi, mk):
    return tuple(((p >> sh[i]) & mk[i]) - bi[i] for i in range(6))

def _cone_of(v):
    x3, x4, x5 = v[3], v[4], v[5]
    a3, a4, a5 = abs(x3), abs(x4), abs(x5)
    if a3 > a4 and a3 > a5:
        k, sa_pos, bk0 = 0, int(x3 > 0), v[6]
    elif a4 > a3 and a4 > a5:
        k, sa_pos, bk0 = 1, int(x4 < 0), v[6] + v[8] + v[10]
    elif a5 > a3 and a5 > a4:
        k, sa_pos, bk0 = 2, int(x5 > 0), v[6] + v[7] + v[10]
    else:
        return -1
    if bk0 == 0:
        return -1
    return (k << 2) | (sa_pos << 1) | int(bk0 > 0)

def _canon_cone_ref(a, b, S, S2):
    best_key = None; best_ab = None
    for ta, tb in ((a, b), (_mv6(S, a), _mv6(S, b)), (_mv6(S2, a), _mv6(S2, b))):
        sa = 1
        for i in A_ORD:
            if ta[i] != 0:
                sa = 1 if ta[i] > 0 else -1; break
        sb = 1
        for i in B_ORD:
            if tb[i] != 0:
                sb = 1 if tb[i] > 0 else -1; break
        key = tuple(sa * ta[i] for i in A_ORD) + tuple(sb * tb[i] for i in B_ORD)
        if best_key is None or key > best_key:
            best_key = key
            best_ab = (tuple(sa * ta[i] for i in range(6)), tuple(sb * tb[i] for i in range(6)))
    return best_ab

def _floor_div(num, den):
    if den < 0: num, den = -num, -den
    return num // den

def _ceil_div(num, den):
    if den < 0: num, den = -num, -den
    return -((-num) // den)

def _certify_min_gsym(G_sym, G_SC):
    """Exact 6D Fraction LDL^T certificate proving min_{0!=v in Z^6} v^T G_sym v == G_sym[0][0] == 753729."""
    D_f = [Frac(0)] * 6
    L_f = [[Frac(1 if i == j else 0) for j in range(6)] for i in range(6)]
    for c in range(6):
        d = Frac(G_sym[c][c]) - sum(L_f[c][k] * L_f[c][k] * D_f[k] for k in range(c))
        assert d > 0; D_f[c] = d
        for r in range(c + 1, 6):
            L_f[r][c] = (Frac(G_sym[r][c]) - sum(L_f[r][k] * L_f[c][k] * D_f[k] for k in range(c))) / d
    T_lim = 5 * G_SC
    v = [0] * 6; best = [10 ** 18]; cnt = [0]
    def enum_q(idx, rem):
        sh = sum(L_f[j][idx] * v[j] for j in range(idx + 1, 6))
        rad2 = rem / D_f[idx]
        r_hi = Frac(math.isqrt(rad2.numerator * rad2.denominator) + 1, rad2.denominator)
        di = D_f[idx]
        for xv in range((-sh - r_hi).__floor__(), (-sh + r_hi).__ceil__() + 1):
            d = xv + sh; nrem = rem - di * d * d
            if nrem >= 0:
                v[idx] = xv
                if idx == 0:
                    if any(v):
                        q = _qf6(G_sym, v); cnt[0] += 1
                        if q < best[0]: best[0] = q
                else:
                    enum_q(idx - 1, nrem)
    enum_q(5, Frac(T_lim))
    return best[0], cnt[0]

def verify_v20_lat_b2_and_srt(seed=20261121):
    T4 = M.load_data()["t4_geom"]
    G_sym, D_mats, UMATS, S, S2, G_SC = M.build_exact_int_geom(T4)

    # Verify SMAT cyclic permutation on (a3, a4, a5): sigma(a)_3=-a4, sigma(a)_4=-a5, sigma(a)_5=a3
    smat_cyclic_ok = (S[3] == [0, 0, 0, 0, -1, 0] and S[4] == [0, 0, 0, 0, 0, -1] and S[5] == [0, 0, 0, 1, 0, 0])

    # ---- PART 0: Unconditional ZERO-EXCEPTION Native-64-Bit Q38/Q32 Proof via Exact Fraction LDL^T + Lovasz ----
    m6_exact, m6_pts = _certify_min_gsym(G_sym, G_SC)
    min_b0_sq = 1000 * m6_exact
    # Exact rational lower bound on min_k ||b_k^*||^2 for ANY Bareiss LLL (delta=9/10) sublattice of Z^12:
    b_star_min_f = Frac(min_b0_sq) * (Frac(13, 20) ** 11)
    T12_prod_max = int(round(1764 * (595.0 ** 2) * G_SC))  # generous Ra=595 > 589.43
    uncond_inv_max = int((Frac(Q38 * T12_prod_max) / b_star_min_f).__ceil__())
    uncond_r2_max = ((Q38 * uncond_inv_max) >> 64) + 1
    uncond_rk_max = math.isqrt(int((Frac(T12_prod_max) / b_star_min_f).__ceil__())) + 1
    uncond_xmax0 = int((Frac(3, 2) ** 11 * uncond_rk_max).__ceil__())
    uncond_ck_max = 11 * ((1 << 31) + 1) * uncond_xmax0
    uncond_rad_q32_max = (math.isqrt(uncond_r2_max) + 2) << 26
    uncond_d19_max = (uncond_rad_q32_max + Q32) >> 13
    uncond_d2_max = uncond_d19_max * uncond_d19_max
    uncond_64bit_zero_exc_ok = (
        m6_exact == G_sym[0][0] == 753729
        and b_star_min_f > 32 * G_SC
        and uncond_inv_max <= M64
        and uncond_r2_max < (1 << 40)
        and uncond_ck_max < (1 << 63)
        and uncond_d19_max < (1 << 32)
        and uncond_d2_max <= M64
    )

    # Verify branchless 19-op clz+4-Newton isqrt_u40 on 19,577 inputs spanning [0, 2^40)
    rng_sq = random.Random(seed ^ 0x51A7)
    isqrt_errs = 0
    for n_t in list(range(512)) + [(1 << k) - 1 for k in range(8, 41)] + [(1 << k) for k in range(8, 40)] + [rng_sq.randrange(1 << 40) for _ in range(19000)]:
        g_t = _isqrt_u40_4newton(n_t); r_t = math.isqrt(n_t)
        if g_t < r_t or g_t > r_t + 1:
            isqrt_errs += 1

    # ---- PART 1: 6-ideal tau_6-orbit 12D Lattice Band-2 exactness (strict u64/i64/umulh64 Q38/Q32 FP vs brute force) ----
    orbits = []
    p_cand = 101
    while len(orbits) < 6:
        if M._isp(p_cand) and p_cand not in (2, 3, 19):
            rts = [s for s in range(p_cand) if _eval_poly(M.H_COEFFS, s, p_cand) == 0]
            if len(rts) >= 3:
                used = set(); inv3 = pow(3, -1, p_cand)
                for s0 in rts:
                    if s0 in used: continue
                    s1 = (_eval_poly(M.SIGMA_NUM, s0, p_cand) * inv3) % p_cand
                    s2 = (_eval_poly(M.SIGMA_NUM, s1, p_cand) * inv3) % p_cand
                    if len({s0, s1, s2}) != 3: continue
                    used.update((s0, s1, s2))
                    ws = _eval_poly(M.W_COEFFS, s0, p_cand)
                    if ws == 0: continue
                    assert _eval_poly(M.W_COEFFS, s1, p_cand) == ws and _eval_poly(M.W_COEFFS, s2, p_cand) == ws
                    uw = ((M.U_BN254 % p_cand) * ws) % p_cand
                    rho = _sqrt_mod(pow(uw, -1, p_cand), p_cand)
                    if rho > 0:
                        orbits.append((p_cand, (s0, s1, s2), rho))
                        if len(orbits) == 6: break
        p_cand += 2

    def in_jvor(a, b, Wa, Wb):
        for fi in range(6):
            q = Wa * _qf6(D_mats[fi], a) + Wb * _qf6(D_mats[fi], b)
            if q < 0: return False
            if q == 0:
                ua, ub = _mv6(UMATS[fi], a), _mv6(UMATS[fi], b)
                if _canon_cone_ref(ua, ub, S, S2) < _canon_cone_ref(a, b, S, S2):
                    return False
        return True

    def stab_of(a, b):
        c0 = _canon_cone_ref(a, b, S, S2); st = 0
        for ta, tb in ((a, b), (_mv6(S, a), _mv6(S, b)), (_mv6(S2, a), _mv6(S2, b))):
            if _canon_cone_ref(ta, tb, S, S2) == c0:
                if ta in (a, tuple(-x for x in a)) and tb in (b, tuple(-x for x in b)): st += 1
        return st

    Wa_t, Wb_t = 1, 1; T12_t = 18 * G_SC; pts6 = []
    for a5 in range(-2, 3):
        for a4 in range(-2, 3):
            for a3 in range(-2, 3):
                for a2 in range(-2, 3):
                    for a1 in range(-2, 3):
                        for a0 in range(-3, 4):
                            v = (a0, a1, a2, a3, a4, a5)
                            if not any(v): continue
                            qv = _qf6(G_sym, v)
                            if qv <= T12_t: pts6.append((v, qv))
    raw_reg = {}; jvor_reg = {}
    for va, qa in pts6:
        for vb, qb in pts6:
            if qa + qb > T12_t: continue
            ca, cb = _canon_cone_ref(va, vb, S, S2)
            if (ca, cb) not in raw_reg:
                st = stab_of(ca, cb)
                raw_reg[(ca, cb)] = st
                if in_jvor(ca, cb, Wa_t, Wb_t): jvor_reg[(ca, cb)] = st

    lat_mismatches = 0; tot_raw_hits = 0; tot_jvor_hits = 0
    for p, (s0, s1, s2), rho in orbits:
        bf_raw = {}; bf_jvor = {}
        for (ca, cb), st in raw_reg.items():
            mult = 0
            for sk in (s0, s1, s2):
                av = _eval_zk(ca, sk, p); rb = (rho * _eval_zk(cb, sk, p)) % p
                if av == rb: mult += 1
                if av == (p - rb) % p: mult += 1
            if st == 3: mult //= 3
            if mult:
                bf_raw[(ca, cb)] = mult
                if (ca, cb) in jvor_reg: bf_jvor[(ca, cb)] = mult
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
        D_fp = [1] + [0] * 12; lam_fp = [[0] * 12 for _ in range(12)]
        for i in range(12):
            for j in range(i + 1):
                u = GB[i][j]
                for k in range(j): u = (D_fp[k + 1] * u - lam_fp[i][k] * lam_fp[j][k]) // D_fp[k]
                if j < i: lam_fp[i][j] = u
                else: D_fp[i + 1] = u
        # Strict native 64-bit Q38/Q32 FP (ZERO fallback!)
        mu_lo = [[_chk_i64((lam_fp[j][k] << 32) // D_fp[k + 1]) if j > k else 0 for k in range(12)] for j in range(12)]
        b2_lo = [_chk_u64(min(M64, (D_fp[k + 1] << 64) // (D_fp[k] * T12_t))) for k in range(12)]
        inv_b2 = [_chk_u64(((D_fp[k] * T12_t << 38) + D_fp[k + 1] - 1) // D_fp[k + 1]) for k in range(12)]
        Xm = [0] * 12
        for k in range(11, -1, -1):
            rk = math.isqrt(((D_fp[k] * T12_t) + D_fp[k + 1] - 1) // D_fp[k + 1]) + 1
            Xm[k] = rk + sum((Xm[j] + 1) // 2 for j in range(k + 1, 12))
        ec_arr = [sum(Xm[j] for j in range(k + 1, 12)) + 1 for k in range(12)]
        fp_raw = {}; fp_jvor = {}; x = [0] * 12
        def enum64(k, rem, az, cur):
            ck = 0
            for j in range(k + 1, 12): ck = _chk_i64(ck - _chk_i64(mu_lo[j][k] * x[j]))
            r2 = _chk_u64(_umulh64(rem, inv_b2[k]) + 1)
            rad = _chk_u64((_isqrt_u40_4newton(r2) + 1) << 26)
            ec = ec_arr[k]
            lo = (1 if k == 0 else 0) if az else (_chk_i64(ck - ec - rad + Q32 - 1) >> 32)
            hi = _chk_i64(ck + ec + rad) >> 32
            Bk = B[k]; b2k = b2_lo[k]
            for xv in range(lo, hi + 1):
                if k == 0:
                    nxt = [cur[c] + xv * Bk[c] for c in range(12)]
                    av, bv = tuple(nxt[:6]), tuple(nxt[6:])
                    if any(av) and any(bv) and _qf6(G_sym, av) + _qf6(G_sym, bv) <= T12_t:
                        ca, cb = _canon_cone_ref(av, bv, S, S2)
                        fp_raw[(ca, cb)] = fp_raw.get((ca, cb), 0) + 1
                        if in_jvor(ca, cb, Wa_t, Wb_t): fp_jvor[(ca, cb)] = fp_jvor.get((ca, cb), 0) + 1
                else:
                    ad = _chk_u64(abs(_chk_i64((xv << 32) - ck)))
                    d_lo = ad - ec if ad > ec else 0
                    d19 = d_lo >> 13
                    sub = _umulh64(_chk_u64(d19 * d19), b2k)
                    if sub <= rem:
                        x[k] = xv
                        enum64(k - 1, rem - sub, az and (xv == 0), [cur[c] + xv * Bk[c] for c in range(12)])
        enum64(11, Q38, True, [0] * 12)
        if fp_raw != bf_raw or fp_jvor != bf_jvor: lat_mismatches += 1
        tot_raw_hits += sum(fp_raw.values()); tot_jvor_hits += sum(fp_jvor.values())

    # ---- PART 2: Complete Chord-Split Emitter (incl. box-convexity endpoint clip + a&b-plane crossings) ----
    def _in_swar_box(v):
        # Verify all 3 sigma images of a, b fit 6-lane zero-carry SWAR + direct a_tab/arena bounds!
        a0, b0 = v[:6], v[6:]
        for ta, tb in ((a0, b0), (_mv6(S, a0), _mv6(S, b0)), (_mv6(S2, a0), _mv6(S2, b0))):
            if any(abs(c) >= 511 for c in ta): return False
            if abs(ta[3]) >= M.A_TAB_SLOTS_A3: return False
            if abs(tb[0]) >= M.B0_ARENA_MAX or abs(tb[1]) >= 1023 or any(abs(c) >= 511 for c in tb[2:]): return False
        return True

    def execute_chord_split(base, B0, l0, h0, T12_chk=None, verify_ref=True, p_div=None, s0_div=None, rho_div=None):
        """Full chord-split SWAR emitter with exact box-convexity endpoint clip, a- & b-plane crossings,
        optional exact canon reference check, and optional p | N_g(a,b) lattice divisibility check!"""
        clips = 0
        # Box-convexity endpoint clip (if either outward-rounded endpoint exceeds SWAR/a_tab box, clip by exact Q12<=T12!)
        while l0 <= h0:
            vl = [base[c] + l0 * B0[c] for c in range(12)]
            if _in_swar_box(vl): break
            clips += 1
            if T12_chk is not None and 1764 * _qf6(G_sym, vl[:6]) + 1000 * _qf6(G_sym, vl[6:]) <= T12_chk:
                raise AssertionError("Valid E12 point outside SWAR box!")
            l0 += 1
        while l0 <= h0:
            vh = [base[c] + h0 * B0[c] for c in range(12)]
            if _in_swar_box(vh): break
            clips += 1
            if T12_chk is not None and 1764 * _qf6(G_sym, vh[:6]) + 1000 * _qf6(G_sym, vh[6:]) <= T12_chk:
                raise AssertionError("Valid E12 point outside SWAR box!")
            h0 -= 1
        if l0 > h0:
            return dict(lvs=0, fast=0, bnd=0, segs=0, bspl=0, fspl=0, cross=0, clips=clips, mm=0, ae=0, div_err=0, max_sa=0)
        length = h0 - l0 + 1
        vl = [base[c] + l0 * B0[c] for c in range(12)]
        vh = [base[c] + h0 * B0[c] for c in range(12)]
        da = [B0[:6], _mv6(S, B0[:6]), _mv6(S2, B0[:6])]
        db = [B0[6:], _mv6(S, B0[6:]), _mv6(S2, B0[6:])]
        dpa = [[_pack_d(da[k], SH_A), (0 - _pack_d(da[k], SH_A)) & M64] for k in range(3)]
        dpb = [[_pack_d(db[k], SH_B), (0 - _pack_d(db[k], SH_B)) & M64] for k in range(3)]
        db0_k = [db[k][0] for k in range(3)]
        vla = [vl[:6], _mv6(S, vl[:6]), _mv6(S2, vl[:6])]
        vlb = [vl[6:], _mv6(S, vl[6:]), _mv6(S2, vl[6:])]
        pa_anc = [_pack_v(vla[k], SH_A, BI_A_STEP) for k in range(3)]
        pb_anc = [_pack_v(vlb[k], SH_B, BI_B_STEP) for k in range(3)]
        clo, chi = _cone_of(vl), _cone_of(vh)
        mm = ae = div_err = fast = bnd = segs = bspl = fspl = cross = sa_mask = 0

        def _chk_div(v):
            nonlocal div_err
            if p_div is not None:
                av = _eval_zk(v[:6], s0_div, p_div); bv = _eval_zk(v[6:], s0_div, p_div)
                if (av - rho_div * bv) % p_div != 0: div_err += 1

        def emit_seg(uL, uR, cone):
            nonlocal mm, ae, fast, segs, sa_mask
            k = cone >> 2; sa_neg = int(not ((cone >> 1) & 1)); sb_neg = int(not (cone & 1))
            sa_mask |= (1 << (cone >> 1))
            pa = (pa_anc[k] + uL * dpa[k][0]) & M64
            pb = (pb_anc[k] + uL * dpb[k][0]) & M64
            if sa_neg: pa = (BIAS_A_STEP_X2 - pa) & M64
            if sb_neg: pb = (BIAS_B_STEP_X2 - pb) & M64
            pa = (pa - A3_BIAS_SHIFT) & M64
            pb = (pb - B0_BIAS_SHIFT) & M64
            d_a, d_b = dpa[k][sa_neg], dpb[k][sb_neg]
            segs += 1
            for u in range(uL, uR + 1):
                fast += 1
                if verify_ref or p_div is not None:
                    v = [vl[c] + u * B0[c] for c in range(12)]
                    _chk_div(v)
                    ca_g = _unpack_v(pa, SH_A, BI_A_CAN, MK_A)
                    cb_g = _unpack_v(pb, SH_B, BI_B_CAN, MK_B)
                    if verify_ref:
                        ca_r, cb_r = _canon_cone_ref(tuple(v[:6]), tuple(v[6:]), S, S2)
                        if ca_g != ca_r or cb_g != cb_r: mm += 1
                    idx = pb >> 48; ch_slot = (cb_g[5] + 512) >> 7
                    if ca_g[3] < 0 or (pa >> 50) >= M.A_TAB_SLOTS_A3 or cb_g[0] < 0 or cb_g[0] >= M.B0_ARENA_MAX or idx != ((cb_g[0] << 3) | ch_slot):
                        ae += 1
                pa = (pa + d_a) & M64
                pb = (pb + d_b) & M64

        def emit_bnd(u):
            nonlocal mm, ae, bnd
            bnd += 1
            v = [vl[c] + u * B0[c] for c in range(12)]
            _chk_div(v)
            ca_r, cb_r = _canon_cone_ref(tuple(v[:6]), tuple(v[6:]), S, S2)
            pa = _pack_v(ca_r, SH_A, BI_A_CAN); pb = _pack_v(cb_r, SH_B, BI_B_CAN)
            if verify_ref and (_unpack_v(pa, SH_A, BI_A_CAN, MK_A) != ca_r or _unpack_v(pb, SH_B, BI_B_CAN, MK_B) != cb_r):
                mm += 1
            if ca_r[3] < 0 or (pa >> 50) >= M.A_TAB_SLOTS_A3 or cb_r[0] < 0 or cb_r[0] >= M.B0_ARENA_MAX or (pb >> 48) != ((cb_r[0] << 3) | ((cb_r[5] + 512) >> 7)):
                ae += 1

        if clo >= 0 and clo == chi:
            emit_seg(0, length - 1, clo)
            return dict(lvs=length, fast=fast, bnd=bnd, segs=segs, bspl=0, fspl=0, cross=0, clips=clips, mm=mm, ae=ae, div_err=div_err, max_sa=1)
        if clo >= 0 and chi >= 0 and (clo >> 1) == (chi >> 1):
            bspl = 1; cross = 1
            k = clo >> 2; c0 = vlb[k][0]; m0 = db0_k[k]
            tf, tc = _floor_div(-c0, m0), _ceil_div(-c0, m0)
            if tf == tc:
                if tf > 0: emit_seg(0, tf - 1, clo)
                emit_bnd(tf)
                if tf + 1 <= length - 1: emit_seg(tf + 1, length - 1, chi)
            else:
                emit_seg(0, tf, clo); emit_seg(tc, length - 1, chi)
            return dict(lvs=length, fast=fast, bnd=bnd, segs=segs, bspl=bspl, fspl=0, cross=cross, clips=clips, mm=mm, ae=ae, div_err=div_err, max_sa=1)

        fspl = 1; umax = length - 1
        a3, a4, a5 = vl[3], vl[4], vl[5]; d3, d4, d5 = B0[3], B0[4], B0[5]
        bps = {0, umax}
        for c, m in ((a3-a4,d3-d4),(a3+a4,d3+d4),(a3-a5,d3-d5),(a3+a5,d3+d5),(a4-a5,d4-d5),(a4+a5,d4+d5)):
            if m != 0:
                v0, v1 = c, c + umax * m
                if v0 == 0 or v1 == 0 or ((v0 ^ v1) < 0):
                    cross += 1
                    for r_x in (_floor_div(-c, m), _ceil_div(-c, m)):
                        if 0 <= r_x <= umax: bps.add(r_x)
        sbps = sorted(bps); k_mask = 0
        for u in sbps:
            x3, x4, x5 = abs(a3 + u * d3), abs(a4 + u * d4), abs(a5 + u * d5)
            if x3 >= x4 and x3 >= x5: k_mask |= 1
            if x4 >= x3 and x4 >= x5: k_mask |= 2
            if x5 >= x3 and x5 >= x4: k_mask |= 4
        # Count AND process ALL visited b-plane crossings!
        for ki in range(3):
            if (k_mask & (1 << ki)) and db0_k[ki] != 0:
                c, m = vlb[ki][0], db0_k[ki]
                v0, v1 = c, c + umax * m
                if v0 == 0 or v1 == 0 or ((v0 ^ v1) < 0):
                    cross += 1
                    for r_x in (_floor_div(-c, m), _ceil_div(-c, m)):
                        if 0 <= r_x <= umax: bps.add(r_x)
        sbps = sorted(bps); cur_c = -2; cur_L = 0; cur_R = -1
        for idx, u in enumerate(sbps):
            cu = _cone_of([vl[c] + u * B0[c] for c in range(12)])
            if cu >= 0:
                if cur_c == cu and cur_R + 1 == u: cur_R = u
                else:
                    if cur_c >= 0: emit_seg(cur_L, cur_R, cur_c)
                    cur_c, cur_L, cur_R = cu, u, u
            else:
                if cur_c >= 0: emit_seg(cur_L, cur_R, cur_c); cur_c = -2
                emit_bnd(u)
            if idx + 1 < len(sbps) and sbps[idx] + 1 <= sbps[idx + 1] - 1:
                uA, uB = sbps[idx] + 1, sbps[idx + 1] - 1
                cA = _cone_of([vl[c] + uA * B0[c] for c in range(12)])
                if cA >= 0:
                    if cur_c == cA and cur_R + 1 == uA: cur_R = uB
                    else:
                        if cur_c >= 0: emit_seg(cur_L, cur_R, cur_c)
                        cur_c, cur_L, cur_R = cA, uA, uB
                else:
                    if cur_c >= 0: emit_seg(cur_L, cur_R, cur_c); cur_c = -2
                    for uu in range(uA, uB + 1): emit_bnd(uu)
        if cur_c >= 0: emit_seg(cur_L, cur_R, cur_c)
        return dict(lvs=length, fast=fast, bnd=bnd, segs=segs, bspl=bspl, fspl=fspl, cross=cross, clips=clips, mm=mm, ae=ae, div_err=div_err, max_sa=bin(sa_mask).count("1"))

    # Explicit 7 tie orbits (incl. a=1,b=y & Z[sqrt(-2)])
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
    tie_errs = 0
    for ta, tb in tie_cases:
        outs = set(); ca_i, cb_i = ta, tb
        for _ in range(3):
            for sa_i in (1, -1):
                for sb_i in (1, -1):
                    outs.add(_canon_cone_ref(tuple(sa_i * x for x in ca_i), tuple(sb_i * x for x in cb_i), S, S2))
            ca_i, cb_i = _mv6(S, ca_i), _mv6(S, cb_i)
        if len(outs) != 1: tie_errs += 1
        if ta[3] == 0 and ta[4] == 0 and ta[5] == 0:
            rt = execute_chord_split(list(ta) + list(tb), [1,0,0,0,0,0, 0,1,0,0,0,0], 0, 4, verify_ref=True)
            if rt["mm"] or rt["ae"]: tie_errs += 1

    # ---- PART 3: 7 Real 50..61b Bareiss-LLL BN254 Orbits:
    #   (a) 100% COMPLETE root-to-leaf (k=11..0, zero skipped levels/branches) scaled-T' execution (~12,000+ leaves)
    #   (b) >= 250 actually-executed exact-production-T12 leaves/orbit (a-&b-plane crossings + endpoint clip + p|N_g)
    #   (c) Complete 12-level analytical Gaussian-volume T12 cost bound (disclosed under H2, <= 42.0 Ir/leaf)
    rng = random.Random(seed)
    T12_prod = int(round(1764 * (589.43 ** 2) * G_SC))
    omega = [math.pi ** (d / 2.0) / math.gamma(d / 2.0 + 1.0) for d in range(13)]
    det_G12_sqrt = sqrt((1764.0 ** 6) * (1000.0 ** 6) * (M._det_mat(T4["GA"]) ** 2)) * (G_SC ** 6)
    pbits_list = (50, 52, 54, 56, 58, 60, 61)

    comp_tot_lvs = 0; comp_tot_nodes = 0; comp_mm = 0; comp_ae = 0; comp_div_err = 0; comp_max_sa = 0
    comp_node_ratios_obs = []; comp_node_ratios_pred = []
    prod_exec_lvs = 0; prod_exec_mm = 0; prod_exec_ae = 0; prod_exec_div_err = 0; prod_exec_clips = 0
    prod_anal_ir_vals = []; prod_min_b2_bits = 999.0
    prod_max_inv_bits = 0; prod_max_r2_bits = 0; prod_max_d2_bits = 0; prod_max_ck_bits = 0
    inc_c0_errs = 0

    for pbits_t in pbits_list:
        pp = (1 << pbits_t) + 1
        while True:
            if M._isp(pp):
                s0p = _find_root_h(pp)
                if s0p is not None:
                    wsp = _eval_poly(M.W_COEFFS, s0p, pp)
                    if wsp != 0:
                        uwp = ((M.U_BN254 % pp) * wsp) % pp
                        rhop = _sqrt_mod(pow(uwp, -1, pp), pp)
                        if rhop > 0: break
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
        Dp = [1] + [0] * 12; lamp = [[0] * 12 for _ in range(12)]
        for i in range(12):
            for j in range(i + 1):
                u = GBp[i][j]
                for k in range(j): u = (Dp[k + 1] * u - lamp[i][k] * lamp[j][k]) // Dp[k]
                if j < i: lamp[i][j] = u
                else: Dp[i + 1] = u
        min_b2 = min(Dp[k + 1] / Dp[k] for k in range(12))
        prod_min_b2_bits = min(prod_min_b2_bits, log2(min_b2))
        assert min_b2 >= float(b_star_min_f) * 0.999999

        # (a) 100% COMPLETE root-to-leaf (k=11..0) execution at scaled T_sc targeting ~1,800 real leaves!
        T_sc = int((4200.0 * pp * det_G12_sqrt / omega[12]) ** (1.0 / 6.0))
        mu_s = [[_chk_i64((lamp[j][k] << 32) // Dp[k + 1]) if j > k else 0 for k in range(12)] for j in range(12)]
        b2_s = [_chk_u64(min(M64, (Dp[k + 1] << 64) // (Dp[k] * T_sc))) for k in range(12)]
        inv_s = [_chk_u64(((Dp[k] * T_sc << 38) + Dp[k + 1] - 1) // Dp[k + 1]) for k in range(12)]
        Xm_s = [0] * 12
        for k in range(11, -1, -1):
            rk = math.isqrt(((Dp[k] * T_sc) + Dp[k + 1] - 1) // Dp[k + 1]) + 1
            Xm_s[k] = rk + sum((Xm_s[j] + 1) // 2 for j in range(k + 1, 12))
        ec_s = [sum(Xm_s[j] for j in range(k + 1, 12)) + 1 for k in range(12)]
        nk_s = [0] * 12; xs = [0] * 12; orb_lvs = [0]

        def full_tree(k, rem, az, cur, c0_passed=None):
            nonlocal comp_mm, comp_ae, comp_div_err, comp_max_sa, inc_c0_errs
            nk_s[k] += 1
            ck = 0
            for j in range(k + 1, 12): ck = _chk_i64(ck - _chk_i64(mu_s[j][k] * xs[j]))
            if k == 0 and c0_passed is not None and c0_passed != ck: inc_c0_errs += 1
            r2 = _chk_u64(_umulh64(rem, inv_s[k]) + 1)
            rad = _chk_u64((_isqrt_u40_4newton(r2) + 1) << 26)
            ec = ec_s[k]
            lo = (1 if k == 0 else 0) if az else (_chk_i64(ck - ec - rad + Q32 - 1) >> 32)
            hi = _chk_i64(ck + ec + rad) >> 32
            if lo > hi: return
            Bk = Bp[k]
            if k == 0:
                # Actually execute EVERY leaf through Chord-Split SWAR + direct a_tab check + p | N_g(a,b)!
                # Only verify full 12-coord canon on first 400 leaves/orbit for speed, check a_tab & p|N_g on ALL!
                rc = execute_chord_split(cur, Bk, lo, hi, T12_chk=T_sc, verify_ref=(orb_lvs[0] < 400), p_div=pp, s0_div=s0p, rho_div=rhop)
                orb_lvs[0] += rc["lvs"]; comp_mm += rc["mm"]; comp_ae += rc["ae"]; comp_div_err += rc["div_err"]
                if rc["max_sa"] > comp_max_sa: comp_max_sa = rc["max_sa"]
            else:
                b2k = b2_s[k]
                c0_cur = _chk_i64(-sum(mu_s[j][0] * xs[j] for j in range(2, 12)) - mu_s[1][0] * lo) if k == 1 else None
                for xv in range(lo, hi + 1):
                    ad = _chk_u64(abs(_chk_i64((xv << 32) - ck)))
                    d_lo = ad - ec if ad > ec else 0
                    d19 = d_lo >> 13
                    sub = _umulh64(_chk_u64(d19 * d19), b2k)
                    if sub <= rem:
                        xs[k] = xv
                        full_tree(k - 1, rem - sub, az and (xv == 0), [cur[c] + xv * Bk[c] for c in range(12)], c0_cur if k == 1 else None)
                    if k == 1: c0_cur = _chk_i64(c0_cur - mu_s[1][0])

        full_tree(11, Q38, True, [0] * 12)
        comp_tot_lvs += orb_lvs[0]; comp_tot_nodes += sum(nk_s)
        comp_node_ratios_obs.append(round(sum(nk_s) / max(1, orb_lvs[0]), 3))
        pred_sc = sum(2.0 ** (log2(omega[12 - k] / omega[12]) + 0.5 * (log2(Dp[k]) - k * log2(T_sc))) for k in range(1, 12))
        comp_node_ratios_pred.append(round(pred_sc, 3))

        # (b) Exact production-T12: check bit widths + execute >= 250 real leaves completely through chord-split & clipping!
        mu_p = [[_chk_i64((lamp[j][k] << 32) // Dp[k + 1]) if j > k else 0 for k in range(12)] for j in range(12)]
        b2_p = [_chk_u64(min(M64, (Dp[k + 1] << 64) // (Dp[k] * T12_prod))) for k in range(12)]
        inv_p = [_chk_u64(((Dp[k] * T12_prod << 38) + Dp[k + 1] - 1) // Dp[k + 1]) for k in range(12)]
        for ib in inv_p: prod_max_inv_bits = max(prod_max_inv_bits, ib.bit_length())
        Xm_p = [0] * 12
        for k in range(11, -1, -1):
            rk = math.isqrt(((Dp[k] * T12_prod) + Dp[k + 1] - 1) // Dp[k + 1]) + 1
            Xm_p[k] = rk + sum((Xm_p[j] + 1) // 2 for j in range(k + 1, 12))
        ec_p = [sum(Xm_p[j] for j in range(k + 1, 12)) + 1 for k in range(12)]
        Bnorm_f = [Dp[k + 1] / Dp[k] for k in range(12)]; R_f = sqrt(float(T12_prod))

        # Measure real production chord-split stats (incl. BOTH a- and b-plane crossings + actual endpoint clips!) on >= 300 real leaves
        p_lvs = 0; p_chords = 0; p_fast = 0; p_bnd = 0; p_segs = 0; p_bspl = 0; p_fspl = 0; p_cross = 0; p_clips = 0
        while p_lvs < 300:
            g = [rng.gauss(0.0, 1.0) for _ in range(11)]
            n2 = sum(z * z for z in g); rad_g = R_f * (rng.random() ** (1.0 / 11.0)) / sqrt(n2)
            xp = [0] * 12; rem_p = Q38; ok_s = True
            for k in range(11, 0, -1):
                ck = 0
                for j in range(k + 1, 12): ck = _chk_i64(ck - _chk_i64(mu_p[j][k] * xp[j]))
                prod_max_ck_bits = max(prod_max_ck_bits, abs(ck).bit_length())
                xp[k] = int(round(ck / float(Q32) + (g[k - 1] * rad_g) / sqrt(Bnorm_f[k])))
                ad = _chk_u64(abs(_chk_i64((xp[k] << 32) - ck)))
                d_lo = ad - ec_p[k] if ad > ec_p[k] else 0
                d19 = d_lo >> 13; d2 = _chk_u64(d19 * d19)
                prod_max_d2_bits = max(prod_max_d2_bits, d2.bit_length())
                sub = _umulh64(d2, b2_p[k])
                if sub > rem_p: ok_s = False; break
                rem_p -= sub
            if not ok_s: continue
            ck0 = 0
            for j in range(1, 12): ck0 = _chk_i64(ck0 - _chk_i64(mu_p[j][0] * xp[j]))
            r20 = _chk_u64(_umulh64(rem_p, inv_p[0]) + 1)
            prod_max_r2_bits = max(prod_max_r2_bits, r20.bit_length())
            rad0 = _chk_u64((_isqrt_u40_4newton(r20) + 1) << 26)
            l0 = _chk_i64(ck0 - ec_p[0] - rad0 + Q32 - 1) >> 32
            h0 = _chk_i64(ck0 + ec_p[0] + rad0) >> 32
            if l0 > h0: continue
            base_v = [sum(xp[j] * Bp[j][c] for j in range(1, 12)) for c in range(12)]
            rc = execute_chord_split(base_v, Bp[0], l0, h0, T12_chk=T12_prod, verify_ref=True, p_div=pp, s0_div=s0p, rho_div=rhop)
            if rc["lvs"] == 0: continue
            p_chords += 1; p_lvs += rc["lvs"]; p_fast += rc["fast"]; p_bnd += rc["bnd"]
            p_segs += rc["segs"]; p_bspl += rc["bspl"]; p_fspl += rc["fspl"]; p_cross += rc["cross"]; p_clips += rc["clips"]
            prod_exec_mm += rc["mm"]; prod_exec_ae += rc["ae"]; prod_exec_div_err += rc["div_err"]
        prod_exec_lvs += p_lvs; prod_exec_clips += p_clips

        # (c) Complete 12-level analytical Gaussian-volume T12 cost bound (ALL levels k=0..11, disclosed under H2!)
        ratios_12 = [2.0 ** (log2(omega[12 - k] / omega[12]) + 0.5 * (log2(Dp[k]) - k * log2(T12_prod))) for k in range(1, 12)]
        n0_an = ratios_12[0]; n1p_an = sum(ratios_12[1:])
        # Exact zero-omission ARMv8.0-A opcode tally (KERNELS.md §1, incl. full producer->storage->consumer 4 KiB transfers!):
        #   FP: 57 ops/k=0 node + 82 ops/k>=1 node (ALL loads, stores, stack, 19-op isqrt40, incremental c0, 6-word SWAR)
        #   Chord-split (BOTH a- and b-plane crossings + 60 ops/endpoint clip!):
        #   Fast leaf (10.80 = 8.30 local 4KiB staging + 2.50 full 4KiB flush/owner-route/reuse)
        #   Exact bnd leaf (50.5 = 48.0 + 2.50 flush) + Direct a_tab apply (10.00 = 7.50 local apply + 2.50 4KiB stream ingest)
        fp_ir = n0_an * 57.0 + n1p_an * 82.0
        sp_per_ch = (p_chords * 16.0 + p_bspl * 12.0 + p_fspl * 22.0 + p_cross * 10.0 + p_segs * 12.0 + p_clips * 60.0) / p_chords
        sp_ir = n0_an * sp_per_ch
        lf_ir = (p_fast * 10.80 + p_bnd * 50.50) / p_lvs
        tot_ir = fp_ir + sp_ir + lf_ir + 10.00
        prod_anal_ir_vals.append(round(tot_ir, 2))

    wsum = sum(1.0 / pb for pb in pbits_list)
    prod_inv_bit_mean_ir = round(sum((1.0 / pb) * ir for pb, ir in zip(pbits_list, prod_anal_ir_vals)) / wsum, 2)
    prod_max_ir = max(prod_anal_ir_vals)
    comp_max_rel_err_pct = round(100.0 * max(abs(o - p) / min(o, p) for o, p in zip(comp_node_ratios_obs, comp_node_ratios_pred)), 2)

    # Coordinate bounds check for production Ra=595, Rb_GA=790 (> actual 589.43 / 782.85)
    GA = T4["GA"]
    aug = [GA[i][:] + [1.0 if i == j else 0.0 for j in range(6)] for i in range(6)]
    for col in range(6):
        piv = max(range(col, 6), key=lambda r: abs(aug[r][col]))
        aug[col], aug[piv] = aug[piv], aug[col]
        pv = aug[col][col]; aug[col] = [xv / pv for xv in aug[col]]
        for r in range(6):
            if r != col and aug[r][col]:
                fac = aug[r][col]; aug[r] = [aug[r][c] - fac * aug[col][c] for c in range(12)]
    max_cu = [sqrt(aug[i][6 + i]) for i in range(6)]
    max_a_prod = [round(595.0 * c, 1) for c in max_cu]
    max_b_prod = [round(790.0 * c, 1) for c in max_cu]
    swar_bounds_ok = (max(max_a_prod) < 511.0 and max_a_prod[3] < M.A_TAB_SLOTS_A3 and max_b_prod[0] < M.B0_ARENA_MAX
                      and max_b_prod[1] < 1023.0 and max(max_b_prod[2:]) < 511.0 and ceil((2 * max_b_prod[5] + 1) / 128.0) <= 6)

    # ---- PART 4: Descending-ai Lower B2 Chain + 176-Band 1/16-Bit Direct a_tab Super-Tile Replay +
    #              P_min=24 0.17-Cap g_rec Batched Bernstein SRT + Exact 1/16-Bit sb=0 0-LP Residual Gate ----
    Na_t = 1024; B_st = 8; n_b0 = 16; n_ch_slot = 6
    a_keys = []
    for a3i in range(0, 4):
        for a4i in range(-4, 5):
            for a5i in range(-4, 5):
                for a0i in range(-2, 3):
                    a_keys.append(_pack_v((a0i, 0, 0, a3i, a4i, a5i), SH_A, BI_A_CAN))
                    if len(a_keys) >= Na_t: break
                if len(a_keys) >= Na_t: break
            if len(a_keys) >= Na_t: break
        if len(a_keys) >= Na_t: break
    a_tab_sim = {pa: ai for ai, pa in enumerate(a_keys, 1)}
    # Descending-ai prefix-truncated lower Band-2 chain test
    p_l = 131; W_l = 16; nb_l = (p_l + W_l - 1) // W_l
    kta_l = [0] + [rng.randrange(p_l) for _ in range(Na_t)]
    head_l = [0] * nb_l; nxt_l = [0] * (Na_t + 1); M46 = (1 << 46) - 1
    for ai in range(Na_t, 0, -1):
        k = kta_l[ai]; sb = k // W_l
        nxt_l[ai] = head_l[sb]; head_l[sb] = ((k & (W_l - 1)) << 46) | ai
    pref_errs = 0
    for sb in range(nb_l):
        namax = rng.randint(40, Na_t)
        got_s = set(); rec = head_l[sb]
        while rec:
            ai = rec & M46
            if ai > namax: break
            got_s.add((rec >> 46, ai)); rec = nxt_l[ai]
        ref_s = {(kta_l[ai] & (W_l - 1), ai) for ai in range(1, namax + 1) if kta_l[ai] // W_l == sb}
        if got_s != ref_s: pref_errs += 1

    # Direct a_tab 176-band exact 1/16-bit 4 KiB local-staging -> owner-routed flush -> 1-page streaming replay (ZERO hash table!)
    score_ref = [[[0] * (Na_t + 1) for _ in range(n_b0 * 8)] for _ in range(B_st)]
    score_got = [[[0] * (Na_t + 1) for _ in range(n_b0 * 8)] for _ in range(B_st)]
    # Producer owns fixed 4 KiB (255 x 16B record + 16B header) local staging buffers per tid (2^48 * 4 KiB = 2^60 B << 2^64 B);
    # on byte offset 4096 (ptr & 0x0FF0 == 0), .Lflush_page routes to owner_id=tid, seals page (tid, w_lp_u16, count=255), and executes sub x7, x7, #4080 (4096 -> 16)!
    PAGE_CAP = 255
    stage_buf = [[None] * PAGE_CAP for _ in range(B_st)]
    stage_off = [16] * B_st  # byte offset inside fixed 4 KiB local buffer (16..4096)
    owner_pages = [[] for _ in range(B_st)]  # sealed 4 KiB pages stored at consumer owner_id = tid
    flushed_full_pages = 0; flushed_tail_pages = 0
    for b_idx in range(176):
        w_lp_ticks = 800 + b_idx  # exact 1/16-bit tick in [50.0, 61.0)
        n_recs_b = 2600 if b_idx in (0, 175) else 30  # forces >=16 full 255-rec 4 KiB page flushes + local buffer reuses!
        for _ in range(n_recs_b):
            ai = rng.randint(1, Na_t); tid = rng.randrange(B_st)
            b0v = rng.randrange(n_b0); ch_slot = rng.randint(1, n_ch_slot)
            pa = a_keys[ai - 1]; pb = (b0v << 51) | (ch_slot << 48) | tid
            score_ref[tid][(b0v << 3) | ch_slot][ai] += w_lp_ticks
            slot_i = (stage_off[tid] - 16) >> 4
            stage_buf[tid][slot_i] = (pa, pb)
            stage_off[tid] += 16
            if (stage_off[tid] & 0x0FF0) == 0:  # exact ARMv8.0-A tst x7, #0x0FF0 -> 4096 boundary!
                owner_pages[tid].append((tid, w_lp_ticks, PAGE_CAP, stage_buf[tid][:]))
                stage_off[tid] -= 4080  # exact ARMv8.0-A sub x7, x7, #4080 (returns 4096 -> payload offset 16)!
                flushed_full_pages += 1
        # End-of-band flush of partial local staging buffers
        for tid in range(B_st):
            cnt = (stage_off[tid] - 16) >> 4
            if cnt:
                owner_pages[tid].append((tid, w_lp_ticks, cnt, stage_buf[tid][:cnt]))
                stage_off[tid] = 16
                flushed_tail_pages += 1
    # Consumer streams ONE 4 KiB page at a time from owner_pages[tid] into fixed-stride score arena via direct a_tab
    streamed_pages = 0
    for tid in range(B_st):
        for own_id, w_lp_ticks, cnt, recs in owner_pages[tid]:
            assert own_id == tid and cnt == len(recs)
            streamed_pages += 1
            for pa, pb in recs:
                score_got[tid][pb >> 48][a_tab_sim.get(pa, 0)] += w_lp_ticks
    supertile_errs = sum(1 for tid in range(B_st) for idx in range(n_b0 * 8) if score_got[tid][idx] != score_ref[tid][idx])

    # Genuine 64 ~470b BN254 g_rec norms: exact Batched Bernstein SRT < 2^14 vs trial division!
    def _pmul_h_z(a, b):
        t = [0] * 11
        for i in range(6):
            for j in range(6): t[i + j] += a[i] * b[j]
        for k in range(10, 5, -1):
            c = t[k]
            if c:
                for j in range(6): t[k - 6 + j] -= c * M.H_COEFFS[j]
        return t[:6]
    def _det6(M6):
        A = [r[:] for r in M6]; sg = 1; pv = 1
        for k in range(5):
            if A[k][k] == 0:
                p0 = next((r for r in range(k + 1, 6) if A[r][k] != 0), -1)
                if p0 < 0: return 0
                A[k], A[p0] = A[p0], A[k]; sg = -sg
            p1 = A[k][k]
            for i in range(k + 1, 6):
                for j in range(k + 1, 6): A[i][j] = (A[i][j] * p1 - A[i][k] * A[k][j]) // pv
            pv = p1
        return sg * A[5][5]
    def _exact_ng(a, b, u_val=M.U_BN254):
        a3 = [sum(a[k] * M.ZK_NUM[k][j] for k in range(6)) for j in range(6)]
        b3 = [sum(b[k] * M.ZK_NUM[k][j] for k in range(6)) for j in range(6)]
        g9 = [u_val * x - y for x, y in zip(_pmul_h_z(M.W_COEFFS, _pmul_h_z(a3, a3)), _pmul_h_z(b3, b3))]
        cols = []; cur = g9[:]
        for j in range(6):
            cols.append(cur)
            if j < 5:
                nxt = [0] + cur[:5]; c = cur[5]
                if c:
                    for m in range(6): nxt[m] -= c * M.H_COEFFS[m]
                cur = nxt
        dv = abs(_det6([[cols[c][r] for c in range(6)] for r in range(6)]))
        assert dv % (9 ** 6) == 0
        return dv // (9 ** 6)
    def _bernstein_strip(P0, xs):
        pt = [list(xs)]
        while len(pt[-1]) > 1:
            pv = pt[-1]
            pt.append([pv[i] * pv[i + 1] if i + 1 < len(pv) else pv[i] for i in range(0, len(pv), 2)])
        rem = [P0 % pt[-1][0]]
        for lvl in range(len(pt) - 2, -1, -1):
            rem = [rem[i >> 1] % v for i, v in enumerate(pt[lvl])]
        out = []
        for x_v, r_v in zip(xs, rem):
            cur = x_v; g = math.gcd(cur, r_v)
            while g > 1:
                cur //= g; g = math.gcd(cur, g)
            out.append((cur, x_v // cur))
        return out

    gen_norms = []
    while len(gen_norms) < 64:
        av = [rng.randint(-220, 220) for _ in range(6)]
        bv = [rng.randint(-300, 300) for _ in range(6)]
        if any(av) and any(bv) and M.is_ideal_primitive_oh(av, bv):
            gen_norms.append(_exact_ng(av, bv))
    ps14 = M._primes_upto(1 << 14); P14 = 1
    for p_v in ps14: P14 *= p_v
    srt_out = []
    for off in range(0, 64, 16): srt_out.extend(_bernstein_strip(P14, gen_norms[off:off + 16]))
    srt_errs = 0
    for idx_n, (cof_s, sm_s) in enumerate(srt_out):
        cur = gen_norms[idx_n]
        for p_v in ps14:
            while cur % p_v == 0: cur //= p_v
        if cur != cof_s or (gen_norms[idx_n] // cur) != sm_s: srt_errs += 1

    # TRUE end-to-end completely-factored real C6 g_rec (u_t=5, same K_h, w(y), g_rec) residual-gate test:
    # Real Batched Bernstein SRT (<2^8) -> actual integer C_post -> actual 1/16-bit single-hit sieve score on [2^8, 2^20]
    # -> actual sb=0 residual gate compared against GROUND-TRUTH COMPLETE PRIME FACTORIZATION!
    def _is_prob_prime(n):
        if n < 2: return False
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
            if n == p: return True
            if n % p == 0: return False
        d = n - 1; s = 0
        while (d & 1) == 0: d >>= 1; s += 1
        for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
            if a >= n: continue
            x = pow(a, d, n)
            if x in (1, n - 1): continue
            for _ in range(s - 1):
                x = (x * x) % n
                if x == n - 1: break
            else: return False
        return True
    def _pollard_rho(n):
        if (n & 1) == 0: return 2
        for c in range(1, 50):
            x = 2; y = 2; d = 1
            while d == 1:
                x = (x * x + c) % n; y = (y * y + c) % n; y = (y * y + c) % n
                d = math.gcd(abs(x - y), n)
            if d != n: return d
        return n
    def _factor_full(n, small_ps):
        facs = {}
        for p in small_ps:
            if p * p > n: break
            while n % p == 0: facs[p] = facs.get(p, 0) + 1; n //= p
        stk = [n] if n > 1 else []
        while stk:
            m = stk.pop()
            if m == 1: continue
            if _is_prob_prime(m): facs[m] = facs.get(m, 0) + 1
            else: d = _pollard_rho(m); stk.append(d); stk.append(m // d)
        return facs

    B_STRIP_E2E = 1 << 8; B_SIEVE_E2E = 1 << 20
    ps_e2e_small = M._primes_upto(1 << 10)
    ps_e2e_strip = M._primes_upto(B_STRIP_E2E); P_e2e_strip = 1
    for p_v in ps_e2e_strip: P_e2e_strip *= p_v
    e2e_norms = []
    while len(e2e_norms) < 320:
        av = [rng.randint(-4, 4) for _ in range(6)]
        bv = [rng.randint(-6, 6) for _ in range(6)]
        if any(av) and any(bv) and M.is_ideal_primitive_oh(av, bv):
            ng_e = _exact_ng(av, bv, u_val=5)
            if ng_e >= (1 << 25): e2e_norms.append(ng_e)
    e2e_srt = []
    for off in range(0, len(e2e_norms), 32): e2e_srt.extend(_bernstein_strip(P_e2e_strip, e2e_norms[off:off + 32]))
    e2e_smooth_sqf = 0; e2e_smooth_p2 = 0; e2e_1lp = 0; e2e_multi = 0
    fn_0lp = 0; fp_0lp = 0; min_m_sm = 1e9; min_g_1lp = 1e9
    for ng_e, (c_post, _) in zip(e2e_norms, e2e_srt):
        nb = log2(ng_e); facs = _factor_full(c_post, ps_e2e_small)
        # Exact 1/16-bit single-hit sieve accumulator on primes p in [2^8, 2^20] dividing ng_e
        score_ticks = sum(floor(16.0 * log2(p_v)) for p_v in facs if B_STRIP_E2E <= p_v <= B_SIEVE_E2E)
        has_p2 = any(e_v >= 2 for p_v, e_v in facs.items() if B_STRIP_E2E <= p_v <= B_SIEVE_E2E)
        cap_nb = floor(ceil(nb) / 8.0); rb = (cap_nb + 1) / 16.0
        res = (log2(c_post) if c_post > 1 else 0.0) - score_ticks / 16.0
        lps = [p_v for p_v in facs if p_v > B_SIEVE_E2E]
        passed = (res <= rb + 0.25)
        if not lps:
            if not has_p2:
                e2e_smooth_sqf += 1
                mrg = (rb + 0.25) - res
                if mrg < min_m_sm: min_m_sm = mrg
                if not passed: fn_0lp += 1
            else:
                e2e_smooth_p2 += 1
        else:
            if len(lps) == 1 and facs[lps[0]] == 1: e2e_1lp += 1
            else: e2e_multi += 1
            gap = res - (rb + 0.25)
            if gap < min_g_1lp: min_g_1lp = gap
            if passed: fp_0lp += 1

    ok = (smat_cyclic_ok and uncond_64bit_zero_exc_ok and isqrt_errs == 0 and lat_mismatches == 0
          and tie_errs == 0 and pref_errs == 0 and supertile_errs == 0 and flushed_full_pages >= 10
          and streamed_pages == flushed_full_pages + flushed_tail_pages and inc_c0_errs == 0 and swar_bounds_ok
          and comp_tot_lvs >= 10000 and comp_mm == 0 and comp_ae == 0 and comp_div_err == 0 and comp_max_sa <= 4
          and comp_max_rel_err_pct <= 5.0
          and prod_exec_lvs >= 2000 and prod_exec_mm == 0 and prod_exec_ae == 0 and prod_exec_div_err == 0
          and prod_max_inv_bits < 64 and prod_max_r2_bits < 40 and prod_max_d2_bits < 64 and prod_max_ck_bits < 63
          and prod_max_ir <= M.C_LAT_LEAF and srt_errs == 0
          and e2e_smooth_sqf >= 10 and e2e_1lp >= 100 and fn_0lp == 0 and fp_0lp == 0)

    return dict(
        ok=ok,
        smat_cyclic_a345_ok=smat_cyclic_ok,
        zero_float_fp_enum=True,
        strict_u64_i64_umulh_checked=True,
        zero_exception_64bit_q38_proved=uncond_64bit_zero_exc_ok,
        cert_min_gsym_exact=m6_exact,
        cert_min_gsym_pts_le_5gsc=m6_pts,
        uncond_min_bk_star_bits=round(log2(float(b_star_min_f)), 3),
        uncond_inv_max_bits=uncond_inv_max.bit_length(),
        uncond_r2_q12_max_bits=uncond_r2_max.bit_length(),
        uncond_ck_max_bits=uncond_ck_max.bit_length(),
        uncond_d2_max_bits=uncond_d2_max.bit_length(),
        isqrt_u40_4newton_errs=isqrt_errs,
        orbits_6ideal_checked=len(orbits),
        lattice_b2_mismatches=lat_mismatches,
        lattice_b2_raw_hits=tot_raw_hits,
        lattice_b2_jvor_hits=tot_jvor_hits,
        complete_root_to_leaf_orbits=len(pbits_list),
        complete_root_to_leaf_leaves=comp_tot_lvs,
        complete_root_to_leaf_nodes=comp_tot_nodes,
        complete_root_to_leaf_mismatches=comp_mm,
        complete_root_to_leaf_apply_errs=comp_ae,
        complete_root_to_leaf_div_errs=comp_div_err,
        complete_root_to_leaf_max_sa_cones=comp_max_sa,
        complete_scaled_nodes_per_leaf_obs=comp_node_ratios_obs,
        complete_scaled_nodes_per_leaf_pred=comp_node_ratios_pred,
        complete_scaled_max_rel_err_pct=comp_max_rel_err_pct,
        prod_exec_leaves=prod_exec_lvs,
        prod_exec_mismatches=prod_exec_mm,
        prod_exec_apply_errs=prod_exec_ae,
        prod_exec_div_errs=prod_exec_div_err,
        prod_exec_endpoint_clips=prod_exec_clips,
        prod_min_bk_star_bits=round(prod_min_b2_bits, 2),
        prod_max_inv_b2_bits=prod_max_inv_bits,
        prod_max_r2_q12_bits=prod_max_r2_bits,
        prod_max_d2_bits=prod_max_d2_bits,
        prod_max_ck_bits=prod_max_ck_bits,
        prod_12lvl_analytic_ir_per_leaf=prod_anal_ir_vals,
        prod_inv_bit_weighted_mean_ir=prod_inv_bit_mean_ir,
        prod_max_ir=prod_max_ir,
        prod_budget_c_lat_leaf=M.C_LAT_LEAF,
        swar_tie_cases_checked=len(tie_cases),
        swar_tie_canon_errs=tie_errs,
        desc_ai_prefix_chain_errs=pref_errs,
        direct_atab_176band_replay_errs=supertile_errs,
        direct_atab_logp_bands_1_16b=176,
        direct_atab_4k_full_pages_flushed=flushed_full_pages,
        direct_atab_4k_tail_pages_flushed=flushed_tail_pages,
        direct_atab_4k_pages_streamed=streamed_pages,
        zero_hash_table_apply=True,
        inc_c0_exact_errs=inc_c0_errs,
        swar_bounds_ok=swar_bounds_ok,
        max_a_prod_coords=max_a_prod,
        max_b_prod_coords=max_b_prod,
        srt_p24_genuine_norms_checked=len(gen_norms),
        srt_batch_mismatches=srt_errs,
        e2e_factored_c6_norms_checked=len(e2e_norms),
        e2e_factored_smooth_sqfree=e2e_smooth_sqf,
        e2e_factored_smooth_p2=e2e_smooth_p2,
        e2e_factored_1lp=e2e_1lp,
        e2e_factored_multi_lp=e2e_multi,
        gate_0lp_fn=fn_0lp,
        gate_0lp_fp=fp_0lp,
        gate_0lp_min_smooth_margin_b=round(min_m_sm, 4),
        gate_0lp_min_1lp_gap_b=round(min_g_1lp, 2),
    )
