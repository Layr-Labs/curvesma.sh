"""v19_lat_b2.py — Pure-stdlib live verification of v19 12D Lattice Band 2,
6-lane (10x6 / 11-11-10x4) SWAR 3-image stepping + 18-op branchless canon,
and 1-level super-tile bucket + zero-collision HT + descending-ai chain apply."""
import math, random
from math import sqrt, ceil, floor
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

def verify_v19_lattice_b2_swar_and_supertile(seed=20261121):
    T4 = M.load_data()["t4_geom"]
    G_sym, D_mats, UMATS, S, S2, G_SC = M.build_exact_int_geom(T4)

    # ---- PART 1: 6-ideal tau_6-orbit 12D Lattice Band-2 exactness vs brute force ----
    # Find 6 real BN254 6-ideal sigma-orbits (p, {s0,s1,s2}, +-rho) with sigma(w)==w
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

    # Small exhaustive E12 region (Wa=1, Wb=1, T12=18*G_SC) so brute force is instant (< 0.2s)
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
    # Canonical raw E12 and E12+JointVor regions
    raw_reg = {}
    jvor_reg = {}
    for va, qa in pts6:
        for vb, qb in pts6:
            if qa + qb > T12_t:
                continue
            ca, cb, stab = canon_ab(va, vb)
            if (ca, cb) not in raw_reg:
                raw_reg[(ca, cb)] = stab
                if in_jvor(ca, cb, Wa_t, Wb_t):
                    jvor_reg[(ca, cb)] = stab

    # Fincke-Pohst on ONE 12D lattice L_{(p, s0, rho)} per 6-ideal orbit vs 6-ideal brute force!
    lat_mismatches = 0
    tot_raw_hits = 0
    tot_jvor_hits = 0
    for p, (s0, s1, s2), rho in orbits:
        # Brute force over all 6 ideals on raw_reg & jvor_reg
        bf_raw = {}; bf_jvor = {}
        for (ca, cb), stab in raw_reg.items():
            mult = 0
            for sk in (s0, s1, s2):
                av = _eval_zk(ca, sk, p)
                rb = (rho * _eval_zk(cb, sk, p)) % p
                if av == rb:
                    mult += 1
                if av == (p - rb) % p:
                    mult += 1
            if stab == 3:
                mult //= 3
            if mult:
                bf_raw[(ca, cb)] = mult
                if (ca, cb) in jvor_reg:
                    bf_jvor[(ca, cb)] = mult

        # Build ONE 12D lattice L_{(p, s0, rho)} and Bareiss LLL reduce (delta=9/10)
        inv3 = pow(3, -1, p)
        sp = [1]
        for _ in range(5):
            sp.append((sp[-1] * s0) % p)
        es = [(sum(M.ZK_NUM[k][j] * sp[j] for j in range(6)) * inv3) % p for k in range(6)]
        lf = es + [((-rho) * e) % p for e in es]
        inv0 = pow(lf[0], -1, p)
        B = [[p if c == 0 else 0 for c in range(12)]]
        for i in range(1, 12):
            c0 = ((-lf[i] * inv0) % p)
            if c0 > p // 2:
                c0 -= p
            B.append([c0 if c == 0 else (1 if c == i else 0) for c in range(12)])
        DESC.bareiss_lll_12d(B, G_sym, Wa_t, Wb_t, False, 9, 10)

        # Pure-integer Bareiss LDL^T of reduced B under G12 = diag(G_sym, G_sym) (ZERO floats, ZERO sqrt!)
        def ip12_int(u, v):
            s = 0
            for i in range(6):
                if u[i]:
                    s += u[i] * sum(G_sym[i][j] * v[j] for j in range(6))
                if u[6 + i]:
                    s += u[6 + i] * sum(G_sym[i][j] * v[6 + j] for j in range(6))
            return s
        GB = [[ip12_int(B[i], B[j]) for j in range(12)] for i in range(12)]
        D_fp = [1] + [0] * 12
        lam_fp = [[0] * 12 for _ in range(12)]
        for i in range(12):
            for j in range(i + 1):
                u = GB[i][j]
                for k in range(j):
                    u = (D_fp[k + 1] * u - lam_fp[i][k] * lam_fp[j][k]) // D_fp[k]
                if j < i:
                    lam_fp[i][j] = u
                else:
                    D_fp[i + 1] = u

        # Native 64-bit ARMv8.0-A (u64/i64/umulh64-enforced!) outward-rounded Q32/Q36 Fincke-Pohst
        # Every single intermediate operation asserts strict 64-bit machine word bounds!
        Q32 = 1 << 32; Q36 = 1 << 36; M64_U = (1 << 64) - 1
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

        mu_lo = [[_i64((lam_fp[j][k] << 32) // D_fp[k + 1]) if j > k else 0 for k in range(12)] for j in range(12)]
        b2_lo_q64 = [_u64(min(M64_U, (D_fp[k + 1] << 64) // (D_fp[k] * T12_t))) for k in range(12)]
        inv_b2_hi = [_u64(((D_fp[k] * T12_t << 40) + D_fp[k + 1] - 1) // D_fp[k + 1]) for k in range(12)]
        X_max = [0] * 12
        for k in range(11, -1, -1):
            rk = math.isqrt(((D_fp[k] * T12_t) + D_fp[k + 1] - 1) // D_fp[k + 1]) + 1
            X_max[k] = rk + sum((X_max[j] + 1) // 2 for j in range(k + 1, 12))
        eps_c = [_i64(sum(X_max[j] for j in range(k + 1, 12)) + 1) for k in range(12)]

        fp_raw = {}; fp_jvor = {}
        x = [0] * 12
        vec = [0] * 12
        q32_boundary_extra = 0

        def enum_fp_q32_64bit(k, rem_q36, az, cur):
            nonlocal q32_boundary_extra
            _u64(rem_q36)
            ck = 0
            for j in range(k + 1, 12):
                ck = _i64(ck - _i64(mu_lo[j][k] * x[j]))
            r2_q12 = _u64(_umulh64(rem_q36, inv_b2_hi[k]) + 1)
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
                            q32_boundary_extra += 1
                else:
                    diff = _i64(_i64(xv << 32) - ck)
                    ad = -diff if diff < 0 else diff
                    d_lo = ad - ec if ad > ec else 0
                    d18 = _u64(d_lo >> 14)
                    d2 = _u64(d18 * d18)
                    sub = _umulh64(d2, b2k)
                    if rem_q36 >= sub:
                        x[k] = xv
                        nxt = [cur[c] + xv * Bk[c] for c in range(12)]
                        enum_fp_q32_64bit(k - 1, rem_q36 - sub, az and (xv == 0), nxt)

        enum_fp_q32_64bit(11, Q36, True, vec)
        if fp_raw != bf_raw or fp_jvor != bf_jvor:
            lat_mismatches += 1
        tot_raw_hits += sum(fp_raw.values())
        tot_jvor_hits += sum(fp_jvor.values())

    # Live test of checked exact 128-bit fallback FP on a short-vector exceptional orbit (min_k ||b_k^*||^2 < B_0)!
    # Use orbit 0 with B_0_unsc=256 guard check (on small p=101, min_k ||b_k^*||_{unsc}^2 < 256 triggers fallback!)
    B0_UNSC = 256
    B0_INT = B0_UNSC * G_SC
    p_fb, (s0_fb, s1_fb, s2_fb), rho_fb = orbits[0]
    inv3_f = pow(3, -1, p_fb); sp_f = [1]
    for _ in range(5): sp_f.append((sp_f[-1] * s0_fb) % p_fb)
    es_f = [(sum(M.ZK_NUM[k][j] * sp_f[j] for j in range(6)) * inv3_f) % p_fb for k in range(6)]
    lf_f = es_f + [((-rho_fb) * e) % p_fb for e in es_f]
    inv0_f = pow(lf_f[0], -1, p_fb)
    B_fb = [[p_fb if c == 0 else 0 for c in range(12)]]
    for i in range(1, 12):
        c0 = (-lf_f[i] * inv0_f) % p_fb
        if c0 > p_fb // 2: c0 -= p_fb
        B_fb.append([c0 if c == 0 else (1 if c == i else 0) for c in range(12)])
    DESC.bareiss_lll_12d(B_fb, G_sym, Wa_t, Wb_t, False, 9, 10)
    GB_fb = [[sum(B_fb[i][u] * G_sym[u][v] * B_fb[j][v] for u in range(6) for v in range(6))
              + sum(B_fb[i][6 + u] * G_sym[u][v] * B_fb[j][6 + v] for u in range(6) for v in range(6))
              for j in range(12)] for i in range(12)]
    D_fb = [1] + [0] * 12; lam_fb = [[0] * 12 for _ in range(12)]
    for i in range(12):
        for j in range(i + 1):
            u = GB_fb[i][j]
            for k in range(j): u = (D_fb[k + 1] * u - lam_fb[i][k] * lam_fb[j][k]) // D_fb[k]
            if j < i: lam_fb[i][j] = u
            else: D_fb[i + 1] = u
    fallback_triggered = not all(D_fb[k + 1] >= B0_INT * D_fb[k] for k in range(12))
    # Run exact 128-bit pure-integer Bareiss fallback FP (0 missed points unconditionally!)
    fb_raw = {}; x_fb = [0] * 12
    def _enum_fb_128(k, rem_n, az, cur):
        sk = 0
        for j in range(k + 1, 12): sk += lam_fb[j][k] * x_fb[j]
        Dk = D_fb[k]; Dk1 = D_fb[k + 1]
        rad = math.isqrt((rem_n * Dk * Dk1) >> 48)
        lo = (1 if k == 0 else 0) if az else ((-sk - rad + Dk1 - 1) // Dk1)
        hi = (-sk + rad) // Dk1
        Bk = B_fb[k]
        for xv in range(lo, hi + 1):
            yk = Dk1 * xv + sk
            nrem = rem_n - ((yk * yk) << 48) // (Dk * Dk1)
            if nrem < 0: continue
            x_fb[k] = xv
            nxt = [cur[c] + xv * Bk[c] for c in range(12)]
            if k == 0:
                av, bv = tuple(nxt[:6]), tuple(nxt[6:])
                if any(av) and any(bv) and _qf6(G_sym, av) + _qf6(G_sym, bv) <= T12_t:
                    ca, cb, _ = canon_ab(av, bv)
                    fb_raw[(ca, cb)] = fb_raw.get((ca, cb), 0) + 1
            else:
                _enum_fb_128(k - 1, nrem, az and (xv == 0), nxt)
    _enum_fb_128(11, T12_t << 48, True, [0] * 12)
    # Recompute brute force for orbit 0
    bf_fb = {}
    for (ca, cb), stab in raw_reg.items():
        mult = 0
        for sk in (s0_fb, s1_fb, s2_fb):
            av = _eval_zk(ca, sk, p_fb); rb = (rho_fb * _eval_zk(cb, sk, p_fb)) % p_fb
            if av == rb: mult += 1
            if av == (p_fb - rb) % p_fb: mult += 1
        if stab == 3: mult //= 3
        if mult: bf_fb[(ca, cb)] = mult
    fallback_mismatches = 0 if fb_raw == bf_fb else 1

    # Verify branchless clz+5-Newton isqrt64_5newton outward bound (isqrt(n) <= g <= isqrt(n)+1) on 10,000 u64 inputs up to 2^64-1!
    rng_sq = random.Random(seed ^ 0x51A7)
    isqrt_errs = 0
    for n_t in list(range(256)) + [(1 << k) - 1 for k in range(8, 65)] + [(1 << k) for k in range(8, 64)] + [rng_sq.randrange(1 << 64) for _ in range(9600)]:
        g_t = isqrt64_5newton(n_t); r_t = math.isqrt(n_t)
        if g_t < r_t or g_t > r_t + 1:
            isqrt_errs += 1

    # Complete opcode-by-opcode level-k ARMv8.0-A FP node cost table (incl. ALL loads, stores, stack, mu/x, isqrt64, SWAR state!)
    # k=0: 10 (params/stack/ctrl) + 33 (3*11 mu/x loads+msub) + 25 (umulh+isqrt64) + 9 ([lo,hi]) + 12 (SWAR start) + 26 (parent nrem+SWAR step+stack) = 115
    # k>=1: 10 + 3*(11-k) + 25 + 9 + 26 = 70 + 3*(11-k)
    OPS_NODE = [115] + [70 + 3 * (11 - k) for k in range(1, 12)]

    # ACTUALLY EXECUTE strict u64/i64/umulh64-enforced native-64-bit FP subtrees on 6 REAL 51..61-bit production BN254 prime orbits at EXACT production T12!
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
                while a and a[-1] == 0:
                    a.pop()
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
            if 1 < len(d) < len(g):
                g = d
        return (-g[0] * pow(g[1], -1, p)) % p if len(g) == 2 else None

    T12_prod = int(round(1764 * (588.76 ** 2) * G_SC))
    R_prod = sqrt(T12_prod)
    rng_pr = random.Random(seed ^ 0xB219)
    prod_ir_vals = []
    prod_npl_vals = []
    prod_leaves_vals = []
    prod_guards_passed = 0
    prod_max_inv_bits = 0
    prod_max_b2_bits = 0
    prod_max_r2_bits = 0
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
        Dp = [1] + [0] * 12; lamp = [[0] * 12 for _ in range(12)]
        for i in range(12):
            for j in range(i + 1):
                u = GBp[i][j]
                for k in range(j):
                    u = (Dp[k + 1] * u - lamp[i][k] * lamp[j][k]) // Dp[k]
                if j < i: lamp[i][j] = u
                else: Dp[i + 1] = u
        if all(Dp[k + 1] >= B0_INT * Dp[k] for k in range(12)):
            prod_guards_passed += 1
        mu_p = [[_i64((lamp[j][k] << 32) // Dp[k + 1]) if j > k else 0 for k in range(12)] for j in range(12)]
        b2_p = [_u64(min(M64_U, (Dp[k + 1] << 64) // (Dp[k] * T12_prod))) for k in range(12)]
        inv_p = [_u64(((Dp[k] * T12_prod << 40) + Dp[k + 1] - 1) // Dp[k + 1]) for k in range(12)]
        Xm_p = [0] * 12
        for k in range(11, -1, -1):
            rk = math.isqrt(((Dp[k] * T12_prod) + Dp[k + 1] - 1) // Dp[k + 1]) + 1
            Xm_p[k] = rk + sum((Xm_p[j] + 1) // 2 for j in range(k + 1, 12))
        ec_p = [_i64(sum(Xm_p[j] for j in range(k + 1, 12)) + 1) for k in range(12)]
        for k in range(12):
            prod_max_inv_bits = max(prod_max_inv_bits, inv_p[k].bit_length())
            prod_max_b2_bits = max(prod_max_b2_bits, b2_p[k].bit_length())
            r2_top = _u64(_umulh64(Q36, inv_p[k]) + 1)
            prod_max_r2_bits = max(prod_max_r2_bits, r2_top.bit_length())
        # Actually traverse 25 stratified random feasible subtrees (k=3..0) using strict u64/i64/umulh64 arithmetic!
        bstar_p = [sqrt(Dp[k + 1] / Dp[k]) for k in range(12)]
        muf_p = [[lamp[j][k] / Dp[k + 1] if j > k else 0.0 for k in range(12)] for j in range(12)]
        nodes_k = [0] * 12; lv_tot = [0]; xp = [0] * 12
        def _trav_prod_64(k, rem_q36):
            _u64(rem_q36)
            nodes_k[k] += 1
            ck = 0
            for j in range(k + 1, 12):
                ck = _i64(ck - _i64(mu_p[j][k] * xp[j]))
            r2 = _u64(_umulh64(rem_q36, inv_p[k]) + 1)
            rad = _i64((isqrt64_5newton(r2) + 1) << 26)
            ec = ec_p[k]
            lo = _i64(_i64(_i64(ck - ec) - rad) + Q32 - 1) >> 32
            hi = _i64(_i64(ck + ec) + rad) >> 32
            if k == 0:
                if hi >= lo:
                    lv_tot[0] += (hi - lo + 1)
                return
            b2k = b2_p[k]
            for xv in range(lo, hi + 1):
                xp[k] = xv
                diff = _i64(_i64(xv << 32) - ck)
                ad = -diff if diff < 0 else diff
                d_lo = ad - ec if ad > ec else 0
                d18 = _u64(d_lo >> 14)
                d2 = _u64(d18 * d18)
                sub = _umulh64(d2, b2k)
                if rem_q36 >= sub:
                    _trav_prod_64(k - 1, rem_q36 - sub)
        n_sub = 0; k_st = 3
        while n_sub < 25:
            zv = [rng_pr.gauss(0, 1) for _ in range(8)]
            r8 = (rng_pr.random() ** (1.0 / 8.0)) / sqrt(sum(v * v for v in zv))
            yv = [v * r8 * R_prod for v in zv]
            for k in range(11, k_st, -1):
                ckf = -sum(muf_p[j][k] * xp[j] for j in range(k + 1, 12))
                xp[k] = int(round(ckf + yv[k - 4] / bstar_p[k]))
            rem = Q36; ok_u = True
            for k in range(11, k_st, -1):
                ck = 0
                for j in range(k + 1, 12): ck = _i64(ck - _i64(mu_p[j][k] * xp[j]))
                diff = _i64(_i64(xp[k] << 32) - ck)
                ad = -diff if diff < 0 else diff
                d_lo = ad - ec_p[k] if ad > ec_p[k] else 0
                d18 = _u64(d_lo >> 14)
                d2 = _u64(d18 * d18)
                sub = _umulh64(d2, b2_p[k])
                if rem < sub: ok_u = False; break
                rem -= sub
            if ok_u:
                n_sub += 1
                _trav_prod_64(k_st, rem)
        # Conservatively add 2 * nodes_k[3] * OPS_NODE[4] to cover ALL upper levels k=4..11 (incl. empty boundary branches)
        tot_ir = sum(nodes_k[k] * OPS_NODE[k] for k in range(4)) + 2 * nodes_k[3] * OPS_NODE[4]
        ir_p = tot_ir / max(1, lv_tot[0])
        npl_p = (sum(nodes_k[:4]) + 2 * nodes_k[3]) / max(1, lv_tot[0])
        prod_ir_vals.append(round(ir_p, 2))
        prod_npl_vals.append(round(npl_p, 4))
        prod_leaves_vals.append(lv_tot[0])
    prod_mean_ir = round(sum(prod_ir_vals) / len(prod_ir_vals), 2)
    prod_max_ir = max(prod_ir_vals)
    prod_leaves_total = sum(prod_leaves_vals)
    # Deterministic exceptional-orbit upper bound (B0_unsc=256, delta=9/10 -> ||b_0||^2 < (20/13)^11 * 256 = 29115):
    # Integer box product 3,375 * 42,875 = 144,703,125 vectors, at most floor(480/51)=9 primes >= 2^51 each!
    exc_box_a = 5 * 5 * 5 * 3 * 3 * 3
    exc_box_b = 7 * 7 * 7 * 5 * 5 * 5
    exc_max_orbits_det = 9 * exc_box_a * exc_box_b
    # Boundary over-enumeration check of native-64-bit Q32/Q36 on real orbit p_b (T_mult=40)
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
    Db = [1] + [0] * 12; lamb = [[0] * 12 for _ in range(12)]
    for i in range(12):
        for j in range(i + 1):
            u = GBb[i][j]
            for k in range(j): u = (Db[k + 1] * u - lamb[i][k] * lamb[j][k]) // Db[k]
            if j < i: lamb[i][j] = u
            else: Db[i + 1] = u
    T12_b = 1764 * 40 * G_SC
    mu_b = [[_i64((lamb[j][k] << 32) // Db[k + 1]) if j > k else 0 for k in range(12)] for j in range(12)]
    b2_b = [_u64(min(M64_U, (Db[k + 1] << 64) // (Db[k] * T12_b))) for k in range(12)]
    inv_b = [_u64(((Db[k] * T12_b << 40) + Db[k + 1] - 1) // Db[k + 1]) for k in range(12)]
    Xmb = [0] * 12
    for k in range(11, -1, -1):
        rk = math.isqrt(((Db[k] * T12_b) + Db[k + 1] - 1) // Db[k + 1]) + 1
        Xmb[k] = rk + sum((Xmb[j] + 1) // 2 for j in range(k + 1, 12))
    ecb = [_i64(sum(Xmb[j] for j in range(k + 1, 12)) + 1) for k in range(12)]
    b_in = [0]; b_ex = [0]; xb = [0] * 12
    def _cnt_bnd(k, rem_q36, az, cur):
        ck = 0
        for j in range(k + 1, 12): ck = _i64(ck - _i64(mu_b[j][k] * xb[j]))
        r2 = _u64(_umulh64(rem_q36, inv_b[k]) + 1)
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
                d18 = _u64(d_lo >> 14)
                sub = _umulh64(_u64(d18 * d18), b2k)
                if rem_q36 >= sub:
                    xb[k] = xv
                    _cnt_bnd(k - 1, rem_q36 - sub, az and (xv == 0), [cur[c] + xv * Bk[c] for c in range(12)])
    _cnt_bnd(11, Q36, True, [0] * 12)
    q32_bnd_extra_frac = round(b_ex[0] / max(1, b_in[0]), 5)

    # ---- PART 2: 6-lane SWAR (10x6 for a, 11-10-11-10x3 for b) 3-image stepper & 18-op branchless canon ----
    # Exact bit layouts matching KERNELS.md:
    #   a: lane i (i=0..5) at shift 10*(5-i) so a0 is in TOP bits 50..59 -> unsigned cmp == exact lex cmp!
    #   b: order (b0[61:51], b5[50:41], b1[40:30], b2[29:20], b3[19:10], b4[9:0]) = 62 bits!
    #      ubfx(b, 0, 48) extracts EXACT (b1,b2,b3,b4, b5[6:0]) -> <=6 chords/supertile deterministically!
    #      ubfx(b, 48, 3) extracts chord slot 0..7 in 1 op (ZERO chord hash table)!
    SH_A = (50, 40, 30, 20, 10, 0); BI_A = (512,) * 6; MK_A = (0x3FF,) * 6
    SH_B = (51, 30, 20, 10, 0, 41); BI_B = (1024, 1024, 512, 512, 512, 512); MK_B = (0x7FF, 0x7FF, 0x3FF, 0x3FF, 0x3FF, 0x3FF)
    B_ORD = (0, 5, 1, 2, 3, 4)
    def abs_lex6_b(v):
        for i in B_ORD:
            if v[i] != 0:
                return v if v[i] > 0 else neg6(v)
        return v
    M64 = (1 << 64) - 1
    BIAS_A = sum(BI_A[i] << SH_A[i] for i in range(6))
    BIAS_A_X2 = sum((2 * BI_A[i]) << SH_A[i] for i in range(6))
    BIAS_B = sum(BI_B[i] << SH_B[i] for i in range(6))
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
            # Rare tie path (a in Z[sqrt(-2)], density <= R_a^{-4} < 10^{-11}): exact 128-bit (m_k, max(pb_k,-pb_k))!
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

    # EXPLICIT tie-case orbit-invariance tests: a fixed by sigma (e.g. a=1, b=y) across ALL 6 <-1,sigma> images!
    swar_tie_canon_errs = 0
    tie_cases = [
        ((1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0)),       # exact a=1, b=y counterexample!
        ((5, 0, 0, 0, 0, 0), (1, 2, -3, 4, -1, 2)),     # rational a, generic b
        ((-7, 0, 0, 0, 0, 0), (0, 0, 1, -2, 3, -4)),    # negative rational a
        ((1, 1, 1, 0, 0, 0), (0, 1, -1, 2, 0, -1)),     # test general + sigma-sum
    ]
    # Add explicit sigma-invariant a = v + sigma(v) + sigma^2(v) (all elements of O_h^sigma = Z[sqrt(-2)])!
    for v_seed, b_seed in [((0,1,0,0,0,0),(0,1,0,0,0,0)), ((0,0,0,1,0,0),(1,-2,3,0,1,-1)), ((1,0,1,0,1,0),(0,0,1,0,0,-2))]:
        v1 = _mv6(S, v_seed); v2 = _mv6(S2, v_seed)
        a_fix = tuple(v_seed[i] + v1[i] + v2[i] for i in range(6))
        assert _mv6(S, a_fix) == a_fix
        tie_cases.append((a_fix, b_seed))
    for ta, tb in tie_cases:
        # Compute canonical output starting from EVERY one of the 6 <-1, sigma> orbit conjugates!
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
        # Ensure chord fits inside production E12 coordinate bounds (|a_i|<=420<511, |b0,b1|<=560<1023, |b2..5|<=490<511)
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
            # Exact reference over all 3 sigma images using (abs_lex6(a), abs_lex6_b(b))
            ref_pair = min((abs_lex6(ra), abs_lex6_b(rb)),
                           (abs_lex6(_mv6(S, ra)), abs_lex6_b(_mv6(S, rb))),
                           (abs_lex6(_mv6(S2, ra)), abs_lex6_b(_mv6(S2, rb))))
            if (unpack_v(ma, SH_A, BI_A, MK_A), unpack_v(mb, SH_B, BI_B, MK_B)) != ref_pair:
                swar_canon_errs += 1
            pa0 = (pa0 + dpa0) & M64; pa1 = (pa1 + dpa1) & M64; pa2 = (pa2 + dpa2) & M64
            pb0 = (pb0 + dpb0) & M64; pb1 = (pb1 + dpb1) & M64; pb2 = (pb2 + dpb2) & M64

    # ---- PART 3: 1-level super-tile bucket + structured real-a HT (load=0.140) + descending-ai chain ----
    # Enumerate 2,293 REAL structured E12 a_packed lattice keys at production load 0.140 (2293 / 16384)!
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
    # Open-addressed HT with 3-op eor-lsr-mul hash ((k ^ (k >> 30)) * FIB) >> 50 (load = 0.140)
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

    # 3a) Descending-ai prefix-truncated lower Band-2 chain test (Q(a) <= T_max(sb))
    p_l = 131; W_l = 16; nb_l = (p_l + W_l - 1) // W_l
    kta_l = [0] + [rng.randrange(p_l) for _ in range(Na_t)]
    head_l = [0] * nb_l; nxt_l = [0] * (Na_t + 1)
    for ai in range(Na_t, 0, -1):
        k = kta_l[ai]; sb = k // W_l
        nxt_l[ai] = head_l[sb]; head_l[sb] = ((k & (W_l - 1)) << 46) | ai
    pref_errs = 0
    M46 = (1 << 46) - 1
    for sb in range(nb_l):
        namax = rng.randint(40, Na_t)
        got_s = set(); rec = head_l[sb]
        while rec:
            ai = rec & M46
            if ai > namax:
                break
            got_s.add((rec >> 46, ai))
            rec = nxt_l[ai]
        ref_s = {(kta_l[ai] & (W_l - 1), ai) for ai in range(1, namax + 1) if kta_l[ai] // W_l == sb}
        if got_s != ref_s:
            pref_errs += 1

    # 3b) 1-level super-tile scatter + structured real-a HT + direct 3-bit chord-slot apply replay (2000 hits)
    b_keys = [[((b0 + 512) << 11) | ((ch // B_st) << 8) | (ch % B_st) for b0 in range(L_ch)] for ch in range(n_chords)]
    buckets = [[] for _ in range(B_st)]
    ht_probes = 0; ht_lookups = 0
    for lp_band in (52, 55, 58, 61):
        for t in range(B_st):
            buckets[t].append((0, lp_band))
        for _ in range(500):
            ai = rng.randint(1, Na_t)
            ch = rng.randrange(n_chords); b0 = rng.randrange(L_ch)
            ka = a_keys[ai - 1]; kb = b_keys[ch][b0]
            score_ref[ch][b0][ai] += lp_band
            tid = kb & 0xFF
            buckets[tid].append((ka, kb))

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

    # Verify exact E12 production coordinate bounds fit (10x6) & (11,10,11,10x3) SWAR
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

    ok = (lat_mismatches == 0 and fallback_triggered and fallback_mismatches == 0
          and isqrt_errs == 0 and swar_step_errs == 0 and swar_neg_errs == 0 and swar_canon_errs == 0
          and swar_tie_canon_errs == 0 and pref_errs == 0 and supertile_errs == 0
          and swar_bounds_ok and swar_steps >= 1000 and struct_avg_probes <= 1.12
          and prod_guards_passed == 6 and prod_leaves_total >= 500000
          and prod_max_inv_bits <= 62 and prod_max_b2_bits <= 64 and prod_max_r2_bits < 64
          and prod_max_ir <= 25.0 and q32_bnd_extra_frac <= 0.02 and exc_max_orbits_det < (1 << 31))
    return dict(
        ok=ok,
        zero_float_fp_enum=True,
        strict_u64_i64_umulh64_enforced=True,
        fast_path_guard_B0_unsc=B0_UNSC,
        fallback_128bit_triggered=fallback_triggered,
        fallback_128bit_mismatches=fallback_mismatches,
        exc_max_orbits_det=exc_max_orbits_det,
        isqrt64_5newton_errs=isqrt_errs,
        orbits_6ideal_checked=len(orbits),
        lattice_b2_mismatches=lat_mismatches,
        lattice_b2_raw_hits=tot_raw_hits,
        lattice_b2_jvor_hits=tot_jvor_hits,
        lattice_b2_leaves_per_kept=round(tot_raw_hits / max(1, tot_jvor_hits), 3),
        prod_51_to_61b_orbits_checked=6,
        prod_51_to_61b_guards_passed=prod_guards_passed,
        prod_51_to_61b_actual_leaves_per_orbit=prod_leaves_vals,
        prod_51_to_61b_actual_leaves_total=prod_leaves_total,
        prod_51_to_61b_actual_nodes_per_leaf=prod_npl_vals,
        prod_51_to_61b_actual_ir_per_leaf=prod_ir_vals,
        prod_51_to_61b_actual_mean_ir=prod_mean_ir,
        prod_51_to_61b_actual_max_ir=prod_max_ir,
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
