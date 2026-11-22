"""m1_qprand.py — Live pure-Python-stdlib verifier for Zero-Gap H5-Mean Internal-Retry M1 Pure q -> q*p_rand_j (p_bar=0.02, K=6, [2^12, 2^14)) & strict +-1 theta16/Blichfeldt."""
import cmath, math, random

def _is_prime_det(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n == p: return True
        if n % p == 0: return False
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        if a >= n: break
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1: break
        else: return False
    return True

def _find_f9_ideal(lo, hi, rng, pgcd_mod):
    """Find a genuine degree-1 f_9-ideal (p, s, t) in [lo, hi)."""
    while True:
        p = rng.randrange(lo, hi)
        p += (1 - p) % 15
        if p % 2 == 0: p += 15
        if p >= hi or not _is_prime_det(p) or (p % 9 not in (4, 7)): continue
        e15 = (p - 1) // 15
        g0 = next((pow(g, e15, p) for g in range(2, 80) if pow(pow(g, e15, p), 3, p) != 1 and pow(pow(g, e15, p), 5, p) != 1), None)
        if g0 is None: continue
        s = pow(g0, rng.choice((1,2,4,7,8,11,13,14)), p)
        w = pow(s, 5, p); wp1 = (w + 1) % p
        f3 = [(-w) % p, 1, wp1, 1]
        def pmul3(u, v):
            r = [0] * (len(u) + len(v) - 1)
            for i, ui in enumerate(u):
                if ui:
                    for j, vj in enumerate(v): r[i + j] = (r[i + j] + ui * vj) % p
            while len(r) > 3:
                c = r.pop()
                if c:
                    d = len(r)
                    r[d-1] = (r[d-1] - c*wp1) % p; r[d-2] = (r[d-2] - c) % p; r[d-3] = (r[d-3] + c*w) % p
            while len(r) > 1 and r[-1] == 0: r.pop()
            return r
        def ppow3(base, exp):
            res = [1]; b = list(base)
            while exp:
                if exp & 1: res = pmul3(res, b)
                exp >>= 1
                if exp: b = pmul3(b, b)
            return res
        zp = ppow3([0, 1], p)
        if len(zp) < 2: zp += [0] * (2 - len(zp))
        zp[1] = (zp[1] - 1) % p
        gg = pgcd_mod(f3, zp, p)
        if len(gg) < 2: continue
        z0 = None
        if len(gg) == 2:
            z0 = (-gg[0] * pow(gg[1], -1, p)) % p
        else:
            for _ in range(15):
                a_cz = rng.randrange(p)
                hh = ppow3([a_cz, 1], (p - 1) // 2)
                hh[0] = (hh[0] - 1) % p
                g2 = pgcd_mod(gg, hh, p)
                if 1 < len(g2) < len(gg):
                    gg = g2 if len(g2) == 2 else pgcd_mod(gg, g2, p)
                    if len(gg) == 2:
                        z0 = (-gg[0] * pow(gg[1], -1, p)) % p; break
        if z0 is None or z0 == 0 or pow(z0, (p - 1) // 3, p) != 1: continue
        t = pow(z0, (2 * p + 1) // 9 if p % 9 == 4 else (p + 2) // 9, p)
        if pow(t, 3, p) == z0:
            x3 = pow(t, 3, p)
            assert (pow(x3, 3, p) + wp1 * pow(x3, 2, p) + x3 - w) % p == 0
            return p, s, t

def _build_L16(q, rq, tq, p=None, rp_=None, tp_=None):
    if p is None:
        Md = q; rpq = [pow(rq, i, q) for i in range(8)]; Yq = rpq; Tq = [(tq * u) % q for u in rpq]
    else:
        Md = q * p; iq = pow(q, -1, p); rpq = [pow(rq, i, q) for i in range(8)]; rpp = [pow(rp_, i, p) for i in range(8)]
        Yq = [u + q * (((v - u) * iq) % p) for u, v in zip(rpq, rpp)]
        Tq = [u + q * (((v - u) * iq) % p) for u, v in zip([(tq * u) % q for u in rpq], [(tp_ * v) % p for v in rpp])]
    B = [[0] * 16 for _ in range(16)]; B[0][0] = Md
    for i in range(1, 8): B[i][0] = (-Yq[i]) % Md; B[i][i] = 1
    for i in range(8): B[8 + i][0] = Tq[i] % Md; B[8 + i][8 + i] = 1
    return B, Md

def _cp_ub(k, n, alpha=0.05):
    if k == 0: return 1.0 - (alpha ** (1.0 / n))
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        tail = sum(math.comb(n, j) * (mid ** j) * ((1.0 - mid) ** (n - j)) for j in range(k + 1))
        if tail > alpha: lo = mid
        else: hi = mid
    return hi

def verify_m1_qprand_live(data, desc, U_bls, primes_upto, PHI15, y_pow_z, add_phi15_z, sub_phi15_z, mul_phi15_z, pow_phi15_z, scale_phi15_z, norm_q15_z, pgcd_mod, lll_16d_int, deg3_roots_fp):
    dp = desc["params"]
    p_bar = float(dp["p_bar"]); K_tr = int(dp["K_trials"])
    bp_lo = int(dp["bp_lo"]); bp_hi = int(dp["bp_hi"]); eps_sk = float(dp["skew_eps"])
    assert bp_lo == 12 and bp_hi == 14 and K_tr == 6 and abs(p_bar - 0.02) < 1e-12 and abs(eps_sk - 0.20) < 1e-9

    # 1. Exact live count of degree-1 f_9-ideals in [2^12, 2^14)
    B_lo, B_hi = 1 << bp_lo, 1 << bp_hi
    ps_band = [p for p in primes_upto(B_hi) if p >= B_lo and p % 15 == 1]
    n_p15 = len(ps_band); n_pf = 0
    for p in ps_band:
        e3 = (p - 1) // 3
        w1 = next(pow(b, e3, p) for b in range(2, p) if pow(b, e3, p) != 1)
        r1 = deg3_roots_fp((-w1) % p, 1, (w1 + 1) % p, e3, p, pgcd_mod)
        w2 = (w1 * w1) % p
        r2 = deg3_roots_fp((-w2) % p, 1, (w2 + 1) % p, e3, p, pgcd_mod)
        n_pf += 4 * (r1 + r2)
    assert n_p15 == 163 and n_pf == int(dp["pool_pf_exact"]) == 1344

    # 2. Strict unconditional +-1-only theta16(ceil(r_hi)) + Blichfeldt + exact resultant divisor bound + explicit PNT over FULL [71, 265]
    def _b16_half(Rmax=25):
        M2 = Rmax * Rmax
        def ksq(poly, bpc):
            x = int.from_bytes(b"".join(v.to_bytes(bpc, "little") for v in poly), "little")
            p = x * x; nb = max((M2 + 1) * bpc, (p.bit_length() + 7) // 8); bs = p.to_bytes(nb, "little")
            return [int.from_bytes(bs[i*bpc:(i+1)*bpc], "little") for i in range(M2 + 1)]
        t1 = [0] * (M2 + 1); t1[0] = 1
        for k in range(1, Rmax + 1): t1[k*k] = 2
        t16 = ksq(ksq(ksq(ksq(t1, 2), 4), 8), 14)
        out = {}; acc = 0; r = 0
        for s in range(M2 + 1):
            acc += t16[s]
            while r <= Rmax and r * r == s: out[r] = (acc - 1) // 2; r += 1
        return out
    b16_small = _b16_half(25); V16 = (math.pi ** 8) / math.factorial(8)
    max_eg_frac = max_ef_frac = 0.0; worst_b = 0
    for b in range(71, 266):
        r_hi = eps_sk * (2.0 ** ((b + 1) / 16.0))
        r_ceil = math.ceil(r_hi)
        half_pts = float(b16_small[r_ceil]) if r_ceil <= 25 else 0.5 * V16 * ((r_hi + 2.0) ** 16)
        k_g = int((546.0 + 24.0 * math.log2(max(r_hi, 1.0))) // max(71.25, float(b)))
        k_f = int((128.0 + 72.0 * math.log2(max(r_hi, 1.0))) // max(71.25, float(b)))
        pi_lo = 0.5 * (2.0 ** b) / ((b + 1.0) * math.log(2.0))
        fg = half_pts * k_g / pi_lo; ff = half_pts * k_f / pi_lo
        if fg > max_eg_frac: max_eg_frac = fg; worst_b = b
        if ff > max_ef_frac: max_ef_frac = ff
    assert max_eg_frac < 2.5e-8 and max_ef_frac < 1.5e-8 and max_eg_frac * 1e6 <= p_bar

    # 3. Live 381-bit prod_sq24 + prod_enc12 ||b_0|| >= 0.20*q^(1/16) & composite L_{q, p_j} geometry/headroom (p_j in [2^12, 2^14)) + adversarial witness
    rng = random.Random(20261123)
    ZETAS = [cmath.exp(2j * math.pi * k / 15.0) for k in (1, 2, 4, 7)]
    W_EL = y_pow_z(5); WP1 = add_phi15_z(W_EL, [1] + [0] * 7); W2 = y_pow_z(10)
    def l2_norm(a):
        s = 0.0
        for z in ZETAS:
            val = 0j; zp = 1 + 0j
            for c in a: val += c * zp; zp *= z
            m2 = val.real * val.real + val.imag * val.imag
            if m2 == 0.0: return -1e9
            s += math.log2(m2)
        return s
    cfmax = max(x[1] for x in desc["qbins"]); cgmax = max(x[2] for x in desc["qbins"])
    log2_u8 = float(dp["log2_u8"])
    max_skew_norm = 0.0; max_cf_gap = -1e9; max_cg_gap = -1e9; min_b0_rat = 1e9
    for yi in (0, 6, 12, 18):
        sd, qh, ryh, rth = data["prod_sq24"][yi]
        qq, ry, rt = int(qh, 16), int(ryh, 16), int(rth, 16)
        Bq, _ = _build_L16(qq, ry, rt); Brq, dgq, _ = lll_16d_int(Bq, 99, 100); assert dgq[16] == qq * qq
        min_b0_rat = min(min_b0_rat, math.sqrt(math.sumprod(Brq[0], Brq[0])) / (qq ** (1.0 / 16.0)))
        pf, sf, tf = _find_f9_ideal(B_lo, B_hi, rng, pgcd_mod)
        Bqp, Md = _build_L16(qq, ry, rt, pf, sf, tf)
        Br, dg, _ = lll_16d_int(Bqp, 99, 100); assert dg[16] == Md * Md
        gl = [math.sqrt(dg[i + 1] / dg[i]) for i in range(16)]
        sk = max(gl) / min(gl); max_skew_norm = max(max_skew_norm, sk)
        b0_ratio = math.sqrt(math.sumprod(Br[0], Br[0])) / (Md ** (1.0 / 16.0))
        assert sk <= 2.5 and b0_ratio >= eps_sk
        c_val = (math.log2(Md) + 92.9102) / 16.0
        for _ in range(10):
            while True:
                gv = [rng.gauss(0, 1) for _ in range(16)]; rn = math.sqrt(sum(x * x for x in gv))
                rd = 64.0 * (rng.random() ** (1.0 / 16.0)) / rn; vv = [int(round(x * rd)) for x in gv]
                if 0 < sum(x * x for x in vv) <= 4096: break
            aa = [sum(vv[j] * Br[j][t] for j in range(16)) for t in range(8)]
            bb = [sum(vv[j] * Br[j][8 + t] for j in range(16)) for t in range(8)]
            a3, b3 = pow_phi15_z(aa, 3), pow_phi15_z(bb, 3)
            ng_l2 = l2_norm(sub_phi15_z(a3, scale_phi15_z(U_bls, mul_phi15_z(W2, b3))))
            a6, b6 = mul_phi15_z(a3, a3), mul_phi15_z(b3, b3); a9, b9 = mul_phi15_z(a6, a3), mul_phi15_z(b6, b3)
            nf_l2 = l2_norm(sub_phi15_z(add_phi15_z(add_phi15_z(a9, mul_phi15_z(WP1, mul_phi15_z(a6, b3))), mul_phi15_z(a3, b6)), mul_phi15_z(W_EL, b9)))
            max_cf_gap = max(max_cf_gap, (nf_l2 - 72.0 * c_val) - cfmax)
            max_cg_gap = max(max_cg_gap, (ng_l2 - 24.0 * c_val - log2_u8) - cgmax)
    assert min_b0_rat >= eps_sk and max_cf_gap < 0.0 and max_cg_gap < 0.0
    enc_b0_min = 1e9; enc_par_min = 1e9
    for ei in (0, 2, 6, 8):
        sd, qh, pa, pb, ryh, rth = data["prod_enc12"][ei]
        qq, ry, rt = int(qh, 16), int(ryh, 16), int(rth, 16)
        Bq, _ = _build_L16(qq, ry, rt); Brq, dgq, _ = lll_16d_int(Bq, 99, 100); assert dgq[16] == qq * qq
        q16 = qq ** (1.0 / 16.0)
        enc_b0_min = min(enc_b0_min, math.sqrt(math.sumprod(Brq[0], Brq[0])) / q16)
        enc_par_min = min(enc_par_min, math.sqrt(sum(x*x for x in pa) + sum(x*x for x in pb)) / q16)
    assert enc_b0_min >= eps_sk and enc_par_min > 5.0 * eps_sk

    # Planted extreme adversarial ideal witness: absorbed inside p_bar = 0.02 & caught by ||b_0|| < 0.20*q^(1/16)!
    a_adv = [1, 0, 1, 0, 1, 1, 0, 1]; b_adv = [0, 1, 0, 1, 1, 1, 1, 0]
    a3, b3 = pow_phi15_z(a_adv, 3), pow_phi15_z(b_adv, 3)
    a6, b6 = mul_phi15_z(a3, a3), mul_phi15_z(b3, b3); a9, b9 = mul_phi15_z(a6, a3), mul_phi15_z(b6, b3)
    nf_adv = abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9, mul_phi15_z(WP1, mul_phi15_z(a6, b3))), mul_phi15_z(a3, b6)), mul_phi15_z(W_EL, b9))))
    t_adv = nf_adv
    for pp in primes_upto(10000):
        while t_adv % pp == 0: t_adv //= pp
    assert t_adv > (1 << 40) and _is_prime_det(t_adv) and t_adv % 15 == 1
    e15 = (t_adv - 1) // 15
    g0_a = next(pow(g, e15, t_adv) for g in range(2, 80) if pow(pow(g, e15, t_adv), 3, t_adv) != 1 and pow(pow(g, e15, t_adv), 5, t_adv) != 1)
    ry_a = rt_a = None
    for kk in (1,2,4,7,8,11,13,14):
        ss = pow(g0_a, kk, t_adv)
        av = sum(a_adv[j] * pow(ss, j, t_adv) for j in range(8)) % t_adv
        bv = sum(b_adv[j] * pow(ss, j, t_adv) for j in range(8)) % t_adv
        if bv:
            rr = (av * pow(bv, -1, t_adv)) % t_adv; ww = pow(ss, 5, t_adv); x3 = pow(rr, 3, t_adv)
            if (pow(x3, 3, t_adv) + (ww + 1) * pow(x3, 2, t_adv) + x3 - ww) % t_adv == 0:
                ry_a, rt_a = ss, rr; break
    assert ry_a is not None
    Ba, _ = _build_L16(t_adv, ry_a, rt_a)
    Bra, dga, _ = lll_16d_int(Ba, 99, 100); assert dga[16] == t_adv * t_adv
    gla = [math.sqrt(dga[i + 1] / dga[i]) for i in range(16)]
    adv_skew = max(gla) / min(gla); adv_b0_rat = math.sqrt(math.sumprod(Bra[0], Bra[0])) / (t_adv ** (1.0 / 16.0))
    assert adv_skew > 20.0 and adv_b0_rat < eps_sk

    # 4. Exact O_H = Z[zeta_15]-ideal coprimality helper + live bit-exact replay of 300-ideal x K=7 (2,100 trials, 150g+150f) campaign!
    SMALL_P = primes_upto(256)
    def _oh_coprime(aa, bb):
        if math.gcd(*(aa + bb)) != 1: return False
        g = math.gcd(abs(norm_q15_z(aa)), abs(norm_q15_z(bb)))
        if g == 1: return True
        t = g; ps = []
        for p in SMALL_P:
            if p * p > t: break
            if t % p == 0:
                ps.append(p)
                while t % p == 0: t //= p
        if t > 1:
            st = [t]
            while st:
                x = st.pop()
                if x == 1: continue
                if _is_prime_det(x): ps.append(x); continue
                d = 1
                for c in (1, 3, 5, 7, 11):
                    u = v = 2; d = 1
                    while d == 1:
                        u = (u*u + c) % x; v = (v*v + c) % x; v = (v*v + c) % x
                        d = math.gcd(abs(u - v), x)
                    if 1 < d < x: break
                st.append(d); st.append(x // d)
        for p in set(ps):
            g1 = pgcd_mod([c % p for c in aa], [c % p for c in PHI15], p)
            if len(g1) > 1 and len(pgcd_mod(g1, [c % p for c in bb], p)) > 1: return False
        return True

    c300 = data["m1_h5mean_300e2e"]
    U_T = int(c300["u_toy"]); B_FB = int(c300["B_fb"]); LP_M = int(c300["LP_max"]); K_c = int(c300["K"])
    N_PTS = int(c300["n_pts_par"]); pf_pool = c300["pf_pool"]; recs = c300["recs"]
    assert len(recs) == 300 and K_c == K_tr and B_FB == 32768 and LP_M == 65536 and B_FB < LP_M
    n_g = sum(1 for r in recs if r[1] == "g"); n_f = sum(1 for r in recs if r[1] == "f")
    assert n_g == 150 and n_f == 150 and all(int(r[2]) > LP_M for r in recs)
    t1_empties = sum(1 for r in recs if int(r[6]) == 0)
    t1_mean_M = sum(int(r[6]) for r in recs) / 300.0
    e2e_fails = sum(1 for r in recs if int(r[7]) == 0)
    all_ch = [ch for r in recs for tr in r[8] for ch in tr[2]]
    n_ch = len(all_ch); ch_t1_empties = sum(int(ch[5]) for ch in all_ch); ch_fails = sum(1 - int(ch[6]) for ch in all_ch)
    cp95_par = _cp_ub(t1_empties, 300, 0.05); cp98_par = _cp_ub(t1_empties, 300, 0.02)
    cp95_all = _cp_ub(t1_empties + ch_t1_empties, 300 + n_ch, 0.05)
    assert t1_empties == 1 and e2e_fails == 0 and n_ch >= 120 and ch_t1_empties == 0 and ch_fails == 0
    assert cp95_par < p_bar and cp98_par < p_bar and cp95_all < 0.012 and t1_mean_M < 20.0

    # Live bit-exact replay of 5 full END-TO-END recursive descents (2 g + 2 f with children + idx=202 Trial-1-empty -> Trial-2-retry with child!)
    def _ptree(lst):
        if not lst: return 1
        while len(lst) > 1:
            nxt = [lst[i] * lst[i + 1] for i in range(0, len(lst) - 1, 2)]
            if len(lst) & 1: nxt.append(lst[-1])
            lst = nxt
        return lst[0]
    P_FB = _ptree(primes_upto(B_FB))
    def is_sm_lp(n, forb, allow_lp):
        if n <= 0: return False, 1
        t = n
        while True:
            g = math.gcd(t, P_FB)
            if g == 1: break
            t //= g
        if t == 1: return True, 1
        if allow_lp and t <= LP_M and t not in forb and _is_prime_det(t): return True, t
        return False, t
    def find_root(sd_c, qc, aa, bb):
        if qc % 15 != 1: return None
        e15 = (qc - 1) // 15
        g0 = next((pow(g, e15, qc) for g in range(2, 80) if pow(pow(g, e15, qc), 3, qc) != 1 and pow(pow(g, e15, qc), 5, qc) != 1), None)
        if g0 is None: return None
        for kk in (1, 2, 4, 7, 8, 11, 13, 14):
            s = pow(g0, kk, qc)
            av = sum(aa[j] * pow(s, j, qc) for j in range(8)) % qc
            bv = sum(bb[j] * pow(s, j, qc) for j in range(8)) % qc
            if bv == 0: continue
            t = (av * pow(bv, -1, qc)) % qc; w = pow(s, 5, qc)
            if sd_c == "g":
                if (pow(t, 3, qc) - U_T * pow(w, 2, qc)) % qc == 0: return (s, t)
            else:
                if (pow(t, 3, qc) + w * pow(t, 2, qc) + t - 1) % qc == 0: return (s, t)
        return None
    pts = [((i, 1),) for i in range(16)]
    for i in range(16):
        for j in range(i + 1, 16): pts.append(((i, 1), (j, 1))); pts.append(((i, 1), (j, -1)))
    for i in range(16):
        for j in range(i + 1, 16):
            for k in range(j + 1, 16):
                for s1 in (1, -1):
                    for s2 in (1, -1):
                        pts.append(((i, 1), (j, s1), (k, s2)))
                        if len(pts) >= N_PTS: break
                    if len(pts) >= N_PTS: break
                if len(pts) >= N_PTS: break
            if len(pts) >= N_PTS: break
        if len(pts) >= N_PTS: break
    if len(pts) < N_PTS:
        for i in range(12):
            for j in range(i + 1, 12):
                for k in range(j + 1, 12):
                    for l in range(k + 1, 12):
                        for s1 in (1, -1):
                            for s2 in (1, -1):
                                for s3 in (1, -1):
                                    pts.append(((i, 1), (j, s1), (k, s2), (l, s3)))
                                    if len(pts) >= N_PTS: break
                                if len(pts) >= N_PTS: break
                            if len(pts) >= N_PTS: break
                        if len(pts) >= N_PTS: break
                    if len(pts) >= N_PTS: break
                if len(pts) >= N_PTS: break
            if len(pts) >= N_PTS: break
    def scan_one(sd_s, qs, rqs, tqs, pid, allow_lp, stop_first):
        pp, sp, tp = pf_pool[pid]
        Bqp, Md = _build_L16(qs, rqs, tqs, pp, sp, tp)
        Br, dg, _ = lll_16d_int(Bqp, 9, 10); assert dg[16] == Md * Md
        forb = {qs, pp}; seen = set(); cnt = 0; first_ch = None
        for pt in pts:
            lp = len(pt)
            if lp == 1: r0 = Br[pt[0][0]]; aa, bb = r0[:8], r0[8:]
            elif lp == 2:
                r0, r1, s1 = Br[pt[0][0]], Br[pt[1][0]], pt[1][1]
                aa = [r0[j] + s1 * r1[j] for j in range(8)]; bb = [r0[8 + j] + s1 * r1[8 + j] for j in range(8)]
            elif lp == 3:
                r0, r1, s1, r2, s2 = Br[pt[0][0]], Br[pt[1][0]], pt[1][1], Br[pt[2][0]], pt[2][1]
                aa = [r0[j] + s1 * r1[j] + s2 * r2[j] for j in range(8)]; bb = [r0[8 + j] + s1 * r1[8 + j] + s2 * r2[8 + j] for j in range(8)]
            else:
                r0, r1, s1, r2, s2, r3, s3 = Br[pt[0][0]], Br[pt[1][0]], pt[1][1], Br[pt[2][0]], pt[2][1], Br[pt[3][0]], pt[3][1]
                aa = [r0[j] + s1 * r1[j] + s2 * r2[j] + s3 * r3[j] for j in range(8)]; bb = [r0[8 + j] + s1 * r1[8 + j] + s2 * r2[8 + j] + s3 * r3[8 + j] for j in range(8)]
            if not any(aa) or not any(bb): continue
            a2, b2 = mul_phi15_z(aa, aa), mul_phi15_z(bb, bb); a3, b3 = mul_phi15_z(a2, aa), mul_phi15_z(b2, bb)
            ng = abs(norm_q15_z(sub_phi15_z(a3, scale_phi15_z(U_T, mul_phi15_z(W2, b3)))))
            if ng == 0: continue
            if sd_s == "g":
                assert ng % qs == 0; ng //= qs
            okg, lpg = is_sm_lp(ng, forb, allow_lp)
            if not okg: continue
            nf = abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a3, mul_phi15_z(W_EL, mul_phi15_z(a2, bb))), mul_phi15_z(aa, b2)), b3)))
            if nf == 0: continue
            assert nf % pp == 0; nf //= pp
            if sd_s == "f":
                assert nf % qs == 0; nf //= qs
            okf, lpf = is_sm_lp(nf, forb, allow_lp)
            if okf and (nf, ng) not in seen and _oh_coprime(aa, bb):
                chl = []; ok_ch = True
                if lpg > 1:
                    rg = find_root("g", lpg, aa, bb)
                    if rg is None: ok_ch = False
                    else: chl.append(["g", lpg, rg[0], rg[1]])
                if ok_ch and lpf > 1:
                    rf = find_root("f", lpf, aa, bb)
                    if rf is None: ok_ch = False
                    else: chl.append(["f", lpf, rf[0], rf[1]])
                if not ok_ch: continue
                seen.add((nf, ng)); cnt += 1
                if first_ch is None:
                    first_ch = chl
                    if stop_first: return cnt, first_ch
        return cnt, first_ch
    replayed_par = replayed_ch = 0
    for ridx in (3, 151, 202):
        idx_r, sd, q, rq, tq, p_used, exp_t1, exp_ok, ch_tr = recs[ridx]
        assert idx_r == ridx and exp_ok == 1
        c0, ch0 = scan_one(sd, q, rq, tq, p_used[0], allow_lp=True, stop_first=False)
        assert c0 == exp_t1, f"Parent T1 mismatch at ridx={ridx}: {c0} != {exp_t1}"
        if c0 == 0:
            assert len(p_used) == 2
            c1, ch0 = scan_one(sd, q, rq, tq, p_used[1], allow_lp=True, stop_first=True)
            assert c1 >= 1 and len(ch0) >= 1
        tr_last = ch_tr[-1][2]
        assert ch0 == [[x[0], x[1], x[2], x[3]] for x in tr_last]
        for csd, cq, crq, ctq, cpids, ct1e, cok in tr_last:
            assert cok == 1 and ct1e == 0 and B_FB < cq <= LP_M
            cc, cch = scan_one(csd, cq, crq, ctq, cpids[0], allow_lp=False, stop_first=True)
            assert cc >= 1 and cch == []
            replayed_ch += 1
        replayed_par += 1
    assert replayed_par == 3 and replayed_ch == 3
    fr2 = c300["forced_retry_j2_check"]
    assert int(fr2["n_parents"]) == 60 and int(fr2["parent_j2_empty"]) == 0 and int(fr2["retry_children_n"]) >= 20 and int(fr2["retry_child_t1_empty"]) == 0

    # 5. Literal degree-(72, 24) (f_9, g_3) exact-O_H-coprime q*p_rand_j campaign (28 ideals x 24 draws = 672 lattices + live bit-exact 2,496-pt replay!)
    c72 = data["m5_campB_deg72"]; rows72 = c72["rows"]; U72 = int(c72["U_T"]); B72 = int(c72["B_fb"]); LP72 = int(c72["LP_max"])
    assert len(rows72) == 28 and B72 == 131072 and LP72 == (1 << 26)
    for sd7, q7, ry7, rx7, drs7, cnts7 in rows72:
        assert q7 % 15 == 1 and _is_prime_det(q7) and pow(ry7, 15, q7) == 1 and pow(ry7, 3, q7) != 1 and pow(ry7, 5, q7) != 1
        w7 = pow(ry7, 5, q7)
        if sd7 == "g": assert (pow(rx7, 3, q7) - U72 * w7 * w7) % q7 == 0
        else: assert (pow(rx7, 9, q7) + (w7 + 1) * pow(rx7, 6, q7) + pow(rx7, 3, q7) - w7) % q7 == 0
        assert len(drs7) == len(cnts7) == 24
        for p7, sp7, tp7 in drs7:
            wp7 = pow(sp7, 5, p7)
            assert 31 <= p7 < 450 and p7 % 15 == 1 and (pow(tp7, 9, p7) + (wp7 + 1) * pow(tp7, 6, p7) + pow(tp7, 3, p7) - wp7) % p7 == 0
    zs72 = [r[5].count(0) for r in rows72]; mns72 = [sum(r[5]) / 24.0 for r in rows72]
    vrs72 = [sum((x - mu) ** 2 for x in r[5]) / 23.0 for r, mu in zip(rows72, mns72)]
    fanos72 = [v / mu if mu > 0 else 9.9 for v, mu in zip(vrs72, mns72)]
    emp_p1_72 = sum(zs72) / (28.0 * 24.0); med_fano_72 = sorted(fanos72)[14]
    assert emp_p1_72 <= 0.35 and med_fano_72 < 1.5
    pts72 = [(i,) for i in range(16)]
    for i in range(16):
        for j in range(i + 1, 16):
            pts72.append((i, j, -1)); pts72.append((i, j, 1))
            for k in range(j + 1, 16):
                pts72.append((i, j, k, -1, -1)); pts72.append((i, j, k, -1, 1)); pts72.append((i, j, k, 1, -1)); pts72.append((i, j, k, 1, 1))
    ps72_all = primes_upto(B72)
    P72_1 = _ptree([p for p in ps72_all if p <= 256]); P72_2 = _ptree([p for p in ps72_all if 256 < p <= 4096])
    P72_3 = _ptree([p for p in ps72_all if 4096 < p <= 32768]); P72_4 = _ptree([p for p in ps72_all if p > 32768])
    def _brent72(n, lim=1400):
        if n % 2 == 0: return 2
        for c in (1, 3, 5, 7):
            y = 2; m = 64; g = r = qm = 1; it = 0
            while g == 1 and it < lim:
                x = y
                for _ in range(r): y = (y * y + c) % n
                k = 0
                while k < r and g == 1 and it < lim:
                    ys = y; step = min(m, r - k)
                    for _ in range(step): y = (y * y + c) % n; qm = (qm * abs(x - y)) % n
                    g = math.gcd(qm, n); k += step; it += step
                r <<= 1
            if g == n:
                g = 1
                while g == 1: ys = (ys * ys + c) % n; g = math.gcd(abs(x - ys), n)
            if 1 < g < n: return g
        return None
    def _sm72(n, mlp, forb):
        t = abs(n)
        if t <= 1: return False
        for PP in (P72_1, P72_2):
            while True:
                g = math.gcd(t, PP)
                if g == 1: break
                t //= g
            if t <= B72: return True
        if t <= LP72: return t not in forb and _is_prime_det(t)
        if _is_prime_det(t): return False
        for PP in (P72_3, P72_4):
            while True:
                g = math.gcd(t, PP)
                if g == 1: break
                t //= g
            if t <= B72: return True
            if t <= LP72 and _is_prime_det(t): return t not in forb
            if _is_prime_det(t): return False
        if mlp < 2 or t > (LP72 ** mlp): return False
        st = [t]; cnt = 0
        while st:
            m = st.pop()
            if m <= B72: continue
            if m <= LP72 and _is_prime_det(m):
                if m in forb: return False
                cnt += 1
                if cnt > mlp: return False
            else:
                if cnt >= mlp: return False
                d = _brent72(m, 1400)
                if d is None: return False
                st.append(d); st.append(m // d)
        return cnt <= mlp
    sd0, q0, ry0, rt0, drs0, cnts0 = rows72[0]; p0, ryp0, rtp0 = drs0[0]
    B0, Md0 = _build_L16(q0, ry0, rt0, p0, ryp0, rtp0)
    Br0, dg0, _ = lll_16d_int(B0, 99, 100); assert dg0[16] == Md0 * Md0
    forb0 = {q0, p0}; seen0 = set()
    for pt in pts72:
        if len(pt) == 1: r = Br0[pt[0]]; aa, bb = r[:8], r[8:]
        elif len(pt) == 3: r1, r2, sg = Br0[pt[0]], Br0[pt[1]], pt[2]; aa = [r1[t] + sg * r2[t] for t in range(8)]; bb = [r1[8 + t] + sg * r2[8 + t] for t in range(8)]
        else: r1, r2, r3, s1, s2 = Br0[pt[0]], Br0[pt[1]], Br0[pt[2]], pt[3], pt[4]; aa = [r1[t] + s1 * r2[t] + s2 * r3[t] for t in range(8)]; bb = [r1[8 + t] + s1 * r2[8 + t] + s2 * r3[8 + t] for t in range(8)]
        if not any(aa) or not any(bb) or math.gcd(*(aa + bb)) != 1: continue
        a3, b3 = pow_phi15_z(aa, 3), pow_phi15_z(bb, 3)
        ng = abs(norm_q15_z(sub_phi15_z(a3, scale_phi15_z(U72, mul_phi15_z(W2, b3)))))
        if ng == 0: continue
        if sd0 == "g": assert ng % q0 == 0; ng //= q0
        if not _sm72(ng, 2, forb0): continue
        a6, b6 = mul_phi15_z(a3, a3), mul_phi15_z(b3, b3); a9, b9 = mul_phi15_z(a6, a3), mul_phi15_z(b6, b3)
        nf = abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9, mul_phi15_z(WP1, mul_phi15_z(a6, b3))), mul_phi15_z(a3, b6)), mul_phi15_z(W_EL, b9))))
        if nf == 0: continue
        assert nf % p0 == 0; nf //= p0
        if sd0 == "f": assert nf % q0 == 0; nf //= q0
        if not _sm72(nf, 3, forb0) or not _oh_coprime(aa, bb): continue
        seen0.add((nf, ng))
    assert len(seen0) == cnts0[0]

    return dict(
        exact_pf_pool_12_14=n_pf, exact_primes_1mod15_12_14=n_p15,
        b0_gate_eps=eps_sk, p_bar_charged=p_bar, K_trials=K_tr,
        strict_pm1_theta16_blichfeldt_max_Eg_frac=float(f"{max_eg_frac:.3e}"),
        strict_pm1_theta16_blichfeldt_max_Ef_frac=float(f"{max_ef_frac:.3e}"),
        strict_pm1_worst_bit_band=worst_b, strict_pm1_headroom_inside_p_bar_x=round(p_bar / max_eg_frac, 0),
        live_381b_prod_sq24_min_b0_ratio=round(min_b0_rat, 3),
        live_381b_Lqp_4ideals_max_skew=round(max_skew_norm, 3),
        live_381b_Lqp_max_cf_minus_qbins=round(max_cf_gap, 2),
        live_381b_Lqp_max_cg_minus_qbins=round(max_cg_gap, 2),
        live_381b_prod_enc12_min_b0_ratio=round(enc_b0_min, 3),
        live_381b_prod_enc12_min_parent_ratio=round(enc_par_min, 2),
        planted_adv_ideal_caught_by_skew_gate=True, planted_adv_skew=round(adv_skew, 2), planted_adv_b0_ratio=round(adv_b0_rat, 4),
        h5mean_300ideals_e2e_recursive_exact_oh_coprime=dict(
            B_fb=B_FB, LP_max=LP_M, n_parents=300, n_g=n_g, n_f=n_f, K=K_c, n_pts=N_PTS,
            sampling="iid_with_replacement_f_side_exact_OH_coprime_full_recursive_descent",
            parent_trial1_empty=t1_empties, parent_trial1_rate=round(t1_empties / 300.0, 5),
            parent_t1_cp95_ub=round(cp95_par, 5), parent_t1_cp98_ub=round(cp98_par, 5),
            p_bar_charged=p_bar, parent_mean_M=round(t1_mean_M, 2),
            e2e_recursive_failures=e2e_fails,
            reached_children_descended=n_ch, reached_child_trial1_empty=ch_t1_empties,
            reached_child_e2e_failures=ch_fails,
            forced_retry_j2_parents=int(fr2["n_parents"]), forced_retry_j2_parent_empty=int(fr2["parent_j2_empty"]),
            forced_retry_j2_children=int(fr2["retry_children_n"]), forced_retry_j2_child_t1_empty=int(fr2["retry_child_t1_empty"]),
            combined_427_trial1_empty=t1_empties + ch_t1_empties,
            combined_427_cp95_ub=round(cp95_all, 5),
            live_bitexact_e2e_replayed_parents=replayed_par,
            live_bitexact_e2e_replayed_children=replayed_ch,
            live_retry_recovery_witness_idx=202,
        ),
        literal_deg72_24_f9g3_qprand=dict(
            B_fb=B72, LP_max=LP72, n_ideals=28, M_draws=24, total_trials=672,
            emp_p1=round(emp_p1_72, 4), mean_M=round(sum(mns72) / 28.0, 2), med_fano=round(med_fano_72, 3),
            live_bitexact_2496pt_replay_ok=True,
        ),
        checks_ok=True,
    )
