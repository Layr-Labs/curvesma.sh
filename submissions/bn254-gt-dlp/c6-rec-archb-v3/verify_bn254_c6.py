"""
verify_bn254_c6.py — Self-contained pure-stdlib live prototypes for BN254 v3:
  1. 3-regime Band-3 7D LLL + weighted-cylinder Fincke-Pohst coset-walk prototype spanning h in [6, 105]
     (0 missed hits + live empirical leaves/hit <= 2.15 & nodes/hit <= exact Gaussian-volume bound);
  2. Exact 12D C6 special-q (6 bases spanning 2^66..2^213) & 2 boot lattice audit: exact Bareiss dets,
     exact 6x6 integer resultant norms, full 6/6 QBINS + 6/6 QBINS_LO + 24/24 QBINS_SUB + R=16..128 +
     3/3 BOOT_BINS + 10/10 BOOT_SUB B_12(4) stochastic dominance + live lattice CV^2 Gamma-Poisson r_eff;
  3. Pure-Python production hexagonal Z[omega] filtering simulator (exact Dickman-kappa occurrence profile,
     singleton + union-find clique purge, Markowitz merge W=200, Boudot recorded pivots, TABM dense-band
     grouping) measuring ALL FOUR cap metrics (N_sw/U', N_l/U', k5_exp, tab_r) live at U'=2^13 & 2^14
     plus exact Z[omega] Markowitz-merge kernel preservation check.
"""
import bisect, heapq, math, random
from math import ceil, exp, floor, log2, pi, sqrt
import model_bn254_c6 as M
import descent_bn254_c6 as DESC

def verify_7d_walk_prototype():
    """Live 3-regime 7D Fincke-Pohst coset walk bracketing the full Band-3 density range h in [6, 105]!"""
    T4 = M.load_data()["t4_geom"]; GB = T4["GB"]
    def eval_poly_mod(cs, s, p):
        v = 0
        for c in reversed(cs): v = (v * s + c) % p
        return v
    def eval_zk_mod(vec, s, p):
        inv3 = pow(3, -1, p); sp = [1]
        for _ in range(5): sp.append((sp[-1] * s) % p)
        return (sum(vec[k] * (sum(M.ZK_NUM[k][j] * sp[j] for j in range(6)) % p) for k in range(6) if vec[k]) * inv3) % p
    def qf6(v): return sum(GB[i][j] * v[i] * v[j] for i in range(6) for j in range(6))
    def is_prime(n): return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))

    L6 = [[0.0] * 6 for _ in range(6)]
    for r in range(6):
        for c in range(r + 1):
            sv = GB[r][c] - sum(L6[r][k] * L6[c][k] for k in range(c))
            L6[r][c] = sqrt(sv) if r == c else sv / L6[c][c]
    RB2, wt_t = 85.0, M.FP_WT_T
    R7_sq = 1.0 + wt_t * wt_t
    b_all, bv = [], [0] * 6
    def enum_b(idx, rem):
        sh = sum(L6[j][idx] * bv[j] for j in range(idx + 1, 6))
        diag = L6[idx][idx]; rad = sqrt(max(0.0, rem))
        for x in range(ceil((-rad - sh) / diag - 1e-12), floor((rad - sh) / diag + 1e-12) + 1):
            d = diag * x + sh; nr = rem - d * d
            if nr >= -1e-9:
                bv[idx] = x
                if idx == 0: b_all.append(tuple(bv))
                else: enum_b(idx - 1, max(0.0, nr))
    enum_b(5, RB2 + 1e-9)

    def run_regime(name, p_lo, p_hi, n_id, L_a, cosets):
        w_ideals = []
        for p in range(p_lo, p_hi, 2):
            if not is_prime(p): continue
            for s in range(p):
                if eval_poly_mod(M.H_COEFFS, s, p) == 0:
                    uw = ((M.U_BN254 % p) * eval_poly_mod(M.W_COEFFS, s, p)) % p
                    if uw == 0: continue
                    inv_uw = pow(uw, -1, p)
                    if pow(inv_uw, (p - 1) // 2, p) == 1:
                        rho_r = next(r for r in range(1, (p + 1) // 2) if (r * r) % p == inv_uw)
                        w_ideals.append((p, s, rho_r)); break
            if len(w_ideals) >= n_id: break
        hits = missed = 0; nodes = [0]; leaves = [0]
        sc_t = (2.0 / L_a) * wt_t; inv_rb = 1.0 / sqrt(RB2)
        def to_metric(v):
            return [sum(L6[j][i] * v[j] for j in range(i, 6)) * inv_rb for i in range(6)] + [v[6] * sc_t]
        for p, s, rho_r in w_ideals:
            inv_rho = pow(rho_r, -1, p)
            zk_s = [eval_zk_mod([1 if j == k else 0 for j in range(6)], s, p) for k in range(6)]
            B7 = [[0] * 7 for _ in range(7)]; B7[0][0] = p
            for k in range(1, 6): B7[k][k] = 1; B7[k][0] = (-zk_s[k]) % p
            B7[6][6] = 1; B7[6][0] = inv_rho % p
            Bs = [[round(1000 * x) for x in to_metric(row)] for row in B7]
            Bi = [row[:] for row in B7]
            def gs7():
                Bo, mu, Bn = [], [[0.0] * 7 for _ in range(7)], [0.0] * 7
                for i in range(7):
                    v = [float(x) for x in Bs[i]]
                    for j in range(i):
                        mu[i][j] = sum(Bs[i][k] * Bo[j][k] for k in range(7)) / Bn[j]
                        v = [v[k] - mu[i][j] * Bo[j][k] for k in range(7)]
                    Bo.append(v); Bn[i] = sum(x * x for x in v)
                return Bo, mu, Bn
            Bo, mu, Bn = gs7(); k = 1
            while k < 7:
                for j in range(k - 1, -1, -1):
                    q = round(mu[k][j])
                    if q:
                        Bi[k] = [Bi[k][i] - q * Bi[j][i] for i in range(7)]
                        Bs[k] = [Bs[k][i] - q * Bs[j][i] for i in range(7)]
                        Bo, mu, Bn = gs7()
                if Bn[k] >= (0.99 - mu[k][k - 1] ** 2) * Bn[k - 1]: k += 1
                else:
                    Bi[k], Bi[k - 1] = Bi[k - 1], Bi[k]; Bs[k], Bs[k - 1] = Bs[k - 1], Bs[k]; Bo, mu, Bn = gs7(); k = max(k - 1, 1)
            Mmat = [to_metric(row) for row in Bi]
            Bo, mu, Bn = [], [[0.0] * 7 for _ in range(7)], [0.0] * 7
            for i in range(7):
                v = Mmat[i][:]
                for j in range(i):
                    mu[i][j] = sum(Mmat[i][k] * Bo[j][k] for k in range(7)) / Bn[j]
                    v = [v[k] - mu[i][j] * Bo[j][k] for k in range(7)]
                Bo.append(v); Bn[i] = sum(x * x for x in v)
            for a0_mod in cosets:
                rhs_mod = (a0_mod * inv_rho) % p
                td = to_metric([-rhs_mod, 0, 0, 0, 0, 0, 0.5 * L_a])
                c_gs = [sum(td[k] * Bo[j][k] for k in range(7)) / Bn[j] for j in range(7)]
                found, x = set(), [0] * 7
                def fp7(idx, rem):
                    nodes[0] += 1
                    sh = c_gs[idx] - sum(mu[j][idx] * x[j] for j in range(idx + 1, 7))
                    rad = sqrt(max(0.0, rem / Bn[idx]))
                    for xi in range(ceil(sh - rad - 1e-9), floor(sh + rad + 1e-9) + 1):
                        d = xi - sh; nr = rem - d * d * Bn[idx]
                        if nr >= -1e-9:
                            x[idx] = xi
                            if idx == 0:
                                leaves[0] += 1
                                w = [rhs_mod + sum(x[i] * Bi[i][0] for i in range(7))] + [sum(x[i] * Bi[i][k] for i in range(7)) for k in range(1, 7)]
                                if 0 <= w[6] < L_a and qf6(w[:6]) <= RB2 + 1e-9:
                                    found.add(tuple(w))
                            else:
                                fp7(idx - 1, max(0.0, nr))
                fp7(6, R7_sq + 1e-7)
                gt = set()
                for b in b_all:
                    bs = sum(b[k] * zk_s[k] for k in range(6)) % p
                    t = ((bs - rhs_mod) * rho_r) % p
                    if t < L_a: gt.add(b + (t,))
                hits += len(gt); missed += len(gt ^ found)
        h_emp = hits / max(1, len(w_ideals) * len(cosets))
        lph = leaves[0] / max(1, hits); nph = nodes[0] / max(1, hits)
        g_bound = M.fp7_gauss_nph(h_emp)
        return dict(
            name=name, n_ideals=len(w_ideals), cosets=len(cosets), hits=hits, missed=missed,
            h_emp=round(h_emp, 2), leaves_per_hit=round(lph, 4), nodes_per_hit=round(nph, 4),
            gauss_nodes_bound=round(g_bound, 4), ratio_to_gauss=round(nph / g_bound, 3),
            ok=(missed == 0 and hits > 0 and lph <= M.FP_LEAVES_PER_HIT and nph <= g_bound),
        )

    regimes = [
        run_regime("sparse_band3_hi", 10007, 18000, 6, 32, (0, 17, 123, 311, 509)),
        run_regime("mid_band3",        2503,  5000, 6, 32, (0, 17, 123)),
        run_regime("dense_band3_lo",   1009,  3000, 6, 68, (0, 17, 123)),
    ]
    return dict(
        ok=all(r["ok"] for r in regimes),
        vr_exact=round(M.fp7_vr(), 6),
        bound_leaves_per_hit=M.FP_LEAVES_PER_HIT,
        total_hits=sum(r["hits"] for r in regimes),
        total_missed=sum(r["missed"] for r in regimes),
        regimes=regimes,
    )

def verify_12d_descent_lattices():
    """Live audit of 6 special-q (2^66..2^213) & 2 boot 12D lattices: exact dets, 6x6 integer norms,
    full QBINS/QBINS_LO/QBINS_SUB + R=16..128 + BOOT_BINS/BOOT_SUB B_12(4) stochastic dominance & CV^2!"""
    D = M.load_data(); T4 = D["t4_geom"]; B12 = D["descent"]["bases_12d"]
    EMB_ZK = [[complex(r, im) for r, im in row] for row in T4["EMB_ZK"]]
    EMB_W = [complex(r, im) for r, im in T4["EMB_W"]]
    EMB_F = {d: [sum(cs[k] * (EMB_ZK[j][1] ** k) for k in range(6)) for j in range(3)] for d, cs in M.F_REC_COEFFS.items()}
    LOG2U6, LOG2P6 = 6.0 * log2(M.U_BN254), 6.0 * log2(M.P_BN254)

    def det_int(Mat):
        n = len(Mat); A = [row[:] for row in Mat]; sign = 1; prev = 1
        for k in range(n - 1):
            if A[k][k] == 0:
                piv = next((r for r in range(k + 1, n) if A[r][k] != 0), -1)
                if piv < 0: return 0
                A[k], A[piv] = A[piv], A[k]; sign = -sign
            piv = A[k][k]
            for i in range(k + 1, n):
                for j in range(k + 1, n):
                    A[i][j] = (A[i][j] * piv - A[i][k] * A[k][j]) // prev
            prev = piv
        return sign * A[n - 1][n - 1]
    def pmul_h_int(a, b):
        res = [0] * 11
        for i in range(6):
            if a[i]:
                for j in range(6): res[i + j] += a[i] * b[j]
        for k in range(10, 5, -1):
            c = res[k]
            if c:
                res[k] = 0
                for j in range(6): res[k - 6 + j] -= c * M.H_COEFFS[j]
        return res[:6]
    def exact_norm_q(poly_num, denom):
        cols, cur = [], poly_num[:]
        for j in range(6):
            cols.append(cur)
            if j < 5:
                nxt = [0] + cur[:5]; c = cur[5]
                if c:
                    for k in range(6): nxt[k] -= c * M.H_COEFFS[k]
                cur = nxt
        d = abs(det_int([[cols[c][r] for c in range(6)] for r in range(6)]))
        assert d % (denom ** 6) == 0
        return d // (denom ** 6)
    def exact_log2_norms(a, b):
        an = [sum(a[k] * M.ZK_NUM[k][j] for k in range(6)) for j in range(6)]
        bn = [sum(b[k] * M.ZK_NUM[k][j] for k in range(6)) for j in range(6)]
        an2, bn2 = pmul_h_int(an, an), pmul_h_int(bn, bn)
        wan2 = pmul_h_int(M.W_COEFFS, an2)
        ng = exact_norm_q([M.U_BN254 * wan2[j] - bn2[j] for j in range(6)], 9)
        a2p = [[1,0,0,0,0,0], an2]; b2p = [[1,0,0,0,0,0], bn2]
        for _ in range(3):
            a2p.append(pmul_h_int(a2p[-1], an2)); b2p.append(pmul_h_int(b2p[-1], bn2))
        fnum = [0] * 6
        for deg_a, cs in M.F_REC_COEFFS.items():
            k = deg_a // 2
            term = pmul_h_int(pmul_h_int(cs, a2p[k]), b2p[4 - k])
            for j in range(6): fnum[j] += term[j]
        nf = exact_norm_q(fnum, 3 ** 8)
        return log2(nf), log2(ng)
    def eval_zk_mod(vec, s, mod):
        inv3 = pow(3, -1, mod); sp = [1]
        for _ in range(5): sp.append((sp[-1] * s) % mod)
        return (sum(vec[k] * (sum(M.ZK_NUM[k][j] * sp[j] for j in range(6)) % mod) for k in range(6) if vec[k]) * inv3) % mod
    def zk_to_poly_mod(vec, mod):
        inv3 = pow(3, -1, mod)
        return [(sum(vec[k] * M.ZK_NUM[k][j] for k in range(6)) * inv3) % mod for j in range(6)]
    def log2_norms_fp(a, b, need_f=True):
        lf = lg = 0.0
        mx = max(max(abs(x) for x in a), max(abs(x) for x in b), 1)
        sh = max(0, mx.bit_length() - 40); sc = float(1 << sh)
        af = [x / sc for x in a]; bf = [x / sc for x in b]
        for j in range(3):
            za = sum(EMB_ZK[j][k] * af[k] for k in range(6))
            zb = sum(EMB_ZK[j][k] * bf[k] for k in range(6))
            gv = M.U_BN254 * EMB_W[j] * za * za - zb * zb
            lg += log2(gv.real**2 + gv.imag**2) + 4.0 * sh
            if need_f:
                za2, zb2 = za * za, zb * zb
                fv = (EMB_F[8][j]*(za2**4) + EMB_F[6][j]*(za2**3)*zb2 + EMB_F[4][j]*(za2**2)*(zb2**2)
                      + EMB_F[2][j]*za2*(zb2**3) + EMB_F[0][j]*(zb2**4))
                lf += log2(fv.real**2 + fv.imag**2) + 16.0 * sh
        return lf, lg

    S_R4 = log2(D["descent"]["cells"]["4"])
    rng = random.Random(42)
    sq_ok = True; max_exact_err = 0.0
    cfs_r4, cgs_r4, cfs_hiR, cgs_hiR, lat_yields = [], [], [], [], []
    q_bits_list = []
    for item in B12["special_q"]:
        q, s, rho_r = int(item["q"]), int(item["s"]), int(item["rho"])
        q_bits_list.append(q.bit_length())
        B = item["B"]
        if abs(det_int(B)) != q: sq_ok = False
        for col in range(12):
            av = [B[r][col] for r in range(6)]; bv = [B[6 + r][col] for r in range(6)]
            if (eval_zk_mod(av, s, q) - rho_r * eval_zk_mod(bv, s, q)) % q != 0: sq_ok = False
        # Exact 6x6 big-int resultant norm vs complex embedding cross-check
        v0 = [B[r][0] + B[r][1] for r in range(12)]
        lf_e, lg_e = exact_log2_norms(v0[:6], v0[6:])
        lf_f, lg_f = log2_norms_fp(v0[:6], v0[6:], True)
        max_exact_err = max(max_exact_err, abs(lf_e - lf_f), abs(lg_e - lg_f))
        c_r4 = (log2(q) + S_R4) / 12.0
        ys_lat = []
        for _ in range(80):
            while True:
                x = [rng.randint(-2, 2) for _ in range(12)]
                if any(x) and sum(z * z for z in x) <= 16: break
            v = [sum(B[r][c] * x[c] for c in range(12)) for r in range(12)]
            lf, lg = log2_norms_fp(v[:6], v[6:], True)
            cf_v = lf - 48.0 * c_r4; cg_v = lg - 12.0 * c_r4 - LOG2U6
            cfs_r4.append(cf_v); cgs_r4.append(cg_v)
            nf_ref = cf_v + 48.0 * ((100.0 + S_R4) / 12.0)
            ng_ref = cg_v + 12.0 * ((100.0 + S_R4) / 12.0) + LOG2U6
            ys_lat.append(M.rho(nf_ref / 59.0) * M.rho(ng_ref / 59.0))
        lat_yields.append(sum(ys_lat) / len(ys_lat))
        # Large-radius R in {16, 64, 128} scale-invariance check
        for R_big in (16, 64, 128):
            S_Rb = log2((pi ** 6 / 720.0) * (R_big ** 12))
            c_rb = (log2(q) + S_Rb) / 12.0
            for _ in range(15):
                gv = [rng.gauss(0.0, 1.0) for _ in range(12)]
                gn = sqrt(sum(z * z for z in gv))
                rad = R_big * (rng.random() ** (1.0 / 12.0))
                x = [int(round(z * rad / gn)) for z in gv]
                if not any(x): x[0] = 1
                v = [sum(B[r][c] * x[c] for c in range(12)) for r in range(12)]
                lf, lg = log2_norms_fp(v[:6], v[6:], True)
                cfs_hiR.append(lf - 48.0 * c_rb); cgs_hiR.append(lg - 12.0 * c_rb - LOG2U6)

    # Full stochastic dominance checks on special-q (ALL 6/6 QBINS, 6/6 QBINS_LO, 24/24 QBINS_SUB, 6/6 R=16..128!)
    cfs_r4.sort(); cgs_r4.sort(); cfs_hiR.sort(); cgs_hiR.sort()
    qbins_ok = all(
        cfs_r4[min(len(cfs_r4)-1, ceil(((k+1)/6.0)*len(cfs_r4))-1)] <= DESC.QBINS[k][1] and
        cgs_r4[min(len(cgs_r4)-1, ceil(((k+1)/6.0)*len(cgs_r4))-1)] <= DESC.QBINS[k][2]
        for k in range(6)
    )
    qbins_lo_ok = all(
        cfs_r4[min(len(cfs_r4)-1, floor((k/6.0)*len(cfs_r4)))] >= DESC.QBINS_LO[k][0] and
        cgs_r4[min(len(cgs_r4)-1, floor((k/6.0)*len(cgs_r4)))] >= DESC.QBINS_LO[k][1]
        for k in range(6)
    )
    qbins_sub_ok = all(
        cfs_r4[min(len(cfs_r4)-1, ceil(((k+1)/24.0)*len(cfs_r4))-1)] <= DESC.QBINS_SUB[k][0] and
        cgs_r4[min(len(cgs_r4)-1, ceil(((k+1)/24.0)*len(cgs_r4))-1)] <= DESC.QBINS_SUB[k][1]
        for k in range(24)
    )
    hiR_ok = all(
        cfs_hiR[min(len(cfs_hiR)-1, ceil(((k+1)/6.0)*len(cfs_hiR))-1)] <= DESC.QBINS[k][1] and
        cgs_hiR[min(len(cgs_hiR)-1, ceil(((k+1)/6.0)*len(cgs_hiR))-1)] <= DESC.QBINS[k][2]
        for k in range(6)
    )

    # Lattice-to-lattice smooth-yield overdispersion CV^2 -> implied Gamma-Poisson shape r_eff
    my = sum(lat_yields) / len(lat_yields)
    vy = sum((y - my) ** 2 for y in lat_yields) / max(1, len(lat_yields) - 1)
    cv2_lat = vy / max(my * my, 1e-300)
    r_eff = 1.0 / max(cv2_lat, 1e-9)

    # Boot lattices: exact det + congruence + 256 radius-4 B_12(4) cell combinations vs ALL BOOT_BINS & 10 BOOT_SUB deciles!
    boot_ok = True; boot_offs = []
    for item in B12["boot"]:
        R_p = [int(x) for x in item["R"]]; B = item["B"]
        if abs(det_int(B)) != M.P_BN254 ** 6: boot_ok = False
        for col in range(12):
            av = [B[r][col] for r in range(6)]; bv = [B[6 + r][col] for r in range(6)]
            if M._pmul_h(R_p, zk_to_poly_mod(av, M.P_BN254), M.P_BN254) != zk_to_poly_mod(bv, M.P_BN254): boot_ok = False
        for _ in range(128):
            while True:
                x = [rng.randint(-2, 2) for _ in range(12)]
                if any(x) and sum(z * z for z in x) <= 16: break
            v = [sum(B[r][c] * x[c] for c in range(12)) for r in range(12)]
            _, lg = log2_norms_fp(v[:6], v[6:], False)
            boot_offs.append(lg - LOG2P6 - S_R4 - LOG2U6)
    boot_offs.sort()
    boot_bins_ok = (
        boot_offs[len(boot_offs) // 2] <= DESC.BOOT_BINS[0][1] and
        boot_offs[int(0.9 * len(boot_offs))] <= DESC.BOOT_BINS[1][1] and
        boot_offs[-1] <= DESC.BOOT_BINS[2][1]
    )
    boot_sub_ok = all(
        boot_offs[min(len(boot_offs)-1, ceil(((k+1)/10.0)*len(boot_offs))-1)] <= DESC.BOOT_SUB[k]
        for k in range(10)
    )

    all_env_ok = (max_exact_err < 1e-8 and qbins_ok and qbins_lo_ok and qbins_sub_ok and hiR_ok and boot_bins_ok and boot_sub_ok and r_eff > 2.0)
    return dict(
        ok=(sq_ok and boot_ok and all_env_ok),
        special_q_count=len(B12["special_q"]), special_q_bits=q_bits_list,
        boot_count=len(B12["boot"]), boot_r4_samples=len(boot_offs),
        special_q_det_and_mem_ok=sq_ok, boot_det_and_mem_ok=boot_ok,
        exact_int_vs_fp_max_err_bits=round(max_exact_err, 12),
        qbins_6_of_6_ok=qbins_ok, qbins_lo_6_of_6_ok=qbins_lo_ok,
        qbins_sub_24_of_24_ok=qbins_sub_ok, hi_radius_R16_128_ok=hiR_ok,
        boot_bins_3_of_3_ok=boot_bins_ok, boot_sub_10_of_10_ok=boot_sub_ok,
        special_q_med_cf=round(cfs_r4[len(cfs_r4)//2], 2),
        special_q_p90_cf=round(cfs_r4[int(0.9*len(cfs_r4))], 2),
        special_q_med_cg=round(cgs_r4[len(cgs_r4)//2], 2),
        boot_r4_med_off=round(boot_offs[len(boot_offs)//2], 2),
        boot_r4_p90_off=round(boot_offs[int(0.9*len(boot_offs))], 2),
        boot_r4_max_off=round(boot_offs[-1], 2),
        lattice_yield_cv2=round(cv2_lat, 4), gamma_poisson_r_eff=round(r_eff, 2),
    )

def verify_live_filtsim_z6(log2U=13, seed=20261120, W_target=200, check_kernel=False):
    """Exact pure-Python production hexagonal Z[omega] filtering simulator:
    Dickman-kappa occurrence profile at (Lf=60, Lg=59, Nf=405.5, Ng=468.7, exc=0.95),
    singleton + union-find clique purge, Markowitz merge to W=200, Boudot recorded pivots,
    and dense-band TABM grouping — measuring ALL FOUR cap metrics (N_sw/U', N_l/U', k5_exp, tab_r) live!"""
    rng = random.Random(seed)
    LN2 = math.log(2.0)
    Lf, Lg, Nf, Ng, excess = 60.0, 59.0, 405.5, 468.7, 0.95
    U6 = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))
    U6_INV = {(1, 0): (1, 0), (0, 1): (1, -1), (-1, 1): (0, -1), (-1, 0): (-1, 0), (0, -1): (-1, 1), (1, -1): (0, 1)}
    def z6_mul(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0] + x[1]*y[1])
    def z6_decomp(a, b):
        if (a >= 0 and b >= 0) or (a <= 0 and b <= 0): return (a, b, 0)
        if abs(a) <= abs(b): return (0, a + b, -a)
        return (a + b, 0, b)
    def z6_npt(v):
        x, y, z = z6_decomp(v[0], v[1])
        return (x != 0) + (y != 0) + (z != 0)
    def z6_n1(v):
        x, y, z = z6_decomp(v[0], v[1])
        return abs(x) + abs(y) + abs(z)
    def z6_mx(v):
        x, y, z = z6_decomp(v[0], v[1])
        return max(abs(x), abs(y), abs(z))
    def z6_is_unit(v): return v in U6_INV

    uf, ug = Nf / Lf, Ng / Lg
    LiBf = (2.0 ** Lf) / (Lf * LN2); LiBg = (2.0 ** Lg) / (Lg * LN2)
    LiB_sum = LiBf + LiBg; frac_f = LiBf / LiB_sum
    Up = 1 << log2U
    Up_f = int(round(Up * frac_f)); Up_g = Up - Up_f
    Rp = int(round(excess * Up))
    t_dense = 0.80
    lsmall = min(t_dense * min(Lf, Lg), max(8.0, log2(Rp / 25.0) + log2(3.0)))

    plim = int(floor(2.0 ** lsmall))
    sieve = bytearray(b'\x01') * (plim + 1)
    sieve[0:2] = b'\x00\x00'
    for i in range(2, int(plim ** 0.5) + 1):
        if sieve[i]: sieve[i*i:plim+1:i] = b'\x00' * (((plim - i*i) // i) + 1)
    ps = [i for i in range(2, plim + 1) if sieve[i]]
    nsmall_side = (len(ps) + 5) // 6
    U = Up + 2 * nsmall_side

    nt = 400
    ts = [t_dense + (1.0 - t_dense) * (i + 0.5) / nt for i in range(nt)]
    cum_f, cum_g = [], []
    tf = tg = 0.0
    for t in ts:
        tf += (2.0 ** (Lf * (t - 1.0))) / t; tg += (2.0 ** (Lg * (t - 1.0))) / t
        cum_f.append(tf); cum_g.append(tg)
    cum_f = [x / tf for x in cum_f]; cum_g = [x / tg for x in cum_g]

    kap_f = [M.kappa_rho(uf, t) for t in ts]; kap_g = [M.kappa_rho(ug, t) for t in ts]
    w = [0.0] * U; col_log2p = [0.0] * U; cid = 0; dt_cell = (1.0 - t_dense) / nt
    for _ in range(Up_f):
        idx = min(nt - 1, bisect.bisect_left(cum_f, rng.random()))
        t = ts[idx] + (rng.random() - 0.5) * dt_cell
        lam = excess * LiB_sum * kap_f[idx] / (2.0 ** (Lf * t))
        w[cid] = min(1.0, lam / Rp); col_log2p[cid] = Lf * t; cid += 1
    for _ in range(Up_g):
        idx = min(nt - 1, bisect.bisect_left(cum_g, rng.random()))
        t = ts[idx] + (rng.random() - 0.5) * dt_cell
        lam = excess * LiB_sum * kap_g[idx] / (2.0 ** (Lg * t))
        w[cid] = min(1.0, lam / Rp); col_log2p[cid] = Lg * t; cid += 1

    small = []
    for u_s, L_s in ((uf, Lf), (ug, Lg)):
        for i in range(0, len(ps), 6):
            q = ps[i]; t = log2(q) / L_s
            pr = min(1.0, 6.0 * M.kappa_rho(u_s, t) / q)
            w[cid] = pr; col_log2p[cid] = log2(q); small.append((cid, pr, q)); cid += 1

    w_mid = w_mid_28 = 0.0
    for u_s, L_s in ((uf, Lf), (ug, Lg)):
        tmin = lsmall / L_s; dt = (t_dense - tmin) / 200.0
        for i in range(200):
            t = tmin + dt * (i + 0.5); dw = M.kappa_rho(u_s, t) * dt / t
            w_mid += dw
            if L_s * t <= 28.0: w_mid_28 += dw

    dense = [i for i in range(Up) if w[i] >= 0.2]
    sparse = [i for i in range(Up) if w[i] < 0.2]
    scum = []; stot = 0.0
    for i in sparse: stot += w[i]; scum.append(stot)

    def pois(lam):
        if lam < 30.0:
            L = exp(-lam); p = 1.0; k = 0
            while True:
                p *= rng.random()
                if p <= L: return k
                k += 1
        return max(0, int(round(lam + sqrt(lam) * rng.gauss(0.0, 1.0))))

    P_CHK = 1000000009
    chi = next(x for g in range(2, 50) for x in [pow(g, (P_CHK - 1) // 6, P_CHK)] if (x * x - x + 1) % P_CHK == 0)
    def z6_eval(v): return (v[0] + v[1] * chi) % P_CHK
    v_sec = [rng.randrange(1, P_CHK) for _ in range(U)] if check_kernel else None

    rows = []; col_rows = [set() for _ in range(U)]
    for ri in range(Rp):
        d = {}
        for c in dense:
            if rng.random() < w[c]:
                u = U6[rng.randrange(6)]
                nv = (d[c][0] + u[0], d[c][1] + u[1]) if c in d else u
                if nv == (0, 0): del d[c]
                else: d[c] = nv
        for _ in range(pois(stot)):
            c = sparse[min(len(sparse) - 1, bisect.bisect_left(scum, rng.random() * stot))]
            u = U6[rng.randrange(6)]
            nv = (d[c][0] + u[0], d[c][1] + u[1]) if c in d else u
            if nv == (0, 0): del d[c]
            else: d[c] = nv
        for c, pr, q in small:
            if rng.random() < pr:
                v = 1; invq = 1.0 / q
                while rng.random() < invq: v += 1
                u = U6[rng.randrange(6)]; mv = (v * u[0], v * u[1])
                nv = (d[c][0] + mv[0], d[c][1] + mv[1]) if c in d else mv
                if nv == (0, 0): del d[c]
                else: d[c] = nv
        rw = w_mid + sum(z6_npt(v) for v in d.values())
        rhs = (-sum(z6_eval(val) * v_sec[c] for c, val in d.items())) % P_CHK if check_kernel else 0
        rows.append([d, 1, w_mid, w_mid, rw, True, rhs])
        for c in d: col_rows[c].add(ri)

    def purge_singletons(nalive):
        q_s = [c for c in range(U) if len(col_rows[c]) == 1]
        while q_s:
            c = q_s.pop()
            if len(col_rows[c]) != 1: continue
            r0 = next(iter(col_rows[c]))
            if not rows[r0][5]: continue
            rows[r0][5] = False; nalive -= 1
            for cc in rows[r0][0]:
                col_rows[cc].discard(r0)
                if len(col_rows[cc]) == 1: q_s.append(cc)
            rows[r0][0].clear()
        return nalive

    nalive = purge_singletons(Rp)
    ncols = sum(1 for c in range(U) if col_rows[c])
    texc = max(2, int(0.0001 * Up))
    while nalive - ncols > texc:
        par = list(range(Rp))
        def find(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for c in range(U):
            if len(col_rows[c]) == 2:
                r1, r2 = col_rows[c]
                a, b = find(r1), find(r2)
                if a != b: par[a] = b
        comp_w = {}
        for i in range(Rp):
            if rows[i][5]:
                rt = find(i); comp_w[rt] = comp_w.get(rt, 0.0) + rows[i][4]
        comps = sorted(comp_w.items(), key=lambda kv: -kv[1])
        need = nalive - ncols - texc
        nb = max(1, min(need // 2 + 1, max(1, len(comps) // 100)))
        kill = {rt for rt, _ in comps[:nb]}
        for i in range(Rp):
            if rows[i][5] and find(i) in kill:
                rows[i][5] = False; nalive -= 1
                for cc in rows[i][0]: col_rows[cc].discard(i)
                rows[i][0].clear()
        nalive = purge_singletons(nalive)
        ncols = sum(1 for c in range(U) if col_rows[c])

    kmax, cmax = 64, 255
    totw = sum(rows[i][4] for i in range(Rp) if rows[i][5])
    def col_key(c):
        wgt = len(col_rows[c])
        if wgt <= 2: return float(wgt - 3)
        best = 1e18
        for r in col_rows[c]:
            if z6_is_unit(rows[r][0][c]) and rows[r][4] < best: best = rows[r][4]
        return float(wgt - 2) * best
    heap = [(col_key(c), c) for c in range(U) if 1 <= len(col_rows[c]) <= kmax]
    heapq.heapify(heap)
    merges = 0
    while heap and (totw / max(1, nalive)) < W_target:
        key, c = heapq.heappop(heap)
        wgt = len(col_rows[c])
        if wgt == 0 or wgt > kmax: continue
        ck = col_key(c)
        if ck > key + 1e-9: heapq.heappush(heap, (ck, c)); continue
        if wgt == 1:
            r0 = next(iter(col_rows[c]))
            rows[r0][5] = False; nalive -= 1; totw -= rows[r0][4]
            for cc in rows[r0][0]:
                col_rows[cc].discard(r0)
                if 1 <= len(col_rows[cc]) <= kmax: heapq.heappush(heap, (col_key(cc), cc))
            rows[r0][0].clear(); continue
        piv = -1; best_rw = 1e18
        for r in col_rows[c]:
            if z6_is_unit(rows[r][0][c]) and rows[r][4] < best_rw: best_rw = rows[r][4]; piv = r
        if piv < 0: continue
        pr = rows[piv]; pd = pr[0]; pc_inv = U6_INV[pd[c]]
        snap = [r for r in col_rows[c] if r != piv]
        ok = True; delta_w = -pr[4]
        for i in snap:
            ri = rows[i]; m = z6_mul(ri[0][c], pc_inv); mnpt = z6_npt(m)
            if ri[1] + z6_n1(m) * pr[1] > cmax: ok = False; break
            delta_w += pr[2] if mnpt == 1 else mnpt * pr[3]
            rd = ri[0]
            for cb, pb in pd.items():
                ov = rd.get(cb); mp = z6_mul(m, pb)
                nv = (-mp[0], -mp[1]) if ov is None else (ov[0] - mp[0], ov[1] - mp[1])
                if z6_mx(nv) > cmax: ok = False; break
                delta_w += z6_npt(nv) - (0 if ov is None else z6_npt(ov))
            if not ok: break
        if not ok: continue
        # Exact deterministic post-merge weight guard: never commit a merge pushing avg_w > W_ROW (201.0)!
        if (totw + delta_w) > M.W_ROW * max(1, nalive - 1):
            break
        for i in snap:
            ri = rows[i]; rd = ri[0]; m = z6_mul(rd[c], pc_inv); mnpt = z6_npt(m)
            ri[1] += z6_n1(m) * pr[1]
            dm = pr[2] if mnpt == 1 else mnpt * pr[3]
            ri[2] += dm; ri[3] += pr[3]; totw += dm
            if check_kernel: ri[6] = (ri[6] - z6_eval(m) * pr[6]) % P_CHK
            for cb, pb in pd.items():
                ov = rd.get(cb); mp = z6_mul(m, pb)
                if ov is None:
                    nv = (-mp[0], -mp[1]); rd[cb] = nv; col_rows[cb].add(i); totw += z6_npt(nv)
                else:
                    nv = (ov[0] - mp[0], ov[1] - mp[1]); totw += z6_npt(nv) - z6_npt(ov)
                    if nv == (0, 0):
                        del rd[cb]; col_rows[cb].discard(i)
                        if 1 <= len(col_rows[cb]) <= kmax: heapq.heappush(heap, (col_key(cb), cb))
                    else: rd[cb] = nv
            ri[4] = ri[2] + sum(z6_npt(v) for v in rd.values())
            for cb in rd:
                if len(col_rows[cb]) <= kmax: heapq.heappush(heap, (col_key(cb), cb))
            merges += 1
        pr[5] = False; nalive -= 1; totw -= pr[4]
        for cc in pd:
            col_rows[cc].discard(piv)
            if 1 <= len(col_rows[cc]) <= kmax: heapq.heappush(heap, (col_key(cc), cc))
        pd.clear()

    alive = [i for i in range(Rp) if rows[i][5]]
    N_sw = len(alive)
    kern_errs = 0
    if check_kernel:
        for r in alive:
            if (sum(z6_eval(val) * v_sec[c] for c, val in rows[r][0].items()) + rows[r][6]) % P_CHK != 0:
                kern_errs += 1
    occ = sorted([(len(col_rows[c]), c) for c in range(Up) if col_rows[c]])
    dc = {}; lab = {}; mn = {}; nrec = 0
    for _, c in occ:
        hi = mn.get(c, 1.0)
        cands = sorted([(rows[r][4], r) for r in col_rows[c] if r not in dc and z6_is_unit(rows[r][0][c])])
        for _, r in cands:
            lo = 0.0
            for x in rows[r][0]:
                if x != c and x in lab and lab[x] > lo: lo = lab[x]
            v = lo + 0.5 * (hi - lo)
            if lo < v < hi:
                dc[r] = c; lab[c] = v; nrec += 1
                for x in rows[r][0]:
                    if x != c and mn.get(x, 1.0) > v: mn[x] = v
                break
    N_l = N_sw - nrec
    U5, G1 = M.C_UNIT_5, M.C_GEN_1
    exp_terms = exp_ones = tab_plain = tab_grp = 0.0
    mid28_frac = w_mid_28 / max(w_mid, 1e-12)
    for r in alive:
        piv_c = dc.get(r, -1); nb4 = nb3 = nb2 = 0; rplain = 0.0
        for c, val in rows[r][0].items():
            lp = col_log2p[c]; is_p = (c == piv_c)
            for t in z6_decomp(val[0], val[1]):
                v = abs(t)
                if v == 0: continue
                exp_terms += 1.0
                if v == 1: exp_ones += 1.0
                tc = (U5 * v) if v <= 4 else G1
                if not is_p:
                    rplain += tc
                    if v <= 4 and lp <= 28.0:
                        if lp <= 14.0: nb4 += v
                        elif lp <= 18.25: nb3 += v
                        else: nb2 += v
        u_cols = rows[r][3]; rplain += u_cols * U5
        if mid28_frac > 0.0 and u_cols > 0.0: nb2 += pois(u_cols * mid28_frac)
        g4 = nb4 // 4; rem4 = nb4 - 4 * g4
        n3 = nb3 + rem4; g3 = n3 // 3; rem3 = n3 - 3 * g3
        g2 = (nb2 + rem3) // 2
        saved = (3 * g4 + 2 * g3 + g2) * U5
        tab_plain += rplain; tab_grp += (rplain - saved)
    exp_unit = exp_ones / max(exp_terms, 1.0)
    k5_exp = exp_unit * U5 + (1.0 - exp_unit) * G1
    tab_r = tab_grp / max(tab_plain, 1e-12)
    avg_w_fin = totw / max(1, N_sw)
    return dict(
        ok=(kern_errs == 0 and N_sw > 0 and N_l > 0 and merges > 0 and avg_w_fin <= M.W_ROW and tab_r < M.TAB_R_CAP),
        log2U=log2U, seed=seed, W_target=W_target, initial_rows=Rp,
        N_sw=N_sw, N_l=N_l, N_sw_U=round(N_sw / Up, 6), N_l_U=round(N_l / Up, 6),
        avg_w=round(avg_w_fin, 2), merges=merges,
        exp_unit=round(exp_unit, 5), k5_exp=round(k5_exp, 4), tab_r=round(tab_r, 5),
        kernel_errors=kern_errs,
    )

def verify_live_filtsim_ladder(seed=20261120):
    """Run exact production filtsim_z6 live at BOTH U'=2^13 (with kernel check) and U'=2^14!"""
    r13 = verify_live_filtsim_z6(log2U=13, seed=seed, check_kernel=True)
    r14 = verify_live_filtsim_z6(log2U=14, seed=seed, check_kernel=False)
    mono = (r14["N_sw_U"] < r13["N_sw_U"] and r14["N_l_U"] < r13["N_l_U"] and
            r14["k5_exp"] < r13["k5_exp"] and r14["tab_r"] < r13["tab_r"])
    return dict(ok=(r13["ok"] and r14["ok"] and mono), monotonic_13_to_14=mono, u13=r13, u14=r14)
