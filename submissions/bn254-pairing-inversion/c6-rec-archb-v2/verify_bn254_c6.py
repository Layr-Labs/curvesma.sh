"""
verify_bn254_c6.py — Self-contained pure-stdlib live prototypes for BN254 v2:
  1. Band-3 7D LLL + weighted-cylinder Fincke-Pohst coset-walk prototype (0 missed hits + nodes/leaves bounds);
  2. Exact 12D C6 special-q (2^65..2^125) & boot lattice determinant, congruence, and norm-envelope audit;
  3. Pure-Python hexagonal Z[omega] singleton + clique filtering simulator with exact kernel verification.
"""
import random
from math import ceil, floor, log2, sqrt
import model_bn254_c6 as M

def verify_7d_walk_prototype():
    """Live 7D LLL + weighted-cylinder Fincke-Pohst coset walk on 6 real C6 ideals x 3 cosets!"""
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

    w_ideals = []
    for p in range(1009, 3000, 2):
        if not is_prime(p): continue
        for s in range(p):
            if eval_poly_mod(M.H_COEFFS, s, p) == 0:
                uw = ((M.U_BN254 % p) * eval_poly_mod(M.W_COEFFS, s, p)) % p
                if uw == 0: continue
                inv_uw = pow(uw, -1, p)
                if pow(inv_uw, (p - 1) // 2, p) == 1:
                    rho_r = next(r for r in range(1, (p + 1) // 2) if (r * r) % p == inv_uw)
                    w_ideals.append((p, s, rho_r)); break
        if len(w_ideals) >= 6: break

    L_a, RB2, wt_t = 32, 85.0, 0.40
    R7_sq = 1.0 + wt_t * wt_t
    L6 = [[0.0] * 6 for _ in range(6)]
    for r in range(6):
        for c in range(r + 1):
            sv = GB[r][c] - sum(L6[r][k] * L6[c][k] for k in range(c))
            L6[r][c] = sqrt(sv) if r == c else sv / L6[c][c]
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
        for a0_mod in (0, 17, 123):
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
    lph = leaves[0] / max(1, hits); nph = nodes[0] / max(1, hits)
    return dict(
        ok=(missed == 0 and hits > 0 and lph <= M.FP_LEAVES_PER_HIT and nph <= M.FP_NODES_PER_HIT),
        n_ideals=len(w_ideals), cosets_per_ideal=3, total_hits=hits, missed=missed,
        leaves_per_hit=round(lph, 4), nodes_per_hit=round(nph, 4),
        bound_leaves_per_hit=M.FP_LEAVES_PER_HIT, bound_nodes_per_hit=M.FP_NODES_PER_HIT,
    )

def verify_12d_descent_lattices():
    """Live big-int det/membership + norm-envelope check on exact C6 12D special-q (2^65..2^125) & boot bases!"""
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
    def eval_zk_mod(vec, s, mod):
        inv3 = pow(3, -1, mod); sp = [1]
        for _ in range(5): sp.append((sp[-1] * s) % mod)
        return (sum(vec[k] * (sum(M.ZK_NUM[k][j] * sp[j] for j in range(6)) % mod) for k in range(6) if vec[k]) * inv3) % mod
    def zk_to_poly_mod(vec, mod):
        inv3 = pow(3, -1, mod)
        return [(sum(vec[k] * M.ZK_NUM[k][j] for k in range(6)) * inv3) % mod for j in range(6)]
    def log2_norms(a, b, need_f=True):
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
    sq_ok = True; cfs, cgs = [], []
    for item in B12["special_q"]:
        q, s, rho_r = int(item["q"]), int(item["s"]), int(item["rho"])
        B = item["B"]
        if abs(det_int(B)) != q: sq_ok = False
        for col in range(12):
            av = [B[r][col] for r in range(6)]; bv = [B[6 + r][col] for r in range(6)]
            if (eval_zk_mod(av, s, q) - rho_r * eval_zk_mod(bv, s, q)) % q != 0: sq_ok = False
        c_r = (log2(q) + S_R4) / 12.0
        for _ in range(60):
            while True:
                x = [rng.randint(-2, 2) for _ in range(12)]
                if any(x) and sum(z * z for z in x) <= 16: break
            v = [sum(B[r][c] * x[c] for c in range(12)) for r in range(12)]
            lf, lg = log2_norms(v[:6], v[6:], True)
            cfs.append(lf - 48.0 * c_r); cgs.append(lg - 12.0 * c_r - LOG2U6)

    boot_ok = True; boot_offs = []
    for item in B12["boot"]:
        R_p = [int(x) for x in item["R"]]; B = item["B"]
        if abs(det_int(B)) != M.P_BN254 ** 6: boot_ok = False
        for col in range(12):
            av = [B[r][col] for r in range(6)]; bv = [B[6 + r][col] for r in range(6)]
            if M._pmul_h(R_p, zk_to_poly_mod(av, M.P_BN254), M.P_BN254) != zk_to_poly_mod(bv, M.P_BN254): boot_ok = False
            _, lg = log2_norms(av, bv, False)
            boot_offs.append(lg - LOG2P6 - LOG2U6)

    cfs.sort(); cgs.sort(); boot_offs.sort()
    med_cf = cfs[len(cfs) // 2]; p90_cf = cfs[int(0.9 * len(cfs))]; med_cg = cgs[len(cgs) // 2]
    med_boot = boot_offs[len(boot_offs) // 2]; min_boot = boot_offs[0]
    env_ok = (med_cf < 31.5 and p90_cf < 47.6 and med_cg < 14.6 and med_boot < -178.78)
    return dict(
        ok=(sq_ok and boot_ok and env_ok),
        special_q_count=len(B12["special_q"]), boot_count=len(B12["boot"]),
        special_q_det_and_mem_ok=sq_ok, boot_det_and_mem_ok=boot_ok,
        special_q_med_cf=round(med_cf, 2), special_q_p90_cf=round(p90_cf, 2), special_q_med_cg=round(med_cg, 2),
        boot_min_off=round(min_boot, 2), boot_med_off=round(med_boot, 2),
        beats_conservative_qbins_and_boot_bins=env_ok,
    )

def verify_live_filtsim_z6(log2U=13, seed=20261120):
    """Live pure-Python hexagonal Z[omega] singleton + clique filtering with exact kernel check!"""
    U_s = 1 << log2U; n_rel = int(0.95 * U_s); rng = random.Random(seed); P_CHK = 1000000009
    chi = next(x for g in range(2, 50) for x in [pow(g, (P_CHK - 1) // 6, P_CHK)] if (x * x - x + 1) % P_CHK == 0)
    UNITS = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
    def z6_eval(c): return (c[0] + c[1] * chi) % P_CHK
    def z6_mul(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0] + x[1]*y[1])
    def z6_add(x, y): return (x[0] + y[0], x[1] + y[1])
    def z6_inv(u): return next(v for v in UNITS if z6_mul(u, v) == (1, 0))

    v_sec = [rng.randrange(1, P_CHK) for _ in range(U_s)]
    rows, col_rows = [], [set() for _ in range(U_s)]
    for ri in range(n_rel):
        deg = rng.randint(5, 11); cols = set()
        while len(cols) < deg:
            c = int(U_s * (rng.random() ** 2.2))
            if c < U_s: cols.add(c)
        d = {c: UNITS[rng.randrange(6)] for c in cols}
        rhs = sum(z6_eval(val) * v_sec[c] for c, val in d.items()) % P_CHK
        rows.append([d, (-rhs) % P_CHK])
        for c in d: col_rows[c].add(ri)

    alive_r = set(range(n_rel))
    queue = [c for c in range(U_s) if 1 <= len(col_rows[c]) <= 2]; in_q = set(queue)
    merges = singletons = 0
    while queue:
        c = queue.pop(); in_q.discard(c)
        rs = [r for r in col_rows[c] if r in alive_r]
        if len(rs) == 1:
            r0 = rs[0]; alive_r.discard(r0); singletons += 1
            for cc in list(rows[r0][0].keys()):
                col_rows[cc].discard(r0)
                if 1 <= len(col_rows[cc]) <= 2 and cc not in in_q: queue.append(cc); in_q.add(cc)
        elif len(rs) == 2:
            r0, r1 = rs; u0, u1 = rows[r0][0][c], rows[r1][0][c]
            if u0 in UNITS and u1 in UNITS:
                fac = z6_mul(u1, z6_inv(u0)); neg_fac = (-fac[0], -fac[1])
                d0, b0 = rows[r0]; d1, b1 = rows[r1]
                if len(d0) + len(d1) <= 120:
                    alive_r.discard(r0); merges += 1
                    for cc, val in d0.items():
                        col_rows[cc].discard(r0)
                        nv = z6_add(d1.get(cc, (0, 0)), z6_mul(neg_fac, val))
                        if nv == (0, 0): d1.pop(cc, None); col_rows[cc].discard(r1)
                        else: d1[cc] = nv; col_rows[cc].add(r1)
                        if 1 <= len(col_rows[cc]) <= 2 and cc not in in_q: queue.append(cc); in_q.add(cc)
                    rows[r1][1] = (b1 + z6_eval(neg_fac) * b0) % P_CHK

    kern_errs = tot_terms = unit_terms = 0
    for r in alive_r:
        d, b = rows[r]
        if (sum(z6_eval(val) * v_sec[c] for c, val in d.items()) + b) % P_CHK != 0: kern_errs += 1
        for val in d.values():
            tot_terms += 1
            if val in UNITS: unit_terms += 1
    u_share = unit_terms / max(1, tot_terms)
    return dict(
        ok=(kern_errs == 0 and merges > 0 and u_share > 0.90),
        log2U=log2U, seed=seed, initial_rows=n_rel, alive_rows=len(alive_r),
        singletons_pruned=singletons, z6_unit_merges=merges,
        unit_share=round(u_share, 5), kernel_errors=kern_errs,
    )
