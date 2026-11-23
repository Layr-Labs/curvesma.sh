"""v25_verify.py — Pure-stdlib live verifier for BN254 C6 v25 upgrades:
1. Exact closed-form K_B1=512 B1_PREF_COV quadrature (<= 1.0060).
2. Band-1 [2^24, 2^29.5) uint32_t delta-encoded runs (APPLY256_DELTA32_LSE, bit-exact).
3. Fixed a0-chord-stepping kta(a) with explicit Pass-0 ai_of_chord -> kta_arr[ai] scatter preserving exact Q(a)-sorted ai order (3,072 genuine 42..48b (p,a) pairs + end-to-end ai-order x4/ll build check).
4. 4-plane cyclic Cone3 algebraic invariance under S[3..5] + 450-point bijection check.
5. k=1->k=0 finite-diff rem_0 one-sided safety + unconditional r2_ub < 2^40 certificate across [48..61]b.
6. 32 KiB-strided single 16 KiB staging slot tbnz #14 exact replay (zero double-buffer misfire).
7. Lower Band 2 x4 (ai << 17) | j_rel + add lsr #16 encoding across variable k_s in {35, 42, 48, 50}, W in {8, 16}, sgn in {+1, -1}.
8. W=64 floor_div b_0 sub-band nesting at b_0 in {-64, -1, 0, 63, 64}, Graham greedy deterministic per-C[s] bound, derived over-radius row-capacity guard (delta = 0.005), dummy sink block disjointness on C[0], and init_a + read-only scan8 bit-exact check.
"""
import math, random
from math import log, log2, sqrt, ceil, exp
import model_bn254_v6 as M
import v20_lat_b2 as V20
import descent_bn254_v6 as DESC

M64 = (1 << 64) - 1

def verify_v25_upgrades(seed=20261122):
    rng = random.Random(seed)
    T4 = M.load_data()["t4_geom"]
    G_sym, _, _, S, S2, G_SC = M.build_exact_int_geom(T4)

    # 1. Exact closed-form piecewise-exact B1_PREF_COV integration for B_a=4, K_B1 in {64, 512}
    def _cov_exact(K):
        tot = 0.0
        n = 4 * K
        inv_n = 1.0 / n
        for m in range(1, n + 1):
            u_lo = 1.0 - (m * inv_n) ** (1.0 / 3.0)
            u_hi = 1.0 - ((m - 1) * inv_n) ** (1.0 / 3.0)
            tot += (m * inv_n) * (u_hi ** 3 - u_lo ** 3)
        return tot / 0.05
    cov_64 = _cov_exact(64)
    cov_512 = _cov_exact(M.K_B1_PREF)
    b1_cov_ok = (cov_64 <= 1.0485 and cov_512 <= M.B1_PREF_COV)

    # 2. Band-1 [2^24, 2^29.5) uint32_t delta-encoded run verification (APPLY256_DELTA32_LSE)
    p_d32 = 29.5
    steps = 500
    dt = (p_d32 - 24.0) / steps
    w_tot = brk_tot = max_qb = 0.0
    for i in range(steps):
        t = 24.0 + (i + 0.5) * dt
        qb = exp(-(2.0 ** 32) / (2.0 ** t))
        w = dt / t
        w_tot += w
        brk_tot += w * qb
        if qb > max_qb:
            max_qb = qb
    avg_qb = brk_tot / w_tot
    c_d32 = 642.0 / 256.0 + 0.525 + avg_qb * 140.1
    frac_d32 = log(p_d32 / 24.0) / log(M.OP["p_cut"] / 24.0)
    c_hit_lo_exact = frac_d32 * c_d32 + (1.0 - frac_d32) * (642.0 / 256.0 + 1.05)
    b1_d32_cost_ok = (c_d32 <= 3.06 and c_hit_lo_exact <= M.C_HIT_LO)

    # Bit-exact APPLY256_DELTA32_LSE replay (including >= 2^32 gap run break)
    d32_errs = 0
    SCORE_N = 1024
    for trial in range(5):
        idxs = sorted(rng.sample(range(1, SCORE_N), 320))
        ref_sc = [rng.randrange(0, 32768) for _ in range(SCORE_N)]
        got_sc = ref_sc[:]
        w_lp = 400 + trial * 17
        for ai in idxs:
            ref_sc[ai] = (ref_sc[ai] + w_lp) & 0xFFFF
        runs = []
        cur_start = idxs[0]
        cur_deltas = []
        prev = idxs[0]
        max_d = 12 if trial == 0 else ((1 << 32) - 1)
        for ai in idxs[1:]:
            d = ai - prev
            if d > max_d:
                runs.append((cur_start, cur_deltas))
                cur_start = ai
                cur_deltas = []
            else:
                cur_deltas.append(d)
            prev = ai
        runs.append((cur_start, cur_deltas))
        for start_ai, deltas in runs:
            cur = start_ai
            V20._staddh_u16(got_sc, cur, w_lp)
            j = 0
            nd = len(deltas)
            while j + 1 < nd:
                d0, d1 = deltas[j], deltas[j + 1]
                assert 0 < d0 < (1 << 32) and 0 < d1 < (1 << 32)
                cur += d0
                V20._staddh_u16(got_sc, cur, w_lp)
                cur += d1
                V20._staddh_u16(got_sc, cur, w_lp)
                j += 2
            if j < nd:
                cur += deltas[j]
                V20._staddh_u16(got_sc, cur, w_lp)
        if got_sc != ref_sc:
            d32_errs += 1

    # 3. Fixed a0-chord-stepping kta(a) with explicit Pass-0 ai_of_chord -> kta_arr[ai] scatter
    # preserving exact Q(a)-sorted ai order across 3,072 genuine 42..48b (p, a) pairs
    chord_errs = 0
    chord_tested = 0
    chord_order_errs = 0
    for pb in (42, 46, 48):
        p = (1 << pb) + 1
        while True:
            if (p & 7) in (1, 3) and M._isp(p):
                s0 = V20._find_root_h(p)
                if s0 is not None:
                    ws = V20._eval_poly(M.W_COEFFS, s0, p)
                    if ws != 0:
                        uw = ((M.U_BN254 % p) * ws) % p
                        rho = V20._sqrt_mod(pow(uw, -1, p), p)
                        if rho > 0:
                            break
            p += 2
        inv3 = pow(3, -1, p)
        sp = [1]
        for _ in range(5):
            sp.append((sp[-1] * s0) % p)
        ek = [(sum(M.ZK_NUM[k][j] * sp[j] for j in range(6)) * inv3) % p for k in range(6)]
        neg_ir = (-pow(rho, -1, p)) % p
        ck = [(e * neg_ir) % p for e in ek]
        c0 = ck[0]
        # Build 32 chords of length 32 = 1,024 a-vectors, assign Q(a)-sorted ai indices,
        # and verify Pass-0 chord walk + scatter into kta_arr[ai] matches exact (a.c) mod p in ai order!
        all_vecs = []
        chords = []
        for ch_i in range(32):
            a_tail = [rng.randint(-300, 300) for _ in range(5)]
            a0_lo = rng.randint(-200, 150)
            ch_vecs = [[a0_lo + step] + a_tail for step in range(32)]
            chords.append(ch_vecs)
            all_vecs.extend(ch_vecs)
        # Sort all_vecs by Q(a) to define exact ai = 0..1023
        qa_key = lambda v: math.sumprod(V20._mv6(G_sym, v), v)
        order = sorted(range(len(all_vecs)), key=lambda idx: (qa_key(all_vecs[idx]), all_vecs[idx]))
        ai_of_flat = [0] * len(all_vecs)
        for ai, fidx in enumerate(order):
            ai_of_flat[fidx] = ai
        # Pass 0: walk a0-chords, step kta by add/cmp/csel (3 ops), scatter to kta_arr[ai]
        kta_arr = [-1] * len(all_vecs)
        flat_pos = 0
        for ch_vecs in chords:
            a_start = ch_vecs[0]
            kta = sum(a_start[k] * ck[k] for k in range(6)) % p
            for v in ch_vecs:
                ai = ai_of_flat[flat_pos]
                flat_pos += 1
                kta_arr[ai] = kta
                ref = sum(v[k] * ck[k] for k in range(6)) % p
                if kta != ref:
                    chord_errs += 1
                chord_tested += 1
                nxt = kta + c0
                kta = nxt - p if nxt >= p else nxt
        # Pass 1/2 check: read kta_arr[ai] in exact ascending ai order
        for ai, fidx in enumerate(order):
            v = all_vecs[fidx]
            ref = sum(v[k] * ck[k] for k in range(6)) % p
            if kta_arr[ai] != ref:
                chord_order_errs += 1

    # 3b. Mean a_0-chord piece length (continuum; Q(a) <= 1 is the ellipsoid 1764 a^T G a <= T12 = 1764 G_SC R_a^2):
    # half-length of a full chord at projected radius r: R_a sqrt(G_SC / G_00) sqrt(1 - r^2); lines uniform in the
    # projected 5-ball (density 5 r^4); B_a = 4 equal-count shells have Q-thresholds (j/4)^(1/3), and a line at r
    # crosses each threshold above r^2 twice: pieces = 1 + 2 #{j : s_j > r^2}.
    half0 = 589.09 * sqrt(196608.0 / 753729.0)
    n_int = 20000
    e_len = sum(5.0 * ((i + 0.5) / n_int) ** 4 * sqrt(1.0 - ((i + 0.5) / n_int) ** 2) for i in range(n_int)) / n_int
    chord_mean_pts = 2.0 * half0 * e_len
    pieces_mean = 1.0 + 2.0 * sum(((j / 4.0) ** (1.0 / 3.0)) ** 2.5 for j in (1, 2, 3))
    piece_mean_pts = chord_mean_pts / pieces_mean
    if not (piece_mean_pts >= 34.0):
        chord_errs += 1

    # 4. 4-plane cyclic Cone3 algebraic invariance under S[3..5] + 450-point bijection check
    # Verify S maps (a3, -a4, a5) -> (-a4, a5, a3) identically (hence sx = a3 - a4 + a5 is S-invariant!)
    s35_alg_ok = (S[3] == [0, 0, 0, 0, -1, 0] and S[4] == [0, 0, 0, 0, 0, -1] and S[5] == [0, 0, 0, 1, 0, 0])
    cone3_errs = 0
    cone3_tested = 0
    for _ in range(450):
        a = [rng.randint(-300, 300) for _ in range(6)]
        b = [rng.randint(-400, 400) for _ in range(6)]
        v = a + b
        orb = []
        cur_a, cur_b = a[:], b[:]
        for _ in range(3):
            for sa in (1, -1):
                for sb in (1, -1):
                    orb.append([sa * x for x in cur_a] + [sb * x for x in cur_b])
            cur_a = V20._mv6(S, cur_a)
            cur_b = V20._mv6(S, cur_b)
        sx0 = a[3] - a[4] + a[5]
        s_a = V20._mv6(S, a)
        if (s_a[3] - s_a[4] + s_a[5]) != sx0:
            cone3_errs += 1
        # Check exact 1-to-1 bijection on non-boundary orbits (sx != 0, strict max, b0 != 0)
        c0 = V20._cone_of(v)
        if c0 >= 0:
            ca0, cb0 = V20._canon_cone_ref(tuple(v[:6]), tuple(v[6:]), S, S2)
            can0 = list(ca0) + list(cb0)
            hits = [u for u in orb if (u[3] - u[4] + u[5] > 0 and u[3] > -u[4] and u[3] > u[5] and u[6] > 0)]
            if len(hits) != 1 or hits[0] != can0:
                cone3_errs += 1
            for u in orb:
                cau, cbu = V20._canon_cone_ref(tuple(u[:6]), tuple(u[6:]), S, S2)
                if (cau, cbu) != (ca0, cb0):
                    cone3_errs += 1
            cone3_tested += 1

    # 5. Unconditional k=1 -> k=0 finite-diff rem_0 upper bound certificate (r2_ub < 2^40) across [48..61]b
    T12_prod = int(round(1764 * (589.43 ** 2) * G_SC))
    Q38 = 1 << 38
    fd_rem0_max_r2_bits = 0.0
    fd_rem0_errs = 0
    fd_A1_min = 1 << 64
    for pb in (48, 52, 56, 61):
        pp = (1 << pb) + 1
        while True:
            if (pp & 7) in (1, 3) and M._isp(pp):
                s0p = V20._find_root_h(pp)
                if s0p is not None:
                    wsp = V20._eval_poly(M.W_COEFFS, s0p, pp)
                    if wsp != 0:
                        uwp = ((M.U_BN254 % pp) * wsp) % pp
                        rhop = V20._sqrt_mod(pow(uwp, -1, pp), pp)
                        if rhop > 0:
                            break
            pp += 2
        inv3p = pow(3, -1, pp)
        spp = [1]
        for _ in range(5):
            spp.append((spp[-1] * s0p) % pp)
        esp = [(sum(M.ZK_NUM[k][j] * spp[j] for j in range(6)) * inv3p) % pp for k in range(6)]
        lfp = esp + [((-rhop) * e) % pp for e in esp]
        inv0p = pow(lfp[0], -1, pp)
        Bp = [[pp if c == 0 else 0 for c in range(12)]]
        for r1 in range(1, 12):
            c0_v = (-lfp[r1] * inv0p) % pp
            if c0_v > pp // 2:
                c0_v -= pp
            Bp.append([c0_v if c == 0 else (1 if c == r1 else 0) for c in range(12)])
        _, Dp, _ = DESC.bareiss_lll_12d(Bp, G_sym, 1764, 1000, False, 9, 10, (99, 100), True)
        inv0 = ((Dp[0] * T12_prod << 38) + Dp[1] - 1) // Dp[1]
        A1_fd = (Dp[2] << 38) // (Dp[1] * T12_prod)
        fd_A1_min = min(fd_A1_min, A1_fd)
        Xm = [0] * 12
        for k in range(11, -1, -1):
            rk = math.isqrt(((Dp[k] * T12_prod) + Dp[k + 1] - 1) // Dp[k + 1]) + 1
            Xm[k] = rk + sum((Xm[j] + 1) // 2 for j in range(k + 1, 12))
        ec1 = sum(Xm[j] for j in range(2, 12)) + 1
        R1_max = 2 * Xm[1] + 2
        pad1_max = ((2 * R1_max * ec1 * A1_fd) >> 32) + 2
        rem0_ub = Q38 + pad1_max + R1_max * R1_max
        r2_ub = ((rem0_ub * inv0) >> 64) + 1
        if r2_ub >= (1 << 40):
            fd_rem0_errs += 1
        fd_rem0_max_r2_bits = max(fd_rem0_max_r2_bits, log2(r2_ub))

    # 6. 32 KiB-strided single 16 KiB staging slot tbnz #14 exact replay (zero double-buffer misfire!)
    # Each s_id has its single 16 KiB staging slot at base + (s_id << 15), starting at +16.
    # For records 1..2045, post-increment wp has offset 24..16376 -> bit 14 == 0.
    # After record 2046, post-increment wp has offset 16384 -> bit 14 == 1 -> .Lflush_s copies 16 KiB to link_E_to_C[s_id] and resets wp to +16!
    tbnz_errs = 0
    base_stage = 1 << 40
    flushes_seen = 0
    for s_id in (0, 1, 17, (1 << 36) - 1):
        slot_base = base_stage + (s_id << 15)
        wp = slot_base + 16
        for rec_i in range(1, 4100):
            wp += 8
            bit14 = (wp >> 14) & 1
            should_flush = ((rec_i % 2046) == 0)
            if bool(bit14) != should_flush:
                tbnz_errs += 1
            if bit14:
                if wp != slot_base + 16384:
                    tbnz_errs += 1
                wp = slot_base + 16
                flushes_seen += 1
    if flushes_seen != 4 * 2:
        tbnz_errs += 1

    # 7. Lower Band 2 x4 (ai << 17) | j_rel + add lsr #16 encoding across variable k_s, W in {8,16}, sgn in {+1,-1}
    sh17_x4_errs = 0
    for W_t in (8, 16):
        for k_s in (35, 42, 48, 50):
            row_stride = 1 << k_s
            namax_t = 4000
            namax_sh17 = (namax_t << 17) | 15
            # Build 200 sorted entries (ai in [1..5000], j_rel in [0..W_t-1])
            pairs = sorted((rng.randint(1, 5000), rng.randrange(W_t)) for _ in range(200))
            entries = [(ai << 17) | jr for ai, jr in pairs]
            for sgn_v in (1, -1):
                got_map = {}
                ref_map = {}
                base_addr = (1 << 55) if sgn_v > 0 else ((1 << 55) + (W_t - 1) * row_stride)
                for e in entries:
                    if e > namax_sh17:
                        break
                    jr = e & 0xF
                    r_addr = (base_addr + (jr << k_s)) if sgn_v > 0 else (base_addr - (jr << k_s))
                    a_addr = r_addr + (e >> 16)
                    got_map[a_addr] = (got_map.get(a_addr, 0) + 720) & 0xFFFF
                for ai, jr in pairs:
                    if ai > namax_t:
                        break
                    jr_eff = jr if sgn_v > 0 else (W_t - 1 - jr)
                    a_addr = (1 << 55) + (jr_eff << k_s) + (ai << 1)
                    ref_map[a_addr] = (ref_map.get(a_addr, 0) + 720) & 0xFFFF
                if got_map != ref_map:
                    sh17_x4_errs += 1
            # Also verify exact boundary ai == namax_t (j_rel = 15) passes and ai == namax_t + 1 (j_rel = 0) breaks!
            if ((namax_t << 17) | 15) > namax_sh17 or (((namax_t + 1) << 17) | 0) <= namax_sh17:
                sh17_x4_errs += 1

    # 8. W=64 floor_div b_0 sub-band nesting, derived row-capacity guard, dummy sink disjointness, and init_a + scan8 check
    sb_nest_errs = 0
    for b0 in (-64, -63, -1, 0, 1, 63, 64, 559):
        sb64 = b0 >> 6
        for W_sub in (8, 16, 32, 64):
            sh = int(round(log2(W_sub)))
            lo_b0 = (b0 >> sh) << sh
            hi_b0 = lo_b0 + W_sub - 1
            if (lo_b0 >> 6) != sb64 or (hi_b0 >> 6) != sb64:
                sb_nest_errs += 1
    # Over-radius bound eps_or (KERNELS.md 1.3 Theorem B): input-uniform, from fd_bounds.py (exact rationals).
    import fd_bounds as FDB
    fdb = FDB.main()
    import exact_levels as EXL
    exl = EXL.main()
    eps_or = fdb["eps_or_bound"]
    guard_errs = 0
    if not (eps_or <= M.SC_GUARD_DELTA / 2.0):
        guard_errs += 1
    # Row capacity and storage: exact-count capacity incl. the Band-1 checkpoint spacing D_B1, outward
    # (input-independent) bounds on rows, per-node RAM and machine total: stor_bound.py (exact rationals).
    import stor_bound as SBND
    sbo = SBND.main()
    if not (sbo["ok"] and fdb["ok"] and eps_or <= float(SBND.EPS_OR_MAX) and float(SBND.EPS_OR_MAX) <= float(SBND.DELTA) / 2
            and fdb["exact_per_orbit_ir_log2"] <= log2(M.C_LAT_EXC) and exl["ok"]):
        guard_errs += 1
    dummy_sink_ok = bool(sbo["sink_ok"])
    # 8-row group scan/init length: every row of an 8-row group is initialised and scanned to the group's
    # longest logical row L_g; extra cells (Sum 8 L_g / Sum N_a(b) - 1) along b_0-chords, N_a(b) ~ (1-Q(b))^3,
    # chords weighted by their 5-dim projection, 120 stratified chords x 8 phases.
    num_g = den_g = 0.0
    R0 = 551.94
    for i_r in range(120):
        qp = ((i_r + 0.5) / 120.0) ** 0.4
        w_c = sqrt(1.0 - qp) * R0
        for i_c in range(8):
            c_c = 8.0 * (i_c + 0.5) / 8.0
            grp = {}
            for b in range(ceil(c_c - w_c), int(c_c + w_c) + 1):
                v = max(0.0, 1.0 - qp - ((b - c_c) / R0) ** 2) ** 3
                grp.setdefault(b >> 3, []).append(v)
            for vs in grp.values():
                num_g += 8.0 * max(vs)
                den_g += sum(vs)
    scan_group_overhead = num_g / den_g
    if not (scan_group_overhead <= 1.05):
        guard_errs += 1
    # Bit-exact init_a flat copy + read-only scan8 ORR check
    scan8_errs = 0
    n_cells_t = 256
    init_a_t = [0x8000 - rng.randint(2000, 8000) for _ in range(n_cells_t)]
    rows_8 = [init_a_t[:] for _ in range(8)]
    expected_surv = set()
    for r in range(8):
        for c in range(n_cells_t):
            if rng.random() < 0.03:
                rows_8[r][c] = (rows_8[r][c] + 9000) & 0xFFFF
            if rows_8[r][c] & 0x8000:
                expected_surv.add((r, c))
    found_surv = set()
    for c0 in range(0, n_cells_t, 4):
        orr_acc = 0
        for r in range(8):
            w = rows_8[r][c0] | (rows_8[r][c0 + 1] << 16) | (rows_8[r][c0 + 2] << 32) | (rows_8[r][c0 + 3] << 48)
            orr_acc |= w
        if orr_acc & 0x8000800080008000:
            for r in range(8):
                for dc in range(4):
                    if rows_8[r][c0 + dc] & 0x8000:
                        found_surv.add((r, c0 + dc))
    if found_surv != expected_surv:
        scan8_errs += 1

    ok = (
        b1_cov_ok and b1_d32_cost_ok and d32_errs == 0
        and chord_errs == 0 and chord_order_errs == 0 and chord_tested == 3072
        and s35_alg_ok and cone3_errs == 0 and cone3_tested >= 300
        and fd_rem0_errs == 0 and fd_rem0_max_r2_bits < 25.0
        and tbnz_errs == 0 and sh17_x4_errs == 0
        and sb_nest_errs == 0 and guard_errs == 0 and dummy_sink_ok and scan8_errs == 0
    )
    return dict(
        ok=ok,
        b1_cov_k64_exact=round(cov_64, 5),
        b1_cov_k512_exact=round(cov_512, 5),
        b1_pref_cov_budget=M.B1_PREF_COV,
        b1_d32_p_cut=p_d32,
        b1_d32_max_qb=round(max_qb, 6),
        b1_d32_avg_qb=round(avg_qb, 7),
        b1_d32_cost_ir=round(c_d32, 4),
        b1_c_hit_lo_exact=round(c_hit_lo_exact, 4),
        b1_c_hit_lo_budget=M.C_HIT_LO,
        b1_delta32_apply_errs=d32_errs,
        a0_chord_kta_tested=chord_tested,
        a0_chord_kta_errs=chord_errs,
        a0_chord_ai_order_errs=chord_order_errs,
        a0_chord_mean_points=round(chord_mean_pts, 2),
        a0_chord_mean_pieces=round(pieces_mean, 4),
        a0_chord_piece_mean_points=round(piece_mean_pts, 2),
        cone3_s35_alg_ok=s35_alg_ok,
        cone3_bijection_tested=cone3_tested,
        cone3_bijection_errs=cone3_errs,
        fd_rem0_max_r2_bits=round(fd_rem0_max_r2_bits, 2),
        fd_rem0_u40_errs=fd_rem0_errs,
        tbnz14_single_slot_errs=tbnz_errs,
        sh17_x4_var_ks_errs=sh17_x4_errs,
        b0_floor_div_nest_errs=sb_nest_errs,
        over_radius_eps_bound=round(eps_or, 7),
        fd_bounds=fdb, exact_levels=exl,
        scan_group_overhead=round(scan_group_overhead, 4),
        stor_bound=sbo,
        fd_A1_min_log2=round(log2(fd_A1_min), 2),
        row_capacity_guard_errs=guard_errs,
        dummy_sink_blocks_c0=1 << 15,
        dummy_sink_disjoint_ok=dummy_sink_ok,
        init_scan8_readonly_errs=scan8_errs,
    )
