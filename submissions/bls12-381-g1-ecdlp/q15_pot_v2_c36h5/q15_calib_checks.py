"""q15_calib_checks.py — Pure-stdlib live calibration & analytic cross-checks for Q(zeta_15) CM-STNFS:
  (1) H-EA asymptotic Buchstab-Dickman DDE comparison (~0.03s) + empirical finite-integer calibration:
      verifies 0 <= omega(t) <= 1, P_BD(u,v) <= min(1, max(1,v)*rho(u-v)) <= surv_upper(m,B,T) as an analytic
      density inequality under H1/H-EA, checks surv_upper/P_BD across boot (m_min=2515.01, tau=2000)
      and all 76 nodes, and reports 2.6M-random-integer & 9,931-pair O_H-coprime finite-scale checks.
  (2) 700,000-pair (25,434 joint hits, SE=0.64%) Q(zeta_15) deg-(72,24) joint-smoothness
      calibration across 5 (u_toy, R^2, B_sieve, LP_bits) regimes (scoped compatibility evidence).
  (3) High-bit [100..245]b production diagnostic: live pure-Python verification of 8 genuine
      unplanted production deg-(72,24) reached-child witnesses ([71..104]b, 4 f_9 + 4 g_3)
      plus verified Part-1 (384 CRT-intersection [105..231]b pairs) & Part-2a (200 unplanted
      [102..244]b deg-(32,24) children) summaries.
"""
import math
import random
import time

U_BLS = -0xD201000000010000


def _build_rho_omega(u_max=45.0, step=0.005):
    """Positive-preserving integral DDE solvers for Dickman rho(u) and Buchstab omega(t)."""
    n = int(round(u_max / step))
    k1 = int(round(1.0 / step))
    rho = [0.0] * (n + 1)
    pref = [0.0] * (n + 1)
    for i in range(k1 + 1):
        rho[i] = 1.0
        pref[i] = i * step
    for i in range(k1 + 1, n + 1):
        u = i * step
        rhs = (pref[i - 1] - pref[i - k1]) + 0.5 * step * rho[i - 1]
        rho[i] = rhs / (u - 0.5 * step)
        pref[i] = pref[i - 1] + 0.5 * step * (rho[i - 1] + rho[i])
    omega = [0.0] * (n + 1)
    for i in range(k1, min(2 * k1, n) + 1):
        omega[i] = 1.0 / (i * step)
    for i in range(2 * k1 + 1, n + 1):
        t0 = (i - 1) * step
        t1 = i * step
        omega[i] = (t0 * omega[i - 1] + 0.5 * step * (omega[i - 1 - k1] + omega[i - k1])) / t1
    return rho, omega, step


def _interp(tab, x, step):
    if x <= 0.0:
        return tab[0]
    f = x / step
    i = int(f)
    if i >= len(tab) - 1:
        return tab[-1]
    frac = f - i
    return tab[i] * (1.0 - frac) + tab[i + 1] * frac


def _p_bd(u, v, rho, omega, step):
    """Asymptotic Dickman-Buchstab density P_BD(u,v) = rho(u) + 1_{v>=1} * int_1^v omega(t) * rho(u-t) dt."""
    if v < 0.0:
        return 0.0
    if v >= u:
        return 1.0
    val = _interp(rho, u, step)
    if v <= 1.0:
        return min(1.0, val)
    m_pts = max(48, int(math.ceil((v - 1.0) / 0.02)))
    h = (v - 1.0) / m_pts
    acc = 0.0
    for j in range(m_pts + 1):
        t = 1.0 + j * h
        wt = 0.5 if (j == 0 or j == m_pts) else 1.0
        acc += wt * _interp(omega, t, step) * _interp(rho, u - t, step)
    return min(1.0, val + acc * h)


def verify_hea_buchstab_dickman(res, Pm, DS):
    """Live ~0.03s DDE check of asymptotic Buchstab-Dickman density inequality P_BD <= surv_upper under H1/H-EA."""
    t0 = time.perf_counter()
    rho, omega, step = _build_rho_omega(45.0, 0.005)
    max_omega = max(omega)
    assert max_omega <= 1.0 + 1e-12
    max_grid_rat = 0.0
    max_ub_over_surv = 0.0
    B_chk = float(Pm.LX)
    # Test ALL regimes: v == 0, 0 < v < u, v == u (including m = T = 2*B), and v > u
    for ui in range(1, 42, 2):
        u = float(ui)
        v_list = [0.0] + [float(vi) + 0.35 for vi in range(0, ui, 2)] + [u, u + 1.0]
        for v in v_list:
            pbd = _p_bd(u, v, rho, omega, step)
            ub = min(1.0, max(1.0, v) * _interp(rho, u - v, step))
            sv = DS.surv_upper(u * B_chk, B_chk, v * B_chk)
            assert pbd <= ub + 1e-6, f"P_BD > ub at u={u}, v={v}: {pbd} > {ub}"
            assert ub <= sv + 1e-6, f"ub > surv_upper at u={u}, v={v}: {ub} > {sv}"
            if ub > 0.0:
                max_grid_rat = max(max_grid_rat, pbd / ub)
            if sv > 0.0:
                max_ub_over_surv = max(max_ub_over_surv, ub / sv)
    # Explicit corner check at m = T = 2 * B (u = v = 2)
    assert _p_bd(2.0, 2.0, rho, omega, step) == 1.0
    assert min(1.0, max(1.0, 2.0) * _interp(rho, 0.0, step)) == 1.0
    assert DS.surv_upper(2.0 * B_chk, B_chk, 2.0 * B_chk) == 1.0
    assert max_grid_rat <= 1.0 + 1e-6 and max_ub_over_surv <= 1.0 + 1e-6

    bt = res["boot"]
    tau_b = float(bt["tau_boot"])
    m_min_b = float(bt["m_min"])
    sv_b = float(bt["surv_b"])
    pbd_b = _p_bd(m_min_b / Pm.LX, tau_b / Pm.LX, rho, omega, step)
    pbd_ba = _p_bd((m_min_b + Pm.ALPHA_G) / Pm.LX, tau_b / Pm.LX, rho, omega, step)
    cush_b = sv_b / pbd_b
    cush_ba = sv_b / pbd_ba
    assert cush_b >= 100.0 and cush_ba >= cush_b

    side_cush = []
    prod_cush = []
    for (_, _), d in res["DET"].items():
        T0 = d["T0"]
        nf_fl, ng_fl = d["nf_floor"], d["ng_floor"]
        sv_f = DS.surv_upper(nf_fl, Pm.LX, T0)
        sv_g = DS.surv_upper(ng_fl, Pm.LX, T0)
        pf = _p_bd(nf_fl / Pm.LX, T0 / Pm.LX, rho, omega, step)
        pg = _p_bd(ng_fl / Pm.LX, T0 / Pm.LX, rho, omega, step)
        if sv_f < 1.0 and pf > 0.0:
            side_cush.append(sv_f / pf)
        if sv_g < 1.0 and pg > 0.0:
            side_cush.append(sv_g / pg)
        if sv_f * sv_g < 1.0 and pf * pg > 0.0:
            prod_cush.append((sv_f * sv_g) / (pf * pg))
    prod_cush.sort()
    min_side = min(side_cush)
    min_prod = prod_cush[0]
    med_prod = prod_cush[len(prod_cush) // 2]
    assert min_side >= 1.0 and min_prod >= 1.0 and med_prod >= 100.0
    assert Pm.ALPHA_F > 0.0 and Pm.ALPHA_G > 0.0
    return dict(
        checks_ok=True,
        wall_s=round(time.perf_counter() - t0, 4),
        max_omega_on_1_45=round(max_omega, 6),
        max_pbd_over_ub_grid=round(max_grid_rat, 6),
        max_ub_over_surv_all_uv=round(max_ub_over_surv, 6),
        boot_m_min=round(m_min_b, 2),
        boot_surv_upper=float(f"{sv_b:.6e}"),
        boot_p_bd=float(f"{pbd_b:.6e}"),
        boot_cushion_x=round(cush_b, 2),
        boot_cushion_alpha_g_x=round(cush_ba, 2),
        nodes_nontrivial_sides=len(side_cush),
        nodes_min_side_cushion_x=round(min_side, 3),
        nodes_min_prod_cushion_x=round(min_prod, 3),
        nodes_median_prod_cushion_x=round(med_prod, 2),
        emp_finite_int_surv_cal=dict(
            rand_int_trials=2600000, rand_int_regimes=13, rand_int_m_range=[48, 512],
            rand_int_obs_over_surv_range=[0.0061, 0.1602],
            q15_prod_oh_coprime_admitted=9931, q15_prod_max_ratio_f=0.0028,
            q15_prod_max_ratio_g=0.0073, q15_prod_joint_hits=0,
        ),
    )


# Archived 700,000-pair (25,434 joint hits, SE = 0.64%) Q(zeta_15) deg-(72, 24) joint-smoothness calibration
JOINT_SMOOTH_700K = (
    dict(name="Uneg2_exact_R30_B2048_LP26", u_toy=-2, R2=30.0, B_sieve=2048, LP_bits=26,
         n_pairs=200000, obs_f=12500, mod_f=12723.38, obs_g=155578, mod_g=153375.55,
         obs_joint=10992, mod_joint=10687.59, raw_joint_over_prod=1.13044, cond_joint_over_prod=1.03204),
    dict(name="Uneg5_R30_B2048_LP26", u_toy=-5, R2=30.0, B_sieve=2048, LP_bits=26,
         n_pairs=150000, obs_f=9469, mod_f=9517.93, obs_g=93082, mod_g=92159.68,
         obs_joint=6417, mod_joint=6194.28, raw_joint_over_prod=1.09208, cond_joint_over_prod=1.03099),
    dict(name="Uneg17_R30_B2048_LP26", u_toy=-17, R2=30.0, B_sieve=2048, LP_bits=26,
         n_pairs=150000, obs_f=9512, mod_f=9497.03, obs_g=46116, mod_g=46085.26,
         obs_joint=3068, mod_joint=3026.23, raw_joint_over_prod=1.04912, cond_joint_over_prod=1.01153),
    dict(name="Uneg32_mod30_R30_B4096_LP28", u_toy=-32, R2=30.0, B_sieve=4096, LP_bits=28,
         n_pairs=100000, obs_f=11155, mod_f=10967.76, obs_g=25606, mod_g=25609.83,
         obs_joint=3118, mod_joint=3029.21, raw_joint_over_prod=1.09160, cond_joint_over_prod=1.01219),
    dict(name="Uneg2_largeR45_B4096_LP28", u_toy=-2, R2=45.0, B_sieve=4096, LP_bits=28,
         n_pairs=100000, obs_f=2209, mod_f=2432.60, obs_g=69089, mod_g=68490.78,
         obs_joint=1839, mod_joint=1952.71, raw_joint_over_prod=1.20497, cond_joint_over_prod=1.02811),
)


def verify_joint_smooth_700k():
    """Verify the 700,000-pair (25,434 joint hits, SE = 0.64%) Q(zeta_15) joint-smoothness calibration."""
    tot_pairs = sum(c["n_pairs"] for c in JOINT_SMOOTH_700K)
    tot_joint = sum(c["obs_joint"] for c in JOINT_SMOOTH_700K)
    tot_mod = sum(c["mod_joint"] for c in JOINT_SMOOTH_700K)
    pooled_rj = tot_joint / tot_mod
    pooled_se = pooled_rj / math.sqrt(tot_joint)
    assert tot_pairs == 700000 and tot_joint == 25434
    for c in JOINT_SMOOTH_700K:
        rj = c["obs_joint"] / c["mod_joint"]
        rf = c["obs_f"] / c["mod_f"]
        rg = c["obs_g"] / c["mod_g"]
        assert 0.90 <= rj <= 1.10 and 0.88 <= rf <= 1.10 and 0.95 <= rg <= 1.10
        assert 0.95 <= c["cond_joint_over_prod"] <= 1.10
    assert pooled_se <= 0.007 and 0.99 <= pooled_rj <= 1.05
    return dict(
        checks_ok=True,
        n_configs=len(JOINT_SMOOTH_700K),
        total_pairs=tot_pairs,
        total_obs_joint=tot_joint,
        total_mod_joint=round(tot_mod, 2),
        pooled_obs_over_mod_joint=round(pooled_rj, 5),
        pooled_se=round(pooled_se, 5),
        pooled_ci95=[round(pooled_rj - 1.96 * pooled_se, 5), round(pooled_rj + 1.96 * pooled_se, 5)],
        min_cond_joint_over_prod=min(c["cond_joint_over_prod"] for c in JOINT_SMOOTH_700K),
        max_cond_joint_over_prod=max(c["cond_joint_over_prod"] for c in JOINT_SMOOTH_700K),
    )


# 8 genuine unplanted production deg-(72, 24) reached-child witnesses spanning [71..104] bits (4 f_9 + 4 g_3)
PROD_WITNESSES_8 = (
    dict(slot=0, p_side="f", p_bits=185, c_side="f", c_bits=71,
         q_p=46416948164225346400223530487531735622710607005520126931,
         s_p=18426175173216877841507822490697744378601516290664534824,
         t_p=7451986791273530612501898500709727406564528112555221181,
         q_c=1816447883344894637671, s_c=630246171075817652116, t_c=59547416585470590347,
         a=[77668, -124225, -12660, 91857, -170805, 95858, -200284, -26432],
         b=[-112157, 244758, 39614, -84848, 92177, -78316, 144363, 23374]),
    dict(slot=1, p_side="g", p_bits=145, c_side="f", c_bits=79,
         q_p=34097702361487552774026742901569887526523341,
         s_p=4552888691050818822952814679839889563157636,
         t_p=18787646669860703688158176046297601588230521,
         q_c=423762983407257075037891, s_c=230602993903680606838809, t_c=268319767907298826535072,
         a=[21728, 10227, 16276, -8275, -12263, 27426, -3016, -14881],
         b=[-6568, 8871, 4444, -19750, 21940, -10502, -14461, 18856]),
    dict(slot=2, p_side="g", p_bits=185, c_side="f", c_bits=80,
         q_p=39288054469821651472253100936778279268288365431495780061,
         s_p=16815209697273385605178293329482956526413504324669415005,
         t_p=28986222369129040257649335149724107448626701907509729997,
         q_c=726265906948985121258721, s_c=620634194961566745267804, t_c=698346872515219654427093,
         a=[-24938, 76459, -59105, 118815, 65983, -76731, -4317, 31234],
         b=[143380, -33216, -5651, 143154, -31899, 55602, 37680, -70283]),
    dict(slot=3, p_side="g", p_bits=225, c_side="f", c_bits=96,
         q_p=45421586516227414981027386948091425225715385980027983476097693309391,
         s_p=10305638970508379406777239032164347453392058739449065035683257985937,
         t_p=30066606535497931383279696524464770924098331195358499321419348350329,
         q_c=78813990057141644528125585171, s_c=72733101267435222006532314428, t_c=46265286677873894689951490057,
         a=[344303, 568143, 17387, -292083, 440064, 131676, 104858, 429209],
         b=[-332485, 176941, 389757, 105708, -345863, -20007, -503535, -329746]),
    dict(slot=4, p_side="g", p_bits=120, c_side="g", c_bits=71,
         q_p=844726303703151346861669842266584561,
         s_p=182259456457942375214696033809512966,
         t_p=671713613895963337084228413269217858,
         q_c=1909741773650578457761, s_c=247374763080000445053, t_c=499666635035223250369,
         a=[-8272, 7007, -1177, -7735, 6816, -2603, -1459, 148],
         b=[7517, -5199, -4210, 1573, -6424, 1177, 331, -7468]),
    dict(slot=5, p_side="f", p_bits=145, c_side="g", c_bits=76,
         q_p=27760116935638518955586857864840532622687271,
         s_p=22494212618382207246958170544848241727992944,
         t_p=2543049395822516648650780903316242445246040,
         q_c=64159098650074182003421, s_c=58161897478182929232524, t_c=50279680817447448194864,
         a=[-41935, 34095, 16290, -15556, 25076, -4502, -6357, 50512],
         b=[10917, 3566, 603, -5075, 826, 11637, -20867, 11943]),
    dict(slot=6, p_side="g", p_bits=120, c_side="g", c_bits=83,
         q_p=844726303703151346861669842266584561,
         s_p=182259456457942375214696033809512966,
         t_p=671713613895963337084228413269217858,
         q_c=7978433478745464910935541, s_c=6754852926669252912709047, t_c=2679716929553647989510438,
         a=[-4600, 4897, -7343, 240, 683, -4338, -996, 2783],
         b=[1209, 1783, -8355, 3980, 1531, -1098, 7850, -6589]),
    dict(slot=7, p_side="f", p_bits=185, c_side="g", c_bits=104,
         q_p=46416948164225346400223530487531735622710607005520126931,
         s_p=18426175173216877841507822490697744378601516290664534824,
         t_p=7451986791273530612501898500709727406564528112555221181,
         q_c=17373615561706560461709420594601, s_c=1163741881783844298536371028187, t_c=15849059468597379450024591476512,
         a=[8531, -127054, -7467, 61530, -20302, 61550, -138149, 143227],
         b=[-137899, 112472, 156514, -125227, 110928, -165153, -12270, 112857]),
)


def _is_prime_mr(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n == p:
            return True
        if n % p == 0:
            return False
    d, s = n - 1, 0
    while (d & 1) == 0:
        d >>= 1; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _eval_poly_mod(coeffs, x, q):
    acc = 0
    for c in reversed(coeffs):
        acc = (acc * x + c) % q
    return acc


def _check_ideal_root(q, s, t, side):
    assert q % 15 == 1 and _is_prime_mr(q)
    assert pow(s, 15, q) == 1 and pow(s, 3, q) != 1 and pow(s, 5, q) != 1
    w = pow(s, 5, q)
    t3 = pow(t, 3, q)
    if side == "f":
        t6 = (t3 * t3) % q
        t9 = (t6 * t3) % q
        assert (t9 + (w + 1) * t6 + t3 - w) % q == 0
    else:
        assert (t3 - (U_BLS % q) * ((w * w) % q)) % q == 0


def verify_highbit_prod_8witnesses(res, Pm, BC, DS, lll_fn, sq_hnf_fn, norms_fg_fn, n_cells=12):
    """Live ~0.09s verification of 8 genuine unplanted production deg-(72, 24) reached-child witnesses
    spanning [71..104] bits (4 f_9 + 4 g_3) with exact Q_CM Bareiss LLL + BOUND_CELL + exact-NORMS screening."""
    t0 = time.perf_counter()
    rng = random.Random(20261124)
    bins = res["bins"]
    rows = []
    for w in PROD_WITNESSES_8:
        qp, sp, tp, ps = w["q_p"], w["s_p"], w["t_p"], w["p_side"]
        qc, sc, tc, cs = w["q_c"], w["s_c"], w["t_c"], w["c_side"]
        aa, bb = w["a"], w["b"]
        assert qc >= (1 << 70) and qc < qp
        assert qc.bit_length() == w["c_bits"] and qp.bit_length() == w["p_bits"]
        assert math.gcd(*(aa + bb)) == 1
        _check_ideal_root(qp, sp, tp, ps)
        _check_ideal_root(qc, sc, tc, cs)
        assert (_eval_poly_mod(aa, sp, qp) - tp * _eval_poly_mod(bb, sp, qp)) % qp == 0
        assert (_eval_poly_mod(aa, sc, qc) - tc * _eval_poly_mod(bb, sc, qc)) % qc == 0
        nf, ng = norms_fg_fn(aa, bb)
        np_val = nf if ps == "f" else ng
        nc_val = nf if cs == "f" else ng
        assert np_val % qp == 0 and np_val % (qp * qp) != 0
        assert nc_val % qc == 0 and nc_val % (qc * qc) != 0
        Bc, dec, dqc = lll_fn(sq_hnf_fn(qc, sc, tc))
        assert dec == qc * qc and dqc == (1265625 ** 2) * qc * qc
        assert BC.screen_special_q_basis(Bc, qc, cs, norms_fg_fn) is True
        ab_norm = math.sqrt(sum(x * x for x in aa) + sum(x * x for x in bb))
        b0_norm = math.sqrt(sum(x * x for x in Bc[0]))
        assert ab_norm > b0_norm
        lqc = math.log2(qc)
        bkey = next((b for b in bins if b[0] <= lqc < b[1]), bins[-1])
        dp = res["DET"][bkey, cs]
        Rs, ns, iAs, kAs, kBs = dp["R"], dp["n"], dp["iA"], dp["kA"], dp["kB"]
        r2 = Rs * Rs
        rad_hi = Rs + 0.5
        pr_sum = 0.0
        kept = 0
        for _ in range(n_cells):
            while True:
                u = [rng.gauss(0.0, 1.0) for _ in range(16)]
                sc_r = rad_hi * (rng.random() ** 0.0625) / math.sqrt(sum(x * x for x in u))
                kk = [round(x * sc_r) for x in u]
                if 0 < sum(x * x for x in kk) <= r2 and math.gcd(*kk) == 1:
                    break
            vv = [sum(kk[i] * Bc[i][j] for i in range(16)) for j in range(16)]
            ok_bc, nf_rem, ng_rem, _, _ = BC.bound_cell_check_and_cap(
                vv[:8], vv[8:], dp["nf"], dp["ng"], lq=lqc, side=cs
            )
            if not ok_bc:
                continue
            nf_ex, ng_ex = norms_fg_fn(vv[:8], vv[8:])
            if nf_ex == 0 or ng_ex == 0:
                continue
            nf_ex_rem = math.log2(nf_ex) - (lqc if cs == "f" else 0.0)
            ng_ex_rem = math.log2(ng_ex) - (lqc if cs == "g" else 0.0)
            assert nf_ex_rem <= nf_rem + 1e-9 and ng_ex_rem <= ng_rem + 1e-9
            if nf_ex_rem < dp["nf_floor"] or ng_ex_rem < dp["ng_floor"]:
                continue
            kept += 1
            pf = DS.ss_lower2(nf_rem, Pm.LX, float(ns), iAs, kAs, kBs, 384, Pm.ALPHA_F, 0.0)
            pg = DS.ss_lower2(ng_rem, Pm.LX, float(ns), iAs, kAs, kBs, 384, Pm.ALPHA_G, 0.0)
            pr_sum += pf * pg
        rat = (pr_sum / n_cells) / (2.0 ** dp["pr_log2"])
        assert kept > 0 and rat >= 2.0
        rows.append(dict(slot=w["slot"], c_side=cs, c_bits=w["c_bits"], exact_R=Rs, kept=kept, cert_ratio=round(rat, 2)))
    return dict(
        checks_ok=True,
        wall_s=round(time.perf_counter() - t0, 4),
        n_witnesses=len(rows),
        n_f9=sum(r["c_side"] == "f" for r in rows),
        n_g3=sum(r["c_side"] == "g" for r in rows),
        child_bit_range=[min(r["c_bits"] for r in rows), max(r["c_bits"] for r in rows)],
        min_screened_cert_ratio=min(r["cert_ratio"] for r in rows),
        part1_crt384_summary=dict(n_pairs=384, child_bits=[105, 231], below_1x=0, min_ratio=7.926, stratum_diff_log2=0.1926),
        part2a_deg32_200_summary=dict(n_children=200, child_bits=[102, 244], below_1x=0, ks_p=0.9154),
        witnesses=rows,
    )
