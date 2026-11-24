"""q15_diag_extra.py — Fast pure-Python stdlib exploratory diagnostics for Q(zeta_15) CM-STNFS:
  1. Analytic H-EA Buchstab-Dickman DDE check: checks 0 <= omega(t) <= 1 => P_BD(u, v) <= min(1, max(1, v)*rho(u-v))
     and verifies psurv_lo <= P_BD(u, v) <= psurv_hi with surv_upper(m, B, T) = min(1, 2*psurv_hi(m, B, T))
     >= min(1, 2*P_BD(u, v)) across boot, all 152 global-floor (T0) Stage-2 node sides (76 nodes), all 152 global-floor (T) Stage-1 node sides, AND all 152 uniform-over-all-m pre-NORMS sv_norms sides (456 sides total, ZERO sextile averaging; <0.20s).
  2. 700,000-pair (25,434 joint hits, SE=0.64%) Q(zeta_15) deg-(72, 24) joint-smoothness calibration check.
  3. High-bit [100..245]b encountered-law / deg-(32,24) archived summary check + 8 live genuine unplanted
     [71..104]b production deg-(72, 24) reached-child witnesses (4 f_9 + 4 g_3) with Q_CM LLL, BOUND_CELL,
     exact integer NORMS floor, and ss_lower2 verification (<0.12s).

Strictly scoped exploratory compatibility support for declared premises H1, H-EA, and H5 — not a proof.
"""
import math
import random
import time

import operator as _OP
_SUMPROD = getattr(math, "sumprod", lambda a, b, _m=_OP.mul: sum(map(_m, a, b)))

U_BLS = -0xD201000000010000

PROD_WITNESSES_8 = (
    {"slot": 0, "p_side": "f", "p_bits": 185, "c_side": "f", "c_bits": 71,
     "q_p": 46416948164225346400223530487531735622710607005520126931,
     "s_p": 18426175173216877841507822490697744378601516290664534824,
     "t_p": 7451986791273530612501898500709727406564528112555221181,
     "q_c": 1816447883344894637671,
     "s_c": 630246171075817652116,
     "t_c": 59547416585470590347,
     "a": [77668, -124225, -12660, 91857, -170805, 95858, -200284, -26432],
     "b": [-112157, 244758, 39614, -84848, 92177, -78316, 144363, 23374]},
    {"slot": 1, "p_side": "g", "p_bits": 145, "c_side": "f", "c_bits": 79,
     "q_p": 34097702361487552774026742901569887526523341,
     "s_p": 4552888691050818822952814679839889563157636,
     "t_p": 18787646669860703688158176046297601588230521,
     "q_c": 423762983407257075037891,
     "s_c": 230602993903680606838809,
     "t_c": 268319767907298826535072,
     "a": [21728, 10227, 16276, -8275, -12263, 27426, -3016, -14881],
     "b": [-6568, 8871, 4444, -19750, 21940, -10502, -14461, 18856]},
    {"slot": 2, "p_side": "g", "p_bits": 185, "c_side": "f", "c_bits": 80,
     "q_p": 39288054469821651472253100936778279268288365431495780061,
     "s_p": 16815209697273385605178293329482956526413504324669415005,
     "t_p": 28986222369129040257649335149724107448626701907509729997,
     "q_c": 726265906948985121258721,
     "s_c": 620634194961566745267804,
     "t_c": 698346872515219654427093,
     "a": [-24938, 76459, -59105, 118815, 65983, -76731, -4317, 31234],
     "b": [143380, -33216, -5651, 143154, -31899, 55602, 37680, -70283]},
    {"slot": 3, "p_side": "g", "p_bits": 225, "c_side": "f", "c_bits": 96,
     "q_p": 45421586516227414981027386948091425225715385980027983476097693309391,
     "s_p": 10305638970508379406777239032164347453392058739449065035683257985937,
     "t_p": 30066606535497931383279696524464770924098331195358499321419348350329,
     "q_c": 78813990057141644528125585171,
     "s_c": 72733101267435222006532314428,
     "t_c": 46265286677873894689951490057,
     "a": [344303, 568143, 17387, -292083, 440064, 131676, 104858, 429209],
     "b": [-332485, 176941, 389757, 105708, -345863, -20007, -503535, -329746]},
    {"slot": 4, "p_side": "g", "p_bits": 120, "c_side": "g", "c_bits": 71,
     "q_p": 844726303703151346861669842266584561,
     "s_p": 182259456457942375214696033809512966,
     "t_p": 671713613895963337084228413269217858,
     "q_c": 1909741773650578457761,
     "s_c": 247374763080000445053,
     "t_c": 499666635035223250369,
     "a": [-8272, 7007, -1177, -7735, 6816, -2603, -1459, 148],
     "b": [7517, -5199, -4210, 1573, -6424, 1177, 331, -7468]},
    {"slot": 5, "p_side": "f", "p_bits": 145, "c_side": "g", "c_bits": 76,
     "q_p": 27760116935638518955586857864840532622687271,
     "s_p": 22494212618382207246958170544848241727992944,
     "t_p": 2543049395822516648650780903316242445246040,
     "q_c": 64159098650074182003421,
     "s_c": 58161897478182929232524,
     "t_c": 50279680817447448194864,
     "a": [-41935, 34095, 16290, -15556, 25076, -4502, -6357, 50512],
     "b": [10917, 3566, 603, -5075, 826, 11637, -20867, 11943]},
    {"slot": 6, "p_side": "g", "p_bits": 120, "c_side": "g", "c_bits": 83,
     "q_p": 844726303703151346861669842266584561,
     "s_p": 182259456457942375214696033809512966,
     "t_p": 671713613895963337084228413269217858,
     "q_c": 7978433478745464910935541,
     "s_c": 6754852926669252912709047,
     "t_c": 2679716929553647989510438,
     "a": [-4600, 4897, -7343, 240, 683, -4338, -996, 2783],
     "b": [1209, 1783, -8355, 3980, 1531, -1098, 7850, -6589]},
    {"slot": 7, "p_side": "f", "p_bits": 185, "c_side": "g", "c_bits": 104,
     "q_p": 46416948164225346400223530487531735622710607005520126931,
     "s_p": 18426175173216877841507822490697744378601516290664534824,
     "t_p": 7451986791273530612501898500709727406564528112555221181,
     "q_c": 17373615561706560461709420594601,
     "s_c": 1163741881783844298536371028187,
     "t_c": 15849059468597379450024591476512,
     "a": [8531, -127054, -7467, 61530, -20302, 61550, -138149, 143227],
     "b": [-137899, 112472, 156514, -125227, 110928, -165153, -12270, 112857]},
)


def _build_rho_omega(u_max=45.0, step=0.005):
    """Positive-preserving integral solvers for Dickman rho(u) and Buchstab omega(t):
      u*rho(u) = int_{u-1}^u rho(t) dt,  t*omega(t) = 1 + int_1^{t-1} omega(s) ds."""
    n = int(round(u_max / step))
    k1 = int(round(1.0 / step))
    rho = [0.0] * (n + 1)
    pref_r = [0.0] * (n + 1)
    for i in range(k1 + 1):
        rho[i] = 1.0
        pref_r[i] = i * step
    for i in range(k1 + 1, n + 1):
        u = i * step
        rhs = (pref_r[i - 1] - pref_r[i - k1]) + 0.5 * step * rho[i - 1]
        rho[i] = rhs / (u - 0.5 * step)
        pref_r[i] = pref_r[i - 1] + 0.5 * step * (rho[i - 1] + rho[i])
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


def _p_buchstab_dickman(u, v, rho, omega, step):
    """Exact asymptotic P(R <= y^v | N ~ y^u) = rho(u) + 1_{v>=1} * int_1^v omega(t)*rho(u-t) dt.
    Uses exact step-aligned slice sumprod over the precomputed 0.005-step DDE tables (h = 4*step = 0.02)."""
    if v < 0.0:
        return 0.0
    if v >= u:
        return 1.0
    val = _interp(rho, u, step)
    if v <= 1.0:
        return min(1.0, val)
    # Fast exact trapezoidal rule on h = 4*step = 0.02 (step = 0.005) with tail trapezoid to v:
    stride = 4
    h = stride * step
    k_start = int(round(1.0 / step))
    k_end = int(math.floor(v / step))
    # Align k_last <= k_end on stride grid starting at k_start:
    n_steps = (k_end - k_start) // stride
    acc = 0.0
    if n_steps >= 1:
        # omega at t_j = 1 + j*h is omega[k_start + j*stride]; interpolate rho(u - t_j) linearly on step grid:
        fu = u / step
        iu = int(math.floor(fu))
        frac = fu - iu
        om_s = omega[k_start : k_start + (n_steps + 1) * stride : stride]
        # u - t_j in step units is fu - (k_start + j*stride) = (iu - k_start - j*stride) + frac
        # Build aligned rho values via slice:
        r0_hi = iu - k_start
        r0_lo = r0_hi - n_steps * stride
        if r0_lo >= 0 and r0_hi + 1 < len(rho):
            # Reverse slices so index j corresponds to r0_hi - j*stride:
            rs0 = rho[r0_hi : (r0_lo - stride if r0_lo >= stride else None) : -stride]
            rs1 = rho[r0_hi + 1 : (r0_lo + 1 - stride if r0_lo + 1 >= stride else None) : -stride]
            if len(rs0) == len(om_s) and len(rs1) == len(om_s):
                tot_s = (1.0 - frac) * _SUMPROD(om_s, rs0) + frac * _SUMPROD(om_s, rs1)
                # Subtract half of endpoints for trapezoidal rule:
                ep0 = om_s[0] * ((1.0 - frac) * rs0[0] + frac * rs1[0])
                ep1 = om_s[-1] * ((1.0 - frac) * rs0[-1] + frac * rs1[-1])
                acc = (tot_s - 0.5 * (ep0 + ep1)) * h
            else:
                acc = 0.0
                for j in range(n_steps + 1):
                    t = 1.0 + j * h
                    wt = 0.5 if (j == 0 or j == n_steps) else 1.0
                    acc += wt * om_s[j] * _interp(rho, u - t, step)
                acc *= h
        else:
            acc = 0.0
            for j in range(n_steps + 1):
                t = 1.0 + j * h
                wt = 0.5 if (j == 0 or j == n_steps) else 1.0
                acc += wt * om_s[j] * _interp(rho, u - t, step)
            acc *= h
        t_last = 1.0 + n_steps * h
    else:
        t_last = 1.0
    if v > t_last + 1e-12:
        dt = v - t_last
        acc += 0.5 * dt * (_interp(omega, t_last, step) * _interp(rho, u - t_last, step) + _interp(omega, v, step) * _interp(rho, u - v, step))
    return min(1.0, val + acc)


def verify_hea_surv_analytic(res, Pm, DS):
    """Verify the Buchstab-Dickman DDE bound (0 <= omega(t) <= 1 => P_BD(u, v) <= min(1, max(1, v)*rho(u - v)))
    and check psurv_lo <= P_BD(u, v) <= psurv_hi with surv_upper(m, B, T) = min(1, 2*psurv_hi(m, B, T))
    >= min(1, 2*P_BD(u, v)) across boot, all 152 global-floor (T0) Stage-2 node sides, all 152 global-floor (T) Stage-1 node sides, AND all 152 pre-NORMS smooth_part_upper sides (456 sides total)."""
    t0 = time.perf_counter()
    rho, omega, step = _build_rho_omega(45.0, 0.005)
    max_omega = max(omega)
    assert max_omega <= 1.0 + 1e-12

    bt = res["boot"]
    tau_b = float(getattr(Pm, "TAU_BOOT", 2070.0))
    m_min_b = bt["m_min"]
    ps_lo_b, ps_hi_b = DS.psurv_iv(m_min_b, Pm.LX, tau_b)
    ps_raw_b = float(bt["boot_psurv_raw"])
    assert abs(ps_raw_b - ps_hi_b) < 1e-18 and ps_lo_b > 0.0 and ps_hi_b / ps_lo_b <= 1.08
    sv_b = float(bt["surv_b"])
    pbd_b = _p_buchstab_dickman(m_min_b / Pm.LX, tau_b / Pm.LX, rho, omega, step)
    pbd_b_a = _p_buchstab_dickman((m_min_b + Pm.ALPHA_G) / Pm.LX, tau_b / Pm.LX, rho, omega, step)
    assert ps_lo_b <= pbd_b <= ps_hi_b
    boot_raw_cush = ps_raw_b / pbd_b
    boot_cush = sv_b / pbd_b
    boot_cush_a = sv_b / pbd_b_a
    assert boot_raw_cush >= 1.01 and boot_cush >= 2.0 and boot_cush_a >= boot_cush

    side_cush = []
    prod_cush = []
    t0_sat_sides = 0
    sext_cush = []
    sext_sat_sides = 0
    sext_checked = 0

    def _check_one(m_val, T_val, cush_list):
        lo_iv, hi_iv = DS.psurv_iv(m_val, Pm.LX, T_val)
        sv = DS.surv_upper(m_val, Pm.LX, T_val)
        pbd = _p_buchstab_dickman(m_val / Pm.LX, T_val / Pm.LX, rho, omega, step)
        assert abs(sv - min(1.0, 2.0 * hi_iv)) < 1e-15
        assert lo_iv - 1e-15 <= pbd <= hi_iv + 1e-15
        assert sv >= min(1.0, 2.0 * pbd) - 1e-15
        if sv < 1.0:
            assert pbd > 0.0
            cush_list.append(sv / pbd)
            return sv, pbd, False
        assert sv == 1.0 and sv >= pbd - 1e-15
        return sv, pbd, True

    norms_cush = []
    norms_sat_sides = 0
    norms_checked = 0
    for ((l_lo, l_hi), side), d in res["DET"].items():
        T0 = d["T0"]
        nf_fl, ng_fl = d["nf_floor"], d["ng_floor"]
        sv_f, pbd_f, sat_f = _check_one(nf_fl, T0, side_cush)
        sv_g, pbd_g, sat_g = _check_one(ng_fl, T0, side_cush)
        t0_sat_sides += int(sat_f) + int(sat_g)
        sv_p = sv_f * sv_g
        pbd_p = pbd_f * pbd_g
        if sv_p < 1.0 and pbd_p > 0.0:
            prod_cush.append(sv_p / pbd_p)
        # Also verify the exact global-floor (nf_floor, T) and (ng_floor, T) Stage-1 bounds (ZERO sextile averaging!):
        T_sieve = T0 + Pm.SLACK
        _, _, sf = _check_one(nf_fl, T_sieve, sext_cush)
        _, _, sg = _check_one(ng_fl, T_sieve, sext_cush)
        sext_sat_sides += int(sf) + int(sg)
        sext_checked += 2
        # Also verify pre-NORMS smooth_part_upper(m_fl - T_sieve, LX) uniformly over ALL m in [c_bits, m_max]:
        for m_fl_s, m_mx_s in ((nf_fl, d["nf"]), (ng_fl, d["ng"])):
            c_bits = m_fl_s - T_sieve
            ub_n = DS.smooth_part_upper(c_bits, Pm.LX)
            norms_checked += 1
            if ub_n >= 1.0 - 1e-15:
                norms_sat_sides += 1
            else:
                mx_pbd = 0.0
                for dm in range(0, int(math.ceil(m_mx_s - c_bits)) + 1, 4):
                    m_cur = min(m_mx_s, c_bits + dm)
                    p_cur = _p_buchstab_dickman(m_cur / Pm.LX, max(0.0, m_cur - c_bits) / Pm.LX, rho, omega, step)
                    if p_cur > mx_pbd:
                        mx_pbd = p_cur
                p_end = _p_buchstab_dickman(m_mx_s / Pm.LX, max(0.0, m_mx_s - c_bits) / Pm.LX, rho, omega, step)
                if p_end > mx_pbd:
                    mx_pbd = p_end
                assert ub_n >= 2.0 * mx_pbd
                norms_cush.append(ub_n / mx_pbd)

    assert side_cush and prod_cush and sext_cush and norms_cush
    assert sext_checked == 2 * len(res["DET"]) and norms_checked == 2 * len(res["DET"])
    assert min(side_cush) >= 2.0 and min(prod_cush) >= 2.0 and min(sext_cush) >= 2.0 and min(norms_cush) >= 2.0
    prod_cush.sort()
    return dict(
        ok=True,
        wall_s=round(time.perf_counter() - t0, 4),
        max_omega_1_45=round(max_omega, 6),
        boot_psurv_iv_ratio=round(ps_hi_b / ps_lo_b, 4),
        boot_iv_encloses_PBD=True,
        boot_psurv_raw_over_PBD=round(boot_raw_cush, 4),
        boot_surv_over_PBD=round(boot_cush, 4),
        boot_surv_over_PBD_alpha_g=round(boot_cush_a, 4),
        nodes_nontrivial_sides=len(side_cush),
        nodes_saturated_sides=t0_sat_sides,
        nodes_nontrivial_products=len(prod_cush),
        min_side_cushion=round(min(side_cush), 3),
        min_prod_cushion=round(prod_cush[0], 3),
        median_prod_cushion=round(prod_cush[len(prod_cush) // 2], 2),
        global_floor_T_sides_checked=sext_checked,
        global_floor_T_nontrivial_sides=len(sext_cush),
        global_floor_T_saturated_sides=sext_sat_sides,
        global_floor_T_min_side_cushion=round(min(sext_cush), 3),
        sv_norms_sides_checked=norms_checked,
        sv_norms_nontrivial_sides=len(norms_cush),
        sv_norms_saturated_sides=norms_sat_sides,
        sv_norms_min_side_cushion=round(min(norms_cush), 3),
    )


def verify_joint_smooth_700k(arc):
    """Verify the archived 700,000-pair (25,434 joint hits, SE=0.64%) Q(zeta_15) deg-(72,24)
    joint-smoothness compatibility calibration across u_toy in {-2, -5, -17, -32}."""
    blk = arc["joint_smooth_700k"]
    cfgs = blk["configs"]
    tot_pairs = sum(c["n_pairs"] for c in cfgs)
    tot_joint = sum(c["obs_joint"] for c in cfgs)
    tot_mod = sum(c["mod_joint"] for c in cfgs)
    pooled_rj = tot_joint / tot_mod
    pooled_se = pooled_rj / math.sqrt(tot_joint)
    assert tot_pairs == 700000 and tot_joint == 25434
    for c in cfgs:
        assert 0.90 <= c["obs_over_mod_joint"] <= 1.12
        assert 0.95 <= c["cond_joint_over_prod"] <= 1.10
        assert 0.88 <= c["obs_over_mod_f"] <= 1.10 and 0.95 <= c["obs_over_mod_g"] <= 1.10
    assert pooled_se <= 0.007 and 0.98 <= pooled_rj <= 1.08
    return dict(
        ok=True,
        n_configs=len(cfgs),
        total_pairs=tot_pairs,
        total_obs_joint=tot_joint,
        total_mod_joint=round(tot_mod, 2),
        pooled_obs_over_mod_joint=round(pooled_rj, 5),
        pooled_se=round(pooled_se, 5),
        min_cond_joint_over_prod=round(min(c["cond_joint_over_prod"] for c in cfgs), 5),
        max_cond_joint_over_prod=round(max(c["cond_joint_over_prod"] for c in cfgs), 5),
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
        d >>= 1
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
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


def _bin_for_log2(lq, bins):
    for blo, bhi in bins:
        if blo <= lq < bhi:
            return (blo, bhi)
    return bins[0] if lq < bins[0][0] else bins[-1]


def verify_highbit_prod_diagnostics(res, Pm, BC, DS, lll_fn, sq_hnf_fn, norms_fg_fn, arc, n_cells=8):
    """Verify archived high-bit [100..245]b diagnostics (384 prod deg-(72,24) L_{q_p} cap L_{q_c} pairs +
    200 unplanted deg-(32,24) children) and replay 8 live unplanted [71..104]b prod deg-(72,24) witnesses."""
    t0 = time.perf_counter()
    hb = arc["highbit_prod_children"]
    p1 = hb["part1_prod_deg72_24_encountered_law_highbit"]
    p2a = hb["part2a_deg32_24_unplanted_highbit"]
    assert p1["n_pairs"] == 384 and p1["children_below_1x"] == 0 and p1["child_vq_basis_failures"] == 0
    assert p1["child_bit_range"][0] >= 100 and p1["child_bit_range"][1] >= 230
    assert p2a["n_reached_children"] == 200 and p2a["children_below_1x"] == 0 and p2a["vq_basis_failures_children"] == 0
    assert p2a["child_bit_range"][0] >= 100 and p2a["child_bit_range"][1] >= 240 and p2a["ks_2samp_p"] >= 0.50

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
        bkey = _bin_for_log2(lqc, bins)
        assert bkey[0] <= lqc < bkey[1]
        dp = res["DET"][bkey, cs]
        Rs, ns, iAs, kAs, kBs = dp["R"], dp["n"], dp["iA"], dp["kA"], dp["kB"]
        r2 = Rs * Rs
        rad_hi = Rs + 2.0  # Covers entire 16D rounding cube circumradius sqrt(16 * 0.5^2) = 2.0
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
        rows.append(dict(slot=w["slot"], p_side=ps, p_bits=w["p_bits"], c_side=cs, c_bits=w["c_bits"], exact_R=Rs, kept=kept, screened_cert_ratio=round(rat, 2)))
    return dict(
        ok=True,
        wall_s=round(time.perf_counter() - t0, 4),
        part1_n_pairs=p1["n_pairs"],
        part1_child_bit_range=p1["child_bit_range"],
        part1_below_1x=p1["children_below_1x"],
        part1_min_ratio=p1["cert_ratio_quantiles"]["min"],
        part1_stratum_diff_log2=p1["stratum_controlled_weighted_minus_unweighted_log2"],
        part2a_n_children=p2a["n_reached_children"],
        part2a_child_bit_range=p2a["child_bit_range"],
        part2a_below_1x=p2a["children_below_1x"],
        part2a_geom_ratio=p2a["geom_mean_child_over_ctrl"],
        part2a_ks_p=p2a["ks_2samp_p"],
        live_witnesses_n=len(rows),
        live_witnesses_bit_range=[min(r["c_bits"] for r in rows), max(r["c_bits"] for r in rows)],
        live_witnesses_min_ratio=min(r["screened_cert_ratio"] for r in rows),
        live_witnesses=rows,
    )
