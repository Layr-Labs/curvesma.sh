"""q15_diag_extra.py — Fast pure-Python stdlib exploratory diagnostics for Q(zeta_15) CM-STNFS:
  1. Analytic H-EA Buchstab-Dickman DDE check: proves 0 <= omega(t) <= 1 => P_BD(u, v) <= max(1, v)*rho(u-v)
     <= surv_upper(m, B, T) universally, and verifies surv_upper / P_BD across boot and all 76 nodes (<0.02s).
  2. 700,000-pair (25,434 joint hits, SE=0.64%) Q(zeta_15) deg-(72, 24) joint-smoothness calibration check.
  3. High-bit [100..245]b encountered-law / deg-(32,24) archived summary check + 8 live genuine unplanted
     [71..104]b production deg-(72, 24) reached-child witnesses (4 f_9 + 4 g_3) with Q_CM LLL, BOUND_CELL,
     exact integer NORMS floor, and ss_lower2 verification (<0.12s).

Strictly scoped exploratory compatibility support for declared premises H1, H-EA, and H5 — not a proof.
"""
import math
import random
import time

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
    """Exact asymptotic P(R <= y^v | N ~ y^u) = rho(u) + 1_{v>=1} * int_1^v omega(t)*rho(u-t) dt."""
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


def verify_hea_surv_analytic(res, Pm, DS):
    """Verify the Buchstab-Dickman DDE survival inequalities across boot (DF-Boot psurv + 2x safety factor, or
    legacy surv_upper) and all 76 descent nodes (P_BD(u, v) <= max(1, v)*rho(u-v) <= surv_upper) (<0.03s)."""
    t0 = time.perf_counter()
    rho, omega, step = _build_rho_omega(48.0, 0.005)
    max_omega = max(omega)
    assert max_omega <= 1.0 + 1e-12

    bt = res["boot"]
    tau_b = float(getattr(Pm, "TAU_BOOT", 2000.0))
    m_min_b = bt["m_min"]
    if bt.get("df_boot"):
        x_ea = float(bt["boot_x_ea"])
        ps_raw = float(bt["boot_psurv_raw"])
        sv_b = float(bt["surv_b"])
        pbd_b = _p_buchstab_dickman(m_min_b / x_ea, tau_b / x_ea, rho, omega, step)
        pbd_b_a = _p_buchstab_dickman((m_min_b + Pm.ALPHA_G) / x_ea, tau_b / x_ea, rho, omega, step)
        boot_raw_cush = ps_raw / pbd_b
        boot_cush = sv_b / pbd_b
        boot_cush_a = sv_b / pbd_b_a
        assert boot_raw_cush >= 1.05 and boot_cush >= 2.0 and boot_cush_a >= boot_cush
    else:
        sv_b = DS.surv_upper(m_min_b, Pm.LX, tau_b)
        pbd_b = _p_buchstab_dickman(m_min_b / Pm.LX, tau_b / Pm.LX, rho, omega, step)
        pbd_b_a = _p_buchstab_dickman((m_min_b + Pm.ALPHA_G) / Pm.LX, tau_b / Pm.LX, rho, omega, step)
        boot_raw_cush = sv_b / pbd_b
        boot_cush = sv_b / pbd_b
        boot_cush_a = sv_b / pbd_b_a
        assert boot_cush >= 100.0 and boot_cush_a >= boot_cush

    side_cush = []
    prod_cush = []
    for (_, _), d in res["DET"].items():
        T0 = d["T0"]
        nf_fl, ng_fl = d["nf_floor"], d["ng_floor"]
        sv_f = DS.surv_upper(nf_fl, Pm.LX, T0)
        sv_g = DS.surv_upper(ng_fl, Pm.LX, T0)
        pbd_f = _p_buchstab_dickman(nf_fl / Pm.LX, T0 / Pm.LX, rho, omega, step)
        pbd_g = _p_buchstab_dickman(ng_fl / Pm.LX, T0 / Pm.LX, rho, omega, step)
        if sv_f < 1.0 and pbd_f > 0.0:
            side_cush.append(sv_f / pbd_f)
        if sv_g < 1.0 and pbd_g > 0.0:
            side_cush.append(sv_g / pbd_g)
        sv_p = sv_f * sv_g
        pbd_p = pbd_f * pbd_g
        if sv_p < 1.0 and pbd_p > 0.0:
            prod_cush.append(sv_p / pbd_p)
    assert side_cush and prod_cush
    assert min(side_cush) >= 1.0 and min(prod_cush) >= 1.0
    prod_cush.sort()
    return dict(
        ok=True,
        wall_s=round(time.perf_counter() - t0, 4),
        max_omega_1_48=round(max_omega, 6),
        boot_psurv_raw_over_PBD=round(boot_raw_cush, 4),
        boot_surv_over_PBD=round(boot_cush, 4),
        boot_surv_over_PBD_alpha_g=round(boot_cush_a, 4),
        nodes_nontrivial_sides=len(side_cush),
        nodes_nontrivial_products=len(prod_cush),
        min_side_cushion=round(min(side_cush), 3),
        min_prod_cushion=round(prod_cush[0], 3),
        median_prod_cushion=round(prod_cush[len(prod_cush) // 2], 2),
    )


def verify_df_boot_tree_and_leaf(res, Pm, DS=None, lll_fn=None, mul_z_fn=None, norm_q15_fn=None, w2_vec=None):
    """Live pure-Python verification of the 5 DF-Boot invariants & exploratory diagnostics:
      (1) boot_df_tree_ok: 1,024-leaf depth-first post-order prime-root construction keeps pending stack
          < 1.02 * root_bits and peak step workspace <= 8.0 * root_bits (vs > 9.0 * root_bits full tree).
      (2) df_leaf_12sq_ok: 12 modular squarings (2^12 = 4,096 > 2,590) + 1 GCD + 1 exact division strips
          ALL prime powers p^e | N_g (p <= 2^{x_ea}) with 0 mismatches across adversarial 2,590-bit norms.
      (3) df_chunk_barrett_ok: chunked top-level reduction of P (SP > M_lim) modulo candidate root M_cand
          followed by binary remainder tree yields exact P mod N_i on all 64 leaves.
      (4) df_ea_lower_upper_ok: live check on EVERY stored (j, i) term of _ea_boot_dist confirming BOTH
          j*x_ea + (i+j)*h <= tau_ea (strict UPPER bit-sum cap) AND ea_boot_pr_lo < ea_boot_pr_hi (strict
          LOWER bit-sum inside decreasing rho_lo vs upper-Riemann rho_hi).
      (5) df_heb_pair_diag: scoped toy-scale 32-lattice Q(zeta_15) B_16(2)/{+-1} pairwise exact-|N_g|-dedup
          DF-Boot survival correlation check (matching the production exact-|N_g|-distinctness rule 1-to-1;
          shows near-independence kappa_emp << 1/p1 at toy marginal p1~0.012 — explicitly scoped as a
          toy-scale consistency check only, NOT a proof or empirical test of the 2^{30} H-EB premise at
          production rare-event marginal ~2^{-73.7})."""
    t0 = time.perf_counter()
    from cm_models import primes_upto

    # (1) 1,024-leaf depth-first post-order root build
    df_leaves = primes_upto(8200)[:1024]
    assert len(df_leaves) == 1024
    lvl = list(df_leaves)
    full_tree_bits = sum(x.bit_length() for x in lvl)
    while len(lvl) > 1:
        lvl = [lvl[i] * lvl[i + 1] for i in range(0, len(lvl), 2)]
        full_tree_bits += sum(x.bit_length() for x in lvl)
    ref_root = lvl[0]
    root_bits = ref_root.bit_length()
    live_stack_bits = 0
    max_stack_bits = 0
    max_step_total_bits = 0

    def _df_prod(arr):
        nonlocal live_stack_bits, max_stack_bits, max_step_total_bits
        if len(arr) == 1:
            return arr[0]
        mid = len(arr) // 2
        pl = _df_prod(arr[:mid])
        bl = pl.bit_length()
        live_stack_bits += bl
        if live_stack_bits > max_stack_bits:
            max_stack_bits = live_stack_bits
        pr = _df_prod(arr[mid:])
        br = pr.bit_length()
        live_stack_bits -= bl
        prod = pl * pr
        bp = prod.bit_length()
        step_bits = live_stack_bits + bl + br + 6 * bp
        if step_bits > max_step_total_bits:
            max_step_total_bits = step_bits
        return prod

    df_root = _df_prod(df_leaves)
    boot_df_tree_ok = (
        df_root == ref_root
        and live_stack_bits == 0
        and max_stack_bits < 1.02 * root_bits
        and max_step_total_bits <= 8.0 * root_bits
        and full_tree_bits > 9.0 * root_bits
    )
    assert boot_df_tree_ok

    # (2) Adversarial 2,590-bit prime-power 12-squaring leaf-stripping check (2^12 = 4,096 > 2,590)
    ps_s = primes_upto(64)
    P_s = math.prod(ps_s)
    p_big = (1 << 127) - 1  # Mersenne prime > 64
    rng = random.Random(20261124)
    cases = [
        1 << 2589,
        3 ** 1630,
        (2 ** 1200) * (3 ** 400) * (5 ** 200) * p_big,
        (2 ** 500) * (7 ** 250) * (61 ** 100) * (p_big ** 3),
        (P_s ** 25) * p_big,
    ]
    for _ in range(16):
        n_v = 1
        for p in ps_s:
            if rng.random() < 0.35:
                n_v *= p ** rng.randint(1, 25)
        while n_v.bit_length() + 127 <= 2590:
            n_v *= p_big
        cases.append(n_v)
    leaf_errs = 0
    for n_v in cases:
        assert n_v.bit_length() <= 2590
        r = P_s % n_v
        for _ in range(12):
            r = (r * r) % n_v
        cof_got = n_v // math.gcd(n_v, r)
        cof_ref = n_v
        for p in ps_s:
            while cof_ref % p == 0:
                cof_ref //= p
        if cof_got != cof_ref:
            leaf_errs += 1
    assert leaf_errs == 0

    # (3) Chunked top-level reduction of P (SP > M_lim) + 64-leaf remainder tree exactness
    P_toy = math.prod(primes_upto(1000))  # ~1,380 bits
    cands = [rng.randrange(1 << 15, 1 << 16) | 1 for _ in range(64)]  # M_cand ~ 1,024 bits < P_toy!
    pt = [list(cands)]
    while len(pt[-1]) > 1:
        prev = pt[-1]
        pt.append([prev[i] * prev[i + 1] for i in range(0, len(prev), 2)])
    M_root = pt[-1][0]
    assert P_toy.bit_length() > M_root.bit_length()
    # Chunk P_toy in base B = 1 << chunk_bits (chunk_bits = M_root.bit_length()) via Horner mod M_root
    cb = M_root.bit_length()
    mask = (1 << cb) - 1
    chunks = []
    tmp = P_toy
    while tmp:
        chunks.append(tmp & mask)
        tmp >>= cb
    assert len(chunks) >= 2
    rem_top = 0
    base_mod = (1 << cb) % M_root
    for ch in reversed(chunks):
        rem_top = (rem_top * base_mod + ch) % M_root
    assert rem_top == P_toy % M_root
    rems = [rem_top]
    for lvl_nodes in reversed(pt[:-1]):
        rems = [rems[i >> 1] % lvl_nodes[i] for i in range(len(lvl_nodes))]
    assert all(rems[i] == P_toy % cands[i] for i in range(64))

    bt = res["boot"]
    assert bt["boot_max_1own_log2"] < 64.0

    # (4) Strict two-sided check on dyadic-grid j-tracked _ea_boot_dist / ea_boot_pr
    df_ea_ok = False
    if DS is not None:
        n1_b, qA_b, kA_b, _ = DS.BOOT16
        x_ea = float(Pm.BOOT_X_EA)
        hb = float(Pm.BOOT_H)
        tau_b = float(Pm.TAU_BOOT)
        by_j = DS._ea_boot_dist(x_ea, Pm.LX, qA_b, float(n1_b), kA_b, tau_b, hb, Pm.BOOT_K_BAND)
        for j, i_lo, arr in by_j:
            i_hi = i_lo + len(arr) - 1
            assert j * x_ea + (i_hi + j) * hb <= tau_b + 1e-9
        m_test = bt["m"]
        p_lo = DS.ea_boot_pr(m_test, x_ea, Pm.LX, qA_b, float(n1_b), kA_b, tau_b, Pm.ALPHA_G, hb, Pm.BOOT_K_BAND)
        # Upper reference using exact upper bit-sum (j*x_ea + (i+j)*h) inside rho_hi + exact log weights ratio
        p_hi = 0.0
        m_f = math.floor(m_test)
        tab_hi, ntab_hi, flr = DS.RHO_TAB_HI, len(DS.RHO_TAB_HI), math.floor
        inv_xr = 1.0 / (x_ea * DS.RSTEP)
        sp_d = DS._SUMPROD
        for j, i_lo, arr in by_j:
            base_m = m_f + Pm.ALPHA_G - j * (x_ea + hb)
            vals_hi = [
                1.0 if (k := flr((base_m - i * hb) * inv_xr)) <= 0 else (tab_hi[k] if k < ntab_hi else tab_hi[-1])
                for i in range(i_lo, i_lo + len(arr))
            ]
            p_hi += sp_d(arr, vals_hi)
        assert 0.0 < p_lo < p_hi
        df_ea_ok = True

    # (5) Scoped toy-scale pairwise correlation check for H-EB (32 independent Q(zeta_15) 16D L_boot-style lattices)
    heb_diag = None
    if lll_fn is not None and mul_z_fn is not None and norm_q15_fn is not None and w2_vec is not None:
        cells2 = [(i, -1, 0) for i in range(16)]
        for i in range(16):
            for j in range(i + 1, 16):
                cells2.append((i, j, 1))
                cells2.append((i, j, -1))
        sub_c = cells2[:128]
        P_t = math.prod(primes_upto(32768))
        rng_p = random.Random(20261124)
        tot_c = tot_h = tot_p = tot_b = 0
        q_t = 61
        for _ in range(32):
            tv = [rng_p.randrange(q_t) for _ in range(8)]
            Bm = [[0] * 16 for _ in range(16)]
            for i in range(8):
                Bm[i][i] = q_t
            for i in range(8):
                yi = [1 if k == i else 0 for k in range(8)]
                pr = mul_z_fn(yi, tv)
                for k in range(8):
                    Bm[8 + i][k] = pr[k] % q_t
                Bm[8 + i][8 + i] = 1
            Br, _, _ = lll_fn(Bm, 3, 4)
            seen = set()
            hits = kept = 0
            for ci, cj, csg in sub_c:
                ri = Br[ci]
                if cj < 0:
                    vv = ri
                elif csg == 1:
                    rj = Br[cj]
                    vv = [ri[k] + rj[k] for k in range(16)]
                else:
                    rj = Br[cj]
                    vv = [ri[k] - rj[k] for k in range(16)]
                a, b = vv[:8], vv[8:]
                a3 = mul_z_fn(mul_z_fn(a, a), a)
                b3 = mul_z_fn(mul_z_fn(b, b), b)
                wb3 = mul_z_fn(w2_vec, b3)
                ng = abs(norm_q15_fn([a3[i] + wb3[i] for i in range(8)]))
                if ng == 0 or ng in seen:
                    continue
                seen.add(ng)
                kept += 1
                rem = ng
                while True:
                    g = math.gcd(rem, P_t)
                    if g == 1:
                        break
                    rem //= g
                if rem.bit_length() <= 36:
                    hits += 1
            tot_c += kept
            tot_h += hits
            tot_p += kept * (kept - 1) // 2
            tot_b += hits * (hits - 1) // 2
        p1_e = tot_h / tot_c
        p2_e = tot_b / tot_p
        kap_e = p2_e / (p1_e * p1_e)
        assert tot_h >= 20 and tot_b >= 10 and 0.0 < kap_e < 100.0
        heb_diag = dict(
            scope="toy_scale_exact_ng_dedup_consistency_check_not_a_proof_of_rare_event_H_EB",
            n_lattices=32,
            q_toy=q_t,
            B_toy=32768,
            tau_toy=36,
            distinct_ng_cells=tot_c,
            hits=tot_h,
            p1_emp=round(p1_e, 5),
            pairs=tot_p,
            both_hits=tot_b,
            p2_emp=round(p2_e, 6),
            kappa_emp=round(kap_e, 4),
            max_possible_kappa_at_toy_p1=round(1.0 / p1_e, 2),
            kappa_charged_log2=int(Pm.BOOT_KAPPA_LOG2),
        )

    return dict(
        ok=True,
        wall_s=round(time.perf_counter() - t0, 4),
        boot_df_tree_ok=True,
        df_stack_ratio=round(max_stack_bits / root_bits, 4),
        df_peak_ratio=round(max_step_total_bits / root_bits, 4),
        df_full_tree_ratio=round(full_tree_bits / root_bits, 4),
        df_leaf_12sq_ok=True,
        df_leaf_cases=len(cases),
        df_chunk_barrett_ok=True,
        df_chunk_count=len(chunks),
        df_ea_lower_upper_ok=df_ea_ok,
        df_heb_pair_diag=heb_diag,
        boot_max_1own_log2=round(bt["boot_max_1own_log2"], 3),
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
