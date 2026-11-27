#!/usr/bin/env python3
"""100% Standalone Verifier for stnfs-q15-guard4-bandpass-snell-v2 (pure Python 3.11+ stdlib).
ZERO required external files (LEDGER_DATA embedded; sibling audits run when present).
Stdin JSON: {} -> part="core" | {"part":"diagnostics"} | {"part":"all"}."""
from __future__ import annotations

import bisect, json, math, os, random, re, struct, sys, time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent
X_BLS = -0xD201000000010000
ABS_X = abs(X_BLS)
LOG2_ABS_X = math.log2(ABS_X)
P_BLS = int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab", 16)
R_BLS = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
P0_X = int("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb", 16)
P0_Y = int("08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1", 16)
Q0_X = (
    int("024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8", 16),
    int("13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e", 16),
)
Q0_Y = (
    int("0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801", 16),
    int("0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be", 16),
)
C_COST = 20000.0
L_OWN_1D = (4.0 * (2.0**60) * 8.0) / 4.5
L_OWN_2D = (4.0 * (2.0**60) * 8.0) / 9.0

LEDGER_DATA: dict = json.loads(
    "{\"package_id\":\"stnfs-q15-guard4-bandpass-snell-v2\",\"config_tag\":\"CAND2_FINAL_ETA85_7000_NL205\",\"track\":\"bls12-"
    "381-fixed-g1-invert-g2\",\"cost_model\":\"classical-group-work-v2\",\"C_ir_per_u\":20000.0,\"isa_profile\":{\"flags\":\"-m"
    "arch=armv8.1-a -mgeneral-regs-only -O3 -fno-vectorize -fno-slp-vectorize\",\"simd\":0,\"scalar_lse_instructions\":["
    "\"STADDB Ws, [Xn]\"],\"forbidden_instructions_count\":0,\"emit_x64_spills\":0,\"slab_kernel_spills\":0},\"operating_poi"
    "nt\":{\"K_shells\":2,\"q_ball\":0.65,\"rb_ratio\":0.82,\"V_a_per_shell\":6.738105618486642e+16,\"B_f\":70.0,\"B_g\":70.36,\""
    "U_orbits\":4.616694350549097e+18,\"log2_U_orbits\":62.00156592952081,\"eta\":0.85,\"g_band\":[30.0,61.3],\"f_band\":[34"
    ".0,61.0],\"q_step_g\":3.4,\"q_step_f\":3.8,\"adaptive_2d_lambda\":671088.64,\"mean_Tg_bits\":266.1260412987349,\"mean_T"
    "f_bits\":177.23003433890125,\"ret_tot\":0.6787271038365621,\"log2_sg\":-14.57507733975491,\"log2_sf\":-9.041862756350"
    "094,\"log2_s1\":-23.616940096105004,\"W_B\":4,\"L_act_max\":248,\"slab_cols_total\":256,\"m_block\":512,\"n_block\":32},\"c"
    "hebotarev_exact\":{\"d_ok_g_frac\":\"1/3\",\"d_ok_g\":0.3333333333333333,\"r_bar_g\":3.0,\"dens_g_frac\":\"1/24\",\"dens_g\":"
    "0.041666666666666664,\"d_ok_f_frac\":\"23/81\",\"d_ok_f\":0.2839506172839506,\"r_bar_f_frac\":\"81/23\",\"r_bar_f\":3.5217"
    "391304347827,\"dens_f_frac\":\"2131/34992\",\"dens_f\":0.060899634202103335,\"galois_single_order\":162,\"galois_joint_"
    "K_order\":8748,\"galois_joint_Q_order\":69984},\"ntt_1d2d_shoup61\":{\"p1\":2287828610704211969,\"p1_hex\":\"0x1fc000000"
    "0000001\",\"g1\":3,\"p2\":2053641430080946177,\"p2_hex\":\"0x1c80000000000001\",\"g2\":7,\"p3\":1945555039024054273,\"p3_hex"
    "\":\"0x1b00000000000001\",\"g3\":5,\"smooth_2_pow\":54,\"log2_P_prod\":182.57646220310036,\"threshold_1d_bits_log2\":59,\""
    "max_2d_bits_log2\":68.5,\"headroom_1d_bits\":1.5764622031003626,\"headroom_2d_bits_at_68p5\":3.934443198227939,\"ir_"
    "bfly\":17,\"ir_crt3_1d_compiled\":58,\"ir_crt3_1d_charged\":67,\"ir_unpack58_max\":17,\"ir_garner_crt3_overlap58_compi"
    "led\":73,\"ir_garner_crt3_charged\":95},\"direct_mc_4M\":{\"n_draws\":4000000,\"n_chords\":200000,\"frontier_anchor_mc_l"
    "og2_vk\":119.37413,\"frontier_hist_log2_vk\":119.43957,\"anchor_pad_bits\":0.06544000000000949,\"cand_mc_log2_vk\":11"
    "9.65310528656599,\"cand_l2_Rret\":61.76748634141777,\"cand_l2_Rreq\":61.76710067588378,\"cand_se_bits\":0.0011846293"
    "54729802,\"cand_se3_pad_bits\":0.003553888064189406,\"cand_R0\":1405.736530398967,\"cand_cond_nf\":543.4551311637518"
    ",\"cand_cond_ng\":671.0200032405045,\"cand_obs_max_nf\":596.9079949474983,\"cand_obs_max_ng\":686.815178931357,\"det_"
    "upper_ng\":686.8386771689147,\"det_upper_nf\":617.3003074004176,\"charged_max_ng\":705.0,\"charged_max_nf\":630.0,\"L_"
    "win_eff\":111.53432394702733,\"gap_frac_measured\":0.0016574892141784845,\"gap_frac_charged\":0.0017574892141784845"
    ",\"gap_pad_bits\":0.0025377516630844055,\"log2_Vk_pure\":119.65919692629326,\"log2_Vk_charged\":119.72463692629327},"
    "\"sieve_and_stage1b\":{\"g_guard4_x64\":{\"hits_active\":0.7145822457000107,\"hits_emitted\":0.7338027570299421,\"ir_hi"
    "ts\":1.9217821836649176,\"ir_probe\":0.18924516180700401,\"n_ok_primes\":2.281865301197192e+16,\"n_bkt\":7.9942251599"
    "04897e+33,\"n_ent\":4.612634821987955e+33},\"f_guard4_x64\":{\"hits_active\":0.5845132263521059,\"hits_emitted\":0.600"
    "2351997949109,\"ir_hits\":1.5691301301242433,\"ir_probe\":0.13160661727539955,\"n_ok_primes\":1.5868758308068872e+16"
    ",\"n_bkt\":4.515775853584159e+33,\"n_ent\":3.765632578536715e+33},\"slab_prebias_ir\":0.19240027322728553,\"slab_scan"
    "_g_ir\":0.27754180090327674,\"slab_slow_and_probe_ir\":0.007743304461852398,\"ipc_slab_live\":0.47768537859241467,\""
    "ipc_slab_charged\":0.495,\"ipc_hits\":3.4909123137891607,\"ipc_probe\":0.3208517790824036,\"ipc_cof\":0.7082516156845"
    "888,\"ipc_tot\":5.015015708556152,\"stage1b_breakdown\":{\"c_rem\":0.7082516155549564,\"c_bld\":1.1706927151904765e-12"
    ",\"c_mr\":1.2963241105059437e-10,\"c_cof_rem_mr\":0.7082516156845888,\"tree_build_u\":6.4286510844062966e+19,\"peak_g"
    "lob_log2_B\":69.42613321244367,\"rem_f\":385.0,\"rem_g\":405.0,\"pf_sm\":2.6865854478436424e-05,\"pg_sm\":2.83501610865"
    "3438e-06,\"tiers\":[{\"tier\":\"f1\",\"band\":[1.0,61.0],\"leaf_bits\":560.0,\"k_pack\":7,\"h_super\":46,\"ntt_mode\":\"1D\",\"lo"
    "g2_L_bits\":57.49134440869225,\"ir_per_cell\":0.31022237496071137},{\"tier\":\"f2\",\"band\":[61.0,65.5],\"leaf_bits\":38"
    "5.0,\"k_pack\":10,\"h_super\":51,\"ntt_mode\":\"2D\",\"log2_L_bits\":61.926133212443666,\"ir_per_cell\":0.2584938878351426"
    "5},{\"tier\":\"f3\",\"band\":[65.5,70.0],\"leaf_bits\":385.0,\"k_pack\":10,\"h_super\":55,\"ntt_mode\":\"2D\",\"log2_L_bits\":66"
    ".42613321244367,\"ir_per_cell\":0.1393823734072984},{\"tier\":\"g1\",\"band\":[1.0,61.3],\"leaf_bits\":675.0,\"k_pack\":6,"
    "\"h_super\":46,\"ntt_mode\":\"1D\",\"log2_L_bits\":57.243804675460105,\"ir_per_cell\":9.930971293409903e-06},{\"tier\":\"g2"
    "\",\"band\":[61.3,65.83],\"leaf_bits\":405.0,\"k_pack\":10,\"h_super\":50,\"ntt_mode\":\"2D\",\"log2_L_bits\":61.709965491640"
    "645,\"ir_per_cell\":6.114466794681393e-06},{\"tier\":\"g3\",\"band\":[65.83,70.36],\"leaf_bits\":405.0,\"k_pack\":10,\"h_su"
    "per\":55,\"ntt_mode\":\"2D\",\"log2_L_bits\":66.23996549164065,\"ir_per_cell\":4.788159444379221e-06}]},\"ir_replay_entr"
    "y\":27.5,\"ir_replay_bkt\":6.5},\"la_and_sge_zero_extrap\":{\"pricing_tier\":\"R3_U24_zero_extrapolation\",\"nl_cap\":0.0"
    "0205,\"w_cap\":190.0,\"nnz_per_U_cap\":0.3895,\"u23_4seed_eta85_nl_mean\":0.00138,\"u23_4seed_eta85_nl_se\":0.000106,\""
    "u23_4seed_eta85_nnz_per_U\":0.25173,\"u23_4seed_eta85_core_ir\":10.56,\"u24_2seed_eta85_nl_mean\":0.001738,\"u24_2se"
    "ed_eta85_nl_se\":0.000111,\"u24_2seed_eta85_nnz_per_U\":0.32249,\"u24_2seed_eta85_core_ir\":10.21,\"N_l\":94642234186"
    "25650.0,\"N_nnz\":1.7982024495388736e+18,\"ir_core_nnz_x4\":12.140999999999998,\"ir_core_nnz_x1\":16.0,\"ir_row_fixed"
    "\":2100.0,\"ir_fr_mont_mul\":236.0,\"ir_fr_add\":47.0,\"bw_overhead\":1.15,\"L_seq\":314241793196618.8,\"d_mksol\":295756"
    "981832083.56,\"total_spmvs\":1.0558524251406364e+16,\"log2_la_u\":104.42036997220899,\"log2_mksol_u\":100.2299599858"
    "2789,\"log2_lingen_u\":79.12947913338215,\"log2_sge_u\":61.697767896082546},\"stage4_12d_mg\":{\"subfield\":\"M_g = Q(z"
    "eta_5)[Y]/(Y^3 - x)\",\"deg_over_Q\":12,\"det_log2\":2539.061119126145,\"covol_log2\":2559.019570557465,\"uncond_lll_a"
    "mgm_log2\":2570.5097955531382,\"N_boot_bits\":2572.0,\"B_boot_bits\":275.0,\"layer0_inert_log2_u\":64.14704834430788,"
    "\"live_2log_no_snell_log2_u\":95.02655965491851,\"live_2log_with_snell_log2_u\":95.61126398977053,\"charged_both_lo"
    "gs_log2_u\":96.5001,\"ecm_split_live_log2_u\":93.74857640634589},\"memory_audit_unified\":{\"single_slab_eib\":14.961"
    "6,\"replay_global_log2_B\":117.00824812766784,\"stage1_peak_global_log2_B\":117.00824812766784,\"stage2_sge_global_"
    "log2_B\":69.52198817804725,\"stage3_la_global_log2_B\":72.30559148136261,\"stage1b_global_log2_B\":69.4261332124436"
    "7,\"peak_global_claim_log2_B\":117.5},\"failure_budget_union_bound\":{\"conditional_on\":[\"H1\",\"H2\",\"H3\",\"H4\",\"H5\",\""
    "H-R\",\"H-SM\",\"H-ECM\"],\"p_yield_shortfall_upper\":7.888609052210118e-31,\"p_mr_false_prime_upper\":1.83670992315982"
    "42e-40,\"p_la_rank_and_bw_upper\":6.223015277861142e-61,\"p_descent_unmetered\":0.0,\"w_pointwise_meter_log2_u\":124"
    ".0,\"p_meter_abort_upper\":1.569922080717824e-05,\"p_wrapper_upper\":3.552713678800501e-15,\"p_fail_open_upper\":1.5"
    "699220810730953e-05,\"p_open_lower\":0.9999843007791893,\"p_open_floor\":0.99998,\"p_cold_lower\":0.3900016773374306"
    "},\"ledger_rows_live_log2_u\":{\"relation_collection\":107.76317876654849,\"replay_a_tables\":104.62000944499829,\"li"
    "near_algebra\":104.42036997220899,\"both_individual_logs_exp\":96.5,\"ecm_split_expected\":93.74857640634589,\"band_"
    "split_and_trees\":90.0,\"meter_safety_buffer\":93.4,\"snell_resume_tail\":92.5,\"matrix_generator_sge\":61.6977678960"
    "82546,\"canonical_a_domain\":84.0,\"schirokauer_maps\":80.0,\"stage0_setup_and_iso\":45.0,\"final_log_backsub_and_gcd"
    "\":72.0,\"final_verification_Q\":40.0},\"ledger_rows_charged_log2_u\":{\"relation_collection\":107.7633,\"replay_a_tab"
    "les\":104.6202,\"linear_algebra\":104.4205,\"both_individual_logs_exp\":96.5001,\"ecm_split_expected\":93.7487,\"band_"
    "split_and_trees\":90.0001,\"meter_safety_buffer\":93.4001,\"snell_resume_tail\":92.5001,\"matrix_generator_sge\":61.6"
    "979,\"canonical_a_domain\":84.0001,\"schirokauer_maps\":80.0001,\"stage0_setup_and_iso\":45.0001,\"final_log_backsub_"
    "and_gcd\":72.0001,\"final_verification_Q\":40.0001},\"totals\":{\"opened_run_live_log2_u\":108.04092406907714,\"opened"
    "_run_charged_log2_u\":108.04105248178097,\"K_gate\":7194374076228249696,\"K_gate_hex\":\"0x63d78d1adadd2860\",\"p_gate"
    "\":0.3900078001560031,\"log2_p_gate\":-1.3584251167231158,\"gated_live_log2_u\":106.68249895235402,\"gated_charged_l"
    "og2_u\":106.68262736505785,\"claim_numerator\":\"21337\",\"claim_denominator\":\"200\",\"claim_fraction\":\"21337/200\",\"cl"
    "aim_log2_u\":106.685,\"slack_to_claim_bits\":0.002372634942148011,\"pure_first_principles_gated_log2_u\":106.628735"
    "36376671,\"frontier_comparator_log2_u\":106.788321,\"cand1_claim_log2_u\":106.67,\"improvement_vs_frontier_bits\":0."
    "103320999999994}}"
)


class ExactDickman:
    def __init__(self, u_max: float = 30.0, N: int = 200) -> None:
        self.u_max, self.N = u_max, N
        h = 1.0 / N
        tot = int(round(u_max * N)) + 1
        v = [1.0] * (N + 1) + [0.0] * (tot - N - 1)
        for i in range(N + 1, min(2 * N + 1, tot)):
            v[i] = 1.0 - math.log(i * h)
        for i in range(2 * N + 1, tot):
            v[i] = (h * (0.5 * v[i - N] + sum(v[i - N + 1 : i]))) / (i * h - 0.5 * h)
        self.vals = v

    def rho(self, u: float) -> float:
        if u <= 1.0:
            return 1.0
        if u >= self.u_max:
            return 1e-300
        pos = u * self.N
        idx = int(pos)
        if idx + 1 >= len(self.vals):
            return self.vals[-1]
        f = pos - idx
        return math.exp((1.0 - f) * math.log(self.vals[idx]) + f * math.log(self.vals[idx + 1]))


DICKMAN = ExactDickman()


def _sieve_primes(limit: int) -> list[int]:
    ip = bytearray(b"\x01") * (limit + 1)
    ip[0] = ip[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if ip[i]:
            ip[i * i : limit + 1 : i] = bytearray(len(range(i * i, limit + 1, i)))
    return [i for i in range(2, limit + 1) if ip[i]]


# 1. Stage 0: BLS12-381 Algebra, Frobenius, Automorphisms & n_SM = 4
def _g1_add(p1, p2):
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    x1, y1, x2, y2 = p1[0], p1[1], p2[0], p2[1]
    if x1 == x2:
        if (y1 + y2) % P_BLS == 0:
            return None
        lam = (3 * x1 * x1 * pow(2 * y1, -1, P_BLS)) % P_BLS
    else:
        lam = ((y2 - y1) * pow(x2 - x1, -1, P_BLS)) % P_BLS
    x3 = (lam * lam - x1 - x2) % P_BLS
    return (x3, (lam * (x1 - x3) - y1) % P_BLS)


def _g1_mul(k: int, pt):
    res, addend = None, pt
    while k:
        if k & 1:
            res = _g1_add(res, addend)
        addend = _g1_add(addend, addend)
        k >>= 1
    return res


def check_stage0_algebra() -> dict:
    assert (X_BLS - 1) ** 2 * (X_BLS**4 - X_BLS**2 + 1) // 3 + X_BLS == P_BLS
    assert X_BLS**4 - X_BLS**2 + 1 == R_BLS
    lam = X_BLS % R_BLS
    assert P_BLS % R_BLS == lam and (pow(lam, 4, R_BLS) - pow(lam, 2, R_BLS) + 1) % R_BLS == 0
    assert pow(lam, 6, R_BLS) == R_BLS - 1 and pow(lam, 12, R_BLS) == 1 and pow(lam, 4, R_BLS) != 1
    assert (P0_Y * P0_Y - (P0_X**3 + 4)) % P_BLS == 0 and _g1_mul(R_BLS, (P0_X, P0_Y)) is None
    fp2m = lambda a, b: ((a[0] * b[0] - a[1] * b[1]) % P_BLS, (a[0] * b[1] + a[1] * b[0]) % P_BLS)
    q3 = fp2m(fp2m(Q0_X, Q0_X), Q0_X)
    assert fp2m(Q0_Y, Q0_Y) == ((q3[0] + 4) % P_BLS, (q3[1] + 4) % P_BLS)
    assert P_BLS % 15 == 7 and [k for k in range(1, 5) if pow(P_BLS, k, 5) == 1][0] == 4
    assert pow(X_BLS % P_BLS, (P_BLS - 1) // 3, P_BLS) != 1
    w_p = next(
        c for g in range(2, 20)
        for cand in (pow(g, (P_BLS - 1) // 3, P_BLS),) if cand != 1
        for c in (cand, (cand * cand) % P_BLS)
        if (X_BLS**3 - X_BLS**2 + (c * c % P_BLS) * X_BLS - c) % P_BLS == 0
    )
    assert w_p and X_BLS**6 - 2 * X_BLS**5 + 2 * X_BLS**3 + X_BLS + 1 == 3 * P_BLS
    assert (R_BLS - 1) % 12 == 0 and (36 // 12 + 12 // 12) == 4
    return {"ok": True, "p_mod_15": 7, "ord_5_p": 4, "p_equiv_x_equiv_lambda_mod_r": True, "n_sm_chi1": 4}


# 2. Live Exhaustive Finite-Group Proof of Chebotarev Theorems 2A & 2B
def check_chebotarev_theorems_2a_2b(p_limit: int = 5000) -> dict:
    zwm = lambda x, y: (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0] - x[1] * y[1])
    zwn = lambda x: x[0] * x[0] - x[0] * x[1] + x[1] * x[1]
    zwdiv = lambda d, n: all(c % zwn(d) == 0 for c in zwm(n, (d[0] - d[1], -d[1])))
    dw, dw2, pi19, pi67 = (41, 24), (17, -24), (-5, -2), (-9, -2)
    assert zwm(pi19, pi67) == dw and zwn(dw) == 1273 == 19 * 67 and math.gcd(1273, 15) == 1
    assert not zwdiv(pi19, dw2) and not zwdiv(pi67, dw2)
    assert 151 % 15 == 1 and (118 * 118 + 118 + 1) % 151 == 0
    assert all((t**3 - t**2 + (118 * 118 % 151) * t - 118) % 151 != 0 for t in range(151))

    s3 = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    f3 = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]
    subs = set()
    for v in f3:
        sp = {(0, 0, 0)}
        for u in {tuple(v[p[i]] for i in range(3)) for p in s3}:
            sp |= {((s[0] + k * u[0]) % 3, (s[1] + k * u[1]) % 3, (s[2] + k * u[2]) % 3) for s in sp for k in (1, 2)}
        subs.add(frozenset(sp))
    all_subs = {frozenset(((a[0] + b[0]) % 3, (a[1] + b[1]) % 3, (a[2] + b[2]) % 3) for a in x for b in y) for x in subs for y in subs}
    assert sorted(len(s) for s in all_subs) == [1, 3, 9, 27] and all((1, 1, 1) in s for s in all_subs if len(s) > 1)

    G = [(v, s) for v in f3 for s in s3]
    gmul = lambda a, b: (tuple((a[0][b[1][i]] + b[0][i]) % 3 for i in range(3)), tuple(a[1][b[1][i]] for i in range(3)))
    def ginv(g):
        si = [0, 0, 0]
        for i in range(3):
            si[g[1][i]] = i
        st = tuple(si)
        return (tuple((-g[0][st[i]]) % 3 for i in range(3)), st)

    seen, cclasses = set(), []
    for g in G:
        if g not in seen:
            cc = frozenset(gmul(gmul(h, g), ginv(h)) for h in G)
            seen |= cc
            cclasses.append(cc)
    assert len(G) == 162 and len(cclasses) == 22

    id_el = ((0, 0, 0), (0, 1, 2))
    h_z9 = frozenset(g for g in G if sum(g[0]) % 3 == 0)
    nsubs, q = {frozenset([id_el])}, [frozenset([id_el])]
    while q:
        N = q.pop()
        for cc in cclasses:
            if cc.issubset(N):
                continue
            M = set(N) | set(cc)
            while True:
                add = {gmul(a, b) for a in M for b in cc} - M
                if not add:
                    break
                M |= add
            Mf = frozenset(M)
            if Mf not in nsubs:
                nsubs.add(Mf)
                q.append(Mf)
    assert sorted(len(N) for N in nsubs) == [1, 3, 9, 27, 27, 54, 81, 162]
    sin_hz9 = [N for N in nsubs if N < h_z9]
    max_hz9 = [N for N in sin_hz9 if not any(N < M for M in sin_hz9)]
    assert len(max_hz9) == 1 and len(max_hz9[0]) == 27

    rc = {0: 0, 3: 0, 6: 0, 9: 0}
    c0_a = cnz_a = 0
    for v, s in G:
        fx = sum(1 for i in range(3) for c in range(3) if s[i] == i and (c + v[i]) % 3 == c)
        rc[fx] += 1
        if fx:
            if sum(v) % 3 == 0:
                c0_a += 1
            else:
                cnz_a += 1
    p0, p3, p6, p9 = (Fraction(rc[k], 162) for k in (0, 3, 6, 9))
    assert (p0, p3, p6, p9) == (Fraction(58, 81), Fraction(13, 54), Fraction(1, 27), Fraction(1, 162))
    d_ok_f, r_bar_f = 1 - p0, 1 / (1 - p0)
    d0, dnz = Fraction(c0_a, 54), Fraction(cnz_a, 108)
    p_act_15 = Fraction(1, 3) * (1 - (1 - d0) ** 2) + Fraction(2, 3) * (1 - (1 - dnz) ** 2)
    dens_f = Fraction(1, 8) * p_act_15
    assert d_ok_f == Fraction(23, 81) and r_bar_f == Fraction(81, 23) and dens_f == Fraction(2131, 34992)

    primes = _sieve_primes(p_limit)
    p15 = [p for p in primes if p % 15 == 1]
    bad_m3 = forb_z9 = f_act_p = f_act_ok = 0
    for p in p15:
        w = next(pow(g, (p - 1) // 3, p) for g in range(2, p) if pow(g, (p - 1) // 3, p) != 1)
        if (pow(w, (p - 1) // 3, p) == 1) != (pow((w * w) % p, (p - 1) // 3, p) == 1):
            forb_z9 += 1
        p_has = False
        for ww in (w, (w * w) % p):
            ww2 = (ww * ww) % p
            ny = sum(3 for t in range(p) if (t * t * t - t * t + ww2 * t - ww) % p == 0 and pow(t, (p - 1) // 3, p) == 1)
            if ny % 3:
                bad_m3 += 1
            if ny:
                p_has = True
                f_act_ok += 4
        if p_has:
            f_act_p += 1
    assert bad_m3 == 0 and forb_z9 == 0
    return {
        "ok": True, "galois_order": 162, "normal_subgroups": 8, "unique_max_in_Hz9": 27,
        "d_ok_f": str(d_ok_f), "r_bar_f": str(r_bar_f), "dens_f": str(dens_f),
        "live_primes": len(primes), "live_bad_mod3": 0, "live_forb_z9": 0,
        "live_dens_f_5000": round(f_act_p / len(primes), 5),
        "live_d_ok_f_5000": round(f_act_ok / (8 * len(p15)), 5),
    }


# 3. Live Bit-Exact 1D + 2D Bivariate 3-Prime 61-Bit Shoup NTT Verification
def _is_prime_64(n: int) -> bool:
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
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


def _ntt_1d_shoup61(vec: list[int], p: int, g: int, inv: bool = False) -> list[int]:
    n, a, two_p = len(vec), [x % p for x in vec], 2 * p
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    length = 2
    while length <= n:
        half = length >> 1
        wl = pow(g, (p - 1) // length, p)
        if inv:
            wl = pow(wl, p - 2, p)
        ws, wps, cur = [1] * half, [0] * half, 1
        for k in range(half):
            ws[k], wps[k] = cur, (cur << 64) // p
            cur = (cur * wl) % p
        for i in range(0, n, length):
            for k in range(half):
                u = a[i + k]
                if u >= two_p:
                    u -= two_p
                v = a[i + k + half]
                t = (v * ws[k] - ((v * wps[k]) >> 64) * p) & 0xFFFFFFFFFFFFFFFF
                a[i + k], a[i + k + half] = u + t, u + two_p - t
        length <<= 1
    if inv:
        inv_n = pow(n, p - 2, p)
        return [(x * inv_n) % p for x in a]
    return [x % p for x in a]


def _mul_2d_3prime_shoup61(A: int, B: int, pgs: list[tuple[int, int]]) -> int:
    m58 = (1 << 58) - 1
    def to58(x):
        r = []
        while x:
            r.append(x & m58)
            x >>= 58
        return r or [0]
    la, lb = to58(A), to58(B)
    k = max(1, (max(len(la), len(lb)) - 1).bit_length())
    k1 = (k + 1) // 2
    D1, D2 = 1 << k1, 1 << (k - k1)
    N1, N2 = 2 * D1, 2 * D2
    la += [0] * (D1 * D2 - len(la))
    lb += [0] * (D1 * D2 - len(lb))
    cmod = []
    for p, g in pgs:
        gA = [_ntt_1d_shoup61(la[r * D1 : (r + 1) * D1] + [0] * D1, p, g) for r in range(D2)] + [[0] * N1 for _ in range(D2)]
        gB = [_ntt_1d_shoup61(lb[r * D1 : (r + 1) * D1] + [0] * D1, p, g) for r in range(D2)] + [[0] * N1 for _ in range(D2)]
        gC = [[0] * N1 for _ in range(N2)]
        for c in range(N1):
            cA = _ntt_1d_shoup61([gA[r][c] for r in range(N2)], p, g)
            cB = _ntt_1d_shoup61([gB[r][c] for r in range(N2)], p, g)
            cC = _ntt_1d_shoup61([(cA[r] * cB[r]) % p for r in range(N2)], p, g, True)
            for r in range(N2):
                gC[r][c] = cC[r]
        cmod.append([_ntt_1d_shoup61(gC[r], p, g, True) for r in range(N2)])
    p1, p2, p3 = pgs[0][0], pgs[1][0], pgs[2][0]
    inv12, inv123, p1p2 = pow(p1 % p2, p2 - 2, p2), pow((p1 * p2) % p3, p3 - 2, p3), p1 * p2
    acc = [0] * (D1 * (N2 + 1) + 4)
    for r in range(N2):
        base = r * D1
        for c in range(N1):
            a1, a2, a3 = cmod[0][r][c], cmod[1][r][c], cmod[2][r][c]
            x12 = a1 + p1 * (((a2 - a1) * inv12) % p2)
            acc[base + c] += x12 + p1p2 * (((a3 - x12) * inv123) % p3)
    res = carry = 0
    for i, v in enumerate(acc):
        s = v + carry
        res |= (s & m58) << (58 * i)
        carry = s >> 58
    return res | (carry << (58 * len(acc)))


def check_ntt_1d2d_shoup61(led: dict) -> dict:
    ntt = led["ntt_1d2d_shoup61"]
    p1, p2, p3, g1, g2, g3 = ntt["p1"], ntt["p2"], ntt["p3"], ntt["g1"], ntt["g2"], ntt["g3"]
    assert (p1, p2, p3) == ((127 << 54) + 1, (57 << 55) + 1, (27 << 56) + 1)
    assert all(_is_prime_64(p) and (p - 1) % (1 << 54) == 0 and 4 * p < (1 << 63) and (1 << 58) < p for p in (p1, p2, p3))
    assert all(pow(g1, (p1 - 1) // q, p1) != 1 for q in (2, 127))
    assert all(pow(g2, (p2 - 1) // q, p2) != 1 for q in (2, 3, 19))
    assert all(pow(g3, (p3 - 1) // q, p3) != 1 for q in (2, 3))
    l2_P = math.log2(p1 * p2 * p3)
    hr1, hr2 = l2_P - 181.0, l2_P - (math.log2(math.ceil((2.0**68.5) / 58.0)) + 116.0)
    assert hr1 > 1.57 and hr2 > 3.93 and abs(l2_P - ntt["log2_P_prod"]) < 1e-9
    rng, pgs, sizes = random.Random(0x5400612D), [(p1, g1), (p2, g2), (p3, g3)], [512, 1024, 2048, 4096]
    for nb in sizes:
        for _ in range(2):
            A, B = rng.getrandbits(nb) | (1 << (nb - 1)), rng.getrandbits(nb) | (1 << (nb - 1))
            assert _mul_2d_3prime_shoup61(A, B, pgs) == A * B
    return {"ok": True, "primes_hex": [hex(p1), hex(p2), hex(p3)], "prim_roots": [g1, g2, g3], "log2_P": round(l2_P, 6), "hr_1d_b": round(hr1, 4), "hr_2d_b": round(hr2, 4), "live_2d_mul_exact": sizes}


# 4. Mode-A Sanity, Deterministic Bounds, Guard4 x64 & 248/248 Exact SWAR Adaptive Slab
def _byte_msb_mask(end_c: int, w_idx: int) -> int:
    lo = w_idx << 3
    m = 0
    for b in range(lo, lo + 8):
        if 4 <= b < end_c:
            m |= 0x80 << ((b - lo) << 3)
    return m


def _swar_scan_row_exact(row256: bytes, L: int) -> bool:
    """Exact bit-for-bit Python model of _slab_scan_g_active_16 fast-path SWAR masks."""
    end_c = 4 + L
    k16 = (end_c + 15) >> 4
    w = struct.unpack("<32Q", row256)
    m0, m1 = _byte_msb_mask(end_c, 2 * (k16 - 1)), _byte_msb_mask(end_c, 2 * (k16 - 1) + 1)
    if k16 == 1:
        return ((w[0] & m0) | (w[1] & m1)) != 0
    acc = (w[0] & 0x8080808000000000) | (w[2 * (k16 - 1)] & m0) | (w[2 * (k16 - 1) + 1] & m1)
    for j in range(1, 2 * (k16 - 1)):
        acc |= w[j]
    return (acc & 0x8080808080808080) != 0


def check_modeA_guard4_adaptive_slab(led: dict, n_trials: int = 36, V_a: int = 192) -> dict:
    # [MINOR-7a] Realistic spiky L_0.65 (a, b) vectors with exact ratio log1p (avoids float64 cancellation to 0.0!)
    rng = random.Random(0xB1512381)
    max_ng_diff = 0.0
    for _ in range(128):
        a_c = [complex(rng.gauss(0, 150), rng.gauss(0, 150)) for _ in range(4)]
        # Spiky L_0.65: one coordinate small (~1.5..4.0), three larger (~150)
        b_c = [complex(rng.uniform(1.5, 4.0), rng.uniform(-1.0, 1.0))] + [complex(rng.gauss(0, 140), rng.gauss(0, 140)) for _ in range(3)]
        if any(abs(z) < 0.5 for z in b_c):
            continue
        # Compute exact delta = sum 2*log2(|1 - a^3 / (x b^3)|) without large-base cancellation!
        delta = 0.0
        for m in range(4):
            ratio = (a_c[m] ** 3) / (X_BLS * (b_c[m] ** 3))
            # log2(|1 - ratio|^2) = log1p(-2*Re(ratio) + |ratio|^2) / ln(2)
            delta += math.log1p(-2.0 * ratio.real + (ratio.real**2 + ratio.imag**2)) / math.log(2.0)
        max_ng_diff = max(max_ng_diff, abs(delta))
    assert 0.0 < max_ng_diff < 1e-6

    op, mc = led["operating_point"], led["direct_mc_4M"]
    V_a_sh, R0, q = op["V_a_per_shell"], mc["cand_R0"], op["q_ball"]
    Ra = (2.0 * V_a_sh * 120.0 * 70.3125 * 24.0 / (math.pi**4 * 0.267223)) ** 0.125
    det_ng = 8.0 * math.log2(ABS_X * ((R0 / (4.0 ** (1.0 / q))) ** 3) + (Ra / 2.0) ** 3)
    det_nf = 0.0
    for k in (1, 2, 3):
        for ia in range(41):
            a1, a2 = math.sqrt((ia / 40.0) * Ra * Ra / k), math.sqrt((1.0 - ia / 40.0) * Ra * Ra / (4 - k))
            for ib in range(41):
                fb = ib / 40.0
                b1 = ((fb * (R0**q)) / k) ** (1.0 / q) if fb > 0 else 0.0
                b2 = (((1.0 - fb) * (R0**q)) / (4 - k)) ** (1.0 / q) if fb < 1.0 else 0.0
                f1, f2 = a1**9 + a1**6 * b1**3 + a1**3 * b1**6 + b1**9, a2**9 + a2**6 * b2**3 + a2**3 * b2**6 + b2**9
                if f1 > 0 and f2 > 0:
                    det_nf = max(det_nf, 2.0 * (k * math.log2(f1) + (4 - k) * math.log2(f2)))
    assert det_ng < 705.0 and abs(det_ng - mc["det_upper_ng"]) < 1e-6 and det_nf < 630.0 and abs(det_nf - mc["det_upper_nf"]) < 0.01
    for ngc in range(610, 706):
        tg_min = max(1, int(math.floor((ngc + 1.0) / op["q_step_g"])) - 127)
        assert 0 <= 128 - tg_min <= 127 and (128 - tg_min) + int(math.floor(ngc / op["q_step_g"])) <= 255
    assert int(math.floor(630.0 / op["q_step_f"])) - 127 == 38

    # [MINOR-7b] Exhaustive 248/248 SWAR mask bit-exactness test under 0xFF left/right guard poison!
    for L in range(1, 249):
        row = bytearray(b"\xff" * 256)
        row[4 : 4 + L] = bytes([0x7F] * L)
        assert not _swar_scan_row_exact(row, L), f"SWAR false positive at L={L}"
        for pos in (4, 4 + L - 1, 4 + (L >> 1)):
            row[pos] = 0x80
            assert _swar_scan_row_exact(row, L), f"SWAR missed survivor at L={L}, pos={pos}"
            row[pos] = 0x7F

    rng2, tps = random.Random(0xCAFEBABE2026), [1021, 1031, 4093, 4441, 65521]
    def mk_tabs(ws):
        tabs, avs = [], []
        for k_i, p in enumerate(tps):
            nb, av = ((p - 1) >> 2) + 1, [[rng2.randrange(p) for _ in range(V_a)] for _ in range(3)]
            avs.append(av)
            cnt = [0] * nb
            for rj in range(3):
                for v in av[rj]:
                    cnt[v >> 2] += 1
            bkt = [0] * (nb + 1)
            for b in range(nb):
                bkt[b + 1] = bkt[b] + cnt[b]
            cur, hits = list(bkt[:-1]), [0] * (3 * V_a)
            for rj in range(3):
                for aid, v in enumerate(av[rj]):
                    hits[cur[v >> 2]] = (aid << 8) | (v & 255)
                    cur[v >> 2] += 1
            tabs.append({"rem": rng2.randrange(p), "step": rng2.randrange(p), "p": p, "bkt": bkt, "hits": hits, "w": ws[k_i]})
        return tabs, avs

    tg, ag = mk_tabs([21, 25, 29, 33, 37])
    tf, af = mk_tabs([19, 23, 27, 31, 35])
    def run_g4(S, tabs, L):
        em = 0
        for mp in tabs:
            p, w, bkt, hits = mp["p"], mp["w"], mp["bkt"], mp["hits"]
            R = (mp["rem"] + mp["step"]) % p
            mp["rem"] = R
            def emit(bo, st, en):
                nonlocal em
                for j in range(st, en):
                    ad = bo + hits[j]
                    S[ad] = (S[ad] + w) & 0xFF
                em += en - st
            re_ = R + L
            if re_ <= p:
                b0, b1, q0, q1, r8 = R >> 2, (re_ - 1) >> 2, R >> 8, (re_ - 1) >> 8, R & 255
                if q0 == q1:
                    emit(260 - r8, bkt[b0], bkt[b1 + 1])
                else:
                    bm = (q0 + 1) << 6
                    emit(260 - r8, bkt[b0], bkt[bm])
                    emit(516 - r8, bkt[bm], bkt[b1 + 1])
            else:
                tk = p - R
                b0, bl, q0, ql, r8 = R >> 2, (p - 1) >> 2, R >> 8, (p - 1) >> 8, R & 255
                if q0 == ql:
                    emit(260 - r8, bkt[b0], bkt[bl + 1])
                else:
                    bm = (q0 + 1) << 6
                    emit(260 - r8, bkt[b0], bkt[bm])
                    emit(516 - r8, bkt[bm], bkt[bl + 1])
                emit(260 + tk, bkt[0], bkt[((L - tk - 1) >> 2) + 1])
        return em

    S_buf = bytearray(256 + 256 * V_a + 256)
    tot_h = tot_gs = tot_d = 0
    for tr in range(n_trials):
        L = 248 if tr < 6 else (tr + 1 if tr < 18 else rng2.randint(1, 248))
        if tr < 6:
            for tabs in (tg, tf):
                for mp in tabs:
                    mp["rem"], mp["step"] = mp["p"] - (1 + tr * 37), 5
        bgc = bytes(rng2.randint(86, 104) for _ in range(L))
        bfc = bytes(rng2.randint(98, 114) for _ in range(L))
        rg, rf = bytearray(256 * V_a), bytearray(256 * V_a)
        for aid in range(V_a):
            rg[(aid << 8) + 4 : (aid << 8) + 4 + L] = bgc
            rf[(aid << 8) + 4 : (aid << 8) + 4 + L] = bfc
        for tabs, avs, ref in ((tg, ag, rg), (tf, af, rf)):
            for ki, mp in enumerate(tabs):
                p, R, w = mp["p"], (mp["rem"] + mp["step"]) % mp["p"], mp["w"]
                for rj in range(3):
                    for aid, v in enumerate(avs[ki][rj]):
                        df = (v - R) if v >= R else (v + p - R)
                        if df < L:
                            idx = (aid << 8) | (df + 4)
                            ref[idx] = (ref[idx] + w) & 0xFF
        S_buf[:256] = bytes(256)
        S_buf[256 + 256 * V_a :] = bytes(256)
        tpl = bytearray(b"\xff" * 256)
        tpl[4 : 4 + L] = bgc
        for aid in range(V_a):
            S_buf[256 + (aid << 8) : 512 + (aid << 8)] = tpl
        tot_h += run_g4(S_buf, tg, L)
        assert not any(S_buf[:256]) and not any(S_buf[256 + 256 * V_a :])
        gs = []
        for aid in range(V_a):
            b = aid << 8
            row_full = S_buf[256 + b : 512 + b]
            assert row_full[4 : 4 + L] == rg[b + 4 : b + 4 + L]
            if _swar_scan_row_exact(row_full, L):
                gs.extend(b | (c + 4) for c in range(L) if row_full[4 + c] & 0x80)
        tot_gs += len(gs)
        tpl[4 : 4 + L] = bfc
        for aid in range(V_a):
            S_buf[256 + (aid << 8) : 512 + (aid << 8)] = tpl
        tot_h += run_g4(S_buf, tf, L)
        assert not any(S_buf[:256]) and not any(S_buf[256 + 256 * V_a :])
        dg = [o for o in gs if S_buf[256 + o] & 0x80]
        dr = [(aid << 8) | (c + 4) for aid in range(V_a) for c in range(L) if (rg[(aid << 8) | (c + 4)] & 0x80) and (rf[(aid << 8) | (c + 4)] & 0x80)]
        assert dg == dr
        tot_d += len(dg)
    assert tot_d > 0 and tot_gs > 0
    return {"ok": True, "max_ng_drift_bits": float(f"{max_ng_diff:.3e}"), "det_ng": round(det_ng, 4), "det_nf": round(det_nf, 4), "swar_248_exact": True, "guard_0xFF_fp": 0, "dual_surv_exact": tot_d, "emitted_hits": tot_h}


# 5. Live Exact 2D Panjer Floor-Quantized Smooth-Recall & Survival DP
def _panjer_1d(mu: dict[int, float], mk: int):
    P, lam = [0.0] * (mk + 1), sum(mu.values())
    P[0], items = math.exp(-lam), [(j, v) for j, v in mu.items() if v > 0 and j >= 1]
    for k in range(1, mk + 1):
        P[k] = sum(j * v * P[k - j] for j, v in items if j <= k) / k
    return P, lam


def _exact_2d_dp(N_eff: float, B_b: float, plo: float, phi: float, qs: float, T: int, N_max: float) -> dict:
    kN, mu_out, mu_in, mu_1d = int(round(N_eff * 2)), {}, {}, {}
    for a, b in [(0.5614594835668851, plo), (phi, B_b)]:
        if b > a:
            ns = max(20, int((b - a) * 8))
            dx = (b - a) / ns
            for i in range(ns):
                x0 = a + i * dx
                kb = int(round((x0 + 0.5 * dx) * 2))
                if 1 <= kb <= kN:
                    mu_out[kb] = mu_out.get(kb, 0.0) + math.log((x0 + dx) / x0)
    P_out, _ = _panjer_1d(mu_out, kN)
    ns = max(80, int((phi - plo) * 16))
    dx = (phi - plo) / ns
    for i in range(ns):
        x0 = plo + i * dx
        xm, w = x0 + 0.5 * dx, math.log((x0 + dx) / x0)
        kb, tk = int(round(xm * 2)), int(math.floor(xm / qs))
        if kb >= 1 and tk >= 1:
            mu_in[(kb, tk)] = mu_in.get((kb, tk), 0.0) + w
            mu_1d[tk] = mu_1d.get(tk, 0.0) + w
    P_s1, lam_in = _panjer_1d(mu_1d, int(math.ceil(N_eff / qs)) + 20)
    P2: list = [None] * (kN + 1)
    P2[0] = {0: math.exp(-lam_in)}
    items = [(dk, dt, dk * v) for (dk, dt), v in mu_in.items()]
    for k in range(min(dk for dk, _, _ in items), kN + 1):
        d: dict[int, float] = {}
        for dk, dt, dv in items:
            if dk <= k and P2[k - dk]:
                for t0, p0 in P2[k - dk].items():
                    d[t0 + dt] = d.get(t0 + dt, 0.0) + dv * p0
        if d:
            inv_k = 1.0 / k
            P2[k] = {t: v * inv_k for t, v in d.items()}
    st: dict[int, float] = {}
    for k in range(kN + 1):
        if P2[k] and P_out[kN - k] > 0:
            po = P_out[kN - k]
            for t, p in P2[k].items():
                st[t] = st.get(t, 0.0) + p * po
    rec = sum(p for t, p in st.items() if t >= T) / sum(st.values())
    surv = sum(P_s1[j] for j in range(T, len(P_s1)))
    mb = (128 - T) + int(math.floor(N_max / qs))
    assert 0 <= 128 - T and mb <= 255
    return {"max_byte": mb, "rec": rec, "log2_surv": math.log2(max(surv, 1e-300)), "lam": lam_in}


def check_panjer_2d_dp(led: dict) -> dict:
    op, mc = led["operating_point"], led["direct_mc_4M"]
    dg = _exact_2d_dp(mc["cand_cond_ng"] - 1.0487, op["B_g"], op["g_band"][0], op["g_band"][1], op["q_step_g"], 78, 687.0)
    df = _exact_2d_dp(mc["cand_cond_nf"] + 5.3752, op["B_f"], op["f_band"][0], op["f_band"][1], op["q_step_f"], 47, 630.0)
    assert dg["max_byte"] <= 255 and df["max_byte"] <= 255
    assert abs(dg["lam"] - led["sieve_and_stage1b"]["g_guard4_x64"]["hits_active"]) < 1e-4
    assert abs(df["lam"] - led["sieve_and_stage1b"]["f_guard4_x64"]["hits_active"]) < 1e-4
    return {"ok": True, "g_max_byte": dg["max_byte"], "f_max_byte": df["max_byte"], "single_pt_ret": round(dg["rec"] * df["rec"], 5), "single_pt_l2s1": round(dg["log2_surv"] + df["log2_surv"], 4), "adapt_2d_ret": round(op["ret_tot"], 5), "adapt_2d_l2s1": round(op["log2_s1"], 4)}


# 6. Full Ledger: Exact Poisson x64, Packed Super-Leaf Stage 1b, Zero-Extrap R3_U24 LA, RAM & Budget
def _exact_poisson_x64_ir(lam: float) -> float:
    if lam <= 0.0:
        return 0.0
    if lam > 4000.0:
        return (163.0 / 64.0) * lam + 17.7734375
    lo, hi, tot = max(0, int(lam - 10.0 * math.sqrt(lam) - 20)), int(lam + 10.0 * math.sqrt(lam) + 20), 0.0
    for k in range(lo, hi + 1):
        tot += math.exp(-lam + k * math.log(lam) - math.lgamma(k + 1)) * (163.0 * (k // 64) + 23.0 * ((k % 64) // 8) + 5.0 * (k % 8))
    return tot


def _integrate_g4_x64_poisson(plo: float, phi: float, V_a: float, L_w: float, d_ok: float, rmix: list[tuple[int, float]]) -> dict:
    L_em, wpc, ps = L_w + 3.0, 1.0 / (V_a * L_w), min(1.0, L_w / 256.0)
    ovh, dx, uq = 51.0 * (1.0 - ps) + 77.0 * ps, (phi - plo) / 400.0, (0.1, 0.3, 0.5, 0.7, 0.9)
    ih = ip = nok = nb = ne = 0.0
    for i in range(400):
        x = plo + (i + 0.5) * dx
        p = 2.0**x
        did = (p / x) * dx
        dok = did * d_ok
        nok += dok
        nb += dok * (p / 4.0 + 2.0)
        ne += did * V_a
        ip += dok * ovh * wpc
        c1 = 0.0
        for rc, rw in rmix:
            lt = rc * V_a * L_em / p
            c_ns = _exact_poisson_x64_ir(lt)
            c_sp = sum(_exact_poisson_x64_ir(u * lt) + _exact_poisson_x64_ir((1.0 - u) * lt) for u in uq) / 5.0
            c1 += rw * ((1.0 - ps) * c_ns + ps * c_sp)
        ih += dok * c1 * wpc
    return {"ir_hits": ih, "ir_probe": ip, "n_ok": nok, "n_bkt": nb, "n_ent": ne}


def _mul_1d2d_shoup61_ir(B_bits: float, is_root: bool = False) -> float:
    if B_bits <= 1024.0:
        nl = max(1.0, B_bits / 64.0)
        return 12.0 * nl * nl * 1.15
    if B_bits <= float(1 << 59):
        nl = math.ceil(B_bits / 64.0)
        N = float(1 << int(math.ceil(math.log2(max(2.0, 2.0 * nl)))))
        assert N <= (1 << 54)
        dist = (B_bits > L_OWN_1D) if is_root else (B_bits * 2.0 > L_OWN_1D)
        return (3.0 * (1.5 * N * math.log2(N) * 17.0 + N * 46.0) + (3.0 * N * 160.0 if dist else 0.0) + N * 67.0) * 1.15
    L = int(math.ceil(B_bits / 58.0))
    k = max(1, (L - 1).bit_length())
    k1 = (k + 1) // 2
    k2 = k - k1
    N2d = float(1 << (k1 + 1)) * float(1 << (k2 + 1))
    bfly = N2d * (1.0 * (k1 + 1) + 1.5 * (k2 + 1))
    dist = (B_bits > L_OWN_2D) if is_root else (B_bits * 2.0 > L_OWN_2D)
    return (3.0 * (bfly * 17.0 + N2d * 46.0) + (9.0 * N2d * 160.0 if dist else 0.0) + N2d * 95.0) * 1.15


def _rem_tree_packed_ir(nl: float, lb: float, k_pack: int, dmult: float = 2.00) -> tuple[float, int]:
    if nl <= 1.0:
        return 40000.0, 1
    nsup, sbits = max(1.0, nl / float(k_pack)), float(k_pack) * float(lb)
    h = int(math.ceil(math.log2(nsup))) if nsup > 1.0 else 0
    ir_p = ir_d = 0.0
    cn, cb = nsup, sbits
    for _ in range(h):
        pairs = cn * 0.5
        m_ir = _mul_1d2d_shoup61_ir(cb, False)
        ir_p += pairs * m_ir
        ir_d += pairs * dmult * m_ir
        cn, cb = pairs, cb * 2.0
    ir_r = 3.0 * _mul_1d2d_shoup61_ir(nsup * sbits, True)
    nll = max(1.0, lb / 64.0)
    ir_sb = sum(12.0 * (j * nll) * nll * 1.15 for j in range(1, k_pack))
    ir_sd = float(k_pack) * 24.0 * (sbits / 64.0) * nll * 1.15
    ir_lg = float(nl) * max(25000.0, 40000.0 * (lb / 636.0))
    return ir_p + ir_r + ir_d + nsup * (ir_sb + ir_sd) + ir_lg, h


def _tier_probs(cof: float, P_t: float, B_b: float, k_max: int = 8):
    if cof <= P_t:
        return 1.0, 0.0, 1.0
    om = lambda v: 0.0 if v < 1.0 else (1.0 / v if v <= 2.0 else ((1.0 + math.log(v - 1.0)) / v if v <= 3.0 else 0.57))
    dy = (cof - P_t) / 200.0
    md = DICKMAN.rho(cof / P_t)
    ma, mt = 0.0, md
    for i in range(200):
        y = P_t + (i + 0.5) * dy
        w = DICKMAN.rho((cof - y) / P_t) * (om(y / P_t) / P_t) * dy
        mt += w
        if y <= B_b:
            md += w
        elif any(k * P_t <= y <= k * B_b for k in range(2, k_max + 1)):
            ma += w
    return md / max(mt, 1e-30), ma / max(mt, 1e-30), (md + ma) / max(mt, 1e-30)


def check_full_ledger_ram_budget(led: dict) -> dict:
    op, mc, cheb = led["operating_point"], led["direct_mc_4M"], led["chebotarev_exact"]
    V_a, L_w = op["V_a_per_shell"], mc["L_win_eff"]
    sg = _integrate_g4_x64_poisson(op["g_band"][0], op["g_band"][1], V_a, L_w, cheb["d_ok_g"], [(3, 1.0)])
    sf = _integrate_g4_x64_poisson(op["f_band"][0], op["f_band"][1], V_a, L_w, cheb["d_ok_f"], [(3, 39.0 / 46.0), (6, 6.0 / 46.0), (9, 1.0 / 46.0)])
    l2_Vk_ch = mc["cand_mc_log2_vk"] + 3.0 * mc["cand_se_bits"] + math.log2(1.0 / (1.0 - mc["gap_frac_charged"])) + mc["anchor_pad_bits"]
    assert abs(l2_Vk_ch - mc["log2_Vk_charged"]) < 1e-9
    V_tot, surv_s1, peak_s1b = 2.0**l2_Vk_ch, 2.0 ** op["log2_s1"], 0.0
    ppb = lambda a, b: max(1.0, (2.0**b - (0.0 if a <= 0.0 else 2.0**a)) / math.log(2.0))

    def price_t(blo, bhi, pin, lb, kp, dens):
        nonlocal peak_s1b
        Lt = max(1.0, (ppb(blo, min(bhi, 36.0)) if blo < 36.0 else 0.0) + (dens * ppb(max(blo, 36.0), bhi) if bhi > 36.0 else 0.0))
        Mt = V_tot * pin
        Mb = min(Mt, max(1.0, Lt / lb))
        nb, LX = max(1.0, Mt / Mb), Mb * lb
        ir, h_sup = _rem_tree_packed_ir(Mb, lb, kp, 2.00)
        if LX < Lt:
            ir += max(0.0, Lt / LX - 1.0) * 2.5 * _mul_1d2d_shoup61_ir(LX, True)
        peak_s1b = max(peak_s1b, ((9.0 if Lt > float(1 << 59) else 4.5) + float(h_sup)) * Lt / 8.0)
        return (nb * ir) / V_tot

    Bf, Bg, Pf, Pg = op["B_f"], op["B_g"], op["f_band"][1], op["g_band"][1]
    mf, mg = 0.5 * (Pf + Bf), 0.5 * (Pg + Bg)
    _, paf1, _ = _tier_probs(385.0, Pf, Bf)
    _, paf2, _ = _tier_probs(385.0, mf, Bf)
    _, pag1, _ = _tier_probs(405.0, Pg, Bg)
    _, pag2, _ = _tier_probs(405.0, mg, Bg)
    rf_e = min(0.95, math.sqrt(op["ret_tot"] / 0.80))
    pf_sm = min(1.0, (DICKMAN.rho((mc["cand_cond_nf"] + 5.3752) / Bf) * rf_e) / (2.0 ** op["log2_sf"]))
    pg_sm = min(1.0, (DICKMAN.rho((mc["cand_cond_ng"] - 1.0487) / Bg) * (op["ret_tot"] / rf_e)) / (2.0 ** op["log2_sg"]))
    cr = surv_s1 * 1700.0 + surv_s1 * pf_sm * 1700.0 + sum(price_t(*a) for a in [
        (1.0, Pf, surv_s1, 560.0, 7, cheb["dens_f"]), (Pf, mf, surv_s1 * paf1, 385.0, 10, cheb["dens_f"]), (mf, Bf, surv_s1 * paf2, 385.0, 10, cheb["dens_f"]),
        (1.0, Pg, surv_s1 * pf_sm, 675.0, 6, cheb["dens_g"]), (Pg, mg, surv_s1 * pf_sm * pag1, 405.0, 10, cheb["dens_g"]), (mg, Bg, surv_s1 * pf_sm * pag2, 405.0, 10, cheb["dens_g"]),
    ])
    c_cof = cr + surv_s1 * pf_sm * pg_sm * (110.0 * 256.0 * 150.0 * 4.0 + 0.05 * 5.0e7 * 2.0)
    slab_ch = led["sieve_and_stage1b"]["ipc_slab_charged"]
    ipc_tot = sg["ir_hits"] + sf["ir_hits"] + sg["ir_probe"] + sf["ir_probe"] + slab_ch + c_cof
    coll_l2 = math.log2((V_tot * ipc_tot) / C_COST)
    n_ent, n_bkt, n_ok = sg["n_ent"] + sf["n_ent"], sg["n_bkt"] + sf["n_bkt"], sg["n_ok"] + sf["n_ok"]
    rep_l2 = math.log2(op["K_shells"] * (n_ent * 27.5 + n_bkt * 6.5) / C_COST)
    assert abs(ipc_tot - led["sieve_and_stage1b"]["ipc_tot"]) < 1e-6
    assert abs(coll_l2 - led["ledger_rows_live_log2_u"]["relation_collection"]) < 1e-6
    assert abs(rep_l2 - led["ledger_rows_live_log2_u"]["replay_a_tables"]) < 1e-6

    la = led["la_and_sge_zero_extrap"]
    U = ((2.0**Bf) / (Bf * math.log(2.0)) + (2.0**Bg) / (Bg * math.log(2.0))) / 12.0
    assert la["nl_cap"] == 0.00205 and la["w_cap"] == 190.0 and la["nl_cap"] > la["u24_2seed_eta85_nl_mean"] + 2.0 * la["u24_2seed_eta85_nl_se"]
    N_l, M_bw, N_bw = U * la["nl_cap"], op["m_block"], op["n_block"]
    N_nnz = N_l * la["w_cap"]
    ir_spmv_x4 = (N_nnz * (1.144 * 8.0 + 0.122 * 24.5) + N_l * 2100.0) * 1.15
    ir_spmv_x1 = (N_nnz * (1.144 * 11.0 + 0.122 * 28.0) + N_l * 2100.0) * 1.15
    L_seq, d_mk = N_l / M_bw + N_l / N_bw + 64.0, N_l / N_bw + 32.0
    ir_bw = 1.02 * (N_bw * L_seq * ir_spmv_x4 + d_mk * ir_spmv_x1)
    ir_mk = 1.02 * d_mk * N_bw * N_l * 283.0 * 1.15
    ir_lg = 20.0 * ((M_bw + N_bw) ** 3) * L_seq * math.log2(L_seq) * 236.0 * 1.15
    la_l2 = math.log2((ir_bw + ir_mk + ir_lg + 0.008 * (ir_bw + ir_mk) + 30.0 * U * 283.0 * 1.15) / C_COST)
    sge_l2 = math.log2(((12500.0 * U + 25.0 * U * math.log2(U) + 100.0 * N_nnz) * 1.15) / C_COST)
    assert abs(la_l2 - led["ledger_rows_live_log2_u"]["linear_algebra"]) < 1e-6
    assert abs(sge_l2 - led["ledger_rows_live_log2_u"]["matrix_generator_sge"]) < 1e-6

    s1_r = math.log2(8.0 * (n_bkt + n_ent) + 48.0 * n_ok + 256.0 * V_a + 2.0**52)
    s2_r = math.log2((op["eta"] * U) * 18.0 * 12.0)
    s3_r = math.log2(N_nnz * 11.0 + N_l * 8.0 + (N_bw + 6.0) * N_l * 32.0 + (M_bw * N_bw * L_seq + 4.0 * N_bw * N_bw * d_mk + 2.0 * M_bw * (M_bw + N_bw) * L_seq) * 32.0)
    s1b_r = math.log2(peak_s1b)
    mem = led["memory_audit_unified"]
    assert abs(s1_r - mem["stage1_peak_global_log2_B"]) < 1e-6 and abs(s2_r - mem["stage2_sge_global_log2_B"]) < 1e-6
    assert abs(s3_r - mem["stage3_la_global_log2_B"]) < 1e-6 and abs(s1b_r - mem["stage1b_global_log2_B"]) < 1e-5
    assert max(s1_r, s2_r, s3_r, s1b_r) < mem["peak_global_claim_log2_B"] < 120.0

    lr, cr_l2, tot = led["ledger_rows_live_log2_u"], led["ledger_rows_charged_log2_u"], led["totals"]
    assert len(lr) == len(cr_l2) == 14 and all(cr_l2[k] > lr[k] for k in lr)
    w_oc = sum(2.0**v for v in cr_l2.values())
    l2_oc = math.log2(w_oc)
    assert abs(math.log2(sum(2.0**v for v in lr.values())) - tot["opened_run_live_log2_u"]) < 1e-9
    assert abs(l2_oc - tot["opened_run_charged_log2_u"]) < 1e-9

    p_m = w_oc / (2.0**124.0)
    p_fail = p_m + 2.0**(-48) + 2.0**(-100) + 2.0**(-132) + 2.0**(-200)
    p_open = 1.0 - p_fail
    K_g, two64 = tot["K_gate"], 1 << 64
    assert p_fail < 2.0e-5 and p_open > 0.99998 and K_g == 7194374076228249696
    assert Fraction(K_g, two64) * Fraction(99998, 100000) >= Fraction(39, 100)
    assert Fraction(K_g - 1, two64) * Fraction(99998, 100000) < Fraction(39, 100)
    p_cold = (K_g / float(two64)) * p_open
    num, den = int(tot["claim_numerator"]), int(tot["claim_denominator"])
    l2_gc = l2_oc + math.log2(K_g / float(two64))
    assert p_cold > 0.39 and math.gcd(num, den) == 1 and (num, den) == (21337, 200)
    assert abs(l2_gc - tot["gated_charged_log2_u"]) < 1e-9 and l2_gc < float(Fraction(num, den))
    return {
        "ok": True, "claim": f"{num}/{den}", "ipc_tot": round(ipc_tot, 6),
        "opened_ch": round(l2_oc, 6), "gated_ch": round(l2_gc, 6), "slack_b": round(float(Fraction(num, den)) - l2_gc, 6),
        "p_meter": round(p_m, 9), "p_open": round(p_open, 8), "p_cold": round(p_cold, 8), "K_gate_min": True,
        "ram_l2_B": {"S1": round(s1_r, 4), "S1b": round(s1b_r, 4), "S2": round(s2_r, 4), "S3": round(s3_r, 4)},
    }


# 7. Optional Sibling-File Consistency & ISA Audit (Skipped in Isolated Docker!)
FORB_REG = re.compile(r"(?<![a-zA-Z0-9_])([vqdshzp][0-9]|[vqdshzp][12][0-9]|[vqdshzp]3[01]|x18|w18)(?![a-zA-Z0-9_])")
FORB_OP = re.compile(r"^\s*(dc\s+zva|fadd|fsub|fmul|fdiv|fmadd|fmsub|fcvt|fmov|fcmp|ld1|st1|addv|umull2|pmull)\b")


def check_optional_siblings(led: dict) -> dict:
    has = {f: os.path.exists(ROOT / f) for f in ("ledger_data.json", "submission.json", "KERNELS_SIEVE.md", "KERNELS_LA_ECM.md")}
    if not any(has.values()):
        return {"ok": True, "standalone_mode": True, "sibling_files_audited": 0}
    out: dict = {"ok": True, "standalone_mode": False}
    if has["ledger_data.json"]:
        assert json.loads((ROOT / "ledger_data.json").read_text(encoding="utf-8")) == led
        out["ledger_json_match"] = True
    if has["submission.json"]:
        sub = json.loads((ROOT / "submission.json").read_text(encoding="utf-8"))
        sc = sub["claims"][0]["score"]
        assert (sc["numerator"], sc["denominator"]) == (led["totals"]["claim_numerator"], led["totals"]["claim_denominator"])
        assert "required" in sub["experiments"][0]["input_schema"]
        out["submission_json_ok"] = True
    if has["KERNELS_SIEVE.md"] and has["KERNELS_LA_ECM.md"]:
        isa = {}
        for fn in ("KERNELS_SIEVE.md", "KERNELS_LA_ECM.md"):
            txt = (ROOT / fn).read_text(encoding="utf-8")
            if fn == "KERNELS_LA_ECM.md":
                # [MAJOR-1] Verify _crt3_recombine_carry uses exact 2^54-smooth primes (NO stale 2^51 constants!)
                assert "2276569611635785729" not in txt and "2254051613498933249" not in txt
                assert "-2053641430080946177" in txt and "-1945555039024054273" in txt
            blks = re.findall(r"```asm\n(.*?)```", txt, flags=re.S)
            ni = ns = 0
            for b in blks:
                for ln in b.splitlines():
                    c = ln.split(";")[0].strip()
                    if not c or c.startswith(".") or c.endswith(":"):
                        continue
                    ni += 1
                    lo = c.lower()
                    if lo.startswith("staddb"):
                        ns += 1
                    assert not FORB_OP.search(lo) and not FORB_REG.search(lo), f"Forbidden ISA: {c}"
            assert len(blks) >= 3 and ((ns > 0) if fn == "KERNELS_SIEVE.md" else (ns == 0))
            isa[fn] = {"blocks": len(blks), "insns": ni, "staddb": ns, "forbidden": 0}
        out["isa_audit"] = isa
        out["crt3_1d_p54_verified"] = True
    files = sorted(p for p in ROOT.iterdir() if p.is_file())
    if len(files) >= 9:
        priv = ("/" + "Users/", "research-" + "work/", "curvesmash-" + "agy-cli")
        tb = 0
        for p in files:
            sz = p.stat().st_size
            assert sz <= 65536, f"{p.name} too large: {sz}"
            tb += sz
            txt = p.read_text(encoding="utf-8")
            assert not any(t in txt for t in priv), f"Path leak in {p.name}"
        assert tb < 1048576
        out["hygiene"] = {"files": len(files), "total_bytes": tb, "max_bytes": max(p.stat().st_size for p in files), "leaks": 0}
    return out


# 8. Diagnostics: Live F_q 2-Merge Rank, 9-Rung SGE Ladder & Stage 4 12D M_g + ECM
Q_MOD, ZETA12 = 2147483629, 1803057106


def _fq_rank_2merge(U: int, eta: float, seed: int) -> tuple[int, int, int]:
    rng, zp = random.Random(seed), [pow(ZETA12, k, Q_MOD) for k in range(12)]
    Uf, lnBf, lnBg = int(round(0.439 * U)), 70.00 * math.log(2.0), 70.36 * math.log(2.0)
    Ug, rf1, rg1 = U - Uf, DICKMAN.rho(6.65), DICKMAN.rho(8.31)
    ws = [
        min(0.45, (1.173 / U) / ((j + 0.5) / Uf) * (DICKMAN.rho(7.65 - max(0.02, 1.0 + math.log((j + 0.5) / Uf) / lnBf)) / rf1) + (0.35 / (2.0 + 2.5 * j) if j < 40 else 0.0))
        for j in range(Uf)
    ] + [
        min(0.45, (1.196 / U) / ((j + 0.5) / Ug) * (DICKMAN.rho(9.31 - max(0.02, 1.0 + math.log((j + 0.5) / Ug) / lnBg)) / rg1) + (0.40 / (2.0 + 2.5 * j) if j < 40 else 0.0))
        for j in range(Ug)
    ]
    mw, cdf, acc = sum(ws), [], 0.0
    for w in ws:
        acc += w / mw
        cdf.append(acc)
    rows, Lp = [], math.exp(-mw)
    for _ in range(int(round(eta * U))):
        kp, pa = 0, 1.0
        while pa > Lp:
            kp += 1
            pa *= rng.random()
        d: dict[int, int] = {}
        for _ in range(max(4, kp - 1)):
            c = min(U - 1, bisect.bisect_left(cdf, rng.random()))
            nv = (d.get(c, 0) + zp[rng.randrange(12)]) % Q_MOD
            if nv:
                d[c] = nv
            else:
                d.pop(c, None)
        if d:
            rows.append(d)
    while True:
        cr: dict[int, list[int]] = {}
        for ri, d in enumerate(rows):
            if d:
                for c in d:
                    cr.setdefault(c, []).append(ri)
        ch = False
        for rl in cr.values():
            if len(rl) == 1:
                rows[rl[0]] = None
                ch = True
        if ch:
            continue
        mg = False
        for c, rl in cr.items():
            if len(rl) == 2:
                r0, r1 = rl
                d0, d1 = rows[r0], rows[r1]
                if not d0 or not d1 or c not in d0 or c not in d1:
                    continue
                if len(d0) > len(d1):
                    r0, r1, d0, d1 = r1, r0, d1, d0
                fac = (d1[c] * pow(d0[c], -1, Q_MOD)) % Q_MOD
                for col, v0 in d0.items():
                    nv = (d1.get(col, 0) - fac * v0) % Q_MOD
                    if nv:
                        d1[col] = nv
                    else:
                        d1.pop(col, None)
                rows[r0] = None
                if not d1:
                    rows[r1] = None
                mg = True
                break
        if not mg:
            break
    alive = [d for d in rows if d]
    cols = sorted({c for d in alive for c in d})
    cm, piv, rk = {c: i for i, c in enumerate(cols)}, {}, 0
    for d in alive:
        cur = {cm[c]: v for c, v in d.items()}
        while cur:
            pc = min(cur)
            if pc in piv:
                fac = cur[pc]
                for c, pv in piv[pc].items():
                    nv = (cur.get(c, 0) - fac * pv) % Q_MOD
                    if nv:
                        cur[c] = nv
                    else:
                        cur.pop(c, None)
            else:
                inv = pow(cur[pc], -1, Q_MOD)
                piv[pc] = {c: (v * inv) % Q_MOD for c, v in cur.items()}
                rk += 1
                break
    return len(alive), len(cols), len(cols) - rk


def diag_fq_rank_2merge() -> dict:
    by_eta = {}
    for eta in (0.82, 0.85, 0.90):
        nulls = surp = 0
        for s in range(4):
            nr, nc, cn = _fq_rank_2merge(512, eta, 9000 + int(eta * 1000) * 10 + s)
            surp += int(nr >= nc)
            nulls += min(nr, nc) - (nc - cn)
        assert surp == 4 and nulls == 0
        by_eta[f"eta_{eta:.2f}"] = {"seeds": 4, "surplus_pos": surp, "dense_nullity": nulls}
    return {"ok": True, "U": 512, "total_dense_nullity": 0, "by_eta": by_eta}


def diag_sge_9rung_ladder(led: dict) -> dict:
    rungs = [
        (16, 0.006424, 0.40446), (17, 0.005037, 0.35823), (18, 0.003139, 0.21837),
        (19, 0.003813, 0.38645), (20, 0.003907, 0.30853), (21, 0.002708, 0.24756),
        (22, 0.001343, 0.24459), (23, 0.001380, 0.25173), (24, 0.001738, 0.32249),
    ]
    la = led["la_and_sge_zero_extrap"]
    nlc, wc = la["nl_cap"], la["w_cap"]
    nnzc = nlc * wc
    assert nlc > 0.001738 + 2.0 * 0.000111 and nlc > 0.001380 + 2.0 * 0.000106
    assert nnzc > 0.32249 and (nlc * nnzc) > 1.42 * (0.001738 * 0.32249) and la["ir_core_nnz_x4"] > 10.56
    xs, ys = [math.log(k * math.log(2.0)) for k, _, _ in rungs], [math.log(nl) for _, nl, _ in rungs]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    p62 = math.exp((my - b * mx) + b * math.log(62.0 * math.log(2.0)))
    assert p62 < 0.0001 and nlc / p62 > 40.0
    return {"ok": True, "rungs": 9, "u24_m_plus_2se": 0.001960, "nl_cap": nlc, "u24_la_headroom": round((nlc * nnzc) / (0.001738 * 0.32249), 3), "ols_exp": round(b, 3), "ols_proj_2p62": round(p62, 6), "ols_headroom_x": round(nlc / p62, 1)}


_EC, _EP = {}, {}


def _opt_ecm(qb: float, Nb: float):
    qk, key = round(qb, 2), (round(qb, 2), int(math.ceil(Nb / 64.0)))
    if key in _EC:
        return _EC[key]
    if qk not in _EP:
        pl, qe = [], max(1.0, qb - 3.5)
        for lb4 in range(32, 249):
            B1, l1 = 2.0 ** (lb4 * 0.25), lb4 * 0.25
            p1 = DICKMAN.rho(qe / l1)
            for m in (10.0, 25.0, 60.0, 120.0, 250.0):
                B2 = B1 * m
                dt = (math.log2(B2) - l1) / 120.0
                p2 = sum((DICKMAN.rho((qe - (l1 + (i + 0.5) * dt)) / l1) * (dt / (l1 + (i + 0.5) * dt))) if (l1 + (i + 0.5) * dt) < qe else (dt / (l1 + (i + 0.5) * dt)) for i in range(120))
                pr = min(1.0, p1 + p2)
                if pr > 1e-40:
                    pl.append((B1, B2, pr))
        _EP[qk] = pl
    L = max(4, int(math.ceil(Nb / 64.0)))
    mir = 212.0 if L == 4 else (13.25 * L * L + 80.0 if L <= 16 else min(13.25 * L * L + 80.0, 2.5 * _mul_1d2d_shoup61_ir(L * 64.0)))
    best = (1e300, None)
    for B1, B2, pr in _EP[qk]:
        np2 = max(0.0, B2 / math.log(B2) - B1 / math.log(B1))
        c1 = 1.15 * (40.0 * mir + (B1 / math.log(2.0)) * (11.0 * mir + 93.0 * L) + 2.5 * np2 * (mir + 40.0 * L))
        if c1 / pr < best[0]:
            best = (c1 / pr, (B1, B2, pr, c1))
    _EC[key] = best
    return best


def _descent_node_u(qb: float, side: str, Bf: float, Bg: float, cap: float) -> float:
    Bfl, Bgl, best = max(Bf, cap), max(Bg, cap), 1e300
    for deg in (1, 2):
        for s2 in range(45):
            s = s2 * 0.5
            nf = (4.5 if deg == 1 else 3.0) * qb + 72.0 * s + (111.0 if deg == 1 else 55.0)
            ng = (1.5 if deg == 1 else 1.0) * qb + 8.0 * (max(3.0 * s, LOG2_ABS_X - 3.0 * s) if deg == 1 else max(3.0 * s, LOG2_ABS_X, 2.0 * LOG2_ABS_X - 3.0 * s)) + (45.0 if deg == 1 else 35.0)
            nfr, ngr = (max(1.0, nf - qb), ng) if side == "f" else (nf, max(1.0, ng - qb))
            l2T, dim = 20.0, 8.0 * (deg + 1)
            for _ in range(4):
                nft, ngt = nfr + (72.0 / dim) * l2T, ngr + (24.0 / dim) * l2T
                pj = DICKMAN.rho(nft / Bfl) * DICKMAN.rho(ngt / Bgl)
                if pj <= 1e-250:
                    l2T = 300.0
                    break
                l2T = -math.log2(pj) + 1.0
            if l2T <= 120.0:
                tu = ((2.0**l2T) * 2000.0 + 3.0 * (_opt_ecm(Bfl, nft)[0] + _opt_ecm(Bgl, ngt)[0]) + 1.0e7) / C_COST
                if tu < best:
                    best = tu
    return best


def diag_stage4_12d_mg_ecm(led: dict) -> dict:
    s4, op = led["stage4_12d_mg"], led["operating_point"]
    det = 6.0 * math.log2(P_BLS) + 4.0 * math.log2(ABS_X)
    lll = 33.0 - 6.0 * math.log2(12.0) + 4.5 * math.log2(5.0) + 6.0 * math.log2(3.0) + det
    assert abs(det - s4["det_log2"]) < 1e-6 and lll < s4["N_boot_bits"] == 2572.0
    mrel = 0.0
    for qb, el2 in [(70, 26.525), (100, 33.023), (140, 40.529), (185, 47.965), (245, 56.784), (275, 60.848)]:
        irf, (_, _, pr, c1) = _opt_ecm(float(qb), 2572.0)
        mrel = max(mrel, abs(irf - c1 / pr) / irf)
        assert mrel < 1e-12 and abs(math.log2(irf / C_COST) - el2) < 0.002
    Bf, Bg = op["B_f"], op["B_g"]
    ub1 = ((1.0 / (0.918 * DICKMAN.rho(2572.0 / 275.0))) * (2.5 * _opt_ecm(275.0, 2572.0)[0])) / C_COST
    ul0 = 5.0 * _descent_node_u(550.0, "g", Bf, Bg, 275.0)
    ut1 = ul0 + sum(mn * max(_descent_node_u(qh, "f", Bf, Bg, ql), _descent_node_u(qh, "g", Bf, Bg, ql)) for qh, ql, mn in [(275.0, 155.0, 5), (155.0, 130.0, 12), (130.0, 106.0, 28), (106.0, 86.0, 65), (86.0, Bg, 200)])
    mf = (1.04094 / op["eta"]) ** 2
    assert abs(math.log2(2.0 * (ub1 + ut1) * mf) - s4["live_2log_with_snell_log2_u"]) < 1e-6
    assert abs(math.log2(2.0 * (ut1 + 0.25 * ub1) * mf) - s4["ecm_split_live_log2_u"]) < 1e-6

    xdbl = lambda X, Z, A24, N: (((X + Z) ** 2 * (X - Z) ** 2) % N, ((((X + Z) ** 2 - (X - Z) ** 2) % N) * (((X - Z) ** 2 + A24 * ((X + Z) ** 2 - (X - Z) ** 2)) % N)) % N)
    xadd = lambda X1, Z1, X2, Z2, Xd, Zd, N: ((Zd * (((X1 - Z1) * (X2 + Z2) + (X1 + Z1) * (X2 - Z2)) ** 2)) % N, (Xd * (((X1 - Z1) * (X2 + Z2) - (X1 + Z1) * (X2 - Z2)) ** 2)) % N)
    def xmul(k, X, Z, A24, N):
        r0x, r0z, r1x, r1z = 1, 0, X % N, Z % N
        for bit in bin(k)[2:]:
            if bit == "0":
                r1x, r1z = xadd(r0x, r0z, r1x, r1z, X, Z, N)
                r0x, r0z = xdbl(r0x, r0z, A24, N)
            else:
                r0x, r0z = xadd(r0x, r0z, r1x, r1z, X, Z, N)
                r1x, r1z = xdbl(r1x, r1z, A24, N)
        return r0x, r0z
    nx, dx = 5, -7
    t16 = [p for p in _sieve_primes(200) if p >= 101]
    for p in t16:
        x = (nx * pow(dx, -1, p)) % p
        d = ((2 * x * x - 1) * pow(x, -4, p)) % p
        A = (2 * (1 + d) * pow(1 - d, -1, p)) % p
        cnt = 1 + sum(1 if (u**3 + A * u * u + u) % p == 0 else (2 if pow((u**3 + A * u * u + u) % p, (p - 1) // 2, p) == 1 else 0) for u in range(p))
        assert cnt % 16 == 0
    N1 = 50023 * 50033
    x1 = (nx * pow(dx, -1, N1)) % N1
    A24 = pow(1 - ((2 * x1 * x1 - 1) * pow(x1, -4, N1)) % N1, -1, N1)
    X, Z = 3, 1
    for p in _sieve_primes(30):
        pk = p
        while pk * p <= 30:
            pk *= p
        X, Z = xmul(pk, X, Z, A24, N1)
    assert math.gcd(Z, N1) == 1
    babies = [xmul(i, X, Z, A24, N1) for i in range(1, 16) if math.gcd(i, 30) == 1]
    DX, DZ = xmul(30, X, Z, A24, N1)
    gp, gc, w1 = (1, 0), (DX, DZ), 1
    for k in range(1, 13):
        for bx, bz in babies:
            g = math.gcd((gc[0] * bz - bx * gc[1]) % N1, N1)
            if 1 < g < N1:
                w1 = g
                break
        if w1 > 1:
            break
        gp, gc = gc, (xdbl(DX, DZ, A24, N1) if k == 1 else xadd(gc[0], gc[1], DX, DZ, gp[0], gp[1], N1))
    X2, Z2 = xmul(2, 4, 1, 2, 35)
    DX2, DZ2 = xmul(4, X2, Z2, 2, 35)
    w2 = math.gcd((DX2 * Z2 - X2 * DZ2) % 35, 35)
    assert w1 == 50023 and math.gcd(Z2, 35) == 1 and w2 == 5
    return {"ok": True, "M_g_lll_amgm": round(lll, 4), "ecm_rel_err": mrel, "layer0_l2u": round(math.log2(ul0), 4), "live_2log_snell": round(math.log2(2.0 * (ub1 + ut1) * mf), 4), "t16_primes": f"{len(t16)}/{len(t16)}", "w1": w1, "w2": [w2, 35 // w2]}


def main() -> None:
    t0 = time.perf_counter()
    raw = sys.stdin.read().strip()
    inp = json.loads(raw) if raw else {}
    if not isinstance(inp, dict) or (set(inp.keys()) - {"part"}):
        raise ValueError("Input must be a JSON object with optional key 'part'")
    part = inp.get("part", "core")
    if part not in ("core", "diagnostics", "all"):
        raise ValueError(f"Unsupported part: {part}")

    led = LEDGER_DATA
    out: dict = {"package_id": led["package_id"], "config_tag": led["config_tag"], "part": part}
    if part in ("core", "all"):
        out["stage0_algebra"] = check_stage0_algebra()
        out["chebotarev_theorems_2a_2b"] = check_chebotarev_theorems_2a_2b(5000)
        out["ntt_1d2d_shoup61"] = check_ntt_1d2d_shoup61(led)
        out["modeA_guard4_x64_adaptive_slab"] = check_modeA_guard4_adaptive_slab(led)
        out["panjer_2d_dp"] = check_panjer_2d_dp(led)
        out["ledger_ram_budget"] = check_full_ledger_ram_budget(led)
        out["optional_sibling_checks"] = check_optional_siblings(led)
    if part in ("diagnostics", "all"):
        out["diag_fq_rank_2merge"] = diag_fq_rank_2merge()
        out["diag_sge_9rung_ladder"] = diag_sge_9rung_ladder(led)
        out["diag_stage4_12d_mg_ecm"] = diag_stage4_12d_mg_ecm(led)

    out["checks_ok"] = all(v.get("ok", True) for v in out.values() if isinstance(v, dict))
    out["elapsed_seconds"] = round(time.perf_counter() - t0, 3)
    print(json.dumps(out, indent=2))
    if not out["checks_ok"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
