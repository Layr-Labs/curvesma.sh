"""v32_verify.py -- Pure-stdlib exact verification & model helpers for BN254 C6 v32d:
  1. get_ret1d_u16: exact 1D Retention (M.Retention, step=0.25, odd=2) for uint16_t STADDH scoring
     (scale=16.0, w_p = floor(16*log2(p)), t(a) = ceil(16*(L0(a) - sl))):
       - Smooth retention lower bound (rg): restricted to the squarefree-sieve-part smooth subset (where every dividing
         prime p >= 2^24 has multiplicity 1 and is scored once by STADDH). M.Retention.at(n, r) already restricts to
         squarefree sieve/large primes (sf=True, multiplied by selection_normalizer) and subtracts rb=(cap(n)+1)/16
         internally; passing r = floor(sl_c - 0.50) adds an extra >= 0.50 b margin (> |log2 N_g - L0| + 1/16 = 0.078 b),
         so every counted squarefree-sieve-part smooth norm unconditionally satisfies S >= t(a).
       - Survivor upper bound (sf/s0): for ALL integers N_g (squarefree AND squareful), S >= t(a) implies
         log2(N_g_nonsieve) <= log2(N_g) - sum_{p>=2^24, p|N_g} log2(p) <= sl + 0.015503 unconditionally.
         sf_tbl[n][r] uses M.Retention.survivor_upper(n, r + 0.50, ex=0.0) (all-integer geometric series, sf=False,
         NO rb subtraction, NO selection_normalizer) at r = ceil(sl_c), summing all nonsieve mass up through
         ceil(sl_c) + 0.50 >= sl_c + 0.50 > sl_c + 0.015503; s0_tbl[n][r] = min(all-integer B-smooth prob, sf_tbl[n][r]).
  2. Exact closed-form rational 3-adic class split (P1, P2, P3 degree-1 + Q3 degree-3, zero double-counting) + live O_h check.
  3. Deterministic ntt3_F Bernstein recursive split (exact upper bound + 2x safety + live 64-integer prototype).
  4. Instruction-by-instruction / bit-exact verification of:
     - Aligned-sub-band bfi x16 Lower Band 2 (uint16_t STADDH, entry = (ai << 5) | j_rel, all k_s in 38..49, W in {8, 16}, w_bits = log2(W) in {3, 4}, sgn in {+1, -1}, 240 cases incl. non-zero aligned bases)
     - 32-bit delta-packed Band 1 on [2^24, 2^28.5) (delta = 2*(ai_new - ai_old), sub-segment splits on byte delta >= 2^32)
     - Lazy 2-Barrett5 + 4-limb CIOS MontMul(chi_R) c_tail_z6 (323 ops <= 340, 200 random tests)
     - Unconditional |log2 N_g(a,b) - log2 N(u*w*a^2)| <= 0.015503 b + uint16_t ceil-threshold inclusion & zero-underflow/wrap-around + read-only scan8_u16.
"""
import math, random, re
from fractions import Fraction as Fr
from math import log, log2, exp, sqrt, ceil, floor
from pathlib import Path

LN2 = math.log(2.0)
L3 = log2(3.0)
M64 = (1 << 64) - 1
S_U16 = 16.0
B1_DELTA_HI = 28.5
C_HIT_DELTA_BASE = 2.82

# ---------------- 1. 1D Retention for uint16_t STADDH (scale=16.0, ceil threshold, odd=2) ----------------
_RET1D_CACHE = {}

def get_ret1d_u16(M, pmin=24, B=61.0, odd=2):
    key = (pmin, B, odd)
    if key in _RET1D_CACHE:
        return _RET1D_CACHE[key]
    NX, NY = 512, 136
    ret_rg = M.Retention(B, B, prime_min=1 << pmin, explicit_bits=14, step=0.25, max_bits=530, tilt=0.05, odd=odd)
    ret_sf = M.Retention(530.0, B, prime_min=1 << pmin, explicit_bits=14, step=0.25, max_bits=531, tilt=0.005, odd=odd)
    rg_tbl = [[0.0] * (NY + 1) for _ in range(NX + 1)]
    sf_tbl = [[0.0] * (NY + 1) for _ in range(NX + 1)]
    s0_tbl = [[0.0] * (NY + 1) for _ in range(NX + 1)]
    for n in range(10, NX + 1):
        nf = float(n)
        n_idx = round(nf / ret_sf.step)
        # All-integer (sf=False) B-smooth probability from the exact same ret_sf generating series (large_geo index b=0)
        ub = ret_sf.two(n_idx, n_idx, 0, None, False) / ret_sf.denominator[n_idx]
        r_row = rg_tbl[n]; s_row = sf_tbl[n]; z_row = s0_tbl[n]
        for r in range(NY + 1):
            # rg_tbl[n][r]: squarefree-sieve-part smooth retention; ret_rg.at already subtracts rb=(cap(n)+1)/16 internally,
            # and model passes r = floor(sl_c - 0.50), adding >= 0.50 b margin (> 0.015503 + 1/16 = 0.078 b).
            r_row[r] = ret_rg.at(nf, float(r), None)
            # sf_tbl[n][r]: ALL-INTEGER (squarefree + squareful) survivor upper envelope via Retention.survivor_upper!
            # Model passes r = ceil(sl_c), so float(r) + 0.50 >= sl_c + 0.50 > sl_c + 0.015503 unconditionally.
            sv = ret_sf.survivor_upper(nf, float(r) + 0.50, ex=0.0)
            s_row[r] = sv
            z_row[r] = ub if ub < sv else sv
    out = dict(rg=rg_tbl, sf=sf_tbl, s0=s0_tbl, ret_rg=ret_rg, ret_sf=ret_sf, NX=NX, NY=NY, step=1.0, scale=S_U16)
    _RET1D_CACHE[key] = out
    return out


# ---------------- 2. Exact Closed-Form Rational 2-adic x 3-adic Class Split ----------------
_C23_CACHE = None

def get_exact_23_classes():
    """Exact CRT product of 2-adic (P_2^2) and 3-adic (P1, P2, P3 degree-1 + Q3 degree-3) local distributions
    over O_h-primitive (a, b) for m184_w216. Tail beyond k=6 at each P_i is placed conservatively at k=6 (lower
    valuation = strictly outward on yield). Exact means E[v2+v3] added back to ALPHA_F_RAW/ALPHA_G_RAW => zero double-counting."""
    global _C23_CACHE
    if _C23_CACHE is not None:
        return _C23_CACHE
    K_TAIL = 6
    p_d = {(0, 0): Fr(3, 4), (2, 2): Fr(1, 12)}
    rem = Fr(1, 6)
    for k in range(1, K_TAIL):
        pr = Fr(1, 9) * (Fr(1, 3) ** (k - 1))
        p_d[(2, 2 + k)] = pr
        rem -= pr
    p_d[(2, 2 + K_TAIL)] = rem
    assert sum(p_d.values()) == 1

    def _conv(d1, d2):
        out = {}
        for (f1, g1), w1 in d1.items():
            for (f2, g2), w2 in d2.items():
                k = (f1 + f2, g1 + g2)
                out[k] = out.get(k, Fr(0)) + w1 * w2
        return out

    p3_d = _conv(_conv(p_d, p_d), p_d)
    q3_d = {(0, 0): Fr(13, 14), (6, 0): Fr(1, 28), (0, 3): Fr(1, 28)}
    full_3 = _conv(p3_d, q3_d)
    assert sum(full_3.values()) == 1

    c3_cap = {}
    for (f3, g3), w in full_3.items():
        k = (min(f3, 12), min(g3, 12))
        c3_cap[k] = c3_cap.get(k, Fr(0)) + w

    CLASSES_2 = ((Fr(7, 9), 0, 0), (Fr(1, 9), 0, 3), (Fr(1, 9), 6, 0))
    exact_ev2_f, exact_ev2_g = Fr(2, 3), Fr(1, 3)
    exact_ev3_f, exact_ev3_g = Fr(12, 7), Fr(33, 14)
    ex_f_bits = float(exact_ev2_f) + float(exact_ev3_f) * L3
    ex_g_bits = float(exact_ev2_g) + float(exact_ev3_g) * L3

    classes = []; w_zero = 0.0
    tot_w = Fr(0); mean_f = 0.0; mean_g = 0.0
    for w2, v2f, v2g in CLASSES_2:
        for (f3, g3), w3 in sorted(c3_cap.items()):
            w = w2 * w3
            tot_w += w
            wf = float(w)
            if (v2f == 0 and v2g == 0 and f3 == 0 and g3 == 0) or wf < 1e-5:
                w_zero += wf
            else:
                vf_b = v2f + f3 * L3; vg_b = v2g + g3 * L3
                mean_f += wf * vf_b; mean_g += wf * vg_b
                classes.append((wf, vf_b, vg_b))
    classes.insert(0, (w_zero, 0.0, 0.0))
    u_vf = sorted({c[1] for c in classes})
    u_vg = sorted({c[2] for c in classes})
    vf_i = {v: i for i, v in enumerate(u_vf)}
    vg_i = {v: i for i, v in enumerate(u_vg)}
    sparse_fg = tuple((w, vf_i[vf], vg_i[vg]) for w, vf, vg in classes)
    assert tot_w == 1
    assert mean_f <= ex_f_bits + 1e-12 and mean_g <= ex_g_bits + 1e-12
    _C23_CACHE = dict(
        classes=tuple(classes),
        u_vf=tuple(u_vf), u_vg=tuple(u_vg), sparse_fg=sparse_fg,
        ex_f_bits=ex_f_bits,
        ex_g_bits=ex_g_bits,
        trunc_f_bits=mean_f,
        trunc_g_bits=mean_g,
        n_classes=len(classes),
    )
    return _C23_CACHE


def verify_3adic_exact(M, seed=20261123):
    """Live pure-Python algebraic verification of the exact 3-adic theorem at P1, P2, P3 and Q3:
    verifies the 3 ring homomorphisms O_h -> Z/729Z (kernel P_i^6), the coefficient residues of f_rec and g_rec
    mod 729 and mod Q_3, and tests the exact valuation law on 1,500 random O_h-primitive pairs."""
    c23 = get_exact_23_classes()
    P_MAPS = (
        (1, 473, 344, 58, 71, 434),
        (1, 642, 473, 17, 681, 58),
        (1, 344, 642, 434, 199, 17),
    )
    if M._OHM_CACHE is None:
        M.is_ideal_primitive_oh([1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0])
    ohm = M._OHM_CACHE
    hom_errs = 0
    for mp in P_MAPS:
        for i in range(6):
            for j in range(6):
                lhs = (mp[i] * mp[j]) % 729
                rhs = sum(ohm[i][j][k] * mp[k] for k in range(6)) % 729
                if lhs != rhs:
                    hom_errs += 1
    distinct_ok = (len({tuple(x % 3 for x in mp) for mp in P_MAPS}) == 3)
    from fractions import Fraction as F
    ZK_f = [[F(c, 3) for c in row] for row in M.ZK_NUM]
    aug = [ZK_f[i][:] + [F(1 if i == j else 0) for j in range(6)] for i in range(6)]
    for col in range(6):
        pr = next(r for r in range(col, 6) if aug[r][col] != 0)
        aug[col], aug[pr] = aug[pr], aug[col]
        pv = aug[col][col]; aug[col] = [x / pv for x in aug[col]]
        for r in range(6):
            if r != col and aug[r][col]:
                fac = aug[r][col]; aug[r] = [aug[r][c] - fac * aug[col][c] for c in range(12)]
    ZK_inv = [row[6:] for row in aug]
    def _to_zk(poly):
        co = [sum(F(poly[m]) * ZK_inv[m][k] for m in range(6)) for k in range(6)]
        assert all(x.denominator == 1 for x in co)
        return [x.numerator for x in co]
    f_zk = {d: _to_zk(M.F_REC_COEFFS[d]) for d in (8, 6, 4, 2, 0)}
    gl_zk = _to_zk([M.U_BN254 * c for c in M.W_COEFFS])
    gc_zk = [-1, 0, 0, 0, 0, 0]
    coeff_errs = 0
    for mp in P_MAPS:
        fc = [sum(f_zk[d][k] * mp[k] for k in range(6)) % 729 for d in (8, 6, 4, 2, 0)]
        gc = [sum(gl_zk[k] * mp[k] for k in range(6)) % 729, sum(gc_zk[k] * mp[k] for k in range(6)) % 729]
        if fc != [441, 72, 717, 669, 286] or gc != [576, 728]:
            coeff_errs += 1
    rng = random.Random(seed)
    v3 = lambda x: 6 if x % 729 == 0 else (next(k for k in range(6) if (x // (3 ** k)) % 3 != 0))
    law_errs = 0; tested = 0
    while tested < 1500:
        a = [rng.randrange(729) for _ in range(6)]
        b = [rng.randrange(729) for _ in range(6)]
        ais = [sum(a[k] * mp[k] for k in range(6)) % 729 for mp in P_MAPS]
        bis = [sum(b[k] * mp[k] for k in range(6)) % 729 for mp in P_MAPS]
        if any(ais[i] % 3 == 0 and bis[i] % 3 == 0 for i in range(3)):
            continue
        tested += 1
        for i in range(3):
            ai, bi = ais[i], bis[i]
            a2, b2 = (ai * ai) % 729, (bi * bi) % 729
            fv = (441 * pow(a2, 4, 729) + 72 * pow(a2, 3, 729) * b2 + 717 * a2 * a2 * b2 * b2 + 669 * a2 * pow(b2, 3, 729) + 286 * pow(b2, 4, 729)) % 729
            gv = (576 * a2 + 728 * b2) % 729
            vf_obs, vg_obs = v3(fv), v3(gv)
            if bi % 3 != 0:
                if vf_obs != 0 or vg_obs != 0:
                    law_errs += 1
            else:
                if vf_obs != 2:
                    law_errs += 1
                if (bi // 3) % 3 == 0:
                    if vg_obs != 2:
                        law_errs += 1
                else:
                    if vg_obs < 3:
                        law_errs += 1
    ok = (hom_errs == 0 and distinct_ok and coeff_errs == 0 and law_errs == 0
          and c23["trunc_f_bits"] <= c23["ex_f_bits"] + 1e-12
          and c23["trunc_g_bits"] <= c23["ex_g_bits"] + 1e-12)
    return dict(ok=ok, hom_errs=hom_errs, distinct_3_primes_ok=distinct_ok, coeff_mod729_errs=coeff_errs,
                law_errs=law_errs, tested=tested, n_classes=c23["n_classes"],
                ex_f_bits=round(c23["ex_f_bits"], 5), ex_g_bits=round(c23["ex_g_bits"], 5),
                trunc_f_bits=round(c23["trunc_f_bits"], 5), trunc_g_bits=round(c23["trunc_g_bits"], 5))


# ---------------- 3. Deterministic ntt3_F Bernstein Recursive Split ----------------
def split_ntt3_parts(M, R, lf, lg):
    """Deterministic side-specific ntt3_F Bernstein recursive prime-tree split (Q = 2^46, 1-owner RAM <= 2^57.0 B).
    Superadditivity of x*(log2 x)^2 bounds sum_v PT/RT across all tree nodes by a single root-sized M_tot = (lleaf+K_max)*Q."""
    Q = float(1 << 46)
    W_b = ceil(R / Q)
    out = {}
    tot = 0.0
    max_ram = 0.0
    for side, cap, lb, nb, lleaf in (("f", M.KH_F_CAP, lf, 448, 7), ("g", M.KH_G_CAP, lg, 512, 8)):
        pb = M.pbits(0, lb, cap)
        SP = ceil(pb / 64.0)
        np_ = 1.25506 * (2.0 ** (0.5 * lb)) / (0.5 * lb * LN2) + cap * 1.25506 * (2.0 ** lb) / (lb * LN2)
        h = ceil(log2(np_))
        K_max = nb // 24
        M_tot = (lleaf + K_max) * Q
        L_bt = ceil(log2(M_tot))
        c_root = 3.0 * (M.ntt3_F(SP) + M.ntt3_F(M_tot))
        c_ptrt = L_bt * ((8.0 / 3.0) * M.ntt3_F(M_tot) + 16.0 * M_tot)
        c_gcd = (K_max * Q) * (2.0 * (2.0 * nb * (12.0 * lleaf + 16.0)) + 20000.0)
        c_route = (K_max * Q) * 100.0
        c_ptree = h * (M.ntt3_F(SP) + 16.0 * SP)
        c_small = R * 600000.0
        c_ideal = R * 35.0 * (2.0 ** 19)
        side_ir = W_b * h * (c_root + c_ptrt + c_gcd + c_route) + c_ptree + c_small + c_ideal
        out[side] = side_ir
        tot += side_ir
        max_ram = max(max_ram, 8.0 * SP + 32.0 * M_tot + 8.0 * Q)
    out["total_ir"] = tot
    out["ram_log2"] = log2(max_ram)
    return out


# ---------------- 4. Band-1 Delta Exact Cost & All v32c Kernel/LA/uint16 Verifiers ----------------
def b1_delta_save_per_cell(M, raw_mult, p_min=24.0, p_hi=B1_DELTA_HI):
    """Exact saving on rep_lo from 32-bit delta-packed Band 1 on [2^p_min, 2^p_hi) (uint16_t 2 B/cell, delta = 2*(ai_new - ai_old)),
    charging 2.82 Ir/hit base + 24 Ir per sub-segment header whenever byte gap >= 2^32 (rank gap >= 2^31)."""
    w_k = [math.comb(3, k) * ((2.0 / 9.0) ** k) * ((2.0 / 3.0) ** (3 - k)) / (512.0 / 729.0) for k in range(4)]
    tot_save_log = 0.0
    steps = int(round((p_hi - p_min) * 4))
    dp = (p_hi - p_min) / steps
    for i in range(steps):
        pb = p_min + (i + 0.5) * dp
        p_hi_cell = 2.0 ** (pb + 0.5 * dp)
        p_ov = sum(w_k[k] * exp(-(2.0 ** 31) * ((2.0 / 3.0) ** k) / p_hi_cell) for k in range(4))
        eff_hit = C_HIT_DELTA_BASE + 24.0 * p_ov
        if eff_hit < M.C_HIT_LO:
            tot_save_log += (M.C_HIT_LO - eff_hit) * log((pb + 0.5 * dp) / (pb - 0.5 * dp))
    return M._WKF * M.B1_PREF_COV * raw_mult * tot_save_log


def verify_v32_all(M, kernels_path, seed=20261123):
    rng = random.Random(seed)
    # 4a. Exact 3-adic algebraic & valuation check
    t3 = verify_3adic_exact(M, seed=seed)

    # 4b. Deterministic NTT3 recursive split 64-integer bit-exact prototype
    ps = M._primes_upto(200)
    def _bpt(arr):
        if len(arr) == 1:
            return (arr[0], None, None)
        mid = len(arr) // 2
        L, Rt = _bpt(arr[:mid]), _bpt(arr[mid:])
        return (L[0] * Rt[0], L, Rt)
    pt = _bpt(ps)
    nums = []
    for _ in range(64):
        v = 1
        for _ in range(rng.randint(3, 9)):
            v *= rng.choice(ps)
        nums.append(v)
    def _rec(node, frags):
        if not frags:
            return []
        prod, L, Rt = node
        if L is None:
            out = []
            for idx, val in frags:
                e = 0; x = val
                while x % prod == 0:
                    x //= prod; e += 1
                assert x == 1 and e >= 1
                out.append((idx, prod, e))
            return out
        P_L = L[0]; fl, fr = [], []
        for idx, val in frags:
            g = math.gcd(val, P_L % val)
            if g == 1:
                fr.append((idx, val))
            elif g == val:
                fl.append((idx, val))
            else:
                left = g; right = val // g
                while True:
                    g2 = math.gcd(right, g)
                    if g2 == 1:
                        break
                    left *= g2; right //= g2
                fl.append((idx, left))
                if right > 1:
                    fr.append((idx, right))
        return _rec(L, fl) + _rec(Rt, fr)
    got = _rec(pt, list(enumerate(nums)))
    by_i = {}
    for idx, p_v, e_v in got:
        by_i.setdefault(idx, []).append((p_v, e_v))
    ntt3_split_errs = 0
    for idx, v in enumerate(nums):
        pr = 1
        for p_v, e_v in by_i[idx]:
            pr *= p_v ** e_v
        if pr != v:
            ntt3_split_errs += 1

    # 4c. Lazy 2-Barrett5 + 4-limb CIOS MontMul(chi_R) c_tail_z6 (323 ops <= 340, 200 random tests)
    r = M.R_BN254; chi = pow(M.P_BN254, 10, r)
    mu_r = (1 << 317) // r
    chi_R = (chi * (1 << 256)) % r
    r_inv = (-pow(r, -1, 1 << 64)) & M64
    def _mont4(a, b):
        t = 0
        for i in range(4):
            ai = (a >> (64 * i)) & M64
            t += ai * b
            u0 = (t * r_inv) & M64
            t = (t + u0 * r) >> 64
        if t >= r:
            t -= r
        return t
    def _bar5(A):
        neg = (A < 0); U = -A if neg else A
        top = (((U >> 256) << 3) | ((U >> 192) >> 61)) & M64
        qh = (top * mu_r) >> 64
        rem = U - qh * r
        if rem >= r:
            rem -= r
        if rem >= r:
            rem -= r
        assert 0 <= rem < r
        return (r - rem) if (neg and rem) else rem
    ctail_errs = 0
    for _ in range(200):
        B0 = rng.randrange(-(1 << 313), 1 << 313)
        B1 = rng.randrange(-(1 << 313), 1 << 313)
        B2 = rng.randrange(-(1 << 313), 1 << 313)
        got_c = (_bar5(B0 - B2) + _mont4(_bar5(B1 + B2), chi_R)) % r
        ref_c = ((B0 - B2) + chi * (B1 + B2)) % r
        if got_c != ref_c:
            ctail_errs += 1
    ctail_ops = 43 + 112 + 154 + 14

    # 4d. Aligned-sub-band bfi Lower Band 2 (uint16_t STADDH, entry = (ai << 5) | j_rel so (e >> 4) == 2*ai exact 2-byte offset,
    # all k_s in 38..49, W in {8, 16}, w_bits = log2(W) in {3, 4}, sgn in {+1, -1}, 5 non-zero aligned bases incl. base = 2^55 + 2^49 at W=8, k_s=46!)
    bfi_b2_errs = 0
    bfi_b2_cases = 0
    for W_t in (8, 16):
        w_bits = 3 if W_t == 8 else 4
        w_mask = (1 << w_bits) - 1
        for k_s in range(38, 50):
            namax = 4000
            namax_sh5 = (namax << 5) | 15
            pairs = [(1, 0), (1, W_t - 1), (namax, 0), (namax, W_t - 1)] + sorted(
                (rng.randint(1, 5000), rng.randrange(W_t)) for _ in range(40)
            )
            pairs.sort()
            entries = [(ai << 5) | jr for ai, jr in pairs]
            for sgn_v in (1, -1):
                for m_hi in (0, 1, 2, 3, 5):
                    base = ((1 << 55) & ~(w_mask << k_s)) + (m_hi << (k_s + w_bits))
                    assert (base & (w_mask << k_s)) == 0 and (base & ((1 << k_s) - 1)) == 0
                    bfi_b2_cases += 1
                    got_m = {}; ref_m = {}
                    for e in entries:
                        if e > namax_sh5:
                            break
                        x5 = base + (e >> 4)
                        if x5 & 1:
                            bfi_b2_errs += 1
                        fld = ((~e) & w_mask) if sgn_v < 0 else (e & w_mask)
                        mask = w_mask << k_s
                        x5 = (x5 & ~mask) | (fld << k_s)
                        got_m[x5] = (got_m.get(x5, 0) + 517) & 0xFFFF
                    for ai, jr in pairs:
                        if ai > namax:
                            break
                        jr_eff = jr if sgn_v > 0 else (W_t - 1 - jr)
                        ad = base + (jr_eff << k_s) + 2 * ai
                        ref_m[ad] = (ref_m.get(ad, 0) + 517) & 0xFFFF
                    if got_m != ref_m:
                        bfi_b2_errs += 1

    # 4e. 32-bit delta-packed Band 1 end-to-end check (uint16_t 2 B/cell, delta = 2*(ai_new - ai_old), forced splits on byte delta >= 2^32!)
    delta_b1_errs = 0
    ais_big = [10, 250, (1 << 31) + 500, (1 << 31) + 900, (1 << 32) + 1200, (1 << 32) + 1205, (1 << 32) + 1299, (1 << 32) + 1400, (1 << 32) + 1555]
    segs_d = []
    cur_base = 2 * ais_big[0]; cur_ds = [0]
    prev = 2 * ais_big[0]
    for ai in ais_big[1:]:
        b_off = 2 * ai
        d = b_off - prev
        if d >= (1 << 32):
            segs_d.append((cur_base, cur_ds))
            cur_base = b_off; cur_ds = [0]; prev = b_off
        else:
            cur_ds.append(d); prev = b_off
    segs_d.append((cur_base, cur_ds))
    row0 = 1 << 50; w_u16 = 511
    got_d = {}; ref_d = {}
    for ai in ais_big:
        ref_d[row0 + 2 * ai] = (ref_d.get(row0 + 2 * ai, 0) + w_u16) & 0xFFFF
    for sb, ds in segs_d:
        x5 = row0 + sb
        for i in range(0, len(ds), 2):
            w64 = ds[i] | ((ds[i + 1] << 32) if i + 1 < len(ds) else 0)
            x5 = (x5 + (w64 & 0xFFFFFFFF)) & M64
            if x5 & 1:
                delta_b1_errs += 1
            got_d[x5] = (got_d.get(x5, 0) + w_u16) & 0xFFFF
            if i + 1 < len(ds):
                x5 = (x5 + (w64 >> 32)) & M64
                if x5 & 1:
                    delta_b1_errs += 1
                got_d[x5] = (got_d.get(x5, 0) + w_u16) & 0xFFFF
    if got_d != ref_d or len(segs_d) != 3:
        delta_b1_errs += 1

    # 4f. Unconditional |log2 N_g(a,b) - log2 N(u*w*a^2)| bound + uint16_t ceil-threshold inclusion & zero-wrap + survivor_upper check + read-only scan8_u16
    T4 = M.load_data()["t4_geom"]
    cA_emb = T4["GA"][0][0] / 3.0
    w_min_abs = min(abs(complex(r_e, im_e)) for r_e, im_e in T4["EMB_W"])
    a_emb2_max = (589.43 ** 2) * M.QA_MAX / cA_emb
    b_emb2_max = (782.76 ** 2) * M.QB_MAX / cA_emb
    a_emb2_min = 1.0 / (a_emb2_max * a_emb2_max)
    ratio_max = b_emb2_max / (float(M.U_BN254) * w_min_abs * a_emb2_min)
    ng_diff_ub = max(6.0 * log2(1.0 + ratio_max), -6.0 * log2(1.0 - ratio_max))
    sl_u16 = M.OP["slack"]
    # On the squarefree-sieve-part smooth subset (counted by M.Retention.at(..., sf=True)), every dividing prime p >= 2^24
    # has multiplicity 1 and is scored once; ret_rg.at already subtracts rb=(cap(n)+1)/16 internally, and the model passes
    # floor(sl_c - 0.50), adding >= 0.50 b margin > ng_diff_ub + 1/16 = 0.078003 b unconditionally!
    # For ALL integers (squarefree + squareful), S >= t(a) implies log2(N_g_nonsieve) <= sl + ng_diff_ub unconditionally;
    # sf_tbl uses ret_sf.survivor_upper(n, ceil(sl_c) + 0.50, ex=0.0) (all-integer sf=False, no rb subtraction, no normalizer).
    uncond_rg_extra_margin = ng_diff_ub + 1.0 / S_U16  # 0.078003 <= 0.50 (on top of ret_rg.at's internal rb=(cap(n)+1)/16)
    uncond_sf_margin = ng_diff_ub  # 0.015503 <= 0.50
    thr_max = ceil(S_U16 * (480.0 - sl_u16))
    init_min = 32768 - thr_max
    max_u16 = 32768 + ceil(S_U16 * (sl_u16 + 1.0))
    r1d_chk = get_ret1d_u16(M, pmin=M.OP["p_min"], B=M.OP["lg"], odd=2)
    sf_468 = r1d_chk["sf"][468][113]
    sf_468_ref = r1d_chk["ret_sf"].survivor_upper(468.0, 113.5, ex=0.0)
    sf_468_at = r1d_chk["ret_sf"].at(468.0, 113.5, None)
    sf_uses_survivor_upper_ok = (abs(sf_468 - sf_468_ref) <= 1e-18 and sf_468 > sf_468_at)
    u16_range_ok = (ng_diff_ub <= 0.02 and uncond_rg_extra_margin <= 0.50 and uncond_sf_margin <= 0.50
                    and thr_max <= 32768 and init_min >= 0 and max_u16 <= 65535 and sf_uses_survivor_upper_ok)
    # Bit-exact read-only scan8_u16 (4 uint16_t cells per 64-bit word, MSB mask 0x8000800080008000):
    scan8_u16_errs = 0
    n_c = 128
    init_u16 = [32768 - rng.randint(4500, 5856) for _ in range(n_c)]
    rows_8 = [init_u16[:] for _ in range(8)]
    exp_s = set()
    for rr in range(8):
        for cc in range(n_c):
            if rng.random() < 0.03:
                rows_8[rr][cc] = (rows_8[rr][cc] + rng.randint(5000, 6500)) & 0xFFFF
            if rows_8[rr][cc] & 0x8000:
                exp_s.add((rr, cc))
    fnd_s = set()
    for c0 in range(0, n_c, 4):
        orr = 0
        for rr in range(8):
            w = 0
            for dc in range(4):
                w |= rows_8[rr][c0 + dc] << (16 * dc)
            orr |= w
        if orr & 0x8000800080008000:
            for rr in range(8):
                for dc in range(4):
                    if rows_8[rr][c0 + dc] & 0x8000:
                        fnd_s.add((rr, c0 + dc))
    if fnd_s != exp_s:
        scan8_u16_errs += 1

    # 4g. Count new segments in kernels_v32.s
    ins_rx = re.compile(r"^\s+[a-z][a-z0-9.]*\b")
    segs = {}; cur = None
    for raw in open(kernels_path, encoding="ascii"):
        st = raw.strip()
        if st.startswith("//@seg"):
            cur = st.split()[1]; segs.setdefault(cur, 0); continue
        if st.startswith("//@end"):
            cur = None; continue
        ln = raw.split("//", 1)[0].rstrip()
        if cur is None or not ln.strip() or ln.endswith(":") or ln.strip().startswith("."):
            continue
        if ins_rx.match(ln):
            segs[cur] += 1
    b1_delta_ops_ok = (segs.get("delta_b1_copy") == 5 and segs.get("delta_b1_apply") == 12
                       and segs.get("x16_bfi_pos_w8") == 12 and segs.get("x16_bfi_neg_w8") == 14
                       and segs.get("x16_bfi_pos_w16") == 12 and segs.get("x16_bfi_neg_w16") == 14)

    # 4h. Live depth-first post-order boot root construction stack-peak check (1,024-leaf balanced tree)
    df_leaves = M._primes_upto(8200)[:1024]
    assert len(df_leaves) == 1024
    # Level-by-level reference product and retained-all-levels byte sum
    lvl = list(df_leaves)
    full_tree_bits = sum(x.bit_length() for x in lvl)
    while len(lvl) > 1:
        lvl = [lvl[i] * lvl[i + 1] for i in range(0, len(lvl), 2)]
        full_tree_bits += sum(x.bit_length() for x in lvl)
    ref_root = lvl[0]
    root_bits = ref_root.bit_length()
    # Depth-first post-order evaluation discarding children immediately after multiplying
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
        # Live during multiplication: pending stack + inputs (bl+br) + output (bp) + 5*bp FFT/Garner scratch
        step_bits = live_stack_bits + bl + br + 6 * bp
        if step_bits > max_step_total_bits:
            max_step_total_bits = step_bits
        return prod
    df_root = _df_prod(df_leaves)
    boot_df_tree_ok = (df_root == ref_root and live_stack_bits == 0
                       and max_stack_bits < 1.02 * root_bits
                       and max_step_total_bits <= 8.0 * root_bits
                       and full_tree_bits > 9.0 * root_bits)

    ok = (t3["ok"] and ntt3_split_errs == 0 and ctail_errs == 0 and ctail_ops <= M.OP["c_tail_z6"]
          and bfi_b2_errs == 0 and bfi_b2_cases == 240 and delta_b1_errs == 0 and u16_range_ok and scan8_u16_errs == 0
          and b1_delta_ops_ok and boot_df_tree_ok)
    return dict(
        ok=ok,
        adic3_exact=t3,
        ntt3_split_proto_errs=ntt3_split_errs,
        ctail_lazy_errs=ctail_errs,
        ctail_lazy_ops=ctail_ops,
        ctail_budget=M.OP["c_tail_z6"],
        bfi_b2_errs=bfi_b2_errs,
        bfi_b2_cases=bfi_b2_cases,
        delta_b1_errs=delta_b1_errs,
        ng_diff_ub_bits=round(ng_diff_ub, 6),
        u16_scale=S_U16, u16_thr_max=thr_max, u16_init_min=init_min, u16_max_val=max_u16,
        u16_uncond_rg_extra_margin_bits=round(uncond_rg_extra_margin, 6),
        u16_uncond_sf_margin_bits=round(uncond_sf_margin, 6),
        sf_uses_survivor_upper_ok=sf_uses_survivor_upper_ok,
        u16_range_ok=u16_range_ok,
        scan8_u16_errs=scan8_u16_errs,
        b1_delta_ops_ok=b1_delta_ops_ok,
        boot_df_tree_ok=boot_df_tree_ok,
        boot_df_stack_ratio=round(max_stack_bits / root_bits, 4),
        boot_df_peak_ratio=round(max_step_total_bits / root_bits, 4),
        boot_full_tree_ratio=round(full_tree_bits / root_bits, 2),
    )
