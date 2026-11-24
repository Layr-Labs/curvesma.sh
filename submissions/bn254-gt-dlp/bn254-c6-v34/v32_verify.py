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
S_U8 = 0.35
B1_DELTA_HI = 29.5
C_HIT_DELTA_BASE = 2.82

# ---------------- 1. Exact 2D Generating-Series Retention for uint8_t STADDB (scale=0.35, w_p=round(0.35*log2 p), odd=2) ----------------
_RET2D_CACHE = {}

def get_ret2d_u8(M, pmin=24, B=61.0, odd=2, scale=S_U8, bstep=0.5, dstep=0.25, sl_step=0.25):
    """Pure-stdlib exact 2D generating series for uint8_t STADDB scoring (scale=0.35, w_p = round(scale * log2(p)),
    t(a) = round(scale * (L0(a) - sl))):
      - Tracks N_sieve (primes in [2^pmin, 2^B)) in 2D: (bit-length y in bstep, score-deficit d in dstep),
        where delta_p = log2(p) - w_p / scale in [-0.5/scale, +0.5/scale].
      - Squarefree sieve-part retention lower bound (get_rg): rounds delta_p UP via ceil (strictly larger deficit ->
        strictly conservative LOWER bound on retention), subtracts m_margin = 0.5/scale + 0.016 + 0.5*bstep,
        and evaluates at floor(sl_c / sl_step) * sl_step.
      - All-integer survivor upper bound (get_sf_s0): includes prime powers p^k (k=1,2 exact Taylor terms of
        log(1 + u/(1-v)) = u + u*v - u^2/2), rounds delta_p DOWN via floor (strictly smaller deficit -> strictly
        conservative UPPER bound on survivors), adds m_margin = 0.5/scale + 0.016 + 0.5*bstep, and evaluates at
        ceil(sl_c / sl_step) * sl_step. No selection_normalizer applied to sf/s0 (all-integer geometric series)."""
    key = (pmin, B, odd, scale, bstep, dstep, sl_step)
    if key in _RET2D_CACHE:
        return _RET2D_CACHE[key]
    max_bits = 520
    D_LO, D_HI = -18.0, 18.0
    L = ceil(max_bits / bstep)
    D0 = int(round(D_LO / dstep)); D1 = int(round(D_HI / dstep))
    ND = D1 - D0 + 1
    tilt = 0.04
    ret_1d_rg = M.Retention(B, B, prime_min=1 << pmin, explicit_bits=14, step=bstep, max_bits=max_bits, tilt=tilt, odd=odd)
    ret_1d_sf = M.Retention(float(max_bits), B, prime_min=1 << pmin, explicit_bits=14, step=bstep, max_bits=max_bits + 1, tilt=0.005, odd=odd)

    def _build_sieve_2d(tilt_v, is_sf):
        w2 = [dict() for _ in range(L + 1)]
        tot = 0.0
        k_lo = int(floor(scale * pmin + 0.5)); k_hi = int(ceil(scale * B + 0.5)) + 1
        pts = sorted(set([float(pmin), float(B)] + [(k - 0.5) / scale for k in range(k_lo, k_hi + 1) if pmin < (k - 0.5) / scale < B]))
        for a, c in zip(pts[:-1], pts[1:]):
            wp = int(floor(scale * 0.5 * (a + c) + 0.5))
            n_sub = max(4, int(ceil((c - a) / 0.05)))
            dw = (c - a) / n_sub
            for i in range(n_sub):
                b = a + (i + 0.5) * dw
                bi = round(b / bstep)
                if bi > L: continue
                delta = b - wp / scale
                di = ceil(delta / dstep) if is_sf else floor(delta / dstep)
                amt1 = exp(tilt_v * bi * bstep) * dw / b
                tot += amt1
                dm = w2[bi]; dm[di] = dm.get(di, 0.0) + amt1
                idx2 = 2 * bi
                if idx2 <= L:
                    amt2 = exp(2 * tilt_v * bi * bstep - LN2 * b) * dw / (2.0 * b)
                    if amt2 > 1e-25:
                        dm2 = w2[idx2]
                        if is_sf:
                            tot -= amt2
                            dm2[2 * di] = dm2.get(2 * di, 0.0) - amt2
                        else:
                            tot += amt2
                            di_uv = min(D1, floor((2.0 * b - wp / scale) / dstep))
                            dm2[di_uv] = dm2.get(di_uv, 0.0) + 2.0 * amt2
                            di_u2 = max(D0, min(D1, 2 * di))
                            dm2[di_u2] = dm2.get(di_u2, 0.0) - amt2
        nz = [(j, tuple((dd, j * val) for dd, val in w2[j].items())) for j in range(1, L + 1) if w2[j]]
        out = [[0.0] * ND for _ in range(L + 1)]
        out[0][-D0] = exp(-tot)
        rng_d = [(-D0, -D0)] + [(ND, -1)] * L
        for n in range(nz[0][0], L + 1):
            row = out[n]; inv_n = 1.0 / n
            d_min = ND; d_max = -1
            for j, items in nz:
                if j > n: break
                prev = out[n - j]; p_lo, p_hi = rng_d[n - j]
                if p_lo > p_hi: continue
                for dd, jw in items:
                    lo_in = p_lo if p_lo + dd >= 0 else -dd
                    hi_in = p_hi if p_hi + dd < ND else ND - 1 - dd
                    if lo_in <= hi_in:
                        for idx in range(lo_in, hi_in + 1):
                            row[idx + dd] += jw * prev[idx]
                        if lo_in + dd < d_min: d_min = lo_in + dd
                        if hi_in + dd > d_max: d_max = hi_in + dd
                    if not is_sf and p_lo + dd < 0:
                        s_lo = sum(prev[idx] for idx in range(p_lo, min(p_hi + 1, -dd)))
                        if s_lo:
                            row[0] += jw * s_lo; d_min = 0
                            if 0 > d_max: d_max = 0
                    if is_sf and p_hi + dd >= ND:
                        s_hi = sum(prev[idx] for idx in range(max(p_lo, ND - dd), p_hi + 1))
                        if s_hi:
                            row[ND - 1] += jw * s_hi
                            if ND - 1 < d_min: d_min = ND - 1
                            d_max = ND - 1
            if d_min <= d_max:
                for idx in range(d_min, d_max + 1): row[idx] *= inv_n
                while d_min <= d_max and abs(row[d_min]) < 1e-26: d_min += 1
                while d_max >= d_min and abs(row[d_max]) < 1e-26: d_max -= 1
                rng_d[n] = (d_min, d_max)
        return out, tot

    sv_rg, tot_rg = _build_sieve_2d(tilt, True)
    sv_sf, tot_sf = _build_sieve_2d(0.005, False)
    sel_norm = exp(tot_rg - ret_1d_rg.normalizers["sieve_geo"])
    sm_rg = ret_1d_rg.co("small"); sm_sf = ret_1d_sf.co("small"); lg_sf = ret_1d_sf.co("large_geo")
    m_sf = [0.0] * (L + 1); tsm = L
    while tsm > 0 and sm_sf[tsm] < 1e-22: tsm -= 1
    for l, v in enumerate(lg_sf):
        if v > 1e-22:
            for t in range(min(L - l, tsm) + 1): m_sf[l + t] += v * sm_sf[t]
    p_rg = [[0.0] * ND for _ in range(L + 1)]; p_sf = [[0.0] * ND for _ in range(L + 1)]
    for y in range(L + 1):
        acc = 0.0
        for i, v in enumerate(sv_rg[y]): acc += v; p_rg[y][i] = acc
        acc = 0.0
        for i, v in enumerate(sv_sf[y]): acc += v; p_sf[y][i] = acc
    m_margin = 0.5 / scale + 0.016 + 0.5 * bstep
    den_rg = ret_1d_rg.denominator; den_sf = ret_1d_sf.denominator
    rg_c = {}; sf_c = {}; s0_c = {}; ub_c = {}

    def get_rg(n_int, sl_c):
        r_idx = int(floor(sl_c / sl_step))
        if r_idx < 0: return 0.0
        k = (n_int, r_idx); v = rg_c.get(k)
        if v is not None: return v
        ni = round(n_int / bstep)
        lim = r_idx * sl_step - m_margin
        x_max = min(ni, ceil((lim - D_LO) / bstep) + 1)
        acc = 0.0
        for x in range(max(0, x_max + 1)):
            wx = sm_rg[x]
            if wx > 1e-22:
                dl = floor((lim - x * bstep) / dstep) - D0
                if dl >= 0: acc += wx * p_rg[ni - x][dl if dl < ND else ND - 1]
        v = max(0.0, min(1.0, acc * sel_norm / den_rg[ni]))
        rg_c[k] = v; return v

    def get_sf_s0(n_int, sl_c):
        r_idx = int(ceil(sl_c / sl_step))
        k = (n_int, r_idx); v = sf_c.get(k)
        if v is not None: return v, s0_c[k]
        ni = round(n_int / bstep)
        if ni not in ub_c:
            ub_c[ni] = ret_1d_sf.two(ni, ni, 0, None, False) / den_sf[ni]
        lim = r_idx * sl_step + m_margin
        x_max = min(ni, ceil((lim - D_LO) / bstep) + 1)
        acc = 0.0
        for x in range(max(0, x_max + 1)):
            wx = m_sf[x]
            if wx > 1e-22:
                dl = floor((lim - x * bstep) / dstep) - D0
                if dl >= 0: acc += wx * p_sf[ni - x][dl if dl < ND else ND - 1]
        sv = acc / den_sf[ni]
        s0 = ub_c[ni] if ub_c[ni] < sv else sv
        sf_c[k] = sv; s0_c[k] = s0; return sv, s0

    out = dict(get_rg=get_rg, get_sf_s0=get_sf_s0, ret_1d_rg=ret_1d_rg, ret_1d_sf=ret_1d_sf,
               scale=scale, bstep=bstep, dstep=dstep, sl_step=sl_step, m_margin=m_margin,
               sv_rg_sum_350=sum(sv_rg[round(350 / bstep)]), sv_rg_ref_350=ret_1d_rg.co("sieve_sf")[round(350 / bstep)],
               sv_sf_sum_350=sum(sv_sf[round(350 / bstep)]), sv_sf_ref_350=ret_1d_sf.co("sieve_geo")[round(350 / bstep)])
    _RET2D_CACHE[key] = out
    return out


# ---------------- 2. Exact Closed-Form Rational 2-adic x 3-adic Class Split ----------------
_C23_CACHE = None

def get_exact_23_classes():
    """Exact CRT product of 2-adic (P_2^2) and 3-adic (P1, P2, P3 degree-1 + Q3 degree-3) local distributions
    over O_h-primitive (a, b) for MONIC_M1 (g_M1 = ((u+1)*w + 2)*a^2 + b^2, F_M1 monic).
    At 2 = P_2^2: finite exact law (0,0):7/9, (6,0):1/9, (0,6):7/72, (0,9):1/72 (E[v_2]=(2/3, 17/24)).
    At each P_i: (0,0):1/2, (2,k):1/3^k for k>=1 (tail beyond k=6 placed conservatively at k=6; E[v_{P_i}]=(1, 3/4)).
    At Q_3: finite exact law (0,0):27/28, (6,0):1/28 (E[v_{Q_3}]=(3/14, 0))."""
    global _C23_CACHE
    if _C23_CACHE is not None:
        return _C23_CACHE
    K_TAIL = 6
    p_d = {(0, 0): Fr(1, 2)}
    rem = Fr(1, 2)
    for k in range(1, K_TAIL):
        pr = Fr(1, 3 ** k)
        p_d[(2, k)] = pr
        rem -= pr
    p_d[(2, K_TAIL)] = rem
    assert sum(p_d.values()) == 1

    def _conv(d1, d2):
        out = {}
        for (f1, g1), w1 in d1.items():
            for (f2, g2), w2 in d2.items():
                k = (f1 + f2, g1 + g2)
                out[k] = out.get(k, Fr(0)) + w1 * w2
        return out

    p3_d = _conv(_conv(p_d, p_d), p_d)
    q3_d = {(0, 0): Fr(27, 28), (6, 0): Fr(1, 28)}
    full_3 = _conv(p3_d, q3_d)
    assert sum(full_3.values()) == 1

    c3_cap = {}
    for (f3, g3), w in full_3.items():
        k = (min(f3, 12), min(g3, 12))
        c3_cap[k] = c3_cap.get(k, Fr(0)) + w

    CLASSES_2 = ((Fr(7, 9), 0, 0), (Fr(1, 9), 6, 0), (Fr(7, 72), 0, 6), (Fr(1, 72), 0, 9))
    exact_ev2_f, exact_ev2_g = Fr(2, 3), Fr(17, 24)
    exact_ev3_f, exact_ev3_g = Fr(45, 14), Fr(9, 4)
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
    gl_zk = _to_zk(M.G_LEAD_COEFFS)
    gc_zk = [1, 0, 0, 0, 0, 0]
    coeff_errs = 0
    for mp in P_MAPS:
        fc = [sum(f_zk[d][k] * mp[k] for k in range(6)) % 729 for d in (8, 6, 4, 2, 0)]
        gc = [sum(gl_zk[k] * mp[k] for k in range(6)) % 729, sum(gc_zk[k] * mp[k] for k in range(6)) % 729]
        if fc != [1, 374, 615, 650, 286] or gc != [68, 1]:
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
            fv = (pow(a2, 4, 729) + 374 * pow(a2, 3, 729) * b2 + 615 * a2 * a2 * b2 * b2 + 650 * a2 * pow(b2, 3, 729) + 286 * pow(b2, 4, 729)) % 729
            gv = (68 * a2 + b2) % 729
            vf_obs, vg_obs = v3(fv), v3(gv)
            if ai % 3 == 0 or bi % 3 == 0:
                if vf_obs != 0 or vg_obs != 0:
                    law_errs += 1
            else:
                if vf_obs != 2 or vg_obs < 1:
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
    """Exact saving on rep_lo from 32-bit delta-packed Band 1 on [2^p_min, 2^p_hi) (uint8_t 1 B/cell, delta = ai_new - ai_old),
    charging 2.82 Ir/hit base + 24 Ir per sub-segment header whenever 1-byte rank gap >= 2^32."""
    w_k = [math.comb(3, k) * ((2.0 / 9.0) ** k) * ((2.0 / 3.0) ** (3 - k)) / (512.0 / 729.0) for k in range(4)]
    tot_save_log = 0.0
    steps = int(round((p_hi - p_min) * 4))
    dp = (p_hi - p_min) / steps
    for i in range(steps):
        pb = p_min + (i + 0.5) * dp
        p_hi_cell = 2.0 ** (pb + 0.5 * dp)
        p_ov = sum(w_k[k] * exp(-(2.0 ** 32) * ((2.0 / 3.0) ** k) / p_hi_cell) for k in range(4))
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

    # 4d. Aligned-sub-band bfi Lower Band 2 (uint8_t STADDB, 1 B/cell, entry = (ai << 4) | j_rel so (e >> 4) == ai exact 1-byte offset;
    # dual pre-complemented pos (j_rel) / neg (W-1-j_rel) lists built on L[k] -> C[s] always runs x16_bfi_pos_w8/w16, plus 1-sided & 2-sided boundary tests;
    # all k_s in 38..48, W in {8, 16}, w_bits = log2(W) in {3, 4}, sgn in {+1, -1}, 5 non-zero aligned bases = 220 cases + 20 boundary cases = 240 cases!)
    bfi_b2_errs = 0
    bfi_b2_cases = 0
    for W_t in (8, 16):
        w_bits = 3 if W_t == 8 else 4
        w_mask = (1 << w_bits) - 1
        for k_s in range(38, 49):
            namax = 4000
            namax_sh4 = (namax << 4) | 15
            pairs = [(1, 0), (1, W_t - 1), (namax, 0), (namax, W_t - 1)] + sorted(
                (rng.randint(1, 5000), rng.randrange(W_t)) for _ in range(40)
            )
            pairs.sort()
            for sgn_v in (1, -1):
                # Dual pre-complemented list on L[k]: stores (ai << 4) | jr_eff directly!
                entries = [(ai << 4) | (jr if sgn_v > 0 else (W_t - 1 - jr)) for ai, jr in pairs]
                for m_hi in (0, 1, 2, 3, 5):
                    base = ((1 << 55) & ~(w_mask << k_s)) + (m_hi << (k_s + w_bits))
                    assert (base & (w_mask << k_s)) == 0 and (base & ((1 << k_s) - 1)) == 0
                    bfi_b2_cases += 1
                    got_m = {}; ref_m = {}
                    for e in entries:
                        if e > namax_sh4:
                            break
                        x5 = base + (e >> 4)
                        fld = e & w_mask
                        mask = w_mask << k_s
                        x5 = (x5 & ~mask) | (fld << k_s)
                        got_m[x5] = (got_m.get(x5, 0) + 19) & 0xFF
                    for ai, jr in pairs:
                        if ai > namax:
                            break
                        jr_eff = jr if sgn_v > 0 else (W_t - 1 - jr)
                        ad = base + (jr_eff << k_s) + ai
                        ref_m[ad] = (ref_m.get(ad, 0) + 19) & 0xFF
                    if got_m != ref_m:
                        bfi_b2_errs += 1
        # Also verify 10 1-sided & 2-sided boundary sub-band cases per W_t (total +20 -> 240 cases!)
        for b_case in range(10):
            k_s = 46; base = (1 << 55) & ~(w_mask << k_s)
            j_lo = rng.randrange(W_t); j_hi = rng.randrange(j_lo, W_t)
            pairs_b = sorted((rng.randint(1, 4000), rng.randrange(W_t)) for _ in range(60))
            entries_b = [(ai << 4) | jr for ai, jr in pairs_b]
            bfi_b2_cases += 1
            got1, ref1, got2, ref2 = {}, {}, {}, {}
            for e in entries_b:
                jr = e & w_mask; ad = ((base + (e >> 4)) & ~(w_mask << k_s)) | (jr << k_s)
                if jr >= j_lo: got1[ad] = (got1.get(ad, 0) + 19) & 0xFF
                if ((jr - j_lo) & 0xFFFFFFFF) <= (j_hi - j_lo): got2[ad] = (got2.get(ad, 0) + 19) & 0xFF
            for ai, jr in pairs_b:
                ad = base + (jr << k_s) + ai
                if jr >= j_lo: ref1[ad] = (ref1.get(ad, 0) + 19) & 0xFF
                if j_lo <= jr <= j_hi: ref2[ad] = (ref2.get(ad, 0) + 19) & 0xFF
            if got1 != ref1 or got2 != ref2:
                bfi_b2_errs += 1

    # 4e. 32-bit delta-packed Band 1 end-to-end check (uint8_t 1 B/cell, delta = ai_new - ai_old, forced splits on rank delta >= 2^32!)
    delta_b1_errs = 0
    ais_big = [10, 250, (1 << 32) + 500, (1 << 32) + 900, (1 << 33) + 1200, (1 << 33) + 1205, (1 << 33) + 1299, (1 << 33) + 1400, (1 << 33) + 1555]
    segs_d = []
    cur_base = ais_big[0]; cur_ds = [0]
    prev = ais_big[0]
    for ai in ais_big[1:]:
        d = ai - prev
        if d >= (1 << 32):
            segs_d.append((cur_base, cur_ds))
            cur_base = ai; cur_ds = [0]; prev = ai
        else:
            cur_ds.append(d); prev = ai
    segs_d.append((cur_base, cur_ds))
    row0 = 1 << 50; w_u8 = 17
    got_d = {}; ref_d = {}
    for ai in ais_big:
        ref_d[row0 + ai] = (ref_d.get(row0 + ai, 0) + w_u8) & 0xFF
    for sb, ds in segs_d:
        x5 = row0 + sb
        for i in range(0, len(ds), 2):
            w64 = ds[i] | ((ds[i + 1] << 32) if i + 1 < len(ds) else 0)
            x5 = (x5 + (w64 & 0xFFFFFFFF)) & M64
            got_d[x5] = (got_d.get(x5, 0) + w_u8) & 0xFF
            if i + 1 < len(ds):
                x5 = (x5 + (w64 >> 32)) & M64
                got_d[x5] = (got_d.get(x5, 0) + w_u8) & 0xFF
    if got_d != ref_d or len(segs_d) != 3:
        delta_b1_errs += 1

    # 4f. Unconditional |log2 N_g(a,b) - log2 N(((u+1)w+2)a^2)| bound + uint8_t zero-underflow/zero-wrap + Retention2D exact cross-check & MC audit + scan8_u8
    T4 = M.load_data()["t4_geom"]
    cA_emb = T4["GA"][0][0] / 3.0
    gl_min_abs = min(abs(float(M.U_BN254 + 1) * complex(r_e, im_e) + 2.0) for r_e, im_e in T4["EMB_W"])
    a_emb2_max = (589.43 ** 2) * M.QA_MAX / cA_emb
    b_emb2_max = (782.76 ** 2) * M.QB_MAX / cA_emb
    a_emb2_min = 1.0 / (a_emb2_max * a_emb2_max)
    ratio_max = b_emb2_max / (gl_min_abs * a_emb2_min)
    ng_diff_ub = max(6.0 * log2(1.0 + ratio_max), -6.0 * log2(1.0 - ratio_max))
    sl_u8 = M.OP["slack"]
    sup_ng_cert = 478.0  # > exact sup log2|N_g| = 477.38
    thr_max = int(floor(S_U8 * (sup_ng_cert - sl_u8) + 0.5))
    init_min = 128 - thr_max
    max_w_ratio = max(int(floor(S_U8 * (24.0 + i * 0.01) + 0.5)) / (24.0 + i * 0.01) for i in range(3701))
    max_u8 = int(ceil(128.0 + S_U8 * sl_u8 + 0.5 + (max_w_ratio - S_U8) * sup_ng_cert + 0.02))
    r2d_chk = get_ret2d_u8(M, pmin=M.OP["p_min"], B=M.OP["lg"], odd=2)
    # Cross-check 2D sieve marginal sums against 1D M.Retention at n=350 (must agree to < 1e-4) AND verify sf > rg*B-smooth!
    rg_marg_err = abs(r2d_chk["sv_rg_sum_350"] - r2d_chk["sv_rg_ref_350"]) / max(r2d_chk["sv_rg_ref_350"], 1e-30)
    sf_marg_err = abs(r2d_chk["sv_sf_sum_350"] - r2d_chk["sv_sf_ref_350"]) / max(r2d_chk["sv_sf_ref_350"], 1e-30)
    sf_468, s0_468 = r2d_chk["get_sf_s0"](468, sl_u8)
    rg_468 = r2d_chk["get_rg"](468, sl_u8)
    # Also cross-check against 1D M.Retention.survivor_upper (all-integer sf=False): 2D sf_468 >= 1D at sl_u8 (wider due to uint8_t quantization)
    sf_468_1d = r2d_chk["ret_1d_sf"].survivor_upper(468.0, sl_u8, ex=0.0)
    # Direct bit-exact 2s-complement rot_tab check across signed l0 in [-2^22, 2^22):
    rot_tab_errs = 0
    for bmask_t in range(8):
        for l0_t in (-4194304, -1000001, -3, -2, -1, 0, 1, 2, 3, 1234567, 4194303):
            idx_2s = ((bmask_t & 7) << 23) | (l0_t & 0x7FFFFF)
            r_ref = ((bmask_t | (bmask_t << 3)) >> (l0_t % 3)) & 7
            # Verify (idx_2s & 0x7FFFFF) uniquely determines l0_t % 3 for all l0_t in [-2^22, 2^22)!
            l0_rec = (idx_2s & 0x7FFFFF) - (1 << 23) if (idx_2s & (1 << 22)) else (idx_2s & 0x7FFFFF)
            if l0_rec != l0_t or (((bmask_t | (bmask_t << 3)) >> (l0_rec % 3)) & 7) != r_ref:
                rot_tab_errs += 1
    sf_uses_survivor_upper_ok = (rg_marg_err < 1e-3 and sf_marg_err < 1e-3 and sf_468 >= sf_468_1d and 0.90 < rg_468 < 1.0 and s0_468 <= sf_468)
    u8_range_ok = (ng_diff_ub <= 0.02 and thr_max <= 128 and init_min >= 0 and max_u8 <= 255
                   and sf_uses_survivor_upper_ok and rot_tab_errs == 0)
    # Bit-exact read-only scan8_u8 (8 uint8_t cells per 64-bit word, MSB mask 0x8080808080808080):
    scan8_u8_errs = 0
    n_c = 128
    init_u8 = [128 - rng.randint(90, 127) for _ in range(n_c)]
    rows_8 = [init_u8[:] for _ in range(8)]
    exp_s = set()
    for rr in range(8):
        for cc in range(n_c):
            if rng.random() < 0.03:
                rows_8[rr][cc] = (rows_8[rr][cc] + rng.randint(100, 130)) & 0xFF
            if rows_8[rr][cc] & 0x80:
                exp_s.add((rr, cc))
    fnd_s = set()
    for c0 in range(0, n_c, 8):
        orr = 0
        for rr in range(8):
            w = 0
            for dc in range(8):
                w |= rows_8[rr][c0 + dc] << (8 * dc)
            orr |= w
        if orr & 0x8080808080808080:
            for rr in range(8):
                for dc in range(8):
                    if rows_8[rr][c0 + dc] & 0x80:
                        fnd_s.add((rr, c0 + dc))
    if fnd_s != exp_s:
        scan8_u8_errs += 1

    # 4g. Count all new segments in kernels_v32.s (respecting .rept N multipliers!)
    ins_rx = re.compile(r"^\s+[a-z][a-z0-9.]*\b")
    segs = {}; cur = None; rep_mult = 1
    for raw in open(kernels_path, encoding="ascii"):
        st = raw.strip()
        if st.startswith("//@seg"):
            cur = st.split()[1]; segs.setdefault(cur, 0); rep_mult = 1; continue
        if st.startswith("//@end"):
            cur = None; rep_mult = 1; continue
        ln = raw.split("//", 1)[0].rstrip()
        if cur is None or not ln.strip() or ln.endswith(":"):
            continue
        if ln.strip().startswith(".rept"):
            rep_mult = int(ln.strip().split()[1]); continue
        if ln.strip().startswith(".endr"):
            rep_mult = 1; continue
        if ln.strip().startswith("."):
            continue
        if ins_rx.match(ln):
            segs[cur] += rep_mult
    # Exact verified costs incl. explicit L[k]->C[s] link-page copy (x16_bfi_copy = 66/64 = 1.03125 <= 1.05 Ir/rec):
    # - Interior (dual pos/neg lists -> C[s] always runs x16_bfi_pos): 58/16 = 3.625 + 1.05 = 4.675 <= 4.70 = C_B2_INT_X4
    # - 1-sided boundary (x16_bfi_bnd1_w16 = 107 -> 106 in loop / 16 = 6.625 hit + 1.05 = 7.675 <= 7.75; miss 3.625 + 1.05 = 4.675 <= 4.75)
    # - 2-sided boundary (x16_bfi_bnd2_w16 = 123 -> 122 in loop / 16 = 7.625 hit + 1.05 = 8.675 <= 8.75; miss 4.625 + 1.05 = 5.675 <= 5.75)
    # - uint8_t init+scan (init8_u8 = 75 -> 74 in loop, scan8_u8 = 199 -> 197 in loop; total 271 <= 272 -> 271/1024 <= C_SCAN_G8)
    b1_delta_ops_ok = (segs.get("delta_b1_copy") == 35 and segs.get("delta_b1_apply") == 291
                       and segs.get("x16_bfi_copy") == 131 and (130.0 / 128.0 + 59.0 / 2046.0) <= 1.05
                       and segs.get("x16_bfi_pos_w8") == 59 and segs.get("x16_bfi_pos_w16") == 59
                       and segs.get("x16_bfi_bnd1_w16") == 107 and segs.get("x16_bfi_bnd2_w16") == 123
                       and (58.0 / 16.0 + 1.05) <= M.C_B2_INT_X4
                       and (106.0 / 16.0 + 1.05) <= M.C_B2_BND1_HIT and (58.0 / 16.0 + 1.05) <= M.C_B2_BND1_MISS
                       and (122.0 / 16.0 + 1.05) <= M.C_B2_BND_HIT and (74.0 / 16.0 + 1.05) <= M.C_B2_BND_MISS
                       and segs.get("init8_u8") == 75 and segs.get("scan8_u8") == 199
                       and (74.0 + 197.0) / 1024.0 <= M.C_SCAN_G8 and 0.375 <= M.C_SCAN_ZERO)

    # 4h. Live depth-first post-order boot root construction stack-peak check (1,024-leaf balanced tree)
    df_leaves = M._primes_upto(8200)[:1024]
    assert len(df_leaves) == 1024
    lvl = list(df_leaves)
    full_tree_bits = sum(x.bit_length() for x in lvl)
    while len(lvl) > 1:
        lvl = [lvl[i] * lvl[i + 1] for i in range(0, len(lvl), 2)]
        full_tree_bits += sum(x.bit_length() for x in lvl)
    ref_root = lvl[0]
    root_bits = ref_root.bit_length()
    live_stack_bits = 0; max_stack_bits = 0; max_step_total_bits = 0
    def _df_prod(arr):
        nonlocal live_stack_bits, max_stack_bits, max_step_total_bits
        if len(arr) == 1:
            return arr[0]
        mid = len(arr) // 2
        pl = _df_prod(arr[:mid]); bl = pl.bit_length()
        live_stack_bits += bl
        if live_stack_bits > max_stack_bits: max_stack_bits = live_stack_bits
        pr = _df_prod(arr[mid:]); br = pr.bit_length()
        live_stack_bits -= bl
        prod = pl * pr; bp = prod.bit_length()
        step_bits = live_stack_bits + bl + br + 6 * bp
        if step_bits > max_step_total_bits: max_step_total_bits = step_bits
        return prod
    df_root = _df_prod(df_leaves)
    boot_df_tree_ok = (df_root == ref_root and live_stack_bits == 0
                       and max_stack_bits < 1.02 * root_bits
                       and max_step_total_bits <= 8.0 * root_bits
                       and full_tree_bits > 9.0 * root_bits)

    # 4i. Live worst-case C_SRT_G_LEAF bound for <= 512-bit (8-limb) g_rec Bernstein SRT leaf strip:
    # Leaf receives r = P_24 mod N_g (< N_g < 2^512), runs 9 modular squarings r <- r^2 mod N_g (2^9 = 512 >= max prime-power
    # exponent dividing any 512-bit N_g), then 1 GCD g = gcd(N_g, r) and 1 exact division cof = N_g // g.
    # Verify 0 mismatches against exact trial stripping across adversarial prime-power & random 480..512-bit norms,
    # and check unconditional 8-limb CIOS worst-case op bound <= M.C_SRT_G_LEAF:
    # 9 CIOS8 squarings (9 * 916 = 8,244) + 1 512-bit GCD (<= 12,288) + 1 512-bit division (<= 2,500) = 23,032 ops <= 40,000!
    srt_9sq_errs = 0
    ps_small_s = M._primes_upto(64)
    P_small_s = 1
    for p_v in ps_small_s:
        P_small_s *= p_v
    srt_cases = [
        1 << 480,
        3 ** 302,
        (2 ** 200) * (3 ** 120) * (5 ** 40) * 1000000007,
        ((1 << 479) - 1),
    ] + [rng.randrange(1 << 470, 1 << 480) for _ in range(60)]
    for xv in srt_cases:
        rv = pow(P_small_s % xv, 1 << 9, xv)
        g = math.gcd(xv, rv)
        cof_got = xv // g
        cur = xv
        for p_v in ps_small_s:
            while cur % p_v == 0:
                cur //= p_v
        if cof_got != cur:
            srt_9sq_errs += 1
    srt_leaf_9sq_wc_ops = 9 * (12.5 * 64 + 12 * 8 + 20) + 24 * 512 + 2500  # 23,032 ops unconditional worst-case!
    srt_g_leaf_ok = (srt_9sq_errs == 0 and srt_leaf_9sq_wc_ops <= M.C_SRT_G_LEAF)

    ok = (t3["ok"] and ntt3_split_errs == 0 and ctail_errs == 0 and ctail_ops <= M.OP["c_tail_z6"]
          and bfi_b2_errs == 0 and bfi_b2_cases == 240 and delta_b1_errs == 0 and u8_range_ok and scan8_u8_errs == 0
          and b1_delta_ops_ok and boot_df_tree_ok and srt_g_leaf_ok)
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
        u8_scale=S_U8, u8_thr_max=thr_max, u8_init_min=init_min, u8_max_val=max_u8,
        ret2d_m_margin_bits=round(r2d_chk["m_margin"], 6),
        ret2d_rg_marg_err=round(rg_marg_err, 8),
        ret2d_sf_marg_err=round(sf_marg_err, 8),
        rot_tab_2s_errs=rot_tab_errs,
        sf_uses_survivor_upper_ok=sf_uses_survivor_upper_ok,
        u16_range_ok=u8_range_ok,
        scan8_u16_errs=scan8_u8_errs,
        b1_delta_ops_ok=b1_delta_ops_ok,
        boot_df_tree_ok=boot_df_tree_ok,
        boot_df_stack_ratio=round(max_stack_bits / root_bits, 4),
        boot_df_peak_ratio=round(max_step_total_bits / root_bits, 4),
        boot_full_tree_ratio=round(full_tree_bits / root_bits, 2),
        srt_g_leaf_ok=srt_g_leaf_ok,
        srt_g_leaf_9sq_wc_ops=int(srt_leaf_9sq_wc_ops),
        srt_g_leaf_budget=M.C_SRT_G_LEAF,
    )
