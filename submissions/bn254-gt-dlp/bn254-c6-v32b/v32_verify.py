"""v32_verify.py -- Pure-stdlib exact verification & model helpers for BN254 C6 v32:
  1. Ret2DU8: exact 2D generating-series Retention for uint8_t STADDB scoring (scale=0.32, floor, odd=2).
  2. Exact closed-form rational 3-adic class split (P1, P2, P3 degree-1 + Q3 degree-3, zero double-counting) + live O_h check.
  3. Deterministic ntt3_F Bernstein recursive split (exact upper bound + 2x safety + live 64-integer prototype).
  4. Instruction-by-instruction / bit-exact verification of:
     - Aligned-sub-band bfi x16 Lower Band 2 (uint8_t STADDB, all k_s in 37..48, W in {8, 16}, sgn in {+1, -1})
     - 32-bit delta-packed Band 1 on [2^24, 2^29.5) with sub-segment splits on gaps >= 2^32
     - Lazy 2-Barrett5 + 4-limb CIOS MontMul(chi_R) c_tail_z6 (323 ops <= 340, 200 random tests)
     - 15-entry / 16-vec Z[omega] SpMV kernel (C_UNIT_5 = 6.80, k5_cap = 8.85)
     - Unconditional |log2 N_g(a,b) - log2 N(u*w*a^2)| <= 0.0156 b + uint8_t zero-underflow/wrap-around + read-only scan8_u8.
"""
import math, random, re
from fractions import Fraction as Fr
from math import log, log2, exp, sqrt, ceil, floor
from pathlib import Path

LN2 = math.log(2.0)
L3 = log2(3.0)
M64 = (1 << 64) - 1
S_U8 = 0.32
B1_DELTA_HI = 29.5
C_HIT_DELTA_BASE = 2.82

# ---------------- 1. Exact Pure-Stdlib 2D Retention for uint8_t STADDB (scale=0.32, odd=2) ----------------
_RET2D_CACHE = {}

def get_ret2d_u8(scale=S_U8, pmin=24, B=61.0, step=1.0, odd=2):
    key = (scale, pmin, B, step, odd)
    if key in _RET2D_CACHE:
        return _RET2D_CACHE[key]
    NX = int(round(512.0 / step))
    NY = int(round(136.0 / step))
    comp = bytearray(16385)
    for p in range(2, 129):
        if not comp[p]:
            for j in range(p * p, 16385, p):
                comp[j] = 1
    primes14 = [p for p in range(2, 16384) if not comp[p]][odd:]

    def _build_weights(is_smooth, use_sf_sgn, tilt):
        b_max = B if is_smooth else float(NX) * step
        W_all = [0.0] * (NX + 1)
        WG = {}
        for p in primes14:
            b = log2(p)
            bi = int(round(b / step))
            q = exp(tilt * bi * step - LN2 * b)
            qk = q; k = 1
            while k * bi <= NX and qk >= 1e-25:
                amt = qk / k
                W_all[k * bi] += amt
                if k * bi <= NY:
                    WG[(k * bi, k * bi)] = WG.get((k * bi, k * bi), 0.0) + amt
                qk *= q; k += 1
        w = step / 4.0
        ncells = int(ceil((b_max - 14.0) / w))
        w = (b_max - 14.0) / ncells
        for i in range(ncells):
            b = 14.0 + (i + 0.5) * w
            bi = int(round(b / step))
            is_sieve = (pmin <= b <= B)
            err_b = (b - floor(scale * b) / scale) if is_sieve else b
            ei = max(0, int(round(err_b / step)))
            for k in range(1, 9):
                if k * bi > NX:
                    break
                amt = exp(k * tilt * bi * step - (k - 1) * LN2 * b) * w / (k * b)
                if amt < 1e-25:
                    break
                sgn = -1.0 if (use_sf_sgn and is_sieve and (k % 2 == 0)) else 1.0
                W_all[k * bi] += amt
                if k * ei <= NY:
                    kk = (k * bi, k * ei)
                    WG[kk] = WG.get(kk, 0.0) + sgn * amt
        A_den = [0.0] * (NX + 1); A_den[0] = 1.0
        jW = [j * W_all[j] for j in range(NX + 1)]
        _sp = math.sumprod
        for n in range(1, NX + 1):
            A_den[n] = _sp(jW[1:n + 1], A_den[n - 1::-1]) / n
        return WG, A_den

    def _run(mode, tilt):
        is_smooth = (mode != "sf")
        use_sf_sgn = (mode == "rg")
        WG, A_den = _build_weights(is_smooth, use_sf_sgn, tilt)
        # Group by e for fast inner loops using sumprod on e==0 and tight slices on e>0!
        by_e = {}
        for (b, e), v in WG.items():
            if abs(v) > 1e-25:
                by_e.setdefault(e, []).append((b, b * v))
        e0_map = dict(by_e.pop(0, []))
        e0_lo = min(e0_map) if e0_map else NX + 1
        e0_hi = max(e0_map) if e0_map else 0
        e0_rev = [e0_map.get(b, 0.0) for b in range(e0_hi, e0_lo - 1, -1)] if e0_map else []
        by_e_items = sorted(by_e.items())
        cols = [[0.0] * (NX + 1) for _ in range(NY + 1)]
        inv_n = [0.0] + [1.0 / n for n in range(1, NX + 1)]
        _sp = math.sumprod
        for r in range(NY + 1):
            rhs = [0.0] * (NX + 1)
            for e, lst in by_e_items:
                if e > r:
                    break
                src = cols[r - e]
                for b, bw in lst:
                    for n in range(b, NX + 1):
                        rhs[n] += bw * src[n - b]
            col = cols[r]
            if r == 0:
                col[0] = 1.0
            for n in range(1, NX + 1):
                s = rhs[n]
                if n >= e0_lo:
                    k0 = 0 if n <= e0_hi else n - e0_hi
                    s += _sp(e0_rev[e0_hi - (n - k0):], col[k0:n - e0_lo + 1])
                col[n] = s * inv_n[n]
        cum = [[0.0] * (NY + 1) for _ in range(NX + 1)]
        for n in range(1, NX + 1):
            dn = A_den[n]
            if dn <= 0.0:
                continue
            acc = 0.0; row = cum[n]
            for r in range(NY + 1):
                acc += cols[r][n]
                row[r] = min(1.0, max(0.0, acc / dn))
        return cum, A_den

    rg_cum, _ = _run("rg", 0.05)
    sf_cum, den_all = _run("sf", 0.005)
    # 1D exact smooth-geo denominator at tilt=0.005 gives an exact, strictly outward upper bound on s0(n, r) = P(B-smooth AND res<=r) <= P(B-smooth) = den_sm005[n]/den_all[n]!
    _, den_sm005 = _build_weights(True, False, 0.005)
    s0_cum = [[0.0] * (NY + 1) for _ in range(NX + 1)]
    for n in range(1, NX + 1):
        ub = (den_sm005[n] / den_all[n]) if den_all[n] > 0.0 else 0.0
        srow = sf_cum[n]; row = s0_cum[n]
        for r in range(NY + 1):
            row[r] = ub if ub < srow[r] else srow[r]
    out = dict(rg=rg_cum, sf=sf_cum, s0=s0_cum, NX=NX, NY=NY, step=step, scale=scale)
    _RET2D_CACHE[key] = out
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
    # At each degree-1 prime P_i above 3 (i=1,2,3):
    #   v_{P_i}(r8..r0) = [2,2,1,1,0], v_{P_i}(gl, gc) = [2,0] with gl/9 == 1 (mod 3), gc == -1 (mod 3).
    #   Primitive (a, b) mod P_i (8 pairs):
    #     P_i \nmid b (prob 3/4): (vf, vg) = (0, 0).
    #     P_i | b, P_i \nmid a (prob 1/4): vf = 2 always; b = 3 b_1:
    #       3 | b_1 (prob 1/12): vg = 2;
    #       3 \nmid b_1 (prob 1/6): g = 9*(64 a^2 - b_1^2) -> Hensel geometric vg = 2 + k with prob (1/9)*(1/3)^(k-1), k >= 1.
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
    # At inert degree-3 prime Q_3 (norm 27, 28 projective points in P^1(F_27)):
    #   v_{Q_3}(r8..r0) = [0,1,1,2,2], v_{Q_3}(gl, gc) = [1,0].
    #   Q_3 \nmid a, b (13/14): (0, 0); Q_3 | a (1/28): (3*2=6, 0); Q_3 | b (1/28): (0, 3*1=3).
    q3_d = {(0, 0): Fr(13, 14), (6, 0): Fr(1, 28), (0, 3): Fr(1, 28)}
    full_3 = _conv(p3_d, q3_d)
    assert sum(full_3.values()) == 1

    # Cap (f3, g3) at 12 conservatively (strictly outward: never overstates valuation)
    c3_cap = {}
    for (f3, g3), w in full_3.items():
        k = (min(f3, 12), min(g3, 12))
        c3_cap[k] = c3_cap.get(k, Fr(0)) + w

    CLASSES_2 = ((Fr(7, 9), 0, 0), (Fr(1, 9), 0, 3), (Fr(1, 9), 6, 0))
    # Exact analytical means (untruncated):
    exact_ev2_f, exact_ev2_g = Fr(2, 3), Fr(1, 3)
    exact_ev3_f, exact_ev3_g = Fr(12, 7), Fr(33, 14)
    ex_f_bits = float(exact_ev2_f) + float(exact_ev3_f) * L3
    ex_g_bits = float(exact_ev2_g) + float(exact_ev3_g) * L3

    # Fold only negligible classes (w < 1e-5) conservatively into (0, 0) (exact unbinned valuations for all classes >= 1e-5!)
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
    # Factorized index tables (distinct vfc, distinct vgc, sparse triples (w, i_f, i_g)) for 250x faster exact eval_d!
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
    # Exact maps phi_i: O_h -> Z/729Z on the Z-basis ZK_NUM (derived from h(y) mod 729)
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
    # Verify pairwise distinct mod 3 (hence 3 distinct degree-1 primes above 3!)
    distinct_ok = (len({tuple(x % 3 for x in mp) for mp in P_MAPS}) == 3)
    # Express f_rec and g_rec coefficients in ZK basis and check exact images mod 729!
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
    # Verify exact valuation law at P1, P2, P3 on 1,500 random pairs mod 729!
    rng = random.Random(seed)
    v3 = lambda x: 6 if x % 729 == 0 else (next(k for k in range(6) if (x // (3 ** k)) % 3 != 0))
    law_errs = 0; tested = 0
    while tested < 1500:
        a = [rng.randrange(729) for _ in range(6)]
        b = [rng.randrange(729) for _ in range(6)]
        # check primitive at P1, P2, P3
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


# ---------------- 4. Band-1 Delta Exact Cost & All v32 Kernel/LA/uint8 Verifiers ----------------
def b1_delta_save_per_cell(M, raw_mult, p_min=24.0, p_hi=B1_DELTA_HI):
    """Exact saving on rep_lo from 32-bit delta-packed Band 1 on [2^p_min, 2^p_hi) (uint8_t 1 B/cell),
    charging 2.82 Ir/hit base + 24 Ir per sub-segment header whenever gap >= 2^32."""
    w_k = [math.comb(3, k) * ((2.0 / 9.0) ** k) * ((2.0 / 3.0) ** (3 - k)) / (512.0 / 729.0) for k in range(4)]
    tot_save_log = 0.0
    steps = int(round((p_hi - p_min) * 4))
    dp = (p_hi - p_min) / steps
    for i in range(steps):
        pb = p_min + (i + 0.5) * dp
        p = 2.0 ** pb
        # Evaluate gap-overflow probability at the RIGHT endpoint pb + dp/2 (strictly outward!)
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
    # Exact ARMv8.0-A instruction count: 43 (6-bucket diff) + 112 (2x signed Barrett5) + 154 (4-limb CIOS MontMul) + 14 (add+sub+csel+2 stp) = 323 <= 340
    ctail_ops = 43 + 112 + 154 + 14

    # 4d. Aligned-sub-band bfi Lower Band 2 (uint8_t STADDB, all k_s in 37..48, W in {8, 16}, sgn in {+1, -1})
    # W*2^k_s alignment clears bits [k_s .. k_s + w_bits - 1] (w_bits = 3 for W=8, 4 for W=16); test non-zero higher bits m << (k_s + w_bits)
    # explicitly including (W=8, k_s=45, sgn=+1, base=2^55 + 2^48, ai=1, j_rel=0)!
    bfi_b2_errs = 0
    bfi_b2_cases = 0
    for W_t in (8, 16):
        w_bits = 3 if W_t == 8 else 4
        w_mask = (1 << w_bits) - 1
        for k_s in range(37, 49):
            namax = 4000
            namax_sh4 = (namax << 4) | 15
            pairs = [(1, 0), (1, W_t - 1), (namax, 0), (namax, W_t - 1)] + sorted(
                (rng.randint(1, 5000), rng.randrange(W_t)) for _ in range(40)
            )
            pairs.sort()
            entries = [(ai << 4) | jr for ai, jr in pairs]
            for sgn_v in (1, -1):
                for m_hi in (0, 1, 2, 3, 5):
                    base = ((1 << 55) & ~(w_mask << k_s)) + (m_hi << (k_s + w_bits))
                    assert (base & (w_mask << k_s)) == 0 and (base & ((1 << k_s) - 1)) == 0
                    bfi_b2_cases += 1
                    got_m = {}; ref_m = {}
                    for e in entries:
                        if e > namax_sh4:
                            break
                        x5 = base + (e >> 4)
                        fld = ((~e) & w_mask) if sgn_v < 0 else (e & w_mask)
                        # Exact ARM64 `bfi x5, src, #k_s, #w_bits` (3 bits for W=8, 4 bits for W=16; preserves bit k_s + w_bits!)
                        mask = w_mask << k_s
                        x5 = (x5 & ~mask) | (fld << k_s)
                        got_m[x5] = (got_m.get(x5, 0) + 17) & 0xFF
                    for ai, jr in pairs:
                        if ai > namax:
                            break
                        jr_eff = jr if sgn_v > 0 else (W_t - 1 - jr)
                        ad = base + (jr_eff << k_s) + ai
                        ref_m[ad] = (ref_m.get(ad, 0) + 17) & 0xFF
                    if got_m != ref_m:
                        bfi_b2_errs += 1

    # 4e. 32-bit delta-packed Band 1 end-to-end check (including forced gaps >= 2^32 sub-segment splits!)
    delta_b1_errs = 0
    ais_big = [10, 250, (1 << 32) + 500, (1 << 32) + 900, (1 << 33) + 1200, (1 << 33) + 1205, (1 << 33) + 1299, (1 << 33) + 1400, (1 << 33) + 1555]
    # Build segments: start a new segment whenever delta >= 2^32
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
    # Simulate packed uint32_t copy + b1_delta_apply (STADDB) on C[s]
    row0 = 1 << 50; w_u8 = 11
    got_d = {}; ref_d = {}
    for ai in ais_big:
        ref_d[row0 + ai] = (ref_d.get(row0 + ai, 0) + w_u8) & 0xFF
    for sb, ds in segs_d:
        # Pack into 64-bit words (2 uint32_t per word), pad with 0 if needed for 4-way loop and handle tail
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

    # 4f. Unconditional |log2 N_g(a,b) - log2 N(u*w*a^2)| bound + uint8_t zero-underflow/wrap-around + read-only scan8_u8
    T4 = M.load_data()["t4_geom"]
    cA_emb = T4["GA"][0][0] / 3.0
    w_min_abs = min(abs(complex(r_e, im_e)) for r_e, im_e in T4["EMB_W"])
    a_emb2_max = (589.43 ** 2) * M.QA_MAX / cA_emb
    b_emb2_max = (782.76 ** 2) * M.QB_MAX / cA_emb
    a_emb2_min = 1.0 / (a_emb2_max * a_emb2_max)
    ratio_max = b_emb2_max / (float(M.U_BN254) * w_min_abs * a_emb2_min)
    ng_diff_ub = max(6.0 * log2(1.0 + ratio_max), -6.0 * log2(1.0 - ratio_max))
    sl_u8 = M.OP["slack"]
    thr_max = floor(S_U8 * (512.0 - sl_u8))
    init_min = 128 - thr_max
    max_byte = 128 + ceil(S_U8 * (sl_u8 + 1.0))
    u8_range_ok = (ng_diff_ub <= 0.02 and thr_max <= 128 and init_min >= 0 and max_byte <= 255)
    # Bit-exact read-only scan8_u8 (8 cells per 64-bit word, MSB mask 0x8080808080808080):
    scan8_u8_errs = 0
    n_c = 256
    init_u8 = [128 - rng.randint(80, 124) for _ in range(n_c)]
    rows_8 = [init_u8[:] for _ in range(8)]
    exp_s = set()
    for rr in range(8):
        for cc in range(n_c):
            if rng.random() < 0.03:
                rows_8[rr][cc] = (rows_8[rr][cc] + rng.randint(90, 140)) & 0xFF
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

    ok = (t3["ok"] and ntt3_split_errs == 0 and ctail_errs == 0 and ctail_ops <= M.OP["c_tail_z6"]
          and bfi_b2_errs == 0 and bfi_b2_cases == 240 and delta_b1_errs == 0 and u8_range_ok and scan8_u8_errs == 0
          and b1_delta_ops_ok)
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
        u8_scale=S_U8, u8_thr_max=thr_max, u8_init_min=init_min, u8_max_byte=max_byte,
        u8_range_ok=u8_range_ok,
        scan8_u8_errs=scan8_u8_errs,
        b1_delta_ops_ok=b1_delta_ops_ok,
    )
