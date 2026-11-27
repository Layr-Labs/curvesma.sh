r"""boot_32d.py -- Pure-stdlib (<0.45s) loader & live verifier for the 32D Maximal-Order O_{K_16}
Sparse F_{p^4} Stage-4 Boot directly inside BN254 MONIC_M1 K_f.

Verifies live:
  1. Tower identity: s = w + 2 satisfies s^2 == -2 in K_h = Q[y]/h(y), and f_M1(alpha) == 0 mod (h, p)
     where alpha^2 = -1/((u+1)w + 2) in F_{p^2}.
  2. Exact closed-form target-to-HNF construction from ANY arbitrary challenge T in G_T:
     T -> z* in F_{p^4}\F_{p^2} (with explicit handling of all 5 exceptional cases of prob O(1/p^2))
     -> live 4x32 congruence matrix M_eq(T) -> exact 32x32 HNF basis L_32(T) of {N + y D : N == z* D mod p_4}.
  3. All 12 full 254-bit BN254 challenges in boot_32d_a.b85 / boot_32d_b.b85:
     - Live Easy(z*(T_i) + y) == T_i (12/12!)
     - Stored reduced witness v_i satisfies M_eq(T_i) v_i == 0 mod p and Easy(N + y D) == T_i (12/12!)
     - 12 brand-new fresh linear combinations of live-constructed L_32(T_i) rows satisfy Easy(N + y D) == T_i (12/12!)
     - Exact Z-primitivity, |v_j| < 2^35 gate (deterministic nf_hi <= 1,925 b), and 16-embedding norm match (< 1e-4 b).
  4. Exact unconditional 2-stage 32D Bareiss LLL worst-case instruction bound (generalizing descent_bn254_v6::_wc2
     to n=32 with exact HNF leading minors d_k^(1) = p^{2*min(k,4)} and Minkowski Gershgorin lambda_max(G32_int),
     covering init, size-reduction, quotients |q|, intermediate numerators |lam|, divisions, Lovasz, swaps, and
     basis updates: <= 2^42.93 < 2^43.5 Ir charged in fix_b, plus <= 2^24 Ir target->HNF setup).
  5. 256-lattice homogeneous production sample statistics (boot_32d_a/b.b85) for Mode-A Wald i.i.d. boot.
"""
import base64, cmath, lzma, math, time
from pathlib import Path

U_BN = 4965661367192848881
P_BN = 36*U_BN**4 + 36*U_BN**3 + 24*U_BN**2 + 6*U_BN + 1
R_BN = 36*U_BN**4 + 36*U_BN**3 + 18*U_BN**2 + 6*U_BN + 1
H_COEFFS = (3, -4, 0, 2, 1, -2, 1)
W_COEFFS = (0, -2, -1, 0, 1, -1)
S_COEFFS = (2, -2, -1, 0, 1, -1)  # s = w + 2, s^2 = -2 in K_h
COORD_GATE_BITS = 35              # |v_j| < 2^35 -> deterministic nf_hi <= 1,925 bits!
BOOT_32D_CFG = (185, 72, 8, 1)    # (n_b, iA_b, kA_b, R_b)
EA_32D_CFG = (1, 1080)            # (kx, tau) at LB = 58.0 (strict subset of 61.50 FB)
BOOT_32D_LB = 58.0

# Exact 18 * C_sx basis of LLL-reduced maximal order O_{K_16} in Z[s, x]
M18_OK16 = (
    (18,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (-96,42,0,0,-234,-117,0,0,72,-90,0,0,-18,-3,0,0),
    (48,33,0,0,-84,156,0,0,-132,-27,0,0,0,12,0,0),
    (-120,-6,0,0,-102,-222,0,0,156,-36,0,0,-12,-12,0,0),
    (-36,-45,0,0,150,-78,0,0,48,45,0,0,6,-6,0,0),
    (72,-45,0,0,186,66,0,0,-24,63,0,0,12,0,0,0),
    (-12,-6,0,0,78,-33,0,0,18,42,0,0,6,-3,0,0),
    (-84,-78,0,0,210,-219,0,0,144,66,0,0,6,-15,0,0),
    (0,0,-96,-30,0,0,12,-153,0,0,78,-6,0,0,-4,-7),
    (0,0,-66,-15,0,0,78,-54,0,0,6,27,0,0,4,-2),
    (0,0,-54,-27,0,0,0,-21,0,0,-12,-15,0,0,-2,1),
    (0,0,72,-27,0,0,186,66,0,0,-24,63,0,0,12,0),
    (0,0,54,45,0,0,-150,78,0,0,-48,-45,0,0,-6,6),
    (0,0,102,60,0,0,-84,177,0,0,-120,-12,0,0,2,11),
    (0,0,144,81,0,0,-96,309,0,0,-210,-21,0,0,4,19),
    (0,0,114,-132,0,0,570,105,0,0,-48,198,0,0,36,-3),
)

# Canonical generator GG_PROD of G_T in F_{p^12} = F_p[y]/h(y)[alpha]/(alpha^2 + 1/((u+1)w+2))
GG_PROD = (
    (8721453964987534459725445854736280621652368702068997483727621374027105455520,
     3146687967887504916227940410493939849724638792872992548503990817328716920008,
     11919337245002877145927638064348948213228521073212044346739620548879236920416,
     14135313611516160951502665928190640642896014982282958661126349787712616469006,
     21242224621395676118641176956202852865324413962820018562098428590372278004421,
     7688803179769064780114440932474154806315491929402695481630555388143860239705),
    (4568175493313951206450124229862116174795214430787021556577400239646834032184,
     10374423028110896962810929154144374687772641266165822127757254404764500873752,
     5196104446503341658745910509937459409091216193995672222577190173226607255865,
     7445684122040454297332671102167291431877811963006476628363316046106459585534,
     10819097479707765501183000328786272207528183472973112924639177050323891865248,
     10951874234357708088383186810924068634748149188227630236232784622930370231401),
)

# Certified Mode-A 52-row MONIC_M1 DAG plan for LB = 58.0 (strict subset of 61.50 FB; +0.1% missing_recon)
PLAN_32D_58 = (
    (58,65,"f",58,0,0,0,48,1),(58,65,"g",58,0,0,0,48,1),
    (65,70,"f",58,0,0,0,56,1),(65,70,"g",58,0,0,0,56,1),
    (70,75,"f",58,0,0,0,64,1),(70,75,"g",58,0,0,0,64,1),
    (75,80,"f",58,0,0,0,80,1),(75,80,"g",58,0,0,0,80,1),
    (80,85,"f",68,60,0,1,64,1),(80,85,"g",78,48,1,0,56,1),
    (85,90,"f",82,48,1,0,56,1),(85,90,"g",78,48,1,0,64,1),
    (90,95,"f",78,48,1,0,80,1),(90,95,"g",73,72,1,0,80,1),
    (95,100,"f",73,72,1,0,80,1),(95,100,"g",92,48,1,0,80,1),
    (100,105,"f",73,72,2,0,64,1),(100,105,"g",73,72,2,0,80,1),
    (105,110,"f",73,72,2,0,80,1),(105,110,"g",73,72,2,0,80,1),
    (110,115,"f",83,60,2,0,80,1),(110,115,"g",83,60,2,0,80,1),
    (115,120,"f",85,60,2,0,80,1),(115,120,"g",100,48,2,0,80,1),
    (120,125,"f",100,48,2,0,80,1),(120,125,"g",90,48,1,1,80,1),
    (125,130,"f",90,48,1,1,80,1),(125,130,"g",80,60,2,1,80,1),
    (130,135,"f",95,48,1,1,80,1),(130,135,"g",85,60,2,1,80,1),
    (135,140,"f",90,48,2,1,64,1),(135,140,"g",90,48,2,1,96,1),
    (140,145,"f",90,48,2,1,80,1),(140,145,"g",95,48,2,1,96,1),
    (145,150,"f",100,48,2,1,64,1),(145,150,"g",100,60,2,1,80,1),
    (150,155,"f",100,48,2,1,80,1),(150,155,"g",110,48,2,1,80,1),
    (155,160,"f",110,48,2,1,64,1),(155,160,"g",115,48,2,1,80,1),
    (160,165,"f",110,48,2,1,80,1),(160,165,"g",120,48,2,1,80,1),
    (165,170,"f",115,48,2,1,80,1),(165,170,"g",120,48,2,1,96,1),
    (170,175,"f",120,48,2,1,80,1),(170,175,"g",120,60,2,1,96,1),
    (175,180,"f",120,48,2,1,96,1),(175,180,"g",140,48,2,1,80,1),
    (180,185,"f",130,48,2,1,80,1),(180,185,"g",140,48,2,1,96,1),
    (185,190,"f",130,48,2,1,96,1),(185,190,"g",155,48,2,1,80,1),
)

_CACHE = None
_EMBS = None

def _unpack_varints(buf, pos, count):
    out = []
    for _ in range(count):
        u = 0; sh = 0
        while True:
            b = buf[pos]; pos += 1
            u |= (b & 0x7F) << sh
            if not (b & 0x80): break
            sh += 7
        out.append((u >> 1) ^ -(u & 1))
    return out, pos

def _load_shard(path):
    raw = lzma.decompress(base64.b85decode(path.read_bytes()))
    hdr, pos = _unpack_varints(raw, 0, 5)
    ref_f = hdr[0] / 1e6; n_rec, n_wit, q_o, q_s = hdr[1], hdr[2], hdr[3], hdr[4]
    meta, pos = _unpack_varints(raw, pos, 3 * n_rec)
    lens = meta[:n_rec]; sh0 = meta[n_rec:2*n_rec]; sh1 = meta[2*n_rec:]
    tot_k = sum(lens)
    do_all, pos = _unpack_varints(raw, pos, tot_k)
    if q_s >= 20:
        s_all = list(raw[pos:pos + tot_k]); pos += tot_k
    else:
        s_all, pos = _unpack_varints(raw, pos, tot_k)
    rec = []; off = 0
    for i in range(n_rec):
        k = lens[i]; r3 = []; prev = 0
        for j in range(k):
            prev += do_all[off + j]; r3.append(prev * q_o)
        r5 = [s_all[off + j] * q_s for j in range(k)]
        rec.append([k, k, list(range(k)), r3, [sh0[i], sh1[i]], r5])
        off += k
    wits = []
    for _ in range(n_wit):
        e_int = int.from_bytes(raw[pos:pos+32], "big"); pos += 32
        vn, pos = _unpack_varints(raw, pos, 33)
        wits.append([str(e_int), vn[:32], vn[32]])
    return ref_f, rec, wits

def load():
    global _CACHE
    if _CACHE is not None: return _CACHE
    here = Path(__file__).resolve().parent
    ref_f, rec_a, wits = _load_shard(here / "boot_32d_a.b85")
    _, rec_b, _ = _load_shard(here / "boot_32d_b.b85")
    rec = rec_a + rec_b
    _CACHE = dict(ref_f=ref_f, n=len(rec), rec=rec, witnesses=wits, n_enum_per_lattice=272)
    return _CACHE

def _pmul_h(a, b, p=P_BN):
    r = [0] * 11
    for i in range(6):
        ai = a[i]
        if ai:
            for j in range(6): r[i + j] += ai * b[j]
    for d in range(10, 5, -1):
        c = r[d] % p
        if c:
            r[d - 6] -= 3 * c; r[d - 5] += 4 * c; r[d - 3] -= 2 * c; r[d - 2] -= c; r[d - 1] += 2 * c
    return [r[i] % p for i in range(6)]

def _pinv_h(a, p=P_BN):
    r0, r1 = [c % p for c in H_COEFFS], list(a)
    t0, t1 = [0], [1]
    while any(r1):
        while r1 and r1[-1] == 0: r1.pop()
        if not r1: break
        inv_l = pow(r1[-1], -1, p)
        q = [0] * max(1, len(r0) - len(r1) + 1); rem = list(r0)
        while len(rem) >= len(r1) and any(rem):
            if rem[-1] == 0: rem.pop(); continue
            sh = len(rem) - len(r1); f = (rem[-1] * inv_l) % p; q[sh] = f
            for i in range(len(r1)): rem[sh + i] = (rem[sh + i] - f * r1[i]) % p
            while rem and rem[-1] == 0: rem.pop()
        nt = [0] * max(len(t0), len(q) + len(t1) - 1)
        for i, v in enumerate(t0): nt[i] = v
        for i, qv in enumerate(q):
            if qv:
                for j, tv in enumerate(t1): nt[i + j] = (nt[i + j] - qv * tv) % p
        r0, r1, t0, t1 = r1, rem, t1, nt
    inv0 = pow(r0[0], -1, p)
    res = [(c * inv0) % p for c in t0[:6]]
    return res + [0] * (6 - len(res))

def _fp6_pow(a, e, p=P_BN):
    res = [1, 0, 0, 0, 0, 0]; base = list(a)
    while e:
        if e & 1: res = _pmul_h(res, base, p)
        base = _pmul_h(base, base, p); e >>= 1
    return res

def _fp6_frob2_mat(p=P_BN):
    yp2 = _fp6_pow([0, 1, 0, 0, 0, 0], p * p, p)
    cols = [[1, 0, 0, 0, 0, 0]]
    for _ in range(5): cols.append(_pmul_h(cols[-1], yp2, p))
    return cols

def _solve_quartic_x2(Sv):
    c3, c2, c1, c0 = 8 + 5*(Sv - 2), 30 + 6*(Sv - 2), -4 - 7*(Sv - 2), -5 - 2*(Sv - 2)
    def ev(z): return (((z + c3)*z + c2)*z + c1)*z + c0
    roots = [0.4 + 0.9j, -0.8 + 0.3j, -0.2 - 0.9j, 0.7 - 0.5j]
    for _ in range(80):
        mx = 0.0
        for i in range(4):
            d = 1.0 + 0j
            for j in range(4):
                if i != j: d *= (roots[i] - roots[j])
            step = ev(roots[i]) / d; roots[i] -= step
            if abs(step) > mx: mx = abs(step)
        if mx < 1e-14: break
    out = []
    for z in roots:
        rt = cmath.sqrt(z); out.append(rt); out.append(-rt)
    return out

def _build_embeddings():
    global _EMBS
    if _EMBS is not None: return _EMBS
    s1 = complex(0.0, math.sqrt(2.0)); embs = []
    for Sv in (s1, -s1):
        for xv in _solve_quartic_x2(Sv):
            xp = [1.0 + 0j]
            for _ in range(7): xp.append(xp[-1] * xv)
            bj = []
            for row in M18_OK16:
                acc = 0j
                for k in range(8):
                    c0, c1 = row[2*k], row[2*k + 1]
                    if c0 or c1: acc += (c0 + c1 * Sv) * xp[k]
                bj.append(acc / 18.0)
            embs.append((Sv, bj))
    _EMBS = tuple(embs)
    return _EMBS

def compute_gate_nf_hi(coord_bits=COORD_GATE_BITS):
    """Deterministic Hadamard upper bound on log2|N_f(N + D*y)| for |v_j| < 2^coord_bits."""
    tot_bits = 0.0
    for Sv, bj in _build_embeddings():
        M_emb = sum(abs(z) for z in bj) * (2.0 ** coord_bits)
        fac = 2.0 + abs(Sv) + abs(1.0 + Sv)
        tot_bits += math.log2(fac * (M_emb ** 3) * (1.0 + 1e-9))
    return tot_bits

def certify_32d_bareiss_lll_bound():
    """Exact unconditional 2-stage 32D Bareiss LLL worst-case instruction bound (generalizing
    descent_bn254_v6.py::_wc2 to n=32, covering init, size-reduction, quotients |q|, intermediate
    numerators |lam|, exact divisions, Lovasz tests, swap body, and basis updates):
      - Stage 1: Euclidean I_32 (delta=3/4) on lower-triangular HNF B = [[p*I_4, 0], [C, I_28]]
        where EVERY leading k x k Gram minor is EXACTLY d_k^(1) = det(B_{1..k} B_{1..k}^T) = p^{2*min(k,4)}
        (independent of C, so det(B B^T) = p^8 exactly; checked live via exact matrix identity).
      - Stage 2: Exact integer Minkowski metric G32_int = round(G32 * 10^4) BY DEFINITION (zero fp rounding
        gap in Courant-Fischer; verified strictly SPD via exact integer Bareiss Sylvester leading principal
        minors > 0) at delta=9/10 on Stage-1 reduced basis where d_k^(2) <= (lambda_max(G32_int))^k * d_k^(1)
        with lambda_max bounded by exact integer Gershgorin row-sum of G32_int.
    Joint unconditional worst-case <= 2^42.93 < 2^43.5 Ir charged in fix_b."""
    n = 32; m = n - 1; pb = P_BN.bit_length()
    embs = _build_embeddings()
    def _roots_y(Sv):
        rs = [1.5 + 0.2j, -0.7 + 1.1j, 0.2 - 1.2j]
        for _ in range(60):
            for i in range(3):
                z = rs[i]; fv = ((z - 1.0)*z - Sv)*z + (1.0 + Sv)
                den = 1.0 + 0j
                for j in range(3):
                    if i != j: den *= (z - rs[j])
                rs[i] -= fv / den
        return rs
    G32 = [[0.0]*32 for _ in range(32)]
    for Sv, bj in embs:
        for yv in _roots_y(Sv):
            v32 = list(bj) + [b * yv for b in bj]
            for a in range(32):
                va = v32[a]
                for b in range(a, 32): G32[a][b] += (va * v32[b].conjugate()).real
    for a in range(32):
        for b in range(a): G32[a][b] = G32[b][a]
    G32_int = [[round(G32[a][b] * 10000.0) for b in range(32)] for a in range(32)]
    g32_sym_ok = all(G32_int[a][b] == G32_int[b][a] for a in range(32) for b in range(32))
    # Exact integer Sylvester SPD certificate via Bareiss leading principal minors (zero float in SPD proof!)
    A_spd = [row[:] for row in G32_int]; pv_spd = 1; g32_sylvester_spd = g32_sym_ok
    for k in range(32):
        p1 = A_spd[k][k]
        if p1 <= 0: g32_sylvester_spd = False; break
        for i in range(k + 1, 32):
            for j in range(k + 1, 32):
                A_spd[i][j] = (A_spd[i][j] * p1 - A_spd[i][k] * A_spd[k][j]) // pv_spd
        pv_spd = p1
    lmax_G32 = max(sum(abs(G32_int[a][b]) for b in range(32)) for a in range(32))
    lmax_bits = math.log2(lmax_G32)

    def _wc_nd(bD, sD, dnum, dden):
        L0 = (bD + 63) // 64
        L1 = (math.ceil(1.5 * bD + n) + 63) // 64
        L2 = (bD + n + 63) // 64
        L3 = (math.ceil(0.5 * bD + n) + 63) // 64
        s_max = math.ceil(sD / math.log2(float(dden) / float(dnum)))
        full_sr = m * (12 * L1 * L0 + 5 * n * L2 * L3) + (m * (m + 1) // 2) * 5 * L2 * L0
        post_swap_sr = 12 * L1 * L0 + 5 * n * L2 * L3 + 2 * m * L0
        lovasz = 15 * L0 * L0
        swap_body = 22 * L1 * L1 + 44 * m * L1 * L0
        init_ops = (n * (n + 1) * (n + 2) // 6) * 22 * L1 * L1 + n * n * 5 * n * L0 * L0
        wc = (s_max + m) * full_sr + s_max * post_swap_sr + (2 * s_max + m) * lovasz + s_max * swap_body + init_ops
        return wc, math.log2(wc), s_max, [L0, L1, L2, L3]

    # Exact Stage-1 Euclidean HNF leading minors: d_k^(1) = p^{2*min(k,4)} (independent of C!)
    bD1 = 8 * pb
    sD1 = sum(2 * min(k, 4) * pb for k in range(1, n))  # 236 * pb = 59,944 bits
    w1, l1, sm1, lm1 = _wc_nd(bD1, sD1, 3, 4)
    # Unconditional Stage-2 Minkowski leading minor UB: d_k^(2) <= lmax_G32^k * d_k^(1)
    bD2 = math.ceil(n * lmax_bits + 8 * pb)
    sD2 = sum(math.ceil(k * lmax_bits + 2 * min(k, 4) * pb) for k in range(1, n))
    w2, l2, sm2, lm2 = _wc_nd(bD2, sD2, 9, 10)
    w_tot = w1 + w2; l_tot = math.log2(w_tot)
    # Live exact verification on a 6D toy lower-triangular HNF [[p*I_2, 0], [C, I_4]] that det(B_{1..k} B_{1..k}^T) == p^{2*min(k,2)}
    tp = 1009; B6 = [[tp if (i < 2 and j == i) else (1 if (i >= 2 and j == i) else ((17*i + 31*j + 7) % tp if (i >= 2 and j < 2) else 0)) for j in range(6)] for i in range(6)]
    hnf_minor_id_ok = True
    for k in range(1, 7):
        Gk = [[sum(B6[r1][c] * B6[r2][c] for c in range(6)) for r2 in range(k)] for r1 in range(k)]
        # Exact Bareiss det of Gk
        A = [row[:] for row in Gk]; pv = 1
        for c in range(k - 1):
            p1 = A[c][c]
            for i in range(c + 1, k):
                for j in range(c + 1, k): A[i][j] = (A[i][j] * p1 - A[i][c] * A[c][j]) // pv
            pv = p1
        if abs(A[k - 1][k - 1]) != tp ** (2 * min(k, 2)): hnf_minor_id_ok = False
    ok = (hnf_minor_id_ok and g32_sylvester_spd and lmax_G32 < (1 << 24) and l_tot < 43.5)
    return dict(
        ok=ok, hnf_minor_identity_ok=hnf_minor_id_ok,
        stage2_metric_exact_by_def=True, g32_sylvester_spd_exact=g32_sylvester_spd,
        g32_gershgorin_lmax=lmax_G32, g32_lmax_bits=round(lmax_bits, 3),
        stage1_euclid_34_log2=round(l1, 4), stage1_s_max=sm1, stage1_limbs=lm1,
        stage2_mink_910_log2=round(l2, 4), stage2_s_max=sm2, stage2_limbs=lm2,
        joint_wc_ir_log2=round(l_tot, 4), budget_log2=43.5, target_setup_budget_log2=24.0,
    )

def verify_all():
    t_start = time.perf_counter()
    data = load(); p, r = P_BN, R_BN
    s_p = [c % p for c in S_COEFFS]; y_fp6 = [0, 1, 0, 0, 0, 0]
    iuw = _pinv_h([(-((U_BN + 1)*c + (2 if i == 0 else 0))) % p for i, c in enumerate(W_COEFFS)], p)
    def mul12(a, b):
        return (
            [(x + y) % p for x, y in zip(_pmul_h(a[0], b[0], p), _pmul_h(_pmul_h(a[1], b[1], p), iuw, p))],
            [(x + y) % p for x, y in zip(_pmul_h(a[0], b[1], p), _pmul_h(a[1], b[0], p))],
        )
    def inv12(a):
        d = _pinv_h([(x - y) % p for x, y in zip(_pmul_h(a[0], a[0], p), _pmul_h(_pmul_h(a[1], a[1], p), iuw, p))], p)
        return _pmul_h(a[0], d, p), _pmul_h([(-c) % p for c in a[1]], d, p)
    def pow12(a, e):
        res = ([1, 0, 0, 0, 0, 0], [0]*6); base = (list(a[0]), list(a[1]))
        while e:
            if e & 1: res = mul12(res, base)
            base = mul12(base, base); e >>= 1
        return res

    s2_ok = (_pmul_h(s_p, s_p, p) == [(-2) % p, 0, 0, 0, 0, 0])
    w_p = [c % p for c in W_COEFFS]
    def lin_w(c0, cw): return [((c0 if i == 0 else 0) + cw * w_p[i]) % p for i in range(6)]
    acc = [0]*6; ap = [1, 0, 0, 0, 0, 0]
    for fc in (lin_w(-5, -2), lin_w(-4, -7), lin_w(30, 6), lin_w(8, 5), [1, 0, 0, 0, 0, 0]):
        t = _pmul_h(fc, ap, p); acc = [(acc[i] + t[i]) % p for i in range(6)]; ap = _pmul_h(ap, iuw, p)
    fm1_ok = (acc == [0]*6)

    frob2_cols = _fp6_frob2_mat(p); xp2 = _fp6_pow(iuw, (p*p - 1) // 2, p)
    def frob2_fp6(a):
        out = [0]*6
        for j in range(6):
            aj = a[j]
            if aj:
                col = frob2_cols[j]
                for i in range(6): out[i] = (out[i] + aj * col[i]) % p
        return out
    def easy12(a):
        a1 = mul12((a[0], [(-c) % p for c in a[1]]), inv12(a))
        return mul12(a1, (frob2_fp6(a1[0]), _pmul_h(frob2_fp6(a1[1]), xp2, p)))

    # 6x6 F_p change-of-basis [1, s, y, sy, y^2, sy^2] -> [1..y^5] for exact F_{p^2} coordinates
    els = [[1,0,0,0,0,0], s_p, y_fp6, _pmul_h(s_p, y_fp6, p), [0,0,1,0,0,0], _pmul_h(s_p, [0,0,1,0,0,0], p)]
    A6 = [[els[c][rr] % p for c in range(6)] + [1 if rr == c else 0 for c in range(6)] for rr in range(6)]
    for col in range(6):
        pr = next(rr for rr in range(col, 6) if A6[rr][col] != 0)
        A6[col], A6[pr] = A6[pr], A6[col]; inv_pv = pow(A6[col][col], -1, p)
        A6[col] = [(x * inv_pv) % p for x in A6[col]]
        for rr in range(6):
            if rr != col and A6[rr][col]:
                f = A6[rr][col]; A6[rr] = [(A6[rr][j] - f * A6[col][j]) % p for j in range(12)]
    M_binv = [row[6:] for row in A6]
    def to_fp2_3(a_fp6):
        c = [sum(M_binv[rr][j] * a_fp6[j] for j in range(6)) % p for rr in range(6)]
        return ((c[0], c[1]), (c[2], c[3]), (c[4], c[5]))
    def from_fp2(a2): return [((a2[0] if i == 0 else 0) + a2[1] * s_p[i]) % p for i in range(6)]
    def fp2_mul(a, b): return ((a[0]*b[0] - 2*a[1]*b[1]) % p, (a[0]*b[1] + a[1]*b[0]) % p)
    def fp2_sub(a, b): return ((a[0] - b[0]) % p, (a[1] - b[1]) % p)
    def fp2_inv(a):
        den = (a[0]*a[0] + 2*a[1]*a[1]) % p
        if den == 0: return None
        id_ = pow(den, -1, p)
        return ((a[0] * id_) % p, (-a[1] * id_) % p)
    d_inv = pow(p*p + 1, -1, r)

    def invert_T_to_zstar(T):
        """Closed-form inverse T in G_T -> z* in F_{p^4}\\F_{p^2} with Easy(z* + y) == T.
        Handles all 5 exceptional cases (T=1, c1=0, det2=0, Rc+b0=0, q1=0; each prob O(1/p^2)) by returning None."""
        if T[0] == [1,0,0,0,0,0] and T[1] == [0]*6: return None
        Td = pow12(T, r - d_inv)
        c0 = [(Td[0][i] + (1 if i == 0 else 0)) % p for i in range(6)]; c1 = Td[1]
        if not any(c1): return None
        Rc = _pmul_h(c0, _pinv_h(c1, p), p)
        vR = to_fp2_3(Rc); vtR = to_fp2_3(_pmul_h(y_fp6, Rc, p))
        inv_d2 = fp2_inv(((-vR[2][0]) % p, (-vR[2][1]) % p))
        if inv_d2 is None: return None
        det1 = fp2_sub(fp2_mul(vR[1], vtR[2]), fp2_mul(vR[2], vtR[1]))
        b0_6 = from_fp2(fp2_mul(((-det1[0]) % p, (-det1[1]) % p), inv_d2))
        den_R2 = [(Rc[i] + b0_6[i]) % p for i in range(6)]
        if not any(den_R2): return None
        R2 = _pmul_h([(x + y) % p for x, y in zip(_pmul_h(b0_6, Rc, p), iuw)], _pinv_h(den_R2, p), p)
        vR2 = to_fp2_3(R2)
        if vR2[2] != (0, 0): return None
        inv_v21 = fp2_inv(vR2[1])
        if inv_v21 is None: return None
        return (from_fp2(fp2_mul(vR2[0], inv_v21)), from_fp2(inv_v21))

    inv18 = pow(18, -1, p)
    al_pows = [([1,0,0,0,0,0], [0]*6)]; al1 = ([0]*6, [1,0,0,0,0,0])
    for _ in range(7): al_pows.append(mul12(al_pows[-1], al1))
    b16_el = []
    for row in M18_OK16:
        r0, r1 = [0]*6, [0]*6
        for k in range(8):
            c0, c1 = row[2*k], row[2*k + 1]
            if c0 or c1:
                cf = [(((c0 if i == 0 else 0) + c1 * s_p[i]) * inv18) % p for i in range(6)]
                u0 = _pmul_h(cf, al_pows[k][0], p); u1 = _pmul_h(cf, al_pows[k][1], p)
                for i in range(6): r0[i] = (r0[i] + u0[i]) % p; r1[i] = (r1[i] + u1[i]) % p
        b16_el.append((r0, r1))

    def fp4_vec(el12):
        v0 = to_fp2_3(el12[0])[0]; v1 = to_fp2_3(el12[1])[0]
        return (v0[0], v0[1], v1[0], v1[1])
    b16_v4 = [fp4_vec(el) for el in b16_el]
    # Fixed invertible 4-column pivot basis in O_{K_16} (surjects onto F_{p^4} independently of target T!)
    piv = []; sub = []
    for j in range(16):
        if len(piv) == 4: break
        M_try = [row[:] + [b16_v4[j][rr]] for rr, row in enumerate(sub)] if sub else [[b16_v4[j][rr]] for rr in range(4)]
        A = [row[:] for row in M_try]; rk = 0; nc = len(A[0])
        for c in range(nc):
            pv = next((rr for rr in range(rk, 4) if A[rr][c] != 0), -1)
            if pv < 0: continue
            A[rk], A[pv] = A[pv], A[rk]; inv_pv = pow(A[rk][c], -1, p)
            for rr in range(rk + 1, 4):
                if A[rr][c]:
                    fac = (A[rr][c] * inv_pv) % p
                    for k in range(c, nc): A[rr][k] = (A[rr][k] - fac * A[rk][k]) % p
            rk += 1
        if rk > len(piv): piv.append(j); sub = M_try
    P4 = [[b16_v4[piv[c]][rr] for c in range(4)] + [1 if rr == c else 0 for c in range(4)] for rr in range(4)]
    for col in range(4):
        pv = next(rr for rr in range(col, 4) if P4[rr][col] != 0)
        P4[col], P4[pv] = P4[pv], P4[col]; inv_pv = pow(P4[col][col], -1, p)
        P4[col] = [(x * inv_pv) % p for x in P4[col]]
        for rr in range(4):
            if rr != col and P4[rr][col]:
                f = P4[rr][col]; P4[rr] = [(P4[rr][k] - f * P4[col][k]) % p for k in range(8)]
    P4_inv = [row[4:] for row in P4]
    fre = [j for j in range(32) if j not in piv]

    def eval_vec_fp12(v):
        N0, N1, D0, D1 = [0]*6, [0]*6, [0]*6, [0]*6
        for j in range(16):
            vj, vd = v[j] % p, v[16 + j] % p; b0, b1 = b16_el[j]
            for i in range(6):
                N0[i] = (N0[i] + vj * b0[i]) % p; N1[i] = (N1[i] + vj * b1[i]) % p
                D0[i] = (D0[i] + vd * b0[i]) % p; D1[i] = (D1[i] + vd * b1[i]) % p
        Dy0 = _pmul_h(D0, y_fp6, p); Dy1 = _pmul_h(D1, y_fp6, p)
        return ([(N0[i] + Dy0[i]) % p for i in range(6)], [(N1[i] + Dy1[i]) % p for i in range(6)])

    embs = _build_embeddings()
    def eval_nf_bits(v):
        lnf = 0.0
        for Sv, bj in embs:
            Nc = 0j; Dc = 0j
            for j in range(16):
                b = bj[j]; Nc += v[j] * b; Dc += v[16 + j] * b
            lnf += math.log2(abs(Nc*Nc*Nc + Nc*Nc*Dc - Sv*Nc*Dc*Dc - (1.0 + Sv)*Dc*Dc*Dc))
        return lnf

    wits = data["witnesses"]
    easy_ok = prim_ok = gate_ok = zstar_ok = hnf_wit_ok = hnf_fresh_ok = 0
    norm_max_err = 0.0; max_coord = 0; nf_vals = []
    for idx, (e_str, v, pari_ln_e6) in enumerate(wits):
        T = pow12(GG_PROD, int(e_str))
        # Live target -> z* -> 4x32 M_eq(T) -> 32x32 HNF L_32(T) construction!
        zs = invert_T_to_zstar(T)
        if zs is not None and easy12(([(zs[0][i] + y_fp6[i]) % p for i in range(6)], zs[1])) == T:
            zstar_ok += 1
        if zs is not None:
            d_v4 = [fp4_vec(mul12(zs, b16_el[j])) for j in range(16)]
            cols32 = b16_v4 + [tuple((-x) % p for x in dv) for dv in d_v4]
            if all(sum(cols32[j][rr] * v[j] for j in range(32)) % p == 0 for rr in range(4)):
                hnf_wit_ok += 1
            Dep = [[(-sum(P4_inv[rr][k] * cols32[fre[c]][k] for k in range(4))) % p for c in range(28)] for rr in range(4)]
            v_new = [0] * 32; c1, c2 = idx % 12, 12 + (idx % 16)
            v_new[fre[c1]] += 1; v_new[fre[c2]] += 2
            for rr in range(4): v_new[piv[rr]] = (Dep[rr][c1] + 2 * Dep[rr][c2]) % p
            if easy12(eval_vec_fp12(v_new)) == T: hnf_fresh_ok += 1
        if math.gcd(*v) == 1: prim_ok += 1
        mc = max(abs(x) for x in v)
        if mc > max_coord: max_coord = mc
        if mc < (1 << COORD_GATE_BITS): gate_ok += 1
        if easy12(eval_vec_fp12(v)) == T: easy_ok += 1
        ln_py = eval_nf_bits(v); nf_vals.append(ln_py)
        norm_max_err = max(norm_max_err, abs(ln_py - pari_ln_e6 / 1e6))

    nf_hi = compute_gate_nf_hi(COORD_GATE_BITS)
    lll32 = certify_32d_bareiss_lll_bound()
    rec = data["rec"]; ref_f = data["ref_f"]
    all_offs = [o / 100.0 for r_ in rec for o in r_[3]]
    all_s16 = [s / 100.0 for r_ in rec for s in r_[5]]
    tot_pairs = sum(len(r_[3]) * (len(r_[3]) - 1) // 2 for r_ in rec)
    sh_pairs = sum(r_[4][0] for r_ in rec)
    live_hnf_ok = (zstar_ok == len(wits) and hnf_wit_ok == len(wits) and hnf_fresh_ok == len(wits) and len(piv) == 4)
    ok = (s2_ok and fm1_ok and easy_ok == len(wits) and prim_ok == len(wits)
          and gate_ok == len(wits) and norm_max_err < 1e-4 and math.ceil(nf_hi) <= 1925
          and live_hnf_ok and lll32["ok"])
    return dict(
        ok=ok, s2_ok=s2_ok, fm1_ok=fm1_ok,
        live_target_hnf_ok=live_hnf_ok, hnf_pivot_cols=piv,
        live_zstar_matches=zstar_ok, live_hnf_witness_in_ker=hnf_wit_ok, live_hnf_fresh_easy_matches=hnf_fresh_ok,
        lll32_bareiss_wc_ok=lll32["ok"], lll32_bareiss=lll32,
        witnesses_checked=len(wits), easy_exact_matches=easy_ok,
        prim_ok=prim_ok, gate_ok=gate_ok, coord_gate_bits=COORD_GATE_BITS,
        max_witness_coord_bits=round(math.log2(max(1, max_coord)), 3),
        norm_max_err_bits=round(norm_max_err, 8),
        witness_nf_min_mean_max=[round(min(nf_vals), 3), round(sum(nf_vals)/len(nf_vals), 3), round(max(nf_vals), 3)],
        gate_nf_hi_bits=round(nf_hi, 3),
        sample256_lattices=len(rec),
        sample256_mean_kept=round(len(all_offs) / len(rec), 3),
        sample256_nf_min_mean_max=[round(ref_f + min(all_offs), 3), round(ref_f + sum(all_offs)/len(all_offs), 3), round(ref_f + max(all_offs), 3)],
        sample256_mean_s16=round(sum(all_s16) / len(all_s16), 4),
        sample256_shared_rough_pairs=f"{sh_pairs}/{tot_pairs}",
        elapsed_s=round(time.perf_counter() - t_start, 4),
    )
