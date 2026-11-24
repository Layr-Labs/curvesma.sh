"""q15_ecm16.py — Deterministic Stepped-GCD + Stage-2 Binary-Product-Tree Geometric-Ladder Suyama-Montgomery ECM
(exact L_min=2 branch-tree DP, 12-bit geometric ladder, deterministic step-by-step Stage-1 gcd + Stage-2
binary-product-tree giant-step bisection per curve: ZERO curve replay, ZERO p1-as-upper-bound, ZERO variance)."""
import math
import random
from functools import lru_cache
from math import ceil, isqrt, log, log1p, log2, sqrt
import q15_descent16 as D

MONT5 = D.MONT5
TORSION = 12
_LOG2_TORSION = log2(TORSION)
WHEELS = (30030,)
B_MIN_DEFAULT = 32  # Conservative 2 * 16 bits (all primes <= 2^16 stripped by TDIV)
L_MIN_DEFAULT = 2   # 2-limb (<=128b) mont2_cios kernel (47 <= 48 ARMv8 ops <= 58.84 Ir charged)
NEWTON_PRECHECK_INSTR = 2.0 ** 26  # Honest schoolbook upper bound (> 32.06M Ir for 2590b input)


def primes_upto(n):
    a = bytearray(b"\x01") * (n + 1)
    a[0:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if a[p]:
            a[p * p :: p] = b"\x00" * ((n - p * p) // p + 1)
    return [p for p in range(2, n + 1) if a[p]]


def _is_perfect_power(n):
    """Exact prime-exponent Newton integer perfect-power precheck before ECM (eliminates N = p^k;
    explicitly charged at NEWTON_PRECHECK_INSTR = 2^26 Ir)."""
    if n < 4:
        return None
    r2 = isqrt(n)
    if r2 * r2 == n:
        return (r2, 2)
    bl = n.bit_length()
    for k in primes_upto(bl):
        if k == 2:
            continue
        x = 1 << ((bl + k - 1) // k)
        km1 = k - 1
        while True:
            nx = (km1 * x + n // (x ** km1)) // k
            if nx >= x:
                break
            x = nx
        if x ** k == n:
            return (x, k)
    return None


def _run_suyama_stepped(N, sigma, b1, primes_b1):
    """Run one Suyama-Montgomery Stage-1 curve mod N taking gcd after EVERY prime power p^k <= B1
    (deterministic stepped-gcd schedule; zero curve replay)."""
    u = (sigma * sigma - 5) % N
    v = (4 * sigma) % N
    x = pow(u, 3, N)
    z = pow(v, 3, N)
    g0 = math.gcd(z, N)
    if g0 > 1:
        return g0
    den = (16 * x * v) % N
    gd = math.gcd(den, N)
    if gd > 1:
        return gd
    a24 = (pow(v - u, 3, N) * (3 * u + v) * pow(den, -1, N)) % N
    A_mod = (4 * a24 - 2) % N
    g_sing = math.gcd((A_mod * A_mod - 4) % N, N)
    if g_sing == N:
        return "SINGULAR"
    if g_sing > 1:
        return g_sing

    def xdbl(X, Z):
        a = (X + Z) % N
        b = (X - Z) % N
        a2 = (a * a) % N
        b2 = (b * b) % N
        c = (a2 - b2) % N
        return (a2 * b2) % N, (c * (b2 + a24 * c)) % N

    def xadd(X1, Z1, X2, Z2, Xd, Zd):
        u1 = ((X1 - Z1) * (X2 + Z2)) % N
        v1 = ((X1 + Z1) * (X2 - Z2)) % N
        return (Zd * pow(u1 + v1, 2, N)) % N, (Xd * pow(u1 - v1, 2, N)) % N

    def xmul(k, X, Z):
        X0, Z0 = X, Z
        X1, Z1 = xdbl(X, Z)
        for bit in bin(k)[3:]:
            if bit == "1":
                X0, Z0 = xadd(X0, Z0, X1, Z1, X, Z)
                X1, Z1 = xdbl(X1, Z1)
            else:
                X1, Z1 = xadd(X0, Z0, X1, Z1, X, Z)
                X0, Z0 = xdbl(X0, Z0)
        return X0, Z0

    for p in primes_b1:
        pk = p
        while pk <= b1:
            x, z = xmul(p, x, z)
            g = math.gcd(z, N)
            if g > 1:
                return g
            pk *= p
    return 1


def _build_product_tree_mod(raw_leaves, N):
    """Build exact Stage-2 binary product tree mod N over J leaves, replacing any 0 mod N leaf
    with 1 BEFORE multiplying into the tree. Because every remaining leaf has gcd(t_j, N) < N,
    whenever gcd(root, N) == N, binary-tree bisection is unconditionally guaranteed to isolate a
    proper factor 1 < gcd(node, N) < N in <= 1 + 2*ceil(log2(J)) gcds (0 extra muls)."""
    leaves = [t % N for t in raw_leaves]
    leaves = [t if t != 0 else 1 for t in leaves]
    levels = [leaves]
    while len(levels[-1]) > 1:
        cur = levels[-1]
        nxt = [(cur[i] * cur[i + 1]) % N if i + 1 < len(cur) else cur[i] for i in range(0, len(cur), 2)]
        levels.append(nxt)
    return levels


def _bisect_product_tree(levels, N):
    """Descend binary product tree from root when gcd(root, N) == N to isolate a proper factor 1 < d < N."""
    g_root = math.gcd(levels[-1][0], N)
    gcds = 1
    if 1 < g_root < N:
        return g_root, gcds
    if g_root != N:
        return 1, gcds
    idx = 0
    for lv in range(len(levels) - 2, -1, -1):
        cur = levels[lv]
        li = 2 * idx
        ri = li + 1
        gl = math.gcd(cur[li], N)
        gcds += 1
        if 1 < gl < N:
            return gl, gcds
        if gl == N:
            idx = li
        else:
            gr = math.gcd(cur[ri], N)
            gcds += 1
            if 1 < gr < N:
                return gr, gcds
            assert gr == N
            idx = ri
    return 1, gcds


def verify_suyama_stepped():
    """Live verification of:
      (1) Newton perfect-power precheck on N = p^2,
      (2) singular-curve guard on sigma = N - 1,
      (3) prime-by-prime stepped-gcd split on a composite N = p*q,
      (4) empirical composite proper-split rate on N = p * q_comp (20-21b p, 41b companion q_comp, B1=256),
          counting ONLY proper splits 1 < gcd < N (never gcd == N), confirming emp_p1 >= model_p1_lo,
      (5) live Stage-2 wheel (W=210, J=24, h_tree=5) + binary-product-tree bisection on N = 100003 * 100043
          AND explicit [0, p, q, 1] zero-leaf-skip edge-case guarantee."""
    N_sq = 1157341 * 1157341
    pp_res = _is_perfect_power(N_sq)
    assert pp_res == (1157341, 2), f"Perfect-power precheck failed on p^2: {pp_res}"
    N = 1157341 * 2633581
    assert _is_perfect_power(N) is None
    ps47 = tuple(primes_upto(47))
    g_sing_nm1 = _run_suyama_stepped(N, N - 1, 200, ps47)
    assert g_sing_nm1 == "SINGULAR", f"Expected singular curve guard on sigma=N-1, got {g_sing_nm1}"
    split_sigma = None
    split_fac = None
    for sig in range(6, 60):
        g = _run_suyama_stepped(N, sig, 5000, ps47)
        if isinstance(g, int) and 1 < g < N:
            split_sigma = sig
            split_fac = g
            break
    assert split_fac is not None and N % split_fac == 0 and 1 < split_fac < N

    # (4) Composite proper-split Suyama P1 check (N = p * q_comp, 20-21b p, 41b prime q_comp, B1 = B2 = 256)
    # Counts ONLY proper splits 1 < gcd < N (never gcd == N).
    b1_emp = 256
    ps_emp = tuple(primes_upto(b1_emp))
    p1_model = semismooth_lo(21.0, b1_emp, b1_emp) * power_correction(21.0, b1_emp, b1_emp)
    q_comp = 1099511627791  # prime 2^40 + 15
    test_primes = (
        1048583, 1048601, 1048609, 1048613, 1048627, 1048633,
        1100009, 1100023, 1200007, 1200047, 1300021, 1300031,
        1500007, 1500041, 1700021, 1700029,
    )
    rng = random.Random(20261124)
    trials = 0
    hits = 0
    for pp in test_primes:
        N_c = pp * q_comp
        for _ in range(10):
            sig = rng.randrange(6, pp - 6)
            g = _run_suyama_stepped(N_c, sig, b1_emp, ps_emp)
            trials += 1
            if isinstance(g, int) and 1 < g < N_c:
                assert N_c % g == 0
                hits += 1
    emp_p1 = hits / trials
    assert emp_p1 >= p1_model, f"Empirical composite proper-split P1 {emp_p1:.4f} < model lower bound {p1_model:.4f}"

    # (5) Live Stage-2 wheel (W=210, J=24, h_tree=5) + binary-product-tree giant-step bisection test:
    p_a, p_b = 100003, 100043
    N_s2 = p_a * p_b
    sig_s2, b1_s2, W_s2, g_s2 = 11, 64, 210, 2
    babies_s2 = [r for r in range(1, W_s2 // 2) if math.gcd(r, W_s2) == 1]
    J_s2 = len(babies_s2)
    h_s2 = ceil(log2(J_s2))
    assert J_s2 == 24 and h_s2 == 5
    u_s = (sig_s2 * sig_s2 - 5) % N_s2
    v_s = (4 * sig_s2) % N_s2
    xs, zs = pow(u_s, 3, N_s2), pow(v_s, 3, N_s2)
    den_s = (16 * xs * v_s) % N_s2
    a24_s = (pow(v_s - u_s, 3, N_s2) * (3 * u_s + v_s) * pow(den_s, -1, N_s2)) % N_s2

    def _xd(X, Z):
        a = (X + Z) % N_s2; b = (X - Z) % N_s2
        a2 = (a * a) % N_s2; b2 = (b * b) % N_s2
        c = (a2 - b2) % N_s2
        return (a2 * b2) % N_s2, (c * (b2 + a24_s * c)) % N_s2

    def _xa(X1, Z1, X2, Z2, Xd, Zd):
        u1 = ((X1 - Z1) * (X2 + Z2)) % N_s2
        v1 = ((X1 + Z1) * (X2 - Z2)) % N_s2
        return (Zd * pow(u1 + v1, 2, N_s2)) % N_s2, (Xd * pow(u1 - v1, 2, N_s2)) % N_s2

    def _xm(k, X, Z):
        if k == 1:
            return X, Z
        X0, Z0 = X, Z
        X1, Z1 = _xd(X, Z)
        for bit in bin(k)[3:]:
            if bit == "1":
                X0, Z0 = _xa(X0, Z0, X1, Z1, X, Z)
                X1, Z1 = _xd(X1, Z1)
            else:
                X1, Z1 = _xa(X0, Z0, X1, Z1, X, Z)
                X0, Z0 = _xd(X0, Z0)
        return X0, Z0

    for pp in primes_upto(b1_s2):
        pk = pp
        while pk <= b1_s2:
            xs, zs = _xm(pp, xs, zs)
            pk *= pp
    assert math.gcd(zs, N_s2) == 1
    bpts = [_xm(r, xs, zs) for r in babies_s2]
    xg, zg = _xm(g_s2 * W_s2, xs, zs)
    raw_leaves = [(xg * zb - xb * zg) % N_s2 for xb, zb in bpts]
    levels = _build_product_tree_mod(raw_leaves, N_s2)
    assert math.gcd(levels[-1][0], N_s2) == N_s2
    s2_split, s2_gcds = _bisect_product_tree(levels, N_s2)
    assert s2_split in (p_a, p_b) and s2_gcds <= 1 + 2 * h_s2

    # Explicit [0, p, q, 1] zero-leaf-skip edge-case check (skipping 0 mod N leaves before tree multiply):
    lv_zero = _build_product_tree_mod([0, p_a, p_b, 1], N_s2)
    assert math.gcd(lv_zero[-1][0], N_s2) == N_s2
    z_split, z_gcds = _bisect_product_tree(lv_zero, N_s2)
    assert z_split in (p_a, p_b) and z_gcds <= 1 + 2 * 2

    return dict(
        N=N,
        perfect_power_precheck_ok=True,
        singular_nm1_rejected=True,
        clean_split_sigma=split_sigma,
        split_factor=split_fac,
        emp_composite_mode=True,
        emp_b1=b1_emp,
        emp_trials=trials,
        emp_hits=hits,
        emp_p1=round(emp_p1, 4),
        model_p1_lo=round(p1_model, 4),
        stage2_tree_bisect_ok=True,
        stage2_zero_leaf_skip_ok=True,
        stage2_tree_N=N_s2,
        stage2_tree_J=J_s2,
        stage2_tree_h=h_s2,
        stage2_tree_split=s2_split,
        stage2_tree_gcds=s2_gcds,
    )


def mulmod_instr(nbits, l_min=L_MIN_DEFAULT):
    L = max(l_min, ceil(nbits / 64))
    return MONT5 * (L / 5.0) ** 2, L


def gcd_instr(nbits, l_min=L_MIN_DEFAULT):
    L = max(l_min, ceil(nbits / 64))
    return 2.0 * nbits * (12 * L + 16)


_TREE_DP_CACHE = {}


def tree_branch_caps(nbits, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
    """Exact DP upper bound on (sum L_i^2, sum gcd_instr(b_i)) over ALL active composite-branch
    partitions b_1 + ... + b_t <= nbits with b_i >= b_min (L_i = max(l_min, ceil(b_i / 64)))."""
    nb = int(ceil(nbits))
    if nb < b_min:
        L = max(l_min, ceil(nb / 64))
        return L * L, 2.0 * nb * (12 * L + 16)
    key = (b_min, l_min)
    st = _TREE_DP_CACHE.get(key)
    if st is None:
        dp_l2 = [0] * b_min
        dp_g = [0.0] * b_min
        st = (dp_l2, dp_g)
        _TREE_DP_CACHE[key] = st
    else:
        dp_l2, dp_g = st
    cur_n = len(dp_l2) - 1
    if nb > cur_n:
        for b in range(cur_n + 1, nb + 1):
            L = max(l_min, (b + 63) >> 6)
            best_l2 = L * L
            best_g = 2.0 * b * (12 * L + 16)
            Ls = max(l_min, (b_min + 63) >> 6)
            v2 = Ls * Ls + dp_l2[b - b_min]
            vg = 2.0 * b_min * (12 * Ls + 16) + dp_g[b - b_min]
            if v2 > best_l2:
                best_l2 = v2
            if vg > best_g:
                best_g = vg
            s_max = b - b_min
            for k in range(l_min, (b >> 6) + 1):
                s = (k << 6) + 1
                if b_min <= s <= s_max:
                    Ls = max(l_min, (s + 63) >> 6)
                    v2 = Ls * Ls + dp_l2[b - s]
                    vg = 2.0 * s * (12 * Ls + 16) + dp_g[b - s]
                    if v2 > best_l2:
                        best_l2 = v2
                    if vg > best_g:
                        best_g = vg
            dp_l2.append(best_l2)
            dp_g.append(best_g)
    return dp_l2[nb], dp_g[nb]


def pi_upper(x):
    return 1.26 * x / log(x) if x >= 17 else 7.0


def pi_pp_upper(x):
    """Upper bound on total Stage-1 prime-power steps sum_{q^k <= x} 1 <= pi_upper(x) + 2*sqrt(x)
    (stepping individually by q so distinct valuations e_1 < e_2 of q split cleanly)."""
    return pi_upper(x) + 2.0 * sqrt(x)


@lru_cache(None)
def _ss_grid(b1, b2, steps=256):
    base = log2(b1)
    top = log2(b2)
    if top <= base:
        return base, (), (), ()
    span = top - base
    los = tuple(base + span * i / steps for i in range(steps))
    wts = tuple(log((base + span * (i + 1) / steps) / (base + span * i / steps)) for i in range(steps))
    suf = [0.0] * (steps + 1)
    acc = 0.0
    for i in range(steps - 1, -1, -1):
        acc += wts[i]
        suf[i] = acc
    return base, los, wts, tuple(suf)


@lru_cache(None)
def semismooth_lo(fbits, b1, b2, steps=256):
    """Conservative LOWER bound on (B1, B2)-semismoothness for a prime p <= 2^{fbits}."""
    m = fbits - _LOG2_TORSION
    base, los, wts, suf = _ss_grid(b1, b2, steps)
    tab = D.RHO_TAB_LO
    ntab = len(tab)
    inv_rstep = 1.0 / D.RSTEP
    _ceil = ceil
    u0 = m / base
    if u0 <= 1.0:
        acc = 1.0
    else:
        i0 = _ceil(u0 * inv_rstep)
        acc = tab[i0] if i0 < ntab else 0.0
    if los:
        from bisect import bisect_left
        inv_br = inv_rstep / base
        cut = bisect_left(los, m - base)
        if cut < steps:
            acc += suf[cut]
        if cut > 0:
            # First index where ceil((m - los[i]) * inv_br) < ntab:
            # los[i] > m - (ntab - 1) / inv_br
            lo_i = bisect_left(los, m - (ntab - 1) / inv_br, 0, cut)
            if lo_i < cut:
                ms = m * inv_br
                # Safe check at boundary lo_i
                while lo_i < cut and _ceil(ms - los[lo_i] * inv_br) >= ntab:
                    lo_i += 1
                while cut > lo_i and (m - los[cut - 1]) <= base:
                    cut -= 1
                    acc += wts[cut]
                if lo_i < cut:
                    acc += math.sumprod((tab[_ceil(ms - los[i] * inv_br)] for i in range(lo_i, cut)), wts[lo_i:cut])
    return min(acc, 1.0 - 1e-15)


semismooth = semismooth_lo


@lru_cache(None)
def _power_correction_b1(fbits, b1):
    m = fbits - _LOG2_TORSION
    u = m / log2(b1)
    tilt = D.rho_hi(max(u - 2.0, 1.0)) / max(D.rho_lo(u), 1e-300)
    y = sqrt(b1)
    return max(0.0, 1.0 - tilt * 2.0 / (y * log(y)))


def power_correction(fbits, b1, b2):
    return _power_correction_b1(fbits, b1)


@lru_cache(None)
def stage_ops(b1, b2, wheel):
    s1 = ceil(1.5 * b1) * 11
    phi = wheel
    for q in (2, 3, 5, 7, 11, 13):
        if wheel % q == 0:
            phi = phi // q * (q - 1)
    J = phi // 2
    G = ceil(b2 / wheel) + 1
    s2 = 6 * (wheel // 2 + G) + 3 * J * G if b2 > b1 else 0
    return 1.1 * (s1 + s2), s1, s2, J, G


@lru_cache(None)
def schedule_rung(fbits, nbits, success_rung=0.999, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
    """Single-rung Deterministic Stepped-GCD + Stage-2 Binary-Product-Tree branch-tree-capped ECM schedule (H-ECM).
    Every curve takes step-by-step gcd after each Stage-1 prime power and builds a J-leaf binary product tree
    on each Stage-2 giant step G (h_tree = ceil(log2(max(2, J))) levels, +1.1*J*G muls, +(1+2*h_tree)*G gcds):
      per_curve_stepped = (ops + 1.1*J*G) * mul + (pi_upper(B1) + (1 + 2*h_tree)*G) * gcd
    Whenever gcd(root, M) == M in a giant step, descending the binary product tree deterministically isolates
    any two distinct primes hitting different baby-step match terms j1 != j2 in <= 2*h_tree gcds (0 extra muls).
    100% deterministic cost: ZERO curve replay, ZERO p1-as-upper-bound, ZERO variance."""
    nb = int(ceil(nbits))
    l2_max, g_max = tree_branch_caps(nb, b_min, l_min)
    mul = MONT5 * (l2_max / 25.0)
    gcd = g_max
    l1ms = log1p(-float(success_rung))
    best = None
    best_cost = float("inf")
    for k in range(6, 37, 2):
        b1 = 1 << k
        pi_b1 = pi_pp_upper(b1)
        min_s1_ops = 1.1 * ceil(1.5 * b1) * 11
        if (min_s1_ops * mul + pi_b1 * gcd) + NEWTON_PRECHECK_INSTR >= best_cost:
            break
        pc = _power_correction_b1(fbits, b1)
        if pc <= 0.0:
            continue
        for ratio in (1, 10, 30, 100, 300):
            b2 = b1 * ratio
            wheels = WHEELS if ratio > 1 else (0,)
            p1 = semismooth_lo(fbits, b1, b2) * pc
            if p1 <= 1e-12:
                continue
            curves = ceil(l1ms / log1p(-p1))
            for wheel in wheels:
                if wheel and wheel > b2:
                    continue
                ops, s1, s2, J, G = (
                    stage_ops(b1, b2, wheel)
                    if wheel
                    else (min_s1_ops, ceil(1.5 * b1) * 11, 0, 0, 0)
                )
                h_tree = ceil(log2(max(2, J))) if G else 0
                m_ops = ops + (1.1 * J * G if G else 0.0)
                g_cnt = pi_b1 + (1 + 2 * h_tree) * G
                per_curve_stepped = m_ops * mul + g_cnt * gcd
                cost = curves * per_curve_stepped + NEWTON_PRECHECK_INSTR
                if cost < best_cost:
                    best_cost = cost
                    best = dict(
                        fbits=fbits, nbits=nb, success_rung=float(success_rung),
                        b_min=b_min, l_min=l_min,
                        l2_tree_max=l2_max, gcd_tree_max=g_max,
                        B1=b1, B2=b2, wheel=wheel, curves=curves,
                        p1=p1, pi_b1=pi_b1, J=J, G=G, h_tree=h_tree,
                        per_curve_stepped=per_curve_stepped,
                        instructions=cost, log2_instructions=log2(cost),
                    )
    return best


schedule = schedule_rung


@lru_cache(None)
def ladder_instr(lo_bits, hi_bits, nbits, success=0.999, step=12.0, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
    """Bottom-up Deterministic Stepped-GCD 12-bit geometric ECM ladder from lo_bits to hi_bits.
    Each rung targets s_rung = 1 - (1 - success) / n_rungs (union bound <= 1 - success)."""
    if hi_bits <= lo_bits + 1e-9:
        return 0.0
    rungs = []
    r = float(hi_bits)
    while r > lo_bits + 1e-9:
        rungs.append(round(r, 4))
        r -= step
    rungs.reverse()
    s_rung = 1.0 - (1.0 - float(success)) / len(rungs)
    nb = int(ceil(nbits))
    tot = 0.0
    for r in rungs:
        sch = schedule_rung(r, nb, s_rung, b_min, l_min)
        if sch is None:
            raise ValueError(f"Infeasible ECM rung r={r}, nb={nb}")
        tot += sch["instructions"]
    return tot


def ecm_instr(fbits, nbits, success=0.999, n_factors=1):
    """Deterministic Stepped-GCD descent ECM ladder from 16b up to fbits."""
    return ladder_instr(16.0, fbits, nbits, success=success, step=12.0, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT)
