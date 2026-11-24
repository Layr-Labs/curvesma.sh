"""q15_ecm16.py — Deterministic Stepped-GCD + Stage-2 Binary-Product-Tree Bisection ECM schedule (STEPPED_LADDER).
Every curve deterministically pays:
  - Stage 1: stepped gcd after every prime power q^e <= B1 (pi_upper(B1) gcds);
  - Stage 2 (wheel W = 30030, J = phi(W)//2 = 2880, G = ceil(B2/W) + 1):
      for each giant step g = 1..G, multiply the J match terms via a binary product tree of height
      h_tree = ceil(log2(J)) (J - 1 muls + explicit +1.1*J*G mul cushion -> m_ops = ops + 1.1*J*G),
      take d = gcd(root, M), and whenever d == M descend the binary product tree (<= h_tree gcds)
      to isolate distinct prime factors dividing distinct baby-step terms j1 != j2 with ZERO replay;
      we deterministically charge g_steps = pi_upper(B1) + (1 + 2 * h_tree) * G gcds on EVERY curve
      over exact L_min = 2 branch-tree DP caps (tree_branch_caps), with 12-bit geometric ladder
      union-bounded across rungs under standard literature H-ECM (Brent / Zimmermann-Dodson)."""
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


def _mont_setup(N, sigma):
    """Initialize Suyama-Montgomery curve parameters and x-only ladder ops mod N."""
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
        if k == 1:
            return X, Z
        if k == 2:
            return xdbl(X, Z)
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

    return x, z, xmul


def _run_suyama_stepped(N, sigma, b1, primes_b1):
    """Run one Suyama-Montgomery Stage-1 curve mod N with step-by-step gcd after each prime power."""
    st = _mont_setup(N, sigma)
    if not isinstance(st, tuple):
        return st
    x, z, xmul = st
    for p in primes_b1:
        pk = p
        while pk * p <= b1:
            pk *= p
        x, z = xmul(pk, x, z)
        g = math.gcd(z, N)
        if g > 1:
            return g
    return 1


def _build_prod_tree(terms, N):
    """Build a balanced binary product tree over J Stage-2 match terms mod N (uses J - 1 muls).
    Returns (node_val, left_subtree, right_subtree) where leaves have left = right = None."""
    nodes = [(t % N, None, None) for t in terms]
    while len(nodes) > 1:
        nxt = []
        i = 0
        n = len(nodes)
        while i + 1 < n:
            L, R = nodes[i], nodes[i + 1]
            nxt.append(((L[0] * R[0]) % N, L, R))
            i += 2
        if i < n:
            nxt.append(nodes[i])
        nodes = nxt
    return nodes[0]


def _stage2_tree_bisect(tree, N):
    """Given a binary product tree with gcd(tree[0], N) == N, descend <= ceil(log2(J)) levels
    using 1 gcd(left, N) per level to isolate a non-trivial factor 1 < d < N with ZERO extra muls."""
    node = tree
    gcds = 0
    while node[1] is not None:
        left, right = node[1], node[2]
        d = math.gcd(left[0], N)
        gcds += 1
        if 1 < d < N:
            return d, gcds
        node = left if d == N else right
    return N, gcds


def _run_suyama_stage1_stage2_bisect(N, sigma, b1, b2, wheel, primes_b1):
    """Execute full Stage-1 (stepped gcd) + Stage-2 wheel with binary-product-tree giant-step
    bisection on composite N. Returns (stage, giant_step, factor, root_gcd, bisect_gcds, J, h_tree)."""
    st = _mont_setup(N, sigma)
    if not isinstance(st, tuple):
        return ("INIT", 0, st, 0, 0, 0, 0)
    x, z, xmul = st
    for p in primes_b1:
        pk = p
        while pk * p <= b1:
            pk *= p
        x, z = xmul(pk, x, z)
        g = math.gcd(z, N)
        if g > 1:
            return ("S1", 0, g, g, 0, 0, 0)
    baby_js = [j for j in range(1, wheel // 2) if math.gcd(j, wheel) == 1]
    J = len(baby_js)
    h_tree = ceil(log2(max(2, J)))
    bpts = [xmul(j, x, z) for j in baby_js]
    G_max = ceil(b2 / wheel) + 1
    for g_step in range(1, G_max + 1):
        Xg, Zg = xmul(g_step * wheel, x, z)
        gz = math.gcd(Zg, N)
        if 1 < gz < N:
            return ("S2_G", g_step, gz, gz, 0, J, h_tree)
        terms = [(Xg * Zj - Zg * Xj) % N for (Xj, Zj) in bpts]
        tree = _build_prod_tree(terms, N)
        d_root = math.gcd(tree[0], N)
        if 1 < d_root < N:
            return ("S2_ROOT", g_step, d_root, d_root, 0, J, h_tree)
        if d_root == N:
            fac, b_gcds = _stage2_tree_bisect(tree, N)
            return ("S2_BISECT", g_step, fac, d_root, b_gcds, J, h_tree)
    return ("MISS", G_max, 1, 1, 0, J, h_tree)


def verify_suyama_ecm_live():
    """Live verification of:
      (1) Newton perfect-power precheck on N = p^2,
      (2) singular-curve guard on sigma = N - 1,
      (3) clean Stage-1 stepped-gcd factor splitting on composite N,
      (4) empirical Stage-1 factor-finding rate on 20-bit primes (B1 = 256) vs semismooth lower bound,
      (5) REAL Stage-2 wheel (W = 210, J = 24, h_tree = 5) + binary-product-tree bisection on a
          composite N_s2 = 1000039 * 1000081 (sigma = 11, B1 = 30, B2 = 1200) whose two prime factors
          collide on the SAME Stage-2 giant step g = 1 at distinct baby steps (j = 61 vs j = 97),
          confirming gcd(root, N_s2) == N_s2 and deterministic bisection split in <= ceil(log2(J)) gcds."""
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

    # Empirical Stage-1 factor-finding rate on 20-bit primes p (B1 = 256, B2 = 256) against semismooth lower bound:
    b1_emp = 256
    ps_emp = tuple(primes_upto(b1_emp))
    p_mod_lo = semismooth(20.0, b1_emp, b1_emp) * _power_correction_b1(20.0, b1_emp)
    q_big = 1000000007
    test_ps = (524309, 580007, 650011, 720007, 790003, 860009, 930017, 1000003)
    rng = random.Random(20261124)
    trials = 0
    found = 0
    for pp in test_ps:
        Nc = pp * q_big
        for _ in range(10):
            sig = rng.randrange(6, Nc - 6)
            g = _run_suyama_stepped(Nc, sig, b1_emp, ps_emp)
            trials += 1
            if isinstance(g, int) and g == pp:
                found += 1
    emp_rate = found / trials
    se = math.sqrt(p_mod_lo * (1.0 - p_mod_lo) / trials)
    assert emp_rate >= p_mod_lo - 3.0 * se and found >= 10, f"Empirical Stage-1 hit rate {emp_rate:.4f} too far below {p_mod_lo:.4f}"

    # Live Stage-2 wheel + binary-product-tree giant-step bisection on a real gcd(root, N_s2) == N_s2 collision:
    p_s2, q_s2, sig_s2 = 1000039, 1000081, 11
    N_s2 = p_s2 * q_s2
    ps30 = tuple(primes_upto(30))
    stage_tag, g_hit, fac_s2, d_root, b_gcds, J_s2, h_s2 = _run_suyama_stage1_stage2_bisect(
        N_s2, sig_s2, 30, 1200, 210, ps30
    )
    assert stage_tag == "S2_BISECT", f"Expected Stage-2 bisection split, got {stage_tag}"
    assert d_root == N_s2, f"Expected simultaneous giant-step collision gcd==N_s2, got {d_root}"
    assert fac_s2 in (p_s2, q_s2) and N_s2 % fac_s2 == 0
    assert 1 <= b_gcds <= h_s2, f"Bisection gcds {b_gcds} exceeded tree height {h_s2}"

    # Also exercise full multi-level binary-product-tree descent (both d==N left and d==1 right branches, depth == h_s2):
    deep_terms = [1] * J_s2
    deep_terms[J_s2 - 2] = p_s2
    deep_terms[J_s2 - 1] = q_s2
    deep_tree = _build_prod_tree(deep_terms, N_s2)
    assert math.gcd(deep_tree[0], N_s2) == N_s2
    deep_fac, deep_gcds = _stage2_tree_bisect(deep_tree, N_s2)
    assert deep_fac in (p_s2, q_s2) and 2 <= deep_gcds <= h_s2

    return dict(
        N=N, perfect_power_precheck_ok=True, singular_nm1_rejected=True,
        clean_split_sigma=split_sigma, split_factor=split_fac,
        emp_b1=b1_emp, emp_trials=trials, emp_found=found,
        emp_hit_rate=round(emp_rate, 4), model_p1_lo=round(p_mod_lo, 4),
        stage2_bisect_ok=True, stage2_N=N_s2, stage2_sigma=sig_s2,
        stage2_wheel=210, stage2_J=J_s2, stage2_h_tree=h_s2,
        stage2_giant_step=g_hit, stage2_root_gcd_eq_N=True,
        stage2_split_factor=fac_s2, stage2_bisect_gcds=b_gcds,
        stage2_deep_bisect_ok=True, stage2_deep_bisect_gcds=deep_gcds,
    )


def mulmod_instr(nbits, l_min=L_MIN_DEFAULT):
    L = max(l_min, ceil(nbits / 64))
    return MONT5 * (L / 5.0) ** 2, L


def gcd_instr(nbits, l_min=L_MIN_DEFAULT):
    L = max(l_min, ceil(nbits / 64))
    return 2.0 * nbits * (12 * L + 16)


def pi_upper(x):
    return 1.26 * x / log(x) if x >= 17 else 7.0


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
def semismooth(fbits, b1, b2, steps=256):
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
        inv_b = 1.0 / base
        for i in range(steps):
            u = (m - los[i]) * inv_b
            if u <= 1.0:
                acc += suf[i]
                break
            idx = _ceil(u * inv_rstep)
            if idx < ntab:
                acc += tab[idx] * wts[i]
    return min(acc, 1.0 - 1e-15)


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
def _rung_candidates_stepped(fbits):
    """Precompute all feasible (b1, b2, wheel, m_ops, gcd_steps, p1, log1m_p1) for stepped-gcd
    Stage-1 + binary-product-tree Stage-2 giant-step bisection ECM.
    When G > 0, each giant step builds a binary product tree of height h_tree = ceil(log2(max(2, J)))
    over its J match terms (explicit +1.1*J*G mul cushion -> m_ops = ops + 1.1*J*G) and pays
    (1 + 2*h_tree) gcds deterministically on every giant step (g_steps = pi_upper(B1) + (1 + 2*h_tree)*G)."""
    out = []
    for k in range(10, 37, 2):
        if fbits / k > 18.0:
            continue
        b1 = 1 << k
        pc = _power_correction_b1(fbits, b1)
        if pc <= 0.0:
            continue
        min_s1_ops = 1.1 * ceil(1.5 * b1) * 11
        pi_b1 = pi_upper(b1)
        for ratio in (1, 10, 30, 100, 300):
            b2 = b1 * ratio
            if ratio > 1 and b2 < 30030:
                continue
            p1 = semismooth(fbits, b1, b2) * pc
            if p1 <= 1e-12:
                continue
            l1m = log1p(-p1)
            if ratio > 1:
                wheel = 30030
                ops, _, _, J, G = stage_ops(b1, b2, wheel)
                h_tree = ceil(log2(max(2, J)))
                m_ops = ops + 1.1 * J * G
                g_steps = pi_b1 + (1 + 2 * h_tree) * G
            else:
                wheel = 0
                m_ops = min_s1_ops
                g_steps = pi_b1
            out.append((b1, b2, wheel, m_ops, g_steps, p1, l1m))
    return tuple(out)


@lru_cache(None)
def schedule_rung(fbits, nbits, success=0.999, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
    """Single-rung deterministic stepped-gcd + Stage-2 binary-product-tree bisection ECM schedule (H-ECM).
    Every curve pays deterministically:
      per_curve = (ops + 1.1*J*G) * mul + (pi_upper(B1) + (1 + 2*ceil(log2(J))) * G) * gcd.
    Curve count follows standard literature H-ECM: curves = ceil(log1p(-success) / log1p(-p1))."""
    nb = int(ceil(nbits))
    l2_max, g_max = tree_branch_caps(nb, b_min, l_min)
    mul = MONT5 * (l2_max / 25.0)
    gcd = g_max
    l1ms = log1p(-success)
    best = None
    best_cost = float("inf")
    for b1, b2, wheel, m_ops, g_steps, p1, l1m in _rung_candidates_stepped(fbits):
        per_curve = m_ops * mul + g_steps * gcd
        if per_curve + NEWTON_PRECHECK_INSTR >= best_cost:
            continue
        curves = ceil(l1ms / l1m)
        cost = curves * per_curve + NEWTON_PRECHECK_INSTR
        if cost < best_cost:
            best_cost = cost
            best = dict(
                fbits=fbits, nbits=nb, success=success, b_min=b_min, l_min=l_min,
                l2_tree_max=l2_max, gcd_tree_max=g_max,
                B1=b1, B2=b2, wheel=wheel, curves=curves, p_curve=p1,
                mul_ops=m_ops, gcd_steps=g_steps, per_curve_stepped=per_curve,
                instructions=cost, log2_instructions=log2(cost),
            )
    return best


schedule = schedule_rung


@lru_cache(None)
def ladder_instr(lo_bits, hi_bits, nbits, success=0.999, step=12.0, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
    """Bottom-up deterministic stepped-gcd + Stage-2 bisection 12-bit geometric ECM ladder.
    Union bound across n_rungs rungs gives s_rung = 1 - (1 - success) / n_rungs -> total failure <= 1 - success."""
    if hi_bits <= lo_bits + 1e-9:
        return 0.0
    rungs = []
    r = float(hi_bits)
    while r > lo_bits + 1e-9:
        rungs.append(round(r, 4))
        r -= step
    rungs.reverse()
    s_rung = 1.0 - (1.0 - success) / len(rungs)
    nb = int(ceil(nbits))
    tot = 0.0
    for r in rungs:
        sch = schedule_rung(r, nb, s_rung, b_min, l_min)
        if sch is None:
            raise ValueError(f"Infeasible stepped ECM rung r={r}, nb={nb}")
        tot += sch["instructions"]
    return tot


def ecm_instr(fbits, nbits, success=0.999, n_factors=1):
    """Deterministic stepped-gcd + Stage-2 bisection descent ECM ladder from 16b up to fbits."""
    return ladder_instr(16.0, fbits, int(ceil(nbits)), success=success, step=12.0, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT)
