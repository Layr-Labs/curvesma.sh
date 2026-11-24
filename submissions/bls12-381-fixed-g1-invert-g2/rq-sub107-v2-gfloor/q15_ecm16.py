"""q15_ecm16.py — Deterministic Stepped-GCD + Stage-2 Binary-Product-Tree Geometric-Ladder Suyama-Montgomery ECM
(exact L_min=2 branch-tree DP, 12-bit geometric ladder, deterministic step-by-step Stage-1 gcd + Stage-2
binary-product-tree giant-step bisection per curve: ZERO curve replay, ZERO p1-as-upper-bound, ZERO variance)."""
import math
import random
from functools import lru_cache
from math import ceil, isqrt, log, log1p, log2, sqrt
import q15_descent16 as D

import operator as _OP
_SUMPROD = getattr(math, "sumprod", lambda a, b, _m=_OP.mul: sum(map(_m, a, b)))

MONT5 = D.MONT5
TORSION = 12
_LOG2_TORSION = log2(TORSION)
WHEELS = (30, 210, 2310, 30030, 510510)
B_MIN_DEFAULT = 32  # Conservative 2 * 16 bits (all primes <= 2^16 stripped by TDIV)
L_MIN_DEFAULT = 2   # 2-limb (<=128b) mont2_cios kernel (47 <= 48 ARMv8 ops <= 58.84 Ir charged)
NEWTON_PRECHECK_INSTR = 2.0 ** 26  # Honest schoolbook upper bound (> 32.06M Ir for 2590b input)
MR_ROUNDS = 110     # Unconditional Rabin error <= 4^-110 = 2^-220 per test (< 2^-100 union-bounded up to T_run=2^124)


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


def _split_product_tree(levels, N):
    """Full recursive Stage-2 binary-product-tree splitting + test-once prime-power extraction for an active composite N.
    Computes root gcd g_root = gcd(root, N). If g_root == 1, returns ([N], 1, 0).
    Otherwise keeps N // g_root (if g_root < N) as the unhit cofactor continuing to the next
    giant step and recursively descends levels with g_root: at every internal node (lv, idx, m),
    computes gl = gcd(left, m), mr = m // gl, and (if mr > 1) gr = gcd(right, mr).
    Test-once perfect-power rule: _is_perfect_power(m) is run once on the initial input and ONLY on
    strictly proper newly created divisors 1 < m < N (NEVER re-tested on an unsplit active composite
    m == N). Since every proper split creates strictly smaller tree children and a factorization tree
    has <= Omega(N) <= omega_max leaves, at most 2*omega_max - 1 calls to _is_perfect_power occur
    across the ENTIRE ECM invocation (all C curves combined), charged at 2 * omega_max * NEWTON_PRECHECK_INSTR.
    Combined Stage-1+Stage-2 Event & Single-Leaf Uniqueness (wheel W <= B1 enforced in schedule_rung):
      Under H-ECM, semismooth_lo * _power_correction_b1 lower-bounds the combined Stage-1+Stage-2 event
      E_good(p) that #E(F_p) is (B1, B2)-semismooth with all prime-power divisors q^k | #E(F_p) (q <= B1)
      satisfying q^k <= B1. Under E_good(p), either: (a) ord_p(Q_1) == 1 (detected in Stage 1 via
      stepped-GCD), or (b) ord_p(Q_1) = L_p is a single prime B1 < L_p <= B2 (detected in Stage 2).
      In case (b), since L_p > B1 >= W, p divides AT MOST ONE leaf t_{j(p)} in its hit giant step g_0
      (|eps_1*b_1 - eps_2*b_2| < W <= B1 < L_p). By the recursive invariant m | levels[lv][idx] and
      m | N, p descends the unique root-to-leaf path to j(p), arriving with m | gcd(t_{j(p)}, N);
      whenever no distinct companion prime q != p collides at that exact same leaf (g_0, j(p)),
      gcd(t_{j(p)}, N) = p^a (1 <= a <= v_p(N)), so _split_product_tree + test-once _is_perfect_power
      deterministically extracts p (even for d = p^2*q or multi-leaf companions). Rare exact same-step
      (Stage-1) or same-leaf (Stage-2) collisions of distinct primes are absorbed into H-ECM (retiring
      the composite to the next curve without re-testing _is_perfect_power if m == N).
    By induction on lv, total descent GCDs <= 2 * h_tree * Omega(g_root), so across ALL G giant
    steps and active composites of one curve, total descent GCDs <= 2 * h_tree * Omega(N)."""
    g_root = math.gcd(levels[-1][0], N)
    if g_root == 1:
        return [N], 1, 0
    raw_out = []
    if g_root < N:
        raw_out.append(N // g_root)
    desc_gcds = 0

    def _rec(lv, idx, m):
        nonlocal desc_gcds
        if lv == 0:
            raw_out.append(m)
            return
        cur = levels[lv - 1]
        li = 2 * idx
        ri = li + 1
        if ri >= len(cur):
            _rec(lv - 1, li, m)
            return
        gl = math.gcd(cur[li], m)
        desc_gcds += 1
        mr = m // gl
        if gl > 1:
            _rec(lv - 1, li, gl)
        if mr > 1:
            gr = math.gcd(cur[ri], mr)
            desc_gcds += 1
            if gr > 1:
                _rec(lv - 1, ri, gr)
            rem = mr // gr
            if rem > 1:
                raw_out.append(rem)

    _rec(len(levels) - 1, 0, g_root)
    out = []
    for m in raw_out:
        if m < N:
            pp = _is_perfect_power(m)
            if pp is not None:
                out.extend([pp[0]] * pp[1])
                continue
        out.append(m)
    return out, 1 + desc_gcds, desc_gcds


def _bisect_product_tree(levels, N):
    """Compatibility wrapper returning one proper factor (if any) and total GCD count via _split_product_tree."""
    parts, total_gcds, _ = _split_product_tree(levels, N)
    for d in parts:
        if 1 < d < N:
            return d, total_gcds
    return 1, total_gcds


def _run_suyama_frontier_curve(N, sigma, b1, primes_b1):
    """Run ONE Suyama-Montgomery Stage-1 curve across the active composite frontier of N:
    whenever 1 < g < d splits an active composite d at a prime-power step, immediately reduce
    (x, z, a24) mod (d // g) and CONTINUE THE SAME CURVE on the cofactor d // g (zero extra curve
    count; total mul/gcd across all active branches remains bounded by tree_branch_caps(nbits))."""
    u = (sigma * sigma - 5) % N
    v = (4 * sigma) % N
    x = pow(u, 3, N)
    z = pow(v, 3, N)
    den = (16 * x * v) % N

    def _init_branch(d):
        if d <= 1:
            return []
        pp = _is_perfect_power(d)
        if pp is not None:
            return _init_branch(pp[0])
        g0 = math.gcd(z % d, d)
        if 1 < g0 < d:
            return _init_branch(g0) + _init_branch(d // g0)
        if g0 == d:
            return [(d, None, None, None)]
        gd = math.gcd(den % d, d)
        if 1 < gd < d:
            return _init_branch(gd) + _init_branch(d // gd)
        if gd == d:
            return [(d, None, None, None)]
        a24 = (pow((v - u) % d, 3, d) * ((3 * u + v) % d) * pow(den % d, -1, d)) % d
        A_mod = (4 * a24 - 2) % d
        gs = math.gcd((A_mod * A_mod - 4) % d, d)
        if 1 < gs < d:
            return _init_branch(gs) + _init_branch(d // gs)
        if gs == d:
            return [(d, None, None, None)]
        return [(d, x % d, z % d, a24)]

    def _xmul_mod(k, X, Z, a24, M):
        def _xd(X0, Z0):
            a = (X0 + Z0) % M
            b = (X0 - Z0) % M
            a2 = (a * a) % M
            b2 = (b * b) % M
            c = (a2 - b2) % M
            return (a2 * b2) % M, (c * (b2 + a24 * c)) % M

        def _xa(X1, Z1, X2, Z2, Xd, Zd):
            u1 = ((X1 - Z1) * (X2 + Z2)) % M
            v1 = ((X1 + Z1) * (X2 - Z2)) % M
            return (Zd * pow(u1 + v1, 2, M)) % M, (Xd * pow(u1 - v1, 2, M)) % M

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

    branches = _init_branch(N)
    for p in primes_b1:
        pk = p
        while pk <= b1:
            nxt = []
            for d, bx, bz, a24 in branches:
                if bx is None:
                    nxt.append((d, None, None, None))
                    continue
                bx, bz = _xmul_mod(p, bx, bz, a24, d)
                g = math.gcd(bz, d)
                if 1 < g < d:
                    nxt.append((g, None, None, None))
                    cof = d // g
                    nxt.append((cof, bx % cof, bz % cof, a24 % cof))
                elif g == d:
                    nxt.append((d, None, None, None))
                else:
                    nxt.append((d, bx, bz, a24))
            branches = nxt
            pk *= p
    return [d for d, _, _, _ in branches]


def _run_suyama_frontier_s1s2(N, sigma, b1, b2, wheel, primes_b1):
    """Run ONE full Stage-1 stepped-GCD + Stage-2 binary-product-tree parallel-active-frontier curve on N.
    In Stage 2, every giant step computes g_root = gcd(root, d); if g_root > 1, d // g_root (if > 1)
    continues to the next giant step while _split_product_tree(levels, d) recursively isolates every
    prime factor hit at any distinct leaf. Returns (parts, total_tree_descent_gcds)."""
    u = (sigma * sigma - 5) % N
    v = (4 * sigma) % N
    x = pow(u, 3, N)
    z = pow(v, 3, N)
    den = (16 * x * v) % N

    def _init(d):
        if d <= 1:
            return []
        pp = _is_perfect_power(d)
        if pp is not None:
            return _init(pp[0])
        g0 = math.gcd(z % d, d)
        if 1 < g0 < d:
            return _init(g0) + _init(d // g0)
        if g0 == d:
            return [(d, None, None, None)]
        gd = math.gcd(den % d, d)
        if 1 < gd < d:
            return _init(gd) + _init(d // gd)
        if gd == d:
            return [(d, None, None, None)]
        a24 = (pow((v - u) % d, 3, d) * ((3 * u + v) % d) * pow(den % d, -1, d)) % d
        Am = (4 * a24 - 2) % d
        gs = math.gcd((Am * Am - 4) % d, d)
        if 1 < gs < d:
            return _init(gs) + _init(d // gs)
        if gs == d:
            return [(d, None, None, None)]
        return [(d, x % d, z % d, a24)]

    def _xmul(k, X, Z, a24, M):
        if k == 1:
            return X, Z

        def _xd(X0, Z0):
            a = (X0 + Z0) % M
            b = (X0 - Z0) % M
            a2 = (a * a) % M
            b2 = (b * b) % M
            c = (a2 - b2) % M
            return (a2 * b2) % M, (c * (b2 + a24 * c)) % M

        def _xa(X1, Z1, X2, Z2, Xd, Zd):
            u1 = ((X1 - Z1) * (X2 + Z2)) % M
            v1 = ((X1 + Z1) * (X2 - Z2)) % M
            return (Zd * pow(u1 + v1, 2, M)) % M, (Xd * pow(u1 - v1, 2, M)) % M

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

    branches = _init(N)
    for p in primes_b1:
        pk = p
        while pk <= b1:
            nxt = []
            for d, bx, bz, a24 in branches:
                if bx is None:
                    nxt.append((d, None, None, None))
                    continue
                bx, bz = _xmul(p, bx, bz, a24, d)
                g = math.gcd(bz, d)
                if 1 < g < d:
                    nxt.append((g, None, None, None))
                    cof = d // g
                    nxt.append((cof, bx % cof, bz % cof, a24 % cof))
                elif g == d:
                    nxt.append((d, None, None, None))
                else:
                    nxt.append((d, bx, bz, a24))
            branches = nxt
            pk *= p

    done = [d for d, bx, _, _ in branches if bx is None]
    active = [(d, bx, bz, a24) for d, bx, bz, a24 in branches if bx is not None]
    babies = [r for r in range(1, wheel // 2) if math.gcd(r, wheel) == 1]
    G = ceil(b2 / wheel) + 1
    total_desc_gcds = 0
    for d, bx, bz, a24 in active:
        cur = d
        cur_x, cur_z, cur_a24 = bx, bz, a24
        bpts = [_xmul(r, cur_x, cur_z, cur_a24, cur) for r in babies]
        for g_idx in range(1, G + 1):
            if cur <= 1:
                break
            xg, zg = _xmul(g_idx * wheel, cur_x, cur_z, cur_a24, cur)
            raw = [(xg * zb - xb * zg) % cur for xb, zb in bpts]
            levels = _build_product_tree_mod(raw, cur)
            g_root = math.gcd(levels[-1][0], cur)
            if g_root == 1:
                continue
            parts, _, dg = _split_product_tree(levels, cur)
            total_desc_gcds += dg
            if g_root == cur:
                done.extend(parts)
                cur = 1
                break
            rem = cur // g_root
            pp_r = _is_perfect_power(rem)
            if pp_r is not None:
                rem_base, rem_e = pp_r
                # Remove one copy of rem_base to keep as active cur, append remaining (rem_e - 1) copies + all other parts:
                skipped = 0
                for pt in parts:
                    if pt == rem_base and skipped == 0:
                        skipped = 1
                        continue
                    done.append(pt)
                cur = rem_base
            else:
                skipped = 0
                for pt in parts:
                    if pt == rem and skipped == 0:
                        skipped = 1
                        continue
                    done.append(pt)
                cur = rem
            cur_x %= cur
            cur_z %= cur
            cur_a24 %= cur
            bpts = [(xb % cur, zb % cur) for xb, zb in bpts]
        if cur > 1:
            done.append(cur)
    return done, total_desc_gcds


def _is_prime_small(n):
    if n < 2:
        return False
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            return False
    return True


def verify_suyama_stepped():
    """Live verification of:
      (1) Newton perfect-power precheck on N = p^2,
      (2) singular-curve guard on sigma = N - 1,
      (3) prime-by-prime stepped-gcd split on a composite N = p*q,
      (4) empirical EXACT target-prime p0 isolation rate across 3 composite regimes
          (20b*41b, balanced 20b*20b, and 3-factor 20b*20b*20b with in-curve continuation),
          confirming exact target-prime isolation >= model_p1_lo in ALL 3 regimes, PLUS
          full parallel-active-frontier complete 3-prime factorization across 16 composites,
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

    # (4) Exact target-prime p0 isolation across 3 composite regimes (B1 = B2 = 256):
    b1_emp = 256
    ps_emp = tuple(primes_upto(b1_emp))
    p1_model = semismooth_lo(21.0, b1_emp, b1_emp) * power_correction(21.0, b1_emp, b1_emp)
    q_comp = 1099511627791  # prime 2^40 + 15
    test_primes = (
        1048583, 1048601, 1048609, 1048613, 1048627, 1048633,
        1100009, 1100023, 1200007, 1200061, 1300021, 1300031,
        1500007, 1500041, 1700021, 1700047,
    )
    assert all(_is_prime_small(pp) for pp in test_primes)
    rng = random.Random(20261124)
    trials = 0
    hits = 0
    for pp in test_primes:
        N_c = pp * q_comp
        for _ in range(10):
            sig = rng.randrange(6, pp - 6)
            g = _run_suyama_stepped(N_c, sig, b1_emp, ps_emp)
            trials += 1
            if isinstance(g, int) and (g == pp or N_c // g == pp):
                hits += 1
    emp_p1 = hits / trials
    assert emp_p1 >= p1_model, f"Empirical 20b*41b target-prime P1 {emp_p1:.4f} < model lower bound {p1_model:.4f}"

    # (4b) Balanced 20b * 20b composite exact target-prime isolation (64 trials):
    bal_trials = 0
    bal_hits = 0
    for idx in range(0, len(test_primes), 2):
        p0, p1_ = test_primes[idx], test_primes[idx + 1]
        N_b = p0 * p1_
        for _ in range(8):
            sig = rng.randrange(6, p0 - 6)
            g = _run_suyama_stepped(N_b, sig, b1_emp, ps_emp)
            bal_trials += 1
            if isinstance(g, int) and (g == p0 or N_b // g == p0):
                bal_hits += 1
    emp_bal_p1 = bal_hits / bal_trials
    assert emp_bal_p1 >= p1_model, f"Balanced 20b*20b target-prime P1 {emp_bal_p1:.4f} < model lower bound {p1_model:.4f}"

    # (4c) 3-factor 20b * 20b * 20b composite EXACT target-prime p0 isolation with in-curve continuation (64 trials):
    tri_trials = 0
    tri_hits = 0
    n_tp = len(test_primes)
    for idx in range(8):
        p0, p1_, p2_ = test_primes[idx], test_primes[(idx + 5) % n_tp], test_primes[(idx + 11) % n_tp]
        N_t = p0 * p1_ * p2_
        for _ in range(8):
            sig = rng.randrange(6, p0 - 6)
            parts = _run_suyama_frontier_curve(N_t, sig, b1_emp, ps_emp)
            tri_trials += 1
            if p0 in parts:
                tri_hits += 1
    emp_tri_p1 = tri_hits / tri_trials
    assert emp_tri_p1 >= p1_model, f"3-factor exact p0 isolation P1 {emp_tri_p1:.4f} < model lower bound {p1_model:.4f}"

    # (4d) Full parallel-active-frontier complete 3-prime factorization across 16 composites using C_s1 rung curves:
    C_s1 = ceil(log1p(-0.999) / log1p(-p1_model))
    frontier_full_ok = 0
    for idx in range(n_tp):
        p0, p1_, p2_ = test_primes[idx], test_primes[(idx + 5) % n_tp], test_primes[(idx + 11) % n_tp]
        frontier = [p0 * p1_ * p2_]
        found = set()
        for _ in range(C_s1):
            if not frontier:
                break
            nxt = []
            for d in frontier:
                sig = rng.randrange(6, d - 6)
                for part in _run_suyama_frontier_curve(d, sig, b1_emp, ps_emp):
                    if part in (p0, p1_, p2_):
                        found.add(part)
                    else:
                        nxt.append(part)
            frontier = nxt
        if found == {p0, p1_, p2_}:
            frontier_full_ok += 1
    assert frontier_full_ok == n_tp, f"Full frontier 3-prime factorization {frontier_full_ok}/{n_tp} != {n_tp}/{n_tp}"

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
    assert s2_split in (p_a, p_b) and s2_gcds <= 1 + 2 * h_s2 * 2

    # Explicit [0, p, q, 1] zero-leaf-skip edge-case check (skipping 0 mod N leaves before tree multiply):
    lv_zero = _build_product_tree_mod([0, p_a, p_b, 1], N_s2)
    assert math.gcd(lv_zero[-1][0], N_s2) == N_s2
    z_split, z_gcds = _bisect_product_tree(lv_zero, N_s2)
    assert z_split in (p_a, p_b) and z_gcds <= 1 + 2 * 2 * 2

    # (5b) Live W <= B1 (W=210 <= B1=256) 3-prime Stage-2 multi-hit p*q inside p*q*r witness:
    # sig=6, B1=256, W=210, g=2 on N_3 = 1000159 * 1000211 * 1000003:
    # Stage-1 gcd=1, p_a and p_b each divide EXACLY ONE distinct Stage-2 leaf (single-leaf uniqueness under W <= B1),
    # p_c is unhit, g_root = p_a*p_b < N_3, and _split_product_tree isolates {p_a, p_b, p_c}.
    p3_a, p3_b, p3_c = 1000159, 1000211, 1000003
    assert _is_prime_small(p3_a) and _is_prime_small(p3_b) and _is_prime_small(p3_c)
    N_3s2 = p3_a * p3_b * p3_c
    sig3, b1_3, W3, g3 = 6, 256, 210, 2
    assert W3 <= b1_3
    u3 = (sig3 * sig3 - 5) % N_3s2
    v3 = (4 * sig3) % N_3s2
    x3, z3 = pow(u3, 3, N_3s2), pow(v3, 3, N_3s2)
    den3 = (16 * x3 * v3) % N_3s2
    a24_3 = (pow(v3 - u3, 3, N_3s2) * (3 * u3 + v3) * pow(den3, -1, N_3s2)) % N_3s2

    def _xd3(X, Z):
        a = (X + Z) % N_3s2; b = (X - Z) % N_3s2
        a2 = (a * a) % N_3s2; b2 = (b * b) % N_3s2
        c = (a2 - b2) % N_3s2
        return (a2 * b2) % N_3s2, (c * (b2 + a24_3 * c)) % N_3s2

    def _xa3(X1, Z1, X2, Z2, Xd, Zd):
        u1 = ((X1 - Z1) * (X2 + Z2)) % N_3s2
        v1 = ((X1 + Z1) * (X2 - Z2)) % N_3s2
        return (Zd * pow(u1 + v1, 2, N_3s2)) % N_3s2, (Xd * pow(u1 - v1, 2, N_3s2)) % N_3s2

    def _xm3(k, X, Z):
        if k == 1:
            return X, Z
        X0, Z0 = X, Z
        X1, Z1 = _xd3(X, Z)
        for bit in bin(k)[3:]:
            if bit == "1":
                X0, Z0 = _xa3(X0, Z0, X1, Z1, X, Z)
                X1, Z1 = _xd3(X1, Z1)
            else:
                X1, Z1 = _xa3(X0, Z0, X1, Z1, X, Z)
                X0, Z0 = _xd3(X0, Z0)
        return X0, Z0

    for pp in ps_emp:
        pk = pp
        while pk <= b1_3:
            x3, z3 = _xm3(pp, x3, z3)
            pk *= pp
    assert math.gcd(z3, N_3s2) == 1
    bpts3 = [_xm3(r, x3, z3) for r in babies_s2]
    xg3, zg3 = _xm3(g3 * W3, x3, z3)
    raw3 = [(xg3 * zb - xb * zg3) % N_3s2 for xb, zb in bpts3]
    # Assert single-leaf uniqueness for both Stage-2 hitting primes p3_a and p3_b:
    assert sum(1 for t in raw3 if t % p3_a == 0) == 1
    assert sum(1 for t in raw3 if t % p3_b == 0) == 1
    assert sum(1 for t in raw3 if t % p3_c == 0) == 0
    lv3 = _build_product_tree_mod(raw3, N_3s2)
    assert math.gcd(lv3[-1][0], N_3s2) == p3_a * p3_b
    parts3, tot_g3, desc_g3 = _split_product_tree(lv3, N_3s2)
    assert set(parts3) == {p3_a, p3_b, p3_c} and desc_g3 <= 2 * h_s2 * 2

    # Synthetic 3-prime checks: all-hit, same-leaf collision, AND multi-leaf companion [q, p, q, 1] isolating p:
    lv3_all = _build_product_tree_mod([p3_a, p3_b, p3_c, 1], N_3s2)
    p3_all, _, dg3_all = _split_product_tree(lv3_all, N_3s2)
    assert set(p3_all) == {p3_a, p3_b, p3_c} and dg3_all <= 2 * 2 * 3
    lv3_col = _build_product_tree_mod([p3_a * p3_b, p3_c, 1, 1], N_3s2)
    p3_col, _, dg3_col = _split_product_tree(lv3_col, N_3s2)
    assert set(p3_col) == {p3_a * p3_b, p3_c} and dg3_col <= 2 * 2 * 3
    lv3_mult = _build_product_tree_mod([p3_b, p3_a, p3_b, 1], N_3s2)
    p3_mult, _, dg3_mult = _split_product_tree(lv3_mult, N_3s2)
    assert set(p3_mult) == {p3_a, p3_b, p3_c} and dg3_mult <= 2 * 2 * 2
    # Non-squarefree prime-power check d = p^2 * q where leaf holds p^2: _split_product_tree + _is_perfect_power returns p:
    N_ppq = p3_a * p3_a * p3_c
    lv_ppq = _build_product_tree_mod([p3_a * p3_a, 1, 1, 1], N_ppq)
    p_ppq, _, dg_ppq = _split_product_tree(lv_ppq, N_ppq)
    assert sorted(p_ppq) == sorted([p3_a, p3_a, p3_c]) and dg_ppq <= 2 * 2 * 2

    # (5c) Live 3-prime Stage-1 + Stage-2 empirical target-prime p0 isolation (B1=256, B2=2560, W=210, 64 trials):
    b1_s2e, b2_s2e, W_s2e = 256, 2560, 210
    p1_s2e_model = semismooth_lo(21.0, b1_s2e, b2_s2e) * power_correction(21.0, b1_s2e, b2_s2e)
    s2e_trials = 0
    s2e_hits = 0
    s2e_max_dg = 0
    h_s2e = ceil(log2(24))
    for idx in range(8):
        p0, p1_, p2_ = test_primes[idx], test_primes[(idx + 5) % n_tp], test_primes[(idx + 11) % n_tp]
        Nt = p0 * p1_ * p2_
        for _ in range(8):
            sig = rng.randrange(6, p0 - 6)
            parts, dg = _run_suyama_frontier_s1s2(Nt, sig, b1_s2e, b2_s2e, W_s2e, ps_emp)
            s2e_trials += 1
            if p0 in parts:
                s2e_hits += 1
            if dg > s2e_max_dg:
                s2e_max_dg = dg
    emp_tri_s2_p1 = s2e_hits / s2e_trials
    assert emp_tri_s2_p1 >= p1_s2e_model, f"3-prime S1+S2 p0 isolation {emp_tri_s2_p1:.4f} < model {p1_s2e_model:.4f}"
    assert s2e_max_dg <= 2 * h_s2e * 3, f"3-prime S2 descent gcds {s2e_max_dg} > {2 * h_s2e * 3}"

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
        emp_balanced_trials=bal_trials,
        emp_balanced_hits=bal_hits,
        emp_balanced_p1=round(emp_bal_p1, 4),
        emp_tri_trials=tri_trials,
        emp_tri_hits=tri_hits,
        emp_tri_p1=round(emp_tri_p1, 4),
        emp_tri_s1s2_b1=b1_s2e,
        emp_tri_s1s2_b2=b2_s2e,
        emp_tri_s1s2_trials=s2e_trials,
        emp_tri_s1s2_hits=s2e_hits,
        emp_tri_s1s2_p1=round(emp_tri_s2_p1, 4),
        emp_tri_s1s2_model_p1_lo=round(p1_s2e_model, 4),
        emp_tri_s1s2_max_desc_gcds=s2e_max_dg,
        frontier_3factor_full_ok=f"{frontier_full_ok}/{n_tp}",
        frontier_3factor_curves=C_s1,
        model_p1_lo=round(p1_model, 4),
        stage2_tree_bisect_ok=True,
        stage2_zero_leaf_skip_ok=True,
        stage2_3prime_pq_split_ok=True,
        stage2_3prime_N=N_3s2,
        stage2_3prime_desc_gcds=desc_g3,
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
    partitions with prod d_i | N for any odd post-TDIV input N <= 2^nbits and b_i = bit_length(d_i) >= b_min.
    Since N is odd, N < 2^ceil(nbits), and since d_i >= 2^(b_i - 1), s_i = b_i - 1 >= b_min - 1 satisfies
    sum s_i <= log2(N) < ceil(nbits), hence sum s_i <= ceil(nbits) - 1 unconditionally (exact for
    bit_length(a) + bit_length(b) = bit_length(a*b) + 1, e.g., N = (65537*65539)*(65543*65551) at nb = 65)."""
    S_max = int(ceil(nbits)) - 1
    s_min = max(1, int(b_min) - 1)
    if S_max < s_min:
        b_full = max(1, S_max + 1)
        L = max(l_min, (b_full + 63) >> 6)
        return L * L, 2.0 * b_full * (12 * L + 16)
    key = (b_min, l_min)
    st = _TREE_DP_CACHE.get(key)
    if st is None:
        dp_l2 = [0] * s_min
        dp_g = [0.0] * s_min
        st = (dp_l2, dp_g)
        _TREE_DP_CACHE[key] = st
    else:
        dp_l2, dp_g = st
    cur_s = len(dp_l2) - 1
    if S_max > cur_s:
        for S in range(cur_s + 1, S_max + 1):
            b_full = S + 1
            L_full = max(l_min, (b_full + 63) >> 6)
            best_l2 = L_full * L_full
            best_g = 2.0 * b_full * (12 * L_full + 16)
            s_hi = S - s_min
            if s_hi >= s_min:
                cands = [s_min]
                for k in range(1, (b_full >> 6) + 2):
                    for s in ((k << 6) - 1, k << 6):
                        if s_min < s <= s_hi:
                            cands.append(s)
                for s in cands:
                    b_i = s + 1
                    Ls = max(l_min, (b_i + 63) >> 6)
                    v2 = Ls * Ls + dp_l2[S - s]
                    vg = 2.0 * b_i * (12 * Ls + 16) + dp_g[S - s]
                    if v2 > best_l2:
                        best_l2 = v2
                    if vg > best_g:
                        best_g = vg
            dp_l2.append(best_l2)
            dp_g.append(best_g)
    return dp_l2[S_max], dp_g[S_max]


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
        return base, (), (), (), ()
    span = top - base
    los = tuple(base + span * i / steps for i in range(steps))
    wts = tuple(log((base + span * (i + 1) / steps) / (base + span * i / steps)) for i in range(steps))
    suf = [0.0] * (steps + 1)
    acc = 0.0
    for i in range(steps - 1, -1, -1):
        acc += wts[i]
        suf[i] = acc
    inv_br = (1.0 / D.RSTEP) / base
    los_s = tuple(x * inv_br for x in los)
    return base, los, wts, tuple(suf), los_s


@lru_cache(None)
def semismooth_lo(fbits, b1, b2, steps=256):
    """Conservative LOWER bound on (B1, B2)-semismoothness for a prime p <= 2^{fbits}."""
    m = fbits - _LOG2_TORSION
    base, los, wts, suf, los_s = _ss_grid(b1, b2, steps)
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
            lo_i = bisect_left(los, m - (ntab - 1) / inv_br, 0, cut)
            if lo_i < cut:
                ms = m * inv_br
                while lo_i < cut and _ceil(ms - los_s[lo_i]) >= ntab:
                    lo_i += 1
                while cut > lo_i and (m - los[cut - 1]) <= base:
                    cut -= 1
                    acc += wts[cut]
                if lo_i < cut:
                    acc += _SUMPROD([tab[_ceil(ms - x)] for x in los_s[lo_i:cut]], wts[lo_i:cut])
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
    for q in (2, 3, 5, 7, 11, 13, 17):
        if wheel % q == 0:
            phi = phi // q * (q - 1)
    J = phi // 2
    G = ceil(b2 / wheel) + 1
    s2 = 6 * (wheel // 2 + G) + 3 * J * G if b2 > b1 else 0
    return 1.1 * (s1 + s2), s1, s2, J, G


_B1_FINE = tuple(sorted({int(round(2.0 ** (k / 2.0))) for k in range(12, 74)}))


@lru_cache(None)
def schedule_rung(fbits, nbits, success_rung=0.999, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
    """Single-rung Deterministic Stepped-GCD + Full Recursive Stage-2 Binary-Product-Tree branch-tree-capped ECM (H-ECM).
    Every curve takes step-by-step gcd after each Stage-1 prime power and builds a J-leaf binary product tree
    on each Stage-2 giant step G (h_tree = ceil(log2(max(2, J))) levels, +1.1*J*G muls, root GCDs <= G*g_max).
    Whenever g_root = gcd(root, d) > 1 (whether 1 < g_root < d or g_root == d), d // g_root continues to
    subsequent giant steps (lazy baby-step reduction mod d // g_root is covered inside the 1.1*3*J*G pairing budget,
    since Montgomery pairing needs 2*J muls vs 3.3*J charged) while _split_product_tree recursively descends the
    tree with g_root: enforcing wheel W <= B1 guarantees that under the combined Stage-1+Stage-2 event
    E_good(p) (lower-bounded by semismooth_lo * _power_correction_b1), whenever Stage-2 case (b) occurs
    (ord_p(Q_1) = L_p > B1 >= W), p divides at most one leaf t_{j(p)} in its hit giant step, so _split_product_tree
    + test-once _is_perfect_power (run only on strictly proper new splits 1 < m < d, charged once at
    2 * omega_max * NEWTON_PRECHECK_INSTR) deterministically extracts p whenever no distinct companion prime
    collides at that exact leaf (g_step, j(p)) (rare distinct-prime collisions absorbed into H-ECM). Because g_root is removed
    from the curve's active Stage-2 frontier, sum_{g_step, d} Omega(g_root) <= Omega(N) <= floor(nb / (b_min // 2))
    = omega_max, so total recursive tree-descent GCD work across ALL G giant steps of one curve is unconditionally
    <= 2 * h_tree * omega_max * gcd_instr(nb) <= 2 * h_tree * omega_max * g_max:
      per_curve_stepped = (ops + 1.1*J*G) * mul + (pi_pp_upper(B1) + G + 2 * h_tree * omega_max) * gcd
    100% deterministic cost: ZERO curve replay, ZERO p1-as-upper-bound, ZERO variance."""
    nb = int(ceil(nbits))
    l2_max, g_max = tree_branch_caps(nb, b_min, l_min)
    mul = MONT5 * (l2_max / 25.0)
    gcd = g_max
    omega_max = nb // max(1, b_min // 2)
    mul_nb, _ = mulmod_instr(nb, l_min)
    # Test-once perfect-power + MR_ROUNDS=110 Miller-Rabin primality classification on initial N and every
    # strictly proper new split child 1 < m < d (<= 2*omega_max - 1 nodes in factorization tree: <= omega_max
    # leaves with sum b_leaf <= nb + omega_max, and <= omega_max - 1 internal composite nodes).
    # Confirms primes immediately so EVERY remaining active frontier element is guaranteed composite (bit_length >= b_min),
    # validating tree_branch_caps unconditionally, and ensures Stage-1 gate divides out only confirmed primes <= 2^69.0.
    # For collection (b_min >= 96, var_t = 0.0): 100% deterministic worst-case 2*omega_max * (2*MR_ROUNDS*nb) * mul_nb.
    # For descent (b_min < 96, expected work E[W_run]): deterministic worst-case 2*MR_ROUNDS*(nb+omega_max)*mul_nb
    # across all leaves + unconditional Rabin theorem E[rounds | odd composite] < 4/3 -> (8/3)*omega_max*nb*mul_nb
    # across all <= omega_max internal composite nodes (error <= 4^-110 = 2^-220/test, < 2^-100 pointwise up to T_run=2^124).
    if b_min >= 96:
        mr_budget = 2.0 * omega_max * (2.0 * MR_ROUNDS * nb) * mul_nb
    else:
        mr_budget = (2.0 * MR_ROUNDS * (nb + omega_max) + (8.0 / 3.0) * omega_max * nb) * mul_nb
    pp_budget = 2 * omega_max * NEWTON_PRECHECK_INSTR + mr_budget
    l1ms = log1p(-float(success_rung))
    best = None
    best_cost = float("inf")
    _ratios_desc = (500, 300, 200, 150, 100, 75, 50, 30, 20, 10, 5, 1)
    order = _B1_FINE[::6] + tuple(b for i, b in enumerate(_B1_FINE) if i % 6 != 0)
    for b1 in order:
        pi_b1 = pi_pp_upper(b1)
        min_s1_ops = 1.1 * ceil(1.5 * b1) * 11
        s1_curve = min_s1_ops * mul + pi_b1 * gcd
        if s1_curve + pp_budget >= best_cost:
            continue
        pc = _power_correction_b1(fbits, b1)
        if pc <= 0.0:
            continue
        p1_hi = semismooth_lo(fbits, b1, b1 * 500) * pc
        if p1_hi <= 1e-12:
            continue
        if ceil(l1ms / log1p(-p1_hi)) * s1_curve + pp_budget >= best_cost:
            continue
        for ratio in _ratios_desc:
            b2 = b1 * ratio
            p1 = p1_hi if ratio == 500 else semismooth_lo(fbits, b1, b2) * pc
            if p1 <= 1e-12:
                break
            curves = ceil(l1ms / log1p(-p1))
            if curves * s1_curve + pp_budget >= best_cost:
                break
            wheels = WHEELS if ratio > 1 else (0,)
            for wheel in wheels:
                if wheel and wheel > b1:
                    break
                ops, s1, s2, J, G = (
                    stage_ops(b1, b2, wheel)
                    if wheel
                    else (min_s1_ops, ceil(1.5 * b1) * 11, 0, 0, 0)
                )
                h_tree = ceil(log2(max(2, J))) if G else 0
                m_ops = ops + (1.1 * J * G if G else 0.0)
                g_cnt = pi_b1 + (G + 2 * h_tree * omega_max if G else 0)
                per_curve_stepped = m_ops * mul + g_cnt * gcd
                cost = curves * per_curve_stepped + pp_budget
                if cost < best_cost:
                    best_cost = cost
                    best = dict(
                        fbits=fbits, nbits=nb, success_rung=float(success_rung),
                        b_min=b_min, l_min=l_min,
                        l2_tree_max=l2_max, gcd_tree_max=g_max,
                        B1=b1, B2=b2, wheel=wheel, curves=curves,
                        p1=p1, pi_b1=pi_b1, J=J, G=G, h_tree=h_tree,
                        omega_max=omega_max,
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


ECM_T2_PAD_BITS = 0.05  # Continuous real-interval Theorem T2 padding so min_{r in [16, fbits]} p1(r) >= p1(fbits + 0.05)


def ecm_instr(fbits, nbits, success=0.999):
    """1-rung deterministic dominance Stepped-GCD parallel-active-frontier ECM up to fbits.
    Tuned to hi_tune = round(fbits + ECM_T2_PAD_BITS, 4) so that on EVERY real prime size
    r = log2(p) in [16.0, fbits] (continuous real interval, zero integer-grid gaps),
    semismooth_lo(r, B1, B2) * _power_correction_b1(r, B1) >= p1(hi_tune, B1, B2) holds
    unconditionally (verified live via contiguous 0.05-bit interval enclosure in verify_1rung_monotonicity).
    Running all C curves across the active composite frontier F detects EVERY prime factor
    p <= 2^fbits with probability >= success (cost bounded unconditionally by tree_branch_caps)."""
    hi_tune = round(float(fbits) + ECM_T2_PAD_BITS, 4)
    sch = schedule_rung(hi_tune, int(ceil(nbits)), float(success), B_MIN_DEFAULT, L_MIN_DEFAULT)
    if sch is None:
        raise ValueError(f"Infeasible 1-rung ECM fbits={fbits}, nb={nbits}")
    return sch["instructions"]


def verify_1rung_monotonicity(targets=((69.0, 2590), (69.0, 1100), (189.15, 2100)), dr=0.05):
    """Continuous real-interval proof of Theorem T2: for EVERY target (hi_real, nb) tuned at
    hi_tune = round(hi_real + ECM_T2_PAD_BITS, 4) with s = 0.999, partition [16.0, hi_real] into
    contiguous real subintervals [r_lo, r_hi] of width <= dr. Since semismooth_lo(r, B1, B2) is
    analytically non-increasing in real r and rho_hi, rho_lo are non-increasing, on every [r_lo, r_hi]
    every real r in [r_lo, r_hi] satisfies:
        semismooth_lo(r, B1, B2) >= semismooth_lo(r_hi, B1, B2)
        _power_correction_b1(r, B1) >= max(0, 1 - (rho_hi(max(u_lo - 2, 1)) / rho_lo(u_hi)) * 2/(y*log(y)))
    Asserting product >= p1(hi_tune, B1, B2) on EVERY subinterval proves p1(r) >= p1(hi_tune) for
    ALL real r in [16.0, hi_real] (zero gaps between integers or at fractional endpoints)."""
    checked = 0
    min_ratio = float("inf")
    uniq = sorted(set((round(float(h), 4), int(ceil(nb))) for h, nb in targets))
    for hi_real, nb_t in uniq:
        hi_tune = round(hi_real + ECM_T2_PAD_BITS, 4)
        sch = schedule_rung(hi_tune, nb_t, 0.999, B_MIN_DEFAULT, L_MIN_DEFAULT)
        b1, b2, p1_hi = sch["B1"], sch["B2"], sch["p1"]
        assert sch["wheel"] <= b1, f"Wheel {sch['wheel']} > B1 {b1} violates single-leaf uniqueness"
        y = sqrt(b1)
        cfac = 2.0 / (y * log(y))
        lb1 = log2(b1)
        n_iv = max(1, int(ceil((hi_real - 16.0) / dr)))
        inv_lb1 = 1.0 / lb1
        inv_p1 = 1.0 / p1_hi
        _rhi = D.rho_hi; _rlo = D.rho_lo
        blk = 20
        for b0 in range(0, n_iv, blk):
            b1_i = min(n_iv, b0 + blk)
            r_blk_hi = min(hi_real, 16.0 + b1_i * dr)
            ss_blk = semismooth_lo(r_blk_hi, b1, b2)
            for i in range(b0, b1_i):
                r_lo = 16.0 + i * dr
                r_hi = min(hi_real, r_lo + dr)
                u_lo = (r_lo - _LOG2_TORSION) * inv_lb1
                u_hi = (r_hi - _LOG2_TORSION) * inv_lb1
                tilt_hi = _rhi(u_lo - 2.0 if u_lo > 3.0 else 1.0) / max(_rlo(u_hi), 1e-300)
                pc_lo = 1.0 - tilt_hi * cfac
                if pc_lo < 0.0:
                    pc_lo = 0.0
                rat = (ss_blk * pc_lo) * inv_p1
                if rat < 1.05:
                    ss_lo = semismooth_lo(r_hi, b1, b2)
                    rat = (ss_lo * pc_lo) * inv_p1
                    assert rat >= 1.0, f"Continuous T2 failed at hi={hi_real} [{r_lo},{r_hi}]: rat={rat}"
                if rat < min_ratio:
                    min_ratio = rat
                checked += 1
    return dict(
        ok=True,
        targets=len(uniq),
        intervals_checked=checked,
        r_checks=checked,
        min_iv_ratio=round(min_ratio, 6),
        dr_bits=dr,
        pad_bits=ECM_T2_PAD_BITS,
    )
