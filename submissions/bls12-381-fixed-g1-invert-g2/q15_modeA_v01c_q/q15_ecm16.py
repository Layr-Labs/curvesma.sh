"""q15_ecm16.py — Deterministic Stepped-GCD Suyama-Montgomery ECM schedule (STEPPED_LADDER).
Every curve pays step-by-step gcd deterministically (per_curve_stepped = ops * mul + (pi_upper(B1) + G) * gcd)
over exact L_min = 2 branch-tree DP caps (tree_branch_caps), with 12-bit geometric ladder union-bounded
across rungs under standard literature H-ECM (Brent / Zimmermann-Dodson) random-integer group orders."""
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
    """Run one Suyama-Montgomery Stage-1 curve mod N with step-by-step gcd after each prime power."""
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
        while pk * p <= b1:
            pk *= p
        x, z = xmul(pk, x, z)
        g = math.gcd(z, N)
        if g > 1:
            return g
    return 1


def verify_suyama_ecm_live():
    """Live verification of:
      (1) Newton perfect-power precheck on N = p^2,
      (2) singular-curve guard on sigma = N - 1,
      (3) clean stepped-gcd factor splitting on composites N = p * q_large (20-bit p, B1 = 256)
          confirming observed Stage-1 hit rate >= semismooth(20.0, 256, 256) * _power_correction_b1(20.0, 256)."""
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
    return dict(
        N=N, perfect_power_precheck_ok=True, singular_nm1_rejected=True,
        clean_split_sigma=split_sigma, split_factor=split_fac,
        emp_b1=b1_emp, emp_trials=trials, emp_found=found,
        emp_hit_rate=round(emp_rate, 4), model_p1_lo=round(p_mod_lo, 4),
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
    """Precompute all feasible (b1, b2, wheel, ops, gcd_steps, p1, log1m_p1) for stepped-gcd ECM."""
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
                ops, _, _, _, G = stage_ops(b1, b2, wheel)
            else:
                wheel = 0
                ops, G = min_s1_ops, 0
            out.append((b1, b2, wheel, ops, pi_b1 + G, p1, l1m))
    return tuple(out)


@lru_cache(None)
def schedule_rung(fbits, nbits, success=0.999, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
    """Single-rung deterministic stepped-gcd branch-tree-capped ECM schedule (H-ECM).
    Every curve pays step-by-step gcd deterministically: per_curve_stepped = ops * mul + (pi_upper(B1) + G) * gcd.
    Curve count follows standard literature H-ECM: curves = ceil(log1p(-success) / log1p(-p1))."""
    nb = int(ceil(nbits))
    l2_max, g_max = tree_branch_caps(nb, b_min, l_min)
    mul = MONT5 * (l2_max / 25.0)
    gcd = g_max
    l1ms = log1p(-success)
    best = None
    best_cost = float("inf")
    for b1, b2, wheel, ops, g_steps, p1, l1m in _rung_candidates_stepped(fbits):
        per_curve = ops * mul + g_steps * gcd
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
                gcd_steps=g_steps, per_curve_stepped=per_curve,
                instructions=cost, log2_instructions=log2(cost),
            )
    return best


schedule = schedule_rung


@lru_cache(None)
def ladder_instr(lo_bits, hi_bits, nbits, success=0.999, step=12.0, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
    """Bottom-up deterministic stepped-gcd 12-bit geometric ECM ladder from lo_bits to hi_bits.
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
    """Deterministic stepped-gcd descent ECM ladder from 16b up to fbits (STEPPED_LADDER)."""
    return ladder_instr(16.0, fbits, int(ceil(nbits)), success=success, step=12.0, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT)
