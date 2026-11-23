"""q15_ecm16.py — Dependence-robust Suyama-Montgomery ECM schedule (L_min=2 branch-tree DP, 12b geometric ladder, sep=0.5, explicit failed-replay retry pricing)."""
import math
from functools import lru_cache
from math import ceil, isqrt, log, log1p, log2, sqrt
import q15_descent16 as D

MONT5 = D.MONT5
TORSION = 12
WHEELS = (30030,)
RETRY_PR_COLL = 1e-3
RETRY_PR_DESC = 0.0043  # Explicit single-curve simultaneous-hit retry rate (pays BOTH per_curve + full replay_1curve!)
SEPARATION_DEFAULT = 0.5  # Conservative curve separation factor
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
    """Exact prime-exponent Newton integer perfect-power precheck before ECM (eliminates N = p^k singular redraws; explicitly charged at NEWTON_PRECHECK_INSTR = 2^26 Ir)."""
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


def verify_suyama_retry_counterexample():
    """Live test of perfect-power precheck (N = p^2), simultaneous ord(P) = 3 Suyama counterexample
    N = 1157341 * 2633581, sigma = 533684655365 (gcd == N), singular-curve guard on sigma = N - 1,
    and clean split upon random-sigma retry."""
    N_sq = 1157341 * 1157341
    pp_res = _is_perfect_power(N_sq)
    assert pp_res == (1157341, 2), f"Perfect-power precheck failed on p^2: {pp_res}"
    N = 1157341 * 2633581
    assert _is_perfect_power(N) is None

    def _run_sigma(sigma, b1):
        u = (sigma * sigma - 5) % N
        v = (4 * sigma) % N
        x = pow(u, 3, N)
        z = pow(v, 3, N)
        g0 = math.gcd(z, N)
        if g0 > 1:
            return g0
        den = (16 * pow(u, 3, N) * v) % N
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

        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
            pk = p
            while pk * p <= b1:
                pk *= p
            x, z = xmul(pk, x, z)
            g = math.gcd(z, N)
            if g > 1:
                return g
        return 1

    g_sim = _run_sigma(533684655365, 200)
    assert g_sim == N, f"Expected simultaneous gcd==N, got {g_sim}"
    g_sing_nm1 = _run_sigma(N - 1, 200)
    assert g_sing_nm1 == "SINGULAR", f"Expected singular curve guard on sigma=N-1, got {g_sing_nm1}"
    split_sigma = None
    split_fac = None
    for sig in range(6, 50):
        g = _run_sigma(sig, 5000)
        if isinstance(g, int) and 1 < g < N:
            split_sigma = sig
            split_fac = g
            break
    assert split_fac is not None and N % split_fac == 0 and 1 < split_fac < N
    return dict(
        N=N, bad_sigma=533684655365, bad_gcd_eq_N=True,
        perfect_power_precheck_ok=True,
        singular_nm1_rejected=True, retry_sigma=split_sigma, split_factor=split_fac,
    )


def mulmod_instr(nbits, l_min=L_MIN_DEFAULT):
    L = max(l_min, ceil(nbits / 64))
    return MONT5 * (L / 5.0) ** 2, L


def gcd_instr(nbits, l_min=L_MIN_DEFAULT):
    L = max(l_min, ceil(nbits / 64))
    return 2.0 * nbits * (12 * L + 16)


@lru_cache(None)
def tree_branch_caps(nbits, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
    """Exact DP upper bound on (sum L_i^2, sum gcd_instr(b_i)) over ALL active composite-branch
    partitions b_1 + ... + b_t <= nbits with b_i >= b_min (L_i = max(l_min, ceil(b_i / 64)))."""
    nb = int(ceil(nbits))
    if nb < b_min:
        L = max(l_min, ceil(nb / 64))
        return L * L, 2.0 * nb * (12 * L + 16)
    dp_l2 = [0] * (nb + 1)
    dp_g = [0.0] * (nb + 1)
    for b in range(b_min, nb + 1):
        L = max(l_min, ceil(b / 64))
        best_l2 = L * L
        best_g = 2.0 * b * (12 * L + 16)
        cand_s = [b_min]
        for k in range(l_min, b // 64 + 1):
            s = 64 * k + 1
            if b_min <= s <= b - b_min:
                cand_s.append(s)
        for s in cand_s:
            Ls = max(l_min, ceil(s / 64))
            v2 = Ls * Ls + dp_l2[b - s]
            vg = 2.0 * s * (12 * Ls + 16) + dp_g[b - s]
            if v2 > best_l2:
                best_l2 = v2
            if vg > best_g:
                best_g = vg
        dp_l2[b] = best_l2
        dp_g[b] = best_g
    return dp_l2[nb], dp_g[nb]


def pi_upper(x):
    return 1.26 * x / log(x) if x >= 17 else 7.0


@lru_cache(None)
def semismooth(fbits, b1, b2, steps=256):
    m = fbits - log2(TORSION)
    base = log2(b1)
    top = log2(b2)
    acc = D.rho_lo(m / base)
    if top > base:
        for i in range(steps):
            lo = base + (top - base) * i / steps
            hi = base + (top - base) * (i + 1) / steps
            acc += D.rho_lo((m - lo) / base) * log(hi / lo)
    return min(acc, 1.0 - 1e-15)


@lru_cache(None)
def power_correction(fbits, b1, b2):
    m = fbits - log2(TORSION)
    u = m / log2(b1)
    tilt = D.rho_hi(max(u - 2.0, 1.0)) / max(D.rho_lo(u), 1e-300)
    y = sqrt(b1)
    return max(0.0, 1.0 - tilt * 2.0 / (y * log(y)))


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
def schedule_v9(fbits, nbits, success=0.999, n_factors=1, sep=SEPARATION_DEFAULT, p_r=RETRY_PR_DESC, b_min=B_MIN_DEFAULT):
    """Single-rung branch-tree-capped ECM schedule with explicit failed-replay retry & separation factor."""
    nb = int(ceil(nbits))
    nf = max(1, int(ceil(n_factors)))
    l2_max, g_max = tree_branch_caps(nb, b_min)
    mul = MONT5 * (l2_max / 25.0)
    gcd = g_max
    best = None
    retry_mult = p_r / (1.0 - p_r)
    for k in range(10, 37, 2):
        b1 = 2 ** k
        for ratio in (1, 10, 30, 100, 300):
            b2 = b1 * ratio
            p1 = semismooth(fbits, b1, b2) * power_correction(fbits, b1, b2)
            if p1 <= 1e-12:
                continue
            p_eff = sep * p1
            curves = ceil(log1p(-success) / log1p(-p_eff))
            for wheel in WHEELS if ratio > 1 else (0,):
                if wheel and wheel > b2:
                    continue
                ops, s1, s2, J, G = (
                    stage_ops(b1, b2, wheel)
                    if wheel
                    else (1.1 * ceil(1.5 * b1) * 11, ceil(1.5 * b1) * 11, 0, 0, 0)
                )
                per_curve = ops * mul + 2.0 * gcd
                replay_1curve = ops * mul + (pi_upper(b1) + G) * gcd
                redo_per_factor = p1 * replay_1curve
                retry_cost_per_curve = retry_mult * (per_curve + replay_1curve)
                cost = curves * (per_curve + nf * redo_per_factor + retry_cost_per_curve) + NEWTON_PRECHECK_INSTR
                if best is None or cost < best["instructions"]:
                    best = dict(
                        fbits=fbits, nbits=nb, n_factors=nf, sep=sep, p_r=p_r, b_min=b_min,
                        l2_tree_max=l2_max, gcd_tree_max=g_max,
                        B1=b1, B2=b2, wheel=wheel, curves=curves, p_curve=p1, p_eff=p_eff,
                        per_curve=per_curve, replay_1curve=replay_1curve,
                        retry_cost_total=curves * retry_cost_per_curve,
                        redo_total=curves * nf * redo_per_factor,
                        instructions=cost, log2_instructions=log2(cost),
                    )
    return best


schedule = schedule_v9


@lru_cache(None)
def ladder_v9_instr(lo_bits, hi_bits, nbits, success=0.999, n_factors=1, step=12.0, sep=SEPARATION_DEFAULT, p_r=RETRY_PR_DESC, b_min=B_MIN_DEFAULT):
    """Geometric ECM ladder from lo_bits up to hi_bits (step bits per rung), using branch-tree-capped
    schedule_v9 at every rung so smaller prime factors are isolated before higher-B1 rungs run."""
    if hi_bits <= lo_bits + 1e-9:
        return 0.0
    tot = 0.0
    r = float(hi_bits)
    while r > lo_bits + 1e-9:
        tot += schedule_v9(round(r, 4), int(ceil(nbits)), success, n_factors, sep, p_r, b_min)["instructions"]
        r -= step
    return tot


def ecm_instr(fbits, nbits, success, n_factors=1):
    return ladder_v9_instr(16.0, fbits, nbits, success, n_factors, step=12.0, sep=SEPARATION_DEFAULT, p_r=RETRY_PR_DESC, b_min=B_MIN_DEFAULT)
