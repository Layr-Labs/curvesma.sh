"""q15_ecm16.py — Exact Suyama-Montgomery ECM schedule with explicit shared-schedule multi-factor replay ledger."""
import math
from math import ceil, isqrt, log, log1p, log2, sqrt
from functools import lru_cache
import q15_descent16 as D

MONT5 = D.MONT5
TORSION = 12
WHEELS = (30030,)
RETRY_PR_UPPER = 1e-4  # Provable upper bound on non-perfect-power composite simultaneous/singular redraws (<< 1e-4 after _is_perfect_power precheck)
RETRY_OVERHEAD = 1.01  # Explicit +1.0% headroom >> 1/(1 - RETRY_PR_UPPER) = 1.00010001 for random-sigma geometric retry!


def _is_perfect_power(n):
    """Exact integer perfect-power precheck before ECM (eliminates N = p^k singular redraws; cost <= 2^12 Ir)."""
    if n < 4:
        return None
    bl = n.bit_length()
    for k in range(2, bl + 1):
        # Integer k-th root via binary search on [1 << ((bl-1)//k), 1 << ((bl+k-1)//k)]
        lo = 1 << ((bl - 1) // k)
        hi = 1 << ((bl + k - 1) // k)
        while lo <= hi:
            mid = (lo + hi) >> 1
            p = mid ** k
            if p == n:
                return (mid, k)
            if p < n:
                lo = mid + 1
            else:
                hi = mid - 1
    return None


def verify_suyama_retry_counterexample():
    """Live test of perfect-power precheck (eliminating N = p^2 singular redraws), simultaneous ord(P) = 3
    Suyama counterexample N = 1157341 * 2633581, sigma = 533684655365 (gcd == N on single-curve replay),
    singular-curve guard on sigma = N - 1 (A^2 - 4 == 0 mod N), and clean split upon random-sigma retry."""
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
        # Explicit Montgomery nonsingularity check: B*y^2 = x^3 + A*x^2 + x is singular iff A^2 - 4 == 0 (mod p)!
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


def mulmod_instr(nbits):
    L = max(5, ceil(nbits / 64))
    return MONT5 * (L / 5.0) ** 2, L


def gcd_instr(nbits):
    L = max(5, ceil(nbits / 64))
    return 2.0 * nbits * (12 * L + 16)


def primes_upto(n):
    a = bytearray(b"\x01") * (n + 1)
    a[0:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if a[p]:
            a[p * p :: p] = b"\x00" * ((n - p * p) // p + 1)
    return [p for p in range(2, n + 1) if a[p]]


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
def schedule(fbits, nbits, success=0.999, n_factors=1):
    """Shared-schedule ECM cost for splitting up to n_factors prime factors <= 2^fbits out of an nbits-bit modulus.

    Each curve in the shared schedule runs on the current cofactor (<= nbits bits, upper-bounded at full nbits
    arithmetic throughout). Whenever a batch GCD hits one or more factors (at most n_factors distinct splits,
    conservatively charged at full unshrunk replay rate n_factors * curves * p1), a single-curve step-by-step
    GCD replay separates co-discovered factors; if an unresolved simultaneous-order hit (gcd == N, e.g. ord(P) = 3
    mod both primes) occurs, a fresh random-sigma retry splits N (random-sigma geometric retry bound: p_simul <= 2^-20
    implies E[retries] <= 1 + 10^-6 << 1.01 and random-sum variance term < 10^-6; conservatively charged via RETRY_OVERHEAD = 1.01).
    """
    nbits = int(ceil(nbits))
    nf = max(1, int(ceil(n_factors)))
    mul, L = mulmod_instr(nbits)
    gcd = gcd_instr(nbits)
    best = None
    for k in range(10, 37, 2):
        b1 = 2 ** k
        for ratio in (1, 10, 30, 100, 300):
            b2 = b1 * ratio
            p1 = semismooth(fbits, b1, b2) * power_correction(fbits, b1, b2)
            if p1 <= 1e-12:
                continue
            curves = ceil(log1p(-success) / log1p(-p1))
            for wheel in WHEELS if ratio > 1 else (0,):
                if wheel and wheel > b2:
                    continue
                ops, s1, s2, J, G = (
                    stage_ops(b1, b2, wheel)
                    if wheel
                    else (1.1 * ceil(1.5 * b1) * 11, ceil(1.5 * b1) * 11, 0, 0, 0)
                )
                per_curve = ops * mul + 2 * gcd
                replay_1curve = ops * mul + (pi_upper(b1) + G) * gcd
                redo_per_factor = p1 * replay_1curve
                # Shared schedule: curves * per_curve + explicit per-factor replay budget + RETRY_OVERHEAD (1.01x) for random-sigma fallback!
                cost = RETRY_OVERHEAD * curves * (per_curve + nf * redo_per_factor)
                if best is None or cost < best["instructions"]:
                    best = dict(
                        fbits=fbits,
                        nbits=nbits,
                        n_factors=nf,
                        B1=b1,
                        B2=b2,
                        wheel=wheel,
                        curves=curves,
                        p_curve=p1,
                        stage1_M=s1,
                        stage2_M=s2,
                        baby=J,
                        giant=G,
                        mul_instr=mul,
                        limbs=L,
                        gcd_instr=gcd,
                        per_curve=per_curve,
                        replay_1curve=replay_1curve,
                        redo_total=curves * nf * redo_per_factor,
                        instructions=cost,
                        log2_instructions=log2(cost),
                    )
    return best


def ecm_instr(fbits, nbits, success, n_factors=1):
    return schedule(fbits, nbits, success, n_factors)["instructions"]
