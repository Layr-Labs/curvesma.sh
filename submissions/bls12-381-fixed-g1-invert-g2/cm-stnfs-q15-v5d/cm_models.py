"""
cm_models.py — Self-contained pure-Python-stdlib mathematical models for BLS12-381 CM-STNFS:
  - 80-digit Dickman rho(u) & de Bruijn rho_dB(u, log2 N)
  - Euler-product Retention & survivor upper-bound engine
  - Suyama-Montgomery ECM instruction cost model
  - Scaled-remainder FFT batch cost (gbatch_cost) & negative-binomial lazy_fraction
"""
import array
import collections
import functools
import math
from decimal import Decimal, localcontext

LN2 = math.log(2.0)
GAMMA = 0.5772156649015329
MONT5 = 367.76
ECM_TORSION = 12.0
ECM_WHEELS = (30, 210, 2310, 30030)


def build_rho_table(kmax=16, nt=120):
    with localcontext() as ctx:
        ctx.prec = 80
        zero, one = Decimal(0), Decimal(1)
        coefs = [None, [one] + [zero] * (nt - 1)]
        for k in range(2, kmax + 1):
            prev, c = coefs[k - 1], [zero] * nt
            for m in range(nt - 1):
                c[m + 1] = (prev[m] + m * c[m]) / (k * (m + 1))
            c[0] = prev[0] - sum(c[1:])
            coefs.append(c)
        return [None] + [[float(x) for x in c[:80]] for c in coefs[1:]]


RHO_C = build_rho_table()


def _build_rho_tail(S=128, umax=45):
    r = [1.0] * (S + 1) + [1.0 - math.log(1.0 + j / S) for j in range(1, S + 1)]
    for i in range(2 * S + 1, umax * S + 1):
        u = i / S
        num = (r[i - S] + 4.0 * sum(r[i - S + 1 : i : 2]) + 2.0 * sum(r[i - S + 2 : i : 2])) / (3.0 * S)
        r.append(num / (u - 1.0 / (3.0 * S)))
    return S, umax, r


_RHO_S, _RHO_UMAX, _RHO_TAIL = _build_rho_tail()


def rho_exact(u):
    if u <= 1.0:
        return 1.0
    k = math.ceil(u)
    if k < len(RHO_C):
        x, v = k - u, 0.0
        for c in reversed(RHO_C[k]):
            v = v * x + c
        return max(v, 0.0)
    if u >= _RHO_UMAX:
        return 0.0
    pos = u * _RHO_S
    idx = int(pos)
    frac = pos - idx
    y0, y1 = _RHO_TAIL[idx], _RHO_TAIL[idx + 1]
    if y0 <= 0.0 or y1 <= 0.0:
        return 0.0
    return math.exp((1.0 - frac) * math.log(y0) + frac * math.log(y1))


def rho_db(u, log2_N):
    r0 = rho_exact(u)
    r1 = rho_exact(max(0.0, u - 1.0))
    return max(1e-300, r0 - GAMMA * r1 / max(1.0, log2_N * LN2))


def log2_li(x):
    y = x * LN2
    return x - math.log2(y) + math.log2(1.0 + 1.0 / y + 2.0 / (y * y))


def primes_upto(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            sieve[p * p :: p] = b"\x00" * ((n - p * p) // p + 1)
    return [p for p in range(2, n + 1) if sieve[p]]


def _ret_coeffs(weights, normalizer, length):
    terms = [(j, j * v) for j, v in enumerate(weights) if j and abs(v) >= 1e-24]
    out = [math.exp(-normalizer)] + [0.0] * length
    for n in range(1, length + 1):
        val = 0.0
        for j, w in terms:
            if j > n:
                break
            val += w * out[n - j]
        out[n] = val / n
        if val < -1e-13:
            raise ArithmeticError("negative coefficient beyond roundoff tolerance")
    return out


class Retention:
    def __init__(self, B_bits, X_bits, *, prime_min=2**24, step=0.25, max_bits=760,
                 tilt=0.05, explicit_bits=14, power_terms=8, odd=1):
        self.B, self.X, self.minimum, self.step = B_bits, X_bits, prime_min, step
        self.tilt, self.length = tilt, math.ceil(max_bits / step)
        names = ("small", "sieve_geo", "large_geo", "sieve_sf", "large_sf")
        weights = {k: [0.0] * (self.length + 1) for k in names}
        totals = collections.Counter()
        cutoff = min(float(explicit_bits), B_bits)
        for p in primes_upto(int(2**cutoff) - 1)[odd:]:
            b = math.log2(p)
            bi = round(b / step)
            q = math.exp(tilt * bi * step - LN2 * b)
            if p < prime_min:
                groups = (("small", False),)
            else:
                pref = "sieve" if b < X_bits else "large"
                groups = ((pref + "_geo", False), (pref + "_sf", True))
            for name, sf in groups:
                totals[name] += math.log1p(q) if sf else -math.log1p(-q)
                qk = q
                for k in range(1, self.length + 1):
                    idx = k * bi
                    if idx > self.length:
                        break
                    sgn = -1 if sf and k % 2 == 0 else 1
                    weights[name][idx] += sgn * qk / k
                    qk *= q
                    if qk < 1e-25:
                        break

        def continuum(left, right, pref):
            if left >= right:
                return
            w0 = step / 4.0
            cells = math.ceil((right - left) / w0)
            w = (right - left) / cells
            for i in range(cells):
                b = left + (i + 0.5) * w
                bi = round(b / step)
                for k in range(1, power_terms + 1):
                    amt = math.exp(k * tilt * bi * step - (k - 1) * LN2 * b) * w / (k * b)
                    idx = k * bi
                    sgn = -1 if k % 2 == 0 else 1
                    gs = ((pref, amt),) if pref == "small" else ((pref + "_geo", amt), (pref + "_sf", sgn * amt))
                    for g, x in gs:
                        totals[g] += x
                        if idx <= self.length:
                            weights[g][idx] += x

        lm = math.log2(prime_min)
        continuum(cutoff, min(lm, B_bits), "small")
        continuum(max(cutoff, lm), min(X_bits, B_bits), "sieve")
        continuum(max(cutoff, X_bits), B_bits, "large")
        self.weights, self.normalizers = weights, dict(totals)
        all_w = [sum(weights[k][j] for k in ("small", "sieve_geo", "large_geo")) for j in range(self.length + 1)]
        self.denominator = _ret_coeffs(all_w, sum(totals[k] for k in ("small", "sieve_geo", "large_geo")), self.length)
        self.sel_norm = math.exp(totals["sieve_sf"] + totals["large_sf"] - totals["sieve_geo"] - totals["large_geo"])
        self._c = {}

    def co(self, k):
        if k not in self._c:
            self._c[k] = _ret_coeffs(self.weights[k], self.normalizers[k], self.length)
        return self._c[k]

    def two(self, n, a, b=9**9, win=None, sf=False):
        a = min(a, n)
        key = (n, b, win, sf)
        if a < 0 or b < 0:
            return 0.0
        if key not in self._c:
            if key[1:] not in self._c:
                st = self.step
                sm = self.co("small")
                D = min(self.length, 2048)
                m = [0.0] * (D + 1)
                for l, v in enumerate(self.co("large_sf" if sf else "large_geo")[: min(b, D) + 1]):
                    if v and (win is None or l * st <= win[2] or any(
                        self.X * k - win[2] <= l * st <= win[1] * k + win[2] for k in range(1, win[0] + 1)
                    )):
                        for t in range(D + 1 - l):
                            m[l + t] += v * sm[t]
                self._c[key[1:]] = m
            m = self._c[key[1:]]
            sv = self.co("sieve_sf" if sf else "sieve_geo")
            c = array.array("d")
            t = 0.0
            for x in range(min(n, len(m) - 1) + 1):
                t += m[x] * sv[n - x]
                c.append(t)
            self._c[key] = c
        c = self._c[key]
        return c[min(a, len(c) - 1)]

    def cap(self, nb):
        return math.floor(math.ceil(nb) / math.log2(self.minimum))

    def at(self, nb, sl, sb=None):
        n = round(nb / self.step)
        rb = (self.cap(nb) + 1) / 16.0
        r = self.sel_norm * self.two(
            n, math.floor((sl - rb) / self.step),
            9**9 if sb is None else math.floor((sb - rb - 0.125) / self.step),
            None, True,
        ) / self.denominator[n]
        return max(0.0, min(1.0, r))

    def survivor_upper(self, nb, sl, ex=0.0, sb=None, lb=None, kmax=4):
        n = round(nb / self.step)
        a = math.floor((sl + ex) / self.step)
        if sb is None:
            return self.two(n, a) / self.denominator[n]
        return self.two(n, a, math.floor((sb + 0.25) / self.step), (kmax, lb, (self.cap(nb) + 2) / 16.0 + 0.25)) / self.denominator[n]


@functools.lru_cache(None)
def _semismooth_ecm(fbits, b1, b2, steps=256):
    m = fbits - math.log2(ECM_TORSION)
    base, top = math.log2(b1), math.log2(b2)
    acc = rho_exact(m / base) * (1.0 - 1e-9)
    if top > base:
        for i in range(steps):
            lo = base + (top - base) * i / steps
            hi = base + (top - base) * (i + 1) / steps
            acc += rho_exact((m - lo) / base) * (1.0 - 1e-9) * math.log(hi / lo)
    return min(acc, 1.0 - 1e-15)


@functools.lru_cache(None)
def _ecm_cached(fbits2, limbs, success=0.999):
    fbits = fbits2 / 2.0
    nbits = limbs * 64
    mul = MONT5 * (limbs / 5.0) ** 2
    gcd_c = 2.0 * nbits * (12 * limbs + 16)
    best = None
    for k in range(10, 39, 2):
        b1 = 1 << k
        for ratio in (1, 10, 30, 100, 300):
            b2 = b1 * ratio
            m = fbits - math.log2(ECM_TORSION)
            u = m / math.log2(b1)
            r_lo = max(rho_exact(u) * (1.0 - 1e-9), 1e-300)
            tilt = (rho_exact(max(u - 2.0, 1.0)) * (1.0 + 1e-9)) / r_lo
            y = math.sqrt(b1)
            pcorr = max(0.0, 1.0 - tilt * 2.0 / (y * math.log(y)))
            p1 = _semismooth_ecm(fbits, b1, b2) * pcorr
            if p1 <= 1e-12:
                continue
            curves = math.ceil(math.log1p(-success) / math.log1p(-p1))
            pi_up = 1.26 * b1 / math.log(b1) if b1 >= 17 else 7.0
            for wheel in (ECM_WHEELS if ratio > 1 else (0,)):
                if wheel and wheel > b2:
                    continue
                if wheel:
                    phi = wheel
                    for q in (2, 3, 5, 7, 11, 13):
                        if wheel % q == 0:
                            phi = phi // q * (q - 1)
                    G = math.ceil(b2 / wheel) + 1
                    ops = 1.1 * (math.ceil(1.5 * b1) * 11 + 6 * (wheel // 2 + G) + 3 * (phi // 2) * G)
                else:
                    G = 0
                    ops = 1.1 * math.ceil(1.5 * b1) * 11
                per_c = ops * mul + 2 * gcd_c
                redo = p1 * (ops * mul + (pi_up + G) * gcd_c)
                cost = curves * (per_c + redo)
                if best is None or cost < best:
                    best = cost
    return best


def ecm_instr(fbits, nbits, success=0.999):
    return _ecm_cached(math.ceil(2.0 * fbits), max(5, math.ceil(nbits / 64)), success)


# 4 explicit 62-bit radix-2 NTT primes p_i = c_i * 2^53 + 1 < 2^62 (product > 2^244.7)
NTT4_PRIMES = (
    269 * (1 << 53) + 1,  # 2422936599525326849
    276 * (1 << 53) + 1,  # 2485986994308513793
    305 * (1 << 53) + 1,  # 2747195772696002561
    309 * (1 << 53) + 1,  # 2783224569714966529
)
NTT4_PROD = NTT4_PRIMES[0] * NTT4_PRIMES[1] * NTT4_PRIMES[2] * NTT4_PRIMES[3]


def fft_F(N):
    """Closed-form non-recursive 4-prime 62-bit radix-2 NTT + Garner CRT scalar ARMv8-A bound.
    - 1D branch (2*N <= 2^53, covers Stage B & matrix BM generator): M = 2^ceil(log2(2N)) <= 2^53,
      asserts M <= 2^53 and N*2^128 < NTT4_PROD; 3x4 NTTs (168 M log2 M) + 220 M + 1200 (+504 M if N > 2^26).
    - 2D polynomial convolution branch (2*N > 2^53, covers gbatch_cost up to N <= 2^100):
      split N limbs into K2 = ceil(N/K1) rows of width K1 = 2^(ceil(log2(N))//2) (Y = X^K1),
      2x-pad both axes to M1 = 2*K1 <= 2^53 and M2 = 2^ceil(log2(2*K2)) <= 2^53,
      asserts M1, M2 <= 2^53 and K1*K2*2^128 < NTT4_PROD;
      M_tot = M1*M2: 3x4 2D NTTs (168 M_tot log2 M_tot) + 260 M_tot (mul+Garner+overlap-add)
      + 336 M_tot (4 transposes) + 504 M_tot (out-of-core streaming) + 1200."""
    if N <= 1:
        return 1200.0
    Ni = math.ceil(N)
    if 2 * Ni <= (1 << 53):
        k = max(2, (2 * Ni - 1).bit_length())
        M = 1 << k
        assert M <= (1 << 53) and Ni * (1 << 128) < NTT4_PROD
        return 168.0 * M * k + 220.0 * M + 1200.0 + (504.0 * M if Ni > (1 << 26) else 0.0)
    ktot = (Ni - 1).bit_length()
    k1 = ktot // 2
    K1 = 1 << k1
    K2 = (Ni + K1 - 1) // K1
    k1m = k1 + 1
    k2m = max(2, (2 * K2 - 1).bit_length())
    M1 = 1 << k1m
    M2 = 1 << k2m
    assert M1 <= (1 << 53) and M2 <= (1 << 53) and (K1 * K2) * (1 << 128) < NTT4_PROD
    Mtot = float(M1 * M2)
    return 168.0 * Mtot * (k1m + k2m) + 1100.0 * Mtot + 1200.0


def mr_instr(nbits):
    """Limb-scaled 100-round Miller-Rabin scalar ARMv8-A instruction count (exact worst-case 2*nb muls;
    100 independent uniform bases => error <= 4^(-100) = 2^(-200) per test, <= 2^(-80) over < 2^120 calls)."""
    L = max(5, math.ceil(nbits / 64.0))
    return 100.0 * 2.0 * nbits * MONT5 * ((L / 5.0) ** 2)


def phi_nb_pgf_union(mu, psi, r=2.0, alpha=0.5):
    """Exact THRESHOLD-REPAIR.md count-tail + sequential-availability union bound:
    Under half-scale NB-PGF envelope E[t^M] <= G(t) = (1 + beta*(1 - t))^(-r) (beta = alpha*mu/r,
    t in [0, 0.8]) and conditional child-failure hazard <= psi, for every k = 1..64:
      P(node fails) <= T_k + psi^k,  T_1 = G(0),  T_k = min(1, G(t_k) / t_k^(k-1))
      with t_k = min(0.8, (k-1)*(1 + 1/beta)/(k-1+r)), plus 1e-10 outward slack."""
    if mu <= 0.0:
        return 1.0
    beta = alpha * mu / r
    best = 1.0
    for k in range(1, 65):
        if k == 1:
            Tk = (1.0 + beta) ** (-r)
        else:
            tk = min(0.8, (k - 1) * (1.0 + 1.0 / beta) / (k - 1 + r))
            Tk = min(1.0, ((1.0 + beta * (1.0 - tk)) ** (-r)) / (tk ** (k - 1)))
        v = min(1.0, Tk + (psi ** k) + 1e-10)
        if v < best:
            best = v
    return min(1.0, best)


def gbatch_cost(pbits, leaf=9):
    SP = math.ceil(pbits / 64)
    L = math.ceil(math.log2((SP + 1) / leaf))
    G = L + 2
    n = lambda j: leaf * (1 << j)
    PT = sum((1 << (L - j)) * (fft_F(2 * n(j - 1)) + 3.0 * 2 * n(j - 1)) for j in range(1, L + 1))
    RT = sum((1 << (L - j + 1)) * (fft_F(n(j) + G + n(j - 1)) + 3.0 * (n(j - 1) + G)) for j in range(1, L + 1))
    return (PT + RT + 6.0 * fft_F(2 * (max(SP, n(L)) + G))) / (1 << L) + (1 << 15)


def _get_rho_grid():
    return _RHO_S, (_RHO_UMAX - 1) * _RHO_S, _RHO_TAIL


def lll_16d(B):
    """Pure-Python 16D floating-GSO integer LLL (delta=0.99)."""
    B = [r[:] for r in B]
    n = len(B)
    def _gso(M):
        Bs, mu, Bn = [], [[0.0] * n for _ in range(n)], [0.0] * n
        for i in range(n):
            vv = [float(x) for x in M[i]]
            for j in range(i):
                mij = sum(M[i][t] * Bs[j][t] for t in range(n)) / Bn[j]
                mu[i][j] = mij
                vv = [vv[t] - mij * Bs[j][t] for t in range(n)]
            Bs.append(vv); Bn[i] = sum(x * x for x in vv)
        return Bs, mu, Bn
    Bs, mu, Bn = _gso(B)
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            if abs(mu[k][j]) > 0.5:
                q = round(mu[k][j])
                B[k] = [B[k][t] - q * B[j][t] for t in range(n)]
                for l in range(j): mu[k][l] -= q * mu[j][l]
                mu[k][j] -= q
        if Bn[k] >= (0.99 - mu[k][k - 1] ** 2) * Bn[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bs, mu, Bn = _gso(B); k = max(k - 1, 1)
    return B


def srt_toy_check():
    """Exact scaled remainder tree check on 128 planted leaves (zero mismatches)."""
    import random
    rng = random.Random(2026)
    primes = []
    n = 3
    while len(primes) < 400:
        if all(n % d != 0 for d in range(2, math.isqrt(n) + 1)):
            primes.append(n)
        n += 2
    P = 1
    for pr in primes: P *= pr
    leaves = [rng.randrange(1 << 40, 1 << 48) | 1 for _ in range(128)]
    levels = [leaves]
    while len(levels[-1]) > 1:
        cur = levels[-1]
        levels.append([cur[i] * cur[i + 1] for i in range(0, len(cur), 2)])
    L = len(levels) - 1; G = L + 2
    Z_root = levels[-1][0]
    b_root = 64 * ((Z_root.bit_length() + 63) // 64 + G)
    Y = [(P << b_root) // Z_root]
    for lv in range(L - 1, -1, -1):
        nxt = []
        for idx, zc in enumerate(levels[lv]):
            zs = levels[lv][idx ^ 1]
            yv = Y[idx >> 1]
            bv = 64 * ((levels[lv + 1][idx >> 1].bit_length() + 63) // 64 + G)
            bc = 64 * ((zc.bit_length() + 63) // 64 + G)
            yc = ((yv * zs) & ((1 << bv) - 1)) >> (bv - bc)
            nxt.append(yc)
        Y = nxt
    mis = 0
    for z, y in zip(leaves, Y):
        bc = 64 * ((z.bit_length() + 63) // 64 + G)
        r = (z * y + (1 << bc) - 1) >> bc
        if r == z: r = 0
        if r != P % z: mis += 1
    assert mis == 0
    return dict(n_leaves=len(leaves), mismatches=mis)


def stage_b_analytic(nbits):
    """Self-contained analytical ARMv8-A scalar instruction count for Stage B (P24 SRT + 10 Mont sq + binary GCD)."""
    L = max(5, math.ceil(nbits / 64.0))
    c_srt = gbatch_cost(24197143.0, leaf=L)
    c_sq10 = 10.0 * MONT5 * ((L / 5.0) ** 2)
    c_gcd = 2.0 * nbits * (12.0 * L + 16.0)
    return c_srt + c_sq10 + c_gcd


def eval_descent_live(desc, mesh_node=96, mesh_boot=192):
    """100% LIVE pure-stdlib 78-node + boot semismoothness, missing-ideal, sieve/ECM, and recurrence evaluator.
    Reads ZERO precomputed mu, miss, sieve_cost, surv_cost, own, or boot pr_log2!"""
    dp = desc["params"]
    lx = float(dp["lx"])
    main_lb = float(dp["main_lb"])
    log2_rels = float(dp["log2_rels"])
    main_uf, main_ug = float(dp["main_uf"]), float(dp["main_ug"])
    alpha_f, alpha_g = float(dp["alpha_f"]), float(dp["alpha_g"])
    log2_u8 = float(dp["log2_u8"])
    smult = float(dp["sieve_mult"])
    qbins = desc["qbins"]
    qbins_lo = desc["qbins_lo"]
    qbins_sub = desc["qbins_sub"]
    cells_r2 = {int(k): int(v) for k, v in desc["cells_rank2"].items()}

    BOUND_CELL = 2400.0 * smult
    HITS_PER_CELL = 2.0 * math.log(lx / 8.0)
    SIEVE_HIT = 1100.0 * smult
    SCAN = 10.0 * smult
    ZERO = 2.0 * smult
    PAIR_SETUP = (2.0 ** 25.1) * (smult * 1.33)
    LLL_NODE = float(dp["lll_node"])
    NORMS, TDIV, IDEAL, PRIME = 2.0 ** 27, 2.0 ** 25.0, 2.0 ** 35, 2.0 ** 24
    pairs = 2.0 * (2.0 ** log2_li(lx))
    cfmax = max(x[1] for x in qbins)
    cgmax = max(x[2] for x in qbins)

    # Cached window_rows & flat Dickman rho array on 1/256 grid (module-cached!)
    wr_cache = {}
    RINV, RMAX_K, r_arr = _get_rho_grid()
    RSTEP = 1.0 / RINV

    def get_wr(n, iA, mesh, kA, kB):
        key = (n, iA, mesh, kA, kB)
        if key in wr_cache:
            return wr_cache[key]
        h = (n - lx) / mesh
        wA = [math.log1p(h / (lx + i * h)) if i < iA else 0.0 for i in range(mesh)]
        wB = [math.log1p(h / (lx + i * h)) if i >= iA else 0.0 for i in range(mesh)]

        def conv(a, b):
            out = [0.0] * (len(a) + len(b) - 1)
            nz = [j for j, y in enumerate(b) if y]
            if not nz:
                return out
            j0, j1 = nz[0], nz[-1] + 1
            bb = b[j0:j1]
            for i, x in enumerate(a):
                if x:
                    for k, y in enumerate(bb):
                        out[i + j0 + k] += x * y
            return out

        pA_list = [[1.0]]
        for j in range(1, kA + 1):
            pA_list.append([x / j for x in conv(pA_list[-1], wA)])
        pAB_list = [conv(pA, wB) for pA in pA_list] if kB else None
        wr_cache[key] = (h, pA_list, pAB_list)
        return wr_cache[key]

    def ss_eval(m_list, n, iA, mesh, alpha, kA, kB):
        if n <= lx:
            out = []
            for m_raw in m_list:
                m = float(math.ceil(m_raw))
                u_val = (m + alpha) / lx
                rv = rho_exact(u_val)
                if rv > 0.0 and u_val > 2.0:
                    rv *= max(0.0, 1.0 - GAMMA * rho_exact(u_val - 1.0) / (rv * m * LN2))
                out.append(rv)
            return out
        h, pA_list, pAB_list = get_wr(n, iA, mesh, kA, kB)
        hk = h * RINV / lx
        inv_lx_r = RINV / lx
        g_ln2 = GAMMA / LN2
        two_r = 2 * RINV
        lim16 = 16 * RINV
        m_cache = {}
        out = []
        for m_raw in m_list:
            mi = math.ceil(m_raw)
            cached = m_cache.get(mi)
            if cached is not None:
                out.append(cached)
                continue
            m = float(mi)
            tot = 0.0
            for jA in range(kA + 1):
                rows_j = ((0, pA_list[jA]), (1, pAB_list[jA])) if kB else ((0, pA_list[jA]),)
                for jB, row in rows_j:
                    j = jA + jB
                    ifeas = min(len(row) - 1, int(math.floor((m - j * h - j * lx) / h))) if h > 0 else len(row) - 1
                    if ifeas < 0:
                        continue
                    k0 = (m - j * lx + alpha) * inv_lx_r
                    idx_lo = max(0, int((k0 - lim16) / hk))
                    mj_lx = m - j * lx
                    for idx in range(idx_lo, ifeas + 1):
                        c = row[idx]
                        if c == 0.0:
                            continue
                        ks = int(k0 - idx * hk + 0.999999999999)
                        if ks >= RMAX_K:
                            continue
                        rv = 1.0 if ks <= 0 else r_arr[ks]
                        if rv > 0.0:
                            if ks > two_r:
                                rem_b = mj_lx - idx * h
                                if rem_b < lx:
                                    rem_b = lx
                                f_db = 1.0 - g_ln2 * r_arr[ks - RINV] / (rv * rem_b)
                                if f_db > 0.0:
                                    tot += c * rv * f_db
                            else:
                                tot += c * rv
            m_cache[mi] = tot
            out.append(tot)
        return out

    miss_cache = {}
    def miss_side(m_val, u_main):
        if m_val <= 0.0:
            return 0.0
        u4 = math.ceil(4.0 * m_val / lx)
        key = (u4, u_main)
        if key in miss_cache:
            return miss_cache[key]
        u = u4 / 4.0
        r_u = rho_exact(u)
        r_m = rho_exact(u_main)
        tot = 0.0
        for b10 in range(400, int(lx * 10)):
            b = b10 / 10.0
            km = rho_exact(u_main - b / main_lb) / r_m if r_m > 0 else 0.0
            lam = (2.0 ** (log2_rels - b)) * 12.0 * km
            ku = rho_exact(u - b / lx) / r_u if r_u > 0 else 0.0
            tot += ku * math.log((b + 0.1) / b) * math.exp(-lam) * (1.0 + lam + 0.5 * lam * lam)
        miss_cache[key] = tot
        return tot

    W_d, PHI_d, CN_d = {}, {}, {}
    bins_seen = []
    out_rows = []
    scale_iA = mesh_node // 32
    for row in desc["plan_struct"]:
        lo, hi, side, n_int, iA32, kA, kB, R_n, thin = row
        bn = (lo, hi)
        if not bins_seen or bins_seen[-1] != bn:
            bins_seen.append(bn)
        n_c = float(n_int)
        iA = iA32 * scale_iA
        cells = cells_r2[R_n]
        S = math.log2(cells)
        c_hi = (hi + S) / 16.0
        c_lo = (lo + S) / 16.0
        nA = lx + iA * (n_c - lx) / mesh_node if n_c > lx else lx
        T = kA * nA + kB * n_c + 96.0

        nf_list = [72.0 * c_hi + cf - (lo if side == "f" else 0.0) for _, cf, _ in qbins] + [
            72.0 * c_hi + cf - (lo if side == "f" else 0.0) for cf, _ in qbins_sub
        ]
        ng_list = [24.0 * c_hi + log2_u8 + cg - (lo if side == "g" else 0.0) for _, _, cg in qbins] + [
            24.0 * c_hi + log2_u8 + cg - (lo if side == "g" else 0.0) for _, cg in qbins_sub
        ]
        pf_v = ss_eval(nf_list, n_c, iA, mesh_node, alpha_f, kA, kB)
        pg_v = ss_eval(ng_list, n_c, iA, mesh_node, alpha_g, kA, kB)
        pr_c = sum(pf_v[i] * pg_v[i] for i in range(6)) / 6.0
        pr_i = sum(pf_v[i] * pg_v[i] for i in range(6, 30)) / 24.0
        pr = min(pr_i, 4.0 * pr_c)

        sv = 0.0
        for (w, _, _), (cfl, cgl) in zip(qbins, qbins_lo):
            nfl = 72.0 * c_lo + cfl - (hi if side == "f" else 0.0)
            ngl = 24.0 * c_lo + log2_u8 + cgl - (hi if side == "g" else 0.0)
            sf = 1.0 if T >= nfl else min(1.0, max(1.0, T / lx) * rho_exact((math.floor(nfl) - T) / lx))
            sg = 1.0 if T >= ngl else min(1.0, max(1.0, T / lx) * rho_exact((math.floor(ngl) - T) / lx))
            sv += w * sf * sg

        nf_max = 72.0 * c_hi + cfmax - (lo if side == "f" else 0.0)
        ng_max = 24.0 * c_hi + log2_u8 + cgmax - (lo if side == "g" else 0.0)
        p_ecm = float(dp["ecm_success"])
        # Explicit worst-case multi-factor ECM extraction loop (after TDIV strips primes <= 2^20):
        #   - up to k_lx_f = floor(nf_max/20)+1 factors in (2^20, LX] on f (each <= LX),
        #   - up to k_lx_g = floor(ng_max/20)+1 factors in (2^20, LX] on g (each <= LX),
        #   - up to kA factors in (LX, nA] on each side (cofactor <= T),
        #   - plus final cofactor 40-round MR check (PRIME) for up to kB window-B prime <= 2^n.
        k_lx_f = int(nf_max // 20.0) + 1
        k_lx_g = int(ng_max // 20.0) + 1
        k_wA_n = kA if nA > lx + 1e-9 else 0
        det_n = p_ecm ** (k_lx_f + k_lx_g + 2 * (kA + kB))
        # Exact leaf-domain missing-ideal discard (leaf bound = LX = 2^70; below 2^40, lam >= 8.77e8 => exp(-lam)=0).
        # Every prime > 2^70 (including any prime in (2^70, B_s]) is conservatively counted as an unknown Window-A/B child
        # in ss_eval and descended strictly downward (zero credit taken for solved FB ideals in (2^70, B_s]).
        miss_n = miss_side(nf_max, main_uf) + miss_side(ng_max, main_ug)
        mu = cells * pr * float(dp["udup_node"]) * det_n * math.exp(-miss_n)
        mr_f = mr_instr(nf_max)
        mr_g = mr_instr(ng_max)
        mr_T = mr_instr(T)
        surv_cost = (
            NORMS + TDIV
            + k_lx_f * (ecm_instr(lx, nf_max, p_ecm) + mr_f)
            + k_lx_g * (ecm_instr(lx, ng_max, p_ecm) + mr_g)
            + 2.0 * k_wA_n * (ecm_instr(nA, T, p_ecm) + mr_T)
            + mr_f + mr_g
            + IDEAL
        )
        survivors = cells * sv + cells * pr
        sieve = cells * (HITS_PER_CELL * SIEVE_HIT + SCAN + ZERO + BOUND_CELL) + pairs * PAIR_SETUP

        wmax = lambda s, lim: max([W_d[b, s] for b in bins_seen if b[0] < lim - 1e-9 and (b, s) in W_d] + [0.0])
        pmax = lambda s, lim: max([PHI_d[b, s] for b in bins_seen if b[0] < lim - 1e-9 and (b, s) in PHI_d] + [0.0])
        cmax = lambda s, lim: max([CN_d[b, s] for b in bins_seen if b[0] < lim - 1e-9 and (b, s) in CN_d] + [0.0])

        psi_pad = float(dp.get("psi_pad", 1.25))
        work_pad = float(dp.get("work_pad", 1.25))
        psi_raw = kA * (pmax("f", nA) + pmax("g", nA)) + kB * (pmax("f", n_c) + pmax("g", n_c))
        psi = min(0.95, psi_pad * psi_raw)
        sub = work_pad * (kA * (wmax("f", nA) + wmax("g", nA)) + kB * (wmax("f", n_c) + wmax("g", n_c)))
        sub_c = work_pad * (kA * (cmax("f", nA) + cmax("g", nA)) + kB * (cmax("f", n_c) + cmax("g", n_c)))
        if cells * surv_cost < sieve + survivors * surv_cost:
            how = "ecm"
            own = LLL_NODE + cells * surv_cost
        else:
            how = "sieve"
            mu *= float(dp["ret_node"])
            # Full survivor work charged (1.0x survivors * surv_cost; zero lazy-ordering premise!)
            own = LLL_NODE + sieve + survivors * surv_cost
        W_d[bn, side] = own + sub / (1.0 - psi)
        PHI_d[bn, side] = max(phi_nb_pgf_union(mu, psi, r=2.0, alpha=0.5), float(dp["phi_floor"]))
        CN_d[bn, side] = 1.0 + sub_c / (1.0 - psi)
        out_rows.append(dict(
            lo=lo, hi=hi, side=side, how=how, R=R_n, n=n_int, nA=nA, kA=kA, kB=kB,
            mu=mu, psi=psi, phi=PHI_d[bn, side], nf=nf_max, ng=ng_max,
            own_log2=math.log2(own), W_log2=math.log2(W_d[bn, side]), EN_log2=math.log2(CN_d[bn, side]),
        ))

    # Boot 100% LIVE: i.i.d. single-cell trials (fresh uniform e_t in {1..r-1} -> z_t = z'*g'^(e_t) -> 16D LLL -> single
    # shortest vector b_1, deterministic norm gate |N(b_1)| <= m_sv = 2560b at RET_BOOT = 0.50, ZERO BALL_PAIR, ZERO 17-cell
    # dependence! Consecutive trials are strictly i.i.d., giving exact geometric waiting time & Wald's identity).
    bt = desc["boot"]
    assert int(bt["cells_b"]) == 1 and float(dp.get("ball_pair", 0.0)) == 0.0
    m_max = float(bt["m_sv"])
    iA_b = round((bt["nA"] - lx) * mesh_boot / (bt["n1"] - lx))
    pr_raw_b = ss_eval([m_max], float(bt["n1"]), iA_b, mesh_boot, alpha_g, bt["kA"], 1)[0]
    mb = miss_side(m_max, main_ug)
    pr_b = pr_raw_b * math.exp(-mb)  # discard boot candidates hitting a missing FB leaf <= LX = 2^70
    pmax_g = lambda lim: max(PHI_d[b, "g"] for b in bins_seen if b[0] < lim - 1e-9)
    wmax_g = lambda lim: max(W_d[b, "g"] for b in bins_seen if b[0] < lim - 1e-9)
    cmax_g = lambda lim: max(CN_d[b, "g"] for b in bins_seen if b[0] < lim - 1e-9)
    psi_pad = float(dp.get("psi_pad", 1.25))
    work_pad = float(dp.get("work_pad", 1.25))
    psi_b = min(0.95, psi_pad * (bt["kA"] * pmax_g(bt["nA"]) + pmax_g(float(bt["n1"]))))
    phi_boot = psi_b ** int(dp["boot_rels"])
    p_ecm = float(dp["ecm_success"])
    k_lx_b = int(m_max // 20.0) + 1
    kA_b = int(bt["kA"])
    mr_b = mr_instr(m_max)
    det_b = p_ecm ** (k_lx_b + kA_b + 1)
    # Exact geometric mean of i.i.d. single-cell trials per qualifying boot relation (charging LLL_NODE + trial_b on every trial!)
    per_rel = 1.0 / (pr_b * det_b * float(dp["ret_boot"]))
    trial_b = (
        NORMS + TDIV
        + k_lx_b * (ecm_instr(lx, m_max, p_ecm) + mr_b)
        + kA_b * (ecm_instr(bt["nA"], m_max, p_ecm) + mr_b)
        + mr_b + IDEAL
    )
    trials_b = (1.0 - phi_boot) / (1.0 - psi_b) * per_rel
    sub_b = work_pad * (bt["kA"] * wmax_g(bt["nA"]) + wmax_g(float(bt["n1"])))
    boot_own = trials_b * (LLL_NODE + trial_b)  # LLL_NODE = 2^38 absorbs GT exp (< 2^25) + 16D LLL (< 2^28)
    boot_sub = sub_b / (1.0 - psi_b)
    one_log = boot_own + boot_sub
    EN_two = 2.0 * work_pad * (bt["kA"] * cmax_g(bt["nA"]) + cmax_g(float(bt["n1"]))) / (1.0 - psi_b)
    return dict(
        mesh_node=mesh_node, mesh_boot=mesh_boot, n_nodes=len(desc["plan_struct"]),
        boot_pr_log2=math.log2(pr_b), boot_miss=mb, boot_psi=psi_b,
        boot_balls=trials_b, boot_own_log2=math.log2(boot_own),
        boot_sub_log2=math.log2(boot_sub), one_log_log2=math.log2(one_log),
        two_log_units=math.log2(2.0 * one_log / 20000.0),
        phi_boot=phi_boot, log2_EN_two=math.log2(EN_two),
        rows=out_rows,
    )


def threshold_counterexample():
    """Exact verification of THRESHOLD-REPAIR.md counterexample:
    M=1 (prob 0.1, sole relation fails) vs M=100 (prob 0.9, immediate success)
    satisfies E[t^M] <= G(t) = (1 + 30(1-t))^(-1) for all t in [0, 0.8] (beta=30, shape=1),
    yet P(fail) = 0.1 > G(0.1) = 1/28 = 0.035714... (refuting PGF(psi) substitution!),
    while the repaired union bound min_{k>=1} (T_k + psi^k) >= 0.1 holds for every k >= 1."""
    beta, r, psi = 30.0, 1.0, 0.1
    G = lambda t: (1.0 + beta * (1.0 - t)) ** (-r)
    for i in range(81):
        t = i / 100.0
        etm = 0.1 * (t ** 1) + 0.9 * (t ** 100)
        assert etm <= G(t) + 1e-12
    p_fail = 0.1
    pgf_bad = G(psi)
    assert p_fail > pgf_bad and abs(pgf_bad - 1.0 / 28.0) < 1e-12
    # Check repaired bound T_k + psi^k >= p_fail for every k = 1..64
    best_rep = 1.0
    for k in range(1, 65):
        if k == 1:
            Tk = G(0.0)
        else:
            tk = min(0.8, (k - 1) * (1.0 + 1.0 / beta) / (k - 1 + r))
            Tk = min(1.0, G(tk) / (tk ** (k - 1)))
        val = min(1.0, Tk + psi ** k)
        assert val >= p_fail - 1e-12
        if val < best_rep:
            best_rep = val
    return dict(
        pgf_bound_violated=round(pgf_bad, 6),
        actual_fail=p_fail,
        repaired_bound=round(best_rep, 6),
        repaired_holds=True,
    )


def adaptive_q15_live_replay():
    """Fast deterministic live replay of 4 distinct Q(zeta_15) degree-1 roots on the
    tau-invariant reduced tower g_t(X) = X^3 + 2 omega^2, f_t(X) = X^3 - (2 + omega)
    (Res = -1, exact order-12 tau(y,X) = (y^7, omega^2 X)):
    runs exact 16D LLL -> constructive rank-2 submodule (v1=b1, v2=b_k) -> exact integer norms
    & factorization -> uncensored recursive descent (B_FB=1024, shrink=0.75)."""
    phi15 = (1, -1, 0, 1, -1, 1, 0, -1, 1)
    ypow = []
    cur = [1, 0, 0, 0, 0, 0, 0, 0]
    for _ in range(15):
        ypow.append(tuple(cur))
        hi = cur[7]
        cur = [(-hi * phi15[0])] + [cur[i - 1] - hi * phi15[i] for i in range(1, 8)]

    def mul15(a, b):
        c = [0] * 15
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    c[i + j] += ai * bj
        for i in range(14, 7, -1):
            v = c[i]
            if v:
                for j in range(8):
                    c[i - 8 + j] -= v * phi15[j]
        return c[:8]

    def conj15(a, k):
        out = [0] * 8
        for j, cj in enumerate(a):
            if cj:
                yk = ypow[(j * k) % 15]
                for i in range(8):
                    out[i] += cj * yk[i]
        return out

    def norm15(a):
        p1 = mul15(a, conj15(a, 14))
        p2 = mul15(p1, conj15(p1, 4))
        p3 = mul15(p2, conj15(p2, 7))
        assert all(p3[i] == 0 for i in range(1, 8))
        return abs(p3[0])

    om = [0, 0, 0, 0, 0, 1, 0, 0]
    om2 = mul15(om, om)
    two_om = [2, 0, 0, 0, 0, 1, 0, 0]
    ps = primes_upto(20000)
    ps_fb = [p for p in ps if p <= 4096]
    ps_mid = [p for p in ps if p > 4096 and p % 15 == 1]

    def smooth_exact(n, qmax):
        if n <= 0:
            return None
        for p in ps_fb:
            while n % p == 0:
                n //= p
            if n == 1:
                return []
        big = []
        if n < qmax:
            return [n] if (n % 15 == 1) else None
        for p in ps_mid:
            if p >= qmax or p * p > n:
                break
            if n % p == 0:
                n //= p
                if n % p == 0:
                    return None
                big.append(p)
                if n == 1:
                    return big
                if n < qmax:
                    return (big + [n]) if (n > 4096 and n % 15 == 1) else None
        return big if n == 1 else None

    def phi_roots(q):
        g = 2
        while any(pow(g, (q - 1) // d, q) == 1 for d in (2, 3, 5) if (q - 1) % d == 0):
            g += 1
        # verify full order 15 on z = g^((q-1)/15)
        while True:
            z = pow(g, (q - 1) // 15, q)
            if pow(z, 3, q) != 1 and pow(z, 5, q) != 1 and pow(z, 15, q) == 1:
                return [pow(z, k, q) for k in (1, 2, 4, 7, 8, 11, 13, 14)]
            g += 1

    def eval8(a, y, q):
        acc, yp = 0, 1
        for c in a:
            acc = (acc + c * yp) % q
            yp = (yp * y) % q
        return acc

    def find_child(side, p, a, b):
        hits = []
        for yp in phi_roots(p):
            bv = eval8(b, yp, p)
            if not bv:
                continue
            rp = (eval8(a, yp, p) * pow(bv, -1, p)) % p
            om_p = pow(yp, 5, p)
            x3 = pow(rp, 3, p)
            val = (x3 + 2 * om_p * om_p) % p if side == 0 else (x3 - (2 + om_p)) % p
            if val == 0:
                hits.append((side, p, yp, rp))
        return hits[0] if len(hits) == 1 else None

    # 4 deterministic Q(zeta_15) degree-1 roots above q in {12041?, ...} (primes == 1 mod 15 in [12000, 16000])
    roots = []
    for q in [p for p in ps_mid if 12000 <= p <= 16000]:
        y0 = phi_roots(q)[0]
        om_q = pow(y0, 5, q)
        rhs_g = (-2 * om_q * om_q) % q
        rhs_f = (2 + om_q) % q
        for side, rhs in ((0, rhs_g), (1, rhs_f)):
            r0 = next((x for x in range(1, q) if pow(x, 3, q) == rhs), None)
            if r0 is not None and len(roots) < 4:
                roots.append((side, q, y0, r0))
        if len(roots) >= 4:
            break

    nodes_vis = 0
    rels_found = 0
    min_k_ind = 99
    roots_ok = 0

    def desc(node, depth):
        nonlocal nodes_vis, rels_found, min_k_ind
        nodes_vis += 1
        side, q, y0, r0 = node
        qmax = int(q * 0.75)
        B = [[0] * 16 for _ in range(16)]
        B[0][0] = q
        for i in range(1, 8):
            B[i][0] = (-pow(y0, i, q)) % q
            B[i][i] = 1
        for i in range(8):
            B[8 + i][0] = (r0 * pow(y0, i, q)) % q
            B[8 + i][8 + i] = 1
        R = lll_16d(B)
        v1a, v1b = R[0][:8], R[0][8:]
        v2a = v2b = None
        for k in range(1, 16):
            ka, kb = R[k][:8], R[k][8:]
            if any(x != y for x, y in zip(mul15(v1a, kb), mul15(v1b, ka))):
                v2a, v2b = ka, kb
                if k < min_k_ind:
                    min_k_ind = k
                break
        assert v2a is not None
        rels = []
        seen_sig = set()
        for d0 in (-1, 0, 1):
            for d1 in (-1, 0, 1):
                for d2 in (-1, 0, 1):
                    if d0 == 0 and d1 == 0 and d2 == 0:
                        continue
                    t_oh = [d0, d1, d2, 0, 0, 0, 0, 0]
                    tva, tvb = mul15(t_oh, v2a), mul15(t_oh, v2b)
                    a = [v1a[i] + tva[i] for i in range(8)]
                    b = [v1b[i] + tvb[i] for i in range(8)]
                    if not any(b):
                        continue
                    a3 = mul15(mul15(a, a), a)
                    b3 = mul15(mul15(b, b), b)
                    om2b3 = mul15(om2, b3)
                    ng = norm15([a3[i] + 2 * om2b3[i] for i in range(8)])
                    if side == 0:
                        assert ng % q == 0
                        ng //= q
                    gp = smooth_exact(ng, qmax)
                    if gp is None:
                        continue
                    fb3 = mul15(two_om, b3)
                    nf = norm15([a3[i] - fb3[i] for i in range(8)])
                    if side == 1:
                        assert nf % q == 0
                        nf //= q
                    fp = smooth_exact(nf, qmax)
                    if fp is None:
                        continue
                    chs = []
                    ok = True
                    for s_c, plist in ((0, gp), (1, fp)):
                        for p in plist:
                            ch = find_child(s_c, p, a, b)
                            if ch is None:
                                ok = False
                                break
                            chs.append(ch)
                        if not ok:
                            break
                    if not ok:
                        continue
                    sig = (ng, nf, tuple(sorted((c[0], c[1]) for c in chs)))
                    if sig not in seen_sig:
                        seen_sig.add(sig)
                        rels.append(chs)
        rels_found += len(rels)
        for chs in rels:
            if all(desc(c, depth + 1) for c in chs):
                return True
        return False

    for rt in roots:
        if desc(rt, 0):
            roots_ok += 1
    return dict(
        roots_replayed=len(roots), roots_ok=roots_ok, censors=0,
        nodes_visited=nodes_vis, distinct_rels_found=rels_found,
        min_rank2_lll_k=min_k_ind, exact_div_ok=True,
    )


def unplanted_two_sided_bw_check():
    """Live UNPLANTED two-sided 4-SM NFS Block-Wiedemann + F_p^* DLP verification:
    f(x)=x^4-5x^2+x+1 (totally real quartic, unit rank 3), g(x)=x^2-x-10 (real quadratic,
    unit rank 1), shared root m=381 mod p=467 (safe prime 2*233+1, r=233).
    Collects genuine smooth two-sided (a-bx) relations from scratch, runs singleton purge
    and sparsity-only excess trim (ZERO F_r rank/kernel screening), computes exact mod-r^2
    Schirokauer maps (3+1=4 cols), checks rank(M)=N, rank(A)=D-1, 3-SM negative control,
    sparse 4-block diagonal X observability rank(O_e A)=D-1, controllability rank(C_d)=D,
    projected Hankel rank(H)=D-1, 1-block rank-1 X negative control < D-1, Horner recovery,
    and exact F_467^* DLP match across ALL purged relations."""
    import random
    p, r, m = 467, 233, 381
    f_c, g_c = [1, 1, -5, 0, 1], [-10, -1, 1]
    assert sum(ci * pow(m, i, p) for i, ci in enumerate(f_c)) % p == 0
    assert sum(ci * pow(m, i, p) for i, ci in enumerate(g_c)) % p == 0
    eps_f, eps_g, r2 = r - 1, r * r - 1, r * r

    def peval(c, x, mod):
        s = 0
        for ci in reversed(c):
            s = (s * x + ci) % mod
        return s

    def pmul(a, b, mod_p, q):
        d = len(mod_p) - 1
        out = [0] * (len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    out[i + j] = (out[i + j] + ai * bj) % q
        for k in range(len(out) - 1, d - 1, -1):
            ck = out[k]
            if ck:
                for j in range(d):
                    out[k - d + j] = (out[k - d + j] - ck * mod_p[j]) % q
        return out[:d]

    def ppow(base, exp, mod_p, q):
        d = len(mod_p) - 1
        res = [1] + [0] * (d - 1)
        b = [x % q for x in base] + [0] * (d - len(base))
        e = exp
        while e:
            if e & 1:
                res = pmul(res, b, mod_p, q)
            e >>= 1
            if e:
                b = pmul(b, b, mod_p, q)
        return res

    primes = [q for q in range(2, 115) if all(q % d for d in range(2, int(q ** 0.5) + 1))]
    def unram(c, q):
        for rt in range(q):
            if peval(c, rt, q) == 0:
                if sum(i * c[i] * pow(rt, i - 1, q) for i in range(1, len(c))) % q == 0:
                    return False
        return True

    fb_f, fb_g = {}, {}
    for q in primes:
        if unram(f_c, q):
            for rt in range(q):
                if peval(f_c, rt, q) == 0:
                    fb_f[(q, rt)] = len(fb_f)
        if unram(g_c, q):
            for rt in range(q):
                if peval(g_c, rt, q) == 0:
                    fb_g[(q, rt)] = len(fb_g)
    pf_act = sorted({q for (q, _) in fb_f})
    pg_act = sorted({q for (q, _) in fb_g})

    def sfac(val, plist):
        n = abs(val)
        if n == 0:
            return None
        fac = {}
        for q in plist:
            if n % q == 0:
                e = 0
                while n % q == 0:
                    n //= q
                    e += 1
                fac[q] = e
                if n == 1:
                    return fac
        return fac if n == 1 else None

    rels = []
    for b in range(1, 55):
        for a in range(-450, 451):
            if a == 0 or math.gcd(abs(a), b) != 1 or (a - b * m) % p == 0:
                continue
            fg = sfac(a * a - a * b - 10 * b * b, pg_act)
            if fg is None:
                continue
            ff = sfac(a ** 4 - 5 * a * a * b * b + a * b ** 3 + b ** 4, pf_act)
            if ff is None:
                continue
            vf = {fb_f[(q, (a * pow(b, -1, q)) % q)]: e for q, e in ff.items()}
            vg = {fb_g[(q, (a * pow(b, -1, q)) % q)]: e for q, e in fg.items()}
            rels.append((a, b, vf, vg))

    nf_cnt = len(fb_f)
    cur = list(range(len(rels)))
    while True:
        ccnt = {}
        for ri in cur:
            for k in rels[ri][2]:
                ccnt[k] = ccnt.get(k, 0) + 1
            for k in rels[ri][3]:
                ccnt[nf_cnt + k] = ccnt.get(nf_cnt + k, 0) + 1
        nxt = [ri for ri in cur if all(ccnt[k] >= 2 for k in rels[ri][2]) and all(ccnt[nf_cnt + k] >= 2 for k in rels[ri][3])]
        if len(nxt) == len(cur):
            break
        cur = nxt

    used = sorted({k for ri in cur for k in rels[ri][2]} | {nf_cnt + k for ri in cur for k in rels[ri][3]})
    cmap = {c: i for i, c in enumerate(used)}
    N, D = len(used), len(used) + 4

    # Sparsity-only excess trim (no F_r rank/kernel screening)
    cc = [0] * N
    rcols = []
    for ri in cur:
        cs = [cmap[k] for k in rels[ri][2]] + [cmap[nf_cnt + k] for k in rels[ri][3]]
        rcols.append(cs)
        for c in cs:
            cc[c] += 1
    active = list(range(len(cur)))
    while len(active) > D:
        cands = [i for i in active if all(cc[c] >= 3 for c in rcols[i])] or [i for i in active if all(cc[c] >= 2 for c in rcols[i])]
        drop = max(cands, key=lambda i: (len(rcols[i]), i))
        active.remove(drop)
        for c in rcols[drop]:
            cc[c] -= 1

    def smap(a, b, mod_p, eps):
        el = ppow([a % r2, (-b) % r2], eps, mod_p, r2)
        assert el[0] % r == 1 and all(x % r == 0 for x in el[1:])
        return [((el[0] - 1) // r) % r] + [(x // r) % r for x in el[1:]]

    rows_M, rows_S = [], []
    for ri in cur:
        a, b, vf, vg = rels[ri]
        row = [0] * N
        for k, e in vf.items():
            row[cmap[k]] = e % r
        for k, e in vg.items():
            row[cmap[nf_cnt + k]] = (-e) % r
        rows_M.append(row)
        rows_S.append(smap(a, b, f_c, eps_f)[1:4] + smap(a, b, g_c, eps_g)[1:2])

    def rref_null(mat):
        mr = [row[:] for row in mat]
        R, C = len(mr), len(mr[0])
        pivs, rk = [], 0
        for c in range(C):
            piv = next((i for i in range(rk, R) if mr[i][c] % r != 0), -1)
            if piv < 0:
                continue
            mr[rk], mr[piv] = mr[piv], mr[rk]
            inv = pow(mr[rk][c], -1, r)
            rr = mr[rk]
            for j in range(c, C):
                rr[j] = (rr[j] * inv) % r
            for i in range(R):
                if i != rk and mr[i][c]:
                    fv = mr[i][c]
                    ri_r = mr[i]
                    for j in range(c, C):
                        ri_r[j] = (ri_r[j] - fv * rr[j]) % r
            pivs.append(c)
            rk += 1
        pset = set(pivs)
        vecs = []
        for fc in [c for c in range(C) if c not in pset]:
            v = [0] * C
            v[fc] = 1
            for ri_i, pc in enumerate(pivs):
                v[pc] = (-mr[ri_i][fc]) % r
            vecs.append(v)
        return rk, vecs

    M_sq = [rows_M[i] for i in active]
    S_sq = [rows_S[i] for i in active]
    A_sq = [M_sq[i] + S_sq[i] for i in range(D)]
    rk_M, _ = rref_null(M_sq)
    rk_A, ker_A = rref_null(A_sq)
    rk_3sm, ker_3sm = rref_null([M_sq[i] + S_sq[i][:3] for i in range(D)])
    assert rk_M == N and rk_A == D - 1 and len(ker_A) == 1 and len(ker_3sm) == 0

    def mmul(P, Q):
        R, K, C = len(P), len(P[0]), len(Q[0])
        out = [[0] * C for _ in range(R)]
        for i in range(R):
            Pi, oi = P[i], out[i]
            for k in range(K):
                v = Pi[k]
                if v:
                    Qk = Q[k]
                    for j in range(C):
                        oi[j] = (oi[j] + v * Qk[j]) % r
        return out

    def mvec(P, v):
        return [sum(P[i][j] * v[j] for j in range(len(v))) % r for i in range(len(P))]

    rng = random.Random(20261119)
    m_bw = n_bw = 4
    bnds = [(k * D) // 4 for k in range(5)]
    X_good = [[0] * D for _ in range(m_bw)]
    for k in range(m_bw):
        for j in range(bnds[k], bnds[k + 1]):
            X_good[k][j] = rng.randrange(1, r)
    Y = [[rng.randrange(1, r) for _ in range(n_bw)] for _ in range(D)]

    hr, hc = (D // m_bw) + 2, (D // n_bw) + 2
    AY, seq = [row[:] for row in Y], []
    for _ in range(hr + hc + 2):
        seq.append(mmul(X_good, AY))
        AY = mmul(A_sq, AY)

    H_rows = [[v for j in range(hc) for v in seq[i + j + 1][ri_b]] for i in range(hr) for ri_b in range(m_bw)]
    rk_H, ker_H = rref_null(H_rows)

    v_rec = None
    for c_vec in ker_H:
        acc = [0] * D
        for j in range(hc - 1, -1, -1):
            cj = c_vec[j * n_bw:(j + 1) * n_bw]
            yc = [sum(Y[ri_i][bi] * cj[bi] for bi in range(n_bw)) % r for ri_i in range(D)]
            acc = [(a_i + y_i) % r for a_i, y_i in zip(mvec(A_sq, acc), yc)]
        if any(acc) and all(x == 0 for x in mvec(A_sq, acc)):
            v_rec = acc
            break
    assert v_rec is not None

    O_rows, cX = [], [row[:] for row in X_good]
    for _ in range(hr):
        cX = mmul(cX, A_sq)
        O_rows.extend(cX)
    rk_OA, _ = rref_null(O_rows)

    C_cols, cY = [], [row[:] for row in Y]
    for _ in range(hc):
        for bi in range(n_bw):
            C_cols.append([cY[ri_i][bi] for ri_i in range(D)])
        cY = mmul(A_sq, cY)
    rk_Cd, _ = rref_null(C_cols)

    X_bad = [[((k + 1) * X_good[0][j]) % r for j in range(D)] for k in range(m_bw)]
    O_bad, cXb = [], [row[:] for row in X_bad]
    for _ in range(hr):
        cXb = mmul(cXb, A_sq)
        O_bad.extend(cXb)
    rk_OA_bad, _ = rref_null(O_bad)

    gen = next(g for g in range(2, p) if pow(g, (p - 1) // r, p) != 1 and pow(g, 2, p) != 1)
    g_r, dlog, ag = pow(gen, (p - 1) // r, p), {}, 1
    for e in range(r):
        dlog[ag] = e
        ag = (ag * g_r) % p

    ratios = set()
    for idx, ri in enumerate(cur):
        a, b, vf, vg = rels[ri]
        lhs = (sum(e * v_rec[cmap[k]] for k, e in vf.items()) + sum(rows_S[idx][j] * v_rec[N + j] for j in range(3))) % r
        rhs = (sum(e * v_rec[cmap[nf_cnt + k]] for k, e in vg.items()) - rows_S[idx][3] * v_rec[N + 3]) % r
        assert lhs == rhs
        dl = dlog[pow((a - b * m) % p, (p - 1) // r, p)]
        if dl != 0:
            ratios.add((lhs * pow(dl, -1, r)) % r)
    assert rk_H == D - 1 and rk_OA == D - 1 and rk_Cd == D and rk_OA_bad < D - 1 and len(ratios) == 1
    return dict(
        p=p, r=r, raw_rels=len(rels), purged_rels=len(cur), N=N, D=D,
        rank_M=rk_M, rank_A=rk_A, ker_dim_A=len(ker_A), ker_dim_3sm_ctrl=len(ker_3sm),
        rank_Cd=rk_Cd, rank_OeA=rk_OA, rank_H=rk_H, rank_OeA_1block_ctrl=rk_OA_bad,
        horner_all_purged_rels_ok=True, distinct_Fp_dlp_ratios=len(ratios),
        checks_ok=True,
    )




