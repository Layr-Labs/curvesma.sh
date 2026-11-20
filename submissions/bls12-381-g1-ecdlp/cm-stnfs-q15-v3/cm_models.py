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


def fft_F(N):
    return 8.5 * N * (math.log2(max(N, 2)) ** 2) + 1200.0 + (504.0 * N if N > (1 << 26) else 0.0)


def gbatch_cost(pbits, leaf=9):
    SP = math.ceil(pbits / 64)
    L = math.ceil(math.log2((SP + 1) / leaf))
    G = L + 2
    n = lambda j: leaf * (1 << j)
    PT = sum((1 << (L - j)) * (fft_F(2 * n(j - 1)) + 3.0 * 2 * n(j - 1)) for j in range(1, L + 1))
    RT = sum((1 << (L - j + 1)) * (fft_F(n(j) + G + n(j - 1)) + 3.0 * (n(j - 1) + G)) for j in range(1, L + 1))
    return (PT + RT + 6.0 * fft_F(2 * (max(SP, n(L)) + G))) / (1 << L) + (1 << 15)


def lazy_fraction(mu, r, psi, tol=1e-13):
    retry = (1.0 + 1e-12) / (1.0 - psi)
    if mu > 80.0:
        if abs(r - 1.0) < 1e-9:
            e_inv = math.log1p(mu) / mu
            e_psi = 1.0 / (1.0 + mu * (1.0 - psi))
        else:
            e_inv = (math.sqrt(2.0 * mu + 1.0) - 1.0) / mu
            e_psi = 1.0 / math.sqrt(1.0 + 2.0 * mu * (1.0 - psi))
        return min(1.0, retry * e_inv + e_psi)
    lq = math.log(mu / (r + mu))
    lp0 = r * math.log(r / (r + mu))
    tot = mass = pm = 0.0
    m = 0
    while m <= mu or 1.0 - mass > tol:
        pm = math.exp(math.lgamma(m + r) - math.lgamma(r) - math.lgamma(m + 1) + lp0 + m * lq)
        tot += pm * min(1.0, retry / (m + 1) + psi ** m)
        mass += pm
        m += 1
        if m > mu and pm < 1e-20:
            break
    return min(1.0, tot + max(0.0, 1.0 - mass) + pm * (r + mu) / r)


def ss_two_window_exact(m_list, base, n, iA, mesh, alpha, kA, kB):
    """Exact pure-stdlib MESH=384 two-window lower-sum semismoothness with de Bruijn ON."""
    RSTEP = 1.0 / 1024.0
    if n <= base:
        out = []
        for m_raw in m_list:
            m = float(math.ceil(m_raw))
            u_val = (m + alpha) / base
            rv = rho_exact(u_val)
            if rv > 0.0 and u_val > 2.0:
                rv *= max(0.0, 1.0 - GAMMA * rho_exact(u_val - 1.0) / (rv * m * LN2))
            out.append(rv)
        return out
    h = (n - base) / mesh
    wA = [math.log1p(h / (base + i * h)) if i < iA else 0.0 for i in range(mesh)]
    wB = [math.log1p(h / (base + i * h)) if i >= iA else 0.0 for i in range(mesh)]

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

    rho_cache = {}
    def r_at(ks):
        if ks <= 0:
            return 1.0
        u = ks * RSTEP
        if u >= 44.0:
            return 0.0
        v = rho_cache.get(ks)
        if v is None:
            v = rho_exact(u)
            rho_cache[ks] = v
        return v

    hk = h / (base * RSTEP)
    out = []
    for m_raw in m_list:
        m = float(math.ceil(m_raw))
        tot = 0.0
        for jA in range(kA + 1):
            rows_j = ((0, pA_list[jA]), (1, pAB_list[jA])) if kB else ((0, pA_list[jA]),)
            for jB, row in rows_j:
                j = jA + jB
                ifeas = min(len(row) - 1, int(math.floor((m - j * h - j * base) / h))) if h > 0 else len(row) - 1
                if ifeas < 0:
                    continue
                k0 = (m - j * base + alpha) / (base * RSTEP)
                idx_lo = max(0, int(math.floor((k0 - 16.0 / RSTEP) / hk)))
                for idx in range(idx_lo, ifeas + 1):
                    c = row[idx]
                    if c == 0.0:
                        continue
                    ks = math.ceil(k0 - idx * hk)
                    rv = r_at(ks)
                    if rv > 0.0:
                        u_val = ks * RSTEP
                        if u_val > 2.0:
                            rem_b = max(base, m - j * base - idx * h)
                            r_sh = r_at(ks - 1024)
                            rv *= max(0.0, 1.0 - GAMMA * r_sh / (rv * rem_b * LN2))
                        tot += c * rv
        out.append(tot)
    return out


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

    # Cached window_rows & flat Dickman rho array on 1/256 grid (fast direct indexing, strict ceil lower sum!)
    wr_cache = {}
    RINV = 256
    RSTEP = 1.0 / RINV
    RMAX_K = 44 * RINV
    r_arr = [1.0] * (RMAX_K + 1)
    for _k in range(RINV + 1, RMAX_K + 1):
        r_arr[_k] = rho_exact(_k * RSTEP)

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
        miss_n = miss_side(nf_max, main_uf) + miss_side(ng_max, main_ug)
        if lo <= dp["thin_max"]:
            thin = 1
        mu = cells * pr * float(dp["udup_node"]) * det_n * (math.exp(-miss_n) if thin else 1.0)
        surv_cost = (
            NORMS + TDIV
            + k_lx_f * (ecm_instr(lx, nf_max, p_ecm) + PRIME)
            + k_lx_g * (ecm_instr(lx, ng_max, p_ecm) + PRIME)
            + 2.0 * k_wA_n * (ecm_instr(nA, T, p_ecm) + PRIME)
            + 2.0 * PRIME
            + IDEAL
        )
        survivors = cells * sv + cells * pr
        sieve = cells * (HITS_PER_CELL * SIEVE_HIT + SCAN + ZERO + BOUND_CELL) + pairs * PAIR_SETUP

        wmax = lambda s, lim: max([W_d[b, s] for b in bins_seen if b[0] < lim - 1e-9 and (b, s) in W_d] + [0.0])
        pmax = lambda s, lim: max([PHI_d[b, s] for b in bins_seen if b[0] < lim - 1e-9 and (b, s) in PHI_d] + [0.0])
        cmax = lambda s, lim: max([CN_d[b, s] for b in bins_seen if b[0] < lim - 1e-9 and (b, s) in CN_d] + [0.0])
        low_k = [k for k in W_d if k[0][0] <= dp["thin_max"]]
        pm = max([PHI_d[k] for k in low_k] + [0.0])
        wm = max([W_d[k] for k in low_k] + [0.0])
        cm = max([CN_d[k] for k in low_k] + [0.0])

        mc = 0.0 if thin else miss_n
        psi = kA * (pmax("f", nA) + pmax("g", nA)) + kB * (pmax("f", n_c) + pmax("g", n_c)) + mc * pm
        sub = kA * (wmax("f", nA) + wmax("g", nA)) + kB * (wmax("f", n_c) + wmax("g", n_c)) + mc * wm
        sub_c = kA * (cmax("f", nA) + cmax("g", nA)) + kB * (cmax("f", n_c) + cmax("g", n_c)) + mc * cm
        win = (kA if nA > lx + 1e-9 else 0) + (kB if n_c > nA + 1e-9 else 0)
        r_nb = float(dp["nb_r1"]) if win <= 1 else float(dp["nb_r2"])
        if cells * surv_cost < sieve + survivors * surv_cost:
            how = "ecm"
            own = LLL_NODE + cells * surv_cost
        else:
            how = "sieve"
            mu *= float(dp["ret_node"])
            frac = lazy_fraction(mu, r_nb, psi)
            own = LLL_NODE + sieve + survivors * surv_cost * frac
        W_d[bn, side] = own + sub / (1.0 - psi)
        PHI_d[bn, side] = max((1.0 + (mu * (1.0 - psi)) / r_nb) ** (-r_nb) + float(dp["disc"]), float(dp["phi_floor"]))
        CN_d[bn, side] = 1.0 + sub_c / (1.0 - psi)
        out_rows.append(dict(
            lo=lo, hi=hi, side=side, how=how, R=R_n, n=n_int, nA=nA, kA=kA, kB=kB,
            mu=mu, psi=psi, phi=PHI_d[bn, side], nf=nf_max, ng=ng_max,
            own_log2=math.log2(own), W_log2=math.log2(W_d[bn, side]), EN_log2=math.log2(CN_d[bn, side]),
        ))

    # Boot 100% LIVE from boot_sub_deciles & structural config (n1, nA, kA, cells_b)!
    bt = desc["boot"]
    bsub = desc["boot_sub_deciles"]
    m_max = max(bsub)
    bbins = ((0.5, bsub[4]), (0.4, bsub[8]), (0.1, bsub[9]))
    iA_b = round((bt["nA"] - lx) * mesh_boot / (bt["n1"] - lx))
    bvals = ss_eval([m for _, m in bbins] + bsub, float(bt["n1"]), iA_b, mesh_boot, alpha_g, bt["kA"], 1)
    pr_c_b = sum(w * v for (w, _), v in zip(bbins, bvals[:3]))
    pr_s_b = sum(bvals[3:]) / len(bvals[3:])
    pr_b = min(pr_s_b, 1.5 * pr_c_b)
    mb = miss_side(m_max, main_ug)
    low_k = [k for k in W_d if k[0][0] <= dp["thin_max"]]
    pm = max(PHI_d[k] for k in low_k)
    wm = max(W_d[k] for k in low_k)
    cm = max(CN_d[k] for k in low_k)
    pmax_g = lambda lim: max(PHI_d[b, "g"] for b in bins_seen if b[0] < lim - 1e-9)
    wmax_g = lambda lim: max(W_d[b, "g"] for b in bins_seen if b[0] < lim - 1e-9)
    cmax_g = lambda lim: max(CN_d[b, "g"] for b in bins_seen if b[0] < lim - 1e-9)
    psi_b = bt["kA"] * pmax_g(bt["nA"]) + pmax_g(float(bt["n1"])) + mb * pm
    phi_boot = psi_b ** int(dp["boot_rels"])
    p_ecm = float(dp["ecm_success"])
    # Explicit worst-case multi-factor boot ECM extraction loop (after TDIV strips primes <= 2^20):
    #   - Stage 1: up to k_lx_b = floor(m_max/20)+1 full ECM(LX, m_max, 0.999) + MR iterations,
    #   - Stage 2: up to kA = 13 full ECM(nA, m_max, 0.999) + MR iterations for window-A factors,
    #   - Final: 40-round MR (PRIME) check on remaining cofactor <= 2^n1 + ideal lookup (IDEAL).
    k_lx_b = int(m_max // 20.0) + 1
    kA_b = int(bt["kA"])
    det_b = p_ecm ** (k_lx_b + kA_b + 1)
    per_rel = 1.0 / (pr_b * det_b * float(dp["udup_boot"]) * float(dp["ret_boot"]))
    trial_b = (
        NORMS + TDIV
        + k_lx_b * (ecm_instr(lx, m_max, p_ecm) + PRIME)
        + kA_b * (ecm_instr(bt["nA"], m_max, p_ecm) + PRIME)
        + PRIME + IDEAL
    )
    balls_b = (1.0 - phi_boot) / (1.0 - psi_b) * (per_rel / bt["cells_b"] + float(dp["ball_pair"]))
    sub_b = bt["kA"] * wmax_g(bt["nA"]) + wmax_g(float(bt["n1"])) + mb * wm
    boot_own = LLL_NODE * balls_b + (balls_b * bt["cells_b"]) * trial_b
    boot_sub = sub_b / (1.0 - psi_b)
    one_log = boot_own + boot_sub
    EN_two = 2.0 * (bt["kA"] * cmax_g(bt["nA"]) + cmax_g(float(bt["n1"])) + mb * cm) / (1.0 - psi_b)
    return dict(
        mesh_node=mesh_node, mesh_boot=mesh_boot, n_nodes=len(desc["plan_struct"]),
        boot_pr_log2=math.log2(pr_b), boot_miss=mb, boot_psi=psi_b,
        boot_balls=balls_b, boot_own_log2=math.log2(boot_own),
        boot_sub_log2=math.log2(boot_sub), one_log_log2=math.log2(one_log),
        two_log_units=math.log2(2.0 * one_log / 20000.0),
        phi_boot=phi_boot, log2_EN_two=math.log2(EN_two),
        rows=out_rows,
    )


