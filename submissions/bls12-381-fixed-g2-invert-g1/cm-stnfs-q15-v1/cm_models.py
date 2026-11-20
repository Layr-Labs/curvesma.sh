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


def rho_exact(u):
    if u <= 1.0:
        return 1.0
    k = math.ceil(u)
    if k >= len(RHO_C):
        return 0.0
    x, v = k - u, 0.0
    for c in reversed(RHO_C[k]):
        v = v * x + c
    return max(v, 0.0)


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
