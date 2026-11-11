"""Stage-1-only ECM with per-prime gcd checkpoints. Pure Python, stdlib only (3.9-3.12).

Curve:  B y^2 = x^3 + A x^2 + x  (Montgomery), XZ coordinates, A24 = (A+2)/4.
Suyama: u = s^2-5, v = 4s, (X0:Z0) = (u^3 : v^3), A24 = (v-u)^3 (3u+v) / (16 u^3 v) mod N.
xDBL:   a=(X+Z)^2, b=(X-Z)^2, c=a-b, X'=a*b, Z'=c*(b+A24*c)       -> 2M + 2S + 1 mul by A24 = 5
xADD:   u=(X0-Z0)(X1+Z1), v=(X0+Z0)(X1-Z1), X'=Zd*(u+v)^2, Z'=Xd*(u-v)^2  -> 4M + 2S = 6
Stage 1: for each prime q <= B1 (increasing) multiply by q^k (q^k <= B1 < q^(k+1)), then
checkpoint d = gcd(Z, N).  q = 2 is done by k plain doublings, odd q^k by a Montgomery
ladder that costs (bits-1) xDBL + (bits-1) xADD.  No stage 2.  The modulus is N throughout.
"""
import json
import math
import random
from math import gcd

# ----------------------------------------------------------------------------- primes


def primes_upto(n):
    if n < 2:
        return []
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i :: i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(2, n + 1) if s[i]]


_PLAN = {}


def plan(B1):
    """(k2, [odd prime powers q^k <= B1 < q^(k+1) in increasing q])."""
    p = _PLAN.get(B1)
    if p is None:
        k2 = B1.bit_length() - 1 if B1 >= 2 else 0
        sc = []
        for q in primes_upto(B1)[1:]:
            m = q
            while m * q <= B1:
                m *= q
            sc.append(m)
        if len(_PLAN) > 64:
            _PLAN.clear()
        p = _PLAN[B1] = (k2, sc)
    return p


def stage1_bits(B1):
    """Exact sum over primes q <= B1 of bit_length(q^k), q^k <= B1 < q^(k+1)."""
    k2, sc = plan(B1)
    return (k2 + 1 if k2 else 0) + sum(m.bit_length() for m in sc)


def stage1_ops(B1):
    """Exact per-curve operation counts of stage 1 (independent of N and sigma)."""
    k2, sc = plan(B1)
    lad = sum(m.bit_length() - 1 for m in sc)
    dbl, add = k2 + lad, lad
    return {"xdbl": dbl, "xadd": add, "mulmod": 5 * dbl + 6 * add,
            "mul": 2 * dbl + 4 * add, "sqr": 2 * dbl + 2 * add, "mul_a24": dbl,
            "gcd_checkpoints": (1 if k2 else 0) + len(sc)}


_MR = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n):
    """Miller-Rabin, bases 2..37: deterministic below 3.3e24, probabilistic above."""
    if n < 2:
        return False
    for p in _MR:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def rand_prime(rng, bits):
    while True:
        n = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(n):
            return n


# ----------------------------------------------------------------------------- Dickman

_RHO = None
_RHO_N = 2000  # grid points per unit


def _rho_table(umax=12):
    n = _RHO_N
    h = 1.0 / n
    r = [1.0] * (n + 1)
    win = sum(r[1:n])  # sum of r[i-n+1 .. i-1] for i = n
    for i in range(n + 1, umax * n + 1):
        win += r[i - 1] - r[i - n]
        # u rho(u) = int_{u-1}^{u} rho  (trapezoid, implicit in r[i])
        r.append(h * (0.5 * r[i - n] + win) / (i * h - 0.5 * h))
    return r


def rho(u):
    """Dickman rho by the integral equation u rho(u) = int_{u-1}^u rho(t) dt."""
    global _RHO
    if u <= 1:
        return 1.0
    if _RHO is None:
        _RHO = _rho_table()
    x = u * _RHO_N
    i = int(x)
    if i >= len(_RHO) - 1:
        return 0.0
    f = x - i
    return _RHO[i] * (1 - f) + _RHO[i + 1] * f


# ----------------------------------------------------------------------------- curve


def suyama(sigma, N):
    """Return ((X0, Z0, A24), 1) or (None, g) when 16 u^3 v is not invertible (g = gcd > 1)."""
    u = (sigma * sigma - 5) % N
    v = 4 * sigma % N
    u3 = u * u * u % N
    den = 16 * u3 * v % N
    g = gcd(den, N)
    if g != 1:
        return None, g
    t = (v - u) % N
    a24 = t * t * t * (3 * u + v) * pow(den, -1, N) % N
    return (u3, v * v * v % N, a24), 1


def xdbl(X, Z, A24, n):
    t = X + Z
    a = t * t % n
    t = X - Z
    b = t * t % n
    c = a - b
    return a * b % n, c * (b + A24 * c) % n


def ladder(m, X, Z, A24, n):
    """[m](X:Z) for odd m >= 3.  (bits-1) xDBL + (bits-1) xADD."""
    X0 = X
    Z0 = Z
    t = X + Z
    a = t * t % n
    t = X - Z
    b = t * t % n
    c = a - b
    X1 = a * b % n
    Z1 = c * (b + A24 * c) % n
    for bit in bin(m)[3:-1]:
        u = (X0 - Z0) * (X1 + Z1) % n
        v = (X0 + Z0) * (X1 - Z1) % n
        t = u + v
        Xa = Z * (t * t % n) % n
        t = u - v
        Za = X * (t * t % n) % n
        if bit == "1":
            t = X1 + Z1
            a = t * t % n
            t = X1 - Z1
            b = t * t % n
            c = a - b
            X1 = a * b % n
            Z1 = c * (b + A24 * c) % n
            X0 = Xa
            Z0 = Za
        else:
            t = X0 + Z0
            a = t * t % n
            t = X0 - Z0
            b = t * t % n
            c = a - b
            X0 = a * b % n
            Z0 = c * (b + A24 * c) % n
            X1 = Xa
            Z1 = Za
    u = (X0 - Z0) * (X1 + Z1) % n
    v = (X0 + Z0) * (X1 - Z1) % n
    t = u + v
    Xa = Z * (t * t % n) % n
    t = u - v
    return Xa, X * (t * t % n) % n


def curve_factors(N, sigma, B1):
    """One curve on the full modulus N.  Returns (recorded factors, inversion_ok)."""
    cv, g = suyama(sigma, N)
    if cv is None:
        return ([g] if g < N else []), False
    X, Z, A24 = cv
    k2, sc = plan(B1)
    out = []
    dprev = 1
    for _ in range(k2):
        X, Z = xdbl(X, Z, A24, N)
    if k2:
        d = gcd(Z, N)
        if d != dprev:
            out.append(d // dprev)
            dprev = d
    for m in sc:
        X, Z = ladder(m, X, Z, A24, N)
        d = gcd(Z, N)
        if d != dprev:
            out.append(d // dprev)
            dprev = d
    return out, True


def capture(l, sigma, B1):
    """Ladder modulo the prime l itself.  1 = captured (Z = 0 mod l at some checkpoint),
    0 = not captured, -1 = degenerate sigma (16 u^3 v = 0 mod l).  Stops at the capture:
    (X:0) stays (X':0) under xDBL/xADD, so the outcome equals the full run."""
    cv, g = suyama(sigma, l)
    if cv is None:
        return -1
    X, Z, A24 = cv
    k2, sc = plan(B1)
    for _ in range(k2):
        X, Z = xdbl(X, Z, A24, l)
    if not Z:
        return 1
    for m in sc:
        X, Z = ladder(m, X, Z, A24, l)
        if not Z:
            return 1
    return 0


# ----------------------------------------------------------------------------- split


def coprime_base(xs):
    """Pairwise-coprime refinement: replace a, b with gcd g > 1 by a/g, g, b/g until stable."""
    base = []
    work = [x for x in xs if x > 1]
    ngcd = 0
    while work:
        x = work.pop()
        if x == 1:
            continue
        for i, b in enumerate(base):
            g = gcd(x, b)
            ngcd += 1
            if g != 1:
                base.pop(i)
                work.extend(t for t in (x // g, g, b // g) if t > 1)
                break
        else:
            base.append(x)
    return sorted(base), ngcd


def split(N, B1, curves, seed):
    """Run `curves` stage-1 curves on N.  Returns dict with the pairwise-coprime base built
    from all recorded checkpoint factors, the remaining cofactor and exact operation counts."""
    rng = random.Random(seed)
    rec = []
    ok = 0
    for _ in range(curves):
        f, good = curve_factors(N, rng.randrange(6, 1 << 64), B1)
        rec.extend(f)
        ok += good
    base, ngcd = coprime_base(rec)
    cof = N
    for b in base:
        while cof % b == 0:
            cof //= b
    o = stage1_ops(B1)
    ops = {"curves": curves, "curves_run": ok, "xdbl": ok * o["xdbl"], "xadd": ok * o["xadd"],
           "mulmod": ok * o["mulmod"], "mul": ok * o["mul"], "sqr": ok * o["sqr"],
           "mul_a24": ok * o["mul_a24"], "setup_mulmod": 11 * ok,
           "gcd": curves + ok * o["gcd_checkpoints"], "gcd_refine": ngcd, "inversions": ok}
    return {"base": base, "cofactor": cof, "recorded": len(rec), "ops": ops}


# ----------------------------------------------------------------------------- measurements


def b1_for(l, u):
    return max(2, int(math.exp(math.log(l) / u)))


def cell_prime(bits, i):
    return rand_prime(random.Random("ecm-prime-%d-%d" % (bits, i)), bits)


def measure_prime(l, u, curves, seed):
    """Capture counts for one prime l at B1 = floor(l^(1/u))."""
    B1 = b1_for(l, u)
    rng = random.Random(seed)
    cap = deg = 0
    for _ in range(curves):
        r = capture(l, rng.randrange(6, 1 << 64), B1)
        if r > 0:
            cap += 1
        elif r < 0:
            deg += 1
    return {"l": l, "B1": B1, "curves": curves, "captures": cap, "degenerate": deg}


def curves_for(p, target=0.98):
    return int(math.ceil(math.log(1 - target) / math.log(1 - p)))


def planted(seed, total_bits=1700, lo=22, hi=48, u=3.5):
    """N = product of distinct primes of lo..hi bits, about total_bits bits in all."""
    rng = random.Random("ecm-planted-%s" % seed)
    k = rng.randint(3, 8)
    sizes = [rng.randint(lo, hi) for _ in range(k)]
    rem = total_bits - sum(sizes)
    while rem > 2 * hi:
        s = rng.randint(lo, hi)
        sizes.append(s)
        rem -= s
    if rem > hi:
        sizes += [rem // 2, rem - rem // 2]
    elif rem >= lo:
        sizes.append(rem)
    ps = set()
    for s in sizes:
        while True:
            p = rand_prime(rng, s)
            if p not in ps:
                ps.add(p)
                break
    ps = sorted(ps)
    N = 1
    for p in ps:
        N *= p
    B1 = b1_for(ps[-1], u)
    ua = math.log(ps[-1]) / math.log(B1)
    return {"N": N, "primes": ps, "k_main": k, "B1": B1, "curves": curves_for(rho(ua))}


def planted_run(seed, total_bits=1700, lo=22, hi=48, u=3.5, curves=None):
    pl = planted(seed, total_bits, lo, hi, u)
    N, ps = pl["N"], pl["primes"]
    c = pl["curves"] if curves is None else curves
    r = split(N, pl["B1"], c, "ecm-split-%s" % seed)
    parts = [sum(1 for p in ps if b % p == 0) for b in r["base"]]
    cofp = sum(1 for p in ps if r["cofactor"] % p == 0)
    prod = r["cofactor"]
    for b in r["base"]:
        prod *= b
    return {"seed": seed, "nbits": N.bit_length(), "nprimes": len(ps), "k_main": pl["k_main"],
            "maxbits": ps[-1].bit_length(), "B1": pl["B1"], "curves": c,
            "parts": sorted(parts), "max_part_bits": max([b.bit_length() for b in r["base"]] or [0]),
            "cof_primes": cofp, "cof_bits": r["cofactor"].bit_length() if r["cofactor"] > 1 else 0,
            "product_ok": prod == N, "recorded": r["recorded"], "ops": r["ops"]}


def summarize_cell(rows):
    """rows: list of measure_prime dicts for one (bits, u) cell."""
    n = sum(r["curves"] for r in rows)
    c = sum(r["captures"] for r in rows)
    p = c / n if n else 0.0
    ref = sum(r["curves"] * rho(math.log(r["l"]) / math.log(r["B1"])) for r in rows) / n
    ref12 = sum(r["curves"] * rho(math.log(r["l"] / 12) / math.log(r["B1"])) for r in rows) / n
    return {"curves": n, "captures": c, "degenerate": sum(r["degenerate"] for r in rows),
            "p": p, "se": math.sqrt(p * (1 - p) / n) if n else 0.0,
            "rho_u": ref, "ratio": p / ref, "rho_u12": ref12, "ratio12": p / ref12,
            "primes": len(set(r["l"] for r in rows))}


def replay(seed=1, budget_curves=2000):
    """Small deterministic re-measurement: 32-bit primes at u = 4 (budget_curves curves over
    6 primes) and 5 planted composites of 600 bits (primes of 10..14 bits, u = 3.5)."""
    per = max(1, budget_curves // 6)
    rows = [measure_prime(cell_prime(32, i), 4.0, per, "replay-%s-%d" % (seed, i)) for i in range(6)]
    s = summarize_cell(rows)
    comp = [planted_run("replay-%s-%d" % (seed, i), 600, 10, 14, 3.5) for i in range(5)]
    for c in comp:
        del c["ops"]
    return {"capture_32bit_u4": s, "capture_rows": rows, "planted_600": comp,
            "planted_all_isolated": sum(1 for c in comp if set(c["parts"]) <= {1} and c["cof_primes"] == 0)}


# ----------------------------------------------------------------------------- self test


def selftest():
    assert abs(rho(2) - 0.3069) < 1e-4 and abs(rho(5) / 3.547e-4 - 1) < 1e-3 and abs(rho(8) / 3.232e-8 - 1) < 1e-3
    # ladder against group orders: [#E * #E'] P = O on the x-line for small primes
    rng = random.Random(7)
    for l in (1009, 10007, 65537):
        for _ in range(5):
            cv, g = suyama(rng.randrange(6, l), l)
            if cv is None:
                continue
            X, Z, A24 = cv
            A = (4 * A24 - 2) % l
            if (A * A - 4) % l == 0:
                continue
            s = 0
            e = (l - 1) // 2
            for x in range(l):
                y = pow((x * x * x + A * x * x + x) % l, e, l)
                s += -1 if y == l - 1 else y
            n1, n2 = l + 1 + s, l + 1 - s
            assert n1 % 4 == 0 and n2 % 4 == 0
            assert n1 % 12 == 0 or n2 % 12 == 0, (l, n1, n2)
            for n in (n1, n2):
                while n % 2 == 0:
                    X, Z = xdbl(X, Z, A24, l)
                    n //= 2
                if n > 1:
                    X, Z = ladder(n, X, Z, A24, l)
            assert Z == 0
    # instrumented op count against stage1_ops
    for B1 in (2, 3, 16, 100, 1000):
        k2, sc = plan(B1)
        d = k2
        a = 0
        for m in sc:
            d += 1 + sum(1 for _ in bin(m)[3:-1])
            a += 1 + sum(1 for _ in bin(m)[3:-1])
        o = stage1_ops(B1)
        assert (d, a) == (o["xdbl"], o["xadd"]) and d <= stage1_bits(B1) and a <= stage1_bits(B1)
    # a real split
    ps = [1000003, 2000003, 4000037, 999983, 15485863]
    N = 1
    for p in ps:
        N *= p
    r = split(N, 300, 80, 3)
    prod = r["cofactor"]
    for b in r["base"]:
        prod *= b
    assert prod == N and len(r["base"]) >= 3, r
    # capture() agrees with curve_factors() on N = l
    l = cell_prime(24, 0)
    rr = random.Random(5)
    for _ in range(300):
        s = rr.randrange(6, 1 << 64)
        f, good = curve_factors(l, s, 64)
        c = capture(l, s, 64)
        assert (c == 1) == (good and f == [l]), (s, f, c)
    return True


if __name__ == "__main__":
    import time
    t = time.time()
    print("selftest", selftest(), round(time.time() - t, 2), "s")
    for e in (8, 12, 16, 20):
        b = stage1_bits(1 << e)
        print("stage1_bits(2^%d) = %d  ratio %.5f" % (e, b, b / (1 << e)), stage1_ops(1 << e))
    t = time.time()
    r = replay()
    print(json.dumps({k: r[k] for k in ("capture_32bit_u4", "planted_all_isolated")}), round(time.time() - t, 2), "s")
