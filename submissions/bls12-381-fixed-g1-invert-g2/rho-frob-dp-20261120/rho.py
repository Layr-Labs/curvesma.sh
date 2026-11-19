"""Executable mathematical reference. Python runtime is NOT the scored machine.

Frobenius-class distinguished-point rho walk on GT (classes {Y^(p^i)}, size 12).
In the attack the base is g=e(P,Q0) (protected py_ecc pairing, computed once),
the target is z, elements are kept in the blst tower representation with
canonical Montgomery coefficients, and the output is Q=[x]Q0.

The walker is generic. A group is given by callbacks:
  mul(u, v), power(n) = g^n, tpower(n) = z^n, frob(i, y) = y^(p^i)
  canon(y) -> (i_star, key)
Here key is a class invariant: an integer that is equal for all 12 members.
i_star is the frame: frob(i_star, y) is the canonical member, and for every
member the rule picks that same canonical element (equivariance).
lam = p mod r is the eigenvalue: log(frob(1, y)) = lam * log(y).
key bits: index = key & 31, DP iff (key >> 64) & (2^D - 1) == 0,
fingerprint = (key >> 192) & (2^128 - 1).
Test coins are reproducible, not a success proof.
"""
import math


def sample_scalar(r, coins):
    for _ in range(1024):
        v = coins.getrandbits(r.bit_length())
        if v < r:
            return v
    raise RuntimeError("bounded randomness rejection")


def walk_cap(r, gap):
    """S = floor(31*(isqrt(r//12)+1)/20) + gap >= 1.55*sqrt(r/12) + gap."""
    return (31 * (math.isqrt(r // 12) + 1)) // 20 + gap


def frob_rho(r, lam, target, identity, mul, power, tpower, frob, canon, coins,
             dbits, gap, ring_len=64, max_entries=1 << 21, max_verify=64,
             step_limit=None, stats=None):
    """Returns (x, steps, outcome, escapes). Steps count every multiplication."""
    if target == identity:
        return 0, 0, "identity", 0
    lp = [pow(lam, i, r) for i in range(12)]
    T = [[None] * 32 for _ in range(12)]   # T[i][j] = frob^-i(M_j), labels scaled
    for j in range(32):
        c, d = sample_scalar(r, coins), sample_scalar(r, coins)
        m = mul(power(c), tpower(d))
        for i in range(12):
            s = lp[(12 - i) % 12]
            T[i][j] = (frob((12 - i) % 12, m), s * c % r, s * d % r)
    a, b = sample_scalar(r, coins), sample_scalar(r, coins)
    y = mul(power(a), tpower(b))
    dmask = (1 << dbits) - 1
    cap = walk_cap(r, gap)
    if step_limit is not None:      # experiments only: truncated walk
        cap = min(cap, step_limit)
    stored, ring = {}, []
    verifications = escapes = steps = 0

    def solve(a1, b1, a2, b2):
        den = (b1 - b2) % r
        x = ((a2 - a1) * pow(den, r - 2, r)) % r
        return x if power(x) == target else None

    while steps < cap:
        i, key = canon(y)
        # short-cycle check against the last ring_len states (full key compare)
        hit = next((n for n, e in enumerate(ring) if e[0] == key), None)
        if hit is not None:
            _, yh, ah, bh, ih = ring[hit]
            ca, cb = lp[i] * a % r, lp[i] * b % r
            ha, hb = lp[ih] * ah % r, lp[ih] * bh % r
            if (ca, cb) != (ha, hb):          # genuine relation
                if cb == hb:
                    return None, steps, "degenerate", escapes
                x = solve(ca, cb, ha, hb)
                verifications += 1
                if x is not None:
                    return x, steps, "success", escapes
                if verifications >= max_verify:
                    return None, steps, "verification_failed", escapes
            # fruitless cycle: square the member with minimal key, clear ring
            cyc = ring[hit:] + [(key, y, a, b, i)]
            _, ym, am, bm, _ = min(cyc, key=lambda e: e[0])
            y = mul(ym, ym)
            a, b = 2 * am % r, 2 * bm % r
            ring = []
            escapes += 1
            steps += 1
            continue
        ring.append((key, y, a, b, i))
        if len(ring) > ring_len:
            ring.pop(0)
        if (key >> 64) & dmask == 0:
            fp = (key >> 192) & ((1 << 128) - 1)
            ca, cb = lp[i] * a % r, lp[i] * b % r
            prev = stored.get(fp)
            if prev is None or prev == (ca, cb):
                if prev is None and len(stored) >= max_entries:
                    return None, steps, "table_full", escapes
                stored[fp] = (ca, cb)
            else:
                if prev[1] == cb:
                    return None, steps, "degenerate", escapes
                x = solve(ca, cb, prev[0], prev[1])
                verifications += 1
                if x is not None:
                    return x, steps, "success", escapes
                if verifications >= max_verify:
                    return None, steps, "verification_failed", escapes
                stored[fp] = (ca, cb)
        m, dc, dd = T[i][key & 31]
        if stats is not None:
            stats.setdefault("index", [0] * 32)[key & 31] += 1
            stats.setdefault("frame", [0] * 12)[i] += 1
        y = mul(y, m)
        a, b = (a + dc) % r, (b + dd) % r
        steps += 1
    if stats is not None:
        stats["final"] = (y, a, b)
    return None, steps, "cap", escapes


# ---- Protected polynomial basis Fp[w]/(w^12 - 2 w^6 + 2) ----
def fp12_mul(a, b, p):
    c = [0] * 23
    for i in range(12):
        for j in range(12):
            c[i + j] = (c[i + j] + a[i] * b[j]) % p
    for k in range(22, 11, -1):
        t = (2 * c[k]) % p
        c[k - 6] = (c[k - 6] + t) % p
        c[k - 12] = (c[k - 12] - t) % p
    return tuple(c[:12])


def fp12_pow(a, n, p, mul=None):
    mul = mul or (lambda x, y: fp12_mul(x, y, p))
    out = None
    while n:
        if n & 1:
            out = a if out is None else mul(out, a)
        a = mul(a, a)
        n >>= 1
    return out


# ---- blst tower: Fp2=Fp[u]/(u^2+1), Fp6=Fp2[v]/(v^3-(1+u)), Fp12=Fp6[w]/(w^2-v)
# Flat tower order: index 6*e + 2*k + i is coefficient c_i of the Fp2 value
# multiplying v^k w^e.  Since w^2=v and w^6=1+u, this is the same w as the
# protected basis, so u = w^6 - 1.
def tower_to_poly(t, p):
    out = [0] * 12
    for e in range(2):
        for k in range(3):
            c0, c1 = t[6 * e + 2 * k], t[6 * e + 2 * k + 1]
            n = 2 * k + e
            out[n] = (out[n] + c0 - c1) % p
            out[n + 6] = (out[n + 6] + c1) % p
    return tuple(out)


def poly_to_tower(f, p):
    out = [0] * 12
    for e in range(2):
        for k in range(3):
            n = 2 * k + e
            out[6 * e + 2 * k] = (f[n] + f[n + 6]) % p
            out[6 * e + 2 * k + 1] = f[n + 6] % p
    return tuple(out)


def _f2mul(a, b, p):
    return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)


def _f2add(a, b, p):
    return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)


def _f2xi(a, p):  # multiply by 1+u
    return ((a[0] - a[1]) % p, (a[0] + a[1]) % p)


def _f6mul(a, b, p):
    m = [[_f2mul(a[i], b[j], p) for j in range(3)] for i in range(3)]
    c0 = _f2add(m[0][0], _f2xi(_f2add(m[1][2], m[2][1], p), p), p)
    c1 = _f2add(_f2add(m[0][1], m[1][0], p), _f2xi(m[2][2], p), p)
    c2 = _f2add(_f2add(m[0][2], m[1][1], p), m[2][0], p)
    return (c0, c1, c2)


def _f6add(a, b, p):
    return tuple(_f2add(a[i], b[i], p) for i in range(3))


def _f6v(a, p):  # multiply by v
    return (_f2xi(a[2], p), a[0], a[1])


def tower_mul(x, y, p):
    def split(t):
        return [tuple((t[6 * e + 2 * k], t[6 * e + 2 * k + 1]) for k in range(3))
                for e in range(2)]
    A0, A1 = split(x)
    B0, B1 = split(y)
    c0 = _f6add(_f6mul(A0, B0, p), _f6v(_f6mul(A1, B1, p), p), p)
    c1 = _f6add(_f6mul(A0, B1, p), _f6mul(A1, B0, p), p)
    return tuple(v for half in (c0, c1) for pair in half for v in pair)




# ---- Frobenius on the blst tower and the production class key / frame ----
def _f2pow(a, n, p):
    out = (1, 0)
    while n:
        if n & 1:
            out = _f2mul(out, a, p)
        a = _f2mul(a, a, p)
        n >>= 1
    return out


def frob_constants(p):
    """gamma[i][n] = xi^(n*(p^i-1)/6), xi = 1+u, for the monomial w^n, n=0..5."""
    xi = (1, 1)
    order = p * p - 1  # xi lies in Fp2*, so exponents reduce mod p^2-1
    return [[_f2pow(xi, (n * (p ** i - 1) // 6) % order, p) for n in range(6)]
            for i in range(12)]


def tower_frob(i, t, p, gam):
    """pi^i(Y) = Y^(p^i): coefficient of v^k w^e maps to sigma^i(c)*gamma[i][2k+e]."""
    out = [0] * 12
    for e in range(2):
        for k in range(3):
            c0, c1 = t[6 * e + 2 * k], t[6 * e + 2 * k + 1]
            if i % 2:
                c1 = (-c1) % p
            d = _f2mul((c0, c1), gam[i][2 * k + e], p)
            out[6 * e + 2 * k], out[6 * e + 2 * k + 1] = d
    return tuple(out)


def mont(c, p):
    return (c << 384) % p


def production_canon(t, p, gam):
    """Class key and frame as specified in machine.md.

    key = canonical Montgomery residue of Re(c_{0,0}). This is invariant
    because pi^i maps c00 to sigma^i(c00) and gamma[i][0] = 1.
    frame: for i in 0..5 let V_i = Mont(Re(sigma^i(c_{1,0}) * gamma[i][2]))
    (coefficient of v). Since pi^6 fixes the v coefficient, i0 = argmin V_i.
    Then the w coefficient decides between i0 and i0+6, because pi^6 negates
    it: with W = Mont(Re(sigma^i0(c_{0,1}) * gamma[i0][1])), choose i0 if
    1 <= W <= (p-1)/2, else i0+6.
    """
    key = mont(t[0], p)
    vc = (t[2], t[3])
    wc = (t[6], t[7])
    best = None
    for i in range(6):
        c = vc if i % 2 == 0 else (vc[0], (-vc[1]) % p)
        val = mont(_f2mul(c, gam[i][2], p)[0], p)
        if best is None or val < best[0]:
            best = (val, i)
    i0 = best[1]
    c = wc if i0 % 2 == 0 else (wc[0], (-wc[1]) % p)
    W = mont(_f2mul(c, gam[i0][1], p)[0], p)
    return (i0 if 1 <= W <= (p - 1) // 2 else i0 + 6), key
