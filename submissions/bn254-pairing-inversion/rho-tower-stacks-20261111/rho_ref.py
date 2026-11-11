"""Executable mathematical reference for the tower-arithmetic, partitioned-stack rho.

Python is NOT the priced machine. machine.md prices the scalar implementation.
This module (standard library only) provides:

* exact Fp2/Fp6/Fp12 tower arithmetic for BN254 with operation counters, and
  the linear basis maps between the protected polynomial basis
  Fp[w]/(w^12-18w^6+82) and the tower Fp2[w]/(w^6-xi), xi=9+u;
* the generic walk (1024-adding multiplicative walk, multiplicative-hash
  partition, K-way partitioned Nivasch stacks, hard step cap), written for
  any group given by (mul, eq-key, power) callbacks, used by the toy runs
  and identical in control flow to the production description in proof.md;
* affine G2 arithmetic for the output Q=[x]Q0.
"""

MASK64 = (1 << 64) - 1


# ---------------------------------------------------------------- counters
class Counter:
    def __init__(self):
        self.mul = 0   # Fp Montgomery products (squares counted as products)
        self.add = 0   # Fp add / sub / neg / double (all priced alike)

    def reset(self):
        self.mul = 0
        self.add = 0


CNT = Counter()
P = 0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47
R = 0x30644E72E131A029B85045B68181585D2833E84879B9709143E1F593F0000001


def fmul(a, b):
    CNT.mul += 1
    return a * b % P


def fadd(a, b):
    CNT.add += 1
    return (a + b) % P


def fsub(a, b):
    CNT.add += 1
    return (a - b) % P


# ---------------------------------------------------------------- Fp2 = Fp[u]/(u^2+1)
def f2add(a, b):
    return (fadd(a[0], b[0]), fadd(a[1], b[1]))


def f2sub(a, b):
    return (fsub(a[0], b[0]), fsub(a[1], b[1]))


def f2mul(a, b):
    # Karatsuba: 3 products, 5 additive operations.
    t0 = fmul(a[0], b[0])
    t1 = fmul(a[1], b[1])
    s = fmul(fadd(a[0], a[1]), fadd(b[0], b[1]))
    return (fsub(t0, t1), fsub(fsub(s, t0), t1))


def f2mul_xi(a):
    # (a0 + a1 u)(9 + u) = (9a0 - a1) + (a0 + 9a1) u ; 9a = 8a + a (3 doublings + 1 add)
    def nine(t):
        d = fadd(t, t)
        d = fadd(d, d)
        d = fadd(d, d)
        return fadd(d, t)
    return (fsub(nine(a[0]), a[1]), fadd(a[0], nine(a[1])))


# ---------------------------------------------------------------- Fp6 = Fp2[v]/(v^3-xi)
def f6add(a, b):
    return tuple(f2add(x, y) for x, y in zip(a, b))


def f6sub(a, b):
    return tuple(f2sub(x, y) for x, y in zip(a, b))


def f6mul(a, b):
    # Karatsuba / Devegili-O hEigeartaigh-Scott-Dahab: 6 Fp2 products.
    a0, a1, a2 = a
    b0, b1, b2 = b
    v0 = f2mul(a0, b0)
    v1 = f2mul(a1, b1)
    v2 = f2mul(a2, b2)
    c0 = f2add(f2mul_xi(f2sub(f2sub(f2mul(f2add(a1, a2), f2add(b1, b2)), v1), v2)), v0)
    c1 = f2add(f2sub(f2sub(f2mul(f2add(a0, a1), f2add(b0, b1)), v0), v1), f2mul_xi(v2))
    c2 = f2add(f2sub(f2sub(f2mul(f2add(a0, a2), f2add(b0, b2)), v0), v2), v1)
    return (c0, c1, c2)


def f6mul_v(a):
    return (f2mul_xi(a[2]), a[0], a[1])


# ---------------------------------------------------------------- Fp12 = Fp6[W]/(W^2-v)
def f12mul(a, b):
    # Karatsuba over Fp6: 3 Fp6 products = 18 Fp2 products = 54 Fp products.
    g0, g1 = a
    h0, h1 = b
    t0 = f6mul(g0, h0)
    t1 = f6mul(g1, h1)
    s = f6mul(f6add(g0, g1), f6add(h0, h1))
    return (f6add(t0, f6mul_v(t1)), f6sub(f6sub(s, t0), t1))


ONE12 = (((1, 0), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0)))


def f12pow(a, n):
    out = ONE12
    while n:
        if n & 1:
            out = f12mul(out, a)
        a = f12mul(a, a)
        n >>= 1
    return out


# ---------------------------------------------------------------- basis maps
def poly_to_tower(c):
    """Protected basis sum c_i w^i  ->  tower ((A0,A2,A4),(A1,A3,A5)), A_k in Fp2.

    Since w^6 = xi = 9+u, c_k w^k + c_{k+6} w^{k+6} = ((c_k + 9 c_{k+6}) + c_{k+6} u) w^k.
    """
    A = [((c[k] + 9 * c[k + 6]) % P, c[k + 6] % P) for k in range(6)]
    return ((A[0], A[2], A[4]), (A[1], A[3], A[5]))


def tower_to_poly(t):
    (A0, A2, A4), (A1, A3, A5) = t
    A = [A0, A1, A2, A3, A4, A5]
    c = [0] * 12
    for k in range(6):
        c[k + 6] = A[k][1]
        c[k] = (A[k][0] - 9 * A[k][1]) % P
    return tuple(c)


def poly_mul(a, b):
    """Protected-basis schoolbook product (independent check of the tower)."""
    c = [0] * 23
    for i in range(12):
        for j in range(12):
            c[i + j] = (c[i + j] + a[i] * b[j]) % P
    for k in range(22, 11, -1):
        c[k - 6] = (c[k - 6] + 18 * c[k]) % P
        c[k - 12] = (c[k - 12] - 82 * c[k]) % P
    return tuple(c[:12])


def tower_words(t):
    """The 12 Fp coordinates of a tower element in fixed order (u-basis of A0..A5)."""
    (A0, A2, A4), (A1, A3, A5) = t
    out = []
    for A in (A0, A1, A2, A3, A4, A5):
        out.extend(A)
    return out


# ---------------------------------------------------------------- generic walk
def sample_scalar(n, coins):
    for _ in range(1024):
        v = coins.getrandbits(n.bit_length())
        if v < n:
            return v
    raise RuntimeError("bounded randomness rejection")


def walk(n, gen_pow, target_pow, mul, key64, order_key, cap, coins,
         log2_adders=10, log2_stacks=8, stack_cap=512):
    """One bounded run of the multiplicative r-adding walk with partitioned Nivasch stacks.

    n          prime group order
    gen_pow(a) g^a ; target_pow(b) z^b
    mul(x,y)   group product
    key64(x)   a fixed 64-bit word of the canonical representation of x
    order_key(x) a total-order key on canonical representations (lexicographic)
    cap        hard bound on walk steps
    Returns (x or None, steps, outcome, stack_peak).
    """
    nadd = 1 << log2_adders
    nst = 1 << log2_stacks
    # Fresh odd 64-bit multipliers (the salt) for the index and class hashes.
    h1 = coins.getrandbits(64) | 1
    h2 = coins.getrandbits(64) | 1
    table = []
    for _ in range(nadd):
        c, d = sample_scalar(n, coins), sample_scalar(n, coins)
        table.append((mul(gen_pow(c), target_pow(d)), c, d))
    a, b = sample_scalar(n, coins), sample_scalar(n, coins)
    y = mul(gen_pow(a), target_pow(b))
    stacks = [[] for _ in range(nst)]
    peak = 0
    for step in range(cap + 1):
        w = key64(y)
        cls = ((w * h2) & MASK64) >> (64 - log2_stacks)
        st = stacks[cls]
        ky = order_key(y)
        # Nivasch: pop strictly larger entries, detect equality, push.
        while st and st[-1][0] > ky:
            st.pop()
        if st and st[-1][0] == ky:
            _, a0, b0, s0 = st[-1]
            den = (b - b0) % n
            if den == 0:
                return None, step, "degenerate", peak
            x = (a0 - a) * pow(den, -1, n) % n
            if gen_pow(x) != target_pow(1):
                return None, step, "verification_failed", peak
            return x, step, "success", peak
        if len(st) >= stack_cap:
            return None, step, "stack_overflow", peak
        st.append((ky, a, b, step))
        peak = max(peak, len(st))
        if step == cap:
            break
        j = ((w * h1) & MASK64) >> (64 - log2_adders)
        m, c, d = table[j]
        y = mul(y, m)
        a = (a + c) % n
        b = (b + d) % n
    return None, cap, "cap", peak


# ---------------------------------------------------------------- G2 output (affine, twist over Fp2)
def _f2inv(a):
    nrm = (a[0] * a[0] + a[1] * a[1]) % P
    ni = pow(nrm, P - 2, P)
    return (a[0] * ni % P, (-a[1]) * ni % P)


def _f2m(a, b):
    return ((a[0] * b[0] - a[1] * b[1]) % P, (a[0] * b[1] + a[1] * b[0]) % P)


def g2_add(Pt, Qt):
    if Pt is None:
        return Qt
    if Qt is None:
        return Pt
    (x1, y1), (x2, y2) = Pt, Qt
    if x1 == x2:
        if y1 == y2 and y1 != (0, 0):
            num = _f2m((3, 0), _f2m(x1, x1))
            lam = _f2m(num, _f2inv(((2 * y1[0]) % P, (2 * y1[1]) % P)))
        else:
            return None
    else:
        lam = _f2m(((y2[0] - y1[0]) % P, (y2[1] - y1[1]) % P),
                   _f2inv(((x2[0] - x1[0]) % P, (x2[1] - x1[1]) % P)))
    l2 = _f2m(lam, lam)
    x3 = ((l2[0] - x1[0] - x2[0]) % P, (l2[1] - x1[1] - x2[1]) % P)
    t = _f2m(lam, ((x1[0] - x3[0]) % P, (x1[1] - x3[1]) % P))
    y3 = ((t[0] - y1[0]) % P, (t[1] - y1[1]) % P)
    return (x3, y3)


def g2_mul(Pt, k):
    out = None
    for bit in bin(k)[2:] if k else "":
        out = g2_add(out, out)
        if bit == "1":
            out = g2_add(out, Pt)
    return out


def g2_on_curve(Pt):
    if Pt is None:
        return True
    x, y = Pt
    b2 = _f2m((3, 0), _f2inv((9, 1)))
    lhs = _f2m(y, y)
    x3 = _f2m(_f2m(x, x), x)
    return lhs == ((x3[0] + b2[0]) % P, (x3[1] + b2[1]) % P)
