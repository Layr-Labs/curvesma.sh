#!/usr/bin/env python3
"""Toy evidence for the additive Cheon integer-walk BSGS (v9b, stdlib only).

Part A  algebra of the production scalar field r and the production curve G1.
Part B  exact |S|/m of the product construction at several X2/d.
Part C  full toy attack: production pseudocode, >=1500 trials, chunk path exercised,
        and the incomplete-addition routine ARMED (raises on P = +-Q).
Part D  the three censuses that v9 only asserted:
        D1 fibre counts per character class over the integer interval (vs 1/d),
        D2 gap histogram between consecutive fibre elements (vs the span/B bound),
        D3 pair multiplicity  b/B = b'/B'  in one (baby, giant) class pair.

Every number printed here is recomputed from the group law; nothing is fitted.
"""
import json
import random
import sys
from math import gcd, isqrt

# ---------------------------------------------------------------- primitives
MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n):
    if n < 2:
        return False
    for q in MR_BASES:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in MR_BASES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factor_small(n, bound=100000):
    fs, q = [], 2
    while q < bound and q * q <= n:
        while n % q == 0:
            fs.append(q)
            n //= q
        q += 1 if q == 2 else 2
    return fs, n


def multiplicative_order(a, r, fac):
    o = r - 1
    for q, e in fac:
        for _ in range(e):
            if o % q == 0 and pow(a, o // q, r) == 1:
                o //= q
            else:
                break
    return o


def sqrt_mod(a, r):
    """Tonelli-Shanks square root mod the prime r (None when a is a non-residue)."""
    a %= r
    if a == 0:
        return 0
    if pow(a, (r - 1) // 2, r) != 1:
        return None
    if r % 4 == 3:
        return pow(a, (r + 1) // 4, r)
    q, s = r - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (r - 1) // 2, r) != r - 1:
        z += 1
    m, c, t, x = s, pow(z, q, r), pow(a, q, r), pow(a, (q + 1) // 2, r)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % r
            i += 1
        b = pow(c, 1 << (m - i - 1), r)
        m, c = i, b * b % r
        t, x = t * c % r, x * b % r
    return x


def dlog_pow2(v, zeta, k, r):
    """e in [0,2^k) with zeta^e = v, for v of order dividing 2^k. Pohlig-Hellman."""
    e = 0
    for i in range(k):
        w = v * pow(zeta, -e, r) % r
        if pow(w, 1 << (k - 1 - i), r) != 1:
            e |= 1 << i
    if pow(zeta, e, r) != v:
        raise AssertionError('dlog_pow2 failed')
    return e


# ---------------------------------------------------------------- toy curve
class Curve:
    """y^2 = x^3 + b over F_p with p = 3 mod 4 (affine, None is the identity).

    incomplete=True models the production routine: it raises instead of taking the
    doubling / opposite branch, so a call that needs the exceptional case is caught.
    """

    def __init__(self, p, b):
        self.p, self.b = p, b

    def on(self, A):
        p, b = self.p, self.b
        return A is None or A[1] * A[1] % p == (A[0] ** 3 + b) % p

    def add(self, A, B, incomplete=False):
        if A is None:
            return B
        if B is None:
            return A
        p = self.p
        x1, y1 = A
        x2, y2 = B
        if x1 == x2:
            if incomplete:
                raise AssertionError('incomplete addition hit P = +-Q')
            if (y1 + y2) % p == 0:
                return None
            lam = 3 * x1 * x1 % p * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)

    def mul(self, k, A):
        R = None
        while k:
            if k & 1:
                R = self.add(R, A)
            A = self.add(A, A)
            k >>= 1
        return R

    def sqrt(self, a):
        s = pow(a % self.p, (self.p + 1) // 4, self.p)
        return s if s * s % self.p == a % self.p else None

    def point(self, x):
        y = self.sqrt((x ** 3 + self.b) % self.p)
        return None if y is None else (x, y)


# ---------------------------------------------------------------- part A
BLS_P = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624'
            '1eabfffeb153ffffb9feffffffffaaab', 16)
BLS_R = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
P0 = (int('17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac58'
          '6c55e83ff97a1aeffb3af00adb22c6bb', 16),
      int('08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3ed'
          'd03cc744a2888ae40caa232946c5e7e1', 16))
BLS_FAC = [(2, 32), (3, 1), (11, 1), (19, 1), (10177, 1), (125527, 1), (859267, 1),
           (906349, 2), (2508409, 1), (2529403, 1), (52437899, 1), (254760293, 2)]
BLS_D = 1 << 20


def curve_g1():
    p = BLS_P

    def add(A, B):
        if A is None:
            return B
        if B is None:
            return A
        x1, y1 = A
        x2, y2 = B
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = 3 * x1 * x1 % p * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)

    def mul(k, A):
        R = None
        while k:
            if k & 1:
                R = add(R, A)
            A = add(A, A)
            k >>= 1
        return R

    return add, mul


def part_a():
    out = {}
    r, d = BLS_R, BLS_D
    prod = 1
    for q, e in BLS_FAC:
        prod *= q ** e
    out['r_prime'] = is_prime(r)
    out['r_minus_1_factorisation'] = prod == r - 1 and all(is_prime(q) for q, _ in BLS_FAC)
    out['v2_r_minus_1'] = 32
    out['d_divides_r_minus_1'] = (r - 1) % d == 0
    m = (r - 1) // d
    out['m_log2_floor'] = m.bit_length() - 1
    out['v2_m'] = (m & -m).bit_length() - 1

    g = 2
    while pow(g, (r - 1) // 2, r) == 1 or pow(g, (r - 1) // 3, r) == 1:
        g += 1
    zeta = pow(g, m, r)
    out['zeta_has_order_d'] = pow(zeta, d, r) == 1 and pow(zeta, d // 2, r) != 1
    rng = random.Random(20261023)
    out['dth_powers_are_mth_roots'] = all(
        pow(pow(rng.randrange(2, r), d, r), m, r) == 1 for _ in range(6))
    out['two_not_in_H'] = pow(2, m, r) != 1
    out['order_of_2_mod_r'] = multiplicative_order(2, r, BLS_FAC)
    out['order_of_2_does_not_divide_m'] = m % out['order_of_2_mod_r'] != 0
    out['order_of_2_pow_m'] = multiplicative_order(pow(2, m, r), r, BLS_FAC)

    def chi(x):
        return dlog_pow2(pow(x % r, m, r), zeta, 20, r)

    add_ok = True
    kern_ok = True
    for _ in range(6):
        x = rng.randrange(2, r)
        y = rng.randrange(2, r)
        add_ok &= chi(x * y % r) == (chi(x) + chi(y)) % d
        kern_ok &= (chi(x) == 0) == (pow(x, m, r) == 1)
    out['chi_additive'] = add_ok
    out['chi_kernel_is_H'] = kern_ok
    out['chi_vanishes_on_dth_powers'] = chi(pow(3, d, r)) == 0
    out['chi_not_identically_zero'] = chi(2) != 0 or chi(3) != 0

    # L1: a*b = a'*b' with a,a' prime > 2^70 and 1 <= b,b' <= 2^70 forces (a,b) = (a',b')
    ok = True
    for _ in range(8):
        a = rng.randrange(2 ** 75, 2 ** 76)
        a2 = rng.randrange(2 ** 75, 2 ** 76)
        b = rng.randrange(1, 2 ** 70)
        b2 = rng.randrange(1, 2 ** 70)
        ok &= (a * b == a2 * b2) == (a == a2 and b == b2)
    out['L1_product_injective_shapes'] = ok

    add, mul = curve_g1()
    out['P0_on_curve'] = (P0[1] * P0[1] - P0[0] ** 3 - 4) % BLS_P == 0
    out['P0_order_r'] = mul(r, P0) is None
    beta = pow(5, (BLS_P - 1) // 3, BLS_P)
    out['beta_cube_root_of_unity'] = pow(beta, 3, BLS_P) == 1 and beta != 1
    lam = pow(g, (r - 1) // 3, r)
    out['lambda_cube_root_of_unity'] = pow(lam, 3, r) == 1 and lam != 1
    Q = mul(rng.randrange(2, r), P0)
    out['glv_map_matches_lambda'] = (beta * Q[0] % BLS_P, Q[1]) == mul(lam, Q)
    orbit = []
    for u in (1, lam, lam * lam % r):
        for sgn in (1, -1):
            orbit.append(mul(u * sgn % r, Q))
    out['orbit_six_distinct'] = len(set(orbit)) == 6
    out['orbit_two_raw_y_values'] = len({X[1] for X in orbit}) == 2
    out['orbit_one_canonical_y_value'] = len({min(X[1], BLS_P - X[1]) for X in orbit}) == 1
    out['orbit_three_distinct_x_values'] = len({X[0] for X in orbit}) == 3
    return out


# ---------------------------------------------------------------- toy field
TOY_P = 268520611
TOY_R = 67122049
TOY_H = 4
TOY_D = 1 << 6


def toy_curve():
    """Sextic twist y^2 = x^3 + b with 4r | ... : #E = TOY_H * TOY_R, plus an order-r point."""
    order = TOY_H * TOY_R
    for b in range(1, 60):
        C = Curve(TOY_P, b)
        pts = []
        ok = True
        for x in range(1, 90):
            Q = C.point(x)
            if Q is None:
                continue
            if C.mul(order, Q) is not None:       # #E does not divide 4r
                ok = False
                break
            pts.append(Q)
        if not ok or len(pts) < 4:
            continue
        for Q in pts:
            if C.mul(TOY_R, Q) is None:           # order r
                return C, Q
            A = C.mul(TOY_H // 2, Q)              # 2r -> r
            if A is not None and C.mul(TOY_R, A) is None:
                return C, A
    raise AssertionError('no twist found')


def toy_setup(seed=7):
    C, A = toy_curve()
    r, d = TOY_R, TOY_D
    m = (r - 1) // d
    fs, cof = factor_small(r - 1)
    if cof != 1:
        assert is_prime(cof)
        fs = fs + [cof]
    fac = [(q, fs.count(q)) for q in sorted(set(fs))]
    prod = 1
    for q, e in fac:
        prod *= q ** e
    assert prod == r - 1 and all(is_prime(q) for q, _ in fac)
    g = 2
    while pow(g, (r - 1) // 2, r) == 1 or pow(g, (r - 1) // 3, r) == 1:
        g += 1
    zeta = pow(g, m, r)
    assert pow(zeta, d, r) == 1 and pow(zeta, d // 2, r) != 1
    lam = pow(g, (r - 1) // 3, r)
    assert pow(lam, 3, r) == 1 and lam != 1
    omega = pow(g % TOY_P, (TOY_P - 1) // 3, TOY_P)
    assert pow(omega, 3, TOY_P) == 1 and omega != 1
    rng = random.Random(seed)
    while True:
        Q = C.point(rng.randrange(2, TOY_P))
        if Q is None:
            continue
        A = C.mul(TOY_H // 2, Q)          # [2]Q kills the order-2 part
        if A is not None and C.mul(r, A) is None:
            break
    if C.mul(lam, A) != (omega * A[0] % TOY_P, A[1]):
        lam = lam * lam % r
        assert C.mul(lam, A) == (omega * A[0] % TOY_P, A[1])
    return dict(C=C, p=TOY_P, r=r, d=d, m=m, zeta=zeta, lam=lam, g=g, P0=A, fac=fac,
                two_in_H=pow(2, m, r) == 1)


def chi_of(env, x):
    return dlog_pow2(pow(x % env['r'], env['m'], env['r']), env['zeta'], 6, env['r'])


def class_lists(env, X2):
    d = env['d']
    buckets = [[] for _ in range(d)]
    for b in range(1, X2 + 1):
        buckets[chi_of(env, b)].append(b)
    return buckets


def primes_in(lo, hi, need):
    out, n = [], lo + 1
    if n % 2 == 0:
        n += 1
    while n <= hi and len(out) < need:
        if is_prime(n):
            out.append(n)
        n += 2
    if len(out) < need:
        raise AssertionError('prime range too small')
    return out


# ---------------------------------------------------------------- part B
def product_coverage(env, ratio, per_side):
    """exact |S| for S = {+-lambda^i x/y} over the product sets."""
    X2 = env['d'] * ratio
    d, r, lam = env['d'], env['r'], env['lam']
    buckets = class_lists(env, X2)
    L = len(buckets[0])
    npr = (per_side + L - 1) // L
    primes = primes_in(X2 + 1, X2 * 400 + 40000, 4 * npr)
    baby, giant = [], []
    for a in primes[:npr]:
        for b in buckets[(-chi_of(env, a)) % d][:max(0, per_side - len(baby))]:
            baby.append(a * b % r)
    for a in primes[npr:npr + npr]:
        for b in buckets[(-chi_of(env, a)) % d][:max(0, per_side - len(giant))]:
            giant.append(a * b % r)
    S = set()
    for y in giant:
        iy = pow(y, -1, r)
        for x in baby:
            t = x * iy % r
            for u in (1, lam, lam * lam % r):
                S.add(u * t % r)
                S.add(-u * t % r)
    pairs = len(baby) * len(giant) * 6
    return dict(X2_over_d=ratio, X2=X2, N=len(baby), S=len(S), pairs=pairs,
                coverage=len(S) / env['m'], model=1 - pow(2.718281828459045, -pairs / env['m']),
                duplicate_ppm=(pairs - len(S)) * 1000000 // pairs)


def part_b(env):
    target = isqrt(env['m'] * 52 // 600)     # 6N^2/m = 0.52
    return [product_coverage(env, ratio, target) for ratio in (4, 16, 64, 256)]


# ---------------------------------------------------------------- part D
def fiber_census(env, ratio):
    """D1: counts per character class over the INTEGER interval [1, X2] vs X2/d.

    The fibre of -chi(a) is a set of integers, not a subgroup: its size is whatever
    the integers give.  Report the exact worst absolute deviation from the mean X2/d
    and the normalised form dev^2 * mean, which exposes the 1/sqrt(L) sampling scale.
    """
    d = env['d']
    X2 = d * ratio
    buckets = class_lists(env, X2)
    counts = [len(b) for b in buckets]
    worst = max(abs(c * d - X2) for c in counts)
    used = [c for c in counts if c > 0]
    return dict(X2=X2, ratio=ratio, mean_num=X2, mean_den=d,
                counts_min=min(counts), counts_max=max(counts),
                counts_nonempty_min=min(used) if used else 0,
                empty_classes=sum(1 for c in counts if c == 0),
                class_count=len(counts),
                max_abs_dev_num=worst, max_abs_dev_den=X2,
                dev_sq_times_mean=[worst * worst, X2 * d])


def part_d1(env):
    return [fiber_census(env, ratio) for ratio in (4, 16, 64, 256, 1024)]


def gap_census(env, ratio, bounds):
    """D2: gap histogram between consecutive integers of one fibre."""
    d = env['d']
    X2 = d * ratio
    buckets = class_lists(env, X2)
    bounds = list(bounds)
    over = {b: 0 for b in bounds}
    per_class = {b: 0 for b in bounds}
    gaps = 0
    max_gap = 0
    total_span = 0
    for cl in buckets:
        here = {b: 0 for b in bounds}
        for i in range(1, len(cl)):
            g = cl[i] - cl[i - 1]
            gaps += 1
            total_span += g
            if g > max_gap:
                max_gap = g
            for b in bounds:
                if g > b:
                    over[b] += 1
                    here[b] += 1
        for b in bounds:
            per_class[b] = max(per_class[b], here[b])
    # the deterministic span bound is PER CLASS (= per prime in production):
    # a gap above B consumes more than B of the class span, and the span is <= X2.
    return dict(X2=X2, ratio=ratio, gaps=gaps, max_gap=max_gap,
                span=total_span,
                over_bound=[[b, over[b]] for b in bounds],
                per_class_max_over=[[b, per_class[b]] for b in bounds],
                det_bound_per_class=[[b, X2 // b] for b in bounds],
                det_bound_total=[[b, (X2 // b) * d] for b in bounds])


def part_d2(env, bounds):
    return [gap_census(env, ratio, bounds) for ratio in (4, 16, 64, 256)]


def pair_multiplicity(env, ratio, a_list, A_list):
    """D3: fraction of (b,B) pairs in one (baby, giant) class pair whose reduced
    ratio b/B repeats.  Exact count, no model."""
    d = env['d']
    X2 = d * ratio
    buckets = class_lists(env, X2)
    out = []
    for a, A in zip(a_list, A_list):
        fb = buckets[(-chi_of(env, a)) % d]
        fg = buckets[(-chi_of(env, A)) % d]
        if not fb or not fg:
            continue
        mult = {}
        for b in fb:
            for B in fg:
                g = gcd(b, B)
                key = (b // g, B // g)
                mult[key] = mult.get(key, 0) + 1
        pairs = len(fb) * len(fg)
        loss = pairs - len(mult)          # pairs absorbed by a repeated ratio
        offdiag = sum(c * (c - 1) // 2 for c in mult.values())
        out.append(dict(ratio=ratio, L_baby=len(fb), L_giant=len(fg), pairs=pairs,
                        distinct_ratios=len(mult), loss=loss,
                        offdiag_pairs=offdiag,
                        max_multiplicity=max(mult.values()),
                        phi_num=loss, phi_den=pairs))
    return out


def part_d3(env, ratios):
    d = env['d']
    res = []
    for ratio in ratios:
        X2 = d * ratio
        primes = primes_in(X2 + 1, X2 * 400 + 40000, 6)
        res.extend(pair_multiplicity(env, ratio, primes[:3], primes[3:6]))
    return res


# ---------------------------------------------------------------- part C
def run_trials(env, trials=1500, B=8, ratio=64, k_num=13, k_den=25, seed=99, trace=None):
    """Production pseudocode on the toy curve.

    N per side is fixed by 6N^2/m = k_num/k_den.  B is the walk table bound, kept
    small so the overflow-chunk path runs constantly.  The single-digit step runs
    through the ARMED incomplete routine; multi-digit chunk steps run through the
    complete routine, because a chunk summand is a sum, not a fibre product.
    """
    C, r, d, m = env['C'], env['r'], env['d'], env['m']
    lam, P0, zeta = env['lam'], env['P0'], env['zeta']
    X2 = d * ratio
    N = isqrt(m * k_num // (k_den * 6))
    buckets = class_lists(env, X2)
    L = len(buckets[0])
    npr = (N + L - 1) // L
    primes = primes_in(X2 + 1, X2 * 400 + 40000, 2 * npr + 2)
    baby_pr, giant_pr = primes[:npr], primes[npr:npr + npr]
    rng = random.Random(seed)

    def make_tabs(base):
        """tabs[e][k] = [k*B^e](base) for e below the digit count needed by X2."""
        tabs, step = [], base
        levels = 1
        while B ** levels <= X2:
            levels += 1
        for _ in range(levels):
            row = [None] * (B + 1)
            row[1] = step
            for k in range(2, B + 1):
                row[k] = C.add(row[k - 1], step)
            tabs.append(row)
            step = row[B]
        return tabs

    def key(A):
        y = A[1]
        return y if y <= TOY_P - y else TOY_P - y

    def walk(cur, gap, tabs, stats):
        """add [gap](base) to cur, digits base B; returns (point, extra adds)."""
        digits = []
        g = gap
        while g:
            digits.append(g % B)
            g //= B
        nz = [e for e, dg in enumerate(digits) if dg]
        extra = 0
        single = (len(nz) == 1 and nz[0] == 0)
        for e in nz:
            if single:
                try:
                    cur = C.add(cur, tabs[e][digits[e]], incomplete=True)
                    stats['incomplete_steps'] += 1
                except AssertionError:
                    stats['incomplete_trips'] += 1
                    cur = C.add(cur, tabs[e][digits[e]])
            else:
                cur = C.add(cur, tabs[e][digits[e]])
                stats['complete_steps'] += 1
                if e:
                    extra += 1
        return cur, extra

    succ = wrong = false_pos = 0
    stats = dict(incomplete_steps=0, incomplete_trips=0, complete_steps=0)
    walk_adds = chunk_adds = overflow_gaps = cand_tests = 0
    for _ in range(trials):
        s = rng.randrange(2, r)
        T = C.mul(s, P0)
        U = C.mul(pow(s, d, r), P0)
        table = {}
        walked = 0
        for a in baby_pr:
            if walked >= N:
                break
            tabs = make_tabs(C.mul(a, P0))
            cl = buckets[(-chi_of(env, a)) % d]
            cur = None
            prev = 0
            for b in cl:
                gap = b - prev
                if gap > B:
                    overflow_gaps += 1
                cur, extra = walk(cur, gap, tabs, stats)
                walk_adds += 1 + extra
                chunk_adds += extra
                table.setdefault(key(cur), []).append((a, b))
                prev = b
                walked += 1
                if walked >= N:
                    break
        found = None
        walked = 0
        for a in giant_pr:
            if walked >= N or found is not None:
                break
            tabs = make_tabs(C.mul(a, U))
            cl = buckets[(-chi_of(env, a)) % d]
            cur = None
            prev = 0
            for b in cl:
                gap = b - prev
                if gap > B:
                    overflow_gaps += 1
                cur, extra = walk(cur, gap, tabs, stats)
                walk_adds += 1 + extra
                chunk_adds += extra
                prev = b
                walked += 1
                row = table.get(key(cur))
                if row:
                    iy = pow(a * b % r, -1, r)
                    for (aa, bb) in row:
                        x = aa * bb % r
                        for u in (1, lam, lam * lam % r):
                            for cand in (x * u % r * iy % r, (-x) * u % r * iy % r):
                                cand_tests += 1
                                if C.mul(cand, P0) == U:
                                    found = cand
                if found is not None or walked >= N:
                    break
        if trace is not None:
            trace.append(found is not None)
        if found is None:
            continue
        rho = found
        for _k in range(6):                       # d = 2^6 successive square roots
            rho = sqrt_mod(rho, r)
            assert rho is not None
        assert pow(rho, d, r) == found
        zinv = pow(zeta, -1, r)
        cand, ok = rho, False
        for _j in range(d):
            if C.mul(cand, P0) == T:
                ok = cand == s
                break
            cand = cand * zinv % r
        if ok:
            succ += 1
        else:
            false_pos += 1
    model = 1 - pow(2.718281828459045, -(6 * N * N) / m)
    sigma = (model * (1 - model) / trials) ** 0.5
    return dict(trials=trials, successes=succ, success_rate=succ / trials,
                model=model, sigma=sigma, deviation_sigma=(succ / trials - model) / sigma,
                false_positives=false_pos, wrong_outputs=wrong, N=N, table_bound=B,
                walk_elements=N, X2_over_d=ratio, K=(6 * N * N) / m,
                walk_additions=walk_adds, chunk_additions=chunk_adds,
                overflow_gaps=overflow_gaps, candidate_tests=cand_tests,
                incomplete_steps=stats['incomplete_steps'],
                incomplete_trips=stats['incomplete_trips'],
                complete_steps=stats['complete_steps'])


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    out = {}
    if which in ('all', 'partA'):
        out['part_a'] = part_a()
    if which in ('all', 'partB', 'partC', 'partD'):
        env = toy_setup()
        out['toy'] = dict(p=env['p'], r=env['r'], d=env['d'], m=env['m'],
                          b=env['C'].b, P0=list(env['P0']), lam=env['lam'],
                          two_in_H=env['two_in_H'])
        if which in ('all', 'partB'):
            out['coverage'] = part_b(env)
        if which in ('all', 'partD'):
            out['census'] = dict(fiber=part_d1(env),
                                 gap=part_d2(env, (env['d'], env['d'] * 4,
                                                   env['d'] * 16)),
                                 pair=part_d3(env, (64, 256)))
        if which in ('all', 'partC'):
            attacks = []
            for ratio, trials in ((64, 1500), (256, 1500)):
                res = run_trials(env, trials=trials, ratio=ratio)
                cov = product_coverage(env, ratio, res['N'])
                res['exact_coverage'] = cov['coverage']
                res['exact_pairs_over_m'] = cov['pairs'] / env['m']
                res['deviation_vs_exact_sigma'] = (
                    (res['success_rate'] - cov['coverage'])
                    / ((cov['coverage'] * (1 - cov['coverage']) / trials) ** 0.5))
                attacks.append(res)
            out['attack'] = attacks
    print(json.dumps(out, indent=1, sort_keys=True))


if __name__ == '__main__':
    main()
