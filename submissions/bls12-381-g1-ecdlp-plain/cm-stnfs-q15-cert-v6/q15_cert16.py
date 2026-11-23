"""q15_cert16.py - per-lattice yield certificate run by every descent trial (DESCENT.md section 3).

After reducing L = L_{Q, p_rand}, the node draws S_CERT = 256 uniform primitive half-cells c of B_16(R)/{+-1}
(fresh coins), computes the exact norms and
    X_c = s_f(N_f(c)) * s_g(N_g(c)) / m_c,
s_f, s_g = the node's two-window semismooth lower sums (q15_descent16.ss_lower2, same windows as the plan), m_c = the
exact number of cells of the unit orbit {eps c : eps in O_H^*} inside B_16(R)/{+-1} (Fincke-Pohst over
eps = sum e_k y^k with ||eps c|| <= R, kept iff det(mult by eps) = +-1). Summing X_c over the ball counts every
orbit once, so mu(L) = Nprim * E_c[X_c] * det * lin_keep * ret is the expected number of distinct usable relations
of the ball. Certificate (layer cake): for levels x_1 < .. < x_J, E[X] >= sum_j (x_j - x_{j-1}) P(X >= x_j);
each P(X >= x_j) is replaced by its one-sided Clopper-Pearson lower bound at level (DELTA/D_MAX)/J. For every
lattice, P(the draw passes and mu(L) < MU_CERT) <= DELTA/D_MAX over the node's own coins; a slot makes at most
D_MAX draws, so its sieved lattice violates mu(L) >= MU_CERT with probability <= DELTA in total. Admission rules of the
search (basis screen, v_Q >= 2 rejection, BOUND_CELL cap) are applied: rejected cells get weight 0.
"""
import math
from fractions import Fraction

S_CERT = 256
LEVEL_EXPS = tuple(k / 2.0 for k in range(-8, 13))  # levels x_plan * 2^(k/2), J = 21

_Y8 = [[0] * 8 for _ in range(8)]
for _i in range(7):
    _Y8[_i][_i + 1] = 1
_Y8[7] = [-1, 1, 0, -1, 1, -1, 0, 1]  # y^8 = y^7 - y^5 + y^4 - y^3 + y - 1

# multiplication table of Z[zeta_15] in the power basis: y^i * y^j
_MT = {}


def _ypow(k):
    v = [0] * 8
    if k < 8:
        v[k] = 1
        return v
    v = _ypow(k - 1)
    out = [0] * 8
    for i in range(8):
        if v[i]:
            for j in range(8):
                out[j] += v[i] * _Y8[i][j]
    return out


_YP = [_ypow(k) for k in range(15)]


def _mul(a, b):
    acc = [0] * 15
    for i in range(8):
        if a[i]:
            for j in range(8):
                if b[j]:
                    acc[i + j] += a[i] * b[j]
    out = [0] * 8
    for k in range(15):
        if acc[k]:
            for j in range(8):
                out[j] += acc[k] * _YP[k][j]
    return out


def _det8(M):
    A = [row[:] for row in M]
    prev, sg = 1, 1
    for k in range(7):
        if A[k][k] == 0:
            r = next((i for i in range(k + 1, 8) if A[i][k]), None)
            if r is None:
                return 0
            A[k], A[r] = A[r], A[k]
            sg = -sg
        for i in range(k + 1, 8):
            for j in range(k + 1, 8):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sg * A[7][7]


def is_unit(e):
    M = [_mul(e, [int(i == r) for i in range(8)]) for r in range(8)]
    return abs(_det8(M)) == 1


def _inv_frac(B):
    n = len(B)
    A = [[Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(B)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def lattice_y(B):
    """Multiplication by y in reduced-basis coordinates: Y = B diag(Y8, Y8) B^-1 (integer since L is an O_H-module)."""
    Bi = _inv_frac(B)
    BM = [[sum(B[i][k] * (_Y8[k % 8][j % 8] if (k < 8) == (j < 8) else 0) for k in range(16)) for j in range(16)] for i in range(16)]
    Y = [[sum(BM[i][k] * Bi[k][j] for k in range(16)) for j in range(16)] for i in range(16)]
    assert all(x.denominator == 1 for row in Y for x in row)
    return [[int(x) for x in row] for row in Y]


def orbit_size(x, Y, R2):
    """m_c: number of cells eps*c (eps a unit of O_H) in B_16(R)/{+-1}, by exact Fincke-Pohst in rational arithmetic:
    the form e^T G e = ||sum_k e_k (c Y^k)||^2 has an integer Gram matrix G; its LDL^T factors are kept as Fractions and
    every enumeration bound is decided by exact comparisons, so the enumeration is complete; each e with
    e^T G e <= R2 is kept iff det(multiplication by eps) = +-1."""
    rows = [list(x)]
    for _ in range(7):
        rows.append([sum(rows[-1][i] * Y[i][j] for i in range(16)) for j in range(16)])
    G = [[sum(a * b for a, b in zip(rows[i], rows[j])) for j in range(8)] for i in range(8)]
    # exact LDL^T: G = L D L^T, L unit lower triangular (mu[i][j] = L[i][j], j < i)
    mu = [[Fraction(0)] * 8 for _ in range(8)]
    Dg = [Fraction(0)] * 8
    for i in range(8):
        for j in range(i):
            mu[i][j] = (Fraction(G[i][j]) - sum(mu[i][k] * mu[j][k] * Dg[k] for k in range(j))) / Dg[j]
        Dg[i] = Fraction(G[i][i]) - sum(mu[i][k] * mu[i][k] * Dg[k] for k in range(i))
    # q(e) = sum_i D_i (e_i + sum_{j>i} mu[j][i] e_j)^2 ; enumerate e_7, e_6, ..., e_0
    e = [0] * 8
    cnt = 0
    R2f = Fraction(R2)

    def rec(i, rem):
        nonlocal cnt
        if i < 0:
            if any(e) and is_unit(e):
                cnt += 1
            return
        c = -sum(mu[j][i] * e[j] for j in range(i + 1, 8))
        x0 = math.floor(c)
        x = x0
        while Dg[i] * (x - c) ** 2 <= rem:
            e[i] = x
            rec(i - 1, rem - Dg[i] * (x - c) ** 2)
            x -= 1
        x = x0 + 1
        while Dg[i] * (x - c) ** 2 <= rem:
            e[i] = x
            rec(i - 1, rem - Dg[i] * (x - c) ** 2)
            x += 1
        e[i] = 0

    rec(7, R2f)
    return max(1, cnt // 2)


_POOL = None


def _prim15(p):
    for g in range(2, p):
        z = pow(g, (p - 1) // 15, p)
        if pow(z, 5, p) != 1 and pow(z, 3, p) != 1:
            return z


def prand_pool(V, RT):
    """The declared p_rand pool: the degree-1 prime ideals (p, y - s, x - t) of the f_9 side with p in [2^14, 2^16)
    (p = 1 mod 15, s one of the 8 primitive 15th roots of unity mod p, f_9(t) = 0 at omega = s^5). Stored as
    (p, z, [r(omega_1), r(omega_2)]) with r = number of roots, from deg gcd(f_9, x^p - x); 4,740 ideals in all."""
    global _POOL
    if _POOL is None:
        pool, tot = [], 0
        for p in range(1 << 14, 1 << 16):
            if p % 15 != 1 or not V._is_prime(p):
                continue
            z = _prim15(p)
            rs = []
            for w in (pow(z, 5, p), pow(z, 10, p)):
                f9 = [(-w) % p, 0, 0, 1, 0, 0, (w + 1) % p, 0, 0, 1]
                xp = RT.ppowmod([0, 1], p, f9, p)
                xp += [0] * max(0, 2 - len(xp))
                xp[1] = (xp[1] - 1) % p
                rs.append(len(RT.pgcd(f9, RT._trim(xp), p)) - 1)
            pool.append((p, z, rs))
            tot += 4 * sum(rs)
        _POOL = (pool, tot)
    return _POOL


def draw_prand(V, RT, rng):
    """Uniform draw from the pool: a uniform index among all 4,740 ideals, decoded to (p, s, t)."""
    pool, tot = prand_pool(V, RT)
    k = rng.randrange(tot)
    for p, z, rs in pool:
        n_p = 4 * sum(rs)
        if k >= n_p:
            k -= n_p
            continue
        # the 8 primitive roots s = z^k, k in {1,2,4,7,8,11,13,14}; s^5 = z^5 for k = 1,4,7,13 and z^10 for 2,8,11,14
        for wi, ks in ((0, (1, 4, 7, 13)), (1, (2, 8, 11, 14))):
            blk = 4 * rs[wi]
            if k >= blk:
                k -= blk
                continue
            s = pow(z, ks[k // rs[wi]], p)
            w = pow(s, 5, p)
            f9 = [(-w) % p, 0, 0, 1, 0, 0, (w + 1) % p, 0, 0, 1]
            roots = sorted(RT.roots_Fp(f9, p))
            assert len(roots) == rs[wi]
            return p, s, roots[k % rs[wi]]
    raise AssertionError("pool index out of range")


def lcb_binom(k, n, delta):
    """One-sided Clopper-Pearson lower bound: the largest p with P(Bin(n, p) >= k) < delta (bisection)."""
    if k == 0:
        return 0.0
    lo, hi = 0.0, k / n
    for _ in range(60):
        mid = (lo + hi) / 2
        tail = 1.0 - sum(math.comb(n, i) * mid ** i * (1 - mid) ** (n - i) for i in range(k))
        if tail < delta:
            lo = mid
        else:
            hi = mid
    return lo


def layer_cake_lcb(xs, levels, delta):
    tot, prev = 0.0, 0.0
    for lv in levels:
        k = sum(1 for x in xs if x >= lv)
        tot += (lv - prev) * lcb_binom(k, len(xs), delta / len(levels))
        prev = lv
    return tot


def _crt(a, m, b, n):
    return (a + m * ((b - a) * pow(m, -1, n) % n)) % (m * n)


def draw_delta(delta, d_max=64):
    """Error budget of every certificate draw: delta / D_MAX (a slot makes at most D_MAX draws, so the union bound
    over its draws is <= delta; all draws use the same rule and are identically distributed)."""
    return delta / d_max


def cell_weight(DS, V, BC, Pm, d, side, n, iA, kA, kB, v, q, lq, lp, hw_env=0.625):
    """X_c numerator for one cell, applying the search's admission rules exactly: zero norm, v_Q(N_side) >= 2
    (rejected at exact norm verification) and the BOUND_CELL cap (nf_max, ng_max) on the certified upper bounds
    with q and p_rand removed from their sides. Rejected cells get weight 0."""
    nf, ng = V._norms_fg(v[:8], v[8:])
    if nf == 0 or ng == 0:
        return 0.0
    ns = nf if side == "f" else ng
    if (abs(ns) // q) % q == 0:
        return 0.0
    nf_ub, ng_ub, _, _ = BC.bound_cell_int_fg(v[:8], v[8:])
    if nf_ub < 0.0 or ng_ub < 0.0:
        return 0.0
    if nf_ub + hw_env - (lq if side == "f" else 0.0) - lp > d["nf"] or ng_ub + hw_env - (lq if side == "g" else 0.0) > d["ng"]:
        return 0.0
    lf = math.log2(abs(nf)) - (lq if side == "f" else 0.0) - lp
    lg = math.log2(abs(ng)) - (lq if side == "g" else 0.0)
    return DS.ss_lower2(lf, Pm.LX, float(n), iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC) * \
        DS.ss_lower2(lg, Pm.LX, float(n), iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)


def certify_entry(DS, V, BC, RT, Pm, res, key, nq, npr, seed, first_only=False, cells_only=None):
    """Production campaign for one plan entry: nq ideals Q (bits = bin top edge) x npr independent first draws of
    p_rand; per lattice the basis screen and the certificate of S_CERT cells at the per-draw budget
    draw_delta(DELTA_CERT, D_MAX). Deterministic in (seed, entry). Returns one record per lattice."""
    import random
    (blo, bhi), side = key
    d = res["DET"][key]
    n, iA, kA, kB, R = DS.PLAN16[blo, side]
    nprim = DS.CELLS[R] * 0.9999
    const = d["det"] * d["lin_keep"] * d["ret"]
    xplan = d["mu"] / (nprim * const)
    levels = [xplan * 2.0 ** k for k in LEVEL_EXPS]
    rd, pref = V._build_dp_ball(rmax=R)
    ei = sorted(res["DET"].keys()).index(key)
    rng = random.Random(seed * 1000003 + ei)
    bits = int(math.ceil(bhi))
    out = []
    for iq in range(nq):
        q, s, t = sample_ideal_uniform(V, bits, side, rng)
        for ip in range(npr):
            p, sp, tp = draw_prand(V, RT, rng)
            B, _, _ = V._lll_16d_qcm(V._sq_hnf(q * p, _crt(s, q, sp, p), _crt(t, q, tp, p)))
            screen = BC.screen_special_q_basis(B, q, side, V._norms_fg)
            Y = lattice_y(B)
            lq, lp = math.log2(q), math.log2(p)
            xs, inv_m, n_rej = [], 0.0, 0
            for _ in range(S_CERT):
                c = V._sample_prim_ball_dp(rd, pref, R, rng)
                v = [sum(c[i] * B[i][j] for i in range(16)) for j in range(16)]
                sc = cell_weight(DS, V, BC, Pm, d, side, n, iA, kA, kB, v, q, lq, lp)
                n_rej += sc == 0.0
                m = orbit_size(c, Y, R * R)
                xs.append(sc / m)
                inv_m += 1.0 / m
                if cells_only and len(xs) >= cells_only:
                    return xs
            mu_hat = nprim * const * sum(xs) / S_CERT
            cert = nprim * const * layer_cake_lcb(xs, levels, draw_delta(Pm.DELTA_CERT, Pm.D_MAX)) if screen else 0.0
            out.append([ei, iq, ip, q.bit_length(), round(d["mu"], 4), round(mu_hat, 4), round(cert, 4), round(inv_m / S_CERT, 4), int(screen), n_rej])
            if first_only:
                return (out, xs) if cells_only == 0 else out
    return out


def sample_ideal_uniform(V, bits, side, rng):
    """Uniform degree-1 prime ideal (q, y - s, x - t) of the given side with q of the given bit length: q uniform
    among primes = 1 mod 15 of that length, s uniform among the 8 primitive 15th roots of unity mod q, and the pair
    (q, s) kept with probability k / deg (k = number of roots of the side polynomial at omega = s^5, deg = 3 for g,
    9 for f), then t uniform among the k roots. Every ideal of that bit length has the same probability."""
    deg = 3 if side == "g" else 9
    while True:
        q = rng.randrange(1 << (bits - 1), 1 << bits)
        q = q - (q % 15) + 1
        if q.bit_length() != bits or not V._is_prime(q):
            continue
        z = _prim15(q)
        s = pow(z, rng.choice((1, 2, 4, 7, 8, 11, 13, 14)), q)
        w = pow(s, 5, q)
        poly = [(-V.U_BLS * w * w) % q, 0, 0, 1] if side == "g" else [(-w) % q, 0, 0, 1, 0, 0, (w + 1) % q, 0, 0, 1]
        rts = V._roots_mod_q(poly, q, rng)
        if rts and rng.randrange(deg) < len(rts):
            return q, s, rts[rng.randrange(len(rts))]
