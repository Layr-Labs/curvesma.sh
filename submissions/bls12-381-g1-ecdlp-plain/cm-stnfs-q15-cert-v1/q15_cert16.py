"""q15_cert16.py - per-lattice yield certificate run by every descent trial (DESCENT.md section 3).

After reducing L = L_{Q, p_rand}, the node draws S_CERT = 256 uniform primitive half-cells c of B_16(R)/{+-1}
(fresh coins), computes the exact norms and
    X_c = s_f(N_f(c)) * s_g(N_g(c)) / m_c,
s_f, s_g = the node's two-window semismooth lower sums (q15_descent16.ss_lower2, same windows as the plan), m_c = the
exact number of cells of the unit orbit {eps c : eps in O_H^*} inside B_16(R)/{+-1} (Fincke-Pohst over
eps = sum e_k y^k with ||eps c|| <= R, kept iff det(mult by eps) = +-1). Summing X_c over the ball counts every
orbit once, so mu(L) = Nprim * E_c[X_c] * det * lin_keep * ret is the expected number of distinct usable relations
of the ball. Certificate (layer cake): for levels x_1 < .. < x_J, E[X] >= sum_j (x_j - x_{j-1}) P(X >= x_j);
each P(X >= x_j) is replaced by its one-sided Clopper-Pearson lower bound at level DELTA/J, so
    P(certificate passes and mu(L) < MU_CERT) <= DELTA
over the node's own coins, for every lattice L. A failing lattice is discarded and p_rand redrawn.
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
    """m_c: number of cells eps*c (eps a unit) in B_16(R)/{+-1}; e^T G e = ||sum e_k (c Y^k)||^2 <= R2."""
    rows = [list(x)]
    for _ in range(7):
        rows.append([sum(rows[-1][i] * Y[i][j] for i in range(16)) for j in range(16)])
    G = [[float(sum(a * b for a, b in zip(rows[i], rows[j]))) for j in range(8)] for i in range(8)]
    L = [[0.0] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(i + 1):
            s = G[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            L[i][j] = math.sqrt(s) if i == j else s / L[j][j]
    e = [0] * 8
    cnt = 0

    def rec(i, rem):
        nonlocal cnt
        if i < 0:
            if any(e) and is_unit(e):
                cnt += 1
            return
        c = sum(L[j][i] * e[j] for j in range(i + 1, 8))
        r = math.sqrt(max(rem, 0.0)) / L[i][i]
        ce = -c / L[i][i]
        for v in range(math.ceil(ce - r - 1e-9), math.floor(ce + r + 1e-9) + 1):
            t = (L[i][i] * v + c) ** 2
            if t <= rem * (1 + 1e-12) + 1e-9:
                e[i] = v
                rec(i - 1, rem - t)
        e[i] = 0

    rec(7, float(R2))
    return max(1, cnt // 2)


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


def certify_entry(DS, V, Pm, res, key, nq, npr, seed, first_only=False):
    """Production campaign for one plan entry: nq ideals Q (bits = bin top edge) x npr fresh p_rand draws; per
    lattice the certificate of S_CERT cells. Deterministic in (seed, entry). Returns one record per lattice."""
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
        q, s, t = V._sample_ideal(bits, side, rng)
        for ip in range(npr):
            p, sp, tp = V._sample_ideal(rng.choice((15, 16)), "f", rng)
            B, _, _ = V._lll_16d_qcm(V._sq_hnf(q * p, _crt(s, q, sp, p), _crt(t, q, tp, p)))
            Y = lattice_y(B)
            lq, lp = math.log2(q), math.log2(p)
            xs, inv_m = [], 0.0
            for _ in range(S_CERT):
                c = V._sample_prim_ball_dp(rd, pref, R, rng)
                v = [sum(c[i] * B[i][j] for i in range(16)) for j in range(16)]
                nf, ng = V._norms_fg(v[:8], v[8:])
                lf = math.log2(abs(nf)) - (lq if side == "f" else 0.0) - lp
                lg = math.log2(abs(ng)) - (lq if side == "g" else 0.0)
                sc = DS.ss_lower2(lf, Pm.LX, float(n), iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC) * \
                    DS.ss_lower2(lg, Pm.LX, float(n), iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
                m = orbit_size(c, Y, R * R)
                xs.append(sc / m)
                inv_m += 1.0 / m
            mu_hat = nprim * const * sum(xs) / S_CERT
            cert = nprim * const * layer_cake_lcb(xs, levels, Pm.DELTA_CERT)
            out.append([ei, iq, ip, q.bit_length(), round(d["mu"], 4), round(mu_hat, 4), round(cert, 4), round(inv_m / S_CERT, 4)])
            if first_only:
                return out
    return out
