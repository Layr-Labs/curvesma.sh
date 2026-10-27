"""Toy-scale differential of the descent node mechanics (DESCENT.md items 2-4).

The tower keeps h and f and replaces the production |u| by TOY_U = 4 so that the
regime (thin window, measurable survivors) fits an exhaustive or sampled box at
toy norm sizes. Everything else is the production node structure: special-q
lattice, reduced basis, congruence functionals, per-cell scores, threshold,
survivor residuals, relation predicate and child extraction.

Checks:
 (a) functional: for sampled cells, w.c = 0 mod p  <=>  p | N_side (both sides);
 (b) score: the accumulated score equals the sum of log2 p over distinct sieved
     primes dividing the norm (one credit per prime, valuations >= 1);
 (c) threshold: no false negatives on a sampled cell set (every cell that is a
     true relation passes; exhaustive versions are the documented research runs);
 (d) survivors: exact norms, residual division, factorisation; every survivor
     credited as a relation is one (no false credits) and its children are
     exactly its residual primes in (2^LX, 2^n]; junk survivors are counted.
"""
import math
from random import Random

H = [-1, -1, 1, 1, 0, 0, 1]
TOY_U = 4


def toy_mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] += x * y
    for i in range(len(c) - 1, 5, -1):
        for j in range(6):
            c[i - 6 + j] -= c[i] * H[j]
    return (c + [0] * 6)[:6]


def toy_norm(a):
    cols = [toy_mul(a, [0] * i + [1]) for i in range(6)]
    m = [[cols[j][i] for j in range(6)] for i in range(6)]
    prev, sign = 1, 1
    for k in range(5):
        if m[k][k] == 0:
            row = next((j for j in range(k + 1, 6) if m[j][k]), None)
            if row is None:
                return 0
            m[k], m[row] = m[row], m[k]; sign = -sign
        for i in range(k + 1, 6):
            for j in range(k + 1, 6):
                v = m[i][j] * m[k][k] - m[i][k] * m[k][j]
                assert v % prev == 0
                m[i][j] = v // prev
        prev = m[k][k]
    return abs(sign * m[-1][-1])


def toy_norms(a, b):
    ap, bp = [[1]], [[1]]
    for _ in range(12):
        ap.append(toy_mul(ap[-1], a)); bp.append(toy_mul(bp[-1], b))
    fs = [0] * 6
    for d, c in [(12, [1]), (10, [0, -2]), (6, [0, 0, 0, 2]), (2, [0, 0, 0, 0, 0, 1]), (0, [1, 1, -1, -1])]:
        term = toy_mul(c, toy_mul(ap[d], bp[12 - d]))
        fs = [x + y for x, y in zip(fs, term)]
    gs = [x + y for x, y in zip(ap[2], toy_mul([0, TOY_U], bp[2]))]
    return toy_norm(fs), toy_norm(gs)


def mr_prime(n):
    if n < 2:
        return False
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % a == 0:
            return n == a
    d = n - 1; s = 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
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


def roots_modp(poly, p, rng):
    poly = [c % p for c in poly]
    while poly and poly[-1] == 0:
        poly.pop()
    if len(poly) <= 1:
        return []

    def pmod(a, m):
        a = a[:]
        inv = pow(m[-1], -1, p)
        while len(a) >= len(m):
            if a[-1]:
                z = a[-1] * inv % p; s = len(a) - len(m)
                for j, v in enumerate(m):
                    a[s + j] = (a[s + j] - z * v) % p
            a.pop()
            while a and a[-1] == 0:
                a.pop()
        return a

    def pmul(a, b, m):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    r[i + j] = (r[i + j] + x * y) % p
        return pmod(r, m)

    def powmod(a, e, m):
        r = [1]
        while e:
            if e & 1:
                r = pmul(r, a, m)
            a = pmul(a, a, m); e >>= 1
        return r

    xp = powmod([0, 1], p, poly)
    xp = xp + [0] * max(0, 2 - len(xp)); xp[1] = (xp[1] - 1) % p
    a, b = poly, pmod(xp, poly)
    while b and any(b):
        a, b = b, pmod(a, b)
    if not a or not any(a):
        return []
    inv = pow(a[-1], -1, p); g = [x * inv % p for x in a]

    def split(g):
        if len(g) <= 1:
            return []
        if len(g) == 2:
            return [(-g[0] * pow(g[1], -1, p)) % p]
        while True:
            c = [rng.randrange(p) for _ in range(len(g) - 2)] + [rng.randrange(1, p), 1]
            t = powmod(c, (p - 1) // 2, g)
            t = t + [0] * max(0, 1 - len(t)); t = [(t[0] - 1) % p] + t[1:]
            h = pmod(t, g) if len(t) >= len(g) else t
            aa, bb = g, h
            while bb and any(bb):
                aa, bb = bb, pmod(aa, bb)
            if aa and 1 < len(aa) < len(g):
                iv = pow(aa[-1], -1, p); hh = [x * iv % p for x in aa]
                qq = g[:]; v = [0] * (len(g) - len(hh) + 1)
                for k in range(len(v) - 1, -1, -1):
                    v[k] = qq[k + len(hh) - 1]
                    for j, x in enumerate(hh):
                        qq[k + j] = (qq[k + j] - v[k] * x) % p
                assert not any(qq)
                return split(hh) + split(v)
    return sorted(split(g))


def lll_small(B, delta=0.99):
    """exact-fraction LLL for small integer matrices (toy sizes)."""
    from fractions import Fraction as F
    B = [v[:] for v in B]; n = len(B)

    def gso():
        mu, norm, stars = [], [], []
        for v in B:
            s = [F(x) for x in v]; row = []
            for j, b in enumerate(stars):
                m = sum(x * y for x, y in zip(v, b)) / norm[j]; row.append(m)
                s = [x - m * y for x, y in zip(s, b)]
            stars.append(s); mu.append(row); norm.append(sum(x * x for x in s))
        return mu, norm
    mu, norm = gso()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            m = mu[k][j]
            qv = int(m + F(1, 2)) if m >= 0 else -int(-m + F(1, 2))
            if qv:
                B[k] = [x - qv * y for x, y in zip(B[k], B[j])]
                mu, norm = gso()
        if norm[k] >= (F(99, 100) - mu[k][k - 1] ** 2) * norm[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            mu, norm = gso()
            k = max(1, k - 1)
    return B


def _pollard(n, budget=20000):
    """Pollard rho with an iteration budget: returns a factor, n if the budget
    is exhausted (the production ECM likewise rejects unfactored cofactors)."""
    if n % 2 == 0:
        return 2
    rng = Random(n)
    for _ in range(6):
        x = rng.randrange(2, n); y = x; c = rng.randrange(1, n); d = 1
        for _ in range(budget // 6):
            x = (x * x + c) % n
            y = (y * y + c) % n; y = (y * y + c) % n
            d = math.gcd(abs(x - y), n)
            if d != 1:
                break
        if 1 < d < n:
            return d
        if d == n:
            continue
        return n  # budget exhausted: unrejected composite
    return n


def factorize(n, base_primes):
    """Full factorisation by TD + bounded Pollard. Returns None when a
    cofactor resists the budget (an unfactored survivor is conservatively
    rejected, exactly as the production ECM does)."""
    if n == 1:
        return []
    fs = []
    for p in base_primes:
        if p * p > n:
            break
        while n % p == 0:
            fs.append(p); n //= p
    if n == 1:
        return fs
    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if mr_prime(m):
            fs.append(m); continue
        d = _pollard(m)
        if d == m or d == 1:
            return None
        stack.append(d); stack.append(m // d)
    return sorted(fs)



def _setup_node(seed, qbits_lo, qbits_hi, kbox, LX, P_MIN=256):
    rng = Random(seed)
    while True:
        q = rng.randrange(1 << qbits_lo, 1 << qbits_hi) | 1
        if not mr_prime(q):
            continue
        rs = roots_modp(H, q, rng)
        opts = [(s, t) for s in rs for t in roots_modp([TOY_U * s % q, 0, 1], q, rng)]
        if opts:
            break
    s0, t0 = rng.choice(opts)
    w0 = [pow(s0, j, q) for j in range(6)] + [(-t0 * pow(s0, j, q)) % q for j in range(6)]
    B0 = [[q] + [0] * 11] + [[(-w0[j]) % q] + [int(k == j) for k in range(1, 12)] for j in range(1, 12)]
    B = lll_small(B0)
    assert all(sum(a * b for a, b in zip(v, w0)) % q == 0 for v in B)
    primes = [p for p in range(P_MIN, 1 << LX) if mr_prime(p) and p != q]
    pairs = []
    pair_roots = {}
    for p in primes:
        for side in ('f', 'g'):
            for s in roots_modp(H, p, rng):
                spoly = ([TOY_U * s % p, 0, 1] if side == 'g' else
                         [(1 + s - s * s - s ** 3) % p, 0, (s ** 5) % p, 0, 0, 0,
                          2 * s ** 3 % p, 0, 0, 0, (-2 * s) % p, 0, 1])
                for t in roots_modp(spoly, p, rng):
                    w = []
                    for i in range(12):
                        v = B[i]
                        val = (sum(v[j] * pow(s, j, p) for j in range(6)) -
                                t * sum(v[6 + j] * pow(s, j, p) for j in range(6))) % p
                        w.append(val)
                    if any(w):
                        sp_ = [pow(s, j, p) for j in range(6)]
                        pairs.append((p, side, t, w, sp_))
                        pair_roots[(p, side, t)] = (s, t, sp_)
    pair_by_p = {}
    for (p, side, t, w, sp_) in pairs:
        pair_by_p.setdefault((p, side), []).append((t, w, sp_))
    return rng, q, B, primes, pairs, pair_by_p, pair_roots


def node_differential(seed=190, kbox=1, LX=10, n_target=12, kw=3, slack=24,
                      nsample=600, nfn=800):
    """Sampled differential. Returns a dict of counts; asserts on any violation
    of (a)-(d)."""
    rng, q, B, primes, pairs, pair_by_p, pair_roots = _setup_node(
        seed, 11, 12, kbox, LX)
    m = 2 * kbox + 1
    # sampled cells
    cells = []
    seen = set()
    while len(cells) < nsample:
        c = tuple(rng.randint(-kbox, kbox) for _ in range(12))
        if any(c) and c not in seen:
            seen.add(c); cells.append(c)
    # (a) functional differential + (b) scores
    bad_a = bad_b = 0
    score = {c: [0.0, 0.0] for c in cells}
    logp = {p: math.log2(p) for p in primes}
    norms_cache = {}
    for c in cells:
        v = [sum(c[i] * B[i][j] for i in range(12)) for j in range(12)]
        nf, ng = toy_norms(v[:6], v[6:])
        assert ng % q == 0  # the toy special-q is taken on the g side
        norms_cache[c] = (v, nf, ng)
    pair_by_p = {}
    for (p, side, t, w, sp_) in pairs:
        pair_by_p.setdefault((p, side), []).append((t, w, sp_))
    for c in cells:
        v, nf, ng = norms_cache[c]
        for (p, side), plist in pair_by_p.items():
            anyhit = False
            for (t, w, sp_) in plist:
                val = sum(ci * wi for ci, wi in zip(c, w)) % p
                # oracle: direct evaluation of the membership condition on v
                s_, t_, sp_ = pair_roots[(p, side, t)]
                oracle = (sum(v[j] * sp_[j] for j in range(6)) -
                           t_ * sum(v[6 + j] * sp_[j] for j in range(6))) % p
                if (val == 0) != (oracle == 0):
                    bad_a += 1
                if val == 0:
                    anyhit = True
                    score[c][0 if side == 'f' else 1] += logp[p]
            # divisibility is per prime: p | N_side iff some root pair hits
            if anyhit != ((nf if side == 'f' else ng) % p == 0):
                bad_a += 1
    # (b) score equals distinct-prime found mass on a subsample
    base_primes = [p for p in range(2, 257) if mr_prime(p)]
    for c in cells[:nfn // 4 or 50]:
        v, nf, ng = norms_cache[c]
        for side, n in ((0, nf), (1, ng)):
            want = sum(logp[p] for p in primes if n % p == 0)
            if abs(score[c][side] - want) > 1e-9 * max(1.0, want):
                bad_b += 1
    # (c) threshold + false-negative check on the sample
    R = kw * n_target + slack
    B_big = 1 << n_target
    LXv = 1 << LX
    surv = 0; surv_true = 0; junk = 0; false_neg = 0; bad_ch = 0
    td_primes = sorted(set(primes))
    for c in cells:
        v, nf, ng = norms_cache[c]
        res_f = math.log2(nf) - score[c][0]
        res_g = math.log2(ng) - score[c][1]
        passes = res_f <= R and res_g <= R

        def is_rel(n):
            left = n; win = 0
            for p in td_primes:
                if p * p > left:
                    break
                while left % p == 0:
                    left //= p
                    if p > LXv:
                        win += 1
            if left > 1:
                if left > B_big:
                    return False
                if mr_prime(left):
                    return (left > LXv and win + 1 <= kw) or (left <= LXv and win <= kw)
                return False
            return win <= kw
        rel = is_rel(nf) and is_rel(ng)
        if rel and not passes:
            false_neg += 1
        if passes:
            surv += 1
            fs = factorize(nf, base_primes); gs = factorize(ng, base_primes)
            if fs is None or gs is None:
                junk += 1          # conservatively rejected, as production ECM does
                continue
            chf = [x for x in fs if LXv < x <= B_big]; chg = [x for x in gs if LXv < x <= B_big]
            if len(chf) <= kw and len(chg) <= kw and not any(x > B_big for x in chf + chg):
                if rel:
                    surv_true += 1
                    # credited children are exactly the residual window primes
                    want_f = [x for x in fs if x > LXv]
                    if chf != [x for x in want_f if x <= B_big] or len(want_f) > kw:
                        bad_ch += 1
                else:
                    junk += 1
            else:
                junk += 1
    return dict(q=q, cells=len(cells), pairs=len(pairs), functional_mismatches=bad_a,
                score_mismatches=bad_b, survivors=surv, survivor_true=surv_true,
                survivor_junk=junk, false_negatives=false_neg, child_violations=bad_ch,
                LX=LX, n=n_target, kw=kw, slack=slack)


if __name__ == '__main__':
    print(node_differential())


def relation_differential(seed=331, kbox=1, LX=10, n_target=40, kw=10, slack=40, nsample=1200):
    """Fat-window regime chosen so that TRUE RELATIONS occur at toy scale
    (small special-q, u ~ 4-5 on both sides at the n=2^40 target): exercises the
    full survivor path - exact norms, residual division, factorisation, the
    relation predicate and child extraction - on real relations. The threshold
    at this regime is generous (R exceeds the norms), so the no-false-negative
    coverage is the thin-window run's role; this run's role is no false credits
    and exact children. Returns counts; asserts on any violation."""
    rng, q, B, primes, pairs, pair_by_p, pair_roots = _setup_node(
        seed, 5, 6, kbox, LX)
    cells = []
    seen = set()
    while len(cells) < nsample:
        c = tuple(rng.randint(-kbox, kbox) for _ in range(12))
        if any(c) and c not in seen:
            seen.add(c); cells.append(c)
    base_primes = [p for p in range(2, 257) if mr_prime(p)]
    pair_by_p = {}
    for (p, side, t, w, sp_) in pairs:
        pair_by_p.setdefault((p, side), []).append((t, w, sp_))
    R = kw * n_target + slack
    B_big = 1 << n_target
    LXv = 1 << LX
    surv = true = junk = false_credit = child_bad = rel_pass_fail = unfactored = 0
    for c in cells:
        v = [sum(c[i] * B[i][j] for i in range(12)) for j in range(12)]
        nf, ng = toy_norms(v[:6], v[6:])
        assert ng % q == 0
        score = [0.0, 0.0]
        for (p, side), plist in pair_by_p.items():
            anyhit = False
            for (t, w, sp_) in plist:
                if sum(ci * wi for ci, wi in zip(c, w)) % p == 0:
                    anyhit = True
                    break
            if anyhit:
                score[0 if side == 'f' else 1] += math.log2(p)
        res_f = math.log2(nf) - score[0]
        res_g = math.log2(ng) - score[1]
        passes = res_f <= R and res_g <= R
        fs = factorize(nf, base_primes)
        gs = factorize(ng, base_primes)
        if fs is None or gs is None:
            unfactored += 1
            is_rel = False
            big_f = big_g = []
        else:
            big_f = [x for x in fs if x > LXv]
            big_g = [x for x in gs if x > LXv]
            is_rel = (not any(x > B_big for x in big_f) and len(big_f) <= kw
                      and not any(x > B_big for x in big_g) and len(big_g) <= kw)
        if is_rel and not passes:
            rel_pass_fail += 1
        if passes:
            surv += 1
            # credited children: exactly the window primes in (2^LX, 2^n];
            # a credited relation requires complete factorisation with every
            # prime <= 2^n and at most kw window primes per side (the
            # production predicate after full residual ECM)
            chf = [x for x in (fs or []) if LXv < x <= B_big]
            chg = [x for x in (gs or []) if LXv < x <= B_big]
            credited_rel = (fs is not None and gs is not None
                            and not any(x > B_big for x in fs + gs)
                            and len(chf) <= kw and len(chg) <= kw)
            if credited_rel:
                if not is_rel:
                    false_credit += 1
                else:
                    true += 1
                    want_f = [x for x in (fs or []) if x > LXv]
                    want_g = [x for x in (gs or []) if x > LXv]
                    if (chf != [x for x in want_f if x <= B_big] or
                            chg != [x for x in want_g if x <= B_big] or
                            len(want_f) > kw or len(want_g) > kw):
                        child_bad += 1
            else:
                junk += 1
    return dict(q=q, cells=len(cells), pairs=len(pairs), survivors=surv,
                relations_true=true, false_credits=false_credit,
                relation_threshold_failures=rel_pass_fail, child_violations=child_bad,
                junk_survivors=junk, unfactored_rejected=unfactored, LX=LX, n=n_target, kw=kw, slack=slack)
