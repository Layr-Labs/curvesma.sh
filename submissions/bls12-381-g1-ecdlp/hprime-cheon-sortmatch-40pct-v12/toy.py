#!/usr/bin/env python3
"""Toy evidence for the H-prime Cheon BSGS with sort-based matching (v10, stdlib only).

Part A  algebra of the production field r (d = 2^20) and the H-prime density on it.
Part B  H-prime density and gap maxima on three toy fields: the per-prime fraction
        against 1/d, with the per-integer density kept as a separate quantity.
Part C  ratio distinctness, exhaustively at toy scale: x = a*w with a, w H-primes
        determines (a, w), and the ratios x/y are pairwise distinct except the
        w = w' diagonal, whose share is measured against (|W|-1)/|W|^2.
Part D  the matching pipeline on a genuine j = 0 curve: gap-addition walk, truncated
        y-key, radix routing, sort, equal-key runs, capped candidate verification,
        d-th root and recovery of s.  Planted challenges force the recovery path;
        random challenges check that no wrong answer is ever emitted.
Part E  coverage: an exact distinct-coverage census at three toy sizes.  At each
        size the full candidate set {x/(u x')} is enumerated, its exact set
        cardinality D is counted, and D is compared with the random-set occupancy
        prediction m(1-exp(-K_cand/m)), K_cand = candidate occurrences / m, both on
        the raw occurrence count and on the count with the toy's shared-walk
        diagonal repetitions removed.  The hit rate on 20000 uniform challenges is
        reported as a separate check on the challenge sampler (D/m), not as a test
        of the distinct-set prediction.

The model is H1/H3occ: the candidate set behaves like a random multiset of
K_cand draws from H, so the expected distinct fraction is m(1-exp(-K_cand/m)).
The occupancy form already contains every random repetition of a multiset; only
systematic (exactly repeated) candidates sit outside it, and those are measured.
"""
import json
import math
import random
import sys
from collections import Counter

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


class Curve:
    """y^2 = x^3 + b over F_p, p = 3 mod 4, affine; None is the identity."""

    def __init__(self, p, b):
        self.p, self.b = p, b

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


def sqrt_mod(a, r):
    a %= r
    if a == 0:
        return 0
    if pow(a, (r - 1) // 2, r) != 1:
        return None
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
        m, c, t, x = i, b * b % r, t * b * b % r, x * b % r
    return x


BLS_P = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624'
            '1eabfffeb153ffffb9feffffffffaaab', 16)
BLS_R = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
BLS_FAC = [(2, 32), (3, 1), (11, 1), (19, 1), (10177, 1), (125527, 1), (859267, 1),
           (906349, 2), (2508409, 1), (2529403, 1), (52437899, 1), (254760293, 2)]
BLS_D = 1 << 20

TOY_P = 268520611
TOY_R = 67122049
TOY_H = 4
TOY_D = 1 << 3            # 8 | r-1 ; v2(r-1) = 7 ; 6 | m
TOY_M = (TOY_R - 1) // TOY_D


# ------------------------------------------------------------------ part A
def part_a():
    out = {}
    r, d = BLS_R, BLS_D
    prod = 1
    for q, e in BLS_FAC:
        prod *= q ** e
    out['r_prime'] = is_prime(r)
    out['r_minus_1_factorisation'] = prod == r - 1 and all(is_prime(q) for q, _ in BLS_FAC)
    out['d_divides_r_minus_1'] = (r - 1) % d == 0
    out['v2_is_exactly_32'] = (r - 1) % (1 << 32) == 0 and (r - 1) % (1 << 33) != 0
    m = (r - 1) // d
    out['m_bit_length'] = m.bit_length()
    out['six_divides_m_mu6_in_H'] = m % 6 == 0
    out['six_d_divides_r_minus_1'] = (r - 1) % (6 * d) == 0
    g = 2
    while any(pow(g, (r - 1) // q, r) == 1 for q, _ in BLS_FAC):
        g += 1
    zeta = pow(g, m, r)
    out['zeta_order_d'] = pow(zeta, d, r) == 1 and pow(zeta, d // 2, r) != 1
    rng = random.Random(20261023)
    xs = [pow(rng.randrange(2, r), d, r) for _ in range(6)]
    out['hprime_set_closed_under_product'] = all(pow(x * y % r, m, r) == 1 for x in xs for y in xs)
    out['dth_powers_in_H'] = all(pow(pow(x, d, r), m, r) == 1 for x in xs)
    out['two_pow_m_is_not_1'] = pow(2, m, r) != 1
    primes = [q for q in range(3, 1 << 16) if is_prime(q)]
    dens = {}
    for dd in (16, 256):
        e = (r - 1) // dd
        cnt = sum(1 for q in primes if pow(q, e, r) == 1)
        dens[str(dd)] = dict(count=cnt, primes=len(primes),
                             fraction=round(cnt / len(primes), 6),
                             one_over_dprime=round(1.0 / dd, 6))
    out['dprime_power_density_below_2_16'] = dens
    return out


# ------------------------------------------------------------------ part B
def toy_field(seed, dd, bits, m_min=5000):
    rng = random.Random(seed)
    while True:
        q = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if not is_prime(q) or (q - 1) % dd:
            continue
        m = (q - 1) // dd
        if m >= m_min and m % 6 == 0:
            return q, m


def part_b():
    rows = []
    for (seed, dd, bits, lo, span) in ((11, 16, 22, 1 << 18, 30000),
                                       (12, 64, 24, 1 << 20, 30000),
                                       (13, 256, 26, 1 << 22, 30000)):
        r, m = toy_field(seed, dd, bits)
        hp = tp = 0
        gaps, last = [], None
        for x in range(lo, lo + span):
            if not is_prime(x):
                continue
            tp += 1
            if pow(x, m, r) == 1:
                hp += 1
                if last is not None:
                    gaps.append(x - last)
                last = x
        rows.append(dict(d=dd, r_bits=bits, window=span, primes=tp, hprimes=hp,
                         per_prime=round(hp / tp, 6), one_over_d=round(1.0 / dd, 6),
                         max_gap=max(gaps) if gaps else 0,
                         mean_gap=round(sum(gaps) / len(gaps), 2) if gaps else 0,
                         gaps_measured=len(gaps)))
    return rows


# ------------------------------------------------------------------ ranges
def hprimes_range(r, m, lo, hi, need):
    out = []
    for x in range(lo, hi):
        if is_prime(x) and pow(x, m, r) == 1:
            out.append(x)
            if len(out) >= need:
                break
    return out


def lists(r, m, nw, nb, ng):
    walk = [w for w in hprimes_range(r, m, 1 << 8, 1 << 12, nw) if True]
    base_b = hprimes_range(r, m, 1 << 12, 1 << 14, nb + 4)
    base_g = hprimes_range(r, m, (1 << 12) + 4096, 1 << 14, ng + 4)
    base_b = [a for a in base_b if all(a * w < r for w in walk)][:nb]
    base_g = [a for a in base_g if all(a * w < r for w in walk)][:ng]
    return walk, base_b, base_g


def lists_c(r, m, nw, nb, ng):
    walk = hprimes_range(r, m, 1 << 8, 1 << 15, nw)
    bb = hprimes_range(r, m, 1 << 15, 1 << 19, nb)
    bg = hprimes_range(r, m, (1 << 19) + 200000, 1 << 20, ng)
    return walk, bb, bg


# ------------------------------------------------------------------ part C
def part_c_field():
    """r about 2^40 with 64 | r-1 so that every product a*w and every cross product
    a*A of the toy lists is an exact integer below r (no modular wrap)."""
    for bits in (40, 41, 42):
        r, m = toy_field(2027, 64, bits, m_min=1 << 20)
        if r > (1 << 19) * (1 << 19) and (r - 1) % 64 == 0:
            return r, m
    raise AssertionError('no part-C field')


def part_c(r, m, nw=37, nb=8, ng=8):
    walk, bb, bg = lists_c(r, m, nw, nb, ng)
    pb = [a * w for a in bb for w in walk]
    pg = [a * w for a in bg for w in walk]
    out = dict(walk_primes=len(walk), base_baby=len(bb), base_giant=len(bg),
               baby_products=len(pb), giant_products=len(pg),
               products_distinct=len(set(pb)) == len(pb) and len(set(pg)) == len(pg),
               walk_below_base=walk[-1] < bb[0] and walk[-1] < bg[0],
               max_product_below_r=max(pb + pg) < r,
               disjoint_base_pools=not (set(bb) & set(bg)))
    inv_w = {w: pow(w, -1, r) for w in walk}
    inv_g = {a: pow(a, -1, r) for a in bg}
    cnt = Counter()
    for a in bb:
        for A in bg:
            ab = a * inv_g[A] % r
            for w in walk:
                for w2 in walk:
                    cnt[ab * (w * inv_w[w2] % r) % r] += 1
    total = len(pb) * len(pg)
    mm1 = sum(v - 1 for v in cnt.values())
    out.update(pairs=total, distinct_ratios=len(cnt), max_multiplicity=max(cnt.values()),
               sum_multiplicity_minus_one=mm1,
               predicted_multiplicity_minus_one=len(bb) * len(bg) * (len(walk) - 1),
               loss_fraction_measured=round(mm1 / total, 9),
               loss_fraction_predicted=round((len(walk) - 1) / len(walk) ** 2, 9),
               collision_classes=sum(1 for v in cnt.values() if v > 1))
    return out


# ------------------------------------------------------------------ part D
def toy_curve():
    order = TOY_H * TOY_R
    for b in range(1, 60):
        C = Curve(TOY_P, b)
        pts, ok = [], True
        for x in range(1, 90):
            Q = C.point(x)
            if Q is None:
                continue
            if C.mul(order, Q) is not None:
                ok = False
                break
            pts.append(Q)
        if not ok or len(pts) < 4:
            continue
        for Q in pts:
            if C.mul(TOY_R, Q) is None:
                return C, Q
            A = C.mul(TOY_H // 2, Q)
            if A is not None and C.mul(TOY_R, A) is None:
                return C, A
    raise AssertionError('no twist found')


def toy_env(seed=7):
    C, A = toy_curve()
    r, d, m = TOY_R, TOY_D, TOY_M
    g = 2
    while pow(g, (r - 1) // 2, r) == 1 or pow(g, (r - 1) // 3, r) == 1:
        g += 1
    zeta = pow(g, m, r)
    lam = pow(g, (r - 1) // 3, r)
    omega = pow(g % TOY_P, (TOY_P - 1) // 3, TOY_P)
    if C.mul(lam, A) != (omega * A[0] % TOY_P, A[1]):
        lam = lam * lam % r
    assert C.mul(lam, A) == (omega * A[0] % TOY_P, A[1])
    assert pow(zeta, d, r) == 1 and pow(zeta, d // 2, r) != 1
    return dict(C=C, p=TOY_P, r=r, d=d, m=m, zeta=zeta, lam=lam, omega=omega, P0=A,
                rng=random.Random(seed))


def mu6_scalars(env):
    r, lam = env['r'], env['lam']
    l2 = lam * lam % r
    return [1, r - 1, lam, r - lam, l2, r - l2]


def key_of(P, p, bits, route_bits):
    y = P[1]
    if y > p - y:
        y = p - y
    k = y & ((1 << bits) - 1)
    return k >> (bits - route_bits), k & ((1 << (bits - route_bits)) - 1)


def walk_chain(env, base, elems, z=None):
    """[base*elems[j](*z)]P0 by gap additions (the production walk)."""
    C, P0, r = env['C'], env['P0'], env['r']
    s = base * elems[0] % r
    if z is not None:
        s = s * z % r
    cur = C.mul(s, P0)
    pts = [cur]
    tab = {}
    for j in range(len(elems) - 1):
        gam = elems[j + 1] - elems[j]
        key = (base * gam % r) * (z if z is not None else 1) % r
        if key not in tab:
            tab[key] = C.mul(key, P0)
        cur = C.add(cur, tab[key], incomplete=True)
        pts.append(cur)
    return pts


def dth_root(env, z):
    """A d-th root of z in F_r^* (the root chain used by the reduction)."""
    r, d = env['r'], env['d']
    s = z
    for _ in range(d.bit_length() - 1):
        s = sqrt_mod(s, r)
        if s is None:
            return None
    return s


def part_d(env, trials=40, key_bits=18, route_bits=6, nb=8, ng=8, nw=37):
    C, p, r, m, d = env['C'], env['p'], env['r'], env['m'], env['d']
    walk, bb, bg = lists(r, m, nw, nb, ng)
    baby_scal, baby_pts = [], []
    for a in bb:
        baby_pts += walk_chain(env, a, walk)
        baby_scal += [a * w % r for w in walk]
    us = mu6_scalars(env)
    n_rec = len(baby_pts)
    cap_pairs = max(64, 4 * (n_rec * n_rec) // (1 << key_bits) + 8)
    cap_ver = 6 * cap_pairs
    st = dict(records=2 * n_rec, key_bits=key_bits, route_bits=route_bits, trials=0,
              planted=0, recovered_z=0, recovered_s=0, false_positive=0,
              max_run_pairs=0, total_run_pairs=0, max_verifications=0,
              cap_pairs=cap_pairs, cap_verifications=cap_ver, cap_exceeded=0,
              incomplete_trips=0, nodes=1 << route_bits, coincident_wrong=0,
              s_checked_against_T=0, unresolved_challenges=0, diagonal_trials=0, expected_false_pairs=round(n_rec * n_rec / 2 / 2 ** key_bits, 3))
    rng = env['rng']
    plans = []
    for _ in range(trials):
        x0 = rng.choice(baby_scal)
        y0 = rng.choice(bg) * rng.choice(walk) % r
        z = x0 * pow(y0, -1, r) % r
        sw = dth_root(env, z)
        plans.append(('planted', z, None if sw is None else C.mul(sw, env['P0'])))
    for _ in range(trials // 4):
        sw = rng.randrange(2, r)
        plans.append(('random', pow(sw, d, r), C.mul(sw, env['P0'])))
    for kind, z, T in plans:
        Ad = C.mul(z, env['P0'])
        st['unresolved_challenges'] += 1 if T is None else 0
        giant_pts, giant_scal = [], []
        for A in bg:
            for w in walk:
                y = A * w % r
                giant_scal.append(y)
                giant_pts.append(C.mul(z * y % r, env['P0']))
        recs = []
        for i, P in enumerate(baby_pts):
            k0, k1 = key_of(P, p, key_bits, route_bits)
            recs.append((k0, k1, 0, baby_scal[i]))
        for i, P in enumerate(giant_pts):
            k0, k1 = key_of(P, p, key_bits, route_bits)
            recs.append((k0, k1, 1, giant_scal[i]))
        recs.sort(key=lambda t: (t[0], t[1]))
        runs, run = [], [recs[0]]
        for rec in recs[1:]:
            if rec[:2] == run[-1][:2]:
                run.append(rec)
            else:
                runs.append(run)
                run = [rec]
        runs.append(run)
        pair_cmp = vers = 0
        found, extra = None, 0
        for grp in runs:
            for i in range(len(grp)):
                for j in range(i + 1, len(grp)):
                    if grp[i][2] == grp[j][2]:
                        continue
                    pair_cmp += 1
                    xb = grp[i][3] if grp[i][2] == 0 else grp[j][3]
                    xg = grp[j][3] if grp[j][2] == 1 else grp[i][3]
                    for u in us:
                        vers += 1
                        cand = xb * pow(u * xg % r, -1, r) % r
                        if C.mul(cand, env['P0']) == Ad:
                            if found is None:
                                found = cand
                            elif cand != found:
                                extra += 1
        st['trials'] += 1
        st['max_run_pairs'] = max(st['max_run_pairs'], pair_cmp)
        st['total_run_pairs'] += pair_cmp
        st['max_verifications'] = max(st['max_verifications'], vers)
        if pair_cmp > cap_pairs or vers > cap_ver:
            st['cap_exceeded'] += 1
        if pair_cmp == len(walk):
            st['diagonal_trials'] += 1
        if kind == 'planted':
            st['planted'] += 1
        if found is not None and found == z:
            st['recovered_z'] += 1
            s0 = found
            for _ in range(d.bit_length() - 1):
                s0 = sqrt_mod(s0, r) if s0 is not None else None
                if s0 is None:
                    break
            ok = False
            if s0 is not None:
                for j in range(d):
                    cand = s0 * pow(env['zeta'], j, r) % r
                    if T is not None and C.mul(cand, env['P0']) == T:
                        ok = True
                        break
            st['recovered_s'] += 1 if ok else 0
            st['s_checked_against_T'] += 1 if ok else 0
        elif found is not None:
            st['false_positive'] += 1
        st['coincident_wrong'] += extra
    return st


# ------------------------------------------------------------------ part E
def coverage_census(env, nb, ng, nw, challenges=20000):
    """Exact distinct-coverage census of the candidate set at one toy size.

    Enumerates every candidate v = x/(u x') (x from the baby side, x' from the
    giant side, u in mu6), counts the exact set cardinality D, and compares it
    with the occupancy prediction for a random multiset of the same occurrence
    count, K_cand = occurrences / m.  The toy's shared-walk diagonal is counted
    separately because it is the one systematic repetition the toy has.
    """
    r, m, d = env['r'], env['m'], env['d']
    walk, bb, bg = lists(r, m, nw, nb, ng)
    inv_w = {w: pow(w, -1, r) for w in walk}
    inv_g = {a: pow(a, -1, r) for a in bg}
    us = mu6_scalars(env)
    invu = [pow(u, -1, r) for u in us]
    cov = set()
    for a in bb:
        for A in bg:
            ab = a * inv_g[A] % r
            for w in walk:
                for w2 in walk:
                    base = ab * (w * inv_w[w2] % r) % r
                    for iu in invu:
                        cov.add(base * iu % r)
    nbs, ngs = len(bb) * len(walk), len(bg) * len(walk)
    occ = 6 * nbs * ngs                       # candidate occurrences
    D = len(cov)                              # exact distinct covered residues
    K = occ / m
    pred = int(round(m * (1 - math.exp(-K))))  # random-set occupancy prediction
    diag_reps = 6 * len(bb) * len(bg) * (len(walk) - 1)   # shared-walk diagonal
    Ke = (occ - diag_reps) / m
    prede = int(round(m * (1 - math.exp(-Ke))))
    rng = random.Random(99)
    hits = 0
    for _ in range(challenges):
        if pow(rng.randrange(2, r), d, r) in cov:
            hits += 1
    return dict(nb=len(bb), ng=len(bg), nw=len(walk), records_per_side=nbs,
                pairs=nbs * ngs, candidate_occurrences=occ, distinct_covered=D, m=m,
                K_cand=round(K, 6), predicted_distinct=pred,
                ratio_measured_over_predicted=round(D / pred, 6),
                shortfall_over_m=round((pred - D) / m, 8),
                diagonal_repetitions=diag_reps, K_eff=round(Ke, 6),
                predicted_distinct_eff=prede,
                ratio_eff=round(D / prede, 6),
                covered_over_m=round(D / m, 6),
                naive_six_nn_over_m=round(6 * nbs * ngs / m, 6),
                measured_hits=hits, challenges=challenges,
                measured_success=round(hits / challenges, 6),
                z_score=round((hits - challenges * D / m) /
                              math.sqrt(challenges * (D / m) * (1 - D / m)), 3))


CENSUS_SIZES = ((12, 12, 37), (12, 12, 57), (12, 12, 73))


def part_e(env, challenges=20000, nb=12, ng=12, nw=37):
    c = coverage_census(env, nb, ng, nw, challenges)
    return dict(records_per_side=c['records_per_side'], pairs=c['pairs'],
                distinct_covered=c['distinct_covered'], m=c['m'],
                candidate_occurrences=c['candidate_occurrences'],
                covered_over_m=c['covered_over_m'],
                naive_six_nn_over_m=c['naive_six_nn_over_m'],
                diagonal_loss_predicted=round((c['nw'] - 1) / c['nw'] ** 2, 6),
                predicted_success=c['covered_over_m'],
                measured_success=c['measured_success'],
                measured_hits=c['measured_hits'], challenges=c['challenges'],
                z_score=c['z_score'])


def main():
    env = toy_env()
    rc, mc = part_c_field()
    out = dict(partA=part_a(), partB=part_b(),
               partC=dict(part_c(rc, mc), field_r_bits=rc.bit_length(), field_d=64),
               partD=part_d(env), partE=part_e(env),
               partEcensus=[coverage_census(env, *s) for s in CENSUS_SIZES],
               toy=dict(p=TOY_P, r=TOY_R, d=TOY_D, m=TOY_M,
                        curve='y^2=x^3+b (j=0), #E=4r, order-r point, mu6 acts as (omega*x, y)'))
    json.dump(out, open('toy-summary.json', 'w'), indent=1, sort_keys=True)
    print(json.dumps(out, indent=1, sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
