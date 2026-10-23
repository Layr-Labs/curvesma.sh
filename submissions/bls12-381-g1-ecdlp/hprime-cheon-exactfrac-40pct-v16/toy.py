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


def key_of(P, p, bits, route_bits):
    """Even representative of y, truncated to `bits`, as 4 limbs of `bits//4` bits.

    Production keeps the 224-bit key as four 56-bit limbs and a 48-byte record;
    the toy keeps a (bits)-bit key as four limbs and a 6-word (48-byte) record,
    so every copy and compare below has the same word count as production.
    """
    y = P[1]
    if y > p - y:
        y = p - y
    key = y & ((1 << bits) - 1)
    limb_bits = bits // 4
    limbs = [(key >> (limb_bits * i)) & ((1 << limb_bits) - 1) for i in range(4)]
    node = key >> (bits - route_bits)              # routed on the top route_bits
    local = key & ((1 << (bits - route_bits)) - 1)  # ordered locally on the rest
    return node, local, limbs, key


# production radix parameters, mirrored by the toy: 3 LSD passes over the
# locally-ordered bits, the routed high bits excluded, and the middle bits that
# no pass orders compared explicitly in the group scan.
SORT_PASSES = 3
KEY_BITS = 24
ROUTE_BITS = 6
SORT_BITS = 12
SORT_DIGIT = SORT_BITS // SORT_PASSES     # 4 bits: toy digit width, 3 passes
# per-record-per-pass instruction weights, identical to the ledger itemisation
PASS_OPS = dict(digit=1, counter=3, address=2, copy=12, loop=2)   # = 21 (UBFX per digit)
# per-pair instruction weights of the full-key group scan, identical to the ledger
SCAN_OPS = dict(key_load=4, compare=3, branch=1, emit=1, loop=1)  # = 10


def radix_pass(recs, shift, width, ops):
    """One LSD counting-sort pass over `width` bits at `shift` (production shape).

    recs are 6-word records [l0,l1,l2,l3,side,scalar]; the counted operations and
    their weights are exactly the ledger's sort row.
    """
    mask = (1 << width) - 1
    cnt = [0] * (1 << width)
    for rec in recs:
        digit = ((rec[0] | (rec[1] << 6)) >> shift) & mask
        ops['digit'] += PASS_OPS['digit']
        cnt[digit] += 1
        ops['counter'] += PASS_OPS['counter']
    run = 0
    for i in range(1 << width):            # prefix sums: per bucket, per-node work
        run += cnt[i]
        cnt[i] = run
        ops['prefix'] += 1
    out = [None] * len(recs)
    for i in range(len(recs) - 1, -1, -1):
        rec = recs[i]
        digit = ((rec[0] | (rec[1] << 6)) >> shift) & mask
        ops['digit'] += PASS_OPS['digit']
        cnt[digit] -= 1
        ops['address'] += PASS_OPS['address']
        out[cnt[digit]] = list(rec)        # 48-byte record: 6 words
        ops['copy'] += PASS_OPS['copy']
        ops['loop'] += PASS_OPS['loop']
    return out


def radix_sort(recs, ops):
    for i in range(SORT_PASSES):
        recs = radix_pass(recs, i * SORT_DIGIT, SORT_DIGIT, ops)
    return recs


def group_runs(recs, ops):
    """Split into groups of equal locally-ordered bits (the sort's key)."""
    mask = (1 << SORT_BITS) - 1
    runs, run = [], [recs[0]]
    for rec in recs[1:]:
        if (rec[0] | (rec[1] << 6)) & mask == ((run[-1][0] | (run[-1][1] << 6)) & mask):
            run.append(rec)
        else:
            runs.append(run)
            run = [rec]
    runs.append(run)
    return runs


def full_key_equal(a, b, ops):
    """The group scan's full-key compare: all 4 limbs, not just the sorted bits."""
    ops['scan_key_load'] += SCAN_OPS['key_load'] // 2
    ops['scan_key_load'] += SCAN_OPS['key_load'] // 2
    ops['scan_compare'] += SCAN_OPS['compare']
    ops['scan_branch'] += SCAN_OPS['branch']
    same = a[0] == b[0] and a[1] == b[1] and a[2] == b[2] and a[3] == b[3]
    ops['scan_loop'] += SCAN_OPS['loop']
    ops['scan_emit'] += SCAN_OPS['emit']        # emit slot of the run scan
    return same


def part_d(env, trials=200, nb=8, ng=8, nw=37):
    C, p, r, m, d = env['C'], env['p'], env['r'], env['m'], env['d']
    walk, bb, bg = lists(r, m, nw, nb, ng)
    baby_scal, baby_pts = [], []
    for a in bb:
        baby_pts += walk_chain(env, a, walk)
        baby_scal += [a * w % r for w in walk]
    us = mu6_scalars(env)
    n_rec = len(baby_pts)
    cap_pairs = max(64, 4 * (n_rec * n_rec) // (1 << KEY_BITS) + 8)
    cap_ver = 6 * cap_pairs
    st = dict(records=2 * n_rec, key_bits=KEY_BITS, route_bits=ROUTE_BITS,
              sort_passes=SORT_PASSES, sort_digit=SORT_DIGIT, sort_bits=SORT_BITS,
              middle_bits_scan_only=KEY_BITS - ROUTE_BITS - SORT_BITS,
              trials=0, planted=0, recovered_z=0, recovered_s=0, false_positive=0,
              max_run_pairs=0, total_run_pairs=0, max_verifications=0,
              cap_pairs=cap_pairs, cap_verifications=cap_ver, cap_exceeded=0,
              incomplete_trips=0, nodes=1 << ROUTE_BITS, coincident_wrong=0,
              s_checked_against_T=0, unresolved_challenges=0, diagonal_trials=0,
              z_from_secret=0, s_pow_d_matches_z=0,
              expected_false_pairs=round(n_rec * n_rec / 2 / 2 ** KEY_BITS, 3))
    rng = env['rng']
    plans = []
    # M5 end-to-end: ONE secret s per trial.  The challenge is built from that
    # secret (z = s^d, A_d = [z]P0, T = [s]P0); the pipeline recovers z from the
    # records alone, then s from z, and the *original* s is checked only at the end.
    for _ in range(trials):
        s_sec = rng.randrange(2, r)
        z = pow(s_sec, d, r)
        plans.append(('endtoend', z, C.mul(s_sec, env['P0']), s_sec))
    ops = dict((k, 0) for k in ('digit', 'counter', 'prefix', 'address', 'copy', 'loop',
                                'scan_key_load', 'scan_compare', 'scan_branch',
                                'scan_emit', 'scan_loop'))
    for kind, z, T, s_sec in plans:
        Ad = C.mul(z, env['P0'])
        st['z_from_secret'] += 1 if (T is not None and pow(s_sec, d, r) == z) else 0
        giant_pts, giant_scal = [], []
        for A in bg:
            for w in walk:
                y = A * w % r
                giant_scal.append(y)
                giant_pts.append(C.mul(z * y % r, env['P0']))
        # records routed to their owner node on the top route_bits, exactly as
        # production routes on the top 85 bits; each node sorts its own block
        nodes = {}
        for i, P in enumerate(baby_pts):
            node, local, limbs, key = key_of(P, p, KEY_BITS, ROUTE_BITS)
            nodes.setdefault(node, []).append(limbs + [0, baby_scal[i]])
        for i, P in enumerate(giant_pts):
            node, local, limbs, key = key_of(P, p, KEY_BITS, ROUTE_BITS)
            nodes.setdefault(node, []).append(limbs + [1, giant_scal[i]])
        pair_cmp = vers = 0
        found, extra = None, 0
        for node in sorted(nodes):
            recs = nodes[node]
            if len(recs) < 1:
                continue
            recs = radix_sort(recs, ops)
            for grp in group_runs(recs, ops):
                for i in range(len(grp)):
                    for j in range(i + 1, len(grp)):
                        if grp[i][4] == grp[j][4]:
                            continue
                        pair_cmp += 1
                        same = full_key_equal(grp[i], grp[j], ops)
                        if not same:
                            continue
                        xb = grp[i][5] if grp[i][4] == 0 else grp[j][5]
                        xg = grp[j][5] if grp[j][4] == 1 else grp[i][5]
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
        if kind == 'endtoend':
            st['planted'] += 1
        if found is not None and found == z:
            st['recovered_z'] += 1
            s0 = found
            for _ in range(d.bit_length() - 1):
                s0 = sqrt_mod(s0, r) if s0 is not None else None
                if s0 is None:
                    break
            ok = False
            ok_pow = False
            if s0 is not None:
                for j in range(d):
                    cand = s0 * pow(env['zeta'], j, r) % r
                    if T is not None and C.mul(cand, env['P0']) == T:
                        ok = True
                        if pow(cand, d, r) == z:
                            ok_pow = True
                        break
            st['recovered_s'] += 1 if ok else 0
            st['s_checked_against_T'] += 1 if ok else 0
            st['s_pow_d_matches_z'] += 1 if ok_pow else 0
        elif found is not None:
            st['false_positive'] += 1
        st['coincident_wrong'] += extra
    passes = st['trials'] * st['records'] * SORT_PASSES
    st['ops_per_record_per_pass'] = (
        ops['digit'] // passes + ops['counter'] // passes + ops['address'] // passes
        + ops['copy'] // passes + ops['loop'] // passes)
    # the digit is extracted once in the count pass and once in the scatter pass
    st['ledger_ops_per_record_per_pass'] = (2 * PASS_OPS['digit'] + PASS_OPS['counter']
                                            + PASS_OPS['address'] + PASS_OPS['copy']
                                            + PASS_OPS['loop'])
    pairs = st['total_run_pairs']
    st['ops_per_run_pair'] = (ops['scan_key_load'] // pairs + ops['scan_compare'] // pairs
                              + ops['scan_branch'] // pairs + ops['scan_emit'] // pairs
                              + ops['scan_loop'] // pairs)
    st['ledger_ops_per_run_pair'] = (SCAN_OPS['key_load'] + SCAN_OPS['compare']
                                     + SCAN_OPS['branch'] + SCAN_OPS['emit'] + SCAN_OPS['loop'])
    st['counted_ops'] = dict((k, v) for k, v in sorted(ops.items()))
    st['prefix_ops_per_node_pass'] = ops['prefix'] // max(1, st['trials'] * SORT_PASSES * st['nodes'])
    # the sampler is an end-to-end test of the coverage model too: z = s^d is a
    # UNIFORM element of H, so the recovery fraction must match the exactly
    # enumerated covered fraction D/m of the very same record set.
    cc = coverage_census(env, nb, ng, nw, 1)
    st['coverage_D_over_m'] = cc['covered_over_m']
    exp_hits = st['trials'] * cc['covered_over_m']
    sd = math.sqrt(st['trials'] * cc['covered_over_m'] * (1 - cc['covered_over_m']))
    st['expected_recoveries'] = round(exp_hits, 3)
    st['recovery_z'] = round((st['recovered_z'] - exp_hits) / sd, 3) if sd else None
    st['recovery_matches_coverage_model'] = bool(
        abs(st['recovery_z']) <= 3 and st['recovered_s'] == st['recovered_z'])
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


# ------------------------------------------------------------------ part F
def prime_field(bits):
    """A prime r about 2^bits with 384 = 64*6 dividing r-1, so every grid d
    (8, 16, 64) divides r-1 and every m = (r-1)/d is a multiple of 6."""
    x = (1 << bits) | 1
    while True:
        if (x - 1) % 384 == 0 and is_prime(x):
            return x
        x += 2


def grid_lists(r, m, nw, nb):
    walk = hprimes_range(r, m, 1 << 8, 1 << 22, nw)
    rest = hprimes_range(r, m, max(1 << 10, walk[-1] + 2), 1 << 22, 2 * nb)
    return walk, rest[:nb], rest[nb:2 * nb]


def grid_census(r, dd, target_K, nb=4):
    """One cell of the (field x d) coverage grid at a target K_cand ~ production K.

    Same candidate enumeration as part E, on a field given by its prime r and a
    divisor d of r-1 (no curve needed: the census is a statement about F_r^*).
    """
    m = (r - 1) // dd
    nw = max(24, int(round(math.sqrt(target_K * m / (6 * nb * nb)))))
    walk, bb, bg = grid_lists(r, m, nw, nb)
    g = 2
    while pow(g, (r - 1) // 2, r) == 1 or pow(g, (r - 1) // 3, r) == 1:
        g += 1
    lam6 = pow(g, (r - 1) // 6, r)                 # generator of the order-6 subgroup
    us = [pow(lam6, j, r) for j in range(6)]       # = {1,-1,omega,-omega,omega^2,-omega^2}
    invu = [pow(u, -1, r) for u in us]
    inv_w = {w: pow(w, -1, r) for w in walk}
    inv_g = {a: pow(a, -1, r) for a in bg}
    cov = set()
    for a in bb:
        for A in bg:
            ab = a * inv_g[A] % r
            for w in walk:
                base0 = ab * w % r
                for w2 in walk:
                    base = base0 * inv_w[w2] % r
                    for iu in invu:
                        cov.add(base * iu % r)
    nbs, ngs = len(bb) * len(walk), len(bg) * len(walk)
    occ = 6 * nbs * ngs
    D = len(cov)
    K = occ / m
    pred = int(round(m * (1 - math.exp(-K))))
    diag = 6 * len(bb) * len(bg) * (len(walk) - 1)
    Ke = (occ - diag) / m
    prede = int(round(m * (1 - math.exp(-Ke))))
    return dict(r=r, r_bits=r.bit_length(), d=dd, m=m, nw=len(walk), nb=len(bb), ng=len(bg),
                records_per_side=nbs, pairs=nbs * ngs, candidate_occurrences=occ,
                distinct_covered=D, K_cand=round(K, 6), predicted_distinct=pred,
                ratio_measured_over_predicted=round(D / pred, 6),
                shortfall_over_m=round((pred - D) / m, 8),
                diagonal_repetitions=diag, K_eff=round(Ke, 6),
                predicted_distinct_eff=prede, ratio_eff=round(D / prede, 6),
                shortfall_eff_over_m=round((prede - D) / m, 8),
                diagonal_share=round(K * (len(walk) - 1) / len(walk) ** 2, 8))


GRID_D = (8, 16, 64)
GRID_K = (0.51,)   # production K, one point per (field, d) cell


def part_f():
    fields = [prime_field(23), prime_field(25)]
    cells = [grid_census(r, dd, K) for r in fields for dd in GRID_D for K in GRID_K]
    eff = [c['shortfall_eff_over_m'] for c in cells]
    raw = [c['shortfall_over_m'] for c in cells]
    # first axis: does the residual fall with field size?  compare the |residual|
    # of the two fields cell by cell at matched (d, K)
    shrink = []
    for i, dd in enumerate(GRID_D):
        for j, K in enumerate(GRID_K):
            a = cells[i * len(GRID_K) + j]
            b = cells[len(GRID_D) * len(GRID_K) + i * len(GRID_K) + j]
            shrink.append(dict(d=dd, K=K, small_field_abs=abs(a['shortfall_eff_over_m']),
                               large_field_abs=abs(b['shortfall_eff_over_m']),
                               shrinks=abs(b['shortfall_eff_over_m']) < abs(a['shortfall_eff_over_m'])))
    return dict(fields=fields, field_bits=[f.bit_length() for f in fields],
                d_values=list(GRID_D), K_values=list(GRID_K), cells=cells,
                cells_count=len(cells),
                max_shortfall_eff=max(eff), min_shortfall_eff=min(eff),
                max_abs_shortfall_eff=max(abs(v) for v in eff),
                max_raw_shortfall=max(raw),
                raw_excess_per_cell=[dict(r_bits=c['r_bits'], d=c['d'],
                                          excess=c['shortfall_over_m'] - c['diagonal_share'])
                                     for c in cells],
                max_raw_excess=max(c['shortfall_over_m'] - c['diagonal_share'] for c in cells),
                raw_minus_diagonal=max(abs(c['shortfall_over_m'] - c['diagonal_share'])
                                       for c in cells),
                shrink_pairs=shrink, shrink_count=sum(1 for s in shrink if s['shrinks']))


# ------------------------------------------------------------------ part G
# 128-bit field.  Two premises the production claim needs are measured here
# rather than assumed:
#   (M2) the conditional uniformity of FULL keys given equal routed prefixes --
#        H3 as originally stated covers only the prefix, while the verification
#        cap is a statement about full-key collisions;
#   (M4) the H-prime gap tail, measured on a real consecutive H-prime list
#        instead of being read off a geometric model.
G_BITS = 128
G_D = 16


def g_field(bits=G_BITS):
    """A prime p near 2^bits with 48 | p-1, so d=16 divides p-1 and 6 | m."""
    x = (1 << (bits - 1)) | 1
    while True:
        if (x - 1) % 48 == 0 and is_prime(x):
            return x
        x += 2


def g_hprimes(p, m, lo, hi, need):
    """Consecutive integers in [lo,hi) that are prime and lie in H (q^m = 1 mod p)."""
    out, q = [], lo | 1
    while len(out) < need and q < hi:
        if is_prime(q) and pow(q, m, p) == 1:
            out.append(q)
        q += 2
    return out


def g_point(P_mod, b=7):
    C = Curve(P_mod, b)
    for x in range(1, 400):
        Q = C.point(x)
        if Q is not None:
            return C, Q
    raise AssertionError('no point')


def g_curve_field(bits=G_BITS):
    """The curve field (p = 3 mod 4), distinct from the scalar field as at production."""
    x = (1 << (bits - 1)) | 3
    while True:
        if x % 4 == 3 and is_prime(x):
            return x
        x += 4


def g_key(P, p):
    """Even representative of y: the production key rule, at full field width."""
    y = P[1]
    return y if y <= p - y else p - y


def part_g(pool=25, nw_rec=240, nw_gap=900, resid_bits=(8, 12, 16), routes=(6, 12)):
    p = g_field()
    d = G_D
    m = (p - 1) // d
    assert (p - 1) % d == 0 and m % 6 == 0
    hp = g_hprimes(p, m, 1 << 40, 1 << 62, 2 * pool + nw_rec + nw_gap)
    walk_gap = hp[:nw_gap]                       # consecutive H-primes: the gap tail
    rest = hp[nw_gap:]
    A, B = rest[:pool], rest[pool:2 * pool]
    W = rest[2 * pool:2 * pool + nw_rec]
    in_H = all(pow(q, m, p) == 1 for q in hp)
    # ---- (M4) measured gap tail, consecutive H-primes
    gaps = [walk_gap[i + 1] - walk_gap[i] for i in range(len(walk_gap) - 1)]
    mean_gap = sum(gaps) / len(gaps)
    max_gap = max(gaps)
    hist = [dict(bucket='2^%d' % j, count=sum(1 for g in gaps if (1 << j) <= g < (1 << (j + 1))))
            for j in range(0, 30) if any((1 << j) <= g < (1 << (j + 1)) for g in gaps)]
    tail = []
    for j in (9, 12, 15, 18, 21, 24, 27):
        T = 1 << j
        obs = sum(1 for g in gaps if g > T)
        tail.append(dict(threshold='2^%d' % j, observed=obs, n=len(gaps),
                         empirical=round(obs / len(gaps), 8),
                         geometric_model=round(math.exp(-T / mean_gap), 8)))
    max_even = max_gap + (max_gap & 1)
    table_bits = max(1, (max_even // 2).bit_length())
    model_under = [t for t in tail if t['empirical'] > t['geometric_model']]
    # ---- (M2) full-key conditional collision measurement at 128-bit keys
    pc = g_curve_field()
    C, P0 = g_point(pc)
    keys = []
    for a in A:
        for w in W:
            keys.append(g_key(C.mul(a * w % p, P0), pc))
    key_bits = max(keys).bit_length()
    out_routes = []
    for rb in routes:
        buckets = {}
        for k in keys:
            buckets.setdefault(k >> (key_bits - rb), []).append(k)
        n_pairs = sum(len(v) * (len(v) - 1) // 2 for v in buckets.values())
        full = 0
        for v in buckets.values():
            cnt = {}
            for k in v:
                cnt[k] = cnt.get(k, 0) + 1
            full += sum(c * (c - 1) // 2 for c in cnt.values())
        # positive control: the counter must fire when a duplicate is planted
        ctrl = 0
        for v in sorted(buckets.values(), key=len, reverse=True)[:4]:
            cnt = {}
            for k in v + [v[0]]:
                cnt[k] = cnt.get(k, 0) + 1
            ctrl += sum(c * (c - 1) // 2 for c in cnt.values())
        resid = []
        for rb2 in resid_bits:
            mask = (1 << rb2) - 1
            obs = 0
            worst = 0.0
            for v in buckets.values():
                cnt = {}
                for k in v:
                    r = k & mask
                    cnt[r] = cnt.get(r, 0) + 1
                obs += sum(c * (c - 1) // 2 for c in cnt.values())
                exp_b = len(v) * (len(v) - 1) / 2 / (1 << rb2)
                if exp_b >= 2.0:
                    obs_b = sum(c * (c - 1) // 2 for c in cnt.values())
                    worst = max(worst, abs(obs_b - exp_b) / exp_b)
            pred = n_pairs / (1 << rb2)
            sd = math.sqrt(pred) if pred > 0 else 0.0
            resid.append(dict(
                residual_bits=rb2, pairs=n_pairs, observed=obs,
                predicted=round(pred, 3),
                ratio=round(obs / pred, 6) if pred else None,
                z=round((obs - pred) / sd, 3) if sd else None,
                max_bucket_rel_dev=round(worst, 6),
                charged_factor=round(max(1.0, (obs + 3 * math.sqrt(obs + 1)) / pred), 6) if pred else None))
        out_routes.append(dict(route_bits=rb, buckets=len(buckets),
                               pairs_with_equal_prefix=n_pairs,
                               full_key_collisions=full,
                               full_key_positive_control=ctrl,
                               uniform_full_key_prediction=round(n_pairs / 2.0 ** key_bits, 12),
                               conditional_rate_upper_95=round(3.0 / n_pairs, 10) if n_pairs else None,
                               residual=resid))
    return dict(field_bits=p.bit_length(), p_mod_48=((p - 1) % 48 == 0), d=d,
                m_bits=m.bit_length(), m_multiple_of_6=(m % 6 == 0),
                hprimes=len(hp), hprimes_all_in_H=in_H,
                key_bits=key_bits, records=len(keys),
                pool=pool, walk_records=nw_rec,
                gap=dict(n=len(gaps), mean=round(mean_gap, 4), max_gap=max_gap,
                         max_gap_bits=max_gap.bit_length(),
                         measured_table_entries_exp=table_bits,
                         measured_table_entries=2 ** table_bits,
                         charged_table_entries_exp=32,
                         histogram=hist, tail=tail,
                         tail_measured_above_model=len(model_under),
                         tail_cells=len(tail)),
                keys=out_routes)


def main():
    env = toy_env()
    rc, mc = part_c_field()
    out = dict(partA=part_a(), partB=part_b(),
               partC=dict(part_c(rc, mc), field_r_bits=rc.bit_length(), field_d=64),
               partD=part_d(env), partE=part_e(env),
               partEcensus=[coverage_census(env, *s) for s in CENSUS_SIZES],
               partFgrid=part_f(),
               partG=part_g(pool=32, nw_rec=500, nw_gap=1200),
               toy=dict(p=TOY_P, r=TOY_R, d=TOY_D, m=TOY_M,
                        curve='y^2=x^3+b (j=0), #E=4r, order-r point, mu6 acts as (omega*x, y)'))
    json.dump(out, open('toy-summary.json', 'w'), indent=1, sort_keys=True)
    print(json.dumps(out, indent=1, sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
