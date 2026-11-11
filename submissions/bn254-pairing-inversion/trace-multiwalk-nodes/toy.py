"""Toy evidence at small scale (imported by experiment.py): the whole 16-walk algorithm in BN-family toy fields,
and the exact cross-walk coverage of the small toy field for fresh coins. Standard library only."""

import math
import random
from fractions import Fraction

from progs import NW, A_W, Tower, bn_toy, charpoly, walk, khash4, walk_params, cross_triples

TOY_FIELDS = {"small": -15, "large": -41}
TOY_SIZES = (20, 40, 60)
COV_DRAWS = 2                                        # fresh coin draws of the exact coverage count


def ceil_dp(v, d):
    return math.ceil(v * 10 ** d) / 10 ** d


def denominators_ok(rp, lam):
    """The toy's own section 1 facts: (w + 1) - lambda^k (v + 1) nonzero mod r' (w != v, or w = v and k != 0),
    and for w < v the inverses for k != k' distinct (lambda has order 12 mod r')."""
    lk = [pow(lam, k, rp) for k in range(12)]
    dens = [(A_W[w] - lk[k] * A_W[v]) % rp for w in range(NW) for v in range(NW) for k in range(12) if w != v or k]
    return all(dens) and len(set(lk)) == 12


def toy_field(name):
    q, rp, a0 = bn_toy(TOY_FIELDS[name])
    n = walk_params(rp)
    mu = Fraction(cross_triples(n), rp)
    return q, rp, a0, n, q % rp, mu - mu * mu / 2


def toy_bounds():
    """Seed-independent figures: each field's walk length and exact Bonferroni bound, and the binomial spread
    bound 0.5/sqrt(n) of the run fraction."""
    out = {}
    for f in TOY_FIELDS:
        q, rp, _, n, _, lb = toy_field(f)
        out[f] = {"r_prime": rp, "walk_length": n, "bonferroni_lower_bound_floor": math.floor(lb * 10 ** 4) / 10 ** 4}
    out["run_spread_at_most"] = {n: ceil_dp(0.5 / math.sqrt(n), 3) for n in TOY_SIZES}
    return out


def exact_coverage(rp, lam, n, deltas, betas):
    """Fraction of all x in [0, r') with some cross-walk triple (w < v, a < n, b < n, k < 12) giving
    (w + 1) x + beta_w - delta_w a = lambda^k ((v + 1) x + beta_v - delta_v b), by listing every triple."""
    S = set()
    lk = [pow(lam, k, rp) for k in range(12)]
    for w in range(NW):
        for v in range(w + 1, NW):
            for k in range(12):
                inv = pow((A_W[w] - lk[k] * A_W[v]) % rp, -1, rp)
                c0, c1, c2 = (lk[k] * betas[v] - betas[w]) * inv % rp, deltas[w] * inv % rp, -lk[k] * deltas[v] * inv % rp
                row = [(c0 + c2 * b) % rp for b in range(n)]
                for a in range(n):
                    s = c1 * a
                    S.update((s + t) % rp for t in row)
    return Fraction(len(S), rp)


def coverage(seed):
    """Exact cross-walk coverage of the small toy field for COV_DRAWS fresh coin draws; each draw's fraction is
    exact over x (no sampling), and must reach the field's Bonferroni bound."""
    q, rp, _, n, lam, lb = toy_field("small")
    if not denominators_ok(rp, lam):
        return {"denominators_nonzero": False, "every_draw_at_least_bound": False}
    rr = random.Random(2 * 10 ** 6 + seed)
    lbf = math.floor(lb * 10 ** 4) / 10 ** 4
    fr = []
    for _ in range(COV_DRAWS):
        deltas, betas = [rr.randrange(1, rp) for _ in range(NW)], [rr.randrange(rp) for _ in range(NW)]
        fr.append(exact_coverage(rp, lam, n, deltas, betas))
    return {"field": "small", "r_prime": rp, "walk_length": n, "cross_triples": cross_triples(n),
            "denominators_nonzero": denominators_ok(rp, lam),
            "coin_draws": COV_DRAWS, "exact_fractions_floor": [math.floor(f * 10 ** 5) / 10 ** 5 for f in fr],
            "bonferroni_lower_bound_floor": lbf, "every_draw_at_least_bound": all(f >= lb for f in fr)}


def toy(seed, targets, field="small", relays=4, homes_per_relay=4, cap=None):
    """The whole algorithm per target: 16 walks with fresh coins, hashed homes with a cap, sort, every
    equal-trace run expanded, every pair rechecked against all 12 conjugates (first match), x kept only if
    g^x = z. Ground truth from exponents: some pair of records (any walks) in one Frobenius class."""
    q, rp, a0, n, lam, lb = toy_field(field)
    assert denominators_ok(rp, lam), "a walk-pair denominator vanishes in the toy field"
    TW = Tower(q, a0)
    rng = random.Random(10 ** 6 + seed)
    while True:
        y = tuple(tuple((rng.randrange(q), rng.randrange(q)) for _ in range(3)) for _ in range(2))
        g = TW.pow(y, (q ** 12 - 1) // rp)
        if g != TW.one:
            break
    M = NW * n
    capd = cap or 2 * -(-M // (relays * homes_per_relay))
    gt = []                                          # g^(j 16^i): a power of g costs one product per 4 bits
    base = g
    for _ in range((rp.bit_length() + 3) // 4):
        row = [TW.one]
        for _ in range(15):
            row.append(TW.mul(row[-1], base))
        gt.append(row); base = TW.mul(row[-1], base)

    def gpow(e):
        acc, i = TW.one, 0
        while e:
            if e & 15:
                acc = TW.mul(acc, gt[i][e & 15])
            e >>= 4; i += 1
        return acc
    lk = [pow(lam, k, rp) for k in range(12)]
    solved = wrong = conj_t = cross_t = overflow = falsep = agree_clean = clean = 0
    within_t = within_solved = agree_within = 0
    for _ in range(targets):
        x = rng.randrange(rp)
        z = gpow(x)
        zp = [TW.one]
        for _ in range(NW):
            zp.append(TW.mul(zp[-1], z))
        deltas, betas = [rng.randrange(1, rp) for _ in range(NW)], [rng.randrange(rp) for _ in range(NW)]
        s1, s2 = rng.getrandbits(64), rng.getrandbits(64)
        homes, flag = {}, False
        for w in range(NW):
            D = gpow(rp - deltas[w])
            init, cur = [], TW.mul(zp[A_W[w]], gpow(betas[w]))
            for _ in range(12):
                init.append(12 * cur[0][0][0] % q); cur = TW.mul(cur, D)
            for idx, t in enumerate(walk(charpoly(TW, D, q), init, n, q)):
                Hl = homes.setdefault((khash4(t, s1) % relays, khash4(t, s2) % homes_per_relay), [])
                if len(Hl) >= capd:
                    flag = True
                    continue
                Hl.append((t, w, idx))
        cls, conj, cross, within = {}, False, False, False   # ground truth at the exponent level
        for w in range(NW):
            for idx in range(n):
                e = (A_W[w] * x + betas[w] - deltas[w] * idx) % rp
                c = min(e * l_ % rp for l_ in lk)
                if c in cls:
                    conj = True
                    cross |= any(w2 != w for w2 in cls[c])
                    within |= w in cls[c]
                cls.setdefault(c, set()).add(w)
        found = found_w = None                           # found_w: negative control, pairs from one walk only
        if not flag:
            for Hl in homes.values():
                Hl.sort()
                i = 0
                while i < len(Hl) - 1:
                    j = i
                    while j + 1 < len(Hl) and Hl[j + 1][0] == Hl[i][0]:
                        j += 1
                    for a in range(i, j + 1):
                        for b in range(a + 1, j + 1):
                            (_, wa, na), (_, wb, nb) = Hl[a], Hl[b]
                            ea, eb = (betas[wa] - deltas[wa] * na) % rp, (betas[wb] - deltas[wb] * nb) % rp
                            u = TW.mul(zp[A_W[wa]], gpow(ea))
                            yv = TW.mul(zp[A_W[wb]], gpow(eb))
                            hit = [c_ for c_ in range(12) if u == TW.frob_k(yv, c_)]
                            if not hit:
                                falsep += 1
                                continue
                            xx = (lk[hit[0]] * eb - ea) * pow((A_W[wa] - lk[hit[0]] * A_W[wb]) % rp, -1, rp) % rp
                            ok = gpow(xx) == z
                            if found is None and ok:
                                found = xx
                            if found_w is None and wa == wb and ok:
                                found_w = xx
                    i = j + 1
        overflow += flag
        conj_t += conj
        cross_t += cross
        solved += found is not None
        wrong += found is not None and found != x
        if not flag:
            clean += 1
            agree_clean += (found == x) == conj
            within_t += within
            within_solved += found_w == x
            agree_within += (found_w == x) == within
    return {"field": field, "q": q, "r_prime": rp, "walks": NW, "walk_length": n, "relays": relays,
            "denominators_nonzero": denominators_ok(rp, lam),
            "homes_per_relay": homes_per_relay, "cap_per_home": capd,
            "algorithm": {"targets": targets, "covered_cross_walk": cross_t, "covered_any_pair": conj_t,
                          "solved": solved, "wrong_outputs": wrong, "runs_without_overflow": clean,
                          "solved_equals_ground_truth_without_overflow": agree_clean == clean,
                          "overflow_runs": overflow, "pairs_failing_all_12_conjugates": falsep,
                          "one_walk_pairs_only": {"covered_inside_one_walk": within_t, "solved": within_solved,
                                                  "solved_equals_ground_truth": agree_within == clean},
                          "run_fraction_cross_walk": round(cross_t / targets, 4),
                          "run_fraction_spread_at_most": ceil_dp(0.5 / math.sqrt(targets), 3)},
            "coverage": coverage(seed)}
