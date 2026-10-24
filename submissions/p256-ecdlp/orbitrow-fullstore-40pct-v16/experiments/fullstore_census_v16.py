"""v16 census: PRIME-order fixture at the production walk shape.

v15 censused the composite-order scan fixture (1000839 = 3*333613) and measured
178/600 = 0.297, below the 0.39 floor.  The v15 note blamed the ~140-step toy
recurrence; that cannot explain an 8-step-per-walker census.  The actual cause is
measurable: the collision RATE matched the birthday law (178 successes + 200
gcd-nonunit denominators = 378 relations vs 384 predicted), but the composite
order destroys every relation whose denominator shares a factor with n.

v16 therefore measures on y^2 = x^3 - 3x + 1 over p = 1000003, whose group order
is the PRIME 999979 (Miller-Rabin in v16kit), so no nonunit denominators exist --
the same arithmetic regime as P-256's prime order.  It censuses three walk shapes
at the same lambda (all V = 800 visits: 100x8, 200x4, 50x16), each with per-walker
seeds and exact Clopper-Pearson 95% bounds, and reports the collision count against
the birthday prediction as a direct uniformity check of this table+selector.
"""
import json, math, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v16kit as K

R = 64
SHAPES = [(100, 8, 4000), (200, 4, 2000), (50, 16, 2000)]
SECRET = 777001
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")


def cp95_lower(k, n):
    if k == 0:
        return 0.0
    def tail(p):
        # stable log-space binomial upper tail
        lp, lq = math.log(p), math.log1p(-p)
        s = 0.0
        for i in range(k, n + 1):
            s += math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                          + i * lp + (n - i) * lq)
        return s
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if tail(mid) < 0.025:
            lo = mid
        else:
            hi = mid
    return lo


def trial(P, Q, tab, ti, W, STEPS):
    store = {}
    succ = False
    hits = noprog = gcds = visits = 0
    first_bucket = [0] * 64
    for w in range(W):
        rng = random.Random((20251024 << 24) ^ (ti << 12) ^ (w * 0x9E3779B97F4A7C15) ^ (STEPS << 40))
        c, d = rng.randrange(1, K.N_PRIME), rng.randrange(1, K.N_PRIME)
        S = K.add(K.mul(c, P), K.mul(d, Q))
        while S is None:
            c, d = rng.randrange(1, K.N_PRIME), rng.randrange(1, K.N_PRIME)
            S = K.add(K.mul(c, P), K.mul(d, Q))
        first_bucket[S[0] % 64] += 1
        trail = ti * W + w
        for t in range(STEPS):
            r = ((S[0] ^ (t * 2654435761) ^ (trail * 7919)) % R) if S else 0
            a, b, V = tab[r]
            T = K.add(S, V)
            if T is None:
                break                            # identity hit: walker restarts, production rule
            c2, d2 = (c + a) % K.N_PRIME, (d + b) % K.N_PRIME
            C = K.canon(T)
            if C != T:
                c2, d2 = (-c2) % K.N_PRIME, (-d2) % K.N_PRIME
            visits += 1
            key = K.cls(C)
            if key in store:
                hits += 1
                pc, pd, py = store[key]
                if py != C[1]:
                    da, db = (c2 + pc) % K.N_PRIME, (d2 + pd) % K.N_PRIME
                else:
                    da, db = (c2 - pc) % K.N_PRIME, (d2 - pd) % K.N_PRIME
                if da == 0 and db == 0:
                    noprog += 1
                elif math.gcd(db, K.N_PRIME) != 1:
                    gcds += 1
                else:
                    cand = (-da * K.inv(db, K.N_PRIME)) % K.N_PRIME
                    if K.mul(cand, P) == Q:
                        succ = True
            else:
                store[key] = (c2, d2, C[1])
            S = C
            c, d = c2, d2
    return succ, visits, hits, noprog, gcds, first_bucket


def main():
    P = K.base_point(K.B_PRIME)
    Q = K.mul(SECRET, P)
    tab = K.row_table(P, Q, R)
    classes = (K.N_PRIME + 1) // 2
    shapes_ev = []
    for W, STEPS, N in SHAPES:
        V = W * STEPS
        lam = V * V / (2.0 * classes)
        wins = relations = noprog = gcds = visits = 0
        bucket = [0] * 64
        for ti in range(N):
            s, v, h, npg, g, fb = trial(P, Q, tab, ti, W, STEPS)
            wins += 1 if s else 0
            relations += h - npg
            noprog += npg
            gcds += g
            visits += v
            bucket = [x + y for x, y in zip(bucket, fb)]
        lo = cp95_lower(wins, N)
        pred = 1.0 - math.exp(-lam)
        exp_coll = N * V * (V - 1) / (2.0 * classes)
        exp_first = visits and 0 or 0
        # chi-square of first-visit x mod 64 vs uniform
        expb = sum(bucket) / 64.0
        chi2 = sum((b - expb) ** 2 / expb for b in bucket)
        shapes_ev.append({
            "walkers": W, "steps_per_walker": STEPS, "trials": N,
            "visits_per_trial": V, "lambda": lam,
            "verified_successes": wins, "measured_success": wins / N,
            "cp95_lower": lo, "birthday_predicted": pred,
            "relation_collisions": relations, "birthday_collisions": exp_coll,
            "collision_ratio": relations / exp_coll if exp_coll else None,
            "no_progress": noprog, "gcd_nonunit": gcds,
            "first_visit_chi2_mod64": chi2, "chi2_df": 63,
            "chi2_99pct_cutoff": 93.0, "meets_floor": lo >= 0.39,
        })
        print(json.dumps(shapes_ev[-1], indent=1))
    ev = {
        "fixture": {"p": K.P_TOY, "b": K.B_PRIME, "order": K.N_PRIME,
                    "order_prime": K.is_probable_prime(K.N_PRIME),
                    "curve": "y^2 = x^3 - 3x + 1 over p=1000003",
                    "why": "prime group order, so the composite-order nonunit loss of "
                           "the v15 scan fixture (1000839=3*333613) is absent; same "
                           "arithmetic regime as P-256's prime order",
                    "classes": classes,
                    "class_semantics": "cls = canonical x (negate iff y odd); P/-P merged"},
        "config": {"rows": R, "shapes": [f"{w}x{s}" for w, s, _ in SHAPES],
                   "detector": "insert-time full store keyed on canonical x; unconditional product; verified scalar",
                   "seeds": "each walker: (20251024<<24) ^ (trial<<12) ^ (w*gamma) ^ (steps<<40)"},
        "shapes": shapes_ev,
        "method": "exact Clopper-Pearson 95% lower bounds; collision count compared with "
                  "the birthday prediction V(V-1)/(2C) as a direct uniformity check of the "
                  "shipped table+selector at the production walk shape",
    }
    os.makedirs(OUT, exist_ok=True)
    json.dump(ev, open(os.path.join(OUT, "census_v16.json"), "w"), indent=1, sort_keys=True)
    print("wrote census_v16.json")


if __name__ == "__main__":
    main()
