"""v15 census: production-shaped shared full store, per-worker seeds, unconditional
product, exact Clopper-Pearson 95% lower bound.

WHY MULTI-WALKER.  Production is 2^20 walkers of 2^107.67 steps sharing one store;
the birthday variable is the TOTAL visit count V.  A single 563-step toy walk would
already sit at its rho-cycle scale (the toy functional graph cycles in ~sqrt(500k)
= 443 steps), which is not the production regime (2^107.67 << 2^127.5 per-walker
rho length).  The census therefore reproduces the production shape: W independent
walkers, short per-walker trails, one shared insert-time store, total V fixed so
lambda = V^2/(2C) = 0.6329.  Each walker has its own seed; every repeat yields a
candidate scalar unconditionally and counts only if cand*P == Q.
"""
import json, math, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v15kit as K

R = 64
N_TRIALS = 600
WALKERS = 100
STEPS = 8
SECRET = 1234567
INV_TRIALS = 649969
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")


def count_points():
    n = 1
    for x in range(K.P_TOY):
        rhs = (x * x * x - 3 * x + K.B_TOY) % K.P_TOY
        if rhs == 0:
            n += 1
        elif pow(rhs, (K.P_TOY - 1) // 2, K.P_TOY) == 1:
            n += 2
    return n


def cp95_lower(k, n):
    """p with P(X>=k | p) = 0.025; the binomial tail is increasing in p."""
    if k == 0:
        return 0.0
    def tail(p):
        return sum(math.comb(n, i) * p ** i * (1.0 - p) ** (n - i) for i in range(k, n + 1))
    lo, hi = 0.0, 1.0
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if tail(mid) < 0.025:
            lo = mid
        else:
            hi = mid
    return lo


def trial(P, Q, tab, ti):
    store = {}
    succ = False
    visits = noprog = gcds = 0
    for w in range(WALKERS):
        rng = random.Random((20251024 << 24) ^ (ti << 12) ^ (w * 0x9E3779B97F4A7C15))
        c, d = rng.randrange(1, K.N_TOY), rng.randrange(1, K.N_TOY)
        S = K.add(K.mul(c, P), K.mul(d, Q))
        trail = ti * WALKERS + w
        for t in range(STEPS):
            r = ((S[0] ^ (t * 2654435761) ^ (trail * 7919)) % R) if S else 0
            a, b, V = tab[r]
            T = K.add(S, V)
            c2, d2 = (c + a) % K.N_TOY, (d + b) % K.N_TOY
            C = K.canon(T)
            if C != T:
                c2, d2 = (-c2) % K.N_TOY, (-d2) % K.N_TOY
            visits += 1
            key = K.cls(C)
            if key in store:
                pc, pd, py = store[key]
                if py != C[1]:                    # stored representative is -C
                    da, db = (c2 + pc) % K.N_TOY, (d2 + pd) % K.N_TOY
                else:
                    da, db = (c2 - pc) % K.N_TOY, (d2 - pd) % K.N_TOY
                if da == 0 and db == 0:
                    noprog += 1
                elif math.gcd(db, K.N_TOY) != 1:
                    gcds += 1
                else:
                    cand = (-da * K.inv(db, K.N_TOY)) % K.N_TOY   # unconditional product
                    if K.mul(cand, P) == Q:                        # verified scalar
                        succ = True
            else:
                store[key] = (c2, d2, C[1])
            S = C
            c, d = c2, d2
    return succ, visits, noprog, gcds


def main():
    P = K.base_point()
    Q = K.mul(SECRET, P)
    tab = K.row_table(P, Q, R)
    npoints = count_points()
    classes = (npoints + 1) // 2 + 1
    V = WALKERS * STEPS
    lam = V * V / (2.0 * classes)
    wins = noprog = gcds = visits = 0
    for ti in range(N_TRIALS):
        s, v, npg, g = trial(P, Q, tab, ti)
        wins += 1 if s else 0
        visits += v
        noprog += npg
        gcds += g
    lo = cp95_lower(wins, N_TRIALS)
    rng = random.Random(99)
    misses = nonunits = 0
    for _ in range(INV_TRIALS):
        a = rng.randrange(1, K.N_TOY)
        if math.gcd(a, K.N_TOY) != 1:
            nonunits += 1
            continue
        try:
            x = K.inv(a, K.N_TOY)
        except ValueError:
            misses += 1
            continue
        if (a * x) % K.N_TOY != 1:
            misses += 1
    m = K.N_TOY
    factor = []
    q = 2
    while q * q <= m:
        while m % q == 0:
            factor.append(q)
            m //= q
        q += 1
    if m > 1:
        factor.append(m)
    ev = {
        "fixture": {"p": K.P_TOY, "b": K.B_TOY, "order": K.N_TOY,
                    "order_prime": K.is_probable_prime(K.N_TOY),
                    "order_factorization": " * ".join(str(f) for f in factor),
                    "v14_order_claim": "1000839 asserted prime; it is 3 * 333613",
                    "inverse_method": "extended Euclid (pow(a,-1,n)); no Fermat claim",
                    "inverse_trials": INV_TRIALS, "inverse_misses": misses,
                    "inverse_nonunits_skipped": nonunits,
                    "class_semantics": "cls = canonical x (y even); P and -P collide, x and p-x are distinct keys",
                    "classes": classes, "points": npoints},
        "config": {"rows": R, "trials": N_TRIALS, "walkers_per_trial": WALKERS,
                   "steps_per_walker": STEPS, "visits_per_trial": V,
                   "visits_total": visits, "lambda_fixture": lam,
                   "detector": "shared insert-time full store keyed on cls; unconditional product",
                   "seeds": "each walker: seed = (20251024<<24) ^ (trial<<12) ^ (w * 0x9E3779B97F4A7C15)"},
        "counts": {"verified_successes": wins, "trials": N_TRIALS,
                   "measured_success": wins / N_TRIALS, "cp95_lower": lo,
                   "no_progress_repeats": noprog, "gcd_nonunit_repeats": gcds},
        "birthday_check": {"predicted": 1.0 - math.exp(-lam),
                           "measured": wins / N_TRIALS,
                           "note": "fixture-scope test of the same birthday law the production V is sized by"},
        "meets_floor": lo >= 0.39,
    }
    os.makedirs(OUT, exist_ok=True)
    json.dump(ev, open(os.path.join(OUT, "census_v15.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({k: ev[k] for k in ("fixture", "config", "counts", "birthday_check", "meets_floor")}, indent=1))


if __name__ == "__main__":
    main()
