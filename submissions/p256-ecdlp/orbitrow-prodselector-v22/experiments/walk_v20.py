"""v20 stability-ladder runner: the PRODUCTION selector measured across period
counts and walker counts, with collision-provenance tagging.

WHY (verdict d32b4f73, flaw 1)
The reviewer accepted the fixture arithmetic but rejected the transfer:
"the measured cell uses Z/q with q=443966840237629, one walker and 2^24 steps,
not P-256 ... even the submission's birthday model gives approximately 0.634
collisions per fixture trail but only (2^107.67)^2/r ~ 2^-40.66 per production
trail.  Production success would therefore depend predominantly on cross-trail
collisions, not the asserted conservative per-walker comparison."

So there are two axes to measure, not one:
  * PERIOD axis -- does success per trial stay flat as the trail covers
    16 -> 1024 full selector periods (per-100-periods hazard)?
  * MIX axis -- production's useful relations come from collisions between
    DISTINCT walkers (2^20 walkers, one shared store, within-walker load
    ~2^-40.66).  This runner tags every collision with the walker that wrote
    the record, so `success_cross` (a useful relation from two distinct
    walkers) is measured separately from `success` (any useful relation) and
    from `success_within` (relations whose two records came from one walker).
    Within-walker load is lambda/W, so W = 64 leaves 1.6% of collisions
    within-walker, two orders of magnitude closer to production than W = 1.

The selector is the PRODUCTION selector, not a representative-based stand-in:
the same source lines (hashed, control (a)) evaluate
    prod = (S ^ (t*2654435761 mod 2^64) ^ (trail*7919)) mod R
on the same true-point walk, at every rung of the ladder.  Nothing about the
selector changes with R, W or the period count; only the modulus R does, which
is what production also does (R = 2^20).

Lambda is held at 0.634 at every cell by sizing q ~ V^2/lambda, exactly as the
v17/v19 production cell did, so p_hat is comparable across cells: if the
per-period hazard is stable, p_hat is flat in coverage and in W.
"""
import json, math, os, random, sys
import multiprocessing as mp
from multiprocessing import Pool

mp.set_start_method("fork", force=True)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")
LAMBDA = 0.634
TMUL = 2654435761
GOLDEN = 0x9E3779B97F4A7C15
TAG = 1024                     # walker tag radix (production W = 2^20 is not run)
CROSS = "cross"                # collision between two distinct walkers
WITHIN = "within"              # collision inside one walker's own trail


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
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


def next_prime(n):
    n |= 1
    while not is_prime(n):
        n += 2
    return n


def cp95_lower(k, n):
    """Exact Clopper-Pearson one-sided 95% lower bound (bisection on the tail)."""
    if k == 0:
        return 0.0
    if k == n:
        return (0.025) ** (1.0 / n)

    def tail(p):
        lp, lq = math.log(p), math.log1p(-p)
        s = 0.0
        lg = math.lgamma
        for i in range(k, n + 1):
            s += math.exp(lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq)
        return s

    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if tail(mid) < 0.025:
            lo = mid
        else:
            hi = mid
    return lo


def cp95_upper(k, n):
    """Exact Clopper-Pearson one-sided 95% upper bound."""
    if k == n:
        return 1.0

    def tail(p):
        lp, lq = math.log(p), math.log1p(-p)
        s = 0.0
        lg = math.lgamma
        for i in range(0, k + 1):
            s += math.exp(lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq)
        return s

    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if tail(mid) > 0.025:
            lo = mid
        else:
            hi = mid
    return hi


def wald_lower(k, n, z=1.959963984540054):
    p = k / n
    return p - z * math.sqrt(p * (1 - p) / n)


def row_coeffs(r, R):
    return 1 + (r * 7 + 3) % R, 1 + (r * 11 + 5) % R


def make_table(q, x, R):
    tab = []
    for r in range(R):
        a, b = row_coeffs(r, R)
        tab.append(((a + b * x) % q, a, b))
    return tab


def trial(args):
    """One trial = one fresh secret x, W walkers of STEPS steps, one GLOBAL store.

    Identical walk, selector and store rules to the v19 production cell; the only
    change is that each stored record also carries the id of the walker that wrote
    it, so a collision is attributed to CROSS or WITHIN.
    """
    q, R, W, STEPS, ti = args
    rng = random.Random((20251024 << 26) ^ (ti * GOLDEN) ^ (R << 3) ^ (W << 50))
    x = rng.randrange(2, q - 1)
    tab = make_table(q, x, R)
    half = q >> 1
    store = {}
    coll = within_coll = cross_coll = 0
    trivial = useful = useful_within = useful_cross = 0
    success_cross = 0
    visits = 0
    ident = 0
    tmul = 2654435761
    for w in range(W):
        c = rng.randrange(1, q)
        d = rng.randrange(1, q)
        S = (c + d * x) % q
        trail = (ti * W + w) & 0xFFFFFF
        for t in range(STEPS):
            prod = (S ^ ((t * tmul) & 0xFFFFFFFFFFFFFFFF) ^ (trail * 7919)) % R
            v, a, b = tab[prod]
            T = S + v
            if T >= q:
                T -= q
            if T == 0:
                ident += 1
                break
            c2 = c + a
            if c2 >= q:
                c2 -= q
            d2 = d + b
            if d2 >= q:
                d2 -= q
            key, kc, kd = T, c2, d2
            if key > half:
                key = q - key
                kc = (-kc) % q
                kd = (-kd) % q
            visits += 1
            # packed record: (kc*q + kd)*TAG + writer id  -- one int, not a tuple
            got = store.get(key)
            if got is None:
                store[key] = (kc * q + kd) * TAG + w
            else:
                coll += 1
                gw = got % TAG
                gr = got // TAG
                pd = gr % q
                pc = gr // q
                da = (kc - pc) % q
                db = (kd - pd) % q
                same = 1 if gw == w else 0
                if same:
                    within_coll += 1
                else:
                    cross_coll += 1
                if da == 0 and db == 0:
                    trivial += 1
                else:
                    cand = (-da * pow(db, -1, q)) % q if db else None
                    if cand == x:
                        useful += 1
                        if same:
                            useful_within += 1
                        else:
                            useful_cross += 1
                            success_cross = 1
            S = T
            c = c2
            d = d2
    return {
        "success": 1 if useful else 0,
        "success_cross": success_cross,
        "success_within": 1 if useful_within else 0,
        "useful_relations": useful,
        "useful_cross": useful_cross,
        "useful_within": useful_within,
        "collisions": coll,
        "cross_collisions": cross_coll,
        "within_collisions": within_coll,
        "trivial": trivial,
        "visits": visits,
        "identity_breaks": ident,
    }


def make_cell(R, W, periods, trials, label, note="", procs=None):
    S = int(round(periods * R))
    V = W * S
    q = next_prime(max(1 << 8, int(round(V * V / LAMBDA)))) | 1
    if q % 2 == 0:
        q = next_prime(q + 1)
    classes = (q + 1) // 2
    lam = V * (V - 1) / (2.0 * classes)
    within_load = (S * (S - 1)) / (2.0 * classes)
    return {
        "label": label, "R": R, "walkers": W, "steps_per_walker": S,
        "periods": periods, "visits_per_trial": V, "trials": trials,
        "q": q, "q_bits": q.bit_length(), "classes": classes, "lambda": lam,
        "within_walker_load": within_load,
        "cross_share_expected": (lam - W * within_load) / lam if lam else None,
        "note": note, "procs": procs,
    }


CELLS = [
    # -- period axis: production selector, fixed lambda, coverage 16..1024 periods
    make_cell(1024,     16, 16,   200, "P_R1024_p16",   "period rung 1", 24),
    make_cell(1024,     16, 64,   160, "P_R1024_p64",   "period rung 2", 24),
    make_cell(1024,     16, 256,   96, "P_R1024_p256",  "period rung 3", 24),
    make_cell(1024,     16, 1024,  48, "P_R1024_p1024", "1024 periods at R=1024", 20),
    make_cell(4096,      4, 16,   200, "P_R4096_p16",   "period rung 1", 24),
    make_cell(4096,      4, 64,   160, "P_R4096_p64",   "period rung 2", 24),
    make_cell(4096,      4, 256,   96, "P_R4096_p256",  "period rung 3", 24),
    make_cell(4096,      4, 1024,  48, "P_R4096_p1024", "1024 periods at R=4096", 20),
    make_cell(65536,     1, 16,   120, "P_R65536_p16",  "period rung 1", 24),
    make_cell(65536,     1, 64,    64, "P_R65536_p64",  "period rung 2", 20),
    make_cell(65536,     1, 256,   40, "P_R65536_p256", "period rung 3", 16),
    make_cell(65536,     1, 1024,  16, "P_R65536_p1024", "2^26 steps/walker, W=1", 8),
    # -- walker/mix axis at R=4096, 256 periods, lambda held
    make_cell(4096,      1, 256,   96, "W_R4096_W1",    "cross share 0", 24),
    make_cell(4096,     16, 256,   48, "W_R4096_W16",   "cross share 15/16", 20),
    make_cell(4096,     64, 256,   24, "W_R4096_W64",   "cross share 63/64", 12),
    # -- cross-only axis: W=64 keeps within-walker load at lambda/64 = 0.0099,
    #    98.4% of collisions between distinct walkers (production: ~100%)
    make_cell(256,      64, 16,  2000, "X_R256_W64_p16",   "cross-tagged, 16 periods", 24),
    make_cell(256,      64, 64,   600, "X_R256_W64_p64",   "cross-tagged, 64 periods", 24),
    make_cell(256,      64, 256,  420, "X_R256_W64_p256",  "cross-tagged, 256 periods", 24),
    make_cell(1024,     64, 16,   800, "X_R1024_W64_p16",  "R axis at W=64", 24),
    make_cell(4096,     64, 16,   400, "X_R4096_W64_p16",  "R axis at W=64", 20),
]


def run_cell(cell):
    procs = cell.get("procs") or 16
    args = [(cell["q"], cell["R"], cell["walkers"], cell["steps_per_walker"], ti)
            for ti in range(cell["trials"])]
    with Pool(procs) as pool:
        res = pool.map(trial, args, chunksize=max(1, len(args) // (procs * 4)))
    n = len(res)
    wins = sum(r["success"] for r in res)
    wins_cross = sum(r["success_cross"] for r in res)
    wins_within = sum(r["success_within"] for r in res)
    coll = sum(r["collisions"] for r in res)
    visits = sum(r["visits"] for r in res)
    out = dict(cell)
    out.update({
        "trials": n,
        "wins": wins, "measured_success": wins / n,
        "cp95_lower": cp95_lower(wins, n),
        "cp95_upper": cp95_upper(wins, n),
        "wald_lower": wald_lower(wins, n),
        "wins_cross": wins_cross, "measured_success_cross": wins_cross / n,
        "cp95_lower_cross": cp95_lower(wins_cross, n),
        "wins_within_only": wins_within - 0,
        "collisions": coll,
        "cross_collisions": sum(r["cross_collisions"] for r in res),
        "within_collisions": sum(r["within_collisions"] for r in res),
        "useful_relations": sum(r["useful_relations"] for r in res),
        "useful_cross": sum(r["useful_cross"] for r in res),
        "useful_within": sum(r["useful_within"] for r in res),
        "trivial_collisions": sum(r["trivial"] for r in res),
        "visits_total": visits,
        "identity_breaks": sum(r["identity_breaks"] for r in res),
        "outcomes": "".join(str(r["success"]) for r in res),
        "outcomes_cross": "".join(str(r["success_cross"]) for r in res),
        "seed_rule": "Random((20251024<<26) ^ ti*0x9E3779B97F4A7C15 ^ R<<3 ^ W<<50); fresh secret per trial",
        "selector": "prod = (S ^ ((t*2654435761) & (2^64-1)) ^ (trail*7919)) % R; row r = ((1+(7r+3) mod R)*P + (1+(11r+5) mod R)*Q)",
    })
    return out


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    os.makedirs(OUT, exist_ok=True)
    order = sorted(CELLS, key=lambda c: c["visits_per_trial"] * c["trials"])
    for cell in order:
        if only and only not in cell["label"]:
            continue
        path = os.path.join(OUT, "ladder20_%s.json" % cell["label"])
        if os.path.exists(path) and os.environ.get("FORCE") != "1":
            print("skip", cell["label"], flush=True)
            continue
        r = run_cell(cell)
        json.dump(r, open(path, "w"), indent=1, sort_keys=True)
        print(json.dumps({k: r[k] for k in ("label", "R", "walkers", "periods",
                                            "measured_success", "cp95_lower",
                                            "measured_success_cross",
                                            "cp95_lower_cross", "collisions",
                                            "cross_collisions", "within_collisions")}), flush=True)
        os.system("rm -rf /tmp/none")  # no-op keeps stdout flush ordering stable
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
