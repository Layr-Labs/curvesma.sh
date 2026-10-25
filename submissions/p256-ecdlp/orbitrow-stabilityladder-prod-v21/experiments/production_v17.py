"""v17 selector-period census ladder -- the decisive floor measurement.

WHY THIS EXISTS
v16's census fixed R=64, one secret, and trails of 4-16 steps.  A selector of R
rows addressed by a step-counter hash is periodic in the step counter, so a
4-16 step trail exercises only 4-16 of 64 rows and does not cross even one
period.  Production trails are 2^107.67 steps at R=2^20, i.e. 2^87.67 periods.
The v16 ticket therefore could not tell whether its 0.36-class success at 64/128
steps was the selector's periodicity or a small-table artifact.

WHAT THIS MEASURES
A ladder in R with the period coverage held at >= 16 full periods at every rung:
S = 16*R steps per walker, R in {64, 256, 1024, 4096, 2^20}, many random
secrets (one per trial), lambda held at 0.634 by sizing the group to
q ~ V^2/lambda.  The question is whether success-per-trial is stable in R at
fixed coverage.  A coverage sweep at R=64 (0.25 .. 64 periods) is run as the
control that coverage is the variable that matters.

FIXTURE
The walk lives in a prime-order cyclic group Z/q.  Any two cyclic groups of
prime order are isomorphic, so the walk's group arithmetic (prefix sums of
table rows) is the same computation as on the P-256 prime-order subgroup.  The
two things a curve supplies are replaced explicitly:
  * sign classes: canon(v) = min(v, q-v) merges v and -v, exactly as the
    shipped rule merges P and -P; classes = (q+1)/2, and the coefficient sign
    flips when the representative is negated;
  * the row selector hashes the element's representative integer instead of the
    x-coordinate -- a hash of the element, as x(S) is.
The useful-relation predicate is NOT assumed: a collision yields (da, db) with
da + db*x = 0 mod q.  It is useful iff (da, db) != (0, 0) and db != 0; then
x = -da/db and the trial is a verified success.  This is the same predicate the
curve census uses, and it is why this census can report "at least one relation
that recovers" rather than a mean collision count.

Cost: steps * ~1.2 us in CPython, so cells are sized to the 55-minute budget and
run on a process pool.
"""
import json, math, os, random, sys
import multiprocessing as mp
mp.set_start_method("fork", force=True)
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")
LAMBDA = 0.634
ZERO_OPT = 0


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
    """Exact Clopper-Pearson 95% one-sided lower bound."""
    if k == 0:
        return 0.0
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


def row_coeffs(r, R):
    return 1 + (r * 7 + 3) % R, 1 + (r * 11 + 5) % R


def make_table(q, x, R):
    """Row r = a*P + b*Q with P=1, Q=x, so the group element is a + b*x mod q."""
    tab = []
    for r in range(R):
        a, b = row_coeffs(r, R)
        tab.append(((a + b * x) % q, a, b))
    return tab


def trial(args):
    """One trial = one fresh secret x, W walkers of STEPS steps, one GLOBAL store.

    The store is global for the whole trial (all walkers, all batches), which is
    the detector the ledger prices: a collision is any key seen before anywhere
    in the challenge, not just inside one batch.
    """
    q, R, W, STEPS, ti = args
    rng = random.Random((20251024 << 26) ^ (ti * 0x9E3779B97F4A7C15) ^ (R << 3) ^ (W << 50))
    x = rng.randrange(2, q - 1)
    tab = make_table(q, x, R)
    half = q >> 1
    store = {}
    coll = trivial = useful = 0
    visits = 0
    ident = 0
    cyc = None
    tmul = 2654435761
    track = STEPS <= (1 << 17)
    for w in range(W):
        c = rng.randrange(1, q)
        d = rng.randrange(1, q)
        S = (c + d * x) % q
        trail = (ti * W + w) & 0xFFFFFF
        seen = {} if track else None          # exact-state repeat = cycle (coefficients repeat)
        for t in range(STEPS):
            prod = (S ^ ((t * tmul) & 0xFFFFFFFFFFFFFFFF) ^ (trail * 7919)) % R
            v, a, b = tab[prod]
            T = S + v
            if T >= q:
                T -= q
            if T == 0:
                ident += 1
                break                       # identity hit: production restart rule
            c2 = c + a
            if c2 >= q:
                c2 -= q
            d2 = d + b
            if d2 >= q:
                d2 -= q
            # The canonical representative is a STORAGE rule, not a walk rule.
            # v16 canonicalised the walk state itself (S = canon(T)); then two
            # consecutive equal rows are a permanent 2-cycle: canon(canon(S+v)+v)=S
            # and the coefficients return exactly, so the walker can never make
            # progress again.  The ladder measures that (every v16 trial cycles at
            # gap 2).  Here the walk steps on the true point and only the KEY and
            # the stored coefficients are canonicalised.
            key, kc, kd = T, c2, d2
            if key > half:
                key = q - key
                kc = (-kc) % q
                kd = (-kd) % q
            visits += 1
            if track:
                k = (c2, d2)
                if k in seen:
                    if cyc is None or (t - seen[k]) < cyc:
                        cyc = t - seen[k]
                    break                   # exact true-state repeat: the walk cycles
                seen[k] = t
            got = store.get(key)
            if got is None:
                store[key] = (kc, kd)
            else:
                coll += 1
                pc, pd = got
                da = (kc - pc) % q
                db = (kd - pd) % q
                if da == 0 and db == 0:
                    trivial += 1
                else:
                    cand = (-da * pow(db, -1, q)) % q if db else None
                    if cand == x:
                        useful += 1
            S = T
            c = c2
            d = d2
    return {
        "success": 1 if useful else 0,
        "useful_relations": useful,
        "collisions": coll,
        "trivial": trivial,
        "visits": visits,
        "identity_breaks": ident,
        "cycle_len": cyc,
    }


def make_cell(R, W, S, trials, label, notes=""):
    V = W * S
    q = next_prime(max(1 << 8, int(round(V * V / LAMBDA)))) | 1
    if q % 2 == 0:
        q = next_prime(q + 1)
    classes = (q + 1) // 2
    lam = V * (V - 1) / (2.0 * classes)
    periods = S / float(R)
    return {
        "label": label, "R": R, "walkers": W, "steps_per_walker": S,
        "visits_per_trial": V, "trials": trials, "q": q, "q_bits": q.bit_length(),
        "classes": classes, "lambda": lam,
        "periods_covered": periods,
        "coverage_ok": periods >= 1.0,
        "coverage_meets_16": periods >= 16.0,
        "coverage_flag": "OK" if periods >= 16.0 else ("FLAGGED_SUBPERIOD" if periods < 1.0 else "PARTIAL"),
        "notes": notes,
    }


CELLS = [
    # ladder: 16 full selector periods at every rung, lambda held at 0.634
    make_cell(64,      64, 16 * 64,        600, "R64_W64_16p", "ladder rung 1"),
    make_cell(256,     64, 16 * 256,       400, "R256_W64_16p", "ladder rung 2"),
    make_cell(1024,    64, 16 * 1024,      300, "R1024_W64_16p", "ladder rung 3"),
    make_cell(4096,    64, 16 * 4096,      240, "R4096_W64_16p", "ladder rung 4"),
    make_cell(1 << 20, 1, 16 * (1 << 20),   64, "R2p20_W1_16p", "ladder rung 5, W=1 (conservative within-trail fraction)"),
    # coverage sweep at R=64: the control that coverage is the operative variable
    make_cell(64, 64, 16,    600, "R64_W64_0p25", "coverage 0.25 period"),
    make_cell(64, 64, 64,    600, "R64_W64_1p",    "coverage 1 period"),
    make_cell(64, 64, 256,   600, "R64_W64_4p",    "coverage 4 periods"),
    make_cell(64, 64, 1024,  600, "R64_W64_16p",   "coverage 16 periods"),
    make_cell(64, 64, 4096,  600, "R64_W64_64p",   "coverage 64 periods"),
    # walker-count control at fixed R and coverage: does W matter?
    make_cell(64, 1,  1024,  600, "R64_W1_16p",    "W=1"),
    make_cell(64, 4,  1024,  600, "R64_W4_16p",    "W=4"),
    make_cell(64, 16, 1024,  600, "R64_W16_16p",   "W=16"),
]


PROCS = {"R2p20_W1_16p": 6, "R4096_W64_16p": 10, "R1024_W64_16p": 16}


def run_cell(cell, procs=None):
    procs = procs or PROCS.get(cell["label"], 24)
    args = [(cell["q"], cell["R"], cell["walkers"], cell["steps_per_walker"], ti)
            for ti in range(cell["trials"])]
    with Pool(procs) as pool:
        res = pool.map(trial, args, chunksize=max(1, len(args) // (procs * 4)))
    wins = sum(r["success"] for r in res)
    cyc = [r.get("cycle_len") for r in res if r.get("cycle_len")]
    coll = sum(r["collisions"] for r in res)
    triv = sum(r["trivial"] for r in res)
    useful = sum(r["useful_relations"] for r in res)
    visits = sum(r["visits"] for r in res)
    n = len(res)
    lo = cp95_lower(wins, n)
    # birthday collision prediction uses the ACTUAL visits, not the nominal V
    exp_coll = sum(r["visits"] * (r["visits"] - 1) / (2.0 * cell["classes"]) for r in res)
    out = dict(cell)
    out.update({
        "verified_successes": wins,
        "measured_success": wins / n,
        "cp95_lower": lo,
        "meets_floor": lo >= 0.39,
        "collisions": coll,
        "trivial_collisions": triv,
        "useful_relations": useful,
        "birthday_collisions": exp_coll,
        "collision_ratio": coll / exp_coll if exp_coll else None,
        "useful_per_collision": useful / coll if coll else None,
        "visits_total": visits,
        "identity_breaks": sum(r["identity_breaks"] for r in res),
        "cycle_lengths_min": min(cyc) if cyc else None,
        "cycle_lengths_max": max(cyc) if cyc else None,
        "trials_with_cycle": len(cyc),
        "seed_rule": "per trial: secret and per-walker seeds from Random((20251024<<26) ^ ti*golden ^ R<<3 ^ W<<50); fresh secret every trial",
    })
    return out


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    os.makedirs(OUT, exist_ok=True)
    if only == "--table":
        print(json.dumps([{k: c[k] for k in ("label", "R", "walkers", "steps_per_walker",
                                             "trials", "q_bits", "lambda",
                                             "periods_covered", "coverage_ok")} for c in CELLS], indent=1))
        return
    results = []
    for cell in CELLS:
        if only and only not in cell["label"]:
            continue
        path0 = os.path.join(OUT, "ladder_%s.json" % cell["label"])
        if os.path.exists(path0) and os.environ.get("FORCE") != "1":
            results.append(json.load(open(path0)))
            continue
        r = run_cell(cell)
        print(json.dumps({k: r[k] for k in ("label", "R", "walkers", "steps_per_walker",
                                            "periods_covered", "measured_success",
                                            "cp95_lower", "collision_ratio",
                                            "useful_per_collision")}, indent=1), flush=True)
        path = os.path.join(OUT, "ladder_%s.json" % cell["label"])
        json.dump(r, open(path, "w"), indent=1, sort_keys=True)
        results.append(r)
    json.dump({"cells": results, "lambda_target": LAMBDA,
               "question": "does success stabilize in R with period coverage held >=16?",
               "config_table": [{k: c[k] for k in ("label", "R", "walkers", "steps_per_walker",
                                                   "trials", "periods_covered", "coverage_ok",
                                                   "coverage_meets_16")} for c in CELLS]},
              open(os.path.join(OUT, "period_ladder_v17.json"), "w"), indent=1, sort_keys=True)
    print("wrote period_ladder_v17.json")


if __name__ == "__main__":
    main()
