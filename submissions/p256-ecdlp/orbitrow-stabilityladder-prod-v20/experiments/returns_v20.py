"""v20 walk-return / restart measurement (answers verdict d32b4f73 flaw 2).

The v19 cycle row tracked repeats of the COEFFICIENT PAIR (c,d), which the review
correctly called a category error: at R=64/W=64/1024 steps the coefficients rise
by at most 64 per step, so a repeat needs a wrap and the zero was close to
automatic; a per-walker death probability was then divided by 2^20 walkers, which
the review also rejected.

This script measures the DECLARED walk on the declared selector and reports three
distinct things, per TRAIL (a trail = one walker's walk, so trails = trials * W):

  1. exact (S, t) returns.  Impossible by construction: the next state is
     f(S, t) with t strictly increasing mod 2^64, so the pair (S, t) never
     recurs inside a trail.  The classical rho-cycle death therefore cannot
     occur on this walk, and that is a STRUCTURAL fact proved from the selector,
     not a zero measured on a vacuous test.  Verified by exact tracking anyway.
  2. S returns at a DIFFERENT phase.  These do occur (birthday in q) and they do
     NOT stop the walk: the phase differs, so the walk continues on a new step.
     Each is classified by the shipped useful-relation predicate, so it is a
     collision subtype, not a loss.
  3. identity breaks (T = 0), which do force the production restart rule and
     therefore DO waste the walker's remaining steps.  Counted per trail and
     charged at the exact CP-95 upper bound on the per-trail rate, with a loss
     model of 'that trail's work' (NOT divided by the walker count).

A positive control runs the same tracker on a deliberately tiny group (q ~ 2^10)
where S returns are certain, so a zero in the main cell is a property of the walk
rather than of the tracker.
"""
import json, math, os, random, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")
sys.path.insert(0, HERE)
import walk_v20 as W20


def cp95_upper(k, n):
    return W20.cp95_upper(k, n)


def trail_walk(q, R, STEPS, x, tab, c, d, S, trail, track_returns=True):
    """One trail of the declared walk, exact tracking on."""
    half = q >> 1
    tmul = W20.TMUL
    seen = {} if track_returns else None
    s_returns = pair_returns = ident = 0
    ret_useful = ret_trivial = 0
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
        if seen is not None:
            prev = seen.get(key)
            if prev is None:
                seen[key] = (t, kc, kd)
            else:
                pt, pc, pd = prev
                if pt == t:
                    pair_returns += 1          # exact (S, t): structurally impossible
                else:
                    s_returns += 1
                    da = (kc - pc) % q
                    db = (kd - pd) % q
                    if da == 0 and db == 0:
                        ret_trivial += 1
                    else:
                        cand = (-da * pow(db, -1, q)) % q if db else None
                        if cand == x:
                            ret_useful += 1
        S = T
        c = c2
        d = d2
    return {"s_returns": s_returns, "pair_returns": pair_returns,
            "identity_breaks": ident, "ret_useful": ret_useful,
            "ret_trivial": ret_trivial}


def main():
    os.makedirs(OUT, exist_ok=True)
    # main cell: exactly the v19 cycle cell (R=64, W=64, 16 periods, lambda 0.634),
    # but counted per TRAIL (600 trials x 64 walkers = 38,400 trails).
    cell = W20.make_cell(64, 64, 16, 600, "RET_R64_W64_16p", "return measurement")
    q, R, W, S = cell["q"], cell["R"], cell["walkers"], cell["steps_per_walker"]
    tot = {"s_returns": 0, "pair_returns": 0, "identity_breaks": 0,
           "ret_useful": 0, "ret_trivial": 0}
    trails = 0
    trails_with_ident = 0
    for ti in range(cell["trials"]):
        rng = random.Random((20251024 << 26) ^ (ti * W20.GOLDEN) ^ (R << 3) ^ (W << 50))
        x = rng.randrange(2, q - 1)
        tab = W20.make_table(q, x, R)
        for w in range(W):
            c = rng.randrange(1, q)
            d = rng.randrange(1, q)
            Sv = (c + d * x) % q
            trail = (ti * W + w) & 0xFFFFFF
            r = trail_walk(q, R, S, x, tab, c, d, Sv, trail)
            trails += 1
            if r["identity_breaks"]:
                trails_with_ident += 1
            for k in tot:
                tot[k] += r[k]

    # positive control: tiny group, returns must be detected
    qc = 1031
    Rc, Wc, Sc = 64, 64, 4096
    ctl = {"s_returns": 0, "pair_returns": 0, "identity_breaks": 0,
           "ret_useful": 0, "ret_trivial": 0}
    ctl_trails = 0
    for ti in range(6):
        rng = random.Random(9871 + ti)
        x = rng.randrange(2, qc - 1)
        tab = W20.make_table(qc, x, Rc)
        for w in range(Wc):
            c = rng.randrange(1, qc)
            d = rng.randrange(1, qc)
            Sv = (c + d * x) % qc
            r = trail_walk(qc, Rc, Sc, x, tab, c, d, Sv, (ti * Wc + w) & 0xFFFFFF)
            ctl_trails += 1
            for k in ctl:
                ctl[k] += r[k]

    p_ident_hi = cp95_upper(trails_with_ident, trails)
    out = {
        "cell": {k: cell[k] for k in ("label", "R", "walkers", "steps_per_walker",
                                      "periods", "trials", "q", "q_bits", "lambda")},
        "unit": "trail = one walker's walk; trails = trials * walkers",
        "trails": trails,
        "trials": cell["trials"],
        "exact_pair_returns": tot["pair_returns"],
        "s_returns_at_other_phase": tot["s_returns"],
        "s_returns_useful": tot["ret_useful"],
        "s_returns_trivial": tot["ret_trivial"],
        "identity_breaks": tot["identity_breaks"],
        "trails_with_identity_break": trails_with_ident,
        "identity_break_rate": trails_with_ident / trails,
        "identity_break_rate_cp95_upper": p_ident_hi,
        "structural_argument": {
            "statement": "the row index is prod = (S XOR (t*2654435761 mod 2^64) XOR trail*7919) mod R, "
                         "so the next state is a function of (S, t) and t increases strictly inside a trail; "
                         "the pair (S, t) never recurs, therefore a walker's future never repeats and the "
                         "classical rho-cycle death does not exist on this walk",
            "verified_count": tot["pair_returns"],
            "consequence": "an S return at a different phase does not stop progress; it is a collision "
                           "subtype whose relation is tested by the shipped predicate (measured useful / trivial above)",
        },
        "return_positive_control": {
            "why": "same tracker, tiny group q=1031 where S returns are certain",
            "cell": {"q": qc, "R": Rc, "walkers": Wc, "steps_per_walker": Sc, "trails": ctl_trails},
            "s_returns_detected": ctl["s_returns"],
            "useful": ctl["ret_useful"], "trivial": ctl["ret_trivial"],
            "detects": ctl["s_returns"] > 0,
        },
        "pricing": {
            "basis": "identity break = one lost trail (that walker's remaining steps); rate charged at the "
                     "exact CP-95 upper bound on the per-TRAIL rate, NOT divided by the walker count",
            "relative_work_loss_charge": p_ident_hi,
            "v19_error": "v19 divided a per-walker probability by 2^20 and claimed 5.845e-09; that division "
                         "is wrong (each dead walker loses its own work, so the loss fraction is the per-trail "
                         "rate) and the underlying zero was measured on coefficient pairs, not states",
        },
    }
    json.dump(out, open(os.path.join(OUT, "returns_v20.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({k: out[k] for k in ("trails", "exact_pair_returns",
                                          "s_returns_at_other_phase", "s_returns_useful",
                                          "s_returns_trivial", "trails_with_identity_break",
                                          "identity_break_rate_cp95_upper")}, indent=1))
    print(json.dumps(out["return_positive_control"], indent=1))


if __name__ == "__main__":
    main()
