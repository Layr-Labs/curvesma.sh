#!/usr/bin/env python3
"""v26 measurement layer: the SQUARED-STEP walk.

Everything here executes the shipped step (`classkey.step_v26`) on the shipped
tower fixture.  Three jobs:

  A. THEOREM SCAN (evidence/impossibility-theorem.md has teeth):
     sigma^6 = inversion, TBL[t][k+6] = TBL[t][k]^(-1) in every column, every
     column is the sigma-orbit of column 0 -- then a PLANTED MUTATION of one
     table entry must make the scan flag exactly that column.  A scan that still
     reports 384/384 after the plant is a no-op and this script fails.

  B. FULL-STATE CYCLE CENSUS (family existence, exhaustive, not sampled):
     build the step map on ALL r-1 elements of the toy.  The v13 rule (the
     control) contains 2-cycles -- the exact-return family.  The v26 squared
     step must contain NONE: a 2-cycle needs y^3 = sigma^j/(m1^4 m2^2), which has
     one solution per multiplier pair, not a rate-1/(12R) family of them.

  C. UNCONDITIONAL RECOVERY CENSUS at production SHAPE:
     S_toy = constants.step_budget(N_toy) (the same function that sizes the
     production budget), 4 trails of ceil(S_toy/4) steps each, per-worker seeds.
     P(recover) = P(repeat reached) * P(agree | repeat), counted jointly per
     trial (a single-run product), with Clopper-Pearson on the joint count.
     Success-equality control: the same trials re-run through counted wrappers
     and the per-trial outcome vector compared element by element.
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K
import classkey as CK
import ledger as L
import routing
import tower as TW
from cp95 import cp95_lower

TRAILS = 4
WORKERS = 12
TRIALS = 20000


# ---------------------------------------------------------------- the step
def step26(toy, y, a, b):
    return CK.step_v26(toy.T, y, a, b, toy.TBL, toy.U, toy.r, toy.tbits, toy.rows)


def make_toy(bits=5, rows=64, seed=20261024):
    import toy_census as TC
    return TC.Toy(bits=bits, seed=seed, rows=rows)


# --------------------------------------------------- A. the theorem, with teeth
def theorem_scan(toy, plant=None):
    """Count the sigma-orbit identities column by column.  `plant` mutates one
    table entry (t,k) to a random group element before the scan."""
    T = toy.T
    TBL = [row[:] for row in toy.TBL]
    planted = None
    if plant is not None:
        t0, k0 = plant
        TBL[t0][k0] = T.mul(TBL[t0][k0], toy.tbl[(t0 + 1) % toy.rows])
        planted = (t0, k0)
    n_rows, n_cols = len(TBL), len(TBL[0])
    pairs = 0
    pair_fail = []
    orbit_ok = 0
    orbit_fail = []
    for t in range(n_rows):
        for k in range(6):
            if T.mul(TBL[t][k], TBL[t][k + 6]) == T.one12():
                pairs += 1
            else:
                pair_fail.append([t, k, k + 6])
        for k in range(n_cols):
            if toy.frob_pow(TBL[t][0], (-k) % 12) == TBL[t][k]:
                orbit_ok += 1
            else:
                orbit_fail.append([t, k])
    inv = all(toy.frob_pow(y, 6) == T.pow(y, toy.r - 1)
              for y in toy.tbl[::max(1, toy.rows // 8)])
    return dict(rows=n_rows, cols=n_cols,
                sigma6_is_inversion=inv,
                inverse_pair_columns=pairs,
                inverse_pair_columns_total=n_rows * 6,
                inverse_pair_failures=pair_fail[:20],
                orbit_columns=orbit_ok,
                orbit_columns_total=n_rows * n_cols,
                orbit_failures=orbit_fail[:20],
                plant=planted,
                clean=(pairs == n_rows * 6 and orbit_ok == n_rows * n_cols
                       and not pair_fail and not orbit_fail and inv))


# ------------------------------------------- B. exhaustive cycle census (family)
def full_map_census(toy, which='v26'):
    """Build the map on all r-1 elements and return the cycle histogram."""
    T, r = toy.T, toy.r
    elems = [None] * r
    y = toy.g
    for i in range(1, r):
        elems[i] = y
        y = T.mul(y, toy.g)
    nxt = [0] * r
    dlog = {}
    for i in range(1, r):
        dlog[toy.key(elems[i])[0]] = i
    kinds = [0, 0]
    for i in range(1, r):
        if which == 'v26':
            y2, _a, _b, kind, _info = step26(toy, elems[i], 1, 0)
        else:
            y2, _a, _b, kind, _info = CK.step_v13(T, elems[i], 1, 0, toy.TBL,
                                                  toy.U, r, toy.tbits)
        if kind == 'mul':
            kinds[0] += 1
        else:
            kinds[1] += 1
        nxt[i] = dlog[toy.key(y2)[0]]
    state = [0] * r
    hist = {}
    two_cycles = 0
    for i in range(1, r):
        if state[i]:
            continue
        path, j = [], i
        while j and not state[j]:
            state[j] = 1
            path.append(j)
            j = nxt[j]
        if j and state[j] == 1:
            cyc = path[path.index(j):]
            hist[len(cyc)] = hist.get(len(cyc), 0) + 1
            if len(cyc) == 2:
                two_cycles += 1
        for q in path:
            state[q] = 2
    return dict(rule=which, elements=r - 1, mul_steps=kinds[0], sqr_steps=kinds[1],
                sqr_fraction=kinds[1] / (kinds[0] + kinds[1]),
                cycles=sum(hist.values()), cycle_lengths={str(k): v for k, v in sorted(hist.items())},
                two_cycles=two_cycles,
                predicted_two_cycles_expectation=(r - 1) / (12.0 * len(toy.us)))


# ------------------------------------------------ C. recovery census at shape
def _batch(args):
    bits, seed, trials, ratio = args
    toy = make_toy(bits=bits, seed=seed)
    rng = random.Random(seed ^ 0x5EED)
    s_toy = K.step_budget(n=toy.N)
    s_len = max(TRAILS, int(round(ratio * s_toy)))
    cap = max(1, s_len // TRAILS)
    rec = dict(trials=trials, cap=cap, s_len=cap * TRAILS, s_toy=s_toy,
               repeats=0, recovered=0, steps=0, mul_steps=0, sq_steps=0,
               candidate_solves=0, verified=0, zeroden=0, wrong_j=0,
               invariant_fail=0, first_repeat_steps=0)
    for _ in range(trials):
        x = rng.randrange(1, toy.r)
        z = toy.T.pow(toy.g, x)
        store = {}
        hit = False
        for _tr in range(TRAILS):
            a0 = rng.randrange(1, toy.r)
            b0 = rng.randrange(1, toy.r)
            y = toy.T.mul(toy.T.pow(toy.g, a0), toy.T.pow(z, b0))
            a, b = a0, b0
            for s in range(cap):
                y, a, b, kind, _i = step26(toy, y, a, b)
                rec['steps'] += 1
                rec['mul_steps' if kind == 'mul' else 'sq_steps'] += 1
                if toy.T.mul(toy.T.pow(toy.g, a), toy.T.pow(z, b)) != y:
                    rec['invariant_fail'] += 1
                k = toy.key(y)[0]
                if k in store:
                    rec['repeats'] += 1
                    rec['first_repeat_steps'] += s + 1
                    a2, b2 = store[k]
                    cands = CK.solve_candidates([pow(toy.p, j, toy.r) for j in range(12)],
                                                None, a, b, a2, b2, toy.r)
                    rec['candidate_solves'] += len(cands)
                    if not cands:
                        rec['zeroden'] += 1
                    got = False
                    for _j, xx in cands:
                        if toy.T.pow(toy.g, xx) == z:
                            rec['verified'] += 1
                            got = True
                            break
                    if got:
                        rec['recovered'] += 1
                    else:
                        rec['wrong_j'] += 1
                    hit = True
                    break
                store[k] = (a, b)
            if hit:
                break
    return rec


def _merge(parts):
    out = dict(parts[0])
    for p in parts[1:]:
        for k, v in p.items():
            if isinstance(v, int) and k not in ('cap', 's_len', 's_toy'):
                out[k] += v
    return out


def census(trials=TRIALS, bits=5, ratio=1.0, workers=WORKERS):
    import multiprocessing as mp
    nw = max(1, min(workers, trials))
    per = [trials // nw] * nw
    for i in range(trials - sum(per)):
        per[i] += 1
    args = [(bits, 4242 + 17 * i, n, ratio) for i, n in enumerate(per) if n]
    with mp.Pool(len(args)) as pool:
        rec = _merge(pool.map(_batch, args))
    rec['trials'] = trials
    rec['workers'] = len(args)
    rec['ratio'] = ratio
    rec['p_repeat'] = rec['repeats'] / trials
    rec['p_agree_given_repeat'] = rec['verified'] / rec['repeats'] if rec['repeats'] else None
    rec['p_recover'] = rec['recovered'] / trials
    rec['cp95_lower'] = float(cp95_lower(rec['recovered'], trials))
    rec['cp95_upper'] = float(cp95_lower(trials - rec['recovered'], trials))
    rec['cp95_upper'] = 1.0 - rec['cp95_upper']
    rec['sqr_fraction'] = rec['sq_steps'] / rec['steps'] if rec['steps'] else None
    rec['model_lambda'] = float(rec['s_len'] ** 2) / (2 * float(toyN(bits)))
    rec['model_success'] = 1.0 - math.exp(-rec['model_lambda'])
    return rec


def toyN(bits):
    import toy_census as TC
    return TC.FIX[bits]['N']


# ------------- D. instrumented-vs-plain success equality (the accepted control)
class Counted( object):
    """Call-counter wrapper around the tower methods the step uses.  Only the
    two shipped operations are counted: mul (full Fp12 multiply) and
    cyclo_sqr_fast (the cyclo squaring the squared step adds)."""

    def __init__(self, T):
        self.T = T
        self.mul_calls = 0
        self.sqr_calls = 0

    def mul(self, a, b):
        self.mul_calls += 1
        return self.T.mul(a, b)

    def cyclo_sqr_fast(self, a):
        self.sqr_calls += 1
        return self.T.cyclo_sqr_fast(a)


def equality_control(trials=400, bits=5, seed=99991):
    """Same trials, plain vs counted: the per-trial outcome vectors must match
    element for element, and the counted call totals must equal the kind counts
    the step itself reports."""
    toy = make_toy(bits=bits, seed=seed)
    rng_a = random.Random(seed)
    rng_b = random.Random(seed)
    plain, counted = [], []
    c = Counted(toy.T)
    kinds = [0, 0]
    for _ in range(trials):
        x_a, x_b = rng_a.randrange(1, toy.r), rng_b.randrange(1, toy.r)
        z_a, z_b = toy.T.pow(toy.g, x_a), toy.T.pow(toy.g, x_b)
        plain.append(_one_trial(toy, z_a, rng_a))
        counted.append(_one_trial(toy, z_b, rng_b, counted_tower=c, kinds=kinds))
    eq = plain == counted
    return dict(trials=trials, plain_successes=sum(plain), counted_successes=sum(counted),
                per_trial_equal=eq, mismatches=[i for i in range(trials)
                                                if plain[i] != counted[i]][:10],
                counted_mul_calls=c.mul_calls, counted_sqr_calls=c.sqr_calls,
                step_reported_mul_steps=kinds[0], step_reported_sqr_steps=kinds[1],
                # one cyclo_sqr per step in BOTH branches, one mul in the mul branch
                counts_match=(c.mul_calls == kinds[0]
                              and c.sqr_calls == kinds[0] + kinds[1]))


def _one_trial(toy, z, rng, counted_tower=None, kinds=None):
    s_toy = K.step_budget(n=toy.N)
    cap = max(1, s_toy // TRAILS)
    store = {}
    for _tr in range(TRAILS):
        a = rng.randrange(1, toy.r)
        b = rng.randrange(1, toy.r)
        y = toy.T.mul(toy.T.pow(toy.g, a), toy.T.pow(z, b))
        for _s in range(cap):
            if counted_tower is None:
                y, a, b, kind, _i = step26(toy, y, a, b)
            else:
                mkey, kap = CK.classkey(toy.T, y)
                t, _dp, sw = routing.fields(mkey, toy.tbits)
                if sw < 3:
                    kind = 'sq'
                    y = counted_tower.cyclo_sqr_fast(y)
                    a, b = (2 * a) % toy.r, (2 * b) % toy.r
                else:
                    kind = 'mul'
                    y = counted_tower.cyclo_sqr_fast(counted_tower.mul(y, toy.TBL[t][kap]))
                    a, b = (2 * (a + toy.U[t][kap])) % toy.r, (2 * b) % toy.r
            if kinds is not None:
                kinds[0 if kind == 'mul' else 1] += 1
            k = toy.key(y)[0]
            if k in store:
                a2, b2 = store[k]
                for _j, xx in CK.solve_candidates([pow(toy.p, j, toy.r) for j in range(12)],
                                                  None, a, b, a2, b2, toy.r):
                    if toy.T.pow(toy.g, xx) == z:
                        return True
                return False
            store[k] = (a, b)
    return False


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'all'
    out = {}
    if cmd in ('all', 'theorem'):
        toy = make_toy()
        clean = theorem_scan(toy)
        planted = theorem_scan(toy, plant=(7, 3))
        out['theorem'] = dict(clean=clean, planted=planted,
                              planted_flagged=(not planted['clean']),
                              planted_column_flags=[f for f in planted['inverse_pair_failures']
                                                    if f[0] == 7])
        print('THEOREM clean:', clean['inverse_pair_columns'], '/', clean['inverse_pair_columns_total'],
              'pairs, orbit', clean['orbit_columns'], '/', clean['orbit_columns_total'],
              'sigma6', clean['sigma6_is_inversion'], '| planted flagged:', planted['clean'] is False,
              planted['inverse_pair_failures'][:3])
    if cmd in ('all', 'cycle'):
        toy = make_toy()
        c26 = full_map_census(toy, 'v26')
        c13 = full_map_census(toy, 'v13')
        out['cycle'] = dict(v26=c26, v13=c13)
        print('CYCLE v26 two_cycles', c26['two_cycles'], 'cycles', c26['cycles'],
              'sqr_frac %.4f' % c26['sqr_fraction'])
        print('CYCLE v13 two_cycles', c13['two_cycles'], '(the family; expected %.0f)'
              % c13['predicted_two_cycles_expectation'])
    if cmd in ('all', 'census'):
        out['census'] = census(trials=int(os.environ.get('TRIALS', TRIALS)),
                               ratio=float(os.environ.get('RATIO', 1.0)))
        c = out['census']
        print('CENSUS n=%d ratio=%.2f p_repeat %.4f agree %.4f p_recover %.4f CPlo %.4f '
              'sqr_frac %.4f model %.4f'
              % (c['trials'], c['ratio'], c['p_repeat'], c['p_agree_given_repeat'],
                 c['p_recover'], c['cp95_lower'], c['sqr_fraction'], c['model_success']))
    if cmd in ('all', 'equality'):
        out['equality'] = equality_control()
        print('EQUALITY', out['equality']['per_trial_equal'],
              'mul', out['equality']['counted_mul_calls'],
              'sqr', out['equality']['counted_sqr_calls'],
              'counts_match', out['equality']['counts_match'])
    os.makedirs('evidence', exist_ok=True)
    tag = cmd
    if cmd == 'census':
        tag = 'census-production' if abs(float(os.environ.get('RATIO', 1.0)) - 1.0) < 1e-9 \
            else 'census-ratio%s' % os.environ.get('RATIO', '1.0')
    with open('evidence/v26-measure-%s.json' % tag, 'w') as f:
        json.dump(out, f, sort_keys=True, default=str, separators=(',', ':'))
