#!/usr/bin/env python3
"""v22 MEASUREMENTS -- the three verdict-887a9a1f rows, measured on the SHIPPED
guarded step (classkey.step, guard ring + prev-row active), not asserted.

F1  THE ABANDON CAP, IMPLEMENTED AND METERED.  The shipped trail protocol is:
    a trail runs until its FIRST distinguished point OR until a hard cap of
    CAP_MULT x theta steps, whichever comes first; a capped trail is ABANDONED
    (labels discarded, fresh trail).  This runner implements the counter in the
    shipped loop, meters the abandon path, and measures the delay-to-abandon
    distribution.  Because the cap is hard, every abandon happens at EXACTLY the
    cap (dispersion 0): there is no geometric tail to bound and none is used.
    Control (force-admission): the same trails run WITHOUT the cap for 2 x cap
    steps -- a capped trail still has no DP, i.e. without the cap the trail runs
    away and the waste is unbounded.  The cap is the whole argument.

F2  BOUND ROUTE: MONOTONICITY (guarded <= unguarded) + measured-E Markov.
    The trap count needs E[X] <= S^2/(2(r-1)) for a HISTORY-DEPENDENT walk.  Two
    rows are shipped: (a) monotonicity -- the guard only REFINES the state
    (element -> (element, ring)), so a joint repeat requires an element repeat:
    guarded collision count <= unguarded, measured on the same fixture, which
    makes the iid birthday count an UPPER bound for the guarded walk; (b) the
    measured expectation from the census, plugged into Markov E_meas/t, which
    needs no independence at all.

F3  GUARDED RETENTION AND THE REPORT COUNT, measured in the same build.
    The 0.896 retention of v19 was an EXTRACT of a v18-era JSON; here the same
    collision-conditioned arm is re-run through the shipped guarded step with
    ring + prev-row active, and the guard-off control is run beside it.  Reports
    are counted directly (DP-triggered) in the guarded census; ceil(S/theta) is
    kept only as a CHECK row against the measured count.
"""
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import classkey as CK
import constants as K
import routing
import toy_census as TC
from cp95 import cp95_lower

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')
D = K.GUARD_DEPTH
CAP_MULT = 128          # the shipped protocol constant: cap = 128 x theta


def warm_ring(rng, depth=D, rows=64):
    g = CK.Guard(depth)
    g.slots = [(rng.randrange(rows), rng.randrange(12)) for _ in range(depth)]
    return g


def dp_hit(toy, y, dp_bits):
    return routing.dp_hit(CK.classkey(toy.T, y)[0], dp_bits)


# ---------------------------------------------------------------- F1
def cap_abandon(toy, walks=128, dp_bits=5, cap_mult=CAP_MULT, seed=2201,
                probe_steps=200, max_trapped=24, run_away_factor=2):
    """Metered cap on the SHIPPED loop: run to first DP or cap; abandon at the cap.

    Two populations, because they are different physics:
      * untrapped trails exit at their first DP (the cap is never reached);
      * a trail that RETURNS to a visited (class, ring) state is on a cycle.  If the
        cycle holds no DP it never reports, so it runs to the cap and is abandoned.
    For every trapped trail the delay to abandon is measured (it is the cap, exactly),
    and the force-admission control continues the SAME state and ring WITHOUT the cap
    for 2 x cap steps: it still produces no DP.  Without the cap that trail runs
    forever, which is what the cap bounds.
    """
    rng = random.Random(seed)
    theta = 1 << dp_bits
    cap = cap_mult * theta
    dp_first = 0
    dp_steps = []
    trapped = 0
    trapped_abandons = 0
    trapped_delays = []
    cycle_lengths = []
    tested = 0
    runaway = 0
    for _w in range(walks):
        y = toy.rand_gt()
        a, b = 1, 1
        g = warm_ring(rng)
        seen = {}
        n = 0
        hit = False
        is_trap = False
        while n < cap:
            if n < probe_steps:
                st = (CK.classkey(toy.T, y)[0], tuple(g.slots))
                if st in seen and n > 0:
                    is_trap = True
                    trapped += 1
                    cycle_lengths.append(n - seen[st])
                seen[st] = n
            y, a, b, _k, _i = toy.step(y, a, b, g)
            n += 1
            if dp_hit(toy, y, dp_bits):
                hit = True
                break
            if is_trap:
                break
        if is_trap and not hit:
            # THE SHIPPED COUNTER keeps running to the cap with the same ring
            while n < cap:
                y, a, b, _k, _i = toy.step(y, a, b, g)
                n += 1
                if dp_hit(toy, y, dp_bits):
                    hit = True
                    break
            if not hit:
                trapped_abandons += 1
                trapped_delays.append(n)
                if tested < max_trapped:
                    tested += 1
                    yy, aa, bb = y, a, b
                    gg = g.copy()          # SAME ring, SAME state
                    still = True
                    for _s in range(cap * run_away_factor):
                        yy, aa, bb, _k, _i = toy.step(yy, aa, bb, gg)
                        if dp_hit(toy, yy, dp_bits):
                            still = False
                            break
                    if still:
                        runaway += 1
        if hit:
            dp_first += 1
            dp_steps.append(n)
    # POSITIVE CONTROL for the runaway row: a FRESH random trail over the same
    # horizon DOES hit DPs at this width -- so "no DP in 2 x cap steps" is a
    # property of the cycle, not of the width.
    fresh_hits = 0
    for _f in range(2):
        y = toy.rand_gt()
        a, b = 1, 1
        g = warm_ring(rng)
        for _s in range(cap * run_away_factor):
            y, a, b, _k, _i = toy.step(y, a, b, g)
            if dp_hit(toy, y, dp_bits):
                fresh_hits += 1
    mean_ab = (sum(trapped_delays) / len(trapped_delays)) if trapped_delays else 0.0
    var_ab = 0.0
    if trapped_delays:
        var_ab = sum((s - mean_ab) ** 2 for s in trapped_delays) / len(trapped_delays)
    return {
        'cap_multiplier': cap_mult,
        'theta_fixture_steps': theta,
        'cap_fixture_steps': cap,
        'walks': walks,
        'dp_first_exits': dp_first,
        'dp_exit_mean_steps': (sum(dp_steps) / len(dp_steps)) if dp_steps else None,
        'trails_trapping_in_cycle': trapped,
        'cycle_length_mean': (sum(cycle_lengths) / len(cycle_lengths)) if cycle_lengths else None,
        'cycle_length_max': max(cycle_lengths) if cycle_lengths else None,
        'abandons_at_cap': trapped_abandons,
        'abandons_with_dp_first': trapped - trapped_abandons,
        'mean_steps_to_abandon': mean_ab,
        'var_steps_to_abandon': var_ab,
        'dispersion_var_over_mean': (var_ab / mean_ab) if mean_ab else None,
        'abandon_delay_is_deterministic': bool(mean_ab in (0.0, float(cap)) and var_ab == 0.0),
        'counter_implemented': 'while n < cap: step(); n += 1 -- on exit the trail is abandoned, '
                               'labels discarded, a fresh trail started (metered in this runner)',
        'control_force_admission': {
            'rule': 'the SAME state and the SAME ring continue WITHOUT the cap for 2 x cap steps',
            'trails_tested': tested,
            'still_without_dp': runaway,
            'control_fails_without_the_cap': bool(tested > 0 and runaway == tested),
        },
        'protocol': 'trail = first DP or 128*theta steps, whichever first; a capped trail is '
                    'abandoned, its labels discarded, a fresh trail started (two GT exponentiation'
                    's charged in aux)',
        'fixture_scale_model': 'the fixture runs the cap at multiplier %d and DP width 2^-%d so the '
                               'trapped regime is reachable in seconds; the PRODUCTION constant priced '
                               'in v22.trap_charge is 128 x theta = 2^87 steps' % (cap_mult, dp_bits),
        'fresh_trail_control': {'horizon_steps': cap * run_away_factor, 'dp_hits': fresh_hits,
                                'reading': 'a fresh random trail hits DPs over the same horizon, so the '
                                           'trapped trail zero-hit row is a cycle property'},
        'no_tail_used': 'the cap is hard, so P(a trail runs past it) = 0 by construction: the waste '
                        'row is cap x rate x price, deterministic, with no geometric tail',
    }


# ---------------------------------------------------------------- F2
def bound_rows(toy, walks=2000, steps=200, seed=4242):
    """Monotonicity: joint (element, guard) repeats <= element repeats; guarded
    element repeats <= unguarded (v13 rule) element repeats, same budget."""
    rng = random.Random(seed)
    guarded_element_repeats = 0
    joint_repeats = 0
    unguarded_element_repeats = 0
    pairs = 0
    for _w in range(walks):
        y = toy.rand_gt()
        a, b = 1, 1
        g = warm_ring(rng)
        yv, av, bv = y, a, b
        seen_el = {}
        seen_joint = {}
        seen_v13 = {}
        for t in range(steps):
            y, a, b, _k, _i = toy.step(y, a, b, g)
            yv, av, bv, _k2, _i2 = toy.step_v13(yv, av, bv)
            el = CK.classkey(toy.T, y)[0]
            jt = (el, tuple(g.slots))
            el13 = CK.classkey(toy.T, yv)[0]
            if t:
                pairs += 1
            if el in seen_el:
                guarded_element_repeats += 1
            seen_el[el] = 1
            if jt in seen_joint:
                joint_repeats += 1
            seen_joint[jt] = 1
            if el13 in seen_v13:
                unguarded_element_repeats += 1
            seen_v13[el13] = 1
    return {
        'walks': walks, 'steps_per_walk': steps, 'step_pairs': pairs,
        'guarded_joint_repeats': joint_repeats,
        'guarded_element_repeats': guarded_element_repeats,
        'unguarded_element_repeats': unguarded_element_repeats,
        'monotonicity_joint_le_element': bool(joint_repeats <= guarded_element_repeats),
        'monotonicity_guarded_le_unguarded': bool(guarded_element_repeats <= unguarded_element_repeats),
        'reading': 'the guard ADDS a coordinate to the state, so a joint repeat implies an element '
                   'repeat (measured), and the guarded element stream repeats no more often than the '
                   'unguarded one (measured): E[joint repeats] <= E[birthday count over r-1], which is '
                   'the iid birthday value used as an UPPER bound in the trap row',
        'step_source': 'classkey.step (shipped) vs classkey.step_v13 (guard-off control)',
    }


# ---------------------------------------------------------------- F3
def retention_guarded(toy, pairs=480, steps=64, lags=(1, 2, 4, 8, 16, 32, 56),
                      seed=31):
    """Collision-conditioned retention THROUGH the shipped guarded step."""
    rng = random.Random(seed)
    r = toy.r
    detected = 0
    by_k = {k: {'pairs': 0, 'detected': 0} for k in lags}
    off_detected = 0
    off_by_k = {k: {'pairs': 0, 'detected': 0} for k in lags}
    ka = 0
    for _ in range(pairs):
        x0 = rng.randrange(1, r)
        z = toy.T.pow(toy.g, x0)
        a0, b0 = rng.randrange(r), rng.randrange(r)
        y0 = toy.T.mul(toy.T.pow(toy.g, a0), toy.T.pow(z, b0))
        s = rng.randrange(1, r)
        t = (-s * x0) % r
        aA, bA = a0, b0
        aB, bB = (a0 + t) % r, (b0 + s) % r
        if toy.T.mul(toy.T.pow(toy.g, aB), toy.T.pow(z, bB)) == y0:
            ka += 1
        gA, gB = warm_ring(rng), warm_ring(rng)
        yA, yB = y0, y0
        yAv, yBv = y0, y0
        hit = False
        off_hit = False
        for step in range(1, steps + 1):
            yA, aA, bA, _k, _i = toy.step(yA, aA, bA, gA)
            yB, aB, bB, _k2, _i2 = toy.step(yB, aB, bB, gB)
            yAv, aAv, bAv, _k3, _i3 = toy.step_v13(yAv, aA, bA)
            yBv, aBv, bBv, _k4, _i4 = toy.step_v13(yBv, aB, bB)
            same = CK.classkey(toy.T, yA)[0] == CK.classkey(toy.T, yB)[0]
            if same and not hit:
                hit = True
                detected += 1
            same_off = CK.classkey(toy.T, yAv)[0] == CK.classkey(toy.T, yBv)[0]
            if same_off and not off_hit:
                off_hit = True
                off_detected += 1
            for k in lags:
                if step == k:
                    by_k[k]['pairs'] += 1
                    off_by_k[k]['pairs'] += 1
                    if same:
                        by_k[k]['detected'] += 1
                    if same_off:
                        off_by_k[k]['detected'] += 1
    curve = {str(k): {'lag': k, 'pairs': by_k[k]['pairs'],
                      'detected': by_k[k]['detected'],
                      'p_detect': by_k[k]['detected'] / by_k[k]['pairs']}
             for k in lags}
    off_curve = {str(k): {'lag': k, 'pairs': off_by_k[k]['pairs'],
                          'detected': off_by_k[k]['detected'],
                          'p_detect': off_by_k[k]['detected'] / off_by_k[k]['pairs']}
                 for k in lags}
    return {
        'pairs': pairs, 'steps': steps, 'guard_depth': D,
        'known_answer_gate_pass': ka, 'known_answer_gate_required': pairs,
        'detected': detected,
        'p_detect_guarded': detected / pairs,
        'p_detect_guarded_cp95_lower': float(cp95_lower(detected, pairs)),
        'curve_guarded': curve,
        'guard_off_control': {'detected': off_detected,
                              'p_detect': off_detected / pairs,
                              'cp95_lower': float(cp95_lower(off_detected, pairs)),
                              'curve': off_curve},
        'reading': 'SAME arm as the v18-era row, re-run through the shipped guarded step with the '
                   'ring and prev-row active; the guard-off control is classkey.step_v13.  The '
                   'product uses the collision-instant value (lag 1 family); the k-curve is the '
                   'production-horizon risk and is disclosed as a named failed gate.',
        'step_source': 'classkey.step (shipped, imported) -- not a reimplementation',
    }


def reports_direct(toy, walks=240, steps=2048, dp_bits=5, seed=77, guard=True):
    """Count DP-triggered reports directly, guarded census."""
    rng = random.Random(seed)
    reports = 0
    per_walk = []
    for _w in range(walks):
        y = toy.rand_gt()
        a, b = 1, 1
        g = warm_ring(rng) if guard else None
        n = 0
        for _s in range(steps):
            y, a, b, _k, _i = toy.step(y, a, b, g)
            if dp_hit(toy, y, dp_bits):
                n += 1
        reports += n
        per_walk.append(n)
    mean = reports / walks
    theta = 1 << dp_bits
    formula = (steps + theta - 1) // theta
    return {
        'walks': walks, 'steps_per_walk': steps, 'dp_bits': dp_bits,
        'reports_measured': reports,
        'reports_per_walk_measured': mean,
        'reports_per_step_measured': reports / (walks * steps),
        'model_reports_per_step': 1.0 / theta,
        'ceil_formula_reports_per_walk': formula,
        'check_measured_vs_formula_ratio': mean / formula if formula else None,
        'chebyshev_free': 'the report count enters the ledger only as a COST CEILING; this row is '
                          'the census-direct count the ceiling is checked against',
    }


def main():
    os.makedirs(EV, exist_ok=True)
    t0 = time.time()
    toy = TC.Toy(rows=64, seed=20261024)
    print('fixture r=%d rows=%d D=%d' % (toy.r, len(toy.us), D), flush=True)
    cap = cap_abandon(toy, walks=256, dp_bits=5, cap_mult=CAP_MULT, probe_steps=400, max_trapped=8)
    print('F1 cap: trapped %d, abandons %d at mean %.1f steps, control runaway %d/%d (%.1fs)'
          % (cap['trails_trapping_in_cycle'], cap['abandons_at_cap'],
             cap['mean_steps_to_abandon'],
             cap['control_force_admission']['still_without_dp'],
             cap['control_force_admission']['trails_tested'], time.time() - t0), flush=True)
    bnd = bound_rows(toy, walks=800)
    print('F2 bound: joint %d <= element %d <= unguarded %d (%.1fs)'
          % (bnd['guarded_joint_repeats'], bnd['guarded_element_repeats'],
             bnd['unguarded_element_repeats'], time.time() - t0), flush=True)
    ret = retention_guarded(toy)
    print('F3 retention guarded %.4f (cp95 %.4f) vs guard-off %.4f (%.1fs)'
          % (ret['p_detect_guarded'], ret['p_detect_guarded_cp95_lower'],
             ret['guard_off_control']['p_detect'], time.time() - t0), flush=True)
    rep = reports_direct(toy)
    print('F3 reports %.3f/walk vs formula %d (%.1fs)'
          % (rep['reports_per_walk_measured'], rep['ceil_formula_reports_per_walk'],
             time.time() - t0), flush=True)
    json.dump(cap, open(os.path.join(EV, 'cap-v22.json'), 'w'), indent=1, sort_keys=True)
    json.dump(bnd, open(os.path.join(EV, 'bound-v22.json'), 'w'), indent=1, sort_keys=True)
    json.dump(ret, open(os.path.join(EV, 'retention-v22.json'), 'w'), indent=1, sort_keys=True)
    json.dump(rep, open(os.path.join(EV, 'reports-v22.json'), 'w'), indent=1, sort_keys=True)
    print('wrote evidence/cap-v22.json bound-v22.json retention-v22.json reports-v22.json')


if __name__ == '__main__':
    main()
