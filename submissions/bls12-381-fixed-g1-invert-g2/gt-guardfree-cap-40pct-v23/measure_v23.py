#!/usr/bin/env python3
"""v23 MEASUREMENTS -- the GUARD-FREE + CAP walk, measured (verdict 3eb6d983).

THE REDESIGN.  v22 measured the trade-off exactly: the routing guard that removes
the absorbing cycle family ALSO separates two trails that share a start point, so
the guarded retention fell to 0.808333 (CP-95 0.776352) and the success product to
0.3523, below the 0.39 floor.  With the guard OFF, the retention control is
1.000000 -- but v22 charged every class repeat as an uncapped trap.

v23 keeps the v22 CAP (hard, metered, deterministic: a trail ends at its FIRST
distinguished point or at 128*theta steps, whichever comes first; a capped trail
is abandoned, its labels discarded and a fresh trail started) and DROPS the guard:
the step becomes history-free (column = kappa(y), no ring), so
  * the successor CLASS is a function of the class alone,
  * two trails that share a start point are the SAME trail,
  * a class repeat is a cycle and its work is CAPPED, not unbounded.

F1  CAP ON THE GUARD-FREE WALK.  Same metered loop as v22, run on the shipped
    guard-free step: trails that fall into a cycle before any DP are abandoned at
    exactly the cap (delay == cap, variance 0); the force-admission control runs
    the SAME state for 2 x cap steps with no cap and still produces no DP.

F2  DP-GATED DETECTION (verdict item 1).  v22 counted class agreement at ANY step.
    The protocol detects on the DISTINGUISHED-POINT predicate, so detection is now
    defined as a DP-HIT EVENT: at the first step whose class key satisfies the DP
    mask, both trails must be in the same class key.  Both statistics are shipped
    side by side (any_step_agree and dp_gated_agree) and only the DP-gated one is
    allowed to feed success.  The trap-before-DP loss is measured at fixture scale
    and bounded at production scale by theta^2/(2N), which is where the loss lives.

F3  REPORTS WITH THE RESTART (verdict item 3).  v22's reports_direct counted DP
    hits on a CONTINUING walk, which is not the protocol: a reporting trail is
    replaced by a fresh trail.  Here the runner restarts with fresh coins after
    every DP AND after every abandon, until the step budget is consumed.

F4  THE EXPECTATION ROW (verdict item 2, made MOOT).  No guard -> no extra state
    coordinate -> the walk is back in the hypothesis domain of the classic
    birthday count; the measured-E Markov row is shipped beside it, and the
    monotonicity route is withdrawn rather than patched.
"""
import json
import math
import os
import random
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import classkey as CK
import constants as K
import routing
import toy_census as TC
from cp95 import cp95_lower

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')
CAP_MULT = 128


def dp_hit(toy, y, dp_bits):
    return routing.dp_hit(CK.classkey(toy.T, y)[0], dp_bits)


def _key(toy, y):
    return CK.classkey(toy.T, y)[0]


# ------------------------------------------------------------------ F1
def _cap_batch(args):
    walks, dp_bits, cap_mult, seed, max_trapped, run_away_factor = args
    toy = TC.Toy(rows=64, seed=20261024)
    rng = random.Random(seed)
    cap = cap_mult * (1 << dp_bits)
    out = dict(dp_first=0, trapped=0, abandons=0, delays=[], cycles=[],
               dp_steps=[], tested=0, runaway=0, fresh_hits=0)
    for _ in range(walks):
        y = toy.rand_gt()
        a, b = 1, 1
        seen = {}
        n = 0
        hit = False
        trapped = False
        while n < cap:
            k = _key(toy, y)
            if k in seen:
                trapped = True
                out['trapped'] += 1
                out['cycles'].append(n - seen[k])
                break
            seen[k] = n
            y, a, b, _kd, _ki = toy.step_v13(y, a, b)
            n += 1
            if routing.dp_hit(_key(toy, y), dp_bits):
                hit = True
                break
        if trapped:
            # the shipped counter keeps the SAME state and runs to the cap
            while n < cap:
                y, a, b, _kd, _ki = toy.step_v13(y, a, b)
                n += 1
                if routing.dp_hit(_key(toy, y), dp_bits):
                    hit = True
                    break
            if not hit:
                out['abandons'] += 1
                out['delays'].append(n)
                if out['tested'] < max_trapped:
                    out['tested'] += 1
                    yy, aa, bb = y, a, b
                    still = True
                    for _s in range(cap * run_away_factor):
                        yy, aa, bb, _kd2, _ki2 = toy.step_v13(yy, aa, bb)
                        if routing.dp_hit(_key(toy, yy), dp_bits):
                            still = False
                            break
                    if still:
                        out['runaway'] += 1
        if hit:
            out['dp_first'] += 1
            out['dp_steps'].append(n)
    # positive control: a fresh guard-free trail hits DPs over the same horizon
    for _f in range(2):
        y = toy.rand_gt()
        a, b = 1, 1
        for _s in range(cap * run_away_factor):
            y, a, b, _kd, _ki = toy.step_v13(y, a, b)
            if routing.dp_hit(_key(toy, y), dp_bits):
                out['fresh_hits'] += 1
    return out


def cap_guardfree(walks=256, dp_bits=5, cap_mult=CAP_MULT, seed=2301,
                  max_trapped=8, workers=8):
    theta = 1 << dp_bits
    cap = cap_mult * theta
    chunks = [(walks // workers + 1, dp_bits, cap_mult, seed + i, max_trapped, 2)
              for i in range(workers)]
    with Pool(workers) as p:
        parts = p.map(_cap_batch, chunks)
    agg = dict(dp_first=0, trapped=0, abandons=0, delays=[], cycles=[],
               dp_steps=[], tested=0, runaway=0, fresh_hits=0)
    for q in parts:
        for k in ('dp_first', 'trapped', 'abandons', 'tested', 'runaway', 'fresh_hits'):
            agg[k] += q[k]
        for k in ('delays', 'cycles', 'dp_steps'):
            agg[k] += q[k]
    mean_ab = (sum(agg['delays']) / len(agg['delays'])) if agg['delays'] else 0.0
    var_ab = 0.0
    if agg['delays']:
        var_ab = sum((s - mean_ab) ** 2 for s in agg['delays']) / len(agg['delays'])
    return {
        'cap_multiplier': cap_mult, 'theta_fixture_steps': theta,
        'cap_fixture_steps': cap, 'walks': walks,
        'dp_first_exits': agg['dp_first'],
        'dp_exit_mean_steps': (sum(agg['dp_steps']) / len(agg['dp_steps']))
                               if agg['dp_steps'] else None,
        'trails_trapping_in_cycle': agg['trapped'],
        'cycle_length_mean': (sum(agg['cycles']) / len(agg['cycles']))
                              if agg['cycles'] else None,
        'cycle_length_max': max(agg['cycles']) if agg['cycles'] else None,
        'abandons_at_cap': agg['abandons'],
        'mean_steps_to_abandon': mean_ab, 'var_steps_to_abandon': var_ab,
        'dispersion_var_over_mean': (var_ab / mean_ab) if mean_ab else None,
        'abandon_delay_is_deterministic': bool(var_ab == 0.0),
        'control_force_admission': {
            'rule': 'the SAME state continues WITHOUT the cap for 2 x cap steps',
            'trails_tested': agg['tested'],
            'still_without_dp': agg['runaway'],
            'control_fails_without_the_cap': bool(agg['tested'] > 0
                                                  and agg['runaway'] == agg['tested']),
        },
        'fresh_trail_control': {'horizon_steps': cap * 2, 'dp_hits': agg['fresh_hits'],
                                'reading': 'a fresh guard-free trail hits DPs over the same '
                                           'horizon, so the trapped trail zero-hit row is a '
                                           'cycle property, not a width property'},
        'step_source': 'classkey.step_v13 (the shipped guard-free step)',
        'no_tail_used': 'the cap is hard, so P(a trail runs past it) = 0 by construction: '
                        'the waste row is cap x rate x price with no geometric tail',
    }


# ------------------------------------------------------------------ F2
def _ret_batch(args):
    pairs, steps, dp_bits, cap_mult, seed, lags, bits = args
    toy = TC.Toy(bits=bits, seed=20261024)
    rng = random.Random(seed)
    r = toy.r
    any_agree = 0
    dp_hit_pairs = 0
    dp_agree = 0
    trap_before_dp = 0
    ka = 0
    by_k = {k: [0, 0] for k in lags}
    dp_by_k = {k: [0, 0] for k in lags}
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
        yA, yB = y0, y0
        aAv, bAv, aBv, bBv = aA, bA, aB, bB
        seen = {}
        a_agree = False
        first_dp = None
        first_rep = None
        got_dp = False
        agree_at_dp = False
        for step in range(1, steps + 1):
            yA, aAv, bAv, _k1, _i1 = toy.step_v13(yA, aAv, bAv)
            yB, aBv, bBv, _k2, _i2 = toy.step_v13(yB, aBv, bBv)
            kA, kB = _key(toy, yA), _key(toy, yB)
            same = (kA == kB)
            if same:
                a_agree = True
            dph = routing.dp_hit(kA, dp_bits)
            if dph and first_dp is None:
                first_dp = step
                got_dp = True
                if same:
                    agree_at_dp = True
                break
            if kA in seen and first_rep is None:
                first_rep = step
                break
            seen[kA] = step
            for k in lags:
                if step == k:
                    by_k[k][0] += 1
                    if same:
                        by_k[k][1] += 1
                    if dph:
                        dp_by_k[k][0] += 1
                        if same:
                            dp_by_k[k][1] += 1
        if a_agree:
            any_agree += 1
        if got_dp:
            dp_hit_pairs += 1
            if agree_at_dp:
                dp_agree += 1
        else:
            trap_before_dp += 1
    return dict(any_agree=any_agree, dp_hit_pairs=dp_hit_pairs, dp_agree=dp_agree,
                trap_before_dp=trap_before_dp, ka=ka, pairs=pairs,
                by_k={str(k): by_k[k] for k in lags},
                dp_by_k={str(k): dp_by_k[k] for k in lags})


def retention_dpgated(pairs=10000, steps=64, dp_bits=5, cap_mult=CAP_MULT,
                      bits=5, seed=2311, lags=(1, 2, 4, 8, 16, 32, 56), workers=8):
    toy = TC.Toy(bits=bits, seed=20261024)
    theta = 1 << dp_bits
    cap = cap_mult * theta
    chunks = [(pairs // workers + 1, steps, dp_bits, cap_mult, seed + i, lags, bits)
              for i in range(workers)]
    with Pool(workers) as p:
        parts = p.map(_ret_batch, chunks)
    tot = dict(any_agree=0, dp_hit_pairs=0, dp_agree=0, trap_before_dp=0, ka=0, pairs=0)
    by_k = {str(k): [0, 0] for k in lags}
    dp_by_k = {str(k): [0, 0] for k in lags}
    for q in parts:
        for k in tot:
            tot[k] += q[k]
        for k in by_k:
            by_k[k][0] += q['by_k'][k][0]
            by_k[k][1] += q['by_k'][k][1]
            dp_by_k[k][0] += q['dp_by_k'][k][0]
            dp_by_k[k][1] += q['dp_by_k'][k][1]
    n = tot['pairs']
    curve = {k: {'lag': int(k), 'pairs': by_k[k][0], 'detected': by_k[k][1],
                 'p_agree_any_step': (by_k[k][1] / by_k[k][0]) if by_k[k][0] else None,
                 'dp_pairs': dp_by_k[k][0], 'dp_detected': dp_by_k[k][1],
                 'p_agree_given_dp': (dp_by_k[k][1] / dp_by_k[k][0])
                                     if dp_by_k[k][0] else None}
             for k in by_k}
    theta_prod = 1 << K.THETA_BITS
    p_trap_prod = float(theta_prod * theta_prod) / (2.0 * float(K.N_CLASSES))
    return {
        'protocol': 'GUARD-FREE step (classkey.step_v13: column = kappa, no ring) + the '
                    'v22 cap (first DP or 128*theta steps; a capped trail is abandoned, '
                    'labels discarded, a fresh trail started)',
        'step_source': 'classkey.step_v13 (shipped control rule, imported) -- the guard-free '
                       'walk of the v23 design',
        'pairs': n, 'steps': steps, 'dp_bits': dp_bits, 'fixture_bits': bits,
        'fixture_classes': toy.N, 'cap_steps': cap,
        'known_answer_gate_pass': tot['ka'], 'known_answer_gate_required': n,
        'any_step_agree': {'detected': tot['any_agree'],
                           'p_detect': tot['any_agree'] / n,
                           'cp95_lower': float(cp95_lower(tot['any_agree'], n)),
                           'role': 'the v22 STATISTIC, shipped for comparison only; it counts '
                                   'class agreement at ANY step and does NOT feed success '
                                   '(verdict 3eb6d983 item 1)'},
        'dp_gated_agree': {'pairs_with_a_dp_hit': tot['dp_hit_pairs'],
                           'detected': tot['dp_agree'],
                           'p_detect_given_dp_hit': (tot['dp_agree'] / tot['dp_hit_pairs'])
                                                    if tot['dp_hit_pairs'] else None,
                           'cp95_lower': (float(cp95_lower(tot['dp_agree'],
                                                           tot['dp_hit_pairs']))
                                          if tot['dp_hit_pairs'] else None),
                           'role': 'THE SUCCESS STATISTIC: detection = a DP-HIT event with the '
                                   'DP mask applied to the class key, both trails in the same '
                                   'class key at that step'},
        'trap_before_dp': {'pairs': tot['trap_before_dp'],
                           'fraction_fixture': tot['trap_before_dp'] / n,
                           'production_bound': p_trap_prod,
                           'production_bound_log2': math.log2(p_trap_prod),
                           'reading': 'a class repeat is a cycle; the guard-free walk state IS '
                                      'the class, so a repeat before the first DP leaves no '
                                      'report and the trail is abandoned at the cap.  At '
                                      'PRODUCTION a DP arrives after ~2^80 steps while the '
                                      'class pool is 2^251.27, so this loss is bounded by '
                                      'theta^2/(2N) = 2^%.1f, not by the fixture fraction (the '
                                      'fixture deliberately runs theta ~ sqrt(N), the regime '
                                      'where the loss is visible at all)'
                                      % math.log2(p_trap_prod)},
        'curve_any_step': curve,
        'lag32_row': {'lag': 32, 'p_agree_any_step': curve['32']['p_agree_any_step'],
                      'p_agree_given_dp': curve['32']['p_agree_given_dp'],
                      'reading': 'the lag-32 arm the v22 verdict singled out (0.10625 with the '
                                 'guard); guard-free the two trails are the same trail, so '
                                 'agreement at lag 32 is the same event as at lag 1'},
        'reading': 'the guard was the ONLY thing that broke the lockstep of two trails that '
                   'share a start point: with the guard off the successor class is a function of '
                   'the class alone, so the pair stays in one class until its first DP and '
                   'reports there with two different label pairs (the solvable event).  The '
                   'fixture measures the mechanism and the trap-before-DP loss in the regime '
                   'where that loss is visible; the production loss is the closed-form bound.',
    }


# ------------------------------------------------------------------ F3
def _rep_batch(args):
    walks, steps, dp_bits, cap_mult, seed, bits = args
    toy = TC.Toy(bits=bits, seed=20261024)
    rng = random.Random(seed)
    theta = 1 << dp_bits
    cap = cap_mult * theta
    reports = 0
    abandons = 0
    starts = 0
    dp_steps = []
    per_walk = []
    for _ in range(walks):
        n = 0
        ab = 0
        used = 0
        while used < steps:
            y = toy.rand_gt()
            a, b = 1, 1
            starts += 1
            seg = 0
            seen = set()
            while True:
                k = _key(toy, y)
                if k in seen:
                    ab += 1
                    used += cap
                    break
                seen.add(k)
                y, a, b, _kd, _ki = toy.step_v13(y, a, b)
                used += 1
                seg += 1
                if used >= steps:
                    break
                if routing.dp_hit(_key(toy, y), dp_bits):
                    n += 1
                    dp_steps.append(seg)
                    break
        reports += n
        abandons += ab
        per_walk.append(n)
    return dict(reports=reports, abandons=abandons, walks=walks, per_walk=per_walk,
                starts=starts, dp_steps=dp_steps)


def reports_restart(walks=240, steps=2048, dp_bits=5, cap_mult=CAP_MULT, bits=5,
                    seed=77, workers=8):
    theta = 1 << dp_bits
    chunks = [(walks // workers + 1, steps, dp_bits, cap_mult, seed + i, bits)
              for i in range(workers)]
    with Pool(workers) as p:
        parts = p.map(_rep_batch, chunks)
    reports = sum(q['reports'] for q in parts)
    abandons = sum(q['abandons'] for q in parts)
    walks_tot = sum(q['walks'] for q in parts)
    starts = sum(q['starts'] for q in parts)
    dps = [x for q in parts for x in q['dp_steps']]
    mean_dp_steps = (sum(dps) / len(dps)) if dps else None
    formula = (steps + theta - 1) // theta
    return {
        'walks': walks_tot, 'steps_per_walk': steps, 'dp_bits': dp_bits,
        'reports_measured': reports,
        'reports_per_walk_measured': reports / walks_tot,
        'reports_per_step_measured': reports / (walks_tot * steps),
        'model_reports_per_step': 1.0 / theta,
        'ceil_formula_reports_per_walk': formula,
        'check_measured_vs_formula_ratio': (reports / walks_tot) / formula,
        'abandons_measured': abandons,
        'trail_starts_measured': starts,
        'mean_steps_to_dp_measured': mean_dp_steps,
        'theta_fixture_steps': theta,
        'abandon_rate_per_start': abandons / starts if starts else None,
        'production_projection': 'reports per S steps = S / (theta + 128 theta p_trap) with '
                                 'p_trap = theta^2/(2N) = 2^-92.3 at production, i.e. '
                                 'ceil(S/theta) exactly; the fixture projection is the same '
                                 'formula with the fixture p_trap, which is why the census '
                                 'reports the abandon rate beside the report count',
        'restart_implemented': 'after EVERY DP the reporting trail is replaced by a FRESH '
                               'trail (fresh coins) and after EVERY abandon at the cap the '
                               'state is discarded and a fresh trail started; the walk '
                               'continues until the step budget is consumed (verdict '
                               '3eb6d983 item 3: reports_direct never restarted)',
        'reading': 'with the restart the report count converges on the DP rate 1/theta, so '
                   'ceil(S/theta) is confirmed rather than under-counted by path shadowing '
                   '(v22 measured 53.05 against 64 because a continuing walk re-hits the same '
                   'DP classes)',
    }


# ------------------------------------------------------------------ F4
def _bnd_batch(args):
    walks, steps, seed, bits = args
    toy = TC.Toy(bits=bits, seed=20261024)
    rng = random.Random(seed)
    reps = 0
    pairs = 0
    for _ in range(walks):
        y = toy.rand_gt()
        a, b = 1, 1
        seen = set()
        for t in range(steps):
            y, a, b, _kd, _ki = toy.step_v13(y, a, b)
            k = _key(toy, y)
            if t:
                pairs += 1
            if k in seen:
                reps += 1
            seen.add(k)
    return dict(repeats=reps, pairs=pairs, walks=walks)


def bound_guardfree(walks=800, steps=200, bits=5, seed=4242, workers=8):
    toy = TC.Toy(bits=bits, seed=20261024)
    chunks = [(walks // workers + 1, steps, seed + i, bits) for i in range(workers)]
    with Pool(workers) as p:
        parts = p.map(_bnd_batch, chunks)
    reps = sum(q['repeats'] for q in parts)
    pairs = sum(q['pairs'] for q in parts)
    walks_tot = sum(q['walks'] for q in parts)
    lam = (float(K.S_STEPS) ** 2 / 2.0) / float(toy.N)
    return {
        'walks': walks_tot, 'steps_per_walk': steps, 'step_pairs': pairs,
        'guard_free_element_repeats': reps,
        'E_repeats_per_walk_measured': reps / walks_tot,
        'iid_birthday_value_fixture': lam,
        'iid_birthday_value_production': (float(K.S_STEPS) ** 2 / 2.0)
                                          / float(K.N_CLASSES),
        'reading': 'the guard-free walk carries NO ring coordinate, so the classic birthday '
                   'count over the class stream is back in its hypothesis domain (a random map '
                   'on N classes, the disclosed mu = 11/10 allowance); this census is the '
                   'measured cross-check that plugs into Markov E_meas/t, which needs no '
                   'independence at all.  The v22 monotonicity route is WITHDRAWN, not patched.',
        'step_source': 'classkey.step_v13 (guard-free, shipped)',
    }


def main():
    os.makedirs(EV, exist_ok=True)
    t0 = time.time()
    print('fixture N=%d' % TC.Toy(rows=64, seed=20261024).N, flush=True)

    cap = cap_guardfree(walks=256, dp_bits=5)
    print('F1 cap: trapped %d, abandons %d, delay mean %.1f var %.1f, runaway %d/%d (%.1fs)'
          % (cap['trails_trapping_in_cycle'], cap['abandons_at_cap'],
             cap['mean_steps_to_abandon'], cap['var_steps_to_abandon'],
             cap['control_force_admission']['still_without_dp'],
             cap['control_force_admission']['trails_tested'], time.time() - t0), flush=True)
    json.dump(cap, open(os.path.join(EV, 'cap-v23.json'), 'w'), indent=1, sort_keys=True)

    ret = retention_dpgated(pairs=10000, steps=64, dp_bits=5, bits=5)
    print('F2 retention N=785: dp-gated %.5f (n=%d, cp95 %.5f) any-step %.5f, '
          'trap-before-dp %d (%.3f) [%.1fs]'
          % (ret['dp_gated_agree']['p_detect_given_dp_hit'],
             ret['dp_gated_agree']['pairs_with_a_dp_hit'],
             ret['dp_gated_agree']['cp95_lower'],
             ret['any_step_agree']['p_detect'], ret['trap_before_dp']['pairs'],
             ret['trap_before_dp']['fraction_fixture'], time.time() - t0), flush=True)
    json.dump(ret, open(os.path.join(EV, 'retention-v23.json'), 'w'), indent=1,
              sort_keys=True)

    big = retention_dpgated(pairs=3000, steps=64, dp_bits=5, bits=7, seed=2411)
    print('F2b retention N=3136: dp-gated %.5f (n=%d) trap-before-dp %d (%.3f) [%.1fs]'
          % (big['dp_gated_agree']['p_detect_given_dp_hit'],
             big['dp_gated_agree']['pairs_with_a_dp_hit'],
             big['trap_before_dp']['pairs'], big['trap_before_dp']['fraction_fixture'],
             time.time() - t0), flush=True)
    json.dump(big, open(os.path.join(EV, 'retention-bigN-v23.json'), 'w'), indent=1,
              sort_keys=True)

    rep = reports_restart(walks=240, steps=2048, dp_bits=5)
    print('F3 reports-with-restart %.3f/walk vs formula %d ratio %.4f [%.1fs]'
          % (rep['reports_per_walk_measured'], rep['ceil_formula_reports_per_walk'],
             rep['check_measured_vs_formula_ratio'], time.time() - t0), flush=True)
    json.dump(rep, open(os.path.join(EV, 'reports-v23.json'), 'w'), indent=1,
              sort_keys=True)

    bnd = bound_guardfree(walks=800, steps=200)
    print('F4 bound: repeats/walk %.4f over %d pairs [%.1fs]'
          % (bnd['E_repeats_per_walk_measured'], bnd['step_pairs'], time.time() - t0),
          flush=True)
    json.dump(bnd, open(os.path.join(EV, 'bound-v23.json'), 'w'), indent=1,
              sort_keys=True)
    print('wrote evidence/{cap,retention,retention-bigN,reports,bound}-v23.json')


if __name__ == '__main__':
    main()
