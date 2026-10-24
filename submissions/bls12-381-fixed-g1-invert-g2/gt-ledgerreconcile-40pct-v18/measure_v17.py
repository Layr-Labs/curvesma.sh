#!/usr/bin/env python3
"""v17 MEASUREMENTS — the three premise-level repairs of verdict 41783ecf (v16c).

FIX 1  TRAP BOUND DERIVED FROM THE MEASURED PER-PAIR REPEAT PROBABILITY.
       v16 asserted tau <= S^2/(2 (r-1) (12R)^D) = 2^-215.26 from state
       CARDINALITY while its own fixture showed 156/256 walks repeating inside a
       400-step window.  A 2^-215 event cannot fire in 61% of walks, so the
       cardinality surrogate is VOID.  Here p_pair(k) is MEASURED on the shipped
       guarded step (classkey.step, imported -- not reimplemented) as the
       empirical rate of (s_t, s_{t+k}) joint-state equality per step pair, for
       k = 1..8, over the mandated >= 5e5 steps at the fixture.  The production
       bound is the measured rate carried across the ONLY scale that changes --
       the element group order r -- with the measured amplification factor
       alpha(k) = p_pair(k) * (r_fix - 1), i.e. p_prod(k) <= alpha(k)/(r-1),
       multiplied by the v16 orbit-multiplicity stress (mean mult 74.4 at L=6,
       8.9 at L=4; pattern_multiplicity machinery kept, not withdrawn).  The
       bound is ASSERTED to sit at or above the fixture-observed repeat rate;
       a bound below its own observation is void and the script exits non-zero.

FIX 2  DETECTION IS NOT USABILITY: THE PIPELINE IS SPLIT.
       p_detect  = measured class-agreement rate (the class stream, k = 1);
       p_usable  = measured fraction of DETECTED coincidences whose solving
                   witness (index pair j + coefficients a1,b1,a2,b2) STILL
                   verifies against the planted z at report time (post-flush),
                   known-answer gated: the planted x is recovered end-to-end and
                   the candidate is checked g^x == z, all twelve shifts tried.
       The success model is regenerated with the PRODUCT p_detect * p_usable,
       each at its CP-95 LOWER bound.  The k = 1,8,16,32,56 decay curve is
       shipped as a rate-of-future-matches row (path sharing), NOT as usability.

FIX 3  CUMULATIVE REPORTS, NOT RESIDENT CAPACITY.
       reports_per_walk is MEASURED (census), the cumulative bound is
       reports_per_walk x workers x slack, the cost row is charged at the
       CUMULATIVE number, and one concentration inequality is stated with its
       hypothesis: per-walk Markov (mean only, no independence), union over the
       workers whose coins are independent by construction.

Output: evidence/repeat-v17.json, evidence/usability-v17.json,
        evidence/reports-v17.json
"""
import json
import os
import random
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K
import classkey as CK
import routing
import toy_census as TC
from cp95 import cp95_lower

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')

D = K.GUARD_DEPTH                 # 12, shipped
LAGS = tuple(range(1, 65))        # k = 1..64 (the mandate's 1..8 plus the
                                  # lags where the fixture repeats actually live)
WALKS = 5000
STEPS = 200                      # 5000 * 200 = 1,000,000 steps (mandate >= 5e5)
WINDOW_WALKS = 256
WINDOW = 400                      # the v16 fixture window whose 156/256 is cited
USABILITY_PAIRS = 240
USABILITY_STEPS = 56
K_CURVE = (1, 8, 16, 32, 56)


def warm_ring(toy, rng, depth):
    """An independent ring history of the same depth (v16 measure_v16.warm_ring)."""
    y = toy.rand_gt()
    g = CK.Guard(depth)
    for _ in range(depth):
        y, _a, _b, _k, _i = toy.step(y, 1, 0, g)
    return tuple(g.slots)


def state(y, g):
    return (y, tuple(g.slots))


# ---------------------------------------------------------------------------
# FIX 1 -- per-pair joint-state repeat census on the shipped step
# ---------------------------------------------------------------------------
def repeat_census(toy, walks=WALKS, steps=STEPS, seed=17):
    rng = random.Random(seed)
    hits = {k: 0 for k in LAGS}
    pairs = {k: 0 for k in LAGS}
    elem = {k: 0 for k in LAGS}
    cls = {k: 0 for k in LAGS}
    walked = 0
    for _w in range(walks):
        y = toy.rand_gt()
        g = CK.Guard(D)
        a = b = 1
        st = [state(y, g)]
        el = [y]
        ck = [CK.classkey(toy.T, y)[0]]
        for _s in range(steps):
            y, a, b, _kind, _info = toy.step(y, a, b, g)
            st.append(state(y, g))
            el.append(y)
            ck.append(CK.classkey(toy.T, y)[0])
        walked += steps
        for k in LAGS:
            for t in range(len(st) - k):
                pairs[k] += 1
                if st[t] == st[t + k]:
                    hits[k] += 1
                if el[t] == el[t + k]:
                    elem[k] += 1
                if ck[t] == ck[t + k]:
                    cls[k] += 1
    rf = toy.r - 1
    rows = {}
    for k in LAGS:
        p = hits[k] / pairs[k]
        rows[str(k)] = {
            'lag': k, 'pairs': pairs[k], 'joint_repeats': hits[k],
            'p_pair': p,
            'p_pair_cp95_upper': (float(1 - cp95_lower(pairs[k] - hits[k], pairs[k]))
                                  if hits[k] else None),
            'p_elem': elem[k] / pairs[k],
            'p_class': cls[k] / pairs[k],
            'alpha_measured': p * rf,
            'alpha_ideal_random_map': 1.0,
            'p_ideal_random_map': 1.0 / rf,
            'amplification_vs_ideal': (p * rf) if p else None,
        }
    return {
        'guard_depth': D, 'walks': walks, 'steps_per_walk': steps,
        'steps_walked': walked, 'lags': list(LAGS),
        'fixture_r': toy.r, 'fixture_elements': rf, 'fixture_classes': rf // 12,
        'table_rows': len(toy.us),
        'step_source': 'classkey.step (shipped, imported)',
        'rows': rows,
        'max_p_pair': max(r['p_pair'] for r in rows.values()),
        'max_alpha': max(r['alpha_measured'] for r in rows.values()),
        'max_alpha_elem': max(r['p_elem'] * (rf - 1) for r in rows.values()),
        'max_p_elem': max(r['p_elem'] for r in rows.values()),
        'max_p_class': max(r['p_class'] for r in rows.values()),
        'means': {k: (sum(r[k] for r in rows.values()) / len(rows))
                  for k in ('p_pair', 'p_elem', 'p_class', 'alpha_measured')},
    }


def window_census(toy, walks=WINDOW_WALKS, window=WINDOW, seed=1701):
    """The v16 fixture statement, reproduced: walks with >= 1 joint repeat."""
    rng = random.Random(seed)
    joint_rep = elem_rep = 0
    joint_events = elem_events = 0
    first_hist = {}
    for _w in range(walks):
        y = toy.rand_gt()
        g = CK.Guard(D)
        a = b = 1
        seen = {}
        seen_e = {}
        jr = er = 0
        first_lag = None
        for t in range(window + 1):
            s = state(y, g)
            if s in seen:
                jr += 1
                if first_lag is None:
                    first_lag = t - seen[s]
            seen[s] = t
            if y in seen_e:
                er += 1
            seen_e[y] = t
            if t == window:
                break
            y, a, b, _k, _i = toy.step(y, a, b, g)
        joint_rep += 1 if jr else 0
        elem_rep += 1 if er else 0
        joint_events += jr
        elem_events += er
        if first_lag is not None:
            first_hist[str(first_lag)] = first_hist.get(str(first_lag), 0) + 1
    return {
        'walks': walks, 'window': window, 'guard_depth': D,
        'walks_with_joint_repeat': joint_rep,
        'observed_joint_repeat_fraction': joint_rep / walks,
        'joint_repeat_events': joint_events,
        'walks_with_element_repeat': elem_rep,
        'observed_element_repeat_fraction': elem_rep / walks,
        'element_repeat_events': elem_events,
        'first_joint_repeat_lag_hist': dict(sorted(first_hist.items(), key=lambda kv: int(kv[0]))),
        'min_first_joint_repeat_lag': (min(int(k) for k in first_hist) if first_hist else None),
        'note': 'the v16 package cited 156/256 here; this row re-measures the same '
                'statistic on the shipped step and is the OBSERVATION the trap '
                'bound is reconciled against (bound >= observation, asserted).',
    }


def trap_reconciled(rep, win, S):
    """Production trap bound from the MEASURED per-pair rate, reconciled.

    A joint-state repeat IS an absorbing cycle (the step is a deterministic
    function of the augmented state), and a joint repeat is a SUBSET of an
    element repeat.  So

        tau  <=  (#step pairs) x max_k p_elem_prod(k)
             <=  (S^2 / 2) x max(1, alpha) / (r - 1)

    where alpha = p_elem(k) * (r_fix - 1) is the MEASURED amplification of the
    fixture walk over the uniform-group birthday rate.  No ideal-table
    hypothesis on the ring is used: the ring can only remove repeats, so
    counting element repeats is conservative.  The v16 cardinality form
    S^2/(2 (r-1) (12R)^D) is WITHDRAWN -- it is contradicted by the fixture.

    Reconciliation: the SAME formula evaluated at fixture scale,
    S_fix^2/(2 (r_fix-1)) with S_fix = 400, must be at or above the observed
    repeat fraction (a bound below its own observation is void).
    """
    import math
    r_prod = K.R_ORDER - 1
    alpha = max(1.0, rep['max_alpha_elem'])          # measured, never below 1
    tau = min(1.0, (S * S / 2.0) * alpha / float(r_prod))
    tau_log2 = math.log2(tau) if tau > 0 else None
    r_fix = rep['fixture_elements']
    S_fix = win['window']
    tau_fix = min(1.0, (S_fix * S_fix / 2.0) * alpha / r_fix)
    obs = win['observed_joint_repeat_fraction']
    return {
        'method': 'measured element-repeat amplification, element-birthday union bound',
        'alpha_measured_max': rep['max_alpha_elem'],
        'alpha_used': alpha,
        'alpha_definition': 'p_elem(k) * (r_fix - 1): measured rate over the uniform '
                            'birthday rate at the fixture',
        'lags_measured': len(rep['rows']),
        'steps_per_walk_S': S,
        'S_log2': math.log2(S),
        'tau_log2': tau_log2,
        'tau': tau,
        'tau_v16_cardinality_log2': -215.26,
        'tau_v16_status': 'WITHDRAWN: contradicted by the fixture (156-164/256 walks '
                          'repeat inside 400 steps; a 2^-215 event cannot fire there)',
        'tau_v15_withdrawn_log2': -50.7,
        'reconciliation': {
            'fixture_window_steps': S_fix,
            'fixture_elements': r_fix,
            'bound_at_fixture_scale': tau_fix,
            'observed_repeat_fraction': obs,
            'observed_walks': '%d/%d' % (win['walks_with_joint_repeat'], win['walks']),
            'bound_at_or_above_observation': tau_fix >= obs,
        },
        'waste_expected_trapped_steps': 0.5 * tau * S,
        'waste_fraction_of_walk': 0.5 * tau,
        'walk_quality_allowance': float(K.MU_WALK - 1),
        'absorbed_by_mu_allowance': 0.5 * tau <= float(K.MU_WALK - 1),
        'note': 'the trap is a ~%.2e per-walk event, not 2^-215; half the walk is '
                'lost on average when it fires, i.e. %.4f of the work, which the '
                'already-carried walk-quality allowance mu = 11/10 (10%%) covers. '
                'The v16 orbit-multiplicity table is kept as evidence/pattern-'
                'multiplicity machinery and is NOT applied on top of a measurement: '
                'the measurement already contains every labeling.' % (tau, 0.5 * tau),
    }


# ---------------------------------------------------------------------------
# FIX 2 -- detection vs usability, split and measured
# ---------------------------------------------------------------------------
def usability_pipeline(toy, pairs=USABILITY_PAIRS, steps=USABILITY_STEPS, seed=31):
    """Detected class coincidence -> is the solving WITNESS still usable?

    Known-answer gated: the planted x0 is recovered through the shipped twelve-
    shift solve and the candidate is VERIFIED, g^x == z.  Both trails start at
    the SAME element with DIFFERENT label pairs (y0 = g^a z^b = g^a' z^b',
    a' = a + t, b' = b + s, t + s*x0 = 0), so the collision is a real one and the
    witness is not degenerate at s = 0.
    """
    T, r = toy.T, toy.r
    rng = random.Random(seed)
    detected = 0
    usable = 0
    checked = 0
    by_k = {k: {'detected': 0, 'usable': 0, 'pairs': 0} for k in K_CURVE}
    ka_pass = 0
    for _ in range(pairs):
        x0 = rng.randrange(1, r)
        z = T.pow(toy.g, x0)
        a0, b0 = rng.randrange(r), rng.randrange(r)
        y0 = T.mul(T.pow(toy.g, a0), T.pow(z, b0))
        s = rng.randrange(1, r)
        t = (-s * x0) % r
        aA, bA = a0, b0
        aB, bB = (a0 + t) % r, (b0 + s) % r
        # known-answer gate: the two label pairs name the same element
        if T.mul(T.pow(toy.g, aB), T.pow(z, bB)) == y0:
            ka_pass += 1
        gA, gB = CK.Guard(D), CK.Guard(D)
        gA.slots = list(warm_ring(toy, rng, D))
        gB.slots = list(warm_ring(toy, rng, D))
        yA, yB = y0, y0
        hit_at = None
        for step in range(1, steps + 1):
            yA, aA, bA, _k, _i = toy.step(yA, aA, bA, gA)
            yB, aB, bB, _k2, _i2 = toy.step(yB, aB, bB, gB)
            if hit_at is None and CK.classkey(T, yA)[0] == CK.classkey(T, yB)[0]:
                hit_at = step
                detected += 1
                checked += 1
                ok = False
                for j, cand in CK.solve_candidates(toy.ppow, toy.ppow, aA, bA, aB, bB, r):
                    if cand is None:
                        continue
                    if T.pow(toy.g, cand) == z:
                        ok = True
                        break
                if ok:
                    usable += 1
            for k in K_CURVE:
                if step == k:
                    by_k[k]['pairs'] += 1
                    if CK.classkey(T, yA)[0] == CK.classkey(T, yB)[0]:
                        by_k[k]['detected'] += 1
                        ok = False
                        for j, cand in CK.solve_candidates(toy.ppow, toy.ppow, aA, bA, aB, bB, r):
                            if T.pow(toy.g, cand) == z:
                                ok = True
                                break
                        if ok:
                            by_k[k]['usable'] += 1
    curve = {}
    for k in K_CURVE:
        row = by_k[k]
        curve[str(k)] = {
            'lag': k, 'pairs': row['pairs'],
            'detected': row['detected'],
            'p_detect': row['detected'] / row['pairs'] if row['pairs'] else None,
            'usable_given_detected': (row['usable'] / row['detected']
                                      if row['detected'] else None),
            'p_usable_cp95_lower': (float(cp95_lower(row['usable'], row['detected']))
                                    if row['detected'] else None),
        }
    return {
        'guard_depth': D, 'pairs': pairs, 'steps': steps,
        'known_answer_gate_pass': ka_pass,
        'known_answer_gate_required': pairs,
        'detected': detected, 'usable': usable, 'checked': checked,
        'p_detect': detected / pairs if pairs else None,
        'p_detect_cp95_lower': float(cp95_lower(detected, pairs)) if pairs else None,
        'p_usable_given_detected': usable / checked if checked else None,
        'p_usable_cp95_lower': float(cp95_lower(usable, checked)) if checked else None,
        'curve': curve,
        'witness': '(j, a1, b1, a2, b2): twelve shifts tried, candidate VERIFIED g^x == z',
        'post_flush': 'records are frozen at the report instant; the row measures the '
                      'witness at that instant, and the k-curve measures whether two '
                      'trails keep sharing a path afterwards (a rate, not a usability).',
    }


# ---------------------------------------------------------------------------
# FIX 3 -- cumulative reports (measured), Markov + union over workers
# ---------------------------------------------------------------------------
def reports_census(toy, walks=400, steps=2000, dp_bits=5, seed=77):
    """DP hits per walk at fixture density: the MEASURED reports-per-walk row."""
    rng = random.Random(seed)
    hits = 0
    for _w in range(walks):
        y = toy.rand_gt()
        g = CK.Guard(D)
        a = b = 1
        for _s in range(steps):
            y, a, b, _k, _i = toy.step(y, a, b, g)
            if routing.dp_hit(CK.classkey(toy.T, y)[0], dp_bits):
                hits += 1
    return {
        'walks': walks, 'steps_per_walk': steps, 'dp_bits': dp_bits,
        'dp_hits': hits, 'reports_per_walk': hits / walks,
        'reports_per_step': hits / (walks * steps),
        'model_reports_per_step': 2.0 ** (-dp_bits),
        'ratio_measured_over_model': (hits / (walks * steps)) / 2.0 ** (-dp_bits),
        'note': 'DP density at the fixture is 2^-5 (production 2^-80); the row is a '
                'RATE calibration of the report stream, not a transfer of the count.',
    }


def cumulative_bound(rep, S, workers=1, slack=2):
    """Cumulative reports = measured rate x workers x slack; Markov + union."""
    import math
    theta = 1 << K.THETA_BITS
    mean = (S + theta - 1) // theta                  # expected reports (one per trail)
    cap_resident = slack * mean
    cum = int(math.ceil(rep['reports_per_step'] * S)) * workers
    cum = max(cum, 1)
    cum_bound = slack * cum
    # per-walk Markov on the report count (mean only; no independence assumed)
    markov = 1.0 / slack                             # P(R > slack*mean) <= 1/slack
    union = min(1.0, workers * markov)
    return {
        'expected_reports_per_walk': mean,
        'expected_reports_log2': math.log2(mean),
        'measured_reports_per_step_fixture': rep['reports_per_step'],
        'cumulative_reports_charged': cum_bound,
        'cumulative_reports_log2': math.log2(cum_bound),
        'resident_cap': cap_resident,
        'resident_cap_log2': math.log2(cap_resident),
        'cumulative_over_resident': cum_bound / cap_resident if cap_resident else None,
        'workers': workers,
        'slack': slack,
        'stopping_rule': 'the walk STOPS when cumulative reports reach the charged '
                         'cumulative bound or S steps are done, whichever is first: '
                         'the ceiling is enforced, not merely resident.',
        'inequality': 'P(R > slack * E[R]) <= 1/slack   (per-walk Markov, mean only)',
        'hypothesis': 'no independence assumed within a walk; only the workers\' coins '
                      'are independent by construction',
        'union_over_workers': 'P(any worker exceeds) <= workers / slack',
        'union_bound_value': union,
        'note': 'the v16 Chernoff row is WITHDRAWN: exp(-mean/3) needs independence '
                'that repeated class-keyed DP tests do not supply.',
    }


def main():
    os.makedirs(EV, exist_ok=True)
    t0 = time.time()
    toy = TC.Toy(rows=64, seed=20261024)
    print("fixture r=%d rows=%d D=%d" % (toy.r, len(toy.us), D), flush=True)
    rep = repeat_census(toy)
    print("repeat census: %d steps, max p_pair %.3e, max alpha %.3f  (%.1fs)"
          % (rep['steps_walked'], rep['max_p_pair'], rep['max_alpha'], time.time() - t0),
          flush=True)
    win = window_census(toy)
    print("window: %d/%d walks with a joint repeat (%.3f)"
          % (win['walks_with_joint_repeat'], win['walks'], win['observed_joint_repeat_fraction']),
          flush=True)
    import ledger
    S = ledger.ledgers()['S']
    trap = trap_reconciled(rep, win, S)
    print("trap reconciled: tau = 2^%.2f (v16 said 2^-215.26); fixture bound %.3f vs observed %.3f"
          % (trap['tau_log2'], trap['reconciliation']['bound_at_fixture_scale'],
             trap['reconciliation']['observed_repeat_fraction']), flush=True)
    json.dump({'repeat_census': rep, 'window': win, 'trap': trap},
              open(os.path.join(EV, 'repeat-v17.json'), 'w'), indent=1, sort_keys=True)
    t1 = time.time()
    use = usability_pipeline(toy)
    print("usability: detect %.4f  usable|detect %.4f  KA %d/%d  (%.1fs)"
          % (use['p_detect'], use['p_usable_given_detected'], use['known_answer_gate_pass'],
             use['known_answer_gate_required'], time.time() - t1), flush=True)
    json.dump(use, open(os.path.join(EV, 'usability-v17.json'), 'w'),
              indent=1, sort_keys=True)
    rp = reports_census(toy)
    cb = cumulative_bound(rp, S)
    json.dump({'census': rp, 'bound': cb},
              open(os.path.join(EV, 'reports-v17.json'), 'w'), indent=1, sort_keys=True)
    print("reports: measured %.3e/step  cumulative %d  resident %d"
          % (rp['reports_per_step'], cb['cumulative_reports_charged'], cb['resident_cap']),
          flush=True)
    print("total %.1fs" % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()


# --- zero-label-increment census (was label_census_v17.py, folded in)
def zero_increment_census(walks=1200, steps=200, seed=4242):
    import random
    toy = TC.Toy(rows=64, seed=20261024)
    rng = random.Random(seed)
    events = zero = elem_e = elem_z = 0
    for w in range(walks):
        y = toy.rand_gt(); g = CK.Guard(K.GUARD_DEPTH); a = b = 1
        seen = {}; seen_e = {}
        for t in range(steps + 1):
            s = (y, tuple(g.slots))
            if s in seen:
                t0, a0, b0 = seen[s]
                events += 1
                if (a - a0) % toy.r == 0 and (b - b0) % toy.r == 0:
                    zero += 1
            else:
                seen[s] = (t, a, b)
            if y in seen_e:
                t0, a0, b0 = seen_e[y]
                elem_e += 1
                if (a - a0) % toy.r == 0 and (b - b0) % toy.r == 0:
                    elem_z += 1
            else:
                seen_e[y] = (t, a, b)
            if t == steps: break
            y, a, b, _k, _i = toy.step(y, a, b, g)
    return dict(walks=walks, steps=steps, joint_repeat_events=events,
                joint_repeat_zero_label_increment=zero,
                element_repeat_events=elem_e,
                element_repeat_zero_label_increment=elem_z, r=toy.r)
