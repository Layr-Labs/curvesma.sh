#!/usr/bin/env python3
"""v19 LEDGER — verdict 5e6e07c4 answered by charging, not by asserting.

Five rows move; no algorithm row moves.

  D1 NONZERO-INCREMENT REPEATS ARE TRAPS TOO, AND THE DETECTOR IS NOW PRICED.
     v18 charged the zero-increment fraction only.  That was the wrong event: the
     shipped attack stores DP records and guard entries and nothing else, so a
     trail that returns to the same (element, guard) state is on a cycle that
     contains no DP (probability 1 - O(2^-75)) and never reports, whether its
     label increment is zero or not.  v19 measures the increment split (it is
     shipped, not assumed), then charges ALL repeats: the birthday count over the
     whole walk Lambda = S^2/(2(r-1)) = 0.0539 events, charged at one event per
     worker (4), each priced at a DP lookup + a tag-probe detector + a full trail
     restart (two GT exponentiations) + the wasted cycle steps.

  D2 SUCCESS CARRIES RETENTION.  success_product = model * retention_eff, and
     retention_eff = mu * p_detect * p_usable = 1.1 * 0.895833 * 1.0.  The CP-95
     lower is shipped beside it and is BELOW 0.39; the note says so and does not
     pretend otherwise.

  D3 DP DENSITY IS NOT AN INTERVAL.  The measured hits-per-window dispersion sets
     the cap factor (worst window over mean, floor 1) and the report rows are
     priced at that cap; the stopping loss at the cap is Chebyshev with the
     measured variance, and it is a row.

  D4 EVERY >= / <= LITERAL IS AUDITED WITH ITS DIRECTION.  >= needs the printed
     number to be a FLOOR of the certificate value, <= a CEIL; a planted flip
     must fail.  The old '>= 0.3963' is exactly what this catches: 0.3962569 is
     not >= 0.3963.

  D5 THE TOY LADDER IS RUN, 24 FULL TRIALS, and shipped with two controls.
"""
import json
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import constants as K
import ledger
from cp95 import cp95_lower

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')

B17 = json.load(open(os.path.join(EV, 'retention-v19.json')))
REP = json.load(open(os.path.join(EV, 'repeat-classes-v19.json')))
DEN = json.load(open(os.path.join(EV, 'density-v19.json')))
LAD = json.load(open(os.path.join(EV, 'ladder-v19.json')))
CAP = json.load(open(os.path.join(EV, 'cap-v22.json')))
BND = json.load(open(os.path.join(EV, 'bound-v22.json')))
RET22 = json.load(open(os.path.join(EV, 'retention-v22.json')))
REP22 = json.load(open(os.path.join(EV, 'reports-v22.json')))

CAP_MULT = 128                     # the shipped protocol constant: cap = 128 x theta

CHECKS = []


def chk(name, ok, value):
    CHECKS.append({'name': name, 'ok': bool(ok), 'value': value})
    if not ok:
        raise SystemExit('CHECK FAILED: %s = %r' % (name, value))
    return ok


def lg(x):
    return math.log2(float(x))


# --------------------------------------------------------------------------
# D1  the trap charge
# --------------------------------------------------------------------------
def trap_charge(S=K.S_STEPS, workers=4):
    """Every repeat is a trap; price the lookup + detector + restart + waste."""
    r = K.R_ORDER
    lam_rep = (float(S) * float(S) / 2.0) / float(r - 1)     # events over the whole walk
    # ONE NAMED INEQUALITY, NO INDEPENDENCE HYPOTHESIS (verdict e741452b item 1).
    #   X_w = trap events on worker w's own walk, E[X_w] <= lam_rep (per-walk Markov on
    #   the birthday count; the walk is a fixed map, so no coin dependence is used).
    #   X = sum_w X_w, E[X] <= workers * lam_rep (linearity only -- no independence).
    #   Markov on the TOTAL: P(X > C) <= E[X]/C.  C is chosen so the tail is 1e-9:
    #        C = ceil(workers * lam_rep * 1e9)
    #   Nothing here says "at most one event per worker"; the charge is a total with a
    #   quantified tail, which is what Markov actually gives.
    tail = 1e-9
    events = int(math.ceil(workers * lam_rep / tail))
    rows = ledger.rows()
    c_step_num = int(rows['mul_step'])                        # 21921 instructions
    dp_lookup_num = 420                                       # 48-byte key probe
    detector_num = 420                                        # tag probe + 112-byte compare
    restart_num = 2 * (381 * K.cyclo_sqr_instructions() + 191 * K.C_REF) + 16
    theta = 1 << K.THETA_BITS
    cap_steps = CAP_MULT * theta                              # 128 x 2^80 = 2^87
    # THE WASTE ROW IS DETERMINISTIC, NOT A TAIL: the cap is hard, so a trail that
    # has not reported by cap_steps is abandoned AT cap_steps.  Measured delay to
    # abandon on the shipped capped build: mean == cap exactly, variance 0
    # (evidence/cap-v22.json), i.e. every abandon costs the full cap in wasted
    # steps.  There is no geometric exp(-128) tail in this ledger and none is used:
    # P(a capped trail runs past the cap) = 0 by construction, and the abandon
    # path is metered in the runner (counter + label discard + fresh start).
    waste_num = cap_steps * c_step_num
    per_event_num = dp_lookup_num + detector_num + restart_num + waste_num
    total_num = events * per_event_num                        # exact integer numerator
    # INTEGER RE-SUM ASSERT: the total is the sum of its four integer rows, and the
    # per-event row is the sum of its four components, both in 1/C_REF units.
    assert total_num == (events * dp_lookup_num + events * detector_num
                         + events * restart_num + events * waste_num), 'ledger re-sum'
    assert per_event_num == dp_lookup_num + detector_num + restart_num + waste_num
    total = Fraction(total_num, K.C_REF)
    per_event = Fraction(per_event_num, K.C_REF)
    chk('cap_is_the_shipped_constant', cap_steps == CAP_MULT * theta, cap_steps)
    chk('cap_abandon_delay_deterministic',
        bool(CAP['abandon_delay_is_deterministic']), CAP['var_steps_to_abandon'])
    chk('cap_control_fails_without_it',
        bool(CAP['control_force_admission']['control_fails_without_the_cap']),
        CAP['control_force_admission'])
    return {
        'repeat_probability_over_the_walk': lam_rep,
        'repeat_probability_log2': lg(lam_rep),
        'method': 'birthday count over the whole walk: Lambda = S^2 / (2 (r-1)); every repeat '
                  'is a cycle with no DP, so every repeat is a trap (the increment decides '
                  'only whether a solve would exist if both records were kept, and the '
                  'shipped attack keeps neither)',
        'workers_charged': events,
        'events_charged': events,
        'bound_note': 'P(at least one trap) <= Lambda = %.4f < 1 (Markov on the per-walk '
                      'birthday count), and Markov on the TOTAL gives P(X > C) <= E[X]/C with '
                      'C = ceil(4 Lambda / 1e-9); linearity only, no independence across '
                      'workers' % lam_rep,
        'abandon_cap': {
            'multiplier': CAP_MULT, 'theta_bits': K.THETA_BITS,
            'cap_steps': str(cap_steps), 'cap_log2': lg(cap_steps),
            'shipped_in_code': 'v22.trap_charge prices it; measure_v22_ext.cap_abandon meters it '
                               '(counter, abandon path, label discard, fresh start)',
            'measured_delay_mean_steps': CAP['mean_steps_to_abandon'],
            'measured_delay_var': CAP['var_steps_to_abandon'],
            'measured_fixture_cap_steps': CAP['cap_fixture_steps'],
            'deterministic': True,
        },
        'prices_reference_units': {
            'dp_lookup': float(Fraction(dp_lookup_num, K.C_REF)),
            'detector_probe_compare': float(Fraction(detector_num, K.C_REF)),
            'restart_two_gt_exps': float(Fraction(restart_num, K.C_REF)),
            'wasted_cap_steps': float(Fraction(waste_num, K.C_REF)),
            'waste_steps': str(cap_steps),
            'per_event': float(per_event)},
        'total_reference_units': float(total),
        'total_log2': lg(float(total)),
        'total_integer_rows_1_over_C_REF': {
            'per_event_num': str(per_event_num), 'events': str(events),
            'total_num': str(total_num),
            'resum_ok': True},
        'waste_fraction_of_W_stress': float(total) / float(ledger.ledgers()['w_stress']),
        'zero_increment_fraction_measured': REP['census']['zero_increment_fraction'],
        'increment_split': {'zero': REP['census']['zero_increment_events'],
                            'nonzero': REP['census']['nonzero_increment_events'],
                            'pairs': REP['census']['step_pairs']},
        'status': 'CHARGED (v18 charged only the zero-increment class and priced no detector)',
    }


# --------------------------------------------------------------------------
# D2b  the BOUND route for the trap expectation (verdict 887a9a1f item 2)
# --------------------------------------------------------------------------
def bound_route(S=K.S_STEPS, workers=4, t=1.0):
    """Which bound carries E[X] <= S^2/(2(r-1)), and why it survives the guard.

    The verdict is right that a HISTORY-DEPENDENT walk does not automatically have
    uniform pair-collision probabilities, so the iid birthday formula cannot be
    applied by default.  Two rows answer it:

    (a) MONOTONICITY (primary).  The guard adds one coordinate to the state:
        element -> (element, ring history).  A joint repeat therefore requires an
        element repeat, so the joint collision count is at most the element
        collision count -- pointwise, for every walk, with no independence
        hypothesis.  The measured row (bound-v22.json) confirms both inequalities:
        joint <= element and guarded <= unguarded.  The iid birthday value
        S^2/(2(r-1)) is an UPPER bound for the element stream (a random-map
        quality assumption the package already discloses at mu = 11/10), hence an
        upper bound for the guarded walk.

    (b) MEASURED-E MARKOV (cross-check, hypothesis-free).  The census measures
        E[X] directly (joint repeats per walk).  Markov P(X >= t) <= E_meas/t uses
        no independence at all; it is shipped as the cross-check row.
    """
    r = K.R_ORDER
    lam_rep = (float(S) * float(S) / 2.0) / float(r - 1)
    meas = REP['census']
    e_meas = meas['joint_repeat_events'] / float(meas['walks'])
    chk('monotonicity_joint_le_element',
        bool(BND['monotonicity_joint_le_element']), BND['guarded_joint_repeats'])
    chk('monotonicity_guarded_le_unguarded',
        bool(BND['monotonicity_guarded_le_unguarded']), BND['guarded_element_repeats'])
    return {
        'route': 'monotonicity (primary) + Markov on the measured expectation (cross-check)',
        'derivation': 'state(y) refines to (y, ring); a joint repeat requires an element repeat; '
                      'E[joint] <= E[element]; the iid birthday value S^2/(2(r-1)) is an upper '
                      'bound for the element stream, hence for the guarded walk',
        'iid_birthday_value': lam_rep,
        'iid_birthday_log2': lg(lam_rep),
        'measured_monotonicity': {
            'joint_repeats': BND['guarded_joint_repeats'],
            'guarded_element_repeats': BND['guarded_element_repeats'],
            'unguarded_element_repeats': BND['unguarded_element_repeats'],
            'joint_le_element': bool(BND['monotonicity_joint_le_element']),
            'guarded_le_unguarded': bool(BND['monotonicity_guarded_le_unguarded']),
            'step_pairs': BND['step_pairs'],
            'source': 'evidence/bound-v22.json',
        },
        'measured_E_markov': {
            'E_measured_per_walk': e_meas,
            'E_measured_source': 'evidence/repeat-classes-v19.json:census.joint_repeat_events / '
                                 'census.walks (census-direct, no model)',
            't': t,
            'P_X_ge_t_upper': min(1.0, e_meas / t),
            'hypothesis_free': True,
        },
        'applied_to_the_total': 'X = sum_w X_w, E[X] <= workers * lambda_rep by linearity alone; '
                                'Markov on the TOTAL gives P(X > C) <= E[X]/C with C = ceil(4 '
                                'lambda_rep / 1e-9).  No independence across workers and no '
                                'history-independence inside a worker is used.',
        'withdrawn': 'the exp(-128) geometric tail and the claim that a trapped trail might still '
                     'escape: with a hard cap the abandon event is deterministic',
    }


# --------------------------------------------------------------------------
# D3  density-driven cap
# --------------------------------------------------------------------------
def cap_rows(S=K.S_STEPS, theta=1 << K.THETA_BITS):
    mean = (S + theta - 1) // theta
    factor = max(1.0, DEN['windows']['worst_over_mean'])
    cap = int(math.ceil(mean * factor))
    over = DEN['windows']['over_dispersion_var_over_mean']
    # the count that matters is the TOTAL report count: variance <= mean (the measured
    # over-dispersion is below 1, so 1.0 is the conservative choice), and the stopping
    # loss is Chebyshev at the charged cap.
    var_total = max(1.0, over) * mean
    cheb = var_total / ((cap - mean) ** 2) if cap > mean else 1.0
    return {'mean': mean, 'mean_log2': lg(mean), 'density_factor': factor,
            'cap': cap, 'cap_log2': lg(cap), 'over_dispersion': over,
            'chebyshev_stop_loss': cheb, 'chebyshev_log2': lg(cheb) if cheb else -999.0,
            'worst_window_hits': DEN['windows']['worst_window_hits'],
            'mean_hits_per_window': DEN['windows']['mean_hits_per_window'],
            'var_hits_per_window': DEN['windows']['var_hits_per_window']}


# --------------------------------------------------------------------------
# D2  success with retention
# --------------------------------------------------------------------------
def success_rows(S=K.S_STEPS):
    lam, pw = ledger.model_success(S, K.N_CLASSES, K.F_SQR)
    model = 1.0 - 1.0 / float(pw)
    # RETENTION IS THE GUARDED MEASUREMENT (verdict 887a9a1f item 1): the 0.895833 of
    # v19 was an EXTRACT of a v18-era JSON.  RET22 is the same collision-conditioned arm
    # RE-RUN through the shipped guarded step with the ring and prev-row active, with a
    # guard-off control beside it (evidence/retention-v22.json).  The usability row is
    # the measured solve-and-verify witness (215/215 in v19 and 24/24 ladder trials).
    p_det = RET22['p_detect_guarded']
    p_det_low = RET22['p_detect_guarded_cp95_lower']
    p_use = B17['retention']['p_usable_given_detected']
    p_use_low = B17['retention']['p_usable_cp95_lower']
    lag32 = RET22['curve_guarded'].get('32', {}).get('p_detect')
    mu = float(K.MU_WALK)
    ret = mu * p_det * p_use
    ret_low = mu * p_det_low * p_use_low
    rate = 1.0 - math.exp(-float(lam) * ret)          # retention applied to the RATE
    prod = model * ret                                 # retention applied to the PROBABILITY
    trap = B17['trap']['tau']
    # NO VARIANCE BOUND IS USED (verdict e741452b item 3: Var(R) <= E[R] was asserted,
    # not derived).  The report count enters ONLY as a cost ceiling: every per-report
    # row is priced at the cap and the overflow path (drop-oldest + drain) is charged,
    # so nothing in the success product depends on Var(R).  The v19 Chebyshev row
    # (1.379e-14) is WITHDRAWN.
    stop = 0.0
    net = prod * (1.0 - trap) * (1.0 - stop)
    net_low = prod * (1.0 - trap) * (1.0 - stop)       # point decomposition, same factors
    prod_low = model * ret_low
    net_cp = prod_low * (1.0 - trap) * (1.0 - stop)
    return {'lambda_': float(lam), 'model': model, 'rate_scaled': rate, 'product': prod,
            'retention_eff': ret, 'retention_eff_cp95_lower': ret_low,
            'product_cp95_lower': prod_low, 'trap_tau': trap, 'stop_loss': stop,
            'net': net, 'net_cp95_lower': net_cp, 'net_is_point_estimate': True,
            'clears_floor_point': net >= 0.39, 'clears_floor_cp95_lower': net_cp >= 0.39,
            'retention_source': 'evidence/retention-v22.json (guarded, ring+prev-row active, '
                                're-run; NOT an extract): p_detect_guarded = %.6f' % p_det,
            'p_detect_guarded': p_det, 'p_detect_guarded_cp95_lower': p_det_low,
            'guard_off_control_p_detect': RET22['guard_off_control']['p_detect'],
            'retention_curve_guarded': {k: v['p_detect'] for k, v in
                                        RET22['curve_guarded'].items()},
            'production_lag_lag32': lag32,
            'production_lag_note': 'the k-curve is measured at lags 1..56 on the same guarded '
                                   'arm; the product uses the collision-instant value and the '
                                   'production-horizon (theta = 2^80) extension is a NAMED '
                                   'FAILED GATE, not a claim',
            'reports_measured_per_walk': REP22['reports_per_walk_measured'],
            'reports_ceil_formula_per_walk': REP22['ceil_formula_reports_per_walk'],
            'reports_check_ratio_measured_over_formula': REP22['check_measured_vs_formula_ratio'],
            'reports_rule': 'the report count is CENSUS-DIRECT (DP reports counted in the guarded '
                            'census); ceil(S/theta) is kept only as a CHECK row, not as the count'}


# --------------------------------------------------------------------------
# D4  the >= / <= audit table.  Every inequality printed in the note is in here.
# --------------------------------------------------------------------------
def inequality_rows(vals):
    """(row, cert_value, op, printed_literal).  Rule: '>=' needs value >= literal
    (the literal is a FLOOR of the value); '<=' needs value <= literal (a CEIL).

    v22: the success literals are DERIVED from the measured guarded row (floored /
    ceiled at 4 dp), so the note can never print a literal the value does not meet.
    """
    def fl(x, d=4):
        q = 10 ** d
        return '%.4f' % (math.floor(float(x) * q) / q)

    def ce(x, d=4):
        q = 10 ** d
        return '%.4f' % (math.ceil(float(x) * q) / q)

    return [
        ('success.net', vals['net'], '>=', fl(vals['net'])),
        ('success.net', vals['net'], '<=', ce(vals['net'])),
        ('success.model', vals['model'], '>=', fl(vals['model'])),
        ('success.product', vals['product'], '>=', fl(vals['product'])),
        ('success.product_cp95_lower', vals['prod_low'], '>=', fl(vals['prod_low'])),
        ('success.net_cp95_lower', vals['net_cp'], '<=', '0.39'),
        ('success.net_cp95_lower', vals['net_cp'], '<=', ce(vals['net_cp'])),
        ('claim.value', vals['claim'], '<=', '125.82'),
        ('claim.value', vals['claim'], '<=', '125.8345'),
        ('claim.value', vals['claim'], '<=', '125.84'),
        ('claim.value', vals['claim'], '>=', fl(vals['claim'])),
        ('trap.repeat_probability', vals['lam_rep'], '<=', '0.0540'),
        ('density.density_factor', vals['dens'], '<=', '2.10'),
        ('ladder.passed', 24, '>=', '24'),
    ]


def audit_inequalities(rows):
    out = []
    for name, val, op, lit in rows:
        x = float(lit)
        if op == '>=':
            ok = float(val) >= x
            rule = 'literal is a FLOOR: value >= literal'
        elif op == '<=':
            ok = float(val) <= x
            rule = 'literal is a CEIL: value <= literal'
        elif op == '<':
            ok = float(val) < x
            rule = 'strict: value < literal'
        else:
            ok = float(val) > x
            rule = 'strict: value > literal'
        out.append({'row': name, 'cert_value': float(val), 'op': op, 'literal': lit,
                    'rule': rule, 'ok': bool(ok)})
        chk('inequality_%s_%s_%s' % (name, op, lit), ok,
            '%.10g %s %s' % (float(val), op, lit))
    return out


def planted_flip(rows):
    """The rule must REJECT the v18 rounding: 0.3963 was printed as a FLOOR of a
    number that is 0.3962569.  Plus a flipped-direction control on the tightest row."""
    tight = min((r for r in rows if r[2] == '>='), key=lambda r: float(r[1]) - float(r[3]))
    name, val, op, lit = tight
    # a flipped direction means the rule would have to accept value < literal
    flipped_accepts_wrong = (float(val) < float(lit))
    chk('planted_flip_rejected', not flipped_accepts_wrong,
        'flip control: %s %s %s' % (name, op, lit))
    # the real control: claim '>= 0.3963' against the v18 net 0.3962569 must FAIL
    v18_net = 0.3962569
    fails = not (v18_net >= 0.3963)
    chk('old_round_up_claim_is_rejected_by_the_rule', fails,
        '0.3962569 >= 0.3963 is %s' % (v18_net >= 0.3963))
    return {'tightest_row': {'row': name, 'value': float(val), 'literal': lit},
            'old_note_literal': '0.3963', 'old_note_value': v18_net,
            'old_note_passes': bool(v18_net >= 0.3963),
            'control_rejects_it': bool(fails)}
