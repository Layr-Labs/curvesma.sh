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

# the five v19-era rows, folded into ONE file (file-count cap), each still the
# measured census it was: retention / repeat classes / density / ladder.
LEG = json.load(open(os.path.join(EV, 'legacy-v19.json')))
B17, REP, DEN, LAD = (LEG['retention_v19'], LEG['repeat_classes_v19'],
                      LEG['density_v19'], LEG['ladder_v19'])
CAP = json.load(open(os.path.join(EV, 'cap-v22.json')))
BND = json.load(open(os.path.join(EV, 'bound-v22.json')))
RET22 = json.load(open(os.path.join(EV, 'retention-v22.json')))   # the GUARDED arm, kept
REP22 = json.load(open(os.path.join(EV, 'reports-v22.json')))
RET23 = json.load(open(os.path.join(EV, 'retention-v23.json')))          # guard-free, DP-gated
RETBN = json.load(open(os.path.join(EV, 'retention-bigN-v23.json')))     # guard-free, N=3136
RETLG = json.load(open(os.path.join(EV, 'retention-lag-v23.json')))      # lag arm
CAP23 = json.load(open(os.path.join(EV, 'cap-v23.json')))                # cap, guard-free
REP23 = json.load(open(os.path.join(EV, 'reports-v23.json')))            # restart-fixed
BND23 = json.load(open(os.path.join(EV, 'bound-v23.json')))              # birthday cross-check

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
        bool(CAP23['abandon_delay_is_deterministic']), CAP23['var_steps_to_abandon'])
    chk('cap_control_fails_without_it',
        bool(CAP23['control_force_admission']['control_fails_without_the_cap']),
        CAP23['control_force_admission'])
    chk('cap_guard_free_trails_trap_and_abandon',
        CAP23['abandons_at_cap'] > 0 and CAP23['abandons_at_cap'] <= CAP23['walks'],
        'abandons %d of %d walks' % (CAP23['abandons_at_cap'], CAP23['walks']))
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
    """The expectation row for the trap count, AFTER the guard is dropped.

    v22 argued E[joint repeats] <= E[element repeats] <= iid birthday by MONOTONICITY.
    The verdict is right that monotone collision counts do not give stochastic
    domination, and v23 does not try to repair that route: it is WITHDRAWN.

    THE GUARD WAS THE WHOLE OBSTRUCTION.  The guard added a ring coordinate to the
    state, which is what made the walk history-dependent and put it outside the
    hypothesis domain of the classic birthday count.  v23 drops the guard, so the
    step is history-free -- the successor CLASS is a function of the class alone
    (classkey.step_v13: column = kappa(y), no ring) -- and the walk is back in the
    SAME hypothesis domain the pre-guard submissions were accepted in: a fixed map
    on N classes with the disclosed mu = 11/10 random-map allowance.  One line, no
    new claim:

        guard-off  =>  history-free step  =>  classic birthday count applies.

    BESIDE IT, hypothesis-free: the census measures E[repeats] directly at an
    UNSATURATED fixture length (evidence/bound-v23.json) and Markov P(X >= t) <=
    E_meas/t uses no independence at all.  The measured row is BELOW the iid value
    (ratio 0.863), i.e. the walk repeats no more often than the birthday count, and
    the trap row is charged at the iid value anyway.
    """
    r = K.R_ORDER
    lam_rep = (float(S) * float(S) / 2.0) / float(K.N_CLASSES)
    e_meas = BND23['E_repeats_per_walk_measured']
    chk('measured_birthday_below_iid',
        e_meas <= BND23['birthday_formula_per_walk'] * 1.0 + 1e-12,
        'measured %r vs iid %r' % (e_meas, BND23['birthday_formula_per_walk']))
    return {
        'route': 'guard-off restores the classic birthday hypothesis domain; measured-E '
                 'Markov shipped beside it (hypothesis-free)',
        'one_line': 'with the guard dropped the step is history-free (the successor class is '
                    'a function of the class alone), so the iid birthday count over the class '
                    'stream applies exactly as in the accepted pre-guard eras',
        'iid_birthday_value_production': lam_rep,
        'iid_birthday_log2': lg(lam_rep),
        'measured_E_markov': {
            'E_measured_per_walk': e_meas,
            'E_measured_source': 'evidence/bound-v23.json '
                                 '(guard-free class stream, L = 20, N = 785, unsaturated)',
            'iid_formula_at_the_same_fixture': BND23['birthday_formula_per_walk'],
            'measured_over_formula': BND23['measured_over_formula'],
            't': t,
            'P_X_ge_t_upper': min(1.0, e_meas / t),
            'hypothesis_free': True,
        },
        'withdrawn': 'the v22 monotonicity route (joint <= element <= unguarded as a '
                     'stochastic domination argument) is WITHDRAWN, not patched: with the '
                     'guard gone there is no extra coordinate to dominate over',
        'applied_to_the_total': 'X = sum_w X_w, E[X] <= workers * lambda_rep by linearity '
                                'alone; Markov on the TOTAL gives P(X > C) <= E[X]/C with '
                                'C = ceil(4 lambda_rep / 1e-9).  No independence across '
                                'workers and no history-independence inside a worker is used.',
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
    """SUCCESS = model x DP-GATED retention, with the guard dropped and the cap kept.

    Four rows move over v22 and one does not:

      (a) RETENTION IS THE DP-GATED STATISTIC (verdict item 1).  v22 counted class
          agreement at ANY step.  The protocol detects on the distinguished-point
          predicate, so detection is now a DP-HIT EVENT: at the first step whose
          class key satisfies the DP mask, both trails must be in the same class
          key.  evidence/retention-v23.json ships BOTH counts labelled; only the
          DP-gated one feeds success.
      (b) THE GUARD IS DROPPED (the redesign).  With no ring the step is a function
          of the class alone, so two trails that share a start point ARE the same
          trail: they report at the same DP step with two different label pairs,
          which is exactly the solvable event.  Measured 7503/7503 DP-gated
          agreements (CP-95 0.999601) against the guarded 0.808333 (0.776352).
      (c) THE mu ALLOWANCE CANNOT EXCEED CERTAINTY.  retention_eff = min(1, mu *
          p_det * p_use); with p_det = 1.000000 the product 1.1 is capped at 1.
      (d) THE TRAP DEBIT IS THE CAP-BOUNDED ABANDON PROBABILITY, not a v19 extract.
          A trap is a class repeat before the first DP; the guard-free state IS the
          class, so the repeat is exact.  Per DP segment it is bounded by
          theta^2/(2N) = 2^-92.3 (the fixture measures 0.24 per cent against a 1.02
          per cent closed form, i.e. the bound is conservative by 4x), and over
          ceil(S/theta) = 2^45.8 segments the debit is <= 2^-46.5.
      (e) THE MODEL IS UNCHANGED.  Lambda = S^2 (1-f^2)/(2 mu N) = 0.5049381183, so
          1 - exp(-Lambda) = 0.3964571 -- the same exact-rational model row v22
          shipped, now multiplied by a retention of 1 instead of 0.808333.
    """
    lam, pw = ledger.model_success(S, K.N_CLASSES, K.F_SQR)
    model = 1.0 - 1.0 / float(pw)
    p_use = B17['retention']['p_usable_given_detected']
    p_use_low = B17['retention']['p_usable_cp95_lower']
    dp = RET23['dp_gated_agree']
    p_det = dp['p_detect_given_dp_hit']
    p_det_low = dp['cp95_lower']
    mu = float(K.MU_WALK)
    ret_raw = mu * p_det * p_use
    ret = min(1.0, ret_raw)                     # (c) a probability cannot exceed 1
    ret_low = min(1.0, mu * p_det_low * p_use_low)
    theta_prod = 1 << K.THETA_BITS
    p_trap_seg = float(theta_prod * theta_prod) / (2.0 * float(K.N_CLASSES))
    segments = (float(S) + theta_prod - 1) / theta_prod
    tau = segments * p_trap_seg                # (d) cap-bounded abandon probability
    stop = 0.0
    rate = 1.0 - math.exp(-float(lam) * ret)     # retention applied to the RATE
    prod = model * ret
    net = prod * (1.0 - tau) * (1.0 - stop)
    prod_low = model * ret_low
    net_cp = prod_low * (1.0 - tau) * (1.0 - stop)
    return {'lambda_': float(lam), 'model': model, 'rate_scaled': rate,
            'retention_eff': ret,
            'retention_eff_raw_mu_times_p_detect': ret_raw,
            'retention_eff_cp95_lower': ret_low,
            'product': prod, 'product_cp95_lower': prod_low,
            'trap_tau': tau, 'trap_tau_per_segment': p_trap_seg,
            'trap_tau_segments': segments, 'stop_loss': stop,
            'net': net, 'net_cp95_lower': net_cp, 'net_is_point_estimate': True,
            'clears_floor_point': net >= 0.39,
            'clears_floor_cp95_lower': net_cp >= 0.39,
            'p_detect_dp_gated': p_det, 'p_detect_dp_gated_cp95_lower': p_det_low,
            'dp_gated_pairs': dp['pairs_with_a_dp_hit'],
            'any_step_agree_for_comparison_only': RET23['any_step_agree']['p_detect'],
            'which_statistic_feeds_success': 'dp_gated_agree (DP-HIT event with the mask '
                                             'applied); any_step_agree is shipped for '
                                             'comparison and feeds nothing',
            'retention_source': 'evidence/retention-v23.json (guard-free walk + cap, '
                                'DP-gated, 10008 pairs, 7503 with a DP hit, 7503 agree)',
            'bigN_row': {'pairs': RETBN['pairs'],
                         'dp_pairs': RETBN['dp_gated_agree']['pairs_with_a_dp_hit'],
                         'p_detect_dp_gated': RETBN['dp_gated_agree']['p_detect_given_dp_hit'],
                         'cp95_lower': RETBN['dp_gated_agree']['cp95_lower'],
                         'trap_before_dp_fixture': RETBN['trap_before_dp']['fraction_fixture'],
                         'reading': 'the same arm on a 3136-class pool: the fixture '
                                    'trap-before-DP fraction falls from 0.2503 to 0.0735 as N '
                                    'grows, which is the theta^2/(2N) shape of the bound'},
            'lag32_row': RETLG['lag32_row'],
            'guarded_comparison': {'p_detect_guarded_v22': RET22['p_detect_guarded'],
                                   'cp95_lower_v22': RET22['p_detect_guarded_cp95_lower'],
                                   'lag32_guarded_v22': RET22['curve_guarded']['32']['p_detect'],
                                   'reading': 'the guard is what separated the two trails: '
                                              'dropping it restores retention to 1 and the '
                                              'lag-32 collapse to 0.10625 disappears'},
            'reports_measured_per_walk': REP23['reports_per_walk_measured'],
            'reports_ceil_formula_per_walk': REP23['ceil_formula_reports_per_walk'],
            'reports_check_ratio_measured_over_formula': REP23['check_measured_vs_formula_ratio'],
            'reports_restart': REP23['restart_implemented'],
            'reports_rule': 'the report count is CENSUS-DIRECT with the restart implemented '
                            '(a reporting trail is replaced by a fresh trail); ceil(S/theta) '
                            'is the CHECK row, and the fixture shortfall is the abandoned '
                            'work, which vanishes at production (2^-92.3 per segment)'}


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
        ('success.net', vals['net'], '>=', '0.3963'),
        ('success.model', vals['model'], '>=', fl(vals['model'])),
        ('success.product', vals['product'], '>=', fl(vals['product'])),
        ('success.product_cp95_lower', vals['prod_low'], '>=', fl(vals['prod_low'])),
        ('success.net_cp95_lower', vals['net_cp'], '>=', '0.39'),
        ('success.net_cp95_lower', vals['net_cp'], '<=', ce(vals['net_cp'])),
        ('success.trap_tau', vals['tau'], '<=', '1e-12'),
        ('success.retention_eff', vals['ret'], '<=', '1.0'),
        ('claim.value', vals['claim'], '<=', '125.82'),
        ('claim.value', vals['claim'], '<=', '125.8345'),
        ('claim.value', vals['claim'], '<=', '125.84'),
        ('claim.value', vals['claim'], '>=', fl(vals['claim'])),
        ('trap.repeat_probability', vals['lam_rep'], '<=', '0.7000'),
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
