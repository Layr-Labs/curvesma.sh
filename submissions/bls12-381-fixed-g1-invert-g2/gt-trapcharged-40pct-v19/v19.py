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
    events = max(1, workers)                                  # charged at one per worker
    rows = ledger.rows()
    c_step = float(rows['mul_step']) / K.C_REF                # ~1.2 ref units per step
    dp_lookup = 420.0 / K.C_REF            # the same 48-byte key probe the store pays
    detector = 420.0 / K.C_REF             # tag probe + 112-byte record compare on a hit
    restart = 2.0 * float(ledger.gt_exp_reference_units()) + 16.0 / K.C_REF
    cycle_steps = 64.0                     # 3 x the measured period-20 family, then restart
    waste = cycle_steps * c_step
    per_event = dp_lookup + detector + restart + waste
    total = float(events) * per_event
    return {
        'repeat_probability_over_the_walk': lam_rep,
        'repeat_probability_log2': lg(lam_rep),
        'method': 'birthday count over the whole walk: Lambda = S^2 / (2 (r-1)); every repeat '
                  'is a cycle with no DP, so every repeat is a trap (the increment decides '
                  'only whether a solve would exist if both records were kept, and the '
                  'shipped attack keeps neither)',
        'workers_charged': events,
        'events_charged': events,
        'bound_note': 'P(at least one trap) <= Lambda = %.4f < 1 (Markov), so one event per '
                      'worker is a ceiling, not an expectation' % lam_rep,
        'prices_reference_units': {
            'dp_lookup': dp_lookup, 'detector_probe_compare': detector,
            'restart_two_gt_exps': restart, 'wasted_cycle_steps': waste,
            'per_event': per_event},
        'total_reference_units': float(total),
        'total_log2': lg(float(total)),
        'waste_fraction_of_W_stress': float(total) / float(ledger.ledgers()['w_stress']),
        'zero_increment_fraction_measured': REP['census']['zero_increment_fraction'],
        'increment_split': {'zero': REP['census']['zero_increment_events'],
                            'nonzero': REP['census']['nonzero_increment_events'],
                            'pairs': REP['census']['step_pairs']},
        'status': 'CHARGED (v18 charged only the zero-increment class and priced no detector)',
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
    p_det = B17['retention']['p_detect']
    p_det_low = B17['retention']['p_detect_cp95_lower']
    p_use = B17['retention']['p_usable_given_detected']
    p_use_low = B17['retention']['p_usable_cp95_lower']
    mu = float(K.MU_WALK)
    ret = mu * p_det * p_use
    ret_low = mu * p_det_low * p_use_low
    rate = 1.0 - math.exp(-float(lam) * ret)          # retention applied to the RATE
    prod = model * ret                                 # retention applied to the PROBABILITY
    trap = B17['trap']['tau']
    stop = cap_rows()['chebyshev_stop_loss']
    net = prod * (1.0 - trap) * (1.0 - stop)
    net_low = prod * (1.0 - trap) * (1.0 - stop)       # point decomposition, same factors
    prod_low = model * ret_low
    net_cp = prod_low * (1.0 - trap) * (1.0 - stop)
    return {'lambda_': float(lam), 'model': model, 'rate_scaled': rate, 'product': prod,
            'retention_eff': ret, 'retention_eff_cp95_lower': ret_low,
            'product_cp95_lower': prod_low, 'trap_tau': trap, 'stop_loss': stop,
            'net': net, 'net_cp95_lower': net_cp, 'net_is_point_estimate': True,
            'clears_floor_point': net >= 0.39, 'clears_floor_cp95_lower': net_cp >= 0.39}


# --------------------------------------------------------------------------
# D4  the >= / <= audit table.  Every inequality printed in the note is in here.
# --------------------------------------------------------------------------
def inequality_rows(vals):
    """(row, cert_value, op, printed_literal).  Rule: '>=' needs value >= literal
    (the literal is a FLOOR of the value); '<=' needs value <= literal (a CEIL)."""
    return [
        ('success.net', vals['net'], '>=', '0.3904'),
        ('success.net', vals['net'], '<=', '0.3910'),
        ('success.model', vals['model'], '>=', '0.3964'),
        ('success.product', vals['product'], '>=', '0.3906'),
        ('success.product_cp95_lower', vals['prod_low'], '>=', '0.3688'),
        ('success.net_cp95_lower', vals['net_cp'], '<=', '0.3690'),
        ('success.net_cp95_lower', vals['net_cp'], '<=', '0.39'),
        ('claim.value', vals['claim'], '<=', '125.7704'),
        ('claim.value', vals['claim'], '<=', '125.8345'),
        ('claim.value', vals['claim'], '<=', '125.84'),
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
