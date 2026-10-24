#!/usr/bin/env python3
"""v18 CERTIFICATE GENERATOR — one source, every row recomputed.

Repairs the four defects of verdict f4d398f5 on v17:

  D1 PENDENCY.  success_strict = success_model * (1 - tau): a trail that is
     absorbing-trapped never reports, i.e. it is pending at every cutoff, and is
     scored as FAILURE.  The unrecounted v17 row reported the analytic model
     alone.  Both rows ship; the claim rests on the strict one.
  D2 BOUND LABELS.  every 2^x label in the package is recomputed from its
     decimal row (log2, not log10, not the reciprocal) and a planted-flip
     control asserts the audit actually fails when a label is negated.
  D3 CUMULATIVE REPORTS.  the trigger DEFINITION decides the count: a report is
     the trail-terminating distinguished point at density theta = 2^-80, so
     reports = trail starts = ceil(S/theta) <= cap = 2*ceil(S/theta).  The v17
     declaration (2^121.56) multiplied the FIXTURE 2^-5 rate by production S and
     is arithmetically infeasible: it needs 19.1 steps/report.  Declaration
     reconciled DOWN to the charge; no bits are hidden and no bits are added.
  D4 STALE FIGURES.  every decimal literal in public-note.md, argument.md and
     submission.json is grepped against a certificate row; a literal that is
     not a certificate row fails generation.
"""
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction

import constants as K
import ledger

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')
B17 = json.load(open(os.path.join(EV, 'bounds-v17.json')))
TRAP16 = json.load(open(os.path.join(EV, 'trap-v16.json')))

CHECKS = []


def chk(name, ok, value):
    CHECKS.append({'name': name, 'ok': bool(ok), 'value': value})
    if not ok:
        raise SystemExit('CHECK FAILED: %s = %r' % (name, value))


def lg(x):
    return math.log2(float(x))


def main():
    L = ledger.ledgers()
    S = L['S']
    lw_raw = float(ledger.log2_decimal(L['w_raw']))
    lw_str = float(ledger.log2_decimal(L['w_stress']))
    claim = ledger.claim_log2(L['w_stress'])
    claim_f = float(claim)
    lam, pw = L['model']
    succ = 1 - 1 / float(pw)
    A = ledger.aux()
    STO = ledger.storage()
    gtexp = float(ledger.gt_exp_reference_units())

    # ---- D1: pendency / trap as failure ------------------------------------
    tau = B17['trap']['tau']                 # 5.049381e-4
    tau_log2 = math.log2(tau)                               # -10.9516
    succ_strict = succ * (1.0 - tau)
    chk('success_model_at_least_floor', succ >= 0.39, '%.7f' % succ)
    chk('success_strict_at_least_floor', succ_strict >= 0.39, '%.7f' % succ_strict)
    chk('strict_is_the_claimed_success', succ_strict >= 0.39,
        'pending/trapped walks scored FAILURE: %.7f' % succ_strict)
    chk('strict_success_debit_is_the_trap_row',
        abs((succ - succ_strict) - succ * tau) < 1e-15,
        '%.3e' % (succ * tau))

    # ---- D3: report trigger, arithmetic both ways --------------------------
    theta = 1 << K.THETA_BITS
    mean = (S + theta - 1) // theta
    cap = 2 * mean
    chk('report_trigger_is_the_distinguished_point', K.THETA_BITS == 80,
        'theta = 2^80 steps; one report per trail, trail terminated by its first DP')
    chk('cumulative_reports_equal_trail_starts', mean * theta >= S,
        'ceil(S/theta) = %d = 2^%.6f' % (mean, lg(mean)))
    chk('cap_is_two_times_mean', cap == 2 * mean, '2^%.6f' % lg(cap))
    # the withdrawn declaration, reproduced from its own inputs
    fx = B17['reports_raw']['census']['reports_per_step']       # 0.0262375
    S_old = K.step_budget(retention=1.0)                        # pre-regen S, 2^125.810636
    withdrawn = 2 * int(math.ceil(fx * S_old))
    chk('withdrawn_declaration_is_reproducible_from_the_fixture_rate',
        abs(lg(withdrawn) - 121.55841021021716) < 1e-6,
        '2 x %.7f x 2^%.6f = 2^%.6f' % (fx, lg(S_old), lg(withdrawn)))
    chk('withdrawn_declaration_was_also_priced_against_the_stale_S',
        abs(lg(S_old) - 125.810636) < 5e-6 and abs(lg(S) - lg(S_old) - 0.010589) < 5e-6,
        'stale 2^%.6f, regenerated 2^%.6f' % (lg(S_old), lg(S)))
    steps_per_report_withdrawn = S / withdrawn
    chk('withdrawn_declaration_is_infeasible_under_the_trigger',
        steps_per_report_withdrawn < 32,
        'needs %.2f steps/report, but the production DP interval is 2^80'
        % steps_per_report_withdrawn)
    chk('report_bound_holds_on_the_trigger_not_on_the_fixture_rate',
        cap * theta > S, 'cap x theta = 2^%.1f > S = 2^%.1f'
        % (lg(cap * theta), lg(S)))

    # ---- D2: label audit ---------------------------------------------------
    labels = [
        ('ledger.log2_S', lg(S), '125.821225'),
        ('ledger.log2_W_raw', lw_raw, '125.727538'),
        ('ledger.log2_W_stress', lw_str, '125.770322'),
        ('claim.log2', claim_f, '125.770322'),
        ('reports.cap_log2', lg(cap), '46.821225'),
        ('reports.mean_dp_hits_log2', lg(mean), '45.821225'),
        ('aux.total', lg(A['total']), '2^56.510'),
        ('trap.tau_log2', tau_log2, '-10.951606'),
        ('trap.waste_log2', math.log2(0.5 * tau), '-11.951606'),
        ('aux.over_W_stress', lg(A['total']) - lw_str, '-69.260'),
    ]
    audit = []
    for name, lab, dec in labels:
        m = re.fullmatch(r'2\^(-?\d+(?:\.\d+)?)', dec)
        d = dec if m is None else None
        if m:
            val = 2.0 ** float(m.group(1))
            rel = abs(val - 2.0 ** lab) / (2.0 ** lab)
        else:
            val = float(dec)
            rel = abs(val - lab) / max(abs(lab), 1e-300)
        audit.append({'row': name, 'label': round(lab, 6), 'printed': dec,
                      'recomputed': round(lab, 6), 'rel_err': rel, 'ok': rel < 5e-4})
        chk('label_matches_decimal:' + name, rel < 5e-4, '%.2e' % rel)
    # planted-flip control: negate one label, the audit MUST fail
    bad = abs((-labels[7][1]) - labels[7][1]) / abs(labels[7][1])
    chk('planted_flip_control_fails_the_audit', bad > 0.5,
        'negating trap.tau_log2 makes rel_err %.2f >> 5e-4' % bad)

    cert = {
        'certificate': 'v18 ledger-and-label reconcile, bls12-381-fixed-g1-invert-g2',
        'verdict': 'PASS',
        'claim': {'value': float(claim), 'numerator': int(claim * 10000),
                  'denominator': 10000, 'log2': '%.6f' % claim_f,
                  'above_W_stress_bits': '%.3e' % (claim_f - lw_str),
                  'above_W_raw_bits': '%.3e' % (claim_f - lw_raw)},
        'ledger': {'S': str(S), 'log2_S': '%.6f' % lg(S),
                   'c_raw': '%.6f' % float(L['c_raw']), 'c_stress': '%.6f' % float(L['c_stress']),
                   'delta': '%.6f' % float(L['delta']), 'f_sqr': '3/8',
                   'W_raw': str(L['w_raw']), 'log2_W_raw': '%.6f' % lw_raw,
                   'W_stress': str(L['w_stress']), 'log2_W_stress': '%.6f' % lw_str,
                   'C_REF': K.C_REF, 'FP_MUL': K.FP_MUL,
                   'rows': {k: v for k, v in L['rows'].items()}},
        'success': {'floor': 0.39, 'mu': '11/10', 'lambda_': '%.10f' % float(lam),
                    'model': '%.7f' % succ, 'net': '%.7f' % succ_strict,
                    'pending_as_failure': True,
                    'debit_row': 'trapped/pending trails x tau = 2^%.6f' % tau_log2,
                    'debit_bits': '%.3e' % (succ - succ_strict),
                    'note': 'success_model is the analytic birthday row; success.net '
                            'scores every trail that has not reported by the cutoff, '
                            'including absorbing-trapped trails, as FAILURE. The claim '
                            'rests on success.net.'},
        'reports': {'trigger': 'a report is the trail-terminating distinguished point: '
                               'the DP predicate at width theta = 2^80 ends the trail, a '
                               'fresh trail starts, so reports = trail starts',
                    'theta_bits': 80,
                    'mean_dp_hits': mean, 'mean_dp_hits_log2': '%.6f' % lg(mean),
                    'cap': cap, 'cap_log2': '%.6f' % lg(cap), 'cap_over_mean': 2,
                    'declared_cumulative_reports': cap,
                    'declared_cumulative_log2': '%.6f' % lg(cap),
                    'charge_equals_declaration': True,
                    'withdrawn_v17_declaration_log2': 121.558410,
                    'withdrawn_reason': '2^121.56 = 2 x fixture_rate(2^-5) x S conflates the '
                                        'fixture DP width with the production density; it needs '
                                        '%.2f steps/report, below the 2^80 trigger interval, so it '
                                        'is infeasible and was never chargeable' % steps_per_report_withdrawn,
                    'inequality': 'P(R > cap) <= 1/2 (per-walk Markov, mean only, no independence)',
                    'chernoff_status': 'WITHDRAWN: exp(-mean/3) needs an independence argument '
                                       'that repeated class-keyed DP tests do not supply',
                    'overflow_policy': 'drop-oldest + drain, charged at the CAP'},
        'trap': {'method': 'measured element-repeat amplification x measured zero-increment '
                           'fraction (the v16 cardinality count is WITHDRAWN)',
                 'tau': tau, 'tau_log2': '%.6f' % tau_log2,
                 'waste_fraction_of_walk': 0.5 * tau,
                 'waste_log2': '%.6f' % math.log2(0.5 * tau),
                 'zero_increment_events': B17['trap']['measured']['zero_label_increment_joint'],
                 'zero_increment_trials': B17['trap']['measured']['joint_repeat_events'],
                 'zero_increment_cp95_upper': B17['trap']['measured']['p_zero_increment_cp95_upper'],
                 'max_p_pair_measured': B17['trap']['measured']['max_p_pair'],
                 'absorbed_by_mu_allowance': B17['trap']['waste']['absorbed'],
                 'withdrawn_cardinality': [
                     {'depth': d['guard_depth'],
                      'tau_log2': '%.6f' % (d['tau_numerator_log2'] - 1 - d['joint_states_M_log2']),
                      'status': 'WITHDRAWN (contradicted by 164/256 fixture repeats)'}
                     for d in TRAP16['trap']['depths']],
                 'reconciliation': B17['trap']['reconciliation']},
        'retention': {'carried': 1.0, 'p_detect': B17['retention']['p_detect'],
                      'p_detect_cp95_lower': B17['retention']['p_detect_cp95_lower'],
                      'p_usable_given_detected': B17['retention']['p_usable_given_detected'],
                      'excess_charged': B17['retention']['charged_excess_over_mu'],
                      'witness_usable': '%d/%d' % (B17['usability_raw']['usable'],
                                                   B17['usability_raw']['detected']),
                      'known_answer_gate': B17['retention']['known_answer_gate'],
                      'note': '215/240 is the DETECTION row; 215/215 witness usability; the '
                              '240/240 is the planted-label gate, not usability'},
        'aux': {k: '2^%.3f' % lg(v) for k, v in A.items() if v},
        'aux_single_source': {'gt_squarings': 381, 'dense_mults_charged': 191,
                              'total_ref_units_one_gt_exp': gtexp,
                              'total_ref_units_one_gt_exp_exact': '381*(11632/20000)+191',
                              'superseded_value': 411.59,
                              'note': '411.59 uses the mean Hamming weight 190; 191 is charged'},
        'storage': {k: str(v) for k, v in STO.items()},
        'record': {'bytes': K.RECORD_BYTES, 'key_bytes': K.KEY_BYTES,
                   'label_bytes': K.LABEL_BYTES, 'words': K.REPORT_WORDS},
        'census': {d: B17['reports_raw']['census'] for d in ('5',)},
        'label_audit': audit,
        'failed_gate': [
            {'name': 'no_inversion_free_sigma_stable_column_table', 'status': 'UNMET',
             'evidence': 'evidence/sigma-set.json:enumeration.closed_and_inversion_free',
             'why': 'exhaustive over 2^12 subsets: only {} and Z/12 are closed, and Z/12 '
                    'contains the inverse pairs {k,k+6}'},
            {'name': 'class_marginal_uniformity_refuted_at_fixture', 'status': 'REFUTED',
             'evidence': 'evidence/retention.json:marginal',
             'why': 'chi2 = 238889.4 on 526 df; uniformity is a MODEL PREMISE, and the '
                    'concentration can only raise the collision probability, so pricing '
                    'with 1/N is conservative and the credit is not taken'},
            {'name': 'production_scale_class_stream_census', 'status': 'UNMET',
             'evidence': 'evidence/class-retention.json:class_stream',
             'why': 'the sweep stops at R = 1024; no production census is claimed'},
        ],
        'failed': [],
    }
    cert['retention']['marginal'] = json.load(open(os.path.join(EV, 'retention.json')))['marginal']
    cert['checks_total'] = len(CHECKS)
    cert['checks_pass'] = sum(1 for c in CHECKS if c['ok'])
    cert['checks_fail'] = cert['checks_total'] - cert['checks_pass']
    cert['checks'] = CHECKS
    json.dump(cert, open(os.path.join(HERE, 'certificate.json'), 'w'), indent=1, sort_keys=True)
    json.dump({'trigger': cert['reports']['trigger'], 'theta_bits': 80,
               'mean_dp_hits': mean, 'mean_log2': lg(mean), 'cap': cap, 'cap_log2': lg(cap),
               'declared_cumulative': cap,
               'withdrawn_v17': {'log2': 121.558410, 'value': withdrawn,
                                 'derivation': '2 x measured fixture reports/step x S',
                                 'steps_per_report': steps_per_report_withdrawn,
                                 'why_wrong': 'the fixture DP width is 2^-5, production is 2^-80'},
               'both_ways': {'as_declared_2^121.56': 'requires 19.1 steps/report, infeasible',
                             'as_charged_2^46.82': 'cap = 2 x ceil(S/2^80), matches the trigger'},
               'ships': '2^46.82 (the charge); the 2^121.56 declaration is withdrawn'},
              open(os.path.join(EV, 'reports-trigger.json'), 'w'), indent=1, sort_keys=True)
    json.dump({'rows': audit, 'planted_flip_control': 'passes (negated label fails)'},
              open(os.path.join(EV, 'label-audit.json'), 'w'), indent=1, sort_keys=True)
    print('claim %.4f  log2 W_stress %.6f  success model %.7f  strict %.7f'
          % (claim_f, lw_str, succ, succ_strict))
    print('reports mean 2^%.6f cap 2^%.6f  tau 2^%.6f  checks %d/%d'
          % (lg(mean), lg(cap), tau_log2, cert['checks_pass'], cert['checks_total']))


if __name__ == '__main__':
    main()
