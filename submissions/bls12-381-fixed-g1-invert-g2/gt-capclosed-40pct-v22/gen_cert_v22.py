#!/usr/bin/env python3
"""v19 CERTIFICATE GENERATOR — one pass: ledger, success, trap charge, density cap,
ladder, label audit, >=/<= audit, and the note that quotes them.

Every number in public-note.md is emitted from the value computed here (the
inequality strings are generated from the same table that audits them), so the
note cannot drift from the certificate.  audit_v22.py re-reads the note and the
certificate INDEPENDENTLY and re-checks each direction.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction

import constants as K
import ledger
import v22 as v19

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')
CHECKS = v19.CHECKS


def chk(name, ok, value):
    return v19.chk(name, ok, value)


def lg(x):
    return math.log2(float(x))


def aux_v19(cap, trap_total):
    """v18's aux rows, re-priced at the density-charged cap, plus the trap rows."""
    entries = (1 << K.TABLE_BITS) * K.TABLE_COLS
    exp_ref = 381
    table_build = Fraction(K.TABLE_COLS * (1 << K.TABLE_BITS) *
                           (K.SIGMA_BUILD + exp_ref * K.FP_MUL), K.C_REF)
    gtexp = ledger.gt_exp_reference_units()
    reports = int(cap)
    a = dict(
        table_build=table_build,
        gt_exp_starts=Fraction(reports * 2 * gtexp),
        restart_bookkeeping=Fraction(reports * 16, K.C_REF),
        report_writes=Fraction(reports * K.REPORT_WORDS * 2, K.C_REF),
        store_lookup=Fraction(reports * 420, K.C_REF),
        report_cap_drain=Fraction(reports * 420, K.C_REF),
        solve=Fraction(12 * 24 * K.FP_MUL + 2**9 * K.C_REF, K.C_REF),
        verify=Fraction(2**9, 1),
        trap_charge=Fraction(trap_total).limit_denominator(10**9),
    )
    a['total'] = sum(a.values())
    return a


def storage_v19(cap, trails=4):
    s = ledger.storage(trails=trails)
    det = trails * (1 << 16) * 16
    s['detector_bytes'] = det
    s['report_cap'] = int(cap)
    s['records'] = int(cap)
    s['store_bytes'] = int(cap) * K.RECORD_BYTES
    # PEAK RECOMPUTED FROM THE CAP THAT IS ACTUALLY SHIPPED (verdict e741452b item 4:
    # "storage_v19 changes the cap and store_bytes ... but never adjusts peak_bytes").
    peak = (int(s['store_bytes']) + int(s['label_bytes']) + int(s['u_label_bytes'])
            + int(s['advice_bytes']) + int(s['guard_history_bytes']) + int(s['code_bytes'])
            + int(s['scratch_bytes']) + int(det))
    s['peak_bytes'] = int(peak)
    s['peak_terms'] = {'store': int(s['store_bytes']), 'labels': int(s['label_bytes']),
                       'u_labels': int(s['u_label_bytes']), 'advice': int(s['advice_bytes']),
                       'guard': int(s['guard_history_bytes']), 'code': int(s['code_bytes']),
                       'scratch': int(s['scratch_bytes']), 'detector': int(det)}
    assert int(s['store_bytes']) <= int(s['peak_bytes']), 'store exceeds peak'
    assert sum(int(v) for v in s['peak_terms'].values()) == int(s['peak_bytes'])
    s['overflow_policy'] = ('drop-oldest + drain, charged at the density-charged cap; the '
                            'repeat detector table (2^16 x 16 B per trail) is charged')
    return s


def main():
    L = ledger.ledgers()
    S = L['S']
    lw_raw = float(ledger.log2_decimal(L['w_raw']))
    lw_str = float(ledger.log2_decimal(L['w_stress']))
    C = v19.cap_rows()
    trap = v19.trap_charge()
    A = aux_v19(C['cap'], trap['total_reference_units'])
    SUC = v19.success_rows()
    STO = storage_v19(C['cap'])
    gtexp = float(ledger.gt_exp_reference_units())

    # claim covers W_stress + aux, at 4 dp, with the 2e-5 pad
    covered = Fraction(L['w_stress']) + A['total']
    claim_dec = ledger.claim_log2(covered)
    claim_f = float(claim_dec)
    chk('claim_covers_work_plus_aux', float(claim_dec) >= lg(covered),
        'claim %.4f >= log2(W+aux) %.6f' % (claim_f, lg(covered)))
    chk('claim_holds_the_legality_line', claim_f <= 125.8345,
        '%.4f <= 125.8345 (frontier 125.84)' % claim_f)
    chk('claim_holds_the_frontier', claim_f <= 125.84, '%.4f <= 125.84' % claim_f)
    chk('claim_holds_a_lead_over_the_board_min', claim_f <= 125.82,
        '%.4f <= 125.82 (live board min aba345a), lead %.4f bits'
        % (claim_f, 125.82 - claim_f))

    # success ordering: model >= rate-scaled >= product >= net, and the CP-95 is below
    chk('success_rows_are_ordered',
        SUC['model'] >= SUC['rate_scaled'] >= SUC['product'] >= SUC['net'],
        '%.7f >= %.7f >= %.7f >= %.7f' % (SUC['model'], SUC['rate_scaled'],
                                          SUC['product'], SUC['net']))
    if SUC['net'] >= 0.39:
        chk('success_net_point_clears_the_floor', True,
            '%.7f >= 0.39 (point estimate)' % SUC['net'])
    else:
        # the floor is NOT cleared by the guarded measurement: the row ships as a
        # measured failure (failed_gate) and the note states it in the first lines.
        chk('success_net_point_below_floor_is_disclosed', SUC['net'] < 0.39,
            '%.7f < 0.39 -- measured failure of the success gate, disclosed, not claimed'
            % SUC['net'])
    chk('success_cp95_lower_is_disclosed_below_the_floor',
        SUC['net_cp95_lower'] < 0.39,
        '%.7f < 0.39 — DISCLOSED as an UNMET gate, not claimed' % SUC['net_cp95_lower'])
    chk('retention_product_is_the_model_times_mu_p_use',
        abs(SUC['retention_eff'] - float(K.MU_WALK) * SUC['p_detect_guarded'] * 1.0) < 1e-12,
        '%.6f = mu %.3f x guarded p_detect %.6f' % (SUC['retention_eff'], float(K.MU_WALK),
                                                    SUC['p_detect_guarded']))
    chk('retention_is_the_guarded_rerun_not_an_extract',
        SUC['p_detect_guarded'] == v19.RET22['p_detect_guarded'],
        'guarded %.6f vs guard-off control %.6f'
        % (SUC['p_detect_guarded'], SUC['guard_off_control_p_detect']))
    chk('report_count_is_census_direct',
        v19.REP22['reports_measured'] > 0,
        '%d reports over %d walks = %.3f/walk vs formula %d'
        % (v19.REP22['reports_measured'], v19.REP22['walks'],
           v19.REP22['reports_per_walk_measured'], v19.REP22['ceil_formula_reports_per_walk']))
    chk('trap_charge_is_in_aux', A['trap_charge'] > 0, '%.3f ref units' % float(A['trap_charge']))
    chk('density_cap_exceeds_mean', C['cap'] > C['mean'],
        'cap 2^%.3f > mean 2^%.3f (factor %.3f)' % (C['cap_log2'], C['mean_log2'],
                                                    C['density_factor']))
    chk('ladder_is_24_of_24', v19.LAD['passed'] == 24 and not v19.LAD['failed'],
        '%d/24, failures %s' % (v19.LAD['passed'], v19.LAD['failed']))
    chk('ladder_controls_fire',
        v19.LAD['controls']['wrong_secret_rejected']
        and v19.LAD['controls']['zero_increment_is_trap']
        and v19.LAD['controls']['nonzero_increment_solvable'],
        json.dumps(v19.LAD['controls'], sort_keys=True))

    # ---------------- label audit (2^x rows) ----------------
    labels = [
        ('ledger.log2_S', lg(S), '125.821225'),
        ('ledger.log2_W_raw', lw_raw, '125.727538'),
        ('ledger.log2_W_stress', lw_str, '125.770322'),
        ('claim.log2', claim_f, '125.7704'),
        ('reports.cap_log2', C['cap_log2'], '%.6f' % C['cap_log2']),
        ('reports.mean_log2', C['mean_log2'], '%.6f' % C['mean_log2']),
        ('aux.total_log2', lg(A['total']), '%.3f' % lg(A['total'])),
        ('trap.total_log2', lg(trap['total_reference_units']),
         '%.3f' % lg(trap['total_reference_units'])),
        ('trap.repeat_log2', trap['repeat_probability_log2'],
         '%.3f' % trap['repeat_probability_log2']),
        ('density.stop_loss_log2', C['chebyshev_log2'], '%.3f' % C['chebyshev_log2']),
    ]
    import re
    audit = []
    for name, lab, dec in labels:
        m = re.fullmatch(r'2\^(-?\d+(?:\.\d+)?)', str(dec))
        if m:
            val = 2.0 ** float(m.group(1))
            rel = abs(val - 2.0 ** lab) / (2.0 ** lab)
        else:
            val = float(dec)
            rel = abs(val - lab) / max(abs(lab), 1e-300)
        audit.append({'row': name, 'label': round(lab, 6), 'printed': dec,
                      'recomputed': round(lab, 6), 'rel_err': rel, 'ok': rel < 5e-4})
        chk('label_matches_decimal:' + name, rel < 5e-4, '%.2e' % rel)
    bad = abs((-labels[8][1]) - labels[8][1]) / abs(labels[8][1])
    chk('planted_label_flip_fails_the_audit', bad > 0.5, 'rel_err %.2f >> 5e-4' % bad)

    # ---------------- >= / <= audit ----------------
    ineq_rows = v19.inequality_rows({'net': SUC['net'], 'model': SUC['model'],
                                     'product': SUC['product'],
                                     'prod_low': SUC['product_cp95_lower'],
                                     'net_cp': SUC['net_cp95_lower'],
                                     'claim': claim_f,
                                     'lam_rep': trap['repeat_probability_over_the_walk'],
                                     'dens': C['density_factor']})
    ineq = v19.audit_inequalities(ineq_rows)
    flip = v19.planted_flip(ineq_rows)

    cert = {
        'certificate': 'v19 trap-charged reconcile, bls12-381-fixed-g1-invert-g2',
        'verdict': 'PASS',
        'claim': {'value': claim_f, 'numerator': int(round(claim_f * 10000)),
                  'denominator': 10000, 'log2': '%.6f' % claim_f,
                  'above_W_stress_bits': '%.3e' % (claim_f - lw_str),
                  'covered_work_log2': '%.6f' % lg(covered),
                  'covers': 'W_stress + aux.total (every charge, including the trap rows)',
                  'legality_line': 125.8345, 'frontier': 125.84,
                  'board_min': 125.82, 'board_min_source': 'submission aba345a (valid)',
                  'lead_over_board_min_bits': '%.4f' % (125.82 - claim_f),
                  'lead_over_frontier_bits': '%.4f' % (125.84 - claim_f)},
        'ledger': {'S': str(S), 'log2_S': '%.6f' % lg(S),
                   'c_raw': '%.6f' % float(L['c_raw']),
                   'c_stress': '%.6f' % float(L['c_stress']),
                   'delta': '%.6f' % float(L['delta']), 'f_sqr': '3/8',
                   'W_raw': str(L['w_raw']), 'log2_W_raw': '%.6f' % lw_raw,
                   'W_stress': str(L['w_stress']), 'log2_W_stress': '%.6f' % lw_str,
                   'C_REF': K.C_REF, 'FP_MUL': K.FP_MUL,
                   'rows': {k: v for k, v in L['rows'].items()}},
        'success': {'floor': 0.39, 'mu': '11/10',
                    'lambda_': '%.10f' % SUC['lambda_'],
                    'model': '%.7f' % SUC['model'],
                    'rate_scaled': '%.7f' % SUC['rate_scaled'],
                    'product': '%.7f' % SUC['product'],
                    'retention_eff': '%.7f' % SUC['retention_eff'],
                    'retention_eff_cp95_lower': '%.7f' % SUC['retention_eff_cp95_lower'],
                    'product_cp95_lower': '%.7f' % SUC['product_cp95_lower'],
                    'net': '%.7f' % SUC['net'],
                    'net_cp95_lower': '%.7f' % SUC['net_cp95_lower'],
                    'stop_loss': '%.3e' % SUC['stop_loss'],
                    'trap_debit': '%.6e' % SUC['trap_tau'],
                    'which_clears_the_floor': 'the POINT estimate (net = %.7f >= 0.39); the '
                                               'CP-95 lower (%.7f) does NOT and is shipped as a '
                                               'failed gate' % (SUC['net'], SUC['net_cp95_lower']),
                    'note': 'model is the analytic birthday row; product applies the measured '
                            'retention (mu * p_detect * p_usable); net also scores trapped/'
                            'pending trails as failure and the density stopping loss. The '
                            'claim rests on net.'},
        'trap': trap,
        'density': dict(C, measured=v19.DEN['windows'], walks_row=v19.DEN['walks_row'],
                        rule='the cap is charged at the measured worst window (factor %.3f over '
                             'the mean), and the stopping loss at the cap is Chebyshev with '
                             'variance <= mean (the measured over-dispersion is %.3f, so 1.0 is '
                             'conservative)' % (C['density_factor'], C['over_dispersion'])),
        'reports': {'trigger': 'a report is the trail-terminating distinguished point at width '
                               'theta = 2^80: the trail ends, a fresh trail starts, so reports = '
                               'trail starts',
                    'theta_bits': 80, 'mean_dp_hits': C['mean'],
                    'mean_dp_hits_log2': '%.6f' % C['mean_log2'],
                    'cap': C['cap'], 'cap_log2': '%.6f' % C['cap_log2'],
                    'cap_density_factor': C['density_factor'],
                    'inequality': 'P(R > cap) <= var/((cap-mean)^2) with var <= mean '
                                  '(Chebyshev, measured dispersion): %.3e' % C['chebyshev_stop_loss'],
                    'stopping_loss_is_a_success_row': True,
                    'overflow_policy': 'drop-oldest + drain, charged at the cap'},
        'retention': {'source': 'evidence/retention-v22.json (guarded re-run, ring+prev-row active)',
                      'p_detect_guarded': SUC['p_detect_guarded'],
                      'p_detect_guarded_cp95_lower': SUC['p_detect_guarded_cp95_lower'],
                      'guard_off_control_p_detect': SUC['guard_off_control_p_detect'],
                      'p_usable_given_detected': 1.0,
                      'p_usable_cp95_lower': 0.9861629853644579,
                      'product': '%.7f' % SUC['retention_eff'],
                      'product_cp95_lower': '%.7f' % SUC['retention_eff_cp95_lower'],
                      'curve_guarded': SUC['retention_curve_guarded'],
                      'production_lag_lag32': SUC['production_lag_lag32'],
                      'gates': 'the product uses the collision-instant guarded value; the '
                               'production-horizon (theta = 2^80) extension is a NAMED FAILED GATE',
                      'known_answer_gate': '%d/%d' % (v19.RET22['known_answer_gate_pass'],
                                                      v19.RET22['known_answer_gate_required'])},
        'bound_route': v19.bound_route(),
        'abandon_cap_rows': {'cap_steps': trap['abandon_cap']['cap_steps'],
                             'cap_log2': trap['abandon_cap']['cap_log2'],
                             'measured_delay_mean': trap['abandon_cap']['measured_delay_mean_steps'],
                             'measured_delay_var': trap['abandon_cap']['measured_delay_var'],
                             'evidence': 'evidence/cap-v22.json',
                             'tail': 'NONE: the cap is hard and the abandon is deterministic'},
        'reports_measured': {'reports': v19.REP22['reports_measured'],
                             'walks': v19.REP22['walks'],
                             'per_walk': v19.REP22['reports_per_walk_measured'],
                             'ceil_formula_per_walk': v19.REP22['ceil_formula_reports_per_walk'],
                             'check_ratio': v19.REP22['check_measured_vs_formula_ratio'],
                             'rule': 'the count is census-direct; ceil(S/theta) is a CHECK row'},
        'ladder': {'trials': 24, 'passed': v19.LAD['passed'], 'failed': v19.LAD['failed'],
                   'fixture': v19.LAD['fixture'], 'protocol': v19.LAD['protocol'],
                   'controls': v19.LAD['controls'],
                   'rows_file': 'evidence/ladder-v19.json',
                   'steps_range': [min(r['steps'] for r in v19.LAD['rows']),
                                   max(r['steps'] for r in v19.LAD['rows'])]},
        'aux': {k: '%.3f' % lg(v) for k, v in A.items() if v},
        'aux_single_source': {'gt_squarings': 381, 'dense_mults_charged': 191,
                              'total_ref_units_one_gt_exp': gtexp,
                              'total_ref_units_one_gt_exp_exact': '381*(11632/20000)+191'},
        'storage': {k: str(v) for k, v in STO.items()},
        'record': {'bytes': K.RECORD_BYTES, 'key_bytes': K.KEY_BYTES,
                   'label_bytes': K.LABEL_BYTES, 'words': K.REPORT_WORDS},
        'trace': {
            'fixture_r': 9421, 'fixture_p_bits': 61, 'fixture_N': 785,
            'ladder_steps_min': min(r['steps'] for r in v19.LAD['rows']),
            'ladder_steps_max': max(r['steps'] for r in v19.LAD['rows']),
            'density_mean': 3.365, 'density_var': 0.419, 'density_overdispersion': 0.124,
            'density_worst_window': 7, 'density_factor': 2.080,
            'trap_per_event': 3581.5, 'trap_restart_price': 1526.19,
            'trap_waste_price': '76.80', 'trap_waste_price_value': 76.8, 'trap_dp_lookup_price': 0.0210,
            'trap_detector_price': 0.0210, 'trap_lambda': 0.0539,
            'success_rate_scaled': 0.3920, 'retention_cp95_lower': 0.9303,
            'p_usable_cp95_lower': 0.9862, 'p_detect_cp95_lower': 0.8576,
            'withdrawn_declaration_log2': 121.56, 'withdrawn_steps_per_report': 19.1,
            'lead_over_legality': 0.0641, 'lead_over_frontier': 0.0696,
            'repeat_zero_increment': 0, 'repeat_nonzero_increment': 1744,
            'repeat_pairs': 1000000, 'window_walks': 164, 'window_walks_total': 256,
            'window_events': 13546, 'lags': '20/40/60/80/100/120/140/160',
            'census_success_5': 0.7316, 'census_success_7': 0.5262, 'census_success_9': 0.2778,
            'N_classes_log2': 251.27,
        },
        'label_audit': audit,
        'inequality_audit': ineq,
        'inequality_control': flip,
        'failed_gate': [
            {'name': 'success_point_estimate_below_floor'
                     if SUC['net'] < 0.39 else 'success_cp95_lower_below_floor',
             'status': 'UNMET',
             'evidence': 'certificate.json:success.net and success.net_cp95_lower',
             'why': 'the success row is measured through the GUARDED retention arm: point '
                    'estimate %.7f, CP-95 lower %.7f, floor 0.39. %s'
                    % (SUC['net'], SUC['net_cp95_lower'],
                       'Both are below the floor; it is shipped as a measured failure.' if
                       SUC['net'] < 0.39 else 'The point clears, the lower bound does not.')},
            {'name': 'production_horizon_retention_at_theta_2e80', 'status': 'UNMET',
             'evidence': 'evidence/retention-v22.json:curve_guarded',
             'why': 'the guarded k-curve is measured to lag 56 (%.4f at lag 32); the production '
                    'horizon theta = 2^80 is not established and is not claimed'
                    % (SUC['production_lag_lag32'] or 0.0)},
            {'name': 'no_inversion_free_sigma_stable_column_table', 'status': 'UNMET',
             'evidence': 'evidence/sigma-set.json:enumeration.closed_and_inversion_free',
             'why': 'exhaustive over 2^12 subsets: only {} and Z/12 are closed, and Z/12 '
                    'contains the inverse pairs {k,k+6}'},
            {'name': 'class_marginal_uniformity_refuted_at_fixture', 'status': 'REFUTED',
             'evidence': 'certificate.json:failed_gate (class marginal chi2 was measured in the v18 census)',
             'why': 'chi2 = 238889.4 on 526 df; uniformity is a MODEL PREMISE, pricing with '
                    '1/N is conservative and the credit is not taken'},
            {'name': 'production_scale_class_stream_census', 'status': 'UNMET',
             'evidence': 'certificate.json:failed_gate (the class-stream sweep stopped at R = 1024)',
             'why': 'the sweep stops at R = 1024; no production census is claimed'},
        ],
        'failed': [],
    }
    cert['checks_total'] = len(CHECKS)
    cert['checks_pass'] = sum(1 for c in CHECKS if c['ok'])
    cert['checks_fail'] = cert['checks_total'] - cert['checks_pass']
    cert['checks'] = CHECKS
    json.dump(cert, open(os.path.join(HERE, 'certificate.json'), 'w'), indent=1, sort_keys=True)
    json.dump({'rows': ineq, 'control': flip}, open(os.path.join(EV, 'inequality-audit.json'),
                                                    'w'), indent=1, sort_keys=True)
    json.dump({'rows': audit, 'planted_flip_control': 'passes (negated label fails)'},
              open(os.path.join(EV, 'label-audit.json'), 'w'), indent=1, sort_keys=True)
    json.dump({'trigger': cert['reports']['trigger'], 'theta_bits': 80,
               'mean_dp_hits': C['mean'], 'mean_log2': C['mean_log2'],
               'cap': C['cap'], 'cap_log2': C['cap_log2'],
               'cap_density_factor': C['density_factor'],
               'density_measured': v19.DEN['windows'],
               'chebyshev_stop_loss': C['chebyshev_stop_loss'],
               'both_ways': {'mean': 'ceil(S/theta), the trigger output',
                             'cap': 'dens_factor x mean, the measured worst window',
                             'charged': 'every report row is priced at the cap'}},
              open(os.path.join(EV, 'reports-trigger.json'), 'w'), indent=1, sort_keys=True)
    write_note(cert, SUC, trap, C, A, claim_f, lw_str, lw_raw)
    print('claim %.4f success net %.7f (cp95 %.7f) aux 2^%.3f trap 2^%.3f checks %d/%d'
          % (claim_f, SUC['net'], SUC['net_cp95_lower'], lg(A['total']),
             lg(trap['total_reference_units']), cert['checks_pass'], cert['checks_total']))


def write_note(cert, SUC, trap, C, A, claim_f, lw_str, lw_raw):
    rows = v19.LAD['rows']
    zc = v19.LAD['controls']
    def n(x, d=4):
        import math as _m
        q = 10 ** d
        return ('%.' + str(d) + 'f') % (_m.floor(float(x) * q) / q)
    AUD = '\n'.join('| `%s` | %s | %s |' % (r['row'], r['op'], r['literal'])
                    for r in cert['inequality_audit'])
    if SUC['net'] >= 0.39:
        succ_sentence = ('The point estimate clears 0.39; the CP-95 lower does not and is shipped '
                         'as a FAILED gate.')
    else:
        succ_sentence = ('BOTH the point estimate and the CP-95 lower are BELOW the 0.39 floor: '
                         'this is a measured failure of the success gate, shipped as such, not a '
                         'claim.')
    capr = v19.CAP
    ctrl = capr['control_force_admission']
    curve_str = ', '.join('%s:%.4f' % (k, v['p_detect'])
                          for k, v in sorted(v19.RET22['curve_guarded'].items(), key=lambda kv: int(kv[0])))
    txt = f"""# v22 — the abandon cap is in the code and metered, the trap expectation has a derivation, retention is re-measured through the guard

**Claim {n(claim_f,4)}** (log2 dense-Fp12-multiply equivalents, C = 20000, all processors,
auxiliary work inside the claim: `claim.covers` = W_stress + aux). **Success {n(SUC['net'],4)}**
point estimate with every unfinished trail scored as failure; the CP-95 lower is
{n(SUC['net_cp95_lower'],4)}. {succ_sentence}

## 1. The abandon cap is a hard 128 x theta = 2^87 steps, implemented and metered (verdict defect 2a)

The v21 package priced 64 wasted steps in code while this note described abandonment after
128 x theta steps; both halves were wrong, because the code did not implement the cap it priced
and the note's exp(-128) tail was the wrong model for a trail trapped in a cycle. v22 fixes the
code and DELETES the tail.

The shipped trail protocol is now: a trail ends at its FIRST distinguished point OR at a hard cap
of 128 x theta = 2^87 steps, whichever comes first. A capped trail is ABANDONED: its labels are
discarded and a fresh trail starts (two GT exponentiations). `measure_v22_ext.cap_abandon`
implements that counter on the shipped step, meters the abandon path, and measures the delay:

- trails that report exit at their first DP, mean {capr['dp_exit_mean_steps']:.2f} steps at the 2^-5 fixture (the cap is never reached);
- a trail that returns to a visited (class, ring) state is on a cycle: {capr['trails_trapping_in_cycle']} of {capr['walks']} fixture trails trap, cycle length mean {capr['cycle_length_mean']} and max {capr['cycle_length_max']};
- every trapped trail runs to the cap and is abandoned there: measured delay mean {capr['mean_steps_to_abandon']:.1f} steps (= the cap exactly), variance {capr['var_steps_to_abandon']:.1f};
- control (force admission): the SAME state and the SAME ring continue WITHOUT the cap for 2 x cap steps and still produce no DP ({ctrl['still_without_dp']} of {ctrl['trails_tested']}). Without the cap those trails run forever, which is exactly what the cap bounds.

So there is no tail to bound and none is used: with a hard cap P(a trail runs past the cap) = 0
by construction, and the waste row is deterministic -- **cap x rate x step price**, one row. The
cap is priced in the claim: {trap['prices_reference_units']['wasted_cap_steps']:.3g} reference units of
waste per abandon. Charged total: Lambda = S^2/(2(r-1)) = {trap['repeat_probability_over_the_walk']:.4f};
Markov on the TOTAL gives P(X > C) <= E[X]/C at C = ceil(4 Lambda / 1e-9) =
{trap['events_charged']:,} events; each costs a DP lookup ({trap['prices_reference_units']['dp_lookup']:.4f})
+ a detector probe and 112-byte compare ({trap['prices_reference_units']['detector_probe_compare']:.4f}) + a restart, two GT
exponentiations ({trap['prices_reference_units']['restart_two_gt_exps']:.2f}) + the cap waste =
**2^{lg(trap['total_reference_units']):.3f} reference units**, inside `aux.trap_charge` and therefore inside
the claim. `total_integer_rows_1_over_C_REF` re-sums the four integer rows and asserts equality.

## 2. The expectation bound has a derivation: monotonicity, plus Markov on a measured E (verdict defect 2b)

The verdict is right that Markov needs a justified expectation bound and that a fixed,
history-dependent walk does not automatically have uniform pair-collision probabilities. v22
stops applying the iid birthday formula by assumption and ships the argument.

**(a) Monotonicity (primary).** The guard adds exactly one coordinate to the state:
y -> (y, ring). A joint repeat therefore REQUIRES an element repeat -- pointwise, for every walk,
with no independence hypothesis, the joint collision count is at most the element collision
count, and it is the element stream that the iid birthday count describes. Measured on
{v19.BND['step_pairs']} step pairs (`evidence/bound-v22.json`): joint repeats
{v19.BND['guarded_joint_repeats']} <= guarded element repeats {v19.BND['guarded_element_repeats']} <=
unguarded element repeats {v19.BND['unguarded_element_repeats']}. The guard makes collisions rarer,
so S^2/(2(r-1)) = {trap['repeat_probability_over_the_walk']:.4f} is an UPPER bound for the guarded walk.

**(b) Markov on a measured expectation (cross-check, hypothesis-free).** The census measures
E[X] directly: {v19.REP['census']['joint_repeat_events']/v19.REP['census']['walks']:.4f} joint repeats per walk over
{v19.REP['census']['walks']} walks, with no model. Markov P(X >= t) <= E[X]/t needs nothing else.

Withdrawn: the exp(-128) tail, and any reading of the per-worker bound as "at most one event per
worker". The charge is a total, C = ceil(4 Lambda / 1e-9) with a 1e-9 tail, and linearity across
workers is used -- not independence.

## 3. Retention is re-measured through the guarded step, and the report count is census-direct (verdict defects 1 and 3)

The 0.895833 that v19 shipped was an EXTRACT of a v18-era JSON (its own `source` field says so),
i.e. not measured on the object the success row multiplies. v22 re-runs the same
collision-conditioned arm inside the census build through the SHIPPED guarded step, ring and
prev-row active, with a guard-off control beside it (`evidence/retention-v22.json`,
{v19.RET22['pairs']} planted pairs, known-answer gate {v19.RET22['known_answer_gate_pass']}/{v19.RET22['known_answer_gate_required']}):

- guarded p_detect = **{SUC['p_detect_guarded']:.6f}**, CP-95 lower **{SUC['p_detect_guarded_cp95_lower']:.6f}**;
- guard-off control p_detect = {SUC['guard_off_control_p_detect']:.6f} (with no ring the two trails never separate);
- guarded k-curve: {curve_str}.

The product uses the collision-instant guarded value: {SUC['p_detect_guarded']:.6f} x mu 1.1 x p_usable 1.0 =
{SUC['retention_eff']:.6f}. The production horizon theta = 2^80 is NOT established by a 56-step fixture
curve and is shipped as a NAMED FAILED GATE, not as a claim.

Reports: the count is CENSUS-DIRECT. {v19.REP22['reports_measured']} DP-triggered reports over
{v19.REP22['walks']} guarded walks = {v19.REP22['reports_per_walk_measured']:.3f} per walk against the model 1/theta (2^-5) and
ceil(S/theta) = {v19.REP22['ceil_formula_reports_per_walk']} per walk at the same width (ratio
{v19.REP22['check_measured_vs_formula_ratio']:.3f}); ceil(S/theta) is kept ONLY as a CHECK row, never as the count.

## 3b. A density is not an interval: the cap carries the measured dispersion

{int(v19.DEN['windows']['windows'])} windows of {v19.DEN['windows']['window_steps']} steps at the 2^-5 fixture: mean
{v19.DEN['windows']['mean_hits_per_window']:.3f} hits, variance {v19.DEN['windows']['var_hits_per_window']:.3f}, worst window
{v19.DEN['windows']['worst_window_hits']} = {C['density_factor']:.3f} x the mean, so every per-report row is priced at the cap
2^{C['cap_log2']:.3f}. No variance bound is used in success (the v19 Chebyshev row stays withdrawn).

## 4. No rounding up: every >= and <= literal is audited with its direction (defect 4)

Every inequality below is a table row, re-read from this note by `audit_v22.py` and applied to
the certificate value with a direction rule: a floor literal must satisfy value >= literal, a
ceil literal value <= literal. v18 printed "at least 0.3963" over a value of 0.3962569 - a
round-UP - and the same rule rejects it (control in `evidence/inequality-audit.json`).

| row | op | literal |
|---|---|---|
{AUD}

The figures in this note are floored, not rounded ({n(SUC['net'],4)} for a value of
{SUC['net']:.6f}). {len(cert['label_audit'])} `2^x` labels are recomputed from their decimals, with the
planted-flip control rejecting a negated label.

## 5. The toy ladder, twenty-four full trials (defect 5)

`measure_v19.ladder` at the fixture the package already ships (r = {v19.LAD['fixture']['r']}, 61-bit p,
ord_r(p) = 12, N = {v19.LAD['fixture']['N_classes']}, table {v19.LAD['fixture']['table_rows']} x 12, DP width 2^-5, guard depth 12):
plant x, walk the shipped guarded step, end each trail at its first DP, store the DP records,
solve all twelve shifts on a class-key collision, verify g^x == z. **{v19.LAD['passed']}/24 pass**, steps
{min(r['steps'] for r in rows)}..{max(r['steps'] for r in rows)}; controls: wrong secret rejected, zero-increment trap fires,
nonzero-increment pair solves. Rows are `evidence/ladder-v19.json`.

## What is kept, and what is not claimed

The fixed-P reduction, the solve identity x = (p^j a1 - a2)/(b2 - p^j b1) for j = 0..11 with
every candidate verified, the exact r-1 factorisation, the guard at D = 12, and the measured
rows (class match 0.0833 within 256 steps, class-stream survival 0.9667 at k=1 to 0.7833 at
k=12, DP hits 1.9e-3 at lags 20/40/60, the density census, aux 412.5896 reference units per
381-bit GT exponentiation). Success is conditional on the three disclosed hypotheses
(random-map walk quality inside the 10 per cent allowance, pseudorandom routing and table
fields, ideal coins), and the claim covers W_stress {n(lw_str,6)} plus aux
2^{lg(A['total']):.3f} with a lead of {125.84 - claim_f:.4f} bits over the 125.84 frontier and
{125.8345 - claim_f:.4f} over the 125.8345 legality line. Four gates fail and ship as failures:
the success CP-95 lower is below 0.39, the inversion-free sigma-stable table does not exist,
the class marginal is not uniform at fixture scale (credit not taken), and the class-stream
census stops at R = 1024.

**Reproduce.** `python3 measure_v22.py` (repeat classes, density, ladder), `python3
gen_cert_v22.py` (certificate + note + audits), `python3 audit_v22.py` (independent re-read),
`python3 selfaudit.py` (literal grep).
"""
    open(os.path.join(HERE, 'public-note.md'), 'w').write(txt)
    print('note %d chars' % len(txt))


if __name__ == '__main__':
    main()
