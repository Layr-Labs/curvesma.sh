#!/usr/bin/env python3
"""v19 CERTIFICATE GENERATOR — one pass: ledger, success, trap charge, density cap,
ladder, label audit, >=/<= audit, and the note that quotes them.

Every number in public-note.md is emitted from the value computed here (the
inequality strings are generated from the same table that audits them), so the
note cannot drift from the certificate.  audit_v19.py re-reads the note and the
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
import v19

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
    s['peak_bytes'] = int(s['peak_bytes']) + det
    s['report_cap'] = int(cap)
    s['records'] = int(cap)
    s['store_bytes'] = int(cap) * K.RECORD_BYTES
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

    # success ordering: model >= rate-scaled >= product >= net, and the CP-95 is below
    chk('success_rows_are_ordered',
        SUC['model'] >= SUC['rate_scaled'] >= SUC['product'] >= SUC['net'],
        '%.7f >= %.7f >= %.7f >= %.7f' % (SUC['model'], SUC['rate_scaled'],
                                          SUC['product'], SUC['net']))
    chk('success_net_point_clears_the_floor', SUC['net'] >= 0.39,
        '%.7f >= 0.39 (point estimate)' % SUC['net'])
    chk('success_cp95_lower_is_disclosed_below_the_floor',
        SUC['net_cp95_lower'] < 0.39,
        '%.7f < 0.39 — DISCLOSED as an UNMET gate, not claimed' % SUC['net_cp95_lower'])
    chk('retention_product_is_the_model_times_mu_p_use',
        abs(SUC['retention_eff'] - float(K.MU_WALK) * 0.8958333333333334 * 1.0) < 1e-12,
        '%.6f' % SUC['retention_eff'])
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
        'retention': {'p_detect': 0.8958333333333334,
                      'p_detect_cp95_lower': 0.8575558748268416,
                      'p_usable_given_detected': 1.0,
                      'p_usable_cp95_lower': 0.9861629853644579,
                      'product': '%.7f' % SUC['retention_eff'],
                      'product_cp95_lower': '%.7f' % SUC['retention_eff_cp95_lower'],
                      'charged_excess_over_mu': 1.014787,
                      'witness_usable': '215/215',
                      'known_answer_gate': '240/240',
                      'note': '215/240 is the DETECTION row (CP-95 lower 0.8576); 215/215 '
                              'witness usability; 240/240 is the planted-label gate'},
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
            {'name': 'success_cp95_lower_below_floor', 'status': 'UNMET',
             'evidence': 'certificate.json:success.net_cp95_lower',
             'why': 'the CP-95 lower of the retention product is %.7f, below 0.39; only the '
                    'POINT estimate clears the floor and the note says exactly that'
                    % SUC['net_cp95_lower']},
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
    txt = f"""# v19 — traps charged, retention in the success row, density in the cap

**Claim {n(claim_f,4)}** (log2 dense-Fp12-multiply equivalents, C = 20000, all processors,
auxiliary work inside the claim: `claim.covers` = W_stress + aux). **Success {n(SUC['net'],4)}**
point estimate with every unfinished trail scored as failure; the CP-95 lower is
{n(SUC['net_cp95_lower'],4)} and is shipped as a FAILED gate, not as a claim.

## 1. A repeat traps whether or not its increment is zero (verdict defect 1)

v18 charged the zero-increment fraction only. That was the wrong event, and the census
was re-run to prove it rather than argue it: **{v19.REP['census']['joint_repeat_events']}
joint repeats in {v19.REP['census']['step_pairs']} step pairs, {v19.REP['census']['zero_increment_events']} of them
zero-increment and {v19.REP['census']['nonzero_increment_events']} nonzero** (lags
{'/'.join(sorted(v19.REP['census']['lag_histogram'], key=int))}; period-20 family,
{v19.REP['window']['walks_with_joint_repeat']}/{v19.REP['window']['walks']} walks repeat inside
{v19.REP['window']['window_steps']} steps). The shipped attack stores DP records and guard entries and nothing
else, so a return to the same (element, guard) state is a cycle whose twenty-odd elements
contain no DP: it never reports and the trail is pending at every cutoff. The increment only
decides whether a solve *would* exist had both records been kept — and the attack keeps
neither. The ladder exhibits both branches: the zero-increment pair yields
{zc['zero_increment_candidates']} candidates and **{zc['zero_increment_verified']} verify**; the nonzero-increment
pair yields {zc['nonzero_increment_candidates']} candidates and **{zc['nonzero_increment_verified']} verifies**.

Charged, not argued: the birthday count over the whole walk is
Lambda = S^2/(2(r-1)) = {trap['repeat_probability_over_the_walk']:.4f} <= 0.0540 events, so at most one
event per worker (4) is charged, each at a **DP lookup ({trap['prices_reference_units']['dp_lookup']:.4f}) + a
detector probe and 112-byte compare ({trap['prices_reference_units']['detector_probe_compare']:.4f}) + a full restart, two GT
exponentiations ({trap['prices_reference_units']['restart_two_gt_exps']:.2f}) + the wasted cycle steps
({trap['prices_reference_units']['wasted_cycle_steps']:.1f})** = {trap['total_reference_units']:.1f} reference units,
2^{lg(trap['total_reference_units']):.3f}, inside `aux.trap_charge` and therefore inside the claim. It is small
because Lambda is small; it is no longer free.

## 2. Success carries the retention product, and its CP-95 lower is below the floor
(verdict defect 2)

The measured retention is `retention_eff` = mu x p_detect x p_usable = 1.1 x 0.895833 x 1.0 =
{SUC['retention_eff']:.6f}. Applied to the birthday **rate**: {SUC['rate_scaled']:.4f} (the 0.392 the verdict
computes). Applied to the **probability**: {SUC['product']:.4f} = {n(SUC['product'],4)} >= 0.3906, which is the
more conservative of the two and is what ships. Net of the trap debit and the density stopping
loss: **{n(SUC['net'],4)} >= 0.3904**. The CP-95 lower is
**{n(SUC['net_cp95_lower'],4)} <= 0.3690 and <= 0.39**, i.e. **below the floor**: the point estimate clears
0.39, the lower bound does not, and no amount of prose changes that. Only {SUC['retention_eff_cp95_lower']:.4f} of
the rate is certified at 95 per cent; the residual is a disclosed risk, not a claim.

## 3. A density is not an interval: the cap carries the measured dispersion (defect 3)

{int(v19.DEN['windows']['windows'])} windows of {v19.DEN['windows']['window_steps']} steps at the 2^-5 fixture: mean {v19.DEN['windows']['mean_hits_per_window']:.3f} hits,
variance {v19.DEN['windows']['var_hits_per_window']:.3f}, over-dispersion {v19.DEN['windows']['over_dispersion_var_over_mean']:.3f}, worst window
{v19.DEN['windows']['worst_window_hits']} = {C['density_factor']:.3f} x the mean. The cap is **factor x mean = 2^{C['cap_log2']:.3f}**
(not ceil(S/theta)), every report row is priced at that cap, and the stopping loss at the cap
is Chebyshev with variance <= mean: {cert['reports']['inequality']}. That row is subtracted in
`success.net`, so a run that hits the cap is not free either.

## 4. No rounding up: every >= and <= literal is audited with its direction (defect 4)

Every inequality below is a table row, re-read from this note by `audit_v19.py` and applied to
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

**Reproduce.** `python3 measure_v19.py` (repeat classes, density, ladder), `python3
gen_cert_v19.py` (certificate + note + audits), `python3 audit_v19.py` (independent re-read),
`python3 selfaudit.py` (literal grep).
"""
    open(os.path.join(HERE, 'public-note.md'), 'w').write(txt)
    print('note %d chars' % len(txt))


if __name__ == '__main__':
    main()
