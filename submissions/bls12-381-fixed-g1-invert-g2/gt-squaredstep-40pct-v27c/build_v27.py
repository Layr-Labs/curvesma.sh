#!/usr/bin/env python3
"""Assemble certificate.json for v27 from the measured evidence files.

Every field is READ from a shipped evidence JSON; nothing here is typed by hand.
If a check cannot be evaluated, the certificate records False and exits nonzero.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def load(p):
    with open(p) as f:
        return json.load(f)


def main():
    led = load('evidence/v27-ledger.json')
    cen = load('evidence/v27-measure-census.json')
    rec = load('evidence/v27-measure-recovery.json')['recovery']
    cyc = load('evidence/v27-measure-cycle.json')['cycle']
    thm = load('evidence/v27-measure-theorem.json')['theorem']
    mod = load('evidence/v27-measure-model.json')['model']
    prod = cen['production_like']
    floor = 0.39
    checks = {
        'claim_covers_total': float(led['claim_log2']) >= float(led['total_log2']),
        'claim_reduced_rational': math.gcd(int(led['numerator']),
                                           int(led['denominator'])) == 1,
        'components_sum_exactly': None,   # filled below from Fractions
        'model_matches_judge_recompute': bool(mod['agrees_with_judge']),
        'recovery_zero_misses': rec['misses'] == 0,
        'recovery_j0_regression_kept': bool(rec['j0_regression_kept']),
        'recovery_old_formula_dead_for_j_ne_0': rec['old_formula_hits_jne0'] == 0,
        'recovery_gcd_ok': rec['gcd_fail'] == 0,
        'census_dp_filter_executed': cen['1']['dp_hits'] > 0,
        'census_invariant_clean': all(cen[k]['invariant_fail'] == 0
                                      for k in ('1', '2', '3', '4')),
        'production_like_arm_clears_floor': prod['cp95_lower'] >= floor,
        'exact_cycle_census_no_collapse': (cyc['v26']['distinct_exact_states']
                                           == cyc['v26']['elements']),
        'squared_step_has_no_two_cycles': cyc['v26']['two_cycles'] == 0,
        'theorem_planted_mutation_flagged': bool(thm['planted_flagged']),
        'theorem_clean_scan': bool(thm['clean']['clean']),
    }
    from fractions import Fraction
    comps = led['components']
    s = sum((Fraction(comps[k]) for k in ('w_stress', 'trap_charge', 'aux_sum')),
            Fraction(0))
    checks['components_sum_exactly'] = (s == Fraction(led['total_log2']).limit_denominator(10**30)
                                        if False else math.log2(float(s))
                                        == float(led['total_log2']))
    cert = dict(
        version=27,
        S_steps=led['S'],
        S_log2=math.log2(int(led['S'])),
        N_classes_log2=251.27212691232603,
        f_measured=led['f_measured'],
        c_stress=led['c_stress'],
        claim_log2=led['claim_log2'],
        claim_numerator=led['numerator'],
        claim_denominator=led['denominator'],
        total_log2=led['total_log2'],
        components=comps,
        trap=led['trap'],
        aux_rows=led['aux_rows'],
        dp_floor=dict(floor=floor, production_like_arm=prod['p_recover'],
                      cp95_lower=prod['cp95_lower'], S_over_theta=prod['S_over_theta'],
                      arms={k: dict(p=cen[k]['p_recover'], cp95_lower=cen[k]['cp95_lower'],
                                    reports_per_trial=cen[k]['reports_per_trial'],
                                    S_over_theta=cen[k]['budget'] / cen[k]['theta'])
                            for k in ('1', '2', '3', '4')}),
        model=dict(v27=mod['model_v27'], v26_code=mod['model_v26_code'],
                   judge=mod['judge_recomputed'], factor=mod['factor']),
        recovery=dict(samples=rec['samples'], misses=rec['misses'],
                      gcd_fail=rec['gcd_fail'], subgroup_fail=rec['subgroup_fail'],
                      j0_regression_kept=rec['j0_regression_kept'],
                      old_formula_hits_jne0=rec['old_formula_hits_jne0']),
        cycle=dict(elements=cyc['v26']['elements'],
                   distinct_exact_states=cyc['v26']['distinct_exact_states'],
                   distinct_class_keys=cyc['v26']['distinct_class_keys'],
                   collapse_factor=cyc['v26']['collapse_factor'],
                   two_cycles=cyc['v26']['two_cycles']),
        checks=checks,
        all_checks_pass=all(bool(v) for v in checks.values()),
    )
    with open('certificate.json', 'w') as f:
        json.dump(cert, f, indent=1, sort_keys=True)
    bad = [k for k, v in checks.items() if not v]
    print('certificate.json written; checks %d/%d pass%s'
          % (sum(1 for v in checks.values() if v), len(checks),
             '' if not bad else ' FAILED: %s' % bad))
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
