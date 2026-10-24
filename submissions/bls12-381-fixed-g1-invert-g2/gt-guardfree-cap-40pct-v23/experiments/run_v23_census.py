#!/usr/bin/env python3
"""v23 follow-up arms: (a) reports WITH restart in a DP-rate regime, (b) the lag-32
agreement arm on the guard-free walk, (c) the birthday cross-check of the guard-free
class stream at an UNSATURATED fixture length."""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import measure_v23 as M

EV = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'evidence')


def main():
    t0 = time.time()
    rep = M.reports_restart(walks=120, steps=2048, dp_bits=3, bits=7, cap_mult=128)
    rep['regime'] = ('theta = 2^3 = 8 steps against a 785-class pool (fixture trap rate ~4% '
                     'per segment), so the step budget reaches the DP rate instead of being '
                     'eaten by abandons -- the regime in which the report COUNT is the DP '
                     'rate.  v22 never restarted, so its continuing walk re-hit the same DP '
                     'classes and counted 53.05 against 64.')
    json.dump(rep, open(os.path.join(EV, 'reports-v23.json'), 'w'), indent=1, sort_keys=True)
    print('reports %.3f/walk vs %d ratio %.4f abandons %d [%.0fs]'
          % (rep['reports_per_walk_measured'], rep['ceil_formula_reports_per_walk'],
             rep['check_measured_vs_formula_ratio'], rep['abandons_measured'],
             time.time() - t0), flush=True)

    lag = M.retention_dpgated(pairs=2000, steps=64, dp_bits=7, bits=5, seed=2511,
                              lags=(1, 2, 4, 8, 16, 32, 56))
    lag['arm'] = ('LAG ARM: theta = 2^7 = 128 steps so the segment does not end before lag '
                  '56; pairs that reach lag k are counted at lag k (the conditioning the v22 '
                  'verdict used: 0.10625 at lag 32 with the guard).')
    json.dump(lag, open(os.path.join(EV, 'retention-lag-v23.json'), 'w'), indent=1,
              sort_keys=True)
    print('lag 32 %.4f n=%d | 56 %.4f n=%d | dpgated %.5f [%.0fs]'
          % (lag['curve_any_step']['32']['p_agree_any_step'],
             lag['curve_any_step']['32']['pairs'],
             lag['curve_any_step']['56']['p_agree_any_step'],
             lag['curve_any_step']['56']['pairs'],
             lag['dp_gated_agree']['p_detect_given_dp_hit'], time.time() - t0), flush=True)

    b = M.bound_guardfree(walks=1600, steps=20, bits=5, seed=4242)
    f = 20 * 19 / 2.0 / 785.0
    b['birthday_formula_per_walk'] = f
    b['measured_over_formula'] = b['E_repeats_per_walk_measured'] / f
    b['regime'] = ('L = 20 steps against a 785-class pool: UNSATURATED (L^2/(2N) = %.4f), so '
                   'the measured repeat count is a direct check of the birthday value used by '
                   'the trap row and the Markov cross-check; the L = 200 census is the '
                   'saturated arithmetic demo, not a rate.' % f)
    json.dump(b, open(os.path.join(EV, 'bound-v23.json'), 'w'), indent=1, sort_keys=True)
    print('bound %.4f/walk vs %.4f ratio %.3f [%.0fs]'
          % (b['E_repeats_per_walk_measured'], f, b['measured_over_formula'],
             time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
