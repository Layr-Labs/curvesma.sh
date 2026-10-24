#!/usr/bin/env python3
"""v16: the three arithmetic repairs, each as an exact row set.

(1) PATTERN MULTIPLICITIES.  v15 priced longer cycles with
        C(L, L/2) * R^{-L/2} * theta
    which counts pattern CLASSES once each.  A class with k labelings enters the
    union bound k times, so the count must be multiplied back.  The table below
    enumerates every word of length L over a small alphabet, groups the words by
    their row-multiset class, and shows that
        sum over classes of  mult(class) * R^{-L}
    is exactly R^L / R^L = 1: the multiplicities restore the WHOLE word space.
    There is no residual "few patterns" saving, so the pattern route is withdrawn
    entirely and replaced by (2).  (This is the reviewer's "pattern
    multiplicities, sum p^2 inference" defect, closed by not using the inference.)

(2) THE TRAP, counted on JOINT STATES instead of row patterns.  The step is a
    deterministic function of the augmented state (element, ring tuple).  An
    absorbing cycle is exactly a repeat of an augmented state, so the cycle
    event is a collision in a state space of size
        M = (r-1) * (12 R)^D
    and the union bound over the S step pairs of one walk is
        tau <= S^2 / (2 M)          (random-map / ideal-table hypothesis,
                                     the same hypothesis class as the rest of
                                     the package, disclosed not hidden).
    No pattern is enumerated, no theta is extrapolated, no 1/(12R) factor is
    multiplied: the multiplicity defect cannot arise.

(3) REPORTS AND THE CAP.  Reports are DP hits, one per distinguished class
    visited.  Mean hits = S * 2^-80; the cap is 2x the mean (a stated slack),
    and a Chernoff row bounds the probability that a run exceeds it.  The
    measured census calibration (records per expected hit at the shipped DP
    width) is shipped beside the bound, not in place of it.

(4) AUX, ONE NUMBER.  381 * (11632/20000) + 191 = 412.5896 reference units per
    381-bit GT exponentiation, written once into evidence/aux-v16.json and read
    by every document.
"""
import json
import math
import os
import sys
from collections import Counter
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')

R_PROD = 1 << K.TABLE_BITS          # 2^14 rows
N_ELEM = K.R_ORDER - 1              # nonidentity elements of the order-r subgroup
COLS = K.TABLE_COLS                 # 12
LOG2_S = 125.828667                 # log2 of the step budget S (ledger.py)
S = None                            # filled from constants


def s_value():
    import ledger
    return ledger.ledgers()['S']


def pattern_multiplicity_table(Ls=(4, 6), R=5):
    """Every word of length L over R letters, grouped by multiset class."""
    rows = {}
    for L in Ls:
        words = [[]]
        for _ in range(L):
            words = [w + [t] for w in words for t in range(R)]
        classes = Counter(tuple(sorted(w)) for w in words)
        total = sum(classes.values())
        # the v15 form: classes counted once, C(L, L//2) * R^-(L/2)
        v15_style = math.comb(L, L // 2) * R ** (-(L / 2))
        rows[str(L)] = {
            'words': total,
            'word_space_R_pow_L': R ** L,
            'classes': len(classes),
            'sum_of_class_multiplicities': total,
            'multiplicities_restore_the_word_space': total == R ** L,
            'max_class_multiplicity': max(classes.values()),
            'mean_class_multiplicity': total / len(classes),
            'v15_style_count_C_L_L2_times_R_pow_minus_L2': v15_style,
            'undercount_factor_vs_mean_multiplicity': (total / len(classes)) / v15_style,
            'sample_classes': [{'multiset': list(k), 'multiplicity': v}
                               for k, v in sorted(classes.items(),
                                                  key=lambda kv: -kv[1])[:5]],
        }
    return rows


def trap_joint(D, S):
    ring = (COLS * R_PROD) ** D
    M = N_ELEM * ring
    tau = Fraction(S * S, 2 * M)
    tau_log2 = 2 * math.log2(S) - 1 - math.log2(M)
    return {
        'guard_depth': D,
        'ring_alphabet': COLS * R_PROD,
        'ring_states_log2': D * math.log2(COLS * R_PROD),
        'joint_states_M_log2': math.log2(N_ELEM) + D * math.log2(COLS * R_PROD),
        'tau_union_bound_log2': tau_log2,
        'tau_numerator_log2': 2 * math.log2(S),
        'tau_expression': 'S^2 / (2 * (r-1) * (12R)^D)',
        'guard_instructions_per_mul_step': 3 * D + 12,
        'v15_spurious_bound_log2': -50.7,
    }


def reports(S, theta_bits=K.THETA_BITS):
    th = 1 << theta_bits
    mean = -(-S // th)                       # ceil(S/theta)
    cap = 2 * mean                           # stated slack, cap >= 2 x observed max
    cher_log2 = -mean / (3.0 * math.log(2.0))   # P(X > 2 mean) <= exp(-mean/3)
    return {
        'theta_bits': theta_bits,
        'mean_dp_hits': mean,
        'mean_dp_hits_log2': math.log2(mean),
        'cap': cap,
        'cap_log2': math.log2(cap),
        'cap_over_mean': 2,
        'chernoff_P_exceed_cap_log2': cher_log2,
        'chernoff_form': 'P(X > 2 mean) <= exp(-mean/3) for independent DP hits',
        'overflow_policy': 'drop-oldest + drain, charged as aux.report_cap_drain at the CAP',
        'measured_calibration': {
            'source': 'evidence/t4-dp5.json',
            'dp_bits': 5,
            'trials': 6400,
            'records': 20240,
            'records_per_trial': 3.1625,
            'expected_hits': None,
            'note': 'records <= cap = 512 in 6400/6400 trials; the ratio row is a '
                    'calibration, the BOUND is the cap with the Chernoff row.',
        },
    }


def aux_single_source():
    cyclo = K.cyclo_sqr_instructions()
    sq = Fraction(381 * cyclo, K.C_REF)
    mul = 191
    total = sq + mul
    return {
        'cyclo_sqr_instructions': cyclo,
        'gt_squarings': 381,
        'dense_mults_charged': mul,
        'gt_squaring_ref_units': float(sq / 381),
        'sq_block_ref_units': float(sq),
        'total_ref_units_one_gt_exp': float(total),
        'total_ref_units_one_gt_exp_exact': '381*(%d/20000)+191' % cyclo,
        'total_ref_units_one_gt_exp_6dp': round(float(total), 6),
        'superseded_value': 411.59,
        'superseded_explanation': '411.59 = 381*(11632/20000) + 190, i.e. the MEAN '
                                  'Hamming weight 190; v16 charges 191 in one field '
                                  'and every document reads this file.',
    }


def main():
    os.makedirs(EV, exist_ok=True)
    S = s_value()
    out = {
        'pattern_multiplicity': pattern_multiplicity_table(),
        'trap': {
            'method': 'joint-state counting; no pattern enumeration',
            'S_log2': LOG2_S,
            'S': str(S),
            'depths': [trap_joint(D, S) for D in (12, 24, 48, 128)],
            'selected_depth': 12,
            'selection_reason': 'D=12 has the lowest guard cost (48 instructions '
                                'per multiply step), the highest measured class-'
                                'stream survival, and a trap bound of 2^-214 that '
                                'is already far beyond the 2^-100 class; larger D '
                                'buys trap margin that is not needed and degrades '
                                'measured class-stream survival (class-retention.json).',
        },
        'reports': reports(S),
        'aux': aux_single_source(),
    }
    with open(os.path.join(EV, 'trap-v16.json'), 'w') as f:
        json.dump(out, f, indent=1, sort_keys=True)
    t = out['trap']['depths'][0]
    print("trap D=12  M=2^%.1f  tau=2^%.2f" % (t['joint_states_M_log2'],
                                               t['tau_union_bound_log2']))
    for d in out['trap']['depths']:
        print("  D=%3d  M=2^%.1f  tau=2^%.1f" % (d['guard_depth'],
                                                 d['joint_states_M_log2'],
                                                 d['tau_union_bound_log2']))
    for L, r in out['pattern_multiplicity'].items():
        print("  L=%s words=%d classes=%d mean_mult=%.1f undercount=%.1fx"
              % (L, r['words'], r['classes'], r['mean_class_multiplicity'],
                 r['undercount_factor_vs_mean_multiplicity']))
    print("reports mean=%.3e cap=%.3e chernoff=2^%.1f"
          % (out['reports']['mean_dp_hits'], out['reports']['cap'],
             out['reports']['chernoff_P_exceed_cap_log2']))
    print("aux one GT exp = %.6f ref units" % out['aux']['total_ref_units_one_gt_exp'])


if __name__ == '__main__':
    main()
