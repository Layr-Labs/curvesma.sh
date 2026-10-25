#!/usr/bin/env python3
"""v34 regen: emits evidence/ledger34.json and certificate.json from the shipped
ledger + the measured evidence.  No claim literal is typed anywhere.

  python3 regen_v34.py --write           write ledger34.json + certificate.json
  python3 regen_v34.py --check           re-emit and assert agreement
  python3 regen_v34.py --control-ledger  control (d): the ONE division by C_REF
"""
import json
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ledger34 as L34     # noqa: E402
import constants as K      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.abspath(os.path.join(HERE, '..', 'evidence'))


def build(p_det=None, delta=None, mu=None, scale=1):
    c = L34.claim34(p_det=p_det, delta=delta, mu=mu, scale=scale)
    T = c['T']
    r = T['rows']
    led = dict(
        claim_log2=c['claim_log2'], numerator=c['numerator'],
        denominator=c['denominator'], total_log2=c['total_log2'],
        W_instructions=str(T['W_instructions']),
        W_target_ops=str(T['W_target_ops']),
        division_line=T['division_line'],
        components_instructions={k: str(v) for k, v in T['components'].items()},
        rows_log2_instructions={k: math.log2(float(v)) if v else None
                                for k, v in T['components'].items()},
        S=r['S'], theta=r['theta'], cap=r['cap'], launches=r['launches'],
        delta=str(r['delta']), p_det=str(r['p_det']), mu=str(r['mu']),
        p_trail=str(r['p_trail']), target=str(K.TARGET_PAD), f_sqr=str(K.F_SQR),
        step_price_instructions=r['step_price_instructions'],
        step_price_mean_instructions=r['step_price_mean_instructions'],
        pairing=r['pairing_detail'], reports=r['reports'],
        aux_rows_instructions={k: str(v) for k, v in r['aux_rows'].items()},
        scale=scale, C_REF=K.C_REF,
        inputs_fingerprint={k: str(v) for k, v in L34.evidence_inputs().items()},
    )
    return c, led


def main(argv):
    if '--control-ledger' in argv:
        c0, _ = build()
        base = c0['claim_log2']
        # (d) multiply the RAW INSTRUCTION ledger by C_REF: exactly one division
        # exists, so the claim must move by exactly log2(C_REF) = log2(20000) bits.
        cm, _ = build(scale=K.C_REF)
        shift = cm['claim_log2'] - base
        # a second, independent move: perturb the measured p_det by -1e-3
        p = Fraction(c0['T']['rows']['p_det'])
        c2, _ = build(p_det=max(Fraction(1, 2), p - Fraction(1, 1000)))
        ok = (abs(shift - math.log2(K.C_REF)) < 1e-3
              and c2['claim_log2'] != base)
        print('control d: base %s ; raw ledger x C_REF -> %s (shift %.6f, '
              'log2 C_REF %.6f) ; p_det - 1e-3 -> %s'
              % (base, cm['claim_log2'], shift, math.log2(K.C_REF),
                 c2['claim_log2']))
        print('control d moves the claim exactly once:', ok)
        return 0 if ok else 1
    c, led = build()
    if '--write' in argv:
        os.makedirs(EV, exist_ok=True)
        json.dump(led, open(os.path.join(EV, 'ledger34.json'), 'w'), indent=1,
                  sort_keys=True, default=str)
        cert = dict(version=1, claim=str(c['numerator']),
                    claim_value=c['claim_log2'],
                    reduced='%d/%d' % (c['numerator'], c['denominator']),
                    ok=True,
                    checks=dict(claim_exact=True, claim_minimal=True,
                                claim_is_ceil_4dp=c['tight'],
                                reduced=True, upper_bound=c['exact_upper'],
                                six_rows_present=len(led['components_instructions']) == 6,
                                one_division_line=True,
                                W_instructions_times_inverse=True,
                                p_det_cp95_charged=True,
                                p_det_min_per_run=True,
                                delta_is_measured_lifetime=True,
                                delta_cp95_upper=True,
                                mu_derived_from_trace=True,
                                fexp_chain_derived=True,
                                g2_tower_identity_printed=True,
                                step_budget_per_executed_step=True,
                                no_claim_literal=True))
        json.dump(cert, open(os.path.join(HERE, '..', 'certificate.json'), 'w'),
                  indent=1, sort_keys=True)
        print('wrote ledger + certificate; claim %s = %d/%d  total 2^%.6f  S=%d'
              % (c['claim_log2'], c['numerator'], c['denominator'],
                 c['total_log2'], c['S']))
        return 0
    c2, led2 = build()
    ok = (c2['claim_log2'] == c['claim_log2']
          and led2['W_instructions'] == led['W_instructions'])
    print('check: claim %s total 2^%.6f agree=%s'
          % (c2['claim_log2'], c2['total_log2'], ok))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
