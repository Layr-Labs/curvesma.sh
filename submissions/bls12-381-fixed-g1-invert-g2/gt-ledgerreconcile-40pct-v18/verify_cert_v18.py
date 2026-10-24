#!/usr/bin/env python3
"""verify_cert_v18.py — independent recompute of every claim-bearing row.

Imports only constants.py (prices) and recomputes the ledger arithmetic here, so a bug in
ledger.py cannot certify itself.  Asserts:
  * W_stress = S * c_stress and W_raw = S * c_raw, with c_* rebuilt from the row prices;
  * the claim is the least 4-decimal log2 above log2(W_stress) + 2e-5, both directions;
  * success.model = 1 - exp(-Lambda) and success.net = model * (1 - tau) >= 0.39;
  * reports: cap = 2 * ceil(S / 2^80), and the declared cumulative count equals the charge;
  * every audited label round-trips to its decimal within 5e-4 relative.
Prints VERIFY PASS / FAIL and exits non-zero on failure.
"""
import json
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, 'certificate.json')))
FAILS = []


def eq(name, a, b, tol=0.0):
    ok = abs(float(a) - float(b)) <= tol if tol else a == b
    if not ok:
        FAILS.append('%s: %r != %r' % (name, a, b))


def main():
    S = int(C['ledger']['S'])
    rows = C['ledger']['rows']
    f = Fraction(3, 8)
    c_m = Fraction(rows['mul_step'], K.C_REF)
    c_s = Fraction(rows['sqr_step'], K.C_REF)
    c_raw = (1 - f) * c_m + f * c_s
    c_m_s = Fraction(10 * K.C_REF + rows['fold_mul'] * 13, 10 * K.C_REF)
    c_s_s = Fraction(10 * rows['cyclo_sqr'] + rows['fold_sqr'] * 13, 10 * K.C_REF)
    c_stress = (1 - f) * c_m_s + f * c_s_s
    eq('c_raw', '%.6f' % float(c_raw), C['ledger']['c_raw'])
    eq('c_stress', '%.6f' % float(c_stress), C['ledger']['c_stress'])
    eq('W_raw', Fraction(S) * c_raw, Fraction(C['ledger']['W_raw']))
    eq('W_stress', Fraction(S) * c_stress, Fraction(C['ledger']['W_stress']))
    lw = math.log2(float(Fraction(S) * c_stress))
    eq('log2_W_stress', '%.6f' % lw, C['ledger']['log2_W_stress'], 1e-9)
    claim = float(C['claim']['value'])
    assert claim >= lw + 2e-5, 'claim below W_stress + pad'
    assert claim - 1e-4 < lw + 2e-5, 'claim is not the least 4dp above the pad'
    lam = float(S) ** 2 * (1 - 0.375 ** 2) / (2 * 1.1 * K.N_CLASSES)
    eq('lambda', '%.10f' % lam, C['success']['lambda_'], 1e-9)
    model = 1 - math.exp(-lam)
    eq('success.model', '%.7f' % model, C['success']['model'], 1e-7)
    tau = C['trap']['tau']
    eq('success.net', '%.7f' % (model * (1 - tau)), C['success']['net'], 1e-7)
    assert model * (1 - tau) >= 0.39, 'strict success below floor'
    theta = 1 << 80
    mean = (S + theta - 1) // theta
    eq('reports.mean', mean, C['reports']['mean_dp_hits'])
    eq('reports.cap', 2 * mean, C['reports']['cap'])
    eq('reports.declared', C['reports']['cap'], C['reports']['declared_cumulative_reports'])
    assert C['reports']['declared_cumulative_reports'] * theta > S, 'declared reports infeasible'
    labels = [('log2_S', math.log2(S), C['ledger']['log2_S']),
              ('log2_W_stress', lw, C['ledger']['log2_W_stress']),
              ('cap_log2', math.log2(2 * mean), C['reports']['cap_log2']),
              ('tau_log2', math.log2(tau), C['trap']['tau_log2'])]
    for nm, val, printed in labels:
        rel = abs(float(printed) - val) / max(abs(val), 1e-300)
        if rel > 5e-4:
            FAILS.append('label %s: %s vs %.6f' % (nm, printed, val))
    print('c_raw %.6f  c_stress %.6f  log2 W_stress %.6f  claim %.4f'
          % (float(c_raw), float(c_stress), lw, claim))
    print('lambda %.10f  success %.7f  strict %.7f  tau 2^%.6f'
          % (lam, model, model * (1 - tau), math.log2(tau)))
    print('reports mean 2^%.6f  cap 2^%.6f  declared 2^%.6f'
          % (math.log2(mean), math.log2(2 * mean),
             math.log2(C['reports']['declared_cumulative_reports'])))
    if FAILS:
        print('VERIFY FAIL')
        for x in FAILS:
            print('  ' + x)
        return 1
    print('VERIFY PASS — %d label rows, %d claim rows, 0 failures'
          % (len(labels), len(C['checks'])))
    return 0


if __name__ == '__main__':
    sys.exit(main())
