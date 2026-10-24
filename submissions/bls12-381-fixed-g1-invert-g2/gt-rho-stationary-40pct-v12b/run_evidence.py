#!/usr/bin/env python3
"""Run the v12b evidence battery: audits at four fixtures, sweeps, equivariance.

Writes a-<bits>.json, d-<bits>.json and e-equivariance.json in the package dir.
Nothing is retyped by hand: every number in the certificate is read from these.
"""
import json
import os
import sys

D = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, D)
import toy12d as T
import equivariance as E

TRIALS = {5: 1600, 7: 3200, 9: 3200, 11: 400}

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('all', 'audit'):
        for b in (5, 7, 9, 11):
            r = T.audit(b, 900 + b)
            json.dump(r, open(os.path.join(D, 'a-%d.json' % b), 'w'), indent=1)
            print('audit', b, r['audit_A_pass'], r['audit_B_pass'], r['audit_E_pass'],
                  r['control_C_old_rule_splits'])
    if which in ('all', 'sweep'):
        for b in (5, 7, 9, 11):
            r = T.sweep(b, TRIALS[b], 900 + b)
            r['trials_requested'] = TRIALS[b]
            json.dump(r, open(os.path.join(D, 'd-%d.json' % b), 'w'), indent=1)
            print('sweep', b, r['recovery'], 'orbit_censor', r['orbit_censor'],
                  'inv_fail', r['inv_failures'])
    if which in ('all', 'equiv'):
        r = E.run_all()
        json.dump(r, open(os.path.join(D, 'e-equivariance.json'), 'w'), indent=1)
        print('equivariance', r['splits'], 'splits of', r['checks'])
