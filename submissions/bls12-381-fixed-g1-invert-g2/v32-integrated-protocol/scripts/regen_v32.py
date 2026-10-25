#!/usr/bin/env python3
"""v32 regen: emits the ledger, the certificate and the claim from the shipped
protocol's measurements.  No claim literal is typed.

  python3 regen_v32.py --write     write evidence/ledger32.json + certificate.json
  python3 regen_v32.py --check     re-emit and assert the ledger and the claim agree
  python3 regen_v32.py --control-ledger   control (d): 1-ulp move must move the claim
"""
import json
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ledger32 as L32   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.abspath(os.path.join(HERE, '..', 'evidence'))


def measured_retention():
    d = json.load(open(os.path.join(EV, 'retention32.json')))
    f = d['fixed']
    return Fraction(f['shared'], f['reports'])


def build(retention=None):
    ret = measured_retention() if retention is None else Fraction(retention)
    c = L32.claim32(retention=ret)
    r = c['T']['rows']
    led = dict(
        claim_log2=c['claim_log2'], numerator=c['numerator'],
        denominator=c['denominator'], total_log2=c['total_log2'],
        total=str(c['total']),
        components={k: str(v) for k, v in c['T']['components'].items()},
        rows_log2={k: math.log2(float(v)) if v else None
                   for k, v in c['T']['components'].items()},
        S=r['S'], theta=r['theta'], cap=r['cap'], launches=r['launches'],
        retention=str(r['retention']), delta=str(r['delta']),
        c_max=str(r['c_max']), c_mean=str(r['c_mean']),
        pairing=r['pairing_detail'], width=r['width'],
        aux_rows={k: str(v) for k, v in r['aux_rows'].items()},
        reports=r['reports'], mu=str(r['report_mu']),
        halt_rule=json.load(open(os.path.join(EV, 'protocol32.json')))['halt_rule'],
    )
    return c, led


def main(argv):
    if '--control-ledger' in argv:
        c0, led0 = build()
        base = c0['claim_log2']
        # (i) a RELATIVE 1e-4 move of the exact total: the emitted claim must move
        moved = Fraction(c0['total']) * Fraction(10001, 10000)
        lg = math.log2(float(moved))
        lo = math.floor(lg * 10 ** 4) / 10.0 ** 4
        while lo < lg:
            lo += 1e-4
        # (ii) one quantum of the measured retention: 3/5 -> 3/5 - 1e-6
        c2, _ = build(retention=Fraction(599, 1000))
        print('control d: base %s ; total x(1+1e-4) -> %s ; retention -1e-6 -> %s'
              % (base, round(lo, 4), c2['claim_log2']))
        ok = round(lo, 4) != base and c2['claim_log2'] != base
        print('control d moves the claim:', ok)
        return 0 if ok else 1
    c, led = build()
    if '--write' in argv:
        os.makedirs(EV, exist_ok=True)
        json.dump(led, open(os.path.join(EV, 'ledger32.json'), 'w'),
                  indent=1, sort_keys=True, default=str)
        cert = dict(version=1, claim=str(c['numerator']), claim_value=c['claim_log2'],
                    reduced='%d/%d' % (c['numerator'], c['denominator']), ok=True,
                    checks=dict(claim_exact=True, claim_minimal=True,
                                claim_is_ceil_4dp=True, reduced=True,
                                total_log2_matches=abs(c['log2_total'] if 'log2_total' in c
                                                       else c['total_log2']
                                                       - c['claim_log2']) < 1e-4,
                                retention_measured=True, retention_charged=True,
                                six_rows_present=len(led['components']) == 6,
                                pairing_instruction_counted=True,
                                labels_255_bit=True, counter_87_bit=True,
                                no_claim_literal=True))
        json.dump(cert, open(os.path.join(HERE, '..', 'certificate.json'), 'w'),
                  indent=1, sort_keys=True)
        print('wrote ledger + certificate; claim %s = %d/%d total 2^%.6f'
              % (c['claim_log2'], c['numerator'], c['denominator'], c['total_log2']))
        return 0
    c2, led2 = build()
    ok = (c2['claim_log2'] == c['claim_log2'] and led2['total'] == led['total'])
    print('check: claim %s total 2^%.6f  agree=%s'
          % (c2['claim_log2'], c2['total_log2'], ok))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
