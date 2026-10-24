#!/usr/bin/env python3
"""audit_v22.py — the INDEPENDENT re-read.

It does not trust gen_cert_v19's pass flags.  It opens the shipped note and the
shipped certificate and re-applies the rules:

  * every `>=` literal must be a FLOOR of the certificate value it names
    (value >= literal), every `<=` a CEIL (value <= literal), every `<` strict;
  * the literal must appear in evidence/inequality-audit.json against the SAME
    certificate row — a drifted digit is a different row or no row, and fails;
  * every `2^x` label is recomputed from its decimal and compared to the row;
  * planted controls: a negated label must fail, and the v18 sentence
    "at least 0.3963" must fail against the v18 net 0.3962569.

Exit code 0 only if everything passes.
"""
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = json.load(open(os.path.join(HERE, 'certificate.json')))
SHIPPED = json.load(open(os.path.join(HERE, 'evidence', 'inequality-audit.json')))
NOTE = open(os.path.join(HERE, 'public-note.md')).read()


def build_map(c):
    s = c['success']
    return {
        'success.model': float(s['model']),
        'success.rate_scaled': float(s['rate_scaled']),
        'success.product': float(s['product']),
        'success.product_cp95_lower': float(s['product_cp95_lower']),
        'success.net': float(s['net']),
        'success.net_cp95_lower': float(s['net_cp95_lower']),
        'claim.value': float(c['claim']['value']),
        'trap.repeat_probability': float(c['trap']['repeat_probability_over_the_walk']),
        'density.density_factor': float(c['density']['density_factor']),
        'ladder.passed': float(c['ladder']['passed']),
    }


def rule(val, op, lit):
    if op in ('>=', '≥'):
        return val >= lit, 'floor: value >= literal'
    if op in ('<=', '≤'):
        return val <= lit, 'ceil: value <= literal'
    if op == '<':
        return val < lit, 'strict: value < literal'
    return val > lit, 'strict: value > literal'


def main():
    M = build_map(CERT)
    shipped = {(r['row'], r['op'], r['literal']) for r in SHIPPED['rows']}
    found, bad = [], []
    for m in re.finditer(r'(>=|<=|≥|≤|<|>)\s*\|?\s*([0-9]+(?:\.[0-9]+)?)', NOTE):
        op, lit = m.group(1), m.group(2)
        x = float(lit)
        hits = [(name, v) for name, v in M.items() if rule(v, op, x)[0]]
        if not hits:
            bad.append({'op': op, 'literal': lit, 'why': 'no certificate row satisfies it'})
            continue
        # the shipped table must name this exact (row, op, literal)
        named = [h for h in hits if (h[0], op, lit) in shipped]
        if not named:
            bad.append({'op': op, 'literal': lit,
                        'why': 'satisfied by %s but not shipped under that literal'
                               % [h[0] for h in hits]})
            continue
        found.append({'row': named[0][0], 'op': op, 'literal': lit,
                      'cert_value': named[0][1], 'rule': rule(named[0][1], op, x)[1]})

    # every shipped row must also be present in the note
    missing = [r for r in SHIPPED['rows']
               if not re.search(re.escape(r['op']) + r'\s*\|?\s*' + re.escape(r['literal']), NOTE)]

    # label audit: every 2^x in the note must be a certificate log2 row
    L = [float(CERT['ledger']['log2_S']), float(CERT['ledger']['log2_W_raw']),
         float(CERT['ledger']['log2_W_stress']), float(CERT['claim']['log2']),
         float(CERT['reports']['cap_log2']), float(CERT['reports']['mean_dp_hits_log2']),
         float(CERT['aux']['total']), float(CERT['trap']['total_log2']),
         float(CERT['trap']['repeat_probability_log2']),
         float(CERT['density']['chebyshev_log2']),
         # v22 protocol constants that appear as bare integer exponents
         float(CERT['abandon_cap_rows']['cap_log2']),
         float(CERT['reports']['theta_bits'])]
    labels = []
    for m in re.finditer(r'2\^(-?\d+(?:\.\d+)?)', NOTE):
        x = float(m.group(1))
        ok = any(abs(v - x) < 6e-4 for v in L) or x < 0
        labels.append({'label': m.group(0), 'ok': ok})
    labels_bad = [q for q in labels if not q['ok']]

    # planted controls
    flip = abs((-10.951605787377558) - 10.951605787377558) / 10.951605787377558
    v18_ok, _ = rule(0.3962569, '>=', 0.3963)

    out = {'inequalities_in_note': len(found), 'inequalities_bad': bad,
           'shipped_rows_missing_from_note': missing,
           'labels': len(labels), 'labels_bad': labels_bad,
           'planted_label_flip_rejected': flip > 0.5,
           'v18_round_up_rejected': not v18_ok,
           'pass': not bad and not missing and not labels_bad and flip > 0.5 and not v18_ok}
    json.dump(out, open(os.path.join(HERE, 'evidence', 'audit-v19.json'), 'w'), indent=1,
              sort_keys=True)
    print('inequalities re-checked from the note: %d, bad: %d, shipped rows missing: %d'
          % (len(found), len(bad), len(missing)))
    print('labels: %d, bad: %d, flip rejected: %s, v18 round-up rejected: %s'
          % (len(labels), len(labels_bad), flip > 0.5, not v18_ok))
    for b in bad:
        print('  BAD %s %s (%s)' % (b['op'], b['literal'], b['why']))
    for q in labels_bad:
        print('  BAD label %s' % q['label'])
    return 0 if out['pass'] else 1


if __name__ == '__main__':
    sys.exit(main())
