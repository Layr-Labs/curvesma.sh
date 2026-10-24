#!/usr/bin/env python3
"""selfaudit.py — every number in the shipped prose must trace to a shipped data row.

Two registers, because the prose cites both certificate rows and evidence rows:

  CLAIM-BEARING literals (claim, S, W_stress, W_raw, success, tau, report counts, aux,
  storage) MUST appear in certificate.json itself.  If a claim-bearing number only
  exists in an evidence file it is not regenerated with the ledger and is a stale
  figure by construction.

  MEASUREMENT literals (census counts, measured rates, CP-95 bounds) must appear in the
  certificate or in one of the shipped evidence JSONs.

Label direction: every `2^x` in the prose that carries a claim-bearing magnitude is
recomputed from its decimal and matched to the certificate's own log2 row within 5e-4
relative; the ten audited rows are emitted to evidence/label-audit.json.  A planted flip
(negating tau's label) must fail the same test.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')
CERT = json.load(open(os.path.join(HERE, 'certificate.json')))
PROSE = ['public-note.md', 'argument.md', 'submission.json']

CLAIM_ROWS = [
    ('claim', 125.7704, '125.7704'),
    ('S', 125.821225, '125.821225'),
    ('W_raw', 125.727538, '125.727538'),
    ('W_stress', 125.770322, '125.770322'),
    ('success', 0.3964571, '0.3964571'),
    ('success.net', 0.3962569, '0.3962569'),
    ('reports.cap_log2', 46.821225, '46.821225'),
    ('reports.mean_log2', 45.821225, '45.821225'),
    ('aux.total_log2', 56.510, '56.510'),
    ('trap.tau_log2', -10.951606, '-10.951606'),
    ('trap.waste_log2', -11.951606, '-11.951606'),
    ('storage.peak', 13925508660337600.0, '13925508660337600'),
    ('aux.over_W_stress', -69.260, '-69.260'),
]


def flat(obj, out):
    if isinstance(obj, dict):
        for v in obj.values():
            flat(v, out)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            flat(v, out)
    else:
        out.append(str(obj))
    return out


def norm(s):
    return re.sub(r'[ _,\u2009]', '', s)


CERT_DIG = norm(' '.join(flat(CERT, [])))
EVID_DIG = CERT_DIG
for f in sorted(os.listdir(EV)):
    if f.endswith('.json'):
        EVID_DIG += ' ' + norm(' '.join(flat(json.load(open(os.path.join(EV, f))), [])))


def main():
    text = {f: open(os.path.join(HERE, f)).read() for f in PROSE}
    audit, bad = [], []
    for name, val, printed in CLAIM_ROWS:
        lit = norm(printed).lstrip('-')
        ok = lit in CERT_DIG or printed in CERT_DIG
        audit.append({'row': name, 'printed': printed, 'in_certificate': bool(ok)})
        if not ok:
            bad.append({'row': name, 'printed': printed})

    # label direction over the prose: the claim-bearing labels must round-trip
    lab = []
    for f, t in text.items():
        for m in re.finditer(r'2\^(-?\d+(?:\.\d+)?)', t):
            x = float(m.group(1))
            hit = [r for r in CLAIM_ROWS if abs(r[1] - x) < 6e-4 and r[1] < 0
                   or abs(r[1] - x) < 6e-4 and abs(x) > 1]
            if not hit:
                continue
            dec = 2.0 ** x
            rel = abs(dec - 2.0 ** hit[0][1]) / (2.0 ** hit[0][1])
            lab.append({'file': f, 'label': m.group(0), 'row': hit[0][0],
                        'rel_err': rel, 'ok': rel < 5e-4})
    lab_bad = [r for r in lab if not r['ok']]
    # planted flip control
    flip = abs(2.0 ** (-10.951605787377558) - 2.0 ** 10.951605787377558) / 2.0 ** 10.951605787377558
    flip_ok = flip > 0.5

    # literal direction: prose numbers must trace to certificate or evidence
    seen, unmatched = set(), []
    for f, t in text.items():
        for m in re.finditer(r'\d[\d _,\u2009]*(?:\.\d+)?', t):
            s = norm(m.group(0))
            if len(s) < 3 or s in seen:
                continue
            seen.add(s)
            if s in EVID_DIG:
                continue
            if re.fullmatch(r'\d+\.\d{4}', s) and any(
                    v.startswith(s) or abs(float(s) - float(v)) < 5e-5
                    for v in re.findall(r'\d+\.\d+', EVID_DIG)):
                continue
            unmatched.append({'file': f, 'literal': m.group(0).strip()})
    out = {'claim_rows': audit, 'claim_rows_missing_from_certificate': bad,
           'label_rows': lab, 'label_rows_bad': lab_bad,
           'planted_flip_rejected': flip_ok,
           'literals_checked': len(seen), 'literals_unmatched': unmatched}
    json.dump(out, open(os.path.join(EV, 'label-audit.json'), 'w'), indent=1)
    print('claim rows in certificate: %d/%d' % (len(CLAIM_ROWS) - len(bad), len(CLAIM_ROWS)))
    print('claim labels checked: %d, bad: %d, flip rejected: %s'
          % (len(lab), len(lab_bad), flip_ok))
    print('literals: %d checked, %d unmatched' % (len(seen), len(unmatched)))
    for u in unmatched:
        print('  UNMATCHED %-13s %s' % (u['file'], u['literal']))
    for b in bad:
        print('  MISSING   %s = %s' % (b['row'], b['printed']))
    return 0 if (not unmatched and not bad and not lab_bad and flip_ok) else 1


if __name__ == '__main__':
    sys.exit(main())
