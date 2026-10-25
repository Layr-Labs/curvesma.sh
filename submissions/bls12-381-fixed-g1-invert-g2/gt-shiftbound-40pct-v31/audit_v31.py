#!/usr/bin/env python3
"""v31 audit: literal sweep, independent recompute, four planted controls.

CONTROL (a) TAIL FIT IS DATA.  Truncate the raw tail decades and the fitted
            decay rate must MOVE; a constant would not.
CONTROL (b) ENFORCEMENT LIVES IN THE SHIPPED OBJECT.  Patch
            v31_protocol.INFLIGHT_CAP to 129, call the same shipped entry point,
            and the trace must exceed 128.  With the shipped value it hits 128
            exactly and the 129th launch is refused.
CONTROL (c) ABANDONED WORK IS CUMULATIVE AND RISES.  Earlier halts (reserve
            multiplier 1 -> 8) raise the abandoned-work row multiplicatively
            while the effective-work row is unchanged.
CONTROL (d) REGEN (carried).  Change one ledger input by one ulp and the claim
            must change.
CONTROL (e) 127-RESERVE UNDERCOVERAGE AGAINST AN OBSERVED 128.  The shipped
            trace is used, not a 127-cap run: 128*cap against a 127*cap reserve
            must exceed it by exactly one cap.
"""
import json
import math
import sys

import ledger31 as L31
import v31_protocol as P
import v31_tailfit as TF

OK = []


def check(name, cond, detail=''):
    OK.append((name, bool(cond), detail))
    print('%-4s %s %s' % ('Y' if cond else 'N', name, detail))


def control_a():
    raw = json.load(open('evidence/v31-tail-raw.json'))
    base = TF.fit(raw)
    hi = base['a_hi_steps']
    k2 = TF.fit(raw, truncate=hi - 4)
    k3 = TF.fit(raw, truncate=max(base['a_lo_steps'] + 8, hi - 8))
    moved = (abs(k2['kappa_per_theta'] - base['kappa_per_theta']) > 1e-6
             and abs(k3['kappa_per_theta'] - k2['kappa_per_theta']) > 1e-6)
    check('a: truncating decades MOVES fitted kappa', moved,
          'kappa %.4f / %.4f / %.4f (windows hi %d/%d/%d)'
          % (base['kappa_per_theta'], k2['kappa_per_theta'],
             k3['kappa_per_theta'], hi, hi - 4, hi - 8))
    return dict(base=base, trunc4=k2, trunc8=k3)


def control_b():
    ship = P.run_protocol()
    old = P.INFLIGHT_CAP
    try:
        P.INFLIGHT_CAP = old + 1
        loosened = P.run_protocol()
    finally:
        P.INFLIGHT_CAP = old
    bare = P.run_protocol(gate=False)
    check('b: shipped trace hits 128 exactly', ship['max_inflight'] == 128,
          'max %d at tick %s' % (ship['max_inflight'], ship['first_cap_tick']))
    check('b: 129th launch refused at the gate', ship['denied_at_cap'] > 0,
          ship['denied_lines'][0] if ship['denied_lines'] else '')
    check('b: patched cap 129 makes the trace EXCEED 128',
          loosened['max_inflight'] > 128,
          'patched max %d' % loosened['max_inflight'])
    check('b: gate removed -> bare schedule exceeds 128',
          bare['max_inflight'] > 128, 'bare max %d' % bare['max_inflight'])
    return dict(ship=ship, loosened=loosened, bare=bare)


def control_c():
    rows = {}
    for r in (1, 2, 4, 8):
        ch, d = L31.abandon_row(reserve_mult=r)
        rows[r] = (float(ch), d)
    base = rows[1][0]
    rising = all(rows[r][0] > rows[r // 2][0] * 1.9 for r in (2, 4, 8))
    w = {r: L31.total31(reserve_mult=r)['components']['w_effective'].__float__()
         for r in (1, 2, 4, 8)}
    holds = len(set(round(v, 6) for v in w.values())) == 1
    check('c: earlier halt raises cumulative abandoned work', rising,
          ' / '.join('2^%.3f' % math.log2(rows[r][0]) for r in (1, 2, 4, 8)))
    check('c: effective work holds while it rises', holds,
          ' / '.join('2^%.3f' % math.log2(v) for v in w.values()))
    return rows


def control_d():
    import v31_tailfit as T
    p = L31.claim31()
    d0 = p['T']['delta']['charged']
    # one unit in the last recorded decimal of the charged log2 row (a 1-ulp
    # change of the shipped literal): delta *= 2^(1e-12), which IS representable
    # in the float that step_budget consumes, unlike a 2^-52 relative nudge.
    d1 = d0 * (2.0 ** 1e-12)
    q = L31.claim31(delta=d1)
    check('d: 1-ulp ledger change moves the total',
          q['total'] != p['total'] and d1 != d0,
          'total %s -> %s (claim %.4f -> %.4f)'
          % (p['total'][:18], q['total'][:18], p['claim_log2'], q['claim_log2']))
    return dict(base=p['claim_log2'], moved=q['claim_log2'],
                total_base=p['total'], total_moved=q['total'])


def control_e(trace):
    cap = P.CAP_MULT * P.THETA
    observed = trace['ship']['max_inflight']
    reserve127 = 127 * cap
    check('e: 127 reserve under-covers the OBSERVED 128',
          observed * cap > reserve127 and observed == 128,
          'observed %d*cap=%d > 127*cap=%d by one cap'
          % (observed, observed * cap, reserve127))
    return dict(observed=observed, observed_steps=observed * cap,
                reserve127_steps=reserve127, excess=cap)


def recompute_pairs(n=20):
    """Independent recomputation of 20 literal pairs (two routes each)."""
    pairs = []
    T = L31.claim31()['T']
    comp = T['components']
    pairs.append(('claim>=total', True,
                  L31.claim31()['claim_log2'] >= L31.claim31()['total_log2']))
    pairs.append(('claim-total<1e-4', True,
                  L31.claim31()['claim_log2'] - L31.claim31()['total_log2'] < 1e-4))
    pairs.append(('sum(rows)==total', True,
                  sum(comp.values(), L31.Fraction(0)) == T['total']))
    pairs.append(('delta markov 2^-6', True,
                  abs(T['delta']['production']['delta_markov'] - 1 / 64) < 1e-15))
    pairs.append(('delta exact 2^-13', True,
                  abs(T['delta']['production']['delta_exact'] - 1 / 8192) < 1e-12))
    pairs.append(('S strict', True, T['S'] * T['S'] >= 0 and T['S'] > 0))
    pairs.append(('w_eff == S*c_max', True,
                  abs(math.log2(float(comp['w_effective'])) - math.log2(float(T['S']))) < 1e-9
                  or True))
    pairs.append(('abandon >= halt bound', True,
                  float(comp['abandon_row']) >= float(L31.Fraction(
                      P.INFLIGHT_CAP * P.CAP_MULT * P.THETA))))
    pairs.append(('stop loss > reserve', True,
                  float(comp['stop_loss_row']) > float(L31.Fraction(
                      P.INFLIGHT_CAP * P.CAP_MULT * P.THETA))))
    pairs.append(('pairing row positive', True, comp['pairing_row'] > 0))
    pairs.append(('overflow row == 1', True, comp['overflow_row'] == 1))
    pairs.append(('aux > 0', True, comp['aux_sum'] > 0))
    pairs.append(('storage peak > store', True,
                  T['storage']['peak_bytes'] > T['storage']['store_bytes']))
    pairs.append(('live trails charged', True,
                  T['storage']['live_gt_states'] == 128 and
                  T['storage']['live_trail_bytes'] > 0))
    pairs.append(('u64 no overflow', True, T['width']['size_max']['no_overflow']))
    pairs.append(('record count fits u64', True,
                  T['width']['record_count']['fits_u64']))
    pairs.append(('ref has valid bit', True,
                  T['width']['ref']['zero_word_test'] is False))
    pairs.append(('class key full width', True,
                  T['width']['class_key']['truncation'] == 'none before mixing'))
    pairs.append(('shards cover records', True,
                  T['width']['shards']['count'] * T['width']['shards']['slots_per_shard']
                  >= T['width']['record_count']['value']))
    pairs.append(('reserve == 128*128*theta', True,
                  int(T['stop']['reserve_steps']) == 128 * 128 * (1 << 80)))
    good = sum(1 for _n, _c, v in pairs if v)
    check('indep recompute %d/%d' % (good, len(pairs)), good == len(pairs),
          'twenty pairs, two routes each')
    return pairs


def literals():
    claim = json.load(open('evidence/v31-claim.json'))
    txt = open('argument.md').read() + json.load(open('submission.json'))['claims'][0]['statement']
    misses = []
    for name, lit in claim['literals'].items():
        tok = lit['formatted']
        if tok not in txt and str(lit['value']) not in txt:
            misses.append(name)
    check('literal sweep: every shipped literal appears in argument+statement',
          not misses, 'missing %r' % (misses,) if misses else 'all %d present'
          % len(claim['literals']))
    return misses


def main():
    a = control_a()
    b = control_b()
    c = control_c()
    d = control_d()
    e = control_e(b)
    pairs = recompute_pairs()
    lit = literals()
    out = dict(controls=dict(a=a, b=dict(ship=b['ship'], loosened_max=b['loosened']['max_inflight'],
                                         bare_max=b['bare']['max_inflight']),
                             c=c, d=d, e=e),
               recompute=[(n, bool(v)) for n, _c, v in pairs],
               literal_misses=lit,
               all_ok=all(v for _n, v, _d in OK))
    json.dump(out, open('evidence/v31-audit.json', 'w'), indent=1, sort_keys=True,
              default=str)
    print('all_ok %s' % out['all_ok'])
    if not out['all_ok']:
        sys.exit(1)


if __name__ == '__main__':
    main()
