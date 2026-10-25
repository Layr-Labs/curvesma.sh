#!/usr/bin/env python3
"""v31 ledger.  The four v30 flaws repaired where they are computed.

FLAW 1  RETENTION TAIL.  v30 charged an asserted exponential envelope evaluated
        at 2^-3.589688 into a claim whose error budget is 2^-100-class, and
        extrapolated the fitted decay past its support.  v31 (a) runs the merge
        simulator over the FULL production range (cap_mult 128, budget 8x the old
        support) and integrates the raw ages directly, (b) closes the only
        unobserved window (a_hi, 128 theta] with a residual whose rate is FITTED
        to the raw tail decades with residuals shipped, and (c) does not use the
        fixture law as a production envelope at all: the production charge is the
        pair-uniform age law + Markov in the report distance, delta = 2/L = 2^-6,
        a COMPUTED row that is 155 bits larger than v30's fixture integral.

FLAW 2  LAUNCH PROTOCOL + CUMULATIVE ABANDONED WORK.  The gate now lives in the
        shipped loop (v31_protocol.RULE 3) and this module imports the constant
        and the halt string from it; there is no separate harness protocol.  The
        abandoned-trail work is a row of its own (abandon_row): it supersedes v30's
        trap_charge and is the cumulative charge n_inflight*cap plus the
        pre-halt population n_launches * P(collide inside cap) * cap, both priced
        at the max single-step price.  W = w_effective + abandon_row + aux +
        pairing + overflow + stop_loss.

FLAW 3  127 CONTROL.  See audit_v31: the shipped protocol is run at max launch
        rate until in-flight EXACTLY HITS 128 and the 129th launch is refused;
        the 127-trail reserve is then charged against that observed 128 and is
        shown to under-cover by exactly one cap.

FLAW 4  PAIRING / STORAGE.  storage31() prices the fixed-P pairing base, the
        returned G2 point, the pinned pairing check and the 128 LIVE GT states
        (v30's inherited storage() defaulted live-trail storage to zero), and
        pairing_production_width() audits every record field for production-width
        representability: multiword counts, full-width refs with an explicit valid
        bit (never a zero-word test), and no 64-bit overflow at production sizes.
"""
from fractions import Fraction
import json
import math

import constants as K
import ledger as L
import ledger26 as L26
import ledger27 as L27
import v31_protocol as P
import v31_tailfit as TF

THETA = 1 << K.THETA_BITS
CAP_MULT = 128
N_INFLIGHT = P.INFLIGHT_CAP
AGE_PATH = 'evidence/v31-tail-raw.json'
MOM_PATH = 'evidence/v30-moments.json'


def _load(path):
    with open(path) as fh:
        return json.load(fh)


def delta_closed_form(cap_mult=CAP_MULT):
    m = float(cap_mult)
    e = math.exp(-m)
    return e * (m ** 3 / 3.0 + m ** 2) / (2.0 * (1.0 - e * (1.0 + m)))


def delta_charge(path=AGE_PATH):
    """Production pair-uniform charge, floored by the closed form and the
    fixture routes; every term computed."""
    p = TF.production_delta()
    d = dict(production=p, fixture=None, analytic=delta_closed_form())
    d['analytic_log2'] = math.log2(d['analytic'])
    terms = {'production_markov': p['charged'], 'closed_form': d['analytic']}
    try:
        raw = _load(path)
        f = TF.fit(raw)
        d['fixture'] = f
        terms['fixture_measured_plus_residual'] = f['delta_fixture']
    except Exception as exc:                     # evidence not yet written
        d['fixture_error'] = repr(exc)
    name = max(terms, key=lambda k: terms[k])
    d.update(terms=terms, charged=terms[name], charged_term=name,
             charged_log2=math.log2(terms[name]),
             rule='max(production pair-uniform 2^-6, closed form, fixture '
                  'measured+slope-fitted residual); no asserted envelope')
    return terms[name], d


def delta_fixture_only(path=AGE_PATH, truncate=None):
    return TF.fit(_load(path), truncate=truncate)


def abandon_row(s=None, cap_mult=CAP_MULT, n_inflight=N_INFLIGHT,
                reserve_mult=1, c_mul=None, n=None, theta=THETA):
    """CUMULATIVE abandoned-trail work, priced at the max single-step price.

    (i)  post-halt in-flight population: the shipped gate holds at most
         n_inflight trails live; at the halt each may have spent up to its cap,
         so <= n_inflight*cap*reserve_mult steps are abandoned after the halt.
    (ii) pre-halt population: a trail that never collides inside its cap is
         abandoned at the cap; over n_launches = ceil(S/theta) launches this is
         n_launches * P(collide inside cap) * cap, with P = cap^2/(2N) from the
         same pair-uniform law the success model uses.
    Both are summed; neither is inside w_effective, which prices only the S steps
    the halt rule admits.
    """
    s = S31() if s is None else s
    n = K.N_CLASSES if n is None else n
    c_mul = Fraction(K.C_REF + L.rows()['fold_mul'] * K.STRESS_FOLD, K.C_REF) \
        if c_mul is None else c_mul
    cap = cap_mult * theta
    n_launch = -(-int(s) // int(theta))
    p_cap = Fraction(cap * cap, 2 * int(n))
    halt_steps = int(n_inflight) * int(cap) * int(reserve_mult)
    pre_steps = Fraction(n_launch) * p_cap * int(cap)
    charge = (Fraction(halt_steps) + pre_steps) * c_mul
    return charge, dict(
        population='cumulative abandoned trails: post-halt in-flight x cap x '
                   'reserve_mult + n_launches x P(cap) x cap',
        n_inflight=int(n_inflight), cap=cap, cap_mult=cap_mult,
        reserve_mult=int(reserve_mult), halt_abandoned_steps=str(halt_steps),
        n_launches=n_launch, p_cap_log2=math.log2(float(p_cap)),
        pre_halt_abandoned_steps=str(float(pre_steps)),
        pre_halt_log2=math.log2(float(pre_steps)) if pre_steps > 0 else None,
        total_steps_log2=math.log2(float(Fraction(halt_steps) + pre_steps)),
        charge_log2=math.log2(float(charge)),
        price='max single-step price c_mul_stress',
        note='supersedes v30 trap_charge (same in-flight bound, now with the '
             'cumulative pre-halt population added)')


def pairing_row(c_mul=None):
    """Work of the final pinned check: fixed-P pairing base, G2 reconstruction,
    one Miller loop + final exponentiation.  Priced in Fp12-mul equivalents."""
    c_mul = Fraction(K.C_REF + L.rows()['fold_mul'] * K.STRESS_FOLD, K.C_REF) \
        if c_mul is None else c_mul
    base = Fraction(2 ** 12)          # 68 fixed-P line pairs, Fp2 muls
    g2 = Fraction(2 ** 11)            # x -> affine G2 point (Fp2 sqrt chain)
    miller = Fraction(2 ** 14)        # Miller loop over Q, Fp12 muls
    final = Fraction(2 ** 13)         # final exponentiation
    steps = base + g2 + miller + final
    charge = steps * c_mul
    return charge, dict(components=dict(pairing_base=str(base), g2_point=str(g2),
                                        miller_loop=str(miller),
                                        final_exp=str(final)),
                        steps=str(steps), charge_log2=math.log2(float(charge)),
                        rule='one pinned e(P,Q)=z check after the 12-shift '
                             'solve; the solve itself is charged in aux.solve')


def storage31(s=None, theta=THETA, n_inflight=N_INFLIGHT):
    """PEAK STORAGE at production width, with the live GT states charged.

    v30 inherited storage() with live-trail storage defaulted to zero.  Here the
    128 live trails each hold one GT element (y, 12 x 48 B = 576 B), the two
    80-bit labels, an own-step counter and a report position; the fixed-P pairing
    base and the returned G2 point are charged; the record store is the
    deterministic cap with a 2x slack.
    """
    s = S31() if s is None else s
    reports = -(-int(s) // int(theta)) * 2
    gt = K.GT_BYTES if hasattr(K, 'GT_BYTES') else 576
    rec = K.RECORD_BYTES
    live = int(n_inflight) * (gt + 2 * K.LABEL_BYTES + 8 + 8)
    u_bytes = (1 << K.TABLE_BITS) * K.TABLE_COLS * K.LABEL_BYTES
    advice = (1 << K.TABLE_BITS) * K.TABLE_COLS * gt
    pair_base = 68 * 2 * 96                      # line coefficient pairs, Fp2
    g2_point = 2 * 96
    target = gt
    code = 262144
    scratch = 8 * gt
    peak = (reports * rec + live + u_bytes + advice + pair_base + g2_point
            + target + code + scratch)
    return dict(records=reports, store_bytes=reports * rec,
                live_gt_states=int(n_inflight), live_trail_bytes=live,
                u_label_bytes=u_bytes, advice_bytes=advice,
                pairing_base_bytes=pair_base, g2_point_bytes=g2_point,
                pinned_target_bytes=target, code_bytes=code,
                scratch_bytes=scratch, peak_bytes=peak,
                peak_log2=math.log2(peak),
                overflow_policy='drop-oldest + drain (charged in aux); capacity '
                                'is the deterministic cap')


def pairing_production_width(s=None, theta=THETA, n_inflight=N_INFLIGHT):
    """PRODUCTION-WIDTH REPRESENTABILITY AUDIT of the record, ref and count types.

    Fleet failure class (g2-invert v36): a calloc'd count overflowed SIZE_MAX and
    a reference was tested by its low word being zero, which aliases a live slot.
    Every field below is checked at the production size S = 2^126.67.
    """
    s = S31() if s is None else s
    reports = -(-int(s) // int(theta)) * 2
    rec = 112
    count_bits = reports.bit_length()
    offset_bits = (reports * rec).bit_length()
    words = (K.KAPPA_BITS + 63) // 64
    gt_words = (12 * 381 + 63) // 64
    rows = dict(
        record_count=dict(value=reports, bits=count_bits, width=64,
                          fits_u64=count_bits <= 64, repr='2-limb multiword '
                          '(hi,lo) with carry assert', sharded=True),
        byte_offset=dict(value=reports * rec, bits=offset_bits, width=64,
                         fits_u64=offset_bits <= 64, repr='(shard_id u32, '
                         'slot u32, limb offset 2x64)'),
        class_key=dict(words=words, bits=K.KAPPA_BITS, truncation='none before '
                       'mixing', repr='%d x u64 little-endian' % words),
        gt_element=dict(words=gt_words, bits=12 * 381,
                        repr='12 Fp limbs x 6 u64 = 576 B, full width'),
        ref=dict(form='(valid:u1, shard:u32, slot:u32)', zero_word_test=False,
                 note='refs carry an explicit valid bit; a slot is never '
                      'identified by a zero word (the aliasing mode that killed '
                      'g2-invert v36)'),
        labels=dict(bits=K.THETA_BITS, words=2,
                    repr='(a,b) 2 x 2 u64, no truncation; every branch input is a '
                         'function of the class key alone'),
        size_max=dict(record_bytes=reports * rec, size_max=2 ** 64 - 1,
                      headroom_log2=math.log2(float(2 ** 64 - 1) / (reports * rec)),
                      no_overflow=(reports * rec < 2 ** 64)),
        shards=dict(count=-(-reports // (1 << 32)), slots_per_shard=1 << 32,
                    bytes_per_shard=(1 << 32) * rec),
    )
    return rows


def production_f31(path='evidence/v30-census.json'):
    d = _load(path)
    sq = sum(d[k]['sq_steps'] for k in ('1', '2', '3', '4'))
    st = sum(d[k]['steps'] for k in ('1', '2', '3', '4'))
    return sq / st


def S31(delta=None):
    d = delta if delta is not None else delta_charge()[0]
    return L27.step_budget27(delta=d)


def aux31(s=None, theta=THETA):
    s = S31() if s is None else s
    reports = -(-int(s) // int(theta))
    rows, _ = L26.aux26(s=s, theta=theta, slack=2)
    return rows, reports


def stop_loss31(s=None, theta=THETA, cap_mult=CAP_MULT, c_stress=None,
                n_inflight=N_INFLIGHT):
    s = S31() if s is None else s
    c_stress = Fraction(K.C_REF + L.rows()['fold_mul'] * K.STRESS_FOLD, K.C_REF) \
        if c_stress is None else c_stress
    cap = cap_mult * theta
    n_inflight = int(n_inflight)
    reserve_steps = n_inflight * cap
    rel = Fraction(2 * reserve_steps, int(s))
    charge = Fraction(2 * reserve_steps) * c_stress
    return charge, dict(
        n_inflight=n_inflight, launch_interval_steps=theta, trail_cap_steps=cap,
        reserve_steps=str(reserve_steps), reserve_over_S=reserve_steps / float(s),
        loss_rel_S2=float(rel),
        loss_bits=-math.log2(1.0 - float(rel)) if rel < 1 else None,
        price='max single-step price (c_mul_stress)', halt_rule=P.HALT_RULE,
        rule='halt is a rule of the shipped loop (v31_protocol.RULE 3 + HALT_RULE)')


def total31(f_measured=None, delta=None, reserve_mult=1):
    if f_measured is None:
        f_measured = production_f31()
    d, ddetail = delta_charge() if delta is None else (delta, dict(charged=delta))
    s = S31(d)
    Lg = L26.ledgers26(f_measured=f_measured, s=s)
    w = Fraction(int(s)) * Lg['c_mul_stress']
    ab, abd = abandon_row(s=s, c_mul=Lg['c_mul_stress'], reserve_mult=reserve_mult)
    rows, reports = aux31(s=s)
    aux_sum = sum(rows.values(), Fraction(0))
    pr, prd = pairing_row(c_mul=Lg['c_mul_stress'])
    overflow_row = Fraction(1)
    sl, sld = stop_loss31(s=s, c_stress=Lg['c_mul_stress'])
    comps = dict(w_effective=w, abandon_row=ab, aux_sum=aux_sum,
                 pairing_row=pr, overflow_row=overflow_row, stop_loss_row=sl)
    total = sum(comps.values(), Fraction(0))
    return dict(components=comps, total=total, log2_total=math.log2(float(total)),
                abandon=abd, aux_rows={k: str(v) for k, v in rows.items()},
                reports=reports, ledgers=Lg, S=s, f=Lg['f'], delta=ddetail,
                pairing=prd, stop=sld, storage=storage31(s=s),
                width=pairing_production_width(s=s),
                diagnostic_average_price=dict(
                    c_stress=str(Lg['c_stress']),
                    w_stress_at_average=float(s) * float(Lg['c_stress']),
                    w_stress_at_average_log2=math.log2(float(s) * float(Lg['c_stress'])),
                    note='mean-price row shown only to expose the undercount; NOT summed'))


def claim31(dp=4, f_measured=None, delta=None, reserve_mult=1):
    T = total31(f_measured=f_measured, delta=delta, reserve_mult=reserve_mult)
    lo = math.floor(T['log2_total'] * 10 ** dp) / 10 ** dp
    while lo < T['log2_total']:
        lo += 10.0 ** (-dp)
    lo = round(lo, dp)
    num = int(round(lo * 10 ** dp))
    den = 10 ** dp
    from math import gcd
    g = gcd(num, den)
    return dict(claim_log2=lo, numerator=num // g, denominator=den // g,
                total_log2=T['log2_total'], total=str(T['total']), T=T)


if __name__ == '__main__':
    dm, ddetail = delta_charge()
    c = claim31()
    T = c['T']
    print('delta charged 2^%.6f (%s)  production 2^%.4f exact 2^%.4f'
          % (ddetail['charged_log2'], ddetail['charged_term'],
             ddetail['production']['delta_markov_log2'],
             ddetail['production']['delta_exact_log2']))
    if ddetail.get('fixture'):
        f = ddetail['fixture']
        print('fixture: window [%.1f,%.1f] theta kappa %.4f +/- %.4f resid max '
              '%.3f bits  I_meas 2^%.2f residual 2^%s'
              % (f['u_lo'], f['u_hi'], f['kappa_per_theta'], f['kappa_se'],
                 f['resid_max_bits'], f['integral_measured_log2'],
                 f['residual_log2']))
    print('S = %d  2^%.6f  f_priced %s' % (T['S'], math.log2(T['S']), T['f']))
    for k, v in T['components'].items():
        print('  %-16s 2^%.6f' % (k, math.log2(float(v))))
    print('total 2^%.9f  claim %.4f = %d/%d'
          % (c['total_log2'], c['claim_log2'], c['numerator'], c['denominator']))
    print('storage peak 2^%.3f B  live trails %d x (%d B) = %d B'
          % (T['storage']['peak_log2'], T['storage']['live_gt_states'],
             576 + 2 * K.LABEL_BYTES + 16, T['storage']['live_trail_bytes']))
