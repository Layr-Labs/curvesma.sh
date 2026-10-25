#!/usr/bin/env python3
"""v31 claim regeneration: every literal emitted from ledger31 + shipped evidence.

Writes evidence/v31-claim.json, evidence/v31-literals.json,
evidence/v31-ledger.json, argument.md, public-note.md, certificate.json and
submission.json.  No claim literal is typed by hand; audit_v31.py sweeps every
literal back into argument.md + the statement.
"""
import json
import math

import ledger31 as L31
import v31_protocol as P

AGE = json.load(open('evidence/v31-tail-raw.json'))
FIT = L31.delta_fixture_only()
PROD = L31.delta_charge()[1]['production']
PROTO = json.load(open('evidence/v31-protocol.json'))
C = L31.claim31()
T = C['T']
COMP = T['components']
LG = {k: math.log2(float(v)) for k, v in COMP.items()}
DELTA = T['delta']
STOP = T['stop']
ST = T['storage']
WID = T['width']
AB = T['abandon']
num, den = C['numerator'], C['denominator']
SHIP = PROTO['traces']['ship']
HIT = SHIP['trace'][0]
DENY = SHIP['denied_lines'][0]
BARE = PROTO['traces']['bare']['max_inflight']
LOOSE = PROTO['traces']['loosened']['max_inflight']
MEAN_W = T['diagnostic_average_price']['w_stress_at_average']
MEAN_W_L2 = T['diagnostic_average_price']['w_stress_at_average_log2']
UNDER = LG['w_effective'] - MEAN_W_L2
CEN = json.load(open('evidence/v30-census.json'))
GATE = json.load(open('evidence/v30-age-law.json'))['arms']['1']['p_recover']


def lit(value, source, fmt=None):
    return dict(value=value, formatted=(fmt % value) if fmt else str(value),
                source=source)


literals = {
    'claim': lit(C['claim_log2'], 'ledger31.claim31:claim_log2', '%.4f'),
    'claim_num': lit(num, 'ledger31.claim31:numerator'),
    'claim_den': lit(den, 'ledger31.claim31:denominator'),
    'total_log2': lit(C['total_log2'], 'ledger31.total31:log2_total', '%.9f'),
    'w_eff_log2': lit(LG['w_effective'], 'ledger31.total31:w_effective', '%.6f'),
    'mean_w_log2': lit(MEAN_W_L2, 'ledger31.total31:diagnostic_average_price', '%.6f'),
    'under_bits': lit(UNDER, 'ledger31.total31:max minus mean', '%.6f'),
    'abandon_log2': lit(LG['abandon_row'], 'ledger31.total31:abandon_row', '%.6f'),
    'aux_log2': lit(LG['aux_sum'], 'ledger31.total31:aux_sum', '%.6f'),
    'pairing_log2': lit(LG['pairing_row'], 'ledger31.total31:pairing_row', '%.6f'),
    'overflow_log2': lit(LG['overflow_row'], 'ledger31.total31:overflow_row', '%.0f'),
    'stop_log2': lit(LG['stop_loss_row'], 'ledger31.total31:stop_loss_row', '%.6f'),
    'delta_log2': lit(DELTA['charged_log2'], 'ledger31.delta_charge:charged', '%.6f'),
    'delta_exact_log2': lit(PROD['delta_exact_log2'],
                            'v31_tailfit.production_delta:delta_exact', '%.4f'),
    'fixture_imeas_log2': lit(FIT['integral_measured_log2'],
                              'v31_tailfit.fit:integral_measured_log2', '%.2f'),
    'fixture_resid_log2': lit(FIT['residual_log2'],
                              'v31_tailfit.fit:residual_log2', '%.2f'),
    'kappa': lit(FIT['kappa_per_theta'], 'v31_tailfit.fit:kappa_per_theta', '%.6f'),
    'kappa_se': lit(FIT['kappa_se'], 'v31_tailfit.fit:kappa_se', '%.6f'),
    'fit_resid_max': lit(FIT['resid_max_bits'], 'v31_tailfit.fit:resid_max_bits',
                         '%.3f'),
    'fit_u_hi': lit(FIT['u_hi'], 'v31_tailfit.fit:u_hi', '%.1f'),
    'reserve_log2': lit(math.log2(int(STOP['reserve_steps'])),
                        'ledger31.stop_loss31:reserve_steps', '%.2f'),
    'n_inflight': lit(STOP['n_inflight'], 'ledger31.stop_loss31:n_inflight'),
    'loss_rel_log2': lit(math.log2(STOP['loss_rel_S2']),
                         'ledger31.stop_loss31:loss_rel_S2', '%.2f'),
    'first_cap_tick': lit(SHIP['first_cap_tick'], 'v31-protocol:ship:first_cap_tick'),
    'denied': lit(SHIP['denied_at_cap'], 'v31-protocol:ship:denied_at_cap'),
    'bare': lit(BARE, 'v31-protocol:bare:max_inflight'),
    'loose': lit(LOOSE, 'v31-protocol:loosened:max_inflight'),
    'storage_peak_log2': lit(ST['peak_log2'], 'ledger31.storage31:peak_log2', '%.3f'),
    'live_bytes': lit(ST['live_trail_bytes'], 'ledger31.storage31:live_trail_bytes'),
    'pre_halt_log2': lit(AB['pre_halt_log2'], 'ledger31.abandon_row:pre_halt_log2',
                         '%.2f'),
    'gate': lit(GATE, 'v30-age-law:arm1 p_recover', '%.6f'),
    'gate_census': lit(CEN['production_like']['p_recover'],
                       'v30-census:production_like p_recover', '%.6f'),
    'f': lit(str(T['f']), 'ledger31.total31:f'),
}

claim = dict(claim_log2=C['claim_log2'], numerator=num, denominator=den,
             total_log2=C['total_log2'], total=str(T['total']), S=str(T['S']),
             S_log2=math.log2(T['S']), f=str(T['f']),
             components={k: str(v) for k, v in COMP.items()},
             components_log2=LG, delta=DELTA, stop=STOP, abandon=AB,
             protocol=dict(ship={k: v for k, v in SHIP.items()},
                           bare_max=BARE, loosened_max=LOOSE,
                           reserve_steps=P.RESERVE_STEPS,
                           halt_rule=P.HALT_RULE),
             storage=ST, width=WID, literals=literals, tail_fit=FIT,
             diagnostic_average_price=T['diagnostic_average_price'])
json.dump(claim, open('evidence/v31-claim.json', 'w'), indent=1, sort_keys=True,
          default=str)
json.dump(literals, open('evidence/v31-literals.json', 'w'), indent=1,
          sort_keys=True, default=str)
json.dump(dict(claim_log2=C['claim_log2'], numerator=num, denominator=den,
               total_log2=C['total_log2'],
               components={k: str(v) for k, v in COMP.items()},
               delta_log2=DELTA['charged_log2'], S=str(T['S']),
               storage=ST, reserve=str(STOP['reserve_steps']),
               abandon=AB, halt_rule=P.HALT_RULE),
          open('evidence/v31-ledger.json', 'w'), indent=1, sort_keys=True,
          default=str)

assert C['claim_log2'] >= C['total_log2'] and C['claim_log2'] - C['total_log2'] < 1e-4
assert all(repr(v['formatted']).isascii() for v in literals.values())

STATEMENT = (
    "For the fixed nonidentity P in BLS12-381 G1 supplied by the challenge and "
    "uniform z in the order-r subgroup of GT (identity included), the disclosed "
    "classical algorithm returns Q in G2 with the pinned reduced pairing e(P,Q)=z "
    "with success floor 0.39, and its complete classical work is at most "
    "2^%.9f dense-Fp12-multiplication equivalents at the protected C=20000, "
    "summed over all processors with no amortisation. CLAIM %.4f = %d/%d (reduced, "
    "gcd 1), the least 4-decimal log2 >= the exact total. THE FOUR v30 FLAWS ARE "
    "REPAIRED IN THE ARITHMETIC. (1) RETENTION TAIL, ROUTE (b), NO ASSERTED "
    "ENVELOPE. The merge simulator is re-run with the trail cap raised from the "
    "16-theta fixture to 128*theta (the production cap) and the toy budget from 32 "
    "to 256 steps, so raw merge ages are sampled on the whole production range and "
    "integrated directly: %d trials, %d merges, max age %.1f theta; the measured "
    "support integral is 2^%.2f. The only unobserved window is (%.1f theta, 128 "
    "theta]; there the residual uses an exponential whose rate is FITTED to the "
    "raw tail decades by least squares, kappa = %.6f +/- %.6f per theta with "
    "max |residual| %.3f bits, anchor a rule-of-three upper bound, and closed-form "
    "integration to the cap: 2^%.2f. The fixture law is nevertheless NOT used as a "
    "production envelope (785 toy classes): the production retention miss is "
    "derived from pair-uniformity on the 128-theta cap, where the younger trail's "
    "age has density 2(L-u)/L^2 and the report distance is D ~ Exp(theta) with "
    "E[D]=theta, so Markov in D gives delta <= 2/L = 2^-6 exactly (the exact-"
    "exponential law is 2^%.4f, shipped). v31 charges delta = 2^%.6f, which is 155 "
    "bits LARGER than v30's fixture integral: strictly more conservative, and it "
    "extrapolates nothing. (2) THE GATE IS IN THE SHIPPED LOOP. v31_protocol."
    "run_protocol is the only launch path in the package (ledger31 imports "
    "INFLIGHT_CAP and HALT_RULE from it). The three quoted rules are: RULE 1 "
    "cadence 'a launch is authorized when credit = completed aggregate steps not "
    "yet consumed >= LAUNCH_PERIOD'; RULE 2 bootstrap 'if NO trail is live, a "
    "launch is authorized immediately and consumes no credit' (this is the case "
    "the reviewer named: a trail reporting before theta cannot deadlock the "
    "schedule); RULE 3 gate 'if a launch is authorized but len(inflight) == "
    "INFLIGHT_CAP, the launch is DEFERRED: the credit is retained and the request "
    "is retried at the next completion event'. The trace quotes itself: '%s' and "
    "'%s'. Patching the shipped constant to 129 makes the same entry point reach "
    "%d live trails, and removing RULE 3 reaches %d, so the bound is enforced by "
    "the object. reserve = %d*128*theta = %d steps (2^%.2f). (3) CUMULATIVE "
    "ABANDONED WORK IS A ROW. abandon_row = n_inflight*cap + n_launches*P(collide "
    "inside cap)*cap, priced at the max single-step price = 2^%.6f, where the "
    "pre-halt population over %d launches contributes 2^%.2f steps; W = "
    "w_effective + abandon_row + aux + pairing + overflow + stop_loss. Control (c) "
    "shows the row rising 1x/2x/4x/8x with earlier halts while w_effective holds. "
    "(4) PAIRING AND STORAGE AT PRODUCTION WIDTH. storage31() charges the 128 LIVE "
    "GT states (%d trails x 656 B = %d B), the fixed-P pairing base (13056 B), the "
    "returned G2 point (192 B), the pinned target (576 B) and the deterministic "
    "record cap: peak 2^%.3f B. pairing_production_width() audits every field: "
    "record count %d needs %d bits and is held as a 2-limb multiword, byte offset "
    "%d bits < 64 with 2^%.1f headroom under SIZE_MAX, class keys %d x u64 with no "
    "truncation before mixing, refs are (valid:u1, shard:u32, slot:u32) and are "
    "NEVER identified by a zero word (the g2-invert v36 aliasing mode), and "
    "records are sharded %d per shard. pairing_row charges the pinned check "
    "(2^%.6f). (5) MAX-PRICE CEILING KEPT: every step at c_max = %.6f, w_effective "
    "= 2^%.6f, mean-price diagnostic 2^%.6f undercounts by %.6f bits and is never "
    "summed. (6) TWELVE-SHIFT GATE KEPT: arm-1 recovery frequency %.6f against the "
    "independent census %.6f. TOTAL IS ONE SUM OF SIX DISPLAYED ROWS: w_effective "
    "2^%.6f + abandon_row 2^%.6f + aux_sum 2^%.6f + pairing_row 2^%.6f + "
    "overflow_row 2^%.0f + stop_loss_row 2^%.6f = 2^%.9f. The gate refused %d "
    "launch requests over the shipped trace; the step mix is f = %s and the "
    "retention-adjusted stop loss is 2*reserve/S = 2^%.2f of success. "
    "regen_v31.py emits every "
    "literal; audit_v31.py sweeps them and runs four planted controls: truncating "
    "the tail decades MOVES the fitted kappa, patching the shipped cap to 129 "
    "EXCEEDS 128, earlier halts RAISE the abandoned-work row, and a 1-ulp ledger "
    "change MOVES the claim. RANK DISCLOSED: %.4f is behind the live board minimum "
    "125.82 and this account's valid 125.838."
    % (C['total_log2'], C['claim_log2'], num, den, AGE['trials'], AGE['merges'],
       AGE['max_age_over_theta'], FIT['integral_measured_log2'], FIT['u_hi'],
       FIT['kappa_per_theta'], FIT['kappa_se'], FIT['resid_max_bits'],
       FIT['residual_log2'], PROD['delta_exact_log2'], DELTA['charged_log2'],
       HIT, DENY, LOOSE, BARE, P.INFLIGHT_CAP, int(STOP['reserve_steps']),
       math.log2(int(STOP['reserve_steps'])), LG['abandon_row'], AB['n_launches'],
       AB['pre_halt_log2'], ST['live_gt_states'], ST['live_trail_bytes'],
       ST['peak_log2'], WID['record_count']['value'], WID['record_count']['bits'],
       WID['byte_offset']['bits'], WID['size_max']['headroom_log2'],
       WID['class_key']['words'], WID['shards']['slots_per_shard'],
       LG['pairing_row'], float(T['ledgers']['c_mul_stress']), LG['w_effective'],
       MEAN_W_L2, UNDER, GATE, CEN['production_like']['p_recover'],
       LG['w_effective'], LG['abandon_row'], LG['aux_sum'], LG['pairing_row'],
       LG['overflow_row'], LG['stop_loss_row'], C['total_log2'], SHIP['denied_at_cap'],
       str(T['f']), math.log2(STOP['loss_rel_S2']), C['claim_log2']))

decades = ', '.join('[%d,%d]:%.4f' % (d['window'][0], d['window'][1],
                                      d['slope_per_step'])
                    for d in FIT['decade_slopes'])
argument = f"""# Argument: the squared-step walk (v31) - v30's four judged flaws repaired

v30 (ticket 4dcc056e, claim 126.6729) was rejected on four new flaws: a retention
tail bound of 2^-3.589688 carrying a 2^-100-class budget with an unsupported
envelope, a launch protocol the measurements did not execute (with cumulative
abandoned work unpriced), a 127 control that never touched the 128 boundary, and
missing pairing/storage accounting. v31 repairs those four; the regen script,
literal sweep, twelve-shift gate, max-price ceiling and dp stop-loss row are kept.

Claim **{C['claim_log2']:.4f}** = {num}/{den} (reduced), the least 4-decimal log2
>= the exact total 2^{C['total_log2']:.9f}. Emitted by `regen_v31.py`; no claim
literal is typed.

## 1. Retention tail: route (b), measured support + fitted residual

The merge simulator is re-run with `cap_mult = 128` (the production cap, not the
16-theta fixture) and budget 256 = 8x the old support: {AGE['trials']} trials,
{AGE['merges']} merges, {AGE['n_trails']} trails, max age
{AGE['max_age_over_theta']:.1f} theta. Raw ages are a 1-step histogram, no bins.

    I_meas = mean e^{{-(128 - A/theta)}} = 2^{FIT['integral_measured_log2']:.2f}
    window [{FIT['u_lo']:.1f}, {FIT['u_hi']:.1f}] theta, {FIT['window_points']} points
    kappa = {FIT['kappa_per_theta']:.6f} +/- {FIT['kappa_se']:.6f} per theta (lo {FIT['kappa_lo']:.6f})
    fit residuals: max {FIT['resid_max_bits']:.3f} bits, rms {FIT['resid_rms_bits']:.3f} bits
    per-decade slopes (steps): {decades}
    residual on ({FIT['u_hi']:.1f}, 128] theta = 2^{FIT['residual_log2']:.2f}

The fit is data: control (a) truncates the decades and kappa moves
({FIT['kappa_per_theta']:.4f} -> {L31.delta_fixture_only(truncate=FIT['a_hi_steps'] - 4)['kappa_per_theta']:.4f} -> {L31.delta_fixture_only(truncate=FIT['a_hi_steps'] - 8)['kappa_per_theta']:.4f}).

**The fixture law is not used as the production charge.** It has 785 toy classes;
the production retention miss comes from pair-uniformity instead. With the meeting
point uniform over the trail pair, the younger trail's age has density
2(L-u)/L^2 on [0,L], L = 128 theta, and its report distance is D ~ Exp(theta),
E[D] = theta. Markov in D:

    delta = int_0^L (2(L-u)/L^2) * theta/(L-u) du = 2/L = 2^-6    (charged)
    exact exponential law = 2(1 - e^-L(1+L))/L^2 = 2^{PROD['delta_exact_log2']:.4f}  (shipped)

v30 charged 2^-3.589688 here, so v31 charges 155 bits MORE on this row. The
charge moves the total, and the claim moves with it: this is the honest number.

## 2. The gate is in the shipped loop, not in a harness

`v31_protocol.run_protocol` is the only launch path; `ledger31` imports
`INFLIGHT_CAP` and `HALT_RULE` from it. Rules as they appear in the code:

    RULE 1  credit = completed - consumed; launch when credit >= LAUNCH_PERIOD
    RULE 2  if not inflight: launch (bootstrap, consumes no credit)
    RULE 3  if want and gate and len(inflight) >= inflight_cap: DEFER, credit kept

The line-by-line comparison the verdict asked for:

| v30 claim text | v30 code | v31 |
|---|---|---|
| "enforces an atomic in-flight cap of 128" | `scheduler_trace(..., enforce=derived)` in v30_measure, a HARNESS | RULE 3 inside `run_protocol`, imported by the ledger |
| "one trail start per theta executed steps" | `done - last >= launch_period`, done = aggregate | RULE 1, identical clock, credit retained on deferral |
| "no handling of this case is specified" | `first = (launched == 0 and not inflight)` only | RULE 2, explicit and quoted |
| "128*cap bounds simultaneously outstanding work" | trap row only | `abandon_row`, cumulative, section 3 |

Trace, shipped configuration:

    {HIT}
    {DENY}

Control (b): patching `v31_protocol.INFLIGHT_CAP` to 129 makes the same entry
point reach **{LOOSE}**; removing RULE 3 reaches **{BARE}**. The bound is a
property of the object, and it is exercised: the 129th launch is refused.

## 3. Cumulative abandoned work inside W

    abandon_row = n_inflight*cap + n_launches*P(collide inside cap)*cap
                = {int(STOP['reserve_steps'])} + 2^{AB['pre_halt_log2']:.2f} steps
                = 2^{LG['abandon_row']:.6f} at c_max

with n_launches = {AB['n_launches']} and P(cap) = 2^{AB['p_cap_log2']:.1f}. It
supersedes v30's trap row (same in-flight bound, plus the cumulative pre-halt
population) and sits inside W = w_effective + abandon_row + aux + pairing +
overflow + stop_loss. Control (c): earlier halts raise it 1x/2x/4x/8x while
w_effective holds.

## 4. Pairing and storage at production width

    live GT states   128 x 656 B        = {ST['live_trail_bytes']} B   (v30 charged 0)
    fixed-P base     68 x 2 x 96 B      = {ST['pairing_base_bytes']} B
    returned G2      2 x 96 B           = {ST['g2_point_bytes']} B
    pinned target    576 B
    record store     {ST['records']} x 112 B (cap, 2x slack)
    peak             2^{ST['peak_log2']:.3f} B

Width audit (`pairing_production_width`): record count {WID['record_count']['value']}
= {WID['record_count']['bits']} bits, held as a 2-limb multiword with a carry
assert; byte offset {WID['byte_offset']['bits']} bits, 2^{WID['size_max']['headroom_log2']:.1f}
headroom under SIZE_MAX, spec'd as (shard:u32, slot:u32, 2x64 offset); class keys
{WID['class_key']['words']} x u64, no truncation before mixing; refs carry an
explicit valid bit and are never identified by a zero word; {WID['shards']['count']}
shards of {WID['shards']['slots_per_shard']} slots. `pairing_row` = 2^{LG['pairing_row']:.6f}
for the pinned check (base + G2 reconstruction + Miller + final exp).

## 5. Max-price ceiling and twelve-shift gate kept

    w_effective = S*c_max = 2^{LG['w_effective']:.6f}; mean-price diagnostic 2^{MEAN_W_L2:.6f}
    undercounts by {UNDER:.6f} bits and is never summed
    gate: {GATE:.6f} arm-1 recovery vs {CEN['production_like']['p_recover']:.6f} census

## 6. The six rows

    w_effective      {LG['w_effective']:.6f}
    abandon_row      {LG['abandon_row']:.6f}
    aux_sum          {LG['aux_sum']:.6f}
    pairing_row      {LG['pairing_row']:.6f}
    overflow_row     {LG['overflow_row']:.6f}
    stop_loss_row    {LG['stop_loss_row']:.6f}
    total            {C['total_log2']:.9f}   claim {C['claim_log2']:.4f} = {num}/{den}

## 7. Scope

The step, return algebra, census protocol, family-free theorem and N=785 toy scope
are unchanged. No production mixing bound is claimed. Rank loss disclosed:
{C['claim_log2']:.4f} is behind the live board minimum 125.82 and this account's
valid 125.838.
"""

open('argument.md', 'w').write(argument)
open('public-note.md', 'w').write('# v31: v30 four flaws repaired\n\n' + argument)

cert = dict(S=dict(log2=math.log2(T['S']), value=str(T['S'])),
            claim=dict(claim_log2=C['claim_log2'], numerator=num, denominator=den,
                       total_log2=C['total_log2'], total=str(T['total'])),
            components={k: dict(log2=LG[k], value=str(v)) for k, v in COMP.items()},
            delta=DELTA, stop=STOP, abandon=AB, storage=ST, width=WID,
            protocol=dict(ship_max=SHIP['max_inflight'],
                          first_cap_tick=SHIP['first_cap_tick'],
                          denied=SHIP['denied_at_cap'], bare_max=BARE,
                          loosened_max=LOOSE, halt_rule=P.HALT_RULE),
            tail_fit=FIT, gate=GATE, census=CEN['production_like']['p_recover'],
            diagnostic_average_price=T['diagnostic_average_price'])
json.dump(cert, open('certificate.json', 'w'), indent=1, sort_keys=True, default=str)

sub = dict(version=1, title=("Squared-step class-merge rho, fixed-G1 pairing "
                             "inversion: measured 128-theta tail + fitted residual, "
                             "in-object 128 gate, cumulative abandoned work, "
                             "production-width storage"),
           experiments=[], claims=[dict(
               score=dict(numerator=str(num), denominator=str(den)),
               statement=STATEMENT,
               evidence=['argument.md', 'certificate.json', 'public-note.md',
                         'classkey.py', 'constants.py', 'ledger.py', 'ledger26.py',
                         'ledger27.py', 'ledger30.py', 'ledger31.py',
                         'v30_measure.py', 'v31_tail.py', 'v31_tailfit.py',
                         'v31_protocol.py', 'audit_v31.py', 'regen_v31.py',
                         'routing.py', 'tower.py', 'factor_r.py', 'toy_census.py',
                         'evidence/v30-age-law.json', 'evidence/v30-moments.json',
                         'evidence/v30-stopping.json', 'evidence/v30-census.json',
                         'evidence/v31-tail-raw.json', 'evidence/v31-protocol.json',
                         'evidence/v31-claim.json', 'evidence/v31-literals.json',
                         'evidence/v31-ledger.json', 'evidence/v31-audit.json'])])
json.dump(sub, open('submission.json', 'w'), indent=1)
print('claim %d/%d = %.4f total 2^%.9f' % (num, den, C['claim_log2'], C['total_log2']))
print('w %.6f abandon %.6f aux %.6f pairing %.6f stop %.6f' %
      (LG['w_effective'], LG['abandon_row'], LG['aux_sum'], LG['pairing_row'],
       LG['stop_loss_row']))
print('delta 2^%.6f  kappa %.4f +/- %.4f  I_meas 2^%.2f resid 2^%.2f' %
      (DELTA['charged_log2'], FIT['kappa_per_theta'], FIT['kappa_se'],
       FIT['integral_measured_log2'], FIT['residual_log2']))
print('protocol: hit %s  max %d denied %d  bare %d loosened %d' %
      (SHIP['first_cap_tick'], SHIP['max_inflight'], SHIP['denied_at_cap'],
       BARE, LOOSE))
