#!/usr/bin/env python3
"""Emit certificate.json, argument.md, public-note.md from the computed ledger."""
import json
import math

import ledger28 as L28

c = L28.claim28()
T = c['T']
C = T['components']
lg = {k: math.log2(float(v)) for k, v in C.items()}
age = json.load(open('evidence/v28-age-law.json'))
mom = json.load(open('evidence/v28-moments.json'))['arms']
sto = json.load(open('evidence/v28-stopping.json'))['arms']
a1 = age['arms']['1']

cert = dict(
    claim=dict(numerator=c['numerator'], denominator=c['denominator'],
               claim_log2=c['claim_log2'], total_log2=c['total_log2'],
               total=str(T['total'])),
    S=dict(value=str(T['S']), log2=math.log2(T['S'])),
    delta=T['delta'], retention_uniform_rejected=age['uniform_rejected'],
    retention_analytic_log2=age['arithmetic_under_declared_law_log2'],
    delta_measured_log2={k: v['delta_measured_production_cap_log2']
                         for k, v in age['arms'].items()},
    age_law={k: dict(merges=v['merges'], trials=v['trials'], theta=v['theta'],
                     age_mean_over_theta=v['age_mean_over_theta'],
                     age_max_over_theta=v['age_max_over_theta'],
                     survival_gt_k_theta={kk: v['survival_gt_k_theta'][kk]
                                          for kk in ('1', '2', '3', '4')},
                     p_recover=v['p_recover']) for k, v in age['arms'].items()},
    moments={k: dict(n=v['n_trails'], sum_l=v['sum_l'], sum_l2=v['sum_l2'],
                     sum_l3=v['sum_l3'], ratio_over_theta=v['ratio_over_theta'])
             for k, v in mom.items()},
    wasted_theta_mult=T['trap']['wasted_theta_mult'],
    components={k: dict(value=str(v), log2=lg[k]) for k, v in C.items()},
    stop_rule=T['stop'],
    stopping_measured={k: dict(S_toy=v['S_toy'],
                               rec_mean_over_S=v['rec_step_mean_over_S'],
                               reserve_experiments=v['halt_experiments'])
                       for k, v in sto.items()},
    f_priced=str(T['f']),
    proof='audit_v28.py recomputes every literal and runs the four controls',
)
json.dump(cert, open('certificate.json', 'w'), indent=1, sort_keys=True, default=str)

arg = """# Argument: the squared-step walk (v28) - v27c's three judged flaws repaired

## 0. What this package is

v27c (ticket fee953a4, claim 126.2672, intake CLEARED at 189,113 B) was judged
rejected on three precise flaws: an unconnected uniform-age premise, a false
size-bias assertion, and a stopping rule that did not exist while the
one-reference-operation overflow charge sat outside total27. v28 changes the three
things the judge named and nothing else. The step, the census protocol, the return
algebra, the price rows and the family-free theorem are v27c's, unchanged and still
shipped.

Claim **%(claim).4f** = %(num)d/%(den)d (reduced; least 4-decimal log2 >= the exact
total 2^%(tot).9f). It is 0.0073 bits BELOW v27c because the measured age law
retires a premise that was costing 0.0075 bits; the other two fixes cost bits and
are charged.

## 1. The object

BLS12-381, fixed P in G1, uniform z in the order-r subgroup of GT. g = e(P,Q0),
solve z = g^x in GT, return Q = [x]Q0. Classical, single machine plus threads, no
amortisation, no cross-challenge precomputation. Step and routing are unchanged
from v27c (classkey -> one table row -> (y*TBL)^2 with labels 2(a+U), 2b).

## 2. FIX R1 - the retention integral is priced on a MEASURED age law

v27c charged `P(miss|merge) = (1-e^-128)/128 = 0.0078125` from an age uniform on
[0,128 theta]. The judge is right that first-DP restarts make the renewal age
nonuniform, so that integral is not connected to the success event. v28 measures
it. Per-trial traces, 64,000 trials, 24,625 merges, four report depths, age at the
colliding report counted from the FOLLOWER's own restart
(`v28_measure.py:_batch`, `evidence/v28-age-law.json`):

| arm (theta) | merges | mean age/theta | max age/theta | P(age>1t) | P(age>2t) | P(age>3t) | P(age>4t) | delta_meas at 128t |
|---|---|---|---|---|---|---|---|---|
%(agerows)s

Two facts, both the judge's point: the age law is NOT uniform on [0,128 theta] (mean
age is 0.55-1.27 theta, i.e. young, because the restart biases it), and it is
arm-dependent - the deeper-DP arms are younger still. The retention functional is
recomputed on the EMPIRICAL law at the production cap, with no uniform assumption:

    delta_meas = (1/merges) * sum_i exp(-(128 - age_i/theta)),

which reads 2^%(d1).2f at the production-like arm (S/theta = 16) and is larger at
every deeper arm, so the production-like arm is the binding one. The uniform value
0.0078125 remains in `ledger28.retention_uniform` as the DISCLOSED-REJECTED method
and is not the charged number.

CHARGED delta = max(measured 2^%(d1).2f, closed form 8192 e^-128 = 2^%(da).2f) =
2^%(da).2f. The closed form is the same age law already declared for trail lengths -
A size-biased by L^2 with L = min(Exp(theta),128theta), so
E[e^(A/theta)] = int_0^128theta a/theta^2 da = 128^2/2 = 8192 - and it is the
larger, so it is what the ledger carries. The one-scope-up hypothesis (a
785-element toy group standing in for the production group) is disclosed, exactly as
for the accepted census rows; the arithmetic value does not depend on it.

Consequence, charged and not banked: S moves from 2^125.8286 (delta = 1/128) to
2^125.8212, and w_stress moves 2^126.267148 -> 2^%(w).6f. The claim moves 0.0073
bits. No other row is resized.

## 3. FIX R3 - a stopping rule on REALIZED work, its measured loss, and the
## missing overflow charge

Declared rule (`ledger28.stop_loss`, `evidence/v28-stopping.json`):

    reserve = one maximal trail charge = 128*theta*c_mul = 2^87.00 step-equivalents
    stop launching new trails when cumulative charged work >= W - reserve

Proof on realized executions, not in expectation: every launch occurs while
cumulative < W - reserve, and each launched trail's charge is at most reserve
(128 theta steps plus its aux share), so when the last in-flight trail completes

    total <= (W - reserve) + reserve = W.

This is arithmetic on the launch loop; no Chernoff tail, no expected charge, and no
cap language is used for the step sum. The rule's SUCCESS LOSS is measured and
charged: the rule stops 128 theta steps short of S, so the unspent tail is
reserve/S = 2^-38.80 of S and, since success ~ S^2, the success share it would have
bought is 2*reserve/S = 2^-37.82; priced as that fraction of W_stress it is
stop_loss_row = 2^%(sl).4f. The halt-point distribution is shipped as well - the
recovery-work fraction histogram per arm, with 0 of 210,000 recovered trials inside
the final unit of the budget at any arm - and the production figure is taken from
the arithmetic because production stops at 1 - 2^-38.8 of its budget, which no toy
budget can realise. That limitation is stated, not glossed.

total28 is ONE SUM of FIVE DISPLAYED ROWS, and the one-reference-operation overflow
charge that v27c claimed separately and left OUT is now inside it:

    w_stress        2^%(w).6f
    trap_charge     2^%(tr).4f
    aux_sum         2^%(aux).6f
    overflow_row    2^0            (the Chernoff overflow tail, one reference op)
    stop_loss_row   2^%(sl).4f
    ------------------------------------------
    total           2^%(tot).9f        claim %(claim).4f = %(num)d/%(den)d

`audit_v28.py` asserts that the total equals the sum of exactly these five rows and
that deleting overflow_row flips that assertion (control d).

## 4. FIX R2 - the wasted-work term is E[L^3]/(2E[L^2]), computed from the sums

v27c asserted that the size-biased prefix mean is E[L^2]/(2E[L]) <= theta. That is
false as written and is WITHDRAWN. Conditioning a collision weights lengths by
L^2, and a uniform prefix of such a trail has mean L/2, so the conditional prefix
mean is

    E[L^3] / (2 E[L^2])           (= 3 theta/2 for untruncated exponential)

The exact integer sums are shipped per arm (`evidence/v28-moments.json`), so the
ratio is computed from data and moves if a sample moves (control b):

| arm | trail starts | sum L | sum L^2 | sum L^3 | E[L^3]/(2E[L^2])/theta |
|---|---|---|---|---|---|
%(momrows)s

The production-like arm gives %(r3).4f theta. The fixture cap is 16 theta (the
32-step toy budget binds first), so this is a TRUNCATED moment ratio taken from the
actual length sample - not the exponential 1.5 theta, and not assumed. It is larger
than both, so the trap row is charged at it:

    trap28 = ceil(S/theta) * E[L^2]/(2N) * (%(r3).4f theta) * c_mul = 2^%(tr).4f

against v27c's 2^35.3257 (which charged 1.0 theta). The row is 2^%(gap).1f below the
step term, so the 4-decimal total does not move from it - but the claim text now
states the corrected quantity, which is what the judge asked for.

## 5. Unchanged, and why

The per-j return exponent `y^(4-p^j) = (m1^4 m2^2)^(-1)` with 40,000 exact solves
and 0 misses; the DP-predicate census (35,211 collisions / 35,207 recoveries =
0.9998863991 at dp=1, 23,366/23,365 = 0.9999572028 at dp=2, 7,449/7,449 and
4,503/4,503); the disclosed model 1-exp(-S^2(1-f^2)/(2 mu N)) = 0.399237 at the
judge's own S=32, N=785, f=3/8, mu=1.1; f = 1734128/4991779 = 0.3473967898 priced
as 14566/41929; the orbit theorem (384/384 inverse-pair columns, planted mutation
flags exactly the mutated column); the class-level two-cycle and the N=785 toy
scope. None of these is restated as a new claim.

## 6. Scope and rank

The algebra proves uniqueness for a FIXED multiplier pair and shift; the census is a
785-element toy at four report depths; the age law is measured there and priced
against the declared production law with the one-scope-up hypothesis disclosed. No
production mixing bound and no minimality theorem is claimed.

126.2599 is 0.4399 bits behind the live board minimum 125.82 and 0.4219 behind this
account's valid 125.838. Disclosed.
""" % dict(claim=c['claim_log2'], num=c['numerator'], den=c['denominator'],
           tot=c['total_log2'], w=lg['w_stress'], tr=lg['trap_charge'],
           aux=lg['aux_sum'], sl=lg['stop_loss_row'],
           d1=a1['delta_measured_production_cap_log2'],
           da=age['arithmetic_under_declared_law_log2'],
           r3=T['trap']['wasted_theta_mult'],
           gap=lg['w_stress'] - lg['trap_charge'],
           agerows='\n'.join(
               '| dp=%s (%d) | %d | %.3f | %.2f | %s | %s | %s | %s | 2^%.2f |' % (
                   k, v['theta'], v['merges'], v['age_mean_over_theta'],
                   v['age_max_over_theta'], v['survival_gt_k_theta']['1'],
                   v['survival_gt_k_theta']['2'], v['survival_gt_k_theta']['3'],
                   v['survival_gt_k_theta']['4'],
                   v['delta_measured_production_cap_log2'])
               for k, v in sorted(age['arms'].items())),
           momrows='\n'.join(
               '| %s | %d | %d | %d | %d | %.4f |' % (
                   k, v['n_trails'], v['sum_l'], v['sum_l2'], v['sum_l3'],
                   v['ratio_over_theta']) for k, v in sorted(mom.items())))
open('argument.md', 'w').write(arg)

note = """# v28: the v27c judge's three flaws, measured

v27c (fee953a4, 126.2672) was rejected on three math/accounting flaws. This package
repairs exactly those and re-derives the claim.

- **Age law, not a uniform premise.** 64,000 toy trials / 24,625 merges at four
  report depths measure the age at merge conditional on the merge: mean 1.27, 1.20,
  0.94, 0.55 theta, max 6.0 theta at the production-like arm, survival
  P(age>k theta) falling to 103/11,580 by k=4. The 1/128 uniform integral is kept
  only as DISCLOSED-REJECTED. Charged delta = max(measured 2^-182.02, closed form
  8192 e^-128 = 2^-171.66) = 2^-171.66; S falls 2^125.8286 -> 2^125.8212.
- **E[L^3]/(2E[L^2])**, from exact integer sums of L, L^2, L^3 (shipped per arm):
  5.0503 steps = 2.5250 theta at the production-like arm, a truncated moment ratio
  from the sample, charged, not the exponential 1.5 theta and not v27c's false
  "<= theta". Trap row 2^35.3257 -> 2^36.6547.
- **Stopping rule on realized work**: stop launching at W - reserve with
  reserve = 128 theta c_mul = 2^87.00; each launched trail charges <= reserve, so
  total <= W deterministically. Loss measured/charged: 2*reserve/S = 2^-37.82 of
  success = 2^88.4386 of work. total28 now sums five displayed rows INCLUDING the
  one-reference-operation overflow charge v27c left out.
- Claim **126.2599 = 1262599/10000** (reduced); total 2^126.259807.

Controls in `audit_v28.py`: the uniform integral is rejected not charged; perturbing
one length sample moves the moment ratio; the launch sim overruns W with the rule
disabled and never with it; deleting the overflow row flips the sum assertion.
"""
open('public-note.md', 'w').write(note)
print('certificate.json, argument.md, public-note.md written; claim %.4f' % c['claim_log2'])
