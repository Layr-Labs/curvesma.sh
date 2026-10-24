# Argument: the squared-step walk (v29) - v28's four judged flaws repaired

## 0. What this package is

v28 (ticket 347806a6, claim 126.2599, commit 331ed8d) was judged rejected on four
precision flaws: a code-vs-text exponent literal (code 2^-115, text 2^-171.66), an
age law conditioned on endpoint report matches, a one-trail reserve that does not
cover concurrent in-flight trails, and an unsupported zero-loss stopping assertion.
v29 repairs those four and re-derives the claim. The step, census protocol, return
algebra, price rows and family-free theorem are v28's, still shipped.

Claim **126.2599** = 1262599/10000 (reduced), the least 4-decimal log2 >=
the exact total 2^126.259807111. Emitted by `regen_v29.py` from
`ledger29`; no claim literal is typed by hand.

## 1. FIX 1 - one literal, computed once

v28's `ledger28.delta_charge` computed `8192*2**-128 = 2^-115`; its text claimed
`8192*exp(-128) = 2^-171.66`. It also dropped the L^2 normalisation and the cap
atom. v29 computes the closed form with both:

    L = min(Exp(theta), 128 theta), collision-weighted by L^2
    E[L^2]          = 2 theta^2 (1 - e^-128(1+128))
    E[e^(A/theta)]  = (128^3/3 + 128^2) / (2(1 - 129 e^-128))
    delta_closed    = e^-128 * E[e^(A/theta)] = 2^-166.216505

The charged value is max(raw-merge measured 2^-168.61,
closed 2^-166.216505) = 2^-166.216505. The
old 2^-171.66 is withdrawn as an unnormalised value; the old 2^-115 is withdrawn as
a typo-class literal. `audit_v29.py` maps every exponent in `submission.json` and
this file to its computing line or to a labelled historical reference: all Y.

## 2. FIX 2 - age at the merge event

v28 stored only distinguished-point report keys, so its age histogram was
conditioned on a detected endpoint match. v29 stores every visited class key and
counts every repeat as a merge event; after that repeat the follower's class-key
sequence is the elder's, so its next report is at a deterministic remaining
distance D. Exact 1-step histograms, no coarse bins, no endpoint substitution:

| arm (theta) | merges | report merges | mean age/theta | max age/theta | raw-merge delta | endpoint-only delta | censored merges |
|---|---|---|---|---|---|---|---|
| dp=1 (2) | 51656 | 16112 | 3.404 | 16.00 | 2^-168.61 | 2^-182.13 | 0/23166 |
| dp=2 (4) | 80243 | 9609 | 3.014 | 8.00 | 2^-177.20 | 2^-182.30 | 0/52719 |
| dp=3 (8) | 53725 | 3827 | 1.904 | 4.00 | 2^-181.11 | 2^-183.12 | 0/40178 |
| dp=4 (16) | 52346 | 2304 | 1.066 | 2.00 | 2^-182.91 | 2^-183.80 | 0/43137 |

The raw-merge residual integral is 2^-168.61 at arm 1
against 2^-182.13 for the endpoint-only law: the
endpoint conditioning understated the charge by 13.5 bits. The fixture cap is 16theta and
23166 of 51656 arm-1 merges are censored there, so a direct
production-cap miss count is not resolvable and is not claimed; the normalised closed
form, not a zero, is charged.

## 3. FIX 3 - concurrency-bound reserve, completed-step rule

The launch schedule is the one ledger27 already declares: one trail start per theta
executed steps (the census measures trail_starts ~ S/theta). A trail lasts at most
cap = 128 theta steps, so at most 128 trails are in flight. The rule is stated on
completed steps, not on realized work at an averaged price:

    reserve = 128 * 128 theta = 2^94.00 step-equivalents
    halt launches when completed steps >= S - reserve

Every launch occurs while completed steps < S - reserve, every in-flight trail
finishes within its 128theta ceiling, so total executed steps <= S deterministically.
The loss is charged, not asserted zero: reserve/S = 2^-31.82
and success ~ S^2, so the lost success share is 2^-30.82,
priced as stop_loss_row = 2^95.438583.

## 4. FIX 4 - unit-exact stopping counters with the verify gate

A repeated report key now counts only after `_solve_verify` recovers
x = (a2-a)/(b-b2) mod r and checks g^x = z. Recovery steps are unit-exact integers
per trial; the halt experiments are exact counts of trials with rec_step > S_toy - r
(no decile histogram arithmetic, no "0 of 210,000"):

| arm | lost trials at r=1,2,4,8,16 | lost fractions |
|---|---|---|
| 1 | 1:84, 16:1023, 2:165, 4:329, 8:600 | 1:0.0042, 16:0.0512, 2:0.0083, 4:0.0164, 8:0.0300 |
| 2 | 1:49, 16:614, 2:104, 4:193, 8:366 | 1:0.0024, 16:0.0307, 2:0.0052, 4:0.0097, 8:0.0183 |
| 3 | 1:22, 16:254, 2:41, 4:83, 8:142 | 1:0.0018, 16:0.0212, 2:0.0034, 4:0.0069, 8:0.0118 |
| 4 | 1:16, 16:188, 2:36, 4:63, 8:112 | 1:0.0013, 16:0.0157, 2:0.0030, 4:0.0053, 8:0.0093 |

The measured loss rises monotonically as the halt moves earlier (control d). The
production halt is 2^-30.82 of S, far below toy resolution, so the production row is
the arithmetic above, charged, not a measured zero.

## 5. The five rows

    w_stress        2^126.259807
    trap_charge     2^36.497842   (E[L^3]/(2E[L^2]) = 2.2650 theta)
    aux_sum         2^56.509866
    overflow_row    1
    stop_loss_row   2^95.438583
    ------------------------------------------
    total           2^126.259807111   claim 126.2599 = 1262599/10000

## 6. Unchanged and scope

The per-j return exponent and its 40,000 exact solves; the DP-predicate census;
the success model 1-exp(-S^2(1-f^2)/(2 mu N)) = 0.399237; f priced as 14566/41929;
the orbit theorem and class-level two-cycle; the N=785 toy scope and the
one-scope-up hypothesis, disclosed. No production mixing bound is claimed.

126.2599 is 0.4399 bits behind the live board minimum 125.82 and
0.4219 behind this account's valid 125.838. Disclosed.
