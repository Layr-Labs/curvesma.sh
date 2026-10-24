# v29: v28's four flaws repaired in the arithmetic, not in prose

v28 (ticket 347806a6, claim 126.2599, commit 331ed8d) cleared intake and was judged
rejected on four precision defects. This note is organised against those four
defects; each repair is executable in the package and each number below is emitted
by `regen_v29.py` from `ledger29`, then swept by `audit_v29.py`.

## 1. The literal that did not match its code

v28's `delta_charge` computed `8192 * 2**-128 = 2^-115`, while its text claimed
`8192 * exp(-128) = 2^-171.66`. Both were wrong: the first is a typo-class exponent,
the second drops the L^2 normalisation and the cap atom. v29 computes one value:

    L = min(Exp(theta), 128 theta), collision-weighted by L^2
    E[L^2]         = 2 theta^2 (1 - e^-128 (1 + 128))
    E[e^(A/theta)] = (128^3/3 + 128^2) / (2 (1 - 129 e^-128))
    delta_closed   = e^-128 E[e^(A/theta)] = 2^-166.216505

Charged delta = max(raw-merge measured 2^-168.61,
closed 2^-166.216505) =
2^-166.216505. The literal sweep in
`audit_v29.py` extracts every `2^...` token from the submission statement and
`argument.md` and maps it to a computing line or a labelled historical reference:
20/20 tokens Y.

## 2. Age at the merge event, not at a report

v28 stored only report keys, so its ages were conditioned on reaching an
already-stored report. v29 stores every visited class key; every repeat is a merge
event; each merge records the follower age A in exact integer steps (1-step
histogram, no coarse bins) and the remaining distance D to the follower's next
report. Raw-merge versus endpoint-only integrals:

| arm (theta) | merges | report merges | mean age/theta | raw-merge delta | endpoint-only delta | censored |
|---|---|---|---|---|---|---|
| dp=1 (2) | 51656 | 16112 | 3.404 | 2^-168.61 | 2^-182.13 | 23166 |
| dp=2 (4) | 80243 | 9609 | 3.014 | 2^-177.20 | 2^-182.30 | 52719 |
| dp=3 (8) | 53725 | 3827 | 1.904 | 2^-181.11 | 2^-183.12 | 40178 |
| dp=4 (16) | 52346 | 2304 | 1.066 | 2^-182.91 | 2^-183.80 | 43137 |

At arm 1 the raw-merge integral is 2^-168.61
against 2^-182.13 endpoint-only: the
endpoint conditioning understated the charge by 13.5 bits. The fixture cap is
16theta and 23166 of 51656
arm-1 merges are censored there, so a direct production-cap miss count is not
resolvable and is not claimed; the normalised closed form, not a zero, is charged.
The 44.8% figure is disclosed as a fixture limitation, not as production evidence.

## 3. Concurrency-bound reserve and a completed-step rule

The launch schedule is the one ledger27 already declares and the census measures:
one trail start per theta executed steps (trail starts ~ S/theta). A trail lasts at
most cap = 128 theta steps, so at most ceil(cap/theta) = 128 trails are in flight.
The rule is stated on COMPLETED STEPS, not on realized work at an averaged price:

    reserve = 128 * 128 theta = 2^94.00 step-equivalents
    halt launches when completed steps >= S - reserve

Every launch occurs while completed steps < S - reserve; every in-flight trail
finishes within its ceiling; total executed steps <= S deterministically. The
success loss is charged, not asserted zero: reserve/S = 2^-31.82,
success ~ S^2, so 2*reserve/S = 2^-30.82 of
success is charged as stop_loss_row = 2^95.438583.

## 4. Unit-exact stopping counters with the verify gate

A repeated report key counts only after the solver recovers
x = (a2-a)/(b-b2) mod r and checks g^x = z. Recovery steps are unit-exact integers
per trial, so each halt experiment is an exact count of trials with
rec_step > S_toy - r, not decile-histogram arithmetic. Measured lost fractions:

| arm | lost fractions at r = 1,2,4,8,16 steps short of S_toy |
|---|---|
| 1 | 1:0.0042, 2:0.0083, 4:0.0164, 8:0.0300, 16:0.0512 |
| 2 | 1:0.0024, 2:0.0052, 4:0.0097, 8:0.0183, 16:0.0307 |
| 3 | 1:0.0018, 2:0.0034, 4:0.0069, 8:0.0118, 16:0.0212 |
| 4 | 1:0.0013, 2:0.0030, 4:0.0053, 8:0.0093, 16:0.0157 |

The measured loss rises monotonically as the halt moves earlier (control d: halving
the budget multiplies the lost fraction by roughly 10x-30x). The production halt is
2^-30.82 of S, far below toy resolution, so the production row is the arithmetic
above, charged, not a measured zero. The v28 sentence "0 of 210,000 recovered trials
in the final unit" is withdrawn: it was decile-histogram arithmetic.

## 5. The five displayed rows

    w_stress        2^126.259807
    trap_charge     2^36.497842   (E[L^3]/(2E[L^2]) = 2.2650 theta)
    aux_sum         2^56.509866
    overflow_row    1
    stop_loss_row   2^95.438583
    ------------------------------------------
    total           2^126.259807111
    claim           126.2599 = 1262599/10000

The moment sums moved from v28's 2.5250 theta to 2.2650 theta because the verify
gate lengthens a trial after a rejected merge; the sums are exact integers in
`evidence/v29-moments.json` and the harsher arm is charged. The trap row is 2^89.8
below the step term, so the 4-decimal claim does not move from it, but the text now
states the quantity the ledger computes.

## 6. Controls

`audit_v29.py` runs four planted controls (4/4 fire): (a) one-ulp claim-grid
perturbation moves the reduced rational by 1/10000 and one priced-fraction ulp
changes the exact total fraction; (b) coarsening the exact age histogram to
[1,2),[2,4),[4,8),[8,16),[16,32) makes the endpoint integrals visibly differ from
the raw integral (2^-175.92 lower, 2^-164.42 upper, raw 2^-168.61 inside) -- the
raw integral is used, and no binning is used for the charge; (c) the launch sim at
maximum rate attains the 128-trail bound, which is exactly the reserve's
N_inflight; (d) forcing the halt earlier raises the measured loss at every arm.
Twenty independent five-row recomputations reproduce the sum exactly (20/20), and
every evidence path is a bare in-tree name with no `..` component.

## 7. Scope and rank

The algebra proves uniqueness for a fixed multiplier pair and shift; the census is
a 785-element toy at four report depths; the age law is measured there and priced
against the declared production law with the one-scope-up hypothesis disclosed. No
production mixing bound and no minimality theorem is claimed. 126.2599 is 0.4399
bits behind the live board minimum 125.82 and 0.4219 behind this account's valid
125.838. Disclosed.
