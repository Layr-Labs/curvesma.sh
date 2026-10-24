# v30: v29's four flaws repaired in the arithmetic

# Argument: the squared-step walk (v30) - v29's four judged flaws repaired

v29 (ticket d7952ec2, claim 126.2599) cleared intake and was rejected on four NEW
precision flaws: a zero-shift-only recovery gate, a concurrency bound that mixed
the aggregate launch clock with the own-step trail cap, a mean-price work ceiling,
and an unsupported retention tail. v30 repairs those four; the regen script,
literal sweep, un-binned raw-age integral and measured stopping rows are kept.

Claim **126.6729** = 1266729/10000 (reduced), the least 4-decimal log2
>= the exact total 2^126.672831633. Emitted by `regen_v30.py` from
`ledger30`; no claim literal is typed.

## 1. Twelve-shift recovery gate

v29's `_solve_verify` solved only `x=(a2-a)/(b-b2)`. A Frobenius-class collision
is `a + b*x = p^j*(a2 + b2*x)` for one j in Z/12, so v29 rejected 11 of 12
legitimate witnesses and its arm-1 recovery frequency was 0.066. v30 tries every
shift and verifies `g^x = z`:

    x_j = (p^j*a2 - a) * (b - p^j*b2)^-1 mod r,  j = 0..11

Measured arm-1 recovery frequency: **0.579100** (11582/20000 trials);
**0.999976** of 41426 merges verify; shift histogram
{'0': 5532, '1': 2209, '10': 3154, '11': 2198, '2': 3526, '3': 3817, '4': 3913, '5': 2051, '6': 5055, '7': 2255, '8': 3610, '9': 4105}. The independent census production-like rate is
**0.586783**; the two literals agree to 0.008 and are both shipped.

## 2. Age at the merge event, support stated

Every visited class key is stored and every repeat is a merge event; ages are
exact integer steps, 1-step histograms, no bins. The support is [1,16] theta (the
fixture cap), not the production [1,128] theta:

| arm | merges | report merges | mean age/theta | max age/theta | raw-merge delta | endpoint-only delta |
|---|---|---|---|---|---|---|
| dp=1 (2) | 41426 | 11583 | 3.939 | 16.00 | 2^-168.25 | 2^-182.02 |
| dp=2 (4) | 70297 | 7768 | 3.155 | 8.00 | 2^-177.09 | 2^-182.28 |
| dp=3 (8) | 49733 | 3170 | 1.972 | 4.00 | 2^-181.05 | 2^-183.01 |
| dp=4 (16) | 49520 | 2109 | 1.095 | 2.00 | 2^-182.87 | 2^-183.78 |

The raw-merge support integral is 2^-168.25 at arm 1
against 2^-182.02 for the endpoint-only law.

## 3. Retention tail: route (a), charged at a conservative constant

The support stops at 16 theta; the production cap is 128 theta. v30 does NOT
present the support integral as a full-range one. Beyond support it charges
`e^-128 * integral_16^128 e^a P(A>a) da` with
`P(A>a) <= s4*exp(-lam*(a-4))`, where s4 is the measured survival at 4 theta and
lam is HALF the measured 1->4 theta decay (slower decay = larger tail):

    s1=0.676580, s4=0.336262,
    lam=0.116527/theta,
    tail row = 2^-3.589688.

Charged delta = max(support 2^-168.25 + tail
2^-3.589688, closed 2^-166.216505) =
2^-3.589688. Route (a) is disclosed; the simulation was
not extended.

## 4. Concurrency: derivation withdrawn, atomic cap enforced

The code lines are: launch cadence `one trail start per theta executed steps`
(aggregate completed steps), trail cap `cap = 128*theta` of OWN steps. Those
clocks do not give `ceil(cap/theta)=128`. The bare-schedule max-launch-rate trace
measures **1023** live trails at a 4-cap horizon. v30 therefore WITHDRAWS the
derivation and ENFORCES an atomic in-flight cap of 128; the same
trace then measures **128**. A cap one smaller measures 127 and is shown to
under-cover by one trail's charge, so the audit's assert trips on it.

    reserve = 128 * 128 theta = 2^94.00 step-equivalents
    halt launches when completed steps >= S - reserve

Every in-flight trail is within its 128 theta cap at the halt, so extra work
<= n_inflight*cap = reserve and total completed steps <= S.

## 5. Max-price work ceiling, not a mean

A mean is not an upper bound. v30 charges every completed step at the dearer step
price `c_max = 1.704385`:

    w_stress = S * c_max = 2^126.672832
    mean-price diagnostic (NOT summed) = 2^126.342163; undercount 0.330669 bits

The trap row is the deterministic max single-trail charge times the enforced
in-flight cap: `128*128 theta` step-equivalents at c_max =
2^94.769251. The measured E[L^3]/(2E[L^2]) =
2.5252 theta is a diagnostic.

## 6. Stopping loss, twelve-shift, unit-exact

| arm | lost trials at r=1,2,4,8,16 | lost fractions |
|---|---|---|
| 1 | 1:487, 16:7968, 2:947, 4:1900, 8:3908 | 1:0.0244, 16:0.3984, 2:0.0474, 4:0.0950, 8:0.1954 |
| 2 | 1:370, 16:5767, 2:769, 4:1541, 8:3066 | 1:0.0185, 16:0.2883, 2:0.0384, 4:0.0770, 8:0.1533 |
| 3 | 1:201, 16:2495, 2:398, 4:786, 8:1426 | 1:0.0168, 16:0.2079, 2:0.0332, 4:0.0655, 8:0.1188 |
| 4 | 1:127, 16:1710, 2:275, 4:528, 8:984 | 1:0.0106, 16:0.1425, 2:0.0229, 4:0.0440, 8:0.0820 |

The production loss is charged: reserve/S = 2^-31.90,
success ~ S^2, so 2*reserve/S = 2^-30.90 of success,
stop_loss_row = 2^95.769251.

## 7. The five rows

    w_stress        126.672832
    trap_charge     94.769251
    aux_sum         56.592221
    overflow_row    0.000000
    stop_loss_row   95.769251
    total           126.672831633   claim 126.6729 = 1266729/10000

## 8. Scope

The step, return algebra, census protocol, family-free theorem and N=785 toy
scope are unchanged. No production mixing bound is claimed. Rank loss disclosed:
126.6729 is behind the live board minimum 125.82 and this account's
valid 125.838.
