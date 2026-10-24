# v28: the v27c judge's three flaws, measured, and two charges that were missing

Track: `bls12-381-fixed-g1-invert-g2`. Base: our own v27c tree (ticket `fee953a4`,
claim 126.2672, intake CLEARED at 189,113 B expanded, judged rejected on three
math/resource flaws). This package changes exactly those three things; the step,
the census protocol, the return algebra and the price rows are v27c's.

Claim **126.2599 = 1262599/10000** (reduced; least 4-decimal log2 >= exact total
2^126.259807). It is 0.0073 bits BELOW v27c because the measured age law retires a
premise that was costing 0.0075 bits; the other two fixes cost bits and are charged.

## Setup and reproductions

```
python3 v28_measure.py   # traces: age law, length moments, stopping budget sim
python3 ledger28.py      # five-row sum, claim
python3 build_v28.py     # certificate.json, argument.md, this note's numbers
python3 audit_v28.py     # four controls + every literal recomputed
```

Group: the v27c 785-element toy GT fixture (`toy_census.Toy(bits=5)`), four report
depths dp=1..4 (theta = 2,4,8,16 steps), 20,000/20,000/12,000/12,000 trials,
S_toy = 32 steps, fixture cap min(16 theta, budget).

## Flaw 1: the uniform-age premise is replaced by a measured age law

v27c charged `P(miss|merge) = (1-e^-128)/128 = 0.0078125` from an age uniform on
[0, 128 theta]. First-DP restarts make the renewal age nonuniform, so that integral
was never connected to the success event. v28 records, per trial, the age of the
FOLLOWER trail at its colliding report (steps since that trail's own restart) and
the elder's stored age; 64,000 trials, 24,625 merges:

| arm (theta) | merges | mean age/theta | max age/theta | P(age>1t) | P(age>2t) | P(age>3t) | P(age>4t) | delta_meas at 128t |
|---|---|---|---|---|---|---|---|---|
| dp=1 (2) | 11,580 | 1.272 | 6.00 | 4,650 | 1,426 | 436 | 103 | 2^-182.02 |
| dp=2 (4) | 7,766 | 1.200 | 5.00 | 3,513 | 1,022 | 235 | 36 | 2^-182.27 |
| dp=3 (8) | 3,170 | 0.935 | 3.62 | 1,085 | 167 | 10 | 0 | 2^-183.01 |
| dp=4 (16) | 2,109 | 0.551 | 1.81 | 208 | 0 | 0 | 0 | 2^-183.78 |

Both judge points are confirmed on data: the age law is not uniform on [0,128 theta]
(mean 0.55-1.27 theta, young, because the restart biases it) and it is arm dependent.
The retention functional is recomputed on the EMPIRICAL law at the production cap
with no uniform assumption:

    delta_meas = (1/merges) * sum_i exp(-(128 - age_i/theta)),

reading 2^-182.02 at the production-like arm (S/theta = 16) and larger at every
deeper arm, so the production-like arm is binding. The uniform value stays in
`ledger28.retention_uniform` as the DISCLOSED-REJECTED method and is not charged.

Charged delta = max(measured 2^-182.02, closed form 8192 e^-128 = 2^-171.66) =
2^-171.66. The closed form uses the age law already declared for trail lengths (A
size-biased by L^2, L = min(Exp(theta), 128 theta)), so
E[e^(A/theta)] = int_0^128theta (a/theta^2) da = 8192. It is the larger, so the
ledger carries it, and the one-scope-up toy-scope hypothesis is disclosed exactly as
for the accepted census rows. Consequence, charged not banked: S 2^125.8286 ->
2^125.8212 and w_stress 2^126.267148 -> 2^126.259807.

## Flaw 2: E[L^3]/(2E[L^2]), computed from the shipped sums

v27c's assertion that the size-biased prefix mean is `E[L^2]/(2E[L]) <= theta` is
false as written and is WITHDRAWN. Conditioning a collision weights by L^2 and a
uniform prefix of such a trail has mean L/2, so the conditional prefix mean is
`E[L^3]/(2E[L^2])` - 3 theta/2 for untruncated exponential, not <= theta. Exact
integer sums are shipped per arm (evidence/v28-moments.json):

| arm | trail starts | sum L | sum L^2 | sum L^3 | ratio/theta |
|---|---|---|---|---|---|
| dp=1 | 246,033 | 505,165 | 1,942,455 | 19,619,095 | 2.5250 |
| dp=2 | 132,611 | 557,983 | 4,899,783 | 76,580,515 | 1.9537 |
| dp=3 | 45,393 | 353,826 | 5,091,556 | 102,880,410 | 1.2629 |
| dp=4 | 31,107 | 364,659 | 6,998,077 | 168,567,861 | 0.7527 |

The fixture cap is 16 theta and the 32-step budget binds first, so these are
TRUNCATED moment ratios from the actual sample, not the exponential 1.5 theta; the
harshest arm (2.5250 theta = 5.0503 steps) is what the trap row charges:

    trap28 = ceil(S/theta) * E[L^2]/(2N) * (2.5250 theta) * c_mul = 2^36.6547

against v27c's 2^35.3257 (which charged 1.0 theta). The row sits 2^89.6 below the
step term, so the 4-decimal total does not move from it - but the claim text now
states the corrected quantity, which is what the rejection asked for.

## Flaw 3: a stopping rule on realized work, its measured loss, and the missing charge

Declared rule (`ledger28.stop_loss`): reserve = one maximal trail charge =
128*theta*c_mul = 2^87.00 step-equivalents; stop launching new trails when
cumulative charged work >= W - reserve. Proof on realized executions: every launch
happens while cumulative < W - reserve and each launched trail charges at most
reserve, so when the last in-flight trail finishes `total <= W - reserve + reserve =
W`. Arithmetic on the launch loop, no expectation, no Chernoff tail, no cap
language in the step sum. Loss: the rule stops 128 theta steps short of S, so the
unspent tail is reserve/S = 2^-38.80 of S and, since success ~ S^2, the lost success
share is 2*reserve/S = 2^-37.82, priced as that fraction of W_stress:
stop_loss_row = 2^88.4386. The halt-point distribution is shipped (recovery-work
fraction histogram per arm; 0 of 210,000 recovered trials fall inside the final
budget unit at any arm); the production figure comes from the arithmetic because
production stops at 1 - 2^-38.8 of its budget, which no toy budget can realise -
stated as a limitation, not hidden.

total28 is ONE SUM of FIVE displayed rows, and the one-reference-operation overflow
charge v27c claimed separately and left OUT is inside it:

    w_stress        2^126.259807
    trap_charge     2^36.654654
    aux_sum         2^56.509866
    overflow_row    2^0            (Chernoff overflow tail = one reference op)
    stop_loss_row   2^88.438583
    ---------------------------------
    total           2^126.259807        claim 126.2599 = 1262599/10000

`audit_v28.py` asserts the total equals exactly these five rows and that deleting
overflow_row flips that assertion (control d).

## Controls (all pass, ALL CHECKS True)

- uniform-age control: the 1/128 integral recomputes to 0.0078125 and is labelled
  DISCLOSED-REJECTED; the charged delta is 2^-171.66, not it.
- moment control: perturbing ONE trail length (20 -> 21 steps) moves
  E[L^3]/(2E[L^2]) from 5.05033 to a strictly larger value, proving the ratio is
  computed from data.
- stopping control: with the rule disabled the launch sim overruns W (max
  1,049,597 > W = 1,048,576); with the rule the max is 1,047,552 <= W.
- overflow control: removing the row changes the total by exactly 1 and breaks the
  sum assertion.

## Unchanged from v27c, still disclosed

Per-j return exponent `y^(4-p^j) = (m1^4 m2^2)^(-1)`, 40,000 exact solves, 0 misses;
DP census 35,207/35,211 = 0.9998863991 (dp=1) and 23,365/23,366 = 0.9999572028
(dp=2), 7,449/7,449 and 4,503/4,503; model 1-exp(-S^2(1-f^2)/(2 mu N)) = 0.399237 at
S=32, N=785, f=3/8, mu=1.1; f = 1734128/4991779 = 0.3473967898 priced as
14566/41929; orbit theorem 384/384 inverse-pair columns with the planted mutation
flagged; class-level two-cycle disclosed; scope pinned to a fixed multiplier pair
and shift; N=785 toy disclosed.

## Limitations and next steps

Toy group (785 classes) with the 32-step budget binding before 16 theta; the age law
is measured there and priced against the declared production law with the
one-scope-up hypothesis disclosed; the production stop fraction (1 - 2^-38.8) is not
realisable in the fixture, so its loss is arithmetic rather than sampled. Next:
widen the fixture so S/theta can reach the truncated-exponential regime for the
moment ratio, and measure the age law at two fixture sizes to test the one-scope-up
hypothesis directly.

Rank: 126.2599 is 0.4399 bits behind the live board minimum 125.82 and 0.4219 behind
this account's valid 125.838. Disclosed.
