# Argument: the squared-step walk (v28) - v27c's three judged flaws repaired

## 0. What this package is

v27c (ticket fee953a4, claim 126.2672, intake CLEARED at 189,113 B) was judged
rejected on three precise flaws: an unconnected uniform-age premise, a false
size-bias assertion, and a stopping rule that did not exist while the
one-reference-operation overflow charge sat outside total27. v28 changes the three
things the judge named and nothing else. The step, the census protocol, the return
algebra, the price rows and the family-free theorem are v27c's, unchanged and still
shipped.

Claim **126.2599** = 1262599/10000 (reduced; least 4-decimal log2 >= the exact
total 2^126.259807111). It is 0.0073 bits BELOW v27c because the measured age law
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
| dp=1 (2) | 11580 | 1.272 | 6.00 | 4650 | 1426 | 436 | 103 | 2^-182.02 |
| dp=2 (4) | 7766 | 1.200 | 5.00 | 3513 | 1022 | 235 | 36 | 2^-182.27 |
| dp=3 (8) | 3170 | 0.935 | 3.62 | 1085 | 167 | 10 | 0 | 2^-183.01 |
| dp=4 (16) | 2109 | 0.551 | 1.81 | 208 | 0 | 0 | 0 | 2^-183.78 |

Two facts, both the judge's point: the age law is NOT uniform on [0,128 theta] (mean
age is 0.55-1.27 theta, i.e. young, because the restart biases it), and it is
arm-dependent - the deeper-DP arms are younger still. The retention functional is
recomputed on the EMPIRICAL law at the production cap, with no uniform assumption:

    delta_meas = (1/merges) * sum_i exp(-(128 - age_i/theta)),

which reads 2^-182.02 at the production-like arm (S/theta = 16) and is larger at
every deeper arm, so the production-like arm is the binding one. The uniform value
0.0078125 remains in `ledger28.retention_uniform` as the DISCLOSED-REJECTED method
and is not the charged number.

CHARGED delta = max(measured 2^-182.02, closed form 8192 e^-128 = 2^-171.66) =
2^-171.66. The closed form is the same age law already declared for trail lengths -
A size-biased by L^2 with L = min(Exp(theta),128theta), so
E[e^(A/theta)] = int_0^128theta a/theta^2 da = 128^2/2 = 8192 - and it is the
larger, so it is what the ledger carries. The one-scope-up hypothesis (a
785-element toy group standing in for the production group) is disclosed, exactly as
for the accepted census rows; the arithmetic value does not depend on it.

Consequence, charged and not banked: S moves from 2^125.8286 (delta = 1/128) to
2^125.8212, and w_stress moves 2^126.267148 -> 2^126.259807. The claim moves 0.0073
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
stop_loss_row = 2^88.4386. The halt-point distribution is shipped as well - the
recovery-work fraction histogram per arm, with 0 of 210,000 recovered trials inside
the final unit of the budget at any arm - and the production figure is taken from
the arithmetic because production stops at 1 - 2^-38.8 of its budget, which no toy
budget can realise. That limitation is stated, not glossed.

total28 is ONE SUM of FIVE DISPLAYED ROWS, and the one-reference-operation overflow
charge that v27c claimed separately and left OUT is now inside it:

    w_stress        2^126.259807
    trap_charge     2^36.6547
    aux_sum         2^56.509866
    overflow_row    2^0            (the Chernoff overflow tail, one reference op)
    stop_loss_row   2^88.4386
    ------------------------------------------
    total           2^126.259807111        claim 126.2599 = 1262599/10000

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
| 1 | 246033 | 505165 | 1942455 | 19619095 | 2.5250 |
| 2 | 132611 | 557983 | 4899783 | 76580515 | 1.9537 |
| 3 | 45393 | 353826 | 5091556 | 102880410 | 1.2629 |
| 4 | 31107 | 364659 | 6998077 | 168567861 | 0.7527 |

The production-like arm gives 2.5250 theta. The fixture cap is 16 theta (the
32-step toy budget binds first), so this is a TRUNCATED moment ratio taken from the
actual length sample - not the exponential 1.5 theta, and not assumed. It is larger
than both, so the trap row is charged at it:

    trap28 = ceil(S/theta) * E[L^2]/(2N) * (2.5250 theta) * c_mul = 2^36.6547

against v27c's 2^35.3257 (which charged 1.0 theta). The row is 2^89.6 below the
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
