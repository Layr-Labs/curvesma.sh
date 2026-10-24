# Argument: the squared-step walk (v27) — the judge's three defects fixed

## 0. What this package is

A validity play, corrected. v26 (claim 126.2508, ticket `1e3fd098`) was rejected on
three MATH defects, not on packaging: the census never executed the distinguished
point predicate, the model dropped the `(1-f^2)/mu` factor, the trap row charged
four trails against an `S/theta` restart schedule, and the return exponent was
cubing for every Frobenius class. This package fixes all four, re-measures, and
re-derives the claim from the corrected ledger.

Claim **126.2672** = 78917/625. It is 0.4472 bits behind the live board minimum
125.82 and 0.4292 behind this account's valid 125.838 row — disclosed, not hidden.
The claim is 0.0164 bits HIGHER than v26's because two corrections cost bits: the
step budget now pays for DP detection retention (+0.0075) and the executed
squaring fraction is priced at the new measured 0.346631 instead of 0.355811
(+0.0089).

## 1. The object

BLS12-381, fixed P in G1, uniform z in the order-r subgroup of GT. `g = e(P,Q0)`,
solve `z = g^x` in GT, return `Q = [x]Q0`. Classical, single machine plus threads,
no amortisation, no cross-challenge precomputation.

## 2. The step

    kappa = classkey(y)                       # the 12-element sigma-orbit class
    t     = routing.fields(kappa).row         # one table row, from the class only
    y'    = (y * TBL[t][kappa])^2             # multiply, then cyclotomic square
    a'    = 2*(a + U[t][kappa]) mod r ; b' = 2*b mod r

No prev_row, no ring, no history: a class repeat is a state repeat, so the cycle
abort is sound. Equivariant: `(sigma^j y · sigma^j m)^2 = sigma^j((y m)^2)`.
Family-free by the orbit theorem (`evidence/impossibility-theorem.md`, 384/384
inverse-pair columns, 768/768 orbit columns, planted mutation flags exactly the
mutated column).

## 3. FIX 3 — the return exponent is per class index, and the recovery is dynamic

The v26 text wrote the two-step return as `y^3 = sigma^j/(m1^4 m2^2)`. That is
wrong for every `j != 0`: `sigma^j(y) = y^(p^j)` in the order-r subgroup, so the
return equation is

    y^(4 - p^j) = (m1^4 m2^2)^(-1)      ⇒      y = [ (m1^4 m2^2)^(-1) ]^(e_j^(-1) mod r)

with `e_j = (4 - p^j) mod r`. Cubing covers only `j = 0` (where `e_0 = 3`).
`recover_y_old()` in `v27_measure.py` is the v26 formula, kept as the control.

Inversion is on the SUBGROUP order r (prime), not on p-1: `gcd(e_j, r) = 1` unless
`e_j ≡ 0 (mod r)`, and it never is for j = 0..11 on the shipped fixture
(`recovery_proof.gcd_fail = 0` over 40,000 solves).

**Dynamic proof, run on the shipped object** (`v27_measure.recovery_proof`):
200 random table rows x random columns x 50 solves, for each j in {0,1,2,3}
(40,000 exact solves). Each recovered y is pushed through the SHIPPED two-step map
`((y·m1)^2 · m2)^2` and compared to `sigma^j(y)` as group elements:

    j=0   10,000 solves   0 misses    old cubing == corrected (regression kept)
    j=1   10,000 solves   0 misses    old cubing 0/10,000 returns
    j=2   10,000 solves   0 misses    old cubing 0/10,000 returns
    j=3   10,000 solves   0 misses    old cubing 0/10,000 returns

Every recovered y also satisfies `y^r = 1` (subgroup membership), 0 failures. So
the uniqueness argument the theorem needs — exactly one y per multiplier pair for
every class index — is now true per j, not just at j=0.

**FIX 3b — the cycle census no longer collapses representatives.** v26 keyed
states by `dlog[min(keys12)]`, one node per class, so twelve sigma-partners shared
a node. v27 keys by the full twelve-value signature:

    distinct exact states 9,420 (was 785)   collapse factor 12.00 (disclosed)
    v27 squared step: 4 cycles, 0 two-cycles out of 9,420 nodes
    v13 control rule: reported as measured, not as expected

The family the plain step carries is a 2-cycle; on the exact 9,420-node map the
squared step has none.

## 4. FIX 1 — the DP predicate runs at census time, and the model has its factor

v26 stored and checked EVERY visited class key and called any key repeat a
success. Production does not: it stores a report only when the distinguished-point
predicate fires (density 2^-80) and restarts the trail at each report. v27 runs
that protocol literally: `routing.dp_hit(k, dp_bits)` is evaluated at census time
(source line in `v27_measure._dp_batch`: `if not routing.dp_hit(k, dp_bits):
continue`), the store receives reports only, a report restarts the trail, and
success is **DP-hit AND agree**, verified against the planted secret.

Count row, same seeds, four report depths (theta = 2^dp_bits), total n = 180,000:

| dp_bits | theta | S/theta | stored keys/trial | DP hits | reports/trial | collisions | recovered | p-hat | CP-95 low |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 2  | 16 | 23.63 | 735,725 | 12.26 | 35,207 | 35,207 | 0.5868 | **0.5835** |
| 2 | 4  | 8  | 24.47 | 374,647 |  6.24 | 23,365 | 23,365 | 0.3894 | 0.3861 |
| 3 | 8  | 4  | 25.43 |  97,261 |  3.24 |  7,449 |  7,449 | 0.2483 | 0.2442 |
| 4 | 16 | 2  | 25.81 |  54,508 |  1.82 |  4,503 |  4,503 | 0.1501 | 0.1467 |

`p_agree_given_collision = 1.000000` at every depth (0 zero-denominator, 0 wrong-j
recoveries, 0 label-invariant failures at reports). The DP predicate costs nothing
in AGREEMENT; it costs REPORTS, and the report count is what the fixture is short
of.

**Why the production regime is the S/theta -> 0 end, and why the floor is charged
for the difference.** A merge is credited only if the follower reaches the common
downstream report. At production the trail cap is 128*theta and the distance to
the next report is Exp(theta), so

    delta = P(miss | merge) <= E[Exp(theta)] / (128*theta) = 1/128 = 0.0078125,

and S is sized so that `(1-delta)*(1-exp(-S^2(1-f^2)/(2 mu N))) >= 0.39`. The
non-equivalence the judge named is CHARGED IN S, not argued away. In the fixture
the budget (32 steps) is the binding cap, not 128*theta, so the small-S/theta arms
truncate trails early — that is a fixture limitation of the 785-element group and
is disclosed as such, not offered as a production rate. The production-like arm
(S/theta = 16, the largest the fixture admits) reads 0.5868 with CP-low 0.5835.

**The model, with the missing factor** (`v27_measure.model_v27`):

    model = 1 - exp( -S^2 (1-f^2) / (2 mu N) )

Reconciliation row, at the fixture shape S=32, N=785, f=3/8, mu=1.1:

    v27 model        0.399237
    judge recompute  0.399       (agrees to 2.4e-4)
    v26 code         0.479117    (the defect: factor omitted)

So the 0.003 "model agreement" v26 claimed was an artifact of computing a
different formula than the one disclosed. The model now matches the judge's own
arithmetic.

All-key control (v26 protocol, n = 30,000, kept for comparison): 0.4761
[0.4713, 0.4809]. The DP-filtered production-like arm EXCEEDS it (0.5868), because
restarts add fresh independent starts.

## 5. FIX 2 — the trap row charges the true trail population

v26 charged four abandon probabilities. The restart schedule starts one trail per
report interval: `n_trails = ceil(S/theta) = 62,484,585,866,263` = 2^45.83 trail
starts, and every one can abandon. Row by row, from the shipped config:

    L = min(Exp(theta), 128 theta)          E[L^2] = 2 theta^2 (1 - e^-128(1+128))
    P(self-collide | L) = E[L^2] / (2N)     N = 2^251.2721
    wasted work per abandon <= E[L] = theta
    trap27 = n_trails * P * theta * c_mul = 2^35.3257      (v26: 2^11.497, 4 trails)

The population was understated by 2^43.83. The corrected charge is still 2^90.9
below the step term, so it does not move the 4th decimal — but it is now the right
population, and the row is written as an integer ledger term, not a rate.

**Overflow semantics.** Because a trail restarts at its first report, a trail
stores AT MOST ONE report; the report count is deterministically
`<= ceil(S/theta)`. v26's "Chernoff overflow bound" was an expected tail labelled
as a cap. It is DELETED; the report rows are charged at the literal cap
`ceil(S/theta)`, and the aux rows keep their conservative slack-2 ceiling.

## 6. The ledger: one sum, every term displayed, cost rounded up

    w_stress          2^126.267148    S * c_stress
    trap_charge       2^35.325690     trail-start population, corrected
    aux_sum (8 rows)  2^56.517206     every non-step row, at the report cap
    --------------------------------------------------------------
    total             2^126.267148    claim 126.2672 = 78917/625 (reduced, gcd 1)

S is the least integer with `(1-1/128)(1-exp(-S^2(1-f^2)/(2 mu N))) >= 0.39`;
f = 0.346631 is the EXECUTED squaring fraction over all four census arms (declared
3/8 = 0.375), priced at the more expensive reading as before. Fractions
throughout, no float in the claim, claim rounded UP to 4 decimals and reduced.

## 7. What is NOT claimed

No rank claim: 126.2672 loses to 125.82. No claim that the fixture measures a
production-width recovery rate: it cannot, and section 4 says exactly which arm is
production-like and why the others are not. No claim that the squared step is
optimal: it is the minimal sound fix for the orbit table. No claim that a
history-free 125.77-class package exists: the v25 theorem says it cannot.
