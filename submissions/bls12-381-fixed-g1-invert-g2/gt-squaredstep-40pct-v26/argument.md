# Argument: the squared-step walk (v26)

## 0. What this package is

A validity play, not a rank play. The claim is **126.2508**, which is 0.4308 bits
behind the live board minimum 125.82 and 0.4128 bits behind this account's own
valid 125.838 row (72a73a40). It is submitted because it is the only sound
history-free family-free walk this track admits: the 125.77-class that v24 tried
is mathematically dead, and the proof is shipped in-tree as
`evidence/impossibility-theorem.md`.

## 1. The object

BLS12-381, fixed P in G1, uniform z in the order-r subgroup of GT. Compute
g = e(P, Q0) with the pinned py_ecc map, solve z = g^x in GT, return Q = [x]Q0.
Everything is classical, single-machine plus threads, no amortisation, no
precomputation across challenges.

## 2. The step, and why it is squared

`classkey.step_v26`, the shipped step, executes the census and the certificate:

    kappa = classkey(y)                       # the 12-element sigma-orbit class
    t     = routing.fields(kappa).row         # one table row from the class only
    y'    = (y * TBL[t][kappa])^2             # multiply, then cyclotomic square
    a'    = 2*(a + U[t][kappa]) mod r ; b' = 2*b mod r

Three properties, each tested in-package:

* **Deterministic function of the class.** No prev_row, no ring, no history. So a
  class repeat IS a state repeat and the cycle abort is sound. (v24's `step_v24`
  used prev_row; its verdict names that as defect 2. The scaffolding is REMOVED
  here, not re-litigated: `step_v24`, the class-repeat abort it invalidated, and
  the sampler that used them are gone from the shipped path.)
* **Equivariant.** `(sigma^j y * sigma^j m)^2 = sigma^j((y m)^2)`, so class
  merging is preserved: partner elements stay merged.
* **Family-free.** A plain multiply step CANNOT be family-free: class merging
  forces `m(k) = sigma^(-k)(A_t)`, hence `m(k+6) = m(k)^(-1)` in every column
  (measured 384/384 on the fixture, with a planted mutation control that makes the
  scan flag exactly the mutated column). Squaring is the minimal fix: the two-step
  return becomes `y^3 = sigma^j/(m1^4 m2^2)`, one solution per multiplier pair.

## 3. Price, measured not estimated

The multiply step gains exactly one Granger-Scott cyclotomic squaring per multiply
step. The census counts the calls on the shipped object: 510,194 mul steps and
510,194 cyclo_sqr calls in the same run (`equality_control.counts_match`), and the
instrumented re-run reproduces the plain run's per-trial outcome vector element for
element.

    mul_step 21889 -> 33521 (+11632)      sqr_step 13445 (already a square)
    c_mul_stress 1.122799 -> 1.704385     c_stress 0.964032 -> 1.346817

The squaring fraction is priced at the EXECUTED 0.355811 (declared 3/8), i.e. the
more expensive of the two readings: with f = 0.355811 the claim is 0.004 bits
HIGHER than with 3/8, and the higher number is the one submitted.

Delta against v24's claim: **+0.4824 bits** (measured), against v25's +0.473
estimate. The measured delta is larger because the executed squaring fraction is
below the declared 3/8 and because the cyclo squaring is priced as a counted row
(24*FP_MUL + 122*FP_ADD) rather than being multiplied by the 1.3 stress factor.

## 4. Success: unconditional product, measured at production shape

The census protocol mirrors production: S_toy = `constants.step_budget(N_toy)`, the
SAME function that sizes the production budget; 4 trails of ceil(S_toy/4) steps;
per-worker seeds (distinct start streams; the v24-accepted fix kept); one shared
store per trial; every recovery verified against the planted secret; UNCONDITIONAL
statistic per trial, i.e. a single-run product P(repeat reached) x P(agree |
repeat), not a conditional.

n = 30,000 trials:

    P(repeat reached) = 0.4761        P(agree | repeat) = 0.999930
    P(recover)        = 0.4761        Clopper-Pearson 95% [0.4713, 0.4809]
    model 1-exp(-S^2(1-f^2)/(2 mu N)) at the same shape = 0.4791

Measured and model agree to 0.003, and the CP lower bound 0.4713 clears the claimed
0.39 by 0.081. The claim is deliberately BELOW both: 0.39 is the floor the track
prices, and the package claims the floor, not the measurement.

CP vs bootstrap: Clopper-Pearson is used (exact, conservative at every n), stated
here as the method; the binomial interval does not rely on a normal approximation,
and the interval is reported on the JOINT count (recovered/trials), so the product
structure is inside the interval rather than multiplied outside it.

Scope, stated not glossed: the fixture is N = 785 with a 32-step budget, so the
census establishes the merge -> solve -> verify mechanism AT THE PRODUCTION
RATIO S/sqrt(N) and calibrates the model there; it cannot establish a
production-width recovery rate, and the claim does not rest on it doing so.

## 5. The ledger: one sum, every term displayed, cost rounded up

`ledger26.total26()` returns the components it sums; `claim26()` asserts the claim
is >= log2 of that same sum, and asserts the components sum exactly to the total
(Fraction equality, no float):

    w_stress             2^126.250778    S * c_stress
    trap_charge          2^11.497        per-trail birthday abandon charge
    aux_sum (8 rows)     2^56.510        every non-step row, charged at the CAP
    aux_overflow_bound   2^0             exp(-2^45.8/3) ceiling, charged 1 ref op
    ---------------------------------------------------------------
    total                2^126.250778    claim 126.2508 (least 4-dp log2 >= total)

* **The v24 reserve is DELETED.** v24 carried 2^114.812 of "reserve" that no
  displayed row accounted for (its verdict called it padding). Here every row is
  enumerated, and the completeness question cannot move the claim: the whole aux
  block is 2^56.5, seventy bits below the step term, so even a factor-2^60 error
  in aux changes the 4th decimal of the claim by nothing. The bound that replaces
  the reserve is constructed, not reserved: per-report rows are charged at the CAP
  (2*ceil(S/theta)), and exceeding the cap has Chernoff probability
  exp(-mean/3) with mean = S/theta = 2^45.8.
* **N = 2^251.2721** = log2((r-1)/12), verified in constants, used everywhere. The
  v24 note's N = 2^371 is the field size and is gone; the trap row here is
  t(t-1)/2N with t = cap = 128*theta, i.e. 2^-78.27, charged PER TRAIL (4 trails),
  because an abandoning trail wastes its own cap independently.
* **Ceilings everywhere.** cost rounds up, claim rounds up, the squaring fraction
  is priced at the cheaper-for-us-forbidden reading, and Python integers carry no
  float into the claim.

## 6. What is NOT claimed

No rank claim: 126.2508 loses to 125.82 by 0.4308 bits, disclosed in the statement
and in the note. No claim that any history-free 125.77-class package exists: the
theorem says it cannot. No claim that the squared step is optimal: it is the
minimal SOUND fix for the orbit table, and the price is what it is.
