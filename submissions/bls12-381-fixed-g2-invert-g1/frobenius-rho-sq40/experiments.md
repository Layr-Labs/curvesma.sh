# Experiments in detail (`experiment.py`, standard library only)

This file carries the per-part detail behind proof.md section 7. Every figure quoted here is
printed by the packaged source, which is one standard-library file; the longer runs behind
`measure/experiment-results.md` use the same source at seed 1 with a larger trial target, and
their counts differ on a rerun.

Every per-step decision is read from the Montgomery-domain value of the class key,
T0M = T0 2^384 mod p and K in the same form, the representation the measured kernels hold.
experiment.py forms each of the twelve keys as key*R mod p before comparing them, takes kappa
over those values, and tests bits 70..197 of the resulting Montgomery-ordered minimum; the toy
fields do the same with their own Montgomery radix R_t. This is not distribution-preserving:
K is a minimum of twelve keys, so its bits are not uniform, and the experiment measures them
rather than asserting uniformity (proof.md section 6).

Part A (BLS12-381 itself): on random GT elements, checks y^r = 1; the Frobenius formula
equals y^p; the four-multiplication key formula equals the keys of all twelve conjugates and
they are distinct; the sign-set minimum equals the twelve-way argmin kappa and the class
key; kappa selects the same canonical element from all twelve conjugates; p^6 is inversion;
ord_r(p) = 12; 1 + p^4 + p^8 = 0 mod r; the factorisation of r-1 and ord_r(2) = (r-1)/16
with the smallest t a fruitless cycle could contain, 2^247.27; T0 equal on all twelve
conjugates of the Montgomery value; and cyclotomic squaring equal to y*y on GT elements.

Part A2 (bit frequencies, read in the Montgomery domain): up to 2,000 time-limited samples,
each a product of three random GT elements; for the low 76 bits of T0M (10 partition, 6
branch, 60 distinguished) and for the report slice, bits 70..197 of the MONTGOMERY-ordered
class key K (the twelve keys formed as key*R mod p and minimised in that form, as the kernels
hold them), it reports one-bit frequencies and their largest deviation from 1/2 in standard
errors, plus the fraction of elements whose branch decision changes between the canonical and
the Montgomery reading. A run at the declared inputs gave largest deviations 3.40 (T0M) and
4.83 (K slice 70..197) standard errors (means 1.18 and 1.12), branch-field fraction 0.5475
against 0.546875, low-4-zero fraction 0.0585 against 0.0625, and 0.3025 of the elements on a
different branch under the two readings. The K-slice figure is a measurement of a nonuniform
quantity (K is a minimum of twelve keys), not a proof that those bits are uniform; proof.md
section 6 takes them as pseudo-random under H1. It cannot measure a 2^-60 event; the toy
measures that instead.

Part P (the production step rule on real GT elements): two walkers start at two conjugates
y0 and y0^{p^5} of one GT element and run the step rule in lockstep with a 16-entry table,
every field read from T0M, keys only on multiply steps and at attempts, and every kappa and
class key taken from the twelve Montgomery keys (key*R mod p), as the kernels produce them.
They must produce the same T0M sequence, keep y = g^alpha z^beta, and report key computations
equal to multiply steps plus attempts. A 4,000-step run gave squaring fraction 0.5313
(2,124 + 1,874), 2,103 key computations against 1,874 multiply steps and 229 attempts,
identical T0M sequences and the invariant holding for both walkers.

Part P also runs two constructed cycles. A two-cycle (all-multiply, built by setting the two
table cells it uses so that y = g^alpha z^beta still holds) is detected by the two-step ring:
one walker holding it in its ring escapes at once, one started fresh at the same class after
two steps, they escape at the same class (the smaller-T0 class squared) and share their
class sequence for the 10 steps compared, with different exponent pairs at the first common
class, which is a usable collision. A three-cycle, invisible to the ring, is detected only by
the checkpoint, and the checkpoint hit runs the declared escape in full: the walker traverses
the cycle with the ordinary step rule, records the class it came from at each step, takes the
class with the smallest T0M, squares that class's element, empties its ring and checkpoint
slot and continues. A walker whose phase puts a checkpoint inside the cycle at once detects
it at iteration 3, a walker that must wait for its next multiple of 64 detects it at
iteration 66, inside the 64 + 3 that section 3.1 allows; both traverse the same cycle, both
escape at the min-T0M class squared, and from there their class sequences agree for the 6
steps the run compares after the escape. The test is detection timing and escape action, not
a production-scale cycle rate.

Part E (the collision law, measured directly). In the toy field q = 67 distinct classes are
sampled uniformly, one step is applied to each, and image coincidences are counted by pair
type. At q = 67 T0 has 7 bits, too few for the fields, so the bits come from a hash of the
class (a function of the class alone, as the production fields are); this is the one place
the toy's bits are not the field's own. The pooled statistic N x coincidences/pairs has
expectation 1 - s^2 - (1-s)^2/(12R); the charged 1 - s^2 - (1-s)^2/R is reported next to it,
and squaring-squaring coincidences must be exactly zero.

Part B (toy, same tower shape): q = 7 mod 12 near 2^20, r' a prime factor of q^4 - q^2 + 1
with ord_r'(q) = 12, N = (r'-1)/12; the production walk of section 3 at toy scale (T0 has 20
bits, so the fields fit, and the bits are read from the toy field's Montgomery form), every
distinguished hit counted as an attempt and charged 3,000 instructions to the slot's meter,
the attempt cap carried per slot, every trial run to its full fixed budget 1.05 S0' with the
section-5 meter and no stop at the first recovery, then the production matcher (full-record
sort, adjacent dedup, charged pairs against a hard budget, canonical recomputation,
verification of g^x = z). Because the toy's distinguished fraction is 2^-4 or 2^-3 rather
than 2^-60, its attempt charge is about 1.2% or 2.4% of the step price: the accounting is
tested roughly 2^56 times more harshly than production. A measurement-only class index
timestamps each trial's first useful collision; no success flag reads it. Stress passes
exercise the attempt cap (four slots stopped at a 4-attempt cap with the log still within the
slots' attempts), the stuck-walk restart with the carried slot counter, the meter, and the
pair budget.

Part C recomputes, without randomness: the zero-sum probabilities, the fruitless-cycle
bounds, N' and S0', the deterministic log capacity, the KL attempt-cap bound, the attempt
charge inside the meter and inside the capacity group, the detection-window offset, the
Montgomery key ordering and the cycle allowance stated as a declared allowance, the memory
bound, the per-record report charge, the pair budget, the
meter's exact integer cap, the dependent-process ledger, the H1+ support figures (headroom in standard deviations and over the cycle mean), the
counter table, and the W ledger with the exact integer claim check and the 2^57-unit
preprocessing bound. It re-derives the measured squaring price from `measure/results.tsv`
and the six raw Callgrind records and fails if any of them differs.

**Part E and part B results** are tabulated in `measure/experiment-results.md`. Part E
(2,500 distinct classes, q = 67, R = 64) measures the pooled statistic 0.9095 +- 0.1697
against 0.6969 predicted at the sample's own branch mix and 0.6977 charged (1.25 standard
errors, a toy-size sample), with
squaring-squaring coincidences exactly zero as injectivity requires; it cannot resolve the
2.624 x 10^-4 factor difference at R = 1024, and does not need to, because the charge is the
larger N'. Part B (1,000 trials per field, 3,000 in all) gives pooled P[T_c <= 1.05 S0'] = 1,224/3,000 =
0.4080 +- 0.0090 against the model 0.42013 at 1.05 S0' (0.971 of it; the toy budget IS the
production 1.05 S0' budget, so this is the model's own scope, not a bound from above; after
the 0.01 allowance the toy leaves 0.3980 >= 0.39), pooled mean T_c/sqrt(N') =
1.2634 +- 0.0120 against sqrt(pi/2) = 1.2533 and 1.0486 for a model without the branch
adjustment (17.9 standard errors away, so the measurement sits on the branch-adjusted scale),
and a recovered-and-verified fraction of 0.4567. The same file has the per-row squaring and
distinguished fractions, the meter's measured per-slot work deviation (0.3933, 0.3933, 0.3974
sqrt(steps) against 0.3877 predicted by independent step costs), the squarings inside the
checkpoint escapes' cycle traversals, the attempt-cap and pair-budget stress rows, and the
toy's 2- and long-cycle rates, with the note that the toy's R/S0' is too large for the
independence regime production satisfies (R/S0 = 2^-115.6).

The packaged run at its initial inputs (a 10 s budget, seed 1) completes smaller samples of
the same measurements, so its part A2 sample, part E class count and part B trial counts
differ on a rerun; they are reported for reproducibility, not as the evidence.

## Source notes (comments and docstrings moved out of `experiment.py` to keep it under the
per-file limit)

* `bits`: (h, squaring?, distinguished?) are disjoint bit fields of T0 (mode t0) or of a
  class hash (toy part E only).
* `walk_pass`: one trial's walk phase in production order (meter, both cycle checks before
  the report append, append rule, mixed step); `instrument` timestamps the first useful
  class collision and no decision reads it.
* Attempt charge: every distinguished hit is an attempt, counted and charged to the slot's
  meter whether the record is stored or suppressed. The toy meter cap carries the worst-case
  attempt charge (cap_walk x I_REP), as production's nominal carries I_REP x R_w.
* `SLOT_METER_MAX = W_CAP - 1 + I_MUL + I_REP` is the largest over-cap update: a multiply
  step that is also a DP attempt. The step counter bound SLOT_STEPS_MAX divides it by the
  cheapest step.
* `SELF_DELTA_LG = -20.96` is the Markov threshold on the wasted fraction of S; the wasted
  slot loss is `n^2/(2N')` (per slot) times the 2^20 slots over that threshold.
* The meter deviation printed by part B removes the attempt charge before taking the
  per-step deviation, so it compares with the independent-step prediction of step costs only.
* Part C's strings and floats are the digits proof.md, README.md, experiments.md and
  submission.json quote: log2 costs and probabilities are rounded up, success figures and
  headroom figures down, and each text literal is equal to the printed value or on the safe
  side of it.
