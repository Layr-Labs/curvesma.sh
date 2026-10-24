# Experiments in detail (`experiment.py`, standard library only)

This file carries the per-part detail behind proof.md section 7. Every figure quoted here is
printed by the packaged source, which is one standard-library file; the longer runs behind
`measure/experiment-results.md` use the same source at seed 1 with a larger trial target, and
their counts differ on a rerun.

Every per-step decision is read from the Montgomery-domain value of the class key,
T0M = T0 2^384 mod p and K in the same form, the representation the measured kernels hold.
The toy fields use their own Montgomery form (a R_t mod q); the map is a bijection on F_q,
so this changes no distribution.

Part A (BLS12-381 itself): on random GT elements, checks y^r = 1; the Frobenius formula
equals y^p; the four-multiplication key formula equals the keys of all twelve conjugates and
they are distinct; the sign-set minimum equals the twelve-way argmin kappa and the class
key; kappa selects the same canonical element from all twelve conjugates; p^6 is inversion;
ord_r(p) = 12; 1 + p^4 + p^8 = 0 mod r; the factorisation of r-1 and ord_r(2) = (r-1)/16
with the smallest t a fruitless cycle could contain, 2^247.27; T0 equal on all twelve
conjugates of the Montgomery value; and cyclotomic squaring equal to y*y on GT elements.

Part A2 (bit frequencies, read in the Montgomery domain): up to 2,000 time-limited samples,
each a product of three random GT elements; for the low 76 bits of T0M (10 partition, 6
branch, 60 distinguished) and the low 128 bits of the class key K it reports one-bit
frequencies and their largest deviation from 1/2 in standard errors, plus the fraction of
elements whose branch decision changes between the canonical and the Montgomery reading. A
run at the declared inputs gave largest deviations 3.40 and 3.13 standard errors (means 1.18
and 1.01), branch-field fraction 0.5475 against 0.546875, low-4-zero fraction 0.0585 against
0.0625, and 0.3025 of the elements on a different branch under the two readings. It cannot
measure a 2^-60 event; the toy measures that instead.

Part P (the production step rule on real GT elements): two walkers start at two conjugates
y0 and y0^{p^5} of one GT element and run the step rule in lockstep with a 16-entry table,
every field read from T0M, keys only on multiply steps and at attempts, cycle checks on T0M.
They must produce the same T0M sequence, keep y = g^alpha z^beta, and report key
computations equal to multiply steps plus attempts. A 4,000-step run gave squaring fraction
0.5480 (2,191 + 1,807), 2,059 key computations against 1,807 multiply steps and 252 attempts,
identical T0M sequences and the invariant holding for both walkers.

Part P also runs two constructed cycles. A two-cycle (all-multiply, built by setting the two
table cells it uses so that y = g^alpha z^beta still holds) is detected by the two-step ring:
one walker holding it in its ring escapes at once, one started fresh at the same class after
two steps, they escape at the same class (the smaller-T0 class squared) and share their
class sequence, with different exponent pairs at the first common class, which is a usable
collision. A three-cycle, invisible to the ring, is detected only by the checkpoint: a walker
whose phase puts a checkpoint inside the cycle at once detects it after one traversal, a
walker that must wait for its next multiple of 64 detects it 63 steps later, inside the
64 + 3 that section 3.1 allows, and both detect the same cycle.

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
the attempt cap carried per slot, every trial run to its full fixed budget 1.03 S0' with the
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
Montgomery bijection, the memory bound, the per-record report charge, the pair budget, the
meter's exact integer cap, the dependent-process ledger with both Hoeffding exponents, the
counter table, and the W ledger with the exact integer claim check and the 2^57-unit
preprocessing bound. It re-derives the measured squaring price from `measure/results.tsv`
and the six raw Callgrind records and fails if any of them differs.

**Part E and part B results** are tabulated in `measure/experiment-results.md`. Part E
(2,500 distinct classes, q = 67, R = 64) measures the pooled statistic 0.6615 +- 0.1713
against 0.7100 predicted at the sample's own branch mix and 0.6977 charged, with
squaring-squaring coincidences exactly zero as injectivity requires; it cannot resolve the
2.62 x 10^-4 factor difference at R = 1024, and does not need to, because the charge is the
larger N'. Part B (600 trials per field) gives pooled P[T_c <= 1.03 S0'] = 745/1,800 =
0.4139 +- 0.0116 against the model 0.40808 at 1.03 S0' (the toy budget is 1% above
production's 1.02 S0', so this bounds production from above), pooled mean T_c/sqrt(N') =
1.2314 +- 0.0154 against sqrt(pi/2) = 1.2533 and 1.0488 for a model without the branch
adjustment (11.9 standard errors away, so the measurement sits on the branch-adjusted scale),
and a recovered-and-verified fraction of 0.462. The same file has the per-row squaring and
distinguished fractions, the meter's measured per-slot work deviation (0.3982, 0.4064, 0.4004
sqrt(steps) against 0.3877 predicted by independent step costs), the attempt-cap and
pair-budget stress rows, and the toy's 2- and long-cycle rates, with the note that the toy's
R/S0' is too large for the independence regime production satisfies (R/S0 = 2^-115.6).

The packaged run at its initial inputs (a 10 s budget, seed 1) completes smaller samples of
the same measurements, so its part A2 sample, part E class count and part B trial counts
differ on a rerun; they are reported for reproducibility, not as the evidence.
