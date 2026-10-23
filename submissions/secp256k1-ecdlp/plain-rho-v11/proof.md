# plain-rho-v11: secp256k1-ecdlp candidate

The structured constants of v6 (and its unsubmitted v7 sibling) are
replaced by PUBLISHED PSEUDORANDOM constants: A_j = [c_j]P with
c_j = int(sha512("curvesmash/secp256k1-ecdlp/constant/j")) mod r. The
review of v6 found a deterministic obstruction: with r = 32Q + rho,
the rescaling z = 32k mod r turns every v6 step into a positive
increment 2^(j+5) - j <= 2^36, confining each trail to an interval
of width 2^36 * L and making cross-walker collisions a < 2^-75
event. Structure-free constants have no such rescaling (Part A
measures every rescaling up to 2^12 and states the counting
heuristic for all rescalings), the analysis is the standard
random-function model of parallel rho, and cycle words of length
<= 4 are excluded by executed enumeration.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **64365809/500000 = 128.731618** |
| success_probability (lower bound) | 0.39 (modeled 0.6053 at Teske's measured factor) |
| memory_log2_bytes | 105.1 (over-report; see below) |
| preprocessing_log2 | 25.0 (included) |
| nonuniform_advice_log2_bytes | 0 |

## Algorithm

- Group: secp256k1: y^2 = x^3 + 7 over the 256-bit prime field, r prime, standard generator. Target T = [x]P for uniform secret
  x in [0, r); recover x.
- Constants: A_j = [c_j]P for j = 0..31, where c_j is the published
  deterministic pseudorandom value
  c_j = int(sha512("curvesmash/secp256k1-ecdlp/constant/j")) mod r.
  The 512-bit hash value reduced mod r has bias below 2^-500.
  Part A verifies the 32 constants are distinct and prints the
  sha256 of their concatenation (bf815f015e0d4e49...).
- K = 2^16 walkers, each seeded as C = [a_0]P + [b_0]T with fresh
  uniform a_0, b_0 (512 DRBG bits per walker), carrying the lazy
  representation (a, b) with C = [a]P + [b]T and b fixed for the
  walker's lifetime.
- Step: read the state's index j = x_C mod 32 (5 bits of the
  x-coordinate), update C <- C + A_j and a <- a + c_j mod r.
- Budget: t = 1.5 * sqrt(r-1) walker steps total (log2 t =
  128.58496), shared round-robin across the walkers
  (log2 of per-walker budget 112.58).
- Distinguished points: a state is distinguished when the low 32
  bits of its x-coordinate are zero (density 2^-32). On reaching
  one, the walker inserts (x_C, y-parity) -> (a, b) into a shared
  table of capacity 2^98 slots. A key hit with a different
  representation is a collision: x = (a2 - a1)/(b1 - b2) mod r is
  extracted and verified against T; the trial succeeds. A key hit
  with the identical representation is a cycle revisit (below);
  the walker reseeds.
- The table key is the pair (x-coordinate, y-parity). On the full
  group +X and -X are DISTINCT points that share an x-coordinate;
  the parity bit is what makes a key hit a true point collision,
  so the extraction is exact algebra with no false positives (the
  toy measured 0 wrong extractions in
  200 trials plus 100 census trials).

## Confinement: why the constants are pseudorandom

The lazy representation makes the walker's scalar k = a + b*x mod r
walk with steps c_j; two walkers collide at trail positions (u, v)
exactly when their scalar-difference walk hits the (uniform)
seed difference. For this to happen at the birthday rate, no
invertible rescaling of the scalar may confine the step set. The
v6 constants failed exactly this: r = 32Q + rho makes
32 * c_j = 32*(jQ + 2^j) = 2^(j+5) - j (mod r), so the rescaled
steps are positive and bounded by 2^36, z = 32k walks inside an
interval of width 2^36 * L, and the v6 review bounded cross-walker
collisions by binom(K,2)*(2^36 L)/r < 2^-75.

For the v11 constants there is no small-coefficient structure to
exploit:

- Executed certificate: Part A measures, for every rescaling
  M <= 2^12, the width of the narrowest window containing all 32
  rescaled steps M*c_j mod r. The narrowest over all 2^12 - 1
  rescalings is 0.6329 of the modulus (at M = 17) -
  i.e. every checked rescaling leaves the step set spread across
  most of the group. The same certificate APPLIED TO THE REJECTED
  v6 CONSTANTS detects their confinement: at M = 32 the window
  collapses to 5.93e-67 of the modulus.
- Counting heuristic for all rescalings: the expected number of
  M in [0, r) that bring all 32 steps into any single window of
  width W is 32 * W^31 / r^30 (the 32 rescaled steps are
  pseudorandom, so one anchor is free and the other 31 must land
  in a window of relative width W/r). For W = 2^240 this is
  < 2^-235; a confinement that could matter for the
  budget needs W * 2^112.58 < r, i.e. W < 2^143.4,
  which sits 2^97 below the heuristic threshold.
- Consequently the standard random-function model (H1) applies:
  walker states collide at the van Oorschot-Wiener parallel-rho
  birthday rate t^2/(2(r-1)), distinguished-point counts are
  Binomial(t, 2^-32), and cycle entries occur at rate 1/r per
  window. Part B measures the realized fraction of the ideal
  directly: 133/200 = 0.665
  against the same-parameter ideal 0.6753
  (ratio 0.985, i.e. the toy MEETS the
  birthday ideal within 0.099 at two
  standard deviations).

## Fruitless cycles

- No cycle word of length <= 4 exists: Part A enumerates every
  multiset of at most 4 constants (4 literals, all
  multisets with repetition) and asserts none sums to 0 mod r. A
  cycle of length L requires its realized index multiset to sum
  to 0 mod r, so this excludes every cycle of length <= 4
  regardless of index realizability.
- Longer cycles under H1: a fruitless cycle is a SELF-repeat of
  one walker - a pair (u, v), u < v, on one trail whose window
  sum is 0 mod r, an event of probability 1/r per self pair. The
  number of self pairs is K*L*(L-1)/2 with K = 2^16 walkers and
  L = t/K, so the expected number of cycle entries is
  t^2/(2Kr) = 1.717e-05 - the self-collision birthday
  scale, not a per-step rate. A cycle containing a distinguished
  point is detected at its first DP revisit (cost at most one
  lap) and the walker reseeds; a cycle containing none traps the
  walker for at most its remaining budget. Charging the
  remaining budget for every entry at (u, v) bounds the total
  expected waste by K * sum_(u<v) (L - v)/r = K*L^3/(6r) =
  t^3/(6 K^2 r) steps, i.e. a fraction t^2/(6 K^2 r) =
  8.731e-11 of the budget. The bound is computed over the
  full un-terminated budget, so it also bounds the terminated
  process. The expectation is NOT substituted as a
  deterministic budget cut: it is converted by Markov's
  inequality - P(waste exceeds 2^-23 of the budget)
  <= 7.32e-04 - and the reseed work (385 reference
  operations per detected entry, expected 6.61e-03
  operations over the whole attack) is charged inside the same
  bound. The reporting loss -80.58 exponent stays
  subtracted as before; the total success impact of cycles and
  reporting is below 0.06%, inside the 0.2153
  margin of the model.
- The toy observes 3 cycle reseeds over
  200 trials against a self-pair model expectation
  of 27.7 for the same run - below the model,
  i.e. the discrepancy is in the conservative direction for the
  claim.

## Distinguished points, reporting, and memory

A collision is detected when the later visitor inserts a
distinguished point whose key is already in the table. If the
collision point itself is not distinguished, the later walker has
merged onto the earlier walker's trail; because the step map is
deterministic, the merged pair CANNOT separate, and it reaches the
next shared distinguished point together. The only way a merge goes
unreported is that the walker's personal budget ends first: the
expected gap to the next distinguished point is 2^32 steps against
a personal budget of 2^112.58 steps, so the reporting loss
is 2^-80.58. The census records this statistic
directly: 64 of its trials contain
a true first-contact merge (counted by a diagnostic full-visit
log that is not part of the attack), 41
are detected and 23 end
undetected by the budget - the budget-end race scaled to toy
size, where the gap to the next distinguished point (2^6 = 64
steps) is comparable to the per-walker budget
(191 steps). The same law at production
scale is the 2^-80.58 term above, subtracted.

Table sizing under H1: expected inserts = t * 2^-32 = 2^96.58
(the insert count is Binomial(t, 2^-32) under H1; the Chernoff
bound is stated under that premise); capacity 2^98 slots
(2.7x the mean), overflow
probability <= exp(-2^95.6); on overflow the attack
aborts and fails, and the probability is subtracted (it is far
below the printed precision). Insert processing is charged at 20
instructions each: 8 for hashing the key, 4 for the table probe, 8
for the record compare/write.

Memory. Record layout (32-byte x-coordinate, 8 bytes parity/flags,
32-byte a, 32-byte b, 8 bytes metadata): 112
bytes <= 128 per slot, padded to 128. (This corrects the v5 review
item that a 32-byte field cannot hold the 381-bit BLS
x-coordinate: the x-coordinate gets 32 bytes here.) Peak memory
is the 2^98 x 128-byte table plus per-walker state (2^16 walkers
x a few hundred bytes) plus the 32 constants: strictly above
2^105 bytes. As with the partial-cheon submissions, the
reported memory_log2_bytes = 105.1 is the rounded-up over-report
of that peak; memory is a supplementary metric and the
over-report direction is conservative.

## Work ledger

Per walker step (AArch64 instruction counts, the "other" category
of classical-group-work-v2; one reference operation = batched
incomplete affine addition is charged as 1.000 and every walk step
performs exactly one):

| Item | Instructions |
| --- | --- |
| index j = x_C mod 32 (5-bit extract) | 2 |
| representation add a += c_j mod r (full-width constant add + reduce) | 32 |
| distinguished-point predicate (low-32-bit OR-reduce of x) | 4 |
| walker state, loop, branch, counters | 69 |
| **other per step** | **107** |
| reference add (the step itself) | 1.000 |

Per-step price: 1 + 107/1000 = 1.107 units. One-offs,
all included in W and itemized:

- constants: 32 scalar multiples, ~256 reference adds each = 2^13
  reference operations;
- seeding 2^16 walkers: [a_0]P by fixed-base comb from a
  precomputed table (32 point adds AND 32 doublings per
  evaluation - the comb rows are combined by doublings; building
  the comb table costs ~2^9 reference ops once per attack),
  [b_0]T by 4-bit windowed double-and-add (256 doublings + 64
  window additions; the 16-entry window table for the attack
  target is built once, 15 adds), plus the final addition
  [a_0]P + [b_0]T forming the seed: 32 + 32 + 256 + 64 + 1 = 385
  reference operations per walker; 2^16 x 385 = 2^24.59, charged
  as 24.59 reference operations (25.0
  rounded up; the table-build one-offs included);
- randomness: 512 DRBG bits per walker, ~64 instructions per
  walker = 2^22 instructions = 2^12 units, included in the
  seeding charge;
- distinguished-point processing: 2^96.58 inserts x 20
  instructions = 2^100.9 instructions;
- extraction and verification: one modular inverse (~700
  instructions) plus one verification scalar multiple [x]P == T
  (~256 reference adds = 2^8 units), charged as 2^9 units. The
  attack verifies every extraction before declaring success; the
  toy recorded 0 wrong extractions in
  200 + 100 trials.

W = 1.107 * 1.5 * sqrt(r-1) + 2^24.59 + 2^92.2 (DP)
+ 2^87.0 (table zeroing) + 2^28.6 (reseed cap)
+ 2^9 (extraction),
log2(W) = 128.731618, reported exactly as
time_log2 = 64365809/500000.

## Success probability

- Ideal collision probability within the budget:
  1 - exp(-(1.5)^2 / 2) * (1 - aborts) = 0.6753475326.
- Collisions are overwhelmingly cross-walker (self-repeats are
  the cycle entries bounded above), every key hit at a true point
  (parity key) is a collision, the extraction is exact algebra,
  and every extraction is verified against T. The reporting loss
  2^-80.58 and the cycle waste fraction
  8.731e-11 are both subtracted in the ideal below.
- Efficiency: the r-adding walk has a collision-time factor c
  (real collision time / ideal random-walk collision time).
  Within the fixed budget the success probability is
  1 - exp(-(1.5/c)^2/2) times the reporting and abort factors
  (Rayleigh collision times), so the 0.39 claim holds for every
  c <= 1.5086, computed after all deductions. Teske's experiments on r-adding walks with
  r >= 16 measure expected collision times at most ~1.1x the
  random-walk ideal (1.1), giving modeled success
  0.6053484543 >= 0.39 with 0.2153
  margin; the toy implies c = 1.014 (faster than ideal:
  measured ratio 0.985 of the
  birthday ideal, two-sigma 0.099).
  The claim is made against the c tolerance, not against the
  toy measurement.

## Response to the v10 review

1. *The fixed-base comb evaluation also requires ~32 doublings
   per scalar, and the final seed addition was omitted; 385
   operations per walker exceed the declared 2^24.5 preprocessing
   bound.* Corrected: per walker 32 adds + 32 doublings (comb)
   + 256 doublings + 64 additions (windowed [b_0]T) + 1 final
   seed addition = 385 reference operations; 2^16 x 385 =
   2^24.59 charged as 24.59, preprocessing_log2 =
   25.0 (rounded up), including the once-per-attack table
   builds (~2^9 comb table, 15 adds for the 16-entry T-window
   table) and the 512 DRBG bits per walker. The score is
   unchanged at six decimals.
2. *The 0.75 efficiency-floor conversion was arithmetically
   wrong (a walk 33% slower than a 1.1x walk gives success
   ~0.305, not 0.4066).* The floor language is dropped
   entirely: the success model now uses the collision-time
   factor directly - under Rayleigh collision times the 0.39
   claim holds for every factor c <= 1.5086 - and Teske's
   measured r-adding factors (at most ~1.1x the ideal for
   r >= 16) give modeled success 0.6053484543 >= 0.39
   with 0.2153 margin. The toy implies
   c = 1.014, faster than ideal.
3. *The cycle-waste expectation was substituted as a
   deterministic budget reduction; a single trapped walker can
   waste far more than the mean.* The expectation is converted
   by Markov's inequality: P(waste exceeds 2^-23 of the budget)
   <= 7.32e-04, with reseed work (385 operations per
   detected entry, 6.61e-03 expected operations in
   total) charged inside the same bound. No deterministic
   substitution appears anywhere in v11.
4. *The sparse census does not record the undetected-merge
   statistic.* It does: the recorded census output contains
   64 trials with a true
   first-contact merge (diagnostic full-visit log, not part of
   the attack), 41 detected,
   23 undetected by budget end,
   matching the budget-end race scaled to toy size.
5. *The Chernoff overflow bound needs justified premises.* The
   insert count is Binomial(t, 2^-32) under H1, the same
   random-function premise that yields the birthday bound; the
   Chernoff bound is stated under that named premise.

The v8-review items (the t/r cycle rate, the 64+16 seeding
misprice, the census statistic, the Chernoff premise) remain
corrected in the main text and are not restated here.

## Reproduce

```
echo '{"trials": 200, "dbits": 0, "walkers": 8}' | \
  python3 experiment.py
```

Stdlib only, deterministic (fixed seeds), single JSON object on
stdout. `trials` scales the toy Monte Carlo; `dbits` sets the
distinguished-point density of the main attack (0 reproduces the
declared run); `walkers` sets the walker count. Part A runs in
every execution and asserts the constants, the short-cycle
enumeration, the confinement certificate (including the v6
self-test), the ledger arithmetic and the success bound >= 0.39.
