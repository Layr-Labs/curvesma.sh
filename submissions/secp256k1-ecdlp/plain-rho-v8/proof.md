# plain-rho-v8: secp256k1-ecdlp candidate

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
| time_log2 | **16058573/125000 = 128.468584** |
| success_probability (lower bound) | 0.39 (modeled 0.4066) |
| memory_log2_bytes | 105.1 (over-report; see below) |
| preprocessing_log2 | 23.0 (included) |
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
- Budget: t = 1.25 * sqrt(r-1) walker steps total (log2 t =
  128.32193), shared round-robin across the walkers
  (log2 of per-walker budget 112.32).
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

For the v8 constants there is no small-coefficient structure to
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
  budget needs W * 2^112.32 < r, i.e. W < 2^143.7,
  which sits 2^96 below the heuristic threshold.
- Consequently the standard random-function model (H1) applies:
  walker states collide at the van Oorschot-Wiener parallel-rho
  birthday rate t^2/(2(r-1)), distinguished-point counts are
  Binomial(t, 2^-32), and cycle entries occur at rate 1/r per
  window. Part B measures the realized fraction of the ideal
  directly: 114/200 = 0.57
  against the same-parameter ideal 0.5420
  (ratio 1.052, i.e. the toy MEETS the
  birthday ideal within 0.129 at two
  standard deviations).

## Fruitless cycles

- No cycle word of length <= 4 exists: Part A enumerates every
  multiset of at most 4 constants (4 literals, all
  multisets with repetition) and asserts none sums to 0 mod r. A
  cycle of length L requires its realized index multiset to sum
  to 0 mod r, so this excludes every cycle of length <= 4
  regardless of index realizability.
- Longer cycles under H1: the realized word of any window is
  uniformly distributed over the 32-letter alphabet, so the
  expected number of cycle entries over the whole walk is
  t/r < 2^-127 (secp) / 2^-126 (BLS). Each entry is detected for
  free by the distinguished-point table as a key hit with an
  identical representation and the walker reseeds; undetected
  entries waste at most the remaining budget, so the total waste
  is below 2^-127 * t - negligible. The claim does not depend on
  the absence of cycles: it survives cycle-entry rates up to
  2^-60 with the success bound unchanged.
- The toy observes 3 cycle reseeds over
  200 trials - below the H1 estimate for the toy
  (about 20), consistent with Teske's analysis that r-adding
  walks with r = 32 partitions realize fewer fruitless cycles
  than the random-function model; the discrepancy is in the
  conservative direction for the claim.

## Distinguished points, reporting, and memory

A collision is detected when the later visitor inserts a
distinguished point whose key is already in the table. If the
collision point itself is not distinguished, the later walker has
merged onto the earlier walker's trail; because the step map is
deterministic, the merged pair CANNOT separate, and it reaches the
next shared distinguished point together. The only way a merge goes
unreported is that the walker's personal budget ends first: the
expected gap to the next distinguished point is 2^32 steps against
a personal budget of 2^112.32 steps, so the reporting loss
is 2^-80.32, measured at toy scale by the census (its
undetected-merge fraction matches the scaled budget-end race) and
subtracted.

Table sizing under H1: expected inserts = t * 2^-32 = 2^96.32
(the insert count is Binomial(t, 2^-32) under H1; the Chernoff
bound is stated under that premise); capacity 2^98 slots
(3.2x the mean), overflow
probability <= exp(-2^95.9); on overflow the attack
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
- seeding 2^16 walkers: [a_0]P by fixed-base comb (32 point adds
  from a precomputed table; the table is part of the 2^105
  memory allowance), [b_0]T by 4-bit windowed double-and-add
  (64 doublings + 16 adds), <= 112 reference operations per
  walker; 2^16 x 112 = 2^22.8, charged as 2^23 reference
  operations;
- randomness: 512 DRBG bits per walker, ~64 instructions per
  walker = 2^22 instructions = 2^12 units, included in the
  seeding charge;
- distinguished-point processing: 2^96.32 inserts x 20
  instructions = 2^100.6 instructions;
- extraction and verification: one modular inverse (~700
  instructions) plus one verification scalar multiple [x]P == T
  (~256 reference adds = 2^8 units), charged as 2^9 units. The
  attack verifies every extraction before declaring success; the
  toy recorded 0 wrong extractions in
  200 + 100 trials.

W = 1.107 * 1.25 * sqrt(r-1) + 2^23 + 2^90.7 + 2^9,
log2(W) = 128.468584, reported exactly as
time_log2 = 16058573/125000.

## Success probability

- Ideal collision probability within the budget:
  1 - exp(-(1.25 * (1 - reporting loss))^2 / 2) = 0.5421666382.
- Every collision is cross-walker at a true point (parity key),
  the extraction is exact algebra, and every extraction is
  verified against T.
- Efficiency of the r-adding walk relative to the ideal random
  walk: charged at a hard floor of **0.75**. Support: Teske's
  analysis of r-adding walks with r = 32 partitions; the toy
  measurement 1.052 of the all-pair
  ideal at two-sigma 0.129; and the
  floor sits ~29% below the measurement.
- Claimed: 0.5421666382 * 0.75 = 0.4066249787 >= 0.39.
  The claim is made at the floor; it does not depend on the toy
  measurement beyond the floor's support.

## Response to the v6 review

1. *The rescaling z = 32k mod r confines the v6 walk (steps
   2^(j+5) - j <= 2^36), making cross-walker collisions a < 2^-75
   event.* The structured constants are replaced by published
   pseudorandom constants; the executed certificate measures every
   rescaling up to 2^12 (narrowest window 0.6329 of the
   modulus), the same certificate detects the v6 confinement
   (5.93e-67), and the counting heuristic excludes confining
   rescalings of width < 2^240 for all M (expected count
   < 2^-235).
2. *The ledger reproduces the score only after assuming its
   prices.* Every one-off is now itemized with its method and
   instruction count: seeding (fixed-base comb 32 adds + windowed
   double-and-add 80 ops per walker, 2^22.8 charged 2^23),
   randomness (2^16 x 64 instructions inside the seeding charge),
   DP processing (20 instructions per insert: hash 8, probe 4,
   compare/write 8), and extraction WITH verification (inverse
   ~700 instructions + 2^8-unit scalar multiple, charged 2^9).
3. *The Chernoff overflow bound lacks justified premises.* The
   insert count is Binomial(t, 2^-32) under H1, the same
   random-function premise that yields the birthday bound; the
   Chernoff bound is stated under that premise.
4. *Span is not a mixing criterion.* Agreed - that is why v8 does
   not argue from span: the confinement analysis is in terms of
   rescaled step windows, measured for all M <= 2^12 and bounded
   by the counting heuristic for all M.
5. *The declared statement said 2^1997276 work and time_log2 =
   15625/127 while the score object said 1997276/15625.* The
   statement generator was corrected; the declared statement now
   matches the score object and proof token for token (this was
   also fixed in the unsubmitted v7).
6. *The toy cannot support the efficiency floor.* The floor is
   0.75 with the toy measuring 1.052;
   the claim is made at the floor, not at the measurement, and the
   floor is supported by Teske's r-adding-walk analysis in
   addition to the toy.

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
