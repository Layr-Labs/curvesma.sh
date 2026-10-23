# plain-rho-v9: secp256k1-ecdlp candidate

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
| preprocessing_log2 | 24.5 (included) |
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

For the v9 constants there is no small-coefficient structure to
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
- Longer cycles under H1: a fruitless cycle is a SELF-repeat of
  one walker - a pair (u, v), u < v, on one trail whose window
  sum is 0 mod r, an event of probability 1/r per self pair. The
  number of self pairs is K*L*(L-1)/2 with K = 2^16 walkers and
  L = t/K, so the expected number of cycle entries is
  t^2/(2Kr) = 1.192e-05 - the self-collision birthday
  scale, not a per-step rate. A cycle containing a distinguished
  point is detected at its first DP revisit (cost at most one
  lap) and the walker reseeds; a cycle containing none traps the
  walker for at most its remaining budget. Charging the
  remaining budget for every entry at (u, v) bounds the total
  expected waste by K * sum_(u<v) (L - v)/r = K*L^3/(6r) =
  t^3/(6 K^2 r) steps, i.e. a fraction t^2/(6 K^2 r) =
  6.063e-11 of the budget. The bound is computed over the
  full un-terminated budget, so it also bounds the terminated
  process, and the fraction is subtracted in the success model
  alongside the reporting loss.
- The toy observes 3 cycle reseeds over
  200 trials against a self-pair model expectation
  of 19.2 for the same run - below the model,
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
a personal budget of 2^112.32 steps, so the reporting loss
is 2^-80.32. The census records this statistic
directly: 56 of its trials contain
a true first-contact merge (counted by a diagnostic full-visit
log that is not part of the attack), 28
are detected and 28 end
undetected by the budget - the budget-end race scaled to toy
size, where the gap to the next distinguished point (2^6 = 64
steps) is comparable to the per-walker budget
(159 steps). The same law at production
scale is the 2^-80.32 term above, subtracted.

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
  from a precomputed table; building the comb table costs ~2^9
  reference ops once per attack), [b_0]T by 4-bit windowed
  double-and-add (256 doublings + 64 window additions; the
  16-entry window table for the attack target is built once,
  15 adds); 32 + 256 + 64 = 352 reference operations per
  walker; 2^16 x 352 = 2^24.46, charged as 24.46
  reference operations (the table-build one-offs included);
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

W = 1.107 * 1.25 * sqrt(r-1) + 2^24.46 + 2^90.7 + 2^9,
log2(W) = 128.468584, reported exactly as
time_log2 = 16058573/125000.

## Success probability

- Ideal collision probability within the budget:
  1 - exp(-(1.25 * (1 - reporting loss))^2 / 2) = 0.5421666382.
- Collisions are overwhelmingly cross-walker (self-repeats are
  the cycle entries bounded above), every key hit at a true point
  (parity key) is a collision, the extraction is exact algebra,
  and every extraction is verified against T. The reporting loss
  2^-80.32 and the cycle waste fraction
  6.063e-11 are both subtracted in the ideal below.
- Efficiency of the r-adding walk relative to the ideal random
  walk: charged at a hard floor of **0.75**. Support: Teske's
  analysis of r-adding walks with r = 32 partitions; the toy
  measurement 1.052 of the all-pair
  ideal at two-sigma 0.129; and the
  floor sits ~29% below the measurement.
- Claimed: 0.5421666382 * 0.75 = 0.4066249786 >= 0.39.
  The claim is made at the floor; it does not depend on the toy
  measurement beyond the floor's support.

## Response to the v8 review

1. *The cycle-entry rate t/r < 2^-127 does not follow from the
   random-function model; self-repeat opportunities scale as
   K*L*(L-1)/(2r) = t^2/(2Kr).* Corrected throughout: cycle
   entries are counted per self pair, giving t^2/(2Kr) =
   1.192e-05 expected entries, and the waste bound
   K*sum_(u<v)(L-v)/r = t^3/(6K^2 r) steps (a fraction
   6.063e-11 of the budget, computed over the
   un-terminated budget so it bounds the terminated process
   too) is subtracted in the success model. Detection is split
   honestly: cycles containing a distinguished point are caught
   at the first DP revisit; cycles containing none waste at most
   the remaining budget; both are inside the same bound.
2. *The 4-bit windowed seeding multiplication is 256 doublings
   + 64 additions, not 64 + 16; the 2^23 preprocessing
   allowance is not justified.* Corrected: per walker 32 (comb)
   + 256 + 64 = 352 reference operations, 2^16 x 352 =
   2^24.46 charged as 24.46 (preprocessing_log2 =
   24.5), including the once-per-attack table builds
   (~2^9 for the comb table, 15 adds for the 16-entry T-window
   table) and the 512 DRBG bits per walker (2^22 instructions
   total). The score is unchanged at six decimals.
3. *The census does not record the undetected-merge statistic.*
   It now does: the recorded census output contains
   56 trials with a true
   first-contact merge, 28 detected,
   28 undetected by budget end
   (a diagnostic full-visit log counts first contacts; it is
   not part of the attack), matching the scaled budget-end race.
4. *The Teske reference supplies no applicable quantitative
   result.* Stated quantitatively: Teske's experiments on
   r-adding walks with r >= 16 measure expected collision times
   at most ~1.1x the random-walk ideal; the floor 0.75 covers a
   production walk 33% slower than that worst measured factor
   and sits 29% below this toy's measured ratio of
   1.052.
3. *The Chernoff overflow bound lacks justified premises.* The
   insert count is Binomial(t, 2^-32) under H1, the same
   random-function premise that yields the birthday bound; the
   Chernoff bound is stated under that premise.
4. *Span is not a mixing criterion.* Agreed - that is why v9 does
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
