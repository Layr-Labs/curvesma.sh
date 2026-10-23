# plain-rho-v13: bls12-381-g1-ecdlp-plain candidate

The walk is the plain additive walk on the full group with published
pseudorandom constants; distinguished points go to an append-only
ARENA and duplicates are found by ONE final merge sort. Every work
term in the ledger is a hard worst-case bound; expected-value
quantities (cycle waste) enter only as probabilities via Markov's
inequality and are included in the abort total before the
collision-factor tolerance is computed. No external collision-time
factor is used anywhere.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **128088699/1000000 = 128.088699** |
| success_probability (lower bound) | 0.39 (modeled 0.6090 at the toy's measured factor) |
| memory_log2_bytes | 114.1 (over-report; see below) |
| preprocessing_log2 | 25.0 (included) |
| nonuniform_advice_log2_bytes | 0 |

## Algorithm

- Group: BLS12-381 G1: y^2 = x^3 + 4 over the 381-bit prime field, r the 255-bit prime order, standard generator. Target T = [x]P for uniform secret
  x in [0, r); recover x.
- Constants: A_j = [c_j]P for j = 0..31, where c_j is the published
  deterministic pseudorandom value
  c_j = int(sha512("curvesmash/bls12-381-g1-ecdlp-plain/constant/j")) mod r.
  The 512-bit hash value reduced mod r has bias below 2^-500.
  Part A verifies the 32 constants are distinct and prints the
  sha256 of their concatenation.
- K = 2^16 walkers, each seeded as C = [a_0]P + [b_0]T with fresh
  uniform a_0, b_0 (512 DRBG bits per walker), carrying the lazy
  representation (a, b) with C = [a]P + [b]T and b fixed for the
  walker's lifetime.
- Step: read the state's index j = x_C mod 32 (5 bits of the
  x-coordinate), update C <- C + A_j and a <- a + c_j mod r.
- Budget: t = 1.5 * sqrt(r-1) walker steps total (log2 t =
  128.01351), shared round-robin across the walkers (log2 of
  per-walker budget 112.01). The walk always consumes the
  full budget.
- Distinguished points: a state is distinguished when the low 32
  bits of its x-coordinate are zero (density 2^-32). On reaching
  one, the walker APPENDS the record ((x_C, y-parity), a, b) to a
  shared append-only arena of capacity 2^cap records; the next
  append aborts the attack (probability bounded by MARKOV on the
  exact expected record count - no independence assumption - and
  subtracted). There is no hash table, no chaining, no probing.
- After the walk: the arena is merge-sorted by key (a buffered
  merge sort: at most 2^cap * log2(2^cap) record moves at 30
  instructions each - key compare 3, two-word address 6, load 8
  and store 8 as four paired transfers, bookkeeping 5 - plus one
  linear scan; the sort uses one arena-sized workspace buffer,
  which is declared in the memory figure) and scanned once. Two
  adjacent records
  with equal keys and DIFFERENT representations are a collision:
  x = (a2 - a1)/(b1 - b2) mod r is extracted and verified against
  T. This is the same detection event as immediate table lookup:
  a collision is detected iff some distinguished point was
  recorded by two walkers.
- The key is the pair (x-coordinate, y-parity). On the full group
  +X and -X are DISTINCT points that share an x-coordinate; the
  parity bit is what makes equal keys true point collisions, so
  the extraction is exact algebra with no false positives (the
  toy measured 0 wrong extractions in
  200 trials plus 100 census trials).

## Confinement: why the constants are pseudorandom

The lazy representation makes the walker's scalar k = a + b*x mod r
walk with steps c_j; two walkers collide at trail positions (u, v)
exactly when their scalar-difference walk hits the (uniform) seed
difference. For this to happen at the birthday rate, no invertible
rescaling of the scalar may confine the step set. The v6 constants
failed exactly this: r = 32Q + rho makes 32*c_j = 32*(jQ + 2^j) =
2^(j+5) - j (mod r), so z = 32k walks inside an interval of width
2^36 * L and cross-walker collisions are a < 2^-75 event at
production scale.

For the v13 constants there is no small-coefficient structure to
exploit:

- Executed certificate: Part A measures, for every rescaling
  M <= 2^12, the width of the narrowest window containing all 32
  rescaled steps M*c_j mod r. The narrowest over all 2^12 - 1
  rescalings is 0.6897 of the modulus (at M = 1266) -
  every checked rescaling leaves the step set spread across most
  of the group. The same certificate APPLIED TO THE REJECTED v6
  CONSTANTS detects their confinement: at M = 32 the window
  collapses to 1.31e-66 of the modulus.
- Counting heuristic for all rescalings: the expected number of M
  in [0, r) that bring all 32 steps into any single window of
  width W is 32 * W^31 / r^30. For W = 2^240 this is
  < 2^-201; a confinement that could matter for the
  budget needs W * 2^112.01 < r, i.e. W < 2^142.8,
  which sits 2^97 below the heuristic threshold.
- Consequently the standard random-function model (H1) applies:
  walker states collide at the van Oorschot-Wiener parallel-rho
  birthday rate t^2/(2(r-1)), distinguished-point counts are
  Binomial(t, 2^-32), and cycle entries occur at rate 1/r per
  self pair. Part B measures the realized fraction of the ideal
  directly: 122/200 = 0.61
  against the same-parameter ideal 0.6753
  (ratio 0.903, two-sigma
  0.102), with the full budget always
  consumed, exactly as the ledger charges it.

## Fruitless cycles

- No cycle word of length <= 4 exists: Part A enumerates every
  multiset of at most 4 constants and asserts none sums to 0 mod r.
  A cycle of length L requires its realized index multiset to sum
  to 0 mod r, so this excludes every cycle of length <= 4
  regardless of index realizability.
- Longer cycles under H1: a fruitless cycle is a SELF-repeat of one
  walker - a pair (u, v), u < v, on one trail whose window sum is
  0 mod r, an event of probability 1/r per self pair. The number
  of self pairs is K*L*(L-1)/2 with K = 2^16, L = t/K, so the
  expected number of cycle entries is t^2/(2Kr) = 1.717e-05.
  A trapped walker wastes at most its remaining budget; charging
  the remaining budget L - v for every entry at (u, v) bounds the
  total expected waste by K * sum_(u<v) (L - v)/r = K*L^3/(6r) =
  t^3/(6K^2 r) steps, a fraction 8.731e-11 of the budget.
  The bound is computed over the full budget, so it also bounds
  the terminated process. The expectation is NOT substituted as a
  work or budget deduction: by Markov's inequality,
  P(waste exceeds 2^-23 of the budget) <= 7.32e-04 =
  2^-10.415, and this probability is included in the
  abort total below.
- The toy observes 11 sort-time
  identical-representation duplicates (cycle evidence) over
  200 trials against a self-pair model expectation of
  27.7 - below the model,
  i.e. the discrepancy is in the conservative direction.

## Distinguished points, reporting, and memory

A collision is detected after the walk iff some distinguished point
was recorded by two walkers. If the collision point itself is not
distinguished, the later walker has merged onto the earlier walker's
trail; because the step map is deterministic, the merged pair CANNOT
separate, and it reaches the next shared distinguished point
together - so both records exist unless the walker's personal
budget ends first. The expected gap to the next distinguished point
is 2^32 steps against a personal budget of 2^112.01 steps,
so the reporting loss is 2^-80.01, measured at toy scale
by the census (the recorded true-merge statistic appears below) and
subtracted.

Arena sizing: the expected record count is EXACTLY t * 2^-32 =
2^96.01 under H1 - each step's state is uniform, so it
appends iff its state is distinguished, regardless of dependence
between steps. The cap is 2^cap records; the next append aborts;
the abort probability is bounded by Markov's inequality on that
exact mean: 0.000986, with no independence assumption.
This replaces the v12 Binomial/Chernoff premise, which the
correlated walk (merged trails repeat distinguished statuses) does
not justify.

Memory. Record layout (48-byte x-coordinate, 8 bytes parity/flags,
32-byte a, 32-byte b, 8 bytes metadata): 128
bytes <= 128 per record, padded to 128. No links, no bucket heads,
no indices: the records are stored directly and sorted at the end,
so nothing needs to address a 2^cap-record pool. (This corrects the
v5 review item that a 32-byte field cannot hold the 381-bit BLS
x-coordinate: the x-coordinate gets 48 bytes here, and the v11
review's 98-bit-index objection is void because no index is stored.)
Peak memory is the 2^cap x 128-byte arena plus the arena-sized sort
workspace buffer plus per-walker state plus the 32 constants:
strictly above 2^114 bytes. As with the
partial-cheon submissions, the reported memory_log2_bytes = 114.1
is the rounded-up over-report of that peak; memory is a
supplementary metric and the over-report direction is conservative.

## Work ledger: every term a hard bound

Per walker step (AArch64 instruction counts, the "other" category of
classical-group-work-v2; one reference operation = batched incomplete
affine addition is charged as 1.000 and every walk step performs
exactly one):

| Item | Instructions |
| --- | --- |
| index j = x_C mod 32 (5-bit extract) | 2 |
| representation add a += c_j mod r (full-width constant add + reduce) | 32 |
| distinguished-point predicate (low-32-bit OR-reduce of x) | 4 |
| walker state, loop, branch, counters | 69 |
| **other per step** | **107** |
| reference add (the step itself) | 1.000 |

Per-step price: 1 + 107/2000 = 1.0535 units. The walk
consumes exactly t = 1.5*sqrt(r-1) steps. One-offs, all hard
worst-case bounds, all included in W:

- constants: 32 scalar multiples (~256 reference adds each) - 2^13
  reference operations;
- seeding 2^16 walkers: [a_0]P by fixed-base comb from a
  precomputed table (32 point adds AND 32 doublings per
  evaluation; building the comb table costs ~2^9 reference ops
  once), [b_0]T by 4-bit windowed double-and-add (256 doublings +
  64 window additions; the 16-entry table for the attack target is
  built once, 15 adds), plus the final seed addition:
  32 + 32 + 256 + 64 + 1 = 385 reference operations per walker,
  2^16 x 385 = 2^24.59 charged as
  24.59, including the 512
  DRBG bits per walker (2^22 instructions);
- arena appends: at most 2^cap appends at 20 instructions each
  (8 paired 64-bit stores for the 128-byte record, 4 two-word
  address instructions, 8 bookkeeping) = 2^110.3
  instructions worst case = 2^110.3 -
  2^99.4 units;
- one final merge sort of at most 2^cap records:
  cap*log2(cap) record moves at 30 instructions each (key compare 3,
  two-word address 6, load 8 and store 8, bookkeeping 5) =
  2^106.6690 units, plus one
  arena-sized workspace buffer declared in the memory figure;
- identity reseeds (state = -A_j with matching index, rate at most
  32/(r-1) per step): at most 2^20 reseeds at 385 reference
  operations each, abort beyond - 2^28.6 units;
- extraction and verification: one modular inverse (~700
  instructions) plus one verification scalar multiple [x]P == T
  (~256 reference adds) - 2^9 units. The attack verifies every
  extraction before declaring success; the toy recorded
  0 wrong extractions.

W = 1.0535 * 1.5 * sqrt(r-1) + 2^24.59 (seeding)
  + 2^99.3561 (appends)
  + 2^106.6690 (sort)
  + 2^9 (extraction),
log2(W) = 128.088699, reported exactly as
time_log2 = 128088699/1000000. Global cutoff: walker steps
<= t (round-robin), appends <= 2^cap (abort, Markov
0.000986 on the exact mean), identity reseeds <= 2^20 (abort, Markov
<= 2^-142), one sort, one extraction; W bounds every
non-aborting path, and all abort probabilities are in the abort
total below.

## Success probability

- Collision-time factor c (real collision time / ideal random-walk
  collision time). Within the fixed budget the success probability
  is q(c) = 1 - exp(-(1.5/c)^2/2) under H1's Rayleigh collision
  times, and bad events are combined by UNION-BOUND SUBTRACTION -
  P(success) >= q(c) - abort_total, since P(A and B) >= P(A) -
  P(not B) - so the 0.39 claim holds for every c <= 1.5043,
  computed AFTER all deductions.
- The deductions, all stated and SUBTRACTED: reporting loss
  2^-80.01; Markov cycle-waste bound
  2^-10.415; identity-reseed abort 2^-142; arena
  overflow Markov bound 0.000986; abort_total =
  0.001718.
- The only measurement is the toy, which runs the same machinery
  (arena, deferred detection, verified extraction) at 20-bit scale
  with the same constants construction: it measures
  122/200 = 0.61
  against ideal 0.6753, i.e. a collision-time
  factor of 1.093. The tolerance 1.5043 covers a
  production walk 138% slower than the toy. No external
  factor enters the model.
- Claimed: 0.39. Modeled at the toy's measured factor:
  0.6089519161 >= 0.39 with 0.2190
  margin.

## Response to the v12 review

1. *The sort charge rounded the recorded figure down (98.8 declared
   vs 98.8189 computed) and lacked instruction itemization and
   workspace accounting.* The sort is now itemized at 30
   instructions per record move (key compare 3, two-word address 6,
   load 8 and store 8, bookkeeping 5), the recorded computed value
   is quoted at four decimals
   (2^106.6690 units), and the
   arena-sized workspace buffer is declared in the memory figure
   (peak strictly above 2^114 bytes, reported 114.1).
2. *Distinguished-point counts are not Binomial under H1 once
   merged trails repeat statuses.* Agreed - the Binomial/Chernoff
   premise is dropped entirely. The expected record count is
   exactly t * 2^-32 under H1 (marginal uniformity per step, no
   independence), and the overflow probability is bounded by
   Markov's inequality on that exact mean: 0.000986 at the
   2^cap-record cap. No correlation assumption is made anywhere.
3. *Marginal Markov bounds justify subtracting, not multiplying.*
   Corrected: the success bound is now P(success) >= q(c) -
   abort_total with all abort terms subtracted, and the
   collision-factor tolerance 1.5043 is computed from the
   subtractive inequality.

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
self-test), the hard-bound ledger arithmetic, and the success
bound >= 0.39 after all deductions.
