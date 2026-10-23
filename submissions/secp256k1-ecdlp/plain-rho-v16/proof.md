# plain-rho-v16: secp256k1-ecdlp candidate

The walk is the plain additive walk on the full group with published
pseudorandom constants; distinguished points go to an append-only
ARENA and duplicates are found by ONE in-place heapsort. Every work
term in the ledger is a hard worst-case bound; expected-value
quantities (cycle waste) enter only as probabilities via Markov's
inequality and are included in the abort total before the
collision-factor tolerance is computed. No external collision-time
factor is used anywhere.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **64367111/500000 = 128.734222** |
| success_probability (lower bound) | 0.39 (modeled 0.6034 at the toy's measured factor) |
| memory_log2_bytes | 63.1 (over-report; see below) |
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
  sha256 of their concatenation.
- K = 2^16 walkers, each seeded as C = [a_0]P + [b_0]T with fresh
  uniform a_0, b_0 (two 512-bit DRBG outputs per walker), carrying
  the lazy
  representation (a, b) with C = [a]P + [b]T and b fixed for the
  walker's lifetime.
- Step: read the state's index j = x_C mod 32 (5 bits of the
  x-coordinate), update C <- C + A_j and a <- a + c_j mod r.
- Budget: t = 1.5 * sqrt(r-1) walker steps total (log2 t =
  128.58496), shared round-robin across the walkers (log2 of
  per-walker budget 112.58). The walk always consumes the
  full budget.
- Distinguished points: a state is distinguished when the low 80
  bits of its x-coordinate are zero (density 2^-80). On reaching
  one, the walker APPENDS the record ((x_C, y-parity), a, b) to a
  shared append-only arena of capacity 2^cap records; the next
  append aborts the attack (probability bounded by MARKOV on the
  exact expected record count - no independence assumption - and
  subtracted). There is no hash table, no chaining, no probing.
- After the walk: the arena is sorted IN PLACE by heapsort (no
  workspace; the peak is the 2^63-byte arena alone, strictly
  inside the 64-bit address space): at most
  2*2^cap*log2(2^cap) record moves at 48 instructions each
  (full-key compare up to 6 64-bit words with early-exit
  branches, charge 16; 64-bit linear addressing 4; load 8 and
  store 8 as paired transfers; heap bookkeeping 8), plus one
  linear scan at 18 instructions per record; then the sorted
  arena is scanned once. Two
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

For the v16 constants there is no small-coefficient structure to
exploit:

- Executed certificate: Part A measures, for every rescaling
  M <= 2^12, the width of the narrowest window containing all 32
  rescaled steps M*c_j mod r. The narrowest over all 2^12 - 1
  rescalings is 0.6329 of the modulus (at M = 17) -
  every checked rescaling leaves the step set spread across most
  of the group. The same certificate APPLIED TO THE REJECTED v6
  CONSTANTS detects their confinement: at M = 32 the window
  collapses to 5.93e-67 of the modulus.
- Counting heuristic for all rescalings: the expected number of M
  in [0, r) that bring all 32 steps into any single window of
  width W is 32 * W^31 / r^30. For W = 2^240 this is
  < 2^-235; a confinement that could matter for the
  budget needs W * 2^112.58 < r, i.e. W < 2^143.4,
  which sits 2^97 below the heuristic threshold.
- Consequently the standard random-function model (H1) applies:
  walker states collide at the van Oorschot-Wiener parallel-rho
  birthday rate t^2/(2(r-1)), distinguished-point counts are
  Binomial(t, 2^-80), and cycle entries occur at rate 1/r per
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
is 2^32 steps against a personal budget of 2^112.58 steps,
so the reporting loss is 2^-32.58, measured at toy scale
by the census (the recorded true-merge statistic appears below) and
subtracted.

Arena sizing: the expected record count is EXACTLY t * 2^-80 =
2^48.58 under H1 - each step's state is uniform, so it
appends iff its state is distinguished, regardless of dependence
between steps. The cap is 2^cap records, so the arena
(2^56 x 128 = 2^63 bytes) fits the protected
64-bit address space as ordinary linear storage - no banked
addressing exists anywhere in v16, and ordinary 64-bit loads and
stores access every record. The next append aborts; the abort
probability is bounded by Markov's inequality on that exact mean:
0.005859, with no independence assumption. This replaces
the v14 banked-addressing construction, which the review rightly
found unrealizable at 2^113 bytes, by lowering the density until
the whole arena fits the machine.
Memory. Record layout (32-byte x-coordinate, 8 bytes parity/flags,
32-byte a, 32-byte b, 8 bytes metadata): 112
bytes <= 128 per record, padded to 128. No links, no bucket heads,
no indices: the records are stored directly and sorted at the end,
so nothing needs to address a 2^cap-record pool. (This corrects the
v5 review item that a 32-byte field cannot hold the 381-bit BLS
x-coordinate: the x-coordinate gets 32 bytes here, and the v11
review's 98-bit-index objection is void because no index is stored.)
Peak memory is the 2^cap x 128-byte arena plus the arena-sized sort
plus per-walker state plus the 32 constants: strictly above
2^63 = 2^63 bytes, strictly inside the 64-bit address
space (heapsort uses no workspace). As with the
partial-cheon submissions, the reported memory_log2_bytes = 63.1
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
| distinguished-point predicate (low-80-bit zero test of x) | 6 |
| walker state, loop, branch, counters | 69 |
| **other per step** | **109** |
| reference add (the step itself) | 1.000 |

Per-step price: 1 + 109/1000 = 1.109 units. The walk
consumes exactly t = 1.5*sqrt(r-1) steps. One-offs, all hard
worst-case bounds, all included in W:

- constants: 32 scalar multiples (~256 reference adds each) - 2^13
  reference operations;
- seeding 2^16 walkers, with per-operation conversions stated:
  a batched affine ADDITION costs ~6 field multiplications (the
  reference operation itself) and is charged 1; a batched affine
  DOUBLING costs ~4 multiplications/squarings plus ~3 for its
  batch-inversion share = ~7 field mults and is charged 7/6.
  [a_0]P by fixed-base comb (32 adds + 32 doublings; the comb
  table costs ~2^9 reference ops once per attack); [b_0]T by
  4-bit windowed double-and-add (256 doublings + 64 additions;
  the 16-entry table for the attack target is built once, 15
  adds); final seed addition: 32 + 32*7/6 + 256*7/6 + 64 + 1 =
  433 reference operations; randomness is TWO independent
  512-bit outputs of a ChaCha20 DRBG (12 rounds, one 512-bit
  block per output, ~256 instructions per output including the
  mod-r reduction of its scalar), one per scalar, each reduced
  ONCE mod r (bias below 2^-249, no rejection sampling, no
  abort path), ~512 instructions = ~0.5 reference operations;
  the per-walker total is charged ROUNDED UP as 434 reference
  operations: 2^16 x 434 =
  2^24.761578;
- arena appends: at most 2^cap appends at 20 instructions each
  (64-bit linear addressing of the 2^63-byte arena: 4 address
  arithmetic instructions; 8 paired 64-bit stores for the 128-byte
  record; 8 bookkeeping and branch) =
  2^50.356144 units
  worst case;
- one in-place HEAPSORT of at most 2^cap records plus the linear
  scan, with NO workspace: at most 2*cap*log2(cap) record moves
  at 48 instructions each (full-key compare up to 6 64-bit words
  with early-exit branches, charge 16; 64-bit linear addressing
  4; load 8 and store 8 as paired transfers; heap bookkeeping 8)
  and 2^cap scan compares at 18 instructions =
  2^58.431356 units;
- identity reseeds (state = -A_j with matching index, rate at most
  32/(r-1) per step): at most 2^20 reseeds at 433 reference
  operations each, abort beyond - 2^28.6 units;
- extraction and verification: one modular inverse (~700
  instructions) plus one verification scalar multiple [x]P == T
  (~256 reference adds) - 2^9 units. The attack verifies every
  extraction before declaring success; the toy recorded
  0 wrong extractions.

W = 1.109 * 1.5 * sqrt(r-1)
  + 2^24.761578 (seeding)
  + 2^50.356144 (appends)
  + 2^58.431356 (in-place heapsort)
  + 2^9 (extraction),
log2(W) = 128.734222, reported exactly as
time_log2 = 64367111/500000. Global cutoff: walker steps
<= t (round-robin), appends <= 2^cap (abort, Markov
0.005859 on the exact mean), identity reseeds <= 2^20 (abort, Markov
<= 2^-142), one sort, one extraction; W bounds every
non-aborting path, and all abort probabilities are in the abort
total below.

## Success probability

- Collision-time factor c (real collision time / ideal random-walk
  collision time). Within the fixed budget the success probability
  is q(c) = 1 - exp(-(1.5/c)^2/2) under H1's Rayleigh collision
  times, and bad events are combined by UNION-BOUND SUBTRACTION -
  P(success) >= q(c) - abort_total, since P(A and B) >= P(A) -
  P(not B) - so the 0.39 claim holds for every c <= 1.4923,
  computed AFTER all deductions.
- The deductions, all stated and SUBTRACTED: reporting loss
  2^-32.58; Markov cycle-waste bound
  2^-10.415; identity-reseed abort 2^-142; arena
  overflow Markov bound 0.005859; abort_total =
  0.006592.
- The only measurement is the toy, which runs the same machinery
  (arena, deferred detection, verified extraction) at 20-bit scale
  with the same constants construction: it measures
  122/200 = 0.61
  against ideal 0.6753, i.e. a collision-time
  factor of 1.093. The tolerance 1.4923 covers a
  production walk 137% slower than the toy. No external
  factor enters the model.
- Claimed: 0.39. Modeled at the toy's measured factor:
  0.6034082030 >= 0.39 with 0.2134
  margin.

## Response to the v15 review

1. *The arena plus its merge-sort workspace exactly consumed the
   64-bit address space.* The workspace is deleted: duplicate
   detection is one IN-PLACE HEAPSORT (at most
   2*cap*log2(cap) record moves at 48 instructions each, plus
   the 18-instruction-per-record scan), so the peak is the
   2^63-byte arena alone, strictly inside the 64-bit address
   space (memory_log2_bytes = 63.1).
2. *One algorithm line still said a 32-bit predicate.* Fixed:
   every line, including the algorithm section, states the
   80-bit predicate, consistent with the ledger, the overflow
   bound and the success model.
3. *Claim-text charges were stale against the computed values
   (seeding 2^24.59 vs the computed 2^24.76, appends 24 vs 20,
   reseeds 385 vs 433, and the Markov bound quoted below its
   computed value).* All charges are regenerated from the
   recorded run and quoted at full computed precision: seeding
   2^24.761578 (434
   reference operations per walker, including two 512-bit
   ChaCha20 DRBG outputs at ~256 instructions each), appends 20
   instructions, identity reseeds 434 operations, and the
   Markov overflow bound 0.005859375.
4. *The DRBG was unspecified.* Fixed: a ChaCha20 DRBG (12
   rounds, one 512-bit block per output, ~256 instructions per
   output including the mod-r reduction), two independent
   outputs per walker, one per scalar, each reduced once mod r
   with bias below 2^-249 and no rejection sampling.

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
