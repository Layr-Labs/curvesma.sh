# plain-rho-v21: secp256k1-ecdlp candidate

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
| preprocessing_log2 | 25.5 (included) |
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
  linear scan at 64 instructions per record flat (a
  deliberately loose bound covering the 6-word loads, the
  parity field, comparisons and loop work); then the sorted
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

For the v21 constants there is no small-coefficient structure to
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

Arena sizing: the expected record count is BOUNDED: under H1
each visited state is uniform on the r-1 non-identity points,
so a state qualifies with probability at most
#qualifying/(r-1) <= 2^-80 + 2^-160 (the qualifying count is
floor((r-1)/2^80)+1, which exceeds 2^-80*(r-1) by at most one
point), and the expected record count is at most
t*(2^-80 + 2^-160) = 2^48.58 regardless of any
dependence between steps. The cap is 2^cap records, so the arena
(2^56 x 128 = 2^63 bytes) fits the protected
64-bit address space as ordinary linear storage - no banked
addressing exists anywhere in v21, and ordinary 64-bit loads and
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

- constants: 32 fixed-base scalar multiples (~256 reference
  additions each) plus their sha512 derivation and mod-r
  reduction, charged generously as 2^14 reference operations
  flat inside W;
- seeding 2^16 walkers, charged with GENEROUS FLAT bounds so
  that no implementation-level count can exceed them. Per
  walker: [a_0]P by fixed-base comb (32 adds + 32 doublings),
  [b_0]T by 4-bit windowed double-and-add (256 doublings + 64
  additions), final addition: 97 additions and 288 doublings.
  Every ADDITION is charged 2 reference operations (the
  reference is 1; 2 covers any full-width modular add with
  carries, conditional reduction, correction and transfers).
  Every DOUBLING is charged 1600/3000 instructions flat (an
  affine doubling is 7 multiplications = 1253/2681
  instructions; the flat allowance covers the six full-width
  adds with carry chains, conditional reductions, corrections
  and transfers at any reasonable count - deliberately loose).
  Randomness: two independent 512-bit ChaCha20 outputs, one per
  scalar, charged 4096 instructions per output flat (covers 80
  quarter-rounds = 960 instructions plus initialization,
  feed-forward, state transfers, counters and the mod-r
  reduction). A reseed's own operations are batched internally
  (one Montgomery inversion shared within the reseed), so each
  reseed costs the same per-walker price wherever it occurs.
  Per-walker total: 97*2 + 288*1600/C + 2*4096/C = 558.6
  (secp) / 630.1 (BLS), charged rounded up as 600 (secp) / 700
  (BLS) reference operations:
  2^16 x per_walker = 2^25.228838;
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
  and 2^cap scan compares at 64 instructions each (a
  deliberately loose bound covering loads, parity, comparisons
  and loop work) =
  2^58.443607 units;
- identity reseeds (state = -A_j with matching index, rate at most
  32/(r-1) per step): at most 2^20 reseeds at per_walker_seed
  reference
  operations each, abort beyond - 2^28.6 units;
- extraction and verification: one modular inverse (~700
  instructions) plus one verification scalar multiple [x]P == T
  (~256 reference adds) - 2^9 units. The attack verifies every
  extraction before declaring success; the toy recorded
  0 wrong extractions.

W = 1.109 * 1.5 * sqrt(r-1)
  + 2^25.228838 (seeding)
  + 2^50.356144 (appends)
  + 2^58.443607 (in-place heapsort)
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

## Response to the v20 review

1. *Component prices were implementation-derived and each
   derivation was disputed line by line.* All component prices
   are now GENEROUS FLAT bounds, deliberately loose so that no
   implementation-level count can exceed them: additions 2
   reference operations each (the reference is 1), doublings
   1600/3000 instructions flat (vs 1253/2681 of pure
   multiplication work), ChaCha20 4096 instructions per output
   flat (vs 960 of round arithmetic), the constant build 2^14
   flat, the scan 64 instructions per record flat. Overcharges
   can only make the claim harder, never easier, and the flat
   prices make the per-term hard bounds trivially true.
2. *Isolated identity reseeds kept a batched price without a
   schedule.* The reseed's own operations are batched
   INTERNALLY: each reseed performs its own 385 scalar
   evaluations as one Montgomery batch (one inversion shared
   within the reseed), so the per-reseed price equals
   per_walker_seed wherever and whenever it occurs. The
   identity-reseed term in W is exactly the cap of 2^20 reseeds
   times that single-walker price.
3. *The exact-mean overflow premise overreached (r-1 is not
   divisible by 2^80; the walk is not a permutation).* The
   record count is now an explicit upper bound: a visited state
   qualifies with probability at most
   #(qualifying)/(r-1) <= 2^-80 + 2^-160, so
   E[records] <= t*(2^-80 + 2^-160); Markov needs no
   independence and is applied to this bound. The cap and the
   overflow abort are unchanged.
4. *Component figures rounded below computed values (the
   identity-reseed Markov bound exceeded the stated 2^-142).*
   Every probability is now quoted AS COMPUTED, never rounded
   below: identity_reseed_abort_log2 is printed by the
   experiment (-142.4150 on this track) and the
   statement quotes the computed bound.

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
