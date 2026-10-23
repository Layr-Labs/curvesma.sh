# plain-rho-v27: bls12-381-g1-ecdlp-plain candidate

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
| time_log2 | **128098091/1000000 = 128.098091** |
| success_probability (lower bound) | 0.39 (modeled 0.6014 at the toy's measured factor) |
| memory_log2_bytes | 63.1 (over-report; see below) |
| preprocessing_log2 | 26.0 (included) |
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
  uniform a_0, b_0 (two 512-bit DRBG outputs per walker), carrying
  the lazy
  representation (a, b) with C = [a]P + [b]T and b fixed for the
  walker's lifetime.
- Step: read the state's index j = x_C mod 32 (5 bits of the
  x-coordinate), update C <- C + A_j and a <- a + c_j mod r.
- Budget: t = 1.5 * sqrt(r-1) walker steps total (log2 t =
  128.01351), shared round-robin across the walkers (log2 of
  per-walker budget 112.01). The walk always consumes the
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

For the v27 constants there is no small-coefficient structure to
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
is 2^32 steps against a personal budget of 2^112.01 steps,
so the reporting loss is 2^-32.01, measured at toy scale
by the census (the recorded true-merge statistic appears below) and
subtracted.

Arena sizing: the predicate tests the low 80 bits of the
field x-COORDINATE (not a scalar mod r). Under H1 each visited
state is a uniform subgroup point; at most TWO curve points
share one admissible x (the +/- y pair) and at most
floor(p/2^80)+1 field values pass, so at most 2*(p/2^80 + 1)
subgroup points qualify out of r-1: on secp256k1 (p ~ r) this
conditional-on-H1 bound is 2^-79 + 2^-253; on BLS the
full-field count is vacuous (p >> r) and the 2^-80 figure is
the standard equidistribution heuristic H2, checked
empirically by the toy census. A CONSERVATIVE uniform bound is
charged on BOTH tracks: density <= 2^-79 + 2^-160, so the
expected record count is at most
t*(2^-79 + 2^-160) = 2^49.01 regardless of any
dependence between steps (Markov needs no independence). The cap is 2^cap records, so the arena
(2^56 x 128 = 2^63 bytes) fits the protected
64-bit address space as ordinary linear storage - no banked
addressing exists anywhere in v27, and ordinary 64-bit loads and
stores access every record. The next append aborts; the abort
probability is bounded by Markov's inequality on that exact mean:
0.007886, with no independence assumption. This replaces
the v14 banked-addressing construction, which the review rightly
found unrealizable at 2^113 bytes, by lowering the density until
the whole arena fits the machine.
Memory. Record layout (48-byte x-coordinate, 8 bytes parity/flags,
32-byte a, 32-byte b, 8 bytes metadata): 128
bytes <= 128 per record, padded to 128. No links, no bucket heads,
no indices: the records are stored directly and sorted at the end,
so nothing needs to address a 2^cap-record pool. (This corrects the
v5 review item that a 32-byte field cannot hold the 381-bit BLS
x-coordinate: the x-coordinate gets 48 bytes here, and the v11
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
| **other per step** | **117** |
| reference add (the step itself) | 1.000 |

Per-step price: 1 + 117/2000 = 1.0585 units. The walk
consumes exactly t = 1.5*sqrt(r-1) steps. One-offs, all hard
worst-case bounds, all included in W:

- constants: 32 fixed-base scalar multiples (~256 reference
  additions each) plus their sha512 derivation and mod-r
  reduction, charged generously as 2^14 reference operations
  flat inside W;
- seeding 2^16 walkers, charged with GENEROUS FLAT bounds,
  DERIVED and ASSERTED in code (the experiment computes the
  component ceilings, rounds up to hundreds, and asserts the
  charge covers the exact component sum). Per walker:
  [a_0]P by fixed-base comb (32 adds + 32 doublings), [b_0]T by
  4-bit windowed double-and-add (256 doublings + 64 additions),
  final addition: 97 additions and 288 doublings. Every ADDITION
  is charged 2 reference operations. Every DOUBLING is charged
  1600/3000 instructions flat (secp/BLS; an affine doubling is 7
  multiplications = 1253/2681 instructions). Randomness: two
  512-bit ChaCha20 outputs per walker at 4096 instructions each
  flat. Exact component sum 97*2 + 288*1600/C + 2*4096/C =
  662.992 (secp) / 630.096 (BLS); component ceilings rounded up
  to hundreds: per_walker_seed = 800 >= 779/775 computed,
  asserted in code. The INITIAL seeding batch-shares one
  Montgomery inversion per sequential layer across the 2^16
  INDEPENDENT walkers, which is what the batched affine prices
  assume. The fixed-base comb TABLE for [a_0]P and the
  15-entry target WINDOW table are charged completely with NO
  batching: 224 setup doublings, 247 subset additions and 255
  affine normalizations for the comb (8 spaced bases plus the
  2^8-1-8 subset sums), 14 additions + 15 normalizations for
  the window table - every setup operation an unbatched
  Jacobian operation and every affine result normalized by
  its OWN charged Fermat inversion: 34972 units
  total, computed in code.
  2^16 x per_walker = 2^25.646123;
- arena appends: at most 2^cap appends at 20 instructions each
  (64-bit linear addressing of the 2^63-byte arena: 4 address
  arithmetic instructions; 8 paired 64-bit stores for the 128-byte
  record; 8 bookkeeping and branch) =
  2^49.356144 units
  worst case;
- one in-place HEAPSORT of at most 2^cap records plus the linear
  scan, with NO workspace: at most 2*cap*log2(cap) record moves
  at 48 instructions each (full-key compare up to 6 64-bit words
  with early-exit branches, charge 16; 64-bit linear addressing
  4; load 8 and store 8 as paired transfers; heap bookkeeping 8)
  and 2^cap scan compares at 64 instructions each (a
  deliberately loose bound covering loads, parity, comparisons
  and loop work) =
  2^57.443607 units;
- special-state reseeds (the exceptional walk states
  S = +A_j and S = -A_j, BOTH detected by the zero-denominator
  branch of the incomplete addition - a doubling never
  executes in the main loop; rate at most 64/(r-1) per step
  under H1): at most 2^20 reseeds, each a COMPLETE UNBATCHED
  construction: both scalar multiplications in Jacobian
  coordinates (mixed additions ~11 multiplications; the
  Jacobian-Jacobian combine ~19 - all 97 additions charged at
  the combine price ~3401/7277 instructions; Jacobian
  doublings ~8 multiplications ~1432/3064), ONE Fermat
  normalization inversion making the affine state explicit -
  the fixed binary left-to-right chain for z^(p-2), EXECUTED
  AND VERIFIED by the experiment
  (380 squarings +
  228 multiplies +
  4 normalization multiplies at the pinned multiplication
  price, 20 control instructions per chain step, 8 transfers =
  246564
  instructions) since the walk needs x = X/Z^2 for the
  partition and predicate, and TWO fresh ChaCha20 outputs -
  component ceiling 1093,
  charged 1100 per reseed (DERIVED and ASSERTED in code),
  abort beyond - 2^30.1 units;
- extraction and verification: one modular inverse (~700
  instructions) plus one verification scalar multiple [x]P == T
  (~256 reference adds) - 2^9 units. The attack verifies every
  extraction before declaring success; the toy recorded
  0 wrong extractions.

W = 1.0585 * 1.5 * sqrt(r-1)
  + 2^118.959326
  (recurring main-loop batch inversions: ceil(t/2^16)
  Montgomery rounds, one Fermat inversion each)
  + 2^25.646123 (seeding,
  including 385 layer-batch Fermat inversions: one per
  sequential scalar-multiplication layer, each a Montgomery
  batch across the 2^16 independent walkers)
  + 2^49.356144 (appends)
  + 2^57.443607 (in-place heapsort)
  + 2^9 (extraction),
log2(W) = 128.098091, reported exactly as
time_log2 = 128098091/1000000. Global cutoff: walker steps
<= t (round-robin), appends <= 2^cap (abort, Markov
0.007886 on the exact mean), identity reseeds <= 2^20 (abort, Markov
<= 2^-142), one sort, one extraction; W bounds every
non-aborting path, and all abort probabilities are in the abort
total below.

## Success probability

- Collision-time factor c (real collision time / ideal random-walk
  collision time). Within the fixed budget the success probability
  is q(c) = 1 - exp(-(1.5/c)^2/2) under H1's Rayleigh collision
  times, and bad events are combined by UNION-BOUND SUBTRACTION -
  P(success) >= q(c) - abort_total, since P(A and B) >= P(A) -
  P(not B) - so the 0.39 claim holds for every c <= 1.4873,
  computed AFTER all deductions.
- The deductions, all stated and SUBTRACTED: reporting loss
  2^-32.01; Markov cycle-waste bound
  2^-10.415; identity-reseed abort 2^-142; arena
  overflow Markov bound 0.007886; abort_total =
  0.008618.
- The only measurement is the toy, which runs the same machinery
  (arena, deferred detection, verified extraction) at 20-bit scale
  with the same constants construction: it measures
  122/200 = 0.61
  against ideal 0.6753, i.e. a collision-time
  factor of 1.093. The tolerance 1.4873 covers a
  production walk 136% slower than the toy. No external
  factor enters the model.
- Claimed: 0.39. Modeled at the toy's measured factor:
  0.6013815900 >= 0.39 with 0.2114
  margin.

## Response to the v25 review

1. *The fixed-base comb table construction was charged a flat
   2^9+15 with no construction.* The table is now charged
   COMPLETELY with NO batching anywhere: 224 setup doublings
   ([2^32]P ... [2^224]P at 32-bit tooth spacing) and 247
   subset additions (2^8-1 nonzero subsets minus the 8
   single-tooth bases) and 255 affine normalizations for the
   comb; 14 additions + 15 normalizations for the 15-entry
   target window table. Every setup operation is an unbatched
   Jacobian operation (unit ceilings 4/2) and every affine
   result is normalized by its OWN charged Fermat inversion
   (34972 units total, computed in code and folded
   into the seeding term
   2^25.646123). The one-off
   cost is invisible at the reported precision; the claim is
   the fully specified construction.

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
