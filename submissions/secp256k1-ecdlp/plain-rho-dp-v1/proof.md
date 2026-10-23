# Plain parallel Pollard rho with distinguished points, conservatively charged

This candidate claims **time_log2 = 12903/100 = 129.03** target-group-operation
equivalents (C = 1000 scalar AArch64 instruction equivalents per reference
operation), with **success_probability = 0.40** over a uniform challenge
x in [0, r) and fresh algorithm coins, **memory_log2_bytes = 63.01**,
**preprocessing_log2 = 25** (included in the total), and
**nonuniform_advice_log2_bytes = 0**. It is a classical, single-target, cold
attack. The success bound rests on the standard, disclosed random-walk
heuristic for r-adding walks (section "Heuristic"); every other step is either
exact or bounded with an explicit failure term.

The design is deliberately the textbook construction of van Oorschot and
Wiener (J. Cryptology 1999) with Teske's r-adding walk (Math. Comp. 2001),
**without** the negation map or the GLV endomorphism. Consequently there are no
fruitless cycles to detect, escape or charge, no class canonicalisation and no
orbit bookkeeping. Nothing from earlier CurveSmash submissions is relied on.

## Target

The protected secp256k1 profile is used unchanged:

```
p = fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
r = fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
```

P is the fixed standard generator; Q = [x]P with x uniform in [0, r). No
reduced challenge, interval promise or advice is assumed.

## Parameters

| Symbol | Value | Meaning |
| --- | --- | --- |
| R | 2^10 = 1024 | number of walk increments (partitions) |
| theta | 76 | distinguished-point (DP) bits |
| W | 2^32 | number of concurrent walks |
| m | 2^12 | walks per processor sharing one batched inversion (2^20 processors) |
| cap | 2^82 = 64 * 2^theta | maximal segment length |
| k1 | ceil(2^128.019214) | phase-1 step budget, see "Success" |
| L1 | ceil(k1 / W) | phase-1 steps per walk |

All work of all processors is summed; parallelism only changes wall-clock time.

## Algorithm

1. **Input validation.** Check that Q is a canonical affine point on
   y^2 = x^3 + 7 and that [r]Q = O. If Q = O output 0 (verified).
2. **Increment table (preprocessing).** Draw c_j, d_j uniformly in [0, r) for
   j < R (rejection sampling of 256-bit strings, at most 64 draws per scalar,
   abort on exhaustion; probability < 2^-8000) and compute the affine points
   M_j = [c_j]P + [d_j]Q with complete Jacobian arithmetic and one inversion each.
   If some M_j = O, redraw it.
3. **Segment start.** A walk starts a segment by drawing a0, b0 uniformly in
   [0, r), computing X = [a0]P + [b0]Q (complete arithmetic, affine output), and
   storing (a0, b0) and a length counter ell = 0 in its state. If X = O, redraw.
   Because a0 is uniform and independent of b0, X is uniform in the group.
   The DP predicate (step 5) is applied to the start point as well: if X is
   distinguished, the record (x, a0, b0, 0) is stored and a new start drawn.
4. **Walk step.** With X = (x, y) in canonical affine form, let
   j = x mod 2^10 (the 10 low bits of x). If x = x(M_j) the step is exceptional:
   the walk discards its segment and starts a new one (step 3). Otherwise
   X <- X + M_j by the incomplete affine formula
   lambda = (y_j - y)/(x_j - x), x' = lambda^2 - x - x_j,
   y' = lambda (x - x') - y, and ell <- ell + 1. The m walks of one processor
   advance together and share one field inversion by Montgomery's trick
   (prefix products, one inversion, backward pass). An exceptional walk
   contributes the factor 1 to the batch so that no zero is ever inverted.
   The walk function X -> X + M_{h(X)} is the same deterministic map for all
   walks; no coefficients (a, b) are updated during the walk.
5. **Distinguished points.** After a step, X is distinguished iff bits
   10..85 of the canonical x are all zero (76 bits, disjoint from the
   partition bits). The walk then appends the 128-byte record
   (x, a0, b0, ell) to its processor's DP log and starts a new segment
   (step 3). If ell reaches cap without a DP, the segment is abandoned and a
   new one started. If the total number of segment starts exceeds 2^56 the
   attack aborts with failure.
6. **Schedule.** Phase 1: every walk performs exactly L1 steps. Phase 2: every
   walk performs exactly cap further steps under the same rules. Then all DP
   logs are concatenated (abort with failure if more than 2^55 records exist).
7. **Collision search.** Merge-sort all records by x and scan neighbours for
   equal x. For at most 16 equal-x pairs, in sorted order: re-walk both segments
   from their stored (a0, b0) for their stored lengths, now tracking
   coefficients (a, b) <- (a + c_j, b + d_j) mod r with complete arithmetic, to
   obtain a1 P + b1 Q = D1 and a2 P + b2 Q = D2 with x(D1) = x(D2). If D1 = -D2
   negate (a2, b2). If b1 != b2 mod r, set
   x* = (a2 - a1)/(b1 - b2) mod r, verify [x*]P = Q and output x* on success.
   Otherwise continue with the next pair; after 16 pairs output failure.

The relation step is exact: D1 = +-D2 and the invariant D = [a]P + [b]Q imply
a1 + b1 x = +-(a2 + b2 x) mod r, so a nonzero denominator yields the correct x,
and the charged final verification makes a wrong output impossible.

## Success

**Detection is exact up to explicit tails.** Suppose that among the points
visited during phase 1 (segment starts and step results), two visits of
different segments hit the same group element Z. Both segments then follow the
identical deterministic path from Z. The earlier-visiting segment ends at the
first DP D on that path unless it is abandoned at the cap; the later one
reaches the same D, because phase 2 gives every walk cap more steps and its
segment length is also capped at cap. Hence two records with x(D) appear,
unless (i) a segment on that path exceeds cap steps without a DP, or (ii) an
exceptional step occurs on the shared path. A segment of a random-like walk
exceeds 64 * 2^76 steps without a DP with probability (1 - 2^-76)^(2^82) < e^-64
< 2^-92. An exceptional step happens with probability at most 2R/r per step
(x = x(M_j) for one of 2 points per partition), so over all at most 2^130
executed steps the probability of any exceptional step is < 2^-115. If both
segments have distinct uniform starts (a0, b0), the resulting relation has
b1 = b2 mod r only with probability about 1/r; we bound the union of
degenerate relations by 2^-100. A segment colliding with itself (a cycle
without DP; its share of repeats is about 2 * 2^theta / k1) requires a self-collision within at most 2^82 steps, probability
< 2^164 / 2^256 = 2^-92 per segment, < 2^-36 over at most 2^56 segments.
The start-count abort (2^56, versus at most about 2^52.1 DP-ended, 2^46.1
cap-abandoned and 2^32 initial segments) is reached with probability
< 2^-1000. The DP record
count exceeds 2^55 (eight times its mean of about 2^52.02) with probability
< 2^-1000 (Chernoff). The 16-pair limit is never binding when the first pair
is non-degenerate. Total explicit failure terms: delta < 2^-35.

**Collision probability (heuristic H).** Let N = r and N_eff = N / (1 - 2/R).
Under H (below), the probability that k visited points are pairwise distinct
is at most exp(-k(k-1) / (2 N_eff)). Phase 1 visits at least k1 points, where

```
k1 = ceil( sqrt( 2 * N_eff * ln(1/0.599) ) ) = 2^128.019214 (rounded up),
```

so a collision occurs in phase 1 with probability at least 1 - 0.599 = 0.401.
Subtracting delta < 2^-35 gives **success >= 0.40**, above the required 0.39
(the 0.39 threshold alone would need about 1.8 % fewer points). Only phase-1
points are credited; phase-2 points can only help. H is applied to the whole
distribution of the first-repeat time, not only to its mean; the toy test
below measures exactly this event (a repeat among phase-1 visited points) at
quantiles around 0.22, 0.40 and 0.68.

## Heuristic H and its support

H: for the r-adding walk with R = 1024 uniformly random increments M_j,
partition index from the low 10 bits of x and DP predicate from the next 76
bits, run in parallel with uniform restarts, the number of visited points
before the first repeated group element is stochastically dominated by the
first-repeat time of i.i.d. uniform sampling from a set of N_eff = N / (1 - 2/R)
elements, i.e. P(no repeat among k points) <= exp(-k(k-1) / (2 N_eff)).

Support:

- Teske (2001) showed experimentally that r-adding walks with r >= 16 behave
  close to random mappings; Brent and Pollard's analysis and Bernstein and
  Lange ("Computing small discrete logarithms faster", 2012, appendix) give
  the inflation factor 1/sqrt(1 - 1/R) for the expected collision time, with
  higher-order corrections of order 1/R^2. We charge 1/(1 - 2/R) on N, which
  dominates 1/(1 - 1/R - 1/R^2 - ...) = 1/(1 - 1/(R-1)) for R = 1024.
- Parallel collision search with distinguished points and fresh uniform
  restarts is the analysed setting of van Oorschot and Wiener (1999); restart
  points are exactly uniform, which can only help the birthday bound.
- The bits of x used for partition and DP selection are standard, not
  structured by the walk; the walk does not use the negation map or the
  endomorphism, so secp256k1's automorphism group creates no fruitless cycles.
- **Toy evidence** (`experiment.py`). It runs the same procedure as above
  (shared r-adding walk, partition from the low log2 R bits of x, DP from the
  next theta bits, DP test also on start points, uniform restarts, redraw of
  M_j = O, segment cap, ceil-divided phase 1 plus a cap-step phase 2,
  sort-and-match, at most 16 pairs, re-walk, verification, input validation)
  on prime-order curves y^2 = x^3 + 7 (a = 0, j = 0, like secp256k1). For
  measurement only it also records every visited group element. That lets it
  test **exactly the event H bounds**: whether a repeated element occurs among
  the k1 phase-1 visited points (steps plus restart points). The prediction
  is 1 - exp(-k1^2 / (2 N_eff)) with N_eff = n / (1 - 2/R). Parameters were
  chosen so that restart points are below 2 % of visits, as in production
  (there about 2^-76). Each row pools the listed seeds (each a separate run
  with the listed trials):

| log2 n | R | theta | walks | cap | c | seeds x trials | mean k1 | restart share | phase-1 repeat rate | H prediction | z | algorithm success | repeats detected and solved |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 28 | 1024 | 7 | 2 | 2048 | 0.700 | 201-202 x 250 | 11572 | 0.8 % | 0.238 | 0.220 | +0.95 | 0.378 | 119/119 |
| 28 | 1024 | 7 | 2 | 2048 | 1.011 | 101-108 x 250 | 16713 | 0.8 % | 0.423 | 0.405 | +1.64 | 0.551 | 839/846 |
| 28 | 1024 | 7 | 2 | 2048 | 1.500 | 301-302 x 250 | 24796 | 0.8 % | 0.662 | 0.681 | -0.92 | 0.794 | 330/331 |
| 32 | 1024 | 9 | 4 | 8192 | 1.011 | 501-504 x 100 | 66457 | 0.2 % | 0.393 | 0.401 | -0.36 | 0.640 | 157/157 |
| 24 | 16 | 6 | 2 | 1024 | 1.011 | 401-404 x 300 | 4498 | 1.6 % | 0.436 | 0.410 | +1.82 | 0.698 | 514/523 |

  Inputs are `{"bits", "trials", "c_milli", "partitions", "theta", "walks",
  "cap_mult": 16, "seed"}` with the values above (cap = cap_mult * 2^theta).
  The toy uses cap_mult = 16 instead of production's 64 so that runs fit the
  time limit; this only shortens the cap tail (e^-16 instead of e^-64) and can
  only lower toy detection.
  The phase-1 repeat rate matches the H prediction within sampling error
  (|z| < 2) at the 0.22, 0.40 and 0.68 quantiles. That includes R = 16, where
  the 1/(1 - 2/R) factor is 12.5 % and therefore visible, so H is tested as a
  distribution, not only as a mean. The 2/R correction is too small to
  resolve at R = 1024 (0.2 %); there it rests on the cited analyses, and the
  toy shows no deviation in either direction. "Algorithm success" is higher
  than the phase-1 repeat rate because phase-2 visits (up to about 50 % extra at toy
  scale, 2^-14 in production) also produce collisions; the production bound
  credits phase-1 visits only.

  The few phase-1 repeats not solved (7/846 and 9/523) are **self-collisions
  inside one segment**, which create a cycle without a DP. Their share of all
  repeats is about 2 * 2^theta / k1 (1.5 % and 2.8 % here). In production it
  is 2 * 2^76 / 2^128 = 2^-51, inside the delta term above.

  No wrong logarithm was returned, and no degenerate relation or exceptional
  step occurred. Equal-x DP pairs with opposite y (counted separately, at most
  11 per 1200 trials) were resolved correctly through the negated relation.
  These runs are evidence for H and for detection only at toy scale
  (n <= 2^32). They say nothing about instruction counts.

Limitation: H is a heuristic, not a theorem; the extrapolation from toy
curves to r ~ 2^256 is the usual one for generic rho analyses.

## Computation ledger

Cost = H + R_instr / 1000, where H counts reference operations. All work of
all processors is included exactly once.

**Per walk step: 1 reference operation + at most 1000 instructions.**
The reference operation is the batched incomplete affine addition with the 6M
convention (3M for lambda, lambda^2, lambda(x - x'), 3M for batch
amortisation). Our step is exactly that operation. Everything else is
charged in instructions. Itemised upper count for one walk in one batch
(5x52 representation, 64-bit instructions):

| Item | Instructions |
| --- | ---: |
| load walk state (x, y: 10 limbs) and ell | 12 |
| full normalisation of x for bit extraction | 40 |
| partition index, table address, load M_j (10 limbs) | 16 |
| exceptional test x = x_j and flag | 12 |
| x_j - x, y_j - y (field subtractions, weak reduction) | 30 |
| batch bookkeeping: store difference and prefix, reload in backward pass | 30 |
| x' = lambda^2 - x - x_j and lambda(x - x') - y subtractions | 45 |
| full normalisation of x' for DP test, DP mask test on two words | 46 |
| store x', y', ell + 1, cap comparison | 16 |
| loop counters, branches, batch control | 40 |
| **amortised inversion**: one Fermat inversion per batch, at most 255 squarings + 256 multiplications + 2^14 control = 2^18 + 2^15 < 2^19 instructions, divided by m = 2^12 | 128 |
| **itemised total** | **415** |

The allowance used is **1000 instructions per step** (more than 2.4 times the
itemised count, covering any omitted copying, spilling and branch overhead).
The inversion is charged here even though the 6M convention already includes
batch amortisation, so it is double-charged on purpose. Each step therefore
costs at most **2** units.

**Number of steps.** Phase 1: W * L1 <= k1 + W. Phase 2: W * cap = 2^114.
Total steps S <= k1 + 2^114 + 2^32.

**Everything else (bounded by 2^100 units in total):**

| Item | Bound |
| --- | --- |
| Input validation ([r]Q with complete Jacobian arithmetic, <= 2^23 instr) | < 2^14 units |
| Increment table: 2R scalar multiplications + R inversions, each <= 2^23 instr | < 2^25 units (preprocessing) |
| Segment starts: at most 2^56 (abort above); each = sampling (<= 2^14 instr) + two scalar multiplications + inversion + DP test (<= 2^24 instr) | < 2^70.1 units |
| DP records: write and transfer 16 words to the global log, <= 2^10 instr each, <= 2^55 records | < 2^55 units |
| Merge sort of <= 2^55 records, 55 passes, <= 256 instr per record per pass, then linear scan | < 2^59 units |
| Collision resolution: <= 16 pairs * 2 segments * 2^82 steps, unbatched complete arithmetic with coefficient updates, <= 2^20 instr per step; plus 16 inversions mod r and 16 verifications | < 2^97.1 units |
| Randomness, output and control | < 2^60 units |

**Total.**

```
T <= 2 * (k1 + 2^114 + 2^32) + 2^100
  = 2 * 2^128.019214 + 2^115 + 2^33 + 2^100
  < 2^129.0194
  < 2^129.03.
```

The claim is the rounded-up bound **129.03**. The rounding slack
(~0.0106 bits, about 0.7 % of the work) is not needed by any item above.

## Memory, preprocessing, advice

- Walk state: 2^32 walks * (x, y, a0, b0, ell, batch scratch) <= 2^32 * 256
  bytes = 2^40 bytes.
- DP log: <= 2^55 records of 128 bytes = 2^62 bytes; merge-sort buffer of the
  same size. Total peak memory < 2^63 + 2^40 + 2^24 < 2^63.01 bytes.
  Each individual store is below 2^64 bytes, so ordinary 64-bit addressing is
  sufficient; all accesses are sequential appends, sorted passes or scans.
- Preprocessing (increment table): < 2^25 units, included in the total.
- Nonuniform advice: none (2^0 bytes).

## Objections considered

- *"Rho is heuristic."* Yes; H is disclosed with literature support, a
  1/(1 - 2/R) safety factor and toy measurements of the exact algorithm. The
  rules allow supported plausible heuristics.
- *"Fruitless cycles."* None: no negation map, no endomorphism, no
  canonicalisation. A walk can only return to a previous point by a genuine
  self-collision, which is covered by the 2^-36 term.
- *"Undetected collisions near the end of the run."* Phase 2 gives every walk
  cap additional steps and is fully charged (2^115 units, 2^-13.0 of the total).
- *"Batch inversion with exceptional or restarting walks."* Exceptional walks
  contribute the factor 1; restarts happen between batch rounds; the batch
  inversion is additionally charged at 128 instructions per step.
- *"Memory addressing beyond 64 bits."* Not needed: every array is < 2^64 bytes.
- *"Toy experiment differs from production."* It runs the same procedure with
  the same predicates (low bits for the partition, the next theta bits for the
  DP), only with small n, theta and walk counts; it is evidence at that scale
  only.
