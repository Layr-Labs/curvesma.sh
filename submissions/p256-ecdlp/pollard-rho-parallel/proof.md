# p256-ecdlp: parallel Pollard rho on group points

Status: ready for independent exploratory review; no qualification or emitted
score is claimed by this file. This is a generic best-known-style attack
accounted under `classical-group-work-v2`, not a cryptanalytic advance.
`claim.json` binds `tracks/p256-ecdlp/target-profile.json` and
`tracks/classical-group-work-v2.json` by SHA-256. No auxiliary input, advice,
oracle or quantum step is used.

## Why this generic method

P-256 has prime order r (cofactor 1), r != p (not anomalous), embedding degree
far above any feasible pairing transfer, and a = -3 gives no known efficiently
computable endomorphism other than negation. This package uses neither:
the walk runs on the r group points directly with an R-adding iteration.

We deliberately do not use the negation map. It would halve the state space
to N = (r+1)/2 classes and save a factor sqrt(2) in expected steps, but the
correctness of the resulting walk depends on the detailed behavior of
canonical-sign collisions: Bernstein, Lange and Schwabe (PKC 2011) show that
negation-map walks enter "fruitless cycles" whose frequency depends on the
partition geometry in a way that resists a short self-contained analysis. The
verifier's shared evidence policy requires that a nonstandard external result
material to the claim be reproduced in the supplied evidence; we cannot
reproduce that analysis compactly and our scoped toy curve is too small to
bind its constants. The plain walk removes the entire question: every state
is a point, every repeat of a point is a solvable relation, and the only
heuristic left is the standard random-mapping idealization whose birthday
calculation is carried out in full below. The cost of this choice is a
factor of at most sqrt(2) in the step budget, about 0.5 in log2 score,
which we pay in exchange for a success argument that needs no external
machinery beyond the explicitly quantified heuristic M1.

## Algorithm

Notation: p, r, P, a, b from the protected profile; input Q. Constants of the
run: partition count R = 2048, walk count w = 2^40 (2^30 processor cores each
advancing a batch of B = 1024 walks), distinguished-point probability
theta = 2^-32, point budget k = 46*2^123, per-walk length
L = k/w + 2^37 = 46*2^83 + 2^37, total steps K = w*L.

1. Input handling. Decode Q canonically and check it lies on the curve
   (on a curve of prime order, an on-curve point automatically has order
   dividing the prime r, i.e. is in the group). If Q = infinity output 0.
   If Q = P output 1; these fast paths are optional and do not affect the
   bound.
2. Partition table. For i in 0..R-1 draw c_i, d_i uniformly in [0,r) by
   rejection sampling 256-bit words (rejection probability < 2^-32 per draw,
   since r > 2^255) and compute M_i = [c_i]P + [d_i]Q in affine form. If any
   M_i is infinity, redraw it.
3. Walk. The partition index is h(X) = (x >> 8) mod R on the affine
   x-coordinate, a well-defined function on points. From X with i = h(X) set
   Y = X + M_i. Look-ahead (Wiener-Zuccherato style, applied to the plain
   walk): if h(Y) = i, instead use Y = X + M_{(i+1) mod R}, and keep
   incrementing the index while it repeats, up to R tries (then restart the
   walk). This makes the step depend on the preceding index and removes
   immediate backtracking at cost one extra addition with probability
   1/(R-1) per step.
4. Start points. Walk j draws seed coefficients (u_j, v_j) uniformly in
   [0,r) and starts at X_j = [u_j]P + [v_j]Q, which is uniform in the group.
   Coefficients are not tracked during the main walk; only the seed is
   retained.
5. Batching. Each core advances its 1024 walks in lockstep, using
   Montgomery's simultaneous-inversion trick so one field inversion serves
   1024 affine additions.
6. Distinguished points. A point is distinguished when the low 32 bits of
   its affine x-coordinate are zero. On reaching one, the record
   (x, walk id, step index), 48 bytes, is inserted into a shared
   open-addressing hash table of capacity ndp_cap defined below, keyed on x:
   equal points have equal x, and distinct points share an x only when they
   are negatives, which the key comparison resolves by accepting a match and
   letting step 7 determine usefulness. If the same key was already stored,
   process the collision (step 7) and restart the newer walk from a fresh
   seed.
7. Collision solving. On a key match between stored walk i and current
   walk j, re-walk both from their seeds to the matched point, now tracking
   coefficients: the state (alpha, beta) with X = [alpha]P + [beta]Q is
   updated by (alpha, beta) += (c_i, d_i) on each step, including the
   look-ahead substitutions, which are ordinary additions with known
   indices. At the match, [alpha_i]P + [beta_i]Q =
   +/- ([alpha_j]P + [beta_j]Q) with the sign determined by comparing the
   y-coordinates (equal or summing to p mod p). Write the relation as
   [alpha_i - s*alpha_j]P = [s*beta_j - beta_i]Q, s = +/-1. If
   s*beta_j - beta_i != 0 mod r, output
   x = (alpha_i - s*alpha_j) / (s*beta_j - beta_i) mod r, verify [x]P = Q,
   and return it. The relation is useless only if both coefficient
   differences vanish mod r. At most 16 collisions are processed; after the
   16th, or when all walks have run L steps, or the DP cap is reached, stop
   with failure.

Every repeat is potentially useful: because the state space is the r group
points themselves rather than sign classes, two walks that reach the same
point carry independent coefficient vectors and the relation fails only when
the two vectors coincide mod r. There is no analogue of a fruitless cycle:
a walk that re-enters a previously visited point produces the same kind of
solvable relation, and a walk trapped in a cycle containing no
distinguished point is a bounded budget loss, not a silent corruption.

## Success argument

Modeling assumption M1 (random-mapping idealization). For collision
purposes the R-adding iteration above, driven by the high-bit partition on
the affine x-coordinate, behaves like a uniform random mapping on the r
points: the joint distribution of the w walks' iterates is modeled as
independent uniform draws until the first repeat. Under M1 the probability
that k = w*L_main points (L_main = k/w per walk) contain no repeat is
exactly the birthday bound prod_{i<k} (1 - i/r) <= exp(-k(k-1)/(2r)), and
the expected number of points to the first repeat is sqrt(pi r/2).

Status of M1. M1 is a heuristic, not a theorem, and it is the only
heuristic the success claim depends on. Its support is: (a) the standard
r-adding-walk analysis (Teske, ANTS 1998; Brent-Pollard) which finds that
adding walks on a fixed random partition with R >= 20 deviate from the
random-mapping idealization by a small constant factor in expected
collision time, quantified as (1 - 1/R)^(-1/2) for the idealized partition
model — cited here as qualitative support only, since the claim does not
rely on the exact factor; (b) the packaged toy experiment, which runs this
exact walk construction end-to-end at its own scope and reports the
measured steps-to-solve against the sqrt(pi r/2) idealization; and
(c) the quantified margin below: the bound survives a degradation of the
collision exponent by a factor of 2.09. M1 is declared as an assumption
whose failure would falsify the success bound, not asserted as proved.

Calculation. k = 46*2^123 gives k^2/(2r) = 2116*2^246/(2r) =
1058*2^246/r. Since r < 2^256, k^2/(2r) > 1058/1024 = 1.0332031 (in fact
1.0332031252...). Under M1 the probability of at least one repeat among
the k main-walk points is at least 1 - e^(-1.0332031) = 0.64413.

Success requires at least 0.39, i.e. a collision exponent of
-ln(1 - 0.39) = 0.4943. The budget therefore carries a degradation margin
of 1.0332031/0.4943 = 2.09: even if the true walk required 2.09 times more
points than the random-mapping idealization for the same collision
probability, success would still exceed 0.39.

From a repeat to a solution. A repeat between two different walks merges
them: after the merge both walks follow the same deterministic trajectory.
Each merged walk reaches a distinguished point within its remaining 2^37
steps except with probability (1 - 2^-32)^(2^37) <= e^(-32) < 2^-46, so the
pair registers a table collision except with probability < 2^-45.
A repeat inside one walk (the walk colliding with itself) needs a
self-intersection of a walk of length at most L < 2^88.6 in a set of
r > 2^255 points: probability < L^2/(2r) < 2^-78 per walk and, by the
union bound over 2^40 walks, < 2^-38 in total. When it does occur it is
not a loss: it is detected at the next distinguished point of the shared
cycle (probability >= 1 - 2^-46 that the cycle or its tail contains one
within the remaining budget) and yields the same solvable relation. The
only true losses are a cycle with no distinguished point in the remaining
budget (< 2^-46 per self-repeat) and the collision being useless.
Uselessness requires both (alpha_i - s*alpha_j, beta_i - s*beta_j) == 0
mod r; since the two seeds and the accumulated index sequences are
independent mod r, this has probability at most 2^-250, and 16 collisions
cover the eventuality that early collisions involve the same walk pair.
The DP cap ndp_cap = 2*floor(K/2^32) + 2 exceeds the binomial(K, 2^-32)
count except with probability below e^(-2^90) (Chernoff at twice the mean).
Hash probing longer than 64 at load factor below 1/2 can only delay a
detection; its probability per insertion is below 2^-50 under the modeled
uniform key distribution, and a delayed detection merely consumes budget
already charged. Summing every loss gives success >=
0.64413 - 2^-39 > 0.39 under M1, over uniform x and fresh coins.

## Instruction ledger (scalar AArch64, 64-bit GPRs)

Field elements are 4 limbs in Montgomery form. For P-256,
-p^{-1} mod 2^64 = 1 and p has words (2^64-1, 2^32-1, 0, 2^64-2^32+1),
so each Montgomery reduction round is shift-and-add work plus one
MUL/UMULH pair.

| Kernel | Bound (instructions) | Justification |
| --- | ---: | --- |
| Field mult/square M | 256 | 4x4 schoolbook: 16 MUL + 16 UMULH + <= 64 ADDS/ADCS; 4 Montgomery rounds <= 16 each (64); final conditional subtract 12; 8 loads + 4 stores; call/return 8: total 192, charged 256 |
| Field add/sub A | 32 | 4 ADDS/ADCS, 4 SUBS/SBCS against p, 4 CSEL, 8 loads, 4 stores, call: <= 28, charged 32 |
| Field inversion I | 76,800 = 300M | Fermat x^(p-2) by a fixed chain of 255 squarings and <= 45 multiplications, plus loop control |
| Unbatched affine add/double | 78,080 = I + 4M + 8A | slope with inversion, x3, y3, infinity/opposite/doubling case tests |
| Scalar multiplication | 39,976,960 = 512 * 78,080 | double-and-add over 256 bits, unbatched |

Per main-walk step (batched), S = 2240:

| Item | Instructions |
| --- | ---: |
| Montgomery batch trick, 3M per element | 768 |
| Affine-add multiplies (lambda^2, lambda*(x1-x3)): 2M | 512 |
| 6 field subtractions | 192 |
| Amortised inversion I/1024 = 75, plus bookkeeping | 128 |
| Partition index and 8-word table load | 32 |
| Walk state and batch scratch loads/stores | 64 |
| Distinguished test | 8 |
| Look-ahead extra addition: probability 1/(R-1) of an unbatched-style add path | 64 (needs 38.2) |
| Subtotal 1768; reserve | 472 |

No cycle detection is required and none is charged: a trapped walk is a
budget loss bounded in the success argument, not a correctness hazard, and
the expected number of lost walks is < 2^-40 of the total. The charged
per-step cost S = 2240 applies to all K steps, including the 2^37-step
detection tail of every walk.

Other components, all charged in full:

- Setup: 2^24 for decoding, curve and constant checks of Q, plus
  2048 * (2 scalar multiplications + 1 addition) for the table and
  2048 * 1024 for its randomness and rejection sampling.
- Walk starts: per walk, 2 scalar multiplications, 1 addition and 1024 for
  seed randomness.
- DP handling: per record, 2^14 covering hashing, at most 64 probes of a
  48-byte record compare, and insertion; charged for the full capacity
  ndp_cap = 2*floor(K/2^32) + 2 records.
- Re-walks: at most 16 collisions, two walks each of at most L steps,
  unbatched with coefficient tracking at <= 2^17 per step (78,080 plus two
  mod-r additions and the look-ahead index bookkeeping).
- Solving, verification and output: 2^28 covering the mod-r inversion, the
  [x]P = Q check and serialization.

## Submitted computation bound

| Component | Instructions |
| --- | ---: |
| main walk K * 2240 | 1095709221485422190851295728026462558617600 |
| setup | 163924410368 |
| walk starts | 87997240494075674624 |
| DP handling | 3731963367071908510927798843506688 |
| re-walks | 1865981683535954255463899421736960 |
| solving, verification, output | 268435456 |

k = 489155902448849041228600998183166803968 points;
K = 489155902448849192344328450011813642240 steps;
L = 444884701618183673730826240 steps.

Total T = **1095709227083367241459246491658819092381696** instructions.
Dividing by the fixed C = 1000 gives log2(T/1000) = 129.6870606...,
submitted as time_log2 = **129687061/1000000 = 129.687061**. Preprocessing
(setup, independent of later walks) is 163924410368 instructions, giving
preprocessing_log2 **27.288456**, included in the total. All point
arithmetic is charged through its instruction count; no addition is priced
as one unit.

## Memory and advice

The DP hash table holds ndp_cap = 227780967228510040950183034882 records of
48 bytes (x, walk id, step index). Walk states, batch scratch and seeds
take 2^40 * 128 bytes. Code, the 2048-entry table and constants fit in a
2^30-byte reserve. Peak memory <= 10933486426968482106347347771488 bytes,
memory_log2_bytes **103.108525**. Nonuniform advice: none (2^0 bytes).

## Evidence and limitations

The cost and success arguments are analytic; the only modeling assumption
is M1, disclosed above with its support and its 2.09x margin.
`experiment.py` runs this same R-adding walk, with look-ahead and
distinguished points, end-to-end on a small prime-order Weierstrass curve
(p = 1999, prime order r = 2089, R = 16, theta = 1/8, 8 concurrent walks)
and reports measured steps-to-solve against the idealized sqrt(pi r/2)
value. With its packaged deterministic seed it recovers 32/32 instances
at a mean of 1.40x the idealized step count and zero unsolvable relations.
It is evidence only at that scope: a finite toy run does not bind
the constant in M1 at r ~ 2^256, which is exactly why the budget carries
an explicit degradation margin. The instruction bounds are analytical
upper bounds, not Callgrind measurements of a P-256 kernel. Compared to
the negation-map variant this design pays ~0.64 bits of score
(129.687 vs ~129.05 for the same confidence construction) to eliminate
the fruitless-cycle analysis entirely.

References for the generic method and its scoped support: van Oorschot-
Wiener, "Parallel collision search with cryptanalytic applications"
(J. Cryptology 1999); Teske, "Speeding up Pollard's rho method for
computing discrete logarithms" (ANTS 1998); Brent-Pollard, "Factorization
of the eighth Fermat number" (1981); Wiener-Zuccherato, "Faster attacks on
elliptic curve cryptosystems" (SAC 1998) for the look-ahead idea applied
here without the negation map; Bernstein-Lange-Schwabe (PKC 2011) cited
only to explain why the negation map was declined, not as support for any
claim.
