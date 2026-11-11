# secp256k1-ecdlp: Pollard rho with negation map and GLV endomorphism

Status: new submission draft; no qualification or emitted score. This is an
analytic cryptanalysis argument, not a claim of practical execution. The
normative inputs/outputs are those of the track's `target-profile.json` and
`classical-group-work-v2.json`. No restricted challenge distribution, advice,
random-oracle premise or quantum operation is used. One disclosed heuristic
is used (the random-walk / random-mapping model for the quotient walk,
Section 4); everything else is exact arithmetic.

## 1. The claim

Let `p` be the secp256k1 field prime, `n` the (prime) group order, `P` the
standard generator from the target profile, and `x` uniform in `[0, n)`.
The input is `Q = [x]P`; the output is `x`. The Pollard rho algorithm
below, with the negation map and the GLV endomorphism (Section 3),
recovers `x` with probability at least `0.544` (hence `>= 0.39`) at the
stated mean step budget, with **expected** total work

- `log2(W) <= 128.140687` (exact rational score `128140687/1000000`),
- preprocessing `<= 2^15.0` target-group-operation equivalents (included),
- peak memory `<= 2^102.62` bytes,
- nonuniform advice `0` bytes.

All logarithms below are rounded **up** to six decimal places; the exact
arithmetic in Sections 5-6 is the audit authority. The reference price is
`C = 1000` (batched incomplete affine addition, 6M convention). The
iteration count is an expectation under the disclosed random-walk
heuristic (Section 4); the score is the log2 of the expected total work
at the mean budget, which is the standard accounting for Pollard rho.

**A correction to a naive reading, stated up front.** The textbook
statement of Pollard rho with a 6-fold symmetry is
`sqrt(pi*n/12) ~= 2^127.03` *group operations*. That figure counts only
the group operations. Under this challenge's cost model *all* work must
be charged once, including the per-step scalar arithmetic (coefficient
updates, endomorphism evaluations, canonicalisation), the fruitless-cycle
allowance, the distinguished-point sorting, and the walk restarts. When
all of that is charged honestly (Sections 5-6), the bound is
`2^128.140687`, not `2^127.03`. This submission claims the honest
`2^128.140687`. It still beats the deterministic BSGS control
(`2^129.000802`, Section 8) by about 0.86 bits, and it beats plain rho
(`2^128.559637`) and negation-only rho (`2^128.206478`).

## 2. Number-theoretic facts (all verified by direct computation)

```
p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
  = 2^256 - 2^32 - 977
n = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
```

`n` is prime. `log2(n) = 256 - eps` with `eps < 2^-127`, so `log2(n)`
rounds to `256.000000` at six decimals and every `sqrt(pi*n/(2m))`
figure below is unaffected at the stated precision (the exact `n` was
used in the script that produced the Section 6 numbers).

The secp256k1 curve `y^2 = x^3 + 7` over `F_p` has the efficiently
computable endomorphism `phi(x, y) = (beta*x, y)` where

```
beta  = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
lambda = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
```

verified by direct computation:

- `beta^3 = 1 (mod p)`, `beta != 1 (mod p)`, so `beta` has exact order 3
  in `F_p^*` and `phi` has exact order 3 as a curve endomorphism;
- `lambda^2 + lambda + 1 = 0 (mod n)`, hence `lambda^3 = 1 (mod n)` and
  `lambda != 1 (mod n)`, so `lambda` has exact order 3 in `F_n^*`;
- `phi(P) = [lambda]P` and `phi^2(P) = [lambda^2]P`, checked by direct
  elliptic-curve arithmetic on the generator `P` from the target
  profile; hence `phi(R) = [lambda]R` for every `R` in `<P>` (both sides
  are endomorphisms agreeing on the generator of the prime-order group).

So `phi` is the usual GLV endomorphism (Gallant-Lambert-Vanstone 2001).
The six-element classes `{+-phi^j(R) : j = 0, 1, 2}` partition `<P>`
except for the identity (whose class is `{O}`) and the points of order
2 (none exist: `n` is odd, so no point satisfies `R = -R` except `O`).

## 3. The algorithm

Notation: `[a]P` is scalar multiplication of a point `P` by an integer
`a` (taken mod `n`). Group operations are secp256k1 additions/doublings.
All field arithmetic is in `F_p`; all scalar arithmetic is mod `n`.

### 3.1 Distinguished-point rho, parallel walks

Fix a 32-entry `r`-adding table: 32 random pairs `(a_k, b_k)` in
`[0, n)` and the points `T_k = [a_k]P + [b_k]Q`. A walk state is
`(R, a, b)` with the invariant `R = [a]P + [b]Q`. A step from a state
whose (canonicalised) point is `R'` with coefficients `(a', b')` is:
`k = partition(R')` (low 5 bits of the x-coordinate of `R'`),
`R <- R' + T_k`, `a <- a' + a_k (mod n)`, `b <- b' + b_k (mod n)`.

A point is *distinguished* if the low 32 bits of its x-coordinate are
zero (fraction `2^-32` of points). Walks run until they reach a
distinguished point, which is stored (with its `(a, b)`); the walk then
restarts from a fresh random start. `2^20` walks run in parallel in
rounds: in each round every active walk takes one step, and the `2^20`
independent additions of the round are batched (Section 5). This is the
standard van Oorschot-Wiener parallelisation.

A collision is two stored distinguished points with the same `R` but
different `(a, b)`: then `[a1]P + [b1]Q = [a2]P + [b2]Q`, i.e.
`(a1 - a2) = (b2 - b1)*x (mod n)`. If `b2 != b1 (mod n)` (the contrary
case is a trivial collision with `a1 = a2`, see Section 4), then
`x = (a1 - a2) * (b2 - b1)^-1 mod n`. The candidate is verified by
computing `[x]P` and accepting only if it equals `Q`.

Edge case: `Q = O` (i.e. `x = 0`) is detected by a free identity check
before the walk starts. For `x != 0` the analysis below applies.

### 3.2 Negation map and GLV quotient walk

Before each step, the current point `R` is replaced by the canonical
representative of its class `{+-phi^j(R) : j = 0, 1, 2}`: the class
element whose x-coordinate is smallest (as an integer in `[0, p)`),
with ties broken by taking the point with the smaller y-coordinate.
Concretely: compute `phi(R) = (beta*x_R, y_R)` and
`phi^2(R) = (beta^2*x_R, y_R)` (one `F_p` multiplication each), compare
the three x-coordinates (and, on an x-tie, the y-coordinates against
their negations), and select the representative `R' = [s]R` with
`s in {+-1, +-lambda, +-lambda^2}`. The coefficients transform the same
way: `(a', b') = (s*a mod n, s*b mod n)`, which needs the two scalar
multiplications by `lambda^j` (two `F_p`-sized multiplications at 179
instructions each, Section 5). The step of Section 3.1 is then taken
from `(R', a', b')`.

The walk therefore lives on the quotient by the 6-element classes: the
effective search space has size `n/6`, and the expected number of steps
to a collision is `sqrt(pi*n/12)` (Section 4). The classes have size
exactly 6 for every `R != O` (Section 2), so no class is over- or
under-counted; the only fixed point of the quotient is `O`, which a
walk reaches with probability `2^-256` per step and which is handled by
treating `O` as its own class (partition index 0, no canonicalisation).

### 3.3 Fruitless cycles

The negation map is known to admit *fruitless cycles*: short cycles of
the iteration function that a walk can enter and never leave, wasting
steps (Bernstein-Lange-Schwabe 2011). The countermeasure used here is
the standard one: each walk keeps the last few visited points, and if a
point repeats without an intervening distinguished point, the walk is
aborted and restarted from a fresh random start. The cost of aborted
walks is charged by the 5% fruitless-cycle allowance in Section 5,
which is generous against the literature's estimate of roughly one to
two percent of steps lost to fruitless cycles under this countermeasure
(see Section 10 for the disclosure). Cycle detection itself (comparing
against a short history) is inside the 32-instruction distinguished-point
check allowance and the per-step load/addressing headroom of Section 5.

### 3.4 Batching schedule for the `C = 1000` price

The cost model defines `C = 1000` as the price of a *batched*
incomplete affine addition under the 6M convention, with no standalone
inversion charged per addition; this proof charges every group
addition/doubling at exactly that price, as the model prescribes. The
batching is concrete: in each round the `2^20` independent additions
(one per active walk) have all their denominators available up front
for one batch inversion, exactly as the cost model defines `C`. The
number of rounds is `2^127.033267 / 2^20 = 2^107.033267`; walks that
reach a distinguished point are replaced by fresh walks at round
boundaries. This is a material premise of the submitted bound: it
relies on the cost model's batched-price convention. If a judge instead
required a standalone affine-addition price (about 75 units per
addition, i.e. ~75,000 instructions for a Fermat inversion), the total
becomes `2^133.352571` — worse than the BSGS control (Section 9). The
submitted bound `128.140687` uses the model's batched price.

## 4. Correctness and success probability

**Heuristic H1 (random-walk / random-mapping model, disclosed).** The
quotient walk of Section 3.2 is modelled as a random mapping on a set
of `n/6` classes: each step lands on a fresh class uniformly at random,
independently of the past. This is the standard heuristic for Pollard
rho (Pollard 1978; Bernstein-Lange-Schwabe 2011 for the negation map;
Galbraith-Ruprai for the GLV quotient). It is *not* a theorem for this
concrete iteration function: the `r`-adding table is fixed, the
partition function is deterministic, and the quotient structure could
in principle bias the walk. The heuristic is supported by the usual
evidence: `r`-adding walks with a 32-entry table are empirically
indistinguishable from random mappings at cryptographic sizes, and the
same heuristic underlies every published rho complexity estimate,
including the textbook `sqrt(pi*n/2)` figure. A reviewer who rejects
H1 rejects the claim; the claim does not pretend otherwise.

Under H1, the number of steps `T` to the first repeated class is the
birthday collision time for a set of size `N = n/6`: `E[T] =
sqrt(pi*N/2) = sqrt(pi*n/12)`, i.e. `log2 = 127.033267`. For a budget of
`T` steps the collision probability is `1 - exp(-T^2/(2N))`; at the mean
budget `T = sqrt(pi*N/2)` this is `1 - exp(-pi/4) ~= 0.544 >= 0.39`.
At a `1.5x` budget the probability is `1 - exp(-(pi/4)*2.25) ~= 0.829`
for `+log2(1.5) = 0.585` bits (`2^128.725650`, still beating the control),
so the claim is not balanced on the success bound.

**Collision resolution.** Given a genuine collision (same `R`, different
`(a, b)`), the recovered `x` is correct by the algebra of Section 3.1,
and the verification `[x]P = Q` confirms it. A *trivial* collision
(`b1 = b2 (mod n)`, hence also `a1 = a2`) carries no information; the
coefficients at a collision are effectively uniform in `[0, n)`, so
this happens with probability `~2^-256` per collision, and the walk
simply continues (charged inside the 5% allowance). The stored
distinguished points are compared exactly (96-byte records), so no
false collision is ever accepted: a wrong `x` fails the verification
and the walk continues.

**Deterministic parts.** Everything except H1 is exact: the class
structure (Section 2), the coefficient invariant `R = [a]P + [b]Q`
(preserved by construction through canonicalisation and table steps),
the distinguished-point fraction `2^-32`, and the final verification.
Success probability `>= 0.544` is therefore the birthday probability
under H1, not a bound on an unmodelled failure mode.

## 5. Per-step cost charges (instructions, C = 1000)

Work is counted in scalar AArch64 instruction equivalents, then divided
by the fixed `C = 1000` to give group-operation equivalents, exactly as
the cost model prescribes. One group addition/doubling is one reference
operation, i.e. 1000 instructions. All other work is charged at
justified instruction counts; nothing is charged twice (a reference
operation is never also charged through its internal instructions).

| Work item | Charge | Justification |
| --- | --- | --- |
| Group addition / doubling | 1000 instructions (1 unit) | Reference operation, fixed `C` (batched incomplete affine addition, 6M convention; Section 3.4). |
| Partition index (low 5 bits of x) | 16 instructions | A few loads, a mask and a table-index computation. |
| Coefficient updates `a, b mod n` (2 x 64) | 128 instructions | Two 256-bit modular additions: 4-5 limb adds plus a conditional subtraction each; 64 instructions each is a generous upper bound. |
| Distinguished-point check (low 32 bits of x zero) | 32 instructions | Load the low limb of x, compare against zero, branch; includes the short-history comparison for cycle detection (Section 3.3). |
| Negation canonicalisation | 64 instructions | Compare y against `p - y` (or equivalently test the sign bit convention) and conditionally negate; a few field operations. |
| `phi`, `phi^2` evaluations (2 `F_p` mults @ 179) + representative selection (96) + `lambda^j` coefficient multiplications (2 x 179) | 812 instructions | `phi(R)` needs one `F_p` multiplication (`x * beta`); `phi^2(R)` one more (`x * beta^2`); each `F_p` multiplication is charged at the organizer's pinned measurement of 179 instructions for secp256k1 field multiplication (5x52). Representative selection compares three x-coordinates (12 limb compares) plus branches: 96 instructions. The coefficient update multiplies `a` and `b` by `lambda^j mod n` (two 256-bit multiplications at 179 instructions each, same bound as a field multiplication). |
| **Per-step instructions (GLV + negation)** | **2052** | Sum of the rows above. |
| Fruitless-cycle allowance (x1.05, negation variants only) | x1.05 | Covers aborted walks and restarts from Section 3.3; the literature puts the expected loss at roughly one to two percent of steps, so 5% is generous headroom. |
| **Per-step units (GLV + negation)** | **2.1546** | `2052 * 1.05 / 1000`. |

For comparison, the per-step charges of the two simpler variants are:
plain rho `1176` instructions (`1.176` units, no canonicalisation, no
`phi` work, no cycle allowance); negation-only rho `1240` instructions
x1.05 = `1302` instructions (`1.302` units).

**Other charged work.** The `r`-adding table build: 32 entries x 2
scalar multiplications each (one for `[a_k]P`, one for `[b_k]Q`), each
scalar multiplication at most 512 group operations = 512 units, i.e.
`32 * 2 * 512 = 32768 = 2^15.0` units. Walk restarts: each restart
costs the same 2 scalar multiplications (1024 units); the expected
number of walks is `2^127.033267 / 2^32 = 2^95.033267`, so restarts cost
`2^105.033267` units in total. Final verification lump: `2^22`
instructions (`2^12.03` units), covering the verification scalar
multiplication `[x]P` (at most 512 group operations = 512 units), the
modular inversion for collision resolution, encoding checks and output.
Distinguished-point sorting: `N * L * 448` instructions with
`N ~= 2^95.033267` stored points and `L = ceil(log2 N) = 96`, following
the secp256k1 baseline's merge-sort formula (Section 6). A hash table
is deliberately **not** used: its total work is only an expected value
(collision-dependent), while the sorting charge is a deterministic
upper bound given `N`.

## 6. Total work, and the score

Expected number of steps to a collision (Section 4):
`T = sqrt(pi*n/12)`, `log2(T) = 127.033267`. With per-step units
`2.1546` (Section 5),

```
W = T * 2.1546                                  (walk steps, dominant)
  + 2^15.0                                      (table build)
  + 2^105.033267                                (walk restarts)
  + 2^100.459800                                (sort distinguished points)
  + 2^12.03                                     (verification lump)
```

in target-group-operation equivalents, all work charged exactly once.
Numerically (log2, rounded up to six decimals):

| Component | log2 (units) |
| --- | --- |
| Walk steps `127.033267 + log2(2.1546)` | **128.140687** |
| Sort distinguished points (`N*L*448/1000`) | 100.459800 |
| Walk restarts | 105.033267 |
| Table build (preprocessing) | 15.000000 |
| Verification lump | 12.030000 |
| **Total `log2(W)`** | **128.140687** |

The non-walk terms together add less than `2^-27` relative to the walk
term, so they do not move the total at six decimals. Direct check:
`2^(128140687/1000000) >= W` (the true `log2(W)` including all small
terms is `128.1406868799...`, so the submitted rational
`128140687/1000000` is a ceiling of it). The exact rational is already
reduced (`gcd(128140687, 1000000) = 1`).

The two simpler variants, for comparison:

| Variant | Steps (log2) | Per-step units | Total (log2) |
| --- | --- | --- | --- |
| Plain rho | 128.325749 | 1.176 | 128.559637 |
| Negation-only rho | 127.825749 | 1.302 | 128.206478 |
| **GLV + negation rho (claimed)** | **127.033267** | **2.1546** | **128.140687** |

The GLV edge over negation-only is only `0.065791` bits: the `sqrt(6)`
iteration saving (1.292 bits) is almost eaten by the `phi`/`lambda`
overhead (~0.73 bits of extra per-step cost). Negation-only
(`128.206478`) is the robust fallback if a reviewer disputes any `phi`
charge; plain rho (`128.559637`) is the simplest fallback. All three
beat the BSGS control (`129.000802`).

## 7. Memory, preprocessing, advice

- **Peak memory** (upper bound): the distinguished-point store holds
  `N ~= 2^95.033267` records of 96 bytes (x-coordinate + `a` + `b`),
  i.e. `2^101.618230` bytes; merge sort needs a second buffer of the
  same size, so the store plus sort buffers are at most `2^102.618230`
  bytes, stated as `2^102.62` bytes. Walk state for `2^20` parallel
  walks (one point + two coefficients + partition state each, at most
  256 bytes per walk) is at most `2^28` bytes. The `r`-adding table is
  32 points (`2^11` bytes). Total peak memory is at most `2^102.62`
  bytes. Memory has no tie-break; it is stated for completeness.
- **Preprocessing**: the `r`-adding table build, `2^15.0` units,
  included in the total above. No preprocessing is reused across
  targets and no credit is taken for it. (The distinguished-point sort
  is online work during the walk, not preprocessing.)
- **Nonuniform advice**: `0` bytes. The `r`-adding table is built by
  the attack itself from the public inputs (`P`, `Q`) and its
  construction is charged in full above; it is not an embedded
  precomputed answer. No other advice is used.

## 8. Comparison with the BSGS control and generic rho

For a group of prime order `n`, plain Pollard rho needs about
`sqrt(pi*n/2)` group operations, i.e. `log2 = 127.754293` here
(`128.559637` with the honest per-step charges of Section 5); with
the negation map about `sqrt(pi*n/4)`, i.e. `log2 = 127.254293`
(`128.206478` honestly charged). The textbook group-operation-only
figures are lower than this submission's bound because they charge
only group operations, not the per-step scalar work — the honest
charging of Section 5 is what moves the bound up, as stated in
Section 1. On the honestly-charged comparison this submission's
`128.140687` beats plain rho by `0.418950` bits and negation-only rho
by `0.065791` bits, and it beats the track's deterministic BSGS
control (`129.000802`) by `0.860115` bits. The control was qualified only under a retired
evaluator, so this submission needs fresh review regardless; it is not
compared against a live best, because there is no prior accepted score
on this track beyond the control itself: the selected template is the
unscored BSGS baseline.

## 9. Sensitivity: what if the charges are higher?

- **Standalone addition price (load-bearing, disclosed).** If every
  group addition were instead charged a standalone affine price of
  about 75 units (~75,000 instructions, a Fermat inversion per
  addition, contradicting the cost model's own definition, which has
  no standalone inversion term), the per-step charge becomes
  `(75 + 1052/1000) * 1.05 = 79.8546` units and the total becomes
  `log2(W) = 133.352571` — **worse** than the BSGS control. The claim
  stands or falls on the cost model's batched 6M convention plus the
  `2^20`-parallel-walk schedule of Section 3.4. (Contrast the Cheon
  submission on the BLS track, which survived standalone pricing.)
  The submitted bound `128.140687` uses the model's batched price.
- **`F_p` multiplication at twice the pinned price.** If each field
  multiplication cost 358 instructions instead of the organizer's
  pinned 179, the per-step charge becomes `2906.4` instructions
  (`2.9064` units) and the total `2^128.572500` — still beating the
  control by 0.43 bits. The break-even point is an `F_p`
  multiplication of about 597 instructions, i.e. about 3.34x the
  pinned measurement; the claim survives any field-multiplication
  price up to that point.
- **Sort at twice the cost.** Doubling the merge-sort charge moves
  the sort term to `2^101.46` units, still negligible against the
  `2^128.14` walk term; the total is unchanged at six decimals.
- **Fruitless-cycle allowance.** Doubling the allowance from 5% to
  10% adds `log2(1.10/1.05) = 0.067` bits, giving `2^128.208` — still
  beating the control. The margin is not balanced on the 5% figure.
- **Success probability.** At a `1.5x` step budget the success
  probability rises to `~0.829` for `+0.585` bits (`2^128.725650`),
  still beating the control; the claim is not balanced on the success
  bound either.

## 10. Heuristics, limitations, and what the experiment shows

- **Heuristic H1 (random-walk model):** the only heuristic in the
  claim (Section 4). The quotient walk is modelled as a random
  mapping on `n/6` classes. Support: the standard rho analysis and
  the empirical behaviour of `r`-adding walks with a 32-entry table.
  Limitation: it is not a theorem for this concrete iteration
  function; a reviewer demanding a proven mixing bound would reject.
  The fruitless-cycle handling (Section 3.3) is a second disclosed
  heuristic-adjacent ingredient: the 5% allowance is a charged bound
  on the loss, not a measured figure, and the literature estimate it
  rests on is roughly one to two percent.
- **Limitations.** The bound is an analytic bound on *expected* work,
  not a measured runtime; no full-size execution is claimed or
  possible. Memory is astronomical (`2^102.62` bytes), as with the
  BSGS control (`2^136` bytes); this is not a practical attack. The
  score is the log2 of expected work at the mean budget, so a single
  run could cost more (or less); the success probability `0.544` is
  the corresponding birthday figure, and the `1.5x`-budget sensitivity
  in Section 9 bounds the cost of buying a higher success
  probability.
- **Experiment.** `experiment.py` implements the exact walk logic of
  Section 3 (negation + GLV canonicalisation, `r`-adding table,
  distinguished points, cycle detection with restart, collision
  resolution) on the additive group `Z/qZ` for the declared input
  `(q, secret) = (10009, 1234)`, where `q` is prime, `q = 1 (mod 3)`,
  and the endomorphism is multiplication by a cube root of unity
  `lambda mod q` (found by search in the experiment). It recovers the
  secret in every run and also checks the class structure (orbits of
  size 6 except 0) and `lambda^3 = 1 (mod q)`. The toy uses a hash
  dict, not merge sort, and a 32-entry table with small distinguished
  threshold — it validates algorithm logic only. It does not establish
  production cost bounds, memory bounds, or secp256k1 arithmetic;
  those rest on the analytic accounting in Sections 5-6.
- **What could still be wrong.** The main risks are (a) heuristic H1
  failing for this concrete walk (Section 4), (b) the batched-price
  premise for `C = 1000` (Section 3.4 states the schedule and Section
  9 bounds the standalone-price consequence at `2^133.352571`), (c) a
  charge the reviewer judges too low (Section 9 bounds the
  consequence of the largest such charges; the field-multiplication
  break-even is ~3.34x the pinned measurement), and (d) an arithmetic
  slip in the Section 6 numbers (they were computed by script from
  the formulas in Section 5; the formulas, not just the numbers, are
  stated so they can be rechecked).

## References

- J. M. Pollard, "Monte Carlo methods for index computation (mod p)",
  Math. Comp. 32 (1978), pp. 918-924. The rho algorithm and the
  `sqrt(pi*n/2)` analysis.
- D. J. Bernstein, T. Lange and P. Schwabe, "On the correct use of the
  negation map in the Pollard rho method", PKC 2011, LNCS 6571,
  pp. 128-146. The negation map, fruitless cycles and countermeasures.
- R. P. Gallant, R. J. Lambert and S. A. Vanstone, "Faster point
  multiplication on elliptic curves with efficient endomorphisms",
  CRYPTO 2001, LNCS 2139, pp. 190-200. The GLV endomorphism.
- S. D. Galbraith and R. S. Ruprai, "Using equivalence classes to
  accelerate solving the discrete logarithm problem in a short
  interval", PKC 2010, LNCS 6056, pp. 368-383. Quotient walks with
  the GLV endomorphism.
- P. C. van Oorschot and M. J. Wiener, "Parallel collision search
  with cryptanalytic applications", J. Cryptology 12 (1999),
  pp. 1-28. Distinguished points and parallel walks.
- A. Menezes, P. van Oorschot and S. Vanstone, Handbook of Applied
  Cryptography, Section 3.6 (generic DLP algorithms).
- The track's own `target-profile.json`,
  `classical-group-work-v2.json` and `tracks/docs/COST_MODEL.md` for
  the target, the reference price `C = 1000`, and the accounting rules
  applied in Sections 5-6. The values of `beta` and `lambda` in
  Section 2 were verified by direct computation from the profile's
  `p` and `n` (including `phi(P) = [lambda]P` by direct
  elliptic-curve arithmetic).
