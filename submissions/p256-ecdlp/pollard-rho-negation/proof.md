# p256-ecdlp: Pollard rho with the negation map

Status: new submission draft; no qualification or emitted score. This is an
analytic cryptanalysis argument, not a claim of practical execution. The
normative inputs/outputs are those of the track's `target-profile.json` and
`classical-group-work-v2.json`. No restricted challenge distribution, advice,
random-oracle premise or quantum operation is used. One disclosed heuristic
is used (the random-walk / random-mapping model for the quotient walk,
Section 4); everything else is exact arithmetic.

## 1. The claim

Let `p` be the P-256 field prime, `r` the (prime) group order, `P` the
standard generator from the target profile, and `x` uniform in `[0, r)`.
The input is `Q = [x]P`; the output is `x`. The Pollard rho algorithm
below, with the negation map (Section 3), recovers `x` with probability
at least `0.544` (hence `>= 0.39`) at the stated mean step budget, with
**expected** total work

- `log2(W) <= 128.724770` (exact rational score `12872477/100000`),
- preprocessing `<= 2^15.0` target-group-operation equivalents (included),
- peak memory `<= 2^103.42` bytes,
- nonuniform advice `0` bytes.

All logarithms below are rounded **to** six decimal places; the exact
arithmetic in Sections 5-6 is the audit authority. The reference price is
the borrowed normalization unit `C = 1000` (Section 3.4): all work is
counted in scalar AArch64 instruction equivalents and divided by 1000.
The iteration count is an expectation under the disclosed random-walk
heuristic (Section 4); the score is the log2 of the expected total work
at the mean budget, which is the standard accounting for Pollard rho.

**A correction to a naive reading, stated up front.** The textbook
statement of Pollard rho with the negation map is `sqrt(pi*r/4) ~=
2^127.83` *group operations*. That figure counts only the group
operations. Under this challenge's cost model *all* work must be charged
once, including the per-step scalar arithmetic (coefficient updates,
canonicalisation), the fruitless-cycle allowance, the distinguished-point
sorting, and the walk restarts, and a P-256 point addition is not
automatically one unit: the rules require actual P-256 arithmetic to be
charged through justified instruction costs (Section 3.4). When all of
that is charged honestly (Sections 5-6), the bound is `2^128.724770`,
not `2^127.83`. This submission claims the honest `2^128.724770`. It
beats honestly-charged plain rho (`2^129.101431`, Section 8) by about
0.377 bits. There is no prior accepted score on this track: the selected
package is an unscored template, so this claim, if accepted, sets the
first score.

**No endomorphism bonus is claimed.** P-256 has no efficiently
computable endomorphism of small order available to this attack (unlike
secp256k1's GLV endomorphism): the curve has `a = -3`, `j != 0, 1728`,
and no small-degree complex-multiplication endomorphism, so the only
symmetry exploited here is the negation map. No automorphism-aware
speedup beyond the factor 2 of Section 3.2 is claimed.

## 2. Number-theoretic facts (all verified by direct computation)

```
p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
  = 2^256 - 2^224 + 2^192 + 2^96 - 1
r = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
a = 0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc
  = p - 3
b = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
```

`p` and `r` are prime. `log2(r) = 256 - eps` with `eps < 2^-31`, so
`log2(r)` rounds to `256.000000` at six decimals and every
`sqrt(pi*r/(2m))` figure below is unaffected at the stated precision
(the exact `r` was used in the script that produced the Section 6
numbers). The generator `P` from the target profile lies on the curve
`y^2 = x^3 + a*x + b` over `F_p`, and `[r]P = O` (checked by direct
elliptic-curve arithmetic), with `[2]P != O`; since `r` is prime and the
cofactor is 1, `<P>` is the full group of prime order `r`.

The negation classes `{R, -R}` partition `<P>` except for the identity
(whose class is `{O}`); there are no points of order 2 (`r` is odd, so
no point satisfies `R = -R` except `O`).

## 3. The algorithm

Notation: `[a]P` is scalar multiplication of a point `P` by an integer
`a` (taken mod `r`). Group operations are P-256 additions/doublings.
All field arithmetic is in `F_p`; all scalar arithmetic is mod `r`.

### 3.1 Distinguished-point rho, parallel walks

Fix a 32-entry `r`-adding table: 32 random pairs `(a_k, b_k)` in
`[0, r)` and the points `T_k = [a_k]P + [b_k]Q`. A walk state is
`(R, a, b)` with the invariant `R = [a]P + [b]Q`. A step from a state
whose (canonicalised) point is `R'` with coefficients `(a', b')` is:
`k = partition(R')` (low 5 bits of the x-coordinate of `R'`),
`R <- R' + T_k`, `a <- a' + a_k (mod r)`, `b <- b' + b_k (mod r)`.

A point is *distinguished* if the low 32 bits of its x-coordinate are
zero (fraction `2^-32` of points). Walks run until they reach a
distinguished point, which is stored (with its `(a, b)`); the walk then
restarts from a fresh random start. `2^20` walks run in parallel in
rounds: in each round every active walk takes one step, and the `2^20`
independent additions of the round are batched (Section 3.4). This is the
standard van Oorschot-Wiener parallelisation.

A collision is two stored distinguished points with the same `R` but
different `(a, b)`: then `[a1]P + [b1]Q = [a2]P + [b2]Q`, i.e.
`(a1 - a2) = (b2 - b1)*x (mod r)`. If `b2 != b1 (mod r)` (the contrary
case is a trivial collision with `a1 = a2`, see Section 4), then
`x = (a1 - a2) * (b2 - b1)^-1 mod r`. The candidate is verified by
computing `[x]P` and accepting only if it equals `Q`.

Edge case: `Q = O` (i.e. `x = 0`) is detected by a free identity check
before the walk starts. For `x != 0` the analysis below applies.

### 3.2 Negation map

Before each step, the current point `R` is replaced by the canonical
representative of its class `{R, -R}`: the class element whose
y-coordinate is smallest (as an integer in `[0, p)`), equivalently the
point with the smaller y-coordinate, with the coefficients negated when
the representative is `-R`: `(a', b') = (-a mod r, -b mod r)`. The step
of Section 3.1 is then taken from `(R', a', b')`.

The walk therefore lives on the quotient by the 2-element classes: the
effective search space has size `r/2`, and the expected number of steps
to a collision is `sqrt(pi*r/4)` (Section 4). The classes have size
exactly 2 for every `R != O` (Section 2), so no class is over- or
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

### 3.4 P-256 addition price and batching schedule

**The P-256 field-multiplication price (justified, not borrowed).**
The rules for this track state that the secp256k1 6M reference formula
is not a P-256 implementation bound and that actual P-256 arithmetic
must be charged through justified instruction counts divided by the
borrowed `C = 1000`. This proof therefore prices P-256 field
multiplication from its own implementation, in scalar AArch64
instruction equivalents:

- A 256-bit by 256-bit schoolbook product is 4x4 limbs: 16 64-bit
  multiply-accumulate operations plus ~16 accumulation adds, about
  40-50 instructions for the 512-bit product.
- Solinas fast reduction for `p = 2^256 - 2^224 + 2^192 + 2^96 - 1`
  expresses the product as a sum of shifted 256-bit terms (the nine
  terms `s1..s9` of the standard formulation), i.e. roughly 8-10
  256-bit additions; each 4-limb addition with carry is about 8
  instructions, so ~70-90 instructions, plus a final conditional
  subtraction of about 20 instructions.
- Loads/stores of operands and result, addressing, and loop/call
  overhead: ~60-80 instructions.

A realistic total is ~180-220 instructions. This proof charges **256
instructions per `F_p` multiplication**, a generous upper bound (about
1.4x the organiser's pinned secp256k1 5x52 measurement of 179
instructions, reflecting the extra reduction additions). All other
field operations are charged on the same basis: a 256-bit modular
addition/subtraction is at most ~12 instructions per limb operation
including the conditional reduction, and the per-step scalar charges in
Section 5 use the same generous limb counts.

**Batched affine addition: 6 field multiplications.** An incomplete
affine addition `R' + T_k` (with `R' != +-T_k`, the generic case; the
doubling case is covered by the same bound, see below) needs, by
formula, 2 field multiplications (`lambda^2` and `lambda*(x1 - x3)`)
plus ~6 field additions/subtractions (~70 instructions) for the
lambda computation and the coordinate updates. The batch-inversion
share is negligible: one Fermat inversion is ~266 field multiplications
(~68,096 instructions at the charged price), spread over a batch of
`2^20` additions it is ~0.07 instructions per addition. The direct
formula cost is therefore ~582 instructions. This proof charges **6
field multiplications = 1536 instructions** per batched addition,
about 2.6x the direct formula cost: the extra 4 multiplications of
headroom cover the batch amortisation share, the formula's
additions/subtractions, loads/stores, loop overhead, and the rare
doubling case (a doubling needs ~2 multiplications plus squarings,
within the same bound). This 6M figure is a P-256-specific upper bound
built from the 256-instruction field multiplication above; it is not
the secp256k1 6M calibration imported as a bound.

**Batching schedule for the `C = 1000` price.** The cost model defines
`C = 1000` as the price of a *batched* incomplete affine addition
under the 6M convention, with no standalone inversion charged per
addition; this proof charges every group addition/doubling at exactly
the batched price justified above, as the model prescribes. The
batching is concrete: in each round the `2^20` independent additions
(one per active walk) have all their denominators available up front
for one batch inversion, exactly as the cost model defines `C`. The
number of rounds is `2^127.825748 / 2^20 = 2^107.825748`; walks that
reach a distinguished point are replaced by fresh walks at round
boundaries. This is a material premise of the submitted bound: it
relies on the cost model's batched-price convention. If a judge instead
required a standalone affine-addition price (about 68.7 units per
addition, i.e. ~68,678 instructions for a Fermat inversion per
addition, contradicting the cost model's own definition, which has no
standalone inversion term), the total becomes `2^134.002946`
(Section 9). The submitted bound `128.724770` uses the model's
batched price.

## 4. Correctness and success probability

**Heuristic H1 (random-walk / random-mapping model, disclosed).** The
quotient walk of Section 3.2 is modelled as a random mapping on a set
of `r/2` classes: each step lands on a fresh class uniformly at random,
independently of the past. This is the standard heuristic for Pollard
rho (Pollard 1978; Bernstein-Lange-Schwabe 2011 for the negation map).
It is *not* a theorem for this concrete iteration function: the
`r`-adding table is fixed, the partition function is deterministic, and
the quotient structure could in principle bias the walk. The heuristic
is supported by the usual evidence: `r`-adding walks with a 32-entry
table are empirically indistinguishable from random mappings at
cryptographic sizes, and the same heuristic underlies every published
rho complexity estimate, including the textbook `sqrt(pi*r/2)` figure.
A reviewer who rejects H1 rejects the claim; the claim does not
pretend otherwise.

Under H1, the number of steps `T` to the first repeated class is the
birthday collision time for a set of size `N = r/2`: `E[T] =
sqrt(pi*N/2) = sqrt(pi*r/4)`, i.e. `log2 = 127.825748`. For a budget of
`T` steps the collision probability is `1 - exp(-T^2/(2N))`; at the mean
budget `T = sqrt(pi*N/2)` this is `1 - exp(-pi/4) ~= 0.544 >= 0.39`.
At a `1.5x` budget the probability is `1 - exp(-(pi/4)*2.25) ~= 0.829`
for `+log2(1.5) = 0.585` bits (`2^129.309732`, still the submitted
family), so the claim is not balanced on the success bound.

**Collision resolution.** Given a genuine collision (same `R`, different
`(a, b)`), the recovered `x` is correct by the algebra of Section 3.1,
and the verification `[x]P = Q` confirms it. A *trivial* collision
(`b1 = b2 (mod r)`, hence also `a1 = a2`) carries no information; the
coefficients at a collision are effectively uniform in `[0, r)`, so
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
by the borrowed `C = 1000` to give group-operation equivalents, exactly
as the cost model prescribes. One group addition/doubling is charged at
the batched price of Section 3.4 (1536 instructions). All other work is
charged at justified instruction counts; nothing is charged twice (a
reference operation is never also charged through its internal
instructions).

| Work item | Charge | Justification |
| --- | --- | --- |
| Group addition / doubling | 1536 instructions (1.536 units) | Batched affine addition at 6 field multiplications x 256 instructions (Section 3.4). |
| Partition index (low 5 bits of x) | 16 instructions | A few loads, a mask and a table-index computation. |
| Coefficient updates `a, b mod r` (2 x 64) | 128 instructions | Two 256-bit modular additions: 4-5 limb adds plus a conditional subtraction each; 64 instructions each is a generous upper bound. |
| Distinguished-point check (low 32 bits of x zero) | 32 instructions | Load the low limb of x, compare against zero, branch; includes the short-history comparison for cycle detection (Section 3.3). |
| Negation canonicalisation | 64 instructions | Compare y against `p - y` (or equivalently test the sign bit convention) and conditionally negate; a few field operations. |
| **Per-step instructions (negation rho)** | **1776** | Sum of the rows above. |
| Fruitless-cycle allowance (x1.05) | x1.05 | Covers aborted walks and restarts from Section 3.3; the literature puts the expected loss at roughly one to two percent of steps, so 5% is generous headroom. |
| **Per-step units (negation rho)** | **1.8648** | `1776 * 1.05 / 1000`. |

For comparison, the per-step charge of plain rho (no negation map, no
canonicalisation, no cycle allowance) is `1536 + 16 + 128 + 32 = 1712`
instructions (`1.712` units).

**Other charged work.** The `r`-adding table build: 32 entries x 2
scalar multiplications each (one for `[a_k]P`, one for `[b_k]Q`), each
scalar multiplication at most 512 group operations = 512 units, i.e.
`32 * 2 * 512 = 32768 = 2^15.0` units. (A P-256 scalar multiplication
needs at most 256 doublings plus ~128 additions, i.e. ~384 group
operations, so 512 is an upper bound.) Walk restarts: each restart
costs the same 2 scalar multiplications (1024 units); the expected
number of walks is `2^127.825748 / 2^32 = 2^95.825748`, so restarts cost
`2^105.825748` units in total. Final verification lump: `2^22`
instructions (`~2^12.034` units), covering the verification scalar
multiplication `[x]P` (at most 512 group operations = 512 units), the
modular inversion for collision resolution, encoding checks and output.
Distinguished-point sorting: `N * L * 448` instructions with
`N ~= 2^95.825748` stored points and `L = ceil(log2 N) = 96`, following
the secp256k1 baseline's merge-sort formula (Section 6). A hash table
is deliberately **not** used: its total work is only an expected value
(collision-dependent), while the sorting charge is a deterministic
upper bound given `N`.

## 6. Total work, and the score

Expected number of steps to a collision (Section 4):
`T = sqrt(pi*r/4)`, `log2(T) = 127.825748`. With per-step units
`1.8648` (Section 5),

```
W = T * 1.8648                                  (walk steps, dominant)
  + 2^15.0                                      (table build)
  + 2^105.825748                                (walk restarts)
  + 2^101.252282                                (sort distinguished points)
  + 2^12.034216                                (verification lump)
```

in target-group-operation equivalents, all work charged exactly once.
Numerically (log2, rounded to six decimals):

| Component | log2 (units) |
| --- | --- |
| Walk steps `127.825748 + log2(1.8648)` | **128.724769** |
| Sort distinguished points (`N*L*448/1000`) | 101.252282 |
| Walk restarts | 105.825748 |
| Table build (preprocessing) | 15.000000 |
| Verification lump | 12.034216 |
| **Total `log2(W)`** | **128.724769** |

The non-walk terms together add less than `2^-22` relative to the walk
term, so they do not move the total at six decimals. Direct check:
`2^(12872477/100000) >= W` (the true `log2(W)` including all small
terms is `128.7247691663...`, so the submitted rational
`12872477/100000` is a ceiling of it). The exact rational is already
reduced (`gcd(12872477, 100000) = 1`).

The simpler variant, for comparison:

| Variant | Steps (log2) | Per-step units | Total (log2) |
| --- | --- | --- | --- |
| Plain rho | 128.325748 | 1.712 | 129.101431 |
| **Negation rho (claimed)** | **127.825748** | **1.8648** | **128.724770** |

The negation map beats plain rho by `0.376662` bits: the `sqrt(2)`
iteration saving (0.5 bits) exceeds the canonicalisation plus
fruitless-cycle overhead (~0.12 bits).

## 7. Memory, preprocessing, advice

- **Peak memory** (upper bound): the distinguished-point store holds
  `N ~= 2^95.825748` records of 96 bytes (x-coordinate + `a` + `b`),
  i.e. `2^102.410711` bytes; merge sort needs a second buffer of the
  same size, so the store plus sort buffers are at most `2^103.410711`
  bytes, stated as `2^103.42` bytes. Walk state for `2^20` parallel
  walks (one point + two coefficients + partition state each, at most
  256 bytes per walk) is at most `2^28` bytes. The `r`-adding table is
  32 points (`2^11` bytes). Total peak memory is at most `2^103.42`
  bytes. Memory has no tie-break; it is stated for completeness.
- **Preprocessing**: the `r`-adding table build, `2^15.0` units,
  included in the total above. No preprocessing is reused across
  targets and no credit is taken for it. (The distinguished-point sort
  is online work during the walk, not preprocessing.)
- **Nonuniform advice**: `0` bytes. The `r`-adding table is built by
  the attack itself from the public inputs (`P`, `Q`) and its
  construction is charged in full above; it is not an embedded
  precomputed answer. No other advice is used.

## 8. Comparison with plain rho

For a group of prime order `r`, plain Pollard rho needs about
`sqrt(pi*r/2)` group operations, i.e. `log2 = 128.325748` here in
group-operation-only terms (`129.101431` with the honest per-step
charges of Section 5); with the negation map about `sqrt(pi*r/4)`,
i.e. `log2 = 127.825748` group-operation-only (`128.724770` honestly
charged). The textbook group-operation-only figures are lower than this
submission's bound because they charge only group operations, not the
per-step scalar work and not the P-256-specific addition price — the
honest charging of Sections 3.4-5 is what moves the bound up, as stated
in Section 1. On the honestly-charged comparison this submission's
`128.724770` beats plain rho by `0.376662` bits. There is no prior
accepted score on this track beyond the unscored template, so the
comparison is against plain rho itself, not against a live best.

## 9. Sensitivity: what if the charges are higher?

- **Standalone addition price (load-bearing, disclosed).** If every
  group addition were instead charged a standalone affine price of
  about 68.7 units (~68,678 instructions: a Fermat inversion of ~266
  field multiplications at the charged 256-instruction price, plus the
  2-multiplication formula and its additions, contradicting the cost
  model's own definition, which has no standalone inversion term),
  the per-step charge becomes `(68.678 + 0.240) * 1.05 = 72.3639`
  units and the total becomes `log2(W) = 134.002946`. The claim stands
  or falls on the cost model's batched-price convention plus the
  `2^20`-parallel-walk schedule of Section 3.4. The submitted bound
  `128.724770` uses the model's batched price.
- **`F_p` multiplication at twice the charged price.** If each field
  multiplication cost 512 instructions instead of the charged 256,
  the per-step charge becomes `(6*512 + 240) * 1.05 = 3477.6`
  instructions (`3.4776` units) and the total `2^129.623840`. The
  break-even points: the claim stays below `2^129` for field
  multiplications up to ~318 instructions, below plain rho's
  `2^129.101431` up to ~344 instructions, and below `2^130` up to
  ~676 instructions (about 2.6x the charged price).
- **Sort at twice the cost.** Doubling the merge-sort charge moves
  the sort term to `2^102.25` units, still negligible against the
  `2^128.72` walk term; the total is unchanged at six decimals.
- **Fruitless-cycle allowance.** Doubling the allowance from 5% to
  10% adds `log2(1.10/1.05) = 0.067` bits, giving `2^128.791883`.
  The margin is not balanced on the 5% figure.
- **Success probability.** At a `1.5x` step budget the success
  probability rises to `~0.829` for `+0.585` bits (`2^129.309732`);
  the claim is not balanced on the success bound either.

## 10. Heuristics, limitations, and what the experiment shows

- **Heuristic H1 (random-walk model):** the only heuristic in the
  claim (Section 4). The quotient walk is modelled as a random
  mapping on `r/2` classes. Support: the standard rho analysis and
  the empirical behaviour of `r`-adding walks with a 32-entry table.
  Limitation: it is not a theorem for this concrete iteration
  function; a reviewer demanding a proven mixing bound would reject.
  The fruitless-cycle handling (Section 3.3) is a second disclosed
  heuristic-adjacent ingredient: the 5% allowance is a charged bound
  on the loss, not a measured figure, and the literature estimate it
  rests on is roughly one to two percent.
- **Limitations.** The bound is an analytic bound on *expected* work,
  not a measured runtime; no full-size execution is claimed or
  possible. Memory is astronomical (`2^103.42` bytes); this is not a
  practical attack. The score is the log2 of expected work at the mean
  budget, so a single run could cost more (or less); the success
  probability `0.544` is the corresponding birthday figure, and the
  `1.5x`-budget sensitivity in Section 9 bounds the cost of buying a
  higher success probability.
- **Experiment.** `experiment.py` implements the exact walk logic of
  Section 3 (negation canonicalisation, `r`-adding table,
  distinguished points, cycle detection with restart, collision
  resolution) on the additive group `Z/qZ` for the declared inputs
  `(q, secret) = (10009, 1234)`, `(127, 42)`, and `(10009, 1)`. It
  recovers the secret in every run (272, 85, and 388 steps; 5, 9, and
  19 cycle restarts respectively; 14, 3, and 13 distinguished points
  stored) and also checks the class structure (orbits of size 2 except
  0). The toy uses a hash dict, not merge sort, and a 32-entry table
  with small distinguished threshold — it validates algorithm logic
  only. It does not establish production cost bounds, memory bounds,
  or P-256 arithmetic; those rest on the analytic accounting in
  Sections 3.4-6.
- **What could still be wrong.** The main risks are (a) heuristic H1
  failing for this concrete walk (Section 4), (b) the batched-price
  premise for `C = 1000` (Section 3.4 states the schedule and Section
  9 bounds the standalone-price consequence at `2^134.002946`), (c) a
  charge the reviewer judges too low (Section 9 bounds the
  consequence of the largest such charges; the field-multiplication
  break-even is ~676 instructions for the total to stay below
  `2^130`), and (d) an arithmetic slip in the Section 6 numbers (they
  were computed by script from the formulas in Section 5; the
  formulas, not just the numbers, are stated so they can be
  rechecked).

## References

- J. M. Pollard, "Monte Carlo methods for index computation (mod p)",
  Math. Comp. 32 (1978), pp. 918-924. The rho algorithm and the
  `sqrt(pi*n/2)` analysis.
- D. J. Bernstein, T. Lange and P. Schwabe, "On the correct use of the
  negation map in the Pollard rho method", PKC 2011, LNCS 6571,
  pp. 128-146. The negation map, fruitless cycles and countermeasures.
- P. C. van Oorschot and M. J. Wiener, "Parallel collision search
  with cryptanalytic applications", J. Cryptology 12 (1999),
  pp. 1-28. Distinguished points and parallel walks.
- A. Menezes, P. van Oorschot and S. Vanstone, Handbook of Applied
  Cryptography, Section 3.6 (generic DLP algorithms).
- The track's own `target-profile.json`,
  `classical-group-work-v2.json` and `tracks/docs/COST_MODEL.md` for
  the target, the borrowed reference price `C = 1000`, and the
  accounting rules applied in Sections 3.4-6. The values of `p`, `r`,
  `a`, `b` and `P` in Section 2 were verified by direct computation
  from the profile (primality of `p` and `r`, `a = p - 3`, `P` on the
  curve, `[r]P = O`).
