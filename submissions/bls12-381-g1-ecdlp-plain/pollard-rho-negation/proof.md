# bls12-381-g1-ecdlp-plain: Pollard rho with the negation map

Status: new submission draft; no qualification or emitted score. This is an
analytic cryptanalysis argument, not a claim of practical execution. The
normative inputs/outputs are those of the track's `target-profile.json` and
`classical-group-work-v2.json`. No restricted challenge distribution, advice,
random-oracle premise or quantum operation is used. One disclosed heuristic
is used (the random-walk / random-mapping model for the quotient walk,
Section 4); everything else is exact arithmetic.

## 1. The claim

Let `p` be the BLS12-381 field prime, `r` the (prime) group order of
`G1 = <P0>`, `P0` the standard generator from the target profile, and
`x` uniform in `[0, r)`. The input is `T = [x]P0`; the output is `x`.
The Pollard rho algorithm below, with the negation map (Section 3),
recovers `x` with probability at least `0.544` (hence `>= 0.39`) at the
stated mean step budget, with **expected** total work

- `log2(W) <= 127.488182` (exact rational score `63744091/500000`),
- preprocessing `<= 2^15.0` target-group-operation equivalents (included),
- peak memory `<= 2^103.07` bytes,
- nonuniform advice `0` bytes.

All logarithms below are rounded **up** to six decimal places; the exact
arithmetic in Sections 5-6 is the audit authority. The reference price is
`C = 2000` (batched incomplete affine G1 addition, 6M convention). The
iteration count is an expectation under the disclosed random-walk
heuristic (Section 4); the score is the log2 of the expected total work
at the mean budget, which is the standard accounting for Pollard rho.

**A correction to a naive reading, stated up front.** The textbook
statement of Pollard rho with the negation map is
`sqrt(pi*r/4) ~= 2^127.25` *group operations*. That figure counts only
the group operations. Under this challenge's cost model *all* work must
be charged once, including the per-step scalar arithmetic
(coefficient updates, canonicalisation), the fruitless-cycle
allowance, the distinguished-point sorting, and the walk restarts. When
all of that is charged honestly (Sections 5-6), the bound is
`2^127.488182`, not `2^127.25`. This submission claims the honest
`2^127.488182`. It still beats honestly-charged plain rho
(`2^127.875971`, Section 8) by about 0.39 bits. This track is unscored
(no qualified baseline); if accepted, this claim sets the first score.

## 2. Number-theoretic facts (all verified by direct computation)

```
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
  = 2^381 - 2^126 + ... (see target-profile.json; log2(p) = 380.700672)
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
P0 = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
      0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
```

The curve is `y^2 = x^3 + 4` over `F_p` (`j = 0`). `G1 = <P0>` has
prime order `r`. `log2(r) = 254.857089`, so `log2(r)` is **not** 256:
every `sqrt(pi*r/(2m))` figure below uses the exact `r` (the script
that produced the Section 6 numbers used the exact integer `r`, not
`2^255`). The cofactor of `G1` in `E(F_p)` is
`h = 0x396c8c005555e1568c00aaab0000aaab`; the walk of Section 3 stays
inside `<P0>` because every table point is a combination of `P0` and
`T = [x]P0`, so the cofactor never enters.

Verified by direct computation (py_ecc, pinned commit in the target
profile):

- `P0` is on the curve `y^2 = x^3 + 4` over `F_p`;
- `[r]P0 = O` (the identity), so `P0` has exact order `r`;
- `-P0 = (x, -y)` is distinct from `P0` and on the curve; `P0 + (-P0)
  = O`. Since `r` is odd, no point of `<P0>` satisfies `R = -R`
  except `O`.

**The GLV endomorphism exists — and is deliberately not claimed.**
The curve has `j = 0` and `p = 1 (mod 3)`, so `F_p` contains a
primitive cube root of unity

```
zeta = 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe
```

with `zeta^3 = 1 (mod p)`, `zeta != 1`. The map `phi(x, y) =
(zeta*x, y)` is an endomorphism of `E` over `F_p`, and on `<P0>` it
acts as scalar multiplication by `lambda^2`, where

```
lambda = 0xac45a4010001a40200000000ffffffff
```

satisfies `lambda^2 + lambda + 1 = 0 (mod r)`: verified by direct
computation, `phi(P0) = [lambda^2]P0` and `phi^2(P0) = [lambda]P0`
(affine coordinates compared exactly). Section 8 shows why this
endomorphism buys essentially nothing here: its honestly-charged
total is `2^127.483789`, only `0.004393` bits better than the
negation-only claim — within the estimation noise of the per-step
charges — so it is disclosed as a sensitivity, not claimed. A
reviewer who knows the BLS12-381 G1 endomorphism is therefore not
catching an omission; the omission is deliberate and quantified.

## 3. The algorithm

Notation: `[a]P` is scalar multiplication of a point `P` by an integer
`a` (taken mod `r`). Group operations are BLS12-381 G1
additions/doublings. All field arithmetic is in `F_p`; all scalar
arithmetic is mod `r`.

### 3.1 Distinguished-point rho, parallel walks

Fix a 32-entry `r`-adding table: 32 random pairs `(a_k, b_k)` in
`[0, r)` and the points `T_k = [a_k]P0 + [b_k]T`. A walk state is
`(R, a, b)` with the invariant `R = [a]P0 + [b]T`. A step from a state
whose (canonicalised) point is `R'` with coefficients `(a', b')` is:
`k = partition(R')` (low 5 bits of the x-coordinate of `R'`),
`R <- R' + T_k`, `a <- a' + a_k (mod r)`, `b <- b' + b_k (mod r)`.

A point is *distinguished* if the low 32 bits of its x-coordinate are
zero (fraction `2^-32` of points). Walks run until they reach a
distinguished point, which is stored (with its `(a, b)`); the walk then
restarts from a fresh random start. `2^20` walks run in parallel in
rounds: in each round every active walk takes one step, and the `2^20`
independent additions of the round are batched (Section 5). This is the
standard van Oorschot-Wiener parallelisation.

A collision is two stored distinguished points with the same `R` but
different `(a, b)`: then `[a1]P0 + [b1]T = [a2]P0 + [b2]T`, i.e.
`(a1 - a2) = (b2 - b1)*x (mod r)`. If `b2 != b1 (mod r)` (the contrary
case is a trivial collision with `a1 = a2`, see Section 4), then
`x = (a1 - a2) * (b2 - b1)^-1 mod r`. The candidate is verified by
computing `[x]P0` and accepting only if it equals `T`.

Edge case: `T = O` (i.e. `x = 0`) is detected by a free identity check
before the walk starts. For `x != 0` the analysis below applies.

### 3.2 Negation map quotient walk

Before each step, the current point `R` is replaced by the canonical
representative of its class `{R, -R}`: the class element whose
y-coordinate is smaller as an integer in `[0, p)` (negation flips `y`
only, so `x` is unchanged; the two y-values are `y` and `p - y`, and
they are equal only for `y = 0`, which never occurs for `R != O`
since `r` is odd). The coefficients transform the same way: if the
representative is `-R`, then `(a', b') = (-a mod r, -b mod r)`.
The step of Section 3.1 is then taken from `(R', a', b')`.

The walk therefore lives on the quotient by the 2-element classes: the
effective search space has size `r/2`, and the expected number of steps
to a collision is `sqrt(pi*r/4)` (Section 4). The classes have size
exactly 2 for every `R != O` (Section 2), so no class is over- or
under-counted; the only fixed point of the quotient is `O`, which a
walk reaches with probability `2^-255` per step and which is handled by
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

### 3.4 Batching schedule for the `C = 2000` price

The cost model defines `C = 2000` as the price of a *batched*
incomplete affine G1 addition under the 6M convention, with no
standalone inversion charged per addition; this proof charges every
group addition/doubling at exactly that price, as the model
prescribes. The batching is concrete: in each round the `2^20`
independent additions (one per active walk) have all their
denominators available up front for one batch inversion, exactly as
the cost model defines `C`. The number of rounds is
`2^127.254293 / 2^20 = 2^107.254293`; walks that reach a distinguished
point are replaced by fresh walks at round boundaries. This is a
material premise of the submitted bound: it relies on the cost
model's batched-price convention. If a judge instead required a
standalone affine-addition price (about 75 units per addition, i.e.
~150,000 instructions for a Fermat inversion, see Section 9), the
total becomes `2^133.555807` — worse than honestly-charged plain rho
(Section 9). The submitted bound `127.488182` uses the model's
batched price.

## 4. Correctness and success probability

**Heuristic H1 (random-walk / random-mapping model, disclosed).** The
quotient walk of Section 3.2 is modelled as a random mapping on a set
of `r/2` classes: each step lands on a fresh class uniformly at random,
independently of the past. This is the standard heuristic for Pollard
rho (Pollard 1978; Bernstein-Lange-Schwabe 2011 for the negation map).
It is *not* a theorem for this concrete iteration function: the
`r`-adding table is fixed, the partition function is deterministic,
and the quotient structure could in principle bias the walk. The
heuristic is supported by the usual evidence: `r`-adding walks with a
32-entry table are empirically indistinguishable from random mappings
at cryptographic sizes, and the same heuristic underlies every
published rho complexity estimate, including the textbook
`sqrt(pi*r/2)` figure. A reviewer who rejects H1 rejects the claim;
the claim does not pretend otherwise.

Under H1, the number of steps `T` to the first repeated class is the
birthday collision time for a set of size `N = r/2`: `E[T] =
sqrt(pi*N/2) = sqrt(pi*r/4)`, i.e. `log2 = 127.254293`. For a budget of
`T` steps the collision probability is `1 - exp(-T^2/(2N))`; at the mean
budget `T = sqrt(pi*N/2)` this is `1 - exp(-pi/4) ~= 0.544 >= 0.39`.
At a `1.5x` budget the probability is `1 - exp(-(pi/4)*2.25) ~= 0.829`
for `+log2(1.5) = 0.585` bits (`2^128.073143`, still beating
honestly-charged plain rho), so the claim is not balanced on the
success bound.

**Collision resolution.** Given a genuine collision (same `R`, different
`(a, b)`), the recovered `x` is correct by the algebra of Section 3.1,
and the verification `[x]P0 = T` confirms it. A *trivial* collision
(`b1 = b2 (mod r)`, hence also `a1 = a2`) carries no information; the
coefficients at a collision are effectively uniform in `[0, r)`, so
this happens with probability `~2^-255` per collision, and the walk
simply continues (charged inside the 5% allowance). The stored
distinguished points are compared exactly (112-byte records), so no
false collision is ever accepted: a wrong `x` fails the verification
and the walk continues.

**Deterministic parts.** Everything except H1 is exact: the class
structure (Section 2), the coefficient invariant `R = [a]P0 + [b]T`
(preserved by construction through canonicalisation and table steps),
the distinguished-point fraction `2^-32`, and the final verification.
Success probability `>= 0.544` is therefore the birthday probability
under H1, not a bound on an unmodelled failure mode.

## 5. Per-step cost charges (instructions, C = 2000)

Work is counted in scalar AArch64 instruction equivalents, then divided
by the fixed `C = 2000` to give group-operation equivalents, exactly as
the cost model prescribes. One group addition/doubling is one reference
operation, i.e. 2000 instructions. All other work is charged at
justified instruction counts; nothing is charged twice (a reference
operation is never also charged through its internal instructions).

| Work item | Charge | Justification |
| --- | --- | --- |
| Group addition / doubling | 2000 instructions (1 unit) | Reference operation, fixed `C` (batched incomplete affine G1 addition, 6M convention; Section 3.4). |
| Partition index (low 5 bits of x) | 16 instructions | A few loads, a mask and a table-index computation. |
| Coefficient updates `a, b mod r` (2 x 64) | 128 instructions | Two 255-bit modular additions: 4 limb adds plus a conditional subtraction each; 64 instructions each is a generous upper bound. |
| Distinguished-point check (low 32 bits of x zero) | 32 instructions | Load the low limb of x, compare against zero, branch; includes the short-history comparison for cycle detection (Section 3.3). |
| Negation canonicalisation | 64 instructions | Compare y against `p - y` (6-limb compare on the 381-bit field) and conditionally negate; a few field operations. |
| **Per-step instructions (negation rho)** | **2240** | Sum of the rows above. |
| Fruitless-cycle allowance (x1.05) | x1.05 | Covers aborted walks and restarts from Section 3.3; the literature puts the expected loss at roughly one to two percent of steps, so 5% is generous headroom. |
| **Per-step units (negation rho)** | **1.176** | `2240 * 1.05 / 2000`. |

For comparison, the per-step charge of plain rho (no canonicalisation,
no cycle allowance) is `2176` instructions (`1.088` units).

**No field multiplication is charged per step.** The claimed variant
uses no `F_p` multiplication outside the reference operation itself:
canonicalisation is a compare plus a conditional negation, the
coefficient updates are modular additions, and the table build and
verification use group operations. The organizer's pinned `F_p`
multiplication price (383 instructions, blst) therefore does not enter
the claimed bound at all; it enters only the unclaimed GLV sensitivity
of Section 8. This makes the claimed bound insensitive to the field
multiplication price.

**Other charged work.** The `r`-adding table build: 32 entries x 2
scalar multiplications each (one for `[a_k]P0`, one for `[b_k]T`), each
scalar multiplication at most 512 group operations = 512 units (a
255-bit double-and-add needs at most 254 doublings + 127 additions =
381 operations, so 512 is headroom), i.e. `32 * 2 * 512 = 32768 =
2^15.0` units. Walk restarts: each restart costs the same 2 scalar
multiplications (1024 units); the expected number of walks is
`2^127.254293 / 2^32 = 2^95.254293`, so restarts cost `2^105.254293`
units in total. Final verification lump: `2^22` instructions
(`2^11.03` units), covering the verification scalar multiplication
`[x]P0` (at most 512 group operations = 512 units), the modular
inversion for collision resolution, encoding checks and output.
Distinguished-point sorting: `N * L * 448` instructions with
`N ~= 2^95.254293` stored points and `L = ceil(log2 N) = 96`, following
the secp256k1 baseline's merge-sort formula (Section 6). A hash table
is deliberately **not** used: its total work is only an expected value
(collision-dependent), while the sorting charge is a deterministic
upper bound given `N`.

## 6. Total work, and the score

Expected number of steps to a collision (Section 4):
`T = sqrt(pi*r/4)`, `log2(T) = 127.254293`. With per-step units
`1.176` (Section 5),

```
W = T * 1.176                                     (walk steps, dominant)
  + 2^15.0                                        (table build)
  + 2^105.254293                                  (walk restarts)
  + 2^99.680826                                   (sort distinguished points)
  + 2^11.03                                       (verification lump)
```

in target-group-operation equivalents, all work charged exactly once.
Numerically (log2, rounded up to six decimals):

| Component | log2 (units) |
| --- | --- |
| Walk steps `127.254293 + log2(1.176)` | **127.488181** |
| Sort distinguished points (`N*L*448/2000`) | 99.680826 |
| Walk restarts | 105.254293 |
| Table build (preprocessing) | 15.000000 |
| Verification lump | 11.034216 |
| **Total `log2(W)`** | **127.488182** |

The non-walk terms together add less than `2^-22` relative to the walk
term, so they do not move the total at six decimals except for the
ceiling: the true `log2(W)` including all small terms is
`127.4881811300...`, so the submitted rational `63744091/500000 =
127.488182` is a ceiling of it (verified: `2^(63744091/500000) >= W`).
The exact rational is already reduced (`gcd(63744091, 500000) = 1`).

The plain-rho variant, for comparison:

| Variant | Steps (log2) | Per-step units | Total (log2) |
| --- | --- | --- | --- |
| Plain rho | 127.754293 | 1.088 | 127.875971 |
| **Negation rho (claimed)** | **127.254293** | **1.176** | **127.488182** |

The negation map buys `0.387790` bits over plain rho on the
honestly-charged comparison.

## 7. Memory, preprocessing, advice

- **Peak memory** (upper bound): the distinguished-point store holds
  `N ~= 2^95.254293` records of 112 bytes (48-byte x-coordinate +
  32-byte `a` + 32-byte `b`), i.e. `2^102.061648` bytes; merge sort
  needs a second buffer of the same size, so the store plus sort
  buffers are at most `2^103.061648` bytes, stated as `2^103.07`
  bytes. Walk state for `2^20` parallel walks (one point + two
  coefficients + partition state each, at most 256 bytes per walk) is
  at most `2^28` bytes. The `r`-adding table is 32 points (`2^11`
  bytes). Total peak memory is at most `2^103.07` bytes. Memory has no
  tie-break; it is stated for completeness.
- **Preprocessing**: the `r`-adding table build, `2^15.0` units,
  included in the total above. No preprocessing is reused across
  targets and no credit is taken for it. (The distinguished-point sort
  is online work during the walk, not preprocessing.)
- **Nonuniform advice**: `0` bytes. The `r`-adding table is built by
  the attack itself from the public inputs (`P0`, `T`) and its
  construction is charged in full above; it is not an embedded
  precomputed answer. No other advice is used.

## 8. Comparison with plain rho, and the GLV sensitivity

For a group of prime order `r`, plain Pollard rho needs about
`sqrt(pi*r/2)` group operations, i.e. `log2 = 127.754293` here
(`127.875971` with the honest per-step charges of Section 5); with
the negation map about `sqrt(pi*r/4)`, i.e. `log2 = 127.254293`
(`127.488182` honestly charged). The textbook group-operation-only
figures are lower than this submission's bound because they charge
only group operations, not the per-step scalar work — the honest
charging of Section 5 is what moves the bound up, as stated in
Section 1. On the honestly-charged comparison this submission's
`127.488182` beats plain rho by `0.387790` bits. There is no prior
accepted score on this track (the selected package is an unscored
template), so this is a first-score claim, not a comparison against a
live best.

**The GLV endomorphism (disclosed, not claimed).** Section 2 shows
`phi(x, y) = (zeta*x, y)` is an endomorphism of `E` over `F_p` with
`phi(P0) = [lambda^2]P0`. A GLV quotient walk on the 6-element classes
`{+-phi^j(R)}` would need `sqrt(pi*r/12)` steps, i.e. `log2 =
126.461812`, but each step would additionally cost two `F_p`
multiplications for `phi`/`phi^2` (2 x 383 instructions at the
organizer's pinned price), a 96-instruction representative selection,
and two `lambda^j` coefficient multiplications (2 x 383): per-step
`3868` instructions x1.05 / 2000 = `2.0307` units, for a total of
`2^127.483789`. That is only `0.004393` bits better than the claimed
`2^127.488182` — the `sqrt(3)` iteration saving (0.792 bits) is almost
entirely eaten by the per-step overhead (0.788 bits), because the BLS
field multiplication (383 instructions) is more than twice the
secp256k1 one (179). The gain is smaller than the estimation noise in
the per-step charges (for example, the 96-instruction representative
selection alone is worth 0.036 bits at this scale), so this submission
does **not** claim the GLV variant. If a reviewer prefers the GLV
accounting, the bound would be `127.483789`; nothing in the claimed
bound depends on ignoring the endomorphism.

## 9. Sensitivity: what if the charges are higher?

- **Standalone addition price (load-bearing, disclosed).** If every
  group addition were instead charged a standalone affine price of
  about 75 units (~150,000 instructions — a Fermat inversion per
  addition, contradicting the cost model's own definition, which has
  no standalone inversion term), the per-step charge becomes
  `(75 + 240/2000) * 1.05 = 78.8760` units and the total becomes
  `log2(W) = 133.555807` — **worse** than honestly-charged plain rho.
  The 75-unit figure is itself conservative: a Fermat inversion in
  `F_p` needs about 380 field multiplications, i.e. `380 * 383 =
  145,540` instructions `= 72.77` units, below the 75 charged here.
  The claim stands or falls on the cost model's batched 6M convention
  plus the `2^20`-parallel-walk schedule of Section 3.4. The submitted
  bound `127.488182` uses the model's batched price.
- **Fruitless-cycle allowance.** Doubling the allowance from 5% to
  10% adds `log2(1.10/1.05) = 0.067` bits, giving `2^127.555295` —
  still beating plain rho. The margin is not balanced on the 5% figure.
- **Sort at twice the cost.** Doubling the merge-sort charge moves
  the sort term to `2^100.68` units, still negligible against the
  `2^127.49` walk term; the total is unchanged at six decimals.
- **Success probability.** At a `1.5x` step budget the success
  probability rises to `~0.829` for `+0.585` bits (`2^128.073143`),
  still beating plain rho; the claim is not balanced on the success
  bound either.
- **Field-multiplication price.** The claimed bound charges no `F_p`
  multiplication per step (Section 5), so it is unaffected by the
  field-multiplication price. The unclaimed GLV variant of Section 8
  does depend on the pinned 383-instruction price; that is one more
  reason it is disclosed rather than claimed.

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
  possible. Memory is astronomical (`2^103.07` bytes), as with the
  BSGS controls on other tracks (`2^136` bytes on secp256k1); this is
  not a practical attack. The score is the log2 of expected work at
  the mean budget, so a single run could cost more (or less); the
  success probability `0.544` is the corresponding birthday figure,
  and the `1.5x`-budget sensitivity in Section 9 bounds the cost of
  buying a higher success probability.
- **Experiment.** `experiment.py` implements the exact walk logic of
  Section 3 (negation-map canonicalisation, `r`-adding table,
  distinguished points, cycle detection with restart, collision
  resolution) on the additive group `Z/qZ` for the declared input
  `(q, secret) = (10009, 1234)`, where `q` is prime. It recovers the
  secret in every run and also checks the class structure (orbits of
  size 2 except 0). The toy uses a hash dict, not merge sort, and a
  32-entry table with small distinguished threshold — it validates
  algorithm logic only. It does not establish production cost bounds,
  memory bounds, or BLS12-381 arithmetic; those rest on the analytic
  accounting in Sections 5-6.
- **What could still be wrong.** The main risks are (a) heuristic H1
  failing for this concrete walk (Section 4), (b) the batched-price
  premise for `C = 2000` (Section 3.4 states the schedule and Section
  9 bounds the standalone-price consequence at `2^133.555807`), (c) a
  charge the reviewer judges too low (Section 9 bounds the
  consequence of the largest such charges; the claimed bound charges
  no field multiplication per step, so that price cannot move it),
  and (d) an arithmetic slip in the Section 6 numbers (they were
  computed by script from the formulas in Section 5; the formulas,
  not just the numbers, are stated so they can be rechecked).

## References

- J. M. Pollard, "Monte Carlo methods for index computation (mod p)",
  Math. Comp. 32 (1978), pp. 918-924. The rho algorithm and the
  `sqrt(pi*n/2)` analysis.
- D. J. Bernstein, T. Lange and P. Schwabe, "On the correct use of the
  negation map in the Pollard rho method", PKC 2011, LNCS 6571,
  pp. 128-146. The negation map, fruitless cycles and countermeasures.
- R. P. Gallant, R. J. Lambert and S. A. Vanstone, "Faster point
  multiplication on elliptic curves with efficient endomorphisms",
  CRYPTO 2001, LNCS 2139, pp. 190-200. The GLV endomorphism (discussed
  in Section 8; not claimed here).
- S. D. Galbraith and R. S. Ruprai, "Using equivalence classes to
  accelerate solving the discrete logarithm problem in a short
  interval", PKC 2010, LNCS 6056, pp. 368-383. Quotient walks.
- P. C. van Oorschot and M. J. Wiener, "Parallel collision search
  with cryptanalytic applications", J. Cryptology 12 (1999),
  pp. 1-28. Distinguished points and parallel walks.
- R. Wahby and D. Boneh, "Fast and simple constant-time hashing to
  the BLS12-381 elliptic curve", ePrint 2019/403. The BLS12-381
  parameters (field prime, group order, generator) used in Section 2.
- A. Menezes, P. van Oorschot and S. Vanstone, Handbook of Applied
  Cryptography, Section 3.6 (generic DLP algorithms).
- The track's own `target-profile.json`,
  `classical-group-work-v2.json` and `tracks/docs/COST_MODEL.md` for
  the target, the reference price `C = 2000`, and the accounting rules
  applied in Sections 5-6. The values of `zeta` and `lambda` in
  Section 2 were verified by direct computation from the profile's
  `p` and `r` (including `phi(P0) = [lambda^2]P0` by direct
  elliptic-curve arithmetic with the pinned py_ecc commit).
