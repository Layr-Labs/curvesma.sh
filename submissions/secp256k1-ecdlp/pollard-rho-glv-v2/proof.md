# secp256k1-ecdlp: Pollard rho over GLV automorphism equivalence classes (v2)

Status: candidate for independent review; no qualification or emitted score.
This is an analytic attack argument in the fixed `classical-group-work-v2`
accounting model, not an executed production attack. `claim.json` binds the
`production-v2` target profile and the `classical-group-work-v2` cost model
by SHA-256. Total claimed work is **2^127.22 target-group-operation
equivalents** with success probability **at least 0.405** over the uniform
challenge distribution and fresh algorithm coins. The only heuristic input
is the standard random-mapping heuristic H1 for the Pollard rho iteration
(Section 5), disclosed with its support and scope; every other step is exact.

Revision note. This v2 candidate keeps the v1 candidate's algorithm,
parameters (except the step budget), walk analysis and ledger structure, and
revises its accounting in exactly two places, both exact: (i) the
canonicalisation no longer prices a second field multiplication — since
`beta^2 = -1 - beta` in `F_p`, the value `beta^2*x` is derived from
`beta*x` by one field addition and one negation (Section 7); (ii) the
success exponent is recomputed for a claimed probability of 0.405 instead of
0.43, still clearing the 0.39 threshold with an absolute margin above 0.015
after all deductions (Section 5). The remaining `beta` multiplication is
still charged at the organiser's full measured price of 179 instructions,
and the scalar-update and miscellaneous per-step lines are unchanged.

## 1. Target and claim

Let `p` be the secp256k1 field prime, `E: y^2 = x^3 + 7` over `F_p`, `P` the
standard generator from the protected profile, and

```
r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
  = 115792089237316195423570985008687907852837564279074904382605163141518161494337
```

the (prime) group order. Input: `Q = [x]P` with `x` uniform in `[0, r)`.
Output: `x`. Claim (first declared claim supplies the score):

* total classical work `<= 2^127.22` group-operation equivalents, all charges
  included (ledger in Section 7);
* success probability `>= 0.405 >= 0.39` within that budget;
* preprocessing `<= 2^27` equivalents (included in the total);
* peak memory `<= 2^113` bytes; nonuniform advice `2^0` bytes (none; the
  adding-point table is built per instance from the public input and its
  construction is charged as preprocessing).

The generic BSGS accounting control (129.000802) is beaten by 1.78 bits.
Section 9 explains the margin: the step count gains `sqrt(6)` from the
equivalence classes; part of the gain is spent honestly on per-step
canonicalisation at the organiser's measured field-multiplication price and
on a success exponent chosen to clear the 0.39 threshold with margin.

## 2. Automorphism classes

secp256k1 has the GLV endomorphism `phi(x, y) = (beta*x, y)`, where `beta`
is a primitive cube root of unity in `F_p` (`beta != 1`, `beta^3 = 1`), and
`phi(P) = [lambda]P` with `lambda^2 + lambda + 1 = 0 (mod r)`. Together with
negation this gives an automorphism group `mu = {+-1, +-lambda, +-lambda^2}`
of order 6 acting on `<P>`. Every orbit ("class") of a non-identity point
has exactly 6 elements; the identity forms a class of its own. Since
`r = 1 (mod 6)`, the number of classes is exactly

```
N'' = 1 + (r - 1)/6
    = 19298681539552699237261830834781317975472927379845817397100860523586360249057
    ~ 2^253.415.
```

The canonical representative of a class is the lexicographically smallest
`(x, y)` among its six elements (smallest `x`; on ties, `y < p/2`). Write
`canon(R)` for the representative and `tau(R) in mu` for the scalar with
`canon(R) = [tau(R)]R`. Both are computed from `R` alone (Section 7 prices
this), so any class-invariant function of the representative is a
well-defined function on classes.

Two structural facts drive the cheap canonicalisation of Section 7:

* `phi` fixes the `y`-coordinate, so the three points `[lambda^k]R`,
  `k = 0, 1, 2`, all share `R`'s `y`-coordinate; the class's six
  `x`-coordinates are `x, x, beta*x, beta*x, beta^2*x, beta^2*x`. Choosing
  the canonical element therefore needs the minimum of the three values
  `x, beta*x, beta^2*x`, after which the sign is decided by `R`'s own `y`
  alone (pick `y` or `-y`, whichever is `< p/2`).
* `beta` is a root of `X^2 + X + 1` in `F_p`, i.e. `beta^2 = -1 - beta`
  exactly. Hence `beta^2 * x = -(x + beta*x) (mod p)`: the third candidate
  costs one field addition and one negation once `beta*x` is known, not a
  second multiplication.

## 3. The walk

Per instance (charged as preprocessing, Section 7): choose `J = 2^16`
adding points `S_j = [u_j]P + [v_j]Q`, `j = 0..J-1`, with public random
`u_j, v_j in [0, r)`. The partition index of a class is
`j(class) = (canonical x-coordinate) mod J`, a class invariant.

One step of the walk on classes: from the canonical representative `R` of
the current class, compute `T = R + S_{j(R)}` (one reference addition),
then output `canon(T)` as the next class representative. This is a function
on classes: two representatives of the same class have the same canonical
representative, hence the same `j`, the same `T` up to the automorphism, and
the same next class. The walk therefore descends to a single deterministic
iteration on the `N''` classes, the setting of heuristic H1 (Section 5).

Scalar bookkeeping. The walk maintains `(A, B)` with
`current point = [A]P + [B]Q`. Stepping adds `(u_j, v_j)`; canonicalisation
multiplies `(A, B)` by `tau`. Scalars are stored in the `lambda`-split form
`A = A1 + A2*lambda`, `B = B1 + B2*lambda (mod r)`. In split coordinates,
multiplication by `lambda` is `(A1, A2) -> (-A2, A1 - A2)` — additions only,
because `lambda^2 = -1 - lambda`. Reductions mod `r` are lazy (values are
kept below `2^258` and reduced when a trail reports).

Trails (van Oorschot–Wiener form). A trail starts at a fresh random point
`[a0]P + [b0]Q`, canonicalises, and walks until the first class whose
canonical x-coordinate satisfies `x = 0 (mod 2^24)` (distinguished, fraction
`theta = 2^-24`), or until `L_max = 2^30` steps (abandoned, Section 6). A
trail that reports stores `(canonical point, A, B)` in a table capped at
`2^104` entries. The algorithm stops when two stored reports coincide in
the canonical point with different `(A, B)`, or when the total step budget
`t` (Section 5) is exhausted (failure). The number of trail starts is capped
at `2^104`; reaching any cap is a failure event whose probability is
absorbed in Section 5.

Collision solving. Two reports of one class give
`[A1]P + [B1]Q = [A2]P + [B2]Q` as points (the canonical representative is
one fixed point of the class, so equality of representatives is equality of
points). Hence `(A1 - A2) + (B1 - B2) x = 0 (mod r)`. If `B1 != B2`,
`x = (A2 - A1)(B1 - B2)^{-1} mod r`, verified by checking `[x]P = Q` before
output (verification charged, Section 7). The case `B1 = B2` forces
`A1 = A2` (same point, same `B` determines `A` mod `r`), i.e. two identical
records; its probability is bounded in Section 5 and it is treated as a
non-event (the trail is discarded, the walk continues).

Exceptional steps. The reference addition is the batched incomplete form
and needs `R != +-S_j`. If `T in {infinity, R, -R}` (detected from the
addition's exceptional branch), the step is retried with
`S_{(j+1) mod J}`, then `S_{(j+2) mod J}`; if all three fail (possible only
if three distinct table points each coincide up to sign with `R`), the
trail is abandoned. Retry and abandonment costs are bounded worst-case in
Section 7.

## 4. Parameters

```
J        = 2^16            adding points
theta    = 2^-24           distinguished fraction
L_max    = 2^30 = 64/theta trail length cap
t        = 141561145530318902406949923548642162214   total step budget
         ~ 2^126.735
table cap    = 2^104 entries
starts cap   = 2^104
```

`t` is the least integer with `t^2 >= 2 * N'' * 0.519195` (exact integer
inequality: `t^2 * 10^6 >= 2 * N'' * 519195`, and `(t-1)` fails it). The
constant `0.519195` is `-ln(1 - 0.405) = 0.5191936...` rounded up, so the
birthday exponent of Section 5 delivers probability at least 0.405 with the
margin computed there.

## 5. Success probability

Heuristic H1 (random mapping). The iteration of Section 3, as a function on
the `N''` classes, behaves as a uniform random mapping for the purpose of
collision statistics: the classes visited by successive trails, before any
class repeats, are distributed as independent uniform samples from the
`N''` classes. This is the standard heuristic underlying all Pollard rho
security estimates (Pollard 1978; van Oorschot–Wiener 1999; the
automorphism-class variant is the standard treatment of rho with the
negation map and with GLV equivalence classes, e.g. Wiener–Zuccherato and
subsequent analyses of rho on GLV curves). It is supported at toy scale by
the included experiments (Section 8), which implement exactly the walk of
Section 3 — canonicalisation (in its v2 form, cross-checked against brute
force), split-coordinate bookkeeping, distinguished trails — and measure
visited-class collision statistics consistent with the birthday model below
(mean visited classes at collision 267.3 observed vs 262.2 predicted;
success rates in Section 8). H1 is the only heuristic in the claim.

Birthday bound. Under H1, after `t` steps the visited classes are (up to
the within-trail distinctness discussed below) `t` samples, and

```
P(no repeated class) <= prod_{i=0}^{t-1} (1 - i/N'')
                     <= exp(-t(t-1) / (2 N'')).
```

With `t^2 >= 2 N'' * 0.519195` and `t/(2N'') < 2^-127`, the exponent obeys
`t(t-1)/(2N'') >= 0.519195 - 2^-127 >= 0.519194`, so

```
P(collision within budget) >= 1 - e^-0.519194 >= 0.4050.
```

(The exact evaluation at the integer `t` of Section 4 gives
`>= 0.4050007`.) A repeated class yields a stored distinguished report
collision as follows: once two trails visit one class they follow the same
deterministic walk and report the same distinguished point, unless the
budget ends inside the residual gap to that report (expected `1/theta`
steps; the standard van Oorschot–Wiener accounting charges this overshoot
inside the step budget, and the exponent above already uses the raw step
count `t`, so no additional factor is due: a collision occurring by step
`t - L` for any realised residual gap `L` is detected by step `t`; the
probability that the first repeated class occurs in the final `2^25` steps
is `< 2^-101` of the total mass and is absorbed below).

Deductions (all under H1, total `< 2^-83`):

* Degenerate report collision `B1 = B2`: needs two distinct trails to
  accumulate identical split scalars at one point; probability
  `<= (#trail pairs) * 2^-256 < 2^-104`.
* Starts cap `2^104` reached before `t` steps: starts needed concentrate at
  `t * theta = 2^102.74` with geometric tails; overflow probability
  `< exp(-2^100)`.
* Table cap `2^104`: distinguished reports likewise concentrate at
  `2^102.74`; overflow `< exp(-2^100)`.
* A trail is lost to a distinguished-free cycle with probability `< 2^-187`
  per trail (Section 6), hence `< 2^-83` over at most `2^104` trails;
  within-trail class repeats short of a cycle waste `< 2^-50` steps in
  total expectation and do not change the exponent at the stated precision.

Hence

```
P(success within the stated budget) >= 0.4050007 - 2^-83 >= 0.405,
```

as claimed. The margin over the 0.39 threshold is therefore at least
`0.015` in absolute probability after every deduction. This is a bound
within the stated worst-case step budget `t`, not an expected-runtime
restatement: the algorithm halts at `t` steps and the probability is over
the challenge `x` and all algorithm coins.

## 6. Fruitless cycles and trail discipline

A deterministic iteration has cyclic classes; a trail that enters a cycle
containing no distinguished class would loop forever. Three facts dispose
of this, with the algorithm's `L_max` cap as the enforcement:

1. Short cycles. Under H1 the expected number of cycles of any given length
   `ell << sqrt(N'')` is `O(1/ell)`, so the expected number of classes lying
   on cycles of length `<= 2^30` is `O(2^30)`, a fraction `<= 2^-222` of all
   classes. A trail of at most `2^30` steps lands on such a class with
   probability `<= 2^30 * 2^-222 = 2^-192`.
2. Long cycles. A trail enters any cyclic class with probability at most
   `(trail length) * (cyclic fraction) <= 2^30 * 2^-126.4 < 2^-96` per
   trail (the cyclic fraction of a random mapping on `N''` classes is
   `sqrt(pi/(2N'')) < 2^-126.4`). On a cycle of length `> 2^30`, any
   `2^30` consecutive classes contain no distinguished class with
   probability at most `(1 - theta)^{2^30} <= e^-64 < 2^-92`, so an entered
   long cycle still produces a report within `L_max` steps except with
   probability `< 2^-92`. Combined per-trail loss `< 2^-187`; over at most
   `2^104` trails, `< 2^-83`.
3. Independently of probability, `L_max = 2^30 = 64/theta` caps every
   trail. A trail that exceeds `L_max` without reporting is abandoned and
   charged inside the budget `t`: under H1 the expected fraction of steps
   spent in over-long trails is `<= (1-theta)^{L_max} * (L_max * theta)
   <= e^-64 * 64 < 2^-86` of all steps, and such steps still visit classes
   and still contribute to the collision process until the cap.

The toy experiments (Section 8) exhibit this trapping at toy scale
(`N'' ~ 2^15.4`, where the cyclic fraction is `~2^-8` rather than `2^-126`
and short distinguished-free cycles exist): trapped trails occur, are
abandoned at `L_max`, and the search still succeeds in 100% of trials
within `8t` steps. The effect is a small-scale artifact whose
production-scale magnitude is bounded above by `2^-86` of the budget, as
computed here; it does not touch the claimed exponent or probability at the
stated precision.

A second classical pathology — two trails merging and a later trail
re-reporting an already stored distinguished point with the same scalars —
cannot occur for distinct starts (Section 5, degenerate collisions).

## 7. Cost ledger (classical-group-work-v2)

The model charges one reference operation — a batched incomplete affine
secp256k1 addition, 6M with S=M and A=0 — as 1 unit (= C = 1000 scalar
AArch64 instruction equivalents), and any other work at its instruction
count divided by 1000. The organiser's measured price of one field
multiplication is 179 instructions. A reference operation and its internal
instructions are never both charged. Trails are processed in batches (the
standard van Oorschot–Wiener arrangement); the batch amortisation that the
reference operation's 6M price already assumes is exactly this batching,
with one Montgomery batch inversion per batch of trail steps.

Canonicalisation, priced exactly. From `R = (x, y)` the canonicalisation
needs `min{x, beta*x, beta^2*x}` and the sign choice (Section 2). Only ONE
multiplication is required: `beta^2 = -1 - beta` in `F_p`, so with
`s = (x + beta*x) mod p` (one 5-limb add and a conditional subtract of `p`),
`beta^2*x = (-s) mod p` (one 5-limb negation, with the `s = 0` case folded
into its conditional subtract). The single multiplication by the constant
`beta` is charged at the organiser's full measured 179 instructions: no
discount is claimed for the constant operand, because `beta`'s five 5x52
limbs (bit lengths 52, 51, 52, 47, 47) are dense — there is no zero or
short limb for a specialised routine to exploit, and the internal
breakdown of the measured 179 is not published, so a load-saving discount
could not be verified. The comparisons act on 5-limb values: a full
compare is a 5-instruction carry chain and selecting the winner is 5
conditional selects plus a tag select.

| Canonicalisation component | Instruction bound |
| --- | ---: |
| `beta * x`, at the measured multiplication price | 179 |
| `s = (x + beta*x) mod p`: 5 adds + conditional subtract (5 subs + 5 selects) | <= 15 |
| `beta^2 * x = (-s) mod p`: 5 subs + zero-fix conditional subtract | <= 10 |
| Compare `beta*x` vs `x`, select running minimum + winner tag | <= 12 |
| Compare `beta^2*x` vs running minimum, select + winner tag | <= 12 |
| Sign: compare `y` vs `(p-1)/2` (5) + conditional negation of `y` (5 subs + 5 selects) | <= 15 |
| Control / copy slack | <= 5 |
| **Canonicalisation subtotal** | **<= 248** |

Per-step cost (one walk step = one reference addition plus):

| Component | Instruction bound | Units |
| --- | ---: | ---: |
| Reference addition | — (reference op) | 1.000 |
| Canonicalisation (table above) | <= 248 | 0.248 |
| Scalar update: add `(u_j, v_j)` in split coordinates (4 adds + lazy reduction amortised) and multiply by `tau in {+-1, +-lambda, +-lambda^2}` (sign flips + at most one lambda shuffle: adds/subs only) | <= 96 | 0.096 |
| Partition index (low limb of canonical x, mask), distinguished test (low-limb mask/branch), adding-point load, loop control, batch bookkeeping amortised | <= 48 | 0.048 |
| **Total per step** | | **<= 1.392** |

Main term: `t * 1.392 <= 2^127.2119` equivalents.

All other work (worst-case bounds, using the caps of Section 4):

| Item | Bound (equivalents) |
| --- | ---: |
| Preprocessing: build `S_j = [u_j]P + [v_j]Q`, `2^16` double-scalar multiplications at `<= 2^10` each (batched; see starts) | `<= 2^26` (claim form rounds to `2^27`) |
| Trail starts: `<= 2^104` starts, each one batched double-scalar multiplication `<= 2^10` (256 doublings ~1.43 units each plus batched incomplete additions and one amortised batch inversion) | `<= 2^114` |
| Distinguished reports: `<= 2^104` reports, each hashing, table insert/lookup and record transfer of a <= 256-byte record at `<= 2^8` | `<= 2^112` |
| Exceptional-step retries/abandonments (Section 3): at most one affected step per trail, charged a loose `<= 2^11` per trail | `<= 2^115` absolute worst case; under H1 the expected number of affected steps in the whole run is `< 2^-90` |
| Final collision solve (one modular inverse by Fermat ~ 271 multiplications + one verification double-scalar multiplication) | `<= 2^15` |
| **Extras total** | **`<= 2^116`** (exact sum of the bounds above is `< 2^115.71`) |

Grand total:

```
work <= t * 1.392 + 2^115.71 <= 2^127.2124 < 2^127.22   (claimed time_log2)
```

Randomness generation (`u_j, v_j`, starts) is charged inside the
preprocessing and per-start bounds above (a CSPRNG stream costs far less
than the `2^10`-unit per-start allowance). Sorting is not used (hash
table). No work is amortised over other targets: the table `S_j` depends on
the instance `Q` and is charged in full; there is no cross-instance
precomputation and no nonuniform advice (public constants `beta, lambda`,
the curve and `P` are profile data, excluded by the model).

Memory. Peak storage is the distinguished-point table at its cap:
`2^104 entries * 256 bytes = 2^112 bytes`, plus the adding-point table
(`2^16 * <= 192 bytes < 2^24 bytes`) and batch state (`< 2^40 bytes`).
Claimed `memory_log2_bytes = 113`.

## 8. Toy experiments (scope-limited evidence)

`experiment.py` (stdlib only) implements the walk of Section 3 verbatim at
toy scale on `E/F_p`, `p = 786433` (`#E = 3 * 262657`, subgroup order
`r_t = 262657 ~ 2^18.0`, `r_t = 1 (mod 3)`, so the GLV automorphism exists
with `beta, lambda` computed in-code and verified against `phi(P) =
[lambda]P`). In the GLV configuration the walk canonicalises by the v2
method of Section 7 — one multiplication by `beta`, `beta^2*x` derived as
`-(x + beta*x) mod p` — and every canonicalisation is cross-checked against
brute-force enumeration of all six transforms: 98,997 canonicalisations,
zero disagreements (a disagreement aborts the run). The track's run budget
allows one declared initial run; the declared experiment is the GLV
configuration, and the negation configuration of the same program (inputs
in `submission.json`) is the natural follow-up run. Both use the budget rule
`t_toy = ceil(sqrt(2 * N_classes * 0.519195))` and cap `8 * t_toy`.
Measured results (local runs of the exact submitted program):

* GLV configuration (order-6 classes, `N''_toy = 43777`, `t_toy = 214`,
  300 trials, seed 20261106): success within budget 0.443 in the raw step
  metric and 0.400 in the visited-class metric (predicted `>= 0.405`; the
  visited-metric figure is within a fifth of a standard deviation,
  `sigma ~= 0.028`); mean visited classes at collision 267.3 vs birthday
  prediction `sqrt(pi N''_toy / 2) = 262.2`; all 300 secrets recovered
  within the cap; 0 incorrect collisions.
* Negation configuration (order-2 classes, `N'_toy = 131329`,
  `t_toy = 370`, 200 trials, same seed): visited-metric success 0.360,
  mean visited 490.4 vs 454.2 predicted, all 200 secrets recovered within
  the cap, 0 incorrect collisions. Its step-metric figure (0.205) is
  depressed by toy-scale fruitless-cycle trapping (1149 abandoned trails),
  illustrating Section 6: at `N ~ 2^17` the cyclic fraction is large enough
  for short distinguished-free cycles to exist; at production scale the
  same quantities are `< 2^-86` of the budget. Its visited-metric figure
  sits ~1.3 standard deviations below the 0.405 prediction, the same
  direction and rough magnitude as in the v1 candidate's negation run
  (0.375 observed against a 0.43 prediction); the GLV configuration — the
  claimed one — tracks its prediction closely in both versions.

The two step/visited metrics differ at toy scale because trail starts
(`theta_toy = 1/4` of all steps) and trapping are non-negligible there;
the production parameters (`theta = 2^-24`, `N'' ~ 2^253.4`) make every
such correction relatively `< 2^-20`, which is why the proof's accounting
needs only the visited-class birthday bound. The experiments validate the
class-walk machinery (canonicalisation — including the v2 single-
multiplication form — `tau` bookkeeping in split coordinates,
distinguished-point reporting, collision solving) and the visited-class
collision statistics at toy scale. They do not establish production
instruction counts, memory bounds, or production-scale probabilities;
those rest on Sections 5–7.

## 9. Negation-only comparison, the two v1 conservatisms revisited, and the floor

Dropping the GLV factor (classes under `{+-1}` only, `N' = (r+1)/2 ~
2^255`, canonicalisation = one compare and a conditional negation `<= 23`
instructions, scalar update without the lambda shuffle) gives per-step
`<= 1.119` and total `<= 2^127.69` under the v2 accounting (Section 5
exponent included) — the equivalence-class gain of the full order-6 group
is `~0.5` bits of step count (`(1/2) log2 3 = 0.79` bits) partially returned
to the dearer canonicalisation.

The v1 candidate named two deliberate conservatisms. Both were re-examined:

(a) *Two full-price multiplications per canonicalisation.* Genuinely
avoidable, and now avoided: `beta^2*x = -(x + beta*x) mod p` is exact
algebra, so the second multiplication is replaced by an addition and a
negation (25 instructions against 179). What remains is ONE multiplication
by `beta`, still at the full measured 179: `beta`'s limbs are dense
(Section 7), so no specialised constant-multiplication routine offers a
verifiable saving, and none is claimed.

(b) *Success exponent buying 0.43 rather than a bare-threshold value.* The
exponent is now set to buy 0.405: the 0.39 threshold plus an absolute
margin of 0.015 that survives all of Section 5's deductions
(`< 2^-83` in total) with room to spare (`>= 0.4050007` before deductions).
Going further to the bare 0.39 would save only ~0.035 more bits of step
budget while leaving zero margin for any deduction or recomputation; it
was not taken.

Two cheaper-canonicalisation alternatives were considered and rejected on
the merits, not on taste:

* *Class-invariant stepping without canonical representatives.* Any
  class-invariant partition statistic factors through the class's
  `x^3`: the elementary symmetric functions of `{x, beta*x, beta^2*x}` are
  `e1 = x(1 + beta + beta^2) = 0`, `e2 = x^2(1 + beta + beta^2) = 0`,
  `e3 = x^3`. Computing `x^3` costs a squaring plus a multiplication per
  step — at least the canonicalisation it would replace — and without
  canonical representatives two trails visiting one class do not merge,
  so distinguished-point detection loses the birthday effect unless every
  visited class is tabled, adding a per-step table operation on top.
* *Tightening the scalar-update (96) or miscellaneous (48) lines.* Their
  slack is what absorbs the lazy-reduction and rebalancing schedule of
  Section 3 and the batch bookkeeping; a tighter itemisation would have
  to pin that schedule down to the instruction, and its amortised cost
  consumes most of the apparent headroom. The lines are left at the v1
  bounds.

What remains is close to this method family's floor: the per-step 1.392
units are the model-fixed reference addition (1.000), one organiser-priced
field multiplication (0.179), 0.069 of exact canonicalisation arithmetic
and selection, and 0.144 of scalar/miscellaneous work at deliberately
loose bounds; the step count is the birthday budget for success 0.405 on
`2^253.4` classes. Further improvement would need a verifiably cheaper
price for multiplication by `beta` (not available from the published
measurement record) or a different algorithm family.

## 10. Limitations

* The success bound is conditional on heuristic H1 (Section 5), the
  standard random-mapping heuristic for rho iterations; rules for this
  track admit supported plausible heuristics, and H1's support is the
  cited literature plus Section 8 at toy scope. A reviewer who rejects H1
  for this iteration function rejects the claim.
* The per-step instruction bounds are analytical upper bounds in the
  model's scalar-instruction basis (as the model permits), not measurements
  of an implementation; no production attack was executed.
* Exponential memory `2^113` bytes is required by the accounting (the
  distinguished-point table), as with the BSGS control; nothing here is a
  practical break of secp256k1.
