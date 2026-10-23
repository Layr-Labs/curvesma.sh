# Adding-walk discrete logarithm attack on BLS12-381 G1

This document specifies a complete classical attack on
`bls12-381-g1-ecdlp-plain` and derives an integer bound `W` on its total
scalar-AArch64-instruction work. The claimed score is the rational
137562769/1000000, a ceiling of `log2(W / 2000)` in target-group-operation
equivalents (`C = 2000` for this track under `classical-group-work-v2`).
`experiment.py` recomputes every line of the ledger from the same constants,
executes the exact walk function on toy groups, and returns nonzero on any
disagreement.

## 1. Target and conventions

The target is `T = [x]P0` in BLS12-381 `G1`, `x` uniform in `[0, r)`, with

    p = 0x1A0111EA397FE6...AAAB   (field modulus, 381 bits)
    r = 0x73EDA753299D7D...0001   (subgroup order, 255 bits)

Both are the constants of `target-profile.json`. The identity target is
handled before the walk: if `T` is the point at infinity the attack returns
`x = 0` and stops; its cost is inside the setup reserve.

Arithmetic is charged in scalar AArch64 instruction-equivalents. The only
field-multiplication price used is `FMUL = 383`, the pinned `blst_fp_mul`
measurement of `scripts/reference_costs/measurement.json` (`bls_field`
difference column: `(3,137,556 - 1,568,788) / (8,192 - 4,096) = 383`,
Valgrind `Ir` instruction counts on AArch64). The same file measures the
reference `6M` BLS G1 addition at 2,298 instructions while the scoring model
rounds the unit to `C = 2,000`; using the rounded `C` makes every converted
instruction count strictly conservative. Field additions and subtractions are
charged at `FADD = 128` instructions each, an analytical ceiling for a 6-limb
modular add or subtract (six `adds`/`subs` limbs, a conditional subtract of
`p`, and loop overhead), which is consistent with the 383-instruction
multiplication's inclusion of six limb multiplications and a full reduction.

## 2. The walk function

The attack uses a Teske `r`-adding walk with `R = 32` jumps and Floyd's
cycle-finding, so memory is `O(1)` group elements rather than a distinguished-
point table. All walk state is canonical: affine coordinates.

### 2.1 Seeds and jumps

One 64-bit `seed` is drawn. A deterministic 512-bit expansion
`block512(seed, tag)` is defined by eight splitmix64 steps keyed by
`seed ^ ((tag+1) * 0x9E3779B97F4A7C15)`. For `j` in `0..31`, jump `j` is the
point `J_j = [c_j]P0 + [d_j]T` with

    c_j = block512(seed, 2j)   mod (r - 1) + 1   (nonzero mod r)
    d_j = block512(seed, 2j+1) mod r

so every jump is a known `[c_j]P0 + [d_j]T` pair of scalars. The bias of a
512-bit uniform block reduced mod `r` is at most `r / 2^512 < 2^-257` per
coefficient, far below the heuristic allowance in section 5.

### 2.2 Canonical word and partition

For a finite affine point `(x, y)` in `G1`, the canonical 64-bit word is

    word(P) = (x xor (y << 1)) mod 2^64

and `word(infinity) = 0`. Because points are stored in affine form, this word
is a well-defined function of the group element: a Jacobian representative
`(k^2 X : k^3 Y : k Z)` of the same point hashes identically only after the
affine reduction, which is why every step of the walk normalizes. The
partition index is

    h(P) = ((mix64(word(P)) * ((2*seed + 1) mod 2^64)) >> 32) & 31

with `mix64` the four-step bijective mix `w ^= w>>30; w *= 0xBF58476D1CE4E5B9;
w ^= w>>27; w *= 0x94D049BB133111EB; w ^= w>>31` (all mod `2^64`). `h` is
deterministic, depends only on the group element and the seed, and selects one
of the 32 jumps.

### 2.3 Iteration and recovery

The walk state is `(P, alpha, beta)` with the invariant
`P = [alpha]P0 + [beta]T`. One adding step replaces `P` by `P + J_{h(P)}` in
affine coordinates and updates `(alpha, beta) += (c_j, d_j) mod r`. Tortoise
and hare both start at `(infinity, 0, 0)`. Each Floyd iteration advances the
tortoise once and the hare twice, so three adding steps per iteration.

At each iteration the states' points are compared. On the first equality with
`(beta_h - beta_t)` nonzero mod `r`,

    [alpha_t]P0 + [beta_t]T = [alpha_h]P0 + [beta_h]T
    => (beta_h - beta_t) * x = alpha_t - alpha_h   (mod r, since T = [x]P0)
    => x = (alpha_t - alpha_h) * (beta_h - beta_t)^(-1)  mod r

The candidate is verified by `[x]P0 = T` before return. If `beta_h = beta_t`,
the iteration is skipped and the walk continues; such a beta collision on one
seed is counted as an ordinary failure inside the success bound of section 5.
The budget is `I = 3 * ceil(sqrt(r))` Floyd iterations, i.e.
`S = 9 * ceil(sqrt(r))` adding steps; exhausting the budget fails the attack.

## 3. Instruction ledger

`ceil(sqrt(r)) = 228,988,810,152,649,578,064,853,576,960,394,133,504`, so

    I = 686,966,430,457,948,734,194,560,730,881,182,400,512   (iterations)
    S = 2,060,899,291,373,846,202,583,682,192,643,547,201,536  (adding steps)

### 3.1 One affine adding step

Each step inverts one field element by Fermat `z^(p-2) mod p` using the plain
binary chain. `p - 2` has `bit_length 381` and Hamming weight `229`, so the
chain is `380` squarings plus `228` additions: `608` Fp multiplications, an
explicit construction, not an estimate. Charged:

    inversion        608 * 383 + 381 * 32 scalar loop bits   =   245,056
    lambda, x3, y3   4 multiplications + 8 adds/subtracts    = 4*383 + 8*128
    per adding step                                          =   247,612 instructions

The affine add `P + J` computes `inv = (x2 - x1)^(p-2)` (the inversion above),
then `lambda = (y2 - y1) * inv`, `x3 = lambda^2 - x1 - x2`,
`y3 = lambda * (x1 - x3) - y1`: three multiplications suffice and four are
charged; the eight additions/subtractions cover `y2 - y1`, `x2 - x1`,
`lambda^2 - x1 - x2`, `lambda * (x1 - x3) - y1` and operand moves. The double
case `P + P` uses the tangent `3 x1^2 / (2 y1)` inside the same 4M + 8A bound.

### 3.2 Per-iteration and per-step auxiliaries

    partition word + mix64 + product + mask   256 instructions per adding step
    (alpha, beta) updates (3 steps worth)     1,536 instructions per iteration
    affine equality test (generous bound)     8 FMUL + 8 FADD = 4,088 per iteration

The equality test's affine `x` compare needs no multiplication; 8M + 8A is a
deliberately inflated bound covering compares, moves and control. The
coefficient updates are two 255-bit modular additions per step; 1,536 per
iteration (512 per step) covers them with margin.

### 3.3 Setup, resolve and reserve

    two 512-bit derivations      2 * 64 * (8*40 + 16,384)   =      2,138,112
    64 scalar multiplications    64 * 510 * (16M + 28A)
                                 = 64 * 510 * 9,712        =    316,999,680
    affine normalize 32 jumps    32 * 245,056               =      7,841,792
    output, resolve, verify overhead                        =      5,500,000
    stream and loop setup                                   =      1,000,000
    setup total                                           ≈  333,479,584

The 64 scalar multiplications build `[c_j]P0` and `[d_j]T` for the 32 jumps at
the charged Jacobian cost `16M + 28A` per group operation and `≤ 510` double-
and-add operations per 255-bit scalar; the same line covers the identity
check, the initial `T` handling, and the final `[x]P0 = T` verification. A
flat `2^40 = 1,099,511,627,776` instruction reserve is added for bookkeeping,
encoding and any residual control flow; it is larger than the setup block by
a factor of three and is not needed anywhere else.

### 3.4 Total

    steps       S * 247,612   = 510,303,395,335,660,805,914,150,715,084,854,009,666,732,032
    partitions  S * 256       =     527,590,218,591,704,627,861,422,641,316,748,083,593,216
    coefficients I * 1,536    =   1,055,180,437,183,409,255,722,845,282,633,496,167,186,432
    equality    I * 4,088     =   2,808,318,767,712,094,425,387,364,267,842,273,653,293,056
    setup                       =                 333,479,584
    reserve                     =         1,099,511,627,776
    ─────────────────────────────────────────────────────────────────────────────
    W = 514,694,484,759,148,014,223,122,347,276,647,627,415,912,096 instructions

All six lines are integers; `experiment.py` recomputes each from the same
constants and asserts the total. In units of `C = 2,000`,
`W / 2000 ≈ 2.5735 × 10^41`, and `log2(W / 2000) = 137.56276…`. The claimed
rational `137562769 / 1000000 = 137.562769` is a ceiling; the program asserts
both `137.562769 ≥ log2(W/2000)` and `137.562768 < log2(W/2000)`.

Preprocessing that does not depend on the challenge is only the public
generator and parameter load, `ceil(2^30 / 2000) = 536,871` units
(`log2 = 19.034216`), and is already inside `W`'s setup line. Peak memory is
bounded by `2^20` bytes: 32 affine jump points (2 Fp elements each, ~3 KiB),
two walk states, the 64 jump coefficients, and code; the bound is reported as
a cap, not as a minimum. Nonuniform advice is zero bytes: jumps, coefficients
and the walk are derived from public constants and a single fresh seed.

## 4. Correctness

On equality with `beta_h ≠ beta_t` the algebra of section 2.3 produces the only
residue consistent with `[alpha_t - alpha_h] = (beta_h - beta_t) x`, and the
explicit `[x]P0 = T` check rejects it if the walk's group arithmetic were ever
wrong. The experiment additionally asserts the returned scalar equals the
drawn secret on every successful draw, exercising the sign of the recovery
formula — `(alpha_t - alpha_h)`, not its negation — against the group law.
The identity target is returned before the walk.

## 5. Success bound and heuristics

Claimed success probability is `0.39`, the model's floor. The support is:

* **Random-mapping heuristic.** A Teske `r`-adding walk with `R = 32` jumps
  and an independent-looking partition is treated as an approximately random
  mapping on the `r` group elements. This is the standard heuristic behind
  every published rho estimate for pairing-friendly curves; Teske's
  experiments show `r`-adding walks with `R ≥ 20` already behave close to
  random mappings (average collision-length multiplier ≈ 1.42–1.45 versus the
  ideal `sqrt(pi/2) ≈ 1.253`), and Bernstein–Lange–Schwabe's negation-map
  measurements confirm the family on elliptic groups. It is not proved that
  this particular `mix64`/word composition is independent; that is the
  declared heuristic.
* **Budget.** Under the random-mapping model the first repeat appears at
  `≈ sqrt(pi r / 2)` steps in expectation, and Floyd detects it within roughly
  three times the collision depth; the `3·ceil(sqrt(r))` iteration budget
  (9 adding steps each) leaves a wide margin over the observed toy
  distributions below.
* **Sample.** `experiment.py` runs the exact production walk function —
  same `block512` jump derivation, same `mix64` partition, same Floyd budget —
  on prime additive groups of order `10,007`, `100,003` and `1,000,003`
  (80 draws each, seeds and secrets derived from `master = 12,648,430`), and
  on the cofactored curve `y^2 = x^3 + 4` over `F_9613` (subgroup of prime
  order `3,187`, real affine curve arithmetic including the inversions that
  dominate `W`; 40 draws). Successes in the shipped initial run are
  `80/80`, `80/80`, `80/80`, `40/40`; the worst observed iteration length is
  `2.9·sqrt(n)`, inside the `3·sqrt(n)` budget. The sample is evidence at its
  actual scope only: it supports the shape of the iteration-length
  distribution and the correctness of the recovery, not the production
  constant itself.
* **Why 0.39 and not higher.** The measured sample is far above the floor;
  a sharper constant would require an independence proof this document does
  not give. Beta collisions, seed degeneracy and the heuristic gap are all
  charged to the `0.61` failure allowance.

## 6. Limitations

* Analytic attack with a stated random-walk heuristic; no production-size run
  exists or is possible. The toy experiments bound the distribution at their
  own scope, not at `2^127.5`.
* `FMUL = 383` is the pinned `blst_fp_mul` measurement; a faster BLS field
  multiplication would only lower `W`. `FADD = 128` is an analytical ceiling
  for a 6-limb modular add; the margin to the next ledger line is large.
* This is a worst-case-work bound at success `0.39`, not a statement that
  BLS12-381 is broken in practice.
