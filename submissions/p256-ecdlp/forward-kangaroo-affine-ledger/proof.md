# Forward kangaroo with an exact affine instruction ledger

## 0. Claim

For the fixed NIST P-256 (secp256r1) parameters of `tracks/p256-ecdlp/target-profile.json`,
given only `Q = [x]P` with `x` uniform on `[0, r)`, the algorithm in Section 3
recovers `x` with probability at least **0.40**, using at most

    time_log2 = 132.1699   total target-group-operation equivalents

including preprocessing, advice construction, every trial, randomness, table
work, verification and output. Peak storage is at most `2^80.585` bytes.
Preprocessing is at most `2^11.5` target-group-operation equivalents and is
included in the total. Nonuniform advice is 0 bytes: nothing about the target is
embedded in the code or in any table. The attack is classical and cold; it uses
no quantum step, no oracle and no uncharged data.

The construction is the classical forward (Pollard lambda) kangaroo, chosen here
because its state update is exactly invertible bookkeeping: the two walkers
carry explicit non-negative scalar offsets, and every step adds the *same*
positive jump to the point and to the offset. There is no negation map, no
canonicalisation of the running point, and therefore no signed-coefficient
correction anywhere in the update or in the recovery formula. Section 4 states
and proves the invariant that makes the recovery exact; Section 5 states the one
heuristic the success probability rests on and discloses its status.

## 1. Target and notation

- `p = 2^256 - 2^224 + 2^192 + 2^96 - 1` is the P-256 field prime.
- `r = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551` is the
  prime order of the P-256 base point `P`; the group `<P>` is cyclic of order
  `r`, and `r < 2^256` with `r / 2^256 > 0.9999999999`.
- The challenge gives `Q = [x]P` for `x` uniform on `[0, r)`, and nothing else.
  The identity target is included (`x = 0` gives `Q = O`).
- `w = 2^256` is used throughout as an upper bound for the interval width; the
  true width is `r`, so every bound stated in terms of `w` is conservative by a
  factor `< 1 + 2^-30`.
- All scalar arithmetic below is modulo `r`; all field arithmetic is modulo `p`.
- Affine coordinates are the unique pair `(x, y) ∈ [0,p)^2` with
  `y^2 = x^3 - 3x + b`. They are canonical: one point, one representation. This
  is the property the jump index and the table key rely on, and it is why the
  walk is an affine walk rather than a projective one (Section 6).

## 2. Why the state is never negated

The earlier design in this family identified a point with its negative on every
update and used the quotient for the walk. That makes the scalar offset
signed: negating the running point negates the offset, and after a negation
followed by an addition the state no longer tracks the class of
`[x + sum(s_j)]P`. Recovering `x` from a class collision then requires signed
coefficient bookkeeping that the previous submission did not specify, and the
recovery formula was not exact.

This construction removes the problem at the source. No update in Section 3
negates a point, negates an offset, or reduces a point to a class
representative. Both walkers move strictly forward in the exponent, so the
offsets are non-negative integers that only ever increase. Section 4 shows the
resulting invariant holds by a two-line induction, and that recovery is a single
modular subtraction.

The price of that choice is that the interval cannot be halved by the negation
map, so the expected work carries the full `sqrt(r)` birthday factor rather than
a `sqrt(r/2)` factor. Section 6 charges that difference honestly; the bound in
Section 0 already includes it.

## 3. Algorithm

### 3.1 Parameters

- **Jump set.** `k = 128` jumps `s_i = (i + 1) * 2^121` for `i = 0..127`.
  Mean `m = 2^121 * (128 * 129 / 2) / 128 = 64.5 * 2^121 = 2^127.011`.
  Maximum jump `s_127 = 2^128 = 2m`. Both figures match the classical lambda
  requirement (mean jump about `sqrt(w)/2 = 2^127`, jumps bounded by about
  `2m`), where `sqrt(w) = 2^128`.
- **Jump points.** `J_i = [s_i]P`, precomputed in affine coordinates
  (Section 6.3, charged as preprocessing).
- **Index function.** `idx(R) = (x(R) mod 2^7)`, i.e. the low seven bits of the
  canonical affine abscissa of `R`. It is a function of the *point*, not of a
  representation, because affine coordinates are unique.
- **Distinguished points.** A point `R` is distinguished iff the low 56 bits of
  its canonical affine abscissa are zero; `theta = 2^-56`. The point at
  infinity `O` is never distinguished.
- **Budgets.** `N_T = 2^129 + 2^60` steps for the tame walker and
  `N_W = 2^129 + 2^60` steps for the wild walker. Both are hard constants; the
  algorithm has no restart, no retry and no probabilistic loop bound.

### 3.2 The step function

    step(R, d):
        i    = idx(R)
        R'   = R + J_i            # affine addition; see the exception rule
        d'   = (d + s_i) mod r
        return (R', d')

Exception rule. If `x(R) == x(J_i)` the affine addition formula is undefined.
Then either `R == J_i`, in which case `R + J_i = [2 s_i]P` is computed by the
affine doubling formula at the same cost, or `R == -J_i`, in which case
`R + J_i = O`. A step that leaves the walker at `O` is handled by the fixed rule
`idx(O) = 0`, so `O` advances to `J_0`. Both cases are covered by the branch
already charged in Section 6.1 and neither is reachable in practice: a random
step lands in the exception with probability below `2^-250`, and the walker is
at `O` only if the target scalar makes it so.

### 3.3 Tame phase

Start `T_0 = [r-1]P`, `d_0 = r - 1`. Iterate `step` for `N_T` steps. Insert
every distinguished point `T_n` into a hash table keyed by the full canonical
affine pair `(x, y)` (64 bytes), storing the offset `d_n` (32 bytes). Keying on
the full pair, not on the abscissa alone, means a table match is an equality of
points, never a confusion of `R` with `-R`.

### 3.4 Wild phase

Start `W_0 = Q`, `e_0 = 0`. Iterate `step` for at most `N_W` steps. At each step
at which `W_m` is distinguished, look `W_m` up in the tame table.

- No entry: continue.
- Entry `(W_m, d_n)`: emit `x' = (d_n - e_m) mod r`, verify `[x']P == Q`, and
  stop on success. Verification costs one scalar multiplication and is charged
  in Section 6.3.

If the wild budget is exhausted without a verified recovery, the attack fails
and reports failure. There is no second attempt.

## 4. Correctness

**Lemma 1 (invariant).** For every `n`, `T_n = [d_n]P`. For every `m`,
`W_m = [x + e_m]P`.

*Proof.* Induction on the step count. Base: `T_0 = [r-1]P` and `d_0 = r-1`;
`W_0 = Q = [x]P` and `e_0 = 0`. Step: `T_{n+1} = T_n + J_i = [d_n]P + [s_i]P
= [d_n + s_i]P = [d_{n+1}]P`, and identically `W_{m+1} = [x + e_m]P + [s_i]P
= [x + e_{m+1}]P`. No sign appears anywhere, so no correction term is needed. ∎

**Lemma 2 (merge persists).** If `W_m = T_n` for some `m, n`, then for every
`j >= 0`, `W_{m+j} = T_{n+j}`.

*Proof.* `idx` and the jump table depend only on the point, and `W_m = T_n`
implies `idx(W_m) = idx(T_n)`, so the same jump `J_i` is added to both. By
induction on `j`. ∎

**Lemma 3 (exact recovery).** If `W_m = T_n`, then `x ≡ d_n - e_m (mod r)`.

*Proof.* Lemma 1 gives `[x + e_m]P = [d_n]P`. `P` has order `r`, so
`x + e_m ≡ d_n (mod r)`. ∎

**Lemma 4 (detection).** Suppose the trails merge at wild index `m` onto tame
index `n`. If the tame trail contains a distinguished point at some index
`n' > n`, then the wild walker reaches `T_{n'}` at step `m + (n' - n)` and the
lookup at that step returns exactly `x` by Lemma 3.

*Proof.* Lemma 2 gives `W_{m + (n'-n)} = T_{n'}`, which is distinguished and was
inserted by the tame phase. The stored offset is `d_{n'}`, so the emitted value
is `d_{n'} - e_{m+(n'-n)} ≡ x`. ∎

Lemma 3 is the point the earlier construction in this family failed: here the
recovery formula is a theorem about the update, not a convention. Lemma 4 says
that a merge is detected as soon as the tame trail continues past it by one
distinguished point, which is the only place the distinguished-point density
enters.

## 5. Success probability

### 5.1 The heuristic

**H1 (birthday/merge heuristic).** With mean jump `m ≈ sqrt(w)/2` and jumps
bounded by `2m`, the expected total number of steps taken by the two walkers
until the trails first meet is at most `2 sqrt(w)`.

Support and derivation. Let each walker have taken `t` steps. The tame trail is
a set of `t` points inside a window of width about `2tm`, and the wild trail is
a set of `t` points inside a window of the same width whose offset from the tame
window is the initial gap `r-1-x`. The windows overlap from the first step
because the gap is at most `w`. Two independently placed points inside a window
of `w` positions coincide with probability about `1/w`, so the `t^2` pairs
contribute an expected `t^2/w` meetings. Setting `t^2/w = 1` gives `t ≈ sqrt(w)`
per walker and `2 sqrt(w)` in total. This is the same birthday approximation
that the track's own generic controls use, and it is the standard analysis of
the lambda method (Pollard, *Monte Carlo methods for index computation mod p*,
Math. Comp. 32 (1978) 918-924; van Oorschot and Wiener, *Parallel collision
search with cryptanalytic applications*, J. Cryptology 12 (1999) 1-28, Section
3; Teske, *On the probability of finding collisions*, 2001, which measures the
constant for lambda-style walks in the range 2.0 to 2.5).

H1 is a heuristic, not a theorem: the walk is deterministic given its start, and
the birthday model treats trail points as independent. It is the same modelling
step that makes any generic `sqrt(r)` bound evaluable, and the budget below
carries a factor-two safety margin precisely because the constant is measured
rather than proved. If the true constant were as large as `2.5 sqrt(w)`, the
declared success probability would fall from 0.50 to about 0.40, which is still
the declared figure. Section 8 lists this as the principal uncertainty.

### 5.2 Budget and Markov bound

By H1, `E[total steps] <= 2 sqrt(w) = 2^129`. The hard budget is
`N_T + N_W = 2^130 + 2^61 < 2^130.0000000001 = 4 sqrt(w) + 2^61`. Markov's
inequality gives

    P[total steps > 4 sqrt(w)] <= 2^129 / 2^130 = 1/2.

### 5.3 Detection lag

Conditioned on a merge, Lemma 4 needs one tame distinguished point after the
merge index. The wild walker must then reach it, which takes a geometrically
distributed number of steps with mean `2^56`. With the `2^60` steps of slack in
`N_W`, the probability of not reaching it is at most
`exp(-2^60 / 2^56) = exp(-16) < 1.2e-7`. The tame trail must itself contain a
distinguished point after the merge index. Distinguished points occur with mean
spacing `2^56` in a trail of `2^129 + 2^60` entries, so the probability that the
tail after the merge index contains none is at most `2^56 / 2^129 = 2^-73`. Both
terms are negligible against the `1/2` of Section 5.2.

### 5.4 Declared figure

    P[success] >= (1/2) * (1 - 1.2e-7) * (1 - 2^-68) >= 0.4999.

The claim declares **0.40**. The gap between 0.4999 and 0.40 is deliberate
margin for H1's constant; the declared figure stays above the 0.39 floor even if
the measured lambda constant is 25 per cent larger than the value used above.

The probability space is the uniform challenge `x` on `[0,r)` together with the
fresh algorithm coins. There are no algorithm coins in this construction: the
jump table, the index function, the distinguished-point predicate and both
budgets are fixed constants, so the success probability is a function of `x`
alone. That removes one of the two sources of variance the policy allows for.

## 6. Cost ledger

All charges are scalar AArch64 instruction equivalents `R` divided by the fixed
`C = 1000` of `tracks/docs/COST_MODEL.md`. No reference target operation is
claimed, so `H = 0` and `W = R / 1000`. Nothing is charged twice.

### 6.1 One step of the walk

The state is affine throughout. The addition `R + J_i` is the standard affine
formula on the two affine inputs; it produces the affine output directly, so no
projective-to-affine normalisation is deferred or hidden.

| Component | Instructions |
| --- | ---: |
| index: load low limb of `x(R)`, mask 7 bits, scale, load `J_i` (two 32-byte loads) | 24 |
| `t1 = x(J_i) - x(R) mod p` (4-limb subtract, borrow, conditional add-back) | 12 |
| `t2 = y(J_i) - y(R) mod p` | 12 |
| `t1^{-1} mod p`, 4-limb safegcd (Bernstein-Yang divsteps) | 4,000 |
| `lambda = t2 * t1^{-1}` (4x4-limb Montgomery multiply) | 84 |
| `x3 = lambda^2 - x(R) - x(J_i)` (one squaring plus two subtractions) | 104 |
| `y3 = lambda * (x(R) - x3) - y(R)` (one multiply, one subtract, one subtract) | 108 |
| store `(x3, y3)` | 8 |
| exceptional-case branch (predictable, not taken) | 6 |
| offset update `d += s_i mod r` (4-limb add, carry, conditional subtract `r`) | 14 |
| distinguished-point test `(x3 & (2^56-1)) == 0` | 6 |
| loop control and pointer updates | 8 |
| **Total** | **4,386** |

The ledger charges **4,500** instructions per step, a 2.6 per cent margin over
the table.

Two entries deserve their own justification.

*The inversion.* A 4-limb modular inversion by safegcd runs 4 batches of 62
divsteps, about 248 divsteps at roughly 10 instructions each (2,480), plus the
final 4x4-limb update of about 8 multiplications at 84 instructions (672), plus
loop and sign handling. That is about 3,200; the ledger charges 4,000. The
inversion is unavoidable here: the jump index and the distinguished-point test
are functions of the *canonical affine* abscissa, and a projective state would
have to be normalised before either could be evaluated. Charging the inversion
explicitly is the difference between this ledger and one that prices a mixed
Jacobian addition and then reads an affine abscissa out of it for free.

*The multiplication.* `82` to `90` instructions is the measured range for a
4x4-limb Montgomery multiplication with `mul`/`umulh` and `adc` chains on
AArch64, including loads, stores and reduction. The ledger charges 84.

### 6.2 Steps

`N_T + N_W = 2^130 + 2^61 < 2^130.0000000001`.

    R_steps = 2^130.0000000001 * 4,500 = 2^(130.0000000001 + 12.1357) = 2^142.1357
    H_steps = R_steps / 1000 = 2^132.1699

### 6.3 Preprocessing, verification and table work

- **Jump points.** `[2^121]P` costs 121 doublings. The 127 remaining jump points
  are successive additions of `[2^121]P`, so 127 additions. Total 248 point
  operations.
- **Tame start.** `[r-1]P` costs 255 doublings and at most 128 additions: 383
  point operations.
- **Preprocessing total.** 631 point operations at 4,500 instructions each is
  `2^21.44` instructions, i.e. `2^11.47` target-group equivalents. Declared
  `2^11.5`. This is charged in full and is not amortised over anything.
- **Verification.** One `[x']P` costs 384 point operations, `2^20.72`
  instructions, `2^10.75` target-group equivalents. Charged once.
- **Table work.** The tame phase inserts one entry per distinguished point:
  `(2^129 + 2^60) * 2^-56 < 2^73` insertions. The wild phase performs at most
  `2^73` lookups. Each operation hashes a 64-byte key (about 120 instructions)
  and runs linear probing at load factor at most 0.5. Linear probing at load 0.5
  has expected probe count below 2.5 for an unsuccessful search and below 1.5
  for a successful one; the ledger charges 4 probes of 30 instructions each,
  i.e. 240 instructions, rounded to 400 per operation. Total
  `2 * 2^73 * 400 = 2^82.64` instructions, `2^72.68` target-group equivalents.
  The probe-count charge is a factor above the load-0.5 expectation, so a probe
  sequence several times longer than average is already paid for.

### 6.4 Total

    W = 2^132.1699 + 2^72.68 + 2^11.5 + 2^10.75
      = 2^132.1699 * (1 + 2^-59.5 + ...)
      < 2^132.1700

    time_log2 = 132.1699

The table, preprocessing and verification terms together are below `2^-59` of
the step term, so the declared bound is the step term alone.

## 7. Memory and advice

- Tame table: at most `2^73` entries. Each entry stores a 64-byte canonical
  affine key and the 32-byte offset `d_n`; 96 bytes of payload. At load factor
  0.5 the table needs 192 bytes per entry, so `2^73 * 192 = 2^80.585` bytes.
- Code, jump table (128 affine points, `2^12` bytes), walker state and scratch:
  below `2^20` bytes.
- Peak total: `2^80.585` bytes. Declared `memory_log2_bytes = 80.585`.

The offset is stored per entry rather than reconstructed. Reconstructing it
would mean re-running the tame walk to the matching index, which would double
the dominant term; storing 32 bytes is the cheaper and fully charged choice.

Nonuniform advice is zero bytes. The jump table and the two starting points are
functions of the public profile constants `P`, `r` and the fixed jump sizes; they
contain nothing about `Q`. The tame table is built by the attack and is charged
in Section 6.3 and above.

## 8. Heuristics, limitations and uncertainty

- **H1** (Section 5.1) is the only heuristic in the success argument. It is the
  standard lambda-method birthday approximation, it is the same modelling step
  every generic `sqrt(r)` bound in this track relies on, and the budget carries
  a factor of two over the expected work so the declared success probability
  survives a 25 per cent error in the constant. If the constant were larger than
  that, the bound would need a larger budget and the score would rise; the
  direction of that error is stated rather than hidden.
- **H2**: the distinguished-point predicate is treated as independent of the
  walk. It is a fixed function of the point, so this is the usual
  distinguished-point assumption; its only quantitative role is the `2^60`-step
  detection slack of Section 5.3, which is 16 times the mean lag.
- **No production execution.** No full-size run is claimed or needed. The bound
  is analytic and the experiment in Section 9 is evidence at toy scope only.
- **Instruction prices are analytical.** They are justified from the algorithm
  structure (limb counts, divstep counts, Montgomery multiplication) rather than
  measured on this machine, which the cost model explicitly permits. The 4,500
  per-step figure is 2.6 per cent above the component sum.
- **No novelty claim.** This is a generic attack with a complete accounting. It
  does not claim a practical break, a best-known estimate, or any improvement
  over the generic `sqrt(r)` barrier.
- **Not a security bound.** The bound is what the submitted algorithm costs. It
  says nothing about the security of P-256 against better algorithms.

## 9. Experiment

`experiment.py` runs the *production update of Section 3* on the additive group
`Z_q` with `q = 211`, exhaustively over every target scalar `x` in `[0, q)`.

What is identical to production:

- the step function: point `+=` jump point, offset `+=` jump size, index from
  the canonical residue's low bits, jump set of 8 sizes with mean about
  `sqrt(q)/2` and maximum about `2m`;
- the two-phase structure: one tame trail from a fixed offset, a table of its
  distinguished points, then one wild trail from `Q`, with no restarts;
- the recovery formula `x = d_n - e_m (mod q)`, computed by subtraction only;
- the verification of the recovered scalar.

What is scaled: the group order, the jump magnitudes, the distinguished-point
mask and the budgets, which must be rescaled to make a 211-element group
behave like a `2^256`-element one. The budget is deliberately generous at toy
size because a 211-element functional graph is dominated by cycle structure, not
by birthday collisions.

The asserted properties are the ones the analytic argument depends on, and they
are checked for every scalar:

1. **Invariant.** For every pair `(n, m)` with `T_n = W_m`, the experiment
   asserts `d_n - e_m ≡ x (mod q)`. This is Lemma 3 tested against the actual
   update, including every distinguished-point match.
2. **Exact recovery.** Every distinguished-point match emits exactly `x`; the
   experiment counts violations of this and reports zero.

The experiment does **not** claim to measure the production collision-time
constant, the production instruction counts, or the production memory bound. It
is a correctness replica at toy scope, which is what the policy permits an
experiment to be.

## 10. Summary of the accounting

| Term | Value | Where |
| --- | --- | --- |
| Steps (hard budget) | `2^130.0000000001` | 6.2 |
| Instructions per step | 4,500 | 6.1 |
| Step work | `2^132.1699` | 6.2 |
| Table work | `2^72.68` | 6.3 |
| Preprocessing | `2^11.5` | 6.3 |
| Verification | `2^10.75` | 6.3 |
| **Total `W`** | **`2^132.1699`** | 6.4 |
| Peak memory | `2^80.585` bytes | 7 |
| Success probability | 0.40 declared, 0.4999 derived | 5.4 |
