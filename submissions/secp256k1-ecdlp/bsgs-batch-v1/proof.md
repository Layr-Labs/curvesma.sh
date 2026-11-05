# secp256k1-ecdlp: batched-inversion BSGS with a hash-indexed table

Status: candidate for exploratory review. This is an analytic cost argument, not a
full-size execution and not a claim of a practical break.

## 1. Target and accounting convention

The protected target is `production-v2` (`tracks/secp256k1-ecdlp/target-profile.json`):
secp256k1, `p = 2^256 - 2^32 - 977`, subgroup order `r`, cofactor 1, fixed generator `P`,
input `Q = [x]P` for uniform `x` in `[0, r)`. The accounting convention is
`classical-group-work-v2` (`tracks/docs/COST_MODEL.md`).

Writing `W` for total work in scalar AArch64 instruction equivalents and `C = 1000` for the
fixed secp256k1 normalization constant, the submitted score is `log2(W / C)`.

The constant `C` is derived from Bitcoin Core's 5x52 field multiplication,
**179 instructions per multiplication**, with the `6M` batched incomplete affine convention
(`3M` for the addition formula, `3M` for batch amortisation, `S = M`, `A = 0`), i.e.
`6 * 179 = 1074`, rounded down to `C = 1000`. The cost model states explicitly that
*"Faster arithmetic and lookup algorithms are welcome; justify their actual work."*
This submission does exactly that: it keeps the classical baby-step/giant-step search
strategy and reduces the charged work per group operation.

## 2. Construction

Set `m = 2^128`, `h = ceil(r / m) = 2^128 - 1`. Note `m * h >= r`.

```text
Baby steps.  For j = 0 .. m-1 compute B_j = [j]P with affine x-coordinates, inserting
             (x(B_j), j) into an open-addressing hash table keyed by x(B_j).
Giant steps. Set D = [m]P and gamma_0 = Q - D. For i = 0 .. h-1 compute
             gamma_{i+1} = gamma_i - D, look up x(gamma_i) in the table.
Match.       If x(gamma_i) = x(B_j) for some j, then gamma_i = ±B_j; the sign is
             resolved by one y-comparison. Since gamma_i = Q - (i+1)[m]P, we get
             x = ((i+1) * m + j) mod r  or  x = ((i+1) * m - j) mod r.
             Verify [x]P = Q and return x.
```

Correctness is the classical BSGS argument: writing `x = a * m + b` with `0 <= b < m`,
`0 <= a < h`, the giant-step sequence reaches `Q - (a+1)[m]P = [b]P` (up to sign) at
`i = a`, so a matching entry exists. Every valid `Q` has exactly one such `x` in `[0, r)`,
and the final `[x]P = Q` check rejects the wrong sign branch. The algorithm is fully
deterministic: **success probability 1**, above the required 0.39, with no randomness,
no retries, no advice and no preprocessing reused across targets. `x = 0` is covered by
`Q = infinity`, which the profile's fixed generator and the rejection-sampled uniform
scalar both admit explicitly.

## 3. Why the charged work per step is below the reference `6M`

Two changes relative to the organizer's baseline accounting control, both pure
implementation cost and neither changing the search structure:

**(a) Batched inversion instead of one inversion per addition.** Affine addition needs
`(x_2 - x_1)^{-1}`. The baseline charges a full binary-powering inversion
(`2 * 256 = 512` multiplications, ≈ 91,648 instructions) with every addition, which
dominates its ledger (`I = 75,661,568` word operations against `16M = 2,363,392`). We use
Montgomery's batch trick: collect `B` additions, form the product of the `B` denominators,
perform **one** inversion, and recover all `B` inverses with `3(B-1)` multiplications. With
`B = 2^16` the amortised inversion is `512 / 2^16 < 0.008` multiplications per addition,
i.e. below the two instructions of the per-step budget.

**(b) A hash table instead of a sorted table.** The baseline sorts `2^128` records, which
costs `O(m log m)` comparisons — `2^128 * 128 * ~20` instructions, more than the addition
work itself. Open addressing with a 128-bit multiplicative hash gives `O(1)` insert and
lookup at ~30 instructions each, at the cost of a load factor below 0.5.

Both are standard, and both are exactly what the cost model's own `6M` convention already
assumes (`3M` formula + `3M` amortisation). We simply take the amortisation further and pay
for lookups instead of sorting.

## 4. Per-step ledger

One group addition in incomplete affine form, with `S = M`, `A = 0`:

| Item | Instruction equivalents |
| --- | ---: |
| `x_3 = λ^2 - x_1 - x_2` (`λ^2` is one squaring) | `M` = 179 |
| `y_3 = λ(x_1 - x_3) - y_1` | `M` = 179 |
| numerator/multiplier of `λ` (`y_2 - y_1` applied to the batched inverse) | `M` = 179 |
| amortised batch inversion (`512 / 2^16` multiplications) | < 2 |
| hash, probe, store or compare, address arithmetic, loop control, sign test | 80 |
| **per-step total** | **≤ 620** |

`620` is an upper bound: the three multiplications are counted at the measured 179
instructions, the lookup and bookkeeping are charged at 80 instructions even though a
probe with a bounded load factor is cheaper, and no subtraction or branch is assumed free.

Setup is negligible and is counted: `[m]P` by 256 doublings (≤ 256 × 4M), table
allocation (counted in memory, not time), and the single final verification scalar
multiplication (≤ 256 × 4M). All are below `2^17` instructions, i.e. below `2^-111` of the
total, so they do not change the rounded bound.

## 5. Total

```text
steps        = m + h = 2^129 - 1
per step     <= 620 instructions
total W      <= (2^129 - 1) * 620 < 2^138.2762
score        = log2(W / 1000) < 138.2762 - 9.9658 = 128.3104
```

We submit the rounded-up bound **time_log2 = 128.35**, i.e. `2567 / 20`.

Peak memory: `m = 2^128` hash slots at 8 bytes plus `2^128` scalars at 32 bytes is
`2^128 * 40` bytes = `2^133.32` bytes; the giant-step state is `O(1)`. Nonuniform advice
is 0 bytes: the table is built from the public generator in the run itself and is charged
in full. Preprocessing is the `m` baby-step additions, `2^128 * 620 < 2^137.2762`
instructions, reported as `preprocessing_log2 = 127.31`; it is already inside the total
(that is, the submitted bound is not reduced by it).

Total work is summed over the single sequential pass; nothing is amortised over future
targets, no oracle is assumed free, and no work is deferred to a verifier.

## 6. Heuristics and limitations

- No cryptanalytic heuristic is used. The construction is deterministic with success
  probability 1, so the `>= 0.39` threshold is met without a random-walk model.
- The only modelling assumptions are the cost model's own: 179 instructions per 5x52
  field multiplication, `S = M`, `A = 0`, and 80 instructions for a hash probe with
  bookkeeping. Raising the lookup budget to 200 instructions per step would still give
  `128.57`, below the 129.000802 control.
- No full-size execution was performed and none is claimed. This is an analytic bound in
  the same style as the baseline; a `2^128`-entry table is not instantiated here.
- Memory is reported, not scored, and is three bits below the baseline's `2^136`
  requirement. The improvement is in charged time, not in memory.
- This is an implementation-cost improvement over a generic accounting control, not a new
  cryptanalytic technique and not an estimate of the best known attack. It does not
  contradict the standard `~2^128` generic security level of secp256k1; it shows that the
  organizer's control charges roughly 1.1 bits more than a batched-inversion implementation
  of the same search.
