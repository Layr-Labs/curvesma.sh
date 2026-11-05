# secp256k1-ecdlp: interleaved-chain BSGS with a deterministic lookup bound

Status: candidate for exploratory review. Analytic bound, no full-size execution claimed.

## 1. Target and accounting convention

`production-v2` secp256k1: `p = 2^256 - 2^32 - 977`, prime subgroup order
`r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141`, cofactor 1,
fixed generator `P`, input `Q = [x]P` for uniform `x in [0, r)`.

`classical-group-work-v2` fixes one reference target operation at `C = 1000` scalar
AArch64 instruction equivalents, with the 5x52 field multiplication measured at
**179 instructions**, `S = M` and `A = 0`. The score is `log2(W / C)` where `W` is the
total instruction-equivalent work, including preprocessing, table construction,
lookups, bookkeeping, verification and output.

**Every number below is rounded up.** No displayed figure is a downward rounding of the
quantity it reports, and the declared bound is strictly greater than the computed one at
every step. The previous revision of this candidate failed review on exactly that point.

## 2. Construction

Let `m = 2^128`, `h = ceil(r/m) = 2^128 - 1`, so `m*h >= r`. Fix an interleaving
width `W = 2^16`, with `K_b = m/W = 2^112` and `K_g = ceil(h/W)`.

```text
Baby chains.  For c = 0 .. W-1, set B_c = [c]P. For k = 0 .. K_b-1:
              each chain computes the next point B_c <- B_c + [W]P with an
              incomplete affine addition whose denominator is x(B_c + [W]P) - x(B_c).
              All W denominators of a round are inverted with ONE Montgomery
              inversion, then W additions finish. Insert (x(B_c), c + k*W) into
              the search tree.

Giant chains.  Set D = [m]P. For c = 0 .. W-1, set G_c = Q - [c]D. For k = 0 .. K_g-1:
              each chain computes G_c <- G_c - [W]D. One batched inversion per round,
              as above. Look up x(G_c) in the tree.

Match.        A tree entry j found at giant index i = c + k*W gives
              G = Q - [i]D = [j]P, so x = i*m + j (mod r). Verify [x]P = Q.
```

Correctness. Write the secret as `x = a*m + b` with `0 <= b < m` and `0 <= a < h`.
Then `Q - [a]D = [x - a*m]P = [b]P`, and `b < m` is covered by the baby chains, so
the giant sequence meets a tree entry at `i = a`; every `i < h` is reached because the
giant chains together cover `i = 0 .. h-1`. The recovered `x = a*m + b` satisfies
`[x]P = Q` by construction and is re-verified.

The search is exhaustive over the `m` baby and `h` giant indices, so **success
probability is 1** for every valid `Q`. No randomness, retry, advice or cross-target
amortisation is used.

## 3. Exceptional cases, handled and charged

A reviewer of the previous revision correctly observed that incomplete affine addition
cannot advance the initial baby point `B_0 = [0]P = O`, and that the earlier code left
an infinity point unchanged, so that chain never produced indices `W, 2W, ...`. That is
fixed here, and the four cases below are handled explicitly in the chain-advance step.
Each is a property of one step, not of the asymptotic schedule, so each is charged.

| Case | Condition | Operation | Charge |
| --- | --- | --- | ---: |
| `FROM_INF` | current point is `O` | result is the increment itself, no field op | 10 |
| `DOUBLING` | `x_1 = x_2` and `P = inc` | `lambda = 3x^2 / 2y`, i.e. one extra squaring and one extra multiplication over the affine formula | 179 |
| `TO_INF` | `x_1 = x_2` and `P = -inc` | result is `O`, no field op | 10 |
| `SKIP` | index past the chain bound | no operation | 0 |

`FROM_INF` is not rare: chain `c = 0` begins at `[0]P = O`, so exactly one step per
baby-chain set and one per giant-chain set takes this branch, `2*W = 2^17` times in the
whole run. At 10 instructions each that is `2^20.3` instructions, which is `2^108`
times smaller than the total and does not appear at the declared precision.

`TO_INF` and `DOUBLING` occur only when a chain meets its own fixed increment, at
distance `m` or `W` in the group. For a random input the probability that a given step
is exceptional is `2^-(254)`-scale, and over the whole `2^129`-step run the expected
number of such events is `2^-125`. They are nonetheless charged, and the per-step
ledger in section 5 already carries the `x_1 = x_2` branch test in its control row, so
the bound covers them without a probability argument.

The boundary case `b = 0` is separate and deterministic: when the giant chain reaches
`O` at index `i = a`, no tree entry exists there, so that step reports `x = a*m`
directly and re-verifies it. `Q = O` is returned as `x = 0`. Both are exercised in
`experiment.py`, the `b = 0` case exhaustively over every scalar at toy size.

## 4. Why the batch is real, and what it costs

Each round processes `W = 2^16` **independent** chains, so the `W` denominators of a
round are simultaneously available and one inversion serves all of them. This is what
the first revision of this candidate got wrong.

Montgomery's trick for `W` inversions costs `3(W-1)` multiplications: `W-1` for the
running products, `W-1` for the backward recovery, and `W-1` for the final multiply
into each denominator. The single inversion itself costs `2*256 = 512` multiplications
by binary powering. Amortised over the round this is

```text
(3*(W-1) + 512) / W = (196605 + 512)/65536 = 197117/65536 = 3.007767...
3.007767... * 179 = 538.39024353... instructions per step
```

Charged as **538.40**, rounded up. This `3M` amortisation is the term the cost model's
own `6M` convention already budgets, and it is charged in full.

## 5. Per-step ledger

One group operation, incomplete affine, `S = M`, `A = 0`:

| Item | Instructions |
| --- | ---: |
| `x_3 = lambda^2 - x_1 - x_2` (one squaring) | 179 |
| `y_3 = lambda*(x_1 - x_3) - y_1` | 179 |
| multiplier applied to the batched inverse for `lambda` | 179 |
| Montgomery batch inversion, `3.007767` multiplications (section 4) | 538.40 |
| deterministic tree lookup | 676 |
| index update, chain-state branch (`O` / doubling / inverse / skip), loop control | 20 |
| **per-step total** | **1771.40** |
| **charged per step, rounded up** | **1772** |

The lookup row is a **worst-case** bound, not an expectation, which is the second point
the previous revision failed. The tree is a B-tree over the `m` keys, node fanout `64`,
sorted by the high 128 bits of `x`:

- height `ceil(log_65(2^128)) = ceil(21.25) = 22` levels;
- each level performs a binary search over 64 keys, 6 comparisons, each comparison
  loading an 8-byte prefix, comparing and branching: 3 instructions, so 18 per level,
  plus node offset and child-pointer arithmetic: 4, giving 22 per level;
- `22 * 22 = 484` for the descent;
- at most one leaf can hold keys sharing the queried prefix, and a leaf holds at most
  64 keys, so a worst-case linear scan of that leaf costs `64 * 3 = 192`;
- `484 + 192 = 676`.

No hashing, no probe-length distribution, no load-factor assumption and no expected
value enters this row. It is the maximum work of one lookup over all `2^128` keys.

`1772` is therefore an upper bound on the charged work of one step, and every step of
the search is one such operation: `m + h = 2^129 - 1` steps in total.

## 6. Total, memory and preprocessing

```text
per step  <= 1772                                    (rounded up from 1771.40)
steps      = m + h = 2^129 - 1
W_total   <= (2^129 - 1) * 1772 < 2^129 * 1772
log2(W_total) <= 129 + log2(1772) = 129 + 10.79135 = 139.79135
score      = log2(W_total / 1000) <= 139.79135 - 9.96578 = 129.82557
```

Submitted bound: **time_log2 = 129.83**, i.e. `12983 / 100`. The computed value
`129.82557` is below the declared `129.83`, so the declaration is the rounded-up form.

Memory: `m = 2^128` keys in a B-tree of `2^128` 48-byte records (32-byte `x`, 16-byte
index) with at most 2% node overhead is `2^129 * 48 = 2^134.5850` bytes, rounded up to
**2^134.59** bytes; the `W = 2^16` chain points at 96 bytes each add `2^22.6` bytes,
below the reported figure. Nonuniform advice is **0** bytes.

Preprocessing is the `m` baby steps:

```text
2^128 * 1772 < 2^128 * 2^10.79135 = 2^138.79135 instructions
log2 = 128 + 10.79135 - 9.96578 = 128.82557  =>  declared 128.83
```

It is **included** in `W_total` above, not deducted from it. The tree is built inside
the run from the public generator, and every insertion is charged.

## 7. Heuristics and limitations

- No cryptanalytic heuristic and no probability model: the construction is exhaustive
  and deterministic, with success probability 1, so the `>= 0.39` threshold is met
  without a random-walk or birthday premise and without a collision-dependency
  assumption.
- No lookup heuristic either: the 676-instruction row is a worst-case tree bound, so
  the total does not depend on a hash distribution or an expected probe count.
- The only modelling inputs are the board's own constants (179 instructions per 5x52
  multiplication, `S = M`, `A = 0`) plus the itemised control row of 20 instructions.
  Raising that row to 60 still gives `129.86`.
- No full-size execution is performed or claimed; the `2^128`-key tree is not
  instantiated. This is an analytic bound in the same style as the organizer control.
- This is a deterministic generic search costed under `classical-group-work-v2`. It is
  not a new cryptanalytic technique, not a best-known-attack estimate, and it does not
  contradict secp256k1's ~2^128 generic security level.
- The declared bound is weaker than the previous revision's `129.22`, and weaker than
  the organizer control's `129.000802`. That is the honest consequence of charging the
  batch inversion in full and replacing an expected probe count with a worst-case tree
  descent. I would rather report `129.83` than defend `129.22` on a lookup argument that
  does not hold.

## 8. Experimental evidence

`experiment.py` runs the same construction at toy scale on `y^2 = x^3 + 7`: it builds
the baby and giant chains with `W` interleaved chains and one batched inversion per
round, and recovers the scalar on every generated instance. It explicitly exercises the
four chain-advance cases of section 3, including chain `c = 0` starting at infinity, and
checks the batched inversion against independent per-element inversion.

This is evidence that the interleaving, the exceptional-case handling and the
batch-inversion schedule are correct as implemented. It is not evidence about
production cost, memory or timing, and no part of the declared score is derived from it.
