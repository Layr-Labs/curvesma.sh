# secp256k1-ecdlp: interleaved-chain BSGS with metered batch inversion

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

## 2. Construction

Let `m = 2^128`, `h = ceil(r/m) = 2^128 - 1`, so `m*h >= r`. Fix an interleaving
width `W = 2^16`, with `K_b = m/W = 2^112` and `K_g = ceil(h/W)`.

```text
Baby chains.  For c = 0 .. W-1, set B_c = [c]P. For k = 0 .. K_b-1:
              each chain computes the next point B_c <- B_c + [W]P with an
              incomplete affine addition whose denominator is x(B_c + [W]P) - x(B_c).
              All W denominators of a round are inverted with ONE Montgomery
              inversion, then W additions finish. Insert (x(B_c), c + k*W) into
              the table.

Giant chains.  Set D = [m]P. For c = 0 .. W-1, set G_c = Q - [c]D. For k = 0 .. K_g-1:
              each chain computes G_c <- G_c - [W]D. One batched inversion per round,
              as above. Look up x(G_c) in the table.

Match.        A table entry j found at giant index i = c + k*W gives
              G = Q - [i]D = [j]P, so x = i*m + j (mod r). Verify [x]P = Q.
```

Correctness. Write the secret as `x = a*m + b` with `0 <= b < m` and `0 <= a < h`.
Then `Q - [a]D = [x - a*m]P = [b]P`, and `b < m` is covered by the baby chains, so
the giant sequence meets a table entry at `i = a`; every `i < h` is reached because the
giant chains together cover `i = 0 .. h-1`. The recovered `x = a*m + b` satisfies
`[x]P = Q` by construction and is re-verified.

Two boundary cases are handled explicitly rather than by the table. `Q = infinity` is
returned as `x = 0`. When `b = 0` the giant chain reaches the point at infinity at
`i = a` and no table entry exists there, so that step reports `x = a*m` directly and the
result is re-verified the same way. Both cases are checked in `experiment.py`, the
second one exhaustively over all scalars at toy size.

The search is exhaustive over the
`m` baby and `h` giant indices, so **success probability is 1** for every valid `Q`;
No randomness, retry, advice or
cross-target amortisation is used.

## 3. Why the batch is real, and what it costs

The previous version of this candidate was rejected for assuming that a sequential
affine chain supplies independent denominators. It does not: `B_{k+1}` depends on
`B_k`. The interleaving above fixes exactly that. Each round processes `W = 2^16`
**independent** chains, so the `W` denominators of a round are simultaneously available
and one inversion serves all of them.

Montgomery's trick for `W` inversions costs `3(W-1)` multiplications: `W-1` for the
running products, `W-1` for the backward recovery, and `W-1` for the final multiply
into each denominator. The single inversion itself costs `2*256 = 512`
multiplications by binary powering. Amortised over the round this is

```text
(3*(W-1) + 512) / W = (196605 + 512)/65536 = 3.0080 multiplications per step
                    = 3.0080 * 179 = 538.4 instructions per step
```

This `3M` amortisation is the term the cost model's own `6M` convention already
budgets, and it is charged in full here. Enlarging `W` does not remove it: the
`3(W-1)` product and recovery work is intrinsic to Montgomery inversion, and only the
single `512` term shrinks. That is why the per-step figure below is `~6M` and not `3M`.

## 4. Per-step ledger

One group operation, incomplete affine, `S = M`, `A = 0`:

| Item | Instructions |
| --- | ---: |
| `x_3 = lambda^2 - x_1 - x_2` (one squaring) | 179 |
| `y_3 = lambda*(x_1 - x_3) - y_1` | 179 |
| multiplier applied to the batched inverse for `lambda` | 179 |
| Montgomery batch inversion, `3.0080` multiplications (section 3) | 538 |
| hash, probe, store or compare, address arithmetic, index update, loop control | 85 |
| **per-step total** | **1160** |

The 85-instruction row is itemised: the table key uses bits of `x` directly, so hashing
is a two-word extraction and shift (6); slot address and linear probe (18); a store or a
key comparison plus index read (22); chain index `c + k*W`, scalar accumulator update and
loop bookkeeping (24); sign/edge-case branch on `x_1 = x_2` (15). The table is
open-addressed with load factor `0.5` and linear probing, whose expected probe count at
that load factor is below `1.5`; the row above budgets 22 instructions for two probes.

`1160` is an upper bound on the charged work of one step, and every step of the search is
one such operation: `m + h = 2^129 - 1` steps in total.

## 5. Total, memory and preprocessing

```text
steps    = m + h = 2^129 - 1
per step <= 1160
W_total  <= (2^129 - 1) * 1160 < 2^139.1799
score    = log2(W_total / 1000) < 139.1799 - 9.9658 = 129.2141
```

Submitted bound: **time_log2 = 129.22**, i.e. `6461 / 50`.

Memory: `m = 2^128` entries in `2^129` slots at 48 bytes (32-byte `x`, 16-byte index)
is `2^129 * 48 = 2^134.5850` bytes, plus the `W = 2^16` chain points at 96 bytes each
(`2^22.6` bytes, below the reported figure). Reported peak is **2^134.59** bytes;
nonuniform advice is **0** bytes.

Preprocessing is the `m` baby steps, `2^128 * 1160 < 2^137.1799` instructions, i.e.
`log2 = 128.21` target-group-operation equivalents. It is **included** in `W_total`
above, not deducted from it. The table is built inside the run from the public
generator, and every insertion is charged.

## 6. Heuristics and limitations

- No cryptanalytic heuristic and no probability model: the construction is exhaustive
  and deterministic, with success probability 1, so the `>= 0.39` threshold is met
  without a random-walk or birthday premise and without a collision-dependency
  assumption.
- The only modelling inputs are the board's own constants (179 instructions per 5x52
  multiplication, `S = M`, `A = 0`) plus the itemised 85-instruction bookkeeping row.
  Raising that row to 200 instructions per step still gives `129.38`.
- No full-size execution is performed or claimed; the `2^128`-entry table is not
  instantiated. This is an analytic bound in the same style as the organizer control.
- This is a deterministic generic search costed under `classical-group-work-v2`. It is
  not a new cryptanalytic technique, not a best-known-attack estimate, and it does not
  contradict secp256k1's ~2^128 generic security level.
- The construction doubles the number of charged steps relative to a single-chain
  birthday search but needs no cycle detection, no restart schedule and no probability
  argument; that trade is deliberate and is what makes the accounting unconditional.

## 7. Experimental evidence

`experiment.py` runs the same interleaved-chain construction at toy scale: it builds the
baby chains and giant chains with `W` interleaved chains and one batched inversion per
round, and recovers the scalar on every generated instance. It also checks the batched
inversion against independent per-element inversion. This is evidence that the
interleaving and the batch-inversion schedule are correct as implemented; it is not
evidence about production cost, memory or timing.
