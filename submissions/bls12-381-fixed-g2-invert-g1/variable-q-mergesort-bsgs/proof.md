# Variable-Q mergesort baby-step giant-step

Classical cold attack on one challenge. `Q` is the supplied nonidentity G2 point. `z` is the supplied GT element. No auxiliary powers, no randomness, no advice table built offline from other targets.

Let `p` and `r` be the protected field modulus and subgroup order. Both are fixed by `target-profile.json`. `r` is an odd prime, `bit_length(r) = 255`, and the Hamming weight of `r` is 134. Set `L = 127`, `m = 2^127` and `h = ceil(r / m)`. The integers are

```text
r = 52435875175126190479447740508185965837690552500527637822603658699938581184513
m = 170141183460469231731687303715884105728
h = 308190375243917312857644309053803048971
```

`h * m >= r`, and `h * m - r = 58830589151200628854506063105544421375 < m`.

## Procedure

Decode `Q` and `z` under the protected encoding. Require `Q` on the G2 curve, `Q` not infinity, and `[r]Q = infinity`. Require `z != 0` and `z^r = 1`. Otherwise reject. The challenge sampler never produces those rejects.

Compute `g = e(P0, Q)` with the pinned pairing: Miller loop of the protected ate count `x = 0xd201000000010000`, then exponentiation by `e = (p^12 - 1) / r`. The base `g` lies in GT. It is not the identity: `Q` has order `r`, `P0` generates G1, and the pairing is non-degenerate. Because `|GT| = r` is prime, `g` generates GT.

Build the baby table by `m` dense Fp12 multiplications. Start at `acc = 1`. For `j = 0 .. m-1`, store the record `(canonical(acc), j)` and set `acc = acc * g`. The stored keys are `g^0, ..., g^(m-1)`, all distinct because `m < r`. The final `acc` is `g^m`.

Set `D = (g^m)^(r-1)` by binary exponentiation, squarings charged as dense multiplications. Since `(g^m)^r = 1`,

```text
(g^m) * (g^m)^(r-1) = (g^m)^r = 1,
```

so `D = g^(-m)`. The exponent is `r-1`, which has bit length 255 and Hamming weight 133: 254 squarings and 132 multiplications, 386 dense multiplications.

Sort the `m` records by the canonical key with bottom-up mergesort, `L` passes, two buffers. The key is the 12 GT coefficients in degree order. Each coefficient is an integer in `[0, p)`, written as six little-endian u64 limbs, and limbs inside a coefficient are compared most-significant first. The key is `K = 72` limbs. The record appends `j` as four u64 limbs, so `W = 76` limbs.

Set `gamma = z`. For `i = 0 .. h-1`, binary-search `gamma` in the sorted table. On a hit at index `j`, set `a = i * m + j` and stop. Otherwise set `gamma = gamma * D`. Then compute `P = [a] P0`, require `[r] P = infinity` and `e(P, Q) = z`, and return `P`.

## Coverage

`z = g^x` for a unique `x` in `[0, r)`. Write `x = i * m + j` with `0 <= j < m` and `0 <= i < h`. This writing exists because `h = ceil(r / m)`, and it is unique because `j` is the remainder modulo `m`. Then

```text
z * D^i = g^x * g^(-m * i) = g^j.
```

The baby table contains that key, so the search returns this `i` and `j`, hence returns `a = x`. Bilinearity gives `e([a] P0, Q) = e(P0, Q)^a = g^a = z`. The identity `z = 1` forces `a = 0` and `P` infinity, which the profile accepts exactly then.

Every valid challenge input is recovered. There are no algorithm coins. Success probability on the uniform challenge is 1, which is at least 0.39. The work bound below is worst-case, not an expectation.

## Reference price

`C = 20000`. One dense Fp12 multiplication, the pinned blst tower product, costs one target-group operation. Its internal instructions are not charged again. Every other AArch64 instruction costs `1/C`. For `H` dense multiplications and `R` other instructions, `W = H + R / C = (H * C + R) / C`.

## Instruction schedule

Limb comparison, at most 16 instructions: two loads, a compare, a conditional branch, index advance, trip test, loop branch, two address formations, and seven reserved moves or spills. Across 72 limbs plus a 64-instruction epilogue,

```text
B = 16 * 72 + 64 = 1216.
```

Record move, at most 16 instructions per limb (load, store, two pointer updates, counter, branch, and reserved spills) plus a 64-instruction epilogue,

```text
D = 16 * 76 + 64 = 1280.
```

A mergesort pass on a power-of-two length does fewer than `m` key comparisons and writes each record once. Charging `m` comparisons and `2 * m` record moves is an upper bound. Loop overhead is 8 instructions per record, and each pass has a 4096-instruction setup. Sorting therefore costs

```text
sort = m * L * (B + 2 * D + 8) + 4096 * L
     = 81764408253230777754833504172134992922013696.
```

Writing the unsorted baby records costs `baby = m * D = 217780714829400616616559748756331655331840`.

Binary search uses at most `L + 1 = 128` probes. Each probe is one comparison plus 16 instructions of index arithmetic. Returning the record costs `D + 64`. One search costs

```text
look = 128 * (B + 16) + D + 64 = 159040
giant = h * look = 49014597278792609436879750911916836908347840.
```

The search still charges all `h` probes when it stops early.

G1 arithmetic is not the reference operation of this track. One affine addition or doubling is charged at 8000 instructions, four times the published G1 price of 2000. A 255-bit binary scalar multiplication is at most `2 * 255` point operations. The attack performs at most four of them (`[a] P0`, `[r] P`, and two conditional infinity cleanups):

```text
g1 = 4 * 510 * 8000 = 16320000.
```

Exponent scanning for `D` and for the input check `z^r` costs `64` instructions per bit of a 255-bit exponent:

```text
exp_ctrl = 2 * 64 * 255 = 32640.
```

Fixed decode, I/O and control sit in `setup = 2^20 = 1048576`.

## Pairing and subgroup checks

`x` has bit length 64 and Hamming weight 6. The Miller loop squares once per bit after the leading bit (63 squarings), multiplies once per extra set bit (5 multiplications), and applies one line accumulation per bit plus the initial line (65). Each line accumulation is charged as a whole dense Fp12 multiplication.

`e = (p^12 - 1) / r` has bit length 4314 and Hamming weight 2124. Binary exponentiation uses 4313 squarings and 2123 multiplications. Squarings are charged as dense multiplications. That exponentiation is a correct final exponent: the Miller value to the power `(p^12 - 1) / r` is the pinned pairing.

One pairing is therefore at most

```text
63 + 5 + 65 + 4313 + 2123 = 6569
```

dense multiplications. The attack computes two pairings, `g = e(P0, Q)` and the output check `e(P, Q)`.

The input check `z^r = 1` is one more binary exponentiation: 254 squarings and 133 multiplications, 387 dense multiplications. The same price covers `[r] Q = infinity` only on the G1/G2 side through `g1`; the GT check is this exponentiation.

## Totals

Dense multiplications:

```text
H = m + h + 386 + 2 * 6569 + 387
  = 478331558704386544589331612769687168610
```

`m` builds the baby table and `g^m`. `h` is the giant walk. `386` builds `D`. `387` checks `z^r`.

Other instructions:

```text
R = sort + baby + giant + g1 + setup + exp_ctrl
  = 130996786246852787808329814832808161503094592
```

```text
H * C + R = 140563417420940518700116447088201904875294592
W = (H * C + R) / C
log2(W) = 132.368344967...
```

The submitted score is this logarithm rounded upward to six decimal places: `132368345/1000000 = 132.368345`. Hence `W <= 2^132.368345`.

Preprocessing is the first pairing, the baby table, the construction of `D`, and the sort. It excludes the giant walk, the second pairing, `z^r`, and G1 output.

```text
H_pre = m + 386 + 6569
      = 170141183460469231731687303715884112683
R_pre = sort + baby + setup + 64 * 255
      = 81982188968060178371450063920891324578410432
```

`log2(H_pre + R_pre / C) = 131.649178560...`, rounded upward to `131.649179`. Preprocessing is inside `W`.

Peak storage is the two `m`-record buffers plus a `2^20`-byte reserve for code, field scratch and the live `gamma`:

```text
bytes = 2 * m * 76 * 8 + 2^20
      = 206891679087930585785731761318515073613824
log2(bytes) = 137.247927513...
```

rounded upward to `137.247928`. Advice is empty: `log2(max(1, 0)) = 0`.

No step is charged both as a reference multiplication and as its internal instructions. Public profile constants are not charged. The bound does not amortise the table across future targets.

## Toy evidence

`experiment.py` runs the same join, sort then binary search, in the multiplicative group of a caller-supplied prime field, and in a fixed subgroup of order 29. The giant-step multiplier is `h^(order - 1)`. The order-29 witness also requires that `h^(order - 2)` does not invert `h`. The toy order is at most 5002. It is not a production execution.
