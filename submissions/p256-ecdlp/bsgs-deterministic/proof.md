# p256-ecdlp: deterministic baby-step giant-step

Status: ready for independent exploratory review; no qualification or emitted
score is claimed by this file. This is a generic accounting control under
`classical-group-work-v2`, not a cryptanalytic advance. `claim.json` binds
`tracks/p256-ecdlp/target-profile.json` and `tracks/classical-group-work-v2.json`
by SHA-256. No auxiliary input, advice, oracle, randomness or quantum step
is used. The success argument does not contain a heuristic.

P-256 has prime order r (cofactor 1) and a = -3. The algorithm does not use
a = -3, endomorphisms, or the negation map. Every bound below is an inequality
in the group law or an instruction tally, and every inequality is in the
direction that increases the charged cost.

## Parameters

Profile constants p, r, P, a, b are those in the protected target profile.
Write

```
m = ceil(sqrt(r)) = 340282366881324382213159967663442624512
```

Checked integers:

```
(m-1)^2 = 115792089210356248762697446949407573529389627945510159738472418292043737989121
r       = 115792089210356248762697446949407573529996955224135760342422259061068512044369
m^2     = 115792089210356248762697446949407573530070192679272808502898738227370623238144
(m-1)*m = 115792089210356248762697446949407573529729910312391484120685578259707180613632
r-1     = 115792089210356248762697446949407573529996955224135760342422259061068512044368
```

So (m-1)^2 < r ≤ m^2, hence m = ceil(sqrt(r)). Also (m-1)*m ≤ r-1 < m*m, so
floor((r-1)/m) = m-1. Every scalar x in 0..r-1 has a unique writing
x = i + j*m with integers i, j satisfying 0 ≤ i < m and 0 ≤ j < m.
Since r > 2^255 and m < 2^128, one has 2m < r, so the m baby scalars are
distinct modulo r.

## Algorithm

1. Decode Q under the profile encoding. Evaluate the curve equation (a few
   field operations, charged in the reserve below). Compute [r]Q by the
   binary method and require the result to be the identity. If Q is the
   identity, output x = 0. The challenge distribution is exactly the set of
   Q = [x]P, so this check accepts every legal input. An illegal input is
   outside the probability space; the algorithm then stops with failure.

2. Baby steps. Compute B_i = [i]P for i = 0..m-1.
   B_0 is the identity. B_1 = P. For i ≥ 2, B_i = B_{i-1} + P.
   If the two summands are equal, use the doubling formula; if they are
   opposite, the sum is the identity. Both formulas are covered by the
   per-step instruction bound below, which charges the more expensive one
   on every index.

3. Uniqueness of baby x-coordinates. Suppose 0 ≤ i < k ≤ m-1 and [i]P and
   [k]P are finite with the same affine x. Then [k]P = ±[i]P. The plus case
   gives k ≡ i (mod r), impossible. The minus case gives k + i ≡ 0 (mod r).
   But 0 < k + i < 2m < r, so k + i = 0, impossible. Therefore the finite
   baby points have distinct x-coordinates, and none is the identity.
   The table key is the 32-byte affine x, together with the 32-byte y and
   the 16-byte index i stored beside it. A later lookup accepts a hit only
   when both coordinates match.

4. Sort. LSD radix sort of the m finite records on the 32-byte x key, 8 bits
   per pass, 32 passes, two buffers. The identity is stored aside as the
   i = 0 record and is not part of the radix key space. Counting-sort
   histograms are 256 counters per pass.

5. Giant generator. Compute G = [m]P by the binary method and replace G by
   -G (one field negation). This is preprocessing: it depends only on P.

6. Giant steps. Set S_0 = Q. For j = 0, 1, ..., m-1: if S_j is the identity,
   the match is i = 0; otherwise binary-search the sorted baby table on x
   and accept the row only when y matches. On a hit at index i, set
   x = i + j*m, check [x]P = Q by one binary scalar multiplication, and
   output x. Otherwise set S_{j+1} = S_j + (-G), using doubling when the
   summands are equal. The loop bound j < m is exactly the range proved
   above, so a legal input produces a hit. The cost bound charges all m
   iterations on every input, which dominates an implementation that returns
   at the first hit.

7. No randomness is drawn. No table depends on a hidden precomputation:
   the baby table is built in this attack and its construction is charged
   in full against this single target. It is not amortized over other
   targets and it is not treated as advice.

## Why the output is the discrete log

On a hit, S_j = [i]P. The recurrence is S_j = Q + [j](-G) = Q - [j*m]P.
Thus Q = [i + j*m]P. The writing i + j*m is the unique scalar in 0..r-1
with that shape, and the cofactor is 1, so it is the unique discrete log.
The final [x]P = Q check rejects a corrupted index. There is no collision
probability, no random-mapping approximation, and no cycle analysis.

## Instruction ledger (scalar AArch64, 64-bit GPRs)

Field elements are four 64-bit Montgomery limbs. For this prime,
-p^{-1} ≡ 1 (mod 2^64), so the CIOS multiply by n0' is a move, not a MUL.
The tally nevertheless charges a generic 4-by-4 CIOS round (four MUL and
four UMULH per round). That is an upper bound on the shift-and-add
reduction this prime actually allows. Squares are charged as general
multiplications.

Field multiplication or square, charged M = 384 instructions:

| Item | Instructions |
| --- | ---: |
| 16 MUL + 16 UMULH schoolbook | 32 |
| 16 partial products, ≤ 8 ADC each into the 8-limb accumulator | 128 |
| 4 CIOS rounds, each ≤ 4 MUL + 4 UMULH + 8 ADC | 64 |
| Final conditional subtract: 4 SUB, 4 CSEL, flags | 16 |
| 8 loads + 4 stores | 12 |
| Address arithmetic and moves | 48 |
| Subtotal | 300 |
| Charged | 384 |

Field addition or subtraction, charged A = 48 instructions: 4 ADDS, 4 SUBS
against p, 4 CSEL, 8 loads, 4 stores, address arithmetic. The executed
prefix is at most 32; the charge is 48.

Inversion, charged I = 154880 instructions. The exponent p-2 has bit length
256 and Hamming weight 128, so the binary method performs 255 squarings and
127 multiplications, 382 field multiplications, plus a loop. 382*M = 146688.
Loop control, bit tests and address arithmetic are charged an additional
8192. Total I = 154880. A shorter addition chain would only reduce this
number; the binary method is the bound.

One batched affine point operation (addition or doubling), batch size
B = 1024, charged S = 3384 instructions:

- Simultaneous inversion of a batch costs 1 inversion and 3(B-1) field
  multiplications. Per element that is at most I/B + 3M. A short final
  batch is repaired below by charging two extra full inversions, so every
  element may be charged ceil(I/B) + 3M.
- The addition formula uses at most 2 further multiplications. The doubling
  formula uses at most 4 (x^2, the multiplication by the inverted
  denominator, λ^2, and λ*(x-x3); multiplication by 3 is an addition and a
  shift, charged as additions). Charge 4M for either formula.
- Additions and subtractions in either formula, including the *3 and the
  *2, are at most 10. Charge 10A.
- Branches, equality tests, infinity tests and address arithmetic: 64.

ceil(I/1024) = ceil(154880/1024) = 152.
S = 7*384 + 10*48 + 152 + 64 = 2688 + 480 + 152 + 64 = 3384.

The *3M term charges 3B multiplications per full batch against an actual
3(B-1), so each full batch is overcharged by 3M. That surplus is not used
to pay for anything else.

Radix-sort move, one record, one 8-bit pass, charged 64 instructions:
12 loads and 12 stores of a 96-byte padded record (32-byte x, 32-byte y,
16-byte index, 16 bytes of pad), plus byte extract, histogram update and
destination-pointer update. The executed prefix is at most 48. Thirty-two
passes cost 2048 instructions per record. Histogram prefix sums are
32*256 counter operations and are absorbed by the per-record pad at this
table size (the pad alone is 16 instructions * 32 passes * m, far above
32*256*20).

Binary search, charged 4096 instructions per giant step. The table length
is m < 2^128, so a binary search performs at most 128 comparisons. Each
comparison is charged 32 instructions (index arithmetic, four limb loads,
four compares, branch) against an executed prefix of about 14. Final y
comparison on the matching x is inside the same 4096.

Unbatched binary scalar multiplication, charged 80363520 instructions:
512 point operations (256 doubles and 256 adds; the real bit length is 256
and the real Hamming weight is far below 256) at I + 4M + 10A + 64 each.
One such operation is 154880 + 1536 + 480 + 64 = 156960 instructions.
512*156960 = 80363520. This pays for a curve-validation [r]Q, for [m]P,
for the output check [x]P, and for one spare scalar multiplication.

## Work

Every one of the m baby indices is charged a full batched point operation
(the identity and the copy of P are strictly cheaper) and a full radix sort:

```
baby point operations  m * 3384 = 1151515529526401709409333330573089841348608
baby radix sort        m * 2048 =  696898287372952334772551613774730495000576
```

Every one of the m giant indices is charged a full batched point operation
and a full binary search:

```
giant point operations m * 3384 = 1151515529526401709409333330573089841348608
giant lookups          m * 4096 = 1393796574745904669545103227549460990001152
```

Fixed terms, charged on every input:

```
4 scalar multiplications             321454080
2 extra inversions (short batches)      309760
decode, equation check, output reserve   65536
```

The two extra inversions cover the short-batch deficit. Partitioning m
into batches of 1024, the formula "ceil(I/B) per element" under-counts the
true inversion count by at most one inversion per side, and by zero when
the side divides evenly. It never under-counts field multiplications.

Total instructions

```
T = 4393725921171660423136321502470371489528320
```

The six lines sum to T. Division by the fixed C = 1000 is the only
conversion into target-group-operation equivalents. No point addition is
priced as one unit, and the secp256k1 6M reference is not used.

log2(T/1000) = 131.690640578... . The submitted score is the ceiling at
six decimal places,

```
time_log2 = 131690641/1000000 = 131.690641
```

so 2^time_log2 > T/1000. The ratio 1000 * 2^time_log2 / T is approximately
1.000000292.

Preprocessing is the baby table, its radix sort, [m]P, and one short-batch
inversion. It does not depend on Q and it is included in T, not added again:

```
P_instr = 1848413816899354044181884944347820416867584
preprocessing_log2 = 130.441484
```

which is log2(P_instr/1000) rounded up at six decimal places.

## Memory and advice

Two radix-sort buffers hold m records of 96 bytes each. Code, field
constants, the identity record, the giant point and scalar scratch are
charged a further 2^20 bytes.

```
peak bytes = 2 * m * 96 + 2^20
           = 65334214441214281384926713791380984954880
memory_log2_bytes = 135.584963
```

again rounded up at six places. Nonuniform advice is empty: nothing is
embedded in the program, and the baby table's construction is the
preprocessing above. nonuniform_advice_log2_bytes = 0.

## Success

On every Q = [x]P with x uniform in 0..r-1 the algorithm returns that x.
The algorithm draws no coins. Success probability on the required
probability space is exactly 1, which is at least 0.39. There is no
heuristic entry.

## Evidence and limitations

`experiment.py` runs this recurrence on the prime-order curve
y^2 = x^3 + x + 8 over F_53, group order 61, generator (1, 13), and checks
every scalar in 0..60, including the identity. It is evidence only that the
recurrence recovers the discrete log on that group. It does not measure
AArch64 instructions or the production memory bound. The production bound
is the ledger above. No full-size attack is executed. The exponential
memory is required for this algorithm. This package is a deterministic
generic control, not a best-known-attack estimate and not a practical break.
