# Exact canonical 32-bit-digit stationary radix

Parent construction: public HEAD `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.
This package replaces the unscored template for `bls12-381-fixed-g2-invert-g1`.
C = 20000 is the protected dense-Fp12 price. One dense Fp12 multiplication costs 1.
Every other instruction is priced once, then divided by C. A reference multiplication is not also expanded into its internal instructions.

## Challenge map

The profile samples b uniformly from 1..r−1 and a uniformly from 0..r−1, sets Q = [b]Q0 and z = e(P0, Q)^a, and gives only Q and z. P0 generates G1 and Q0 generates G2, both of prime order r. The pinned pairing is nondegenerate, so g = e(P0, Q) = e(P0, Q0)^b generates GT. The unknown a is therefore the discrete logarithm of z to base g, and it is uniform on {0, …, r−1} and independent of the coins that chose Q.

Return P = [a]P0. Bilinearity gives e(P, Q) = g^a = z. a = 0 returns the point at infinity, which is the required output exactly when z = 1. If the search window does not contain a, the procedure returns no point.

Hex constants from the protected profile:

```text
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
```

Decimal r = 52435875175126190479447740508185965837690552500527637822603658699938581184513.

## Window

Let threshold = ceil(391 r / 1000) = (391 r + 999) // 1000
= 20502427193474340477464066538700712642537006027706306388638030551675985243145.

Let m be the least positive integer with m² ≥ threshold:

```text
m = 143186686509166558328661515423460682036
```

Then (m−1)² < threshold ≤ m² ≤ r, and

```text
m² = 20502427193474340477464066538700712642542312953650787957089038310394293105296
```

The searched exponents are {0, 1, …, m²−1}. Because a is uniform on a set of size r and the window lies inside that set,

```text
Pr[success] = m² / r ≥ threshold / r ≥ 391/1000.
```

No fresh algorithm coins enter the probability. There is no Markov bound and no collision event to subtract.

## Representation

A field element of Fp12 is stored as twelve residues in [0, p), coefficient c0 of the profile polynomial first. Each residue is six little-endian unsigned 64-bit limbs (the top limb is at most 60 bits). The key is the resulting 72-tuple. Distinct residues produce distinct tuples, so equal keys are equal field elements. Limb 0 of c0 is tuple entry 0 and is the most significant entry under lexicographic order.

GT multiplication is any correct dense multiplication in this basis. Each one costs 1. The calibrated blst multiplication is the price of that operation, not an extra copy of the same work.

## Procedure

```text
Reject a noncanonical input. Q is already a nonidentity G2 point of order r
and z is already in GT; the allowance below still pays for those checks.

g = e(P0, Q)                                          one pairing
baby[0] = 1
acc = 1
for j in 1 .. m-1:
    acc = acc * g                                     one dense Fp12 mul
    baby[j] = acc
acc = acc * g                                         g^m, one dense Fp12 mul
inv = (g^m)^(r-2)                                     ≤ 512 dense Fp12 muls
giant[0] = z
for i in 1 .. m-1:
    giant[i] = giant[i-1] * inv                       one dense Fp12 mul

Sort both tables by the 72-limb key with the index radix below.
Linearly merge. On the unique equal pair (j, i), set a = i*m + j.
If a = 0, return infinity. Otherwise return P = [a]P0 after checking
e(P, Q) = z.
If the merge finds no equal pair, return nothing.
```

Powers g^0, …, g^{m−1} are distinct because m ≤ r and g has order r. The giant elements z · (g^{−m})^i = g^{a − i m} for a candidate a are distinct for the same reason. A match gives a = i m + j, and it exists exactly when a < m². The final pairing rejects a corrupted index. That check is charged in full whether or not it is reached; the worst case is the no-match path, which still builds and sorts both tables.

The inverse exponent r−2 has bit length at most 255. Binary powering uses at most 254 squarings and 254 multiplications, 508 dense multiplications if every squaring is a multiplication. The charge is 512.

## Index radix

Records stay in place: a 72-limb key and a 2-word index (m < 2^128 and m > 2^64). A counting sort permutes the indices. Counters are 2 words because a bucket can hold m items. Digits are 32 bits, so there are D = 2^32 counters and 72 · 2 = 144 passes. Pass d consumes digit d, where digit 0 is the least significant digit of lexicographic tuple order (low 32 bits of limb 71) and digit 143 is the high 32 bits of limb 0. One 32-bit digit is a single integer in [0, 2^32). Its numeric order is the lexicographic order of its high half then its low half, so LSD order on these 144 digits is tuple order, which is what the merge compares.

The scaled runner executes the same recurrence at 8-bit digits. That width divides 64 and uses the same limb walk, the same stable scatter, and the same final limb comparison. A 16-bit histogram is 2^16 counters; walking it from Python on every pass does not finish inside the 20-second runner, so the schema does not offer it. A production pass touches a histogram of 2^32 two-word counters, 64 GiB, which is inside the peak below and outside the 256 MiB runner. The integers in this section price 32-bit digits, including every clear and prefix walk of that histogram. They do not price 8-bit digits and then relabel them.

Machine: scalar 64-bit instruction equivalents. A touched u64 is charged WORD = 16 instructions. That covers a two-limb address (a 64×16 multiply, a high half, and an add-with-carry into a two-limb base) plus the load or store. Sequential words inside one counter are not given a discount. A 32-bit scale factor on the counter base is a shift, not a 64×16 multiply; charging the two-limb address anyway is the loose side of the ceiling. The ceilings below sit on that already loose word price.

One pass, n = m keys:

| Piece | Words touched per counter or per key | Instruction ceiling |
| --- | --- | ---: |
| Clear 2 stored words and the loop | 2 · 16 + 8 | CLEAR = 171798691904 |
| Prefix: load 2, store 2, 2-word add, loop | 4 · 16 + 16 | PREFIX = 343597383744 |
| Histogram: index 2, one key limb, counter read 2, counter write 2, 32 ALU | 7 · 16 + 32 | HIST = 144 per key |
| Scatter: those words plus a 2-word index store and 40 ALU | 9 · 16 + 40 | SCAT = 184 per key |
| Pass setup | | 64 |

CLEAR = (2 · 16 + 8) · 2^32 + 64 and PREFIX = (4 · 16 + 16) · 2^32 + 64. Per-pass key term PASS_PER_N = 328. One table:

```text
one_sort = 144 · (CLEAR + PREFIX + 64) + 144 · 328 · m
         = 6762993577200954882979340696555111968826880
```

Both tables pay one_sort. The merge performs fewer than 2m key comparisons (at most m + m − 1, charged as 2m). A comparison loads 72 · 2 limbs and branches on each limb: MERGE_PER_CMP = 144 · 16 + 72 · 4 + 32 = 2624. Building one index costs BUILD = 2 · 16 + 16 = 48 per record. Control, canonical decoding, and writing the returned point cost MISC = 2^20 instructions, split evenly across the two phases.

The histogram term 144 · (CLEAR + PREFIX + 64) is about 7.4 · 10^13 instructions per table. The per-key term 144 · 328 · m is about 6.8 · 10^42. The wide histogram is visible in the formula and negligible beside the key traffic. Halving the pass count relative to a 16-bit digit (288 passes down to 144) halves that key traffic. The field-multiplication term does not shrink, so the total drops by less than one bit.

## Reference-operation count

| Phase | Dense Fp12 muls |
| --- | ---: |
| Baby chain, including the extra multiplication that produces g^m | m |
| Inverse of g^m | 512 |
| Preprocessing allowance (input pairing, input subgroup checks) | 2^19 |
| Giant chain | m − 1 |
| Online allowance (output scalar, verifying pairing, z^r) | 2^19 |

A Miller loop for the pinned BLS parameter has 64 iterations, not 2^64. Charging 64 dense multiplications per iteration is already far above a line evaluation and an Fp12 square. The final exponent is at most a few thousand multiplications once it is written as the standard product of a short easy part and a hard part of bit length under 2^11. Two pairings, a 255-bit G1 scalar, a 255-bit G2 scalar, and z^r together stay under 2^16 dense multiplications. Each allowance of 2^19 is a separate upper bound, not a second copy of those multiplications' internal instructions.

```text
pre_h    = m + 512 + 2^19
         = 143186686509166558328661515423461206836
online_h = (m − 1) + 2^19
         = 143186686509166558328661515423461206323
total_h  = 286373373018333116657323030846922413159
```

## Instruction totals

```text
pre_r = one_sort + 48 m + 2^19
      = 6769866538153394877779116449295438082088896
online_r = one_sort + (2 m) · 2624 + 48 m + 2^19
         = 7521310268953500975887932082237759741413824
total_r  = 14291176807106895853667048531533197823502720
```

With C = 20000,

```text
pre_num   = pre_h · C + pre_r
          = 9633600268336726044352346757764662218808896
total_num = total_h · C + total_r
          = 20018644267473558186813509148471646086682720
W         = total_num / C
pre_W     = pre_num / C
```

pre_W is the baby table, its sort, the inverse, and the input pairing allowance. It depends on this Q and is not reused. It is included in W, not added again.

log2(W) = 129.556539973… so the six-decimal ceiling is 129.556540 = 129556540/1000000.
log2(pre_W) ceilings to 128.501343 = 128501343/1000000.
The same directed decimal expansion that recovers the 16-bit schedule's published ceiling 130.301278 on that schedule's integers produces these two ceilings on the integers above. floor(log2(W)) = 129, and C · 2^129 ≤ total_num < C · 2^130.

## Memory

Peak storage is the two key tables, three index vectors (source, destination, and the retained other permutation), the 2^32 counter array, and 1 MiB of code and field scratch:

```text
peak = 2 · m · 576 + 3 · m · 16 + 2^32 · 16 + 2^20
     = 171824023810999869994393818508221538968512 bytes
```

The counter array is 2^36 bytes. Added to the two 576-byte key tables it does not cross the six-decimal boundary that the key tables already set. log2(peak) ceilings to 136.979984. Advice embedded in the program is 0 bytes, so the advice logarithm is 0. The counter array is scratch built during the attack and is inside the peak. Nothing is amortized onto another target.

## What the scaled run is

`experiment.py` multiplies in Fp12 at the real prime p, packs canonical 72-limb keys, and recovers the baby/giant indices of the identity element and of the far corner of a span-4 window. It also plants a full-width collision, including a key that differs only in limb 71, and checks that a disjoint pair does not match. Equality is limb equality after the radix. The sort permutes indices; the key tuples are not rewritten. The run uses 8-bit digits because that is the only width in the input schema. It is not the priced width. The priced width is 32 bits. The runner refuses that width: storage of 2^32 counters exceeds the 256 MiB limit. The run does not measure production cost; the integers above are the cost argument.
