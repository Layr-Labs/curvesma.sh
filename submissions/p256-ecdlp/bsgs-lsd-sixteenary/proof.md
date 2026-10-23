# Odd-lift BSGS with global LSD distribution and implicit sixteen-ary lookup

Status: ready for exploratory review. No qualification or emitted score is assumed.
Normative inputs and outputs are the protected P-256 profile. `claim.json` binds
`tracks/p256-ecdlp/target-profile.json` and `tracks/classical-group-work-v2.json`.
One classical cold attack. No quantum step, no uncharged table, no restricted
challenge distribution.

The score unit is the track constant C = 1000 scalar AArch64 instruction
equivalents. P-256 point addition is not charged as one unit. Total work is
W = T/C, where T is an instruction upper bound derived below. The submitted
time_log2 is 1329/10 because W < 2^(1329/10).

## Parameters

Let p and r be the profile field modulus and subgroup order, and let P be the
profile generator. Both p and r are odd, the cofactor is 1, and a = -3.
Write M = 2^128, B = 2^127, G = 2^128 and block length L = 2^12.
Field elements are four little-endian 64-bit limbs. Points used in the
search chain are Jacobian (X : Y : Z) in the Montgomery domain, except the
stored baby records and the fixed giant stride, which are affine Montgomery
residues. The Montgomery map x ↦ xR mod p is a bijection on the field, so
equality of residues is equality of integers in [0, p).

## Odd-lift coverage

The identity Q = ∞ is the scalar 0 and is returned immediately.

For x in 1..r-1 let z = x when x is odd and z = x+r when x is even.
Then z is odd, 1 ≤ z ≤ 2r-1, and [z]P = [x]P because rP = ∞.

Every odd positive integer has one representation
z = (2i+1) M ± (2j+1) with i ≥ 0 and 0 ≤ j < B.
The two windows for a fixed i are the odd integers in
[2i M + 1, (2i+1) M - 1] and [(2i+1) M + 1, (2i+2) M - 1].
Consecutive i tile the odd positives. The last odd covered by i < G is
(2G - 1) M - 1 = 2^257 - 2^128 - 1.

This reaches every admissible z. The profile order satisfies
r ≤ 2^256 - 2^224 + 2^192 - 1,
because its top 32 bits are 0xffffffff, the next 32 bits are 0, and every
lower limb is at most 2^32-1. Then
2^256 - r ≥ 2^224 - 2^192 + 1 > 2^127,
so r < 2^256 - 2^127 and 2r - 1 < 2^257 - 2^128 - 1.
Every odd representative is therefore represented by some i < G and some
sign and baby index. Windows are disjoint, so the integer representation
is unique. A group match can still occur for a congruent integer; reduction
modulo r returns the same scalar.

No two baby scalars collide on x. Baby scalars lie in {1, 3, ..., 2^128-1},
hence are smaller than r. If (2j1+1) ≡ (2j2+1) (mod r) then j1 = j2.
If (2j1+1) ≡ -(2j2+1) (mod r) then their sum is divisible by r. The sum is
even, positive and strictly less than 2^129 < 2r, while r is odd, which is
impossible. Distinct babies have distinct x-coordinates. Negation does not
change x, and the stored y selects the sign.

## Algorithm

Reject a non-canonical or off-curve Q, and reject a Q with [r]Q ≠ ∞.
If Q = ∞, return 0.

Baby table. Double P to 2P and convert 2P to affine Montgomery form.
Starting from Jacobian P, add the affine point 2P repeatedly to obtain
Jacobian points (2j+1)P for j = 0..B-1. Batch-invert their Z coordinates
and write affine records (x, y, 2j+1), 80 bytes each: 32 + 32 + 16.
The Jacobian array and the prefix-product array are released before the
distribution buffer is allocated.

Distribution. Sort the B records by the 32-byte Montgomery x with a stable
global LSD radix-256 pass. Byte positions run from the least significant
byte to the most significant byte. Each pass builds one 256-bin histogram
of the whole table, the exclusive prefix of that histogram, and a stable
scatter into a second buffer of B records. The destination index is global,
so the destination page index//32 may differ from the source page. A pass
whose byte is constant on every record is not rewritten. There is no
recursion, no per-bucket planted directory, no merge heap and no bounded
port queue. The only auxiliary state is the 256 counters, reused each pass,
and the ping-pong of two record arrays.

Giant search. Precompute the affine stride S = 2^129 P and the start
R = Q - 2^128 P in Jacobian coordinates. The G successive points
R, R-S, R-2S, ... are exactly Q - (2i+1) 2^128 P. They are produced in
blocks of L. Each block is L mixed Jacobian-plus-affine additions, one
batch inversion of the block Z values, and the affine maps
x = X Z^{-2}, y = Y Z^{-3}.

Lookup is an implicit sixteen-ary descent on the sorted array. While the
live index interval has length n ≥ 16, split it into 16 contiguous parts
whose lengths differ by at most one. The part starts are
origin + k·width + min(k, remainder). Binary-search the part index in
0..15 by comparing the query with the record at the part start. No tree
node is stored. When the interval has length at most 15, finish with a
binary lower bound and test equality of x. On an x hit, compare y with the
baby y and with its negation. The sign selects
z = (2i+1) M ± (2j+1). Return z mod r after checking [z mod r]P = Q.

The true odd representative appears as some giant point, so a hit exists.
Any hit is a scalar of Q. The loops are charged in full even though the
implementation may return at the first hit. No randomness and no restart.

Success probability on the uniform profile challenge is 1.

## Instruction kernels

AArch64 general-purpose instructions only. Loads, stores, arithmetic,
addressing and branches count. Cache latency and issue width do not.
Thirty-one registers are available. The kernels below are straight-line
enough to keep the stated operands in registers.

Field addition or subtraction, at most A = 32 instructions. Four-limb
add-with-carry or sub-with-borrow (4), comparison against p (4), four
conditional selects, two `ldp` pairs per operand and one `stp` pair, plus
a fixed control margin that fills the cap.

Field multiplication, at most M = 768 instructions, Montgomery CIOS with
public n0 = -p^{-1} mod 2^64. That constant is profile arithmetic, not
searched advice. One CIOS inner step `(C, s) = s + u·v + C` is

```text
mul; umulh; adds; adc; adds; adc; mov
```

which is 7 instructions. Two loops of 4 steps and 4 outer iterations are
224 instructions. Limb loads that are not already resident, the Montgomery
multiplier `m = t0·n0`, the shift-down of the accumulator and one
conditional subtraction of p are at most 256 further instructions.
The kernel cap is 768. An inner step padded from 7 instructions to 16
is 4·2·4·16 = 512 instructions, and 512+256 = 768 still fits. Squaring
uses the same cap. S ≤ M.

Inversion of a nonzero residue by square-and-multiply on p-2. The exponent
is shorter than 256 bits, so the chain is at most 256 squares and 256
multiplications: I = 512 M = 393216 instructions.

Mixed Jacobian-plus-affine addition, the EFD `madd-2007-bl` schedule:
Z1², X2·Z1², Z1·Z1², Y2·(Z1·Z1²), H², H·(4H²), X1·(4H²), r²,
Y1·J, r·(V-X3), (Z1+H)². That is 7M+4S ≤ 11M, charged as 12M.
Limb additions, doublings by 2 and 4, and the infinity/opposite cases are
at most 16 field additions. One mixed addition is
12 M + 16 A = 9728 instructions. The same cap covers doubling.

Doubling with a = -3, EFD `dbl-2001-b`: 2M+5S ≤ 7M, charged as 8M+16A
= 6656 instructions. Used only for the power-of-two setup points.

One lookup probe, at most 32 instructions. The query x occupies four
registers. One probe forms a record address (`add`, `lsl` or `madd` by the
constant 80) and loads four limbs of the stored x, then compares them and
branches. Cap 32 against an enumerated body of about 15. The same cap
covers a pivot probe and a final lower-bound probe.

Sixteen-ary control on this table. Start with n_0 = 2^127. While n ≥ 16,
the next interval has length ceil(n/16). Here 2^127 is divisible by
16^31 = 2^124, so the lengths are n_i = 2^{127-4i}. The body runs for
i = 0..30, which is 31 stages, and the next length is n_31 = 2^3 = 8.
Each stage binary-searches one of 16 part indexes and uses at most 4
probes, because the candidate interval starts with 16 answers and each
probe cuts it. The final lower bound on an interval of length 8 uses at
most 4 probes: a probe on a span s replaces it by at most max(floor(s/2),
s-floor(s/2)-1), and that recurrence is 4 on span 8. Total probes per
giant are at most 31·4 + 4 = 128. Charged as 128·32 = 4096 instructions.

One LSD byte pass, at most 64 instructions per record, paid on all 32
bytes even when a unary byte skips the scatter. Histogram update of one
key byte is an address, a byte load, a counter load, an add and a store
(16 with loop overhead). The scatter copies 80 bytes as ten 8-byte loads
and ten stores, plus destination address arithmetic (48). Cap 64. The
256-bin prefix is 256 iterations of a few instructions per pass,
32·256·8 = 65536 instructions for the whole attack, inside the fixed
budget below. Across 32 key bytes and B records the paid body is
32·64·B = 1024·2^128 instructions.

## Work tally

Units below are coefficients of 2^128 instructions.

| Component | Coefficient of 2^128 |
| --- | ---: |
| Baby mixed additions, B·9728 | 4864 |
| Baby batch affine, 8·B·M | 3072 |
| Global LSD distribution, 32·64·B | 1024 |
| Giant mixed additions, G·9728 | 9728 |
| Giant batch affine, 8·G·M | 6144 |
| Giant sixteen-ary lookup, G·128·32 | 4096 |
| Giant batch inversions, (G/L)·I | 96 |
| Listed sum | 29024 |

The baby inversion is one I, not one per baby. Giant inversions are one per
block of L = 4096 points: G/L · 393216 = 96·2^128. Setup doubling (at most
130 doublings), the single baby inversion, histogram installation, decoding,
the final scalar reduction, [r]Q and the verifying scalar multiplication
are together below 2^30 instructions.

Unlisted attack work, charged in full and not folded into the kernels above:

- Baby record stores and address arithmetic, at most 16·B = 8·2^128.
- Baby loop control, at most 8·B = 4·2^128.
- Giant loop control and the stride pointer update, at most 16·G = 16·2^128.
- Final y and sign comparisons, at most 16·G = 16·2^128.
- Histogram clears, prefixes, setup doublings, the baby inversion and the
  verifying multiplication, together below 2^30 instructions, which is less
  than one coefficient of 2^128.

The sum of those coefficients is 8+4+16+16 = 44, and the fixed tail is
smaller than 2^128. Unlisted work is therefore strictly less than
176·2^128 instructions. Total
T < (29024 + 176)·2^128 = 29200·2^128.
W = T/1000 < 29.2·2^128.

Now 2^(49/10) = 2^4.9 > 29.2, because 292^10 < 100^10 · 2^49
(292^10 is 82 bits and 100^10 · 2^49 is 116 bits, and the integer
comparison holds). Therefore
29.2·2^128 < 2^(49/10)·2^128 = 2^(1329/10),
so W < 2^(1329/10) and log2(W) < 1329/10.

Preprocessing is the baby chain, the baby affine conversion and the LSD
distribution: coefficients 4864+3072+1024 = 8960, plus the baby portion of
the unlisted stores and loop control (12·2^128) and the fixed baby setup
(below 2^30). T_pre < 9000·2^128 and W_pre < 9·2^128 < 2^132.
Preprocessing is included in T. Advice is empty: the baby table is attack
work, not an embedded table. Nonuniform advice bytes are 0, and
log2(max(1, 0)) = 0.

## Memory

Peak is the batch-inversion phase of the baby table, before the LSD
second buffer exists.

- Jacobian babies: B·96 = 3·2^132 bytes.
- Prefix products of the batch inversion: B·32 = 2^132 bytes.
- Affine output records: B·80 = 2.5·2^132 bytes.

Sum 6.5·2^132. The LSD ping-pong buffers are 2·B·80 = 5·2^132 and are
live only after the Jacobian and prefix arrays are released. Giant blocks
hold L points and one 256-bin histogram, under 2^20 bytes with the code
and field scratch. Peak bytes < 7·2^132 < 8·2^132 = 2^135.
memory_log2_bytes = 135.

## Objections

The bound is an instruction count, not a timing model. A load is one
instruction whether or not it misses. The track scores instruction
equivalents divided by 1000, not memory latency, and it forbids treating a
P-256 addition as one reference operation.

The CIOS cap is not a measured cycle count. The 7-instruction inner step
uses 19 live registers (four multiplicand limbs, four prime limbs, six
accumulator limbs, and the running operand, product, carry and Montgomery
limb), which fits in the 31 general-purpose registers, so the bound does
not rely on an uncounted spill loop. Padding that step to 16 instructions
still fits in M = 768.

The odd-lift does not drop x = 0 and does not drop even scalars. Even x
is lifted by adding r, which flips parity because r is odd, and the range
check uses the profile upper bound on r rather than an approximation
2^256. Baby x-coordinates are unique for the parity reason above, so the
lookup has no collision case to charge.

The distribution is one global histogram per byte. It does not recurse into
buckets, does not plant a base per bucket, and does not keep a port queue.
A constant byte skips the scatter; the cost bound still pays for 32 full
passes. The sixteen-ary lookup does not allocate a heap: every pivot is an
index in the sorted record array. The 128-probe cap is the exact control
count on B = 2^127, not an average.

The experiment is exhaustive on the toy curve y^2 = x^3 + 2x + 2 over F_17
(order 19, the same odd-lift with the toy modulus 2^4) and checks global
LSD on 4096 synthetic records, including stability, unary skips and
cross-page scatters, plus agreement of the sixteen-ary lower bound with
the sorted order. It does not measure production instruction counts. The
production bound is the tally above.

## Submitted bound

| Quantity | Bound |
| --- | --- |
| time_log2 | 1329/10 |
| preprocessing_log2 | 132/1 |
| memory_log2_bytes | 135/1 |
| success_probability | 1 |
| nonuniform_advice_log2_bytes | 0 |

Heuristics: none. The group law, the tiling, the LSD permutation, the
sixteen-ary partition and the instruction caps are exact. The exponential
memory is required by the baby table of this algorithm.
