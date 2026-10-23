# P-256 field instruction listing (the counts behind proof.md section 5)

Every line in the three listings below is one scalar AArch64 instruction.
`experiment.py` Part C builds the same sequences (functions `listing_mul`,
`listing_add`, `listing_sub`; this file was printed from them), executes them
on a 64-bit register and memory model and checks every result. The counts are
the lengths of the executed straight-line sequences: there are no branches or
loops inside a kernel.

**Flag rule (AArch64).** `adds`, `adcs`, `subs` and `sbcs` write the carry flag
C. `adc` and `sbc` read C (for `sbc` the borrow is 1 − C) and write no flag.
Every instruction whose carry-out the next instruction reads is therefore the
flag-setting form; a plain `adc`/`sbc` ends a chain. The model implements this
rule exactly, and Part C checks that demoting the chain-internal `adcs`/`sbcs`
of the multiplication to `adc`/`sbc` makes its results wrong (the model notices
the flag rule).

**What is counted.**

- Operand loads (`ldp`) and result stores (`stp`) are in every kernel. Each
  kernel reads its operands from memory and writes its result to memory, so no
  value is kept in registers between kernels and kernels compose by address
  only. `x0`, `x1` and `x2` stand for base register plus immediate offset (the
  walk's element block, the table entry and a scratch block, each addressed by
  one base register); no per-kernel pointer arithmetic is needed.
- The limbs of p are held in reserved registers for the whole run: `P0`, `P1`,
  `P3`, and `P2` = 0 is `xzr`. No constant is rematerialised per kernel.
- The kernels are inlined in the step loop: no call, `ret`, prologue or spill.
  Registers: `A0..A3`, `B0..B3`, `t0..t8`, `lo`, `hi` (with `mlo`/`mhi` reusing
  `lo`/`hi` and `r0..r3`, `d0..d3`, `cy` reusing dead operand registers), the
  three p limbs and three base registers make 25 of the 28 general registers
  outside x18, x29 and x30. The loop's own instructions are not in these counts;
  proof.md section 5 itemises them in the control row.
- A squaring is charged as a multiplication (S = M).

p = 2^256 − 2^224 + 2^192 + 2^96 − 1, limbs (little-endian)
P0 = 0xFFFFFFFFFFFFFFFF, P1 = 0x00000000FFFFFFFF, P2 = 0, P3 = 0xFFFFFFFF00000001.
Field elements are in Montgomery form (x·2^256 mod p), always fully reduced
(< p): every kernel's output is below p, and all inputs are kernel outputs or
reduced constants.

## 1. Montgomery multiplication: 139 instructions

**Product.** The 16 limb products are summed by output column k into a
three-word accumulator (t_k, t_k+1, t_k+2): `mul`, `umulh`, then `adds`/`adcs`
into the low two words and `adc` into the third. The first term of a column
writes the third word from `xzr`, so no column needs a separate clear. Column
0 is a bare `mul`/`umulh` pair: 2 + 15 × 5 = 77.

**Reduction.** p ≡ −1 (mod 2^64), so −p^-1 mod 2^64 = 1 and the Montgomery
digit of round r is the current low limb m = t_r. Then
m·p = m·2^256 − m·2^224 + m·2^192 + m·2^96 − m, with m·2^96 split as
mlo = m << 32 at limb r+1 and mhi = m >> 32 at limb r+2. The −m term cancels
limb r exactly (no instruction; limb r is never read again). Each round adds
the positive terms with one `adds`/`adcs` chain, then subtracts m·2^224 with
one `subs`/`sbcs` chain, both carried up to the top limb. Rounds 0 to 2 carry
up to t7 and round 3 also into t8: 14 + 12 + 10 + 10 = 46.

**No carry is lost.** For reduced inputs the product T ≤ (p − 1)^2 <
2^512 − 2^480. Before round r the accumulated value is below
T + 2^(64r)·p, and the positive terms of round r add less than 2^(64r + 321).
So the value stays below 2^512 through round 2 (t0..t7 suffice) and below
2^514 in round 3 (t8 holds the top bits). Each round's subtraction chain
leaves the pre-round value plus m·(p + 1)·2^(64r), which is non-negative (the
skipped −m stays in limb r, which is never read again), so no borrow leaves the
top limb. The quotient after round 3 is below
(T + 2^256·p)/2^256 < 2p, so one conditional subtraction (11 instructions,
including the two stores) gives the reduced result.

By opcode: `mul` 16, `umulh` 16, `adds` 19, `adcs` 30, `adc` 19, `subs` 5,
`sbcs` 11, `sbc` 4, `lsr` 4, `lsl` 4, `csel` 4, `ldp` 4, `stp` 2, `mov` 1.
Part C prints the same histogram.

```
; setup: operand loads and one zeroed accumulator word (5)
  1  ldp    A0, A1, [x0]
  2  ldp    A2, A3, [x0, #16]
  3  ldp    B0, B1, [x1]
  4  ldp    B2, B3, [x1, #16]
  5  mov    t2, xzr
; 4 x 4 product in column order; one product term = mul, umulh, adds, adcs, adc (77)
  6  mul    t0, A0, B0
  7  umulh  t1, A0, B0
  8  mul    lo, A0, B1
  9  umulh  hi, A0, B1
 10  adds   t1, t1, lo
 11  adcs   t2, t2, hi
 12  adc    t3, xzr, xzr
 13  mul    lo, A1, B0
 14  umulh  hi, A1, B0
 15  adds   t1, t1, lo
 16  adcs   t2, t2, hi
 17  adc    t3, t3, xzr
 18  mul    lo, A0, B2
 19  umulh  hi, A0, B2
 20  adds   t2, t2, lo
 21  adcs   t3, t3, hi
 22  adc    t4, xzr, xzr
 23  mul    lo, A1, B1
 24  umulh  hi, A1, B1
 25  adds   t2, t2, lo
 26  adcs   t3, t3, hi
 27  adc    t4, t4, xzr
 28  mul    lo, A2, B0
 29  umulh  hi, A2, B0
 30  adds   t2, t2, lo
 31  adcs   t3, t3, hi
 32  adc    t4, t4, xzr
 33  mul    lo, A0, B3
 34  umulh  hi, A0, B3
 35  adds   t3, t3, lo
 36  adcs   t4, t4, hi
 37  adc    t5, xzr, xzr
 38  mul    lo, A1, B2
 39  umulh  hi, A1, B2
 40  adds   t3, t3, lo
 41  adcs   t4, t4, hi
 42  adc    t5, t5, xzr
 43  mul    lo, A2, B1
 44  umulh  hi, A2, B1
 45  adds   t3, t3, lo
 46  adcs   t4, t4, hi
 47  adc    t5, t5, xzr
 48  mul    lo, A3, B0
 49  umulh  hi, A3, B0
 50  adds   t3, t3, lo
 51  adcs   t4, t4, hi
 52  adc    t5, t5, xzr
 53  mul    lo, A1, B3
 54  umulh  hi, A1, B3
 55  adds   t4, t4, lo
 56  adcs   t5, t5, hi
 57  adc    t6, xzr, xzr
 58  mul    lo, A2, B2
 59  umulh  hi, A2, B2
 60  adds   t4, t4, lo
 61  adcs   t5, t5, hi
 62  adc    t6, t6, xzr
 63  mul    lo, A3, B1
 64  umulh  hi, A3, B1
 65  adds   t4, t4, lo
 66  adcs   t5, t5, hi
 67  adc    t6, t6, xzr
 68  mul    lo, A2, B3
 69  umulh  hi, A2, B3
 70  adds   t5, t5, lo
 71  adcs   t6, t6, hi
 72  adc    t7, xzr, xzr
 73  mul    lo, A3, B2
 74  umulh  hi, A3, B2
 75  adds   t5, t5, lo
 76  adcs   t6, t6, hi
 77  adc    t7, t7, xzr
 78  mul    lo, A3, B3
 79  umulh  hi, A3, B3
 80  adds   t6, t6, lo
 81  adcs   t7, t7, hi
 82  adc    t8, xzr, xzr
; reduction round 0 (14)
 83  lsr    mhi, t0, #32
 84  lsl    mlo, t0, #32
 85  adds   t1, t1, mlo
 86  adcs   t2, t2, mhi
 87  adcs   t3, t3, t0
 88  adcs   t4, t4, t0
 89  adcs   t5, t5, xzr
 90  adcs   t6, t6, xzr
 91  adc    t7, t7, xzr
 92  subs   t3, t3, mlo
 93  sbcs   t4, t4, mhi
 94  sbcs   t5, t5, xzr
 95  sbcs   t6, t6, xzr
 96  sbc    t7, t7, xzr
; reduction round 1 (12)
 97  lsr    mhi, t1, #32
 98  lsl    mlo, t1, #32
 99  adds   t2, t2, mlo
100  adcs   t3, t3, mhi
101  adcs   t4, t4, t1
102  adcs   t5, t5, t1
103  adcs   t6, t6, xzr
104  adc    t7, t7, xzr
105  subs   t4, t4, mlo
106  sbcs   t5, t5, mhi
107  sbcs   t6, t6, xzr
108  sbc    t7, t7, xzr
; reduction round 2 (10)
109  lsr    mhi, t2, #32
110  lsl    mlo, t2, #32
111  adds   t3, t3, mlo
112  adcs   t4, t4, mhi
113  adcs   t5, t5, t2
114  adcs   t6, t6, t2
115  adc    t7, t7, xzr
116  subs   t5, t5, mlo
117  sbcs   t6, t6, mhi
118  sbc    t7, t7, xzr
; reduction round 3 (10)
119  lsr    mhi, t3, #32
120  lsl    mlo, t3, #32
121  adds   t4, t4, mlo
122  adcs   t5, t5, mhi
123  adcs   t6, t6, t3
124  adcs   t7, t7, t3
125  adc    t8, t8, xzr
126  subs   t6, t6, mlo
127  sbcs   t7, t7, mhi
128  sbc    t8, t8, xzr
; final conditional subtraction of p, select, store (11)
129  subs   r0, t4, P0
130  sbcs   r1, t5, P1
131  sbcs   r2, t6, P2
132  sbcs   r3, t7, P3
133  sbcs   t8, t8, xzr
134  csel   r0, r0, t4, cs
135  csel   r1, r1, t5, cs
136  csel   r2, r2, t6, cs
137  csel   r3, r3, t7, cs
138  stp    r0, r1, [x2]
139  stp    r2, r3, [x2, #16]
```

## 2. Modular addition: 20 instructions

a + b < 2p. The 4-limb sum and its carry word `cy` form a 257-bit value. One
5-limb subtraction of p sets C exactly when the sum is at least p, and `csel`
keeps the difference in that case.

```
  1  ldp    A0, A1, [x0]
  2  ldp    A2, A3, [x0, #16]
  3  ldp    B0, B1, [x1]
  4  ldp    B2, B3, [x1, #16]
  5  adds   r0, A0, B0
  6  adcs   r1, A1, B1
  7  adcs   r2, A2, B2
  8  adcs   r3, A3, B3
  9  adc    cy, xzr, xzr
 10  subs   d0, r0, P0
 11  sbcs   d1, r1, P1
 12  sbcs   d2, r2, P2
 13  sbcs   d3, r3, P3
 14  sbcs   cy, cy, xzr
 15  csel   r0, d0, r0, cs
 16  csel   r1, d1, r1, cs
 17  csel   r2, d2, r2, cs
 18  csel   r3, d3, r3, cs
 19  stp    r0, r1, [x2]
 20  stp    r2, r3, [x2, #16]
```

## 3. Modular subtraction: 19 instructions

`cy` is the borrow mask (all ones when a < b). p is masked with it and added
back.

```
  1  ldp    A0, A1, [x0]
  2  ldp    A2, A3, [x0, #16]
  3  ldp    B0, B1, [x1]
  4  ldp    B2, B3, [x1, #16]
  5  subs   r0, A0, B0
  6  sbcs   r1, A1, B1
  7  sbcs   r2, A2, B2
  8  sbcs   r3, A3, B3
  9  sbc    cy, xzr, xzr
 10  and    d0, P0, cy
 11  and    d1, P1, cy
 12  and    d2, P2, cy
 13  and    d3, P3, cy
 14  adds   r0, r0, d0
 15  adcs   r1, r1, d1
 16  adcs   r2, r2, d2
 17  adcs   r3, r3, d3
 18  stp    r0, r1, [x2]
 19  stp    r2, r3, [x2, #16]
```

proof.md charges A = 20 for every addition and every subtraction.

## 4. Inversion: 53,098 instructions, charged I = 2^17

z^-1 = z^(p−2). p − 2 has 256 bits, 128 of them set. The left-to-right chain
for this fixed public exponent is straight-line: 255 squarings and 127
multiplications, 382 kernels × 139 = 53,098 instructions. Part C runs this
chain through the listing and checks z·z^-1 = 1. proof.md charges
I = 2^17 = 131,072.

## 5. One batched affine step

Montgomery's trick shares one inversion over a batch of n elements. Element i
has denominator d_i = x(T_i) − x(P_i) for an addition, or 2y(P_i) for a
doubling. The batch then runs:

- forward: d_i (one SUB, or one ADD for a doubling) and the prefix products
  c_i = c_(i−1)·d_i (one M each, c_0 = d_0);
- one inversion of c_(n−1);
- backward: inv_i = u·c_(i−1) and u ← u·d_i (two M each; inv_0 = u).

Per addition: λ = (y_T − y_P)·inv_i, x3 = λ² − x_P − x_T,
y3 = λ(x_P − x3) − y_P. With the batch that is 3 + 3 = 6 M and 6 SUB:
6·139 + 6·20 = 954 instructions charged (948 executed).

Per doubling (a = −3): numerator 3(x² − 1) (one M, one SUB, two ADD),
λ = num·inv_i, x3 = λ² − 2x, y3 = λ(x − x3) − y. With the batch that is 7 M
and 8 additions or subtractions: 7·139 + 8·20 = 1,133 charged.

Part C runs a batch of production P-256 points (additions and doublings mixed)
through the listing. The affine results equal the reference arithmetic. The
executed total equals 139·(M calls) + 20·(ADD calls) + 19·(SUB calls), with
6 M per addition, 7 M per doubling, 382 M for the inversion and 3 fewer M
because c_0 and inv_0 need no product.

## 6. Consistency with the reference measurement

`tracks/docs/COST_MODEL.md` measures the secp256k1 5×52 field multiplication
at 179 instructions, call overhead included. A 4×64 Montgomery multiplication
has 16 limb products instead of 25, and its reduction has no multiplications
because p's low limb is −1. Inlined, it is 139. These are instruction counts,
which the benchmark scores; they make no claim about cycles.
