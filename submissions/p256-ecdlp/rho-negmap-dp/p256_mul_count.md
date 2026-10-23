# Instruction listing: P-256 field multiply, add and subtract (4 x 64-bit limbs)

This file fixes M and A used in `proof.md` by listing the actual instruction sequence and
counting it. Every line below is one scalar AArch64 instruction. The listing is executed and
checked line by line by `check_mul_mont.py` in this folder (a Python model of the same
instruction sequence, run on 20,000 random and 20,000 carry-rich inputs for the multiplication
and 20,000 random and edge-case inputs for each of the addition and the subtraction); the counts
in the tables are produced by that model from the same sequence.

**Flag semantics (AArch64).** `adds`, `adcs`, `subs` and `sbcs` write the carry flag C; `adc` and
`sbc` read C (for `sbc` the borrow is 1 - C) but do not write any flag. Hence every instruction
in a carry or borrow chain whose carry-out is read by the next instruction is the flag-setting
form; the plain `adc`/`sbc` appears only as the last instruction of a chain (after which the
next carry reader is preceded by a flag-writing `adds`/`subs`, or none follows). The model in
`check_mul_mont.py` implements exactly this rule, and it also asserts that demoting the
chain-internal `adcs`/`sbcs` to `adc`/`sbc` makes the multiplication wrong.

## 0. Parameters, registers, and what is counted

p = 2^256 - 2^224 + 2^192 + 2^96 - 1 (SEC 2 P-256), in 4 little-endian 64-bit limbs:

```
P0 = 0xFFFFFFFFFFFFFFFF
P1 = 0x00000000FFFFFFFF
P2 = 0x0000000000000000
P3 = 0xFFFFFFFF00000001
```

p ≡ -1 (mod 2^64), so (-p^-1 mod 2^64) = 1 and the Montgomery quotient digit of a round is
exactly the current low limb. The same fact makes m*p a sum of shifts of m, so the four
reduction rounds need no multiplication at all: 16 MUL + 16 UMULH for the 4 x 4 product, and
the reduction is adds, adcs, subs and sbcs.

Notation: `A0..A3`, `B0..B3` are the limbs of the two operands, `t0..t8` the 8-limb product and
its 9th (carry) limb, `P0..P3` the limbs of p, `lo`/`hi` the words of one 64x64->128 product,
`mlo`/`mhi` the 32-bit split of a reduction digit, `r0..r3`/`d0..d3` results and their
conditional difference, `mhi`/`mlo` reused in the final reduction.

**Counting convention.** Each operation is counted as the straight-line sequence that a walk
step inlines:

* operand loads from the walker and table arrays (`ldp`), and result stores (`stp`), are in the
  count;
* the limbs of p are held in four reserved registers for the whole walk, so no per-operation
  constant materialisation is charged here (the walk-loop setup line in `proof.md` covers it);
* there is no per-call prologue, no register spill and no `ret`, because these sequences are
  inlined in the walk loop; the loop's own state saving and control flow stay in the separate
  walk-loop line of `proof.md`;
* M is charged at S = M: a squaring is priced as a multiplication.

With that convention the exact count is **M = 139 instructions** for one Montgomery
multiplication, **20** for a modular addition and **19** for a modular subtraction.
`proof.md` charges M = 146 and A = 21, i.e. the counted values plus about 5 %
(139 * 1.05 = 145.95 -> 146, 20 * 1.05 = 21) as a reserve for scheduling, operand routing and
the difference between an analytic listing and a compiled or hand-scheduled kernel.

## 1. Multiplication, block A: setup (5 instructions)

```
    ldp    A0, A1, x0          ; a: four limbs of the walker coordinate
    ldp    A2, A3, x0, 16
    ldp    B0, B1, x1          ; b: four limbs of the second operand
    ldp    B2, B3, x1, 16
    mov    t2, xzr             ; third accumulator of column 1 starts at zero
```

Note: `t2` is the only accumulator word that has to be zeroed explicitly; `t0`, `t1` are written
by column 0, and every later column's high word is written by the previous column's `adc`.

## 2. Multiplication, block B: the 4 x 4 product, column (Comba) order (77 instructions)

The 16 limb products are summed by output column k, each column in a three-word accumulator
(t_k, t_k+1, t_k+2). One product term is five instructions: `mul`, `umulh`, and `adds`/`adcs`/`adc`
into the accumulator. The first term of a column writes the third accumulator word from the zero
register (`adc t_k+2, xzr, xzr`), which is why no column needs a separate clear. Column 0 holds
one product only, so it is a bare `mul`/`umulh` pair.

```
column 0: a0*b0
    mul    t0, A0, B0
    umulh  t1, A0, B0
column 1: a0*b1 + a1*b0                                    (10)
    mul    lo, A0, B1
    umulh  hi, A0, B1
    adds   t1, t1, lo
    adcs   t2, t2, hi
    adc    t3, xzr, xzr
    mul    lo, A1, B0
    umulh  hi, A1, B0
    adds   t1, t1, lo
    adcs   t2, t2, hi
    adc    t3, t3, xzr
column 2: a0*b2 + a1*b1 + a2*b0                            (15)
    mul    lo, A0, B2
    umulh  hi, A0, B2
    adds   t2, t2, lo
    adcs   t3, t3, hi
    adc    t4, xzr, xzr
    mul    lo, A1, B1
    umulh  hi, A1, B1
    adds   t2, t2, lo
    adcs   t3, t3, hi
    adc    t4, t4, xzr
    mul    lo, A2, B0
    umulh  hi, A2, B0
    adds   t2, t2, lo
    adcs   t3, t3, hi
    adc    t4, t4, xzr
column 3: a0*b3 + a1*b2 + a2*b1 + a3*b0                    (20)
    mul    lo, A0, B3
    umulh  hi, A0, B3
    adds   t3, t3, lo
    adcs   t4, t4, hi
    adc    t5, xzr, xzr
    mul    lo, A1, B2
    umulh  hi, A1, B2
    adds   t3, t3, lo
    adcs   t4, t4, hi
    adc    t5, t5, xzr
    mul    lo, A2, B1
    umulh  hi, A2, B1
    adds   t3, t3, lo
    adcs   t4, t4, hi
    adc    t5, t5, xzr
    mul    lo, A3, B0
    umulh  hi, A3, B0
    adds   t3, t3, lo
    adcs   t4, t4, hi
    adc    t5, t5, xzr
column 4: a1*b3 + a2*b2 + a3*b1                            (15)
    mul    lo, A1, B3
    umulh  hi, A1, B3
    adds   t4, t4, lo
    adcs   t5, t5, hi
    adc    t6, xzr, xzr
    mul    lo, A2, B2
    umulh  hi, A2, B2
    adds   t4, t4, lo
    adcs   t5, t5, hi
    adc    t6, t6, xzr
    mul    lo, A3, B1
    umulh  hi, A3, B1
    adds   t4, t4, lo
    adcs   t5, t5, hi
    adc    t6, t6, xzr
column 5: a2*b3 + a3*b2                                    (10)
    mul    lo, A2, B3
    umulh  hi, A2, B3
    adds   t5, t5, lo
    adcs   t6, t6, hi
    adc    t7, xzr, xzr
    mul    lo, A3, B2
    umulh  hi, A3, B2
    adds   t5, t5, lo
    adcs   t6, t6, hi
    adc    t7, t7, xzr
column 6: a3*b3                                            (5)
    mul    lo, A3, B3
    umulh  hi, A3, B3
    adds   t6, t6, lo
    adcs   t7, t7, hi
    adc    t8, xzr, xzr
```

The seven columns hold 2 + 10 + 15 + 20 + 15 + 10 + 5 = **77 instructions**: column 0 is a bare
`mul`/`umulh` pair, and the other 15 limb products are five instructions each (15 x 5 = 75). The
carry out of the last column is zero (the product of two 256-bit numbers fits in 8 limbs), which
is why `t8` is initialised for free.

## 3. Multiplication, block C: four Montgomery reduction rounds (46 instructions)

Round r uses the digit m = t_r (the current low limb, which is eliminated) and adds
m*p*2^(64r) to the accumulator. With p = 2^256 - 2^224 + 2^192 + 2^96 - 1:

```
m*p = m*2^256 - m*2^224 + m*2^192 + m*2^96 - m
```

and with mlo = (m & 0xffffffff) << 32 = m << 32, mhi = m >> 32:

* `-m` at limb r cancels the accumulator's low limb exactly (zero borrow, no instruction),
* `+m*2^96` puts `mlo` at limb r+1 and `mhi` at limb r+2,
* `+m*2^192` and `+m*2^256` put `m` at limbs r+3 and r+4,
* `-m*2^224` subtracts `mlo` at limb r+3 and `mhi` at limb r+4.

The positive terms are added first (`adds` then four chained `adcs`, then `adcs` into each
higher limb with the last one a plain `adc`), then the negative terms are subtracted (`subs`,
`sbcs`, then `sbcs` into each higher limb with the last one a plain `sbc`). The additions run before the subtractions so that every intermediate value is
non-negative; the total (T + sum of m*p*2^(64r)) is a non-negative number below 2^513, so no
borrow leaves the top.

```
round 0 (14 instructions)
    lsr    mhi, t0, 32
    lsl    mlo, t0, 32
    adds   t1, t1, mlo
    adcs   t2, t2, mhi
    adcs   t3, t3, t0
    adcs   t4, t4, t0
    adcs   t5, t5, xzr
    adcs   t6, t6, xzr
    adc    t7, t7, xzr
    subs   t3, t3, mlo
    sbcs   t4, t4, mhi
    sbcs   t5, t5, xzr
    sbcs   t6, t6, xzr
    sbc    t7, t7, xzr
round 1 (12 instructions)
    lsr    mhi, t1, 32
    lsl    mlo, t1, 32
    adds   t2, t2, mlo
    adcs   t3, t3, mhi
    adcs   t4, t4, t1
    adcs   t5, t5, t1
    adcs   t6, t6, xzr
    adc    t7, t7, xzr
    subs   t4, t4, mlo
    sbcs   t5, t5, mhi
    sbcs   t6, t6, xzr
    sbc    t7, t7, xzr
round 2 (10 instructions)
    lsr    mhi, t2, 32
    lsl    mlo, t2, 32
    adds   t3, t3, mlo
    adcs   t4, t4, mhi
    adcs   t5, t5, t2
    adcs   t6, t6, t2
    adc    t7, t7, xzr
    subs   t5, t5, mlo
    sbcs   t6, t6, mhi
    sbc    t7, t7, xzr
round 3 (10 instructions)
    lsr    mhi, t3, 32
    lsl    mlo, t3, 32
    adds   t4, t4, mlo
    adcs   t5, t5, mhi
    adcs   t6, t6, t3
    adcs   t7, t7, t3
    adc    t8, t8, xzr
    subs   t6, t6, mlo
    sbcs   t7, t7, mhi
    sbc    t8, t8, xzr
```

14 + 12 + 10 + 10 = **46 instructions**. Round r carries into limbs r+5..t7 (and, in the last
round, into `t8`) through a chain of `adcs` ended by one `adc` (and likewise `sbcs`/`sbc` for the
borrow), because a carry out of the reduced window can propagate into limbs that still
hold product bits.

## 4. Multiplication, block D: final conditional subtraction and store (11 instructions)

After the four rounds the quotient (t4..t7, t8) is below 2p < 2^257, so t8 is 0 or 1 and one
conditional subtraction of p finishes the reduction. The 5-limb subtract is
`subs`/`sbcs` on t4..t7 followed by `sbcs t8, t8, xzr`; the carry flag of that last `sbcs` is
set exactly when the value was at least p (t8 = 1 with no borrow, or t8 = 0 and t4..t7 >= p).

```
    subs   r0, t4, P0
    sbcs   r1, t5, P1
    sbcs   r2, t6, P2
    sbcs   r3, t7, P3
    sbcs   t8, t8, xzr
    csel   r0, r0, t4, cs
    csel   r1, r1, t5, cs
    csel   r2, r2, t6, cs
    csel   r3, r3, t7, cs
    stp    r0, r1, x2          ; result: four limbs, Montgomery form
    stp    r2, r3, x2, 16
```

**11 instructions.**

## 5. Multiplication: exact count

| Block | Instructions |
| --- | ---: |
| A setup: 4 `ldp`, 1 `mov` | 5 |
| B product: 2 + 15 x 5 | 77 |
| C reduction: 4 rounds | 46 |
| D final subtract, select, 2 `stp` | 11 |
| **Counted M** | **139** |
| Reserve 5 % (139 x 1.05 = 145.95) | 7 |
| **Charged M** | **146** |

By opcode: `mul` 16, `umulh` 16, `adds` 19, `adcs` 30, `adc` 19, `subs` 5, `sbcs` 11, `sbc` 4,
`lsr` 4, `lsl` 4, `csel` 4, `ldp` 4, `stp` 2, `mov` 1 = 139 (the model prints the same line).

## 6. Modular addition a + b mod p (20 instructions)

a + b < 2p, so the sum needs one conditional subtraction. `cy` holds the carry out of the
4-limb add; the subtraction `d = r - p` and `sbcs cy, cy, xzr` combine the 257-bit value's
borrow, and its carry flag is set exactly when the sum was at least p.

```
    ldp    A0, A1, x0
    ldp    A2, A3, x0, 16
    ldp    B0, B1, x1
    ldp    B2, B3, x1, 16
    adds   r0, A0, B0
    adcs   r1, A1, B1
    adcs   r2, A2, B2
    adcs   r3, A3, B3
    adc    cy, xzr, xzr
    subs   d0, r0, P0
    sbcs   d1, r1, P1
    sbcs   d2, r2, P2
    sbcs   d3, r3, P3
    sbcs   cy, cy, xzr
    csel   r0, d0, r0, cs
    csel   r1, d1, r1, cs
    csel   r2, d2, r2, cs
    csel   r3, d3, r3, cs
    stp    r0, r1, x2
    stp    r2, r3, x2, 16
```

**20 instructions** (4 loads, 5 add-with-carry including the carry word, 5 subtract-with-borrow,
4 select, 2 stores).

## 7. Modular subtraction a - b mod p (19 instructions)

`cy` is the borrow mask (all ones when a < b); p is masked with it and added back.

```
    ldp    A0, A1, x0
    ldp    A2, A3, x0, 16
    ldp    B0, B1, x1
    ldp    B2, B3, x1, 16
    subs   r0, A0, B0
    sbcs   r1, A1, B1
    sbcs   r2, A2, B2
    sbcs   r3, A3, B3
    sbc    cy, xzr, xzr
    and    d0, P0, cy
    and    d1, P1, cy
    and    d2, P2, cy
    and    d3, P3, cy
    adds   r0, r0, d0
    adcs   r1, r1, d1
    adcs   r2, r2, d2
    adcs   r3, r3, d3
    stp    r0, r1, x2
    stp    r2, r3, x2, 16
```

**19 instructions** (4 loads, 5 subtract-with-borrow including the mask, 4 masked limbs, 4
add-with-carry, 2 stores). A walk step uses both kinds; the ledger charges the more expensive
one (A = 21, i.e. the 20-instruction addition plus the same 5 % reserve) for all six add/sub
operations of a step, which over-charges three of them by two instructions.

## 8. Why this is below the earlier estimate and consistent with the reference prices

The earlier text charged M = 220 from "about 185" instructions; the exact count of the sequence
above is 139, so 220 was 58 % above the listed code, not 18 %. The count is also consistent
with the reference measurements in `tracks/docs/COST_MODEL.md`: the measured secp256k1 5x52
field multiplication is 179 instructions including its loads and stores, and a 4 x 64-limb
Montgomery multiplication has 16 partial products instead of 25 and needs no second reduction
pass, because p's low limb is -1 and the reduction is adds and shifts.

What the reserve and the excluded overheads do not cover, and what a reviewer should weigh:

* instruction count is not schedule quality: a Cortex-A78/A710 core issues about 4 scalar
  instructions per cycle, but the 16 UMULH and the long add-with-carry chains are dependency
  bound, so wall-clock time per multiplication is well above 139/4 cycles; the benchmark scores
  instructions, not cycles, so this does not change the claim, but it is why the 5 % reserve is a
  scheduling reserve and not a timing claim;
* the walk step's own work (table index, canonicalisation, `A`/`B` updates, address
  arithmetic, loop control) is not in this file: `proof.md` Section 4 itemises each such line
  as an explicit upper bound under its charged value;
* operands are P-256 field elements below p; the sequence is also correct for any 256-bit
  inputs, and the check script tests both.
