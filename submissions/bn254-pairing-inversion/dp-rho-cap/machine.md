# Scalar AArch64 implementation and instruction counts

All counts are scalar AArch64 instruction counts on 64-bit general registers.
They are not timings and are not derived from the normalization constant C.
The priced walk step is a concrete straight-line instruction sequence. The
generator `a64.py` emits it, and `experiment.py` executes the emitted sequence
on production BN254 values. Its length is the executed count, because the step
has no internal loops. There is one processor.

## Instruction set used

Each emitted tuple in `a64.py` is one A64 instruction:

- MUL and UMULH for the low and high words of a 64x64 product;
- ADDS, ADCS and ADC (carry chains), SUBS and SBCS (borrow chains);
- CSEL on CS, CSETM on CC, and register AND, for selection and masking;
- AND with the logical immediate 0x3ff, ORR with a shifted register, and ADD
  and SUB with a shifted register (LSL #9 and LSL #6);
- LDP and STP with a signed 7-bit immediate scaled by 8 (offsets -512..504),
  LDR and STR with an unsigned scaled offset, ADD with a 12-bit immediate;
- CBZ, B.CS and B.

`Prog.e` asserts at generation time that every immediate is encodable and that
no LDP has two equal destinations. The experiment checks that these guards
reject seven malformed instructions: an out-of-range LDP offset (512), an LDP
with equal destinations, a misaligned LDR offset (4), and AND masks 0x5,
2^64-1, 2^65-1 and -1, none a contiguous low mask of 1..63 bits. The emulator
models carries as A64 does
(C = NOT borrow after SUBS/SBCS). No vector unit, 128-bit unit-cost operation,
division instruction or table lookup outside ordinary loads is used.

Registers: x0-x17 are temporaries. x19-x24 hold the bases of the current
element, the selected table record, the next element, two scratch windows and
the control block. x25-x28 hold p and x29 holds -p^-1 mod 2^64. x18 and the
stack are not used by the step.

## Field elements

An Fp element is four little-endian words in Montgomery form v*2^256 mod p,
always canonical (< p). Fp12 uses the tower
Fp2 = Fp[u]/(u^2+1), Fp6 = Fp2[v]/(v^3-(9+u)), Fp12 = Fp6[w]/(w^2-v), where w is
the pinned basis variable, so v = w^2 and u = w^6 - 9. Coefficient i < 6 of
the pinned basis pairs with coefficient i+6 as A_i = (c_i + 9 c_{i+6}) +
c_{i+6} u, and A_{k+2m} sits at Fp2 slot m of Fp6 part k. The experiment
checks this map and its inverse, and checks that the tower product equals the
pinned-basis product on 8 inputs, including g*g and the all-(p-1) input.

## Kernel counts (emitted and executed)

| Sequence | Instructions | Content |
|---|---:|---|
| Fp Montgomery product (CIOS) | 149 | 4 LDP, 16 MUL + 16 UMULH product, 4 rounds of (1 + 4 MUL, 4 UMULH) reduction with carries, 4-word conditional subtraction, 2 STP |
| Fp addition | 18 | 4 LDP, 4-word add, 4-word trial subtraction, 4 CSEL, 2 STP (16 when doubling one operand) |
| Fp subtraction | 19 | 4 LDP, 4-word subtract, CSETM, 4 AND, 4-word add-back, 2 STP |
| Fp2 product (Karatsuba) | 540 | 3 products, 2 additions, 3 subtractions |
| Fp2 product by 9+u | 169 | 9a via three doublings and one addition, for both parts |
| Fp6 product | 4,130 | 6 Fp2 products, 2 products by 9+u, Fp2 additions and subtractions |
| Fp12 product | 13,111 | 3 Fp6 products, 2 Fp6 additions, 2 Fp6 subtractions, 1 product by 9+u, 3 Fp2 additions |

The product needs no reduction beyond the CIOS rounds. Each partial sum fits in
five words, and the CIOS result is below 2p < 2^255, so one conditional
subtraction makes it canonical. Additions cannot overflow four words because
2p < 2^255; the same holds modulo r for the label additions. The experiment
checks the product, addition and subtraction on all 25 pairs of edge values
{0, 1, 2, m-2, m-1} and on 400 random pairs (addition also modulo r).

## The walk step

One step, as emitted by `step_program(cur, nxt)`:

| Part | Instructions |
|---|---:|
| partition: LDR word 0 of A0.im of the current element, AND 0x3ff | 2 |
| record address: LDR table base, ADD j LSL 9, SUB j LSL 6 (448 j) | 3 |
| Fp12 product into the next element | 13,111 |
| load r (2 LDP), a += c_j and b += d_j modulo r | 2 + 36 |
| distinguished-point test: LDP, ORR (word1 LSL 34), CBZ | 3 |
| 2-word step counter: LDP, ADDS, ADC, STP; compare with the cap: LDP, SUBS, SBCS, B.CS; loop branch B | 9 |
| **total emitted / charged** | **13,166 / 13,312** |

The loop body holds two copies of the step with the roles of the two element
buffers exchanged, so no copying is needed. The charge of 13,312 = 13 x 1,024
leaves 146 instructions of slack per step. Scratch for the product uses at
most 42 of the 64 slots of 32 bytes in the two 1 KiB windows.

The distinguished-point path (`dp_program`) runs only when the test succeeds:
it forms the address of record n of the distinguished-point table, copies the
384-byte element and both 32-byte labels (28 LDP and 28 STP), increments n and
compares it with D_max. It emits 66 instructions and is charged 128.

## Bounded work outside the walk

Fp12 exponentiation by a scalar below 2^254 uses at most 508 Fp12 products
(below 2^14 instructions each) and 2^16 instructions of control: below 2^23.

| Item | Bound (instructions) |
|---|---:|
| decode 12 coefficients, canonical checks, Montgomery conversion, z^r = 1 | 2^24 |
| 2,050 scalars by 254-bit rejection sampling, at most 1,024 draws each, 64 instructions per random bit | 2^36 |
| 1,025 elements g^c z^d (two exponentiations and one product each), records | 2^35 |
| sort of at most 2^36 eight-byte record indices by full element (merge sort, 36 passes, 512 instructions per comparison and move), then one adjacent scan | 2^51 |
| at most 16 recovery attempts: inversion modulo r by exponentiation (508 products modulo r), one product, and the check g^x = z | 2^29 |
| [x]P0 by 508 Jacobian doublings and additions (each under 2^12), one Fp inversion, the curve check and [r]P = infinity | 2^23 |
| attack-side reduced pairing check e(P, Q0) = z: optimal ate Miller loop over 6x+2 (65 bits, under 400 line evaluations and Fp12 products), G2 arithmetic, and the final exponentiation by binary powering of (p^12-1)/r (under 2^12 bits, so under 2^13 Fp12 products) | 2^28 |
| allocation, input and output, all remaining fixed control | 2^32 |
| **fixed allowance F** | **2^52** |

Random bits are the algorithm's ordinary fresh independent coins, read 64 bits per
load; no generator expands a seed, so E1, E5 and H1 use ideal uniform coins as in
the protected probability space. They are not an uncharged group or advice oracle.
The allowance of 64 instructions per bit (4,096 per 64-bit word) covers the load,
assembling each 254-bit candidate and the comparison with r.

The listed items sum to below 2^51.1, so F = 2^52 leaves slack.

## Memory

| Store | Bytes |
|---|---:|
| distinguished-point table, D_max = 2^36 records of 448 bytes | 2^44.8074 |
| two index arrays for the merge sort, 2^36 eight-byte entries each | 2^40 |
| adder table, 1,024 records of 448 bytes | 458,752 |
| element buffers, scratch windows, control block, code, stacks | below 2^20 |
| **peak total** | **below 2^45** |

All of it is addressed with ordinary 64-bit pointers; no segmented, banked or
multi-word addressing is used.
