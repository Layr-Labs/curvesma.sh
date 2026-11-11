# Scalar AArch64 implementation and instruction counts

All work is counted in scalar AArch64 instructions on 64-bit general-purpose
registers. Each emitted instruction is one architectural instruction. LDP and STP
count as one instruction each, as do BL and RET. There is one processor, and
nothing uses SIMD or vector registers.

The arithmetic kernels are not estimated. `experiment.py` (mode `kernels`)
**generates** them as straight-line programs and **executes** them on a small
interpreter that updates the 64-bit registers and the carry and zero flags. It
counts every instruction executed and compares every output with big-integer
arithmetic. The counts below are those printed values.

## Conventions

- **Residues.** Each residue is four little-endian words in Montgomery form
  (x * 2^256 mod m), canonical in [0, m), for m = p or m = r. Equality and order
  comparisons act on these words; Montgomery form is a bijection.
- **Constants.** The constants pool (p, -p^-1 mod 2^64, r, -r^-1 mod 2^64,
  Montgomery zeta, d, e and zero) is addressed through x28. The walk sets x28 once
  per run.
- **Encodability.** The generator asserts encodable immediates:
  - ADD immediate: 0..4095.
  - LDP/STP: -512..504, multiples of 8.
  - LDR: 0..32760, multiples of 8.
  - SUB/ADD sp: 0 to 4095 x 4096, in multiples of 4096.
  - LSR and EOR shift amounts: 0..63.
  - MOVZ: 16 bits.
  - x18 (the platform register) is never used.
- **Registers.** Used: x0-x17 (caller-saved), plus x19-x30. Each
  function saves and restores the callee-saved registers it uses as scratch with
  STP (pre-index) and LDP (post-index), and those instructions are counted; `glue`
  updates the walk-state registers x21-x23 on purpose. Peak
  stack use is 4,208 bytes.

## Kernels (executed counts, including RET; a call site adds 3 address moves and BL)

| Kernel | Instructions | Content |
| --- | ---: | --- |
| fp_mul (CIOS Montgomery, mod p) | 155 | 7 loads (a, b first half, p, n0), one save pair; per word of b: 4 MUL + 4 UMULH with an ADDS/ADCS chain into a 5-word accumulator, q = t0 * n0, 4 MUL + 4 UMULH for q*p with its carry chain; second-half b load; final SUBS/SBCS and 4 CSEL; 2 STP; restore; RET |
| fp_add (mod p) | 21 | 4 LDP, ADDS/ADCS x4, 2 LDP of p, SUBS/SBCS x4, CSEL x4, 2 STP, RET |
| fp_sub (mod p) | 21 | 4 LDP, SUBS/SBCS x4, 2 LDP of p, CSEL x4 (p or zero by borrow), ADDS/ADCS x4, 2 STP, RET |
| fr_mul, fr_add (mod r) | 155, 21 | the same programs with r's constants |
| inline reduced add/sub (in f12_mul, selector) | 20 | 4 LDP, 4-word add or subtract, 2 LDP of p, 4-word correction and CSEL, 2 STP |
| inline unreduced add (Fp2 Karatsuba sums) | 10 | 4 LDP, ADDS/ADCS x4, 2 STP |

**Range argument for fp_mul.** For canonical inputs (a, b < p): p < 2^254. Before each word of b, the accumulator t
is below 2p < 2^255. Adding a * b_i < p 2^64 and q p < 2^64 p keeps it below 2^320,
so five words never overflow. Division by 2^64 then leaves t < 2p, the standard
CIOS bound. One conditional subtraction gives the canonical result. For inputs below 2p, used by the
Fp12 product, see below: the accumulator stays at most 3p - 1. The
interpreter's random and edge tests cover all these operations: 0, 1, 2, p-1, p-2,
2^255 mod p, 2^256 mod p, 2^512 mod p, plus destination-equals-source.

## Fp12 product: 13,843 executed instructions

`f12_mul(out, a, b)` is Karatsuba at each tower level:
- Fp2: 3M + 5A, with its two Karatsuba sums left unreduced.
- Fp6: 6 Fp2 products, xi = 9+u by additions only.
- Fp12: 3 Fp6 products.

The experiment counts the operations: **54 fp_mul calls, 244 inline reduced
additions or subtractions, and 36 inline unreduced additions**.
- A call costs 159: the 155-instruction kernel, plus x0, x1, x2 address formation
  and BL.
- An inline reduced addition costs 20: 4 LDP for the operands, ADDS/ADCS x4, 2 LDP
  of p, SUBS/SBCS x4, CSEL x4, 2 STP. A subtraction is the same shape, with 4 CSEL
  selecting p or zero after the borrow.
- An inline unreduced addition costs 10: 4 LDP, ADDS/ADCS x4, 2 STP.

Operands sit in 32-byte slots at offsets of at most 480 from a base register:
- inputs from x19/x20 and output at x21;
- scratch at sp, sp+512 and sp+1024 (x22-x24), with at most 43 live slots.

So every LDP/STP immediate is encodable. Prologue and epilogue save x19-x25 and
x30, and allocate 4 KiB of stack. The output never aliases an input.

**Static bounds.** The generator tracks a bound B (value < B p) for every slot and
asserts it before each operation:
- reduced additions and subtractions take canonical inputs (B = 1);
- an unreduced sum has B = 2 (2p < 2^256);
- fp_mul requires (B_a p)(B_b p) < p 2^256, which holds for B_a = B_b = 2 because
  4p < 2^256.

Unreduced values only ever feed fp_mul.

**fp_mul with inputs below 2p.** This case is used by the Fp12 product. Each CIOS step gives
t' = (t + a b_i + q p)/2^64 with b_i, q <= 2^64 - 1 and a <= 2p - 1.
- *Accumulator:* if t <= 3p - 1, then t + a b_i + q p <= (3p - 1) + (2p - 1)(2^64 - 1) + p(2^64 - 1) = (3p - 1) 2^64,
  so t' <= 3p - 1. The accumulator therefore stays at most 3p - 1, and five words never overflow, because
  (3p - 1) 2^64 < 2^320.
- *Final value:* (a b + Q p)/2^256 <= ((2p - 1)^2 + (2^256 - 1) p)/2^256 < 2p, because 4p < 2^256. So one
  conditional subtraction still gives the canonical result.
- *Checks:* the generator asserts (B_a + 1) p < 2^256 for every multiplication, and
  the experiment tests inputs up to 2p - 1 against big integers.

Arithmetic check: 54 x 159 + 244 x 20 + 36 x 10 = 13,826, plus 17 instructions of
prologue, epilogue and base setup, = 13,843. The experiment checks the product
against a Python tower implementation and against schoolbook multiplication modulo
w^12 - 18 w^6 + 82 in the pinned basis.

## Class selector: 1,067 executed instructions

`selector(record, y)` reads e_1 = (a, b) from slots 6 and 7 of y.
- **Keys:** a*d, b*e, B = a*d + b*e, a*zeta, B*zeta, then a(zeta-1), B(zeta-1) and the
  six negations (4 fp_mul calls at 159, 9 inline reduced operations at 20).
- **Minimum:** for each of the 11 later keys, a 4-word SUBS/SBCS difference,
  3 ORR and CSET for "less" and "equal", and CSEL updates of the 4 minimum words,
  the index and the tie flag (21 instructions each).
- **Output:** stores K (4 words), k* and the tie flag.

The experiment compares K, k* and the tie flag on random GT elements with the
minimum taken directly over the 12 explicitly computed Frobenius conjugates
(Montgomery representatives).

## Per-slot glue: 357 executed instructions, the same on every path

The rest of a slot is a generated program, `glue`. It is branch-free: the same 357
instructions run whether the probe is accepted, hits the window, ties, or reaches a
distinguished key. The outcome is committed with CSEL, and only the final trail-exit
branch depends on data; it is counted, and the exit work is charged in the report and
trail-start rows. The program does the following, in order:
- reloads the walk constants through x26 = x28 + 512;
- forms the table entry ((h + a) AND (R-1)) x 12 + k*, plus the TW and EX addresses
  (MADD);
- calls f12_mul, then the selector into a candidate record (address moves and BL counted
  here, bodies counted above);
- computes both splitmix64 hash words of the candidate key (EOR with shifted register,
  MUL);
- compares the first word with all 24 stored window words (LDR, CMP, CSET, ORR each);
- forms the accept flag (no window hit and no tie);
- calls fr_add twice into temporaries and commits alpha and beta with CSEL;
- commits the candidate record, the window-ring entry and index, the next index h and
  the state-buffer swap with CSEL;
- sets the probe index to a + 1 on rejection and 0 on acceptance;
- computes the 87-bit distinguished flag (only when accepted);
- adds 1 to the slot and trail-slot counters and the accept flag to the step counter,
  each 2 words, compared with its 2-word cap;
- tests the 8-probe limit and branches once on the combined exit flag.

The experiment executes it on two accepted probes and two planted window hits per run, and checks every effect against Python:
- *both paths:* the counters, the flags and the probe index;
- *on acceptance:* the next state y * TW entry, the committed record, the ring, h, the
  exponents mod r and the swap;
- *on a hit:* that state, record, ring, h and exponents are unchanged.

It also asserts that the instruction count is identical on both paths.

**Slots chain without a call.** The walk sets x30 to the address of `glue` once per
trail. `glue` saves and restores x30 and sp, so its final RET branches to its own first
instruction and starts the next slot. No call, loop or branch instruction runs between
slots; a trail ends only through the counted CBNZ. The experiment models BL writing x30,
checks every RET target against the call stack, and checks that x30 and sp are
unchanged after each executed slot.

**Per slot:** 13,843 + 1,067 + 357 = **15,267** instructions.

Table reads are ordinary LDPs inside f12_mul from the 75,497,472-byte TW array. Its
address is the pointer computed above. No access is priced separately, and none
needs more than 64-bit addresses.

## Non-slot work

| Item | Charge (instructions) | Basis |
| --- | ---: | --- |
| exponentiation by a < 2^254 exponent | 508 products | square-and-multiply; squarings charged as products |
| Frobenius map | 1 product | 5 Fp2 products by gamma^k (15 fp_mul, 25 fp_add/fp_sub) plus 6 negations: under 3,200 |
| fresh coin word | 16 | one 64-bit word of the algorithm's ordinary fresh independent coins: one RNDR (FEAT_RNG) with its status branch, store, loop |
| random scalar | 128 draws x (4 x 16 + 64) = 16,384 | 4 coin words per draw, mask to 254 bits, compare with r; acceptance r/2^254 > 0.756 |
| table entry h | 1,028 products + 2^13 + 2 x 16,384 | two exponentiations, one product, 11 Frobenius maps; 22 fr_mul plus copies for EX; two scalars |
| trail start | 1,017 products + 1,067 + 2 x 16,384 + 2^10 | two exponentiations, one product, selector, two scalars, window and counter reset |
| report | 2^15 | canonical exponents (2 fr_mul by stored lambda^k), then red-black tree lookup and insert on 4-word keys |
| recovery | 508 products + 512 fr_mul + 2^10 | inversion mod r by exponentiation; check g^x = z (12-coefficient compare) |
| input check | 508 products + 2^16 | z^r = 1; decoding and Montgomery conversion |
| derived constants | 2^22 | zeta, gamma^k, lambda^k, Montgomery constants (exponentiations in Fp2 and mod r) |
| output and fixed overhead | 2^16 + 2^32 | encoding x; parsing, allocation, code loading |

**Reports in detail.** A red-black tree holding at most 2^40 nodes has height at
most 2 log2(2^40 + 1) < 82.
- Each level costs about 17 instructions: 2 LDP of the key, a 4-word compare, child
  pointer load and branch.
- Rebalancing is at most 82 recolourings and 2 rotations.
- Node allocation is a bump pointer.

That is under 4,000 instructions; 2^15 is charged.

## Storage

| Item | Bytes |
| --- | ---: |
| store: 2^40 nodes x 128 (K 32, alpha_c 32, beta_c 32, two children 16, parent and colour 16) | 140,737,488,355,328 |
| TW table 12 x 2^14 x 384 | 75,497,472 |
| EX table 12 x 2^14 x 64 | 12,582,912 |
| stack, scratch, window, states, code, constants | < 1,048,576 |
| **total** | **140,737,577,484,288 = 2^47.000001** |

All of it lies in one address space below 2^48, so ordinary 64-bit pointers reach it.
