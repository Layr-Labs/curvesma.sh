# Scalar AArch64 implementation and instruction counts

All work is counted in scalar AArch64 instructions on 64-bit general-purpose
registers. Each emitted instruction is one architectural instruction. LDP and STP
count as one instruction each, as do BL and RET. There is one processor, and
nothing uses SIMD or vector registers.

The priced programs are not estimated. `experiment.py` (mode `kernels`, with `kern.py`)
**generates** them and **executes** them on a small interpreter that updates the 64-bit
registers and the carry and zero flags and follows real branches. It counts every instruction
executed and compares every effect with big-integer arithmetic. The counts below are those
printed values. A memory word that no program or test wrote reads as a scramble of its address,
never as zero, so no listing can rely on zeroed memory.

## Conventions

- **Residues.** Each residue is four little-endian words in Montgomery form
  (x * 2^256 mod m), canonical in [0, m), for m = p or m = r. Equality and order
  comparisons act on these words; Montgomery form is a bijection.
- **Constants.** The constants pool (p, -p^-1 mod 2^64, r, -r^-1 mod 2^64, and the walk
  constants: 6, THR, the trip budget 1,025, the 43-bit digit mask, the exit cap and the walk
  trip cap) is addressed through x28. The walk sets x28 once per run.
- **Encodability.** The generator asserts encodable immediates:
  - ADD immediate: 0..4095.
  - LDP/STP: -512..504, multiples of 8.
  - LDR: 0..32760, multiples of 8.
  - SUB/ADD sp: 0 to 4095 x 4096, in multiples of 4096.
  - LSR and LSL shift amounts, EXTR lsb, ADD (shifted register) amount: 0..63.
  - SUBS immediate: 0..4095. MOVZ: 16 bits.
  - Branch and call targets are labels.
  - x18 (the platform register) is never used.
- **Registers.** Used: x0-x17 (caller-saved), plus x19-x30. Each
  function saves and restores the callee-saved registers it uses as scratch with
  STP (pre-index) and LDP (post-index), and those instructions are counted; the step loop
  keeps its state in x19, x20 and x23 to x25 on purpose, and `rebase` swaps x21 and x22. Peak
  stack use is 4,208 bytes.

## Kernels (executed counts, including RET; a call site adds 3 address moves and BL)

| Kernel | Instructions | Content |
| --- | ---: | --- |
| fp_mul (CIOS Montgomery, mod p) | 155 | 7 loads (a, b first half, p, n0), one save pair; per word of b: 4 MUL + 4 UMULH with an ADDS/ADCS chain into a 5-word accumulator, q = t0 * n0, 4 MUL + 4 UMULH for q*p with its carry chain; second-half b load; final SUBS/SBCS and 4 CSEL; 2 STP; restore; RET |
| fp_add (mod p) | 21 | 4 LDP, ADDS/ADCS x4, 2 LDP of p, SUBS/SBCS x4, CSEL x4, 2 STP, RET |
| fp_sub (mod p) | 21 | 4 LDP, SUBS/SBCS x4, 2 LDP of p, CSEL x4 (p or zero by borrow), ADDS/ADCS x4, 2 STP, RET |
| fr_mul, fr_add (mod r) | 155, 21 | the same programs with r's constants |
| inline reduced add/sub (in f12_mul) | 20 | 4 LDP, 4-word add or subtract, 2 LDP of p, 4-word correction and CSEL, 2 STP |
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


## Step loop: **153** executed instructions per walk step

`tt_run` is the loop of proof.md step 3. x19 points at the 72 digits of Y0 and x20 at the entry e
(both pointers sit 288 bytes into their 576-byte arrays, so every LDP offset lies in -288..280).
One trip:

| Part | Instructions |
| --- | ---: |
| 36 x (LDP of two digits, LDP of two table words, two MADD into x8; the first product is a MUL) | 144 |
| certificate: MUL x9 = 6 A mod 2^64, CMP with THR, B.CS to the uncertified exit | 3 |
| selector and pointer: UMULH x10 = floor(6 A / 2^64), ADD x11 = e + 8 x10, LDR of the child, CBZ to the window-end exit | 4 |
| trip budget: SUBS x23, B.NE back to the top | 2 |
| **one trip** | **153** |

- **Entry** (3 instructions): the pool address, LDP of 6 and THR, LDR of the trip budget 1,025.
- **Exits** (MOVZ of a code and RET, 2 instructions): not certified (code 1, the tag is left in x8),
  window end (code 2, the selector in x10 and e + 8 x10 in x11, so the entry is x11 - 8 x10), budget
  used (code 3).
- **Bounded.** x23 starts at 1,025 and one trip takes 1 from it, so a call makes at most 1,025 trips,
  whatever the table holds.
- **No multiplication in the field.** The trip uses 73 MUL/MADD and 1 UMULH on 64-bit words and no
  field product. MADD accumulates modulo 2^64, which is the tag of proof.md section 2.

The experiment executes the loop on three real windows of 3, 4 and 5 steps (a target z, six multipliers
g^c z^d, the table entries along the walk's own path written with their real words, labels and T_e).
It checks, against Python:
- every tag against the lemma with the exact rational 2^64 t/p, the selector against floor(6 t/p) and
  that no element is distinguished;
- the exit code, the selector and entry registers and the trips left;
- that every full trip executes the same count.

It also executes crafted entries whose tag makes 6 A mod 2^64 equal to THR - 1 and 0 (certified) and
THR + 1 and 2^64 - 6 (not certified), and a ring of two entries on which a budget of 5 trips ends the call.

## Re-base: 14,115 executed instructions

`rebase` runs after every exit. With x12 = e it does, in order:
- Y0 <- Y0 * T_e by `f12_mul` into the spare buffer, then swaps the two buffer pointers;
- (alpha, beta) <- (alpha + alpha_e, beta + beta_e) mod r by two `fr_add` calls in place;
- the 72 digits of the new Y0: per coordinate 2 LDP, then AND, EXTR and AND (lsb 43), EXTR and AND (22),
  LSR and AND (1), EXTR and AND (44), LSR (23), and 3 STP: 15 instructions;
- the run exit counter: 3 words (ADDS, ADCS, ADC), compared with the 3-word cap X_cap (SUBS, SBCS,
  SBCS), flag stored;
- the walk trip counter: 2 words, increased by the trips of the call (1,025 - x23 + 1), compared with
  2^97, flag stored.

The experiment checks after each of the three windows: the new Y0 against the Python product and against
the independently multiplied walk element; g^alpha z^beta = Y0 for the new labels; every digit; both
counters (starting from values that carry into the second and third word and that sit just below each
cap) and both flags; the swapped buffer pointers; sp and x27 unchanged. The count is the same in all runs.

**One exit is charged 15,144 instructions**: the re-base (14,115), the loop entry and exit (5), and
2^10 for the dispatch on the exit code, the exact evaluation of an uncertified element, the load of the
root's child pointer and the flag tests. The exact evaluation reads the four stored words of t, forms 6 t on
5 words (4 MUL, 4 UMULH, 4 ADDS/ADCS), compares it with p, 2p, ..., 5p (5 x 5 words) to find the selector
and with (s + 1) p - floor(p/2^90) for the distinguished test: under 150 instructions. Every exit is charged
the full amount, also those that need no product (an uncertified tag at the root).

## Other work

| Item | Charge (instructions) | Basis |
| --- | ---: | --- |
| run gate | 64 | one coin word (16), the compare with K (LDR, CMP, B.CS) and the return; charged on every run, open or closed |
| exponentiation by a < 2^254 exponent | 508 products | square-and-multiply; squarings charged as products |
| fresh coin word | 16 | one 64-bit word of the algorithm's ordinary fresh independent coins: one RNDR (FEAT_RNG) with its status branch, store, loop |
| random scalar | 128 draws x (4 x 16 + 64) = 16,384 | 4 coin words per draw, mask to 254 bits, compare with r; acceptance r/2^254 > 0.756; a draw that fails 128 times fails the run |
| multiplier | 1,017 products + 2 x 16,384 | two exponentiations, one product; two scalars |
| table entry | 2^20 | see below |
| walk start | 1,017 products + 2 x 16,384 + 2^10 | two exponentiations, one product, two scalars; digits, counter reset, compare with 1 |
| report | 2^18 | red-black tree lookup and insert on 48-word keys, see below |
| recovery | 508 products + 512 fr_mul + 2^10 | inversion mod r by exponentiation; check g^x = z (12-coefficient compare) |
| input check | 508 products + 2^16 | z^r = 1; decoding and Montgomery conversion |
| derived constants | 2^22 | Montgomery constants, the 144 structure constants of the row, 2^(43 j) mod p |
| output and fixed overhead | 2^16 + 2^32 | encoding x; parsing, allocation, code loading |

**Table entry in detail.** The entries are built in order of multiset size; an entry is made from its
parent (one multiplier fewer) and one multiplier:
- T_e: one Fp12 product (13,843) and two `fr_add` for the labels;
- the row v(T_e): 12 coordinates, each a fixed integer combination of the 12 coordinates of T_e with
  coefficients from the tower constants (products by 9 and sign changes, done by additions): under 12 x 12 x 4
  additions of 21 instructions, 12,096;
- the 72 residues v_k 2^(43 j) mod p: 72 `fp_mul` calls by stored constants, 72 x 159 = 11,448;
- the 72 words floor(2^64 u/p): each by 64 steps of shift, compare and subtract on 4 words (at most 16
  instructions a step), 72 x 1,024 = 73,728;
- the six child pointers: the index of a multiset in the table by the combinatorial number system
  (6 additions of stored binomial coefficients) and one MADD: under 6 x 64;
- stores of the 136 words of the entry: 68 STP.

That is under 112,000 instructions; 2^20 = 1,048,576 is charged for each of the 1,634,386,955,681,025 entries.

**Reports in detail.** A red-black tree holding at most 2^50 nodes has height at
most 2 log2(2^50 + 1) < 102.
- Each level costs at most 160 instructions: 48 LDP for the two 48-word keys, a 48-word compare
  (CMP and 47 CCMP), child pointer load and branch.
- Rebalancing is at most 102 recolourings and 2 rotations: under 2,000.
- Node allocation is a bump pointer; storing the node (the element, both labels, pointers) is 30 STP.

That is under 102 x 160 + 2,000 + 100 = 18,420 instructions; 2^18 = 262,144 is charged.

## Storage

| Item | Bytes |
| --- | ---: |
| table: 1,634,386,955,681,025 entries x 1,088 (72 words 576, six child pointers and padding 64, alpha_e 32, beta_e 32, T_e 384) | 1,778,213,007,780,955,200 |
| store: 2^50 nodes x 480 (element 384, a 32, b 32, two children 16, parent and colour 16) | 540,431,955,284,459,520 |
| stack, buffers, digits, state, code, constants | < 1,048,576 |
| **total** | **2,318,644,963,066,463,296 = 2^61.007988** |

All of it lies in one address space below 2^63, so ordinary 64-bit pointers reach it; table pointers are
plain addresses stored in the entries. The store is sized for the report cap 2^50; the expected number of
reports is below 2^38.
