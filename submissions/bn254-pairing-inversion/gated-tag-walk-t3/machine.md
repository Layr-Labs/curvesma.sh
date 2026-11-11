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
  constants: 5, THR, the 51-bit digit mask, the exit cap and the walk
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
  keeps its state in x19, x20 and x23 to x25 on purpose (and saves and restores x21, x22 and x27-x30, which hold digits), and `rebase` swaps x21 and x22. Peak
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


## Step loop: **117** executed instructions per walk step, **175** more for a second tag

`tt_run` is the loop of proof.md step 3. x19 points at the 60 digits of Y0 (it sits 240 bytes into the
480-byte digit array, so digit offsets lie in -240..232; the second threshold THR2 is stored at x19 + 240
once per run) and x20 at the entry e (the 60 first words lie at e - 480..e - 8 and the 60 second words at
e..e + 472, so every LDP offset lies in -480..464). The digits of Y0 do not change inside a call, so the 20
digits of pairs 20..29 are loaded into registers once per call and stay there: x4-x7, x9-x17, x21, x22 and
x26-x30 (x18 is never used). One trip:

| Part | Instructions |
| --- | ---: |
| pairs 0..19: 20 x (LDP of two digits, LDP of two first words, two MADD into x8; the first product is a MUL) | 80 |
| pairs 20..29: 10 x (LDP of two first words, two MADD with the digits held in registers) | 30 |
| certificate: MUL x0 = 5 A mod 2^64, CMP with THR, B.CS to the second tag | 3 |
| selector and pointer: UMULH x1 = floor(5 A / 2^64), ADD x2 = e + 8 x1, LDR of the child, CBNZ back to the top (a null child is the window end) | 4 |
| **one trip** | **117** |

When the first tag is not certified, the trip runs the second tag before the selector:

| Part | Instructions |
| --- | ---: |
| pairs 0..19: 20 x (LDP of two digits, LDP of two second words, two UMULH, two ADD into x8) | 120 |
| pairs 20..29: 10 x (LDP of two second words, two UMULH with the digits held in registers, two ADD) | 50 |
| certificate: LDR of THR2, MUL x0 = 5 A2 mod 2^64, CMP, B.CS to the uncertified exit | 4 |
| B back to the selector | 1 |
| **second tag** | **175** |

- **Entry** (16 instructions): three STP (pre-index) saving x21/x22, x27/x28 and x29/x30, the pool
  address, LDP of 5 and THR, LDR of the start entry's multiset size (kept in x23), and 10 LDP of the 20
  register-held digits.
- **Exits** restore the three pairs (three LDP, post-index), set the code (MOVZ) and return (RET):
  second tag not certified (code 1, the second tag is left in x8; 5 instructions) or window end (code 2,
  when CBNZ falls through on a null child; two MOV first put the selector in x10 and e + 8 x10 in x11, so the
  entry is x11 - 8 x10; 7 instructions). Entry and the longest exit are 23 instructions per call.
- **Registers preserved.** x19, x21, x22, x27-x30 and sp are unchanged by a call; x28 still points at
  the constants pool. The experiment sets x21, x22, x27, x29 and x30 to random values before each call on a
  real window and checks them, x19, x28 and sp after it.
- **Bounded without a trip counter.** Each entry stores the size of its multiset; the child of an entry of
  size k for any selector is the multiset with one more multiplier, of size k + 1, and entries of size L =
  3,500 have no children (proof.md step 2). So every trip of a call moves to an entry one larger, and a call
  makes at most L + 1 = 3,501 trips. The trips of a call are counted at the exit as the size of the exit
  entry minus the start size in x23, plus 1 (the re-base adds them to the walk trip counter).
- **No multiplication in the field.** A trip uses 61 MUL/MADD and 1 UMULH on 64-bit words (60 more UMULH for a
  second tag) and no field product. MADD and ADD accumulate modulo 2^64, which is the tag arithmetic of
  proof.md section 2.

The experiment executes the loop on three real windows of 3, 4 and 5 steps (a target z, five multipliers
g^c z^d, the table entries along the walk's own path written with their real first and second words, labels
and T_e); the first window is chosen to contain a step whose first tag is not certified, so the second tag
runs inside a real walk. It checks, against Python:
- every first tag against the lemma's bound E1 and every second tag against E2 with the exact rational
  2^64 t/p, the selector against floor(5 t/p) and that no element is distinguished;
- the exit code, the selector and entry registers, the start size and the preserved registers;
- that every full trip executes the same count once the second tags (175 each) are taken out, and the entry
  and exit counts above.

It also executes crafted entries whose tags make 5 A mod 2^64 equal to THR - 1 and 0 (first tag
certified), THR, THR + 1 and THR2 - 1 (second tag certified), THR2 and 2^64 - 5 (second tag not certified,
code 1 with A2 in x8), and 2^64 - 5 with a second word that adds 1 to the tag (second tag 0, certified,
selector 0); the difference of the counts of the two certified kinds is the 175 of the second tag.

## Re-base: 14,085 executed instructions

`rebase` runs after every exit. With x12 = e it does, in order:
- Y0 <- Y0 * T_e by `f12_mul` into the spare buffer, then swaps the two buffer pointers;
- (alpha, beta) <- (alpha + alpha_e, beta + beta_e) mod r by two `fr_add` calls in place;
- the 60 digits of the new Y0, two coordinates at a time: per coordinate 2 LDP, then AND, EXTR and AND
  (lsb 51), EXTR and AND (38), EXTR and AND (25), LSR (12); then 5 STP for the ten digits of the two: 25
  instructions per two coordinates;
- the run exit counter: 3 words (ADDS, ADCS, ADC), compared with the 3-word cap X_cap (SUBS, SBCS,
  SBCS), flag stored;
- the walk trip counter: 2 words, increased by the trips of the call (the size stored in the exit entry
  minus the start size x23, plus 1), compared with 2^97, flag stored.

The experiment checks after each of the three windows: the new Y0 against the Python product and against
the independently multiplied walk element; g^alpha z^beta = Y0 for the new labels; every digit and THR2;
both counters (starting from values that carry into the second and third word and that sit just below each
cap) and both flags; the swapped buffer pointers; sp and x27 unchanged. The count is the same in all runs.

**One exit is charged 15,132 instructions**: the re-base (14,085), the loop entry and exit (23), and
2^10 for the dispatch on the exit code, the exact evaluation of an uncertified element, the load of the
root's child pointer and the flag tests. The exact evaluation reads the four stored words of t, forms 5 t on
5 words (4 MUL, 4 UMULH, 4 ADDS/ADCS), compares it with p, 2p, ..., 4p (4 x 5 words) to find the selector
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
| derived constants | 2^22 | Montgomery constants, the 144 structure constants of the row, 2^(51 j) mod p; THR2 stored next to the digits |
| output and fixed overhead | 2^16 + 2^32 | encoding; parsing, allocation, code loading |

**Table entry in detail.** The entries are built in order of multiset size; an entry is made from its
parent (one multiplier fewer) and one multiplier:
- T_e: one Fp12 product (13,843) and two `fr_add` for the labels;
- the row v(T_e): 12 coordinates, each a fixed integer combination of the 12 coordinates of T_e with
  coefficients from the tower constants (products by 9 and sign changes, done by additions): under 12 x 12 x 4
  additions of 21 instructions, 12,096;
- the 60 residues u = v_k 2^(51 j) mod p: 60 `fp_mul` calls by stored constants, 60 x 159 = 9,540;
- the first and second words, floor(2^128 u/p) split into two words: each u by 128 steps of shift, compare
  and subtract on 4 words (at most 16 instructions a step), 60 x 2,048 = 122,880;
- the five child pointers and the size: the index of a multiset in the table by the combinatorial number
  system (5 additions of stored binomial coefficients) and one MADD: under 5 x 64;
- stores of the 184 words of the entry: 92 STP.

That is under 160,000 instructions; 2^20 = 1,048,576 is charged for each of the 4,395,611,121,935,076 entries.

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
| table: 4,395,611,121,935,076 entries x 1,472 (60 first words 480, 60 second words 480, five child pointers, the multiset size and padding 64, alpha_e 32, beta_e 32, T_e 384) | 6,470,339,571,488,431,872 |
| store: 2^50 nodes x 480 (element 384, a 32, b 32, two children 16, parent and colour 16) | 540,431,955,284,459,520 |
| stack, buffers, digits, state, code, constants | < 1,048,576 |
| **total** | **7,010,771,526,773,939,968 = 2^62.604279** |

All of it lies in one address space below 2^63, so ordinary 64-bit pointers reach it; table pointers are
plain addresses stored in the entries. The store is sized for the report cap 2^50; the expected number of
reports is below 2^38.

## Output step

- **Point output** (`outstep.py`, executed in every mode
  on 1, r - 1 and a seeded x): Jacobian doubling for a = 0 (7 products, 14 sums) and mixed addition of P0
  (11 products, 14 sums), counted operation by operation; the affine result must equal an independent
  affine double-and-add and lie on y^2 = x^3 + 3, and its operation count times the `fp_mul` count must stay
  within the point-output bound of proof.md section 5. The affine reference also checks [r]P0 = infinity.
- **Pairing check.** Priced from the executed Fp12 product and `fp_mul` counts (proof.md section 5); it is
  not executed.
