# Machine model and instruction ledger (lazy-reduction kernels)

All work is charged in scalar AArch64 instruction equivalents R and divided by C = 20,000
(`tracks/docs/COST_MODEL.md`, BN254 section). No BN254 operation is charged as one unit. Each
routine priced here is a **generated straight-line program** of scalar AArch64 instructions on the
64-bit general registers (no SIMD, no vector registers). `lazy.py` part L generates it, executes it
instruction by instruction on a simulator of the subset used (LDP/STP/LDR, MUL/UMULH,
ADDS/ADCS/ADC, SUBS/SBCS/SBC, EXTR/LSL/LSR, CSEL/CSETM/AND, UDIV, ADD-immediate, RET), compares
every output with big-integer arithmetic on the protected group, and reports the executed count.
Straight-line programs execute every instruction exactly once, so the executed count equals the
listing length on every input; the runs confirm this (`executed_instructions`). The listing of the
dense product is printed with input `{"seed":1,"listings":true}`.

The design follows the lazy-reduction method for BN pairing towers of Aranha, Karabina, Longa,
Gebotys and Lopez (EUROCRYPT 2011, eprint 2010/526): products of field elements are accumulated
unreduced in double width and Montgomery-reduced once per output coefficient.

## 1. Representation and register use

* An Fp element is stored as four little-endian 64-bit words in Montgomery form x~ = x 2^256 mod p,
  **canonical** (0 <= x~ < p) in every walk state and every routine output. p < 2^253.6, so 5p < 2^256.
* Double-width intermediates are 9-word two's-complement integers (80-byte slots, the tenth word
  unused); raw products are 8 nonnegative words.
* Tower: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 9 + u, the tower the profile basis
  maps to (`experiment.py` part A checks the map). For the product the coefficients are grouped
  as Fp4 = Fp2[s]/(s^2 - xi), s = w^3, and Fp12 = Fp4[t]/(t^3 - s), t = w:
  y = Y0 + Y1 t + Y2 t^2 with Y_j = c_j + c_{j+3} s. This is a regrouping of the same twelve stored
  coefficients, not a change of basis.
* Base registers: x0 output, x1 input y, x2/x3 the two 1 KiB windows of the table entry, x4 the
  constant pool, x5..x8 four 1 KiB scratch windows on the stack. Every LDP/STP offset is in
  [-512, 504] and a multiple of 8 (asserted by the generator). Computation uses x9..x17 and
  x19..x30 (21 registers): a 9-word accumulator, 9 operand temporaries and 3 more.
* Each routine saves and restores x19..x30 (6 STP, 6 LDP), allocates and frees its frame and sets its
  base registers (ADD-immediate), and ends with RET; all of this is inside the counted length.

## 2. Primitives (executed counts, measured on the generator)

| Primitive | Instructions | Content |
| --- | ---: | --- |
| 256x256 -> 512-bit product into the accumulator | 69 | 2 LDP (b), 4 LDR (a_i), 16 MUL, 16 UMULH, 31 ADDS/ADCS/ADC (operand scanning) |
| store 8-word product / 9-word value | 4 / 5 | STP pairs (the tenth word of a 9-word slot is xzr) |
| accumulator +- product (8 words) | 13 | 4 LDP + 9-word ADDS/ADCS/ADC (or SUBS/SBCS/SBC) chain |
| accumulator +- stored 9-word value | 14 | 5 LDP + 9-word chain |
| accumulator + 2^k * value, per extra signed digit | +18 | 9 shifted words (LSL, EXTR) interleaved with the chain (EXTR and LSL leave the carry flag intact) |
| accumulator + 9 * value / + 2 * value | 32 / 23 | 9 = 8 + 1, 2 = 2^1 |
| accumulator +- 2 c 2^256 (single-width term, squaring only) | 12 | 2 LDP, 5 shifted words into words 4..8 |
| add the routine constant K (a multiple of p) | 14 | 5 LDP + chain |
| unreduced single-width add / reduced add | 10 / 20 | 4 LDP, 4-word add, 2 STP (+ 2 LDP p, 4 SUBS/SBCS, 4 CSEL) |
| a - b + k p (k = 1, 2), unreduced | 16 | |
| REDC + canonicalisation + store | 107 | 3 LDP constants; 4 word rounds of Montgomery reduction on 9 words (78); q = UDIV(t >> 200, floor(p / 2^200) + 1), t - q p (4 MUL, 3 UMULH, 3 adds, 4 subtracts); one conditional subtraction of p (8); 2 STP |

**Static bounds (asserted while generating, worst case, not sampled).** Every stored value carries an
exact integer interval [mn, mx]. Operands of a product are below 2^256 (y-side sums at most 4p, the
squaring's s0 + s1 and s0 - s1 + 2p at most 4p); products are below 2^512; every 9-word value lies in
(-2^575, 2^575), so its 576-bit two's-complement image is exact and every accumulator chain computes
the true value modulo 2^576. Before reduction the routine adds K, the least multiple of p with
K >= -mn over all its outputs, so the reduced input T is in [0, 2^575). Montgomery reduction of a
9-word T in [0, Tmax] gives t < Tmax/2^256 + p, asserted below 2^264, so t >> 200 fits 64 bits. With
d = floor(p/2^200) + 1 and q = floor(floor(t/2^200)/d): q p <= t, and t - q p < t 2^200/p + p/d + p,
asserted below 2p; one conditional subtraction gives the canonical result. REDC is the identity on
congruence classes mod p of T 2^-256, and K is a multiple of p, so the output is
(sum of the products' true values) 2^-256 mod p = the Montgomery form of the field result.

## 3. The three routines (part L: executed counts, checks)

**Dense Fp12 product y * M, M a table entry (8,785 instructions).** Cubic Karatsuba over Fp4 (6 Fp4
products: v0, v1, v2, (Y1+Y2)(M1+M2), (Y0+Y1)(M0+M1), (Y0+Y2)(M0+M2)), each a Karatsuba Fp4 product of
3 Fp2 products, each a Karatsuba Fp2 product of 3 Fp products: 54 products. Fp2 results
(re = P0 - P1, im = P2 - P0 - P1) and Fp4 results (u0 + xi u1, u2 - u0 - u1) are kept unreduced as
9-word values; each of the 12 output coefficients is one accumulator over the Fp4-level values
(Z0 = v0 + s(t12 - v1 - v2), Z1 = t01 - v0 - v1 + s v2, Z2 = t02 - v0 + v1 - v2, s(A0 + A1 s) =
xi A1 + A0 s), plus K, then one REDC: **12 reductions instead of 54**. The y-side sums at the Fp4
level are reduced (12 reduced adds), deeper sums unreduced (< 4p). The table side is fixed: for each
of the 12,288 table entries setup stores the 54 canonical table-side operands (the 12 coefficients
and their Karatsuba sums), 1,728 bytes per entry, 21.2 MB in all, computed once with additions mod p
(charged in setup, below 2^20 units). Opcode histogram: 1,808 ADCS, 1,330 LDP, 1,104 MUL,
1,092 UMULH, 920 SBCS, 590 ADC, 558 STP, 542 ADDS, 216 LDR, 200 SUBS, 172 EXTR, 116 SBC, 96 CSEL,
20 LSL, 12 UDIV, 8 ADD, 1 RET. Peak scratch: 34 slots.
Checks: 30 cases (random, edge-valued coefficients {0, 1, 2, p-1, p-2, ...}, all-limbs-large), all
outputs canonical and equal to the schoolbook product mod (w^6 - xi).

**Granger-Scott cyclotomic squaring (4,323 instructions).** For each Fp4 pair (a, b): Fp2 squares
from two products each (x^2 = ((x0+x1)(x0-x1+kp), 2 x0 x1)), for a, b and a+b: 18 products in all.
X0 = a^2 + xi b^2 and X1 = (a+b)^2 - a^2 - b^2 (stored as half) are 9-word values; each output
coefficient is 3X -+ 2c (c the input coefficient, added as c 2^256 so that REDC returns it in Montgomery
form), plus K, one REDC: 12 reductions. Checks: 4 GT elements g^k of the protected group per run, all
outputs canonical and equal to y^2.

**Twelve Frobenius keys and their argmin (898 instructions).** B = a Re(d1) + b Im(d1) (two products,
one reduction), a zeta and B zeta (one product and reduction each), the six values
x, x zeta, x zeta - x, -x, -x zeta, x - x zeta mod p for x = a and x = B, and a running minimum over the
twelve stored keys with a full four-word comparison per candidate (SUBS/SBCS into xzr, 4 CSEL for
the minimum, MOV + CSEL for the index; the first minimum wins on equality). Keys are distinct except
with probability at most 66/p per element (proof.md section 2), a failure event already charged.
Checks: the returned index equals the argmin of the twelve stored keys computed from the **twelve
explicitly computed conjugates** y^(p^k), on GT elements and on elements with edge-valued w^1
coefficients.

## 4. Per-step bookkeeping (COMMON, charged 300)

| Work | Instructions |
| --- | ---: |
| load limbs 0 and 1 of stored T = Re(c_0); branch bits, table index h, DP test across two limbs | 12 |
| table address ((12h + t) * 1728 + base) and twist index t = (12 - kappa) mod 12 | 7 |
| exponent update mod r: two modular additions (base package fp_add listing, 26 each with call site; `experiment.py` part K executes it for r) | 52 |
| 2-cycle test: four limbs of T against T two steps back | 10 |
| checkpoint test: four limbs of T against the checkpoint T | 10 |
| history: store T and rotate the three state-buffer pointers | 12 |
| step counter of the walk, checkpoint period test, amortised checkpoint store | 6 |
| stuck-walk counter (add, cmp, branch) | 3 |
| **ordinary-step counter**: 192-bit (3-word ADDS/ADCS/ADC) increment, 3-word compare with S_cap + E (SUBS/SBCS/SBCS) and branch | 7 |
| **all-step counter** (also incremented by cycle-escape steps): 192-bit increment, compare with 2(S_cap + E), branch | 7 |
| branch dispatch | 3 |
| walk loop: next walk-state pointer, index increment, wrap test, branch | 8 |
| **itemised** | **137** |
| **charged** | **300** |

The caps are at most 2^132.3 < 2^192, so the 192-bit counters never wrap. Call sites: 4 instructions
per routine call (BL and up to three pointer set-ups); a multiply step calls the key routine and the
product (8), a squaring step one routine (4).

| Step | Routines | Call sites | Bookkeeping | Instructions | Units |
| --- | ---: | ---: | ---: | ---: | ---: |
| multiply | 8,785 + 898 | 8 | 300 | **9,991** | **0.49955** |
| squaring | 4,323 | 4 | 300 | **4,627** | **0.23135** |

With squaring probability s = 1/2 the average step costs **c = 0.36545 units** (base package 0.6289).

## 5. Reports, setup and output

A distinguished-point report computes the keys if the step did not (898), multiplies alpha and beta by
the stored constant p^kappa 2^256 mod r (2 x 162 with the base fp_mul listing), extracts the 128-bit
fingerprint, hashes it into the chained multimap and writes an 80-byte record: under 2,000
instructions, **charged 10,000** (0.5 unit). Setup (pairing-free: g is a profile constant): conversion
of z, z^r = 1, the 1,024 multipliers g^a z^b with their 11 Frobenius images and their 54 precomputed
operands each, 2^20 walk starts, each exponentiation at most 254 squarings and 254 products (under 400
units even at the base prices), all under 2^31 units. Verification of a candidate (one inverse mod r
and g^x = z) is under 600 units.
