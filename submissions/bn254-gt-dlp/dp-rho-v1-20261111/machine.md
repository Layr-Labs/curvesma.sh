# Priced scalar implementation and instruction bounds

All counts are **executed scalar AArch64 instructions** (64-bit general
registers; arithmetic, logic, loads/stores, address formation, calls and
returns all count one each).  They are analytical bounds read off supplied
straight-line code, not timings.  One processor.

## Evidence chain for the per-step count

* `arithmetic.c` is the complete source of the step: `fp_add`, `fr_add`,
  `fp_sub`, `fp_mul` (leaf functions), `fp2_mul`, `fp2_add`, `fp2_sub`,
  `fp2_mul_xi`, `fp6_mul`, `fp6_add`, `fp6_sub`, `fp12_mul`, `walk_step`.
  No function contains a loop or a conditional branch; conditional selection
  is by CSEL/masks.  Hence each function executes each of its instructions
  exactly once per call, and its executed count is its static count plus the
  executed counts of its callees, independent of the data.
* `arithmetic.disasm.txt` is the disassembly (`objdump -d
  --no-show-raw-insn`) of that file compiled with Apple clang 17.0.0, target
  arm64, flags `-O2 -fno-vectorize -fno-slp-vectorize`.  It uses only
  general-purpose registers (the recount asserts that no v/q/d/s/z register
  occurs) and its only control transfers are `bl`, tail-call `b` to another
  function, and `ret` (asserted).  `fp_mul` contains MRS/MSR NZCV pairs
  emitted by the compiler to save/restore the carry flag; they are ordinary
  scalar instructions and are counted.
* `experiment.py` section 6 parses that text inside the verifier container
  and recomputes the table below; `development-evidence.json` records the
  compile command and the same numbers.  The compiler is a convenience for
  producing a concrete listing: the bound is a property of the listing, which
  is itself a valid implementation, not of any optimisation assumption.

| function | static | executed | composition |
|---|---:|---:|---|
| fp_add | 35 | 35 | leaf: 4 limb add, 4 limb subtract p (constants by MOV/MOVK), 4 CSEL, loads/stores, RET |
| fr_add | 34 | 34 | same modulo r |
| fp_sub | 35 | 35 | leaf: 4 limb subtract, masked add of p |
| fp_mul | 218 | 218 | leaf: 4 rows x (4 MUL+4 UMULH+carry adds; q = t0*n0; 4 MUL+4 UMULH+carry adds), final subtract and CSEL |
| fp2_mul | 44 | 873 | 3 fp_mul + 2 fp_add + 3 fp_sub + 44 |
| fp2_add / fp2_sub | 15 | 85 | 2 leaf + 15 |
| fp2_mul_xi | 49 | 399 | 9 fp_add + 1 fp_sub + 49 |
| fp6_mul | 101 | 7,412 | 6 fp2_mul + 9 fp2_add + 6 fp2_sub + 2 fp2_mul_xi + 101 |
| fp6_add / fp6_sub | 19 | 274 | 3 fp2 add/sub + 19 |
| fp12_mul | 58 | 24,044 | 3 fp6_mul + 2 fp6_add + 2 fp6_sub + 1 fp2_mul_xi + 3 fp2_add + 58 |
| walk_step | 44 | **24,156** | fp12_mul + 2 fr_add + 44 |

Check of the compositions: 3*218+2*35+3*35+44 = 873; 9*35+35+49 = 399;
6*873+9*85+6*85+2*399+101 = 7,412; 3*7,412+2*274+2*274+399+3*85+58 = 24,044;
24,044+2*34+44 = 24,156.  In field operations one Fp12 multiplication is 54
Montgomery multiplications (11,772 instructions) and 286 Fp additions or
subtractions (10,010), plus 2,262 instructions of calls, argument set-up
and returns.

### Why 218 is a credible Montgomery multiplication count

Independent sanity check not relying on the compiler: the organiser's cost
model measures 179 instructions for a 4x64-bit-class secp256k1 multiplication
with special-form reduction and 383 for the six-limb blst Montgomery
multiplication.  A four-limb coarsely-integrated Montgomery product needs
16+16 word products for a*b and q*p (64 MUL/UMULH) plus 4 MUL for the q
values, i.e. 68 multiply instructions, about 2-3 carry instructions per word
product, 8 loads, 4 stores and the final conditional subtraction; 218 is in
the expected range between the two organiser measurements scaled by limb
count (383 * (4/6)^2 = 170 is a lower reference for the quadratic part).

### Correctness of the leaf arithmetic

p < 2^254 and r < 2^254, R = 2^256.  `fp_mul`: after each row the running
value t (five words plus one carry bit) equals (previous + a*b_i + q*p)/2^64
with q chosen so the low word vanishes.  If the previous value is < 2p then
the new one is < (2p + (2^64-1)p + (2^64-1)p)/2^64 < 2p, so t < 2p after
every row by induction (a < p, b_i,q <= 2^64-1);
so the final value is < 2p < 2^255, fits
four words, and one conditional subtraction makes it canonical.
`rho.mont_mul_words` is the word-for-word Python mirror with these bounds as
assertions; the experiment checks it on 544 full-width edge and random cases
for both moduli.  Additions cannot overflow four words because 2p < 2^255.
The Karatsuba tower formulas are the textbook ones; all reads of the inputs
precede the first write to the output, and the experiment checks them
against the protected polynomial basis.

## The step loop

Around `walk_step` the loop does: call (1), two-word step counter increment
and comparison with N (ADDS/ADC/CMP/CCMP/branch, <= 8), branch on the
distinguished flag (<= 2), swap of the in/out state pointers (<= 3), passing
four arguments (<= 4).  That is below 20 instructions; we charge

    S = 24,576 = 24,156 + 420 per step,

a uniform worst-case bound including table-record addressing (the record is
448 bytes; index*448 is one multiply-add inside the 44) and all memory access.
The 1024-record table (458,752 bytes) and two states are ordinary 64-bit
addressed memory.  Non-distinguished steps do nothing else.

## Distinguished-point table

A record is a state: 384 + 32 + 32 = 448 bytes.  At a distinguished step the
algorithm scans all stored records (at most T = 2^29), comparing 48 words
each: <= 48*(2 loads + compare + branch) + 16 < 512 instructions per record
(early exit is allowed but not assumed).  At most T appends happen, each
preceded by one scan, plus at most one final scan: fewer than (T+1)*T*512 <
2^67.1 instructions in the worst case, plus T copies of 56 words (< 2^37).
This crude quadratic scan is affordable because distinguished states are
2^-100 rare; a hash table is unnecessary.  Addresses are below 2^38, so
64-bit addressing suffices and no wide-address cost arises.

## Fixed work F = 2^70

| item | bound (instructions) |
|---|---:|
| parse the fixed 12 x 64-hex JSON form, canonicity checks, basis change, conversion to Montgomery form (24 fp_mul and ~36 additions) | < 2^20 |
| subgroup check z^r = 1: <= 508 fp12_mul | < 2^24 |
| salt and 2,050 scalars, <= 1024 draws of 254 bits each, 64 instructions per bit | < 2^36 |
| 1,025 records: two exponentiations (<= 508 fp12_mul each) and one product, i.e. <= 1,017 fp12_mul each; conversion of scalars | 1,025*1,017*24,044 + 2^24 < 2^35 |
| all distinguished-table scans and appends (worst case) | < 2^67.2 |
| recovery: one inverse mod r by exponentiation (<= 508 Montgomery products mod r, each <= 218+slack 256), two products, subtraction | < 2^18 |
| verification g^x = z (<= 508 fp12_mul), 48-word compare, scalar output | < 2^24 |
| allocation, code/data loading, shutdown, all remaining control | < 2^32 |
| **sum** | < 2^67.3 < **2^70** |

Montgomery arithmetic modulo r for the recovery uses the same `fp_mul`
algorithm with the constants of r (n0 and R^2 mod r are public profile-derived
constants; computing them takes < 2^16 instructions, inside the allowance).

## Peak storage

Distinguished table 2^29 * 448 = 2^37.81 bytes = 240,518,168,576.  Increment
table 458,752 bytes.  States, scratch, stack (the deepest call chain keeps
fewer than 4 KiB of temporaries), constants, code (< 64 KiB), randomness
buffer: < 1 MiB together.  Total < 240,520,000,000 < **2^38 bytes**.
Advice bytes: zero.  Memory is reported only; it does not enter the score.

## Total

    N = ceil(13*ceil(sqrt(r))/10)
    R_total <= N*24,576 + 2^70
    W = R_total/20,000,  log2 W <= 127.474108.

The bound is a cap valid for every outcome (success, cap reached, table
full, degenerate, verification failure, rejection-sampling failure).
