# Finite scalar implementation and work bounds

These are conservative **analytical scalar AArch64 instruction-equivalent
bounds**, not timings and not prices inferred from measurements. The
implementation below uses only 64-bit scalar registers and ordinary
addresses. Python is an executable mathematical reference, not the priced
implementation. There is one processor. The walk, partition, and cap are
the baseline's; the detection method (distinguished points instead of
Floyd) and the Fp12 arithmetic (tower Karatsuba instead of schoolbook
polynomial convolution) are new, and every changed number is derived
below from explicit cell counts.

## Word operations and modular arithmetic

Let B=2^64 and R=B^4. Store a residue as four little-endian words. For
either modulus m=p or m=r, m is odd and <2^254. Store multiplication
operands in Montgomery form aR mod m. Public constants may be loaded;
generating n0=-m^-1 mod B and R^2 mod m instead takes six Newton word
iterations and at most 512 modular doublings, included in setup.

`arithmetic.c` specifies the product and REDC algorithms with fixed
bounded loops and alias-safe scratch, following the baseline's loop
structure. A 64x64 product is MUL plus UMULH (two scalar instructions);
128-bit sums/carries use ADDS/ADC, comparisons/subtractions use SUBS/SBC
and CSET, and logical selection uses ordinary word masks. No 128-bit
unit-cost arithmetic is assumed. Indexes are at most 9 and can be
unrolled into immediate offsets. Loads, stores, address formation, loops
and calls are included below. The implementation can be emitted directly
with these expansions; no particular C compiler optimization is a
premise. The multiply is register-resident: the accumulator, the three
operands, and all temporaries live in registers for the whole multiply
(the explicit 28-register allocation is stated below), so no load or
store occurs inside a cell. Operand loads and result stores are charged
once, and spills, calls, and loop control once, in the table below.

Per modular multiplication, the following cell counts suffice:

| Work | Cells | Maximum scalar instructions each | Bound |
|---|---:|---:|---:|
| Schoolbook product cell (register-resident: 5 core, charged at 8) | 16 | 8 | 128 |
| Initialise 9-word accumulator (one MOV per word) | 1 | 9 | 9 |
| Product carry-propagation cell | 14 | 4 | 56 |
| REDC q = t[i]*n0 + pass setup | 4 | 4 | 16 |
| REDC product cell (register-resident: 5 core, charged at 8) | 16 | 8 | 128 |
| REDC carry-propagation cell | 14 | 4 | 56 |
| Final conditional subtract + select | 1 | 32 | 32 |
| Operand loads + result stores | 1 | 24 | 24 |
| Call/return, spills, loop control | 1 | 48 | 48 |
| **Total** | | | **497 < M = 512** |

A product cell, with all operands resident in registers, executes five
core instructions: MUL and UMULH for the 64x64 product, then ADDS (low
word plus accumulator word), ADCS (plus carry), and ADC (high-word
carry-out) to add the 128-bit product, the accumulator word, and the
incoming carry. No load or store occurs inside a cell; the table
charges 8 per cell, leaving 3 per cell of headroom for index and
addressing. The register allocation is explicit: accumulator t[0..8]
(9 registers), operands a, b, m (12), n0, q, carry (3), product
temporaries (2), index/addressing (2), totalling 28 registers, within
the 30 usable general-purpose registers of a leaf function. The 9-word
accumulator is zero-initialised once (one MOV per word, 9). A
propagation cell adds one accumulator word with carry and controls its
fixed loop: 4 suffice.
The final conditional subtract compares four words, subtracts with
borrow, and selects with masks: 32 suffice. There are 14 propagation
cells per stage (positions i+4..8 across the four passes, counted with
their actual bounds), not 28: each pass propagates at most 5-k cells
for pass k, and the table charges 14 per stage as an upper envelope.

Correctness/overflow: each product-plus-two-words is at most B^2-1. Full
product a*b<R^2 is represented in 9 words. REDC cancels the next low
word at each pass; its final U=(a*b+q*m)/R is <2m since a,b<m<R and
q<R. Thus U<2^255, the ninth word is zero and one subtraction of m makes
the answer canonical. Intermediate carry propagation is retained through
word 8. Four-limb addition cannot overflow 256 bits because 2m<2^255.

Independent cross-check (support, not pricing): the repo's own pinned
measurement (tracks/docs/COST_MODEL.md) is **383 instructions for
`blst_fp_mul`, a 6-limb multiply**. A 4-limb schoolbook multiply has
16/36 = 44% of the multiply cells, so the blst-scaled estimate is
~200-250 instructions; M <= 512 carries ~2x margin over that too. The
bound is an analytical cell count; no measurement prices the arithmetic.

Modular addition, subtraction, negation and comparison each cost at most
**64** instructions: four carry/borrow cells (4 ADDS/ADC or SUBS/SBC), a
four-word conditional correction/selection (4 SUBS/SBC + 4 CSEL/mask),
loads/stores (~12) and control (~8) total ~32; the claimed 64 is 2x
margin. `add_mod` is supplied in arithmetic.c; subtraction subtracts b
with borrow and conditionally adds m, and negation subtracts from m and
masks the zero case. A four-word copy costs <=32. These bounds also
cover scalar coefficients modulo r. Inverse of a nonzero residue is
exponentiation to m-2 with <=508 modular multiplications and <=2^16
scalar overhead, hence <2^19 instructions. Representation conversions
each cost one Montgomery multiplication; explicit setup/cleanup allows
for them.

## Tower Fp12 multiplication

Use the tower Fp subset Fp2 = Fp[u]/(u^2+1) subset Fp6 = Fp2[v]/(v^3-xi),
xi = 9+u, subset Fp12 = Fp6[t]/(t^2-v), with canonical Montgomery
coefficients throughout. The embedding u = w^6-9, v = w^2, t = w into
the protected polynomial basis is proved in proof.md. The input z is
converted once: the conversion is a fixed 12x12 linear map over Fp,
costing at most 144 modular multiplications and 132 modular additions,
hence <2^17 instructions, charged in the fixed budget. (The direct
formula in rho.py is cheaper; the bound charges the generic map.)

Karatsuba formulas (all counts verified by the instrumented reference
in rho.py, re-checked by experiment.py):

Fp2 multiply (complex Karatsuba): for (a0+a1*u)(b0+b1*u), compute
m0 = a0*b0, m1 = a1*b1, m2 = (a0+a1)(b0+b1); result
(m0-m1) + (m2-m0-m1)*u. **3 Fp muls + 5 Fp add/subs.**

Multiply by xi = 9+u: (a0+a1*u)(9+u) = (9*a0-a1) + (a0+9*a1)*u, with
9x = 8x+x via three doubles and one add. **10 Fp add/subs, 0 muls.**

Fp6 multiply: for a = a0+a1*v+a2*v^2, b = b0+b1*v+b2*v^2, v^3 = xi,
compute v0 = a0*b0, v1 = a1*b1, v2 = a2*b2,
v3 = (a0+a1)(b0+b1)-v0-v1, v4 = (a0+a2)(b0+b2)-v0-v2,
v5 = (a1+a2)(b1+b2)-v1-v2; result c0 = v0+xi*v5, c1 = v3+xi*v2,
c2 = v4+v1. **6 Fp2 muls + 15 Fp2 add/subs + 2 xi-muls.**

Multiply by v: (a0+a1*v+a2*v^2)*v = xi*a2 + a0*v + a1*v^2. **1 xi-mul.**

Fp12 multiply (Karatsuba over Fp6): for (a0+a1*t)(b0+b1*t), t^2 = v,
compute v0 = a0*b0, v1 = a1*b1,
v2 = (a0+a1)(b0+b1)-v0-v1; result (v0+v*v1, v2).
**3 Fp6 muls + 5 Fp6 add/subs + 1 mul-by-v.**

Total per Fp12 multiply, by instrumented count: **54 Fp multiplications
+ 280 Fp add/subs.** The data-movement lump below covers call/return for
the 54 Fp-multiply calls, stack traffic for the ~40 live Fp2/Fp6
temporaries, and addressing. Rough itemisation of the ~4,000 estimate:
54 Fp-multiply calls at ~20 instructions each ≈ 1,100; spill/fill of
the ~40 live temporaries ≈ 650; addressing and loop control ≈ 500;
the remaining ~1,750 covers argument/result shuffling and unitemised
movement. The charged lump of 8,192 is 2x the estimate:

    G = 54*512 + 280*64 + 8192 = 27,648 + 17,920 + 8,192 = 53,760.

Squaring may use this same multiplication. There is no double counting
of a normalization price. Tower scratch occupies <2 KiB and is reused,
never recursively allocated. Cross-check: the pinned blst dense Fp12
multiply (6-limb base field, assembly) measures 22,996 instructions;
our bound of 53,760 for a 4-limb field is 2.3x that measurement, i.e.
conservative, not aggressive. No measurement prices the arithmetic.

Fp12 exponentiation by a <2^254 scalar uses <=508 such multiplications,
plus <=2^20 instructions for bits, copies and loop overhead (the
baseline's allowance, kept).

## Hashing and walk state

Hash SHA-256(salt || encode(Y)), taking the low five bits of digest byte
zero as the table index. Salt is 32 bytes. encode(Y) is the 384-byte
concatenation of the twelve tower coefficients as **canonical Montgomery
residues**, each in fixed big-endian 32-byte form, in the fixed order
(a0[0][0], a0[0][1], a0[1][0], a0[1][1], a0[2][0], a0[2][1],
a1[0][0], a1[0][1], a1[1][0], a1[1][1], a1[2][0], a1[2][1]).
This internal encoding is injective (the tower representation is unique
and Montgomery form is canonical), and explicitly differs from the
external standard-coefficient input encoding. No conversion out of
Montgomery representation is needed on each step. Scalar software
SHA-256 uses seven 64-byte padded blocks for this 416-byte message.

For clarity, this is ordinary SHA-256: initialize the standard eight IV
words, expand the 16 input words to 64 with the standard small-sigma
rotates/shifts and additions, then perform 64 rounds of the standard Ch,
Maj, large-sigma and modular-2^32 additions per block and feed forward.
All its 32-bit instructions are scalar AArch64 W-register operations.
Constants are the public standard SHA-256 IV and 64 K words. A rotate
costs one ROR.

Per block: message schedule 48 words x 16 = 768 (sigma0 5 + sigma1 5 +
3 adds + store/control ~3 per word, with the 16-word ring buffer held
in registers); 64 rounds x 32 = 2,048 (Ch 3 + Maj 3 + Sigma0 5 + Sigma1
5 + adds 6 + W/K loads 2 + rotation/control ~8); feed-forward + state
32; total 2,848 per block. Seven blocks give 19,936. Serialization
(104 words) + padding ≈ 1,024. Setup/extraction ≈ 512. Total ≈ 21,472
< **H = 2^15 = 32,768** (1.5x margin). A single input byte is never an uncharged oracle query. This
is a finite deterministic implementation; random-mapping behavior of
the walk is a disclosed heuristic, not a free ideal hash implementation.
The distinguished test (word ANDs/compares on digest bytes 1..14) and
the slot-index extraction (digest bytes 15..31) use the same digest and
are charged inside A below, not as extra hashes.

A step reads one of 32 records, multiplies Y by its stored group
element, adds the two 4-word coefficients modulo r, tests the
distinguished condition, and updates the state. The table record is 448
bytes (384+32+32). Record read (112 words, loads + addressing) ≈ 224;
two mod-r adds 2x64 = 128; distinguished test ≈ 32; state write-back
≈ 128; three-word step-counter increment/compare and control/spills
≈ 256. Total ≈ 768 < **A = 1,024** beyond G+H. The distinguished-point
hash-table lookup is not a per-step cost: it happens only at
distinguished states and is charged in the fixed budget below.

Therefore

    S = G + H + A = 53,760 + 32,768 + 1,024 = 87,552

is a uniform worst-case bound per step, including memory access. All
loops have concrete integer bounds; the result is a cap, not an expected
work estimate. No cost depends on an unrepresented enormous address.

## Distinguished-point table

Open-addressed hash table: 2^21 slots x 512-byte records = 2^30 bytes.
A record holds the 384-byte group element, the two 32-byte scalar labels,
and 64 bytes of slot metadata including an 8-byte fingerprint (digest
bytes 15..22) used to reject non-matching probes without a full
384-byte compare; a fingerprint match always triggers the full compare,
so correctness never depends on the fingerprint. The slot index comes
from the same digest as the partition index and distinguished test.
At most 2^20 records are stored, so the load never exceeds 0.5.

Each lookup probes linearly with a cap of 256 probes; exceeding the cap
reports failure (probability quantified in proof.md: negligible at load
<= 0.5). Each lookup costs at most 2^13 instructions: 256 fingerprint
probes x 16 instructions, one full element compare (<=512), one record
write (<=512), and control (<=512), totaling <=5,632 < 8,192. The number
of lookups is uniformly bounded: every lookup either hits (recovery or
degenerate failure, ending the walk) or misses (inserting a record,
bounded by the 2^20 record cap), so there are at most 2^20+1 < 2^21
lookups. Total lookup work <= 2^21 x 2^13 = 2^34, charged in the fixed
budget. Lookup work is 2^-110 of the total R; it cannot move the score.

## One-target setup, termination and output

Use exactly one trial and at most K steps. No retry is implicit. Failure
by the cap, rejection-sampling failure, degenerate collision, probe-cap
or record-cap overflow, and failed output verification terminate within
the same cap and count as failure.

The following additional fixed **F=2^40** instructions cover both tracks:

| Item | Bound |
|---|---:|
| Input validation + z^r = 1 check in tower form after conversion (<=508 GT products + parsing) | < 2^26 |
| Fresh randomness (baseline's charge: 66 scalars, <=1024 draws, 64 instr/bit) | < 2^31 |
| Table build: 33 exponentiations, <=33*1017 GT products x G | 1.81x10^9 < 2^31 |
| Basis conversion of z (144 Fp muls + 132 Fp adds, one-time) | < 2^17 |
| Recovery: inversion mod r + final g^x = z verification | < 2^26 |
| Pairing output (pairing track; baseline's <2^35 bound kept) | < 2^35 |
| Distinguished-point table lookups (<=2^21 x 2^13) | 2^34 |
| Parsing/encoding/allocation/shutdown (baseline's allowance) | <= 2^32 |
| **Total** | **< 2^36 < F = 2^40** |

Table-build detail: build all 32 M_j=g^c_j*z^d_j, plus the initial
g^a*z^b, charging two exponentiations and one multiplication per
record, <=33*1017 GT products. Charge all searched precomputation to
this target, never to future inputs. Generate/read at most
66*1024*256+256 fresh bits, charging even 64 instructions per bit,
<2^31. Only one candidate and scalar at a time are retained. Physical
entropy timing is not the instruction metric; these are the algorithm's
ordinary fresh independent coins, not uncharged group/advice oracles.
The above disjoint quantities sum to <6.0x10^10 < 2^36; F=2^40 leaves
additional slack. Four further small items are not given their own rows
and are absorbed by the slack in the <2^36 total (F = 2^40 is a further
16x above it): exponentiation loop overhead (<=2^20 per exponentiation
over ~68 exponentiations), Montgomery-form conversions of the g/z
coefficients and R^2 generation (~24 multiplications plus setup), the
one-time tower conversion of the public generator g (~2^17, as for z),
and zeroing the 2^21-slot hash table (~2^28 word stores); even with
generous allowances for all four, the total stays below 2^36. This
allowance is derived from finite algorithms, not from Python timings.

## Pairing output and its explicit cost

After verifying recovered x with g^x=z, return P0 and Q=[x]Q0. Use binary
double-and-add in affine Fp2; at most 254 doubles and 254 additions. No
pairing computation is needed by the attack: its verified scalar and
bilinearity establish the output relation. The organizer's pairing check
is excluded by the protected rules. The derivation below is the
baseline's, kept with its bounds: under this package's tighter
per-operation bounds (M <= 512, add <= 64) every step of it holds a
fortiori.

For u^2=-1, multiplication (a+bu)(c+du)=(ac-bd)+(ad+bc)u takes four Fp
multiplications and two additions, <2^15 instructions including
loads/control. Fp2 addition/subtraction takes two Fp operations plus
copies, <=1024. Inversion is (a-bu)/(a^2+b^2), using two squarings, one
Fp inversion by p-2, two scaling multiplications and adds/negation,
<2^23 instructions. For nonzero a+bu the norm is nonzero: p=3 mod 4, so
-1 is not a square.

For nonexceptional affine point addition, compute
lambda=(y2-y1)/(x2-x1), X=lambda^2-x1-x2,
Y=lambda*(x1-X)-y1.
For doubling lambda=3*x1^2/(2*y1), X=lambda^2-2*x1,
Y=lambda*(x1-X)-y1.
Infinity, equal x/opposite y and y=0 are checked first and handled by
the usual infinity or doubling case; they do not call inverse at zero.
At most one Fp2 inverse, 8 Fp2 multiplications and 30 Fp2 add/sub/copy
operations plus <2^16 scalar control suffice for any point operation:
< **2^24** instructions.

508 operations for [x]Q0 and 508 for the optional check [r]Q=infinity use
<2^34 instructions. Curve-equation evaluation, canonical conversion and
serialization, infinity and P0 handling fit in <2^30. The total is <2^35,
as included in F. An order-r multiple of Q0 is on the required subgroup;
the optional explicit checks add evidence without weakening that
argument. x=0 returns Q=infinity, which the protected profile permits
exactly for z=1.

## Peak storage

Distinguished-point table: 2^21 slots x 512 bytes = 2^30 bytes. 32
table records: 14,336 bytes. Walk states and group scratch: <8192. SHA
scratch/encoding: <4096. Scalar/curve temporaries, random draw and loop
counters: <16,384. Profiles/constants: <8192. The scalar implementation,
SHA, arithmetic, parser and output code are finite loops; 128 KiB
suffices even at 32,768 scalar instructions of 4 bytes each. Reserve
another 256 KiB for stacks, code-loading scratch, allocation metadata
and output. Non-table storage sums to <2^20 bytes (the baseline's 1 MiB
bound, reused). Peak total < 2^30 + 2^20 < **2^31 bytes**. Table slot
addresses are ordinary addresses below 2^30; no enormous address is
assumed. Iteration counts do not grow live storage beyond the capped
table; no trace or random tape is retained. This is the algorithm's
total code/data/processor storage on a scalar reference machine, not
the excluded organizer or Python test interpreter. Advice bytes are
zero (log2(max(1,advice))=0).

## Total and score

    K = 2*ceil(sqrt(r))
    R_total = K*87,552 + 2^40
    W = R_total/20,000
    submitted score = ceil(10^6*log2(W))/10^6

Memory log2 bytes <=31. Preprocessing is <=2^32 instruction equivalents
(table build plus randomness generation), or 2^32/20,000 reference
units, already included in F. No history, earlier score, former rejected
estimate, or existing table is used as a premise.
