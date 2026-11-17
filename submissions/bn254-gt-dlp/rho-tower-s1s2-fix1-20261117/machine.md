# Finite scalar implementation and work bounds

These are conservative **analytical scalar AArch64 instruction-equivalent
bounds**, not timings. The implementation below uses only 64-bit scalar
registers and ordinary addresses within a 1 MiB allocation. Python is an
executable mathematical reference, not the priced implementation. There is
one processor. Where the accepted `rho-baseline-v2-20261110` used
deliberately loose bounds, this file derives tighter ones from explicit
per-operation counts; every count below is either derived line by line or
anchored to a published measurement cited in the track cost model.

## Word operations and modular arithmetic

Let B=2^64 and R=B^4. Store a residue as four little-endian words. For either
modulus m=p or m=r, m is odd and <2^254. Store multiplication operands in
Montgomery form aR mod m. Public constants may be loaded; generating n0=-m^-1
mod B and R^2 mod m instead takes six Newton word iterations and at most 512
modular doublings, included in setup.

`arithmetic.c` specifies the product and REDC algorithms with fixed bounded
loops and alias-safe scratch. A 64x64 product is MUL plus UMULH (two scalar
instructions); 128-bit sums/carries use ADDS/ADC, comparisons/subtractions
use SUBS/SBC and CSET, and logical selection uses ordinary word masks. No
128-bit unit-cost arithmetic is assumed. Indexes are at most 9 and can be
unrolled into immediate offsets. Loads, stores, address formation, loops and
calls are included below. The implementation can be emitted directly with
these expansions; no particular C compiler optimization is a premise.
All loops below have fixed public bounds (4 limbs, 9 accumulator words),
so the priced expansion is fully unrolled: every count below is a
straight-line count, with loop control charged only through the explicit
frame allowances.

Per modular multiplication, the following counts suffice, derived cell by
cell from the fixed loop structure of `arithmetic.c` (mont_mul). A product
cell (one (i,j) pair in either the schoolbook or the REDC phase) performs:
load the operand word (LDR), multiply low (MUL), multiply high (UMULH),
load the accumulator word (LDR), add the low product with carry (ADDS),
add the high product with carry (ADC), add the incoming carry (ADDS),
propagate it (ADC), store the accumulator word (STR) — 9 instructions.
There are 16 schoolbook cells and 16 REDC cells. The four a[i] operand
words are loaded once and held in registers (4). Each q = t[i]*n0 costs
one load plus one multiply (2 each, 8 total); n0 is held in a register,
loaded once inside the frame allowance. Carry propagation after each row
runs over the remaining accumulator words: in the schoolbook phase
5+4+3+2 = 14 cells, and the same 14 in the REDC phase, 28 cells in total;
each propagation cell is one load, one add-with-carry and one store
(3 each, 84 total). Zeroing the 9-word accumulator costs 5 (four STP
pairs plus one STR). The final conditional subtraction and selection costs
16 (4 loads of the candidate, 4 subtract-with-borrow, 2 mask formation,
4 word selects, 2 stores). Frame, spill and call allowance: 16.

| Work | Count | Scalar instructions each | Bound |
|---|---:|---:|---:|
| Schoolbook product cell | 16 | 9 | 144 |
| REDC product cell | 16 | 9 | 144 |
| Hoisted a[i] loads | 4 | 1 | 4 |
| q = t[i]*n0 | 4 | 2 | 8 |
| Carry-propagation cell | 28 | 3 | 84 |
| Accumulator zeroing | 1 | 5 | 5 |
| Final conditional subtraction and select | 1 | 16 | 16 |
| Frame, spill, call allowance | 1 | 16 | 16 |
| Total derived | | | 421 |
| Charged | | | **448** |

The charged bound of 448 carries 27 instructions (6%) of margin over the
derived 421. The blst measurement in the track cost model (383 instructions
for a 72-cell 6-limb Montgomery multiply, about 5.3 per cell) is a
hand-scheduled assembly figure; our bound charges every load and carry
separately at 9 per product cell, which is why it sits above that
measurement. The gap is in the safe direction for an upper bound:
over-counting cannot inflate the score.

Correctness/overflow: each product-plus-two-words is at most B^2-1. Full
product a*b<R^2 is represented in 9 words. REDC cancels the next low word at
each pass; its final U=(a*b+q*m)/R is <2m since a,b<m<R and q<R. Thus U<2^255,
the ninth word is zero and one subtraction of m makes the answer canonical.
Intermediate carry propagation is retained through word 8. Four-limb
addition cannot overflow 256 bits because 2m<2^255.

Modular addition and subtraction each cost at most **32** instructions,
derived at the register level. For add_mod (a, b < m, so s = a+b < 2m):
load the two operands (8 words, 4 LDP); add chain (ADDS + 3 ADC, 4);
load the modulus (4 words, 2 LDP); subtract chain (SUBS + 3 SBC, 4)
computing d = s-m and leaving the borrow in the flags; select per word
(4 CSEL on the borrow flag, choosing d where s >= m and s otherwise);
store the result (4 words, 2 STP). Total 20 instructions. sub_mod is
symmetric (subtract chain, add chain, select on the borrow flag), also
20. The charged bound of 32 carries 12 instructions (60%) of margin for
register pressure and spill. Negation (subtract from m, mask the zero
case), comparison (4 CMP plus control) and a four-word copy (2 LDP +
2 STP) are all bounded by 32 as well. These bounds also cover scalar
coefficients modulo r. Inverse of a nonzero residue is exponentiation to
m-2 with <=508 modular multiplications and <=2^16 scalar overhead, hence
508*448 + 2^16 = 293,120 < 2^19 instructions. Representation conversions
each cost one Montgomery multiplication; explicit setup/cleanup allows
for them.

Aliasing. Every function in arithmetic.c permits `out` to alias any
input operand: each function reads all of its inputs (or completes all
input reads) before writing any output word. In particular fp2_mul_xi
computes both result components into scratch before writing either, so
the fp6_mul call fp2_mul_xi(t, t, m) is correct. An earlier revision
wrote out[0] before reading a[0] for out[1], which corrupted exactly
that aliased call; the defect is fixed here and is covered by the
experiment's aliased-call test (64 random aliased inputs, all matching
the unaliased result) and by the v*v^2 = xi vector below.

## The tower: Fp2, Fp6, Fp12

Work in the tower basis, not the profile's polynomial basis:

- Fp2 = Fp[u]/(u^2+1). Irreducible because p = 3 mod 4, so -1 is a
  non-square in Fp (Legendre symbol -1, verified: (-1)^((p-1)/2) = -1).
- Fp6 = Fp2[v]/(v^3-xi) with xi = 9+u. Irreducible because xi is a
  non-cube in Fp2, verified by direct exponentiation:
  xi^((p^2-1)/3) = 21888242871839275220042445260109153167277707414472061641714758635765020556616
  (an element of Fp, a primitive cube root of unity), which is not 1.
  Note xi = 1+u would NOT work: 1+u is a cube in Fp2 for this prime, so
  the BLS12-381 tower constant cannot be reused here.
- Fp12 = Fp6[w]/(w^2-v). Irreducible because v is a non-square in Fp6,
  verified by direct exponentiation: v^((p^6-1)/2) = -1 in Fp6.
  Independently, xi is a non-square in Fp2, and the norm N(v) = xi, so
  v = y^2 would force xi = N(y)^2 to be a square in Fp2 — a second,
  structural proof that v is a non-square.

The tower is a field with p^12 elements, hence isomorphic to the profile's
polynomial-basis Fp12 by the uniqueness of finite fields. The basis
conversion is a 12x12 Fp-linear map (see below); the walk, the hash
partition and the collision statistics are unaffected by the choice of
basis, since the group operation and the encoding are both basis-covariant
and the encoding remains injective.

### Fp2 multiplication (Karatsuba)

For a = a0+a1*u, b = b0+b1*u with u^2 = -1:
t0 = a0*b0, t1 = a1*b1, t2 = (a0+a1)*(b0+b1);
c0 = t0-t1, c1 = t2-t0-t1. That is 3 Fp multiplications and 5 Fp
additions/subtractions. Cost: 3*448 + 5*32 = 1,344 + 160 = 1,504.
Squaring costs no more than multiplication (charged as a multiplication
throughout). An Fp2 addition/subtraction costs 2*32 = 64.

### Fp6 multiplication (Karatsuba)

For a = a0+a1*v+a2*v^2, b = b0+b1*v+b2*v^2 with v^3 = xi:
v0 = a0*b0, v1 = a1*b1, v2 = a2*b2,
v3 = (a0+a1)*(b0+b1), v4 = (a1+a2)*(b1+b2), v5 = (a0+a2)*(b0+b2);
c0 = v0 + xi*(v4-v1-v2), c1 = v3-v0-v1 + xi*v2, c2 = v5-v0-v2+v1.
That is 6 Fp2 multiplications, 15 Fp2 additions/subtractions, and 2
multiplications by xi. Multiplication by xi = 9+u:
(x0+x1*u)*(9+u) = (9*x0-x1) + (x0+9*x1)*u, where 9*x = 8*x+x costs
3 doubles + 1 add = 4 Fp additions; each xi-multiplication costs
10 Fp additions. Cost per Fp6 multiplication:
6*1,504 + 15*64 + 2*320 = 9,024 + 960 + 640 = 10,624.
An Fp6 addition/subtraction costs 6*32 = 192.

### Fp12 multiplication (Karatsuba)

For a = a0+a1*w, b = b0+b1*w with w^2 = v:
t0 = a0*b0, t1 = a1*b1, t2 = (a0+a1)*(b0+b1);
c0 = t0 + v*t1, c1 = t2-t0-t1. That is 3 Fp6 multiplications, 5 Fp6
additions/subtractions, and 1 multiplication by v. Multiplication by v
permutes the Fp6 coefficients and applies one xi-multiplication:
v*(d0+d1*v+d2*v^2) = xi*d2 + d0*v + d1*v^2, costing 10 Fp additions
plus coefficient moves. Cost per Fp12 multiplication:
3*10,624 + 5*192 + 320 = 31,872 + 960 + 320 = 33,152.
Coefficient moves, copies and call overhead fit within a further 512
instructions (there are 54 Fp multiplications, 280 Fp additions and under
100 word moves in the whole routine; the per-operation bounds above already
include their local loads/stores). The total is 33,664.

Cross-checks on the count, two independent ways:
- Flat count: 3*6*3 = 54 Fp multiplications, and Fp additions:
  18 Fp2 multiplications * 5 = 90, 3 Fp6 multiplications * 15 Fp2 adds
  = 45 Fp2 adds = 90 Fp adds, 5 Fp12-level Fp6 adds = 15 Fp2 adds
  = 30 Fp adds, 7 xi/v multiplications * 10 = 70; total 280 Fp
  additions. 54*448 + 280*32 = 24,192 + 8,960 = 33,152, agreeing
  exactly with the layered count above; plus 512 overhead = 33,664.
  The experiment's harness measures the same counts from the compiled
  code (54 mont_mul, 171 add_mod, 109 sub_mod) and derives this total
  from them; see "Direct validation of the C code" below.
- Published anchor: blst measures 22,996 instructions for a dense Fp12
  multiply on BLS12-381 (381-bit field, same tower structure, same
  Karatsuba counts), per the track cost model. Our derived bound of
  33,664 sits about 46% above that measurement for the larger field.
  The gap is explained, not hidden: blst is hand-scheduled assembly
  that fuses loads and carries (about 5.3 instructions per Montgomery
  cell), while our bound charges every load and carry separately
  (9 per product cell, 3 per propagation cell) and prices a fully
  scalar, un-fused expansion. The gap is in the safe direction for an
  upper bound: over-counting cannot inflate the score, and the charged
  G' = 40,000 below is higher still.

Charge **G' = 40,000** per Fp12 multiplication: 19% above the derived
bound of 33,664. Squaring may use this same multiplication. Polynomial
scratch occupies <2 KiB and is reused, never recursively allocated.

Fp12 exponentiation by a <2^254 scalar uses <=508 such multiplications,
plus <=2^20 instructions for bits, copies and loop overhead.

## Basis conversion (one-time)

The profile supplies z and g in polynomial-basis coefficients; the walk
uses tower-basis coefficients. The conversion in either direction is a
fixed 12x12 Fp-linear map: at most 144 Fp multiplications and 132 Fp
additions, under 144*448 + 132*32 = 68,736 instructions, once per input
(and once for g). This is charged inside the setup allowance F below; it
is not a per-step cost. No conversion is needed on any walk step: the
table, the states and the hash encoding all live in the tower basis.

## Hashing and walk state

Hash SHA-256(salt || encode(Y)), taking the low five bits of digest byte zero
as the table index. Salt is 32 bytes. encode(Y) is the 384-byte concatenation
of its twelve **canonical Montgomery residues in the tower basis**, each in
fixed big-endian 32-byte form. This internal encoding is injective, and
explicitly differs from the external standard-coefficient input encoding.
No conversion out of Montgomery representation is needed on each step.
Scalar software SHA-256 uses seven 64-byte padded blocks for this 416-byte
message.

For clarity, this is ordinary SHA-256: initialize the standard eight IV
words, expand the 16 input words to 64 with the standard small-sigma rotates/
shifts and additions, then perform 64 rounds of the standard Ch, Maj,
large-sigma and modular-2^32 additions per block and feed forward. All its
32-bit instructions are scalar AArch64 W-register operations. Constants are
the public standard SHA-256 IV and 64 K words. A rotate costs one ROR.
State-variable rotation between rounds is free by register renaming.

Per round: Sigma1 (3 ROR + 2 EOR = 5), Ch (1 BIC + 2 AND + 1 EOR = 4),
Sigma0 (3 ROR + 2 EOR = 5), Maj (3 AND + 2 EOR = 5), T1 (4 ADD),
T2 (1 ADD), new a and new e (2 ADD), K and W loads (2) — 28 instructions.
Per expanded schedule word: sigma0 (2 ROR + 1 LSR + 2 EOR = 5), sigma1 (5),
3 ADD, 4 loads, 1 store — 18 instructions, for 48 words. Feed-forward
8 ADD + loads/stores — 24. Per block: 64*28 + 48*18 + 24 = 2,680.
Seven blocks: 18,760. Serialization (104 words, load + byte-swap + store)
312; padding and setup 64. Total 19,136. Charge **H' = 20,000**.

A step reads one of 32 records, multiplies Y by its stored group element,
adds the two 4-word coefficients modulo r, and updates the state. The table
record is 448 bytes (384+32+32). The baseline's bound **A' = 4,096** is
retained unchanged for table access, coefficient adds, function control
and spills; with the tighter Fp bounds above the true cost is under 512,
so 4,096 carries more than 8x slack and needs no re-derivation.

Floyd performs three steps per iteration, one 48-word group equality test,
and increments/compares a three-word counter against N=3*ceil(sqrt(r))/2.
N is less than 2^128; three words avoid boundary conventions. The equality
test (48 CMP + control, <=128), counter work (<=16), state copies (<=80)
and loop control (<=16) total under 240 instructions. Charge **L' = 1,000**.
Therefore

    I' = 3*(G'+H'+A')+L' = 3*(40,000+20,000+4,096)+1,000 = 193,288

is a uniform worst-case bound per Floyd iteration, including memory access.
All loops have concrete integer bounds; the result is a cap, not an expected
work estimate. No cost depends on an unrepresented enormous address.

## One-target setup, termination and output

Use exactly one trial and at most N Floyd iterations. No retry is implicit.
Failure by the cap, rejection-sampling failure, degenerate collision, and
failed output verification terminate within the same cap and count as failure.

The following additional fixed **F=2^40** instructions cover both tracks:

* Validate bounded canonical input and convert 12 coefficients (including
  the one-time basis conversion, <2^17 instructions); check z^r=1.
  Invalid input is outside the stated success distribution. Valid JSON is
  bounded by the fixed 12x64-hex profile format; parse only that bounded form.
* Generate a 256-bit salt and 66 scalars in [0,r). Each scalar uses a
  254-bit candidate, rejecting >=r, at most 1024 times. Generate/read at most
  66*1024*256+256 fresh bits, charging even 64 instructions per bit, <2^31.
  Only one candidate and scalar at a time are retained. Physical entropy
  timing is not the instruction metric; these are the algorithm's ordinary
  fresh independent coins, not uncharged group/advice oracles.
* Build all 32 M_j=g^c_j*z^d_j, plus the initial g^a*z^b, charging two
  exponentiations and one multiplication per record, <=33*1017 GT products
  at G'=40,000: 33*1017*40,000 = 1,342,440,000 < 2^31. Charge all searched
  precomputation to this target, never to future inputs.
* Include subgroup testing, final g^x=z, coefficient inversion modulo r,
  modular division and the constant generation. Fewer than 35,000 GT
  products and 2^32 other instructions cover this part (<2^33).
* For pairing output and optional subgroup/curve checks use the affine G2
  procedure below, <2^31 instructions. GT scalar serialization is smaller.
* Up to 2^32 additional instructions cover bounded parsing/encoding, initial
  memory allocation, code/data loading, shutdown and all fixed overhead.

The above disjoint quantities sum to <2^36; F=2^40 leaves additional slack.
This allowance is derived from finite algorithms, not from Python timings.

## Pairing output and its explicit cost

After verifying recovered x with g^x=z, return P0 and Q=[x]Q0. Use binary
double-and-add in affine Fp2; at most 254 doubles and 254 additions. No pairing
computation is needed by the attack: its verified scalar and bilinearity
establish the output relation. The organizer's pairing check is excluded by
the protected rules.

The Fp2 here is the same Fp2 as the tower base (u^2=-1), so the tighter
bounds above apply: Fp2 multiplication costs 1,504 instructions
(3*448+5*32), Fp2 addition/subtraction 64, and Fp2 inversion
(a-bu)/(a^2+b^2) costs two squarings, one Fp inversion by p-2 (<=508
multiplications, 508*448 = 227,584 < 2^19), two scaling multiplications
and adds/negation, under 2^19 instructions. For nonzero a+bu the norm is
nonzero: p=3 mod 4, so -1 is not a square.

For nonexceptional affine point addition, compute
lambda=(y2-y1)/(x2-x1), X=lambda^2-x1-x2,
Y=lambda*(x1-X)-y1.
For doubling lambda=3*x1^2/(2*y1), X=lambda^2-2*x1,
Y=lambda*(x1-X)-y1.
Infinity, equal x/opposite y and y=0 are checked first and handled by the
usual infinity or doubling case; they do not call inverse at zero.
At most one Fp2 inverse (<2^19), 8 Fp2 multiplications (12,032) and 30
Fp2 add/sub/copy operations (1,920) plus <2^16 scalar control suffice for
any point operation: 524,288 + 12,032 + 1,920 + 65,536 = 603,776
< **2^20** instructions.

508 operations for [x]Q0 and 508 for the optional check [r]Q=infinity use
<2^30 instructions. Curve-equation evaluation, canonical conversion and
serialization, infinity and P0 handling fit in <2^30. The total is <2^31,
inside the <2^35 allowed in F above. An order-r multiple of Q0 is on the
required subgroup; the optional explicit checks add evidence without weakening
that argument. x=0 returns Q=infinity, which the protected profile permits
exactly for z=1.

## Peak storage

32 table records: 14,336 bytes. Three walk states and group scratch: <8192.
SHA scratch/encoding: <4096. Scalar/curve temporaries, random draw and loop
counters: <16,384. Profiles/constants: <8192. The scalar implementation,
SHA, arithmetic, parser and output code are finite loops; 128 KiB suffices
even at 32,768 scalar instructions of 4 bytes each. Reserve another 256 KiB
for stacks, code-loading scratch, allocation metadata and output.
Sum <2^20 bytes. Iteration counts do not grow live storage; no trace or
random tape is retained. This is the algorithm's total code/data/processor
storage on a scalar reference machine, not the excluded organizer or Python
test interpreter. Advice bytes are zero (log2(max(1,advice))=0).

## Direct validation of the C code

experiment.py validates arithmetic.c directly rather than restating its
constants. It compiles arithmetic.c together with a test harness
(cc/gcc/clang, -O2 -DCOUNT_OPS), runs it, and checks the printed vectors
against the independent schoolbook tower reference in rho.py: v*v^2 =
xi = 9+u (the vector the first revision got wrong), w*w = (v, 0),
64 aliased fp2_mul_xi calls against their unaliased twins, 32 random
Fp6 products and 32 random Fp12 products. The harness also reports the
operation counts for exactly one fp12_mul: 54 mont_mul, 171 add_mod,
109 sub_mod (280 Fp additions), matching the call-graph count derived
above. The Fp12 instruction total is computed from those measured
counts (54*448 + 280*32 + 512 = 33,664), and G' = 40,000 is charged
above it. When no C compiler exists in the experiment environment (the
shared verifier's slim container has none), experiment.py executes a
faithful Python transcription of arithmetic.c's exact call graph
instead and reports that honestly as simulated; the compiled validation
was executed during development and is recorded in
development-evidence.json, and both paths produce identical counts.

## Total and score

    N = 3*ceil(sqrt(r))/2
    R_total = N*193,288 + 2^40
    W = R_total/20,000
    submitted score = ceil(10^6*log2(W))/10^6 = 130.655989

Memory log2 bytes <=20. Preprocessing is <=2^36 instruction equivalents,
or 2^36/20,000 reference units, already included in F. No history, earlier
score, former rejected estimate, or existing table is used as a premise.
