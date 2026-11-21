# Finite scalar implementation and work bounds

These are conservative **analytical scalar AArch64 instruction-equivalent
bounds**, not timings and not prices inferred from the previous slow C kernel.
The implementation below uses only 64-bit scalar registers and ordinary
addresses within a 1 MiB allocation. Python is an executable mathematical
reference, not the priced implementation. There is one processor.

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

Per modular multiplication, the following deliberately loose counts suffice:

| Work | Count | Maximum scalar instructions each | Bound |
|---|---:|---:|---:|
| Schoolbook product cell | 16 | 48 | 768 |
| REDC product cell | 16 | 48 | 768 |
| Carry propagation cell, fixed bounds | 28 | 16 | 448 |
| q=t[i]*n0, loop setup for eight outer passes | 8 | 32 | 256 |
| Initialize 9 words, subtract m, select/copy 4-word result | 1 | 512 | 512 |
| Calls, spills, remaining public loop control | 1 | 512 | 512 |
| Total | | | 3264 < **4096** |

A product cell loads two factors, t and carry (at most four address/load pairs),
uses two multiplies and at most four carry instructions, stores its low word
and carry, and advances/tests its index. Even allowing separate address
instructions and stack-resident temporaries, this is below 48. A propagation
cell loads two words, adds with carry, stores, and controls its fixed loop,
below 16. There are 5+4+3+2 cells in each propagation stage, twice.

Correctness/overflow: each product-plus-two-words is at most B^2-1. Full
product a*b<R^2 is represented in 9 words. REDC cancels the next low word at
each pass; its final U=(a*b+q*m)/R is <2m since a,b<m<R and q<R. Thus U<2^255,
the ninth word is zero and one subtraction of m makes the answer canonical.
Intermediate carry propagation is retained through word 8. Four-limb
addition cannot overflow 256 bits because 2m<2^255.

Modular addition, subtraction, negation and comparison each cost at most
**256** instructions. They use at most four carry/borrow cells, a four-word
conditional correction/selection and fixed loads/stores. `add_mod` is supplied;
subtraction subtracts b with borrow and conditionally adds m, and negation
subtracts from m and masks the zero case. A four-word copy costs <=32.
These bounds also cover scalar coefficients modulo r. Inverse of a nonzero
residue is exponentiation to m-2 with <=508 modular multiplications and
<=2^16 scalar overhead, hence <2^22 instructions. Representation conversions
each cost one Montgomery multiplication; explicit setup/cleanup allows for them.

## Fp12 multiplication

Use the protected polynomial basis and canonical Montgomery coefficients.
Form the degree-22 convolution with 144 modular multiplications and 144
modular additions. Reduce k=22,...,12 using w^12=18w^6-82.
Multiplication by 18 uses <=7 modular additions (five doubles/two adds);
by 82 uses <=10 (seven doubles/three adds); two accumulator additions and
one negation suffice per high coefficient. Thus <=400 modular
add/sub/neg operations cover the convolution and reduction.

144*4096 + 400*256 = 692224.
Zeroing/copying the 23-coefficient scratch and 12-coefficient output and all
the small loops, indexing, calls and spills fit within a further 2^16
instructions: there are 144 convolution cells, 11 reduction cells and under
200 other word-copy/zero cells; 128 overhead instructions per cell suffices.
The total is < **G=2^20**. Squaring may use this same multiplication.
There is no double counting of a normalization price. Polynomial scratch
occupies <2 KiB and is reused, never recursively allocated.

Fp12 exponentiation by a <2^254 scalar uses <=508 such multiplications,
plus <=2^20 instructions for bits, copies and loop overhead.

## Hashing and walk state

Hash SHA-256(salt || encode(Y)), taking the low five bits of digest byte zero
as the table index. Salt is 32 bytes. encode(Y) is the 384-byte concatenation
of its twelve **canonical Montgomery residues**, each in fixed big-endian
32-byte form. This internal encoding is injective, and explicitly differs
from the external standard-coefficient input encoding. No conversion out of
Montgomery representation is needed on each step. Scalar software SHA-256
uses seven 64-byte padded blocks for this 416-byte message.

For clarity, this is ordinary SHA-256: initialize the standard eight IV
words, expand the 16 input words to 64 with the standard small-sigma rotates/
shifts and additions, then perform 64 rounds of the standard Ch, Maj,
large-sigma and modular-2^32 additions per block and feed forward. All its
32-bit instructions are scalar AArch64 W-register operations. Constants are
the public standard SHA-256 IV and 64 K words. A rotate costs one ROR.

At most 100 scalar instructions per expanded-word/round pair allow separate
loads, stores and addressing for eight working words, K, the ring-buffer W,
all sigma/Ch/Maj operations and control (each round has fewer than 40 ALU
operations). Charging all 64 pairs even for the first 16 unexpanded words
gives 7*64*100=44800. A further 16384 covers 416-byte serialization, padding,
block decoding, state setup, feed-forward and extraction. Thus **H=2^16**
includes hashing and encoding. A single input byte is never an uncharged oracle
query. This is a finite deterministic implementation; random-mapping behavior
of the walk is a disclosed heuristic, not a free ideal hash implementation.

A step reads one of 32 records, multiplies Y by its stored group element,
adds the two 4-word coefficients modulo r, and updates the state. The table
record is 448 bytes (384+32+32). Multiplying the five-bit index by 448 and
adding a normal pointer, reading/copying at most 112 words, two modular
adds, function control and spills fit in **A=4096** beyond G+H.

Floyd performs three steps per iteration, one 48-word group equality test,
and increments/compares a three-word counter against
N=162741432569968250906303761388095885571 = ceil(1.1*sqrt(r)).
N is less than 2^128; three words avoid boundary conventions. An additional
**L=65536** instructions per iteration more than covers these comparisons,
copies and control. Therefore

    I = 3*(G+H+A)+L = 3,420,160

is a uniform worst-case bound per Floyd iteration, including memory access.
All loops have concrete integer bounds; the result is a cap, not an expected
work estimate. No cost depends on an unrepresented enormous address.

## One-target setup, termination and output

Use exactly one trial and at most N Floyd iterations. No retry is implicit.
Failure by the cap, rejection-sampling failure, degenerate collision, and
failed output verification terminate within the same cap and count as failure.

The following additional fixed **F=2^40** instructions cover both tracks:

* Validate bounded canonical input and convert 12 coefficients; check z^r=1.
  Invalid input is outside the stated success distribution. Valid JSON is
  bounded by the fixed 12x64-hex profile format; parse only that bounded form.
* Generate a 256-bit salt and 66 scalars in [0,r). Each scalar uses a
  254-bit candidate, rejecting >=r, at most 1024 times. Generate/read at most
  66*1024*256+256 fresh bits, charging even 64 instructions per bit, <2^31.
  Only one candidate and scalar at a time are retained. Physical entropy
  timing is not the instruction metric; these are the algorithm's ordinary
  fresh independent coins, not uncharged group/advice oracles.
* Build all 32 M_j=g^c_j*z^d_j, plus the initial g^a*z^b, charging two
  exponentiations and one multiplication per record, <=33*1017 GT products.
  Charge all searched precomputation to this target, never to future inputs.
* Include subgroup testing, final g^x=z, coefficient inversion modulo r,
  modular division and the constant generation. Fewer than 35,000 GT
  products and 2^32 other instructions cover this part (<2^36).
* For pairing output and optional subgroup/curve checks use the affine G2
  procedure below, <2^35 instructions. GT scalar serialization is smaller.
* Up to 2^32 additional instructions cover bounded parsing/encoding, initial
  memory allocation, code/data loading, shutdown and all fixed overhead.

The above disjoint quantities sum to <2^37; F=2^40 leaves additional slack.
This allowance is derived from finite algorithms, not from Python timings.

## Pairing output and its explicit cost

After verifying recovered x with g^x=z, return P0 and Q=[x]Q0. Use binary
double-and-add in affine Fp2; at most 254 doubles and 254 additions. No pairing
computation is needed by the attack: its verified scalar and bilinearity
establish the output relation. The organizer's pairing check is excluded by
the protected rules.

For u^2=-1, multiplication (a+bu)(c+du)=(ac-bd)+(ad+bc)u takes four Fp
multiplications and two additions, <2^15 instructions including loads/control.
Fp2 addition/subtraction takes two Fp operations plus copies, <=1024.
Inversion is (a-bu)/(a^2+b^2), using two squarings, one Fp inversion by
p-2, two scaling multiplications and adds/negation, <2^23 instructions.
For nonzero a+bu the norm is nonzero: p=3 mod 4, so -1 is not a square.

For nonexceptional affine point addition, compute
lambda=(y2-y1)/(x2-x1), X=lambda^2-x1-x2,
Y=lambda*(x1-X)-y1.
For doubling lambda=3*x1^2/(2*y1), X=lambda^2-2*x1,
Y=lambda*(x1-X)-y1.
Infinity, equal x/opposite y and y=0 are checked first and handled by the
usual infinity or doubling case; they do not call inverse at zero.
At most one Fp2 inverse, 8 Fp2 multiplications and 30 Fp2 add/sub/copy
operations plus <2^16 scalar control suffice for any point operation:
< **2^24** instructions.

508 operations for [x]Q0 and 508 for the optional check [r]Q=infinity use
<2^34 instructions. Curve-equation evaluation, canonical conversion and
serialization, infinity and P0 handling fit in <2^30. The total is <2^35,
as included in F. An order-r multiple of Q0 is on the required subgroup;
the optional explicit checks add evidence without weakening that argument.
x=0 returns Q=infinity, which the protected profile permits exactly for z=1.

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

## Total and score

    N = 162741432569968250906303761388095885571 = ceil(1.1*sqrt(r))
    R_total = N*3,420,160 + 2^40
    W = R_total/20,000
    submitted score = ceil(10^6*log2(W))/10^6

Memory log2 bytes <=20. Preprocessing is <=2^36 instruction equivalents,
or 2^36/20,000 reference units, already included in F. No history, earlier
score, former rejected estimate, or existing table is used as a premise.
