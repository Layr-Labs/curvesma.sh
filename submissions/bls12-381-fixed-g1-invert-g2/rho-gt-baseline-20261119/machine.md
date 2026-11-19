# Finite scalar implementation and work bounds

These are conservative **analytical scalar AArch64 instruction-equivalent
bounds**, not timings. The implementation uses only 64-bit scalar registers
and ordinary addresses within a 1 MiB allocation. Python is an executable
mathematical reference, not the priced implementation. One processor.

## Word operations and modular arithmetic

Let B=2^64. The field prime p < 2^381 uses NL=6 limbs (R_p=B^6); the
subgroup order r < 2^255 uses NL=4 limbs (R_r=B^4). Both moduli are odd with
2m < B^NL. Multiplication operands are in Montgomery form. Public constants
may be loaded; generating n0=-m^-1 mod B and R^2 mod m instead takes six
Newton word iterations and at most 768 modular doublings, included in setup.

`arithmetic.c` specifies add and Montgomery multiply/REDC for both limb
counts with fixed bounded loops and alias-safe scratch. A 64x64 product is
MUL plus UMULH (two instructions); 128-bit sums use ADDS/ADC; comparisons and
subtractions use SUBS/SBC and CSET; selection uses word masks. No 128-bit
unit-cost arithmetic is assumed. Indexes are below 13 and can be unrolled.
Loads, stores, addressing, loops and calls are included below.

Per **six-limb** (mod p) Montgomery multiplication, deliberately loose:

| Work | Count | Max instructions each | Bound |
|---|---:|---:|---:|
| Schoolbook product cell | 36 | 48 | 1728 |
| REDC product cell | 36 | 48 | 1728 |
| Carry propagation cell (7+6+5+4+3+2 per stage, two stages) | 54 | 16 | 864 |
| q=t[i]*n0 and loop setup, twelve outer passes | 12 | 32 | 384 |
| Initialize 13 words, subtract m, select/copy result | 1 | 512 | 512 |
| Calls, spills, remaining loop control | 1 | 512 | 512 |
| Total | | | 5728 < **8192** |

The **four-limb** (mod r) instance has 16+16 product cells, 28 propagation
cells, 8 passes: 3264 < **4096**.

A product cell loads two factors, t and carry, uses two multiplies and at most
four carry instructions, stores its low word and carry, and advances its index:
below 48 even with separate address instructions and stack temporaries. A
propagation cell (load, add with carry, store, loop control) is below 16.

Correctness/overflow: each product plus two words is at most B^2-1. The full
product is held in 2NL+1 words. REDC cancels one low word per pass; the result
U=(a*b+q*m)/R_m < 2m < B^NL since a,b<m<R_m, q<R_m. Thus the top word is zero
and one conditional subtraction of m makes the result canonical. Addition
cannot overflow NL words because 2m < B^NL.

Modular addition, subtraction, negation and comparison each cost at most
**512** instructions (at most 6 carry/borrow cells, a 6-word conditional
correction/selection, fixed loads/stores). A 6-word copy costs <= 48.
Inverse of a nonzero residue mod p is exponentiation to p-2: <= 762 modular
multiplications plus <= 2^16 overhead, < 2^23 instructions. Inverse mod r:
<= 510 four-limb multiplications, < 2^22. Representation conversions cost one
Montgomery multiplication each and are covered by setup/cleanup allowances.

## Fp12 multiplication (protected basis)

Fp12 = Fp[w]/(w^12 - 2w^6 + 2), coefficients canonical Montgomery residues.
Form the degree-22 convolution with 144 modular multiplications and 144
modular additions. Reduce k=22,...,12 using w^12 = 2w^6 - 2: t=2*c_k (one
modular addition), c_{k-6} += t, c_{k-12} -= t: 3 add/sub per high
coefficient, 33 in total. So <= 400 modular add/sub/neg operations cover the
convolution and reduction:

    144*8192 + 400*512 = 1,384,448.

Zeroing/copying the 23-coefficient scratch and 12-coefficient output, loops,
indexing, calls and spills fit within a further 2^16 instructions (144
convolution cells, 11 reduction cells, under 200 copy/zero cells, at most 128
overhead instructions per cell). Total < **G = 2^21 = 2,097,152**. Squaring
uses the same routine. Scratch is < 2 KiB, reused, not recursive. (The pinned
blst tower multiplication measured 22,996 instructions; our deliberately
simple bound is ~91 times larger and is not a normalization substitute.)

Fp12 exponentiation by a scalar < 2^k uses <= 2k such multiplications plus
<= 2^20 instructions of bit extraction, copies and control.

## Hashing and walk state

Hash SHA-256(salt || encode(Y)) and take the low five bits of digest byte 0 as
the table index. Salt is 32 bytes; encode(Y) is the 576-byte concatenation of
the twelve **canonical Montgomery residues**, each in fixed big-endian 48-byte
form. This internal encoding is injective (canonical residues) and needs no
conversion out of Montgomery form per step. The 608-byte message pads to ten
64-byte SHA-256 blocks.

Ordinary scalar SHA-256: standard IV, message expansion to 64 words with the
standard small-sigma rotates/shifts and additions, 64 rounds of Ch, Maj,
large-sigma and mod-2^32 additions per block, feed-forward; all 32-bit
W-register operations; ROR is one instruction. At most 100 instructions per
expanded-word/round pair (fewer than 40 ALU operations, plus loads, stores and
addressing for eight working words, K and the W ring buffer, and control).
Charging all 64 pairs per block: 10*64*100 = 64,000. A further 16,384 covers
the 608-byte serialization (76 byte-reversed word stores), padding, block
decoding, state setup, feed-forward and extraction. Thus **H = 2^17** covers
hashing and encoding. No byte of the hash is an uncharged oracle query; the
random-mapping behaviour of the walk is a disclosed heuristic.

A step reads one of 32 table records, multiplies Y by its stored GT element
(one G), adds two four-word coefficients mod r, and updates the state. A record
is 640 bytes (576 + 32 + 32). Multiplying the five-bit index by 640 and adding
a base pointer, reading/copying at most 160 words, two modular adds mod r,
function control and spills fit in **A = 2^13 = 8192** beyond G+H.

Floyd performs three steps per iteration, one 72-word group equality test, and
increments/compares a three-word counter against N=2*ceil(sqrt(r)) < 2^129.
An extra **L = 65,536** per iteration more than covers comparisons, copies and
control. Therefore

    I = 3*(G+H+A) + L = 3*(2,097,152+131,072+8,192) + 65,536 = 6,774,784

is a uniform worst-case bound per Floyd iteration, including memory access.
All loops have concrete integer bounds; this is a cap, not an expectation.

## One-target setup, pairing, termination and output

Exactly one trial and at most N Floyd iterations. Failure by the cap,
rejection-sampling failure, degenerate collision and failed verification all
terminate within the same cap and count as failure. The following fixed
**F = 2^40** instructions cover everything outside the Floyd loop:

* Validate bounded canonical input (fixed 96-hex-character coefficients),
  check P on curve and z != 0, convert to Montgomery form; z^r = 1 check
  (<= 510 Fp12 products); optional [r]P = infinity check (affine G1, below).
* **Pairing g = e(P,Q0).** Evaluate the protected function
  py_ecc.optimized_bls12_381.pairing(Q0,P) (Miller loop over the protected
  64-bit ate loop count d201000000010000 = |x|, then exponentiation by the
  protected (p^12-1)/r). We implement its value with the dense Fp12 routine
  above, following the library's projective-in-Fp12 algorithm: twist Q0 and
  cast P into Fp12 (<= 64 Fp12-equivalent operations); then at most 64 loop
  iterations, each with one tangent line function, the two numerator/
  denominator square-and-multiply updates, one projective point doubling and,
  on the 5 set bits after the leading bit, one chord line function, two
  multiplications and one projective addition. Each line function, doubling
  or addition uses fewer than 32 Fp12 multiplications/subtractions/scalar
  multiplications, so an iteration uses fewer than 200 Fp12 operations, each
  bounded by one product G: < 64*200 = 12,800 products. Dividing f_num/f_den
  uses a Fermat inverse f_den^(p^12-2) (< 2^4572 exponent): <= 9,144
  products. The final exponent (p^12-1)/r has fewer than 4,320 bits:
  <= 8,640 products. Total < 31,000 < 2^15 products, < 2^36 instructions.
  The value, not the internal branch path, is what is used; it equals the
  protected function value on this exact P because it is the same
  deterministic computation.
* Generate a 256-bit salt and 66 scalars in [0,r): 255-bit candidates,
  rejecting values >= r, at most 1024 draws each; at most 66*1024*255+256
  fresh bits, even at 64 instructions per bit < 2^31. These are the
  algorithm's ordinary fresh coins, not group/advice oracles.
* Build the 32 M_j = g^c_j * z^d_j and the initial g^a0*z^b0: two
  exponentiations and one multiplication per record, <= 33*1021 products.
* Final g^x = z (<= 510 products), coefficient subtraction and inversion
  mod r (< 2^22), constant generation (< 2^24).
  All GT-product items above total < 80,000 < 2^17 products (< 2^38
  instructions) plus < 2^33 other instructions.
* Output Q=[x]Q0 and optional [r]Q=infinity check in affine G2 (below),
  < 2^36. Parsing, encoding, allocation, code loading, shutdown and all other
  fixed overhead < 2^32.

These disjoint quantities sum to < 2^39; F = 2^40 leaves slack. The allowance
is derived from finite algorithms, not from Python timings. It is charged in
full to this single target.

## G2 output and its explicit cost

After verifying g^x = z, return Q = [x]Q0 by binary double-and-add in affine
Fp2 coordinates: at most 255 doublings and 255 additions. No further pairing
is needed: the verified scalar and bilinearity establish e(P,Q)=z (the
organizer's acceptance pairing is excluded work).

With u^2 = -1, (a+bu)(c+du) = (ac-bd)+(ad+bc)u: four Fp multiplications and
two add/sub, < 2^16 instructions with loads/control. Fp2 add/sub < 2^11.
Inverse (a-bu)/(a^2+b^2): two squarings, one Fp inverse, two multiplications,
negation, < 2^24. For nonzero a+bu the norm is nonzero since p = 3 mod 4.
Affine addition lambda=(y2-y1)/(x2-x1), X=lambda^2-x1-x2,
Y=lambda(x1-X)-y1; doubling lambda=3x1^2/(2y1), X=lambda^2-2x1,
Y=lambda(x1-X)-y1. Infinity, equal x / opposite y and y=0 are tested first
and handled as infinity/doubling, never inverting zero. One Fp2 inverse, at
most 8 Fp2 multiplications, 30 add/sub/copies and < 2^16 control: < 2^25 per
point operation. 510 operations for [x]Q0 plus 510 for the optional
[r]Q=infinity check: < 2^35. The same affine formulas over Fp (cheaper) serve
the optional [r]P check on G1 (< 2^34). x=0 returns infinity, valid exactly
when z=1. Serialization and curve-equation checks < 2^30.

## Peak storage

32 table records: 20,480 bytes. Three walk states and Fp12 scratch: < 8192.
SHA scratch/encoding: < 4096. Pairing state (a few Fp12 values and
Fp12-projective points): < 16,384. Scalar/curve temporaries, random draw and
counters: < 16,384. Profile/constants (including the 4,314-bit final exponent
and p^12-2 as bit strings): < 8192. Code for these finite loops: < 128 KiB.
Another 256 KiB for stacks, loading scratch, allocation metadata and output.
Sum < 2^20 bytes. Iterations do not grow live storage; no trace or random tape
is retained. Advice bytes are zero.

## Total and score

    N = 2*ceil(sqrt(r))            (r is not a square; N = 2*(isqrt(r)+1))
    R_total = N*6,774,784 + 2^40
    W = R_total/20,000
    submitted score = ceil(10^6*log2(W))/10^6 = 136.832576

Memory log2 bytes <= 20. Preprocessing < 2^40 instruction equivalents, all
included in F. No history, earlier score or external table is a premise.
