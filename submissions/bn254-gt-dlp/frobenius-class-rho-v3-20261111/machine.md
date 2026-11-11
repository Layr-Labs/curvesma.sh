# Finite scalar implementation and work bounds

Conservative **analytical scalar AArch64 instruction-equivalent bounds**, not
timings. One processor, 64-bit scalar registers, ordinary addresses inside a
4 MiB allocation. Python (`rho.py`) is the executable mathematical reference,
not the priced implementation.

**Provenance.** The primitive bounds in the first three sections (Montgomery
word arithmetic, the polynomial-basis Fp12 product G, SHA-256 hashing H and
the table record step A) are taken verbatim in substance from the organizer's
reviewed baseline `submissions/bn254-gt-dlp/rho-baseline-v2-20261110`
(`machine.md`, `arithmetic.c`), with credit. They are restated here in full so
that this submission is self-contained; nothing below relies on that baseline
having been accepted. The canonicalisation K, stack S, control, setup F and
storage sections are new.

## Word operations and modular arithmetic (from the baseline)

Let B=2^64 and R=B^4. A residue is four little-endian words; for m=p or m=r,
m is odd and <2^254. Multiplication operands are in Montgomery form aR mod m.
Generating n0=-m^-1 mod B and R^2 mod m takes six Newton word iterations and
at most 512 modular doublings, charged in setup.

`arithmetic.c` specifies product and REDC with fixed bounded loops and
alias-safe scratch. A 64x64 product is MUL plus UMULH; 128-bit sums use
ADDS/ADC, comparisons/subtractions SUBS/SBC and CSET, selection uses word
masks. No 128-bit unit-cost arithmetic is assumed. Per Montgomery product:

| Work | Count | Max scalar instructions each | Bound |
|---|---:|---:|---:|
| Schoolbook product cell | 16 | 48 | 768 |
| REDC product cell | 16 | 48 | 768 |
| Carry propagation cell, fixed bounds | 28 | 16 | 448 |
| q=t[i]*n0 and loop setup, 8 outer passes | 8 | 32 | 256 |
| Initialise 9 words, subtract m, select/copy result | 1 | 512 | 512 |
| Calls, spills, remaining loop control | 1 | 512 | 512 |
| Total | | | 3264 < **MM=4096** |

A product cell loads two factors, t and carry (at most four address/load
pairs), uses two multiplies and at most four carry instructions, stores its low
word and carry and advances its index: below 48. A propagation cell is below 16;
there are 5+4+3+2 cells per propagation stage, twice. Each product-plus-two-
words is at most B^2-1; REDC output is <2m<2^255, so one conditional
subtraction makes it canonical. `rho.py:montgomery_word_reference` executes the
same word loops with assertions on every carry; `experiment.py` checks it on
544 edge/random full-width cases for both p and r.

Modular addition, subtraction, negation and four-word comparison each cost at
most **AS=256** instructions; a four-word copy at most 32 (8 per word). An
inverse modulo m is an exponentiation to m-2 with <=508 products and <=2^16
overhead, <2^22 instructions.

## Fp12 product and squaring (from the baseline)

Protected polynomial basis, canonical Montgomery coefficients, w^12=18w^6-82.
The degree-22 convolution uses 144 Montgomery products and 144 additions;
reducing k=22..12 uses <=400 add/sub/neg operations in total (x18 by five
doubles and two adds, x82 by seven doubles and three adds, two accumulations
and a negation per high coefficient). 144*4096+400*256 = 692,224; zeroing,
copying, indexing, calls and spills fit in a further 2^16 (155 product/
reduction cells and <200 copy cells at 128 overhead each). Total < **G=2^20**.
Squaring uses the same routine. An exponentiation by a <2^254 scalar uses <=508
products plus <=2^20 instructions of bit handling.

## Hashing and table step (from the baseline, index widened to 10 bits)

j is taken from SHA256(salt || encode(Y)), salt 32 bytes, encode(Y) the
injective 384-byte concatenation of the 12 canonical Montgomery residues in
big-endian 32-byte form (no conversion out of Montgomery form). Ordinary
scalar SHA-256, seven 64-byte blocks for the 416-byte message: at most 100
instructions per expanded-word/round pair, 7*64*100=44,800, plus 16,384 for
serialisation, padding, state setup, feed-forward and extraction: **H=2^16**.
Here j = (256*d[0]+d[1]) mod 1024 for digest bytes d[0],d[1]: a load, shift,
OR and AND inside H.

A table record is 448 bytes (384-byte element, two 32-byte scalars). Forming
j*448 plus a base pointer, reading/copying at most 112 words, two modular
additions of labels, control and spills fit in **A=4096**, unchanged from the
baseline (the address arithmetic does not depend on the table having 32 or
1024 records).

## Canonical class representative (new): K = 2^18

Frobenius phi(y)=y^p is Fp-linear on Fp12. In the protected basis it has a
fixed matrix; phi^k has matrix M_k (k=0..11), computed in setup and stored as
Montgomery constants (mont_mul(yR, mR) = (ym)R). `experiment.py` derives M_k
from the profile and verifies: every nonzero entry of every M_k lies in a
2x2 block on coefficient pair {i,i+6}; even k are diagonal (12 nonzeros), odd
k have 23 nonzeros; M_6=diag((-1)^i); M_1*M_11=I; and phi^k(g)=g^(p^k mod r)
on the generator for all k.

Order: elements are compared lexicographically on their canonical Montgomery
coefficients in priority (1,7,0,2,3,4,5,6,8,9,10,11). Each coefficient is
compared as a four-word integer.

Given Z (the product just formed), the representative is the minimum of the
12 conjugates phi^k(Z):

| Work | Count | Each | Bound |
|---|---:|---:|---:|
| Keys (c1,c7) of phi^k(Z), k=1..5, rows 1 and 7 (read only z1,z7) | 16 MM + 6 AS | | 67,072 |
| Keys for k=6..11: negate the key of k-6 (rows 1,7 are odd) | 12 AS | 256 | 3,072 |
| Running minimum over 12 keys, two coefficient compares each | 22 AS | 256 | 5,632 |
| Apply the winning M_k (<=23 products, <=11 two-term row sums) | 23 MM + 11 AS | | 97,024 |
| Scale labels a,b by lambda_k = p^k mod r (Montgomery constants mod r) | 2 MM | 4096 | 8,192 |
| Copies of 48 four-word values | 48 | 32 | 1,536 |
| Loop control, indexing, spills | 1 | 16,384 | 16,384 |
| Total (tally printed by experiment.py: 198,912) | | | < **K=262,144** |

If the minimum two-coefficient key is attained by two or more k, the run
takes the *tie path*: it applies all tied M_k fully and compares full keys
(<=12 applications and <=11 full compares, each full compare <=12*256). This
path is budgeted separately: at most TIE_MAX=1024 tie events per run (a 1025th
aborts as failure), each <2^22 instructions, charged to F. The identity is the
only GT element with equal conjugates (r does not divide p^k-1 for proper
divisors k of 12, verified via the order of p mod r); expected ties on the
walk are ~2^-121.7 (printed), so the cap is a cost guarantee, not a likely
event. With the tie path the representative is a genuine class function.

## Stack cycle detector (new): S = 2^13

The stack holds entries (Y, a, b): 384+64 bytes, padded to 512. On each new
element: while the top is larger than Y in the order, pop; if the top equals
Y, test B-b mod r (detection); otherwise push. A full order compare is at most
12 coefficient compares, 12*256=3,072. Every compare either pops (at most once
per pushed entry) or ends the loop (at most once per step), so over a run the
number of compares over a run of s steps is at most 2s+1: charge two per
step (6,144) and the extra one to F. This is an amortised, not a per-step,
bound: a single step may pop many entries, but the run total is bounded. Push
copies 56 words (<=448), the degeneracy test is one four-word compare (256),
depth check and pointer update 256. Tally 7,104 < **S=8192**. The depth cap
is D_MAX=4096 entries (2 MiB); exceeding it aborts as failure.

## Per-step bound

Every produced element is one step. A normal step: hash (H), table record
and label additions (A), one Fp12 product (G), canonicalisation including
label scaling (K), stack processing (S), and control (a three-word step
counter compared with E, the escape flag and loop branches, <=**CTRL=4096**).
An escape step replaces H+A by label doubling (two additions) and the product
by a squaring (G): it costs less. Hence any run of s <= E steps costs at
most s*I instructions beyond F, where

    I = G + H + A + K + S + CTRL
      = 1,048,576 + 65,536 + 4,096 + 262,144 + 8,192 + 4,096 = 1,392,640.

All terms except S are worst-case per step; S is the amortised share above
(pops <= pushes <= s+1), so the statement is about run totals.

The step counter has three words; E<2^126.

## One-target setup, termination and output: F = 2^42

All of the following are charged once, deterministically bounded
(`experiment.py` itemises them with these exact formulas and asserts the sum
is below 2^42; the printed tally is about 2^40.08):

* Bounded input parse/validation and z^r=1 check: <=1017 GT products plus
  2^32 for parsing/conversion (<2^33).
* Randomness: 256-bit salt and 2*1025 scalars by 254-bit rejection, at most
  1024 draws each: <2^30 fresh bits at 64 instructions per bit, <2^36.
* Table of 1024 records plus the starting element: 1025*(two exponentiations
  and one product) <= 1025*1017 GT products plus 2050*2^20 overhead, <2^40.1.
* Frobenius constants: phi(w)=w^p by one exponentiation, 11 12x12 matrix
  products modulo p (1728 MM each), lambda_k by 11 products mod r, Montgomery
  conversion of 12*23 entries: <2^30.
* Tie budget: 1024 events, each <= 12 full applications (23 MM + 11 AS)
  plus 11 full compares (12 AS each) plus 2^16 control: < 2^31.
* Start element: one canonicalisation (K), its push and the extra stack
  compare (S).
* Final work: inverse of B-b mod r, one product, verification g^x=z (<=1017
  GT products), scalar output; for the pairing track, P0 and [x]Q0 by the
  affine G2 procedure below (<2^35).
* Stack/table allocation and initialisation, code loading, fixed overhead: <2^32.

The sum is below 2^41; F=2^42 leaves slack. Failure by cap, stack cap, tie
budget, sampling failure or failed verification terminates inside the same
bounds.

## Pairing output and its explicit cost (from the baseline)

After verifying x with g^x=z, return P0 and Q=[x]Q0 by binary double-and-add
in affine Fp2 coordinates: at most 254 doublings and 254 additions. No pairing
is computed by the attack; the verified scalar and bilinearity give
e(P0,[x]Q0)=g^x=z.

For u^2=-1, (a+bu)(c+du)=(ac-bd)+(ad+bc)u takes four Fp products and two
additions, <2^15 instructions with loads/control. Fp2 add/sub: two Fp
operations plus copies, <=1024. Inversion (a-bu)/(a^2+b^2): two squarings,
one Fp inverse by p-2, two scaling products and adds/negation, <2^23; the
norm is nonzero for nonzero input since p=3 mod 4.

Nonexceptional addition: lambda=(y2-y1)/(x2-x1), X=lambda^2-x1-x2,
Y=lambda(x1-X)-y1; doubling: lambda=3x1^2/(2y1), X=lambda^2-2x1,
Y=lambda(x1-X)-y1. Infinity, equal x/opposite y and y=0 are tested first and
handled by the infinity or doubling case, never inverting zero. Any point
operation: at most one Fp2 inverse, 8 Fp2 products, 30 Fp2 add/sub/copy and
<2^16 control, < **2^24**.

508 operations for [x]Q0 and 508 for the optional [r]Q=infinity check are
<2^34; curve-equation checks, canonical serialisation, infinity and P0
handling <2^30; total <2^35, charged in F. x=0 (only for z=1) returns the
permitted infinity output.

## Storage

Table 1024*448 = 458,752 bytes; stack 4097*512 = 2,097,664; twelve Frobenius
matrices (<=23 Montgomery constants each) 8,832; walk state, scratch, SHA,
scalar/curve temporaries, profile/constants, code and stacks: <2^20 (baseline
allocation). Total < 2^22 bytes (tallied and asserted in experiment.py). No trace or random tape is retained. Advice
is zero.

## Total

    n  = (r-1)/12 + 1                                 (classes)
    E  = isqrt(floor((9/5)^2 * n * 1024/1023)) + 1    (step cap, ~2^125.854566)
    R_total = E*1,392,640 + 2^42
    W  = R_total/20,000,  log2 W <= 131.976245 (rounded up)
    claimed score = 131.976245 + 0.020000 margin = 131.996245

`experiment.py` recomputes every integer above from the protected profile.
