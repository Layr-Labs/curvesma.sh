# Deterministic fixed-G1 inversion with explicit conservative accounting

This is an analytic classical cold attack for the exact production profile. The
claimed total work is at most 2^146 dense-Fp12 reference units, at C=20000
scalar AArch64 instruction equivalents per unit. Success is one on every valid
input. Peak storage is at most 2^142 bytes. Preprocessing, including the complete
baby table, is at most 2^146 units and is included in the total. There is no
searched nonuniform advice. Only the published profile constants are embedded.
The bound is intentionally loose. It is not a practical attack, a novelty claim,
or a claim to improve the best published heuristic number-field-sieve estimate.

## 1. Algebra and exact input

Use p, r, P0, Q0, the field polynomial and the pinned pairing from the protected
profile. In particular r is
0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001.
Thus 2^254 < r < 2^255. Let m=2^128. The supplied P is nonidentity in G1.
Set g=e(P,Q0), using the library argument order pairing(Q0,P), with final
exponentiation. The standard nondegenerate bilinear prime-order pairing implies
that g has order r: writing P=[b]P0 with b nonzero modulo r gives
 g=e(P0,Q0)^b, a generator. This uses the production pairing's defining
nondegeneracy property, not knowledge of b. We never recover b or replace P.
The target z is an arbitrary member of GT, including 1. There is a unique
0<=a<r with z=g^a. Once it is found, Q=[a]Q0 solves the required instance by
bilinearity. Q is returned as the canonical Fp2 point object, or null for infinity.

## 2. Finite deterministic algorithm

Construct all m pairs (g^j,j) for 0<=j<m, starting from (1,0) and multiplying
by g between consecutive records. Store all twelve canonical Fp coefficients
of each field element. Sort the table lexicographically by this complete key.
No fingerprint, hash assumption, collision heuristic or random-oracle model is
used. Canonical reduction modulo p is done after every field operation.

Compute h=g^(-m) by binary powering, with g^(-1)=g^(r-1). Put t=z.
For i=0,...,ceil(r/m)-1, binary-search the sorted table for t. If t=g^j,
set a=im+j. Accept this match only if a<r and g^a=z. Return [a]Q0 after
checking e(P,[a]Q0)=z. If no accepted match occurs, multiply t by h and
continue. The final unused update may be omitted. Binary exponentiation uses
ordinary integer exponents, never discrete logarithms or an oracle.

For the unique true a, Euclidean division gives a=im+j with 0<=j<m and
0<=i<ceil(r/m). At that iteration t=z g^(-im)=g^j. The baby powers are
distinct because m<r. Thus binary search finds j. The range guard and the
two verification checks pass, and output is the required Q. There is no failure
case on valid inputs. In particular a=0 is found with i=j=0 and returns infinity.
This proves success one, exceeding 0.39 for the specified independent uniform
P,z distribution without excluding any subset of inputs or requiring coins.

## 3. Explicit arithmetic model, with deliberately generous upper bounds

The attack can use the profile's polynomial representation
Fp[w]/(w^12-2w^6+2). It does not require a field isomorphism or smoothness.
All following work is scalar instruction-equivalent work. A 64-bit add,
subtract, comparison, Boolean operation, shift, multiplication yielding a low
or high word, scalar load/store, branch and address operation each costs one.
For loops we allow at least several operations per primitive for control.
No 381-bit operation is counted as one instruction.

For a conservative constructive base-field multiply, store each input in six
64-bit limbs. Compute the 12-limb product using the schoolbook 36 limb products,
using low/high multiplies and add-with-carry. Each partial product accumulation,
including loop/control and a carry ripple through all 12 positions, costs at
most 256 scalar instructions; all 36 cost below 2^14. Reduce the resulting
nonnegative integer of at most 768 bits by binary long division. Maintain a
7-limb remainder below 2p. For each input bit, shift the remainder by one,
insert the bit, compare with p, and conditionally subtract p. Fewer than
256 scalar instructions suffice per bit, including loading the input bit,
seven-limb shift, comparison, subtraction, selection, stores and loop control.
For all 768 bits this is below 2^18. Together with entry/exit and canonical
six-limb output, fewer than 2^20 scalar instructions suffice per Fp multiply.
The generous bound also covers squaring with the same algorithm. Addition or
subtraction modulo p uses at most 256 instructions: operate on six limbs with
carry/borrow, compare, subtract/add the modulus and select. Reduction by a
small integer multiple is implemented by repeated additions, not a free wide op.

A dense Fp12 multiply is implemented by the 144 schoolbook Fp coefficient
multiplications, and accumulation modulo p in the degree-22 polynomial. All
144 products cost at most 144*2^20 < 2^28 instructions. Their accumulations
cost fewer than 144*256 < 2^16. Reduce degrees from 22 down through 12 by
w^k=2w^(k-6)-2w^(k-12), a consequence of the exact profile polynomial.
There are 11 high degrees; each requires at most six modular additions or
subtractions (doubling and two updates, with scratch allowance). This is below
2^15 instructions. Initializing, loading, storing and loop control over these
fixed arrays is below 2^18 more. Hence a complete canonical dense Fp12 product
costs less than 2^30 instructions, including all arithmetic and memory work.
The margin between the displayed subtotal and 2^30 is more than a factor of
three. This bound is derived from elementary algorithms, not a timing claim or
an assertion that the fixed reference price equals this slow implementation.

Fp12 inversion where needed uses x^(p^12-2) by binary powering. Since p<2^381,
the exponent has fewer than 4572 bits. At most 9144 dense products suffice.
Group GT inverse is cheaper using r-1, but either bound is harmless. Public
integer exponent construction uses bounded limb arithmetic; p^12 has 72 limbs.
The target's coefficients are stored in 72 words (576 bytes). These are exact
full-production sizes throughout. No toy curve substitutes for the target.

## 4. Work for the search and table

Because r<2^255 and m=2^128, ceil(r/m)<=2^127. There are fewer than 2^128
baby multiplications and at most 2^127 giant updates. A successful exponent
check needs at most 510 additional dense products. The other fixed powers
and pairing operations are covered separately below. Thus search-related
Fp12 arithmetic, including this check, uses fewer than 2^129 products. By the
bound above this is less than 2^159 instructions.

Store each table record in 1024 bytes. It contains the complete 576-byte field
key, 16-byte baby index and padding/scratch. Sort by a stable bottom-up mergesort
using two buffers, retaining only these m baby records, not all giant points.
There are exactly 128 merge levels. Each merge level copies every record once
and performs at most m comparisons. A comparison of two full keys scans at most
72 words, with fewer than 1024 scalar instructions including address/control.
A complete 1024-byte record copy costs fewer than 1024 instructions including
load/store, address and control, using a scalar loop. Using 2^13 instructions
per record per merge level is a safe joint bound, including the merge-index
arithmetic and loop management. Sorting therefore costs less than
2^128 * 128 * 2^13 = 2^148 instructions. Allocating and zeroing both buffers
with scalar loops costs fewer than 2^140 instructions and is included below.

A binary search on m records takes at most 129 comparisons. Each uses the full
canonical field key, and costs at most 2^13 instructions jointly with table
address construction and updates. At most 2^127 searches therefore cost less
than 2^127*2^8*2^13=2^148 instructions. There is no collision-resolution cost
beyond this exact comparison. Sequential generation, counters, recording baby
indices, input parsing and other search bookkeeping are assigned an additional
2^148 instruction budget, which exceeds 2^20 instructions for each iteration
of the larger loop. This also covers table allocation/initialization.

Addresses and table indices here exceed a single 64-bit word. They are not
assumed to fit: an offset below 2^139 is represented by three 64-bit limbs.
Addition, subtraction, comparison, halving, and multiplying a record index by
1024 take fewer than 128 instructions with these fixed-width multiword values.
The pointer overhead is included in the 2^13 per-operation and 2^20 per-loop
allowances above. One can realize the abstract huge store as banks of at most
2^60 bytes, with bank identifier and within-bank address; the bank directory
has fewer than 2^80 entries and is included in the storage/work allowance.
This is a conventional analytic memory-access accounting model: each touched
64-bit word, its addressing and control are charged; no distance-dependent
physical latency is claimed. As with exponential BSGS itself, an actual VM
cannot store this table. No production execution is asserted or required by
this analytic challenge. If physical geographic routing rather than the stated
scalar load/store convention were imposed, that would require a different
cost analysis. This conditional interpretation is disclosed explicitly.

## 5. Pairing, conversion, output and all fixed costs

The pinned pairing can be evaluated directly over the profile's Fp12 field.
The fixed ate loop has fewer than 65 bits. A Miller iteration comprises fixed
curve addition/doubling and line formulas and optional affine inversions;
use even the slow affine formulas with inversion by binary powering. The
final exponent has fewer than 4572 bits and takes fewer than 9144 products.
A deliberately generous allowance of 2^24 dense-field multiplication-equivalent
calls covers both pairings, their curve formulas, inversions, input validation
and scalar multiplication. For specificity, an affine point operation has
at most ten extension-field adds, ten multiplies and one inversion, inversion
costing fewer than 9144 products; fewer than 2048 point/line operations suffice
for two 65-bit loops, two subgroup checks and 256-bit double-and-add output.
Their subtotal is below 2^25 products; use 2^27 to cover all such work instead.
Fp and Fp2 operations can be embedded in Fp12 at this conservative bound.
Together with final exponentiation, decoding, encoding, integer arithmetic,
and public-constant handling, a fixed allowance of 2^80 scalar instructions
is ample (2^27 * 2^30 = 2^57, leaving more than 2^22 multiplicative slack).
This includes construction and verification of Q, and an optional full input
validation pass. No cost of verifier service or challenge generation is charged
because that is organizer work, not an attack oracle.

Summing strictly bounded search work and the explicit allowances gives
 I < 2^159 + 3*2^148 + 2^80 < 2^160 scalar instructions.
No primitive is additionally charged as a reference operation. Therefore
 W = I/20000 < 2^160/20000 < 2^146,
because 20000>2^14. This establishes the submitted bound unchanged; it is
unnecessary to infer a tighter score. All preprocessing and failures are
included: deterministic success has no failed restart or amplification.
Preprocessing alone is less than the total and hence also bounded by 2^146.
Random sampling is absent. Constants such as r and Q0 are public profile data;
no table or external advice is supplied for free.

## 6. Peak storage and evidence scope

Two sorting buffers occupy exactly 2*m*1024 = 2^139 bytes. The directory,
implementation, public profile, pointers, counters, binary-search scratch,
field temporaries, pairing state, and output fit within another 2^139 bytes;
this is extremely loose, even allowing a directory entry of 64 bytes for each
2^60-byte bank. Thus total storage is below 2^140 bytes; we declare the looser
2^142-byte bound. The sorting scratch is reused thereafter. We do not retain
all giant steps or coins. Advice size is zero; the convention log2(max(1,0))=0
applies. No executable experiment is declared, and no computational evidence is
claimed. This is solely a finite algorithm, algebraic completeness argument,
and conservative analytical instruction bound. It requires enormous hypothetical
memory and work, so it establishes neither feasibility on this VM nor a practical
cryptographic break. Arithmetic correctness can be implemented by the elementary
algorithms specified above without relying on an external measurement.

## 7. Attribution and limitations

Baby-step/giant-step is the standard Shanks generic discrete logarithm method.
Its use here follows the permitted GT-logarithm reduction in the protected rules.
The construction and elementary conservative arithmetic bound were written for
this candidate. The repository's original historical fixed-Q/recover-G1 baseline
manifest was inspected, but its proof, kernels, constants and experiments were
not used. Recent public notes include ethvs's deterministic Shanks control at
162 and kaankolcu's much stronger heuristic STNFS estimate near 119.215; these
are comparisons only, not accepted premises or source code for this argument.
We make no claim that 146 beats those heuristic estimates, nor that a valid
review causes manual promotion. Independent critique and official configured
review are separate from this mathematical claim.
