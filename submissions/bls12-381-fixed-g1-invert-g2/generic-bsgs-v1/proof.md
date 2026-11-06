# Deterministic cold fixed-P inversion, conservative accounting

This is an analytic generic control, not a practical production run or a new
cryptanalytic technique. No heuristic or accepted prior submission is a premise.
The exact target is the protected production profile, including its pairing,
field modulus p, prime subgroup order r, and generators. C remains 20000.

## Algorithm and correctness

Let m=2^128. Compute g=e(P,Q0) using the specified reduced pairing. Nondegeneracy
and the prime order of the groups imply that g has order r for every allowed
nonidentity P. If z=1 return the G2 identity. Otherwise construct the table
(g^j,j), 0<=j<m, by successive multiplication, starting at the identity. Use the
canonical twelve-coefficient representation, compare lexicographically and sort
by bottom-up mergesort. Compute h=g^(-m) by binary powering and field inversion.
For i=0,...,ceil(r/m)-1, binary-search z*h^i in the table. On finding j, put
 a=(i*m+j) mod r, construct Q=[a]Q0 by binary double-and-add, and return Q after
checking the specified pairing equality e(P,Q)=z. Update the search value by
multiplication by h between iterations. Failure of the final verification is an
internal algorithm error, not an allowed successful output.

Every a in [0,r) has a unique quotient i=floor(a/m) and remainder j=a-i*m.
The loop includes that quotient. Thus z*h^i=g^j at that iteration. Conversely
any equality found implies z=g^(im+j), giving the required Q for the original
supplied P by bilinearity. The scalar is an internal intermediate; the output
is the actual canonically encoded G2 point. This establishes success 1, including
identity targets, over the full production distribution. P is not replaced by
P0 and its discrete logarithm is never needed. No powers input is used.

## Explicit conservative arithmetic implementation

The following deliberately slow implementation gives instruction upper bounds,
not timings or a claim to reproduce the organizer's dense multiplication cost.
Store each Fp coefficient in 381 one-byte bits. Implement addition, subtraction,
comparison and bit shifts with fixed-length loops over bits. A full-adder bit
step (loads, stores, boolean operations, address increments, carry and control)
uses at most 64 scalar instructions. An unsigned multiplication computes a
762-bit product by 381 shifted conditional additions of 762 bits. Reduce by
binary long division: for each of at most 762 positions compare the 763-bit
remainder with p and conditionally subtract. Each position uses fewer than
4*763*64 instructions. Together these stages take less than 2^29 scalar
instructions, including initialization and conversion. All branches can be
implemented literally; no constant-time requirement is imposed. Shorter
operations fit the same bound.

Represent Fp12 in the protected polynomial basis modulo w^12-2*w^6+2.
Convolve two degree-11 polynomials using 144 Fp products and at most 144
coefficient additions. Reduce high monomials from degree 22 downward using
w^12=2*w^6-2. At most 11 high monomials each cause two coefficient updates,
with doubling by addition. Fewer than 512 modular coefficient operations suffice.
The modular-coefficient bound above therefore gives less than 2^38 instructions
per complete Fp12 multiplication, including canonicalization and loop control.
An Fp12 inverse is binary powering to p^12-2, with fewer than 2^14 such
multiplications (p^12 has fewer than 4572 bits). Fp and Fp2 inverses use the
same construction with smaller exponents. Multiplication is implemented rather
than charged again as a reference operation: all instructions are divided by C.

A stored GT value uses 4572 bytes in this representation. Allow an 8192-byte
record including the exponent, alignment and tags. Comparing, copying or
indexing a record uses less than 2^22 instructions. Comparisons scan all bits
if needed. Table indices are two 64-bit words; the table has fewer than 2^142
bytes and addresses use three words. Addition, comparison and shifts of these
fixed widths are explicitly bounded by the same allowance. The memory model is
the challenge's scalar-instruction-equivalent accounting with total storage
reported separately, not an assertion that a literal current machine can allocate
the full table. One may partition storage into addressed pages; routing a request
through at most 142 binary decisions and returning a record remains below the
2^22 allowance. All data movement, addressing and control count. No table exists
before the cold attack, and its construction and sorting are charged.

## Setup, verification and output upper bound

Use the profile's rational Miller line formulas with Fp12 arithmetic and final
exponentiation by the exact integer (p^12-1)/r. There are fewer than 128 Miller
steps, each with a fixed number less than 1024 of field operations, even if every
line division is implemented with a full inverse. Including such inversions,
curve arithmetic, conversion, nonzero checks and final exponentiation costs
less than 2^32 Fp12 multiplications. This is intentionally excessive. Handle
exceptional projective line cases by the ordinary rational-function cases; the
specified pairing is well-defined on the prime-order inputs. Binary scalar
multiplication for the output uses fewer than 512 point additions/doublings;
using affine formulas with exceptional identity/inverse cases and powering for
inversion fits the same allowance. Two pairings, output scalar multiplication,
input validation (including subgroup powering), inversion to obtain h and all
fixed setup together are bounded by 2^74 instructions. The bound allows more
than four times 2^32*2^38 and includes canonical JSON input/output, with only
fixed-size production objects. This is polynomial setup, not an uncharged oracle.

## Complete work bound

The order satisfies r<2^255, so ceil(r/m)<=2^127. Table construction has fewer
than 2^128 Fp12 multiplications. Giant steps use fewer than 2^127. Both together
are less than 2^129 multiplications, hence fewer than 2^167 instructions using
the conservative 2^38 price. Table sorting takes 128 passes over m records.
Each pass performs at most m comparisons and at most 2m copies. Even charging
2^22 to each operation, with index/control allowances, sorting is below 2^160
instructions. Binary searches have at most 129 comparisons each, costing below
2^158 instructions in all; giant-loop scalar/control work is included by that
allowance. Setup and verification add less than 2^74. Other allocation,
initialization, canonical encoding, table reads/writes and cleanup are covered
by an extra 2^160 instructions, including freeing the table if required.
Therefore total R is strictly less than 2^168 scalar instruction equivalents.
Since 20000>2^14, W=R/20000<2^154. The submitted bound is 154 exactly. No
reference operation is additionally charged, so internal work is not double counted.

Peak storage: two arrays for mergesort, each m*8192=2^141 bytes, plus metadata,
code, counters, scratch and output, fit strictly within 2^143 bytes. There is no
random tape. Preprocessing (the table and sorting and setup) is a subset of the
charged total, and is bounded by 2^154 normalized equivalents. There is no
nonuniform searched advice: public profile constants and the algorithm are the
only embedded information. Report advice bound 1 byte (log bound 0); public
constants and code are included in peak storage. No future-target amortization,
parallel elapsed-time discount or physical allocation is assumed.

## Evidence scope and limitations

This proof-only package declares no computational experiments. In particular it
has not allocated the exponential table, executed the production attack, measured
an arithmetic kernel or demonstrated practical inversion. The correctness is the
deterministic BSGS argument and nondegenerate bilinearity of the specified pairing.
The instruction ledger is an intentionally conservative analytic implementation
bound. Structural checks and independent advisory review cannot certify platform
acceptance. Only the configured hosted verifier can accept and assign this claim.
