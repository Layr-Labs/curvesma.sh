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

Store each Fp coefficient in six 64-bit limbs. Modular addition/subtraction
uses at most six passes over seven limbs (including carry) at sixteen
instructions per limb, plus 1376 for correction, copies, calls, addresses
and control. Thus A=2048. Canonical operands need one modular correction.
Multiplication uses 384 conditional-add/double iterations costing 2A+64,
plus 2048 setup: M=384*(4096+64)+2048=1599488<2^21. This is a finite
limb implementation, with no arbitrary-width unit-cost multiplication.

Fp12 polynomial convolution uses 144 M and 188 A including descending
reduction of degrees 22 through 12 using w^12=2*w^6-2. Add 65536 local
instructions for copies, control and scratch. Total is below 2^28. Use the
looser E=2^30 for each complete multiplication, including three remote
record transfers at less than 2^22 each and conversion. Table keys can
remain 381 one-byte bits per coefficient: conversion of 12*384 bits costs
below 2^20 instructions, and records remain 8192 bytes. Fp12 inversion by
powering to p^12-2 uses fewer than 2^14 products. Fp2 multiplication/addition
fit E; inverse (a-bi)/(a^2+b^2) costs less than 2^33 instructions.

The explicit routing tree in routing.md supplies multiword addresses,
initialization, transfer, code replication and memory bounds. Comparison
or copying costs less than 2^22 instructions per record. This is an analytic
instruction-equivalent machine, not a production cluster.

This arithmetic and bounded routing-interface design materially use public
submission f9b7a390-0898-4369-a97c-3d752133ef4d. The required arithmetic is
reproduced here; its acceptance supplies no mathematical premise. Our earlier
drafts used slower bit arithmetic. This revision retains our larger records
and routing schedule and covers the complete exponent interval instead of a
partial rectangle. No NFS heuristic is used.

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
than four times 2^32*2^30 and includes canonical JSON input/output, with only
fixed-size production objects. This is polynomial setup, not an uncharged oracle.

## Complete work bound

The order satisfies r<2^255, so ceil(r/m)<=2^127. Table construction has fewer
than 2^128 Fp12 multiplications. Giant steps use fewer than 2^127. Both together
are less than 2^129 multiplications, hence fewer than 2^159 instructions using
the conservative 2^30 price. Table sorting takes 128 passes over m records.
Each pass performs at most m comparisons and at most 2m copies. Even charging
2^22 to each operation, with index/control allowances, sorting is below 2^160
instructions. Binary searches have at most 129 comparisons each, costing below
2^158 instructions in all; giant-loop scalar/control work is included by that
allowance. Setup and verification add less than 2^74. Other allocation,
initialization, canonical encoding, table reads/writes and cleanup are covered
by an extra 2^160 instructions, including freeing the table if required.
Therefore total R is strictly less than 2^162 scalar instruction equivalents.
Since 20000>2^14, W=R/20000<2^148. The submitted bound is 148 exactly. No
reference operation is additionally charged, so internal work is not double counted.

Peak storage: two arrays for mergesort, each m*8192=2^141 bytes, plus metadata,
code, counters, scratch and output, fit strictly within 2^147 bytes. There is no
random tape. Preprocessing (the table and sorting and setup) is a subset of the
charged total, and is bounded by 2^148 normalized equivalents. There is no
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
