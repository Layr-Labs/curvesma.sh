# bls12-381-fixed-g2-invert-g1: deterministic bsgs accounting control

## Target and reduction

Given Q, compute h=e(P0,Q), solve h^x=z and return [x]P0 as a G1 point. Nondegeneracy makes h a generator of GT because the supplied point is nonidentity. Bilinearity proves the returned point pairs to z for the exact supplied point. For z=1 return infinity. The supplied point is never replaced. The recovered exponent is an intermediate value, not the required final output.

## Algorithm: deterministic sorted baby-step/giant-step

Set m=ceil(sqrt(r)), computed by integer arithmetic. Form all m records
(j*G,j), 0<=j<m, including the identity; in GT use G^j. Construct iteratively,
starting at identity. Sort by the full canonical representation with bottom-up
mergesort and a second array. Compute D=-(m*G) (or G^(-m)). Set V=target.
For i=0..ceil(r/m)-1, binary-search V in the sorted baby records. A match j
supplies x=i*m+j. Require x<r and verify the original equation; output x or
the required reconstructed pairing preimage. Otherwise update V by D.
For x in [0,r), division with remainder gives x=i*m+j in exactly this range,
so the algorithm always finds it. A wraparound match with x>=r is discarded.

Since r<2^256, m<=2^128 and the number of giants is <=m. Table construction
and giant steps use at most 2^129 group operations, each charged 2^48
instructions. Together with setup/finish they cost <2^178 instructions.
Mergesort uses <=128*m comparisons and <=256*m record copies, bounded together
by 2^137 record operations; binary searches add <=129*m comparisons, and
initialization adds <=m copies. Charge the combined record operations fewer
than 2^139*2^24=2^163 instructions. Including all costs gives <2^179
instructions, or <2^170 reference equivalents since C>2^9. Declare the looser
upper bound 2^180, not a measured running time. Preprocessing, including all
baby arithmetic and sorting, is <2^170 reference equivalents and included in
the total. Two arrays of <=2^128 records of <=8192 bytes require <=2^142 bytes;
code, scratch, input and counters fit within the declared peak <2^145 bytes.
Nonuniform advice is zero bytes. No sorting, lookup or memory-transfer cost is
omitted; this is deterministic and uses no ideal hash table.

## Arithmetic and instruction accounting

This is an analytic algorithm specification, not a timing claim about Python.
Use fixed-length bit arrays for residues; one bit occupies one byte, so even
schoolbook arithmetic has a simple scalar implementation. Each bit loop uses
loads, stores, shifts, boolean operations, carry arithmetic, comparisons and
branches. Charge at most 64 scalar AArch64 instructions per bit per pass, including
address computation. A 384-bit add/subtract/compare, with conditional reduction,
uses fewer than 16 such passes and hence fewer than 2^19 instructions. A modular
multiply uses 384 rounds of double-and-add, each with fewer than eight such
passes: fewer than 2^31 instructions. These bounds intentionally waste orders of
magnitude rather than assume optimized Montgomery arithmetic. Every residue is
canonical after the operation. Modular inversion of a nonzero residue uses binary
exponentiation to p-2, fewer than 768 modular multiplications, hence <2^41
instructions. To keep the later bound simple we charge **2^48 instructions per
complete curve step**, including inversion, special cases, comparisons and output
encoding. Use affine formulas: for distinct x, lambda=(y2-y1)/(x2-x1), x3=lambda^2-x1-x2,
y3=lambda*(x1-x3)-y1. For doubling use lambda=(3*x1^2+a)/(2*y1).
Handle infinity, inverse pairs, and y=0 explicitly before division. This applies
to both a=0 and P-256's a=-3; no secp256k1 microbenchmark is used to price P-256.

For GT multiplication use polynomial convolution (144 base-field products) and
reduce degrees 12..22 with w^12=2*w^6-2. This uses fewer than 2^10 base-field
multiplications/additions, below 2^42 instructions, again within the 2^48 step
allowance. Fp2 arithmetic uses u^2=-1. Its inversion is (a-b*u)/(a^2+b^2),
so G2 scalar multiplication can use the same generous per-step allowance.

Input decoding and subgroup checks, pairing computation when needed, computing
one inverse/stride, final scalar multiplication, output construction and final
verification are covered by a separate **2^80 instruction setup/finish allowance**.
A subgroup test or scalar multiplication uses fewer than 512 curve steps. An
Fp12 inverse can be computed by exponentiation to p^12-2 using fewer than 9216
Fp12 products. For pairings use exactly the pinned optimal-ate map and final
exponentiation, not a different pairing: a 64-bit Miller loop with bounded
line-function arithmetic, followed by exponentiation by (p^12-1)/r (less than
4572 bits). Even allowing 2^20 extension operations, each bounded by 2^58
instructions including inversions, this is below 2^78; the 2^80 allowance covers
two pairings plus the other fixed costs. These are executable arithmetic schedules,
not free oracles. There is no precomputed table beyond the public profile constants.

Tables use canonical full point coordinates with an infinity tag, or all twelve
GT coefficients, followed by the scalar index. A record fits in 8192 bytes even
with the deliberately inefficient bit-array arithmetic representation. Scalar
indices and addresses use multiple 64-bit words when necessary: never assume a
128-bit table fits in a real machine's 64-bit address space. The model is an
abstract collection of storage blocks, with logical addresses explicitly encoded
in up to 256 bits. Address calculation and memory transfers are charged. This is
the same theoretical large-memory interpretation required by generic BSGS;
physical availability of the memory is not claimed. Each record compare or copy,
including all address words and control, costs less than 2^24 instructions.
The claim counts total logical transfers in the instruction-equivalent model,
not latency or a claim about a physically realized interconnect.

## Correctness and scope

All groups used have the prime subgroup order r in the protected profile. The
public base is nonidentity, hence generates the relevant subgroup. Complete
arithmetic and exact comparison mean no hash collision assumptions are needed.
No probability heuristic, birthday approximation, free preprocessing, external
advice, parallel wall-clock discount or uncharged failed trial is used. The
algorithm succeeds on every valid challenge, giving success probability 1, above
0.39. Malformed inputs are outside the honest challenge distribution; the
specified canonical decoding and subgroup checks reject them.

This submission is a conservative baseline, not a new cryptanalytic result or a
best-known security estimate. The finite toy experiment checks group-search
control flow only. It cannot establish the cost bound or security of BLS12-381,
P-256 or secp256k1, and does not execute a production-size discrete logarithm.
The analytic instruction upper bound is intentionally loose. Memory is reported
separately and is not an extra objective. A production implementation would need
substantial engineering and, for BSGS, infeasible memory.
