# p256-ecdlp: deterministic enumeration accounting control

## Target and reduction

Use the exact standard generator and subgroup order from target-profile.json. Solve the supplied point as [x]G and return x in canonical 32-byte lowercase hexadecimal encoding. Return zero exactly for infinity.

## Algorithm: deterministic enumeration

Start V at the identity and j=0. Compare V to the target. If equal, return j
(and reconstruct the output point for pairing inversion). Otherwise update
V=V*G, or V=V+G in additive notation, increment j, and repeat until j=r-1 has
been compared. Thus every exponent 0..r-1 is covered exactly once. No table is
stored. Verify the returned answer independently before encoding it.

At most r<2^256 group steps/comparisons occur. Charge each complete iteration
2^48 instructions as above; this includes scalar counter update and comparison.
Add setup/finish 2^80. Total instructions are <2^305. Every reference price C
is at least 1000>2^9, so work is <2^296 reference equivalents. The declared
log2 upper bound is 296. Setup/preprocessing is at most 2^80/C<2^71 reference
equivalents and is included. Peak storage is <2^30 bytes including code,
input window, arithmetic temporaries and fixed constants. No large table and
no stored random coins are needed; nonuniform advice is zero bytes.

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
