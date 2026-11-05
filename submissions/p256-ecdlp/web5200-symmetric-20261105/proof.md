# P-256: deterministic symmetric interval search with a complete cold-work bound

## Claim and scope

For the exact full-size P-256 target in the protected profile, the algorithm
below returns the discrete logarithm with probability strictly greater than
1/2 over the uniform scalar in [0,r). Its worst-case total classical work,
including unsuccessful inputs, is at most 2^151 target-group-operation
equivalents under classical-group-work-v2 with C=1000. Peak memory is at
most 2^144 bytes. Preprocessing is included in total work and is itself
bounded by 2^151 units. There is no searched nonuniform advice and no
randomness, heuristic, quantum step, reusable free table, or production run.
Zero searched advice is conservatively reported by the conventional
log2(max(1,advice_bytes))=0.

This is an analytic, generic baby-step/giant-step construction, using equality
up to elliptic-curve negation and a deterministic partial-coverage interval.
It is not a new leading-exponent attack, not a practical break, and not a
claim to improve the best displayed research result. The deliberately loose
cost certificate uses scalar 64-bit instructions, not the 256-bit unit-word
ledger in the repository's old organizer control. Memory is astronomical
and explicitly part of this hypothetical construction. The attack has not
been executed on the participant's computer or elsewhere.

## 1. Fixed target and input validation

All numbers in this paragraph are hexadecimal. Use

    p = ffffffff00000001000000000000000000000000ffffffffffffffffffffffff
    r = ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
    a = ffffffff00000001000000000000000000000000fffffffffffffffffffffffc
    b = 5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
    P.x = 6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
    P.y = 4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5

The coefficient a is -3 modulo p, NOT zero. These are the protected
SEC 2 v2.0 Section 2.4.2 parameters. P has prime order r on
y^2=x^3+a*x+b over F_p, and the target is Q=[x]P for uniform 0<=x<r.
The specified group structure and public constants are the problem's
premises, not secret advice. In particular

    2^255+1 < r < 2^256,     r odd,     2^128 < r.

The first inequality follows already from the most significant hex digits:
r begins with ffffffff whereas 2^255+1 begins with 80000000. The other
inequalities are immediate from the bit length. No approximate order is
substituted into the correctness or success proof.

Read the exact UTF-8 JSON input with key Q. Infinity is null; other points
have exactly x,y, each 64 lowercase hex characters for an integer less
than p. Reject duplicate keys, noncanonical coordinates and off-curve
points, and verify [r]Q=O. A parser can operate on the fixed-size canonical
profile syntax, rejecting overlength or noncanonical input early. The work
claim concerns well-formed generated challenges. No attacker-chosen
unbounded whitespace or arbitrary input length is part of that distribution.
The input contains no powers, field oracle, scalar, or sampling coins.

If Q=O, return the canonical 32-byte scalar zero immediately. This is a
successful challenge, not a discarded sample. Otherwise perform the search
below. All comparisons use canonical field representatives.

## 2. Algorithm

Set m=2^127 and h=2^127. Every preprocessing step is performed afresh for
this one target. The scalar constants m,h have no hidden construction cost.

1. Starting with P, generate A_j=[j]P for j=1,...,m by successive ordinary
   affine additions of P. Store (A_j.x,A_j.y,j) as a fixed 128-byte record.
   Use canonical affine coordinates; no hash-function assumption is used.
   Retain A_m. Compute D=-[2]A_m.
2. Sort the m records by their x-coordinate with bottom-up merge sort and
   two m-record arrays. Keys have 256 bits. Store the full y-coordinate and
   index j so that the sign of a recovered scalar is unambiguous.
3. Put gamma=Q. For each i=0,...,h-1, do the following:
   - If gamma=O, set candidate c=2mi modulo r.
   - Otherwise binary-search the sorted x-keys. If none equals gamma.x,
     this iteration has no candidate.
   - On a match with record (A_j.x,A_j.y,j), choose epsilon=+1 when
     gamma.y=A_j.y, and epsilon=-1 when gamma.y=-A_j.y modulo p. If
     neither equality holds, reject the candidate. Set
     c=(2mi+epsilon*j) modulo r.
   - For any candidate, independently check [c]P=Q and 0<=c<r. On success
     emit the exact JSON output key x and a 64-character lowercase hex
     scalar. If the check unexpectedly fails, terminate with failure.
   - If no candidate was found and i<h-1, set gamma=gamma+D and continue.
4. If every iteration finishes without success, return failure. No restart,
   coins, amplification or second table is invoked.

The residue reduction in step 3 uses the SIGNED integer 2mi+epsilon*j.
Its range is [-m,2mh-m]. As 2mh=2^255<r, a negative value is reduced
by adding r and a nonnegative value is already less than r. Thus reduction
needs at most one correction, not a large division or free oracle.

The identity A_0 is not stored in the table; the explicit gamma=O branch
covers it. No A_j for 1<=j<=m is the identity because j<r.

## 3. Why x-coordinate lookup is exact

Two affine points on this short Weierstrass curve having the same
x-coordinate obey y_1^2=y_2^2. In the field F_p,
(y_1-y_2)(y_1+y_2)=0, so they are equal or negatives. This does not assume
random point encodings. Within our table the x-keys are distinct: if
A_j=+A_k then r divides j-k, forcing j=k in [1,m]; if A_j=-A_k then r
divides j+k, impossible since 2<=j+k<=2m=2^128<r.

A finite nonidentity element of the odd-order subgroup cannot equal its
negative: that would give order two. Therefore table y-coordinates are
nonzero and the signs in step 3 cannot both hold. The y-coordinate is
retained so there is no assumed free square root, no ambiguity from
point compression, and no probability of a spurious hash match.

By induction on i, gamma=Q-[2mi]P at the start of iteration i. Every hit
with j and sign epsilon implies Q=[2mi+epsilon*j]P. The returned scalar
is its unique representative in [0,r), and the explicit final scalar
multiplication checks the result. The infinity branch is the case j=0.

## 4. Exact support and success probability

At iteration i the set of scalar residues recovered is precisely the
image modulo r of the integer interval

    I_i = [2mi-m, 2mi+m].

Its two halves are epsilon=-1 and +1 with 1<=j<=m, and the center is
the infinity branch. Adjacent intervals meet at one endpoint because
the centers are 2m apart. Their union for 0<=i<h is the consecutive
integer interval

    I = [-m, 2mh-m],      |I|=2mh+1=2^255+1.

The difference between its largest and smallest integer is 2^255<r, so
reduction modulo r is injective on I. Thus the search succeeds on exactly
2^255+1 distinct scalars out of the r possible ones, including scalar 0.
Its success probability is exactly (2^255+1)/r, which is strictly greater
than 1/2 since r<2^256, and therefore exceeds 0.39. These are uniform
full-production challenges; no distribution conditioning or key-range
promise is made. The other scalars simply produce failure after the fixed
budget. All such failures are included in the worst-case work below.

## 5. Executable arithmetic specification

See arithmetic.md for the detailed scalar-instruction and memory ledger.
Every field element is FOUR 64-bit limbs, with a fifth limb used for carry
and subtraction correction. A 256-bit integer is never treated as a
single register or a single-cost instruction.

Use the following complete affine group law for U=(x,y), V=(u,v):
return the other point for an identity operand; if x=u and y+v=0 modulo p,
return O; otherwise take lambda=(v-y)/(u-x) for x!=u, or
lambda=(3x^2+a)/(2y) for doubling. Then return
(lambda^2-x-u, lambda*(x-x')-y). The possible zero denominator at doubling
has already returned O. The numerator includes the actual a=-3.
Infinity is a separate tag and is not passed to coordinate inversion.
The formulas are the ordinary chord/tangent group law from the target
curve, not an incomplete affine formula counted as a complete operation.

Modular multiplication uses fixed 256-bit shift-and-add; inversion of a
nonzero z uses binary powering z^(p-2). Fermat's identity follows by
permuting the nonzero field elements under multiplication by z and
cancelling their product, so z^(p-1)=1. Scalar multiplication is ordinary
binary double-and-add with at most 512 point operations for a 256-bit
nonnegative scalar. This also specifies both subgroup checking and final
candidate verification. The implementation is slow but concrete.

## 6. All work, memory, preprocessing and advice

The ledger proves an upper bound of fewer than 2^160 scalar AArch64
instruction equivalents for every generated challenge, including a
full unsuccessful run. With no separately charged reference point
operations, H=0 and W=R/1000. Since 1000>2^9,

    W < 2^160/1000 < 2^151.

We submit the conservative exact rational score 151, meaning W<=2^151;
we do not claim equality with a measured runtime. The upper bound is
deliberately rounded upward, as the review checks the submitted bound.

Preprocessing is a subset of this computation and is also bounded by
2^151 equivalents. It includes baby generation, step construction,
both table buffers and their initialization, sorting, directory/routing
initialization, and all advice-like structures. None is reused across
future challenges. Peak memory, including distributed routing state,
code, tables and scratch, is at most 2^144 bytes. The materialized baby
and sorting tables really do require exponential storage.

The algorithm uses no random bits. There is no quantum computation,
external advice, uncharged oracle, free global lookup, cryptanalytic
heuristic, or uncharged precomputation. Uniform input generation is
organizer work, while input parsing, subgroup testing, scalar recovery,
the final attack-side multiplication and output are all charged here.
No asymptotic big-O statement is used as a numerical coefficient.

## 7. Evidence, comparisons and provenance

This proof and arithmetic.md are the entire mathematical evidence. No
experiments are declared. No full-size scalar recovery was performed, and
passing structural validation cannot establish any practical production
attack. Independent AI critique and platform exploratory review are not
formal machine proofs or independent human approval.

Baby-step/giant-step and elliptic-curve negation are classical methods.
The organizer secp256k1 control at the same repository commit was read
for orientation; it is not an accepted prior participant result and its
numerical ledger is not adopted. In particular that old control used a
different word-width discussion and cannot price P-256 addition. The
64-bit implementation here is explicitly rederived, includes a=-3,
and charges global addressing rather than assuming 256-bit unit-cost RAM.
The contribution is a complete exact-support and cold-cost certificate
for this particular deterministic half-coverage construction, not a claim
of a new generic DLP algorithm. No prior participant proof, implementation,
or submission-specific optimization was used. Public board metadata was
viewed only to confirm that stronger results already exist.

Protected context: Layr-Labs/curvesma.sh commit
6423c2728dc03d5f2d647c7f22adf1dac0643a66, tracks/p256-ecdlp/{rules.md,
target-profile.json,verifier.json}, tracks/classical-group-work-v2.json,
and tracks/docs/COST_MODEL.md. The complete construction and reasoning
above do not require the reviewer to retrieve external references.
