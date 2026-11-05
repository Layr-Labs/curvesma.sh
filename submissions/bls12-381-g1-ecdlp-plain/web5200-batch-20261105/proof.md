# BLS12-381 G1: deterministic symmetric interval search with a complete cold-work bound

## Claim and scope

For the exact full-size BLS12-381 G1 target in the protected profile, the algorithm
below returns the discrete logarithm with probability strictly greater than
1/2 over the uniform scalar in [0,r). Its worst-case total classical work,
including unsuccessful inputs, is at most 2^152 target-group-operation
equivalents under classical-group-work-v2 with C=2000. Peak memory is at
most 2^144 bytes. Preprocessing is included in total work and is itself
bounded by 2^152 units. There is no searched nonuniform advice and no
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

All numbers in the following parameter display are hexadecimal:

    p = 1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
    r = 73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
    a = 0
    b = 4
    P.x = 17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb
    P.y = 08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1

Use the exact protected BLS12-381 G1 profile. P denotes P0
of prime order r on y^2=x^3+a*x+b over F_p. The input point is
T=[x]P for uniform 0<=x<r; we use Q as its mathematical
name throughout the search below. The G1 cofactor is not one; subgroup validation is mandatory.
We use no auxiliary powers, SRS information, pairing or private sampling coins.
The prime field and specified subgroup order are protected problem premises.
In particular

    2^254+1 < r < 2^255,     r odd,     2^128 < r.

These inequalities follow directly from the displayed hexadecimal order:
its leading 0x73 lies strictly between 0x40 and 0x80.
No approximate order is used in the correctness or support argument.

For the finite I/O cost certificate, fix a canonical serializer for the
mathematical point input: the outer key is T; infinity is null;
otherwise the inner object has x followed by y, each exactly 96
lowercase hexadecimal characters, with no redundant whitespace or escapes.
This is a representation convention, not a restriction on x or on the
uniform group element: every production point, including infinity, has
exactly one such representation, of fewer than 512 bytes. Decode the
coordinates as canonical integers below p, reject duplicate or extra keys,
off-curve values and noncanonical coordinate strings, and verify [r]Q=O.
The target point is denoted Q in the mathematics even when its input key
is T. A finite-state parser for this bounded grammar makes at most two
passes at 64 scalar instructions per byte plus 2^16 fixed control work;
this is inside the fixed setup reserve in arithmetic.md.

The scored bound is for that explicitly specified bounded representation
of uniformly sampled production points. The protected profile does not
specify a maximum amount of redundant JSON whitespace. No assertion is
made that every arbitrarily long JSON spelling can be consumed at the
same constant cost: accepting a spelling of L bytes by first normalizing
its redundant syntax adds at most 128L+2^16 parsing instructions, which
must be charged if such an expanded wire format is part of an interface.
The submitted score here concerns the fixed canonical representation;
its use by an actual platform transport remains an interface assumption,
not an authenticated fact inferred from the profile. The mathematical
support proof below holds for all of the stated group elements regardless
of how their encoding is written.

If Q=O, return the canonical 32-byte scalar zero immediately. This is a
successful challenge, not a discarded sample. Otherwise perform the search
below. All comparisons use canonical field representatives.

## 2. Algorithm

Set m=2^127 and h=2^126. Every preprocessing step is performed afresh for
this one target. The scalar constants m,h have no hidden construction cost.

1. Starting with P, generate A_j=[j]P for j=1,...,m by successive ordinary
   affine additions of P. Store (A_j.x,A_j.y,j) as a fixed 128-byte record.
   Use canonical affine coordinates; no hash-function assumption is used.
   Retain A_m. Compute D=-[2]A_m.
2. Sort the m records by their x-coordinate with bottom-up merge sort and
   two m-record arrays. Keys have 381 significant bits in six 64-bit limbs. Store the full y-coordinate and
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
Its range is [-m,2mh-m]. As 2mh=2^254<r, a negative value is reduced
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

    I = [-m, 2mh-m],      |I|=2mh+1=2^254+1.

The difference between its largest and smallest integer is 2^254<r, so
reduction modulo r is injective on I. Thus the search succeeds on exactly
2^254+1 distinct scalars out of the r possible ones, including scalar 0.
Its success probability is exactly (2^254+1)/r, which is strictly greater
than 1/2 since r<2^255, and therefore exceeds 0.39. These are uniform
full-production challenges; no distribution conditioning or key-range
promise is made. The other scalars simply produce failure after the fixed
budget. All such failures are included in the worst-case work below.

## 5. Executable arithmetic specification

See arithmetic.md for the detailed scalar-instruction and memory ledger.
Every field element is SIX 64-bit limbs, with a seventh limb used for carry
and subtraction correction. A 381-bit field integer is never treated as a
single register or a single-cost instruction.

Use the following complete affine group law for U=(x,y), V=(u,v):
return the other point for an identity operand; if x=u and y+v=0 modulo p,
return O; otherwise take lambda=(v-y)/(u-x) for x!=u, or
lambda=(3x^2+a)/(2y) for doubling. Then return
(lambda^2-x-u, lambda*(x-x')-y). The possible zero denominator at doubling
has already returned O. The numerator includes the actual a=0.
Infinity is a separate tag and is not passed to coordinate inversion.
The formulas are the ordinary chord/tangent group law from the target
curve, not an incomplete affine formula counted as a complete operation.

Modular multiplication uses fixed 384-bit shift-and-add with leading zero bits included; inversion of a
nonzero z uses binary powering z^(p-2). Fermat's identity follows by
permuting the nonzero field elements under multiplication by z and
cancelling their product, so z^(p-1)=1. Scalar multiplication is ordinary
binary double-and-add with at most 512 point operations for a 256-bit
nonnegative scalar. This also specifies both subgroup checking and final
candidate verification. The implementation is slow but concrete.

## 6. All work, memory, preprocessing and advice

The ledger proves an upper bound of fewer than 2^162 scalar AArch64
instruction equivalents for every generated challenge, including a
full unsuccessful run. With no separately charged reference point
operations, H=0 and W=R/2000. Since 2000>2^10,

    W < 2^162/2000 < 2^152.

We submit the conservative exact rational score 152, meaning W<=2^152;
we do not claim equality with a measured runtime. The upper bound is
deliberately rounded upward, as the review checks the submitted bound.

Preprocessing is a subset of this computation and is also bounded by
2^152 equivalents. It includes baby generation, step construction,
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

This proof and arithmetic.md are the complete mathematical evidence. No
experiments are declared and no full-size scalar recovery was performed.
Mechanical checks and advisory critique are not formal proof or human approval.

This entry adapts our own earlier P-256 submission
e428ce82-2ee3-48a1-a863-cf71d01a52da, commit
5e2a880fa867b6ac772c7661c768c462aa336ce9, files
submissions/p256-ecdlp/web5200-symmetric-20261105/proof.md and arithmetic.md.
We reuse its exact signed-interval BSGS argument and charged local-mailbox
storage implementation and declare that dependency. Its accepted judgment
supplies no mathematical premise; every needed argument is reproduced here.

The adaptation binds the actual BLS12-381 G1 prime field, subgroup, generator,
input syntax, a=0, b=4, and normalization C=2000.
The field grows to six limbs, the nontrivial cofactor makes subgroup checking essential, h is halved to 2^126, and all field arithmetic and total-work bounds are rederived.
Public board notes were viewed for comparison, not as accepted premises or
as the source of this algorithm. Lower submitted costs already exist; this
generic analysis does not claim a frontier improvement, new attack or novelty.

Protected source: Layr-Labs/curvesma.sh commit
6423c2728dc03d5f2d647c7f22adf1dac0643a66, tracks/bls12-381-g1-ecdlp-plain/rules.md,
target-profile.json and verifier.json, tracks/classical-group-work-v2.json,
and tracks/docs/COST_MODEL.md. The hosted verifier is pinned to
8c1f098875fae0750c29ecfd641f7c9fc53a0f41. No web retrieval is needed to
assess the complete argument in the supplied files.
