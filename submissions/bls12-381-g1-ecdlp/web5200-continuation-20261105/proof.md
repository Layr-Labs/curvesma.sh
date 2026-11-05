# BLS12-381 G1 with auxiliary powers: deterministic symmetric interval search with a complete cold-work bound

## Claim and scope

For the exact full-size BLS12-381 G1 target in the protected profile, the algorithm
below returns the discrete logarithm with probability strictly greater than
1/2 over the uniform nonzero secret s in [1,r). Its worst-case total classical work,
including unsuccessful inputs, is bounded by W(n) in Section 6 and, at the required scoring n=2^20, is at most 2^152 target-group-operation
equivalents under classical-group-work-v2 with C=2000. Peak memory is at
most 2^145 bytes at n=2^20, including the compact supplied input. Preprocessing is included in total work and is itself
bounded by W(n), and by 2^152 units at n=2^20. There is no searched nonuniform advice and no
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
T=[s]P for uniform 1<=s<r; we use Q as its mathematical
name throughout the search below. The G1 cofactor is not one; subgroup validation is mandatory.
The additional public inputs are n and A_i=[s^i]P0 for i=1,...,n, with T=A_1. The same secret s is used throughout. We read the array as described below but do not exploit its powers in the search; the protected access rule permits using only a subset. No pairing or private sampling coins are used.
The prime field and specified subgroup order are protected problem premises.
In particular

    2^254+1 < r < 2^255,     r odd,     2^128 < r.

These inequalities follow directly from the displayed hexadecimal order:
its leading 0x73 lies strictly between 0x40 and 0x80.
No approximate order is used in the correctness or support argument.

For the finite I/O certificate, use compact canonical UTF-8 JSON with
exact outer keys n,A,T in that order, no redundant whitespace/escapes,
and inner point keys x,y with 96 lowercase hex characters each. The
integer n satisfies 1<=n<r and has at most 77 decimal digits. A has exactly
n points in increasing index order; T equals A_1 on the honest distribution.
Every permissible mathematical instance has this representation, of less
than 512(n+1) bytes. This is an explicit serialization assumption, not an
authenticated fact about a deployed transport or a restriction on s.

Stream the input one point at a time. Retain n, a four-limb index and A_1;
check each coordinate, curve equation and [r]A_i=O and discard it after
processing, except A_1. Read T, validate its subgroup and nonidentity,
and check T=A_1. The powers relation is an honest-input promise; we neither
reconstruct s to validate it nor treat any supplied exponent as secret advice.
The profile explicitly permits not reading or validating all powers; this
construction chooses to do the inexpensive-in-n validation above for a
complete stream-I/O bound. It then uses only T for its search.

At most two byte visits at 64 instructions each plus 2^16 per point cover
the fixed grammar, hex decoding and counters. Group validation uses at
most 512 complete G1 operations plus a fixed number of field operations.
The conservative per-point allowance 2^45 instructions covers all these
steps, ordinary stream transfers, local parsing buffers and loop control.
No n-point working table is retained. Incoming compact input itself is
charged as at most 512(n+1) bytes in the total-memory bound. Indexing needs
at most four/five limbs for the entire permitted n range. We use the same
bounded local stream/receive convention as the declared machine below.

An expanded JSON spelling of L bytes needs an additional bounded streaming
normalization charge 128L+2^16 if that spelling, rather than the specified
compact one, is an interface input. The fixed scoring bound does not assert
constant-cost handling of arbitrary redundant whitespace. All representation
and local-event assumptions are explicitly part of the declared claim.

If T=O, reject: the sampled secret is nonzero, so this case is outside the
honest challenge distribution. Set Q=T for the mathematics below. The
output is a scalar under the JSON key s, not a point or a new power.

## 2. Algorithm

Set m=2^127 and h=2^126. Every preprocessing step is performed afresh for
this one target. The scalar constants m,h have no hidden construction cost.

1. Starting with P, generate B_j=[j]P for j=1,...,m by successive ordinary
   affine additions of P. Store (B_j.x,B_j.y,j) as a fixed 128-byte record.
   Use canonical affine coordinates; no hash-function assumption is used.
   Retain B_m. Compute D=-[2]B_m.
2. Sort the m records by their x-coordinate with bottom-up merge sort and
   two m-record arrays. Keys have 381 significant bits in six 64-bit limbs. Store the full y-coordinate and
   index j so that the sign of a recovered scalar is unambiguous.
3. Put gamma=Q. For each i=0,...,h-1, do the following:
   - If gamma=O, set candidate c=2mi modulo r.
   - Otherwise binary-search the sorted x-keys. If none equals gamma.x,
     this iteration has no candidate.
   - On a match with record (B_j.x,B_j.y,j), choose epsilon=+1 when
     gamma.y=B_j.y, and epsilon=-1 when gamma.y=-B_j.y modulo p. If
     neither equality holds, reject the candidate. Set
     c=(2mi+epsilon*j) modulo r.
   - For any candidate, independently check [c]P=Q and 1<=c<r. On success
     emit the exact JSON output key s and a 64-character lowercase hex
     scalar. If the check unexpectedly fails, terminate with failure.
   - If no candidate was found and i<h-1, set gamma=gamma+D and continue.
4. If every iteration finishes without success, return failure. No restart,
   coins, amplification or second table is invoked.

The residue reduction in step 3 uses the SIGNED integer 2mi+epsilon*j.
Its range is [-m,2mh-m]. As 2mh=2^254<r, a negative value is reduced
by adding r and a nonnegative value is already less than r. Thus reduction
needs at most one correction, not a large division or free oracle.

The baby identity B_0=O is not stored in the table; the explicit gamma=O branch
covers it. No B_j for 1<=j<=m is the identity because j<r.

## 3. Why x-coordinate lookup is exact

Two affine points on this short Weierstrass curve having the same
x-coordinate obey y_1^2=y_2^2. In the field F_p,
(y_1-y_2)(y_1+y_2)=0, so they are equal or negatives. This does not assume
random point encodings. Within our table the x-keys are distinct: if
B_j=+B_k then r divides j-k, forcing j=k in [1,m]; if B_j=-B_k then r
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

For giant index i, equality up to sign with baby j=1,...,m, together with
the gamma=O center, recovers the residue image of

    I_i=[2mi-m,2mi+m].

Adjacent intervals share their endpoints. Their union for 0<=i<h is

    I=[-m,2mh-m],   |I|=2mh+1=2^254+1.

Its span 2^254 is less than r, so residues are injective. Exactly one of
these residues is zero. Since the target secret is uniform on 1,...,r-1,
the supported NONZERO residue count is exactly 2^254, and

    Pr(success)=2^254/(r-1)>1/2>0.39.

The argument holds for every 1<=n<r without any divisor or smoothness
condition. The public powers are correlated with s exactly as specified,
but the algorithm uses no independent-randomness premise about them.
Removing scalar zero changes the numerator and denominator; we do not
carry over the plain-G1 probability unchanged. Every hit on an honest
nonzero target produces a nonzero canonical candidate and passes final
verification. All other secrets produce an explicit fixed-budget failure.

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

## 6. Work as a function of n

The search arithmetic and routed-storage certificate is reproduced in
arithmetic.md. Define the following explicit instruction upper bound:

    R(n)=3*2^160+2^149+2^140+2^52+(n+1)*2^45,
    W(n)=R(n)/2000,                      1<=n<r.

W(n) bounds total classical work in the protected reference units; no
separate reference group operation is charged (H=0). The first four
terms bound the entire cold BSGS construction, including failures,
bootstrap and final verification. The last term covers streaming and
validation of n powers and T, including all accesses. Input generation
is organizer work and is not charged again. The formula has no omitted
n-dependent assumption and needs no number-theoretic condition on n.

At the REQUIRED scoring n=2^20, (n+1)*2^45<2^66. Therefore

    R(2^20)<2^162,  W(2^20)<2^162/2000<2^152.

The first claim submits the exact conservative rational 152/1 for that
scoring point. This is an upper bound, not an exact runtime or an evaluated
score at a convenient smaller n. Preprocessing is charged cold and is
bounded by W(n), and by 2^152 at the scoring point.

Search working storage, including replicated routing code and controller
state, is below 2^144 bytes. Counting the entire compact input as well,

    M(n)<=2^144+512(n+1) bytes,
    M(2^20)<2^145 bytes.

Input can be consumed sequentially with one point buffer; no extra n-point
copy is needed. Any external input storage is counted by the displayed
addend, while working addresses are routed as in arithmetic.md. The
construction assumes the explicitly declared compact input stream and
bounded local mailbox/reset interface, not free global array addressing.
Searched nonuniform advice is zero. No random bits, restarts, amplification,
quantum computation, real attack or uncharged precomputation are used.

## 7. Attribution and scope

This adapts our own plain-G1 submission e25b1379-d936-4940-9739-d4d90fe30718,
commit 6c39d6f28139142dfadd08ddc5fc957ab062a7a7. We reuse and reproduce
its signed-interval argument, six-limb arithmetic and charged local-mailbox
storage construction. Earlier acceptance is not a mathematical premise.
The changed target uses a nonzero s and an honest correlated powers array.
This manuscript rederives exact nonzero coverage, outer input/output keys,
all input-dependent work W(n), and the mandated n=2^20 evaluation.

This is a generic accounting control which deliberately ignores the
algebraic acceleration offered by the auxiliary powers. It does not
claim a Cheon improvement, new cryptanalysis, a competitive record or a
practical break. No new subgroup/input promises on s are introduced.
The source is curvesma.sh commit 6423c2728dc03d5f2d647c7f22adf1dac0643a66,
the auxiliary-powers track profile/rules and classical-group-work-v2.
The finite input encoding and hypothetical machine assumptions remain
explicit, and exploratory review is distinct from human acceptance.
