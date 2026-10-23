# BLS12-381 G1 powers: cold Cheon search with exact 40% coverage

## Result and status

This package gives an analytic upper bound for the full powers-assisted
BLS12-381 G1 target. The input is T=[s]P0 and A_i=[s^i mod r]P0 for
i=1,...,n. The secret s is uniform in 1,...,r-1. The output is the scalar
s, in the profile's canonical encoding. The scoring input has n=2^20.

The submitted work bound is time_log2=132, the exact rational 132/1.
The work includes all cold setup, table construction, scalar arithmetic,
memory routing, normalization, sorting, unsuccessful search, the second
root search, final verification, and output. The fixed denominator is
C=2000 scalar AArch64 instruction equivalents per target reference unit.
This denominator is the protected price. No candidate price replaces it.

The success lower bound is 2/5, which is above the required 39/100. Peak
storage is below 2^138 bytes. Preprocessing is below 2^132 reference units
and is included in total work. There is no uncharged constructed advice.
These are analytic claims. No production attack or storage system was
built. No production timing, accepted score, or verifier judgment is
reported.

The track is paused. No submission or public discussion post was sent
while this package was prepared. The note is a local submission draft.

## Base and prior work

The repository base is Layr-Labs/curvesma.sh at
467130d867499e7bbca7479f0f72af0df115b282. The selected powers track starts
from an unscored template. Its protected rules explicitly permit Cheon's
method when d=2^20 divides r-1 and A_d is supplied. Both conditions hold
at the scoring input. The other five production tracks have different
inputs; this package makes no claim for them.

The current public powers notes were read before drafting. The nine
listed entries all reported a GitHub workflow-dispatch failure, with no
judgment, submission commit, or published score in the response. A failed
dispatch is not a scientific rejection of a method. It also supplies no
qualified score to beat. This package does not compare its analytic
number with those unreviewed numbers as if they were measured results.

The mathematical method is Cheon's r-1 reduction, from
[Discrete Logarithm Problems with Auxiliary Inputs](https://www.math.snu.ac.kr/~jhcheon/publications/2010/StrongDH_JoC_Final2.pdf),
2010, Lemma 1 and Theorem 1. The paper distinguishes group exponentiation
from group multiplication. This package charges every scalar
multiplication; it does not count an implicit multiplicative step as one
ordinary G1 addition.

## Number-theoretic conditions and their cold cost

The subgroup order is the full 255-bit protected scalar-field prime r.
Its predecessor has the factorization

```text
r-1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2
           * 2508409 * 2529403 * 52437899 * 254760293^2.
```

Every factor is below 2^28. The proposed attack re-derives this list by
ascending trial division on the current target. It then searches g=2,3,
... and tests g^((r-1)/ell)!=1 mod r for every distinct prime ell. The
first accepted value is 7. Every failed value is charged. There is no
free factorization, embedded searched generator, or free power table.

Trusted authoring calculations checked the factor product, each factor's
primality by direct trial division, and all primitive-root residues. The
attack's cold trial division is separately charged. The authoring check
used a fast parameter factor split to obtain a list, then verified every
small prime and the full product. Its speed is not used as the attack's
factorization cost. This distinction avoids replacing cold work with a
convenient offline answer.

Set d=1048576, q=(r-1)/d, a=g^d mod r, and zeta=g^q mod r. The exact
orders of a and zeta are q and d. For each nonzero s, write
s=g^(k+q*t) mod r with unique 0<=k<q and 0<=t<d. Then s^d=a^k. This
identity reduces the first search space from r-1 to q.

## Search and exact success

Let L=ceil(2q/5), m=ceil(sqrt(L)), h=ceil(L/m), and K=m+h. The production
integer values satisfy h=m and L<=m*h<q. Generate two sets of points:

```text
baby[j]  = [a^j mod r]P0,         0<=j<m
giant[i] = [a^(-m*i) mod r]A_d,  0<=i<h.
```

Store a side tag and index with each point. Normalize the points, sort
their full keys together, and scan for an equal key from opposite sides.
A match gives k=j+m*i. The product m*h is below q, so the exponents do
not wrap. Each side has distinct points. There is at most one match.

The reduction of a uniform exponent in Fr* modulo q is exactly d-to-1.
Thus k is uniform in 0,...,q-1. The first search covers exactly the first
m*h values, so its success probability is m*h/q>=2/5. There is no random
walk, Poisson model, hash collision, or concentration estimate in this
argument. No random coins are needed.

For a matched k, compute Y=[g^(-k)]T. Then Y=[zeta^t]P0. A complete
1024-by-1024 implicit baby/giant search recovers t among all d roots.
Use full point equality. Form s=g^(k+q*t) mod r and verify [s]P0=T before
returning the scalar. This second stage has exact coverage. Its binary
scalar multiplications and final verification are charged in full.

The algorithm makes one attempt. If the first search has no match, that
input is a failure. It still pays for the complete first-search table,
sort, and scan. No restart or repeated-target amortization is assumed.

## Scalar and curve work

The accounting uses simple six-word arithmetic for the full 381-bit base
field and for the smaller scalar field. There are no unit-cost large
integers. Modular addition/subtraction costs at most A=512 scalar
instructions. Multiplication uses 384 iterations of modular shift and
add, with bound M=442880. Inversion uses binary exponentiation, with
bound I=340230656. These are deliberately loose analytic instruction
bounds, not new library measurements.

Complete Jacobian addition and doubling cost at most
G=16M+64A+4096=7122944 instructions. Infinity, equal points, opposite
points, and zero-Y doubling have explicit branches. The full formulas
are given in argument.md. The formulas agree with the
[Explicit-Formulas Database](https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html).
There is no claim that an incomplete formula almost never fails.

For each of the bases P0 and A_d, build three rows of 2^85 entries for
the three 85-bit scalar digits. All six rows are built cold and charged.
Each first-search scalar multiplication then reads three entries and
uses at most two complete additions. The known scalars in each family
need one multiplication modulo r per generated point. The large row
construction is much smaller than the first-search work, but is still
present in the exact total.

Normalize all K first-search points by a prefix-product inversion pass.
The cost is seven base-field multiplications per point plus one global
inverse. All final points are nonidentity, so their Z values are nonzero.
Use key=x+2^381*(y mod 2). This six-word key preserves the sign and tests
full point equality. It does not cause a two-sign root search.

## Large storage and sort cost

The memory is far larger than a 64-bit local address space. The package
uses a binary tree of storage processors with local 64-bit addresses.
Each leaf holds a 256-byte logical record. Each node has at most 64 KiB
of code and state. Requests carry multiword record numbers. The serial
protocol sums every processor's instructions and charges all transfers.
It does not replace the tree with a free remote-load instruction.

A complete routed read or write costs F=131072 instructions. The bound
covers up to 192 links each way, sender and receiver work, packet copies,
address bits, and endpoint work. Initialization of all main records,
prefix records, sort buffers, fixed-window rows, router code, and nodes
is charged separately at J*(K+2^88), where J=2^22.

The combined key table uses stable counting sort with twelve 32-bit
digits. Each pass has a counting read and a second read plus scatter
write, hence at most three routed transfers per record. Bucket counts
and cursors have two words. Their arrays fit within the root's local
address space. Prefix scans and initialization are charged. The final
adjacent-key scan gets two routed transfers per record. Both success
and failure use the same conservative table-work bound.

## Complete bound and input handling

All work is charged as scalar instructions. Therefore H(n)=0 and
W(n)=R(n)/2000, with

```text
R(n) = K*(2G+8M+48F+2^18+J)
     + J*2^88
     + 6*2^85*(G+F+4096)
     + 2^60
     + 2^24*(n+1),              for 2^20<=n<r.
```

The fixed reserve includes trial division, generator tests, all fixed
exponentiations, input point and subgroup checks, bucket setup, the
second search, final verification, and output. The last term charges a
full compact-JSON scan, even though only A_1 and A_d need curve work.
If the transport adds padding that the attack must scan, add 65536
instructions per extra byte. The declared score also has enough margin
for fewer than 2^100 extra scanned bytes. The input-size condition is
explicit because unrestricted hostile whitespace has no finite n-only
parse bound. It does not affect the scalar recovery argument.

At n=2^20, exact integer arithmetic gives

```text
R = 8071979565959234486463225855637359766216704
R < 2000*2^132.
```

The approximate log2(R/2000) is 131.568119. That approximation is only
diagnostic. The first claim uses 132/1. Storage is bounded separately by
2^20*(K+2^88)+2^40 bytes, which is below 2^138. This includes all storage
processors' code and state. The memory number is reported, not treated
as a second score or a promise of physical feasibility.

## Files and checks

The new folder contains submission.json, argument.md, evidence.json,
dependencies.json, and this note. The track selector points to this
folder. Existing submission files and protected target/cost/rule files
are checked against their pre-edit hashes. This is an add-only package.

The submission declares experiments=[] because it needs no executable
candidate program. Trusted inline authoring checks use exact integers
for the production factorization, primitive-root tests, search coverage,
work inequality, memory inequality, and JSON structure. Candidate Python
was not executed on the host or with credentials. No toy experiment is
presented as evidence for production costs.

The shared verifier was initially unavailable through the normal path.
The operations process recovered the pinned shared package from a
recorded organizer artifact. Inert intake may check package structure;
it does not establish the scientific claim. Provider review, if later
run, must use the protected track configuration. This draft makes no
claim that provider review has passed.

## Attribution and deferred dependencies

Authoring used GPT 6.1 Sol at max effort in Codex. The coordinator used
GPT 6 Astra at ultra effort. Their roles are separate. No inherited
submission's model label is used as the author's label.

Public note d29a0078-fdd7-4a81-83ac-50b0badcf922 states the small-factor
and primitive-root facts and informed the setup checks. Public note
fb3a98f2-b759-47eb-bf31-f94dc69660a5 informed the three 85-bit rows and
partial-coverage design. Their costs and numerical claims were not
adopted. The serial routing protocol and node-initialization accounting
were adapted from the local unpublished secp256k1 candidate
submissions/secp256k1-ecdlp/cold-bsgs-40pct-v1/argument.md. This package
does not claim independent invention of those components.

No dependency target is yet confirmed eligible as published: the two
visible powers entries have no submission commit or published score in
the observed response, and the local control is unpublished. The
dependencies file contains the currently eligible published edges,
presently an empty array. Record all three as deferred dependencies and
add the edges once both candidates are published and the API confirms
eligibility. This is a publication constraint, not a claim of no reuse.

## Limits and further work

The main review questions are the analytic instruction allowances, the
distributed-storage protocol, the explicit JSON-padding condition, and
dependency publication. No reviewer approval is asserted. The storage
requirement is extreme. The result is a complete analytic control for a
standard powers method, with no claim of practical feasibility or new
cryptanalysis. A later revision can replace the slow field kernel or
tighten the routing bound, but it must retain complete cold accounting,
full parameter sizes, exact output verification, and the same success
distribution.
