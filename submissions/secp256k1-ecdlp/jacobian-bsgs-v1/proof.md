# secp256k1-ecdlp: Jacobian-coordinate deterministic BSGS with one batch inversion

Status: ready for independent review; no qualification or emitted score.
This is a generic accounting construction, not a best-known attack claim.
The normative inputs/outputs are those of `production-v2`; `claim.json`
binds that profile and `classical-group-work-v2` by SHA-256. No restricted
challenge distribution, advice, random-oracle premise or quantum operation
is used.

Parent: the affine-coordinate control submitted at `129.000802`. The
algorithm, the success argument, the coverage, the table layout, the
sorting and the lookup are unchanged. What changes is the coordinate
system used for the group steps: the affine formulas charged a base-field
inversion per step, the projective formulas below charge none, and the
single inversion the table needs is paid once through Montgomery's batch
trick.

## Construction and success

Let p and r be the protected field modulus and subgroup order. Set
b=bit_length(r)=256, L=floor(b/2)=128, m=2^L and h=ceil(r/m).
All table indices/addresses and i*m+j fit 256 bits for these fixed targets.

```text
Use P from the protected profile; input Q. Starting at infinity, store (X,Y,Z,j)
in Jacobian coordinates for jP, j=0..m-1, by repeated mixed addition of the
affine P; retain mP after the final update. Batch-normalise every record with
Montgomery's trick: prefix-product the m Z-coordinates, invert the product
once, and walk the prefix products backwards, so each record holds affine (x,y,j).
Set D=-(mP) in affine form and gamma=Q in Jacobian form. Sort the table by
(infinity_tag,x,y), tag 1 for infinity and 0 otherwise. For i=0..h-1,
binary-search gamma against the affine table, testing X=x*Z^2 and Y=y*Z^3
modulo p (and Z=0 against the infinity record); on a match with j, set
a=(i*m+j) mod r, verify [a]P=Q, and return a. Otherwise set gamma=gamma+D
by mixed addition. Reject malformed inputs; exhaust the fixed loops with failure.
```

For every valid Q there is a unique x in [0,r) with Q=[x]P. Write
x=i*m+j with 0<=j<m and 0<=i<h. At giant step i, gamma=[j]P; the
complete table therefore contains a match. The comparison above is an
exact projective equality test, so a match is a match of points and not
of representatives. Any match yields [a]P=Q after reduction modulo r. The
independent attack-side equality check accepts it. Infinity gives x=0 and
is represented explicitly; a giant step that reaches infinity has Z=0 and
is compared against the infinity record only.

The algorithm succeeds for every valid input: probability 1, hence >=0.39
over the uniform challenge distribution. No randomness, retries or
amplification are needed. A finite worst-case work budget is used; it is
not an expected-runtime or group-operation exponent standing in for total work.

## Word-operation ledger

Integers are little-endian arrays of s=ceil(k/256) limbs, where
k=bit_length(p)=256 and s=1. Each listed cost-model primitive costs one.
Operands of a field operation are canonical integers in [0,p).

| Kernel | Upper bound in word operations |
| --- | ---: |
| Modular addition/subtraction | A=128*(s+1)=256 |
| Base-field multiplication | M=k*(2*A+64)+128*(s+1)=147712 |
| Nonzero base-field inverse | I=2*k*(M+64)+128*(s+1)=75661568 |
| Affine point addition/doubling (reference price) | C=I+16*M+32*A+4096=78037248 |

For A, zero-extend to s+1 limbs, add/subtract with carry/borrow, compare
against p, and make one conditional correction. Each limb pass needs at
most 16 operations for loads/stores, arithmetic and carry/borrow handling;
allow six passes plus 32*(s+1) for comparison, selection and control.
For M, scan the k bits of the multiplier; conditionally add the running
multiplicand to the accumulator and double it modulo p each iteration.
Two A calls and 64 bit-extraction/loop operations per bit suffice; reserve
128*(s+1) for initialization/output. Inversion uses binary powering by p-2:
at most 2*k multiplications plus 64 operations per bit. All instructions,
including memory transfers, are counted; multiplication/inversion are not
unit-cost primitives. C is retained only as the organizer's fixed reference
price; no step of this construction is charged C.

### Projective kernels

Jacobian doubling for a=0 (the `dbl-2009-l` schedule) uses
A=X1^2, B=Y1^2, C=B^2, D=2*((X1+B)^2-A-C), E=3*A, F=E^2,
X3=F-2*D, Y3=E*(D-X3)-8*C, Z3=2*Y1*Z1: five squarings, two
multiplications and at most nine additions, since multiplication by
2, 3 and 8 is repeated addition. Jacobian mixed addition with an affine
second operand (the `madd-2007-bl` schedule) uses
Z1Z1=Z1^2, U2=X2*Z1Z1, S2=Y2*Z1*Z1Z1, H=U2-X1, HH=H^2, I=4*HH,
J=H*I, r=2*(S2-Y1), V=X1*I, X3=r^2-J-2*V,
Y3=r*(V-X3)-2*Y1*J, Z3=(Z1+H)^2-Z1Z1-HH: four squarings, seven
multiplications and at most eight additions. Both schedules are covered by

| Kernel | Upper bound in word operations |
| --- | ---: |
| Jacobian doubling / mixed addition | J=12*M+16*A+4096=1780736 |
| Jacobian-vs-affine equality probe | Q=6*M+4*A+4096=891392 |
| Montgomery batch normalisation, per point | BN=3*M+64=443200 |

J allows twelve multiplication-class calls where the two schedules need at
most nine, sixteen additions where they need at most nine, and 4096 for
loads, stores, index arithmetic and branches. Q covers
Z^2, Z^3, x*Z^2 and y*Z^3 (four multiplication-class calls), the two
comparisons, the Z=0 test and control, with margin. BN covers one
multiplication for the forward prefix product, two for the backward pass
and the per-point copy/selection; the single inversion of the batch is
charged separately as I, exactly once.

Point addition uses affine formulas only for the reference price C and for
the affine operands P, mP and the normalised table. The stored table holds
(X,Y,Z,index) in W=4 words while the walk is in progress; after batch
normalisation Z=1 and the same record holds affine (x,y,index). No
per-step inversion is performed anywhere in the baby build or the giant
walk.

Setup reserves 2^24 operations for decoding/constants and the step negation.
Verification/output reserves 4*b*C+2^24: this covers [r]Q and [a]P checks,
curve/encoding checks, scalar reconstruction by shift/add and serialization.
That reserve is deliberately kept at the reference price so the claimed
total does not depend on a cheaper verification schedule.

Use bottom-up merge sort with two m-record buffers and exactly L passes.
A lexicographic key comparison costs B=16*K+64; moving a record costs
D=16*W+64. These bounds allow all operand loads, stores, index arithmetic,
comparisons and branches, up to 16 operations per word plus 64 fixed work.
Each merge pass has at most m comparisons and 2*m record moves, with
128*m additional loop work. Thus sorting costs
S=m*L*(B+2*D+128)+4096*L. Binary search performs at most L+1 key
comparisons; V=(L+1)*(B+64)+D+256 covers lookup and the returned record.
No constant-time hashing/table oracle is assumed.

For group-step cost, charge m*(J+BN+D) to build and normalise the baby
table plus the single I, S to sort it, and h*(J+Q+V) for the full giant
search. Setup and verification/output are added in full. The accounting
below gives the exact integer components and their sum T. Preprocessing
includes the whole baby table, its normalisation and its sorting, with no
reuse credit.

Peak memory is 2*m*W*32+2^20 bytes. The two arrays include all keys,
projective coordinates and indices; the backward Montgomery pass rewrites
the Z slots in place, so no third array exists. The 1 MiB reserve includes
code, public constants, fixed field/point scratch, loop state and stacks.
No table is uncharged advice. Input/output and all record lifetimes are
included. This memory is astronomical and mandatory for this construction;
it is not a practical attack.

## Evidence and comparison

The correctness argument and cost calculation are analytic; no full-size
execution is claimed. Organizer arithmetic fixtures check the exact
point/GT convention; organizer exhaustive finite-group checks exercise
table coverage, last-block handling and modular output recovery. They do
not certify production cost, scientific approval or an executed production
break. Heuristic premises: none.

Pollard rho and automorphism-aware methods must be compared under this same
cost model before selecting a best-known ECDLP starting score. This package
does not implement those optimizations. Beating the affine control alone
does not establish novel cryptanalysis or a weakened hardness assumption;
it establishes that the inversion-bearing affine step was not the cheapest
correct accounting of the same deterministic search.

References: [generic DLP algorithms, HAC section 3.6](https://cacr.uwaterloo.ca/hac/about/chap3.pdf);
[Guide to Elliptic Curve Cryptography, Jacobian coordinates](https://link.springer.com/book/10.1007/b97644);
[automorphism-aware Pollard rho](https://www.iacr.org/archive/pkc2014/83830154/83830154.pdf).
The needed construction and accounting are given above; review need not
fetch citations to reconstruct them.

## Submitted computation bound

The following integer word-operation bounds follow from the implementation
formulas above. They are ordinary proof evidence, not a required submission
format.

| Component | Word operations |
| --- | ---: |
| setup | 16777216 |
| baby table and normalisation | 2^128*2224064 + 75661568 |
| merge sort | 2^128*63488 + 524288 |
| giant search | (2^128-1)*2695216 |
| verification and output | 79926919168 |

Summing the coefficients: 2224064+2695216+63488 = 4982768, and the
constants are -2695216+75661568+524288+16777216+79926919168 = 80017187024.
Hence

```text
T <= 2^128*4982768 + 80017187024 < 2^124.1 * 78037248.
```

Dividing the exact integer T by the organizer reference price 78037248
gives submitted time_log2 **124.1**, rounded upward. Preprocessing is
`2^128*2287552 + 76185856 < 2^122 * 78037248`, so preprocessing_log2 is
**122.0**, rounded upward, and is included in total work.
Memory_log2_bytes remains **136.000001**. All table access, normalisation,
sorting, failure paths, verification and output costs are included;
reference-operation normalization changes units, not the work or algorithm.
