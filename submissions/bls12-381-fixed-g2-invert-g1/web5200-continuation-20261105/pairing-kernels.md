# Exact protected pairing formulas, lowered to scalar 64-bit work

The formulas below are reproduced/adapted from the organizer appendix
submissions/bls12-381-invert-g1/pairing-kernels.md at curvesma.sh commit
6423c2728dc03d5f2d647c7f22adf1dac0643a66 and py_ecc 8.0.0 commit
04151f01f59f902ab932a51e0ca0ebce3883fc51. They compute the exact protected
reduced e(P,Q) with library argument order pairing(Q,P). Here P and Q
are arbitrary valid points in the respective prime-order subgroups;
neither is silently fixed to its generator. Source hashes were checked
against the protected profile; code was read as text and never executed.

Use the arithmetic.md implementation: six 64-bit limbs per Fp value,
A=2048, M=1599488, and E=2^30 as a common finite envelope for Fp12/Fp2
arithmetic and movement. Canonical Fp12 multiplication costs
144M+188A+65536=230776832<E. Fp12 additions use at most 12A+4096<E;
Fp2 multiplication uses four M and two A plus 4096, also below E.
Use padded k=384 in exponent-length bounds. Constants are not unit-cost
arbitrary-precision primitives. The appendix's old wide-word numeric
ledger is not reused; all calls below use these scalar bounds.

## Projective point operations over Fp2

Coordinates (x,y,z) represent (x/z,y/z); z=0 represents infinity.
The following are straight-line formulas in Fp2.

For doubling:

```text
W = 3*x*x
S = y*z
B = x*y*S
H = W*W - 8*B
S2 = S*S
X = 2*H*S
Y = W*(4*B-H) - 8*y*y*S2
Z = 8*S*S2
return (X,Y,Z)
```

This has 17 multiplications (constants included) and three subtractions.
For addition of (x1,y1,z1) and (x2,y2,z2):

```text
if z1=0: return (x2,y2,z2)
if z2=0: return (x1,y1,z1)
U1=y2*z1; U2=y1*z2; V1=x2*z1; V2=x1*z2
if V1=V2:
    if U1=U2: return double(x1,y1,z1)
    else: return (1,1,0)
U=U1-U2; V=V1-V2
Vsq=V*V; VsqV2=Vsq*V2; Vcube=V*Vsq; W=z1*z2
A0=U*U*W - Vcube - 2*VsqV2
X=V*A0
Y=U*(VsqV2-A0) - Vcube*U2
Z=Vcube*W
return (X,Y,Z)
```

The ordinary branch has 15 products and six subtractions. The doubling
branch has the four preliminary products plus the doubling kernel: 21
products and three subtractions. Charge at most 32 arithmetic calls for
any branch. Equality tests compare the twelve 64-bit limbs of an Fp2 value.
Argument/result copying, these tests and branch/loop control together are
bounded by 65536 word operations, less than E. Thus 64*E bounds either
point operation, including all exceptional cases, loads and stores.

## Twist and base-point embedding

The flat Fp12 basis satisfies w^12=2*w^6-2 and u=w^6-1, so u^2=-1.
For each Fp2 coordinate t=a+bu compute c=a-b. The twist is:

```text
nx=(ax-bx)*w + bx*w^7
ny=(ay-by)   + by*w^6
nz=(az-bz)*w^3 + bz*w^9
return (nx,ny,nz)
```

All omitted coefficients are zero. This is exactly the pinned `twist`
coefficient placement. It needs three A calls and construction/copying
of three 72-word Fp12 arrays. At most 65536 word operations cover all array
initialization, coefficient loads/stores and control; three A plus this
reserve is below 4*E. Charge 64*E. To cast a base-field projective P into
Fp12, put each coordinate in coefficient zero, setting coefficients
1..11 to zero. The same 65536 reserve bounds this one-time embedding.

## Rational line evaluation over Fp12

For projective points R1=(x1,y1,z1), R2=(x2,y2,z2) and evaluation point
T=(xt,yt,zt), compute numerator n and denominator d without division:

```text
mn=y2*z1-y1*z2
md=x2*z1-x1*z2
if md != 0:
    n=mn*(xt*z1-x1*zt) - md*(yt*z1-y1*zt)
    d=md*zt*z1
elif mn=0:
    mn=3*x1*x1
    md=2*y1*z1
    n=mn*(xt*z1-x1*zt) - md*(yt*z1-y1*zt)
    d=md*zt*z1
else:
    n=xt*z1-x1*zt
    d=z1*zt
return (n,d)
```

Counting the preliminary mn and md in every branch, the ordinary branch
needs 12 products and five subtractions, the tangent branch 16 products
and five subtractions, and the vertical branch seven products and three
subtractions. Fp12 zero tests inspect 72 limbs each. A 65536-word reserve
bounds coordinate movement, tests and control, less than E. Charge 64*E
for a complete line call. Multiplication by 2 or 3 is not free.

## Exact Miller schedule and wrapper

The positive loop parameter is 0xd201000000010000, 64 bits long. Initialize
R=Q and fnum=fden=1; compute twist(Q) and cast(P). Process bits 62..0
in descending order, exactly as the reference's `pseudo_binary_encoding[62::-1]`:

```text
(n,d)=line(twist(R),twist(R),cast(P))
fnum=fnum*fnum*n; fden=fden*fden*d
R=double(R); twistR=twist(R)
if bit=1:
    (n,d)=line(twistR,twist(Q),cast(P))
    fnum=fnum*n; fden=fden*d
    R=add(R,Q); twistR=twist(R)
```

There are 63 doubling substeps and at most 63 addition substeps. Each
substep costs at most:

| Component | Charged E units | Justification |
| --- | ---: | --- |
| Fp2 point operation | 64 | Counts and exceptional branches above |
| Fp12 line evaluation | 64 | At most 16 products, five subtractions and bounded movement |
| Twist | 64 | Explicit coefficient placement above |
| Accumulators | 16 | Four products for doubling or two for addition, plus copying |
| Schedule/control | 1 | At most 65536 word operations |
| Total | 209 | Strictly less than the declared 512-E allowance |

The two unused slots in 128*512*E cover initialization/wrapper and final
quotient multiplication. Curve validation uses y^2*z=x^3+b*z^3 in the
appropriate field, at most eight products and one addition per point,
with zero-z tests. This is the pinned project's `is_on_curve` condition.
Initialization, validation of both points, embedding, initial twist,
constant/accumulator construction, quotient multiplication and output
copying together fit 256*E plus 65536 word operations, below these two
slots. If P or Q is infinity the wrapper returns 1, as the pinned map does.

For nonidentity valid inputs the loop's Fp2 points are [n]Q with
0<n<=0xd201000000010000<r. They never reach infinity or a point of order
two. An addition step cannot have R=-Q: its positive scalar n is far
below r-1. Consequently the projective z values and tangent/ordinary line
denominators in this schedule are nonzero. With cast(P).z=1, fden is a
product of nonzero factors. The final inverse below is well defined.

Compute f=fnum*(fden^(p^12-2)), then return f^((p^12-1)/r).
These are the reference quotient and its exact final exponent; Fermat
inversion changes the implementation, not the pairing value. Each exponent
has at most 12*k bits. Binary powering takes at most 24*k products per
exponent, at most 48*k in total. The final quotient multiplication is
already in the two wrapper slots; no product is omitted. Scalar exponent
bits are scanned once; at most 64 operations per bit plus 4096 per power
cover index/control and accumulator initialization, below 2^28.

Construct p^12 and divide p^12-1 by r using scalar limb arithmetic,
not a big-integer oracle. At most 73 64-bit limbs hold the 4608-bit
padded numerator and carry. Twelve bit-serial products, each with at
most 384 multiplier bits and 128 instructions per limb per bit, cost
at most 12*384*73*128=43057152. Binary long division uses at most
4608 shift/compare/subtract stages, each at most 128*73 instructions,
also 43057152. Add 2^27 instructions for constant decoding, powering
bit control, copies and initialization. The sum is below 2^28.

Thus J=(128*512+48*384)*E+2^28<2^47 prices the exact reduced map.
The rational quotient and the full final exponent are included. No raw
Miller comparison, free pairing, precomputed field table or unstated
numeric cost premise is used. Input and output subgroup tests are
separately in arithmetic.md's fixed reserve.

## Retained upstream license


The MIT License (MIT)

Copyright (c) 2015 Vitalik Buterin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
