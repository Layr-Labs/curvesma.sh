# Point-valued output from a verified affine intersection

For this track the recovered scalar is an internal intermediate only. After
the exact test g^t=z, return P=P0 and Q=[t]Q0 on the protected twist over Fp2.
For t=0 return Q=null; otherwise output affine x,y, each as [c0,c1], and each
Fp component as 64 lowercase hex digits. P0 is the protected G1 generator.
The JSON keys are exactly P and Q. Neither output is an unverified scalar.

Represent Fp2 as pairs with u^2=-1. Addition is componentwise; multiplication
uses (a+bu)(c+du)=(ac-bd)+(ad+bc)u. For a nonzero pair, inverse is
(a-bu)/(a^2+b^2). The denominator cannot vanish for nonzero pairs in the
protected Fp2 field. Compute its Fp inverse by a^(p-2), with bounded binary
exponentiation. This gives exact finite operations using the halfword kernel.

Use affine binary scalar multiplication starting at infinity, at most 256
doublings and 256 additions. Handle infinity by copying the other point;
same x/opposite y gives infinity, including doubling a point of y=0. Otherwise
use slope 3*x^2/(2*y) for doubling or (y2-y1)/(x2-x1) for distinct x, then
x3=slope^2-x1-x2 and y3=slope*(x1-x3)-y1. All arithmetic is in Fp2, and the
twist has zero x coefficient, so these are its exact group formulas. A
nonzero denominator is inverted only after exceptional cases are handled.

The loop invariant proves the result is [t]Q0; group closure ensures curve
membership and order-r subgroup membership. Explicitly check the curve
equation and [r]Q=infinity if desired; all those bounded checks are included
in the fixed budget. P0 and Q0 are public profile constants, not searched
advice. Each affine operation needs fewer than 2^16 base-field products/adds
including inversion, control and conversions, so <2^36 instructions. Fewer
than 1024 point operations and all remaining serialization/curve tests are
below 2^48, within implementation.md's 2^52 fixed allowance.

The protected reduced pairing is bilinear. Hence e(P0,[t]Q0)=g^t=z, by the
scalar check already performed. No uncharged full pairing computation is
claimed or needed to prove this equality. This is one admissible method,
not a claim that arbitrary pairing inversion has the same exact cost as GT DLP.
