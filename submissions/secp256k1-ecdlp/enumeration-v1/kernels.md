# Pinned pairing and raw Miller formulas

Adapted from the original repository baseline pairing-kernels.md and py_ecc 8.0.0, commit 04151f01f59f902ab932a51e0ca0ebce3883fc51. MIT notice in py-ecc-LICENSE.txt. All costs are derived anew in proof.md. Coordinates are (x/z,y/z); z=0 is infinity. Extension elements use the target profile flat polynomial basis. Arithmetic is reduced after each operation.

## Projective doubling over Fp2

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

## Projective addition over Fp2

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

## Twist (t=a+b*u, u=w^6-1)

```text
nx=(ax-bx)*w + bx*w^7
ny=(ay-by)   + by*w^6
nz=(az-bz)*w^3 + bz*w^9
return (nx,ny,nz)
```

## Rational line, returning numerator and denominator

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

## Miller substep

```text
(n,d)=line(twist(R),twist(R),cast(P))
fnum=fnum*fnum*n; fden=fden*fden*d
R=double(R); twistR=twist(R)
if bit=1:
    (n,d)=line(twistR,twist(Q0),cast(P))
    fnum=fnum*n; fden=fden*d
    R=add(R,Q0); twistR=twist(R)
```

## Wrapper and normalization

For input infinity return Fp12 one. Otherwise take the AFFINE G1 input (x,y) as (x,y,1), and the fixed affine Q0 as (xQ,yQ,1). Cast each G1 coordinate to coefficient zero in Fp12. Start R=Q0 and fnum=fden=1. In the Miller substep block above, process bits 62 down to 0 of the positive integer 0xd201000000010000; the leading bit 63 is represented by initial R=Q0. No signed recoding or conjugation is applied. Return fnum * fden^(p^12-2) for raw M. For reduced e, raise that quotient to (p^12-1)/r. This is the profile's pinned quotient, including all denominators.

Throughout the schedule the Fp2 points are positive multiples [j]Q0 with 0<j<=0xd201000000010000<r. They are neither infinity nor order two; an addition never encounters R=-Q0. With cast(P).z=1, the line denominators in this schedule are nonzero. Fermat inversion is well defined. Point curve checks use y^2*z=x^3+b*z^3, with b=4 in G1 and b=4(1+u) in G2. The profile inputs are honestly generated; enumeration itself preserves membership in the prime-order subgroup.
