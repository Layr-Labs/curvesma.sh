# Exact three-element trace fibers at full production size

Let sigma(X)=X^p in protected prime-order GT of embedding degree12. For nonidentity X,the roots X,sigma4(X),sigma8(X) are distinct because xi=p modr has order12. Their product is1:1+p^4+p^8=0 modr. If c is their sum,the second elementary symmetric function is c^(p^2):inversion is sigma6 and{6,10,2}={0,4,8}+2 mod12. Their exact monic polynomial is therefore
T^3-c*T^2+c^(p^2)*T-1.
Equal traces identify precisely these three roots,with no additional GT fibers. This is the XTR trace mechanism;the needed polynomial argument is supplied here. Material external source:Lenstra/Verheul,The XTR Public Key System,CRYPTO2000,Lemma2.2.1,https://www.cs.ru.nl/E.Verheul/papers/crypto2000/crypto2000.pdf . No source performance/qualification is inherited.

Use the pinned native tower Fp2=Fp[u]/(u^2+1),Fp6=Fp2[v]/(v^3-1-u),Fp12=Fp6[w]/(w^2-v). Write raw X=a0+a1*w^2+a2*w^4+b0*w+b1*w^3+b2*w^5. sigma4 fixes w^3 and multiplies w by a nontrivial cube root;its trace kills powers not divisible by3. Thus c=3(a0+b1*w^3). Since3 is invertible,use S=c/3 without computing a divide:copy native a0 at0..95 and b1 at384..479. These four full Fp coefficients give192 canonical native Montgomery bytes,an injective encoding of S in Fp4. Native canonicalization within the reference product is retained;no truncated field or free arithmetic primitive is used.

No four-image minimum,full-twelve-orbit quotient,orientation bit or field multiplication is used. This key intentionally identifies only the three-element automorphism orbit on the FULL order-r group. Equality of its canonical native words is equivalent to trace equality. The mathematical known coefficient A=1 in X=h*z^beta is unrelated to raw a0;unused scalar-A metadata is never read.

b1=0 iff X=1 in GT. If b1=0,then c=3a0 is in Fp2 and the polynomial factors(T-1)(T^2+(1-c)T+1). A nonidentity root would be in Fp4,impossible for an order-r point at embedding degree12. Thus testing all twelve b1 limbs detects identity. The key operation needs12 LDP/12 STP=24 for the two96-byte blocks,at most18 for the complete B-zero reduction/branch,and fewer than32 for pointers/frame/status/RET;retain128 instructions including every argument/address/output/interface. No arithmetic,exponentiation or Frobenius is performed per key.

On equal S and different canonical beta,raw points differ because z generates GT. Their relation is therefore sigma4 or sigma8. Try both public scalar equations as proof.md specifies,skipping a zero wrong denominator and verifying h^a=z. The true relation has nonzero denominator by prime-order h. Both fixed attempts and canonical G2/pairing verification are paid cold. No stored orientation or raw-point reconstruction is needed.

Pinned tower definition:https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fp12_tower.c . The same full protected p/r and native basis are used throughout.
