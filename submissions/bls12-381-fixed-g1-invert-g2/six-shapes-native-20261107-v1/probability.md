# Deterministic singular-denominator certificate

K=ceil2061R/2000 and all symbols are as in proof.md. The full-order public inequalities below are checked at this exact cap.

## Complete finite Frobenius interval lemma

Let U=0xd201000000010000,the positive magnitude of protected parameter_x. Full public profile arithmetic verifies

 r=U^4-U^2+1, xi=p modr=-U modr.

Let M=K+1 and A=floor(M/U). All the following exact integer inequalities hold for the full protected r and chosen K:

 U^2-1>M,
 r-U^2*M>M,
 r-A*(U^2-1)-(U-1)*U^3>M,
 r-(U-1)*(U^3-U)>M,
 U^3-U-A>M,
 r-A>M,
 U-1>T=128.

These are finite public certificates,not asymptotic size arguments. Also M<U^2/3. We prove that for integers1<=s<t<=M and ell0..11,

 t=xi^ell*s modr implies ell=7 and t=U*s.

Write s=aU+b,0<=a<=A,0<=b<U. Reduced power representatives are

|ell|xi^ell modr represented as an integer polynomial in U|
|---|---|
|0|1|
|1|-U|
|2|U^2|
|3|-U^3|
|4|U^2-1|
|5|-U^3+U|
|6|-1|
|7|U|
|8|-U^2|
|9|U^3|
|10|1-U^2|
|11|U^3-U|

They follow by repeated multiplication using U^4=U^2-1 modulo r and xi=-U. For ell1/6/8/10,the negative product -c*s has0<c*s<=U^2*M<r-M;its canonical residue r-c*s exceeds M. For ell2/4,the positive product lies between U^2-1>M and U^2*M<r. Ell0 gives t=s,excluded. Ell7 gives exactly t=U*s without wrapping because U*M<r.

For ell3/9,reduce using s=aU+b. The products are respectively -H and H,where H=a(U^2-1)+bU^3. For s>0,H>=U^2-1>M,and H<=A(U^2-1)+(U-1)U^3<r-M. Thus both canonical residues exceed M.

For ell5/11,the reduced products are respectively a-b(U^3-U) and b(U^3-U)-a. If b=0,the former is a<s=aU and the latter is -a with canonical residue r-a>M. If b>=1,set J=b(U^3-U)-a. Its lower bound U^3-U-A exceeds M;its upper bound(U-1)(U^3-U) is<r-M. Thus both J and -J have residues>M. This covers every exponent and every s;only ell7 remains. No unproved lattice,short-relation,complex-norm or random-alpha assumption is used.

## Singular denominator occurs only on very few long pairs

A true raw orbit match X_t=sigma^ell X_s obeys

 h^[(t+1)-xi^ell*(s+1)] * z^[B_t-xi^ell*B_s]=1.

If D=B_t-xi^ell*B_s were0,then h's prime order forces t+1=xi^ell*(s+1) modulo r. Applying the interval lemma at1<=s+1<t+1<=M leaves only ell7,t+1=U(s+1). There are at most A=floor(M/U) such chronological pairs. Every such lag is (U-1)(s+1)>=U-1>T.

On a good generator table,the pointwise conditional kernel at these lags is at most(1+eta)/r for any specified raw target,eta=2^-128. In particular it bounds the fixed twist7 target sigma7(X_s),conditional on X_s. Union over the at most A exceptional time pairs gives

 Pr[any singular-denominator orbit match]<=A*(1+eta)/r<2^-180.

This bounds the first actual terminal match by extending an ideal iid-index walk through K only for analysis. No adaptive conditioning,independent pair assumption or uniform-denominator claim is needed. The exact last inequality uses the full public order. All other genuine nonidentity orbit matches have D!=0 deterministically. Identity visits and fingerprint/dictionary/coupling failures are separately charged.

The dictionary returns matched logical p and native payload. Since it allocates initial index1 and exactly one index per earlier transition,the saved time is p-1 and A_s=p;current A_t=t+1 is the preserved accepted counter. The old raw/payload A slots are unused initialized zero padding and must not be read. The unique actual raw twist is recovered by reconstructing the saved raw point from these exact coefficients and stored B and testing ks-kc+4m,m0..2. This returns the group invariant and never chooses an unchecked scalar-dependent alternative. With D!=0,compute a=(xi^ell*A_s-A_t)/D,verify h^a=z,form canonical Q=[a]Q0,check curve/subgroup and exact pairing(Q,P)=z before output. All reconstruction,Frobenius constants,three tests,inversion,scalar multiplication,verification and output are bounded cold paid. Scalar-only or unchecked output is never accepted.

## Complete failure union

Raw generator/initial state stationarity is proved in success.md;hash.md bounds NH coupling2^-55 and scalar cap2^-970. Good-table failure2^-66,identity visits(K+1)/r<2^-100,storage.md dictionary cap2^-100,and the singular bound2^-180 above have a union bound,without independence. success.md gives full finite repeat>0.3900031430168671;subtracting all these gives success>0.39000314301686>0.39. No lazy/parity/rare-zero cost cap is used. Standard full BLAKE3 H and explicit finite banked storage remain material. No official score/qualification is claimed.
