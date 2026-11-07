# Exact production reduction and canonical b-half key

## Exact production reduction and constants

Use the protected full BLS12-381 p,r,Q0 andpinned pairing. Validate/convert canonical P,z,including nonidentity/subgroup checks. If z=1 return canonical G2 infinity. Otherwise compute h=pairing(Q0,P) in the library's reversed argument order. Nondegeneracy andprime order make h a GT generator for every allowed P. There is a unique a withz=h^a;neither a,nor challenge coins,nor powers are supplied. Final output is Q=[a]Q0 with independent h^a=z andpairing(Q,P)=z checks,not a scalar orreplacement input point.

Let sigma(X)=X^p on GT andxi=p mod r. At this protected profile xi has order12 andPhi12(xi)=xi^4-xi^2+1=0 mod r. p=7 mod12,r=1 mod12. Every nonidentity GT element has a free12-element sigma orbit. Work in the standard native Montgomery tower Fp2=Fp[u]/(u²+1),Fp6=Fp2[v]/(v³-1-u),Fp12=Fp6[w]/(w²-v),matching the public polynomial model through w^6=1+u. Conversion is paid cold. Write X=a+bw,with a=a0+a1v+a2v²,b=b0+b1v+b2v²,each coefficient in Fp2.

For k=0..11,d=0..5 construct cold c[k,d]=(1+u)^(d*(p^k-1)/6) in Fp2. Then sigma^k(x*w^d)=c[k,d]*conj^k(x)*w^d. Native constants andtheir sums for Karatsuba are cold constructed,not free advice. The public degree/order identities can be independently checked exactly.

## Unique canonical orbit key without twelve full point transforms

Profile-specific exact arithmetic refinement: p=19 mod24. Since (1+u)^4=-4 is real, c[k,d] is real,purely imaginary,orhas equal/opposite real/imaginary coefficients. For even k the exponent d(p^k-1)/6 is0 mod4. For odd k its residues at d=0..5 are0,3,2,1,0,3. Thus every needed constant product can use at most two Fp multiplications,not three. For c=(s,0) or(0,s),scale/swap/negate two coordinates;for c=(s,s),use s(x-y),s(x+y);for c=(s,-s),use s(x+y),s(y-x). Conjugate the input when required. These are exact field identities,with additions/negations charged below. No new assumption orhash output change is introduced.

If b0=b2=0,then vb²=b1²v³ lies in Fp2 and a²=1+vb² lies in Fp2. Since a belongs to the odd-degree extension Fp6/Fp2,its degree divides both2 and3,so a is in Fp2. Hence X is in Fp2[w³]=Fp4. Its intersection with prime-order GT is identity because gcd(r,p^4-1)=1. Thus b0=b2=0 characterizes X=1 on GT. Except z=1 handled above,abort on such a visited state;probability.md charges this rare event.

Choose x=b0 if b0!=0,otherwise x=b2,with d=1 or5 respectively. Let T(x)=c[1,d]*conj(x),beta=c[2,d]. beta is a primitive sixth root in Fp,and T²(x)=beta*x. Its order is6 because w has degree6 over Fp2 and d is coprime to6. For x!=0 the12 T-images are distinct:an even stabilizer requires beta^m=1;an odd stabilizer would give beta^k=1 for odd k,impossible. Define y=T(x). The12 images are beta^m*x (applied exponent k=2m) and beta^m*y (k=2m+1),m=0..5.

Compute y by one constant Fp2 multiplication,then beta*x and beta*y byfour Fp multiplications. For either z=x ory,generate its six scalar rotations using beta²=beta-1:
 z,beta*z,beta*z-z,-z,-beta*z,z-beta*z.
Choose the unique minimum of these12 Fp2 values,comparing real then imaginary canonical native Montgomery residues as unsigned integers. Let k be its applied exponent;C=sigma^k X is the canonical point. This choice is invariant under starting orientation,so k(sigma X)=k(X)-1 mod12.

To build the exact key of C, reuse the selected b pair and transform the other two b pairs by c[k,d]*conj^k(b_i), at most four further Fp multiplications. The b half alone is injective on GT. If two norm-one points X=a+bw,Y=c+bw share b, then a²=c². Distinct points require c=-a andY=-X^-1, soXY=-1 would belong toGT. This is impossible in prime odd-order GT. Validated subgroup membership is preserved by every jump. The claim is restricted to GT, not the full even-order norm-one group.

Thus no a slot, sign, square root or inverse is needed for canonical-key generation. Store36 b words and a constant-zero37th word, with applied k separate. The twelve selected-pair images remain distinct, so canonical orientation and equivariance are unchanged. Full37-word equality is redundant but exact. Cold setup initializes word36 in both working records; key-copy writes only36 words and reference/scalar/prefix/k writes never touch it. Saved copies retain zero.


Current working records retain rawX at0..575,A/B at576..639,canonicalb words640..927,constantzero928..935 andappliedk at944. U andlabel/contextbuffers are separate shared data. Saved tree payloads are512 bytes with the layout in tree.md;they do not store full rawpoints. Canonical b equality plus saved/current k is sufficient to establish the raw Frobenius twist. No a/sign recovery is required.

The current implementation uses the explicit cheap-class branches andhardweighted counter in canonical-budget.md. It still computes exactly the same canonicalb/key/k;copies andnegations replace only multiplicationby+1/-1. Worst-case7696 remainsvalid,andthehardtotalcapimprovesaggregatework.
