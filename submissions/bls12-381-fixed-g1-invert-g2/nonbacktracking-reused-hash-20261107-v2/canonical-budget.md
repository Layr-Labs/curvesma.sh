# Exact orientation-dependent canonical cost and a hard total cap

This refines the same canonical key without changing its value. There is no additional heuristic. The full worst-case helper bound remains7696, but cheaper Frobenius classes are implemented explicitly and a hard counter bounds their total work. Any budget exhaustion is an algorithmic failure whose probability is proved below.

## Cheap classes and exact upper cost

The selected pair is b0 or b2, so the other two pairs always include b1. The first stage still computes y=T(x),beta*x,beta*y andall twelve candidates,using six Fp multiplications. The minimum determines k before the other two b pairs are transformed.

In the native degree-six extension over Fp2, w³ has degree two: its square lies in Fp2 and it cannot itself lie there without reducing the degree of w below six. Thus w³ generates Fp4 and sigma²(w³)=-w³. Consequently sigma^(2m) acts on the b1*w³ term by(-1)^m, with no Fp multiplication. For k=0, all b pairs are copied. For k=6, the quadratic Fp12/Fp6 involution negates every b pair. At other even k, only the remaining b0/b2 pair may need two Fp multiplications. Odd k may need four for both remaining pairs.

Copies and canonical negations use at most the existing four48-instruction field-operation allowances per structured pair. They do not introduce an uncharged copy. All rotation,comparison,key/control andwrapper allowances are retained. Define cost class c(k)=0 for k=0,6; c(k)=1 for other even k; c(k)=2 for odd k. The helper is bounded by

 f(k)=6032+832*c(k).

Here6032=6*416+3536, where3536 is the unchanged nonmultiplication allowance. Each class unit adds two Fp multiplications,832. Thus f is6032,6864,7696 respectively. This is an upper instruction count for an explicit copy/negate/general-product branch, not a new field-operation price.

## Hard counter and bounded work on every path

After initial cold canonicalization,initialize a128-bit remaining class budget B=floor(75K/52). It fits because B<2K<2^127. After finding each new k andbefore remaining b transforms,subtract c(k),failing on underflow. LDP,SUBS,SBCS,STP andborrow branch,with class-table lookup andinterfaces,fit16 extra instructions per transition. Cold computation of B can use ordinary fixed-limb arithmetic;it is within2^120.

For any completed prefix,extra class units consumed are at most B. Common work is at most6032K. Hence total canonical work is at most6032K+832B<=7232K. A partially executed aborting helper is covered by one fixed7696 allowance in cold2^120;placing the check before optional products makes this conservative. This is a hard work bound,not just an expectation. Other algorithm caps still apply.

## Failure probability from proven stationary mixing

Applied k runs through all twelve values once on every nonidentity orbit, since k(sigma X)=k(X)-1. Uniform raw GT therefore gives each class frequency exactly from the twelve values. Extend g(X)=c(k(X)) by g(1)=2. Its stationary mean is

 mu=4/3+2/(3r), 0<=g<=2.

This is algebraic orbit counting, not an assumption about sampled orientations. The ideal nonbacktracking chain starts stationary. For lag at least T=32,the proven full-state relative kernel error epsilon=2^-127 implies

 |E[g(X_j)|state_i]-mu|<=2epsilon,
 |Cov(g(X_i),g(X_j))|<=4epsilon.

For shorter gaps use absolute covariance<=4. Thus S=sum_(t=1)^K g(X_t) obeys

 Var(S)<=4(2T+1)K+4epsilon K²=260K+4epsilon K².

At the protected K,r, B/K-mu>1/10, including floor rounding. Chebyshev gives

 Pr[S>B] <=26000/K+400epsilon<2^-100.

The bound holds for every good generator table because stationarity and mixing do. The actual algorithm follows the same ideal trajectory until a prior stop or this cap; g is nonnegative, so exhaustion on a prefix implies S>B on the full ideal extension. No conditioning on already adaptive fingerprints or successful paths is used. Fingerprint,phase,identity andtable failures remain separately subtracted. Budget decisions use group k only,not coefficients,so denominator uniformity is preserved.

## Updated complete bound

Canonicalbudget7232 isretained unchanged. Thecurrentdictionary-specificledgeris inscalar-audit.md:7232canonical+807tabulation+1594hash+46scalar+848controller+192head+224payload+32cappedtraversals+1retry=10976,retainedas11008. Cold2^130 includeseverybucket-headzero store. Thus W<=(31008K+2^130)/20000. Theexactinequality(31008K+2^130)^1000<20000^1000*2^126314 supports126.314. Addthecanonical-budgetfailure2^-100 anddictionary-capfailure2^-100 separately;successremains>0.39. Thecanonicalfieldvalues andgroup/hashtrajectoryareunchangeduntilanexplicitabort. Peakmemory<2^136bytes;H andfinitebanked-memoryadmissibilityremainmaterial.
