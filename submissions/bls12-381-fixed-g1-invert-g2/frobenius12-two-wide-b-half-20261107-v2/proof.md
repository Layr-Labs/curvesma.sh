# Frobenius12 quotient walk with two wide jump components

Proposed claim: cold classical work below2^(127033/1000)=2^127.033 protected dense-Fp12 reference equivalents,success>=0.39 conditional on the disclosed standard BLAKE3 idealization H. Peak storage<2^87 bytes,zero advice;preprocessing andall failures are included. This is analytic exploratory evidence,not an unconditional fixed-hash theorem,production execution,practical break,official score orpromotion. C remains20000. The banked-storage realization andits charged interface are explicit in storage.md;it is not installed orallocated on the verification VM.

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

A complete mutable/checkpoint record is960 bytes:576 raw X,64 scalar A/B,288 b-key bytes plus8 constantzero,8 prefix,8 applied k,8 padding. Key starts640,constantzero928,prefix936,k944,padding952. Padding is cold initialized;all other bytes are written by displayed arithmetic/canonicalization/preparation. Key equality excludes prefix,k,padding.

## Cold tables, standard hashes and labels

Set w=2^72 and use two independent component tables b=0,1. Tags are the exact16-byte ASCII strings FAPI-RHO-INIT-v1,FAPI-RHO-SCAL-v1,FAPI-RHO-TABL-v1,FAPI-RHO-JUMP-v1,FAPI-RHO-PHAS-v1. Seed is standard BLAKE3-256 of INIT||P_bytes||z_bytes. P_bytes is x||y,each48-byte big-endian canonical Fp coordinate;z_bytes is public polynomial coefficients0..11,each48-byte big-endian canonical Fp value. Thus INIT has688 bytes as in the source. A global unsigned128-bit counter starts0 and encodes as16 big-endian bytes. A scalar candidate is the low255 bits of B3(SCAL||seed||counter_16), interpreted big-endian; accept iff<r, with1024 attempts per scalar and failure otherwise. SCAL now has64 bytes, so this cold single-block standard hash uses len64,flags11,counter0. This explicitly replaces the old eight-byte counter, which would overflow at the new table width.

Generate(u[b,j],v[b,j]) for b=0,1,j=0..w-1 and alpha0, at most4w+1<2^75 accepted scalars and<2^85 attempts. The128-bit counter cannot wrap. At protected r>7*2^252, rejection probability is<1/8, so all capped failures have union probability<2^75*2^-3072<2^-970. Accepted values remain independent uniform inside H.

Construct every M[b,j]=h^u*z^v and its12 variants sigma^e M, with coefficients xi^e*u,xi^e*v modulo r. Each record is640 bytes:576 native GT andtwo32-byte canonical residues. Every power,variant,hash,write andbank selection is charged cold. No table is supplied as advice or reused over targets.

Retain the same73-position simple tabulation on canonical orbit b words and a constant-zero final digit:72 tables of2^32 plus one table of2 entries,each40 bytes. TABLE messages are54 bytes, with the same one-byte position,four-byte big-endian digit andone-byte part0/1. For digit i=0..71,use the low32-bit half of exact-key word floor(i/2) when i is even andits high32-bit half when odd;digit72 is0. Entry words0..3 are four big-endian numerical64-bit words of the part0 digest;word4 is the first big-endian64-bit word of part1. Store those five native words. These supply320 independent modeled bits. U is T72[0] XOR the72 selected b-digit entries, serialized asfive big-endian64-bit words. Every entry is constructed.

The two online standard hashes are J=B3(JUMP||seed||zero16||U) and T=B3(PHAS||seed||zero16||U). Their64-byte prefixes cache a paid chaining value with len64,flags1; final40-byte blocks use len40,flags10,counter0. All seven standard rounds and full required outputs remain. Let D_i be little-endian64-bit digest words. Define j0=D0(J)+2^64*(D1(J)&255),j1=D2(J)+2^64*(D3(J)&255). These are two independent uniform72-bit indices under H; unused bits are ignored. For V=sum_i2^(64i)D_i(T), set t=floor(12V/2^256). The phase differs from uniform12 by total variation<=12/2^256, separately charged over all labels. Four-limb multiplication by12 computes it without a wide unit-cost operation.

H is the disclosed unproved bridge treating these actual full standard domain-separated BLAKE3 computations as random functions. It includes the new64-byte SCAL messages and72-bit index extraction. No free hash oracle, reduced-round hash, independent random mapping, collision independence or extra relation heuristic is supplied. INIT remains688 bytes with eleven standard blocks; TABLE54,SCAL64 use flags11; all generic cold hash work is paid. The changed counter and index extraction change actual jump values and paths, not the standard primitive.

## Equivariant transition, checkpoints and exact output

Initialize(X,A,B)=(h^alpha0,alpha0,0), prepare canonical C=sigma^k X, its key,prefix andcached(j0,j1,t),and save it. For each transition set e=t-k mod12. Multiply raw X by M[0,j0,e] andM[1,j1,e],adding each rotated coefficient pair to A,B modulo r. There are exactly two reference dense Fp12 multiplications per transition. Intermediate factors are not canonicalized. After both factors recompute canonical key,k,prefix andlabels andsearch checkpoints.

Since canonical orbit labels agree at X andsigma X while k decreases by1, e increases by1 and f(sigma X)=sigma f(X). Thus the quotient transition is deterministic on group orbits and independent of coefficient histories. Abort on identity visits,whose loss is included.

Let R be the least integer with12R²>=r,L=ceil(1031R/1000),s=ceil(L/1022),K=L+s. Prefix lower_bound followed byfull37-word equality scans every equal-prefix record. The global confirmation cap is floor(K/256),decremented before each full comparison;abort before exceeding it. No prefix-only match is accepted. Save after unsuccessful lookup at s,2s,...,1022s;there are at most1023 records and1023s>K. Copies/inserts are paid cold.

At a key match set ell=k_saved-k_current mod12,so X_current=sigma^ell X_saved. Let D=B_current-xi^ell*B_saved. If D=0 fail;otherwise a=(xi^ell*A_saved-A_current)/D. Verify h^a=z,compute canonical Q=[a]Q0,check curve/subgroup andpairing(Q,P)=z,andoutput. Failed finishes are charged.

For quotient tail mu andcycle lambda,first repeat tau<=L,the saved cycle index b=ceil(mu/s)*s satisfies mu<=b<=mu+s-1 andb<=1022s. Equivariance gives equal orbit keys at b andb+lambda bytau+s-1<K. Corresponding raw twists rotate each later coefficient increment,so the formal relation is unchanged when both endpoints advance. The first detected checkpoint relation therefore has the same formal vector as the first orbit repeat,even across repeated raw rotations. probability.md bounds fruitless relations andproves denominator uniformity for the remaining cases.

## Complete computation and memory

Charge at most11008 ordinary scalar instructions per transition, including the two-factor controller and the reduced canonical helper. The explicit current list is in scalar-audit.md;this is a conservative upper count,not a changed reference price. Exactly two reference multiplications give

 W <=2K+(11008K+2^120)/20000=(51008K+2^120)/20000.

At the protected constants,(51008K+2^120)^1000<20000^1000*2^127033. Thus the proposed bound is127033/1000=127.033,above diagnostic127.03224282321547. All trials,failures,hashes,transfers,verification andoutput are included without target amortization.

Payload is2*12*640*2^72 bytes. storage.md specifies2^36 banks of capacity2^50 bytes,allocated capacity2^86 plus<2^45 shared data,reported peak<2^87 including padding,code andmetadata. The bank id is36 bits;all algorithmic arithmetic still uses ordinary64-bit limbs andexplicit charged selector stores/barriers. This finite banked-RAM realization is a material machine-model assumption requiring review,not hardware installed on this VM. If inadmissible,the claim must be rejected rather than silently qualified. No actual exponential arrays are allocated bylocal validation.

Cold2^120 is deliberately loose. At most2^85 scalar attempts at<2^20 instructions cost<2^105. At most2^73 base jumps times<1024 dense products at<2^20 scalar instructions cost<2^103. Variant/coefficient generation,all data writes andbank configuration remain<2^105. Fully writing2^86 bytes andpadding takes<2^90 ordinary instructions even at a conservative per-word allowance. Fixed pairings,input/output checks,constants,checkpoint copies andgeneric hashes are far smaller. Their sum is below2^120. Public profile constants andalgorithm design are excluded as in the rules;no searched advice or free table is hidden.

## Provenance and limitations

Published23d3d4cd-768d-4cac-a9a5-872fc010580d directly supplies the original reduction,norm-one representation andfinite interval method. Published13f0490c-4535-4ee8-8112-6a72a5254024 supplies standard BLAKE3/helper contracts. Exact citations are declared. Source13f was AI-valid127.756,not promoted;that doesnot qualify this orbit proposal. Its BLAKE2s ancestor doesnot establish BLAKE3 H.

Material local predecessors include orbit58640c1a-877a-4a4f-8f50-2c69f8eca147,explicit field64c0a8bf-aee5-4496-9c88-ac60614fc09c,paired orbit e6237cd4-1bee-48d6-935f-a7fa76763f48 andpaired plain ee01563a-29fe-4bd8-81d7-76d8f3ff68d3. They were unpublished/unqualified at preparation andare acknowledged,not called independent research. Current change trades much larger stated memory for two components,with a corrected128-bit scalar counter andfully rederived spectral,fruitless andstorage bounds. No novelty,practical break,scientific qualification orpromotion is claimed.

field-ops.S contains explicit canonical add/sub leaves36/37;48 withinterfaces is retained. Fixed native Fp multiplication416 includes the pinned383 and33 interface/slack. Canonicalization uses10 multiplications andcharged additions/control,with7696 total. The paired baseline uses encodable local scalar offsets afterADD#576 andvalid paired tabulation offsets;it is an accounting template,not a claim that the old one-factor controller is the complete new two-factor program. All added controller work is separately charged. Candidate files are untrusted static evidence,never executed on thehost orwith credentials. Shared llm-verifier owns execution/review/scoring;experiments are empty. Fullparameters,0.39 threshold,referenceprices andsingle computation metric remain unchanged.

## B-half key and constant-digit partial evaluation

Odd subgroup order makes b injective, so the canonical helper drops all a/sign work: selected-pair6 field multiplications plus remaining-b4 gives10, not12. At416 each this costs4160. Three structured Fp2 products have at most12 field operations at48=576; rotations and up tothree conjugations retain19 operations=912. Comparisons768,key/output control768 andwrapper512 give7696 total. No hidden recovery ofa is needed: full rawpoints remain stored for reference multiplication.

The tabulation leaf initializes accumulators from the cold-builtT72[0] constant in four loads and omits its twelve-instruction constant-digit body, giving807. The one-factor accounting template drops three sign-interface instructions and the71-instruction sign leaf:55main+46scalar+37copy+807tabulation+797fullhash+128ordinaryprefix=1870, plus confirmations<=K gives1871K, retained as1888K. Extra controller1280 is unchanged. Total1888+7696+1280=10864, bounded by11008 with144 per-step slack. Two reference multiplications give the stated51008 coefficient. The oldsign table's unusedentry1 remains fully constructed andpaid.

Actual fingerprints/labels change when sign becomesconstantzero, but the paired universality and all probability arguments remain valid: distinct canonicalGTpoints differ in a b digit andthe commonT72[0] cancels. H still concerns actualstandardhash calls,not an unchanged source trajectory. This exact representation refinement introduces no new heuristic. The wide72-bit tables,counter16 and36-bit banked-memory assumption remain fully disclosed andmaterial toreview.
