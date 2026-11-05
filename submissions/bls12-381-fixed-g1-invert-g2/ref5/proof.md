# Reference-primitive supplied-P inversion
Insns mean scalar AArch64 instructions.
Claim:log2 cold work<=128.99,success>=0.39 conditional on H,memory<2^35 bytes,advice zero.

Algorithm uses protected p,r,Q0 and pinned pairing(Q,P). z=1 returns infinity. Compute h=e(P,Q0),a generator by prime order/nondegeneracy. Solve z=h^a,return[a]Q0 with subgroup/exact pairing verification. C=20000 unchanged.

INIT,SCAL,WALK are16-byte tags FAPI-RHO-INIT-v1,FAPI-RHO-SCAL-v1,FAPI-RHO-WALK-v1. Hashes are actual standard SHA-512. seed=SHA512(INIT||P_bytes||z_bytes),using protected96/576-byte encodings. Sample2w+1 scalars,w=2^24:low255bits SHA512(SCAL||seed||counter_be64),reject>=r,cap1024 attempts,exhaustion fails. Counter globally increases. Store w records M_j=h^u_j*z^v_j,u_j,v_j;start(X,A,B)=(h^alpha,alpha,0). Counter<2^36.

Use canonical Montgomery tower u^2=-1,v^3=1+u,w^2=v,equivalent to protected u=w^6-1,v=w^2. Write X=a+bw. Norm one implies a^2-vb^2=1. Encode six canonical48-byte coeffs of b and one sign byte comparing a lexicographically with -a. The289-byte encoding is injective:given b,roots +/-a are distinguished by sign,with zero designated. Select j=low24bits SHA512(WALK||seed||zero48||encoding(X)),zero48 being48 zero bytes;update X<-X*M_j,A<-A+u_j,B<-B+v_j mod r. Preserve X=h^A*z^B.

R=isqrt(r-1)+1,L=(1167R+999)//1000,s=(L+1023)//1024,K=L+s. Save640-byte(X,A,B) at0,s,...,1024s in sorted-pointer order by72 canonical64-bit X words. For t=1..K advance,search equality before saving. First equality gives a=(A_saved-A)/(B-B_saved);zero denominator fails. Verify/return Q. Otherwise insert pointer when saving. At most1025 distinct immutable keys;decisions inspect group values only. No equality by K fails.

H alone idealizes computed domain-separated SHA512 as indep. random functions on disjoint inputs over uniform challenges. Fixed SHA behavior is unproved;work bounded regardless. No random-mapping walk assumption.

Under H samplers are uniform indep.,failure<2^-990;M_j uniform indep. GT. Real/imaginary Hoeffding gives Pr(|w^-1 sum chi(M_j)|>rho=1/128)<=4exp(-w/65536) per nontrivial character. Union over r-1 gives bad-table Pr<2^-100:4r*3^256*2^100<8^256,using e>8/3. No table test.

On good tables take indep.-index surrogate increments. Fourier/Parseval give sup and unweighted l2 norms of epsilon_d=p_d-1/r at most rho^d. Put T=271,eta=1/65536,D=(L-T+1)(L-T+2)/2. Count collisions C with interval length>=T. Then E C>=m=(1-eta)D/r since rho^T<=eta/r.

Second moment:identical intervals contribute(1+eta)D/r. Give distinct ordered pairs baseline(1+3eta)/r^2. Disjoint/touching increments are independent. Nested gaps a,c>=0,b>=T,a+c>0 have joint p_b(1)p_(a+c)(1)<=(1+eta)/r^2+(1+eta)rho^(a+c)/r. Geometric sum over a,c is255/16129;both orders,start/b choices give excess510(1+eta)(L+1)^2/(16129r). Crossing positive gaps a,b,c,a+b,b+c>=T have joint sum_x p_a(x)p_b(x^-1)p_c(x)<=1/r^2+(rho^(a+b)+rho^(a+c)+rho^(b+c))/r+rho^(a+b+c),by expansion/Cauchy-Schwarz. Summing rho^(a+c) gives2(L+1)^2/(16129r);sum_(a,b>=1,a+b>=T)rho^(a+b)<=2T*2^-T bounds the other terms by8T eta(L+1)^2/r^2. Constrained sum rho^(a+b+c)<6*2^-T:dominate by rho=1/2,split b>=T-1 (<=4*2^-T) and b<T-1 (<2*2^-T). Both orders/start give12eta(L+1)/r. These exhaust ordered interval patterns,including shared endpoints.

Thus E C^2<=U=(1+3eta)D^2/r^2+(1+eta)D/r+((512+510eta)/16129)(L+1)^2/r+8T eta(L+1)^2/r^2+12eta(L+1)/r. Exact protected-r rational arithmetic gives m^2>(3903/10000)U;Cauchy-Schwarz yields repeat Pr>0.3903 by L.

Couple actual ideal-hash iteration to the indep.-index surrogate until first repeat:every distinct group element has a distinct encoding,so its previously unused complete hash message receives a fresh uniform label. First-repeat laws agree.

For tail mu,cycle lambda,tau=mu+lambda<=L,the saved anchor a=ceil(mu/s)s satisfies mu<=a<=mu+s-1,a<=1024s. Equality recurs at a+lambda<=tau+s-1<K and binary search detects it. Conditional on group jumps,partition and initial group value,each v_j stays uniform independently,as exactly one u_j corresponds to each v_j at fixed M_j. Detected denominator=sum n_j v_j,sum n_j=d,0<d<=K<r,so some n_j!=0 mod r. Singularity Pr=1/r,with no coeff-based selection. Hence success>0.3903-2^-100-2^-990-1/r>0.39 under H.

Cost uses the protected H_ref+R/C rule:each X*M_j is exactly one pinned blst_fp12_mul dense reference operation,charged1 including its internals. External work follows.

Normalize every384-bit output coeff by16 unrolled conditional p subtractions: p>2^380 implies input<16p. Each six-limb subtraction/selection costs12 scalar insns;12 coeffs*(16*12+12 loads/stores)=2448. Charge4096 including constants,pointers,call-site/control. Negating six Fp coeffs and lexicographically comparing a with -a fits2048 insns.

Seed has64 bytes;SCAL88-byte messages pad to one128-byte block. WALK417-byte message has fixed128-byte prefix(tag16,seed64,zero48),precomputed once;289-byte suffix,0x80,78 zeros,16-byte big-endian3336 pad to three blocks. For each block,unroll64 SHA512 schedule words:two small sigmas(10 insns),four loads,three adds,one store=18<=20. Each of80 rounds:two big sigmas10,Ch3,Maj5,seven adds,two loads,eight moves=35<=40. Use scalar64-bit operations per FIPS180-4 sections5.1.2,6.4.2. Expand schedule before loading state;registers suffice,all offsets<640. Word loads/reversals32,chaining loads8,feed-forward24,prologue/pointers/return64 give128 overhead. Three blocks cost3*(80*40+64*20+128)=13824;2048 formatting/state copies/padding/extraction gives15872. No free hash.

Binary search:unroll72 words*(two LDR,CMP,B.ne)=288,plus96 caller/pointer/midpoint/lo/hi/control=384 per key. Eleven comparisons=4224. Add1024 two four-limb scalar updates,256 jump index/address,512 multiword counters/save control,512 record copy,256 dispatch,1024 spills/flags:7808<=8192. At most1025*1024 pointer shifts*16 insns fit2^30 fixed work.

Per step=1 reference multiply+(4096+2048+15872+8192)/20000. Fixed2^90 insns cover table,samplers,basis/input/subgroup checks,seed/prefix,insertions,coeff inversion,pairings,verification,output,failures. Bit-serial384-bit product/reduction<2^24 scalar operations gives Fp12 product<2^34,inverse by p^12-2<2^48. Table<2^35 products costs<2^69;samplers<2^36 hashes*65536<2^52. Two pairings each<2^20 extension operations even at2^58 each cost<2^79 total. Complete affine point formulas test identity,equal/opposite,y=0 before division;Fp2 inverse=(a-bu)/(a^2+b^2). Even4096 point steps at2^52 plus fixed work<2^90. No amortization.

Exact check:(K*50208+2^90)^100<20000^100*2^12899,so score128.99. Storage2^24*640=10737418240 plus1026 states,1025 pointers,code/constants/scratch/seed/output<2^35 bytes;offsets<2^34. Advice zero;no host execution.

Uses declared reduction/encoding and checkpoint/horizon sources. New:direct SHA512 and reference primitive. Shared verifier scores.
