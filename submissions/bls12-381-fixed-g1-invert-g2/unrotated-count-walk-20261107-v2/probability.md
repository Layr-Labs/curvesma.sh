# Complete coupling, count relations and extraction

All statements are conditional on the disclosed standard full BLAKE3 random-function bridge H,with actual hashes computed/paid. Protected p/r remain unchanged. Let W=2^76,q=12,N=r/q,R=min{R:12R^2>=r},K=L=ceil(1243R/1250),T=128,eta=2^-128. Ordinary raw increments are simply M_j=h^u_j*z^v_j with independently uniform index j. No Frobenius phase is applied to a jump. Full12-orbit detection still uses the exact canonical Fp4trace key. mixing.md/success.md prove stationary finite correlated-walk repeat bounds; raw endpoints and collision events are not assumed independent.

## Exact hash and dictionary coupling

nh-fingerprint.md computes U40 from192-byte exact orbit keys with fixed-pair collision<2^-312,using65 independently domain-separated paid BLAKE3 seed calls. Extend an ideal group walk with fresh independent indices,independent of those fingerprint keys,through its first full-orbit repeat. Before that repeat its exact keys are distinct. Since K+1<2^126,the fixed-pair union bound is<2^251*2^-312=2^-61;retain2^-55. No adaptive-universality claim is made.

JUMP is the full standard BLAKE3-256 message JUMPtag16||seed32||zero16||U40. Prefix64 is compressed/prepared cold;final40 has counter0,len40,flags10 and all56 G/seven rounds. Interpret digest Y little-endian. Its low76 bits select j;the NEXT128 bits give bucket=(Y>>76)mod2^128. These disjoint output bits are exactly independent uniform under H,with no bias/coupling approximation. The remaining52 digest bits are unused. There is no phase,orientation correction,previous edge,inverse guard,fictitious generator,retry hash or rejection cap per transition.

Before the first repeat,distinct U inputs yield fresh independent JUMP digests. Thus indices are iid uniform on W. Conditional on the complete ideal group path/indices,the bucket draws for distinct orbit keys remain exactly independent uniform on2^128 buckets. The final repeated query reuses its old bucket,as explicitly counted in dictionary.md. Its initialized capped-chaining proof has failure<2^-100;all transfers,head zeros,payload initialization and comparisons are charged. Uniform x0 and stationarity bound identity visits by(K+1)/r<2^-100.

## Scalar sampler

SCAL uses the original global unsigned128-bit counter encoded16-byte big-endian in a standard64-byte BLAKE3 message. Accept low255 digest bits interpreted big-endian iff<r;cap1024 attempts per scalar. There are at most2W+1<2^78 accepted scalars,r>7*2^252 implies rejection<1/8,and union scalar-cap loss<2^78*(1/8)^1024<2^-970. All attempted hashes,counters and target-specific group construction are cold paid. No scalar sampler is called on an ordinary transition.

## Every terminal relation has a nonzero count coefficient

Let N_j(t) be the integer number of times index j occurs through prefix t. These are nonnegative monotone counts. The ideal iid index law gives N_j(K) binomial with mean K/W<2^50. Set Hc=2^52. Chernoff bounds Pr[N_j(K)>=Hc]<=(e*(K/W)/Hc)^Hc<(3/4)^Hc. Since(3/4)^4<1/2,union over W indices is<2^76*2^(-2^50)<2^-200. On its complement every prefix count is strictly<Hc. This is an analysis event,not an algorithmic count test or extra heuristic.

For an actual orbit match x_t=xi^ell*x_s,t>s,ell0..11,the denominator coefficient at index j is

 c_j=N_j(t)-xi^ell*N_j(s) modr.

If ell=0,at least one interval count N_j(t)-N_j(s) is a strictly positive integer<Hc<r,so that coefficient is nonzero modulo r.

For ell!=0,fix an index used by time t,so a=N_j(t)>0,b=N_j(s)>=0. Let zeta=exp(pi*i/6),the formal complex primitive twelfth root in Z[zeta],Phi12(X)=X^4-X^2+1. The algebraic integer a-zeta^ell*b is nonzero:if b=0 then a>0;otherwise vanishing would make a/b=zeta^ell positive rational,impossible for a nontrivial twelfth root (the other real root is-1). Its four complex conjugates have magnitude<=a+b<2Hc,so its nonzero integer norm has absolute value<(2Hc)^4=2^212<r. Evaluating zeta at protected xi modulo r cannot yield0,as that would make r divide this norm/resultant. Therefore c_j!=0 for that index. This proves a nonzero coefficient for every twist and every t>s,outside the count exception.

No finite-field value is identified with a complex number. The complex root is used solely to prove a formal polynomial nonzero and bound its integer norm;reduction zeta->xi is justified by Phi12(xi)=0. No small-relation,short-return,long-return,positive-phase or initial-alpha exception is assumed. Actual group returns can occur and are precisely the useful collisions being detected;their coefficient vector remains nonzero.

## Uniform denominator and verified G2 output

Fix all group jumps M_j,initial alpha0,and ideal label/hash randomness through the first orbit match. Since h generates GT,the mapping(u_j,v_j)->(M_j,v_j) is bijective on Fr^2,so v_j remain independent uniform conditional on M_j. The group path,counts,exact canonical keys,NH outputs and unique true twist are group-determined,independent of that remaining scalar freedom. Reconstructing a saved raw point from its stored histories returns exactly its group invariant and hence does not select an alternative scalar-dependent twist.

The terminal denominator D=B_t-xi^ell*B_s=sum_j c_j*v_j is therefore uniform Fr,conditional on that whole group path,because at least one c_j is nonzero. The first terminal match has singular probability1/r,not an unsupported independent-trial argument. Otherwise a=(xi^ell*A_s-A_t)/D satisfies h^a=z;explicitly verify it,form canonical Q=[a]Q0 and check exact pairing/curve/subgroup before output. Reconstruction,three possible raw twists,inversion,scalar multiplication and all verification are bounded cold work. No unchecked or scalar-only output is accepted.

## Complete finite success

success.md gives repeat>0.3900538335845625 at full protected r. Subtract fingerprint2^-55,bad generator table2^-98,scalar cap2^-970,identity2^-100,count exception2^-200,dictionary cap2^-100,and1/r. Success remains>0.39005383358456>0.39. No loss independence is assumed. No phase bias,previous-edge,small-return,lattice martingale,long-return,orientation or post-repeat assumption remains. Standard full BLAKE3 H and explicit finite banked-memory admissibility remain material;this is conditional analytic evidence,not qualification,an official score,a practical attack execution or a fixed-hash theorem.
