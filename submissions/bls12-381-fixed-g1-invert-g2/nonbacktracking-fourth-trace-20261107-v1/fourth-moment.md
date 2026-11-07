# Finite fourth-order success bound for the nonbacktracking twelve-orbit walk

Keep the same full r,q=12,W=2^76,D=qW,N=r/q,T=32,eta=2^-16 and good-table conditions. Change only L=K=ceil(199R/200),where R=min{R:12R²>=r}. The exact trace key,raw transitions,hashes,reference prices,dictionary and terminal extraction are unchanged. This proof adds no heuristic. All statements are inside the already disclosed full standard BLAKE3 bridge H and the same good-table event.

## Multi-time Fourier bound, including short correlated gaps

Let delta=2^-30,b=1/(D-1),a=D/(D-1). For nontrivial raw-group character chi,the previous-edge transition matrix is

 T_chi=aP_chi-bQ_chi,
 P_chi(u,v)=chi(v)/D,
 Q_chi(u,v)=chi(-u)1[v=-u].

The good-table event gives |lambda_chi|<=delta for lambda_chi=D^-1 sum_v chi(v). P_chi is rank one and norm1; Q_chi is unitary. P_chi P_psi=lambda_chi P_psi. At character0,T0 has norm1,preserves the all-ones column and uniform row,and all powers do likewise. In particular P_chi T0^d P_psi=lambda_chi P_psi and P_chi T0^d 1=lambda_chi*1.

For any sequence of endpoint frequencies,the stationary raw path's joint characteristic function is zero unless their sum is zero. Otherwise it is

 w^T T_theta1^d1 ... T_theta_(ell-1)^d_(ell-1) 1,

where w is uniform on D previous edges and theta are the partial-sum frequencies across consecutive time gaps. Let F be the nonzero theta gaps and M=sum_(i in F)d_i. Leave all zero-character stretches as exact norm1 T0 powers. Expand ONLY the M nonzero-character single steps into aP-bQ.

A term with k Q factors has M-k P factors. In each P block,all factors except possibly its last before a Q contribute their small lambda. A final P block ending at the all-ones column also contributes its last lambda. Thus at most k P factors can fail to contribute,and at least M-2k small factors remain. Zero-character stretches do not change this count by the identities above. Remaining operator norms are<=1,and ||w||*||1||=1,so there is no hidden D factor.

Its absolute value is at most a^(M-k)b^k delta^max(M-2k,0),which is also bounded by a^(M-k)b^k delta^(M-2k) even when the exponent is negative. Sum over subsets of Q positions to obtain

 |joint characteristic function| <=(a*delta+b/delta)^M <=rho^M,
 rho=2^-29.

The final strict inequality is checked exactly at protected D. This applies to arbitrarily short/adjacent nonzero gaps and arbitrary changes of character. It does not replace the previous-edge chain by iid increments or claim the arc chain is reversible.

## Collision moments and pointwise Bonferroni

Let C count pairs0<=i<j<=L with j-i>=T and equal FULL Frobenius12 orbits. The exact trace key tests precisely this equivalence. Set B=(L-T+1)(L-T+2)/2. The previous finite first/second-moment proof still gives

 m=(1-eta)B/N <=E C,
 U=(1+3eta)B²/N²+(1+eta)B/N
   +2(q+1)(1+eta)T²(L+1)²/[N(D-1)] >=E C².

At the shorter L,the quadratic term alone is about0.3724730974;it is NOT asserted to exceed0.39.

For the third factorial moment,use only triples of collision edges with six distinct endpoints whose five consecutive gaps are all>=T. Joint-state lower mixing gives factor(1-epsilon)^5,epsilon<2^-127. Uniform independent raw positions have pair orbit probability q/r*(1-(q-1)/(qr));three disjoint pairs give its cube. Both deficits are safely bounded by eta. The number of endpoint choices is binom(L-5(T-1)+1,6),and there are15 matchings. Therefore

 E binom(C,3) >=V3=15(1-eta)^8 binom(L-5(T-1)+1,6)/N³.

For integer c>=0,the EVEN four-term Bonferroni expression is pointwise a lower bound:

 1[c>0] >= c-binom(c,2)+binom(c,3)-binom(c,4).

This follows from the alternating binomial identity,not event independence. Thus

 Pr[first orbit repeat by L] >=(3m-U)/2+V3-U4

for any upper bound U4>=E binom(C,4). We derive that bound below;no Poisson or independent-pair approximation is used.

## Weighted orbit-constraint Fourier graphs

A term in binom(C,4) is an unordered set H of four distinct collision edges on ell=4..8 ordered time endpoints,covering every endpoint. Let s be the number of equality components and rank=ell-s<=4. Inside each component,fix its root orientation and assign each other vertex a multiplier xi^h,h0..11. There are q^rank assignments. The orbit-equality indicator is bounded above by their sum of linear-constraint indicators;identity vertices may be overcounted,which is harmless for an upper bound.

For a fixed assignment,the constraints x_v=zeta_v*x_root have rank rank and Fourier normalizer r^-rank. The dual frequencies satisfy sum_(v in component) zeta_v*psi_v=0. Uniform raw starting position kills total-frequency-nonzero terms. Parameterize the remaining endpoint frequencies by partial sums theta across adjacent time gaps,with nonzero support F.

Contract H components. Each selected adjacent gap becomes a weighted edge between its two component labels,with nonzero incidence weights zeta at its endpoints. A loop has the corresponding difference weight and may impose extra constraints. If some active vertex has degree1,its equation forces that gap's theta to zero,so the term vanishes. Otherwise let f=|F|,v=number of incident vertices,c=number of active connected components,and bF=f-v+c.

The weighted incidence rank is at least v-c:spanning-forest columns are independent by leaf elimination,since all nonloop endpoint coefficients are nonzero. Loops or inconsistent rotations can only increase rank. Hence there are at most r^bF possible frequency assignments. Combined with the multi-time lemma,the contribution for a fixed endpoint configuration is bounded by

 q^rank * r^(-rank+bF) *rho^(sum selected gap lengths).

This explicitly handles phase rotations and previous-edge dependence. It does not import the iid plain-walk graph bound without proof.

## Finite endpoint templates and all error classes

templates.md specifies the complete finite enumeration,with no secret input or candidate execution. For each ell4..8,enumerate covering four-edge H and every subset F of its ell-1 adjacent gaps;discard degree1 terms. Mark F forced if it contains the whole time path of any collision edge of H. Such an edge has length>=T,so total selected gap length is>=T. Define

 E=-rank+bF+(ell-f)/2.

The verified finite certificate has900 covering graphs,22246 retained H/F terms including empty F,no positive-E nonempty free term,and all negative free/empty terms have E<=-1/2. Forced terms have E<=2. Empty ell8 terms are the105 perfect matchings;all other empty terms are negative. Free E0 counts by f1..6 are

 90,139,122,121,72,24.

The actual coefficient table and full language-independent pseudocode are in templates.md. Two independently written public combinatorial checks (union-find and component sets) agree. These checks concern only finite graph integers,not an attack execution,production success measurement or extra assumption.

At the protected shorter L,we have L+1<sqrt(N),checked exactly. For fixed selected gap lengths,there are at most(L+1)^(ell-f) choices of start and remaining gaps. Thus the contribution after summing endpoint choices is bounded by a factor

 q^(rank-(ell-f)/2) *r^E

times its selected-gap geometric sum.

For free E0 terms,the q exponent equals bF<=f. With t=q*rho/(1-rho),their total is at most

 90t+139t²+122t³+121t⁴+72t⁵+24t⁶.

For negative terms,the q factor is<=q^rank<=q^4,the r factor<=r^-1/2<2^-127,and the gap sum<=1. All such terms together are bounded by22246*q^4*2^-127<2^-90.

For forced terms,sum of selected gap lengths S>=T. Write rho^S=rho^(3S/4)*rho^(S/4). Since rho^(1/4)<1/2,the sum over all positive selected lengths is<=rho^(3T/4)=rho^24;the remaining geometric factor is<1. The q factor<=q^4 and r^E<=r². All forced terms together are bounded by22246*q^4*r²*rho^24<2^-140. These bounds include every positive-E diagram;none is hidden in an asymptotic remainder.

Combine the two remainder classes conservatively as2^-80. The leading empty-F ell8 term is105 binom(L+1,8)/N⁴. Therefore the complete upper bound is

 U4=105 binom(L+1,8)/N⁴
    +90t+139t²+122t³+121t⁴+72t⁵+24t⁶+2^-80.

All inequalities use full protected r,D and finite integer/rational arithmetic.

## Full finite success and complete computation

Exact protected evaluation gives

 (3m-U)/2 ≈0.3724730973788035,
 V3≈0.020213626291806727,
 U4≈0.002503816457812061,
 repeat lower bound≈0.39018290721279814 >0.39018.

Subtract the unchanged fingerprint2^-63,table2^-98,scalar/local2^-970,initial/retry-phase2^-110,ordinary-phase2^-50,identity2^-100,short-fruitless2^-16 and2^-32,other-fruitless2^-199,global-retry2^-45,dictionary-cap2^-100 andsingular1/r losses. No canonical cost cap exists in this trace version. The resulting exact lower bound is >0.3901676481908756 >0.39. No failure independence is claimed.

All instruction ledgers,full parameter sizes,reference priceC20000,one computation metric,raw GT buffers,exact Fp4trace key,cold raw reconstruction/three-twist recovery,terminal G2 verification,physical capacity<2^136 andzero advice are unchanged. Only the maximum transitions is shortened toK=ceil199R/200. Complete coldD=540*2^119 remains paid per target,not amortized. Therefore

 W<=(23776K+540*2^119)/20000.

Exact integer certificate(23776K+540*2^119)^1000<20000^1000*2^125879 establishes proposed125879/1000=125.879,above diagnostic125.87866896466751. This is analytic conditional evidence under H and explicit finite banked-memory admissibility,not an official score,qualification,production execution,fixed-hash theorem or practical hardware claim.
