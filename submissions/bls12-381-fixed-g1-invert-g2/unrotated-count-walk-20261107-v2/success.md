# Finite six-term orbit-collision bound for the one-sided walk

Keep full protected r,q=12,W=2^76,D=W,N=r/q,T=128,eta=2^-128,rho=2^-33 and K=L=ceil(1243R/1250),R=min{R:12R^2>=r}. mixing.md proves an exact independent-increment Fourier product bound rho^S for any nonzero partial-frequency support, including arbitrary short gaps. No previous-edge matrix, endpoint independence, reversibility, birthday approximation or extra heuristic is used.

Let C count pairs0<=i<j<=L with j-i>=T and equal FULL twelve-element Frobenius orbits. The exact trace key identifies this equivalence. Write B=(L-T+1)(L-T+2)/2. A stationary raw value is uniform; its mean orbit size is q-(q-1)/r. Mixing relative error<(r-1)rho^T<2^-3900 and identity deficit(q-1)/(qr)<2^-250 together are<eta. Hence

 m=(1-eta)B/N <=E C.

For odd factorial moments retain only disjoint edges with all consecutive endpoint gaps>=T. The exact adjusted endpoint counts and numbers of perfect matchings give

 V3=15(1-eta)^8 binom(L-5(T-1)+1,6)/N^3 <=E binom(C,3),
 V5=945(1-eta)^14 binom(L-9(T-1)+1,10)/N^5 <=E binom(C,5).

Each full-state kernel deficit and each pair identity deficit is<eta; five/nine consecutive kernel factors and three/five pair factors explain exponents8/14. This is a finite joint-kernel lower bound, not independence of actual endpoints.

## General diagram bound without enumerating large graphs

For k=2,4 or6, expand E binom(C,k) by unordered sets H of k distinct collision edges on ell covered, chronologically ordered endpoints, ell<=2k. H has s components and rank h=ell-s<=k. Root each component and assign its other vertices one of q nonzero orbit multipliers. Summing these q^h assignments bounds the orbit-constraint indicator above even at identity. Each fixed assignment has Fourier normalization r^-h.

Let F be the nonzero partial-frequency support on the ell-1 time gaps, f=|F|. Contract H components. Selected gaps become weighted edges on these components. If an active contracted vertex has degree1, its frequency equation forces that selected frequency zero, so discard that term. Let v be the active vertex count and c its nonempty component count; bF=f-v+c. Weighted spanning-forest columns are independent by leaf elimination, so the number of allowable frequencies is at most r^bF, regardless of loop weights or rotations. The contribution at fixed endpoint times is at most q^h r^(-h+bF) rho^S.

Call F forced if the whole adjacent-gap path between the endpoints of some H edge lies in F. Then S>=T because every H edge is one of the counted long-lag collisions. Otherwise call it free. Fixed selected gap lengths leave at most (L+1)^(ell-f) choices of start and remaining gaps. The protected cap satisfies L+1<sqrt(N), exactly. Thus the prefactor is bounded by q^(h-(ell-f)/2) r^E, with E=-h+bF+(ell-f)/2.

Here is a structural proof of E<=0 for every free diagram, replacing a graph enumeration. Let u be the component count of H union F. Since H starts with s components, u=s-v+c, so E=u-(ell-f)/2. Contract the f edges of F in the ordered time path, leaving ell-f blocks. Freeness means every H edge joins different blocks. Since H covers every endpoint, every block has an incident external H edge. Consequently each connected component of the block graph has at least two blocks, u<=(ell-f)/2 and E<=0. Negative exponents are half-integers <=-1/2. For every diagram, including forced ones, u<=s<=ell/2, giving E<=f/2<k. The safe upper bound E<=k suffices.

Empty F has E=s-ell/2. Equality means H is a perfect matching on ell=2k vertices; all other empty diagrams have E<=-1/2. Thus the leading terms are (2k-1)!! binom(L+1,2k)/N^k.

For free E=0 and f=1, the selected gap must be a loop after contracting H, since a nonloop gives degree1. Its two adjacent time vertices belong to the same H component but are not joined by an H edge (which would force the diagram). E=0 implies s=(ell-1)/2: exactly one H component has three vertices and all others have two. The three-vertex component must be a two-edge path with the selected adjacent pair as its leaves. A triangle would include the forbidden leaf edge. Hence ell=2k-1. Choose the adjacent leaf pair in 2k-2 ways, its center in 2k-3 ways, and match the remaining vertices in (2k-5)!! ways, with(-1)!!=1. This gives exactly2 diagrams for k2,90 for k4 and9450 for k6.

Every other free E0 term has f>=2. Since E0 makes h-(ell-f)/2=bF<=f, summing positive selected gap lengths bounds each such diagram by x^f<=x^2, where x=q*rho/(1-rho)<1. Total diagram counts, including all supports and discarded terms, are safely bounded by

 M2<=3*binom(6,2)*2^3<2^9,
 M4<=5*binom(28,4)*2^7<2^24,
 M6<=9*binom(66,6)*2^11<2^41.

Indeed ell ranges2..4,4..8 or4..12; each H is a k-subset of at most6,28 or66 vertex pairs and each F is a subset of at most3,7 or11 gaps. Overcounting only enlarges an upper bound.

All negative free/empty terms together are <=M_k q^k r^-1/2 <=M_k q^k 2^-127, since the selected-gap geometric sums are <=1. For forced terms S>=128; split rho^S as rho^(3S/4)*rho^(S/4). Since rho^(1/4)<2^-8, the remaining geometric sums are <1 for nonempty F. These terms total at most M_k q^k r^k rho^96. Exact protected integer checks give a combined remainder <2^-100 for k2,<2^-80 for k4 and <2^-50 for k6. These bounds cover every positive-exponent diagram and every short correlated gap.

Therefore valid complete upper bounds are

 U2=3 binom(L+1,4)/N^2+2x+2^9*x^2+2^-100,
 U4=105 binom(L+1,8)/N^4+90x+2^24*x^2+2^-80,
 U6=10395 binom(L+1,12)/N^6+9450x+2^41*x^2+2^-50.

## Finite success and complete computation

For every integer c>=0, the EVEN six-term Bonferroni polynomial satisfies

 1[c>0]>=c-binom(c,2)+binom(c,3)-binom(c,4)+binom(c,5)-binom(c,6).

Insert the lower odd and upper even factorial-moment bounds in their correct signed directions:

 Pr[first full orbit repeat by L]>=m-U2+V3-U4+V5-U6
   >0.3900538335845625.

All quantities are evaluated with exact protected integers/rationals, not a Poisson or asymptotic collision formula. probability.md proves complete coupling/extraction failures: fingerprint2^-55,bad table2^-98,scalar cap2^-970,identity2^-100,count/norm exception2^-200,dictionary cap2^-100,singular denominator1/r. Their union yields success>0.39005383358456>0.39. The integer-count/formal-norm argument eliminates fruitless-return losses,phase selection and inverse retries,without asserting that actual finite-field group returns cannot occur.

Ordinary scalar work<=1408canonical+460fingerprint+776fullBLAKE3+46Frupdates+72controller+192head+224payload+32cappedvisits=3210,retained3224 with14 extra spare instructions. Every transition performs one protected reference C20000. Complete per-target coldD=540*2^119 is unchanged and includes the initial key/hash,all paid raw jumps,head zeros,seed/generator/scalar construction,code/constants and verified G2 finishing. Therefore work<=(23224K+540*2^119)/20000,diagnosticlog2=125.84391726376873. Exact certificate(23224K+540*2^119)^250<20000^250*2^31461 supports proposed31461/250=125.844. No reference price, full parameter, success threshold, six-track structure or computation metric is changed. Memory<2^136,zero advice. H and explicit finite banked-memory admissibility remain material. This is conditional analytic evidence,not an official score,qualification,promotion,practical execution or fixed-hash theorem.

Material local sources7a74d390-af9f-4837-bd36-c9ff9678d232/04b8156b-b862-4793-bbc2-36c94bf83bed/ae4d5871-142e-4aa9-8b9b-cb904b2eea31 supplied the earlier graph/weighted-flow/even-Bonferroni method. The present exact independent-increment product,second-factorial upper bound,tighter eta and one-sided denominator argument are derived explicitly. Unpublished-source formal links are deferred,not attribution omitted or qualification inherited.
