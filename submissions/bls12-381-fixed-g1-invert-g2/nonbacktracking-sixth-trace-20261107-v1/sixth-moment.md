# Finite sixth-order success bound

All protected field/group parameters, reference price C=20000, hashes, exact Fp4 trace key, physical dictionary, cold initialization and output verification remain unchanged. Only the accepted-transition cap becomes L=K=ceil(1989R/2000), where R=min{R:12R^2>=r}. Analysis uses q=12, W=2^76, D=qW, N=r/q, eta=2^-16, T=128. T is an analysis lag, not an extra runtime operation or checkpoint.

## Stronger good-table event

Use delta=2^-33. The same orbit-average Hoeffding argument in mixing.md gives failure <=4r exp(-256)<2^-100. A wholly integer certificate is 4r*3^256*2^100<8^256, since e>8/3. Distinct/nonzero table failure remains <2^-99, so their union remains <2^-98. No search, test or free conditioning operation is performed by the algorithm. For a good table, ||T_chi^2||<=a^2 delta+2ab+b^2<2^-32, where a=D/(D-1), b=1/(D-1). Full-state relative mixing error is <2^-127 already at lag32 and hence at T=128.

The multi-frequency operator lemma is as follows: rho=a delta+b/delta<2^-32. In its statement we retain the upper bound rho0=2^-32. For any collection of ordered endpoint frequencies with total zero and nonzero partial-sum support F, absolute joint characteristic function <=rho0^S, S=sum of selected positive time-gap lengths. This follows by expanding only nonzero-character steps as aP-bQ, leaving zero-character stretches exact. Here P_chi(u,v)=chi(v)/D is rank one of norm1, Q_chi(u,v)=chi(-u)1[v=-u] is unitary, and w is the uniform previous-edge row. Uniform raw starting position kills total-frequency-nonzero terms; otherwise the characteristic function is w^T product(T_theta_i^d_i)1. P_chi P_psi=lambda_chi P_psi and P_chi T0^d P_psi=lambda_chi P_psi. At most k of M-k P factors lose their small lambda in a term with k Q factors. Summing gives (a delta+b/delta)^M. There is no endpoint dimension factor, independence assumption, or replacement of the previous-edge chain by a raw-position Markov chain.

## Collision count and odd lower moments

Let C count equal full twelve-orbit pairs among times0..L with lag>=T. Write B=(L-T+1)(L-T+2)/2. The finite first/second-moment argument in probability.md, valid for any T>=32, gives

 m=(1-eta)B/N <=E C,
 U=(1+3eta)B^2/N^2+(1+eta)B/N
    +2(q+1)(1+eta)T^2(L+1)^2/[N(D-1)] >=E C^2.

For j=3 or5 retain only j disjoint collision edges with all 2j-1 consecutive endpoint gaps>=T. Full-state kernel lower bounds and free orbit size q for nonidentity points give

 V3=15(1-eta)^8 binom(L-5(T-1)+1,6)/N^3 <=E binom(C,3),
 V5=945(1-eta)^14 binom(L-9(T-1)+1,10)/N^5 <=E binom(C,5).

The factors account for respectively5 or9 kernel deficits and3 or5 pair identity deficits; each deficit is much smaller than eta. Endpoint subsets have the displayed exact gap-adjusted binomial counts; the factors15 and945 count perfect matchings. No pair-independence claim is used for the actual path.

## General diagram bound without enumerating large graphs

For k=4 or6, expand E binom(C,k) by unordered sets H of k distinct collision edges on ell covered, chronologically ordered endpoints, ell<=2k. H has s components and rank h=ell-s<=k. Root each component and assign its other vertices one of q nonzero orbit multipliers. Summing these q^h assignments bounds the orbit-constraint indicator above even at identity. Each fixed assignment has Fourier normalization r^-h.

Let F be the nonzero partial-frequency support on the ell-1 time gaps, f=|F|. Contract H components. Selected gaps become weighted edges on these components. If an active contracted vertex has degree1, its frequency equation forces that selected frequency zero, so discard that term. Let v be the active vertex count and c its nonempty component count; bF=f-v+c. Weighted spanning-forest columns are independent by leaf elimination, so the number of allowable frequencies is at most r^bF, regardless of loop weights or rotations. The contribution at fixed endpoint times is at most q^h r^(-h+bF) rho0^S.

Call F forced if the whole adjacent-gap path between the endpoints of some H edge lies in F. Then S>=T because every H edge is one of the counted long-lag collisions. Otherwise call it free. Fixed selected gap lengths leave at most (L+1)^(ell-f) choices of start and remaining gaps. The protected cap satisfies L+1<sqrt(N), exactly. Thus the prefactor is bounded by q^(h-(ell-f)/2) r^E, with E=-h+bF+(ell-f)/2.

Here is a structural proof of E<=0 for every free diagram, replacing a graph enumeration. Let u be the component count of H union F. Since H starts with s components, u=s-v+c, so E=u-(ell-f)/2. Contract the f edges of F in the ordered time path, leaving ell-f blocks. Freeness means every H edge joins different blocks. Since H covers every endpoint, every block has an incident external H edge. Consequently each connected component of the block graph has at least two blocks, u<=(ell-f)/2 and E<=0. Negative exponents are half-integers <=-1/2. For every diagram, including forced ones, u<=s<=ell/2, giving E<=f/2<k. The safe upper bound E<=k suffices.

Empty F has E=s-ell/2. Equality means H is a perfect matching on ell=2k vertices; all other empty diagrams have E<=-1/2. Thus the leading terms are (2k-1)!! binom(L+1,2k)/N^k.

For free E=0 and f=1, the selected gap must be a loop after contracting H, since a nonloop gives degree1. Its two adjacent time vertices belong to the same H component but are not joined by an H edge (which would force the diagram). E=0 implies s=(ell-1)/2: exactly one H component has three vertices and all others have two. The three-vertex component must be a two-edge path with the selected adjacent pair as its leaves. A triangle would include the forbidden leaf edge. Hence ell=2k-1. Choose the adjacent leaf pair in 2k-2 ways, its center in 2k-3 ways, and match the remaining vertices in (2k-5)!! ways. This gives exactly90 diagrams for k4 and9450 for k6.

Every other free E0 term has f>=2. Since E0 makes h-(ell-f)/2=bF<=f, summing positive selected gap lengths bounds each such diagram by x^f<=x^2, where x=q*rho0/(1-rho0)<1. Total diagram counts, including all supports and the discarded terms, are safely bounded by

 M4<=5*binom(28,4)*2^7<2^24,
 M6<=9*binom(66,6)*2^11<2^41.

Indeed ell ranges4..8 or4..12; each H is a k-subset of at most28 or66 vertex pairs and each F is a subset of at most7 or11 gaps. Overcounting only enlarges an upper bound.

All negative free/empty terms together are <=M_k q^k r^-1/2 <=M_k q^k 2^-127, since the selected-gap geometric sums are <=1. For forced terms S>=128; split rho0^S as rho0^(3S/4)*rho0^(S/4). Since rho0^(1/4)=2^-8, the remaining geometric sums are <1 for nonempty F. These terms total at most M_k q^k r^k rho0^96. Exact protected integer checks give a combined remainder <2^-80 for k4 and <2^-50 for k6. These bounds cover every positive-exponent diagram and every short correlated gap.

Therefore valid complete upper bounds are

 U4=105 binom(L+1,8)/N^4+90x+2^24*x^2+2^-80,
 U6=10395 binom(L+1,12)/N^6+9450x+2^41*x^2+2^-50.

## Finite success and computation certificate

The pointwise EVEN six-term Bonferroni identity gives, for every integer C>=0,

 1[C>0]>=C-binom(C,2)+binom(C,3)-binom(C,4)+binom(C,5)-binom(C,6).

Using the signed bounds in their correct directions,

 Pr[first full orbit repeat by L]>=(3m-U)/2+V3-U4+V5-U6
   >0.3900650822981647.

This is a finite rational evaluation at the full protected r, not an asymptotic birthday or Poisson approximation. Subtract all disclosed failure losses: fingerprint2^-63, bad table2^-98, scalar/local caps2^-970, initial/retry phase2^-110, ordinary phase2^-50, identity2^-100, short fruitless2^-16 and2^-32, other fruitless2^-199, global retry cap2^-45, dictionary cap2^-100 and singular denominator1/r. Their union yields success>0.3900498232762422>0.39. No independence of losses is needed.

All instruction counts and costs are unchanged. Work <=(23776K+540*2^119)/20000, with diagnostic log2=125.87794397649496. The integer certificate (23776K+540*2^119)^500<20000^500*2^62939 establishes the proposed bound62939/500=125.878. The reference price, one computation metric, full field/order, six production tracks, all per-target cold work, zero advice, peak memory<2^136 and verified G2 finishing remain intact. The existing full standard BLAKE3 bridge H and finite banked-memory admissibility remain explicit material conditions. This is analytic conditional evidence, not official scoring, qualification, production execution or a fixed-hash theorem.

Material local sources: 04b8156b-b862-4793-bbc2-36c94bf83bed supplies the previous-edge multi-frequency/weighted-rotation graph method. Earlier ae4d5871-142e-4aa9-8b9b-cb904b2eea31 supplied endpoint-template/even-Bonferroni methods. These are acknowledged material sources, not independent discoveries or inherited qualifications. Their formal links are deferred while unpublished. The structural free-graph proof and six-term finite refinement are explicitly derived here.
