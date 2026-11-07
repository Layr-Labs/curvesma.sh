# Direct six-shape second-moment proof

Full protected parameters:q=12,N=r/q,W=2^94,rho=2^-42,T=128,J=2T=256,eta=2^-128,c=1+eta,R=min{integer R:12R^2>=r},K=ceil2061R/2000,n=K+1. This proof needs only stationary raw increments,pointwise long-gap mixing and a short-gap maximum. It uses no weighted frequency diagrams,higher moments,random mapping,Poisson approximation or independent collision events.

## Good table and two elementary kernel bounds

For nonidentity z=h^a,a!=0,the sampled jumps h*z^v_j and initial h*z^beta0 have independent uniform exponents because multiplication by a is bijective. Table zero/duplicate probability<=W/r+W(W-1)/(2r)<2^-67. For each nontrivial additive character,lambda_chi=W^-1 sum_j chi(1+a*v_j) averages independent zero-mean complex unit values. Hoeffding on real/imaginary parts and union over r-1 characters gives spectral failure<=4r exp(-W*rho^2/4)=4r exp(-256)<2^-100. Integer certificate4r*3^256*2^100<8^256 uses e>8/3. Combined good-table failure<2^-66.

Fix ANY good table. Extend an ideal iid-index walk through K only for analysis. Raw exponents x_t have uniform stationary start and independent increment BLOCKS on disjoint chronological gaps. Let p_d(y) be the d-step increment probability. Fourier inversion gives relative error<= (r-1)rho^d;for d>=T,this is<2^-3900. In particular p_d(y)<=c/r. Since the W one-step increments are distinct,p_1(y)<=1/W. Convolution of probability distributions cannot increase their maximum,so p_d(y)<=1/W for EVERY d>=1,including short gaps. No reversibility or independent endpoints are assumed.

## First moment

Let C count chronological pairs0<=i<j<=K with j-i>=J and equal FULL12-element Frobenius orbits. There are B=(K-J+1)(K-J+2)/2 pairs. Long-gap relative error<2^-3900 and identity orbit deficit11/(12r)<2^-250 together are<eta. Thus m=(1-eta)B/N<=E C.

For an upper bound,an orbit equality can be expanded into q equations x_later=alpha*x_earlier with alpha in A={xi^ell:0<=ell<12}. At identity this overcounts,which is safe. Two distinct collision edges cover three or four distinct ordered times. Exactly three two-edge trees on three vertices and three perfect matchings on four vertices are possible. For each shape sum q^2 phase assignments. Write consecutive time gaps as d1,d2[,d3]. All counted edge lags are>=J.

## The three shared-endpoint shapes

For the middle-center tree(edges01,12),both gaps are long,so each fixed-phase joint probability<=c^2/r^2 by conditioning successive kernels.

For the first-center tree(edges01,02),write x1=alpha*x0,x2=beta*x0. Gap d1>=J is long. If d2>=T,both constraints again cost<=c^2/r^2. If d2<T,the exact probability is r^-1 sum_x p_d1((alpha-1)x)*p_d2((beta-alpha)x). For beta!=alpha,the second argument is a bijection in x,so bounding the first kernel by c/r gives<=c/r^2. For beta=alpha,use p_d2(0)<=1/W to get<=c/(rW). There are only q such equal-phase assignments.

For the last-center tree(edges02,12),write x2=alpha*x0,x1=beta*x0 by a bijective reparametrization of its q^2 assignments. Gap d2>=J is long. If d1>=T,use c^2/r^2. Otherwise its probability is r^-1 sum_x p_d1((beta-1)x)*p_d2((alpha-beta)x). For beta!=1,averaging the first kernel gives<=c/r^2;for beta=1,use its short maximum1/W,at most q bad assignments.

Therefore all three-endpoint shapes contribute at most3q^2*c^2*binom(n,3)/r^2 plus2q*c*n^2*T/(rW). Indeed for each exceptional tree there are at most n^2*T placements with its indicated short adjacent gap. We deliberately overcount times and phase assignments.

## The three four-endpoint shapes

Separated edges01,23 have d1,d3>=J. Condition on the intermediate state and bound these two kernels by c/r,for<=c^2/r^2 per assignment,regardless of d2.

Crossing edges02,13 impose x2=alpha*x0,x3=beta*x1. If d2<T,the two counted lags imply d1,d3>T because J=2T. Bound both outer kernels by c/r and sum the remaining normalized middle kernel;the result is<=c^2/r^2. Explicitly sum r^-1 p_d1(x1-x0)*p_d2(alpha*x0-x1)*p_d3(beta*x1-alpha*x0) over x0,x1.

If crossing d2>=T,let U,V be the independent outer increment blocks. The long middle kernel is bounded by c/r,and the remaining constraint is V=(beta-alpha)*x0+beta*U. For beta!=alpha,each U,V determines exactly one uniform x0,so the bound is<=c/r^2. For beta=alpha,it is V=beta*U. If either outer gap>=T,one outer distribution has maximum c/r,giving<=c^2/r^2. Only when BOTH outer gaps<T do we use maximum1/W,giving<=c/(rW),with q equal-phase assignments.

Nested edges03,12 always have d2>=J. Write x2=beta*x1,x3=alpha*x0 and outer blocks U=x1-x0,V=x3-x2. Bound the long middle kernel by c/r;the remaining constraint is V=(alpha-beta)*x0-beta*U. If alpha!=beta,average the unique uniform x0 for<=c/r^2. If alpha=beta,it is V=-beta*U. One long outer gap gives maximum c/r;only both short outer gaps require1/W. Again at most q assignments cost c/(rW).

Thus all four-endpoint shapes contribute<=3q^2*c^2*binom(n,4)/r^2 plus2q*c*n^2*T^2/(rW). For each exceptional crossing/nested shape choose the two outer short gaps in at most T^2 ways and the remaining two positions in at most n^2 ways. Every use of uniform x0 is valid because it is independent of subsequent increment blocks on a stationary raw Cayley walk.

## Complete finite certificate and costs

Combining the six shapes gives

 U2=3q^2*c^2*(binom(n,4)+binom(n,3))/r^2
       +2q*c*n^2*T/(rW)+2q*c*n^2*T^2/(rW)
 >=E binom(C,2).

Pointwise1[C>0]>=C-binom(C,2),so repeat probability>=m-U2. Exact full-order rationals at the stated K give>0.3900031430168671 before the separate failure union. Subtract fingerprint2^-55,table2^-66,sampler2^-970,identity2^-100,dictionary2^-100,singular2^-180 to retain success>0.39000314301686>0.39. The actual algorithm detects the first repeat at ANY lag. Extending the ideal path for C is only analysis:P[first repeat<=K]>=P[C>0]. Actual hashing is coupled only through that first repeat;no actual post-repeat quotient-map behavior is assumed.

The exact interval and dictionary certificates in probability.md/storage.md are checked at this K. Native ordinary scalar work is12288,with one C20000 reference per transition and complete coldD540*2^119. Thus work<=(32288K+540*2^119)/20000. Exact(32288K+540*2^119)^10000<20000^10000*2^1263707 supports reduced proposed1263707/10000=126.3707,above diagnostic126.37064143144556. Full parameters,threshold,reference prices,metric,tracks and initialization/memory conditions are unchanged. Material conditions remain full-standard-BLAKE3 H and finite banked-storage admissibility. This is conditional analytic evidence,not qualification or official scoring.
