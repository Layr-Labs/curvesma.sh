# Finite repeat bound inside the sole disclosed hash model

Write N=r, w=2^40, rho=1/30000, T=271, eta=1/65536, R=ceil(sqrt(N)),
L=ceil(995R/1000), s=ceil(L/1022), K=L+s. This is the inherited finite proof,
not an assertion that the adding walk is a uniform random mapping. H is disclosed
in proof.md; fixed BLAKE3 is not mathematically certified by this calculation.

## Good table

For nontrivial cyclic-group character chi, lambda_chi=(1/w)sum_j chi(M_j).
Independent uniform jumps give mean zero real and imaginary summands in[-1,1].
If |lambda_chi|>rho, some component exceeds rho/sqrt(2) in magnitude.
Hoeffding and union over both components and N-1 characters give failure
 beta <=4N exp(-w*rho²/4)<4N exp(-305)<2^-174.
A finite certificate is4N*3^305*2^174<8^305, using e>8/3 and
w>305*4*30000². Duplicate jumps are permitted.
The attack does not test or select good tables; every table and failure is paid.

Condition on a good table. For the independent-uniform-index addition chain let
p_d(x) be d-step increment distribution and epsilon_d(x)=p_d(x)-1/N. Fourier
inversion and Parseval give sup|epsilon_d|<=rho^d and unweighted l2 norm<=rho^d,
including d=0. This is valid for complex eigenvalues/nonreversible chains.
Since rho^T<=eta/N, returns after d>=T have probability in[(1-eta)/N,(1+eta)/N].

## All interval correlations

Among Y_0,...,Y_L count C=sum I_(i,j) for i<j,j-i>=T with Y_i=Y_j.
Let D=(L-T+1)(L-T+2)/2. Then m=(1-eta)D/N<=E C. Equality events depend on
independent increment intervals, so the starting point does not matter.
We bound ordered indicator pairs for E C². Identical intervals contribute
at most(1+eta)D/N. A baseline(1+3eta)D²/N² covers the uniform terms of all
distinct ordered intervals; the following excess terms cover their dependence.

Disjoint intervals, including those touching at an endpoint, are independent,
with joint probability<=((1+eta)/N)², already within that baseline.
For nested intervals, write consecutive lengths a>=0,b>=T,c>=0,a+c>0.
The inner return and remaining outer increments are independent, so joint
probability=p_b(1)*p_(a+c)(1)
 <=(1+eta)/N²+(1+eta)rho^(a+c)/N.
There are at most(L+1)² choices of start and b. The geometric sum over a,c
with a+c>0 is(1-rho)^(-2)-1. Both indicator orders give
nested excess<=2(1+eta)((1-rho)^(-2)-1)(L+1)²/N. Cases a=0 or c=0 include
shared-start/shared-end intervals; identical intervals were separated.

For crossing intervals, successive positive gaps a,b,c satisfy a+b>=T,b+c>=T.
Joint probability=sum_x p_a(x)*p_b(x^-1)*p_c(x). Expand p as uniform+epsilon.
Linear terms vanish. Cauchy-Schwarz bounds quadratic terms and sup/l2 bounds
bound the cubic term, giving
 joint<=1/N²+(rho^(a+b)+rho^(a+c)+rho^(b+c))/N+rho^(a+b+c).
The rho^(a+c) sum over positive a,c is(rho/(1-rho))²; free start
and b, with both indicator orders, give excess<=2(rho/(1-rho))²(L+1)²/N.
For the other two quadratic terms,
 sum_(a,b>=1,a+b>=T) rho^(a+b)<=2T*2^-T.
Both choices and indicator orders with free start/remaining gap give
<=8T*eta*(L+1)²/N². For the cubic constrained sum, dominate rho by1/2:
b>=T-1 contributes<=4*2^-T; b=1..T-2 with a,c>=T-b contributes<2*2^-T.
After free starting position and both orders the cubic excess is
<12eta*(L+1)/N. All distinct endpoint orders are disjoint, nested or crossing;
there are no omitted correlations or assumptions of collision independence.

Thus E C²<=U, where
 U=(1+3eta)D²/N²+(1+eta)D/N
   +[2(1+eta)((1-rho)^(-2)-1)+2(rho/(1-rho))²](L+1)²/N
   +8T eta (L+1)²/N²+12eta(L+1)/N.
The first two Bonferroni terms satisfy E[C-choose(C,2)] >=(3m-U)/2. At the shorter L=ceil(995R/1000) this alone is insufficient. We now bound the next two factorial moments explicitly.

## Third factorial moment: lower bound

Choose six distinct ordered time indices, with each consecutive gap at least T. Their count is binom(L-5(T-1)+1,6): subtract i(T-1) from the i-th index. Each of the15 perfect matchings gives three distinct eligible equality events. Every intervening transition has pointwise density at least(1-eta)/N. Conditional on the first endpoint value, equality leaves two free values, so the probability is at least(1-eta)^5/N^3. These are distinct terms in choose(C,3). Thus

 V3 =15(1-eta)^5*binom(L-5(T-1)+1,6)/N^3 <= E choose(C,3).

No assumption about independence of equality events is used.

## Fourth factorial moment: equality graphs and Fourier flows

A term of choose(C,4) is a set of four distinct eligible interval edges. Sort its ell distinct time endpoints,4<=ell<=8, and put a simple graph G with four edges on these ordered vertices. Every vertex is incident. Let q be its number of connected components and rank=ell-q. Equality identifies all endpoint values inside each component; fix the first endpoint value and sum over the other q-1 class values. The ell-1 time gaps have transition kernels p_d=1/N+epsilon_d.

Expand this product. For a subset F of the ell-1 gap positions, use epsilon on F and1/N elsewhere. Contract vertices by the equality classes of G. Let f=|F|,v the number of incident contracted vertices,c their number of nonempty connected components, and b=f-v+c (the cycle-space dimension, including loops). If the contracted F graph has a degree-one vertex, its term is zero: summing its epsilon kernel gives zero. Loops count twice toward degree.

Otherwise Fourier expansion gives the exact useful bound

 |term_F| <= N^(-rank+b)*rho^(sum_(i in F) d_i).

Indeed epsilon_d(x)=N^-1*sum_(chi!=1)lambda_chi^d*chi(x). Summing class values enforces zero character divergence at each vertex. There are at most N^(f-v+c) allowed character flows; every factor has magnitude<=rho^d. With q-1 free class values and ell-1-f uniform factors this is exactly the displayed power. This argument works for complex eigenvalues and arbitrary orientations, parallel edges and loops. No reversibility is assumed.

For empty F take v=c=b=0. The endpoint variable ell is distinct from the first-moment lower bound m above.

Let A=1002/1000 and t=rho/(1-rho)=1/29999. At protected N and L=ceil(995R/1000), (L+1)/sqrt(N)<A. For fixed F-gap lengths, the start and other gaps have at most(L+1)^(ell-f) choices. Summing positive F gaps gives t^f. The net power of N after this count is

 E =-rank+b+(ell-f)/2.

The finite graph certificate in templates.md establishes the following, exhaustively over G and F after removing degree-one terms. There are22246 retained graph/subset terms, fewer than32768.

* F empty with ell=8 consists of105 perfect matchings and contributes at most105*binom(L+1,8)/N^4. Every other empty term has E<=-1/2.
* If no original G edge has all its intervening gap positions in F (the free case), then E<=0. The counts for E=0, aggregated by f=1..6, are90,139,122,121,72,24. Other free terms have E<=-1/2.
* Otherwise some original equality edge is wholly inside F gaps. Eligibility forces sum_(i in F)d_i>=T. Every such term has E<=2.

Negative-exponent terms, including empty shared-endpoint terms, contribute at most32768*A^8/sqrt(N)<2^-110 (using N>2^254 and A^8<2). Long forced terms contribute less than32768*A^8*N^2*2^-1897<2^-1300: split rho^D into rho^(D/2)rho^(D/2), use D>=T=271 and rho^(T/2)<2^-1897, and sum the remaining positive gaps with sqrt(rho)<1/128. Their geometric factors are at most1. Both error classes together are below2^-100.

## Centering the90 single-loop terms over the jump table

The90 free E=0,f=1 terms all have ell=7,q=3,rank=4 and a single contracted loop. For fixed endpoint times and gap d their exact expansion term is epsilon_d(0)/N^3. The graph certificate in templates.md identifies these terms explicitly. We now average over random jump tables, conditional on the good-table event, instead of taking a worst-case rho^d bound for them. All other bounds remain uniform on good tables.

For every1<=d<=L<N, E_table p_d(0)=1/N exactly. To see this, first fix the d iid uniform jump indices. Their multiplicities n_j are nonnegative, sum to d, and some n_j is nonzero modulo prime N. Independent uniform jump exponents make sum n_j*M_j uniform. Averaging indices preserves this. Hence E_table epsilon_d(0)=0, with no extra hash or relation heuristic.

Since |epsilon_d(0)|<=1 for every table, conditional centering gives

 E[epsilon_d(0) | good] =-E[1_bad*epsilon_d(0)]/Pr[good] <= beta/(1-beta).

The interval eligibility weights and endpoint counts are table-independent. Summing the90 templates over at most(L+1)^7 ordered endpoint assignments bounds their entire conditional contribution by

 90*beta/(1-beta)*(L+1)^7/N^3 <2^-36.

Here beta<2^-174,(L+1)/sqrt(N)<A,A^7<2 and sqrt(N)<2^128. This is an expectation bound over good tables, not a false uniform bound for each individual good table. No rare event is multiplied by an unbounded count without its full count: the displayed(L+1)^7/N^3 factor is explicitly included.

Thus a valid conditional fourth-moment upper bound is

 U4c =105*binom(L+1,8)/N^4
       +A^8*(139t^2+122t^3+121t^4+72t^5+24t^6)
       +2^-36+2^-100 >= E[choose(C,4) | good].

The remaining free/shared-endpoint terms are still retained. Only the exactly centered90 single-loop terms replace their prior90t allowance. First/second/third moment bounds apply uniformly on every good table, so they also apply after averaging good tables. This adds no new heuristic.

For every integer C>=0,

 1_(C>0) >= C-choose(C,2)+choose(C,3)-choose(C,4).

For C>=1 the right side is1-choose(C-1,4), so this is pointwise. Therefore

 Pr[C>0 | good] >= (3m-U)/2+V3-U4c >3901/10000.

The exact diagnostic value of this lower bound is0.39011968016756776. This is finite rational arithmetic at the full production order, not a Poisson approximation or empirical extrapolation.

## Compression and relation extraction

For full-K adaptive coupling use the independent group-state-label surrogate J from proof.md: labels are uniform independent on first visits and cached on revisits. Before first repeat it is the independent-index chain above; its full path is independent of tabulation entries. Each fixed unequal pair has full-fingerprint collision probability2^-320, so union over at mostK+1<2^129 states loses <2^-63. Couple every new actual fingerprint input to J(X), retaining equality on revisits. Without a distinct-state fingerprint collision the full K-step actual path agrees, including its cycle. Sorting on prefixes never changes group transitions.

For the64-bit cached prefix, fixed unequal keys collide with probability2^-64. Condition on the full table-independent surrogate path, scheduled checkpoint set and ideal group-equality stopping time. At most1023K unequal query/checkpoint pairs can cause false full comparisons. Prefix sorting/scanning chooses a subset; no independent-prefix assumption is needed. Hence E F<=1023K/2^64+1 including one eventual true comparison. With B_cap=floor(K/256)>=K/512 and K>2^64, Pr[F>B_cap]<=512*(1023/2^64+1/K)<2^-43. This is a union/expectation/Markov argument before conditioning on no full-fingerprint collisions; losses are subtracted separately. It is not an adaptive union bound on the actual prefix-dependent trajectory. Every full comparison decrements the global cap and costs<=256; total extra instructions<=256B_cap<=K. The extra failure loss2^-43 is part of the success claim.

Checkpoint detection is the deterministic tail/cycle lemma in proof.md. Given
jump group values and the trajectory, v_j stay independent uniform and the first
relation's D=sum n_j*v_j is uniform since0<sum n_j<=K<N. Its zero probability
is1/N even conditional on meeting. Main scalar/key rejection failures sum below
2^-980 under H. Accepted values are uniform independent of their rejection
counts. Subtract compression2^-63, bad table2^-174, sampling2^-980, prefix-budget2^-43 and singular
relation1/N from0.3901; the result is strictly above0.39. Identity targets succeed.

## Exact public checks

Set N to the protected subgroup order,R=isqrt(N-1)+1,L=(995R+999)//1000,s=(L+1021)//1022,K=L+s,T=271,eta=1/65536,rho=1/30000,A=1002/1000,t=1/29999. Compute m,U,V3,U4c as exact rationals. The checks are:

 (3m-U)/2+V3-U4c >3901/10000;
 3901/10000-1/2^63-1/2^174-1/2^980-1/N-1/2^43 >39/100;
 (21936K+2^90)^1000 <20000^1000*2^127557;
 (L+1)^2*1000^2 <1002^2*N;
 N>2^254,N<2^255,K<N,K+1<2^129,1023s>K;
 4N*3^305*2^174<8^305;
 2^40>305*4*30000^2;
 90*(L+1)^7*2^-174/(N^3*(1-2^-174)) <2^-36.

The diagnostic work logarithm is127.55602456233001; the conservative submitted bound is127557/1000=127.557. The finite success lower bound after all stated losses exceeds0.3901196. The fingerprint subtraction uses the proved full320-bit union bound<2^-63, rather than the older unnecessarily loose2^-13. This does not change the algorithm, hash model or coupling. Public arithmetic does not establish fixed-hash H, production execution or qualification.
