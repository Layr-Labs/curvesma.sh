# Finite repeat bound inside the sole disclosed hash model

Write N=r, w=2^40, rho=1/32768, T=271, eta=1/65536, R=ceil(sqrt(N)),
L=ceil(1031R/1000), s=ceil(L/1022), K=L+s. This is the inherited finite proof,
not an assertion that the adding walk is a uniform random mapping. H is disclosed
in proof.md; fixed BLAKE3 is not mathematically certified by this calculation.

## Good table

For nontrivial cyclic-group character chi, lambda_chi=(1/w)sum_j chi(M_j).
Independent uniform jumps give mean zero real and imaginary summands in[-1,1].
If |lambda_chi|>rho, some component exceeds rho/sqrt(2) in magnitude.
Hoeffding and union over both components and N-1 characters give failure
<=4N exp(-w*rho²/4)=4N exp(-256)<2^-100. A finite certificate is
4N*3^256*2^100<8^256, using e>8/3. Duplicate jumps are permitted.
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
with a+c>0 is(1-rho)^(-2)-1=65535/1073676289. Both indicator orders give
nested excess<=(131070/1073676289)(1+eta)(L+1)²/N. Cases a=0 or c=0 include
shared-start/shared-end intervals; identical intervals were separated.

For crossing intervals, successive positive gaps a,b,c satisfy a+b>=T,b+c>=T.
Joint probability=sum_x p_a(x)*p_b(x^-1)*p_c(x). Expand p as uniform+epsilon.
Linear terms vanish. Cauchy-Schwarz bounds quadratic terms and sup/l2 bounds
bound the cubic term, giving
 joint<=1/N²+(rho^(a+b)+rho^(a+c)+rho^(b+c))/N+rho^(a+b+c).
The rho^(a+c) sum over positive a,c is(rho/(1-rho))²=1/1073676289; free start
and b, with both indicator orders, give excess<=(2/1073676289)(L+1)²/N.
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
   +((131072+131070eta)/1073676289)(L+1)²/N
   +8T eta (L+1)²/N²+12eta(L+1)/N.
For every integer c>=0, 1_(c>0)>=(3c-c²)/2: it equals0,1,1 at c=0,1,2,
zero at3, nonpositive afterwards. Hence Pr[C>0]>=(3m-U)/2>39015/100000.
This finite two-term Bonferroni bound uses no higher-moment or mapping heuristic.

## Compression and relation extraction

For full-K adaptive coupling use the independent group-state-label surrogate J from proof.md: labels are uniform independent on first visits and cached on revisits. Before first repeat it is the independent-index chain above; its full path is independent of tabulation entries. Each fixed unequal pair has full-fingerprint collision probability2^-320, so union over at mostK+1<2^129 states loses <2^-63<2^-13. Couple every new actual fingerprint input to J(X), retaining equality on revisits. Without a distinct-state fingerprint collision the full K-step actual path agrees, including its cycle. Sorting on prefixes never changes group transitions.

For the64-bit cached prefix, fixed unequal keys collide with probability2^-64. Condition on the full table-independent surrogate path, scheduled checkpoint set and ideal group-equality stopping time. At most1023K unequal query/checkpoint pairs can cause false full comparisons. Prefix sorting/scanning chooses a subset; no independent-prefix assumption is needed. Hence E F<=1023K/2^64+1 including one eventual true comparison. With B_cap=floor(K/256)>=K/512 and K>2^64, Pr[F>B_cap]<=512*(1023/2^64+1/K)<2^-43. This is a union/expectation/Markov argument before conditioning on no full-fingerprint collisions; losses are subtracted separately. It is not an adaptive union bound on the actual prefix-dependent trajectory. Every full comparison decrements the global cap and costs<=256; total extra instructions<=256B_cap<=K. The extra failure loss2^-43 is part of the success claim.

Checkpoint detection is the deterministic tail/cycle lemma in proof.md. Given
jump group values and the trajectory, v_j stay independent uniform and the first
relation's D=sum n_j*v_j is uniform since0<sum n_j<=K<N. Its zero probability
is1/N even conditional on meeting. Main scalar/key rejection failures sum below
2^-980 under H. Accepted values are uniform independent of their rejection
counts. Subtract compression2^-13, bad table2^-100, sampling2^-980, prefix-budget2^-43 and singular
relation1/N from0.39015; the result is strictly above0.39. Identity targets succeed.

## Exact public checks

Use integer N from protected target-profile.json and rational arithmetic with
R=isqrt(N-1)+1,L=(1031R+999)//1000,s=(L+1021)//1022,K=L+s,
D=(L-T+1)*(L-T+2)//2,T=271,eta=1/65536. The exact inequalities are:

 (3m-U)/2 > 39015/100000;
 39015/100000-1/8192-1/2^100-1/2^980-1/N-1/2^43 > 39/100;
 (22112K+2^90)^50 < 20000^50*2^6381;
 K<N and K+1<2^129;
 4N*3^256*2^100<8^256.

They are public finite arithmetic, not production challenge execution, empirical
success rates or reviewer confidence. All hold for the protected production N.
The sole heuristic is H's bridge to actual fixed BLAKE3, disclosed in proof.md.
