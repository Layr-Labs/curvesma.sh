# Finite repeat bound inside the sole disclosed hash model

Write N=r, w=2^40, rho=1/32768, T=271, eta=1/65536, R=ceil(sqrt(N)),
L=ceil(1031R/1000), s=ceil(L/62), K=L+s. This is the inherited finite proof,
not an assertion that the adding walk is a uniform random mapping. H is disclosed
in proof.md; fixed SHA-512 is not mathematically certified by this calculation.

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

On the independent-index surrogate the group states are independent of the five
compression keys. Conditional on any surrogate, every fixed distinct-state pair
has compression-collision probability<=(38/q)^5. With at mostK+1<2^129 states,
union loss<2^257*(2^-54)^5=2^-13, using q>2^60,38<64. Before first group repeat
or compression collision, the ideal hash walk can be coupled exactly to that
surrogate: a new compressed input has a fresh uniform40-bit label. This avoids
invalid conditioning on an already key-dependent trajectory. After the first
repeat the two chains need not agree; only first-repeat distribution transfers.

Checkpoint detection is the deterministic tail/cycle lemma in proof.md. Given
jump group values and the trajectory, v_j stay independent uniform and the first
relation's D=sum n_j*v_j is uniform since0<sum n_j<=K<N. Its zero probability
is1/N even conditional on meeting. Main scalar/key rejection failures sum below
2^-980 under H. Accepted values are uniform independent of their rejection
counts. Subtract compression2^-13, bad table2^-100, sampling2^-980 and singular
relation1/N from0.39015; the result is strictly above0.39. Identity targets succeed.

## Exact public checks

Use integer N from protected target-profile.json and rational arithmetic with
R=isqrt(N-1)+1,L=(1031R+999)//1000,s=(L+61)//62,K=L+s,
D=(L-T+1)*(L-T+2)//2,T=271,eta=1/65536. The exact inequalities are:

 (3m-U)/2 > 39015/100000;
 39015/100000-1/8192-1/2^100-1/2^980-1/N > 39/100;
 (33088K+2^90)^500 < 20000^500*2^64111;
 K<N and K+1<2^129;
 4N*3^256*2^100<8^256.

They are public finite arithmetic, not production challenge execution, empirical
success rates or reviewer confidence. All hold for the protected production N.
The sole heuristic is H's bridge to actual fixed SHA-512, disclosed in proof.md.
