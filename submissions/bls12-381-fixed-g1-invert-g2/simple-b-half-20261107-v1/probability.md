# Finite repeat bound inside the sole disclosed hash model

Write N=r, w=2^40, rho=1/32768, T=271, eta=1/65536, R=ceil(sqrt(N)),
L=ceil(1031R/1000), s=ceil(L/126), K=L+s. This is the inherited finite proof,
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

The tables for simple tabulation are independent uniform320-bit values under the expanded disclosed H. On the independent-index surrogate, all group states are independent of these tables. For every fixed distinct pair of injective keys, choose a differing digit and condition on all entries except its entry from the first key. Since each position has its own independent table and the two digits at that position are different, the one unconditioned entry appears once in their XOR difference. Thus Pr[U(X)=U(Y)]=2^-320. Union over at mostK+1<2^129 surrogate states gives loss<2^257*2^-320=2^-63<2^-13. This uses pairwise universality only, not a claim that fingerprints across all states are independent. Couple the actual walk and independent-label surrogate only until first group repeat or distinct-state fingerprint collision. Until then every fingerprint is new, so each WALK hash output provides a fresh uniform label under H. The inherited conservative2^-13 loss is retained below, and no adaptive union bound on the actual table-dependent walk is asserted.

Checkpoint detection is the deterministic tail/cycle lemma in proof.md. Given
jump group values and the trajectory, v_j stay independent uniform and the first
relation's D=sum n_j*v_j is uniform since0<sum n_j<=K<N. Its zero probability
is1/N even conditional on meeting. Main scalar/key rejection failures sum below
2^-980 under H. Accepted values are uniform independent of their rejection
counts. Subtract compression2^-13, bad table2^-100, sampling2^-980 and singular
relation1/N from0.39015; the result is strictly above0.39. Identity targets succeed.

## Exact public checks

Use integer N from protected target-profile.json and rational arithmetic with
R=isqrt(N-1)+1,L=(1031R+999)//1000,s=(L+125)//126,K=L+s,
D=(L-T+1)*(L-T+2)//2,T=271,eta=1/65536. The exact inequalities are:

 (3m-U)/2 > 39015/100000;
 39015/100000-1/8192-1/2^100-1/2^980-1/N > 39/100;
 (22432K+2^90)^20 < 20000^20*2^2553;
 K<N and K+1<2^129;
 4N*3^256*2^100<8^256.

They are public finite arithmetic, not production challenge execution, empirical
success rates or reviewer confidence. All hold for the protected production N.
The sole heuristic is H's bridge to actual fixed BLAKE3, disclosed in proof.md.

The current exact key is36canonical b words plus a constantzero word. proof.md proves b alone injective on odd-order GT. Its constantT72[0] fingerprint term cancels between distinct points,so the same320-bit pairwise-universality and surrogate coupling argument applies. Actualfingerprints/labels change from the sign-based source,buttheir distribution underH andallfinite moments above are unchanged. No added prefix cap,orbit model,nonbacktracking hypothesis orbanked memory is used.
