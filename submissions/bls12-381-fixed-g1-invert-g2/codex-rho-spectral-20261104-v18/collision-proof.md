# Finite collision bound in the ideal-hash model

Let N be the protected prime subgroup order, w=65536, and M_1,...,M_w independent uniform elements of GT. The implementation generates these from uniformly idealized initialization scalars; no test for a good table or discrete-log oracle is performed. This note proves a finite collision bound in that idealization. It does not assert that fixed SHA-256 is a mathematically random function.

## Good-table probability

For every nontrivial character chi of the cyclic group, define lambda_chi=(1/w)*sum_j chi(M_j). Its real and imaginary summands lie in[-1,1],have mean zero,and are independent over j. If |lambda_chi|>1/2,one component has magnitude>1/(2sqrt(2)). Hoeffding's inequality on each component and a union bound give probability<=4exp(-w/16). Union over N-1 characters gives bad-table probability<=4Nexp(-4096)<2^-5000 for the protected N<2^256. For an exact check, e>8/3 implies this is below4N*(3/8)^4096; the claimed dyadic bound is the integer inequality4N*3^4096*2^5000<8^4096. The experiment checks this public constant inequality; it does not sample all production tables.

Condition henceforth on a good table, so every nontrivial eigenvalue of its uniform-addition Cayley transition matrix has modulus<=rho=1/2. Duplication among the sampled M_j is allowed and does not invalidate Hoeffding or the transition weights. The attack pays for all samples and accepts bad tables too, counting their possible failures.

## Markov surrogate and coupling

Consider independent uniform indices J_t in[1,w] and the Markov walk Y_(t+1)=Y_t*M_(J_t). A random-function partition assigns independent uniform indices to previously unvisited elements. Therefore its deterministic iterated walk and this independent-index walk can be coupled exactly until the first repeated group element. They need not agree after that repeat. Only the distribution of the first repeat is transferred below. Floyd runs on the actual deterministic map, not on a fresh-index walk.

Write p_d(x) for the d-step increment distribution and epsilon_d(x)=p_d(x)-1/N. Cyclic convolution is diagonalized by the characters,including complex eigenvalues for this possibly nonreversible chain. Fourier inversion and Parseval imply

sup_x |epsilon_d(x)|<=rho^d, and (sum_x |epsilon_d(x)|^2)^(1/2)<=rho^d.

These also hold at d=0. Set L=5ceil(sqrt(N)), T=ceil(log2(N))+16=271 and eta=2^-16. Since rho^T<=eta/N, for d>=T the return probability p_d(1) lies between(1-eta)/N and(1+eta)/N.

## Long-pair second moment

Among states Y_0,...,Y_L,count collisions C=sum I_(i,j) over i<j and j-i>=T. The number of such pairs is D=(L-T+1)(L-T+2)/2,so E C>=(1-eta)D/N. All equality events depend only on independent increments over their intervals; the starting group value does not affect the argument.

We bound ordered pairs of indicators in C^2. Identical intervals contribute at most(1+eta)D/N. For distinct intervals,include a baseline(1+3eta)/N^2 per ordered pair,whose sum is at most(1+3eta)D^2/N^2. The following bounds account for every excess above that baseline.

Disjoint increment intervals,including intervals touching at one endpoint,are independent and have joint probability<=((1+eta)/N)^2,already covered. Intervals sharing their first or last endpoint belong to the nested case below; identical intervals have already been separated.

For nested intervals with consecutive segment lengths a>=0,b>=T,c>=0 and a+c>0,the inner return and the outer residual return depend on disjoint increments. Their joint probability is p_b(1)*p_(a+c)(1),at most(1+eta)/N^2+(1+eta)*rho^(a+c)/N. There are at most(L+1)^2 choices for starting position and b. Sum over a,c with a+c>0 gives sum rho^(a+c)=3. Both orders of the two indicators give total excess<=6(1+eta)(L+1)^2/N. This includes a=0 or c=0,so shared-endpoint nested intervals are counted.

For crossing intervals,the four ordered endpoints have positive successive gaps a,b,c with a+b>=T,b+c>=T. Independence of these three increment blocks gives joint probability sum_x p_a(x)*p_b(x^-1)*p_c(x). Expand each p as uniform+epsilon. Linear terms vanish. Cauchy-Schwarz on quadratic terms and sup/l2 bounds on the cubic term give

joint <=1/N^2+(rho^(a+b)+rho^(a+c)+rho^(b+c))/N+rho^(a+b+c).

Sum the rho^(a+c) term over positive a,c:the sum is1. Free starting position and b,plus the two indicator orders,give excess<=2(L+1)^2/N. For the other two quadratic terms,sum_(a,b>=1,a+b>=T)rho^(a+b)=2T*rho^T. The free remaining segment and start,with both orientations and both terms,give excess<=8T*eta*(L+1)^2/N^2.

For the cubic term,its constrained sum over a,b,c is below6rho^T. Indeed,b>=T-1 contributes4rho^T; for1<=b<=T-2,the conditions require a,c>=T-b and the geometric sum contributes below2rho^T. After both indicator orders and the starting position,the contribution is below12eta*(L+1)/N. No interval pattern is omitted:ordered distinct intervals are disjoint,nested,or crossing.

Combining terms proves

E C^2 <= (1+3eta)D^2/N^2 +(1+eta)D/N +(8+6eta)(L+1)^2/N +8T*eta*(L+1)^2/N^2 +12eta*(L+1)/N.

This is a finite bound,not an assumption of independent collision indicators. In particular the short-return correlations of an adding walk are explicitly included.

## Production numerical bounds

Here N>2^254. The definition of L and T gives12.499<D/N<12.501 and(L+1)^2/N<25.001. For example,D/N>12.5-5T/sqrt(N)>12.499,and L+1<=5sqrt(N)+6 gives the upper bounds. The last two terms in the displayed second moment sum to less than1 by these same bounds and N>2^254.

Using eta=1/65536,((1-eta)*12.499)^2>156,and

(1+3eta)*12.501^2 +(1+eta)*12.501 +(8+6eta)*25.001 +1 <370.

Thus(E C)^2>156 and E C^2<370. Cauchy-Schwarz on C times its positivity indicator yields P(C>0)>156/370>0.42. Any such collision implies the first repeat occurred by L. Coupling until the first repeat transfers that bound to the ideal-partition deterministic adding walk. For a deterministic map with tail mu and cycle lambda,Floyd's first positive meeting time is a multiple of lambda no greater than mu+lambda,the first-repeat index. Hence it also meets by L with probability>0.42 conditioned on a good table. The analysis never uses independence after a repeated state.

## Useful relation and total success

Let z=h^a. Conditional on the group values M_j,each v_j is still uniform and independent:for every v_j there is exactly one u_j giving that M_j=h^u_j*z^v_j. The trajectory,partition and initial group value do not inspect these coefficient representations. Conditioning on the good-table property,or on a particular meeting,does not change this fact.

At a Floyd meeting t<=L<N,the difference B_fast-B_slow equals sum_j n_j v_j,where n_j is the nonnegative count of jump j between trajectory positions t and2t. Their sum is t>0 and<N,so some n_j is nonzero modulo the prime N. The linear form is therefore uniform modulo N and vanishes with probability exactly1/N. Otherwise the invariant X=h^A*z^B yields a=(A_slow-A_fast)/(B_fast-B_slow) and the algorithm can verify and output[a]Q0. This is a conditional distribution proof,not an extra nonsingular-relation heuristic. No automorphism/negation normalization is applied to the walk.

With2w+1 initialization scalars and at most1024 independent255-bit rejection attempts each,the sampling-failure probability under the idealized hash is below(2w+1)*2^-1024<2^-1000,since N>2^254. Accepted values remain independent uniform N-residues conditioned on success; the rejection count is independent of the accepted value. Domain separation makes the partition independent of that initialization in the idealized analysis.

Subtracting bad-table probability,sampling failure and singular-relation probability gives idealized success>0.42-2^-5000-2^-1000-1/N>0.41. The submitted success lower bound is the more conservative0.4. The sole remaining computational heuristic is transferring this random-function analysis to the fully computed,fixed SHA-256 construction over uniform challenges. No unconditional success theorem for fixed SHA is asserted. All work is bounded and charged even if that heuristic fails.
