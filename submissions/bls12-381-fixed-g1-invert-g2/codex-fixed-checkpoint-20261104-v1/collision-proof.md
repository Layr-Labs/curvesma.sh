# Finite collision bound in the ideal-hash model

Let N be the protected prime subgroup order, w=16777216, and M_1,...,M_w independent uniform elements of GT. The implementation generates these from uniformly idealized initialization scalars; no test for a good table or discrete-log oracle is performed. This note proves a finite collision bound in that idealization. It does not assert that fixed SHA-256 is a mathematically random function.

## Good-table probability

For every nontrivial character chi of the cyclic group, define lambda_chi=(1/w)*sum_j chi(M_j). Its real and imaginary summands lie in[-1,1],have mean zero,and are independent over j. If |lambda_chi|>1/128,one component has magnitude>1/(128sqrt(2)). Hoeffding's inequality on each component and a union bound give probability<=4exp(-w/65536). Union over N-1 characters gives bad-table probability<=4Nexp(-256)<2^-100 for the protected N<2^256. For an exact check, e>8/3 implies this is below4N*(3/8)^256; the claimed dyadic bound is the integer inequality4N*3^256*2^100<8^256. This is an analytic public-integer inequality; no experimental result is claimed.

Condition henceforth on a good table, so every nontrivial eigenvalue of its uniform-addition Cayley transition matrix has modulus<=rho=1/128. Duplication among the sampled M_j is allowed and does not invalidate Hoeffding or the transition weights. The attack pays for all samples and accepts bad tables too, counting their possible failures.

## Markov surrogate and coupling

Consider independent uniform indices J_t in[1,w] and the Markov walk Y_(t+1)=Y_t*M_(J_t). A random-function partition assigns independent uniform indices to previously unvisited elements. Therefore its deterministic iterated walk and this independent-index walk can be coupled exactly until the first repeated group element. They need not agree after that repeat. Only the distribution of the first repeat is transferred below. Fixed checkpoint detection runs on the actual deterministic map,not on a fresh-index walk.

Write p_d(x) for the d-step increment distribution and epsilon_d(x)=p_d(x)-1/N. Cyclic convolution is diagonalized by the characters,including complex eigenvalues for this possibly nonreversible chain. Fourier inversion and Parseval imply

sup_x |epsilon_d(x)|<=rho^d, and (sum_x |epsilon_d(x)|^2)^(1/2)<=rho^d.

These also hold at d=0. Set R=ceil(sqrt(N)), L=ceil(119R/100), T=ceil(log2(N))+16=271 and eta=2^-16. Since rho^T<=eta/N, for d>=T the return probability p_d(1) lies between(1-eta)/N and(1+eta)/N.

## Long-pair second moment

Among states Y_0,...,Y_L,count collisions C=sum I_(i,j) over i<j and j-i>=T. The number of such pairs is D=(L-T+1)(L-T+2)/2,so E C>=(1-eta)D/N. All equality events depend only on independent increments over their intervals; the starting group value does not affect the argument.

We bound ordered pairs of indicators in C^2. Identical intervals contribute at most(1+eta)D/N. For distinct intervals,include a baseline(1+3eta)/N^2 per ordered pair,whose sum is at most(1+3eta)D^2/N^2. The following bounds account for every excess above that baseline.

Disjoint increment intervals,including intervals touching at one endpoint,are independent and have joint probability<=((1+eta)/N)^2,already covered. Intervals sharing their first or last endpoint belong to the nested case below; identical intervals have already been separated.

For nested intervals with consecutive segment lengths a>=0,b>=T,c>=0 and a+c>0,the inner return and the outer residual return depend on disjoint increments. Their joint probability is p_b(1)*p_(a+c)(1),at most(1+eta)/N^2+(1+eta)*rho^(a+c)/N. There are at most(L+1)^2 choices for starting position and b. Sum over a,c with a+c>0 gives sum rho^(a+c)=1/(1-rho)^2-1=255/16129. Both orders of the two indicators give total excess<=(510/16129)(1+eta)(L+1)^2/N. This includes a=0 or c=0,so shared-endpoint nested intervals are counted.

For crossing intervals,the four ordered endpoints have positive successive gaps a,b,c with a+b>=T,b+c>=T. Independence of these three increment blocks gives joint probability sum_x p_a(x)*p_b(x^-1)*p_c(x). Expand each p as uniform+epsilon. Linear terms vanish. Cauchy-Schwarz on quadratic terms and sup/l2 bounds on the cubic term give

joint <=1/N^2+(rho^(a+b)+rho^(a+c)+rho^(b+c))/N+rho^(a+b+c).

Sum the rho^(a+c) term over positive a,c:the sum is(rho/(1-rho))^2=1/16129. Free starting position and b,plus the two indicator orders,give excess<=(2/16129)(L+1)^2/N. For the other two quadratic terms,sum_(a,b>=1,a+b>=T)rho^(a+b)<=2T*2^-T. The free remaining segment and start,with both orientations and both terms,give excess<=8T*eta*(L+1)^2/N^2.

For the cubic term,its constrained sum over a,b,c is below6*2^-T (dominate termwise by rho=1/2). Indeed,b>=T-1 contributes at most4*2^-T; for1<=b<=T-2,the conditions require a,c>=T-b and the geometric sum contributes below2*2^-T. After both indicator orders and the starting position,the contribution is below12eta*(L+1)/N. No interval pattern is omitted:ordered distinct intervals are disjoint,nested,or crossing.

Combining terms proves

E C^2 <= (1+3eta)D^2/N^2 +(1+eta)D/N +((512+510eta)/16129)(L+1)^2/N +8T*eta*(L+1)^2/N^2 +12eta*(L+1)/N.

This is a finite bound,not an assumption of independent collision indicators. In particular the short-return correlations of an adding walk are explicitly included.

## Production numerical bounds

Here N>2^254, T=271, eta=1/65536, L=ceil(119ceil(sqrt(N))/100), and D=(L-T+1)(L-T+2)/2. Let

m=(1-eta)D/N,
U=(1+3eta)D^2/N^2+(1+eta)D/N+((512+510eta)/16129)(L+1)^2/N+8T eta (L+1)^2/N^2+12eta(L+1)/N.

The exact rational inequality m^2/U>799/2000=0.3995 holds for the protected public N. A simple coarse check uses 0.708049<D/N<0.708051, (L+1)^2/N<1.416101, and the final two terms together less than0.000001. Those give m^2>0.50131 and U<1.2544, with0.50131/1.2544>0.3995. Cauchy-Schwarz yields P(C>0)>0.3995. Thus the first repeat is by L with that probability on a good table. The deterministic fixed checkpoint lemma detects it by K=2L-1. All error-tail bounds use rho<=1/2 and2^-T<=eta/N. No independence after a repeated state is used.

## Useful relation and total success

Let z=h^a. Conditional on the group values M_j,each v_j is still uniform and independent:for every v_j there is exactly one u_j giving that M_j=h^u_j*z^v_j. The trajectory,partition and initial group value do not inspect these coefficient representations. Conditioning on the good-table property,or on a particular meeting,does not change this fact.

At a fixed checkpoint equality,the current and anchor indices differ by d>0,with d<=L<N. The difference B_current-B_anchor equals sum_j n_j v_j for the nonnegative jump counts between those indices. Their sum is d,so some count is nonzero modulo N and the difference is uniform modulo N. Its zero probability is exactly1/N. Otherwise a=(A_anchor-A_current)/(B_current-B_anchor) and verification returns[a]Q0. Conditioning on the checkpoint schedule or an equality does not inspect coefficient representations. This uses no additional nonsingular-relation heuristic.

With2w+1 initialization scalars and at most1024 independent255-bit rejection attempts each,the sampling-failure probability under the idealized hash is below(2w+1)*2^-1024<2^-990,since N>2^254. Accepted values remain independent uniform N-residues conditioned on success; the rejection count is independent of the accepted value. Domain separation makes the partition independent of that initialization in the idealized analysis.

Subtracting bad-table probability,sampling failure and singular-relation probability gives idealized success>0.3995-2^-100-2^-990-1/N>0.399. The submitted success lower bound is the more conservative0.395. The sole remaining computational heuristic is transferring this random-function analysis to the fully computed,fixed SHA-256 construction over uniform challenges. No unconditional success theorem for fixed SHA is asserted. All work is bounded and charged even if that heuristic fails.
