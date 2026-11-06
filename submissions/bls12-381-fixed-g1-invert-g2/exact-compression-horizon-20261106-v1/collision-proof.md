# Finite collision bound in the ideal-hash model

Let N be the protected prime subgroup order, w=1099511627776, and M_1,...,M_w independent uniform elements of GT. The implementation generates these from uniformly idealized initialization scalars; no test for a good table or discrete-log oracle is performed. This note proves a finite collision bound in that idealization. It does not assert that fixed SHA-512 is a mathematically random function.

## Good-table probability

For every nontrivial character chi of the cyclic group, define lambda_chi=(1/w)*sum_j chi(M_j). Its real and imaginary summands lie in[-1,1],have mean zero,and are independent over j. If |lambda_chi|>1/32768,one component has magnitude>1/(32768sqrt(2)). Hoeffding's inequality on each component and a union bound give probability<=4exp(-w/(4*32768^2)). Union over N-1 characters gives bad-table probability<=4Nexp(-256)<2^-100 for the protected N<2^256. For an exact check, e>8/3 implies this is below4N*(3/8)^256; the claimed dyadic bound is the integer inequality4N*3^256*2^100<8^256. This is an analytic public-integer inequality; no experimental result is claimed.

Condition henceforth on a good table, so every nontrivial eigenvalue of its uniform-addition Cayley transition matrix has modulus<=rho=1/32768. Duplication among the sampled M_j is allowed and does not invalidate Hoeffding or the transition weights. The attack pays for all samples and accepts bad tables too, counting their possible failures.

## Markov surrogate and coupling

Consider independent uniform indices J_t in[1,w] and the Markov walk Y_(t+1)=Y_t*M_(J_t). A random-function partition assigns independent uniform indices to previously unvisited elements. Therefore its deterministic iterated walk and this independent-index walk can be coupled exactly until the first repeated group element. They need not agree after that repeat. Only the distribution of the first repeat is transferred below. This revision first compresses encoded states with an independently keyed polynomial family. Coupling also stops on a collision of distinct states under that compression. universal-proof.md bounds the probability of that extra stopping event by less than2^-13 on the independent-index surrogate. Subtract this event after the first-repeat bound. Equally spaced checkpoint detection runs on the actual deterministic map,not on a fresh-index walk.

Write p_d(x) for the d-step increment distribution and epsilon_d(x)=p_d(x)-1/N. Cyclic convolution is diagonalized by the characters,including complex eigenvalues for this possibly nonreversible chain. Fourier inversion and Parseval imply

sup_x |epsilon_d(x)|<=rho^d, and (sum_x |epsilon_d(x)|^2)^(1/2)<=rho^d.

These also hold at d=0. Set R=ceil(sqrt(N)), L=ceil(25767R/25000), T=ceil(log2(N))+16=271 and eta=2^-16. Since rho^T<=eta/N, for d>=T the return probability p_d(1) lies between(1-eta)/N and(1+eta)/N.

## Long-pair second moment

Among states Y_0,...,Y_L,count collisions C=sum I_(i,j) over i<j and j-i>=T. The number of such pairs is D=(L-T+1)(L-T+2)/2,so E C>=(1-eta)D/N. All equality events depend only on independent increments over their intervals; the starting group value does not affect the argument.

We bound ordered pairs of indicators in C^2. Identical intervals contribute at most(1+eta)D/N. For distinct intervals,include a baseline(1+3eta)/N^2 per ordered pair,whose sum is at most(1+3eta)D^2/N^2. The following bounds account for every excess above that baseline.

Disjoint increment intervals,including intervals touching at one endpoint,are independent and have joint probability<=((1+eta)/N)^2,already covered. Intervals sharing their first or last endpoint belong to the nested case below; identical intervals have already been separated.

For nested intervals with consecutive segment lengths a>=0,b>=T,c>=0 and a+c>0,the inner return and the outer residual return depend on disjoint increments. Their joint probability is p_b(1)*p_(a+c)(1),at most(1+eta)/N^2+(1+eta)*rho^(a+c)/N. There are at most(L+1)^2 choices for starting position and b. Sum over a,c with a+c>0 gives sum rho^(a+c)=1/(1-rho)^2-1=65535/1073676289. Both orders of the two indicators give total excess<=(131070/1073676289)(1+eta)(L+1)^2/N. This includes a=0 or c=0,so shared-endpoint nested intervals are counted.

For crossing intervals,the four ordered endpoints have positive successive gaps a,b,c with a+b>=T,b+c>=T. Independence of these three increment blocks gives joint probability sum_x p_a(x)*p_b(x^-1)*p_c(x). Expand each p as uniform+epsilon. Linear terms vanish. Cauchy-Schwarz on quadratic terms and sup/l2 bounds on the cubic term give

joint <=1/N^2+(rho^(a+b)+rho^(a+c)+rho^(b+c))/N+rho^(a+b+c).

Sum the rho^(a+c) term over positive a,c:the sum is(rho/(1-rho))^2=1/1073676289. Free starting position and b,plus the two indicator orders,give excess<=(2/1073676289)(L+1)^2/N. For the other two quadratic terms,sum_(a,b>=1,a+b>=T)rho^(a+b)<=2T*2^-T. The free remaining segment and start,with both orientations and both terms,give excess<=8T*eta*(L+1)^2/N^2.

For the cubic term,its constrained sum over a,b,c is below6*2^-T (dominate termwise by rho=1/2). Indeed,b>=T-1 contributes at most4*2^-T; for1<=b<=T-2,the conditions require a,c>=T-b and the geometric sum contributes below2*2^-T. After both indicator orders and the starting position,the contribution is below12eta*(L+1)/N. No interval pattern is omitted:ordered distinct intervals are disjoint,nested,or crossing.

Combining terms proves

E C^2 <= (1+3eta)D^2/N^2 +(1+eta)D/N +((131072+131070eta)/1073676289)(L+1)^2/N +8T*eta*(L+1)^2/N^2 +12eta*(L+1)/N.

This is a finite bound,not an assumption of independent collision indicators. In particular the short-return correlations of an adding walk are explicitly included.

## Integer-valued Bonferroni probability bound

For every nonnegative integer c, 1_{c>0} >= (3c-c^2)/2. Indeed the polynomial is 0 at c=0, 1 at c=1 and c=2, 0 at c=3, and nonpositive for c>=3. Hence for the collision count C,

P(C>0) >= (3 E C - E C^2)/2 >= (3m-U)/2,

where m=(1-eta)D/N is the lower first-moment bound and U is the upper second-moment bound proved above. No independence of collision indicators is assumed, and no third or higher moment is required. This is the two-term inclusion-exclusion (Bonferroni) bound written in terms of raw moments. It supplements the earlier Cauchy-Schwarz estimate; it does not alter any nested or crossing interval accounting.

Use the protected full-production N=r, R=ceil(sqrt(N)), L=ceil(25767R/25000), T=271 and eta=1/65536. Define D=(L-T+1)(L-T+2)/2 and the same U as above:

U=(1+3eta)D^2/N^2+(1+eta)D/N+((131072+131070eta)/1073676289)(L+1)^2/N+8T eta (L+1)^2/N^2+12eta(L+1)/N.

An exact rational calculation gives (3m-U)/2 > 39015/100000. Subtracting compression loss 1/8192, bad-table loss 2^-100, bounded-sampler loss 2^-980 and singular denominator loss 1/N leaves more than 39/100. All these are exact finite arithmetic comparisons specified in horizon-audit.md, not reported production experiments. The checkpoint schedule remains K=L+ceil(L/62). A first repeat by L is still detected by K deterministically. The old longer-horizon Cauchy-Schwarz numerical examples are superseded; only this shorter-horizon Bonferroni bound supplies the current success claim.

## Useful relation and total success

Let z=h^a. Conditional on the group values M_j,each v_j is still uniform and independent:for every v_j there is exactly one u_j giving that M_j=h^u_j*z^v_j. The trajectory,partition and initial group value do not inspect these coefficient representations. Conditioning on the good-table property,or on a particular meeting,does not change this fact.

At a saved checkpoint equality,the current and anchor indices differ by d>0,with d<=K<N. The difference B_current-B_anchor equals sum_j n_j v_j for the nonnegative jump counts between those indices. Their sum is d,so some count is nonzero modulo N and the difference is uniform modulo N. Its zero probability is exactly1/N. Otherwise a=(A_anchor-A_current)/(B_current-B_anchor) and verification returns[a]Q0. Conditioning on the checkpoint schedule or an equality does not inspect coefficient representations. This uses no additional nonsingular-relation heuristic.

With2w+1 initialization scalars and at most1024 independent255-bit rejection attempts each,the sampling-failure probability under the idealized hash is below(2w+1)*2^-1024<2^-980,since N>2^254. Accepted values remain independent uniform N-residues conditioned on success; the rejection count is independent of the accepted value. Domain separation makes the partition independent of that initialization in the idealized analysis.

Subtracting all four losses from (3m-U)/2 gives idealized success strictly greater than 0.39. The submission declares the conservative 0.39 conditional on the sole disclosed SHA-512 idealization. Input identity targets succeed deterministically. No unconditional fixed-SHA theorem or measured production success is asserted. All cold work is bounded independently of whether the idealization holds.

## General spectral coefficient check

For rho=1/d the sum of nested excesses is 2(1+eta)((1-rho)^(-2)-1)(L+1)^2/N, and the crossing rho^(a+c) excess is 2(rho/(1-rho))^2(L+1)^2/N. Their total coefficient is (4d+(4d-2)eta)/(d-1)^2. Here d=32768, giving (131072+131070eta)/1073676289. All other tail bounds above used only rho<=1/2 and rho^T<=eta/N, so remain valid. With w=2^40, w/(4d^2)=256, and the same 4N exp(-256)<2^-100 table-failure bound applies. No condition is tested by the attack; all sampled tables and failures are charged.
