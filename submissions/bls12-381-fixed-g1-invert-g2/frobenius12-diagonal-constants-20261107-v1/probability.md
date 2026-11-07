# Finite success argument for the12-orbit,three-component walk

Everything here is conditional on H's actual-hash-to-random-function bridge,not a theorem for fixed BLAKE3. Once inside H,no independent-random-mapping,independent-collision,Poisson orrelation heuristic is used. All bounds apply to the full protected r andevery allowed h,z with z!=1. Identity z succeeds separately. Write q=12,w=2^47,xi=p mod r,N=r/q (a rational effective size,not a replacement subgroup),rho=1/32768,T=271,eta=1/65536. Let R be the least integer with qR²>=r,L=ceil(1031R/1000),s=ceil(L/1022),K=L+s.

## Labels,sampling andstationarity

The three j indices use disjoint47-bit parts of J,so they are independent uniform w choices inside H. PHAS is independent of J. The map t=floor(qV/2^256) from uniform256-bit V has distribution at total-variation distance<=q/2^256 from uniform q:each bucket size differs from2^256/q by at most1. Couple each new orbit's t to uniform;over at mostK+1 labels,the loss is<2^-120 at the protected K. No rejection work is hidden online. Scalar rejection caps have union failure<2^-970:there are<2^50 scalars,each1024-attempt failure<2^-1024 because r>2^254. Accepted scalars remain independent uniform conditional on acceptance;counter/rejection counts do not bias their values.

Condition on the base group tables M[b,j]. Before first quotient repeat,new ideal orbit labels are independent;given history,applied k is fixed,and e=t-k is uniform q independent of the three j values. Hence raw exponents follow a translation-invariant chain with independent increments
 mu= (1/q) sum_e sigma^e(mu0*mu1*mu2),
where mu_b is uniform over the w base entries of component b. Initial X0=h^alpha0 is uniform GT,including identity,independent of increments. We analyze this uniform chain first andsubtract identity visits rather than condition it away.

For each nontrivial character chi andcomponent b,lambda_b(chi) is the average of w independent mean-zero unit complex variables. Hoeffding on real/imaginary parts andunion over3(r-1) characters gives
 Pr[some |lambda_b|>1/32] <=12r exp(-w/4096)=12r exp(-2^35)<2^-100.
On a good table,each nontrivial Fourier eigenvalue of mu is an average of products of three component eigenvalues at rotated nontrivial characters,so its modulus<=1/32³=rho. mu is sigma-invariant (andincludes negation through sigma^6);its d-step distribution p_d satisfies p_d(sigma x)=p_d(x). With epsilon_d=p_d-1/r,Fourier inversion/Parseval give sup|epsilon_d|<=rho^d andunweighted l2 norm<=rho^d. For d>=T,point probabilities lie in[(1-eta)/r,(1+eta)/r],since rho^T<=eta/r. T is deliberately conservative.

## Orbit return indicators

For a stationary uniform raw state U andindependent d-step increment S,orbit equality means U+S=xi^t U for some t. The sum over q labels overcounts only U=S=0,q times instead ofonce. Therefore
 Pr[orbit(U+S)=orbit(U)]
 =(1-(q-1)/r)*p_d(0)+(q-1)/r.
For d>=T this is between q(1-eta)/r andq(1+eta)/r;the lower inequality uses r*eta>=1. Conditional on any fixed starting state,the return probability is at most q(1+eta)/r. These statements handle the smaller identity orbit explicitly;the quotient is not falsely assumed uniform over its orbit set.

Among Y0,...,YL count C=sum I_(i,j) over i<j,j-i>=T,where I is orbit equality. Let D=(L-T+1)(L-T+2)/2. Then m=(1-eta)D/N<=E C. We upper-bound the second moment by considering all endpoint orders,with no independence claim for overlapping intervals.

Identical intervals contribute<=(1+eta)D/N. Disjoint intervals,including touching endpoints,have joint probability<=q²(1+eta)²/r² by the Markov conditional return bound. This fits the baseline(1+3eta)D²/N² because eta<1. Independence of those orbit events is not needed.

For nested intervals write independent consecutive increments A,B,C of lengths a>=0,b>=T,c>=0,a+c>0. Let the outer andinner rotation labels be u,v. After summing over the uniform starting state:
- u=v=0 gives p_b(0)*p_(a+c)(0).
- v=0,u!=0 gives p_b(0)/r.
- u=0,v!=0 gives p_(a+b+c)(0)/r.
- u=v!=0 gives p_(a+c)(0)/r,using sigma-invariance.
- u,v nonzero andunequal gives at most(1+eta)/r²:eliminating the starting state leaves xi^u A+C=kappa B for a nonzero kappa,and B is pointwise mixed because b>=T.
Thus the actual joint probability,which is at most the sum over rotation labels,is bounded by
 q²(1+eta)/r² +(q+eta)*rho^(a+c)/r
 <= q²(1+eta)/r² +q(1+eta)*rho^(a+c)/r.
The latter is exactly the inherited nested bound with N=r/q. Summing both indicator orders andall a,c uses(1-rho)^(-2)-1=65535/1073676289.

For crossing intervals i<j<k<l let the positive gap lengths be a,b,c,with a+b>=T,b+c>=T. Rotation labels u,v satisfy
 A+B=(xi^u-1)U,
 B+C=(xi^v-1)(U+A).
The q² cases are:
- u=v=0:the usual joint return bound1/r²+[rho^(a+b)+rho^(a+c)+rho^(b+c)]/r+rho^(a+b+c),by expanding the three densities into uniform+epsilon andusing l2/sup bounds.
- exactly one label zero,2(q-1) cases:at most(1+eta)/r²,from p_(a+b)(0)/r orp_(b+c)(0)/r.
- u=v!=0,q-1 cases:1/r²+rho^(a+c)/r,from(1/r)sum_x p_a(x)p_c(xi^u x).
- u,v nonzero unequal,(q-1)(q-2) cases:1/r²+rho^(a+b+c)/r. Elimination yields C=alpha A+beta B with both coefficients nonzero;scaling permutes nontrivial characters,so Fourier products bound the deviation by rho^(a+b+c).
The sum is covered by a q² baseline,the a+c excess multiplied by q,andthe original cubic excess multiplied by at most q because1+(q-1)(q-2)/r<=q. The other two quadratic excesses are no larger than the inherited bound at N=r/q. Summing positive a,c uses(rho/(1-rho))²=1/1073676289. For the constrained a+b/b+c sums one may use rho<=1/2 and2^-T<=eta/r;the cubic sum is covered by the same12eta(L+1)/N term after the q factor. This covers all distinct endpoint orders andall shared-start/shared-end cases through the nested a=0/c=0 cases.

Consequently E C²<=U where exactly the source formula applies at this effective N:
 U=(1+3eta)D²/N²+(1+eta)D/N
   +((131072+131070eta)/1073676289)*(L+1)²/N
   +8T eta*(L+1)²/N²+12eta*(L+1)/N.
For every nonnegative integer c,1_(c>0)>=(3c-c²)/2. Thus Pr[C>0]>=(3m-U)/2>39015/100000 at the protected r,q,R,L. This is a finite second-moment/Bonferroni calculation,not a Poisson oruniform-mapping approximation.

## Adaptive hash coupling andcheckpoint delays

Define a surrogate independent of fingerprint tables:assign an ideal independent (j0,j1,j2,t) to each newly visited orbit andcache it on revisits,then use the same equivariant transition. Before first orbit repeat it is the iid-increment chain above. Its entire K-step quotient tail/cycle is independent of the tabulation entries. Fixed unequal canonical orbit keys have full-fingerprint collision probability2^-320. Union over at mostK+1<2^129 states loses<2^-63,below the conservative2^-13 used in the claim. Couple actual JUMP/PHAS calls to the surrogate whenever a new U appears;revisited orbits have the same U andsame labels. Without a distinct-orbit fingerprint collision andthe separately charged phase-bias event,the entire K-step paths agree. This avoids conditioning on an already table-dependent path.

The deterministic checkpoint lemma in proof.md detects a first quotient repeat byL no later thanK. Let its raw twist be ell,so X_(mu+lambda)=sigma^ell X_mu. Equivariance makes later paired jumps differ by xi^ell. The formal A/B coefficient relation V_(b+lambda)-xi^ell V_b is unchanged as b advances along the cycle. Thus the first checkpoint-detected relation has the same formal vector as the first quotient repeat,even if its first saved cycle point occurs after several raw traversals. No new singularity is introduced by checkpoint delay.

Prefix cap:conditional on the full table-independent surrogate andideal first exact-orbit/checkpoint stop,there are<=1023K unequal query/checkpoint pairs. Each64-bit prefix collision has probability2^-64;sorting/scanning chooses a subset. Hence expected full confirmations<=1023K/2^64+1. With B_cap=floor(K/256)>=K/512 andK>2^64,Markov bounds exhaustion<2^-43. The full comparison cap pays all its work,including the exhausted check;its loss is subtracted separately from full-fingerprint collisions,not conditioned on their absence.

## Fruitless relations:free cyclotomic walk

Condition on group table values,orbit labels andinitial alpha0. The base v[b,j] remain independent uniform Fr because for every chosen v there is exactly one u giving its group value M=h^u z^v. The trajectory,orientations andequality selection use only group states,not coefficient histories. Therefore the extracted denominator is uniform Fr if its formal coefficient vector is nonzero modulo r. We now bound the exception,rather than assume all quotient collisions are useful.

Write zeta for a complex primitive12th root anduse the basis1,zeta,zeta²,zeta³ of Z[zeta],withPhi12=zeta^4-zeta²+1. Each zeta^e has coefficients in{-1,0,1},andopposite phases negate. In the iid pre-repeat process,a formal step is
 zeta^e*(E_(0,j0)+E_(1,j1)+E_(2,j2))
in Z[zeta]^(3w). The common phase is uniform12 independent of history;the three indices are independent uniform w. This follows from t-k uniformity,not a new walk assumption. Extend these independent choices to K steps solely for bounding bad prefixes.

First exclude modular coincidences not visible in the integer lattice. Each coordinate increment is0 or±1,symmetric,nonzero with probability<=1/w. For0<=lambda<=1,its mgf<=exp(lambda²/w),using cosh(lambda)-1<=lambda². Set B0=2^55,lambda=B0*w/(2K)<1. Chernoff for each coordinate/time gives
 Pr[|partial sum|>=B0]<=2 exp(-B0²*w/(4K))<2 exp(-2^29).
Union over K times and12w coordinates is<2^178 exp(-2^29)<2^-200. On its complement every interval coefficient has four coordinates of magnitude<2B0=2^56. Its algebraic norm is bounded by(sum|c_i|)^4<=(4*2^56)^4=2^232<r. If its evaluation at xi in Fr is zero,r divides this integer norm (the norm/resultant modulo r has a zero factor). Nonzero algebraic integers have nonzero integer norm. Hence a coefficient vanishes modulo r only when all four integer coordinates vanish. No short-relation assumption about xi is made.

For raw twist ell=0,a zero formal relation therefore requires a zero interval of this integer walk. Length1 cannot return. Length2 returns exactly when all three indices repeat andphases are opposite,probability1/(12w³). Union over starts gives K/(12w³)<2^-18.

For3<=d<1024,zero return implies no singleton index in each component table. In one component there are at most d^d partitions into repeated-index blocks,andat most floor(d/2) distinct indices,so probability<=d^d*w^(-ceil(d/2)). The three index arrays are independent,even though phase is shared. Thus return probability<=d^(3d)*w^(-3ceil(d/2))<=2^(30d-141ceil(d/2)). Summing the odd/even geometric bounds over this range is at most(2^-192+2^-162)/(1-2^-81)<2^-161. Union over K starts is<2^-35;we retain2^-32.

For d>=1024,use Fourier orthogonality of the free integer lattice,not finite-field mixing. For uniform character angles theta,let a_(b,e) be the average over w indices of exp(i<zeta^e E_(b,j),theta>). For fixed b,e these are averages of independent mean-zero unit complex variables. Hoeffding and36 cases give Pr_theta[some |a|>1/2]<=144 exp(-w/16). On its complement the step characteristic function phi=(1/12)sum_e product_b a_(b,e) has modulus<=1/8. Therefore return probability<=E|phi|^d<=8^-d+144 exp(-w/16). Union over all<=K² intervals is<2^-200 at w=2^47,K<2^126. Combining modular-coincidence failure,long returns andshort returns,the fruitless raw-repeat loss is<2^-18+2^-32+2^-199.

For ell!=0,a zero formal coefficient vector in Fr would make the group relation reduce to(1-xi^ell)*alpha0=0. As xi has order12 andr is prime,this forces alpha0=0. Initial identity is included in the identity-visit loss below. Thus nonzero twists need no extra free-lattice estimate. On the complement of all these exceptions,a formal coefficient is nonzero,andthe denominator is uniform Fr with singular probability1/r. Otherwise extraction/output follows algebraically from X=h^A z^B andis verified.

## Remaining losses andexact public checks

Before first orbit repeat,the iid chain is uniform at each raw state. Pr[any identity visit]<= (K+1)/r<2^-100. After a quotient repeat,no new orbit occurs before checkpoint detection,so no additional identity can appear if the pre-repeat part avoided it. This also handles the only orbit with nontrivial stabilizer andthe otherwise undefined equivariant map there.

All events are union-bounded without asserting mutual independence. Success is greater than
 39015/100000
 -2^-13 fingerprint -2^-43 cap -2^-100 group-table
 -2^-970 sampling -2^-120 phase bias -2^-100 identity
 -2^-18 fruitless2 -2^-32 short fruitless -2^-199 remaining fruitless
 -1/r singular denominator
 >39/100.
The proof uses conservative2^-13 despite the stronger full-fingerprint bound. Identity target z=1 succeeds deterministically. The reduction holds for every allowed supplied P;therefore the success statement has the required independent uniform P,z scope.

Public finite certificates use exact rational N=r/12 andinteger R=min{R:12R²>=r},L=(1031R+999)//1000,s=(L+1021)//1022,K=L+s:
 (3m-U)/2>39015/100000;
 the displayed post-loss bound>39/100;
 (74336K+2^120)^1000<20000^1000*2^127577;
 K<2^126,K>2^124,K²<r,1023s>K;
 K/(12w³)<2^-18;
 (1023K/2^64+1)/floor(K/256)<2^-43;
 12(K+1)/2^256<2^-120;
 2^232<r.
These are independently checkable public arithmetic,not production experiments orproof of H. No candidate file has executed on the host orwith credentials.
