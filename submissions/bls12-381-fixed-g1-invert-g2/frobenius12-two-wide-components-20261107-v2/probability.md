# Finite success bound for two wide-component Frobenius12 walk

Everything is conditional on H's actual-standard-hash random-function bridge,not a fixed-hash theorem. Inside H there is no random-mapping,collision-independence,Poisson orrelation heuristic. Use full protected r,q=12,w=2^72,xi=p mod r,N=r/q,rho=1/32768,T=271,eta=1/65536,R=min{R:12R²>=r},L=ceil(1031R/1000),s=ceil(L/1022),K=L+s.

## Labels, sampling andstationary increments

The two72-bit j indices use disjoint128-bit digest halves,with only their low72 bits selected. They are independent uniform w choices under H. PHAS is independent. t=floor(12V/2^256) differs from uniform12 in total variation<=12/2^256;over<=K+1 labels the coupling loss is<2^-120. Scalar generation uses128-bit counter_16 and64-byte SCAL messages. At most4w+1<2^75 scalars and1024 attempts each fit<2^85 counter values. Protected r>7*2^252 makes rejection<1/8,so all scalar failures total<2^75*2^-3072<2^-970. Accepted scalars are uniform independent oftheir rejection counts.

Given the base group tables,before first orbit repeat fresh labels give independent j0,j1 anduniform phase e=t-k. Thus raw exponents follow independent increments
 mu=(1/12)sum_e sigma^e(mu0*mu1).
Initial h^alpha0 is uniform raw GT,including identity,independent of increments;identity loss is subtracted rather than conditioned away.

For each component andnontrivial character,lambda_b is an average ofw independent mean-zero unit complex variables. Hoeffding andunion over2(r-1) characters gives
 Pr[some|lambda_b|>1/256] <=8r exp(-w/262144)=8r exp(-2^54)<2^-100.
On a good table,every nontrivial eigenvalue ofmu has magnitude<=1/256²=1/65536<=rho. Symmetrization is sigma-invariant. Fourier inversion/Parseval give sup|epsilon_d|<=rho^d andunweightedl2<=rho^d. For d>=T,point probabilities lie in[(1-eta)/r,(1+eta)/r]. All following orbit second-moment steps are unchanged anddepend only onthese properties,not the number ofcomponents.

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

Define a surrogate independent of fingerprint tables:assign an ideal independent (j0,j1,t) to each newly visited orbit andcache it on revisits,then use the same equivariant transition. Before first orbit repeat it is the iid-increment chain above. Its entire K-step quotient tail/cycle is independent of the tabulation entries. Fixed unequal canonical orbit keys have full-fingerprint collision probability2^-320. Union over at mostK+1<2^129 states loses<2^-63,below the conservative2^-13 used in the claim. Couple actual JUMP/PHAS calls to the surrogate whenever a new U appears;revisited orbits have the same U andsame labels. Without a distinct-orbit fingerprint collision andthe separately charged phase-bias event,the entire K-step paths agree. This avoids conditioning on an already table-dependent path.

The deterministic checkpoint lemma in proof.md detects a first quotient repeat byL no later thanK. Let its raw twist be ell,so X_(mu+lambda)=sigma^ell X_mu. Equivariance makes later paired jumps differ by xi^ell. The formal A/B coefficient relation V_(b+lambda)-xi^ell V_b is unchanged as b advances along the cycle. Thus the first checkpoint-detected relation has the same formal vector as the first quotient repeat,even if its first saved cycle point occurs after several raw traversals. No new singularity is introduced by checkpoint delay.

Prefix cap:conditional on the full table-independent surrogate andideal first exact-orbit/checkpoint stop,there are<=1023K unequal query/checkpoint pairs. Each64-bit prefix collision has probability2^-64;sorting/scanning chooses a subset. Hence expected full confirmations<=1023K/2^64+1. With B_cap=floor(K/256)>=K/512 andK>2^64,Markov bounds exhaustion<2^-43. The full comparison cap pays all its work,including the exhausted check;its loss is subtracted separately from full-fingerprint collisions,not conditioned on their absence.

## Fruitless relations in the free cyclotomic walk

Given group tables,orbit labels andalpha0,the base v[b,j] remain independent uniform Fr:for each v there is exactly oneu giving M=h^u z^v. Trajectory/equality selection uses group states only. Therefore a nonzero formal coefficient vector modulo r makes the extracted denominator uniform. We bound zero formal relations explicitly.

In Z[zeta12] with basis1,zeta,zeta²,zeta³ andPhi12=zeta^4-zeta²+1,a pre-repeat step is zeta^e*(E_(0,j0)+E_(1,j1)). Phase is uniform12 andindices are independent uniformw,as proved above. Extend iid choices toK steps only to bound bad pre-repeat prefixes. Each coordinate increment is symmetric0/±1,nonzero with probability<=1/w. For0<=lambda<=1,mgf<=exp(lambda²/w). Set B0=2^44,lambda=B0*w/(2K)<1 using K>2^124. Chernoff bounds each coordinate/time by2exp(-B0²*w/(4K))<2exp(-2^32). Union overK times and8w coordinates is<2^202 exp(-2^32)<2^-200.

Outside that event,an interval coefficient has four coordinates<2B0=2^45. Its integer algebraic norm is at most(4*2^45)^4=2^188<r. If evaluation atxi vanishes modulo r,the integer norm/resultant is divisible byr. A nonzero algebraic integer has a nonzero integer norm,so this bound forces the coefficient itself tozero in the lattice. No short algebraic-relation heuristic is used.

For raw twistell=0,a zero formal relation therefore requires a zero integer interval. Length1 cannot return. Length2 requires both indices torepeat andopposite phases,probability1/(12w²),with unionK/(12w²)<2^-18 at full constants.

For3<=d<1024,each component's index sequence must have no singleton. Its probability is<=d^d*w^(-ceil(d/2)). The two sequences are independent,so joint return probability<=d^(2d)*w^(-2ceil(d/2))<=2^(20d-144ceil(d/2)). Oddd starts at3with2^-228,evend at4with2^-208;each increase byd2 multiplies byatmost2^-104. Their sum<(2^-228+2^-208)/(1-2^-104)<2^-207. Union overK<2^126 starts is<2^-81;retain the much looser2^-32.

For d>=1024,use Fourier orthogonality on the free integer lattice. With uniform character angles,each of24 means a_(b,e) is an average ofw independent mean-zero unit complex variables. Hoeffding gives exceptional angle measure<=96exp(-w/16). Otherwise |phi|=|(1/12)sum_e a_(0,e)*a_(1,e)|<=1/4. Thus return probability<=4^-d+96exp(-w/16). Union over<=K² intervals is<2^-200 atw=2^72,K<2^126. Together with modular-prefix failure this contributes<2^-199. These are finite inequalities,not return heuristics.

For raw twistell!=0,a zero formal coefficient vector modulo r would reduce the group relation to(1-xi^ell)*alpha0=0. Sincexi hasorder12,this forcesalpha0=0,already included in the identity-visit loss. Otherwise some formal coefficient is nonzero andthe denominator is uniform,singular with probability1/r. Checkpoint delay preserves the first repeat's formal vector byequivariance,as proved above,so it adds no later singularity.

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
 (52032K+2^120)^500<20000^500*2^63531;
 K<2^126,K>2^124,K²<r,1023s>K;
 K/(12w²)<2^-18;
 (1023K/2^64+1)/floor(K/256)<2^-43;
 12(K+1)/2^256<2^-120;
 2^188<r.
These are independently checkable public arithmetic,not production experiments orproof of H. No candidate file has executed on the host orwith credentials.
