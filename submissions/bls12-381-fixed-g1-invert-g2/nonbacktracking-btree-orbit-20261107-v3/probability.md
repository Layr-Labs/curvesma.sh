# Finite nonbacktracking orbit success bound

All statements are inside the disclosed standard BLAKE3 bridge H. Use the full protected prime r, q=12, W=2^76, generator degree D=12W, effective size N=r/q, delta=2^-30, T=32 and eta=2^-16. Set R=min{R:12R²>=r}, L=ceil(1031R/1000), K=L. No random-mapping, collision-independence, Poisson or heuristic mixing assumption is added. mixing.md proves the required full-state kernel bounds.

## Uniform choices and bounded retries

Before the first orbit repeat, fresh ideal hash labels give independent uniform76-bit indices and phases whose distance from uniform12 is at most12/2^256. Couple all ordinary/retry/initial phases to uniform. At most1024K+2 are used, so coupling loss is below2^-110. Given history, k is fixed, hence e=t-k is uniform. Rejecting exactly the inverse previous label selects uniformly among D-1 generators. This is the nonbacktracking Markov chain of mixing.md.

The1024-attempt local cap has failure at most K*D^-1024. Scalar rejection has failure at most2^78*(1/8)^1024. Their sum is below2^-970. Expected additional label attempts are at most K/(D-1). The global cap floor(K/2^32)>=K/2^33 therefore exhausts with probability at most2^33/(D-1)<2^-45 by Markov. Each extra attempt costs less than2^20 scalar instructions, including generic48-byte final hashes, nonce/copy, guards and counters. Global extra work is below K/4096 plus fixed exit work; one K allowance pays it. No group multiplication occurs on a rejected label. Caps bound work on every path, even outside H.

## First and second orbit-collision moments

On a good table, the extended state (raw position, previous generator) is stationary uniform. mixing.md proves pointwise joint mixing with relative error below2^-127 for any gap at least T, and raw endpoint mass at most1/(D-1) at every positive gap. It also proves reversal symmetry of the stationary raw path.

Let C count pairs i<j among0..L with j-i>=T and equal sigma-orbits. Put B=(L-T+1)(L-T+2)/2. Every long conditional orbit return has probability at most q(1+eta)/r=(1+eta)/N. The stationary raw starting value has average orbit size q-(q-1)/r, including the identity orbit. The mixing lower bound gives

 E C >= m=(1-eta)B/N,

because2^-127+(q-1)/(qr)<eta. This does not assume a uniform distribution over the orbit set.

Identical indicator pairs contribute at most(1+eta)B/N. Disjoint intervals, including touching endpoints, have joint probability at most(1+eta)^2/N² by the conditional Markov return bound. No independence is asserted.

For nested or crossing intervals, call the two outside gaps a,c. If c>=T, condition on the earlier event and apply the long-gap bound to the last return. The joint probability is at most(1+eta)^2/N². If a>=T, stationary raw-path reversal gives the same result. These distinct terms fit baseline(1+3eta)B²/N².

Only a,c<T remain. For nested intervals, inner length b>=T and a+c>0. The pointwise full-state upper bound over b permits replacing its end value and previous edge by uniform, with factor1+eta. Let U0 be the uniform initial raw value, A the a-step preceding increment, and C0 the c-step following increment after this uniform edge reset. Sum outer and inner rotation labels u,v. For u!=v, the equation

 C0=(xi^u-xi^v)U0-xi^v A

has probability1/r, giving at most(1+eta)/r² per label pair. For u=v, it requires C0=-xi^v A. At least one of a,c is positive; its endpoint distribution has maximum mass1/(D-1). Thus each equal-label case contributes at most(1+eta)/(r(D-1)). The initial raw value is independent of generator increments by the translation-invariant chain law; after the long-gap domination the reset edge is independent of A. Identity rotation overcount only enlarges this upper bound. The total is bounded by a q² baseline plus q(1+eta)/(r(D-1)). Counting both indicator orders, starts and inner lengths, with at most T² choices of outside gaps, gives nested excess

 2(1+eta)T²(L+1)²/[N(D-1)].

Shared starts/ends are a=0/c=0; identical intervals were separated.

For crossing intervals, positive a,b,c satisfy a+b>=T and b+c>=T. When a,c<T, the first event has probability at most q(1+eta)/r. Given history through its end, the final c-step orbit return to the known earlier value has probability at most q/(D-1). Hence the joint is at most q²(1+eta)/(r(D-1)). Both orders, free start/b and the short outside gaps give crossing excess

 2q(1+eta)T²(L+1)²/[N(D-1)].

All endpoint orders are disjoint, nested or crossing. Therefore

 U=(1+3eta)B²/N²+(1+eta)B/N
   +2(q+1)(1+eta)T²(L+1)²/[N(D-1)] >=E C².

For integer c>=0,1_(c>0)>=(3c-c²)/2. Consequently

 Pr[first orbit repeat by L] >=Pr[C>0]>=(3m-U)/2>3902/10000.

The exact diagnostic lower bound is0.3902220543021415. This is finite Markov moment arithmetic at the full subgroup order, not a birthday approximation.

## Hash coupling and identity

Extend an ideal nonbacktracking chain with fresh choices after repeats only to bound its first repeat and bad prefixes. The actual algorithm stops at the first orbit match, so it needs coupling only before that repeat. The ideal group trajectory is independent of tabulation tables. Distinct orbit keys have fingerprint collision probability2^-320 by the b-only injectivity lemma. Union over K+1<2^129 states loses below2^-63. The common T72[0] XOR constant cancels in differences.

Before the first repeat or fingerprint collision, all ordinary/retry inputs are fresh under their separate domains. They supply exactly the ideal choices after phase coupling. Tree ordering affects storage only, not transitions. On equal U but different canonical b, the algorithm fails rather than accepting a false match. No adaptive independence of actual fingerprints is assumed.

Uniform stationary raw marginals give identity-visit probability at most(K+1)/r<2^-100. This includes initial alpha0=0 and the sole non-free orbit. Actual identity visits abort; their loss is subtracted. Identity target z=1 succeeds separately.

## Formal coefficient coincidences

Conditional on group jump values and group trajectory, v_j remain independent uniform Fr: each chosen v has exactly one u giving M=h^u z^v. Transition, previous-generator exclusion and equality selection inspect group values, not coefficient histories. A nonzero formal coefficient vector therefore makes the extracted denominator uniform, with singular probability1/r. We bound zero vectors explicitly.

In Z[zeta12]^W, basis1,zeta,zeta²,zeta³, step Z_t=zeta^e E_j is uniform over12W units except -Z_(t-1). Initial fictitious Z0 is uniform. Consider one scalar coordinate and write b=1/(D-1). Then E[Z_t|history]=bZ_(t-1). With S_t=sum_(s=1)^t Z_s,

 M_t=sum_(s=1)^t(Z_s-bZ_(s-1))
    =(1-b)S_t-bZ0+bZ_t

is a martingale. Its increments have magnitude at most2 and conditional variance at most2/W. For0<=lambda<=1/2, conditional mgf is at most exp(2lambda²/W), since e^y<=1+y+y² for |y|<=1.

Set Bm=2^40 and lambda=Bm*W/(4K)<1/2. Each coordinate/time has tail at most2exp(-Bm²W/(8K))<2exp(-2^27). Union over4W coordinates and K times is below2^205 exp(-2^27)<2^-200. Outside that event, |S_t|<2^42 and interval coordinates are below2^43. Their algebraic norm is at most(4*2^43)^4=2^180<r. Evaluation at xi vanishing modulo r makes r divide the integer norm/resultant; a nonzero algebraic integer has nonzero integer norm. Thus modular vanishing forces an actual zero lattice coefficient. No short-relation heuristic is used.

Length1 zero return is impossible. Length2 is precisely an immediate inverse and is forbidden. For d<1024, stationary nonbacktracking sequence likelihood relative to iid uniform generators is at most(D/(D-1))^(d-1)<2. A zero sum requires no singleton index, with iid probability at most d^d W^-ceil(d/2). Hence length3/4 union is at most2K(3³+4⁴)/W²<2^-16. For5<=d<1024, use2*2^(10d-76ceil(d/2)). Odd5 starts at2^-178, even6 at2^-168, and each two-step increase multiplies by2^-56. Union is below2K(2^-178+2^-168)/(1-2^-56)<2^-32.

For d>=1024 use free-lattice Fourier orthogonality. The twelve phase averages have magnitude at most1/2 except angle measure48exp(-W/16), by Hoeffding. The uniform unit-step characteristic function phi is real because the set is symmetric. Exact nonbacktracking characteristic functions satisfy

 F0=1,F1=phi,F2=(D phi²-1)/(D-1),
 Fd=(D/(D-1))phi F_(d-1)-(1/(D-1))F_(d-2), d>=3.

This follows by counting regular-graph walks and removing immediate backtracks. If |phi|<=1/2, induction gives |Fd|<=2(2/3)^d: its inductive factor(3/4)a+(9/4)b is below1, and the first three cases hold. Return probability is at most2(2/3)^d+48exp(-W/16). Union over at most K² intervals is below2^-200;3^1024>2^1577 and K²<2^252 give a finite certificate. Martingale-prefix and long exceptions together fit2^-199.

For raw twist ell=0, a zero modular formal relation requires a zero integer interval just bounded. For ell!=0, a zero formal vector would reduce the group relation to(1-xi^ell)*alpha0=0, forcing alpha0=0, already in identity loss. Otherwise some coefficient is nonzero and the denominator is uniform. There is no checkpoint delay: the global tree stops on the first orbit repeat. Extraction is verified for the supplied P.

## Total success and exact checks

Subtract from3902/10000: fingerprint2^-63, bad table2^-98, scalar/local caps2^-970, phase2^-110, identity2^-100, short3/4 fruitless2^-16, short5..1023 fruitless2^-32, other fruitless2^-199, global retry cap2^-45 and singular1/r. The result is strictly above39/100. Using the actual finite moment bound gives more than0.3902067952802199. No independence of failure events is asserted.

All certificates use exact protected r and rational arithmetic: the displayed moment and post-loss inequalities,(41504K+2^120)^200<20000^200*2^25347,K+1<2^126,K>2^124,K²<r,2^180<r,r>7*2^252, generator/mixing bounds, and all retry/short-return bounds above. They are public mathematical checks, not production attack execution, fixed-hash security proof, hardware availability or qualification.
