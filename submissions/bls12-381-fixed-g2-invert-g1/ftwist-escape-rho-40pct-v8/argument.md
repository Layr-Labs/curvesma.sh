# Orientation-guarded Frobenius R-adding rho: fixed-G2 inversion

## Claim, input and accounting scope

The candidate score is **12881/100 = 128.81** in log2 dense-Fp12-multiplication equivalents, for success at least **0.39** on the exact supplied nonidentity Q in G2 and uniform z in GT. Protected accounting is H+R/20000, summed over every processor, failed attempt, setup, lookup, transfer, verification and output. This is a conditional exploratory cryptanalytic estimate, not a production execution or a theorem of production mixing. The named success hypotheses below are part of the claim.

Compute g=e(P0,Q), with library argument order pairing(Q,P0,True). Nondegeneracy and prime order make g a generator for every allowed Q. Recover x from z=g^x, return P=[x]P0 and verify e(P,Q)=z for that same Q. z=1 returns infinity. No auxiliary powers, replaceable Q, future-target amortization or external advice are used.

## Stationary class map and the paid orientation guard

Let p,r be the protected BLS12-381 constants, lambda=p mod r, and N=(r-1)/12. The exact certificate checks ord(lambda)=12, lambda^6=-1 and 1+lambda^4+lambda^8=0. For a nonidentity y, rep(y)=y^(p^kappa(y)) is its unique canonical Frobenius conjugate. In the tower Fp2[w]/(w^6-(1+u)), u^2=-1, write c1=a+b*u, d=(1+u)^((p-1)/6), zeta=(1+u)^((p^2-1)/6) in Fp, and Bc=a*d0+b*d1. The twelve real-coordinate keys are the interleaving of

    even=(a,a*zeta,a*zeta-a,-a,-a*zeta,a-a*zeta)
    odd =(Bc,Bc*zeta,Bc*zeta-Bc,-Bc,-Bc*zeta,Bc-Bc*zeta).

There are four fixed-constant Fp products and nine add/sub/neg operations. K is the minimum of the twelve canonical integer residues. Resolve ties by full conjugate lexicographic comparison, then index. K is not injective: equality candidates require full canonical-point confirmation. Full-point equality, not a coordinate or truncated fingerprint alone, authorizes solving the scalar equation.

Retain one random 1024-entry table for the whole run: independent uniform a_h,b_h mod r, M_h=g^a_h*z^b_h, and all twelve inverse-Frobenius twists and corresponding scalar labels. For the current noncanonical state y=g^A*z^D, keep its already computed K,k=kappa(y). Put h=K mod1024 and form the tentative R-add

    v=y*M_h^(p^-k)
    Av=A+lambda^-k*a_h; Dv=D+lambda^-k*b_h  (mod r).

Compute k_v=kappa(v). If (k_v-k) mod12=0, emit v and its labels. Otherwise emit v^2, double both labels modulo r, and compute the metadata of that squared output. A square costs one dense product, without any discounted cyclotomic-squaring assumption. Every call is charged for the extra square and metadata, whether the guard fires or not. There is no history-dependent escape, unpriced fresh offset, restart, or assumption that escapes are rare.

Indeed y^(p^k)*M_h=(y*M_h^(p^-k))^(p^k). Relative orientation k_v-k is therefore the canonical alignment of rep(y)*M_h. On canonical input U, define V=U*M_h and alpha=lambda^kappa(V). The class map emits class(V) when alpha=1 and class(V^2) otherwise. All quantities in this decision are functions of the input class and the retained table. Conjugate inputs make the same decision and give the same next class. The raw-coordinate partition is unchanged between the baseline and guard controls: only this guard predicate differs.

## Why the forced basin disappears, not merely gets detected

The unguarded same-branch inverse pair has labels -(-(A+a_h)+a_h)=A and -(-(D+b_h)+b_h)=D. Its two-cycle relation is identically zero. Bounded Floyd alone does not repair it. The inverse-only ablation tests alpha=-1 plus equal next branch; it breaks most such pairs but leaves longer structural cancellations and some residual two-cycles. The full guard is deliberately broader.

Track canonical labels and condition on the initial group element and all group jumps j_h=a_h+x*b_h. The initial target label D0 and all b_h remain independent uniform mod r. Each corrected step has label multiplier c*lambda^j, c=1 or2. An unguarded step has c=1 and j=0. Consider the segment between Floyd positions t and2t:

* If it contains s>0 guarded squares, its multiplier is 2^s*lambda^j. The certificate computes ord_r(2)/gcd(ord_r(2),12)=273103516537115575413790315146801905404638294273581446992727389062180110336 >2T. Thus 2^s*lambda^j is not1 for 1<=s<=t<=T. The denominator has a nonzero coefficient of the still-uniform D0.
* If it contains no square, every alignment in that segment is1. Its denominator is minus a sum of n_h*b_h, where n_h are nonnegative visit counts, sum n_h=t, and 0<t<=T<r. Some coefficient is strictly between0 andr. That independent b_h makes the denominator uniform.

Thus the aligned denominator vanishes with probability exactly1/r conditional on the point trajectory, not with a structural short-cycle basin probability. This is a label argument, not a mixing argument. It applies at the first Floyd meeting without a speculative one-lap or multilap repair.

Both Floyd states start from the same independently labelled random group element. One iteration advances the tortoise once and hare twice: three map calls, not one. At full confirmed class equality y_t=y_h^(p^j), recover

    x=(lambda^j*A_h-A_t)/(D_t-lambda^j*D_h) mod r.

Verify g^x=z and the final supplied-Q pairing. A zero denominator, exhausted sampling/exception budget, or the hard T-iteration limit fails; no free retries are allowed. If the trajectory has entry mu and cycle length ell, first positive Floyd meeting is ell*max(1,ceil(mu/ell))<=mu+ell.

## Faithful paired evidence, including failed ablations

The final fixed seeds are 2026102309+i for the seven configurations, 25,000 independently sampled table/start/control sets each: 175,000 sets. Each set runs the baseline, inverse-only guard and full orientation guard on the same retained jumps and start, plus one fresh uniform-map control. The resulting 525,000 design trajectories are coupled within each set, not falsely called independent samples. The final 109/100 stretch was fixed before this larger captured run.

The tower is identical to the production construction. Enumeration gives every actual coordinate minimum, exact tie representative and complete nonidentity-class branch histogram. All branches are live: histogram ranges for B=64,256,1024,1024,8,32,128 are respectively9..25,3..29,6..30,36..89,8146..8258,1947..2140,456..561. Each setting includes24 literal tower metadata checks and32 literal adding/twisted/guard checks per design. The exponent oracle accelerates exactly the same class transition; it is not a production operation or a hash substitute. The fresh map memoizes each state's destination, so repeat visits do not refresh its randomness.

All following counts have denominator25,000. Structural failures of the baseline and inverse-only guard are retained, not removed to fit a birthday model.

|N,B|baseline structural|inverse-only structural|full guard structural|baseline structural 2-cycles|full guard 2-cycles (all)|
|---|---:|---:|---:|---:|---:|
|1024,64|2879|1657|0|1414|1016|
|4098,256|2319|1716|0|711|479|
|16388,1024|2088|1774|0|362|232|
|65536,1024|2479|1882|0|669|131|
|65536,8|22783|3435|0|21894|123|
|65536,32|14317|1753|0|13354|113|
|65536,128|6027|1708|0|4578|136|

At N=65536,B=8 the fruitless basin is22783/25000 before and0 after; its structural two-cycle component is21894/25000 before and0 after. The fixed map still has123 ordinary two-cycles there, all nonstructural and useful; eliminating fruitless two-cycles does not mean eliminating every point two-cycle. Full-guard denominator zeros are5,2,0,0,0,0,0, with no structural cancellation in any setting. These accidental toy zeros are retained as failures.

|N,B|control repeat<=T0|guard useful repeat<=T|guard useful Floyd<=T|guard denominator zero|paired lower difference|
|---|---:|---:|---:|---:|---:|
|1024,64|10611|12044|16122|5|683/25000|
|4098,256|10172|11650|15322|2|91/3125|
|16388,1024|10213|11531|15251|0|71/3125|
|65536,1024|10110|11516|15248|0|82/3125|
|65536,8|11192|12505|16135|0|563/25000|
|65536,32|10303|11701|15476|0|81/3125|
|65536,128|10145|11486|15286|0|591/25000|

The CDF statistic is the paired difference I(guard useful first repeat<=T)-I(control first repeat<=T0), not an inverse-mean-time proxy. It lies in[-1,1]. Subtracting radius3/100 gives the positive lower differences shown. For25,000 independent toy sets, the one-sided Hoeffding tail is exp(-45/4) per row, and the seven-row union is below1/100. Endpoint radius3/200 has the same exponent; the fourteen two-sided endpoint tails have union28*exp(-45/4)<1/100. The certificate proves these inequalities by a rational lower series for exp(45/4). They are toy-only ideal-coin confidence statements, not a production theorem or a guarantee about a fixed deterministic PRNG seed.

Exact mu, cycle length, first-repeat and Floyd distributions, all branch visits, repeated-table accesses, raw-data SHA256 digests and DP censor counts are retained in losslessly split evidence. Actual DP bits immediately above the branch bits use widths2,4,6; unhit trajectories are counted as censored, not discarded. Sparse-production DP behavior is not inferred from these dense predicates: core success depends on Floyd, not DP transport. Toy class spaces span about10..16 bits, far below production; high-reuse B=8,32,128 stress does not prove fixed-1024 asymptotics.

## Explicit production hypotheses and success reserve

**Hmix-guard:** at T=ceil(109*T0/100), the production stationary guarded map with the actual coordinate partition and shared1024-entry table has first-repeat CDF at least that of a uniform map on Np states at T0. Here

    Np=ceil(N*1024/1023)
      =4373927678342230290889091420689347427321206073683634826518910598626352813045
    T0=ceil(sqrt(2*Np*ln(5/3)))
      =66847802276462331670412349757297191499
    T =72864104481343941520749461235453938734.

The branch correction is a disclosed modeling convention. The paired CDF and reuse stress support this particular 9% stretch at toy scope; they do not establish a quantitative transfer-error bound. Controls have N+1 class states including identity, while quotient outputs of uniform group exponents give identity different weight. Np is a conservative comparison parameter, not an exact theorem for this map.

**Hkey:** identity/exceptional-key encounters obey field-scale rates sufficient for the conservative unions6T/N and6T*1024/p. The latter covers within-element key ties and different-class key-confirmation candidates; full exact slow paths are implemented in the tower audit and metered in the cost model. The certificate shows these unions below2^-120 and2^-240. Hkey is separate from Hmix-guard and not implied by it.

A uniform-map no-repeat probability through t transitions is product(j=1..t)(1-j/Np)<=exp(-t*(t+1)/(2*Np)). An80-term rational atanh series encloses ln(5/3) with width below2^-300; both enclosure endpoints give the same T0 integer. Hence the comparison target is at least2/5. Hmix-guard and Floyd's mu+ell bound transfer that comparison to detection by T. The label proof gives1/r degeneracy. Sample every scalar by at most1024 independent256-bit rejection attempts; with fewer than4096 sampled scalars and acceptance above1/4, total exhaustion is below2^-400. The named **0.01 success reserve** from0.40 to0.39 covers these negligible conditional losses. It is not an empirical bound on extrapolation error. Correlations violating either production hypothesis can invalidate the claimed success.

## Complete worst-case work and native representation

Walk elements and twisted multipliers remain native Montgomery-encoded GT elements. The coordinate key nevertheless uses canonical residues: decode the two Fp2 c1 components with two383-instruction Fp multiplications, then compute four mixed-representation fixed-constant products, whose constants are already Montgomery encoded. A standard-residue input multiplied by an encoded constant through Montgomery reduction returns a standard-residue product. `montgomery_keys` checks this entire semantic pipeline on24 production conjugates. No coordinate decoding is silently free.

Per call, the primary rows are: two decodes766; four Fp products1532; nine six-limb add/sub/neg576; minimum/index600; two scalar add/reduce updates256; twisted loads/address256; DP dispatch128; Floyd comparison/control256; loop/meters128; state/scratch512. The guard adds orientation test128, a second complete key (766+1532+576+600), two scalar doublings128 and extra traffic256. Sum **8996**, ceiling **9000**. These are analytical fixed-limb allowances, not newly measured assembly. Scalar carries, reductions, constant loads, addressing and branch control are covered by their rows. The polynomial toy's conversion of its flattened coefficients is not an extra production operation: the production Fp2 tower stores c1 directly.

At most two dense products and9000 other instructions per call, three calls per iteration, give

    Wexec   =(147/20)*T+2^100
    Wscored =8*T+2^100
            =582912837118402132394225091380334715248
    Wstress =(381/50)*T+2^101.

The scored envelope has log2 in(128.7765,128.7766), versus submitted128.81: slack0.0334..0.0335 bits. Repricing every nonreference ceiling by6/5 gives log2 in(128.7063,128.7064), at least0.1036 bits below the claim. Dense internal instructions are not double-charged; all key Fp multiplication rows are nonreference and are stress-repriced. Unused8T engineering slack is not fictitiously executed and then repriced. Executed work itself lies in(128.6542,128.6543). There is no statistical escape-overhead term: a guarded step's complete cost is paid even when every step is guarded.

Ancillary caps: at most129*M reports at2^26 instructions, clearing at129*M*128*16 instructions,2^32 candidate checks at2^15 units, setup2^50 units and exceptions/pairing/output2^90 units; M=2^69. Their exact rational sum is below2^100 and their1.2 repricing below2^101. This covers target-dependent table construction, all random sampling, processors, memory traffic, finite candidates and final verification. Public fixed profile constants are not target advice.

## Optional fingerprints, hard storage and termination

DP reporting generations s=0..128 use zero bits10..69+s, capacity M=2^69 and a hard2^32 candidate limit. Exhaustion clears/raises density or disables reports, without changing or stopping Floyd. No expected table size is presented as a peak. Nodes occupy128 bytes: fingerprint16, labels32+32, two128-bit links32 and metadata16. Code, matrix, tables and scratch fit in the additional2^40 bytes, so peak128*M+2^40<2^78 bytes. Advice is zero.

The full canonical twelve-coordinate, six-limb point encoding has4608 bits. Apply an independent128-by4608 binary matrix,73728 bytes, to obtain a128-bit fingerprint independent of the DP predicate. Every nonzero point-encoding difference has collision probability2^-128, by conditioning on one differing column. The exhaustive4-input/3-output analogue checks all15 nonzero differences over4096 matrices. A finite false-pair expectation is at most129*M*(M-1)/2^129<2^17; Markov gives a false-pair-cap probability below2^-15, not an exact Poisson law. True duplicates are handled by the same hard candidate cap. Both caps only disable optional reporting, so neither tail is a core-success premise. Five deliberately colliding tiny fingerprints are rejected, and no wrong answer is accepted.

## Reproduction

Exact commands, captured runtime and certificate scope are in `evidence/reproduction.md`. The manifest declares no organizer experiments.
