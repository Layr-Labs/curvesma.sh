# Recover-P 117.619455262253: f-only collection, certified counts and audited norm transport

Track `bls12-381-fixed-g2-invert-g1`. The first manifest claim gives the exact rational score. The complete argument is this file together with DESCENT.md, KERNELS.md, STORAGE.md and field-map.md. The bound is analytic and conditional on the disclosed heuristics, including GRH, H5n and H5p.

## Fixed supplied-Q reduction

Validate the supplied nonidentity Q in G2 and the GT target z. For z=1 return infinity. Otherwise compute g=e(P0,Q), using the pinned `pairing(Q,P0)`. Map both g,z to the STNFS tower using the cold isomorphism construction in field-map.md. One relation/LA phase supports TWO individual logarithms, of g and z. Recover a=log_g(z) modulo r and return [a]P0 only if e([a]P0,Q)=z. Nondegeneracy makes g a generator. Uniform Q gives uniform nonidentity g independently of uniform z, exactly the two-log input law. No replacement Q or free advice is used.

The wrapper has a fixed2^80-instruction cap, including all random field-factorization attempts, pairings, mapping, scalar multiplication, output and checks. At C=20000 it fits2^66 units; failure allowance2^-48 (field-map.md). Every failed check or exhausted cap is a failure.

## Inherited STNFS collection and arithmetic

Use h=y^6+y^3+y^2-y-1, f=x^12-2yx^10+2y^3x^6+y^5x^2+1+y-y^2-y^3, g=x^2+U y, U=15132376222941642752; x->-x is the odd sigma quotient. config.json fixes the exact rational region: A²=2^21, pooltheta3/8, regiontheta2719/8000, lam17/20, kap1. Exact canonical unit/HNF and primitive tests are implemented in ideal_canonical.py/runs16.py. lattice_count_bounds.py supplies rational yield/work geometry factors0.970276956/1.028809574.

Only f is sieved; the exact SRT g batch decides g smoothness. model16.py, KERNELS.md and STORAGE.md specify the line/record/run walks, certified thresholds and warm starts, contamination fallback, routing/recovery, division, bandstrip/split, cofactor-window ECM and g splitting. Every setup, tree, failed trial, transfer and fallback is charged there; DESCENT.md gives both logarithms. S/G eachhave656000draws; their recorded work sums/squares and yield sums in pool_d.json are empirical inputs. The supplied anchors are recorded inputs. The fast certificate does not regenerate that historical campaign.

Credit gbotrel's transferred f-only scan/bound deletion: g-array operations are removed with unchanged f predicate. The shipped raw kernel log yields exact slope260281025/189262038<=1.376, leaf excess148329356571/3942959125<=37.62, bound cost3416387/2959<=1154.6 and residual2271913/23788165<=.125<.5. The certificate decodes/rechecks every point. S/G work SE=.004213089760336546 and1.00843>=1+2SE support the explicit H2 capped mean<=2^116.92454; the actual hard meter is2^116.925. G yield after.009-bit bandloss is2^63.62660550353605; its2-SE lower2^63.576124844219954 supports descent input2^63.575. The fixed inherited region remains feasible with the corrected U; the historical S minimum-grid rule is not reasserted. H1/H2/H4/H-FFT assumptions and their measured scope remain explicit in the associated files.

## Direct prime-ideal counts (reproduced proof)

The truncated li expansion is a heuristic density only. Dimensions and repeat bounds use the direct derivation reproduced at the end of this file, conditional on GRH for the two fixed Dedekind zeta functions. It derives, for all integer69<=b<=225,

|pi_K(2^b)-integral(3,2^b,dt/logt)| < E_b=2^(ceil(3b/4)+3).

In particular E_69=2^55. No quantitative Grenie-Molteni estimate or external zero-count/Laurent-constant bound is used. The standard analytic continuation, functional equation and Hadamard product are stated with their hypotheses; the zero-weight estimate, every explicit-formula residue, unsmoothing and partial-summation constants are derived below. Numerical checks support the reproduced argument rather than purporting to prove an assumed theorem.

The absolute resultant polynomials have degrees72/12 and are irreducible over Q (Wolfram code/raw in stats16.b85). Since |y|<=2, at |x|=3 the f leading magnitude531441 exceeds the other upper magnitudes248163, so all f roots satisfy|x|<3. Thus logDelta_f<=72*71*log6<10000. For g, |x|<2^33 because U<2^64, giving logDelta_g<=4488log2<10000. Field discriminants divide their monic polynomial discriminants.

For F(x)=x/(logx-1.1), F'(x)>=1/logx when logx>=12.1. Splitting at2^18 gives integral(3,B)<=F(B)+2^18. A sigma-fixed ideal has zero odd log; other orbits have size two. The direct bound counts ALL absolute prime ideals and upper-bounds the relative-degree-one subset. The field-specific endpoint refinement proved below bounds the average of f/g count errors by17523154306377711616/1175<15000000000000000. This is sufficient for the half-sum; individual repeat counts still use E_b. With redundant boundary84,

U_actual<=ceil[B/(69log2-11/10)+2^18+17523154306377711616/1175+84]=12647734742496954934.

counts16.py uses100positive atanh terms as a rational lower log2 bound. It checks the direct-proof constants, root/discriminant bounds and exact ceiling. Ncap=.03330Ucap<2^58.55<6*2^56, so six owners suffice. The complete main-phase appendix below defines actual-U production moments, the protected filter, rank/unit prerequisites and every cap/projection tax. The larger cap gives LA<=2^115.822475078 and generator<=2^110.610733127.

Independently on EACH side in every actual bin[2^a,2^b], subtract E_a+E_b, at most72*2^ceil(b/2) higher-degree ideals, and2^20 ideals above bad polynomial primes. Simple nonzero roots form free sigma pairs. Each orbit pool is at least

[(2^b-2^a)/(b log2)-E_a-E_b-72*2^ceil(b/2)-2^20]/2 >2^68.

All31exact rational checks pass; minimum floor357084587023637308834. Under the explicit H5 arrival-weight allowance2^6, Nmax=2^27.3 gives repeat tax<=2^(2*27.3-1+6-68)=2^-8.4. No abort is added by discarding the exceptional primes in this lower count.

## Runtime norm-CDF audit and empirical transport (new)

The historical norm bins are proposed thresholds, not universal envelopes. Before using a node basis, audit its finite integer ball with fresh, exact uniform samples. For the24node thresholds(a_i,b_i), replace b_i by min_{j>=i}b_j. These nested rectangles E_i have monotone thresholds. Require true mass P(E_i)>=eta*i/24 with eta=19/20; boot uses its10scalar thresholds, checked by exact whole-ball enumeration for each fresh basis. A layer-cake sum proves that these CDF bounds imply at least eta times the original sub-bin semismooth lower sum for every coordinatewise decreasing smoothness function. Since the old charge was the minimum of that sum and a capped sextile estimate, multiplying the old charge by eta remains conservative. No sampled maximum is treated as a universal bound.

At ordinary nodes use2^46samples, epsilon=2^-15, accepting only if every empiricalCDF>=eta*i/K+epsilon. Hoeffding gives total false acceptance<=48exp(-2^17), and false rejection on a basis whose trueCDF>=eta*i/K+2epsilon is equally negligible. A uniform integer-box proposal in[-R,R]^12, rejecting points outside the ball and zero and folding signs, is exactly uniform on the half ball. Its acceptance probability is>=2^-16 for all shipped radii (checked from integer CELLS). Cap proposals at2^64 per audit; the binomial lower-tail abort and capped uniform-integer generation are charged. Exact norm-threshold comparisons use precomputed integer100th-root thresholds; no unpriced arbitrary-width unit operation is used. Total audit work, including rejected bases, setup and randomness, fits a separate2^105-unit row (detailed in DESCENT.md).

A bad node basis is rejected and fails. H5n is the explicit quantitative transport heuristic: under the fresh special-q distribution actually reached in each production bin, the mass of bases violating the strongerCDF condition is<=delta=1/200. This is NOT deduced as a theorem from pooled samples. A fresh held-out replication uses2000production-sized special-q lattices covering62side/size strata (each sampled in the top1bit of its5/6bit bin) and32fresh uniform-GT boot inputs; recorded minimum-constraint witnesses and exact statistical summaries are in stats16.b85. All2000node and32boot tests pass at eta=.95,2048norms per lattice, conservative eps=.041/.036. The certificate recomputes these empirical summaries and the pooled binomial/Chernoff evidence, not every lattice or the historical full campaign.

Scope: witness primes/root choices and top1bit samples differ from the production arrival law. H5n extrapolates across root multiplicity, weighting, within-bin shapes and exact-norm to certified-bound classification. Offline Gaussian/rounding sampling and floating logs are approximations, not exact certificates. The positive margins and coverage are evidence for this testable assumption, not a proof or an algorithmic95%success probability. The complete generator sources and fixed seeds are embedded in stats16.b85 with the minimum-constraint witnesses; those witnesses are empirical inputs, not a replay of every recorded count. The runtime audit removes reliance on untested CDF envelopes for accepted bases; H5n supplies its rejection-rate tax. A62-lattice/992point calibration gives1984positive gaps, zeroCDF/cap disagreements and maxgap5.161e-10bits, at that scope only.

The conditional negative-binomial/independent-subtree model H5 remains, supported at its measured scope by the existing toy fits and exact unit-orbit checks (DESCENT.md). Its0.02floor is additional padding, not the norm-tail allowance. The corrected recurrence is phi=delta+(1-delta)*max(phi_NB,0.02). Every bad-lattice failure is therefore ADDED. Work uses the corresponding retry factors; charging a full node even when its audit fails is conservative. Boot resamples fresh exponents after failedCDF audits, capped at128 attempts. H5n also assumes boot bad mass<=1/2;32positive input tests support this. Fewer than2^40good balls are attempted under the charged meter, so the union of128-attempt audit exhaustion events is<=2^-88. All attempts are charged. The model uses safe LOG2U6=382.286, above the exact382.28563411957; LOG2P6 remains the existing safe upper2285.982.

## Meter, success and complete ledger

The recurrence, node count and boot are recomputed by audit_refinement.py under the20-second hosted experiment limit. Ordered zero-only convolution trimming preserves terms. Moving upper omitted-mass tails above512bits to infinite mass is conservative: all denominator cutoffs are below512bits, all increments are nonnegative, and the loss numerator is bounded by1. Lower omitted laws remain uncapped. No recorded descent-output cache substitutes for the recurrence.

Expected two-log work<=2^113.329, boot phi<=.001755, expected node log2<=16.682, node-abort charge2^-10.61, repeat charge2^-8.4, new audit-error charge2^-80 and prefix-abort charge2^-13 are in refinement.json. The meter for BOTH logarithms is2^114.085 units. Markov charges E[W]/meter; node overflow is charged from the recomputed expected count; repeats use2^-8.4; all audit, boot and wrapper failures are charged. All failure charges are combined by the explicit UNION bound in the main-phase appendix, including production matrix bad mass.01 and the new descent-map abort2^-45; no main/descent independence is assumed. The older+10%work-sensitivity claim is removed. The final certificate requires success>.39, not reviewer confidence.

Ledger rows, all in dense-Fp12 units at C=20000: collection116.925; LA115.822475078; generator110.610733127; both-log meter114.085; canonical list91; ALL Schirokauer maps103; filtering/misc80; factor base73; wide lattices70.6; rec7bases73.5; certificate construction50; verification11; g-tree74; g-split95; band-split101.54; newCDFaudits105; main administration74. Sum all rows across all processors once. The wrapper2^66 is additional. The exact first claim is an outward-rounded logarithm covering this sum; the certificate checks it with80-digit arithmetic and Wolfram verifies the rational inputs independently.

Rows include all cold advice, failures, arithmetic, memory/lookup/transfer, processors, randomness, filtering, output and amplification; no amortization or free oracle. Broad fixed factor-base/root construction remains2^73 by the inherited aggregate per-rational-prime price and unconditional pi_Q(B)<1.25506B/logB. Memory remains below2^106bytes: the2^105 core plus wrapper and streamedCDFcounts/thresholds. Samples are streamed, not stored en masse.

## Remaining hypotheses and provenance

H1/H1' are Dickman/Murphy/retention and smoothness-local-law assumptions; H-EA is the survivor-law ECM share; H2 is the uniform-pool work/yield model; H3-algebra, H3-filter/H3c and H3-reg are precisely defined in the main-phase appendix; H4/H4u are kernel/shape/root-density extrapolations and entry-cost/unit shares; H-X is the formal Krylov regularity part of H3-reg, with explicit.01 bad-mass transport and independently derived projection-coin tax; H-FFT is the distributed multiplication extrapolation; H5 is the descendant independence, negative-binomial, missing-log and unit-orbit model; H5p is the marked common-intensity extension below. New GRH and H5n are explicit above. Small experiments only support their actual scale; all extrapolations remain material.

Credit gbotrel explicitly: recover-P base e1e57459-af88-4c69-b5ec-4e86829fd81d and sibling f-only lever ff9f9487-e306-40e3-b6cb-feb0a199515f (after9524d511). Earlier provenance remains the credited lineage in those source packages. This new folder addresses21ec3732's exact main-phase objection, retaining the reproduced count/tail corrections. No prior verdict is a substitute for a missing premise. Protected files and old folders remain untouched.

The core meter117.619455262252 plus wrapper2^66 fits the exact declared score117.619455262253. The certificate recomputes every ledger row and the conservative union success bound0.39063>.39.

Audit arithmetic bounds: the fixed12-dimensional certifiedLLL input has at most381-bit coordinates; the stated Lovasz reduction keeps basis coordinates below2^420 (a loose12-dimensional LLL bound). Multiplying by a ball vector and summing12terms keeps coordinates below2^434. The f/g norm calculations have total degree72/12 and small fixed tower coefficients, so intermediate integer lengths stay below2^15bits using fraction-free resultants. A schoolbook product/division at this width costs below2^26scalar instructions; fewer than2^16such bounded operations, comparisons and control per point fit2^44instructions. Threshold radicands have fewer than2^19bits; bounded integer-root binary search/exponentiation for24thresholds fits2^61instructions per basis.  Node audits therefore cost<=2^28(2^78+2^90+2^61)instructions. Boot whole-ball audits cost<=2^69per basis and at most2^47bases, hence<=2^116instructions. Their combined normalized work is<2^105. This larger row replaces the initial93estimate; it changes only the reported conservative score.

## Bounded survivor prefix (work correction)

A lower yield alone does not bound survivor work. The old QBINS_LO array now DESIGNS a public survivor budget K=ceil(S*(sv_model+pr_model)); it is not asserted to be an actual norm envelope or survivor upper count. Before the sieve's survivor scan, discard every cell whose certified UPPER norm bound exceeds the proposed outer rectangle. The runtime CDF audit evaluates that SAME bound/filter predicate, not a hypothetical exact-norm-only predicate. This makes both retained cofactor bit-size caps valid.

Require each computed archimedean central magnitude M and certified error E to satisfy M>=512E; otherwise discard the cell before thresholding. Then its upper/true ratio is<=513/511, and (513/511)^84<2 proves the total norm-bound gap<1bit for every credited cell. There are43distinct places after conjugation; scaling E by512, comparing and branching costs3instructions/place plus<=64control, within the inherited bound path's1500-to2400instruction margin. Thus the retention proof uses a CERTIFIED1bit margin, not the false universal empirical2e-11 number. The CDF audit includes this guard. Existing floating-bound calibration is evidence only at its stated scope; H4/H5n cover extrapolation.

Retain at most K survivor IDs in uniform order by inside-out reservoir sampling. At j<K append the newID, draw r uniformly0..j, and swap slots r,j only if r<j. At j>=K draw r uniformly0..j and replace slotr only when r<K. Induction gives a uniform ordered Ksample; no uncharged second shuffle is needed. Process that prefix lazily. Generation, addressing and128-bitIDs cost<=64expected instructions/survivor (two64bit uniform words at<=8each, <2expected rejection trials plus control/stores), charged64 per EVERY cell: SCAN10->74. Capped256rejections contribute<2^-125 globally and are in the2^-80error allowance. Allocation, ceiling and one-record endpoint slack fit an additional2^110instructions inside the broad105row. Largest K has log2<70.597;128bitIDs use<2^74.597bytes, below the inherited peak.

For each side and every64bit norm interval[L,U], the certificate lower-bounds its semismooth law by retaining only quadrature cells feasible at L and evaluating rho at ceil(U). It groups16positive weights at their FIRSTindex (the smallest rho), with a2^-30rounding guard. The survivor upper law is evaluated at max(0,L-1). Multiplying the two one-sided ratio maxima plus1 bounds survivor/relation density for ANY retained norm distribution. All62nodes clear raw coverage42; the outward rounding guards discard uncertain boundary terms. This replaces the unsupported survivor-mean use of empirical lower norm bins.

H5p: conditional on a norm profile, survivors have intensity lambda_s*G, G Gamma of shape r; distinct-good marks are independent thinnings of rate p. The full good mean is>=mu_model. Assume weighted distinct-class fraction>=1/64; the interval rate ratio>=42 and model factor .65 give p>=mu_model/K. The 8200 sampled orbits have multiplicity<=62; this supports only that sample, not a production maximum. Marking and production transport remain heuristics.

Under H5p, the full-ball NB failure bound persists. If S>K and the prefix exhausts, abort the whole attack. Coupling to the full process gives extra failure<=E[N]exp(-min(mu_model(1-psi)))<=2^-13, independently recomputed and charged. Low-S full-node failures follow the old recurrence. Conditional marking also gives E[processed]<=1/[p(1-psi)], hence fraction<=1/[mu_model(1-psi)]+1/K. The certificate checks this against every old charged lazy fraction. No invalid conditional-NB-given-S assertion is used. Bad-lattice failures still enter additively.

The exact runtime thresholds use integer U^6 and q (or p^6 at boot), not rounded LOG2U6/LOG2P6 decimals. For hundredth-offset a_i,b_i, use floor100throot(floor((q*CELLS)^600*2^(100a_i))) and floor100throot(floor((q*CELLS)^100*U^600*2^(100b_i))); clear negative-power denominators first. The rounded model logs only enlarge conservative norm caps.

# A direct GRH count bound, with its numerical constants derived

This appendix replaces the previously quoted quantitative Grenie–Molteni theorem. No numerical prime-ideal theorem, zero-count bound, bound on a zeta Laurent coefficient, or unevaluated constant from that theorem is used below. The weaker bound suffices at the benchmark's endpoints. GRH for the two fixed Dedekind zeta functions remains an explicit conditional premise.

Throughout, K is either of the fixed number fields; n=[K:Q]<=72, Delta is its absolute field discriminant, d=log Delta<10000, and r1+2r2=n. Earlier elementary defining-polynomial root bounds give this discriminant inequality, because the field discriminant divides the monic integral polynomial discriminant. Let Lambda_K be the ideal von Mangoldt function and

    psi(x)=sum_{N(p)^m<=x} log N(p),
    psi1(x)=integral_0^x psi(t)dt,
    theta(x)=sum_{N(p)<=x} log N(p),
    pi_K(x)=#{p:N(p)<=x}.

All logarithms in this appendix are natural.

## 1. Standard analytic structure used, with hypotheses

We use the standard analytic continuation, functional equation, Euler product and Hadamard product of a number-field Dedekind zeta function. In the precise normalization needed here, put

    Gamma_R(s)=pi^(-s/2) Gamma(s/2),
    Gamma_C(s)=2(2pi)^(-s) Gamma(s),
    Xi(s)=s(s-1) Delta^(s/2) Gamma_R(s)^r1 Gamma_C(s)^r2 zeta_K(s).

Xi is entire of order one, is nonzero at 0 and 1, and satisfies Xi(s)=Xi(1-s). Its zeros, counted with multiplicity, are the nontrivial zeros rho of zeta_K. Euler-product logarithmic differentiation is valid absolutely for Re(s)>1. Under the stated GRH, rho=1/2+i gamma. Complex conjugation pairs the zeros; real central zeros are included. These structural facts are not a quantitative prime-ideal estimate.

For clarity, the numerical consequences needed from the Hadamard product are derived, not quoted. Write

    Xi(s)=exp(A+Bs) product_rho [(1-s/rho) exp(s/rho)].

The order-one product gives sum_rho 1/(1+|gamma|^2)<infinity. On the real axis, pair conjugate zeros, so all real-part sums below converge absolutely. Logarithmic differentiation at 0 and 1, followed by the functional equation, gives

    B=-sum_rho Re(1/rho),
    Xi'(sigma)/Xi(sigma)=sum_rho Re[1/(sigma-rho)]

for real sigma outside the zeros. Indeed Xi'/Xi(0)=B; Xi'/Xi(1)=B+2 sum Re(1/rho), while the functional equation makes these negatives. This fixes the otherwise ambiguous linear term in the product.

## 2. A convergent zero sum bounded at s=2

Let L2=Xi'(2)/Xi(2). Differentiating the displayed definition of Xi gives

    L2=3/2+d/2+r1 Gamma_R'(2)/Gamma_R(2)
         +r2 Gamma_C'(2)/Gamma_C(2)+zeta_K'(2)/zeta_K(2).

The Euler product makes its last term nonpositive. Both gamma terms are negative:

    Gamma_R'/Gamma_R(2)=-(gamma_E+log pi)/2<0,
    Gamma_C'/Gamma_C(2)=1-gamma_E-log(2pi)<0.

Here 0<gamma_E<1, 2<e<3 and 3<pi<4 suffice; the gamma identities follow from the recurrence and duplication formulas for Gamma. Consequently

    0<L2<=d/2+3/2.

No estimated distribution of zero ordinates enters this step. For u=gamma^2>=0,

    1/sqrt[(u+1/4)(u+9/4)] <= 3/(u+9/4),

because u+9/4<=9(u+1/4). Thus the absolutely convergent zero weight obeys

    S=sum_rho 1/|rho(rho+1)| <= 2L2 <= d+3 <10003 <101^2.       (Z)

A second termwise comparison gives

    sum_rho Re[1/(1-rho)] <=3L2,

using (1/2)/(u+1/4)<=3(3/2)/(u+9/4). Hence |B|<=3L2.

## 3. The smoothed explicit formula and its residues

This identity is included to identify every lower-order term and avoid importing a bound on a zeta Laurent coefficient. For x>1,

    psi1(x)=x^2/2-sum_rho x^(rho+1)/[rho(rho+1)] +R(x),         (EF)

where

    R(x)=-x r_K+r'_K-a(x log x-x)+r2(log x+1)
         -(r1+r2) f1(x)-r2 f2(x),
    a=r1+r2-1,
    f1(x)=sum_{j>=1} x^(1-2j)/[2j(2j-1)],
    f2(x)=sum_{j>=2} x^(2-2j)/[(2j-1)(2j-2)].

The coefficients r_K and r'_K are defined by the Laurent expansions

    zeta_K'/zeta_K(s)=a/s+r_K+O(s)             near 0,
    zeta_K'/zeta_K(s)=r2/(s+1)+r'_K+O(s+1)   near -1.

Here is the contour derivation of (EF). On Re(s)=2 the Euler-product series and the integrable s^-2 kernel give, by absolute convergence,

    psi1(x)=-(1/(2pi i)) integral_(2-i infinity)^(2+i infinity)
                  (zeta_K'/zeta_K)(s) x^(s+1)/[s(s+1)] ds.

The elementary Mellin kernel integral of y^s/[s(s+1)] is 1-1/y for y>1 and 0 for 0<y<=1. Multiplying by x and summing Lambda_K produces sum Lambda_K(I)(x-NI)_+=psi1(x), so there is no endpoint convention problem.

Shift this contour to negative half-integers, avoiding zeros on the horizontal edges. The usual limiting shift is justified without any numerical zero-count constant: the functional equation and Stirling bounds, together with the absolutely convergent Euler product at real part 2, imply via Jensen's formula that a bounded-height strip contains O_K(log(T+3)) zeros. Choose edge heights in successive unit intervals at distance >>_K 1/log(T+3) from their ordinates. The Hadamard logarithmic derivative there is O_K(log^2(T+3)); the horizontal integrals tend to zero because of the quadratic denominator. On left edges the functional equation and x^(Re(s)+1), with x>1, make the integrals tend to zero. These qualitative limiting statements establish an identity; their implicit constants do not enter any numerical bound below. Equivalently the same shift may first be performed on a finite rectangle and then passed through this convergent sequence.

The enclosed residues are explicit:

- At 1, the simple pole of zeta_K gives x^2/2.
- A nontrivial zero rho of multiplicity m contributes -m x^(rho+1)/[rho(rho+1)]. The resulting series is absolutely convergent by (Z).
- At 0, expand x^(s+1) and 1/[s(s+1)]; the residue is -x r_K-a(x log x-x).
- At -1, the residue is r'_K+r2(log x+1).
- At -2j, j>=1, the trivial-zero multiplicity is r1+r2 and contributes -(r1+r2)x^(1-2j)/[2j(2j-1)].
- At -(2j-1), j>=2, the trivial-zero multiplicity is r2 and contributes -r2 x^(2-2j)/[(2j-1)(2j-2)].

These are precisely the terms displayed in R(x). The two trivial-zero series and their derivatives converge uniformly for x>=2.

## 4. Elementary bounds on the Laurent terms and R'

The finite part of the logarithmic derivative of Xi at 0 gives

    r_K=B+1-d/2-g0,
    g0=-r1(gamma_E+log pi)/2-r2(gamma_E+log(2pi)).

Since |B|<=3L2 and |g0|<=2n using 0<gamma_E<1, log pi<2 and log(2pi)<3,

    |r_K| <=2d+6+2n.                                        (R0)

For completeness, the same calculation at -1 gives

    r'_K=Xi'(-1)/Xi(-1)+3/2-d/2-gminus1,
    gminus1=r1(1-gamma_E/2-log2-(log pi)/2)
                 +r2(1-gamma_E-log(2pi)).

The functional equation gives Xi'(-1)/Xi(-1)=-L2, and |gminus1|<=2n. Thus |r'_K|<=d+3+2n. This last constant cancels in every finite difference, so its magnitude is not required for the bound.

Differentiate the explicitly displayed R(x). For x>=2,

    |f1'(x)| <=1/[2(x^2-1)],
    |f2'(x)| <=1/[3x(x^2-1)].

These follow by differentiating their positive term series and dropping denominator factors. Therefore

    |R'(x)| <= |r_K|+n log x+n
              <=2d+6+3n+n log x.                            (R1)

For any x in [2^44,2^225], put h=202 x^(3/4). We have h/x<=202/2048<1/10. All intermediate points x-h through x+h exceed 2 and are below 2^226. Since log2<7/10, (R1) is at most

    20000+6+3*72+72*(226*7/10) <32768=2^15.                 (R2)

All constants in this inequality are rational or integral.

## 5. Unsmoothing: a uniform numerical psi bound

Nonnegativity of Lambda_K makes psi nondecreasing. Thus

    [psi1(x)-psi1(x-h)]/h <= psi(x)
           <=[psi1(x+h)-psi1(x)]/h.

Insert (EF). For the upper difference, (Z) gives

    |sum_rho[(x+h)^(rho+1)-x^(rho+1)]/[h rho(rho+1)]|
       <= S[(x+h)^(3/2)+x^(3/2)]/h
       <= (221/100)*101/2 *x^(3/4).

The last step uses x+h<1.1x and (1.1)^(3/2)<(1.1)^2=121/100. The main-term difference is x+h/2, with h/2=101x^(3/4). The lower difference has the smaller zero bound 101x^(3/4). By the mean value theorem and (R2), either difference of R contributes at most2^15. Since x^(3/4)>=2^33,

    101+(221/200)*101+2^15/2^33 <213.

Consequently, with every constant now justified,

    |psi_K(x)-x| <=213 x^(3/4)
                       for 2^44<=x<=2^225.                 (P)

This is intentionally weaker than a sharp GRH prime-ideal theorem, but strong enough here. No Trudgian zero-count constants or Grenie–Molteni quantitative estimates occur in its proof.

## 6. A direct endpoint bound for pi_K

Let T=2^44 and li_3(x)=integral_3^x dt/log t. First,

    0<=psi_K(x)-theta_K(x)<=n sqrt(x) log x.                 (PP)

Each proper ideal-power term lies over a rational prime p<=sqrt(x). There are at most n prime ideals over p; for any one ideal, the number of exponents times log N(p) is at most log x. Counting the underlying rational primes by sqrt(x) proves (PP), without a prime-counting estimate. The function (log x)/x^(1/4) decreases once log x>4. Hence for x>=T its contribution divided by x^(3/4) is at most72*(44 log2)/2^11<72*44*(7/10)/2^11<2. From (P),

    |theta_K(x)-x|<=256 x^(3/4), T<=x<=2^225.               (TH)

The exact Stieltjes partial-summation identity, using the atomic mass log N(p) at each prime-ideal norm, gives

    pi_K(x)-pi_K(T)-integral_T^x dt/log t
      =[theta_K(x)-x]/log x-[theta_K(T)-T]/log T
       +integral_T^x [theta_K(t)-t]/(t log^2 t)dt.           (PS)

At endpoints use (T,x] for both differences, equivalently subtract the two cumulative identities with <=. This handles atoms exactly.

The low boundary is bounded elementarily: pi_K(T)<=nT<=72T because at most n ideals lie over each rational prime<=T; 0<li_3(T)<T because log3>1. Thus |pi_K(T)-li_3(T)|<73T. Also 0<=theta_K(T)<=psi_K(T)<=T+213T^(3/4)<2T, as213<2^11. Therefore |theta_K(T)-T|/logT<=T. Apply (TH) inside the remaining integral and use log t>=logT:

    |pi_K(x)-li_3(x)|
       <74T+256 x^(3/4)/log x
             +[256/(logT)^2] integral_T^x t^(-1/4)dt
       <=74T+256 x^(3/4)[1/log x+4/(3(logT)^2)].             (PC)

Every constant and integral in this step is explicit. In particular the negative lower endpoint of the integral was dropped in the direction enlarging the bound; no asymptotic series for li is substituted.

The positive series log2=2 sum_(j>=0)1/[(2j+1)3^(2j+1)] gives log2>2/3+2/81=56/81 and log2<2/3+1/36<7/10. For the upper inequality, bound1/(2j+1)<=1/3 forj>=1 and sum the resulting geometric tail. Thus69log2>47 and44log2>30. Also e<3 follows from k!>=2^(k-1) for k>=2, strictly for some terms, in the exponential series; hence log3>1.

For each integer69<=b<=225 put C=ceil(3b/4). Then C>=52, (2^b)^(3/4)<=2^C andT/2^C<=1/256. Divide (PC) by2^C. Its coefficient is strictly below

    74/256+256*(1/47+4/2700)<8.

Consequently the derived endpoint bound is

    |pi_K(2^b)-li_3(2^b)| < E_b,
    E_b=2^(ceil(3b/4)+3), 69<=b<=225.                       (COUNT)

Atb=69 this isE_69=2^55. For any two such endpoints the interval count differs from integral_(2^a)^(2^b)dt/logt by at mostE_a+E_b, by the triangle inequality. Higher-degree ideals, ignored polynomial-bad primes and symmetry quotients used by the attack must then be counted separately. This absolute prime-ideal count theorem has now been derived here from the displayed analytic structure and elementary inequalities, rather than quoted from Grenie–Molteni or another quantitative counting theorem.

## 7. Details of the limiting contour shift (no quantitative zero-count input)

This paragraph explains the qualitative estimates used to pass from the Perron integral to its residues. Its field-dependent constants disappear in limits and are not used to price the algorithm or bound psi.

On Re(s)>=2 the Euler product gives both |zeta_K(s)|<=zeta_K(2) and |zeta_K(s)|>=1/zeta_K(2): for the lower bound apply |1-z|<=1+|z|<=1/(1-|z|) to each local factor. The functional equation, Stirling's formula, and this right-half-plane bound give polynomial growth in |Im(s)| on every fixed strip to the left. The Phragmen–Lindelof principle in a strip (equivalently, a damped rectangle maximum-principle argument), applied after removing the pole at1, gives the same polynomial-growth statement between the left and right edges.

Apply Jensen's formula to Xi on disks centered at2+it, with outer radius4 and inner radius3, for t sufficiently large. The lower Euler-product bound at the center and the polynomial-growth bound throughout the disk, together with Stirling for the gamma factors, show that the logarithm of the maximum modulus divided by the center modulus is O_K(log(t+3)). The exponential gamma decay changes by only a bounded field-dependent factor on this fixed disk. Jensen therefore bounds the number of zeros in the inner disk by O_K(log(t+3)). The disk contains the portion of the critical strip with ordinates between t-1 and t+1. Thus every bounded-height strip contains O_K(log(t+3)) zeros; summing such strips also gives O_K(T log(T+3)) zeros up to heightT.

For each sufficiently large integer j, remove from [j,j+1] intervals of radius c_K/log(j+3) around the zero ordinates in [j-2,j+3]. Since there are O_K(log(j+3)) such ordinates, choosing a fixed sufficiently small positive c_K leaves a height T_j. At that height the nearest-zero terms of the Hadamard logarithmic derivative contribute O_K(log^2(T_j+3)). For ordinates at distance between k and k+1, the local zero count is O_K(log(T_j+k+3)); summing the reciprocal-distance terms, with the Hadamard 1/rho regularizer, gives another O_K(log^2(T_j+3)). Ordinates beyond2T_j give O_K(log(T_j+3)) by the convergent second-power tail. This proves zeta_K'/zeta_K(s)=O_K(log^2(T_j+3)) on any fixed-width horizontal edge. The gamma logarithmic derivatives only add O_K(log(T_j+3)).

Fix a negative half-integer left edge -M-1/2, M>=2, and first let T_j tend to infinity. On each horizontal edge, x^(Re(s)+1) is bounded for this fixed M and x, while |s(s+1)| is of order T_j^2. Its integral therefore tends to zero. On the left vertical edge, the logarithmic derivative of the functional equation relates zeta_K'/zeta_K(s) to its absolutely convergent Euler-product derivative at1-s, plus gamma derivatives and d. The edge stays a fixed half-integer distance from every trivial zero. Its magnitude is O_K(log(M+|Im(s)|+3)). The integral is consequently

    O_K[x^(-M+1/2) integral_R
                   log(M+|t|+3)/(M^2+t^2) dt]
      =O_K[x^(-M+1/2) log(M+3)/M],

which tends to zero as M tends to infinity because x>1. The nontrivial-zero residues already converge absolutely by (Z), and the trivial-zero residues form the two convergent power series displayed in (EF). The residue theorem now gives precisely (EF), rather than an identity with an uncontrolled truncation error. Only the familiar analytic continuation, functional equation, gamma identities, Hadamard factorization, Jensen formula, Stirling formula and the Phragmen–Lindelof strip principle were used; no numerical estimate from another prime-ideal theorem was substituted.


## Sharper average count error at the factor-base endpoint only

The preceding direct proof also permits a sharper endpoint calculation for B=2^69, without changing its generic E_b bound for repetitions. This refinement uses the two degrees and the already proved polynomial-root discriminant bounds separately:

    n_f=72, d_f<72*71*(9/5)=9201.6;
    n_g=12, d_g<4488*(7/10)=3141.6.

The derived zero-sum inequality S_K<=d_K+3 therefore gives S_f<96^2 and S_g<57^2. Repeat exactly the finite-difference proof with A_f=96 and A_g=57 instead of101, and h=2A_K*x^(3/4). For x>=T=2^44 the ratios h/x are at most192/2048 and114/2048, respectively, both below1/10. The same uniform remainder derivative bound2^15 applies. The resulting upper coefficients, including that remainder, are bounded by

    (421/200)*96+2^-18<203,
    (421/200)*57+2^-18<121.

The lower differences have smaller coefficients, so throughout [T,2^225],

    |psi_f(x)-x|<203*x^(3/4),
    |psi_g(x)-x|<121*x^(3/4).

The elementary proper-power correction is below2*x^(3/4) on f and below1*x^(3/4) on g: its maxima are bounded by72*44*(7/10)/2^11<2 and12*44*(7/10)/2^11<1. Thus the theta coefficients205 and122 are valid. The low boundary in the previously proved partial-summation inequality is (n_K+2)T; retain74T for f and14T for g. This proves the following bound on the AVERAGE of the two absolute count errors:

    (|pi_f(B)-li_3(B)|+|pi_g(B)-li_3(B)|)/2
      <44T+(327/2)*B^(3/4)
                       *[1/log B+4/(3(log T)^2)].

Use logB>47 and logT>30 as proved above. Also B^(3/4)=2^51*8^(1/4)<(27/16)*2^51, since(27/16)^4>8. All remaining arithmetic is rational, and the displayed right side is less than

    44*2^44+(327/2)*(27/16)*2^51*(1/47+4/2700)
      =17523154306377711616/1175
      <15000000000000000.

This last constant bounds the average, not each field's individual error. For sigma-quotient factor-base counting, the average is exactly what is required. Therefore one may use

    U_cap=ceil[B/(69 log2-11/10)+2^18+17523154306377711616/1175+84]
         =12647734742496954934,

with the same exact atanh lower interval for log2. The84 allowance may be retained conservatively. For repeat probabilities on an individual field, continue to use the generic proved endpoint radius E_b=2^(ceil(3b/4)+3), then account for non-degree-one ideals, bad-prime subsets and free sigma orbits as before.

Every numerical step in this refinement is independently checked by `counts16.py` and the recorded Wolfram code/raw under `wolfram-field69.wl`/`wolfram-field69.txt` in stats16.b85.


The exact-rational endpoint refinement uses17523154306377711616/1175 directly, rather than rounding it to15*10^15. A distinct lower bound on the ACTUAL factor universe is proved in the main-phase appendix and checked separately by counts16.py. It is used only in failure probabilities.

## Main-phase contracts, complete probability space, and construction

This appendix replaces all inherited shorthand H3/H3c/H-X arguments. These are explicit conditional production-model hypotheses, with their actual support and extrapolations separated from deductions. Let U be the ACTUAL sigma-orbit factor-base size, Ucap its proved work upper bound, and Ulo=12206442875139313724 its independently proved lower bound. In particular Ulo is not inferred from Ucap. The count appendix gives

Ulo=floor((2^69-3)/(69*7/10)-15*10^15-72*2^35-2^20-84)>10^19.

It follows by retaining only simple absolute-degree-one free sigma orbits, subtracting the average endpoint error, higher-degree ideals and polynomial-bad primes. All these are a subset of the actual universe. The extra84 is conservative. Upper and lower universes have distinct roles below.

### Fixed algebraic prerequisites, not projection probabilities

H3-algebra requires r prime to both field class numbers and squarefree reductions of BOTH defining absolute power-basis polynomials modulo r. The latter is checked by gcd(F,F')=gcd(G,G')=1 modulo r; field-unramifiedness alone would not exclude a bad power-basis index. Signature checks give Kf=(0,36), Kf^sigma=(0,18), Kg=(2,5), Kh=(2,2), so Dirichlet odd unit ranks are18 and3. These facts and squarefreeness are recomputed by main_phase.py. Neither signatures nor the gcd test computes class numbers or proves a nonzero unit regulator.

Work in each etale power-basis algebra modulo r^2. Choose the UNIVERSAL exponent E_K=product_(d=1)^n(r^d-1), not an unspecified fast factorization-based exponent. It is divisible by every residue-field multiplicative order and nonzero modulo r. For an r-adic unit alpha define lambda(alpha)=(alpha^E_K-1)/r modulo r, coefficientwise. Since(1+ra)(1+rb)=1+r(a+b) modulo r^2, lambda is additive and annihilates rth powers. Sigma commutes with it; use lambda^-=(lambda-sigma(lambda))/2. H3-algebra assumes this full odd map is injective on odd units modulo rth powers. This is the precise unverified regulator condition, not an inference from a planted-kernel toy.

Sample eighteen independent uniform F_r linear functionals on the72 coefficients for f and three on the12 for g, and apply them to lambda^-. Conditional on full odd injectivity, their restrictions to unit images are uniform18-by-18 and3-by-3 matrices. Their determinant degrees give failure at most21/r. The same polynomial-root argument proved below applies. There are1332 sampled coefficients. This is a coin error; full unit injectivity and the class-number conditions remain fixed-field conditional premises. Accepted B-smooth relation elements are invertible modulo r because B<r.

These fixed algebraic conditions imply compatible virtual logs, and the following exact anchor proves their retained vector nonzero without adding compatibility to a random-matrix hypothesis. Let p be the production base prime, h=(p^12-1)/r and k=(h modr)^-1 modr. Exact arithmetic gives h modr=13860579304294609274647100502445671956830243711357126442998 and k=24831094541369581385436675797109464480232629635954898292486708404843532461745. Define pi(z)=(z^h)^k. Its image lies in GT, hk=1 modr makes pi the IDENTITY on GT, and pi^2=pi. Since r divides p^6+1, h is a multiple of p^6-1; pi kills nonzero F_(p^6) multipliers in boot lifts.

Take alpha0=1+x in Kf. Exact resultants give Norm(alpha0)=1399, a prime below B. The prime is unramified (the absolute polynomial discriminant is288 mod1399), and gcd(h(y),f(-1,y)) mod1399 is y+40. Thus (alpha0)=I=(1399,y+40,x+1), a single relative-degree-one ideal of absolute norm1399; sigma(I), with x-1, is distinct. Small primes remain factor-base ideals even where their sieve-score mass is omitted. Independent exact Python and Wolfram field arithmetic prove omega=pi(eval_f(alpha0)) is nonidentity, omega^r=1 and sigma(omega)*omega=1. For example its x coefficient is1823330539042453828028069116870546032896926197478728071363869126167756345638209516204631600716902807865027915471103, not zero. Hence omega generates GT. All resultants, prime-ideal, exponent and finite-field code/raw outputs are retained with the proof.

Define the TRUE virtual logarithm by L_K(alpha)=log_omega(pi(eval_K(alpha))). This is an existence definition, not a free DLOG oracle. It is already sigma-odd, because p^6=-1 modr and pi commutes with sigma. In particular L_f(alpha0)=1 exactly. Do NOT apply an additional unnormalized L-Lsigma subtraction; if odd-projection notation is used its convention is (L-Lsigma)/2. On GT, L is the ordinary logarithm to omega: blinding z*omega^b adds b, not hb, and boot lifts differing by an F_(p^6) multiplier have the same logarithm.

Unit-map injectivity and the successful random unit-coordinate maps let ell_K agree with L_K on odd units. Then delta_K=L_K-ell_K dot SM_K annihilates all units. Choose a class-group exponent h_K prime to r and generators beta_J of J^h_K. Set v_J=h_K^-1 delta_K(beta_J). Unit independence makes this well-defined/additive. For smooth alpha, alpha^h_K/product beta_J^e_J(alpha) is a unit; applying delta and dividing by h_K proves L_K(alpha)=sum e_J(alpha)v_J+ell_K dot SM_K(alpha). No class exponent, unit basis or unknown log is computed as free advice: existence is proved here, and the charged matrix solve recovers the coordinates.

Orient the stored protected orbit coordinate as I, with sigma(I) carrying the negative coefficient in every rewritten row. Keep that semantic I-orbit coordinate and ALL21 Schirokauer coordinates throughout filtering: do not purge, eliminate, or merge any of these coordinates away. Other exact eliminations preserve the restriction of the true advice. With the signed row convention (e_f,-e_g,SM_f,-SM_g) and advice v*=(v_f,v_g,ell_f,ell_g), equal common-field evaluations make EVERY checked row vanish on v*. The protected anchor gives eta(v*)=v*_I+ell_f dot SM_f(alpha0)=L_f(alpha0)=1. Thus the retained v* is provably nonzero; this does not assert that v*_I alone is nonzero. If the selected square has rank D-1, its kernel is automatically this true line. For any extracted nonzero w=c*v*, eta(w)=c!=0; dividing w by eta(w) normalizes the true advice. The anchor map packet, omega construction, inverse/normalization, checks and counters are charged below. The construction and normalization use no matrix-compatibility heuristic.

### Actual row ensemble, gates, square selection and production moments

The candidate universe is the fixed canonical primitive region in config.json. Rational factor-base primes and root labels are constructed by the charged cold preprocessing. Every accepted relation is exactly checked: factor-base membership, rational primality by the constructed prime table, factorization/product identity and ideal roots. Probable-prime tests may screen candidates but cannot authorize a false row. Curves use independent per-candidate random streams. All tests, rejected trials and local counters are included in the work model. Each candidate's stochastic cofactor work is capped at2^72 instructions. R denotes the retained-row count of the full virtual collection, after all tests and that local cap. Early actual aborts only reduce the work relative to this process.

The explicit H1/H2 true-yield premise is

1.087 Ucap <= E[R] <=2^64, with independent nonidentical Bernoulli per-candidate retention, hence Var(R)<=E[R].

The supplied G mean and two-SE summaries are evidence for the mean premise, not an algorithmic confidence factor. The checker recomputes their values. Enforce Rneed=ceil((679/625)Ucap) and abort on fewer rows; it checks exactly Rneed^40>=2^2543, so Rneed>=2^63.575, the actual descent input, and also Rneed>=1.08U. Abort before mapping row2^65+1. This cold retained-row cap is Rmax=2^65, enforced before its map packet. With mu0=1.087Ucap, the lower-tail tax is mu0/(mu0-Rneed)^2 by Chebyshev and monotonicity. The upper-tail tax is at most2^-64, because Rmax-E[R]>=2^64 and Var(R)<=E[R]<=2^64. These are separate from filtering concentration. The old weaker1.081Ucap premise did not guarantee the descent relation gate and is not used here.

Let W be total collection instructions before its hard meter, including deterministic scan/layout and the independent, locally capped cofactor costs. H2 is the quantitative production premise E[W]<=20000*2^116.92454. The empirical work sums, squares and1.00843 allowance support this finite-population extrapolation; they do not prove it. Independent Xi in[0,2^72] give Var(W)<=2^72 E[W]. The NEW collection meter is20000*2^116.925. If mu_max=20000*2^116.92454, then mu_max>2^130 and the gap exceeds mu_max/4096, because(116.925-116.92454)*(2/3)>1/4096 and exp(t)>=1+t. Thus the work-abort probability is less than2^-34. Fixed scans use bounded blocks/descriptors; stochastic counters update at curve/batch boundaries inside charged stage/control allowances, not with uncharged per-cell overhead. The planned69-bit ECM schedules at cofactor sizes241/310 have full-redo costs far below the local cap; the mean premise nevertheless applies AFTER all stated caps.

Use ALL retained rows up to the hard row cap, in fresh uniform record order, in singleton purge, clique removal and the specified Markowitz merge. Filtering and shape construction have a hard2^79-unit meter, leaving the other half of the existing misc80 row for its other work. Row/column keys and coefficients are preserved exactly. Keep the protected semantic I-orbit column and the21 dense coordinates throughout; only other columns may be purged or eliminated. Let N_* be retained ROW count, a conservative dimension statistic, and set N_*=U on a filtering-work or shape failure. Shape includes coefficient representation, enough rows for the augmented unknowns, and the selected square's sparse mean row weight<=301. If ideal-column count is N and D=N+21, require at least D reduced rows. Uniformly shuffle reduced rows, choose the first D as square A, and retain the rest for verification. No independent rows are selected by free Gaussian elimination. Canonically sort columns by ideal key then append21 dense columns. There is no fresh column-support shuffle.

H3-filter/H3c are moment assumptions on THIS ACTUAL PROTECTED-COLUMN production variant, conditional on the preceding collection gates, including the abort encoding. Protecting I may alter an elimination cascade; no bound that its effect is merely +1 column is asserted:

E[N_*/U]<=.03329, Var(N_*/U)<=1/U.

Abort when N_*>.03330U. Therefore its probability is at most10^10/Ulo. The work/memory cap uses Dcap=ceil(.03330Ucap)+21. The row statistic dominates ideal-column count on usable overdetermined systems; the extra21 is explicitly retained. The shipped ideal-only surrogate at excess1.08 is support, not an identical production ensemble at the actual Rneed gate. Sixteen/eight seeds at scales2^19/2^20 have mean+2SE below.03329; sample sd*sqrt(U_sample)<.644 on the shipped adequately sampled scales. Changed-norm paired data also support the chosen envelope. These facts do not prove the production moments, dependent-row transport, shape-abort mass or a10^-9 tail.

H4u is likewise stated on the ACTUAL SELECTED square A. For fractions f_c of its ORIGINAL sparse nonzeros with absolute coefficient c, c=1..4, the measured five-vector price functional is

k5=9.2 f1+18.4 f2+27.6 f3+36.8 f4+40(1-f1-f2-f3-f4).

Assume E[k5]<=9.966 and Var(k5)<=256/U, conditional on the preceding shape gate. Abort if k5>9.97. Chebyshev bounds this event by16*10^6/Ulo. The20 preregistered scaled runs support the moments: all sample variances times U_sample are below256 and the two larger-scale mean+2SE values are below9.966. These statistics are recomputed. Transport to the selected production square is unproved and material. A mere86% share of +/-1 entries would permit k5=13.512 and is not the claimed justification. On passing, the same coefficients are reused in every matrix pass; no per-iteration independence is needed. The k5 bound also implies at least(18.4-9.97)/9.2>.86 unit entries, sufficient for the more conservative single-vector entry price.

The sum of shortage, row overflow, collection work, filter-size and k5 taxes is recomputed exactly and is below10^-9. The whole declared P_CAP=10^-9 now covers these specified events, not a two-SE confidence statement.

### Matrix regularity is an explicit probability hypothesis

For the actual selected A after the preceding passes, define G_reg ONLY by rank(A)=D-1 and the two formal Krylov ranks specified next equal D. The compatible nonzero true line now follows from the fixed class/unit construction, checked relation rows and protected eta(v*)=1 anchor; it is not a further random regularity predicate. H3-reg is the explicit PRODUCTION TRANSPORT hypothesis P(not G_reg | preceding passes, fixed H3-algebra)<=.01. This bound is not inferred from rank dimensions, unit shares or an assumed independent-entry matrix. It is a new quantitative assumption about the actual arithmetic/ECM/filter/row-selection ensemble. Unlike an assertion that random A is always regular, its possible failures are charged.

Fresh preregistered evidence is600 synthetic matrices with D=80,m=72,n=5 over the actual scalar field, seeds40951573..40952172, with a planted row-sum-zero line, sparse rows, the fixed coordinate partition and independently sampled X/Y. Every rank(A), rank(O), rank(R), H and H' identity test passed. The checker recomputes all recorded summaries and a disjoint fresh witness. Under an independent surrogate-draw model, a bad mass.01 gives zero-failure probability(.99)^600<.05. This supports only the synthetic planted-kernel family. It does not prove production rank/formal-rank transport, class numbers or unit injectivity. True-line compatibility/nonzeroness is instead established by the separate exact protected-anchor construction. The raw records, fixed seeds, frozen generator and limitations are supplied. No confidence level is multiplied into attack success.

### Fixed-support projection and its derived coin error

Fix m=72,n=5,a=ceil(D/m),b=ceil(D/n),L=a+b. Partition canonical coordinates by block(j)=floor(72j/D). X is D-by72 with X[j,block(j)] independently uniform in{0,...,2^128-1} and all other entries zero. Y is D-by5 with independent uniform F_r entries. Coins are fresh after A is fixed. Define O_a by stacking X^T A^i for0<=i<a, and R_b=[Y,AY,...,A^(b-1)Y]. G_reg requires their FORMAL ranks over F_r(X) and F_r(Y) to be D. Rank(A)=D-1 alone is insufficient, as diag(0,1,...,1) demonstrates. The premise is nonempty: distinct eigenvalues including zero give block Vandermonde witnesses for the fixed partition and an appropriate five-column Y. This example is not a proof about production matrices.

Choose formal nonzero D-minors. Each has degree at most D in its independent X or Y variables. The finite-set polynomial bound is proved by induction: view a polynomial as univariate in its last variable, bound the chance its leading coefficient vanishes inductively, and otherwise use the univariate root bound. Total degree d gives at most d/|S| zero probability. Hence, conditional on G_reg,

P(rank O<D or rank R<D)<=D/2^128+D/r.

Uniform field generation rejects255-bit candidates>=r and aborts after256 attempts. Since r>2^254, the combined sampler tax for Y and the1332 map coefficients is at most(5D+1332)*2^-256. Add21/r for unit-coordinate rank. With D<=Dcap the sum is below2^-66, checked exactly. There is one projection/solve attempt; any sampler or reconstruction failure aborts and is charged. No probability of a bad support permutation is silently covered: the partition is fixed and G_reg explicitly includes its formal ranks.

Compute S_i=X^T A^iY for0<=i<L. H with blocks S_(i+j),0<=i<a,0<=j<b, satisfies H=O_a R_b, so has rank D on a good draw. H' with blocks S_(i+j+1) equals O_a A R_b and has rank D-1. Therefore c in kerH' but not kerH gives nonzero R_b c in kerA. The algorithm uses the following structured version, not dense D-cubed elimination.

Let Phi(p)=sum A^jY p_j on F_r[t]^5. Right spanning makes Phi onto. Its kernel has a column-reduced polynomial basis F(t) with column degrees<=b and determinant degree D: eliminate spanning chains Y,AY,... by their first dependence; the quotient has dimension D and t acts as A. For degree<=b, the finite equations sum S_(i+j)p_j=0 are equivalent to Phi(p)=0 because O_a is injective. Reverse degree-b polynomials, set G(t)=[S(t),-I_72], and use approximant order L with shift(0 repeated5,1 repeated72). Selecting the five reduced columns of shifted degree<=b encodes degP<=b,degQ<b and yields the true relation basis, not a projected spurious relation.

The approximant algorithm is explicit: order1 uses Gaussian elimination with minimal-shift pivots and multiplies pivot directions by t. At ell=ell1+ell2 compute B1 for the first ell1 terms, residual G2=(G B1)/t^ell1 modulo t^ell2, recursively compute B2 with the new degree shift, and return B1B2. Every first-stage approximant is uniquely B1q; the remaining equations are exactly G2q=0. This proves the recurrence and basis correctness. F(0) has nullity one because the quotient modulo t has dimension dim(coker A)=1. For c!=0 in kerF(0), set w=sum_(j>=1)A^(j-1)YF_jc. Then Aw=0. If w=0, Fc/t belongs to kerPhi, so Fc/t=Fq and q=c/t, impossible for polynomial q. Thus w!=0 without any semisimple-zero assumption. Check w!=0,Aw=0 and all original/dropped relations after exact back-substitution; failed checks abort.

### Explicit work, including new map and administration charges

One map PACKET computes both fields' universal maps and all21 coordinates (a one-field packet costs no more). The universal E has fewer than255*n(n+1)/2<2^20 bits. Binary powering uses fewer than2^21 ring products. A schoolbook product plus monic reduction at n<=72 uses fewer than2n^2<2^14 scalar operations modulo r^2; eight-limb Montgomery loops give the conservative2^14-instruction price per scalar operation. Both field maps, coefficientwise divisions, odd projection, random functionals, global counter increment/comparison and control fit2^51 instructions per packet.

Cold relation packets have hard cap2^65, enforced before packet2^65+1. In descent, defer maps until a candidate is semismooth and deduplicated; do not map every sieve survivor or an already rejected candidate. On the actual conditional descent law, each of the EN expected ordinary nodes and the two boot nodes needs at most1/(1-psi)<=5 map packets in expectation, using the checked uniform psi<=.8 bound and the lazy candidate order. This counts packets on failed candidate branches too. With EN<2^17, E[descendant packets]<=5(EN+2)<2^20. Enforce a separate descendant-packet cap2^65; Markov gives an ADDITIONAL failure tax2^-45. This is a global counter abort, not uncharged repeated work, and it is added separately in the success union below.

At most2^65 cold plus2^65 descendant plus ONE alpha0 anchor packet are ever computed, so total packets<=2^66+1. Absolute power-basis conversion y(x), reduced polynomials, exponent construction and other fixed setup cost below2^70 instructions are included. Therefore total map/setup work is at most(2^66+1)*2^51+2^70 instructions, below20000*2^103. This REPLACES the old unsupported SM88 row by103. The already conservative packet price includes every counter update; no unmetered map can bypass a cap. Omega construction uses fewer than2^14 tower products at the stated implementation price, below2^50 instructions and inside the cold wrapper. Advice normalization and its checks fit administration74. No class-number/unit-regulator oracle is used.

On a passing selected matrix the sparse row price is c=301*9.97+1215+21*97=6252.97. The existing Krylov/Horner/projection expression is applied conservatively to Dcap, with the extra21 coordinates included. Supported X gives nDL products at60, not72 times that count. A128-by255-bit product accumulated over fewer than2^59 coordinates fits seven limbs (<2^442); each of mn block coefficients reduces at420. Horner recovery splits into five chunks, keeps the20 additional shifted Y checkpoints, and evaluates five chains in lockstep, giving the charged L D*(c+5*97+420) term. At most four endpoint/checkpoint steps are charged separately.

For the generator k=77,C=97+420=517,Lcap=(Dcap+21)(1/72+1/5)+2. Padding to a power of two P gives a scalar polynomial-product bound3C P^log2(3): leaf products plus at most16P^log2(3) additions fit because517+16*64<=3*517. At each recursive basis node, two polynomial matrix products plus residual/pivot bookkeeping fit8C k^3 P^log2(3); the recursion ratio is2/3 and Gaussian/inversion bases are covered by4C k^3. Total<=32C k^3 P^log2(3)<=96C k^3 Lcap^log2(3). Since log2Lcap>24, this is at most4C k^3 Lcap^log2(3) log2Lcap, the charged generator expression. Small actual matrices are bounded by this large cap, not by extrapolating the production inequality backwards.

A NEW administration row74 covers at most2^88 instructions for all capped original-row checks (at most2^65 rows,2^22 each),2^81 for random words/order/index preparation, and2^76 for ceiling/end-step, anchor eta evaluation, vector normalization and Aw checks. Their sum is below20000*2^74. This row is separate from misc80. The other misc work and the2^79 filter meter fit misc80. All retained-row verification, random attempts, failed checks, owner transfers and output work are included. Memory remains below the existing2^106-byte cap: added25 vectors, record indices, exponent/coordinate data and generator buffers are explicitly bounded by the earlier storage dimensions.

### Success is a union bound, not a confidence or independence product

The fixed class/unit/power-basis conditions remain conditional algebraic premises. The stochastic production moment and transport hypotheses above have precisely the stated scopes. Main cap failures, including the revised2^-64 cold-row overflow, total at most10^-9; main matrix rank/formal-Krylov irregularity at most.01; coordinate/projection/sampler errors at most2^-66. The separate descendant-map counter abort is at most2^-45. H5/H5n/H5p refer to the actual distribution after all regular main gates, with R>=Rneed and true advice. Their conditional descent failure/work bounds control P(descent failure AND passed regular main); no independence of main and descent is asserted.

The global lower success is therefore

1-.01-10^-9-2^-66-[2*.001755+2^-10.61+2^-8.4+2^-80+2^-13+2^-45+2/r+2^(113.329-114.085)]-2^-48.

The BOTH-log meter is deliberately raised to114.085 to retain success above.39 under this conservative union. The exact certificate and independent Wolfram query evaluate this expression. It does not charge a false structural premise by a Schwartz-Zippel number, nor interpret two SE or synthetic confidence as attack success. The complete ledger includes collection116.925, SM103 and administration74 as actual new costs.

### Lemma: the actual universe has at least Ulo columns

Let pi_f(B),pi_g(B) count all prime ideals of norm<=B=2^69. The proved average count error gives (pi_f(B)+pi_g(B))/2>li_3(B)-15*10^15. Since1/logt>=1/logB on[3,B] andlog2<7/10, li_3(B)>(B-3)/(69*(7/10)).

Count only the following subset. In a degree-n field, ideals of absolute degree>=2 lie above rational primes<=sqrtB; there are at mostn above each prime. Their average over degrees72 and12 is at most42sqrtB; charge the larger72*2^35. Next omit all ideals above primes dividing the monic absolute defining-polynomial discriminant. The logarithm of each polynomial discriminant is<10000, so each field loses at mostn*10000/log2 such ideals. Even replacing both degrees by72 bounds their average by72*10000/log2<2^20; the last inequality follows exactly fromlog2>56/81 and72*10000*81<2^20*56.

At every remaining prime the defining power basis is maximal locally, and absolute-degree-one ideals correspond to simple F_p roots. Their contractions to K_h also have residue field F_p, so they belong to the required relative-degree-one universe. This good subset is invariant under sigma:x->-x. Its roots form free pairs:2 is a discriminant prime of these even polynomials; the f polynomial has constant1; and a prime dividing the g polynomial's constant dividesU_BLS, where its reduction isx^12 and is not squarefree. Thus no remaining root is zero and the odd-characteristic rootsx and-x are distinct. Half the number of these good ideals is therefore a genuine subset of the actual nonfixed sigma-orbit column universe U; fixed ideals are not being counted as columns.

Consequently, retaining an additional conservative84 subtraction,

    U>(2^69-3)/(69*(7/10))-15*10^15-72*2^35-2^20-84
     =5895711908692288528742/483.

Hence U>=Ulo=12206442875139313724>10^19. The checker verifies the floor and every rational substitution. This lemma is separate from Ucap; only Ulo can be used in the denominator of the concentration bounds. It is conditional on the already stated GRH count derivation, with no additional stochastic hypothesis.

Reproduction locations: losslessly decode stats16.b85 for historical empirical inputs, source_archive(source,SHA256) and the600 preregistered matrix records/frozen source. Decode filt16.b85/wolfram_main_phase_records for every new exact query and its untouched Wolfram output, including recoverp-main-final.wl/.txt. The live checker derives main caps/projection and all62node+boot recurrences; it does not substitute recorded descent outputs.

At the boot interface the polynomial element isalpha=a-b*x, while a finite-field plus-coordinate pair(A,B) representsA+B*x. The latticea*z1+b*z0=0 givesalpha=-(b/z1)z forz=T^e, so pi(alpha)=z; anchor1+x uses(a,b)=(1,-1). Use multiplicative boot blindingT^e and divide the normalized log bye. The discussion of additiveomega^b offsets only states the normalization identity, not a second algorithmic step.

Exact signatures use integer-polynomial Sturm sequences at infinity and zero: h has one positive and one negative root, while P(t)=t^6-2t^5+2t^3+t+1 has no real roots. Since f=y^6P(x^2/y)-h and fplus=y^6P(s/y)-h, both have no real embedding. The unique negative h-root supplies two real g-roots. Thus the fixed signatures and odd unit ranks are exact; no approximate root isolation is used. Protect theI orientation explicitly, even if a default orbit key would choose sigmaI.
