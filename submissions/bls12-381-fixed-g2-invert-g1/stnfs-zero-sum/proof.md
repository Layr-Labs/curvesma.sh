# Recover-P 117.618259017806: f-only collection, certified counts and audited norm transport

Track `bls12-381-fixed-g2-invert-g1`. The first manifest claim gives the exact rational score. The complete argument is this file together with DESCENT.md, KERNELS.md, STORAGE.md and field-map.md. The bound is analytic and conditional on the disclosed heuristics, including GRH, H5n and H5p.

## Fixed supplied-Q reduction

Validate the supplied nonidentity Q in G2 and the GT target z. For z=1 return infinity. Otherwise compute g=e(P0,Q), using the pinned `pairing(Q,P0)`. Map both g,z to the STNFS tower using the cold isomorphism construction in field-map.md. One relation/LA phase supports TWO individual logarithms, of g and z. Recover a=log_g(z) modulo r and return [a]P0 only if e([a]P0,Q)=z. Nondegeneracy makes g a generator. Uniform Q gives uniform nonidentity g independently of uniform z, exactly the two-log input law. No replacement Q or free advice is used.

The wrapper has a fixed2^80-instruction cap, including all random field-factorization attempts, pairings, mapping, scalar multiplication, output and checks. At C=20000 it fits2^66 units; failure allowance2^-48 (field-map.md). Every failed check or exhausted cap is a failure.

## Inherited STNFS collection and arithmetic

Use h=y^6+y^3+y^2-y-1, f=x^12-2yx^10+2y^3x^6+y^5x^2+1+y-y^2-y^3, g=x^2+U y, U=15132376222941642752; x->-x is the odd sigma quotient. config.json fixes the exact rational region: A²=2^21, pooltheta3/8, regiontheta2719/8000, lam17/20, kap1. Exact canonical unit/HNF and primitive tests are implemented in ideal_canonical.py/runs16.py. lattice_count_bounds.py supplies rational yield/work geometry factors0.970276956/1.028809574.

Only f is sieved; the exact SRT g batch decides g smoothness. model16.py, KERNELS.md and STORAGE.md specify the line/record/run walks, certified thresholds and warm starts, contamination fallback, routing/recovery, division, bandstrip/split, cofactor-window ECM and g splitting. Every setup, tree, failed trial, transfer and fallback is charged there; DESCENT.md gives both logarithms. S/G eachhave656000draws; their recorded work sums/squares and yield sums in pool_d.json are empirical inputs, with selected anchor replay in evidence16.py. The fast certificate does not regenerate that historical campaign.

Credit gbotrel's transferred f-only scan/bound deletion: g-array operations are removed with unchanged f predicate. The shipped raw kernel log yields exact slope260281025/189262038<=1.376, leaf excess148329356571/3942959125<=37.62, bound cost3416387/2959<=1154.6 and residual2271913/23788165<=.125<.5. The certificate decodes/rechecks every point. S/G work SE=.004213089760336546 and1.00843>=1+2SE give collection<=2^116.924540000. G yield after.009-bit bandloss is2^63.62660550353605; its2-SE lower2^63.576124844219954 supports descent input2^63.575. The fixed inherited region remains feasible with the corrected U; the historical S minimum-grid rule is not reasserted. H1/H2/H4/H-FFT assumptions and their measured scope remain explicit in the associated files.

## Direct prime-ideal counts (reproduced proof)

The truncated li expansion is a heuristic density only. Dimensions and repeat bounds use the direct derivation reproduced at the end of this file, conditional on GRH for the two fixed Dedekind zeta functions. It derives, for all integer69<=b<=225,

|pi_K(2^b)-integral(3,2^b,dt/logt)| < E_b=2^(ceil(3b/4)+3).

In particular E_69=2^55. No quantitative Grenie-Molteni estimate or external zero-count/Laurent-constant bound is used. The standard analytic continuation, functional equation and Hadamard product are stated with their hypotheses; the zero-weight estimate, every explicit-formula residue, unsmoothing and partial-summation constants are derived below. Numerical checks support the reproduced argument rather than purporting to prove an assumed theorem.

The absolute resultant polynomials have degrees72/12 and are irreducible over Q (Wolfram code/raw in stats16.b85). Since |y|<=2, at |x|=3 the f leading magnitude531441 exceeds the other upper magnitudes248163, so all f roots satisfy|x|<3. Thus logDelta_f<=72*71*log6<10000. For g, |x|<2^33 because U<2^64, giving logDelta_g<=4488log2<10000. Field discriminants divide their monic polynomial discriminants.

For F(x)=x/(logx-1.1), F'(x)>=1/logx when logx>=12.1. Splitting at2^18 gives integral(3,B)<=F(B)+2^18. A sigma-fixed ideal has zero odd log; other orbits have size two. The direct bound counts ALL absolute prime ideals and upper-bounds the relative-degree-one subset. The field-specific endpoint refinement proved below bounds the average of f/g count errors by15000000000000000. This is sufficient for the half-sum; individual repeat counts still use E_b. With redundant boundary84,

U_actual<=ceil[B/(69log2-11/10)+2^18+15000000000000000+84]=12647821419683016456.

counts16.py uses100positive atanh terms as a rational lower log2 bound. It checks the direct-proof constants, root/discriminant bounds and exact ceiling. Ncap=.03330Ucap<2^58.55<6*2^56, so six owners suffice. Filtering H3/H3c refers to actual U, with abort allowance<=1e-9. The larger cap gives LA<=2^115.822494853 and generator<=2^110.610749051.

Independently on EACH side in every actual bin[2^a,2^b], subtract E_a+E_b, at most72*2^ceil(b/2) higher-degree ideals, and2^20 ideals above bad polynomial primes. Simple nonzero roots form free sigma pairs. Each orbit pool is at least

[(2^b-2^a)/(b log2)-E_a-E_b-72*2^ceil(b/2)-2^20]/2 >2^68.

All31exact rational checks pass; minimum floor357084587023637308834. Under the explicit H5 arrival-weight allowance2^6, Nmax=2^27.3 gives repeat tax<=2^(2*27.3-1+6-68)=2^-8.4. No abort is added by discarding the exceptional primes in this lower count.

## Runtime norm-CDF audit and empirical transport (new)

The historical norm bins are proposed thresholds, not universal envelopes. Before using a node basis, audit its finite integer ball with fresh, exact uniform samples. For the24node thresholds(a_i,b_i), replace b_i by min_{j>=i}b_j. These nested rectangles E_i have monotone thresholds. Require true mass P(E_i)>=eta*i/24 with eta=19/20; boot uses its10scalar thresholds, checked by exact whole-ball enumeration for each fresh basis. A layer-cake sum proves that these CDF bounds imply at least eta times the original sub-bin semismooth lower sum for every coordinatewise decreasing smoothness function. Since the old charge was the minimum of that sum and a capped sextile estimate, multiplying the old charge by eta remains conservative. No sampled maximum is treated as a universal bound.

At ordinary nodes use2^46samples, epsilon=2^-15, accepting only if every empiricalCDF>=eta*i/K+epsilon. Hoeffding gives total false acceptance<=48exp(-2^17), and false rejection on a basis whose trueCDF>=eta*i/K+2epsilon is equally negligible. A uniform integer-box proposal in[-R,R]^12, rejecting points outside the ball and zero and folding signs, is exactly uniform on the half ball. Its acceptance probability is>=2^-16 for all shipped radii (checked from integer CELLS). Cap proposals at2^64 per audit; the binomial lower-tail abort and capped uniform-integer generation are charged. Exact norm-threshold comparisons use precomputed integer100th-root thresholds; no unpriced arbitrary-width unit operation is used. Total audit work, including rejected bases, setup and randomness, fits a separate2^105-unit row (detailed in DESCENT.md).

A bad node basis is rejected and fails. H5n is the explicit quantitative transport heuristic: under the fresh special-q distribution actually reached in each production bin, the mass of bases violating the strongerCDF condition is<=delta=1/200. This is NOT deduced as a theorem from pooled samples. A fresh held-out replication uses2000production-sized special-q lattices covering62side/size strata (each sampled in the top1bit of its5/6bit bin) and32fresh uniform-GT boot inputs; recorded minimum-constraint witnesses and exact statistical summaries are in stats16.b85. All2000node and32boot tests pass at eta=.95,2048norms per lattice, conservative eps=.041/.036. The certificate recomputes these empirical summaries and the pooled binomial/Chernoff evidence, not every lattice or the historical full campaign.

Scope: witness primes/root choices and top1bit samples differ from the production arrival law. H5n extrapolates across root multiplicity, weighting, within-bin shapes and exact-norm to certified-bound classification. Offline Gaussian/rounding sampling and floating logs are approximations, not exact certificates. The positive margins and coverage are evidence for this testable assumption, not a proof or an algorithmic95%success probability. Full counts/reproduction code are retained in the public note. The runtime audit removes reliance on untested CDF envelopes for accepted bases; H5n supplies its rejection-rate tax. A62-lattice/992point calibration gives1984positive gaps, zeroCDF/cap disagreements and maxgap5.161e-10bits, at that scope only.

The conditional negative-binomial/independent-subtree model H5 remains, supported at its measured scope by the existing toy fits and exact unit-orbit checks (DESCENT.md). Its0.02floor is additional padding, not the norm-tail allowance. The corrected recurrence is phi=delta+(1-delta)*max(phi_NB,0.02). Every bad-lattice failure is therefore ADDED. Work uses the corresponding retry factors; charging a full node even when its audit fails is conservative. Boot resamples fresh exponents after failedCDF audits, capped at128 attempts. H5n also assumes boot bad mass<=1/2;32positive input tests support this. Fewer than2^40good balls are attempted under the charged meter, so the union of128-attempt audit exhaustion events is<=2^-88. All attempts are charged. The model uses safe LOG2U6=382.286, above the exact382.28563411957; LOG2P6 remains the existing safe upper2285.982.

## Meter, success and complete ledger

The recurrence, node count and boot are recomputed by audit_refinement.py under the20-second hosted experiment limit. Ordered zero-only convolution trimming preserves terms. Moving upper omitted-mass tails above512bits to infinite mass is conservative: all denominator cutoffs are below512bits, all increments are nonnegative, and the loss numerator is bounded by1. Lower omitted laws remain uncapped. No recorded descent-output cache substitutes for the recurrence.

Expected two-log work<=2^113.329, boot phi<=.001755, expected node log2<=16.682, node-abort charge2^-10.61, repeat charge2^-8.4, new audit-error charge2^-80 and prefix-abort charge2^-13 are in refinement.json. The meter for BOTH logarithms is2^114.075 units. Markov charges E[W]/meter; node overflow is charged from the recomputed expected count; repeats use2^-8.4; all audit, boot and wrapper failures are charged. Core success is the product of main-phase success and one minus these failure charges. The older+10%work-sensitivity claim is removed. The final certificate requires success>.39, not reviewer confidence.

Ledger rows, all in dense-Fp12 units at C=20000: collection116.924540000; LA115.822494853; generator110.610749051; both-log meter114.075; canonical list91; Schirokauer maps88; filtering/misc80; factor base73; wide lattices70.6; rec7bases73.5; certificate construction50; verification11; g-tree74; g-split95; band-split101.54; newCDFaudits105. Sum all rows across all processors once. The wrapper2^66 is additional. The exact first claim is an outward-rounded logarithm covering this sum; the certificate checks it with80-digit arithmetic and Wolfram verifies the rational inputs independently.

Rows include all cold advice, failures, arithmetic, memory/lookup/transfer, processors, randomness, filtering, output and amplification; no amortization or free oracle. Broad fixed factor-base/root construction remains2^73 by the inherited aggregate per-rational-prime price and unconditional pi_Q(B)<1.25506B/logB. Memory remains below2^106bytes: the2^105 core plus wrapper and streamedCDFcounts/thresholds. Samples are streamed, not stored en masse.

## Remaining hypotheses and provenance

H1/H1' are Dickman/Murphy/retention and smoothness-local-law assumptions; H-EA is the survivor-law ECM share; H2 is the uniform-pool work/yield model; H3/H3c are filtering, matrix rank, unit conditions and concentration; H4/H4u are kernel/shape/root-density extrapolations and entry-cost/unit shares; H-X is Block-Wiedemann projection nondegeneracy; H-FFT is the distributed multiplication extrapolation; H5 is the descendant independence, negative-binomial, missing-log and unit-orbit model; H5p is the marked common-intensity extension below. New GRH and H5n are explicit above. Small experiments only support their actual scale; all extrapolations remain material.

Credit gbotrel explicitly: recover-P base e1e57459-af88-4c69-b5ec-4e86829fd81d and sibling f-only lever ff9f9487-e306-40e3-b6cb-feb0a199515f (after9524d511). Earlier provenance remains the credited lineage in those source packages. The present new folder corrects4d15ccc5's count/tail objections; no prior verdict is a premise. Protected files and old folders remain untouched.

The core meter117.618259017805 plus wrapper2^66 fits the exact declared score117.618259017806. The certificate recomputes every ledger row and the conservative success bound0.396>.39.

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

    U_cap=ceil[B/(69 log2-11/10)+2^18+15000000000000000+84]
         =12647821419683016456,

with the same exact atanh lower interval for log2. The84 allowance may be retained conservatively. For repeat probabilities on an individual field, continue to use the generic proved endpoint radius E_b=2^(ceil(3b/4)+3), then account for non-degree-one ideals, bad-prime subsets and free sigma orbits as before.

Every numerical step in this refinement is independently checked by `field69_bound.py` and the recorded Wolfram query/output `wolfram-field69.wl`/`wolfram-field69.txt`.

