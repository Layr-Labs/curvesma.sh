# Probabilistic sieve descent, with independently randomized boot

finalmodel4.py contains the full 82-row plan and mesh96 evaluator. It replaces the fixed-auxiliary-family algorithm of BN3. There is NO assertion that every individual ideal has a nonempty fixed family. The target, field, full parameter sizes and C=20,000 remain unchanged.

## Nodes and qualification

Use the original rank12 special-q lattice L_q, reduced under QN, determinant q. Enumerate the row's radius-R coefficient ball uniformly modulo sign; full enumeration is deterministic. The cheaper of sieving and finite ECM on every cell is chosen by node_cost. Sieve rows use all degree-one root ideals in [2^8,2^56], certified norm upper bounds (2400 instructions/cell), the original threshold T=kA*nA+kB*n+96, exact surviving norms, small-prime division and finite cofactorization. A relation obeys both two-window prime bounds, primitivity, ordinary-ideal, unit-duplicate and known-log rules. Missing small logs are rejected on discard rows or handled by the original <=77-bit auxiliary helper rows; those helpers discard missing logs, preventing a helper cycle. All other prime windows strictly decrease the input bin.

Every ball's complete sieve/bound pass and ALL modeled cofactor survivors are charged. There is no lazy-work discount. The old DISC certificate-discard step is removed: it is not used to turn pooled zero-failure samples into a uniform per-type guarantee. Geometric certificates remain scoped diagnostics, not actual relation-count certificates.

## H-PGF: reached-ideal ensemble envelope

Let M be the number of distinct qualifying relations before child descents. The raw scale mu is node_cost's Dickman/Murphy, retention, detection and norm-envelope prediction. It is NOT asserted to equal the measured mean of M. For t in [0,0.8], H-PGF uses

    E[t^M] <= (1 + 0.5*mu*(1-t)/r)^(-r),

with r=0.5 for effective windows of at most one prime, r=1 otherwise. The factor0.5 is a calibrated conservative envelope parameter. No exact negative-binomial distribution or validated production mean is claimed. The law concerns the reached-ideal ensembles used by the recurrence, including parent selection and earlier relation/child failures, before individual labels are conditioned away. It is not a random-count law conditional on a fully fixed q.

For one relation, psi is the union bound over its child failure envelopes. H-PGF and the conditional next-relation failure envelope after prior failures give phi=max(PGF(psi),0.02), at most0.04 in this plan, and expected attempted relations at most1/(1-psi). This conditional transfer is an explicit empirical hypothesis; no multilevel independence theorem is asserted. Work is own_full_ball + child_work/(1-psi). finalmodel4 evaluates this increasing-bin recurrence. It also recomputes expected node invocations. The run aborts at the declared node budget or a repeated ideal; the original size/weight repeat bound and Markov node-budget bound are both charged. The global work meter additionally covers every failed attempt and abort.

## Relevant evidence and limits

The BLS-style port tests actual complete radius2 balls over the BN base ring, using toy f=x^4+w*x^2+w^2 and g=x^2+4w. Finite windows are fully factored and root-qualified; these are not unrestricted-prime-residual proxies. All192 fixed-ball runs and all335 encountered-child/control pairs are retained. Encountered/control means were8.510/8.740, zero counts7/12;140 parent relations gave failure0.050 and sibling-conditioned0.0522. The parent-cluster difference interval was[-0.0483,0.0360], descriptive at this toy scope.

Raw-mean NB lazy pricing FAILED some tests (up to1.508 times its bound); we therefore charge full cofactor work. Raw PGF also failed a small-mean toy. With the fixed0.5 envelope parameter, all862 count observations pass PGF comparisons on t=0..0.8 in steps0.01, worst empirical/model ratio0.824. The child/control functional check is a transfer check, not a blind holdout. pgf_data.b85 retains all counts and raw scales rounded DOWN to0.001, conservatively increasing the model PGF. No observations are clipped. Replay part=pgf recomputes the table; part=porttoy regenerates full toy/parent/child data from seeds.

The production encountered-child campaign adds12 parents,3072 parent cells and48 matched children/controls. Mean conditioned offset shifts were cf+1.212,cg+1.171bits; raw maxima73.73/6.82 were below charged98.3/14.6. Individual paired mean shifts exceeded the old nominal10/4 margins, so those are NOT universal bounds. Original production norm samples and per-ideal predicted-yield ratios support the transfer from toy per-ideal scales to production bin scales. Neither geometry nor toy data proves that transfer. H-PGF extrapolates across polynomial degree, prime size and adaptive depth; it includes production norm heterogeneity and survivor-law uncertainty. The argument is empirical exploratory support, not proof from a positive predicted mean.

## ECM hit-to-extraction bridge

All schedules are finite. Keep the original composite modulus throughout each curve, collect stage1 prime-power-prefix and stage2 giant-block product gcds, and refine the divisor partition. No gcd on each baby/giant term is priced. Charge unconditional replay plus five full-modulus gcd budgets per prefix/block output for gcd, remainder and refinement. extraction_price charges K*(2*ops*mul+2*gcd+5*(pi_upper(B1)+G)*gcd).

The standard random-order model supplies corrected hit probability P per curve. H-SEP separately lower-bounds complete target-factor isolation conditional on its hit by1/2; it is not inferred from P alone. Fresh curve coins then give extraction probability at least P/2. K=ceil(log(0.001)/log(1-P/2)) bounds each required prime's miss by0.001, followed by the explicit factor-count union bound. This amplification is included in the price, unlike BN3's unsupported complete-partition delta.

stage2-summary.json records480 curves,36/40-bit target factors and original moduli200/800/1746bits:133/135 hit factors isolated by same-curve prefix/block replay;35/36 small factors separated overall, one had no hit. Timeouts and bounded reruns are disclosed. This tests real stage2 extraction, not merely order hits. H-SEP still extrapolates to production factors and multiplicities; the measured rate is not a uniform theorem. Scalar loop/output accounting and the per-node construction allowance include control work.

## Fresh-exponent boot and memory

Each trial chooses a fresh uniform nonzero e, forms T^e, builds/certifies the boot lattice under QB and samples ONE radius4 cell. This is an iid single-cell target population for nonidentity T. Use (n,iA,kA)=(220,48,16), nA=138, plus the original shifted norm-bin model, counted-cell/duplicate factors and an extra half-rate allowance for boot gates. Factorization uses the amplified finite schedules above. After a relation, try its children; if they fail, start a fresh exponent. Expected boot work is (trial_cost/p_trial+child_work)/(1-psi). No BALL_PAIR, whole-ball waiting assumption or universal node positivity is used. The global meter and node/repeat events give the declared success bound.

The model charges all rows, including unreachable high-bin table entries; the reachable child primes are below2^220. Conservative peak storage is below2^91bytes (full-plan sieve bound2^89.281 plus main storage and scratch). Memory is separate and unscored. No work is amortized over another target.

The fixed node budget is2^22: expected calls2^12.993756, node-abort bound2^-9.00624 and repeat bound2^-8.60928. These losses are included.

Survivors follow deterministic cell traversal: factor and descend one at a time, stopping at first success. Full survivor work is still charged. H5 applies to this order; complete factorizations are not retained together.

Primality uses128 fresh MR bases; total error<2^-125 by2^131 calls. Checks fit node TDIV2^32 and boot2^40 budgets.
