# Probabilistic sieve descent, with independently randomized boot

finalmodel4.py loads the full 82-row plan from model4core.py in lib.b85 and mesh96 evaluator. There is NO assertion that every individual ideal has a nonempty fixed family. The target, field, full parameter sizes and C=20,000 remain unchanged.

## Nodes and qualification

Use the original rank12 special-q lattice L_q, reduced under QN, determinant q. Enumerate the row's radius-R coefficient ball uniformly modulo sign; full enumeration is deterministic. The cheaper of sieving and finite ECM on every cell is chosen by node_cost. Sieve rows use all degree-one root ideals in [2^8,2^56], certified norm upper bounds (2400 instructions/cell), the original threshold T=kA*nA+kB*n+96, exact surviving norms, small-prime division and finite cofactorization. A relation obeys both two-window prime bounds, primitivity, ordinary-ideal, unit-duplicate and known-log rules. Missing small logs are rejected on discard rows or handled by the original <=77-bit auxiliary helper rows; those helpers discard missing logs, preventing a helper cycle. All other prime windows strictly decrease the input bin.

Every ball's complete sieve/bound pass and ALL modeled cofactor survivors are charged. There is no lazy-work discount. The old DISC certificate-discard step is removed: it is not used to turn pooled zero-failure samples into a uniform per-type guarantee. Geometric certificates remain scoped diagnostics, not actual relation-count certificates.

## H-PGF: reached-ideal ensemble envelope

Let M be the number of distinct qualifying relations before child descents. The raw scale mu is node_cost's Dickman/Murphy, retention, detection and norm-envelope prediction. It is NOT asserted to equal the measured mean of M. For t in [0,0.8], H-PGF uses

    E[t^M] <= (1 + 0.5*mu*(1-t)/r)^(-r),

with r=0.5 for effective windows of at most one prime, r=1 otherwise. The factor0.5 is a calibrated conservative envelope parameter. No exact negative-binomial distribution or validated production mean is claimed. The law concerns the reached-ideal ensembles used by the recurrence, including parent selection and earlier relation/child failures, before individual labels are conditioned away. It is not a random-count law conditional on a fully fixed q.

For each relation psi is the sum of its child failure envelopes. The sequential hypothesis is P(F_i|M>=i,F_1,...,F_(i-1),pre-entry history)<=psi; the per-attempt child-work envelope feeds the expected descent work, which the claim uses only through the Markov cap loss (proof.md section 7). THRESHOLD-REPAIR.md proves phi<=P(M<k)+psi^k, with the low-count tail bounded by G(t)/t^(k-1), and uses a fixed finite certificate grid k1..64,t<=.8. It does not substitute PGF(psi) or assume count/solvability independence. Expected attempted relations remain at most1/(1-psi), giving own_full_ball+child_work/(1-psi). The fixed82-row plan increases only radii needed to restore phi<=.04 under this corrected bound. Expected invocations, node/repeat losses, phase meter and total work are recomputed. The run aborts on the declared node budget or repeated ideal.


## Relevant evidence and limits

The BLS-style port tests actual complete radius2 balls over the BN base ring, using toy f=x^4+w*x^2+w^2 and g=x^2+4w. Finite windows are fully factored and root-qualified; these are not unrestricted-prime-residual proxies. All192 fixed-ball runs and all335 encountered-child/control pairs are retained. Encountered/control means were8.510/8.740, zero counts7/12;140 parent relations gave failure0.050 and sibling-conditioned0.0522. The parent-cluster difference interval was[-0.0483,0.0360], descriptive at this toy scope.

Raw-mean NB lazy pricing FAILED some tests (up to1.508 times its bound); we therefore charge full cofactor work. Raw PGF also failed a small-mean toy. With the fixed0.5 envelope parameter, all862 count observations pass PGF comparisons on t=0..0.8 in steps0.01, worst empirical/model ratio0.824. The child/control functional check is a transfer check, not a blind holdout. pgf_data.b85 retains all counts and raw scales rounded DOWN to0.001, conservatively increasing the model PGF. No observations are clipped. Replay part=pgf recomputes the table; part=porttoy regenerates full toy/parent/child data from seeds.

The production encountered-child campaign adds12 parents,3072 parent cells and48 matched children/controls. Mean conditioned offset shifts were cf+1.212,cg+1.171bits; raw maxima73.73/6.82 were below charged98.3/14.6. Individual paired mean shifts exceeded the old nominal10/4 margins, so those are NOT universal bounds. Original production norm samples and per-ideal predicted-yield ratios support the transfer from toy per-ideal scales to production bin scales. Neither geometry nor toy data proves that transfer. H-PGF extrapolates across polynomial degree, prime size and adaptive depth; it includes production norm heterogeneity and survivor-law uncertainty. The argument is empirical exploratory support, not proof from a positive predicted mean.

## ECM hit-to-extraction bridge

All schedules are finite. Keep the original composite modulus throughout each curve. Store stage1 prime-power-prefix Z coordinates and stage2 giant-block products on the first pass. The stored-prefix lemma below replaces unconditional arithmetic replay and charges only the necessary endpoint-search prefix gcds plus every stage2 block; it preserves the complete factor partition. No gcd on each baby/giant term is priced.

The standard random-order model supplies corrected hit probability P per curve. H-SEP separately lower-bounds complete target-factor isolation conditional on its hit by1/2; it is not inferred from P alone. Fresh curve coins then give extraction probability at least P/2. K=ceil(log(0.001)/log(1-P/2)) bounds each required prime's miss by0.001, followed by the explicit factor-count union bound. This amplification is included in the price, unlike BN3's unsupported complete-partition delta.

stage2-summary.json records480 curves,36/40-bit target factors and original moduli200/800/1746bits:133/135 hit factors isolated by same-curve prefix/block replay;35/36 small factors separated overall, one had no hit. Timeouts and bounded reruns are disclosed. This tests real stage2 extraction, not merely order hits. H-SEP still extrapolates to production factors and multiplicities; the measured rate is not a uniform theorem. Scalar loop/output accounting and the per-node construction allowance include control work.

## Fresh-exponent boot and memory

Each trial chooses a fresh uniform nonzero e, forms T^e, builds/certifies the boot lattice under QB and samples ONE radius4 cell. This is an iid single-cell target population for nonidentity T. Use (n,iA,kA)=(192,71,12), nA=156+7/12, plus the original shifted norm-bin model, counted-cell/duplicate factors and an extra half-rate allowance for boot gates. Factorization uses the amplified finite schedules above. After a relation, try its children; if they fail, start a fresh exponent. Expected boot work is (trial_cost/p_trial+child_work)/(1-psi). No BALL_PAIR, whole-ball waiting assumption or universal node positivity is used. The logarithm-phase meter and node/repeat events, together with split-only preprocessing metering, give the declared success bound; no numeric whole-attack work cap is asserted.

The model charges all rows, including unreachable high-bin table entries; the reachable child primes are below2^192. Conservative peak storage is below2^100bytes (full-plan sieve bound2^92.436788 plus main storage and scratch). Memory is separate and unscored. No work is amortized over another target.

The fixed node budget is2^21: expected calls at most2^12.687665913, node-abort bound2^-8.312334087 and repeat bound2^-10.609281891. These losses are included.

Survivors follow deterministic cell traversal: factor and descend one at a time, stopping at first success. Full survivor work is still charged. H5 applies to this order; complete factorizations are not retained together.

Primality uses128 fresh MR bases; total error<2^-125 by2^131 calls. Checks fit node TDIV2^32 and boot2^40 budgets.

The fixed plan includes threshold-feasibility radius repairs and retains all82 production node types, alpha=0.5, H-SEP=0.5, full cofactor work and the stated qualification rules. The complete exact plan is the PLAN constant in model4core.py inside lib.b85, loaded by finalmodel4.py. Curve p, order r, input distribution and the0.39 success threshold are unchanged.

# Stored-prefix ECM extraction

Before each norm's first ECM schedule, completely divide all primes <=2^20. This exact step is covered by the existing node TDIV=2^32 instruction row; boot construction's2^40 allowance also contains it. Nested cofactor schedules inherit the absence of small primes. Keep the resulting composite modulus N fixed throughout every curve, as before. If N=1 or prime, skip needless curves; pricing them anyway is conservative. Let b=ceil(log2 N), L=max(5,ceil(b/64)), Q an upper bound on pi(B1), and G the charged number of stage2 giant-block products.

Store the N-bit Z coordinate after every stage1 prime-power multiplication and each stage2 giant-block product during the FIRST curve pass. No curve-arithmetic replay is needed. These are precisely the values used by the old unconditional replay, with the same sigma, modulus and arithmetic. Storage is contiguous limb arrays, not per-object allocations.

For every divisor d|N, a Montgomery doubling preserves d|Z: its new Z is4XZ*((X-Z)^2+4a24XZ). For a differential addition, if both operand Z coordinates are divisible by d, then u-v=2(zy-xw) is divisible by d, and the new Z=diffX*(u-v)^2 is divisible by d. A ladder starts with P and2P, so induction gives d|Z(P) => d|Z([k]P). Reduction modulo N preserves this assertion. No field, squarefreeness, nonsingularity or hit-distribution assumption is used.

Consequently d_j=gcd(Z_j,N) form a divisibility chain. Since N has no prime <=2^20, every strict chain increment increases log2(d_j) by more than20. There are at most F=ceil(b/20) such increments, including prime-power valuation increments.

Compute the endpoint gcds of a prefix-array interval. If they agree, every interior gcd agrees, because it lies between the endpoints in the divisibility order. Otherwise bisect, evaluate the midpoint gcd, and recurse, carrying endpoint values down the recursion. This recovers every distinct prefix gcd. At most min(Q,2F*ceil(log2 Q)+2) gcd queries are needed: at each depth every unequal interval contains a strict increment, and depth is at most ceil(log2 Q). Endpoint queries and integer-rounding slack are included. Distinct divisors suffice: reapplying an already-refined divisor cannot further change the integer factor partition.

Scan ALL G stage2 block products with the original gcd/refinement procedure. The complete set of distinct stage1 divisors and all stage2 divisors is unchanged. The same common-refinement partition therefore results, including prime powers. H-SEP and the sigma law are unchanged; no stronger conditional extraction probability is claimed. The work bound holds for unsuccessful and nonsmooth norms too, without a hit-conditioned amortization premise.

Per curve charge

    ops*mul + 8L(Q+G) + 2g + 5*(min(Q,2F*ceil(log2 Q)+2)+G)*g.

Here ops already contains the original1.1 arithmetic margin; mul and g are the original ecm16 multiplication/gcd prices. Eight scalar instructions per stored limb budget initial writes, later operand copies/reads, endpoint handling and address/control work; each stored prefix is queried at most once. The original five full-modulus gcd/remainder/refinement budgets remain per queried output. Curve setup/final tests retain2g. Curves remain ceil(log(0.001)/log(1-0.5*Pcurve)). The same original finite B1/B2/wheel grid is minimized under this revised deterministic price.

Stored output memory is8L(Q+G) bytes, plus the baby points, O((F+log Q)L) endpoint/divisor limbs and existing curve buffers. The combined fixed-plan replay observed maximum stored-output allocation2^29.57454 bytes; production node arrays dominate the run's stated memory.

Exact structural verification: six original-modulus curves at200/800/1746-bit N containing1009^2*1013*65537^3 and odd composite cofactors, stage1 B1=1024 and stage2 B2=51200, wheel210. All172 prefix gcds were checked, all244 block outputs retained, and old/new refined partitions compared exactly. Four to six distinct prefix gcds were recovered in14–17 queries. All six partitions agreed, including multiplicities. This is an exact scaled algorithm-equivalence check, not an extraction-probability calibration. An initial fixture with a cofactor divisible by3 yielded no valid Suyama setup; that failed artifact is retained, and the revised fixture chooses cofactors coprime to primes<=97.

The fixed plan has61 sieve rows, raw mu>=374.960544394, phi<=0.039483332, psi<=0.550175912, retention>=0.828360856 and boot retry2.367893567. Construction, failures and full cofactor work remain charged. Peak arrays are below2^92.436788 bytes each; at most84 simultaneous frames plus main/scratch fit2^100bytes (STORAGE.md). Memory is unscored.

## ECM-cell qualification and bounded preprocessing

The full-cell ECM branch applies the identical exact base-primitivity and ordinary-degree-one-ideal predicates used by the sieve branch. Its raw mu is therefore multiplied by retain16.counted_fraction(), the union bound for the two ordinary predicates times base_primitive(). It does not pay the sieve threshold-loss retention factor, because no approximate score threshold selects these cells. The fixed (242,g) row uses radius12 after this correction; all82 types still satisfy the checked PGF envelopes. That row lies above the boot192 reachable range.

For every full-cell ECM trial, construct both exact norms, remove the special-q factor exactly, and run the finite <=2^56 stripping schedule, including small-prime trial division. For each resulting residual z, compute its exact bitlength and reject the cell unless z.bit_length()<=floor(T). This test is performed BEFORE invoking the nA-bit extraction schedule priced at modulus bound T. The nA schedule is never run on a larger residual. A qualifying two-window relation has residual product at most2^(kA*nA+kB*n), so T adds the existing96-bit slack and this exact gate removes no qualifying relation in the modeled population. Both sides must pass. Fully factor and root-qualify accepted residues, verify all exponent counts and windows, then apply the ordinary/primitivity, unit-duplicate, known-log and auxiliary-log predicates exactly. Charge all cells and both stripping prices even if an earlier rejection would save work.

Perfect powers are handled explicitly before asking for prime isolation: after exact trial division through2^20, for every prime exponent k<=ceil(input_bitlength/20), compute an integer kth root by bounded binary search and verify its kth power exactly; recurse on a successful power decomposition and retain multiplicity. The tested integer and intermediates are capped by the original norm bound (less than2000 bits in this plan); multiplication compares against the bound and terminates once exceeded. A prime-power block is represented by its base and exponent, so extraction need not ever produce gcd=p rather than p^e. The exact polynomial-time root checks and reconstruction are allocated at most2^31 scalar instructions within node TDIV=2^32 and boot construction2^40; trial division, bounded residual comparisons and other TD checks use the remaining node allowance. The trial-division invariant implies every nonunit power base exceeds2^20, so its exponent is less than input_bitlength/20; no possible power is omitted. For a 2000-bit operand there are25 candidate prime exponents<=100. The sum of ceil(2000/k)*2ceil(log2 k) over those exponents is22058 power-multiply steps. At64 scalar instructions per32x32 limb pair this is at most1,445,593,088 instructions; comparisons, root-search control and reconstruction fit the remaining2^31 allocation. Smaller recursive blocks reduce the operand lengths and are covered by that slack. These checks are deterministic and bounded. The original H-SEP bridge still concerns the resulting prime/prime-power blocks.

Keep BOTH original final gcds as well as the stored output arrays: one gcd for the final stage1 Z coordinate and one for the TOTAL stage2 product, together charged by the unchanged +2g term. Scan every stored stage2 giant-block product and refine by the total-product divisor too. This preserves the old final-gcd partition, including valuation information when prime powers occur. The endpoint search replaces only repeated stage1 prefix queries; the deterministic price retains all block/refinement work and both final gcds.

See ADAPTIVE.md for28new preregistered multilevel recursive runs and all171entered-node records, including censors and pre-entry history strata. They supplement the earlier862-count calibration; they do not prove uniform production adaptive transfer.
