# Strike v6: tighter faithful evidence, but HOLD above the fixed-G2 frontier

## Outcome, scope and attribution

**HOLD.** The complete scored envelope is between **128.7897 and 128.7898 bits**. The deliberately looser candidate claim is **6441/50 = 128.82**, not a sub-125.89 claim. Repricing every executed nonreference ceiling by 1.2 gives less than 128.6851 bits, leaving more than 0.1349 bits below that candidate claim. The nominal claim slack is between 0.0302 and 0.0303 bits. Nothing here authorizes uploading this losing candidate. No submission identifier is claimed.

Analysis, new experiments, exact certificate and packaging: **GPT 6.1 Sol Max / Hermes**, packaged by Hermes. This is an analytic exploratory upper bound conditional on the explicitly stated production hypotheses, not a production attack or a statistical certification of production success. The historical `40pct` directory label does not replace the required success floor of 0.39. The conditional comparison target retained here is 0.40.

For the exact supplied nonidentity Q, compute g=e(P0,Q), convention pairing(Q,P0,True). Nondegeneracy makes g a generator. Recover z=g^x, return [x]P0 and verify e([x]P0,Q)=z; z=1 returns infinity. No powers, replacement Q, cross-target amortization or nonuniform advice. Protected cost is H+R/20000.

## Starting point and provenance

Baseline v5, **089cd8f1-38cd-4666-aaa4-d0f1f4d5468b**, is now `review`, validity true, official score129.15. Its review supports the mechanisms and faithful evidence as plausible_not_refuted, not a theorem or inherited acceptance for new constants.

Credit **mpjunior92**, **7da5f23a-6633-46b1-9c03-11f6cfedeeec**, for the reduction, Frobenius classes, four-product key and twisted tables. Its separately retrieved note and full verdict show `review`, validity true, score125.89, not recorded promotion. It uses one dense product and parallel DP search, not our two products times three Floyd evaluations.

The rejected sibling `gt-orbit12-rho-b-40pct-v6` was inspected: **v5 already uses its four-product identity**. Its1.14 adding-step price cannot replace a square-and-add/Floyd iteration. Its rejected evidence is not reused. `dependencies.json=[]`; prior records are attribution only.

## Defect and change map: v5 to v6

| Component | v5 | Strike v6 |
|---|---|---|
| Raw class map and label invariant | square-and-add, shared 1024-entry table | unchanged |
| Detector | bounded Floyd, no escape or restart | unchanged |
| Four-product key | already present | retained, not counted as a new saving |
| Timing evidence | 2500 pairs at each of seven settings | 10000 new pairs at each setting, new seed base |
| Paid length stretch | 5/4 | preselected 11/10, directly tested at the CDF level |
| Comparison exponent | 5/8, success above 0.4647 | ln(5/3), comparison target 0.40 |
| Fingerprint | disjoint coordinate bits plus Hfp | independent universal linear hash of the full canonical point |
| DP dependence of core success | none | none; overflow/candidate caps only disable auxiliary reports |
| Executed/scored ceilings | (36/5)T and 8T, plus ancillary reserve | same coefficients, smaller T |

The tighter stretch and success comparison reduce the envelope without changing the walk family. The fingerprint change repairs a distribution premise, not a dominant-cost term. No unsupported 1.14-per-iteration substitution or uncharged detector work is taken.

## Unchanged exact class mechanism

Let lambda=p mod r. The verifier checks ord_r(lambda)=12, lambda^6=-1 and 1+lambda^4+lambda^8=0. Thus nonidentity GT has N=(r-1)/12 twelve-element Frobenius classes. In the pinned Fp2[w]/(w^6-(1+u)) tower, write c1=a+b*u, d=(1+u)^((p-1)/6), zeta=(1+u)^((p^2-1)/6) in Fp, B=a*d0+b*d1 and zeta^2=zeta-1. Interleave

    even = (a, a*zeta, a*zeta-a, -a, -a*zeta, a-a*zeta)
    odd  = (B, B*zeta, B*zeta-B, -B, -B*zeta, B-B*zeta).

Products are a*d0,b*d1,a*zeta,B*zeta. K=min(keys), with index kappa and full-conjugate tie-break. K is not injective; full equality confirmation and hard exception caps remain.

Sample independent uniform a_h,b_h modulo r; set J_h=g^a_h*z^b_h and precompute all twelve twists and scalar labels once. For y=g^A*z^D, h=K mod1024 and k=kappa(y), use

    y' = y^2 * J_h^(p^-k)
    A' = 2*A + lambda^-k*a_h mod r
    D' = 2*D + lambda^-k*b_h mod r.

This induces [rep(y)^2*J_h]. It is class-consistent and history-independent. The table remains shared for the complete run. Squaring is charged as a dense product; multiplying the jump is a second dense product. There is no free cyclotomic squaring assumption.

Floyd advances the tortoise once and hare twice. Full class equality gives x=(k*A_h-A_t)/(D_t-k*D_h), k=lambda^(kappa(hare)-kappa(tortoise)); verify the logarithm and original output. Exhaustion or zero denominator fails without restart.

For trajectory entry mu and cycle length ell, first positive meeting is ell*max(1,ceil(mu/ell)) <= mu+ell. Conditioning on initial/jump group elements leaves D0 uniform; its coefficient in den is 2^t*(1-k*2^t). The exact order calculation verifies ord_r(2)/gcd(ord_r(2),12)>2T, excluding coefficient closure throughout this budget. Conditional denominator degeneracy is therefore 1/r. This algebra does not establish the distribution of repeat times.

## Larger faithful evidence: all declared branches live

The new seed base is **2026102306**, with a distinct seed per setting, **10000 independent shared-table/fresh-map pairs per setting**, 70000 total. The 11/10 stretch and ln(5/3) comparison were fixed before the captured run. No seed search, selection of a favorable stretch after observing results, or substitution of refreshed jumps for a retained table was used.

The same tower arithmetic enumerates each small subgroup and its exact Frobenius classes, retaining raw K and full-tie canonical representatives. The exact toy acceleration class_of(2*canonical_exponent+jump[h] mod r) is cross-checked against24 literal metadata and32 literal transition calculations per setting. This is not a hash partition or a production oracle. Cached universes are re-censused for each B.

| q | r | N | B | branch min..max | shared/control repeat-time ratio |
|---:|---:|---:|---:|---:|---:|
|14009539|12289|1024|64|9..25|197321/200344|
|2677963|49177|4098|256|3..29|800713/798942|
|159508651|196657|16388|1024|6..30|1586793/1607663|
|19066387|786433|65536|1024|36..89|797211/806279|
|19066387|786433|65536|8|8146..8258|3204725/3221623|
|19066387|786433|65536|32|1947..2140|1070533/1070804|
|19066387|786433|65536|128|456..561|1074547/1074503|

Every declared branch is reachable; the exhaustive histograms are printed in `evidence/faithful_run.log` and retained in the per-setting JSON files. The exact maximum empirical mean ratio is 800713/798942 <103/100. This is a finite-sample ratio bound on the observed sums, **not a confidence bound on a production mean or a CDF coefficient**.

Let n'=ceil(N*B/(B-1)), t0=ceil(sqrt(2*n'*ln(5/3))) using the rational upper enclosure, and t1=ceil(11*t0/10). Relevant endpoint counts, each out of 10000, are:

| N,B | control repeat<=t0 | shared repeat<=t1 | shared Floyd<=t1 | paired lower difference |
|---|---:|---:|---:|---:|
|1024,64|4180|5123|6663|543/10000|
|4098,256|4101|4723|6316|111/5000|
|16388,1024|4039|4789|6232|7/200|
|65536,1024|3972|4663|6212|291/10000|
|65536,8|4446|5071|6584|9/400|
|65536,32|4151|4697|6175|73/5000|
|65536,128|3959|4681|6136|161/5000|

The paired variable I(shared<=t1)-I(control<=t0) lies in [-1,1]. Subtracting radius1/25 yields the last column; its one-sided Hoeffding failure is exp(-10000*(1/25)^2/2)=exp(-8) per row. The seven-row union is below1/100. Separately, endpoint radius1/50 has two-sided tail2*exp(-8); the fourteen endpoint means have union failure28*exp(-8)<1/100. All lower differences are positive. These finite confidence statements concern independent toy coins only, not production or a guarantee about a deterministic PRNG seed.

B=8 stress reuses entries on624963/640945 transitions (>97%). Complete summaries/digests retain mu, ell, Floyd and branch reuse. Controls are fixed uniform maps, lazily sampled once per state. Identity is separate: uniform N+1-state controls differ slightly from shared uniform-exponent quotient weights.

DP bits immediately above the branch bits were tested at widths2,4,6. Exact census densities, hit times, control statistics and every no-hit-before-repeat censor count are retained; no DP-free cycle is discarded. Full per-setting counts appear in argument.md.

## Honest production extrapolation and constants

**Hmix-strike:** at T=ceil(11*T0/10), the production shared-table raw-coordinate map's first-repeat CDF is at least the uniform-Np-map CDF at T0. Here Np=ceil(N*1024/1023), T0 is defined below. This is the tighter 10% length-stretch hypothesis supported at toy scope by the direct paired CDF tests. It is not proved by the small observed mean ratios. The branch correction is a retained modeling convention, not a theorem for square-and-add walks. Hkey still assumes field-scale exceptional-key frequencies, with exact slow paths and hard caps.

    Np = 4373927678342230290889091420689347427321206073683634826518910598626352813045
    T0 = 66847802276462331670412349757297191499
    T  = 73532582504108564837453584733026910649

For a uniform map on n points, no-repeat through t transitions is product(j=1..t)(1-j/n), bounded above by exp(-t*(t+1)/(2n)). The verifier encloses ln(5/3) using 80 positive atanh-series terms and an explicit geometric remainder, all Fractions. It chooses the exact integer ceiling square root of 2*Np times the upper endpoint. Consequently the comparison success is at least2/5. Taking exactly ln(100/61) would target0.39 but leave no allowance for exceptional losses; that is a diagnostic only, not the candidate success target.

Toy log2(N) spans about10..16, versus production near251.27. Small-B stress increases finite-table reuse but does not prove the asymptotic fixed-1024 family. The 0.40-versus-0.39 difference is a success allowance, **not a bound on transfer error**. Conditional on Hmix-strike/Hkey, negligible identity, denominator and bounded rejection-sampling losses leave success at least0.39. Production correlations violating Hmix can invalidate that success statement. Verifying answers prevents wrong outputs; it does not prove frequent termination.

## Complete cost, fingerprint repair and decision

The ordinary nonreference rows are: four Fp products1532; nine six-limb add/sub/neg576; minimum/index600; two scalar double/add/reduce384; twisted loads/address256; DP predicate/dispatch128; Floyd comparison/control256; loop/meters128; state/scratch512. Sum4372, charged ceiling8000. Fixed-limb carry/subtract/select and four-limb scalar reductions include carry handling; twists are precomputed. These are inherited analytical allowances, not newly measured assembly.

Each Floyd iteration executes three map calls: **six dense products plus24000 instructions**, or36/5 units. The scored envelope deliberately retains8T, distinct from executed work. Ancillary caps are129*M insertions at2^26 instructions, clearing at129*M*128*16 instructions, 2^32 candidate checks at2^15 units, setup2^50 units, and exceptions/output/pairings2^90 units, M=2^69. Their exact sum is below2^100 and their1.2 repricing below2^101. Table construction is charged once within setup, not once per processor and not amortized over targets.

    Wexec = (36/5)*T + 2^100
    W     = 8*T + 2^100
          = 588260661300519118927858079360918490568
    Wstress = (186/25)*T + 2^101.

Optional DP generations s=0..128 use zero bits10..69+s. Fingerprints now apply an independently sampled128-by4608 binary matrix to the complete canonical12-coordinate, six-limb encoding. For distinct encodings, conditioning on all other columns leaves one independent uniform128-bit column, so collision probability is exactly2^-128, independently of the DP predicate and core walk. The73728-byte matrix, full canonicalization, scalar-label normalization and hashing are paid within setup/report ceilings. An exhaustive4-input/3-output analogue verifies all15 nonzero differences over4096 matrices.

Nodes remain128 bytes: fingerprint16, two full scalar labels64, two128-bit links32, metadata16. Hard memory128*M+2^40<2^78 bytes includes matrix, table, code and scratch. Capacity/verification saturation disables DP reporting but never stops or changes Floyd. Distinct-point false-pair expectation is at most129*M*(M-1)/2^129<2^17; Markov bounds false-pair-cap risk by2^-15, not an exact Poisson law or a bound on true duplicate pairs. Neither auxiliary tail is needed for core success. All verification is charged.

The rival ledger reconstructs to125.8845..125.8846. Ours cannot inherit its one-product parallel birthday constant: even zero overhead, zero stretch and the0.39 ideal comparison give a six-reference mu+ell schedule of128.2134..128.2135 bits. This is a floor for that retained sufficient schedule, not a universal lower bound on pairing inversion. Even the unrealistically cheap2800-instruction transplant remains above128.4722 bits at our T. The honest candidate128.82 therefore loses; it is held, not submitted.

## Reproduction and captured scope

From this package: `python3 -B scripts/faithful_audit.py --write`; `python3 -B scripts/toy_audit.py --quick > evidence/algebra_results.json`; `python3 -B scripts/verify_v6.py --write`. The certificate embeds its exact verifier source and individual executed checks; it certifies arithmetic/data consistency only. The rerun checks two production GT samples,2995488 toy alignment cases and495 DP fault comparisons. Old-partition pipeline runs are mechanism evidence only.

Per-setting evidence respects the65536-byte file and262144-byte payload caps. V5 and all pre-existing selectors are preserved. Manifest `experiments: []` claims no organizer run. The blocker is architectural cost: a competitive design must change detector/transition economics and supply new evidence.
