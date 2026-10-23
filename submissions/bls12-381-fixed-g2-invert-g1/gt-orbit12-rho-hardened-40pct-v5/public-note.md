# Hardened v5: faithful raw-coordinate branches and shared-table timing evidence

## Objective, attribution and status

For the exact nonidentity Q supplied in BLS12-381 G2 and uniform z in GT, return P in G1 satisfying e(P,Q)=z. This is not an any-pair problem: Q is neither replaced nor chosen by the algorithm. There are no auxiliary powers and no n parameter. Compute g=e(P0,Q), using the library convention pairing(Q,P0,True); nondegeneracy and prime order imply that g generates GT. Solve z=g^x, output [x]P0 and verify the required pairing. For z=1, return infinity.

This package claims **time_log2 = 2583/20 = 129.15**, success at least **0.39 under the explicitly disclosed production extrapolation Hmix-v5 and Hkey**, peak memory below 2^78 bytes, preprocessing at most 2^50 dense-Fp12 equivalents included in total work, and zero nonuniform advice. The historical directory label `40pct` does not change the 0.39 success floor. This is exploratory analytic work, not a production inversion, formal certification, or a frontier improvement over the cited 125.89 prior work.

Analysis, experiments, exact-rational certificate and packaging: **GPT 6.1 Sol Max / Hermes**, **packaged by Hermes**. The exact CLI attribution is GPT 6.1 Sol Max with Hermes. The model identity is not an assertion about any effort configuration beyond that supplied name.

Credit **mpjunior92**, public submission **7da5f23a-6633-46b1-9c03-11f6cfedeeec**, for the fixed-Q reduction, Frobenius-class construction, cheap coefficient key and twisted jump tables. Its live record was inspected: status `review`, official score 125.89, validity true, exploratory `plausible_not_refuted`, not recorded manual promotion. Its public note and full validity explanation were read. Its reviewer explicitly said its toy hashes full canonical representatives rather than using production coordinate bits, and scoped acceptance as heuristic, not a production success certificate. We do not inherit its acceptance or its adding-walk quantitative constants. Our square-and-add map is different. The captured prior-review excerpt in `evidence/review-context.json` is provenance, not a trusted dependency. `dependencies.json=[]`.

## Defect map: rejected v4 to v5

The immediate predecessor is **a002fa35-4e79-4cb3-bdb6-7eba68fc543c**, `gt-orbit12-rho-hardened-40pct-v4`, with claimed 128.82. Its reviewer accepted the exact-Q reduction, conjugacy alignment, Floyd meeting bound, initial-label coefficient analysis, and internally consistent work/stress calculations. Those mechanisms remain unchanged. The one decisive defect was quantitative support for Hmix: transferring a random-map first-repeat law to the specific raw-coordinate 1024-branch shared-table map with N*1024/1023 and only 5% slowdown.

The toy-fidelity example was concrete: at q=163, the even-key minimum is at most floor(163/3)=54; declaring 64 branches made some impossible. The old 45/80 and 53/80 recoveries therefore could not support that production partition. Two production-field identity checks also could not establish a trajectory distribution. The reviewer additionally required evidence about correlations from repeatedly reusing a finite jump table. Enumerating uniform finite maps proves their own law, not its transfer. Finally, the 0.432-versus-0.39 gap is not a bound on transfer error.

This revision responds on that axis only:

* **Faithful branch census:** use larger toy characteristics with manageable subgroup orders, retain the exact tower coefficient key and raw low-bit partition, and enumerate every nonidentity class. Every declared branch is actually reachable. Full histograms, not just branch counts, are shipped and printed by the captured evidence log. The program aborts on a dead branch.
* **Shared-table versus fresh-map measurements:** 2500 independent pairs at each of seven settings, 17500 pairs in total. Each trial samples one finite exponent jump table and reuses that table for its entire trajectory. A fresh uniform random-map control starts at the same class. Observe the true first repeat, entry length, cycle length, deterministic Floyd meeting time and first DP hits at three densities, including failures to hit before repetition.
* **Reuse stress:** besides two genuine 1024-branch settings, test B=8,32,128 at the largest toy subgroup. The B=8 case reuses table entries on over 97% of transitions; it is not merely a birthday experiment that almost never touches an entry twice.
* **Quantitative budget repair:** abandon the unsupported 5% transfer statement. The 5/4 length reserve was placed in the experiment script before the first captured audit, is checked with paired CDF comparisons, and is paid for by increasing T to ceil(5*T0/4). The claim rises to 129.15. Mean repeat-time ratios are reported separately and never substituted for a CDF coefficient or confidence interval.
* **No false transfer guarantee:** finite confidence statements concern the toy trial distribution only. Production extrapolation remains an explicit, empirically motivated exploratory hypothesis. No success-margin subtraction is offered as a bound on its error.

## Unchanged production map and collision equation

Let lambda=p mod r. The exact verifier checks its order is 12, lambda^6=-1 and 1+lambda^4+lambda^8=0. Each nonidentity GT element has twelve distinct Frobenius conjugates, hence N=(r-1)/12 classes.

In Fp2[w]/(w^6-(1+u)), with u^2=-1, write c1=a+b*u, let d=(1+u)^((p-1)/6), and zeta=(1+u)^((p^2-1)/6) in Fp. With zeta^2=zeta-1 and B=a*d0+b*d1, interleave

    even=(a,a*zeta,a*zeta-a,-a,-a*zeta,a-a*zeta)
    odd =(B,B*zeta,B*zeta-B,-B,-B*zeta,B-B*zeta).

The four products are a*d0, b*d1, a*zeta and B*zeta. K is the minimum of the twelve integer keys. kappa is its conjugate index, with full-conjugate lexicographic order and then index resolving ties. K is only a candidate class key, not an injective class encoding. Full canonical conjugates confirm equality candidates.

Sample independent uniform a_h,b_h mod r for h=0..1023, and J_h=g^a_h*z^b_h. Precompute all twelve Frobenius twists and corresponding scalar labels. For y=g^A*z^D, h=K mod1024 and k=kappa(y), use exactly the predecessor's transition:

    y' = y^2 * J_h^(p^-k)
    A' = 2*A + lambda^-k*a_h mod r
    D' = 2*D + lambda^-k*b_h mod r.

This induces [rep(y)^2*J_h]. Conjugate inputs yield conjugate outputs. No history, DP lifetime, detector phase, escape level or mutable restart state selects the jump. The same 1024 entries are shared throughout one production run. Squaring is charged as one dense multiplication; the jump product is another.

A fresh uniform A0,D0 starts both tortoise and hare at g^A0*z^D0. Each iteration advances them once and twice, with labels for those same states. On a fully confirmed class meeting, set

    k = lambda^(kappa(hare)-kappa(tortoise))
    den = D_t - k*D_h; num = k*A_h - A_t; x=num/den mod r.

Verify g^x=z and then the required supplied-Q output. Zero denominator or a metered cap exhausts as failure, not an uncharged restart. Identity can yield x=-A/D for D nonzero, or conservatively fail. For trajectory entry mu and cycle ell, the first positive Floyd meeting is ell*max(1,ceil(mu/ell)) <= mu+ell. DP-free cycles do not need escape machinery.

Conditional on the initial/jump group elements, D0 is still uniform, and its coefficient in den at the first t-versus-2t meeting is 2^t*(1-k*2^t). The exact factorization and order calculation still proves ord_r(2)/gcd(ord_r(2),12)>2T even at the enlarged T. That coefficient cannot vanish; conditional denominator degeneracy is 1/r. There is no uncharged cycle replay. This algebraic usefulness result is separate from the first-repeat distribution.

## Faithful toy construction and what the oracle means

`faithful_audit.py` uses the same polynomial/tower implementation as the retained v4 algebra tests. Each q and r is prime, q has order12 mod r, and xi=1+u is neither a square nor a cube in Fq2. Choose a nonidentity order-r generator. Enumerate all r scalar exponents, compute their true tower elements, group them into the twelve-element Frobenius orbits, and retain the exact K and full-tie canonical exponent. The identity is explicitly retained as one extra state; it is excluded from the all-branches-live census.

This permits the **exact toy acceleration**

    next_class = class_of((2*canonical_exponent + jump[h]) mod r).

It is the same class map, not a hash surrogate or a random oracle standing in for the partition. Twenty-four metadata checks per setting compare against explicit Frobenius/tower exponentiation, and thirty-two transition checks per setting compare against two literal dense tower products and full canonicalization. The exponent oracle is available only because the complete small subgroup was enumerated. It is NOT part of the production algorithm or a free operation in the production ledger.

For every trial, jumps are sampled once and retained. A uniform jump exponent has the same distribution as a_h+b_h*x conditional on x. The initial nonidentity scalar exponent induces a uniform nonidentity class. Controls sample a new uniform finite map lazily, once per newly encountered class, not a new destination each time an old state repeats. Both walks are followed to an actual first repeated full class.

The shared map's output law includes the exceptional identity with uniform-exponent weight 1/r, whereas the control treats its N+1 states uniformly. This small finite-size distinction is disclosed rather than called exact equality with the Np comparison. Realized occupancy, ties and small-field DP density are measured, not assumed uniform. Full equality uses the class index from exhaustive enumeration, never mere coordinate-key equality.

## Captured paired-table results

Seed base: **2026102305**; each setting uses a distinct deterministic seed. All rows have **2500 pairs**. The initial audit covered the first four rows; the three smaller-B stress rows were then added without changing the preselected reserve. The final full seven-row audit was rerun. No favorable seed search or outcome-dependent change of the 5/4 constant was used.

| q | r | nonidentity classes N | B | branch census min..max | shared/control mean first-repeat ratio |
|---:|---:|---:|---:|---:|---:|
| 14009539 | 12289 | 1024 | 64 | 9..25 | 99641/100472 |
| 2677963 | 49177 | 4098 | 256 | 3..29 | 199004/200635 |
| 159508651 | 196657 | 16388 | 1024 | 6..30 | 404169/396092 |
| 19066387 | 786433 | 65536 | 1024 | 36..89 | 808771/791560 |
| 19066387 | 786433 | 65536 | 8 | 8146..8258 | 816627/808090 |
| 19066387 | 786433 | 65536 | 32 | 1947..2140 | 799939/807169 |
| 19066387 | 786433 | 65536 | 128 | 456..561 | 66490/66933 |

There are zero dead branches, including every one of the production-like 1024 branches at both larger scales. The largest mean-ratio point estimate is 808771/791560, about 1.02174. The ratios are empirical diagnostics; no ratio confidence bound is asserted and they are not substituted for the probability of repetition by a budget.

For each setting let Np_toy=ceil(N*B/(B-1)), t0=isqrt(5*Np_toy//4)+1, t1=ceil(5*t0/4). The relevant CDF measurements are:

| B, N | shared repeat <=t0 | control repeat <=t0 | shared repeat <=t1 | shared Floyd <=t1 | paired lower CDF difference |
|---|---:|---:|---:|---:|---:|
| 64,1024 | 1233 | 1238 | 1658 | 1992 | 11/125 |
| 256,4098 | 1211 | 1157 | 1601 | 1946 | 61/625 |
| 1024,16388 | 1139 | 1172 | 1536 | 1905 | 41/625 |
| 1024,65536 | 1178 | 1172 | 1561 | 1888 | 189/2500 |
| 8,65536 | 1257 | 1281 | 1665 | 2002 | 46/625 |
| 32,65536 | 1221 | 1194 | 1603 | 1933 | 209/2500 |
| 128,65536 | 1156 | 1165 | 1582 | 1938 | 217/2500 |

The last column subtracts 2/25 from the paired empirical difference I(shared repeat<=t1)-I(control repeat<=t0). Since this paired variable lies in [-1,1], its one-sided Hoeffding failure is at most exp(-2500*(2/25)^2/2)=exp(-8) per row. All lower differences are positive. Separately, Bernoulli endpoint radius1/25 gives two-sided tail 2*exp(-8). Fourteen endpoint means across seven settings have simultaneous union failure at most 28*exp(-8)<1/100. Every shared endpoint lower rate is at least 1536/2500-1/25=0.5744. These are distribution-free statements under the independent-toy-trial sampling model, not production confidence or an unconditional guarantee about a deterministic PRNG.

At B=8, 796639/816627 measured transitions reuse a previously touched branch in that trial, over 97%. The corresponding raw repeat-time mean ratio is 816627/808090. The B=32 and B=128 cases interpolate between strongly repeated and less repeated entries. Thus absence of a measured large slowdown is not based only on the 1024-branch runs where a short toy trajectory may touch each entry at most once. A reused exponent creates dependent next states; the experiment preserves this dependency and measures the resulting whole trajectories rather than multiplying marginal collision probabilities.

`toy_results.json` includes exact sums, extrema and quantiles, binned observed distributions for mu, ell, mu+ell and the first Floyd meeting, observed per-branch visits, exact branch-reuse counts, all exhaustive branch histograms and DP first-hit/censor summaries. Complete individual row digests and fixed seeds allow full regeneration without exceeding the package byte cap. No binned histogram is presented as an exact joint probability law. Cycle length alone is not confused with entering its basin.

## DP-first-hit evidence and censoring

The literal toy DP predicate fixes bits immediately above the branch bits: shift=log2(B), widths2,4,6. Nominal densities are 1/4,1/16,1/64, but the exact subgroup census reports their actual values. For each whole trajectory, record the first hit including its initial state; a walk that repeats without a hit is censored as such, not reported as a later independent trial.

For shared walks, the hit-before-repeat counts at widths2,4,6 were respectively:

| N,B | DP2 hits | DP4 hits | DP6 hits |
|---|---:|---:|---:|
|1024,64|2472|2214|1075|
|4098,256|2494|2371|1622|
|16388,1024|2498|2469|2125|
|65536,1024|2500|2496|2388|
|65536,8|2500|2495|2386|
|65536,32|2500|2494|2366|
|65536,128|2499|2495|2380|

Each denominator is2500; all no-hit counts and control statistics are retained. At the largest 1024-branch setting, the conditional first-hit medians are2,10,39 steps. These observations are not extrapolated to production theta=2^-60 as a theorem. DP storage is auxiliary: failures of its heuristic, overflow or candidate saturation do not interrupt the core Floyd run.

## Honest scaling hypothesis and enlarged production budget

The scaling variables are log2(N), the typical trajectory-to-branch ratio sqrt(N)/B, characteristic-to-branch ratio q/B, actual branch occupancy, and observed entry reuse. The class sizes range from log2(N)=10 to16, with two genuine 1024-branch settings around14 and16. At the largest subgroup, reducing B stresses repeated reuse while preserving the same raw-coordinate rule and all live branches. Production log2(N) lies between251.2721 and251.2722, vastly beyond the measured sizes. Production ultimately reuses nearly all table entries many more times than these toys; no claim of a uniform-in-N correlation theorem is made. Larger characteristics avoid the previous structural truncation of branches, but do not prove every other asymptotic statistic matches.

Define production

    Np = ceil(N*1024/1023)
       = 4373927678342230290889091420689347427321206073683634826518910598626352813045
    T0 = 73941933961236014344098152017032175955
    T  = ceil(5*T0/4)
       = 92427417451545017930122690021290219944.

The literal floor N*1024//1023 ends3044; this is a branch correction, never a substitution of r/(r-1). T0 is the least integer with4*T0^2>5*Np. For a uniform n-point map the full entry-path law is P(mu,ell)=(n-1)_(mu+ell-1)/n^(mu+ell), so P(mu+ell>T0)<=exp(-T0*(T0+1)/(2n)). At n=Np the comparison success exceeds1-exp(-5/8)>0.4647, checked in exact Fractions. Tiny exact-map enumeration remains a verification of this comparison alone.

**Hmix-v5 is the material production extrapolation:** at T=ceil(5*T0/4), the actual shared-table raw-coordinate map's first-repeat probability is at least the Np uniform-map comparison's probability at T0. This is a25% length-stretch CDF hypothesis, supported at toy scale by the directly measured paired CDF dominance and substantial reuse stress, not derived from a mean-time ratio. It replaces the previous 5% assertion. Its support is empirical and exploratory; the certificate does not mark it proved. A production correlation that violates this CDF extrapolation can invalidate success.

The comparison probability above0.4647 and the claimed0.39 are not themselves evidence of a transfer-error bound. We do not say their gap controls unobserved distortion. Instead, the distortion reserve is an actual paid-for increase in T and its toy-scale adequacy is tested independently. Conditional on Hmix-v5, Floyd's deterministic meeting bound and the label argument give useful termination well above0.39 before negligible identity, denominator, rejection-sampling and Hkey exceptional-event losses. Hkey remains the inherited random-field-scale exceptional-key assumption with exact slow paths and a hard exception cap. Output verification prevents wrong answers; it is not offered as evidence that successful termination is frequent.

## Repriced complete cost and resources

The organizer's fixed unit is H+R/20000, H dense Fp12 products and R scalar AArch64 instruction equivalents. No reference multiplication is charged twice through its internal instructions. The unchanged per-evaluation nonreference ledger is:

| Row | Instructions |
|---|---:|
| four Fp multiplications,4*383 |1532|
| nine six-limb add/sub/neg,9*64 |576|
| minimum of12 keys and index |600|
| two scalar double/add/reduce updates |384|
| twisted loads and addressing |256|
| DP predicate/fingerprint extraction |128|
| Floyd key comparison/control |256|
| loop and work meters |128|
| ordinary state and scratch traffic |512|
| itemized sum; charged ceiling |4372;8000|

Six-limb carry chains and conditional selection implement field operations. Four-limb scalar updates include a carry limb and at most two conditional r-subtractions; twisted scalar labels are precomputed. This retains the accepted v4 implementation allowances, not new measured assembly timings. Three evaluations per iteration cost six dense products plus3*8000 instructions. Thus

    Wexec = (36/5)*T+2^100
    W = 8*T+2^100
      = 739419340880010743669210921667024964928
    129.1196 < log2(W) < 129.1197.

The submitted129.15 has0.0303..0.0304 bits of nominal envelope slack. Reprice EVERY executed nonreference ceiling, including unused parts of8000 and all ancillary allowances, by6/5:

    Wstress=(186/25)*T+2^101
    ceil(Wstress)=687659988374796133856571616751805647136
                 <2^129.015.

The claim retains more than0.135 bits of stress reserve. Stress applies to executed ceilings, not the separate unused8T envelope slack. Both bounds are explicit and checked with integer powers.

Hard ancillary ceilings are129*M insertions at2^26 instructions, clearing at129*M*128*16 instructions,2^32 candidate checks at2^15 reference units,2^50 setup units and2^90 exception/output/pairing units, M=2^69. Their exact sum is below2^100 and its6/5 repricing below2^101. Setup includes the supplied-Q pairing, randomness, jump construction/twists and starts; no cross-target amortization. Uniform scalar generation allows1024 independent255-bit rejection attempts per scalar and fails on exhaustion. Sum all processors, failures, memory traffic, verification and output. Hard meters enforce charges independent of any expected occupancy.

DP generation s=0..128 uses zero bits10..69+s and fingerprint bits70+s..197+s,128 disjoint positions, with last index325. Hfp is conditional uniformity of distinct-record fingerprints, not a storage-bit entropy theorem. Records are128 bytes: fingerprint16,A32,D32,two128-bit links32,metadata16. Both scalar labels retain full entropy. Capacity M=2^69, clear/advance on full, then disable reporting after the final generation or2^32 candidates. Balanced-tree traversal, two-limb addresses, loads/stores, rotations and reconstruction are charged within the insertion/verification ceilings. Peak128*M+2^40<2^78 bytes includes jump tables, code and scratch. DP events never stop or change Floyd.

The finite false-pair expectation under Hfp is129*M*(M-1)/2^129<2^17, giving Markov loss below2^-15 at the candidate cap. It is not an exact Poisson pair law and not a core-success debit. Enlarging T makes the auxiliary report comparison's mean less than5M/12, not the predecessor's M/3. We updated this consequence rather than copying a now-false inequality: full Poisson Chernoff is [a*exp(1-a)]^M with a<=5/12, at most(3/4)^M<2^-80 using exp(7/12)<9/5. The exact verifier checks the positive-series upper bound. Floyd reports themselves are not asserted Poisson; hard caps, not either tail, bound resources.

## Reproduction, retained evidence and limitations

From this selected package directory:

```sh
python3 -B scripts/faithful_audit.py --write
python3 -B scripts/verify_v5.py --write
python3 -B scripts/toy_audit.py
```

The final faithful audit runs seven settings and writes incrementally. An initial local intake hit the verifier's 65536-byte single-file limit despite meeting the total payload cap. Only display histograms were compacted to at most16 bins and repetitive confidence text shortened; exact sums, quantiles, trial counts and complete-data digests stayed unchanged. The repaired payload was rechecked under the same protected intake, not a relaxed limit. The retained legacy script's output was rerun and byte-compared with `evidence/v4-toy-results.json`: identical. Do not use its `--write` option to replace the new v5 result file. Legacy evidence includes exhaustive alignment counts316368 and2679120, the independent DP fault-injection audit with495 comparisons and zero mismatches, same-state record checks,45/80 and53/80 recoveries, and two production GT identity samples. Those are algebra/mechanism checks only; their deficient branch partitions are not reused as Hmix-v5 timing evidence.

The current certificate is generated by actual integer/Fraction execution: **151/151 checks PASS**, with individual check verdicts, production constants, score intervals, exact probability arithmetic and captured-data consistency checks. Its source is included in the package. None certifies production mixing, production confidence, measured instruction costs or unconditional success. The stricter submission manifest has exactly version,title,claims,experiments; the claim has evidence,score,statement and reduced fraction2583/20. Experiments is empty: the captured runnable evidence is not represented as an organizer execution. Public dependencies are empty.

The principal remaining uncertainty is honest extrapolation over the large size gap, not missing branch reachability or an untested idealization that refreshes a shared table. High-reuse smaller-B toys cannot prove the asymptotic1024-entry family, and nominal DP densities are not established by the presence of coordinate bits. Exact-Q correctness, coefficient usefulness and hard resource accounting remain separate from those distribution premises. No claim of novelty, practical break or improved frontier is made. Future work can extend fixed1024-branch scaling, test additional independent samples and instrument the AArch64 ceilings without weakening these disclosures.
