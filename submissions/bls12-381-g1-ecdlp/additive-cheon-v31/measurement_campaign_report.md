# P2 coprime-defect measurement campaign

**Status: reproducible finite measurements; production P2 remains unproved.** Seven original saved controls and one direct dense low-load complement have exact original, unpruned coprime defects. A denser actual-order-d=2^20 construction now has its complete critical-load prefix cache. Whole-row Horvitz–Thompson (HT) samples estimate defects of fixed finite critical-load targets; their exact sample counts are not exact full-target counts. All three registered HT summaries are complete below; no uncomputed target count is imputed.

Snapshot UTC: 2026-10-25T07:53:39.539698+00:00. All scientific sources, computations, records and large artifacts remain on EC2.

## 1. Endpoint and interpretation

Write P=x*y_prime and Q=x_prime*y for ordered product indices, with x,x_prime in B and y,y_prime in G. M0 counts all ordered choices with gcd(P,Q)=1, without a congruence. E_coprime counts their ordered coverage-index collisions, summed over the six roots. The endpoint is

~~~text
T=n^2, N=6*n^2, m=(r-1)/d, h=m/6
B0=36*n^4/m=N^2/m
Delta0=E_coprime-36*M0/m
Delta0/B0=(m*E_coprime-36*M0)/(36*n^4)
theta=M0/n^4
eta_coprime=E_coprime/(36*M0/m)-1
Delta0/B0=theta*eta_coprime
~~~

The exact counter groups integer products by U-coset modulo r, then tests ordinary integer gcds. Each unordered distinct coprime product pair contributes 12 to E_coprime: two orientations and six root pairs. Equal integer products cannot be coprime here. This defect is neither total collision excess nor the complete P2 epsilon. Other residual strata remain separate.

The campaign records a production coprime allowance of 53/1250=0.0424 of B0. Every saved exact control lies below that number. This is a finite observation, not a transfer theorem. The fixed production epsilon<=0.057 goal remains unresolved.

## 2. Seven exact controls

Every constructor retains full49: NA=floor(196*n*d/X2)+2 and M=ceil(NA/2). Baby rows are the first M physical interval primes, giant rows the next M. Genuine character classes determine complementary increasing prime lists; complete rows are consumed except the final required prefix. Empty physical rows are skipped without moving or replacing primes.

Quartic mini-BLS means r=z^4-z^2+1 and L=z^2-1. Quadratic C3 means r=L^2+L+1; it preserves the small cube-root identity but L+1 need not be square. These families remain separate.

| Control | Family | r | d | X1 | X2 | n | Achieved K |
|---|---|---:|---:|---:|---:|---:|---:|
| quartic-d4-full | quartic | 116975041 | 16 | 262144 | 256 | 819 | 0.550485436893 |
| quartic-d8-full | quartic | 2517580801 | 256 | 262144 | 8192 | 950 | 0.550623836979 |
| quartic-d12-full | quartic | 167771750401 | 4096 | 524290 | 262144 | 1938 | 0.550173804135 |
| quartic-d16-full | quartic | 164995450798081 | 65536 | 8388610 | 4194304 | 15192 | 0.550033465987 |
| quartic-d20-low-K-screen | quartic | 175921860024729601 | 1048576 | 536870912 | 16777216 | 11825 | 0.00500073285103 |
| C3-d20-sparse-full | C3 | 2326566556139521 | 1048576 | 67108864 | 16777216 | 14262 | 0.550042879515 |
| C3-d20-rich-full | C3 | 12354112538542081 | 1048576 | 134217730 | 67108864 | 32864 | 0.550022498218 |

Six points use n=ceil(sqrt((11/20)*m/6)); rounding explains loads slightly above 0.55. Only quartic-d20-low-K-screen uses K=0.005 before rounding. It is not pooled with critical-load observations.

### Exact counts and signed defects

| Control | E_coprime | M0 | Exact Delta0/B0 | Percent of B0 | theta |
|---|---:|---:|---|---:|---:|
| quartic-d4-full | 2034036 | 413360542550 | -285121610/449920319121 | -0.0633715789% | 0.918741663763 |
| quartic-d8-full | 2938932 | 802856677406 | -6964153/407253125000 | -0.001710030586% | 0.985697380967 |
| quartic-d12-full | 12296724 | 14005039299990 | -87247845/87076321928 | -0.1001969802% | 0.992817082329 |
| quartic-d16-full | 761142600 | 53227852649654126 | 1018643634437/26633596216117248 | 0.003824656746% | 0.999261463186 |
| quartic-d20-low-K-screen | 4189776 | 19545852694263790 | -4016273518038/3910520737578125 | -0.1027043145% | 0.999654726617 |
| C3-d20-sparse-full | 671020764 | 41361550488069378 | -1137666750037/10343362300191684 | -0.0109990032% | 0.999712407041 |
| C3-d20-rich-full | 3563314368 | 1166151243988044800 | 54208838765/2804066810494976 | 0.001933222082% | 0.999708056945 |

The largest positive value is quartic d16, approximately 0.00382466% of B0; rich C3 d20 is approximately 0.00193322%. Signed negative values are retained, not clipped. The low-d theta deficits are material: small Delta0/B0 need not mean equally small conditional eta_coprime. Both quantities are saved.

### Density, concentration and wrap

| Control | X2/d | pi(X2)/d | Max weighted side delta | Nominal log2 wrap | Actual-max log2 wrap |
|---|---:|---:|---:|---:|---:|
| quartic-d4-full | 16 | 3.375 | 0.2407399953 | 25.198374507 | 23.741171310 |
| quartic-d8-full | 32 | 4.015625 | 1.755439335 | 30.770609064 | 29.152558027 |
| quartic-d12-full | 64 | 5.61523438 | 14.96590274 | 36.712302150 | 34.930035009 |
| quartic-d16-full | 64 | 4.51579285 | 23.96308109 | 42.770581112 | 40.916973613 |
| quartic-d20-low-K-screen | 16 | 1.02793789 | 182.0106765 | 48.712287624 | 46.721179326 |
| C3-d20-sparse-full | 16 | 1.02793789 | 150.7049656 | 48.952876118 | 47.033920119 |
| C3-d20-rich-full | 64 | 3.77446079 | 154.0802301 | 52.544159147 | 50.569295930 |

Weighted delta is d*sum_c W(c)^2/n^2-1 for the actual emitted side. Both sides pass the finite P1 threshold 500 in all seven controls. This does not establish production P1. Nominal wrap is (X1*X2)^2/r; actual wrap uses the largest emitted element on each side.

All seven satisfy X1*X2<r and X2<X1/2. All seven fail both 3*X1^2<r and 2*X1^2*X2<r. Production exclusions using these inequalities cannot be imported. Coprime collisions were instead counted directly. The old d20 points have X2/d only 16 or64 and therefore do not meet the dense-class request by themselves.

## 3. Dense actual-mask critical-load target

The new dense target is quartic and uses the complete full49 construction:

~~~text
z=288768
r=6953384735058167857153
L=z^2-1
g=30
d=1048576
m=6631264433916252
h=1105210738986042
X2=34359738368=2^35
X1=68719476738=2*X2+2
n=24654937
T=n^2=607865918473969
NA=147474, M=73737
~~~

Here n is the least integer satisfying 120*n^2>=11*m. It is a full critical-load target. Its inherited configuration identifier ends in LOWLOAD, but the explicit n override and invocation flags govern the actual experiment. The earlier suggested n=54772 low-load screen is a different construction.

The exact stream computes chi(p)=p^m mod r for every list prime p<=X2, retaining complete lists for 137601 requested complementary classes from the fixed physical block rows. The cache uses exact class offsets rather than full n-element B/G arrays. The 3 GB cache guard is operational, not an equidistribution theorem.

The dense-critical-full invocation is now PASS. It completed in 241.697769seconds. Its cache size is 1553945400bytes and SHA-256 is f3bf24652cc47a9a7a4a197b59f4a8e839b79f656e8e35c58f3bb2281328d8a6. The completed invocation and cache metadata are authoritative. Cache construction is not a full coprime energy count.

The full target contains about 6.079e14 B×G products, so its full pair count is not enumerated at the declared RAM budget. HT samples whole rows of this actual critical target; it does not substitute a low-load construction. Intermediate actual-mask full-critical targets at X2=2^27 and2^31 are also built for the density ladder.

### Direct exact dense low-load complement

A separate predeclared top-density input uses the same r,d,X1,X2 but n=54772, T=2999971984, NA=329 and M=165. Its exact stream-derived prefixes use 39 baby rows and 40 giant rows, with 330 distinct physical block classes and no integer shared-prefix swap. It has K approximately 2.71439e-6 and does not satisfy P1<=500. This is a direct low-load exact count, distinct from the critical target and its HT estimates. It cannot repair the load/concentration mismatch by normalization.

The direct exact count is complete at runs/dense-lowload-X2-35/dense.result.json:

| Quantity | Exact saved value |
|---|---|
| E_coprime | 46032 |
| M0 | 8543838322330758792 |
| Delta0/B0 | -8082691603730571/1124978988098112032 |
| Decimal Delta0/B0 | -0.00718474895019609 |
| Percent of B0 | -0.71847489502% |
| theta | 1067979790291344849/1124978988098112032 |
| Total coverage collision energy | 48456 |
| Integer-product L2 | 2999971984 |

Actual side weighted deltas are 5046875390199/187498249 and 5043654557943/187498249 (approximately 26916.9201 and 26899.7422), both far above 500. The direct point therefore complements the critical-load HT estimates but is not a P1-qualified critical-load control. Its exact integer-product L2=T confirms no repeated integer products in this input. Sorting took 200.408695seconds and peak RSS was 35160444KiB.

The compact package retains the extraction record and its input/row hashes, the exact counter result and a heavy rebuild recipe. The 109544 emitted prime factors, direct .in/.rows files, full prime caches and event streams are not required by the mandatory checker and are omitted from the compact bundle. Saved-result arithmetic is checked without those heavy inputs.

## 4. Exact compressed M0 and HT design

m0_compressed.py computes exact M0 from distinct large-row IDs, class keys and initial-prefix lengths in O(R log R) time and O(R) memory. It uses incidence/Möbius identities for shared small and large primes and reads no modular collision counts. Primality, true classes and prefix semantics remain the exact constructor’s contract. Explicit-prime and brute-force tests check the compressed formula, including its four-row kernel.

For R_B/R_G full emitted rows, independently select k_B/k_G rows uniformly without replacement and retain every selected row’s original prefix. Coprimality forces two distinct baby and two distinct giant rows. Every eligible tuple therefore has common inclusion probability

~~~text
pi_22 = k_B(k_B-1)/(R_B(R_B-1))
      * k_G(k_G-1)/(R_G(R_G-1))
Delta_hat=(E_coprime_sample-36*M0_sample/m)/pi_22
Normalized estimate=Delta_hat/B0_full
~~~

Under that ideal uniform sampling design, both counts are sums over the same full-target coprime tuples, so the estimator is design-unbiased for Delta0. Dividing by fixed full-target B0 preserves this identity. The checked invalid hybrid is E_coprime_sample/pi_22-36*M0_sample/m: its random sample baseline was left unscaled. In contrast, E_coprime_sample/pi_22-36*M0_full/m is also unbiased, but is a different estimator. Sample-dependent B0 and pi_11 also give the wrong identities. ht_package_check.py exhausts 101 small genuine-prefix row samples across five designs and verifies this identity plus those three negative controls. Its small oracle is not a full49/critical-load production rung.

The actual implementation uses predeclared fixed seeds with random.Random(seed).sample. Sizes depend on full-mask metadata and a predeclared pair budget, before collision outcomes. Exact row positions, inclusion probability, Python version and hashes are registered. Each sample count is exact; the full-target value remains an estimate.

Pohlig–Hellman quotient logarithms are independently verified; the native H/U counter sorts uint64 quotient keys plus uint32 base indices, and every emitted collision is checked again by its U ratio and integer gcd. Six predeclared row-mask replicates are used per target. Replicate mean/SE are empirical variability descriptions, not a finite-sample rigorous one-sided guarantee or a production confidence statement.

### Completed finite HT summaries

All three registered finite critical-load targets have six completed samples. Exact sample counts do not make the target estimates exact. No replicate was excluded or clipped.

| Target | n_full | X2/d | R_B,R_G | k_B,k_G | Replicates | Mean HT Delta0/B0 | Empirical SE |
|---|---:|---:|---|---|---:|---:|---:|
| quartic-d20-X2-2pow27-LOWLOAD | 124013 | 128 | (17035, 17091) | (381, 382) | 6 | -0.00457170207346 | 0.00921443607757 |
| quartic-d20-X2-2pow31-LOWLOAD | 972259 | 2048 | (9690, 9713) | (121, 121) | 6 | 0.00969608593372 | 0.00407272178336 |
| quartic-d20-X2-2pow35-LOWLOAD | 24654937 | 32768 | (17464, 17465) | (35, 35) | 6 | 0.00752891866178 | 0.00923015852009 |

Exact full-target mass and concentration metadata from preregistrations:

| Target | theta=M0/n^4 | Weighted delta baby | Weighted delta giant | P1 both<=500 |
|---|---:|---:|---:|---|
| quartic-d20-X2-2pow27-LOWLOAD | 0.999866066015 | 70.1246526693 | 70.0861983977 | True |
| quartic-d20-X2-2pow31-LOWLOAD | 0.999791733132 | 109.176505674 | 109.199148025 | True |
| quartic-d20-X2-2pow35-LOWLOAD | 0.999885398769 | 60.173044293 | 60.1351985131 | True |

#### quartic-d20-X2-2pow27-LOWLOAD

Exact mean fraction: -63290887473514466768166180/13844053364909566832393871101.
Exact inclusion probability: 117064282/470864352329145.
Exact full-target M0: 236488857902354513482.
Sample normalized HT values in preregistered seed order:

| Seed | Exact normalized HT fraction | Decimal |
|---:|---|---:|
| 2026102501 | -168379200719851715981078040/13844053364909566832393871101 | -0.0121625651304 |
| 2026102502 | -185761104330677698964010780/13844053364909566832393871101 | -0.0134181153044 |
| 2026102503 | -207856953633776262988874310/13844053364909566832393871101 | -0.0150141687665 |
| 2026102504 | -239198436820778256165434400/13844053364909566832393871101 | -0.0172780637661 |
| 2026102505 | 571371166601872477838181000/13844053364909566832393871101 | 0.0412719563802 |
| 2026102506 | -149920795937875344347780550/13844053364909566832393871101 | -0.0108292558535 |

Source summary: ht-runs/quartic-d20-X2-27/summary.json.

#### quartic-d20-X2-2pow31-LOWLOAD

Exact mean fraction: 34595822013818306820857257/3568019327624052907970756073.
Exact inclusion probability: 39930/1677377730607.
Exact full-target M0: 893382476186366321402168.
Sample normalized HT values in preregistered seed order:

| Seed | Exact normalized HT fraction | Decimal |
|---:|---|---:|
| 2026102501 | 8655348419968701545906560/1189339775874684302656918691 | 0.00727743963125 |
| 2026102502 | 75474855548403170616505216/17840096638120264539853780365 | 0.00423063042087 |
| 2026102503 | 450141042249263735141373056/17840096638120264539853780365 | 0.0252319845223 |
| 2026102504 | -1424473660054637762793104/540608989033947410298599405 | -0.00263494260908 |
| 2026102505 | 62079116425918916254843298/3568019327624052907970756073 | 0.0173987612526 |
| 2026102506 | 196761297470347505086276/29487763038215313289014513 | 0.00667264238441 |

Source summary: ht-runs/quartic-d20-X2-31/summary.json.

#### quartic-d20-X2-2pow35-LOWLOAD

Exact mean fraction: 28139351270618563532866087359072/3737502360528872489599505150100515.
Exact inclusion probability: 10115/664426359505488.
Exact full-target M0: 369458629575683560877465018672.
Sample normalized HT values in preregistered seed order:

| Seed | Exact normalized HT fraction | Decimal |
|---:|---|---:|
| 2026102501 | 103537383795272285115959368512/747500472105774497919901030020103 | 0.000138511462747 |
| 2026102502 | 86780842012854075200312210728896/3737502360528872489599505150100515 | 0.0232189397201 |
| 2026102503 | -86094046102497937590910945074048/3737502360528872489599505150100515 | -0.0230351817331 |
| 2026102504 | 72423066887714309794486432886592/3737502360528872489599505150100515 | 0.0193773969624 |
| 2026102505 | 136297447442932581886666771867968/3737502360528872489599505150100515 | 0.0364675214342 |
| 2026102506 | -41088889536268009518937743097536/3737502360528872489599505150100515 | -0.0109936758757 |

Source summary: ht-runs/quartic-d20-X2-35/summary.json.

### Top-target design-variance diagnostic

The posthoc top-target diagnostic partitions the coprime count into groups with two baby and two giant rows. There are 354025 such groups in each 35-by-35 row sample. Joint inclusion probabilities depend on the overlap of two groups on each side. The source uses all row cells, including zero observed collision cells, with their exact M0 centering. This is not an iid occupancy variance model.

Its individual design-variance estimators are unbiased under the ideal uniform row design, and can be negative. One of six signed estimates is negative and is retained without clamping. For independent replicate row designs, the estimated variance of the six-replicate mean is exactly sum(V_i)/36. Independence is the design idealization; fixed pseudorandom seeds do not prove a mathematical RNG-independence theorem.

~~~text
sum(V_i) = 4479336357989333134838979197013551439875690334795110426020180243224/1690239791290026168714013120047421250693893645736746033458066995092225
estimated variance(mean) = 373278029832444427903248266417795953322974194566259202168348353602/5070719373870078506142039360142263752081680937210238100374200985276675
descriptive design SE(mean) = 0.00857988420250914
empirical SE(mean) = 0.0092301585200935
~~~

Both SEs are descriptive; no normal, Studentized, distribution-free or production confidence interval is certified. The saved per-seed variance JSON files, their hashes and exact aggregation are included in the arithmetic check. The variance algebra is audited from saved reports; the mandatory check does not redo the all-cell variance computation.

### Root-order decomposition

The following entries are contributions to the total normalized HT mean, rather than relative excesses against separate family baselines. The sample baselines are6*c_order*M0_sample/m with c_order=(1,1,2,2) for orders(1,2,3,6). Their exact sum equals the corresponding total mean.

| Target | Order1 contribution | Order2 contribution | Order3 contribution | Order6 contribution |
|---|---:|---:|---:|---:|
| quartic-d20-X2-2pow27-LOWLOAD | -0.00982329838601 | 4.16945959054e-06 | -0.00189375163158 | 0.00714117848454 |
| quartic-d20-X2-2pow31-LOWLOAD | 0.00106388978059 | 0.00343398537422 | 0.0011043291912 | 0.00409388158771 |
| quartic-d20-X2-2pow35-LOWLOAD | -0.00434475868299 | 0.00921211548321 | 0.00427792748865 | -0.00161636562709 |

The full diagnostic includes every sample’s exact order fractions, coprime-event hashes and verified sums. These profiles do not establish a persistent production order-family mechanism.

## 5. Remaining production scale gaps

| Quantity | Production | Dense critical top |
|---|---:|---:|
| log2 r | 254.857089413047 | 72.558205409578 |
| d | 2^20 | 2^20 |
| X2/d | 2^50 | 2^15 |
| log(X1)/log(r) | 0.321748946395 | 0.496153395703 |
| log(X2)/log(r) | 0.274663734728 | 0.482371356932 |
| nominal log2 wrap | 49.142910586953 | 69.441794590506 |
| K before ceiling | 0.55 | 0.55 |
| 3*X1^2<r | true | false |
| 2*X1^2*X2<r | true | false |

The dense full-critical top repairs the order-d/load mismatch at its finite modulus, and raises X2/d by 9 bits beyond the old rich d20 control. It still leaves 35 bits in X2/d and about 182.30 bits in field size, substantial relative-length exponent gaps, and about 20.30 extra nominal-wrap bits. Actual observed maxima must remain separate from cap-based wrap.

The tier ratio also differs: production X1/X2=4096, whereas the dense target has X1/X2=34359738369/17179869184, just above 2. Its physical 147474 block primes span 34359738421 through 34363323997, giving normalized total block width 5.21769979953e-05. The emitted baby/giant spans are only about 6.22e-6 and 6.15e-6 of X1. In production, even the elementary spacing of distinct odd physical block primes forces total physical span/X1 at least2(2M-1)/X1=0.00487468618686, over 93.426times the dense physical width. No relative-tier or block-width bridge has been proved.

The completed dense sieve records pi(X2)=1480206279 and mean class length 1480206279/1048576. The actual full-target maximum-product wrap is log2=67.441902784133; the nominal cap wrap remains 69.441794590506. These are separate measurements. Production’s49.142910586953 value in the table is its nominal cap, not a measured actual maximum.

No count or sample here proves that a new correlation cannot appear across these gaps. First-eligible deterministic moduli are controls, not random production-modulus replicates. Changing primitive generator only relabels classes; conjugating the cube root leaves U unchanged. Those are not independent observations. Different d, class means, family, load and wrap are not pooled as one fitted distribution.

At K=11/20 the elementary full49/full-load geometry obstruction in family_configs.json gives, assuming enough physical block primes, X1>2X2 and2X1^2X2<r:

~~~text
n < X1*X2/(392*d)
X2/d > 422576/15
n^2/d^2 > 830055571086708736/50625
At d=2^20: n^2 > 18027767942915247477684137
~~~

This explains why practical exact enumeration cannot match all production geometry. HT alters the computational budget, not the geometry itself. A reviewer may assess an explicitly named empirical extrapolation premise from these records; the premise remains a premise. Production P2, including epsilon<=0.057, is not proved.

The iid-N collision null is unsuitable: U creates six-point orbits. Even for independent uniform quotient draws, Var(E)=72*T*(T-1)/h*(1-1/h), approximately 12*B0, rather than the iid-N approximation 2*B0. Actual Cartesian products, class masks and prefixes add dependence. Exact deterministic-control counts have no Monte Carlo counting error, and HT design variability is distinct from either occupancy null. No iid-N z scores or p-values are reported.

### Suggested fixed-instance extrapolation premise

For the single declared production instance

~~~text
r=52435875175126190479447740508185965837690552500527637822603658699938581184513
d=2^20, X1=2^82, X2=2^70, K=11/20,
with the original full49 block rule and increasing class-list prefixes,
assume E_coprime <= 36*M0/m + (53/1250)*B0.
~~~

E_coprime and M0 refer to the original unpruned coprime tuples. This is a fixed-instance empirical extrapolation premise supported by exact finite controls and design-unbiased estimators of full critical-load toy targets. It is not a theorem, a uniform modulus/character bound, or a statistical confidence statement about production. The other declared P1 and P3/P4 antecedents and proved reductions remain prerequisites of the Route B coverage floor; this premise does not replace them.

In particular, no claim is made that a new mechanism or correlation cannot arise above X2/d=2^15. The production ratio remains 2^50, and the geometric and relative-field-size gaps above remain substantive.

## 6. Mandatory lightweight reproduction

~~~sh
python3 measurement_arithmetic_check.py
~~~

This command prints JSON to stdout and leaves the frozen package files unchanged. If saving output, choose a fresh path outside the verified package.

This stdlib wrapper reuses package_checker.py for four genuine small full49 constructions, exact expected counts and negative controls. It runs ht_package_check.py for exhaustive small HT design identities. It calls existing heavy_arithmetic on all seven saved rungs and arithmetic-audits any additional completed dense log-counter result, then checks exact load rounding, block formulas, geometry, theta*eta normalization, native saved counts and allowance flags. Completed HT summaries receive exact arithmetic, inclusion-probability and seeded-position checks.

Latest fresh combined audit: PASS in 0.290321 seconds on this EC2, including eight saved exact rungs (seven original plus the direct dense point), three completed HT summaries and one all-six estimated-variance aggregate. This is an observed timing; the command enforces a five second default limit. Only small oracles are freshly counted. Saved heavy counts are arithmetic checked, not rerun. The mandatory path does not scan a dense list, allocate a pair array or hash the multi-GB cache.

## 7. Optional heavy offline reproduction

Use a separate sufficiently provisioned EC2 working directory containing the package; these commands write results and must not overwrite frozen evidence. They are not part of the mandatory check.

### Original seven exact controls

~~~sh
python3 family_search.py
python3 run_u64_campaign.py
~~~

The runner compiles construction.cpp, rebuilds exact masks/prefixes, checks emitted H membership and invokes native_coprime_runner.py/native_coprime_counter.cpp. The largest saved 24-byte counter array is25,921,019,904bytes. Per-run wall time, CPU, peak RSS, binary/source and input hashes are retained under runs/.

### Dense full-critical cache

~~~sh
python3 dense-family-search.py
python3 dense_builder_runner.py --case quartic-d20-X2-2pow35-LOWLOAD \
  --n 24654937 --output dense-critical-full --mode cache --threads 8 \
  --max-cache-bytes 3000000000 --compile
~~~

The runner compiles GMP/OpenMP using the recorded unprivileged SDK headers and existing libgmp.so.10. Preserve version/package hashes. This visits about 1.48 billion primes; benchmark forecasts do not replace the completed cache invocation. Intermediate full-critical caches use the same runner with their own registered r/X1/X2 and critical n.

### Registered whole-row HT campaign

~~~sh
python3 ht_campaign.py --prefix dense-critical-full \
  --config dense-family-search.json --case quartic-d20-X2-2pow35-LOWLOAD \
  --out ht-runs/quartic-d20-X2-35 --pair-cap 3000000000 \
  --baseline-target 200000 --reps 6
~~~

The actual saved registration, sample-size rule, output path and baseline target take precedence over this example; intermediate targets use their own declared baseline target. Registrations must precede new counts. Packed 12-byte sample records need at most 36 GB at3b products; 16-byte records need 48 GB. Sorting must account for every additional full buffer. The full-target cache and all heavy arrays remain on EC2.

## 8. Source and artifact manifest

Freeze a manifest after all registered work is complete. Each entry should contain SHA-256, relative path, size and role. Include constructor/counter sources and runners; family-search sources/configs/certificates; m0_exact and compressed-M0 sources/tests; package/HT/checker sources; dense builder and pinned SDK details; Pohlig–Hellman/log-counter sources/tests; registrations; this report; and completed saved results/summaries. Retain compiler commands and environment/library versions.

Every heavy run already records source, binary, input and row hashes. Dense invocations record metadata/cache hashes. HT registrations pin metadata/config/source/cache-reader hashes, exact row positions, seeds and Python version. Preserve the large cache as an optional offline artifact with its hash in the small manifest. A full cache digest verification is an explicit heavy integrity operation, not a subsecond mandatory step.

Fresh heavy reproduction is performed in a new output directory. Frozen invocation/metadata/result files contain timing, machine/path and historical source/metadata hashes, so their full JSON hashes need not match a fresh run. Deterministic mathematical outputs to compare are the emitted prime prefixes, cache bytes where format/order are fixed, integer counts, exact fractions and registered row selections for the specified seeds and environment. Historical registration hashes are evidence for that frozen run; reusing them with freshly generated metadata would correctly fail immutability checks. A fresh registration is written before its new collision counts.

The manifest must not contain its own digest internally; store its digest externally or in an archive receipt. Paths stay within the package root. Existing package_checker manifest verification rejects escaping paths and altered files, with a digest negative control. Hash agreement establishes artifact identity, not production extrapolation validity.

Portable parameter-search sources now prefer the byte-identical local prime_count_inputs.py helper. parameter_portability_audit.json records the only path-resolution edits, before/after hashes, and fresh certificate equality for all7original and6dense parameter records. Frozen campaign configurations and active HT sources were not changed. The complete reverse unified diffs recover both historical parameter-source files byte-for-byte with their original provenance hashes; no scientific result was silently recalculated under changed parameters.
