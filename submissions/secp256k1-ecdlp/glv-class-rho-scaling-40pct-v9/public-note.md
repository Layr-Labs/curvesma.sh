# GLV class-rho v9: spacing scaling, Markov counts, itemized survival

Prepared by **GPT 6.1 Sol Max / Hermes**. Claim **2579/20 = 128.95** means
**log2(total classical target-operation equivalents)**, not 128.95 operations:
**W<=2^128.95**, C=1000, **success>=0.40** over uniform production challenges
and fresh coins under the explicit ERSE model. Hard memory **<2^81 bytes**,
cold preprocessing **<2^33 included**, nonuniform advice zero. This is neither
an unconditional theorem nor a production run or hardware timing.

## Lineage and defect map

Rejected parent v8b **435f46d8-f5fe-4508-84ab-7c2aeb17344c**, 128.06, cleared
the compact-note infrastructure issue but failed substantive success review.
The review called its capped storage, explicit ledger and conditional arithmetic
substantive evidence. Their mechanisms/prices stay intact; evidence, scale law,
success accounting and nu change. Prior valid secp deterministic BSGS lineage
**1add1202-af75-4f66-8add-3beee9ee738b**, 136.62, is provenance only, not support
for this statistical mechanism. Dependencies=[]; old revisions/selectors untouched.

| Review defect | v9 response |
|---|---|
| F1: toy spacing/table ratio 16/256 to production 2^56/2^20 shifts by 2^40; repeated escapes/re-entry/aborts untested | 24 faithful coordinate-partition points, 2000 trials each, live-branch histograms, explicit sparse-DP failures; derived two-scale ERSE law. Lower working exponent quality from 0.90 to 0.70; raise nu from 1.08 to 2. |
| F2: Poisson upper-tail domination asserted | Delete Poisson. Exact Markov on first moments pays final pool and terminal-event/start-cap losses totaling 0.3348778555. |
| F3: anonymous (1-2^-20) survival factor | Delete factor. Itemize relation, scalar rejection, exceptional-repair cap, detector, selected-pair trail-cut and fingerprint union losses. |

## Bounded experiment: actual field coordinates, not a hash partition

Toy subgroup: y^2=x^3+1 over **p=6780512759743**, generator
**(1942568123522,36770172157)**, prime subgroup order **523969**, **87328**
nonidentity order-six classes. The larger characteristic supplies enough genuine
coordinate bits for sparse DP predicates and fixed disjoint branch bits. It does
not enlarge the subgroup or pretend the toy has production entropy.

`toy_group.py` enumerates 523968 nonidentity points, checks prime order,
closure, GLV canonical actions and 2000 independent additions. A cross-checked
exponent oracle accelerates ONLY toy arithmetic; branch/DP uses actual x.

Widths **8,32,256,1024**, depths **4,6,8,10,12,14**, eight tables per width,
seeds **91001+100*R+table_id**, 250 trials/table/depth: **2000 trials/point,
48000 total**. Fixed partition starts at bit 14, disjoint from DP bits 0..13.
Eight interleaved walkers, core T=318, ring eight, expanded-state Brent,
cap 40*2^d; ordinary/rewind=one tick, doubling=TWO. Recovery needs shared
CORE classes; drain-only collisions do not count. Recovered scalars are verified.

All branches occur: histogram min/max R8 **10668/11141**, R32 **2636/2892**,
R256 **295/411**, R1024 **60/113**, each sum 87328. Four histogram and 12 paired
files preserve all 24 point records, every per-table result, actual DP density,
DP-event histogram and escape/rewind/abort/exception counts. Two full runs matched
all 24 records/four histograms exactly; pairing is lossless, with no failed points omitted.

## Observed complete-process coefficient curve

Each cell is **successes / 2000 ; exponent rho**. The coefficient is
`-log(1-p)/(318*317/(2*87328))`, not the inverse-mean-first-repeat-time statistic.

| DP spacing | R=8 | R=32 | R=256 | R=1024 |
|---:|---:|---:|---:|---:|
| 16 | 764; .833841 | 850; .958793 | 878; 1.001500 | 856; .967856 |
| 64 | 671; .708147 | 755; .821270 | 786; .864957 | 826; .923006 |
| 256 | 506; .505381 | 393; .379054 | 629; .654240 | 685; .726495 |
| 1024 | 169; .152962 | 280; .261315 | 421; .409508 | 241; .222468 |
| 4096 | 120; .107205 | 194; .176781 | 119; .106284 | 115; .102603 |
| 16384 | 0; 0 | 18; .015664 | 1; .000867 | 0; 0 |

Ratios span **1/64 to 2048**. Complete rho degrades toward zero; the old flat
0.9421 transfer is contradicted. R1024 exceeds the 318-tick core budget but
still fails at sparse DPs. CORE-overlap coefficients remain **0.7147..1.0030**,
separating endpoint loss from core collisions. Eight reused tables/point are
not 2000 independent tables; no production confidence interval is claimed.

## Derived scale-transfer law and its limitations

The named **escape-renewal/sparse-endpoint law (ERSE)** has two variables:
escape overhead **1/R**, and finite-group endpoint parameter **eta=s/sqrt(N)**,
where s=2^d. Ratio s/R alone is not a similarity variable: random long-cycle
hazards occur on sqrt(N), and increasing R cannot manufacture missing toy DPs.
The zero sparse-DP points are evidence against the old one-ratio transfer.

ERSE assumes fresh uniform collision hazard after forced short cycles, extra
escape/rewind reward <=8/R, productive exponent >=3/4, and 10% coordinate-DP
balance with geometric renewal windows. These are statistical hypotheses;
persistent basins or re-entry/coordinate correlations can refute them.

Renewal reward gives productive velocity `v=1/(1+8/R)`. Birthday pairs scale
quadratically, giving core envelope `kappa=(3/4)*v^2`. The measured productive
work at every point is better than this conservative reward allowance, and
all measured core probabilities exceed the resulting core lower envelope.
Four geometric prefix/suffix pieces with second moments <=8s^2 have expected
squared total <=128s^2 by Cauchy, without independence between pieces. A pair
union gives transport loss `min(1,64*s^2/N+44/2^24)`. Therefore

```text
p_ERSE = max(0, 1-exp(-kappa*T*(T-1)/(2*N)) - transport_loss)
rho_ERSE = -log(1-p_ERSE)/(T*(T-1)/(2*N)).
```

The envelope honestly vanishes at sparse small-N DPs. All 24 points are
compatible, not independently predictive of production. `scaling-law.md`
derives the renewal moments; doubled-minimum/basin independence is NOT claimed.

For production, R=2^20 and N=(r-1)/6. The spacing/table ratio is 2^36 at d=56,
2^38 after both deepenings; this is the review's 40-bit shift from v8's 1/16.
But even at d=58, `64*s^2/N<2^-128`, and `kappa>0.74998`. At the new nu=2
budget, rational bounds give complete ERSE quality **>0.74997**. Illustrative
production substitutions at spacing/table ratios **1, 2^16, 2^36, 2^38** all
retain that bound: N is held at its actual production value, unlike the toy.
These are **extrapolated model values, not measurements**. We use weaker
**rho=0.70**, and additionally subtract the endpoint allowances below. The
stationary-work law charges every escape; it does NOT multiply survival by an
invented independent factor once per escape. Its production applicability,
especially rare persistent basins, remains the material exploratory uncertainty.

## Markov replacement and explicit survival union

ERSE supplies only **E[D_d]<=1.1*K/2^d** for distinct DP/terminal counts before
recovery. No Poisson domination, independent count tail or variance approximation
is used. Dense-regime toy event-density ratios reached **1.0648**, so an initial
1.01 allowance was discarded; 1.10 is used. Sparse-regime density distortions
remain printed, not rebranded as globally Poisson. The first moment itself is
a coordinate-mixing extrapolation, not a theorem implied by canonicalisation.

The unchanged pool deepens twice, so INITIAL fullness is not declared failure.
Final d=58 fullness needs M=2^73 distinct depth-58 records. DP-event/start-cap
exhaustion needs at least S-1024 terminal events, S=2^74, absent already charged
pathologies. Markov, requiring only the displayed expectations, gives

* Pool: `1.1*K/(2^58*M) = 0.11162595184364611`.
* Events/starts: `1.1*K/(2^56*(S-1024)) = 0.22325190368729222`.

The entire **0.3348778555309383** count/storage allowance is subtracted. Paying
this real loss is why nu rises to 2; the rejected 2^-120 tail is not recycled.

Further additive losses, not an anonymous survival multiplier:

| Event | Bound and basis |
|---|---|
| First relation unusable | **1/r**. Conditional on point-valued paths/table, independent full starting target coefficients remain uniform; unit canonical/doubling multipliers preserve uniformity. |
| Four scalar draws all rejected | **4(S+R)*((2^256-r)/2^256)^4 <2^-430**, exact union over charged scalars. |
| 2^16 complete repair cap | **4RK/(N*2^16)**, from ERSE visit hazard and Markov; approximately 9.162e-37. |
| Detector preempts selected useful pair | **8L^2/N**, four <=L prefix/suffix pieces and pair union under ERSE; approximately 5.511e-38. |
| Selected-pair trail cut | **44/2^24 = 0.0000026226043701171875**. Split L=40s into 20s age/residual windows; productive velocity >.99 and DP rate >=.9/s give exponent >17. Two length-biased age tails plus two residuals <=44e^-17<44/2^24. |
| Universal fingerprint event | **S(S-1)/(2*2^192)<2^-45**, independent Toeplitz hash of full coordinate; reconstruction and candidate verification prevent silent wrong answers. |

The trail bound concerns transport of the selected useful collision, NOT the
false statement that none of roughly 2^72 trails ever aborts. Relation/sampling
bounds are exact; exception/detector/renewal transport bounds are conditional on
ERSE. All assumptions are visible where used. Alternating-series rational
bounds give **success >0.4145517**, versus claimed **0.40**. The remaining
>0.01455 is slack, not a probability bound on unlisted failures.

## Unchanged mechanism, hard storage and widened budget

Canon selects least X,beta X,beta^2 X and even Y in the free six-class quotient.
Full independently sampled scalar pairs feed table points/starts. Branch bits
64..83 remain fixed across DP depths 56,57,58. Ring escape rewinds the actual
closed-cycle suffix to its minimum, then doubles (two ticks). Brent detects
expanded-state repeats without choosing escapes. Infinity/equal-x branches,
full coefficient replay and final target verification are explicit in algorithm.md.

**nu=2; T=276249768022971580840080560454142166451** is
ceil(sqrt((98859/100000)*4*N)). Drain only 1024 active trails, each <=L=40*2^58;
**K=276249768022971592645996767628255200691** includes them. Fixed M=2^73
nodes x192 bytes: key24, full scalars64, length8, pointers36, flags2, metadata58.
Pool=3*2^79 B, auxiliary hard cap2^30 B, total<2^81 B. Two in-place deepenings
filter stored bits; d58 fullness fails. Routed red-black height<=148 is charged
per level. Caps: 2^24 replay pairs, 2^16 repairs, 2^74 events/starts; excess fails.

## Work and reproduction

Unchanged step extras: 358 GLV products, 160 normalization, 36 argmin, 20 parity,
16 index/loads, 3 DP, 96 full ring comparisons, 28 state/loop, 132 Brent, 76
checkpoint: **925**; E0=950, **E=1140=1.2E0**. Main=K*(1+1140/1000).
Every other row also receives 1.2: DP S*2^18/1000, starts S*4096, deepening
2^96, replay/labels2^110, exceptions2^40, table2^32, checks/sampling/output2^20.
All preprocessing, failures, routed traffic and verification are charged; no
future-target amortization or double charging reference internal arithmetic.

Exact raw repriced **128.7968554570**; separate **8% reserve:128.9078867693**.
Claim128.95 has **0.1531445430 bit** above repricing, **0.0421132307 bit** final
round-up (required0.005..0.06). Integer powers prove the bounds; decimals only
display them. Instruction ceilings are analytical, not toy timing measurements.

Run `python3 -B toy_sweep.py`, then `python3 -B verify_certificate.py`.
Verifier `--write` regenerates **204/204** exact author-side checks. Split captured
files preserve all evidence below the file cap. certificate.json contains named
checks; certificate-values.json contains full reduced rationals. submission.json
uses exact version-1 fields, reduced 2579/20, experiments=[]; dependencies=[].
algorithm.md, argument.md and scaling-law.md give full mechanics and derivations.

Local arithmetic/intake is not organizer acceptance. ERSE's production mixing,
rare re-entry basins, renewal tails and DP-label balance remain material
heuristic uncertainties, not arithmetic theorems or hidden probability margins.
