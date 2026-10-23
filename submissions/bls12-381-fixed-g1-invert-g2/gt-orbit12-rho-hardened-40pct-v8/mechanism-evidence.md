# V8 actual-field mechanism evidence and deviation audit

Attribution: **GPT 6.1 Sol Max / Hermes**. All measurements below are from the
captured seed-83 run, not predictions or production measurements. The raw capture
is split losslessly between `mechanism-output.json` and `mechanism-pipeline.json`
to satisfy the 64-KiB per-file intake limit.

The completed run took **441.380 seconds**, observed **488 cycles**
and used **240 independently sampled shared tables**. Eight additional
basin-restart probes each reached exactly 32 cycles on their fixed table.
The 248 extra cycles are not counted as independent maps.

## Faithfulness contract

The implemented field is Fp2[w]/(w^6-(1+u)), u^2=-1. Each p is prime, p=19 mod 24,
and xi=1+u is nonsquare and noncube in Fp2. Each r is prime, divides p^4-p^2+1,
and p mod r has order exactly twelve. The projected nonidentity generator has
order r. This is an ACTUAL order-12 Frobenius field experiment, not an order-six
stand-in and not a minimum over a SHA-256 hash of visible exponents.

The representative is chosen by the minimum real w coefficient among all
twelve Frobenius conjugates, with the exact even/odd four-product construction
and the same full-encoding tie rule. Coordinate ordering uses reduced one-limb
Montgomery values. R=1024, h=K mod 1024. Marking is precisely
`((K >> 10) & ((1 << d)-1)) == 0`, at d=8,10,12; production substitutes d=60.
The 46--50-bit p values avoid the vacuous high-zero-bit mistake that a tiny p
would cause. The same jump table is shared by both walkers and all density
comparisons in each table trial. Square-and-add, inverse Frobenius twists,
label normalization, three-evaluation Floyd steps and restart decisions are
the same mechanism as production. Visible exponents provide an instrumentation
oracle for equality/recovery ONLY; they never choose the partition or mark.

The toy differs in p,r, word width, number of workers (two), density and budget.
It does not implement production 112-bit fingerprint storage, G2 output or
instruction prices. Those remain in the exact v7 mechanism/ledger. Reproducible
Python random coins are for measurement, not the claimed production coin source.

## Parameter and first-cycle distribution table

| log2 N | p | r | N | independent tables | lambda min / median / max | rho min / max |
|---:|---:|---:|---:|---:|---|---|
| 16.978643 | 70368891220219 | 1549741 | 129146 | 128 | 2 / 177 / 912 | 10 / 1100 |
| 22.109456 | 70377839719483 | 54298861 | 4524906 | 64 | 47 / 1099 / 6797 | 248 / 7394 |
| 26.197334 | 70381155765259 | 923346397 | 76945534 | 32 | 479 / 5260 / 18344 | 2901 / 26889 |
| 33.044421 | 73986475718707 | 106302407641 | 8858533971 | 12 | 500 / 51336 / 166734 | 21832 / 206866 |
| 38.552248 | 1044763899747811 | 4836869906233 | 403072492187 | 4 | 60562 / 327986 / 342138 | 501427 / 904895 |

Each first-cycle observation runs until an actual repeated class, recording
mu, lambda and rho=mu+lambda. The largest observed rho is 904895, not 81.
The small-size long trials run to 32 detected cycles and fresh-start out of
their basins. No purported escape edge is added to the stationary map.

## Exact-rational empirical CDF grid

For each of the 240 independently sampled maps, normalize by its own sqrt(N).
The comparator is a lazily generated uniformly random function on that same N,
not a fitted asymptotic curve. Each displayed comparison is an integer count.
Membership at threshold a/b is tested as t^2*b^2 <= a^2*N, with no floating
rounding. The sample sizes are unequal across fields; this is the disclosed
fixed mixture, not 240 observations at the largest size.

| quantity | threshold / sqrt(N) | concrete count / 240 | random-function count / 240 | absolute difference |
|---|---:|---:|---:|---:|
| rho | 1/2 | 29/240 | 15/240 | 7/120 |
| rho | 1/1 | 95/240 | 89/240 | 1/40 |
| rho | 3/2 | 157/240 | 162/240 | 1/48 |
| rho | 2/1 | 206/240 | 207/240 | 1/240 |
| rho | 3/1 | 236/240 | 237/240 | 1/240 |
| lambda | 1/2 | 117/240 | 112/240 | 1/48 |
| lambda | 1/1 | 181/240 | 190/240 | 3/80 |
| lambda | 3/2 | 219/240 | 224/240 | 1/48 |
| lambda | 2/1 | 229/240 | 238/240 | 3/80 |
| lambda | 3/1 | 239/240 | 240/240 | 1/240 |
| mu | 1/2 | 128/240 | 124/240 | 1/60 |
| mu | 1/1 | 194/240 | 188/240 | 1/40 |
| mu | 3/2 | 226/240 | 221/240 | 1/48 |
| mu | 2/1 | 237/240 | 230/240 | 7/240 |
| mu | 3/1 | 240/240 | 239/240 | 1/240 |

The maximum empirical GRID differences are **7/120 for rho, 3/80 for lambda,
and 7/240 for mu**. These bound only this captured grid discrepancy; they are
not suprema over every threshold and NOT population/production probability
bounds. Large-field counts 12 and 4 are too small for a rare-basin conclusion.
The conservative v8 success ledger reserves a larger **1/10 absolute deficit**
on the production collision-and-delivery event, explicitly as an extrapolation
premise H, not as a theorem deduced from 7/120.

## Partition occupancy and density scaling

Uniform subgroup elements are sampled independently of the walk tables. The
exact Pearson occupancy statistic is R*sum(c_j^2)/m-m, with R=1024 and m=32768.
At the smallest subgroup, repeated class sampling contributes finite-population
overdispersion; no claim of an exact chi-square null law is made. In particular
the largest value 1253.8125 is disclosed rather than described as a clean
p-value. Counts below are raw, so sparse-count uncertainty remains visible.

| log2 N | uniform marks at d=8 / 10 / 12 | nominal counts | Pearson statistic (1023 nominal df) | first-path distinct states | first-path marks at d=8 / 10 / 12 |
|---:|---|---|---:|---:|---|
| 16.978643 | 111 / 34 / 12 | 128 / 32 / 8 | 20061/16 | 57829 | 241 / 49 / 15 |
| 22.109456 | 127 / 35 / 6 | 128 / 32 / 8 | 4123/4 | 173309 | 702 / 166 / 47 |
| 26.197334 | 127 / 30 / 7 | 128 / 32 / 8 | 16319/16 | 375040 | 1428 / 347 / 91 |
| 33.044421 | 121 / 36 / 4 | 128 / 32 / 8 | 1069/1 | 1444809 | 5725 / 1432 / 390 |
| 38.552248 | 129 / 38 / 8 | 128 / 32 / 8 | 8497/8 | 2902457 | 11114 / 2770 / 725 |

Pooled independent-element counts: **[615, 173, 37]** out of **163840**, versus nominal
640,160,40. These are compatible in scale but do not certify exact independence.
First distinct paths contain **4953444 states**, with mark counts **[19210, 4764, 1268]**.
On the largest field, mark-free first cycles are 0/4 at all three densities;
this is a four-trial observation, not a rare-event upper bound.

## Two-worker pipeline, bounded windows and basin costs

Each window has S_toy=ceil(sqrt(6N/5)), the smaller predecessor budget, and
begins with fresh starts; the table remains fixed across windows. The windows
are not independent map trials. Reported per-table spread is retained. All
three evaluations and two dense-operation analogues per evaluation are counted.
DP predicates and class/Floyd tests use the actual minimum-coordinate data.

| log2 N | d | success windows | reports | Floyd detections | cycle/basin iterations / all | verified relations per (T/sqrt(N)) | verified relations per (T/sqrt(R)) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 16.978643 | 8 | 168/512 | 976 | 120 | 11819/201728 | 0.372324 | 209/6304 |
| 16.978643 | 10 | 128/512 | 170 | 147 | 15805/201728 | 0.286814 | 161/6304 |
| 16.978643 | 12 | 146/512 | 22 | 187 | 20486/201728 | 0.333132 | 187/6304 |
| 22.109456 | 8 | 117/256 | 2637 | 60 | 24459/596736 | 0.573916 | 23/2664 |
| 22.109456 | 10 | 103/256 | 456 | 105 | 56368/596736 | 0.434893 | 61/9324 |
| 22.109456 | 12 | 99/256 | 126 | 109 | 63694/596736 | 0.417069 | 13/2072 |

Totals: **2395392 tortoise iterations**, **7186176 F evaluations**,
**14372352 dense-operation analogues**, **4387 reports**, **728 Floyd detections**,
zero timeout restarts, zero equivariance failures and zero incorrect verified
recoveries. The exact checks recompute these accounting relationships.

R in sqrt(R) here means jump-table width 1024, so sqrt(R)=32 exactly; the
birthday normalization is sqrt(N), displayed separately. Relation counts can
contain multiple recoveries in one window and are NOT independent birthday
collisions. No inverse-mean-cycle statistic is substituted for a success exponent.

**Unsuccessful configurations are material:** at the smallest N, successful
window fractions are only 168/512,128/512,146/512; the d=10 per-table counts
range from 1/128 to 69/128. At the next N the fractions are 117/256,103/256,
99/256. We do not hide the last value below 0.39 or claim every toy trial
meets the production floor. Two workers, sparse marks and a short end-of-window
cutoff cause unrecovered/coalesced paths and table/basin variation. This is why
the pipeline is mechanism and regime evidence, not an estimator of production
success. The mu/lambda census separately samples the stationary map beyond
those finite windows.

The cycle/basin row pessimistically calls every iteration before a Floyd
detection lost, even when its relation actually recovers x. This makes the
restart burden visible without treating useful meetings as free work.

## Analytic interpretation and explicit extrapolation

The excluded inverse-undo family would require 2*lambda^k=-1 and hence
2^12=1 mod r, which fails at production and every toy r. This excludes that
family ONLY; random two-cycles still exist (minimum observed lambda=2).
There is no claim that deleting one family proves a random-map theorem.

With R fixed, the observed cycle medians grow 177,1099,5260,51336,327986 as
N grows. Counts with lambda<=32 are 12/128,0/64,0/32,0/12,0/4. Counts with
lambda<=1024 are 128/128,31/64,2/32,1/12,0/4. Thus the captured behavior is
not a persistent O(R) lifetime obstruction across the tested sizes. The CDF
grid and mark-free-cycle counts quantify that statement without extrapolating
a false factorial tail or confusing a point-on-cycle law with basin mass.

The transfer variables are lambda/sqrt(N), rho/sqrt(N), q*sqrt(N), and
Lmax/sqrt(N). Absolute q is not the only sparsity scale. For the largest toy,
q*sqrt(N) is approximately 2480,620,155 at d=8,10,12. Production has
log2(N)=251.2721269, q*sqrt(N)=2^65.6360635 and Lmax/sqrt(N)=2^-58.6360635.
If the scaled cycle law transfers, production marks are much MORE frequent
relative to cycle time, despite their much lower absolute density. The
experiment reaches log2(N)=38.5522484, not 40 exactly or production 251.272.

The new budget sqrt(8N/5) has ideal success 1-exp(-4/5) >11/20. Allow a
GLOBAL 1/10 probability deficit for concrete-map/mark/delivery extrapolation,
the exact finite-fingerprint 1/32 loss, and the 2^-30 ideal-model remainder:
`11/20 - 1/10 - 1/32 - 2^-30 > 2/5`. This permits substantially more empirical
deviation than the observed 7/120 grid difference while charging it explicitly.
H asserts the production event deficit is at most 1/10; the samples cannot
prove this or quantify every rare basin at 2^67 steps. It is the remaining
material exploratory premise, supported now by a faithful mechanism and
multiple scaling regimes rather than by an unrelated whole-exponent hash.

A full production statistical-distance or tail bound cannot be inferred from
finite toy data. No such derivation is claimed. The empirical deviation bound,
its finite-sample scope, the 0.10 transfer tolerance and the gap in log2(N) are
all explicit so the exploratory reviewer can evaluate the complete premise.
