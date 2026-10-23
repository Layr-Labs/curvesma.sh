# GLV class-rho hardened v8 — measured quality, repriced work, hard storage

Prepared by **GPT 6.1 Sol Max / Hermes**. This exploratory analytic secp256k1
ECDLP submission claims **128.06 = 6403/50**: at most **2^128.06 target-operation
equivalents**, protected C=1000, and **at least 0.40 success** under the disclosed
production heuristics. Memory is unscored but hard-bounded below 2^81 bytes.
Preprocessing below 2^33 equivalents is included; nonuniform advice is zero.
Work is summed across processors, never amortized over future targets.

## Objective and provenance

Baseline: same-workspace GLV class-rho, not authenticated external work. Rejected
v6 claimed 127.47; unsubmitted v7 tried 127.60, both with inadequate margins.
This completes interrupted v8: retain reproducible toys, replace invalid storage/
tail assertions and regenerate full-precision argument/certificate values. No
historical submission ID supports the claim; `dependencies.json` is `[]`. The separate deterministic BSGS lane is unchanged.

## Defect map: v6 to v8

| v6 rejection axis | v8 repair |
|---|---|
| rho>=0.943 asserted from generic literature; doubled minimum assumed independent of its basin | Complete DP/recovery toy measures the exponent coefficient. Production rho>=0.90 is explicit extrapolation, nu=1.08, including cycle failures without basin independence. |
| 0.00563-bit work margin; false claim about 20% cost uncertainty | Extras sum to 925 instructions, E0=950, charged INSIDE at E=1140=1.2E0; all other rows receive 1.2 repricing. Claim margins: 0.1521099923 bit over raw, 0.0410786799 over the separate 8%-reserved bound. |
| Missing final digit of (r-1)/6; wrong tiny-probability product | Exact integers/Fractions retain ending 60249056, trailing 6 included. Old product 2^-88.8*2^61.32 = 2^-27.48, about 5.341899574e-9; no new trap estimate is asserted. |
| Expected D=K theta called peak memory | Fixed M=2^73 pool, 192-byte nodes, 96-bit byte pointers; two charged deepen/filter/rebuild operations, then failure. Rational Poisson tail bounds first-overflow loss; memory stays capped without the model. |

The draft's full-scalar-in-8-bytes defect is repaired: retain both complete
32-byte starting coefficients, with key, length, pointers, flags and padding.

## Group reduction and capped procedure

On the protected group, phi(X,Y)=(beta X,Y)=[lambda](X,Y), with matching cube
roots in `argument.md`/`certificate.json`. The six units
{+/-1,+/-lambda,+/-lambda^2} act freely outside the identity. Class count:

`N = 19298681539552699237261830834781317975472927379845817397100860523586360249056`.

Canonicalisation selects least X,beta X,beta^2 X, then even Y, retaining its unit
for labels. Infinity/equal-x branches are explicit; identity target returns x=0.
Every output is challenge-verified.

A fresh 2^20-entry table T_j=s_j P+t_j Q feeds 1024 interleaved affine walkers
with transition canon(C+T_h(C)). Partition bits are always canonical bits 64..83.
DPs require 56 low zero bits, deepening to 57 then 58 without changing partition.
Full scalar pairs use bounded rejection sampling; table creation/starts are charged,
not supplied by a free short-seed random oracle.

The core budget is exactly
`T = 149174874732404653653643502645236769884`, obtained as
ceil(sqrt((98859/100000)*(27/25)^2*N)). Stop new trails after T and drain only
in-flight trails. Each trail has at most 40*2^58 ticks; global main ticks are
bounded by K=T+1024*40*2^58. These are integer ceilings, not rounded-log bounds.

The last-eight ring identifies a closed short-cycle suffix: rewind to its least
CYCLE member, double, canonicalise, clear the ring. Ordinary steps cost one tick;
doubling TWO. No unrelated tail minimum is used. Merged trails have a deterministic
eventual DP unless capped, without doubled-minimum/basin independence.

Brent compares the full expanded point/ring state, aborting repeats, never choosing
escapes. It catches escape-to-escape loops. It uses one evaluation stream versus
Floyd's three evaluations per iteration. The hard trail cap ensures termination;
no divergent-series long-cycle assertion remains.

DP records store complete initial scalar pairs rather than streaming full labels.
Matching keys trigger exact scalar replay: ordinary updates are u(a+s_j,b+t_j),
doubling updates 2u(a,b), modulo r. Equal full endpoints give
x=(a1-a2)/(b2-b1) mod r, with denominator and target checks. Unequal endpoints
are spurious key matches, never answers; their replay is capped and charged.

## Captured toys: two distinct rho statistics

Both programs ran under host `python3 -B`, only at toy scale: **p=524347**, not
approximately 2^40; curve `y^2=x^3+13`, order **523969**, **87328** order-six
classes. The fresh pipeline checks point count, cube-root action, canonicalisation,
labels and recovered scalars with its toy group law/log table. No production
search code ran.

### Fresh DP/recovery pipeline: exponent coefficient

`toy_pipeline.py` uses 8 independent 256-entry tables, seeds 8007 through 8014,
2000 trials/table, DP bits=4, key bits=12, ring length 8, expanded-state Brent,
core budget 318, trail cap 640. Doubling costs two ticks; only the current trail
drains. Counted recoveries require core-budget overlap, excluding drain-only
collisions. Scalars are verified, not sentinel successes.

Captured `toy_pipeline_results.json`:

| seed | core-overlap recoveries / 2000 | exponent rho |
|---|---:|---:|
| 8007 | 821 | 0.915643 |
| 8008 | 858 | 0.970888 |
| 8009 | 809 | 0.898098 |
| 8010 | 825 | 0.921531 |
| 8011 | 866 | 0.983068 |
| 8012 | 855 | 0.966342 |
| 8013 | 867 | 0.984596 |
| 8014 | 810 | 0.899553 |

Total **6711/16000 = 0.4194375**. From
`-log(1-p)/(318*317/(2*87328))`, pooled exponent rho is **0.9421122297**;
table range **0.8980975046 to 0.9845961256**. There were **481 Brent aborts** and
**13038 spurious-key comparisons**. Descriptive Wilson interval: approximately
**0.41181 to 0.42710**, not production confidence or 16000 independent table draws.
Sub-0.90 tables are disclosed, not selected away.

### Legacy repeat-time statistic

Rerunning `toy_experiment.py` reproduced the complete `toy_results.json`.
For 8 tables times 2000 starts, pooled inverse-mean-repeat-time quality is
**1.0514**, mean-of-table values **1.1085**, worst table **0.7003**. Per table:
1.2824, 1.3832, 0.9766, 1.4208, 0.9232, 0.9442, 1.2377, 0.7003.
Pooling uses mean TIMES, not averaged inverted ratios. This statistic,
sqrt(pi*N/2)/E[first-repeat-time], is **NOT** the birthday-exponent coefficient;
1.0514 is never substituted in the success formula.

## Conditional success and named margin

H_rate extrapolates the complete production process, averaged over fresh
coins/tables and uniform challenges: rho>=0.90 in
`1-exp(-rho*T*(T-1)/(2*N))`, including wasted cycle/escape/restart work within T.
It is not a per-table guarantee. H_count separately assumes a Poisson upper tail
for un-deepened DP counts; both are explicit where used in the argument.

Multiply the modeled collision term by (1-2^-20), then subtract exact overflow
and universal-fingerprint union allowances. Alternating-series rational bounds
give **success >0.4048 at rho=0.90**, versus claimed **0.40**: named absolute
margin **>0.0048**, not proof confidence. At rho=0.943 the floor exceeds 0.4193;
at rho=1 it exceeds 0.4381. Break-even for rules floor 0.39 is strictly between
**0.857 and 0.858**, below **0.86**. Arbitrary failure of H_rate invalidates the
success claim; local PASS flags do not settle this scientific uncertainty.

## Hard storage and exact-rational tail

The pool is **2^73 nodes x 192 bytes = 3*2^79 bytes**. Adding table, worker states,
code and routing scratch share a hard 2^30-byte auxiliary cap; peak **<2^81 bytes**.
A red-black tree has height <=148 with actual 96-bit byte addresses; loads and
transfers are charged at each visited level.

Record layout: 24-byte fingerprint, two 32-byte initial scalars, 8-byte trail
length, three 12-byte pointers, 1-byte color, 1-byte next-DP-bit flags, 58 bytes
metadata/padding. Its independent 192-bit binary Toeplitz universal hash covers
the full coordinate, not low 192 bits whose entropy is depleted by 56 DP bits.

On fullness at d=56 or 57, filter on stored next-DP bits and rebuild in place
at increased depth, with no replay. Fullness at d=58 fails. DP events and starts
are capped separately by 2^74; key-match replay pairs by 2^24. Every cap has a
failure branch. Up to three pool scans and two rebuilds fit the deepening row.
These memory/work caps do not depend on the Poisson model.

H_count bounds the DP-count upper tail by Pois(mu), mu=K*2^-56. A virtual split
into 2^65 subcounts, mean below 57 each, makes total >=2^73 imply some subcount
>=256. Dropping e^-57 gives the conservative exact rational
`2^65 * 57^256 * 257 / (256! * 200) < 2^-120`.
The certificate checks the reduced fraction. First overflow earns zero success,
even if deepening works. At most 2^74 distinct inputs give fingerprint pair union
bound <2^-45. Capped spurious reconstruction is charged.

## Complete work bound and uncertainty reserve

Ordinary extras: **925 scalar-64-bit instruction equivalents**, including eight
full 256-bit ring comparisons, Brent comparison/control and full checkpoint copy
even when unnecessary. Set E0=950 and charge **1140=1.2*950**
per tick beyond one reference affine operation; doubling debits two ticks.
The reference addition's six internal multiplications are not double-counted;
extra GLV products use the protected 179-instruction field-multiply calibration.

Separate rows: DP events, 2^18 instructions each for hashing and up to 148
multiword-address tree levels; starts, 4096 equivalents each; deepening, 2^96;
replay/labels, 2^110; exceptions, 2^40; table creation, 2^32; checks/sampling/output,
2^20. **Every row is multiplied by 1.2**. `verify_certificate.py` includes full
256-bit modular coefficient multiplication/reduction in replay.

Exact raw repriced **log2=127.9078900077**; an explicit **8% whole-ledger reserve**
gives **log2 W=128.0189213201**. These are analytic bounds, not measured work.
Submitted **128.06** is **0.0410786799 bit looser** than W, within the deliberate
**0.005–0.06-bit** round-up band, and **0.1521099923 bit above raw 1.2x repricing**,
exceeding **0.10-bit** required headroom. The distinct margins are not double-counted.

Old E=650 to 750 costs log2(1.75/1.65)=**0.0848888976 bit**, not slack absorbed
by nine instructions. The corrected old product is not an independent trap
probability. The audit recomputes integers, reduced fractions and directed bounds.

## Reproduction and validation

Package files: `submission.json`, `argument.md`, `certificate.json`,
`verify_certificate.py`, both toy programs/captured JSONs, `manifest.json`, this
note and empty `dependencies.json`. Manifest hashes bind captured artifacts.
From the package directory:

```text
python3 -B toy_experiment.py
python3 -B toy_pipeline.py
python3 -B verify_certificate.py
```

Toys are finite runs. The verifier compares exact arithmetic to captured values:
**65 named checks**, including full class integer, minimal integer square-root
budget, 1.2 repricing, reduced score, layout, Poisson rational, fingerprint union
bound and directed success bounds. Inequalities are integer/rational; decimal
logs are annotations only.

## Limitations

Neither toy certifies production mixing, 56-bit DP behavior, Poisson tails or
instruction prices: group, table size and DP density differ. Direct exponent
measurement improves on generic adding-walk citations, but transfer stays heuristic.
Candidate-authored checks certify arithmetic, not H_rate/H_count or instruction
estimates. This is no practical attack, hardware timing or formal production
certification. No private candidate program, protected-rule change, external advice
or inaccessible historical dependency was used. Caps stay hard if heuristics fail;
success need not. Remote review differs from local arithmetic/intake acceptance.
