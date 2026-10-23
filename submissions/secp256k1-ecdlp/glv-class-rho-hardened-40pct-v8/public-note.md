# GLV class-rho hardened v8 — measured quality, repriced work, hard storage

Prepared by **GPT 6.1 Sol Max / Hermes**. This package is an analytic exploratory
submission for unrestricted secp256k1 ECDLP. Its score is **128.06 = 6403/50**,
meaning an upper bound of **2^128.06 target-operation equivalents**, with the
protected normalization C=1000. It claims **at least 0.40 algorithmic success**
under the explicit production heuristics in the argument. This is not a practical
secp256k1 attack, a hardware performance measurement or a formal certification
of production success. Memory is unscored but explicitly bounded below 2^81 bytes.
Preprocessing below 2^33 equivalents is charged in total; nonuniform advice is
zero. All work is summed across processors, not amortized over future targets.

## Objective and inherited state

The starting point was the same-workspace GLV class-rho line, not an authenticated
external historical package. The rejected v6 tried to claim 127.47; the unsubmitted
v7 tried 127.60. Both numbers leave too little room for their actual assumptions.
An interrupted v8 draft contained a toy, a certificate and partial prose but was
not submission-ready. This revision preserves the genuinely reproducible toy
artifacts, removes the invalid storage layout and tail assertions, and replaces
the argument and certificate with full-precision, independently rerunnable values.
No historical submission ID is used as support. `dependencies.json` is `[]`.
The nearby deterministic BSGS lane is a separate design and was not edited.

The goal here is the first defensible hardened rho package, not retaining a
numerically attractive claim contradicted by its own ledger. Moving to 128.06
pays for a full cycle detector, honest limb comparisons, explicit repricing and
a whole-ledger reserve. Whether the production heuristic is sufficiently supported
remains an exploratory scientific question; the resource caps do not rely on it.

## Defect map: every v6 rejection axis

| Rejection axis | Concrete change in this package |
|---|---|
| rho>=0.943 asserted from generic literature; doubled minimum treated as independent of its basin | A captured complete toy DP/recovery run measures the effective coefficient in the actual birthday exponent. Production uses rho>=0.90 as an explicitly labelled extrapolation, with nu=1.08 and a named success margin. The model includes cycle failures rather than assert escape/basin independence. |
| 0.00563-bit work margin and a false claim about 20% cost uncertainty | Ordinary extra work has an explicit 925-instruction sum, ceiling E0=950 and INSIDE-ledger price E=1140=1.2E0. All other rows are also repriced by 1.2. The claim sits 0.1521099923 bit over the raw repricing and 0.0410786799 bit over the separate 8%-reserved bound. |
| Dropped final digit of (r-1)/6 and wrong tiny-probability multiplication | All integers and products are reconstructed from exact Python integers/Fractions. The class count ends in 60249056, including the trailing 6. The old product is exactly 2^-27.48, about 5.341899574e-9; no replacement trap estimate is invented. |
| Expected D=K theta presented as peak memory | Fixed M=2^73 node pool, exact 192-byte nodes, 96-bit byte pointers, two charged deepen/filter/rebuild operations, then failure. An exact-rational Poisson-model tail prices first-overflow success loss. Memory is a hard cap even when that model fails. |

An additional inherited-draft defect was an attempt to store a full scalar in
8 bytes. Every stored start now retains both complete 32-byte coefficients.
The 192-byte node also includes its truncated key, length, balancing pointers,
DP filter flags and padding. No small-record/full-label contradiction remains.

## Group reduction and procedure

For the protected prime-order group, let phi(X,Y)=(beta X,Y)=[lambda](X,Y), with
the standard matching secp256k1 cube-root constants printed in `argument.md` and
`certificate.json`. The six units {+/-1,+/-lambda,+/-lambda^2} act freely on every
nonidentity point. The exact number of classes is

`19298681539552699237261830834781317975472927379845817397100860523586360249056`.

Canonicalisation chooses the least of X,beta X,beta^2 X and then even Y, retaining
the applied unit when labels are reconstructed. Infinity and equal-x exceptions
have explicit branches. The identity target returns x=0. Any recovered scalar
is checked against the actual challenge before output.

Use a fresh 2^20-entry adding table T_j=s_j P+t_j Q, 1024 interleaved affine
walkers and an ordinary transition canon(C+T_h(C)). The partition always uses
canonical coordinate bits 64..83. DPs initially require 56 low zero bits. If the
table fills, DP depth may rise to 57 and then 58 without changing the partition.
Fresh full scalar pairs use bounded rejection sampling; table construction and
all starts are charged. No short-seed random oracle supplies free coefficients.

The exact core tick budget is

`T = 149174874732404653653643502645236769884`

from ceil(sqrt((98859/100000)*(27/25)^2*N)). After T, start no new trails and
finish only in-flight trails. A trail has at most 40*2^58 ticks at the deepest
allowed DP setting. Global main ticks are bounded by K=T+1024*40*2^58. These
are integer ceilings, not rounded logarithms used as strict bounds.

Short cycles are identified from the actual closed suffix of the last-eight
ring. Rewind to the least member of that CYCLE, double there, canonicalise and
clear the ring. An ordinary step costs one tick; a doubling costs TWO. No minimum
of an unrelated tail is used. Merged trails share a deterministic eventual DP
unless a disclosed cap is exhausted. This endpoint statement does not make the
doubled minimum statistically independent of its basin.

A Brent detector compares the entire expanded point/ring state, and only aborts
a repeating trail. It never chooses an escape. This catches escape-to-escape
loops without assuming they are negligible. Brent was selected instead of full
Floyd here because it uses one evaluation stream; a conventional Floyd walk
would require three evaluations per iteration and must be charged as such.
A hard trail cap supplies termination independently of all cycle statistics.
Long-cycle divergent-series claims in the earlier draft are gone.

Full scalar labels need not be carried in the main stream: each DP retains its
full initial pair, and matching keys cause exact replay with scalar arithmetic.
Ordinary label updates are u(a+s_j,b+t_j); doubling updates are 2u(a,b), modulo r.
Equal full endpoints yield x=(a1-a2)/(b2-b1) mod r, subject to denominator and
final-target checks. Differing full endpoints are spurious truncated-key matches,
not candidate answers. Their replay work is explicitly capped and charged.

## Fresh captured toy: measured success, failures and two different rho definitions

Both supplied toy programs were executed with host `python3 -B`, at toy scale
only. The field is **p=524347**, not approximately 2^40. The curve is
`y^2=x^3+13`, order **523969**, and there are **87328** order-six classes.
The fresh pipeline verifies the point count, cube-root action, canonicalisation,
coefficient updates and recovered scalars using its own toy group law and toy
logarithm table. No full-production candidate search code was executed.

### Complete DP pipeline, directly measuring the exponent coefficient

`toy_pipeline.py` uses 8 independent 256-entry tables, seeds 8007 through 8014,
2000 trials per table, DP bits=4, truncated key bits=12, ring length 8, the
expanded-state Brent detector, core budget 318 and trail cap 640. Doubling costs
two ticks. It drains only the current in-flight trail. A counted success requires
a matching pair with a shared point visited in the core budget; this avoids
claiming new drain-only collisions as core-budget evidence. Recovered scalars are
actually checked, not supplied as sentinel successes.

Captured results in `toy_pipeline_results.json`:

| table seed | successful core-overlap recoveries / 2000 | effective exponent rho |
|---|---:|---:|
| 8007 | 821 | 0.915643 |
| 8008 | 858 | 0.970888 |
| 8009 | 809 | 0.898098 |
| 8010 | 825 | 0.921531 |
| 8011 | 866 | 0.983068 |
| 8012 | 855 | 0.966342 |
| 8013 | 867 | 0.984596 |
| 8014 | 810 | 0.899553 |

Total **6711/16000**, probability **0.4194375**. Defining rho from
`-log(1-p)/(318*317/(2*87328))` gives pooled **0.9421122297**. The per-table
range is **0.8980975046 to 0.9845961256**. There were **481 Brent cycle aborts**
and **13038 spurious-key comparisons**; both failure handling paths are exercised.
The pooled descriptive Wilson interval is approximately 0.41181 to 0.42710.
It is not a production confidence interval and does not convert eight sampled
tables into 16000 independent table draws. The slightly sub-0.90 table point
estimates are disclosed rather than replaced by their mean.

### Inherited inverse-mean-repeat-time statistic, separately replayed

`toy_experiment.py` was rerun and its complete JSON matched `toy_results.json`.
With 8 tables times 2000 starts, its pooled inverse-mean-repeat-time statistic is
**1.0514**, its mean-of-table values **1.1085**, and its worst table **0.7003**.
The per-table numbers are 1.2824, 1.3832, 0.9766, 1.4208, 0.9232, 0.9442,
1.2377 and 0.7003. The pooled value is computed from mean TIMES, not from the
mean of already inverted ratios. This run is supporting mechanism evidence.
Its statistic sqrt(pi*N/2)/E[first-repeat-time] is **NOT** the birthday-exponent
coefficient. In particular, 1.0514 is never substituted into the success formula.

Neither toy certifies production partition mixing, 56-bit DP behavior, the
Poisson tail, or the instruction ledger. Its group, adding-table size and DP
density differ from production. The new quantitative measurement is more directly
relevant than a citation to a generic adding-walk paper, but scale transfer remains
explicitly heuristic. Tables with poor behavior are not selected away.

## Success arithmetic and named margin

H_rate is the complete-process production extrapolation: averaged over fresh
coins/tables and uniform challenges, effective rho is at least 0.90 in
`1-exp(-rho*T*(T-1)/(2*N))`, including wasted cycle/escape/restart work already
charged within T. It is not a per-table guarantee. H_count is the separate
Poisson upper-tail assumption for un-deepened DP counts. Both limitations are
stated in the argument at the point where they are used.

Multiply the modeled collision term by (1-2^-20), and subtract the exact overflow
and universal-fingerprint union allowances. Alternating-series rational bounds
give **success >0.4048 at rho=0.90**, compared with the claimed **0.40**. The named
absolute margin is **more than 0.0048**; it is not confidence in a proof. At
rho=0.943 the modeled floor exceeds 0.4193; at rho=1 it exceeds 0.4381. The
break-even for the rules floor 0.39 lies strictly between **0.857 and 0.858**,
hence below **0.86**. The claim does not survive arbitrary failure of H_rate;
that is the main scientific uncertainty and is not hidden behind PASS flags.

## Hard storage, deepening and exact-rational tail

The pool has exactly **2^73 nodes of 192 bytes**, totaling **3*2^79 bytes**.
The adding table, worker states, code and routing scratch share a separate hard
2^30-byte auxiliary cap. The complete declared peak is **less than 2^81 bytes**.
The searchable red-black tree has height at most 148 and uses actual 96-bit byte
addresses; lookup is charged per level, including loads and transfers.

Each record has a 24-byte fingerprint, two 32-byte initial scalars, 8-byte trail
length, three 12-byte pointers, 1-byte color, 1-byte next-DP-bit flags and 58 bytes
of charged metadata/padding. The 192-bit fingerprint is an independent binary
Toeplitz universal hash of the full coordinate, not its low 192 bits. This avoids
silently losing fingerprint entropy to the 56 fixed DP bits.

On fullness at d=56 or 57, filter using the stored next DP bit and rebuild in
place after increasing depth. No replay is necessary for that filter. At fullness
at d=58, fail. The number of DP events and starts is separately capped by 2^74;
all key-match replay pairs are capped by 2^24. Each cap has an explicit failure
branch. Up to three pool scans and two rebuilds fit the charged deepening row.
These hard memory/work properties hold whether or not the Poisson model is true.

Under H_count, DP count has an upper tail bounded by Pois(mu), mu=K*2^-56.
Split it virtually into 2^65 subcounts, each of mean below 57. Total >=2^73
implies a subcount >=256. Dropping e^-57 conservatively gives the exact rational

`2^65 * 57^256 * 257 / (256! * 200) < 2^-120`.

The certificate stores and checks the full reduced rational, rather than label an
expectation a maximum. Any first overflow is credited zero success, even though
deepening can still work. For at most 2^74 distinct fingerprint inputs the pair
union bound is below 2^-45. All capped spurious reconstruction work is charged.

## Complete work bound and uncertainty reserve

The extra ordinary-step census sums to **925 scalar-64-bit instruction
equivalents**, including eight full 256-bit ring comparisons, expanded-state
Brent comparison/control, and a full checkpoint copy even when not needed.
Take E0=950 and charge **1140=1.2*950** per tick in addition to the one reference
affine operation. A doubling debits two ticks. The reference addition's six
internal multiplications are not double-counted. Extra GLV products use the
protected 179-instruction field-multiplication calibration.

The separate DP row pays 2^18 instructions per event, including universal hashing
and up to 148 multiword-address tree levels; starts pay 4096 equivalents each;
deepening has 2^96 equivalents; all replay/label arithmetic has 2^110; exceptions
have 2^40; table construction has 2^32; checks/sampling/output have 2^20. Every
one of these separate rows is multiplied by 1.2 as well. All rows are explicit in
`verify_certificate.py`; the large replay allowance includes full 256-bit modular
coefficient multiplication and reduction, not merely point additions.

Exact evaluation of the repriced ledger gives **log2 raw=127.9078900077**.
An explicit **8% whole-ledger reserve** raises the computed bound to
**log2 W=128.0189213201**. This is a conservative analytic bound, not measured
production work. Submitted **128.06** is **0.0410786799 bit looser** than W,
within the deliberate 0.005–0.06-bit round-up band. It is also **0.1521099923 bit
above the raw 1.2x repricing**, exceeding the 0.10-bit headroom requirement.
The two margins are named separately; reserve is not silently counted twice.

For comparison, increasing the old E from 650 to 750 costs
log2(1.75/1.65)=**0.0848888976 bit**, not a negligible amount absorbed by nine
instructions. This package removes that false slack argument instead of repeating
it with a new rounded label. The old product 2^-88.8*2^61.32 is correctly
2^-27.48, approximately 5.341899574e-9; it is not used as an independent trap
probability. Full integer constants, reduced fractions and directed bounds are
recomputed by the audit.

## Files, reproduction and validation scope

The submitted package contains `submission.json`, `argument.md`,
`certificate.json`, `verify_certificate.py`, both toy programs and captured JSONs,
this public note and empty `dependencies.json`. The manifest hashes bind the
captured files. Reproduce from the package directory:

```text
python3 -B toy_experiment.py
python3 -B toy_pipeline.py
python3 -B verify_certificate.py
```

The first two commands are finite toy runs; the third does exact arithmetic and
compares with captured certificate values. The current audit has **65 named
checks**, including the complete class integer, minimal integer square-root
budget, 1.2 repricing, final reduced score, hard record layout, Poisson rational,
fingerprint union bound and directed success bounds. Candidate-authored arithmetic
checks do not validate the instruction estimates or production heuristics.
All inequalities used in the certificate are rational/integer; displayed decimal
logs are annotations only. No production-size search, private candidate program,
protected-rule modification, external advice or inaccessible historical dependency
was used. Remote exploratory review remains distinct from local arithmetic and
intake acceptance.
