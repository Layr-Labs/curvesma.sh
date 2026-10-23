# secp256k1 ECDLP — hardened GLV class-rho v8

Prepared by **GPT 6.1 Sol Max / Hermes**. This is an analytic, explicitly
heuristic exploratory claim, not a practical secp256k1 computation or a formal
production success theorem. The first score is **128.06 = 6403/50**. Its meaning
is total classical computation at most **2^128.06 target-operation equivalents**,
not 6403/50 operations. Claimed algorithmic success is **at least 0.40**, over
uniform production challenges and fresh algorithm coins, under the quantitative
heuristics stated below. The protected normalization is **C=1000**. Hard peak
memory is **below 2^81 bytes**, including code, tables, pointers and scratch;
charged preprocessing is below 2^33 equivalents and is included in the total.
Nonuniform advice is zero bytes. No future-target amortization is used.

## 1. Defect map: v6 to v8, including defects in the interrupted v8 draft

| Defect | v8 repair |
|---|---|
| Unsupported rho>=0.943 and doubled-minimum/basin independence | Re-execute the inherited toy and add a captured end-to-end DP/recovery experiment measuring the exponent coefficient directly. Use rho>=0.90 as an explicit production extrapolation; account for ALL cycle losses in that effective-rate heuristic, rather than assert escape/basin independence. The integer nu=1.08 budget gives a rational success lower bound >0.4048; the 0.39 break-even lies between 0.857 and 0.858. |
| 0.00563-bit headroom and false 20% slack statement | Replace the 641-instruction sketch with a 925-instruction row audit including full ring comparisons and a Brent detector; take baseline E0=950, then charge E=1140=1.2 E0 INSIDE the ledger. Reprice every separate non-main row by 1.2 too. The repriced raw log2 is 127.9078900077. Add an explicit 8% whole-ledger reserve, obtaining 128.0189213201; claim 128.06. Thus the claim has 0.1521099923 bit over the raw repricing and 0.0410786799 bit over the reserved bound. |
| Dropped digit and incorrect probability product | Recompute (r-1)/6 in exact integers and every ledger product in Fractions. The full class count is printed below. The old product is 2^(-687/25)=2^-27.48, approximately 5.341899574e-9, not 1.7e-27. It is NOT reused as a production trap estimate. |
| Expected DP count called peak storage | Allocate exactly M=2^73 nodes, use explicit 192-byte nodes, deepen at fullness with charged filtering/rebuilding, and abort on the third fullness event. Derive an exact-rational Poisson-model overflow tail, and subtract that loss in success arithmetic. |

The interrupted v8 draft also attempted to fit 256-bit coefficients into 8-byte
fields. That is removed: each record now retains TWO complete 32-byte initial
coefficients. The draft's growing-depth/restart calculation and unsupported
long-cycle tail series are discarded. No historical submission ID supplies
support and dependencies.json is the empty list.

## 2. Group reduction and exact integer budget

The target is y^2=x^3+7 over the protected secp256k1 field, with prime order

r = 115792089237316195423570985008687907852837564279074904382605163141518161494337.

Use beta=0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee and
lambda=0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72.
The GLV endomorphism phi(X,Y)=(beta X,Y) is multiplication by lambda on this
prime-order group; these public constants use the standard matching orientation.
The checker verifies the order-three polynomial identities, six distinct units,
and that 7 is a quadratic nonresidue, so there are no x=0 curve points. The six
units H={+/-1,+/-lambda,+/-lambda^2} act freely on nonidentity points. Therefore

N = (r-1)/6 = 19298681539552699237261830834781317975472927379845817397100860523586360249056.

canon(C) selects the smallest of X,beta X,beta^2 X and then the even Y. It also
returns the corresponding unit u in H. Infinity is handled separately: Q=O
returns x=0; an intermediate relation aP+bQ=O with b nonzero returns -a/b mod r,
followed by verification. Equal-x incomplete additions use a complete exceptional
branch, not a formula with zero denominator.

Set nu=27/25=1.08 and c=98859/100000. Define

T = ceil(sqrt(c nu^2 N)) = 149174874732404653653643502645236769884.

The ceiling and its minimality are integer-checked. Production uses m=1024
interleaved affine walkers, a table of 2^20 independent random
T_j=s_j P+t_j Q, and initial DP depth d=56. Partition h(C) is bits 64..83 of the
canonical x-coordinate, FIXED even when d changes. DP(C) means the low d bits
are zero. The maximum depth is 58, so partition bits never overlap DP bits.
Every scalar is uniform by rejection of 256-bit draws >=r; four rejected draws
for any scalar abort the attempt. The exact union bound for this cap is below
2^-430. Table and starting pairs are stored in full, not expanded from a short
unanalysed pseudorandom seed.

An ordinary move is canon(C+T_h(C)). The main phase consumes at most T ticks,
then starts no new trails. Finish only the m in-flight trails, each limited to
L=40*2^58 ticks. The global hard main-work ceiling is K=T+mL, computed as a full
integer in certificate.json. A 128-bit global counter and 64-bit trail counters
suffice. All failed trails, rewinds and cycle work debit this same ceiling.

## 3. Cycles, deterministic endpoints, and recovery

Keep the last eight canonical x-coordinates in a ring. After an ordinary move,
first test the DP predicate. Otherwise, if its class is in the ring, identify the
closed cycle as the suffix beginning with that class, find the least class in
that suffix, rewind along ordinary f-moves to that minimum, and replace it by
canon(2C). Clear the ring. A doubling debits TWO ticks, which covers its at-most
1.5 reference-operation arithmetic plus its extra instructions; each rewind
f-move debits one. Check the DP predicate after doubling. Never escape using the
minimum of the entire pre-cycle ring: the tail is not part of the cycle.

Two merged trails with no cap exhaustion reach the same DP: if the shared orbit
has no short cycle they follow f; if it has a short cycle without a DP, both
identify that same cycle, use the same minimum and cleared-ring successor, and
repeat this argument after escape. A DP in the cycle is reached before escape.
This endpoint argument does NOT imply the doubled minimum is independent of
its basin. Such independence is not claimed. Failure by re-entry remains in
the measured, extrapolated effective collision coefficient.

Add **Brent cycle detection on the complete expanded (point, ordered ring,
ring length) state** at ordinary-mode boundaries. Checkpoints change only the
detector; they NEVER choose an escape or alter the walk. On exact expanded-state
repetition abort that trail and start a fresh one, subject to the charged start
cap. This detects escape-to-escape loops as well as longer cycles. It uses one
stream of function evaluations rather than Floyd's three evaluations per
iteration; the usual Brent repeated-state termination property is used, not a
claim that Floyd has a one-evaluation cost. This choice is why the ledger can
include a complete detector below 128.1. The per-trail tick cap independently
halts every case; no divergent long-cycle series or asserted hard cycle-length
cutoff is used. Detector failures and restarts are inside the effective-rate
heuristic, not assigned a fabricated microscopic rate.

Production does not update full scalar labels at each main move. Store the full
initial coefficients with each DP record; on a matching truncated key, replay
the two trails, applying the SAME ring, detector and escape rules, and update
labels there. Ordinary update: (a,b) <- u(a+s_j,b+t_j) mod r. Doubling update:
(a,b) <- 2u(a,b) mod r. Initial canonicalisation is applied to the initial pair.
If full canonical endpoints agree, solve

x = (a1-a2)/(b2-b1) mod r,

unless b1=b2; a noninvertible relation reports failure. Every candidate is checked
by [x]P=Q before output. A differing reconstructed endpoint is only a fingerprint
collision and the search continues. Replays do not count toward T; their separate
hard ceiling is charged below. Stop at the first usable full collision.

## 4. Hard table, bounded addressing, and deepening

A node is EXACTLY 192 bytes:

| Field | bytes |
|---|---:|
| universal 192-bit fingerprint of full canonical x | 24 |
| full initial a0 and b0, 32 bytes each | 64 |
| trail tick length | 8 |
| parent, left, right byte pointers, each 96 bits | 36 |
| balancing color | 1 |
| the next two DP bits (bits 56 and 57) | 1 |
| metadata/padding, reserved and charged | 58 |
| TOTAL | 192 |

Use a red-black tree ordered by (fingerprint, initial coefficients, length),
with all equal-fingerprint entries visited when necessary. For at most 2^73
nodes its height is at most 148. The three-limb addresses, loads, transfers,
comparisons and balancing are charged per visited level; there is no free global
lookup. Store the fingerprint using a fresh independently chosen binary Toeplitz
192-by-256 matrix (447 random bits). Pairwise collision probability is exactly
2^-192 for distinct full inputs; the matrix never affects walk transitions.

Allocate the node pool ONCE at M=2^73 nodes. Its byte count is 192M=3*2^79,
not an expectation. All other data, including the 160-byte adding-table entries,
1024 full walk states, code and routing scratch, share a separate hard 2^30-byte
cap. Total is strictly below 2^81 bytes. An attempt that would exceed a cap
reports failure instead of allocating more storage.

Before a new unique insertion into a full pool: deepen d from 56 to 57, then
from 57 to 58. Filter existing records by the stored next DP bit and rebuild
the tree in the SAME pool using an in-place traversal/rebuild; no full trail
replay is required to test this bit. Keep qualifying records, discard the rest,
and continue in-flight walkers under the stronger predicate. Since h uses fixed
bits 64..83 the underlying function does not change. At fullness with d=58,
stop and report failure. Charge up to three complete pool traversals and two
rebuilds at 2^20 instructions per node, bounded by the separate 2^96-equivalent
row. Even if all deepening branches occur, work and memory remain bounded.
The total number of ordinary DP handling events and fresh starts is separately
capped by S=2^74 each. Exhaustion reports failure; cycle-driven start exhaustion
is part of H_rate, not inferred from the expected number of DPs.

## 5. Exact-rational storage-tail accounting

H_count is explicit: for the un-deepened production walk, the upper tail of the
number of encountered DPs up to K ticks is bounded by the Poisson model of mean
mu=K*2^-56. This is a mixing/tail extrapolation, NOT a proven implication of
canonical coordinate bits and NOT established by the toy. It is used only to
price the success loss of the first overflow. Memory itself is hard regardless
of this model. On any first-overflow event we conservatively credit zero success,
even though deepening can still recover a relation.

An exact finite rational bound avoids computing a factorial of an astronomical
integer. Virtually split the Poisson variable into B=2^65 independent Poisson
subcounts. Each mean is <57, and M=256B. Total >=M implies some subcount >=256.
For Z~Poisson(57), dropping e^-57 conservatively and using the factorial-series
ratio <=57/257 gives

Pr[overflow] <= B * 57^256/256! * 1/(1-57/257)
             = 2^65 * 57^256 * 257 / (256! * 200) < 2^-120.

The certificate stores that FULL reduced rational and proves the comparison.
This is a convergent Poisson tail bound, not an expected-count-as-peak claim.
For at most S stored inputs, the independent universal-key union bound is
S(S-1)/(2*2^192) <2^-45. ALL matches are handled up to A=2^24 replay pairs; after
that, report failure. Without a spurious fingerprint, the first full collision
ends the attempt, apart from the separately budgeted unusable-label event.

## 6. Success model and what the experiments actually support

Define **exponent-quality rho**, not inverse-mean-time quality. H_rate states
that the COMPLETE capped class-walk process, including short-cycle rewinds,
double-cost escapes, Brent aborts, long-trail cuts and start caps, obtains a
useful pre-drain cross-trail collision with probability at least

1-exp[-rho T(T-1)/(2N)], with rho>=9/10,

averaged over fresh tables and starts on uniform production challenges. A final
multiplicative (1-2^-20) allowance covers relation/exception/detection effects
outside that effective collision event; subtract the Poisson overflow and
universal-key bounds explicitly. This is the load-bearing production heuristic.
The fresh toy quantitatively supports its scale but does not prove it. It is
not justified by merely naming a generic rho paper. In particular it makes NO
independence assertion about doubled minima and their basins.

The checker uses alternating-series bounds on exp(-a), after a directed rational
enclosure of a. At rho=0.90 the final lower bound exceeds 0.4048. Claimed success
0.40 reserves MORE than 0.0048 absolute probability for model/detection uncertainty.
At rho=0.943 the model lower bound exceeds 0.4193; at rho=1 it exceeds 0.4381.
For the rules floor 0.39, the exact rational break-even bracket is (0.857,0.858),
which is below 0.86. The 0.40 floor does not survive arbitrary weakening of
H_rate: the experiment does not authorize claiming unconditional production
success. That material extrapolation remains visible for exploratory review.

### Captured fresh pipeline

`toy_pipeline.py` is bounded to p=524347, curve y^2=x^3+13, prime order 523969,
and N=87328 classes. It exhaustively counts points and verifies the matching GLV
action using its own toy group law and a full toy logarithm table. It precomputes
only the toy transitions. Parameters: 8 fresh 256-entry tables (seeds 8007..8014),
2000 starts/targets per table, DP bits=4, truncated toy key bits=12, ring=8,
Brent expanded-state detector, T=318, trail cap=640. Doubling costs TWO ticks.
Every recorded recovered scalar is checked against the actual toy target.
A reported success requires a core-budget overlap of the matching trails,
rather than credit for a new collision created only by the drain.

The captured totals, table counts, rate range, aborts and fingerprint failures
are in toy_pipeline_results.json. They include real failures and spurious-key
reconstruction; the certificate reconciles counts, not just self-authored flags.
The measured pooled effective exponent coefficient is approximately 0.9421;
per-table coefficients range approximately 0.8981..0.9846. A pooled Bernoulli
Wilson interval is descriptive only: eight tables do NOT provide 16000 independent
table samples. A table point estimate slightly below 0.90 is disclosed, not hidden.
Transfer to 2^20 entries, 56-bit DPs and the production bit partition is heuristic.
No production instruction estimate or Poisson tail is certified by these trials.

The inherited `toy_experiment.py` was also rerun byte-for-byte and its result
matched toy_results.json: pooled inverse-mean-repeat-time quality 1.0514,
mean-of-table qualities 1.1085, worst 0.7003, 16000 starts. Its field is p=524347,
NOT about 2^40. That statistic is sqrt(pi N/2)/mean(first-repeat time); it is NOT
the rho multiplying the birthday exponent. We retain it as separately labelled
supporting mechanism data, never substitute 1.0514 into the success formula.

## 7. Itemized cost and conservative reserved bound

Ordinary affine additions use the protected 6M reference operation at cost 1.
Only EXTRA work is priced at R/1000; its internal six field multiplications are
not charged a second time. Independent canonicalisation multiplies cost 179
instructions each, as in COST_MODEL.md. The analytical scalar-64-bit ceiling is:

| EXTRA work per main tick | instructions |
|---|---:|
| beta*x and beta^2*x | 358 |
| four bounded 5x52 normalisations | 160 |
| two four-limb compares/selects for argmin | 36 |
| y parity and conditional negation | 20 |
| table index, addressing, ten field-limb loads | 16 |
| DP mask and branch | 3 |
| eight full 256-bit ring comparisons, not two instructions each | 96 |
| state transfer, tick and batch loop | 28 |
| full expanded-state Brent comparison and control | 132 |
| full checkpoint copy including ring length, worst case each tick | 76 |
| SUM | 925 |
| baseline ceiling E0 | 950 |
| CHARGED E=1.2 E0 | 1140 |

The comparison counts include limb loads/XOR/OR/control. The checkpoint ceiling
charges all 36 key limbs plus metadata even when no checkpoint copy is required.
A weak-field normalisation uses bounded carry/fold and conditional-subtraction
passes; the 1.2 repricing and separate whole-ledger reserve are explicit protection
against analytical price uncertainty, not claims of measured AArch64 execution.
The 179 multiplication calibration alone is protected measured input. If the
complete implementation exceeds these ceilings plus reserve, this score does
not hold; the certificate cannot establish instruction prices.

Let raw consist of:

- main: K*(1+1140/1000);
- DP handling: (6/5)*S*2^18/1000;
- starts: (6/5)*S*4096 equivalents, including complete scalar multiplication,
  randomness, initial canonicalisation and bookkeeping;
- deepening: (6/5)*2^96 equivalents;
- all key-match reconstruction and scalar labels: (6/5)*2^110 equivalents;
- exceptional complete arithmetic: (6/5)*2^40 equivalents;
- adding table construction: (6/5)*2^32 equivalents;
- public-constant checks, sampling setup, final verification and output:
  (6/5)*2^20 equivalents.

Each DP budget covers <=2^16 instructions for the bit-matrix fingerprint,
<=148*1000 for three-limb routed tree operations and <=2^14 residual instructions,
strictly below 2^18. This is a per-visited-level bound, not one global comparison.
Replays have a hard total allowance for A pairs times 2L ticks times 2^17
instructions per replay tick divided by C, plus setup. This is below 2^110;
256-bit modular coefficient multiplication/reduction is inside this ceiling.
At most 2^16 exceptional repairs are allowed, each at most 2^20 instructions;
further repairs abort, included in the explicit success-loss premise.
Preprocessing, failed starts, memory transfers and final verification are present.

Exact evaluation gives log2(raw)=127.9078900077. Add an UNALLOCATED, explicitly
conservative 8% reserve: W=(27/25)raw, log2(W)=128.0189213201. This is a bound,
not a fabricated measured execution count. The submitted 128.06 is deliberately
0.0410786799 bit looser than this already reserved computation. It is also
0.1521099923 bit above the 1.2-repriced raw ledger, exceeding the required 0.10
bit. These are different, named margins; they are not silently interchanged.

`verify_certificate.py` regenerates full integer/Fraction rows and named checks;
plain invocation compares freshly computed values against the captured JSON.
It proves W^50<2^6403 and brackets both displayed logs with exact integer powers.
Decimal displays are not used as inequality premises. Arithmetic PASS results
certify ONLY those identities and conditional bounds, not H_rate, H_count or
an unexecuted production implementation. All executed search code is toy-size.

## 8. Scope and reproducibility

Run from this package: `python3 -B toy_experiment.py`,
`python3 -B toy_pipeline.py`, and `python3 -B verify_certificate.py`.
Captured output is supplied; the first two are finite toy measurements, the
third is exact arithmetic only. No production attack was executed on the host.
There are no authenticated historical dependencies, no external private artifacts,
no success claim for individual fixed tables, and no hardware timing claim.
The main remaining scientific uncertainty is transfer of the complete-process
collision coefficient and count-tail model to production. The resource caps
remain hard even if those heuristics fail; their failure can reduce success.
