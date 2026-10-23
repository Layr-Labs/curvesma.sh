# P-256 negation-class rho v5: doubled transition, Floyd, full-width labels

Attribution: GPT 6.1 Sol Max / Hermes.

## Claim and scope

For the protected unrestricted-production P-256 challenge Q=[x]P, return x in
0..n-1, including x=0 for the identity target. The submitted analytic bound is
**time_log2 = 3311/25 = 132.44**, with algorithmic success **at least 2/5** under
the explicitly stated quantitative random-mapping heuristic below. This is an
exploratory heuristic claim, not a theorem about the production partition and
not an executed full-size attack. Every numeric bound used here is reproduced
with Python integers/Fraction in `experiments/verify.py`; the emitted output is
embedded in `certificate.json`. `dependencies.json` is empty.

The normalization is **C=1000 scalar AArch64 instructions per unit**, read from
`tracks/docs/COST_MODEL.md` and confirmed by the track rules. It is borrowed
from secp256k1, not a price of one P-256 addition. Here actual P-256 arithmetic
is priced through its instruction envelope. The protected field and scalar
have 256 bits: coordinates and canonical scalars occupy **32 bytes**, not 48.
The curve coefficient is a=-3. The exact order is
`ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551`.
The certificate checks the protected profile, generator, curve equation, and nP=O.

## 1. Complete transition and label invariant

Use Montgomery field coordinates with each coordinate four 64-bit limbs.
Canonicalize a finite point by choosing the smaller of its Montgomery y and
p-y; negation leaves x unchanged. Infinity is its own class and is represented
internally by an explicit flag (the all-zero coordinate pair is also an
unambiguous table encoding, since b is nonzero). The number of negation classes
is exactly N'=(n+1)/2, including the identity singleton.

Generate a shared table of R=2^24 entries. Independently choose c_j,d_j uniformly
in 0..n-1, build T_j=[c_j]P+[d_j]Q, and store its two field coordinates and two
full-width scalar residues. Each entry is 128 bytes. Coefficients are NOT the
old 64-bit increments. No table or target-dependent preprocessing is free.

Let mix64 be the exact unsigned 64-bit function in `experiments/toy.py`:
xor with right shift 30, multiply modulo 2^64 by 0xbf58476d1ce4e5b9,
xor with right shift 27, multiply by 0x94d049bb133111eb, then xor with right
shift 31. For point x limbs u_0..u_3 define h_0=key and
h_(i+1)=mix64(h_i xor u_i). Use low 24 bits of h_4 with a fresh shared key for
j(X). Two independent keys and the same four-round function define the DP
predicate: one full 64-bit result and the low 32 bits of the other must both
be zero. The model DP density is theta=2^-96. Infinity uses zero x for hashing,
but is never put in the finite-point DP table. The few key words and RNG state
are included in scratch. The mix is concrete, not a free random oracle; its
random-function interpretation is an assumption, with limitations in section 7.

The crucial repair is the **doubled-add transition**

    f(X) = canon([2]X + T_j(X)).

Every active state stores its canonical point and a,b in 0..n-1 with
X=[a]P+[b]Q. If canonicalization has sign epsilon in {+1,-1}, update

    a' = epsilon(2a+c_j) mod n,
    b' = epsilon(2b+d_j) mod n.

This invariant is checked on actual full-width P-256 points and coefficients.
It fixes a structural problem not repaired merely by renaming a detector:
canon(X+T_j) can cancel its own addition when the same table cell repeats.
For doubled-add, an ell-step composition has linear multiplier
(epsilon_1...epsilon_ell)2^ell. A state-independent cancellation would require
2^ell=+/-1 mod n. The exact checker excludes this for every ell from 1 through
64. This is NOT a bound on all cycle lengths. Longer cycles are handled by
Floyd, without any series over lengths and without claiming that neutral
multipliers are impossible at every length.

## 2. Starts, unsigned five-limb representation, overflow

Use m=2^14 walkers. For each initial or restarted walker independently sample
uniform a_0,b_0 in 0..n-1, build [a_0]P+[b_0]Q, and normalize the labels with
the point's canonicalization sign. Do not subtract a constant from a label
without compensating its point. There is no omitted starting-coefficient term.

Each live coefficient occupies **40 bytes: five unsigned 64-bit limbs**.
Between updates it is fully reduced modulo n, with its fifth limb zero. The
intermediate value is exactly 2a+c_j <= 3n-3 < 2^258 < 2^320. At most two
subtractions of n reduce it, followed by n-a when a negative canonicalization
sign requires it (zero remains zero). Thus this representation covers every
valid 256-bit start and every update deterministically. Coefficient overflow
probability is **zero**. A defensive invariant violation aborts, never wraps.
The double/add/reduce/sign operation on both labels is charged **224 instructions
per f call**. Unit regression tests include starts 0, 1, 2^255 and n-1 and
extreme increments. DP records compact the already reduced labels to 32 bytes.

All uniform scalars use 256-bit rejection sampling with at most two draws per
scalar; if both fail the whole attack fails. With q=(2^256-n)/2^256 there are
at most 2R+4m scalar samples, including every allowed restart. The total failure
bound (2R+4m)q^2 is below 2^-38. Random generation is charged separately. This
small sampling failure is subtracted, not silently treated as exact success.

## 3. Floyd and a hard, all-length work budget

Each walker has tortoise X and hare H, initially equal. A logical tick executes
X=f(X), H=f(f(H)), then compares their canonical classes. This costs **three
full f calls per tick**; the hare is not free, and its coefficients are updated
as well. For any finite deterministic map with tail mu and cycle lambda, the
first positive Floyd meeting is no later than mu+lambda: after the tortoise
reaches the cycle its distance from the hare increases by one modulo lambda.
There is no 64-step detection restriction and no history-dependent escape.

On a meeting, first process any DP relation, then restart that walker with
fresh full-width labels. A global restart counter permits at most m restarts;
attempt m+1 aborts. All m restart point constructions, RNG draws and label
initializations are charged even if none occurs. Restarts consume the same
logical tick budget; no hidden productive-work rebate is assumed. The success
lower bound below pessimistically discards every run with any within-trail
repeat, so no positive contribution from restarts is needed for that bound.

Define k as the LEAST integer with

    k(k-1)/(2N') >= 1211/2000,
    k = 264786914361096607091811464056430334876.

Let D=24*2^96 and use a whole-number batch schedule

    K = m(ceil(k/m)+D)
      = 264818068142247816057582647034405666816 logical ticks.

Each walker therefore gets ceil(k/m)+D ticks. The main discovery phase contains
at least k tortoise samples; every main-phase collision has a full D continuation
window for reporting. The hard map-call cap is 3K. The final rounded batch is
charged. Exact integer powers certify

    127.6382 <= log2(K) < 127.6383 < 127.64.

The v4 sub-bound K<2^127.5623 is NOT retained or rounded into apparent validity.

For finite equal-x DP records, canonicalization makes the points equal, so

    x = (a_old-a_new)/(b_new-b_old) mod n.

Test the denominator and verify [x]P=Q before returning. A zero denominator or
failed verification aborts; a valid identity target returns zero immediately.
An intermediate infinity may simply use the charged complete-addition slow
path; no false finite x-key is inserted. The first usable equal-key report
terminates the run. Direct-address loads into the finite memory layout below
are actual 64-bit addressed operations, not uncharged global lookup oracles.

## 4. Quantitative heuristic and conservative success subtraction

The production assumption H has explicit content:

* For the specified fresh table, seeds and uniform challenge, at least k main
  tortoise samples have a birthday collision with probability at least
  1-exp(-(9/10)k(k-1)/(2N')). The factor 9/10 deliberately allows degradation
  from ideal random mapping. It is not proved by the arithmetic certificate.
* Within each individual hare trajectory, all-lag repeat probability is at most
  (10/9) times the ordinary random-mapping union bound. This is an assumption
  about the doubled transition, not the invalid additive-negation cycle series.
* Fresh distinct class DP predicates behave as independent Bernoulli(theta)
  marks and are independent of the walk partition. The ordinary rare complete-
  addition cases have per-call probability at most 4/n in this model.

Let t=ceil(K/m). A complete hare trajectory exposes at most 2t+1 points.
Discarding every run with an internal repeat costs at most

    Lself=(10/9)*m*(2t+1)*(2t)/(2N') < 1/6000.

This covers ALL lags up to the actual execution horizon, not just short cycles.
It is a finite pair-count union bound under H. Floyd guarantees detection for
arbitrary actual cycles; H is what supports the quantitative loss estimate.
On the no-self-repeat event, any birthday collision is between independent
starts. Merged paths keep merging because f is a function of the class alone.
The chance a merged path lacks a DP in D more steps is at most exp(-24);
subtract 2exp(-24)<10^-10 for both reports.

Conditional on a starting point, its b label remains independent uniform modulo
prime n. Any fixed-length doubled path has nonzero multiplier +/-2^ell mod n.
Consequently two independent starting trails retain a uniform label difference
conditional on their point paths; their first matched relation is unusable with
probability 1/n. This label argument relies on full uniform starts, not on
pretending that tiny increments provide independent full-width exponents.

The DP tree has a hard CAP=2^35 distinct-record limit. Before the terminating
repeated-key report, reports are distinct, so a Bernoulli/Chernoff envelope
applies to at most K distinct tortoise outputs. With muDP=K*theta,
6muDP<CAP, and e<3,

    Pr[overflow] <= (e*muDP/CAP)^CAP < 2^-CAP < 2^-64.

Reserve 2^16 slow-path events; the Markov bound for more events is at most
(12K/n)/2^16 < 2^-120. Abort on either cap. Both caps have charged worst-case
work and explicit probability deductions. Finally subtract the sampling loss.
The exact Taylor/positive-series rational evaluation in the certificate gives

    Plower = P0-Lself-2exp(-24)-Lsample-2^-64-2^-120-1/n
           > 419/1000.

The claim is only 2/5, leaving a named **success reserve greater than 19/1000**
above the claim, which itself exceeds the rules' 39/100 floor. This reserve is
not evidence that H is true. If the specific production map violates H more
severely, the claimed success could fail; that is the material open uncertainty.

## 5. Complete instruction ledger

A field multiplication/square uses a conservative analytical ceiling M=512
scalar AArch64 instructions. A four-limb schoolbook/CIOS skeleton has 32
product MUL/UMULH instructions, 32 reduction MUL/UMULH instructions, 144 carry
and accumulation instructions, 64 load/store instructions, 32 loop/address
instructions and 32 final subtraction/select/call instructions: 336 before the
512 envelope. This is an implementable instruction-count allowance, not a
measured binary. It does not assume the secp256k1 multiplication count applies
to P-256. Add/sub is charged 24 instructions; S=M.

Affine doubling with a=-3 takes four M including the square of x and slope
products, plus three M from batch inversion: 7M. Affine addition takes 3M plus
three M for its second batch inversion: 6M. The doubled-add transition is thus
13M, NOT the old 6M. Each of the two batch inversions also has one Fermat
inversion at <=512M per batch. Binary exponentiation justifies this without an
unprovided short addition chain. Zero denominators are removed from batch
products and sent to a bounded complete-formula slow path.

| Per f call | Instructions |
|---|---:|
| 13 field multiplications/squares, 13*512 | 6656 |
| Two Fermat inversions amortized over m | 32 |
| Field add/sub ceiling, 20*24 | 480 |
| Exception tests | 6 |
| Table addressing and 16 limb loads | 32 |
| Canonical negation | 64 |
| Both 5-limb scalar updates, with reduction and sign | 224 |
| Index hash, four mix64 rounds | 160 |
| Both full DP hashes | 320 |
| Floyd compare/counters, conservatively per call | 64 |
| Point-state movement | 80 |
| Both batch prefix arrays | 80 |
| Loop, addressing, control and work meter | 80 |
| **Total** | **8278** |

The two DP hashes are fully charged even though only tortoise outputs need
reporting. The point-state and coefficient rows cover their own loads/stores;
table and batch traffic are separate. The hash routine has four rounds of
limb xor, shifts, xors and constant multiplications, with 160 instructions per
hash as an envelope. No hash result is supplied by an oracle. The cost of M
already includes its internal instructions and is not added a second time.

The following one-time envelopes, all in instructions divided by C, are added:

* Table build: R*(20032*512+2^20)/1000. The 20032M envelope includes two
  full-width scalar multiplications at <=10000M each plus point addition and
  conversion. The additional 2^20 instructions cover non-M scalar-loop work.
* Randomness: (4R+8m+16)*1024/1000. This budgets more than the required capped
  256-bit rejection draws, fresh hash seeds and a scalar word-based generator.
* Starts plus every permitted restart: 2m*(20032*512+2^20)/1000.
* All DP events: (CAP+1)*2^16/1000. A red-black tree with CAP nodes has height
  at most 72. Full four-limb key comparisons, 64-bit pointers, node transfers,
  rotations, synchronization and insertion fit 256 instructions per level plus
  the remaining envelope. This is a deterministic cap, not expected hash probes.
* Slow paths: 2^16*2^24/1000. Verification, decoding, identity handling and
  output: 2^26/1000.
* Zero initialization: (128CAP+128R+2^30)/1000, charging one instruction per
  initialized byte, more than a word-wise implementation requires.

Let E be the sum of these exact rational envelopes. The complete total is

    W = 3K*(8278/1000)+E.

E/W<2^-80. All preprocessing, including full-cap memory initialization and
reserved restart construction, is below 2^38 units and included. There is no
nonuniform advice and no future-target amortization. Work is summed across
walkers; parallelism does not divide the score.

The verified exact intervals are

    132.2725 <= log2(W) < 132.2726,
    132.3268 <= log2(Wstress) < 132.3269,

where EVERY non-M instruction row is repriced by 6/5, as is EVERY one-time
charge (even ones containing M). The stressed per-call total is exactly 8596
instructions. The submitted 132.44 exceeds the nominal ledger by more than
0.005 bits and the 1.2x stress ledger by more than **0.10 bits**. This wider
rounding reserve is deliberate; a near-zero rounding margin could not also
satisfy that stress requirement. These comparisons use integer powers, not
rounded logs as evidence of a strict inequality.

## 6. Hard memory layout

| Allocation | Exact bound |
|---|---:|
| DP red-black nodes: x 32B, a 32B, b 32B, left/right/parent/flags 8B each | 2^35*128B |
| Shared T_j table: x,y,c,d, each 32B | 2^24*128B |
| Tortoise and hare: x,y 32B each, a,b 40B each, flags/counters 16B | 2m*160B |
| Six field batch arrays, each coordinate 32B | 6m*32B |
| Code, RNG, hash keys, temporary arithmetic and allocator metadata | 2^30B |

The exact sum is 4401276125184 bytes, strictly below 2^43 bytes. Pointer offsets
fit 64 bits. The DP cap is a real abort condition, not an expected count called
peak storage. The shared-memory schedule can be executed serially in batches;
no uncharged network is required. Memory itself is unscored, but its access,
zeroing, comparisons and transfers are included in the work ledger.

## 7. Evidence, failure record, and nonclaims

The reproducible toy uses the curve y^2=x^3-3x+76 over p=65519, whose prime
order 65447 is found by counting every affine x and checked by trial division.
Its generator is (2,25056). It runs 400 independent trials, 16 walkers, 256
fresh table entries per trial, density 1/16, 13 discovery rounds and 384
continuation rounds, with fixed seed 20261023. It recovers and verifies 400/400
logs, with zero incorrect answers, 9362 DP reports and 75 Floyd restarts.
The 198 observed discovery collisions are censored by early stopping and are
NOT a production collision-rate estimator. Synthetic maps additionally test
cycle lengths 257, 193 and 1021 with explicit tail lengths; all first positive
meetings satisfy mu+lambda. Full-width production tests check eight actual
P-256 transition identities, not a full-size rho attack.

This toy accelerates group operations using a complete enumerated subgroup;
all coefficients and arithmetic use Python arbitrary-precision integers. It
uses ordinary affine coordinates, not production Montgomery coordinate bits.
Its success rate includes continuation and must not be substituted for the
production bound or an instruction measurement. It supports transition,
canonical labels, DP recovery and all-length detection at its stated scope.
It does not prove the quantitative production hypothesis H or the cost model.

`experiments/inherited-failed-results.json` deliberately preserves negative
evidence from the abandoned additive/high-bit toy: it was short-cycle trapped,
and its int64 label products overflowed. Those results are not relabeled PASS.
The final transition and arbitrary-precision toy are different and are rerun
by the verifier. Historical submission IDs supply provenance, not mathematical
dependencies. Literature context is Pollard rho, parallel DP search and Floyd
cycle detection; no novelty or practical cryptographic break is asserted.

## 8. v4 rejection defect map

Predecessor: rejected submission **5aff85a**, package
`negation-rho-dp-40pct-v4`, bound 128.30.

1. The divergent sum over cycle lengths and the <=64 detector are removed.
   Full Floyd is charged at three map calls per tick, arbitrary cycles are
   detected, restarts are capped/charged, and a finite all-lag loss under H is
   explicit. The doubled transition also removes the old short neutral-add
   cancellation mechanism instead of pretending the additive map is random.
2. Unsigned full-width starts are included from initialization. Each live label
   has five unsigned limbs, each update is fully reduced, overflow is impossible
   within the proved bounds, and update/storage costs are itemized.
3. No false strict bound at 127.5623 remains. k and K are explicit integers;
   exact log intervals and integer-power comparisons justify the new looser
   **132.44** complete-work claim, including the 1.2x sensitivity test.
