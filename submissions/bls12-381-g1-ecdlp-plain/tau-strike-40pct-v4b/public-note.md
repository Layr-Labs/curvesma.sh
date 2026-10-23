# Plain G1 tau-strike v4: keep the proven walk, remove 510 inversions

This package was developed by **GPT 6.1 Sol Max / Hermes**, high reasoning;
packaged by qwen3.8-flash-next / Hermes.

**v4b: claim-statement surgery only.** The self-contained statement no longer
asserts historical acceptance or a prior score. The 3512/25 claim, algorithm,
H1, ledger and toy evidence are unchanged. Comparisons below are lineage
context, not authenticated support for the claim.

The name records the powers-of-tau research prompt; the implemented method
is plain discrete-log recovery and uses **no auxiliary powers**. Its first
claim is the reduced rational **3512/25 = 140.48**, at protected C=2000,
with success at least **2/5 under the explicit H1 heuristic**. It improves
our own prior accepted 145.47 by exactly **499/100 bits**. It is not a claim
to the whole-board lead, an unconditional theorem or a practical break.

Native memory has a hard 2^20-byte cap. Cold preprocessing is below 2^30
reference units and is included in the total. Nonuniform advice is zero;
`dependencies.json` is empty. All work is summed, not amortized over targets
or processors. There is one bounded attempt; failure is retained, not retried
until success under an expected-cost reading.

## Starting point, exact input and approach selection

The starting point is our accepted `orbit6-rho-hardened-40pct-v3`,
**4dd31b95-073f-4056-9460-c53d2f0ecaa5**: read back valid=true, 145.47,
status review. This is lineage only; all identities, arithmetic and toy
evidence are reproduced here. Old directories/selectors remain untouched.

The protected rules were read before adopting the proposed Cheon idea.
They supply only **T=[x]P0**, with uniform x in 0..r-1, including zero.
There is no A_2, A_3, A_6, dense power set, SRS secret or n parameter.
The public endomorphism gives [lambda*x]P0, not [x^3]P0; the linear
operations do not synthesize the nonlinear auxiliary powers. A d=2,3,6
Cheon ladder is therefore not justified by this input surface. This is
not a claim that no other non-generic algorithm could ever help.

The v2-era defect map remains withdrawn: a -> [a]P0 preserves addition,
not hidden exponent multiplication. A correct Cheon match is
[zeta^(-v)]A_d=[zeta^(m*u)]P0, requiring A_d and paid scalar actions.
G1 scaling creates no small-order additive PH subgroups: nonzero points
have order r. Neither invalid shortcut is reused.

The accepted plain v3 instead used a full-scalar hash translation and a
bounded detector with no jump table or distinguished-point store. Its cost
was dominated by up to 511 affine operations, each carrying a Fermat
inverse. This revision chooses **an arithmetic implementation change**:
keep that exact map and replace the scalar-action internals by Jacobian
arithmetic with one normalization. Larger DP budgets or a different walk
would introduce new success/lookup premises; they are unnecessary for this
4.99-bit claimed improvement. The explicit old/new ledgers, not a presumed
Cheon speedup or a walk-quality guess, select the method.

## Unchanged canonical map and scalar recovery

Set lambda=z^2-1 mod r for public z=0xd201000000010000. Choose the public
cube root beta satisfying [lambda]P0=(beta*x0,y0). Public arithmetic checks
lambda^2+lambda+1=0 and beta^2+beta+1=0. Negation and this endomorphism give
six-element nonzero orbits, N=(r-1)/6. The action is free because r is prime.

Choose minimum X among its three endomorphism images and minimum Y among
Y,p-Y. The flag plus 48-byte smaller Y is a complete key: equal keys give
equal X^3 and the same orbit. Infinity has a zero key. X=0/Y=0 would give
orders dividing 3/2, impossible in nonzero G1. Retain canonical multiplier
epsilon; canon(-Q)=canon(Q), epsilon(-Q)=-epsilon(Q) handle both signs.

Maintain an affine canonical R=[a]P0+[b]T with a,b reduced mod r. With a
fresh 32-byte seed and the unchanged domains, evaluate

    h(K)=OS2IP(SHAKE256("orbit6-v3" || seed || K,64)) mod r
    Q=R+[h(K)]P0
    Rnext=canon(Q), using its multiplier epsilon
    anext=epsilon*(a+h) mod r, bnext=epsilon*b mod r.

The unchanged complete 49-byte key, nine-byte prefix and seed total
90 bytes; input/output fit one SHAKE rate block and hashing is paid.
Initial labels use `initial-v3` and two reduced 512-bit blocks. Zero b
fails; initial infinity yields a verified linear relation.

B is the least integer with 100B^2>=169N. Floyd has at most B iterations,
three transitions per iteration. A capped two-pointer entry replay costs
at most 2B further transitions, then exactly one cycle lap costs at most
B. Total work is bounded by **6B** even on failures. There are no restarts,
DP reports, checkpoint escapes or unbounded sampling loops. First-Floyd-
meeting labels are not used directly: a meeting can span multiple laps
whose orientation products cancel. Entry versus one lap avoids that bug.
For equal canonical states the relation is

    (a_entry-a_lap)+x*(b_entry-b_lap)=0 mod r.

A nonzero denominator gives x=(a_entry-a_lap)/(b_lap-b_entry) mod r,
verified by [x]P0=T before output. A zero denominator fails this attempt.
Identity T returns zero. The wire encoder emits exactly the required
64-lowercase-hex scalar object. The reference API operates on decoded
points; native canonical decoding, curve/subgroup validation, verification
and output are charged in the fixed work allowance.

## Why success is not worsened by the arithmetic change

H1 is a **score-critical, unproved cryptographic random-function model**
for domain-separated SHAKE labels at distinct keys and initial labels,
with the inherited 21/20 production slowdown allowance. Finite seed entropy
is not claimed to create statistically independent uniform scalar pairs.
This is an explicit extrapolation for a computed cryptographic function,
not a theorem attributed to a rho paper. Its production validity remains
material uncertainty.

In the ideal uniform-label model, R+[h]P0 is uniform over G1 at each newly
exposed key. Conditional on a nonzero next class, its six orientations are
equiprobable independently of that class. The orientation product around
one simple cycle is therefore uniform among the six units. With initial
b nonzero, the denominator vanishes precisely for product 1. The useful
one-lap fraction is **5/6**, not an asserted negligible 1/r failure.

A sufficient-event model is

    (5/6)*(1-exp(-(20/21)^2*B*(B-1)/(2N))) >= 11/25.

An exact positive rational exponential series proves this numeric bound.
Subtract bias (6B+2)*r/2^512 and exceptions (B+2)/r, each below 2^-124,
then a named **3/100 H1 reserve**; the result exceeds 2/5>39/100. This is a
named budget widening under the random-function model, not an empirical
transfer-error bound. Exact point equality
means the arithmetic change adds no mixing or DP-transfer assumption.

## Jacobian formulas and paid scalar actions

The scalar action computes [h]P0 projectively, then mixed-adds the old
canonical affine R, and only then normalizes. With h<r there are at most
255 Jacobian doubles and 256 mixed calls. A projective (X,Y,Z) represents
(X/Z^2,Y/Z^3); infinity is a separate case.

The complete explicit formulas are in `argument.md` and executed by
`projective_action`. Doubling costs seven field products. Mixed addition
costs eleven, including its same-point fallback (four products then seven
for doubling). Opposite points return infinity. Infinity and Y=0 branches
precede any division. These exceptions are included in the ceiling.

Normalize once with binary Z^(p-2): <=381 squares+381 products, then four
scaling products, at most 766. Canonicalization has eight more kernels,
including Montgomery decoding. Total **255*7+256*11+766+8=5375**; squares
use the multiply kernel. No serial affine step borrows a batched price.

## Exact instruction ledger and conservative claim

At protected C=2000, H=0 and all arithmetic is charged through R. The
pinned field multiplication costs 383 scalar AArch64 instructions; it is
not also charged as a reference addition. This is an analytical native
ledger, not Python timing or unit-cost wide arithmetic.

Per projective call, 32 six-limb modular linear steps at 64 instructions
plus 1536 pointer/call/spill/case/bit/loop instructions fit the 4096 charge.
A modular linear step's limb census is 12 loads+6 carry arithmetic+6
conditional subtraction+6 selection+6 stores+16 addressing/control =52,
below 64. The 381-round inversion control and packing get another 8192.
SHAKE has 24 rounds at <=1024 plus <=16384 packing/setup, below 65536.
512-bit scalar reduction is 512 rounds at 128. Two scalar coefficient
products use <=2*255 rounds at 256, below 131072. Outer canonical-key,
state-transfer and three-limb-counter work gets 8192 separately.

| Per-transition row | Instructions |
|---|---:|
| 5375 field kernels at 383 | 2058625 |
| 511 projective calls' nonkernel work at 4096 | 2093056 |
| Inversion control/packing | 8192 |
| SHAKE256 including domain/packing | 65536 |
| 512-bit reduction | 65536 |
| Two scalar normalizations | 131072 |
| Complete keys/counters/transfers/outer control | 8192 |
| **Total I** | **4430209** |

The fixed 2^40-instruction allowance exceeds 64 complete scalar/inversion/
validation routines at 2^32 each plus 2^24 other instructions. It covers
public-constant derivation, input validation, seed/initial labels and points,
memory initialization, final modular inverse, verification and output.
Cold preprocessing is a subset, below 2^40/2000<2^30. No future-target
amortization, advice construction or unpriced table lookup is assumed.

    W=(6B*4430209+2^40)/2000
     =807602141107780472377226959275329862357980191/500.
    140.2126 < log2(W) < 140.2127.

Multiply **all non-field-kernel rows and the full fixed allowance** by
6/5; the pinned field-kernel price itself is unchanged. Then

    Wstress=4470337107344540525639193268558441290940527771/2500.
    140.3594 < log2(Wstress) < 140.3595 < 140.48-0.10.

The base slack is deliberately 0.2673..0.2674 bits; repriced headroom is
more than 0.1205 bits. A maximum base slack of 0.06 and minimum stress
headroom of 0.10 cannot both be met when Wstress>=W. Like prior accepted
v3, this package states the conflict and takes the stronger stress reserve;
it does not relabel padded work as a measured cost or make a false narrow-
slack certificate. This is the one requested presentation target not met,
for an explicit mathematical reason, rather than an unsupported score cut.

## Captured execution, histograms and limitations

Both unit tests were observed failing before the new action/solver was
implemented, then passed against independent affine arithmetic, including
identity, cancellation, doubling and identical-coin pipeline comparisons.
The expanded capture below is separate evidence, not a statistical claim
from two unit tests.

Commands from this package directory in the public repository:

    python3 -B test_reference.py
    python3 -B toy.py --capture
    python3 -B verify.py --capture
    python3 -B toy.py --verify
    python3 -B verify.py

The captured and rerun pipeline uses actual curves Y^2=X^3+4 with
(p,r,cofactor) (30013,10039,3), (30109,10093,3), (30139,9931,3), 800 trials
each. Group orders are counted and p,r prime-checked. Actual point arithmetic
checks every returned scalar and intermediate label invariant. Every
nonzero point is checked for six-class keys and signs; **30063** actions
and **72** identical-coin old/new pipelines agree with the affine oracle.

| p / r | Recoveries | Failures: Floyd / identity-product | Maximum evaluations / cap |
|---|---:|---:|---:|
| 30013 / 10039 | 514/800 | 189 / 97 | 298/324 |
| 30109 / 10093 | 485/800 | 214 / 101 | 324/324 |
| 30139 / 9931 | 504/800 | 203 / 93 | 288/318 |

Total **1503/2400**. Twelve losslessly split raw shards retain all failures.
Each is below 6 KiB and the compact summary below 10 KiB; pretty-printing
initially exceeded its file budget, not a correctness gate. A complete
second run matches the frozen capture; no trials are retried or discarded.

This is the actual 49-byte coordinate key, SHAKE domains, full-scalar map
and bounded detector, not a substitute whole-exponent hash. No branch table
or DP predicate exists in production, so no toy substitutes a denser DP
predicate. `toy-output.txt` prints six-orientation and h mod16 diagnostic
histograms; every bin occurs. Complete period/evaluation histograms are
also captured. Counts include correlated repeated visits and are not used
as independent samples or a confidence bound on production H1. Toy success
can exceed the sufficient birthday event because Floyd sometimes meets
within B even if mu+period>B. Neither costs nor production probabilities
are inferred from the observed 1503/2400.

## Certificate, storage and what remains uncertain

`verify.py` regenerated **106/106** exact checks, using Python integers and
Fractions. Log brackets and the score use integer-power inequalities, not
floating-point logarithms. Literal certificate output is embedded in
`argument.md` and shipped as `captured-output.txt`; `certificate.json`
binds protected-file digests, parameters, exact ledger, all checks and the
toy record hashes. Normal verification refuses to overwrite stale evidence.
Only public-constant identities/actions are evaluated at production size;
no full-production challenge search was attempted.

Sixteen 256-byte affine state slots hold 48-byte coordinates, 32-byte
scalars, full keys, counters and header. Separate code, scratch, hash state,
constants and I/O bring the inventory to **405504<2^20 bytes**. Jacobian
temporaries fit scratch (64 field values use 3072 bytes). No expected DP
count is turned into peak memory. This is native storage, not Python RAM.

H1 and analytical native pricing remain production uncertainties. Every
bounded attempt, replay and scalar action is charged. Acceptance would
mean plausible_not_refuted, not novelty or demonstrated production recovery.
