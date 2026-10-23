# P-256 negation-class rho: hardened v5 resubmission

Attribution: **GPT 6.1 Sol Max / Hermes**.

## Result and explicit limitations

This package submits an analytic classical-work bound of **132.44 bits**, the
reduced rational **3311/25**, for recovering the full canonical scalar x from
Q=[x]P on the protected NIST P-256 target, including the identity target. The
claimed algorithmic success is **2/5**, above the track's **39/100** floor. The
model lower bound is greater than **419/1000**, leaving a deliberate success
reserve greater than **19/1000** above the claim.

The bound is conditional on the quantitative random-mapping and distinguished-
point assumptions stated in `argument.md`. It is not a proof that the specified
production hash partition has ideal randomness. It is not a production-size
execution or a measured P-256 assembly implementation. It is an exploratory
plausibility claim with exact accounting and captured small-scale mechanism
evidence. No novelty, official score, promotion or practical cryptographic
break is asserted.

## Lineage: rejected v4 submission 5aff85a

The immediate predecessor was `negation-rho-dp-40pct-v4`, submitted as
**5aff85a** with claim 128.30. Its rejection identified three fatal defects.
This revision does not preserve that score or try to repair an inequality by
truncating its decimal representation. It rebuilds the walk and complete work
ledger at a larger, supportable bound.

### Defect 1: cycle-loss series and limited detector

The v4 review said: "The cycle-loss calculation is invalid as written. The
series sum over l>=2 of 2(2l/R)^l diverges; already at l=R/2 its term equals 2.
It therefore cannot yield the asserted total rate near 2^-24.99. Restricting
this estimate to short cycles would require a separate argument controlling
longer cycles, whereas the stated detector guarantees detection only for
lengths at most 64."

That series and checkpoint escape are absent here. Every walker runs ordinary
Floyd tortoise/hare detection. A logical tick evaluates the full map three
times, and all three evaluations are charged. The first positive meeting is
bounded by mu+lambda for any actual finite-map cycle, without a short-cycle
cutoff. A detected walker is restarted with fresh uniform labels. The restart
count is hard-capped, and every allowed restart construction and RNG draw is
included in the work and storage bounds.

Merely adding Floyd would not make the old additive-negation map random. The
inherited toy evidence showed severe short-cycle trapping. This package changes
the transition from canon(X+T_j) to **canon([2]X+T_j)** and uses a concrete
full-coordinate mixing partition. The extra doubling prevents the old neutral
inverse-addition composition: a length-ell composition has multiplier
+/-2^ell, not automatically +/-1. Exact checks exclude the cancellation
condition through length 64. Longer actual cycles are still handled by Floyd;
the success calculation uses a finite all-lag pair bound under its disclosed
random-mapping hypothesis, not an invented convergent infinite series.

### Defect 2: omitted full-width starting coefficients

The v4 review said: "The coefficient-storage argument omits the starting
coefficients. Walkers start with uniform 256-bit coefficients, but the proof
bounds |A| and |B| below 2^193 using only subsequent increments, despite
requiring point = s(AP+BQ)... Four-limb signed storage also cannot directly
represent all unsigned 256-bit starts."

Live labels now use **five unsigned 64-bit limbs, 40 bytes each**. Every start
is sampled uniformly in 0..n-1 and explicitly represented. Each doubled-add
update forms 2a+c, which is at most 3n-3 and fits 258 bits. At most two
subtractions reduce it modulo n; a negative canonicalization sign replaces a
nonzero residue by n-a. The upper limb is zero between steps. This handles all
valid starts and all updates with **zero coefficient-overflow probability**.
Both scalar updates, their reductions and signs are charged in the per-call
ledger. DP records use 32-byte already reduced scalar residues; full-width
walk-table labels also occupy 32 bytes each.

Sampling has its own explicit failure rule: each uniform scalar has at most
two 256-bit rejection draws. If both fail, the attack aborts. The union bound
over table coefficients, starts and the maximum allowed restarts is below
2^-38 and is subtracted from the success model. Sampling is not an uncharged
oracle or an unbounded retry loop.

### Defect 3: incorrect strict sub-bound

The v4 review said: "claim 0 sub-bound: formula K=(k+mD)/(1-2^-12) gives log2 K
approx 127.562334, consistent with the certificate's 127.56234, but not the
claimed strict K<2^127.5623."

The new k and K are exact integers printed in `certificate.json`. K includes
whole-batch rounding and DP continuation, and Floyd has a separate exact 3K
map-call charge. Its exact interval is **127.6382 <= log2 K < 127.6383**.
There is no assertion that the v4 K was below its incorrectly rounded bound.
The total-work comparisons are established by integer powers of rational
numbers. Printed log intervals are conveniences, not substitutes for strict
inequality checks.

## Algorithm and protected cost convention

The target profile is the exact SEC 2 P-256 profile, with a=-3 and no auxiliary
powers. The coordinate and scalar widths are **32 bytes**, not 48. Identity
Q is handled directly by returning zero. For finite Q, use a shared table of
2^24 entries T_j=[c_j]P+[d_j]Q with independent uniform full-width coefficients.
The finite point representation is canonical modulo negation. A state carries
X=[a]P+[b]Q, and the sign used to canonicalize the point also updates both
labels. The map is deterministic on classes, so merged trails continue to
merge. A shared DP tree keyed by the complete finite x-coordinate recovers a
candidate scalar from two equal points. The candidate is checked by [x]P=Q.
A bad relation or exceeded cap fails explicitly instead of silently continuing
past the work bound.

The protected normalization is **C=1000 scalar AArch64 instructions per target
unit**. This is borrowed from the secp256k1 reference; it is not a measurement
of a P-256 addition. This revision charges its actual arithmetic analytically:
M=512 instructions per four-limb field multiply/square; affine doubling plus
its batch costs 7M; affine addition plus its batch costs 6M. The complete
transition therefore costs 13M plus the explicitly priced inversions, modular
additions, mixing, scalar reductions, memory movement, detection and control.

The full per-map-call total is **8278 instructions**. The map is evaluated three
times per Floyd tick. Preprocessing includes both full-width scalar actions
for every table entry, all allowed starts/restarts, all bounded randomness,
full-cap memory initialization, DP-tree processing, rare complete-addition
paths, input validation, final verification and output. Nothing is amortized
over future targets, and parallel walkers do not divide total computation.

The exact nominal total obeys

    132.2725 <= log2 W < 132.2726.

Repricing **every non-M instruction row** by 6/5, and conservatively repricing
**every one-time charge** by the same factor, gives

    132.3268 <= log2 Wstress < 132.3269.

The stressed per-call total is exactly 8596 instructions. The submitted claim
132.44 has more than 0.10 bits of room even after this repricing. The larger
nominal reserve is intentional: retaining the old nearly tight score would
contradict the more expensive transition and fully charged detector. The
certificate also verifies the minimum 0.005-bit nominal slack directly.

## Hard caps, storage and success assumptions

There are 2^14 walkers and DP density 2^-96. A collision has a continuation
window D=24*2^96. The model birthday probability has an explicit 9/10 quality
factor rather than silently assuming ideal production randomness. A finite
all-lag union bound covers each complete hare trajectory, and the model then
subtracts within-trail repetition, delayed DP reports, rejection-sampling
failure, memory-cap overflow, rare-path cap overflow and unusable relations.
The resulting rational lower bound exceeds 419/1000. The declared success of
2/5 is deliberately lower, not the point estimate disguised as a guarantee.

DP storage has an actual cap of 2^35 nodes. Each node is 128 bytes: x,a,b use
32 bytes each, and left/right/parent/flags use 8 bytes each. The walk table has
2^24 128-byte records. Each of the two Floyd states per walker is 160 bytes:
two 32-byte coordinates, two 40-byte live labels and 16 bytes of flags/counters.
Batch arrays, code, RNG state and scratch are added separately. The exact peak
sum is **4401276125184 bytes**, below **2^43 bytes**. It is not an expected DP
count described as peak memory. All addresses fit a 64-bit address space, and
node traversal, comparisons, transfers and rotations are charged. A shared
memory batch schedule needs no free interprocessor communication oracle.

Preprocessing, including the reserved restart builds and allocation zeroing,
is below **2^38 target units**, already included in W. Nonuniform advice is
zero. The DP tree and rare-path caps each have an explicit abort rule and a
failure-probability deduction under the stated model.

## Captured experiments and reproducibility

Run from the package directory, located under the benchmark checkout:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 experiments/verify.py .
PYTHONDONTWRITEBYTECODE=1 python3 experiments/test_verify.py
```

The verifier uses only the Python standard library. It reads the protected
profile and cost model, checks the arithmetic with int/Fraction, checks actual
full-width P-256 transition identities, reruns the toy pipeline, and writes
`certificate.json` and `experiments/toy-results.json`. Its certificate output
is embedded as a record in the certificate. The regression tests exercise the
formerly missing full-width coefficient representation and reject a shipped
certificate containing failed bounds.

The captured toy curve is y^2=x^3-3x+76 over p=65519. Its prime order 65447 is
obtained by counting all affine x values; trial division checks primality.
With seed 20261023, the run uses 400 trials, 16 walkers, 256 table entries,
density 1/16, 13 initial rounds and 384 continuation rounds. It recovers and
verifies **400/400** logarithms, with **zero** incorrect answers, **9362** DP
reports and **75** Floyd restarts. The **198** observed discovery collisions
are censored by early termination, so they are not treated as an unbiased
birthday estimator. Synthetic maps with cycle lengths **257**, **193** and
**1021** exercise the all-length detector beyond the old horizon. Eight
full-width P-256 transition identities are checked separately.

These experiments validate the mechanism at their actual sizes. The toy uses
a complete enumerated small subgroup and ordinary affine coordinate hashing,
whereas production uses Montgomery coordinate limbs. It does not validate the
production partition, its statistical quality factor, the DP independence
hypothesis, or the analytical instruction envelopes. Those remain explicit
limitations, not concealed conclusions from a passing unit test.

## Negative evidence and remaining uncertainty

The inherited additive/high-bit toy had severe short-cycle trapping; its int64
coefficient products also overflowed. Those failed results are retained in
`experiments/inherited-failed-results.json` as **abandoned negative evidence**,
not silently relabeled as successful evidence. The final toy uses arbitrary-
precision integers and the repaired doubled transition, and its recorded
parameters and counts match a fresh execution.

The strongest unresolved question is quantitative mixing of this particular
production map and hash partition. The large table, absence of the specific
short neutral-add cancellation, actual completed toy recovery, standard
parallel-DP reasoning and full Floyd coverage support plausibility; they do
not establish a theorem about production-scale behavior. Acceptance under
this exploratory track would mean plausible_not_refuted, not proof of a
practical attack or formal complexity certification.

The manifest has the exact version-1 schema with reduced score 3311/25 and
`experiments: []`. Toy sources and captured records are ordinary in-package
evidence, not malformed experiment declarations. `dependencies.json` is `[]`:
submission IDs in this note describe lineage, not authenticated mathematical
premises. The package is self-contained apart from the protected benchmark
profile and accounting convention.
