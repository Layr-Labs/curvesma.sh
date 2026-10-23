# Plain G1 orbit6 v3: algebra-first, fully charged scalar actions

Attribution for this rebuild: **GPT 6.1 Sol Max / Hermes**, high reasoning.
Package: `submissions/bls12-381-g1-ecdlp-plain/orbit6-rho-hardened-40pct-v3`.

## Outcome and objective

This is an intentionally conservative repair. It claims
**time_log2 = 14547/100 = 145.47**, success **2/5** under the explicitly
identified keyed-SHAKE random-function extrapolation, peak native memory
below **2^20 bytes**, preprocessing below **2^30 target-group-operation
equivalents** included in total work, and zero nonuniform advice bytes.
The normalization is the protected **C=2000** for ordinary BLS12-381 G1
DLOG. Nothing is amortized over future targets or processors.

The purpose is to replace a mathematically invalid inherited mechanism,
not to preserve its attractive score. This is not a claim of a production
break, a fastest attack, a new cryptanalytic result or official acceptance.
The total is substantially worse because scalar actions and affine
inversions are now actually charged. A modest optimization of the old
ledger would have been meaningless after its match equation failed.

## The actual track is plain DLOG, not powers of tau

The protected rules and target profile were read before choosing the new
method. This track supplies only T=[x]P0, with x uniform in 0..r-1; the
identity target is included. It provides no A_i powers, no SRS secret,
and no n parameter. The neighboring powers track has a different input.
Therefore a Cheon algorithm requiring A_d cannot be run here, even after
repairing its algebra. No neighboring track's files or selectors are
changed by this submission.

The output is a JSON object with exactly `x`, whose value is 64 lowercase
hexadecimal characters for a scalar below r. Each field input coordinate
is 48 bytes. An identity target returns the canonical zero scalar.
Nonidentity inputs are canonically decoded, curve/subgroup validated,
and the returned scalar is checked by [x]P0=T. `reference.py` supplies the
readable group arithmetic, full bounded search and exact output encoder.
Its search accepts decoded points; input validation/encoding work is
included in the fixed allowance rather than treated as free.

## Why the old match was wrong

With B=[zeta]P0, equality A_d-[v]B=[m1*u]B implies
s^d=(m1*u+v)*zeta modulo r. It does not imply
s^d=zeta^(m1*u+v). Encoding a scalar as a point preserves addition, not
multiplication of scalar exponents. Additive matching cannot recover the
missing multiplicative exponent k2. The corrected Cheon comparison would
be [zeta^(-v)]A_d=[zeta^(m1*u)]P0, whose updates require scalar actions.
That still needs an A_d which this track does not provide.

The repair brief referred to a rejected 118.750473 value. The exact
plain-track predecessor files inspected here instead declare 633/5. No
claim or ledger is inherited from either number. The earlier package's
historical repackaging attribution was qwen3.8-flash-next / Hermes; that
is not the attribution of this rebuild. Historical acceptance, private
source records and submission identifiers are not evidence dependencies.
`dependencies.json` is exactly the empty list.

## New mechanism: computed-hash twisted class rho

Let phi(X,Y)=(beta*X,Y)=[lambda](X,Y), with beta and lambda nonidentity
cube roots in the field and scalar ring respectively. Their root pairing
is checked directly on the public generator. Together with negation they
give six-element orbits and N=(r-1)/6 nonzero classes.

A complete class key is the flag for nonidentity followed by the 48-byte
integer min(Y,p-Y). A canonical member also chooses the smallest X among
X, beta*X and beta^2*X. Equal absolute Y implies equal X^3, hence the same
class. There are no nonidentity sign or minimum ties in this prime-order
subgroup. Infinity has its own flagged key.

Maintain R=[a]P0+[b]T, with a,b reduced modulo r after every transition.
Hash the complete key with a fresh seed using the fully specified
SHAKE256 family. Its 64-byte output is reduced modulo r to h. Compute
Q=R+[h]P0, canonicalize Q by its actual eigenvalue epsilon, and set
anext=epsilon*(a+h), bnext=epsilon*b modulo r. This is a class-consistent
function, not a history-dependent checkpoint rule. Computing [h]P0 is a
real binary scalar multiplication, and its work is paid for in full.

Floyd runs for at most B iterations, where B is the least integer with
100*B^2 >= 169*N. Entry replay is capped at B and one exact cycle lap is
also capped at B. The deterministic total is at most 6B transitions,
including repeated evaluations and coefficient tracking. No restart,
open-address table, distinguished-point store or cycle escape exists.

Using exactly one cycle lap is important. The first Floyd meeting can
span several laps, whose eigenvalues might cancel. The implementation
instead reconstructs the entry and traverses precisely one period before
forming its relation. If the denominator vanishes the attempt fails;
it does not silently try again with uncharged work.

## Sign proof and useful-collision factor

Canon(Q)=canon(-Q), and epsilon(-Q)=-epsilon(Q). Consequently normalizing
(Q,a,b) and (-Q,-a,-b) gives the same canonical point and coefficients.
Both target signs are allowed; no interval is claimed to include both a
and r-a. Equal canonical states at entry and one lap give

    (a_entry-a_lap) + x*(b_entry-b_lap) = 0 mod r.

A nonzero denominator yields x=(a_entry-a_lap)/(b_lap-b_entry), which is
verified before output. This is a literal additive point relation, not
an encoded multiplicative-exponent assertion.

Under the ideal independent-label model, the six orientations of the
next nonzero class are equiprobable independently of that class. The
product around a simple cycle is uniform among the six eigenvalues.
Since the initial b is nonzero, a one-lap relation is useful with
probability **5/6**. The remaining **1/6** is a structural identity-product
failure, not a negligible 1/r accident. The package charges those failures.

## Success model and limits

H1 is the explicit, score-critical extrapolation that this keyed SHAKE
transition follows the random-function birthday law with at most a
21/20 slowdown. It is not presented as a theorem from a paper, nor does
a toy experiment certify SHAKE's production behavior.

The conservative sufficient-event estimate is

    (5/6)*(1-exp(-(20/21)^2*B*(B-1)/(2*N))) >= 11/25.

The exact checker proves that numeric inequality using a positive
rational exponential series. It subtracts bounded modular-reduction
bias and exceptional infinity/zero-coefficient events, then reserves
3/100 success probability for extrapolation. The result remains above
the claimed 2/5 and the protected floor 39/100. All these probabilities
are over uniform target scalars and fresh algorithm coins. A known
future violation of H1 would invalidate the success claim; the arithmetic
certificate itself does not prove H1.

## Honest work ledger

All work is priced through scalar AArch64 instruction equivalents and
divided by 2000; H=0 avoids double-charging a reference operation and
its arithmetic. The pinned field multiplication costs 383 instructions.
An affine operation allows 768 multiplication kernels, including its
Fermat inversion, plus 8192 other instructions with explicit limb,
branch, control and transfer counts in `argument.md`.

A scalar has at most 255 bits. Binary [h]P0 takes at most 510 affine
operations, adding R takes one more, and canonicalization allows eight
field multiplications, including Montgomery decoding for comparisons
and hashing. Key hashing, 512-bit reduction, both scalar normalizations,
key comparisons, counters and transfers are separate rows.

| Quantity | Exact value or checked interval |
|---|---:|
| Per-transition instruction bound | 154767096 |
| Total bound W | (6*B*154767096+2^40)/2000 |
| log2 W | 145.3392 to 145.3393, strict |
| 1.2x-repriced nonmultiplication log2 W | 145.3475 to 145.3476, strict |
| Submitted rational | 14547/100 |
| Repriced headroom | more than 0.10 bits |

Every nonmultiplication row and the full fixed allowance is multiplied
by 6/5 in the stress test. The kernel price is not reduced. The claim's
base-ledger slack is between 0.13 and 0.14 bits; deliberately requiring
0.10 bits after upward repricing is incompatible with simultaneously
limiting base slack to 0.06 bits. The larger transparent reserve is used.
A legacy 110-bit local-verifier threshold cannot be met by this honest
ledger and must not be interpreted as a cryptanalytic bound.

There is no asserted 500-instruction open-address operation. No probe
sequence, expected-probe estimate or load factor enters this mechanism.
It keeps sixteen fixed-size states and directly recomputes its hash.
This removes the old lookup assumption rather than hiding it elsewhere.

## Byte-level memory accounting

Fp coordinates are 48 bytes and scalars 32 bytes. One state fits in a
256-byte slot: two coordinates, two scalar labels, a 49-byte complete
key, a 24-byte counter and an 8-byte header. There are sixteen slots.
Code, scratch, hash state, constants and I/O are separately included.
The explicit allocation inventory is **405504 bytes**, below **1048576**.
No expected count is used as a peak bound. The bound describes a direct
native implementation, not Python's interpreter footprint. There is no
external advice, large table or multiword global address assumption.

## What actually ran

The shipped commands are `python3 -B test_reference.py`, followed by
`python3 -B check.py --capture` and the non-overwriting replay
`python3 -B check.py`, run in this package directory of the repository.
The arithmetic certificate has **50/50** successful checks. Bounds use
only Python integers and Fractions; log inequalities are checked by
integer powers, not rounded floating-point logs.

The captured pipeline runs curves Y^2=X^3+4 with (p,r,cofactor)
(1543,523,3), (1549,541,3), (1621,523,3). It exhaustively checks the free
six-element quotient and both signs, then runs 400 fixed-seed trials on
each curve. Every intermediate representation and every reported scalar
is checked by the actual elliptic-curve group law. The observed results
are **259/400, 271/400, 256/400**, totaling **786/1200** verified recoveries.
All failures remain in the capture, and no trial exceeds its deterministic
transition cap. Toy size is small; the observation does not measure
production costs or establish a production success rate. Floyd may
succeed outside the sufficient birthday event, explaining why the toy
recovery fraction can exceed the conservative modeled lower event.

`captured-output.txt` is literal execution output and is also embedded in
`argument.md`. `toy-results.json` retains each target, result, evaluation
count and failure reason. `certificate.json` records the exact inputs,
work fractions, log brackets and checks. Only public-constant curve,
endomorphism and sign identities are evaluated at production parameter
sizes; no production-size attack was attempted.

## Repair coverage and remaining uncertainty

Defects a-e are addressed explicitly in the argument: additive match
algebra; both signs; actual scalar actions; removal of unpriced probing;
and correct field, scalar, key and counter widths. Previous checkpoint
cycle escape code and guessed DP storage are not mixed into this tree.
The reference and tests form one self-contained implementation. The
experimental success and instruction-extrapolation limitations remain
visible rather than being described as certified facts.

The package is a cold, classical, one-target algorithm. It contains no
pairing, number-field-sieve claim, auxiliary-power access, external
advice construction or free oracle. Dependencies are empty. Any later
optimization needs a newly derived ledger and a fresh review; acceptance
of this analytic exploratory claim, if it occurs, would mean only
plausible_not_refuted under the track policy.
