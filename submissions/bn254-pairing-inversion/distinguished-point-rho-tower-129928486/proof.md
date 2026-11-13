# Distinguished-point Pollard rho with tower arithmetic for BN254 GT

This is a generic 32-adding Pollard-rho attack on the actual protected
254-bit prime-order GT, with two changes to the track's
rho-baseline-v2-20261110 package: Floyd tortoise-and-hare detection is
replaced by distinguished-point detection, and the loose schoolbook
polynomial-basis Fp12 multiply is replaced by a tight tower-Karatsuba
multiply. It is not a best-known finite-field attack, a novelty claim, a
practical break, or a security-bit estimate. Its work bound is analytical
and its success justification uses the standard, explicitly disclosed
Pollard-rho random-walk heuristic. This submission seeks the protected
exploratory `plausible_not_refuted` standard, not formal proof certification.
No prior result is assumed accepted.

Attribution, stated honestly: the 32-adding walk, the salted SHA-256
partition, the table construction, the step cap, the success-heuristic
structure, and the degeneracy argument below are the baseline's, reused
with the same claim level; the distinguished-point detection lemma, the
tower representation and its isomorphism argument, and the tighter cost
derivations in machine.md are new to this package. The walk itself is not
claimed as original.

## Complete algorithm

Read only the public profile and z. Decode its 12 canonical coefficients,
check z!=0 and z^r=1, and convert once to tower form (see below) with
four-limb Montgomery coefficients. If z=1 return scalar zero (or P0 and
infinity in the pairing track). Let K=2*ceil(sqrt(r)), with the protected
r. There is one bounded trial.

Sample fresh 256-bit salt and 66 independent uniform scalars modulo r,
using bounded 254-bit rejection sampling (at most 1024 draws per scalar).
For each of 32 pairs (c_j,d_j), construct M_j=g^c_j*z^d_j. From the
remaining pair form Y0=g^a0*z^b0. A state is (Y,a,b), maintaining
Y=g^a*z^b. The deterministic iteration is

    j = SHA256(salt || fixed_encoding(Y))[0] & 31
    (Y,a,b) <- (Y*M_j, a+c_j mod r, b+d_j mod r).

fixed_encoding is the injective 384-byte internal encoding specified in
machine.md (the twelve tower coefficients as canonical Montgomery
residues, in fixed order), not an unpriced conversion or hash oracle.

A state is distinguished when the low 110 bits of the big-endian integer
formed by SHA-256 digest bytes 1..14 are zero (theta = 2^-110). The test
uses the same digest already computed for the partition index, so it
costs no additional hash. Each visited state is tested before stepping,
so exactly one digest is computed per step. Distinguished states are
stored in an open-addressed hash table of 2^21 slots holding 512-byte
records (384-byte Y, 64-byte coefficients, 64-byte slot metadata
including a fingerprint); at most 2^20 records are stored, so the load
never exceeds 0.5. The slot index is the low 21 bits of the big-endian
integer formed by digest bytes 15..31, again from the same digest.
Probing is linear with a cap of 256 probes per lookup; exceeding the
probe cap or the record cap reports failure within the step cap.

At a distinguished state already present in the table, let (a1,b1) be
the stored labels and (a,b) the current labels. If b-b1=0 mod r, report
failure. Otherwise compute

    x = (a1-a)/(b-b1) mod r.

Verify g^x=z before output; verification failure is a failed trial, never
an accepted wrong answer. If no repeated distinguished point occurs within
K steps, report failure. There is no hidden restart, uncharged
amplification or amortization. `rho.py` gives the executable mathematical
reference and precise loop order.

The GT track returns x. The pairing track returns P=P0 and Q=[x]Q0 using
the explicitly charged affine procedure in machine.md. This is allowed
by the arbitrary-preimage-pair rules. Bilinearity gives
e(P0,[x]Q0)=g^x=z; both points are in the prescribed subgroups.
For the identity, the permitted infinity output is used.

## The tower representation

Fp12 arithmetic is performed in the tower

    Fp subset Fp2 = Fp[u]/(u^2+1) subset Fp6 = Fp2[v]/(v^3-xi), xi = 9+u
    subset Fp12 = Fp6[t]/(t^2-v),

which is the same field as the protected polynomial basis
Fp[w]/(w^12-18w^6+82) via the explicit embedding u = w^6-9, v = w^2,
t = w. The embedding identities are verified algebraically from the
defining relation w^12 = 18w^6-82: u^2 = (w^6-9)^2 = w^12-18w^6+81 =
-1; v^3 = w^6 = 9+u = xi; t^2 = w^2 = v. The conversion from the twelve
protected coefficients to tower form is the fixed linear map written in
`rho.py` (to_tower/to_poly), charged once in the fixed budget; the
generator g's tower form is a public profile constant. The output is a
scalar (or, on the pairing track, source-group points built from the
scalar), so no output conversion exists. The tower form is internal
representation only: the group, the problem, and the walk are unchanged.
Precedent: the baseline itself computes in Montgomery representation
internally. `experiment.py` verifies the embedding identities, the
conversion roundtrip, and tower-multiply == polynomial-basis-multiply on
200 random full-width pairs, and re-checks the priced Karatsuba operation
counts (54 Fp multiplications, 280 Fp add/subs per Fp12 multiply) by
instrumentation rather than assertion.

## Correctness given a nondegenerate repeated distinguished point

The invariant follows initially and after every step by exponent addition.
If g^a1*z^b1=g^a*z^b and z=g^x, then a1-a=x(b-b1) mod r. Since r is
prime, nonzero b-b1 has a unique inverse. The resulting x is the unique
discrete logarithm modulo r. All scalar operations use four words and are
reduced modulo r. Neither a hidden sampled scalar nor any source-group
preimage is an input. The final exponentiation validates every returned
scalar.

For a deterministic state map on group elements, write mu for the tail
length and lambda for the cycle length. Detection lemma: if the cycle
contains a distinguished point, the walk detects a repeated distinguished
point at step index at most mu+lambda+gap, where gap is the distance
from the cycle entry to the first distinguished point on the cycle. The
first distinguished cycle state is stored on its first visit (its first
occurrence is on the cycle, so no earlier visit stored it), and on its
next visit, lambda steps later, it is found in the table. Scalar labels
are not part of the hash key; detection is on group elements only. A
first-repeat length mu+lambda and a gap with mu+lambda+gap<K are
therefore sufficient for detection within the charged K steps.

## Success heuristic, support and limits

For an ideal uniform random map on r points with uniform starting point,
the no-repeat probability after t distinct transitions is the product
of (1-j/r), and is approximately exp(-t^2/(2r)) at t=Theta(sqrt(r)).
A 32-adding map is not an independent random map. The conventional
r-adding correction accounts for the probability of selecting the same
increment; its effective rate is approximately 1-1/32. Applying it at
t=2*ceil(sqrt(r)) predicts repeat probability

    1-exp(-(31/32)*2) ≈ 0.856.

This is the baseline's heuristic, reused with the same stated limits.
Distinguished-point detection adds three quantified terms, all negligible
at production size:

* Overshoot. Detection needs mu+lambda+gap<K, where gap is geometric-like
  with mean 1/theta = 2^110. P(gap > 2^120) = exp(-2^10) is zero for all
  practical purposes, so with overwhelming probability the walk needs
  mu+lambda <= K-2^120, giving P ≈ 1-exp(-(31/32)*1.991^2/2) ≈ 0.853.
* Empty cycle. The cycle holds lambda*theta distinguished points in
  expectation; lambda is typically ~2^127, so ~2^17. The failure
  probability E[exp(-lambda*theta)] over the Rayleigh lambda-distribution
  is ≈ 10^-10, dominated by the short-cycle tail.
* Record overflow. Stored records are ~Binomial(K, 2^-110) with mean
  2^17.8; P(exceeding the 2^20 record cap) is Chernoff-negligible (the
  cap is 4.6x the mean). Probe-cap overflow is likewise negligible: at
  load <= 0.5, the crude per-lookup bound gives <= 2^21 lookups x 2^-256
  (P(more than 256 probes) <= load^256 = 2^-256 per lookup). The load^256
  estimate assumes uniform hashing: digest-derived slot indices behave
  as uniform independent draws, which is part of the stated heuristic
  model. The table never exceeds its stated size; overflow is a counted
  failure, not an unbounded allocation.

The coefficient degeneracy probability is exactly 1/r by the baseline's
argument, transferred unchanged. Fix the input's true x (without giving
it to the algorithm) and condition on the salt, the 32 increment
logarithms s_j=c_j+x*d_j mod r, and the starting logarithm
s0=a0+x*b0 mod r. Conditional on these, all d_j remain independent
uniform scalars: for each d_j there is exactly one c_j realizing the
conditioned s_j. The group walk, its partition choices, its distinguished
tests, and the selected first repeated distinguished point depend only
on s_j,s0,salt, not on the remaining random d_j. At a repeated
distinguished point, the coefficient difference b-b1 is the sum of the
d_j over the lambda intervening steps, namely sum(n_j*d_j) mod r. Any
two visits to one element of a deterministic functional graph are a
positive multiple of lambda apart, and the first detected pair — the
first distinguished cycle state's second visit — is exactly lambda
apart; the argument needs only a positive intervening count below r.
The nonnegative counts satisfy sum(n_j)=lambda with 0<lambda<=K<r. At least
one n_j is nonzero modulo r. This nonzero linear combination of
independent uniform field elements is exactly uniform, so degeneracy
probability conditional on any such group walk is exactly 1/r. Bounded
rejection sampling is uniform conditional on successful sampling. This
argument applies to the actual deterministic SHA partition as well as to
an idealized partition. Thus only the probability of a group collision
within the cap is heuristic.

The joint heuristic prediction is ≈ 0.85, leaving substantial margin
over the declared **success >=0.50** (and the protected minimum 0.39).
Rejection-sampling failure is negligible: one 254-bit draw accepts with
probability r/2^254>1/2, so union bound gives failure <66*2^-1024 for
the scalar generation. Identity succeeds directly.

The material assumption is that the salted SHA partition and randomly
generated 32 increments do not cause anomalous collision-time behavior or
slow mixing at production size, and that the distinguished-point
statistics follow the same random-map model (gap geometric with mean
2^110, cycle occupancy as above). SHA collision resistance alone does
not prove them. They concern walk statistics across fresh increments,
salt, starting labels and uniform full-production input, not a claim
that a fixed seed is representative or that the whole walk is IID.
The score's work and memory caps hold even if these heuristics fail; the
success lower bound is exploratory rather than a proven unconditional
theorem.

Support is the standard r-adding construction, its random-map analysis
above, the sizeable probability margin, and directly supplied
experiments. Teske, "On random walks for Pollard's rho method",
Mathematics of Computation 70(234), 809-825 (2001),
DOI 10.1090/S0025-5718-00-01213-8, studies this class of walks; we use
32 additions rather than Pollard's less-random original three-part walk.
Distinguished-point detection is the standard van Oorschot-Wiener
method. These references are context, not external unproved theorems
imported to supply our precise probability. The relevant algorithm and
heuristic reasoning are reproduced here; no web retrieval is required.

`experiment.py` tests the actual distinguished-point walk (same 32-way
hash, 384-byte encoding length, coefficient update, cap, recovery and
failure rules, with theta = 2^-4 so distinguished points occur at toy
scale) on four prime additive groups, reporting all outcomes including
failures, and checks Montgomery multiplication on full-width edge/random
inputs for both p and r, the tower/polynomial equality, the protected
BN254 generator order, and full-width BN254 walk invariants in tower
arithmetic. At the declared seed the toy success fractions rise with
group size (0.49, 0.63, 0.84, 0.82 across the four groups), consistent
with the disclosed finite-size account: at toy scale the mean gap
(16 steps) is a large fraction of the cap and small cycles often hold
no distinguished point, while at production size the gap is 2^-18 of
the cap and cycles hold ~2^17 distinguished points. These are
correctness and finite walk-statistics evidence only. They neither
solve production DLP nor prove a probability bound at r. Additive-group
encodings differ from the field representation; extrapolation to GT and
2^127-scale walks is part of the explicitly disclosed heuristic.

This is an ordinary generic-attack improvement over the baseline.
Finite-field index calculus and NFS can be stronger and remain eligible;
nothing here proves this score is tight or is the appropriate reward
threshold for such algorithms.

## Computation, memory and output

`machine.md` gives a complete bounded scalar implementation: explicit
four-word Montgomery multiplication/REDC, tower-Karatsuba Fp12
arithmetic, scalar SHA-256, fixed table access, the distinguished-point
hash table with probe and record caps, three-word step counters, bounded
randomness, setup and failure handling, scalar recovery, and affine Fp2
output conversion with exceptional cases. Peak storage is <=2^31 bytes:
the distinguished-point table (2^21 slots x 512 bytes = 2^30 bytes) plus
at most 2^20 bytes for everything else. Memory creates no extra
objective or tie-break under the protected rules; 2 GiB is a stated,
physically plausible bound, not an unbounded or absurd-memory design.
No unbounded table, wide memory address or network transfer is assumed.

Every step costs <=87,552 scalar instruction equivalents. An additional
2^40 instructions covers all setup and final work, including the pairing
track's conversion and optional on-curve/subgroup checks, the one-time
basis conversion, and all distinguished-point table lookups. Thus for
both tracks

    R_total = 2*ceil(sqrt(r))*87,552 + 2^40,
    W = R_total/20,000.

Submit the upward-rounded six-decimal log2(W) as an exact rational in
the manifest. `experiment.py` recomputes it from the protected order; a
shared bound for both tracks is intentional, not a claim of equal actual
work. Memory log2 bytes <=31; preprocessing <=2^32 scalar instructions
is charged in full; nonuniform advice is zero. Public profile/SHA
constants and algorithm design are the only excluded public data. No
Python timings, historical measurements or normalization-unit
substitution price the arithmetic.
