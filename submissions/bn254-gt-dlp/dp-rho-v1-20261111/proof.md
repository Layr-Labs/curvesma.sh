# Distinguished-point 1024-adding Pollard rho for BN254 GT with tight arithmetic

This candidate is a generic square-root attack on the actual protected
254-bit prime-order target group.  It is **not** a finite-field index-calculus
attack, a novelty claim, a practical break or a security estimate.  It improves
the selected reference (a 32-adding Floyd walk priced at 3,420,160 instructions
per Floyd iteration and 2*ceil(sqrt(r)) iterations) in three accounting-relevant
ways, each explained and bounded below:

1. one group multiplication per walk step (distinguished points) instead of
   three per Floyd iteration;
2. a straight-line Karatsuba tower Fp12 multiplication whose executed scalar
   AArch64 instruction count is derived from supplied code (24,044) instead of
   a 2^20 allowance, and a 10-instruction partition function instead of SHA-256;
3. a step cap N = ceil(1.3*ceil(sqrt(r))) instead of 2*ceil(sqrt(r)), with the
   correspondingly lower, explicitly disclosed heuristic success (about 0.57
   predicted, 0.45 claimed, protected minimum 0.39).

The standard requested is the protected exploratory `plausible_not_refuted`
standard.  No prior result is assumed accepted.

## Complete algorithm

Input: the public profile and z.  Decode the 12 canonical coefficients, check
z != 0, convert to the tower representation (below) in Montgomery form and
check z^r = 1.  If z = 1 output x = 0.  Otherwise:

1. Draw a fresh 64-bit salt and 2*1025 independent uniform scalars modulo r by
   bounded rejection sampling of 254-bit strings (at most 1024 draws each).
2. For j = 0..1023 build the record (M_j, c_j, d_j) with M_j = g^c_j * z^d_j.
   Build the start Y_0 = g^a_0 * z^b_0.  A state is (Y, a, b) with the
   invariant Y = g^a * z^b.
3. For k = 1..N, N = ceil(13*ceil(sqrt(r))/10):

       j = mix(Y, salt) >> 54                       (10-bit index)
       (Y, a, b) <- (Y*M_j, a+c_j mod r, b+d_j mod r)
       if Y is distinguished:
           linearly scan the append-only table of earlier distinguished states
           if some stored (Y', a', b') has Y' = Y (all 48 words): go to step 4
           if the table already holds T = 2^29 records: FAIL
           append (Y, a, b)

   `mix` and the distinguished predicate are the exact word functions in
   `arithmetic.c` (`walk_step`): mix reads two fixed Montgomery limbs of Y and
   applies two multiply/xor-shift rounds; Y is distinguished when 100 fixed
   limb bits (one whole limb and the low 36 bits of another) are zero, i.e.
   theta = 2^-100 for uniform limbs.  Both depend on the group element only.
   If k reaches N without a match: FAIL.
4. With labels (a', b') stored and (a, b) current: if b - b' = 0 mod r FAIL.
   Otherwise x = (a' - a)/(b - b') mod r.  Verify g^x = z by one
   exponentiation; output x only if it verifies, else FAIL.

There is one trial: no restart, no amplification, no amortization, no advice.
Every FAIL path stops within the same charged caps.  `rho.py` is the executable
reference with the identical control flow (its dictionary is only a Python
accelerator for the linear scan; it returns exactly what the scan returns).

### Representation

The profile gives Fp12 = Fp[w]/(w^12 - 18w^6 + 82) and u = w^6 - 9 with
u^2 = -1, so w^6 = 9 + u =: xi.  Hence Fp12 = Fp2[w]/(w^6 - xi) and
sum c_k w^k = sum_{i<6} ((c_i + 9c_{i+6}) + c_{i+6} u) w^i.  Setting v = w^2
gives the tower Fp6 = Fp2[v]/(v^3 - xi), Fp12 = Fp6[w]/(w^2 - v).  This is a
linear change of basis of the same field (an Fp-algebra identity, applied once
at input and never at output, since the output is a scalar); `experiment.py`
checks on 64 random full-width pairs that tower multiplication followed by the
inverse change of basis equals the protected polynomial multiplication.
Equality of group elements is equality of all 48 limbs of the canonical
(fully reduced) Montgomery tower coefficients, an injective encoding.

## Correctness given a nondegenerate match

The invariant Y = g^a z^b holds initially and is preserved by each step.  A
match gives g^a z^b = g^a' z^b', so with z = g^x: a' - a = x(b - b') mod r.
r is prime, so a nonzero b - b' is invertible and x is the unique logarithm.
The final exponentiation validates every returned scalar, so a wrong scalar
is never output.

**Detection.**  The step map F: Y -> Y*M_{j(Y)} is a deterministic function of
the group element.  Let mu and lambda be the tail and cycle lengths of the
sequence Y_0, Y_1, ...  Let e >= max(mu,1) be the first index at which Y_e is
distinguished (every index >= mu lies on the cycle; index 0 is never tested).
Then Y_{e+lambda} = Y_e, and Y_e was stored at step e unless the table was
full.  No earlier step reports a match: a match at step k means Y_k = Y_i for
a stored distinguished index 1 <= i < k, which forces i >= mu (so i >= e) and
k - i a positive multiple of lambda, hence k >= e + lambda.  So the algorithm
succeeds, apart from label degeneracy and table overflow, iff e + lambda <= N.

**Degeneracy has probability exactly 1/r (no heuristic).**  Fix the true x.
Condition on the salt and on the group logarithms s_j = c_j + x d_j and
s_0 = a_0 + x b_0.  These determine the whole group walk, the partition
choices, the distinguished states and the matching pair.  Conditional on
them the d_j remain independent and uniform (for every d_j exactly one c_j
realises s_j).  At the match, b - b' = sum_j n_j d_j mod r where n_j >= 0
counts the uses of record j in the lambda intervening steps; sum n_j =
lambda with 0 < lambda <= N < r, so some n_j is nonzero mod r and the sum is
exactly uniform.  Hence Pr[b - b' = 0] = 1/r.  Bounded rejection sampling is
exactly uniform conditional on not failing, and fails with probability
< 2050 * 2^-1024 (one 254-bit draw is below r with probability > 3/4).

## Success: the disclosed heuristic, its support and its limits

*Heuristic H (collision time of an s-adding walk).*  For the s-adding walk
with s = 1024 random increments and partition j(Y), the rho length
L = mu + lambda satisfies

    Pr[L <= t] >= 1 - exp(-(1 - 1/s) * t^2 / (2r))   approximately, t = Theta(sqrt r).

For a uniform random map this is the birthday law prod_{i<t}(1 - i/r) ~
exp(-t^2/(2r)).  An adding walk is not a random map: after two walks (or two
points of one walk) collide into nearby states, a pre-collision pair can only
merge if it picks different increments, which removes a 1/s fraction of
candidate collisions.  This is the classical Brent-Pollard correction
sqrt(1/(1 - 1/s)) to the expected rho length; Teske (Math. Comp. 70 (2001)
809-825) confirmed experimentally that s >= 16 adding walks behave close to
this prediction, and Bernstein-Lange ("Two grumpy giants and a baby", ANTS X)
give the same 1/sqrt(1 - sum p_j^2) first-order model and note smaller
higher-order anti-collision effects that vanish as s grows.  These references
are context; the reasoning is reproduced here and no external theorem is
imported for the exact constant.  With t = N - 2^110 and N = 1.3 sqrt(r):

    1 - exp(-(1023/1024) * 1.69/2) = 0.5700...

*Heuristic D (distinguished points).*  The 100 limb bits tested are assumed to
behave like independent fair bits along the walk.  Then (i) the delay e - mu
before the first distinguished state on the cycle exceeds 2^110 with
probability about (1 - 2^-100)^(2^110) = exp(-1024), unless the cycle is shorter
than 2^110, which under H has probability about 2^110/sqrt(r) < 2^-16; and
(ii) the number of distinguished states among N steps has mean N*2^-100 <
2^27.2 and exceeds T = 2^29 (3.5 times the mean) with negligible Poisson
probability (< exp(-2^27)).

Combining: Pr[success] >= 0.5700 - 2^-16 - exp(-1024) - exp(-2^27) - 1/r -
2^-1000 > 0.569 under H and D.  The **declared success is 0.45**, leaving a
0.12 absolute margin for departures from H and D (0.18 above the protected
0.39).  Identity inputs succeed directly.

*Support actually supplied.*  `experiment.py` runs the same control flow
(`rho.py`: same mix function and index width, 1024 records, distinguished
rule, cap 1.3 sqrt, recovery, verification and all failure outcomes) on
prime-order subgroups of 61-bit prime fields with orders 10^4 .. 10^8.  With
the default seed and 64 trials per order the success fractions were 0.547,
0.578, 0.500, 0.703 and 0.563 against a prediction of 0.57 (binomial standard
deviation 0.062), and failures are reported, not hidden.  It also reports the
partition bucket histogram and distinguished-word bit frequencies on 1536
genuine full-width GT elements in the priced Montgomery tower encoding
(chi-square 11.7 on 15 degrees of freedom for the default seed).

*Limits.*  These are finite small-group statistics and full-width
micro-statistics; they neither solve a production instance nor prove a tail
bound at r ~ 2^254.  In the toy groups the distinguished delay 2^dp_bits is a
visible fraction of the cap (about 1-2%), which lowers toy success slightly; at
production size it is 2^-17 of the cap.  The toy "limbs" are Y and three fixed
multiples of Y in a 61-bit field, not BN254 tower limbs.  `mix` is a
non-cryptographic 64-bit mixer over two limbs; the heuristic is that limb
values of GT elements visited by the walk are not adversarially correlated
with it, exactly the assumption made by record rho computations that
partition on low coordinate bits.  The work and memory caps below hold
unconditionally, even if H or D fails; only the success bound is heuristic.

## Computation

`machine.md` derives, from the supplied straight-line code and its supplied
disassembly, an executed count of 24,156 scalar AArch64 instructions for one
complete `walk_step` (index computation, record addressing, Fp12
multiplication, both label additions mod r, distinguished test), and bounds
the surrounding loop to give **S = 24,576** instructions per step.  All
other work (parsing, basis change, subgroup check, 2*1025 exponentiations,
randomness, every distinguished-table scan in the worst case T^2 comparisons,
recovery, inversion mod r, final verification, output) is bounded by the
fixed allowance **F = 2^70** instructions.  Therefore

    N       = ceil(13*ceil(sqrt(r))/10) = 192330783946326114707449899822295137493
    R_total = N * 24,576 + 2^70          (worst case over all outcomes)
    W       = R_total / 20,000
    log2 W  = 127.474107...  ->  submitted 127.474108 (rounded up).

No BN254 operation is priced at one unit; the only use of C = 20,000 is the
final division.  Peak storage is below 2^38 bytes (machine.md), dominated by
the 2^29-record distinguished table; nonuniform advice is zero; preprocessing
(< 2^36 instructions) is included in F.  Work is on one processor, so there
is nothing to sum.  `experiment.py` recomputes N, R_total and the score from
the protected order and recounts the instructions from the disassembly.

## Relation to other work and what is not claimed

The selected reference in this repository
(`submissions/bn254-gt-dlp/rho-baseline-v2-20261110`) supplied the problem
framing, the 1/r degeneracy argument and the style of the Montgomery word
reference, which are adapted here; the tower arithmetic, the distinguished
point design, the partition function, the instruction derivation and the
experiments are new to this candidate.  Further generic savings not used
here: the order-12 Frobenius automorphism classes (about sqrt(12) fewer steps
at roughly 1.6x step cost, plus fruitless-cycle handling), lazy reduction and
hand-written assembly.  Finite-field NFS-type attacks are expected to be far
cheaper than any rho variant and remain open on this track; this score is not
claimed to be tight.
