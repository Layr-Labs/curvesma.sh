# Hardened v7: stationary square-and-add Frobenius-class rho with Floyd detection

Attribution: **GPT 6.1 Sol Max / Hermes**. This is a new analytic candidate, not
an accepted result, production execution, practical break, or novelty claim.
The first claim is **128.65 = 2573/20** dense-Fp12-operation log2 units, success
at least **0.40 conditional on the explicitly disclosed exploratory model**.
C is **20,000**, independently checked in `tracks/classical-group-work-v2.json`.
Peak bytes are below 2^73; nonuniform advice is zero. All setup is charged.

## 1. Exact fixed-P reduction and output

For the supplied nonidentity P, compute g=pairing(Q0,P,final_exponentiate=True)
using the pinned reduced pairing. Nondegeneracy and prime order imply g != 1
and ord(g)=r, for every allowed P, not just the profile generator. If z=1,
return {"Q":null}. Otherwise find x with g^x=z, form Q=[x]Q0, and check the
reduced pairing against z for the original P. Return exactly {"Q":Q}, with
Q.x and Q.y each represented by two 96-character lowercase Fp hexadecimal
coefficients as required by the protected profile. Returning x alone is not
an output. Group/order validation and all encoding are included in setup/output.
The uniform P changes the generator, not the argument or the class count.

## 2. Class action, four products, complete tie rule

Use Fp2=Fp[u]/(u^2+1), xi=1+u, Fp12=Fp2[w]/(w^6-xi). This is the profile's
flat representation via u=w^6-1. Put lambda=p mod r. Exact modular checks in
`verify_exact.py` give ord_r(lambda)=12 and lambda^6=-1. Thus nonidentity GT
splits into (r-1)/12 free orbits; include identity to get

    N=(r-1)/12+1
     =4369656264593849206620645042348830486474212708377303151883638224994881765377.

For c1=a+b*u, d1=xi^((p-1)/6), and zeta=xi^((p^2-1)/6), exact Fp2 arithmetic
checks N(d1)=zeta, zeta in Fp, zeta^2=zeta-1, zeta^3=-1. The twelve real
w-coefficient keys are a*zeta^j (even Frobenius indices 2j), and B*zeta^j
(odd indices 2j+1), where B=a*Re(d1)+b*Im(d1). Four Fp products suffice:
a*zeta, the two terms forming B, and B*zeta; use additions/subtractions for
the other powers. `verify_exact.py` checks all twelve gamma identities against
xi^((p^k-1)/6), not a copied certificate.

Let K be the minimum key and k the corresponding Frobenius index. On a tie,
compare the full twelve reduced tower encodings of the tied conjugates and
choose the lexicographically smallest full element. Identity uses k=0. A
nonidentity full-orbit tie cannot occur. These rules make rep(y)=phi^k(y)
and (class(y),k) unambiguous. Reduced Montgomery coordinates may define the
lexicographic ordering: it is a bijective encoding, not an integer ordering
claim about decoded field elements. No ordinary-step flat-basis export is
needed. A key equality is only a possible class equality: full canonical
encodings must agree before treating a local Floyd meeting as real.

The exceptional-key component of the model below treats the coordinate
linear coincidences like uniform Fp coincidences, yielding the conditional
union estimate 3*S*66/p. This is a disclosed distribution heuristic, not a
proved bound for the production subgroup. Ties have a complete algorithmic
rule and a charged work meter; the claim does not silently identify K with a
complete invariant. The arithmetic checks do not validate this distribution.

## 3. One public table, not one random map per processor

R=1024. The public seed is a fresh, genuinely uniform **2R-block coin tape**,
with 512 bits per block. Block 2j, reduced modulo r, is alpha_j; block 2j+1,
reduced modulo r, is beta_j. This is the full seed, not a short uncharged PRNG
seed. The deterministic derivation is therefore fixed-width block extraction
and modular reduction. All processors receive identical exponents. Sampling
and building M_j=g^alpha_j*z^beta_j and all twelve twists is done once; copying
identical tables is charged separately. This avoids an added pseudorandom
exponent-generator premise. A 512-bit reduction has statistical distance at
most r/2^512 from uniform; the certificate unions this over all bounded starts
and table entries. Independent hash/mark salts and fingerprint matrix coins
are sampled and charged too; no sampling oracle or advice is free.

## 4. Why Floyd alone was not a sufficient patch

The old R-adding quotient map has algebraic inverse two-cycles. Merely adding
Floyd detects those cycles but need not preserve a random-mapping lifetime of
sqrt(N); at theta=2^-60, a walk dying after about R steps would almost never
report. Therefore v7 changes the *stationary* transition, not only its detector:

    F([y]) = [rep(y)^2 * M_h(K)],    h(K)=low_10_bits(K).

The implemented twisted update and labels are

    y <- y^2 * phi^(-k)(M_h)
    a <- 2a + lambda^(-k)*alpha_h (mod r)
    b <- 2b + lambda^(-k)*beta_h  (mod r).

This preserves y=g^a*z^b exactly and costs **two dense multiplications**, counting
squaring conservatively as a multiplication. If y'=phi^j(y), then k'=k-j;
the two successor elements differ by phi^j. Thus F is a pure class function.
There is no checkpoint, minimum-member escape, alternate index, ladder, or
history-dependent modification of F. In particular the old inverse undo
requires 2*lambda^k=-1; raising to the twelfth power contradicts 2^12 !=1 mod r.
This removes that specific degeneracy; it is not a theorem that every longer
cycle is random. The joint model below is still essential.

Use P_workers=2^20 independent starting pairs (a,b). Each active segment has
Floyd states t,h initially equal. A tortoise iteration evaluates F(t), F(F(h)):
**three evaluations, six dense multiplications**. Compare full classes if their
cheap keys coincide. Floyd's first positive meeting occurs by mu+lambda
iterations: on the eventual cycle, an iteration i>=mu with i=0 mod lambda
exists in [mu,mu+lambda], taking lambda when mu=0. At a usable meeting, align
labels, verify x, and return. Otherwise restart from a freshly charged pair.
At a distinguished tortoise point, append a report and restart; after Lmax=
128/theta=2^67 tortoise iterations without a report, restart. Distinguished
means bits 10..69 of K are zero, theta=2^-60 under the joint model. Fresh
starts, including *every* DP restart, are charged, not just initial workers.

Coalesced tortoise paths have identical relative future classes and the same
first subsequent DP, hence report/restart there in relative lockstep. Local
Floyd phases from different starts need not agree; no false claim of globally
synchronized Floyd counters is made. A premature local meeting can interfere
only on a segment containing a cycle; the bound below charges/loses all such
segments conservatively. Detector history never changes F. That distinction
is what restores the shared-map cross-processor birthday calculation.

## 5. Exact cycle law, including a valid geometric tail

For a uniformly random function on N points, the exact probability that a
specified point lies on a cycle of length L is

    c_L=(N-1)_(L-1)/N^L,    1<=L<=N,

because the L-1 other ordered cycle vertices have (N-1)_(L-1) choices and the
L edges have probability N^-L. More generally

    Pr(mu=m,lambda=L)=(N-1)_(m+L-1)/N^(m+L), m>=0, L>=1, m+L<=N.

`verify_exact.py` exhausts all 4^4 functions to check the joint law. For
independent DP marking q=theta, put u_L=(1-q)^L*c_L. Then

    u_(L+1)/u_L=(1-q)*(N-L)/N <=1-q,
    sum_(L>L0) u_L <= u_(L0+1)/q.

This is a genuine finite random-mapping law with a decreasing geometric
majorant. It is NOT the rejected, eventually growing L!*R^-ceil(L/2) series.
The point-on-cycle law is not misused as a basin-entry probability. For paths
with tails, use the joint law or the conservative birthday union bound: over
at most 2^67 starts and at most 2*Lmax states explored by a hare per segment,
the probability of any local cycle is at most

    2^67*(2*Lmax)*(2*Lmax+1)/(2N) <2^-47.

All evaluations before detection or timeout are already within the three-
evaluations-per-iteration meter. There is no separate asserted cycle-waste
factor. A no-DP segment has probability (1-q)^Lmax <=exp(-128); unioning starts
gives <2^-100. These estimates are conditional on the joint model, not unconditional
claims about the concrete deterministic F.

## 6. Aligned recovery, brute-force evidence, zero silent false answers

For raw y=g^a*z^b, y'=g^a'*z^b', and y=phi^k(y'), the correct equation is

    x=(lambda^k*a'-a)/(b-lambda^k*b') mod r.

V7 normalizes a report with its own (class,k) before discarding its field state:
A=lambda^k*a, B=lambda^k*b. Two reports whose full canonical elements match
then use x=(A'-A)/(B-B'). It is incorrect to use the unnormalized ratio on
merely equal classes. Handle zero denominators by trying the next candidate;
verify **g^x=z**, then construct and verify Q. An incorrect ratio or fingerprint
match cannot silently succeed. Labels of coalesced paths stay related by
invertible factors 2*lambda^j, so a nonzero aligned difference stays nonzero.
The ideal independent-label model gives a bad-relation union <=2^24/r<2^-230.

`toy.py --alignment` exhausts x in F13, y in F13*, both b labels, and every k
for lambda=2 of order twelve; it excludes exactly zero aligned denominators.

| k | usable pairs | aligned correct | normalized correct | raw wrong | raw undefined |
|---|---:|---:|---:|---:|---:|
| 0 | 24,336 | 24,336 | 24,336 | 0 | 0 |
| each of 1..11 | 24,336 | 24,336 | 24,336 | 22,464 | 1,872 |
| total | 292,032 | 292,032 | 292,032 | 247,104 | 20,592 |

Those are exhaustive toy identities, not a claim that a real DLP was broken.
`alignment-output.json` is the actual captured run, with all twelve rows.

## 7. Fingerprints, capacity, exact Poisson formula and finite-model safety

A report is **14-byte fingerprint +32-byte A +32-byte B +2-byte padding =80 B**,
not a fictitious 80-byte full-field (K,A,B). Its fingerprint is an independent
pairwise-universal linear hash over the 4608-bit complete reduced canonical
GT encoding: fresh random 112x4608 binary matrix plus 112-bit offset. Each
output bit is the parity of a row AND the encoding, XOR the offset. This
hash is not used by F or DP marking, so for any distinct canonical elements
its collision probability is exactly 2^-112 over its fresh coins. Hashing,
full canonicalization at a report, and matrix sampling/storage are charged.

Store at most Dcap=2^66 records in an allocated 3*2^65-slot array, nominal
load <=2/3. Abort on overflow; do not call an expectation a peak bound. Append
records during the walks and use in-place heapsort by fingerprint afterwards.
Enumerate all unordered pairs in each equal-fingerprint bucket. Each candidate
normalizes labels already in the record, solves the ratio, and checks g^x=z.
No need to store the 48-byte K or a full field element. A true collision has
the same fingerprint and is never silently discarded. Stop with counted
failure if more than Bcap=2^24 candidates are required. Every candidate costs
at most 2^15 units, so total checking is <=2^39<2^40 units.

At the nominal report count S*theta, expected distinct-element false pairs
are of order 2^18.5. At the HARD cap, their exact expectation upper bound is
Dcap*(Dcap-1)/2^113 <2^19. Pair indicators are NOT independent. For the actual
universal-hash construction we use Markov, not a fake exact Poisson law:

    Pr(false pairs >=2^24) <1/32.

This intentionally spends a visible part of the success margin. For comparison,
IF one additionally Poissonizes with X~Pois(mu), its exact tail is

    Pr(X>=B)=1-exp(-mu)*sum_(j=0)^(B-1) mu^j/j!
            =exp(-mu)*sum_(j=B)^infinity mu^j/j!.

The ratio of successive terms beyond B is at most mu/(B+1); hence the first
term divided by 1-mu/(B+1) is a valid geometric upper bound. For mu<=2^19 and
B=2^24, Chernoff gives exp(-mu)*(e*mu/B)^B <2^(-3B), since e<3 and 3*mu/B<1/8.
The formula is exact *for the stated Poisson variable only*. This much smaller
Poisson tail is not needed for the finite-fingerprint success floor.

For report overflow, mu_D=S*theta <(7/8)Dcap. The binomial mgf at 9/8 gives
Pr(D>=Dcap) <=[exp(7/64)*(8/9)]^Dcap <(255/256)^Dcap <2^-100. This is a bound
in the independent-marking model, not a measured production table occupancy.

## 8. Model, budget, and deliberate success margin

**Joint exploratory heuristic H:** the stationary square-and-add class map,
with a common fresh table and the concrete coordinate partition, has the
random-mapping first-collision/joint-cycle statistics on N classes; DP key bits
behave as independent marks along distinct classes; coordinate ties have the
uniform-field union behavior; exponent-label relations are unbiased subject
to the exposed class relation. This is an unproved production extrapolation.
It is neither an authenticated dependency nor a theorem from generic-rho
literature. The four-product identity and shared-map identity are exact, but
do not by themselves establish H. Toy support and its limits are given below.

Let

    S=ceil(sqrt(6N/5))
     =72412619877426193341163708563374841327.

Sum **S tortoise iterations over all workers**, not S per worker. Only these
paths supply the birthday opportunities; both hare evaluations are charged
but not double-counted as independent birthday samples. The ideal prediction
is 1-exp(-S*(S-1)/(2N)) approximately 0.451188. The finite S-1 correction is
negligible and included in the reserved loss. Remove local-cycle events,
mark-free timeouts, sampling bias, key ties, and unfinished segments (at most
2^20*Lmax states, relative <2^-38). The certificate allows 2^-30 total for these
small terms, including a conservative factor for perturbing the birthday
exponent. Candidate-cap loss is separately <=1/32. Thus a rational lower
bound is

    9/20 -1/32 -2^-30 >2/5.

The gap is deliberate. We do not promote the 0.451 point estimate to a 0.45
success claim, nor treat a detector/escape conjecture as an exact rate. The
quantitative conclusions follow from H and the displayed exact inequalities;
H itself is an exploratory assumption requiring reviewer judgment.

## 9. Complete resource ledger and 20% nonmultiplication repricing

One dense Fp12 multiplication costs exactly one reference unit; its internals
are not charged again. An ordinary evaluation performs two such operations.
The external instruction ledger is an analytical unrolled allowance, not a
hardware measurement of this implementation:

| external row per F evaluation | instruction ceiling / derivation |
|---|---:|
| four Fp products | 4*383=1532 (protected calibration) |
| eight 6-limb sub/neg operations | 8*24=192 |
| combine two products for B | 24 |
| reduce/normalize twelve key values | 12*24=288 |
| eleven six-limb comparisons | 11*18=198 |
| two four-limb double/add/reductions | 160 |
| table scalar/field-address and metadata loads | 192 |
| index, DP extraction, twist index | 128 |
| external state/key copies | 256 |
| cached Floyd key comparisons/control | 192 |
| two-word counters, metering, loop | 128 |
| residual branch/call/load allowance | 710 |
| total analytical ceiling | **4000** |

The conservative *scored envelope* charges **8000**, twice that ceiling, per
evaluation: 2+8000/C=12/5=2.4 units. The unused 4000 instructions are explicit
engineering slack, not asserted executed instructions. The three-evaluation
Floyd factor is fully charged. A hard meter bounds actual work and aborts on
overspend. Exceptional full ties, canonical report processing and local full
meeting checks use the additive reserve, not an unbounded free branch.

| additive row (included, never amortized) | upper bound |
|---|---:|
| every initial/DP/cycle/timeout start, cap 2^67 | 2^67 *2^10 units |
| report canonicalization, universal hash, storage, sort and routing | 2^66 *2^25/C units |
| all candidate checks, cap 2^24 | 2^24 *2^15 units |
| table construction, table/seed replication, initial setup | 2^42 units |
| exceptional equality processing and remaining output/meters | 2^70 units |
| generous total additive envelope | **2^90 units** |

A fresh start needs two <=255-bit exponentiations plus multiplication and
coin reduction, at most 2^10 units under binary square/multiply counting.
Report sorting is in-place heapsort: at most 2 log2(Dcap)+O(1) comparisons and
record moves per record. Global positions use **two 64-bit words** (67-bit
indices), not overflowing native pointers. The array is sharded over 2^16
native-addressed shards, each at most 3*2^49 slots, under 2^64 bytes. Each
random record move/compare routes a bounded 80-byte packet across <=16 binary
edges; 2^25 scalar instructions per record leaves over an order of magnitude
over the resulting <=134 comparisons/moves times packet/edge operations.
There is no free giant-RAM lookup or assumed constant-cost 73-bit address.
The universal hash uses 112*72 64-bit row operations and scalar SWAR parity,
well under this allowance. Two <=255-bit label inversions/exponent checks and
a final pairing/output fit their separate 2^15-unit check envelope. Reduced
pairing setup/output can conservatively fit 2^20 units by binary final
exponentiation and bounded Miller loops, much below the additive reserve.

The exact rational scored envelope is

    W=(36/5)*S+2^90,
    log2(W)=128.6155775656382889...,
    claim=2573/20=128.65.

The claim is **0.03442243436... bits above this envelope**, inside 0.005..0.06.
The truly itemized ordinary-row bound is smaller. Reprice EVERY nonmultiplication
row at 1.2x, including the entire additive envelope as if all of it were
nonmultiplication, obtaining

    W_stress=(168/25)*S+(6/5)*2^90,
    log2(W_stress)=128.5160418920876192...,
    claim-log2(W_stress)=0.1339581079... >=0.10 bits.

This comparison reprices the actual 4000-instruction row ceilings, NOT a
second repricing of the already doubled 8000-instruction slack. Both quantities
and that distinction are explicit; no uncertainty larger than the available
headroom is hidden. Integer cross-products in `verify_exact.py` check the
claimed power bounds without floating-point logarithms.

## 10. Peak memory, preprocessing, advice

Slot array: 80*(3*2^65) bytes, including unused slots. Identical table replicas:
2^20*12*1024*(576+32+32) bytes. Reserve an additional 2^50 bytes for all walk
states (two GT states, four scalars and cached keys/counters per worker),
code, random seed tape, the universal matrix, route queues, directories and
in-place-sort stacks. The computed hard total is

    8,854,438,289,526,828,826,624 bytes <2^73 bytes.

At overflow the algorithm fails rather than allocating beyond this cap.
No hidden second full report array is used. Preprocessing is below 2^42
operation-equivalents, included in W; the seed/table is challenge-dependent
and not amortized over future targets. External nonuniform advice: zero.
This memory is not a physically attainable resource claim; it is supplementary
and unscored under the protected convention, with its access computation charged.

## 11. Captured evidence, provenance and limits

Own standard-library Python, not executed inherited candidate code:

    python3 -I toy.py --alignment
    python3 -I toy.py --trials 1000 --seed 71
    python3 -I test_pipeline.py
    python3 -I verify_exact.py

The pipeline toy uses r=65521, lambda=53059 (order 12), 5461 classes including
identity, R=64, two workers, DP fraction 1/2, eight-bit fingerprints, and 81
aggregate tortoise iterations per trial. Actual captured output: **629/1000
verified successes**, 81,000 tortoise iterations, 243,000 F evaluations,
486,000 counted dense-operation analogues, 40,242 reports, 36 Floyd restarts,
zero timeouts, 3,013 false-fingerprint pairs tested, and zero class-equivariance
failures. It uses visible exponents and a whole-canonical-exponent SHA-256
partition; it does NOT exercise BLS field arithmetic, the four-product
production partition, production DP sparsity, instruction costs, or certify H
at production size. Its role is to expose normalization, shared-map, detector,
restart, counter and false-pair bugs. All three F evaluations are counted;
there are no hidden retry multiplications.

The exact checker returns **45/45 arithmetic checks** and tests Floyd on all
1024 (four-vertex function,start) pairs. This certifies only arithmetic and
conditional implications. It is not an independent evaluator or a proof of H.

Public prior: mpjunior92 submission
`e444842b-8632-4faa-a36e-dd455a89cfb7`, note retrieved via `yukon submission-note`.
It supplied the useful Floyd mu+lambda framing and a model-vs-claim margin
(about 0.45 vs 0.39), and used roughly 640 MiB under a 1 GiB native-memory cap.
That is context, not an authenticated accepted dependency. Its plain GT map
has no automorphism cycles; transplanting Floyd into the old quotient map
without changing the stationary map was not a justified fix. Its public note
also contains stale intermediate numerical annotations; v7 recomputes all
numbers instead of copying them. `dependencies.json` is [].

The rejected local v5 and b-v6 inform the defect map: escape re-entry -> no
escape; divergent factorial tail -> exact joint random-mapping law; raw ratio
-> (class,k) alignment; private maps -> public full-coin tape/table; false
80-byte full key -> explicit 14+32+32+2-byte record and charged false pairs;
tight/missing-cost claim -> fully counted Floyd/squaring and explicit slack.
The honest cost is higher than 125.89/126.16 because it removes unsubstantiated
savings. It remains below the cited 129.20 public Floyd bound, but no
production qualification is claimed before hosted review.
