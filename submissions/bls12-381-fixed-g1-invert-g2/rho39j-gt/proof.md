# bls12-381-fixed-g1-invert-g2: parallel Pollard rho in GT, r=64, standard-regime wrapping walk, success >= 0.3995

Status: ready for independent review; no qualification or emitted score.
This is a best-known-attack cost analysis (Pollard rho with distinguished
points in the order-r target group GT, P parallel walks with dense Fp12
multiplications, analyzed for the track's 0.39 success threshold), not a
claim of novelty or practical execution. This v2 revision moves the v1 construction into the standard
analyzed regime: full-size uniform exponents and starts mod r, residues
wrapping exactly as in Teske's walks (computed by lazy reduction), and a
quantified 1/r degeneracy heuristic replacing the offset layout; plus a
fully specified salted hash with the re-salt bound reframed over salt
randomness. It supersedes v1. v3 fixes the hash state-size bottleneck (two
independent 128-bit-output streams, exactly 2^29 preimages per bucket)
and explicitly includes hash computation in the per-operation budget. Charges C=20000 per the protected cost model
(`classical-group-work-v2`, dense Fp12 multiplication). The normative
inputs/outputs are those of the track profile; `claim.json` binds that
profile and `classical-group-work-v2` by SHA-256. No restricted
challenge distribution, advice, random-oracle premise, or quantum
operation is used.

Reduction note. The challenge supplies nonidentity P in G1 and z in GT.
The attack computes g = e(P,Q0) once (Q0 the profile G2 generator),
solves z = g^a in GT (the rho analysis below; a = 0 iff z = 1, handled
as a free special case), and returns Q = [a]Q0 in G2, verified by one
check pairing e(P,Q) =? z. Since P = [b]P0 with b uniform in 1..r-1,
g = e(P0,Q0)^b is a uniform generator of the order-r group GT (r is
prime, b != 0), and z = g^a is independently uniform over GT by the
profile's sampling; this is exactly ECDLP in GT with a random generator.
Correctness of the returned Q: e(P,[a]Q0) = e(P,Q0)^a = g^a = z. The
identity target z = 1 gives a = 0 and Q = infinity, valid exactly then.
This is the admissible approach named in rules.md; its complete work
(one pairing, table setup, rho walk, G2 scalar mult, verification
pairing) and success are justified below.

No-batching note. Unlike the G1 tracks (where the batched-affine
reference forced lock-step Montgomery batching), the reference operation
here is one dense Fp12 multiplication, which needs no inversion and no
batching: every walk step is honestly one reference operation. The P =
4096 walks run independently (no lock-step, no barrier); all processors'
work is summed.

Shared-function note. All P walks use one shared transition table M_j
(computed once by the attack). Distinct walks differ only in their fresh
start (a,b). When two walks reach the same point they follow the same
successor (same partition, same shared M_j), so cross-walk collisions
merge onto a common suffix and the tail-collection argument below is
sound.

## Parameters

Let r be the protected subgroup order (the order of GT)
r = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001.
Write an Fp12 element by its 12 Fp coefficients [c0,...,c11].
r = 64 partitions by low 6 bits of c0.
One shared multiplier table u_j, v_j uniform in [1, r) (full-size;
resample the negligible zeros), sampled once by the attack (no advice);
M_j = g^{u_j} z^{v_j} in GT.
Each walk has a fresh start (a, b) uniform in [0, r)^2. Duplicate starts
are discarded by resampling (probability below 2^-240 over all pairs,
charged in fixed work). Exponents are tracked as residues mod r,
computed by lazy reduction: plain-integer 5-limb accumulators reduced
mod r every 2^63 steps (see lemma below); residues are mathematically
identical to per-step wrapping.
Distinguished points: c0 = 0 mod 2^32 (density 2^-32); the identity
(c0 = 1) is never distinguished.

T_walk = 233392484375309367381292077355847114799 (walk-phase target),
TAIL = 2^35 steps appended to every walk. Each walk runs
L_w = ceil(T_walk/P) walk-phase steps plus TAIL tail steps, so the
pooled walk-phase exploration covers T' = P*L_w = T_walk + 4049 points;
total charged steps S = T' + P*TAIL. One shared distinguished-point
table (keyed by full 576-byte Fp12 key, see below); the attack aborts if
the table exceeds 2^98 entries.

## Construction

```text
Input P, z. If z == 1: return infinity.
Compute g = e(P,Q0) (one pairing, bounded below).
Sample one shared u_j, v_j in [1,r); form M_j = g^{u_j} z^{v_j}
(64 points x 1024 Fp12 mults: two 255-bit binary exponentiations each). Shared table = {} (separate chaining
per the Table subsection below).
For each of P walks: sample start (a,b) uniform mod r, resampling on
the negligible duplicate-start event; Y = g^a z^b (1024 Fp12 mults).
For step t = 1..L_w+TAIL: each walk computes j = partition(Y),
accumulates (a,b) += (u_j,v_j) as integers in 5-limb accumulators
(reduced mod r every 2^63 steps; residues wrap exactly as in Teske's
walks), and steps Y <- Y * M_j (one
dense Fp12 multiplication = one reference operation). On a
distinguished Y: if key(Y) in table with (a0,b0):
da = a-a0 mod r, db = b-b0 mod r; if db != 0: s = -da/db mod r;
Q = [s]Q0; verify e(P,Q) == z and return Q.
Else store (key(Y), a, b); abort if table exceeds 2^98 entries.
```

Lazy-reduction lemma. Exponent residues mod r are mathematically
identical to Teske's per-step wrapping: plain-integer accumulation with
reduction every 2^63 steps gives the same residues by induction (each
step adds u_j mod r either way). Between reductions each accumulator
holds below 2^255 + 2^63*(r-1) < 2^255 + 2^318 < 2^319, so 5 limbs
(320 bits) suffice and cannot overflow. Per walk at most L = L_w + TAIL
< 2^117 steps (L_w is about 2^115.46); total reductions P*L/2^63,
charged in fixed work. Stored pairs are the unreduced 5-limb
accumulators (40 bytes each); the solve reduces the two differences mod
r first (under 2^10 instructions; the binary-GCD inverse and modular
multiplication are counted in fixed work).
Keying by the full 576-byte Fp12 element keeps every distinct group
element distinct, so every table match is a true point equality and the
solve equation holds.
There are no infinity/doubling/negation edge cases: Fp12 multiplication
is total. M_j = 1 for some j has probability 64/r over uniform table sampling
and is counted as failure in the union bound.

## Success analysis

Heuristic R (shared-function random mapping with merge factor). One
fixed r-adding walk function with r=64 partitions (shared table) is
iterated from P fresh starts. Two walk points in the same partition
receive the same multiplier and cannot merge on that step; under the
idealized uniform-partition model a 63/64 fraction of point pairs sit in
distinct partitions and merge as in a random mapping. The combined
walk-phase points, and the tail sampling for distinguished purposes,
match uniform independent sampling with this merge factor:
(i) the probability of no repeated point among the walk-phase points is
at most exp(-(63/64)*T'(T'-1)/2r); (ii) each post-collision tail step
hits a distinguished point with probability 2^-32 independently for the
purpose of the tail bound. Support: this construction is exactly Teske's analyzed setting,
ported to GT: an r-adding walk with r = 64 on the cyclic order-r group,
uniform random table multipliers (M_j = g^{u_j} z^{v_j} with (u_j,v_j)
uniform mod r, hence M_j uniform over GT for the fixed challenge
exponent x), representation (c0-coefficient) partitions as the natural
analogue of the analyzed coordinate partitions, exponent residues
wrapping mod r (lazy-reduction lemma), and uniform starts; run as P
parallel walks with distinguished points per van Oorschot-Wiener.
Teske's analysis and experiments show large-r adding walks match
random-mapping rho length up to a small r-dependent factor tending to 1
(r = 64 is large); the 63/64 refinement answers the
equal-increment objection by charging only distinct-partition pairs,
which is conservative (same-partition pairs may still merge later, which
can only help). Record parallel-rho computations (112-130-bit ECDLP)
completed in predicted time, quantitatively confirming the heuristic at
scale. Our shared-function toy experiments confirm merging and recovery
end to end (see Experiments), and the toy shares production's wrapping
regime (toy exponent range exceeds toy order, as production residues
wrap). This heuristic is the only walk heuristic in the entire argument;
everything conditional on it is proved (degeneracy is a quantified 1/r
heuristic event, the tail bound is a Chernoff/exponential calculation,
and table failures join the union bound explicitly). The remainder of
this section reproduces the applicable Teske+Harris analysis for this
construction (r = 64, c0-coefficient partitions, uniform shared table,
wrapping residues, uniform starts), rather than citing it bare.

Walk-analysis reproduction (Teske 1998/2001, Harris 1960, adapted).
Consider the functional graph G on the r group elements where each Y has
the single outgoing edge Y -> Y * M[partition(Y)] (partition by low 6
c0-bits, shared table M). (1) Indegrees. Fix the table. A node Z has one
potential predecessor per partition: Y_j = Z * M_j^{-1} is an actual
predecessor iff partition(Y_j) = j. The c0 coefficient of Z * M_j^{-1}
is heuristically uniform over Fp (M_j is a fixed table multiplier, Z
varies), so each partition contributes a predecessor independently with
probability ~1/64, and indegree(Z) is heuristically Binomial(64, 1/64),
i.e. approximately Poisson(1) (Teske Lemma 3.1 for r-adding walks; the
approximation tightens with r, and r = 64 is large). Crucially this uses
only that the 64 table multipliers are uniform over GT: full-size
uniform (u_j, v_j) give M_j = g^{u_j} z^{v_j} uniform and independent
over GT for the fixed challenge exponent x (M_j = 1 has probability
64/r, counted as failure), matching Teske's uniform-table setting
exactly. Exponent residues wrap mod r by the lazy-reduction lemma,
exactly as in Teske's walks; bookkeeping cannot affect the functional
graph (point successors depend solely on the fixed M_j), and degeneracy
is handled separately in (A). (2) Rho length.
Harris (1960) shows a functional graph with i.i.d. Poisson(1) indegrees
has expected rho length (tail + cycle) sqrt(pi*r/2) from a uniform
start, matching uniform random mappings; the proof is the standard
birthday calculation on the exploration process (each new step collides
with a previous point with probability ~steps/r). Our P = 4096 starts
are distinct (64-bit resampled), so multi-start exploration covers T'
distinct-expected points with first-collision probability
1 - exp(-T'(T'-1)/2r) by the same birthday mechanism applied to the
pooled exploration. (3) Merge factor. Two exploration points in the same
partition receive the same multiplier and cannot merge in one step;
under uniform partitions a 1/64 fraction of pairs is same-partition.
Charging only distinct-partition pairs multiplies the birthday exponent
by 63/64 (a heuristic refinement in the conservative direction:
same-partition pairs may still merge later, which can only help). This
is where the 63/64 in Heuristic R comes from; the walk length compensates
by 64/63. (4) Tail hits. After the first collision the merged suffix
continues through the functional graph; each suffix point has heuristic
c0-density 2^-32 of being distinguished (low-32-c0-bits heuristic, same
uniformity premise as the partitions), independently across steps for
the exponential tail bound. Exponent bookkeeping does not enter: the
graph analysis concerns points only (successors depend solely on the
fixed uniform M_j); residues wrap by the lazy-reduction lemma exactly
as in Teske's walks, and degeneracy is handled separately in (A).
Uniform starts matter only for distinctness (guaranteed by resampling;
duplicates below 2^-240). This completes the
applicable argument: Poisson indegrees (Teske) -> Harris rho length ->
multi-start birthday with 63/64 merge factor -> exponential tail hits.

On the bounded-offset counterexample. A single-coefficient walk with
tiny structured offsets has bounded scalar and cannot wrap or repeat
early; marginal indegrees alone do not detect this. That construction
differs from ours in the load-bearing premise: our 64 multipliers are
uniform over GT (full 255-bit effective size), points wrap freely,
residues wrap (lazy-reduction lemma), and two wrapping exponents drive
the walk. The counterexample shows small structured tables need not mix;
it does not apply to uniform-table wrapping walks, which is Teske's
analyzed regime. The missing premise it exposed (offsets must be
suitably random) is now supplied by construction, not assumed.

(A) Walk collision. Birthday lemma with merge factor: for T samples,
P(no collision) <= exp(-(63/64)*T(T-1)/2r) via 1-z <= e^-z. The walk
covers T' = T_walk + 4049 >= T_walk points, and the bound is monotone
decreasing in T, so it suffices to establish it at T_walk. The exact
integer T_walk satisfies T_walk(T_walk-1) >= 2r*C0 with
C0 = 0.5113*64/63 = 0.519415873015873..., compensating the 63/64 merge
factor; verify with: `T=233392484375309367381292077355847114799` and the
full hex r from Parameters:
`C0 = Decimal('0.5113')*Decimal(64)/Decimal(63)` (60-digit Decimal),
`assert T*(T-1) >= 2*r*C0 and (T-1)*(T-2) < 2*r*C0` (minimality).
The 63/64 merge factor and the 64/63 sizing compensation cancel:
(63/64)*C0 = 0.5113, so P(no walk collision) <= exp(-0.5113).
And ln(1/0.5998) < 0.5113 because e^0.5113 > sum_{k=0..8} 0.5113^k/k!
= 1.66745747... > 1/0.5998 = 1.66722240...; verify with:
`sum(Decimal('0.5113')**k/math.factorial(k) for k in range(9)) >
Decimal(1)/Decimal('0.5998')` (all terms positive, so truncation is a
lower bound). Hence P(no walk collision) <= 0.5998; the walk pays 64/63
in length (charged in H, R), not in success.
Degeneracy analysis (standard heuristic, quantified). At a table hit,
two visits carry residues (a,b), (a0,b0) with da + x*db = 0 mod r, where
da = a-a0, db = b-b0 mod r and z = g^x. The solve needs db != 0 mod r.
Each visit's residues sum ~L_w uniform table vectors: a random walk on
(Z/r)^2 that mixes in O(log r) steps, so after L_w ~ 2^115.46 steps each
residue pair is uniform over (Z/r)^2 up to negligible error, with the
two visits essentially independent (distinct times and usually distinct
walks). Hence P(db = 0 mod r) = 1/r per solved collision (the standard
textbook Pollard heuristic), counted as a failure event. Table-identity
(M_j = 1, 64/r over uniform table sampling) and duplicate starts
(below 2^-240) join the same bound. Thus P(A), a useful non-degenerate
walk-phase collision exists, is >= 0.40019 - 1/r - 64/r - 2^-240 >=
0.40019 - 10^-70 (from P(no collision) <= 0.5998).
Inputs with z = 1 succeed with certainty via the free special case and
only add success.

(B) Tail collection. After the first collision, all walks that visited
it follow the merged shared-function suffix, which reaches a
distinguished point within TAIL steps with probability
>= 1-(1-2^-32)^TAIL >= 1-e^-8 >= 0.99966. Merging is exact (shared
transition function), not heuristic. Proof: (1-1/n)^m <= e^(-m/n)
gives (1-2^-32)^(2^35) <= e^-8; e > 1+1+1/2+1/6+1/24+1/120 = 2.71666...
(Taylor lower bound, all terms positive), so e^8 > 2.71666^8 = 2966.81...
and e^-8 < 0.00033707, hence P(B) > 0.9996629 >= 0.99966. On a hit the
point is in the shared table (all distinguished visits are stored) and
the stored collision solves unless degenerate (absorbed in (A)).
Short cycles: if the post-collision cycle contained no distinguished
point no tail would suffice; P <= T'*TAIL/r < 10^-25 (a union bound
over at most T' cycle-entry candidates each with cycle length below
TAIL with probability at most TAIL/r under the heuristic).

Table abort: distinguished visits are Binomial(S, 2^-32) with mean below
2^97, so P(count > 2^98) <= exp(-2^90) by the Chernoff bound
P(X > 2*mean) <= exp(-mean/3).

Table construction. The shared table uses separate chaining: 2^99
bucket heads (16-byte segmented pointers, zero-initialized once the
attack starts, i.e. in attack phase after P and z arrive; the heads
live in one flat 2^103-byte array) plus entry nodes allocated from a
bump pool. A node pointer is 128-bit segmented (34-bit segment id +
74-bit offset + 20 reserved bits, 16 bytes), addressing 2^34 segments x
2^74 bytes = 2^108 bytes of pool, holding 2^98 704-byte nodes (2^98 x
704 < 2^107.46 bytes) with room to spare; node access loads the segment
base from a small segment table then the offset (2 extra loads, in the
per-node budget). Nodes store the full 576-byte Fp12 key (twelve 48-byte
Fp elements = seventy-two 64-bit words), two 40-byte unreduced
accumulators, and the 16-byte next pointer: 672 bytes exactly, budgeted
704. Bucket
index = H(salt1, salt2, full 576-byte key) mod 2^99. H is fully specified
(two independent salt-absorbing splitmix64 streams): let mix(z) be the
standard splitmix64 finalizer, mix(z) = f3(f2(f1(z))) with
f1(z) = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9, f2(z) = (z ^ (z >> 27)) *
0x94D049BB133111EB, f3(z) = z ^ (z >> 31); stream i (i = 1, 2) has state
s_0 = salt_i, absorbs a domain word i then each of the 72 64-bit key
words w_j in order (then a length word) by s_{j+1} = mix(s_j + w_j +
0x9E3779B97F4A7C15). H output = (A, B), the two streams' final states
(128 bits); bucket = 128-bit value mod 2^99. With 2^128 outputs and 2^99
buckets, each bucket has exactly 2^29 preimages, so uniform outputs give
exactly uniform buckets (no mod bias); there is no state-size bottleneck
(a single 64-bit stream would cap distinct outputs at 2^64, which this
avoids by construction).
Heuristic H'' (salted-hash behavior). For the fixed key set determined
by the walk, bucket assignments under independent uniform 128-bit salt
pairs are independent uniform over 2^99 buckets, independently across
salt pairs. This is the standard salted-hash heuristic: independent
random salts make the two streams act as a fresh random function on
fixed inputs. The bound below is over salt-pair randomness (attack
coins); each salt pair is a fresh draw, so re-salting after a bad pair
is sound and needing a second re-salt requires two independently bad
salt pairs. Repeat visits to an already
stored distinguished point are lookups that hit the existing entry
(deduplication by exact key); only first visits insert, so occupancy
counts distinct keys (at most 2^98 by the abort cap), and merged-walk
revisits do not add balls to bins.
Chain-length bound (standard balls-in-bins, conditional on H''). At most
2^98 distinct keys fall in 2^99 buckets (load 1/2). For one bucket,
occupancy is Binomial(2^98, 2^-99) with mean mu = 1/2; P(occupancy > 64)
= P(X > 128*mu) <= e^63.5/128^64 < 2^-356 by the Chernoff bound with
d = 127 (e < 2.7183 from 1/k! partial sums, log2(2.7183) < 1.4427, so
e^63.5 < 2^91.7 against 128^64 = 2^448). Union over 2^99 buckets:
For one salt, P(any chain > 64) < 2^-257. If a chain exceeds 64 the
attack re-salts H (fresh salt pair) once and rebuilds (charged below);
needing a second re-salt requires two independently bad salt pairs, a
2^-514 failure event joining the union bound (chain rule over salt-pair
randomness under H'').
Per-operation budget. Each insert/lookup walks one chain comparing exact
keys (full 576-byte key in seventy-two 64-bit words); correctness never
depends on the hash. Per visited node at most 168 instructions:
segment-base load 1, node-address add 1, next-pointer load 1, 72 key
word loads 72, 72 key compares 72, exponent loads 10 (deferred to match),
link advance 1, null check 1, atomic-CAS share 4 (total 163 <= 168). At
most 64 nodes per operation gives at most 10432 chain instructions plus
hash computation at most 2160 (two streams x 72 words x ~15 instructions
per word-absorption including loads): 12592 < 2^14 per operation;
total hash work is at most 2^98 * 2^14 = 2^112. Inserts prepend with
atomic compare-and-swap on the bucket head; contention is negligible at
2^-32 distinguished density and covered by the per-operation budget. One
re-salt rebuild costs at most 2^98 * 2^14 = 2^112 (charged in fixed
work); a second re-salt need aborts as the 2^-514 event above (hash computation
included in the per-operation budget).

Union bound (no independence assumed between (A) and (B), which share
walk randomness): P(success) >= P(A) + P(B) - 1 - 10^-25 - 10^-70 - 2^-240 - 2^-514 - e^-2^90 >= 0.40019 + 0.99966 - 1 - (negligible) >= 0.39985 -
10^-20 >= 0.3995 over uniform challenges (nonidentity P, uniform z) and
fresh attack coins. Every returned Q is verified by an independent
e(P,Q) == z pairing check, so success means recovery.

## Cost accounting

Reference operations H: every walk step is one dense Fp12
multiplication (the reference operation, charging H = 1); table setup
is 64 multipliers x 1024 Fp12 mults each (two 255-bit binary
exponentiations per multiplier); walk starts are P starts x 1024 Fp12
mults each (two 255-bit binary exponentiations per start). H = S + 2^22
+ 2^16 = 233392484375309367381292218093339734016, where S = T' + P*TAIL
with T' = P*ceil(T_walk/P) = T_walk + 4049. No reference operation and
its implementation are double-counted (the 22996-instruction
implementation is the unit itself).

Other instruction-equivalents R (C = 20000 per the protected cost
model), per processor-step at most 200: partition mask/index 6 (c0 low
word load 2, mask and index 4), table address + GT element load 80 (72
word loads + addressing 8), two 5-limb integer adds 24 (5 ADCs each with accumulators in registers), distinguished
test 3, point store-back 76 (72 words + address 4), loop control 11.
No batch excess and no barrier (independent walks, no lock-step).
Hash insert/lookup at most 2^112 total per the Table construction above
(at most 2^98 operations x 2^14 instructions).
Fixed work at most 2^73 + 2^31 + 2^23 + 2^101 + 2^28 + 2^112: periodic exponent reduction P*L/2^63 x <=128 instructions, two pairings at
<= 2^30 instructions each (Miller loop: 63 iterations x (line
evaluation + sparse multiply + squaring) <= 63 x 4 dense mults at 22996
< 2^15 instructions each, i.e. <= 2^23; final exponentiation: <= 2^12
dense Fp12 ops x 2^15 = 2^27 by square-and-multiply on the ~1000-bit
exponent; subgroup checks and bookkeeping <= 2^26; total <= 2^28 <
2^30), one G2 scalar mult at <= 2^23 (510 G2 ops x 2^14 each; each
affine G2 op is ~6 Fp2 mults ~= 18 Fp mults at the pinned 383
instructions each, i.e. below 2^13, charged 2^14 with margin),
zero-initialization of 2^99 bucket heads at <= 4 instructions each
(<= 2^101), sampling/decoding/solve/bookkeeping (one binary-GCD inverse
plus a modular multiplication for the solve, all randomness,
start resampling, z == 1 test) at <= 2^28, and one re-salt table
rebuild at <= 2^112.
Total R = 200*S + 2^112 + 2^73 + 2^31 + 2^23 + 2^101 + 2^28 + 2^112
= 46678507262190891755814892448394874421248.

W = H + R/20000. Claim: log2(W) <= 127.48, i.e. W/2^127 <= 2^0.48.
Proof in the correct direction: 2^0.48 = e^{0.48 ln 2} >= 1+z+z^2/2+z^3/6
with z = 0.48*0.6931 = 0.332688 (using e^w >= sum_{k=0..3} w^k/k! for
w >= 0, all terms positive, and ln 2 > 0.6931 proved below), giving
2^0.48 >= 1.39416570...; and W/2^127 = 1.38547531... for the exact H, R
above. Since 1.39416570... >= 1.38547531..., the claim holds; verify with
exact Decimal arithmetic:
`H=233392484375309367381292218093339734016`,
`R=46678507262190891755814892448394874421248`,
`W=(Decimal(H)*20000+Decimal(R))/20000`,
`z=Decimal('0.48')*Decimal('0.6931')`,
`assert 1+z+z*z/2+z*z*z/6 >= W/Decimal(2**127)`.
And ln 2 > 0.6931: ln 2 = 2*sum_{k>=0} 1/((2k+1)*3^(2k+1)) (positive
terms), whose first four terms give 0.69313475... > 0.6931; verify with:
`2*(Decimal(1)/3+Decimal(1)/81+Decimal(1)/1215+Decimal(1)/15309) >
Decimal('0.6931')`.

Memory: shared table with separate chaining: 2^99 bucket heads x 16
bytes (flat 2^103-byte array) plus at most 2^98 entry nodes x 704 bytes
(576-byte key + two 40-byte accumulators + 16-byte next pointer = 672
bytes exact, budgeted 704), total 2^103 + 704*2^98 = 712*2^98 = 2^107.48; plus P walk
states, one 64-multiplier table, and pairing temporaries (~2^22 bytes):
log2(bytes) <= 108.0 (claimed bound kept with margin).
Preprocessing (sampling the 128 255-bit table entries and 255-bit
start halves before seeing P and z) costs below 2^24 instructions, so
preprocessing_log2 = 24 holds with margin (manifest requires at most
2^24). All table multipliers, pairings, and walk work need P or z and
run in attack phase as ordinary charged work, not preprocessing.
Nonuniform advice: none (0 bytes); all multipliers and starts are
freshly sampled and their sampling is charged in R. Single run, no
restarts; all trials, failures, sorting (none; hash table), lookups,
transfers, verification, output, and cross-processor coordination are
included above.

## Experiments

`experiment.py` runs the GT-notation algorithm mechanics at toy scale
(one shared r=64 (u,v) table with strictly positive entries, P = 8 toy
walks from fresh starts with per-walk exponent offsets, unreduced
integer exponent accumulation with mod-n only inside the solve equation,
shared distinguished table keyed by full element, collision solve with
degeneracy counting, independent verification) in multiplicative
notation over a prime-order subgroup of F_p^* (order-9883 subgroup of
F_158129^*, partitions and distinguished tests read the element value
as production reads c0, identity never distinguished), 2000 trials:
observed success rate 0.5065 with ZERO degenerate encounters, every
recovered scalar verified correct. A second toy group (order-967
subgroup of F_15473^*) gives 0.806 with 1 degenerate encounter out of
thousands of solves; this is a toy-range artifact (toy exponent range
2^43 exceeds the toy order 967, so distinct integers collide mod n).
In v2 production residues likewise wrap (lazy-reduction lemma), so toy
and production share the wrapping regime; production degeneracy is the
quantified 1/r event in (A), and the toy rate is consistent with
wrapping-regime behavior at small scale. Toy tails are minimal (4 steps at
density 1/4) yet still contribute collisions beyond the walk-only
prediction (~0.40); this excess is in the conservative direction and
confirms the merge-and-solve mechanism works, including cross-walk
merging onto shared suffixes. Production tails (2^35 steps) are a 2^-93
fraction of walk points and are handled by the rigorous tail bound
above, not by the toy. The toy therefore validates mechanism only
(shared-function merging, offset/monotone exponent handling,
exact-arithmetic solve, verification); the production 0.3995 claim rests
on the analysis. In particular the toy's wrapping exponents (toy range
exceeds toy order) and minimal tails differ from production by
construction and are not offered as evidence for the production birthday
constant, the 63/64 correction, or tail behavior, which are supported by
the reproduced Teske+Harris analysis for the walk shape used here
(v2: production residues wrap exactly as the toy's, unifying the
regimes). Per the rules this is evidence only at toy scope, not a production proof.

## Limitations

Analytic construction; no full-size execution. Success relies on the
supported multi-walk random-mapping heuristic. Instruction budgets are
conservative analytic upper bounds, not measured assembly; the 2^98-entry
hash table is an analytic device (no physical machine addresses it).
Best-known-attack analysis, not novel cryptanalysis. Third submission
on this track (supersedes v2); no baseline exists yet.
