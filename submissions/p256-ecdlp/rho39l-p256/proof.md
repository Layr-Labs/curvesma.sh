# p256-ecdlp: batched parallel Pollard rho, r=64, reproduced walk analysis, success >= 0.3995

Status: ready for independent review; no qualification or emitted score.
This is a best-known-attack cost analysis (Pollard rho with distinguished
points, P parallel walks with Montgomery-batched affine additions,
analyzed for the track's 0.39 success threshold), not a claim of novelty
or practical execution. This v12 revision restructures the P-256 construction per the cost
model (C = 1000 is a borrowed normalization; all P-256 arithmetic is
charged as justified instructions): the walk runs in the Montgomery
domain with an analytic CIOS field-multiplication bound (F_p256 <= 208
instructions), partitions and keys use Montgomery-form bits, H = 0, and
the claim is 2^128.44. It also carries the v12 hash coincidence
exclusion and single-bridge walk argument. It supersedes v11. Charges C=1000 per the protected p256 cost model. The normative inputs/outputs are those of the
track profile; `claim.json` binds that profile and
`classical-group-work-v2` by SHA-256 (the "production-v2" label follows
the baseline package's naming). No restricted challenge distribution,
advice, random-oracle premise, or quantum operation is used.

Batching and pricing note. On this track C = 1000 is a borrowed
normalization, not a reference operation (COST_MODEL.md): every P-256
field operation is charged as justified scalar instructions divided by
C, and H = 0. The walk runs in the Montgomery domain (all coordinates
stored as xR mod p); partitions, distinguished tests, and table keys
use Montgomery-form bits (multiplication by R permutes F_p, so this is
an exact bijection for keys and the same heuristic shape for
partitions). This construction runs P = 262144 walks in lock-step and
invert-batches the P Montgomery-domain additions of each round (one
shared Fermat inversion per round over Montgomery-form operands); each
addition costs 3M formula + 3M batch products in P-256 field
multiplications at F_p256 <= 208 instructions each (bounded below), all
charged in R.

Shared-function note. All P walks use one shared transition table M_j
(sampled once by the attack). Distinct walks differ only in their fresh
start (a,b). When two walks reach the same point they follow the same
successor (same partition, same shared M_j), so cross-walk collisions
merge onto a common suffix and the tail-collection argument below is
sound. Per-walk tables would not merge; this construction does not use
them.

## Parameters

Let q be the protected subgroup order
q = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551.
r = 64 partitions by low 6 Montgomery-form x-coordinate bits (partition 0 for infinity).
One shared multiplier table u_j, v_j uniform in [1, q) (full-size;
resample the negligible zeros), sampled once by the attack (no advice).
Each walk has a fresh start (a, b) uniform in [0, q)^2. Duplicate starts
are discarded by resampling (probability below 2^-220 over all pairs,
charged in fixed work). Coefficients are tracked as residues mod q,
computed by lazy reduction: plain-integer 5-limb accumulators reduced
mod q every 2^63 steps (see lemma below); residues are mathematically
identical to per-step wrapping.
Distinguished points: Montgomery-form x = 0 mod 2^32 (density 2^-32); infinity never is.

T_walk = 346826322834726722136455861790175095936 (combined walk-phase target),
TAIL = 2^35 steps appended to every walk. Each walk runs
L_w = ceil(T_walk/P) walk-phase steps plus TAIL tail steps, so the
pooled walk-phase exploration covers T' = P*L_w = T_walk + 37760 points;
total charged steps S = T' + P*TAIL. One shared
distinguished-point table (keyed by x AND y-parity, see below); the attack
aborts if the table exceeds 2^98 entries.

## Construction

```text
Input Q. Sample one shared u_j, v_j in [1,q); form
M_j = [u_j]P + [v_j]Q (batched Montgomery-domain mults below). Shared table = {} (separate
chaining per the Table subsection below).
For each of P walks: sample start (a,b) uniform mod q, resampling on
the negligible duplicate-start event.
For round t = 1..L_w+TAIL: each walk computes j = partition(Y),
accumulates (a,b) += (u_j,v_j) as integers in 5-limb accumulators
(reduced mod q every 2^63 steps; residues wrap exactly as in Teske's
walks), and offers Y + M_j; the P additions are performed as one
invert-batch (one shared inversion per round over Montgomery-form operands). On a distinguished Y:
if (xbar, ybar-parity) in table with (a0,b0): da = a-a0 mod q, db = b-b0 mod q;
if db != 0: x = -da/db mod q; verify [x]P == Q and return x.
Else store (xbar, ybar-parity, a, b); abort if table exceeds 2^98 entries.
```

Lazy-reduction lemma. Coefficient residues mod q are mathematically
identical to Teske's per-step wrapping: plain-integer accumulation with
reduction every 2^63 steps gives the same residues by induction (each
step adds u_j mod q either way). Between reductions each accumulator
holds below 2^256 + 2^63*(q-1) < 2^256 + 2^319 < 2^320, so 5 limbs
(320 bits) suffice and cannot overflow. Per walk at most L = L_w + TAIL
< 2^117 steps (L_w is about 2^110.03); total reductions P*L/2^63,
charged in fixed work. Stored pairs are the unreduced 5-limb
accumulators (40 bytes each); the solve reduces the two differences mod
q first (under 2^10 instructions; the binary-GCD inverse and modular
multiplication are counted in fixed work).
Keying by (Montgomery-form x, Montgomery-form y-parity) keeps P and -P
distinct (same x, opposite parity; the Montgomery map is a bijection, and
p odd gives opposite parities for y and -y), so every table match is a
true point equality and the solve equation holds. Coefficient residues
live mod q (the group order, lazy-reduction lemma); field arithmetic
lives mod p (the P-256 field, Montgomery domain): two different moduli,
each handled in its own subsection.
Affine edge cases (infinity/doubling/negation inputs) occur
with total probability at most S/q < 2^-127 and are handled by the complete
addition routine; their cost is inside the per-step budget and their failure
contribution is absorbed in the success bound. Table-identity events
(M_j = infinity, probability 64/q over uniform table sampling) join the
union bound in (A).

## Success analysis

Heuristic R (shared-function random mapping with merge factor). One
fixed r-adding walk function with r=64 partitions (shared table) is
iterated from P fresh starts. Two walk points in the same partition
receive the same increment and cannot merge on that step; under the
idealized uniform-partition model a 63/64 fraction of point pairs sit in
distinct partitions and merge as in a random mapping. The combined
T' walk-phase points, and the tail sampling for distinguished
purposes, match uniform independent sampling with this merge factor:
(i) the probability of no repeated point among the T' walk-phase
points is at most exp(-(63/64)*T'(T'-1)/2q); (ii) each
post-collision tail step hits a distinguished point with probability
2^-32 independently for the purpose of the tail bound. Support: a single heuristic bridge plus rigorous downstream. Heuristic T
(walk-as-random-mapping): the r = 64 adding walk with uniform table,
Montgomery-form-x partitions (a permuted representation; same heuristic
shape), wrapping residues, and uniform starts matches a uniform random
mapping for exploration statistics. Support for T: Teske's analysis and
experiments for large-r adding walks (r-dependent correction tending to
1; r = 64 is large); record parallel-rho computations (112-130-bit
ECDLP) completing in predicted time; our shared-function toy experiments
confirming merging and recovery end to end (see Experiments), in the
same wrapping regime as production. Everything downstream of T is
rigorous: multi-start pooling follows van Oorschot-Wiener's parallel-rho
analysis (shared iteration function from P starts, the exact setting
here), giving the birthday inequality below; the 63/64 refinement
conservatively charges only distinct-partition pairs; tail hits are
geometric at density 2^-32 on continued exploration; degeneracy is the
quantified 1/q event. On finite budgets: the asserted inequality
P(no collision) <= exp(-(63/64)T(T-1)/2q) is the birthday bound on pooled
exploration (VW's parallel model), not an extrapolation from Harris's
expected rho length; Harris and Flajolet-Odlyzko's random-mapping
analyses corroborate the quantitative form. P-256 instruction costs are
bounded analytically in the cost section, not via the toy. This heuristic
is the only walk heuristic in the entire argument; everything conditional
on it is proved (degeneracy is a quantified 1/q heuristic event, the
tail bound is a Chernoff/exponential calculation, and table failures
join the union bound explicitly). The remainder of this section
reproduces the applicable Teske+Harris analysis for this construction
(r = 64, Montgomery-form-x partitions, uniform shared table, wrapping
residues, uniform starts), rather than citing it bare.

Walk-analysis reproduction (Teske 1998/2001, Harris 1960, adapted).
Consider the functional graph G on the q group elements where each Y has
the single outgoing edge Y -> Y + M[partition(Y)] (partition by low 6
Montgomery-form-x bits, shared table M). (1) Indegrees. Fix the table. A node Z has one
potential predecessor per partition: Y_j = Z - M_j is an actual
predecessor iff partition(Y_j) = j. The Montgomery-form x-coordinate of Z - M_j is
heuristically uniform over the field (multiplication by R permutes F_p) (M_j is a fixed table point, Z
varies), so each partition contributes a predecessor independently with
probability ~1/64, and indegree(Z) is heuristically Binomial(64, 1/64),
i.e. approximately Poisson(1) (Teske Lemma 3.1 for r-adding walks; the
approximation tightens with r, and r = 64 is large). Crucially this uses
only that the 64 table points are uniform over G: full-size uniform
(u_j, v_j) give M_j = [u_j + x v_j]P uniform and independent over G for
the fixed challenge x (M_j = infinity has probability 64/q, counted as
failure), matching Teske's uniform-table setting exactly. Coefficient
residues wrap mod q by the lazy-reduction lemma, exactly as in Teske's
walks; bookkeeping cannot affect the functional graph (point successors
depend solely on the fixed M_j), and degeneracy is handled separately
in (A). (2) Rho length.
Harris (1960) shows a functional graph with i.i.d. Poisson(1) indegrees
has expected rho length (tail + cycle) sqrt(pi*q/2) from a uniform
start, matching uniform random mappings; the proof is the standard
birthday calculation on the exploration process (each new step collides
with a previous point with probability ~steps/q). Our P = 262144 starts
are distinct (uniform resampled; duplicates below 2^-220), so multi-start exploration covers T
distinct-expected points with first-collision probability
1 - exp(-T(T-1)/2q) by the same birthday mechanism applied to the pooled
exploration. (3) Merge factor. Two exploration points in the same
partition receive the same increment and cannot merge in one step; under
uniform partitions a 1/64 fraction of pairs is same-partition. Charging
only distinct-partition pairs multiplies the birthday exponent by 63/64
(a heuristic refinement in the conservative direction: same-partition
pairs may still merge later, which can only help). This is where the
63/64 in Heuristic R comes from; the walk length compensates by 64/63.
(4) Tail hits. After the first collision the merged suffix continues
through the functional graph; each suffix point has heuristic Montgomery-form-x density
2^-32 of being distinguished (low-32-bits heuristic, same uniformity
premise as the partitions), independently across steps for the
exponential tail bound. Coefficient bookkeeping does not enter: the
graph analysis concerns points only (successors depend solely on the
fixed uniform M_j); residues wrap by the lazy-reduction lemma exactly
as in Teske's walks, and degeneracy is handled separately in (A).
Uniform starts matter only for distinctness (guaranteed by resampling;
duplicates below 2^-220). This completes the
applicable argument: Poisson indegrees (Teske) -> Harris rho length ->
multi-start birthday with 63/64 merge factor -> exponential tail hits.

On the bounded-offset counterexample. A single-coefficient walk with
tiny structured offsets (P, 2P, ..., 64P) has scalar bounded by 64T and
cannot wrap or repeat before q/64 steps; marginal indegrees alone do
not detect this. That construction differs from ours in the load-bearing
premise: our 64 offsets are uniform over G (full 256-bit effective
size), points wrap freely, residues wrap (lazy-reduction lemma), and two
wrapping coefficients drive the walk. The counterexample shows small
structured tables need not mix; it does not apply to uniform-table
wrapping walks, which is Teske's analyzed regime. The missing premise
it exposed (offsets must be suitably random) is now supplied by
construction, not assumed.

(A) Walk collision. Birthday lemma with merge factor: for T samples,
P(no collision) <= exp(-(63/64)*T(T-1)/2q) via 1-z <= e^-z. The walk covers T' = T_walk + 37760 >= T_walk points and the bound decreases in T, so establishing it at T_walk suffices. The exact
integer T_walk satisfies T_walk(T_walk-1) >= 2q*C0 with
C0 = 0.5113*64/63 = 0.519415873015873..., compensating the 63/64 merge
factor; verify with: `T=346826322834726722136455861790175095936` and the
full hex q from Parameters:
`C0 = Decimal('0.5113')*Decimal(64)/Decimal(63)` (60-digit Decimal),
`assert T*(T-1) >= 2*q*C0 and (T-1)*(T-2) < 2*q*C0` (minimality).
The 63/64 merge factor and the 64/63 sizing compensation cancel:
(63/64)*C0 = 0.5113, so P(no walk collision) <= exp(-0.5113).
And ln(1/0.5998) < 0.5113 because e^0.5113 > sum_{k=0..8} 0.5113^k/k!
= 1.66745747... > 1/0.5998 = 1.66722240...; verify with:
`sum(Decimal('0.5113')**k/math.factorial(k) for k in range(9)) >
Decimal(1)/Decimal('0.5998')` (all terms positive, so truncation is a
lower bound). Hence P(no walk collision) <= 0.5998; the walk pays 64/63
in length (charged in H, R), not in success.
Degeneracy analysis (standard heuristic, quantified). At a table hit,
two visits carry residues (a,b), (a0,b0) with da + x*db = 0 mod q, where
da = a-a0, db = b-b0 mod q. The solve needs db != 0 mod q. Each visit's
residues sum ~L_w uniform table vectors: a random walk on (Z/q)^2 that
mixes in O(log q) steps, so after L_w ~ 2^110.03 steps each residue pair
is uniform over (Z/q)^2 up to negligible error, with the two visits
essentially independent (distinct times and usually distinct walks).
Hence P(db = 0 mod q) = 1/q per solved collision (the standard textbook
Pollard heuristic), counted as a failure event. Table-identity (M_j =
infinity, 64/q over uniform table sampling) and duplicate starts
(below 2^-220) join the same bound. Thus P(A), a useful non-degenerate
walk-phase collision exists, is >= 0.40019 - 1/q - 64/q - 2^-220 >=
0.40019 - 10^-70 (from P(no collision) <= 0.5998).

(B) Tail collection. After the first collision, all walks that visited
it follow the merged shared-function suffix, which
reaches a distinguished point within TAIL steps with probability
>= 1-(1-2^-32)^TAIL >= 1-e^-8 >= 0.99966. Merging is exact (shared
transition function), not heuristic. Proof: (1-1/n)^m <= e^(-m/n)
gives (1-2^-32)^(2^35) <= e^-8; e > 1+1+1/2+1/6+1/24+1/120 = 2.71666...
(Taylor lower bound, all terms positive), so e^8 > 2.71666^8 = 2966.81...
and e^-8 < 0.00033707, hence P(B) > 0.9996629 >= 0.99966. On a hit the
point is in the shared table (all distinguished visits are stored) and
the stored collision solves unless degenerate (absorbed in (A)).
Short cycles: if the post-collision cycle contained no distinguished
point no tail would suffice; P <= T'*TAIL/q < 10^-25 (a union bound
over at most T' cycle-entry candidates each with cycle length below
TAIL with probability at most TAIL/q under the heuristic).

Table abort: distinguished visits are Binomial(S, 2^-32) with mean below
2^97, so P(count > 2^98) <= exp(-2^90) by the Chernoff bound
P(X > 2*mean) <= exp(-mean/3).

Table construction. The shared table uses separate chaining: 2^99
bucket heads (16-byte segmented pointers, zero-initialized once the
attack starts, i.e. in attack phase after Q arrives) plus entry nodes
allocated from a bump pool. A node pointer is 128-bit segmented
(34-bit segment id + 72-bit offset + 22 reserved bits, 16 bytes),
addressing 2^34 segments x 2^72 bytes = 2^106 bytes of pool, holding
2^98 144-byte nodes (2^98 x 144 = 2^105.17 bytes) with room to spare;
node access loads the segment base from a small segment table then the
offset (2 extra loads, in the per-node budget). Bucket
index = H(salt1, salt2, Montgomery-form x-coordinate) mod 2^99. H is fully specified
(two independent salt-absorbing splitmix64 streams): let mix(z) be the
standard splitmix64 finalizer, mix(z) = f3(f2(f1(z))) with
f1(z) = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9, f2(z) = (z ^ (z >> 27)) *
0x94D049BB133111EB, f3(z) = z ^ (z >> 31); stream i (i = 1, 2) has state
s_0 = salt_i, absorbs a domain word i then each 64-bit key word w_j in
order (x-words little-endian, then a length/parity word) by
s_{j+1} = mix(s_j + w_j + 0x9E3779B97F4A7C15). H output = (A, B), the
two streams' final states (128 bits); bucket = 128-bit value mod 2^99.
With 2^128 outputs and 2^99 buckets, each bucket has exactly 2^29
preimages, so uniform outputs give exactly uniform buckets (no mod
bias); there is no state-size bottleneck (a single 64-bit stream would
cap distinct outputs at 2^64, which v10 avoids by construction).
Heuristic H'' (salted-hash behavior). For the fixed key set determined
by the walk, bucket assignments under independent uniform 128-bit salt
pairs are independent uniform over 2^99 buckets, independently across
salt pairs. This is the standard salted-hash heuristic: independent
random salts make the two streams act as a fresh random function on
fixed inputs. The bound below is over salt-pair randomness (attack
coins); each salt pair is a fresh draw, so re-salting after a bad pair
is sound and needing a second re-salt requires two independently bad
salt pairs.
Salt-pair resampling protocol (coincidence exclusion). After drawing a
salt pair, the attack compares the two streams' post-domain states
s^(1) = mix(salt_1 + 1 + C) and s^(2) = mix(salt_2 + 2 + C) (under 100
instructions, inside fixed work): if equal the pair is discarded and
redrawn, up to 3 draws per live pair. This is exact, not heuristic:
mix is a permutation, so equal post-domain states stay equal under
identical subsequent absorption (capping outputs at 2^64), while
unequal states, transformed by a permutation at every step, never
collide. Three consecutive coincidences (probability 2^-192) aborts as
a union-bound event. Non-coincident pairs exclude only the 2^64-output
diagonal (a 2^-64 fraction of outputs), absorbed in H''. Repeat visits to an already
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
attack re-salts H (fresh non-coincident salt pair) once and rebuilds
(charged below); needing a second re-salt requires two independently bad
non-coincident salt pairs, a 2^-514 failure event joining the union bound
(chain rule over salt-pair randomness under H'', plus the 2^-192
coincidence-abort event).
Per-operation budget. Each insert/lookup walks one chain comparing exact
keys (full 256-bit x in four 64-bit words, y-parity bit); correctness
never depends on the hash. Per visited node at most 25 instructions:
segment-base load 1, node-address add 1, next-pointer load 1, four key
word loads 4, four key compares 4, parity load+compare 2, coefficient
loads 10 (deferred to match), link advance 1, null check 1. At most 64
nodes per operation gives at most 1600 chain instructions plus hash
computation at most 128 (two streams x 4 words x ~15 instructions
per word-absorption including loads): 1728 < 2^11 per operation; total
hash work is at most 2^98 * 2^11 = 2^109. Inserts prepend with atomic
compare-and-swap on the bucket head; contention is negligible at 2^-32
distinguished density and covered by the per-operation budget. One
re-salt rebuild costs at most 2^98 * 2^11 = 2^109 (charged in fixed
work); a second re-salt need aborts as the 2^-514 event above (hash computation
included in the per-operation budget).

Union bound (no independence assumed between (A) and (B), which share
walk randomness): P(success) >= P(A) + P(B) - 1 - 10^-25 - 10^-70 - 2^-192 - 2^-220 - 2^-514 - e^-2^90
>= 0.40019 + 0.99966 - 1 - (negligible) >= 0.39985 - 10^-20 >= 0.3995
over uniform challenges and fresh attack coins. Every returned scalar is
verified by an independent [x]P == Q check, so success means recovery.

## Cost accounting

P-256 field multiplication (CIOS Montgomery, 4 limbs, analytic bound).
All walk values live in the Montgomery domain (aR mod p); partitions,
tests, and keys use Montgomery-form bits, so no domain conversion ever
occurs on walk paths. One CIOS pass, operands/accumulators in registers
(23 of 31 GPRs: a(4), p(4), T(8), b[i], m, C, temps): for i = 0..3:
m = T[0] + a[0]*b[i] mod 2^64 with p' = 1 (p equiv -1 mod 2^64, so no
multiply for m): 1 MUL + 1 ADD = 2; for j = 0..3: (C,S) = T[j] +
a[j]*b[i] + m*p[j] + C: 2 MUL + 2 UMULH + 5 ADD/ADC = 9; per i at most
2 + 36 + 2 (loop/indexing) = 40; four passes 160. Final conditional
subtract of p (4 CMP + 4 SUB + 4 CSEL = 12) and operand loads (12).
Total at most 160 + 12 + 12 = 184 <= 208 budgeted (F_p256 <= 208, with
margin). Squarings cost the same (S=M convention, conservative).
Fermat inversion I <= 512 P-256-M via left-to-right square-and-multiply
(256 squarings + at most 256 multiplies).

Reference operations H: none. C = 1000 is a borrowed normalization on
this track, not a reference operation; all P-256 arithmetic, including
point additions, is charged as justified instructions in R. H = 0.

Other instruction-equivalents R (C = 1000 borrowed normalization), per
processor-step at most 1321: six P-256 field multiplications x 208
(3M addition formula + 3M invert-batch products per member) = 1248;
partition mask/index 4, table address + point load 10, two 5-limb
integer adds 24 (5 ADCs each with accumulators in registers),
distinguished test 3, point store-back 4, loop control 7 (= 52 walk
operations); round barrier/sync share 8 (lock-step SIMT execution:
synchronization free by construction, 8 retained as margin); batch
control and addressing 12; inverse share 1 (I <= 512 P-256-M divided
across P = 262144 walks: 512/262144 M = 0.0020M x 208 < 1
instruction).
Hash insert/lookup at most 2^109 total per the Table construction above
(at most 2^98 operations x 2^11 instructions, hash computation
included).
Fixed work at most 2^73 + 2^38 + 2^45 + 2^102 + 2^109 (periodic
coefficient reduction P*L/2^63 x <=128 instructions, precompute loop
control over P*64*1024 + P*1024 adds at <=8 instructions each (<=
2^38), precompute and verification P-256 additions ((2^34 + 2^28)
batched adds x 1268 + 2^16 batch inversions x 512 x 208, all <= 2^45),
one binary-GCD inverse plus a modular multiplication for the solve, all
sampling, decoding/verification bookkeeping, start resampling,
coincidence draws, zero-initialization of 2^99 bucket heads at <= 4
instructions each, i.e. <= 2^101 inside the 2^102 term, and one re-salt
table rebuild at <= 2^109).
Total R = 1321*S + 2^109 + 2^73 + 2^38 + 2^45 + 2^102 + 2^109
= " + R256 + ".

W = R/1000 (H = 0). Claim: log2(W) <= 128.44, i.e. W/2^128 <= 2^0.44.
Proof in the correct direction: 2^0.44 = e^{0.44 ln 2} >= 1+z+z^2/2+z^3/6
with z = 0.44*0.6931 = 0.304964 (using e^w >= sum_{k=0..3} w^k/k! for
w >= 0, all terms positive, and ln 2 > 0.6931 proved below), giving
2^0.44 >= 1.35619261...; and W/2^128 = 1.34640410... for the exact R
above. Since 1.35619261... >= 1.34640410..., the claim holds; verify with
exact Decimal arithmetic:
`R=458157573767818816986334649674846300012544`,
`W=Decimal(R)/1000`,
`z=Decimal('0.44')*Decimal('0.6931')`,
`assert 1+z+z*z/2+z*z*z/6 >= W/Decimal(2**128)`.
And ln 2 > 0.6931: ln 2 = 2*sum_{k>=0} 1/((2k+1)*3^(2k+1)) (positive
terms), whose first four terms give 0.69313475... > 0.6931; verify with:
`2*(Decimal(1)/3+Decimal(1)/81+Decimal(1)/1215+Decimal(1)/15309) >
Decimal('0.6931')`.

Memory: shared table with separate chaining: 2^99 bucket heads x 16
bytes plus at most 2^98 entry nodes x 144 budgeted bytes (actual 32-byte
x + 1 parity bit + two 40-byte unreduced accumulators + 16-byte segmented next
pointer = 129 bytes), total 2^103 + 144*2^98 = 2^105.46; plus P walk
states and one batch buffer (~2^26 bytes): log2(bytes) <= 105.59
(claimed bound kept with margin).
Preprocessing (sampling the 128 256-bit table entries and 256-bit
starts before seeing Q) costs below 2^24 instructions, so
preprocessing_log2 = 24 holds with margin (manifest requires at most
2^24). All table multipliers (Montgomery-domain double-scalar mults),
inversions, and walk work need Q and run in attack phase as ordinary
charged work, not preprocessing. Nonuniform advice: none (0 bytes); all
multipliers and starts are freshly sampled and their sampling is charged
in R. Single run, no restarts; all trials, failures, sorting (none; hash
table), lookups, transfers, verification, output, and cross-processor
coordination are included above.

## Experiments

`experiment.py` runs the algorithm mechanics at toy scale (unchanged in v9; the toy already wraps) (one shared
r=64 (u,v) table with strictly positive entries, P = 8 toy walks from
fresh starts with per-walk coefficient offsets, unreduced integer
coefficient accumulation with mod-n only inside the solve equation,
shared distinguished table, collision solve with degeneracy counting,
independent verification) on y^2=x^3+x+21/F_9973 (n = 9883), 2000 trials:
observed success rate 0.509 with ZERO degenerate encounters, every
recovered scalar verified correct. A second toy curve (n = 967) gives
0.8035 with 2 degenerate encounters out of thousands of solves; these are
toy-range artifacts (toy coefficient range 2^43 exceeds the toy order 967,
so distinct integers collide mod n). In v9 production residues likewise
wrap (lazy-reduction lemma), so toy and production share the wrapping
regime; production degeneracy is the quantified 1/q event in (A), and the
toy rate is consistent with wrapping-regime behavior at small scale. Toy tails are minimal (4 steps at density 1/4) yet still
contribute collisions beyond the walk-only prediction (~0.40); this excess
is in the conservative direction and confirms the merge-and-solve
mechanism works, including cross-walk merging onto shared suffixes.
Production tails (2^35 steps) are a 2^-93 fraction of walk points and are
handled by the rigorous tail bound above, not by the toy. The toy therefore
validates mechanism only (shared-function merging, offset/monotone
coefficient handling, exact-arithmetic solve, verification); the
production 0.3995 claim rests on the analysis. In particular the toy's
wrapping coefficients (toy range exceeds toy order) and minimal tails
differ from production by construction and are not offered as evidence
for the production birthday constant, the 63/64 correction, or tail
behavior, which are supported by the reproduced Teske+Harris analysis for the
walk shape used here (production residues wrap exactly as the toy's,
unifying the regimes). The toy validates walk mechanics only, not P-256
instruction costs, which are bounded analytically in the cost section. Per the rules this is evidence only at toy scope,
not a production proof.

## Limitations

Analytic construction; no full-size execution. Success relies on the
supported multi-walk random-mapping heuristic. Instruction budgets are
conservative analytic upper bounds, not measured assembly; the 2^98-entry
hash table is an analytic device (no physical machine addresses it; the
same holds for the baseline's larger table). Best-known-attack analysis,
not novel cryptanalysis. Degeneracy is a quantified 1/q heuristic event; coefficient residues wrap per the lazy-reduction lemma; all P-256 arithmetic is instruction-charged (H = 0). This supersedes the sequential-walk candidate in
submissions/secp256k1-ecdlp/rho39/, which charged sequential additions at
the batched reference price.
