# bls12-381-g1-ecdlp (auxiliary powers): batched parallel Pollard rho, r=64, posited exploration law, success >= 0.3995

Status: ready for independent review; no qualification or emitted score.
This is a best-known-attack cost analysis (Pollard rho with distinguished
points, P parallel walks with Montgomery-batched affine additions,
analyzed for the track's 0.39 success threshold), not a claim of novelty
or practical execution. This v14 revision (a) bounds all sampling with explicit cutoffs (preprocessing table scalars: 24 attempts each; attack-phase starts: two rounds of sample, bitonic sort, dedup scan; abort on exhaustion with accounted failure mass), (b) replaces the per-step abort law with a dependence-aware bound (entries subsumed by lookups; tails trivially bounded; walk-phase lookups via no-cycling conditioning plus per-distinct-point Bernstein), (c) reframes the walk heuristic as a posited exploration law with labeled support (no derivation claimed), and (d) raises C0 0.5113 -> 0.65 for success-margin robustness (modeled 0.4775 vs claimed 0.3995). R, H, and the claim are recomputed (claim 2^127.69). It supersedes v13. Charges C=2000 per the protected BLS12-381 cost model. The normative inputs/outputs are those of the
track profile; `claim.json` binds that profile and
`classical-group-work-v2` by SHA-256 (the "production-v2" label follows
the baseline package's naming). No restricted challenge distribution,
advice, random-oracle premise, or quantum operation is used.

Batching note. The reference operation is batched incomplete affine
addition, so a single sequential walk cannot honestly charge 1 unit per
step (each dependent affine add needs its own inversion). This construction
therefore runs P = 262144 walks in lock-step and Montgomery-batches the P
additions of each round (one shared inversion per round), making every
charged addition a genuine batched addition. All batch overhead beyond
the reference 6M convention is charged explicitly in R.

Shared-function note. All P walks use one shared transition table M_j
(sampled once by the attack). Distinct walks differ only in their fresh
start (a,b). When two walks reach the same point they follow the same
successor (same partition, same shared M_j), so cross-walk collisions
merge onto a common suffix and the tail-collection argument below is
sound. Per-walk tables would not merge; this construction does not use
them.

## Parameters

Let q be the protected subgroup order
q = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001.
r = 64 partitions by low 6 x-coordinate bits (partition 0 for infinity).
One shared multiplier table u_j, v_j uniform in [1, q) (full-size), sampled once by the attack (no advice) in preprocessing: each scalar allows 24 attempts (255-bit draw, reject >= q or zero; abort the attack on exhaustion). Each walk has a fresh start (a, b) uniform in [0, q)^2, sampled in attack phase by two rounds of (sample all P starts with 24 attempts per scalar, bitonic-sort, dedup-scan), aborting on any remaining duplicate or range failure. Total sampling failure is at most 2^-59 (round-two duplicates below 2^-219 plus 255-bit range exhaustion below 2^-60; verify `assert 262144*(262144-1)*2**219 < 2*q` and `assert 1048832*(2**255-q)**24*2**60 < 2**(255*24)` with the full hex q from Parameters; sum below 2^-59), joining the union bound; sampling work is charged in fixed work (below). Coefficients are tracked as residues mod q,
computed by lazy reduction: plain-integer 5-limb accumulators reduced
mod q every 2^63 steps (see lemma below); residues are mathematically
identical to per-step wrapping.
Distinguished points: x = 0 mod 2^32 (density 2^-32); infinity never is.

T_walk = 263151378203355381386577988888241678999 (combined walk-phase target),
TAIL = 2^35 steps appended to every walk. Each walk runs
L_w = ceil(T_walk/P) walk-phase steps plus TAIL tail steps, so the
pooled walk-phase exploration covers T' = P*L_w = T_walk + 218473 points;
total charged steps S = T' + P*TAIL. One shared
distinguished-point table (keyed by x AND y-parity, see below); the attack
aborts if the table exceeds 2^98 entries.

## Construction

```text
Input Q. Sample one shared u_j, v_j in [1,q) (preprocessing: 24 attempts per scalar, abort on exhaustion); form
M_j = [u_j]P + [v_j]Q (batched below). Shared table = {} (separate
chaining per the Table subsection below).
For each of P walks: sample start (a,b) uniform mod q by the two-round setup (sample, bitonic sort, dedup scan; abort on remaining duplicates).
For round t = 1..L_w+TAIL: each walk computes j = partition(Y),
accumulates (a,b) += (u_j,v_j) as integers in 5-limb accumulators
(reduced mod q every 2^63 steps; residues wrap exactly as in Teske's
walks), and offers Y + M_j; the P additions are performed as one
Montgomery batch (one inversion per round). On a distinguished Y:
if (x, y-parity) in table with (a0,b0): da = a-a0 mod q, db = b-b0 mod q;
if db != 0: x = -da/db mod q; verify [x]P == Q and return x.
Else store (x, y-parity, a, b); abort if table exceeds 2^98 entries.
```

Lazy-reduction lemma. Coefficient residues mod q are mathematically
identical to Teske's per-step wrapping: plain-integer accumulation with
reduction every 2^63 steps gives the same residues by induction (each
step adds u_j mod q either way). Between reductions each accumulator
holds below 2^256 + 2^63*(q-1) < 2^256 + 2^319 < 2^320, so 5 limbs
(320 bits) suffice and cannot overflow. Per walk at most L = L_w + TAIL
< 2^117 steps (L_w is about 2^109.63); total reductions P*L/2^63,
charged in fixed work. Stored pairs are the unreduced 5-limb
accumulators (40 bytes each); the solve reduces the two differences mod
q first (under 2^10 instructions; the binary-GCD inverse and modular
multiplication are counted in fixed work).
Keying by (x, y-parity) keeps P and -P distinct (same x, opposite parity),
so every table match is a true point equality and the solve equation holds.
Affine edge cases (infinity/doubling/negation inputs) occur
with total probability at most S/q < 2^-127 and are handled by the complete
addition routine; their cost is inside the per-step budget and their failure
contribution is absorbed in the success bound. Table-identity events
(M_j = infinity, probability 64/q over uniform table sampling) join the
union bound in (A).

## Success analysis

Heuristic R (posited exploration law). One fixed r-adding walk function with r=64 partitions (shared table) is
iterated from P fresh starts. The analysis posits the following exploration law, with support (not derivation) below: (i) pooled birthday: P(no repeated point among the T' walk-phase points) <= exp(-(63/64)*T'(T'-1)/2q), charging only distinct-partition pairs (a 63/64 fraction); same-partition pairs are excluded conservatively (their later merges can only add collisions); (ii) single-walk self-collisions: one walk self-collides within L steps with probability at most L^2/2q (standard birthday upper bound); (iii) per-point labels: distinguished indicators are independent Bernoulli(2^-32) across distinct exploration points, independent of the walk exploration, and key bits [32,131) are i.i.d. uniform across distinct distinguished keys (bucketing); (iv) degeneracy: a solved collision has db = 0 with probability at most 2/q (mixing support in (A)). The v13 per-step-independence framing was false even in the ideal model (repeated visits share status); (iii) holds exactly there.
Support for R (evidence for the posit). Teske's analysis and experiments for large-r adding walks (r-dependent correction tending to 1; r = 64 is large; uniform-table wrapping regime with x-coordinate partitions, matching this construction); Harris's rho-length analysis corroborating the quantitative form; van Oorschot-Wiener's parallel-rho pooling model (shared iteration function from P starts, the exact multi-start setting here); record parallel-rho computations (112-bit ECDLP, Bos et al., completing in predicted time, validating the birthday constant near scale); our shared-function toy experiments confirming merging and recovery end to end (see Experiments), in the same wrapping regime as production. The posit needs only coarse accuracy: C0 = 0.65 sizes the walk for modeled success 0.4775 against the claimed 0.3995 (19% margin), so constant-factor approximation error is absorbed. This heuristic is the only walk heuristic in the entire argument; everything conditional on it is proved (the tail bound is a Chernoff/exponential calculation given (ii)-(iii), table occupancy is balls-in-bins given (iii), and table/cycling/sampling failures join the union bound explicitly).

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

(A) Walk collision. By posited exploration law (i) with the T' pooled walk-phase points: P(no collision) <= exp(-(63/64)*T'(T'-1)/2q). The exact integer T_walk satisfies T_walk(T_walk-1) >= 2q*C0 with C0 = 0.65*64/63 = 0.660317460317..., compensating the 63/64 merge factor; verify with: `T=263151378203355381386577988888241678999` and the full hex q from Parameters:
`C0 = Decimal('0.65')*Decimal(64)/Decimal(63)` (60-digit Decimal),
`assert T*(T-1) >= 2*q*C0 and (T-1)*(T-2) < 2*q*C0` (minimality).
Since T' = T_walk + 218473 >= T_walk and the bound decreases in T, (63/64)*T'(T'-1)/2q >= (63/64)*C0 = 0.65, so P(no walk collision) <= exp(-0.65).
And ln(1/0.5221) < 0.65 because e^0.65 > sum_{k=0..8} 0.65^k/k! = 1.91554076... > 1/0.5221 = 1.91534188...; verify with:
`sum(Decimal('0.65')**k/math.factorial(k) for k in range(9)) > Decimal(1)/Decimal('0.5221')` (all terms positive, so truncation is a lower bound). Hence P(no walk collision) <= 0.5221; the walk pays 64/63 in length (charged in H, R), not in success.
Degeneracy analysis (mixing support, quantified). At a table hit, two visits carry residues (a,b), (a0,b0) with da + x*db = 0 mod q, where da = a-a0, db = b-b0 mod q. The solve needs db != 0 mod q. The difference db equals the signed sum of table vectors over the symmetric difference of the two root-to-visit paths: the two distinct pre-merge segments (totaling ~L_w table vectors) plus a shared-suffix correction that cancels. Pre-merge accumulation is a random walk on (Z/q)^2 with 64 uniform generators, which mixes in O(log q) steps (standard rapid mixing of random Cayley walks; the generators span with overwhelming probability), so after L_w ~ 2^109.63 steps db is 2^-100-close to uniform over Z/q in total variation; we charge P(db = 0 mod q) <= 2/q per solved collision (covering the mixing error with room). Table-identity (M_j = infinity, 64/q over uniform table sampling) joins the same bound. Thus P(A), a useful non-degenerate walk-phase collision exists, is >= 0.4779 - 2/q - 64/q >= 0.4779 - 10^-70 (from P(no collision) <= 0.5221; verify `66*10**70 < q` with the full hex q from Parameters).

(B) Tail collection. After the first collision, all walks that visited it follow the merged shared-function suffix, which reaches a distinguished point within TAIL steps with probability >= 1-e^-8-2^-180 >= 0.99966. Merging is exact (shared transition function), not heuristic. Given exploration law (iii), distinguished indicators are independent across distinct suffix points at density 2^-32; conditional on the suffix being simple (no self-collision; failure <= TAIL^2/2q < 2^-180 by posit (ii), verify `2**70*2**180 < 2*q` with the full hex q from Parameters), its 2^35 distinct points each hit independently, so P(no hit) = (1-2^-32)^(2^35) <= e^-8. Proof: (1-1/n)^m <= e^(-m/n) gives (1-2^-32)^(2^35) <= e^-8; e > 1+1+1/2+1/6+1/24+1/120 = 2.71666... (Taylor lower bound, all terms positive), so e^8 > 2.71666^8 = 2966.81... and e^-8 < 0.00033707, hence P(B) > 0.9996629 - 2^-180 >= 0.99966. On a hit the point is in the shared table (all distinguished visits are stored) and the stored collision solves unless degenerate (absorbed in (A)). The old per-step-independence framing was false even in the ideal model (repeated visits share status); the per-distinct-point framing above holds exactly there and is posited of the construction.

Table abort (entries or total lookups exceed 2^98). Entries are subsumed by lookups (each entry needs a visit: N_entries <= L), so one bound covers both. Split total lookups L = L_walk + L_tail over walk-phase and tail steps. Tails contribute at most one lookup per step: L_tail <= P*TAIL = 2^53 deterministically. For L_walk, condition on the walk-phase exploration: given exploration law (iii), distinguished indicators are independent Bernoulli(2^-32) across distinct walk-phase points, so L_walk = sum_v m_v I_v with visit multiplicities m_v summing to T'. Conditional on no walk-phase cycling (each walk's path simple, so m_v <= P; cycling failure <= P*L_w^2/2q < 2^-18 by posit (ii), verify `assert 262144*1003842842877789998575508075287788*1003842842877789998575508075287788*2**18 < 2*q` with the full hex q from Parameters), Bernstein's inequality gives P(L_walk > 2*mu) <= exp(-3*mu/(8*P)) = exp(-3*L_w/2^35) <= e^-2^70 (mu = T'/2^32; sigma^2 <= P*mu; verify `assert 3*1003842842877789998575508075287788 >= 2**105`). The threshold fits: 2*mu + P*TAIL = T'/2^31 + P*TAIL <= 2^98 (verify `assert 263151378203355381386577988888241897472 + 262144*2**35*2**31 <= 2**129`), so {L > 2^98} implies {L_walk > 2*mu} or cycling. Total abort failure <= 2^-18 + e^-2^70, joining the union bound. Every walk runs exactly L_w + TAIL steps then halts (S deterministic); a lookup counter aborts past 2^98 total table operations and an entry counter past 2^98 entries, so table work is deterministically at most 2^98 x 2^11 on all paths including degenerate and failure paths. No cycle detection is needed: cycling is analysis-only failure mass, not an algorithmic event.

Table construction. The shared table uses separate chaining: 2^99
bucket heads (16-byte segmented pointers, zero-initialized once the
attack starts, i.e. in attack phase after Q arrives) plus entry nodes
allocated from a bump pool. A node pointer is 128-bit segmented
(34-bit segment id + 72-bit offset + 22 reserved bits, 16 bytes),
addressing 2^34 segments x 2^72 bytes = 2^106 bytes of pool, holding
2^98 160-byte nodes (2^98 x 160 = 2^105.32 bytes) with room to spare;
node access loads the segment base from a small segment table then the
offset (2 extra loads, in the per-node budget). Bucket index = key bits [32,131) (the 99 x-coordinate bits directly
above the 32-bit distinguished window; the distinguished low 32 bits are
all zero on table keys and carry no entropy). No hash function is used:
there are no salts, no re-salting, and no hash heuristic. Correctness
(exact matching, deduplication) holds by comparing full keys; occupancy
is a walk-randomness event bounded given Heuristic R; all table work is
deterministically capped (below).
Occupancy given R. Under Heuristic R, distinguished exploration points
are uniform conditioned on the distinguished window, so key bits
[32,131) are i.i.d. uniform over 2^99 buckets. With at most n <= 2^98
distinct keys (entry-cap abort below), each bucket's occupancy is
Binomial(n, 2^-99) with mean mu <= 1/2. Repeat visits to an already
stored distinguished point are lookups that hit the existing entry
(deduplication by exact key); only first visits insert, so occupancy
counts distinct keys (at most 2^98 by the abort cap), and merged-walk
revisits do not add balls to bins.
Chain-length bound (balls-in-bins over walk randomness given R, load <= 1/2). For one bucket, P(occupancy > 64)
= P(X > 128*mu) <= e^63.5/128^64 < 2^-356 by the Chernoff bound with
d = 127 (e < 2.7183 from 1/k! partial sums, log2(2.7183) < 1.4427, so
e^63.5 < 2^91.7 against 128^64 = 2^448). Union over 2^99 buckets given R:
P(any chain > 64) < 2^-257. A chain exceeding 64 aborts the attack (no
re-salt exists); this is a 2^-257 failure event joining the union bound.
Per-operation budget. Each insert/lookup walks one chain comparing exact
keys (full 256-bit x in four 64-bit words, y-parity bit); correctness
never depends on the hash. Per visited node at most 29 instructions:
segment-base load 1, node-address add 1, next-pointer load 1, six key
word loads 6, six key compares 6, parity load+compare 2, coefficient
loads 10 (deferred to match), link advance 1, null check 1. At most 64
nodes per operation gives at most 1856 chain instructions plus bucket
bit-slice addressing (absorbed in link advance): 1856 < 2^11 per
operation; total table work is at most 2^98 operations x 2^11 = 2^109
(deterministic: lookup counter aborts past 2^98 operations). Inserts prepend with atomic
compare-and-swap on the bucket head; contention is negligible at 2^-32
distinguished density and covered by the per-operation budget. No re-salt or rebuild exists. Overflow, entry-cap, and lookup-cap aborts
are the failure events above.

Union bound (no independence assumed between (A) and (B), which share
walk randomness): P(success) >= P(A) + P(B) - 1 - 2^-180 - 10^-70 - 2^-59 - 2^-18 - 2^-257 - e^-2^70
>= 0.4779 + 0.99966 - 1 - (negligible) >= 0.4775 - 10^-5 >= 0.3995 (e^-2^70 < 10^-300, absorbed; verify `assert sum(Decimal(2)**-e for e in (180,59,18,257)) + Decimal(10)**-70 < Decimal('0.00001')`).
over uniform challenges and fresh attack coins. Every returned scalar is
verified by an independent [x]P == Q check, so success means recovery.

## Cost accounting

Reference operations H: every walk-step addition is one of a Montgomery
batch of P = 262144 (3M formula each, one shared inversion per round), i.e.
a genuine batched incomplete affine addition charging H = 1; precompute
adds are likewise batched across walks (width P); verification adds are
plain reference operations. H = S + P*64*1024 + P*1024 + 2^10
= 263151378203355381386586996104944944128, where S = T' + P*TAIL with T' = T_walk + 218473,
precompute is 64 multiples x <=1024 adds each per walk (256-bit double-scalar
mults: 256 doublings + 512 adds = 768 <= 1024), walk starts are P starts x <=1024
adds each, and verification
is two 256-bit scalar mults (<= 512 adds each). The reference 6M convention
covers each batch member's 3M formula plus 3M of Montgomery batch products;
all remaining batch overhead is charged in R below, never double-counted.

Other instruction-equivalents R (C = 2000 per the protected BLS cost model), per processor-step at most 73:
partition mask/index 4, precompute table address + point load 10, two
5-limb integer adds 24 (5 ADCs each with accumulators in registers), distinguished test 3, point store-back 4, loop
control 7 (= 52 walk operations); round barrier/sync share 8 (lock-step SIMT execution: synchronization free by construction, 8 retained as margin); batch excess
13 (inverse share I/P <= 1 instruction with I <= 600M via Fermat exponentiation to p-2 for this 381-bit field under S=M (381 squarings plus under 219 multiplies), divided across P = 262144 walks: 600/262144 M = 0.0023M at the pinned 383 instructions per M < 1 instruction; Montgomery control
and addressing 12). Table insert/lookup at most 2^109 total per the Table
construction above (at most 2^98 operations x 2^11 instructions, deterministically capped).
Fixed work at most 2^73 + 2^38 + 2^102 + 2^109 (periodic coefficient reduction P*L/2^63 x <=128 instructions, precompute batch-excess and loop
control over P*64*1024 + P*1024 adds at <=8 instructions each (= 2^37+2^31), one binary-GCD inverse plus a modular multiplication for the solve (<= 2^21), start-setup sampling (3305209728 <= 2^32: two rounds of P*2*24*16 sampling + 22413312-CAS bitonic sort x 64 + dedup scan, plus 256*24*16 table scalars), decoding/verification bookkeeping and lookup counter (<= 2^21), zero-initialization of 2^99 bucket heads
at <= 4 instructions each, i.e. <= 2^101 inside the 2^102 term (2^38-term total: (2^37+2^31) + 2^32 + 2^21 + 2^21 < 2^38; verify `assert (2**34+2**28)*8 + 3305209728 + 2**21 + 2**21 < 2**38`).
Total R = 73*S + 2^109 + 2^73 + 2^38 + 2^102
= 19210051262952652568431954852700276916224.

W = H + R/2000. Claim: log2(W) <= 127.69, i.e. W/2^127 <= 2^0.69.
Proof in the correct direction: 2^0.69 = e^{0.69 ln 2} >= 1+z+z^2/2+z^3/6
with z = 0.69*0.6931 = 0.478239 (using e^w >= sum_{k=0..3} w^k/k! for
w >= 0, all terms positive, and ln 2 > 0.6931 proved below), giving
2^0.69 >= 1.61082514...; and W/2^127 = 1.60311805... for the exact H, R
above. Since 1.61082514... >= 1.60311805..., the claim holds; verify with
exact Decimal arithmetic:
`H=263151378203355381386586996104944944128`,
`R=19210051262952652568431954852700276916224`,
`W=(Decimal(H)*2000+Decimal(R))/2000`,
`z=Decimal('0.69')*Decimal('0.6931')`,
`assert 1+z+z*z/2+z*z*z/6 >= W/Decimal(2**127)`.
And ln 2 > 0.6931: ln 2 = 2*sum_{k>=0} 1/((2k+1)*3^(2k+1)) (positive
terms), whose first four terms give 0.69313475... > 0.6931; verify with:
`2*(Decimal(1)/3+Decimal(1)/81+Decimal(1)/1215+Decimal(1)/15309) >
Decimal('0.6931')`.

Memory: shared table with separate chaining: 2^99 bucket heads x 16
bytes plus at most 2^98 entry nodes x 160 budgeted bytes (actual 48-byte
BLS x-coordinate + 1 parity bit + two 40-byte unreduced accumulators + 16-byte
segmented next pointer = 145 bytes), total 2^103 + 160*2^98 = 168*2^98 = 2^105.4;
plus P walk states and one batch buffer (~2^26 bytes): log2(bytes) <=
106.0 (claimed bound kept with margin).
Preprocessing (all precomputed multiples with batch excess) <= 2^24
group-equivalents, included in H (attack phase: M_j needs Q, so all
17,448,304,640 precompute additions run after Q arrives and are ordinary
attack work, not preprocessing). True preprocessing (sampling the 128 table-scalar pairs before seeing Q: 256 scalars x 24 attempts x <=16 instructions = 98304 <= 2^17, abort on exhaustion) costs below
2^24 instructions, so preprocessing_log2 = 24 holds with margin
(manifest requires at most 2^24). Starts are sampled in attack phase (charged in R above), so no sampling sort enters preprocessing. Nonuniform
advice: none (0 bytes); all multipliers and starts are freshly sampled
and their sampling is charged in R. Single run, no restarts; all trials,
failures, sorting (start-setup bitonic sorts only, charged above; the rho table needs none), lookups, transfers, verification,
output, and cross-processor coordination are included above.

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
regime; production degeneracy is the quantified (<= 2/q) event in (A), and the
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
walk shape used here (production residues wrap exactly as the toy's, unifying the regimes). The toy's dictionary faithfully implements production matching (exact keys, dedup); production buckets are key bits (no hash function), so occupancy is analysis given R, not a second heuristic. Per the rules this is evidence only at toy scope,
not a production proof.

## Limitations

Analytic construction; no full-size execution. Success relies on the
supported multi-walk random-mapping heuristic. Instruction budgets are
conservative analytic upper bounds, not measured assembly; the 2^98-entry
bit-sliced table is an analytic device (no physical machine addresses it; the
same holds for the baseline's larger table). Best-known-attack analysis,
not novel cryptanalysis. Degeneracy is a quantified (<= 2/q) heuristic event; coefficient residues wrap per the lazy-reduction lemma. This supersedes the sequential-walk candidate in
submissions/secp256k1-ecdlp/rho39/, which charged sequential additions at
the batched reference price.

## Auxiliary-powers input (this track)

The track supplies public powers A_i = [s^i]P0 for i = 1..n. The attack
reads A_1 = [s]P0 (one point decode, charged inside the fixed-work
budget) and runs the identical parallel rho above with Q = A_1, ignoring
A_2..A_n. Since A_1 = [s]P0 with s uniform (the track's secret), recovery
of the discrete log of A_1 is exactly recovery of s; success probability
and cost are as analyzed. Total work W(n) = W is independent of n (one
extra point read regardless of n), so log2(W(2^20)) = 127.69. No Cheon
step is used; this is a conservative generic-group upper bound.
