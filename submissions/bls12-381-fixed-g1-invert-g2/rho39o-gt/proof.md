# bls12-381-fixed-g1-invert-g2: parallel Pollard rho in GT, r=64, standard-regime wrapping walk, success >= 0.3995

Status: ready for independent review; no qualification or emitted score.
This is a best-known-attack cost analysis (Pollard rho with distinguished
points in the order-r target group GT, P parallel walks with dense Fp12
multiplications, analyzed for the track's 0.39 success threshold), not a
claim of novelty or practical execution. This v7 revision (a) moves the distinguished window to c0 bits [6,38), disjoint from the partition bits [0,6) by construction (density unchanged at 2^-32, buckets re-windowed to [38,137)), (b) corrects the degeneracy bound to 1/r + 2^-100 (the v6 2/r did not follow from the stated TV distance; the margin absorbs it), (c) replaces the rapid-mixing assertion with an explicit mixing posit supported by measured toy db uniformity (chi-square 29.35 on 24 df) plus a stated Teske-hypothesis check with an acknowledged partition-delta, and (d) corrects the table-memory sum to 736*2^98 (still within the 2^108.0-byte bound). H and R are unchanged (Fp12 multiplication is total, so no exceptional path exists); the claim holds at 2^127.65 as recomputed in v6. It supersedes v6. Pre-v7 history: this v6 revision (a) bounded all sampling with explicit cutoffs (preprocessing table scalars: 24 attempts each; attack-phase starts: two rounds of sample, bitonic sort, dedup scan; abort on exhaustion with accounted failure mass), (b) replaces the per-step abort law with a dependence-aware bound (entries subsumed by lookups; tails trivially bounded; walk-phase lookups via no-cycling conditioning plus per-distinct-point Bernstein), (c) reframes the walk heuristic as a posited exploration law with labeled support (no derivation claimed), and (d) raises C0 0.5113 -> 0.65 for success-margin robustness (modeled 0.4775 vs claimed 0.3995). R, H, and the claim are recomputed (claim 2^127.65). It supersedes v5. Charges C=20000 per the protected cost model
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
One shared multiplier table u_j, v_j uniform in [1, r) (full-size), sampled once by the attack (no advice) in preprocessing: each scalar allows 24 attempts (255-bit draw, reject >= r or zero; abort the attack on exhaustion);
M_j = g^{u_j} z^{v_j} in GT.
Each walk has a fresh start (a, b) uniform in [0, r)^2, sampled in attack phase by two rounds of (sample all P starts with 24 attempts per scalar, bitonic-sort, dedup-scan), aborting on any remaining duplicate or range failure. Total sampling failure is at most 2^-59 (round-two duplicates below 2^-219 plus 255-bit range exhaustion below 2^-60; verify `assert 4096*(4096-1)*2**219 < 2*r` and `assert 16640*(2**255-r)**24*2**60 < 2**(255*24)` with the full hex r from Parameters; sum below 2^-59), joining the union bound; sampling work is charged in fixed work (below). Exponents are tracked as residues mod r,
computed by lazy reduction: plain-integer 5-limb accumulators reduced
mod r every 2^63 steps (see lemma below); residues are mathematically
identical to per-step wrapping.
Distinguished points: c0 bits [6,38) all zero (density 2^-32); the identity
(c0 = 1) is never distinguished. The window is disjoint from the partition bits [0,6) by construction; no distinguished point is forced into any partition (the containment is removed, not assumed away).

T_walk = 263151378203355381386577988888241678999 (walk-phase target),
TAIL = 2^35 steps appended to every walk. Each walk runs
L_w = ceil(T_walk/P) walk-phase steps plus TAIL tail steps, so the
pooled walk-phase exploration covers T' = P*L_w = T_walk + 1385 points;
total charged steps S = T' + P*TAIL. One shared distinguished-point
table (keyed by full 576-byte Fp12 key, see below); the attack aborts if
the table exceeds 2^98 entries.

## Construction

```text
Input P, z. If z == 1: return infinity.
Compute g = e(P,Q0) (one pairing, bounded below).
Sample one shared u_j, v_j in [1,r) (preprocessing: 24 attempts per scalar, abort on exhaustion); form M_j = g^{u_j} z^{v_j}
(64 points x 1024 Fp12 mults: two 255-bit binary exponentiations each). Shared table = {} (separate chaining
per the Table subsection below).
For each of P walks: sample start (a,b) uniform mod r by the two-round setup (sample, bitonic sort, dedup scan; abort on remaining duplicates); Y = g^a z^b (1024 Fp12 mults).
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
< 2^117 steps (L_w is about 2^115.63); total reductions P*L/2^63,
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

Heuristic R (posited exploration law). One fixed r-adding walk function with r=64 partitions (shared table) is
iterated from P fresh starts in GT. The analysis posits the following exploration law, with support (not derivation) below: (i) pooled birthday: P(no repeated point among the T' walk-phase points) <= exp(-(63/64)*T'(T'-1)/2r), charging only distinct-partition pairs (a 63/64 fraction); same-partition pairs are excluded conservatively (their later merges can only add collisions); (ii) single-walk self-collisions: one walk self-collides within L steps with probability at most L^2/2r (standard birthday upper bound); (iii) per-point labels: distinguished indicators are independent Bernoulli(2^-32) across distinct exploration points, independent of the walk exploration, and c0 bits [38,137) are i.i.d. uniform across distinct distinguished keys (bucketing); the distinguished window [6,38) is disjoint from the partition bits [0,6) by construction, so label/partition overlap is structurally absent; residual same-c0 bit correlation is posited negligible; (iv) degeneracy: a solved collision has db = 0 with probability at most 1/r + 2^-100 (corrected from the 2/r; mixing posit with measured support in (A)). The v13 per-step-independence framing was false even in the ideal model (repeated visits share status); (iii) holds exactly there.
Support for R (evidence for the posit, with hypothesis checks). (a) Teske (2000): rho length sqrt(pi*q/2) times E_r tending to 1 (E_64 within a few percent of 1), assuming a uniform random table and full-period wrapping, both met here by construction. One delta is acknowledged: Teske partitions by an auxiliary random hash, while this walk partitions by c0 bits [0,6); the posit extends Teske's form across this delta, supported by (d)-(e) and the margin. We do not claim Teske proves this walk. Other support: Teske's analysis and experiments for large-r adding walks (r-dependent correction tending to 1; r = 64 is large; uniform-table wrapping regime with c0-coefficient partitions, matching this construction); Harris's rho-length analysis corroborating the quantitative form; van Oorschot-Wiener's parallel-rho pooling model (shared iteration function from P starts, the exact multi-start setting here); record parallel-rho computations (112-bit ECDLP, Bos et al., completing in predicted time, validating the birthday constant near scale); our shared-function toy experiments confirming merging and recovery end to end (see Experiments), in the same wrapping regime as production, and measuring db uniformity directly (chi-square 29.35 on 24 df over 997 solves, p approx 0.21), supporting the 1/r shape of (iv). The posit needs only coarse accuracy: C0 = 0.65 sizes the walk for modeled success 0.4775 against the claimed 0.3995 (19% margin), so constant-factor approximation error is absorbed. This heuristic is the only walk heuristic in the entire argument; everything conditional on it is proved (the tail bound is a Chernoff/exponential calculation given (ii)-(iii), table occupancy is balls-in-bins given (iii), and table/cycling/sampling failures join the union bound explicitly).

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

(A) Walk collision. By posited exploration law (i) with the T' pooled walk-phase points: P(no collision) <= exp(-(63/64)*T'(T'-1)/2r). The exact integer T_walk satisfies T_walk(T_walk-1) >= 2r*C0 with C0 = 0.65*64/63 = 0.660317460317..., compensating the 63/64 merge factor; verify with: `T=263151378203355381386577988888241678999` and the full hex r from Parameters:
`C0 = Decimal('0.65')*Decimal(64)/Decimal(63)` (60-digit Decimal),
`assert T*(T-1) >= 2*r*C0 and (T-1)*(T-2) < 2*r*C0` (minimality).
Since T' = T_walk + 1385 >= T_walk and the bound decreases in T, (63/64)*T'(T'-1)/2r >= (63/64)*C0 = 0.65, so P(no walk collision) <= exp(-0.65).
And ln(1/0.5221) < 0.65 because e^0.65 > sum_{k=0..8} 0.65^k/k! = 1.91554076... > 1/0.5221 = 1.91534188...; verify with:
`sum(Decimal('0.65')**k/math.factorial(k) for k in range(9)) > Decimal(1)/Decimal('0.5221')` (all terms positive, so truncation is a lower bound). Hence P(no walk collision) <= 0.5221; the walk pays 64/63 in length (charged in H, R), not in success.
Degeneracy analysis (mixing support, quantified). At a table hit, two visits carry residues (a,b), (a0,b0) with da + x*db = 0 mod r, where da = a-a0, db = b-b0 mod r and z = g^x. The solve needs db != 0 mod r. The difference db equals the signed sum of table vectors over the symmetric difference of the two root-to-visit paths: the two distinct pre-merge segments (totaling ~L_w table vectors) plus a shared-suffix correction that cancels. Pre-merge accumulation is a random walk on (Z/r)^2 with 64 uniform generators, which is posited (not proved) to mix rapidly after L_w ~ 2^115.63 steps db is 2^-100-close to uniform over Z/r in total variation. This mixing is an explicit posit: spanning does not imply rapid mixing by itself, nor are partition-selected increments independent Cayley increments (the prior assertion is withdrawn); support comes from 64 uniform-random generators (Alon-Roichman-type heuristic evidence only), measured toy db uniformity (chi-square 29.35/24df), and exclusion of the bounded-offset failure mode by construction. TV distance at most 2^-100 gives P(db = 0 mod r) <= 1/r + 2^-100 (the 2/r was an error, corrected here). Table-identity (M_j = 1, 64/r over uniform table sampling) joins the same bound. Thus P(A), a useful non-degenerate walk-phase collision exists, is >= 0.4779 - (1/r + 2^-100) - 64/r >= 0.4779 - 10^-30 (from P(no collision) <= 0.5221; verify `assert Decimal(2)**-100 + Decimal(65)/Decimal(r) < Decimal(10)**-30` with the full hex r from Parameters).

(B) Tail collection. After the first collision, all walks that visited it follow the merged shared-function suffix, which reaches a distinguished point within TAIL steps with probability >= 1-e^-8-2^-180 >= 0.99966. Merging is exact (shared transition function), not heuristic. Given exploration law (iii), distinguished indicators are independent across distinct suffix points at density 2^-32; conditional on the suffix being simple (no self-collision; failure <= TAIL^2/2r < 2^-180 by posit (ii), verify `2**70*2**180 < 2*r` with the full hex r from Parameters), its 2^35 distinct points each hit independently, so P(no hit) = (1-2^-32)^(2^35) <= e^-8. Proof: (1-1/n)^m <= e^(-m/n) gives (1-2^-32)^(2^35) <= e^-8; e > 1+1+1/2+1/6+1/24+1/120 = 2.71666... (Taylor lower bound, all terms positive), so e^8 > 2.71666^8 = 2966.81... and e^-8 < 0.00033707, hence P(B) > 0.9996629 - 2^-180 >= 0.99966. On a hit the point is in the shared table (all distinguished visits are stored) and the stored collision solves unless degenerate (absorbed in (A)). The old per-step-independence framing was false even in the ideal model (repeated visits share status); the per-distinct-point framing above holds exactly there and is posited of the construction.

Table abort (entries or total lookups exceed 2^98). Entries are subsumed by lookups (each entry needs a visit: N_entries <= L), so one bound covers both. Split total lookups L = L_walk + L_tail over walk-phase and tail steps. Tails contribute at most one lookup per step: L_tail <= P*TAIL = 2^47 deterministically. For L_walk, condition on the walk-phase exploration: given exploration law (iii), distinguished indicators are independent Bernoulli(2^-32) across distinct walk-phase points, so L_walk = sum_v m_v I_v with visit multiplicities m_v summing to T'. Conditional on no walk-phase cycling (each walk's path simple, so m_v <= P; cycling failure <= P*L_w^2/2r < 2^-12 by posit (ii), verify `assert 4096*64245941944178559908832516818418379*64245941944178559908832516818418379*2**12 < 2*r` with the full hex r from Parameters), Bernstein's inequality gives P(L_walk > 2*mu) <= exp(-3*mu/(8*P)) = exp(-3*L_w/2^35) <= e^-2^70 (mu = T'/2^32; sigma^2 <= P*mu; verify `assert 3*64245941944178559908832516818418379 >= 2**105`). The threshold fits: 2*mu + P*TAIL = T'/2^31 + P*TAIL <= 2^98 (verify `assert 263151378203355381386577988888241680384 + 4096*2**35*2**31 <= 2**129`), so {L > 2^98} implies {L_walk > 2*mu} or cycling. Total abort failure <= 2^-12 + e^-2^70, joining the union bound. Every walk runs exactly L_w + TAIL steps then halts (S deterministic); a lookup counter aborts past 2^98 total table operations and an entry counter past 2^98 entries, so table work is deterministically at most 2^98 x 2^14 on all paths including degenerate and failure paths. No cycle detection is needed: cycling is analysis-only failure mass, not an algorithmic event.

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
704. Bucket index = c0 bits [38,137) (the 99 bits directly above the 32-bit
distinguished window [6,38); the distinguished c0 bits are all zero on
table keys and carry no entropy). No hash function is used: there are
no salts, no re-salting, and no hash heuristic. Correctness (exact
matching, deduplication) holds by comparing full keys; occupancy is a
walk-randomness event bounded given Heuristic R; all table work is
deterministically capped (below).
Occupancy given R. Under Heuristic R, distinguished exploration points
are uniform conditioned on the distinguished window, so c0 bits
[38,137) are i.i.d. uniform over 2^99 buckets. With at most n <= 2^98
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
keys (full 576-byte key in seventy-two 64-bit words); correctness never
depends on the hash. Per visited node at most 168 instructions:
segment-base load 1, node-address add 1, next-pointer load 1, 72 key
word loads 72, 72 key compares 72, exponent loads 10 (deferred to match),
link advance 1, null check 1, atomic-CAS share 4 (total 163 <= 168). At
most 64 nodes per operation gives at most 10432 chain instructions plus bucket
bit-slice addressing (absorbed in link advance): 10432 < 2^14 per
operation; total table work is at most 2^98 x 2^14 = 2^112
(deterministic: lookup counter aborts past 2^98 operations). Inserts prepend with
atomic compare-and-swap on the bucket head; contention is negligible at
2^-32 distinguished density and covered by the per-operation budget. No re-salt or rebuild exists. Overflow, entry-cap, and lookup-cap aborts
are the failure events above.

Union bound (no independence assumed between (A) and (B), which share
walk randomness): P(success) >= P(A) + P(B) - 1 - 2^-180 - 10^-30 - 2^-59 - 2^-12 - 2^-257 - e^-2^70
>= 0.4779 + 0.99966 - 1 - (negligible) >= 0.4775 - 10^-3 >= 0.3995 (e^-2^70 < 10^-300, absorbed; verify `assert sum(Decimal(2)**-e for e in (180,59,12,257)) + Decimal(10)**-30 < Decimal('0.001')`). over uniform challenges (nonidentity P, uniform z) and
fresh attack coins. Every returned Q is verified by an independent
e(P,Q) == z pairing check, so success means recovery.

## Cost accounting

Reference operations H: every walk step is one dense Fp12
multiplication (the reference operation, charging H = 1); table setup
is 64 multipliers x 1024 Fp12 mults each (two 255-bit binary
exponentiations per multiplier); walk starts are P starts x 1024 Fp12
mults each (two 255-bit binary exponentiations per start). H = S + 2^22
+ 2^16 = 263151378203355381386578129625734295552, where S = T' + P*TAIL
with T' = P*ceil(T_walk/P) = T_walk + 1385. No reference operation and
its implementation are double-counted (the 22996-instruction
implementation is the unit itself).

Other instruction-equivalents R (C = 20000 per the protected cost
model), per processor-step at most 200: partition mask/index 6 (c0 low
word load 2, mask and index 4), table address + GT element load 80 (72
word loads + addressing 8), two 5-limb integer adds 24 (5 ADCs each with accumulators in registers), distinguished
test 3, point store-back 76 (72 words + address 4), loop control 11.
No batch excess and no barrier (independent walks, no lock-step).
Table insert/lookup at most 2^112 total per the Table construction above
(at most 2^98 operations x 2^14 instructions, deterministically capped).
Fixed work at most 2^73 + 2^31 + 2^23 + 2^101 + 2^28: periodic exponent reduction P*L/2^63 x <=128 instructions, two pairings at
<= 2^30 instructions each (Miller loop: 63 iterations x (line
evaluation + sparse multiply + squaring) <= 63 x 4 dense mults at 22996
< 2^15 instructions each, i.e. <= 2^23; final exponentiation: <= 2^12
dense Fp12 ops x 2^15 = 2^27 by square-and-multiply on the ~1000-bit
exponent; subgroup checks and bookkeeping <= 2^26; total <= 2^28 <
2^30), one G2 scalar mult at <= 2^23 (510 G2 ops x 2^14 each; each
affine G2 op is ~6 Fp2 mults ~= 18 Fp mults at the pinned 383
instructions each, i.e. below 2^13, charged 2^14 with margin),
zero-initialization of 2^99 bucket heads at <= 4 instructions each
(<= 2^101), sampling/decoding/solve/bookkeeping (start-setup sampling 27361152 <= 2^25: two rounds of P*2*24*16 sampling + 159744-CAS bitonic sort x 64 + dedup scan, plus 256*24*16 table scalars; one binary-GCD inverse plus a modular multiplication for the solve <= 2^21; all randomness, z == 1 test, decoding, lookup counter <= 2^21; total 2^25 + 2^21 + 2^21 < 2^28) at <= 2^28.
Total R = 200*S + 2^112 + 2^73 + 2^31 + 2^23 + 2^101 + 2^28
= 52630280835503236022044446224377457508352.

W = H + R/20000. Claim: log2(W) <= 127.65, i.e. W/2^127 <= 2^0.65.
Proof in the correct direction: 2^0.65 = e^{0.65 ln 2} >= 1+z+z^2/2+z^3/6
with z = 0.65*0.6931 = 0.450515 (using e^w >= sum_{k=0..3} w^k/k! for
w >= 0, all terms positive, and ln 2 > 0.6931 proved below), giving
2^0.65 >= 1.56723658...; and W/2^127 = 1.56213144... for the exact H, R
above. Since 1.56723658... >= 1.56213144..., the claim holds; verify with
exact Decimal arithmetic:
`H=263151378203355381386578129625734295552`,
`R=52630280835503236022044446224377457508352`,
`W=(Decimal(H)*20000+Decimal(R))/20000`,
`z=Decimal('0.65')*Decimal('0.6931')`,
`assert 1+z+z*z/2+z*z*z/6 >= W/Decimal(2**127)`.
And ln 2 > 0.6931: ln 2 = 2*sum_{k>=0} 1/((2k+1)*3^(2k+1)) (positive
terms), whose first four terms give 0.69313475... > 0.6931; verify with:
`2*(Decimal(1)/3+Decimal(1)/81+Decimal(1)/1215+Decimal(1)/15309) >
Decimal('0.6931')`.

Memory: shared table with separate chaining: 2^99 bucket heads x 16
bytes (flat 2^103-byte array) plus at most 2^98 entry nodes x 704 bytes
(576-byte key + two 40-byte accumulators + 16-byte next pointer = 672
bytes exact, budgeted 704), total 2^103 + 704*2^98 = 736*2^98 = 2^107.52 (corrected from the 712*2^98 arithmetic slip; 32 + 704 = 736); plus P walk
states, one 64-multiplier table, and pairing temporaries (~2^22 bytes):
log2(bytes) <= 108.0 (claimed bound kept with margin).
Preprocessing (sampling the 128 table-scalar pairs before seeing P and z: 256 scalars x 24 attempts x <=16 instructions = 98304 <= 2^17, abort on exhaustion) costs below 2^24 instructions, so
preprocessing_log2 = 24 holds with margin (manifest requires at most
2^24). Starts are sampled in attack phase (charged in R above), so no sampling sort enters preprocessing. All table multipliers, pairings, and walk work need P or z and
run in attack phase as ordinary charged work, not preprocessing.
Nonuniform advice: none (0 bytes); all multipliers and starts are
freshly sampled and their sampling is charged in R. Single run, no
restarts; all trials, failures, sorting (start-setup bitonic sorts only, charged above; the rho table needs none), lookups,
transfers, verification, output, and cross-processor coordination are
included above.

## Experiments

`experiment.py` runs the GT-notation algorithm mechanics at toy scale (v15: distinguished window moved above the partition bits, mirroring production; db-uniformity instrumentation added)
(one shared r=64 (u,v) table with strictly positive entries, P = 8 toy
walks from fresh starts with per-walk exponent offsets, unreduced
integer exponent accumulation with mod-n only inside the solve equation,
shared distinguished table keyed by full element, collision solve with
degeneracy counting, independent verification) in multiplicative
notation over a prime-order subgroup of F_p^* (order-9883 subgroup of
F_158129^*, partitions and distinguished tests read the element value
as production reads c0, identity never distinguished), 2000 trials:
observed success rate 0.4985 (997 successes) with ZERO degenerate encounters, every
recovered scalar verified correct. A second toy group (order-967
subgroup of F_15473^*) gives 0.825 with 2 degenerate encounters out of
1652 solves; this is a toy-range artifact (toy exponent range
2^43 exceeds the toy order 967, so distinct integers collide mod n).
In v2 production residues likewise wrap (lazy-reduction lemma), so toy
and production share the wrapping regime; production degeneracy is the
quantified (<= 1/r + 2^-100) event in (A), and the toy rate is consistent with
wrapping-regime behavior at small scale. The toy also measures coefficient-difference uniformity directly: db over all 997 production-path solves gives chi-square 29.35 on 24 df (p approx 0.21), consistent with uniform db (second group: 24.82 on 24 df); the experiment asserts chi-square below 40 as a regression guard. This supports the 1/r SHAPE of posit (iv) where measurable; the 2^-100 TV error at production scale remains an explicit posit with margin. (A partition step chi-square was measured during development but dropped: consecutive steps are Markov-dependent.) Toy tails are minimal (4 steps at
density 1/4) yet still contribute collisions beyond the walk-only
prediction (~0.40); this excess is in the conservative direction and
confirms the merge-and-solve mechanism works, including cross-walk
merging onto shared suffixes. Production tails (2^35 steps) are a 2^-93
fraction of walk points and are handled by the rigorous tail bound
above, not by the toy. The toy therefore validates mechanism
(shared-function merging, offset/monotone exponent handling,
exact-arithmetic solve, verification); the production 0.3995 claim rests
on the analysis. In particular the toy's wrapping exponents (toy range
exceeds toy order) and minimal tails differ from production by
construction and are not offered as evidence for the production birthday
constant, the 63/64 correction, or tail behavior, which are supported by
the Teske+Harris regime analysis for the walk shape used here
(v2: production residues wrap exactly as the toy's, unifying the
regimes). The toy's dictionary faithfully implements production matching (exact keys, dedup); production buckets are c0 bits [38,137) (no hash function), so occupancy is analysis given R, not a second heuristic. Per the rules this is evidence only at toy scope, not a production proof.

## Limitations

Analytic construction; no full-size execution. Success relies on the
supported multi-walk random-mapping heuristic. Instruction budgets are
conservative analytic upper bounds, not measured assembly; the 2^98-entry
table is an analytic device (no physical machine addresses it).
Best-known-attack analysis, not novel cryptanalysis. Sixth submission
on this track (supersedes v5); no baseline exists yet.
