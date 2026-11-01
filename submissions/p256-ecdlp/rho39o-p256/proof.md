# p256-ecdlp: batched parallel Pollard rho, r=64, posited exploration law, success >= 0.3995

Status: ready for independent review; no qualification or emitted score.
This is a best-known-attack cost analysis (Pollard rho with distinguished
points, P parallel walks with Montgomery-batched affine additions,
analyzed for the track's 0.39 success threshold), not a claim of novelty
or practical execution. This v15 revision (a) moves the distinguished window to Montgomery-form x-coordinate bits [6,38), disjoint from the partition bits [0,6) by construction (the v14 window contained the partition bits; the overlap is now structurally removed; density unchanged at 2^-32, buckets re-windowed to [38,137)), (b) corrects the degeneracy bound to 1/q + 2^-100 (the v14 2/q did not follow from the stated TV distance; the success margin absorbs the corrected term), (c) replaces the rapid-mixing assertion with an explicit mixing posit supported by measured toy coefficient-difference uniformity (db chi-square 19.89 on 24 df over 1069 production-path solves) plus a stated Teske-hypothesis check with an acknowledged partition-delta, (d) specifies the complete-addition exceptional path with instruction counts, batch-safe zero-denominator detection, an explicit E* = 2^20 work cap with Chernoff failure mass, and its R term. P-256 field arithmetic is unchanged (SOS+REDC with explicit carries, six modular subtractions, per-step 1969, H = 0); R gains only the tiny exceptional-path cap terms and the claim holds at 2^129.20 (recomputed below). It supersedes v14. Charges C=1000 per the protected p256 cost model. The normative inputs/outputs are those of the
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
One shared multiplier table u_j, v_j uniform in [1, q) (full-size), sampled once by the attack (no advice) in preprocessing: each scalar allows 24 attempts (256-bit draw, reject >= q or zero; abort the attack on exhaustion). Each walk has a fresh start (a, b) uniform in [0, q)^2, sampled in attack phase by two rounds of (sample all P starts with 24 attempts per scalar, bitonic-sort, dedup-scan), aborting on any remaining duplicate or range failure. Total sampling failure is at most 2^-220 (round-two duplicate pairs below 2^-221 plus range exhaustion below 2^-700; verify `assert 262144*(262144-1)*2**221 < 2*q` and `assert 1048832*(2**256-q+1)**24*2**700 < 2**(256*24)` with the full hex q from Parameters), joining the union bound; sampling work is charged in fixed work (below). Coefficients are tracked as residues mod q,
computed by lazy reduction: plain-integer 5-limb accumulators reduced
mod q every 2^63 steps (see lemma below); residues are mathematically
identical to per-step wrapping.
Distinguished points: Montgomery-form x-coordinate bits [6,38) all zero (density 2^-32); infinity never is. The window is disjoint from the partition bits [0,6) by construction (multiplication by R permutes F_p, so Montgomery-form bits carry the same heuristic shape with an exact bijection for keys); no distinguished point is forced into any partition.

T_walk = 391048688201955837422046370934246227631 (combined walk-phase target),
TAIL = 2^35 steps appended to every walk. Each walk runs
L_w = ceil(T_walk/P) walk-phase steps plus TAIL tail steps, so the
pooled walk-phase exploration covers T' = P*L_w = T_walk + 115025 points;
total charged steps S = T' + P*TAIL. One shared
distinguished-point table (keyed by x AND y-parity, see below); the attack
aborts if the table exceeds 2^98 entries.

## Construction

```text
Input Q. Sample one shared u_j, v_j in [1,q) (preprocessing: 24 attempts per scalar, abort on exhaustion); form
M_j = [u_j]P + [v_j]Q (batched Montgomery-domain mults below). Shared table = {} (separate
chaining per the Table subsection below).
For each of P walks: sample start (a,b) uniform mod q by the two-round setup (sample, bitonic sort, dedup scan; abort on remaining duplicates).
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
< 2^117 steps (L_w is about 2^110.20); total reductions P*L/2^63,
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
Exceptional path (complete addition, specified). Each round's P additions are normally incomplete Montgomery-domain adds inside one invert-batch. An addition is exceptional exactly when its affine denominator vanishes. Detection is batch-safe and costs 1 instruction per round: a single zero-check on the batch's running denominator product detects any zero denominator without per-step comparisons. On a zero product the round bisects to the culprit walks (at most log2(P) = 18 re-batched sub-rounds, charged to the fallback below) and each culprit walk is redone by the scalar complete routine: convert the two points out of Montgomery domain (2 modular multiplications by R^-1 mod p at F_p256 <= 304 each), then if P1 is infinity return P2; elif P2 is infinity return P1; elif x1 == x2: if y1 + y2 == 0 mod p return infinity else the doubling formula; else the incomplete affine formula; convert back (1 multiplication); at most 1024 instructions per redone walk including the bisection share. At most E* = 2^20 exceptional walks are ever redone; past E* the attack aborts (deterministic work cap). Per-step exceptional probability is at most 4/q given Heuristic R uniformity; P(Binomial(S, 4/q) > 2^20) <= exp(-2^20) by the Chernoff bound with mean 4S/q < 2^-124 (verify `assert 4*(391048688201955837422046370934246342656 + 9007199254740992)*2**124 < q` with S = T' + P*TAIL and the full hex q from Parameters), so the cap-abort failure mass is below e^-2^20, joining the union bound. The R ledger gains E* x 1024 = 2^30 for the fallback cap plus one zero-check per round over L_w + TAIL <= 2^111 rounds (below). Table-identity events
(M_j = infinity, probability 64/q over uniform table sampling) join the
union bound in (A).

## Success analysis

Heuristic R (posited exploration law). One fixed r-adding walk function with r=64 partitions (shared table) is
iterated from P fresh starts. The analysis posits the following exploration law, with support (not derivation) below: (i) pooled birthday: P(no repeated point among the T' walk-phase points) <= exp(-(63/64)*T'(T'-1)/2q), charging only distinct-partition pairs (a 63/64 fraction); same-partition pairs are excluded conservatively (their later merges can only add collisions); (ii) single-walk self-collisions: one walk self-collides within L steps with probability at most L^2/2q (standard birthday upper bound); (iii) per-point labels: distinguished indicators are independent Bernoulli(2^-32) across distinct exploration points, independent of the walk exploration, and key bits [38,137) are i.i.d. uniform across distinct distinguished keys (bucketing); the distinguished window [6,38) is disjoint from the partition bits [0,6) by construction, so label/partition overlap is structurally absent; residual dependence is only same-x bit correlation, posited negligible; (iv) degeneracy: a solved collision has db = 0 with probability at most 1/q + 2^-100 (corrected from the v14 2/q, which did not follow from the stated TV distance; mixing posit with measured support in (A)). The v13 per-step-independence framing was false even in the ideal model (repeated visits share status); (iii) holds exactly there.
Support for R (evidence for the posit, with hypothesis checks). (a) Teske (2000): expected rho length is sqrt(pi*q/2) times an r-dependent correction E_r tending to 1 (E_64 within a few percent of 1), assuming a uniform random precomputed table and full-period wrapping, both met here by construction (full-size uniform u_j, v_j; lazy-reduction wrapping lemma). One delta is acknowledged: Teske partitions by an auxiliary random hash, while this walk partitions by Montgomery-form x-coordinate bits [0,6) (a permutation of F_p, same heuristic shape); the posit extends Teske's quantitative form across this delta, supported by (d)-(e) and the margin below. We do not claim Teske proves this walk. (b) Harris's rho-length analysis corroborates the quantitative birthday form for adding walks. (c) van Oorschot-Wiener's parallel-rho pooling model (shared iteration function from P starts) is the exact multi-start setting here. (d) Record parallel-rho computations (112-bit ECDLP, Bos et al., completing in predicted time) validate the birthday constant near scale. (e) Our shared-function toy experiments (see Experiments) confirm merging and recovery end to end in the same wrapping regime, and measure coefficient-difference uniformity directly: db over 1069 production-path toy solves gives chi-square 19.89 on 24 df (p approx 0.70), consistent with uniform db, supporting the 1/q shape of (iv). The posit needs only coarse accuracy: C0 = 0.65 sizes the walk for modeled success 0.4775 against the claimed 0.3995 (19% margin), so constant-factor approximation error is absorbed. This heuristic is the only walk heuristic in the entire argument; everything conditional on it is proved (the tail bound is a Chernoff/exponential calculation given (ii)-(iii), table occupancy is balls-in-bins given (iii), and table/cycling/sampling failures join the union bound explicitly).

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

(A) Walk collision. By posited exploration law (i) with the T' pooled walk-phase points: P(no collision) <= exp(-(63/64)*T'(T'-1)/2q). The exact integer T_walk satisfies T_walk(T_walk-1) >= 2q*C0 with C0 = 0.65*64/63 = 0.660317460317..., compensating the 63/64 merge factor; verify with: `T=391048688201955837422046370934246227631` and the full hex q from Parameters:
`C0 = Decimal('0.65')*Decimal(64)/Decimal(63)` (60-digit Decimal),
`assert T*(T-1) >= 2*q*C0 and (T-1)*(T-2) < 2*q*C0` (minimality).
Since T' = T_walk + 115025 >= T_walk and the bound decreases in T, (63/64)*T'(T'-1)/2q >= (63/64)*C0 = 0.65, so P(no walk collision) <= exp(-0.65).
And ln(1/0.5221) < 0.65 because e^0.65 > sum_{k=0..8} 0.65^k/k! = 1.91554076... > 1/0.5221 = 1.91534188...; verify with:
`sum(Decimal('0.65')**k/math.factorial(k) for k in range(9)) > Decimal(1)/Decimal('0.5221')` (all terms positive, so truncation is a lower bound). Hence P(no walk collision) <= 0.5221; the walk pays 64/63 in length (charged in H, R), not in success.
Degeneracy analysis (mixing posit, quantified honestly). At a table hit, two visits carry residues (a,b), (a0,b0) with da + x*db = 0 mod q, where da = a-a0, db = b-b0 mod q. The solve needs db != 0 mod q. The difference db equals the signed sum of table vectors over the symmetric difference of the two root-to-visit paths: the two distinct pre-merge segments (totaling ~L_w table vectors) plus a shared-suffix correction that cancels. We posit db is 2^-100-close to uniform over Z/q in total variation after L_w ~ 2^110.20 steps. This is a posit, not a theorem: we do not claim spanning implies rapid mixing (it does not by itself), nor that partition-selected increments are independent Cayley increments (they are not; the v14 "standard rapid mixing" assertion is withdrawn). Support: (1) 64 uniform-random generators (Alon-Roichman-type heuristics suggest expansion, though 64 is below the C*log|G| the theorem needs, so heuristic evidence only); (2) measured toy uniformity: db over 1069 production-path toy solves gives chi-square 19.89 on 24 df (p approx 0.70), consistent with uniform db, i.e. the 1/q shape holds where measurable; (3) the bounded-offset failure mode is excluded by construction (full-size uniform table, wrapping residues). A TV distance of at most 2^-100 from uniform gives P(db = 0 mod q) <= 1/q + 2^-100 (the v14 2/q was an error and is corrected here). Table-identity (M_j = infinity, 64/q over uniform table sampling) joins the same bound. Thus P(A), a useful non-degenerate walk-phase collision exists, is >= 0.4779 - (1/q + 2^-100) - 64/q >= 0.4779 - 10^-30 (from P(no collision) <= 0.5221; verify `assert Decimal(2)**-100 + Decimal(65)/Decimal(q) < Decimal(10)**-30` with the full hex q from Parameters).

(B) Tail collection. After the first collision, all walks that visited it follow the merged shared-function suffix, which reaches a distinguished point within TAIL steps with probability >= 1-e^-8-2^-180 >= 0.99966. Merging is exact (shared transition function), not heuristic. Given exploration law (iii), distinguished indicators are independent across distinct suffix points at density 2^-32; conditional on the suffix being simple (no self-collision; failure <= TAIL^2/2q < 2^-180 by posit (ii), verify `2**70*2**180 < 2*q` with the full hex q from Parameters), its 2^35 distinct points each hit independently, so P(no hit) = (1-2^-32)^(2^35) <= e^-8. Proof: (1-1/n)^m <= e^(-m/n) gives (1-2^-32)^(2^35) <= e^-8; e > 1+1+1/2+1/6+1/24+1/120 = 2.71666... (Taylor lower bound, all terms positive), so e^8 > 2.71666^8 = 2966.81... and e^-8 < 0.00033707, hence P(B) > 0.9996629 - 2^-180 >= 0.99966. On a hit the point is in the shared table (all distinguished visits are stored) and the stored collision solves unless degenerate (absorbed in (A)). The old per-step-independence framing was false even in the ideal model (repeated visits share status); the per-distinct-point framing above holds exactly there and is posited of the construction.

Table abort (entries or total lookups exceed 2^98). Entries are subsumed by lookups (each entry needs a visit: N_entries <= L), so one bound covers both. Split total lookups L = L_walk + L_tail over walk-phase and tail steps. Tails contribute at most one lookup per step: L_tail <= P*TAIL = 2^53 deterministically. For L_walk, condition on the walk-phase exploration: given exploration law (iii), distinguished indicators are independent Bernoulli(2^-32) across distinct walk-phase points, so L_walk = sum_v m_v I_v with visit multiplicities m_v summing to T'. Conditional on no walk-phase cycling (each walk's path simple, so m_v <= P; cycling failure <= P*L_w^2/2q < 2^-18 by posit (ii), verify `assert 262144*1491732361610244130790887340294824*1491732361610244130790887340294824*2**18 < 2*q` with the full hex q from Parameters), Bernstein's inequality gives P(L_walk > 2*mu) <= exp(-3*mu/(8*P)) = exp(-3*L_w/2^35) <= e^-2^70 (mu = T'/2^32; sigma^2 <= P*mu; verify `assert 3*1491732361610244130790887340294824 >= 2**105`). The threshold fits: 2*mu + P*TAIL = T'/2^31 + P*TAIL <= 2^98 (verify `assert 391048688201955837422046370934246342656 + 262144*2**35*2**31 <= 2**129`), so {L > 2^98} implies {L_walk > 2*mu} or cycling. Total abort failure <= 2^-18 + e^-2^70, joining the union bound. Every walk runs exactly L_w + TAIL steps then halts (S deterministic); a lookup counter aborts past 2^98 total table operations and an entry counter past 2^98 entries, so table work is deterministically at most 2^98 x 2^11 on all paths including degenerate and failure paths. No cycle detection is needed: cycling is analysis-only failure mass, not an algorithmic event.

Table construction. The shared table uses separate chaining: 2^99
bucket heads (16-byte segmented pointers, zero-initialized once the
attack starts, i.e. in attack phase after Q arrives) plus entry nodes
allocated from a bump pool. A node pointer is 128-bit segmented
(34-bit segment id + 72-bit offset + 22 reserved bits, 16 bytes),
addressing 2^34 segments x 2^72 bytes = 2^106 bytes of pool, holding
2^98 144-byte nodes (2^98 x 144 = 2^105.17 bytes) with room to spare;
node access loads the segment base from a small segment table then the
offset (2 extra loads, in the per-node budget). Bucket index = key bits [38,137) (the 99 Montgomery-form x-coordinate
bits directly above the 32-bit distinguished window [6,38); the distinguished
bits are all zero on table keys and carry no entropy; the bucket range is disjoint
from both the partition bits [0,6) and the distinguished window). No hash
function is used: there are no salts, no re-salting, and no hash
heuristic. Correctness (exact matching, deduplication) holds by comparing
full keys; occupancy is a walk-randomness event bounded given Heuristic
T; all table work is deterministically capped (below).
Occupancy given R. Under Heuristic R, distinguished exploration points
are uniform conditioned on the distinguished window, so key bits
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
keys (full 256-bit x in four 64-bit words, y-parity bit); correctness
never depends on the hash. Per visited node at most 25 instructions:
segment-base load 1, node-address add 1, next-pointer load 1, four key
word loads 4, four key compares 4, parity load+compare 2, coefficient
loads 10 (deferred to match), link advance 1, null check 1. At most 64
nodes per operation gives at most 1600 chain instructions plus bucket
bit-slice addressing (absorbed in link advance): 1600 < 2^11 per
operation; total table work is at most 2^98 operations x 2^11 = 2^109
(deterministic: lookup counter aborts past 2^98 operations). Inserts prepend with atomic
compare-and-swap on the bucket head; contention is negligible at 2^-32
distinguished density and covered by the per-operation budget. No re-salt or rebuild exists. Overflow, entry-cap, and lookup-cap aborts
are the failure events above.

Union bound (no independence assumed between (A) and (B), which share
walk randomness): P(success) >= P(A) + P(B) - 1 - 2^-180 - 10^-30 - 2^-220 - 2^-18 - 2^-257 - e^-2^70 - e^-2^20
>= 0.4779 + 0.99966 - 1 - (negligible) >= 0.4775 - 10^-5 >= 0.3995 (e^-2^70, e^-2^20 < 10^-300000, absorbed; verify `assert sum(Decimal(2)**-e for e in (180,220,18,257)) + Decimal(10)**-30 < Decimal('0.00001')`).
over uniform challenges and fresh attack coins. Every returned scalar is
verified by an independent [x]P == Q check, so success means recovery.

## Cost accounting

P-256 field multiplication (schoolbook product + Montgomery REDC, 4
limbs, analytic bound with explicit carries). All walk values live in
the Montgomery domain (aR mod p); partitions, tests, and keys use
Montgomery-form bits, so no domain conversion ever occurs on walk paths.
Schoolbook 512-bit product: 16 limb pairs, each MUL + UMULH + ADDS + ADC
+ ADC with explicit carry propagation (5 instructions): 80; combining
the 16 partials into 8 words (32 ADD/ADC) plus operand loads (8): total
product 120. REDC: for i = 0..3: m = T[i] (p' = 1 since p equiv -1 mod
2^64, free); for j = 0..3: (C,S) = T[i+j] + m*p[j] + C with the 65-bit
carry kept explicitly (MUL + UMULH + ADDS + ADC + ADC + ADC + addressing
= 7): 28 per i; T[i+4] += C with carry to T[i+5] (3; T has 9 words by
REDC invariants) plus loop overhead (2): 33 per i, 132 total. Final
conditional subtract of p (4 SUB + 4 CMP + 4 CSEL = 12). Total at most
120 + 132 + 12 = 264 <= 304 budgeted (F_p256 <= 304, with margin).
Register allocation fits 31 GPRs (product: a(4)+b(4)+T(8)+temps(4) = 20;
REDC: T(9)+m(1)+p(4)+C(1)+temps(4) = 19). Squarings cost the same (S=M
convention, conservative). Fermat inversion I <= 512 P-256-M via
left-to-right square-and-multiply (256 squarings + at most 256
multiplies).

Reference operations H: none. C = 1000 is a borrowed normalization on
this track, not a reference operation; all P-256 arithmetic, including
point additions, is charged as justified instructions in R. H = 0.

Other instruction-equivalents R (C = 1000 borrowed normalization), per
processor-step at most 1969: six P-256 field multiplications x 304
(3M addition formula + 3M invert-batch products per member) = 1824;
six modular 4-limb subtractions x 12 (numerator/denominator differences
and coordinate completions: 4 SUB + 4 CMP + 4 CSEL each) = 72;
partition mask/index 4, table address + point load 10, two 5-limb
integer adds 24 (5 ADCs each with accumulators in registers),
distinguished test 3, point store-back 4, loop control 7 (= 52 walk
operations); round barrier/sync share 8 (lock-step SIMT execution:
synchronization free by construction, 8 retained as margin); batch
control and addressing 12; inverse share 1 (I <= 512 P-256-M divided
across P = 262144 walks: 512/262144 M = 0.0020M x 304 < 1
instruction).
Table insert/lookup at most 2^109 total per the Table construction above
(at most 2^98 operations x 2^11 instructions, deterministically capped).
Fixed work at most 2^73 + 2^38 + 2^46 + 2^102 (periodic
coefficient reduction P*L/2^63 x <=128 instructions, precompute loop
control over P*64*1024 + P*1024 adds at <=8 instructions each (= 2^37+2^31), precompute and verification P-256 additions ((2^34 + 2^28)
batched adds x 1896 + 2^16 batch inversions x 512 x 304, all <= 2^46),
one binary-GCD inverse plus a modular multiplication for the solve (<= 2^21), start-setup sampling (3305209728 <= 2^32: two rounds of P*2*24*16 sampling + 22413312-CAS bitonic sort x 64 + dedup scan, plus 256*24*16 table scalars), decoding/verification bookkeeping and lookup counter (<= 2^21),
zero-initialization of 2^99 bucket heads at <= 4
instructions each, i.e. <= 2^101 inside the 2^102 term (2^38-term total: (2^37+2^31) + 2^32 + 2^21 + 2^21 < 2^38; verify `assert (2**34+2**28)*8 + 3305209728 + 2**21 + 2**21 < 2**38`).
Total R = 1969*S + 2^109 + 2^73 + 2^38 + 2^46 + 2^102 + 2^30 + 2^111
= 769974870319907182878651958018794638540800 (the last two terms are the v15 exceptional-path fallback cap E* x 1024 and the per-round batch zero-check over L_w + TAIL <= 2^111 rounds; verify `assert 1491732361610244130790887340294824 + 2**35 <= 2**111`).

W = R/1000 (H = 0). Claim: log2(W) <= 129.20, i.e. W/2^128 <= 2^1.20.
Proof in the correct direction: 2^1.20 = e^{1.20 ln 2} >= 1+z+z^2/2+z^3/6
with z = 1.20*0.6931 = 0.831720 (using e^w >= sum_{k=0..3} w^k/k! for
w >= 0, all terms positive, and ln 2 > 0.6931 proved below), giving
2^1.20 >= 2.27349059...; and W/2^128 = 2.26275278... for the exact R
above. Since 2.27349059... >= 2.26275278..., the claim holds; verify with
exact Decimal arithmetic:
`R=769974870319907182878651958018794638540800`,
`W=Decimal(R)/1000`,
`z=Decimal('1.20')*Decimal('0.6931')`,
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
Preprocessing (sampling the 128 table-scalar pairs before seeing Q: 256 scalars x 24 attempts x <=16 instructions = 98304 <= 2^17, abort on exhaustion) costs below 2^24 instructions, so
preprocessing_log2 = 24 holds with margin (manifest requires at most
2^24). Starts are sampled in attack phase (charged in R above), so no sampling sort enters preprocessing. All table multipliers (Montgomery-domain double-scalar mults),
inversions, and walk work need Q and run in attack phase as ordinary
charged work, not preprocessing. Nonuniform advice: none (0 bytes); all
multipliers and starts are freshly sampled and their sampling is charged
in R. Single run, no restarts; all trials, failures, sorting (start-setup bitonic sorts only, charged above; the rho table needs none), lookups, transfers, verification, output, and cross-processor
coordination are included above.

## Experiments

`experiment.py` runs the algorithm mechanics at toy scale (v15: distinguished window moved above the partition bits, mirroring production; db-uniformity instrumentation added) (one shared
r=64 (u,v) table with strictly positive entries, P = 8 toy walks from
fresh starts with per-walk coefficient offsets, unreduced integer
coefficient accumulation with mod-n only inside the solve equation,
shared distinguished table, collision solve with degeneracy counting,
independent verification) on y^2=x^3+x+21/F_9973 (n = 9883), 2000 trials:
observed success rate 0.5345 (1069 successes) with ZERO degenerate encounters, every
recovered scalar verified correct. A second toy curve (n = 967) gives
0.81 with 2 degenerate encounters out of 1622 solves; these are
toy-range artifacts (toy coefficient range 2^43 exceeds the toy order 967,
so distinct integers collide mod n). In v9 production residues likewise
wrap (lazy-reduction lemma), so toy and production share the wrapping
regime; production degeneracy is the quantified (<= 1/q + 2^-100) event in (A), and the
toy rate is consistent with wrapping-regime behavior at small scale. The toy also measures coefficient-difference uniformity directly: db over all 1069 production-path solves gives chi-square 19.89 on 24 df (p approx 0.70), consistent with uniform db over Z/n (second curve: 20.29 on 24 df); the experiment asserts chi-square below 40 as a regression guard. This supports the 1/q SHAPE of posit (iv) where measurable; the 2^-100 TV error at production scale remains an explicit posit with margin. (A partition step chi-square was measured during development but dropped: consecutive steps are Markov-dependent, so it confounds dependence with non-uniformity.) Toy tails are minimal (4 steps at density 1/4) yet still
contribute collisions beyond the walk-only prediction (~0.40); this excess
is in the conservative direction and confirms the merge-and-solve
mechanism works, including cross-walk merging onto shared suffixes.
Production tails (2^35 steps) are a 2^-93 fraction of walk points and are
handled by the rigorous tail bound above, not by the toy. The toy therefore
validates mechanism (shared-function merging, offset/monotone
coefficient handling, exact-arithmetic solve, verification) plus the db-uniformity shape; the
production 0.3995 claim rests on the analysis. In particular the toy's
wrapping coefficients (toy range exceeds toy order) and minimal tails
differ from production by construction and are not offered as evidence
for the production birthday constant, the 63/64 correction, or tail
behavior, which are supported by the Teske+Harris regime analysis for the
walk shape used here (production residues wrap exactly as the toy's,
unifying the regimes). The toy's dictionary faithfully implements production matching (exact keys, dedup); production buckets are key bits (no hash function), so occupancy is analysis given R, not a second heuristic. The toy validates walk mechanics only, not P-256
instruction costs, which are bounded analytically in the cost section. Per the rules this is evidence only at toy scope,
not a production proof.

## Limitations

Analytic construction; no full-size execution. Success relies on the
supported multi-walk random-mapping heuristic. Instruction budgets are
conservative analytic upper bounds, not measured assembly; the 2^98-entry
hash table is an analytic device (no physical machine addresses it; the
same holds for the baseline's larger table). Best-known-attack analysis,
not novel cryptanalysis. Degeneracy is a quantified (<= 1/q + 2^-100) heuristic event (corrected from v14's 2/q); coefficient residues wrap per the lazy-reduction lemma; all P-256 arithmetic is instruction-charged (H = 0). Teske's analysis is invoked for its regime with the partition-delta explicitly acknowledged, not as a theorem proving this walk. This supersedes the sequential-walk candidate in
submissions/secp256k1-ecdlp/rho39/, which charged sequential additions at
the batched reference price.
