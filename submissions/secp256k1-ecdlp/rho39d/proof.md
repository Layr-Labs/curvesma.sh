# secp256k1-ecdlp: batched parallel Pollard rho, r=64 shared function, success >= 0.3995

Status: ready for independent review; no qualification or emitted score.
This is a best-known-attack cost analysis (Pollard rho with distinguished
points, P parallel walks with Montgomery-batched affine additions,
analyzed for the track's 0.39 success threshold), not a claim of novelty
or practical execution. This v4 revision uses one shared transition table with r=64 partitions
across all walks (cross-walk collisions merge), 64-bit fresh starts
(duplicate-start probability below 2^-100), sizes the walk with the
63/64 distinct-partition merge factor, specifies the distinguished-point
table (open addressing, load, probe cap, rehash failure analysis), and
charges C=1000 per the protected secp256k1 cost model. It supersedes v3. The normative inputs/outputs are those of the
track profile; `claim.json` binds that profile and
`classical-group-work-v2` by SHA-256 (the "production-v2" label follows
the baseline package's naming). No restricted challenge distribution,
advice, random-oracle premise, or quantum operation is used.

Batching note. The reference operation is batched incomplete affine
addition, so a single sequential walk cannot honestly charge 1 unit per
step (each dependent affine add needs its own inversion). This construction
therefore runs P = 4096 walks in lock-step and Montgomery-batches the P
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
q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141.
r = 64 partitions by low 6 x-coordinate bits (partition 0 for infinity).
One shared multiplier table u_j, v_j in [0, 2^32), sampled once by the
attack (no advice). Each walk has a fresh start (a, b) uniform in
[0, 2^64)^2. Duplicate-start probability over C(4096,2) pairs is below
2^24/2^128 = 2^-104; duplicate starts are discarded by resampling (at
most a handful expected, charged in fixed work), so all walks begin at
distinct states and no identical-path degenerate collisions occur.
Distinguished points: x = 0 mod 2^32 (density 2^-32); infinity never is.

T_walk = 344106074421942661294185815287096969442 (combined walk points),
TAIL = 2^35 steps appended to every walk. Each walk runs
L = T_walk/P + TAIL steps; total steps S = T_walk + P*TAIL. One shared
distinguished-point table (keyed by x AND y-parity, see below); the attack
aborts if the table exceeds 2^98 entries.

## Construction

```text
Input Q. Sample one shared u_j, v_j in [0,2^32); form
M_j = [u_j]P + [v_j]Q (batched below). Shared table = {} (open
addressing per the Table subsection below).
For each of P walks: sample start (a,b) in [0,2^64), resampling on the
negligible duplicate-start event.
For round t = 1..L: each walk computes j = partition(Y), accumulates
(a,b) += (u_j,v_j) as integers (values stay below 2^150 << q; no modular
reduction ever needed), and offers Y + M_j; the P additions are performed
as one Montgomery batch (one inversion per round). On a distinguished Y:
if (x, y-parity) in table with (a0,b0): da = a-a0 mod q, db = b-b0 mod q;
if db != 0: x = -da/db mod q; verify [x]P == Q and return x.
Else store (x, y-parity, a, b); abort if table exceeds 2^98 entries.
```

Coefficient invariant. Per walk at most L < 2^117 steps (L is about
2^116.03) times multipliers below 2^32, plus start below 2^64: a, b <
2^149 + 2^64 < 2^150 always, accumulated as plain integers (never
reduced). Stored pairs are canonical mod-q representatives; 3-limb
(192-bit) accumulators cannot overflow. Since 2^150 << q, two
coefficient vectors are congruent mod q iff they are equal as integers;
this is used in the non-degeneracy proof.
Keying by (x, y-parity) keeps P and -P distinct (same x, opposite parity),
so every table match is a true point equality and the solve equation holds.
Affine edge cases (infinity/doubling/negation inputs, M_j = infinity) occur
with total probability at most S/q < 2^-127 and are handled by the complete
addition routine; their cost is inside the per-step budget and their failure
contribution is absorbed in the success bound.

## Success analysis

Heuristic R (shared-function random mapping with merge factor). One
fixed r-adding walk function with r=64 partitions (shared table) is
iterated from P fresh starts. Two walk points in the same partition
receive the same increment and cannot merge on that step; under the
idealized uniform-partition model a 63/64 fraction of point pairs sit in
distinct partitions and merge as in a random mapping. The combined
T_walk walk-phase points, and the tail sampling for distinguished
purposes, match uniform independent sampling with this merge factor:
(i) the probability of no repeated point among the T_walk walk-phase
points is at most exp(-(63/64)*T_walk(T_walk-1)/2q); (ii) each
post-collision tail step hits a distinguished point with probability
2^-32 independently for the purpose of the tail bound. Support: Pollard's standard heuristic; r-adding walks are
shown close to random mapping and confirmed experimentally on elliptic
curves (Teske; walks here use x-coordinate partitions as in the analyzed
implementations); a single shared iteration function from multiple starts
is the standard parallel-rho-with-distinguished-points setting, so the
heuristic applies more directly than with per-walk functions; the 63/64
refinement answers the equal-increment objection by charging only
distinct-partition pairs; our shared-function toy experiments confirm
merging and recovery end to end (see Experiments). This heuristic is the
only non-rigorous step for the walk phase; everything conditional on it
is proved.

(A) Walk collision. Birthday lemma with merge factor: for T samples,
P(no collision) <= exp(-(63/64)*T(T-1)/2q) via 1-z <= e^-z. The exact
integer T_walk satisfies T_walk(T_walk-1) >= 2q*C0 with
C0 = 0.5113*64/63 = 0.519415873015873..., compensating the 63/64 merge
factor; verify with: `T=346826322875102620124962777434519038877` and the
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
Non-degeneracy (Heuristic R', conservative). A collision with integer
coefficient difference (da,db) = (0,0) does not solve. Starts are fresh
64-bit (duplicates resampled; residual duplicate-start probability below
2^-104), so distinct walks never share an identical path. Integer
coefficients are unreduced 64-bit-start sums of up to 2^117 uniform
32-bit table entries, so each coordinate has spread >= 2^90; two
colliding positions carry independent-start sums, and exact equality of
both 150-bit sums is estimated at ~2^-180 per collision. Heuristic R'
asserts the conservative bound P(degenerate | collision) <= 2^-80 (a
2^100 margin over the estimate), supported by the coefficient spread,
the 64-bit-start duplicate exclusion, and by zero degenerate solves in
4000+ toy recoveries. Expected collisions are T_walk^2/2q < 0.53 (exact: 25*T_walk^2 < 26*q
holds as integers, so T_walk^2/2q < 0.52), so P(any degenerate) <=
0.52 * 2^-80 < 2^-79. Non-degenerate collisions solve: (da,db) != (0,0)
as integers with |da|,|db| < 2^150 < q; if
db = 0 then da != 0 mod q contradicts da + x*db = 0, so db != 0 mod q
and x = -da/db mod q is well defined and correct (verified by
[x]P == Q). Thus P(A), a useful walk-phase collision exists, is
>= 0.40019 - 2^-79 (from P(no collision) <= 0.5998).

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
point no tail would suffice; P <= T_walk*TAIL/q < 10^-25 (a union bound
over at most T_walk cycle-entry candidates each with cycle length below
TAIL with probability at most TAIL/q under the heuristic).

Table abort: distinguished visits are Binomial(S, 2^-32) with mean below
2^97, so P(count > 2^98) <= exp(-2^90) by the Chernoff bound
P(X > 2*mean) <= exp(-mean/3).

Table construction. The shared table is open addressing with linear
probing over 2^99 slots (load at most 1/2 given the 2^98-entry abort
cap). Slot index = splitmix64(x-coordinate) mod 2^99; stored entries
carry the full x, y-parity, and both 150-bit coefficients, so probing
compares exact keys and correctness never depends on the hash. Each
insert/lookup probes at most 64 slots then triggers a rehash (fresh
splitmix salt, table rebuilt); at load <= 1/2 the probability any of
the <= 2^98 operations needs a rehash is below 2^-50 by the standard
linear-probing tail bound, and a rehash costs at most 2^100
instructions (rebuild 2^98 entries). Expected rehash work is below
2^-50 * 2^100 = 2^50, charged as 2^50 in fixed work below. Inserts use
atomic compare-and-swap; contention is negligible at 2^-32 distinguished
density and covered by the per-operation budget. Non-rehash operations
cost at most 200 instructions (64 probes x ~3 instructions: hash-word
compare, key compare, pointer advance), for at most 2^98 * 200 < 2^106
total hash work.

Union bound (no independence assumed between (A) and (B), which share
walk randomness): P(success) >= P(A) + P(B) - 1 - 10^-25 - 2^-79 - e^-2^90
>= 0.40019 + 0.99966 - 1 - (negligible) >= 0.39985 - 10^-20 >= 0.3995
over uniform challenges and fresh attack coins. Every returned scalar is
verified by an independent [x]P == Q check, so success means recovery.

## Cost accounting

Reference operations H: every walk-step addition is one of a Montgomery
batch of P = 4096 (3M formula each, one shared inversion per round), i.e.
a genuine batched incomplete affine addition charging H = 1; precompute
adds are likewise batched across walks (width P); verification adds are
plain reference operations. H = S + P*64*129 + 2^10
= 346826322875102620124962918172041211805, where S = T_walk + P*TAIL,
precompute is 64 multiples x <=129 adds each per walk, and verification
is two 256-bit scalar mults (<= 512 adds each). The reference 6M convention
covers each batch member's 3M formula plus 3M of Montgomery batch products;
all remaining batch overhead is charged in R below, never double-counted.

Other instruction-equivalents R (C = 1000 per the protected secp256k1 cost model), per processor-step at most 85:
partition mask/index 4, precompute table address + point load 10, two
3-limb integer adds 24, distinguished test 3, point store-back 4, loop
control 7 (= 52 walk operations); round barrier/sync share 8; batch excess
25 (inverse share I/P <= 13 with I <= 300M via Fermat exponentiation to p-2 for this 256-bit field under S=M, Montgomery control
and addressing 12). Hash insert/lookup at most 2^106 total per the Table
construction above (at most 2^98 operations x <= 200 instructions;
rehash expectation charged separately).
Fixed work at most 2^28 + 2^50 (precompute batch-excess and loop control
over P*64*129 adds, one binary-GCD inverse plus a modular multiplication
for the solve, all sampling, decoding/verification bookkeeping, start
resampling, expected table-rehash work 2^50).
Total R = 85*S + 2^106 + 2^28 + 2^50
= 29480237525513361125228530866309808929569.

W = H + R/1000. Claim: log2(W) <= 128.146, i.e. W/2^128 <= 2^0.146.
Proof in the correct direction: 2^0.146 = e^{0.146 ln 2} >= 1+z+z^2/2+z^3/6
with z = 0.146*0.6931 = 0.1011926 (using e^w >= sum_{k=0..3} w^k/k! for
w >= 0, all terms positive, and ln 2 > 0.6931 proved below), giving
2^0.146 >= 1.10648527...; and W/2^128 = 1.10586558... for the exact H, R
above. Since 1.10648527... >= 1.10586558..., the claim holds; verify with
exact Decimal arithmetic:
`H=346826322875102620124962918172041211805`,
`R=29480237525513361125228530866309808929569`,
`W=(Decimal(H)*1000+Decimal(R))/1000`,
`z=Decimal('0.146')*Decimal('0.6931')`,
`assert 1+z+z*z/2+z*z*z/6 >= W/Decimal(2**128)`.
And ln 2 > 0.6931: ln 2 = 2*sum_{k>=0} 1/((2k+1)*3^(2k+1)) (positive
terms), whose first four terms give 0.69313475... > 0.6931; verify with:
`2*(Decimal(1)/3+Decimal(1)/81+Decimal(1)/1215+Decimal(1)/15309) >
Decimal('0.6931')`.

Memory: shared table 2^99 open-addressing slots x 96 budgeted bytes
(actual 32-byte x + 1 parity bit + two sub-192-bit scalars = 81 bytes,
load <= 1/2) plus P walk states and one batch buffer (~2^22 bytes):
log2(bytes) <= 105.59.
Preprocessing (all precomputed multiples with batch excess) <= 2^24
group-equivalents, included in H: preprocessing_log2 = 24. Nonuniform
advice: none (0 bytes); all multipliers and starts are freshly sampled
and their sampling is charged in R. Single run, no restarts; all trials,
failures, sorting (none; hash table), lookups, transfers, verification,
output, and cross-processor coordination are included above.

## Experiments

`experiment.py` runs the v4 algorithm mechanics at toy scale (one shared
(u,v) table with r=64 partitions, P = 8 toy walks from fresh 32-bit
starts, unreduced integer coefficient accumulation with mod-n only inside
the solve equation, shared distinguished table, collision solve,
independent verification) on y^2=x^3+x+21/F_9973 (n = 9883), 2000 trials:
observed success rate 0.5135, every recovered scalar verified correct. A
second toy curve (n = 967) gives 0.807. Toy tails are minimal (4 steps at
density 1/4) yet still contribute collisions beyond the walk-only
prediction (~0.40); this excess is in the conservative direction (more
steps find more collisions) and confirms the merge-and-solve mechanism
works, including cross-walk merging onto shared suffixes. Production tails
(2^35 steps) are a 2^-93 fraction of walk points and are handled by the
rigorous tail bound above, not by the toy. The toy therefore validates
mechanism only; the production 0.3995 claim rests on the analysis. Per the
rules this is evidence only at toy scope, not a production proof.

## Limitations

Analytic construction; no full-size execution. Success relies on the
supported multi-walk random-mapping heuristic. Instruction budgets are
conservative analytic upper bounds, not measured assembly; the 2^98-entry
hash table is an analytic device (no physical machine addresses it; the
same holds for the baseline's larger table). Best-known-attack analysis,
not novel cryptanalysis. This supersedes the sequential-walk candidate in
submissions/secp256k1-ecdlp/rho39/, which charged sequential additions at
the batched reference price.
