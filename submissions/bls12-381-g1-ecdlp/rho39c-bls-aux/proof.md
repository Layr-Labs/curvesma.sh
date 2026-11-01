# bls12-381-g1-ecdlp (auxiliary powers): batched parallel Pollard rho with shared walk function, success >= 0.3995

Status: ready for independent review; no qualification or emitted score.
This is a best-known-attack cost analysis (Pollard rho with distinguished
points, P parallel walks with Montgomery-batched affine additions,
analyzed for the track's 0.39 success threshold), not a claim of novelty
or practical execution. The normative inputs/outputs are those of the
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
q = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001.
r = 16 partitions by low x-coordinate bits (partition 0 for infinity).
One shared multiplier table u_j, v_j in [0, 2^32), sampled once by the
attack (no advice). Each walk has a fresh start (a, b) uniform in
[0, 2^32)^2.
Distinguished points: x = 0 mod 2^32 (density 2^-32); infinity never is.

T_walk = 231561926823223774936527941607558918593 (combined walk points),
TAIL = 2^35 steps appended to every walk. Each walk runs
L = T_walk/P + TAIL steps; total steps S = T_walk + P*TAIL. One shared
distinguished-point table (keyed by x AND y-parity, see below); the attack
aborts if the table exceeds 2^98 entries.

## Construction

```text
Input Q. Sample one shared u_j, v_j in [0,2^32); form
M_j = [u_j]P + [v_j]Q (batched below). Shared table = {}.
For each of P walks: sample start (a,b) in [0,2^32).
For round t = 1..L: each walk computes j = partition(Y), accumulates
(a,b) += (u_j,v_j) as integers (values stay below 2^149 << q; no modular
reduction ever needed), and offers Y + M_j; the P additions are performed
as one Montgomery batch (one inversion per round). On a distinguished Y:
if (x, y-parity) in table with (a0,b0): da = a-a0 mod q, db = b-b0 mod q;
if db != 0: x = -da/db mod q; verify [x]P == Q and return x.
Else store (x, y-parity, a, b); abort if table exceeds 2^98 entries.
```

Coefficient invariant. Per walk at most L < 2^117 steps (L is about
2^115.45) times multipliers below 2^32, plus start below 2^32: a, b < 2^148
always, accumulated as plain integers (never reduced). Stored pairs are
canonical mod-q representatives; 3-limb (192-bit)
accumulators cannot overflow. Since 2^148 << q, two coefficient vectors
are congruent mod q iff they are equal as integers; this is used in the
non-degeneracy proof.
Keying by (x, y-parity) keeps P and -P distinct (same x, opposite parity),
so every table match is a true point equality and the solve equation holds.
Affine edge cases (infinity/doubling/negation inputs, M_j = infinity) occur
with total probability at most S/q < 2^-127 and are handled by the complete
addition routine; their cost is inside the per-step budget and their failure
contribution is absorbed in the success bound.

## Success analysis

Heuristic R (shared-function random mapping). One fixed r-adding walk
function (shared table) is iterated from P fresh starts; the combined
T_walk walk-phase points, and the tail sampling for distinguished
purposes, match uniform independent sampling: (i) the probability of no
repeated point among the T_walk walk-phase points is at most
exp(-T_walk(T_walk-1)/2q); (ii) each post-collision tail step hits a
distinguished point with probability 2^-32 independently for the purpose of
the tail bound. Support: Pollard's standard heuristic; r-adding walks are
shown close to random mapping and confirmed experimentally on elliptic
curves (Teske; walks here use x-coordinate partitions as in the analyzed
implementations); a single shared iteration function from multiple starts
is the standard parallel-rho-with-distinguished-points setting, so the
heuristic applies more directly than with per-walk functions; our
shared-function toy experiments confirm merging and recovery end to end
(see Experiments). This heuristic is the only non-rigorous
step; everything conditional on it is proved.

(A) Walk collision. Birthday lemma: for T uniform samples from q elements,
P(no collision) <= exp(-T(T-1)/2q) via 1-z <= e^-z. The exact integer
T_walk satisfies T_walk(T_walk-1) >= 2q(0.5113); verify with:
`T=231561926823223774936527941607558918593` and the full hex q from
Parameters: `assert T*(T-1) >= 2*q*Decimal('0.5113')` (60-digit Decimal).
And ln(1/0.5998) < 0.5113 because e^0.5113 > sum_{k=0..8} 0.5113^k/k!
= 1.66745747... > 1/0.5998 = 1.66722240...; verify with:
`sum(Decimal('0.5113')**k/math.factorial(k) for k in range(9)) >
Decimal(1)/Decimal('0.5998')` (all terms positive, so truncation is a
lower bound). Hence P(no walk collision) <= 0.5998.
Non-degeneracy (Heuristic R', conservative). A collision with integer
coefficient difference (da,db) = (0,0) does not solve. Integer
coefficients are unreduced sums of up to 2^117 uniform 32-bit table
entries, so each coordinate has spread >= 2^90; two colliding positions
carry independent-start sums, and exact equality of both 148-bit sums is
estimated at ~2^-180 per collision. Heuristic R' asserts the
conservative bound P(degenerate | collision) <= 2^-80 (a 2^100 margin
over the estimate), supported by the coefficient spread and by zero
degenerate solves in 4000+ toy recoveries. Expected collisions are
T_walk^2/2q < 0.52 (exact: 25*T_walk^2 < 26*q holds as integers), so
P(any degenerate) <= 0.52 * 2^-80 < 2^-79. Non-degenerate collisions
solve: (da,db) != (0,0) as integers with |da|,|db| < 2^148 < q; if
db = 0 then da != 0 mod q contradicts da + x*db = 0, so db != 0 mod q
and x = -da/db mod q is well defined and correct (verified by
[x]P == Q). Thus P(A), a useful walk-phase collision exists, is
>= 0.40019 - 2^-79 (from P(no collision) <= 0.5998).

(B) Tail collection. After the first collision, all walks that visited
it follow the merged shared-function suffix, which
reaches a distinguished point within TAIL steps with probability
>= 1-(1-2^-32)^TAIL >= 1-e^-8 >= 0.99966. Proof: (1-1/n)^m <= e^(-m/n)
gives (1-2^-32)^(2^35) <= e^-8; e > 1+1+1/2+1/6+1/24+1/120 = 2.71666...
(Taylor lower bound, all terms positive), so e^8 > 2.71666^8 = 2966.81...
and e^-8 < 0.00033707, hence P(B) > 0.9996629 >= 0.99966. Merging is exact (shared
transition function), not heuristic. On a hit the
point is in the shared table (all distinguished visits are stored) and
the stored collision solves unless degenerate (absorbed in (A)).
Short cycles: if the post-collision cycle contained no distinguished
point no tail would suffice; P <= T_walk*TAIL/q < 10^-25 (a union bound
over at most T_walk cycle-entry candidates each with cycle length below
TAIL with probability at most TAIL/q under the heuristic).

Table abort: distinguished visits are Binomial(S, 2^-32) with mean below
2^97, so P(count > 2^98) <= exp(-2^90) by the Chernoff bound
P(X > 2*mean) <= exp(-mean/3).

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
plain reference operations. H = S + P*16*129 + 2^10
= 231561926823223774936528082345055729089, where S = T_walk + P*TAIL,
precompute is 16 multiples x <=129 adds each per walk, and verification
is two 256-bit scalar mults (<= 512 adds each). The reference 6M convention
covers each batch member's 3M formula plus 3M of Montgomery batch products;
all remaining batch overhead is charged in R below, never double-counted.

Other instruction-equivalents R (C = 2000 per the protected BLS cost model), per processor-step at most 85:
partition mask/index 4, precompute table address + point load 10, two
3-limb integer adds 24, distinguished test 3, point store-back 4, loop
control 7 (= 52 walk operations); round barrier/sync share 8; batch excess
25 (inverse share I/P <= 0.15 with I <= 600M via Fermat exponentiation to p-2 for this 381-bit field under S=M (381 squarings plus O(10) multiplies), Montgomery control
and addressing 12). Hash insert/lookup at most 2^106 total (at most 2^98
entries times <= 200 instructions per hash operation, worst case; the
table is shared, inserts use atomic compare-and-swap whose contention is
negligible at 2^-32 density and covered by the 200-instruction budget).
Fixed work at most 2^28 (precompute batch-excess and loop control over
P*16*129 adds, one binary-GCD inverse plus a modular multiplication for
the solve, all sampling, decoding/verification bookkeeping).
Total R = 19682763861103659284211568695118291862805.

W = H + R/2000. Claim: log2(W) <= 127.51, i.e. W/2^127 <= 2^0.51.
Proof in the correct direction: 2^0.51 = e^{0.51 ln 2} >= 1+z+z^2/2+z^3/6
with z = 0.51*0.6931 = 0.353481 (using e^w >= sum_{k=0..3} w^k/k! for
w >= 0, all terms positive, and ln 2 > 0.6931 proved below), giving
2^0.51 >= 1.42331658...; and W/2^127 = 1.41884113... for the exact H, R
above. Since 1.42331658... >= 1.41884113..., the claim holds; verify with
exact Decimal arithmetic:
`H=231561926823223774936528082345055729089`,
`R=19682763861103659284211568695118291862805`,
`W=(Decimal(H)*2000+Decimal(R))/2000`,
`z=Decimal('0.51')*Decimal('0.6931')`,
`assert 1+z+z*z/2+z*z*z/6 >= W/Decimal(2**127)`.
And ln 2 > 0.6931: ln 2 = 2*sum_{k>=0} 1/((2k+1)*3^(2k+1)) (positive
terms), whose first four terms give 0.69313475... > 0.6931; verify with:
`2*(Decimal(1)/3+Decimal(1)/81+Decimal(1)/1215+Decimal(1)/15309) >
Decimal('0.6931')`.

Memory: shared table at most 2^98 entries x 128 budgeted bytes (actual
48-byte BLS x-coordinate + 1 parity bit + two sub-192-bit scalars = 97 bytes) plus
P walk states and one batch buffer (~2^22 bytes): log2(bytes) <= 105.0.
Preprocessing (all precomputed multiples with batch excess) <= 2^24
group-equivalents, included in H: preprocessing_log2 = 24. Nonuniform
advice: none (0 bytes); all multipliers and starts are freshly sampled
and their sampling is charged in R. Single run, no restarts; all trials,
failures, sorting (none; hash table), lookups, transfers, verification,
output, and cross-processor coordination are included above.

## Experiments

`experiment.py` runs the identical shared-function multi-walk attack path
(one shared (u,v) table, P = 8 toy walks from fresh starts, r = 16
x-partition, shared distinguished table, collision solve, independent
verification) on the fixed prime-order toy curve y^2=x^3+x+21/F_9973
(n = 9883), 2000 trials: observed success rate 0.9955 against the 0.4062
walk-phase birthday prediction, every recovered scalar verified correct;
coefficients are unreduced 32-bit-start sums exactly as in production
(reduced mod n only inside the solve equation). Toy parameters
necessarily differ from production (distinguished density 1/4 instead of
2^-32; P = 8 walks instead of 4096); the shared-function structure,
unreduced-coefficient handling, and the 0.5113 walk-sizing constant are
identical. A second toy curve (n = 967) gives 0.998 against 0.4013. The
excess over the walk-only prediction is merged-suffix tail collection
working as designed (cross-walk collisions visibly merge onto shared
suffixes in the toy logs). This confirms the merge mechanism and solve
logic end-to-end; per the rules it is evidence only at toy scope, not a
production proof.

## Limitations

Analytic construction; no full-size execution. Success relies on the
supported multi-walk random-mapping heuristic. Instruction budgets are
conservative analytic upper bounds, not measured assembly; the 2^98-entry
hash table is an analytic device (no physical machine addresses it; the
same holds for the baseline's larger table). Best-known-attack analysis,
not novel cryptanalysis. This ports the batch-corrected parallel-rho analysis to the BLS12-381 scalar field; walk/table/solve structure is identical.

## Auxiliary-powers input (this track)

The track supplies public powers A_i = [s^i]P0 for i = 1..n. The attack
reads A_1 = [s]P0 (one point decode, charged inside the 2^28 fixed-work
budget) and runs the identical parallel rho above with Q = A_1, ignoring
A_2..A_n. Since A_1 = [s]P0 with s uniform (the track's secret), recovery
of the discrete log of A_1 is exactly recovery of s; success probability
and cost are as analyzed. Total work W(n) = W is independent of n (one
extra point read regardless of n), so log2(W(2^20)) = 127.51. No Cheon
step is used; this is a conservative generic-group upper bound.
