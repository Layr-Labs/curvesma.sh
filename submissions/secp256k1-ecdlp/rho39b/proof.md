# secp256k1-ecdlp: batched parallel Pollard rho, success >= 0.3995

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

## Parameters

Let q be the protected subgroup order
q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141.
r = 16 partitions by low x-coordinate bits (partition 0 for infinity).
Each walk has fresh uniform multipliers u_j, v_j in [0, 2^32) and a fresh
start (a, b) uniform in [0, 2^32)^2, sampled by the attack (no advice).
Distinguished points: x = 0 mod 2^32 (density 2^-32); infinity never is.

T_walk = 344106074421942661294185815287096969442 (combined walk points),
TAIL = 2^35 steps appended to every walk. Each walk runs
L = T_walk/P + TAIL steps; total steps S = T_walk + P*TAIL. One shared
distinguished-point table (keyed by x AND y-parity, see below); the attack
aborts if the table exceeds 2^98 entries.

## Construction

```text
Input Q. For each of P walks: sample u_j, v_j, start (a,b) in [0,2^32);
form M_j = [u_j]P + [v_j]Q (batched below). Shared table = {}.
For round t = 1..L: each walk computes j = partition(Y), accumulates
(a,b) += (u_j,v_j) as integers (values stay below 2^149 << q; no modular
reduction ever needed), and offers Y + M_j; the P additions are performed
as one Montgomery batch (one inversion per round). On a distinguished Y:
if (x, y-parity) in table with (a0,b0): da = a-a0 mod q, db = b-b0 mod q;
if db != 0: x = -da/db mod q; verify [x]P == Q and return x.
Else store (x, y-parity, a, b); abort if table exceeds 2^98 entries.
```

Coefficient invariant. Per walk at most L < 2^117 steps (L is about
2^116.02) times multipliers below 2^32, plus start below 2^32: a, b < 2^149
always. Stored pairs are canonical mod-q representatives; 3-limb (192-bit)
accumulators cannot overflow.
Keying by (x, y-parity) keeps P and -P distinct (same x, opposite parity),
so every table match is a true point equality and the solve equation holds.
Affine edge cases (infinity/doubling/negation inputs, M_j = infinity) occur
with total probability at most S/q < 2^-127 and are handled by the complete
addition routine; their cost is inside the per-step budget and their failure
contribution is absorbed in the success bound.

## Success analysis

Heuristic R (multi-walk random mapping). The combined T_walk walk-phase
points across the P independent walks, and the tail sampling for
distinguished purposes, match uniform independent sampling: (i) the
probability of no repeated point among the T_walk walk-phase points is at
most exp(-T_walk(T_walk-1)/2q); (ii) each post-collision tail step hits a
distinguished point with probability 2^-32 independently for the purpose of
the tail bound. Support: Pollard's standard heuristic; r-adding walks are
shown close to random mapping and confirmed experimentally on elliptic
curves (Teske; walks here use x-coordinate partitions as in the analyzed
implementations); independent walks with fresh randomness are the standard
parallel-rho setting; our multi-walk toy experiments confirm the predicted
success rates (see Experiments). This heuristic is the only non-rigorous
step; everything conditional on it is proved.

(A) Walk collision. Birthday lemma: for T uniform samples from q elements,
P(no collision) <= exp(-T(T-1)/2q) via 1-z <= e^-z. The exact integer
T_walk satisfies T_walk(T_walk-1) >= 2q(0.5113); verify with:
`T=344106074421942661294185815287096969442` and the full hex q from
Parameters: `assert T*(T-1) >= 2*q*Decimal('0.5113')` (60-digit Decimal).
And ln(1/0.5998) < 0.5113 because e^0.5113 > sum_{k=0..8} 0.5113^k/k!
= 1.66745747... > 1/0.5998 = 1.66722240...; verify with:
`sum(Decimal('0.5113')**k/math.factorial(k) for k in range(9)) >
Decimal(1)/Decimal('0.5998')` (all terms positive, so truncation is a
lower bound). Hence P(no walk collision) <= 0.5998.
Degenerate collisions: conditional on a collision, (da,db) lies on the
line da + x*db = 0, so it is (0,0) with heuristic probability 1/q;
expected collisions <= T_walk^2/2q < 0.52 (exact: 25*T_walk^2 < 26*q
holds as integers), so P(any degenerate) <= 0.52/q < 10^-75.
Thus P(A), a useful walk-phase collision exists, is >= 0.40019.

(B) Tail collection. After the first collision the merged walk suffix
reaches a distinguished point within TAIL steps with probability
>= 1-(1-2^-32)^TAIL >= 1-e^-8 >= 0.99966. Proof: (1-1/n)^m <= e^(-m/n)
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

Union bound (no independence assumed between (A) and (B), which share
walk randomness): P(success) >= P(A) + P(B) - 1 - 10^-25 - 10^-75 - e^-2^90
>= 0.40019 + 0.99966 - 1 - (negligible) >= 0.39985 - 10^-20 >= 0.3995
over uniform challenges and fresh attack coins. Every returned scalar is
verified by an independent [x]P == Q check, so success means recovery.

## Cost accounting

Reference operations H: every walk-step addition is one of a Montgomery
batch of P = 4096 (3M formula each, one shared inversion per round), i.e.
a genuine batched incomplete affine addition charging H = 1; precompute
adds are likewise batched across walks (width P); verification adds are
plain reference operations. H = S + P*16*129 + 2^10
= 344106074421942661294185956024593779938, where S = T_walk + P*TAIL,
precompute is 16 multiples x <=129 adds each per walk, and verification
is two 256-bit scalar mults (<= 512 adds each). The reference 6M convention
covers each batch member's 3M formula plus 3M of Montgomery batch products;
all remaining batch overhead is charged in R below, never double-counted.

Other instruction-equivalents R (C = 1000), per processor-step at most 85:
partition mask/index 4, precompute table address + point load 10, two
3-limb integer adds 24, distinguished test 3, point store-back 4, loop
control 7 (= 52 walk operations); round barrier/sync share 8; batch excess
25 (inverse share I/P <= 13 with I <= 300M via Fermat, Montgomery control
and addressing 12). Hash insert/lookup at most 2^106 total (at most 2^98
entries times <= 200 instructions per hash operation, worst case; the
table is shared, inserts use atomic compare-and-swap whose contention is
negligible at 2^-32 density and covered by the 200-instruction budget).
Fixed work at most 2^28 (precompute batch-excess and loop control over
P*16*129 adds, one binary-GCD inverse plus a modular multiplication for
the solve, all sampling, decoding/verification bookkeeping).
Total R = 29249016406994764624612487957879026184970.

W = H + R/1000. Claim: log2(W) <= 128.135, i.e. W/2^128 <= 2^0.135.
Proof in the correct direction: 2^0.135 = e^{0.135 ln 2} >= 1+z+z^2/2+z^3/6
with z = 0.135*0.6931 = 0.0935685 (using e^w >= sum_{k=0..3} w^k/k! for
w >= 0, all terms positive, and ln 2 > 0.6931 proved below), giving
2^0.135 >= 1.09808256...; and W/2^128 = 1.09719200... for the exact H, R
above. Since 1.09808256... >= 1.09719200..., the claim holds; verify with
exact Decimal arithmetic:
`H=344106074421942661294185956024593779938`,
`R=29249016406994764624612487957879026184970`,
`W=(Decimal(H)*1000+Decimal(R))/1000`,
`z=Decimal('0.135')*Decimal('0.6931')`,
`assert 1+z+z*z/2+z*z*z/6 >= W/Decimal(2**128)`.
And ln 2 > 0.6931: ln 2 = 2*sum_{k>=0} 1/((2k+1)*3^(2k+1)) (positive
terms), whose first four terms give 0.69313475... > 0.6931; verify with:
`2*(Decimal(1)/3+Decimal(1)/81+Decimal(1)/1215+Decimal(1)/15309) >
Decimal('0.6931')`.

Memory: shared table at most 2^98 entries x 96 budgeted bytes (actual
32-byte x + 1 parity bit + two sub-192-bit scalars = 81 bytes) plus
P walk states and one batch buffer (~2^22 bytes): log2(bytes) <= 104.59.
Preprocessing (all precomputed multiples with batch excess) <= 2^24
group-equivalents, included in H: preprocessing_log2 = 24. Nonuniform
advice: none (0 bytes); all multipliers and starts are freshly sampled
and their sampling is charged in R. Single run, no restarts; all trials,
failures, sorting (none; hash table), lookups, transfers, verification,
output, and cross-processor coordination are included above.

## Experiments

`experiment.py` runs the identical multi-walk attack path (P = 8 toy walks,
r = 16 x-partition, shared distinguished table, collision solve,
independent verification) on the fixed prime-order toy curve
y^2=x^3+x+21/F_9973 (n = 9883), 2000 trials: observed success rate 0.858
against the 0.4062 walk-phase birthday prediction, every recovered scalar
verified correct. Toy parameters necessarily differ from production
(u,v sampled in [0,n) with mod reduction instead of [0,2^32) without;
distinguished density 1/4 instead of 2^-32; P = 8 walks instead of 4096);
the walk/table/solve structure and the 0.5113 walk-sizing constant are
identical. A second toy curve (n = 967) run locally with the same code
gives 1.0 against 0.4013. The excess over the walk-only prediction is
tail-phase collisions the bound conservatively ignores. This confirms the
algorithm logic end-to-end and the heuristic's order-of-magnitude fit at
small scale; per the rules it is evidence only at toy scope, not a
production proof.

## Limitations

Analytic construction; no full-size execution. Success relies on the
supported multi-walk random-mapping heuristic. Instruction budgets are
conservative analytic upper bounds, not measured assembly; the 2^98-entry
hash table is an analytic device (no physical machine addresses it; the
same holds for the baseline's larger table). Best-known-attack analysis,
not novel cryptanalysis. This supersedes the sequential-walk candidate in
submissions/secp256k1-ecdlp/rho39/, which charged sequential additions at
the batched reference price.
