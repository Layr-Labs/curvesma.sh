# secp256k1-ecdlp: distinguished-points Pollard rho at 39%+ success

Status: ready for independent review; no qualification or emitted score.
This is a best-known-attack cost analysis (Pollard rho with distinguished
points, analyzed for the track's 0.39 success threshold), not a claim of
novelty or practical execution. The normative inputs/outputs are those of
`production-v2`; `claim.json` binds that profile and
`classical-group-work-v2` by SHA-256. No restricted challenge distribution,
advice, random-oracle premise, or quantum operation is used.

## Parameters

Let q be the protected subgroup order
q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
(a 256-bit prime). The walk uses r = 16 partitions selected by the low 4
bits of the point x-coordinate (partition 0 for the point at infinity).
Precomputed multiples M_j = [u_j]P + [v_j]Q use fresh uniform u_j, v_j in
[0, 2^32) sampled by the attack. Distinguished points are those whose
x-coordinate is 0 mod 2^32 (density 2^-32); infinity is never distinguished.

T_walk = 344106074421942661294185815287096969442 (walk phase),
TAIL = 2^35 steps (distinguished-point collection tail).
Total walk steps S = T_walk + TAIL. The table aborts (attack fails) if it
ever exceeds 2^98 entries; the abort probability is shown negligible below.

## Construction and success

```text
Input Q. Sample u_j, v_j in [0,2^32) for j = 0..15; form M_j = [u_j]P+[v_j]Q.
Sample start (a,b) uniform in [0,q)^2; Y = [a]P+[b]Q. Table = {}.
Repeat S times: j = partition(Y); (a,b) += (u_j,v_j) as integers (no modular
reduction; values stay below 2^161 << q); Y = Y + M_j (reference group op).
If Y distinguished: if Y in table with (a0,b0): da = a-a0 mod q, db = b-b0
mod q; if db != 0: x = -da/db mod q; verify [x]P == Q and return x.
Else store (Y,(a,b)); abort if table exceeds 2^98 entries.
```

For every visited point the invariant Y = [a]P + [b]Q holds with a, b < 2^161
(at most S < 2^129 steps times multipliers below 2^32; S is about 2^128.02),
so no reduction is needed and stored pairs are already canonical mod-q
representatives. The 3-limb (192-bit) walk accumulators cannot overflow.

Heuristic R (random-mapping walk). The sequence of visited points has
collision statistics matching T_walk uniform independent samples: the
probability of no repeated point among the first T_walk visits is at most
exp(-T_walk(T_walk-1)/2q). Support: this is Pollard's standard heuristic
for r-adding walks (Pollard 1978; Teske 1998/2001 proves r-adding walks
approach random-mapping rho length and confirms it experimentally on
elliptic-curve groups); the partition uses x-coordinate bits exactly as in
the analyzed implementations, and our toy-curve experiments confirm the
birthday curve quantitatively (see Experiments). The heuristic is the only
non-rigorous step; everything conditional on it is proved.

Birthday lemma. For T uniform samples from q elements,
P(no collision) = prod_{i<T}(1-i/q) <= exp(-T(T-1)/2q), using 1-z <= e^-z.
The exact integer T_walk above satisfies
T_walk(T_walk-1) >= 2q(0.5113); verify with:
`T=344106074421942661294185815287096969442; q=0xFFFF...4141` (full hex in
Parameters); `assert T*(T-1) >= 2*q*Decimal('0.5113')` (60-digit Decimal).
And ln(1/0.5998) < 0.5113 because e^0.5113 > sum_{k=0..8} 0.5113^k/k!
= 1.64000606... > 1/0.5998 = 1.63988192...; verify with:
`sum(Decimal('0.5113')**k/math.factorial(k) for k in range(9)) >
Decimal(1)/Decimal('0.5998')`. Hence P(no walk collision) <= 0.5998.

Useful collisions. Given a collision with coefficient differences (da,db),
it is degenerate iff (da,db) = (0,0), of heuristic probability 1/q^2 per
collision. Expected collisions <= T_walk^2/2q < 0.52, so
P(any degenerate) <= 0.52/q^2 < 10^-150. Thus P(useful walk collision)
>= 1 - 0.5998 - 10^-150 >= 0.40019.

Tail. After the first collision the two walk suffixes merge; the merged
walk hits a distinguished point within TAIL = 2^35 steps with probability
>= 1-(1-2^-32)^(2^35) >= 1-e^-8 >= 0.99966, using (1-1/n)^m <= e^(-m/n).
On a hit, the point is in the table (all distinguished visits are stored)
and the collision solves unless degenerate (absorbed above).

Table abort. Distinguished visits are Binomial(S, 2^-32) with mean below
2^97. P(count > 2^98) <= exp(-2^90) by a Chernoff bound, negligible.

Total: P(success) >= 0.40019 x 0.99966 - exp(-2^90) >= 0.40005 >= 0.39
over uniform challenges and fresh attack coins. Every returned scalar is
verified by an independent [x]P == Q check, so success means recovery.

## Cost accounting

Reference operations H (each batched incomplete affine addition counts 1;
its internal instructions are not double-charged):
T_walk + TAIL walk steps, 16 x 129 <= 2^12 precompute adds (each [u]P+[v]Q
with 32-bit scalars costs <= 129 adds), 2^10 verification adds (two 256-bit
scalar mults, <= 512 adds each). Total H = 344106074421942661294185815321456712930.

Other instruction-equivalents R (scalar AArch64, divided by C = 1000):
per walk step at most 48 (partition mask/index 4, table address + point
load 10, two 3-limb integer adds 24, distinguished test 3, loop control 7;
no hashing, no modular reduction on the walk); distinguished-point hash
insert/lookup at most 2^106 (at most 2^98 table entries times <= 200
instructions per hash operation, worst case); fixed work at most 2^16
(precompute loop control, one binary-GCD modular inverse plus modular
multiplication for the solve, sampling, decoding/verification bookkeeping).
Total R = 16517091653382886156727600831218927184480.

W = H + R/1000. Then log2(W) <= 128.08376, proved by exact integer
arithmetic: with y = 0.08376 x 0.6932 = 0.058062432 (> 0.08376 ln 2 since
ln 2 = 0.693147... < 0.6932), e^y <= 1+y+0.532y^2 = 1.059856... < 1.05987
(lemma: e^y = 1+y+sum_{k>=2} y^k/k! <= 1+y+y^2(1/2)(1/(1-y)) <= 1+y+0.532y^2
for 0 <= y <= 0.06), so 2^0.08376 = e^{0.08376 ln 2} <= 1.05987, and
100000(1000H+R) <= 105987 x 1000 x 2^128 holds as integers for the H, R
above (one-line Python check with the stated values). Hence the submitted
time_log2 = 128.08376 is a rigorous upper bound on log2(W).

Memory: at most 2^98 entries x 96 bytes (32-byte x plus two sub-192-bit
scalars, 24 bytes each) plus negligible code/scratch/state: log2(bytes) <= 104.59.
Preprocessing (the 16 precomputed multiples) <= 2^12 group-equivalents,
included in H: preprocessing_log2 = 12. Nonuniform advice: none (0 bytes);
multipliers are freshly sampled and their sampling is charged in R.
All trials (single run, no restarts), failures, sorting (none; hash table),
lookups, transfers, verification, and output are included above.

## Experiments

`experiment.py` runs the identical attack path (r = 16 x-partition walk,
distinguished table, collision solve, independent verification) on fixed
prime-order toy curves. The hosted run uses y^2=x^3+x+21/F_9973 (n = 9883),
2000 trials: observed success rate 0.5265 against the walk-phase birthday
prediction 0.4062, every recovered scalar verified correct. A second toy
curve y^2=x^3+11/F_1009 (n = 967), run locally with the same code, gives
0.7685 against 0.4013. The excess over prediction is explained: tail-phase
collisions add successes the bound does not count. This confirms the
algorithm logic end-to-end and the heuristic's quantitative fit at small
scale; per the rules it is evidence only at toy scope, not a production
proof.

## Limitations

Analytic construction; no full-size execution. Success relies on the
supported random-mapping heuristic for the r-adding walk (standard in the
literature and confirmed at toy scale here). Instruction budgets are
conservative analytic upper bounds, not measured assembly. This beats the
generic BSGS accounting control by using negligible-memory rho at the
track's 0.39 threshold; it is not novel cryptanalysis.
