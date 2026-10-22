# secp256k1-ecdlp: orbit-collision r-adding walk (v2, 1.05 sqrt(N) budget)

Status: ready for independent exploratory review. This is an analytic classical
attack bound with bounded supporting experiments; it is not a production run, a
practical break, or a claim of novelty. The protected target profile,
`classical-group-work-v2`, and the current `C = 1000` secp256k1 reference price are
used unchanged. `claim.json` binds the protected files by SHA-256.

v2 is a direct revision of the published v1 submission
(`submissions/secp256k1-ecdlp/orbit-rho-v1`, submission UUID
`4cd8a282-d747-4e3e-b67e-8bfdedba5ccd`). The algorithm, correctness argument and
orbit-count proof are unchanged. The changes are (a) the success budget is
`1.05 sqrt(N)` instead of `1.30 sqrt(N)`, supported by a 2,000-trial
success-fraction measurement, and (b) the per-step increment/table charge is
itemized as 200 instructions instead of the v1 conservative 320. The resulting
bound is `127.644312` instead of `128.043099`.

## Target and inputs

Let `p` be the protected field modulus, `r` the prime subgroup order, and `P` the
protected standard generator of secp256k1. The group `E(F_p)` has order `r`
(cofactor 1), so `<P> = E(F_p)`. The input is a single point `Q = [x]P` with `x`
uniform in `[0, r)`. The output is the canonical 32-byte scalar `x`. Acceptance is
exactly `[x]P = Q`, with `0 <= x < r`; `P` has prime order, so a verified answer
is unique.

The attack receives only public inputs: `p`, `r`, `P`, and `Q`. It does not use the
challenge sampling coins, hidden advice, an oracle, or a quantum step. The
challenge is cold: nothing is amortized over future targets.

## Verified six-element automorphism

The public secp256k1 constants are

```text
beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
lam  = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
phi(x, y) = (beta * x mod p, y)
```

For the protected curve and generator, direct integer checks give
`beta^3 = 1 mod p`, `beta != 1`, `lam^3 = 1 mod r`, `lam != 1`, and
`phi(P) = [lam]P`. Hence `phi` is an efficiently computable automorphism of
`E(F_p)` of order 3, and `phi([a]P) = [lam a]P`. It costs one base-field
multiplication; `phi^2(x, y) = (beta^2 x, y)` costs a second one.

Let

```text
H = {1, -1, lam, -lam, lam^2, -lam^2}  (mod r).
```

The six scalars are distinct because `r` is odd and `lam` has exact order 3. The
identity point is fixed by every element of `H`. For every nonzero point `X`, the
orbit `H X` has exactly six elements: if `uX = X` for `u in H`, then
`(u - 1)X = 0`; since `X` has prime order `r`, this forces `u = 1 mod r`, and the
six scalars are distinct modulo `r`.

The remaining possible fixed points would have `x = 0`, because `beta x = x` and
`beta != 1` imply `x = 0`. Euler's criterion gives
`7^((p-1)/2) = -1 mod p`, so `y^2 = 7` has no solution and no point has `x = 0`.
Thus the only non-generic orbit is the identity, and the number of orbits is
exactly

```text
N = (r - 1)/6 + 1
  = 19298681539552699237261830834781317975472927379845817397100860523586360249057.
```

## Algorithm

Parameters: `K = 2^20` increments, `c = 1.05`, and
`T = ceil(c * sqrt(N)) = 145865679298993601362043284063376113664`
(`log2 T = 126.77790807753082`). All scalars below are modulo `r`.

1. **Increments.** Choose `alpha_k, beta_k` uniformly in `[0, r)` for
   `k = 0..K-1`. Compute `M_k = [alpha_k]P + [beta_k]Q` and store
   `(M_k, alpha_k, beta_k)`.
2. **Walk state.** Choose `a_0` uniformly in `[0, r)` and set `b_0 = 0`,
   `X_0 = [a_0]P`, so `X_0 = [a_0]P + [b_0]Q`.
3. **Canonical class.** For the current `X`, compute the six points
   `X, -X, phi(X), -phi(X), phi^2(X), -phi^2(X)` and select the
   lexicographically least `(infinity_tag, x, y)` as `C = canon(X)`. Record the
   unique `u in H` with `X = u C`.
4. **Collision table.** Probe a hash table keyed by `C`. If `C` is absent,
   insert `(C, u, a, b)`; if present as `(C, u', a', b')`, go to step 6.
5. **Step.** Let `k = H64(C) mod K` for a fixed 64-bit mixer `H64`. Set
   `X <- X + M_k`, `a <- a + alpha_k`, `b <- b + beta_k`, and repeat from step 3.
   Stop after `T` steps and report failure.
6. **Solve.** The current point is `X = u C` and the stored point is `X' = u' C`,
   so `X = w X'` with `w = u (u')^{-1} in H`. Because
   `X = [a]P + [b]Q = [a + b x]P`, the collision gives
   `a + b x = w (a' + b' x) mod r`. If `d = b - w b' != 0 mod r`, compute
   `x = (w a' - a) d^{-1} mod r`, verify `[x]P = Q`, and output `x`. If `Q` is
   the identity, output `0` directly. If `d = 0`, report failure.

The table has `2^128` slots of 128 bytes (canonical class, orbit index, two
256-bit scalars), with a one-bit-per-slot occupancy bitmap. The record area is
written only on insertion; the bitmap is the only array that must be cleared.
The walk is a single logical sequence; with parallel walkers, each walker keeps
its 96-byte state in registers/local scratch, and the state accesses are covered
by the arithmetic and control charges below.

## Correctness

Every valid `Q` has a unique `x in [0, r)`. If the walk visits two points in the
same orbit, step 6 recovers a scalar satisfying `[x]P = Q` whenever the
denominator is nonzero. The attack-side check `[x]P = Q` accepts only the correct
scalar. A wrong or trivial relation cannot produce a wrong output; it produces
failure. The only failure modes are budget exhaustion and a zero denominator;
both are counted as failure, not as false success.

## Success probability

The class sequence is an r-adding walk whose step index is a deterministic hash
of the canonical orbit representative. Under the standard random-walk heuristic,
the visited classes behave like independent uniform draws from the `N` classes
until the first repeat. The first-repeat time has

```text
E[T_rho] = sqrt(pi N / 2) = 1.253314... * sqrt(N),
P(no repeat after t) ~ exp(-t^2 / (2N)).
```

With `T = 1.05 sqrt(N)`, the ideal success probability is
`1 - exp(-1.05^2 / 2) = 0.4237709`. The claim conservatively states at least
`0.39`, which leaves an 8.2% margin.

Supporting evidence, measured on the same full-storage orbit walk at
`N = 166668` classes (`order = 1000003`, `K = 1024`, 2,000 trials, seed
20261022), is a direct success-fraction curve rather than only a mean:

| budget `c` | observed `P(T <= c sqrt(N))` | ideal `1-exp(-c^2/2)` |
| ---: | ---: | ---: |
| 1.00 | 0.3945 | 0.3935 |
| 1.02 | 0.4040 | 0.4056 |
| 1.05 | 0.4245 | 0.4238 |
| 1.10 | 0.4565 | 0.4539 |

The observed mean collision time was `1.2474 sqrt(N)` against the ideal
`1.2533 sqrt(N)`. The full-storage design detects any repeat, including a short
cycle, so it does not suffer the distinguished-point fruitless-cycle loss; this
is why its curve tracks the ideal. The 2,000-trial standard error at `p = 0.42`
is about 0.011, and the observed `0.4245` is well above `0.39`. Even applying
the 0.924 toy efficiency measured by the low-memory distinguished-point variant
would give `0.4238 * 0.924 = 0.3916`, still above `0.39`.

Fresh algorithm coins are the `K` increment pairs and the starting point `a_0`.
The challenge scalar `x` is uniform. A zero denominator occurs with probability
`O(1/r)` per collision and is included in the failure probability; no restart or
amplification is claimed, so the work bound is for one run.

## Complete work ledger

The protected normalization is `C = 1000` scalar AArch64 instruction
equivalents per secp256k1 reference operation, the batched incomplete affine
addition (6M, S = M, A = 0). The walk may be run as `M` independent walks on
`M` processors, all sharing the increment table, the collision table and one
batch inversion per round; the total work is the sum over processors, and the
first collision between any two walks is used. The single-walk description is
only for the collision count; parallelization does not change the total work
bound. One walk step performs exactly one such reference addition. Exceptional
additions (doubling, opposite, infinity) occur with probability `O(1/r)` for
random points and are covered by the control reserve; their expected work is
negligible.

Per-step non-reference work, as an upper bound in scalar AArch64 instructions:

| Component | Instructions |
| --- | ---: |
| `phi(X)` and `phi^2(X)`: two pinned 5x52 field multiplications at 179 each | 358 |
| canonical selection: negations, six comparisons, tie-break, index | 96 |
| step-index mixer `H64(C)` | 32 |
| two 4-limb modular additions for `a` and `b` | 64 |
| increment fetch and collision-table probe/compare/write/bitmap | 200 |
| loop, branch, address and state control | 64 |
| **non-reference total R** | **814** |
| reference addition | 1000 |
| **per-step total** | **1814** |

The 200-instruction increment/table charge is itemized as follows. One step
touches 16 words for the 128-byte increment record (point plus two scalars),
4+1 words for the probed 32-byte key plus tag, and 16 words for the 128-byte
collision record write: 37 words. At 4 instructions per 64-bit word (load or
store plus address and loop work) this is 148 instructions; the remaining 52
cover the table hash, key comparison, bitmap test/set, and branch. No memory
latency is charged (the cost model counts instructions), and the single-walk
state is register-resident. This is an upper bound, not a measured trace.

Other charged work:

| Component | Work (target-group equivalents) |
| --- | ---: |
| increment construction: `2^20` increments at 1024 reference ops each | 1,073,741,824 |
| increment randomness: `2^20` increments at 512 instructions each | 536,870,912 |
| occupancy bitmap: `2^125` bytes = `2^122` words, clear at about 2 instructions/word | 2^120 |
| collision solve, final scalar multiplication, verification, output | 1024 |
| **preprocessing total** | **1,610,612,736 = 2^30.584962500721158** |

The total instruction count is

```text
T * 1814 + (1,610,612,736 + 2^120 + 1024) * 1000
= 265929570244159308743650324352855228522496.
```

Dividing by `C = 1000` gives `2^127.64431181342555` target-group-operation
equivalents. The submitted score is the rounded-up six-decimal value
`127.644312`, i.e. the exact rational `15955539 / 125000`. This includes the
complete walk budget (successful or failing), all increment construction and
randomness, every table insertion and lookup, memory transfers, verification,
and output. Nothing is amortized over future targets.

## Memory

```text
collision table: 2^128 slots * 128 bytes = 2^135 bytes
occupancy bitmap: 2^125 bytes
increment table: 2^20 * 128 bytes = 2^27 bytes
code, constants, stacks and scratch: 2^24 bytes
peak = 43598678261745240631244871577195453087744 bytes
     = 2^135.0014081943928 bytes.
```

This is below the organizer baseline's `2^136.000001` bytes. It is astronomical
and mandatory for this construction; the attack is an accounting bound, not a
practical computation. Nonuniform advice is 0 bytes: the increment table is
built from the public input `Q` at run time and is fully charged.

## Heuristics, objections and limitations

1. **Random-walk heuristic.** The deterministic r-adding class sequence is
   assumed to behave like independent uniform draws until its first repeat. This
   is the standard Pollard-rho heuristic. The 2,000-trial success-fraction curve
   above is direct evidence at 166,668 classes; it is not a production proof.
2. **Random-hash heuristic.** `H64` is a fixed deterministic 64-bit mixer, not a
   cryptographic hash; it is assumed to randomize the step index. The same mixer
   is used in the toy check.
3. **Orbit freeness is not a heuristic.** The six scalars are distinct, no
   nonzero point is fixed, and `x = 0` is impossible; the orbit count is exact.
4. **Fruitless cycles.** Because the attack stores every visited orbit class and
   stops at the first repeat, any cycle, including a 2-cycle, is itself a
   collision and yields a relation. The cost is paid in `2^135` bytes of memory
   and in the 200-instruction table charge, not in lost steps. The published
   1.29x negation-map loss for distinguished-point walks does not directly
   apply; it is a low-memory trade this attack does not make.
5. **No production execution.** The bound is analytic. The profile check in
   `experiment.py` verifies the production endomorphism constants; the toy
   check validates the class-walk algebra and success curve at 166,668 classes.
   Neither is a production run or a proof of the production constant.
6. **Normalization.** The 179-instruction field multiplication and `C = 1000` are
   the organizer's pinned normalization constants, not hardware timings. A
   different field representation would change the instruction count and must be
   re-justified.

Sensitivity, disclosed rather than hidden:

- If the increment/table charge were the v1 conservative 320 instructions
  instead of 200, the same budget gives `2^127.736278`.
- If the budget were `1.10 sqrt(N)` instead of `1.05 sqrt(N)` at 200, the bound
  is `2^127.711099`.
- If the budget were `1.02 sqrt(N)`, the bound is `2^127.602704`.

The submitted point is `1.05 sqrt(N)` with the itemized 200 charge. The
2,000-trial curve supports the budget; the itemized transfer count supports the
charge.

## Evidence and provenance

- `experiment.py` — bounded toy orbit-walk and production profile check; reads
  JSON on stdin and writes JSON. It reports the success fraction at the declared
  budget, the mean collision ratio, the exact ledger integers and inequality
  checks.
- `accounting.json` — exact integers, hashes and sensitivity figures.
- `claim.json` — machine-readable claim bound to the current protected files.
- `README.md` — package summary.

Base checkout: `Layr-Labs/curvesma.sh` at
`467130d867499e7bbca7479f0f72af0df115b282`. Protected hashes: target profile
`18c3338de676f1909715814325bffdb1513f1621d678b942e8a7a9563b8fd086`, cost
model `e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294`.

Prior work: Pollard (1978); van Oorschot and Wiener, parallel collision search
(1999); Wiener and Zuccherato, negation map (1998); Duursma, Gaudry and Morain,
automorphisms (1999); Teske, r-adding walks (2001); Bos, Kleinjung and Lenstra,
negation-map cycles (2010); Gallant, Lambert and Vanstone, endomorphisms (2001);
SEC 2 for the curve parameters. The construction and every charged quantity are
stated above, so review does not need to fetch a citation to reconstruct the
claim.

v2 is a revision of v1: the algorithm, correctness argument and orbit count are
reused, and the changes are the budget, the itemized table charge and the new
success curve. v1 is currently unpublished because the track runner rejected its
dispatch; per the dependency rules, v2 submits without a dependency link and the
link can be added through the correction API after both submissions are
published. No other prior submission's code or distinctive technique is
incorporated. The organizer baseline and public notes were read for context and
comparison; existing submission files are untouched.
