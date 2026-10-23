# bls12-381-g1-ecdlp: Cheon subgroup walk with fixed-base comb tables

Claim: for s uniform in Z_r^* and inputs T = [s]P0, A_i = [s^i]P0 (i = 1..n, n = 2^20),
the algorithm below recovers s with probability at least **0.40**, using at most
**2^120.945788** target-group-operation equivalents at n = 2^20 (score 30236447/250000).
The total includes preprocessing, distinguished-point storage and lookup, trail starts and
restarts, the stopping rule's full budget, the square-root chain, verification and output.
Peak memory is at most 2^34.014135 bytes. Preprocessing is at most 2^21.032477 units, and
nonuniform advice is 0 bytes.

This is a classical cold attack on one target. It is not a practical break. The success
bound rests on the disclosed heuristic H1 (random-mapping behaviour of the walk) and on the
fresh toy runs in section 6. The per-step instruction count is a Callgrind measurement of
the full walk step as run (section 4).

Files: `accounting.py` recomputes every number here with exact integers (uses the vendored
`exact.py`; output in `evidence/accounting.json`). `experiment.py` is the bundled toy experiment
(stdlib only, from `cscore.cheon`). `evidence/` holds the fresh toy runs, the earlier toy runs
and the step measurement. `claim.json` and `submission.json` follow the track's schema.

## 1. Reduction

Let r be the group order, d = 2^20 (which divides r−1 since r−1 = 2^32·odd), and
H = {x^d : x in Z_r^*} the subgroup of d-th powers, |H| = N = (r−1)/d. Since
s is uniform in Z_r^*, h = s^d is uniform in H.

- A_d = [h]P0 with h in H. This is a discrete logarithm in H with base P0.
- Given h, recover s by twenty square roots in Z_r. With t_20 = h, for j = 19 down to 0 take
  a root t_j of t_(j+1) with [t_j]P0 = A_(2^j). The true s^(2^j) is one of the two roots, and
  the other fails the check. Finally s = t_0, checked against T = A_1.

The chain costs 41 comb evaluations and twenty Tonelli–Shanks square roots in Fr. Since
2-adicity(r−1) = 32, one square root costs at most 2·32² modular multiplications
(section 4). The chain is therefore not negligible in the bound, and it is charged.

## 2. The walk on H modulo mu6

mu6 = {1, −1, λ, −λ, λ², −λ²}, where λ is a cube root of unity mod r. Since 6 | N, mu6 lies
in H. On G1, −1 is negation and λ acts as the endomorphism (x, y) → (βx, y). Let
M = N/6 = 2^232.272 be the number of classes.

**Canonical form and key.** Among the six images of a point P under mu6, the canonical
representative has the least x-coordinate (comparing Montgomery representations), with y
normalised to even parity. A canonical point is determined by its x-coordinate: y is recovered
from x up to sign, and the parity fixes the sign. The **key** of a canonical point is its
x-coordinate in Montgomery form, a 381-bit integer. Two canonical points are equal exactly when
their keys are equal.

**Multipliers.** Sixty-four elements g_1..g_64 of H, each x^d for a random x in Z_r^*. Because
d = 2^20, x^d is twenty squarings mod r, which is charged in section 4.

**Step.** For a canonical point with key k and tracked exponent c:
- j = (k ≫ 90) mod 64 selects the multiplier;
- c' = c·g_j mod r;
- the new point is [c']B, evaluated by a window-16 signed fixed-base comb on B (P0 for tame
  trails, A_d for wild trails), then canonicalised, and c ← μ·c' for the canonicalising mu6
  element μ.

Both families apply the same rule and the same class function, so a tame and a wild walk that
reach the same canonical point merge.

**Trails and records.** A trail is a sequence of steps from a random start c_0 in H. Each
trail is capped at 2^100 steps and then restarted from a fresh random start. A point is
**distinguished** when its key has its low 90 bits equal to zero, which happens with
probability 2^−90 per point. Each distinguished point is stored as a record of 72 bytes: the
high 291 bits of the key (the other 90 bits are zero by definition, so the key is recovered
exactly; 37 bytes), the exponent c (32 bytes), the family (1 byte), and padding. Lookup of a
distinguished point compares the stored 291-bit prefix with the new key. A match is therefore an
exact key match, with no fingerprint and no false match: the omitted bits are zero for every
distinguished point. The walk does not infer point equality from any shorter fingerprint.

**Collision.** A tame point with exponent c_t and a wild point with exponent c_w sharing a key
satisfy [c_t]P0 = [c_w]A_d = [c_w·h]P0, so h = c_t·c_w^{−1} mod r exactly. The candidate is
verified by [h]P0 = A_d before the chain runs. Tame–tame and wild–wild matches carry no
information about h and are ignored.

**Stopping rule.** The walk stops after a fixed total of B = S(τ) + 2^100 points over both
families, where S(τ) = ceil(2·sqrt(τ·M)) and τ = 13/20. It returns the first verified
collision. Section 3 gives the loss from this rule.

## 3. Success probability

Let S be the number of walk points, tame and wild, performed before a collision is detected.

**H1 (random mapping, disclosed).** The class walk behaves like a random mapping on the M
classes, so the S points contain a tame–wild collision with probability
F(S) = 1 − exp(−S²/(4M)).

**H2 (detection, cycles and caps).**
- *Detection delay.* When two walks merge, they then follow the same class sequence, so the first
  distinguished point on the merged path is reached by both trails. The delay G to that point is
  geometric with mean 2^90 steps per trail, and the pair consumes about 2G points in total. A
  collision at point count t ≤ S(τ) is detected within the budget B if 2G ≤ B − t. With
  B − t ≥ 2^100, this fails with probability at most exp(−2^100/(2·2^90)) = e^−512.
- *Trail caps.* A trail that reaches its 2^100-step cap without a distinguished point has
  probability at most e^−1024 of doing so (no distinguished point in 2^100 steps), and it is
  restarted.
- *Cycles.* A cycle inside one trail has probability at most (2^100)²/(2M) ≈ 2^−33.3 per trail.
  Over the 116,127 trails of the run this is at most 2^−16.4.

So the success probability is at least (1 − e^−512)·F(S(τ)) − 2^−16.4 − e^−1024 under H1, with
F(S(τ)) = 1 − e^−τ = 0.477954 for τ = 13/20. The losses from the stopping rule, caps and cycles
are below 2^−16.

**Claim value.** The model value is above the claimed bound by about 0.078, and the observed rate
at the same budget is lower (section 6). The claimed lower bound 0.40 is the Wilson 95% lower
bound of the fresh toy runs at S(τ), which is 0.4227 over 4000 trials (section 6), minus the
losses above. The claim does not rely on the model value.

**Budget choice.** An earlier budget τ = 3/5 gave a Wilson lower bound of 0.3986 on the fresh
runs, which supports only 0.39. The budget was raised to τ = 13/20 for a claim of 0.40.

## 4. Computation

Work is W = H + R/C with C = 2000 (G1 price, `tracks/docs/COST_MODEL.md`). H counts batched
incomplete affine additions and doublings (reference units, one each, with the 6M convention
that includes the batch inversion). R counts all other scalar AArch64 instructions, including
loads, stores and control.

**Per walk step, measured.** `evidence/measure-full-step/cheon_step_full.c` runs the full walk
step as the walk performs it, over 256 walk records in round-robin order (half tame on the P0
table, half wild on the A_d table). Each step does, with every executed instruction counted:
the family table selection, the multiplier index and exponent multiply, the 16 signed digits
with one table entry copied per window and conditional negation, the 16 batched additions
represented by `ref_add` (a copy, which is an upper bound on the bookkeeping around one
addition; the addition itself is a reference unit), the canonicalisation (β·x, the third
rotation, least-of-three, parity), the exponent update, the walk-record store, the
distinguished-point test and the trail counter. A distinguished point calls a store routine
that is not reached in the measured range; trail restarts are not reached either. Both are
charged analytically below.

| Quantity (Callgrind Ir, organizer harness) | Value |
| --- | ---: |
| 256 steps | 714,288 |
| 512 steps | 1,427,632 |
| per step (differenced) | **2786.5** |
| of which the reference functions (`ref_add`, `ref_fp_mul`, `ref_fr_mul`) | 925.0 |

`ref_fr_mul` and `ref_fp_mul` are counted in R, because only the additions are reference units.
The kernel's own reference pass names the functions with `noipa`, so the toggle captures all
three. Each call of `ref_fr_mul` costs at most 925/2 = 462 instructions, since two calls per
step fit inside the 925 of the reference functions.

**Native and replica checks.** The kernel was built natively and compared with
`replica_full.py` over 512 and 2000 steps; both dumps (all 256 records) match exactly. The
container dump for 512 steps matches the replica too.

**Analytic bounds (not executed in the measured range).**
- *Distinguished points.* Each distinguished point costs at most DP_ACTION = 4000 instructions
  (a hash of the 37-byte key, a probe of at most a few slots, a 72-byte record store, and a
  key comparison). The expected count is B/2^90 = 2^26.825, so the charge is 2^26.825 × 4000.
- *Trail starts.* A start samples c_0 (at most 4000 instructions for the random element) and
  computes twenty squarings (each at most 462 instructions), then one comb evaluation (16
  additions). A start also costs at most one measured step's worth of instructions. There are
  116,127 starts.
- *Random sampling.* At most 4000 instructions per random element, covering a PRNG squeeze and a
  rejection test. This is a generous analytic bound, not a measurement.
- *Square roots.* Each Tonelli–Shanks root in Fr with 2-adicity 32 takes at most 2·32²
  modular multiplications, each at most 462 instructions. Twenty roots are charged.
- *Inversion for the collision solve.* At most 2^18 instructions.

**Preprocessing.** Preprocessing is everything before the walk starts: the comb tables for P0
and A_d, the 64 multipliers, the zeroed distinguished-point table, and the input parse.

| Item | H (units) | R (instructions) |
| --- | ---: | ---: |
| comb tables: batched additions, 2 × 16 × 2^15 | 2^20 | — |
| comb tables: base shifts, 2 × 15 × 16 doublings | 480 | — |
| comb tables: loads and stores of 192-byte entries (2 × 24 words each) | — | 50,331,648 |
| multipliers: 64 × (sampling + 20 squarings) | — | 847,360 |
| distinguished-point table zeroing (237,826,082 slots × 9 words) | — | 2,140,434,738 |
| input parse (T and the powers used) | — | 65,536 |

Preprocessing totals H = 1,049,056 and R/C ≈ 1,095,840 units, so its log2 is 21.032477.

**Walk and finish.**

| Item | H (units) | R (instructions) |
| --- | ---: | ---: |
| walk points: 16 additions each, B = S(τ) + 2^100 points | 16·B = 2^120.825 | — |
| walk points: measured step, 2786.5 each | — | 2786.5·B |
| distinguished points: 2^26.825 × 4000 (bound) | — | 475,652,164,000 |
| trail starts: 116,127 × 16 additions | 1,858,032 | — |
| trail starts: sampling and 20 squarings | — | 1,537,521,480 |
| trail starts: one step's worth each (bound) | — | 323,587,886 |
| chain: 41 comb checks | 656 | — |
| chain: one step's worth per check | — | 114,246 |
| chain: 20 Tonelli–Shanks roots (bound) | — | 18,923,520 |
| collision solve (bound) | — | 262,144 |

Totals: H = 2^120.825332, R = 2^128.269570 instructions, R/C = 2^117.303786, and
W = H + R/C, so **time_log2 = log2 W(2^20) = 120.945788** (= 30236447/250000).

**Range of n.** The bound uses A_1, A_(2^j) for j ≤ 20 and A_d = A_(2^20). It therefore holds for
every n with 2^20 ≤ n < r, with the same W, and the scored value is W(2^20). For n < 2^20 nothing is
claimed.

**Accounting convention.** Every item above is charged once. Loads, stores and control are
counted as instructions wherever they execute; table memory traffic beyond these counts is not
charged, since the score counts instructions. The 16 comb additions of each walk point, and
the doublings of the table build, are reference units under the 6M convention, which includes
batch inversion. The attacker's verification is one comb evaluation.

## 5. Memory and advice

The peak storage is the sum of four items:

| Item | Size (bytes) |
| --- | ---: |
| distinguished-point table: 2·2^26.825 slots (load factor 2) × 72 bytes | 17,123,477,904 |
| comb tables: 2·16·2^15 affine points × 192 bytes | 201,326,592 |
| trail state: 116,127 trails × 64 bytes | 7,432,128 |
| code, stacks and scratch | 16,777,216 |
| **total** | **17,349,013,840 (2^34.014135)** |

Nonuniform advice is 0 bytes, since the comb tables are built from the challenge's public points
at run time and charged in section 4.

## 6. Evidence

**Scope.** The toy (`experiment.py`, from `cscore/cheon.py` in curvesmash-core) runs the full walk
and the chain on a generic prime-order group. Points are opaque keyed labels: the attacker sees
equality, multiplication by known integers and the mu6 action, and never the exponent. It is not
a BLS12-381 curve, so it tests the collision statistics and the correctness of the reduction, not
the per-step cost.

**Fresh runs (the basis of the claim).** `evidence/toy-fresh/` holds 4000 fresh trials, with the
per-trial step counts and the driver `toy_stats.py`. All 4000 secrets were recovered and no
collision gave a wrong h.

| Batch | bits, k | trials | successes at S(13/20) | rate | Wilson 95% |
| --- | --- | ---: | ---: | ---: | --- |
| b32-k6-3000-s4101 | 32, 6 | 3000 | 1314 | 0.438 | [0.4203, 0.4558] |
| b36-k10-1000-s4102 | 36, 10 | 1000 | 438 | 0.438 | [0.4075, 0.4689] |
| **pooled** | | **4000** | **1752** | **0.438** | **[0.4227, 0.4534]** |

Success at other budgets (pooled, same 4000 trials):

| τ | model 1 − e^−τ | successes | rate | Wilson 95% |
| --- | ---: | ---: | ---: | --- |
| 3/5 | 0.4512 | 1655 | 0.4138 | [0.3986, 0.4291] |
| 31/50 | 0.4621 | 1703 | 0.4258 | [0.4105, 0.4411] |
| **13/20** | **0.4780** | **1752** | **0.4380** | **[0.4227, 0.4534]** |
| 7/10 | 0.5034 | 1846 | 0.4615 | [0.4461, 0.4770] |

The claim of 0.40 is supported at τ = 13/20 by the lower bound 0.4227. At τ = 3/5 the lower bound
0.3986 supports only 0.39, which is why the budget was raised. The observed rate is about 0.04
below the model at every τ, which is the finite-size gap that the claim absorbs.

**Earlier runs.** `evidence/toy-initial-b32-k6-60.json`, `toy-k6-400.json`, `toy-k8-400.json`,
`toy-b36-k10-150.json`, `toy-b38-k10-200.json` and `toy-mult16-b32-k6-400.json` were run before the
fresh batch, at τ = 3/5 with 1210 trials on the 64-multiplier settings. Their pooled rate, 503/1210
= 0.416, has a Wilson 95% interval of [0.388, 0.444] (not [0.40, 0.43] as earlier stated). They
agree with the fresh runs at τ = 3/5 and are kept for transparency; the claim does not use them.
The 16-multiplier ablation (`toy-mult16-b32-k6-400.json`) gives 0.383 at τ = 3/5, which is why the
package uses 64 multipliers.

**Declared experiment.** `toy-cheon-walk` (initial input 32 bits, k = 6, 60 trials, seed 2601)
recovers 60/60 and reports its own budget statistics at τ = 3/5 and τ = 0.39. The claim budget
τ = 13/20 is evidenced by the fresh runs above, not by the declared experiment.

**Step measurement.** Section 4 and `evidence/measure-full-step/`: the kernel, the replica, the
input generator, the run log, the summary, and the replica's stderr. The input file (402,218
bytes, SHA-256 `f9d7dda243ff41a6842a09d56c0e43577a2f4e4876f8989f84927f119cc45c7f`) and the
container dump (70,656 bytes, SHA-256 `72d5fb76cc43e6219b6de4cc4a150afec721db4e0e6308fcccf173670b7f0649`)
exceed the 64 KiB per-file limit and are not in the package. `gen_inputs_full.py <target-profile.json> 7`
regenerates the input file byte for byte.

## 7. Objections from the review of submission 8798126c, and how they are answered

1. **Preprocessing was charged as 2^20 additions with no other work.** Section 4 charges every
   preprocessing item: the table loads and stores, the multiplier sampling and exponentiation, the
   zeroing of the distinguished-point table, and the input parse. The total is 2^21.032477.
2. **The 2703.56-instruction kernel was a surrogate, and the 100-instruction allowance had no
   support.** The measured kernel now runs the walk loop over 256 records, both families, the
   record store, the trail counter and the distinguished-point test, with additions as unit copies.
   The allowance is removed. The analytic items that the loop does not reach (DP storage, restarts,
   sampling, Tonelli–Shanks, the solve) carry explicit bounds in section 4.
3. **The memory ledger reserved 8 bytes for the label, which could not hold the key.** Each
   distinguished point is stored as an exact 291-bit high part of the key plus the exponent, 72
   bytes in all. No fingerprint is used, and a match is an exact key match. The memory total rises
   to 2^34.014135 bytes.
4. **H2 had an unbounded "at most" detection allowance.** The walk now runs a fixed budget
   B = S(τ) + 2^100. The loss from the stopping rule, the trail caps and the cycles is bounded in
   section 3 and is below 2^−16.
5. **The pooled interval was misstated.** Section 6 gives the earlier pooled interval as
   [0.388, 0.444], and the claim now rests on 4000 fresh trials at the claim budget, with a Wilson
   lower bound of 0.4227.

## 8. Limitations and heuristics

- **H1** is the random-mapping heuristic for the class walk. Its support is the fresh toy runs,
  the Teske r-adding walk analysis and the random-mapping argument, not a proof. The toy groups
  have orders of 32 to 38 bits, far below the production class count 2^232.27. The claim keeps a
  margin of 0.0227 above 0.40 at the lower Wilson bound.
- **H2** is the detection and trail model of section 3, which uses the merged-path argument for
  detection delay. It is not tested by the toy runs at the production delay scale (the toy uses
  2^3-point distinguished points).
- **Measurement scope.** One compiler (GCC 14.2.0), one pinned blst, Valgrind 3.24.0, 256
  differenced steps on one input set. The 512-step container dump matches the replica.
- **Analytic bounds.** The distinguished-point, sampling, restart and Tonelli bounds are analytic
  instruction bounds, not measurements. Together they are about 9 × 10^−29 of W, so even a large
  error in any of them leaves the printed score unchanged.
- **Memory is about 17 GB** (2^34.014 bytes). The earlier 10 GB figure belonged to the rejected
  package and is superseded.
- **Window trade-off.** Larger comb windows trade memory for steps at the same budget (the table
  is for the earlier τ and is not recomputed here): window 20 needs about 2.6 GB of comb tables,
  window 24 about 35 GB. This candidate uses window 16.

## 9. Reproduction

```sh
python3 -B -I accounting.py
echo '{"bits":32,"k":6,"trials":60,"seed":2601}' | python3 -B -I experiment.py
```

The fresh toy runs: `CURVESMASH_CORE=<curvesmash-core> python3 -B evidence/toy-fresh/toy_stats.py 32 6 3000 4101`
and `... 36 10 1000 4102`. The step measurement uses the organizer's reference-count image (Docker
linux/arm64) on `evidence/measure-full-step/cheon_step_full.c`, with `N1=256 N2=512` and
`REF_FUNCS="ref_fp_mul ref_fr_mul ref_add"`. The replica is `replica_full.py <target-profile.json>
inputs.txt 512 dump-512.txt`.
