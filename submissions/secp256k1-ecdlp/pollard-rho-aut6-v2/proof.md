# secp256k1-ecdlp: Pollard rho with order-6 automorphisms, honest accounting (v2)

Status: ready for independent review; no qualification or emitted score.
This is a revision addressing the review of submission
`pollard-rho-aut6` (rejected for unsupported per-step cost, vague
cycle-escape, and non-peak memory bound). All three are fixed here with
explicit constructions and conservative bounds. This remains an analytic
cryptanalytic claim under exploratory review, not a practical execution.

The normative inputs/outputs are those of `production-v2`; `claim.json`
binds that profile and `classical-group-work-v2` by SHA-256. No restricted
challenge distribution, advice, random-oracle premise, quantum operation or
unpriced oracle is used. All work is charged under the fixed cost model
(C = 1000, reference = batched incomplete affine secp256k1 addition, 6M).

## 1. Target and automorphism group

Same as v1 (see proof.md in `pollard-rho-aut6` for the full setup; we
repeat the essentials). p = 2^256 - 2^32 - 977, n prime as in the
protected profile, E: y^2 = x^3 + 7, generator P, input Q = [x]P.

G_aut = {1, lambda, lambda^2, -1, -lambda, -lambda^2}, order 6, where
psi(x,y) = (beta*x, y), beta^3 = 1 mod p, beta != 1, and psi acts on <P>
as multiplication by lambda with lambda^2 + lambda + 1 = 0 mod n.

## 2. Algorithm: batched parallel Pollard rho on orbits

### 2.1 Parallel walks with batched affine additions (the key construction)

We run P = 2^20 parallel walks (P fixed below). Each walk holds
(R_i, a_i, b_i) with R_i in affine coordinates (x_i, y_i), and the
invariant R_i = a_i*P + b_i*Q as group elements.

**Per round** (all P walks advance one step together):

1. **Partition**: for each walk w, j_w = s(R_w) in {0..31} from the
   affine x-coordinate (low 5 bits after a fixed public mixing;
   10 instructions).

2. **Batch the P affine additions** R_w + T_{j_w} via Montgomery's
   simultaneous inversion:
   - For each w, compute the denominator d_w for the affine addition
     slope: d_w = (x_{T} - x_w) mod p, or for doubling (negligible
     probability; handled by a fallback below) d_w = 2*y_w.
     Cost per w: one modular subtraction (~256 instructions? No—
     modular subtraction is ~256 word ops in the baseline's model,
     but in straight AArch64 it's ~8 instructions for 4 limbs plus
     conditional correction, say 20 instructions). We bound the
     non-inversion part of each affine addition at 3M = 537
     instructions (see below).
   - Batch-invert (d_1, ..., d_P): Montgomery's trick costs
     3*(P-1) field multiplications + 1 field inversion.
     Per walk amortized: 3M + (1 inversion)/P.
   - Complete each addition: slope = num * d_w^{-1}, then
     x' = slope^2 - x_w - x_T, y' = slope*(x_w - x') - y_w.
     This is the standard incomplete affine addition, ~3M beyond the
     inversion (which we already counted via batching).

   **Cost per walk-step**: 3M (pre-inversion) + 3M (amortized batch) +
   3M (post-inversion) = 9M? No—let's be careful.

   Actually, the reference "batched incomplete affine addition, 6M"
   is *defined* as: 3M for the addition formula + 3M for batch
   amortization. So one batched affine addition = 6M = 1000
   instructions = 1 reference unit, *by definition of the reference*.

   Our construction realizes exactly this: P parallel affine additions,
   batched via Montgomery, each costing 6M amortized (for P large
   enough that the single inversion is negligible; with P = 2^20,
   inversion/P = 68600/2^20 ≈ 0.065 instructions, negligible).

   **Therefore each walk-step's addition costs exactly 1 reference
   unit.** This addresses the v1 review objection: we now *specify*
   the batching construction that realizes the reference price, rather
   than asserting it.

   Edge cases (doubling, infinity, opposite points): the partition
   steps T_j are fixed public points; R_w is (heuristically) uniform.
   P(R_w = ±T_j or infinity) is negligible (< 2^-250 per step). We
   handle these by a deterministic fallback (skip the addition, take
   a fixed extra step), charged at 1 unit. This does not affect the
   asymptotics and is included in the cycle-escape overhead.

3. **Canonicalize** each R'_w to rep(R'_w) (lex smallest of 6 images):
   - Compute x1 = beta*x mod p, x2 = beta*x1 mod p: 2 field mults.
     Using the pinned 179-instruction multiply: 358 instructions.
   - Find min x among {x, x1, x2} (negations share x): 2 comparisons
     of 256-bit ints, ~32 instructions. Then min y (y vs p-y):
     ~16 instructions. Total ~48.
   - Coefficient scaling: (a,b) -> mu*(a,b) mod n. mu in
     {±1, ±lambda, ±lambda^2}. Multiplying by ±1 is free (conditional
     negation). Multiplying by lambda or lambda^2 mod n: one 256-bit
     modular multiplication each. Generously 300 instructions each
     (the baseline's M is for field mult; scalar mult mod n is
     similar). For 2 coefficients: 600 instructions. But we only pay
     this when mu is not ±1 (probability 2/3). Expected: 400.
     For an upper bound, charge 600 always.
   - Total canonicalization: 358 + 48 + 600 = 1006 instructions
     = 1.006 units. **We round up to 1.1 units** for safety.

   *Note: this is the dominant overhead. The v1 submission used 0.7;
   the reviewer correctly demanded justification. The 1.1 here is
   itemized and conservative.*

4. **Distinguished check**: t = 40 leading zero bits of affine x.
   ~20 instructions = 0.02 units. If distinguished, store
   (rep, a, b) in the shared hash table (see Section 4 for the
   peak-memory bound). Store cost amortized negligible.

5. **Cycle escape** (concrete rule, addressing v1 objection):
   - After computing rep_{w,i+1}, compare to rep_{w,i} (1-cycle) and
     rep_{w,i-1} (2-cycle). Each compare ~48 instructions, total
     ~96 = 0.096 units, **rounded to 0.1**.
   - If a 1-cycle or 2-cycle is detected, *replace* the next step for
     that walk with: R := R + E, where E is a fixed public point
     (e.g. E = [3]P), then canonicalize normally. This costs 1 extra
     batched addition (already counted in the per-round batch) plus
     the canonicalization (counted). No extra asymptotic cost; the
     walk is perturbed deterministically and continues.
   - Longer fruitless cycles (length ≥ 3): these require the walk to
     return to a previous rep via a non-trivial path. Under the
     random-walk heuristic with r = 32, the probability per step of
     entering a k-cycle for k ≥ 3 is O(1/N_orb), negligible. We do
     not add a special rule; the walk's randomness (via the partition
     function on fresh reps) breaks such cycles with overwhelming
     probability. This is disclosed as a heuristic (Section 5).
   - **Cost accounting**: the 1-cycle/2-cycle check (0.1 units) is
     charged *every step*. The escape step itself replaces a normal
     step (no extra cost beyond the already-counted addition). We add
     0.1 units/step for the checks.

**Total per walk-step**: 1.0 (batched affine addition) + 1.1
(canonicalization) + 0.02 (distinguished) + 0.01 (partition) + 0.1
(cycle checks) = **2.23 units**. We use **2.25** as a clean upper bound.

### 2.2 Walk parameters and collision use

- r = 32 partitions, P = 2^20 parallel walks, t = 40 distinguished bits.
- Precompute T_j = c_j*P + d_j*Q for j = 0..31 (32 scalar mults,
  ≤ 2^15 units, negligible; charged in preprocessing).
- Each walk starts at random (a_{w,0}, b_{w,0}); R_{w,0} =
  rep(a_{w,0}*P + b_{w,0}*Q) via one scalar mult each (charged;
  P * 512 additions = 2^20 * 2^9 = 2^29 units total, negligible
  vs 2^128).
- On a table collision (same rep, different (a,b)): solve
  x = -(a1-a2)/(b1-b2) mod n, verify [x]P = Q. If b1 = b2 (negligible
  probability), continue.
- The algorithm is Las Vegas (terminates w.p. 1); we convert to fixed
  budget via Markov below.

## 3. Work analysis and success probability

### 3.1 Heuristic

Same as v1: the r-adding walk on orbit representatives is modeled as a
uniform random walk on a set of size ≤ n/6, giving expected steps to
first collision E[T] ≤ sqrt(pi*n/12). Supported by the Pollard-rho
literature (Teske; van Oorschot–Wiener; Bos et al. for automorphisms).
Disclosed in Section 5. The P parallel walks contribute to a shared
table; the *total* steps across all walks to a collision is what
matters, and the birthday bound applies to the aggregate.

*Clarification on "distinguished-collision time" (v1 objection):* with
distinguished points, we need a collision *among stored (distinguished)
points*. The expected total steps is E[T] * (1 + ε), where ε accounts
for the fact that the colliding pair must both be distinguished. For
2^-t distinguishing probability, standard analysis (van
Oorschot-Wiener, Section 3) shows ε is negligible when t is chosen so
that 2^t << sqrt(E[T]) (here 2^40 << 2^64). We absorb this in the
explicit 2.25x per-step cost and the Markov 2x budget, and disclose it.

### 3.2 Fixed budget (Markov)

B_total = 2 * E[T] total walk-steps across all P walks.
P(T > B_total) ≤ 1/2, so success ≥ 1/2 (minus negligible bad-collision
probability, giving ≥ 0.49 as in v1).

Using n < 2^256:
  E[T] < sqrt(pi * 2^256 / 12) = 2^128 * sqrt(pi/12).
  sqrt(pi/12) = 0.5116629... (computed below).
  log2(sqrt(pi/12)) = -0.966733...
  E[T] < 2^127.033267.
  B_total < 2^128.033267 walk-steps.

Each walk does B_total / P steps; the *total work* sums across walks
(per the cost model: "Sum work across processors").

### 3.3 Total work and score

W = B_total * 2.25 + W_pre + W_post, where:
- W_pre: 32 T_j + P initial scalar mults ≤ 2^15 + 2^29 < 2^30.
  (P = 2^20 walks * 512 additions each = 2^29.)
- W_post: collision solving (1 inversion) + verification (1 scalar
  mult) < 2^20.
- Both negligible vs 2^128.

W < 2^128.033267 * 2.25 + 2^30
  < 2^128.033267 * 2^1.169925 + 2^30    (log2(2.25) = 1.169925...)
  = 2^129.203192 + 2^30
  < 2^129.203193.

**Score = log2(W) < 129.203193.**

We declare **129.2032** (rounded up), i.e. 1292032/10000 = 80752/625.

*Honest comparison:* this is *worse* than the BSGS baseline
(129.000802) by about 0.20 bits. The v1 submission's 128.84 was
achieved by undercharging the per-step cost; honest accounting for
affine-coordinate canonicalization (1.1 units) and cycle checks (0.1)
erodes the rho+automorphism advantage. We report this transparently:
**this revision does not beat the baseline.** It is submitted as a
corrected, defensible accounting of the rho+automorphism approach,
not as an improvement.

*Why submit a non-improving score?* Per Thang's instruction to revise
with honest accounting. The value is in the corrected methodology
(explicit batching construction, concrete cycle-escape, peak memory
bound), which future work can build on. We do not claim a leaderboard
improvement.

### 3.4 Success probability

≥ 0.49, as in v1 (Markov 1/2 minus negligible bad-collision term).
Meets the 0.39 threshold.

## 4. Peak memory bound (deterministic, addressing v1 objection)

The v1 submission reported an *expected* table size (2^95 bytes). Here
we give a *deterministic worst-case* peak bound.

- The table stores at most one entry per walk-step that yields a
  distinguished point. In the worst case, *every* step could be
  distinguished (probability 2^{-t*B_total}, negligible, but we bound
  deterministically).
- Worst-case entries: B_total < 2^128.04.
- Per entry: rep (64 bytes: x,y) + a (32) + b (32) = 128 bytes
  (we use 128, not 96, for alignment conservatism).
- Worst-case table: 2^128.04 * 2^7 = 2^135.04 bytes.
- Walk states: P = 2^20 walks * 192 bytes (R,a,b + scratch) = 2^27.6.
- Code + constants + batch buffers: < 2^30 bytes (generous).
- **Deterministic peak: < 2^135.1 bytes. We report
  memory_log2_bytes = 136** (rounded up).

*Note:* the *expected* table size is B_total/2^40 ≈ 2^88.0 entries ≈
2^95 bytes, but we report the deterministic worst-case 2^136 as
required. Memory is not scored. This bound is extremely loose; a
tighter high-probability bound via Chernoff would give ~2^96, but we
prioritize determinism per the review.

*Pruning policy (for completeness):* we do not prune; the worst-case
bound above assumes no pruning. If memory were a constraint, one could
evict oldest entries at the cost of a quantified success-probability
reduction; we do not need this.

## 5. Heuristics disclosed

1. **Random-walk heuristic** (same as v1): r=32 adding walk on orbits
   ≈ uniform random walk; E[T] ≤ sqrt(pi*n/12). Standard, literature-
   supported. The 2x Markov margin absorbs small deviations.
2. **Automorphism quotient**: effective space n/6; smaller orbits only
   help. Conservative.
3. **Fruitless cycles**: 1- and 2-cycles handled by the explicit
   escape rule (Section 2.1, step 5) with charged cost (0.1/step for
   checks; escape replaces a normal step). Longer cycles (k≥3) are
   heuristically negligible (probability O(1/N_orb) per step;
   supported by Bos et al.). Disclosed.
4. **Distinguished-collision overhead**: the (1+ε) factor for needing
   *stored* collisions is negligible for t=40 << 64; absorbed in
   Markov 2x. Disclosed (this addresses the v1 "ideal random-map"
   objection by naming the ε and bounding it).
5. **Coefficient usability**: b1≠b2 with overwhelming probability;
   standard.
6. **Batching efficiency**: Montgomery batch inversion achieves the
   3M amortized cost for P=2^20 (inversion/P negligible). This is
   textbook, not heuristic, but we list it for completeness.

No subexponential conjectures, oracles, or unproven number theory.

## 6. Limitations

- Analytic claim; no full-size execution.
- **Does not beat the BSGS baseline** (129.2032 vs 129.000802).
  The rho+automorphism step savings are outweighed by honest per-step
  canonicalization and checking costs under this cost model.
- Score is an upper bound; tighter implementation constants could
  reduce it, but not below ~129.1 given the 1.1-unit canonicalization
  floor.
- Memory bound 2^136 is a loose deterministic worst-case; expected
  is ~2^95.
- Success 0.49 is a lower bound.

## 7. What changed from v1 (for the reviewer)

| Issue (v1 review) | Fix in v2 |
|---|---|
| Per-step cost used batched reference price with no construction | Section 2.1 specifies P=2^20 parallel walks + Montgomery batch inversion, realizing exactly 1 unit/addition; canonicalization itemized at 1.1 units |
| Cycle escape was "a deterministic extra step" with no rule | Section 2.1 step 5: concrete 1-/2-cycle detection (charged 0.1/step) and escape via fixed point E; k≥3 cycles heuristically bounded |
| Memory was expected size, not peak | Section 4: deterministic worst-case 2^136 bytes, no pruning assumed |
| Score 128.84 unsupported | Honest score 129.2032, transparently worse than baseline; submitted per author instruction as a corrected accounting |

## 8. Numerical verification of constants

We verify the key constants by direct computation (Python, exact):

- pi = 3.141592653589793
- pi/12 = 0.2617993877991494
- sqrt(pi/12) = 0.511662902... (see below)
- log2(sqrt(pi/12)) = -0.966733...

To reproduce:
```python
import math
v = math.sqrt(math.pi/12)
print(v)                    # 0.511662902...
print(math.log2(v))         # -0.966733...
print(math.log2(2.25))      # 1.169925...
print(128 + math.log2(v) + 1 + math.log2(2.25))  # 129.20319...
```

Score: 129 + (-0.966733) + 1.169925 = 129.203192 → declared 129.2032.

As a rational: 129.2032 = 1292032/10000 = 80752/625.
Check: 80752/625 = 129.2032. Yes.
