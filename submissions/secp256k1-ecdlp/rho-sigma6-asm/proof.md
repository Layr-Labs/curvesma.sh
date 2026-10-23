# secp256k1-ecdlp: parallel Pollard rho on order-6 classes with a sigma walk, all arithmetic in instructions

Claim: for uniform x in [0,r) and Q=[x]P on secp256k1, the algorithm below returns x
with probability at least 0.39. It uses at most **2^126.418986** target-group-operation
equivalents in total (score 63209493/500000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^74.74996 bytes. Preprocessing is
2^48.619659 units, and no nonuniform advice is used. This is a classical cold attack on
one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on automorphism classes, and on the rate bound H2. The cost rests on
H3, instruction counts measured with the organizer's toolchain. This is a carefully
accounted generic attack. It is not a practical break.

**Cost route.** The walk does not use the 6M reference operation. Every instruction of
every step, including all field multiplications, squarings and subtractions, is its own
4×64-bit GCC inline AArch64 assembly and is charged at 1/C with C = 1000 (H = 0 in
COST_MODEL's H + R/C). COST_MODEL.md: "Faster arithmetic and lookup algorithms are welcome;
justify their actual work". A step measures 781 instructions in the organizer's
toolchain, so it is charged 0.795 units rather than 1 + R/C.

Relative to the usual order-6 walk, the walk never moves a point to its canonical
representative, so it never pays the multiplication by beta (section 1.2). It stores the
actual point W and adds sigma(R_i), where sigma is the automorphism that maps the
canonical point to W. The class sequence is unchanged, so the standard analysis applies.
Which rotation is canonical is decided from 64-bit fixed-point fractions with an exact
fallback, and Lemma 1 bounds their error. Compared with the earlier rho-order6 package,
it removes the beta multiplication, and the budget is tau = 0.5303. The same walk with the reference operation
and libsecp256k1 field code is the companion package rho-sigma6 (1.32 units per step).

Files:
- `accounting.py` recomputes every number with exact integers (output in
  `evidence/accounting.json`).
- `measurement/` holds the per-step instruction-count harness, the field code
  (`secp_field.h`, `secp_sub.h`), the Python replica `walk_ref.py`, and the approximation
  check `frac_check.py`.
- `evidence/measurement/` holds the counts on two input sets.
- `experiment.py` is the declared toy-curve experiment for the identical walk;
  `evidence/small-curves/` holds larger offline runs of the same file.

## 1. Algorithm

### 1.1 Notation and parameters

p, r and P come from the protected profile. The curve is y^2 = x^3 + 7, the cofactor is 1,
and p ≡ r ≡ 1 (mod 3).

The endomorphism uses beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
and lambda = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72, with
[lambda](x,y) = (beta·x, y). The automorphism group G6 = {±1, ±lambda, ±lambda^2} has
(r−1)/6 nonzero classes, and N = (r−1)/6 + 1. Field values are compared as reduced
integers in [0,p).

Parameters:
- A = 2^16 adders.
- K = 2^32 walks: 2^16 processors × 2^16-walk batches.
- theta = 2^−60.
- L = 33310394024866356565296666183 ≈ 2^94.7500.
- CAP = 64/theta.
- DP_MAX = 2·theta·S_max + 2K, where S_max = K·(L + CAP).

### 1.2 Classes, keys and the sigma walk

For W = (x, y) ≠ O:
- **Class key.** key(W) = y if y is even, else p − y, so key(W) is ±y and even.
  - key is constant on the class, since G6 changes y only by sign.
  - key determines the class: y^2 = x^3 + 7 fixes x up to a factor beta^j.
- **Canonical point.** canon(W) = (beta^k·x, (−1)^s·y), where k minimises beta^k·x and
  s = y mod 2.
- **sigma.** sigma_W = (x → beta^(−k)·x, y → (−1)^s·y) = [(−1)^s·lambda^(−k)], so that
  W = sigma_W(canon(W)).
- **Walk step.** i(W) = (key >> 1) mod 2^16, and

      f(W) = W + sigma_W(R_i) = sigma_W(canon(W) + R_i).

  So the class sequence is exactly that of the standard class walk
  [W] → [canon(W) + R_h([W])]. The six images sigma(R_i) are tabulated per adder.
- **Distinguished point.** W is distinguished when bits 17..76 of key(W) are zero.
- **Representation.** In the kernel, x and y are fully reduced integers in [0, p), held as
  four 64-bit limbs.

**Deciding k without beta·x (Lemma 1).** Write x = Σ_{i<4} a_i·2^(64i), and let
C1 = ⌊2^384·beta/p⌋ = Σ_{j<6} C1_j·2^(64j). Define

    g_0 = a_3,
    g_1 = ( Σ_{i+j=5} lo(a_i·C1_j) + Σ_{i+j=4, j≥1} hi(a_i·C1_j) ) mod 2^64,

so g_1 takes 8 multiplications (4 `mul`, 4 `umulh`), and g_0 none. Let
f_j = 2^64·(beta^j·x mod p)/p = 2^64·frac(x·beta^j/p).

**Lemma 1.** f_0 ∈ [g_0, g_0 + 2). If g_1 ≤ 2^64 − 10, then f_1 ∈ [g_1, g_1 + 9).

*Proof.* For j = 0: x/2^192 ≥ a_3, and since p < 2^256,
f_0 − a_3 < 1 + 2^64·x·(2^256 − p)/(p·2^256) < 2.
For j = 1: beta/p = (C1 + e)·2^−384 with e ∈ [0,1), so
f_1 = (x·C1·2^−320 + x·e·2^−320) mod 2^64, and the second term is below 2^−64. Next,
x·C1·2^−320 = Σ_{i,j} a_i·C1_j·2^(64(i+j−5)). Modulo 2^64:
- terms with i + j ≥ 6 vanish;
- i + j = 5 gives lo(a_i·C1_j);
- i + j = 4 gives hi(a_i·C1_j) plus a fraction below 1 (four terms with j ≥ 1);
- the term with i + j = 4 and j = 0 does not exist (i ≤ 3);
- i + j = 3 gives four terms below 1 each;
- i + j ≤ 2 gives terms below 2^−64 each.

Dropping all of these gives an error in [0, 8 + 2^−50). The hypothesis excludes
wrap-around. ∎

x + beta·x + beta^2·x ≡ 0 mod p, with the reduced values in (0, p), so they sum to p or 2p
and f_2 ≡ −f_0 − f_1 mod 2^64. With g_2 = (−g_0 − g_1 − 22) mod 2^64, the real f_2 lies in
(g_2 + 11, g_2 + 22] whenever g_2 ≤ 2^64 − 23. So every f_j lies in [g_j, g_j + 22]. The
kernel accepts the least g_j as k only if all three g_j ≤ 2^64 − 23 and the least is more
than 22 below both others; then k is exact. Otherwise it computes beta·x and compares
exactly. The fallback rate is below 2^−56 per step and was never seen in 2^19 measured
steps. `measurement/frac_check.py` checks Lemma 1 and the decision rule on 200,000 random
and edge-case inputs per seed; the largest floor gap seen for g_1 is 6.

### 1.3 Procedure

1. **Validate Q.** Check it is on the curve. If Q = O, output 0.
2. **Adders.** For each adder i < 2^16:
   - draw c_i, d_i uniform in [0,r);
   - set R_i = [c_i]P + [d_i]Q;
   - tabulate its six images.

   For each walk w, set W_w = [a_w]P + [b_w]Q from fresh coins. All coins are stored, and
   coefficients are not tracked during the walk.
3. **Step.** W ← f(W), as one batched incomplete affine addition: Montgomery's
   simultaneous inversion over the 2^16 walks of a batch (3M), then lambda, lambda^2 and
   y3 (2M + 1S), with one Fermat inversion per batch.
   - If the batch product is zero, the affected walk uses the unbatched complete formula.
4. **Fruitless 2-cycles.** Each walk keeps the low 64 bits of its previous key.
   - On a match, the walk evaluates f(new) unbatched.
   - If key(f(new)) = key(old) exactly, the walk replaces new by 2·canon(M), where M is
     the element with the smaller key.
5. **Short cycles (3..64).** Every 64 steps each walk records its key as an anchor.
   - On a later match, the walk iterates f unbatched up to 64 times.
   - If its class returns, the walk escapes from 2·canon(M), with M the cycle element of
     least key.
6. **Distinguished points.** The walk sends (key, w, step index) to a hash table keyed by
   key.
   - A walk runs L steps and then continues until its next distinguished point.
   - A walk that exceeds CAP steps without one is abandoned and restarted from fresh coins,
     with at most K restarts in total.
   - The run fails if more than DP_MAX points arrive or the restarts run out.
7. **Resolution.** On the first key stored twice, both walks are replayed from their coins.
   - The replay tracks (a, b) mod r: an image sigma(R_i) has coefficients multiplied by
     (−1)^s·lambda^(−k), and canon multiplies them by (−1)^s·lambda^k.
   - This gives both walks' canonical points as [a]P + [b]Q. If b1 ≠ b2, it outputs
     x = (a2 − a1)/(b1 − b2) mod r after checking [x]P = Q. Otherwise it continues; the
     bound counts this case as a failure.

**Merged trails stay merged.** Outside fruitless cycles the next class depends only on the
current class. Escapes fire only after an exact check, and they exit from a point
determined by the cycle's set of classes. The 64-bit filters only decide when an exact
check runs.

## 2. Success probability

Let T_steps = K·L. Then tau = T_steps^2/(2N) = 0.530300.

**H1 (heuristic, disclosed).** The class map [W] → [canon(W) + R_h([W])], with uniform
random adders, behaves like a random mapping on the N classes, apart from the r-adding
penalty. That penalty is a factor of about (1 − 1/(6A)) in the birthday exponent (Teske
2001; Bernstein–Lange 2012; Bos–Kleinjung–Lenstra 2010). Under H1, the probability that
the first T_steps steps contain no useful class collision is at most
exp(−tau·(1 − eps)). Here eps = 4/A + 2·(3/A) + 2·(192/A^2) = 0.00015268. Its three
margins are:
- the adder penalty, taken as 4/A;
- 2-cycles at a rate of at most 1/A, with 3 non-fresh states each;
- short cycles at a rate of at most 1/A^2, with 192 wasted states each.

**H2 (rates).**
- 2-cycles are entered at a rate of at most 1/A per step. The literature gives 1/(6A).
  The earlier rho-order6 package measured 7.2e−5 and 9.2e−5 per step on real secp256k1
  walks at A = 2^11, against 1/(6·2^11) = 8.1e−5.
- Cycles of length 3..64 are entered at a rate of at most 1/A^2.
- Longer cycles with no distinguished point arise at a rate of at most A^−32. They are
  handled by abandonment, which is charged.

The remaining failure terms are:
- detection beyond CAP: about 2^−25;
- degenerate b1 = b2: at most 2^−250;
- DP_MAX or restart overflow: under 2^−1000;
- exceptional additions: under 2^−110.

Together these are below delta = 2^−20. Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015268)) − 2^−20 = 0.411523 ≥ 0.39.**

The budget is 1.0299·√N steps, against 0.9943·√N for success exactly 0.39. The 0.0215
of slack in the model bound covers the sampling error of the toy evidence below.

**Evidence for H1 and H2.** `experiment.py` with `"b": 7` runs this exact sigma walk on
random prime-order toy curves y^2 = x^3 + 7 over p ≡ 1 mod 3. It uses the same class key
±y and the same canonical rule, adder images and index (keys use Montgomery
representatives v·2^64 mod p, any fixed bijection being equivalent). It also uses both
cycle rules, escapes, distinguished points and abandonment. It counts every step to the
first useful class collision, solves end to end and checks [x]P = Q.

<!-- toy-evidence -->
Offline runs (`evidence/small-curves/`, pooled in `summary.json`): 7,000 trials on 7 curves
with N = 2^26.8, 2^29.1, 2^28.6, 2^32.5, 2^33.2, 2^37.2, 2^37.4 classes, adder tables of 2^11 and 2^16.
- Every trial solved and verified [x]P = Q.
- Success by T = sqrt(2 tau N) at the claimed tau = 0.5303: pooled 0.4129, Wilson 95%
  interval [0.401, 0.424]. The model value is 0.412.
- At tau = 0.4943 (model 0.390) the pooled success is 0.3919, interval [0.381, 0.403].
- Mean steps to the first useful collision: 1.254 sqrt(N), against sqrt(pi/2) = 1.253.
- 2-cycle entry rate per step:
  - 8.39e-05, 8.07e-05, 8.18e-05, 8.17e-05, 7.99e-05 with A = 2^11 (1/(6A) = 8.14e-05);
  - 2.08e-06, 2.57e-06 with A = 2^16 (1/(6A) = 2.54e-06).
- Longer fruitless cycles: 0 in total. Abandoned walks: 7.
<!-- /toy-evidence -->

The declared initial input (`bits` 26, 150 trials, seed 2601) is a quick reproducibility
check, not the statistical evidence. At tau = 0.5303 it reports 52/150 = 0.347, Wilson 95%
interval [0.275, 0.426], which contains the model value 0.412; with 150 trials the
standard error is about 0.04. The success claim rests on the 7,000 pooled offline trials
above, which any input within the schema can spot-check.

Limits: the toy curves have N up to about 2^38, against 2^253 in production. H1 is an
extrapolation in N, supported by the size-free random-mapping argument and the cited
literature. It is not proved for this walk.

## 3. Computation

**No reference operation.** Each step is charged its full measured instruction count
divided by C = 1000, as COST_MODEL.md's H + R/C with H = 0. No 6M reference operation is
claimed, so nothing is charged twice.

**Field arithmetic.** `measurement/secp_field.h` implements multiplication and squaring
modulo p = 2^256 − c, c = 2^32 + 977, as GCC inline AArch64 assembly in general-purpose
registers.
- Products are 16 `mul` + 16 `umulh` (squaring: 10 + 10 with doubling). The high half is
  folded twice by c. A final carry selects between y and y − c, which gives a fully
  reduced output without a separate comparison with p.
- `measurement/secp_sub.h` holds branchy subtractions that subtract c (that is, add p
  mod 2^256) only on a borrow, and a fused a − b − c.
- Standalone costs in the same toolchain, per call including about 7 instructions of
  call and loop overhead: mul 108, sqr 85, sub 22 on average. For comparison, the pinned
  libsecp256k1 5x52 multiplication is 179. Inside the kernel the routines are inlined.
- Correctness:
  - The routines matched Python big-integer arithmetic on about 194,000 reduced and
    199,000 arbitrary 256-bit operand pairs. These included constructed cases for the
    rare second-fold carry.
  - The replica check below compares every final point of 2^16 walks after 2^19 steps
    with Python big-integer arithmetic, on two input sets.

Inline assembly is ordinary scalar AArch64 code. Callgrind counts each executed
instruction exactly as for compiled C, which is how the organizer measured blst's assembly.

**Per-step measurement (H3).** `measurement/step.c` implements steps 3–6 for a batch of
2^16 walks with 2^16 adders.
- It was compiled with `-O3 -fno-tree-vectorize -mgeneral-regs-only` in the organizer's
  Docker image (GCC 14.2.0, Valgrind 3.24.0), and counted with Callgrind `Ir`.
- Counts are differenced between 8 and 4 batch-steps (2^18 steps).
- Every instruction inside the measured loop is charged.

| Input set | total Ir/step | 2-cycles in 2^19 steps |
| ---: | ---: | ---: |
| seed 1 | 780.992 | 0 |
| seed 2 | 781.008 | 0 |

This covers:
- 5 multiplications and 1 squaring;
- 5 subtractions, one of them fused;
- the fraction approximations and the k decision;
- key, parity and adder-image selection;
- the 2-cycle, anchor and distinguished-point filters, and any confirmations and
  escapes;
- the batch's prefix stores and loads;
- one Fermat inversion per 2^16 walks (libsecp256k1's addition chain, 255 squarings and
  15 multiplications, about 22,000 Ir, or 0.34 per step);
- loop control.

**795 instructions** are charged per step, so each step costs **0.795 units**.
`measurement/walk_ref.py check` replays all 2^16 walks in Python integers and matches the C
final points exactly on both input sets. Neither set entered a 2-cycle, so the cycle paths
were also checked on a build with 2^6 adders and 2^12 walks
(`evidence/measurement/cycles-a6.json`). Over 2^19 steps it had 1,337 2-cycles and 5 longer
cycles; 1/(6A) predicts about 1,365 2-cycles. The replica matched every final point and
both cycle counts.

Rare paths fit in the 14-instruction margin:
- **2-cycles.** At most 1/A per step × 3 unbatched operations × 2^15 = 1.5 instructions
  per step. The expected rate 1/(6A) gives 0.25.
- **Short cycles and exact-k fallbacks.** Below 2^−15 instructions per step.

A separately measured unbatched step (affine add with its own Fermat inversion, plus
adder selection) costs 22,913.1 and 22,914.0 Ir. Every unbatched operation below is
charged 2^15.

**Total.** S_max = K·(L + CAP) ≤ 2^126.749960 steps.

| Component | Bound (units) | log2 |
| --- | --- | ---: |
| walk steps | 0.795·S_max | 126.418986 |
| adder table | 2^16 × 2 scalar mults × 768 unbatched ops × 2^15 / C | 31.62 |
| walk starts and restarts | 2K × 2 scalar mults × 768 × 2^15 / C | 48.62 |
| distinguished points | DP_MAX = 2^67.75 entries × 2^20 instructions / C | 77.78 |
| collision replay | 2·(L + CAP) unbatched steps × (2^15 + 2^10) / C | 100.83 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 16.03 |
| input validation | 2^24 instructions | 14.03 |

Rounded up, **time_log2 = 126.418986**. Preprocessing is included. All processors' work
is summed, and elapsed time plays no part.

## 4. Memory and advice

Peak memory is at most the sum of:
- 2·DP_MAX·64 bytes for the distinguished-point table;
- K·448 bytes of walk state;
- 2^16·6·64 bytes of adder images;
- 16 MiB.

That is 2^74.74996 bytes. Nonuniform advice is 0 bytes.

## 5. Limitations

- H1 is the standard but unproved random-mapping heuristic. It is supported by the
  literature and toy-curve experiments only. The model bound is 0.4115, 0.0215 above
  0.39.
- Instruction counts come from one implementation and compiler. The sample is 2^18 steps
  on each of two input sets. The branchy subtractions and comparisons were measured, not
  bounded analytically; each subtraction's worst case is 6 instructions above its
  no-borrow path.
- The field routines were tested, not formally verified. The replica check covers 2^19
  production steps per input set.
- About 2^67.7 distinguished points and 2^126.4 work are needed. This is not a practical
  attack and does not weaken secp256k1.

## References

- Pollard (1978).
- van Oorschot–Wiener (J. Cryptology 1999).
- Teske (Math. Comp. 2001).
- Duursma–Gaudry–Morain (Asiacrypt 1999).
- Bos–Kleinjung–Lenstra (ANTS 2010).
- Bernstein–Lange–Schwabe (PKC 2011).
- Bernstein–Lange (ANTS 2012).
- Bos–Costello–Miele (PKC 2014).

## Reproduction

```sh
python3 -I submissions/secp256k1-ecdlp/rho-sigma6-asm/accounting.py
python3 -I submissions/secp256k1-ecdlp/rho-sigma6-asm/measurement/frac_check.py 1
bash submissions/secp256k1-ecdlp/rho-sigma6-asm/measurement/reproduce.sh 1   # Docker, linux/arm64
echo '{"bits":26,"trials":150,"seed":2601,"b":7}' | python3 -I submissions/secp256k1-ecdlp/rho-sigma6-asm/experiment.py
```
