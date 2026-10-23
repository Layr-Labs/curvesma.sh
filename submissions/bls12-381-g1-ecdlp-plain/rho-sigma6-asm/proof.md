# bls12-381-g1-ecdlp-plain: parallel Pollard rho on order-6 classes with a sigma walk, all arithmetic in instructions

Claim: for uniform x in [0,r) and T=[x]P0 in G1, the algorithm below returns x with
probability at least 0.39. It uses at most **2^126.312725** target-group-operation
equivalents in total (score 5052509/40000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^74.178504 bytes. Preprocessing is
2^45.619692 units, and no nonuniform advice is used. This is a classical cold attack on
one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on automorphism classes, and on the rate bound H2. The cost rests on
the instruction counts in H3, measured with the organizer's toolchain. This is a
carefully accounted generic attack. It is not a practical break.

**Cost route.** The walk does not use the 6M reference operation. Its five field
multiplications and one squaring per step are the package's own 6×64-bit Montgomery
multiplication in GCC inline AArch64 assembly (`measurement/bls_field.h`, the same
convention as blst's `mul_mont_384`). Every instruction of the step is charged at 1/C
with C = 2000 (H = 0 in COST_MODEL's H + R/C). COST_MODEL.md: "Faster arithmetic and
lookup algorithms are welcome; justify their actual work". A step measures 2,177
instructions, so it is charged 1.0975 units rather than 1 + R/C = 1.14. The same walk with
the reference operation is the companion package rho-sigma6 (126.367538).

What is new relative to the usual order-6 walk is in section 1.2. The walk never moves
a point to its canonical representative, so it never pays the multiplication by beta.
It stores the actual point W and adds sigma(R_i), where sigma is the automorphism that
maps the canonical point to W. The class sequence is unchanged, so the standard
analysis applies unchanged. Which rotation is canonical is decided from 64-bit
fixed-point approximations with an exact fallback, and Lemma 1 bounds their error.

Files:
- `accounting.py` recomputes every number with exact integers (output in
  `evidence/accounting.json`).
- `measurement/` holds the per-step instruction-count harness, the field code
  `bls_field.h`, the Python replica `walk_ref.py`, and the approximation check
  `frac_check.py`.
- `evidence/measurement/` holds the counts on two input sets.
- `experiment.py` is the declared toy-curve experiment for the identical walk;
  `evidence/small-curves/` holds larger offline runs of the same file.

## 1. Algorithm

### 1.1 Notation and parameters

p, r, P0 and the curve y^2 = x^3 + 4 come from the protected profile. r is prime, and
p ≡ r ≡ 1 (mod 3).

The endomorphism uses beta = 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe,
a primitive cube root of unity mod p, and lambda = 0x73eda753299d7d483339d80809a1d804a7780001fffcb7fcfffffffe00000001
= −z^2 mod r, where z = −0xd201000000010000. Then [lambda](x,y) = (beta·x, y) on G1. This was
checked on P0 (`measurement/walk_ref.py` uses it; [lambda]P0 = (beta·x0, y0)).

The automorphism group G6 = {±1, ±lambda, ±lambda^2} acts on G1. Each nonzero point lies
in a class of exactly 6 points, so there are (r−1)/6 nonzero classes. Write N = (r−1)/6 + 1.

M(v) = v·2^384 mod p is blst's Montgomery representative, a fixed bijection of Fp onto
[0,p). Every comparison and bit test below is on M-values, read as integers.

Parameters:
- A = 2^16 adders.
- K = 2^30 walks: 2^16 processors × 2^14-walk batches.
- theta = 2^−60.
- L = 89663270682989036334838442605 ≈ 2^96.179.
- CAP = 64/theta.
- DP_MAX = 2·theta·S_max + 2K, where S_max = K·(L + CAP).

### 1.2 Classes, keys and the sigma walk

For an affine point W = (x, y) ≠ O:

- **Class key.** key(W) = M(y) if M(y) is even, else p − M(y), so key(W) = M(±y) and is even.
  - key is constant on the class, since G6 changes y only by sign.
  - key determines the class. y^2 = x^3 + 4 fixes x^3, hence x up to a factor beta^j, and
    ±y fixes the sign.
- **Canonical point.** canon(W) = (beta^k·x, (−1)^s·y), where k ∈ {0,1,2} minimises
  M(beta^k·x) and s = M(y) mod 2.
  - The three values are distinct for x ≠ 0, because beta^j·x differ and M is a bijection.
  - x = 0 does not occur in G1: (0, ±2) have order 3.
  - p is odd, so M(−y) = p − M(y) has the other parity. canon(W) therefore has an even
    M(y) and is a fixed point of the class.
- **sigma.** sigma_W = (x → beta^(−k)·x, y → (−1)^s·y), so that W = sigma_W(canon(W)).
  sigma_W is a group automorphism, equal to [(−1)^s·lambda^(−k)].
- **Walk step.** i(W) = (key(W) >> 1) mod 2^16, and

      f(W) = W + sigma_W(R_i(W)).

  Because sigma_W is an automorphism, f(W) = sigma_W(canon(W) + R_i). So the class of f(W)
  is the class of canon(W) + R_i. That is exactly the standard class walk
  [W] → [canon(W) + R_h([W])] with the class function h = i. Every statement about the
  class sequence below is a statement about that standard walk.
- **Images of each adder.** The six images sigma(R_i) = (beta^(−k)·X_i, ±Y_i) are tabulated
  per adder, so a step reads sigma_W(R_i) from T6[i][k][s] and no point is ever rotated
  or negated.
- **Distinguished point.** W is distinguished when bits 17..76 of key(W) are zero, which
  is a class property.

**Deciding k without beta·x (Lemma 1).** Let a = M(x) = a_0 + a_1·2^64 + … + a_5·2^320.
Let C_0 = ⌊2^512/p⌋ and C_1 = ⌊2^512·beta/p⌋, with 64-bit limbs C_{j,m}. Define

    g_j = ( Σ_{i+m=7} lo(a_i·C_{j,m}) + Σ_{i+m=6} hi(a_i·C_{j,m}) ) mod 2^64,   j = 0, 1,

where lo and hi are the low and high words of the 128-bit product. For C_0 only m = 1, 2
are nonzero, so g_0 takes 3 multiplications and g_1 takes 12.

Let f_j = 2^64·M(beta^j·x)/p, a real number in (0, 2^64). Note M(beta·x) = beta·a mod p.

**Lemma 1.** If g_j ≤ 2^64 − 14, then f_j ∈ [g_j, g_j + 13) for j = 0, 1.

*Proof.* f_j = 2^64·frac(a·c_j) with c_0 = 1/p and c_1 = beta/p. Here a·c_j differs from
a·C_j/2^512 by less than a/2^512 < 2^−131. Limb 7 of a·C_j mod 2^512 equals g_j plus the
carry from the dropped terms. Those terms are the low words with i+m = 6 and all terms
with i+m ≤ 5. Their total is below 6·2^448 + 6·2^448·(1 + 2^−64 + …) < 13·2^448, so the
carry is at most 12. Products with i+m ≥ 8, and high words with i+m = 7, are multiples of
2^512 and vanish mod 1. The wrap-around case is excluded by g_j ≤ 2^64 − 14. ∎

M(x) + M(beta·x) + M(beta^2·x) ≡ 0 mod p, with each term in (0, p), so the sum is p or 2p
and f_2 = 2^64·m − f_0 − f_1 for m ∈ {1, 2}. With g_2 = (−g_0 − g_1 − 26) mod 2^64, the
real f_2 lies in (g_2, g_2 + 26] whenever g_2 ≤ 2^64 − 28. The kernel accepts the index of
the least g_j as k only if all three g_j ≤ 2^64 − 28 and the least is more than 26 below
both others. Then f_k ≤ g_k + 26 < g_l ≤ f_l, so k is exact. Otherwise it computes beta·x
and compares exactly.

The fallback probability per step is at most about 3·2·27/2^64 + 3·28/2^64 < 2^−55. It
never occurred in 2^19 measured steps. `measurement/frac_check.py` checks the floor gap
⌊f_j⌋ − g_j ∈ [0,12] on 200,000 random and edge-case inputs per seed; the largest gap seen
is 6.

### 1.3 Procedure

1. **Validate T.** Check it is on the curve and [r]T = O. If T = O, output 0.
2. **Adders.** For each adder i < 2^16:
   - draw c_i, d_i uniform in [0,r);
   - set R_i = [c_i]P0 + [d_i]T;
   - tabulate its six images, using one multiplication by beta and one negation.

   For each walk w, draw a_w, b_w and set W_w = [a_w]P0 + [b_w]T. All coins are stored.
   Coefficients are not tracked during the walk.
3. **Step.** W ← f(W), as one batched incomplete affine addition: Montgomery's
   simultaneous inversion over the 2^14 walks of a batch (3M), then lambda, lambda^2 and
   y3 (3M), with one inversion per batch.
   - If the batch product is zero, then W.x = sigma(R).x for some walk, with probability
     at most 2^−238 per step. That walk then uses the unbatched complete formula.
   - Each walk also records key(W) mod 2^64 and selects sigma_W(R_i) for its next step.
4. **Fruitless 2-cycles.** Each walk keeps the low 64 bits of its previous key.
   - If the new state's key matches, the walk evaluates f(new) unbatched.
   - If key(f(new)) = key(old) exactly, then {old, new} is a 2-cycle of the class walk.
     The walk replaces new by 2·canon(M), where M is the element with the smaller key.
5. **Short cycles (3..64).** Every 64 steps each walk records its key as an anchor.
   - On a later match, the walk iterates f unbatched up to 64 times.
   - If its class returns, the walk escapes from 2·canon(M), with M the cycle element of
     least key. Otherwise nothing changes.
6. **Distinguished points.** The walk sends (key, w, step index) to a hash table keyed by
   key.
   - A walk runs L steps and then continues until its next distinguished point.
   - A walk that exceeds CAP steps without one is abandoned and restarted from fresh coins,
     with at most K restarts in total.
   - The run fails if more than DP_MAX points arrive or the restarts run out.
7. **Resolution.** On the first key stored twice, both walks are replayed from their
   coins to the recorded step indices.
   - The replay tracks (a, b) mod r through sigma: an image sigma(R_i) has coefficients
     multiplied by (−1)^s·lambda^(−k). Escapes double canon(M), whose coefficients are
     multiplied by (−1)^s·lambda^k.
   - It obtains the canonical points canon(W) = [a]P0 + [b]T of both walks. These are
     equal, since the keys are equal.
   - If b1 ≠ b2, it outputs x = (a2 − a1)/(b1 − b2) mod r after checking [x]P0 = T.
     Otherwise it continues; the bound counts this case as a failure.

**Merged trails stay merged.** Outside fruitless cycles the next class depends only on
the current class. An escape fires only after an exact check of a genuine cycle of the
class walk. The exit point 2·canon(min) depends only on the cycle's set of classes, not
on the entry point or on when the walk noticed. The 64-bit filters only decide when an
exact check runs. So two trails that reach the same class follow the same classes
afterwards and reach the same next distinguished class.

## 2. Success probability

Let T_steps = K·L be the main-phase steps, each walk's first L steps. Then
tau = T_steps^2/(2N) = 0.530300.

**H1 (heuristic, disclosed).** The class map [W] → [canon(W) + R_h([W])], with uniform
random adders, behaves like a random mapping on the N classes, apart from the known
r-adding penalty. An r-adding walk loses about a factor (1 − 1/(m·A)) in the birthday
exponent, with A the number of adders and m = 6 the class size (Teske 2001;
Bernstein–Lange 2012; Bos–Kleinjung–Lenstra 2010 eq. (2) for m = 2).

Under H1, the probability that the first T_steps steps contain no useful class collision
is at most exp(−tau·(1 − eps)). Here eps = 4/A + 2·(3/A) + 2·(192/A^2) = 0.00015268. It
collects three margins:

- **Adder penalty, 4/A.** This is 24 times the 1/(6A) estimate.
- **Fruitless 2-cycles, 2·(3/A).** H2 bounds the entry rate by 1/A per step, and each
  event yields at most 3 non-fresh states.
- **Short cycles, 2·(192/A^2).** H2 bounds the rate by 1/A^2, and detection wastes at most
  192 states per event.

**H2 (rates).**
- 2-cycles are entered at a rate of at most 1/A per step. The literature gives 1/(6A)
  for order-6 classes (Bos–Costello–Miele 2014).
- Cycles of length 3..64 are entered at a rate of at most 1/A^2.
- Cycles longer than 64 with no distinguished point arise at a rate of at most A^−32.
  They cost at most CAP steps each through abandonment, which is charged.

A detected collision needs both trails to reach the shared next distinguished class.
That fails only if the class is more than CAP steps away: probability e^−64 per segment,
about 2^−25 over all segments. The remaining failure terms are:

- degenerate b1 = b2: at most 2^−250;
- DP_MAX or restart overflow: under 2^−1000 by Chernoff;
- exceptional additions: under 2^−110.

Together with the abandonment term these are below delta = 2^−20. Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015268)) − 2^−20 = 0.411523 ≥ 0.39.**

The budget is 1.0299·√N steps, against 0.9943·√N for success exactly 0.39. The 0.0215
of slack in the model bound covers the sampling error of the toy evidence below. The
bound is exact arithmetic on H1, not an empirical estimate.

**Evidence for H1 and H2.** `experiment.py` runs this exact sigma walk on toy curves
y^2 = x^3 + 4 over p ≡ 1 mod 3, in a prime-order subgroup of small cofactor (the BLS12-381
G1 shape). The toy uses the same pieces as production:
- Montgomery-representative keys, M(v) = v·2^64 mod p;
- the same canonical rule, adder images and index;
- both cycle rules and escapes;
- distinguished points.

It counts every step, including fruitless, confirming and escape steps, up to the first
useful class collision. It then solves end to end through the distinguished-point table
and checks [x]P = Q.

<!-- toy-evidence -->
Offline runs (`evidence/small-curves/`, pooled in `summary.json`): 5,200 trials on 6 curves
with N = 2^28.8, 2^27.5, 2^35.5, 2^35.5, 2^35.5, 2^39.7 classes, adder tables of 2^11 and 2^16.
- Every trial solved and verified [x]P = Q.
- Success by T = sqrt(2 tau N) at the claimed tau = 0.5303: pooled 0.4094, Wilson 95%
  interval [0.396, 0.423]. The model value is 0.412.
- At tau = 0.4943 (model 0.390) the pooled success is 0.3881, interval [0.375, 0.401].
- Mean steps to the first useful collision: 1.248 sqrt(N), against sqrt(pi/2) = 1.253.
- 2-cycle entry rate per step:
  - 8.09e-05, 8.36e-05, 8.35e-05, 8.14e-05, 8.15e-05 with A = 2^11 (1/(6A) = 8.14e-05);
  - 2.50e-06 with A = 2^16 (1/(6A) = 2.54e-06).
- Longer fruitless cycles: 5 in total. Abandoned walks: 10.
<!-- /toy-evidence -->

Limits: the toy curves have N up to about 2^40, against 2^252 in production. H1 is an
extrapolation in N, supported by the size-free random-mapping argument and the cited
literature, including the 55- and 60-bit negation-map challenge runs in
Bernstein–Lange–Schwabe 2011. It is not proved for this walk.

## 3. Computation

**No reference operation.** Each step is charged its full measured instruction count
divided by C = 2000, as COST_MODEL.md's H + R/C with H = 0. No 6M reference operation is
claimed, so nothing is charged twice.

**Field arithmetic.** `measurement/bls_field.h` implements Montgomery multiplication and
squaring with R = 2^384 as GCC inline AArch64 assembly in general-purpose registers.
- Multiplication interleaves six word-by-word Montgomery reduction rows with the
  product. The cancelled low limb of each row is replaced by a carry test, and one
  final conditional subtraction (`csel`) makes the output fully reduced.
- Squaring computes the 21 distinct products with the doubling fused into the diagonal,
  then a 12-word Montgomery reduction.
- Standalone costs in the same toolchain, per call including 6 instructions of loop and
  call overhead: mul 334 and sqr 279, against blst's 383 and 325 in the same harness.
  Inlined in the kernel, mul is about 325.
- Correctness:
  - The routines matched Python big-integer arithmetic on 380,000 cases, with seven
    outputs each, including aliasing and products landing on 0, 1, p − 1 and p − 2.
  - The replica check below compares every final point of 2^14 walks after 2^18 steps
    on two input sets.
- blst is still used for the one inversion per batch (`blst_fp_inverse`) and for input
  conversion. Those instructions are counted like any others.

Inline assembly is ordinary scalar AArch64 code. Callgrind counts each executed
instruction exactly as for compiled C, which is how the organizer measured blst's assembly.

**Per-step measurement (H3).** `measurement/step.c` implements steps 3–6 for a batch of
2^14 walks with 2^16 adders.
- It links the organizer-pinned blst (f7a88067…), built exactly as in the organizer's
  harness (server.c with `-O3 -fno-builtin -fPIC -fno-tree-vectorize
  -mgeneral-regs-only`, plus the upstream AArch64 assembly).
- The kernel itself is compiled with `-O3 -fno-tree-vectorize -mgeneral-regs-only` in the
  organizer's Docker image (GCC 14.2.0, Valgrind 3.24.0, Callgrind `Ir`).
- Counts are differenced between 16 and 8 batch-steps (2^18 steps).
- Every instruction inside the measured loop is charged.

| Input set | total Ir/step | 2-cycles in 2^18 steps |
| ---: | ---: | ---: |
| seed 1 | 2177.255 | 0 |
| seed 2 | 2177.228 | 1 |

This covers:
- 5 multiplications and 1 squaring;
- all field subtractions (5 per step, one of them fused a−b−c), 6-limb inline-assembly
  sequences that add p back only on a borrow;
- the two fraction approximations and the k decision;
- key, parity and adder-image selection;
- the 2-cycle, anchor and distinguished-point filters, and any 2-cycle confirmations and
  escapes that occurred;
- the batch's prefix stores and loads and the one inversion per 2^14 walks (26,246 Ir,
  or 1.6 per step);
- loop control.

**2,195 instructions** are charged per step, so each step costs **1.0975 units**.
`measurement/walk_ref.py check` replays all 2^14 walks in Python integers and matches the C
final points exactly on both input sets, so the measured binary computes the specified
walk.

The rare paths are charged within the 17.7-instruction margin:
- **2-cycles.** At most 1/A per step, each with at most 3 unbatched operations of at most
  2^15 instructions: at most 1.5 instructions per step.
- **Short cycles.** Below 2^−21 instructions per step.
- **Exact k fallbacks.** Below 2^−45.

A separately measured unbatched step (affine add with its own `blst_fp_inverse`, plus
image selection) costs 28,553.8 Ir on both sets. Every unbatched operation below is
charged 2^15.

**Total.** S_max = K·(L + CAP) ≤ 2^126.178504 steps.

| Component | Bound (units) | log2 |
| --- | --- | ---: |
| walk steps | 1.0975·S_max | 126.312725 |
| adder table | 2^16 adders × 2 scalar mults × 768 unbatched ops × 2^15 / C (plus images) | 30.62 |
| walk starts and restarts | 2K × 2 scalar mults × 768 × 2^15 / C | 45.62 |
| distinguished points | DP_MAX = 2^67.179 entries × 2^20 instructions / C | 76.21 |
| collision replay | 2·(L + CAP) unbatched steps × (2^15 + 2^10) / C | 101.26 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 15.03 |
| input validation ([r]T = O check) | 2 scalar mults + 2^24 instructions | 15.03 |

Rounded up, **time_log2 = 126.312725**. Preprocessing (adders,
starts and restarts) is 2^45.619692 units and is included above. The 2^20 instructions
charged per distinguished point cover hashing, insertion and lookup, transfer to the
table host and bookkeeping. All processors' work is summed, and elapsed time plays no
part.

## 4. Memory and advice

Peak memory is at most the sum of:
- 2·DP_MAX·64 bytes for the distinguished-point table at a 2× load factor;
- K·512 bytes of walk state (two point buffers, batch scratch, keys);
- 2^16·6·96 bytes of adder images;
- 16 MiB of code, stacks and scratch.

That is 2^74.178504 bytes. No advice or precomputed table exists before T is given, so
nonuniform advice is 0 bytes.

## 5. Limitations

- H1 is the standard but unproved random-mapping heuristic for additive class walks. It
  is supported here by the literature and by toy-curve experiments only. The model bound
  is 0.4115, 0.0215 above 0.39.
- Instruction counts come from one implementation, compiler and convention. The sample
  is 2^18 steps on each of two input sets.
- The field multiplication is the package's own assembly, tested against Python big
  integers but not formally verified.
- The variable-time paths (branchy subtractions, comparisons, the exact fallback) were
  measured, not bounded analytically. Lemma 1 makes the fallback's rarity analytic.
- About 2^67 distinguished points and 2^126.3 work are needed. This is not a practical
  attack and does not weaken BLS12-381's security assumptions.
- The pairing (MOV) route to GT is not used.

## References

- Pollard (1978).
- van Oorschot–Wiener, *Parallel collision search* (J. Cryptology 1999).
- Teske, *On random walks for Pollard's rho method* (Math. Comp. 2001).
- Duursma–Gaudry–Morain (Asiacrypt 1999).
- Bos–Kleinjung–Lenstra, *On the use of the negation map in Pollard rho* (ANTS 2010).
- Bernstein–Lange–Schwabe, *On the correct use of the negation map in the Pollard rho
  method* (PKC 2011).
- Bernstein–Lange, *Two grumpy giants and a baby* (ANTS 2012).
- Bos–Costello–Miele, *Elliptic and hyperelliptic curves: a practical security analysis*
  (PKC 2014).

The construction and accounting above are complete and do not depend on fetching these.

## Reproduction

```sh
python3 -I submissions/bls12-381-g1-ecdlp-plain/rho-sigma6-asm/accounting.py
python3 -I submissions/bls12-381-g1-ecdlp-plain/rho-sigma6-asm/measurement/frac_check.py 1
bash submissions/bls12-381-g1-ecdlp-plain/rho-sigma6-asm/measurement/reproduce.sh 1   # Docker, linux/arm64
echo '{"bits":26,"trials":150,"seed":2601}' | python3 -I submissions/bls12-381-g1-ecdlp-plain/rho-sigma6-asm/experiment.py
# offline evidence runs used larger inputs, e.g. {"bits":40,"trials":700,"seed":4001}
```
