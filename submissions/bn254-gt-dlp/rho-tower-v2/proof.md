# BN254 GT discrete logarithm: the rho-tower baseline with a cheaper measured step

Track `bn254-gt-dlp`. Input: `z = g^x` in the prime-order group GT of BN254,
with `g` the profile generator `e(P0,Q0)` and `x` uniform in `0..r-1`. Output:
the scalar `x`. One unit is `C = 20,000` scalar AArch64 instruction equivalents.
This is an analytic exploratory claim: no production-size attack was run, and
the success assessment uses the standard random-map heuristic, stated and
bounded below.

This package is the `rho-tower-v1` algorithm with a **reworked arithmetic
kernel**; it reduces the measured per-step price from 47,239.44 to 31,429.48
instructions, a 33.5% cut, and moves the declared score from 130.851718 (our
resubmitted v1) to 130.131912. The algorithm, ledger structure and heuristic are
unchanged, so the gain is purely measured work, not a new attack.

This package also carries the two fixes required after our sibling submission
`bf1247b4` was rejected (it shares the uncapped sampler): the scalar sampler now
has a hard cap of K = 64 attempts per scalar with its failure probability
charged against the success bound, and the per-step price is the measured
worst-case path over every production table entry and data-dependent branch,
not an average with a 10% pad. Both are changes to the declared algorithm and
are described in sections 2 and 4.

## Claim

For the protected BN254 target, the single-trial algorithm of this document
returns `x` with `g^x = z`, with success probability at least **0.85** over
uniform challenges and fresh algorithm coins, and every run performs at most

    W <= 2^130.131912 < 2^130.131913   units,

where one unit is 20,000 scalar AArch64 instructions. The score is
`130131913/1000000`. This is a per-run deterministic work cap: every loop has a
data-independent trip count (`LOOP_BOUNDS.txt`), the sampler has a hard cap of
K = 64 attempts per scalar, the walk stops within `N = 2*ceil(sqrt(r))`
iterations, and the table build, recovery, verification and output work are
fixed. Success is at least the protected 0.39.

The claim is a **baseline**, not a best-known attack. It does not use NFS or
index calculus. Its point is that the arithmetic kernel is measured, not
assumed, and that the measurement is reproducible from the shipped code.

## 1. Representation and conversion

The profile defines

    Fp12 = Fp[w]/(w^12 - 18 w^6 + 82),   u = w^6 - 9,   Fp2 = Fp[u]/(u^2+1).

The attack uses the isomorphic tower

    Fp2 = Fp[u]/(u^2+1),  Fp6 = Fp2[v]/(v^3 - (9+u)),  Fp12 = Fp6[w_t]/(w_t^2 - v).

`w_t^6 = v^3 = 9+u` and `w_t^12 = (9+u)^2 = 18 w_t^6 - 82`, so `w_t` satisfies
the profile polynomial; `w -> w_t` is an Fp-isomorphism. In the tower basis the
profile-basis coefficients `t_0..t_11` map to

    B=t6, A=t0+9B;  H=t7, G=t1+9H;  D=t8, C=t2+9D;
    J=t9, I=t3+9J;  F=t10, E=t4+9F; L=t11, K=t5+9L,

with tower coordinates `(A+Bu) + (C+Du)v + (E+Fu)v^2 +
((G+Hu) + (I+Ju)v + (K+Lu)v^2) w_t`. This is a fixed public linear map, cost
twelve Fp multiplications and a bounded number of additions, included in the
setup row. The shipped `selftest` prints `g_is_one=0 g^r_is_one=1`, a direct
check of the conversion and of the tower arithmetic together.

## 2. Algorithm

Public input: the profile `p, r, g` and `z`.

1. Decode `z` canonically, check `z != 0` and `z^r = 1`; if `z = 1` output
   `x = 0`. Convert `g` and `z` to the tower (section 1).
2. Draw one 256-bit seed `sigma` from the OS CSPRNG and one fresh 64-bit salt.
   The word stream is `w_i = SHA-256(sigma || LE64(i))` truncated to 256 bits.
   Every uniform scalar in `0..r-1` is produced by the capped rejection sampler

       S:  for k = 1..K,  K = 64:
               v <- next word          (uniform in 0..2^256-1)
               if v < r: return v
           abort: the run stops and reports failure

   Draw the 64 table scalars `c_j, d_j` and then `a_0, b_0` with independent
   calls to `S`, in that fixed order (at most 66*K = 4,224 words).
3. Build a 32-entry table with

       M_j = g^{c_j} z^{d_j},   j = 0..31,

   each entry stored with its two scalar labels. Two fixed-base/fixed-exponent
   exponentiations per entry; see the setup row.
4. Choose `a_0, b_0` uniform in `0..r-1` and set `Y_0 = g^{a_0} z^{b_0}`.
5. The walk state is `(Y, a, b)`. One step is

       j      = H(salt || enc(Y)) mod 32
       Y      = Y * M_j
       a      = a + c_j mod r
       b      = b + d_j mod r

   where `enc(Y)` is the 384-byte concatenation of the twelve Montgomery
   coefficients and `H` is the 64-bit two-accumulator mixer in `rho.c`. The
   invariant `Y = g^a z^b` holds after every step.
6. Floyd: initialise tortoise and hare at `(Y_0,a_0,b_0)`. For
   `k = 1..N`, advance the tortoise one step and the hare two steps, then test
   the two group elements for equality (twelve Fp words).
7. At the first equality let the tortoise label be `(a,b)` and the hare label
   be `(A,B)`. If `B - b = 0 mod r`, report failure. Otherwise

       x = (a - A) * (B - b)^{-1} mod r.

8. Verify `g^x = z` in the tower; a mismatch is failure, never an accepted
   answer. Return `x` (32-byte scalar).

There is exactly one trial, one table, one salt and no restart. `N` is fixed
before the run. Every failure path spends at most the same bounded work.

**Capped sampling and its failure accounting.** Conditional on acceptance the
sampler returns a scalar exactly uniform on `0..r-1` under the PRF assumption
stated in section 3, so the walk analysis of
section 6 is unchanged. One call of `S` exhausts its budget with probability at
most `delta^K`, `delta = 1 - r/2^256 = 0.810969445...`, `delta^64 <= 2^-19.35`;
the union bound over the 66 calls is at most `66*delta^64 <= 2^-13.30`
(recomputed from `r` in `experiment.py`). Success is therefore at least
`0.8559 - 2^-13.30 = 0.855801 > 0.85`, and an aborted run stops after at most
the sampler's capped cost, so the per-run cap holds on the failure path too.

## 3. Correctness

The invariant is `Y = g^a z^b` initially and after a step, because
`Y M_j = g^a z^b * g^{c_j} z^{d_j} = g^{a+c_j} z^{b+d_j}`. If
`g^a z^b = g^A z^B` and `z = g^x`, then `a + xb = A + xB mod r`, hence
`a - A = x(B - b) mod r`. Since `r` is prime, a nonzero `B - b` has a unique
inverse and the recovered `x` is the unique discrete logarithm. The final
`g^x = z` check rejects any wrong scalar.

Floyd's equality is a genuine collision of two walk states with different step
counts. If the deterministic walk has tail `mu` and cycle `lambda`, an equality
is observed at some `k` with `mu <= k <= mu + lambda`, so `N >= mu + lambda`
suffices for detection. The scalar labels are carried along, not derived from
the hash, so the collision yields the linear relation above.

Degenerate collisions (`B - b = 0`) are negligible. Conditional on the group
walk and on the salt, the jump logarithms `s_j = c_j + x d_j mod r` are fixed;
the individual `d_j` remain independent uniform because each `c_j` is uniquely
determined by `s_j`. The coefficient difference at a collision is
`sum_j n_j d_j` with `sum_j n_j = k`, `0 < k <= N < r`, and at least one
`n_j != 0 mod r`; a nonzero integer combination of independent uniform field
elements is exactly uniform, so the degeneracy probability is `1/r`
conditioned on any such walk.

**One idealisation is explicit here.** The word stream is a deterministic
SHA-256 expansion of one 256-bit seed, so the sampler's output is exactly
uniform only under the standard pseudorandomness assumption that the expanded
words are computationally indistinguishable from independent uniform 256-bit
strings (equivalently: SHA-256 in counter mode is a PRF under a secret seed).
Conditional on that assumption the rejection sampler is exactly uniform, the
`1/r` degeneracy argument is exact, and the `delta^K` abort bound is exact;
without it, both are heuristic. The seed is fresh from the OS CSPRNG on every
run and is never reused, so the assumption is the usual PRF one and not a claim
about any fixed public value.

## 4. Arithmetic and the measured step cost

All arithmetic is scalar 64-bit. `Fp` elements are four little-endian words in
Montgomery form; `fp_mul` is a 4-limb CIOS Montgomery multiplication with a
six-word accumulator. Three changes from `rho-tower-v1`:

1. **`Fp6` Karatsuba.** `Fp6 = Fp2[v]/(v^3 - xi)` multiplication now uses six
   `Fp2` products instead of nine:

       v0=x0y0, v1=x1y1, v2=x2y2,
       s01=(x0+x1)(y0+y1)-v0-v1,  s02=(x0+x2)(y0+y2)-v0-v2,
       s12=(x1+x2)(y1+y2)-v1-v2,
       c0 = v0 + xi*s12,  c1 = s01 + xi*v2,  c2 = s02 + v1.

   The three saved `Fp2` products cost 9 fewer `fp_mul` per `Fp12` product
   (three `Fp6` products per `Fp12` product) and add six `Fp2` subtractions.
2. **`Fp2` multiplication by `xi = 9+u`.** The old routine computed `2a`, `3a`,
   `6a` and then recomputed `8a` from scratch, leaving dead arithmetic. The new
   routine computes `8a+a` and `8b+b` directly and drops the dead adds.
3. **Branchless inline `Fp` addition and subtraction.** `fp_add` and `fp_sub`
   are `static inline` in `fp12.h`; the canonical reduction after addition is a
   mask select instead of a comparison loop plus `memcpy`, and subtraction adds
   the modulus under a borrow mask. The tower makes far more calls to these
   than to `fp_mul`, so the per-call overhead mattered.

`fp12_mul` remains three `Fp6` products plus two `v` multiplications
(Karatsuba over `w_t^2 = v`). The binary `rho` was compiled with

    gcc 13.3.0 -O3 -fno-tree-vectorize -mgeneral-regs-only

on the same machine class as the cost model (AArch64, 64-bit general-purpose
registers), and measured with Valgrind 3.22.0 Callgrind (`Ir`).

**The bound is a static composition, not a sampled average.** The review of
submission bf1247b4 rejected the previous sweep on exactly this point, and the
criticism was right: differencing two runs on one trajectory measures the mean
cost of its later steps, not the worst cost of a step. The bound below is
composed instead, and it holds for every step of every run:

* **Call counts are fixed by construction.** One step executes exactly one
  `fp12_mul`, one mixer, two `scalar_add_r` and one table load; the tower
  expansion of `fp12_mul` is fixed by the source (3 `fp6_mul`, 3 `fp2_mul` per
  `fp6_mul`, 3 `fp_mul` per `fp2_mul`), so per step
  81 `fp_mul`, 0 `fp_add`, 0 `fp_sub`, plus the
  fixed mixer/loop code. Nothing about the data changes these counts.
* **Every loop inside those functions has a data-independent trip count**
  (`fp_mul` 4x4 limbs, the mixer 24 rounds over 48 words, `scalar_add_r` 4
  limbs), so a call's instruction count is its straight-line code plus the
  arms of its data-dependent decisions. The arms were measured one by one with
  `percall.c` (`evidence/percall-arms.txt`): `fp_add` 105.0 vs 89.0
  instructions, `fp_sub` 92.0 vs 60.0, `fp_mul` 282.0 vs 282.0 (the final
  reduction comparison is charged 8 extra instructions per call anyway).
  In this kernel `fp_add` and `fp_sub` are `static inline` and branchless (mask
  selects); `percall.c` measures the same cost for both operand extremes, so
  their arm delta is 0 and only `fp_mul` contributes.

* **The step was measured on the algorithm's own table, not a synthetic one.**
  The review also objected that the harness built `M_j = g^(j+1)`. That mode is
  gone: `./rho mktable 20261112 table.bin` runs the sampler of section 2
  (one 256-bit seed, SHA-256 counter-mode words, capped rejection, canonical
  scalars) and builds `M_j = g^c_j z^d_j` for `z = g^x`, and every measurement
  below runs `./rho stepwc K j table.bin` on that table. The seed draws 67
  scalars in 359 words; the table is written with both labels.

With `c_f` the fixed call count and `[min_f, max_f]` the measured per-call
range, the step obeys

    cost_per_step <= straight_line + sum_f c_f * max_f
                  = mean_per_step + sum_f c_f * (max_f - min_f)

because the straight-line part is exactly `mean_per_step - sum_f c_f * min_f`
for every step (the measured mean is an average of quantities that are all
`>= sum_f c_f * min_f + straight_line`). The measured mean on the production
table, over the 33 trajectories `j = 0..31` plus the unpinned production mixer,
is **31,589.48 Ir/step** (`evidence/step-measurement.txt`; per-trajectory spread
under 25 instructions). The composition adds

    81 x 8  +  0 x 16  +  0 x 32  =  648 instructions,

so

    cost_per_step <= 31,589.48 + 648 = 32,237.48    ->    C_STEP = 32,238 instructions.

The charge is the composed bound rounded up to the next instruction; there is
no percentage padding, and no part of it is an average over unmeasured steps.
The measured loop includes the hash, the two scalar additions, the table load
and the multiplication; it excludes the Floyd comparison and loop control,
charged separately at 4,096 instructions per iteration, an analytic over-charge
(48 word compares at ~3 instructions each plus a counter and a branch).

## 5. Complete work ledger

`r` is the profile subgroup order, `N = 2*ceil(sqrt(r)) =
295893513763578638011461384341992519220`.

| component | work | log2 units (C=20,000) |
|---|---:|---:|
| 3N walk steps at 32,238 instructions | 3N * 32,238 / 20,000 | 130.072070 |
| N Floyd iterations, 4,096 instructions each | N * 4,096 / 20,000 | 125.510633 |
| setup: 64 exponentiations, table, conversion | <= 2^36 instructions | 21.712992 |
| capped sampler: 66 * 64 draws at <= 2,048 instructions | <= 8,650,752 | 8.756682 |
| verification `g^x = z`, canonical output | <= 2^32 instructions | 17.712288 |
| output and final scalar inversion margin | <= 2^20 instructions | 5.712288 |
| **total W** | | **130.131912** |

The sampler row is charged at its cap (66 scalars x 64 attempts x 2,048
instructions = 8,650,752 < 2^24), not at an observed attempt count, so it is
input-independent.

The setup row: 64 exponentiations with 254-bit exponents, each at most 508
`Fp12` multiplications, is 32,512 multiplications; even at the loose analytic
bound of 2^20 instructions per multiplication used by the earlier published
BN254 baseline, that is below 2^36. We charge 2^36, about 30x the measured cost
of that many multiplications. The verification row is two exponentiations plus
decoding, charged at 2^32, over 60x the same estimate. The output row covers
the single scalar inversion and the 32-byte encoding with an additional margin.

No row is amortised over future targets, no oracle is free, and no BN254
operation is priced at one unit. Memory: 32 table records of 448 bytes, two
walk states, the hash scratch and code, below 2^20 bytes; the salt and
increments are fresh per run and charged in setup.

## 6. Success probability: heuristic, support and limits

Model the walk as a random map on the `r` points of GT, with the salt and the
64 increments as the fresh coins. For `t` distinct transitions the no-repeat
probability is `prod_{i=1}^{t-1}(1-i/r) <= exp(-t(t-1)/(2r))`. At
`t = N = 2 ceil(sqrt(r))` this is at most `exp(-2) = 0.1353`, so a repeat occurs
with probability at least `0.8647` even before the 32-adding correction. The
conventional r-adding correction multiplies the effective collision rate by
`1 - 1/32`; applying it at `t = 2 sqrt(r)` gives

    1 - exp(-(31/32)*2) = 0.8559.

We claim the smaller, conservative **0.85**, with the capped sampler's abort
probability `2^-13.30` subtracted from the collision model in section 2, so the
claimed floor is `0.8559 - 2^-13.30 = 0.85580`. This leaves a large margin over
the protected 0.39.

**Measured success of the production walk.** `experiment.py` runs the shipped
walk: the shipped 384-byte encoding, the mixer transliterated from `rho.c`
(cross-checked index for index against the binary's `hashdump` mode,
`evidence/mixer-crosscheck.txt`) and the shipped tower multiplication, inside
the order-`Q` subgroup of the *real* BN254 Fp12, `Q = 43,913` a prime factor
of `p^12 - 1`. The subgroup order is the only difference from production.

With seed 20261112 and 64 trials the measured result is **60 of 64 runs
recovered `x` exactly**, 4 runs had no collision inside `N =
2*ceil(sqrt(Q)) = 420`, and no run produced a wrong scalar or a degenerate
`B - b = 0`. The Wilson 95% interval is [0.850, 0.975]; the model's 0.8559 lies
inside it. Frozen output: `evidence/production-walk.txt`.

Limits, stated plainly. This is the standard random-map heuristic for r-adding
walks; it is not a theorem for this fixed field and this fixed mixer. The
`0.8647` union-bound step is unconditional for a truly uniform random map, but
a hash walk is not independent. The 32-adding construction and its analysis
are the ones in Teske, "On random walks for Pollard's rho method", Math. Comp.
70 (2001), 809-825; Pollard's original method is Pollard, "Monte Carlo methods
for index computations", Math. Comp. 32 (1978), 918-924. The material
assumption is that a salted 32-way partition of a uniform GT element does not
produce anomalous collision statistics at production size. Fresh salts and
increments, the uniform challenge `z`, and the sizeable probability margin are
the support; a fixed seed or a fixed challenge is not claimed to be
representative. If the heuristic failed badly the work cap would still hold,
but the success lower bound would not; that is the disclosed risk.

The mixer `H` is not cryptographic. It is a fixed deterministic function of the
injective 384-byte encoding with a fresh 64-bit salt. Correctness of the attack
does not depend on its unpredictability; only the random-map collision
heuristic uses it, and no adversary chooses `z` after seeing the salt (the salt
is fresh per run).

## 7. What was tried and rejected

The obvious structural next step on this board is to quotient GT by the
sign/inversion automorphism or by the full Frobenius orbit. In our `bn254-gt-dlp`
package `rho-tower-v1` we ran a separate experiment (a toy implementation in
Python over `Z/p`, not shipped) with a walk folded by `y -> min(y, -y)`, the
labels negated whenever the representative is inverted, and
`N = 2*ceil(sqrt((r-1)/2))`. It failed: on primes around 2^17 the folded walk
found a solvable collision only about 8-11% of the time, and 83-89% of the
first collisions had `B - b = 0`. The cause is structural: folding creates
directed two-cycles `y -> y' -> y` whenever the folded image of `y` hashes to
the same jump index, and on such a cycle the labels return exactly, so the
collision is useless. Adding a fixed escape jump when that index repeats raised
the toy success rate to about 47%, still far below the 93% of the unfolded walk,
so we did not build it. The full Frobenius-orbit quotient has the same
two-cycle problem plus a per-step cost of about eleven `Fp12` Frobenius maps
(our estimate: roughly 2x the step), which cancels the `sqrt(12)` gain. Neither
is in this package; both remain open for a solver who handles the fruitless
cycles properly. This is why `rho-tower-v2` improves only the measured
arithmetic.

## 8. Relation to prior work and credits

The algorithm is Pollard rho with a 32-adding walk and Floyd detection. The
tower and the profile basis are the ones published in the target profile; the
tower arithmetic is standard (Karatsuba over `Fp2/Fp6/Fp12`; the six-product
`Fp6` form is the standard Karatsuba/Toom decomposition for `v^3 = xi`). Teske
(2001) is the source of the r-adding walk analysis; Pollard (1978) of the
method. The track's accounting rules (`rules.md`, `COST_MODEL.md`) define the
unit, the inclusion list and the `C = 20,000` normalisation. The implementation,
measurement, arithmetic and text here are ours; no claim, parameter or text
was copied from another solver's submission. We compared the accounting
conventions with the accepted baseline `rho-baseline-v2-20261110` but did not
reuse its code or numbers; our measured step cost and our tower are different,
and our score is an independent bound.

## 9. Reproduction

    cd submissions/bn254-gt-dlp/rho-tower-v2
    gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o rho fp12.c sha256.c rho.c
    ./rho selftest                       # g_is_one=0 g^r_is_one=1 / assoc=1
    ./rho hashdump 4                     # mixer indices 22 2 22 12
    gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o crosscheck crosscheck.c fp12.c
    ./crosscheck                         # 2000 Fp6/Fp12 pairs ok / 2000 add-sub ok
    ./rho mktable 20261112 /tmp/table.bin
    for j in $(seq 0 32); do for K in 1000 3000; do
      valgrind --tool=callgrind --callgrind-out-file=cg.$K.$j.out ./rho stepwc $K $j /tmp/table.bin
    done; done
    python3 experiment.py '{"seed": 20261112, "trials": 24}'

`experiment.py` is standard library only and needs no file outside the package.
It recomputes every row of the ledger from `r`, checks that the declared
rational is the rounded-up `log2(W)`, recomputes the capped sampler's failure
bound, measures the sampler's acceptance rate and attempt counts, verifies that
its mixer transliteration reproduces the binary's `hashdump` indices, and runs
the production walk in the order-43,913 subgroup of the real BN254 Fp12.
Those runs are evidence at their stated scope; they are not a proof of the
collision law at order `r`.
