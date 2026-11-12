# BN254 GT discrete logarithm: a measured 32-adding rho baseline in the tower representation

Track `bn254-gt-dlp`. Input: `z = g^x` in the prime-order group GT of BN254, with
`g` the profile generator `e(P0,Q0)` and `x` uniform in `0..r-1`. Output: the
scalar `x`. One unit is `C = 20,000` scalar AArch64 instruction equivalents.
This is an analytic exploratory claim: no production-size attack was run, and
the success assessment uses the standard random-map heuristic, stated and
bounded below.

## Claim

For the protected BN254 target, the single-trial algorithm of this document
returns `x` with `g^x = z`, with success probability at least **0.85** over
uniform challenges and fresh algorithm coins, and every run performs at most

    W <= 2^130.665164 < 2^130.665165   units,

where one unit is 20,000 scalar AArch64 instructions. This is a per-run
deterministic work cap, not an expectation: every loop in the algorithm has a
data-independent trip count (listed, one by one, in `LOOP_BOUNDS.txt`), the
scalar sampler has a hard cap of K = 64 attempts per scalar, the walk stops
within `N = 2*ceil(sqrt(r))` iterations, and the table build, recovery,
verification and output work are fixed. The score is `log2(W)` by the track
rules, submitted as the rational `130665165/1000000`. Success is at least the
protected 0.39.

**This is a changed algorithm, not the rejected submission.** Submission
`bf1247b4` drew its scalars by rejection sampling with no attempt limit, so for
any finite budget some run could exceed it; the review's first objection was
exactly that, and it was correct. The algorithm below is different in two
declared ways: (i) sampling is a capped rejection sampler over a SHA-256
counter-mode word stream, with a hard cap of K = 64 attempts per scalar, an
explicit failure outcome, its abort probability charged against the success
bound and its capped cost in the ledger; (ii) the per-step price is no longer a
~10%-padded average measurement but the measured worst-case path over every
production table entry and every data-dependent branch (section 4). Everything
else -- the tower, the basis map, the 32-adding walk, Floyd detection, the
recovery equation and the random-map heuristic -- is unchanged.

The claim is a **baseline**, not a best-known attack. It does not use NFS or
index calculus. Its point is that the arithmetic kernel is measured, not
assumed: the per-step cost below is the Valgrind instruction count of the
shipped scalar C implementation.

## 1. Representation and conversion

The profile defines

    Fp12 = Fp[w]/(w^12 - 18 w^6 + 82),   u = w^6 - 9,   Fp2 = Fp[u]/(u^2+1).

The attack uses the isomorphic tower

    Fp2 = Fp[u]/(u^2+1),  Fp6 = Fp2[v]/(v^3 - (9+u)),  Fp12 = Fp6[w_t]/(w_t^2 - v).

`w_t^6 = v^3 = 9+u` and `w_t^12 = (9+u)^2 = 80 + 18u = 18 w_t^6 - 82`, so `w_t`
satisfies the profile polynomial; `w -> w_t` is an Fp-isomorphism. In the tower
basis `{1,u,v,uv,v^2,uv^2} + {1,u,v,uv,v^2,uv^2} w_t`, the profile-basis
coefficients `t_0..t_11` map to

    B=t6, A=t0+9B;  H=t7, G=t1+9H;  D=t8, C=t2+9D;
    J=t9, I=t3+9J;  F=t10, E=t4+9F; L=t11, K=t5+9L,

with tower coordinates `(A+Bu) + (C+Du)v + (E+Fu)v^2 +
((G+Hu) + (I+Ju)v + (K+Lu)v^2) w_t`. This is a fixed public linear map. Its
cost is twelve Fp multiplications and a bounded number of additions; it is
included in the setup row. The profile generator converted this way satisfies
`g^r = 1` in our implementation (the shipped `selftest` prints this), which is a
direct check of the conversion and of the tower arithmetic.

The attack never needs a pairing. It needs only GT multiplication, which is the
tower multiplication in section 4.

## 2. Algorithm

Public input: the profile `p, r, g` and `z`.

1. Decode `z` canonically, check `z != 0` and `z^r = 1`; if `z = 1` output
   `x = 0`. Convert `g` and `z` to the tower (section 1).
2. Draw one 256-bit seed `sigma` from the OS CSPRNG and one fresh 64-bit salt.
   The word stream is `w_i = SHA-256(sigma || LE64(i))` truncated to 256 bits,
   `i = 0, 1, 2, ...`. Every uniform scalar in `0..r-1` is produced by the
   capped rejection sampler

       S:  for k = 1..K,  K = 64:
               v <- next word          (uniform in 0..2^256-1)
               if v < r: return v
           abort: the run stops and reports failure

   Draw the 64 table scalars `c_j, d_j` and then `a_0, b_0` with independent
   calls to `S`, in that fixed order. At most `66*K = 4,224` words are drawn and
   the sampler consumes no other randomness.
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
sampler returns a scalar exactly uniform on `0..r-1`, so the walk analysis of
section 6 is unchanged. One call of `S` exhausts its budget with probability at
most

    delta^K,   delta = 1 - r/2^256 = 0.810969445...,   delta^64 <= 2^-19.35,

and the union bound over the 66 calls is at most `66*delta^64 <= 2^-13.30`
(the bound is recomputed from `r` in `experiment.py`). Success is therefore at
least `0.8559 - 2^-13.30 = 0.855801 > 0.85`. An aborted run has spent only the
sampler's capped cost and stops; the per-run cap therefore holds on the failure
path as well, which is the point the review raised.

## 3. Correctness

The invariant is `Y = g^a z^b` initially (`Y_0 = g^{a_0} z^{b_0}`) and after a
step, because `Y M_j = g^a z^b * g^{c_j} z^{d_j} = g^{a+c_j} z^{b+d_j}`.
If `g^a z^b = g^A z^B` and `z = g^x`, then `a + xb = A + xB mod r`, hence
`a - A = x(B - b) mod r`. Since `r` is prime, a nonzero `B - b` has a unique
inverse and the recovered `x` is the unique discrete logarithm. The final
`g^x = z` check rejects any wrong scalar.

Floyd's equality is a genuine collision of two walk states with (in general)
different step counts. If the deterministic walk has tail `mu` and cycle
`lambda`, an equality is observed at some `k` with `mu <= k <= mu + lambda`, so
`N >= mu + lambda` suffices for detection. The scalar labels are carried along,
not derived from the hash, so the collision yields the linear relation above.

Degenerate collisions (`B - b = 0`) are negligible. Conditional on the group
walk and on the salt, the jump logarithms `s_j = c_j + x d_j mod r` are fixed;
the individual `d_j` remain independent uniform because each `c_j` is uniquely
determined by `s_j`. The coefficient difference at a collision is
`sum_j n_j d_j` with `sum_j n_j = k`, `0 < k <= N < r`, and at least one
`n_j != 0 mod r`; a nonzero integer combination of independent uniform field
elements is exactly uniform, so the degeneracy probability is `1/r`
conditioned on any such walk. This is an exact argument, not a heuristic.

## 4. Arithmetic and measured step cost

All arithmetic is scalar 64-bit. `Fp` elements are four little-endian words in
Montgomery form; `fp_mul` is a 4-limb CIOS Montgomery multiplication with a
six-word accumulator (`fp12.c`). The tower multiplication is

* `Fp2`: three `Fp` multiplications (Karatsuba) and additions;
* `Fp6`: six `Fp2` multiplications plus three multiplications by `xi = 9+u`
  (each four `Fp` additions/doublings);
* `Fp12`: three `Fp6` multiplications and two `Fp6` multiplications by `v`
  (Karatsuba over `w_t^2 = v`).

One step performs exactly one `Fp12` multiplication, one hash of 48 words, two
4-word additions modulo `r`, one table record load and the loop bookkeeping.
The shipped binary `rho` was compiled with

    gcc 13.3.0 -O3 -fno-tree-vectorize -mgeneral-regs-only

on the same machine class as the cost model (AArch64, 64-bit general-purpose
registers). Valgrind 3.22.0 Callgrind (`Ir`, the instruction event) was run on
the step loop.

**Worst case, not average.** The count below is a measured worst-case path, not
an average with a percentage margin. The reason a sweep suffices is structural:
every loop in the step has a data-independent trip count (`fp_mul` 4x4 limbs,
the mixer's 24 rounds over 48 words, two 4-limb scalar additions), so the
instructions executed by a step are determined by the step's inputs only through
the arms of its data-dependent decisions -- the carry and reduction arms of
`fp_add`, the borrow arm of `fp_sub`, the final comparison in `fp_mul` and the
carry/borrow select in `scalar_add_r`. Those arms select a sequence; they never
change an iteration count.

We therefore measured the step over the full product of

* all **32 production table entries** `j = 0..31` (the index is pinned by the
  audit mode `./rho worstcase K j pat`, which still computes the mixer and only
  overrides its result, so the measured code is the shipped step);
* **22 operand patterns** per entry: the production start `g`; all limbs
  `0xffff...` (forces the carry arm in every `Fp` addition); all limbs `1`
  (forces the borrow arm in every `Fp` subtraction); all limbs
  `FP_ONE_MONT`; all limbs `p_limb - 1` (carry plus the full comparison chain);
  all limbs `p_limb` (every addition overflows the modulus); and 16
  deterministic pseudo-random 384-byte elements;
* two loop lengths, K = 1,000 and K = 3,000, differenced to cancel the table
  build and conversion.

That is **704 configurations, 1,408 Callgrind runs**. The per-step count is
`(Ir(3000) - Ir(1000)) / 2000`. Results (full table in
`evidence/worstcase-sweep.txt`):

| statistic | Ir per step |
|---|---:|
| minimum over the 704 configurations | 46,319.00 |
| median | 47,248.14 |
| **maximum** (j = 13, pat = 12) | **47,264.24** |

We charge

    C_STEP = 47,265 instructions per step,

the maximum rounded up to the next instruction. There is no percentage padding:
47,264.24 is the largest value observed on any production table entry under any
of the swept operand patterns, and the patterns include the extremes that force
each arm to its longer direction. `measure.json` records the summary and the
command; `machine.md` gives the reproduction loop.

The measured loop includes the hash, the two scalar additions, the table load
and the multiplication; it excludes only the Floyd comparison and loop control,
charged separately below at 4,096 instructions per iteration. That charge is
analytic: twelve 4-word comparisons are 48 word compares (~3 instructions
each), plus a counter and a branch, so 4,096 is more than 25x the operation
count of the comparison, and it is charged for all `N` iterations.

## 5. Complete work ledger

`r` is the profile subgroup order, `N = 2*ceil(sqrt(r)) =
295893513763578638011461384341992519220`.

| component | work | log2 units (C=20,000) |
|---|---:|---:|
| 3N walk steps at 47,265 instructions | 3N * 47,265 / 20,000 | 130.624080 |
| N Floyd iterations, 4,096 instructions each | N * 4,096 / 20,000 | 125.510633 |
| setup: 64 exponentiations, table, conversion | <= 2^36 instructions | 21.712992 |
| capped sampler: 66 * 64 draws at <= 2,048 instructions | <= 8,650,752 | 8.756682 |
| verification `g^x = z`, canonical output | <= 2^32 instructions | 17.712288 |
| output and final scalar inversion margin | <= 2^20 instructions | 5.712288 |
| **total W** | | **130.665164** |

The setup row: 64 exponentiations with 254-bit exponents, each at most 508
`Fp12` multiplications, is 32,512 multiplications; even at the loose analytic
bound of 2^20 instructions per multiplication used by the earlier published
BN254 baseline, that is below 2^36. We charge 2^36, about 30x the measured
cost of that many multiplications.

The sampler row is charged at its **cap**, not at the observed number of
attempts: 66 scalars times K = 64 attempts times 2,048 instructions per
SHA-256-compression-sized attempt is 8,650,752 < 2^24 instructions, which is
what the ledger adds. The charge is therefore independent of the input and of
the observed rejection counts. (`experiment.py` also measures the sampler: over
20,000 draws the acceptance rate is within 0.02 of `r/2^256`, and over 2,000
scalars the largest observed attempt count is well below the cap.) The verification row is two exponentiations
plus decoding, charged at 2^32, over 60x the same estimate. The output row
covers the single scalar inversion (`(B-b)^{-1}`, at most 2^22 instructions by
Fermat) and the 32-byte encoding with an additional margin.

No row is amortised over future targets, no oracle is free, and no BN254
operation is priced at one unit. Memory: 32 table records of 448 bytes, two
walk states, the hash scratch and code, below 2^20 bytes; the claim declares
`memory_log2_bytes = 20` and `nonuniform_advice_log2_bytes = 0` (the salt and
increments are fresh per run and charged in setup).

## 6. Success probability: heuristic, support and limits

Model the walk as a random map on the `r` points of GT, with the salt and the
64 increments as the fresh coins. For `t` distinct transitions the no-repeat
probability is `prod_{i=1}^{t-1}(1-i/r) <= exp(-t(t-1)/(2r))`. At
`t = N = 2 ceil(sqrt(r))` this is at most `exp(-2) = 0.1353`, so a repeat occurs
with probability at least `0.8647` even before the 32-adding correction. The
conventional r-adding correction multiplies the effective collision rate by
`1 - 1/32`; applying it at `t = 2 sqrt(r)` gives

    1 - exp(-(31/32)*2) = 0.8559.

We claim the smaller, conservative **0.85**, and the capped sampler's abort
probability `2^-13.30` is subtracted from the collision model in section 2, so
the claimed floor is `0.8559 - 2^-13.30 = 0.85580`.

**Measured success of the production walk.** The review objected that the
earlier experiment checked the formula and ran a small group with a *different*
mixer. Both points are addressed directly. `experiment.py` now runs the shipped
walk: the shipped 384-byte encoding, the mixer transliterated from `rho.c`
(cross-checked index for index against the binary's `hashdump` mode, see
`evidence/mixer-crosscheck.txt`) and the shipped tower multiplication, inside
the order-`Q` subgroup of the *real* BN254 Fp12, where `Q = 43,913` is a
prime factor of `p^12 - 1`. The subgroup order is the only thing that differs
from production; the field, the tower, the encoding, the mixer and the walk are
the production ones.

With seed 20261112 and 64 trials the measured result is **60 of 64 runs
recovered `x` exactly**, 4 runs had no collision inside `N = 2*ceil(sqrt(Q)) = 420`,
and no run produced a wrong scalar or a degenerate `B - b = 0`. The Wilson 95%
interval for the success rate is [0.850, 0.975]; the model's 0.8559 sits inside
it. The frozen output is `evidence/production-walk.txt`. This is evidence at a
toy subgroup order, not a proof at order `r`: it tests the walk and the mixer at
production structure and small size, which is exactly the scope the review said
was missing.

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
is fresh per run). We state this rather than claim properties the mixer does
not have.

## 7. Relation to prior work and credits

The algorithm is Pollard rho with a 32-adding walk and Floyd detection. The
tower and the profile basis are the ones published in the target profile; the
tower arithmetic is standard (Karatsuba over `Fp2/Fp6/Fp12`). Teske (2001) is
the source of the r-adding walk analysis; Pollard (1978) of the method. The
track's accounting rules (`rules.md`, `COST_MODEL.md`) define the unit, the
inclusion list and the `C = 20,000` normalisation. The implementation,
measurement, arithmetic and text here are ours; no claim, parameter or text
was copied from another solver's submission. We compared the accounting
conventions with the accepted baseline `rho-baseline-v2-20261110` but did not
reuse its code or numbers; our measured step cost and our tower are different,
and our score is an independent bound.

## 8. Reproduction

    cd submissions/bn254-gt-dlp/rho-tower-v1
    gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o rho fp12.c rho.c
    ./rho selftest                # g_is_one=0 g^r_is_one=1 / assoc=1
    ./rho hashdump 4              # mixer indices 22 2 22 12 (cross-checked in experiment.py)

    # worst-case step sweep (704 configurations; evidence/worstcase-sweep.txt)
    for j in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31; do
      for pat in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21; do
        for K in 1000 3000; do
          valgrind --tool=callgrind --callgrind-out-file=cg.$K.$j.$pat.out \
            ./rho worstcase $K $j $pat
        done
      done
    done

    python3 experiment.py '{"seed": 20261112, "trials": 24}'

`experiment.py` (standard library only, no files outside this directory)
recomputes every ledger row from `r`, checks that the declared rational is the
rounded-up `log2(W)`, recomputes the capped sampler's failure bound, measures
the sampler's acceptance rate and attempt counts, verifies that its mixer
transliteration reproduces the shipped binary's `hashdump` indices, and runs the
production walk (shipped encoding, mixer and tower) in the order-43,913
subgroup of the real BN254 Fp12. Those runs are evidence at their stated scope;
they are not a proof of the collision law at order `r`.

Files: `LOOP_BOUNDS.txt` (every loop and its hard bound), `evidence/worstcase-sweep.txt`
(704 measured configurations), `evidence/production-walk.txt` (frozen 64-trial
output), `evidence/mixer-crosscheck.txt` (C vs Python mixer).
