# BN254 GT discrete logarithm: the rho-tower baseline with a constant-time step and a static instruction bound

Track `bn254-gt-dlp`. Input: `z = g^x` in the prime-order group GT of BN254,
with `g` the profile generator `e(P0,Q0)` and `x` uniform in `0..r-1`. Output:
the scalar `x`. One unit is `C = 20,000` scalar AArch64 instruction equivalents.
This is an analytic exploratory claim: no production-size attack was run, and
the success assessment is the standard random-map heuristic, stated and bounded
below.

This package is the `rho-tower-v2` kernel (six-product Karatsuba `Fp6`,
branchless inline `Fp` add/sub) with the **last data-dependent control flow
removed from the walk step** and the per-step price turned into a **static
instruction count taken from the compiled step**, not a measurement. It moves
the declared score from 130.131912 (v2) and 130.851718 (v1) to 130.091157.

Changes relative to our two earlier submissions on this track, both rejected
(`bf1247b4`, `b75580bd`), are exactly these:

1. **The step contains no conditional branch at all.** The final reduction
   comparison of `fp_mul` (`lt4` plus `if`) became a mask select, the `if` in
   `scalar_add_r` became a mask select, the audit pin `if (PIN_J < 32)` and the
   measurement checksum were removed from the shipped step, and the 24-round
   mixer loop is unrolled. `fp_add`/`fp_sub` are the branchless inline forms
   already used by v2.
2. **The price is a static count, not a measurement.** `static_count.py`
   disassembles the shipped binary, proves by construction that every function
   in the step's call graph is straight-line (no conditional branch, no loop,
   no early return), and sums instruction counts times the fixed call
   multiplicities: **31,302 instructions per step**, exactly.
3. **Wording accepted from both reviews.** The uniformity, independence,
   degeneracy and abort statements are now stated as *exact conditional on the
   SHA-256 counter-mode word stream being indistinguishable from independent
   uniform words* (a standard PRF assumption, fresh seed per run), not as
   unconditional facts. The success probability is stated as a random-map
   heuristic supported by small-subgroup evidence at `Q = 43,913`, not as a
   production measurement.

## Claim

For the protected BN254 target, the single-trial algorithm of this document
returns `x` with `g^x = z`, with success probability at least **0.85** over
uniform challenges and fresh algorithm coins, and every run performs at most

    W < 2^130.091157   units,   log2 W = 130.091156869,

where one unit is 20,000 scalar AArch64 instructions. The score is
`130091157/1000000`. This is a per-run deterministic work cap: every loop has a
data-independent trip count (`LOOP_BOUNDS.txt`), the sampler has a hard cap of
K = 64 attempts per scalar, the walk stops within `N = 2*ceil(sqrt(r))`
iterations, the step costs exactly 31,302 instructions on every input
(`STATIC_BOUND.txt`), and the table build, recovery, verification and output
work are fixed. Success is at least the protected 0.39.

The claim is a **baseline**, not a best-known attack. It does not use NFS or
index calculus. Its point is that the arithmetic is constant-time, so the
per-step price is an exact static number reproducible from the shipped binary.

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

**Capped sampling and its failure accounting.** The sampler accepts a word
strictly when `v < r` and rejects `v >= r`, including `v = r`; in the ideal
independent-word model of section 3 it therefore returns a scalar exactly
uniform on `0..r-1`, and under the PRF assumption this holds up to the
distinguishing advantage `eps` stated there. One call of `S` exhausts its budget with probability at
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
elements is exactly uniform in the ideal independent-word model, so the
degeneracy probability there is `1/r` conditioned on any such walk; under the
PRF assumption this is the same statement up to `eps` (section 3).

**Two models, stated separately.** In the *ideal model* the 66*64 words are
independent uniform 256-bit strings. There the rejection sampler is exactly
uniform on `0..r-1` (it accepts strictly `v < r` and rejects `v >= r`,
including `v = r`), the `1/r` degeneracy argument is exact, and the `delta^K`
abort bound is exact.

The shipped implementation is not that model: it expands one 256-bit seed by
SHA-256 in counter mode (FIPS 180-4), which produces at most `2^256` possible
66-scalar tuples while independent uniform draws range over `r^66` tuples. The
two distributions are therefore **not identical**, and we do not claim they
are. We assume only that the expanded stream is computationally
indistinguishable from independent uniform words (the standard PRF assumption
for SHA-256 in counter mode under a secret, fresh, never-reused seed). Under
that assumption the ideal-model statements hold up to a distinguishing
advantage `eps`, and the sampler-derived bounds in this document are
**heuristic with additive error `eps`**, not exact facts about the
deterministic expansion. `eps` is not quantified; the claim is that it is
negligible for a fresh random seed, which is an assumption and is labelled as
one. The work cap does not depend on this at all: the step price is exact
(section 4) and the sampler is capped either way.

## 4. Constant-time arithmetic and the exact static step price

All arithmetic is scalar 64-bit. `Fp` elements are four little-endian words in
Montgomery form; `fp_mul` is a 4-limb CIOS Montgomery multiplication with a
six-word accumulator. Relative to `rho-tower-v1` the kernel keeps the v2
improvements (six `Fp2` products per `Fp6` product instead of nine; no dead
`2a/3a/6a` terms in multiplication by `xi = 9+u`; inline branchless
`fp_add`/`fp_sub`) and removes the remaining data-dependent control flow:

1. **`fp_mul` final reduction (fp12.c lines 47-58).** v1/v2 computed
   `d = t - p` and then chose between `d` and `t` with `if (t[4] ||
   !lt4(t, P)) memcpy(...)`. Now `sub4` returns the borrow, `cond = t[4] |
   (borrow ^ 1)` is 1 exactly when `t >= p`, and the four output limbs are a
   mask select `(d & mask) | (t & ~mask)`. Four `and`, four `orr`, one `neg`,
   one `sub`, no branch.
2. **`scalar_add_r` (rho.c line 86).** Same change: `cond = s[4] | (borrow ^ 1)`
   selects `s - r` or `s` under a mask instead of `if (...) memcpy(...)`.
3. **The shipped step (rho.c, `rho_step`).** It computes `j = hash_index(y,
   salt)` with the 24 rounds unrolled, indexes `tab[j]` directly, multiplies
   `y` by `tab[j].m` and adds the two labels. The audit pin `if (PIN_J < 32)`
   and the measurement checksum were removed: they were instruments, not
   algorithm. What remains has no conditional branch and no data-dependent
   index arithmetic beyond the direct `tab + j` load.
4. **`fp6_add`/`fp6_sub` are written out limb-wise** so the compiler cannot
   keep a loop back-edge in `fp12_mul`; all fixed loops in the kernel
   (4-limb CIOS, 4-limb add/sub, 24-round mixer) are unrolled by the compiler
   at `-O3`, and the disassembly assertion in `static_count.py` checks that
   none survived.

The binary was compiled with

    gcc 13.3.0 -O3 -fno-tree-vectorize -fno-stack-protector -mgeneral-regs-only

on the same machine class as the cost model (AArch64, 64-bit general-purpose
registers).

**The bound is a static count of the compiled step, not a measurement.**
`static_count.py` (shipped) parses `objdump -d --no-show-raw-insn ./rho`,
takes the transitive call graph of the root symbol `rho_step` (marked
`noinline, noipa`, so it is the generic step taking an arbitrary 64-bit salt),
asserts that every reached function contains no conditional branch
(`b.cond`, `cbz`, `cbnz`, `tbz`, `tbnz`), no backwards branch and no indirect
branch, and that no call leaves the binary; then it counts each function's
instructions and multiplies by the number of times the call graph enters it.
`STATIC_BOUND.txt` is the recorded run. The result:

| function | instructions | calls per step | instructions per step |
|---|---:|---:|---:|
| `rho_step` (inlined mixer, 2x `scalar_add_r`, call site, `ret`) | 317 | 1 | 317 |
| `fp12_mul` | 2,581 | 1 | 2,581 |
| `fp6_mul` | 2,982 | 3 | 8,946 |
| `fp2_mul` | 274 | 18 | 4,932 |
| `fp_mul` | 269 | 54 | 14,526 |
| **total** | | | **31,302** |

The call multiplicities are fixed by the source: `fp12_mul` makes three
`fp6_mul` calls (Karatsuba over `w_t^2 = v`), each `fp6_mul` makes six
`fp2_mul` calls (six-product Karatsuba over `v^3 = xi`), and each `fp2_mul`
makes three `fp_mul` calls. Nothing about the data changes a count, because
there is no conditional control transfer to change it.

    C_STEP = 31,302 instructions per step, exact.

No padding, no average, no worst-case arm: the same 31,302 instructions execute
for every table entry, every walk state and every challenge.

**Corroboration, not justification.** `evidence/callgrind-identity.txt` runs
Valgrind Callgrind 3.22.0 with `--toggle-collect='*rho_step*'` on ten
independent production tables built by `./rho mktable <seed>` (the capped
sampler of section 2, `M_j = g^c_j z^d_j`), 1,000 steps each, plus 2,000 and
5,000 steps on one table; `stepwc` varies the salt on every step, so this
audits the arbitrary-salt step rather than a constant-specialized clone.
Collected `Ir` is 31,302,000 for 1,000 steps and
156,510,000 for 5,000 on every table: per-step `Ir` is **31,302 in every run**,
identical to the static count. A difference between any two steps would mean
the step is not data-independent; there is none.

The Floyd comparison and loop control are outside the step and are charged
separately at 4,096 instructions per iteration (48 word comparisons, three
`bl` call instructions, counter and branch), an analytic over-charge that is
not used to justify anything inside the step.

## 5. Complete work ledger

`r` is the profile subgroup order, `N = 2*ceil(sqrt(r)) =
295893513763578638011461384341992519220`.

| component | work | log2 units (C=20,000) |
|---|---:|---:|
| 3N walk steps at 31,302 instructions (static) | 3N * 31,302 / 20,000 | 130.029563 |
| N Floyd iterations, 4,096 instructions each | N * 4,096 / 20,000 | 125.510633 |
| setup: 64 exponentiations, table, conversion | <= 2^36 instructions | 21.712992 |
| capped sampler: 66 * 64 draws at <= 2,048 instructions | <= 8,650,752 | 8.756682 |
| verification `g^x = z`, canonical output | <= 2^32 instructions | 17.712288 |
| output and final scalar inversion margin | <= 2^20 instructions | 5.712288 |
| **total W** | | **130.091157 -> declared 130.091157** |

The sampler row is charged at its cap (66 scalars x 64 attempts x 2,048
instructions = 8,650,752 < 2^24), not at an observed attempt count, so it is
input-independent. The walk row is the product of the exact static step price
and `3N`; `N` itself is fixed before the run.

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

## 6. Success probability: a random-map heuristic, with small-subgroup evidence

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

**Evidence at small order, not a production measurement.** `experiment.py`
runs the shipped walk: the shipped 384-byte encoding, the mixer transliterated from `rho.c`
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
cycles properly. This is why this package improves only the arithmetic and
the accounting of it.

## 8. Relation to prior work and credits

The algorithm is Pollard rho with a 32-adding walk and Floyd detection. The
review response (constant-time kernel, static disassembly count, PRF wording)
follows the two decisions on our earlier submissions `bf1247b4`/`b75580bd` and
`d887db59`/`d54d32f6`. The
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

    cd submissions/bn254-gt-dlp/rho-tower-v3
    gcc -O3 -fno-tree-vectorize -fno-stack-protector -mgeneral-regs-only \
        -o rho fp12.c sha256.c rho.c
    ./rho selftest                       # g_is_one=0 g^r_is_one=1 / assoc=1
    ./rho hashdump 4                     # mixer indices 22 2 22 12
    gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o crosscheck crosscheck.c fp12.c
    ./crosscheck                         # 2000 Fp6/Fp12 pairs ok / 2000 add-sub ok
    python3 static_count.py ./rho        # asserts straight-line, prints 31302
    ./rho mktable 20261112 /tmp/table.bin
    valgrind --tool=callgrind --toggle-collect='*rho_step*' \
        --callgrind-out-file=/tmp/cg.out ./rho stepwc 1000 /tmp/table.bin
    python3 experiment.py '{"seed": 20261112, "trials": 24}'

`static_count.py` is the static audit: it fails (exit 1) if any function in the
step's call graph contains a conditional branch, a loop or an indirect branch,
and otherwise prints the per-function counts and the 31,302 total.
`experiment.py` is standard library only and needs no file outside the package.
It recomputes every row of the ledger from `r` and the static step price,
checks that the declared rational is the rounded-up `log2(W)`, recomputes the
capped sampler's failure bound, measures the sampler's acceptance rate and
attempt counts, verifies that its mixer transliteration reproduces the
binary's `hashdump` indices, and runs the 32-adding walk in the order-43,913
subgroup of the real BN254 Fp12. Those runs are evidence at their stated scope;
they are not a proof of the collision law at order `r`.
