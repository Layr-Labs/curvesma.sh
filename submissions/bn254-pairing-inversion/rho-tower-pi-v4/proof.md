# BN254 pairing inversion by GT discrete logarithm, tower rho and a G1 scalar multiplication

Track `bn254-pairing-inversion`. The protected profile defines `g = e(P0,Q0)` of
order `r` in the BN254 target group GT and, for uniform `x` in `0..r-1`, gives
only `z = g^x`. The required output is any `P` in G1 and `Q` in G2 with
`e(P,Q) = z`.

The track rules state the input and success conditions directly:

> "Given only z=g^x for uniform x in 0..r-1, return any P in G1 and Q in G2
> such that e(P,Q)=z. Neither point is supplied or fixed by the challenge;
> choosing a public generator for one point is allowed."

and name this reduction explicitly:

> "A GT discrete-logarithm solver gives a strategy: recover x and return P0 and
> [x]Q0, charging scalar multiplication and all other work."

We take the transposed form of that strategy: recover `x`, return `P = [x]P0`
and `Q = Q0`. By bilinearity of the reduced pairing and `g = e(P0,Q0)`,

    e(P, Q0) = e([x]P0, Q0) = e(P0, Q0)^x = g^x = z,

so the pair is valid exactly when `x` is the discrete logarithm of `z` to base
`g`. This is the same reduction the rules describe with the factors exchanged;
the G1 scalar multiplication is cheaper than the G2 one and is charged in full.

One unit is `C = 20,000` scalar AArch64 instruction equivalents and the score is
`log2(W)` with `W = R/C` for total scalar instruction-equivalent work `R`, as in
`rules.md` and `../docs/COST_MODEL.md`. This is an analytic exploratory claim:
no production-size attack was run.

## Claim

For the protected BN254 target, the single-trial algorithm of this document
outputs `(P, Q)` with `e(P,Q) = z` for the profile pairing with success
probability at least **0.85** over uniform challenges and fresh algorithm coins,
and every run performs at most

    W < 2^129.988202   units,   log2 W = 129.988201857,

The score is `129988202/1000000`. Success is far above the protected minimum
0.39. The claim is a per-run deterministic work cap: one trial, one table, one
salt, one **capped** sampler (K = 64 attempts per scalar), a fixed number of
walk steps of exactly 29,052 instructions each, one scalar multiplication
charged at a fixed 2^24, and bounded verification and output. Every loop in the
algorithm outside the walk step has a data-independent trip count; the list is
`LOOP_BOUNDS.txt`.

**This is a changed algorithm, not the rejected submission.** Submission
`d887db59` specified rejection sampling without an attempt limit, so no finite
per-run cap followed for the algorithm as written; the review was right about
that, and adding the cap and its failure accounting is a change to the declared
algorithm. The algorithm below draws its scalars from a SHA-256 counter-mode
word stream through a capped rejection sampler that aborts the run after K = 64
attempts per scalar; the abort probability is charged against the success
bound and the sampler's capped cost is in the ledger. The per-step price is also
replaced: it is now an exact static instruction count of the compiled step,
with no data-dependent branch left to vary it.

This is a **baseline, not a best-known attack**, and it is not an improvement
over the best public claims on this board (see section 9). Its purpose is to
give the inversion track the same exact floor that our `bn254-gt-dlp` package
has. Relative to the two rejected submissions on this track (`d887db59`,
`d54d32f6`) this package changes exactly three things:

1. **The GT walk step is constant-time.** The final reduction comparison of
   `fp_mul` and the conditional add of `r` in `scalar_add_r` are mask selects,
   the audit pin and the measurement checksum are gone from the shipped step,
   `fp6_add`/`fp6_sub` are written limb-wise and the mixer is unrolled. The
   compiled step and every function it calls contain no conditional branch, no
   loop and no early return.
2. **The price is a static instruction count, not a measurement.**
   `static_count.py` disassembles the shipped binary, asserts the straight-line
   property, and sums instruction counts times fixed call multiplicities:
   **29,052 instructions per step**, exactly (`STATIC_BOUND.txt`), corroborated
   by identical Callgrind `Ir` on every step of ten production tables.
3. **Wording.** Sampler uniformity, independence, degeneracy and abort
   probabilities are exact in the ideal independent-word model and heuristic up
   to eps under the standard PRF assumption for the
   SHA-256 counter-mode word stream (fresh seed per run), and the 0.85 success
   figure is a random-map heuristic supported by small-subgroup evidence at
   `Q = 43,913`, not a production measurement.

The kernel is the one shipped in our `bn254-gt-dlp` `rho-tower-v3` package
(six-product Karatsuba `Fp6`, branchless inline `Fp` add/sub); the tower, the
walk, the bilinear reduction and the G1 scalar multiplication are unchanged.

## 1. Representation and conversion

The profile defines

    Fp12 = Fp[w]/(w^12 - 18 w^6 + 82),   u = w^6 - 9,   Fp2 = Fp[u]/(u^2+1).

We use the isomorphic tower

    Fp2 = Fp[u]/(u^2+1),  Fp6 = Fp2[v]/(v^3 - (9+u)),  Fp12 = Fp6[w_t]/(w_t^2 - v),

with the fixed public linear basis map of the `bn254-gt-dlp` package we wrote
(`w_t^6 = 9+u`, `w_t^12 = 18 w_t^6 - 82`, so `w -> w_t` is an Fp-isomorphism).
In profile coefficients `t_0..t_11` the map is

    B=t6, A=t0+9B;  H=t7, G=t1+9H;  D=t8, C=t2+9D;
    J=t9, I=t3+9J;  F=t10, E=t4+9F; L=t11, K=t5+9L,

with tower coordinates `(A+Bu) + (C+Du)v + (E+Fu)v^2 +
((G+Hu) + (I+Ju)v + (K+Lu)v^2) w_t`. It costs twelve Fp multiplications and a
bounded number of additions, charged in the setup row. The shipped `selftest`
checks `g^r = 1` after conversion, and also `P0` on the curve and `[r]P0 = O`.

G1 is the profile curve `y^2 = x^3 + 3` over Fp with `P0 = (1,2)`. G1 has
cofactor 1, so every on-curve point is in the order-`r` subgroup; we still check
`[r]P = O` for the output point because the acceptance rule requires it.

## 2. Algorithm

Public input: the profile and `z`.

1. Decode `z` canonically; check `z != 0` and `z^r = 1`; if `z = 1`, output
   `P = null` (infinity) and `Q = Q0`, which is accepted since
   `e(O,Q0) = 1 = z`.
2. Convert `g` and `z` to the tower (section 1). Draw one 256-bit seed `sigma`
   from the OS CSPRNG and one fresh 64-bit salt.
3. The word stream is `w_i = SHA-256(sigma || LE64(i))` truncated to 256 bits.
   Every uniform scalar in `0..r-1` is produced by the capped rejection sampler

       S:  for k = 1..K,  K = 64:
               v <- next word          (uniform in 0..2^256-1)
               if v < r: return v
           abort: the run stops and reports failure

   Draw the 64 table scalars `c_j, d_j` and then `a_0, b_0` with independent
   calls to `S`, in that fixed order (at most 66*K = 4,224 words). Build the
   32-entry table `M_j = g^{c_j} z^{d_j}` with both scalar labels. Two
   fixed-base/fixed-exponent exponentiations per entry.
4. Choose `a_0, b_0` uniform in `0..r-1`, set `Y_0 = g^{a_0} z^{b_0}`.
5. Walk `(Y,a,b) -> (Y*M_j, a+c_j, b+d_j)` with
   `j = H(salt || enc(Y)) mod 32`, where `enc(Y)` is the 384-byte concatenation
   of the twelve Montgomery coefficients and `H` is the shipped 64-bit mixer.
   The invariant `Y = g^a z^b` holds after every step.
6. Floyd detection for `k = 1..N`, `N = 2*ceil(sqrt(r))`; at the first equality
   with tortoise label `(a,b)` and hare label `(A,B)`, if `B-b = 0` report
   failure, else

       x = (a - A) * (B - b)^{-1} mod r.

7. **Inversion step.** Compute `P = [x]P0` on G1 by left-to-right
   double-and-add over the 256-bit scalar `x`, using Jacobian doubling, mixed
   addition with the affine fixed base `(1,2)`, and one Fp inversion for the
   affine output. Set `Q = Q0` (the profile's G2 generator, converted as the
   profile defines; no arithmetic on it is needed).
8. Verify `g^x = z` in the tower. Check `P` is on `y^2 = x^3+3` and
   `[r]P = O`, i.e. the two acceptance conditions for the output points, and
   that `Q0` is the profile generator. A mismatch is failure, never an accepted
   answer. Output `{P: {x,y} or null, Q: {x,y}}` in the profile encoding.

There is exactly one trial, one table, one salt, and no restart. `N` is fixed
before the run; every failure path spends at most the same bounded work. The
`P` computation runs once, after `x` is known; it is not amortised.

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

### Why verification without a pairing is enough

The algorithm does not compute a pairing. It checks `g^x = z` in the tower,
which is the same field and generator the profile pins, and the output relation
then follows from bilinearity of the profile's reduced pairing:

    e(P,Q0) = e([x]P0,Q0) = e(P0,Q0)^x = g^x = z.

Bilinearity and non-degeneracy of the reduced ate pairing are the profile's
definition; the algorithm uses no property beyond that identity. If the shipped
tower conversion or the G1 arithmetic were wrong, the `g^x = z` and
`[r]P = O` checks are what fails; the two checks together with the identity
above cover the acceptance conditions `[r]P = [r]Q = O` and `e(P,Q) = z`, with
`Q = Q0` fixed public data.

## 3. Correctness

The GT part is the argument of our `bn254-gt-dlp` package, restated. The
invariant is `Y = g^a z^b`; a collision gives `a + xb = A + xB mod r`, so
`a - A = x(B-b) mod r`; `r` is prime, so a nonzero `B-b` has a unique inverse
and yields the unique `x`. Floyd's equality is a genuine collision of walk
states; the scalar labels travel with the states, so it yields the linear
relation. The degenerate case `B-b = 0` is a failure and is charged.

**Two models, stated separately.** In the *ideal model* the 66*64 words are
independent uniform 256-bit strings; there the rejection sampler is exactly
uniform (strict `v < r`), the `1/r` degeneracy argument is exact and the
`delta^K` abort bound is exact. The shipped implementation expands one 256-bit
seed by SHA-256 in counter mode, producing at most `2^256` possible 66-scalar
tuples while independent uniform draws range over `r^66`, so the distributions
are **not identical** and we do not claim they are. We assume only
computational indistinguishability from independent uniform words (standard PRF
assumption, fresh never-reused seed). Under it the ideal-model statements hold
up to a distinguishing advantage `eps`, and the sampler-derived bounds are
heuristic with additive error `eps`, not exact facts about the expansion. The
work cap does not depend on this: the step price is exact and the sampler is
capped either way.

For the G1 part, standard Jacobian formulas with `a = 0`: the mixed addition is
`madd-2007-bl` for the affine base, doubling is `dbl-2009-l`, and the
double-and-add loop computes `[x]P0` for the 256-bit scalar (the loop starts at
bit 255; `x < r < 2^254`, so leading zero bits only cost doublings of the point
at infinity). The affine conversion uses one Fermat inversion `Z^{p-2}`. The
`[r]P = O` and on-curve checks are independent of the scalar multiplication
routine. On the shipped vectors in `vectors.json`, the C output equals an
independent Python affine implementation for `x = 1, 2, 3, 7,
0x123456789abcdef` and `r-1`.

## 4. Constant-time arithmetic and the exact static step price

All arithmetic is scalar 64-bit. `Fp` elements are four little-endian words in
Montgomery form; `fp_mul` is a 4-limb CIOS Montgomery multiplication with a
six-word accumulator (`fp12.c`). The tower multiplication is `Fp2` Karatsuba
(three `Fp` multiplications), `Fp6` with six `Fp2` products (six-product
Karatsuba over `v^3 = xi`, `xi = 9+u`), and `Fp12` with three `Fp6` products
and two multiplications by `v`.

One GT step performs one `Fp12` multiplication, one 48-word mixer, two 4-word
additions modulo `r` and one direct table-record load. The binary `rho_pi` was
compiled with

    gcc 13.3.0 -O3 -fno-tree-vectorize -fno-stack-protector -mgeneral-regs-only

on AArch64 (Cortex-X925, Linux 6.8.0), the machine class of the cost model.

**Every data-dependent branch is gone.** The conditional subtract at the end of
`fp_mul` (`sub4` + mask select, `fp12.c` lines 47-58) and the conditional add
of `r` in `scalar_add_r` (mask select, `rho_pi.c`) are instruction sequences,
not branches. The audit pin `if (PIN_J < 32)` and the measurement checksum were
removed from the shipped step; `fp6_add`/`fp6_sub` are written limb-wise and the
24-round mixer is unrolled. The only control transfers in the step are the
fixed `bl` call sites.

**The bound is a static count of the compiled step, not a measurement.**
`static_count.py` (shipped) parses `objdump -d --no-show-raw-insn ./rho_pi`,
takes the transitive call graph of the root symbol `rho_step.constprop.0`,
asserts that every reached function contains no conditional branch, no
backwards branch and no indirect branch (and that no call leaves the binary),
then sums each function's instruction count times its fixed call multiplicity.
`STATIC_BOUND.txt` is the recorded run:

| function | instructions | calls per step | instructions per step |
|---|---:|---:|---:|
| `rho_step` (inlined mixer, 2x `scalar_add_r`, call site, `ret`) | 312 | 1 | 312 |
| `fp12_mul` | 2,581 | 1 | 2,581 |
| `fp6_mul` | 2,982 | 3 | 8,946 |
| `fp2_mul` (fused lazy reduction) | 956 | 18 | 17,208 |
| **total** | | | **29,052** |

The multiplicities are fixed by the source (3 `fp6_mul` per `fp12_mul`, 6
`fp2_mul` per `fp6_mul`), and no conditional control transfer exists to change
them. `fp2_mul` is a fused leaf: three 4x4 products per Fp2 product and one
Montgomery reduction per output coefficient, replacing three canonical
Montgomery multiplications (`fp_mul` survives only in setup, outside the step).

    C_STEP = 29,052 instructions per step, exact.

**Corroboration, not justification.** `evidence/callgrind-identity.txt` runs
Callgrind with `--toggle-collect='*rho_step*'` on ten independent production
tables built by `./rho_pi mktable <seed>` (capped sampler, `M_j = g^c_j z^d_j`),
1,000 steps each, plus 2,000 and 5,000 steps on one table: `Ir` per step is
29,052 in every run, identical to the static count. Any difference between two
steps would mean the step is not data-independent; there is none.

The Floyd comparison and loop control are outside the step and charged
separately at 4,096 instructions per iteration. The G1 part (Jacobian doubling,
mixed addition, one 254-bit left-to-right scalar multiplication, affine
inversion by Fermat) is deterministic once `x` is known; each of its functions
has a fixed instruction count and the loop runs 256 times, so it is charged at
a fixed 2^24, far above its actual cost, with no dependence on the Hamming
weight of `x` entering the total. The `g^x = z` and `[r]P0 = O` checks are
charged in the verification row.

## 5. Complete work ledger

`r` is the profile subgroup order, `N = 2*ceil(sqrt(r)) =
295893513763578638011461384341992519220`.

| component | work | log2 units (C=20,000) |
|---|---:|---:|
| 3N GT walk steps at 29,052 instructions (static) | 3N * 29,052 / 20,000 | 129.921946 |
| N Floyd iterations, 4,096 instructions each | N * 4,096 / 20,000 | 125.510633 |
| setup: 64 GT exponentiations, table, conversion | <= 2^36 instructions | 21.712992 |
| capped sampler: 66 * 64 draws at <= 2,048 instructions | <= 8,650,752 | 8.756682 |
| GT verification `g^x = z`, canonical output | <= 2^32 instructions | 17.712288 |
| G1 `[x]P0`, `[r]P` and on-curve checks, affine inversion | <= 2^24 instructions | 13.712288 |
| output and final scalar inversion margin | <= 2^20 instructions | 5.712288 |
| **total W** | | **129.988202 -> declared 129.988202** |

The total is dominated by the walk row: adding the G1 row changes `log2(W)` by
less than `2^-79`. The sampler row is charged at its cap (66 scalars x 64
attempts x 2,048 instructions per SHA-256-compression-sized attempt =
8,650,752 < 2^24), not at an observed attempt count, so it is
input-independent. Setup, verification and output rows use the same
conservative charges as our GT package (2^36 for 64 exponentiations, 2^32 for
verification, 2^20 for output), each far above the measured cost. No row is amortised over future targets, no oracle
is free, and no BN254 operation is priced at one unit.

Memory: 32 GT table records of 448 bytes, two GT walk states, the hash scratch,
G1 temporaries and code, below 2^20 bytes. No advice; the salt and increments
are fresh per run and charged in setup.

## 6. Success probability: a random-map heuristic, with small-subgroup evidence

The only probabilistic part is the GT rho collision. Model the walk as a random
map on the `r` points of GT with the salt and 64 increments as fresh coins. For
`t` distinct transitions the no-repeat probability is at most
`exp(-t(t-1)/(2r))`; at `t = N = 2 ceil(sqrt(r))` this is `exp(-2) = 0.1353`, so
a repeat occurs with probability at least `0.8647`. The conventional r-adding
correction multiplies the collision rate by `1-1/32`, giving
`1 - exp(-(31/32)*2) = 0.8559`; we claim the smaller **0.85**, with the capped
sampler's abort probability `2^-13.30` subtracted from the collision model in
section 2, so the claimed floor is `0.8559 - 2^-13.30 = 0.85580`.

**Evidence at small order, not a production measurement.** The review noted that the earlier
experiment exited with `FileNotFoundError` for `vectors.json` and that its small
walk used a different mixer. Both are fixed. `experiment.py` needs no file
outside this directory: it regenerates the G1 vectors from `y^2 = x^3 + 3` and
`P0 = (1,2)`, and additionally compares them with `vectors.json` when that file
is present (it is shipped, and the comparison passes). For the GT part it runs
the shipped walk: the shipped 384-byte encoding, the mixer transliterated from
`rho_pi.c` (cross-checked index for index against the binary's `hashdump` mode,
`evidence/mixer-crosscheck.txt`) and the shipped tower multiplication, inside
the order-`Q` subgroup of the *real* BN254 Fp12, `Q = 43,913` a prime factor
of `p^12 - 1`. The subgroup order is the only difference from production.

With seed 20261112 and 64 trials the measured result is **60 of 64 runs
recovered `x` exactly**, 4 runs had no collision inside `N =
2*ceil(sqrt(Q)) = 420`, and no run produced a wrong scalar or a degenerate
`B - b = 0`. The Wilson 95% interval is [0.850, 0.975]; the model's 0.8559 lies
inside it. Frozen output: `evidence/production-walk.txt`.

Limits, stated plainly. This is the standard random-map heuristic for r-adding
walks (Teske 2001), not a theorem for this fixed field and mixer. The `0.8647`
union-bound step is unconditional for a truly uniform random map; a hash walk
is not independent. The material assumption is that a salted 32-way partition
of a uniform GT element does not produce anomalous collision statistics at
production size. The G1 step is deterministic and always correct once `x` is;
the identity case `z=1` is handled exactly. If the heuristic failed badly, the
work cap would still hold but the success lower bound would not; that is the
disclosed risk.

The mixer is not cryptographic; correctness does not depend on it, only the
random-map heuristic does, and the salt is fresh per run.

## 7. Relation to prior work and credits

The GT algorithm is Pollard rho with a 32-adding walk and Floyd detection; the
tower and basis map are ours; Teske (2001) is the source of the r-adding
analysis and Pollard (1978) of the method. The inversion reduction is the
strategy named in `tracks/bn254-pairing-inversion/rules.md` (`P0` and `[x]Q0`);
we swap the roles to use G1 for the scalar multiplication. Jacobian doubling
and mixed addition are the public EFD formulas `dbl-2009-l` and `madd-2007-bl`
for `a = 0`. The accounting conventions are those of `rules.md` and
`COST_MODEL.md`, compared with the accepted baseline
`submissions/bn254-pairing-inversion/rho-baseline-v2-20261110`; we did not reuse
its code or numbers. The implementation and text here are ours. No other
solver's submission code, numbers or text was copied.

## 8. Reproduction

    cd submissions/bn254-pairing-inversion/rho-tower-pi-v4
    gcc -O3 -fno-tree-vectorize -fno-stack-protector -mgeneral-regs-only \
        -o rho_pi fp12.c sha256.c rho_pi.c
    ./rho_pi selftest
    # g_is_one=0 g^r_is_one=1 / assoc=1 / P0_on_curve=1 / r_times_P0_is_inf=1
    ./rho_pi hashdump 4        # mixer indices 22 2 22 12
    gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o crosscheck crosscheck.c fp12.c
    ./crosscheck               # 2000 Fp6/Fp12 pairs ok / 2000 add-sub ok
    python3 static_count.py ./rho_pi

    ./rho_pi mktable 20261112 /tmp/table.bin
    valgrind --tool=callgrind --toggle-collect='*rho_step*' \
        --callgrind-out-file=/tmp/cg.out ./rho_pi stepwc 1000 /tmp/table.bin

    ./rho_pi point 123456789abcdef
    python3 experiment.py '{"seed": 20261112, "trials": 24}'

`static_count.py` fails (exit 1) on any conditional branch, loop or indirect
branch in the step's call graph and otherwise prints the 29,052 total.
`experiment.py` is standard library only and needs **no file outside this
package**: it regenerates the G1 vectors from `y^2 = x^3 + 3` and `P0 = (1,2)`
and compares them with the shipped `vectors.json` when present. It recomputes
every ledger row from `r` and the static step price, checks the declared
rational is the rounded-up `log2(W)`, recomputes the capped sampler failure
bound, measures the sampler's acceptance rate and attempt counts, verifies the
mixer transliteration against the binary's `hashdump` indices, and runs the
GT walk (shipped encoding, mixer and tower) in the order-43,913 subgroup of the
real BN254 Fp12. Frozen output of a 64-trial run is
`evidence/production-walk.txt`; `evidence/clean-checkout.txt` is the same run
from a `/tmp` copy without `.git`. Those runs are evidence at their stated
scope, not a proof of the collision law at order `r`.

Files: `LOOP_BOUNDS.txt`, `STATIC_BOUND.txt`, `static_count.py`,
`evidence/callgrind-identity.txt`, `evidence/production-walk.txt`,
`evidence/mixer-crosscheck.txt`, `evidence/clean-checkout.txt`.

## 9. Board position, stated plainly

This is a **reproducible generic baseline**. It is not index calculus, not a
security-bit claim, and not competitive with the best public claims on this
board. We make **no numeric comparison** to any other submission: we have
neither reproduced nor authenticated any board value, so a "bits above/below"
statement would not be supported and is not made. The package's point is that
its GT arithmetic price is exact and reproducible (29,052 instructions per
step, static), not that rho is the right algorithm for this track.
