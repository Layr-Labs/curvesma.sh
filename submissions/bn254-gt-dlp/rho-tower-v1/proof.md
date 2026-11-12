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

    W <= 2^130.799211 < 2^130.799212   units,

where one unit is 20,000 scalar AArch64 instructions. This is a per-run
deterministic work cap, not an expectation: the loop always stops within
`N = 2*ceil(sqrt(r))` iterations, and the table build, recovery, verification
and output work are fixed. The score is `log2(W)` by the track rules. Success
is at least the protected 0.39.

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
2. Draw a fresh 64-bit salt and 64 independent scalars `c_j, d_j` uniform in
   `0..r-1`, by rejection sampling 256-bit words (accept a draw below `r`).
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
the step loop at three lengths:

| steps K | total Ir | Ir per step |
|---:|---:|---:|
| 5,000 | 237,821,066 | |
| 10,000 | 474,013,018 | 47,238.39 |
| 20,000 | 946,407,371 | 47,239.44 |

The two differences give 47,238.39 and 47,239.44; we charge
**52,000 instructions per step**, a factor 1.1007 over the larger measured
value, to cover compiler scheduling, cache/address variation and the transition
between the measured loop and the full run. The raw Callgrind output is
`evidence/callgrind-summary.txt`. `machine.md` gives the reproduction commands.

The measured loop includes the hash, the two scalar additions, the table load
and the multiplication; it excludes only the Floyd comparison and loop control,
which are charged separately below at 4,096 instructions per iteration (twelve
word comparisons plus a counter and branch; the bound is far above the
measured cost of those operations).

## 5. Complete work ledger

`r` is the profile subgroup order, `N = 2*ceil(sqrt(r)) =
295893513763578638011461384341992519220`.

| component | work | log2 units (C=20,000) |
|---|---:|---:|
| 3N walk steps at 52,000 instructions | 3N * 52,000 / 20,000 | 130.761820 |
| N Floyd iterations, 4,096 instructions each | N * 4,096 / 20,000 | 125.510633 |
| setup: 64 exponentiations, table, conversion, RNG | <= 2^36 instructions | 21.712288 |
| verification `g^x = z`, canonical output | <= 2^32 instructions | 17.712288 |
| output and final scalar inversion margin | <= 2^20 instructions | 5.712288 |
| **total W** | | **130.799211** |

The setup row: 64 exponentiations with 254-bit exponents, each at most 508
`Fp12` multiplications, is 32,512 multiplications; even at the loose analytic
bound of 2^20 instructions per multiplication used by the earlier published
BN254 baseline, that is below 2^36. We charge 2^36, about 30x the measured
cost of that many multiplications. The verification row is two exponentiations
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

We claim the smaller, conservative **0.85**. This leaves a large margin over
the protected 0.39.

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
    ./rho selftest          # prints g_is_one=0 g^r_is_one=1
    valgrind --tool=callgrind --callgrind-out-file=cg.out ./rho 10000
    callgrind_annotate cg.out | grep 'PROGRAM TOTALS'
    python3 experiment.py '{"seed": 20261112, "trials": 24}'

`experiment.py` recomputes every row of the ledger from `r`, checks the
declared rational against the rounded-up `log2(W)`, checks the 0.85 success
lower bound, and runs the same 32-adding walk on small prime groups as a
finite correctness/statistics check. Those small runs are evidence at their
actual scope only; they do not prove production collision statistics.
