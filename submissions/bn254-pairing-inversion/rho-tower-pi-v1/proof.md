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

    W <= 2^130.799211 < 2^130.799212   units.

The score is `130799212/1000000`, which the manifest preprocessor reduces
exactly to `32699803/250000`. Success is far above the protected minimum 0.39.
The claim is a per-run deterministic work cap: one trial, one table, one salt,
a fixed number of walk steps, one scalar multiplication, and bounded
verification and output.

This is a **baseline, not a best-known attack**, and it is not an improvement
over the best public claims on this board (see section 9). Its purpose is to
give the inversion track the same measured floor that `bn254-gt-dlp` has: the
per-step price below is the Valgrind instruction count of the shipped scalar C
implementation, and the G1 part is measured and cross-checked against an
independent Python implementation.

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
2. Convert `g` and `z` to the tower (section 1), and draw a fresh 64-bit salt.
3. Draw 64 independent scalars `c_j, d_j` uniform in `0..r-1` by rejection
   sampling; build the 32-entry table `M_j = g^{c_j} z^{d_j}` with both scalar
   labels. Two fixed-base/fixed-exponent exponentiations per entry.
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

For the G1 part, standard Jacobian formulas with `a = 0`: the mixed addition is
`madd-2007-bl` for the affine base, doubling is `dbl-2009-l`, and the
double-and-add loop computes `[x]P0` for the 256-bit scalar (the loop starts at
bit 255; `x < r < 2^254`, so leading zero bits only cost doublings of the point
at infinity). The affine conversion uses one Fermat inversion `Z^{p-2}`. The
`[r]P = O` and on-curve checks are independent of the scalar multiplication
routine. On the shipped vectors in `vectors.json`, the C output equals an
independent Python affine implementation for `x = 1, 2, 3, 7,
0x123456789abcdef` and `r-1`.

## 4. Arithmetic and measured step cost

All arithmetic is scalar 64-bit. `Fp` elements are four little-endian words in
Montgomery form; `fp_mul` is a 4-limb CIOS Montgomery multiplication with a
six-word accumulator (`fp12.c`). The tower multiplication is `Fp2` Karatsuba
(three `Fp` multiplications), `Fp6` with six `Fp2` products plus three
multiplications by `xi = 9+u` (each four `Fp` additions), and `Fp12` with three
`Fp6` products and two multiplications by `v`.

One GT step performs one `Fp12` multiplication, one hash of 48 words, two 4-word
additions modulo `r`, one table record load and loop bookkeeping. The binary
`rho_pi` was compiled with

    gcc 13.3.0 -O3 -fno-tree-vectorize -mgeneral-regs-only

on AArch64 (Cortex-X925, Linux 6.8.0), the machine class of the cost model.
Valgrind 3.22.0 Callgrind (`Ir`) on the step loop:

| steps K | total Ir | difference Ir/step |
|---:|---:|---:|
| 5,000 | 237,818,071 | |
| 10,000 | 474,005,018 | 47,237.39 |
| 20,000 | 946,389,371 | 47,238.44 |

We charge **52,000 instructions per step**, a factor 1.1007 over the larger
difference, to cover compiler scheduling, cache/address variation and the
transition between measured loop and full run. The annotated raw totals are in
`evidence/callgrind-summary.txt`; `machine.md` has the reproduction commands.

The GT step loop includes the hash, the two scalar additions, the table load
and the multiplication. It excludes the Floyd comparison and loop control,
charged separately at 4,096 instructions per iteration.

The G1 scalar multiplication was measured separately by Callgrind on the
shipped `point` path with a 254-bit scalar: `1,062,076 Ir` for the whole process
including the profile-generator conversion and output formatting. The
double-and-add itself is therefore well below `2^21` instructions; the ledger
charges `2^24` for the scalar multiplication, the `[r]P` and on-curve checks,
and the affine inversion, a margin of more than 15x over the measured process.

## 5. Complete work ledger

`r` is the profile subgroup order, `N = 2*ceil(sqrt(r)) =
295893513763578638011461384341992519220`.

| component | work | log2 units (C=20,000) |
|---|---:|---:|
| 3N GT walk steps at 52,000 instructions | 3N * 52,000 / 20,000 | 130.761820 |
| N Floyd iterations, 4,096 instructions each | N * 4,096 / 20,000 | 125.510633 |
| setup: 64 GT exponentiations, table, conversion, RNG | <= 2^36 instructions | 21.712288 |
| GT verification `g^x = z`, canonical output | <= 2^32 instructions | 17.712288 |
| G1 `[x]P0`, `[r]P` and on-curve checks, affine inversion | <= 2^24 instructions | 13.712288 |
| output and final scalar inversion margin | <= 2^20 instructions | 5.712288 |
| **total W** | | **130.799211** |

The total is dominated by the walk row: adding the G1 row changes `log2(W)` by
less than `2^-79`, so the declared rational `130799212/1000000` is unchanged
from the GT-only bound and remains a rounded-up cap. Setup, verification and
output rows use the same conservative instruction charges as our GT package
(2^36 for 64 exponentiations, 2^32 for verification, 2^20 for output), each
far above the measured cost. No row is amortised over future targets, no oracle
is free, and no BN254 operation is priced at one unit.

Memory: 32 GT table records of 448 bytes, two GT walk states, the hash scratch,
G1 temporaries and code, below 2^20 bytes. No advice; the salt and increments
are fresh per run and charged in setup.

## 6. Success probability: heuristic, support and limits

The only probabilistic part is the GT rho collision. Model the walk as a random
map on the `r` points of GT with the salt and 64 increments as fresh coins. For
`t` distinct transitions the no-repeat probability is at most
`exp(-t(t-1)/(2r))`; at `t = N = 2 ceil(sqrt(r))` this is `exp(-2) = 0.1353`, so
a repeat occurs with probability at least `0.8647`. The conventional r-adding
correction multiplies the collision rate by `1-1/32`, giving
`1 - exp(-(31/32)*2) = 0.8559`; we claim the smaller **0.85**.

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

    cd submissions/bn254-pairing-inversion/rho-tower-pi-v1
    gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o rho_pi fp12.c rho_pi.c
    ./rho_pi selftest
    # g_is_one=0 g^r_is_one=1 / assoc=1 / P0_on_curve=1 / r_times_P0_is_inf=1
    valgrind --tool=callgrind --callgrind-out-file=cg.out ./rho_pi 10000
    callgrind_annotate cg.out | grep 'PROGRAM TOTALS'
    valgrind --tool=callgrind --callgrind-out-file=cg.point.out \
      ./rho_pi point 8000000000000000000000000000000000000000000000000000000000000000
    ./rho_pi point 123456789abcdef
    python3 experiment.py '{"seed": 20261112, "trials": 24}'

`experiment.py` (standard library only) recomputes every ledger row from `r`,
checks the declared rational is the rounded-up `log2(W)`, checks the success
bound, recomputes the G1 vectors in `vectors.json` with an independent affine
implementation and the exact scalar multiplication of this document, checks
`[r]P0 = O`, and runs the same 32-adding walk on small prime additive groups as
scope-limited correctness/statistics evidence.

## 9. Board position, stated plainly

At preparation time the public `bn254-pairing-inversion` board shows the
promoted baseline `135.216266` and several unreviewed `review`-status claims at
119.6-124.76. Our `130.799212` is about 4.4 bits below the promoted baseline
(an arithmetic improvement, not an algorithmic one) and roughly 11 bits above
the best public claim in review. We have not verified those claims; they are
listed here as position, not as support. This package is a measured, fully
charged baseline of the reduction the rules name, suitable as a floor for
better inversion methods.
