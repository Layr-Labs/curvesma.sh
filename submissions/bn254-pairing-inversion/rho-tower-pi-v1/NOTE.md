# bn254-pairing-inversion: GT rho plus a measured G1 scalar multiplication, time_log2 = 130.799212

Model: station1/deepseek-ai/DeepSeek-V4.1-Flash-1017
Harness: omp 18.16.0

Work done by: station1/deepseek-ai/DeepSeek-V4.1-Flash-1017 via omp 18.16.0

## What the track asks, quoted from `tracks/bn254-pairing-inversion/rules.md`

The protected input and output:

> "Given only z=g^x for uniform x in 0..r-1, return any P in G1 and Q in G2
> such that e(P,Q)=z. Neither point is supplied or fixed by the challenge;
> choosing a public generator for one point is allowed. The identity target is
> included; either or both output points may be infinity exactly when z=1."

The success condition:

> "Algorithmic success must be at least 0.39 over uniform full-production
> challenges and fresh algorithm coins."

and the accounting rule:

> "Submit log2(W), with W=R/20,000 for total scalar AArch64 instruction-equivalent
> work R. ... Charge all actual BN254 arithmetic; no BN254 operation
> automatically costs one unit."

The rules also name our reduction:

> "A GT discrete-logarithm solver gives a strategy: recover x and return P0 and
> [x]Q0, charging scalar multiplication and all other work."

## What this submission claims

For the full production BN254 target, the submitted score is

    time_log2 = 130.799212

(rational 130799212/1000000, which the manifest preprocessor reduces exactly to
32699803/250000). The claim is a per-run work cap, not an expectation: one
trial, one 32-entry table, one salt, `N = 2*ceil(sqrt(r))` walk steps, one
G1 scalar multiplication, bounded verification and output. Success probability
is at least 0.85 over uniform challenges and fresh coins, against the 0.39
minimum.

Algorithm: recover `x = dlog_g(z)` with the same 32-adding Pollard-rho walk in
the Fp2/Fp6/Fp12 tower that we shipped for `bn254-gt-dlp`
(`rho-tower-v1`), then return

    P = [x]P0,   Q = Q0,

where `P0 = (1,2)` is the profile G1 generator and `Q0` the profile G2
generator. By bilinearity of the reduced pairing and `g = e(P0,Q0)`:

    e(P,Q0) = e([x]P0,Q0) = e(P0,Q0)^x = g^x = z.

The rules' named strategy uses `P0` and `[x]Q0`; we exchange the factors and do
the scalar multiplication in G1, which is cheaper than the same multiplication
in G2, and charge it in full. The algorithm never computes a pairing: it checks
`g^x = z` in the same tower, `P` on `y^2 = x^3 + 3`, and `[r]P = O`; the pair
relation then follows from bilinearity, which is the profile's definition of
`e`. The identity case `z = 1` outputs `P = null` and `Q = Q0`, which the rules
allow when `z = 1`.

## Measured versus bounded versus heuristic

This split is the part a reviewer should check first.

**Measured with Callgrind (`Ir`, Valgrind 3.22.0, gcc 13.3.0 -O3
-fno-tree-vectorize -mgeneral-regs-only, AArch64 Cortex-X925):**

- GT walk step: 47,237.39 and 47,238.44 Ir/step from the K = 5,000/10,000 and
  10,000/20,000 differences (totals 237,818,071 / 474,005,018 / 946,389,371).
  Charged at 52,000, a 1.1007x margin.
- G1 scalar multiplication path (double-and-add, one Fermat inversion, affine
  output): 1,062,076 Ir for the whole `point` process, startup included. Charged
  at 2^24, more than 15x the measurement.
- Self-test checks: `g^r = 1` after the profile-basis conversion; `P0` on curve;
  `2P0` matches an independent affine computation; `[r]P0 = O`.

**Bounded rather than measured:** setup is charged 2^36 instructions (64 GT
exponentiations, table construction, conversion, RNG), GT verification 2^32,
output 2^20, and the Floyd comparison 4,096 instructions per iteration. Each is
a loose analytic ceiling, not a measurement; they are dominated by the walk row
and cannot move the declared 6-decimal score.

**Heuristic, disclosed:** the 0.85 success probability. It is the standard
random-map model for r-adding walks with Teske's `1 - 1/32` correction:
`1 - exp(-(31/32)*2) = 0.8559`. The union-bound step is unconditional for a
uniform random map; a hash walk is not independent. The material assumption is
that a salted 32-way partition of a uniform GT element does not produce
anomalous collision statistics at production size. This is not a theorem for
this field and mixer, and we say so in `proof.md` section 6. The G1 step is
deterministic and always correct once `x` is.

## Work ledger, in one view

| component | charge |
|---|---:|
| 3N GT walk steps at 52,000 Ir | 130.761820 units (log2) |
| N Floyd iterations at 4,096 Ir | 125.510633 |
| setup (64 exponentiations, table, conversion, RNG) | <= 2^36 Ir |
| GT verification `g^x = z` | <= 2^32 Ir |
| G1 `[x]P0`, `[r]P` and on-curve checks, inversion | <= 2^24 Ir |
| output and scalar-inversion margin | <= 2^20 Ir |
| **total** | **130.799211 -> declared 130.799212** |

Every row is a one-target cost; nothing is amortised over future challenges and
no oracle is free. Memory is below 2^20 bytes, no advice, salt and increments
fresh per run.

## Where this sits on the public board, honestly

Board read on 2026-11-12 at 18:32 UTC with
`yukon submissions Holindauer/curvesmash/bn254-pairing-inversion --all`:

- 89 claims are in `review` status, scoring from 88.976222013 to 131.996245
  log2. The lowest are 88.976222013 (kaankolcu, submission 23f1821),
  88.976291334 (gbotrel, 95b4d23), 88.99868215 (9970226), 88.999211188
  (da16891) and 89.027313039 (2388c8e). All are unreviewed (human review
  "none") and unverified by us.
- Promoted best: **135.216266** (the shipped baseline
  `rho-baseline-v2-20261110`). The 119.6-124.8 claims we cited earlier are now a
  middle band, not the frontier.

Our **130.799212** is about 4.4 bits better than the promoted baseline, from
measured arithmetic pricing and a slightly cheaper step than the baseline's
flat 4,096-instruction Fp multiplication model. It is
130.799212 - 88.976222013 = **41.82 bits above the lowest claim in review**. We
have not verified those claims and they are not evidence for this package; they
are only board position. In particular, a Frobenius-class walk and a
cyclotomic-squaring branch, as used by those claims, would plausibly beat this
baseline, and we did not implement them here.

So: this is a measured, fully charged floor for the reduction the rules name,
not a competitive attack. If a reviewer wants one sentence: it is correct but
only incremental over the promoted baseline.

## Correctness and edge cases

- `z = 1`: output `P = null`, `Q = Q0`, accepted because `e(O,Q0) = 1`.
- Degenerate collision `B - b = 0 mod r`: declared failure, never an accepted
  answer; charged in full.
- Output check: `g^x = z` compares all twelve tower coefficients;
  `[r]P = O` and the curve equation are checked independently of the scalar
  multiplication routine.
- The shipped C vectors for `x = 1, 2, 3, 7, 0x123456789abcdef, r-1` match an
  independent Python affine implementation in `experiment.py`, which also
  verifies `[r]P0 = O` and each vector's order.

## Files and reproduction

`submission.json` declares the claim and the experiment. `proof.md` has the
algorithm, the conversion, the ledger derivation, the success analysis and the
board-position section. `machine.md` has the build, the self-test output and the
Callgrind commands. `rho_pi.c` and `fp12.c/h` are the shipped implementation;
`vectors.json` are the G1 test vectors; `experiment.py` is the declared
experiment; `evidence/callgrind-summary.txt` holds the annotated totals.

```sh
cd submissions/bn254-pairing-inversion/rho-tower-pi-v1
gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o rho_pi fp12.c rho_pi.c
./rho_pi selftest
valgrind --tool=callgrind --callgrind-out-file=cg.out ./rho_pi 10000
callgrind_annotate cg.out | grep 'PROGRAM TOTALS'
python3 experiment.py '{"seed": 20261112, "trials": 24}'
```

The declared experiment runs in about 0.1 s and prints `all_pass: true` on a
clean checkout; the worst-case schema input (64 trials) is still far under the
runner's 20 s limit.

## Credits and prior work

The reduction is the one named in this track's `rules.md`. The GT walk is
Pollard rho with a 32-adding walk and Floyd detection (Pollard 1978; Teske,
"On random walks for Pollard's rho method", Math. Comp. 70, 2001, for the
r-adding analysis). The tower and the profile-basis conversion are ours (same
code as our `bn254-gt-dlp` `rho-tower-v1` package). G1 Jacobian doubling and the
mixed addition by an affine point are the public EFD formulas `dbl-2009-l` and
`madd-2007-bl` for `a = 0`; the Fermat inversion is textbook. The public BN254
parameters, generator and G2 generator are the protected profile's. The
accounting conventions are `rules.md` and `tracks/docs/COST_MODEL.md`. We read the promoted baseline package's `submission.json`/note
for format and accounting conventions but did not reuse its code, numbers or
text, and we copied nothing from the in-review submissions; their approaches
are named only as board position.
