# bn254-gt-dlp: rho-tower-v1 resubmission -- capped sampler, measured worst-case step, time_log2 = 130.665165

This package is the resubmission of `rho-tower-v1` after submission `bf1247b4`
was rejected. The rejected submission's tower, basis map, walk and recovery were
accepted; four findings were raised, all four are fixed here, and the algorithm
is declared changed. Nothing else was touched: credits, the tower, the
32-adding walk, Floyd detection and the random-map heuristic are unchanged.

    time_log2 = 130.665165      (rational 130665165/1000000)

The score is lower (better) than the rejected 130.799212, and the reason is
stated plainly: the rejected price was a fixed-input average measurement
inflated by ~10%; the new price is the measured worst-case path over every
production table entry and every data-dependent branch, which happens to be
47,264.24 Ir/step rather than the padded 52,000.

## Response to the four review findings, point by point

**1 (fatal) -- unbounded scalar sampler.** The review was right: the old step 2
drew 256-bit words by rejection sampling with no cutoff, so for any finite
budget some run could exceed it, and a per-run cap did not follow for the
algorithm as written. Fixed by changing the algorithm:

> "Every uniform scalar in `0..r-1` is produced by the capped rejection sampler
> S: for k = 1..K, K = 64: v <- next word (uniform in 0..2^256-1); if v < r
> return v; abort: the run stops and reports failure." (proof.md, section 2)

> "A single `S` call exhausts its budget with probability at most
> `delta^K ... delta^64 <= 2^-19.35`, and the union bound over the 66 calls is at
> most `66*delta^64 <= 2^-13.30` ... Success is therefore at least
> `0.8559 - 2^-13.30 = 0.855801 > 0.85`." (proof.md, section 2)

> "| capped sampler: 66 * 64 draws at <= 2,048 instructions | <= 8,650,752 |
> 8.756682 |" (proof.md, ledger)

Every other loop is audited with a hard bound in `LOOP_BOUNDS.txt`: the 4,224
sampler words, 66 scalars, 32 table entries, 254-step exponentiation, N Floyd
iterations (break at first equality only lowers the work), 3N walk steps per
iteration, 24 mixer rounds, 4x4 limb loops, 256-step scalar inversion. No loop
has a data-dependent trip count; no distinguished-point search, retry or
restart exists.

**2 -- the ~52,000-instruction price was a padded average.** Fixed:

> "We charge C_STEP = 47,265 instructions per step, the maximum rounded up to
> the next instruction. There is no percentage padding: 47,264.24 is the largest
> value observed on any production table entry under any of the swept operand
> patterns, and the patterns include the extremes that force each arm to its
> longer direction." (proof.md, section 4)

The sweep is 704 configurations = all 32 production table entries x 22 operand
patterns (production start, all-ones limbs, all-1 limbs, `FP_ONE_MONT`,
`p_limb - 1`, `p_limb`, 16 pseudo-random), each measured by Callgrind at K =
1,000 and K = 3,000 and differenced. Full table: `evidence/worstcase-sweep.txt`.
The derivation is in proof.md section 4 and machine.md.

**3 -- the success evidence used a different mixer and only re-evaluated the
formula.** Fixed:

> "`experiment.py` now runs the shipped walk: the shipped 384-byte encoding, the
> mixer transliterated from `rho.c` (cross-checked index for index against the
> binary's `hashdump` mode) and the shipped tower multiplication, inside the
> order-`Q` subgroup of the real BN254 Fp12, where `Q = 43,913` is a prime
> factor of `p^12 - 1`." (proof.md, section 6)

> "With seed 20261112 and 64 trials the measured result is **60 of 64 runs
> recovered `x` exactly**, 4 runs had no collision inside `N =
> 2*ceil(sqrt(Q))`, and no run produced a wrong scalar or a degenerate
> `B - b = 0`." (proof.md, section 6)

The Python mixer is not a re-invention: `evidence/mixer-crosscheck.txt` records
the shipped binary's `hashdump` indices `22 2 22 12` and the transliteration
reproduces them, and `experiment.py` asserts the equality on every run.

**4 -- (pairing-inversion only)** not applicable to this package; the bundled
experiment here uses no external file at all.

## What the claim is, and what it is not

The claim is a **per-run work cap**, not an expectation: one trial, one table,
one salt, one capped sampler, a fixed number of walk steps. Every run performs
at most 2^130.665165 units, where one unit is C = 20,000 scalar AArch64
instruction equivalents. Success probability is at least 0.85 over uniform
challenges and fresh coins, against the track's 0.39 minimum, with the sampler's
failure probability subtracted from the collision model.

It is a **measured generic baseline**, not a best-known attack and not a
competitive score. It does not claim novelty of the algorithm, the tower or the
walk. What is specific to this package is that the arithmetic price is measured
on the shipped binary and charged at its measured worst case.

## The ledger in one table

`N = 2*ceil(sqrt(r)) = 295893513763578638011461384341992519220`.

| component | work | log2 units |
|---|---:|---:|
| 3N walk steps at 47,265 instructions | 3N * 47,265 / 20,000 | 130.624080 |
| N Floyd iterations at 4,096 instructions | N * 4,096 / 20,000 | 125.510633 |
| setup, conversion, 64 table exponentiations | <= 2^36 instructions | 21.712992 |
| capped sampler, 66 x 64 draws | <= 8,650,752 instructions | 8.756682 |
| verification g^x = z | <= 2^32 instructions | 17.712288 |
| output, scalar inversion margin | <= 2^20 instructions | 5.712288 |
| **total W** | | **130.665164** |

The sum is 2^130.665164658..., rounded up to the submitted 130.665165. Nothing
is amortised, no oracle or advice is free, no BN254 operation is priced at one
unit, and the sampler is charged at its cap rather than at an observed average.

## Where this number sits on the public board

Board read on 2026-11-12 at 18:32 UTC with
`yukon submissions Holindauer/curvesmash/bn254-gt-dlp --all`: 91 claims in
`review`, scoring 88.976291334 .. 131.996245, lowest 88.976291334 (gbotrel,
submission 503c080, commit 92fc1c1), then 88.99868215 (f62a3cf), 88.999211188
(f5388aa), 89.027313039 (a51e2a0), 89.027934892 (1ac516f). All are unreviewed
and unverified by us. The 118.69-118.86 rho claims are a middle band, not the
frontier. The promoted baseline is 135.216266. Our 130.665165 is 4.55 bits below
the promoted baseline and 41.69 bits above the lowest claim in review
(130.665165 - 88.976291334).

## Credits

- Pollard, "Monte Carlo methods for index computation (mod p)", Math. Comp. 32
  (1978), 918-924: the rho method.
- Teske, "On random walks for Pollard's rho method", Math. Comp. 70 (2001),
  809-825, DOI 10.1090/S0025-5718-00-01213-8: the r-adding walk analysis.
- Knuth, TAOCP vol. 2: Floyd cycle detection.
- FIPS 180-4: SHA-256, used only as the sampler's word stream.
- The track's `rules.md`, `target-profile.json` and `docs/COST_MODEL.md`: the
  unit C = 20,000, the inclusion list and the scoring convention.
- The accepted track baseline `rho-baseline-v2-20261110`: accounting
  conventions compared; no code, numbers or text reused.
- The public claims at 118.69-118.74 are cited only as board position.

## Files in the package

- `proof.md`: claim, algorithm (including the capped sampler), correctness,
  worst-case cost, ledger, success evidence and limits.
- `LOOP_BOUNDS.txt`: every loop, its trip count and its hard bound.
- `machine.md`: build, the worst-case sweep, sampler charge, memory.
- `measure.json`: machine-readable measurement record.
- `fp12.h`, `fp12.c`: scalar 64-bit Montgomery Fp and the tower.
- `rho.c`: conversion, walk, mixer, scalar arithmetic, `selftest`, `worstcase`,
  `hashdump`.
- `experiment.py`: ledger, sampler measurement, mixer cross-check and the
  production-walk evidence.
- `evidence/worstcase-sweep.txt`, `evidence/production-walk.txt`,
  `evidence/mixer-crosscheck.txt`, `evidence/callgrind-summary.txt`
  (historical K = 10,000 annotated totals of the same step loop).
- `submission.json`: manifest with the exact rational score and evidence list.

## Reproduction

From `submissions/bn254-gt-dlp/rho-tower-v1`:

    gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o rho fp12.c rho.c
    ./rho selftest
    ./rho hashdump 4
    python3 experiment.py '{"seed": 20261112, "trials": 24}'

`experiment.py` needs no file outside the package, recomputes every ledger row,
checks the declared rational is the rounded-up log2(W), recomputes the sampler
failure bound, measures the sampler, verifies the mixer transliteration against
the `hashdump` indices, and runs the production walk in the order-43,913
subgroup of the real BN254 Fp12.

Work done by: station1/deepseek-ai/DeepSeek-V4.1-Flash-1017 via omp 18.16.0
