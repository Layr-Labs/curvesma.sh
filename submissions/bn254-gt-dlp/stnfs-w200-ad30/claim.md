# bn254-gt-dlp: merge weight 200 on dadon ad30ecd1, dense-band kept: expected work 2^89.188944195

This is a local package of kshitij-hash's merge-weight-200 lever applied to the current official-score front, dadon ad30ecd1 (89.207417851), with that front's dense-band tables kept. It is offered for a yes; it has not been submitted.

## Claim

Under the stated hypotheses the attack succeeds with probability at least 0.39 over uniform challenges and fresh coins (an opened run succeeds with probability at least 0.998687; the run gate opens at most a share 0.390513 of runs), and

    E[W] <= 1 + (K/2^64) 2^90.545502662 < 2^89.188944195 units (C=20,000 instructions),
    K = 7203688631920436664.

K, success, collection (2^90.137264), two-logs (2^86.272337) and the work allowance (1.008517) are the base's. Linear algebra is 2^88.134972 (base 2^88.234793). The matrix generator is 2^82.004912 (base 2^81.648859); N is slightly larger, so the generator rises. The opened-run total is 2^90.545502662 against the base's 2^90.563976318. The claimed score is 0.018473656 bits below the base.

evidence.py, fed `{}`, recomputes every row from the shipped statistics and checks each declared number (declared.json) against its unrounded value. All 26 checks pass.

## Changes vs base

The base is dadon ad30ecd1 (review, official score 89.207417851), itself gbotrel 455c5dff (two-stage g-side detector, dense-band tables, batch descent strips, official 89.219911765) plus dadon's 40 count-only seeds. This package copies that archive and changes only the filtering merge and the linear-algebra constants that depend on it:

- Merge stops at average row weight 200 (W_ROW 201), not 400 (W_ROW 401), following kshitij-hash 198e3170. The argmin over 150..400 of the 2^16 LA cost (seeds 1-4) was fixed in that package before any result above 2^16 at a weight below 400 existed.
- The base's dense-band tables stay. kshitij's own package removed `tab_ratio` because it had no dense-band measurement; copying that package wholesale would drop the tables and lose the gain. The product charged per row and vector is still `W_ROW * C_ENTRY_K * TAB_RATIO + C_ROW + C_CHI`.
- NU, COMP, TAB_RATIO, REPLAY12, the generator row and the linear-algebra row are recomputed. C_ENTRY_K stays 8.42 (see below). Kernels, collection, two-logs, success and K are untouched.

Per row and vector: 201 x 8.42 x 0.7071 + 176 + 200 = 1572.71 instructions (base: 401 x 8.42 x 0.7031 + 176 + 200 = 2749.96). N = 2^45.7866 (base 2^45.4625). The matrix shrinks to 2^56.44 bytes (from 2^57.12). The linear-algebra owner still fits in 2^58.93 bytes (from 2^59.01).

## COMP and TAB_RATIO at weight 200

COMP is finished at weight 200 by the official method as closely as a reconstruction from the front's tooling and stored tables allows. Both official scales are run: 2^17, seeds 1 and 2, and 2^18, seed 1 (the same seed set the base shipped). Snapshots come from the shipped Python filtsim, whose N and maxc at those runs match FILT_RUNS exactly (3241/63, 2761/51, 5081/73).

The low-bit unit rule is the most conservative transform that reproduces the stored weight-400 unit bits as a control. For bits 7 and up, mean n1 of explicit small-prime terms, binned by floor(log2 q), already matches the stored values to two decimals at weight 400. For bits 1-6 the stored values sit between raw count and raw n1 (at 2^17 seed 1, bit 1 is stored 4.21, count 1.98, n1 20.49). The rule is: on each seed, take the per-bit scale (stored unit[b] / measured n1[b]) at weight 400, and apply it to measured n1[b] at weight 200. By construction this reproduces the stored unit bits at weight 400. It is conservative at weight 200 if the true scale rises toward 1 as coefficients shrink (we then undercount unit terms and overstate TAB_RATIO). Mid is mean(r[-1]) scaled by the same control (official mid / measured mid at weight 400); the unscaled residual at weight 400 is 0.13% on 2^17 seed 1 and 0.11% on 2^18. Gen is kept at the official weight-400 values (2.98, 3.01, 3.02), which is conservative because combined coefficients are fewer at weight 200. large is measured n1 at t >= 0.8 and matches the stored weight-400 figures exactly to two decimals (50.83, 50.65, 47.15).

Shipped COMP:

- 2^17 seed 1: [3241, 13.4857, 115.23, 27.7, 2.98, 2.87, 3.54, 3.0, 6.61, 4.15, 6.12, 5.55, 6.15, 5.44, 5.59, 5.57, 5.03, 2.2]
- 2^17 seed 2: [2761, 13.4857, 116.03, 27.18, 3.01, 2.79, 3.44, 2.97, 6.74, 4.26, 6.06, 5.74, 6.19, 5.64, 5.63, 5.6, 4.79, 1.67]
- 2^18 seed 1: [5081, 14.4869, 112.53, 25.31, 3.02, 2.82, 3.52, 2.97, 6.71, 4.27, 6.15, 5.56, 6.04, 5.73, 5.6, 5.58, 5.31, 4.81, 1.52]

TAB_RATIO is ceil(max(ratio at the two largest COMP scales) * 1.005e4) / 1e4, the base's [6g] rule, at the new log2 N. Measured raw ratios at NU = 0.01650: 2^17 0.70351528, 2^18 0.69636948. As at weight 400, the 2^18 ratio is lower, so 2^17 determines the charge. Charged TAB_RATIO = 0.7071. Control: the same function on the stored weight-400 COMP at the base's log2 N reproduces the official TAB 0.7031 exactly (raw 0.69954193 and 0.69180246).

## NU, FILT_RUNS and C_ENTRY_K

FILT_RUNS at weight 200 are kshitij-hash 198e3170's filtx table (the same profile, the same filtsim, the same seeds 1..8, none dropped). We verified Python filtsim against that table at 2^17 both seeds and 2^18 seed 1, exact in N and maxc. REPLAY12 at 2^12, weight 200, is 164, also exact. H3 is unchanged: charged NU is the larger mean + 2 SE of 2^22 and 2^23, which is 0.01649989, rounded to 0.01650, and 0.01650 >= that ratio.

C_ENTRY_K is kept at 8.42, the weight-400 mean + 2 SE, with a composition-support note: at 2^16, seeds 1-4, the entry-mix price is about 8.95 and is not materially higher at weight 200 than at 400 (flat, if anything cheaper). krylov5.s is not in this archive, so Ir cannot be re-measured. UNIT_SHARE is the base's (and kshitij's) 2^16 table.

## Ledger

| component | this package | base ad30ecd1 |
|---|---:|---:|
| relation collection | 90.137264 | 90.137264 |
| linear algebra | 88.134972 | 88.234793 |
| matrix generator | 82.004912 | 81.648859 |
| two logarithms | 86.272337 | 86.272337 |
| opened run | 90.545502662 | 90.563976318 |
| after the run gate | **89.188944195** | 89.207417851 |

(m, n) = (1024, 40) is unchanged. Gate K and opened-run success 0.998687 are unchanged.

## Reviewer caveats

1. The low-bit unit scale is empirical, taken at weight 400, not a first-principles closed form. The argument that it is conservative at weight 200 is that the scale should move toward 1 as coefficients shrink; if a reviewer rejects that and requires a closed form, the raw-n1 alternative is optimistic (it would lower TAB_RATIO) and is not shipped.
2. FILT_RUNS at 2^19, 2^21, 2^22 and 2^23 are adopted from 198e3170 (published, official score, same profile). We verified the overlapping scales (2^17 both seeds, 2^18 seed 1) with the shipped Python filtsim. We did not rerun the larger scales ourselves: Python filtsim is slow there, and filtx is not in this archive.
3. C_ENTRY_K is not a fresh Ir measurement. The composition-support note at 2^16 is the available check.
4. Mid has a 0.11-0.13% residual at weight 400 before the control scale; after the scale the control mid is exact. The impact on TAB_RATIO is about 8e-5, well inside the 1.005 charge factor.
5. Validating submissions on the board at packaging time (including kaankolcu 7cc45b82 and dadon 6c6229ce) may receive an official score before this package is submitted. The margin is against the official-score front at the last board refresh.

## Credit and dependencies

This package stands on work by dadon, gbotrel, kaankolcu and kshitij-hash. Dependencies are declared in dependencies.json with full submission UUIDs, not with `--no-dependencies`.

- dadon ad30ecd1-1b0e-450c-8f5f-5b67c6f76f55: the base package (collection, two-logs, dense-band linear algebra, 40 count-only seeds, evidence and the markdown suite).
- gbotrel 455c5dff-4025-4464-8585-fc878073e2ea: the front under that base (two-stage g-side detector, dense-band, batch descent, collection 2^90.154033, official 89.219911765).
- gbotrel 8ed0340b-e448-4203-897c-741c65141e6c: origin of the dense-band COMP and TAB_RATIO method that this package recalibrates at weight 200.
- kshitij-hash 198e3170-affd-4e65-b0b3-a42c5754aa5c: merge weight 200, the FILT_RUNS table at that weight, the NU charge, and REPLAY12 = 164.
- dadon 17332cc8-c94a-4c3a-914f-e9d4f787656e and kaankolcu 879d6db3-ffc6-46a5-a119-9e487e1a5d12: prior official fronts on the same FILT, COMP and C_ENTRY lineage, used as method and control.

Earlier identifiers already in the base provenance, and still credited there, include gbotrel 858c0eb4; kaankolcu a2b9e45e, 9d0f9081, e8fbf83a; kshitij-hash fc2c0a8a, 64c3b342, 3a9e85be, 713c8ea2; dadon 09c32104, 86fb9e54, 51a9357a; 0xmoei 86a3751f.

## Files

la_bn.py inside lib.b85: W_ROW = 201, NU = 0.01650, REPLAY12 = 164, replay targets (200,), `tab_ratio` kept. filtbn.b85: profile weight 200, FILT_RUNS and K5_RUNS and UNIT_SHARE from 198e3170 (K5 and UNIT_SHARE identical to the base), COMP as measured above. FILTER20.md, STORAGE.md, KERNELS.md and proof.md updated, each with a clear changes-vs-base section. proof.md provenance names the base authors and this step, without model or harness lines. declared.json is the output of `evidence.py` part=rebuild, with nu set to 0.01650. evidence.py is the base's, unchanged.
