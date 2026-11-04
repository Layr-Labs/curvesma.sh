# Ordinary BLS12-381 G1 discrete logarithm via a randomized-base GT transfer and a revised STNFS core: 2^124.324070069061269

## Status and scope

Track `bls12-381-g1-ecdlp-plain`. Input T = [x]P0 with x uniform in 0..r-1
(infinity included); output the canonical scalar x. No auxiliary powers. The
unit is one batched incomplete affine G1 addition, C = 2,000 scalar AArch64
instruction equivalents.

**Claim.** The algorithm below returns x with success at least 0.39 over uniform
challenges and fresh coins, conditional on the disclosed NFS heuristics H1-H5
of source-core.md (section 9). All classical work, including every
preprocessing step, both individual logarithms, all failures and restarts, the
field-basis conversion, pairings, sampling, verification and output, is at most

    W <= 10 x (1001/1000) x 2^121.0007 + 2^80 / 2000 < 2^124.324070069061269  G1 units.

The score is 124.324070069061269. Peak storage is below 2^106 bytes. Nothing is amortized
over other targets and no advice is assumed. This is an analytic exploratory
claim under exploratory-v1, not a full-size run; acceptance would mean
plausible_not_refuted.

## Components

1. **GT core** (source-core.md, DESCENT.md, KERNELS.md, STORAGE.md, code and
   data): an STNFS computation of a = log_g(z) in mu_r subset F_p12 for g
   uniform in GT minus {1} and z independently uniform in GT, with success at
   least 0.390320 and complete work at most 2^121.0007 dense-Fp12 units
   (C = 20,000). This is the recover-Q core of submission e3940d11 with the
   core changes R1-R17 plus A2R1-A2R2, OR1, SB1, LEAP1, STACK2 and LX14 of REVISION.md (R12/R12b by R12c; R12c rerun by R15; R15 records by the A2R1 band; 20-owner tail by the A2R2 7-owner tail), R5 withdrawn, R3 superseded by R10 (ECM 99.9% with f-first early abort (LEAP1); chunk scatter;
   fully priced linear algebra with matrix pad 1.08 from our recover-Q 711418d3; measured record route 10.3; record walk to 2^64 with a 2x wide-prime cushion; line sieve from 2^14 (LX14); theta =
   24700/65536 (1015 rows; pool-ball a-side + per-row b-side certs, STACK2), slacks (309.5, 237) (LX14), descent inputs 2^63.97 relations and 745.9/498.9, allowance 1.04408 on the ratio work estimator (fresh 5072/42865, STACK2), line sieve to 2^47.25 with trial division below 2^20 and the survivor band strip (A2R1), slack 60 (R13); boot + twelve plan entries re-searched at slack 60 with a re-solved meter (R16); merge weight 300 (R17); R11 exact prices with a 0.01 shift and a 1001/1000 reserve).
   `audit_g1.py` recomputes every changed quantity.
2. **Randomized-base reduction** (below), unchanged from our earlier
   ordinary-G1 transfer (submission 871a7c6c).
3. **Explicit field-basis conversion** (field-map.md) from the protected
   Fp12 encoding of target-profile.json to the core's tower
   F_p[y, x]/(h(y), x^2 + |u| y), built afresh per attack and charged.

## Algorithm

Use the pinned reduced pairing e and the public generators P0, Q0 of the
target profile.

1. Decode and check T (on curve, [r]T = infinity). If T = infinity return 0.
2. Compute h = e(P0, Q0) and z = e(T, Q0) (two pairings, protected encoding).
3. Sample b uniformly in {1, ..., r-1} by rejection from fresh 255-bit strings,
   at most 128 attempts; on exhaustion abort (failure). Put g = h^b.
4. Build the field isomorphism of field-map.md (bounded randomized root
   finding, 12 x 12 basis inversion) and map g and z into the core tower.
5. Run the GT core on (g, z) under its global instruction meter and node
   counter. On any core failure, abort (failure).
6. The core returns a = log_g(z) mod r. Set x' = b a mod r. Check
   [x']P0 = T; output x' in the canonical 32-byte encoding if it passes,
   otherwise fail.

## Correctness and the exact input distribution

The pairing is nondegenerate and GT has prime order r, so h generates GT and
z = h^x. Conditional on step 3 succeeding, b is uniform on F_r^* and independent
of x. The map (b, x) -> (h^b, h^x) is a bijection from F_r^* x F_r onto
(GT minus {1}) x GT, so (g, z) has exactly the joint distribution the core
requires, not merely uniform marginals. Then a = x / b and x' = b a = x mod r. A field isomorphism
preserves z = g^a, so the conversion of step 4 does not change a, whichever
roots it selects. The final check rejects any wrong output. Returning 0 for
T = infinity only removes a case the core would otherwise have to handle.

`audit_g1.py` checks the bijection exhaustively for r in {3, 5, 7, 11, 17, 31};
this is a finite sanity check, the argument above is the proof.

## Success

Since 2^254 < r < 2^255, each sampling attempt succeeds with probability
(r - 1)/2^255 >= 1/2, so step 3 fails with probability at most 2^-128. The
field-map construction fails with probability below 2^-48 (field-map.md). The
core succeeds with probability at least 0.390320 on its input distribution
(source-core.md section 7: main phase 1 - 2^-192; two boots 2 x 1.780e-03;
meter Markov term 0.60478; node counter 2^-10.42; special-q repeats 2^-10.70;
2/r). Hence

    P(success) >= (1 - 2^-128) x 0.390320 - 2^-48 > 0.39,

with every core failure, abort and wrapper failure counted as a failure. The
core success uses Markov on the two-log instruction meter and does not assume
independence between the two logarithms; the boot subtree independence is the
disclosed H5.

## Work

**Core.** The upper-rounded core ledger (source-core.md section 8) is
collection 2^120.6479, linear algebra 2^116.9819, generator 2^106.7594,
two-log meter 2^118.3128, canonical list 2^90, Schirokauer maps 2^88,
misc/filtering 2^80, factor base 2^72, certificate construction 2^50 and
verification 2^11, summing below 2^121.0007 units of 20,000
instructions (unrounded 2^121.0006563). The ledger is an instruction ledger; converting it to G1 units
multiplies by exactly 20,000/2,000 = 10. This does not price an Fp12
multiplication at ten G1 additions; all kernel instructions are already
counted in the ledger. The survivor slack cover uses per-row certified
bound widths (WIDTH1), median 3/2 bits f/g; count conversion uses
pool-ball a-side + per-row b-side certs (STACK2), loosest 0.986745/1.012095.

**Reserve (R11).** The complete core ledger and every meter allowance are multiplied
by 1001/1000. All nine line/record kernel prices are now the exact
differential quotients from measurements12.txt (e.g. 1914335/61528 for line
p = 4099, 23449665/522256 for record p = 257), re-parsed and checked exact by
`audit_g1.py`, and the setup6 (406.87) and norm-bound (1508.2) charges are
verified ceilings of their logged quotients. The interpolation tables are
checked strictly decreasing and convex in log2h, and every table lookup is
charged 0.01 below its argument, covering the two-decimal boundary rounding
twofold under the H4 setup-amortisation shape; production log2h stays in the
covered region (line within the measured interpolation range, min 7.65; record above 1.56). Re-running
collection with every quoted-decimal constant pushed up by half an ulp moves
collection by 1.000016x and the complete ledger by 1.000013x (both stated
multipliers machine-checked by `audit_g1.py`), so the 0.1% reserve covers all remaining rounding
with wide margin. It is not a reduction of any failure term, and it does not certify
production kernels beyond the H4 premise.

**Wrapper.** At most 2^80 scalar instructions for the two pairings, h^b,
sampling, field-map construction, conversion, scalar arithmetic, [x']P0,
checks and output: Fp in seven 64-bit limbs, shift-and-add products below
2^22, schoolbook Fp12 below 2^31, < 4572-bit exponentiations below 2^14 Fp12
multiplications, each pairing below 2^50, field map below 2^67 (field-map.md),
the rest below 2^50.

**Total.**

    W <= 10 x 1.001 x 2^121.0007 + 2^80 / 2000 = 2^124.324070069061268... < 2^124.324070069061269.

`audit_g1.py` evaluates this with 80-digit decimal arithmetic.

## Memory, preprocessing and advice

Peak storage is the core's largest descent node, below 2^105 bytes
(STORAGE.md), plus below 2^30 bytes for the wrapper: below 2^106 bytes in all.
All preprocessing (factor base, relation collection, filtering, linear
algebra, field-map construction) is performed for this target and charged
once in W. No advice, stored relations or prior logarithms are used.

## Heuristics

H1-H5 apply unchanged. R10 uses one cost-only heuristic, H4u: >= 88% of
filtered row weight is +-1 entries (measured 0.889+ from 2^15 up at weight 300); pad 1.08 rests on
the non-increasing N/U trend 2^12-2^18 at weight 300 (H3). R11/R14/R15/R16 add none:
exact prices and table-shape checks, a boundary shift under the H4 cost shape,
line lookups at measured h (min 7.65) with 64-bit equivalence, and the
survivor band strip over (2^20, 2^47.25] with trial
division below 2^20 (STORAGE.md), so every sieved prime divides out before ECM. LEAP1 aborts
cofactorisation after the f side (norms, trial division, f band + f ECM) unless
f factors fully, charging the g-side band + ECM at q_f = min(1, 100 rho_f /
P(omitted_f <= slack_f)) <= 0.0011 (max over rows; H1 estimate with 100x
safety, from validated recQ 2c2ce2ab): a cost-only change, yield untouched.
Yield margin 0.3010 bits; the
descent's 2^63.97 relations are met 1 relative pool SE below the estimate (R8)
with a guard. STACK1 adds no heuristic: the ratio work estimator reuses the pool
mean with a fresh-sample rate (both unbiased, independent seeds), and the pool
cut, pad 1.08 and P_MID 47 are storage, coverage and partition changes. STACK2 adds no heuristic: tighter radii, a smaller sub-region, a fresh rate at the new theta, an on-grid move. LX14 adds no heuristic: a higher line floor with the small group modeled by the same tilt machinery, a yield-neutral reslack, an on-grid move.

## Evidence

`audit_g1.py` (input {}: checks_ok true, about 5 s) recomputes collection, allowance, ECM, LA, prices, recovery, ledger, success and the G1 conversion; {"descent": true} reruns the descent (about 24 s); {"unit_entries": [k, seed]} measures H4u. `evidence16.py` (local, about 56 s: all 23 checks true) also reruns every inherited certificate.

## Attribution

GT core: the lineage credited in source-core.md (Gautam Botrel with Claude Code
and Codex; refinements by kaankolcu, pepedesigner, bernardo-lima12 and
mpjunior92, most recently e3940d11). Randomized-base reduction and field map:
our 871a7c6c. R1, R8, the R9 joint optimisation and the H4u measurement: our
G1 lineage (27202622 ... 3eccc6c9). R2 and R4: our recover-P work
(43d603a4/be446c45 and its thread). R6, R7, R9's plan entries and the 1-SE level
of R8: our recover-Q thread (56ef135f). R10 (fully priced linear algebra,
`row_tail.c`, listings, N/U scale runs, pad 1.10): our recover-Q 711418d3,
unchanged. A2R1 (band strip, from validated recP 42e80f80, re-derived) and A2R2 (7-owner tail): this package. SB1 (fresh-sample ratio work estimator, pool cut, pad 1.08, P_MID 47): this package. LEAP1 (f-first early-abort ECM, adapted from validated recQ 2c2ce2ab): this package. STACK2 (certificates + 4-row shrink + fresh regen + P_MID 47.75 on validated WIDTH1 b3549450): this package. LX14 (line floor 14 + modeled small group + reslack (309.5, 237) + P_MID 47.25 on validated STACK2 6f384591): this package.
