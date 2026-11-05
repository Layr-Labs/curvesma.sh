# Ordinary BLS12-381 G1 discrete logarithm via a randomized-base GT transfer and a revised STNFS core: 2^124.278670069061269

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

    W <= 10 x (1001/1000) x 2^120.9553 + 2^80 / 2000 < 2^124.278670069061269  G1 units.

The score is 124.278670069061269. Peak storage is below 2^106 bytes. Nothing is amortized
over other targets and no advice is assumed. This is an analytic exploratory
claim under exploratory-v1, not a full-size run; acceptance would mean
plausible_not_refuted.

## Components

1. **GT core** (source-core.md, DESCENT.md, KERNELS.md, STORAGE.md, code and
   data): an STNFS computation of a = log_g(z) in mu_r subset F_p12 for g
   uniform in GT minus {1} and z independently uniform in GT, with success at
   least 0.390320 and complete work at most 2^120.9553 dense-Fp12 units
   (C = 20,000). This is the recover-Q core of submission e3940d11 with the
   core changes of REVISION.md (R1-R17, A2R1-A2R2, OR1, SB1, LEAP1, STACK1, WIDTH1, STACK2, WIDE1; R5 withdrawn, R3 superseded by R10):
   ECM 99.9% with f-first early abort; chunk scatter; fully priced linear algebra with matrix pad 1.08 (our recover-Q 711418d3);
   record walk to 2^64 at measured wide-prime prices (WIDE1); line sieve from 2^11 to 2^47.75 with trial division below 2^20 and the
   survivor band strip; theta = 24700/65536 (1015 rows; pool-ball a-side + per-row b-side certs); slacks (297, 233); descent inputs
   2^63.97 relations and 745.9/498.9; allowance 1.04432 on the ratio work estimator (fresh 5072/42865); slack 60; boot + twelve plan
   entries re-searched (R16); merge weight 300; exact prices with a 0.01 shift and a 1001/1000 reserve.
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
collection 2^120.5896, linear algebra 2^116.9819, generator 2^106.7594,
two-log meter 2^118.3128, canonical list 2^90, Schirokauer maps 2^88,
misc/filtering 2^80, factor base 2^72, wide-prime lattices 2^70.6, certificate construction 2^50 and
verification 2^11, summing below 2^120.9553 units of 20,000
instructions (unrounded 2^120.9551758). The ledger is an instruction ledger; converting it to G1 units
multiplies by exactly 20,000/2,000 = 10. This does not price an Fp12
multiplication at ten G1 additions; all kernel instructions are already
counted in the ledger. The survivor slack cover uses per-row certified
bound widths (WIDTH1), median 3/2 bits f/g; count conversion uses
pool-ball a-side + per-row b-side certs (STACK2), loosest 0.986745/1.012095.

**Reserve (R11).** The complete core ledger and every meter allowance are multiplied
by 1001/1000. All line/record kernel prices are the exact differential quotients of
measurements12.txt (e.g. 23449665/522256 for record p = 257), re-parsed by `audit_g1.py`;
setup6 (406.87) and norm-bound (1508.2) charges are verified ceilings; the tables are checked
decreasing and convex in log2h and every lookup is charged 0.01 below its argument, covering
the two-decimal boundary rounding twofold (production log2h: line min 7.15, record above 1.56).
Every quoted-decimal constant pushed up by half an ulp moves collection by 1.000016x and the
ledger by 1.000013x (machine-checked), well inside the 0.1% reserve. The reserve reduces no
failure term and certifies no kernel beyond H4.

**Wrapper.** At most 2^80 scalar instructions for the two pairings, h^b,
sampling, field-map construction, conversion, scalar arithmetic, [x']P0,
checks and output: Fp in seven 64-bit limbs, shift-and-add products below
2^22, schoolbook Fp12 below 2^31, < 4572-bit exponentiations below 2^14 Fp12
multiplications, each pairing below 2^50, field map below 2^67 (field-map.md),
the rest below 2^50.

**Total.**

    W <= 10 x 1.001 x 2^120.9553 + 2^80 / 2000 = 2^124.278670069061268... < 2^124.278670069061269.

`audit_g1.py` evaluates this with 80-digit decimal arithmetic.

## Memory, preprocessing and advice

Peak storage is the core's largest descent node, below 2^105 bytes
(STORAGE.md), plus below 2^30 bytes for the wrapper: below 2^106 bytes in all.
All preprocessing (factor base, relation collection, filtering, linear
algebra, field-map construction) is performed for this target and charged
once in W. No advice, stored relations or prior logarithms are used.

## Heuristics

H1-H5 apply unchanged. Cost-only: H4u (>= 88% of filtered row weight is +-1 entries;
measured 0.889+ from 2^15 up at weight 300) and the non-increasing N/U trend behind pad 1.08
(H3). Exact prices, the lookup shift, line lookups at measured h, trial division below 2^20 and
the band strip (STORAGE.md) add none. LEAP1 runs the g side only if f factors fully, charged at
q_f = min(1, 100 rho_f / P(omitted_f <= slack_f)) <= 0.0017 (H1 estimate, 100x safety): yield
untouched. The descent's 2^63.97 relations sit 1 relative pool SE below the estimate (R8) with a
guard. STACK1/STACK2 add none (unbiased ratio estimator on independent seeds; tighter radii, a
smaller sub-region, a fresh rate, partition moves). WIDE1 adds none: measured prices of the
same kernels replace the R6 2x cushion under H4 (KERNELS.md).

## Evidence

`audit_g1.py` (input {}: checks_ok true, about 5 s) recomputes collection, allowance, ECM, LA, prices, recovery, ledger, success and the G1 conversion; {"descent": true} reruns the descent (about 24 s); {"unit_entries": [k, seed]} measures H4u. `evidence16.py` (local, about 56 s: all 23 checks true) also reruns every inherited certificate.

## Attribution

GT core: the lineage credited in source-core.md (Gautam Botrel with Claude Code
and Codex; refinements by kaankolcu, pepedesigner, bernardo-lima12 and mpjunior92, most
recently e3940d11). Randomized-base reduction and field map: our 871a7c6c. R1, R8, R9 and H4u:
our G1 lineage (27202622 ... 3eccc6c9). R2, R4: our recover-P work (43d603a4/be446c45). R6, R7,
R9's plan entries and the 1-SE level of R8: our recover-Q thread (56ef135f). R10: our recover-Q
711418d3. A2R1 (band strip from validated recP 42e80f80), A2R2, SB1, LEAP1 (from validated recQ
2c2ce2ab), STACK1, WIDTH1, STACK2 (validated 6f384591, this package's file base): ours. WIDE1
(`prepare_pair_w`/`setup6w` appended to gbotrel's `slice12.c`, curvesma.sh 0bbca0b2): this package.
