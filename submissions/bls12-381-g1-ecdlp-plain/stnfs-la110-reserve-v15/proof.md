# Ordinary BLS12-381 G1 discrete logarithm via a randomized-base GT transfer and a revised STNFS core: 2^124.84318

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

    W <= 10 x (1001/1000) x 2^121.5198 + 2^80 / 2000 < 2^124.84318  G1 units.

The score is 124.84318. Peak storage is below 2^106 bytes. Nothing is amortized
over other targets and no advice is assumed. Analytic exploratory claim, not a full-size run.

## Components

1. **GT core** (source-core.md, DESCENT.md, KERNELS.md, STORAGE.md, code and
   data): an STNFS computation of a = log_g(z) in mu_r subset F_p12 for g
   uniform in GT minus {1} and z independently uniform in GT, with success at
   least 0.390320 and complete work at most 2^121.5198 dense-Fp12 units
   (C = 20,000): the recover-Q core of e3940d11 with REVISION.md R1-R10
   (R5 withdrawn, R3 replaced by R10).
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
requires (uniform nonidentity base, independently uniform target), not merely
uniform marginals. Then a = x / b and x' = b a = x mod r. A field isomorphism
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
(source-core.md section 7: main phase 1 - 2^-192; two boots 2 x 7.746e-04;
meter Markov term 0.60593; node counter 2^-9.29; special-q repeats 2^-10.70;
2/r). Hence

    P(success) >= (1 - 2^-128) x 0.390320 - 2^-48 > 0.39,

with every core failure, abort and wrapper failure counted as a failure. The
core success uses Markov on the two-log instruction meter and does not assume
independence between the two logarithms; the boot subtree independence is the
disclosed H5.

## Work

**Core.** The upper-rounded core ledger (source-core.md section 8) is
collection 2^121.1117, linear algebra 2^117.7327, generator 2^107.2082,
two-log meter 2^118.996, canonical list 2^90, Schirokauer maps 2^88,
misc/filtering 2^80, factor base 2^72, certificate construction 2^50 and
verification 2^11, summing to 2^121.5197816 < 2^121.5198 units of 20,000
instructions. The ledger is an instruction ledger; converting it to G1 units
multiplies by exactly 20,000/2,000 = 10. This does not price an Fp12
multiplication at ten G1 additions; all kernel instructions are already
counted in the ledger.

**Reserve.** The complete core ledger and every meter allowance are multiplied
by 1001/1000 (earlier packages: 101/100) to cover the rounding of measured kernel
prices. `audit_g1.py` (`rounding_ratios`, check `reserve_covers_rounding`)
compares each price with the exact Ir ratio of its differential row in the
measurement log (12 rows, at most 1.000104: 1914335/61528 against 31.11) or, for
prices logged only as rounded values (apply 7.00, recovery 7.003, layout 12.325,
zero 0.625, scan 9.313, Montgomery 367.76), with the top of the rounding
interval (at most 1.0008, zero). Every measured price is within 0.08% of its
charge. The reserve is not a reduction of any failure term, and it does not
certify production kernels beyond the H4 premise.

**Wrapper.** At most 2^80 scalar instructions in total for the two pairings,
h^b, the rejection sampling, the field-map construction (randomized
factorization over F_p12, 12 x 12 inversion, verification), input conversion,
scalar arithmetic mod r, the G1 scalar multiplication [x']P0, subgroup and
encoding checks and output. The bound uses elementary scalar integer
implementations only: Fp in seven 64-bit limbs, a 381-step shift-and-add
modular product below 2^22 instructions, schoolbook Fp12 products with
reduction below 2^31 instructions, exponentiations of fewer than 4572 bits
below 2^14 Fp12 multiplications, each pairing below 2^50 instructions, the
field-map construction below 2^67 instructions (field-map.md), and the
remaining scalar and point work below another 2^50. This is a deliberately
loose ceiling, not a benchmark.

**Total.**

    W <= 10 x 1.001 x 2^121.5198 + 2^80 / 2000 = 2^124.8431701... < 2^124.84318.

`audit_g1.py` evaluates this with 80-digit decimal arithmetic.

## Memory, preprocessing and advice

Peak storage is the core's largest descent node, below 2^105 bytes
(STORAGE.md), plus below 2^30 bytes for the wrapper: below 2^106 bytes in all.
All preprocessing (factor base, relation collection, filtering, linear
algebra, field-map construction) is performed for this target and charged
once in W. No advice, stored relations or prior logarithms are used.

## Heuristics

H1-H5 of source-core.md section 9 apply unchanged, with H4u (at least 88% of
filtered row weight is +-1 entries; measured 0.90-0.95) and, in H3, the measured
non-increasing N/U trend behind pad 1.10 (source-core.md section 4). The yield
margin is 0.3012 bits; the descent's 2^63.97 relations hold at 1 relative pool
SE below the estimate (R8), with a collection guard.

## Evidence

`audit_g1.py` (declared initial experiment, input {}, about 2.5 s) recomputes
every changed quantity (see its docstring); follow-ups {"descent": true},
{"unit_entries": [k, seed]}, {"source": name}. `evidence16.py` (local, about
19 s, all 21 checks true) reruns every inherited certificate.

## Attribution

GT core: the lineage credited in source-core.md (Gautam Botrel with Claude Code
and Codex; refinements by kaankolcu, pepedesigner, bernardo-lima12 and
mpjunior92, most recently e3940d11). Randomized-base reduction and field map:
our 871a7c6c. R1, R8, the R9 joint optimisation and the H4u measurement: our
G1 lineage (27202622 ... 3eccc6c9). R2 and R4: our recover-P work
(43d603a4/be446c45 and its thread). R6, R7, R9's plan entries and the 1-SE level
of R8: our recover-Q thread (56ef135f). R10 (fully priced linear algebra,
`row_tail.c`, listings, N/U scale runs, pad 1.10): our recover-Q 711418d3,
unchanged.
