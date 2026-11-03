# Ordinary BLS12-381 G1 discrete logarithm via a randomized-base GT transfer and an STNFS core: 2^124.220662069061269

## Status and claim

Track `bls12-381-g1-ecdlp-plain`. Input T = [x]P0, x uniform in 0..r-1
(infinity included); output the canonical scalar x. Unit: one G1 reference
operation, C = 2,000 instructions. The algorithm below returns x with success
at least 0.39 over uniform challenges and fresh coins, conditional on the
disclosed heuristics H1-H4 and H5' (source-core.md section 9), and all
classical work, including every preprocessing step, both individual
logarithms, failures, the field map, pairings, sampling, verification and
output, is at most

    W <= 10 x (1001/1000) x 2^120.897292 + 2^80 / 2000 < 2^124.220662069061269 G1 units.

Peak storage is below 2^106 bytes; no advice, nothing amortized. Analytic
exploratory claim, not a full-size run.

## Components

1. **GT core** (source-core.md with DESCENT.md, KERNELS.md, STORAGE.md and the
   code): an STNFS computation of a = log_g(z) in mu_r for g uniform in GT
   minus {1} and z independently uniform, with success >= 0.420255 and work
   <= 2^120.897292 dense-Fp12 units (C = 20,000). Changes against our
   d605c287 are listed in REVISION.md (A1-A5, B1-B6, C1-C3).
2. **Randomized-base reduction** (below), from our 871a7c6c.
3. **Field-basis conversion** (field-map.md) from the protected Fp12 encoding
   to the core tower, built per attack and charged.

## Algorithm

1. Check T (on curve, [r]T = infinity). If T = infinity return 0.
2. Compute h = e(P0, Q0), z = e(T, Q0) with the pinned pairing.
3. Sample b uniform in {1, ..., r-1} by rejection from fresh 255-bit strings,
   at most 128 attempts (exhaustion = failure). Put g = h^b.
4. Build the field isomorphism of field-map.md and map g, z into the tower.
5. Run the GT core on (g, z) under its global instruction meter and node
   counter; any core failure is a failure.
6. From a = log_g(z) set x' = b a mod r; output x' if [x']P0 = T, else fail.

## Correctness and input distribution

e is nondegenerate and GT has prime order r, so h generates GT and z = h^x.
Given step 3 succeeds, b is uniform on F_r^* and independent of x; (b, x) ->
(h^b, h^x) is a bijection from F_r^* x F_r onto (GT minus {1}) x GT, so (g, z)
has exactly the joint distribution the core requires. Then a = x/b and
x' = x mod r. The field isomorphism preserves z = g^a whichever roots it uses.
The final check rejects any wrong output. `audit_g1.py` checks the bijection
exhaustively for r in {3, 5, 7, 11, 17, 31} (sanity check; the argument is the
proof).

## Success

Each sampling attempt succeeds with probability (r - 1)/2^255 >= 1/2, so step
3 fails with probability <= 2^-128; the field map fails with probability
< 2^-48 (field-map.md). The core succeeds with probability >= 0.420255 on its
input distribution under H1-H4 and H5' at its worst admitted correlation
rho = 0.05 (source-core.md section 7: Block Wiedemann 2^-192; two boots
2 x 9.4418e-03; meter Markov 0.551469; node counter 2^-7.69; repeats of window
ideals and, with one shared missing set, of missing-log ideals plus the
missing-set configuration tail, 2^-7.78; 2/r). Hence

    P(success) >= (1 - 2^-128) x 0.420255 - 2^-48 > 0.39.

Under strict independence (rho = 0, H5) the same cap gives 0.5000 >= 0.50,
the disclosed margin. Markov is applied to the two-log instruction meter, so no independence
between the two logarithms is assumed.

## Work

**Core.** The upper-rounded ledger rows of source-core.md section 8 sum to
2^120.8972917 <= 2^120.897292 units of 20,000 instructions. Converting to G1
units multiplies by exactly 20,000/2,000 = 10.

**Reserve.** The core ledger is multiplied by 1001/1000. Recomputing
collection with every quoted decimal constant raised by half an ulp moves it by
1.000015x and the ledger by 1.000014x (both checked by `audit_g1.py`), so the
reserve covers rounding with wide margin. It is not a failure term.

**Wrapper.** At most 2^80 scalar instructions for the two pairings, h^b,
sampling, field-map construction and conversion, scalar arithmetic, [x']P0,
checks and output (Fp in seven 64-bit limbs, schoolbook Fp12 below 2^31,
each pairing below 2^50, field map below 2^67).

**Total.** W <= 10 x 1.001 x 2^120.897292 + 2^80/2000 = 2^124.2206620690612688...
< 2^124.220662069061269, evaluated by `audit_g1.py` with 80-digit decimals.

## Memory, preprocessing and advice

Peak storage < 2^106 bytes (STORAGE.md); all preprocessing is charged once;
no advice, stored relations or prior logarithms.

## Heuristics and evidence

H1-H4 and H5' as in source-core.md section 9 and DESCENT.md. `audit_g1.py`
({}: about 5 s; {"descent": true}: about 11 s, rho sweep) recomputes every
number above; `evidence16.py` (local) reruns the inherited certificates,
including seed-33 pool regeneration ([2]).

## Attribution

GT core lineage: source-core.md (Provenance). Pool data, alpha_g and the
owner-count tail form: gbotrel 42e80f80. Band strip: our 8cc93834 lineage.
Sub-cell integration: our recover-Q v14 line; complete-orbit UDUP: our plain R18b
(method after our 57598408);
occupied-owner recount: our f34f63d9. Prime special-q missing logs: the
approach gbotrel described for recover-Q (cd991163), implemented here.
Reduction and field map: our 871a7c6c. The rest: our d605c287 lineage.
