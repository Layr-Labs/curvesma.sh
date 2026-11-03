# Ordinary BLS12-381 G1 discrete logarithm via a randomized-base GT transfer and an STNFS core: 2^124.216030069061269

## Status and claim

Track `bls12-381-g1-ecdlp-plain`. Input T = [x]P0, x uniform in 0..r-1
(infinity included); output the canonical scalar x. Unit: one G1 reference
operation, C = 2,000 instructions. The algorithm below returns x with success
at least 0.39 over uniform challenges and fresh coins, conditional on the
disclosed heuristics H1-H5 (section Heuristics and limitations), and all
classical work, including every preprocessing step, both individual
logarithms, failures, the field map, pairings, sampling, verification and
output, is at most

    W <= 10 x (1001/1000) x 2^120.892660 + 2^80 / 2000 < 2^124.216030069061269 G1 units.

Peak storage is below 2^106 bytes; no advice, nothing amortized. Analytic
exploratory claim, not a full-size run.

## Components

1. **GT core** (source-core.md with DESCENT.md, KERNELS.md, STORAGE.md and the
   code): an STNFS computation of a = log_g(z) in mu_r for g uniform in GT
   minus {1} and z independently uniform, with success >= 0.410016 and work
   <= 2^120.892660 dense-Fp12 units (C = 20,000). Changes against our
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
< 2^-48 (field-map.md). The core succeeds with probability >= 0.410016 on its
input distribution under H1-H5 (section Heuristics): Block Wiedemann 2^-192
(H3); two boots 2 x 1.2578e-03 (independent kept relations given the main-phase
table, DESCENT.md step 2); meter Markov 0.578902 on the unrounded expected
two-log work; node counter 2^-7.96; repeats 2^-7.78 (special-q 2^-13.90,
missing-log ideals under the shared missing set 2^-7.81, missing-set
Chebyshev tail charged 2^-30, boot lattice shared 2^-25.95); 2/r. Hence

    P(success) >= (1 - 2^-128) x 0.410016 - 2^-48 > 0.39.

The meter is set for modelled success 0.41: the expected two-log work may
exceed its model value by 0.048 bits before the bound reaches 0.39. Markov
applies to the two-log meter, so no independence between the logs is used.

## Work

**Core.** The upper-rounded ledger rows of source-core.md section 8 sum to
2^120.8926595 <= 2^120.892660 units of 20,000 instructions. Converting to G1
units multiplies by exactly 20,000/2,000 = 10.

**Reserve.** The core ledger is multiplied by 1001/1000. Recomputing
collection with every quoted decimal constant raised by half an ulp moves it by
1.000015x and the ledger by 1.000014x (both checked by `audit_g1.py`
{"residual": true}), so the
reserve covers rounding with wide margin. It is not a failure term.

**Wrapper.** At most 2^80 scalar instructions for the two pairings, h^b,
sampling, field-map construction and conversion, scalar arithmetic, [x']P0,
checks and output (Fp in seven 64-bit limbs, schoolbook Fp12 below 2^31,
each pairing below 2^50, field map below 2^67).

**Total.** W <= 10 x 1.001 x 2^120.892660 + 2^80/2000 = 2^124.2160300690612688...
< 2^124.216030069061269, evaluated by `audit_g1.py` with 80-digit decimals.

## Memory, preprocessing and advice

Peak storage < 2^106 bytes (STORAGE.md); all preprocessing is charged once;
no advice, stored relations or prior logarithms.

## Side budgets (derived)

Fixed rows of `model16.total_log2` (units of 20,000 instructions), each above
an explicit count:

- **Canonical list, 2^90.** Per eligible a (< 2^63.8 in the 19/50 region):
  the HNF key of a O_h (6 x 6, entries < 2^11, exact: < 2^16 instructions) and a
  segmented sort of 16-byte keys (< 2^9 per key): < 2^80 instructions = 2^65.7 units.
- **Schirokauer maps, 2^88.** Per retained relation (< 2^64.2) and side,
  (a - b x)^E mod (r^2, h, S), E < r^72: < 18,400 multiplications of
  72-coefficient polynomials with 510-bit coefficients (< 2^20 each): < 2^99.4
  instructions = 2^85.1 units.
- **Filtering and miscellaneous, 2^80.** <= 64 entries per relation, < 2^70.2
  entries per pass, < 2^8 instructions per entry, < 2^6 passes for purge,
  clique removal, merge to weight 300 and back-substitution: < 2^84.2
  instructions = 2^70 units.
- **Factor base, 2^72.** Primes below 2^69 by a segmented sieve (< 2^72
  instructions), roots of h, f and g per prime (< 2^15): < 2^78.4 = 2^64.1 units.
- **Verification 2^11**: the core's check g^a = z, one 255-bit exponentiation
  in F_p12 (< 2^10 units). **Certificates 2^50**: the one-time exact
  certificates ([0]-[0c], [3], [13]) are < 2^40 instructions. The G1 wrapper
  (pairings, field map, [x']P0, output) is charged separately: 2^80 instructions.

## Heuristics and limitations

H1 smoothness (Dickman/Murphy) and H1' threshold retention, with random ECM
group orders (collection, band strip, descent); H2 uniform pool and the
two-standard-error work allowance and descent inputs; H3 the virtual-logarithm
system: r unramified in K_f and K_g ([13]), r dividing neither class number,
the 21 sigma-odd Schirokauer maps injective on units modulo r-th powers, so the
kernel of A = [M | S] is the line of virtual logarithms (Block Wiedemann then
fails with probability <= 4N/r); N/U at production not above the 2^15-2^18
values; H4 measured kernel counts and H4u (+-1 share >= 0.88); H5 the descent
model (DESCENT.md): the Poisson occurrence model, H5c (ideals occurring at
least three times keep their logarithm), the negative-binomial relation
counts, the measured unit-duplicate charge, the worst-sextile norm envelopes
and H5a (independent factorisations of distinct cells within a node). These
are extrapolated from measurements at smaller scale (filtering to 2^18, toy
balls, sampled production lattices); none is proved. `audit_g1.py` ({}: about
10 s; {"descent": true}: about 9 s) recomputes every number above;
`evidence16.py` (local) reruns the inherited certificates.

## Attribution

GT core lineage: source-core.md (Provenance). Pool data, alpha_g and the
owner-count tail form: gbotrel 42e80f80. Band strip: our 8cc93834 lineage.
Sub-cell integration: our plain R18 (245a6b54) after our recover-Q v14 line;
complete-orbit UDUP: our plain R18b (eb01512c, method after 57598408); success
apparatus and side-budget derivations: our recover-P 273b1519;
occupied-owner recount: our f34f63d9. Prime special-q missing logs: the
approach gbotrel described for recover-Q (cd991163), implemented here.
Reduction and field map: our 871a7c6c. The rest: our d605c287 lineage.
