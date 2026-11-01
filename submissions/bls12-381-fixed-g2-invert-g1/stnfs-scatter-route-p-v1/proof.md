# Fixed-G2 pairing inversion into G1 (recover P) via the revised GT STNFS core with direct chunk scatter: 2^121.8542

## Status and scope

Track `bls12-381-fixed-g2-invert-g1`. Input: a nonidentity Q in G2, sampled
uniformly from the nonidentity points of G2, and z independently uniform in GT
(identity included). Output: P in G1 with e(P, Q) = z for that exact Q, where
e(P, Q) = py_ecc `pairing(Q, P)` at the pinned commit. No auxiliary powers. The
unit is one dense F_p12 multiplication, C = 20,000 scalar AArch64 instruction
equivalents, which is also the unit of the GT core ledger.

**Claim.** The algorithm below returns a valid P with success at least 0.39 over
uniform challenges and fresh coins, conditional on the disclosed NFS heuristics
H1-H5 of source-core.md (section 9). All classical work, including every
preprocessing step, both individual logarithms, all failures and restarts, the
field-basis conversion, both pairings, verification and output, is at most

    W <= (101/100) x 2^121.8398 + 2^80 / 20000 < 2^121.8542  Fp12 units.

The score is 121.8542. Peak storage is below 2^106 bytes. Nothing is amortized
over other targets and no advice is assumed. This is an analytic exploratory
claim under exploratory-v1, not a full-size run; acceptance would mean
plausible_not_refuted.

## Components

1. **GT core** (source-core.md, REVISION.md, REFINEMENT.md, DESCENT.md,
   KERNELS.md, STORAGE.md, code and data): an STNFS computation of
   a = log_g(z) in mu_r inside F_p12 for g uniform in GT minus {1} and z
   independently uniform in GT, with success at least 0.391949 and complete
   work at most 2^121.8398 dense-Fp12 units. It is the recover-Q core of
   submission e3940d11 with the two revisions of REVISION.md: R1, the
   main-phase change first published in 27202622 (99.9% per-prime ECM,
   theta = 207/500, slacks (300, 228), descent inputs 747.5/498.9, meter
   2^118.635), and R2, new here: records are written at emission directly to
   their chunk's page (the STORAGE.md listing, 21 instead of 36 instructions per record),
   allowance 1.0593. `audit_p.py`
   recomputes every changed quantity from the shipped code and data.
2. **Recover-P reduction** (below). Both pairing orientations need the same
   two logarithms in GT; only the base and the output group differ from
   recover-Q.
3. **Explicit field-basis conversion** (field-map.md) from the protected F_p12
   encoding of target-profile.json to the core's tower
   F_p[y, x]/(h(y), x^2 + |u| y), built afresh for each attack and charged.

## Algorithm

Use the pinned reduced pairing e and the public generator P0 of the target
profile.

1. Decode and check Q and z as in the acceptance rule (Q on the G2 curve,
   Q != infinity, [r]Q = infinity; z != 0, z^r = 1). If z = 1, output
   P = infinity (valid exactly for z = 1) and stop.
2. Compute g = e(P0, Q), i.e. `pairing(Q, P0)` (one pairing, protected
   encoding).
3. Build the field isomorphism of field-map.md (bounded randomized root
   finding, 12 x 12 basis inversion) and map g and z into the core tower.
4. Run the GT core on (g, z) under its global instruction meter and node
   counter. On any core failure, abort (failure).
5. The core returns a = log_g(z) mod r. Compute P = [a]P0 and check
   e(P, Q) = z with all 12 coefficients compared (second pairing). Output P
   in the canonical encoding if the check passes; otherwise fail.

## Correctness and the exact input distribution

G2 has prime order r and Q is a uniform nonidentity point, so Q = [b]Q0 with b
uniform in {1, ..., r-1}. Bilinearity and nondegeneracy give
g = e(P0, Q) = e(P0, Q0)^b, which is uniform on GT minus {1} because
b -> e(P0, Q0)^b is a bijection from {1, ..., r-1} onto GT minus {1}. The
target z is independent of Q and uniform on GT. So (g, z) has exactly the joint
distribution the core requires (uniform nonidentity base, independently
uniform target), with no re-randomization step. Since g generates GT,
z = g^a for a unique a mod r, and e([a]P0, Q) = e(P0, Q)^a = g^a = z, so
P = [a]P0 is a valid output for the supplied Q; it is the unique valid output,
because e(., Q) is a bijection from G1 onto GT. A field isomorphism preserves
z = g^a, so the conversion of step 3 does not change a, whichever roots it
selects. Returning P = infinity for z = 1 is the unique valid answer in that
case and only removes an input the core would otherwise have to handle. The
final pairing check rejects any wrong output. The returned point is a G1 point
for the supplied Q; no substitute Q and no bare scalar is returned.

## Success

The field-map construction fails with probability below 2^-48
(field-map.md). The core succeeds with probability at least 0.391949 on its
input distribution (source-core.md section 7: main phase 1 - 2^-192; two boots
2 x 1.7864e-04; meter Markov term 0.60589; node counter 2^-9.70; special-q
repeats 2^-10.70; 2/r). Hence

    P(success) >= 0.391949 - 2^-48 > 0.39,

with every core failure, abort and wrapper failure counted as a failure. The
z = 1 shortcut can only increase success. The core success uses Markov on the
two-log instruction meter and does not assume independence between the two
logarithms; the boot subtree independence is the disclosed H5.

## Work

**Core.** The upper-rounded core ledger (source-core.md section 8) is
collection 2^121.3481, linear algebra 2^119.36845, generator 2^106.4172,
two-log meter 2^118.635, canonical list 2^90, Schirokauer maps 2^88,
misc/filtering 2^80, factor base 2^72, certificate construction 2^50 and
verification 2^11, summing to 2^121.839785 < 2^121.8398 units of 20,000
instructions. The ledger is already in this track's unit; no conversion factor
is applied.

**Reserve.** The complete core ledger and every meter allowance are multiplied
by 101/100: some two-decimal kernel prices in measurements12.txt are rounded
down (e.g. 1914335/61528 > 31.11), every checked row by less than 1%. It is
not a reduction of any failure term.

**Wrapper.** At most 2^80 scalar instructions in total for the input checks,
the two pairings, the field-map construction (randomized factorization over
F_p12, 12 x 12 inversion, verification), input conversion, scalar arithmetic
mod r, the G1 scalar multiplication [a]P0, subgroup and encoding checks and
output. The elementary scalar bounds are those of field-map.md (Fp in seven 64-bit
limbs, schoolbook Fp12 below 2^31 instructions, each pairing below 2^50,
field map below 2^67). This is a deliberately loose ceiling, not a
benchmark; it contributes 2^80/20000 < 2^65.72 units.

**Total.**

    W <= 1.01 x 2^121.8398 + 2^80 / 20000 = 2^121.8541553... < 2^121.8542.

`audit_p.py` evaluates this with 80-digit decimal arithmetic.

## Memory, preprocessing and advice

Peak storage is the core's largest descent node, below 2^105 bytes
(STORAGE.md), plus below 2^30 bytes for the wrapper: below 2^106 bytes in all.
All preprocessing (factor base, relation collection, filtering, linear
algebra, field-map construction) is performed for this target and charged
once in W; the preprocessing is at most the whole bound W. No advice, stored
relations or prior logarithms are used.

## Heuristics

H1-H5 of source-core.md section 9 apply unchanged; this package adds none. R1
uses a stronger point of the same ECM per-curve model (`ecm16.py`); R2 is an
instruction count of an explicit listing (H4). The yield margin is 0.3012 bits
(about 2.2 pool SE) above the fixed 0.30-bit floor.

## Evidence

`audit_p.py` (declared initial experiment; no network, one CPU, 256 MiB:
checks_ok true, about 14 s) recomputes every changed core quantity,
the STORAGE.md listing hot-path count and the recover-P total. `evidence16.py` (optional; samples=16, seed=1) also reruns
every inherited certificate. REVISION.md lists each changed number.

## Attribution

GT core and all its machinery: the lineage credited in source-core.md (Gautam
Botrel with Claude Code and Codex; refinements by kaankolcu, pepedesigner,
bernardo-lima12 and mpjunior92, most recently e3940d11). This track's earlier
lean package (cbab3319, d08f8621) is the origin of the core. Main-phase
ECM/region revision and field map: kaankolcu's 27202622 and 871a7c6c (same
account). New here: the direct chunk scatter (REVISION.md R2), the recover-P reduction and success argument for this
track at the revised core bound, and the re-verification of the core and its
audit in this package.
