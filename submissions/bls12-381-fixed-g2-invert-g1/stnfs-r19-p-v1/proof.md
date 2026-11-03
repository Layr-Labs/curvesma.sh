# Fixed-G2 pairing inversion into G1 (recover P) via the revised GT STNFS core, filter-pad re-tune (R19): 2^121.260301167243655

## Status and scope

Track `bls12-381-fixed-g2-invert-g1`. Input: Q uniform over the nonidentity
points of G2 and z independently uniform in GT (identity included). Output: P in G1 with e(P, Q) = z for that exact Q, where e(P, Q) =
py_ecc `pairing(Q, P)` at the pinned commit. No auxiliary powers. Unit: one dense
F_p12 multiplication, C = 20,000 scalar AArch64 instructions, also the unit of
the GT core ledger.

**Claim.** The algorithm below returns a valid P with success at least 0.39 over
uniform challenges and fresh coins, conditional on the NFS heuristics H1-H5 of
source-core.md (section 9). All classical work, including every preprocessing
step, both individual logarithms, all failures and restarts, the field-basis
conversion, both pairings, verification and output, is at most

    W <= (2001/2000) x 2^121.25958 + 2^80 / 20000 < 2^121.260301167243655  Fp12 units.

The score is the tight fifteen-decimal ceiling 121.260301167243655. Peak storage is below 2^106 bytes. Nothing is amortized
and no advice is assumed. Analytic exploratory claim under exploratory-v1, not a
full-size run.

## Components

1. **GT core** (source-core.md with REVISION.md, REFINEMENT.md, DESCENT.md,
   KERNELS.md, STORAGE.md, code and data): a = log_g(z) in mu_r inside F_p12 for g
   uniform in GT minus {1} and z independently uniform, success >= 0.3917, work
   <= 2^121.28540 units: the e3940d11 core with R1-R13 (exact prices, 10.005
   route, 20/13-entry LA with 22x4 blocks, 0.05% reserve, refiltered N/U 0.0614,
   kernels 850/2000, rels allowance, corrected loss,
   re-optimised plan), R14 (2^64 sieve, theta 19/50, slacks 293/222, meter 118.50533, numerical fix),
   R15 (split 2^20 -> 2^48), R16 (band-strip for (2^20, 2^48], allowance 1.0700),
   R17 (plan re-opt, odd-k grid), R18 (excess/filter/descent re-tune), and new here R19 (pad 1.2).
2. **Recover-P reduction** (below). Both pairing orientations need the same two
   logarithms in GT; only the base and the output group differ from recover-Q.
3. **Explicit field-basis conversion** (section below) from the protected F_p12
   encoding to the core tower F_p[y, x]/(h(y), x^2 + |u| y), built afresh for
   each attack and charged.

## Algorithm

Use the pinned reduced pairing e and the public generator P0.

1. Decode and check Q and z as in the acceptance rule (Q on the G2 curve,
   Q != infinity, [r]Q = infinity; z != 0, z^r = 1). If z = 1, output
   P = infinity (valid exactly for z = 1) and stop.
2. Compute g = e(P0, Q), i.e. `pairing(Q, P0)`.
3. Build the field isomorphism below and map g and z into the core tower.
4. Run the GT core on (g, z) under its global instruction meter and node
   counter. On any core failure, abort (failure).
5. The core returns a = log_g(z) mod r. Compute P = [a]P0, check e(P, Q) = z
   on all 12 coefficients (second pairing) and output P canonically if the
   check passes; otherwise fail.

## Correctness and the exact input distribution

G2 has prime order r and Q is uniform nonidentity, so Q = [b]Q0 with b uniform in
{1, ..., r-1} and g = e(P0, Q) = e(P0, Q0)^b is uniform on GT minus {1}
(b -> e(P0, Q0)^b is a bijection onto GT minus {1}). z is independent of Q and
uniform on GT, so (g, z) has exactly the core's input distribution, with no
re-randomization. Since g generates GT, z = g^a for a unique a mod r, and
e([a]P0, Q) = g^a = z, so P = [a]P0 is valid for the supplied Q and unique
(e(., Q) is a bijection G1 -> GT). A field isomorphism preserves z = g^a, so
step 3 does not change a, whichever roots it selects. P = infinity is the unique
answer for z = 1. The final check rejects any wrong output; no substitute Q and no
bare scalar is returned.

## Field-basis conversion

Let F be the protected representation of F_(p^12). Factor h(Y) = Y^6 + Y^3 +
Y^2 - Y - 1 in F[Y] into linear factors and choose a root y (h is irreducible of
degree 6 over Fp and 6 | 12, so it splits in F); then find a root x in F of
X^2 + |u| y (irreducible over F_(p^6), so it splits in F with x outside that
subfield; irreducibility is recomputed by evidence16.py [0], [0b]). Y -> y,
X -> x is a field isomorphism from the core tower onto F (a unital map of a field
of dimension 12 into F). The 12 vectors y^i x^j form an Fp basis; invert their
12 x 12 coordinate matrix by Gaussian elimination and map g and z with it.

Splitting: for a squarefree completely split f of degree d <= 6, pick R of
degree < d with uniform F coefficients and take gcd(R, f) and
gcd(R^((p^12-1)/2) - 1, f). Values at distinct roots are jointly uniform
(Vandermonde), so two roots get opposite quadratic characters with probability
((|F|-1)/|F|)^2/2 > 1/4: each trial splits with probability > 1/4. Five splits
for h and one for the quadratic, each capped at 128 trials: failure
<= 6 (3/4)^128 < 2^-50. Uniform F elements use 381-bit rejection sampling of 12
Fp coordinates (128 attempts each, failure < 2^-128 per coordinate, < 2^16
coordinates: < 2^-112). Total field-map failure allowance: 2^-48.

Cost ceiling: a multiplication mod a degree <= 6 polynomial is < 2^10 field
operations, powering with < 4572 bits < 2^24, gcds < 2^24; charge 2^26 Fp12
operations per trial, < 2^10 trials: < 2^36 Fp12 operations at < 2^31 scalar
instructions each (Fp in seven 64-bit limbs, schoolbook): < 2^67 instructions.
Each pairing is below 2^50 instructions. The maps occupy < 2^30 bytes.

## Success

The field map fails with probability below 2^-48. The core succeeds with
probability at least 0.3917 on its input distribution (source-core.md
section 7: main phase 1 - 2^-192; two boots 2 x 2.4860e-03; meter Markov term
0.60119; node counter 2^-9.50; special-q repeats 2^-10.70; 2/r). Hence

    P(success) >= 0.3917 - 2^-48 > 0.39,

with every core failure, abort and wrapper failure counted as a failure. Markov on
the two-log meter needs no independence between the two logarithms; boot subtree
independence is the disclosed H5. R3-R8 change no success term; R9 re-derives the descent terms; R10 corrects the loss proof; R13 re-derives the allowances; R14 re-tunes sieve/region; R15 moves the sieve split; R16 adds the band-strip (success terms unchanged); R17 re-opts plan/grid (boot phi down, success up); R18 re-tunes excess/filter/descent (success 0.39191); R19 touches LA only.

## Work

**Core.** Upper-rounded core ledger (source-core.md section 8): collection
2^120.61756, linear algebra 2^118.794879, generator 2^109.009417, meter
2^118.76810, canonical list 2^90, Schirokauer maps 2^88, misc/filtering 2^80,
factor base 2^72, certificate construction 2^50, verification 2^11: sum
2^121.25957361 < 2^121.25958 units.

**Reserve.** The ledger and meter allowances are multiplied by 2001/2000.
Raw-counted prices are exact (`prices_exact`); the rest is printed precision
of five constants within half an ulp each, a li-term allowance and binary64
error. The audit (`reserve_residual`) weights each by core share: undercharge
1.32e-5, covered tenfold. Not a failure-term reduction.

**Wrapper.** At most 2^80 scalar instructions for the input checks, two pairings,
field-map construction (above, < 2^67), input conversion, scalar arithmetic mod
r, [a]P0, subgroup and encoding checks and output: a deliberately loose ceiling,
2^80/20000 < 2^65.72 units.

**Total.**

    W <= 1.0005 x 2^121.25958 + 2^80 / 20000 = 2^121.260301... < 2^121.260301167243655.

`audit_p.py` evaluates this with 80-digit decimal arithmetic.

**Memory, preprocessing and advice**

Peak storage is the core's largest descent node, below 2^105 bytes (STORAGE.md),
plus below 2^30 bytes for the wrapper. All preprocessing is performed for this
target and charged once in W. No advice or stored relations are used.

## Heuristics and evidence

H1-H5 of source-core.md section 9 apply unchanged; R3 is a measured count of an
explicit program, R4 an explicit charge, R5 exact arithmetic on shipped raw counts
plus a repriced route, finishing budget and reserve, R6/R11/R18 re-simulated filters
with trend-backed pads, R7 budgeted transfers, R8 tighter descent kernels, R9 rels allowance, R10 corrected loss, R13 union allowances, R14 re-tuned sieve/region. Yield margin 0.30119 bits (about 1.9 pool
SE) above the 0.30-bit floor. `audit_p.py` (initial experiment; one CPU, 256 MiB,
about 13 s) recomputes every changed quantity, the R3 count from the
measurements12.txt R3 lines, the R4 budget and the recover-P total;
`evidence16.py` (optional; samples=16, seed=1) reruns every inherited certificate.

## Attribution

GT core and all its machinery: the lineage credited in source-core.md (Gautam
Botrel with Claude Code and Codex; refinements by kaankolcu, pepedesigner,
bernardo-lima12 and mpjunior92, most recently e3940d11). This track's lean
package (cbab3319, d08f8621) is the origin of the core. Main-phase ECM/region
revision and field map: kaankolcu's 27202622 and 871a7c6c. R2: be446c45. R3 and
R4, the recover-P argument at the revised core bound, and the declaration
refinement: 115f3fd9, 5f39af0 and 339cb2b7 (this package's base). Earlier: R5 (exact
kernel prices, 10.01 route, 13-per-entry LA finishing, 23x4 blocks, 0.05%
reserve), R6 (fresh filtering basis, 1.2 pad, 100.7/entry), R7 (transfer
budget, 10.005 route, 22x4 blocks), R8 (descent kernels), R9 (rels allowance, meter), R10 (loss proof), R11
(2^19 filter basis, 1.25 pad, retention allowances), all by muse-spark-1.3-contributor (Muse Code via T3 Code);
R12-R18 as in REVISION.md, new here R19 (pad 1.2).
