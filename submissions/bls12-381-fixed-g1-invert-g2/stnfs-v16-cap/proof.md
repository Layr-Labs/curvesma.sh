# recover-Q (bls12-381) through two GT logarithms with the plain-track STNFS core (plain P19): 2^119.191213809501 reference units

Track `bls12-381-fixed-g1-invert-g2`. Input: a nonidentity P in G1 and z in GT; output Q in G2 with e(P,Q) = z for
the exact supplied P (Q = infinity for z = 1). The score unit is one dense F_p12 multiplication, C = 20,000 scalar
AArch64 instruction equivalents, the unit the core below counts. No advice, smaller targets or amortization are used.
Predecessor work is reused as credited at the end; no prior verdict is a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied
nonidentity P (Q = infinity for z = 1) with success at least 0.39 over uniform
challenges and fresh coins under the declared heuristics H1-H5 (section 9). All classical work, including
preprocessing, both individual logarithms, failures, restarts, the field map, the pairing, the scalar multiplication,
verification and output, is at most

    W <= 2^119.1912138095 + 2^66 <= 2^119.191213809501  dense-Fp12 units (C = 20,000).

The score is 119.191213809501. 2^119.1912138095 is the GT core's global meter (sections 2-8): upper-rounded ledger
rows sum to at most 2^119.19121380943; 2^66 is the wrapper (section 1). Storage is below 2^106 bytes. The
claim is analytic and exploratory; no full-size run is made.

**Core.** Our 0f679202 (122.537206169587363) = our a05affb1 (122.540321783687363: plain P10b 9c904850 re-optimised
by our REOPT-2000 study at fixed kernels: line sieve from 2^14 with trial division below 2^20, slacks 302/233, region
Q_a <= (2895/8000) A^2, Q_b <= (17/20) A^2 - Q_a, pools S/G (PREREG-G), clean-run rec7, band charge 0.009 bits,
descent hardening, canonical row 2^91; coupled parameters: end of section 3) with mpjunior92's seven matrix owners
(6734bd5a): collection **2^118.895525055**, meter **2^113.279**. **This package (P19)** adds three linear-algebra
changes over P18 (Horner chains, section 4): the **PREREG-H** filtering extension (N(300)/U 0.03672 -> 0.03577),
**Horner chains** from Krylov checkpoints (section 4), **Block Wiedemann m = 27** (the interior minimiser of the
package's own LA formula at n = 5; P18 kept m = 28 from 0f679202, its own proof noting m = 27 lower by 4e-5 bits) and
the **descent meter set to the least value its own success budget supports** (2^113.278864326817; P18's 2^113.279 was
rounded up from the same constraint): linear algebra 2^116.621811206 -> **2^116.621741985**, generator
2^107.328830867 -> **2^107.208650296**, descent meter 2^113.279 -> **2^113.278864326817**. Every other row is
0f679202's.

**Conventions.** Bounds round safely; `evidence16.py` recomputes the checks cited as [n] and every stated bound;
`audit_refinement.py` is the hosted audit. Model: `config.json`, `stats16.b85`, `pool_d.json`, `model16.py`,
`descent16.py` (`owner10.install` sets C_ROW 1,320).

## 1. Reduction (fixed P: two GT logarithms) and field map

With the pinned reduced pairing e and the public generator Q0 of G2:

1. Decode and check P (on curve, [r]P = infinity, P != infinity) and z (z^r = 1); z = 1 returns Q = infinity.
2. g = e(P, Q0) (one pairing, pinned `pairing(Q0, P)`, protected encoding).
3. Build the field isomorphism of field-map.md (bounded randomized root finding, 12 x 12 basis inversion) and map
   g and z into the core's tower.
4. Run the core on (g, z) under its meter and node counter; any core failure is a failure; it returns
   a = log_g(z) mod r.
5. Q = [a]Q0 (one G2 scalar multiplication); output Q if e(P, Q) = z, otherwise fail.

GT has prime order r and e is nondegenerate, so for the exact supplied P (never replaced) g generates GT, z = g^a for
a unique a mod r and e(P, [a]Q0) = z. For uniform nonidentity P, g is uniform on GT minus {1} and independent of the
uniform z: the core's input law. The field isomorphism preserves z = g^a; the final check rejects any wrong output.

One relation/LA phase serves log z and log g (section 6); each target t is randomised as t^e, e uniform in [1, r)
(2/r charged); z = 1 never reaches the core.

Wrapper work below 2^80 instructions (Fp in seven limbs, products < 2^22, Fp12 < 2^31, exponentiations < 2^14 Fp12
products; pairing below 2^50, field map below 2^67 by field-map.md, input checks, scalar multiplication and check
below 2^50): 2^80 / 20,000 = 2^65.71 units, a **2^66** ledger row (section 8); wrapper failure: field map below 2^-48.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is squarefree, so O_h = Z[y]; h is irreducible
mod p. f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks (u y)^((p^6-1)/2) = -1 in F_p6, so F_p[y,x]/(h, g) = F_p12;
[0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P);
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at (5.15, 0.85) in collection and descent;
[14] recomputes it from the exact Kummer-Dedekind structure: (4.868096, 0.843643) over p <= 2000 and (5.085035,
0.731543) over p < 50000, both below the charge.

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The factor base is the relative-degree-1
prime ideals over K_h of norm <= B on both sides; U = 2^63.451 sigma-orbits.

**Region.** Exact rational forms Q_a, Q_b (`config.json`), A^2 = 2^21. Pool: uniform on theta = 3/8 (a5 > 0,
Q_a(a) <= (3/8) A^2; b uniform in the half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a)). Region: a5 > 0, Q_a(a) <=
(2895/8000) A^2; b5 >= 1, Q_b(b) <= (17/20) A^2 - Q_a(a) (invariant under b -> -b and the sign fold). The restriction
of the uniform sample is uniform on the region; G has 35244 region rows (32197 primitive) of 78960 canonical draws,
scaled by the region volume over the effective sample 656000 x (2895/3000)^3 = 589502.674. Two weights, never mixed: the density w = (1 - qa)^3 (yield)
and vb (1 - qa)^3 x factor (work, population); every slice-geometry cost uses S' = vb (17/20 - qa)^3 x factor. The
canonical a of each unit class is the least of its exact HNF key class in (Q_a, lex) order; non-primitive pairs are
excluded exactly. `lattice_count_bounds.py` certifies (exact
rational Machin-pi, determinant/sqrt intervals) the count ratio at every radius R >= 871, below the smallest b-slice
radius 1011 ([3]): [0.985511, 1.013857] (Q_a), [0.985464, 1.013857] (Q_b); the model uses the exact products
**0.971187233 (yield)** and **1.027905204 (work)**. The canonical list costs 2^91 units.

**Pools (PREREG-G).** 656,000 uniform eligible pairs each on theta = 3/8 (indices 0..163999 per seed): S
(seeds 261-264) sizes the region, G (251-254, drawn after the configuration and theta were fixed; else only in [10y],
[11c], H4u) is judged once, S and G give the work; every draw was tested with the exact canonical predicate.
`pool_d.json` carries per-seed counts and list hashes, per-pool sums, 12 anchor rows of G, each G seed's canonical
list over draws 0..49 ([2] regenerates and replays them) and the contamination sample. theta is the smallest t on
the 1/8000 grid with margin_S >= 0.06 + sqrt(2) SE_S: 2895/8000 ([4s]); the band charge applies to G's count.

**Sieve, thresholds, cofactorisation.** Line sieve (scores only) for p in [2^14, 2^48); record walk for p in [2^48,
2^64) (`slice12`, KERNELS.md; X2 = 2^64; [2^63, 2^64) on the measured wide path; [2^56, 2^63) per a-run); primes in [2^64, 2^69] by ECM within the slacks (302, 233). Primes below 2^14 are not sieved: their mass is the retention models' omitted mass (`prime_min` 2^14,
`explicit_bits` 14); their factors and those of [2^14, 2^20] are found at the survivors by **trial division by every
prime below 2^20 on both cofactors** (82,025 primes, 13 + 9 limbs at the 832/576-bit caps): 2^23.649 <=
**C_TDIV = 2^24.5** per survivor ([5t]). Primes in (2^20, 2^48] only add line scores and emit no records (>= 2^5.9
hits per line coset at p = 2^48); their factors are recovered at the survivors by a band strip: two rounds per side of the `ecm16.py` 48-bit schedule
at 99.9% on the unstripped cofactor, charged at ceil(log2|N|) + 1 bits (2 x 2^36.4488 f, 2 x 2^35.3903 g); a band
prime missed twice loses the relation: yield x (1 - 41 x 0.001^2)(1 - 28 x 0.001^2). Survivors'
cofactors are bounded by 310 (f) and 241 (g) bits; per-survivor ECM is the `ecm16.py` schedule at **99.9%** per prime
([12]): 2^38.0909 instructions at 310 bits, 2^38.0882 at 241 bits. **Early abort (2c2ce2ab):** the g side's band strip and ECM run only after the f side factors completely (yield
unchanged): per survivor C_NORMS + C_TDIV + band_f + 2^38.0909 + q_f (band_g + 2^38.0882), q_f =
min(1, 100 rho(u_f)/P(omitted_f <= slack_f)) per pool row (H1, mean ~1e-5). Certified per-block
bounds (`embedding_bounds.py`, 64-cell blocks of width 128, bisected until both bounds are within 8 bits) give
thresholds T = lower - 16 x slack; F_v is Taylor-shifted to the block centre per place, the remainder over the block
disk plus a 2^-32 relative term, rounded outward; `nbound31` replays one evaluation bit for bit at 1,508.2 Ir
(charged as measured). **Box bounds (recover-Q v8).** Rows of one line in an aligned 2 x 2 x 4 box of (b2, b3, b4) (`model16.box_rows`)
share one bisection per block over their contexts' hull (`EB.hull`, `EB.box_leaves`); an unresolved column uses each
row's own context. Boxes x blocks x (b1, b5) partition each slice, so C_NB_EVAL x the box's evaluations per box cell
has the right expectation; each row pays its leaves, context and hull share (C_NB_ROW 4,600). [2b] checks 1,536 exact
norms of six rebuilt boxes.

**Record walk per canonical a-run (rec7; KERNELS.md).** A canonical a's run is the maximal t-interval containing 0
with a + t e_0 in the region and no image of smaller Q_a under the 11,325 listed units (`runs16.py`; [2] replays the
anchors' run ends); a prime ideal's hits in all its slices form one coset of a 7-dim lattice, walked once per (run,
ideal) with one setup6. **Clean-run rule:** rec7 walks a run only if every a-value of it is in the canonical list (whose exact test the 2^91
row pays); otherwise all its slices take the slice walk. One implementation per run and band: a clean run of 32 to
1,024 a-values uses the run walk in each band p in [2^56, 2^63) with h' = run cells / p >= 2^4.99 in place of the
slice walk and empty cosets. Run cells are exact (`runs16.lrs`) and each row pays its cells' share. **Contamination (as recP v12gp 5aad7567).** A contaminated run's slices cost at most the slice walk: per slice
<= D_SUP = `contam_sup` = 2^72.211 on [2^56, 2^63) (each band's record walk at the exact max of h x REC_TAB(h), up to
the largest region slice, plus empty cosets and setup6 in full). Every row pays RUN_ALLOW x D_SUP per slice,
RUN_ALLOW = 0.003: pool draws are uniform over canonical a, and in an a-uniform sample of 1,000 runs the exact
predicate found no non-canonical a-value among 404,564, so contaminated runs hold under 1 - 0.05^(1/1000) <= 0.002992
of the canonical a (95%, a modeling margin). [k7c] replays 3 sampled runs, [2r] 3 anchor runs. Count envelopes ([3]) leave the per-a charge >= 0.002876 D_SUP;
a slice's work is <= 0.9373 D_SUP and the need <= 0.002804 D_SUP (2.5% to spare). The rule's own work is < 2^69.3
units, in the 2^80 misc row.

**Yield ([4], pool G).** rels = (va_region/589502.674) vb (0.971187233) sum over G's canonical primitive region
rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g): Dickman rho ([1]), retentions r from the finite
Euler-product threshold model (`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14;
prime_min 2^14), the ECM miss and the band-strip loss; cells with qa < 0.02 are sieved and charged, their relations
discarded. Estimate 2^63.63225 (SE 0.0234 bits) less the band charge 0.009: retained relations >= **2^63.6232**
against 1.05 U <= 2^63.52173: margin >= 0.1015 bits = 4.3 pool SE; feasibility requires max(0.06, 2 SE).
**Band charge (H1', [11c]).** [11c] factors the exact norms of 6,000 region pairs of G's seeds over [2^8, 2^14)
(units keep norms, so canonicity is immaterial): f 5.876 +- 0.105 bits, g 6.130 +- 0.108 against the per-prime law's
5.923; above 40 bits 12 (f) and 19 (g) pairs against the law's 0.00194. The bound takes the mean excess at 2 SE as a
slack guard through the retention model at G's relation-weighted means and the excess above 40 bits at its Poisson
95% limit as lost: 0.0075 <= REL_BAND_LOSS = 0.009 bits, deducted in [4]. **Descent input.** `LOG2_RELS = 63.575` is at or below the charged count lowered by two relative pool SE, 2^63.57570
([4] `rels_input`); the mean norms 739.1/498.2 truncate G's relation-weighted 739.167/498.217 (lower means: covered by
[10y], pooled ratio + 2 SE <= 1).

**Cost ([5]).** Sampled mean 190.2 instructions/cell over pools S and G, charged 191.75 with the allowance, over
2^125.600 cells (bounds, layout, setups, line path, record walk or rec7 with its contamination charge, routing,
apply, recovery, scan2 with fused zeroing, band strip, cofactorisation). Collection <= 2^118.895525055 units.

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo estimate of a fixed, challenge-independent
domain total over the 1,312,000 draws of pools S and G; its relative SE has a conditional row term 0.1828% and an
incidence term 0.3664%, total <= 0.4095%. The charge multiplies the estimate by **1.0082 >= 1 + 2.00 x 0.004095**
(rounded up); [4a] fails below 1 + 2.00 SE.

**Coupled parameters (vs P10c = P10b + descent hardening; not submitted).** Shape constants: `config.json`; ledger
constants: `model16.py`.

Coupled parameters vs P10c (line-sieve start 2^14 with trial division below 2^20, retention prime_min 2^14 / explicit_bits 14, slacks 302/233 with cofactor bounds 310/241 and ECM 2^38.0909/2^38.0882, region 2895/8000 with Q_b <= (17/20) A^2 - Q_a, pools S/G sized/judged/work = S/G/S+G at theta_pool 3/8, count factors 0.971187233/1.027905204 at radius 871, exact run cells with clean-run rec7 + 0.003 D_SUP/slice contamination, filtering to 2^20 (PREREG-G/H) with NU300 0.03577, relation charge 0.009 bits at means 739.1/498.2 (G), LOG2_RELS 63.575, meter cap 113.278864326817, canonical row 2^91, largest half-slice 2^66.768).

Unchanged: B = 2^69, X2 = 2^64, P_MID 2^48, wide path, band strip, alpha, RUN_T, purge, the LA kernels and
formula, the descent plan and its own sieve from 2^8 (descent16.P.PMIN).

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U, `filtsim.PEXC`), clique removal and Markowitz merge (coefficients <= 255), simulated on
the production occurrence profile (`filtsim.py`; [6] reruns 2^12) as the signed sigma-quotient system, charged at
excess 1.05 (every retained relation is used, which only over-states N). A campaign at P10c's judged means 740.8/499.4 (`FILT_RUNS`: 10, 6, 6, 6, 12, 16, 16, 16, 8 seeds at 2^12..2^20; 2^18 seeds 9-16 and 2^19 seeds 1-8 are
5170c7d8's PREREG-G, 2^19 seeds 9-16 and 2^20 seeds 1-8 our recover-Q lane's PREREG-H, each rule written before its
runs) measures x = N(300)/U, merged weight-300 rows over U: 0.06812 (SE 0.00233) at 2^12 down to 0.03394 (0.00022)
at 2^20. H3 assumes the production expectation does not rise above
the largest simulated scales; the pre-registered charge is the larger of mean + 2 SE at the two largest scales,
max(0.035763 at 2^19, 0.034382 at 2^20) -> **0.03577** ([10]; the eight new 2^19 seeds averaged above the first
eight, mean 0.03498 -> 0.03540). At G's means 739.1/498.2 the same seeds give pooled paired ratio 0.9806 +- 0.0076
(`FILT_Y`, [10y]; 0.9957 <= 1 at 2 SE): the charge covers G's means. Merge to mean weight 300
(charged 301 per row), largest merged coefficient bound 253 <= 255. N = 0.03577 U <= **2^58.6463** < 7 x 2^56:
seven owners (row tail 1,317, charged 1,320).

The square system A = [M | S] has N + 21 rows, N ideal columns and 21 dense Schirokauer columns (sigma-odd unit
ranks 18 + 3 from signatures K_f (0, 36), K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3 its kernel is the
line of virtual logs. Block Wiedemann: m = 27 left vectors X with entries uniform in [0, 2^128), n = 5 right vectors Y uniform
in F_r; Krylov N/m + N/n block steps. Per entry, row and element (STORAGE.md, [10]): Krylov products 13.21 per vector-entry (H4u: >= 87% +-1 entries,
0.8809-0.9366 in 44 runs at G's means), Horner (solution) products in five interleaved chains, chain j from the
Krylov checkpoint A^(jL) Y (L = Horner length / 5; 20 checkpoint vectors copied in the Krylov pass, 4L < N/m + N/n,
< 2^68 bytes), stored interleaved, so they run the five-vector pass: 301 x 13.21 + 1,320 (row tail, seven owners) +
21 x 97 per row and chain step plus n dense terms and one reduction (420) per element, projections (m + n) N^2 terms
at 60. Linear algebra **<= 2^116.621741985 units** (m = 27, the interior minimiser of this formula at n = 5 over m in
[20, 40]: m = 26 is 5.2e-5 bits higher, m = 28 is 6.9e-5 higher; the same formula and measured kernels as P18, which
kept m = 28 from 0f679202 and noted m = 27 lower by 4e-5 bits); generator of length (N + 21)(1/m + 1/n) + 2:
**<= 2^107.208650296 units**.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent, so the relation count is one fixed number; H1 puts it at
>= 2^63.6232 (pool G, band charge deducted), 4.3 SE above 1.05 U; the descent uses 2^63.575, two SE below that. Under H3 (r unramified in K_f and K_g, checked by [13]; r divides neither class number; the 21 sigma-odd
Schirokauer maps are injective on units modulo r-th powers; kernel of A = the line of virtual logarithms) Block
Wiedemann fails with probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel, degree <= 2N, entries uniform on >= 2^128
values); a failed solve (checked on the dropped rows) is a failure.

**Reconstruction.** Purged ideals' logs are recovered by back-substitution; the descent treats every ideal with at
most two occurrences as log-less (section 6); [6b]: the model's count of log-less ideals is >= 2.4x the measured one
in every lambda band with a loss.

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6 lattices, then special-q nodes sieving
an exact lattice ball (`dsieve16`) and recursing on child ideals; counted cells base-primitive (0.9100), retention >=
0.8592 at every sieve node ([7]), `ecm16.py` cofactoring ([12]), negative-binomial relation counts (toy-fit shapes),
log-less ideals discarded or re-descended in a fresh lattice, sub-bin yields capped at 4x (1.5x at the boot) the
sextile pr. The plan (62 nodes + boot) is recover-Q v7's, re-searched under this model; min mu 44.7, max phi 0.0268,
max psi 0.536 (boot). **Lazy descent:** the first relation whose subtree succeeds ends a node (`lazy_fraction`; [7l] simulates the three
costliest nodes: 0.0227/0.0161/0.0368 vs bounds 0.0243/0.0173/0.0420); the boot descends one relation at a time, at
most 18 (2.153 expected relations, 2^56.68 trials of 2^67.02 instructions, pr 2^-54.5); phi, psi and E[N] are the
eager recurrence's.

Expected two-log work **<= 2^112.4260 units** ([7]; unrounded 2^112.42596; one log 2^125.714 instructions). The
meter caps both logs at **2^113.278864326817 units** (cap reached = failure): P(cap) <= E[W]/cap <=
2^(112.4260 - 113.278864326817) <= 0.553684359093004 (Markov); with 10% more expected work, <= 0.6090527949. The meter is the least value its own success budget supports: success at +10% expected work, (1 - 2^-66) x
(1 - 2 phi - 1.1 x markov - P(abort) - P(repeat) - 2/r) >= 0.39 + 2^-48, binds exactly at
markov = 0.553684359093004, i.e. cap = 112.426 + 0.852864326817 = 113.278864326817 (`audit_refinement.py`); P18's
2^113.279 was this constraint rounded up.

## 7. Success

Core: success >= (1 - 2^-66) x (1 - 2 x 1.307e-05 - 0.553684359093004 - 2^-11.61 - 2^-10.70 - 2/r) >= 0.4453684
(audit >= 0.445368435909; model >= 0.4453848); with the Markov term at 10% more expected descent work (0.6090527949)
it is still >= 0.39000 (audit >= 0.3900000000002): main phase, two boots, the meter's Markov bound, the node counter
and special-q repeats.

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX = 2^27.3; E[N] uses the same
  recurrence with unit node cost (`counts16.py`): E[N] <= 2^15.6813 for both logs, P(abort) <= 2^-11.61.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for the 1/p weighting:
  P <= 2^(2 x 27.3 - 1 + 6 - 70.30) <= 2^-10.70, charged.

No lower bound on node work is used. Subtrees are treated as independent (H5). With the wrapper (section 1):
P(success) >= 0.4453684 - 2^-48 > 0.39; every core failure, abort and wrapper failure counts as a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (line sieve from 2^14, region Y, clean-run rec7, scan2, wide path, pool G) | <= 118.895525055 |
| linear algebra (weight 301, N <= 2^58.6463, seven owners, Horner chains, (27, 5)) | <= 116.621741985 |
| matrix generator | <= 107.208650296 |
| both individual logarithms (metered cap; expected 112.4260; success >= 0.39 at +10%) | 113.278864326817 |
| canonical a list | 91 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 72 |
| wide factor-base lattices ([2^63, 2^64), KERNELS.md) | 70.6 |
| rec7 bases (2 per ideal, p in [2^52, 2^64), KERNELS.md) | 68.3 |
| numerical certificate construction | 50 |
| verification | 11 |
| wrapper: pairing, field map, scalar multiplication, check (section 1, <= 2^80 instructions) | 66 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^91 (2^41.8 per eligible a, 2^63.451 of them:
2^90.963 rounded up; P10c's 2^90 row was under-rounded at 2^90.016); Schirokauer maps 2^88 (2^38.5 per relation
against < 2^34.5); factor base 2^72 (2^23.4 per prime against < 2^19); filtering and back-substitution 2^80 (2^24.4
per entry against < 2^14). The upper-rounded rows sum to at most 2^119.19121380943, so the core meter is
**2^119.1912138095** (audit-checked, as is the stated sum).

**Total.** The core meter 2^119.1912138095 bounds every core row above; with the wrapper row, W <= 2^119.1912138095
+ 2^66 <= 2^119.191213809501 (80-digit decimals, `audit_refinement.py`).

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) against section 2's computed values; exact sampled norms; sub-cell
maxima from sampled norms, their effect bounded by the caps). H1' retention (Euler-product model calibrated on exact smooth integers; the band-content excess of [2^8, 2^14)
charged, 0.009 bits). H2 region (uniform pools and region, exact flags, certified counts; sized on S, judged once on the fresh G; the
S+G work mean charged with the 0.82% allowance, a modeling margin). H3 matrix (r unramified, [13]; r not dividing the class numbers, the unit condition and a one-dimensional kernel
assumed; signed filtering on the production occurrence profile at excess 1.05, purge target 0.0001; the expected
N(300)/U not rising above the charged 0.03577 (mean + 2 SE at 2^19, 2^20, PREREG-H; G's means no higher, [10y]) at
production). H4 kernel cost (KERNELS.md, STORAGE.md: static counts of loop-free code, executed counts for routing); H4u (>= 87%
of production row weight is +-1; 0.8809..0.9366 at G's means, `unit_entries.py`); the wide path's per-hit price is
0bbca0b2's (same walk code); rec7's per-hit price depends on h' only; contaminated runs hold at most 0.3% of the
canonical a (the a-uniform sample's 95% bound; a shortfall changes cost only); f/g independence given q_f (charged
at 100x); band-strip ECM at the `ecm16` per-prime rate.
H5 descent (DESCENT.md): lower sums; lazy order independent of relation outcomes; the negative-binomial recurrence
with independent subtrees; missing logs (Poisson occurrence model at the 2-SE input; log-less ideals re-descended in
fresh lattices L_(Q q') or discarded); local omitted-mass laws; random-integer ECM orders; UDUP = 0.65
(relation-weighted E[s/m_true] over COMPLETE unit orbits over the CHARGED pr; live [11] 0.896 / 2.05).

## 10. Memory

Core storage <= 2^105 bytes (STORAGE.md); wrapper maps and polynomials below 2^30 bytes: below 2^106 bytes in all. All
preprocessing is performed for this target and charged once in W.

## Provenance

GT core: recover-Q v6 9332a74d (mpjunior92: lazy descent, five-vector Krylov passes, weight-300 merge); recover-Q 3486094d,
our v7 4c817536 and v8 4199f5f2; lineage 9ff4fc85, 57598408, 90d99567, 2c2ce2ab (gbotrel), 6282ab79 (bernardo-lima12),
8cc93834/a31a5f1d, 42e80f80, mesh 3a2e69ca/e3940d11; after Gautam Botrel, Guillevic-Singh and CADO-NFS. Ours: reduction
and field map 871a7c6c, 6f384591, ed51fb3d; ae060b47; P4 2899d54f; WIDE1 on 0bbca0b2; plain P5-P10b (line table of 734a9723,
per-entry layout of the overhead lane); rec7 and runs16.py (big-lever and core lanes); scan2 / apply / recover (overhead
lane, Muse Code, 0fb4319a); re-optimisation, pools S/G, band study, contamination sample and descent hardening: a05affb1
(REOPT-2000 lane, Claude Opus 5.5 in Claude Code via T3 Code), clean-run rule after 5aad7567. Filtering: 5170c7d8
(PREREG-G) and our recover-Q lane's PREREG-H runs (q13 7182431e); seven owners from 6734bd5a (mpjunior92) via 0f679202
(REOPT lane); Horner chains: Claude Opus 5.5 in Claude Code (yukon automining). P19 (this package): Block Wiedemann m = 27
and the meter cap 113.278864326817, re-optimised inside the package's own formulas. recover-Q wrapper (section 1, 2^66
row, audit): Claude Opus 5.5 in Claude Code via T3 Code; reduction as in our 9ac7a9b4 / 0fbd40b3.
