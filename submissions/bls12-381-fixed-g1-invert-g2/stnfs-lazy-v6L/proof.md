# STNFS with lazy descent and five-vector Krylov passes: 2^120.2457 reference units

Track bls12-381-fixed-g1-invert-g2. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used.
Predecessor work is reused as credited below; no prior verdict is a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied
nonidentity P in G1 (Q = infinity for z = 1) with success at least 0.39 over
uniform full-production challenges and fresh coins under the declared heuristics
H1-H5; the modelled success is at least 0.3901 (model 0.3902123; audit with
rounded inputs 0.3901785). A global instruction meter bounds all work, including
preprocessing, failures, restarts, transfers, verification and output, by
**2^120.2457 units**; the score is **120.2457**. Upper-rounded ledger rows
(section 8) sum to 2^120.2456705 <= 2^120.2457 (unrounded model total
2^120.24563). Storage is at most 2^105 bytes (STORAGE.md). The claim is
analytic and exploratory; no full-size run is made.

**What v6 changes (base: v5 9ff4fc85, 120.3702).** Relation collection,
filtering inputs, polynomials, descent plan and boot parameters are unchanged.
(a) **Lazy descent** (section 6, DESCENT.md): a sieve node cofactors its
survivors in random order and stops at the first relation whose subtree
succeeds, and the boot finds and descends one relation at a time (at most 18)
instead of collecting all 18 first. Expected two-log work 2^116.4941 ->
**2^115.6941**, meter 2^117.221 -> **2^116.421** at the same Markov budget.
(b) **Linear algebra** (section 4, STORAGE.md): the five Krylov sequences share
each pass over the matrix (`krylov5.s`), the merge continues to weight 300
(N(300)/N(100) <= 0.711 measured at 2^12..2^18), eight owners, blocks (28, 5), and
(v6L) block selection and sign once per entry per pass: 2^117.6351 -> **2^117.2564**.

**Numerical conventions.** Decimals are the experiment's values rounded as
shown; bounds are rounded in the safe direction and `declared_bounds` compares
every stated bound with the unrounded value. The model is `model16.py` and
`descent16.py` with `config.json` and `stats16.b85`; `evidence16.py` recomputes
the ledger and the checks cited as [n].

## 1. Reduction

Let g = e(P, Q0) (pinned `pairing(Q0, P)`), one pairing; Q0 is the fixed G2
generator. P has prime order r and e is nondegenerate, so g generates mu_r,
z = g^a for unique a mod r (a = 0 iff z = 1), and Q = [a]Q0 (one G2 scalar
mult), so e(P, Q) = z for that exact P (never replaced). log z and log g share
one relation/LA phase; a = log z / log g mod r.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks that
(u y)^((p^6-1)/2) = -1 in F_p6, so g is irreducible over F_p6 and
F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. The map
x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P), so
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at
(5.15, 0.85) in collection and descent; [14] recomputes it from the exact
Kummer-Dedekind structure: (4.868096, 0.843643) over p <= 2000, (5.085035,
0.731543) over p < 50000 (g side with the primitive-correct p = 2 double-root
term and q/(q+1) factor); margins 0.065, 0.118.

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21. The pool is uniform on theta = 5/16 (a5 > 0,
Q_a(a) <= (5/16) A^2); the collected region is its sub-region theta = 2330/8000
(`production_theta_region`): a5 > 0, Q_a(a) <= (2330/8000) A^2; half-slice
b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The restriction of the uniform sample to
qa <= 2330/8000 is uniform on the sub-region: 16033 canonical sub-region rows
(14709 primitive) of 19765, scaled by the sub-region volume over the effective
sample 164000 x (2330/2500)^3 = 132767.441. The canonical a of each unit class is
the least of its exact HNF key class in (Q_a, lex) order, so each relation orbit
of <units, sigma> has one point in the domain; non-primitive pairs are excluded
by an exact HNF test. `lattice_count_bounds.py` certifies (exact rational
Machin-pi, determinant/sqrt intervals) the count ratio at every radius
R >= 781 = isqrt((2330/8000) 2^21) (b-radii >= 1219): [0.983851, 1.015464] for
Q_a and [0.983798, 1.015464] for Q_b, inside [0.98, 1.02] (asserted); the model
uses the exact products **0.967911517 (yield, floor)** and **1.031166365
(work, ceiling)**. Pool sampling error is handled by the canonical-work
allowance below. The canonical list costs 2^90 units.

**Pool.** 164,000 uniform eligible pairs on theta = 5/16 from seeded indices
(seed 5: 2,000, seed 6: 162,000; `pool12.b85`), 19,765 canonical rows (12.05%);
`load_pool` regenerates every row from its (seed, index) (exact norms, HNF
primitivity, certified nbound statistics); [2] checks list membership both ways
and a 64-row anchor table. Every draw was tested (no selection; one draw needs a
10^6-node orbit budget); seeds and draw count were fixed before any yield was
computed.

**Sieve, thresholds, cofactorisation.** Line sieve (scores only) for p in [2^8, 2^48), record
walk for p in [2^48, 2^64) (`slice12`, KERNELS.md; X2 = 2^64), primes in [2^64, 2^69] taken by ECM
within the slacks (295, 226). **Line sieve to 2^48 with band strip (v3; 8cc93834 /
42e80f80 lineage):** primes in (2^20, 2^48] only add scores on the line path
(every line keeps >= 2^7.11 hits per coset at p = 2^48, inside the measured line
table) and emit no records: records/cell 2.33 -> 0.58. Their factors are recovered at the
survivors by a band strip: two rounds per side of the `ecm16.py` 48-bit schedule at
99.9% on the unstripped cofactors (caps 832 / 576 bits: 2 x 2^36.4488 f side,
2 x 2^35.3903 g side). A band prime missed in both rounds loses the relation: yield x
(1 - 41 x 0.001^2)(1 - 28 x 0.001^2) (at most floor(832/20) = 41 and 28 band
factors per side). The f-side slack is 295, so survivors'
f-side cofactors are bounded by 303 bits; the g-side slack is 226, so survivors'
g-side cofactors are bounded by 234 bits (the slacks are re-optimised for the
sub-region, section 8). The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at **99.9%** per prime (ECM_FIND = 0.999), which [12] checks: the cheapest
grid point for 69-bit primes in a 303-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves,
five-limb Montgomery multiplication at 367.76 Ir, `ecm12`) prices one side at 2^38.0906
instructions and the other (234 bits) at 2^38.0879. **Early abort (2c2ce2ab):** the g-side band strip and ECM run
only after the f side factors completely (yield unchanged). Per survivor: C_NORMS + C_TDIV + 2 x 2^36.4488 + 2^38.0906 +
q_f (2 x 2^35.3903 + 2^38.0879), q_f = min(1, 100 rho(u_f)/P(omitted_f <= slack_f))
per pool row (heuristic under H1, mean ~1e-5, factor 100 a margin). The 99.9% is the same
ECM_FIND used in the yield's miss model.
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20, the band strip and the ECM
schedule above (f side first, g side only after a complete f factorisation).

**Certified block bounds.** F_v Taylor-shifted to the block centre per place,
remainder over the block disk plus a 2^-32 relative term, rounded outward;
`nbound31` replays the Python leaves bit for bit (v13 pool, 17,874 rows) at
1,508.2 Ir per evaluation (charged as measured); [2] checks bounds against exact
norms. Every block of every row is enumerated, sieved and charged (no selection).

**Yield ([4]).** rels = (va_region/132767.441) vb (0.967911517) sum over canonical primitive sub-region
rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
the ECM miss and the band-strip loss. **Central core (v4):** relations from
cells with qa = Q_a(a)/A^2 < Q_CORE = 0.02 are discarded (never passed to
filtering, linear algebra or the descent; those cells are still sieved and
charged), so count and relation-weighted mean norms are those the estimator
measures over the rows with qa >= 0.02; the discard only lowers the count. Retained relations >= 2^63.6617
against 1.1 U = 2^63.589 (excess 1.1, section 4): margin >= 0.0728 bits (unrounded
0.07289; pool SE 0.0364, so 2.00 SE). Feasibility requires a margin of at least
max(`MARGIN` = 0.06, `MARGIN_SE` = 2 SE). The descent's missing-log input `LOG2_RELS =
63.661` is the floor (to 0.001 bit) of this retained count (unrounded
2^63.66172; certified lower product already applied). [4] checks `LOG2_RELS <=
log2_rels` (`rels_input` check). The mean-norm inputs 739.9/498.9 are rounded
down from 739.98/498.95 and checked in [4].

**Cost ([5]).** Sampled mean 306.3 instructions/cell; charged mean 311.6
with the allowance below, over 2^125.934 cells: all-cell bounds, layout,
setups, line path and record walk at the re-measured recentred `slice12` costs
(KERNELS.md; primes in [2^63, 2^64) charged twice), routing
(10.01/record, 8/cell; STORAGE.md), apply 7.00, recovery 7.003/record, scan,
zero, band strip and f-first cofactorisation. Collection: **<= 2^119.9297 units** (estimate
2^119.9048 times allowance 1.0174).

**Canonical-work allowance (H2, [4a]).** The charged work is (va_region/132767.441)
times the per-row sum over in-region canonical rows: a Monte Carlo estimate of
a fixed, challenge-independent domain total. The estimator sums all 164,000
draws (zero outside the sub-region or non-canonical), so its relative SE has a
conditional row-mean term and an incidence term for the 16033 contributing
draws: conditional 0.434%, incidence sqrt((1-16033/164000)/16033) = 0.750%, total
<= 0.867%. The charge multiplies the estimate by **1.0174 >= 1 + 2.00 x 0.00867
= 1.01734** (rounded up; unrounded need 1.0173370); [4a] fails if the allowance falls below 1 + 2.00 SE.
This is a declared modeling margin under H2, not a probability term, and no
sampling tail is charged against success.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge (coefficients <= 255, including a
collision-safe bound on implicit middle-range entries), simulated on the
production occurrence profile (`filtsim.py`; [6] reruns 2^12) as the signed
sigma-quotient system. **Required excess 1.1.** The whole region is sieved and
every retained relation (>= 2^63.6617 = 1.157 U) is filtered and used; the
matrix is charged at excess 1.1, which over-states N at the larger actual
excess. Data are measured at excess 1.1 with the mean norms 739.9/498.9. Four
seeds at 2^17 give N/U <= 0.05797 at weight 100; N/U falls with scale (max over
seeds): 0.1174 (2^12), 0.0748, 0.0675, 0.0618, 0.0592, 0.0580, 0.0529 (2^18)
(`SCALE_NU`, [10] checks non-increasing); pad 1.066 covers every value from 2^15
up. **Merge weight 300 (v6).** The merge continues to mean weight 300 (charged
301 per row); N(300)/N(100) maxima per scale are 0.7098 (2^12, 6 seeds), 0.7034
(2^13, 6), 0.7043 (2^14, 3), 0.6971 (2^15, 3), 0.6997 (2^16, 3), 0.6897 (2^17,
2), 0.6913 (2^18, 1), charged **0.711** (`W300`, [10] checks max <= 0.711); the
weight-200 ratio was 0.776-0.795 (v5 charge 0.80). The coefficient cap never
binds (largest explicit or implicit bound 192 < 255 at weight 300). N = 0.711 x
1.066 x 0.05797 U = **2^58.9429** < 8 x 2^56. N/U and the ratio not rising at
production scale are part of H3.

The system is the square A = [M | S]: N + 21 rows, N ideal columns and 21 dense
Schirokauer columns (sigma-odd unit ranks 18 + 3 from signatures K_f (0, 36),
K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3 its kernel is the line of
virtual logs. Block Wiedemann: m = 30 left vectors X with entries uniform in
[0, 2^128), n = 5 right vectors Y uniform in F_r. The Krylov phase takes N/m +
N/n block steps, N(n/m + 1) vector-products.

- **Krylov products (v6, STORAGE.md).** The five sequences advance in lockstep;
  element j of the five current vectors shares one 256-byte slot, and one pass
  over a row serves all five: `unit_run5` handles a +-1 entry for five vectors
  in 39 instructions, `gen_run1` a coefficient entry in 33 per vector (clang 18
  listings, counted by [10]). Block selection (4) and sign (3) are per entry, once
  per pass: per vector-entry (4 + 3)/5 + 0.88 x 39/5 + 0.12 x 33 = 12.22 -> **12.3** (>=
  88% +-1 by H4u, measured 0.905-0.948 at weight 300, 2^12..2^15). Per row and
  vector: 301 x 12.3 + C_ROW + 21 x 97.
- **Row tail.** Eight active 2^56-element owners: seven sends and seven
  combinations, `_finish`, headers, routing: counted 1,424 by [10]
  (`owner10.row_tail`), charged C_ROW = 1,430.
- **Horner and projections.** Horner (w <- A w + Y F_i, deg F <= N/m + N/n) uses
  single-vector products at 301 x 20.7 + 1,430 + 21 x 97 plus n dense terms and
  one reduction (420) per element. Projections: (m + n) N^2 terms of a 128-bit X
  entry times a field element at 60 (`bw_dot2`, counted 59) plus a reduction
  per entry.

With c_K and c the Krylov and single-vector row costs: N^2 (c_K (n/m + 1) +
(1/m + 1/n)(c + 97 n + 420) + 60 (m + n)) + 420 (m + n) N instructions =
**<= 2^117.2564 units** at N = 2^58.9429 ((28, 5) minimises the formula). The
(m + n) x (m + n) generator of the length L = (N + 21)(1/m + 1/n) + 2 sequence:
<= 4 (m + n)^3 L^log2(3) log2(L) field mults at 97 + 420: **<= 2^107.8065 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** Region and polynomials are coin/challenge-independent, so the
  count is one fixed number; H1 puts it at >= 2^63.6617, >= 0.0728 bits above
  1.1 U (pool SE 0.0364, so 2.00 SE above requirement); all of them are
  filtered and the descent uses 63.661, at or below that count.
- **Linear algebra.** H3 (standard NFS-DL): r is unramified in K_f and K_g
  (checked by [13]), divides neither class number, the 21 sigma-odd Schirokauer
  maps are injective on the units modulo r-th powers, and the kernel of A is the
  line of virtual logarithms. Block Wiedemann then fails with probability
  <= 4N/2^128 < 2^-66 (Schwartz-Zippel on the generator's determinant conditions,
  degree <= 2N in the block entries, each entry uniform on a set of size >= 2^128).
  A failed solve (checked on the dropped rows) counts as failure.
- **Reconstruction.** Purged ideals' logs are recovered by back-substitution.
  The descent counts every ideal with <= 2 occurrences as log-less; [6b]
  measures ideals left without a log on the [6] system, and the model's count
  is >= 12.5x measured in every band with loss (DESCENT.md).

Under H1 and H3 the main phase succeeds with probability >= 1 - 2^-66 over the coins.

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6
lattices, then special-q nodes sieving an exact lattice ball (`dsieve16`) and
recursing on child ideals. Counted cells are base-primitive (0.9100);
survivor-test retention is >= 0.8593 at every sieve node ([7]); cofactors use
the `ecm16.py` schedule ([12]). Relation counts are negative binomial (toy-fit
shapes). Log-less ideals are discarded or re-descended in a fresh lattice. The
plan (62 nodes + boot (215, 344, 14, 4)) is v5's.

**Lazy descent (v6; DESCENT.md).** The first relation whose subtree succeeds
ends a node, so per-relation work is done on demand. (i) After the full sieve, a
node cofactors its s survivors in uniformly random order, descending each distinct
relation at once: E[processed]/s <= E_M min(1, (1 + 1/s)/((1 - psi)(M + 1)) +
psi^M), M ~ NB(mu, r) as in the failure model (`lazy_fraction`; [7l] simulates the
process at the three costliest nodes: 0.0596/0.0405/0.0878 vs bounds
0.0627/0.0449/0.0948). (ii) The boot finds and descends one relation at a time,
at most 18: by Wald (1 - psi^18)/(1 - psi) = 3.31 expected relations (v5: 18),
trials 2^60.69 -> 2^58.25. Failure events, phi, psi, E[N] and repeats are v5's.

Expected two-log work **<= 2^115.6941 units** ([7]; one log 2^128.982
instructions: boot own 2^126.20 + subtree 2^127.02 x 3.318). The meter caps both
logs at **2^116.421 units** (cap reached = failure): P(cap) <= E[W]/cap <=
2^(115.6941 - 116.421) <= 0.6042008 (Markov).

## 7. Success

Success >= (1 - 2^-66) x (1 - 2 x 1.573e-03 - 0.6042008 - 2^-9.06 - 2^-10.70 - 2/r)
>= **0.3901** >= 0.39 (audit 0.3901785; model 0.3902123): main phase, two boots
(phi <= 1.573e-03 each), the instruction-meter Markov bound, the node-counter
Markov bound and special-q repeats.

- Instruction meter: cap 2^116.421 units, P(cap) <= 0.6042008 (section 6).
- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX
  = 2^27.3. E[N] uses the same recurrence with unit node cost (`counts16.py`):
  E[N] <= 2^18.24 for both logs, P(abort) <= 2^-9.06 <= 0.0018734.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for
  1/p weighting: P <= 2^(2 x 27.3 - 1 + 6 - 70.30) <= 2^-10.70, charged.

No lower bound on node work is used. Subtrees are treated as independent (H5).
Targets are uniform nonidentity by power randomisation (2/r charged); z = 1 is
handled directly; the final check e(P, Q) = z rejects any incorrect Q.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (excess 1.1, band strip, 164k pool; estimate 119.9048 + 1.74% canonical-work allowance) | <= 119.9297 |
| linear algebra (weight 301, N = 2^58.9429, eight owners, five-vector Krylov passes, blocks (28, 5)) | <= 117.2564 |
| matrix generator | <= 107.8065 |
| both individual logarithms (metered cap; expected 115.6941) | 116.421 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

**Fixed rows (upper bounds; units x 2^14.29 = instructions).** Canonical list
2^90: 2^41.8 instructions per eligible a (2^62.5 of them) against one LLL-reduced
HNF key, the `canonical_a` orbit enumeration and a segmented key sort (< 2^30).
Schirokauer maps 2^88: 2^38.5 per relation (<= 2^63.8) against one exponentiation
per side mod r^2 (< 2^34.5). Factor base 2^72: 2^23.4 per prime p < 2^69 against
root finding (< 2^19). Filtering, purge, merge bookkeeping and back-substitution
2^80: 2^24.4 per relation entry (2^69.9) against < 2^14. Certificates 2^50;
pairing, [a]Q0, check 2^11.

The upper-rounded rows sum to 2^120.2456705, so W <= **2^120.2457**
(audit-checked).

**Configuration.** Excess 1.1; pool theta 5/16 (164,000 draws); region theta
2330/8000, Q_CORE 0.02, radius 781, allowance 1.0174, X2 = 2^64, slacks (295,
226), alpha (5.15, 0.85); descent inputs 63.661, 739.9/498.9; plan and boot as
v5, UDUP 0.65, THIN_MAX 85, BOOT_RELS 18, lazy nodes and boot; meter 2^116.421,
node counter 2^27.3. LA: weight 301, ratio 0.711, pad 1.066, C_ROW 1,430 (eight
owners), 12.3 (Krylov) / 20.7 (Horner) per entry, 97 per dense term, 420 per
reduction, blocks (28, 5) with 128-bit X entries.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) in collection and descent against the computed values of section 2 ([14]); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool and sub-region, yield lower bound without the qa < 0.02 core, exact flags,
certified counts, and the canonical work mean charged with the 1.74% allowance,
a modeling margin, not a probability term); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
signed filtering on the production occurrence profile at relation excess 1.1,
pad 1.066 over the 2^17 maximum, N/U (max over seeds) non-increasing 2^12..2^18
and assumed not to rise at production; the weight-300 ratio N(300)/N(100) <= 0.711
not rising at production scale is part of H3, measured 0.690-0.710 at
2^12..2^18 at excess 1.1); H4 kernel cost (KERNELS.md, STORAGE.md; LA tail/reductions are
statically counted loop-free code, the route an executed-instruction count),
with H4u: >= 88% of production row weight is +-1 (measured 0.905..0.949 at
weights 100-300, 2^12..2^15, `unit_entries.py`; a shortfall changes cost only); f/g independence given
survival (H1) bounds the early-abort probability q_f (charged at 100x the
estimate); band-strip ECM finds band primes at the `ecm16` per-prime rate (random
group orders, as in collection);
H5 descent (DESCENT.md: lower sums, lazy order independent of relation outcomes, negative-binomial recurrence, missing logs,
independent subtrees, fresh-lattice re-descents, local omitted-mass laws,
random-integer ECM orders, UDUP = 0.65 (relation-weighted E[shat/m]/pr >=
1.06 boot / >= 4.70 nodes, uniform cells, at v13's node regimes; the boot lattice
(R = 4) and the binding (220, g) node are re-checked live by [11] at the shipped
plan (boot kA 14: 1.27); spot cells at re-optimised nodes in DESCENT.md)).

## 10. Memory

Storage <= 2^105 bytes (STORAGE.md).

## Provenance

Algorithm, models, kernels: Gautam Botrel (Claude Code, Codex); polynomials
and alpha after Guillevic-Singh, filtering after CADO-NFS. Refinements:
kaankolcu, pepedesigner, mpjunior92, DeepSeek and muse-spark-1.3-contributor
lineages (v13 57598408); 2c2ce2ab (gbotrel), 6282ab79 (bernardo-lima12); band
strip 8cc93834/a31a5f1d (kaankolcu), 42e80f80 (gbotrel); v1-v5 (mpjunior92,
Claude Opus 5.5) ending in 9ff4fc85 (120.3702). v6 (lazy descent, five-vector
Krylov passes, weight-300 merge, eight owners): Claude Opus 5.5 in Claude Code.
No prior review verdict is a premise.
