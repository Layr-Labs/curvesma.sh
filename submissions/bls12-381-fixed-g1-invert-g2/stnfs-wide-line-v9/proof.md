# STNFS v9: measured line-path table, wide-prime kernels and a freshly judged region on v8: 2^119.889929619243 reference units

Track bls12-381-fixed-g1-invert-g2. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used.
Predecessor work is reused as credited below; no prior verdict is a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied
nonidentity P in G1 (Q = infinity for z = 1) with success at least 0.39 over
uniform full-production challenges and fresh coins under the declared heuristics
H1-H5 (model 0.3902307; audit with rounded inputs 0.3902000). A global
instruction meter bounds all work, including preprocessing, failures, restarts,
transfers, verification and output, by **2^119.889929619243 units**; the score is
**119.889929619243**. Upper-rounded ledger rows (section 8) sum to
2^119.889929619242911 <= 2^119.889929619243. Storage is at most 2^105 bytes (STORAGE.md). The claim is
analytic and exploratory; no full-size run is made.

**What v9 changes (base: v8 4199f5f2, valid 119.933288883129).** (a) **Line path** (KERNELS.md): v8 charged 24.10 Ir per
update for every h >= 2^13.82, while most line work runs at h = 2^20-2^48; v9 charges the measured public `slice12`
line walk on random (p, s) pairs (cap 17.2 for h >= 2^18.81; below 2^13.82 the larger of v8's points and the
random-pair bounds). (b) **Primes in [2^63, 2^64)** run the measured `wide12.c` path instead of v8's 2x cushion. (c)
**Region judged on a fresh pool** (section 3): pool A over-stated v8's yield (an independent 656k pool found 0.05 bits
fewer relations), so the region is sized on pool B at excess 1.05 and judged once on a fresh pool D (pre-registered);
work is priced on B, C, D at v8's allowance. (d) Filtering re-measured at excess 1.05 (unit share 0.87),
ten matrix owners. (e) Descent inputs from D (estimate - 2 SE). Collection 2^119.669318937 -> **2^119.444013904**,
LA 2^117.279842877 -> **2^117.925074783**, meter 2^112.975291 -> **2^113.156048**. v8's box bounds and v7's descent
and LA kernels are kept.

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
Q_a(a) <= (5/16) A^2); the collected region is its sub-region theta = 2326/8000
(`production_theta_region`): a5 > 0, Q_a(a) <= (2326/8000) A^2; half-slice
b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The restriction of a uniform pool to qa <= 2326/8000 is uniform on the sub-region;
each pool is scaled by the sub-region volume over its effective sample 656000 x (2326/2500)^3 = 528339.334. The
canonical a of each unit class is the least of its exact HNF key class in (Q_a, lex) order, so each relation orbit of
<units, sigma> has one point in the domain; non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` certifies (exact rational Machin-pi, determinant/sqrt intervals) the count ratio at every
radius R >= 780 = isqrt((2326/8000) 2^21): [0.983830, 1.015484] for Q_a and [0.983778, 1.015484] for Q_b inside [0.98, 1.02] (asserted); the
model uses the exact products **0.967870931 (yield, floor)** and **1.031206887 (work, ceiling)**.
The canonical list costs 2^90 units.

**Pools (v9; `pool_d.json`).** Every pool draws 164,000 uniform eligible pairs per seed on theta = 5/16 from seeded
indices of the package sampler; every draw is tested with the exact canonical predicate and every canonical row
kept (no selection). B = seeds 7-10, C = 11-14, D = 41-44 (656,000 draws each). Pre-registered before any statistic of
D (PREREG-D): theta is the smallest 1/8000 step with margin_B >= 0.06 + sqrt(2) SE_B at excess 1.05 (2326: margin
0.086607 >= 0.086070; 2325: 0.085747 < 0.086084, [3s]); D judges it once
(margin_D >= max(0.06, 2 SE_D)); the work is estimated on B, C and D; v8's pool A is no longer used. `pool_d.json`
ships each seed's region sums and index-list hash; [2] regenerates D seed 41 draws [0, 1000) and must reproduce that
block's sums. D: 78,782 canonical rows, 63,651 in the region.

**Sieve, thresholds, cofactorisation.** Line sieve (scores only) for p in [2^8, 2^48), record
walk for p in [2^48, 2^64) (KERNELS.md: public `slice12` line walk; record walk `slice12g` below
2^63 and `slice12` + `wide12.c` in [2^63, 2^64); X2 = 2^64), primes in [2^64, 2^69] taken by ECM
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
remainder over the block disk plus a 2^-32 relative term, rounded outward; one
evaluation (`NormBounds.block`) is what `nbound31` replays bit for bit (v13 pool,
17,874 rows) at 1,508.2 Ir (charged as measured). **Box bounds (v8).** The rows of
one line (same a, b1, b5) whose (b2, b3, b4) lie in one aligned 2 x 2 x 4 box
(`model16.box_rows`, up to 16 rows) share one bisection per 64-cell block.
`EB.hull` gives per place the centre of the rows' context centres and a radius
>= |c_k - centre| + r_k for every row k, each step rounded outward: row k's b(v)
lies within r_k + r0 of c_k + m, hence within the hull radius + r0 of centre + m,
so one `block` on the hull bounds every row's cells. `EB.box_leaves` bisects the
union of the rows' b0 ranges until both bounds are within 8 bits; a single column
the hull leaves unresolved is evaluated with each overlapping row's own context,
and a row's cell is rejected only if that fails (no pool row is). Boxes x blocks
x (b1, b5) partition each slice (every box row keeps b5 >= 1, the sampled half),
so charging each cell C_NB_EVAL x its box's evaluations per box cell has
expectation total evaluations / cells over the uniform pool; each row is charged
every leaf it uses (C_SCAN_LEAF x leaves / blen) and its context plus a hull share
(C_NB_ROW 4,300 -> 4,600; the hull update is ~25 instructions per place and row).
Work-pool (B, C, D) region means: evaluations per cell 0.0073, leaves per cell 0.0332.
[2b] rebuilds six boxes of replayed D rows and checks 1,536 exact norms at sampled cells of every box row
against their leaf bounds; [2] replays every field of the replayed rows. Every block of every row is enumerated, sieved and
charged (no selection).

**Yield ([4], judged pool D).** rels = (va_region/528339.334) vb (0.967870931) sum over canonical primitive sub-region
rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
the ECM miss and the band-strip loss. **Central core (v4):** relations from
cells with qa = Q_a(a)/A^2 < Q_CORE = 0.02 are discarded (never passed to
filtering, linear algebra or the descent; those cells are still sieved and
charged), so count and relation-weighted mean norms are those the estimator
measures over the rows with qa >= 0.02; the discard only lowers the count. Retained relations >= 2^63.6016
against 1.05 U = 2^63.5218 (excess 1.05, section 4): margin >= 0.0798 bits (pool SE 0.0177,
so 4.5 SE). Feasibility requires a margin of at least max(`MARGIN` = 0.06, `MARGIN_SE` = 2 SE). The
descent's missing-log input `LOG2_RELS = 63.566` is the judged estimate minus 2 SE, floored to 0.001 bit (unrounded
2^63.56621; red line (a)); [4] checks it. The mean-norm inputs 740.8/499.4 are D's
relation-weighted means 740.86/499.43 rounded down, checked in [4].

**Cost ([5], work pools B, C, D).** Sampled mean 221.0 instructions/cell (v8 255.7); charged mean
224.8 with the allowance below, over 2^125.919 cells: all-cell bounds, layout, setups, line path
at the v9 table, record walk (`slice12g`; wide primes at `REC_TAB_W` x 1.01), routing (10.01/record, 8/cell;
STORAGE.md), apply 7.00, recovery 7.003/record, scan, zero, band strip and f-first cofactorisation. Collection:
**<= 2^119.444013904 units** (estimate 2^119.4191 times allowance 1.0174).

**Canonical-work allowance (H2, [4a]).** The charged work is (va_region/n_eff) times the per-row sum over in-region
canonical rows of B, C and D: a Monte Carlo estimate of a fixed, challenge-independent domain total over 1,968,000
draws (190,262 contributing): total relative SE 0.232%. The charge keeps v8's **1.0174** as a floor ([4a]
checks it is >= max(1.0174, 1 + 2 SE = 1.00464)); a declared modeling margin under H2, not a probability term.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge (coefficients <= 255, including a
collision-safe bound on implicit middle-range entries), simulated on the
production occurrence profile (`filtsim.py`; [6] reruns 2^12) as the signed
sigma-quotient system. **Required excess 1.05.** The whole region is sieved and every retained relation
(>= 2^63.6016) is filtered and used; the matrix is charged at excess 1.05, which over-states N at the
larger actual excess. Data are measured at excess 1.05 with D's mean norms 740.8/499.4 (purge target
0.002 U, v8's): 8 seeds at 2^17 give N/U <= 0.06617 at weight 100; the maxima over seeds by scale
are 0.1235 (2^12, 10 seeds), 0.0856 (2^13, 10 seeds), 0.0818 (2^14, 10 seeds), 0.0738 (2^15, 10 seeds), 0.0679 (2^16, 10 seeds), 0.0662 (2^17, 8 seeds), 0.0597 (2^18, 6 seeds) (`SCALE_NU`); pad 1.116 makes r17 x pad >= every maximum from 2^15 up ([10]). The maxima fall with scale ([10] checks non-increasing).
**Merge weight 300 (v6).** The merge continues to mean weight 300 (charged 301 per row); N(300)/N(100) maxima per
scale are 0.7318 (2^12), 0.7085 (2^13), 0.6932 (2^14), 0.6987 (2^15), 0.6910 (2^16), 0.6912 (2^17), 0.6864 (2^18), charged **0.732** (`W300`, [10] checks the maximum over every scale). The coefficient cap
255 is reached in 1 of 64 runs (largest bound per scale 160, 254, 246, 255, 230, 235, 206); the simulated merge enforces it,
so the measured N include its effect. N = 0.732 x 1.116 x 0.06617 U = **2^59.2420** < 10 x 2^56.
N/U and the ratio not rising at production scale are part of H3.

The system is the square A = [M | S]: N + 21 rows, N ideal columns and 21 dense
Schirokauer columns (sigma-odd unit ranks 18 + 3 from signatures K_f (0, 36),
K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3 its kernel is the line of
virtual logs. Block Wiedemann: m = 28 left vectors X with entries uniform in
[0, 2^128), n = 5 right vectors Y uniform in F_r. The Krylov phase takes N/m +
N/n block steps, N(n/m + 1) vector-products.

- **Krylov products (STORAGE.md).** The five sequences advance in lockstep;
  element j of the five current vectors shares one 256-byte slot, and one pass
  over a row serves all five: `unit_run5` handles a +-1 entry for five vectors
  in 39 instructions, `gen_run1` a coefficient entry in 33 per vector (clang 18
  listings, counted by [10]). The inherited 4 + 3 (block selection, sign) is
  charged once per +-1 entry (one pass for five vectors) and per vector for
  other entries: per vector-entry 0.87 x (4 + 3 + 39)/5 + 0.13 x (4 + 3 + 33) =
  13.204 -> **13.21** (>= 87% +-1 by H4u, measured 0.884-0.935 at 2^10 (20 seeds), 0.877-0.937 at 2^11 (20 seeds), 0.942-0.946 at 2^15 (4 seeds) at weight 300).
  Per row and vector: 301 x 13.21 + C_ROW + 21 x 97.
- **Row tail.** Ten active 2^56-element owners: nine sends and nine
  combinations, `_finish`, headers, routing: counted 1,638 by [10]
  (`owner10.row_tail`), charged C_ROW = 1,650.
- **Horner and projections.** Horner (w <- A w + Y F_i, deg F <= N/m + N/n) uses
  single-vector products at 301 x 20.9 + 1,650 + 21 x 97 plus n dense terms and
  one reduction (420) per element. Projections: (m + n) N^2 terms of a 128-bit X
  entry times a field element at 60 (`bw_dot2`, counted 59) plus a reduction
  per entry.

With c_K and c the Krylov and single-vector row costs: N^2 (c_K (n/m + 1) +
(1/m + 1/n)(c + 97 n + 420) + 60 (m + n)) + 420 (m + n) N instructions =
**<= 2^117.925074783 units** at N = 2^59.2420. The
(m + n) x (m + n) generator of the length L = (N + 21)(1/m + 1/n) + 2 sequence:
<= 4 (m + n)^3 L^log2(3) log2(L) field mults at 97 + 420: **<= 2^108.288165823 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** Region and polynomials are coin/challenge-independent, so the
  count is one fixed number; H1 puts it at >= 2^63.6016, >= 0.0798 bits above
  1.05 U (judged pool D, SE 0.0177); all of them are filtered and the descent uses 63.566, 2 SE below it.
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
  is >= 7.1x measured in every band with loss at excess 1.05 (DESCENT.md).

Under H1 and H3 the main phase succeeds with probability >= 1 - 2^-66 over the coins.

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6
lattices, then special-q nodes sieving an exact lattice ball (`dsieve16`) and
recursing on child ideals. Counted cells are base-primitive (0.9100);
survivor-test retention is >= 0.8592 at every sieve node ([7]); cofactors use
the `ecm16.py` schedule ([12]). Relation counts are negative binomial (toy-fit
shapes). Log-less ideals are discarded or re-descended in a fresh lattice. The
plan (62 nodes + boot (215, 328, 16, 4)) is re-searched at the v7 yields.

**Lazy descent (v6; DESCENT.md).** The first relation whose subtree succeeds
ends a node, so per-relation work is done on demand. (i) After the full sieve, a
node cofactors its s survivors in uniformly random order, descending each distinct
relation at once: E[processed]/s <= E_M min(1, (1 + 1/s)/((1 - psi)(M + 1)) +
psi^M), M ~ NB(mu, r) as in the failure model (`lazy_fraction`; [7l] simulates the
process at the three costliest nodes: 0.0217/0.0175/0.0462 vs bounds
0.0244/0.0173/0.0423, within the 0.01 Monte Carlo allowance; 40,000 runs: 0.0158, 0.0399). (ii) The boot finds and descends one relation at a time,
at most 18: by Wald (1 - psi^18)/(1 - psi) = 2.170 expected relations, 2^56.70
trials. Failure events are those of the recurrence (DESCENT.md).

Expected two-log work **<= 2^112.4402 units** ([7]; one log 2^125.728
instructions: boot own about 2^123.71 + subtree about 2^124.20 x 2.169). The
meter caps both logs at **2^113.156048 units** (cap reached = failure): P(cap) <=
E[W]/cap <= 2^(112.4402 - 113.156048) <= 0.6088472 (Markov; unrounded 0.6088165).

## 7. Success

Success >= (1 - 2^-66) x (1 - 2 x 1.475e-05 - 0.6088472 - 2^-11.60 - 2^-10.70 - 2/r)
>= **0.3902** >= 0.39 (audit 0.3902000; model 0.3902307): main phase, two boots
(phi <= 1.475e-05 each), the instruction-meter Markov bound, the node-counter
Markov bound and special-q repeats.

- Instruction meter: cap 2^113.156048 units, P(cap) <= 0.6088472 (section 6).
- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX
  = 2^27.3. E[N] uses the same recurrence with unit node cost (`counts16.py`):
  E[N] <= 2^15.70 for both logs, P(abort) <= 2^-11.60 <= 0.0003222.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for
  1/p weighting: P <= 2^(2 x 27.3 - 1 + 6 - 70.30) <= 2^-10.70, charged.

No lower bound on node work is used. Subtrees are treated as independent (H5).
Targets are uniform nonidentity by power randomisation (2/r charged); z = 1 is
handled directly; the final check e(P, Q) = z rejects any incorrect Q.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (excess 1.05, band strip, box bounds, v9 line table, wide-prime path; pools B, C, D; estimate 119.4191 + 1.74% allowance) | <= 119.444013904 |
| linear algebra (weight 301, N = 2^59.2420, ten owners, five-vector Krylov passes, blocks (28, 5)) | <= 117.925074783 |
| matrix generator | <= 108.288165823 |
| both individual logarithms (metered cap; expected 112.4402) | 113.156048 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |
| wide-prime factor-base lattices (`wide12.c`, KERNELS.md) | 70.6 |

**Fixed rows (upper bounds; units x 2^14.29 = instructions).** Canonical list
2^90: 2^41.8 instructions per eligible a (2^62.5 of them) against one LLL-reduced
HNF key, the `canonical_a` orbit enumeration and a segmented key sort (< 2^30).
Schirokauer maps 2^88: 2^38.5 per relation (<= 2^63.8) against one exponentiation
per side mod r^2 (< 2^34.5). Factor base 2^72: 2^23.4 per prime p < 2^69 against
root finding (< 2^19). Filtering, purge, merge bookkeeping and back-substitution
2^80: 2^24.4 per relation entry (2^69.9) against < 2^14. Certificates 2^50;
pairing, [a]Q0, check 2^11.

The upper-rounded rows sum to 2^119.889929619242911, so W <=
**2^119.889929619243** (audit-checked).

**Configuration.** Excess 1.05; pools B (sizing), D (judged), B, C, D (work), 164,000 draws per seed on theta
5/16; region theta 2326/8000, Q_CORE 0.02, radius 780, allowance 1.0174, X2 = 2^64, slacks (295, 226), alpha (5.15, 0.85);
line table (18.81, 17.2), (13.82, 24.1), (9.82, 38.99), (5.82, 88.32), (1.83, 352.1); wide primes `REC_TAB_W` x 1.01; descent inputs 63.566, 740.8/499.4; sub-bin yields (caps
4x/1.5x), plan and boot of v7, UDUP 0.65, THIN_MAX 85, BOOT_RELS 18, lazy nodes and boot; meter 2^113.156048, node counter
2^27.3. LA: weight 301, ratio 0.732, pad 1.116, C_ROW 1,650 (ten owners), 13.21 (Krylov) / 20.9 (Horner) per
entry, 97 per dense term, 420 per reduction, blocks (28, 5) with 128-bit X entries.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) in collection and descent against the computed values of section 2 ([14]); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pools and sub-region sized on B and judged once on a fresh pool D, yield lower bound without the qa < 0.02 core,
exact flags, certified counts, and the canonical work mean on B, C, D charged with the 1.74% allowance,
a modeling margin, not a probability term); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
signed filtering on the production occurrence profile at relation excess 1.05,
pad 1.116 over the 2^17 maximum covering every maximum from 2^15 up, N/U assumed not to rise at production
(maxima 2^12..2^18 in section 4); the weight-300 ratio N(300)/N(100) <= 0.732 not rising at production scale is
part of H3, measured 0.686-0.732 at 2^12..2^18 at excess 1.05); H4 kernel cost (KERNELS.md, STORAGE.md; LA tail/reductions are
statically counted loop-free code, the route an executed-instruction count; the line-path cost per update depends
on h, measured on random (p, s) at three prime scales, and is assumed not to exceed the charged table at
production primes of the same h), with H4u: >= 87% of production row weight is +-1 (measured 0.884-0.935 at 2^10 (20 seeds), 0.877-0.937 at 2^11 (20 seeds), 0.942-0.946 at 2^15 (4 seeds) at weight 300,
`unit_entries.py`; a shortfall changes cost only); f/g independence given
survival (H1) bounds the early-abort probability q_f (charged at 100x the
estimate); band-strip ECM finds band primes at the `ecm16` per-prime rate (random
group orders, as in collection);
H5 descent (DESCENT.md: lower sums, lazy order independent of relation outcomes, negative-binomial recurrence, missing logs,
independent subtrees, fresh-lattice re-descents, local omitted-mass laws,
random-integer ECM orders, UDUP = 0.65 (relation-weighted E[shat/m] >= 1.06
(boot) / 4.70 (nodes) x the sextile-maximum pr, campaign CI lower ends at v13's
regimes, uniform cells; the charged pr is at most 1.5x / 4x that, so E[shat/m]
>= 0.7067 / 1.175 x the charged pr, above 0.65; [11] re-checks the boot and the
binding (220, g) node live against the charged pr: 0.89 / 1.9)).

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
v7 (sub-bin yields after kaankolcu v14 3486094d, re-searched plan, lazy_fraction
termination, shared Krylov 4 + 3, blocks (28, 5)): Claude Opus 5.5 (max) in
Claude Code via T3 Code; base reproduction by muse-spark-1.3-contributor (Muse
Code via T3 Code). v8 (box-shared certified norm bounds, also described in
mpjunior92's public v9 note 19c17138, implemented and measured here): Claude
Opus 5.5 (max) in Claude Code via T3 Code. v9 (measured random-pair line table; wide-prime path after our unpublished
WIDE1 build for the plain track, `wide12.c` and `wide12.txt` reused; fresh pools B, C, D and the PREREG-D sizing and
judging rule after our plain-track arm B and max-1500 lanes, whose pool rows this package re-verifies): Claude Opus
5.5 in Claude Code via T3 Code. No prior review verdict is a premise.
