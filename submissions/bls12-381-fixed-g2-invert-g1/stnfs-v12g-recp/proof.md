# STNFS v12: run-consistent rec7 and an excess-1.0 region judged on a fresh pool, on the v10 core: 2^119.464026123585 reference units

Track bls12-381-fixed-g2-invert-g1. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used.
Predecessor work is reused as credited below; no prior verdict is a premise.

## Claim

The algorithm below returns P in G1 with e(P,Q) = z for the exact supplied
nonidentity Q in G2 (P = infinity for z = 1) with success at least 0.39 over
uniform full-production challenges and fresh coins under the declared heuristics
H1-H5 (model 0.4455539; audit with rounded inputs 0.4455416; 0.3902002 at +10% expected descent work). A global
instruction meter bounds all work, including preprocessing, failures, restarts,
transfers, verification and output, by **2^119.464026123585 units**; the score is
**119.464026123585**. Upper-rounded ledger rows (section 8) sum to
2^119.464026123584594 <= 2^119.464026123585. Storage is at most 2^105 bytes (STORAGE.md). The claim is
analytic and exploratory; no full-size run is made.

**What v12 changes (base: v10p 0fbd40b3, valid 119.611470529661 on this track; its core is recover-Q v10 9ac7a9b4).** (a) **Excess 1.0** (pre-registered, `core-1500-PREREG-E`): theta 2253/8000 sized on pool
B, judged once on a fresh pool E (seeds 71-74, 656,000 draws): margin 0.1221 >= max(0.06, 2 SE); descent inputs
from E (63.531, 738.9/498.9); filtering re-measured at excess 1.0 (section 4); rec7 runs recomputed in the smaller
region. (b) **rec7, one implementation per (run, band)**: for p in [2^52, 2^63) every slice of a clean run (every a-value canonical) of 32 to 1,024
a-values is walked by the record walk per canonical a-run when log2(run cells) minus the band's midpoint log2 p is >= 4.99 (quarter-bit bands), by v10's slice kernel otherwise (v11, never judged because of provider timeouts
and now withdrawn, chose the cheaper walk per row and band, which no single schedule attains; that rule got our plain-track
P9 26d7bdc1 rejected). (c) **Kernels**: scan2 with fused zeroing, apply and recovery at 5.76 per record (`slice12oh.c`, our overhead lane).
(d) **Success headroom**: the meter cap is the smallest with success >= 0.39 at +10% expected descent work (section 7).
(e) **Clean runs** (answers the run-allowance objection of bccaaeba): rec7 walks only runs whose every a-value is canonical;
a contaminated run takes the slice kernel, and every row pays RUN_ALLOW x D_SUP, a bound on that extra work (section 9, H2).
From v10: the plain-track P7/P8 core (filtering at mean + 2 SE, per-entry layout, allowance 1 + 2 SE, band strip at the row norm). Collection
**2^119.088326225**, LA **2^117.238031284**, generator **2^107.761930585**, meter **2^113.425867**.

**Numerical conventions.** Decimals are the experiment's values rounded as
shown; bounds are rounded in the safe direction and `declared_bounds` compares
every stated bound with the unrounded value. The model is `model16.py` and
`descent16.py` with `config.json` and `stats16.b85`; `evidence16.py` recomputes
the ledger and the checks cited as [n].

## 1. Reduction

Let g = e(P0, Q) (pinned `pairing(Q, P0)`), one pairing; P0 is the fixed G1
generator. Q has prime order r and e is nondegenerate, so g generates mu_r,
z = g^a for unique a mod r (a = 0 iff z = 1), and P = [a]P0 (one G1 scalar
mult), so e(P, Q) = z for that exact Q (never replaced). log z and log g share
one relation/LA phase; a = log z / log g mod r. g and z enter the tower by
the field isomorphism (roots of h and of X^2 + |u| y in F_p12, a 12 x 12 basis
change), charged with the pairings, [a]P0 and the check in the 2^21 row
(section 8: < 2^20 F_p12 multiplications in all).

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
Q_a(a) <= (5/16) A^2); the collected region is its sub-region theta = 2253/8000
(`production_theta_region`): a5 > 0, Q_a(a) <= (2253/8000) A^2; half-slice
b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The restriction of a uniform pool to qa <= 2253/8000 is uniform on the sub-region;
each pool is scaled by the sub-region volume over its effective sample 656000 x (2253/2500)^3 = 480139.448. The
canonical a of each unit class is the least of its exact HNF key class in (Q_a, lex) order, so each relation orbit of
<units, sigma> has one point in the domain; non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` certifies (exact rational Machin-pi, determinant/sqrt intervals) the count ratio at every
radius R >= 768 = isqrt((2253/8000) 2^21): [0.983579, 1.015728] for Q_a and [0.983525, 1.015728] for Q_b inside [0.98, 1.02] (asserted); the
model uses the exact products **0.967375771 (yield, floor)** and **1.031701500 (work, ceiling)**.
The canonical list costs 2^90 units.

**Pools (v9; `pool_d.b85` (base85 of LZMA JSON)).** Every pool draws 164,000 uniform eligible pairs per seed on theta = 5/16 from seeded
indices of the package sampler; every draw is tested with the exact canonical predicate and every canonical row
kept (no selection). B = seeds 7-10, C = 11-14, D = 41-44, E = 71-74 (656,000 draws each). Pre-registered before any statistic of
E (PREREG-E): theta is the smallest 1/8000 step with margin_B >= 0.06 + sqrt(2) SE_B at excess 1.0 (2253: margin
0.087893 >= 0.087257; 2252: 0.086774 < 0.087277, [3s]); E judges it once
(margin_E >= max(0.06, 2 SE_E)); the work is estimated on B, C and D (D judged the excess-1.05 region of v9-v11). `pool_d.b85`
ships each seed's region sums and index-list hash; [2] regenerates E seed 71 draws [0, 1000) and must reproduce that
block's sums. E: 78,665 canonical rows, 57,491 in the region.

**Sieve, thresholds, cofactorisation.** Line sieve (scores only) for p in [2^8, 2^48), record
walk for p in [2^48, 2^64) (KERNELS.md: public `slice12` line walk; record walk `slice12g` below
2^52, per canonical a-run (`rec7.c`) in [2^52, 2^63), `slice12` + `wide12.c` in [2^63, 2^64); X2 = 2^64), primes in [2^64, 2^69] taken by ECM
within the slacks (295, 226). **Line sieve to 2^48 with band strip (v3; 8cc93834 /
42e80f80 lineage):** primes in (2^20, 2^48] only add scores on the line path
(every line keeps >= 2^7.11 hits per coset at p = 2^48, inside the measured line
table) and emit no records: records/cell 2.33 -> 0.58. Their factors are recovered at the
survivors by a band strip: two rounds per side of the `ecm16.py` 48-bit schedule at
99.9% on the unstripped cofactor; it divides the norm, so v10 charges each row at its own norm size,
min(cap, ceil(log2 norm) + 1) bits (in-model: the `ecm16` schedule's count at that size, no separate measurement; caps 832 / 576 bits; 2 x 2^36.4488 f side and 2 x 2^35.3903 g side at the caps). A band prime missed in both rounds loses the relation: yield x
(1 - 41 x 0.001^2)(1 - 28 x 0.001^2) (at most floor(832/20) = 41 and 28 band
factors per side). Survivors' cofactors are bounded by slack + 8 bits: 303 (f, slack 295) and 234 (g, slack 226); slacks are
v10's (not re-optimised at theta 2253/8000; yield and cost both use them). The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at **99.9%** per prime (ECM_FIND = 0.999), which [12] checks: the cheapest
grid point for 69-bit primes in a 303-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves,
five-limb Montgomery multiplication at 367.76 Ir, `ecm12`) prices one side at 2^38.0906
instructions and the other (234 bits) at 2^38.0879. **Early abort (2c2ce2ab):** the g-side band strip and ECM run
only after the f side factors completely (yield unchanged). Per survivor: C_NORMS + C_TDIV + band_f + 2^38.0906 +
q_f (band_g + 2^38.0879), band_s the two rounds at the row's norm size, q_f = min(1, 100 rho(u_f)/P(omitted_f <= slack_f))
per pool row (heuristic under H1, mean ~1e-5, factor 100 a margin). The 99.9% is the same
ECM_FIND used in the yield's miss model.
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.

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
(C_NB_ROW 4,600 = the measured row context, `row_setup` 4,292.1 rounded to 4,300, + 300 for the hull update: about 25
instructions per place and row over at most 4 places, <= 100 counted, charged 300).
[2b] rebuilds six boxes of replayed E rows and checks 1,536 exact norms at sampled cells of every box row
against their leaf bounds; [2] replays every field of the replayed rows. Every block of every row is enumerated, sieved and
charged (no selection).

**Yield ([4], judged pool E).** rels = (va_region/480139.448) vb (0.967375771) sum over canonical primitive sub-region
rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
the ECM miss and the band-strip loss. **Central core (v4):** relations from
cells with qa = Q_a(a)/A^2 < Q_CORE = 0.02 are discarded (never passed to
filtering, linear algebra or the descent; those cells are still sieved and
charged), so count and relation-weighted mean norms are those the estimator
measures over the rows with qa >= 0.02; the discard only lowers the count. Retained relations >= 2^63.5734
against U = 2^63.4514 (excess 1.0, section 4): margin >= 0.1221 bits (pool SE 0.0212,
so 5.7 SE). Feasibility requires a margin of at least max(`MARGIN` = 0.06, `MARGIN_SE` = 2 SE). The
descent's missing-log input `LOG2_RELS = 63.531` is the judged estimate minus 2 SE, floored to 0.001 bit (unrounded
2^63.53107; red line (a)); [4] checks it. The mean-norm inputs 738.9/498.9 are E's relation-weighted means 738.94/498.91 rounded down, checked in [4].

**Cost ([5], work pools B, C, D).** Sampled mean 187.7 instructions/cell (v10 207.8); charged mean
188.6 with the allowance below, over 2^125.817 cells: all-cell bounds, slice layout (0.109 per cell:
28.17 per suffix-box entry, 2^21 per slice, 2^17 per segment; KERNELS.md), setups, line path at the v9 table, record walk
(`slice12g`; rec7 per (run, band) in [2^52, 2^63) on clean runs of 32-1,024 a-values, + 0.011 x D_SUP per row; wide primes at `REC_TAB_W` x 1.01), routing (10.01/record, 8/cell; STORAGE.md), apply 5.76, recovery
5.76/record and scan2 with fused zeroing (`slice12oh.c`), band strip at the row's norm size and f-first cofactorisation. Collection:
**<= 2^119.088326225 units** (estimate 2^119.0811 times allowance 1.005).

**Canonical-work allowance (H2, [4a]).** The charged work is (va_region/n_eff) times the per-row sum over in-region
canonical rows of B, C and D: a Monte Carlo estimate of a fixed, challenge-independent domain total over 1,968,000
draws (172,740 contributing); its relative SE has a conditional row term 0.0933% and an incidence term
sqrt((1 - 172740/1968000)/172740) = 0.2299%, total <= 0.2481%. The estimate is multiplied by
**1.005 >= 1 + 2 SE = 1.004961** (rounded up; [4a] fails below 1 + 2 SE). v9 kept v8's 1.0174 as a floor;
with three 656k pools the SE is measured directly, so v10 (as the plain-track P7) drops it. A modeling margin, not a probability.

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U, `filtsim.PEXC`), clique removal and Markowitz merge (coefficients <= 255, including a
collision-safe bound on implicit middle-range entries), simulated on the production occurrence profile (`filtsim.py`;
[6] reruns 2^12) as the signed sigma-quotient system. **Required excess 1.0.** The whole region is sieved and every
retained relation (>= 2^63.5734) is filtered and used; the matrix is charged at excess 1.0, which only
over-states N. A campaign at pool B's means 739.9/499.0 (fixed before E existed, PREREG-E; charge rule of the plain track's, PREREG-F; every run shipped
in `FILT_RUNS`, 10, 6, 6, 6, 12, 16, 8 seeds at 2^12..2^18) measures x = N(300)/U, the merged weight-300 matrix's rows over U:
mean (SE, seeds) 0.07402 (0.00265, 10), 0.05101 (0.00273, 6), 0.04770 (0.00113, 6), 0.04508 (0.00134, 6), 0.04523 (0.00060, 12), 0.04241 (0.00035, 16), 0.04110 (0.00043, 8) at 2^12..2^18. The means fall with scale overall. At production the realised x concentrates
on its expectation, which H3 assumes does not rise above the largest simulated scales; the charge is the larger of
mean + 2 SE at 2^17 and 2^18, 0.043098, rounded up to **0.0431** (2^17: 0.0424056 + 2 x 0.0003459 = 0.0430974; 2^18: 0.0410953 + 2 x 0.0004260 = 0.0419473) ([10]; the means at 2^17 and 2^18 lie below
the 2^15 mean; the 2^16 mean sits above the 2^15 one, and the falls from 2^16 to 2^17 and from 2^17 to 2^18 are
4.0 and 2.3 SE). v9's maximum over seeds rises with the number of seeds and is no longer used.
The merge stops at mean weight 300 (charged 301 per row); the largest merged coefficient bound is 255 <= 255.
N = 0.0431 U <= **2^58.9152** < 8 x 2^56.

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
  13.204 -> **13.21** (>= 87% +-1 by H4u; measured >= 0.8867 at excess 1.0, STORAGE.md).
  Per row and vector: 301 x 13.21 + C_ROW + 21 x 97.
- **Row tail.** Eight active 2^56-element owners: seven sends and seven
  combinations, `_finish`, headers, routing: counted 1,424 by [10]
  (`owner10.row_tail`), charged C_ROW = 1,430.
- **Horner and projections.** Horner (w <- A w + Y F_i, deg F <= N/m + N/n) uses
  single-vector products at 301 x 20.9 + 1,430 + 21 x 97 plus n dense terms and
  one reduction (420) per element. Projections: (m + n) N^2 terms of a 128-bit X
  entry times a field element at 60 (`bw_dot2`, counted 59) plus a reduction
  per entry.

With c_K and c the Krylov and single-vector row costs: N^2 (c_K (n/m + 1) +
(1/m + 1/n)(c + 97 n + 420) + 60 (m + n)) + 420 (m + n) N instructions =
**<= 2^117.238031284 units** at N = 2^58.9152. The
(m + n) x (m + n) generator of the length L = (N + 21)(1/m + 1/n) + 2 sequence:
<= 4 (m + n)^3 L^log2(3) log2(L) field mults at 97 + 420: **<= 2^107.761930585 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** Region and polynomials are coin/challenge-independent, so the
  count is one fixed number; H1 puts it at >= 2^63.5734, >= 0.1221 bits above
  1.0 U (judged pool E, SE 0.0212); all of them are filtered and the descent uses 63.531, 2 SE below it.
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
  is >= 5.0x measured in every band with loss at excess 1.0 (DESCENT.md).

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
ends a node, so per-relation work is done on demand. (i) After the full sieve a node cofactors its survivors in random order, descending each
distinct relation at once (bound `lazy_fraction`, DESCENT.md; [7l] simulates the
process at the three costliest nodes: 0.0222/0.0175/0.0382 vs bounds
0.0250/0.0178/0.0456, within the 0.01 Monte Carlo allowance; 40,000 runs: 0.0158, 0.0399). (ii) The boot finds and descends one relation at a time,
at most 18: by Wald (1 - psi^18)/(1 - psi) = 2.335 expected relations, 2^56.80
trials. Failure events are those of the recurrence (DESCENT.md).

Expected two-log work **<= 2^112.5723 units** ([7]; one log 2^125.860
instructions: boot own about 2^123.81 + subtree about 2^124.24 x 2.335). The
meter caps both logs at **2^113.425867 units** (cap reached = failure): P(cap) <=
E[W]/cap <= 2^(112.5723 - 113.425867) <= 0.5534148 (Markov; unrounded 0.5534025).

## 7. Success

Success >= (1 - 2^-66) x (1 - 2 x 4.249e-05 - 0.5534148 - 2^-11.45 - 2^-10.70 - 2/r)
>= **0.4455** >= 0.39 (audit 0.4455416; model 0.4455539): main phase, two boots
(phi <= 4.249e-05 each), the instruction-meter Markov bound, the node-counter
Markov bound and special-q repeats. **Headroom:** the cap is the smallest for which the same bound with the
Markov term x 1.1 (+10% expected descent work) stays >= 0.3902: audit 0.3902002 (`audit_refinement.py` success_headroom),
model 0.3902136 (evidence16 success_head).

- Instruction meter: cap 2^113.425867 units, P(cap) <= 0.5534148 (section 6).
- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX
  = 2^27.3. E[N] uses the same recurrence with unit node cost (`counts16.py`):
  E[N] <= 2^15.85 for both logs, P(abort) <= 2^-11.45 <= 0.0003575.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for
  1/p weighting: P <= 2^(2 x 27.3 - 1 + 6 - 70.30) <= 2^-10.70, charged.

No lower bound on node work is used. Subtrees are treated as independent (H5).
Targets are uniform nonidentity by power randomisation (2/r charged); z = 1 is
handled directly; the final check e(P, Q) = z rejects any incorrect P.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (excess 1.0, band strip at row norm, box bounds, v9 line table, per-entry layout, rec7 run walk, wide-prime path; pools B, C, D; estimate 119.0811 + 0.50% allowance) | <= 119.088326225 |
| linear algebra (weight 301, N(300)/U 0.0431, N = 2^58.9152, eight owners, five-vector Krylov passes, blocks (28, 5)) | <= 117.238031284 |
| matrix generator | <= 107.761930585 |
| both individual logarithms (metered cap; expected 112.5723) | 113.425867 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| field map, verification and output (pairing, G1 scalar multiplication, check) | 21 |
| wide-prime factor-base lattices (`wide12.c`, KERNELS.md) | 70.6 |
| rec7 bases (two reduced 7-dim bases per ideal, p in [2^52, 2^63); KERNELS.md) | 67.3 |

**Fixed rows (upper bounds; units x 2^14.29 = instructions).** Canonical list
2^90: 2^41.8 instructions per eligible a (2^62.5 of them) against one LLL-reduced
HNF key, the `canonical_a` orbit enumeration and a segmented key sort (< 2^30).
Schirokauer maps 2^88: 2^38.5 per relation (<= 2^63.8) against one exponentiation
per side mod r^2 (< 2^34.5). Factor base 2^72: 2^23.4 per prime p < 2^69 against
root finding (< 2^19). Filtering, purge, merge bookkeeping and back-substitution
2^80: 2^24.4 per relation entry (2^69.9) against < 2^14. Certificates 2^50;
field map, pairings, [a]P0, check 2^21 (< 2^20 F_p12 multiplications).

The upper-rounded rows sum to 2^119.464026123584594, so W <=
**2^119.464026123585** (audit-checked).

**Configuration.** Excess 1.0; pools B (sizing), E (judged), B, C, D (work), 164,000 draws per seed on theta
5/16; region theta 2253/8000, Q_CORE 0.02, radius 768, allowance 1.005, layout 28.17 per entry, X2 = 2^64, slacks (295, 226), alpha (5.15, 0.85);
line table (18.81, 17.2), (13.82, 24.1), (9.82, 38.99), (5.82, 88.32), (1.83, 352.1); wide primes `REC_TAB_W` x 1.01; rec7 for p in [2^52, 2^63) on runs of 32-1,024 a-values (`REC7_TAB` steps, +2 per record, clean runs only, 0.011 x D_SUP per row); descent inputs 63.531, 738.9/498.9; sub-bin yields (caps
4x/1.5x), plan and boot of v7, UDUP 0.65, THIN_MAX 85, BOOT_RELS 18, lazy nodes and boot; meter 2^113.425867, node counter
2^27.3. LA: purge 1e-4, weight 301, N(300)/U 0.0431, C_ROW 1,430 (eight owners), 13.21 (Krylov) / 20.9 (Horner) per
entry, 97 per dense term, 420 per reduction, blocks (28, 5) with 128-bit X entries.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) in collection and descent against the computed values of section 2 ([14]); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pools and sub-region sized on B and judged once on a fresh pool E, yield lower bound without the qa < 0.02 core,
exact flags, certified counts, and the canonical work mean on B, C, D charged with the 0.50% allowance; rec7 runs end at the first listed unit
associate (11,325 units of log-radius <= 9) or the region edge. Production walks a run with rec7 only if every a-value of it
is canonical (canonical list); otherwise the slice kernel takes every slice of it, so no non-canonical slice is ever walked.
A computed run can end before the exact canonical interval (runs16 also stops at listed associates with a5 = 0, which the
exact test skips, and at the floored region bound), down to a single a-value; the a-values cut off form other runs, and
runs below 32 a-values take the slice kernel (`RUN_T`). The model charges every run in the window at the rec7 price and
every row RUN_ALLOW x D_SUP (KERNELS.md). A row of a contaminated run costs at most D_SUP beyond its rec7 charge (a
deterministic bound over all slices, whatever the contamination), and the share of canonical rows in contaminated runs is
at most RUN_ALLOW = 0.011 >= 1 - 0.05^(1/272) = 0.0110 at 95% confidence: the exact canonical test of all 97,652 a-values of the
runs of 272 randomly drawn in-region pool rows (plain P10b's 300-row exact sample, restricted to this region) found no contaminated run (a modeling margin, not a probability term). The
same sample checks the truncated unit list: no unit beyond log-radius 9 changed a canonical decision on it; under the
clean-run rule the unit list only changes how often a run is contaminated); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
signed filtering on the production occurrence profile at relation excess 1.0 (purge target 0.0001),
the expected N(300)/U not rising from 2^17 and 2^18 (charged mean + 2 SE there, 0.0431) to production); H4 kernel cost (KERNELS.md, STORAGE.md; LA tail/reductions are
statically counted loop-free code, the route an executed-instruction count; the line-path cost per update depends
on h, measured on random (p, s) at three prime scales, and is assumed not to exceed the charged table at
production primes of the same h; the rec7 cost per hit is assumed not to exceed `REC7_TAB`'s per-level
worst case at production runs (32-1,024 a-values: 32-128 and 181-1,024 measured, 129-180 in-model) and primes below 2^63; rec7 work at h' above the
campaign's 2^22 (up to about 2^25.4) and on slices above 2^66 (up to 2^67.5) is charged at the top level 63.72 measured
at 2^22, assuming the cost per hit keeps falling with h'), with H4u: >= 87% of production row weight is +-1 (measured >= 0.8867 at excess 1.0,
weight 300, 2^10..2^11, `unit_entries.py`; a shortfall changes cost only); f/g independence given
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
Claude Opus 5.5) ending in 9ff4fc85 (120.3702). Claude Opus 5.5 in Claude Code (via T3 Code from v7) built
v6 (lazy descent, five-vector Krylov passes, weight-300 merge, eight owners); v7 (sub-bin yields after kaankolcu v14
3486094d, re-searched plan, lazy_fraction termination, shared Krylov 4 + 3, blocks (28, 5); base reproduction by
muse-spark-1.3-contributor, Muse Code via T3 Code); v8 (box-shared certified norm bounds, also described in mpjunior92's
public v9 note 19c17138); v9 (random-pair line table; the wide-prime path and `wide12.c`/`wide12.txt` of our unpublished
plain-track WIDE1 build; pools B, C, D and PREREG-D of our plain-track arm B and max-1500 lanes, re-verified here); v10
(the plain-track max-1500 lane's P7/P8 core: filtering campaign and mean + 2 SE charge, allowance rule, band strip at
the row norm; layout measured by our overhead lane); v11 (`rec7.c`, its campaign and `REC7_TAB` by our BIGLEVER-1500
research lane, also Claude Opus 5.5); v12 (scan2, apply and recover rows: our overhead lane, Muse Code with
muse-wonder-full-1.0, measured in its v10oh 0fb4319a; per-leaf-length scan points by our plain-track max-1500 lane). No prior review verdict is a premise.
