# STNFS over a special tower: 2^123.5 reference units

Track bls12-381-fixed-g2-invert-g1. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
advice, smaller targets, future-target amortization or prior results are used.

## Claim

The algorithm below returns P with e(P,Q) = z for the exact supplied nonidentity
Q (P = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins; the modelled success is 0.81. A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^123.5 units**; the score
is **247/2**. The ledger (section 8) totals 2^123.093. Distributed storage is at
most 2^105 bytes (descent score arrays, STORAGE.md) in individually bounded,
explicitly addressed objects. The
claim is analytic and exploratory: the standard NFS heuristics H1-H5 (section 9)
are used with the measurements and checks cited, and no full-size run is made.

Numbers in this package have one source: `model16.py` and `descent16.py` with
`config.json` and `stats16.json`; `evidence16.py` recomputes every figure below.

## 1. Reduction

Let g = e(P0, Q) (pinned convention `pairing(Q, P0)`). Q has prime order r and e
is nondegenerate, so g generates mu_r in F_p12^*, z = g^a for a unique a mod r,
and P = [a]P0. Both logarithms log z and log g are computed against one shared
relation and linear-algebra phase; a = log z / log g mod r. Orbit-12 Pollard rho
would need 2^125.96 group operations.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial composed
with x^2/y) and g = x^2 + |u| y with u = -0xd201000000010000. [0b] checks that
(u y)^((p^6-1)/2) = -1 in F_p6, so g is irreducible over F_p6 and
F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. The map
x -> -x is the Frobenius^6 conjugation sigma; it pairs (a, b) with (a, -b) and
log sigma(P) = -log P, so candidates and unknowns are both quotiented by sigma.
Murphy alpha is charged at (5, 1) (rounded up from 4.868, 0.844).

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21, theta = 1/2; eligible a: a5 > 0, Q_a(a) <= A^2/2;
half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The canonical a of each unit class
is the least of its exact HNF key class in (Q_a, lex) order, so each relation
orbit of <units, sigma> has exactly one point in the domain; non-primitive pairs
are excluded by an exact HNF test. `lattice_count_bounds.py` certifies counts
within [0.98, 1.02] of the volumes; the model scales work by 2601/2500 and yield
by 2401/2500. The canonical list costs 2^90 units.

**Pool.** 20,000 uniform eligible pairs from seeded indices (`pool12-*.b85`),
2,401 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics. [2] regenerates sampled indices from their seeds and
verifies every stored field.

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20), record
walk for p in [2^20, 2^63) (`slice12.c`), primes in [2^63, 2^69] taken by ECM
within the slacks (292, 220). Certified per-block bounds (`embedding_bounds.py`:
rational Rouche root disks, outward-rounded intervals, 64-cell blocks of width
128) give thresholds T = lower - 16 x slack, which imply retention of every
B-smooth cell with omitted mass below slack, and bound survivors' cofactors by
slack + 8 bits. Survivors get exact norms, trial division below 2^20 and ECM
(B1 = 5 x 10^4, 99% per prime, measured 367.76-instruction Montgomery
multiplication in `ecm12.c`): 2^40.1 instructions per survivor.

**Block selection.** Blocks with certified bound log2 N_f + log2 N_g > 1285 are
not sieved; their bound evaluation, layout and setup are still charged. Over all
2,401 rows (113,342 blocks) the cut keeps 0.8078 of cells and 0.9617 of yield.

**Yield ([4]).** rels = (va/20000) vb (2401/2500) sum over canonical primitive
rows of w rho(u_f) rho(u_g) r_f r_g 0.99^(lambda_f + lambda_g), with Dickman rho
(checked against an independent 300-digit integrator to u = 64), retentions r
from the finite Euler-product threshold model (`threshold_model.py`; calibrated
on exact smooth integers at B = 2^10..2^14) and the ECM miss. Retained relations
2^64.162 against 1.3 U = 2^63.830: margin 0.332 bits (pool SE 0.118).

**Cost ([5]).** 541.7 instructions per cell over 2^127.389 cells: all-cell bounds,
layout and setups, plus 0.8078 x the kept-cell line path, record walk, routing
(36 per record, 8 per cell; STORAGE.md), apply 7.00, recovery 7.003, scan, zero
and cofactorisation. Relation collection: **2^122.1826 units**.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge with coefficients bounded by 255,
simulated on the production occurrence profile (`filtsim.py`; scale 2^17 table
embedded, [6] reruns scale 2^12): N/U = 0.0538 at weight 100; the model charges
N = 2 x 0.0538 x U = 2^60.2341 (N/U falls with scale: 0.0696 in the sandbox
rerun at 2^12, 0.0538 at 2^17; the factor-2 pad also covers the 2^12 value).
Block Wiedemann with <= 82 Schirokauer columns
as right-hand sides: 3N products at 96 instructions per entry (STORAGE.md;
`la_mac.s` inner loop of 33, [10]): **2^120.9943 units**. Dropped excess rows
validate the solution.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** The region and polynomials do not depend on coins or the
  challenge, so the relation count is one fixed number; H1 puts it at
  2^64.162, 0.332 bits above 1.3 U. The pool standard error (0.118 bits, 2.8
  per margin) is uncertainty in H1, not a success probability.
- **Linear algebra.** Block Wiedemann over F_r with random projections fails
  with probability about N/r < 2^-190 per run for a nonsingular system; the
  dimension carries a factor-2 pad (x4 on the charged cost) and filtering at
  excess 1.2-1.8 is stable in the simulation ([6]). The solve is verified on the
  dropped rows; a failed verification repeats with fresh projections inside the
  padded budget. H3 (heuristic, as in every NFS-DL analysis) allots <= 0.097
  to rank defects after filtering; the filtering simulation ([6]) is its
  support.
- **Reconstruction.** Logs of ideals removed by singleton purge are recovered by
  back-substitution through the purged relations (as in CADO-NFS
  reconstructlog). The descent counts every ideal with at most two occurrences
  as having no log and never assumes such a log (section 6). [6b] runs the
  cascading purge and back-substitution on the [6] system and measures the
  ideals actually left without a log, including those with >= 3 occurrences.
  The model's count is 9-17 times the measured one in every band
  (DESCENT.md).

Together the main phase succeeds with probability >= 1 - 0.097 = 0.903 under H1 and H3.

## 6. Individual logarithms

DESCENT.md specifies the descent: a boot by ECM trials on balls of certified p^6
lattices, then special-q nodes that sieve an exact lattice ball with the measured
`dsieve16.c` kernel, take relations with a two-tier window of children and
recurse with backtracking over the exact child ideals. Relation counts are
negative binomial with shapes below every toy fit; survivors are bounded at each
bin's lowest norms. Factor-base ideals without a logarithm are handled
explicitly: their expected number per relation follows from the main-phase
occurrence model, and each node either discards relations containing one or
descends them as further children (DESCENT.md). The shipped plan (96 nodes and
the boot) gives expected two-log work **2^117.632 units** ([7]). Node failure
uses the negative-binomial shape of the effective window: the empty-window base
nodes take r = 0.5 and meet phi <= 0.045 on larger balls. The descent meter caps
both logs at 2^121 units, so Markov gives P(cap reached) <= 2^-3.37.

## 7. Success

Success >= 0.903 x (1 - 2 x 2.2e-3 - 2^-3.37) = **0.812** >= 0.39 (main phase,
two boots, descent cap). The boot term treats the 16 relations' subtrees as
independent (H5); success stays >= 0.39 for any per-log boot failure <= 0.23,
100 times the modelled 2.2e-3. Uniform power randomization of the targets makes both
individual-log tasks distributed as the sampled production targets for every Q
and z. The final pairing check rejects any incorrect P.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection | 122.1826 |
| linear algebra | 120.9943 |
| both individual logarithms (metered cap) | 121 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| filtering streams, generator, misc | 80 |
| verification and output | 11 |

total_log2W_full_ledger = **123.092838**, under the 123.5 claim by 0.407 bits.
With every sieve constant x1.5, relation collection is 2^122.7675 and the total
2^123.4332 ([5]), still below the claim. Nothing is amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (5, 1), exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool, exact flags, certified counts); H3 matrix (filtering on the
production occurrence profile, pad 2); H4 kernel cost (KERNELS.md, STORAGE.md);
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs).

## 10. Memory and attribution

Storage <= 2^105 bytes, set by the descent score arrays (R = 384 balls); the main phase uses
<= 2^81 bytes (STORAGE.md). Algorithm, models, kernels and package: Gautam
Botrel with Claude Code and Codex; polynomial selection and alpha follow
Guillevic-Singh, filtering CADO-NFS.
