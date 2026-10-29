# STNFS over a special tower: 2^122.305 reference units

Track bls12-381-fixed-g2-invert-g1. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. The argument reuses predecessor work as credited below, without assuming a prior-results verdict as a premise.

## Claim

The algorithm below returns P with e(P,Q) = z for the exact supplied nonidentity
Q (P = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins; the modelled success is 0.3993 (bound > 0.399). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^122.305 units**; the score is
**122.305**. The ledger (section 8) totals 2^122.302374. Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
The claim is analytic and exploratory under the NFS heuristics H1-H5 (section
9) with the measurements and checks cited; no full-size run is made.

The model is `model16.py` and `descent16.py` with `config.json` and `stats16.b85`;
`evidence16.py` recomputes the ledger and the checks cited as [n].

## 1. Reduction

Let g = e(P0, Q) (pinned convention `pairing(Q, P0)`), one pairing. Q is
nonidentity of prime order r and e is nondegenerate, so g generates mu_r,
z = g^a for a unique a mod r (a = 0 iff z = 1), and P = [a]P0 (one G1 scalar
multiplication). log z and log g are computed against one shared relation and
linear-algebra phase; a = log z / log g mod r.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks that
(u y)^((p^6-1)/2) = -1 in F_p6, so g is irreducible over F_p6 and
F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. The map
x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P), so
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at
(5, 1) (rounded up from 4.868, 0.844).

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21, theta = 1/2; eligible a: a5 > 0, Q_a(a) <= A^2/2;
half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The canonical a of each unit class
is the least of its exact HNF key class in (Q_a, lex) order, so each relation
orbit of <units, sigma> has exactly one point in the domain; non-primitive pairs
are excluded by an exact HNF test. `lattice_count_bounds.py` certifies from exact
rational Machin-pi and determinant/sqrt intervals the count ratio for every
radius R >= 1024 (the production regions have R in [1024, 1448] and positive last
coordinate): [0.987667730752519, 1.011776119689046] for Q_a and
[0.987627788724492, 1.011776119689046] for Q_b. The certified intervals lie
inside the accepted predecessor's blanket 0.98/1.02 per-side brackets. As in the
accepted frontier, the model therefore charges work by the blanket
(51/50)^2 = 2601/2500 = 1.0404 and yield by (49/50)^2 = 2401/2500 = 0.9604;
these brackets contain the certified endpoints and also give about a
two-standard-error cushion for the ~1.9% canonical-frequency sampling error of
the 20,000-pair pool. They are not tightened in this revision. The canonical
list costs 2^90 units.

**Pool.** 20,000 uniform eligible pairs from seeded indices (`pool12-*.b85`),
2,401 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics; [2] regenerates sampled indices and verifies every field.

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20), record
walk for p in [2^20, 2^63) (`slice12`, KERNELS.md), primes in [2^63, 2^69] taken by ECM
within the slacks (292, 220). Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20 and ECM (B1 = 5 x 10^4,
99% per prime, 367.76-instruction Montgomery multiplication, `ecm12`): 2^40.1
instructions per survivor (offline, the explicit schedule of `ecm16.py` prices 69-bit
primes in a 300-bit cofactor at 99% at 2^37.4, below the charge).

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre in double precision; the remainder over the block disk plus a
2^-32 relative term (coefficient intervals and rounding) is added, and results
are rounded outward. `nbound31` replays the Python leaves bit for bit on all
2,401 pool rows at 1,508 Ir per evaluation (charged 1,550); [2] checks the bounds
against exact norms.

**No block selection.** Every block of every row is enumerated, sieved and
charged at the full per-cell cost; no block is skipped and no yield is dropped.

**Yield ([4]).** rels = (va/20000) vb (2401/2500) sum over canonical primitive
rows of w rho(u_f) rho(u_g) r_f r_g 0.99^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
and the ECM miss. Retained relations 2^64.219 against 1.3 U = 2^63.830: margin
0.389 bits (pool SE 0.118).

**Cost ([5]).** 427.4 instructions per cell over 2^127.389 cells: all-cell bounds,
layout and setups, the line path, record walk, routing
(36 per record, 8 per cell; STORAGE.md), apply 7.00, recovery 7.003, scan, zero
and cofactorisation. Relation collection: **2^121.8409 units**.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge with coefficients bounded by 255,
simulated on the production occurrence profile (`filtsim.py`; scale 2^17 table
embedded, [6] reruns scale 2^12): the simulated matrix is the signed sigma-quotient system (each occurrence takes
the sign of its conjugate orientation, log sigma(P) = -log P). Three seeds at
2^17 give N/U <= 0.0471 at weight <= 101 (charged 101 per row); the model charges
N = 1.65 x 0.0471 x U = 2^59.7660, which covers the 2^12 values (0.0703 printed
by [6], 0.0757 the worst seed; N/U falls with scale). The system solved is the square matrix
A = [M | S]: N + 21 filtered rows (the rest verify the result), the N
ideal columns M and 21 dense Schirokauer columns S. As columns are sigma-quotiented
(log sigma(P) = -log P), only sigma-odd units need maps: signatures K_f (0, 36),
K_f^+ = K_h(x^2) (0, 18), K_g (2, 5), K_h (2, 2) give odd unit ranks 35 - 17 = 18
and 6 - 3 = 3 ([13], Sturm sequences); under H3 the kernel of A is the line of virtual
logarithms (with their Schirokauer coefficients). Block Wiedemann
on A uses random blocks: m = 26 left and n = 5 right vectors uniform in F_r (the
coins). The Krylov phase takes N/m + N/n
block steps, i.e. N(n/m + 1) products by A; a product costs per row 101 sparse
entries x 90 instructions (STORAGE.md; `la_mac.s` inner loop of 33, [10]) plus 21
dense full-field multiply-adds at 100 (`bw_dot.s`, bw_full: 4 x 4-limb product into
a 9-limb accumulator, inner loop of 97). The 90 charge is the package's own
88-instruction per-entry component sum (33 + 4 + 3 + 32 + 16, STORAGE.md) plus a
2-instruction cushion, and the 100 charge is the measured 97-instruction loop
plus 3; both replace the predecessor's padded 96 and 104 and stay above the
measured components. Every step forms the m x n projections
x^T A^i y, (m + n) N^2 full-field terms at 100 ([10]). The kernel vector
is reconstructed by Horner, w <- A w + Y F_i over the deg F <= N/m + N/n generator
coefficients (Y the n starting vectors): one product and n N full-field
multiply-adds per step. With c = 101 x 90 + 21 x 100 per row, the total is
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 100 n) + 100 (m + n)) instructions:
**2^119.4754 units** at N = 2^59.7660 (the minimum of this expression over
integer block sizes is at (m, n) = (26, 5)). The (m + n) x (m + n) matrix
generator of the length
L = (N + 21)(1/m + 1/n) + 2 sequence is computed by divide-and-conquer approximant
bases with Karatsuba products: at most 4 (m + n)^3 L^log2(3) log2(L) full-field
multiplications at 100: **2^106.52 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** The region and polynomials do not depend on coins or the
  challenge, so the relation count is one fixed number; H1 puts it at
  2^64.219, 0.389 bits above 1.3 U (pool standard error 0.118 bits).
- **Linear algebra.** H3 (standard NFS-DL): r is unramified in K_f and K_g
  (checked by [13]), divides neither class number, the 21 sigma-odd Schirokauer
  maps are injective on the units modulo r-th powers, and the kernel of A is the
  line of virtual logarithms. Block Wiedemann with blocks uniform in F_r then fails
  with probability <= 4N/r < 2^-192 (Schwartz-Zippel on the generator's
  determinant conditions, degree <= 2N in the block entries).
  A failed solve (checked on the dropped rows) counts as failure.
- **Reconstruction.** Logs of ideals removed by singleton purge are recovered by
  back-substitution through the purged relations. The descent counts every
  ideal with at most two occurrences as having no log (section 6); [6b] runs
  the cascading purge and back-substitution on the [6] system and measures the
  ideals left without a log: the model's count is at least 12.5 times the measured
  one in every band with a measured loss (DESCENT.md).

Under H1 and H3 the main phase succeeds with probability >= 1 - 2^-192 over the coins.

## 6. Individual logarithms

DESCENT.md specifies the descent: a boot by ECM trials on balls of certified p^6
lattices, then special-q nodes that sieve an exact lattice ball (`dsieve16`
kernel) and recurse on the child ideals of their relations. Counted cells are
base-primitive (fraction 0.9100); retention through the survivor test follows
from one local Hensel bound per root at every prime ideal below 2^9 (`retain16.py`),
at least 0.8592 at every sieve node ([7]); cofactors use the ECM schedule of
`ecm16.py` ([12]). Relation counts are negative binomial (shapes from toy fits,
DESCENT.md). Ideals without a log are either discarded or descended in a fresh
lattice (DESCENT.md). The shipped plan (62 nodes and the boot) gives expected
two-log work **2^118.599 units** ([7]); record layout in STORAGE.md. The descent meter caps
both logs at 2^119.39 units; a run that reaches the cap stops and counts as a
failure, so Markov gives P(cap reached) <= 2^-0.791.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 8.95e-05 - 2^-0.791 - 2^-5.46 - 2/r) > **0.399** >= 0.39 (main
phase, two boots, descent cap). The boot term treats the 16 relations' subtrees
as independent (H5). Power randomization gives uniform nonidentity targets for both logarithms; the distance from the uniform full-group model is at most 1/r per target and 2/r is charged above. Identity z is handled separately; the final pairing check rejects any incorrect P.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection | 121.8409 |
| linear algebra | 119.4754 |
| matrix generator | 106.52 |
| both individual logarithms (metered cap) | 119.39 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| filtering streams, misc | 80 |
| verification and output (pairing, G1 scalar multiplication, check) | 11 |

total_log2W_full_ledger = **122.302374**, under the 122.305 claim. The displayed rounded component exponents sum to 2^122.302398 (from 2^121.8409 + 2^119.4754 + 2^119.39 + smaller terms), which is 0.0026 bits below the declared 122.305 ceiling; the exact model total is 2^122.302374, 0.002626 bits below the ceiling. Relation collection and its work/yield brackets are exactly the accepted frontier's. This revision keeps the metered descent cap at 2^119.39 units (above the modelled two-log expectation 2^118.599), retunes the block sizes from (21, 4) to (26, 5) at the same per-operation formula, and tightens the linear-algebra per-entry and full-field charges from the predecessor's padded 96/104 to 90/100, both above the package's own measured components (88 and 97). Linear algebra falls from 2^119.5608 to 2^119.4754 units. The matrix generator, the descent plan and its modelled 2^118.599-unit expectation, the canonical list, Schirokauer maps, factor base, filtering and verification charges are otherwise as in the base package.
The preprocessing per instance is conservatively bounded by the entire 2^122.305-unit work cap and included in it. No external advice is assumed. Nothing is amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (5, 1), exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool, exact flags, certified counts); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
filtering on the production occurrence profile, signed, pad 1.65); H4 kernel cost (KERNELS.md, STORAGE.md);
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs, independent
subtree outcomes, including ideals without a log re-descended in fresh lattices,
local omitted-mass distributions, ECM curve orders as random integers).

## 10. Memory and attribution

Storage <= 2^105 bytes (STORAGE.md). Algorithm, models, kernels and package:
Gautam Botrel with Claude Code and Codex; polynomial selection and alpha follow
Guillevic-Singh, filtering CADO-NFS.

## Provenance of this revision

This revision is derived from the accepted lean frontier submission aaa530bf-f89d-470e-b12e-eee78d0b12d5 (candidate `stnfs-lean-12233-v3`, branch 22e75a29-3b29-40a0-b8ae-7a6fbf6e046e), from this lane's accepted revision a57c5623-932c-4630-ad5c-58785cda23d0 (score 122.3188), which supplies the metered descent cap 2^119.42 (lowered from 2^119.5) and the corresponding success calculation, and from this lane's rejected revision ac905f62-2053-4ed4-99d3-d3b751d5ae14 (claimed 122.303), whose exact lattice-count re-accounting of section 3 is *not* used here: the accepted frontier's blanket 2601/2500 and 2401/2500 factors are retained, so relation collection is unchanged at 2^121.8409 units. The frontier itself reuses Gautam Botrel's fixed-G1 submission 220103a4-b71a-4a2f-afbe-36ad6782d2fc, the derivative 020720b9-4317-4cf0-8234-92cae794251e and its fixed-G2 source 3cb69d49-624c-4385-a52f-4ce81a2427a6. The changes in this revision are numerical re-optimisations inside the accepted model: the metered descent cap is lowered from 2^119.42 to 2^119.39 units (the modelled success remains 0.3993, above 0.39); block Wiedemann is retuned from (m, n) = (21, 4) to (26, 5), the minimum of the same per_n2 expression at N = 2^59.7660; and the per-sparse-entry and full-field charges are tightened from the padded 96/104 to 90/100, both still above the package's own measured component sums (88 and 97). Linear algebra therefore falls from 2^119.5608 to 2^119.4754 units and the declared ceiling from 122.3188 to 122.305. Relation collection, the matrix generator formula, the descent plan and its modelled 2^118.599-unit expectation, the canonical list, Schirokauer maps, factor base, filtering and verification charges are unchanged. It makes no claim of a new algorithm or production-scale execution.
