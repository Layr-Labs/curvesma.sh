# STNFS over a special tower: 2^122.192 reference units

Track bls12-381-fixed-g2-invert-g1. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. The argument reuses predecessor work as credited below, without assuming a prior-results verdict as a premise.

## Claim

The algorithm below returns P with e(P,Q) = z for the exact supplied nonidentity
Q (P = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins; the modelled success is 0.4200 and the
directed-rounding bound of section 7 is 0.4199 > 0.39. A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^122.192 units**; the score is
**122.192**. The ledger (section 8) totals 2^122.191732. Distributed storage is at most
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
within the slacks (295, 225). The f-side slack is 295, so survivors'
f-side cofactors are bounded by 303 bits; the g-side slack is 225, so survivors'
g-side cofactors are bounded by 233 bits. The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at 99% per prime, the same schedule the
descent uses and [12] checks: the cheapest grid point for 69-bit primes in a
303-bit (resp. 233-bit) cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 271 curves,
five-limb Montgomery multiplication at 367.76 Ir, `ecm12`) prices one side at 2^37.5177
instructions and the other at 2^37.5150, so the per-survivor charge
C_NORMS + C_TDIV + 2^37.5177 + 2^37.5150 is 2^38.5165 (the predecessor's crude
ECM table charged 2^40.06 per survivor, 2^1.545 above the package's own
schedule; the 99% figure is the same ECM_FIND used in the yield's miss model).
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20 and ECM (B1 = 2^16 ~ 5 x 10^4,
99% per prime, 367.76-instruction Montgomery multiplication, `ecm12`): 2^38.5165
instructions per survivor for the two sides together (C_NORMS + C_TDIV + the
schedule's own per-side counts at the same 99%).

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre in double precision; the remainder over the block disk plus a
2^-32 relative term (coefficient intervals and rounding) is added, and results
are rounded outward. `nbound31` replays the Python leaves bit for bit on all
2,401 pool rows at 1,508.2 Ir per evaluation (charged 1,508.2, the measured marginal
per-block cost excluding the separately charged per-row and per-slice setups); [2] checks the bounds
against exact norms.

**No block selection.** Every block of every row is enumerated, sieved and
charged at the full per-cell cost; no block is skipped and no yield is dropped.

**Yield ([4]).** rels = (va/20000) vb (2401/2500) sum over canonical primitive
rows of w rho(u_f) rho(u_g) r_f r_g 0.99^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
and the ECM miss. Retained relations 2^64.252 against 1.3 U = 2^63.830: margin
0.422 bits (pool SE 0.118), above the model's 0.30-bit feasibility floor
(`MARGIN` in model16.py). The descent's missing-log input is not the pool mean
but a supported lower value: `LOG2_RELS = 64.134 = floor(2^64.2517 - 2^0.1175)`,
the pool estimate minus one standard error, in `descent16.py`; the 2401/2500
yield bracket is already applied inside the estimate, so the input is 0.176 bits
below the raw pool mean. [4] checks `LOG2_RELS <= log2_rels - yield_se_log2`
explicitly (the `rels_input` check). The descent therefore prices every
missing-log work and failure figure at a relation count 1 SE below the
2,401-canonical-row pool mean, answering the review of 2e582f95 that flooring
the mean itself is not a supported lower bound.

**Cost ([5]).** 407.0 instructions per cell over 2^127.389 cells: all-cell bounds,
layout and setups, the line path, record walk, routing
(36 per record, 8 per cell; STORAGE.md), apply 7.00, recovery 7.003, scan, zero
and cofactorisation. Relation collection: **2^121.7703 units**.
The record walk is charged at the package's own measured argument: section B of
`measurements12.txt` states `h = 2 cells / p` for the record path and tabulates
the per-hit charge against `log2 h`; the production charge is therefore
`log2h = log2(2S) - lp`, one bit above the predecessor's `log2 S - lp`, which
overcharged every record-walk point below the table top. [5b] checks the five
measured anchors.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge with coefficients bounded by 255,
simulated on the production occurrence profile (`filtsim.py`; scale 2^17 table
embedded, [6] reruns scale 2^12): the simulated matrix is the signed sigma-quotient system (each occurrence takes
the sign of its conjugate orientation, log sigma(P) = -log P). Three seeds at
2^17 give N/U <= 0.0471 at weight <= 101 (charged 101 per row); the model charges
N = 1.61 x 0.0471 x U = 2^59.7306; the charged ratio 1.61 x 0.0471 = 0.07585
still covers the 2^12 values (0.0703 printed by [6], 0.0757 the worst seed;
N/U falls with scale). The system solved is the square matrix
A = [M | S]: N + 21 filtered rows (the rest verify the result), the N
ideal columns M and 21 dense Schirokauer columns S. As columns are sigma-quotiented
(log sigma(P) = -log P), only sigma-odd units need maps: signatures K_f (0, 36),
K_f^+ = K_h(x^2) (0, 18), K_g (2, 5), K_h (2, 2) give odd unit ranks 35 - 17 = 18
and 6 - 3 = 3 ([13], Sturm sequences); under H3 the kernel of A is the line of virtual
logarithms (with their Schirokauer coefficients). Block Wiedemann
on A uses random blocks: m = 26 left and n = 5 right vectors uniform in F_r (the
coins). The Krylov phase takes N/m + N/n
block steps, i.e. N(n/m + 1) products by A; a product costs per row 101 sparse
entries x 88 instructions (STORAGE.md; `la_mac.s` inner loop of 33, [10]) plus 21
dense full-field multiply-adds at 97 (`bw_dot.s`, bw_full: 4 x 4-limb product into
a 9-limb accumulator, inner loop of 97). The 88 charge is exactly the package's
88-instruction per-entry component sum (33 + 4 + 3 + 32 + 16, STORAGE.md) and the
97 charge is exactly the measured bw_full loop; both replace the predecessor's
padded 96 and 104 and this revision's cushioned 90 and 100, and neither exceeds
its measured component. Every step forms the m x n projections
x^T A^i y, (m + n) N^2 full-field terms at 97 ([10]). The kernel vector
is reconstructed by Horner, w <- A w + Y F_i over the deg F <= N/m + N/n generator
coefficients (Y the n starting vectors): one product and n N full-field
multiply-adds per step. With c = 101 x 88 + 21 x 97 per row, the total is
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n) + 97 (m + n)) instructions:
**2^119.3684 units** at N = 2^59.7306 (the minimum of this expression over
integer block sizes is at (m, n) = (26, 5)). The (m + n) x (m + n) matrix
generator of the length
L = (N + 21)(1/m + 1/n) + 2 sequence is computed by divide-and-conquer approximant
bases with Karatsuba products: at most 4 (m + n)^3 L^log2(3) log2(L) full-field
multiplications at 97: **2^106.42 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** The region and polynomials do not depend on coins or the
  challenge, so the relation count is one fixed number; H1 puts it at
  2^64.252, 0.422 bits above 1.3 U (pool standard error 0.118 bits); the descent
  uses 64.134, the estimate minus one standard error, so its figures are not
  taken at the pool mean.
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
lattice (DESCENT.md). The shipped plan (62 nodes and the boot, 9 kept boot
relations) gives expected two-log work **2^118.171 units** ([7]), priced at the
relation input `LOG2_RELS = 64.134` of section 3; record layout in STORAGE.md. The descent meter caps
both logs at 2^119.034 units; a run that reaches the cap stops and counts as a
failure, so Markov gives P(cap reached) <= 2^-0.862 (the experiment's directed
rounding of the exact model value 2^-0.8625 = 0.5499893).

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 8.101e-03 - 2^-0.862 - 2^-6.18 - 2/r) > **0.4199** >= 0.39 (main
phase, two boots, descent cap); the exact model success is 0.4200. The boot term
treats the 9 relations' subtrees as independent (H5). The 0.030 margin over the
0.39 floor absorbs a 2x increase in the charged boot failure term (success
0.4038), a 5.4% increase in the modelled two-log expectation, a 2x error in the
markov exponent's complement, or any single small term being wrong by a factor
of two. Power randomization gives uniform nonidentity targets for both logarithms; the distance from the uniform full-group model is at most 1/r per target and 2/r is charged above. Identity z is handled separately; the final pairing check rejects any incorrect P.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection | 121.7703 |
| linear algebra | 119.3684 |
| matrix generator | 106.42 |
| both individual logarithms (metered cap) | 119.034 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| filtering streams, misc | 80 |
| verification and output (pairing, G1 scalar multiplication, check) | 11 |

total_log2W_full_ledger = **122.191732**, under the 122.192 claim. The displayed
rounded component exponents sum to 2^122.191788 (from 2^121.7703 + 2^119.3684 +
2^119.034 + smaller terms), 0.000212 bits below the declared 122.192 ceiling;
the exact model total is 2^122.191732, 0.000268 bits below the ceiling. The work
and yield brackets (2601/2500 and 2401/2500) are exactly the accepted frontier's.
Relative to the accepted 122.2065 frontier (cbab3319), this revision changes:
(i) the ECM cofactor slacks are retuned from (291, 222) to (295, 225); the
survivors' bounds are 303/233 bits and the retention at the mean norms is
0.7917/0.9087; the per-survivor ECM charge stays the explicit
`ecm16.py` schedule at 99% (the same grid point as the base): 2^38.5165 per
survivor including the norm and trial-division charges;
(ii) the retained relation count is 2^64.252, 0.422 bits above 1.3 U; the
descent missing-log input is the pool estimate minus one standard error,
`LOG2_RELS = 64.134 = floor(2^64.2517 - 2^0.1175)`, and [4] checks
`LOG2_RELS <= log2_rels - yield_se_log2` exactly. This answers the binding
review of 2e582f95, which held that flooring the pool mean itself is not a
supported lower bound: every missing-log figure is now priced at a relation
count 0.1175 bits below the 2,401-row pool mean (0.176 bits below the raw mean
before the 2401/2500 yield bracket);
(iii) the boot keeps 9 relations instead of 16, so its failure probability is
psi^9 = 8.101e-03 with psi = 0.5856 (charged twice, once per boot) instead of
psi^16 = 9.700e-05; the 62-node plan is unchanged (its minimum per-node relation
count is 35.4, far above either value). With the relation input of (ii) the
shipped plan and boot give a two-log expectation of 2^118.171 units (one log
2^131.459; the k=16 boot at the same input gives 2^118.27), and the metered cap
is 2^119.034 units (directed-rounded Markov abort 2^-0.862, exact
2^-0.8625); the modelled success is 0.4200, a 0.030 margin over the 0.39 floor;
(iv) relation collection is 2^121.7703 units at 407.0 instructions per cell;
the ledger is 2^122.191732 and the declared ceiling is 122.192. The other
components are the accepted 122.2065 tree's: the record path is charged at its
measured argument `h = 2 cells / p` (`measurements12.txt` section B),
`log2h = log2(2S) - lp`; the per-survivor ECM is the explicit `ecm16.py`
schedule (2^38.5165); the norm-bound charge is the measured marginal 1,508.2;
linear algebra uses block sizes (26, 5), per-entry charges 88/97 and pad 1.61.
The norm-bound charge is the measured marginal 1,508.2 Ir
((bench3-bench1)-(setup3-setup1))/(2 x 6509) = 1,508.19; the LA per-entry
charges are the package's 88-instruction component sum and the measured
97-instruction bw_full loop; block sizes (26, 5) remain the minimum of the same
per_n2 expression, and linear algebra is 2^119.3684 units. The pad 1.61, the
62-node plan, the generator formula and the remaining charges are as in the base
package. Preprocessing is bounded by the entire 2^122.192-unit cap and included
in it; no external advice, nothing amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (5, 1), exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool, exact flags, certified counts); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
filtering on the production occurrence profile, signed, pad 1.61); H4 kernel cost (KERNELS.md, STORAGE.md);
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs, independent
subtree outcomes, including ideals without a log re-descended in fresh lattices,
local omitted-mass distributions, ECM curve orders as random integers). H5 is a
disclosed heuristic, not derived: section 7 charges its two consequences
explicitly (the psi^9 boot failure and the cap) and keeps a 0.030 success margin
over the 0.39 floor for it.

## 10. Memory and attribution

Storage <= 2^105 bytes (STORAGE.md). Algorithm, models, kernels and package:
Gautam Botrel with Claude Code and Codex; polynomial selection and alpha follow
Guillevic-Singh, filtering CADO-NFS.

## Provenance of this revision

This revision is derived from this lane's accepted 122.2065 tree
`stnfs-lean-122206-v1` (submitted as cbab3319-552d-4143-955c-e7fc8f3056de,
official 122.2065, ledger 2^122.205983), which repriced the per-survivor ECM
cofactor charge to the package's explicit `ecm16.py` schedule at 99% (2^38.5163
instead of the crude 2^40.06) and retuned the cofactor slacks to (291, 222).
The intermediate revision submitted as 2e582f95-0626-40a5-9a05-d42d1f36c3fd
(claimed 122.172) was rejected: its review held that flooring the pool mean for
LOG2_RELS is not a supported lower bound, and that a modelled success of 0.3912
leaves no margin for the H5 descent extrapolation. This revision answers both
points: the descent input is the pool estimate minus one standard error
(0.1175 bits), and the success target is 0.4200, a 0.030 margin over the 0.39
floor; the boot keeps 9 rather than 16 relations and the cofactor slacks are
retuned to (295, 225). Relation collection is 2^121.7703 units and the ledger is
2^122.191732.
The lineage of the 122.2065 tree runs through this lane's accepted revisions
0ac63018-97d5-4b08-a774-f9f666926f02 (122.224, which fixed the relation-input
rejection of 97b80666 by re-deriving LOG2_RELS from the model), fc9fd9a0
(122.238, corrected Markov text), a73abc31 (122.305, block sizes (26, 5)) and
a57c5623 (122.3188, metered descent-cap method), with the lean base
aaa530bf-f89d-470e-b12e-eee78d0b12d5 (122.33, branch
22e75a29-3b29-40a0-b8ae-7a6fbf6e046e). They supply the measured marginal
charges 1,508.2 and 88/97, the pad 1.61, the boot (215, 43, 15, 4), the
62-node plan, the work/yield brackets and the metered-cap method.
The rejected lattice-count re-accounting ac905f62 is not used: the accepted
frontier's blanket 2601/2500 and 2401/2500 factors are retained and only the
cofactor slacks, the descent relation input and the measured charges are
retuned. The frontier itself reuses Gautam Botrel's fixed-G1 submission
220103a4-b71a-4a2f-afbe-36ad6782d2fc, the derivative
020720b9-4317-4cf0-8234-92cae794251e and its fixed-G2 source
3cb69d49-624c-4385-a52f-4ce81a2427a6. The declared ceiling falls from the
accepted 122.2065 frontier to 122.192. The revision makes no claim of a new
algorithm or production-scale execution: it retunes the charged cofactor slacks
and the boot redundancy, and it replaces the descent's point-estimate relation
input by a one-standard-error lower value, with every figure recomputed by
evidence16.py.
