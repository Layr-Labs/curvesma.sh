# STNFS over a special tower: 2^121.993 reference units

Track bls12-381-fixed-g2-invert-g1. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. The argument reuses predecessor work as credited below, without assuming a prior-results verdict as a premise.

## Claim

The algorithm below returns P with e(P,Q) = z for the exact supplied nonidentity
Q (P = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins. The declared success is **0.5025**,
the descent and linear-algebra conditional bound at the relation input 2^64.130:
the main phase fails only through the Wiedemann projections (<= 2^-192), the two
boots through psi^8 = 2.1438e-02 each, the descent cap with directed Markov
2^-1.189, and the special-q repeats with 2^-5.969. No yield-survival probability
is claimed: the relation input is the epistemic haircut
floor(2^64.2541 - 1.05 x 0.1175) of section 3, and the success is conditional on
the fixed production supply reaching it (disclosed in H6). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^121.993 units**; the score is
**121.993**. The ledger (section 8) totals 2^121.992667. Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
The claim is analytic under H1-H5 and the disclosed H6 input (section 9); no full-size run is made.

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

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^40), record
walk for p in [2^40, 2^63) (`slice12`, KERNELS.md), primes in [2^63, 2^69] taken by ECM
within the slacks (295, 226). The f-side slack is 295, so survivors'
f-side cofactors are bounded by 303 bits; the g-side slack is 226, so survivors'
g-side cofactors are bounded by 234 bits. The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at 99% per prime, the same schedule the
descent uses and [12] checks: for 69-bit primes in a 303-bit (resp. 234-bit)
cofactor it prices one side at 2^37.5177 instructions and the other at 2^37.5150,
so the per-survivor charge
C_NORMS + 600 pi(2^34) + 2^37.5177 + 2^37.5150 is 2^39.632; trial division covers
every prime below 2^34 at a charged 600 instructions per prime, with
pi(2^34) <= 2^34/(ln 2^34 - 1.1) = 2^29.510 by the Rosser-Schoenfeld bound, so the
trial-division term is 2^38.739 per survivor. This is not weaker than the
predecessor's charge: at the 2^20 bound the same rate is 2^25.55, above the
accepted 2^24.5 lump sum. (The predecessor's crude
ECM table charged 2^40.06 per survivor, 2^1.545 above the package's own
schedule; the 99% figure is the same ECM_FIND used in the yield's miss model).
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre in double precision, with an outward-rounded 2^-32 relative term;
`nbound31` replays the Python leaves bit for bit on all 2,401 pool rows at
1,508.2 Ir per evaluation (the measured marginal per-block cost; per-row and
per-slice setups are charged separately), and [2] checks the bounds against exact
norms. Every block of every row is enumerated and charged at the full per-cell cost.

**Yield ([4]).** rels = (va/20000) vb (2401/2500) sum over canonical primitive
rows of w rho(u_f) rho(u_g) r_f r_g 0.99^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
and the ECM miss. Retained relations 2^64.254 against 1.3 U = 2^63.830: margin
0.424 bits (pool SE 0.1175), above the model's 0.30-bit feasibility floor
(`MARGIN` in model16.py). The descent is not run at the point estimate. Its
missing-log input is `LOG2_RELS = 64.130 = floor(2^64.2541 - 1.05 x 0.1175)`,
i.e. 1.05 estimated standard errors below the pool mean; [4] checks this
inequality against the estimate and its standard error. This is an epistemic
haircut, not a probability: the 20,000-pair pool estimates a fixed production
quantity, so the 0.1175-bit standard error measures estimation uncertainty and
is not assigned a production tail probability anywhere in the success bound.
Section 7's declared success is conditional on the fixed supply meeting the
input; the residual risk in that estimate is disclosed as H6 and is not
multiplied into the success; the normal-tail multiplier Phi(1.00) = 0.8413 is
withdrawn.

**Cost ([5]).** 332.10 instructions per cell over 2^127.389 cells: all-cell bounds,
layout and setups, the line path, record walk, routing
(36 per record, 8 per cell; STORAGE.md), apply 7.00, recovery 7.003, scan, zero
and cofactorisation. Relation collection: **2^121.4767 units**. The line/record
split is retuned from P_MID = 20 to 40. At the new boundary the production
line-coset size is S/(Kp) >= 2^14.03 cells (minimum over the 2401 pool rows,
[5c]), i.e. 0.21 bits above the measured table's top anchor h = 2^13.82, so
every charged line hit is at the measured minimum 24.10 instructions and no
point is interpolated below the table. The record walk is charged at the
package's own measured argument `h = 2 cells / p` (section B of
`measurements12.txt`), i.e. `log2h = log2(2S) - lp`; [5b] checks the five
measured anchors. Its top-table value is 44.90 per hit; rows whose argument
falls below the table top are interpolated upward by the h-monotone table, and
the packaged mean is 77.9 instructions per cell for this component, so the
charge is at least the top value per hit.

**Line-factor recovery.** Trial division covers every prime p < 2^34 at 600
instructions per prime, with pi(2^34) <= 2^34/(ln 2^34 - 1.1) (Rosser-Schoenfeld),
i.e. 2^38.739 per survivor. A line-sieved hit with p >= 2^34 is emitted as a
record (cell, prime) and charged the package's measured 36 routing + 7.00 apply +
7.003 recovery; the prime is divided out of the norm before ECM. The total record
count is NBAR ln(2^63/2^34) = 1.2335 per cell (line 0.3250 + record walk 0.9085),
charged 1.2335 x 50.003 + 8 = 69.68 instructions per cell; [5d] checks that
apply, recovery and routing are charged on that full count. No line-sieved factor
above the trial-division bound is credited without the record that recovers it.

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
  challenge, so the relation count has one fixed mean; H1 puts the retained
  count at 2^64.254, 0.424 bits above 1.3 U (pool standard error 0.1175 bits).
  The descent input is 64.130 = floor(estimate - 1.05 SE); its yield tail is
  charged by multiplication in section 7 (H6).
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
lattice (DESCENT.md). The shipped plan (62 nodes and the boot, 8 kept boot
relations) gives expected two-log work **2^117.951 units** ([7]), priced at the
relation input `LOG2_RELS = 64.130` of section 3; record layout in STORAGE.md. The descent meter caps
both logs at 2^119.14 units; a run that reaches the cap stops and counts as a
failure, so Markov gives P(cap reached) <= 2^-1.189 (the experiment's directed
rounding of the exact model value 2^-1.18908 = 0.438583).

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 2.1438e-02 - 2^-1.189 - 2^-5.969) = **0.5025** >= 0.39
(main phase, two boots with psi^8 = 2.1438e-02 each, descent cap with directed
Markov 2^-1.189, special-q repeat bound 2^-5.969). No yield-survival factor is
multiplied in: the relation input of section 3 is an epistemic haircut (H6) and
the bound is conditional on the fixed production supply reaching it. The boot
term treats the 8 relations' subtrees as independent (H5). The 0.1125 margin
over the 0.39 floor is the difference between the charged bound and the floor;
this revision claims no sensitivity percentage. Power randomization gives uniform nonidentity targets; the distance from the
uniform full-group model is at most 1/r per target and 2/r is charged. Identity
z is handled separately; the final pairing check rejects any incorrect P.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection | 121.4767 |
| linear algebra | 119.3684 |
| matrix generator | 106.42 |
| both individual logarithms (metered cap) | 119.14 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| filtering streams, misc | 80 |
| verification and output (pairing, G1 scalar multiplication, check) | 11 |

total_log2W_full_ledger = **121.992667**, under the 121.993 claim. The displayed
rounded component exponents sum to 2^121.992626 (from 2^121.4767 + 2^119.3684 +
2^119.14 + smaller terms), 0.000374 bits below the declared 121.993 ceiling;
the exact model total is 2^121.992667, 0.000333 bits below the ceiling. The work
and yield brackets (2601/2500 and 2401/2500) are exactly the accepted frontier's.
Relative to the accepted 122.2065 frontier (cbab3319): the record path uses the
measured argument `h = 2 cells / p` ([5b]); P_MID = 40 with the line-sieved
factors above the trial-division bound recovered by charged records of section 3
([5c],[5d]); the trial-division bound is 2^34 at 600 Ir/prime (section 3), so
the charged record count is 1.2335 per cell and relation collection is
2^121.4767 units at 332.10 per cell; slacks (295, 226) give yield 2^64.254,
0.424 bits above 1.3 U (SE 0.1175); the descent input is
`LOG2_RELS = 64.130` ([4]) with no yield probability charged (H6); boot
(210, 41, 16, 4) gives two-log 2^117.951, cap 2^119.14 (directed Markov
2^-1.189, exact 2^-1.18908), psi^8 = 2.1438e-02, declared success 0.5025.
Ledger 2^121.992667, ceiling 121.993. The other components are the accepted
122.2065 tree's (ECM schedule, norm-bound charge 1,508.2 Ir, linear algebra
(26,5) at 88/97 with pad 1.61, generator, 62-node plan). Preprocessing is inside
the 2^121.993-unit cap; nothing is amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (5, 1), exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool, exact flags, certified counts); H3 matrix (Schirokauer
premises: r unramified checked by [13]; unit condition and one-dimensional
kernel assumed, section 5; signed filtering on the production profile, pad 1.61); H4 kernel cost (KERNELS.md, STORAGE.md);
H6 relation-supply input (section 3: the pool mean minus 1.05 estimated standard
errors, an epistemic haircut, not a probability; the declared success is
conditional on the fixed supply reaching it; no normal-tail factor is charged);
H5 descent (DESCENT.md:H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing
logs, independent subtree outcomes, re-descended ideals, omitted-mass
distributions, ECM curve orders as random integers). H5 is a
disclosed heuristic, not derived: section 7 charges its two consequences
explicitly (the psi^8 boot failure and the cap) and keeps a 0.1125 success margin
over the 0.39 floor for it.

## 10. Memory and attribution

Storage <= 2^105 bytes (STORAGE.md). Algorithm and package: Gautam Botrel with
Claude Code and Codex; polynomial selection and alpha follow Guillevic-Singh,
filtering CADO-NFS.

## Provenance of this revision

This revision derives from the accepted 122.2065 tree `stnfs-lean-122206-v1`
(cbab3319) and its ancestors aaa530bf, 0ac63018, fc9fd9a0, a73abc31 and a57c5623,
which supply the measured charges, the work/yield brackets and the metered-cap
method. It supplies the line-factor recovery whose absence rejected 78f869d4
([5d]) and replaces the point-estimate relation input by a charged 1.05-SE
quantile (H6); the rejected lattice-count re-accounting ac905f62 is not used.
Every figure is recomputed by evidence16.py; no production-scale run is claimed.
