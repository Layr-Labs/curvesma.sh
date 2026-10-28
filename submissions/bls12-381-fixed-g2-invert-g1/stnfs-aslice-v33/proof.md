# STNFS over a special tower: 2^122.35 reference units

Track bls12-381-fixed-g2-invert-g1. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
advice, smaller targets, future-target amortization or prior results are used.

## Claim

The algorithm below returns P with e(P,Q) = z for the exact supplied nonidentity
Q (P = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins; the modelled success is 0.425. A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^122.35 units**; the score is
**122.35**. The ledger (section 8) totals 2^122.316. Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
The claim is analytic and exploratory under the NFS heuristics H1-H5 (section
9) with the measurements and checks cited; no full-size run is made.

Every number has one source: `model16.py` and `descent16.py` with `config.json`
and `stats16.b85` (LZMA, base85 JSON); `evidence16.py` recomputes every figure.

## 1. Reduction

Let g = e(P0, Q) (pinned convention `pairing(Q, P0)`). Q has prime order r and e
is nondegenerate, so g generates mu_r, z = g^a for a unique a mod r, and
P = [a]P0. log z and log g are computed against one shared relation and
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
are excluded by an exact HNF test. `lattice_count_bounds.py` certifies counts
within [0.98, 1.02] of the volumes; the model scales work by 2601/2500 and yield
by 2401/2500. The canonical list costs 2^90 units.

**Pool.** 20,000 uniform eligible pairs from seeded indices (`pool12-*.b85`),
2,401 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics; [2] regenerates sampled indices and verifies every field.

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20), record
walk for p in [2^20, 2^63) (`slice12.c`), primes in [2^63, 2^69] taken by ECM
within the slacks (292, 220). Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20 and ECM (B1 = 5 x 10^4,
99% per prime, 367.76-instruction Montgomery multiplication, `ecm12.c`): 2^40.1
instructions per survivor (the explicit schedule of `ecm16.py` prices 69-bit
primes in a 300-bit cofactor at 99% at 2^37.4, below the charge).

**Certified block bounds.** At each of the four places v of K_h (rational Rouche
root disks), F_v(w) = sum c_k w^k with w = b_v^2, and b_v ranges over a disk of
radius rho around the block centre (half length, row-centre radius and rounding).
Then |w - w~| <= r_w = rho (2s + rho) + 2^-50 s^2 with s >= |x~| (x~ the computed centre), and the Taylor
shift F_v(w~ + d) = sum T_k d^k gives |F_v - T_0| <=
sum_{k>=1} |T_k| r_w^k + kappa F^(t), with F^ the majorant (|c_k|), t >= |w~| + r_w
and kappa = 2^-32 covering the coefficient intervals (relative width <= 2^-34,
checked) and every rounding of the evaluation (Higham's bounds for the shift are below
2^-45 F^(t)); final values are rounded outward (factors 1 +- 2^-50). G_v = a^2 + u y w is linear in w. The C port (`nbound31.c`,
measurements12.txt) replays the Python leaves bit for bit on all 2,401 pool rows
and costs 1,508 Ir per evaluation (charged 1,550); [2] recomputes the stored
statistics and checks the bound against exact norms.

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
embedded, [6] reruns scale 2^12): N/U = 0.0538 at weight 100; the model charges
N = 1.5 x 0.0538 x U = 2^59.8191 (N/U falls with scale: 0.0696 at 2^12, 0.0538 at
2^17; 1.5 x 0.0538 = 0.0807 covers the 2^12 value). The system solved is the square matrix
A = [M | S]: N + 21 filtered rows (the rest verify the result), the N
ideal columns M and 21 dense Schirokauer columns S. As columns are sigma-quotiented
(log sigma(P) = -log P), only sigma-odd units need maps: signatures K_f (0, 36),
K_f^+ = K_h(x^2) (0, 18), K_g (2, 5), K_h (2, 2) give odd unit ranks 35 - 17 = 18
and 6 - 3 = 3 ([13], Sturm sequences); under H3 the kernel of A is the line of virtual
logarithms (with their Schirokauer coefficients). Block Wiedemann
on A uses random blocks: m = 46 left vectors with entries uniform in [0, 2^64) and
n = 6 right vectors uniform in F_r (the coins). The Krylov phase takes N/m + N/n
block steps, i.e. N(n/m + 1) products by A; a product costs per row 100 sparse
entries x 96 instructions (STORAGE.md; `la_mac.s` inner loop of 33, [10]) plus 21
dense full-field multiply-adds at 104 (`bw_dot.s`, bw_full: 4 x 4-limb product into
a 9-limb accumulator, inner loop of 97). Every step forms the m x n projections
x^T A^i y, (m + n) N^2 dot terms at 40 instructions (`bw_dot.s`: 64-bit x 4-limb
product into a 7-limb lazy accumulator, inner loop of 35, [10]). The kernel vector
is reconstructed by Horner, w <- A w + Y F_i over the deg F <= N/m + N/n generator
coefficients (Y the n starting vectors): one product and n N full-field
multiply-adds per step. The matrix generator (block Berlekamp-Massey,
O((m + n)^3 N log N)) is below 2^90 instructions. With c = 9600 + 21 x 104 per
row, the total is N^2 (c (n/m + 1) + (1/m + 1/n)(c + 104 n) + 40 (m + n))
instructions: **2^119.4650 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** The region and polynomials do not depend on coins or the
  challenge, so the relation count is one fixed number; H1 puts it at
  2^64.219, 0.389 bits above 1.3 U (pool standard error 0.118 bits).
- **Linear algebra.** The matrix is fixed for the target, so its rank is a
  structural property. By Schirokauer's construction the logarithm of
  phi(a - bx) in F_p12 is the sum of the virtual logarithms of the factor-base
  ideals of a - bx and of its Schirokauer-map coordinates, provided (i) r is
  unramified in K_f and K_g, (ii) r divides neither class number and (iii) the
  Schirokauer map is injective on the sigma-odd units of each field modulo r-th powers.
  [13] checks (i): N(f(X)) (degree 72) and N(g(X)) (degree 12) are squarefree
  modulo r and r does not divide disc h = 58673, so r is unramified in Q(x), in
  K_h and in their compositum. H3 is the standard NFS-DL hypothesis: r divides neither h(K_f) nor h(K_g) (a single
  divisibility of a fixed class number by a 255-bit prime; Cohen-Lenstra
  probability of order 1/r), the unit condition (iii) holds (the same kind of
  divisibility), and the filtered matrix over F_r with the Schirokauer columns
  has the line of virtual logarithms as its whole kernel: the excess of 0.3 U
  rows over the N columns leaves no other dependency (a larger kernel would be
  caught by the final pairing check and counted as failure). Under H3, Block Wiedemann on
  A with random blocks fails only if a minor of the projected Krylov matrices
  vanishes. With y uniform in F_r the right Krylov space is generic except with
  probability O(N/r); the left minor has degree <= N in the entries of x (each of
  its rows is linear in one row of x), so with x uniform in [0, 2^64)
  Schwartz-Zippel bounds the failure by N/2^64. The ledger charges 2N/2^64 =
  2^-3.18 = 0.110.
  The solve is verified on the dropped rows; a failure is counted, not retried. The factor-1.5 dimension pad (x2.25 on the charged
  cost) is for N/U alone.
- **Reconstruction.** Logs of ideals removed by singleton purge are recovered by
  back-substitution through the purged relations. The descent counts every
  ideal with at most two occurrences as having no log (section 6); [6b] runs
  the cascading purge and back-substitution on the [6] system and measures the
  ideals left without a log: the model's count is 9-17 times the measured one
  in every band (DESCENT.md).

Under H1 and H3 the main phase succeeds with probability >= 1 - 2^-3.18 = 0.890 over the coins.

## 6. Individual logarithms

DESCENT.md specifies the descent: a boot by ECM trials on balls of certified p^6
lattices, then special-q nodes that sieve an exact lattice ball with the
measured `dsieve16.c` kernel, take relations with a two-tier window of children
and recurse with backtracking over the exact child ideals; survivors are
selected by per-cell certified norm bounds. Only base-primitive
cells are counted (fraction 0.9100; the others repeat smaller relations), and
the retention of relations through the survivor test is computed, not measured:
the omitted factor-base mass of a counted cell has an explicit local
distribution (one Hensel bound per root of the side polynomial at every prime
ideal below 2^9, exact for simple roots), and the loss among
relations weights each mass x by the semismooth probability at size m - x
(`retain16.py`); every sieve node retains at least 0.8592 of its relations
([7]). Cofactors are factored by the explicit ECM schedule of `ecm16.py`
(B1, B2, wheel, curves, per-curve operation count and per-curve probability
from Dickman's rho on the curve order; the same procedure is run on random
36- and 40-bit primes in [12]). Relation counts are negative binomial with
shapes set by the toy counts of ideals without a relation, including empty
windows (DESCENT.md); survivors are bounded at each bin's lowest norms.
Factor-base ideals without a logarithm are handled explicitly: their expected
number per relation follows from the main-phase occurrence model, and each node
either discards relations containing one or descends them as further children
(DESCENT.md). The shipped plan (62 nodes and the boot) gives expected two-log
work **2^118.562 units** ([7]). Node failure uses the negative-binomial shape of
the effective window: the empty-window base nodes take r = 0.5. Records carry a
64-bit owner (index >> 40) and a 40-bit offset; the largest ball has half count
2^102.44, so owners stay below 2^62.44 ([7], STORAGE.md). The descent meter caps
both logs at 2^119.5 units; a run that reaches the cap stops and counts as a
failure, so Markov gives P(cap reached) <= 2^-0.94.

## 7. Success

Success >= (1 - 2^-3.18) x (1 - 2 x 8.94e-05 - 2^-0.94) = **0.425** >= 0.39 (main
phase, two boots, descent cap). The boot term treats the 16 relations' subtrees
as independent (H5); success stays >= 0.39 for any per-log boot failure up to
0.019. Power randomization of the targets makes both
individual-log tasks distributed as the sampled production targets for every Q
and z; the final pairing check rejects any incorrect P.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection | 121.8409 |
| linear algebra | 119.4650 |
| both individual logarithms (metered cap) | 119.5 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| filtering streams, generator, misc | 80 |
| verification and output | 11 |

total_log2W_full_ledger = **122.316011**, under the 122.35 claim.
With every sieve constant x1.5, relation collection is 2^122.4258 and the total
2^122.7593 ([5]); a sensitivity remark, not part of the claim. Nothing is amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (5, 1), exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool, exact flags, certified counts); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
filtering on the production occurrence profile, pad 1.5); H4 kernel cost (KERNELS.md, STORAGE.md);
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs,
local omitted-mass distributions, ECM curve orders as random integers).

## 10. Memory and attribution

Storage <= 2^105 bytes (STORAGE.md). Algorithm, models, kernels and package:
Gautam Botrel with Claude Code and Codex; polynomial selection and alpha follow
Guillevic-Singh, filtering CADO-NFS.
