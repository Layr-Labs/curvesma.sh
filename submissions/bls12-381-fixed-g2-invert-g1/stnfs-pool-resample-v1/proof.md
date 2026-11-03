# Recover P by STNFS on a directly sampled region: 2^121.400000 reference units

Track bls12-381-fixed-g2-invert-g1. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. No
prior verdict is a premise; reused work is credited under Provenance.

## Claim

The algorithm below returns P in G1 with e(P, Q) = z for the exact supplied
nonidentity Q in G2 (P = infinity for z = 1) with success at least 0.39 over
uniform full-production challenges and fresh coins under the declared
heuristics H1-H5; the modelled success is at least 0.392624. A global
instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, field conversion, verification and output, by
**2^121.400000 units**; the score is **121.400000**. The upper-rounded ledger
rows (section 8) sum to less than 2^121.400000. Distributed storage is at most
2^105 bytes (STORAGE.md). The claim is analytic and exploratory under H1-H5
(section 9) with the measurements and checks cited; no full-size run is made.

**Numerical conventions.** Figures used as bounds carry <= or >= and are
rounded in the safe direction; the `declared_bounds` check of `evidence16.py`
compares every stated bound with the unrounded value. The model is
`model16.py` and `descent16.py` with `config.json` and `stats16.b85`;
`evidence16.py` recomputes the ledger and the checks cited as [n].

## 1. Reduction

Let g = e(P0, Q) (pinned `pairing(Q, P0)`), one pairing; P0 is the fixed G1
generator. Q has prime order r and e is nondegenerate, so g generates mu_r,
z = g^a for a unique a mod r (a = 0 iff z = 1), and P = [a]P0 (one G1 scalar
multiplication) satisfies e(P, Q) = z for that exact Q. log z and log g share
one relation/LA phase; a = log z / log g mod r. Q uniform nonidentity makes g
uniform on mu_r minus {1}, independent of z.

**Field conversion.** The core works in F_p[y, x]/(h(y), x^2 + |u| y)
(section 2). In the protected F_p12 encoding, split h into linear factors and
take a root y, then a root x of X^2 + |u| y (equal-degree splitting with random
gcds; each trial splits with probability > 1/4, 128 trials per split, failure
< 2^-48 in total); Y -> y, X -> x is a field isomorphism, and inverting the
12 x 12 coordinate matrix of the basis y^i x^j maps g and z. A field
isomorphism preserves z = g^a. This, both pairings, the input checks, [a]P0
and the final check e(P, Q) = z cost < 2^80 instructions, charged as the
2^66-unit row of section 8; the 2^-48 failure is inside the success slack.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks that
(u y)^((p^6-1)/2) = -1 in F_p6, so g is irreducible over F_p6 and
F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. The map
x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P), so
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at
(5, 0.85) in collection: [14] recomputes it from the exact Kummer-Dedekind
structure, (4.868096, 0.843643) over p <= 2000 and 0.731544 for g over
p < 50000 (p = 2 double roots of g above q = 4 and q = 16: v2(U) = 16 and y is
a unit, so for b a unit, probability q/(q+1) under base-primitivity, the term
is at most 2q/(q^2-1)(1+q^-8)). The descent keeps (5, 1).

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69.
The factor base is the relative-degree-1 prime ideals over K_h of norm <= B on
both sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21, theta = 2/5: eligible a: a5 > 0,
Q_a(a) <= (2/5) A^2; half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The canonical
a of each unit class is the least of its exact HNF key class in (Q_a, lex)
order, so each relation orbit of <units, sigma> has exactly one point in the
domain; non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` certifies from exact rational Machin-pi and
determinant/sqrt intervals the count ratio at every radius R >= 915 =
isqrt((2/5) 2^21), the smallest a-radius (b-radii are at least 1121):
[0.986205, 1.013186] for Q_a and [0.986160, 1.013186] for Q_b. The model uses
the exact products, **0.972556868 (yield, floor)** and **1.026547145 (work,
ceiling)**. The canonical list costs 2^90 units.

**Pool.** The pool is sampled directly on this region: 32,000 uniform eligible
pairs from seeded indices (seed 1: 2,000; seed 2: 30,000; `pool12-*.b85`),
of which 3,907 are canonical rows (3,589 primitive), stored with exact norms,
primitivity flags and certified norm-bound statistics. Every draw lies in the
collected region, so the estimator scales by va/32,000 with no sub-region
incidence loss; [2] regenerates sampled indices and verifies every field. The
region parameters (theta, slacks) were fixed from separate exploratory pools
before this pool was drawn, and the shipped pool was evaluated once.

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20),
record walk for p in [2^20, 2^64) (`slice12`, KERNELS.md; X2 = 2^64), primes in
[2^64, 2^69] taken by ECM within the slacks (295, 224): survivors' cofactors
are bounded by 303 (f) and 232 (g) bits. The per-survivor ECM charge is the
`ecm16.py` schedule at 99.9% per prime (ECM_FIND = 0.999, the same figure as
the yield's miss model), checked by [12]: B1 = 2^16, B2 = 30 B1, wheel 30030,
403 curves, five-limb Montgomery multiplication at 367.76 Ir (`ecm12`):
2^38.0906 and 2^38.0878 instructions per side, per-survivor charge
C_NORMS + C_TDIV + both = 2^39.0892. Certified per-block bounds
(`embedding_bounds.py`, 64-cell blocks of width 128, bisected until both
bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining
every B-smooth cell with omitted mass below slack and bounding survivors'
cofactors by slack + 8 bits.

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre in double precision; the remainder over the block disk plus a
2^-32 relative term is added and results are rounded outward. `nbound31`
replays the Python leaves bit for bit (offline, 2,401 rows of the base pool)
at 1,508.2 Ir per evaluation (charged 1,508.2); [2] checks the bounds against
exact norms. Every block of every row is enumerated, sieved and charged.

**Yield ([4]).** rels = (va/32000) vb (0.972556868) sum over canonical
primitive rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g),
with Dickman rho ([1]), retentions r from the finite Euler-product threshold
model (`threshold_model.py`, calibrated on exact smooth integers at
B = 2^10..2^14) and the ECM miss. Retained relations >= 2^64.3065 against
1.3 U = 2^63.830: margin >= 0.4766 bits (pool SE 0.083), above the 0.30-bit
floor. The descent input `LOG2_RELS = 64.306` is at or below this count and
the descent mean norms 744.2/499.2 are at or below the relation-weighted means
(744.3/499.3); [4] checks both.

**Cost ([5]).** Sampled mean 395.2 instructions/cell; charged mean 408.4 with
the allowance below, over 2^126.863 cells: all-cell bounds, layout, setups,
line path, record walk (primes in [2^63, 2^64) charged twice), routing
(10.01/record, 8/cell; STORAGE.md), apply 7.00, recovery 7.003/record, scan,
zero, cofactorisation. Collection: **<= 2^121.2492 units**.

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo
estimate of a fixed, challenge-independent domain total. Its relative SE has a
conditional row-mean term (0.730%) and an incidence term for the 3,907
canonical draws of 32,000 (sqrt((1 - 3907/32000)/3907) = 1.499%): total
<= 1.668%. The charge multiplies the estimate by **1.0334 >= 1 + 2.00 x
0.01668** ([4a] fails otherwise). This is a declared modelling margin under
H2, not a failure probability.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge (coefficients <= 255), simulated as
the signed sigma-quotient system (`filtsim.py`, kmax 64, weight 100; [6] reruns
2^12). The embedded N/U table (`SCALE_NU`, max over seeds) is 0.0850 (2^12),
0.0635 (2^13), 0.0563 (2^14), 0.0516 (2^15), 0.0469 (2^16), 0.0461 (2^17),
0.0425 (2^18) at relation norm means 748.39/498.93 (`unit_entries.py` uses the
same inputs); at this region's relation-weighted means 744.3/499.3 five seeds
per scale give at most 0.0491 (2^15), 0.0478 (2^16), 0.0460 (2^17) (offline,
same code). The model charges N = 1.096 x 0.0471 x U = 2^59.1758: ratio
0.05163 exceeds every value from 2^15 up; N/U not rising at production scale
is part of H3. Mean row weight is charged at 101.

The square matrix A = [M | S] has N + 21 filtered rows (the rest verify), N
ideal columns and 21 dense Schirokauer columns: signatures K_f (0, 36),
K_f^+ (0, 18), K_g (2, 5), K_h (2, 2) give sigma-odd unit ranks 18 and 3
([13], Sturm); under H3 the kernel of A is the line of virtual logs. Block
Wiedemann uses random blocks, m = 19 left and n = 4 right vectors uniform in
F_r: N(n/m + 1) products by A. Per row and product: 101 x 20.7 sparse
(STORAGE.md: 4 + 3 + 0.88 x 11 + 0.12 x 33 = 20.64; +-1 entries add without
multiplying, >= 88% of row weight is +-1 by H4u, measured 0.917..0.946), 21
dense terms at 97, and row tail C_ROW = 1,640. Projections x^T A^i y cost
(m + n) N^2 dense terms at 97 plus one reduction (C_RED9 = 420, counted 385)
per entry; Horner reconstruction costs one product plus, per element, n dense
terms and one reduction. With c = 101 x 20.7 + 1,640 + 21 x 97:
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n + 420) + 97 (m + n)) + 420 (m + n) N
instructions: **<= 2^117.5157 units**. The (m + n) x (m + n) generator of the
length L = (N + 21)(1/m + 1/n) + 2 sequence (divide-and-conquer approximant
bases, Karatsuba): <= 4 (m + n)^3 L^log2(3) log2(L) field multiplications at
97 plus reduction 420: **<= 2^107.1997 units**.

**Row tail (ten active owners).** N < 2^59.176 < 10 x 2^56: column j goes to
block j >> 56, and only blocks 0..9 are populated. Each row uses a fixed tree
on ten owners (nine sends and nine combinations; owners without entries send
zero partial sums): 9(15 + 18) + 9(62 + 2) + (92 + 180 + 183 + 2) + 100 + 40
+ 168 = 1,638 <= 1,640 instructions, from the counted loop-free kernels of
`la_mac.s` ([10]; `owner10.py` checks the allocation and formula). Column
blocks keep width 2^56, so the five-limb group and six-limb message bounds
hold; the row owner evaluates all 21 dense terms and the final reduction.

## 5. Main-phase success

- **Relations.** Region and polynomials are coin/challenge-independent, so the
  count is one fixed number; H1 puts it at >= 2^64.3065, 0.4766 bits above
  1.3 U.
- **Linear algebra.** H3 (standard NFS-DL): r is unramified in K_f and K_g
  ([13]), divides neither class number, the 21 sigma-odd Schirokauer maps are
  injective on the units modulo r-th powers, and the kernel of A is the line of
  virtual logarithms. Block Wiedemann with uniform blocks then fails with
  probability <= 4N/r < 2^-192 (Schwartz-Zippel). A failed solve (checked on
  the dropped rows) counts as failure.
- **Reconstruction.** Purged ideals' logs are recovered by back-substitution.
  The descent counts every ideal with <= 2 occurrences as log-less; [6b]
  measures ideals left without a log, below the model's count in every band.

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6
lattices, then special-q nodes sieving an exact lattice ball (`dsieve16`) and
recursing on child ideals. Counted cells are base-primitive (0.9100);
survivor-test retention is >= 0.8592 at every sieve node ([7]); cofactors use
the `ecm16.py` schedule ([12]). Relation counts are negative binomial.
Log-less ideals are discarded or re-descended in a fresh lattice. The 62-node
plan and boot (REFINEMENT.md) give expected two-log work **<= 2^115.6681
units** ([7]; one log <= 2^128.956 instructions). The meter caps both logs at
2^116.3915 units (cap reached = failure), so P(cap) <= E[W]/cap <= 0.605669
(Markov).

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 6.7554e-05 - 0.605669 - 2^-10.01 -
2^-10.70 - 2/r) >= **0.392624** > 0.39 + 2^-48: main phase, two boots
(phi <= 6.7554e-05 each), the instruction-meter Markov bound, the node-counter
Markov bound and special-q repeats.

**Two explicit meters.** Both logs run under (i) the instruction meter (cap
2^116.3915 units) and (ii) a special-q node counter (cap N_MAX = 2^27.3) that
increments before each special-q node and aborts as a failure past N_MAX.

- Node counter: P(abort) <= E[N]/N_MAX. E[N] uses the finite-tree recurrence
  of E[W] with unit node cost (`counts16.py`): E[N] <= 2^17.29 for both logs,
  so P(abort) <= 2^-10.01.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for
  1/p weighting share a special-q with probability <= 2^(2 x 27.3 - 1 + 6 -
  70.30) <= 2^-10.70, charged in full.

The boot term treats the boot relations' subtrees as independent (H5).
Identity z is handled directly; the field-conversion failure 2^-48 and the
2/r distance from the uniform-target model are charged; the final pairing
check rejects any incorrect P. The H2 allowance is a resource margin, not a
probability term.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (estimate x 1.0334 allowance) | <= 121.2492 |
| linear algebra, ten active vector owners | <= 117.5157 |
| matrix generator | <= 107.1997 |
| both individual logarithms (metered cap) | 116.3915 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 72 |
| verification, field conversion, pairings, output | 66 |
| numerical certificate construction | 50 |

The rows sum to 2^121.39941 < 2^121.400000 (`audit_refinement.py` recomputes
the sum from refinement.json; `evidence16.py` gives the unrounded
2^121.399351). Preprocessing is inside the cap; nothing is amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; collection alpha (5, 0.85) against the computed
(4.868096, 0.843643) and 0.731544 for g ([14]); descent alpha (5, 1); exact
sampled norms); H1' retention (Euler-product model calibrated on exact smooth
integers); H2 region (uniform pool sampled on the theta = 2/5 region, exact
flags, certified counts with exact rational products, and the canonical work
mean charged with a 3.34% allowance = 2.00 incidence-inclusive standard errors,
a modelling margin, not a probability term); H3 matrix (r unramified checked
by [13]; r not dividing h(K_f) or h(K_g), the unit condition and the
one-dimensional kernel assumed; signed filtering, pad 1.096 over the 2^17
table, N/U assumed not to rise at production); H4 kernel cost (KERNELS.md,
STORAGE.md; LA tail and reductions statically counted loop-free code, the
route an executed-instruction count), with H4u: >= 88% of production row
weight is +-1 (measured 0.917..0.946 at 2^12..2^15); H5 descent (DESCENT.md:
lower sums, negative-binomial recurrence, missing logs, independent subtrees,
fresh-lattice re-descents, local omitted-mass laws, random-integer ECM orders,
UDUP = 0.5 against measured 0.779).

## 10. Memory

Storage <= 2^105 bytes (STORAGE.md), plus < 2^30 bytes for the field maps.

## Provenance

GT discrete-log machinery: Gautam Botrel with Claude Code and Codex (this
track's 3cb69d49 lineage); polynomial selection and alpha follow
Guillevic-Singh, filtering CADO-NFS. ECM 99.9%, X2 = 2^64, slacks, alpha_g
0.85, record route, LA kernels, pad, ten-owner tail, quadrature and plan:
sister-track bls12-381-fixed-g1-invert-g2 lineage (6282ab79, 3a00ffc8,
e3940d11, 04073dd1, 8c1635f2). New here (Claude Opus 5.5 in Claude Code): the
directly sampled theta = 2/5 pool, the 1.0334 allowance, re-derived descent
inputs and meter, and the recover-P wrapper with field conversion.
