# Recover P by STNFS on a directly sampled region: 2^121.400000 units

Track bls12-381-fixed-g2-invert-g1. Unit: one dense F_p12 multiplication,
C = 20,000 scalar AArch64 instructions. No auxiliary powers, advice, smaller
targets or amortization.

## Claim

The algorithm returns P in G1 with e(P, Q) = z for the exact supplied
nonidentity Q (P = infinity for z = 1) with success >= 0.392624 > 0.39 over
uniform challenges and fresh coins, under heuristics H1-H5 (section 9). A
global instruction meter bounds all work (preprocessing, failures, field
conversion, verification, output) by **2^121.400000 units**; score 121.400000.
Storage <= 2^105 bytes. Analytic exploratory claim, not a full-size run.
Stated bounds are rounded in the safe direction; `evidence16.py` recomputes
every figure (check `declared_bounds`), `audit_refinement.py` the ledger sum
and success from refinement.json.

## 1. Reduction and field conversion

g = e(P0, Q) (pinned `pairing(Q, P0)`), P0 the G1 generator. g generates mu_r,
z = g^a for a unique a mod r, and P = [a]P0 satisfies e(P, Q) = z. Q uniform
nonidentity makes g uniform on mu_r minus {1}. Both logs share one relation and
linear-algebra phase; a = log z / log g mod r. The protected F_p12 encoding is
mapped to the core tower by a field isomorphism: roots y of h and x of
X^2 + |u| y in F_p12 by random-gcd splitting (128 trials per split, failure
< 2^-48), then the inverse of the 12 x 12 matrix of the basis y^i x^j. This,
both pairings, [a]P0 and the final check e(P, Q) = z cost < 2^80 instructions
(ledger row 66).

## 2. Polynomials ([0], [0b], [0c], [14])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1 (Disc 58673 squarefree, irreducible
mod p). f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3),
g = x^2 + |u| y, u = -0xd201000000010000; F_p[y,x]/(h, g) = F_p12 ([0b]) and
Res_y(h, f(x^2 -> u y)) = 0 mod p ([0]). sigma: x -> -x (log sigma(P) =
-log P) quotients candidates and unknowns. Murphy alpha charged (5, 0.85)
against computed (4.868096, 0.843643), and 0.731544 for g over p < 50000 ([14]).

## 3. Relation collection ([2]-[5])

(a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69; U = 2^63.451
sigma-orbits. Region (`config.json`): A^2 = 2^21, theta = 2/5, a5 > 0,
Q_a(a) <= (2/5) A^2, b5 >= 1, Q_b(b) <= A^2 - Q_a(a); canonical a per unit
class (exact HNF key), non-primitive pairs excluded. Certified count products
at every radius >= 915: 0.972556868 (yield, floor), 1.026547145 (work,
ceiling) (`lattice_count_bounds.py`). Canonical list 2^90 units.

**Pool.** 32,000 uniform draws on this region from seeded indices (seed 1:
2,000, seed 2: 30,000), 3,907 canonical rows (3,589 primitive) with exact norms
and certified norm-bound statistics (`pool12-*.b85`; [2] regenerates seed-1
draws). Theta and slacks were fixed on separate exploratory pools before this
pool was drawn.

**Sieve.** Line sieve p in [2^8, 2^20), record walk p in [2^20, 2^64), ECM for
[2^64, 2^69] within slacks (295, 224) at 99.9% per prime (`ecm16.py`, [12]):
2^39.0892 instructions per survivor. Certified per-block norm bounds
(`embedding_bounds.py`) set thresholds; every block is enumerated and charged
(KERNELS.md, STORAGE.md).

**Yield ([4]).** Dickman rho, Euler-product retentions (`threshold_model.py`)
and the ECM miss give >= 2^64.3065 relations against 1.3 U = 2^63.8299: margin
>= 0.4766 bits (pool SE 0.083). Descent inputs LOG2_RELS = 64.306 and mean
norms 744.2/499.2 are at or below the pool's values.

**Cost ([5], [4a]).** Charged 408.4 instructions/cell over 2^126.863 cells
(sampled 395.2 times allowance 1.0334): **<= 2^121.2492 units**. The allowance
is >= 1 + 2.00 x the work estimator's relative SE (conditional 0.730%,
incidence 1.499%, total <= 1.668%), a modelling margin (H2), not a probability.

## 4. Linear algebra ([6], [10], [13])

Signed sigma-quotient filtering (`filtsim.py`): the embedded N/U maxima are
0.0516 (2^15), 0.0469 (2^16), 0.0461 (2^17), 0.0425 (2^18), falling with
scale; charged N = 1.096 x 0.0471 x U < 2^59.176 (ratio <= 0.05163), mean weight
101. A = [M | S] with 21 Schirokauer columns (sigma-odd unit ranks 18 + 3,
[13]). Block Wiedemann, blocks (19, 4) uniform in F_r; per row and product
101 x 20.7 sparse (>= 88% of weight is +-1, H4u) + 21 x 97 dense + row tail
1,640 (ten owners, counted 1,638, `owner10.py`); projections 97 per term plus
420 per reduction. LA **<= 2^117.5157 units**; generator (approximant bases,
Karatsuba) **<= 2^107.1997 units** (STORAGE.md).

## 5-7. Descent and success ([7], DESCENT.md)

Boot plus 62-node special-q plan: expected two-log work **<= 2^115.6681
units**; meter cap 2^116.3915, so P(cap) <= 0.605669 (Markov). Boot failure
<= 6.7554e-05 per log; node counter cap 2^27.3 with E[N] <= 2^17.29, abort
<= 2^-10.01; special-q repeats <= 2^-10.70. Main phase fails only by the Block
Wiedemann projections (<= 2^-192) under H3. Success >= (1 - 2^-192)(1 -
2 x 6.7554e-05 - 0.605669 - 2^-10.01 - 2^-10.70 - 2/r) - 2^-48 >= 0.392624.

## 8. Ledger (log2 units)

| row | bound |
|---|---:|
| relation collection | 121.2492 |
| linear algebra | 117.5157 |
| generator | 107.1997 |
| two-log meter | 116.3915 |
| canonical list / Schirokauer / misc / factor base | 90 / 88 / 80 / 72 |
| field conversion, pairings, verification, output | 66 |
| certificate construction | 50 |

Sum 2^121.39941 < 2^121.400000 (unrounded 2^121.399351).

## 9. Heuristics

H1 smoothness and Murphy alpha; H1' threshold retention model; H2 region,
uniform pool and the work allowance; H3 matrix (r unramified [13]; class
numbers, unit maps and one-dimensional kernel assumed; N/U not rising at
production); H4 kernel counts, H4u +-1 share (measured 0.917..0.946 offline);
H5 descent model (DESCENT.md).

## Provenance

Core: Gautam Botrel with Claude Code and Codex (3cb69d49 lineage) and the
sister-track refinements of 6282ab79 / 3a00ffc8 / e3940d11 / 04073dd1 /
8c1635f2. New here (Claude Opus 5.5, Claude Code): direct theta = 2/5 pool,
1.0334 allowance, re-derived descent inputs and meter, recover-P wrapper.
