# STNFS: 2^120.522000 reference units

Track bls12-381-fixed-g1-invert-g2. Unit: one dense F_p12 multiplication,
C = 20,000 scalar AArch64 instructions. No auxiliary powers, advice, smaller
targets or amortization; no prior verdict is a premise.

## Claim

The algorithm returns Q in G2 with e(P,Q) = z for the exact supplied
nonidentity P in G1 (Q = infinity for z = 1) with success at least 0.39
(modelled >= 0.390236) under heuristics H1-H5 (section 9). A global meter
bounds all work, preprocessing included, by **2^120.522000 units**; the score
is **120.522000**. Storage <= 2^105 bytes (STORAGE.md). Analytic exploratory
claim, not a full-size run. `evidence16.py` (initial run, quick mode) recomputes the pool, yield, cost,
linear algebra, descent, ledger and success and compares each stated bound with
its unrounded value (`declared_bounds`); the filtering rerun, missing logs,
retention and ECM checks ([6], [11], [11b], [12]) run with quick=false.

## 1. Reduction

g = e(P, Q0) (pinned `pairing(Q0, P)`), Q0 the G2 generator: g generates mu_r,
z = g^a for a unique a mod r, and Q = [a]Q0. log z and log g share one
relation and linear-algebra phase; a = log z / log g mod r. g and z enter the
tower of section 2 by the standard field isomorphism (roots of h and of
X^2 + |u| y in F_p12, then a 12 x 12 basis change), charged with the pairings,
[a]Q0 and the final check e(P, Q) = z in the 2^66-unit row.

## 2. Tower and polynomials ([0], [0b], [0c], [13], [14])

h = y^6 + y^3 + y^2 - y - 1 (Disc 58673 squarefree, irreducible mod p);
f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3);
g = x^2 + |u| y, u = -0xd201000000010000; F_p[y,x]/(h, g) = F_p12. sigma:
x -> -x quotients candidates and unknowns. Murphy alpha charged (5.15, 0.85) in
collection and descent; computed (4.868096, 0.843643) over p <= 2000 and
(5.085036, 0.731544) over p < 50000 ([14]).

## 3. Relation collection ([2]-[5], [4a], [12])

Pairs (a, b) in O_h^2, B = 2^69, U = 2^63.451 sigma-orbits. Region
(`config.json`): A^2 = 2^21, a5 > 0, Q_a(a) <= (237/800) A^2, b5 >= 1,
Q_b(b) <= A^2 - Q_a(a); one canonical a per unit class (exact HNF key);
non-primitive pairs excluded. Certified count products at radius >= 788:
0.968192772 (yield), 1.030885630 (work).

**Pool.** 164,000 uniform draws on theta = 5/16 (seeds 5/6; pool of mpjunior92's
dd3785f5, regenerated row by row from (seed, index) by `common12.regenerate_pair`;
[2] checks list membership and 64 field anchors). Its 16,887 canonical rows with
qa <= 237/800 are uniform on the region (effective sample 139,723). Rows with
Q_a(a)/A^2 < 0.02 contribute no yield and their relations are discarded (cells
still sieved and charged).

**Sieve.** Line sieve (scores only) for p in [2^8, 2^48), record walk for
[2^48, 2^64) (KERNELS.md), ECM for [2^64, 2^69] within slacks (295, 226).
Survivors: exact norms, trial division below 2^20, per side two rounds of the
48-bit `ecm16.py` schedule at 99.9% (band strip) and the 69-bit schedule at
99.9%; the g side runs only if f factors (q_f = 100 rho(u_f) / P(omitted_f <=
slack_f), a heuristic estimate).

**Yield and cost.** All relations of the region are kept: >= 2^63.6955 against
the required 1.1 U = 2^63.589, margin 0.1066 >= max(0.06, 2 SE) (pool SE 0.036
bits). Charged 311.7 instructions/cell (allowance **1.0169 >= 1 + 2 x 0.00843**
included) over 2^125.988 cells: collection **<= 2^119.9844 units**. Descent
inputs LOG2_RELS = 63.62 <= estimate - 2 SE (63.6243) and mean norms
739.9/498.9 <= the relation-weighted means 740.3/499.0 ([4]).

## 4. Linear algebra ([6], [10], [13])

Signed sigma-quotient filtering at excess 1.1 (`filtsim.py`; every merge keeps
all explicit coefficients and a bound on the implicit middle-prime entries
<= 255) merged to mean weight 200 (charged 201.1). Runs at 2^12..2^18 (40
systems, `SCALE_NU`, `W200`): weight-100 N/U maxima are non-increasing from
2^15 (0.0618, 0.0592, 0.0580, 0.0529); N(200)/N(100) <= 0.7949. Charged N =
0.80 x 1.066 x 0.05797 x U < 2^59.12 (covers every value from 2^15).
A = [M | S] with 21 Schirokauer columns. Block Wiedemann, blocks (19, 4); per
row and product 201.1 x 20.7 sparse (H4u: >= 88% of weight is +-1) + 21 x 97
dense + row tail 1,640 (ten owners; STORAGE.md). LA **<= 2^117.7461 units**,
generator **<= 2^107.0987 units** (printed 107.10).

## 5-7. Descent and success ([7], DESCENT.md)

62-node plan at LOG2_RELS = 63.62 (calibrated missing-log model at excess 1.1;
missing ideals as prime special-q nodes; dependence margins KAPPA_N = 1.2,
KAPPA_T = 1.27, `toy_descent.py`): two-log work **<= 2^117.1982 units**, meter
2^117.918, P(cap) <= 0.6071577 (Markov). Boot failure <= 1.7171e-04 per log;
node counter 2^27.3, E[N] <= 2^16.88 (abort <= 2^-10.42); repeats <= 2^-9.35;
Block Wiedemann <= 2^-192 (H3); field map < 2^-48. Success >= 1 - 2 x
1.7171e-04 - 0.6071577 - 2^-10.42 - 2^-9.35 - 2^-48 - 2^-192 - 2/r = 0.3902366 - 2/r
>= 0.390236.

## 8. Ledger (log2 units)

| row | bound |
|---|---:|
| relation collection | 119.9844 |
| linear algebra / generator | 117.7461 / 107.0987 |
| two-log meter | 117.918 |
| canonical list / Schirokauer / misc / factor base / certificates | 90 / 88 / 80 / 72 / 50 |
| field map, pairings, verification, output | 66 |

Sum 2^120.5213 < **2^120.522000** (unrounded 2^120.521208).

Fixed rows (1 unit = 2^14.29 instructions): canonical list 2^90 = 2^41.7
instructions per eligible a (2^62.59), against < 2^30 for its exact HNF key and
orbit walk; Schirokauer maps 2^88 = 2^38.6 per relation, against < 2^34.5 for one
exponentiation mod r^2 per side; factor base 2^72 = 2^22.8 per ideal (2^63.45),
against < 2^19 for root finding; misc 2^80 = 2^24.4 per relation entry
(filtering streams, purge and merge bookkeeping, back-substitution), against
< 2^14; certificates 2^50; verification 2^66.

## 9. Heuristics

H1 smoothness (Dickman/Murphy), H1' threshold retention, including random ECM
group orders in collection and the early-abort estimate; H2 uniform pool, 2-SE
yield margin and work allowance; H3 matrix (class numbers, unit maps,
one-dimensional kernel; N/U and the weight-200 ratio not rising at
production); H4 measured kernel counts and H4u +-1 share; H5 descent model
(DESCENT.md), including the missing-log calibration and the dependence margins
KAPPA_N, KAPPA_T measured on the toy descent.

## Provenance

Gautam Botrel with Claude Code and Codex; refinements of 6282ab79 / 3a00ffc8 /
8cc93834 / a31a5f1d lineages. Excess 1.1, 2-SE margin, core cut, the 164k
pool and alpha 5.15 follow mpjunior92's 90d99567 / dd3785f5 (filtering re-run
here with the coefficient cap); Claude Opus 5.5 in Claude Code.
