# Multilevel descents under the production node rules

49200 recursive descents on the reduced tower over the production base ring (158215 entered nodes, 107502 adaptively
reached prime-ideal nodes, depth up to 5), run by the node procedure of DESCENT.md. It adds three measurements to the
171-node recursive campaign (part=adaptive_stats) and supports the sequential hazard and the conditional child-work envelope
at toy scale. It does not establish them at production size.

## 1. Attempts after an earlier failed relation

| campaign | prior failed relations | completed | failed | censored | failure rate | Clopper-Pearson 95% upper |
|---|---|---|---|---|---|---|
| A | 0 | 155293 | 2078 | 0 | 0.0134 | 0.0139 |
| A | 1 | 2025 | 40 | 0 | 0.0198 | 0.0256 |
| A | 2 | 40 | 0 | 0 | 0.0000 | 0.0722 |

Campaign A: 2065 completed attempts after an earlier failed relation, 40 failed, 0 censored. At the first-attempt failure rate of the same parent type 42.5 failures were expected against 40 observed: ratio 0.94 [0.65, 1.24].

## 2. Deterministic-order roots

Campaign A: all 49200 runs use the production order; none is shuffled. 0 runs hit the node budget and 0 attempts are unresolved: nothing is censored. Roots are uniformly sampled ideals of every type and both sides (38000 prime-ideal roots, 11200 helper-control roots). Raw run files are not shipped (sha256 a57288bff60b8258); part=toy regenerates any run.

## 3. Conditional child work after a failure

Children entered after a failed relation of their parent against (a) fresh ideals of the same class (bits, side, half-bit) and (b) the slot envelope that child_work charges (largest fresh mean over the types admitted in the child's window slot); ratios of sums, 95% interval or upper limit:

| campaign | children | subtree nodes vs fresh | subtree work vs fresh | failure vs fresh | subtree nodes vs envelope | subtree work vs envelope |
|---|---|---|---|---|---|---|
| A, after a failure | 2373 | 1.03 [1.01, 1.06] | 1.02 [1.00, 1.05] | 0.85 [0.57, 1.12] | 0.86 (<= 0.89) | 0.62 (<= 0.68) |
| A, no failure on the path | 103452 | 1.00 [1.00, 1.01] | 1.00 [0.99, 1.00] | 0.93 [0.81, 1.05] | 0.85 (<= 0.85) | 0.59 (<= 0.60) |

Work is the number of survivors factored in the whole subtree.

## 4. Correspondence with production

| rule | production (DESCENT.md) | campaign (mtoy.py in libextra.b85) |
|---|---|---|
| polynomials | f of degree 8 and g of degree 2 in x over the degree-6 base O_h | f_t = x^4 + w x^2 + w^2 and g_t = x^2 + 4w over the same O_h (the tower of the 862-ball calibration) |
| lattice, ball, traversal | special-q lattice, reduced; radius-R ball, one cell per +- pair, cell order; first success ends the node | same routines and order; squared radius 3 to 4 |
| windows | <= kA primes in (2^56, 2^nA], <= kB in (2^nA, 2^n], row per (bin, side) | the same predicate above 2^lb (lb = 14), row (n, iA, kA, kB) per (bits, side); model scale by the same ss_lower2 |
| known-log rule | factor base to 2^56; ideals whose column occurs <= 2 times have no logarithm; thin rows reject relations containing one, other rows descend it by a helper node on Q q' under the thin rule | factor base to 2^lb; a hash marks ideals of norm >= 2^8 as missing, more often when larger; same thin rule and helper nodes |
| repeat handling | nothing solved earlier is reused; a repeated ideal aborts the run | nothing is reused; a repeated ideal is entered as a flagged node with its first outcome, so the run up to its first repeat is the production process |

## 5. What remains different, and what this does not show

- The toy f has degree 4 in x where production's has degree 8; g has degree 2 in both. Ideals have 15 to 26 bits against 56 to
  192. A toy ball holds about 3 to 15 relations against hundreds, so toy failure rates are far above the production bounds: a
  stress test of selection, not a measurement of production rates.
- Factorisation is exact (no ECM). The set of ideals without logarithm is a hash; with one fixed set there are few of them, so
  reached helper nodes repeat them across runs; they are kept out of the comparisons (part=multilevel lists them).
- Nothing here proves that the ratios of section 3 hold at production size; at this scale, under the production rules,
  reaching an ideal through earlier failures changes its failure rate and subtree work by a few percent at most.

## 6. The recurrence on toy inputs against the toy's measured values

Campaign A: the DESCENT.md recurrence, evaluated on the toy plan with inputs from fresh ideals only (95% upper limits), never predicts less than the measured failure rate, subtree nodes or subtree work of reached nodes at any of the 14 types with at least 500 reached nodes (also at the 95% limits); measured/predicted work is 0.04 to 0.96 (table: part=multilevel). On the toy it errs on the safe side, more so with depth, and is vacuous at the largest types: a check of direction at toy scale, not of production accuracy.

By address stratum (part=multilevel prints every ratio of reached to fresh): no stratum of at least 200 nodes has a 95% lower limit above 1.01; among reached nodes with M >= 2 relations 0 of 97269 failed; the count envelope holds in every stratum (largest ratio 0.8271).
