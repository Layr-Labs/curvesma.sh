# Multilevel descents under the production node rules

40100 recursive descents on the reduced tower over the production base ring (130965 entered nodes, 89621 adaptively
reached prime-ideal nodes, depth up to 5), run by the node procedure of DESCENT.md. It adds three measurements to the
171-node recursive campaign (part=adaptive_stats) and supports the sequential hazard and the conditional child-work envelope
at toy scale. It does not establish them at production size.

## 1. Attempts after an earlier failed relation

| campaign | prior failed relations | completed | failed | censored | failure rate | Clopper-Pearson 95% upper |
|---|---|---|---|---|---|---|
| A | 0 | 128683 | 1769 | 0 | 0.0137 | 0.0143 |
| A | 1 | 1724 | 33 | 0 | 0.0191 | 0.0255 |
| A | 2 | 33 | 0 | 0 | 0.0000 | 0.0868 |

Campaign A: 1757 completed attempts after an earlier failed relation, 33 failed, 0 censored. At the first-attempt failure rate of the same parent type 37.85 failures were expected against 33 observed: ratio 0.87 [0.57, 1.17].

## 2. Deterministic-order roots

Campaign A: all 40100 runs use the production order (relations in cell order, children in the fixed order, first success ends the node); none is shuffled. 0 runs hit the node budget and 0 attempts are unresolved: nothing is censored. Roots are uniformly sampled ideals of every type and both sides (38000 prime-ideal roots, 2100 helper-control roots).

## 3. Conditional child work after a failure

Children entered after a failed relation of their parent, against (a) fresh ideals of the same class (bits, side, half of the bit interval) and (b) the slot envelope: the largest fresh mean over the types the recurrence admits in the child's window slot, which is what child_work charges. Ratios of sums with 95% interval or upper limit:

| campaign | children | subtree nodes vs fresh | subtree work vs fresh | failure vs fresh | subtree nodes vs envelope | subtree work vs envelope |
|---|---|---|---|---|---|---|
| A, after a failure | 2063 | 1.04 [1.01, 1.06] | 1.03 [1.00, 1.06] | 0.81 [0.52, 1.09] | 0.86 (<= 0.89) | 0.61 (<= 0.68) |
| A, no failure on the path | 86030 | 1.01 [1.00, 1.01] | 1.00 [0.99, 1.00] | 0.96 [0.84, 1.08] | 0.84 (<= 0.85) | 0.58 (<= 0.59) |

Work is the number of survivors factored in the whole subtree. Children reached after a failure cost a few percent more than
fresh ideals of their class and stay below the slot envelope.

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
  192. A toy ball holds about 3 to 15 relations on average against hundreds, so toy failure rates are far above the production
  bounds: a stress test of selection, not a measurement of production rates.
- Factorisation is exact (no ECM). The set of ideals without logarithm is a hash; with one fixed set there are few of them, so
  reached helper nodes repeat them across runs and are listed apart (section 7).
- Nothing here proves that the ratios of section 3 hold at production size. They show that, under the production rules at this
  scale, reaching an ideal through earlier failures changes its failure rate and subtree work by a few percent at most.

## 6. Reached versus fresh by address stratum

Sum over the reached nodes of a stratum divided by the sum expected from fresh ideals of each node's class (95% interval
including the error of the fresh controls); 1.00 means no selection. P(M<1), own survivors and helpers are inputs of the
recurrence, the last three its outcomes. part=multilevel also prints the strata after a successful sibling, at depth >= 3 and
after two failures.

Campaign A, strata in the order no failure on the path (86030 nodes); after a failed relation of the parent (2063 nodes); below an ancestor entered after a failure (1528 nodes): P(M<1) 0.96 [0.84, 1.08], 0.83 [0.54, 1.12], 0.93 [0.61, 1.24]; own survivors 0.99 [0.99, 1.00], 1.02 [1.01, 1.04], 0.98 [0.95, 1.01]; helpers 1.02 [1.00, 1.04], 1.05 [1.00, 1.11], 1.01 [0.94, 1.07]; failure 0.96 [0.84, 1.08], 0.81 [0.52, 1.09], 0.94 [0.63, 1.25]; subtree nodes 1.01 [1.00, 1.01], 1.04 [1.01, 1.06], 1.04 [1.02, 1.06]; subtree work 1.00 [0.99, 1.00], 1.03 [1.00, 1.06], 1.00 [0.97, 1.04].

Node failure by relation count M among reached nodes: A: M=1 29/6849, M>=2 0/81020.

Count envelope (alpha = 1/2): the largest ratio of mean t^M to the mean envelope over t = 0..0.8 and all strata is 0.827 (A).

## 7. Helper nodes, repeats

- A: 1238 reached helper nodes on 876 distinct special-q, 21 failed; 2100 fresh controls, 121 failed. 6 repeat nodes. Raw run files are not shipped (sha256 6d832dfba462b083).