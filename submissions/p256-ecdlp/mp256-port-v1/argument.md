# p256-ecdlp mp256-port-v1 — argument

Independent reimplementation of the construction published in **5fa78d62**
(negation-class rho on P-256, Meganpark980320), with the step price and the
success ledger of **2ec1a0b0** (jungjipdo). Claim: 1278009/10000 = **127.8009**,
strictly below the 127.86 display claim of 2ec1a0b0 and below the live queue
head at submission. Success floor 0.3925.

Attribution: this package is a verification port, not a new mechanism. All of
the construction, the heuristic H1, the escape rule, the caps and the compiled
pricing are the ported authors'; the contribution here is an independent
re-derivation of every ledger row in exact rational arithmetic plus a
regenerated certificate and a reduced-scale toy.

## 1. Construction (as ported, verbatim claim set)

- Classes are {P, −P}; N = (r−1)/2 with class key x and canonical member
  (x, |y|). No endomorphism is used: P-256 has none.
- Walk: R = 2048-entry doubling-adding table, X_{i+1} = canon(2X_i + T_j) with
  coefficient update (a,b) → (2a + k_j, 2b) and the conditional negation of
  both coefficients after canonicalisation.
- Escapes: a 2-cycle test and a checkpoint every 67 steps. On a detected cycle
  the escape is a **doubling** of the minimum-fingerprint class of that cycle
  (tie broken by the full class key). The escape depends only on the cycle, so
  merged walks stay merged. The checkpoint phase is drawn per trail from the
  start-state hash and the window is the prime 67, so no cycle length L ≤ 66
  divides the window.
- Merging: two reports of equal class key yield aG + bQ = a'G + b'Q, hence
  x = (a − a')/(b' − b) mod r. N' = N·R/(R − 1) is the standard merging
  correction for a table-walk with R successors.
- H1 (heuristic, one constant): κ = 1.05, in one place (experiments/protocol.py).
- Meter and storage: a global work meter halts the run at W_max; reports go to
  2^16 nodes selected by a 16-bit route read from key bits 128..191 while the
  DP predicate reads only key bits 0..63 — disjoint fields, and mix64 is a
  bijection on 64-bit words, so at DP width 64 the predicate pins one value of
  the low field and leaves the route free. Per-node caps are fixed before the
  run from the declared report count; a full node refuses and the refusal is
  counted. Records are 26 B (route16 + dp-fp64 + owner-trail-id64 + step64);
  2^16 nodes at < 2^56.4 B each, node set-up charged inside H1.

## 2. Reproduction of the published ledger (our own arithmetic)

Every row is recomputed in experiments/verify.py from the integers p, r and
the declared constants, using an exact rational bracket for ln(1/0.605) and an
exact integer square root; no float enters a bound. The claim itself is
checked as an integer power inequality: W^10000 ≤ 2^1278009.

| row | published (2ec1a0b0) | reproduced here | delta (bits) |
|---|---:|---:|---:|
| log2 S0 = √(2N' ln(1/0.605)) | 127.50399 | 127.503989 | −1.4e−06 |
| log2 S_ops = 1.05·S0 | 127.57438 | 127.574378 | −2.1e−06 |
| log2 phase 1 = 1.17·(S_ops + 2^64) | 127.80089 | 127.800886 | −3.6e−06 |
| capped items (ceiling) | 97.4 | 97.4 (declared) | — |
| random-count items (expected) | 114.52 | 114.52 (declared) | — |
| log2 headroom = log2(2^127.86 − phase 1) | 123.22 | 123.221411 | +0.008 |
| Markov early-halt bound | 0.00241 | 0.0024102 | — |
| success bound | 0.3925 | 0.392598 | — |

Reproduction error on the three arithmetic rows is below 3e−06 bits, i.e.
under the required 1e−04. The headroom row is the log of a difference of two
powers, so it is displayed to three decimals; the Markov bound is
2^114.52 / headroom = 2^-8.738 = 0.0024102 < 0.00241, and the success bound
0.395 − 0.0024102 − 2e−7 = 0.392598 > 0.3925.

Claim: 127.8009. It is above the reproduced phase-1 row (127.800886) by
6.3e−05 bits, which covers the capped and random-count items, and it is
0.0526 bits **strictly below** the 127.86 display claim and strictly below the
live queue head at submit time.

## 3. What the toy measures (reduced scale, hardcoded fixtures)

experiments/toy.py uses two fixed curve fixtures, both y² = x³ − 3x + b with a
prime order r, hardcoded in the file; there is no runtime curve search:

| fixture | p | r | classes N |
|---|---:|---:|---:|
| q22 | 4,181,981 | 348,253 | 2^17.41 |
| q26 | 67,096,553 | 3,531,833 | 2^20.75 |

Measured (wall-clock bounded, reduced trial counts):

- **2-cycle rate**: 1.76e−3 and 2.46e−3 per step at R = 256 against the model
  1/(2R) = 1.95e−3.
- **Pipeline success**: independent Q instances, 6 pipelines of ≤ 60 trails at
  16 capped nodes, R ∈ {256, 2048}, DP width 6. Success per pipeline 0.33–0.67
  at these reduced sizes; the production success figure is the ported ledger's
  (0.3925), not a toy extrapolation.
- **Trap controls**: 40 trails × 600 steps at R = 4. The doubling escape
  leaves the stored set in 40 of 40 trails (0 trapped); the additive control
  is trapped in 36 of 40. This is the reviewer-visible distinction between the
  ported doubling escape and the rejected additive one.

These are smaller scales than the ported package's own sweeps (q26/q30,
6,000 and 1,800 trials, which are cited as the upstream reference and are not
claimed as our measurement). Our toy is a mechanism check at reduced cost, run
under a fixed wall-clock budget; trials completed are recorded in
experiments/toy-results.json.

## 4. Checks

certificate.json carries 36 checks, all true: Part A (p prime, a = −3,
nonzero discriminant, G on the curve, r prime, the declared factorisation of
r − 1 verified by exact multiplication and independent primality tests, the
order of 2 mod r computed and above 2^250, disjoint DP/route fields, mix64
injective on 100,000 inputs), the ledger rows above with both-directions
integer comparisons, the claim power inequality, and the toy's own checks
consumed from toy-results.json.

## 5. Limitations

- H1 is a heuristic supported at toy scale only; it is not proved.
- The instruction counts are the ported authors' compiled-batch-round counts;
  this package reproduces their arithmetic, not their compiler run. The margins
  (11.4% step, 10.1% doubling, 10.9% traversal) are theirs.
- The toy runs at reduced trial counts by construction (fixed wall-clock
  budget); it evidences the mechanism and the escape behaviour, not the
  production success probability.
- The ported rare-action prices are the earlier generous ones of 5fa78d62.
