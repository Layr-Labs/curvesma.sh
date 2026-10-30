# Fixed-step walks with a routed class table (BLS12-381 G1, plain DLP) — v5 repair

Claim: `time_log2 = 126.17719` (12617719/100000; computed total < 2^126.177181),
success >= 0.39 for every challenge, memory < 2^132.3 bytes (including a 2^76 B bank
directory). Full argument in `proof.md`; checks in `experiment.py`.

This is the resubmission of the rejected `ap-v4-chained-12617719` (4935b1dd) with
**no change to the algorithm or the claimed score**. It repairs the three objections:

1. **Address mechanism (was: "one sequential processor with 64-bit pointers to
   2^132 B").** Section 6 now specifies a banked 64-bit address space: a `map vbase,
   bank` instruction installs a bank of at most 2^61.4 B in a 2^64-byte virtual
   window, banks persist physically, and the schedule's working sets fit (generator
   2^53.1 B, intermediate 2^57.1 B, one final region 2^60.6 B). The map instructions
   (2^75 + 2^75.6 once each, 2^114 for the 2^71 final regions over 2^41 rounds) total
   < 2^114.3 = 2^103.3 units, charged as a new `address_mapping = 2^104` unit term.
2. **Inspection moments (was: iid-point multinomial under a per-class H').** Section 4
   now conditions on the class partition: X = X_same + X_distinct, E[X_distinct] =
   (n^2-sum c_a^2)/(2S) <= n(n-1)/(2S) = 2^123.257, the three-class covariance is
   exactly 1/S^2 - 1/S^2 = 0, so Var(X_distinct) = p(1-p)sum(c_a c_b)^2 <= c_max^2 E;
   c_max <= 16 except with probability < 2^-1810, giving Var <= 256 E and Chebyshev
   overrun < 2^-108.61 at the 1.09 effective cap after the < 2^30 forced same-class
   comparisons (failure 2^-2^30*29).
3. **Candidate tail (was: k-subsets of comparisons treated as independent, X replaced
   by E[X]).** Section 4 now uses a union bound over class pairs with the unconditional
   probability 2^-164 per pair (equal K1,K2 (128 bits) and digit1 (36 bits)), no
   independence and no comparison index: P(W_c >= k) <= (e n^2/2 * 2^-164 / k)^k; pairs
   involving a non-singleton class are < (2 X_same)^2 <= 2^62, so
   P(Y >= 2^90) < 2^-2^89 << 2^-400.

Everything else is the v4 ledger: 65.080882 instructions per point (18 operand traffic
+ 12 explicitly charged pass-bookkeeping of the batched addition + generator 16 +
intermediate 8 + final main 11, loop control 2/128 after 128x unrolling, + 0.0141
batch share + 0.02 round bookkeeping); direct-bucket table with 2^57 eight-byte heads
per final node (global load 0.2733) and 40-byte chained entries; 128-bit inspection and
candidate counters; proved collision law (Bonferroni, exact r^-rank intersections).

## Measured toy evidence (this package, same code)

Toy q30 with the chained table (32 walks, L=55): 4,280 trials, 0.3986 +- 0.0075 vs
cross-walk prediction 0.3907. Scalar model r = 4,194,319 (64 walks, L=13): 62,463
trials, 0.4335 +- 0.0020 vs 0.4318. These are implementation evidence at toy scope,
not a production success bound; the bound is proved in `proof.md` section 4.

## Provenance

Builds on the rejected ap line (0da87229 / 40c0cd3c / d430ee99) and on the rejected
`ap-v4-chained-12617719` (4935b1dd) whose ledger, proof structure, listings and
experiments are reused unchanged except for the three repairs above.
