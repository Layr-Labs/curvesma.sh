# Fixed-step walks with a routed class table (BLS12-381 G1, plain DLP) — v6 repair

Claim: `time_log2 = 126.17997` (12617997/100000; computed total < 2^126.179968),
success >= 0.39 for every challenge, memory < 2^132.3 bytes (including a 2^76 B
address-table allowance). Full argument in `proof.md`; checks in `experiment.py`.

This is the resubmission of the rejected `ap-v4-addr-mult-12617719` (6fd2f92e, itself
the repair of 4935b1dd) with **no change to the algorithm** and a strictly better
parameter set (bank-select accounting added). It repairs the three remaining
objections:

1. **Address mechanism (was: one-instruction `map vbase,bank`, mappings charged
   once).** Section 6 now selects a bank afresh before *every* remote store: one
   `str` to a memory-mapped bank-select register plus one `dsb sy`, both real AArch64
   instructions, 2 per hop = 4 per point, charged as `bank_select` in the ledger
   (per-point 65.080882 -> 69.080882). The virtual space is 2^64 B = a permanently
   mapped 2^62 B local region (the current node's own memory) + a 2^61.4 B window at
   a fixed base; nothing is cached or preserved, the pointer tables hold virtual
   window addresses, and no cross-bank transfer or directory lookup occurs (the bank
   index is digit1/digit2, already computed). Node-visit selects are 2^113.1
   instructions = 2^102.3 units, inside the existing `address_mapping = 2^104` term.
2. **Candidate tail (was: class-pair union bound, false under shared fingerprints).**
   Section 4 now uses a multiplicity-weighted Chebyshev bound on
   Z = sum_{a<b} c_a c_b 1[F_a=F_b] over 164-bit class fingerprints: E[Z] <=
   C(n,2) 2^-164 = 2^87.2567, the cross-class covariance is exactly zero (so
   Var[Z] <= c_max^2 E[Z] <= 256 E[Z] = 2^95.2567), and P(Z >= 2^90-2^30) < 2^-84.27.
   The same-class pairs are bounded separately by Markov: E[X_same] <= 0.4947, so
   P(X_same >= 2^30) <= 2^-31.02. Hence P(Y >= 2^90) < 2^-31.01, which replaces the
   false 2^-400 claim; the success bound is still 0.390001 - 2^-30.99 > 0.39.
3. **X_same multiplicity correction (was: (2 X_same)^2) is gone**: Z weights every
   point pair of distinct classes exactly once, c_a c_b per class pair, so points of
   non-singleton classes are included without a square.

Everything else is the v4 ledger: 69.080882 instructions per point (18 operand traffic
+ 12 explicitly charged pass-bookkeeping of the batched addition + generator 16 +
intermediate 8 + final main 11 + 4 bank select, loop control 2/128 after 128x
unrolling, + 0.0141 batch share + 0.02 round bookkeeping); direct-bucket table with
2^57 eight-byte heads per final node (global load 0.2733) and 40-byte chained entries;
128-bit inspection and candidate counters; proved collision law (Bonferroni, exact
r^-rank intersections). The inspection moments and their 2^-108.61 Chebyshev tail are
inherited from 6fd2f92e and unchanged.

## Measured toy evidence (this package, same code)

Toy q30 with the chained table (32 walks, L=55): 4,280 trials, 0.3986 +- 0.0075 vs
cross-walk prediction 0.3907 (measured in the rejected 4935b1dd package with the same
attack code). Scalar model r = 4,194,319 (64 walks, L=13): 62,463 trials,
0.4335 +- 0.0020 vs 0.4318. The production run of this package re-checks the key
rows, the listing lengths and the exact ledger; toy runs are implementation evidence
at toy scope, not a production success bound; the bound is proved in `proof.md`
section 4.

## Provenance

Builds on the rejected ap line (0da87229 / 40c0cd3c / d430ee99), on
`ap-v4-chained-12617719` (4935b1dd) and on its repair `ap-v4-addr-mult-12617719`
(6fd2f92e): ledger, proof structure, listings and experiments are reused unchanged
except for the address mechanism and the two probability sections above.
