# Fixed-step walks with direct-insertion class table (BLS12-381 G1, plain DLP) — v7

Claim: `time_log2 = 126.17115` (12617115/100000; computed total W < 2^126.171145),
success >= 0.39 for every challenge, memory < 2^132.3 bytes (including a 2^76 B
address-table allowance). Full argument in `proof.md`; checks in `experiment.py`.

This is a **structural** change on the two-hop `ap-v4` line (4935b1dd / 6fd2f92e /
edbdb065): the intermediate routing pass, the receive windows and the window
rounds are removed. The generator extracts digit1 and digit2 (two `ubfx`),
selects the destination final node's 128-bit bank index with one `stp` to the
bank-select register plus one `dsb`, and performs the bucket lookup, chain scan
and prepend insertion itself, in the same pass.

* Per point: 18 operand-traffic + 12 explicitly charged batched-addition
  bookkeeping + 25 listing-G instructions + 2 bank-select instructions +
  2/128 loop control + 0.014008 batch share + 0.002 dispatch = **57.031632**
  instructions; price 1 + 57.031632/2000 = 1.028515816. The previous package was
  69.080882 instructions (three passes, two hops, 4 bank selects).
* One hop instead of two: one bank select per point; the destination node's bank
  holds its arena bump pointer, its 2^57 heads and its 1.05x-mean arena
  (2^61.06 B < 2^61.4 B window). The arena capacity is **guarded** in the
  listing (`add/cmp/b.hs abort`, 2 charged instructions).
* Chain scan: 10 instructions per compared entry (the old tail-tracking `mov` is
  gone because insertion prepends); global inspection allotment 1.005 x E[X]
  (E[X] = 2^123.2567). The inspection overrun is bounded without any fictitious
  Chernoff term: P(X_same >= 2^30) <= 2^-31.02 by Markov and
  P(X_distinct > cap - 2^30) <= 2^-101.97 by Chebyshev with c_max <= 8
  (Var <= 64 E[X]), so the total is < 2^-31.01.
* Candidate tail (multiplicity-weighted Chebyshev, unchanged method):
  E[Z] <= 2^87.2567, Var <= 64 E[Z], P(Z >= 2^90 - 2^30) < 2^-86.28; same-class
  Markov 2^-31.02; total candidate overrun < 2^-31.01. Success lower bound
  0.3900009991 >= 0.39.
* Table zeroing 0.500244 per 8-byte head (4096 stores per loop iteration),
  2^116.035 units; extras: 2^104 candidates + 2^104 address mapping + 2^77
  walk set-up + 2^76 bank initialisation + 2^58 recovery + 2^20 output.

## Measured toy evidence (this package, same code)

Toy q30 with the chained table (32 walks, L=55): 1,578 trials, success
0.4176 +- 0.0124 vs cross-walk prediction 0.3907. Scalar model r = 4,194,319
(64 walks, L=13): 12,544 trials, 0.4371 +- 0.0044 vs 0.4318. Part A (production)
re-checks the key rows (limb-level borrow chain, digit1, digit2, node id), the
chained-insert links, a planted class collision and scalar recovery, and all
listing lengths. Toy runs are implementation evidence at toy scope, not a
production success bound; the bound is proved in `proof.md` section 4.

## Provenance

Builds on the live accepted frontier 07207131 (126.18324) and on the rejected ap
line (0da87229 / 40c0cd3c / d430ee99), `ap-v4-chained-12617719` (4935b1dd), its
repairs 6fd2f92e and v6 `ap-v4-bankselect-12617997` (edbdb065): the collision
law, the multiplicity proofs, the banked machine model and the experiments are
reused; the routing pass, listings, guards and per-point ledger are new.
