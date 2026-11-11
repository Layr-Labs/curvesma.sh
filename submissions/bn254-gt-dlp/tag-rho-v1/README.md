# Tag-traced Frobenius-class rho on BN254 GT

Claim: time_log2 = **123.02** (W = R/20,000), success probability at least 0.39, peak memory
below 2^63 bytes, the precomputed multiset product table (2^47.35 units) included, no advice.

- [Algorithm, success argument and cost](proof.md)
- [Machine model, kernel listings and instruction ledger](machine.md)
- [The c1dot listing](c1dot.listing.txt)
- [Executable checks, production lock-step walk, toy runs, cycle law, ledger](experiment.py)

Built on the accepted BN254 package 7e58d165 (e-man07, 124.23), itself built on 2ece9f04
(internettrashh, 124.65): same class map, keys, twisted table and kernels. New: tag tracing
(Cheon-Hong-Kim) with c_1-only tags from a fused AArch64 listing, a plain 32-adding walk
whose bits come from the least key, a multiset product table, and an exact useless-cycle
law pricing the escapes. This is an analytic generic-attack bound, not a practical break.
