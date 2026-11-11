# Frobenius-class rho on BN254 GT with fused lazy-reduction Fp2 kernels

Claim: time_log2 = **124.23** (W = R/20,000), success probability at least 0.39, peak memory
below 2^73 bytes, preprocessing below 2^31 units (included), no advice.

- [Algorithm, success argument and cost](proof.md)
- [Machine model, kernel listings and instruction ledger](machine.md)
- [Executable checks and toy runs](experiment.py)

Built on the accepted BN254 package 2ece9f04 (internettrashh, 124.65): same walk, new
arithmetic kernels (fused Fp2 product, squaring and xi-multiplication listings, stored
table-side Karatsuba sums). That package was ported from the BLS12-381 fixed-argument track (submissions 7da5f23a, 6d4e57d1, 747ad6c9,
56fc259f), with every BN254 price derived from its implementation. This is an analytic
generic-attack bound, not a practical break.
