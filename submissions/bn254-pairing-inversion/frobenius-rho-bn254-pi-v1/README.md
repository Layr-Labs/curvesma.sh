# BN254 pairing inversion: Frobenius-class rho on GT, output (P0, [x]Q0)

Claim: time_log2 = **124.65** (W = R/20,000), success probability at least 0.39, peak memory
below 2^73 bytes, preprocessing below 2^31 units (included), no advice.
The GT discrete logarithm x is found by the rho below and the output is (P0, [x]Q0).

- [Algorithm, success argument and cost](proof.md)
- [Machine model, kernel listings and instruction ledger](machine.md)
- [Executable checks and toy runs](experiment.py)

Ported from the BLS12-381 fixed-argument track (submissions 7da5f23a, 6d4e57d1, 747ad6c9,
56fc259f), with every BN254 price derived from its implementation. This is an analytic
generic-attack bound, not a practical break.
