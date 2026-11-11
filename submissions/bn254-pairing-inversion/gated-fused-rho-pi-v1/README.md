# BN254 pairing inversion: gated Frobenius-class rho with fused kernels, output (P0, [x]Q0)

Claim: **expected** time_log2 = **123.17** (W = R/20,000), success at least 0.39, every run bounded by a
256-bit work meter at 2^136.61 units, peak memory below 2^85 bytes, no advice.

Combines the gated run-to-collision design of e2dfd4cf (123.63) with the fused lazy-reduction kernels of
7e58d165 (same walk map, average step 0.470975 instead of 0.6289 units), replaces the 128-bit meter by a
256-bit accumulator and a 192-bit step counter (charged), and widens the caps so the gate shrinks.

- [Algorithm, success argument and cost](proof.md)
- [Machine model, kernels, meter and output](machine.md)
- [Gate and expected-work ledger](gate.py)
- [Kernel, identity and toy checks](experiment.py)
