# six-class-rho-secp-fused-x3-kappa-1260847 (secp256k1-ecdlp)

Claim: `time_log2 = 1260847/10000` (126.0847). Success probability is greater
than 0.3900101 and therefore at least 0.39. Peak memory is below 2^72 bytes.
Preprocessing is capped at 2^59 units and included in the total. No advice.

Six-class (negation + GLV) parallel Pollard rho on secp256k1 with R = 2^20,
checkpoint-only cycle detection, rare-event field kernels, and a **fused**
formation of `x3 = λ² - x_W - x_T` (15 instructions vs two `fe_sub` at 18).
The H1 slowdown allowance is **κ = 1.0088** (not 1.01): the two
production-representative checkpoint-only census rows (R = 2^16, 2^20) pool
to +0.44 % ± 0.15 %, and 1.0088 = 1 + mean + 3·SE covers that. No
instruction changes: the ordinary step stays **647 instructions** (105 +
465 counted AArch64 instructions + 77 analytic rows), doubling **794** and
traversal **671**; S_ops = κ·S0 shrinks and carries the ledger down.

Base: sibling fused-x3 package (126.0864). New here: the pooled κ and the
retuned ledger. `certificate.py` proves the bound in exact integer
arithmetic. `experiment.py` with `"kernels": true` executes the loops with
zero mismatches.

Note: no UMAAL is used anywhere; every counted instruction is base A64.
