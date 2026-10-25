# six-class-rho-secp-fused-x3-1260864 (secp256k1-ecdlp)

Claim: `time_log2 = 1260864/10000` (126.0864). Success probability is greater
than 0.3900101 and therefore at least 0.39. Peak memory is below 2^72 bytes.
Preprocessing is capped at 2^59 units and included in the total. No advice.

Six-class (negation + GLV) parallel Pollard rho on secp256k1 with R = 2^20,
checkpoint-only cycle detection, rare-event field kernels, and a **fused**
formation of `x3 = λ² - x_W - x_T` (15 instructions vs two `fe_sub` at 18).
The ordinary step is **647 instructions** (105 + 465 counted AArch64
instructions + 77 analytic rows), doubling **794** and traversal **671**.

Base: bradcush `736c001a` (126.0931). New here: `fe_sub2` in `kernels/gen_loops.py`,
regenerated loops (pass 2 = 465), the final-borrow failure row, and the retuned
ledger. `certificate.py` proves the bound in exact integer arithmetic.
`experiment.py` with `"kernels": true` executes the loops with zero mismatches.
