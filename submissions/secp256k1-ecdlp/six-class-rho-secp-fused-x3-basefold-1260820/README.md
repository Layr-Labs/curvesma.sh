# six-class-rho-secp-fused-x3-basefold-1260820 (secp256k1-ecdlp)

Claim: `time_log2 = 1260820/10000` (126.0820). Success probability is greater
than 0.3900101 and therefore at least 0.39. Peak memory is below 2^72 bytes.
Preprocessing is capped at 2^59 units and included in the total. No advice.

Six-class (negation + GLV) parallel Pollard rho on secp256k1 with R = 2^20,
checkpoint-only cycle detection, rare-event field kernels, a **fused**
formation of `x3 = λ² - x_W - x_T` (15 instructions vs two `fe_sub` at 18),
and the arrival-pass base fold: the two base-setup memory instructions are
already the counted loop heads, so arrival traffic drops from 3 to 1.
The ordinary step is **645 instructions** (105 + 465 counted AArch64
instructions + 75 analytic rows), doubling **792** and traversal **669**.

Base: sibling `six-class-rho-secp-fused-x3-1260864` (126.0864, fused x3 on the
R=2^20 six-class walk, 647-instruction step). New here: the arrival traffic
fold (proof.md section 3), the 645/792/669 prices and the retuned ledger.
`certificate.py` proves the bound in exact integer arithmetic.
`experiment.py` with `"kernels": true` executes the loops with zero mismatches.
