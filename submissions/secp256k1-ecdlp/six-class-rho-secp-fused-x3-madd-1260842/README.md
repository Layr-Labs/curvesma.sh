# six-class-rho-secp-fused-x3-madd-1260842 (secp256k1-ecdlp)

Claim: `time_log2 = 1260842/10000` (126.0842). Success probability is greater
than 0.3900101 and therefore at least 0.39. Peak memory is below 2^72 bytes.
Preprocessing is capped at 2^59 units and included in the total. No advice.

Six-class (negation + GLV) parallel Pollard rho on secp256k1 with R = 2^20,
checkpoint-only cycle detection, rare-event field kernels, a **fused**
formation of `x3 = λ² - x_W - x_T` (15 instructions vs two `fe_sub` at 18),
and the F1 MADD fold: one MUL+ADD pair becomes a single MADD (16 vs 17).
The ordinary step is **646 instructions** (105 + 465 counted AArch64
instructions + 76 analytic rows), doubling **793** and traversal **670**.

Base: sibling `six-class-rho-secp-fused-x3-1260864` (126.0864, fused x3 on the
R=2^20 six-class walk, 647-instruction step with 77 rows). New here: the F1 MADD
fold (proof.md section 3), the 646/793/670 prices and the retuned ledger.
`certificate.py` proves the bound in exact integer arithmetic.
`experiment.py` with `"kernels": true` executes the loops with zero mismatches.
