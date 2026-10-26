# six-class-rho-secp-fused-x3-addrbase-1260842 (secp256k1-ecdlp)

Claim: `time_log2 = 1260842/10000` (126.0842). Success probability is greater
than 0.3900101 and therefore at least 0.39. Peak memory is below 2^72 bytes.
Preprocessing is capped at 2^59 units and included in the total. No advice.

Six-class (negation + GLV) parallel Pollard rho on secp256k1 with R = 2^20,
checkpoint-only cycle detection, rare-event field kernels, a **fused**
formation of `x3 = λ² - x_W - x_T` (15 instructions vs two `fe_sub` at 18),
and a dual-base table address (CSEL between base and base + 256, row 6 → 5).
The ordinary step is **646 instructions** (105 + 465 counted AArch64
instructions + 76 analytic rows), doubling **793** and traversal **670**.

Base: sibling `six-class-rho-secp-fused-x3-1260864` (126.0864). New here: the
h/DP-test/table-address row cut and the retuned ledger. `certificate.py`
proves the bound in exact integer arithmetic. `experiment.py` with
`"kernels": true` executes the loops with zero mismatches.
