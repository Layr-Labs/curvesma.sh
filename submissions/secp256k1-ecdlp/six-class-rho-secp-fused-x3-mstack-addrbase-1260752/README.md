# six-class-rho-secp-fused-x3-mstack-addrbase-1260752 (secp256k1-ecdlp)

Claim: `time_log2 = 1260752/10000` (126.0752). Success probability is greater
than 0.3900101 and therefore at least 0.39. Peak memory is below 2^72 bytes.
Preprocessing is capped at 2^59 units and included in the total. No advice.

Six-class (negation + GLV) parallel Pollard rho on secp256k1 with R = 2^20,
checkpoint-only cycle detection, rare-event field kernels, a **fused**
formation of `x3 = λ² - x_W - x_T` (15 instructions vs two `fe_sub` at 18),
the full F1 MADD fold (row 17 → 14), the K_cp/length co-load (row 9 → 8),
and the dual-base table address (row 6 → 5: base + 256 precomputed once per
batch in `x14`). The ordinary step is **642 instructions** (105 + 465
counted AArch64 instructions + 72 analytic rows), doubling **789** and
traversal **666**.

Base: the reviewed mergestack (643-step, 126.0775). New here: the dual-base
address and the retuned ledger. The prior 642-step tree was rejected only on
quoted ledger numerals (phase-1/W display values above the true bounds);
this tree quotes the exact bounds the initial run reports. `certificate.py`
proves the bound in exact integer arithmetic. `experiment.py` with
`"kernels": true` executes the loops with zero mismatches.

Note: no UMAAL is used anywhere; every counted instruction is base A64.
