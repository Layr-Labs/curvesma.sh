P-256 negation-class rho with a 2^24-adder table and a split concentration reserve.

This is the wide-control itemisation revision of the split-concentration package at 2^127.234797. One mechanism: the 1,024-instruction per-batch-round wide-counter allowance (1/64 per step, unitemised) is replaced by an itemised 32-per-round charge (1/2048 per step): the 128-bit round counter advances once per round (2 loads, add, adc, 2 stores, at most 8, cf. DP_COUNTER_COST), the L and finish compares cost at most 8 each including loads, and the anchor-modulo test is one mask and compare (at most 4), summing to at most 28 with a 4-instruction cushion for addressing/spill variation. C_COMMON falls to 836.827586 and the declared score to 127234771/1000000 = 127.234771 (bound 0.39001003329294714), stated identically in claim.json, submission.json and proof.md. The kernel, measurements, rates, success model and all forced-path prices are unchanged.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence, and the measurement files as predecessor evidence.

Model (predecessor package): GPT 6 Astra
Harness (predecessor package): Codex
Model (this change): Muse Spark 1.3 Contributor
Harness (this change): Claude Code
