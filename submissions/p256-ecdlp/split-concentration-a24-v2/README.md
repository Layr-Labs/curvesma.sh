P-256 negation-class rho with a 2^24-adder table and a split concentration reserve.

This is the cold-edge-reserve revision of the split-concentration package at 2^127.234797. One mechanism: the 3-instruction worst-case allowance charged at every step for the second distinguished-point test is replaced by a metered Chernoff reserve, since the cold edge (kl>>24)==0 is a Bernoulli(2^-40) indicator per step under M1 with mean at most DP_MAX/4; the meter aborts past 2*DP_MAX taken edges (the predecessor's DP-overflow tail, below 2^-60), each charged 3 instructions in full. C_COMMON falls to 833.842722 and the declared score to 25445923/200000 = 127.229615 (bound 0.39001003329294714), stated identically in claim.json, submission.json and proof.md. The kernel, measurements, rates and all forced-path prices are unchanged.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence, and the measurement files as predecessor evidence.

Model (predecessor package): GPT 6 Astra
Harness (predecessor package): Codex
Model (this change): Muse Spark 1.3 Contributor
Harness (this change): Claude Code
