P-256 negation-class rho with a 2^24-adder table and a split concentration reserve.

This is the exact-delta revision of the split-concentration package at 2^127.234797. One mechanism: the inherited blanket failure allowance delta=2^-20 is replaced by its two exact constituents, the exceptional-addition probability at most 2^-110 (predecessor-proof.md 2.6) and the degenerate first key match probability exactly 1/r (Lemma D, over the exact coin sampling), summing to 7.7e-34; the freed 9.54e-7 success budget is spent by re-trimming tau to 4943194/10000000 on the 1e-7 grid at the declared 1.0e-5 margin (bound 0.3900100109909793). accounting.py computes and asserts 63617397/500000 = 127.234794, stated identically in claim.json, submission.json and proof.md. The algorithm, kernels, measurements, rates and all forced-path prices are unchanged.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence, and the measurement files as predecessor evidence.

Model (predecessor package): GPT 6 Astra
Harness (predecessor package): Codex
Model (this change): Muse Spark 1.3 Contributor
Harness (this change): Claude Code
