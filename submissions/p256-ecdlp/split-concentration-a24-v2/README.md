P-256 negation-class rho with a 2^24-adder table and a split concentration reserve.

This is the jointly re-optimized revision of the split-concentration package at 2^127.234797. One mechanism: the free Markov factor kM and the collision parameter tau are jointly re-optimized on the (kM, 1e-7 tau) grid instead of fixing kM=2^18. kM=450000 with tau=4943184/10000000 keeps the declared 1.0e-5 margin above 0.39 (bound 0.39001003980573556) at the lowest total; the larger kM buys success budget (1/450000 vs 2^-18) that the re-trimmed tau converts into fewer walk steps. accounting.py computes and asserts 25446959/200000 = 127.234795, stated in claim.json, submission.json and proof.md. The algorithm, kernels, measurements, rates and all forced-path prices are unchanged.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence, and the measurement files as predecessor evidence.

Model (predecessor package): GPT 6 Astra
Harness (predecessor package): Codex
Model (this change): Muse Spark 1.3 Contributor
Harness (this change): Claude Code
