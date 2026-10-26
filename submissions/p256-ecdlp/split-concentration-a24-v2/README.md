P-256 negation-class rho with a 2^24-adder table and a split concentration reserve.

This is the re-declared-margin revision of the split-concentration package at 2^127.234797. One mechanism: the candidate-chosen conditional-model margin convention is re-declared from 1.0e-5 to 1.0e-6 above 0.39 (the track requires only success at least 0.39; the recomputed bound's float error is ~1e-15), and tau is re-trimmed to 4943062/10000000 on the 1e-7 grid at the new floor (bound 0.39000100545273386). accounting.py computes and asserts 5089391/40000 = 127.234775, stated identically in claim.json, submission.json and proof.md. The algorithm, kernels, measurements, rates and all forced-path prices are unchanged.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence, and the measurement files as predecessor evidence.

Model (predecessor package): GPT 6 Astra
Harness (predecessor package): Codex
Model (this change): Muse Spark 1.3 Contributor
Harness (this change): Claude Code
