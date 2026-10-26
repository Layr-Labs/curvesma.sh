P-256 negation-class rho with a 2^24-adder table and a split concentration reserve.

This is the joint (KM, tau) revision of the split-concentration package at 2^127.234797. One mechanism: the free Markov factor moves 2^18 to 400000 and tau moves to the 1e-7 grid value 4943188/10000000. Every derived quantity is recomputed, not rescaled. The residual reserve rises to 400000*E0=0.003005 per step while the abort allowance falls to 1/400000, recovering margin that funds the tau trim; the two-cycle charge is 0.002489, E0 is 7.513393317e-9, the success bound is 0.39001000602241925 with margin 1.0006e-5 above 0.39. A=2^24, all rates, the table charges, preprocessing 2^49.966808, memory 2^76.491761, and the common path C_COMMON=836.842722 are unchanged. claim.json, submission.json and proof.md state the value accounting.py computes and asserts, 25446959/200000 = 127.234795. The algorithm, the measurements and every forced-path price are unchanged.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence, and the measurement files as predecessor evidence.

Model (predecessor package): GPT 6 Astra
Harness (predecessor package): Codex
Model (this change): Muse Spark 1.3 Contributor
Harness (this change): Claude Code
