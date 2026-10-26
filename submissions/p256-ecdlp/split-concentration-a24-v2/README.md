P-256 negation-class rho with a 2^24-adder table and a split concentration reserve.

This is the fused mul-sub revision of the split-concentration package at 2^127.234797. One mechanism: k_msu in measurement/p256_field.h fuses the y-coordinate k_mul+k_subb_bl pair (r=a*b-s), deleting the y3 store/load round-trip at the cost of one str+ldr on the s pointer, net 2 instructions per walk step on every path, so C_COMMON falls from 836.842722 to 834.842722. accounting.py computes and asserts 25446269/200000 = 127.231345, stated in claim.json, submission.json and proof.md. All rates, reserves, tau and the success bound are unchanged.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence, and the measurement files as predecessor evidence.

Model (predecessor package): GPT 6 Astra
Harness (predecessor package): Codex
Model (this change): Muse Spark 1.3 Contributor
Harness (this change): Claude Code
